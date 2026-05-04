"""
hgvs_normalize.py
HGVS 표준화 모듈 + gold 비교 재계산

Step 1: Rule-based lite normalization (항상 적용)
  - transcript prefix 제거 (NM_xxx.x:c. → c.)
  - 삭제/삽입 염기 strip (c.724_727delTGCT → c.724_727del)
  - 단일 위치 dup 정리 (c.5dupA → c.5dup)
  - 대문자 표준화

Step 2: Mutalyzer normalization (transcript ID 있는 경우)
  - 3' rule shifting
  - API: GET https://mutalyzer.nl/api/normalize/{description}

Usage:
  python hgvs_normalize.py          # normalization 후 gold 비교
  python hgvs_normalize.py --dry    # normalization 결과만 출력 (비교 생략)
"""

import asyncio
import csv
import json
import re
import sys
import urllib.parse
from pathlib import Path

import aiohttp

_HERE      = Path(__file__).resolve().parent.parent
GOLD_CSV   = _HERE / "results" / "evaluation_detail.csv"
EXTR_JSONL = _HERE / "data" / "extracted_variants.jsonl"
AVAIL_PATH = _HERE / "data" / "pmc_availability.jsonl"
MUTALYZER  = "https://mutalyzer.nl/api/normalize/{}"


# ── Rule-based normalization ───────────────────────────────────────────────────
_TRANSCRIPT_RE  = re.compile(r"^[A-Z]+_[\d.]+\s*:\s*", re.IGNORECASE)
_DEL_BASES_RE   = re.compile(r"(del)[ACGT]+$", re.IGNORECASE)
_DUP_BASES_RE   = re.compile(r"(dup)[ACGT]+$", re.IGNORECASE)
_INS_BASES_RE   = re.compile(r"(ins)([ACGT]+)$", re.IGNORECASE)   # keep inserted seq

def hgvs_lite(hgvs: str) -> str:
    """빠른 rule-based 정규화. transcript prefix와 삭제염기 제거."""
    if not hgvs:
        return ""
    h = hgvs.strip()
    # 1. transcript prefix 제거: NM_xxx.x:c. → c.
    h = _TRANSCRIPT_RE.sub("", h)
    # 2. del 뒤 염기 제거: c.724_727delTGCT → c.724_727del
    h = _DEL_BASES_RE.sub(r"\1", h)
    # 3. dup 뒤 단순 염기 제거: c.5dupA → c.5dup
    h = _DUP_BASES_RE.sub(r"\1", h)
    # 4. 표준화: ACGT 대문자, c/p 소문자
    h = re.sub(r"([cp])\.", lambda m: m.group(1).lower() + ".", h)
    return h.strip()


# ── Protein HGVS normalization ────────────────────────────────────────────────
_AA1 = {
    "a":"ala","r":"arg","n":"asn","d":"asp","c":"cys","q":"gln","e":"glu",
    "g":"gly","h":"his","i":"ile","l":"leu","k":"lys","m":"met","f":"phe",
    "p":"pro","s":"ser","t":"thr","w":"trp","y":"tyr","v":"val","x":"xaa",
}
# 3-letter AA 집합 (소문자)
_AA3 = set(_AA1.values())

def _expand_aa1(aa: str) -> str:
    """단일 1-letter AA → 3-letter. 이미 3-letter이면 그대로."""
    lo = aa.lower()
    if lo in _AA1:
        return _AA1[lo]
    if lo in _AA3:
        return lo
    return lo

def hgvs_protein_norm(prot: str) -> str:
    """
    p.HGVS 표준화 (비교용 canonical form)
    변환 목록:
      1. 외부 괄호 제거: p.(Q568fs*3) → p.Q568fs*3
      2. p. prefix 소문자
      3. * → ter
      4. fs*N / fsterN → fsterN
      5. 1-letter AA → 3-letter AA (위치 앞뒤)
      6. 프레임시프트 중간 AA 제거 (비교용):
         p.gln568argfster3 → p.gln568fster3
         (gold에 중간 AA 있고 pred에 없어도 매칭되도록)
    """
    if not prot:
        return ""
    h = prot.strip()

    # 1. 외부 p.( ) 괄호 제거
    h = re.sub(r"^[Pp]\.\((.+)\)$", r"p.\1", h)
    # p. prefix 정규화
    h = re.sub(r"^[Pp]\.", "p.", h)
    if not h.startswith("p."):
        h = "p." + h

    body = h[2:]  # "p." 이후

    # 2. * → ter (fs* 포함)
    body = body.replace("*", "ter")

    # 3. fster → fster (이미 ok), fs → 그대로
    # fs[숫자] → fster[숫자]
    body = re.sub(r"fs(\d)", r"fster\1", body)

    # 4. 대소문자 통일 → 소문자
    body = body.lower()

    # 5. 1-letter AA 확장
    # 패턴: [단일 AA 1글자][숫자] 앞부분
    def expand_leading(m):
        aa, pos = m.group(1), m.group(2)
        return _expand_aa1(aa) + pos
    body = re.sub(r"^([a-z])(\d+)", expand_leading, body)

    # 패턴: [숫자][단일 AA 1글자] 뒷부분 (missense 두 번째 AA)
    def expand_trailing(m):
        pos, aa = m.group(1), m.group(2)
        return pos + _expand_aa1(aa)
    body = re.sub(r"(\d+)([a-z])(?=fster|ter|$|\d)", expand_trailing, body)

    # 6. 프레임시프트 중간 AA 제거 (비교용 canonical):
    #    p.gln568argfster3 → p.gln568fster3
    #    패턴: [3-letter AA][pos][3-letter AA][fster|fs][숫자]
    #    → [3-letter AA][pos][fster][숫자]
    def strip_fs_intermediate(m):
        leading_aa, pos, _, fs_part = m.group(1), m.group(2), m.group(3), m.group(4)
        return leading_aa + pos + fs_part
    _aa3_pat = "|".join(sorted(_AA3, key=len, reverse=True))
    body = re.sub(
        rf"({_aa3_pat})(\d+)({_aa3_pat})(fster\d+)",
        strip_fs_intermediate, body
    )

    return "p." + body


def hgvs_norm_for_cmp(hgvs: str) -> str:
    """비교용: lite normalization + 모두 소문자"""
    return hgvs_lite(hgvs).lower()


# ── Mutalyzer async normalization ─────────────────────────────────────────────
async def mutalyzer_normalize(
    session: aiohttp.ClientSession,
    hgvs: str,
    semaphore: asyncio.Semaphore,
) -> str:
    """
    transcript ID가 있는 HGVS → Mutalyzer API로 표준화.
    실패하면 lite normalization 결과 반환.
    """
    lite = hgvs_lite(hgvs)
    # transcript ID 없으면 API 호출 불가
    if not re.match(r"[A-Z]+_[\d.]", hgvs, re.IGNORECASE):
        return lite

    encoded = urllib.parse.quote(hgvs, safe="")
    url     = MUTALYZER.format(encoded)
    async with semaphore:
        try:
            async with session.get(
                url, timeout=aiohttp.ClientTimeout(total=15)
            ) as resp:
                if resp.status != 200:
                    return lite
                data = await resp.json(content_type=None)
                norm = data.get("normalized_description", "")
                if norm and not data.get("errors"):
                    # transcript prefix 제거 후 반환
                    return hgvs_lite(norm)
        except Exception:
            pass
    return lite


async def normalize_all(hgvs_list: list[str]) -> dict[str, str]:
    """list of HGVS → {original: normalized} dict"""
    sem = asyncio.Semaphore(5)
    async with aiohttp.ClientSession() as session:
        tasks = {h: mutalyzer_normalize(session, h, sem) for h in set(hgvs_list)}
        results = await asyncio.gather(*tasks.values())
    return dict(zip(tasks.keys(), results))


# ── 평가 ──────────────────────────────────────────────────────────────────────
def run_eval(dry: bool = False):
    # gold 로드
    gold_by_doi: dict[str, list[dict]] = {}
    gold_by_pmid: dict[str, list[dict]] = {}
    with open(GOLD_CSV) as f:
        for row in csv.DictReader(f):
            doi  = row["doi"].strip().lower()
            pmid = row["pmid"].strip()
            if doi:
                gold_by_doi.setdefault(doi, []).append(row)
            if pmid:
                gold_by_pmid.setdefault(pmid, []).append(row)

    # extracted 로드
    extracted = []
    with open(EXTR_JSONL) as f:
        for line in f:
            extracted.append(json.loads(line))

    # pmid→doi 매핑
    pmid2doi: dict[str, str] = {}
    with open(AVAIL_PATH) as f:
        for line in f:
            r = json.loads(line)
            if r.get("doi") and r.get("pmid"):
                pmid2doi[r["pmid"]] = r["doi"].strip().lower()

    # 모든 HGVS 수집 (Mutalyzer 호출 대상)
    all_hgvs = set()
    for e in extracted:
        for f in e.get("findings", []):
            for k in ("hgvs_coding", "hgvs_protein"):
                if f.get(k):
                    all_hgvs.add(f[k])
    for rows in list(gold_by_doi.values()) + list(gold_by_pmid.values()):
        for row in rows:
            for k in ("HGVS_coding_gold", "HGVS_protein_gold"):
                if row.get(k):
                    all_hgvs.add(row[k])

    print(f"HGVS 정규화 대상: {len(all_hgvs)}개 (Mutalyzer 호출 포함)")
    norm_map = asyncio.run(normalize_all(list(all_hgvs)))

    # normalization 변환 요약
    changed = [(o, n) for o, n in norm_map.items() if hgvs_lite(o) != n]
    print(f"  Mutalyzer로 추가 변환: {len(changed)}개")
    for o, n in changed[:5]:
        print(f"    {o!r:40s} → {n!r}")

    if dry:
        print("\n--dry 모드: 비교 생략.")
        return

    # 필드별 매칭
    FIELDS    = ["HGVS_coding", "HGVS_protein"]
    GOLD_COL  = {"HGVS_coding": "HGVS_coding_gold", "HGVS_protein": "HGVS_protein_gold"}
    PRED_KEY  = {"HGVS_coding": "hgvs_coding",      "HGVS_protein": "hgvs_protein"}

    def nrm(s):
        return norm_map.get(s, hgvs_lite(s)).lower()

    stats = {f: {"tp": 0, "fp": 0, "fn": 0} for f in FIELDS}
    paper_details = []

    matched = 0
    for e in extracted:
        pmid = e["pmid"]
        doi  = pmid2doi.get(pmid, "")
        gold_rows = gold_by_pmid.get(pmid) or gold_by_doi.get(doi)
        if not gold_rows:
            continue
        matched += 1

        for field in FIELDS:
            gold_vals = {nrm(r[GOLD_COL[field]]) for r in gold_rows if r[GOLD_COL[field]].strip()}
            pred_vals = {nrm(f.get(PRED_KEY[field], "")) for f in e.get("findings", []) if f.get(PRED_KEY[field], "").strip()}

            tp = len(gold_vals & pred_vals)
            fp = len(pred_vals - gold_vals)
            fn = len(gold_vals - pred_vals)
            stats[field]["tp"] += tp
            stats[field]["fp"] += fp
            stats[field]["fn"] += fn

        paper_details.append({
            "pmid":   pmid,
            "title":  e.get("title", "")[:55],
            "gold_hgvs": sorted({nrm(r["HGVS_coding_gold"]) for r in gold_rows if r["HGVS_coding_gold"].strip()}),
            "pred_hgvs": sorted({nrm(f.get("hgvs_coding","")) for f in e.get("findings",[]) if f.get("hgvs_coding","")}),
        })

    print(f"\nMatched papers: {matched}/{len(extracted)}")
    print(f"\n{'Field':<16} {'Prec':>7} {'Rec':>7} {'F1':>7}  TP  FP  FN")
    print("-" * 55)
    for field in FIELDS:
        s = stats[field]
        tp, fp, fn = s["tp"], s["fp"], s["fn"]
        p  = tp / (tp + fp) if tp + fp else 0
        r  = tp / (tp + fn) if tp + fn else 0
        f1 = 2*p*r/(p+r) if p+r else 0
        print(f"{field:<16} {p:>7.1%} {r:>7.1%} {f1:>7.1%} {tp:>3} {fp:>3} {fn:>3}")

    print(f"\n{'='*55}")
    print("논문별 HGVS_coding 상세:")
    for p in paper_details:
        g, pr = set(p["gold_hgvs"]), set(p["pred_hgvs"])
        tp = g & pr
        print(f"\n  PMID {p['pmid']}: {p['title']}")
        print(f"    gold: {sorted(g)}")
        print(f"    pred: {sorted(pr)}")
        if tp:  print(f"    ✓ match: {sorted(tp)}")
        if g-pr: print(f"    ✗ 누락: {sorted(g-pr)}")
        if pr-g: print(f"    + 추가: {sorted(pr-g)}")


def apply_normalization():
    """
    --apply 모드: extracted_variants.jsonl 전체에 HGVS normalization 적용 후
    extracted_variants_normalized.jsonl 로 저장
    """
    out_path = _HERE / "data" / "extracted_variants_normalized.jsonl"

    records = []
    with open(EXTR_JSONL) as f:
        for line in f:
            records.append(json.loads(line))

    # 모든 HGVS 수집
    all_hgvs = set()
    for e in records:
        for finding in e.get("findings", []):
            for k in ("hgvs_coding", "hgvs_protein"):
                if finding.get(k):
                    all_hgvs.add(finding[k])

    print(f"HGVS normalization: {len(all_hgvs)}개...")
    norm_map = asyncio.run(normalize_all(list(all_hgvs)))

    with open(out_path, "w", encoding="utf-8") as f:
        for e in records:
            for finding in e.get("findings", []):
                for k in ("hgvs_coding", "hgvs_protein"):
                    if finding.get(k):
                        finding[k] = norm_map.get(finding[k], hgvs_lite(finding[k]))
            f.write(json.dumps(e, ensure_ascii=False) + "\n")

    print(f"저장 완료: {out_path} ({len(records):,}건)")


if __name__ == "__main__":
    if "--apply" in sys.argv:
        apply_normalization()
    else:
        dry = "--dry" in sys.argv
        run_eval(dry=dry)
