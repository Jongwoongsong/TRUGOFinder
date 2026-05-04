"""
11_extract_variants.py
PMC full-text에서 GOF variant 정보 구조화 추출

입력: data/fulltext/{PMCID}.json
출력: data/extracted_variants.jsonl

추출 필드 (per finding):
  gene, disease, omim, hgvs_coding, hgvs_protein,
  gof_effect, class (H/N), note, evidence_text

Usage:
  python 11_extract_variants.py --test-gold   # gold 24건 테스트
  python 11_extract_variants.py               # 전체 3,754건
  python 11_extract_variants.py --concurrency 8
"""

import argparse
import asyncio
import json
import os
import re
import logging
import re
from pathlib import Path

import aiohttp

_HERE     = Path(__file__).resolve().parent.parent
FT_DIR    = _HERE / "data" / "fulltext"
AVAIL_PATH= _HERE / "data" / "pmc_availability.jsonl"
OUT_PATH  = _HERE / "data" / "extracted_variants.jsonl"
LOG_PATH  = _HERE / "data" / "step2_extraction.log"
NEEDS_REVIEW_PATH = _HERE / "data" / "needs_review.jsonl"

VLLM_BASE = os.environ.get("VLLM_BASE", "http://localhost:8100/v1")
MODEL     = "Qwen/Qwen3-32B-AWQ"

# gold PMIDs
GOLD_PMIDS = {
    "18760763","39066985","28904206","21565291","29514032","36856110",
    "39028950","33816491","12692554","27588951","27059040","35484142",
    "25488980","39623139","30514661","35122023","28428218","28338294",
    "35792400","25817014","34667213","24325359","17403716","37675773",
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(), logging.FileHandler(LOG_PATH)],
)
log = logging.getLogger(__name__)

# ── JSON Schema ────────────────────────────────────────────────────────────────
FINDING_SCHEMA = {
    "type": "object",
    "properties": {
        "gene":         {"type": "string", "description": "Gene symbol (e.g. TP53, BRCA1)"},
        "disease":      {"type": "string", "description": "Disease or phenotype name"},
        "omim":         {"type": ["string", "null"], "description": "OMIM ID (e.g. 601399), null if not found"},
        "hgvs_coding":  {"type": ["string", "null"], "description": "HGVS coding variant (e.g. c.1234A>G), null if not found"},
        "hgvs_protein": {"type": ["string", "null"], "description": "HGVS protein change (e.g. p.Arg123Gly), null if not found"},
        "gof_effect":   {"type": "string",  "description": "Mechanism of gain-of-function (1-2 sentences)"},
        "class":        {"type": "string",  "enum": ["H", "N"],
                         "description": "H = direct GOF (constitutive activation, enhanced signaling, increased stability); N = toxic aggregation, dominant negative, uncertain, or inferred-only GOF"},
        "is_human":     {"type": "boolean",
                         "description": "True if this is a human HGNC gene AND evidence involves human patients/cells/tissue. False if evidence is ONLY from model organisms (mouse, yeast, C. elegans, plant, bacteria, etc.)"},
        "note":         {"type": ["string", "null"], "description": "Additional notes (e.g. dominant negative, C-terminal truncation)"},
        "evidence_text":{"type": "string",  "description": "Direct verbatim quote from paper supporting GOF classification (max 200 chars)"},
    },
    "required": ["gene", "disease", "omim", "hgvs_coding", "hgvs_protein",
                 "gof_effect", "class", "is_human", "note", "evidence_text"],
    "additionalProperties": False,
}

RESPONSE_SCHEMA = {
    "type": "object",
    "properties": {
        "findings": {
            "type": "array",
            "items": FINDING_SCHEMA,
            "minItems": 0,
            "maxItems": 10,
        }
    },
    "required": ["findings"],
    "additionalProperties": False,
}

# ── 프롬프트 ──────────────────────────────────────────────────────────────────
SYSTEM_PROMPT = """/no_think
You are a biomedical curation expert specializing in gain-of-function (GOF) mutations.
Extract ALL GOF variant findings from the provided paper sections.

Return ONLY valid JSON in this exact format (no markdown, no extra text):
{
  "findings": [
    {
      "gene": "GENE_SYMBOL",
      "disease": "disease name",
      "omim": "123456 or empty string",
      "hgvs_coding": "c.XXX or empty string",
      "hgvs_protein": "p.XXX or empty string",
      "gof_effect": "mechanism description",
      "class": "H or N",
      "is_human": true,
      "note": "additional notes or empty string",
      "evidence_text": "verbatim quote from paper"
    }
  ]
}

Rules:
- gene: HGNC symbol only (e.g. TP53, BRCA1). Use the approved human gene symbol.
- hgvs_coding: c.XXXX format only; empty string if not stated
- hgvs_protein: p.XXXX format only; empty string if not stated
- omim: numeric ID only (e.g. "601399"); empty string if not stated
- is_human: true ONLY if (1) the gene is an HGNC-approved human gene symbol AND (2) the experimental evidence directly involves human patients, human-derived cells, or human tissue. Set false if the evidence is exclusively from model organisms (mouse, rat, yeast, C. elegans, Drosophila, zebrafish, plant, bacteria, or other non-human systems). If a human gene is studied only in a mouse model with no human patient data, set false.

- class "H": DIRECT gain-of-function evidence. Assign H if ANY of the following apply:
    * Constitutive activation, enhanced/hyperactive signaling, ligand-independent receptor activation
    * Increased protein stability, resistance to normal degradation pathways
    * Novel enzymatic activity, neomorphic biochemical function
    * C-terminal extension or truncation that creates new function (e.g. removes PEST/regulatory domain)
    * Hypermorphic allele — same function as WT but at higher level
    * NMD escape: truncating variant explicitly stated to escape nonsense-mediated decay,
      producing a stable truncated protein with demonstrated GOF activity
    * Resistance to proteasomal or autophagic degradation leading to protein accumulation
      with a demonstrated downstream gain-of-function consequence

- class "N": assign N if the mechanism is ANY of the following:
    * Toxic protein aggregation — mutant protein forms cytotoxic aggregates or inclusions
      (e.g. "toxic aggregation", "aggregate", "inclusion body", "amyloid-like")
    * Cytoplasmic mislocalization — mutant mislocalizes from nucleus/organelle to cytoplasm
      or vice versa, causing toxicity by sequestering WT or disrupting normal localization
      (e.g. "mislocalization", "cytoplasmic accumulation", "nuclear retention")
    * Toxic GOF — protein accumulates and is toxic but mechanism is indirect/non-specific
    * Dominant negative — mutant inhibits wild-type protein function
    * Uncertain / ambiguous — paper does not conclusively demonstrate GOF mechanism
    * Loss-of-function with secondary GOF effects not directly demonstrated
    * Frameshift/truncation where GOF is inferred but not experimentally confirmed

- IMPORTANT: Do NOT return empty findings just because the mechanism is toxic/indirect.
  Toxic aggregation, mislocalization, and dominant negative variants ARE valid findings — extract
  them and assign class "N". Only return {"findings": []} if the paper contains NO frameshift
  or stop-gain variants, or if the variants are clearly unrelated to GOF/dominant disease.

- When in doubt between H and N, assign N.
- evidence_text: direct quote from paper, max 200 chars

Examples:
[H] NOTCH2 c.7504C>T (p.Gln2502*): C-terminal truncation removes the PEST domain,
    preventing receptor degradation → constitutive activation of Notch signaling.
    → class "H" (direct constitutive activation, experimentally confirmed)

[H] CXCR4 c.959_960del (p.Val320GlufsTer23): C-terminal frameshift removes
    β-arrestin binding site → prolonged receptor signaling and ligand hypersensitivity.
    → class "H" (direct enhanced receptor signaling, confirmed in functional assay)

[H] ATM c.7271T>G (p.Val2424Gly*): paper explicitly states the truncated protein
    escapes NMD and accumulates, showing constitutive kinase activity.
    → class "H" (NMD escape + demonstrated GOF activity)

[H] OSBPL2 c.158_159del (p.Gln53ArgfsTer100): frameshift causes intracellular
    accumulation of mutant OSBPL2 that impairs autophagy flux with demonstrated
    downstream toxicity in transgenic mouse model.
    → class "H" is borderline; assign "H" only if autophagy impairment is shown as
    direct GOF consequence; assign "N" if mechanism is pure toxic accumulation.

[N] HSPB8 c.515dup (p.Pro173SerFsTer43): frameshift creates aggregation-prone
    C-terminal extension that sequesters wild-type chaperones into toxic inclusions.
    → class "N" (toxic aggregation / indirect cell death, not direct GOF)

[N] HNRNPA2B1 c.992del (p.Gly331GluFsTer28): frameshift reduces karyopherin β2
    affinity → cytoplasmic accumulation and recruitment of wild-type into stress
    granules, inhibiting normal RNA processing.
    → class "N" (dominant negative + toxic mislocalization, not direct GOF)

[N] CAV1 c.474del (p.Leu159SerFsTer22): frameshift introduces a de novo ER-retention
    motif; mutant sequesters in ER, disrupting caveolae formation. GOF is indirect
    (dominant negative over wild-type caveolae).
    → class "N" (dominant negative / indirect mechanism)

[N] OSBPL2 c.180_181del (p.His60GlnfsTer93): mutant OSBPL2 accumulates
    intracellularly and sequesters autophagy proteins into toxic inclusions, causing
    hearing loss. Mechanism is toxic accumulation, not direct GOF.
    → class "N" (toxic protein accumulation / indirect mechanism)
"""

_HGVS_SENT_RE = re.compile(
    r"[^.!?\n]{0,200}(?:c\.\d[\w_>*+\-]+|p\.[A-Za-z]{1,3}\d[\w*]+)[^.!?\n]{0,200}[.!?\n]",
    re.IGNORECASE,
)

def build_prompt(title: str, sections: dict) -> str:
    # Step 1: 모든 섹션에서 HGVS 포함 문장 수집 → [VARIANTS] 앞에 삽입
    hgvs_sents = []
    seen = set()
    for key in ["results", "methods", "discussion", "abstract", "introduction"]:
        for m in _HGVS_SENT_RE.finditer(sections.get(key, "")):
            s = m.group().strip()
            if s not in seen:
                hgvs_sents.append(s)
                seen.add(s)

    parts = []
    budget = 5000

    if hgvs_sents:
        snippet = "\n".join(hgvs_sents)[:1500]
        parts.append(f"[VARIANTS]\n{snippet}")
        budget -= len(snippet)

    # Step 2: 섹션 순서대로 나머지 budget 채우기
    for key in ["abstract", "results", "methods", "discussion"]:
        if key in sections and budget > 0:
            text = sections[key][:budget]
            parts.append(f"[{key.upper()}]\n{text}")
            budget -= len(text)

    return f"Title: {title}\n\n" + "\n\n".join(parts)


# ── vLLM 호출 ─────────────────────────────────────────────────────────────────
async def extract(
    session: aiohttp.ClientSession,
    record: dict,
    semaphore: asyncio.Semaphore,
    errors: list,
) -> dict | None:
    pmcid   = record["pmcid"]
    pmid    = record["pmid"]
    title   = record.get("title", "")
    sections= record.get("sections", {})

    user_msg = build_prompt(title, sections)
    payload  = {
        "model": MODEL,
        "messages": [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user",   "content": user_msg},
        ],
        "max_tokens": 1500,
        "temperature": 0.0,
        "response_format": {"type": "json_object"},
        "chat_template_kwargs": {"enable_thinking": False},
    }

    async with semaphore:
        for attempt in range(3):
            try:
                async with session.post(
                    f"{VLLM_BASE}/chat/completions",
                    json=payload,
                    timeout=aiohttp.ClientTimeout(total=120),
                ) as resp:
                    if resp.status != 200:
                        raise RuntimeError(f"HTTP {resp.status}: {(await resp.text())[:200]}")
                    data    = await resp.json()
                    content = data["choices"][0]["message"]["content"]
                    # think 블록 및 마크다운 코드블록 제거
                    content = re.sub(r"<think>.*?</think>", "", content, flags=re.DOTALL)
                    content = re.sub(r"```(?:json)?\s*|\s*```", "", content).strip()
                    parsed  = json.loads(content)
                    findings = parsed.get("findings", [])
                    # needs_review: 핵심 필드 누락 여부 판단
                    needs_review = any(
                        not f.get("gene") or not f.get("hgvs_coding")
                        for f in findings
                    ) if findings else False
                    return {
                        "pmid":         pmid,
                        "pmcid":        pmcid,
                        "doi":          record.get("doi", ""),
                        "title":        title,
                        "findings":     findings,
                        "n_findings":   len(findings),
                        "needs_review": needs_review,
                    }
            except Exception as e:
                if attempt == 2:
                    log.warning(f"실패 pmcid={pmcid}: {e}")
                    errors.append(pmcid)
                    return None
                await asyncio.sleep(2 ** attempt)
    return None


# ── 메인 ──────────────────────────────────────────────────────────────────────
async def run(args):
    # full-text 파일 목록
    ft_files = {f.stem: f for f in FT_DIR.glob("*.json")}

    # pmc_availability에서 pmid 매핑 (PMC*.json + PMID*.json 모두 커버)
    pmcid2pmid = {}
    with open(AVAIL_PATH) as f:
        for line in f:
            r = json.loads(line)
            pmid_val = str(r.get("pmid", ""))
            if r.get("pmcid"):
                pmcid2pmid[r["pmcid"]] = pmid_val
            if pmid_val:
                pmcid2pmid[f"PMID{pmid_val}"] = pmid_val

    # 대상 결정
    if args.test_gold:
        targets = [
            f for pmcid, f in ft_files.items()
            if pmcid2pmid.get(pmcid) in GOLD_PMIDS
        ]
        log.info(f"Gold test 모드: {len(targets)}/24개 full-text 있음")
    else:
        targets = list(ft_files.values())
        log.info(f"전체 모드: {len(targets):,}건")

    # 이미 처리된 것 스킵
    done_pmcids = set()
    if OUT_PATH.exists():
        with open(OUT_PATH) as f:
            for line in f:
                done_pmcids.add(json.loads(line)["pmcid"])
    targets = [f for f in targets if f.stem not in done_pmcids]
    log.info(f"신규 추출 대상: {len(targets)}건")

    if not targets:
        log.info("모두 처리 완료.")
        return

    # 서버 헬스 체크
    async with aiohttp.ClientSession() as sess:
        try:
            async with sess.get(f"{VLLM_BASE}/models",
                                timeout=aiohttp.ClientTimeout(total=10)) as r:
                if r.status != 200:
                    raise RuntimeError(f"HTTP {r.status}")
            log.info(f"vLLM 서버 OK: {VLLM_BASE}")
        except Exception as e:
            log.error(f"vLLM 서버 연결 실패: {e}")
            return

    semaphore = asyncio.Semaphore(args.concurrency)
    errors: list = []
    saved = 0

    connector = aiohttp.TCPConnector(limit=args.concurrency + 5)
    async with aiohttp.ClientSession(connector=connector) as session:
        with OUT_PATH.open("a", encoding="utf-8") as out_f, \
             NEEDS_REVIEW_PATH.open("a", encoding="utf-8") as rev_f:

            async def process_one(fp: Path):
                nonlocal saved
                record = json.loads(fp.read_text())
                # 파일 내 pmid 우선, 없으면 pmc_availability 매핑으로 보완
                pmid_in_file = str(record.get("pmid", ""))
                record["pmid"] = pmid_in_file or pmcid2pmid.get(fp.stem, "")
                result = await extract(session, record, semaphore, errors)
                if result is None:
                    return
                out_f.write(json.dumps(result, ensure_ascii=False) + "\n")
                out_f.flush()
                if result.get("needs_review"):
                    rev_f.write(json.dumps(result, ensure_ascii=False) + "\n")
                    rev_f.flush()
                saved += 1
                if args.test_gold:
                    n = result["n_findings"]
                    log.info(
                        f"  pmid={result['pmid']} | {n}개 finding | "
                        f"{result['title'][:55]}"
                    )
                elif saved % 1000 == 0:
                    pct = saved / len(targets) * 100
                    log.info(f"  [{saved:,}/{len(targets):,} {pct:.1f}%] 완료 | 오류: {len(errors)}")

            await asyncio.gather(*[process_one(f) for f in targets])

    nr_count = sum(1 for line in open(NEEDS_REVIEW_PATH) if line.strip()) if NEEDS_REVIEW_PATH.exists() else 0
    log.info(f"\n완료: {saved:,}건 저장 / {len(errors)}건 오류 / needs_review: {nr_count}건")
    log.info(f"  → {OUT_PATH}")
    log.info(f"  → needs_review: {NEEDS_REVIEW_PATH}")

    # 전체 실행이면 Mutalyzer normalization 자동 실행
    skip_mut = getattr(args, "skip_mutalyzer", False)
    if not args.test_gold and saved > 0 and not skip_mut:
        log.info("\nMutalyzer normalization 시작...")
        import subprocess, sys
        subprocess.run(
            [sys.executable, str(_HERE / "hgvs_normalize.py"), "--apply"],
            check=False
        )

    # gold 테스트 상세 출력
    if args.test_gold and OUT_PATH.exists():
        log.info("\n=== 추출 결과 상세 ===")
        with open(OUT_PATH) as f:
            for line in f:
                r = json.loads(line)
                if r["pmid"] not in GOLD_PMIDS:
                    continue
                log.info(f"\nPMID {r['pmid']}: {r['title'][:60]}")
                for i, f_ in enumerate(r["findings"], 1):
                    log.info(f"  [{i}] Gene={f_.get('gene')} | "
                             f"HGVS={f_.get('hgvs_coding')} | "
                             f"Class={f_.get('class')} | "
                             f"Effect={f_.get('gof_effect','')[:60]}")
                    log.info(f"       Evidence: {f_.get('evidence_text','')[:80]}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--test-gold", action="store_true",
                        help="gold 24건 테스트만 실행")
    parser.add_argument("--concurrency", type=int, default=4)
    parser.add_argument("--output", type=str, default=None,
                        help="출력 jsonl 경로 (기본: data/extracted_variants.jsonl)")
    parser.add_argument("--log", type=str, default=None,
                        help="로그 파일 경로 (기본: data/step2_extraction.log)")
    parser.add_argument("--skip-mutalyzer", action="store_true",
                        help="Mutalyzer normalization 건너뜀")
    args = parser.parse_args()

    if args.output:
        OUT_PATH = Path(args.output)
    if args.log:
        LOG_PATH = Path(args.log)
        logging.getLogger().handlers[1].baseFilename = str(LOG_PATH)  # FileHandler 교체
        for h in logging.getLogger().handlers:
            if isinstance(h, logging.FileHandler):
                logging.getLogger().removeHandler(h)
        logging.getLogger().addHandler(logging.FileHandler(LOG_PATH))

    asyncio.run(run(args))
