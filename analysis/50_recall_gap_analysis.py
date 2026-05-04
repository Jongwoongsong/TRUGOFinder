"""
50_recall_gap_analysis.py
역방향 recall 분석 — PubTator3/tmVar가 찾았는데 우리 모델이 놓친 변이

접근법:
  PubTator3 biocjson 캐시 (cache_pubtator/) → 각 PMID에서 tmVar가 인식한 Variant 목록
  우리 모델 출력 (gof_task_reclassified_v1.xlsx) → 각 PMID에서 추출된 변이 목록
  → 차집합 = "tmVar 인식 but 모델 미추출" 변이

우선순위 필터:
  Tier-A: frameshift/stop-gain/stop-loss 관련 표기 (*, X, Ter, fs, del, dup, ins)
           → 우리 타겟 variant type이므로 반드시 점검
  Tier-B: 나머지 변이 (missense 등) → 참고용

Output: data/recall_gap_analysis.xlsx
"""
import json
import re
from pathlib import Path
from collections import defaultdict

import pandas as pd

HERE = Path(__file__).resolve().parent.parent
CACHE_DIR  = HERE / "cache_pubtator"
TASK_XL    = HERE / "data/outputs/gof_task_reclassified_v1.xlsx"
OUT_XL     = HERE / "data/outputs/recall_gap_analysis.xlsx"

_ILLEGAL_RE = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f\ud800-\udfff]')
def _clean_df(d):
    return d.map(lambda v: _ILLEGAL_RE.sub("", v) if isinstance(v, str) else v)

# ── 1. 모델 출력 로드 ─────────────────────────────────────────────────────────
print("Loading model output...")
task_df = pd.read_excel(TASK_XL, sheet_name="All_2486", dtype=str).fillna("")

# PMID별 모델 추출 변이 set: (hgvs_coding, hgvs_protein, gene) 튜플
model_by_pmid = defaultdict(set)
for _, row in task_df.iterrows():
    pmid = row["pmid"].strip()
    if not pmid:
        continue
    model_by_pmid[pmid].add((
        row["hgvs_coding"].strip(),
        row["hgvs_protein"].strip(),
        row["gene"].strip(),
    ))

print(f"  Task PMIDs: {len(model_by_pmid)}")
print(f"  Total model variants: {len(task_df)}")

# ── 2. PubTator3 캐시 로드 ────────────────────────────────────────────────────
print("\nLoading PubTator3 cache...")
pt_by_pmid = defaultdict(list)   # pmid → list of annotation dicts

for f in sorted(CACHE_DIR.glob("*.json")):
    with open(f) as fp:
        data = json.load(fp)
    for pmid, annots in data.items():
        for a in annots:
            if a.get("annotation_type") == "Variant":
                pt_by_pmid[pmid].append(a)

print(f"  PMIDs with Variant annotations: {len(pt_by_pmid)}")
total_pt_vars = sum(len(v) for v in pt_by_pmid.values())
print(f"  Total tmVar Variant annotations: {total_pt_vars}")

# ── 3. 정규화 유틸 ────────────────────────────────────────────────────────────

# tmVar normalized_id 형식: "#gene_id#p.XXXX" 또는 "rs12345##" 또는 "#gene_id#c.XXXX"
def _parse_norm_id(norm_id: str) -> tuple[str, str]:
    """Return (coding, protein) extracted from PubTator normalized_id."""
    if not norm_id:
        return "", ""
    # 형식: #gene_id#p.XXX 또는 #gene_id#c.XXX
    m = re.search(r"#(c\.[^#\s]+)", str(norm_id))
    coding = m.group(1) if m else ""
    m2 = re.search(r"#(p\.[^#\s]+)", str(norm_id))
    protein = m2.group(1) if m2 else ""
    return coding, protein


def _extract_positions(text: str) -> set[str]:
    """Extract numeric positions from a variant mention (e.g. '123' from 'p.R123G')."""
    return set(re.findall(r"\d+", text))


def _is_ptc_relevant(mention: str, norm_id: str) -> bool:
    """True if mention looks like a frameshift / stop-gain / stop-loss / PTC variant."""
    combined = (mention + " " + str(norm_id)).lower()
    return bool(re.search(
        r"\bfs\b|frameshift|fsTer|fsX|"
        r"\bter\b|nonsense|\bstop\b|\btrunc|"
        r"[\*x]\d+|x\b|"           # e.g. R123X, p.*123
        r"\bdel\b|\bdup\b|\bins\b|"
        r"delins|ext\*|\bext[a-z]|"
        r"splice|splicing",
        combined,
    ))


def _model_covers(pt_mention: str, pt_norm_id: str, model_variants: set) -> bool:
    """
    Returns True if any model variant for this PMID 'matches' the PubTator annotation.
    Match logic (lenient):
      1. Normalized coding/protein from PubTator3 == model hgvs_coding/hgvs_protein
      2. Numeric positions overlap AND mention substring in model string
    """
    pt_coding, pt_protein = _parse_norm_id(pt_norm_id)
    pt_positions = _extract_positions(pt_mention)

    for m_coding, m_protein, m_gene in model_variants:
        # Exact HGVS match
        if pt_coding and pt_coding == m_coding:
            return True
        if pt_protein and pt_protein == m_protein:
            return True

        # Position-based fuzzy match: same numeric positions AND mention substring
        if pt_positions:
            m_positions = _extract_positions(m_coding + " " + m_protein)
            if pt_positions & m_positions:
                # At least one position overlaps → likely same variant
                return True

    return False


# ── 4. Gap 계산 ───────────────────────────────────────────────────────────────
print("\nComputing recall gap...")

gap_rows     = []   # PubTator found but model missed
covered_rows = []   # PubTator found and model covered

task_pmids = set(model_by_pmid.keys())

for pmid in sorted(task_pmids):
    pt_annots = pt_by_pmid.get(pmid, [])
    model_vars = model_by_pmid[pmid]

    for a in pt_annots:
        mention  = str(a.get("mention_text", ""))
        norm_id  = str(a.get("normalized_id", ""))
        context  = str(a.get("context_text", ""))
        passage  = str(a.get("passage_type", ""))

        is_ptc = _is_ptc_relevant(mention, norm_id)
        covered = _model_covers(mention, norm_id, model_vars)

        row = {
            "pmid":         pmid,
            "pt_mention":   mention,
            "pt_norm_id":   norm_id,
            "passage_type": passage,
            "is_ptc_relevant": is_ptc,
            "context_text": context[:300],
            "model_covered": covered,
            "model_variants_for_pmid": len(model_vars),
        }
        if covered:
            covered_rows.append(row)
        else:
            gap_rows.append(row)

print(f"  PubTator variants covered by model: {len(covered_rows)}")
print(f"  PubTator variants NOT in model (gap): {len(gap_rows)}")

gap_df     = pd.DataFrame(gap_rows)
covered_df = pd.DataFrame(covered_rows)

# ── 5. Gap 분석 ───────────────────────────────────────────────────────────────
ptc_gap = gap_df[gap_df["is_ptc_relevant"] == True].copy() if len(gap_df) > 0 else pd.DataFrame()
ptc_cov = covered_df[covered_df["is_ptc_relevant"] == True].copy() if len(covered_df) > 0 else pd.DataFrame()

print(f"\nPTC-relevant variants (frameshift/stop-gain/stop-loss):")
print(f"  Covered by model:  {len(ptc_cov)}")
print(f"  MISSED by model:   {len(ptc_gap)}")
if len(ptc_gap) + len(ptc_cov) > 0:
    recall = len(ptc_cov) / (len(ptc_cov) + len(ptc_gap)) * 100
    print(f"  Recall (PTC):      {recall:.1f}%")

# ── 6. 놓친 PTC 변이 — GOF 관련 컨텍스트 있는지 체크 ─────────────────────────
GOF_RE = re.compile(
    r"gain.of.function|constitutiv.{0,20}activ|dominant.negative|"
    r"toxic|aggregat|neomorphic|\bGOF\b|NMD.escape|hypermorph|"
    r"activating.mutation|gain.of.function",
    re.I,
)

if len(ptc_gap) > 0:
    ptc_gap["has_gof_context"] = ptc_gap["context_text"].apply(
        lambda t: bool(GOF_RE.search(t))
    )
    ptc_gap_gof = ptc_gap[ptc_gap["has_gof_context"] == True].copy()
    print(f"\nMissed PTC variants WITH GOF context: {len(ptc_gap_gof)}")
else:
    ptc_gap_gof = pd.DataFrame()
    print("\nNo missed PTC variants.")

# ── 7. 모델이 아예 추출하지 않은 PMID ────────────────────────────────────────
# PubTator에 PTC 변이가 있는데 모델 output이 0개인 경우
pmids_with_pt_ptc = set(ptc_gap["pmid"].unique()) if len(ptc_gap) > 0 else set()
pmids_model_zero  = {
    pmid for pmid in pmids_with_pt_ptc
    if len(model_by_pmid.get(pmid, set())) == 0
}
print(f"\nPMIDs where PubTator found PTC variants but model extracted 0 variants: {len(pmids_model_zero)}")

# ── 8. 통계 요약 ─────────────────────────────────────────────────────────────
pt_pmid_count = len([p for p in task_pmids if pt_by_pmid.get(p)])
pt_with_ptc   = len(set(ptc_gap["pmid"].tolist() + (ptc_cov["pmid"].tolist() if len(ptc_cov) > 0 else [])))

summary_rows = [
    ["=== PubTator3/tmVar vs 모델 Recall Gap 분석 ===", ""],
    ["", ""],
    ["[ Coverage ]", ""],
    ["Task dataset PMIDs",                  len(task_pmids)],
    ["PubTator 캐시에 있는 PMIDs",          pt_pmid_count],
    ["PubTator Variant 어노테이션 총계",    total_pt_vars],
    ["", ""],
    ["[ 전체 변이 (missense 포함) ]", ""],
    ["PubTator 커버 (모델 일치)",           len(covered_rows)],
    ["PubTator 갭 (모델 미추출)",           len(gap_rows)],
    ["", ""],
    ["[ PTC 관련 변이 (fs/stop/del/dup/ins) ]", ""],
    ["PTC: 모델이 커버",                    len(ptc_cov)],
    ["PTC: 모델이 놓침",                    len(ptc_gap)],
    ["PTC recall",                          f"{recall:.1f}%" if len(ptc_gap) + len(ptc_cov) > 0 else "N/A"],
    ["", ""],
    ["[ 놓친 PTC 중 GOF 컨텍스트 있음 ]",  ""],
    ["GOF 컨텍스트 있는 놓친 PTC 변이",    len(ptc_gap_gof)],
    ["해당 PMID 수",                        len(ptc_gap_gof["pmid"].nunique()) if len(ptc_gap_gof) > 0 else 0],
    ["", ""],
    ["[ 해석 주의사항 ]", ""],
    ["매칭 방법",  "위치 번호 기반 fuzzy + HGVS 직접 매칭 (완전 일치 아님)"],
    ["오버카운트", "PubTator는 동일 변이를 여러 패시지에서 중복 어노테이션 가능"],
    ["언더카운트", "모델은 fulltext 기반, PubTator는 abstract 중심 → 모델이 더 많이 찾을 수 있음"],
]

summary_df = pd.DataFrame(summary_rows, columns=["항목", "값"])

# ── 9. Excel 저장 ─────────────────────────────────────────────────────────────
with pd.ExcelWriter(OUT_XL, engine="openpyxl") as writer:
    _clean_df(summary_df).to_excel(      writer, sheet_name="Summary",            index=False)
    if len(ptc_gap_gof) > 0:
        _clean_df(ptc_gap_gof).to_excel( writer, sheet_name="Missed_PTC_GOF",     index=False)
    if len(ptc_gap) > 0:
        _clean_df(ptc_gap).to_excel(     writer, sheet_name="Missed_PTC_All",     index=False)
    if len(ptc_cov) > 0:
        _clean_df(ptc_cov).to_excel(     writer, sheet_name="Covered_PTC",        index=False)
    if len(gap_df) > 0:
        _clean_df(gap_df).to_excel(      writer, sheet_name="All_Gap",            index=False)

    for sn in writer.sheets:
        ws = writer.sheets[sn]
        for col_cells in ws.columns:
            ml = max((len(str(c.value)) if c.value else 0 for c in col_cells), default=0)
            ws.column_dimensions[col_cells[0].column_letter].width = min(ml + 2, 60)

print(f"\nSaved: {OUT_XL}")
if len(ptc_gap_gof) > 0:
    print(f"\n=== GOF 컨텍스트 있는 놓친 PTC 변이 TOP 30 ===")
    show_cols = ["pmid", "pt_mention", "pt_norm_id", "context_text"]
    print(_clean_df(ptc_gap_gof[show_cols].head(30)).to_string())
