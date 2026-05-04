"""
48c_stop_loss_merge.py
Stop-loss 증강 파이프라인 — Phase C: 추출 결과 병합 + 중복 제거

Input:
  data/extracted_variants_stoploss_v1.jsonl  (48b 추출 결과)
  data/gof_task_reclassified_v1.xlsx         (기존 task dataset, All_2486)

Pipeline:
  1. JSONL findings 정규화 (HGVS 형식 검사, 비인간 제거, is_human 필터)
  2. stop_loss 관련 variants 우선 추출 (ext* 패턴 포함 ALL variants도 보존)
  3. 기존 task dataset과 중복 제거 (gene+hgvs_coding+pmid key)
  4. primary_variant_class 등 reclassification 컬럼 추가 (47_task_reclassification.py 로직 재사용)
  5. 출력: data/gof_stop_loss_augmented_candidates.xlsx

Output:
  data/gof_stop_loss_augmented_candidates.xlsx
    - Sheet: Summary
    - Sheet: New_StopLoss      (stop_loss primary_variant_class만)
    - Sheet: New_AllVariants   (전체 신규 variants)
    - Sheet: Dedup_Skipped     (기존 dataset과 중복이라 제외된 것)
"""
import json
import re
from pathlib import Path
from collections import Counter

import pandas as pd

HERE = Path(__file__).resolve().parent.parent
SL_JSONL = HERE / "data/extracted_variants_stoploss_v1.jsonl"
TASK_XL  = HERE / "data/outputs/gof_task_reclassified_v1.xlsx"
OUT_XL   = HERE / "data/outputs/gof_stop_loss_augmented_candidates.xlsx"

_ILLEGAL_RE = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f\ud800-\udfff]')

def _clean_df(d: pd.DataFrame) -> pd.DataFrame:
    return d.map(lambda v: _ILLEGAL_RE.sub("", v) if isinstance(v, str) else v)


# ── 1. 기존 task dataset 로드 ─────────────────────────────────────────────────
print("Loading existing task dataset...")
task_df = pd.read_excel(TASK_XL, sheet_name="All_2486", dtype=str).fillna("")
print(f"  Existing rows: {len(task_df)}")

# Dedup key: gene + hgvs_coding + pmid (same as existing pipeline)
existing_keys = set(
    zip(task_df["gene"], task_df["hgvs_coding"], task_df["pmid"])
)
print(f"  Existing unique (gene, hgvs_coding, pmid): {len(existing_keys)}")

# ── 2. 추출 JSONL 로드 및 정규화 ─────────────────────────────────────────────
print(f"\nLoading extraction results: {SL_JSONL}")
if not SL_JSONL.exists():
    print("  ERROR: JSONL file not found. Run 48b first.")
    raise SystemExit(1)

HGVS_CODING_RE = re.compile(r"^c\.[0-9_+\-*]", re.I)
HGVS_PROT_RE   = re.compile(r"^p\.", re.I)
EXT_RE          = re.compile(r"ext[\*\dA-Za-z]|extTer", re.I)

raw_findings = []
paper_count  = 0
error_count  = 0

with open(SL_JSONL) as f:
    for line in f:
        doc = json.loads(line)
        if doc.get("error"):
            error_count += 1
            continue
        paper_count += 1
        pmid   = str(doc.get("pmid", "")).strip()
        pmcid  = str(doc.get("pmcid", "")).strip()
        title  = str(doc.get("title", "")).strip()
        source = str(doc.get("source", "stoploss_extraction")).strip()

        for fnd in doc.get("findings", []):
            gene    = str(fnd.get("gene", "")).strip()
            disease = str(fnd.get("disease", "")).strip()
            omim    = str(fnd.get("omim", "")).strip()
            hgvs_c  = str(fnd.get("hgvs_coding", "")).strip()
            hgvs_p  = str(fnd.get("hgvs_protein", "")).strip()
            gof_eff = str(fnd.get("gof_effect", "")).strip()
            cls     = str(fnd.get("class", "N")).strip().upper()
            is_hum  = fnd.get("is_human", True)
            note    = str(fnd.get("note", "")).strip()
            ev_text = str(fnd.get("evidence_text", "")).strip()

            # Basic filters
            if not gene or len(gene) < 2:
                continue
            if not is_hum:
                continue  # exclude non-human
            if cls not in ("H", "N"):
                cls = "N"

            # HGVS format guard
            if hgvs_c and not HGVS_CODING_RE.match(hgvs_c):
                hgvs_c = ""
            if hgvs_p and not HGVS_PROT_RE.match(hgvs_p):
                hgvs_p = ""

            raw_findings.append({
                "gene":         gene,
                "disease":      disease,
                "omim":         omim,
                "hgvs_coding":  hgvs_c,
                "hgvs_protein": hgvs_p,
                "gof_effect":   gof_eff,
                "class":        cls,
                "is_human":     True,
                "note":         note,
                "evidence_text": ev_text,
                "pmid":         pmid,
                "pmcid":        pmcid,
                "title":        title,
                "source":       source,
            })

print(f"  Papers with findings: {paper_count}")
print(f"  Papers with errors:   {error_count}")
print(f"  Raw findings:         {len(raw_findings)}")

# ── 3. 중복 제거 (기존 dataset 대비) ─────────────────────────────────────────
new_rows  = []
dup_rows  = []

for row in raw_findings:
    key = (row["gene"], row["hgvs_coding"], row["pmid"])
    if key in existing_keys:
        dup_rows.append(row)
    else:
        new_rows.append(row)
        existing_keys.add(key)   # prevent intra-batch duplicates too

print(f"\nAfter dedup:")
print(f"  New variants:     {len(new_rows)}")
print(f"  Already in task:  {len(dup_rows)}")

# ── 4. primary_variant_class 할당 (47_task_reclassification.py 로직) ────────
NMD_FS_RE = re.compile(r"fs|fsTer|fsX|fs\*",        re.I)
NMD_SG_RE = re.compile(r"Ter\d|[\*X]\d|Ter$|[\*X]$", re.I)

def _primary_variant_class(hgvs_c: str, hgvs_p: str) -> str:
    """Classify based on HGVS notation."""
    hc = hgvs_c.lower()
    hp = hgvs_p

    # Stop-loss: extension beyond stop codon
    if EXT_RE.search(hp) and ("ext" in hp.lower()):
        return "STOP_LOSS"
    # Frameshift
    if "fs" in hp.lower() or "del" in hc or "dup" in hc or "ins" in hc:
        if not EXT_RE.search(hp):
            return "FRAMESHIFT"
    # Stop-gain
    if re.search(r"[\*X]$|Ter$|Ter\d", hp):
        return "STOP_GAIN"
    # Fallback heuristics
    if re.search(r"del|dup|ins|fs", hc):
        return "FRAMESHIFT"
    if re.search(r"[\*]", hc):
        return "STOP_GAIN"
    return "OTHER"


def _stop_loss_subtype(hgvs_p: str, pvc: str) -> str:
    if pvc != "STOP_LOSS":
        return "NOT_APPLICABLE"
    if EXT_RE.search(hgvs_p):
        return "GENETIC_STOP_CODON_LOSS"
    return "NONSTOP_EXTENSION"


new_df = pd.DataFrame(new_rows)
if len(new_df) == 0:
    print("\nNo new variants found.")
else:
    new_df["primary_variant_class"] = new_df.apply(
        lambda r: _primary_variant_class(r["hgvs_coding"], r["hgvs_protein"]), axis=1
    )
    new_df["stop_loss_subtype"] = new_df.apply(
        lambda r: _stop_loss_subtype(r["hgvs_protein"], r["primary_variant_class"]), axis=1
    )
    new_df["augmentation_source"] = "stoploss_extraction_v1"

    pvc_counts = new_df["primary_variant_class"].value_counts()
    print(f"\nprimary_variant_class distribution:")
    for k, v in pvc_counts.items():
        print(f"  {k}: {v}")

    sl_df  = new_df[new_df["primary_variant_class"] == "STOP_LOSS"].copy()
    print(f"\nNew stop-loss variants: {len(sl_df)}")
    if len(sl_df) > 0:
        print(sl_df[["gene","hgvs_coding","hgvs_protein","disease","gof_effect"]].to_string())

dup_df = pd.DataFrame(dup_rows)

# ── 5. Summary ─────────────────────────────────────────────────────────────
summary_rows = [
    ["=== Stop-Loss Augmentation Candidates ===", ""],
    ["", ""],
    ["[ Extraction stats ]", ""],
    ["Papers processed",    paper_count],
    ["Papers with errors",  error_count],
    ["Raw findings",        len(raw_findings)],
    ["", ""],
    ["[ Deduplication ]", ""],
    ["Already in task dataset",  len(dup_rows)],
    ["New variants",              len(new_rows)],
    ["", ""],
]
if len(new_df) > 0:
    summary_rows += [["[ primary_variant_class (new) ]", ""]]
    for k, v in new_df["primary_variant_class"].value_counts().items():
        summary_rows.append([f"  {k}", int(v)])
    sl_count = len(new_df[new_df["primary_variant_class"] == "STOP_LOSS"])
    summary_rows += [
        ["", ""],
        ["[ Stop-loss new variants ]", sl_count],
    ]

summary_df = pd.DataFrame(summary_rows, columns=["항목", "값"])

# ── 6. Excel 저장 ─────────────────────────────────────────────────────────────
with pd.ExcelWriter(OUT_XL, engine="openpyxl") as writer:
    _clean_df(summary_df).to_excel(writer, sheet_name="Summary",         index=False)
    if len(new_df) > 0:
        sl_new = new_df[new_df["primary_variant_class"] == "STOP_LOSS"].copy()
        _clean_df(new_df).to_excel(writer, sheet_name="New_AllVariants", index=False)
        _clean_df(sl_new).to_excel(writer, sheet_name="New_StopLoss",    index=False)
    if len(dup_df) > 0:
        _clean_df(dup_df).to_excel(writer, sheet_name="Dedup_Skipped",   index=False)

    for sn in writer.sheets:
        ws = writer.sheets[sn]
        for col_cells in ws.columns:
            ml = max((len(str(c.value)) if c.value else 0 for c in col_cells), default=0)
            ws.column_dimensions[col_cells[0].column_letter].width = min(ml + 2, 60)

print(f"\nSaved: {OUT_XL}")
if len(new_df) > 0:
    print(f"  New_AllVariants: {len(new_df)} rows")
    print(f"  New_StopLoss:    {len(sl_df)} rows")
    print(f"  Dedup_Skipped:   {len(dup_df)} rows")
