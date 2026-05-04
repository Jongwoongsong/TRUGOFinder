"""
49_export_task_split.py
GOF Pipeline v10 최종 export — task-split sheets

Input:
  data/gof_task_reclassified_v1.xlsx          (기존 2,486 rows + 8 reclassification 컬럼)
  data/gof_stop_loss_augmented_candidates.xlsx (48c 출력 — 신규 stop-loss 등 추가 variants)

Pipeline:
  1. 기존 2,486 rows 로드
  2. 신규 augmentation 후보 중 include_in_core 판단 후 append
  3. 최종 dedup
  4. Task-split sheets 생성

Output: data/gof_pipeline_final_v10_task_split.xlsx
  - Summary
  - All_variants         (기존 + 신규 전체)
  - FRAMESHIFT_GOF_core  (include_in_core=YES, primary_gof_task=FRAMESHIFT_GOF)
  - STOP_GAIN_GOF_core   (include_in_core=YES, primary_gof_task=STOP_GAIN_GOF)
  - STOP_LOSS_GOF_core   (include_in_core=YES, primary_gof_task=STOP_LOSS_GOF)
  - Core_REVIEW          (include_in_core=REVIEW)
  - Core_NO              (include_in_core=NO)
  - New_augmented        (신규 추가된 rows only)
"""
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent.parent
TASK_XL   = HERE / "data/outputs/gof_task_reclassified_v1.xlsx"
AUG_XL    = HERE / "data/outputs/gof_stop_loss_augmented_candidates.xlsx"
OUT_XL    = HERE / "data/outputs/gof_pipeline_final_v10_task_split.xlsx"

_ILLEGAL_RE = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f\ud800-\udfff]')

def _clean_df(d: pd.DataFrame) -> pd.DataFrame:
    return d.map(lambda v: _ILLEGAL_RE.sub("", v) if isinstance(v, str) else v)


# ── 1. 기존 task dataset 로드 ─────────────────────────────────────────────────
print("Loading existing task dataset...")
base_df = pd.read_excel(TASK_XL, sheet_name="All_2486", dtype=str).fillna("")
print(f"  Base rows: {len(base_df)}")

# ── 2. 신규 augmentation 로드 ─────────────────────────────────────────────────
print("Loading augmentation candidates...")
if not AUG_XL.exists():
    print("  WARNING: augmentation file not found — using base dataset only")
    aug_new = pd.DataFrame()
else:
    try:
        aug_new = pd.read_excel(AUG_XL, sheet_name="New_AllVariants", dtype=str).fillna("")
        print(f"  Augmentation new variants: {len(aug_new)}")
    except Exception as e:
        print(f"  WARNING: could not load augmentation ({e}) — using base only")
        aug_new = pd.DataFrame()

# ── 3. 신규 variants reclassification 컬럼 추가 (없는 경우 기본값) ─────────────
NEW_TASK_COLS = [
    "primary_variant_class", "stop_loss_subtype", "primary_gof_task",
    "gof_due_to_variant", "gof_evidence_strength", "nmd_escape_status",
    "hn_call", "include_in_core_task_dataset", "task_reclassification_reason",
]

if len(aug_new) > 0:
    # primary_variant_class is already set by 48c
    # Fill remaining reclassification columns with defaults
    if "primary_gof_task" not in aug_new.columns:
        def _task(pvc):
            m = {"FRAMESHIFT": "FRAMESHIFT_GOF", "STOP_GAIN": "STOP_GAIN_GOF",
                 "STOP_LOSS": "STOP_LOSS_GOF", "QC_FAIL": "REVIEW", "OTHER": "REVIEW"}
            return m.get(pvc, "REVIEW")
        aug_new["primary_gof_task"] = aug_new["primary_variant_class"].map(_task)

    if "gof_due_to_variant" not in aug_new.columns:
        cls_map = {"H": "YES", "N": "REVIEW"}
        aug_new["gof_due_to_variant"] = aug_new.get("class", pd.Series("N", index=aug_new.index)).map(cls_map).fillna("REVIEW")

    if "gof_evidence_strength" not in aug_new.columns:
        aug_new["gof_evidence_strength"] = "STRONG_ASSERTION"  # from literature

    if "nmd_escape_status" not in aug_new.columns:
        def _nmd(pvc):
            return "NOT_APPLICABLE" if pvc == "STOP_LOSS" else "UNKNOWN"
        aug_new["nmd_escape_status"] = aug_new["primary_variant_class"].map(_nmd)

    if "hn_call" not in aug_new.columns:
        cls_map2 = {"H": "HYPERMORPH", "N": "REVIEW"}
        aug_new["hn_call"] = aug_new.get("class", pd.Series("N", index=aug_new.index)).map(cls_map2).fillna("REVIEW")

    if "include_in_core_task_dataset" not in aug_new.columns:
        def _include(pvc):
            return "YES" if pvc in ("FRAMESHIFT", "STOP_GAIN", "STOP_LOSS") else "REVIEW"
        aug_new["include_in_core_task_dataset"] = aug_new["primary_variant_class"].map(_include)

    if "task_reclassification_reason" not in aug_new.columns:
        aug_new["task_reclassification_reason"] = "augmented_stoploss_extraction_v1"

    if "stop_loss_subtype" not in aug_new.columns:
        aug_new["stop_loss_subtype"] = "NOT_APPLICABLE"

    # Mark as augmented
    aug_new["augmentation_batch"] = "stoploss_v1"

    # Align columns with base_df (add missing cols as empty)
    for col in base_df.columns:
        if col not in aug_new.columns:
            aug_new[col] = ""
    aug_new = aug_new[[c for c in base_df.columns if c in aug_new.columns] +
                      [c for c in aug_new.columns if c not in base_df.columns]]
else:
    aug_new = pd.DataFrame(columns=base_df.columns)

# ── 4. Combine + final dedup ──────────────────────────────────────────────────
base_df["augmentation_batch"] = "original"

combined = pd.concat([base_df, aug_new], ignore_index=True)
n_before_dedup = len(combined)
combined = combined.drop_duplicates(
    subset=["gene", "hgvs_coding", "hgvs_protein", "pmid"], keep="first"
)
n_after_dedup = len(combined)
print(f"\nCombined: {n_before_dedup} → {n_after_dedup} (dedup removed {n_before_dedup - n_after_dedup})")

# ── 5. Task-split subsets ─────────────────────────────────────────────────────
inc_col  = "include_in_core_task_dataset"
task_col = "primary_gof_task"

df_yes    = combined[combined[inc_col] == "YES"].copy()
df_review = combined[combined[inc_col] == "REVIEW"].copy()
df_no     = combined[combined[inc_col] == "NO"].copy()
df_new    = combined[combined.get("augmentation_batch", pd.Series("original", index=combined.index)) == "stoploss_v1"].copy() if len(aug_new) > 0 else pd.DataFrame()

df_fs = df_yes[df_yes[task_col] == "FRAMESHIFT_GOF"].copy()
df_sg = df_yes[df_yes[task_col] == "STOP_GAIN_GOF"].copy()
df_sl = df_yes[df_yes[task_col] == "STOP_LOSS_GOF"].copy()

print(f"\n=== v10 Final Summary ===")
print(f"Total rows:         {len(combined)}")
print(f"  Original base:    {len(base_df)}")
print(f"  New augmented:    {len(aug_new)}")
print(f"\nCore YES:           {len(df_yes)}")
print(f"  FRAMESHIFT_GOF:   {len(df_fs)}")
print(f"  STOP_GAIN_GOF:    {len(df_sg)}")
print(f"  STOP_LOSS_GOF:    {len(df_sl)}")
print(f"Core REVIEW:        {len(df_review)}")
print(f"Core NO:            {len(df_no)}")

if len(df_sl) > 0:
    print(f"\nStop-loss core YES variants:")
    print(df_sl[["gene","hgvs_coding","hgvs_protein","disease","hn_call"]].to_string())

# ── 6. Summary sheet rows ─────────────────────────────────────────────────────
sl_hn = df_sl["hn_call"].value_counts().to_dict() if len(df_sl) > 0 else {}

summary_rows = [
    ["=== GOF Pipeline v10 — Task Split Final ===", ""],
    ["", ""],
    ["[ Dataset composition ]", ""],
    ["Original base (v9 reclassified)", len(base_df)],
    ["New augmented (stop-loss v1)",    len(aug_new)],
    ["After dedup",                     n_after_dedup],
    ["", ""],
    ["[ include_in_core_task_dataset ]", ""],
    ["YES",    len(df_yes)],
    ["REVIEW", len(df_review)],
    ["NO",     len(df_no)],
    ["", ""],
    ["[ Core task split (YES) ]", ""],
    ["FRAMESHIFT_GOF", len(df_fs)],
    ["STOP_GAIN_GOF",  len(df_sg)],
    ["STOP_LOSS_GOF",  len(df_sl)],
    ["", ""],
    ["[ hn_call distribution (YES) ]", ""],
]
for k, v in df_yes["hn_call"].value_counts().items():
    summary_rows.append([f"  {k}", int(v)])
summary_rows += [
    ["", ""],
    ["[ Stop-loss GOF hn_call ]", ""],
]
for k, v in sl_hn.items():
    summary_rows.append([f"  {k}", int(v)])

summary_df = pd.DataFrame(summary_rows, columns=["항목", "값"])

# ── 7. Excel export ──────────────────────────────────────────────────────────
print(f"\nWriting {OUT_XL}...")
with pd.ExcelWriter(OUT_XL, engine="openpyxl") as writer:
    _clean_df(summary_df).to_excel(writer, sheet_name="Summary",            index=False)
    _clean_df(combined).to_excel(  writer, sheet_name="All_variants",       index=False)
    _clean_df(df_fs).to_excel(     writer, sheet_name="FRAMESHIFT_GOF_core",index=False)
    _clean_df(df_sg).to_excel(     writer, sheet_name="STOP_GAIN_GOF_core", index=False)
    _clean_df(df_sl).to_excel(     writer, sheet_name="STOP_LOSS_GOF_core", index=False)
    _clean_df(df_review).to_excel( writer, sheet_name="Core_REVIEW",        index=False)
    _clean_df(df_no).to_excel(     writer, sheet_name="Core_NO",            index=False)
    if len(df_new) > 0:
        _clean_df(df_new).to_excel(writer, sheet_name="New_augmented",      index=False)

    for sn in writer.sheets:
        ws = writer.sheets[sn]
        for col_cells in ws.columns:
            ml = max((len(str(c.value)) if c.value else 0 for c in col_cells), default=0)
            ws.column_dimensions[col_cells[0].column_letter].width = min(ml + 2, 60)

print(f"Saved: {OUT_XL}")
print(f"  All_variants:          {len(combined)}")
print(f"  FRAMESHIFT_GOF_core:   {len(df_fs)}")
print(f"  STOP_GAIN_GOF_core:    {len(df_sg)}")
print(f"  STOP_LOSS_GOF_core:    {len(df_sl)}")
