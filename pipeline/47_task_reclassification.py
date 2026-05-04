"""
47_task_reclassification.py
기존 2,486 rows에 task-specific reclassification 컬럼 추가

새 컬럼 9개:
  primary_variant_class    : FRAMESHIFT / STOP_GAIN / STOP_LOSS / OTHER / QC_FAIL
  stop_loss_subtype        : GENETIC_STOP_CODON_LOSS / READTHROUGH_ASSOCIATED /
                             NONSTOP_EXTENSION / UNKNOWN / NOT_APPLICABLE
  primary_gof_task         : FRAMESHIFT_GOF / STOP_GAIN_GOF / STOP_LOSS_GOF /
                             OUT_OF_SCOPE_BROAD_GOF / NON_GOF_EXCLUDE / REVIEW
  gof_due_to_variant       : YES / REVIEW / NO / CONFLICTING
  gof_evidence_strength    : DIRECT_FUNCTIONAL / STRONG_ASSERTION / INFERRED_ONLY /
                             NO_GOF_SHOWN / CONTRADICTS_GOF / UNCERTAIN
  nmd_escape_status        : NMD_ESCAPE_REPORTED / NMD_ESCAPE_INFERRED /
                             NMD_SENSITIVE_OR_DEGRADED / UNKNOWN / NOT_APPLICABLE
  hn_call                  : HYPERMORPH / NEOMORPH / TOXIC_GOF_SEPARATE /
                             DOMINANT_NEGATIVE_SEPARATE / NMD_ESCAPE_GOF_REVIEW /
                             NON_GOF / REVIEW
  include_in_core_task_dataset : YES / REVIEW / NO
  task_reclassification_reason : free text

Input:  data/gof_curated_after_manual_qc_round1.xlsx (Curated_All_With_QC_Action)
Output: data/gof_task_reclassified_v1.xlsx
"""
import re
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
IN_XL  = HERE / "data/outputs/gof_curated_after_manual_qc_round1.xlsx"
OUT_XL = HERE / "data/outputs/gof_task_reclassified_v1.xlsx"

df = pd.read_excel(IN_XL, sheet_name="Curated_All_With_QC_Action", dtype=str).fillna("")
print(f"Loaded: {len(df)} rows, {len(df.columns)} columns")

# ── helpers ───────────────────────────────────────────────────────────────────
_ILLEGAL_RE  = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f\ud800-\udfff]')
NMD_FS_RE    = re.compile(r"fs|fsTer|fsX|fs\*",     re.I)
NMD_SG_RE    = re.compile(r"Ter\d|[\*X]\d|Ter$|[\*X]$", re.I)
NMD_TEXT_RE  = re.compile(r"\bNMD\b|nonsense.mediated.decay|nonsense-mediated", re.I)
NEOMORPH_RE  = re.compile(r"neomorp",               re.I)
EXT_RE       = re.compile(r"ext[\*\dA-Za-z]|extTer", re.I)

NON_GOF_MECHS = {"LOSS_OF_FUNCTION_NON_GOF"}
OOS_MECHS     = {"LOF_MEDIATED_PATHWAY_ACTIVATION"}


def _clean_df(d: pd.DataFrame) -> pd.DataFrame:
    def _c(v):
        return _ILLEGAL_RE.sub("", v) if isinstance(v, str) else v
    return d.map(_c)


# ── 1. primary_variant_class ──────────────────────────────────────────────────
def _primary_variant_class(row):
    fvt  = row.get("final_variant_type", "")
    mech = row.get("curated_mechanism_class", "")
    prot = row.get("hgvs_protein", "")

    if mech == "QC_FAIL_VARIANT_ID":
        return "QC_FAIL"
    if fvt == "frameshift":
        return "FRAMESHIFT"
    if fvt == "stop_gain":
        return "STOP_GAIN"
    if fvt == "stop_loss":
        return "STOP_LOSS"
    if fvt == "NMD_escape":
        # Remap to underlying variant class based on protein annotation.
        # fs/fsTer → FRAMESHIFT; Ter/*/X (without fs) → STOP_GAIN
        if NMD_FS_RE.search(prot):
            return "FRAMESHIFT"
        if NMD_SG_RE.search(prot):
            return "STOP_GAIN"
        return "FRAMESHIFT"  # most NMD_escapes originate from frameshifts
    return "OTHER"


df["primary_variant_class"] = df.apply(_primary_variant_class, axis=1)
print("primary_variant_class:", df["primary_variant_class"].value_counts().to_dict())

# ── 2. stop_loss_subtype ──────────────────────────────────────────────────────
# All 9 stop_loss rows in this dataset carry stop-codon mutations that create
# an extension peptide (hgvs_protein contains "ext").  None are readthrough-
# associated (translational suppression without a DNA-level stop-codon change).
def _stop_loss_subtype(row):
    if row["primary_variant_class"] != "STOP_LOSS":
        return "NOT_APPLICABLE"
    prot = row.get("hgvs_protein", "")
    if EXT_RE.search(prot) or "ext" in prot.lower():
        return "GENETIC_STOP_CODON_LOSS"   # mutation IN stop codon → extension
    # fallback for any future stop_loss without clear ext annotation
    return "NONSTOP_EXTENSION"


df["stop_loss_subtype"] = df.apply(_stop_loss_subtype, axis=1)

# ── 3. primary_gof_task ───────────────────────────────────────────────────────
# Assignment is mechanism/variant-based; EXCLUDE rows → NON_GOF_EXCLUDE.
# LOF-mediated pathway activation is out of scope: the variant is LOF; the
# pathway-level GOF effect is secondary and not a direct variant property.

TASK_MAP = {
    "FRAMESHIFT": "FRAMESHIFT_GOF",
    "STOP_GAIN":  "STOP_GAIN_GOF",
    "STOP_LOSS":  "STOP_LOSS_GOF",
    "QC_FAIL":    "REVIEW",
    "OTHER":      "OUT_OF_SCOPE_BROAD_GOF",
}


def _primary_gof_task(row):
    pvc  = row["primary_variant_class"]
    mech = row.get("curated_mechanism_class", "")
    qc   = row.get("qc_round1_action", "")

    # Manually excluded rows → not GOF for our purposes
    if qc.startswith("EXCLUDE"):
        return "NON_GOF_EXCLUDE"

    # Clear non-GOF mechanism
    if mech in NON_GOF_MECHS:
        return "NON_GOF_EXCLUDE"

    # Variant is LOF; pathway activation is an indirect secondary effect
    if mech in OOS_MECHS:
        return "OUT_OF_SCOPE_BROAD_GOF"

    return TASK_MAP.get(pvc, "REVIEW")


df["primary_gof_task"] = df.apply(_primary_gof_task, axis=1)
print("primary_gof_task:", df["primary_gof_task"].value_counts().to_dict())

# ── 4. gof_due_to_variant ─────────────────────────────────────────────────────
def _gof_due_to_variant(row):
    mech = row.get("curated_mechanism_class", "")
    qc   = row.get("qc_round1_action", "")

    # Excluded rows: no valid GOF evidence or proven non-GOF
    if qc in ("EXCLUDE_NO_GOF_EVIDENCE", "EXCLUDE_NON_GOF_LOF",
              "EXCLUDE_NONHUMAN_EVIDENCE", "EXCLUDE_GENE_EVIDENCE_MISMATCH"):
        return "NO"
    if mech in ("LOSS_OF_FUNCTION_NON_GOF", "QC_FAIL_VARIANT_ID"):
        return "NO"

    # Conflicting signals (LOF and GOF keywords both present)
    if mech == "UNCERTAIN_OR_CONFLICTING" or qc == "REVIEW_MECHANISM":
        return "CONFLICTING"

    # Clear GOF mechanisms
    if mech in ("HYPERMORPH_TRUE_GOF", "DOMINANT_NEGATIVE_ANTIMORPH", "TOXIC_GOF"):
        return "YES"

    # NMD-escape: strong presumption of GOF but review recommended
    # LOF-mediated pathway: GOF is indirect, needs review
    if mech in ("NMD_ESCAPE_EFFECT_REVIEW", "LOF_MEDIATED_PATHWAY_ACTIVATION"):
        return "REVIEW"

    if qc.startswith("REVIEW"):
        return "REVIEW"

    return "REVIEW"


df["gof_due_to_variant"] = df.apply(_gof_due_to_variant, axis=1)
print("gof_due_to_variant:", df["gof_due_to_variant"].value_counts().to_dict())

# ── 5. gof_evidence_strength ──────────────────────────────────────────────────
EVSTR_MAP = {
    "direct_functional":     "DIRECT_FUNCTIONAL",
    "literature_claim_only": "STRONG_ASSERTION",
    "inferred_from_HGVS":    "INFERRED_ONLY",
    "indirect_pathway":      "INFERRED_ONLY",
    "insufficient":          "UNCERTAIN",
    "qc_fail":               "NO_GOF_SHOWN",
}


def _gof_evidence_strength(row):
    es   = row.get("evidence_strength", "")
    mech = row.get("curated_mechanism_class", "")
    qc   = row.get("qc_round1_action", "")

    # Proven non-GOF → CONTRADICTS
    if mech == "LOSS_OF_FUNCTION_NON_GOF":
        return "CONTRADICTS_GOF"

    # Explicitly excluded for lack of GOF evidence
    if qc == "EXCLUDE_NO_GOF_EVIDENCE":
        return "NO_GOF_SHOWN"

    return EVSTR_MAP.get(es, "UNCERTAIN")


df["gof_evidence_strength"] = df.apply(_gof_evidence_strength, axis=1)
print("gof_evidence_strength:", df["gof_evidence_strength"].value_counts().to_dict())

# ── 6. nmd_escape_status ──────────────────────────────────────────────────────
# NMD_escape is a secondary property of frameshift/stop_gain variants.
# Stop-loss variants are not subject to NMD (no premature stop → no degradation).
def _nmd_escape_status(row):
    fvt  = row.get("final_variant_type", "")
    mech = row.get("curated_mechanism_class", "")
    pvc  = row["primary_variant_class"]

    if pvc == "STOP_LOSS":       # extension proteins, NMD not relevant
        return "NOT_APPLICABLE"
    if mech == "QC_FAIL_VARIANT_ID":
        return "NOT_APPLICABLE"

    # Original pipeline explicitly classified as NMD_escape
    if fvt == "NMD_escape":
        return "NMD_ESCAPE_REPORTED"

    # LLM curated the mechanism as NMD-escape even though the original
    # pipeline didn't flag the variant as NMD_escape
    if mech == "NMD_ESCAPE_EFFECT_REVIEW":
        return "NMD_ESCAPE_INFERRED"

    # Non-GOF LOF: the truncated/frameshifted protein is degraded via NMD
    if mech == "LOSS_OF_FUNCTION_NON_GOF":
        return "NMD_SENSITIVE_OR_DEGRADED"

    # Check free-text for NMD keywords even when mechanism is not NMD_ESCAPE
    combined = " ".join([row.get("gof_effect", ""), row.get("evidence_text", "")])
    if NMD_TEXT_RE.search(combined):
        return "NMD_ESCAPE_INFERRED"

    return "UNKNOWN"


df["nmd_escape_status"] = df.apply(_nmd_escape_status, axis=1)
print("nmd_escape_status:", df["nmd_escape_status"].value_counts().to_dict())

# ── 7. hn_call ────────────────────────────────────────────────────────────────
# H/N not forced on all GOF; DN and Toxic remain separate categories.
# NEOMORPH is assigned when gof_effect explicitly uses "neomorphic" language AND
# the curated mechanism is HYPERMORPH_TRUE_GOF (meaning the original H/N call
# labelled it H but the mechanism text reveals true neomorphic character).

MECH_HN_MAP = {
    "HYPERMORPH_TRUE_GOF":          "HYPERMORPH",              # may → NEOMORPH below
    "DOMINANT_NEGATIVE_ANTIMORPH":  "DOMINANT_NEGATIVE_SEPARATE",
    "TOXIC_GOF":                    "TOXIC_GOF_SEPARATE",
    "NMD_ESCAPE_EFFECT_REVIEW":     "NMD_ESCAPE_GOF_REVIEW",
    "LOSS_OF_FUNCTION_NON_GOF":     "NON_GOF",
    "LOF_MEDIATED_PATHWAY_ACTIVATION": "NON_GOF",
    "QC_FAIL_VARIANT_ID":           "REVIEW",
    "UNCERTAIN_OR_CONFLICTING":     "REVIEW",
}


def _hn_call(row):
    mech   = row.get("curated_mechanism_class", "")
    qc     = row.get("qc_round1_action", "")
    effect = row.get("gof_effect", "")

    if qc.startswith("EXCLUDE"):
        return "NON_GOF"

    hn = MECH_HN_MAP.get(mech, "REVIEW")

    # Neomorphic override: gof_effect explicitly describes a new/novel function
    if hn == "HYPERMORPH" and NEOMORPH_RE.search(effect):
        return "NEOMORPH"

    return hn


df["hn_call"] = df.apply(_hn_call, axis=1)
print("hn_call:", df["hn_call"].value_counts().to_dict())

# ── 8. include_in_core_task_dataset ──────────────────────────────────────────
INCLUDE_TASKS = {"FRAMESHIFT_GOF", "STOP_GAIN_GOF", "STOP_LOSS_GOF"}


def _include_in_core(row):
    task    = row["primary_gof_task"]
    r1_incl = row.get("include_in_final_dataset_round1", "")

    # Out-of-scope or excluded → NO
    if task in ("NON_GOF_EXCLUDE", "OUT_OF_SCOPE_BROAD_GOF"):
        return "NO"

    # Round-1 QC said NO → preserve that decision
    if r1_incl == "NO":
        return "NO"

    # Any remaining REVIEW from round-1 QC or task assignment
    if r1_incl == "REVIEW" or task == "REVIEW":
        return "REVIEW"

    # Core task type + round-1 included → YES
    if task in INCLUDE_TASKS and r1_incl == "YES":
        return "YES"

    return "REVIEW"


df["include_in_core_task_dataset"] = df.apply(_include_in_core, axis=1)
print("include_in_core_task_dataset:", df["include_in_core_task_dataset"].value_counts().to_dict())

# ── 9. task_reclassification_reason ──────────────────────────────────────────
def _task_reason(row):
    pvc     = row["primary_variant_class"]
    task    = row["primary_gof_task"]
    mech    = row.get("curated_mechanism_class", "")
    qc      = row.get("qc_round1_action", "")
    fvt     = row.get("final_variant_type", "")
    gof_due = row["gof_due_to_variant"]
    nmd     = row["nmd_escape_status"]
    hn      = row["hn_call"]
    sl_sub  = row["stop_loss_subtype"]

    parts = []

    # NMD_escape remapping note
    if fvt == "NMD_escape":
        parts.append(f"NMD_escape hgvs remapped→{pvc}")

    # Task assignment reason
    if task == "NON_GOF_EXCLUDE":
        if qc.startswith("EXCLUDE"):
            parts.append(f"QC_excluded:{qc}")
        elif mech in NON_GOF_MECHS:
            parts.append("mech=LOF_non_GOF")
        else:
            parts.append("excluded_non_GOF")
    elif task == "OUT_OF_SCOPE_BROAD_GOF":
        if mech in OOS_MECHS:
            parts.append("LOF_mediated_pathway:indirect_GOF_out_of_scope")
        else:
            parts.append(f"variant_type_{pvc}_not_primary_task")

    # Mechanism
    parts.append(f"mech={mech}")

    # hn_call override note
    if hn in ("DOMINANT_NEGATIVE_SEPARATE", "TOXIC_GOF_SEPARATE",
              "NMD_ESCAPE_GOF_REVIEW", "NEOMORPH"):
        parts.append(f"hn_call={hn}")

    # Stop-loss subtype
    if sl_sub != "NOT_APPLICABLE":
        parts.append(f"stop_loss_subtype={sl_sub}")

    # NMD note when relevant
    if nmd in ("NMD_ESCAPE_REPORTED", "NMD_ESCAPE_INFERRED", "NMD_SENSITIVE_OR_DEGRADED"):
        parts.append(f"nmd={nmd}")

    # GOF determination caveat
    if gof_due == "NO":
        parts.append("no_GOF_evidence")
    elif gof_due == "CONFLICTING":
        parts.append("conflicting_GOF_signals")

    # QC review flag
    if qc.startswith("REVIEW"):
        parts.append(f"qc_flag:{qc}")

    return "; ".join(parts) if parts else "standard_assignment"


df["task_reclassification_reason"] = df.apply(_task_reason, axis=1)

# ── 10. Column ordering ───────────────────────────────────────────────────────
NEW_COLS = [
    "primary_variant_class",
    "stop_loss_subtype",
    "primary_gof_task",
    "gof_due_to_variant",
    "gof_evidence_strength",
    "nmd_escape_status",
    "hn_call",
    "include_in_core_task_dataset",
    "task_reclassification_reason",
]
QC_FRONT = [c for c in df.columns
            if c not in NEW_COLS and
            any(c.startswith(p) for p in
                ("qc_round1", "include_in_final", "include_in_strict", "hn_view",
                 "curated_mechanism_class_round1", "strict_gof_status_round1",
                 "grobid_context"))]
REST = [c for c in df.columns if c not in NEW_COLS and c not in QC_FRONT]
df = df[NEW_COLS + QC_FRONT + REST]

# ── 11. Subsets ───────────────────────────────────────────────────────────────
df_yes    = df[df["include_in_core_task_dataset"] == "YES"].copy()
df_review = df[df["include_in_core_task_dataset"] == "REVIEW"].copy()
df_no     = df[df["include_in_core_task_dataset"] == "NO"].copy()

df_fs = df_yes[df_yes["primary_gof_task"] == "FRAMESHIFT_GOF"].copy()
df_sg = df_yes[df_yes["primary_gof_task"] == "STOP_GAIN_GOF"].copy()
df_sl = df_yes[df_yes["primary_gof_task"] == "STOP_LOSS_GOF"].copy()

# ── 12. Print summary ─────────────────────────────────────────────────────────
print(f"\n{'='*55}")
print(f"Task Reclassification Summary  (total={len(df)})")
print(f"{'='*55}")
print(f"\ninclude_in_core_task_dataset:")
print(f"  YES={len(df_yes)}  REVIEW={len(df_review)}  NO={len(df_no)}")
print(f"\nCore task YES breakdown:")
print(f"  FRAMESHIFT_GOF : {len(df_fs)}")
print(f"  STOP_GAIN_GOF  : {len(df_sg)}")
print(f"  STOP_LOSS_GOF  : {len(df_sl)}")
print(f"\nhn_call (YES only):")
for k, v in df_yes["hn_call"].value_counts().items():
    print(f"  {k}: {v}")
print(f"\ngof_evidence_strength (YES only):")
for k, v in df_yes["gof_evidence_strength"].value_counts().items():
    print(f"  {k}: {v}")
print(f"\nnmd_escape_status (YES only):")
for k, v in df_yes["nmd_escape_status"].value_counts().items():
    print(f"  {k}: {v}")
print(f"\ngof_due_to_variant (YES only):")
for k, v in df_yes["gof_due_to_variant"].value_counts().items():
    print(f"  {k}: {v}")

# ── 13. Excel export ──────────────────────────────────────────────────────────
summary_rows = [
    ["=== GOF Pipeline — Task Reclassification v1 ===", ""],
    ["", ""],
    ["[ include_in_core_task_dataset ]", ""],
    ["YES",    len(df_yes)],
    ["REVIEW", len(df_review)],
    ["NO",     len(df_no)],
    ["", ""],
    ["[ Core task split (YES only) ]", ""],
    ["FRAMESHIFT_GOF", len(df_fs)],
    ["STOP_GAIN_GOF",  len(df_sg)],
    ["STOP_LOSS_GOF",  len(df_sl)],
    ["", ""],
    ["[ hn_call (YES only) ]", ""],
]
for k, v in df_yes["hn_call"].value_counts().items():
    summary_rows.append([f"  {k}", int(v)])
summary_rows += [["", ""], ["[ gof_evidence_strength (YES only) ]", ""]]
for k, v in df_yes["gof_evidence_strength"].value_counts().items():
    summary_rows.append([f"  {k}", int(v)])
summary_rows += [["", ""], ["[ nmd_escape_status (YES only) ]", ""]]
for k, v in df_yes["nmd_escape_status"].value_counts().items():
    summary_rows.append([f"  {k}", int(v)])
summary_rows += [["", ""], ["[ gof_due_to_variant (YES only) ]", ""]]
for k, v in df_yes["gof_due_to_variant"].value_counts().items():
    summary_rows.append([f"  {k}", int(v)])
summary_rows += [
    ["", ""],
    ["[ primary_gof_task distribution (all 2486) ]", ""],
]
for k, v in df["primary_gof_task"].value_counts().items():
    summary_rows.append([f"  {k}", int(v)])

summary_df = pd.DataFrame(summary_rows, columns=["항목", "값"])

with pd.ExcelWriter(OUT_XL, engine="openpyxl") as writer:
    _clean_df(summary_df).to_excel( writer, sheet_name="Summary",       index=False)
    _clean_df(df).to_excel(         writer, sheet_name="All_2486",      index=False)
    _clean_df(df_yes).to_excel(     writer, sheet_name="Core_YES",      index=False)
    _clean_df(df_review).to_excel(  writer, sheet_name="Core_REVIEW",   index=False)
    _clean_df(df_no).to_excel(      writer, sheet_name="Core_NO",       index=False)
    _clean_df(df_fs).to_excel(      writer, sheet_name="FRAMESHIFT_GOF", index=False)
    _clean_df(df_sg).to_excel(      writer, sheet_name="STOP_GAIN_GOF",  index=False)
    _clean_df(df_sl).to_excel(      writer, sheet_name="STOP_LOSS_GOF",  index=False)

    for sn in writer.sheets:
        ws = writer.sheets[sn]
        for col_cells in ws.columns:
            ml = max((len(str(c.value)) if c.value else 0 for c in col_cells), default=0)
            ws.column_dimensions[col_cells[0].column_letter].width = min(ml + 2, 60)

print(f"\nSaved: {OUT_XL}")
