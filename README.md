# GOF Variant Extraction Pipeline

LLM-based pipeline for extracting **gain-of-function (GOF) variants** (frameshift, stop-gain, stop-loss) from PubMed/PMC full-text literature.

## Overview

This pipeline identifies variants where truncating or extension mutations confer gain-of-function effects (NMD escape, dominant negative, toxic aggregation, constitutive activation), as opposed to the canonical loss-of-function outcome.

**Model**: Qwen3-32B-AWQ via vLLM (self-hosted, SLURM/gsomics server)  
**Final dataset**: 2,566 GOF variants from 1,214 papers (v10)

---

## Pipeline Structure

```
pipeline/
├── 05_pubmed_fetch.py            Step 1: PubMed/PMC search (~111,424 candidates)
├── 04_screen_abstracts_v3.py     Step 2: LLM abstract screening (score_overall)
├── 10c_download_pmc.py           Step 3: PMC full-text download (top 21,042)
├── 11_extract_variants.py        Step 4: LLM variant extraction (core)
├── 46_qc_round1.py               Step 5: QC round 1 review
├── 47_task_reclassification.py   Step 6: Add 9 classification columns
├── 48a_stop_loss_audit.py        Step 7a: Stop-loss augmentation audit
├── 48c_stop_loss_merge.py        Step 7c: Merge stop-loss findings
├── 49_export_task_split.py       Step 8: Final export (v10 task split)
├── 51_expand_download.py         Step 9: Expand to score<0.7 papers (~74k)
├── hgvs_normalize.py             Utility: HGVS normalization
├── omim_local.py                 Utility: local OMIM DB lookup
└── slurm/
    ├── 48b_stop_loss_extract.sbatch   SLURM: stop-loss LLM extraction
    ├── 52_expand_extract.sbatch       SLURM: expanded LLM extraction
    └── start_vllm.sh                  vLLM server startup

analysis/
└── 50_recall_gap_analysis.py     tmVar/PubTator3 vs model recall comparison
```

---

## Key Outputs (`data/outputs/`)

| File | Description |
|------|-------------|
| `gof_pipeline_final_v10_task_split.xlsx` | **Final dataset** — 2,566 variants, task-split sheets |
| `gof_task_reclassified_v1.xlsx` | Base 2,486 variants with classification columns |
| `gof_curated_after_manual_qc_round1.xlsx` | Manual QC source (gold standard) |
| `gof_stop_loss_augmented_candidates.xlsx` | Stop-loss augmentation candidates (80 new) |
| `recall_gap_analysis.xlsx` | tmVar vs LLM model bidirectional comparison |
| `gof_extraction_benchmark_*.xlsx` | Benchmark vs external tools |

---

## Classification Schema

Each variant is annotated with:

| Column | Values | Description |
|--------|--------|-------------|
| `primary_variant_class` | FRAMESHIFT / STOP_GAIN / STOP_LOSS / OTHER | Variant type |
| `primary_gof_task` | FRAMESHIFT_GOF / STOP_GAIN_GOF / STOP_LOSS_GOF | Task label |
| `gof_due_to_variant` | YES / NO / REVIEW | GOF confirmed |
| `hn_call` | HYPERMORPH / HYPOMORPH / DOMINANT_NEGATIVE_SEPARATE / ... | Mechanism |
| `nmd_escape_status` | CONFIRMED / LIKELY / UNLIKELY / UNKNOWN | NMD escape |
| `include_in_core_task_dataset` | YES / NO / REVIEW | Core task inclusion |

---

## v10 Dataset Summary

```
Total variants:       2,566
  FRAMESHIFT_GOF:     1,627  (core YES)
  STOP_GAIN_GOF:        476  (core YES)
  STOP_LOSS_GOF:         28  (core YES, augmented from 9 → 28)
  REVIEW:               435
```

---

## Comparison with tmVar/PubTator3

- tmVar recall of our PTC variants: **94.6%** (791/836)
- Our model-unique variants: **~63%** (not found by tmVar)
- Our model covers **fulltext** (21,042 papers); tmVar is abstract-based

---

## Setup

```bash
conda env create -f environment.yml
conda activate gof_pipeline
```

**Requirements**: SLURM cluster with GPU (A100 or similar), NCBI Entrez account.

---

## Data Availability

Raw data (full-text corpus, screening results) is not included due to size (~1.3 GB).  
Final annotated outputs are in `data/outputs/`.
