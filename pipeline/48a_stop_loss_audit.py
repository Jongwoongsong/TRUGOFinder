"""
48a_stop_loss_audit.py
Stop-loss 증강 파이프라인 — Phase A: 기존 corpus 감사 + PubMed 신규 PMID 발굴

Phase 1  21,042 개 fulltext JSON 스캔 → stop-loss 키워드 포함 문서 목록
Phase 2  PubMed E-utilities 검색 → stop-loss GOF 관련 신규 PMID
Phase 3  Gap 분석 → 기존 corpus에 없는 신규 PMID 목록 (48b SLURM 입력)

Output files:
  data/audit_stop_loss_in_corpus.tsv   - 기존 corpus hit 목록 (pmid, 매칭 섹션, 매칭 패턴)
  data/audit_stop_loss_pubmed_raw.tsv  - PubMed 검색 전체 결과 (pmid, query, title)
  data/audit_stop_loss_new_pmids.txt   - 신규 PMID (48b LLM 추출 입력)
  data/audit_stop_loss_summary.txt     - 전체 요약
"""
import json
import re
import time
import sys
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent.parent
FT_DIR = HERE / "data/fulltext"

OUT_CORPUS   = HERE / "data/audit_stop_loss_in_corpus.tsv"
OUT_PUBMED   = HERE / "data/audit_stop_loss_pubmed_raw.tsv"
OUT_NEW      = HERE / "data/audit_stop_loss_new_pmids.txt"
OUT_PRECISE  = HERE / "data/audit_stop_loss_precise_pmids.txt"  # filtered — 48b input
OUT_SUMMARY  = HERE / "data/audit_stop_loss_summary.txt"

# HGVS_ext query is too broad (2000 hits, 0 overlap with other queries).
# Exclude PMIDs that appear ONLY from this query.
BROAD_LABELS = {"HGVS_ext"}

NCBI_BASE  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = "jongwoongsong@gmail.com"   # required by NCBI policy

# ── Stop-loss keyword patterns ────────────────────────────────────────────────
# Tier-1: specific stop-loss terms (high precision)
TIER1_RE = re.compile(
    r"\bstop[- ]loss\b|"
    r"\bnonstop[- ](?:extension|mutation|variant|change)\b|"
    r"\bnon-?stop\s+(?:codon\s+)?(?:mutation|variant|extension)\b|"
    r"\bloss\s+of\s+(?:the\s+)?stop\s+codon\b|"
    r"\bstop\s+codon\s+(?:loss|elimination)\b",
    re.I,
)

# Tier-2: extension proteins from stop codon mutations (HGVS-based)
TIER2_RE = re.compile(
    r"\bext\*\d+\b|"               # ext*38, ext*19
    r"\bextTer\d+\b|"              # extTer93
    r"\bext[A-Z][a-z]{2}\*\d+\b|" # extLeu*42
    r"p\.\*\d+[A-Za-z]+ext\b|"    # p.*278Glyext
    r"p\.Ter\d+[A-Za-z]+ext\b|"   # p.Ter1153Lysext
    r"\bp\.\([A-Z][a-z]{2}\d+[A-Z][a-z]{2}ext",  # p.(Ter402Leuext
    re.I,
)

# Tier-3: C-terminal extension near stop codon context
TIER3_RE = re.compile(
    r"C[- ]terminal\s+extension.{0,60}stop\s+codon|"
    r"stop\s+codon.{0,60}C[- ]terminal\s+extension|"
    r"protein\s+extension.{0,40}stop\s+codon|"
    r"stop\s+codon\s+mutation.{0,60}extension|"
    r"translational\s+readthrough\b|"
    r"\breadthrough\s+(?:suppression|mutation|variant)\b",
    re.I,
)


def get_all_text(sections: dict) -> str:
    """Concatenate all section texts into a single string for regex matching."""
    return " ".join(str(v) for v in sections.values() if v)


def find_matches(text: str) -> dict:
    """Return dict of tier → matched substring (first match only, for logging)."""
    hits = {}
    m1 = TIER1_RE.search(text)
    if m1:
        hits["tier1"] = m1.group(0)
    m2 = TIER2_RE.search(text)
    if m2:
        hits["tier2"] = m2.group(0)
    m3 = TIER3_RE.search(text)
    if m3:
        hits["tier3"] = m3.group(0)
    return hits


# ══════════════════════════════════════════════════════════════════════════════
# Phase 1: Scan existing corpus
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 60)
print("Phase 1: Scanning existing corpus for stop-loss terms")
print("=" * 60)

all_json = sorted(FT_DIR.glob("*.json"))
print(f"  Files to scan: {len(all_json):,}")

corpus_hits = []       # rows for OUT_CORPUS
corpus_pmid_set = set()  # all PMIDs in corpus (for Phase 3 gap analysis)

for i, fpath in enumerate(all_json):
    if i % 2000 == 0:
        print(f"  [{i:5d}/{len(all_json)}] scanning...")
        sys.stdout.flush()

    with open(fpath) as fp:
        doc = json.load(fp)

    pmid  = str(doc.get("pmid", "")).strip()
    pmcid = str(doc.get("pmcid", "")).strip()
    title = str(doc.get("title", "")).strip()

    if pmid:
        corpus_pmid_set.add(pmid)

    sections = doc.get("sections", {})
    text     = get_all_text(sections)
    hits     = find_matches(text)

    if hits:
        tier = "tier1" if "tier1" in hits else ("tier2" if "tier2" in hits else "tier3")
        matched_term = " | ".join(f"{k}:{v}" for k, v in hits.items())
        hit_sections = [k for k, v in sections.items()
                        if v and (TIER1_RE.search(str(v)) or
                                  TIER2_RE.search(str(v)) or
                                  TIER3_RE.search(str(v)))]
        corpus_hits.append({
            "pmid":          pmid,
            "pmcid":         pmcid,
            "filename":      fpath.name,
            "title":         title[:120],
            "best_tier":     tier,
            "matched_terms": matched_term[:200],
            "hit_sections":  "|".join(hit_sections),
        })

print(f"\n  Corpus scan done.")
print(f"  Total corpus PMIDs: {len(corpus_pmid_set):,}")
print(f"  Stop-loss hits in corpus: {len(corpus_hits)}")

tier_counts = {}
for h in corpus_hits:
    tier_counts[h["best_tier"]] = tier_counts.get(h["best_tier"], 0) + 1
for t, c in sorted(tier_counts.items()):
    print(f"    {t}: {c}")

corpus_hits_df = pd.DataFrame(corpus_hits)
corpus_hits_df.to_csv(OUT_CORPUS, sep="\t", index=False)
print(f"\n  Saved: {OUT_CORPUS}  ({len(corpus_hits_df)} rows)")

# ══════════════════════════════════════════════════════════════════════════════
# Phase 2: PubMed E-utilities search
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("Phase 2: PubMed E-utilities search")
print("=" * 60)

# Queries designed to find stop-loss / nonstop GOF variants in human disease.
# Intentionally separate human-disease and mechanism terms to maximize recall.
PUBMED_QUERIES = [
    # Direct stop-loss terminology
    ('"stop-loss" variant disease',            "SL_direct"),
    ('"stop-loss" mutation pathogenic',        "SL_direct"),
    # Nonstop extension
    ('"nonstop" mutation disease',             "nonstop"),
    ('nonstop extension mutation pathogenic',  "nonstop"),
    # C-terminal extension from stop codon mutation
    ('stop codon loss gain of function human', "SL_GOF"),
    ('stop codon mutation extension pathogenic gain function', "SL_GOF"),
    ('c-terminal extension stop codon dominant pathogenic',   "SL_ext"),
    # Specific to our known genes/diseases (for recall expansion)
    ('protein extension "loss of stop codon" disease',        "SL_ext"),
    # Readthrough in genetic context
    ('readthrough mutation gain function disease',            "readthrough"),
    ('translational readthrough pathogenic variant',          "readthrough"),
    # HGVS extension notation
    ('extTer pathogenic variant disease',                     "HGVS_ext"),
]


def esearch_all(query: str, retmax: int = 500) -> list[str]:
    """Return all PMIDs for a PubMed query (paginated)."""
    params = {
        "db":       "pubmed",
        "term":     query,
        "retmax":   retmax,
        "retmode":  "json",
        "email":    NCBI_EMAIL,
        "tool":     "gof_pipeline_stoploss_audit",
    }
    r = requests.get(f"{NCBI_BASE}/esearch.fcgi", params=params, timeout=20)
    r.raise_for_status()
    result = r.json()["esearchresult"]
    total  = int(result.get("count", 0))
    pmids  = result.get("idlist", [])

    # Paginate if necessary
    retstart = len(pmids)
    while retstart < total and retstart < 2000:  # hard cap 2000 per query
        params2 = {**params, "retstart": retstart}
        r2 = requests.get(f"{NCBI_BASE}/esearch.fcgi", params=params2, timeout=20)
        r2.raise_for_status()
        batch = r2.json()["esearchresult"].get("idlist", [])
        if not batch:
            break
        pmids.extend(batch)
        retstart += len(batch)
        time.sleep(0.4)

    return pmids


def efetch_titles(pmids: list[str]) -> dict[str, str]:
    """Fetch titles for a batch of PMIDs (max 200 at a time)."""
    title_map = {}
    for i in range(0, len(pmids), 200):
        batch = pmids[i:i + 200]
        params = {
            "db":       "pubmed",
            "id":       ",".join(batch),
            "rettype":  "abstract",
            "retmode":  "xml",
            "email":    NCBI_EMAIL,
        }
        r = requests.get(f"{NCBI_BASE}/efetch.fcgi", params=params, timeout=30)
        r.raise_for_status()
        # Simple regex extraction — no lxml dependency required
        for pmid, title in zip(
            re.findall(r"<PMID[^>]*>(\d+)</PMID>", r.text),
            re.findall(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>", r.text, re.S),
        ):
            title_map[pmid] = re.sub(r"<[^>]+>", "", title).strip()
        time.sleep(0.4)
    return title_map


pubmed_rows = []   # {pmid, query_label, query, title}
seen_pmids  = set()

for query_str, label in PUBMED_QUERIES:
    print(f"  [{label}] {query_str}")
    try:
        pmids = esearch_all(query_str)
        print(f"    → {len(pmids)} PMIDs")
        for pid in pmids:
            if pid not in seen_pmids:
                seen_pmids.add(pid)
                pubmed_rows.append({
                    "pmid":        pid,
                    "query_label": label,
                    "query":       query_str,
                    "title":       "",
                })
    except Exception as e:
        print(f"    ERROR: {e}")
    time.sleep(0.5)

print(f"\n  Unique PMIDs from all queries: {len(seen_pmids)}")

# Fetch titles for all unique PMIDs
print("  Fetching titles...")
all_pids = list(seen_pmids)
title_map = efetch_titles(all_pids)
for row in pubmed_rows:
    row["title"] = title_map.get(row["pmid"], "")

pubmed_df = pd.DataFrame(pubmed_rows)
pubmed_df.to_csv(OUT_PUBMED, sep="\t", index=False)
print(f"  Saved: {OUT_PUBMED}  ({len(pubmed_df)} rows)")

# ══════════════════════════════════════════════════════════════════════════════
# Phase 3: Gap analysis — new PMIDs not in corpus
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("Phase 3: Gap analysis")
print("=" * 60)

# Load existing task dataset PMIDs (already processed + classified)
task_df   = pd.read_excel(
    HERE / "data/outputs/gof_task_reclassified_v1.xlsx",
    sheet_name="All_2486",
    dtype=str,
    usecols=["pmid"],
)
task_pmids = set(task_df["pmid"].dropna().str.strip())

# New PMIDs: from PubMed search, not already in fulltext corpus
pubmed_pmids  = seen_pmids
in_corpus     = pubmed_pmids & corpus_pmid_set
not_in_corpus = pubmed_pmids - corpus_pmid_set
in_task       = pubmed_pmids & task_pmids
new_to_fetch  = not_in_corpus - task_pmids  # need fulltext + LLM extraction

print(f"  PubMed search PMIDs:       {len(pubmed_pmids)}")
print(f"  Already in fulltext corpus: {len(in_corpus)}")
print(f"  Already in task dataset:    {len(in_task)}")
print(f"  NEW (need fetch+extract):  {len(new_to_fetch)}")

# Also flag corpus hits that are not yet in the task dataset
corpus_hit_pmids = set(corpus_hits_df["pmid"].dropna().str.strip())
corpus_new_pmids = corpus_hit_pmids - task_pmids  # in corpus, has stop-loss term, not in task
print(f"\n  Corpus stop-loss hits not yet in task dataset: {len(corpus_new_pmids)}")
print(f"  (These already have fulltext — can run LLM extraction directly)")

# Write new PMID list for 48b (all, including broad HGVS_ext)
all_new = new_to_fetch | corpus_new_pmids
with open(OUT_NEW, "w") as f:
    for pid in sorted(all_new):
        f.write(pid + "\n")
print(f"\n  Total PMIDs for LLM extraction (all): {len(all_new)}")
print(f"  Saved: {OUT_NEW}")

# Write PRECISE list: exclude PMIDs that come only from broad HGVS_ext query
pub_df_all      = pubmed_df.copy()
precise_pmids_pub = set(pub_df_all[~pub_df_all["query_label"].isin(BROAD_LABELS)]["pmid"].unique())
precise_new_pub   = precise_pmids_pub - corpus_pmid_set - task_pmids

# Corpus tier1/tier2 hits not in task (tier3 readthrough papers are low-priority)
corpus_precise_new = set(
    corpus_hits_df[
        (corpus_hits_df["best_tier"].isin(["tier1", "tier2"])) &
        (~corpus_hits_df["pmid"].isin(task_pmids))
    ]["pmid"].dropna()
)

precise_all = precise_new_pub | corpus_precise_new
with open(OUT_PRECISE, "w") as f:
    for pid in sorted(precise_all):
        f.write(pid + "\n")
print(f"\n  PRECISE list (48b input):  {len(precise_all)}")
print(f"    Precise PubMed new:       {len(precise_new_pub)}")
print(f"    Corpus tier1/2 not in task:{len(corpus_precise_new)}")
print(f"  Saved: {OUT_PRECISE}")

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════

summary_lines = [
    "=== Stop-Loss Augmentation Audit ===",
    "",
    f"[Corpus scan]",
    f"  Files scanned:          {len(all_json):,}",
    f"  Unique PMIDs in corpus: {len(corpus_pmid_set):,}",
    f"  Stop-loss hits (total): {len(corpus_hits)}",
    f"  tier1 (direct terms):   {tier_counts.get('tier1', 0)}",
    f"  tier2 (HGVS ext):       {tier_counts.get('tier2', 0)}",
    f"  tier3 (C-term ext):     {tier_counts.get('tier3', 0)}",
    "",
    f"[PubMed search]",
    f"  Queries run:            {len(PUBMED_QUERIES)}",
    f"  Unique PMIDs found:     {len(pubmed_pmids)}",
    f"  Already in corpus:      {len(in_corpus)}",
    f"  Already in task dataset:{len(in_task)}",
    f"  NEW (need fetch+extract):{len(new_to_fetch)}",
    "",
    f"[Gap analysis]",
    f"  Corpus hits not in task: {len(corpus_new_pmids)}",
    f"  PubMed new (no fulltext):{len(new_to_fetch)}",
    f"  Total for 48b LLM job:  {len(all_new)}",
    "",
    f"[Gap analysis — precise (48b input)]",
    f"  Precise PubMed new (no fulltext):  {len(precise_new_pub)}",
    f"  Corpus tier1/2 not in task:        {len(corpus_precise_new)}",
    f"  Total precise for 48b:             {len(precise_all)}",
    f"",
    f"[Output files]",
    f"  {OUT_CORPUS}",
    f"  {OUT_PUBMED}",
    f"  {OUT_NEW}",
    f"  {OUT_PRECISE}   ← 48b input",
]

summary_text = "\n".join(summary_lines)
print("\n" + summary_text)
with open(OUT_SUMMARY, "w") as f:
    f.write(summary_text + "\n")
print(f"\nSaved: {OUT_SUMMARY}")
