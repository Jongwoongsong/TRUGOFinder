#!/bin/bash
#SBATCH --job-name=sl_extract
#SBATCH --partition=debug
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/home/jihoon_yoon/jwsong/logs/slurm_sl_extract_%j.log
#SBATCH --error=/home/jihoon_yoon/jwsong/logs/slurm_sl_extract_%j.log

# ─────────────────────────────────────────────────────────────────────────────
# 48b_stop_loss_extract.sbatch
# Stop-loss 증강 파이프라인 — Phase B: fulltext 취득 + LLM 추출
#
# Input:  data/audit_stop_loss_precise_pmids.txt  (230 PMIDs)
# Step 1: 기존 corpus에 있는 23개 → 직접 추출
# Step 2: 새 207개 PMID → PMC EFetch로 fulltext 취득 (없으면 abstract)
#         → data/fulltext_stoploss/ 에 저장
# Step 3: 전체 대상 대해 LLM 추출 (Qwen3-32B-AWQ, port 8004)
# Output: data/extracted_variants_stoploss_v1.jsonl
# ─────────────────────────────────────────────────────────────────────────────

set -e
cd /home/jihoon_yoon/jwsong/projects/gof_pipeline
source /home/jihoon_yoon/jwsong/miniconda3/etc/profile.d/conda.sh
conda activate gof_pipeline
PYTHON=python3

export VLLM_BASE=http://localhost:8004/v1
PMID_LIST=data/audit_stop_loss_precise_pmids.txt
OUT_JSONL=data/extracted_variants_stoploss_v1.jsonl
FT_SL_DIR=data/fulltext_stoploss
CORPUS_DIR=data/fulltext

echo "=== Stop-Loss Augmentation: Phase B ==="
echo "Input PMID list: $(wc -l < $PMID_LIST) PMIDs"
date
mkdir -p $FT_SL_DIR
mkdir -p logs

# ── Step 1: vLLM 시작 ────────────────────────────────────────────────────────
echo ""
echo "[vLLM] 기존 port 8004 정리 후 시작..."
fuser -k 8004/tcp 2>/dev/null || true
sleep 3

VLLM_USE_V1=0 /home/jihoon_yoon/jwsong/miniconda3/envs/gof_pipeline/bin/vllm serve \
    Qwen/Qwen3-32B-AWQ \
    --quantization awq \
    --max-model-len 12288 \
    --gpu-memory-utilization 0.92 \
    --port 8004 \
    --dtype half \
    --enforce-eager \
    > /home/jihoon_yoon/jwsong/logs/vllm_sl_${SLURM_JOB_ID}.log 2>&1 &

VLLM_PID=$!
echo "  vLLM PID: $VLLM_PID"

echo -n "  vLLM 대기"
for i in $(seq 1 120); do
    sleep 5
    MODEL_CHECK=$(curl -s http://localhost:8004/v1/models 2>/dev/null || echo "")
    if echo "$MODEL_CHECK" | python3 -c "import sys,json; d=json.load(sys.stdin); exit(0 if d.get('data') else 1)" 2>/dev/null; then
        echo " OK (${i}×5s)"
        break
    fi
    echo -n "."
done

# ── Step 2: Fulltext 취득 (PMC EFetch + fallback to abstract) ─────────────
echo ""
echo "[Step 2] Fulltext 취득..."
$PYTHON - <<'PYEOF'
import json
import re
import time
import sys
from pathlib import Path
import requests

HERE     = Path(".")
CORPUS   = HERE / "data/fulltext"
FT_SL    = HERE / "data/fulltext_stoploss"
PMID_F   = HERE / "data/audit_stop_loss_precise_pmids.txt"

EMAIL    = "jongwoongsong@gmail.com"
EBASE    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Load target PMIDs
pmids = [l.strip() for l in open(PMID_F) if l.strip()]
print(f"Target PMIDs: {len(pmids)}")

# Build corpus PMID→file map
corpus_map = {}
for f in CORPUS.glob("*.json"):
    with open(f) as fp:
        d = json.load(fp)
    if d.get("pmid"):
        corpus_map[str(d["pmid"])] = f

already_in_corpus = sum(1 for p in pmids if p in corpus_map)
need_fetch = [p for p in pmids if p not in corpus_map]
print(f"Already in corpus: {already_in_corpus}")
print(f"Need to fetch:     {len(need_fetch)}")


def _clean(text):
    """Strip XML tags and normalize whitespace."""
    text = re.sub(r"<[^>]+>", " ", text)
    text = re.sub(r"&[a-z]+;", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def fetch_pmc_fulltext(pmid: str) -> dict | None:
    """Try PMC EFetch for full-text XML; parse into our JSON schema."""
    try:
        r = requests.get(f"{EBASE}/elink.fcgi",
            params={"dbfrom":"pubmed","db":"pmc","id":pmid,
                    "retmode":"json","email":EMAIL},
            timeout=30)
        if r.status_code != 200:
            return None
        try:
            linkdata = json.loads(r.text, strict=False)
        except Exception:
            return None
        try:
            pmcids = linkdata["linksets"][0]["linksetdbs"][0]["links"]
        except (KeyError, IndexError):
            return None
        if not pmcids:
            return None
        pmcid = str(pmcids[0])

        r2 = requests.get(f"{EBASE}/efetch.fcgi",
            params={"db":"pmc","id":pmcid,"rettype":"full","retmode":"xml","email":EMAIL},
            timeout=45)
        if r2.status_code != 200:
            return None

        xml = r2.text
        title_m = re.search(r"<article-title[^>]*>(.*?)</article-title>", xml, re.S)
        title = _clean(title_m.group(1)) if title_m else ""

        sections = {}
        abs_m = re.search(r"<abstract[^>]*>(.*?)</abstract>", xml, re.S)
        if abs_m:
            sections["abstract"] = _clean(abs_m.group(1))
        for sec_m in re.finditer(
            r"<sec[^>]*>.*?<title[^>]*>(.*?)</title>(.*?)</sec>", xml, re.S
        ):
            sec_title = _clean(sec_m.group(1)).lower()[:30]
            sec_text  = _clean(sec_m.group(2))
            if len(sec_text) > 50:
                sections[sec_title] = sec_text

        if not sections:
            return None

        chars = sum(len(v) for v in sections.values())
        return {
            "pmid":          pmid,
            "pmcid":         f"PMC{pmcid}",
            "doi":           "",
            "score_overall": 0.7,
            "title":         title,
            "source":        "pmc_efetch_stoploss",
            "sections":      sections,
            "n_sections":    len(sections),
            "chars_total":   chars,
        }
    except Exception:
        return None


def fetch_abstract(pmid: str) -> dict | None:
    """Fallback: fetch abstract only via PubMed EFetch."""
    try:
        r = requests.get(f"{EBASE}/efetch.fcgi",
            params={"db":"pubmed","id":pmid,"rettype":"abstract","retmode":"xml","email":EMAIL},
            timeout=30)
        if r.status_code != 200:
            return None

        xml = r.text
        title_m = re.search(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>", xml, re.S)
        title   = _clean(title_m.group(1)) if title_m else ""
        abs_m   = re.search(r"<AbstractText[^>]*>(.*?)</AbstractText>", xml, re.S)
        if not abs_m:
            return None
        abstract = _clean(abs_m.group(1))
        if len(abstract) < 50:
            return None

        return {
            "pmid":          pmid,
            "pmcid":         "",
            "doi":           "",
            "score_overall": 0.4,
            "title":         title,
            "source":        "pubmed_abstract_stoploss",
            "sections":      {"abstract": abstract},
            "n_sections":    1,
            "chars_total":   len(abstract),
        }
    except Exception:
        return None


fetched_ok    = 0
fetched_abs   = 0
fetch_failed  = 0

for i, pmid in enumerate(need_fetch):
    out_path = FT_SL / f"PMID{pmid}.json"
    if out_path.exists():
        fetched_ok += 1
        continue

    if i % 20 == 0:
        print(f"  [{i}/{len(need_fetch)}] fetching...", flush=True)

    doc = fetch_pmc_fulltext(pmid)
    time.sleep(0.4)

    if doc is None:
        doc = fetch_abstract(pmid)
        time.sleep(0.4)
        if doc:
            fetched_abs += 1
        else:
            fetch_failed += 1
            continue
    else:
        fetched_ok += 1

    with open(out_path, "w") as f:
        json.dump(doc, f, ensure_ascii=False)

print(f"\nFetch results:")
print(f"  PMC fulltext:  {fetched_ok}")
print(f"  Abstract only: {fetched_abs}")
print(f"  Failed:        {fetch_failed}")
PYEOF

# ── Step 3: LLM 추출 ────────────────────────────────────────────────────────
echo ""
echo "[Step 3] LLM 추출 (Qwen3-32B-AWQ)..."
$PYTHON - <<'PYEOF'
import asyncio
import json
import re
import sys
import time
from pathlib import Path

import aiohttp

HERE     = Path(".")
CORPUS   = HERE / "data/fulltext"
FT_SL    = HERE / "data/fulltext_stoploss"
PMID_F   = HERE / "data/audit_stop_loss_precise_pmids.txt"
OUT_F    = HERE / "data/extracted_variants_stoploss_v1.jsonl"

import os
VLLM_BASE = os.environ.get("VLLM_BASE", "http://localhost:8004/v1")
MODEL     = "Qwen/Qwen3-32B-AWQ"
CONCURRENCY = 4

# Load target PMIDs
target_pmids = set(l.strip() for l in open(PMID_F) if l.strip())

# Already extracted?
done_pmids = set()
if OUT_F.exists():
    with open(OUT_F) as f:
        for line in f:
            d = json.loads(line)
            if d.get("pmid"):
                done_pmids.add(str(d["pmid"]))
print(f"Already extracted: {len(done_pmids)}")

# Build file list: prefer fulltext_stoploss, fall back to corpus
corpus_map = {}
for fdir in [CORPUS, FT_SL]:
    for fp in fdir.glob("*.json"):
        with open(fp) as f:
            d = json.load(f)
        pmid = str(d.get("pmid",""))
        if pmid and pmid not in corpus_map:
            corpus_map[pmid] = fp

todo = [(pmid, corpus_map[pmid]) for pmid in target_pmids
        if pmid in corpus_map and pmid not in done_pmids]
print(f"To extract: {len(todo)}")

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
- is_human: true ONLY if (1) the gene is an HGNC-approved human gene symbol AND
  (2) evidence involves human patients/cells/tissue.

- class "H": DIRECT gain-of-function evidence. Assign H if ANY of the following:
  * Constitutive activation, enhanced/hyperactive signaling
  * Novel function from C-terminal extension (stop codon loss/readthrough creating new domain)
  * Protein extension beyond natural stop codon causing dominant disease
  * Increased protein stability, resistance to degradation
  * NMD escape producing stable truncated/extended protein with demonstrated GOF activity

- class "N": assign N if:
  * Toxic protein aggregation or cytoplasmic mislocalization
  * Dominant negative — mutant inhibits wild-type function
  * Uncertain / ambiguous mechanism
  * Frameshift/truncation where GOF is inferred but not experimentally confirmed

IMPORTANT: This batch focuses on stop-loss (nonstop/extension) variants.
Prioritize finding:
  - Stop codon mutations creating protein C-terminal extensions (p.*NNNextXX notation)
  - Nonstop extension mutations causing dominant disease
  - Stop codon loss leading to gain-of-function (novel domain, ligand-independent activity, etc.)
  Return findings for ALL frameshift/stop-gain/stop-loss GOF variants — do not limit to stop-loss only.
- evidence_text: direct verbatim quote, max 200 chars
"""

def build_user_msg(doc: dict) -> str:
    sections = doc.get("sections", {})
    title    = doc.get("title", "")
    parts    = [f"Title: {title}"] if title else []
    for sec_name, text in sections.items():
        if text and len(str(text)) > 20:
            parts.append(f"[{sec_name.upper()}]\n{str(text)[:4000]}")
    return "\n\n".join(parts)[:11000]


async def extract_one(session, pmid: str, fpath: Path, sem: asyncio.Semaphore) -> dict:
    async with sem:
        with open(fpath) as fp:
            doc = json.load(fp)

        user_msg = build_user_msg(doc)
        if len(user_msg) < 100:
            return {"pmid": pmid, "pmcid": doc.get("pmcid",""), "findings": [],
                    "error": "text_too_short"}

        payload = {
            "model": MODEL,
            "messages": [
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user",   "content": user_msg},
            ],
            "temperature":     0.0,
            "max_tokens":      2048,
            "response_format": {"type": "json_object"},
        }

        for attempt in range(3):
            try:
                async with session.post(
                    f"{VLLM_BASE}/chat/completions",
                    json=payload, timeout=aiohttp.ClientTimeout(total=120)
                ) as resp:
                    resp.raise_for_status()
                    data    = await resp.json()
                    content = data["choices"][0]["message"]["content"]
                    parsed  = json.loads(content)
                    findings = parsed.get("findings", [])
                    return {
                        "pmid":    pmid,
                        "pmcid":   doc.get("pmcid",""),
                        "title":   doc.get("title",""),
                        "source":  doc.get("source",""),
                        "findings": findings,
                    }
            except Exception as e:
                if attempt == 2:
                    return {"pmid": pmid, "pmcid": doc.get("pmcid",""),
                            "findings": [], "error": str(e)}
                await asyncio.sleep(2 ** attempt)


async def main():
    sem = asyncio.Semaphore(CONCURRENCY)
    connector = aiohttp.TCPConnector(limit=CONCURRENCY + 2)
    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = [extract_one(session, pmid, fpath, sem) for pmid, fpath in todo]
        results = []
        for i, coro in enumerate(asyncio.as_completed(tasks)):
            result = await coro
            results.append(result)
            if (i + 1) % 20 == 0:
                n_findings = sum(len(r.get("findings",[])) for r in results)
                print(f"  [{i+1}/{len(todo)}] findings so far: {n_findings}", flush=True)

    with open(OUT_F, "a") as f:
        for r in results:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    total_findings = sum(len(r.get("findings",[])) for r in results)
    errors         = sum(1 for r in results if r.get("error"))
    print(f"\nExtraction done:")
    print(f"  Papers processed: {len(results)}")
    print(f"  Findings extracted: {total_findings}")
    print(f"  Errors: {errors}")
    print(f"  Output: {OUT_F}")


asyncio.run(main())
PYEOF

# ── Step 4: 결과 요약 ─────────────────────────────────────────────────────────
echo ""
echo "[Step 4] 추출 결과 요약..."
$PYTHON - <<'PYEOF'
import json
from pathlib import Path
from collections import Counter

OUT_F = Path("data/extracted_variants_stoploss_v1.jsonl")
if not OUT_F.exists():
    print("Output file not found!")
    exit(0)

total_papers   = 0
total_findings = 0
class_counter  = Counter()
is_human_counter = Counter()
stop_loss_hits = 0

with open(OUT_F) as f:
    for line in f:
        d = json.loads(line)
        total_papers += 1
        findings = d.get("findings", [])
        total_findings += len(findings)
        for fnd in findings:
            class_counter[fnd.get("class","?")] += 1
            is_human_counter[str(fnd.get("is_human","?"))] += 1
            # Stop-loss specific check
            prot = str(fnd.get("hgvs_protein",""))
            if "ext" in prot.lower() and ("Ter" in prot or "*" in prot):
                stop_loss_hits += 1

print(f"Papers processed:   {total_papers}")
print(f"Total findings:     {total_findings}")
print(f"Stop-loss variants: {stop_loss_hits}")
print(f"Class H:            {class_counter.get('H',0)}")
print(f"Class N:            {class_counter.get('N',0)}")
print(f"is_human=True:      {is_human_counter.get('True',0)}")
PYEOF

# ── 종료 ─────────────────────────────────────────────────────────────────────
echo ""
echo "[vLLM] 종료..."
kill $VLLM_PID 2>/dev/null || true

echo ""
echo "=== 완료 ==="
echo "다음 단계: python pipeline/48c_stop_loss_merge.py"
date
