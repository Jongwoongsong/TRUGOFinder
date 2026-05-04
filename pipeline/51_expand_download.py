"""
51_expand_download.py
Abstract 스크리닝 필터(score_overall < 0.7)에서 탈락한 논문 fulltext 다운로드

Rationale:
  현재 파이프라인은 abstract LLM 스코어 >= 0.7 논문만 fulltext를 취득했으나,
  abstract에 GOF 정보가 없어도 fulltext에 있을 수 있음 (false negative).
  → score < 0.7인 77,798개 논문도 PMC fulltext 취득 후 추출 대상에 추가.

Input:
  data/screening_results.jsonl  (111,424건, score_overall 포함)
  data/fulltext/                (기존 21,042개 — 중복 다운로드 방지용)

Output:
  data/fulltext/{PMCID}.json   (기존과 동일 포맷, 기존 디렉토리에 추가)
  data/expand_download_progress.json  (체크포인트 — 재실행 시 이어받기)

Usage:
  python3 51_expand_download.py
  python3 51_expand_download.py --dry-run   # 다운로드 없이 대상 수만 출력
  python3 51_expand_download.py --batch-size 500  # ELink 배치 크기 조정
  nohup python3 51_expand_download.py > logs/expand_download.log 2>&1 &
"""
import argparse
import asyncio
import json
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import aiohttp

HERE = Path(__file__).resolve().parent.parent
SCREEN_JSONL = HERE / "data/screening_results.jsonl"
FT_DIR       = HERE / "data/fulltext"
PROGRESS_F   = HERE / "data/expand_download_progress.json"

EUTILS_BASE  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL        = "jongwoongsong@gmail.com"
TOOL         = "gof_pipeline_expand"

# NCBI 속도 제한: API key 없이 3 req/s → 세마포어 3, 최소 간격 0.34s
MAX_CONCURRENT = 3
MIN_INTERVAL   = 0.34   # seconds between requests


# ── XML 파싱 (10c_download_pmc.py와 동일 로직) ─────────────────────────────
def parse_pmc_xml(xml_text: str) -> dict:
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return {}

    def get_text(elem) -> str:
        return " ".join(elem.itertext()).strip()

    sections: dict[str, str] = {}

    for t in root.findall(".//article-title"):
        sections["title"] = get_text(t)
        break

    for ab in root.findall(".//abstract"):
        sections["abstract"] = get_text(ab)
        break

    target = {
        "introduction":            "introduction",
        "intro":                   "introduction",
        "background":              "background",
        "methods":                 "methods",
        "materials and methods":   "methods",
        "methods and materials":   "methods",
        "methodology":             "methods",
        "experimental procedures": "methods",
        "experimental section":    "methods",
        "patients and methods":    "methods",
        "subjects and methods":    "methods",
        "results":                 "results",
        "results and discussion":  "results",
        "findings":                "results",
        "discussion":              "discussion",
        "conclusion":              "conclusion",
        "conclusions":             "conclusion",
        "summary":                 "conclusion",
    }
    for sec in root.findall(".//sec"):
        title_el = sec.find("title")
        if title_el is None:
            continue
        raw = (title_el.text or "").strip().lower()
        key = target.get(raw)
        if not key:
            for pat, mapped in target.items():
                if pat in raw:
                    key = mapped
                    break
        if key and key not in sections:
            parts = [get_text(c) for c in sec if c.tag != "title"]
            text  = " ".join(parts).strip()
            if text:
                sections[key] = text

    return sections


# ── Rate limiter ──────────────────────────────────────────────────────────────
class RateLimiter:
    def __init__(self, min_interval: float):
        self._min_interval = min_interval
        self._last = 0.0
        self._lock = asyncio.Lock()

    async def wait(self):
        async with self._lock:
            now = time.monotonic()
            wait_time = self._min_interval - (now - self._last)
            if wait_time > 0:
                await asyncio.sleep(wait_time)
            self._last = time.monotonic()


# ── ELink: PMIDs → PMCIDs (배치) ──────────────────────────────────────────────
async def batch_elink(
    session: aiohttp.ClientSession,
    pmids: list[str],
    semaphore: asyncio.Semaphore,
    rate: RateLimiter,
) -> dict[str, str]:
    """Returns {pmid: pmcid} for those that have PMC fulltext."""
    result = {}
    async with semaphore:
        await rate.wait()
        try:
            params = {
                "dbfrom": "pubmed", "db": "pmc",
                "id": ",".join(pmids),
                "retmode": "json",
                "tool": TOOL, "email": EMAIL,
            }
            async with session.get(
                f"{EUTILS_BASE}/elink.fcgi",
                params=params,
                timeout=aiohttp.ClientTimeout(total=45),
            ) as r:
                if r.status != 200:
                    return result
                text = await r.text()
                try:
                    data = json.loads(text, strict=False)
                except Exception:
                    return result

            for ls in data.get("linksets", []):
                input_ids = [str(i) for i in ls.get("ids", [])]
                for ldb in ls.get("linksetdbs", []):
                    links = ldb.get("links", [])
                    if links and input_ids:
                        # links[0] is the primary PMC article
                        pmcid = f"PMC{links[0]}"
                        for pid in input_ids:
                            result[pid] = pmcid
        except Exception:
            pass
    return result


# ── EFetch: PMCID → parsed sections ──────────────────────────────────────────
async def fetch_and_save(
    session: aiohttp.ClientSession,
    pmid: str,
    pmcid: str,
    score: float,
    doi: str,
    semaphore: asyncio.Semaphore,
    rate: RateLimiter,
) -> str:
    """Fetch XML, parse, save. Returns 'saved', 'skip' (already exists), or 'fail'."""
    out_path = FT_DIR / f"{pmcid}.json"
    if out_path.exists():
        return "skip"

    async with semaphore:
        await rate.wait()
        try:
            params = {
                "db": "pmc",
                "id": pmcid.replace("PMC", ""),
                "rettype": "full", "retmode": "xml",
                "tool": TOOL, "email": EMAIL,
            }
            async with session.get(
                f"{EUTILS_BASE}/efetch.fcgi",
                params=params,
                timeout=aiohttp.ClientTimeout(total=60),
            ) as r:
                if r.status != 200:
                    return "fail"
                xml_text = await r.text()
        except Exception:
            return "fail"

    sections = parse_pmc_xml(xml_text)
    if not sections:
        return "fail"

    title = sections.pop("title", "")
    record = {
        "pmid":        pmid,
        "pmcid":       pmcid,
        "doi":         doi,
        "score_overall": score,
        "title":       title,
        "sections":    sections,
        "n_sections":  len(sections),
        "chars_total": sum(len(v) for v in sections.values()),
    }
    out_path.write_text(json.dumps(record, ensure_ascii=False), encoding="utf-8")
    return "saved"


# ── 메인 ──────────────────────────────────────────────────────────────────────
async def main(batch_size: int, dry_run: bool):
    FT_DIR.mkdir(exist_ok=True)

    # ── 1. 기존 fulltext PMID 목록 구축 ─────────────────────────────────────
    print("Building existing fulltext PMID index...", flush=True)
    existing_pmids: set[str] = set()
    existing_pmcids: set[str] = set()
    for f in FT_DIR.glob("*.json"):
        existing_pmcids.add(f.stem)
        try:
            d = json.loads(f.read_text())
            if d.get("pmid"):
                existing_pmids.add(str(d["pmid"]))
        except Exception:
            pass
    print(f"  Existing fulltext PMIDs: {len(existing_pmids):,}", flush=True)

    # ── 2. 체크포인트 로드 ────────────────────────────────────────────────────
    progress = {"completed_pmids": [], "no_pmc_pmids": [], "failed_pmids": []}
    if PROGRESS_F.exists():
        try:
            progress = json.loads(PROGRESS_F.read_text())
            print(f"  Checkpoint: {len(progress['completed_pmids']):,} done, "
                  f"{len(progress['no_pmc_pmids']):,} no-PMC, "
                  f"{len(progress['failed_pmids']):,} failed", flush=True)
        except Exception:
            pass
    done_pmids = set(progress["completed_pmids"]) | set(progress["no_pmc_pmids"])

    # ── 3. 대상 PMID 추출 (score < 0.7, 미다운로드) ──────────────────────────
    print("Loading screening results...", flush=True)
    targets: list[dict] = []   # [{pmid, doi, score}]
    with open(SCREEN_JSONL) as f:
        for line in f:
            d = json.loads(line)
            pmid  = str(d.get("pmid", "")).strip()
            score = float(d.get("score_overall", 0))
            doi   = str(d.get("doi", "")).strip()
            if not pmid:
                continue
            if score >= 0.7:
                continue
            if pmid in existing_pmids:
                continue
            if pmid in done_pmids:
                continue
            targets.append({"pmid": pmid, "doi": doi, "score": score})

    print(f"  Target PMIDs (score < 0.7, not yet downloaded): {len(targets):,}", flush=True)

    if dry_run:
        print("\n[dry-run] 종료 — 실제 다운로드 없음")
        return

    if not targets:
        print("모두 완료. 다운로드할 논문 없음.")
        return

    # ── 4. Batch ELink: PMID → PMCID 매핑 ───────────────────────────────────
    print(f"\nPhase 1: Batch ELink ({len(targets):,} PMIDs, batch={batch_size})...", flush=True)
    t0 = time.time()

    pmid_to_meta = {t["pmid"]: t for t in targets}
    all_pmids    = [t["pmid"] for t in targets]
    batches      = [all_pmids[i:i+batch_size] for i in range(0, len(all_pmids), batch_size)]

    semaphore  = asyncio.Semaphore(MAX_CONCURRENT)
    rate       = RateLimiter(MIN_INTERVAL)
    pmid_pmcid: dict[str, str] = {}

    connector = aiohttp.TCPConnector(limit=MAX_CONCURRENT + 5)
    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = [batch_elink(session, b, semaphore, rate) for b in batches]
        done = 0
        for coro in asyncio.as_completed(tasks):
            mapping = await coro
            pmid_pmcid.update(mapping)
            done += 1
            if done % 50 == 0 or done == len(batches):
                pct = done / len(batches) * 100
                found = len(pmid_pmcid)
                print(f"  ELink: {done}/{len(batches)} batches ({pct:.0f}%) — "
                      f"{found:,} PMC IDs found", flush=True)

    # PMC 없는 PMID 기록
    no_pmc = [p for p in all_pmids if p not in pmid_pmcid]
    progress["no_pmc_pmids"] = list(set(progress["no_pmc_pmids"]) | set(no_pmc))
    PROGRESS_F.write_text(json.dumps(progress))

    elapsed = time.time() - t0
    print(f"  ELink done: {len(pmid_pmcid):,} have PMC ({len(pmid_pmcid)/len(all_pmids)*100:.1f}%), "
          f"{len(no_pmc):,} no PMC — {elapsed:.0f}s", flush=True)

    # ── 5. EFetch + 저장 ─────────────────────────────────────────────────────
    fetch_targets = [
        (pmid, pmcid)
        for pmid, pmcid in pmid_pmcid.items()
        if pmcid not in existing_pmcids
    ]
    already_exists = len(pmid_pmcid) - len(fetch_targets)
    print(f"\nPhase 2: EFetch {len(fetch_targets):,} papers "
          f"({already_exists} already in fulltext/)...", flush=True)

    saved = failed = skipped = 0
    t0 = time.time()

    save_lock = asyncio.Lock()

    async def process_one(pmid: str, pmcid: str):
        nonlocal saved, failed, skipped
        meta = pmid_to_meta.get(pmid, {})
        status = await fetch_and_save(
            session, pmid, pmcid,
            meta.get("score", 0.0),
            meta.get("doi", ""),
            semaphore, rate,
        )
        async with save_lock:
            if status == "saved":
                saved += 1
                progress["completed_pmids"].append(pmid)
                if saved % 500 == 0:
                    PROGRESS_F.write_text(json.dumps(progress))
                    elapsed = time.time() - t0
                    rate_pps = saved / elapsed if elapsed > 0 else 0
                    remaining = len(fetch_targets) - saved - failed
                    eta_min = remaining / rate_pps / 60 if rate_pps > 0 else 0
                    print(f"  Saved {saved:,} | failed {failed:,} | "
                          f"{rate_pps:.1f} papers/s | ETA {eta_min:.0f}m", flush=True)
            elif status == "fail":
                failed += 1
                progress["failed_pmids"].append(pmid)
            else:
                skipped += 1

    async with aiohttp.ClientSession(
        connector=aiohttp.TCPConnector(limit=MAX_CONCURRENT + 5)
    ) as session:
        tasks = [process_one(pmid, pmcid) for pmid, pmcid in fetch_targets]
        await asyncio.gather(*tasks)

    PROGRESS_F.write_text(json.dumps(progress))

    elapsed = time.time() - t0
    print(f"\n=== 완료 ===")
    print(f"  Saved:   {saved:,}  (data/fulltext/ 에 추가)")
    print(f"  Skipped: {skipped:,} (이미 존재)")
    print(f"  Failed:  {failed:,}")
    print(f"  No PMC:  {len(no_pmc):,}")
    print(f"  Time:    {elapsed/60:.1f}분")
    print(f"\n다음 단계: sbatch 52_expand_extract.sbatch")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true",
                        help="다운로드 없이 대상 수만 확인")
    parser.add_argument("--batch-size", type=int, default=200,
                        help="ELink 배치 크기 (기본 200)")
    args = parser.parse_args()

    asyncio.run(main(batch_size=args.batch_size, dry_run=args.dry_run))
