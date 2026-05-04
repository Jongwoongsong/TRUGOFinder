"""
10c_download_pmc.py
pmc_availability.jsonl에서 pmc_oa=True인 3,789건 PMC XML 다운로드 + 섹션 파싱

저장: data/fulltext/{PMCID}.json
로그: data/pmc_download.log
"""
import asyncio, json, logging
from pathlib import Path
import aiohttp
import xml.etree.ElementTree as ET

_HERE      = Path(__file__).resolve().parent.parent
AVAIL_PATH = _HERE / "data" / "pmc_availability.jsonl"
OUT_DIR    = _HERE / "data" / "fulltext"
LOG_PATH   = _HERE / "data" / "pmc_download.log"

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL       = "jongungsong@yuhs.ac"
TOOL        = "gof_pipeline"
RATE        = 0.34   # 3 req/s (NCBI without API key)
WORKERS     = 5

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(), logging.FileHandler(LOG_PATH)],
)
log = logging.getLogger(__name__)


def parse_pmc_xml(xml_text: str) -> dict:
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return {}

    def get_text(elem) -> str:
        return " ".join(elem.itertext()).strip()

    sections: dict[str, str] = {}

    # title
    for t in root.findall(".//article-title"):
        sections["title"] = get_text(t)
        break

    # abstract
    for ab in root.findall(".//abstract"):
        sections["abstract"] = get_text(ab)
        break

    # body sections
    target = {
        "introduction":           "introduction",
        "intro":                  "introduction",
        "background":             "background",
        "methods":                "methods",
        "materials and methods":  "methods",
        "methods and materials":  "methods",
        "methodology":            "methods",
        "experimental procedures":"methods",
        "experimental section":   "methods",
        "patients and methods":   "methods",
        "subjects and methods":   "methods",
        "results":                "results",
        "results and discussion": "results",
        "findings":               "results",
        "discussion":             "discussion",
        "conclusion":             "conclusion",
        "conclusions":            "conclusion",
        "summary":                "conclusion",
    }
    for sec in root.findall(".//sec"):
        title_el = sec.find("title")
        if title_el is None:
            continue
        raw   = (title_el.text or "").strip().lower()
        key   = target.get(raw)
        if not key:
            # 부분 매칭 (e.g. "2. methods", "materials & methods")
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


async def fetch_xml(session: aiohttp.ClientSession, pmcid: str) -> str | None:
    params = {
        "db": "pmc", "id": pmcid.replace("PMC", ""),
        "rettype": "full", "retmode": "xml",
        "tool": TOOL, "email": EMAIL,
    }
    for attempt in range(3):
        try:
            async with session.get(
                f"{EUTILS_BASE}/efetch.fcgi",
                params=params,
                timeout=aiohttp.ClientTimeout(total=90),
            ) as resp:
                if resp.status == 200:
                    return await resp.text()
                if resp.status == 429:
                    await asyncio.sleep(5)
        except Exception as e:
            if attempt == 2:
                log.warning(f"efetch 실패 {pmcid}: {e}")
        await asyncio.sleep(2 ** attempt)
    return None


async def run():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # PMC OA 대상 로드
    oa_recs = []
    with open(AVAIL_PATH) as f:
        for line in f:
            r = json.loads(line)
            if r.get("pmc_oa") and r.get("pmcid"):
                oa_recs.append(r)
    log.info(f"PMC OA 대상: {len(oa_recs):,}건")

    # 이미 완료된 것 스킵
    done_set = {p.stem for p in OUT_DIR.glob("*.json")}
    todo     = [r for r in oa_recs if r["pmcid"] not in done_set]
    log.info(f"이미 완료: {len(done_set):,} / 신규 다운로드: {len(todo):,}건")

    if not todo:
        log.info("모두 완료됨.")
        return

    sem   = asyncio.Semaphore(WORKERS)
    saved = 0
    errors = 0
    total = len(todo)

    connector = aiohttp.TCPConnector(limit=WORKERS + 5)
    async with aiohttp.ClientSession(connector=connector) as session:

        async def download_one(rec):
            nonlocal saved, errors
            async with sem:
                pmcid    = rec["pmcid"]
                xml_text = await fetch_xml(session, pmcid)
                if not xml_text:
                    errors += 1
                    return
                sections = parse_pmc_xml(xml_text)
                out = {
                    "pmid":          rec["pmid"],
                    "pmcid":         pmcid,
                    "doi":           rec.get("doi", ""),
                    "score_overall": rec["score_overall"],
                    "title":         rec["title"],
                    "sections":      sections,
                    "n_sections":    len(sections),
                    "chars_total":   sum(len(v) for v in sections.values()),
                }
                (OUT_DIR / f"{pmcid}.json").write_text(
                    json.dumps(out, ensure_ascii=False, indent=2)
                )
                saved += 1
                if saved % 50 == 0:
                    pct = saved / total * 100
                    log.info(f"  [{saved:,}/{total:,} {pct:.1f}%] 저장 완료 | 오류: {errors}")
                await asyncio.sleep(RATE)

        await asyncio.gather(*[download_one(r) for r in todo])

    log.info(f"\n{'='*50}")
    log.info(f"=== PMC 다운로드 완료 ===")
    log.info(f"  저장: {saved:,}건")
    log.info(f"  오류: {errors}건")
    log.info(f"  저장 경로: {OUT_DIR}/")

    # 섹션 통계
    all_files = list(OUT_DIR.glob("*.json"))
    sec_counts = {"title":0,"abstract":0,"introduction":0,
                  "methods":0,"results":0,"discussion":0,"conclusion":0}
    for fp in all_files:
        secs = json.loads(fp.read_text()).get("sections", {})
        for k in sec_counts:
            if k in secs:
                sec_counts[k] += 1
    n = len(all_files)
    log.info(f"\n=== 섹션 파싱 통계 ({n:,}건) ===")
    for sec, cnt in sec_counts.items():
        log.info(f"  {sec:15s}: {cnt:,} ({cnt/n*100:.0f}%)")


asyncio.run(run())
