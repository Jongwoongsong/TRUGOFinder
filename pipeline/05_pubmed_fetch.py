"""
05_pubmed_fetch.py
PubMed E-utilities API로 키워드 검색 후 abstract 배치 다운로드.

검색 전략:
  - keyword_list.md 기반 3개 개념 그룹 (Frameshift / GOF / DN)
  - 1단계: (Frameshift 그룹) AND (GOF 그룹 OR DN 그룹)  → 핵심 민감도
  - 2단계: 고특이도 메커니즘 표현 단독 (stop-loss, NMD escape, penultimate exon 등)
  - 전체 OR 조합
  - 필터: humans[MeSH], English[lang], 1995:2025[dp]

Usage:
  python 05_pubmed_fetch.py               # 전체 실행
  python 05_pubmed_fetch.py --dry-run     # 검색식만 출력, 히트 수 확인
  python 05_pubmed_fetch.py --resume      # 이미 받은 PMID 스킵
"""

import argparse
import json
import logging
import os
import time
from pathlib import Path
from typing import Iterator

from Bio import Entrez

# ── 설정 ──────────────────────────────────────────────────────────────────────
_HERE       = Path(__file__).resolve().parent.parent
OUT_PATH    = _HERE / "data" / "pubmed_abstracts.jsonl"
QUERY_LOG   = _HERE / "data" / "pubmed_query.txt"

Entrez.email  = "jongungsong@yuhs.ac"
BATCH_SIZE    = 500      # efetch 한 번에 가져올 PMID 수
RATE_SLEEP    = 0.35     # API 호출 간격 (초) — API key 없이 3req/s 허용
MAX_RECORDS   = 100_000  # 안전 상한선

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# ── 검색식 구성 ────────────────────────────────────────────────────────────────
# keyword_list.md 1-A~1-F / 2-A~2-F / 3-A~3-H 전체 반영
# 전략: 전체 OR (sensitivity 최우선) → LLM 스크리닝에서 GOF/DN 판단
# PubMed [tiab] = title + abstract 필드
#
# ※ 제외 이유 명시
#   - 1-C (small insertion/deletion, microindel): 너무 광범위, specificity 0
#   - 1-D HGVS 패턴: PubMed tiab 검색 불가 (p.Xfs*N 등 정규식 미지원)
#   - 1-F "indel" 단독: 2만+ 건, frameshift와 직접 연관 없는 경우 다수
#   - 3-C "enhanced activity", "increased activity" 등: 너무 광범위
#   - 2-E "acts as a decoy": 비특이적

# ── 1. FRAMESHIFT ─────────────────────────────────────────────────────────────

# 1-A: 직접 명칭
FS_DIRECT = [
    'frameshift[tiab]',                          # 가장 넓은 net (단독)
    '"frameshift change"[tiab]',
    '"coding frameshift"[tiab]',
]

# 1-B: Reading frame disruption 묘사
FS_MECH = [
    '"out-of-frame"[tiab]',
    '"reading frame disruption"[tiab]',
    '"reading frame alteration"[tiab]',
    '"reading frame shift"[tiab]',
    '"disruption of the reading frame"[tiab]',
    '"alters the reading frame"[tiab]',
    '"shifts the reading frame"[tiab]',
    '"changes the open reading frame"[tiab]',
]

# 1-E: 결과 표현 (PTC, truncation, novel peptide 등)
FS_CONSEQUENCE = [
    '"premature termination codon"[tiab]',
    '"premature stop codon"[tiab]',
    '"truncating variant"[tiab]',
    '"truncating mutation"[tiab]',
    '"truncating mutations"[tiab]',
    '"C-terminal truncation"[tiab]',
    '"C-terminally truncated"[tiab]',
    '"truncated C terminus"[tiab]',
    'neopeptide[tiab]',
    '"novel peptide"[tiab]',
    '"cryptic ORF"[tiab]',
    '"novel open reading frame"[tiab]',
]

# 3-F: Frameshift 특이적 GOF 표현
FS_GOF_SPECIFIC = [
    '"penultimate exon"[tiab]',
    '"last exon frameshift"[tiab]',
    '"frameshift extension"[tiab]',
    '"frameshift-extended"[tiab]',
    '"frameshift gain-of-function"[tiab]',
    '"gain-of-function frameshift"[tiab]',
    '"frameshift-induced gain-of-function"[tiab]',
    '"activating frameshift"[tiab]',
    '"activating indel"[tiab]',
    '"oncogenic frameshift"[tiab]',
    '"gain-of-function truncation"[tiab]',
    '"truncation-mediated activation"[tiab]',
    '"truncated protein with gain-of-function"[tiab]',
]

# 3-F-2: Stop-loss / Read-through GOF
STOPLOSS = [
    '"stop-loss mutation"[tiab]',
    '"stop-loss variant"[tiab]',
    '"stop codon read-through"[tiab]',
    '"read-through translation"[tiab]',
    '"C-terminal extension"[tiab]',
    '"novel C-terminal sequence"[tiab]',
    '"novel C-terminal peptide"[tiab]',
    '"protein elongation"[tiab]',
]

# ── 2. DOMINANT NEGATIVE ──────────────────────────────────────────────────────

# 2-A: 직접 명칭
DN_DIRECT = [
    '"dominant negative"[tiab]',
    '"dominant-negative"[tiab]',
    '"dominant negative variant"[tiab]',
    '"dominant negative allele"[tiab]',
    '"dominant negative activity"[tiab]',
    '"dominant negative function"[tiab]',
    '"dominant negative mechanism"[tiab]',
    '"dominant-negative manner"[tiab]',
    '"dominant-negative interference"[tiab]',
    '"acts in a dominant negative manner"[tiab]',
]

# 2-B: Muller 용어
DN_MULLER = [
    'antimorph[tiab]',
    '"antimorphic allele"[tiab]',
    '"antimorphic mutation"[tiab]',
    '"dominant antimorphic"[tiab]',
]

# 2-C/D: 기계적 묘사 + 다량체
DN_MECH = [
    '"interferes with wild-type"[tiab]',
    '"sequesters wild-type"[tiab]',
    '"trans-dominant inhibition"[tiab]',
    '"trans-dominant repression"[tiab]',
    '"poison subunit"[tiab]',
    '"poison polypeptide"[tiab]',
    '"poison peptide"[tiab]',
    '"dominant negative heterodimer"[tiab]',
    '"nonfunctional heterodimer"[tiab]',
    '"dominant negative transcription factor"[tiab]',
]

# ── 3. GAIN-OF-FUNCTION ───────────────────────────────────────────────────────

# 3-A: 직접 명칭
GOF_DIRECT = [
    '"gain-of-function"[tiab]',
    '"gain of function"[tiab]',
    '"gain-of-function mutation"[tiab]',
    '"gain-of-function variant"[tiab]',
    '"gain-of-function allele"[tiab]',
    '"gain-of-function mechanism"[tiab]',
    'GOF[tiab]',
    'GoF[tiab]',
]

# 3-B: Muller GOF 용어
GOF_MULLER = [
    'hypermorph[tiab]',
    'neomorph[tiab]',
    '"hypermorphic allele"[tiab]',
    '"neomorphic allele"[tiab]',
    '"hypermorphic mutation"[tiab]',
    '"neomorphic mutation"[tiab]',
]

# 3-C: 활성화 표현
GOF_ACTIVATION = [
    '"constitutive activation"[tiab]',
    '"constitutively active"[tiab]',
    '"constitutive signaling"[tiab]',
    '"constitutive receptor activation"[tiab]',
    '"constitutively activating"[tiab]',
    '"activating mutation"[tiab]',
    '"activating mutations"[tiab]',
    '"activating variant"[tiab]',
    '"activating allele"[tiab]',
    '"ligand-independent activation"[tiab]',
    '"ligand-independent signaling"[tiab]',
    '"loss of autoinhibition"[tiab]',
    '"loss of inhibitory domain"[tiab]',
    '"loss of autoinhibitory domain"[tiab]',
    '"removal of inhibitory domain"[tiab]',
    '"loss of C-terminal regulatory domain"[tiab]',
    '"deletion of autoinhibitory domain"[tiab]',
]

# 3-D: Neomorph 표현
GOF_NEOMORPH = [
    '"toxic gain-of-function"[tiab]',
    '"dominant gain-of-function"[tiab]',
    'neofunction[tiab]',
    '"aberrant activity"[tiab]',
    '"aberrant function"[tiab]',
    '"toxic aggregation"[tiab]',
    '"cryptic exon"[tiab]',
]

# 3-E: NMD escape / 안정성
GOF_STABILITY = [
    '"NMD escape"[tiab]',
    '"NMD-resistant"[tiab]',
    '"escape from NMD"[tiab]',
    '"escapes NMD"[tiab]',
    '"escapes nonsense-mediated"[tiab]',
    '"nonsense-mediated decay escape"[tiab]',
    '"stabilized mRNA"[tiab]',
    '"stabilized transcript"[tiab]',
]

# 3-G: 임상/생물학적 결과
GOF_CLINICAL = [
    '"constitutive receptor signaling"[tiab]',
    '"receptor hyperactivation"[tiab]',
    '"loss of inhibitory domain"[tiab]',
    '"aberrant C-terminal function"[tiab]',
    '"loss of C-terminal regulation"[tiab]',
]

# ── v7 확장: 누락 gold 포착을 위한 추가 키워드 ────────────────────────────────
# 대상: NPM1(mislocalization), TEK(autoinhibit), CCND3(stabilization/MeSH)
#       + 전반적 커버리지 확대

# 4-A: 단백질 미스로컬라이제이션 (NPM1 등)
GOF_MISLOCAL = [
    '"cytoplasmic mislocalization"[tiab]',
    '"nuclear mislocalization"[tiab]',
    '"aberrant localization"[tiab]',
    '"protein mislocalization"[tiab]',
    '"mislocalized protein"[tiab]',
    '"cytoplasmic relocalization"[tiab]',
    '"aberrant nuclear export"[tiab]',
]

# 4-B: 자가억제 해소 / 키나아제 활성 증가 (TEK 등)
GOF_AUTOINHIBIT = [
    '"autoinhibitory mechanism"[tiab]',
    '"relieved autoinhibition"[tiab]',
    '"loss of autoinhibitory"[tiab]',
    '"enhanced kinase activity"[tiab]',
    '"kinase hyperactivation"[tiab]',
    '"constitutive kinase activity"[tiab]',
    '"kinase activation by deletion"[tiab]',
]

# 4-C: 단백질 안정화 (CCND3 등)
GOF_PROTEIN_STAB = [
    '"protein stabilization"[tiab]',
    '"stabilized mutant"[tiab]',
    '"resistance to proteasomal"[tiab]',
    '"escape from proteasomal"[tiab]',
    '"proteasome-resistant"[tiab]',
    '"ubiquitin-resistant"[tiab]',
    '"increased protein stability"[tiab]',
]

# 4-D: Tier 2 — 종양 활성화 표현
GOF_ONCOGENIC = [
    '"oncogenic mutation"[tiab]',
    '"somatic activating mutation"[tiab]',
    '"driver mutation"[tiab]',
    '"oncogenic variant"[tiab]',
]

# 4-E: MeSH term (2015년 이후 논문에 체계적으로 적용됨)
GOF_MESH = [
    '"Gain of Function Mutation"[Mesh]',
]

# ── 검색식 조합 ────────────────────────────────────────────────────────────────
def build_query() -> str:
    def OR(terms):
        return "(" + " OR ".join(terms) + ")"

    all_terms = (
        FS_DIRECT + FS_MECH + FS_CONSEQUENCE + FS_GOF_SPECIFIC + STOPLOSS +
        DN_DIRECT + DN_MULLER + DN_MECH +
        GOF_DIRECT + GOF_MULLER + GOF_ACTIVATION + GOF_NEOMORPH +
        GOF_STABILITY + GOF_CLINICAL +
        # v7 확장
        GOF_MISLOCAL + GOF_AUTOINHIBIT + GOF_PROTEIN_STAB +
        GOF_ONCOGENIC + GOF_MESH
    )

    concept_query = OR(all_terms)

    filters = (
        'English[Language]'
        ' AND 1995:2025[dp]'
    )

    return f"{concept_query} AND {filters}"


# ── PMID 검색 ─────────────────────────────────────────────────────────────────
def search_pmids(query: str) -> list[str]:
    # PubMed API retstart 상한 = 9999 (esearch/efetch 모두 해당)
    # → 연도별 분할로 우회: 연당 평균 ~3,500건, 9999 이내 보장
    base = query.replace(" AND 1995:2025[dp]", "")

    all_pmids: set[str] = set()
    log.info("연도별 PMID 수집 (1995–2025)...")
    for year in range(1995, 2026):
        year_q = f"{base} AND {year}:{year}[dp]"
        h = Entrez.esearch(db="pubmed", term=year_q, retmax=9999)
        r = Entrez.read(h); h.close()
        count = int(r["Count"])
        all_pmids.update(r["IdList"])
        log.info(f"  {year}: {count:,}건 (누적 {len(all_pmids):,})")
        if count >= 9999:
            log.warning(f"  ⚠ {year}년 {count}건 ≥ 9999 — 초과분 누락")
        time.sleep(RATE_SLEEP)

    log.info(f"총 unique PMID: {len(all_pmids):,}건")
    return list(all_pmids)


# ── Abstract 배치 수집 ────────────────────────────────────────────────────────
def fetch_batch(pmids: list[str]) -> list[dict]:
    ids = ",".join(pmids)
    handle = Entrez.efetch(db="pubmed", id=ids,
                            rettype="xml", retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    results = []
    for article in records.get("PubmedArticle", []):
        try:
            med = article["MedlineCitation"]
            art = med["Article"]

            pmid = str(med["PMID"])
            title = str(art.get("ArticleTitle", ""))

            # Abstract
            abstract_texts = art.get("Abstract", {}).get("AbstractText", [])
            if isinstance(abstract_texts, list):
                abstract = " ".join(str(t) for t in abstract_texts)
            else:
                abstract = str(abstract_texts)

            # Publication year
            pub_date = art.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
            year = str(pub_date.get("Year", pub_date.get("MedlineDate", "")))[:4]

            # DOI
            doi = ""
            for loc in article.get("PubmedData", {}).get("ArticleIdList", []):
                if str(loc.attributes.get("IdType", "")) == "doi":
                    doi = str(loc)
                    break

            # MeSH terms
            mesh = [str(m["DescriptorName"])
                    for m in med.get("MeshHeadingList", [])]

            results.append({
                "pmid":     pmid,
                "doi":      doi,
                "year":     year,
                "title":    title,
                "abstract": abstract,
                "mesh":     mesh,
            })
        except Exception as e:
            log.debug(f"파싱 오류: {e}")

    return results


# ── 이미 처리된 PMID 로드 ────────────────────────────────────────────────────
def load_done(path: Path) -> set[str]:
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                try:
                    done.add(json.loads(line)["pmid"])
                except Exception:
                    pass
    return done


# ── 메인 ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true",
                        help="검색식 출력 + 히트 수만 확인, 다운로드 안 함")
    parser.add_argument("--resume", action="store_true", default=True,
                        help="이미 받은 PMID 스킵 (기본 True)")
    args = parser.parse_args()

    OUT_PATH.parent.mkdir(exist_ok=True)

    query = build_query()

    # 쿼리 저장
    QUERY_LOG.write_text(query, encoding="utf-8")
    log.info(f"검색식 저장: {QUERY_LOG}")
    log.info(f"검색식 길이: {len(query)} chars")

    if args.dry_run:
        print("\n" + "=" * 60)
        print("PUBMED 검색식:")
        print("=" * 60)
        # 읽기 쉽게 줄바꿈
        readable = query.replace(" AND humans", "\nAND humans")
        print(readable[:3000])
        print("\n" + "=" * 60)
        # 히트 수 확인
        log.info("dry-run: 히트 수 조회 중...")
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        total = int(record["Count"])
        print(f"\n총 히트: {total:,}건")
        print(f"필터: English + 1995-2025")
        return

    # PMID 검색
    pmids = search_pmids(query)
    log.info(f"수집할 PMID: {len(pmids)}개")

    # Resume: 이미 받은 것 제외
    done = load_done(OUT_PATH) if args.resume else set()
    todo = [p for p in pmids if p not in done]
    log.info(f"기완료: {len(done)}개 / 신규: {len(todo)}개")

    if not todo:
        log.info("모두 처리 완료.")
        return

    # 배치 efetch
    saved = 0
    errors = 0
    with OUT_PATH.open("a", encoding="utf-8") as out_f:
        for batch_start in range(0, len(todo), BATCH_SIZE):
            batch = todo[batch_start:batch_start + BATCH_SIZE]
            try:
                records = fetch_batch(batch)
                for rec in records:
                    out_f.write(json.dumps(rec, ensure_ascii=False) + "\n")
                out_f.flush()
                saved += len(records)
                pct = (batch_start + len(batch)) / len(todo) * 100
                log.info(f"  [{pct:5.1f}%] 누적 저장: {saved}건 "
                         f"(배치 {batch_start//BATCH_SIZE + 1}/"
                         f"{(len(todo)-1)//BATCH_SIZE + 1})")
                time.sleep(RATE_SLEEP)
            except Exception as e:
                log.error(f"  배치 오류 (start={batch_start}): {e}")
                errors += 1
                time.sleep(3)

    log.info(f"\n완료: 저장 {saved}건 / 오류 {errors}배치")
    log.info(f"결과: {OUT_PATH}")

    # 요약
    total_saved = sum(1 for l in OUT_PATH.read_text().splitlines() if l.strip())
    log.info(f"누적 총 레코드: {total_saved}건")


if __name__ == "__main__":
    main()
