"""
로컬 OMIM DB 조회 모듈
- /home/jyoon/FSGOF/OMIM_Gene_Pheno_database_2401.txt 사용
- Gene + Disease(선택) → OMIM ID 반환
- Entrez API 불필요, 즉시 조회

사용법:
    from omim_local import lookup_omim_local
    omim_id = lookup_omim_local("CXCR4", "WHIM syndrome")
"""

import re
import pandas as pd

OMIM_DB_PATH = "/home/jyoon/FSGOF/OMIM_Gene_Pheno_database_2401.txt"

# 모듈 로드 시 한 번만 읽음
_db: dict[str, list[tuple[str, str]]] = {}   # Gene → [(pheno, omim_id), ...]


def _load_db(path: str = OMIM_DB_PATH) -> None:
    global _db
    if _db:
        return
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    for _, row in df.iterrows():
        gene   = row["Gene"].strip()
        phenos = [p.strip() for p in row["OMIM_Pheno"].split(";")]
        ids    = [i.strip() for i in row["OMIM_ID"].split(";")]
        # 길이 맞추기
        max_len = max(len(phenos), len(ids))
        phenos += [""] * (max_len - len(phenos))
        ids    += [""] * (max_len - len(ids))
        _db.setdefault(gene, []).extend(zip(phenos, ids))


def _normalize(s: str) -> str:
    return re.sub(r"\s+", " ", s.lower().strip())


def _is_valid_id(omim_id: str) -> bool:
    """숫자 5-6자리이고 공백/빈값이 아닌 경우만 유효."""
    return bool(re.match(r"^\d{5,6}$", omim_id.strip()))


def _pheno_score(pheno: str, disease: str) -> float:
    """disease 단어가 pheno에 얼마나 겹치는지 (0~1)."""
    if not disease:
        return 0.0
    # {susceptibility} 계열은 패널티
    if pheno.startswith("{") or pheno.startswith("?"):
        return -1.0
    g_words = set(_normalize(disease).split())
    p_words = set(_normalize(pheno).split())
    if not g_words:
        return 0.0
    overlap = g_words & p_words
    return len(overlap) / len(g_words)


def lookup_omim_local(gene: str, disease: str = "") -> tuple[str | None, str]:
    """Gene(+Disease) → (OMIM ID | None, 선택 근거) 반환.

    근거 문자열:
      "exact"         - disease 완전 일치
      "partial:{N}"   - N개 단어 겹침
      "first_valid"   - disease 불일치, 유효한 첫 번째 ID
      "not_found"     - DB에 Gene 없음
    """
    _load_db()
    gene = gene.strip()

    # Gene 없는 경우
    if gene not in _db:
        return None, "not_found"

    entries = _db[gene]   # [(pheno, omim_id), ...]

    # 유효한 ID만 필터
    valid = [(p, i) for p, i in entries if _is_valid_id(i)]
    if not valid:
        return None, "no_valid_id"

    # Disease 매칭 시도
    if disease:
        scored = [(p, i, _pheno_score(p, disease)) for p, i in valid]
        scored.sort(key=lambda x: -x[2])
        best_pheno, best_id, best_score = scored[0]

        if best_score >= 0.5:
            tag = "exact" if best_score >= 0.99 else f"partial:{best_score:.2f}"
            return best_id, tag

    # Disease 매칭 실패 → 첫 번째 유효 ID (susceptibility 제외 우선)
    non_susceptibility = [(p, i) for p, i in valid
                          if not p.startswith("{") and not p.startswith("?")]
    chosen = non_susceptibility[0] if non_susceptibility else valid[0]
    return chosen[1], "first_valid"
