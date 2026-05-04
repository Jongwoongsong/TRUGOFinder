#!/bin/bash
#SBATCH --job-name=expand_extract
#SBATCH --partition=debug
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/home/jihoon_yoon/jwsong/logs/slurm_expand_extract_%j.log
#SBATCH --error=/home/jihoon_yoon/jwsong/logs/slurm_expand_extract_%j.log

# ─────────────────────────────────────────────────────────────────────────────
# 52_expand_extract.sbatch
# 확장 fulltext 추출 — score < 0.7 논문 (51_expand_download.py 다운로드 후 실행)
#
# 51_expand_download.py가 data/fulltext/ 에 새 파일을 추가한 후 이 job을 실행.
# 11_extract_variants.py는 이미 추출된 pmcid를 자동으로 스킵하므로
# 6시간 내에 처리 못한 경우 동일 sbatch를 재제출하면 이어받기 가능.
#
# Input:  data/fulltext/{PMCID}.json  (51 스크립트가 추가한 파일 포함)
# Output: data/extracted_variants_expanded_v1.jsonl
# ─────────────────────────────────────────────────────────────────────────────

set -e
cd /home/jihoon_yoon/jwsong/projects/gof_pipeline
source /home/jihoon_yoon/jwsong/miniconda3/etc/profile.d/conda.sh
conda activate gof_pipeline
PYTHON=python3

export VLLM_BASE=http://localhost:8005/v1
OUT_JSONL=data/extracted_variants_expanded_v1.jsonl

echo "=== 확장 LLM 추출 (52) ==="
echo "대상: data/fulltext/ 전체 중 expanded_v1에 미추출"
date

# ── vLLM 시작 ─────────────────────────────────────────────────────────────────
echo ""
echo "[vLLM] port 8005 시작..."
fuser -k 8005/tcp 2>/dev/null || true
sleep 3

VLLM_USE_V1=0 /home/jihoon_yoon/jwsong/miniconda3/envs/gof_pipeline/bin/vllm serve \
    Qwen/Qwen3-32B-AWQ \
    --quantization awq \
    --max-model-len 12288 \
    --gpu-memory-utilization 0.92 \
    --port 8005 \
    --dtype half \
    --enforce-eager \
    > /home/jihoon_yoon/jwsong/logs/vllm_expand_extract_${SLURM_JOB_ID}.log 2>&1 &

VLLM_PID=$!
echo "  vLLM PID: $VLLM_PID"

echo -n "  vLLM 대기"
for i in $(seq 1 120); do
    sleep 5
    MODEL_CHECK=$(curl -s http://localhost:8005/v1/models 2>/dev/null || echo "")
    if echo "$MODEL_CHECK" | python3 -c "import sys,json; d=json.load(sys.stdin); exit(0 if d.get('data') else 1)" 2>/dev/null; then
        echo " OK (${i}×5s)"
        break
    fi
    echo -n "."
done

# ── 현황 확인 ─────────────────────────────────────────────────────────────────
echo ""
echo "=== 추출 전 현황 ==="
$PYTHON - <<'PYEOF'
import json
from pathlib import Path

base   = Path(".")
ft_dir = base / "data/fulltext"
out_v9 = base / "data/extracted_variants_v9.jsonl"
out_ex = base / "data/extracted_variants_expanded_v1.jsonl"

# 기존 추출 완료 pmcid 수집 (v9 + expanded)
done_pmcids = set()
for out_path in [out_v9, out_ex]:
    if out_path.exists():
        with open(out_path) as f:
            for line in f:
                d = json.loads(line)
                done_pmcids.add(d.get("pmcid", ""))

all_ft        = list(ft_dir.glob("*.json"))
not_extracted = [p for p in all_ft if p.stem not in done_pmcids]
already_done  = len(all_ft) - len(not_extracted)

print(f"  전체 fulltext:    {len(all_ft):,}건")
print(f"  기추출 (v9+exp):  {already_done:,}건")
print(f"  신규 추출 대상:   {len(not_extracted):,}건")

# score 분포 확인
try:
    low_score = 0
    for p in not_extracted:
        d = json.loads(p.read_text())
        if float(d.get("score_overall", 1.0)) < 0.7:
            low_score += 1
    print(f"  그 중 score<0.7: {low_score:,}건 (신규 확장분)")
except Exception:
    pass
PYEOF

# ── LLM 추출 실행 ─────────────────────────────────────────────────────────────
echo ""
echo "=== LLM 추출 시작 ==="

# 11_extract_variants.py는 --output 파일에 이미 있는 pmcid를 자동 스킵.
# v9에 있는 것도 스킵해야 하므로, v9 결과를 expanded에 soft-link하지 않고
# 별도 --skip-done-from 옵션이 없으므로 아래 방식으로 처리:
# expanded_v1에 v9의 pmcid를 미리 기록 (빈 findings로) → 스킵 유도
$PYTHON - <<'PYEOF'
import json
from pathlib import Path

base   = Path(".")
out_v9 = base / "data/extracted_variants_v9.jsonl"
out_ex = base / "data/extracted_variants_expanded_v1.jsonl"

# expanded_v1에 이미 있는 pmcid 수집
existing_ex = set()
if out_ex.exists():
    with open(out_ex) as f:
        for line in f:
            d = json.loads(line)
            existing_ex.add(d.get("pmcid", ""))

# v9에 있는 항목을 expanded_v1에 stub로 추가 (아직 없는 것만)
added = 0
with open(out_ex, "a", encoding="utf-8") as out_f:
    if out_v9.exists():
        with open(out_v9) as f:
            for line in f:
                d = json.loads(line)
                pmcid = d.get("pmcid", "")
                if pmcid and pmcid not in existing_ex:
                    # stub: pmcid 기록만, findings 보존
                    out_f.write(json.dumps({
                        "pmcid": pmcid,
                        "pmid":  d.get("pmid",""),
                        "findings": d.get("findings", []),
                        "_source": "v9"
                    }, ensure_ascii=False) + "\n")
                    existing_ex.add(pmcid)
                    added += 1

print(f"  v9 stubs added to expanded_v1: {added:,}")
PYEOF

$PYTHON pipeline/11_extract_variants.py \
    --output      data/extracted_variants_expanded_v1.jsonl \
    --log         logs/expand_extract_${SLURM_JOB_ID}.log \
    --concurrency 4 \
    --skip-mutalyzer

# ── 결과 확인 ─────────────────────────────────────────────────────────────────
echo ""
echo "=== 추출 완료 현황 ==="
$PYTHON - <<'PYEOF'
import json
from pathlib import Path
from collections import Counter

base    = Path(".")
out_ex  = base / "data/extracted_variants_expanded_v1.jsonl"
out_v9  = base / "data/extracted_variants_v9.jsonl"

# expanded_v1 통계 (v9 stubs 제외)
new_papers = 0
new_findings = 0
is_human_cnt = Counter()

if out_ex.exists():
    with open(out_ex) as f:
        for line in f:
            d = json.loads(line)
            if d.get("_source") == "v9":
                continue   # v9 stub 제외
            new_papers += 1
            findings = d.get("findings", [])
            new_findings += len(findings)
            for fnd in findings:
                is_human_cnt[str(fnd.get("is_human"))] += 1

print(f"  expanded_v1 신규 논문: {new_papers:,}건")
print(f"  신규 findings:         {new_findings:,}건")
print(f"  is_human 분포:         {dict(is_human_cnt)}")

# 잔여 추출 대상
ft_dir = base / "data/fulltext"
done_pmcids = set()
with open(out_ex) as f:
    for line in f:
        done_pmcids.add(json.loads(line).get("pmcid",""))
remaining = sum(1 for p in ft_dir.glob("*.json") if p.stem not in done_pmcids)
print(f"  잔여 추출 대상:        {remaining:,}건")
if remaining > 0:
    print("  → sbatch 52_expand_extract.sbatch 재제출 필요")
PYEOF

# ── vLLM 종료 ─────────────────────────────────────────────────────────────────
echo ""
echo "[vLLM] 종료..."
kill $VLLM_PID 2>/dev/null || true

echo ""
echo "=== 완료 ==="
echo "다음 단계: 53_expand_merge.py 실행 (새 findings 통합)"
date
