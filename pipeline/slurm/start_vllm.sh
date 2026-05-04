#!/bin/bash
# start_vllm.sh
# Qwen3-30B-A3B-Instruct vLLM 서버 시작
#
# Usage:
#   bash start_vllm.sh          # 서버 시작
#   bash start_vllm.sh --port 8001  # 포트 변경

MODEL="Qwen/Qwen3-32B-AWQ"
PORT=${2:-8100}
LOG_FILE="/home/jihoon_yoon/jwsong/projects/gof_pipeline/logs/vllm_server.log"

mkdir -p "$(dirname "$LOG_FILE")"

echo "=== vLLM 서버 시작 ==="
echo "모델: $MODEL"
echo "포트: $PORT"
echo "로그: $LOG_FILE"
echo ""

VLLM_USE_V1=0 /home/jihoon_yoon/jwsong/miniconda3/envs/gof_pipeline/bin/vllm serve "$MODEL" \
    --port "$PORT" \
    --tensor-parallel-size 1 \
    --max-model-len 8192 \
    --gpu-memory-utilization 0.92 \
    --max-num-seqs 32 \
    --quantization awq \
    --dtype half \
    --enforce-eager \
    2>&1 | tee "$LOG_FILE"
