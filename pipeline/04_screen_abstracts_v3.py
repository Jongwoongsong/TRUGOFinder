"""
04_screen_abstracts_v3.py  — v3: stop-loss / mixed 케이스 보완
v2 대비 변경:
  - score_frameshift 기준 확장: stop-loss (C-terminal extension), nonsense/stop-gain
    causing C-terminal truncation GOF도 frameshift-equivalent로 인정 (0.6–0.8)
  - mixed 케이스: 보고된 변이 중 하나라도 frameshift/truncating 성분이 있으면
    그 성분 기준으로 score_frameshift 평가
  - false positive rule 1 완화: 순수 missense/SNV만 LOW 강제,
    truncating variants (nonsense last-exon, stop-loss, stop-gain) + GOF는 HIGH 가능
  - 출력: results/screening_scores_v3.jsonl

Usage:
  python 04_screen_abstracts_v3.py          # 전체 (resume 자동)
  python 04_screen_abstracts_v3.py --test   # 첫 30개만
"""

import argparse
import json
import logging
import os
import time
from pathlib import Path

import anthropic
import pandas as pd

# ── 설정 ──────────────────────────────────────────────────────────────────────
_HERE      = Path(__file__).resolve().parent.parent
XLSX_PATH  = Path("/home/jyoon/FSGOF/meta/FSGOF_meta_250512.xlsx")
OUT_PATH   = _HERE / "results" / "screening_scores_v3.jsonl"
MODEL      = "claude-haiku-4-5-20251001"
TEST_N     = 30
RATE_SLEEP = 0.3

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# ── Abstract_screen 로드 ───────────────────────────────────────────────────────
def load_abstracts(test: bool) -> list[dict]:
    df = pd.read_excel(XLSX_PATH, sheet_name="Abstract_screen")
    df = df[df["Abstract"].notna() & (df["Abstract"].str.strip() != "")].copy()
    df["DOI"] = df["DOI"].fillna("").str.strip()
    if test:
        df = df.head(TEST_N)
    records = []
    for _, row in df.iterrows():
        records.append({
            "no":               int(row["No."]) if pd.notna(row["No."]) else None,
            "doi":              row["DOI"],
            "source":           str(row.get("Source", "") or ""),
            "title":            str(row.get("Title",    "") or ""),
            "abstract":         str(row.get("Abstract", "") or ""),
            "abstract_include": int(row["Abstract_include"]) if pd.notna(row.get("Abstract_include")) else None,
        })
    return records

def load_done(path: Path) -> set[str]:
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                try:
                    done.add(json.loads(line)["doi"])
                except Exception:
                    pass
    return done

# ── 시스템 프롬프트 (v3) ──────────────────────────────────────────────────────
SYSTEM_PROMPT = """You are a biomedical literature classifier for a systematic review on **truncating variants (frameshift or equivalent) that cause gain-of-function (GOF) or dominant-negative (DN) effects**.

## TARGET PAPERS
Papers reporting a truncating variant — where the protein reading frame is disrupted OR the protein is C-terminally altered — leading to GOF or DN consequence.

### Qualifying variant types (score_frameshift ≥ 0.6):
  1. **Classical frameshift**: insertion/deletion/duplication of 1–few bp causing a reading frame shift (explicit "frameshift", HGVS p.Xxxfs notation, "out-of-frame")
  2. **Stop-loss / read-through**: a stop codon mutated to a sense codon, creating a C-terminal extension with novel amino acid sequence (e.g. p.Ter2648SerextTer101, "C-terminal extension", "novel C-terminal sequence via read-through")
  3. **Nonsense/stop-gain causing GOF**: a premature stop codon that removes a C-terminal inhibitory domain, resulting in constitutive activation (NOT classic LOF/NMD). Recognizable by "loss of [inhibitory/regulatory] domain", "constitutive activation", "truncated protein with gain-of-function".
  4. **Truncating indel**: small in-frame or out-of-frame deletion/insertion removing a C-terminal regulatory region and causing GOF.

### NON-qualifying variant types (score_frameshift ≤ 0.2):
  - Missense / amino acid substitution / point mutation / SNV — even if they cause GOF
  - Large chromosomal deletions, CNV, gene fusions
  - Splice variants (unless they cause frameshift or novel exon inclusion)

## SCORING RULES

### score_frameshift  (0.0–1.0)
  - 0.9–1.0 : explicit "frameshift", HGVS fs notation, "out-of-frame indel", or last-exon frameshift
  - 0.7–0.85 : stop-loss / read-through creating novel C-terminal sequence; OR nonsense variant explicitly described as causing GOF via C-terminal domain loss
  - 0.6–0.7  : "truncating variant/mutation" causing C-terminal loss + GOF/constitutive activation (mechanism matches even if variant type ambiguous)
  - 0.3–0.5  : unclear variant type; may be frameshift/nonsense/splice
  - 0.0–0.2  : missense, SNV, large deletion, CNV — confirmed non-truncating

  **For mixed papers**: evaluate score_frameshift based on the truncating/frameshift component,
  not the missense component. If the paper reports BOTH a frameshift and missense variants,
  and the frameshift variant is part of the primary finding, score_frameshift accordingly.

### score_gof  (0.0–1.0)
  - 0.8–1.0 : "gain-of-function", "constitutive activation", "NMD escape", "novel peptide/sequence",
               "C-terminal truncation → activation", "loss of inhibitory domain", "increased stability/expression"
  - 0.5–0.7 : GOF strongly implied (receptor hyperactivation, dominant overgrowth syndrome, ligand-independent signaling)
  - 0.2–0.4 : GOF present but for a different variant than the truncating one; or mixed GOF+LOF
  - 0.0–0.1 : pure LOF / haploinsufficiency / no GOF evidence

### score_dominant_negative  (0.0–1.0)
  - 0.8–1.0 : "dominant negative", "dominant-negative", "poison subunit", "interferes with wild-type"
  - 0.5–0.7 : DN implied by heterozygous truncation disrupting an oligomeric complex
  - 0.0–0.2 : no DN evidence

### score_overall  (0.0–1.0)  ← MOST IMPORTANT
**Rule**: score_overall reflects co-occurrence of a qualifying truncating variant AND (GOF or DN).
  - ≥ 0.7 (HIGH) : score_frameshift ≥ 0.6 AND (score_gof ≥ 0.7 OR score_dominant_negative ≥ 0.7)
  - 0.4–0.69 (MID): truncating variant present but GOF/DN weak; OR GOF/DN present but variant type uncertain
  - < 0.4 (LOW)  : variant is missense/SNV only; OR neither GOF nor DN present

**FALSE POSITIVE RULES:**
  1. Pure missense/SNV GOF (no truncation component) → score_overall < 0.4, score_frameshift ≤ 0.1.
  2. "gain-of-function" in title for a confirmed missense variant → score_overall ≤ 0.35.
  3. Frameshift/truncation mentioned only as background/comparison, not as the primary finding → score_frameshift ≤ 0.25.
  4. Truncating variant causing pure LOF/NMD with no GOF/DN consequence → score_overall < 0.4.

**EXAMPLES:**
  - p.Ter2648SerextTer101 (stop-loss, 101aa novel C-terminal extension) + platelet GOF → score_frameshift=0.8, score_gof=0.9, score_overall=0.85 ✓
  - 6 TID-truncating variants causing constitutive TP63 activation + POI → score_frameshift=0.7, score_gof=0.9, score_overall=0.80 ✓
  - p.Arg175His missense in TP53 causing GOF (no truncation) → score_frameshift=0.05, score_overall=0.2 ✗
  - STAT1 p.Thr385Met missense GOF → score_frameshift=0.05, score_overall=0.15 ✗
  - DVL3 frameshift in last exon → Robinow GOF → score_frameshift=0.95, score_gof=0.85, score_overall=0.90 ✓

### variant_type
  "frameshift"   — reading-frame-disrupting indel (explicit or HGVS fs)
  "stop_loss"    — stop codon mutated to sense, C-terminal extension
  "truncating"   — C-terminal truncating variant (nonsense/stop-gain causing GOF via domain loss)
  "missense"     — amino acid substitution only
  "nonsense_lof" — premature stop → LOF/NMD (not GOF)
  "splice"       — splice site variant
  "mixed"        — paper reports frameshift/truncating AND missense variants as primary findings
  "other"        — CNV, fusion, chromosomal
  "unclear"      — cannot determine

## OUTPUT FORMAT
Return ONLY a valid JSON object, no extra text:
{
  "score_frameshift":        <float 0–1>,
  "score_gof":               <float 0–1>,
  "score_dominant_negative": <float 0–1>,
  "score_overall":           <float 0–1>,
  "variant_type":            <string from list above>,
  "matched_concepts":        [<list of matched concept strings>],
  "reasoning":               "<2–4 sentences: what variant type, what GOF/DN mechanism, why score_overall is this value>"
}"""

USER_PROMPT_TEMPLATE = "Title: {title}\n\nAbstract:\n{abstract}"

# ── LLM 호출 ─────────────────────────────────────────────────────────────────
def call_llm(client: anthropic.Anthropic, title: str, abstract: str) -> dict:
    msg = client.messages.create(
        model=MODEL,
        max_tokens=512,
        system=SYSTEM_PROMPT,
        messages=[{"role": "user",
                   "content": USER_PROMPT_TEMPLATE.format(title=title, abstract=abstract[:3000])}]
    )
    raw = msg.content[0].text.strip()
    if raw.startswith("```"):
        raw = raw.split("```")[1]
        if raw.startswith("json"):
            raw = raw[4:]
    return json.loads(raw.strip())

# ── 메인 ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store_true", help=f"첫 {TEST_N}개만 처리")
    args = parser.parse_args()

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise RuntimeError("ANTHROPIC_API_KEY 환경변수 미설정")
    client = anthropic.Anthropic(api_key=api_key)

    records = load_abstracts(test=args.test)
    done    = load_done(OUT_PATH)

    todo = [r for r in records if not (r["doi"] and r["doi"] in done)]
    log.info(f"대상: {len(records)}개 | 완료(스킵): {len(done)}개 | 신규: {len(todo)}개")

    OUT_PATH.parent.mkdir(exist_ok=True)
    processed = errors = 0

    with OUT_PATH.open("a", encoding="utf-8") as out_f:
        for i, rec in enumerate(todo, 1):
            try:
                result = call_llm(client, rec["title"], rec["abstract"])
                output = {
                    "no":               rec["no"],
                    "doi":              rec["doi"],
                    "source":           rec["source"],
                    "title":            rec["title"],
                    "abstract_include": rec["abstract_include"],
                    "score_frameshift":        result.get("score_frameshift"),
                    "score_gof":               result.get("score_gof"),
                    "score_dominant_negative": result.get("score_dominant_negative"),
                    "score_overall":           result.get("score_overall"),
                    "variant_type":            result.get("variant_type", "unclear"),
                    "matched_concepts":        result.get("matched_concepts", []),
                    "reasoning":               result.get("reasoning", ""),
                }
                out_f.write(json.dumps(output, ensure_ascii=False) + "\n")
                out_f.flush()

                s = output["score_overall"]
                vt = output["variant_type"]
                tag = "★ HIGH" if s >= 0.7 else ("△ MID" if s >= 0.4 else "  low")
                log.info(f"[{i:4d}/{len(todo)}] {tag} {s:.2f} [{vt:<11}] {rec['title'][:55]}")
                processed += 1
                time.sleep(RATE_SLEEP)

            except json.JSONDecodeError as e:
                log.warning(f"  JSON 파싱 실패 ({rec['doi']}): {e}")
                errors += 1
            except Exception as e:
                log.error(f"  API 오류 ({rec['doi']}): {e}")
                errors += 1
                time.sleep(2)

    log.info(f"\n완료: 처리 {processed}개 / 스킵 {len(done)}개 / 오류 {errors}개")
    log.info(f"결과: {OUT_PATH}")

    # 최종 분포
    all_scores = []
    vtype_counts = {}
    for line in OUT_PATH.read_text().splitlines():
        if line.strip():
            try:
                r = json.loads(line)
                s = r.get("score_overall")
                if s is not None:
                    all_scores.append(s)
                vt = r.get("variant_type", "unclear")
                vtype_counts[vt] = vtype_counts.get(vt, 0) + 1
            except Exception:
                pass

    if all_scores:
        high = sum(1 for s in all_scores if s >= 0.7)
        mid  = sum(1 for s in all_scores if 0.4 <= s < 0.7)
        low  = sum(1 for s in all_scores if s < 0.4)
        log.info(f"\n=== 스코어 분포 ({len(all_scores)}개) ===")
        log.info(f"  HIGH (≥0.7): {high}개")
        log.info(f"  MID  (0.4–): {mid}개")
        log.info(f"  LOW  (<0.4): {low}개")
        log.info(f"  평균 score:  {sum(all_scores)/len(all_scores):.3f}")
        log.info(f"\n=== variant_type 분포 ===")
        for vt, cnt in sorted(vtype_counts.items(), key=lambda x: -x[1]):
            log.info(f"  {vt:<12}: {cnt}개")

if __name__ == "__main__":
    main()
