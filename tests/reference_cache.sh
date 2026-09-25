#!/usr/bin/env bash
# Custom-reference index cache: a different GTF in the same outdir must build a new index,
# and STAR genomeGenerate must be given less RAM than the container limit.
# Needs Docker. Usage: bash tests/reference_cache.sh
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cd "$TMP"
run() { nextflow run "$REPO/main.nf" -profile test,docker -ansi-log false -w "$TMP/work" --outdir "$TMP/out" "$@" > "$TMP/run.log" 2>&1 || { tail -20 "$TMP/run.log"; exit 1; }; }

run
# Same genome, edited annotation (drops the last 200 GTF lines) under another name.
head -n -200 "$REPO/assets/test/genes.gtf" > "$TMP/genes_v2.gtf"
run --gtf "$TMP/genes_v2.gtf"

n_index=$(find "$TMP/out/reference" -name SA -path '*star_index*' | wc -l)
[[ $n_index -eq 2 ]] || { echo "FAIL: expected 2 STAR indexes (one per GTF), found $n_index"; exit 1; }
echo "ok: different GTF built its own index ($n_index indexes)"

# The test profile caps memory at 6 GB; genomeGenerate must get 90% of it.
expected=$(( 6442450944 / 10 * 9 ))
cmd=$(grep -l genomeGenerate "$TMP"/work/*/*/.command.sh | head -1)
grep -q -- "--limitGenomeGenerateRAM $expected" "$cmd" || { echo "FAIL: limitGenomeGenerateRAM is not $expected:"; grep -o -- '--limitGenomeGenerateRAM [0-9]*' "$cmd"; exit 1; }
echo "ok: limitGenomeGenerateRAM = $expected (90% of task memory)"
