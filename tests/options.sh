#!/usr/bin/env bash
# --skip_downstream stops after featureCounts (no report); --bam_dir copies BAMs to a chosen folder.
# Needs Docker. Usage: bash tests/options.sh
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cd "$TMP"
nextflow run "$REPO/main.nf" -profile test,docker -ansi-log false -w "$TMP/work" --outdir "$TMP/out" \
    --skip_downstream --bam_dir "$TMP/my_bams" > "$TMP/run.log" 2>&1 || { tail -20 "$TMP/run.log"; exit 1; }

[[ -s "$TMP/out/featurecounts/counts.tsv" ]] || { echo "FAIL: counts.tsv missing"; exit 1; }
[[ -s "$TMP/out/multiqc/multiqc_report.html" ]] || { echo "FAIL: MultiQC report missing"; exit 1; }
[[ ! -e "$TMP/out/report" ]] || { echo "FAIL: report was built despite --skip_downstream"; exit 1; }
echo "ok: --skip_downstream gives counts + MultiQC, no report"

n_bam=$(ls "$TMP"/my_bams/*.bam 2>/dev/null | wc -l)
[[ $n_bam -eq 6 ]] || { echo "FAIL: expected 6 BAMs in --bam_dir, found $n_bam"; exit 1; }
[[ -z "$(ls "$TMP"/out/star/*.bam 2>/dev/null)" ]] || { echo "FAIL: BAMs also copied to outdir"; exit 1; }
echo "ok: --bam_dir holds the 6 BAMs"
