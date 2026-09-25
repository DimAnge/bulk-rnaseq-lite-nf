#!/usr/bin/env bash
# Input-validation tests. They need no Docker: every check in main.nf runs before any
# process starts, and passing cases use -preview (build the workflow, run nothing).
# Usage: bash tests/validation.sh
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
READS="$REPO/assets/test/reads"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cd "$TMP"
pass=0; fail=0

run() { nextflow run "$REPO/main.nf" -ansi-log false -w "$TMP/work" "$@" 2>&1; }

expect_error() {  # expect_error "<message fragment>" <nextflow args...>
    local msg="$1"; shift
    local out code
    out="$(run "$@")"; code=$?
    if [[ $code -ne 0 && "$out" == *"$msg"* ]]; then echo "PASS  error: $msg"; pass=$((pass + 1))
    else echo "FAIL  expected error '$msg' (exit $code)"; echo "$out" | tail -5; fail=$((fail + 1)); fi
}

expect_ok() {  # expect_ok "<label>" <nextflow args...>
    local label="$1"; shift
    local out code
    out="$(run "$@" -preview)"; code=$?
    if [[ $code -eq 0 ]]; then echo "PASS  ok: $label"; pass=$((pass + 1))
    else echo "FAIL  expected success: $label"; echo "$out" | tail -5; fail=$((fail + 1)); fi
}

sheet() {  # sheet <name> <csv body...>  -> writes $TMP/<name>.csv
    local name="$1"; shift
    printf '%s\n' "$@" > "$TMP/$name.csv"
}

HEADER="sample,fastq_1,fastq_2,condition,batch"
OK1="ctrl_1,$READS/ctrl_1/R1.fastq.gz,$READS/ctrl_1/R2.fastq.gz,control,A"
OK2="trt_1,$READS/trt_1/R1.fastq.gz,$READS/trt_1/R2.fastq.gz,treated,A"
OK3="ctrl_2,$READS/ctrl_2/R1.fastq.gz,,control,B"
OK4="trt_2,$READS/trt_2/R1.fastq.gz,,treated,B"

sheet ok "$HEADER" "$OK1" "$OK2" "$OK3" "$OK4"
sheet no_condition "sample,fastq_1,fastq_2" "ctrl_1,$READS/ctrl_1/R1.fastq.gz,"
sheet duplicate "$HEADER" "$OK1" "$OK1" "$OK2"
sheet missing_fastq "$HEADER" "$OK1" "trt_1,$READS/nope/R1.fastq.gz,,treated,A"
sheet plain_fastq "$HEADER" "$OK1" "trt_1,$READS/trt_1/R1.fastq,,treated,A"
sheet one_condition "$HEADER" "$OK1" "$OK3"
sheet bad_level "$HEADER" "$OK1" "trt_1,$READS/trt_1/R1.fastq.gz,,KO-1,A"
sheet empty_batch "$HEADER" "$OK1" "$OK2" "ctrl_2,$READS/ctrl_2/R1.fastq.gz,,control,"
sheet semicolon "sample;fastq_1;fastq_2;condition;batch" "ctrl_1;$READS/ctrl_1/R1.fastq.gz;;control;A"
sheet blank_rows "$HEADER" "$OK1" "$OK2" ",,,," ",,,,"
# Excel-style: UTF-8 BOM + CRLF line endings
printf '\xef\xbb\xbf%s\r\n%s\r\n%s\r\n' "$HEADER" "$OK1" "$OK2" > "$TMP/bom_crlf.csv"
# Relative paths resolve against the samplesheet's own folder
mkdir -p "$TMP/rel" && ln -s "$READS" "$TMP/rel/reads"
sheet rel/relative "$HEADER" "ctrl_1,reads/ctrl_1/R1.fastq.gz,,control,A" "trt_1,reads/trt_1/R1.fastq.gz,,treated,A"

expect_error "Please provide a samplesheet"                 
expect_error "Samplesheet not found"                        --input "$TMP/nope.csv"
expect_error "missing required column(s): condition"        --input "$TMP/no_condition.csv"
expect_error "Duplicate sample name(s) in samplesheet: ctrl_1" --input "$TMP/duplicate.csv"
expect_error "FASTQ file not found"                         --input "$TMP/missing_fastq.csv"
expect_error "must be gzipped"                              --input "$TMP/plain_fastq.csv"
expect_error "at least two conditions"                      --input "$TMP/one_condition.csv"
expect_error "condition values must start with a letter"    --input "$TMP/bad_level.csv"
expect_error "empty design values for sample(s): ctrl_2"    --input "$TMP/empty_batch.csv" --design '~ batch + condition'
expect_error "--strandedness must be one of"                --input "$TMP/ok.csv" --strandedness revrse
expect_error "--genome must be one of"                      --input "$TMP/ok.csv" --genome hg19
expect_error "--fasta and --gtf must be given together"     --input "$TMP/ok.csv" --fasta "$REPO/assets/test/genome.fa"
expect_error "column(s) not in the samplesheet: batch2"     --input "$TMP/ok.csv" --design '~ batch2 + condition'
expect_error "--design must end with 'condition'"           --input "$TMP/ok.csv" --design '~ condition + batch'
expect_error "interaction terms"                            --input "$TMP/ok.csv" --design '~ batch * condition'
expect_error "Contrast 'treated_vs_ctrl' is not valid"      --input "$TMP/ok.csv" --contrasts treated_vs_ctrl
expect_error "compares a condition with itself"             --input "$TMP/ok.csv" --contrasts control_vs_control
expect_error "separated by semicolons"                      --input "$TMP/semicolon.csv"
expect_error "--padj_cutoff must be between 0 and 1"        --input "$TMP/ok.csv" --padj_cutoff 5

expect_ok "valid samplesheet"            --input "$TMP/ok.csv" --design '~ batch + condition' --contrasts treated_vs_control
expect_ok "Excel BOM + CRLF samplesheet" --input "$TMP/bom_crlf.csv" --design '~ batch + condition'
expect_ok "blank trailing rows (Excel)"  --input "$TMP/blank_rows.csv"
expect_ok "relative FASTQ paths"         --input "$TMP/rel/relative.csv"
expect_ok "test profile"                 -profile test
expect_ok "mouse genome"                 --input "$TMP/ok.csv" --genome GRCm39

echo "---- $pass passed, $fail failed"
[[ $fail -eq 0 ]]
