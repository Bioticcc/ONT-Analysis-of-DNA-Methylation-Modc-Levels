#!/usr/bin/env bash
set -Eeuo pipefail


#---------------------------------------------->
# Initializes the project and loads the shared paths, settings, logging, locking, and validation functions.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=../config/project.env
source "${PROJECT_ROOT}/config/project.env"
# shellcheck source=lib/pipeline_common.sh
source "${SCRIPT_DIR}/lib/pipeline_common.sh"

#---------------------------------------------->



#---------------------------------------------->
# Defines the help text describing the methylation analysis and its check-only and forced-rebuild modes.

usage() {
    cat <<'EOF'
Usage: 04_methylation_analysis.sh [--check-only] [--force]

Validate the finalized modBAM, create a strand-combined CpG pileup retaining
separate 5mC and 5hmC records, and calculate weighted methylation summaries.

  --check-only  Validate inputs, tools, tags, and configuration only.
  --force       Archive existing Stage 04 outputs and rebuild them.
  -h, --help    Show this help text.
EOF
}

#---------------------------------------------->



#---------------------------------------------->
# Reads command-line arguments and rejects any option other than --check-only, --force, or help.

CHECK_ONLY=0
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --check-only)
            CHECK_ONLY=1
            ;;
        --force)
            FORCE=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'Unknown argument: %s\n\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
    shift
done

#---------------------------------------------->



#---------------------------------------------->
# Sets the finalized BAM, reference, Stage 04 output, temporary, manifest, log, and Python helper paths.

ALIGNMENT_DIR="${SAMPLE_OUTPUT}/02_alignment_qc"
OUTPUT_DIR="${SAMPLE_OUTPUT}/04_methylation"
MANIFEST_DIR="${SAMPLE_OUTPUT}/00_manifests"
TMP_DIR="${SAMPLE_OUTPUT}/tmp"
LOG_DIR="${SAMPLE_OUTPUT}/logs"
ALIGNMENT_BAM="${ALIGNMENT_DIR}/${SAMPLE_ID}.aligned.sorted.bam"
ALIGNMENT_BAI="${ALIGNMENT_BAM}.bai"
REFERENCE_FAI="${REFERENCE_FASTA}.fai"
SUMMARY_HELPER="${SCRIPT_DIR}/lib/summarize_cpg_bedmethyl.py"

#---------------------------------------------->



#---------------------------------------------->
# Creates required directories, starts logging and monitoring, locks the stage, and verifies every required program.

mkdir -p "${OUTPUT_DIR}" "${MANIFEST_DIR}" "${TMP_DIR}" "${LOG_DIR}"
start_pipeline_log "${LOG_DIR}" "04_methylation_analysis"
acquire_pipeline_lock "${TMP_DIR}/04_methylation_analysis.lock"
start_pipeline_heartbeat 300 "${DATA_ROOT}"

require_executable "${SAMTOOLS_BIN}"
require_executable "${MODKIT_BIN}"
require_executable "${BGZIP_BIN}"
require_executable "${TABIX_BIN}"
require_executable python3

#---------------------------------------------->



#---------------------------------------------->
# Confirms that the Python helper used to validate and summarize CpG bedMethyl output is readable.

if [[ ! -r "${SUMMARY_HELPER}" ]]; then
    log "ERROR: missing Stage 04 summary helper: ${SUMMARY_HELPER}"
    exit 1
fi

#---------------------------------------------->



#---------------------------------------------->
# Validates Modkit thread counts, sampling size, high-depth cap, methylation window size, and probability-filter percentile.

if ! [[ "${MODKIT_THREADS}" =~ ^[1-9][0-9]*$ \
    && "${MODKIT_IO_THREADS}" =~ ^[1-9][0-9]*$ \
    && "${MODKIT_BGZF_THREADS}" =~ ^[1-9][0-9]*$ \
    && "${MODKIT_SAMPLING_THREADS}" =~ ^[1-9][0-9]*$ \
    && "${MODKIT_TAG_CHECK_READS}" =~ ^[1-9][0-9]*$ \
    && "${MODKIT_SAMPLE_READS}" =~ ^[1-9][0-9]*$ \
    && "${MODKIT_MAX_DEPTH}" =~ ^[1-9][0-9]*$ \
    && "${METHYLATION_WINDOW_SIZE}" =~ ^[1-9][0-9]*$ \
    && "${MODKIT_FILTER_PERCENTILE}" =~ ^(0([.][0-9]+)?|1([.]0+)?)$ ]]; then
    log "ERROR: invalid numeric Stage 04 configuration"
    exit 1
fi
if ! awk -v value="${MODKIT_FILTER_PERCENTILE}" 'BEGIN {exit !(value >= 0 && value < 1)}'; then
    log "ERROR: MODKIT_FILTER_PERCENTILE must be at least 0 and less than 1"
    exit 1
fi
if (( MODKIT_MAX_DEPTH > 60000 )); then
    log "ERROR: MODKIT_MAX_DEPTH must not exceed 60000; a lower cap prevents Modkit's 65,535-count saturation"
    exit 1
fi

#---------------------------------------------->



#---------------------------------------------->
# Records the sample, paths, tool versions, methylation method, filtering, and thread settings in the run log.

log "Sample: ${SAMPLE_ID}"
log "Input BAM: ${ALIGNMENT_BAM}"
log "Reference FASTA: ${REFERENCE_FASTA}"
log "Output directory: ${OUTPUT_DIR}"
log "Modkit: $("${MODKIT_BIN}" --version 2>&1 | head -n 1)"
log "Samtools: $("${SAMTOOLS_BIN}" --version | head -n 1)"
log "CpG analysis: separate 5mC/5hmC; combine_strands=yes; window=${METHYLATION_WINDOW_SIZE}"
log "Filtering: lowest ${MODKIT_FILTER_PERCENTILE} probability quantile; sampling_reads=${MODKIT_SAMPLE_READS}"
log "High-depth protection: enabled; maximum pileup observations per position=${MODKIT_MAX_DEPTH}"
log "Tag preflight: approximately ${MODKIT_TAG_CHECK_READS} primary mapped reads sampled across the indexed reference"
log "Threads: compute=${MODKIT_THREADS}; io=${MODKIT_IO_THREADS}; sampling=${MODKIT_SAMPLING_THREADS}; bgzf=${MODKIT_BGZF_THREADS}"

#---------------------------------------------->



#---------------------------------------------->
# Stops if the Stage 02 BAM, its index, the reference FASTA, or the FASTA index is missing or empty.

for required_file in "${ALIGNMENT_BAM}" "${ALIGNMENT_BAI}" "${REFERENCE_FASTA}" "${REFERENCE_FAI}"; do
    if [[ ! -s "${required_file}" ]]; then
        log "ERROR: required input is missing or empty: ${required_file}"
        log "Stage 02 must finish successfully before Stage 04 can run"
        exit 1
    fi
done

#---------------------------------------------->



#---------------------------------------------->
# Validates BAM structure and indexing, confirms coordinate sorting, and checks a mapped read for joint 5mC/5hmC tags.

log "Validating finalized BAM, BAI, reference, and modified-base tags"
"${SAMTOOLS_BIN}" quickcheck -v "${ALIGNMENT_BAM}"
"${SAMTOOLS_BIN}" idxstats -@ "${MODKIT_IO_THREADS}" "${ALIGNMENT_BAM}" >/dev/null
sort_order="$("${SAMTOOLS_BIN}" view -H "${ALIGNMENT_BAM}" \
    | awk -F '\t' '$1 == "@HD" {for (i=1; i<=NF; i++) if ($i == "SO:coordinate") found=1} END {print found+0}')"
if [[ "${sort_order}" -ne 1 ]]; then
    log "ERROR: BAM header does not declare coordinate sorting"
    exit 1
fi

set +o pipefail
first_mapped_record="$("${SAMTOOLS_BIN}" view -F 4 "${ALIGNMENT_BAM}" 2>/dev/null | head -n 1)"
set -o pipefail
if [[ -z "${first_mapped_record}" ]]; then
    log "ERROR: finalized BAM contains no mapped records"
    exit 1
fi
if [[ "${first_mapped_record}" != *$'\tMM:Z:'* \
    || "${first_mapped_record}" != *$'\tML:B:C'* \
    || "${first_mapped_record}" != *$'\tMN:i:'* \
    || "${first_mapped_record}" != *'C+h'* \
    || "${first_mapped_record}" != *'C+m'* ]]; then
    log "ERROR: sampled mapped record lacks joint C+h/C+m MM, ML, or MN tags"
    exit 1
fi
log "Confirmed coordinate sort, readable index, and joint 5mC/5hmC tags"

#---------------------------------------------->



#---------------------------------------------->
# Ends successfully after prerequisite checks when --check-only was requested.

if [[ ${CHECK_ONLY} -eq 1 ]]; then
    log "CHECK-ONLY complete: Stage 04 prerequisites validate successfully"
    exit 0
fi

#---------------------------------------------->



#---------------------------------------------->
# Lists every required final output and defines checks for completion, freshness, indexes, and both modification types.

COMPLETION_MARKER="${OUTPUT_DIR}/.stage04.complete"
REQUIRED_OUTPUTS=(
    "${OUTPUT_DIR}/modkit_tag_check.txt"
    "${OUTPUT_DIR}/modkit_valid_mm_headers.tsv"
    "${OUTPUT_DIR}/modkit_modified_bases.tsv"
    "${OUTPUT_DIR}/modkit_all_context_summary.txt"
    "${OUTPUT_DIR}/modkit_cpg_summary.txt"
    "${OUTPUT_DIR}/cpg_5mc_5hmc.bed.gz"
    "${OUTPUT_DIR}/cpg_5mc_5hmc.bed.gz.tbi"
    "${OUTPUT_DIR}/global_methylation_levels.tsv"
    "${OUTPUT_DIR}/chromosome_methylation_levels.tsv"
    "${OUTPUT_DIR}/window_methylation_levels.tsv.gz"
    "${OUTPUT_DIR}/window_methylation_levels.tsv.gz.tbi"
    "${OUTPUT_DIR}/methylation_call_coverage.tsv"
    "${OUTPUT_DIR}/methylation_analysis_report.txt"
)

outputs_validate() {
    local output
    [[ -s "${COMPLETION_MARKER}" ]] || return 1
    [[ "${COMPLETION_MARKER}" -nt "${ALIGNMENT_BAM}" ]] || return 1
    [[ "${COMPLETION_MARKER}" -nt "${REFERENCE_FAI}" ]] || return 1
    for output in "${REQUIRED_OUTPUTS[@]}"; do
        [[ -s "${output}" ]] || return 1
    done
    "${TABIX_BIN}" -l "${OUTPUT_DIR}/cpg_5mc_5hmc.bed.gz" >/dev/null 2>&1 || return 1
    "${TABIX_BIN}" -l "${OUTPUT_DIR}/window_methylation_levels.tsv.gz" >/dev/null 2>&1 || return 1
    awk -F '\t' '
        NR == 1 {next}
        $2 == "5mC" && $5 > 0 {m=1}
        $2 == "5hmC" && $5 > 0 {h=1}
        END {exit !(m && h)}
    ' "${OUTPUT_DIR}/global_methylation_levels.tsv" || return 1
}

#---------------------------------------------->



#---------------------------------------------->
# Skips a current valid result unless --force requests a recoverable rebuild.

if outputs_validate && [[ ${FORCE} -ne 1 ]]; then
    log "Existing Stage 04 result is complete, current, and indexed; skipping rebuild"
    exit 0
fi
if outputs_validate && [[ ${FORCE} -eq 1 ]]; then
    log "A valid Stage 04 result exists, but --force requested a recoverable rebuild"
fi

#---------------------------------------------->



#---------------------------------------------->
# Refuses to overwrite incomplete outputs, or archives the entire old directory first when --force was supplied.

if [[ -n "$(find "${OUTPUT_DIR}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    if [[ ${FORCE} -ne 1 ]]; then
        log "ERROR: Stage 04 contains incomplete, stale, or unvalidated outputs"
        log "Inspect ${OUTPUT_DIR}; rerun with --force to archive it and rebuild"
        exit 1
    fi
    STALE_ROOT="${TMP_DIR}/stale_stage04_outputs"
    STALE_DIR="${STALE_ROOT}/${PIPELINE_RUN_STAMP}"
    mkdir -p "${STALE_ROOT}"
    if [[ -e "${STALE_DIR}" ]]; then
        log "ERROR: stale-output archive target already exists: ${STALE_DIR}"
        exit 1
    fi
    mv "${OUTPUT_DIR}" "${STALE_DIR}"
    mkdir -p "${OUTPUT_DIR}"
    log "Archived prior Stage 04 directory without deleting it: ${STALE_DIR}"
fi

#---------------------------------------------->



#---------------------------------------------->
# Builds a fingerprint from the BAM and analysis files so compatible interrupted work can be safely resumed.

work_fingerprint="$({
    printf '%s\n' "$(stat -c '%n:%s:%Y' "${ALIGNMENT_BAM}")"
    printf '%s\n' "$(sha256sum "${REFERENCE_FAI}" "${PROJECT_ROOT}/config/project.env" "$0" "${SUMMARY_HELPER}")"
} | sha256sum | awk '{print $1}')"
if [[ ${FORCE} -eq 1 ]]; then
    work_fingerprint="${work_fingerprint}.force.${PIPELINE_RUN_STAMP}"
fi
WORK_DIR="${TMP_DIR}/stage04_work/${work_fingerprint}"
PARTIAL_TOKEN="${PIPELINE_RUN_STAMP}"
mkdir -p "${WORK_DIR}"
log "Restartable working directory: ${WORK_DIR}"

#---------------------------------------------->



#---------------------------------------------->
# Names the Modkit tag-validation tables and defines a check requiring exactly cytosine 5mC and 5hmC codes.

TAG_CHECK_TEXT="${WORK_DIR}/modkit_tag_check.txt"
TAG_HEADERS="${WORK_DIR}/modkit_valid_mm_headers.tsv"
TAG_BASES="${WORK_DIR}/modkit_modified_bases.tsv"
TAG_DEBUG="${WORK_DIR}/modkit_tag_check.debug.log"

tag_tables_validate() {
    [[ -s "${TAG_CHECK_TEXT}" && -s "${TAG_HEADERS}" && -s "${TAG_BASES}" ]] || return 1
    awk -F '\t' '
        NR == 1 {next}
        $2 == "C" && $3 == "m" {m=1; next}
        $2 == "C" && $3 == "h" {h=1; next}
        {bad=1}
        END {exit !(m && h && !bad)}
    ' "${TAG_BASES}"
}

#---------------------------------------------->



#---------------------------------------------->
# Reuses a valid prior tag check or samples mapped records across the indexed reference and atomically saves its validation tables.

if tag_tables_validate; then
    log "Reusing validated sampled Modkit tag check from restartable work"
else
    tag_partial_dir="${WORK_DIR}/tag_check.partial.${PARTIAL_TOKEN}"
    tag_text_partial="${TAG_CHECK_TEXT}.partial.${PARTIAL_TOKEN}"
    tag_debug_partial="${TAG_DEBUG}.partial.${PARTIAL_TOKEN}"
    mkdir -p "${tag_partial_dir}"
    log "Running sampled Modkit MM/ML tag validation on mapped primary records"
    "${MODKIT_BIN}" modbam check-tags \
        --threads "${MODKIT_THREADS}" \
        --num-reads "${MODKIT_TAG_CHECK_READS}" \
        --mapped-only \
        --suppress-progress \
        --log-filepath "${tag_debug_partial}" \
        --out-dir "${tag_partial_dir}" \
        --prefix tag_check \
        "${ALIGNMENT_BAM}" \
        > "${tag_text_partial}" 2>&1
    mv "${tag_text_partial}" "${TAG_CHECK_TEXT}"
    mv "${tag_debug_partial}" "${TAG_DEBUG}"
    mv "${tag_partial_dir}/tag_check_valid_mm_headers.tsv" "${TAG_HEADERS}"
    mv "${tag_partial_dir}/tag_check_modified_bases.tsv" "${TAG_BASES}"
    rmdir "${tag_partial_dir}"
    if ! tag_tables_validate; then
        log "ERROR: tag check did not confirm exactly C:m and C:h modification codes"
        exit 1
    fi
fi

#---------------------------------------------->



#---------------------------------------------->
# Reuses or creates an all-context Modkit summary using reference-matched calls and the configured probability filter.

ALL_CONTEXT_SUMMARY="${WORK_DIR}/modkit_all_context_summary.txt"
ALL_CONTEXT_DEBUG="${WORK_DIR}/modkit_all_context_summary.debug.log"
if [[ -s "${ALL_CONTEXT_SUMMARY}" ]] \
    && grep -q $'^total_reads_used\t' "${ALL_CONTEXT_SUMMARY}"; then
    log "Reusing completed all-context Modkit summary"
else
    log "Generating matched, reference-aligned all-context Modkit summary"
    "${MODKIT_BIN}" summary \
        --tsv \
        --threads "${MODKIT_THREADS}" \
        --io-threads "${MODKIT_IO_THREADS}" \
        --reference "${REFERENCE_FASTA}" \
        --matched-only \
        --filter-quantile "${MODKIT_FILTER_PERCENTILE}" \
        --suppress-progress \
        --log-filepath "${ALL_CONTEXT_DEBUG}.partial.${PARTIAL_TOKEN}" \
        "${ALIGNMENT_BAM}" \
        > "${ALL_CONTEXT_SUMMARY}.partial.${PARTIAL_TOKEN}"
    mv "${ALL_CONTEXT_SUMMARY}.partial.${PARTIAL_TOKEN}" "${ALL_CONTEXT_SUMMARY}"
    mv "${ALL_CONTEXT_DEBUG}.partial.${PARTIAL_TOKEN}" "${ALL_CONTEXT_DEBUG}"
fi

#---------------------------------------------->



#---------------------------------------------->
# Converts the reference FASTA index into a whole-contig BED file used for chromosome-level Modkit statistics.

GENOME_REGIONS="${WORK_DIR}/reference_contigs.bed"
awk -F '\t' 'BEGIN {OFS="\t"} {print $1, 0, $2}' "${REFERENCE_FAI}" \
    > "${GENOME_REGIONS}.partial.${PARTIAL_TOKEN}"
mv "${GENOME_REGIONS}.partial.${PARTIAL_TOKEN}" "${GENOME_REGIONS}"

#---------------------------------------------->



#---------------------------------------------->
# Reuses a valid CpG bedMethyl or runs high-depth-safe Modkit pileup with combined strands and separate 5mC and 5hmC records.

CPG_BED="${WORK_DIR}/cpg_5mc_5hmc.bed.gz"
CPG_BED_TBI="${CPG_BED}.tbi"
PILEUP_DEBUG="${WORK_DIR}/modkit_pileup.debug.log"
if [[ -s "${CPG_BED}" && -s "${CPG_BED_TBI}" ]] \
    && "${TABIX_BIN}" -l "${CPG_BED}" >/dev/null 2>&1; then
    log "Reusing validated CpG bedMethyl from restartable work"
else
    cpg_partial="${WORK_DIR}/cpg_5mc_5hmc.partial.${PARTIAL_TOKEN}.bed.gz"
    pileup_debug_partial="${PILEUP_DEBUG}.partial.${PARTIAL_TOKEN}"
    log "Running Modkit CpG pileup with separate 5mC and 5hmC records"
    log "RUN: modkit pileup --cpg --combine-strands --modified-bases 5mC 5hmC --filter-percentile ${MODKIT_FILTER_PERCENTILE} --num-reads ${MODKIT_SAMPLE_READS} --high-depth --max-depth ${MODKIT_MAX_DEPTH}"
    "${MODKIT_BIN}" pileup \
        --threads "${MODKIT_THREADS}" \
        --io-threads "${MODKIT_IO_THREADS}" \
        --sampling-threads "${MODKIT_SAMPLING_THREADS}" \
        --bgzf-threads "${MODKIT_BGZF_THREADS}" \
        --reference "${REFERENCE_FASTA}" \
        --cpg \
        --combine-strands \
        --modified-bases 5mC 5hmC \
        --filter-percentile "${MODKIT_FILTER_PERCENTILE}" \
        --num-reads "${MODKIT_SAMPLE_READS}" \
        --high-depth \
        --max-depth "${MODKIT_MAX_DEPTH}" \
        --bgzf \
        --suppress-progress \
        --log-filepath "${pileup_debug_partial}" \
        "${ALIGNMENT_BAM}" \
        "${cpg_partial}"
    if [[ ! -s "${cpg_partial}" ]]; then
        log "ERROR: Modkit CpG pileup is empty"
        exit 1
    fi
    "${TABIX_BIN}" -f -p bed "${cpg_partial}"
    "${TABIX_BIN}" -l "${cpg_partial}" >/dev/null
    mv "${cpg_partial}" "${CPG_BED}"
    mv "${cpg_partial}.tbi" "${CPG_BED_TBI}"
    mv "${pileup_debug_partial}" "${PILEUP_DEBUG}"
fi

#---------------------------------------------->



#---------------------------------------------->
# Reuses or calculates chromosome-level statistics from the indexed CpG bedMethyl file.

CPG_SUMMARY="${WORK_DIR}/modkit_cpg_summary.txt"
CPG_STATS_DEBUG="${WORK_DIR}/modkit_cpg_stats.debug.log"
if [[ -s "${CPG_SUMMARY}" ]] && grep -q $'^#chrom\t' "${CPG_SUMMARY}"; then
    log "Reusing completed CpG Modkit stats summary"
else
    log "Aggregating CpG pileup by chromosome with Modkit stats"
    "${MODKIT_BIN}" stats \
        --threads "${MODKIT_THREADS}" \
        --io-threads "${MODKIT_IO_THREADS}" \
        --regions "${GENOME_REGIONS}" \
        --out-table "${CPG_SUMMARY}.partial.${PARTIAL_TOKEN}" \
        --log-filepath "${CPG_STATS_DEBUG}.partial.${PARTIAL_TOKEN}" \
        "${CPG_BED}"
    mv "${CPG_SUMMARY}.partial.${PARTIAL_TOKEN}" "${CPG_SUMMARY}"
    mv "${CPG_STATS_DEBUG}.partial.${PARTIAL_TOKEN}" "${CPG_STATS_DEBUG}"
fi

#---------------------------------------------->



#---------------------------------------------->
# Uses the Python helper to validate paired records and calculate weighted global, chromosome, window, and call-coverage tables.

GLOBAL_LEVELS="${WORK_DIR}/global_methylation_levels.tsv"
CHROMOSOME_LEVELS="${WORK_DIR}/chromosome_methylation_levels.tsv"
WINDOW_LEVELS_TSV="${WORK_DIR}/window_methylation_levels.tsv"
CALL_COVERAGE="${WORK_DIR}/methylation_call_coverage.tsv"

log "Validating paired bedMethyl records and calculating weighted results"
python3 "${SUMMARY_HELPER}" \
    --fasta "${REFERENCE_FASTA}" \
    --fai "${REFERENCE_FAI}" \
    --bedmethyl "${CPG_BED}" \
    --modkit-stats "${CPG_SUMMARY}" \
    --window-size "${METHYLATION_WINDOW_SIZE}" \
    --global-out "${GLOBAL_LEVELS}.partial.${PARTIAL_TOKEN}" \
    --chromosome-out "${CHROMOSOME_LEVELS}.partial.${PARTIAL_TOKEN}" \
    --window-out "${WINDOW_LEVELS_TSV}.partial.${PARTIAL_TOKEN}" \
    --call-coverage-out "${CALL_COVERAGE}.partial.${PARTIAL_TOKEN}"
mv "${GLOBAL_LEVELS}.partial.${PARTIAL_TOKEN}" "${GLOBAL_LEVELS}"
mv "${CHROMOSOME_LEVELS}.partial.${PARTIAL_TOKEN}" "${CHROMOSOME_LEVELS}"
mv "${WINDOW_LEVELS_TSV}.partial.${PARTIAL_TOKEN}" "${WINDOW_LEVELS_TSV}"
mv "${CALL_COVERAGE}.partial.${PARTIAL_TOKEN}" "${CALL_COVERAGE}"

#---------------------------------------------->



#---------------------------------------------->
# BGZF-compresses and Tabix-indexes the window-level methylation table for fast genomic-region queries.

"${BGZIP_BIN}" -@ "${MODKIT_BGZF_THREADS}" -c "${WINDOW_LEVELS_TSV}" \
    > "${WORK_DIR}/window_methylation_levels.tsv.gz.partial.${PARTIAL_TOKEN}"
mv "${WORK_DIR}/window_methylation_levels.tsv.gz.partial.${PARTIAL_TOKEN}" \
    "${WORK_DIR}/window_methylation_levels.tsv.gz"
"${TABIX_BIN}" -f -p bed "${WORK_DIR}/window_methylation_levels.tsv.gz"
"${TABIX_BIN}" -l "${WORK_DIR}/window_methylation_levels.tsv.gz" >/dev/null

#---------------------------------------------->



#---------------------------------------------->
# Defines helpers that retrieve named global methylation and call-coverage values for the readable report.

global_value() {
    local category="$1"
    local column="$2"
    awk -F '\t' -v category="${category}" -v column="${column}" \
        '$2 == category {print $column; exit}' "${GLOBAL_LEVELS}"
}
coverage_value() {
    local column_name="$1"
    awk -F '\t' -v name="${column_name}" '
        NR == 1 {for (i=1; i<=NF; i++) if ($i == name) column=i; next}
        $1 == "genome" {print $column; exit}
    ' "${CALL_COVERAGE}"
}

#---------------------------------------------->



#---------------------------------------------->
# Builds a human-readable report describing the method, weighted 5mC/5hmC levels, call coverage, and interpretation.

REPORT_FILE="${WORK_DIR}/methylation_analysis_report.txt"
{
    printf 'ONT CpG 5mC/5hmC methylation analysis report\n'
    printf 'Sample: %s\n' "${SAMPLE_ID}"
    printf 'Reference: %s\n' "${REFERENCE_NAME}"
    printf 'Input BAM: %s\n' "${ALIGNMENT_BAM}"
    printf 'Generated: %s\n\n' "$(date --iso-8601=seconds)"
    printf 'Method\n'
    printf '  Tool: %s\n' "$("${MODKIT_BIN}" --version 2>&1 | head -n 1)"
    printf '  Context: CpG; complementary strands combined into one CpG dyad\n'
    printf '  Modifications remain separate: m=5mC and h=5hmC\n'
    printf '  Primary alignments only; positions above %s observations are downsampled by Modkit high-depth mode\n' "${MODKIT_MAX_DEPTH}"
    printf '  Calls in the lowest %s probability quantile are classified as failed/filtered\n' "${MODKIT_FILTER_PERCENTILE}"
    printf '  Global percentages are coverage-weighted: sum(modified calls) / sum(valid calls)\n\n'
    printf 'Global CpG levels\n'
    printf '  Valid cytosine calls: %s\n' "$(global_value 5mC 5)"
    printf '  5mC: %s calls; %s%% of valid calls; %s%% of modified calls\n' \
        "$(global_value 5mC 4)" "$(global_value 5mC 6)" "$(global_value 5mC 8)"
    printf '  5hmC: %s calls; %s%% of valid calls; %s%% of modified calls\n' \
        "$(global_value 5hmC 4)" "$(global_value 5hmC 6)" "$(global_value 5hmC 8)"
    printf '  Total modified: %s calls; %s%% of valid calls\n' \
        "$(global_value total_modified 4)" "$(global_value total_modified 6)"
    printf '  Canonical C: %s calls; %s%% of valid calls\n\n' \
        "$(global_value canonical_C 4)" "$(global_value canonical_C 6)"
    printf 'Call and site coverage\n'
    printf '  Reference CpG sites: %s\n' "$(coverage_value reference_cpg_sites)"
    printf '  CpG sites with >=1 valid call: %s (%s%%)\n' \
        "$(coverage_value cpg_sites_with_valid_calls)" \
        "$(coverage_value pct_reference_cpg_sites_with_valid_calls)"
    printf '  Failed/filtered calls: %s\n' "$(coverage_value failed_filtered_calls)"
    printf '  No-call observations: %s\n' "$(coverage_value no_call_calls)"
    printf '  Deletion observations: %s\n' "$(coverage_value deletion_calls)"
    printf '  Different-base observations: %s\n' "$(coverage_value different_base_calls)"
    printf '  Valid calls among all observations: %s%%\n\n' "$(coverage_value pct_valid_calls_of_observations)"
    printf 'Interpretation notes\n'
    printf '  The 5mC and 5hmC percentages share the same valid-call denominator and can\n'
    printf '  therefore be compared directly. The all-context summary is informational;\n'
    printf '  final reported quantities come from the validated CpG bedMethyl counts.\n'
    printf '  Missing CpG-site calls are distinct from low sequencing depth; Stage 03\n'
    printf '  coverage should be consulted when investigating a regional gap.\n'
    printf '  The high-depth cap prevents Modkit count saturation at localized extreme\n'
    printf '  coverage; reported counts at those positions represent the capped pileup.\n'
} > "${REPORT_FILE}"

#---------------------------------------------->



#---------------------------------------------->
# Lists and verifies every file that must be complete and nonempty before Stage 04 results are promoted.

FINAL_WORK_FILES=(
    modkit_tag_check.txt
    modkit_valid_mm_headers.tsv
    modkit_modified_bases.tsv
    modkit_all_context_summary.txt
    modkit_cpg_summary.txt
    cpg_5mc_5hmc.bed.gz
    cpg_5mc_5hmc.bed.gz.tbi
    global_methylation_levels.tsv
    chromosome_methylation_levels.tsv
    window_methylation_levels.tsv.gz
    window_methylation_levels.tsv.gz.tbi
    methylation_call_coverage.tsv
    methylation_analysis_report.txt
)
for final_work_file in "${FINAL_WORK_FILES[@]}"; do
    if [[ ! -s "${WORK_DIR}/${final_work_file}" ]]; then
        log "ERROR: validated Stage 04 output is missing or empty: ${final_work_file}"
        exit 1
    fi
done

#---------------------------------------------->



#---------------------------------------------->
# Moves the validated analysis products from restartable work into the official Stage 04 output directory.

log "Promoting validated Stage 04 outputs atomically within the 20 TB filesystem"
for final_work_file in "${FINAL_WORK_FILES[@]}"; do
    mv "${WORK_DIR}/${final_work_file}" "${OUTPUT_DIR}/${final_work_file}"
done

#---------------------------------------------->



#---------------------------------------------->
# Preserves useful Modkit debug logs in the main log directory with timestamped names.

DEBUG_LOGS=(
    modkit_tag_check.debug.log
    modkit_all_context_summary.debug.log
    modkit_pileup.debug.log
    modkit_cpg_stats.debug.log
)
for debug_log in "${DEBUG_LOGS[@]}"; do
    if [[ -s "${WORK_DIR}/${debug_log}" ]]; then
        debug_destination="${LOG_DIR}/04_${debug_log%.log}.${PIPELINE_RUN_STAMP}.log"
        mv "${WORK_DIR}/${debug_log}" "${debug_destination}"
        log "Retained Modkit debug log: ${debug_destination}"
    fi
done

#---------------------------------------------->



#---------------------------------------------->
# Writes a timestamped reproducibility manifest with settings, tool versions, checksums, and output file sizes.

MANIFEST_FILE="${MANIFEST_DIR}/${SAMPLE_ID}.stage04.${PIPELINE_RUN_STAMP}.manifest.tsv"
{
    printf 'record_type\tname\tvalue\n'
    printf 'metadata\tsample_id\t%s\n' "${SAMPLE_ID}"
    printf 'metadata\treference\t%s\n' "${REFERENCE_NAME}"
    printf 'metadata\tinput_bam\t%s\n' "${ALIGNMENT_BAM}"
    printf 'metadata\tinput_bam_bytes\t%s\n' "$(stat -c '%s' "${ALIGNMENT_BAM}")"
    printf 'metadata\tcontext\tCpG\n'
    printf 'metadata\tcombine_strands\ttrue\n'
    printf 'metadata\tmodified_bases\t5mC,5hmC\n'
    printf 'metadata\tfilter_percentile\t%s\n' "${MODKIT_FILTER_PERCENTILE}"
    printf 'metadata\ttag_check_sample_reads\t%s\n' "${MODKIT_TAG_CHECK_READS}"
    printf 'metadata\tthreshold_sample_reads\t%s\n' "${MODKIT_SAMPLE_READS}"
    printf 'metadata\thigh_depth_mode\ttrue\n'
    printf 'metadata\tmax_pileup_depth\t%s\n' "${MODKIT_MAX_DEPTH}"
    printf 'metadata\twindow_size\t%s\n' "${METHYLATION_WINDOW_SIZE}"
    printf 'tool\tmodkit\t%s\n' "$("${MODKIT_BIN}" --version 2>&1 | head -n 1)"
    printf 'tool\tsamtools\t%s\n' "$("${SAMTOOLS_BIN}" --version | head -n 1)"
    printf 'checksum\tstage04_script_sha256\t%s\n' "$(sha256sum "$0" | awk '{print $1}')"
    printf 'checksum\tsummary_helper_sha256\t%s\n' "$(sha256sum "${SUMMARY_HELPER}" | awk '{print $1}')"
    printf 'checksum\tproject_config_sha256\t%s\n' "$(sha256sum "${PROJECT_ROOT}/config/project.env" | awk '{print $1}')"
    for output in "${REQUIRED_OUTPUTS[@]}"; do
        printf 'output_bytes\t%s\t%s\n' "$(basename -- "${output}")" "$(stat -c '%s' "${output}")"
    done
} > "${MANIFEST_FILE}.partial.${PARTIAL_TOKEN}"
mv "${MANIFEST_FILE}.partial.${PARTIAL_TOKEN}" "${MANIFEST_FILE}"
ln -sfn "$(basename -- "${MANIFEST_FILE}")" "${MANIFEST_DIR}/${SAMPLE_ID}.stage04.latest.manifest.tsv"

#---------------------------------------------->



#---------------------------------------------->
# Writes the completion marker last, then reruns the complete final-output validation.

{
    printf 'status\tcomplete\n'
    printf 'sample_id\t%s\n' "${SAMPLE_ID}"
    printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
    printf 'input_bam\t%s\n' "${ALIGNMENT_BAM}"
    printf 'manifest\t%s\n' "${MANIFEST_FILE}"
} > "${COMPLETION_MARKER}.partial.${PARTIAL_TOKEN}"
mv "${COMPLETION_MARKER}.partial.${PARTIAL_TOKEN}" "${COMPLETION_MARKER}"

if ! outputs_validate; then
    log "ERROR: final Stage 04 validation failed after promotion"
    exit 1
fi

#---------------------------------------------->



#---------------------------------------------->
# Safely removes only the expected restartable working directory after all final validation succeeds.

case "${WORK_DIR}" in
    "${TMP_DIR}"/stage04_work/*)
        rm -rf -- "${WORK_DIR}"
        ;;
    *)
        log "ERROR: refusing to clean unexpected working directory: ${WORK_DIR}"
        exit 1
        ;;
esac

#---------------------------------------------->



#---------------------------------------------->
# Reports successful completion and the locations of the principal Stage 04 results.

log "Stage 04 complete"
log "Methylation report: ${OUTPUT_DIR}/methylation_analysis_report.txt"
log "Global levels: ${OUTPUT_DIR}/global_methylation_levels.tsv"
log "CpG bedMethyl: ${OUTPUT_DIR}/cpg_5mc_5hmc.bed.gz"
log "Manifest: ${MANIFEST_FILE}"
#---------------------------------------------->
