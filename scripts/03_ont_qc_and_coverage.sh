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
# Defines the help text describing Stage 03 and its optional check-only and forced-rebuild modes.

usage() {
    cat <<'EOF'
Usage: 03_ont_qc_and_coverage.sh [--check-only] [--force]

Build alignment QC, exact run-length depth, chromosome coverage, and 100 kb
window coverage from the finalized Stage 02 modBAM.

  --check-only  Validate inputs, indexes, tools, tags, and configuration only.
  --force       Archive existing Stage 03 outputs and rebuild them.
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
# Sets the finalized BAM input, reference index, output, temporary, manifest, log, and Python helper paths.

ALIGNMENT_DIR="${SAMPLE_OUTPUT}/02_alignment_qc"
OUTPUT_DIR="${SAMPLE_OUTPUT}/03_ont_qc_coverage"
MANIFEST_DIR="${SAMPLE_OUTPUT}/00_manifests"
TMP_DIR="${SAMPLE_OUTPUT}/tmp"
LOG_DIR="${SAMPLE_OUTPUT}/logs"
ALIGNMENT_BAM="${ALIGNMENT_DIR}/${SAMPLE_ID}.aligned.sorted.bam"
ALIGNMENT_BAI="${ALIGNMENT_BAM}.bai"
REFERENCE_FAI="${REFERENCE_FASTA}.fai"
SUMMARY_HELPER="${SCRIPT_DIR}/lib/summarize_mosdepth.py"

#---------------------------------------------->



#---------------------------------------------->
# Creates required directories, starts logging and monitoring, locks the stage, and verifies every required program.

mkdir -p "${OUTPUT_DIR}" "${MANIFEST_DIR}" "${TMP_DIR}" "${LOG_DIR}"
start_pipeline_log "${LOG_DIR}" "03_ont_qc_and_coverage"
acquire_pipeline_lock "${TMP_DIR}/03_ont_qc_and_coverage.lock"
start_pipeline_heartbeat 300 "${DATA_ROOT}"

require_executable "${SAMTOOLS_BIN}"
require_executable "${MOSDEPTH_BIN}"
require_executable "${BGZIP_BIN}"
require_executable "${TABIX_BIN}"
require_executable python3

#---------------------------------------------->



#---------------------------------------------->
# Confirms that the Python helper used to validate and summarize Mosdepth output is readable.

if [[ ! -r "${SUMMARY_HELPER}" ]]; then
    log "ERROR: missing Stage 03 summary helper: ${SUMMARY_HELPER}"
    exit 1
fi

#---------------------------------------------->



#---------------------------------------------->
# Parses and validates all coverage settings, including thresholds, thread counts, flags, MAPQ, and low-coverage cutoffs.

IFS=',' read -r -a COVERAGE_THRESHOLD_ARRAY <<< "${COVERAGE_THRESHOLDS}"
if [[ ${#COVERAGE_THRESHOLD_ARRAY[@]} -lt 4 ]] \
    || [[ ",${COVERAGE_THRESHOLDS}," != *,1,* ]] \
    || [[ ",${COVERAGE_THRESHOLDS}," != *,5,* ]] \
    || [[ ",${COVERAGE_THRESHOLDS}," != *,10,* ]] \
    || [[ ",${COVERAGE_THRESHOLDS}," != *,20,* ]]; then
    log "ERROR: COVERAGE_THRESHOLDS must include 1,5,10,20"
    exit 1
fi
if ! [[ "${QC_THREADS}" =~ ^[1-9][0-9]*$ \
    && "${COVERAGE_WINDOW_SIZE}" =~ ^[1-9][0-9]*$ \
    && "${COVERAGE_EXCLUDE_FLAGS}" =~ ^[0-9]+$ \
    && "${COVERAGE_MIN_MAPQ}" =~ ^[0-9]+$ \
    && "${LOW_COVERAGE_DEPTH_THRESHOLD}" =~ ^[2-9][0-9]*$ \
    && "${LOW_COVERAGE_MIN_1X_BREADTH_PERCENT}" =~ ^([0-9]+([.][0-9]+)?|[.][0-9]+)$ ]]; then
    log "ERROR: invalid numeric Stage 03 configuration"
    exit 1
fi
if ! awk -v value="${LOW_COVERAGE_MIN_1X_BREADTH_PERCENT}" \
    'BEGIN {exit !(value >= 0 && value <= 100)}'; then
    log "ERROR: LOW_COVERAGE_MIN_1X_BREADTH_PERCENT must be between 0 and 100"
    exit 1
fi
previous_threshold=0
for threshold in "${COVERAGE_THRESHOLD_ARRAY[@]}"; do
    if ! [[ "${threshold}" =~ ^[1-9][0-9]*$ ]] || (( threshold <= previous_threshold )); then
        log "ERROR: COVERAGE_THRESHOLDS must be unique positive integers in ascending order"
        exit 1
    fi
    previous_threshold="${threshold}"
done

#---------------------------------------------->



#---------------------------------------------->
# Records the sample, inputs, tool versions, and exact coverage definitions in the run log.

log "Sample: ${SAMPLE_ID}"
log "Input BAM: ${ALIGNMENT_BAM}"
log "Reference FAI: ${REFERENCE_FAI}"
log "Output directory: ${OUTPUT_DIR}"
log "Samtools: $("${SAMTOOLS_BIN}" --version | head -n 1)"
log "Mosdepth: $("${MOSDEPTH_BIN}" --version 2>&1 | head -n 1)"
log "BGZF/tabix: $("${TABIX_BIN}" --version 2>&1 | head -n 1)"
log "Coverage: window=${COVERAGE_WINDOW_SIZE}; thresholds=${COVERAGE_THRESHOLDS}; min_mapq=${COVERAGE_MIN_MAPQ}; exclude_flags=${COVERAGE_EXCLUDE_FLAGS}"
log "Low coverage definitions: exact depth <${LOW_COVERAGE_DEPTH_THRESHOLD}x; low window <${LOW_COVERAGE_MIN_1X_BREADTH_PERCENT}% breadth at >=1x"

#---------------------------------------------->



#---------------------------------------------->
# Stops if the finalized Stage 02 BAM, BAM index, or reference FASTA index is missing or empty.

for required_file in "${ALIGNMENT_BAM}" "${ALIGNMENT_BAI}" "${REFERENCE_FAI}"; do
    if [[ ! -s "${required_file}" ]]; then
        log "ERROR: required input is missing or empty: ${required_file}"
        log "Stage 02 must finish successfully before Stage 03 can run"
        exit 1
    fi
done

#---------------------------------------------->



#---------------------------------------------->
# Validates BAM structure and indexing, confirms coordinate sorting, and checks a mapped read for MM, ML, and MN tags.

log "Validating finalized BAM and BAI"
"${SAMTOOLS_BIN}" quickcheck -v "${ALIGNMENT_BAM}"
"${SAMTOOLS_BIN}" idxstats -@ "${QC_THREADS}" "${ALIGNMENT_BAM}" >/dev/null
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
    || "${first_mapped_record}" != *$'\tMN:i:'* ]]; then
    log "ERROR: MM, ML, and MN modified-base tags were not found in the sampled mapped record"
    exit 1
fi
log "Confirmed coordinate sort, readable index, and MM/ML/MN tags"

#---------------------------------------------->



#---------------------------------------------->
# Ends successfully after prerequisite checks when --check-only was requested.

if [[ ${CHECK_ONLY} -eq 1 ]]; then
    log "CHECK-ONLY complete: Stage 03 prerequisites validate successfully"
    exit 0
fi

#---------------------------------------------->



#---------------------------------------------->
# Lists every required final output and defines how to verify completion, freshness, nonempty files, and genomic indexes.

COMPLETION_MARKER="${OUTPUT_DIR}/.stage03.complete"
REQUIRED_OUTPUTS=(
    "${OUTPUT_DIR}/alignment_flagstat.txt"
    "${OUTPUT_DIR}/alignment_stats.txt"
    "${OUTPUT_DIR}/alignment_idxstats.tsv"
    "${OUTPUT_DIR}/coverage_summary.tsv"
    "${OUTPUT_DIR}/chromosome_coverage.tsv"
    "${OUTPUT_DIR}/window_coverage.tsv.gz"
    "${OUTPUT_DIR}/window_coverage.tsv.gz.tbi"
    "${OUTPUT_DIR}/per_base_coverage.bed.gz"
    "${OUTPUT_DIR}/per_base_coverage.bed.gz.csi"
    "${OUTPUT_DIR}/zero_low_coverage_intervals.bed.gz"
    "${OUTPUT_DIR}/zero_low_coverage_intervals.bed.gz.tbi"
    "${OUTPUT_DIR}/mosdepth.global.dist.txt"
    "${OUTPUT_DIR}/mosdepth.region.dist.txt"
    "${OUTPUT_DIR}/mosdepth.summary.txt"
    "${OUTPUT_DIR}/ont_qc_coverage_report.txt"
)

outputs_validate() {
    local output
    [[ -s "${COMPLETION_MARKER}" ]] || return 1
    [[ "${COMPLETION_MARKER}" -nt "${ALIGNMENT_BAM}" ]] || return 1
    [[ "${COMPLETION_MARKER}" -nt "${REFERENCE_FAI}" ]] || return 1
    for output in "${REQUIRED_OUTPUTS[@]}"; do
        [[ -s "${output}" ]] || return 1
    done
    "${TABIX_BIN}" -l "${OUTPUT_DIR}/window_coverage.tsv.gz" >/dev/null 2>&1 || return 1
    "${TABIX_BIN}" -l "${OUTPUT_DIR}/per_base_coverage.bed.gz" >/dev/null 2>&1 || return 1
    "${TABIX_BIN}" -l "${OUTPUT_DIR}/zero_low_coverage_intervals.bed.gz" >/dev/null 2>&1 || return 1
}

#---------------------------------------------->



#---------------------------------------------->
# Skips a current valid result unless --force requests a recoverable rebuild.

if outputs_validate && [[ ${FORCE} -ne 1 ]]; then
    log "Existing Stage 03 result is complete, current, and indexed; skipping rebuild"
    exit 0
fi
if outputs_validate && [[ ${FORCE} -eq 1 ]]; then
    log "A valid Stage 03 result exists, but --force requested a recoverable rebuild"
fi

#---------------------------------------------->



#---------------------------------------------->
# Refuses to overwrite incomplete outputs, or archives them first when --force was explicitly supplied.

if [[ -n "$(find "${OUTPUT_DIR}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    if [[ ${FORCE} -ne 1 ]]; then
        log "ERROR: Stage 03 contains incomplete, stale, or unvalidated outputs"
        log "Inspect ${OUTPUT_DIR}; rerun with --force to archive it and rebuild"
        exit 1
    fi
    STALE_ROOT="${TMP_DIR}/stale_stage03_outputs"
    STALE_DIR="${STALE_ROOT}/${PIPELINE_RUN_STAMP}"
    mkdir -p "${STALE_ROOT}"
    mv "${OUTPUT_DIR}" "${STALE_DIR}"
    mkdir -p "${OUTPUT_DIR}"
    log "Archived prior Stage 03 directory without deleting it: ${STALE_DIR}"
fi

#---------------------------------------------->



#---------------------------------------------->
# Creates a unique working directory so unfinished files remain separate from validated final outputs.

WORK_DIR="${TMP_DIR}/stage03.${PIPELINE_RUN_STAMP}"
MOSDEPTH_PREFIX="${WORK_DIR}/${SAMPLE_ID}"
mkdir -p "${WORK_DIR}"
log "Working directory: ${WORK_DIR}"

#---------------------------------------------->



#---------------------------------------------->
# Generates three complementary Samtools QC reports and confirms that two independent mapped-read counts agree.

log "Generating Samtools alignment statistics"
"${SAMTOOLS_BIN}" flagstat -@ "${QC_THREADS}" "${ALIGNMENT_BAM}" \
    > "${WORK_DIR}/alignment_flagstat.txt"
"${SAMTOOLS_BIN}" stats -@ "${QC_THREADS}" "${ALIGNMENT_BAM}" \
    > "${WORK_DIR}/alignment_stats.txt"
"${SAMTOOLS_BIN}" idxstats -@ "${QC_THREADS}" "${ALIGNMENT_BAM}" \
    > "${WORK_DIR}/alignment_idxstats.tsv"

flagstat_mapped="$(awk '/ mapped \(/ {print $1 + $3; exit}' "${WORK_DIR}/alignment_flagstat.txt")"
idxstats_mapped="$(awk '$1 != "*" {mapped += $3} END {printf "%.0f", mapped}' "${WORK_DIR}/alignment_idxstats.tsv")"
if [[ -z "${flagstat_mapped}" || "${flagstat_mapped}" != "${idxstats_mapped}" ]]; then
    log "ERROR: mapped-record counts disagree (flagstat=${flagstat_mapped:-missing}; idxstats=${idxstats_mapped:-missing})"
    exit 1
fi
log "Mapped-record count agrees between flagstat and idxstats: ${idxstats_mapped}"

#---------------------------------------------->



#---------------------------------------------->
# Runs Mosdepth once to calculate exact run-length depth, chromosome summaries, fixed windows, and breadth thresholds.

log "Running one Mosdepth pass for exact run-length depth, chromosome totals, and fixed windows"
log "RUN: MOSDEPTH_PRECISION=6 mosdepth --threads ${QC_THREADS} --by ${COVERAGE_WINDOW_SIZE} --thresholds ${COVERAGE_THRESHOLDS} --flag ${COVERAGE_EXCLUDE_FLAGS} --mapq ${COVERAGE_MIN_MAPQ} <prefix> <BAM>"
MOSDEPTH_PRECISION=6 "${MOSDEPTH_BIN}" \
    --threads "${QC_THREADS}" \
    --by "${COVERAGE_WINDOW_SIZE}" \
    --thresholds "${COVERAGE_THRESHOLDS}" \
    --flag "${COVERAGE_EXCLUDE_FLAGS}" \
    --mapq "${COVERAGE_MIN_MAPQ}" \
    "${MOSDEPTH_PREFIX}" \
    "${ALIGNMENT_BAM}"

#---------------------------------------------->



#---------------------------------------------->
# Names and verifies every expected Mosdepth file, including the indexed per-base run-length coverage track.

MOSDEPTH_REGIONS="${MOSDEPTH_PREFIX}.regions.bed.gz"
MOSDEPTH_THRESHOLDS="${MOSDEPTH_PREFIX}.thresholds.bed.gz"
MOSDEPTH_PER_BASE="${MOSDEPTH_PREFIX}.per-base.bed.gz"
MOSDEPTH_SUMMARY="${MOSDEPTH_PREFIX}.mosdepth.summary.txt"
for generated_file in \
    "${MOSDEPTH_REGIONS}" \
    "${MOSDEPTH_THRESHOLDS}" \
    "${MOSDEPTH_PER_BASE}" \
    "${MOSDEPTH_PER_BASE}.csi" \
    "${MOSDEPTH_SUMMARY}" \
    "${MOSDEPTH_PREFIX}.mosdepth.global.dist.txt" \
    "${MOSDEPTH_PREFIX}.mosdepth.region.dist.txt"; do
    if [[ ! -s "${generated_file}" ]]; then
        log "ERROR: Mosdepth did not produce expected file: ${generated_file}"
        exit 1
    fi
done
"${TABIX_BIN}" -l "${MOSDEPTH_PER_BASE}" >/dev/null

#---------------------------------------------->



#---------------------------------------------->
# Uses the Python helper to validate interval consistency and calculate genome, chromosome, window, and low-depth summaries.

log "Validating Mosdepth intervals and building coverage summaries"
python3 "${SUMMARY_HELPER}" \
    --fai "${REFERENCE_FAI}" \
    --regions "${MOSDEPTH_REGIONS}" \
    --thresholds-bed "${MOSDEPTH_THRESHOLDS}" \
    --per-base "${MOSDEPTH_PER_BASE}" \
    --mosdepth-summary "${MOSDEPTH_SUMMARY}" \
    --coverage-summary-out "${WORK_DIR}/coverage_summary.tsv" \
    --chromosome-out "${WORK_DIR}/chromosome_coverage.tsv" \
    --window-out "${WORK_DIR}/window_coverage.tsv" \
    --low-interval-out "${WORK_DIR}/zero_low_coverage_intervals.bed" \
    --thresholds "${COVERAGE_THRESHOLDS}" \
    --window-size "${COVERAGE_WINDOW_SIZE}" \
    --low-depth-threshold "${LOW_COVERAGE_DEPTH_THRESHOLD}" \
    --low-window-min-1x-breadth "${LOW_COVERAGE_MIN_1X_BREADTH_PERCENT}" \
    --exclude-flags "${COVERAGE_EXCLUDE_FLAGS}" \
    --min-mapq "${COVERAGE_MIN_MAPQ}"

#---------------------------------------------->



#---------------------------------------------->
# Compresses and indexes genomic tables, then gives retained Mosdepth outputs stable, descriptive filenames.

log "BGZF-compressing and indexing derived genomic tables"
"${BGZIP_BIN}" -@ "${QC_THREADS}" -c "${WORK_DIR}/window_coverage.tsv" \
    > "${WORK_DIR}/window_coverage.tsv.gz"
"${TABIX_BIN}" -f -p bed "${WORK_DIR}/window_coverage.tsv.gz"
"${BGZIP_BIN}" -@ "${QC_THREADS}" -c "${WORK_DIR}/zero_low_coverage_intervals.bed" \
    > "${WORK_DIR}/zero_low_coverage_intervals.bed.gz"
"${TABIX_BIN}" -f -p bed "${WORK_DIR}/zero_low_coverage_intervals.bed.gz"

"${TABIX_BIN}" -l "${WORK_DIR}/window_coverage.tsv.gz" >/dev/null
"${TABIX_BIN}" -l "${WORK_DIR}/zero_low_coverage_intervals.bed.gz" >/dev/null

mv "${MOSDEPTH_PER_BASE}" "${WORK_DIR}/per_base_coverage.bed.gz"
mv "${MOSDEPTH_PER_BASE}.csi" "${WORK_DIR}/per_base_coverage.bed.gz.csi"
mv "${MOSDEPTH_PREFIX}.mosdepth.global.dist.txt" "${WORK_DIR}/mosdepth.global.dist.txt"
mv "${MOSDEPTH_PREFIX}.mosdepth.region.dist.txt" "${WORK_DIR}/mosdepth.region.dist.txt"
cp "${MOSDEPTH_SUMMARY}" "${WORK_DIR}/mosdepth.summary.txt"
"${TABIX_BIN}" -l "${WORK_DIR}/per_base_coverage.bed.gz" >/dev/null

#---------------------------------------------->



#---------------------------------------------->
# Defines a helper that retrieves one named value from the calculated genome-wide coverage summary.

metric_value() {
    local metric="$1"
    awk -F '\t' -v metric="${metric}" '$1 == metric {print $2; exit}' "${WORK_DIR}/coverage_summary.tsv"
}

#---------------------------------------------->



#---------------------------------------------->
# Builds a human-readable report describing the coverage method, headline results, and interpretation notes.

REPORT_FILE="${WORK_DIR}/ont_qc_coverage_report.txt"
{
    printf 'ONT alignment QC and sequencing coverage report\n'
    printf 'Sample: %s\n' "${SAMPLE_ID}"
    printf 'Reference: %s\n' "${REFERENCE_NAME}"
    printf 'Input BAM: %s\n' "${ALIGNMENT_BAM}"
    printf 'Generated: %s\n\n' "$(date --iso-8601=seconds)"
    printf 'Coverage method\n'
    printf '  Mosdepth: %s\n' "$("${MOSDEPTH_BIN}" --version 2>&1 | head -n 1)"
    printf '  Window size: %s bases\n' "${COVERAGE_WINDOW_SIZE}"
    printf '  Included MAPQ: >=%s\n' "${COVERAGE_MIN_MAPQ}"
    printf '  Excluded SAM flag mask: %s (unmapped, secondary, QC-fail, duplicate)\n' "${COVERAGE_EXCLUDE_FLAGS}"
    printf '  Supplementary alignments are included, matching Mosdepth default flag policy.\n'
    printf '  Exact depth is retained as a BGZF run-length BED track; zero bases are explicit.\n\n'
    printf 'Headline coverage results\n'
    printf '  Reference bases: %s\n' "$(metric_value reference_bases)"
    printf '  Mean depth: %sx\n' "$(metric_value mean_depth)"
    for threshold in "${COVERAGE_THRESHOLD_ARRAY[@]}"; do
        printf '  Breadth >=%sx: %s%% (%s bases)\n' \
            "${threshold}" \
            "$(metric_value "pct_ge_${threshold}x")" \
            "$(metric_value "bases_ge_${threshold}x")"
    done
    printf '  Zero-coverage bases: %s (%s%%)\n' \
        "$(metric_value zero_coverage_bases)" \
        "$(metric_value pct_zero_coverage)"
    printf '  Exact zero-coverage intervals: %s\n' "$(metric_value zero_coverage_interval_count)"
    printf '  Exact low-depth intervals (1-%sx): %s\n' \
        "$((LOW_COVERAGE_DEPTH_THRESHOLD - 1))" \
        "$(metric_value low_depth_interval_count)"
    printf '  Windows with no coverage: %s\n' "$(metric_value zero_coverage_window_count)"
    printf '  Nonzero windows below %s%% 1x breadth: %s\n\n' \
        "${LOW_COVERAGE_MIN_1X_BREADTH_PERCENT}" \
        "$(metric_value low_1x_breadth_window_count)"
    printf 'Interpretation notes\n'
    printf '  Sequencing depth is not modified-base call depth. Stage 04 will report the\n'
    printf '  separate valid-call denominators needed for 5mC and 5hmC percentages.\n'
    printf '  The exact depth track and window table can later be queried at SixBase gap\n'
    printf '  regions, but no SixBase comparison is performed in this project.\n'
} > "${REPORT_FILE}"

#---------------------------------------------->



#---------------------------------------------->
# Confirms that every intended Stage 03 result exists and is nonempty before anything is promoted.

for final_work_file in \
    alignment_flagstat.txt \
    alignment_stats.txt \
    alignment_idxstats.tsv \
    coverage_summary.tsv \
    chromosome_coverage.tsv \
    window_coverage.tsv.gz \
    window_coverage.tsv.gz.tbi \
    per_base_coverage.bed.gz \
    per_base_coverage.bed.gz.csi \
    zero_low_coverage_intervals.bed.gz \
    zero_low_coverage_intervals.bed.gz.tbi \
    mosdepth.global.dist.txt \
    mosdepth.region.dist.txt \
    mosdepth.summary.txt \
    ont_qc_coverage_report.txt; do
    if [[ ! -s "${WORK_DIR}/${final_work_file}" ]]; then
        log "ERROR: validated Stage 03 output is missing or empty: ${final_work_file}"
        exit 1
    fi
done

#---------------------------------------------->



#---------------------------------------------->
# Moves all validated results from the temporary working directory into the official Stage 03 output directory.

log "Promoting validated outputs atomically within the 20 TB filesystem"
for final_work_file in \
    alignment_flagstat.txt \
    alignment_stats.txt \
    alignment_idxstats.tsv \
    coverage_summary.tsv \
    chromosome_coverage.tsv \
    window_coverage.tsv.gz \
    window_coverage.tsv.gz.tbi \
    per_base_coverage.bed.gz \
    per_base_coverage.bed.gz.csi \
    zero_low_coverage_intervals.bed.gz \
    zero_low_coverage_intervals.bed.gz.tbi \
    mosdepth.global.dist.txt \
    mosdepth.region.dist.txt \
    mosdepth.summary.txt \
    ont_qc_coverage_report.txt; do
    mv "${WORK_DIR}/${final_work_file}" "${OUTPUT_DIR}/${final_work_file}"
done

#---------------------------------------------->



#---------------------------------------------->
# Writes a timestamped reproducibility manifest with settings, tool versions, checksums, and output file sizes.

MANIFEST_FILE="${MANIFEST_DIR}/${SAMPLE_ID}.stage03.${PIPELINE_RUN_STAMP}.manifest.tsv"
{
    printf 'record_type\tname\tvalue\n'
    printf 'metadata\tsample_id\t%s\n' "${SAMPLE_ID}"
    printf 'metadata\treference\t%s\n' "${REFERENCE_NAME}"
    printf 'metadata\tinput_bam\t%s\n' "${ALIGNMENT_BAM}"
    printf 'metadata\tinput_bam_bytes\t%s\n' "$(stat -c '%s' "${ALIGNMENT_BAM}")"
    printf 'metadata\twindow_size\t%s\n' "${COVERAGE_WINDOW_SIZE}"
    printf 'metadata\tcoverage_thresholds\t%s\n' "${COVERAGE_THRESHOLDS}"
    printf 'metadata\tminimum_mapq\t%s\n' "${COVERAGE_MIN_MAPQ}"
    printf 'metadata\texclude_flag_mask\t%s\n' "${COVERAGE_EXCLUDE_FLAGS}"
    printf 'tool\tsamtools\t%s\n' "$("${SAMTOOLS_BIN}" --version | head -n 1)"
    printf 'tool\tmosdepth\t%s\n' "$("${MOSDEPTH_BIN}" --version 2>&1 | head -n 1)"
    printf 'checksum\tstage03_script_sha256\t%s\n' "$(sha256sum "$0" | awk '{print $1}')"
    printf 'checksum\tsummary_helper_sha256\t%s\n' "$(sha256sum "${SUMMARY_HELPER}" | awk '{print $1}')"
    printf 'checksum\tproject_config_sha256\t%s\n' "$(sha256sum "${PROJECT_ROOT}/config/project.env" | awk '{print $1}')"
    for output in "${REQUIRED_OUTPUTS[@]}"; do
        [[ "${output}" == "${COMPLETION_MARKER}" ]] && continue
        printf 'output_bytes\t%s\t%s\n' "$(basename -- "${output}")" "$(stat -c '%s' "${output}")"
    done
} > "${MANIFEST_FILE}.partial.${BASHPID}"
mv "${MANIFEST_FILE}.partial.${BASHPID}" "${MANIFEST_FILE}"
ln -sfn "$(basename -- "${MANIFEST_FILE}")" "${MANIFEST_DIR}/${SAMPLE_ID}.stage03.latest.manifest.tsv"

#---------------------------------------------->



#---------------------------------------------->
# Writes the completion marker last, then reruns the complete final-output validation.

{
    printf 'status\tcomplete\n'
    printf 'sample_id\t%s\n' "${SAMPLE_ID}"
    printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
    printf 'input_bam\t%s\n' "${ALIGNMENT_BAM}"
    printf 'manifest\t%s\n' "${MANIFEST_FILE}"
} > "${COMPLETION_MARKER}.partial.${BASHPID}"
mv "${COMPLETION_MARKER}.partial.${BASHPID}" "${COMPLETION_MARKER}"

if ! outputs_validate; then
    log "ERROR: final Stage 03 validation failed after promotion"
    exit 1
fi

#---------------------------------------------->



#---------------------------------------------->
# Safely removes only the expected temporary working directory after all final validation succeeds.

case "${WORK_DIR}" in
    "${TMP_DIR}"/stage03.*)
        rm -rf -- "${WORK_DIR}"
        ;;
    *)
        log "ERROR: refusing to clean unexpected working directory: ${WORK_DIR}"
        exit 1
        ;;
esac

#---------------------------------------------->



#---------------------------------------------->
# Reports successful completion and the locations of the principal Stage 03 results.

log "Stage 03 complete"
log "Coverage report: ${OUTPUT_DIR}/ont_qc_coverage_report.txt"
log "Coverage summary: ${OUTPUT_DIR}/coverage_summary.tsv"
log "Exact depth track: ${OUTPUT_DIR}/per_base_coverage.bed.gz"
log "Manifest: ${MANIFEST_FILE}"
#---------------------------------------------->
