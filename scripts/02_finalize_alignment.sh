#!/usr/bin/env bash
set -Eeuo pipefail


#---------------------------------------------->
# Initializes the project, finding where the existing files are.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=../config/project.env
source "${PROJECT_ROOT}/config/project.env"
# shellcheck source=lib/pipeline_common.sh
source "${SCRIPT_DIR}/lib/pipeline_common.sh"

#---------------------------------------------->



#---------------------------------------------->
# Sets the locations of the Stage 01 chunks, Stage 02 outputs, log files, etc.

INPUT_BAM_DIR="${INPUT_ROOT}/reads/${SAMPLE_ID}/bam_pass"
CHUNK_DIR="${SAMPLE_OUTPUT}/01_alignment/aligned_unsorted_chunks"
MANIFEST_DIR="${SAMPLE_OUTPUT}/00_manifests"
QC_DIR="${SAMPLE_OUTPUT}/02_alignment_qc"
TMP_DIR="${SAMPLE_OUTPUT}/tmp"
LOG_DIR="${SAMPLE_OUTPUT}/logs"

#---------------------------------------------->



#---------------------------------------------->
# Creates the required directories, starts the log and resource heartbeat, prevents duplicate runs, and checks Samtools.

mkdir -p "${MANIFEST_DIR}" "${QC_DIR}" "${TMP_DIR}" "${LOG_DIR}"
start_pipeline_log "${LOG_DIR}" "02_finalize_alignment"
acquire_pipeline_lock "${TMP_DIR}/02_finalize_alignment.lock"
start_pipeline_heartbeat 300 "${DATA_ROOT}"
require_executable "${SAMTOOLS_BIN}"

#---------------------------------------------->



#---------------------------------------------->
# Names the final sorted BAM, its index, the temporary BAM list, and the Stage 02 reports.

ALIGNMENT_BAM="${QC_DIR}/${SAMPLE_ID}.aligned.sorted.bam"
ALIGNMENT_BAI="${ALIGNMENT_BAM}.bai"
BAM_LIST_FILE="${TMP_DIR}/${SAMPLE_ID}.aligned_bams.list"
MANIFEST_FILE="${MANIFEST_DIR}/${SAMPLE_ID}.alignment_inputs.tsv"
QC_SUMMARY="${QC_DIR}/${SAMPLE_ID}.alignment_summary.txt"

#---------------------------------------------->



#---------------------------------------------->
# Records the sample ID, Aligned chunk directory, final BAMs, etc.

log "Sample: ${SAMPLE_ID}"
log "Aligned chunk directory: ${CHUNK_DIR}"
log "Final aligned BAM: ${ALIGNMENT_BAM}"
log "Samtools: $("${SAMTOOLS_BIN}" --version | head -n 1)"
log "Sort threads: ${SORT_THREADS}; memory per sort thread: ${SORT_MEMORY_PER_THREAD}"
df -h "${DATA_ROOT}"
free -h

#---------------------------------------------->



#---------------------------------------------->
# Stops if the aligned chunk directory from Stage 01 does not exist.

if [[ ! -d "${CHUNK_DIR}" ]]; then
    log "ERROR: aligned chunk directory does not exist: ${CHUNK_DIR}"
    log "Run ${PROJECT_ROOT}/scripts/01_align_modbam.sh first"
    exit 1
fi

#---------------------------------------------->



#---------------------------------------------->
# Discovers the original source BAMs and confirms that the expected number of source files is present.

mapfile -t SOURCE_BAMS < <(find "${INPUT_BAM_DIR}" -maxdepth 1 -type f -name '*.bam' | sort -V)
if [[ ${#SOURCE_BAMS[@]} -ne ${EXPECTED_PASS_BAM_COUNT} ]]; then
    log "ERROR: expected ${EXPECTED_PASS_BAM_COUNT} source BAMs, found ${#SOURCE_BAMS[@]}"
    exit 1
fi

#---------------------------------------------->



#---------------------------------------------->
# Matches every source BAM to its expected aligned chunk and stops if any chunk is missing or unexpected.

INPUT_BAMS=()
missing=0
for source_bam in "${SOURCE_BAMS[@]}"; do
    source_name="$(basename -- "${source_bam}" .bam)"
    aligned_bam="${CHUNK_DIR}/${source_name}.aligned.unsorted.bam"
    if [[ ! -s "${aligned_bam}" ]]; then
        if [[ ${missing} -lt 10 ]]; then
            log "MISSING: ${aligned_bam}"
        fi
        missing=$((missing + 1))
        continue
    fi
    INPUT_BAMS+=("${aligned_bam}")
done

mapfile -t DISCOVERED_BAMS < <(find "${CHUNK_DIR}" -maxdepth 1 -type f -name '*.aligned.unsorted.bam' | sort -V)
if [[ ${missing} -ne 0 || ${#INPUT_BAMS[@]} -ne ${EXPECTED_PASS_BAM_COUNT} || ${#DISCOVERED_BAMS[@]} -ne ${EXPECTED_PASS_BAM_COUNT} ]]; then
    if [[ ${missing} -gt 10 ]]; then
        log "MISSING: ${missing} total chunks (first 10 shown above)"
    fi
    log "ERROR: alignment is incomplete or contains unexpected chunks."
    log "Expected=${EXPECTED_PASS_BAM_COUNT}; matched=${#INPUT_BAMS[@]}; discovered=${#DISCOVERED_BAMS[@]}; missing=${missing}"
    exit 1
fi

#---------------------------------------------->



#---------------------------------------------->
# Fully reads every aligned chunk, counts its records, and writes a validation log for each chunk before sorting.

log "Fully decompressing all ${#INPUT_BAMS[@]} aligned chunks before the expensive sort"
CHUNK_VALIDATION_FILE="${MANIFEST_DIR}/${SAMPLE_ID}.alignment_chunk_validation.tsv"
chunk_validation_partial="${CHUNK_VALIDATION_FILE}.partial.${PIPELINE_RUN_STAMP}"
printf 'chunk_index\tbam\tbytes\trecord_count\tvalidation\n' > "${chunk_validation_partial}"
total_input_records=0
chunk_index=0
for input_bam in "${INPUT_BAMS[@]}"; do
    chunk_index=$((chunk_index + 1))
    log "FULL CHECK [${chunk_index}/${#INPUT_BAMS[@]}]: $(basename -- "${input_bam}")"
    "${SAMTOOLS_BIN}" quickcheck -v "${input_bam}"
    if ! record_count="$("${SAMTOOLS_BIN}" view -@ "${VALIDATION_THREADS}" -c "${input_bam}")"; then
        log "ERROR: full-stream BAM validation failed: ${input_bam}"
        log "Rerun Stage 01; it will quarantine and regenerate the invalid chunk"
        exit 1
    fi
    if ! [[ "${record_count}" =~ ^[1-9][0-9]*$ ]]; then
        log "ERROR: aligned chunk has zero or invalid record count: ${input_bam}"
        exit 1
    fi
    total_input_records=$((total_input_records + record_count))
    printf '%d\t%s\t%s\t%s\tfull_stream_pass\n' \
        "${chunk_index}" "${input_bam}" "$(stat -c '%s' "${input_bam}")" "${record_count}" \
        >> "${chunk_validation_partial}"
done
mv "${chunk_validation_partial}" "${CHUNK_VALIDATION_FILE}"
log "All chunks passed full-stream validation: total_records=${total_input_records}"

#---------------------------------------------->



#---------------------------------------------->
# Checks the first record of every chunk to confirm that the modified-base tags are correct.

log "Checking MM, ML, and MN modified-base tags in every aligned chunk"
for input_bam in "${INPUT_BAMS[@]}"; do
    first_record="$("${SAMTOOLS_BIN}" view "${input_bam}" 2>/dev/null | head -n 1 || true)"
    if [[ "${first_record}" != *$'\tMM:Z:'* \
        || "${first_record}" != *$'\tML:B:C'* \
        || "${first_record}" != *$'\tMN:i:'* ]]; then
        log "ERROR: MM/ML/MN modified-base tags were not found in ${input_bam}"
        exit 1
    fi
done
log "Confirmed MM, ML, and MN modified-base tags in all aligned chunks"

#---------------------------------------------->



#---------------------------------------------->
# Writes the ordered chunk list used by Samtools and checks to make sure each sorted chunk has a new timestamp compared to its older, unsorted chunk.

printf '%s\n' "${INPUT_BAMS[@]}" > "${BAM_LIST_FILE}"

final_is_current=1
for input_bam in "${INPUT_BAMS[@]}"; do
    if [[ "${input_bam}" -nt "${ALIGNMENT_BAM}" ]]; then
        final_is_current=0
        break
    fi
done

#---------------------------------------------->



#---------------------------------------------->
# Fully validates a current final BAM and only reuses it when its record count exactly matches the sum of all chunks.

existing_final_valid=0
if [[ ${final_is_current} -eq 1 ]] \
    && [[ -s "${ALIGNMENT_BAM}" ]] \
    && [[ -s "${ALIGNMENT_BAI}" ]] \
    && "${SAMTOOLS_BIN}" quickcheck "${ALIGNMENT_BAM}" \
    && "${SAMTOOLS_BIN}" idxstats "${ALIGNMENT_BAM}" >/dev/null; then
    log "Fully validating existing final BAM before accepting it"
    if existing_final_record_count="$("${SAMTOOLS_BIN}" view -@ "${VALIDATION_THREADS}" -c "${ALIGNMENT_BAM}")" \
        && [[ "${existing_final_record_count}" == "${total_input_records}" ]]; then
        existing_final_valid=1
    else
        log "WARNING: existing final BAM failed full validation or record-count comparison"
    fi
fi

#---------------------------------------------->



#---------------------------------------------->
# Reuses a valid final BAM or concatenates all chunks, coordinate-sorts once, validates the result, and builds its index.

if [[ ${existing_final_valid} -eq 1 ]]; then
    log "Existing final BAM and index validate successfully; skipping rebuild"
else
    final_partial="${QC_DIR}/${SAMPLE_ID}.aligned.sorted.${PIPELINE_RUN_STAMP}.partial.bam"
    sort_prefix="${TMP_DIR}/${SAMPLE_ID}.samtools_sort.${PIPELINE_RUN_STAMP}.${BASHPID}"

    log "Concatenating ${#INPUT_BAMS[@]} unsorted chunks and coordinate-sorting once"
    log "RUN: samtools cat -b <list> | samtools sort -@ ${SORT_THREADS} -m ${SORT_MEMORY_PER_THREAD}"
    "${SAMTOOLS_BIN}" cat -b "${BAM_LIST_FILE}" \
        | "${SAMTOOLS_BIN}" sort \
            -@ "${SORT_THREADS}" \
            -m "${SORT_MEMORY_PER_THREAD}" \
            -T "${sort_prefix}" \
            -o "${final_partial}" \
            -

    "${SAMTOOLS_BIN}" quickcheck -v "${final_partial}"
    log "Fully validating sorted partial BAM before promotion"
    final_record_count="$("${SAMTOOLS_BIN}" view -@ "${VALIDATION_THREADS}" -c "${final_partial}")"
    if [[ "${final_record_count}" != "${total_input_records}" ]]; then
        log "ERROR: final BAM record count (${final_record_count}) does not match chunk total (${total_input_records})"
        exit 1
    fi
    sort_order="$("${SAMTOOLS_BIN}" view -H "${final_partial}" | awk -F '\t' '$1 == "@HD" {for (i=1; i<=NF; i++) if ($i == "SO:coordinate") found=1} END {print found+0}')"
    if [[ "${sort_order}" -ne 1 ]]; then
        log "ERROR: finalized BAM header is not coordinate sorted"
        exit 1
    fi

    mv "${final_partial}" "${ALIGNMENT_BAM}"
    index_partial="${ALIGNMENT_BAI}.partial.${BASHPID}"
    "${SAMTOOLS_BIN}" index -@ "${SORT_THREADS}" -o "${index_partial}" "${ALIGNMENT_BAM}"
    mv "${index_partial}" "${ALIGNMENT_BAI}"
fi

#---------------------------------------------->



#---------------------------------------------->
# Writes a manifest describing the sample, inputs, validation record counts, and final BAM and index.

manifest_partial="${MANIFEST_FILE}.partial.${BASHPID}"
{
    printf 'sample_id\t%s\n' "${SAMPLE_ID}"
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'alignment_chunk_dir\t%s\n' "${CHUNK_DIR}"
    printf 'input_bam_count\t%d\n' "${#INPUT_BAMS[@]}"
    printf 'input_record_count\t%s\n' "${total_input_records}"
    printf 'chunk_validation_manifest\t%s\n' "${CHUNK_VALIDATION_FILE}"
    printf 'merged_bam\t%s\n' "${ALIGNMENT_BAM}"
    printf 'merged_bai\t%s\n' "${ALIGNMENT_BAI}"
} > "${manifest_partial}"
mv "${manifest_partial}" "${MANIFEST_FILE}"

#---------------------------------------------->



#---------------------------------------------->
# Creates an alignment summary containing the BAM header overview and Samtools flag statistics.

summary_partial="${QC_SUMMARY}.partial.${BASHPID}"
{
    printf 'Aligned BAM summary for %s\n' "${SAMPLE_ID}"
    printf 'Input chunks: %d\n' "${#INPUT_BAMS[@]}"
    printf 'Merged BAM: %s\n' "${ALIGNMENT_BAM}"
    printf 'Indexed BAM: %s\n' "${ALIGNMENT_BAI}"
    printf '\nHeader overview:\n'
    "${SAMTOOLS_BIN}" view -H "${ALIGNMENT_BAM}" | sed -n '1,40p'
    printf '\nFlagstat:\n'
    "${SAMTOOLS_BIN}" flagstat -@ "${SORT_THREADS}" "${ALIGNMENT_BAM}"
} > "${summary_partial}"
mv "${summary_partial}" "${QC_SUMMARY}"

#---------------------------------------------->



#---------------------------------------------->
# Reports successful completion and the locations of the final Stage 02 outputs.

log "Alignment finalization complete"
log "Manifest: ${MANIFEST_FILE}"
log "Summary: ${QC_SUMMARY}"
log "Final BAM: ${ALIGNMENT_BAM}"
#---------------------------------------------->
