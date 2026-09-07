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
# Checks for command line args, specifically just for help or checking whether required files exist and if they are valid.

CHECK_ONLY=0
case "${1:-}" in
    "") ;;
    --check-only) CHECK_ONLY=1 ;;
    -h|--help)
        printf 'Usage: %s [--check-only]\n' "$0"
        exit 0
        ;;
    *)
        printf 'Unknown argument: %s\n' "$1" >&2
        exit 2
        ;;
esac
#---------------------------------------------->



#---------------------------------------------->
# Sets paths and gets variable names for each

DORADO_BIN="${PROJECT_ROOT}/bin/dorado"
MINIMAP2_BIN="${PROJECT_ROOT}/bin/minimap2"
ALIGNMENT_DIR="${SAMPLE_OUTPUT}/01_alignment"
CHUNK_DIR="${ALIGNMENT_DIR}/aligned_unsorted_chunks"
MANIFEST_DIR="${SAMPLE_OUTPUT}/00_manifests"
LOG_DIR="${SAMPLE_OUTPUT}/logs"
STALE_DIR="${SAMPLE_OUTPUT}/tmp/stale_alignment_partials"
INPUT_BAMS="${INPUT_ROOT}/reads/${SAMPLE_ID}/bam_pass"
#---------------------------------------------->



#---------------------------------------------->
# Creates the required directories if they dont already exist, amd calls the 3 shared functions.

mkdir -p "${REFERENCE_DIR}" "${CHUNK_DIR}" "${MANIFEST_DIR}" \
    "${LOG_DIR}" "${STALE_DIR}"

# Starts a timestamped log so that if we get errors, we have a filled log and will know what happened.
start_pipeline_log "${LOG_DIR}" "01_align_modbam"

# Prevents two copies of script 1 from running at the same time for this sample
acquire_pipeline_lock "${SAMPLE_OUTPUT}/tmp/01_align_modbam.lock"

# Every 300 seconds outputs system load, memory, and free disk space for error checking purposes. (can remove later, but I like it.)
start_pipeline_heartbeat 300 "${DATA_ROOT}"
#---------------------------------------------->



#---------------------------------------------->
# Verify that we have the required programs

require_executable "${DORADO_BIN}"
require_executable "${MINIMAP2_BIN}"
require_executable "${SAMTOOLS_BIN}"
require_executable gzip
#---------------------------------------------->



#---------------------------------------------->
# Begin logging run information

log "Project: ${PROJECT_ROOT}"
log "Sample: ${SAMPLE_ID}"
log "Input BAM directory: ${INPUT_BAMS}"
log "Output chunk directory: ${CHUNK_DIR}"
log "Dorado: $("${DORADO_BIN}" --version 2>&1 | tail -n 1)"
log "Minimap2: $("${MINIMAP2_BIN}" --version 2>&1 | head -n 1)"
log "Samtools: $("${SAMTOOLS_BIN}" --version | head -n 1)"
log "Alignment threads: ${ALIGNMENT_THREADS}"
df -h "${DATA_ROOT}"
free -h
#---------------------------------------------->



#---------------------------------------------->
# --check-only argument behaivour. Just checks if each required file exists and is valid.

if [[ ${CHECK_ONLY} -eq 1 ]]; then
    if [[ ! -s "${REFERENCE_FASTA}" ]]; then
        log "CHECK: reference FASTA will be created from ${REFERENCE_SOURCE}"
    else
        log "CHECK: reference FASTA exists (${REFERENCE_FASTA})"
    fi
    if [[ ! -s "${REFERENCE_MMI}" ]]; then
        log "CHECK: correct lr:hq index will be built at ${REFERENCE_MMI}"
    else
        log "CHECK: lr:hq index exists (${REFERENCE_MMI})"
    fi
    if [[ ! -d "${INPUT_BAMS}" ]]; then
        log "ERROR: missing input BAM directory: ${INPUT_BAMS}"
        exit 1
    fi
    source_bam_count="$(find "${INPUT_BAMS}" -maxdepth 1 -type f -name '*.bam' | wc -l)"
    if [[ "${source_bam_count}" -ne "${EXPECTED_PASS_BAM_COUNT}" ]]; then
        log "ERROR: expected ${EXPECTED_PASS_BAM_COUNT} pass BAMs, found ${source_bam_count}"
        exit 1
    fi
    log "CHECK: found all ${source_bam_count} source pass BAMs"
    log "CHECK-ONLY complete; no alignment was started"
    exit 0
fi
#---------------------------------------------->



#---------------------------------------------->
# Decompress ref genome, from .gz to

if [[ ! -s "${REFERENCE_FASTA}" ]]; then
    reference_partial="${REFERENCE_FASTA}.partial.${BASHPID}"
    log "Decompressing reference to ${reference_partial}"
    gzip -dc "${REFERENCE_SOURCE}" > "${reference_partial}"
    mv "${reference_partial}" "${REFERENCE_FASTA}"
fi
#---------------------------------------------->



#---------------------------------------------->
# Build index(plural).

if [[ ! -s "${REFERENCE_FASTA}.fai" ]]; then
    log "Creating FASTA index ${REFERENCE_FASTA}.fai"
    "${SAMTOOLS_BIN}" faidx "${REFERENCE_FASTA}"
fi

if [[ ! -s "${REFERENCE_MMI}" ]]; then
    index_partial="${REFERENCE_MMI}.partial.${BASHPID}"
    log "Building Minimap2 lr:hq index ${index_partial}"
    "${MINIMAP2_BIN}" -x lr:hq -d "${index_partial}" "${REFERENCE_FASTA}"
    mv "${index_partial}" "${REFERENCE_MMI}"
fi
#---------------------------------------------->



#---------------------------------------------->
# Discover and validate input list

if [[ ! -d "${INPUT_BAMS}" ]]; then
    log "ERROR: missing input BAM directory: ${INPUT_BAMS}"
    exit 1
fi

mapfile -t SOURCE_BAMS < <(find "${INPUT_BAMS}" -maxdepth 1 -type f -name '*.bam' | sort -V)
if [[ ${#SOURCE_BAMS[@]} -ne ${EXPECTED_PASS_BAM_COUNT} ]]; then
    log "ERROR: expected ${EXPECTED_PASS_BAM_COUNT} pass BAMs, found ${#SOURCE_BAMS[@]}"
    exit 1
fi
#---------------------------------------------->



#---------------------------------------------->
# Creates a list of each source BAM file and its size.

manifest_partial="${MANIFEST_DIR}/${SAMPLE_ID}.source_pass_bams.tsv.partial.${BASHPID}"
{
    printf 'source_bam\tbytes\n'
    for source_bam in "${SOURCE_BAMS[@]}"; do
        printf '%s\t%s\n' "${source_bam}" "$(stat -c '%s' "${source_bam}")"
    done
} > "${manifest_partial}"
mv "${manifest_partial}" "${MANIFEST_DIR}/${SAMPLE_ID}.source_pass_bams.tsv"
#---------------------------------------------->



#---------------------------------------------->
# Finds partial BAMs left behind by interrupted runs and moves them to a timestamped archive instead of deleting them.

mapfile -t STALE_PARTIALS < <(find "${CHUNK_DIR}" -maxdepth 1 -type f -name '*.partial.*.bam' | sort -V)
if [[ ${#STALE_PARTIALS[@]} -gt 0 ]]; then
    stale_run_dir="${STALE_DIR}/${PIPELINE_RUN_STAMP}"
    mkdir -p "${stale_run_dir}"
    log "Archiving ${#STALE_PARTIALS[@]} stale partial BAM(s) to ${stale_run_dir}"
    mv "${STALE_PARTIALS[@]}" "${stale_run_dir}/"
fi

#---------------------------------------------->



#---------------------------------------------->
# Counters for completed, skipped, and total chunks.

completed=0
skipped=0
total=${#SOURCE_BAMS[@]}

#---------------------------------------------->



#---------------------------------------------->
# Function that defines a BAM and validates it.

FULL_VALIDATED_RECORD_COUNT=""
full_validate_aligned_bam() {
    local bam="$1"
    local record_count
    "${SAMTOOLS_BIN}" quickcheck -v "${bam}" || return 1
    record_count="$("${SAMTOOLS_BIN}" view -@ "${VALIDATION_THREADS}" -c "${bam}")" || return 1
    if ! [[ "${record_count}" =~ ^[1-9][0-9]*$ ]]; then
        log "ERROR: full-stream validation found no records in ${bam}"
        return 1
    fi
    FULL_VALIDATED_RECORD_COUNT="${record_count}"
}

#---------------------------------------------->



#---------------------------------------------->
# Function that helps us log completed BAM chunks with corresponding info.

write_chunk_marker() {
    local marker="$1"
    local bam="$2"
    local record_count="$3"
    local marker_partial="${marker}.partial.${PIPELINE_RUN_STAMP}"
    {
        printf 'status\tcomplete\n'
        printf 'validation_method\tsamtools_quickcheck_plus_full_view_count_v1\n'
        printf 'record_count\t%s\n' "${record_count}"
        printf 'bam_bytes\t%s\n' "$(stat -c '%s' "${bam}")"
        printf 'bam_mtime_epoch\t%s\n' "$(stat -c '%Y' "${bam}")"
        printf 'validated_at\t%s\n' "$(date --iso-8601=seconds)"
    } > "${marker_partial}"
    mv "${marker_partial}" "${marker}"
}

#---------------------------------------------->



#---------------------------------------------->
# Function that checks if a BAM file matches its log from the above function.

marker_matches_bam() {
    local marker="$1"
    local bam="$2"
    local marker_bytes marker_mtime marker_method
    [[ -s "${marker}" ]] || return 1
    marker_method="$(awk -F '\t' '$1 == "validation_method" {print $2}' "${marker}")"
    marker_bytes="$(awk -F '\t' '$1 == "bam_bytes" {print $2}' "${marker}")"
    marker_mtime="$(awk -F '\t' '$1 == "bam_mtime_epoch" {print $2}' "${marker}")"
    [[ "${marker_method}" == "samtools_quickcheck_plus_full_view_count_v1" \
        && "${marker_bytes}" == "$(stat -c '%s' "${bam}")" \
        && "${marker_mtime}" == "$(stat -c '%Y' "${bam}")" ]]
}

#---------------------------------------------->



#---------------------------------------------->
# MAIN ALIGNMENT LOOP: Aligns each source BAM to the ref genome and validates before moving onto the next.

# Primary loop, we begin processing each BAM stored in our SOURCE_BAMS array.
for source_bam in "${SOURCE_BAMS[@]}"; do
    source_name="$(basename -- "${source_bam}" .bam)"
    output_bam="${CHUNK_DIR}/${source_name}.aligned.unsorted.bam"
    complete_marker="${output_bam}.complete"

    # Reuses an existing BAM when it is nonempty, its marker still matches it, and Samtools confirms its validity.
    if [[ -s "${output_bam}" ]] \
        && marker_matches_bam "${complete_marker}" "${output_bam}" \
        && "${SAMTOOLS_BIN}" quickcheck "${output_bam}"; then
        skipped=$((skipped + 1))
        log "SKIP [$((completed + skipped))/${total}]: fully validated existing $(basename -- "${output_bam}")"
        continue
    fi

    # If BAM exists already, but marker is invalid then we remake it and validate.
    if [[ -s "${output_bam}" ]]; then
        log "FULL CHECK [$((completed + skipped + 1))/${total}]: legacy or changed $(basename -- "${output_bam}")"
        if full_validate_aligned_bam "${output_bam}"; then
            write_chunk_marker "${complete_marker}" "${output_bam}" "${FULL_VALIDATED_RECORD_COUNT}"
            skipped=$((skipped + 1))
            log "SKIP [$((completed + skipped))/${total}]: full stream passed (${FULL_VALIDATED_RECORD_COUNT} records)"
            continue
        fi

        # If we fail to remake the marker, move the possibly corrupted BAM into quarentine.
        invalid_dir="${STALE_DIR}/${PIPELINE_RUN_STAMP}/invalid_completed_chunks"
        mkdir -p "${invalid_dir}"
        log "QUARANTINE: full-stream validation failed for ${output_bam}"
        mv "${output_bam}" "${invalid_dir}/"
        if [[ -e "${complete_marker}" ]]; then
            mv "${complete_marker}" "${invalid_dir}/"
        fi
        log "Archived invalid chunk without deleting it: ${invalid_dir}"
    fi

    # Aligns the current modBAM to the ref genome with dorado.
    output_partial="${output_bam}.partial.${BASHPID}.bam"
    log "START [$((completed + skipped + 1))/${total}]: $(basename -- "${source_bam}")"
    log "RUN: nice -n 5 dorado aligner --no-sort -t ${ALIGNMENT_THREADS} --mm2-opts '-x lr:hq -Y' <mmi> <bam>"

    nice -n 5 "${DORADO_BIN}" aligner \
        --no-sort \
        -t "${ALIGNMENT_THREADS}" \
        --mm2-opts "-x lr:hq -Y" \
        "${REFERENCE_MMI}" \
        "${source_bam}" \
        > "${output_partial}"

    # Reads and validates the entire newly aligned BAM before allowing it to receive its final filename.
    log "FULL CHECK: decompressing the completed partial BAM before promotion"
    if ! full_validate_aligned_bam "${output_partial}"; then
        log "ERROR: newly aligned BAM failed full-stream validation: ${output_partial}"
        exit 1
    fi

    # Checks the first alignment record to confirm that the MM, ML, and MN modified-base tags survived alignment.
    first_record="$("${SAMTOOLS_BIN}" view "${output_partial}" 2>/dev/null | head -n 1 || true)"
    if [[ "${first_record}" != *$'\tMM:Z:'* \
        || "${first_record}" != *$'\tML:B:C'* \
        || "${first_record}" != *$'\tMN:i:'* ]]; then
        log "ERROR: newly aligned BAM lacks MM/ML/MN tags: ${output_partial}"
        exit 1
    fi

    # Promotes the validated temporary BAM to its final name, writes its completion marker, and updates progress.
    mv "${output_partial}" "${output_bam}"
    write_chunk_marker "${complete_marker}" "${output_bam}" "${FULL_VALIDATED_RECORD_COUNT}"
    completed=$((completed + 1))
    log "DONE [$((completed + skipped))/${total}]: $(basename -- "${output_bam}") ($(stat -c '%s' "${output_bam}") bytes; ${FULL_VALIDATED_RECORD_COUNT} records)"
done

#---------------------------------------------->



#---------------------------------------------->
# Reports finall completion, then outputs script 02's location and the next step.

log "Alignment chunk stage complete: newly completed=${completed}, previously complete=${skipped}, total=${total}"
log "Next step: ${PROJECT_ROOT}/scripts/02_finalize_alignment.sh"
#---------------------------------------------->
