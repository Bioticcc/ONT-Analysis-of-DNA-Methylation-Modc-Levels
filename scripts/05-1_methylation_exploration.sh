#!/usr/bin/env bash
set -Eeuo pipefail


#---------------------------------------------->
# Initializes the project and loads shared paths, settings, logging, locking, and validation functions.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=../config/project.env
source "${PROJECT_ROOT}/config/project.env"
# shellcheck source=lib/pipeline_common.sh
source "${SCRIPT_DIR}/lib/pipeline_common.sh"

#---------------------------------------------->



#---------------------------------------------->
# Defines Stage 05 usage and its non-destructive check-only and recoverable forced-rebuild modes.

usage() {
    cat <<'EOF'
Usage: 05-1_methylation_exploration.sh [--check-only] [--force]

Create single-sample ONT CpG QC, genome-wide, regional, and genomic-feature
figures with exact plot-data tables. Compatible SixBase figure families are
reproduced; multi-sample PCA, correlation, DMR, and enrichment plots are
explicitly documented as unavailable for one sample.

  --check-only  Validate Stage 03/04 results, annotations, R, and packages only.
  --force       Archive existing Stage 05 outputs and rebuild them.
  -h, --help    Show this help text.
EOF
}

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
# Declares Stage 03/04 inputs, Stage 05 outputs, annotations, logs, manifests, and the R implementation.

STAGE03_DIR="${SAMPLE_OUTPUT}/03_ont_qc_coverage"
STAGE04_DIR="${SAMPLE_OUTPUT}/04_methylation"
OUTPUT_DIR="${SAMPLE_OUTPUT}/05_methylation_exploration"
MANIFEST_DIR="${SAMPLE_OUTPUT}/00_manifests"
TMP_DIR="${SAMPLE_OUTPUT}/tmp"
LOG_DIR="${SAMPLE_OUTPUT}/logs"
R_SCRIPT="${SCRIPT_DIR}/05-2_methylation_exploration.R"
REFERENCE_FAI="${REFERENCE_FASTA}.fai"

STAGE03_MARKER="${STAGE03_DIR}/.stage03.complete"
STAGE04_MARKER="${STAGE04_DIR}/.stage04.complete"
GLOBAL_LEVELS="${STAGE04_DIR}/global_methylation_levels.tsv"
CHROMOSOME_METHYLATION="${STAGE04_DIR}/chromosome_methylation_levels.tsv"
CALL_COVERAGE="${STAGE04_DIR}/methylation_call_coverage.tsv"
WINDOW_METHYLATION="${STAGE04_DIR}/window_methylation_levels.tsv.gz"
CPG_BEDMETHYL="${STAGE04_DIR}/cpg_5mc_5hmc.bed.gz"
CHROMOSOME_COVERAGE="${STAGE03_DIR}/chromosome_coverage.tsv"
WINDOW_COVERAGE="${STAGE03_DIR}/window_coverage.tsv.gz"
ALIGNMENT_STATS="${STAGE03_DIR}/alignment_stats.txt"
ALIGNMENT_FLAGSTAT="${STAGE03_DIR}/alignment_flagstat.txt"
COVERAGE_SUMMARY="${STAGE03_DIR}/coverage_summary.tsv"
STAGE02_SUMMARY="${SAMPLE_OUTPUT}/02_alignment_qc/${SAMPLE_ID}.alignment_summary.txt"

shopt -s nullglob
ONT_RUN_REPORT_CANDIDATES=("${RUN_ROOT}"/report_*.html)
shopt -u nullglob
if [[ ${#ONT_RUN_REPORT_CANDIDATES[@]} -ne 1 ]]; then
    log "ERROR: expected exactly one ONT instrument HTML report in ${RUN_ROOT}; found ${#ONT_RUN_REPORT_CANDIDATES[@]}"
    exit 1
fi
ONT_RUN_REPORT="${ONT_RUN_REPORT_CANDIDATES[0]}"

#---------------------------------------------->



#---------------------------------------------->
# Creates the output roots, starts the standard log and heartbeat, locks Stage 05, and verifies executables.

mkdir -p "${OUTPUT_DIR}" "${MANIFEST_DIR}" "${TMP_DIR}" "${LOG_DIR}"
start_pipeline_log "${LOG_DIR}" "05_methylation_exploration"
acquire_pipeline_lock "${TMP_DIR}/05_methylation_exploration.lock"
start_pipeline_heartbeat 300 "${DATA_ROOT}"

require_executable "${RSCRIPT_BIN}"
require_executable "${TABIX_BIN}"
require_executable gzip
require_executable awk
require_executable sha256sum

if [[ ! -r "${R_SCRIPT}" ]]; then
    log "ERROR: missing Stage 05 R script: ${R_SCRIPT}"
    exit 1
fi

#---------------------------------------------->



#---------------------------------------------->
# Validates visualization thresholds and confirms every required R package is installed in the configured environment.

if ! [[ "${EXPLORATION_MIN_VALID_COVERAGE}" =~ ^[1-9][0-9]*$ \
    && "${EXPLORATION_MIN_FEATURE_CPGS}" =~ ^[1-9][0-9]*$ \
    && "${EXPLORATION_TOP_WINDOW_COUNT}" =~ ^[1-9][0-9]*$ \
    && "${EXPLORATION_PLOT_DPI}" =~ ^[1-9][0-9]*$ ]]; then
    log "ERROR: invalid numeric Stage 05 configuration"
    exit 1
fi
if (( EXPLORATION_PLOT_DPI < 72 )); then
    log "ERROR: EXPLORATION_PLOT_DPI must be at least 72"
    exit 1
fi

required_r_packages=(
    data.table
    ggplot2
    ggrepel
    scales
    patchwork
    GenomicRanges
    IRanges
    GenomeInfoDb
    rtracklayer
    jsonlite
)
r_package_list="$(IFS=,; printf '%s' "${required_r_packages[*]}")"
"${RSCRIPT_BIN}" -e '
packages <- strsplit(commandArgs(TRUE)[1], ",", fixed = TRUE)[[1]]
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("missing R package(s): ", paste(missing, collapse = ", "))
' "${r_package_list}"

#---------------------------------------------->



#---------------------------------------------->
# Requires completed, indexed Stage 03/04 products and the exact annotation files already used by SixBase.

required_inputs=(
    "${STAGE03_MARKER}"
    "${STAGE04_MARKER}"
    "${REFERENCE_FAI}"
    "${GLOBAL_LEVELS}"
    "${CHROMOSOME_METHYLATION}"
    "${CALL_COVERAGE}"
    "${WINDOW_METHYLATION}"
    "${WINDOW_METHYLATION}.tbi"
    "${CPG_BEDMETHYL}"
    "${CPG_BEDMETHYL}.tbi"
    "${CHROMOSOME_COVERAGE}"
    "${WINDOW_COVERAGE}"
    "${WINDOW_COVERAGE}.tbi"
    "${ALIGNMENT_STATS}"
    "${ALIGNMENT_FLAGSTAT}"
    "${COVERAGE_SUMMARY}"
    "${STAGE02_SUMMARY}"
    "${ONT_RUN_REPORT}"
    "${METHYLATION_GENCODE_GFF3}"
    "${METHYLATION_CPG_ISLANDS_JSON}"
    "${METHYLATION_CCRE_TABLE}"
    "${METHYLATION_INTERGENIC_BED}"
)
for required_input in "${required_inputs[@]}"; do
    if [[ ! -s "${required_input}" ]]; then
        log "ERROR: required Stage 05 input is missing or empty: ${required_input}"
        exit 1
    fi
done
if ! grep -q $'^status\tcomplete$' "${STAGE03_MARKER}" \
    || ! grep -q $'^status\tcomplete$' "${STAGE04_MARKER}"; then
    log "ERROR: Stage 03 and Stage 04 must both have complete status"
    exit 1
fi
"${TABIX_BIN}" -l "${WINDOW_COVERAGE}" >/dev/null
"${TABIX_BIN}" -l "${WINDOW_METHYLATION}" >/dev/null
"${TABIX_BIN}" -l "${CPG_BEDMETHYL}" >/dev/null

#---------------------------------------------->



#---------------------------------------------->
# Records the exact analysis scope, thresholds, annotation provenance, and R environment in the log.

log "Sample: ${SAMPLE_ID}"
log "Stage 03 input: ${STAGE03_DIR}"
log "Stage 04 input: ${STAGE04_DIR}"
log "ONT instrument report: ${ONT_RUN_REPORT}"
log "Output directory: ${OUTPUT_DIR}"
log "R: $("${RSCRIPT_BIN}" --version 2>&1 | head -n 1)"
log "Feature thresholds: valid_call_depth>=${EXPLORATION_MIN_VALID_COVERAGE}; CpGs_per_plotted_region>=${EXPLORATION_MIN_FEATURE_CPGS}"
log "Extreme-window table size per ranking: ${EXPLORATION_TOP_WINDOW_COUNT}"
log "Figure format: PDF (the HTML dashboard embeds relative PDF previews)"
log "GENCODE annotation: ${METHYLATION_GENCODE_GFF3}"
log "CpG islands: ${METHYLATION_CPG_ISLANDS_JSON}"
log "cCREs: ${METHYLATION_CCRE_TABLE}"
log "Intergenic regions: ${METHYLATION_INTERGENIC_BED}"
log "Single-sample boundary: correlations, PCA, DMR statistics, DMR heatmaps, and enrichment are not generated"

if [[ ${CHECK_ONLY} -eq 1 ]]; then
    log "CHECK-ONLY complete: Stage 05 inputs, indexes, annotations, R, and packages validate successfully"
    exit 0
fi

#---------------------------------------------->



#---------------------------------------------->
# Lists the complete Stage 05 contract and defines checks for freshness, figure formats, tables, and report content.

COMPLETION_MARKER="${OUTPUT_DIR}/.stage05.complete"
FIGURE_STEMS=(
    01_Global_CpG_Modification_Composition
    02_Global_CpG_Methylation_Levels
    03_CpG_Valid_Coverage_Distribution
    04_CpG_Coverage_Thresholds
    05_CpG_Site_Methylation_Distributions
    06_Chromosome_CpG_Methylation_Levels
    07_Chromosome_Sequencing_and_Call_Coverage
    08_Sequencing_vs_Modification_Call_Coverage
    09_Genomewide_5mC_100kb_Windows
    10_Genomewide_5hmC_100kb_Windows
    11_Window_5mC_vs_5hmC_Scatter
    12_Window_Methylation_Heatmap
    13_Extreme_Window_Regional_Profiles
    13A_Regional_Profiles_Ranks_3_to_4
    13B_Regional_Profiles_Ranks_5_to_6
    13C_Regional_Profiles_Ranks_7_to_8
    13D_Regional_Profiles_Ranks_9_to_10
    14_Feature_Methylation_5mC_Boxplot
    15_Feature_Methylation_5hmC_Boxplot
    16_Feature_Methylation_5hmC_Boxplot_Zoom
    17_CpG_Feature_Coverage_Violin
    18_CpG_Feature_Length_vs_Measured_CpGs
    19_Feature_Weighted_Methylation_Summary
)
TABLE_FILES=(
    global_plot_data.tsv
    cpg_coverage_distribution.tsv
    cpg_coverage_thresholds.tsv
    cpg_site_plot_sample.tsv.gz
    high_depth_cpg_sites.tsv
    chromosome_plot_data.tsv
    window_plot_data.tsv.gz
    extreme_methylation_windows.tsv
    extreme_window_profile_data.tsv.gz
    13A_Regional_Profiles_Ranks_3_to_4_plot_data.tsv.gz
    13B_Regional_Profiles_Ranks_5_to_6_plot_data.tsv.gz
    13C_Regional_Profiles_Ranks_7_to_8_plot_data.tsv.gz
    13D_Regional_Profiles_Ranks_9_to_10_plot_data.tsv.gz
    feature_annotation_counts.tsv
    feature_region_methylation.tsv.gz
    feature_methylation_summary.tsv
    cpg_feature_length_plot_data.tsv.gz
    top_feature_regions.tsv
    sixbase_figure_mapping.tsv
)

outputs_validate() {
    local stem
    local table_file
    [[ -s "${COMPLETION_MARKER}" ]] || return 1
    [[ "${COMPLETION_MARKER}" -nt "${STAGE03_MARKER}" ]] || return 1
    [[ "${COMPLETION_MARKER}" -nt "${STAGE04_MARKER}" ]] || return 1
    [[ "${COMPLETION_MARKER}" -nt "${R_SCRIPT}" ]] || return 1
    [[ -s "${OUTPUT_DIR}/methylation_exploration_report.txt" ]] || return 1
    [[ -s "${OUTPUT_DIR}/ont_methylation_qc_report.html" ]] || return 1
    [[ -s "${OUTPUT_DIR}/R_session_info.txt" ]] || return 1
    grep -q 'single-sample descriptive analysis' "${OUTPUT_DIR}/methylation_exploration_report.txt" || return 1
    grep -q '<title>ONT sequencing and methylation QC report</title>' \
        "${OUTPUT_DIR}/ont_methylation_qc_report.html" || return 1
    for stem in "${FIGURE_STEMS[@]}"; do
        [[ -s "${OUTPUT_DIR}/figures/${stem}.pdf" ]] || return 1
        [[ "$(head -c 4 "${OUTPUT_DIR}/figures/${stem}.pdf")" == '%PDF' ]] || return 1
    done
    for table_file in "${TABLE_FILES[@]}"; do
        [[ -s "${OUTPUT_DIR}/tables/${table_file}" ]] || return 1
    done
    gzip -t "${OUTPUT_DIR}/tables/cpg_site_plot_sample.tsv.gz" || return 1
    gzip -t "${OUTPUT_DIR}/tables/window_plot_data.tsv.gz" || return 1
    gzip -t "${OUTPUT_DIR}/tables/feature_region_methylation.tsv.gz" || return 1
    awk -F '\t' 'NR > 1 && $2 == "Not scientifically available" {n++} END {exit !(n == 6)}' \
        "${OUTPUT_DIR}/tables/sixbase_figure_mapping.tsv" || return 1
}

#---------------------------------------------->



#---------------------------------------------->
# Skips a current validated result, or archives it before a forced rebuild without deleting prior analysis.

if outputs_validate && [[ ${FORCE} -ne 1 ]]; then
    log "Existing Stage 05 result is complete, current, and validated; skipping rebuild"
    exit 0
fi
if outputs_validate && [[ ${FORCE} -eq 1 ]]; then
    log "A valid Stage 05 result exists, but --force requested a recoverable rebuild"
fi

if [[ -n "$(find "${OUTPUT_DIR}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    if [[ ${FORCE} -ne 1 ]]; then
        log "ERROR: Stage 05 contains incomplete, stale, or unvalidated outputs"
        log "Inspect ${OUTPUT_DIR}; rerun with --force to archive it and rebuild"
        exit 1
    fi
    STALE_ROOT="${TMP_DIR}/stale_stage05_outputs"
    STALE_DIR="${STALE_ROOT}/${PIPELINE_RUN_STAMP}"
    mkdir -p "${STALE_ROOT}"
    if [[ -e "${STALE_DIR}" ]]; then
        log "ERROR: stale-output archive target already exists: ${STALE_DIR}"
        exit 1
    fi
    mv "${OUTPUT_DIR}" "${STALE_DIR}"
    mkdir -p "${OUTPUT_DIR}"
    log "Archived prior Stage 05 directory without deleting it: ${STALE_DIR}"
fi

#---------------------------------------------->



#---------------------------------------------->
# Builds a stable work fingerprint so expensive feature aggregation can resume after an interrupted plotting run.

work_fingerprint="$({
    stat -c '%n:%s:%Y' \
        "${STAGE03_MARKER}" "${STAGE04_MARKER}" "${CPG_BEDMETHYL}" \
        "${WINDOW_METHYLATION}" "${WINDOW_COVERAGE}" \
        "${METHYLATION_GENCODE_GFF3}" "${METHYLATION_CPG_ISLANDS_JSON}" \
        "${METHYLATION_CCRE_TABLE}" "${METHYLATION_INTERGENIC_BED}"
    sha256sum "$0" "${R_SCRIPT}" "${PROJECT_ROOT}/config/project.env"
} | sha256sum | awk '{print $1}')"
if [[ ${FORCE} -eq 1 ]]; then
    work_fingerprint="${work_fingerprint}.force.${PIPELINE_RUN_STAMP}"
fi
WORK_DIR="${TMP_DIR}/stage05_work/${work_fingerprint}"
mkdir -p "${WORK_DIR}"
log "Restartable working directory: ${WORK_DIR}"

#---------------------------------------------->



#---------------------------------------------->
# Runs the R analysis, which writes exact plot data before generating matching PDF and PNG figures.

log "Running Stage 05 R methylation exploration and SixBase-compatible visualization"
TMPDIR="${WORK_DIR}" "${RSCRIPT_BIN}" "${R_SCRIPT}" \
    --sample-id "${SAMPLE_ID}" \
    --reference-name "${REFERENCE_NAME}" \
    --fai "${REFERENCE_FAI}" \
    --global-levels "${GLOBAL_LEVELS}" \
    --chromosome-methylation "${CHROMOSOME_METHYLATION}" \
    --call-coverage "${CALL_COVERAGE}" \
    --window-methylation "${WINDOW_METHYLATION}" \
    --chromosome-coverage "${CHROMOSOME_COVERAGE}" \
    --window-coverage "${WINDOW_COVERAGE}" \
    --alignment-stats "${ALIGNMENT_STATS}" \
    --alignment-flagstat "${ALIGNMENT_FLAGSTAT}" \
    --coverage-summary "${COVERAGE_SUMMARY}" \
    --stage02-summary "${STAGE02_SUMMARY}" \
    --instrument-report "${ONT_RUN_REPORT}" \
    --bedmethyl "${CPG_BEDMETHYL}" \
    --gff3 "${METHYLATION_GENCODE_GFF3}" \
    --cpg-islands "${METHYLATION_CPG_ISLANDS_JSON}" \
    --ccre "${METHYLATION_CCRE_TABLE}" \
    --intergenic "${METHYLATION_INTERGENIC_BED}" \
    --output-dir "${WORK_DIR}" \
    --min-valid-coverage "${EXPLORATION_MIN_VALID_COVERAGE}" \
    --min-feature-cpgs "${EXPLORATION_MIN_FEATURE_CPGS}" \
    --top-window-count "${EXPLORATION_TOP_WINDOW_COUNT}" \
    --plot-dpi "${EXPLORATION_PLOT_DPI}"

#---------------------------------------------->



#---------------------------------------------->
# Verifies the complete work product before promoting any Stage 05 result.

for stem in "${FIGURE_STEMS[@]}"; do
    if [[ ! -s "${WORK_DIR}/figures/${stem}.pdf" ]]; then
        log "ERROR: missing Stage 05 figure: ${stem}.pdf"
        exit 1
    fi
done
for table_file in "${TABLE_FILES[@]}"; do
    if [[ ! -s "${WORK_DIR}/tables/${table_file}" ]]; then
        log "ERROR: missing Stage 05 plot-data table: ${table_file}"
        exit 1
    fi
done
for work_file in methylation_exploration_report.txt ont_methylation_qc_report.html R_session_info.txt; do
    if [[ ! -s "${WORK_DIR}/${work_file}" ]]; then
        log "ERROR: missing Stage 05 work product: ${work_file}"
        exit 1
    fi
done

#---------------------------------------------->



#---------------------------------------------->
# Promotes the validated figures, plot data, report, and R session record to the final Stage 05 directory.

log "Promoting validated Stage 05 figures and plot-data tables"
mv "${WORK_DIR}/figures" "${OUTPUT_DIR}/figures"
mv "${WORK_DIR}/tables" "${OUTPUT_DIR}/tables"
mv "${WORK_DIR}/methylation_exploration_report.txt" "${OUTPUT_DIR}/methylation_exploration_report.txt"
mv "${WORK_DIR}/ont_methylation_qc_report.html" "${OUTPUT_DIR}/ont_methylation_qc_report.html"
mv "${WORK_DIR}/R_session_info.txt" "${OUTPUT_DIR}/R_session_info.txt"

#---------------------------------------------->



#---------------------------------------------->
# Writes a reproducibility manifest covering sources, thresholds, code, annotation checksums, and output sizes.

MANIFEST_FILE="${MANIFEST_DIR}/${SAMPLE_ID}.stage05.${PIPELINE_RUN_STAMP}.manifest.tsv"
{
    printf 'record_type\tname\tvalue\n'
    printf 'metadata\tsample_id\t%s\n' "${SAMPLE_ID}"
    printf 'metadata\treference\t%s\n' "${REFERENCE_NAME}"
    printf 'metadata\tanalysis_type\tsingle_sample_methylation_exploration\n'
    printf 'metadata\tmin_valid_coverage\t%s\n' "${EXPLORATION_MIN_VALID_COVERAGE}"
    printf 'metadata\tmin_feature_cpgs\t%s\n' "${EXPLORATION_MIN_FEATURE_CPGS}"
    printf 'metadata\ttop_window_count\t%s\n' "${EXPLORATION_TOP_WINDOW_COUNT}"
    printf 'metadata\tfigure_formats\tPDF\n'
    printf 'metadata\tplot_dpi\t%s\n' "${EXPLORATION_PLOT_DPI}"
    printf 'input\tstage03_marker\t%s\n' "${STAGE03_MARKER}"
    printf 'input\tstage04_marker\t%s\n' "${STAGE04_MARKER}"
    printf 'input\tont_instrument_report\t%s\n' "${ONT_RUN_REPORT}"
    printf 'input\tcpg_bedmethyl\t%s\n' "${CPG_BEDMETHYL}"
    printf 'annotation\tgencode_gff3\t%s\n' "${METHYLATION_GENCODE_GFF3}"
    printf 'annotation\tcpg_islands\t%s\n' "${METHYLATION_CPG_ISLANDS_JSON}"
    printf 'annotation\tccre_table\t%s\n' "${METHYLATION_CCRE_TABLE}"
    printf 'annotation\tintergenic_bed\t%s\n' "${METHYLATION_INTERGENIC_BED}"
    printf 'tool\tR\t%s\n' "$("${RSCRIPT_BIN}" --version 2>&1 | head -n 1)"
    printf 'checksum\tstage05_shell_sha256\t%s\n' "$(sha256sum "$0" | awk '{print $1}')"
    printf 'checksum\tstage05_r_sha256\t%s\n' "$(sha256sum "${R_SCRIPT}" | awk '{print $1}')"
    printf 'checksum\tproject_config_sha256\t%s\n' "$(sha256sum "${PROJECT_ROOT}/config/project.env" | awk '{print $1}')"
    printf 'checksum\tgencode_gff3_sha256\t%s\n' "$(sha256sum "${METHYLATION_GENCODE_GFF3}" | awk '{print $1}')"
    printf 'checksum\tcpg_islands_sha256\t%s\n' "$(sha256sum "${METHYLATION_CPG_ISLANDS_JSON}" | awk '{print $1}')"
    printf 'checksum\tccre_table_sha256\t%s\n' "$(sha256sum "${METHYLATION_CCRE_TABLE}" | awk '{print $1}')"
    printf 'checksum\tintergenic_bed_sha256\t%s\n' "$(sha256sum "${METHYLATION_INTERGENIC_BED}" | awk '{print $1}')"
    printf 'output_bytes\tfigures\t%s\n' "$(du -sb "${OUTPUT_DIR}/figures" | awk '{print $1}')"
    printf 'output_bytes\ttables\t%s\n' "$(du -sb "${OUTPUT_DIR}/tables" | awk '{print $1}')"
    printf 'output_bytes\tmethylation_exploration_report.txt\t%s\n' "$(stat -c '%s' "${OUTPUT_DIR}/methylation_exploration_report.txt")"
    printf 'output_bytes\tont_methylation_qc_report.html\t%s\n' "$(stat -c '%s' "${OUTPUT_DIR}/ont_methylation_qc_report.html")"
    printf 'output_bytes\tR_session_info.txt\t%s\n' "$(stat -c '%s' "${OUTPUT_DIR}/R_session_info.txt")"
} > "${MANIFEST_FILE}.partial.${PIPELINE_RUN_STAMP}"
mv "${MANIFEST_FILE}.partial.${PIPELINE_RUN_STAMP}" "${MANIFEST_FILE}"
ln -sfn "$(basename -- "${MANIFEST_FILE}")" "${MANIFEST_DIR}/${SAMPLE_ID}.stage05.latest.manifest.tsv"

#---------------------------------------------->



#---------------------------------------------->
# Writes the completion marker last, revalidates the final contract, and removes only the expected work directory.

{
    printf 'status\tcomplete\n'
    printf 'sample_id\t%s\n' "${SAMPLE_ID}"
    printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
    printf 'stage03_marker\t%s\n' "${STAGE03_MARKER}"
    printf 'stage04_marker\t%s\n' "${STAGE04_MARKER}"
    printf 'manifest\t%s\n' "${MANIFEST_FILE}"
} > "${COMPLETION_MARKER}.partial.${PIPELINE_RUN_STAMP}"
mv "${COMPLETION_MARKER}.partial.${PIPELINE_RUN_STAMP}" "${COMPLETION_MARKER}"

if ! outputs_validate; then
    log "ERROR: final Stage 05 validation failed after promotion"
    exit 1
fi

case "${WORK_DIR}" in
    "${TMP_DIR}"/stage05_work/*)
        rm -rf -- "${WORK_DIR}"
        ;;
    *)
        log "ERROR: refusing to clean unexpected Stage 05 work directory: ${WORK_DIR}"
        exit 1
        ;;
esac

log "Stage 05 complete"
log "Exploration report: ${OUTPUT_DIR}/methylation_exploration_report.txt"
log "HTML QC report: ${OUTPUT_DIR}/ont_methylation_qc_report.html"
log "Figures: ${OUTPUT_DIR}/figures"
log "Plot-data tables: ${OUTPUT_DIR}/tables"
log "SixBase figure mapping: ${OUTPUT_DIR}/tables/sixbase_figure_mapping.tsv"
log "Manifest: ${MANIFEST_FILE}"

#---------------------------------------------->
