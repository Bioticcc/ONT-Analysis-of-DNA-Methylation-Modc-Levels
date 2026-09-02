#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../config/project.env
source "${SCRIPT_DIR}/../config/project.env"

READ_LINK="${INPUT_ROOT}/reads/${SAMPLE_ID}"

mkdir -p \
    "${INPUT_ROOT}/reads" \
    "${REFERENCE_DIR}" \
    "${SAMPLE_OUTPUT}/00_manifests" \
    "${SAMPLE_OUTPUT}/01_alignment" \
    "${SAMPLE_OUTPUT}/02_alignment_qc" \
    "${SAMPLE_OUTPUT}/03_ont_qc_coverage" \
    "${SAMPLE_OUTPUT}/04_methylation" \
    "${SAMPLE_OUTPUT}/05_methylation_exploration" \
    "${SAMPLE_OUTPUT}/logs" \
    "${SAMPLE_OUTPUT}/tmp"

if [[ -e "${READ_LINK}" && ! -L "${READ_LINK}" ]]; then
    printf 'Refusing to replace non-symlink path: %s\n' "${READ_LINK}" >&2
    exit 1
fi
ln -sfn "${RUN_ROOT}" "${READ_LINK}"

printf 'Created data layout under %s\n' "${DATA_ROOT}"
printf 'Source run remains at %s\n' "${RUN_ROOT}"
printf 'Readable input link: %s -> %s\n' "${READ_LINK}" "${RUN_ROOT}"
