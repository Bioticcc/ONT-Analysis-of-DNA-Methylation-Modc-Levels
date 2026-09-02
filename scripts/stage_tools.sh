#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
ARCHIVE_DIR="${PROJECT_ROOT}/tools/archives"
BIN_DIR="${PROJECT_ROOT}/bin"

DORADO_VERSION="2.0.0"
MINIMAP2_VERSION="2.31"
MODKIT_VERSION="0.6.4"
MOSDEPTH_VERSION="0.3.14"

DORADO_ARCHIVE="dorado-${DORADO_VERSION}-linux-x64.tar.gz"
MINIMAP2_ARCHIVE="minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2"
MODKIT_ARCHIVE="modkit_v${MODKIT_VERSION}_u16_x86_64.tar.gz"
MOSDEPTH_ARCHIVE="mosdepth-v${MOSDEPTH_VERSION}-linux-x86_64"

DORADO_URL="https://cdn.oxfordnanoportal.com/software/analysis/${DORADO_ARCHIVE}"
MINIMAP2_URL="https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/${MINIMAP2_ARCHIVE}"
MODKIT_URL="https://github.com/nanoporetech/modkit/releases/download/v${MODKIT_VERSION}/${MODKIT_ARCHIVE}"
MOSDEPTH_URL="https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth"

MINIMAP2_SHA256="300bc287f05eb890c6211fa7db043ce98320a401621fadd1cfdbeabd1a6e4ab5"
MODKIT_SHA256="fb332c691431bd336eb0a81cbca17d2a35caf442ac48277ed3e296c2fe061d80"
DORADO_SHA256="311c4be8fe5177ee2ffe08e2b6ae08bd7f4f1c11caf273e6161562bd9ba48b49"
MOSDEPTH_SHA256="c5182b74a8f1b66710efa16e122cbc8a197834874b103e7c5c0bd9a6265ae7b6"

MODKIT_INSTALL_ROOT="${PROJECT_ROOT}/tools/modkit-v${MODKIT_VERSION}"
MODKIT_BINARY="${MODKIT_INSTALL_ROOT}/dist/modkit"
MODKIT_RELEASE_BINARY="${MODKIT_INSTALL_ROOT}/dist_modkit_v${MODKIT_VERSION}_cd85862/modkit"

mkdir -p "${ARCHIVE_DIR}" "${BIN_DIR}" "${PROJECT_ROOT}/tools"

download_if_missing() {
    local url="$1"
    local destination="$2"
    if [[ ! -s "${destination}" ]]; then
        curl -fL --retry 3 --retry-delay 5 --continue-at - \
            "${url}" -o "${destination}"
    fi
}

check_sha256() {
    local expected="$1"
    local archive="$2"
    printf '%s  %s\n' "${expected}" "${archive}" | sha256sum --check --status
}

download_if_missing "${DORADO_URL}" "${ARCHIVE_DIR}/${DORADO_ARCHIVE}"
download_if_missing "${MINIMAP2_URL}" "${ARCHIVE_DIR}/${MINIMAP2_ARCHIVE}"
download_if_missing "${MODKIT_URL}" "${ARCHIVE_DIR}/${MODKIT_ARCHIVE}"
download_if_missing "${MOSDEPTH_URL}" "${ARCHIVE_DIR}/${MOSDEPTH_ARCHIVE}"

check_sha256 "${DORADO_SHA256}" "${ARCHIVE_DIR}/${DORADO_ARCHIVE}"
check_sha256 "${MINIMAP2_SHA256}" "${ARCHIVE_DIR}/${MINIMAP2_ARCHIVE}"
check_sha256 "${MODKIT_SHA256}" "${ARCHIVE_DIR}/${MODKIT_ARCHIVE}"
check_sha256 "${MOSDEPTH_SHA256}" "${ARCHIVE_DIR}/${MOSDEPTH_ARCHIVE}"

if [[ ! -x "${PROJECT_ROOT}/tools/dorado-${DORADO_VERSION}-linux-x64/bin/dorado" ]]; then
    tar -xzf "${ARCHIVE_DIR}/${DORADO_ARCHIVE}" -C "${PROJECT_ROOT}/tools"
fi

if [[ ! -x "${PROJECT_ROOT}/tools/minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2" ]]; then
    tar -xjf "${ARCHIVE_DIR}/${MINIMAP2_ARCHIVE}" -C "${PROJECT_ROOT}/tools"
fi

if [[ -x "${MODKIT_RELEASE_BINARY}" ]]; then
    MODKIT_BINARY="${MODKIT_RELEASE_BINARY}"
elif [[ ! -x "${MODKIT_BINARY}" ]]; then
    mkdir -p "${MODKIT_INSTALL_ROOT}"
    tar -xzf "${ARCHIVE_DIR}/${MODKIT_ARCHIVE}" -C "${MODKIT_INSTALL_ROOT}"
fi

if [[ -x "${MODKIT_INSTALL_ROOT}/dist/modkit" ]]; then
    MODKIT_BINARY="${MODKIT_INSTALL_ROOT}/dist/modkit"
elif [[ -x "${MODKIT_RELEASE_BINARY}" ]]; then
    MODKIT_BINARY="${MODKIT_RELEASE_BINARY}"
else
    printf 'Unable to locate Modkit executable under %s\n' "${MODKIT_INSTALL_ROOT}" >&2
    exit 1
fi

ln -sfn "${PROJECT_ROOT}/tools/dorado-${DORADO_VERSION}-linux-x64/bin/dorado" \
    "${BIN_DIR}/dorado"
ln -sfn "${PROJECT_ROOT}/tools/minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2" \
    "${BIN_DIR}/minimap2"
ln -sfn "${MODKIT_BINARY}" "${BIN_DIR}/modkit"
chmod 0755 "${ARCHIVE_DIR}/${MOSDEPTH_ARCHIVE}"
ln -sfn "${ARCHIVE_DIR}/${MOSDEPTH_ARCHIVE}" "${BIN_DIR}/mosdepth"

"${BIN_DIR}/dorado" --version
"${BIN_DIR}/minimap2" --version
"${BIN_DIR}/modkit" --version
"${BIN_DIR}/mosdepth" --version

printf 'All staged tool artifacts passed their pinned SHA-256 checks.\n'
