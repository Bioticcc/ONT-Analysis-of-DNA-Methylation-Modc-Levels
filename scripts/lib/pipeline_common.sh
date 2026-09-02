#!/usr/bin/env bash

PIPELINE_RUN_STAMP="$(date '+%Y%m%d_%H%M%S')_${BASHPID}"
PIPELINE_LOG_FILE=""
PIPELINE_HEARTBEAT_PID=""

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*"
}

pipeline_exit() {
    local exit_code=$?
    if [[ -n "${PIPELINE_HEARTBEAT_PID}" ]]; then
        kill "${PIPELINE_HEARTBEAT_PID}" 2>/dev/null || true
        wait "${PIPELINE_HEARTBEAT_PID}" 2>/dev/null || true
    fi
    if [[ ${exit_code} -eq 0 ]]; then
        log "EXIT: success"
    else
        log "EXIT: failure (code ${exit_code})"
    fi
}

start_pipeline_heartbeat() {
    local interval_seconds="$1"
    local data_path="$2"
    (
        # The heartbeat and its sleep child must not inherit the stage lock.
        # Otherwise the lock can remain held briefly after the main script exits.
        exec 9>&-
        while sleep "${interval_seconds}"; do
            local available_kib
            local data_free_kib
            available_kib="$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)"
            data_free_kib="$(df -Pk "${data_path}" | awk 'NR == 2 {print $4}')"
            log "HEARTBEAT: load=$(cut -d ' ' -f 1-3 /proc/loadavg); mem_available_kib=${available_kib}; data_free_kib=${data_free_kib}"
        done
    ) &
    PIPELINE_HEARTBEAT_PID=$!
    log "Started ${interval_seconds}s resource heartbeat (PID ${PIPELINE_HEARTBEAT_PID})"
}

start_pipeline_log() {
    local log_dir="$1"
    local stage_name="$2"
    mkdir -p "${log_dir}"
    PIPELINE_LOG_FILE="${log_dir}/${stage_name}.${PIPELINE_RUN_STAMP}.log"
    ln -sfn "$(basename -- "${PIPELINE_LOG_FILE}")" "${log_dir}/${stage_name}.latest.log"
    exec > >(tee -a "${PIPELINE_LOG_FILE}") 2>&1
    trap pipeline_exit EXIT
    log "START: ${stage_name}"
    log "Log file: ${PIPELINE_LOG_FILE}"
    log "Host: $(hostname)"
    log "Script: $0"
}

acquire_pipeline_lock() {
    local lock_file="$1"
    mkdir -p "$(dirname -- "${lock_file}")"
    exec 9>"${lock_file}"
    if ! flock -n 9; then
        log "ERROR: another instance holds lock ${lock_file}"
        exit 1
    fi
    log "Acquired lock: ${lock_file}"
}

require_executable() {
    local executable="$1"
    if [[ "${executable}" == */* ]]; then
        if [[ ! -x "${executable}" ]]; then
            log "ERROR: executable not found: ${executable}"
            exit 1
        fi
    elif ! command -v "${executable}" >/dev/null 2>&1; then
        log "ERROR: executable not found on PATH: ${executable}"
        exit 1
    fi
}
