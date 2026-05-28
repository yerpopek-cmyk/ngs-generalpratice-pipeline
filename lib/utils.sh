#!/usr/bin/env bash
# =============================================================================
#  lib/utils.sh — Pipeline Utilities: logging, colors, checkpoints
# =============================================================================
#  Sourced via: source "$(dirname "${BASH_SOURCE[0]}")/lib/utils.sh"
# =============================================================================

# --- COLORS AND FORMATTING --------------------------------------------------
# Use tput for cross-terminal compatibility
if [[ -t 1 ]] && command -v tput &>/dev/null; then
    CLR_RESET=$(tput sgr0)
    CLR_BOLD=$(tput bold)
    CLR_RED=$(tput setaf 1)
    CLR_GREEN=$(tput setaf 2)
    CLR_YELLOW=$(tput setaf 3)
    CLR_BLUE=$(tput setaf 4)
    CLR_MAGENTA=$(tput setaf 5)
    CLR_CYAN=$(tput setaf 6)
    CLR_WHITE=$(tput setaf 7)
else
    # Fallback: no colors (e.g., when redirected to file)
    CLR_RESET="" CLR_BOLD="" CLR_RED="" CLR_GREEN=""
    CLR_YELLOW="" CLR_BLUE="" CLR_MAGENTA="" CLR_CYAN="" CLR_WHITE=""
fi
export CLR_RESET CLR_BOLD CLR_RED CLR_GREEN CLR_YELLOW CLR_BLUE CLR_MAGENTA CLR_CYAN CLR_WHITE

# --- HELPER FUNCTION: TIMESTAMP ---------------------------------
_timestamp() {
    date "+%Y-%m-%d %H:%M:%S"
}

# --- HELPER FUNCTION: DURATION ------------------------------------
# Usage: start=$(_time_start); ...; _time_elapsed "$start"
_time_start() { date +%s; }
_time_elapsed() {
    local start="$1"
    local end; end=$(date +%s)
    local delta=$(( end - start ))
    local mins=$(( delta / 60 ))
    local secs=$(( delta % 60 ))
    printf "%dm %02ds" "$mins" "$secs"
}

# --- MAIN LOGGING FUNCTIONS --------------------------------------------

# log_step: Header for a major pipeline stage
# Usage: log_step "01" "Data Preparation and Alignment"
log_step() {
    local num="$1"
    local title="$2"
    echo ""
    echo "${CLR_BOLD}${CLR_BLUE}╔══════════════════════════════════════════════════════════════╗${CLR_RESET}"
    printf "${CLR_BOLD}${CLR_BLUE}║  STEP %-2s: %-52s║${CLR_RESET}\n" "$num" "$title"
    echo "${CLR_BOLD}${CLR_BLUE}╚══════════════════════════════════════════════════════════════╝${CLR_RESET}"
    echo ""
    _log_to_file "=== STEP ${num}: ${title} ==="
}

# log_substep: Sub-step inside a major stage
log_substep() {
    local title="$1"
    echo ""
    echo "${CLR_BOLD}${CLR_CYAN}  ▶ ${title}${CLR_RESET}"
    _log_to_file "  >> ${title}"
}

# log_info: Standard informational message
log_info() {
    local msg="$1"
    echo "${CLR_WHITE}[$(_timestamp)] ${CLR_RESET}${msg}"
    _log_to_file "[$(_timestamp)] INFO: ${msg}"
}

# log_theory: EDUCATIONAL BLOCK — explains biology/algorithm
# Highlighted with colors so the user sees theory in the console
log_theory() {
    local title="$1"
    shift
    echo ""
    echo "${CLR_MAGENTA}${CLR_BOLD}  ┌─ 📚 THEORY: ${title}${CLR_RESET}"
    while [[ $# -gt 0 ]]; do
        echo "${CLR_MAGENTA}  │  $1${CLR_RESET}"
        shift
    done
    echo "${CLR_MAGENTA}  └──────────────────────────────────────────${CLR_RESET}"
    echo ""
    _log_to_file "  [THEORY] ${title}"
}

# log_cmd: Display command before running (educational)
log_cmd() {
    local cmd="$1"
    echo "${CLR_YELLOW}  ┌─ CMD:${CLR_RESET}"
    # Multi-line wrap for readability
    echo "$cmd" | fold -s -w 75 | sed 's/^/  │  /'
    echo "${CLR_YELLOW}  └──────────${CLR_RESET}"
    _log_to_file "  [CMD] ${cmd}"
}

# log_ok: Step completed successfully
log_ok() {
    local msg="$1"
    local elapsed="${2:-}"
    if [[ -n "$elapsed" ]]; then
        echo "${CLR_GREEN}${CLR_BOLD}  ✔ OK: ${msg} ${CLR_RESET}${CLR_WHITE}(${elapsed})${CLR_RESET}"
    else
        echo "${CLR_GREEN}${CLR_BOLD}  ✔ OK: ${msg}${CLR_RESET}"
    fi
    _log_to_file "  [OK] ${msg} ${elapsed}"
}

# log_warn: Warning (non-fatal)
log_warn() {
    local msg="$1"
    echo "${CLR_YELLOW}${CLR_BOLD}  ⚠  WARN: ${msg}${CLR_RESET}" >&2
    _log_to_file "  [WARN] ${msg}"
}

# log_error: Error (fatal)
log_error() {
    local msg="$1"
    echo "${CLR_RED}${CLR_BOLD}  ✘ ERROR: ${msg}${CLR_RESET}" >&2
    _log_to_file "  [ERROR] ${msg}"
}



# log_hardware_tip: Info about hardware acceleration (Parabricks / DRAGEN)
log_hardware_tip() {
    echo ""
    echo "${CLR_MAGENTA}  ┌─ ⚡ HARDWARE ACCELERATION (Info):${CLR_RESET}"
    echo "${CLR_MAGENTA}  │${CLR_RESET}"
    echo "${CLR_MAGENTA}  │  NVIDIA Parabricks (GPU):${CLR_RESET}"
    echo "${CLR_MAGENTA}  │    • Migrates BWA + MarkDup + BQSR + HaplotypeCaller to GPU${CLR_RESET}"
    echo "${CLR_MAGENTA}  │    • Speedup: 30–60x compared to CPU${CLR_RESET}"
    echo "${CLR_MAGENTA}  │    • WGS 30x: ~40 mins instead of 30+ hours${CLR_RESET}"
    echo "${CLR_MAGENTA}  │    • Command: pbrun fq2bam / pbrun haplotypecaller${CLR_RESET}"
    echo "${CLR_MAGENTA}  │${CLR_RESET}"
    echo "${CLR_MAGENTA}  │  Illumina DRAGEN (FPGA):${CLR_RESET}"
    echo "${CLR_MAGENTA}  │    • BCL→FASTQ + alignment + variant calling = hardware accelerated${CLR_RESET}"
    echo "${CLR_MAGENTA}  │    • WGS 30x: ~25 mins directly on sequencer${CLR_RESET}"
    echo "${CLR_MAGENTA}  │    • FPGA = hardwired algorithms, no CPU/GPU overhead${CLR_RESET}"
    echo "${CLR_MAGENTA}  └──────────────────────────────────────────${CLR_RESET}"
    echo ""
}

# --- INTERNAL FUNCTION: LOG TO FILE ----------------------------------
_log_to_file() {
    if [[ -n "${LOG_FILE:-}" ]]; then
        echo "$1" >> "${LOG_FILE}"
    fi
}

# --- LOG INITIALIZATION ------------------------------------------------------
# Called at the start of each step script
init_log() {
    local step_name="$1"
    export LOG_FILE="${LOG_DIR:-/tmp}/pipeline_${step_name}.log"
    mkdir -p "$(dirname "$LOG_FILE")"
    echo "=== Pipeline Log: ${step_name} ===" > "$LOG_FILE"
    echo "Started: $(_timestamp)" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"
}

# --- CHECKPOINT SYSTEM ------------------------------------------------------
# Resumes the pipeline without repeating completed steps
# Status is stored in the .checkpoints/ directory

# Check: Was the step already completed successfully?
is_done() {
    local step="$1"
    local checkpoint="${CHECKPOINT_DIR:-${BASE_DIR}/.checkpoints}/${step}.done"
    [[ -f "$checkpoint" ]]
}

# Mark a step as completed
mark_done() {
    local step="$1"
    local checkpoint_dir="${CHECKPOINT_DIR:-${BASE_DIR}/.checkpoints}"
    mkdir -p "$checkpoint_dir"
    echo "$(_timestamp)" > "${checkpoint_dir}/${step}.done"
    log_ok "Checkpoint saved: ${step}"
}

# Reset a specific checkpoint (to re-run it)
reset_checkpoint() {
    local step="$1"
    local checkpoint="${CHECKPOINT_DIR:-${BASE_DIR}/.checkpoints}/${step}.done"
    [[ -f "$checkpoint" ]] && rm "$checkpoint"
    log_info "Checkpoint reset: ${step}"
}

# Reset ALL checkpoints (start over)
reset_all_checkpoints() {
    rm -rf "${CHECKPOINT_DIR:-${BASE_DIR}/.checkpoints}"
    log_warn "All checkpoints reset. Pipeline will restart."
}

# --- FILE AND TOOL CHECKS ------------------------------------------

# Ensure file exists
require_file() {
    local f="$1"
    local desc="${2:-file}"
    if [[ ! -f "$f" ]]; then
        log_error "${desc} not found: ${f}"
        return 1
    fi
}

# Ensure directory exists; if not — create it
ensure_dir() {
    local d="$1"
    mkdir -p "$d"
}

# Ensure tool is available in PATH
require_tool() {
    local tool="$1"
    local hint="${2:-}"
    if ! command -v "$tool" &>/dev/null; then
        log_error "Tool not found: '${tool}'. ${hint}"
        return 1
    fi
    return 0
}

# --- HELPER FUNCTIONS -------------------------------------------------

# Run command with logging: in demo — display only, in real — execute
run_cmd() {
    local cmd="$1"
    local desc="${2:-}"
    local start; start=$(_time_start)

    log_cmd "$cmd"

    if ! eval "$cmd"; then
        log_error "Command failed with error: ${cmd}"
        return 1
    fi

    local elapsed; elapsed=$(_time_elapsed "$start")
    [[ -n "$desc" ]] && log_ok "$desc" "$elapsed"
    return 0
}

# Run GNU Parallel with logging
run_parallel() {
    local jobs="$1"
    local cmd_template="$2"
    local inputs="$3"
    local desc="${4:-parallel execution}"

    local joblog="${LOG_DIR:-/tmp}/parallel_$(date +%s).log"

    log_info "GNU Parallel: -j ${jobs} | Inputs: $(echo "$inputs" | wc -w | tr -d ' ') files"
    log_cmd "parallel --progress --joblog ${joblog} -j ${jobs} '${cmd_template}' ::: ${inputs}"

    # shellcheck disable=SC2086
    parallel --progress --joblog "$joblog" -j "$jobs" "$cmd_template" ::: $inputs

    log_ok "$desc"
}

# Count files in a directory by pattern
count_files() {
    local pattern="$1"
    # shellcheck disable=SC2012
    ls $pattern 2>/dev/null | wc -l | tr -d ' '
}

# Final summary at the end of the pipeline
print_summary() {
    echo ""
    echo "${CLR_BOLD}${CLR_GREEN}╔══════════════════════════════════════════════════════════════╗${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║           🎉  PIPELINE COMPLETED SUCCESSFULLY               ║${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}╠══════════════════════════════════════════════════════════════╣${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║  Results:${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║    VCF:          ${OUT_DIR}/chrY_MT_final.vcf.gz${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║    Annotation:   ${OUT_ANNOTATION}/chrY_MT_final.hg38_multianno.txt${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║    Haplogroups:  ${OUT_HAPLOGROUPS}/haplogroups.txt${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║    Kinship:      ${OUT_KINSHIP}/kinship_summary.tsv${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║    QC Report:    ${OUT_MULTIQC}/multiqc_report.html${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}╚══════════════════════════════════════════════════════════════╝${CLR_RESET}"
    echo ""
}
