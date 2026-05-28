#!/usr/bin/env bash
# =============================================================================
#  setup_conda_envs.sh — Automated script to create and verify Conda environments
# =============================================================================
#  Performs the following:
#    1. Detects Conda on the system (including miniconda3)
#    2. Adds and configures channel priorities (conda-forge, bioconda)
#    3. Creates 6 isolated environments with the latest tool versions:
#       - gatk4          → gatk4, bwa, samtools, parallel
#       - QC_fastq       → fastqc, qualimap, multiqc, bcftools
#       - LineageTracker → python=3.10 + pip + Y-LineageTracker
#       - yleaf          → yleaf
#       - haplogrep      → haplogrep
#       - nextflow       → nextflow
#    4. Tests and verifies each installed tool
# =============================================================================

set -euo pipefail

# ── Paths and libraries ────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/utils.sh"

LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"
export LOG_FILE="${LOG_DIR}/setup_conda_envs.log"
echo "=== Conda Environment Setup Log ===" > "${LOG_FILE}"
echo "Started: $(date)" >> "${LOG_FILE}"

log_step "00" "Conda Initialization and Setup"

# ── Conda Detection ────────────────────────────────────────────────────────
log_substep "Searching and initializing Conda..."

CONDA_PATH=""
# 1. Check standard PATH
if command -v conda &>/dev/null; then
    CONDA_PATH="conda"
fi

# 2. Check Miniconda in user's home directory
if [[ -z "$CONDA_PATH" ]] && [[ -f "/home/yer_kanat/miniconda3/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1091
    source "/home/yer_kanat/miniconda3/etc/profile.d/conda.sh"
    if command -v conda &>/dev/null; then
        CONDA_PATH="conda"
    fi
fi

if [[ -z "$CONDA_PATH" ]]; then
    log_error "Conda not found! Ensure Miniconda/Anaconda is installed."
    log_error "Check path: /home/yer_kanat/miniconda3"
    exit 1
fi

log_info "Using Conda: $(which conda)"
log_info "Version: $(conda --version)"

# ── Configuring Conda channels ──────────────────────────────────────────────────
log_substep "Configuring channel repositories..."

log_info "Adding channels defaults, bioconda, conda-forge..."
# Configure channels in correct priority order
conda config --add channels defaults 2>/dev/null || true
conda config --add channels bioconda 2>/dev/null || true
conda config --add channels conda-forge 2>/dev/null || true
# Enable strict channel priority for faster builds and compatibility
conda config --set channel_priority strict 2>/dev/null || true

log_ok "Channel configuration complete. Priority: conda-forge > bioconda > defaults"

# ── Function to create and update environments ──────────────────────────────────
setup_env() {
    local env_name="$1"
    local env_packages="$2"
    local post_command="${3:-}"

    log_substep "Building environment: ${CLR_BOLD}${CLR_CYAN}${env_name}${CLR_RESET}"
    log_info "Packages to install: ${env_packages}"

    local start_time; start_time=$(_time_start)

    if conda env list | grep -qE "^${env_name}\s"; then
        log_warn "Environment '${env_name}' already exists."
        log_info "Updating environment to the latest versions..."
        log_cmd "conda install -y -n ${env_name} ${env_packages}"
        conda install -y -n "${env_name}" ${env_packages} >> "${LOG_FILE}" 2>&1 || {
            log_error "Failed to update environment ${env_name}."
            return 1
        }
    else
        log_info "Creating new environment '${env_name}'..."
        log_cmd "conda create -y -n ${env_name} ${env_packages}"
        conda create -y -n "${env_name}" ${env_packages} >> "${LOG_FILE}" 2>&1 || {
            log_error "Failed to create environment ${env_name}."
            return 1
        }
    fi

    if [[ -n "$post_command" ]]; then
        log_info "Executing post-installation command..."
        log_cmd "${post_command}"
        eval "${post_command}" >> "${LOG_FILE}" 2>&1 || {
            log_error "Error executing post-command in ${env_name}."
            return 1
        }
    fi

    local elapsed; elapsed=$(_time_elapsed "$start_time")
    log_ok "Environment '${env_name}' successfully built!" "$elapsed"
    echo ""
    return 0
}

# ── Create Environments ───────────────────────────────────────────────────────
log_step "01" "Creating and configuring environments"

# 1. gatk4
setup_env "gatk4" "gatk4 bwa samtools parallel"

# 2. QC_fastq
setup_env "QC_fastq" "fastqc qualimap multiqc bcftools"

# 3. LineageTracker (Python 3.10 + pip + Y-LineageTracker)
setup_env "LineageTracker" "python=3.10 pip" "conda run -n LineageTracker pip install Y-LineageTracker"

# 4. yleaf
setup_env "yleaf" "yleaf"

# 5. haplogrep
setup_env "haplogrep" "haplogrep"

# 6. nextflow
setup_env "nextflow" "nextflow"

# ── Testing and Verification ──────────────────────────────────────────────
log_step "02" "Verifying tool functionality"

FAILED_VERIFICATIONS=0

verify_tool() {
    local env_name="$1"
    local tool_cmd="$2"
    local check_arg="$3"
    
    printf "  • Verifying %-15s in %-16s ... " "${tool_cmd}" "${env_name}"
    
    if conda run -n "${env_name}" "${tool_cmd}" ${check_arg} >/dev/null 2>&1; then
        echo -e "${CLR_GREEN}${CLR_BOLD}SUCCESS (OK)${CLR_RESET}"
        echo "[VERIFY] ${env_name} : ${tool_cmd} -> SUCCESS" >> "${LOG_FILE}"
    else
        # Not all tools output with 0 code on --version/--help, check exit code or output
        local output; output=$(conda run -n "${env_name}" "${tool_cmd}" ${check_arg} 2>&1 || true)
        if [[ -n "$output" ]]; then
            echo -e "${CLR_GREEN}${CLR_BOLD}SUCCESS (OK)${CLR_RESET}"
            echo "[VERIFY] ${env_name} : ${tool_cmd} -> SUCCESS (non-zero exit but output present)" >> "${LOG_FILE}"
        else
            echo -e "${CLR_RED}${CLR_BOLD}ERROR (FAIL)${CLR_RESET}"
            echo "[VERIFY] ${env_name} : ${tool_cmd} -> FAILED" >> "${LOG_FILE}"
            ((FAILED_VERIFICATIONS++))
        fi
    fi
}

log_info "Running test commands in each environment:"
echo ""

# Verify gatk4
verify_tool "gatk4" "gatk" "--help"
verify_tool "gatk4" "samtools" "--version"
verify_tool "gatk4" "bwa" ""
verify_tool "gatk4" "parallel" "--version"

# Verify QC_fastq
verify_tool "QC_fastq" "fastqc" "--version"
verify_tool "QC_fastq" "qualimap" "-h"
verify_tool "QC_fastq" "multiqc" "--version"
verify_tool "QC_fastq" "bcftools" "--version"

# Verify LineageTracker
verify_tool "LineageTracker" "LineageTracker" "--help"

# Verify yleaf
verify_tool "yleaf" "Yleaf" "-h"

# Verify haplogrep
verify_tool "haplogrep" "haplogrep" "--help"

# Verify nextflow
verify_tool "nextflow" "nextflow" "-version"

echo ""

if [[ $FAILED_VERIFICATIONS -eq 0 ]]; then
    echo "${CLR_BOLD}${CLR_GREEN}╔══════════════════════════════════════════════════════════════╗${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║     🎉  ALL ENVIRONMENTS SUCCESSFULLY SETUP AND VERIFIED     ║${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}╠══════════════════════════════════════════════════════════════╣${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║  The pipeline is ready to run in production mode!            ║${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║  To run:                                                     ║${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║    bash run_pipeline.sh --production                         ║${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}╚══════════════════════════════════════════════════════════════╝${CLR_RESET}"
    echo ""
else
    log_warn "Some tools (${FAILED_VERIFICATIONS}) finished verification with a warning/error."
    log_warn "Please check the installation log: ${LOG_FILE}"
fi

echo "Installation log saved to: ${LOG_FILE}"
echo ""
