#!/usr/bin/env bash
# =============================================================================
#  run_pipeline.sh — Main execution script for the NGS Y/MT Analysis pipeline
# =============================================================================
#  USAGE:
#    bash run_pipeline.sh                   # Run the full pipeline
#    bash run_pipeline.sh --step 02         # Run only one step
#    bash run_pipeline.sh --from 03         # Run from step 03 to the end
#    bash run_pipeline.sh --reset           # Reset checkpoints, start over
#    bash run_pipeline.sh --threads 8       # Override number of threads
# =============================================================================

set -euo pipefail

# ── Absolute path to the pipeline directory ───────────────────────────────────
PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BASE_DIR="${PIPELINE_ROOT}"

# ── Source config and utilities ──────────────────────────────────────────────
source "${PIPELINE_ROOT}/config.sh"
source "${PIPELINE_ROOT}/lib/utils.sh"

# ── Global log ────────────────────────────────────────────────────────────
export LOG_FILE="${LOG_DIR:-${PIPELINE_ROOT}/logs}/pipeline_main.log"
mkdir -p "$(dirname "$LOG_FILE")"
mkdir -p "${CHECKPOINT_DIR}"

# =============================================================================
#  ARGUMENT PARSING
# =============================================================================
ONLY_STEP=""
FROM_STEP=""
DO_RESET=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --step)
            ONLY_STEP="$2"
            shift 2 ;;
        --from)
            FROM_STEP="$2"
            shift 2 ;;
        --reset)
            DO_RESET=true
            shift ;;
        --threads)
            export THREADS="$2"
            export PARALLEL_JOBS=1
            export SPARK_CORES="$2"
            shift 2 ;;
        --outdir)
            export OUT_DIR="$2"
            shift 2 ;;
        -h|--help)
            echo "Usage: bash run_pipeline.sh [OPTIONS]"
            echo ""
            echo "  --step N         Run only step N (01–06)"
            echo "  --from N         Run from step N to the end"
            echo "  --reset          Reset all checkpoints"
            echo "  --threads N      Number of threads (default: ${THREADS})"
            echo "  --outdir PATH    Override output directory"
            echo "  -h, --help       Show this help"
            echo ""
            exit 0 ;;
        *)
            log_warn "Unknown argument: $1"
            shift ;;
    esac
done

# =============================================================================
#  RESET CHECKPOINTS
# =============================================================================
if [[ "$DO_RESET" == "true" ]]; then
    reset_all_checkpoints
fi

# =============================================================================
#  START HEADER
# =============================================================================
if [[ -t 1 ]] && [[ -n "${TERM:-}" ]] && command -v clear &>/dev/null; then
    clear || true
fi
echo ""
echo "${CLR_BOLD}${CLR_CYAN}╔══════════════════════════════════════════════════════════════════╗${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}║    🧬  NGS Y/MT ANALYSIS PIPELINE  v${PIPELINE_VERSION}                   ║${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}╠══════════════════════════════════════════════════════════════════╣${CLR_RESET}"
echo "${CLR_BOLD}${CLR_GREEN}║    MODE: PRODUCTION  (real data, real tools)                     ║${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}╠══════════════════════════════════════════════════════════════════╣${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}║  Pipeline stages:                                                ║${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}║   01. Data Preparation      FASTQ → uBAM → BWA alignment         ║${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}║   02. Variant Calling       MarkDup → BQSR → HaplotypeCaller     ║${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}║   03. Quality Control       FastQC → QualiMap → MultiQC          ║${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}║   04. Annotation            ANNOVAR → ClinVar → gnomAD           ║${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}║   05. Haplogroups           LineageTracker → Haplogrep           ║${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}║   06. Kinship               PLINK (IBD/PI_HAT) → KING            ║${CLR_RESET}"
echo "${CLR_BOLD}${CLR_CYAN}╠══════════════════════════════════════════════════════════════════╣${CLR_RESET}"
printf "${CLR_BOLD}${CLR_CYAN}║  %-64s║${CLR_RESET}\n" "Threads: ${THREADS} | Spark: ${SPARK_CORES} cores | Memory: ${SPARK_MEMORY}"
printf "${CLR_BOLD}${CLR_CYAN}║  %-64s║${CLR_RESET}\n" "Build: ${GENOME_BUILD} | Log: logs/pipeline_main.log"
echo "${CLR_BOLD}${CLR_CYAN}╚══════════════════════════════════════════════════════════════════╝${CLR_RESET}"
echo ""

PIPELINE_START=$(_time_start)

# =============================================================================
#  FUNCTION TO RUN A SINGLE STEP
# =============================================================================
run_step() {
    local step_num="$1"
    local step_script="${PIPELINE_ROOT}/steps/${step_num}_*.sh"

    # Find the script (glob)
    local found_script
    found_script=$(ls ${step_script} 2>/dev/null | head -1)

    if [[ -z "$found_script" ]]; then
        log_error "Script for step ${step_num} not found: ${step_script}"
        return 1
    fi

    local step_start; step_start=$(_time_start)
    local step_status=0
    bash "$found_script" || step_status=$?
    local step_elapsed; step_elapsed=$(_time_elapsed "$step_start")
    log_info "⏱  Step ${step_num} took: ${step_elapsed}"
    echo ""
    return "$step_status"
}

# =============================================================================
#  DETERMINE LIST OF STEPS TO RUN
# =============================================================================
ALL_STEPS=("01" "02" "03" "04" "05" "06")

if [[ -n "$ONLY_STEP" ]]; then
    STEPS_TO_RUN=("$ONLY_STEP")
elif [[ -n "$FROM_STEP" ]]; then
    STEPS_TO_RUN=()
    found_from=false
    for s in "${ALL_STEPS[@]}"; do
        if [[ "$s" == "$FROM_STEP" ]]; then found_from=true; fi
        if $found_from; then STEPS_TO_RUN+=("$s"); fi
    done
    if [[ ${#STEPS_TO_RUN[@]} -eq 0 ]]; then
        log_error "Step --from ${FROM_STEP} not found. Valid: ${ALL_STEPS[*]}"
        exit 1
    fi
else
    STEPS_TO_RUN=("${ALL_STEPS[@]}")
fi

log_info "Steps to run: ${STEPS_TO_RUN[*]}"
echo ""

# =============================================================================
#  MAIN EXECUTION LOOP
# =============================================================================
FAILED_STEPS=()

for step in "${STEPS_TO_RUN[@]}"; do
    if run_step "$step"; then
        :
    else
        log_error "Step ${step} failed with an error!"
        log_error "Aborting pipeline. Fix the error and run with --from ${step}"
        exit 1
    fi
done

# =============================================================================
#  FINAL SUMMARY
# =============================================================================
# Synchronize results with Windows FS on successful completion (if OUT_DIR is in internal FS)
if [[ "${OUT_DIR}" == "/home/yer_kanat/ngs-pipeline/out" ]]; then
    log_info "Copying all results and intermediate files to Windows working directory..."
    mkdir -p "${BASE_DIR}/out"
    cp -r "${OUT_DIR}"/* "${BASE_DIR}/out/" 2>/dev/null || true
    log_ok "All results and intermediate files (including BAM and uBAM) successfully synced to Windows!"
fi

TOTAL_ELAPSED=$(_time_elapsed "$PIPELINE_START")

echo ""
if [[ ${#FAILED_STEPS[@]} -gt 0 ]]; then
    echo "${CLR_BOLD}${CLR_YELLOW}╔══════════════════════════════════════════════════════════════════╗${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_YELLOW}║  ⚠  PIPELINE COMPLETED WITH WARNINGS                            ║${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_YELLOW}║  Errors in steps: ${FAILED_STEPS[*]}$(printf '%*s' $((47 - ${#FAILED_STEPS[*]})) '')║${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_YELLOW}╚══════════════════════════════════════════════════════════════════╝${CLR_RESET}"
else
    echo "${CLR_BOLD}${CLR_GREEN}╔══════════════════════════════════════════════════════════════════╗${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}║  ✔  PIPELINE COMPLETED SUCCESSFULLY  ⏱ ${TOTAL_ELAPSED}$(printf '%*s' $((26 - ${#TOTAL_ELAPSED})) '')║${CLR_RESET}"
    echo "${CLR_BOLD}${CLR_GREEN}╠══════════════════════════════════════════════════════════════════╣${CLR_RESET}"
    printf "${CLR_BOLD}${CLR_GREEN}║  %-64s║${CLR_RESET}\n" "📁 Main results:"
    printf "${CLR_BOLD}${CLR_GREEN}║    %-62s║${CLR_RESET}\n" "VCF: out/chrY_MT_final.vcf.gz"
    printf "${CLR_BOLD}${CLR_GREEN}║    %-62s║${CLR_RESET}\n" "Annotation: out/annotation/chrY_MT_final.hg38_multianno.txt"
    printf "${CLR_BOLD}${CLR_GREEN}║    %-62s║${CLR_RESET}\n" "Pathogenic: out/annotation/pathogenic.tsv"
    printf "${CLR_BOLD}${CLR_GREEN}║    %-62s║${CLR_RESET}\n" "Haplogroups: out/haplogroups/haplogroups_summary.tsv"
    printf "${CLR_BOLD}${CLR_GREEN}║    %-62s║${CLR_RESET}\n" "Kinship: out/kinship/kinship_summary.tsv"
    printf "${CLR_BOLD}${CLR_GREEN}║    %-62s║${CLR_RESET}\n" "QC Report: out/multiqc/multiqc_report.html"
    echo "${CLR_BOLD}${CLR_GREEN}╚══════════════════════════════════════════════════════════════════╝${CLR_RESET}"
fi
echo ""
echo "${CLR_WHITE}  Log: ${LOG_FILE}${CLR_RESET}"
echo ""
