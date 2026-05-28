#!/usr/bin/env bash
# =============================================================================
#  steps/03_qc.sh — Step 3: Quality Control (QC)
# =============================================================================
#  Performs the following:
#    1. FastQC           → read quality (duplicates, GC-content, adapters)
#    2. QualiMap bamqc   → coverage, alignment, insert sizes
#    3. MultiQC          → summary HTML report
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${ROOT_DIR}/config.sh"
source "${ROOT_DIR}/lib/utils.sh"

STEP_ID="03_qc"
init_log "$STEP_ID"

if is_done "$STEP_ID"; then
    log_info "Step ${STEP_ID} is already completed. Skipping."
    exit 0
fi

log_step "03" "Quality Control: FastQC → QualiMap → MultiQC"

# Activate QC environment
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CONDA_ENV_QC}" 2>/dev/null || \
    log_warn "Failed to activate ${CONDA_ENV_QC}"

# =============================================================================
#  3.1 FASTQC
# =============================================================================
log_substep "FastQC: quality of final BAM files"

log_theory "FastQC — what are we analyzing?" \
    "FastQC is run on the FINAL BAMs (after BQSR) to ensure" \
    "that alignment and processing were performed correctly." \
    "" \
    "Key metrics:" \
    "  Per-base quality: Q30+ on >80% of positions — normal for modern data" \
    "  Duplicate rate:   already marked by MarkDup, FastQC shows total %" \
    "  GC-content:       should match reference (~41% for hg38)" \
    "  Adapter content:  should be minimal after trimming" \
    "  Insert size:      peak ≈150–250 bp for standard WGS/WES libraries" \
    "" \
    "For the Y-chromosome, we expect lower GC (~37%) — it is AT-rich."

log_cmd "fastqc ${OUT_FINAL_BAM}/*.bam -o ${OUT_FASTQC} -t ${PARALLEL_JOBS_LIGHT}"

    ${FASTQC_BIN} "${OUT_FINAL_BAM}"/*.bam -o "${OUT_FASTQC}" -t "${PARALLEL_JOBS_LIGHT}"

log_ok "FastQC completed"

# =============================================================================
#  3.2 QUALIMAP
# =============================================================================
log_substep "QualiMap bamqc: coverage and alignment statistics"

log_theory "QualiMap — detailed coverage analysis" \
    "QualiMap analyzes BAM and generates detailed reports:" \
    "" \
    "  Coverage distribution: coverage uniformity across the chromosome" \
    "    → for Y: expect peaks in repeat regions (segmental duplications)" \
    "  Genome fraction coverage: % of positions with coverage ≥ X" \
    "    → normal WGS: >90% of positions with coverage ≥10x" \
    "  Insert size: distribution of fragment lengths" \
    "    → peak at 150–350 bp for WGS, 100–200 bp for WES" \
    "  Mapping quality: % of reads with MAPQ≥30" \
    "    → PAR regions of Y have MAPQ=0 (align to both X and Y)" \
    "  GC bias: relationship between coverage and GC content" \
    "" \
    "Flag -gff chrY_MT.bed: analyze ONLY target regions" \
    "Flag -ip (inside provided regions): statistics inside BED"

# QualiMap requires unset DISPLAY (otherwise crashes on servers without X11)
unset DISPLAY 2>/dev/null || true

QUALIMAP_CMD='qualimap bamqc \
    --java-mem-size='"${SPARK_MEMORY}"' \
    --bam {} \
    --genome-gc-distr hg38 \
    -outdir '"${OUT_QUALIMAP}"'/{/.} \
    -outfile {/.} \
    -ip \
    -gff '"${CHRY_MT_BED}"' \
    -outformat HTML'

log_cmd "parallel --progress -j ${PARALLEL_JOBS} '${QUALIMAP_CMD}' ::: ${OUT_FINAL_BAM}/*.bam"

    ${GNU_PARALLEL} --progress -j "${PARALLEL_JOBS}" \
        "qualimap bamqc \
            --java-mem-size=${SPARK_MEMORY} --bam {} \
            --genome-gc-distr hg38 \
            -outdir ${OUT_QUALIMAP}/{/.} -outfile {/.} \
            -ip -gff ${CHRY_MT_BED} -outformat HTML" \
        ::: "${OUT_FINAL_BAM}"/*.bam

log_ok "QualiMap completed"

# =============================================================================
#  3.3 MULTIQC — summary report
# =============================================================================
log_substep "MultiQC: summary report for all samples"

log_theory "MultiQC — QC metric aggregation" \
    "MultiQC automatically finds all QC reports in a directory (FastQC, QualiMap," \
    "bcftools stats, GATK metrics, Picard) and aggregates them into a single HTML report." \
    "" \
    "Especially useful when analyzing >10 samples: allows quick" \
    "identification of outliers across any metric." \
    "" \
    "Run from out/ — MultiQC will recursively find all QC files."

log_cmd "cd ${OUT_DIR} && multiqc ./ -o ${OUT_MULTIQC}"

    cd "${OUT_DIR}" && ${MULTIQC_BIN} ./ -o "${OUT_MULTIQC}"
    cd "${BASE_DIR}"

log_ok "Step 03 completed. QC reports → ${OUT_MULTIQC}/"
mark_done "$STEP_ID"
