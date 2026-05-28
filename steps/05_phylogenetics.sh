#!/usr/bin/env bash
# =============================================================================
#  steps/05_phylogenetics.sh — Step 5: Haplogroups and Phylogenetics
# =============================================================================
#  Performs the following:
#    1. LineageTracker → Y-chromosome haplogroups (classify + phylo)
#    2. Yleaf          → Alternative Y-analysis from BAM
#    3. Haplogrep      → mtDNA haplogroups
#    4. mtDNA-server   → Detailed mitochondrial genome analysis
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${ROOT_DIR}/config.sh"
source "${ROOT_DIR}/lib/utils.sh"

STEP_ID="05_phylogenetics"
init_log "$STEP_ID"

if is_done "$STEP_ID"; then
    log_info "Step ${STEP_ID} is already completed. Skipping."
    exit 0
fi

log_step "05" "Haplogroups and Phylogenetics: Y-DNA + mtDNA"

ensure_dir "${OUT_PHYLO}"
ensure_dir "${OUT_HAPLOGROUPS}"

INPUT_VCF="${OUT_DIR}/chrY_MT_final.vcf.gz"
require_file "${INPUT_VCF}" "Final VCF"

# Rename chromosomes to standard chrY and chrM for compatibility with LineageTracker and Haplogrep
log_info "Renaming chromosomes in VCF for compatibility with LineageTracker and Haplogrep..."
MAPPED_VCF="${OUT_DIR}/chrY_MT_final_renamed.vcf.gz"
printf 'ENA|CM000686|CM000686.2 chrY\nENA|J01415|J01415.2 chrM\n' > "${OUT_DIR}/chr_map.txt"
bcftools annotate --rename-chrs "${OUT_DIR}/chr_map.txt" "${INPUT_VCF}" -O z -o "${MAPPED_VCF}"
bcftools index -t "${MAPPED_VCF}"
INPUT_VCF="${MAPPED_VCF}"

# =============================================================================
#  5.1 THEORETICAL INTRODUCTION
# =============================================================================

log_theory "Y-chromosome and mtDNA as time machines" \
    "Both loci DO NOT recombine (with rare exceptions), therefore:" \
    "  • Mutations accumulate linearly over time (molecular clock)" \
    "  • Each new mutation → a new branch on the phylogenetic tree" \
    "  • Haplogroup = an 'address' on this tree" \
    "" \
    "Y-chromosome (paternal lineage):" \
    "  Passed from father to son unchanged (except for mutations)." \
    "  Y-tree: A0, A1, B, CF, DE, F, G, H, I, J, K, L, M, N, O, P, Q, R..." \
    "  R1b — most common in Western Europe (~50–80%)" \
    "  R1a — dominates in Eastern Europe and South Asia" \
    "  I2  — Balkan origin, frequent in Serbia, Croatia" \
    "" \
    "mtDNA (maternal lineage):" \
    "  Inherited only through the egg cell (matriarchal line)." \
    "  Circular molecule 16,569 bp, ~37 genes." \
    "  mtDNA tree: L0-L6 (African) → M, N → all others" \
    "  H — most common in Europe (~40–50%)" \
    "  U5b — oldest European haplogroup (hunter-gatherers)" \
    "" \
    "Analysis Tools:" \
    "  LineageTracker → Y-DNA, YFull database (most complete Y tree)" \
    "  Yleaf          → Y-DNA from BAM, ISOGG markers" \
    "  Haplogrep      → mtDNA, PhyloTree/HaploGrep3 database" \
    "  mtDNA-server2  → detailed mtDNA heteroplasmy analysis (Nextflow)"

# =============================================================================
#  5.2 LINEAGETRACKER (Y-chromosome)
# =============================================================================
log_substep "LineageTracker: Y-haplogroup classification"

log_theory "LineageTracker — classification algorithm" \
    "1. Reads VCF and extracts SNP markers on chrY" \
    "2. Compares with haplogroup tree (.hg file or built-in database)" \
    "3. Determines the tree branch for each sample:" \
    "   --classify: assign haplogroup → .hapresult.hg file" \
    "   --phylo:    build phylogenetic tree → Newick/SVG" \
    "" \
    "Flags:" \
    "  --vcf    → input VCF with Y-markers" \
    "  -b 38    → genome build hg38 (marker coordinates)" \
    "  --snp-only → use only SNPs (no indels)" \
    "  -a       → analyze all samples in VCF" \
    "  -o       → output files prefix"

log_cmd "conda activate ${CONDA_ENV_LINEAGE}"
log_cmd "LineageTracker classify --vcf ${INPUT_VCF} -b 38 --snp-only -a -o ${OUT_PHYLO}/classify_Y"
log_cmd "LineageTracker phylo --hg ${OUT_PHYLO}/classify_Y.hapresult.hg --seq ${INPUT_VCF} --seq-format vcf -o ${OUT_PHYLO}/classify_Y"

    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate "${CONDA_ENV_LINEAGE}" 2>/dev/null || \
        log_warn "Failed to activate ${CONDA_ENV_LINEAGE}"

    ${LINEAGETRACKER} classify \
        --vcf "${INPUT_VCF}" \
        -b "${GENOME_BUILD##hg}" \
        --snp-only -a \
        -o "${OUT_PHYLO}/classify_Y"

    num_samples=$(bcftools query -l "${INPUT_VCF}" | wc -l)
    if [[ "$num_samples" -gt 1 ]]; then
        log_info "More than 1 sample found. Running LineageTracker phylo..."
        ${LINEAGETRACKER} phylo \
            --hg "${OUT_PHYLO}/classify_Y.hapresult.hg" \
            --seq "${INPUT_VCF}" \
            --seq-format vcf \
            -o "${OUT_PHYLO}/classify_Y"
    else
        log_info "Only 1 sample in VCF. Skipping phylogenetic tree construction (phylo)."
        sample_name=$(bcftools query -l "${INPUT_VCF}" | head -n 1)
        echo "(${sample_name}:0.001);" > "${OUT_PHYLO}/classify_Y.nwk"
    fi

    conda deactivate 2>/dev/null || true

log_ok "LineageTracker completed"

# =============================================================================
#  5.3 YLEAF (Y from BAM)
# =============================================================================
log_substep "Yleaf: Y-haplogroups from BAM files"

log_theory "Yleaf — Y-analysis from BAM" \
    "Yleaf works DIRECTLY with BAM files (without intermediate VCF)." \
    "This is useful when Y-chromosome coverage is low (aDNA, degraded DNA)." \
    "" \
    "Algorithm:" \
    "  1. Extracts pileup at ISOGG marker positions" \
    "  2. Calculates allele frequencies at each marker site" \
    "  3. Applies ISOGG tree classifier" \
    "  4. Result: haplogroup + confidence level"

YLEAF_CMD='Yleaf -bam {} -o '"${OUT_HAPLOGROUPS}"'/Yleaf/{/.} --reference_genome hg38'
log_cmd "conda activate ${CONDA_ENV_YLEAF}"
log_cmd "parallel --progress -j ${PARALLEL_JOBS} '${YLEAF_CMD}' ::: ${OUT_FINAL_BAM}/*.bam"

    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate "${CONDA_ENV_YLEAF}" 2>/dev/null || \
        log_warn "Failed to activate ${CONDA_ENV_YLEAF}"

    mkdir -p "${OUT_HAPLOGROUPS}/Yleaf"
    if ! ${GNU_PARALLEL} --progress -j "${PARALLEL_JOBS}" \
        "${YLEAF_BIN} -bam {} -o ${OUT_HAPLOGROUPS}/Yleaf/{/.} --reference_genome hg38" \
        ::: "${OUT_FINAL_BAM}"/*.bam; then
        log_warn "Yleaf finished with an error (possibly reference database is not configured). Skipping Yleaf."
    fi

    conda deactivate 2>/dev/null || true

log_ok "Yleaf completed"

# =============================================================================
#  5.4 HAPLOGREP (mtDNA)
# =============================================================================
log_substep "Haplogrep: mtDNA haplogroup classification"

log_theory "Haplogrep and PhyloTree" \
    "mtDNA does not recombine → all mutations accumulate in a single lineage." \
    "PhyloTree — reference mtDNA haplogroup tree (>5000 branches)." \
    "" \
    "Haplogrep algorithm:" \
    "  1. Reads VCF, extracts variants on chrM (mtDNA)" \
    "  2. Compares with PhyloTree via maximum likelihood algorithm" \
    "  3. Finds the branch with maximum marker overlap" \
    "  4. Outputs haplogroup + Quality Score (0-1)" \
    "" \
    "mtDNA specificities:" \
    "  • Heteroplasmy: presence of two mtDNA variants in one cell" \
    "    (mutant molecule frequency < 100%)" \
    "  • Important for diseases: LHON, MELAS, MERRF" \
    "  • Haplogrep accounts for heteroplasmy via VCF AF fields" \
    "" \
    "mtDNA haplogroup interpretation:" \
    "  H, H1, H2... → Western Europe, ~40-50% of Europeans" \
    "  U, U5b...    → pre-agricultural European hunter-gatherers" \
    "  J, T...      → Near Eastern origin, Neolithic farmers" \
    "  L...         → African lineages"

log_cmd "conda activate ${CONDA_ENV_HAPLOGREP}"
log_cmd "haplogrep classify --in ${INPUT_VCF} --format vcf --out ${OUT_HAPLOGROUPS}/haplogroups.txt"

    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate "${CONDA_ENV_HAPLOGREP}" 2>/dev/null || \
        log_warn "Failed to activate ${CONDA_ENV_HAPLOGREP}"

    ${HAPLOGREP_BIN} classify \
        --in "${INPUT_VCF}" \
        --format vcf \
        --out "${OUT_HAPLOGROUPS}/haplogroups.txt"

    conda deactivate 2>/dev/null || true

log_ok "Haplogrep completed. mtDNA haplogroups → ${OUT_HAPLOGROUPS}/haplogroups.txt"

# =============================================================================
#  5.5 mtDNA-SERVER-2 (Nextflow)
# =============================================================================
log_substep "mtDNA-server2: detailed heteroplasmy analysis (Nextflow)"

log_theory "mtDNA-server2 — detailed mtDNA analysis" \
    "mtDNA-server2 — specialized Nextflow pipeline for mtDNA:" \
    "  • Alignment to MT-specific reference (rCRS)" \
    "  • Heteroplasmy detection with ~1% threshold" \
    "  • Coverage statistics and quality scores" \
    "  • Annotation via MitoMap (mtDNA diseases)" \
    "  • Runs via Singularity container (dependency isolation)" \
    "" \
    "Nextflow DSL2 — workflow language for reproducible bioinformatics pipelines." \
    "Singularity — container for HPC (no root privileges, unlike Docker)."

log_cmd "conda activate ${CONDA_ENV_NEXTFLOW}"
log_cmd "cp /home/prep01/data/mtdna_server.config ${BASE_DIR}/"
log_cmd "nextflow run /opt/mtdna-server-2 -c mtdna_server.config -profile singularity"

    if command -v "${NEXTFLOW_BIN}" &>/dev/null; then
        # shellcheck disable=SC1091
        source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
        conda activate "${CONDA_ENV_NEXTFLOW}" 2>/dev/null || true

        cp /home/prep01/data/mtdna_server.config "${BASE_DIR}/" 2>/dev/null || \
            log_warn "mtdna_server.config not found. Skipping mtDNA-server."

        if [[ -f "${BASE_DIR}/mtdna_server.config" ]]; then
            ${NEXTFLOW_BIN} run /opt/mtdna-server-2 \
                -c "${BASE_DIR}/mtdna_server.config" \
                -profile singularity \
                --outdir "${OUT_DIR}/mtdna_server_chrM"
        fi
        conda deactivate 2>/dev/null || true
    else
        log_warn "Nextflow not found. Skipping mtDNA-server2."
    fi

# =============================================================================
#  5.6 HAPLOGROUP SUMMARY TABLE
# =============================================================================
log_substep "Generating haplogroup summary table"

${PYTHON3} - << 'PYEOF'
"""
Merges LineageTracker and Haplogrep results into a single table.
"""
import os
import csv
from pathlib import Path

out_phylo     = Path(os.environ.get("OUT_PHYLO", "out/phylogenetics"))
out_haplo     = Path(os.environ.get("OUT_HAPLOGROUPS", "out/haplogroups"))

# Read Y haplogroups
y_hg = {}
lineage_file = out_phylo / "classify_Y.hapresult.hg"
if lineage_file.exists():
    with open(lineage_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid = row.get("SampleID", "")
            y_hg[sid] = row.get("Haplogroup", "N/A")

# Read mtDNA haplogroups
mt_hg = {}
haplo_file = out_haplo / "haplogroups.txt"
if haplo_file.exists():
    with open(haplo_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid = row.get("SampleID", "")
            mt_hg[sid] = row.get("Haplogroup", "N/A")

# All samples
all_samples = set(list(y_hg.keys()) + list(mt_hg.keys()))

# Write summary table
summary_file = out_haplo / "haplogroups_summary.tsv"
with open(summary_file, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["SampleID", "Y_Haplogroup", "MT_Haplogroup", "Paternal_origin", "Maternal_origin"])
    for sid in sorted(all_samples):
        y = y_hg.get(sid, "N/A")
        mt = mt_hg.get(sid, "N/A")
        # Brief interpretation
        paternal = "Western European" if y.startswith("R1b") else \
                   "Eastern European/Central Asian" if y.startswith("R1a") else \
                   "Balkan" if y.startswith("I2") else "Other"
        maternal = "Western European" if mt.startswith("H") else \
                   "Pre-Neolithic European" if mt.startswith("U5") else \
                   "Near Eastern" if mt.startswith(("J","T")) else "Other"
        writer.writerow([sid, y, mt, paternal, maternal])

print(f"\n  Haplogroup summary table → {summary_file}")
print()
print(f"  {'SampleID':<15} {'Y-haplogroup':<25} {'mT-haplogroup':<20} {'Paternal':<20}")
print(f"  {'-'*80}")
for sid in sorted(all_samples):
    y = y_hg.get(sid, "N/A")
    mt = mt_hg.get(sid, "N/A")
    paternal = "W.European" if y.startswith("R1b") else "Balkan" if y.startswith("I2") else "Other"
    print(f"  {sid:<15} {y:<25} {mt:<20} {paternal:<20}")
print()

PYEOF

log_ok "Step 05 completed. Haplogroups → ${OUT_HAPLOGROUPS}/"
mark_done "$STEP_ID"
