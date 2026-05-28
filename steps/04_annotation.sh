#!/usr/bin/env bash
# =============================================================================
#  steps/04_annotation.sh — Step 4: Variant Annotation (ANNOVAR)
# =============================================================================
#  Performs the following:
#    1. ANNOVAR table_annovar.pl → multi-layer VCF annotation
#    2. Python filtering → selection of clinically significant variants
#    3. Final report with ClinVar classification and gnomAD frequencies
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${ROOT_DIR}/config.sh"
source "${ROOT_DIR}/lib/utils.sh"

STEP_ID="04_annotation"
init_log "$STEP_ID"

if is_done "$STEP_ID"; then
    log_info "Step ${STEP_ID} is already completed. Skipping."
    exit 0
fi

log_step "04" "Variant Annotation: ANNOVAR + ClinVar + gnomAD"

ensure_dir "${OUT_ANNOTATION}"
INPUT_VCF="${OUT_DIR}/chrY_MT_final.vcf.gz"

# =============================================================================
#  4.1 ANNOVAR ANNOTATION
# =============================================================================
log_substep "ANNOVAR: multi-layer variant annotation"

log_theory "ANNOVAR — annotation layers" \
    "Every variant in VCF is just a position and alleles. To interpret" \
    "it clinically, we need to overlay several layers of information:" \
    "" \
    "  1. refGene (gene-based):" \
    "     Where is the variant? Exon, intron, UTR, intergenic, splice-site?" \
    "     If in coding region: synonymous (silent) or nonsynonymous?" \
    "" \
    "  2. avsnp150 (filter-based):" \
    "     Is there an rsID in dbSNP150? Helps reference known variants." \
    "" \
    "  3. clinvar_20200316 (filter-based):" \
    "     Clinical significance according to ClinVar:" \
    "       Pathogenic | Likely_pathogenic | VUS | Likely_benign | Benign" \
    "     VUS = Variant of Uncertain Significance." \
    "" \
    "  4. gnomad30_genome / gnomad211_exome (filter-based):" \
    "     Allele frequency in population database gnomAD (125k+ genomes)." \
    "     A variant with MAF > 1% is likely a polymorphism, not disease-causing." \
    "" \
    "  5. dbnsfp35c (filter-based):" \
    "     Functional predictors for amino acid substitutions:" \
    "       SIFT    → tolerated / deleterious (evolution-based)" \
    "       PolyPhen-2 → benign / possibly_damaging / probably_damaging" \
    "       CADD score → combined damage score (>20 = top 1% deleterious)" \
    "       REVEL   → ensemble predictor for missense" \
    "" \
    "  6. dbscsnv11 (filter-based):" \
    "     Prediction of impact on splicing."

require_file "${INPUT_VCF}" "Final VCF"

ANNOVAR_CMD="${PERL} ${ANNOVAR_DIR}/table_annovar.pl \
    ${INPUT_VCF} \
    ${HUMANDB_DIR}/ \
    -buildver ${GENOME_BUILD} \
    -out ${OUT_ANNOTATION}/chrY_MT_final \
    -protocol ${ANNOVAR_PROTOCOL} \
    -operation ${ANNOVAR_OPERATION} \
    -nastring . \
    -vcfinput"

log_cmd "${ANNOVAR_CMD}"

if [[ ! -f "${ANNOVAR_DIR}/table_annovar.pl" ]]; then
    log_error "ANNOVAR not found at: ${ANNOVAR_DIR}/table_annovar.pl"
    exit 1
fi

${ANNOVAR_CMD}

log_ok "ANNOVAR annotation completed"

# =============================================================================
#  4.2 FILTERING BY gnomAD AND ClinVar (Python)
# =============================================================================
log_substep "Filtering: gnomAD MAF + ClinVar pathogenicity"

log_theory "Variant Filtering Strategy" \
    "After annotation we have thousands of variants. We select clinically significant ones:" \
    "" \
    "Step 1 — frequency filtering (gnomAD MAF):" \
    "  Remove variants with MAF > ${GNOMAD_MAF_THRESHOLD} (${GNOMAD_MAF_THRESHOLD} × 100 = $(echo "${GNOMAD_MAF_THRESHOLD} * 100" | bc 2>/dev/null || echo "threshold")%)" \
    "  Logic: disease-causing mutations are rare in the population." \
    "  Exception: recessive diseases can be more frequent (but Y/MT are haploid)." \
    "" \
    "Step 2 — functional filtering:" \
    "  Keep only exonic, splicing, UTR5 variants." \
    "  Intronic and intergenic variants require separate analysis." \
    "" \
    "Step 3 — ClinVar:" \
    "  Pathogenic / Likely_pathogenic → to the final list." \
    "  VUS → to a separate list for manual review." \
    "" \
    "For mtDNA variants, HETEROPLASMY is important:" \
    "  Homoplasmy (=100% mutant molecules): symptoms often pronounced" \
    "  Heteroplasmy (<100%): severity is proportional to mutation level"

${PYTHON3} - << 'PYEOF'
"""
Inline Python filter for annotated variants.
Reads ANNOVAR output (.hg38_multianno.txt) and applies filters.
"""
import os
import sys
import csv
from pathlib import Path

# ─── Parameters (from bash environment variables) ────────────────────────────────
out_annotation = os.environ.get("OUT_ANNOTATION", "out/annotation")
gnomad_threshold = float(os.environ.get("GNOMAD_MAF_THRESHOLD", "0.01"))

anno_file = Path(out_annotation) / "chrY_MT_final.hg38_multianno.txt"

if not anno_file.exists():
    print(f"  [WARN] Annotation file not found: {anno_file}")
    print(f"  [WARN] Skipping filtering (run ANNOVAR first)")
    sys.exit(0)

# ─── Filtering ──────────────────────────────────────────────────────────────
all_variants = []
pathogenic = []
vus = []
rare_functional = []

with open(anno_file, "r", newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        all_variants.append(row)

        # Extract filtering metrics
        func = row.get("Func.refGene", ".")
        clinvar = row.get("CLNSIG", row.get("ClinVar_SIG", "."))
        gnomad_str = row.get("AF_popmax", row.get("gnomAD_genome_ALL", "."))

        # gnomAD Frequency
        try:
            gnomad_af = float(gnomad_str) if gnomad_str not in (".", "", "NA") else 0.0
        except ValueError:
            gnomad_af = 0.0

        # Filter 1: rare variant (MAF < threshold)
        is_rare = gnomad_af < gnomad_threshold

        # Filter 2: functional (not intergenic/deep intronic)
        functional_funcs = {"exonic", "splicing", "exonic;splicing", "UTR5", "UTR3"}
        is_functional = func.strip() in functional_funcs

        # Filter 3: ClinVar pathogenicity
        is_pathogenic = any(term in clinvar for term in ["Pathogenic", "pathogenic"])
        is_likely_path = "Likely_pathogenic" in clinvar
        is_vus = "Uncertain" in clinvar or "VUS" in clinvar

        if is_pathogenic or is_likely_path:
            pathogenic.append(row)
        elif is_vus and is_rare:
            vus.append(row)
        elif is_rare and is_functional:
            rare_functional.append(row)

# ─── Writing results ───────────────────────────────────────────────────
def write_filtered(variants, filepath, label):
    if not variants:
        print(f"  [INFO] {label}: 0 variants")
        return
    fieldnames = list(variants[0].keys())
    with open(filepath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(variants)
    print(f"  [OK] {label}: {len(variants)} variants → {filepath}")

out_dir = Path(out_annotation)
write_filtered(pathogenic,     out_dir / "pathogenic.tsv",     "Pathogenic/Likely_pathogenic")
write_filtered(vus,            out_dir / "vus.tsv",            "VUS + rare")
write_filtered(rare_functional, out_dir / "rare_functional.tsv", "Rare functional")

# ─── Summary ──────────────────────────────────────────────────────────────────
print()
print("  ┌─── ANNOTATION SUMMARY ───────────────────────────────")
print(f"  │  Total variants:           {len(all_variants)}")
print(f"  │  Pathogenic/LP:            {len(pathogenic)}")
print(f"  │  VUS (rare):               {len(vus)}")
print(f"  │  Rare functional:          {len(rare_functional)}")
print(f"  │  MAF threshold (gnomAD):   {gnomad_threshold:.1%}")
print(f"  └──────────────────────────────────────────────────")

PYEOF

log_ok "Filtering completed. Results → ${OUT_ANNOTATION}/"

# =============================================================================
#  4.3 TOP PATHOGENIC VARIANTS OUTPUT
# =============================================================================
log_substep "Top pathogenic variants"

if [[ -f "${OUT_ANNOTATION}/pathogenic.tsv" ]]; then
    count=$(wc -l < "${OUT_ANNOTATION}/pathogenic.tsv")
    if [[ "$count" -gt 1 ]]; then
        log_info "Pathogenic variants found: $((count - 1))"
        echo ""
        echo "  Pathogenic variants:"
        awk -F'\t' 'NR<=6 {printf "  %-12s %-10s %-15s %-15s %-20s\n",
            $1,$2,$5,$6,$NF}' "${OUT_ANNOTATION}/pathogenic.tsv" 2>/dev/null || true
    fi
else
    log_info "Pathogenic variants file not found (possibly no pathogenic variants)"
fi

log_ok "Step 04 completed"
mark_done "$STEP_ID"
