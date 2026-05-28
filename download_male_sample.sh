#!/usr/bin/env bash
# =============================================================================
#  download_male_sample.sh — Fetch lightweight real male reads using prefetch & fasterq-dump
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/lib/utils.sh"

# Two verified GBR male samples from the 1000 Genomes Project
# SRR062634 = HG00096 (British Male)
# SRR062654 = HG00101 (British Male)
SAMPLES=("SRR062634" "SRR062654")

# Lightweight spot limit (2 million reads) to keep downloads and alignment fast!
# (Set to empty "" if you want to download the entire full-size dataset)
MAX_READS="2000000"

echo "================================================================="
echo "  🧬 Downloading Male Samples via SRA Toolkit (prefetch + fastq-dump)"
echo "================================================================="

# 1. Activate conda environment and ensure sra-tools is installed
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate "${CONDA_ENV_QC}" 2>/dev/null || {
    echo "[WARN] Could not activate ${CONDA_ENV_QC}. Attempting with current environment..."
}

if ! command -v prefetch &>/dev/null || ! command -v fastq-dump &>/dev/null; then
    echo "SRA Toolkit (prefetch / fastq-dump) not found."
    echo "Installing 'sra-tools' inside your Conda environment automatically..."
    conda install -y -n "${CONDA_ENV_QC}" -c bioconda sra-tools || {
        echo "[ERROR] Failed to install sra-tools. Please install manually or run: conda install -c bioconda sra-tools"
        exit 1
    }
fi

echo "✔ SRA Toolkit verified successfully."

# 2. Reference Genome verification
REF_FA="${REFERENCE_DIR}/GCA_000001405.29.fasta"
if [ ! -f "$REF_FA" ]; then
    echo "[ERROR] Reference genome not found at: ${REF_FA}"
    exit 1
fi

for accession in "${SAMPLES[@]}"; do
    echo ""
    echo "-----------------------------------------------------------------"
    echo " Downloading SRA Run: ${accession}"
    echo "-----------------------------------------------------------------"

    # 3. Prefetch the SRA file
    echo "Running prefetch to fetch SRA archive..."
    prefetch "${accession}" -O "${TMP_DIR}"

    # 4. Extract FASTQ files using fastq-dump with spot limit and compression
    echo "Extracting reads using fastq-dump..."
    # Note: fastq-dump with --gzip --split-files directly produces accession_1.fastq.gz and accession_2.fastq.gz
    fastq-dump "${TMP_DIR}/${accession}/${accession}.sra" \
        -O "${FASTQ_DIR}" \
        --split-files \
        --gzip \
        -X "${MAX_READS}"

    # 5. Rename files to match the pipeline's expected format (*_1_trimmed.fastq.gz)
    mv "${FASTQ_DIR}/${accession}_1.fastq.gz" "${FASTQ_DIR}/${accession}_1_trimmed.fastq.gz"
    mv "${FASTQ_DIR}/${accession}_2.fastq.gz" "${FASTQ_DIR}/${accession}_2_trimmed.fastq.gz"

    # 7. Cleanup cached SRA archive
    echo "Cleaning up SRA cache..."
    rm -rf "${TMP_DIR}/${accession}"

    echo "✔ Successfully finished downloading & extracting ${accession}!"
done

echo ""
echo "================================================================="
echo "  ✔ SUCCESS! Two lightweight male datasets are ready: "
for accession in "${SAMPLES[@]}"; do
    echo "    - fastq/${accession}_1_trimmed.fastq.gz"
    echo "    - fastq/${accession}_2_trimmed.fastq.gz"
done
echo "================================================================="
