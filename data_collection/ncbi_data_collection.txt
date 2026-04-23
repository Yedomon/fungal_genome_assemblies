#!/bin/bash
#SBATCH --job-name=download_extract
#SBATCH --output=logs/%x.%A.out
#SBATCH --error=logs/%x.%A.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --qos=long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ange.zoclanclounon@rothamsted.ac.uk

# ─────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────
ACCESSION_FILE="/home/data/fungal_longread_review_rp10796-14/00.data/ncbi_genome_accession_list.txt"
ZIPDIR="/home/data/fungal_longread_review_rp10796-14/00.data"
OUTDIR="/home/data/fungal_longread_review_rp10796-14/00.data/extracted"
CONDA_BASE="/home/data/cs_genomes/analysis/Ange/00.Tools/miniconda3/bin/activate"
CONDA_ENV="ncbi_datasets_env"

mkdir -p "$OUTDIR" logs

# ─────────────────────────────────────────────
# STEP 1: DOWNLOAD GENOME ASSEMBLIES FROM NCBI
# ─────────────────────────────────────────────
echo ">>> Activating conda environment: $CONDA_ENV"
source "$CONDA_BASE"
conda activate "$CONDA_ENV"

echo ">>> Downloading genome assemblies..."
while read -r accession; do
    echo "  Downloading: $accession"
    datasets download genome accession "$accession" --filename "$ZIPDIR/$accession.zip"
done < "$ACCESSION_FILE"

conda deactivate
echo ">>> Download complete."

# ─────────────────────────────────────────────
# STEP 2: EXTRACT .fna FILES FROM ZIP ARCHIVES
# ─────────────────────────────────────────────
echo ">>> Extracting .fna files from ZIP archives..."

for zip in "$ZIPDIR"/*.zip; do
    echo "  Processing: $zip"

    # Detect the .fna file path inside the zip
    fna_path=$(unzip -Z1 "$zip" | grep "_genomic.fna$")

    if [[ -z "$fna_path" ]]; then
        echo "  WARNING: No *_genomic.fna found in $zip — skipping."
        continue
    fi

    # Extract using accession as the output filename
    base=$(basename "$zip" .zip)
    unzip -p "$zip" "$fna_path" > "$OUTDIR/${base}.fna"
    echo "  Extracted: $OUTDIR/${base}.fna"
done

echo ">>> Extraction complete. All files saved to: $OUTDIR"