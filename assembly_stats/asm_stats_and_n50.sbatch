#!/bin/bash
#SBATCH --job-name=assembly_stats
#SBATCH --output=logs/%x.%A.out
#SBATCH --error=logs/%x.%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ange.zoclanclounon@rothamsted.ac.uk

# ─────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────
CONDA_BASE="/home/data/cs_genomes/analysis/Ange/00.Tools/miniconda3/bin/activate"
CONDA_ENV="assembly_stats_env"
WORKDIR="/home/data/fungal_longread_review_rp10796-14/00.data/jgi"
DOWNLOAD_DIR="$WORKDIR/assemblies"
STATS_DIR="$WORKDIR/assembly_stats"
N50_SUMMARY="$WORKDIR/n50_summary.txt"

mkdir -p "$STATS_DIR" logs

# ─────────────────────────────────────────────────────────────────
# ACTIVATE CONDA
# ─────────────────────────────────────────────────────────────────
source "$CONDA_BASE"
conda activate "$CONDA_ENV"

# ─────────────────────────────────────────────────────────────────
# INITIALISE N50 SUMMARY FILE WITH HEADER
# ─────────────────────────────────────────────────────────────────
echo -e "assembly\tN50" > "$N50_SUMMARY"

# ─────────────────────────────────────────────────────────────────
# LOOP OVER ASSEMBLIES
# ─────────────────────────────────────────────────────────────────
echo ">>> Running assembly-stats on all FASTA files in $DOWNLOAD_DIR"

for fasta in "$DOWNLOAD_DIR"/*.fasta; do

    # Get base name e.g. Zoorad1_unmasked-assembly
    base=$(basename "$fasta" .fasta)

    echo "  Processing: $base"

    # Run assembly-stats and save full output
    assembly-stats "$fasta" > "$STATS_DIR/assemblystats.${base}.txt"

    if [[ $? -ne 0 ]]; then
        echo "  [ERROR] assembly-stats failed for $fasta"
        continue
    fi

    # Extract N50 value using awk
    n50=$(awk '/^N50 =/ {print $3}' "$STATS_DIR/assemblystats.${base}.txt" | tr -d ',')

    if [[ -z "$n50" ]]; then
        echo "  [WARNING] N50 not found for $base"
        n50="NA"
    fi

    # Append to N50 summary
    echo -e "${base}\t${n50}" >> "$N50_SUMMARY"

    echo "  [OK] N50 = $n50 → $STATS_DIR/assemblystats.${base}.txt"

done

conda deactivate

# ─────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────
TOTAL=$(( $(wc -l < "$N50_SUMMARY") - 1 ))
echo ""
echo "════════════════════════════════════════════"
echo "  Assembly stats complete."
echo "  Genomes processed : $TOTAL"
echo "  Individual stats  : $STATS_DIR/"
echo "  N50 summary table : $N50_SUMMARY"
echo "════════════════════════════════════════════"