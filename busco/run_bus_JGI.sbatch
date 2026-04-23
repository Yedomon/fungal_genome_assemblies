#!/bin/bash
#SBATCH -p standard
#SBATCH --qos=medium
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -a 1-2776%10
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=bus_JGI
#SBATCH --output=logs/bus_%A_%a.out
#SBATCH --error=logs/bus_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ange.zoclanclounon@rothamsted.ac.uk

# ---------------------------
# Environment setup
# ---------------------------
source /home/data/cs_genomes/analysis/Ange/00.Tools/miniconda3/bin/activate
conda activate buscov6_env

working_dir=/home/data/fungal_longread_review_rp10796-14/01.BUSCO/JGI
genome_data_dir=/home/data/fungal_longread_review_rp10796-14/00.data/JGI

cd $working_dir

# ---------------------------
# Load genome name from list
# ---------------------------
GENOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $genome_data_dir/genome_list_2776.txt)
GENOME=$(echo "$GENOME" | xargs)   
GENOME="$genome_data_dir/$GENOME"


# ---------------------------
# Prepare filenames
# ---------------------------
BNAME=$(basename "$GENOME" .fna)
OUTDIR="busco_${BNAME}_results"

echo "Running BUSCO for: $GENOME"
echo "Output to: $OUTDIR"

# ---------------------------
# Run BUSCO
# ---------------------------
busco \
    -i "$GENOME" \
    --offline \
    -m genome \
    -l /home/data/fungal_longread_review_rp10796-14/00.data/busco_fungi_database/fungi_odb12 \
    -c 16 \
    -o "$OUTDIR"

echo "BUSCO completed for $GENOME"

