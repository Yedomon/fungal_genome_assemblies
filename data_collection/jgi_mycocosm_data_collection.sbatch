#!/bin/bash
#SBATCH --job-name=jgi_download
#SBATCH --output=logs/%x.%A.out
#SBATCH --error=logs/%x.%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --qos=long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ange.zoclanclounon@rothamsted.ac.uk

# ─────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────
WORKDIR="/home/data/fungal_longread_review_rp10796-14/00.data/jgi"
CATALOG_FILE="$WORKDIR/jgi_fungi_catalog.tsv"   # Your saved catalog file
PREFIX_FILE="$WORKDIR/prefix_list.txt"
DOWNLOAD_DIR="$WORKDIR/assemblies"
LOG_FILE="$WORKDIR/logs/jgi_download.log"
SLEEP_COUNT=200       # Pause after every N downloads
SLEEP_DURATION=30     # Pause duration in seconds

mkdir -p "$WORKDIR" "$DOWNLOAD_DIR" "$WORKDIR/logs"

# ─────────────────────────────────────────────────────────────────
# STEP 1: EXTRACT PORTAL IDs FROM THE SAVED CATALOG
# Column structure (tab-separated):
# ##  name  portal  NCBI Taxon  assembly length  #of genes  ...
#  1   2      3         4             5               6
# Portal ID is in column 3. Skip header lines starting with '#'.
# ─────────────────────────────────────────────────────────────────
echo ">>> Extracting portal IDs from catalog: $CATALOG_FILE" | tee -a "$LOG_FILE"

if [[ ! -f "$CATALOG_FILE" ]]; then
    echo "ERROR: Catalog file not found at $CATALOG_FILE. Exiting." | tee -a "$LOG_FILE"
    exit 1
fi

awk -F'\t' '!/^#/ && $3 != "" && $3 != "portal" {print $3}' "$CATALOG_FILE" > "$PREFIX_FILE"

TOTAL=$(wc -l < "$PREFIX_FILE")
echo ">>> $TOTAL portal IDs extracted to: $PREFIX_FILE" | tee -a "$LOG_FILE"

if [[ $TOTAL -eq 0 ]]; then
    echo "ERROR: No portal IDs extracted. Check the catalog file format." | tee -a "$LOG_FILE"
    exit 1
fi

# ─────────────────────────────────────────────────────────────────
# STEP 2: DOWNLOAD GENOME ASSEMBLIES IN A POLITE LOOP
# ─────────────────────────────────────────────────────────────────
echo ">>> Starting bulk download from JGI MycoCosm..." | tee -a "$LOG_FILE"
echo ">>> Will pause ${SLEEP_DURATION}s after every ${SLEEP_COUNT} downloads." | tee -a "$LOG_FILE"

DOWNLOAD_COUNT=0
ERROR_COUNT=0
SKIP_COUNT=0

while IFS= read -r prefix; do

    # Skip empty lines and comments
    if [[ -z "$prefix" || "$prefix" =~ ^# ]]; then
        continue
    fi

    OUTPUT_FILE="$DOWNLOAD_DIR/${prefix}_unmasked-assembly.fasta.gz"

    # Skip already downloaded files
    if [[ -s "$OUTPUT_FILE" ]]; then
        echo "  [SKIP] Already exists: $prefix" | tee -a "$LOG_FILE"
        SKIP_COUNT=$(( SKIP_COUNT + 1 ))
        continue
    fi

    echo "  [$(( DOWNLOAD_COUNT + 1 ))/$TOTAL] Downloading: $prefix" | tee -a "$LOG_FILE"

    URL="https://mycocosm.jgi.doe.gov/ext-api/mycocosm/download/${prefix}/unmasked-assembly"

    curl --get -L -s "$URL" --output "$OUTPUT_FILE"

    if [[ $? -eq 0 && -s "$OUTPUT_FILE" ]]; then
        echo "  [OK] Saved: $OUTPUT_FILE" | tee -a "$LOG_FILE"
        DOWNLOAD_COUNT=$(( DOWNLOAD_COUNT + 1 ))
    else
        echo "  [ERROR] Failed to download: $prefix" | tee -a "$LOG_FILE"
        rm -f "$OUTPUT_FILE"   # Remove empty/corrupt file
        ERROR_COUNT=$(( ERROR_COUNT + 1 ))
    fi

    # Polite pause every SLEEP_COUNT successful downloads
    if [[ $(( DOWNLOAD_COUNT % SLEEP_COUNT )) -eq 0 && $DOWNLOAD_COUNT -gt 0 ]]; then
        echo "--- Pause: $DOWNLOAD_COUNT downloads completed. Sleeping ${SLEEP_DURATION}s... ---" | tee -a "$LOG_FILE"
        sleep "$SLEEP_DURATION"
        echo "--- Resuming. ---" | tee -a "$LOG_FILE"
    fi

done < "$PREFIX_FILE"

# ─────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "════════════════════════════════════════════" | tee -a "$LOG_FILE"
echo "  Download complete." | tee -a "$LOG_FILE"
echo "  Successfully downloaded : $DOWNLOAD_COUNT" | tee -a "$LOG_FILE"
echo "  Skipped (exist)         : $SKIP_COUNT" | tee -a "$LOG_FILE"
echo "  Errors                  : $ERROR_COUNT" | tee -a "$LOG_FILE"
echo "  Output directory        : $DOWNLOAD_DIR" | tee -a "$LOG_FILE"
echo "════════════════════════════════════════════" | tee -a "$LOG_FILE"