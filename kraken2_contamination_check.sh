#!/usr/bin/env bash

set -euo pipefail

########################################
# DEFAULTS
########################################

INPUT=""
KRAKEN_DB=""
OUTDIR=""

THREADS=1
MAX_JOBS=1

########################################
# HELP
########################################

usage() {
cat <<EOF

Kraken2 assembly species/contamination screening

Usage:
  $(basename "$0") -i INPUT_DIR -db KRAKEN_DB -o OUTPUT_DIR [options]

Required:
  -i DIR       Input directory containing .fa, .fasta or .fna assemblies
  -db DIR      Path to Kraken2 database
  -o DIR       Output directory

Optional:
  -t INT       Threads per Kraken2 job [default: 1]
  -jobs INT    Number of Kraken2 jobs to run in parallel [default: 1]
  -h           Show this help message

Example:
  $(basename "$0") \
    -i fasta \
    -db /path/to/kraken2_database \
    -o kraken_output \
    -t 8 \
    -jobs 4

Output:
  OUTPUT_DIR/reports/                  Kraken2 reports
  OUTPUT_DIR/raw/                      Raw Kraken2 classifications
  OUTPUT_DIR/kraken_species_summary.tsv
  OUTPUT_DIR/kraken_confirmed_species.tsv
  OUTPUT_DIR/kraken_run.log

Notes:
  Kraken2 and a compatible Kraken2 database must already be installed.

EOF
}

########################################
# ARGUMENTS
########################################

while [[ $# -gt 0 ]]; do
    case "$1" in

        -i)
            INPUT="$2"
            shift 2
            ;;

        -db)
            KRAKEN_DB="$2"
            shift 2
            ;;

        -o)
            OUTDIR="$2"
            shift 2
            ;;

        -t)
            THREADS="$2"
            shift 2
            ;;

        -jobs)
            MAX_JOBS="$2"
            shift 2
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        *)
            echo "ERROR: Unknown option: $1"
            echo ""
            usage
            exit 1
            ;;
    esac
done

########################################
# CHECK INPUTS
########################################

if [[ -z "$INPUT" || -z "$KRAKEN_DB" || -z "$OUTDIR" ]]; then
    echo "ERROR: -i, -db and -o are required."
    echo ""
    usage
    exit 1
fi

if [[ ! -d "$INPUT" ]]; then
    echo "ERROR: Input directory does not exist:"
    echo "  $INPUT"
    exit 1
fi

if [[ ! -d "$KRAKEN_DB" ]]; then
    echo "ERROR: Kraken2 database directory does not exist:"
    echo "  $KRAKEN_DB"
    exit 1
fi

if ! command -v kraken2 >/dev/null 2>&1; then
    echo "ERROR: kraken2 was not found in PATH."
    echo "Install/activate Kraken2 before running this script."
    exit 1
fi

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: -t must be a positive integer."
    exit 1
fi

if ! [[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: -jobs must be a positive integer."
    exit 1
fi

########################################
# OUTPUT DIRECTORIES
########################################

mkdir -p "$OUTDIR/reports"
mkdir -p "$OUTDIR/raw"

LOG="$OUTDIR/kraken_run.log"

{
    echo "Kraken2 assembly screening"
    echo "Date: $(date)"
    echo "Input: $INPUT"
    echo "Database: $KRAKEN_DB"
    echo "Output: $OUTDIR"
    echo "Threads per job: $THREADS"
    echo "Parallel jobs: $MAX_JOBS"
    echo ""
} > "$LOG"

########################################
# FIND FASTA FILES
########################################

mapfile -d '' FASTAS < <(
    find "$INPUT" \
        -maxdepth 1 \
        -type f \
        \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
        -print0 |
    sort -z
)

N=${#FASTAS[@]}

if [[ "$N" -eq 0 ]]; then
    echo "ERROR: No .fa, .fasta or .fna files found in:"
    echo "  $INPUT"
    exit 1
fi

echo "Found $N assemblies."
echo "Threads per job: $THREADS"
echo "Parallel jobs: $MAX_JOBS"
echo ""

echo "Assemblies found: $N" >> "$LOG"

########################################
# CHECK FOR DUPLICATE SAMPLE NAMES
########################################

declare -A SEEN

for fasta in "${FASTAS[@]}"; do

    base=$(basename "$fasta")
    base=${base%.*}

    if [[ -n "${SEEN[$base]:-}" ]]; then
        echo "ERROR: Duplicate sample basename detected: $base"
        echo "Rename files so every assembly has a unique basename."
        exit 1
    fi

    SEEN["$base"]=1
done

########################################
# LIMIT PARALLEL JOBS
########################################

limit_jobs() {

    while [[ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]]; do
        sleep 1
    done
}

########################################
# RUN KRAKEN2
########################################

count=0

for fasta in "${FASTAS[@]}"; do

    limit_jobs

    count=$((count + 1))

    base=$(basename "$fasta")
    base=${base%.*}

    echo "[$count/$N] Starting $base"

    (
        kraken2 \
            --use-names \
            --db "$KRAKEN_DB" \
            --threads "$THREADS" \
            --report "$OUTDIR/reports/${base}.report" \
            --output "$OUTDIR/raw/${base}.kraken" \
            "$fasta"

        echo "[$(date)] Finished: $base" >> "$LOG"

    ) &
done

wait

echo ""
echo "All Kraken2 jobs finished."

########################################
# CREATE DETAILED SUMMARY
########################################

SUMMARY="$OUTDIR/kraken_species_summary.tsv"

printf \
"ID\tTOP_SPECIES\tTOP_SPECIES_PCT\tSECOND_SPECIES\tSECOND_SPECIES_PCT\tC_JEJUNI_PCT\tC_COLI_PCT\tUNCLASSIFIED_PCT\n" \
> "$SUMMARY"

for fasta in "${FASTAS[@]}"; do

    base=$(basename "$fasta")
    base=${base%.*}

    report="$OUTDIR/reports/${base}.report"

    awk -F'\t' -v id="$base" '

    {
        pct  = $1 + 0
        rank = $4
        name = $6

        sub(/^[[:space:]]+/, "", name)
        sub(/[[:space:]]+$/, "", name)

        if (rank == "U") {
            unclassified = pct
        }

        if (rank == "S") {

            if (name == "Campylobacter jejuni") {
                jejuni = pct
            }

            if (name == "Campylobacter coli") {
                coli = pct
            }

            if (pct > top1) {

                top2 = top1
                species2 = species1

                top1 = pct
                species1 = name

            } else if (pct > top2) {

                top2 = pct
                species2 = name
            }
        }
    }

    END {

        if (species1 == "")
            species1 = "-"

        if (species2 == "")
            species2 = "-"

        printf "%s\t%s\t%.4f\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n",
               id,
               species1,
               top1,
               species2,
               top2,
               jejuni,
               coli,
               unclassified
    }

    ' "$report" >> "$SUMMARY"

done

########################################
# CREATE SIMPLE SPECIES SUMMARY
########################################

CONFIRMED="$OUTDIR/kraken_confirmed_species.tsv"

printf "ID\tKRAKEN_SPECIES\n" > "$CONFIRMED"

awk -F'\t' '
NR > 1 {
    print $1 "\t" $2
}
' "$SUMMARY" >> "$CONFIRMED"

########################################
# FINISH
########################################

{
    echo ""
    echo "Completed: $(date)"
    echo "Assemblies analysed: $N"
    echo "Detailed summary: $SUMMARY"
    echo "Species summary: $CONFIRMED"
} >> "$LOG"

echo ""
echo "========================================"
echo " Kraken2 screening complete"
echo "========================================"
echo ""
echo "Assemblies analysed: $N"
echo ""
echo "Detailed summary:"
echo "  $SUMMARY"
echo ""
echo "Simple species summary:"
echo "  $CONFIRMED"
echo ""
echo "Individual Kraken2 reports:"
echo "  $OUTDIR/reports/"
echo ""
echo "Raw Kraken2 output:"
echo "  $OUTDIR/raw/"
echo ""
echo "Run log:"
echo "  $LOG"
echo ""
