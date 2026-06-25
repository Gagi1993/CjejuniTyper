#!/bin/bash
set -euo pipefail
shopt -s nullglob

########################################
# HELP
########################################

show_help() {
  cat <<EOF

CjejuniTyper — Campylobacter jejuni Typing Pipeline
══════════════════════════════════════════════════
Performs multi-layer genotypic characterisation of C. jejuni assemblies:

  • Penner (HS) capsular serotyping   — Two-step BLAST + in silico PCR system:
      Step 1: specific marker genes (blastn, 99.5%/80cov) — fast, precise
      Step 2: full-locus blastn fallback for unresolved samples
      Step 3: complex subtype bracketing (AA blastx / poly-G / PCR)
      Step 4: in silico PCR confirmation (ipcress)
  • LOS class typing                  — BLAST-based, BEST HIT selection (no mosaics)
  • Virulence factor detection        — 90% identity / 60% coverage
  • cstII allele typing               — Asn51 (bifunctional) vs Thr51 (monofunctional)
  • GBS risk prediction               — ganglioside mimic type + combined HS/LOS score
  • MEOP / T4SS / CDT status          — capsule modification + invasion factors
  • AMR detection                     — AMRFinder Plus (acquired genes + point mutations)
                                        gyrA T86I/D90N/P104S, 23S A2075G, tet(O), blaOXA, etc.
  • Assembly QC                       — genome size, GC%, contig count

HS Serotype Complexes
  Complexes are always reported as complex (e.g. HS4c, HS8_17c).
  Where subtype can be predicted, it appears in brackets (e.g. HS8_17c(HS8*)).
  The asterisk (*) denotes a predicted subtype, not confirmed by serology.

  HS4c     — HS4 complex: HS4, HS13, HS16, HS43, HS50, HS62, HS63, HS64, HS65
               Bracket: HS4c(HS16/64*) if Mu_HS4B primer hits
  HS6_7c   — HS6/7 complex: HS6, HS7
               No subtype bracket possible (>99% identical loci)
  HS8_17c  — HS8/17 complex: HS8, HS17
               Bracket: HS8_17c(HS8*) or HS8_17c(HS17*) via AA blastx
  HS5_31c  — HS5/31 complex: HS5, HS31
               Bracket: HS5_31c(HS31*) via Mu_HS5 in silico PCR
  HS23_36c — HS23/36 complex: HS23, HS36
               Bracket: HS23_36c(HS23*) or HS23_36c(HS36*) via poly-G/C tract
  HS1c     — HS1 complex: HS1, HS44
               No subtype bracket (single marker covers both)
  HS15c    — HS15 complex: HS15, HS31, HS58
               Bracket: HS15c(HS58*) via Mu_HS58 in silico PCR
  HS45c    — HS45 complex: HS45, HS32, HS60
               No subtype bracket possible

Output: final_los_hs_gbs.tsv with columns:
  ID  LOS_CLASS  SIALYLATED  HS_TYPE  GBS_RISK
  genome_size_bp  n_contigs  gc_pct  qc_status  qc_warnings
  CSTII_FUNCTION  GANGLIOSIDE_MIMIC  GBS_MIMIC_RISK  HS_GBS_MARKER  COMBINED_GBS
  SERUM_RESISTANCE  MEOP_STATUS  T4SS_STATUS  CDT_STATUS
  SIALYLATION_GENES  LOS_CORE_STATUS  CAPSULE_KPS  WLAN_STATUS
  <AMR gene columns — one per unique gene, %ID %COV>
  <raw VF gene columns>

Usage:
  $0 -i assemblies/ -o results/ [options]

Required:
  -i   Input directory containing .fa / .fasta / .fna assemblies
  -o   Output directory

Optional:
  -hs         HS full-locus BLAST database directory [default: databases/hs_db_new]
  -hs_aa      HS specific marker gene DB (blastn)    [default: databases/hs_specific_aa_per_hs/hs/hs_aa]
  -los        LOS class BLAST database directory     [default: databases/los_db_new]
  -los_b_windows LOS class B window database directory [default: databases/los_db_specific_markers/LOS_B_windows/]
  -vf         Virulence factor BLAST database dir    [default: databases/vfdb_new]
  -t          Threads per job                        [default: 32]
  -j          Parallel jobs                          [default: 1]
  -r          Resume — skip completed samples
  --skip_hs   Skip all HS (Penner) typing — HS_TYPE will be "-" for all samples
  --skip_qc    Skip assembly QC
  --skip_mlst  Skip MLST typing
  --skip_los   Skip LOS typing
  --skip_vf    Skip virulence factor typing
  --skip_amr   Skip AMRFinder detection
  -I          Install conda environment and exit
  -U          Update conda environment (blast + amrfinder) and exit
  -amrfinder_db  Force-update AMRFinder database and exit
  -h          Show this help

Examples:
  $0 -i assemblies/ -o results/
  $0 -i assemblies/ -o results/ -t 16 -j 4
  $0 -i assemblies/ -o results/ -r
  $0 -i assemblies/ -o results/ --skip_hs   # LOS + VF + AMR only
  $0 -I                    # install environment
  $0 -U                    # update blast + amrfinder packages
  $0 -amrfinder_db         # force-update AMRFinder database

Environment:
  Conda/mamba environment: cjejuni_typer
  Packages: blast, ncbi-amrfinderplus, exonerate (ipcress), seqkit
  Python 3 (system) used for interpretation — no conda needed

T4SS database setup (run once):
  grep -A1 "virB\|virD4" campy_gene.fas > virB_virD4.fasta
  mkdir -p databases/vfdb/blast/T4SS
  cp virB_virD4.fasta databases/vfdb/blast/T4SS/T4SS.fasta
  makeblastdb -in databases/vfdb/blast/T4SS/T4SS.fasta \\
              -dbtype nucl -out databases/vfdb/blast/T4SS/T4SS

LOS Typing Logic:
  - Uses BEST HIT scoring (identity × coverage × bitscore)
  - Thresholds: STRICT (99.5% ID, 93% cov), STRICT_2 (97% ID, 93% cov), STRICT_3 (95% ID, 93% cov), MID (94% ID, 80% cov), MID_2 (92% ID, 80% cov), LOOSE (90% ID, 60% cov)
  - True mosaics retained: I/J, I/S (both hits independently >90% coverage)
  - Thresholds: STRICT (99.5% ID, 99% cov), MID (95% ID, 80% cov), LOOSE (80% ID, 60% cov)

AMR Detection:
  - Detects acquired resistance genes (tet, bla, etc.)
  - Detects point mutations (gyrA, 23S rRNA, etc.)
  - Excludes STRESS genes (arsenic, etc.) from final output
  - Requires AMRFinder database update: $0 -amrfinder_db

HS Penner Serotype Complexes and Standalone Types:
  ┌─────────────┬──────────────────────────────────────────────────────────────────┬─────────────────────────────────────┐
  │ Complex     │ Member serotypes                                                  │ Bracket method                      │
  ├─────────────┼──────────────────────────────────────────────────────────────────┼─────────────────────────────────────┤
  │ HS4c        │ HS4, HS13, HS16, HS43, HS50, HS62, HS63, HS64, HS65             │ Mu_HS4B PCR → HS4c(HS16/64*)        │
  │             │                                                                   │ Blast winner → HS4c(HS4*) etc.      │
  │ HS6_7c      │ HS6, HS7                                                          │ None (>99% identical)               │
  │ HS8_17c     │ HS8, HS17                                                         │ AA blastx → HS8_17c(HS8/HS17*)      │
  │ HS5_31c     │ HS5, HS31                                                         │ Mu_HS5 PCR → HS5_31c(HS5/HS31*)     │
  │ HS23_36c    │ HS23, HS36                                                        │ Poly-G/C tract → HS23_36c(HS23/36*) │
  │ HS1c        │ HS1, HS44                                                         │ Blast winner → HS1c(HS1/HS44*)      │
  │ HS15c       │ HS15, HS31, HS58                                                  │ Mu_HS58 PCR → HS15c(HS58*)          │
  │ HS45c       │ HS45, HS32, HS60                                             │ Blast winner → HS45c(HS32/60*)      │
  └─────────────┴──────────────────────────────────────────────────────────────────┴─────────────────────────────────────┘

  Standalone types (not in a complex):
  HS2, HS3, HS9, HS10, HS11, HS12, HS18, HS19, HS21, HS22, HS27, HS29,
  HS33_35, HS37, HS38, HS40, HS41, HS42, HS52, HS53, HS55, HS57

  Note: * = predicted subtype, not confirmed by serology.
        Brackets are best-effort predictions from genomic data.
        HS23/36c subtype prediction requires accurate poly-G tract counting
        which may be unreliable in short-read draft assemblies.

EOF
}

########################################
# JOB CONTROL
########################################

limit_jobs() {
  local MAX_JOBS=$1
  while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do
    sleep 1
  done
}

########################################
# DEFAULTS
########################################

HS_DB="databases/new_hs_db/blast"
HS_MARKER="databases/hs_specific_aa_per_hs/hs/hs_aa"
LOS_DB="databases/los_db_new"
LOS_B_WINDOW_DB="databases/los_db_specific_markers/LOS_B_windows/LOS_B_WINDOWS_DB"
VF_DB="databases/vfdb_new"
ORF11_DB="databases/vfdb/orf11/orf11_db"
THREADS=32
JOBS=1
INSTALL=0
UPDATE=0
AMR_DB_UPDATE=0
RESUME=0
SKIP_HS=0
SKIP_QC=0
SKIP_MLST=0
SKIP_LOS=0
SKIP_VF=0
SKIP_AMR=0
INPUT=""
OUTDIR=""
ENV_NAME="cjejuni_typer"

# Detect package manager
if command -v mamba >/dev/null 2>&1; then
  PKG_MANAGER="mamba"
elif command -v conda >/dev/null 2>&1; then
  PKG_MANAGER="conda"
else
  PKG_MANAGER=""
fi

########################################
# ARGS
########################################

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i)      INPUT="$2";    shift 2 ;;
    -o)      OUTDIR="$2";   shift 2 ;;
    -hs)     HS_DB="$2";    shift 2 ;;
    -hs_aa)  HS_MARKER="$2"; shift 2 ;;
    -los)    LOS_DB="$2";   shift 2 ;;
    -los_b_windows) LOS_B_WINDOW_DB="$2"; shift 2 ;;
    -vf)     VF_DB="$2";    shift 2 ;;
    -orf11_db) ORF11_DB="$2"; shift 2 ;;
    -t)      THREADS="$2";  shift 2 ;;
    -j)      JOBS="$2";     shift 2 ;;
    -r)      RESUME=1;      shift ;;
    --skip_hs) SKIP_HS=1;  shift ;;
    --skip_qc)   SKIP_QC=1; shift ;;
    --skip_mlst) SKIP_MLST=1; shift ;;
    --skip_los)  SKIP_LOS=1; shift ;;
    --skip_vf)   SKIP_VF=1; shift ;;
    --skip_amr)  SKIP_AMR=1; shift ;;
    -I)      INSTALL=1;     shift ;;
    -U)      UPDATE=1;      shift ;;
    -amrfinder_db) AMR_DB_UPDATE=1; shift ;;
    -h|--help) show_help; exit 0 ;;
    *)    echo "Unknown option: $1"; show_help; exit 1 ;;
  esac
done

########################################
# INSTALL
########################################

if [ "$INSTALL" -eq 1 ]; then
  if [ -z "$PKG_MANAGER" ]; then
    echo "ERROR: conda or mamba not found."
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
  fi
  echo "Installing conda environment: $ENV_NAME"
  echo "Using package manager: $PKG_MANAGER"

  if $PKG_MANAGER env list | grep -q "^${ENV_NAME}"; then
    echo "Environment '$ENV_NAME' already exists. Skipping."
    echo "To reinstall: conda env remove -n $ENV_NAME && $0 -I"
  else
    $PKG_MANAGER create -y -n "$ENV_NAME" \
      -c conda-forge -c bioconda \
      blast ncbi-amrfinderplus exonerate seqkit
    echo ""
    echo "✓ Environment '$ENV_NAME' installed."
    echo "  Packages: blast, ncbi-amrfinderplus, exonerate (ipcress), seqkit"
    echo ""
    echo "Next: update the AMRFinder database:"
    echo "  $0 -amrfinder_db"
  fi
  exit 0
fi

########################################
# UPDATE ENV
########################################

if [ "$UPDATE" -eq 1 ]; then
  if [ -z "$PKG_MANAGER" ]; then
    echo "ERROR: conda or mamba not found."; exit 1; fi
  echo "Updating conda environment: $ENV_NAME ..."
  $PKG_MANAGER update -y -n "$ENV_NAME" \
    -c conda-forge -c bioconda \
    blast ncbi-amrfinderplus exonerate seqkit
  echo "✓ Environment updated."
  exit 0
fi

########################################
# UPDATE AMRFINDER DATABASE
########################################

if [ "$AMR_DB_UPDATE" -eq 1 ]; then
  echo "Force-updating AMRFinder database..."
  if [ -n "$PKG_MANAGER" ] && $PKG_MANAGER env list 2>/dev/null | grep -q "^${ENV_NAME}"; then
    $PKG_MANAGER run -n "$ENV_NAME" amrfinder -U
  elif command -v amrfinder >/dev/null 2>&1; then
    amrfinder -U
  else
    echo "ERROR: amrfinder not found. Run '$0 -I' first."; exit 1
  fi
  echo "✓ AMRFinder database updated."
  exit 0
fi

########################################
# VALIDATE INPUT
########################################

if [ "$#" -lt 1 ] && [ -z "$INPUT" ]; then
  show_help; exit 1
fi

[ -z "$INPUT" ]  && echo "ERROR: Missing -i (input directory)"  && exit 1
[ -z "$OUTDIR" ] && echo "ERROR: Missing -o (output directory)" && exit 1
[ -d "$INPUT" ]  || { echo "ERROR: Input directory not found: $INPUT"; exit 1; }

########################################
# BLAST COMMAND WRAPPER
# Runs blastn either via conda env or directly if already in PATH
########################################

run_blastn() {
  if [ -n "$PKG_MANAGER" ] && $PKG_MANAGER env list 2>/dev/null | grep -q "^${ENV_NAME}"; then
    $PKG_MANAGER run -n "$ENV_NAME" blastn "$@"
  elif command -v blastn >/dev/null 2>&1; then
    blastn "$@"
  else
    echo "ERROR: blastn not found. Run '$0 -I' to install the conda environment."
    exit 1
  fi
}

# ── Query sanitiser (FIX) ─────────────────────────────────────────────────
# seqkit is detected GLOBALLY here (not inside the HS block, which is skipped
# under --skip_hs). It is used by clean_query() to normalise assemblies before
# BLAST: strips CR (CRLF -> LF), drops blank lines, and — when seqkit is
# present — re-wraps/validates records so BLAST's CFastaReader stops throwing
# "Near line 2 ... doesn't look like plausible data" on malformed FASTA.
SEQKIT_CMD=""
if [ -n "${PKG_MANAGER:-}" ] && $PKG_MANAGER env list 2>/dev/null | grep -q "^${ENV_NAME}[[:space:]]"; then
  SEQKIT_CMD="$PKG_MANAGER run -n $ENV_NAME seqkit"
elif command -v seqkit >/dev/null 2>&1; then
  SEQKIT_CMD="seqkit"
fi

# Normalise a query FASTA to stdout for BLAST.
clean_query() {
  local qf="$1"
  if [ -n "${SEQKIT_CMD:-}" ]; then
    tr -d '\r' < "$qf" | $SEQKIT_CMD seq -w 0 2>/dev/null
  else
    tr -d '\r' < "$qf" | awk 'NF'
  fi
}

run_amrfinder() {
  if [ -n "$PKG_MANAGER" ] && $PKG_MANAGER env list 2>/dev/null | grep -q "^${ENV_NAME}[[:space:]]"; then
    $PKG_MANAGER run -n "$ENV_NAME" amrfinder "$@"
  elif command -v amrfinder >/dev/null 2>&1; then
    amrfinder "$@"
  else
    return 1
  fi
}

# Test amrfinder availability once at startup (not during run)
AMR_AVAILABLE=0
if [ -n "$PKG_MANAGER" ] && $PKG_MANAGER env list 2>/dev/null | grep -q "^${ENV_NAME}[[:space:]]"; then
  if $PKG_MANAGER run -n "$ENV_NAME" amrfinder --version >/dev/null 2>&1; then
    AMR_AVAILABLE=1
  fi
elif command -v amrfinder >/dev/null 2>&1; then
  AMR_AVAILABLE=1
fi

THREADS_PER_JOB=$(( THREADS / JOBS ))
[ "$THREADS_PER_JOB" -lt 1 ] && THREADS_PER_JOB=1

########################################
# DIRS
########################################

mkdir -p "$OUTDIR/hs_blast_strict" "$OUTDIR/hs_blast_mid" "$OUTDIR/hs_blast_loose"
mkdir -p "$OUTDIR/los_blast_strict" "$OUTDIR/los_blast_mid" "$OUTDIR/los_blast_loose"
mkdir -p "$OUTDIR/vf_blast"
mkdir -p "$OUTDIR/.done"   # tracks completed samples for resume

echo ""
echo "CampyTyper"
echo "══════════════════════════════════════"
echo "  Input       : $INPUT"
echo "  Output      : $OUTDIR"
echo "  HS DB       : $HS_DB"
echo "  HS Marker DB: $HS_MARKER"
echo "  LOS DB      : $LOS_DB"
echo "  VF DB       : $VF_DB"
echo "  ORF11 DB    : $ORF11_DB"
echo "  Threads     : $THREADS (${THREADS_PER_JOB}/job × ${JOBS} jobs)"
echo "  Resume      : $([ "$RESUME" -eq 1 ] && echo YES || echo NO)"
echo "  Skip HS     : $([ "$SKIP_HS" -eq 1 ] && echo YES || echo NO)"

echo "══════════════════════════════════════"
echo ""
echo "  Skip QC     : $([ "$SKIP_QC" -eq 1 ] && echo YES || echo NO)"
echo "  Skip MLST   : $([ "$SKIP_MLST" -eq 1 ] && echo YES || echo NO)"
echo "  Skip LOS    : $([ "$SKIP_LOS" -eq 1 ] && echo YES || echo NO)"
echo "  Skip VF     : $([ "$SKIP_VF" -eq 1 ] && echo YES || echo NO)"
echo "  Skip AMR    : $([ "$SKIP_AMR" -eq 1 ] && echo YES || echo NO)"
TIME_START=$(date +%s)
echo "  Started     : $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Helper: check if a sample is already done
is_done() { [ -f "$OUTDIR/.done/${1}.done" ]; }

# Helper: mark a sample as done and print progress
PROGRESS_FILE="$OUTDIR/.done/.progress"
mark_done() {
  touch "$OUTDIR/.done/${1}.done"
  # atomic increment using a lock file
  (
    flock 9
    count=$(cat "$PROGRESS_FILE" 2>/dev/null || echo 0)
    count=$(( count + 1 ))
    echo "$count" > "$PROGRESS_FILE"
    echo "  [${count}/${N_SAMPLES}] done: ${1}"
  ) 9>"$OUTDIR/.done/.lock"
}

########################################
# SAMPLE LIST
########################################

SAMPLES="$OUTDIR/samples.txt"
> "$SAMPLES"
for f in "$INPUT"/*.{fa,fasta,fna}; do
  [ -f "$f" ] || continue
  n=$(basename "$f"); n=${n%.*}
  echo "$n" >> "$SAMPLES"
done

N_SAMPLES=$(wc -l < "$SAMPLES")
echo "Samples found: $N_SAMPLES"
[ "$N_SAMPLES" -eq 0 ] && echo "ERROR: No .fa/.fasta/.fna files in $INPUT" && exit 1

if [ "$RESUME" -eq 1 ]; then
  N_DONE=$(find "$OUTDIR/.done" -maxdepth 1 -name "*.done" -type f | wc -l)
  N_TODO=$(( N_SAMPLES - N_DONE ))
  echo "Resume mode: $N_DONE already done, $N_TODO remaining"
  echo "$N_DONE" > "$PROGRESS_FILE"
else
  echo "0" > "$PROGRESS_FILE"
fi


########################################
# QC
########################################
if [ "$SKIP_QC" -eq 1 ]; then
  echo "[1/6] QC... SKIPPED"

  QC_FILE="$OUTDIR/qc_report.tsv"
  printf "ID\tgenome_size_bp\tn_contigs\tgc_pct\tqc_status\tqc_warnings\n" > "$QC_FILE"
  while IFS= read -r sample; do
    printf "%s\t-\t-\t-\t-\t-\n" "$sample"
  done < "$SAMPLES" >> "$QC_FILE"

else

echo "[1/6] QC... RUNNING"

########################################
# QC — assembly metrics
# -----------------------------------------------------------------------
# Checks each assembly for:
#   genome_size_bp  — total nucleotide length (excl. headers)
#   n_contigs       — number of sequences in FASTA
#   gc_pct          — GC content (%)
#
# Reference values for C. jejuni / C. coli:
#   GC content  : ~30.5% ± 3%  → warn outside 27.5–33.5%
#   Genome size : 1.3–2.0 Mb   → warn outside range
#   Contigs     : >500          → warn (fragmented assembly)
#
# Output: results/qc_report.tsv
#   ID  genome_size_bp  n_contigs  gc_pct  qc_status  qc_warnings
#
# qc_status: PASS | WARN
# qc_warnings: pipe-separated list of triggered warnings, e.g.
#   HIGH_CONTIGS|SMALL_GENOME
########################################

echo "  computing assembly metrics (genome size / GC% / contigs)..."

QC_FILE="$OUTDIR/qc_report.tsv"

awk '
BEGIN {
    OFS="\t"
    print "ID","genome_size_bp","n_contigs","gc_pct","qc_status","qc_warnings"

    # thresholds
    GC_MIN  = 27.5
    GC_MAX  = 33.5
    SZ_MIN  = 1300000
    SZ_MAX  = 2000000
    CTG_MAX = 500
}

FNR==1 {
    # new file — save previous sample if any
    if (sample != "") {
        gc_pct = (total > 0) ? (gc / total * 100) : 0

        warnings = ""
        if (gc_pct < GC_MIN) warnings = warnings (warnings?"|":"") "LOW_GC"
        if (gc_pct > GC_MAX) warnings = warnings (warnings?"|":"") "HIGH_GC"
        if (total < SZ_MIN)  warnings = warnings (warnings?"|":"") "SMALL_GENOME"
        if (total > SZ_MAX)  warnings = warnings (warnings?"|":"") "LARGE_GENOME"
        if (contigs > CTG_MAX) warnings = warnings (warnings?"|":"") "HIGH_CONTIGS"

        status = (warnings == "") ? "PASS" : "WARN"
        printf "%s\t%d\t%d\t%.1f\t%s\t%s\n", \
            sample, total, contigs, gc_pct, status, (warnings==""?"-":warnings)
    }

    # start new sample
    n = split(FILENAME, a, "/")
    sample = a[n]
    sub(/\.[^.]+$/, "", sample)
    total   = 0
    gc      = 0
    contigs = 0
}

/^>/ { contigs++; next }

{
    len = length($0)
    total += len
    for (i=1; i<=len; i++) {
        c = substr($0,i,1)
        if (c=="G"||c=="C"||c=="g"||c=="c") gc++
    }
}

END {
    if (sample != "") {
        gc_pct = (total > 0) ? (gc / total * 100) : 0

        warnings = ""
        if (gc_pct < GC_MIN) warnings = warnings (warnings?"|":"") "LOW_GC"
        if (gc_pct > GC_MAX) warnings = warnings (warnings?"|":"") "HIGH_GC"
        if (total < SZ_MIN)  warnings = warnings (warnings?"|":"") "SMALL_GENOME"
        if (total > SZ_MAX)  warnings = warnings (warnings?"|":"") "LARGE_GENOME"
        if (contigs > CTG_MAX) warnings = warnings (warnings?"|":"") "HIGH_CONTIGS"

        status = (warnings == "") ? "PASS" : "WARN"
        printf "%s\t%d\t%d\t%.1f\t%s\t%s\n", \
            sample, total, contigs, gc_pct, status, (warnings==""?"-":warnings)
    }
}
' "$INPUT"/*.fa "$INPUT"/*.fasta "$INPUT"/*.fna 2>/dev/null > "$QC_FILE"

# Summary
N_WARN=$(awk 'NR>1 && $5=="WARN"' "$QC_FILE" | wc -l)
echo "  QC complete: $N_SAMPLES samples — $(( N_SAMPLES - N_WARN )) PASS, $N_WARN WARN"
echo "  QC report → $QC_FILE"
echo ""

fi

########################################
# MLST TYPING — 7-gene allele-based (pubMLST scheme)
# Genes: aspA, glnA, gltA, glyA, pgm, tkt, uncA
# Each gene BLASTed separately → allele number from header (>gene_N)
# Allele profile → ST + clonal complex via mlst.txt lookup table
# Thresholds: ≥99% identity, ≥90% coverage per gene
# Output: mlst.tsv — ID, MLST_ST, MLST_CC
########################################
if [ "$SKIP_MLST" -eq 1 ]; then
  echo "[2/6] MLST... SKIPPED"

  MLST_OUT="$OUTDIR/mlst.tsv"
  printf "ID\tMLST_ST\tMLST_CC\n" > "$MLST_OUT"
  while IFS= read -r sample; do
    printf "%s\t-\t-\n" "$sample"
  done < "$SAMPLES" >> "$MLST_OUT"

else


MLST_GENE_DIR="databases/mlst"
MLST_PROFILE="databases/mlst/mlst.txt"
MLST_OUT="$OUTDIR/mlst.tsv"
MLST_RAW="$OUTDIR/mlst_raw"
MLST_PID=99
MLST_GENE_COV=0.90
MLST_GENES="aspA glnA gltA glyA pgm tkt uncA"

echo "[2/6] MLST typing... RUNNING"
mkdir -p "$MLST_RAW"

for f in "$INPUT"/*.{fa,fasta,fna}; do
  [ -f "$f" ] || continue
  n=$(basename "$f"); n=${n%.*}
  out_mlst="$MLST_RAW/${n}_mlst.txt"
  [ "$RESUME" -eq 1 ] && [ -s "$out_mlst" ] && continue

  limit_jobs "$JOBS"
  (
    alleles=""
    for gene in $MLST_GENES; do
      allele=$(blastn \
        -query "$f" \
        -db "$MLST_GENE_DIR/${gene}" \
        -outfmt "6 sseqid pident length slen" \
        -max_target_seqs 1 \
        -num_threads "$THREADS_PER_JOB" \
        -perc_identity "$MLST_PID" \
        2>/dev/null \
        | awk -v cov="$MLST_GENE_COV" '
          BEGIN{best=0; al="-"}
          {
            c=$3/$4; score=$2*c
            if(c>=cov && score>best){best=score; gsub(/.*_/,"",$1); al=$1}
          }
          END{print al}')
      alleles="$alleles $allele"
    done
    echo "$alleles" > "$out_mlst"
  ) &
done
wait

# Collect alleles and look up ST + CC from profile table
printf "ID\tMLST_ST\tMLST_CC\n" > "$MLST_OUT"
for f in "$INPUT"/*.{fa,fasta,fna}; do
  [ -f "$f" ] || continue
  n=$(basename "$f"); n=${n%.*}
  out_mlst="$MLST_RAW/${n}_mlst.txt"

  st="-"; cc="-"
  if [ -f "$out_mlst" ]; then
    read -r a1 a2 a3 a4 a5 a6 a7 < "$out_mlst"
    # only look up if all 7 alleles found
    if [ "$a1" != "-" ] && [ "$a2" != "-" ] && [ "$a3" != "-" ] && \
       [ "$a4" != "-" ] && [ "$a5" != "-" ] && [ "$a6" != "-" ] && [ "$a7" != "-" ]; then
      result=$(awk -v a1="$a1" -v a2="$a2" -v a3="$a3" -v a4="$a4" \
                   -v a5="$a5" -v a6="$a6" -v a7="$a7" \
        'NR>1 && $2==a1 && $3==a2 && $4==a3 && $5==a4 && $6==a5 && $7==a6 && $8==a7 \
         {print $1"\t"$9; exit}' "$MLST_PROFILE")
      if [ -n "$result" ]; then
        st=$(echo "$result" | cut -f1)
        cc=$(echo "$result" | cut -f2)
        [ -z "$cc" ] && cc="-"
      else
        st="novel"; cc="-"
      fi
    fi
  fi
  printf "%s\t%s\t%s\n" "$n" "$st" "$cc"
done >> "$MLST_OUT"

N_MLST=$(awk 'NR>1 && $2!="-" && $2!="novel"' "$MLST_OUT" | wc -l)
N_NOVEL=$(awk 'NR>1 && $2=="novel"' "$MLST_OUT" | wc -l)
echo "  MLST complete: $N_MLST/$N_SAMPLES typed, $N_NOVEL novel → $MLST_OUT"
echo ""

fi

########################################
# HS TYPING — TWO-STEP SYSTEM
#
# STEP 1: Marker gene blastn → precise call for known types
#   — Three tiers: A=99.5%/80cov, B=95%/75cov, C=90%/60cov
#   — Most marker genes wins; tie → higher identity wins
#   — If hit → final call, full locus never runs for this sample
#   — If no hit → go to Step 2
#
# STEP 2: Full-locus blastn → fallback for unmarked types
#   — Only runs for samples with no marker hit
#   — Covers: HS9, HS11, HS12, HS18, HS21, HS22, HS27, HS29,
#             HS31, HS32, HS33_35, HS37, HS38, HS40, HS52,
#             HS55, HS57, HS60, HS63, HS5_31 etc.
#
# STEP 3: Bracket subtypes
#   — HS23_36 → poly-G/C tract count → (HS23*) or (HS36*)
#   — HS8_17  → AA blastx → (HS8*) or (HS17*)
########################################

if [ "$SKIP_HS" -eq 1 ]; then

  echo "[3/6] HS (Penner) typing... SKIPPED (--skip_hs)"

  # Write stub hs_type.tsv — all samples get "-"
  printf "ID\tHS_TYPE\n" > "$OUTDIR/hs_type.tsv"
  while IFS= read -r sample; do
    printf "%s\t-\n" "$sample"
  done < "$SAMPLES" >> "$OUTDIR/hs_type.tsv"
  echo "  HS_TYPE set to '-' for all $N_SAMPLES samples → $OUTDIR/hs_type.tsv"

else

echo "[3/6] HS (Penner) typing... RUNNING"

mkdir -p "$OUTDIR/hs_blast_raw"

HS_MARKER_DB="$HS_MARKER"
HS_FULLLOCUS_DB="$HS_DB/hs_db_new"

########################################
# STEP 1 — Marker gene BLAST (all samples)
########################################

echo "  [Step 1] Marker gene BLAST (all samples)..."

for f in "$INPUT"/*.{fa,fasta,fna}; do
  [ -f "$f" ] || continue
  n=$(basename "$f"); n=${n%.*}
  out_marker="$OUTDIR/hs_blast_raw/${n}_marker.blast"
  [ "$RESUME" -eq 1 ] && [ -s "$out_marker" ] && continue
  limit_jobs "$JOBS"
  (
    run_blastn -query "$f" -db "$HS_MARKER_DB" \
      -outfmt "6 sseqid pident length slen bitscore" \
      -evalue 1e-10 -num_threads "$THREADS_PER_JOB" \
      > "$out_marker"
  ) &
done
wait

########################################
# STEP 1 CALLING — marker gene
########################################

awk -v SAMPLES="$SAMPLES" '
BEGIN{
  FS=OFS="\t"
  while((getline line < SAMPLES)>0) samples[line]=1
  type_map["HS1"]    = "HS1c"
  type_map["HS2"]    = "HS2"
  type_map["HS3"]    = "HS3"
  type_map["HS4"]    = "HS4c"
  type_map["HS6/7"]  = "HS6_7c"
  type_map["HS8/17"] = "HS8_17c"
  type_map["HS10"]   = "HS10"
  type_map["HS15"]   = "HS15c"
  type_map["HS19"]   = "HS19"
  type_map["HS23/36"]= "HS23_36c"
  type_map["HS41"]   = "HS41"
  type_map["HS42"]   = "HS42"
  type_map["HS53"]   = "HS53"
  tid_A=99.5; tcov_A=0.80
  tid_B=95.0; tcov_B=0.75
  tid_C=90.0; tcov_C=0.60
}
function get_tier(pid, cov){
  if(pid >= tid_A && cov >= tcov_A) return 1
  if(pid >= tid_B && cov >= tcov_B) return 2
  if(pid >= tid_C && cov >= tcov_C) return 3
  return 0
}
{
  split(FILENAME,a,"/"); file=a[length(a)]
  sub(/_marker\.blast$/,"",file); sample=file
  if(!(sample in samples)) next
  sseqid=$1; pid=$2+0; cov=($4>0?$3/$4:0)
  idx=index(sseqid,"_")
  if(idx==0) next
  hstype=substr(sseqid,1,idx-1)
  if(!(hstype in type_map)) next
  call=type_map[hstype]
  t=get_tier(pid,cov)
  if(t==0) next
  # HS19 requires Tier A only (99.5% min) — cj1419 cross-reacts at 97-98%
  if(call=="HS19" && t>1) next
  gene_key=sample SUBSEP t SUBSEP call SUBSEP sseqid
  if(!(gene_key in seen)){
    seen[gene_key]=1
    gene_count[sample,t,call]++
    sum_pid[sample,t,call] += pid
    if(pid > best_pid[sample,t,call]) best_pid[sample,t,call]=pid
  }
}
END{
  print "ID","MARKER_TYPE","MARKER_TIER"
  for(s in samples){
    call="-"; best_tier=0
    for(t=1;t<=3;t++){
      best_call="-"; best_count=0; best_avg=0
      second_call="-"; second_count=0; second_avg=0
      for(c in type_map){
        ct=type_map[c]
        cnt=gene_count[s,t,ct]+0
        if(cnt==0) continue
        avg=(cnt>0 ? sum_pid[s,t,ct]/cnt : 0)
        if(cnt > best_count || (cnt==best_count && avg > best_avg)){
          second_call=best_call; second_count=best_count; second_avg=best_avg
          best_call=ct; best_count=cnt; best_avg=avg
        } else if(cnt > second_count || (cnt==second_count && avg > second_avg)){
          second_call=ct; second_count=cnt; second_avg=avg
        }
      }
      if(best_call != "-"){
        if(second_call != "-" && second_count == best_count){
          ba=int(best_avg*10+0.5)/10; sa=int(second_avg*10+0.5)/10
          if(ba == sa){
            call=(best_call < second_call ? best_call"_"second_call : second_call"_"best_call)
          } else { call=best_call }
        } else { call=best_call }
        best_tier=t; break
      }
    }
    tier_name=(best_tier==1?"A":(best_tier==2?"B":(best_tier==3?"C":"NONE")))
    print s, call, tier_name
  }
}
' "$OUTDIR"/hs_blast_raw/*_marker.blast \
  > "$OUTDIR/hs_marker_calls.tsv"

N_MARKER=$(awk 'NR>1 && $3!="NONE"' "$OUTDIR/hs_marker_calls.tsv" | wc -l)
N_FALLBACK=$(awk 'NR>1 && $3=="NONE"' "$OUTDIR/hs_marker_calls.tsv" | wc -l)
echo "  Marker calls: $N_MARKER | Needs full-locus: $N_FALLBACK"

# Write list of samples needing full-locus
FULLLOCUS_SAMPLES="$OUTDIR/fulllocus_samples.txt"
awk 'NR>1 && $3=="NONE" {print $1}' "$OUTDIR/hs_marker_calls.tsv" > "$FULLLOCUS_SAMPLES"

########################################
# STEP 2 — Full-locus BLAST (ALL samples)
# Needed for conflict resolution and fallback
########################################

echo "  [Step 2] Full-locus BLAST (all samples)..."

for f in "$INPUT"/*.{fa,fasta,fna}; do
  [ -f "$f" ] || continue
  n=$(basename "$f"); n=${n%.*}
  out_full="$OUTDIR/hs_blast_raw/${n}_fulllocus.blast"
  [ "$RESUME" -eq 1 ] && [ -s "$out_full" ] && continue
  limit_jobs "$JOBS"
  (
    run_blastn -query "$f" -db "$HS_FULLLOCUS_DB" \
      -outfmt "6 sseqid pident length slen bitscore" \
      -evalue 1e-10 -num_threads "$THREADS_PER_JOB" \
    | awk '$2 >= 78 {
        sseqid=$1; slen=$4; aln=$3; bs=$5; pident=$2
        total_aln[sseqid] += aln
        if(slen > slen_map[sseqid]) slen_map[sseqid]=slen
        if(bs > best_bs[sseqid]){
          best_bs[sseqid]  = bs
          best_pid[sseqid] = pident
        }
      }
      END {
        for(s in total_aln)
          printf "%s\t%.3f\t%d\t%d\t%.0f\n",
            s, best_pid[s], total_aln[s], slen_map[s], best_bs[s]
      }' > "$out_full"
  ) &
done
wait

# Multi-gene BLAST for HS2/HS6/HS53 (all samples)
for mg in HS2 HS6 HS53; do
  mg_db="$HS_DB/$mg/$mg"
  [ -f "${mg_db}.nin" ] || continue
  for f in "$INPUT"/*.{fa,fasta,fna}; do
    [ -f "$f" ] || continue
    n=$(basename "$f"); n=${n%.*}
    out_mg="$OUTDIR/hs_blast_raw/${n}_${mg}.blast"
    [ "$RESUME" -eq 1 ] && [ -s "$out_mg" ] && continue
    limit_jobs "$JOBS"
    (
      run_blastn -query "$f" -db "$mg_db" \
        -outfmt "6 sseqid pident length slen bitscore" \
        -evalue 1e-10 -num_threads "$THREADS_PER_JOB" \
      | awk '{ cov=($4>0?$3/$4:0); if($2>=78 && cov>=0.80) print }' \
        > "$out_mg"
    ) &
  done
  wait
done

########################################
# STEP 2 CALLING — full-locus best-hit
########################################

awk -v SAMPLES="$SAMPLES" '
BEGIN{
  FS=OFS="\t"
  while((getline line < SAMPLES)>0) samples[line]=1
  req["HS2"]=15; req["HS6"]=8; req["HS53"]=7
  split("HS2 HS6 HS53", tmp, " ")
  for(i in tmp) multigene[tmp[i]]=1
  complex_norm["HS23"]    = "HS23_36c"
  complex_norm["HS36"]    = "HS23_36c"
  complex_norm["HS23_36"] = "HS23_36c"
  complex_norm["HS6"]     = "HS6_7c"
  complex_norm["HS7"]     = "HS6_7c"
  complex_norm["HS6_7"]   = "HS6_7c"
  complex_norm["HS8"]     = "HS8_17c"
  complex_norm["HS17"]    = "HS8_17c"
  complex_norm["HS8_17"]  = "HS8_17c"
  complex_norm["HS5"]     = "HS5_31c"
  complex_norm["HS31"]    = "HS5_31c"
  complex_norm["HS5_31"]  = "HS5_31c"
  complex_norm["HS4"]     = "HS4c"
  complex_norm["HS4c"]    = "HS4c"
  complex_norm["HS13"]    = "HS4c"
  complex_norm["HS16"]    = "HS4c"
  complex_norm["HS43"]    = "HS4c"
  complex_norm["HS50"]    = "HS4c"
  complex_norm["HS62"]    = "HS4c"
  complex_norm["HS63"]    = "HS4c"
  complex_norm["HS64"]    = "HS4c"
  complex_norm["HS65"]    = "HS4c"
  complex_norm["HS4_2"]   = "HS4c"
  complex_norm["HS1"]     = "HS1c"
  complex_norm["HS44"]    = "HS1c"
  complex_norm["HS15"]    = "HS15c"
  complex_norm["HS58"]    = "HS15c"
  complex_norm["HS45"]    = "HS45c"
  complex_norm["HS32"]    = "HS45c"
  complex_norm["HS60"]    = "HS45c"
  complex_norm["HS33_2"]  = "HS33_35"
  complex_norm["HS55_2"]  = "HS55"
  no_mosaic["HS4_HS31"]=1;     no_mosaic["HS31_HS4"]=1
  no_mosaic["HS12_HS55"]=1;    no_mosaic["HS55_HS12"]=1
  no_mosaic["HS41_HS55"]=1;    no_mosaic["HS55_HS41"]=1
  no_mosaic["HS42_HS55"]=1;    no_mosaic["HS55_HS42"]=1
  no_mosaic["HS23_36_HS5_31"]=1; no_mosaic["HS5_31_HS23_36"]=1
  no_mosaic["HS3_HS4"]=1;      no_mosaic["HS4_HS3"]=1
  no_mosaic["HS12_HS41"]=1;    no_mosaic["HS41_HS12"]=1
  no_mosaic["HS12_HS42"]=1;    no_mosaic["HS42_HS12"]=1
  no_mosaic["HS41_HS42"]=1;    no_mosaic["HS42_HS41"]=1
  no_mosaic["HS1_HS19"]=1;     no_mosaic["HS19_HS1"]=1
  no_mosaic["HS1_HS37"]=1;     no_mosaic["HS37_HS1"]=1
}
function norm(hs){ return (hs in complex_norm ? complex_norm[hs] : hs) }
function tier_from_cov_pid(pid, cov, is_mg){
  if(!is_mg){
    if(pid >= 99.5 && cov >= 0.90) return 1
    if(pid >= 90   && cov >= 0.60) return 2
    if(pid >= 80   && cov >= 0.40) return 3
  } else {
    if(pid >= 99.5) return 1
    if(pid >= 95)   return 2
    if(pid >= 80)   return 3
  }
  return 0
}
function gene_min(hs,tier){
  if(tier==1) return req[hs]
  if(tier==2) return int(req[hs]*0.8+0.999)
  if(tier==3) return int(req[hs]*0.6+0.999)
  return req[hs]
}
{
  split(FILENAME,a,"/"); file=a[length(a)]
  if(file ~ /_fulllocus\.blast$/){
    sub(/_fulllocus\.blast$/,"",file); sample=file
    if(!(sample in samples)) next
    hs=norm($1); pid=$2+0; cov=($4>0?$3/$4:0); score=pid*cov*$5
    raw_hs=$1
    tier=tier_from_cov_pid(pid,cov,0)
    if(tier==0) next
    if($3>=3000){
      # Track raw bitscore per type for HS3/HS29 conflict resolution
      bs=$5+0
      if(bs > raw_bs[sample,$1]) raw_bs[sample,$1]=bs
      if(score > full_score[sample,tier]){
        full_second_score[sample,tier] = full_score[sample,tier]
        full_second_call[sample,tier]  = full_call[sample,tier]
        full_second_cov[sample,tier]   = full_cov[sample,tier]
        full_second_pid[sample,tier]   = full_pid[sample,tier]
        full_score[sample,tier]        = score
        full_call[sample,tier]         = hs
        full_raw_call[sample,tier]     = raw_hs
        full_cov[sample,tier]          = cov
        full_pid[sample,tier]          = pid
        full_cov_hs[sample,tier,hs]    = cov
      } else if(score > full_second_score[sample,tier] && hs != full_call[sample,tier]){
        full_second_score[sample,tier] = score
        full_second_call[sample,tier]  = hs
        full_second_cov[sample,tier]   = cov
        full_second_pid[sample,tier]   = pid
        full_cov_hs[sample,tier,hs]    = cov
      } else {
        full_cov_hs[sample,tier,hs]    = cov
      }
    }
  } else {
    split(file,b,"_HS"); sample=b[1]; hs="HS"b[2]; sub(/\.blast$/,"",hs)
    if(!(sample in samples)) next
    pid=$2+0; cov=($4>0?$3/$4:0); score=pid*cov*$5
    tier=tier_from_cov_pid(pid,cov,1)
    if(tier==0) next
    if(hs in multigene){
      gene=$1; key=sample SUBSEP tier SUBSEP hs SUBSEP gene
      if(!(key in seen_gene)){ seen_gene[key]=1; gene_count[sample,tier,hs]++ }
      if(score > gene_score[sample,tier,hs]) gene_score[sample,tier,hs]=score
      if(gene=="rfbE") seen_marker[sample,tier,"rfbE"]=1
      if(gene=="siaA") seen_marker[sample,tier,"siaA"]=1
      if(gene=="ispD") seen_marker[sample,tier,"ispD"]=1
    }
  }
}
END{
  print "ID","HS_TYPE","HS_BLAST_WINNER"
  for(s in samples){
    call="-"; raw_winner="-"
    for(t=1;t<=3;t++){
      pcov = (t==1 ? 0.90 : (t==2 ? 0.60 : 0.40))
      has_full = ((s,t) in full_call)
      if(has_full){
        b_hs = full_call[s,t]; b_cov = full_cov[s,t]; b_pid = full_pid[s,t]
        s_hs = full_second_call[s,t]; s_cov = full_second_cov[s,t]+0; s_pid = full_second_pid[s,t]+0
        # Store raw winner BEFORE normalisation
        raw_winner = full_raw_call[s,t]
        if(raw_winner == "") raw_winner = b_hs
        if(s_hs != "" && s_hs != "-" && s_cov >= pcov && t <= 2){
          pid_diff = b_pid - s_pid; if(pid_diff<0) pid_diff=-pid_diff
          cov_diff = b_cov - s_cov; if(cov_diff<0) cov_diff=-cov_diff
          if(b_hs == s_hs){ call=b_hs; break }
          pair_key = (b_hs < s_hs ? b_hs"_"s_hs : s_hs"_"b_hs)
          if(pid_diff <= 0.5 && cov_diff <= 0.10 && !(pair_key in no_mosaic)){
            call = (b_hs < s_hs ? b_hs"_"s_hs : s_hs"_"b_hs)
          } else { call=b_hs }
        } else { call=b_hs }
        if(t <= 2){
          for(mg in req){
            best_gene_count=0; marker_ok=0
            for(tt=1;tt<=3;tt++){
              if(gene_count[s,tt,mg]+0 > best_gene_count)
                best_gene_count=gene_count[s,tt,mg]+0
              if(mg=="HS2"  && seen_marker[s,tt,"rfbE"]) marker_ok=1
              if(mg=="HS6"  && seen_marker[s,tt,"siaA"]) marker_ok=1
              if(mg=="HS53" && seen_marker[s,tt,"ispD"]) marker_ok=1
            }
            min_genes=int(req[mg]*0.6+0.999)
            if(best_gene_count >= min_genes && marker_ok)
              call=(mg=="HS6" ? "HS6_7" : mg)
          }
        }
        break
      }
      best_gene="-"; best_ratio=0; best_gscore=0
      for(hs in req){
        c=gene_count[s,t,hs]+0; minc=gene_min(hs,t)
        if(c>=minc){
          ratio=c/req[hs]; gs=gene_score[s,t,hs]+0
          if(ratio>best_ratio||(ratio==best_ratio&&gs>best_gscore)){
            best_ratio=ratio; best_gscore=gs; best_gene=hs
          }
        }
      }
      if(best_gene!="-"){ call=best_gene; raw_winner=best_gene; break }
    }
    # HS3 vs HS29: if HS29 bitscore > HS3 bitscore by more than 10% → call HS29
    if(call == "HS3"){
      bs3  = raw_bs[s,"HS3"]+0
      bs29 = raw_bs[s,"HS29"]+0
      if(bs29 > 0 && bs3 > 0 && bs29 > bs3 * 1.10){
        call = "HS29"; raw_winner = "HS29"
      }
    }
    # Minimum threshold check — if best hit < 90% identity or < 60% coverage → untypable
    if(call != "-"){
      best_pid_val = full_pid[s,1]+0
      best_cov_val = full_cov[s,1]+0
      if(best_pid_val > 0 && (best_pid_val < 90.0 || best_cov_val < 0.60)){
        call = "-"; raw_winner = "-"
      }
    }
    print s, call, raw_winner
  }
}
' "$OUTDIR"/hs_blast_raw/*_fulllocus.blast \
  "$OUTDIR"/hs_blast_raw/*_HS2.blast \
  "$OUTDIR"/hs_blast_raw/*_HS6.blast \
  "$OUTDIR"/hs_blast_raw/*_HS53.blast \
  2>/dev/null \
  > "$OUTDIR/hs_fulllocus_calls.tsv"

echo "  Full-locus calls done → $OUTDIR/hs_fulllocus_calls.tsv"

########################################
# MERGE Step 1 + Step 2 → hs_type.tsv
# Marker call wins if tier != NONE
# Otherwise keep full-locus call
########################################

python3 - "$OUTDIR/hs_marker_calls.tsv" \
           "$OUTDIR/hs_fulllocus_calls.tsv" \
           "$OUTDIR/hs_type.tsv" <<'PYEOF'
import sys, csv

marker_file, fulllocus_file, out_file = sys.argv[1:4]

COMPLEX_NORM = {
    "HS23": "HS23_36c", "HS36": "HS23_36c", "HS23_36": "HS23_36c", "HS23_36c": "HS23_36c",
    "HS6":  "HS6_7c",   "HS7":  "HS6_7c",   "HS6_7":   "HS6_7c",   "HS6_7c":   "HS6_7c",
    "HS8":  "HS8_17c",  "HS17": "HS8_17c",  "HS8_17":  "HS8_17c",  "HS8_17c":  "HS8_17c",
    "HS5":  "HS5_31c",  "HS31": "HS5_31c",  "HS5_31":  "HS5_31c",  "HS5_31c":  "HS5_31c",
    "HS4":  "HS4c",     "HS4c": "HS4c",     "HS13": "HS4c", "HS16": "HS4c",
    "HS43": "HS4c",     "HS50": "HS4c",     "HS62": "HS4c", "HS63": "HS4c",
    "HS64": "HS4c",     "HS65": "HS4c",     "HS4_2": "HS4c",
    "HS1":  "HS1c",     "HS44": "HS1c",     "HS1c": "HS1c",
    "HS15": "HS15c",    "HS58": "HS15c",    "HS15c": "HS15c",
    "HS45": "HS45c",    "HS32": "HS45c",    "HS60": "HS45c", "HS45c": "HS45c",
    "HS33_2": "HS33_35",
    "HS55_2": "HS55",
}

def norm(call):
    if call in COMPLEX_NORM:
        return COMPLEX_NORM[call]
    parts = call.split("_HS")
    if len(parts) == 1:
        return call
    components = [parts[0]] + ["HS" + p for p in parts[1:]]
    seen = []
    for p in components:
        n = COMPLEX_NORM.get(p, p)
        if n not in seen:
            seen.append(n)
    return "_".join(seen)

# Load marker calls
marker = {}
with open(marker_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        marker[row["ID"]] = (row["MARKER_TYPE"], row["MARKER_TIER"])

# Load full-locus calls — now has 3 columns: ID, HS_TYPE, HS_BLAST_WINNER
fulllocus = {}
fulllocus_winner = {}
with open(fulllocus_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        fulllocus[row["ID"]] = norm(row["HS_TYPE"])
        fulllocus_winner[row["ID"]] = row.get("HS_BLAST_WINNER", "-")

# Complex members in DB — used to add blast-based bracket
# Only add bracket if raw winner is a known member of that complex
COMPLEX_MEMBERS = {
    "HS4c":     ["HS4", "HS63", "HS4_2"],
    "HS6_7c":   ["HS6_7", "HS6", "HS7"],
    "HS8_17c":  ["HS8_17", "HS8", "HS17"],
    "HS5_31c":  ["HS5_31", "HS5", "HS31"],
    "HS23_36c": ["HS23_36", "HS23", "HS36"],
    "HS1c":     ["HS1", "HS44"],
    "HS15c":    ["HS15", "HS58"],
    "HS45c":    ["HS45", "HS32", "HS60"],
    "HS33_35":  ["HS33_35", "HS33_2"],
    "HS55":     ["HS55", "HS55_2"],
}

def add_blast_bracket(complex_call, raw_winner):
    """Add bracket from blast winner if it's a specific member of the complex.
    Single * = blast-predicted subtype within complex."""
    members = COMPLEX_MEMBERS.get(complex_call, [])
    if raw_winner in members and raw_winner not in [complex_call, "-", ""]:
        return f"{complex_call}({raw_winner}*)"
    return complex_call

# Merge — marker wins EXCEPT for known cross-reactive conflict pairs
with open(out_file, "w") as fo:
    fo.write("ID\tHS_TYPE\tHS_BLAST_WINNER\n")
    for sid, (m_call, m_tier) in marker.items():
        fl_call = fulllocus.get(sid, "-")
        fl_winner = fulllocus_winner.get(sid, "-")
        if m_tier != "NONE" and m_call != "-":
            conflict = False
            if m_call == "HS15c" and fl_call == "HS15c(HS58*)": conflict = True
            if m_call == "HS41"  and fl_call == "HS55": conflict = True
            if conflict:
                fo.write(f"{sid}\t{fl_call}\t{fl_winner}\n")
            else:
                final = add_blast_bracket(m_call, fl_winner)
                fo.write(f"{sid}\t{final}\t{fl_winner}\n")
        elif fl_call != "-":
            final = add_blast_bracket(fl_call, fl_winner)
            fo.write(f"{sid}\t{final}\t{fl_winner}\n")
        else:
            fo.write(f"{sid}\t-\t-\n")

PYEOF

echo "  HS typing complete → $OUTDIR/hs_type.tsv"

########################################
# STEP 3 — Bracket subtypes
#   HS8_17  → AA blastx → HS8_17(HS8*) or HS8_17(HS17*)
#   HS23_36 → poly-G/C tract count → HS23_36(HS23*) or HS23_36(HS36*)
#   All other complexes → no bracket, final as-is
########################################

echo "  [Step 3] Adding subtype brackets..."

HS_AA_DB_DIR="$HS_DB/../hs_blast_aa"
[ -d "$HS_AA_DB_DIR" ] || HS_AA_DB_DIR="databases/hs_blast_aa"

BLASTX_CMD=""
if [ -n "$PKG_MANAGER" ] && $PKG_MANAGER env list 2>/dev/null | grep -q "^${ENV_NAME}[[:space:]]"; then
  BLASTX_CMD="$PKG_MANAGER run -n $ENV_NAME blastx"
elif command -v blastx >/dev/null 2>&1; then
  BLASTX_CMD="blastx"
fi

SEQKIT_CMD=""
if [ -n "$PKG_MANAGER" ] && $PKG_MANAGER env list 2>/dev/null | grep -q "^${ENV_NAME}[[:space:]]"; then
  SEQKIT_CMD="$PKG_MANAGER run -n $ENV_NAME seqkit"
elif command -v seqkit >/dev/null 2>&1; then
  SEQKIT_CMD="seqkit"
fi

python3 - \
  "$OUTDIR/hs_type.tsv" \
  "$INPUT" \
  "$HS_AA_DB_DIR" \
  "${BLASTX_CMD:-}" \
  "$OUTDIR/hs_type.tsv" \
  "${SEQKIT_CMD:-}" \
  "${HS_MARKER_DB:-databases/hs_specific_aa_per_hs/hs/hs_aa}" \
<<'PYEOF'
import sys, os, subprocess, re

tsv_in, asm_dir, aa_db_dir, blastx_cmd, tsv_out, seqkit_cmd, marker_db = sys.argv[1:8]

########################################
# HS8_17 bracketing via AA blastx
########################################
def count_aa_matches(asm, db, blastx):
    if not blastx or not os.path.exists(f"{db}.phr"): return 0
    cmd = blastx.split() + [
        "-query", asm, "-db", db,
        "-outfmt", "6 sseqid sstart qseq",
        "-evalue", "1e-5", "-word_size", "2",
        "-matrix", "PAM30", "-seg", "no"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    matches = 0
    for line in result.stdout.strip().split("\n"):
        if not line: continue
        p = line.split("\t")
        if len(p) < 3: continue
        try: sstart = int(p[1])
        except: continue
        qseq = p[2]
        parts = p[0].split("|")
        if len(parts) < 4: continue
        expected = parts[3]
        idx = 11 - sstart
        if 0 <= idx < len(qseq) and qseq[idx] == expected:
            matches += 1
    return matches

########################################
# HS23_36c bracketing via poly-G/C tract
# Concatenates all contig hits in reference order
# then finds longest poly-G/C run
# divisible by 3 → HS23* (gene OFF/frameshifted)
# NOT divisible by 3 → HS36* (gene ON/in-frame)
########################################
def get_polyg_bracket(asm, hs_marker_db, seqkit):
    try:
        cmd = [
            "blastn", "-query", asm,
            "-db", hs_marker_db,
            "-outfmt", "6 sseqid pident length slen qseqid qstart qend sstart send",
            "-evalue", "1e-10", "-perc_identity", "90",
            "-num_alignments", "20"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if not result.stdout.strip(): return ""

        # Collect all hits sorted by reference position
        hits = []
        for line in result.stdout.strip().split("\n"):
            if not line: continue
            p = line.split("\t")
            if len(p) < 9: continue
            if not p[0].startswith("HS23"): continue
            try:
                aln_len = int(p[2])
                qstart  = int(p[5])
                qend    = int(p[6])
                sstart  = int(p[7])
                qseqid  = p[4]
            except:
                continue
            if aln_len < 50: continue
            hits.append((sstart, qseqid, qstart, qend))

        if not hits: return ""

        # Sort by reference position and concatenate
        hits.sort(key=lambda x: x[0])
        seqkit_parts = seqkit.split() if seqkit else ["seqkit"]
        combined = ""
        for sstart, qseqid, qstart, qend in hits:
            cmd2 = seqkit_parts + ["grep", "-p", qseqid, asm]
            r2 = subprocess.run(cmd2, capture_output=True, text=True)
            cmd3 = seqkit_parts + ["subseq", "--quiet", "-r", f"{qstart}:{qend}"]
            r3 = subprocess.run(cmd3, input=r2.stdout, capture_output=True, text=True)
            seq = "".join(l.strip() for l in r3.stdout.split("\n")
                         if not l.startswith(">")).upper()
            combined += seq

        if not combined: return ""

        all_runs = re.findall(r'G{5,}|C{5,}', combined)
        if not all_runs: return ""

        longest = max(all_runs, key=len)
        n = len(longest)

        # divisible by 3 → HS23* (gene OFF)
        # NOT divisible by 3 → HS36* (gene ON)
        return "HS23*" if n % 3 == 0 else "HS36*"

    except Exception:
        return ""

# Read current calls
rows = []
with open(tsv_in) as f:
    header = f.readline()
    for line in f:
        rows.append(line.rstrip().split("\t"))

hs8_17_done = 0
hs23_36_done = 0

with open(tsv_out, "w") as fo:
    fo.write(header)
    for row in rows:
        sample = row[0]
        call   = row[1] if len(row) > 1 else "-"

        asm = None
        for ext in ["fa","fasta","fna"]:
            p = os.path.join(asm_dir, f"{sample}.{ext}")
            if os.path.exists(p):
                asm = p; break

        if asm:
            # HS8_17c bracket via AA blastx
            if call == "HS8_17c" and blastx_cmd:
                scores = {}
                for sub in ["HS8","HS17"]:
                    db = os.path.join(aa_db_dir, sub, sub)
                    scores[sub] = count_aa_matches(asm, db, blastx_cmd)
                best = max(scores, key=scores.get)
                best_s = scores[best]
                other_s = min(scores.values())
                if best_s >= 3 and best_s > other_s:
                    row[1] = f"HS8_17c({best}*)"
                    print(f"  {sample}: HS8_17c → HS8_17c({best}*)")
                    hs8_17_done += 1

            # HS23_36c bracket via poly-G/C tract
            elif call == "HS23_36c":
                if marker_db:
                    bracket = get_polyg_bracket(asm, marker_db, seqkit_cmd)
                    if bracket:
                        row[1] = f"HS23_36c({bracket})"
                        print(f"  {sample}: HS23_36c → HS23_36c({bracket})")
                        hs23_36_done += 1

        fo.write("\t".join(row) + "\n")

print(f"  HS8_17c subtype brackets: {hs8_17_done}")
print(f"  HS23_36c subtype brackets: {hs23_36_done}")
PYEOF

echo "  HS typing with brackets complete → $OUTDIR/hs_type.tsv"

########################################
# STEP 4 — In silico PCR confirmation + subtype bracketing
# Uses ipcress (exonerate package)
# Primers from databases/primers.txt
#
# Confirmation rules (override blast call):
#   HS29_target (150-250bp)                          → HS29
#   HS55_target (300-400bp)                          → HS55
#   HS58_target (70-120bp) + HS15_target (300-400bp) → HS15c(HS58*)
#
# Bracket rules (add bracket to existing complex call):
#   HS4c + HS4B_target (625-680bp)   → HS4c(HS16/64*)
#   HS5_31c + HS5_31_target (830-890bp) + HS15_target (300-400bp) → HS5_31c(HS31*)
#   HS5_31c + HS5_31_target (830-890bp), no HS15 → HS5_31c(HS5*)
#
# All other primers → confirmation only (no call change, logged)
########################################

echo "  [Step 4] In silico PCR (HS29/HS55/HS58)..."

PRIMERS_FILE="databases/primers.txt"
IPCRESS_CMD=""
if [ -n "$PKG_MANAGER" ] && $PKG_MANAGER env list 2>/dev/null | grep -q "^${ENV_NAME}[[:space:]]"; then
  IPCRESS_CMD="$PKG_MANAGER run -n $ENV_NAME ipcress"
elif command -v ipcress >/dev/null 2>&1; then
  IPCRESS_CMD="ipcress"
fi

if [ -z "$IPCRESS_CMD" ]; then
  echo "  WARNING: ipcress not found — skipping HS29/HS55/HS58 in silico PCR"
elif [ ! -f "$PRIMERS_FILE" ]; then
  echo "  WARNING: primers file not found at $PRIMERS_FILE — skipping"
else
  python3 - \
    "$OUTDIR/hs_type.tsv" \
    "$INPUT" \
    "$PRIMERS_FILE" \
    "$IPCRESS_CMD" \
    "$OUTDIR/hs_type.tsv" \
  <<'PYEOF'
import sys, os, subprocess, tempfile

tsv_in, asm_dir, primers_file, ipcress_cmd, tsv_out = sys.argv[1:6]
pcr_log = tsv_out.replace("hs_type.tsv", "hs_pcr_corrections.tsv")

def run_ipcress(asm, primers_file, ipcress_cmd):
    cmd = ipcress_cmd.split() + [
        primers_file, asm,
        "--mismatch", "2",
        "--pretty", "false"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    hits = {}
    for line in result.stdout.strip().split("\n"):
        if not line or line.startswith("#"): continue
        parts = line.split()
        if len(parts) < 4: continue
        if parts[0] != "ipcress:": continue
        name = parts[2]
        try:
            size = int(parts[3])
        except:
            continue
        if name not in hits:
            hits[name] = []
        hits[name].append(size)
    return hits

def check_size(sizes, min_s, max_s):
    return any(min_s <= s <= max_s for s in sizes)

def get_size(sizes, min_s, max_s):
    return next((s for s in sizes if min_s <= s <= max_s), None)

# Read current calls
rows = []
with open(tsv_in) as f:
    header = f.readline()
    for line in f:
        rows.append(line.rstrip().split("\t"))

pcr_done = 0
pcr_records = []

with open(tsv_out, "w") as fo:
    fo.write(header)
    for row in rows:
        sample = row[0]
        call   = row[1] if len(row) > 1 else "-"

        asm = None
        for ext in ["fa","fasta","fna"]:
            p = os.path.join(asm_dir, f"{sample}.{ext}")
            if os.path.exists(p):
                asm = p; break

        if asm:
            hits = run_ipcress(asm, primers_file, ipcress_cmd)
            new_call = None
            pcr_size = None

            # ── Override rules ──────────────────────────────────────────
            # HS29: confirmed by HS29_target 150-250bp
            if "HS29_target" in hits and check_size(hits["HS29_target"], 150, 250):
                new_call = "HS29"
                pcr_size = get_size(hits["HS29_target"], 150, 250)

            # HS55: confirmed by HS55_target 300-400bp
            elif "HS55_target" in hits and check_size(hits["HS55_target"], 300, 400):
                new_call = "HS55"
                pcr_size = get_size(hits["HS55_target"], 300, 400)

            # HS15c(HS58*): HS58_target 70-120bp AND HS15_target 300-400bp
            elif ("HS58_target" in hits and check_size(hits["HS58_target"], 70, 120) and
                  "HS15_target" in hits and check_size(hits["HS15_target"], 300, 400)):
                new_call = "HS15c(HS58*)"
                pcr_size = get_size(hits["HS58_target"], 70, 120)

            # ── Bracket rules (add bracket to complex call) ──────────────
            # HS4c + HS4B → HS4c(HS16/64*)
            elif call == "HS4c" and "HS4B_target" in hits and check_size(hits["HS4B_target"], 625, 680):
                new_call = "HS4c(HS16/64*)"
                pcr_size = get_size(hits["HS4B_target"], 625, 680)

            # HS4c + HS4A (no HS4B, no blast bracket already) → HS4c(HS4A**)
            # ** = PCR confirms HS4 complex only, no specific subtype from blast
            elif (call == "HS4c" and
                  "HS4A_target" in hits and check_size(hits["HS4A_target"], 340, 400) and
                  not ("HS4B_target" in hits and check_size(hits["HS4B_target"], 625, 680)) and
                  "(" not in call):  # only if no bracket already
                new_call = "HS4c(HS4A**)"
                pcr_size = get_size(hits["HS4A_target"], 340, 400)

            # HS5_31c + HS5_31_target + HS15_target → HS5_31c(HS31*)
            elif (call == "HS5_31c" and
                  "HS5_31_target" in hits and check_size(hits["HS5_31_target"], 830, 890) and
                  "HS15_target" in hits and check_size(hits["HS15_target"], 300, 400)):
                new_call = "HS5_31c(HS31*)"
                pcr_size = get_size(hits["HS5_31_target"], 830, 890)

            # HS5_31c + HS5_31_target only → HS5_31c(HS5*)
            elif (call == "HS5_31c" and
                  "HS5_31_target" in hits and check_size(hits["HS5_31_target"], 830, 890) and
                  not ("HS15_target" in hits and check_size(hits["HS15_target"], 300, 400))):
                new_call = "HS5_31c(HS5*)"
                pcr_size = get_size(hits["HS5_31_target"], 830, 890)

            if new_call and new_call != call:
                print(f"  {sample}: {call} → {new_call} ({pcr_size}bp)")
                pcr_records.append((sample, call, new_call, pcr_size))
                row[1] = new_call
                pcr_done += 1

        fo.write("\t".join(row[:2]) + "\n")

# Write PCR log
pcr_log = tsv_out.replace("hs_type.tsv", "hs_pcr_corrections.tsv")
with open(pcr_log, "w") as fl:
    fl.write("ID\tPREV_TYPE\tNEW_TYPE\tPRODUCT_SIZE_BP\n")
    for rec in pcr_records:
        fl.write("\t".join(str(x) for x in rec) + "\n")

print(f"  In silico PCR corrections: {pcr_done} → {pcr_log}")
PYEOF
fi

fi  # end --skip_hs check

########################################
# LOS BLAST
########################################
if [ "$SKIP_LOS" -eq 1 ]; then
  echo "[4/6] LOS class typing... SKIPPED"

  printf "ID\tLOS_CLASS\tGBS_RISK\tSIALYLATED\tLOS_TIER\tLOS_COVERAGE\tLOS_IDENTITY\n" > "$OUTDIR/los_final.tsv"
  while IFS= read -r sample; do
    printf "%s\t-\t-\t-\t-\t-\t-\n" "$sample"
  done < "$SAMPLES" >> "$OUTDIR/los_final.tsv"

else

echo "[4/6] LOS class typing... RUNNING"

mkdir -p "$OUTDIR/los_blast_raw"
LOS_COMBINED="$LOS_DB/los"

# One BLAST per sample — save ALL hits at minimum identity 78%
# AWK calling logic applies strict/mid/loose thresholds afterwards
for f in "$INPUT"/*.{fa,fasta,fna}; do
  [ -f "$f" ] || continue
  n=$(basename "$f"); n=${n%.*}
  out_los="$OUTDIR/los_blast_raw/${n}_LOS.blast"
  [ "$RESUME" -eq 1 ] && [ -s "$out_los" ] && continue
  limit_jobs "$JOBS"
  (
    clean_query "$f" \
    | run_blastn -db "$LOS_COMBINED" \
      -outfmt "6 sseqid pident length slen bitscore" \
      -evalue 1e-10 -num_threads "$THREADS_PER_JOB" 2>/dev/null \
    | awk '$2 >= 78 {
        sseqid=$1; slen=$4; aln=$3; bs=$5; pident=$2
        los=sseqid; sub(/.*_LOS_/,"",los)
        # collapse extracted allele suffix: A_1 -> A, S_3 -> S, etc.
        # mosaic alleles: IJ_1 -> I_J, IS_1 -> I_S (match existing mosaic labels)
        sub(/_[0-9]+$/,"",los)
        if(los=="IJ") los="I_J"; else if(los=="IS") los="I_S"
        # Aggregate per ALLELE (full sseqid), not per class, so multiple
        # alleles of the same class are NOT summed together (which would
        # inflate coverage past 100%). We sum segments within one allele,
        # then pick the single best allele per class in END.
        aln_a[sseqid]  += aln
        slen_a[sseqid]  = slen
        cls_a[sseqid]   = los
        if(bs > bs_a[sseqid]){ bs_a[sseqid]=bs; pid_a[sseqid]=pident }
      }
      END {
        # for each class keep the allele with the most aligned bases
        for(a in aln_a){
          c = cls_a[a]
          if(aln_a[a] > total_aln[c]){
            total_aln[c] = aln_a[a]
            slen_map[c]  = slen_a[a]
            best_bs[c]   = bs_a[a]
            best_pid[c]  = pid_a[a]
          }
        }
        for(l in total_aln)
          printf "%s\t%.3f\t%d\t%d\t%.0f\n",
            l, best_pid[l], total_aln[l], slen_map[l], best_bs[l]
      }' > "$out_los"
  ) &
done
wait

# FIX: diagnose malformed assemblies. If most samples produced no hits, the
# files in $INPUT are almost certainly not valid nucleotide FASTA.
N_EMPTY=$(find "$OUTDIR/los_blast_raw" -name "*_LOS.blast" -size 0 2>/dev/null | wc -l)
if [ "$N_EMPTY" -gt 0 ]; then
  echo "  WARNING: $N_EMPTY/$N_SAMPLES samples produced no LOS BLAST hits."
  if [ "$N_EMPTY" -eq "$N_SAMPLES" ]; then
    echo "           ALL samples empty -> assemblies in '$INPUT' are not being parsed as FASTA."
    echo "           Inspect one file:   head -3 \"\$(ls $INPUT/*.{fa,fasta,fna} 2>/dev/null | head -1)\""
    echo "           (look for: gzip/binary content, coordinate-numbered sequence lines,"
    echo "            UTF-16/BOM encoding, or a non-FASTA file with a .fa/.fasta/.fna name)"
  fi
fi

########################################
# LOS B-WINDOW DIAGNOSTIC COUNT
# Used ONLY for close A/B rescue.
#
# Diagnostic B windows from validation:
#   Cj81-176_LOS_B_sliding:5001-6000
#   Cj81-176_LOS_B_sliding:6001-7000
#   Cj81-176_LOS_B_sliding:8001-9000
#   Cj81-176_LOS_B_sliding:9001-10000
#   Cj81-176_LOS_B_sliding:16001-17000
#
# A strong diagnostic window hit requires:
#   identity >= 95%
#   aligned length >= 800 bp
#
# Output:
#   $OUTDIR/los_b_windows.tsv
#   columns: ID  B_DIAGNOSTIC_WINDOW_COUNT
########################################

echo "  Counting LOS B diagnostic windows for A/B close-case rescue..."

mkdir -p "$OUTDIR/los_b_window_hits"
> "$OUTDIR/los_b_windows.tsv"

# FIX: missing B-window DB is no longer fatal. Without it we simply skip the
# A/B close-case rescue (counts default to 0); the LOS step still finishes
# and los_final.tsv is still written.
if [ ! -f "${LOS_B_WINDOW_DB}.nin" ] && [ ! -f "${LOS_B_WINDOW_DB}.ndb" ]; then
  echo "  NOTE: LOS B-window DB not found ($LOS_B_WINDOW_DB) — skipping A/B rescue."
  while IFS= read -r s; do printf "%s\t0\n" "$s"; done < "$SAMPLES" \
    > "$OUTDIR/los_b_windows.tsv"
else

for f in "$INPUT"/*.{fa,fasta,fna}; do
  [ -f "$f" ] || continue
  n=$(basename "$f"); n=${n%.*}
  out_bwin="$OUTDIR/los_b_window_hits/${n}_BWIN.blast"

  # FIX: route through clean_query; '|| true' so one bad sample can never
  # abort the whole pipeline under 'set -e' (this was why los_final.tsv was
  # never written). Empty $out_bwin just means 0 diagnostic windows.
  clean_query "$f" \
  | run_blastn -db "$LOS_B_WINDOW_DB" \
    -outfmt "6 sseqid pident length" \
    -evalue 1e-10 -num_threads "$THREADS_PER_JOB" 2>/dev/null \
    > "$out_bwin" || true

  # FIX (mawk): single-line condition. mawk does not allow a line break after
  # '&&' before '(' (gawk does), which was an unreached parse bug. The regex
  # matches exactly the same 5 diagnostic windows at id>=95 / aln>=800.
  count=$(awk '$2>=95 && $3>=800 && $1 ~ /^Cj81-176_LOS_B_sliding:(5001-6000|6001-7000|8001-9000|9001-10000|16001-17000)$/ {print $1}' "$out_bwin" | sort -u | wc -l)

  printf "%s\t%d\n" "$n" "$count" >> "$OUTDIR/los_b_windows.tsv"
done

fi  # end B-window DB present

########################################
# ORF11 BLAST — A vs R disambiguation
# Runs for ALL samples in parallel alongside the other LOS blasts.
# Results used in the LOS calling awk below to override A/R calls
# that landed at MID or LOOSE tier only (STRICT tier BLAST is trusted as-is).
#
# DB: $ORF11_DB  (default: databases/vfdb/orf11/orf11_db)
# Threshold: 90% identity, 80% coverage
# Output: $OUTDIR/los_blast_raw/orf11_hits.tsv  (ID <TAB> YES|NO)
########################################

echo "  Running orf11 BLAST for A/R disambiguation..."

ORF11_HITS="$OUTDIR/los_blast_raw/orf11_hits.tsv"

if [ ! -f "${ORF11_DB}.nin" ] && [ ! -f "${ORF11_DB}.ndb" ]; then
  echo "  NOTE: orf11 DB not found ($ORF11_DB) — A/R orf11 correction skipped."
  while IFS= read -r s; do printf "%s\tUNKNOWN\n" "$s"; done < "$SAMPLES" > "$ORF11_HITS"
else
  > "$ORF11_HITS"
  for f in "$INPUT"/*.{fa,fasta,fna}; do
    [ -f "$f" ] || continue
    n=$(basename "$f"); n=${n%.*}
    limit_jobs "$JOBS"
    (
      hit=$(clean_query "$f" \
        | run_blastn -db "$ORF11_DB" \
          -outfmt "6 pident length slen" \
          -evalue 1e-10 -num_threads "$THREADS_PER_JOB" 2>/dev/null \
        | awk '{ cov=($3>0?$2/$3:0); if($1>=90 && cov>=0.80){ print "YES"; exit } }
               END{ if(!found) print "NO" }' found=0)
      [ -z "$hit" ] && hit="NO"
      # atomic write
      (
        flock 7
        printf "%s\t%s\n" "$n" "$hit" >> "$ORF11_HITS"
      ) 7>"$OUTDIR/los_blast_raw/.orf11_lock"
    ) &
  done
  wait
  echo "  orf11 hits: $(awk '$2=="YES"' "$ORF11_HITS" | wc -l) / $N_SAMPLES samples positive"
fi
# Scoring logic:
#   BEST HIT:   pid * coverage           (no bitscore — removes reference
#                                          length bias; fixes A/B and H/P)
#   SECOND HIT: pid * coverage * bitscore (keeps bitscore to correctly rank
#                                          J above S for I/J mosaic detection)
#
# Mosaic rules:
#   1. WHITELIST — only I/J and I/S can ever be called as mosaics
#      (confirmed true mosaics per Parker 2008 and BLAST validation)
#      Everything else: best hit wins, no mosaic called
#   2. Both classes must independently pass 80% coverage threshold
#   3. Identity difference must be <= 1.0%
#   4. Only at STRICT or MID tier — never LOOSE
########################################

awk -v SAMPLES="$SAMPLES" -v BWIN="$OUTDIR/los_b_windows.tsv" -v ORF11="$ORF11_HITS" '
BEGIN{
  FS=OFS="\t"
  while((getline line < SAMPLES)>0) samples[line]=1
  while((getline line < BWIN)>0){
    split(line,bw,OFS)
    bwin[bw[1]]=bw[2]+0
  }
  # orf11 lookup: sample -> YES / NO / UNKNOWN
  while((getline line < ORF11)>0){
    split(line,orf,OFS)
    orf11[orf[1]]=orf[2]
  }
  # WHITELIST — only these pairs can be called as true mosaics
  true_mosaic["I_J"]=1; true_mosaic["J_I"]=1
  true_mosaic["I_S"]=1; true_mosaic["S_I"]=1
}
function tier_from_cov_pid(pid, cov){
  if(pid >= 99.5 && cov >= 0.96) return 1
  if(pid >= 99.0 && cov >= 0.95) return 2
  if(pid >= 97.0 && cov >= 0.93) return 3
  if(pid >= 95.0 && cov >= 0.93) return 4
  if(pid >= 94.0 && cov >= 0.80) return 5
  if(pid >= 92.0 && cov >= 0.80) return 6
  if(pid >= 90.0 && cov >= 0.60) return 7
  return 0
}
function is_sialylated(x){
  if(x=="A"||x=="B"||x=="C"||x=="M"||x=="R") return 1
  return 0
}
{
  split(FILENAME,a,"/"); file=a[length(a)]
  sub(/_LOS\.blast$/,"",file); sample=file
  los=$1; pid=$2+0; aln=$3+0; cov=($4>0?$3/$4:0); bs=$5+0

  tier=tier_from_cov_pid(pid,cov)
  if(tier==0) next

  # BEST HIT: ranked by pid*cov (no bitscore — removes reference length bias)
  # Special case: A vs B uses pid*aln as tiebreaker — A ref (11474bp) is much
  # smaller than B ref (17001bp) so pid*cov always favors A even for true B strains.
  # pid*aln correctly picks whichever class has more total bases aligned.
  sc_best = pid * cov

  # SECOND HIT: ranked by pid*cov*bitscore (keeps bitscore to rank J above S)
  sc_2nd = pid * cov * bs

  # Track best hit
  if(sc_best > best_score[sample,tier]){
    best_score[sample,tier] = sc_best
    best_los[sample,tier]   = los
    best_cov[sample,tier]   = cov
    best_pid[sample,tier]   = pid
    best_aln[sample,tier]   = aln
  } else if(sc_best == best_score[sample,tier]){
    # tiebreaker: pid*aln
    if(pid*aln > best_pid[sample,tier]*best_aln[sample,tier]){
      best_score[sample,tier] = sc_best
      best_los[sample,tier]   = los
      best_cov[sample,tier]   = cov
      best_pid[sample,tier]   = pid
      best_aln[sample,tier]   = aln
    }
  }

  # Track second hit independently (pid*cov*bs), must be different class from best
  # We do this in END block after best hits are finalized
  # Store all hits per sample/tier for post-processing
  hits[sample,tier,los]     = sc_2nd
  hits_cov[sample,tier,los] = cov
  hits_pid[sample,tier,los] = pid
  hits_aln[sample,tier,los] = aln
}
END{
  print "ID","LOS_CLASS","GBS_RISK","SIALYLATED","LOS_TIER","LOS_COVERAGE","LOS_IDENTITY"
  for(s in samples){
    los="-"; risk="-"; sial="-"; tier_name="-"; cov_out="-"; pid_out="-"

    # Store all hits for special pair rules (A/B, H/P, S/J)
    # Build global best pid*cov score and pid*cov*(bs/slen) score per class
    delete _pcov; delete _pnorm; delete _pt; delete _pcovv; delete _ppid; delete _paln; delete _pbs
    for(key in hits_aln){
      n=split(key,kp,SUBSEP)
      if(kp[1]!=s) continue
      cc=kp[3]; tt=kp[2]+0
      if(cc!="A"&&cc!="B"&&cc!="H"&&cc!="P"&&cc!="S"&&cc!="J") continue
      cv=hits_cov[key]; pp=hits_pid[key]; al=hits_aln[key]; bv=hits[key]/(cv>0?cv*pp:1)
      # bv is bitscore — recover from sc_2nd = pid*cov*bs
      bs_val = (cv>0 && pp>0) ? hits[key]/(pp*cv) : 0
      slen_val = (cv>0 && al>0) ? al/cv : 1
      sc_pc   = pp * cv
      sc_norm = pp * cv * (bs_val / slen_val)
      if(!(cc in _pcov) || sc_pc > _pcov[cc]){
        _pcov[cc]=sc_pc; _pnorm[cc]=sc_norm; _pt[cc]=tt
        _pcovv[cc]=cv; _ppid[cc]=pp; _paln[cc]=al
      }
    }

    for(t=1;t<=7;t++){
      if(!((s,t) in best_los)) continue

      b_los = best_los[s,t]
      b_cov = best_cov[s,t]
      b_pid = best_pid[s,t]
      b_aln = best_aln[s,t]+0

      # ── A/B close-call structural classifier ─────────────────────
      # Original BLAST caller remains primary.
      #
      # ONLY activate this logic when:
      #   A and B are close competitors (diff <=5%)
      #
      # Then use statistically learned structural separation:
      #   TRUE B usually has much larger B alignment span.
      #
      # Learned from validation set:
      #   ALN_DELTA = B_ALN - A_ALN
      #
      # Rule:
      #   ALN_DELTA > 5000  -> B
      #   otherwise         -> A

      if(("A" in _pcov) && ("B" in _pcov)){

        a_pc=_pcov["A"]
        b_pc=_pcov["B"]

        diff_pct=(a_pc>b_pc ?
                  (a_pc-b_pc)/a_pc :
                  (b_pc-a_pc)/b_pc) * 100

        # activate ONLY for close A/B competition
        if(diff_pct <= 5.0){

          aln_delta = _paln["B"] - _paln["A"]

          # statistically learned B structural extension
          if(aln_delta > 5000){

            b_los="B"
            b_cov=_pcovv["B"]
            b_pid=_ppid["B"]
            b_aln=_paln["B"]
            t=_pt["B"]

          } else {

            b_los="A"
            b_cov=_pcovv["A"]
            b_pid=_ppid["A"]
            b_aln=_paln["A"]
            t=_pt["A"]
          }
        }
      }

      # ── H/P close-case rule (same logic) ─────────────────────────
      if((b_los=="H"||b_los=="P") && ("H" in _pcov) && ("P" in _pcov)){
        other = (b_los=="H" ? "P" : "H")
        h_pc=_pcov["H"]; p_pc=_pcov["P"]
        diff_pct=(h_pc>p_pc?(h_pc-p_pc)/h_pc:(p_pc-h_pc)/p_pc)*100
        same_tier=(_pt["H"]==_pt["P"])
        if(same_tier || diff_pct <= 5.0){
          if(_pnorm[other] > _pnorm[b_los]){
            b_los=other; b_cov=_pcovv[other]; b_pid=_ppid[other]
            b_aln=_paln[other]; t=_pt[other]
          }
        } else if(_pt[other] < t){
          b_los=other; b_cov=_pcovv[other]; b_pid=_ppid[other]
          b_aln=_paln[other]; t=_pt[other]
        }
      }

      # ── S/J close-case rule (same logic) ─────────────────────────
      if((b_los=="S"||b_los=="J") && ("S" in _pcov) && ("J" in _pcov)){
        other = (b_los=="S" ? "J" : "S")
        s_pc=_pcov["S"]; j_pc=_pcov["J"]
        diff_pct=(s_pc>j_pc?(s_pc-j_pc)/s_pc:(j_pc-s_pc)/j_pc)*100
        same_tier=(_pt["S"]==_pt["J"])
        if(same_tier || diff_pct <= 5.0){
          if(_pnorm[other] > _pnorm[b_los]){
            b_los=other; b_cov=_pcovv[other]; b_pid=_ppid[other]
            b_aln=_paln[other]; t=_pt[other]
          }
        } else if(_pt[other] < t){
          b_los=other; b_cov=_pcovv[other]; b_pid=_ppid[other]
          b_aln=_paln[other]; t=_pt[other]
        }
      }

      # Find second best hit by pid*cov*bs, excluding the best class
      s_los=""; s_cov=0; s_pid=0; s_sc=0
      for(key in hits){
        n=split(key,parts,SUBSEP)
        if(parts[1]==s && parts[2]==t && parts[3]!=b_los){
          sc=hits[key]
          if(sc > s_sc){
            s_sc  = sc
            s_los = parts[3]
            s_cov = hits_cov[key]
            s_pid = hits_pid[key]
          }
        }
      }

      # Build candidate pair name (alphabetical order)
      pair = (b_los < s_los ? b_los"_"s_los : s_los"_"b_los)

      # Tier name mapping
      tname[1]="STRICT"; tname[2]="STRICT_2"; tname[3]="STRICT_3"; tname[4]="STRICT_4"
      tname[5]="MID"; tname[6]="MID_2"; tname[7]="LOOSE"

      los       = b_los
      tier_name = tname[t]

      cov_out = sprintf("%.1f", b_cov*100)
      pid_out = sprintf("%.1f", b_pid)

      # ── QC FLOOR ──────────────────────────────────────────────────
      # cov < 60%            -> untypable (-)
      # 60% <= cov < 68%     -> keep only if identity >= 97%, else (-)
      # cov >= 68%           -> keep
      if(los != "-"){
        if(b_cov < 0.60){
          los="-"; tier_name="-"; risk="-"; sial="-"
        } else if(b_cov < 0.68 && b_pid < 97.0){
          los="-"; tier_name="-"; risk="-"; sial="-"
        }
      }

      if(los=="-"){
        # untypable: leave risk/sial as "-"
      } else if(is_sialylated(los)){ risk="HIGH"; sial="YES" }
      else                  { risk="LOW";  sial="NO"  }
      break
    }

    # ── orf11 A/R correction ──────────────────────────────────────────────
    # Only applies at non-STRICT tiers (MID=5, MID_2=6, LOOSE=7).
    # At STRICT tiers (1-4) the BLAST score is trusted as-is.
    # orf11 present -> A  |  orf11 absent -> R
    # UNKNOWN (DB missing) -> no change
    if((los=="A"||los=="R") && t>=5){
      o=orf11[s]
      if(o=="YES")    { los="A"; sial="YES"; risk="HIGH" }
      else if(o=="NO"){ los="R"; sial="YES"; risk="HIGH" }
      # o=="UNKNOWN" -> leave call as-is
    }

    print s,los,risk,sial,tier_name,cov_out"%",pid_out"%"
  }
}
' "$OUTDIR"/los_blast_raw/*_LOS.blast \
  > "$OUTDIR/los_final.tsv"

fi

########################################
# VF BLAST
########################################
if [ "$SKIP_VF" -eq 1 ]; then
  echo "[5/6] Virulence factor detection... SKIPPED"

  VF_GENES=("cstII" "T4SS")
  printf "ID\tcstII\tT4SS\n" > "$OUTDIR/vf_type.tsv"
  while IFS= read -r sample; do
    printf "%s\t-\t-\n" "$sample"
  done < "$SAMPLES" >> "$OUTDIR/vf_type.tsv"

else

echo "[5/6] Virulence factor detection... RUNNING"

VF_PID=90
VF_PCOV=0.60
T4SS_TOTAL=11
VF_COMBINED="$VF_DB/vf/vf"
CSTII_DB="$VF_DB/cstII/cstII"
T4SS_DB="$VF_DB/T4SS/T4SS.fasta"

# Derive gene list from combined FASTA headers (for column order in output)
VF_GENES_COMBINED=($(grep "^>" "$VF_DB/vf/vf.fasta" \
  | sed 's/^>//' | sed 's/_[0-9]*$//' | sort -u))

for f in "$INPUT"/*.{fa,fasta,fna}; do
  [ -f "$f" ] || continue
  n=$(basename "$f"); n=${n%.*}

  if [ "$RESUME" -eq 1 ] && is_done "$n"; then
    count=$(cat "$PROGRESS_FILE" 2>/dev/null || echo "?")
    echo "  [${count}/${N_SAMPLES}] skip: $n"
    continue
  fi

  limit_jobs "$JOBS"
  (
    vf_out="$OUTDIR/vf_blast/${n}.vf"
    > "$vf_out"

    # ── 1. Combined VF DB — one BLAST for all standard genes ────────────
    run_blastn -query "$f" -db "$VF_COMBINED" \
      -outfmt "6 sseqid pident length slen bitscore" \
      -evalue 1e-10 -num_threads "$THREADS_PER_JOB" \
    | awk -v pid="$VF_PID" -v pcov="$VF_PCOV" '
        {
          cov=($4>0?$3/$4:0)
          if($2>=pid && cov>=pcov && $5>best_score[$1]){
            best_score[$1]=$5
            best_line[$1]=$0
            # Extract gene name by stripping trailing _N allele number
            gene=$1; sub(/_[0-9]+$/,"",gene)
            if($5>gene_score[gene]){
              gene_score[gene]=$5
              best_seqid[gene]=$1
              best_pid[gene]=$2
              best_cov[gene]=cov
            }
          }
        }
        END{
          for(g in gene_score){
            pid_r=int(best_pid[g]+0.5)
            cov_r=int(best_cov[g]*100+0.5)
            printf "%s\t%s_%s (%d%%ID %d%%COV)\n",
              g, g, best_seqid[g], pid_r, cov_r
          }
        }' >> "$vf_out"

    # ── 2. cstII standalone — Asn51/Thr51 allele typing ─────────────────
    hit=$(run_blastn -query "$f" -db "$CSTII_DB" \
      -outfmt "6 qseqid qstart qend sstrand pident length slen bitscore sstart" \
      -evalue 1e-10 -num_threads "$THREADS_PER_JOB" \
      -max_hsps 1 -max_target_seqs 1 2>/dev/null \
      | awk -v pid="$VF_PID" -v pcov="$VF_PCOV" \
          '{ cov=($7>0?$6/$7:0); if($5>=pid && cov>=pcov){ print; exit } }')

    if [ -z "$hit" ]; then
      printf "cstII\t-\n" >> "$vf_out"
    else
      qseqid=$(echo "$hit" | awk '{print $1}')
      qstart=$( echo "$hit" | awk '{print $2}')
      qend=$(   echo "$hit" | awk '{print $3}')
      strand=$( echo "$hit" | awk '{print $4}')
      best_pid_v=$(echo "$hit" | awk '{print $5}')
      aln_len=$( echo "$hit" | awk '{print $6}')
      slen=$(    echo "$hit" | awk '{print $7}')
      sstart=$(  echo "$hit" | awk '{print $9}')
      pid_r=$(echo "$best_pid_v" | awk '{printf "%d", int($1+0.5)}')
      cov_r=$(echo "$aln_len $slen" | awk '{printf "%d", int($1/$2*100+0.5)}')

      codon51=$(awk \
        -v target="$qseqid" -v qstart="$qstart" -v qend="$qend" \
        -v strand="$strand"  -v sstart="$sstart" \
        'BEGIN{ codon_offset=151-sstart; collecting=0; seq="" }
         /^>/{ name=substr($0,2); split(name,a," "); name=a[1]
               collecting=(name==target); next }
         collecting{ seq=seq $0 }
         END{
           gsub(/[^ACGTacgtNn]/,"",seq)
           if(strand=="plus"){
             pos=qstart+codon_offset; codon=toupper(substr(seq,pos,3))
           } else {
             pos=qend-codon_offset-2; raw=toupper(substr(seq,pos,3))
             n=split(raw,c,""); rc=""
             for(i=n;i>=1;i--){
               if(c[i]=="A") rc=rc"T"; else if(c[i]=="T") rc=rc"A"
               else if(c[i]=="G") rc=rc"C"; else if(c[i]=="C") rc=rc"G"
               else rc=rc"N" }
             codon=rc }
           print codon }' "$f")

      allele="unk51"
      case "$codon51" in
        AAT|AAC)         allele="asn51" ;;
        ACT|ACC|ACA|ACG) allele="thr51" ;;
      esac

      printf "cstII\tcstII_%s (%d%%ID %d%%COV)\n" \
        "$allele" "$pid_r" "$cov_r" >> "$vf_out"
    fi

    # ── 3. T4SS standalone — count distinct subunits ────────────────────
    result=$(run_blastn -query "$f" -db "$T4SS_DB" \
      -outfmt "6 sseqid pident length slen bitscore" \
      -evalue 1e-10 -num_threads "$THREADS_PER_JOB" 2>/dev/null \
      | awk -v pid="$VF_PID" -v pcov="$VF_PCOV" -v total="$T4SS_TOTAL" '
          {
            cov=($4>0?$3/$4:0)
            if($2>=pid && cov>=pcov){
              gene_name=$1; sub(/_[0-9]+$/,"",gene_name)
              seen[gene_name]=1
            }
          }
          END{
            n=0; for(k in seen) n++
            if(n>0) printf "T4SS_%d/%d",n,total
            else    printf "-"
          }')

    printf "T4SS\t%s\n" "$result" >> "$vf_out"

    mark_done "$n"
  ) &
done
wait

########################################
# VF CALLING
########################################

# Gene list = combined genes + cstII + T4SS (sorted)
VF_GENES=()
for g in "${VF_GENES_COMBINED[@]}"; do VF_GENES+=("$g"); done
VF_GENES+=("cstII" "T4SS")
IFS=$'\n' VF_GENES=($(printf '%s\n' "${VF_GENES[@]}" | sort -u))
unset IFS

{
  printf "ID"
  for g in "${VF_GENES[@]}"; do printf "\t%s" "$g"; done
  printf "\n"
  while IFS= read -r sample; do
    printf "%s" "$sample"
    vf_file="$OUTDIR/vf_blast/${sample}.vf"
    for g in "${VF_GENES[@]}"; do
      if [ -f "$vf_file" ]; then
        val=$(awk -F'\t' -v gene="$g" '$1==gene{print $2; exit}' "$vf_file")
        [ -z "$val" ] && val="-"
      else val="-"; fi
      printf "\t%s" "$val"
    done
    printf "\n"
  done < "$SAMPLES"
} > "$OUTDIR/vf_type.tsv"

fi

########################################
# neuA1 → LOS C CORRECTION
# If --skip_vf was used there is no VF data, so this block is skipped
# and the raw BLAST LOS call is kept as-is.
# If VF was run and neuA1 is detected for a sample, override the LOS
# class to C in los_final.tsv (C is sialylated → SIALYLATED=YES, GBS_RISK=HIGH).
########################################

if [ "$SKIP_VF" -ne 1 ]; then
  echo "  Applying neuA1 → LOS C correction..."
  python3 - "$OUTDIR/los_final.tsv" "$OUTDIR/vf_type.tsv" <<'PYEOF'
import sys, csv

los_file, vf_file = sys.argv[1], sys.argv[2]

# Find neuA1 column index in vf_type.tsv
neua1_samples = set()
try:
    with open(vf_file, newline="") as fv:
        reader = csv.DictReader(fv, delimiter="\t")
        for row in reader:
            val = row.get("neuA1", "-").strip()
            if val not in ("-", "", "absent", "ABSENT"):
                neua1_samples.add(row["ID"])
except Exception as e:
    print(f"  WARNING: could not read vf_type.tsv: {e}")

if not neua1_samples:
    print("  neuA1 → C: no samples with neuA1 detected, nothing to override")
    sys.exit(0)

# Rewrite los_final.tsv in place
rows = []
changed = 0
with open(los_file, newline="") as fl:
    reader = csv.DictReader(fl, delimiter="\t")
    fieldnames = reader.fieldnames
    for row in reader:
        if row["ID"] in neua1_samples and row["LOS_CLASS"] != "C":
            print(f"  {row['ID']}: LOS {row['LOS_CLASS']} → C (neuA1 present)")
            row["LOS_CLASS"]  = "C"
            row["SIALYLATED"] = "YES"
            row["GBS_RISK"]   = "HIGH"
            changed += 1
        rows.append(row)

with open(los_file, "w", newline="") as fl:
    writer = csv.DictWriter(fl, fieldnames=fieldnames, delimiter="\t",
                            lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

print(f"  neuA1 → C: {changed} sample(s) overridden")
PYEOF
fi

########################################
# AMR — AMRFinder Plus
# -----------------------------------------------------------------------
# Runs amrfinder --organism Campylobacter --plus on each assembly.
# Detects:
#   Acquired genes : tet(O), blaOXA variants, ant(6)-Ia, aph(3')-IIIa,
#                    cmeABC efflux, sat4, etc.
#   Point mutations: gyrA T86I/D90N/P104S  (fluoroquinolone resistance)
#                    23S A2075G/A2074C      (macrolide resistance)
#                    rplV A103V             (macrolide resistance)
#                    blaOXA-61 promoter G-57T
#
# Output per sample: results/amr_raw/<sample>_amr.tsv
# Combined        : results/amr_combined.tsv
# For final table : one column per unique AMR gene across all samples
#                   value = "<gene> (XX%ID XX%COV)" or "-"
#                   Only "core" AMR genes included (Type=AMR), STRESS skipped
########################################
if [ "$SKIP_AMR" -eq 1 ]; then
  echo "[6/6] AMR detection (AMRFinder Plus)... SKIPPED"

  AMR_COMBINED="$OUTDIR/amr_combined.tsv"
  printf "Name\tElement symbol\tElement name\tClass\tType\t%% Coverage of reference\t%% Identity to reference\n" > "$AMR_COMBINED"

else

echo "[6/6] AMR detection (AMRFinder Plus)... RUNNING"

# ---------------------------------------------------------------------------
# AMR robustness guard.
# The whole script runs under `set -euo pipefail`. AMRFinder can legitimately
# exit non-zero (DB quirks, one bad contig, transient failures) which, under
# set -e, would silently kill the ENTIRE pipeline before the merge step ever
# runs — leaving no final table. We relax set -e for this stage only, then
# restore it before MERGE. AMR failures now degrade gracefully (empty AMR
# columns) instead of aborting the run.
# ---------------------------------------------------------------------------
set +e

AMR_RAW="$OUTDIR/amr_raw"
mkdir -p "$AMR_RAW"
AMR_COMBINED="$OUTDIR/amr_combined.tsv"
AMR_SKIP=0

if [ "$AMR_AVAILABLE" -eq 0 ]; then
  echo "  WARNING: amrfinder not available — AMR columns will be empty"
  echo "  To install: $0 -I"
  echo "  To update DB: $0 -amrfinder_db"
  AMR_SKIP=1
fi

if [ "$AMR_SKIP" -eq 0 ]; then

  THREADS_PER_JOB_AMR=$(( THREADS / JOBS ))
  [ "$THREADS_PER_JOB_AMR" -lt 1 ] && THREADS_PER_JOB_AMR=1

  # Per-sample progress for AMR (this is the slowest stage — show it working).
  # Uses its own counter/lock so it does NOT touch the global mark_done counter.
  AMR_PROG="$OUTDIR/.done/.amr_progress"
  echo 0 > "$AMR_PROG"
  echo "  (AMRFinder is the slowest stage; ${N_SAMPLES} genomes at ${JOBS} jobs — please let it run)"

  for f in "$INPUT"/*.{fa,fasta,fna}; do
    [ -f "$f" ] || continue
    n=$(basename "$f"); n=${n%.*}
    out_amr="$AMR_RAW/${n}_amr.tsv"

    # resume: skip if already done
    [ "$RESUME" -eq 1 ] && [ -s "$out_amr" ] && continue

    limit_jobs "$JOBS"
    (
      run_amrfinder \
        -n "$f" \
        --name "$n" \
        --organism Campylobacter \
        --plus \
        --threads "$THREADS_PER_JOB_AMR" \
        -o "$out_amr" \
        --quiet 2>/dev/null || true
      # atomic progress print (does not affect the global done counter)
      (
        flock 8
        c=$(cat "$AMR_PROG" 2>/dev/null || echo 0)
        c=$(( c + 1 ))
        echo "$c" > "$AMR_PROG"
        echo "  [AMR ${c}/${N_SAMPLES}] $n"
      ) 8>"$OUTDIR/.done/.amr_lock"
    ) &
  done
  wait

  # Combine all per-sample AMR files into one TSV
  # Header from first file, data from all
  first=$(ls "$AMR_RAW"/*_amr.tsv 2>/dev/null | head -n1)
  if [ -n "$first" ]; then
    head -n1 "$first" > "$AMR_COMBINED"
    for f in "$AMR_RAW"/*_amr.tsv; do
      tail -n +2 "$f" >> "$AMR_COMBINED"
    done
    N_AMR=$(tail -n +2 "$AMR_COMBINED" | wc -l)
    echo "  AMR: $N_AMR total hits across all samples → $AMR_COMBINED"
  else
    echo "  AMR: no output files found"
    AMR_SKIP=1
  fi

fi

fi

# Restore strict error handling now that the AMR stage is complete.
set -e

########################################
# MERGE
########################################

echo "[merge] Merging results..."

{
  printf "ID\tLOS_CLASS\tGBS_RISK\tSIALYLATED\tHS_TYPE\tMLST_ST\tMLST_CC"
  for g in "${VF_GENES[@]}"; do printf "\t%s" "$g"; done
  printf "\n"

  join -t $'\t' -a1 -a2 -e "-" -o '0,1.2,1.3,1.4,2.2' \
    <(tail -n +2 "$OUTDIR/los_final.tsv" | sort -k1,1) \
    <(tail -n +2 "$OUTDIR/hs_type.tsv"   | sort -k1,1) \
  | sort -k1,1 \
  | join -t $'\t' -a1 -a2 -e "-" -o '0,1.2,1.3,1.4,1.5,2.2,2.3' \
      - \
      <(tail -n +2 "$OUTDIR/mlst.tsv" | sort -k1,1) \
  | sort -k1,1 \
  | join -t $'\t' -a1 -a2 -e "-" \
      -o "0,1.2,1.3,1.4,1.5,1.6,1.7$(printf ',2.%d' $(seq 2 $((${#VF_GENES[@]}+1))))" \
      - \
      <(tail -n +2 "$OUTDIR/vf_type.tsv" | sort -k1,1)
} > "$OUTDIR/merged.tsv"

########################################
# INTERPRETATION ENGINE
########################################

echo "[6/6] Interpreting results..."

python3 - \
  "$OUTDIR/merged.tsv" \
  "$OUTDIR/qc_report.tsv" \
  "${AMR_COMBINED:-/dev/null}" \
  "$OUTDIR/final_los_hs_gbs.tsv" \
<< 'PYEOF'
import sys, csv, re

in_file, qc_file, amr_file, out_file = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

GBS_HS = {
    "HS1","HS1C","HS1/44","HS1/44C","HS44","HS44C",
    "HS2","HS2C",
    "HS4","HS4C",
    "HS19","HS19C",
    "HS23","HS23C","HS36","HS36C","HS23/36","HS23/36C",
    "HS41","HS41C",
}

def normalize_hs(hs):
    # Canonicalize an HS call for GBS-marker lookup ONLY.
    # Does NOT modify the HS_TYPE output column.
    #   "HS23_36c(HS23_36*)" -> "HS23/36C"   "HS4c(HS4*)" -> "HS4C"   "HS1c" -> "HS1C"
    hs = (hs or "-").strip().upper()
    hs = hs.split("(")[0]      # drop bracketed subtype
    hs = hs.replace("_", "/")  # underscore complex form -> slash
    return hs
CDT_GENES  = ["cdtA","cdtB","cdtC"]
SIAL_GENES = ["cstII","cstIII","neuA1","neuB","neuB1","neuC","neuC1"]
CORE_GENES = ["gmhA","gmhA2","gmhB","hddA","hddC","hldD","hldE",
              "rfaD","rfaE1","waaC","waaF"]
KPS_GENES  = ["kpsC","kpsD","kpsE","kpsF","kpsM","kpsS","kpsT"]
T4SS_TOTAL = 11

def present(v):
    return v not in ("-","",None)

def interpret(row):
    def val(g): return row.get(g,"-")

    cstii_val = val("cstII")
    if   not present(cstii_val):   cstii_func = "ABSENT"
    elif "asn51" in cstii_val:     cstii_func = "BIFUNCTIONAL"
    elif "thr51" in cstii_val:     cstii_func = "MONOFUNCTIONAL"
    else:                          cstii_func = "UNKNOWN"

    found_sial      = [g for g in SIAL_GENES if present(val(g))]
    cgtA            = present(val("cgtA"))
    cgtB            = present(val("ctgB")) or present(val("cgtB"))
    wlan            = present(val("wlaN"))
    has_transferase = present(val("cstII")) or present(val("cstIII"))
    neub            = present(val("neuB")) or present(val("neuB1"))

    if not found_sial:
        mimic, mimic_risk = "NON-SIALYLATED", "LOW"
    elif not has_transferase:
        mimic, mimic_risk = "NO_TRANSFERASE", "LOW_UNCERTAIN"
    elif cstii_func == "BIFUNCTIONAL":
        if cgtA and cgtB:   mimic, mimic_risk = "GT1a", "HIGH"
        else:               mimic, mimic_risk = "GD3",  "MODERATE"
    elif cstii_func == "MONOFUNCTIONAL":
        if cgtA and cgtB:   mimic, mimic_risk = "GM1/GD1a", "HIGH"
        else:               mimic, mimic_risk = "GM3", "MODERATE"
    elif cstii_func == "UNKNOWN":          # cstII present, codon 51 unreadable
        if cgtA and cgtB:   mimic, mimic_risk = "GM1/GD1a|GT1a (unresolved)", "HIGH_UNCERTAIN"
        else:               mimic, mimic_risk = "sialylated (unresolved)", "MODERATE"
    elif present(val("cstIII")):
        if wlan:            mimic, mimic_risk = "GM2/GM1", "MODERATE"
        else:               mimic, mimic_risk = "GM3",     "MODERATE"
    elif found_sial:        mimic, mimic_risk = "UNKNOWN", "MODERATE"
    else:                   mimic, mimic_risk = "NON-SIALYLATED", "LOW"

    # HS_GBS_MARKER: match against canonical names (handles complex/bracketed calls)
    hs_gbs = "YES" if normalize_hs(row.get("HS_TYPE","-")) in GBS_HS else "NO"

    # COMBINED_GBS: mimic risk is primary; HS marker corroborates; SIALYLATED gates HIGH
    los_high = row.get("SIALYLATED","-").strip().upper() == "YES"
    if   mimic_risk == "HIGH"           and los_high:                   combined = "HIGH"
    elif mimic_risk == "HIGH":                                          combined = "MODERATE"
    elif mimic_risk == "HIGH_UNCERTAIN" and hs_gbs=="YES" and los_high:  combined = "HIGH"
    elif mimic_risk == "HIGH_UNCERTAIN":                                combined = "MODERATE"
    elif mimic_risk == "MODERATE"       and hs_gbs=="YES" and los_high: combined = "HIGH"
    elif mimic_risk == "MODERATE":                                      combined = "MODERATE"
    elif hs_gbs=="YES" and los_high:                                    combined = "MODERATE"
    else:                                                               combined = "LOW"

    found_kps  = [g for g in KPS_GENES  if present(val(g))]
    found_core = [g for g in CORE_GENES if present(val(g))]
    kps_ok     = len(found_kps)  == len(KPS_GENES)
    core_ok    = len(found_core) == len(CORE_GENES)
    if   kps_ok and core_ok:                       serum = "YES"
    elif len(found_kps)==0 and len(found_core)==0: serum = "NO"
    else:                                          serum = "PARTIAL"

    meop_status = "PRESENT" if present(val("MEOP")) else "ABSENT"

    t4ss_val = val("T4SS")
    if not present(t4ss_val):
        t4ss_status = "ABSENT"
    else:
        m = re.search(r'(\d+)/(\d+)', t4ss_val)
        if m:
            nf, nt = int(m.group(1)), int(m.group(2))
            t4ss_status = f"COMPLETE({nf}/{nt})" if nf==nt else f"PARTIAL({nf}/{nt})"
        else:
            t4ss_status = "ABSENT"

    sial_genes = "+".join(sorted(found_sial)) if found_sial else "NONE"

    found_cdt = [g for g in CDT_GENES if present(val(g))]
    n = len(found_cdt)
    if   n==3: cdt_status = "COMPLETE"
    elif n==0: cdt_status = "ABSENT"
    else:      cdt_status = f"PARTIAL({n}/3:{'+'.join(sorted(found_cdt))})"

    missing_core = [g for g in CORE_GENES if not present(val(g))]
    nc, tc = len(found_core), len(CORE_GENES)
    if   nc==tc: core_status = f"COMPLETE({tc}/{tc})"
    elif nc==0:  core_status = "ABSENT"
    else:        core_status = f"PARTIAL({nc}/{tc}:missing={'+'.join(missing_core)})"

    return {
        "CSTII_FUNCTION":    cstii_func,
        "GANGLIOSIDE_MIMIC": mimic,
        "GBS_MIMIC_RISK":    mimic_risk,
        "HS_GBS_MARKER":     hs_gbs,
        "COMBINED_GBS":      combined,
        "SERUM_RESISTANCE":  serum,
        "MEOP_STATUS":       meop_status,
        "T4SS_STATUS":       t4ss_status,
        "CDT_STATUS":        cdt_status,
        "SIALYLATION_GENES": sial_genes,
        "LOS_CORE_STATUS":   core_status,
        "CAPSULE_KPS":       f"{len(found_kps)}/{len(KPS_GENES)}",
        "WLAN_STATUS":       "PRESENT" if wlan else "ABSENT",
    }

INTERP_COLS = [
    "CSTII_FUNCTION","GANGLIOSIDE_MIMIC","GBS_MIMIC_RISK",
    "HS_GBS_MARKER","COMBINED_GBS",
    "SERUM_RESISTANCE","MEOP_STATUS","T4SS_STATUS","CDT_STATUS",
    "SIALYLATION_GENES","LOS_CORE_STATUS","CAPSULE_KPS","WLAN_STATUS",
]
QC_COLS    = ["genome_size_bp","n_contigs","gc_pct","qc_status","qc_warnings"]
FRONT_COLS = ["ID","LOS_CLASS","SIALYLATED","HS_TYPE","GBS_RISK"]

# ── Load QC ────────────────────────────────────────────────────────────
qc = {}
try:
    with open(qc_file, newline="") as fq:
        for row in csv.DictReader(fq, delimiter="\t"):
            qc[row["ID"]] = row
except Exception:
    pass

# ── Load AMR data ───────────────────────────────────────────────────────
# From AMRFinder combined TSV:
#   Column "Name"           = sample ID
#   Column "Element symbol" = gene name e.g. gyrA_T86I, tet(O), blaOXA-193
#   Column "Element name"   = full description
#   Column "Class"          = drug class e.g. QUINOLONE, TETRACYCLINE
#   Column "Type"           = AMR | STRESS | VIRULENCE  (we keep only AMR)
#   Column "% Coverage of reference" and "% Identity to reference"
#
# Build per-sample dict: {sample: {gene_symbol: (pct_id, pct_cov, element_name, class)}}
# Also collect all unique AMR gene symbols across all samples.

amr_data   = {}   # {sample: {symbol: (pct_id, pct_cov, name, drug_class)}}
amr_genes  = []   # ordered unique gene symbols

try:
    with open(amr_file, newline="") as fa:
        reader = csv.DictReader(fa, delimiter="\t")
        seen_genes = {}   # preserve insertion order
        for row in reader:
            # skip non-AMR rows (STRESS, VIRULENCE etc)
            if row.get("Type","").strip().upper() not in ("AMR","POINT"):
                # also include POINT type — that's how point mutations are typed
                # actually AMRFinder uses "AMR" type for both acquired and point
                # check Subtype instead
                subtype = row.get("Subtype","").strip().upper()
                typ     = row.get("Type","").strip().upper()
                if typ not in ("AMR",):
                    continue

            sample  = row.get("Name","").strip()
            symbol  = row.get("Element symbol","").strip()
            name    = row.get("Element name","").strip()
            dclass  = row.get("Class","").strip()
            pct_id  = row.get("% Identity to reference","").strip()
            pct_cov = row.get("% Coverage of reference","").strip()

            if not sample or not symbol:
                continue

            if sample not in amr_data:
                amr_data[sample] = {}

            # format value: gene (ID% COV%)
            try:
                val = f"{symbol} ({float(pct_id):.0f}%ID {float(pct_cov):.0f}%COV)"
            except ValueError:
                val = symbol

            amr_data[sample][symbol] = val

            if symbol not in seen_genes:
                seen_genes[symbol] = (name, dclass)

        amr_genes = list(seen_genes.keys())

except Exception:
    amr_genes = []

# ── Build final output ──────────────────────────────────────────────────
with open(in_file, newline="") as fi:
    reader   = csv.DictReader(fi, delimiter="\t")
    raw_cols = reader.fieldnames or []
    vf_cols  = [c for c in raw_cols if c not in set(FRONT_COLS)]
    # Column order:
    #   ID | typing | QC | interpretation | AMR genes | VF genes
    out_cols = FRONT_COLS + QC_COLS + INTERP_COLS + amr_genes + vf_cols
    rows     = list(reader)

with open(out_file, "w", newline="") as fo:
    writer = csv.DictWriter(fo, fieldnames=out_cols, delimiter="\t",
                            extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        sid = row["ID"]

        # QC columns
        q = qc.get(sid, {})
        for c in QC_COLS:
            row[c] = q.get(c, "-")

        # interpretation columns
        row.update(interpret(row))

        # AMR columns — one per unique gene, "-" if absent
        sample_amr = amr_data.get(sid, {})
        for g in amr_genes:
            row[g] = sample_amr.get(g, "-")

        writer.writerow(row)

PYEOF

########################################
# DONE
########################################

TIME_END=$(date +%s)
ELAPSED=$(( TIME_END - TIME_START ))
HOURS=$(( ELAPSED / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

echo ""
echo "══════════════════════════════════════"
echo "DONE → $OUTDIR/final_los_hs_gbs.tsv"
echo "  Samples processed : $N_SAMPLES"
N_DONE_FINAL=$(find "$OUTDIR/.done" -maxdepth 1 -name "*.done" -type f | wc -l)
echo "  Samples tracked   : $N_DONE_FINAL / $N_SAMPLES"
echo "  QC report         : $OUTDIR/qc_report.tsv"
echo "  AMR combined      : ${AMR_COMBINED:-N/A}"
echo "  Started           : $(date -d "@$TIME_START" '+%Y-%m-%d %H:%M:%S' 2>/dev/null || date -r "$TIME_START" '+%Y-%m-%d %H:%M:%S')"
echo "  Finished          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Duration          : ${HOURS}h ${MINS}m ${SECS}s"
echo "══════════════════════════════════════"
echo ""
echo "Notes:"
echo "  LOS typing: pid*cov best hit (no length bias) + whitelist mosaics (I/J and I/S only)"
echo "  T4SS ABSENT + chromosome-only assembly = possible false negative"
echo "  MEOP ABSENT = gene not detected (phase-variable, may still express)"
echo "  AMR columns = acquired genes + point mutations (AMRFinder Plus)"
echo "  STRESS genes (arsenic etc) excluded from final table"
[ "$SKIP_HS" -eq 1 ] && echo "  HS typing was SKIPPED (--skip_hs) — HS_TYPE, HS_GBS_MARKER, COMBINED_GBS reflect no HS data"
echo ""
