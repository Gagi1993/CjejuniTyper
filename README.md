# CjejuniTyper

**CjejuniTyper** is a bioinformatics pipeline for multi-layer genotypic characterisation of *Campylobacter jejuni* genome assemblies.

It performs in silico Penner HS capsular serotyping, LOS class typing, MLST typing, virulence factor detection, AMR detection, assembly QC, and GBS risk prediction from whole-genome sequencing assemblies.

---

## 🔬 Overview

CjejuniTyper analyses *C. jejuni* assemblies using BLAST-based marker detection, full-locus fallback typing, in silico PCR confirmation, and rule-based interpretation.

The pipeline is intended for reproducible genomic characterisation of *Campylobacter jejuni* from draft genome assemblies and can be integrated into surveillance workflows.

---

## ⚙️ Features

* Penner HS capsular serotyping
* Two-step HS typing: specific marker BLAST followed by full-locus BLAST fallback
* Subtype bracketing for selected HS complexes
* LOS class typing using best-hit BLAST logic
* MLST typing using the 7-gene *Campylobacter* scheme
* Virulence factor detection
* `cstII` allele typing
* GBS risk prediction based on LOS/HS markers and ganglioside mimic status
* MEOP, T4SS, CDT, capsule and LOS-core status reporting
* AMR detection using AMRFinder Plus
* Assembly QC: genome size, GC content and contig count
* Resume mode and parallel processing
* Optional skip flags for HS, QC, MLST, LOS, VF and AMR modules

---

## ⚙️ Requirements

* Linux-based operating system
* Conda or Mamba
* Python 3
* BLAST
* AMRFinder Plus
* Exonerate / `ipcress`
* SeqKit

The automatic installer creates a conda/mamba environment named:

```bash
cjejuni_typer
```

---

## 🚀 Installation

### Automatic install

```bash
bash CjejuniTyper.sh -I
```

After installation, update the AMRFinder database:

```bash
bash CjejuniTyper.sh -amrfinder_db
```

### Update installed environment

```bash
bash CjejuniTyper.sh -U
```

---

## 🧬 Database setup

CjejuniTyper expects the required databases to be available in the `databases/` folder.

Default database paths used by the script:

```text
databases/new_hs_db/blast
databases/hs_specific_aa_per_hs/hs/hs_aa
databases/los_db_new
databases/los_db_specific_markers/LOS_B_windows/LOS_B_WINDOWS_DB
databases/vfdb_new
databases/vfdb/orf11/orf11_db
databases/mlst/
```

### T4SS database setup

Run once, if using T4SS detection:

```bash
grep -A1 "virB\|virD4" campy_gene.fas > virB_virD4.fasta
mkdir -p databases/vfdb/blast/T4SS
cp virB_virD4.fasta databases/vfdb/blast/T4SS/T4SS.fasta

makeblastdb \
  -in databases/vfdb/blast/T4SS/T4SS.fasta \
  -dbtype nucl \
  -out databases/vfdb/blast/T4SS/T4SS
```

---

## ❓ Help

Display all available options:

```bash
bash CjejuniTyper.sh -h
```

---

## 🧪 Input

Input directory containing genome assemblies.

Accepted extensions:

```text
.fa
.fasta
.fna
```

---

## 📦 Output

Main output file:

```text
final_los_hs_gbs.tsv
```

Main output columns include:

* genome ID
* LOS class
* sialylation status
* HS type
* GBS risk
* genome size
* contig count
* GC percentage
* QC status and warnings
* MLST result
* `cstII` function
* ganglioside mimic prediction
* combined GBS score
* serum resistance status
* MEOP status
* T4SS status
* CDT status
* capsule and LOS-core markers
* AMR genes and mutations
* raw virulence factor gene calls

Intermediate output files include module-specific BLAST results, QC reports, MLST output, HS typing tables, LOS typing tables, VF results and AMRFinder outputs.

---

## ▶️ Example usage

Run the full pipeline:

```bash
bash CjejuniTyper.sh -i assemblies/ -o results/
```

Run with 16 total threads and 4 parallel jobs:

```bash
bash CjejuniTyper.sh -i assemblies/ -o results/ -t 16 -j 4
```

Resume a previous run:

```bash
bash CjejuniTyper.sh -i assemblies/ -o results/ -r
```

Run LOS, VF and AMR only:

```bash
bash CjejuniTyper.sh -i assemblies/ -o results/ --skip_hs
```

---

## ▶️ Usage

```bash
bash CjejuniTyper.sh \
  -i assemblies/ \
  -o results/ \
  -t 16 \
  -j 4
```

### Required parameters

* `-i` : input directory containing `.fa`, `.fasta` or `.fna` assemblies
* `-o` : output directory

### Optional parameters

* `-hs` : HS full-locus BLAST database directory
* `-hs_aa` : HS specific marker gene BLAST database
* `-los` : LOS class BLAST database directory
* `-los_b_windows` : LOS class B window database directory
* `-vf` : virulence factor BLAST database directory
* `-t` : total threads, default `32`
* `-j` : parallel jobs, default `1`
* `-r` : resume previous run and skip completed samples
* `--skip_hs` : skip HS Penner typing
* `--skip_qc` : skip assembly QC
* `--skip_mlst` : skip MLST typing
* `--skip_los` : skip LOS typing
* `--skip_vf` : skip virulence factor typing
* `--skip_amr` : skip AMRFinder detection
* `-I` : install conda environment
* `-U` : update conda environment
* `-amrfinder_db` : force-update AMRFinder database
* `-h` : show help

---

## 🧬 HS Penner typing

CjejuniTyper reports HS serotype complexes as complexes.

Examples:

```text
HS4c
HS8_17c
HS8_17c(HS8*)
```

The asterisk `*` means the subtype is predicted from genomic data and is not confirmed by serology.

### HS complexes

| Complex | Member serotypes | Subtype/bracket method |
|---|---|---|
| HS4c | HS4, HS13, HS16, HS43, HS50, HS62, HS63, HS64, HS65 | Mu_HS4B PCR / BLAST winner |
| HS6_7c | HS6, HS7 | no subtype bracket |
| HS8_17c | HS8, HS17 | AA blastx |
| HS5_31c | HS5, HS31 | Mu_HS5 in silico PCR |
| HS23_36c | HS23, HS36 | poly-G/C tract |
| HS1c | HS1, HS44 | BLAST winner |
| HS15c | HS15, HS31, HS58 | Mu_HS58 in silico PCR |
| HS45c | HS45, HS32, HS60 | BLAST winner |

Standalone HS types include:

```text
HS2, HS3, HS9, HS10, HS11, HS12, HS18, HS19, HS21, HS22, HS27,
HS29, HS33_35, HS37, HS38, HS40, HS41, HS42, HS52, HS53, HS55, HS57
```

---

## 🧫 LOS typing

LOS class typing is BLAST-based and uses best-hit selection.

The script uses strict, mid and loose thresholds based on identity, coverage and bitscore. True mosaic calls such as `I/J` and `I/S` are retained when both hits meet the required coverage.

---

## 💊 AMR detection

AMR detection is performed with AMRFinder Plus.

The pipeline detects acquired resistance genes and point mutations, including markers such as:

* `gyrA`
* 23S rRNA mutations
* `tet(O)`
* `blaOXA`

STRESS genes are excluded from the final output.

Before AMR analysis, update the AMRFinder database:

```bash
bash CjejuniTyper.sh -amrfinder_db
```

---

## 📁 Repository structure

```text
CjejuniTyper/
├── CjejuniTyper.sh
├── README.md
├── LICENSE
├── .gitignore
├── databases/
├── example_data/
└── results/
```

---

## 📊 Example

```bash
bash CjejuniTyper.sh -i example_data/ -o results/
```

---

## 📖 Citation

If you use this tool, please cite:

Jurić D. et al.  
*In silico typing and genomic characterisation of Campylobacter jejuni from whole-genome sequencing data.*  
(Manuscript in preparation)

---

## 📜 License

This project is licensed under the MIT License.

---

## 👤 Author

Dragan Jurić  
Croatian Institute of Public Health  
Zagreb, Croatia

---

## ⚠️ Notes

* Designed for genomic characterisation of *Campylobacter jejuni* assemblies.
* HS subtype brackets are best-effort genomic predictions, not serological confirmation.
* HS23/36c subtype prediction depends on accurate poly-G/C tract counting and may be unreliable in short-read draft assemblies.
* AMR detection requires an updated AMRFinder database.
* Performance depends on genome quality, database completeness and assembly continuity.
