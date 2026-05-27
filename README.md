# 🦠 ViraLong-ID v5.7

**ViraLong-ID** is a long-read viral identification and phylogeny pipeline for multi-sample batches. It takes ONT-style long-read data from raw reads to viral contigs, global alignments, identity heatmaps, and maximum-likelihood trees.

✨ **New in v5.7:** optional **BEAST 2 preparation and execution** for time-scaled and map-ready phylogeographic analyses.

---

## 🚀 What It Does

| Step | Output |
|---|---|
| 🧬 Target reference download | Complete viral genomes for the selected NCBI TaxID |
| ✨ Read QC | Filtered long reads |
| 🧱 Assembly | Sample-level assemblies with Flye |
| 🎯 Viral detection | BLAST-based target contig identification |
| 🌳 Global phylogeny | MAFFT alignment, trimAl filtering, IQ-TREE tree |
| 🎨 Identity analysis | Pairwise identity matrices and heatmaps |
| 🕰️ BEAST 2 preparation | Dates, metadata, traits, coordinates, BEAST-safe alignment |
| 🚀 BEAST 2 execution | Optional BEAST run plus MCC tree summarization |

---

## 🧰 Requirements

### Core tools

```text
datasets
dataformat
unzip
fastplong
flye
makeblastdb
blastn
mafft
trimal
iqtree
```

### Optional BEAST 2 tools

```text
beast
treeannotator
```

### Python dependency

```text
biopython
```

---

## ⚡ Quick Start

```bash
python ViraLong-ID_v5.7.py \
  --taxid 3433772 \
  --reads sample1.fastq.gz sample2.fastq.gz \
  --outdir /path/to/output \
  --refseq-virus-fasta /path/to/refseq_virus.fasta \
  --threads 32
```

---

## 📥 Required Inputs

| Argument | Description |
|---|---|
| `--taxid` | Target NCBI Taxonomy ID |
| `--reads` | One or more input FASTQ/FASTQ.GZ files |
| `--outdir` | Output directory |
| `--refseq-virus-fasta` | Local RefSeq virus FASTA used to build the BLAST database |

---

## ⚙️ Useful Options

| Argument | Default | Description |
|---|---:|---|
| `--threads` | `8` | Number of threads |
| `--sample-names` | none | Optional sample names, same order as `--reads` |
| `--min-q` | `15` | Minimum mean read quality for `fastplong` |
| `--flye-mode` | `meta` | Flye mode: `normal` or `meta` |
| `--flye-iterations` | `1` | Number of Flye polishing iterations |
| `--min-pident` | `70.0` | Minimum BLAST identity for target contig selection |
| `--min-qcov` | `40.0` | Minimum BLAST query coverage for target contig selection |
| `--min-contig-len-phylo` | `300` | Minimum contig length retained for phylogeny |
| `--trimal-gap-threshold` | `0.8` | trimAl gap threshold |
| `--mafft-adjust-direction` | `on` | MAFFT strand correction: `off`, `on`, or `accurate` |

---

## 📂 Output Structure

```text
output/
├── 00_logs/
├── 01_references/
├── 05_blast_database/
├── 06_combined_target_contigs/
├── 07_phylogeny_alignment/
├── 07b_pairwise_identity/
├── 08_phylogeny_tree/
├── 09_report/
├── 10_beast2_preparation/
├── 11_beast2_run/
├── samples/
└── tmp/
```

| Folder | Contents |
|---|---|
| `00_logs/` | Hidden output from external tools |
| `01_references/` | NCBI target-virus references and metadata |
| `06_combined_target_contigs/` | Combined target contigs from all samples |
| `07_phylogeny_alignment/` | MAFFT alignments and trimmed alignment |
| `07b_pairwise_identity/` | Pairwise identity matrices and heatmaps |
| `08_phylogeny_tree/` | IQ-TREE tree files and rendered tree PDF |
| `09_report/` | Per-sample and batch summaries |
| `10_beast2_preparation/` | Optional BEAST 2 input files and editable templates |
| `11_beast2_run/` | Optional BEAST 2 output and MCC tree |

---

## 🕰️ BEAST 2 Workflow

BEAST 2 support is optional and split into clear stages.

```text
1. Prepare BEAST 2 files
2. Fill missing dates and coordinates
3. Rebuild BEAST 2 tables
4. Create/review final XML in BEAUti
5. Run BEAST 2
6. Summarize the MCC tree
7. Use BEAST output for map visualization
```

This prevents long MCMC runs from starting before the temporal and geographic metadata are ready.

---

## 🧪 Stage 1: Prepare BEAST 2 Files

Add:

```bash
--prepare-beast2
```

Example:

```bash
python ViraLong-ID_v5.7.py \
  --taxid 3433772 \
  --reads sample1.fastq.gz sample2.fastq.gz \
  --outdir /path/to/output \
  --refseq-virus-fasta /path/to/refseq_virus.fasta \
  --threads 32 \
  --prepare-beast2
```

This creates:

```text
10_beast2_preparation/
```

### Files created

| File | Description |
|---|---|
| `alignment_beast2_safe_ids.fasta` | Trimmed alignment with BEAST-compatible sequence IDs |
| `alignment_beast2_safe_ids.nexus` | Same alignment in NEXUS format |
| `metadata_beast2.tsv` | Automatically extracted metadata per sequence |
| `tip_dates_beast2.tsv` | Sampling dates formatted for BEAST/BEAUti |
| `manual_dates_template.tsv` | Editable file for missing sampling dates |
| `traits_beauti.tsv` | Optional discrete traits: country, host, region, sample type |
| `map_locations_coordinates_template.tsv` | Editable location table for latitude/longitude |
| `sequence_coordinates_template.tsv` | Per-sequence coordinate mapping table |
| `sequence_id_map.tsv` | Original headers mapped to BEAST-safe IDs |
| `CYVCV_BEAST2_template.xml` | Non-runnable XML template/documentation file |
| `README_BEAST2_preparacion.md` | BEAST 2 preparation notes |

---

## ✍️ Stage 2: Fill Manual Metadata

BEAST 2 needs sampling time. Map-ready phylogeography also needs coordinates.

Edit these files:

```text
manual_dates_template.tsv
map_locations_coordinates_template.tsv
```

Use exact dates when possible:

```text
YYYY-MM-DD
```

If exact dates are unavailable, use the `year` column.

---

## 🔁 Stage 3: Rebuild BEAST 2 Tables

After editing the TSV files, re-run with:

```bash
--prepare-beast2 \
--beast2-manual-dates /path/to/manual_dates_template.tsv \
--beast2-coordinates /path/to/map_locations_coordinates_template.tsv
```

This incorporates the edited dates and coordinates into the BEAST 2 tables.

---

## 🚀 Stage 4: Run BEAST 2

Before running BEAST 2, create and review the final runnable XML in **BEAUti** using the prepared alignment, dates, traits, and coordinates.

Then run:

```bash
--prepare-beast2 \
--beast2-manual-dates /path/to/manual_dates_template.tsv \
--beast2-coordinates /path/to/map_locations_coordinates_template.tsv \
--run-beast2 \
--beast2-xml /path/to/CYVCV_BEAST2.xml
```

Optional burn-in for TreeAnnotator:

```bash
--beast2-burnin 10
```

The BEAST 2 run output is written to:

```text
11_beast2_run/
```

The summarized MCC tree is:

```text
CYVCV_BEAST2.MCC.tree
```

---

## 🗺️ Notes on Maps

The `10_beast2_preparation/` folder does **not** generate a final map by itself.

To produce a map, you need:

| Requirement | Why it matters |
|---|---|
| ✅ Complete sampling dates | Required for time-scaled BEAST analysis |
| ✅ Latitude/longitude coordinates | Required for continuous geographic mapping |
| ✅ Final BEAST 2 XML | Defines the model, clock, priors, and traits |
| ✅ Completed BEAST 2 run | Produces posterior trees/logs |
| ✅ Visualization step | Converts BEAST output into a geographic map |

Coordinates are required for continuous maps. Discrete traits such as country or region can also be used, but plotting them still requires representative coordinates.

---

## ♻️ Resume Behavior

ViraLong-ID skips completed steps when their expected outputs already exist.

When edited BEAST 2 date or coordinate files are provided, the BEAST 2 preparation step is regenerated so the updated metadata are incorporated.

---

## 🧬 Example: CYVCV Batch With BEAST 2 Preparation

```bash
python ViraLong-ID_v5.7.py \
  --taxid 3433772 \
  --reads /path/to/*.fastq.gz \
  --threads 32 \
  --outdir /path/to/CYVCV_output \
  --trimal-gap-threshold 0.95 \
  --refseq-virus-fasta /path/to/sequences.fasta \
  --prepare-beast2
```

---

## 📑 Reporting Checklist

When reporting results, document:

- 🧾 ViraLong-ID version
- 🧬 Target TaxID
- ✨ Read QC thresholds
- 🧱 Assembly parameters
- 🎯 BLAST identity/query coverage thresholds
- 🌳 Alignment and trimming parameters
- 📊 IQ-TREE model/results
- 🕰️ BEAST 2 model, clock, tree prior, chain length, burn-in, and convergence diagnostics if BEAST 2 was used

