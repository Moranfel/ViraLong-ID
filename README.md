# 🦠 ViraLong-ID v6.1

**Long-read viral identification, genome reconstruction, comparative genomics and phylogenetic analysis for multi-sample ONT datasets.**

**ViraLong-ID** is a long-read viral identification, genome reconstruction, comparative genomics, and phylogenetic analysis pipeline designed for multi-sample batches.

It processes ONT-style long-read sequencing data from raw FASTQ reads to target viral contigs, global sequence alignments, pairwise identity matrices, clustered heatmaps, and maximum-likelihood phylogenetic trees.

ViraLong-ID also provides optional **BEAST 2 preparation and execution** for time-scaled phylogenetic and phylogeographic analyses, including automatic generation of geographic maps and time-scaled trees.

## ✨ What's new in v6.1

ViraLong-ID v6.1: Larger text in trees and heatmaps.

### 🚀 Major additions

- Pairwise nucleotide identity analysis for all retained genomes.
- Separate pairwise identity analysis for assembled isolates only.
- Standard and similarity-clustered identity heatmaps.
- Publication-ready phylogenetic tree rendering.
- Optional phylogenetic tree containing assembled target genomes plus a selected reference genome.
- Configurable phylogram or cladogram tree rendering.
- Improved filtering of alignment positions based on sequence occupancy.
- Assembly-specific read selection and coverage control.
- Parallel processing of multiple samples.
- Expanded BEAST 2 workflow.
- Automatic generation of a basic runnable BEAST 2 XML when a custom XML is not supplied.
- Automatic generation of BEAST 2 phylogeographic maps.
- Interactive HTML phylogeographic visualization.
- Automatic rendering of time-scaled phylogenetic trees.

### 🕰️ Optional BEAST 2 tools

Required only when using `--run-beast2`:

```text
beast
treeannotator
```

---
## 🚀 Quick Start

```bash
python ViraLong-ID.py \
  --taxid 3433772 \
  --reads sample1.fastq.gz sample2.fastq.gz \
  --outdir /path/to/output \
  --refseq-virus-fasta /path/to/refseq_virus.fasta \
  --threads 32
```

The pipeline processes all supplied samples and performs the shared comparative and phylogenetic analyses after sample-level target identification.

---

## 📥 Required Inputs

| Argument | Description |
|---|---|
| `--taxid` | Target NCBI Taxonomy ID |
| `--reads` | One or more FASTQ or FASTQ.GZ input files |
| `--outdir` | Output directory |
| `--refseq-virus-fasta` | Local RefSeq virus FASTA used to build the BLAST database |

---

## 💻 Sample and Computational Options

| Argument | Default | Description |
|---|---:|---|
| `--threads` | `8` | Total number of computational threads |
| `--sample-parallel` | `1` | Number of samples processed concurrently |
| `--sample-names` | none | Optional sample names supplied in the same order as `--reads` |

When several samples are processed simultaneously, `--sample-parallel` can reduce total batch processing time. Available CPU and memory should be considered when increasing this value.

---

## ✨ Read QC and Assembly Options

| Argument | Default | Description |
|---|---:|---|
| `--min-q` | `15` | Minimum mean read quality used by fastplong |
| `--flye-mode` | `meta` | Flye mode: `normal` or `meta` |
| `--flye-iterations` | `1` | Number of Flye polishing iterations |
| `--assembly-min-q` | `20.0` | Minimum read quality for assembly-read selection |
| `--assembly-min-len` | `auto` | Minimum read length retained for assembly |
| `--assembly-max-len` | `auto` | Maximum read length retained for assembly |
| `--assembly-target-cov` | `300` | Target read coverage used for assembly |
| `--assembly-retry-all-qc` / `--no-assembly-retry-all-qc` | enabled | Retry assembly using all QC-passed reads when required |

The assembly-specific filtering step allows the pipeline to reduce excessive sequencing depth while preferentially retaining high-quality reads appropriate for viral genome reconstruction.

---

## 🎯 Target Contig Selection

| Argument | Default | Description |
|---|---:|---|
| `--min-pident` | `70.0` | Minimum BLAST nucleotide identity |
| `--min-qcov` | `40.0` | Minimum BLAST query coverage |
| `--min-contig-len-phylo` | `300` | Minimum target contig length retained for phylogenetic analysis |

Target viral contigs are identified using BLAST against the reference database associated with the selected viral TaxID.

---

## 🌳 Alignment and Phylogenetic Analysis

After sample-level target identification, retained viral sequences are combined with the downloaded reference genomes.

The global phylogenetic workflow uses:

```text
MAFFT
  |
  v
strand correction
  |
  v
alignment occupancy filtering
  |
  v
trimAl
  |
  v
IQ-TREE
```

### ⚙️ Main options

| Argument | Default | Description |
|---|---:|---|
| `--trimal-gap-threshold` | `0.8` | trimAl gap threshold |
| `--min-alignment-occupancy` | `0.7` | Minimum occupancy required for retained global alignment positions |
| `--min-assembled-core-occupancy` | `0.7` | Minimum occupancy for assembled-genome core alignment positions |
| `--mafft-adjust-direction` | `on` | MAFFT strand correction: `off`, `on`, or `accurate` |
| `--tree-render-mode` | `phylogram` | Tree representation: `phylogram` or `cladogram` |

`phylogram` preserves branch-length information in the rendered tree, whereas `cladogram` emphasizes tree topology.

---

## 🎨 Pairwise Identity Analysis

Unless disabled, ViraLong-ID calculates pairwise nucleotide sequence identity from the aligned genomes.

Two datasets are analyzed:

1. All retained genomes, including assembled sequences and reference genomes.
2. Assembled target isolates only.

For each dataset, ViraLong-ID generates:

- Pairwise identity matrix.
- Standard identity heatmap.
- Similarity-clustered identity heatmap.
- PDF output.
- PNG output.

Similarity clustering automatically reorders genomes according to their sequence similarity, facilitating visualization of genetic groups and closely related isolates.

### ⚙️ Options

| Argument | Default | Description |
|---|---:|---|
| `--identity-plot-min` | automatic | Minimum value used for the heatmap color scale |
| `--skip-identity-plot` | disabled | Skip pairwise identity matrices and heatmap generation |

Principal outputs include:

```text
pairwise_identity.tsv
pairwise_identity_heatmap.pdf
pairwise_identity_heatmap.png
pairwise_identity_heatmap_clustered.pdf
pairwise_identity_heatmap_clustered.png

pairwise_identity_assembled_only.tsv
pairwise_identity_assembled_only_heatmap.pdf
pairwise_identity_assembled_only_heatmap.png
pairwise_identity_assembled_only_heatmap_clustered.pdf
pairwise_identity_assembled_only_heatmap_clustered.png
```

---

## 🧬 Assembled Genomes + Reference Tree

ViraLong-ID can optionally generate an additional maximum-likelihood tree containing only:

```text
assembled target genomes
          +
one reference genome
```

Enable this analysis with:

```bash
--assembled-reference-tree
```

By default:

```bash
--assembled-reference-id auto
```

ViraLong-ID preferentially selects an appropriate `NC_` RefSeq accession when available.

A specific reference accession can instead be supplied:

```bash
--assembled-reference-tree \
--assembled-reference-id NC_XXXXXXXX.X
```

This analysis is useful when the main objective is to visualize relationships among newly reconstructed genomes without displaying the complete collection of public reference sequences.

---

## 📂 Output Structure

```text
output/
├── 00_logs/
├── 01_references/
├── 02_blast_database/
├── 03_samples/
├── 04_combined_target_contigs/
├── 05_phylogeny_alignment/
├── 06_pairwise_identity/
├── 07_phylogeny_tree/
├── 08_phylogeny_tree_assembled_plus_reference/
├── 09_report/
├── 10_beast2_preparation/
├── 11_beast2_run/
└── 12_tmp/
```

### 📁 Main folders

| Folder | Contents |
|---|---|
| `00_logs/` | Logs generated during pipeline execution |
| `01_references/` | NCBI target-virus reference genomes and metadata |
| `02_blast_database/` | Local BLAST database |
| `03_samples/` | Sample-specific QC, assembly, BLAST and target-contig results |
| `04_combined_target_contigs/` | Combined target contigs from all samples |
| `05_phylogeny_alignment/` | MAFFT and filtered/trimmed alignments |
| `06_pairwise_identity/` | Identity matrices and standard/clustered heatmaps |
| `07_phylogeny_tree/` | Global IQ-TREE phylogenetic analysis and rendered tree |
| `08_phylogeny_tree_assembled_plus_reference/` | Optional assembled-genomes + reference tree |
| `09_report/` | Per-sample and batch summaries |
| `10_beast2_preparation/` | Optional BEAST 2 input files and metadata templates |
| `11_beast2_run/` | Optional BEAST 2 results, MCC tree and visualizations |
| `12_tmp/` | Temporary working files |

---

## 🧱 Sample-Level Output

Each sample receives its own directory:

```text
03_samples/
└── sample_name/
    ├── 00_logs/
    ├── 01_reads_qc/
    ├── 02_reads_renamed/
    ├── 03_reads_for_assembly/
    ├── 04_assembly_flye/
    ├── 05_blast_identification/
    ├── 06_taxon_filtered_contigs/
    ├── 07_report/
    └── 08_tmp/
```

This structure separates sample-specific processing from batch-level comparative analyses.

---

# 🕰️ Optional BEAST 2 Analysis

> **BEAST 2 is optional.** The core ViraLong-ID workflow runs independently of BEAST 2.

The recommended BEAST 2 workflow is:

```text
1. Prepare BEAST 2 files
2. Review and complete sampling dates
3. Review and complete geographic coordinates
4. Rebuild BEAST 2 metadata
5. Run BEAST 2
6. Summarize posterior trees
7. Generate time-scaled and phylogeographic visualizations
```

---

## 🧪 Stage 1 - Prepare BEAST 2 Files

Add:

```bash
--prepare-beast2
```

Example:

```bash
python ViraLong-ID.py \
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

### 📄 Main files

| File | Description |
|---|---|
| `alignment_beast2_safe_ids.fasta` | Trimmed alignment using BEAST-compatible sequence identifiers |
| `alignment_beast2_safe_ids.nexus` | Same alignment in NEXUS format |
| `metadata_beast2.tsv` | Automatically extracted sequence metadata |
| `tip_dates_beast2.tsv` | Sampling dates formatted for BEAST |
| `manual_dates_template.tsv` | Editable table for missing or corrected sampling dates |
| `traits_beauti.tsv` | Discrete traits such as country, host, region or sample type |
| `map_locations_coordinates_template.tsv` | Editable geographic location and coordinate table |
| `sequence_coordinates_template.tsv` | Sequence-to-coordinate mapping |
| `sequence_id_map.tsv` | Original sequence headers mapped to BEAST-safe identifiers |
| `CYVCV_BEAST2_template.xml` | BEAST 2 template/documentation file |
| `README_BEAST2_preparacion.md` | Notes describing the generated BEAST 2 files |

---

## 📅 Stage 2 - Complete Sampling Metadata

Time-scaled phylogenetic reconstruction requires sampling dates.

Edit:

```text
manual_dates_template.tsv
```

Exact dates should preferably use:

```text
YYYY-MM-DD
```

When exact sampling dates are unavailable, the corresponding year can be supplied when appropriate for the intended analysis.

---

## 🗺️ Stage 3 - Complete Geographic Metadata

Phylogeographic visualization requires geographic coordinates.

Edit:

```text
map_locations_coordinates_template.tsv
```

Provide latitude and longitude for the relevant sampling locations.

Coordinates are required for continuous geographic visualization.

Discrete geographic traits such as country or region may also be represented in the metadata.

---

## 🔁 Stage 4 - Rebuild BEAST 2 Tables

After editing the metadata files, run:

```bash
python ViraLong-ID.py \
  --taxid 3433772 \
  --reads sample1.fastq.gz sample2.fastq.gz \
  --outdir /path/to/output \
  --refseq-virus-fasta /path/to/refseq_virus.fasta \
  --prepare-beast2 \
  --beast2-manual-dates /path/to/manual_dates_template.tsv \
  --beast2-coordinates /path/to/map_locations_coordinates_template.tsv
```

The BEAST 2 preparation files are regenerated using the supplied dates and coordinates.

---

## 🚀 Stage 5 - Run BEAST 2

Enable BEAST execution with:

```bash
--run-beast2
```

When `--run-beast2` is used, ViraLong-ID requires completed files supplied through:

```text
--beast2-manual-dates
--beast2-coordinates
```

### 🤖 Automatic XML generation

If `--beast2-xml` is omitted, ViraLong-ID automatically creates a basic runnable XML:

```text
CYVCV_BEAST2.xml
```

Example:

```bash
python ViraLong-ID.py \
  --taxid 3433772 \
  --reads sample1.fastq.gz sample2.fastq.gz \
  --outdir /path/to/output \
  --refseq-virus-fasta /path/to/refseq_virus.fasta \
  --threads 32 \
  --prepare-beast2 \
  --beast2-manual-dates /path/to/manual_dates_template.tsv \
  --beast2-coordinates /path/to/map_locations_coordinates_template.tsv \
  --run-beast2
```

### 🛠️ Custom BEAST XML

For analyses requiring a specifically configured evolutionary model, molecular clock, tree prior or other BEAST settings, a custom XML generated or reviewed externally can be supplied:

```bash
--beast2-xml /path/to/final_BEAST2.xml
```

For example:

```bash
python ViraLong-ID.py \
  --taxid 3433772 \
  --reads sample1.fastq.gz sample2.fastq.gz \
  --outdir /path/to/output \
  --refseq-virus-fasta /path/to/refseq_virus.fasta \
  --threads 32 \
  --prepare-beast2 \
  --beast2-manual-dates /path/to/manual_dates_template.tsv \
  --beast2-coordinates /path/to/map_locations_coordinates_template.tsv \
  --run-beast2 \
  --beast2-xml /path/to/final_BEAST2.xml
```

---

## ⏱️ BEAST 2 MCMC Options

| Argument | Default | Description |
|---|---:|---|
| `--beast2-chain-length` | `10000000` | MCMC chain length used for automatically generated XML |
| `--beast2-log-every` | `10000` | Logging interval |
| `--beast2-burnin` | `10` | Percentage burn-in used by TreeAnnotator |
| `--beast2-xml` | automatic | Optional custom BEAST 2 XML |

The chain length and sampling frequency should be adjusted according to dataset complexity and convergence diagnostics.

---

## 📊 BEAST 2 Output

BEAST 2 results are written to:

```text
11_beast2_run/
```

Principal outputs include:

```text
CYVCV_BEAST2.MCC.tree

CYVCV_BEAST2_phylogeography_map.pdf
CYVCV_BEAST2_phylogeography_map.png
CYVCV_BEAST2_phylogeography_map.html

CYVCV_BEAST2_time_tree.pdf
CYVCV_BEAST2_time_tree.png

CYVCV_BEAST2_phylogeography_edges.tsv
```

### 🌳 MCC tree

```text
CYVCV_BEAST2.MCC.tree
```

contains the maximum clade credibility tree summarized with TreeAnnotator.

### 🗺️ Phylogeographic map

ViraLong-ID automatically generates static phylogeographic maps in:

```text
PDF
PNG
```

and an interactive visualization in:

```text
HTML
```

### ⏳ Time-scaled tree

The BEAST MCC tree is also rendered as a time-scaled phylogenetic tree in:

```text
PDF
PNG
```

### 🔗 Phylogeographic edges

```text
CYVCV_BEAST2_phylogeography_edges.tsv
```

contains the inferred geographic connections extracted for visualization and downstream analysis.

---

## ♻️ Resume Behavior

ViraLong-ID is designed to resume partially completed analyses.

Completed stages are skipped when their expected output files are already present.

This allows interrupted or computationally expensive runs to continue without repeating completed analyses.

When edited BEAST 2 date or coordinate files are explicitly supplied, the corresponding BEAST 2 preparation stage is regenerated so that the updated metadata are incorporated.

---

## 🧪 Example - Standard Multi-Sample Analysis

```bash
python ViraLong-ID.py \
  --taxid 3433772 \
  --reads /path/to/*.fastq.gz \
  --threads 32 \
  --sample-parallel 4 \
  --outdir /path/to/analysis \
  --refseq-virus-fasta /path/to/sequences.fasta
```

---

## 🌳 Example - Assembled Genomes + Reference Tree

```bash
python ViraLong-ID.py \
  --taxid 3433772 \
  --reads /path/to/*.fastq.gz \
  --threads 32 \
  --outdir /path/to/analysis \
  --refseq-virus-fasta /path/to/sequences.fasta \
  --assembled-reference-tree
```

To select a particular reference:

```bash
--assembled-reference-tree \
--assembled-reference-id NC_XXXXXXXX.X
```

---

## 🍊 Example - CYVCV Analysis

```bash
python ViraLong-ID.py \
  --taxid 1214459 \
  --reads /path/to/*.fastq.gz \
  --threads 32 \
  --outdir /path/to/CYVCV_output \
  --trimal-gap-threshold 0.95 \
  --refseq-virus-fasta /path/to/sequences.fasta \
  --assembled-reference-tree
```

---

## 🕰️ Example - CYVCV With BEAST 2 Preparation

```bash
python ViraLong-ID.py \
  --taxid 1214459 \
  --reads /path/to/*.fastq.gz \
  --threads 32 \
  --outdir /path/to/CYVCV_output \
  --trimal-gap-threshold 0.95 \
  --refseq-virus-fasta /path/to/sequences.fasta \
  --prepare-beast2
```

---

## 📑 Reporting Checklist

For reproducible reporting of ViraLong-ID analyses, document at least:

- ViraLong-ID version.
- Target NCBI TaxID.
- Sequencing platform and library preparation.
- Read QC threshold.
- Assembly read-selection parameters.
- Flye mode and polishing iterations.
- BLAST identity threshold.
- BLAST query-coverage threshold.
- Minimum target-contig length.
- MAFFT strand-correction setting.
- Alignment occupancy threshold.
- trimAl gap threshold.
- IQ-TREE model and relevant tree statistics.
- Pairwise identity analysis settings.
- Reference accession used for the assembled + reference tree, when applicable.

For BEAST 2 analyses, additionally report:

- Sampling-date information.
- Evolutionary substitution model.
- Molecular clock model.
- Tree prior.
- MCMC chain length.
- Logging interval.
- Burn-in.
- Convergence diagnostics.
- Effective sample sizes (ESS).
- Geographic traits or coordinates used for phylogeographic reconstruction.

---

## 🔬 Reproducibility

ViraLong-ID records intermediate and final outputs in structured directories to facilitate inspection, troubleshooting and reproducibility.

For scientific analyses, the exact ViraLong-ID version, Conda environment, command-line parameters and reference database version should be retained together with the sequencing data.

Raw sequencing reads and reconstructed viral genomes should be deposited in appropriate public repositories when results are published.

---

## 📚 Citation

If you use ViraLong-ID in your research, please cite the software as: Morán, F. (2026). *ViraLong-ID: A long-read viral identification and phylogeny pipeline*

A formal citation and DOI will be provided in a future release.

---

## 🔬 Workflow

```text
Raw FASTQ / FASTQ.GZ
        |
        v
    fastplong
        |
        v
Assembly read selection
        |
        v
      Flye
        |
        v
 BLAST target identification
        |
        v
 Target viral contigs
        |
        v
Combined target dataset
        |
        v
      MAFFT
        |
        v
     trimAl
        |
        +---------------------------+
        |                           |
        v                           v
 Pairwise identity              IQ-TREE
        |                           |
        v                           v
 Identity matrices          ML phylogenetic tree
        |
        +--> Standard heatmaps
        |
        +--> Clustered heatmaps
        |
        +--> Assembled-only analysis

Optional extensions:
        |
        +--> Assembled genomes + reference ML tree
        |
        +--> BEAST 2 preparation
                 |
                 v
           Sampling metadata
                 |
                 v
             BEAST 2
                 |
        +--------+---------+
        |                  |
        v                  v
   MCC time tree     Phylogeographic maps
```

---

See the repository `LICENSE` file for the complete license terms.
