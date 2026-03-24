# ViraLong-ID 🦠🧬

![alt text](<Screenshot 2026-03-16 at 6.57.19 PM.png>)

**Long-read viral identification, identity heatmap generation, and global phylogeny pipeline**

ViraLong-ID is a command-line workflow for **viral identification from long-read sequencing data** followed by **global phylogenetic reconstruction** and **pairwise genome identity visualization**. It is designed for **batch analysis of multiple samples in a single run**, while reusing shared steps such as reference download, reference preparation, and BLAST database construction.

## Why use ViraLong-ID? 🚀

- 🦠 Detect viral contigs from long-read data
- 📦 Process multiple samples in one run
- 🧬 Reuse shared references for the whole batch
- 🔎 Build the BLAST database only once
- 🔄 Correct contig orientation automatically before final alignment
- ✂️ Trim alignments to retain well-supported homologous regions
- 🌳 Infer one global IQ-TREE phylogeny
- 🎨 Generate a pairwise identity heatmap ordered according to the final tree
- 📄 Export the final tree as PDF with assembled isolates highlighted
- 📁 Keep logs and outputs clearly organized

## What does the pipeline do? ⚙️

### Shared steps for the whole batch

1. 📥 Download complete viral genomes for the target NCBI TaxID
2. 🧾 Build metadata tables
3. 🏷️ Rename complete target references
4. 🗃️ Build a shared local BLAST database from your RefSeq virus FASTA

### Per-sample steps

1. ✨ Filter reads with `fastplong`
2. ✂️ Shorten read headers for Flye
3. 📏 Preselect reads for assembly
4. 🧱 Assemble reads with `Flye`
5. 🔎 BLAST contigs against the shared virus database
6. 🎯 Retain contigs matching the selected target taxon

### Global analysis steps

1. 📦 Combine retained contigs from all samples
2. 🔄 Correct sequence orientation with `MAFFT`
3. 🧬 Build a global alignment with `MAFFT`
4. ✂️ Trim the alignment with `trimAl`
5. 🌳 Infer a maximum-likelihood tree with `IQ-TREE`
6. 🎨 Build a pairwise identity matrix and heatmap from the trimmed alignment
7. 📄 Render the final phylogeny as PDF
8. 📝 Write per-sample and batch summary reports

## What is new in the current version? ✨

The latest version of ViraLong-ID includes:

- robust strand correction before the final fragment alignment
- trimmed phylogenetic alignments using `trimAl`
- pairwise identity heatmaps generated from the trimmed alignment
- heatmap ordering based on the final IQ-TREE topology
- PDF tree rendering with **assembled isolates automatically highlighted**

## Requirements 🛠️

### Python packages

- `python >= 3.9`
- `biopython`
- `matplotlib`

### External tools

Make sure these executables are available in your environment:

- `datasets`
- `dataformat`
- `unzip`
- `fastplong`
- `flye`
- `makeblastdb`
- `blastn`
- `mafft`
- `trimal`
- `iqtree`

## Input 📥

You need:

- 🧬 A target NCBI Taxonomy ID
- 📂 One or more input `FASTQ` or `FASTQ.GZ` files
- 🗃️ A local RefSeq virus FASTA file for BLAST database creation
- 📁 An output directory

Optional arguments let you define:

- 🏷️ Sample names
- 🧵 Number of threads
- ✨ QC thresholds
- 🎯 BLAST filtering thresholds
- 📏 Assembly read selection parameters
- ✂️ trimAl stringency
- 🔄 MAFFT strand correction mode
- 🎨 heatmap color scale minimum

## Example command 💻

```bash
python ViraLong-ID_v4.8.py \
  --taxid 1214459 \
  --reads /path/to/Selection.fastq.gz /path/to/Test.fastq.gz \
  --sample-names Selection Test \
  --outdir resultados \
  --refseq-virus-fasta /path/to/viral.1.1.genomic.fna \
  --threads 16 \
  --mafft-adjust-direction on \
  --trimal-gap-threshold 0.8 \
  --identity-plot-min 96.5
```

## Main parameters 🧾

- `--taxid`: target NCBI Taxonomy ID
- `--reads`: one or more input `FASTQ` / `FASTQ.GZ` files
- `--sample-names`: optional sample names in the same order as `--reads`
- `--outdir`: output directory
- `--refseq-virus-fasta`: local FASTA used to build the shared BLAST database
- `--threads`: number of threads
- `--min-q`: minimum mean read quality for `fastplong`
- `--flye-mode`: `Flye` mode (`normal` or `meta`)
- `--min-pident`: minimum BLAST percent identity for target contig retention
- `--min-qcov`: minimum BLAST query coverage for target contig retention
- `--min-contig-len-phylo`: minimum contig length retained for phylogeny
- `--assembly-min-q`: minimum mean read quality retained for assembly
- `--assembly-min-len`: minimum read length for assembly preselection or `auto`
- `--assembly-max-len`: maximum read length for assembly preselection or `auto`
- `--assembly-target-cov`: target effective coverage used to cap Flye input
- `--trimal-gap-threshold`: trimAl gap threshold
- `--mafft-adjust-direction`: MAFFT strand correction mode (`off`, `on`, `accurate`)
- `--identity-plot-min`: lower bound for heatmap color scale

## Output structure 📂

The batch run creates **shared outputs** at the top level and **sample-specific outputs** inside `samples/`.

```text
outdir/
├── 00_logs/
├── 01_references/
├── 05_blast_database/
├── 06_combined_target_contigs/
├── 07_phylogeny_alignment/
├── 07b_pairwise_identity/
├── 08_phylogeny_tree/
├── 09_report/
└── samples/
    ├── Sample1/
    └── Sample2/
```

### Shared outputs 🌍

- `01_references/` → downloaded target references and metadata
- `05_blast_database/` → shared BLAST database files
- `06_combined_target_contigs/` → combined retained contigs
- `07_phylogeny_alignment/alignment_mafft.fasta` → global alignment
- `07_phylogeny_alignment/alignment_mafft.trimmed.fasta` → trimmed alignment
- `07b_pairwise_identity/pairwise_identity.tsv` → pairwise identity matrix
- `07b_pairwise_identity/pairwise_identity_heatmap.pdf` → identity heatmap in PDF
- `07b_pairwise_identity/pairwise_identity_heatmap.png` → identity heatmap in PNG
- `08_phylogeny_tree/alignment_mafft.trimmed.treefile` → final Newick tree
- `08_phylogeny_tree/alignment_mafft.trimmed.tree.pdf` → rendered phylogeny PDF
- `09_report/multi_sample_summary.tsv` → batch summary table
- `09_report/multi_sample_summary.txt` → batch summary report

### Per-sample outputs 🧪

Each sample folder contains:

- filtered reads
- renamed reads
- assembly input subset
- `Flye` assembly output
- BLAST identification results
- retained target contigs
- per-sample report

## Visual outputs 🎨

ViraLong-ID produces two publication-friendly visual outputs:

### 🌳 Phylogenetic tree PDF

- based on the trimmed alignment
- inferred with `IQ-TREE`
- assembled isolates highlighted automatically
- reference genomes shown separately in the legend

### 🟨 Pairwise identity heatmap

- built from the trimmed alignment
- pairwise identity values calculated between all sequences
- ordered according to the final tree topology
- exported as both `PDF` and `PNG`

## Terminal experience 🖥️

ViraLong-ID includes a styled terminal interface with:

- 🦠 startup logo
- 📊 progress bars
- 🎯 clear section blocks
- 📁 compact summaries of outputs and logs
- `-h` help message with styled descriptions

## Important notes ⚠️

- Shared steps are run **only once per batch**
- The final phylogeny is **global**, not per sample
- The heatmap is generated from the **trimmed alignment**
- The heatmap is ordered according to the **final IQ-TREE topology**
- If no target contigs are retained across the batch, the phylogeny cannot be built
- External tool output is hidden from the terminal and saved in log files
- NCBI download steps may fail temporarily if the remote service is unavailable

## Citation 📚

Morán F. 2026. ViraLong-ID: long-read viral identification and global phylogeny pipeline. GitHub repository. Available at: https://github.com/Moranfel/ViraLong-ID
