![alt text](<Screenshot 2026-03-16 at 6.57.19 PM.png>)
**Long-read viral identification and global phylogeny pipeline**

ViraLong-ID is a command-line pipeline for **fast viral identification from long-read sequencing data** with a final **global phylogeny** built from all retained target contigs across the batch.

It is designed to process **multiple samples in one run**, while reusing shared steps such as reference download, reference preparation, and BLAST database construction.

## Why use ViraLong-ID? 🚀

- 🦠 Detect viral contigs from long-read data
- 📦 Process multiple samples in a single run
- 🧬 Reuse shared references for the full batch
- 🔎 Build the BLAST database only once
- 🌿 Generate one global MAFFT alignment
- 🌳 Infer one global IQ-TREE phylogeny
- 📄 Export the final tree as both `treefile` and `PDF`
- 📁 Keep logs and outputs well organized

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
4. 🧱 Assemble reads with Flye
5. 🔎 BLAST contigs against the shared virus database
6. 🎯 Retain contigs matching the selected target taxon

### Global phylogeny steps

1. 📦 Combine retained contigs from all samples
2. 🧬 Build a global MAFFT alignment
3. 🌳 Infer a global maximum-likelihood tree with IQ-TREE
4. 📄 Render the final tree as PDF
5. 📝 Write per-sample and batch summary reports

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

## Example command 💻

```bash
python ViraLong-ID_v4.2.py \
  --taxid 1214459 \
  --reads /path/to/Selection.fastq.gz /path/to/Test.fastq.gz \
  --sample-names Selection Test \
  --outdir resultados \
  --refseq-virus-fasta /path/to/viral.1.1.genomic.fna \
  --threads 16
```

## Main parameters 🧾

- `--taxid`: target NCBI Taxonomy ID
- `--reads`: one or more input `FASTQ` / `FASTQ.GZ` files
- `--sample-names`: optional sample names in the same order as `--reads`
- `--outdir`: output directory
- `--refseq-virus-fasta`: local FASTA used to build the shared BLAST database
- `--threads`: number of threads
- `--min-q`: minimum mean read quality for `fastplong`
- `--flye-mode`: Flye mode (`normal` or `meta`)
- `--min-pident`: minimum BLAST percent identity for target contig retention
- `--min-qcov`: minimum BLAST query coverage for target contig retention
- `--min-contig-len-phylo`: minimum contig length retained for phylogeny
- `--assembly-min-q`: minimum mean read quality retained for assembly
- `--assembly-min-len`: minimum read length for assembly preselection or `auto`
- `--assembly-max-len`: maximum read length for assembly preselection or `auto`
- `--assembly-target-cov`: target effective coverage used to cap Flye input

## Output structure 📂

The batch run creates **shared outputs** at the top level and **sample-specific outputs** inside `samples/`.

```text
outdir/
├── 00_logs/
├── 01_references/
├── 05_blast_database/
├── 07_phylogeny_alignment/
├── 08_phylogeny_tree/
├── 09_report/
└── samples/
    ├── Sample1/
    └── Sample2/
```

### Shared outputs 🌍

- `01_references/` → downloaded target references and metadata
- `05_blast_database/` → shared BLAST database files
- `07_phylogeny_alignment/alignment_mafft.fasta` → global alignment
- `08_phylogeny_tree/alignment_mafft.treefile` → final Newick tree
- `08_phylogeny_tree/alignment_mafft.tree.pdf` → rendered tree in PDF format
- `09_report/multi_sample_summary.tsv` → batch summary table
- `09_report/multi_sample_summary.txt` → batch summary report

### Per-sample outputs 🧪

Each sample folder contains:

- filtered reads
- renamed reads
- assembly input subset
- Flye assembly output
- BLAST identification results
- retained target contigs
- per-sample report

## Terminal experience 🎨

ViraLong-ID includes a more visual terminal interface with:

- 🦠 a startup logo
- 📊 progress bars
- 🎯 clear step sections
- 📁 compact output summaries
- `-h` help message with styled descriptions

## Important notes ⚠️

- Shared steps are run **only once per batch**
- The final phylogeny is **global**, not per sample
- If no target contigs are retained across the batch, the phylogeny cannot be built
- External tool output is hidden from the terminal and saved in log files

## Citation 📚

Morán F. 2026. ViraLong-ID: long-read viral identification and global phylogeny pipeline. GitHub repository. Available at: https://github.com/Moranfel/ViraLong-ID
 

