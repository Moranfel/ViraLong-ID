#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ViraLong-ID v3.1
Long-read viral identification and phylogeny pipeline.

Main features
- English-only output
- Clean terminal progress bar by pipeline step
- External tool output hidden from terminal and saved in 00_logs/
- Resume support: completed steps are skipped automatically
- Short final report
- Assembly preselection step to reduce excessive ONT input:
  * mean Q filter
  * read length interval filter
  * target assembly coverage cap

Pipeline
1. Download complete viral genomes for a target TaxID from NCBI
2. Build metadata and rename target references
3. Filter long reads with fastplong
4. Shorten read headers for Flye
5. Preselect reads for assembly
6. Assemble with Flye
7. BLAST contigs against local RefSeq Virus database
8. Keep contigs matching the selected target TaxID
9. Align target contigs against target references with MAFFT
10. Build ML phylogeny with IQ-TREE
11. Write a short report
"""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
import textwrap
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Tuple

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


PIPELINE_NAME = "ViraLong-ID"
PIPELINE_VERSION = "5.0-batch"


# ---------------------------------------------------------------------
# Basic utilities
# ---------------------------------------------------------------------

def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def eprint(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def die(msg: str, code: int = 1) -> None:
    eprint(f"[{now()}] ERROR: {msg}")
    sys.exit(code)


def warn(msg: str) -> None:
    eprint(f"[{now()}] WARNING: {msg}")


def require_executable(name: str) -> None:
    if shutil.which(name) is None:
        die(f"Required executable not found: {name}")


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def open_maybe_gz(path: Path, mode: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def sanitize_field(text: str | None) -> str:
    if text is None:
        return "NA"
    text = str(text).strip()
    if not text:
        return "NA"
    for old, new in [
        ("/", "-"),
        ("\\", "-"),
        ("|", "-"),
        (":", "-"),
        (";", "-"),
        (",", "-"),
        ("(", "-"),
        (")", "-"),
        ("[", "-"),
        ("]", "-"),
        ("{", "-"),
        ("}", "-"),
        ('"', "-"),
        ("'", "-"),
        (" ", "_"),
    ]:
        text = text.replace(old, new)
    while "__" in text:
        text = text.replace("__", "_")
    while "--" in text:
        text = text.replace("--", "-")
    return text.strip("_-") or "NA"


def run_logged(cmd: List[str], log_file: Path, cwd: Path | None = None) -> None:
    mkdir(log_file.parent)
    with open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"\n[{now()}] CMD: {' '.join(cmd)}\n")
        logh.flush()
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=logh,
            stderr=logh,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"Command failed: {' '.join(cmd)}")


def fasta_lengths(fasta_path: Path) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        out[rec.id] = len(rec.seq)
    return out


def parse_flye_genome_size(bp: int) -> str:
    if bp >= 1_000_000:
        return f"{bp/1_000_000:.2f}m"
    if bp >= 1_000:
        return f"{bp/1_000:.2f}k"
    return str(bp)


def count_fastq_reads(path: Path) -> int:
    if not path.exists():
        return 0
    n_lines = 0
    try:
        with open_maybe_gz(path, "rt") as fh:
            for n_lines, _ in enumerate(fh, start=1):
                pass
    except (EOFError, OSError):
        return 0
    if n_lines == 0:
        return 0
    return n_lines // 4


def fastq_output_usable(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    return count_fastq_reads(path) > 0


def ensure_no_corrupt_fastq(path: Path, log_file: Path | None = None) -> None:
    if not path.exists() or path.stat().st_size == 0:
        return
    if fastq_output_usable(path):
        return
    warn(f"Detected an unreadable or truncated FASTQ output at {path}. It will be regenerated.")
    if log_file is not None:
        with open(log_file, "a", encoding="utf-8") as logh:
            logh.write(
                f"[{now()}] WARNING: Existing FASTQ output is unreadable or truncated and will be regenerated: {path}\n"
            )
    path.unlink(missing_ok=True)


def mean_phred(qualities: List[int]) -> float:
    if not qualities:
        return 0.0
    return sum(qualities) / len(qualities)


# ---------------------------------------------------------------------
# Terminal UI
# ---------------------------------------------------------------------

USE_COLOR = sys.stdout.isatty()


def style(text: str, color: str = "", bold: bool = False) -> str:
    if not USE_COLOR:
        return text
    colors = {
        "blue": "34",
        "cyan": "36",
        "green": "32",
        "yellow": "33",
        "red": "31",
        "magenta": "35",
        "white": "37",
    }
    codes = []
    if bold:
        codes.append("1")
    if color in colors:
        codes.append(colors[color])
    if not codes:
        return text
    return f"\033[{';'.join(codes)}m{text}\033[0m"


def print_banner(title: str, subtitle: str | None = None) -> None:
    width = 78
    print(style("=" * width, "cyan", bold=True))
    print(style(title.center(width), "cyan", bold=True))
    if subtitle:
        print(style(subtitle.center(width), "white"))
    print(style("=" * width, "cyan", bold=True))


def print_logo() -> None:
    logo = r"""
🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬
██╗   ██╗██╗██████╗  █████╗ ██╗      ██████╗ ███╗   ██╗ ██████╗       ██╗██████╗
██║   ██║██║██╔══██╗██╔══██╗██║     ██╔═══██╗████╗  ██║██╔════╝       ██║██╔══██╗
██║   ██║██║██████╔╝███████║██║     ██║   ██║██╔██╗ ██║██║  ███╗🧬🦠🧬██║██║  ██║
╚██╗ ██╔╝██║██╔══██╗██╔══██║██║     ██║   ██║██║╚██╗██║██║   ██║      ██║██║  ██║
 ╚████╔╝ ██║██║  ██║██║  ██║███████╗╚██████╔╝██║ ╚████║╚██████╔╝      ██║██████╔╝
  ╚═══╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝       ╚═╝╚═════╝
🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬 By Felix Morán 🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬
"""
    print(style(logo, "cyan", bold=True))
    print(style("  Long-read viral identification and global phylogeny for multi-sample batches", "white", bold=True))
    print()


class LogoHelpFormatter(argparse.RawTextHelpFormatter):
    pass


class LogoArgumentParser(argparse.ArgumentParser):
    def format_help(self) -> str:
        from io import StringIO

        buffer = StringIO()
        stdout = sys.stdout
        try:
            sys.stdout = buffer
            print_logo()
        finally:
            sys.stdout = stdout
        return buffer.getvalue() + super().format_help()


def print_section(title: str) -> None:
    print()
    print(style(f"[ {title} ]", "blue", bold=True))


def print_status_line(label: str, value: str, color: str = "white") -> None:
    print(f"  {style(label + ':', bold=True)} {style(value, color)}")


def draw_progress(step_index: int, total_steps: int, label: str, status: str) -> None:
    width = 30
    filled = int(width * step_index / total_steps)
    bar = "█" * filled + "·" * (width - filled)
    pct = int(100 * step_index / total_steps)
    color = {
        "RUNNING": "yellow",
        "DONE": "green",
        "SKIPPED": "blue",
        "FAILED": "red",
    }.get(status, "white")
    line = f"\r[{bar}] {pct:3d}% | Step {step_index}/{total_steps} | {label} | {status:<7}"
    print(style(line, color, bold=(status != "RUNNING")), end="", flush=True)
    if step_index == total_steps and status in {"DONE", "SKIPPED"}:
        print()


def summarize_path(path: Path) -> str:
    return str(path)


def strip_fastq_suffix(path: Path) -> str:
    name = path.name
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return path.stem


def sample_name_from_reads(path: Path) -> str:
    return sanitize_field(strip_fastq_suffix(path)) or "sample"


def mafft_direction_flag(mode: str) -> List[str]:
    if mode == "accurate":
        return ["--adjustdirectionaccurately"]
    if mode == "on":
        return ["--adjustdirection"]
    return []


# ---------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------

def make_shared_layout(base: Path) -> Dict[str, Path]:
    layout = {
        "logs": base / "00_logs",
        "refs": base / "01_references",
        "blast_db": base / "05_blast_database",
        "combined": base / "06_combined_target_contigs",
        "aln": base / "07_phylogeny_alignment",
        "identity": base / "07b_pairwise_identity",
        "tree": base / "08_phylogeny_tree",
        "report": base / "09_report",
        "samples": base / "samples",
        "tmp": base / "tmp",
    }
    for p in layout.values():
        mkdir(p)
    return layout


def make_sample_layout(shared_layout: Dict[str, Path], sample_name: str) -> Dict[str, Path]:
    base = shared_layout["samples"] / sample_name
    layout = {
        "base": base,
        "logs": base / "00_logs",
        "qc": base / "02_reads_qc",
        "renamed_reads": base / "03_reads_renamed",
        "assembly_reads": base / "03b_reads_for_assembly",
        "assembly": base / "04_assembly_flye",
        "blast": base / "05_blast_identification",
        "target": base / "06_taxon_filtered_contigs",
        "report": base / "09_report",
        "tmp": base / "tmp",
    }
    for p in layout.values():
        mkdir(p)
    return layout


# ---------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = LogoArgumentParser(
        prog="ViraLong-ID.py",
        formatter_class=LogoHelpFormatter,
        description=(
            "🧬 Long-read viral identification and phylogeny pipeline for one or more samples\n"
            "🦠 Shared references and BLAST database across the batch\n"
            "🌳 One combined alignment with robust strand correction, trimAl filtering, identity heatmap generation, and one global IQ-TREE phylogeny for all retained contigs"
        )
    )
    p.add_argument("--taxid", required=True, help="🧬 Target NCBI Taxonomy ID")
    p.add_argument("--reads", required=True, nargs="+", type=Path,
                   help="📥 One or more input FASTQ / FASTQ.GZ files")
    p.add_argument("--outdir", required=True, type=Path, help="📂 Output directory")
    p.add_argument("--refseq-virus-fasta", required=True, type=Path,
                   help="🗃️ Local RefSeq virus FASTA for BLAST database")
    p.add_argument("--sample-names", nargs="*",
                   help="🏷️ Optional sample names, same order as --reads")
    p.add_argument("--threads", type=int, default=8, help="🧵 Threads")
    p.add_argument("--min-q", type=int, default=15,
                   help="✨ Minimum mean read quality for fastplong")
    p.add_argument("--flye-mode", choices=["normal", "meta"], default="normal",
                   help="🧱 Flye mode")
    p.add_argument("--min-pident", type=float, default=70.0,
                   help="🎯 Minimum BLAST identity for target contig selection")
    p.add_argument("--min-qcov", type=float, default=40.0,
                   help="📏 Minimum BLAST query coverage for target contig selection")
    p.add_argument("--min-contig-len-phylo", type=int, default=300,
                   help="🌿 Minimum contig length retained for phylogeny")
    p.add_argument("--assembly-min-q", type=float, default=20.0,
                   help="🔬 Minimum mean Q for reads retained for assembly")
    p.add_argument("--assembly-min-len", default="auto",
                   help='📐 Minimum read length for assembly preselection, or "auto"')
    p.add_argument("--assembly-max-len", default="auto",
                   help='📏 Maximum read length for assembly preselection, or "auto"')
    p.add_argument("--assembly-target-cov", type=int, default=300,
                   help="🧮 Target effective coverage used to cap input for Flye")
    p.add_argument("--trimal-gap-threshold", type=float, default=0.8,
                   help="✂️ trimAl gap threshold used to keep well-aligned columns")
    p.add_argument("--mafft-adjust-direction", choices=["off", "on", "accurate"], default="on",
                   help="🔄 Automatic strand correction in MAFFT")
    p.add_argument("--identity-plot-min", type=float, default=None,
                   help="🎨 Minimum value for pairwise identity heatmap color scale")
    return p


# ---------------------------------------------------------------------
# Step 1 - Download NCBI data
# ---------------------------------------------------------------------

def step1_outputs(layout: Dict[str, Path], taxid: str) -> Tuple[Path, Path, Path]:
    zip_path = layout["refs"] / f"taxid_{taxid}.zip"
    dataset_dir = layout["refs"] / f"taxid_{taxid}_dataset"
    raw_fasta = dataset_dir / "ncbi_dataset" / "data" / "genomic.fna"
    jsonl = dataset_dir / "ncbi_dataset" / "data" / "data_report.jsonl"
    return zip_path, raw_fasta, jsonl


def step1_done(layout: Dict[str, Path], taxid: str) -> bool:
    _, raw_fasta, jsonl = step1_outputs(layout, taxid)
    return raw_fasta.exists() and jsonl.exists()


def step1_download_ncbi(layout: Dict[str, Path], taxid: str) -> None:
    log_file = layout["logs"] / "step01_download_ncbi.log"
    zip_path, raw_fasta, jsonl = step1_outputs(layout, taxid)
    dataset_dir = layout["refs"] / f"taxid_{taxid}_dataset"

    if dataset_dir.exists():
        shutil.rmtree(dataset_dir)

    run_logged(
        ["datasets", "download", "virus", "genome", "taxon", taxid, "--filename", str(zip_path)],
        log_file
    )
    mkdir(dataset_dir)
    run_logged(["unzip", "-o", str(zip_path), "-d", str(dataset_dir)], log_file)

    if not raw_fasta.exists() or not jsonl.exists():
        raise RuntimeError("NCBI dataset download did not produce expected files")


# ---------------------------------------------------------------------
# Step 2 - Metadata + renamed references
# ---------------------------------------------------------------------

def step2_outputs(layout: Dict[str, Path]) -> Tuple[Path, Path, Path]:
    meta_tsv = layout["refs"] / "virus_metadata.tsv"
    complete_meta = layout["refs"] / "virus_metadata_complete.tsv"
    renamed_refs = layout["refs"] / "complete_genomes_renamed.fna"
    return meta_tsv, complete_meta, renamed_refs


def step2_done(layout: Dict[str, Path]) -> bool:
    meta_tsv, complete_meta, renamed_refs = step2_outputs(layout)
    return meta_tsv.exists() and complete_meta.exists() and renamed_refs.exists()


def build_metadata_table(jsonl: Path, meta_tsv: Path, complete_meta: Path, log_file: Path) -> None:
    cmd = [
        "dataformat", "tsv", "virus-genome",
        "--inputfile", str(jsonl),
        "--fields", "accession,completeness,isolate-lineage,host-name,geo-location,length",
    ]
    with open(meta_tsv, "w", encoding="utf-8") as out, open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"\n[{now()}] CMD: {' '.join(cmd)}\n")
        proc = subprocess.run(cmd, stdout=out, stderr=logh, text=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError("dataformat failed")

    with open(meta_tsv, "r", encoding="utf-8") as inp, open(complete_meta, "w", encoding="utf-8", newline="") as out:
        reader = csv.DictReader(inp, delimiter="\t")
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            if (row.get("Completeness") or "").strip().lower() == "complete":
                writer.writerow(row)


def rename_reference_headers(raw_fasta: Path, complete_meta: Path, renamed_refs: Path) -> int:
    meta: Dict[str, Dict[str, str]] = {}
    with open(complete_meta, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            acc = row["Accession"].strip()
            meta[acc] = {
                "isolate": sanitize_field(row.get("Isolate Lineage")),
                "host": sanitize_field(row.get("Host Name")),
                "origin": sanitize_field(row.get("Geographic Location")),
            }

    records = []
    for rec in SeqIO.parse(str(raw_fasta), "fasta"):
        acc = rec.id.split()[0]
        if acc not in meta:
            continue
        rec.id = f"{acc}/{meta[acc]['isolate']}/{meta[acc]['host']}/{meta[acc]['origin']}"
        rec.name = rec.id
        rec.description = ""
        records.append(rec)

    if not records:
        raise RuntimeError("No complete target genomes were retained after metadata filtering")

    SeqIO.write(records, str(renamed_refs), "fasta")
    return len(records)


def step2_build_refs(layout: Dict[str, Path], taxid: str) -> None:
    log_file = layout["logs"] / "step02_build_refs.log"
    _, raw_fasta, jsonl = step1_outputs(layout, taxid)
    meta_tsv, complete_meta, renamed_refs = step2_outputs(layout)

    build_metadata_table(jsonl, meta_tsv, complete_meta, log_file)
    rename_reference_headers(raw_fasta, complete_meta, renamed_refs)


# ---------------------------------------------------------------------
# Step 3 - Read QC
# ---------------------------------------------------------------------

def step3_outputs(layout: Dict[str, Path], min_q: int) -> Tuple[Path, Path, Path]:
    fq = layout["qc"] / f"reads.Q{min_q}.fastq.gz"
    html = layout["qc"] / "fastplong.html"
    js = layout["qc"] / "fastplong.json"
    return fq, html, js


def step3_done(layout: Dict[str, Path], min_q: int) -> bool:
    fq, html, js = step3_outputs(layout, min_q)
    return fastq_output_usable(fq)


def step3_qc_reads(layout: Dict[str, Path], reads: Path, min_q: int, threads: int) -> None:
    log_file = layout["logs"] / "step03_qc_reads.log"
    fq, html, js = step3_outputs(layout, min_q)
    ensure_no_corrupt_fastq(fq, log_file)
    cmd = [
        "fastplong",
        "-i", str(reads),
        "-o", str(fq),
        "-m", str(min_q),
        "-w", str(threads),
        "-h", str(html),
        "-j", str(js),
    ]

    try:
        run_logged(cmd, log_file)
    except Exception:
        if fastq_output_usable(fq):
            warn(
                f"fastplong reported an error for {reads.name}, but a usable filtered FASTQ was produced. "
                "Continuing without HTML/JSON QC reports."
            )
            with open(log_file, "a", encoding="utf-8") as logh:
                logh.write(
                    f"[{now()}] WARNING: fastplong exited with an error code, "
                    "but reads.Q*.fastq.gz was produced and will be used.\n"
                )
            return
        raise

    if not fastq_output_usable(fq):
        raise RuntimeError("fastplong output FASTQ not found or empty")


# ---------------------------------------------------------------------
# Step 4 - Short read headers
# ---------------------------------------------------------------------

def step4_outputs(layout: Dict[str, Path]) -> Path:
    return layout["renamed_reads"] / "reads.short_headers.fastq.gz"


def step4_done(layout: Dict[str, Path]) -> bool:
    return step4_outputs(layout).exists()


def shorten_ont_headers(in_fastq: Path, out_fastq: Path) -> int:
    n_reads = 0
    with open_maybe_gz(in_fastq, "rt") as fin, gzip.open(out_fastq, "wt") as fout:
        block: List[str] = []
        for line in fin:
            block.append(line.rstrip("\n"))
            if len(block) == 4:
                n_reads += 1
                _, seq, _, qual = block
                fout.write(f"@read_{n_reads:08d}\n")
                fout.write(seq + "\n+\n" + qual + "\n")
                block = []
    return n_reads


def step4_rename_reads(layout: Dict[str, Path], min_q: int) -> None:
    log_file = layout["logs"] / "step04_rename_reads.log"
    in_fastq, _, _ = step3_outputs(layout, min_q)
    out_fastq = step4_outputs(layout)
    n = shorten_ont_headers(in_fastq, out_fastq)
    with open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"[{now()}] Renamed reads: {n}\n")
    if not out_fastq.exists():
        raise RuntimeError("Short-header FASTQ was not created")


# ---------------------------------------------------------------------
# Step 5 - Read preselection for assembly
# ---------------------------------------------------------------------

def step5_outputs(layout: Dict[str, Path]) -> Tuple[Path, Path]:
    subset_fastq = layout["assembly_reads"] / "reads.assembly_subset.fastq.gz"
    stats_tsv = layout["assembly_reads"] / "reads.assembly_subset.stats.tsv"
    return subset_fastq, stats_tsv


def step5_done(layout: Dict[str, Path]) -> bool:
    subset_fastq, stats_tsv = step5_outputs(layout)
    return subset_fastq.exists() and stats_tsv.exists()


def parse_len_arg(value: str, default_value: int) -> int:
    if isinstance(value, int):
        return value
    if str(value).lower() == "auto":
        return default_value
    return int(value)


def step5_select_reads_for_assembly(
    sample_layout: Dict[str, Path],
    shared_layout: Dict[str, Path],
    assembly_min_q: float,
    assembly_min_len: str,
    assembly_max_len: str,
    assembly_target_cov: int
) -> None:
    log_file = sample_layout["logs"] / "step05_preselect_assembly_reads.log"
    in_fastq = step4_outputs(sample_layout)
    renamed_refs = step2_outputs(shared_layout)[2]
    subset_fastq, stats_tsv = step5_outputs(sample_layout)

    ref_lengths = fasta_lengths(renamed_refs)
    if not ref_lengths:
        raise RuntimeError("Reference FASTA is empty")
    genome_size = max(ref_lengths.values())

    auto_min_len = max(1000, int(genome_size * 0.15))
    auto_max_len = int(genome_size * 2.5)

    min_len = parse_len_arg(assembly_min_len, auto_min_len)
    max_len = parse_len_arg(assembly_max_len, auto_max_len)

    if min_len > max_len:
        raise RuntimeError("assembly-min-len cannot be greater than assembly-max-len")

    target_bases = genome_size * assembly_target_cov
    candidates: List[Tuple[float, int, object]] = []

    total_reads = 0
    total_bases = 0
    qpass_reads = 0
    lenpass_reads = 0

    with open_maybe_gz(in_fastq, "rt") as handle:
        for rec in SeqIO.parse(handle, "fastq"):
            total_reads += 1
            read_len = len(rec.seq)
            total_bases += read_len
            read_q = mean_phred(rec.letter_annotations["phred_quality"])

            if read_q < assembly_min_q:
                continue
            qpass_reads += 1

            if read_len < min_len or read_len > max_len:
                continue
            lenpass_reads += 1

            candidates.append((read_q, read_len, rec))

    if not candidates:
        raise RuntimeError("No reads passed the assembly preselection filters")

    candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)

    selected = []
    selected_reads = 0
    selected_bases = 0

    for read_q, read_len, rec in candidates:
        selected.append(rec)
        selected_reads += 1
        selected_bases += read_len
        if selected_bases >= target_bases:
            break

    if not selected:
        raise RuntimeError("No reads were selected for assembly")

    with gzip.open(subset_fastq, "wt") as out:
        SeqIO.write(selected, out, "fastq")

    estimated_cov = selected_bases / genome_size if genome_size > 0 else 0.0

    with open(stats_tsv, "w", encoding="utf-8", newline="") as out:
        fields = [
            "genome_size_bp",
            "total_reads_input",
            "total_bases_input",
            "min_q",
            "min_len",
            "max_len",
            "target_bases",
            "reads_q_pass",
            "reads_len_pass",
            "selected_reads",
            "selected_bases",
            "estimated_selected_coverage",
        ]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "genome_size_bp": genome_size,
            "total_reads_input": total_reads,
            "total_bases_input": total_bases,
            "min_q": assembly_min_q,
            "min_len": min_len,
            "max_len": max_len,
            "target_bases": target_bases,
            "reads_q_pass": qpass_reads,
            "reads_len_pass": lenpass_reads,
            "selected_reads": selected_reads,
            "selected_bases": selected_bases,
            "estimated_selected_coverage": f"{estimated_cov:.2f}",
        })

    with open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"[{now()}] Genome size: {genome_size} bp\n")
        logh.write(f"[{now()}] Input reads: {total_reads}\n")
        logh.write(f"[{now()}] Input bases: {total_bases}\n")
        logh.write(f"[{now()}] Assembly min Q: {assembly_min_q}\n")
        logh.write(f"[{now()}] Assembly min len: {min_len}\n")
        logh.write(f"[{now()}] Assembly max len: {max_len}\n")
        logh.write(f"[{now()}] Target bases: {target_bases}\n")
        logh.write(f"[{now()}] Reads passing Q: {qpass_reads}\n")
        logh.write(f"[{now()}] Reads passing length: {lenpass_reads}\n")
        logh.write(f"[{now()}] Selected reads: {selected_reads}\n")
        logh.write(f"[{now()}] Selected bases: {selected_bases}\n")
        logh.write(f"[{now()}] Estimated selected coverage: {estimated_cov:.2f}x\n")


# ---------------------------------------------------------------------
# Step 6 - Flye assembly
# ---------------------------------------------------------------------

def step6_outputs(sample_layout: Dict[str, Path]) -> Path:
    return sample_layout["assembly"] / "assembly.fasta"


def step6_done(sample_layout: Dict[str, Path]) -> bool:
    return step6_outputs(sample_layout).exists()


def step6_assemble(sample_layout: Dict[str, Path], shared_layout: Dict[str, Path], threads: int, flye_mode: str) -> None:
    log_file = sample_layout["logs"] / "step06_flye.log"
    renamed_refs = step2_outputs(shared_layout)[2]
    subset_fastq = step5_outputs(sample_layout)[0]
    assembly_fasta = step6_outputs(sample_layout)

    ref_lengths = fasta_lengths(renamed_refs)
    if not ref_lengths:
        raise RuntimeError("Reference FASTA is empty")
    max_ref = max(ref_lengths.values())

    cmd = [
        "flye",
        "--nano-raw", str(subset_fastq),
        "--out-dir", str(sample_layout["assembly"]),
        "--threads", str(threads),
        "--genome-size", parse_flye_genome_size(max_ref),
    ]
    if flye_mode == "meta":
        cmd.append("--meta")

    run_logged(cmd, log_file)
    if not assembly_fasta.exists():
        raise RuntimeError("Flye assembly.fasta not found")


# ---------------------------------------------------------------------
# Step 7 - BLAST
# ---------------------------------------------------------------------

def step7_outputs(sample_layout: Dict[str, Path]) -> Tuple[Path, Path]:
    blast_tsv = sample_layout["blast"] / "assembly_vs_refseq_virus.tsv"
    top_hits = sample_layout["blast"] / "top_hits_per_contig.tsv"
    return blast_tsv, top_hits


def step7_done(sample_layout: Dict[str, Path]) -> bool:
    blast_tsv, top_hits = step7_outputs(sample_layout)
    return blast_tsv.exists() and top_hits.exists()


def step7_db_outputs(shared_layout: Dict[str, Path]) -> Tuple[Path, List[Path]]:
    db_prefix = shared_layout["blast_db"] / "refseq_virus"
    db_files = [
        db_prefix.with_suffix(".nhr"),
        db_prefix.with_suffix(".nin"),
        db_prefix.with_suffix(".nsq"),
    ]
    return db_prefix, db_files


def step7_db_done(shared_layout: Dict[str, Path]) -> bool:
    _, db_files = step7_db_outputs(shared_layout)
    return all(path.exists() for path in db_files)


def best_hit_per_contig(blast_tsv: Path, top_hits_tsv: Path) -> Dict[str, Dict[str, str]]:
    best: Dict[str, Dict[str, str]] = {}
    with open(blast_tsv, "r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            rec = {
                "qseqid": row[0],
                "sacc": row[1],
                "stitle": row[2],
                "pident": float(row[3]),
                "length": int(row[4]),
                "qlen": int(row[5]),
                "qstart": int(row[6]),
                "qend": int(row[7]),
                "sstart": int(row[8]),
                "send": int(row[9]),
                "evalue": row[10],
                "bitscore": float(row[11]),
                "qcovs": float(row[12]),
            }
            qseqid = rec["qseqid"]
            if qseqid not in best or (rec["bitscore"], rec["pident"], rec["qcovs"]) > (
                best[qseqid]["bitscore"], best[qseqid]["pident"], best[qseqid]["qcovs"]
            ):
                best[qseqid] = rec

    with open(top_hits_tsv, "w", encoding="utf-8", newline="") as out:
        fields = ["qseqid", "sacc", "stitle", "pident", "length", "qlen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(best):
            writer.writerow(best[key])

    return best


def step7_build_blast_db(shared_layout: Dict[str, Path], refseq_virus_fasta: Path) -> None:
    log_file = shared_layout["logs"] / "step03_build_blast_db.log"
    db_prefix, _ = step7_db_outputs(shared_layout)

    run_logged(
        ["makeblastdb", "-in", str(refseq_virus_fasta), "-dbtype", "nucl", "-out", str(db_prefix)],
        log_file
    )


def step7_blast(sample_layout: Dict[str, Path], shared_layout: Dict[str, Path], threads: int) -> None:
    log_file = sample_layout["logs"] / "step07_blast.log"
    assembly_fasta = step6_outputs(sample_layout)
    blast_tsv, top_hits = step7_outputs(sample_layout)
    db_prefix, _ = step7_db_outputs(shared_layout)

    outfmt = "6 qseqid sacc stitle pident length qlen qstart qend sstart send evalue bitscore qcovs"
    run_logged(
        [
            "blastn",
            "-task", "megablast",
            "-query", str(assembly_fasta),
            "-db", str(db_prefix),
            "-out", str(blast_tsv),
            "-num_threads", str(threads),
            "-max_target_seqs", "20",
            "-evalue", "1e-10",
            "-outfmt", outfmt,
        ],
        log_file
    )
    best_hit_per_contig(blast_tsv, top_hits)

    if not blast_tsv.exists() or not top_hits.exists():
        raise RuntimeError("BLAST outputs not created")


# ---------------------------------------------------------------------
# Step 8 - Target contig selection
# ---------------------------------------------------------------------

def step8_outputs(sample_layout: Dict[str, Path]) -> Tuple[Path, Path]:
    fasta = sample_layout["target"] / "contigs_target_taxon.fasta"
    tsv = sample_layout["target"] / "contigs_target_taxon.tsv"
    return fasta, tsv


def step8_done(sample_layout: Dict[str, Path]) -> bool:
    fasta, tsv = step8_outputs(sample_layout)
    return fasta.exists() and tsv.exists()


def load_top_hits(top_hits_tsv: Path) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    with open(top_hits_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            row["pident"] = float(row["pident"])
            row["qcovs"] = float(row["qcovs"])
            row["bitscore"] = float(row["bitscore"])
            out[row["qseqid"]] = row
    return out


def get_target_accessions(renamed_refs: Path) -> set:
    accs = set()
    for rec in SeqIO.parse(str(renamed_refs), "fasta"):
        accs.add(rec.id.split("/")[0])
    return accs


def step8_select_target_contigs(sample_layout: Dict[str, Path], shared_layout: Dict[str, Path], min_pident: float, min_qcov: float,
                                min_contig_len_phylo: int, sample_name: str) -> None:
    log_file = sample_layout["logs"] / "step08_target_contigs.log"
    renamed_refs = step2_outputs(shared_layout)[2]
    top_hits_tsv = step7_outputs(sample_layout)[1]
    assembly_fasta = step6_outputs(sample_layout)
    out_fasta, out_tsv = step8_outputs(sample_layout)

    target_accessions = get_target_accessions(renamed_refs)
    top_hits = load_top_hits(top_hits_tsv)

    keep_hits = []
    keep_ids = set()

    for qseqid, hit in top_hits.items():
        sacc = str(hit["sacc"]).split()[0]
        if sacc in target_accessions and float(hit["pident"]) >= min_pident and float(hit["qcovs"]) >= min_qcov:
            keep_ids.add(qseqid)
            keep_hits.append(hit)

    out_records = []
    for rec in SeqIO.parse(str(assembly_fasta), "fasta"):
        if rec.id in keep_ids and len(rec.seq) >= min_contig_len_phylo:
            rec.id = f"{sample_name}__{rec.id}"
            rec.name = rec.id
            rec.description = ""
            out_records.append(rec)

    SeqIO.write(out_records, str(out_fasta), "fasta")

    with open(out_tsv, "w", encoding="utf-8", newline="") as out:
        fields = ["qseqid", "sacc", "stitle", "pident", "length", "qlen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in keep_hits:
            writer.writerow(row)

    with open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"[{now()}] Retained target contigs: {len(out_records)}\n")


# ---------------------------------------------------------------------
# Step 9 - Combined contigs and MAFFT
# ---------------------------------------------------------------------

def step9_outputs(shared_layout: Dict[str, Path]) -> Tuple[Path, Path, Path, Path, Path]:
    combined = shared_layout["combined"] / "all_target_contigs.fasta"
    guide = shared_layout["aln"] / "direction_guide_plus_contigs.fasta"
    oriented = shared_layout["aln"] / "all_target_contigs.oriented.fasta"
    dataset = shared_layout["aln"] / "sequences_for_phylogeny.fasta"
    aln = shared_layout["aln"] / "alignment_mafft.fasta"
    return combined, guide, oriented, dataset, aln


def step9_done(shared_layout: Dict[str, Path]) -> bool:
    combined, _, oriented, dataset, aln = step9_outputs(shared_layout)
    return combined.exists() and oriented.exists() and dataset.exists() and aln.exists()


def concatenate_fastas(files: List[Path], out_fasta: Path) -> int:
    records = []
    for f in files:
        if not f.exists() or f.stat().st_size == 0:
            continue
        for rec in SeqIO.parse(str(f), "fasta"):
            records.append(rec)
    SeqIO.write(records, str(out_fasta), "fasta")
    return len(records)


def step9_collect_and_align(shared_layout: Dict[str, Path], sample_layouts: List[Dict[str, Path]], threads: int,
                            adjust_direction: str) -> None:
    log_file = shared_layout["logs"] / "step09_mafft.log"
    renamed_refs = step2_outputs(shared_layout)[2]
    combined_contigs, direction_guide, oriented_contigs, dataset, aln = step9_outputs(shared_layout)
    ref_aln = shared_layout["aln"] / "reference_alignment.fasta"

    target_fastas = [step8_outputs(sample_layout)[0] for sample_layout in sample_layouts]
    n_contigs = concatenate_fastas(target_fastas, combined_contigs)
    if n_contigs == 0:
        raise RuntimeError("No target contigs were retained across the batch")

    cmd1 = [
        "mafft",
        "--thread", str(threads),
        *mafft_direction_flag(adjust_direction),
        "--auto",
        str(renamed_refs)
    ]
    with open(ref_aln, "w", encoding="utf-8") as out, open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"\n[{now()}] CMD: {' '.join(cmd1)}\n")
        proc = subprocess.run(cmd1, stdout=out, stderr=logh, text=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError("MAFFT reference alignment failed")

    if adjust_direction == "off":
        shutil.copyfile(str(combined_contigs), str(oriented_contigs))
    else:
        first_ref = None
        for rec in SeqIO.parse(str(renamed_refs), "fasta"):
            first_ref = rec
            break
        if first_ref is None:
            raise RuntimeError("Reference FASTA is empty")

        with open(direction_guide, "w", encoding="utf-8") as out:
            SeqIO.write([first_ref], out, "fasta")
            for rec in SeqIO.parse(str(combined_contigs), "fasta"):
                SeqIO.write([rec], out, "fasta")

        cmd_orient = [
            "mafft",
            "--thread", str(threads),
            *mafft_direction_flag(adjust_direction),
            "--auto",
            str(direction_guide)
        ]
        oriented_full = shared_layout["aln"] / "direction_guide_plus_contigs.aligned.fasta"
        with open(oriented_full, "w", encoding="utf-8") as out, open(log_file, "a", encoding="utf-8") as logh:
            logh.write(f"\n[{now()}] CMD: {' '.join(cmd_orient)}\n")
            proc = subprocess.run(cmd_orient, stdout=out, stderr=logh, text=True, check=False)
            if proc.returncode != 0:
                raise RuntimeError("MAFFT strand correction failed")

        oriented_records = []
        for rec in SeqIO.parse(str(oriented_full), "fasta"):
            if rec.id == first_ref.id or rec.id.startswith("_R_" + first_ref.id):
                continue
            seq = str(rec.seq).replace("-", "")
            rec.seq = rec.seq.__class__(seq)
            rec.id = rec.id.removeprefix("_R_")
            rec.name = rec.id
            rec.description = ""
            oriented_records.append(rec)

        if not oriented_records:
            raise RuntimeError("No oriented contigs were produced for MAFFT")
        SeqIO.write(oriented_records, str(oriented_contigs), "fasta")

    concatenate_fastas([renamed_refs, oriented_contigs], dataset)

    cmd2 = [
        "mafft",
        "--thread", str(threads),
        "--reorder",
        "--addfragments", str(oriented_contigs),
        str(ref_aln)
    ]
    with open(aln, "w", encoding="utf-8") as out, open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"\n[{now()}] CMD: {' '.join(cmd2)}\n")
        proc = subprocess.run(cmd2, stdout=out, stderr=logh, text=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError("MAFFT addfragments step failed")

    if not aln.exists() or aln.stat().st_size == 0:
        raise RuntimeError("MAFFT produced an empty final alignment file")


# ---------------------------------------------------------------------
# Step 10 - trimAl
# ---------------------------------------------------------------------

def step10_outputs(shared_layout: Dict[str, Path]) -> Tuple[Path, Path]:
    trimmed = shared_layout["aln"] / "alignment_mafft.trimmed.fasta"
    html = shared_layout["aln"] / "alignment_mafft.trimmed.html"
    return trimmed, html


def step10_done(shared_layout: Dict[str, Path]) -> bool:
    trimmed, _ = step10_outputs(shared_layout)
    return trimmed.exists() and trimmed.stat().st_size > 0


def step10_trimal(shared_layout: Dict[str, Path], gap_threshold: float) -> None:
    log_file = shared_layout["logs"] / "step10_trimal.log"
    aln = step9_outputs(shared_layout)[4]
    trimmed, html = step10_outputs(shared_layout)

    run_logged(
        [
            "trimal",
            "-in", str(aln),
            "-out", str(trimmed),
            "-htmlout", str(html),
            "-gt", str(gap_threshold),
        ],
        log_file
    )

    if not trimmed.exists() or trimmed.stat().st_size == 0:
        raise RuntimeError("trimAl did not produce a trimmed alignment")


# ---------------------------------------------------------------------
# Step 11 - Pairwise identity heatmap
# ---------------------------------------------------------------------

def step11_identity_outputs(shared_layout: Dict[str, Path]) -> Tuple[Path, Path, Path]:
    matrix_tsv = shared_layout["identity"] / "pairwise_identity.tsv"
    pdf = shared_layout["identity"] / "pairwise_identity_heatmap.pdf"
    png = shared_layout["identity"] / "pairwise_identity_heatmap.png"
    return matrix_tsv, pdf, png


def step11_identity_done(shared_layout: Dict[str, Path]) -> bool:
    matrix_tsv, pdf, png = step11_identity_outputs(shared_layout)
    return matrix_tsv.exists() and pdf.exists() and png.exists()


def load_alignment_records(alignment_fasta: Path):
    return list(SeqIO.parse(str(alignment_fasta), "fasta"))


def pairwise_identity(seq_a: str, seq_b: str) -> float:
    matches = 0
    compared = 0
    for a, b in zip(seq_a.upper(), seq_b.upper()):
        if a in "-?" or b in "-?":
            continue
        compared += 1
        if a == b:
            matches += 1
    if compared == 0:
        return 0.0
    return 100.0 * matches / compared


def compute_identity_matrix(records) -> Tuple[List[str], List[List[float]]]:
    labels = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]
    matrix: List[List[float]] = []
    for seq_a in seqs:
        row = []
        for seq_b in seqs:
            row.append(pairwise_identity(seq_a, seq_b))
        matrix.append(row)
    return labels, matrix


def reorder_matrix(labels: List[str], matrix: List[List[float]], ordered_labels: List[str]) -> Tuple[List[str], List[List[float]]]:
    index = {label: i for i, label in enumerate(labels)}
    keep = [label for label in ordered_labels if label in index]
    if not keep:
        return labels, matrix
    reordered = []
    for row_label in keep:
        i = index[row_label]
        reordered.append([matrix[i][index[col_label]] for col_label in keep])
    return keep, reordered


def render_identity_heatmap(labels: List[str], matrix: List[List[float]], pdf_path: Path, png_path: Path,
                            title: str, vmin: float | None) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    off_diag = [matrix[i][j] for i in range(len(matrix)) for j in range(len(matrix)) if i != j]
    auto_vmin = min(off_diag) if off_diag else 95.0
    if vmin is None:
        vmin = max(0.0, round((auto_vmin - 0.2) * 2) / 2)

    n = max(1, len(labels))
    fig_w = max(12, min(28, 4 + n * 0.42))
    fig_h = max(10, min(28, 4 + n * 0.34))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(matrix, cmap="viridis", vmin=vmin, vmax=100.0, interpolation="nearest", aspect="auto")
    ax.set_title(title, fontsize=14, pad=12)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xticks([x - 0.5 for x in range(1, n)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n)], minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Genome identity (%)")
    fig.tight_layout()
    fig.savefig(str(pdf_path), format="pdf", bbox_inches="tight")
    fig.savefig(str(png_path), format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def step11_identity_plot(shared_layout: Dict[str, Path], plot_min: float | None) -> None:
    from Bio import Phylo

    log_file = shared_layout["logs"] / "step11_identity_heatmap.log"
    trimmed_alignment = step10_outputs(shared_layout)[0]
    matrix_tsv, pdf, png = step11_identity_outputs(shared_layout)
    treefile = step12_tree_outputs(shared_layout)[0]
    records = load_alignment_records(trimmed_alignment)
    if not records:
        raise RuntimeError("Trimmed alignment is empty")
    labels, matrix = compute_identity_matrix(records)
    if treefile.exists():
        tree = Phylo.read(str(treefile), "newick")
        tree_order = [term.name for term in tree.get_terminals()]
        labels, matrix = reorder_matrix(labels, matrix, tree_order)

    with open(matrix_tsv, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sequence"] + labels)
        for label, row in zip(labels, matrix):
            writer.writerow([label] + [f"{value:.4f}" for value in row])

    render_identity_heatmap(labels, matrix, pdf, png, "CYVCV genome identity heatmap", plot_min)

    with open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"[{now()}] Pairwise identity matrix written to: {matrix_tsv}\n")
        logh.write(f"[{now()}] Heatmap PDF written to: {pdf}\n")
        logh.write(f"[{now()}] Heatmap PNG written to: {png}\n")


# ---------------------------------------------------------------------
# Step 12 - IQ-TREE + PDF
# ---------------------------------------------------------------------

def step12_tree_outputs(shared_layout: Dict[str, Path]) -> Tuple[Path, Path, Path]:
    treefile = shared_layout["tree"] / "alignment_mafft.trimmed.treefile"
    iqtree = shared_layout["tree"] / "alignment_mafft.trimmed.iqtree"
    pdf = shared_layout["tree"] / "alignment_mafft.trimmed.tree.pdf"
    return treefile, iqtree, pdf


def step12_tree_done(shared_layout: Dict[str, Path]) -> bool:
    treefile, iqtree, pdf = step12_tree_outputs(shared_layout)
    return treefile.exists() and iqtree.exists() and pdf.exists()


def is_assembled_tip(label: str) -> bool:
    return "__" in label


def render_tree_pdf(treefile: Path, pdf_path: Path) -> None:
    from Bio import Phylo
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    tree = Phylo.read(str(treefile), "newick")
    terminals = max(1, len(tree.get_terminals()))
    fig_width = max(12, min(36, terminals * 0.35))
    fig_height = max(8, min(42, terminals * 0.45))
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False)

    for text in ax.texts:
        label = text.get_text().strip()
        if not label:
            continue
        if is_assembled_tip(label):
            text.set_color("#c0392b")
            text.set_fontweight("bold")
        else:
            text.set_color("#1f1f1f")

    legend_handles = [
        Line2D([0], [0], color="#c0392b", lw=0, marker="o", markersize=7, label="Assembled isolates"),
        Line2D([0], [0], color="#1f1f1f", lw=0, marker="o", markersize=7, label="Reference genomes"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(str(pdf_path), format="pdf", bbox_inches="tight")
    plt.close(fig)


def step12_iqtree(shared_layout: Dict[str, Path], threads: int) -> None:
    log_file = shared_layout["logs"] / "step12_iqtree.log"
    aln = step10_outputs(shared_layout)[0]
    prefix = shared_layout["tree"] / "alignment_mafft.trimmed"

    run_logged(
        [
            "iqtree",
            "-s", str(aln),
            "-m", "MFP",
            "-bb", "1000",
            "-nt", "AUTO",
            "--prefix", str(prefix),
        ],
        log_file
    )

    treefile, iqtree_txt, pdf = step12_tree_outputs(shared_layout)
    if not treefile.exists() or not iqtree_txt.exists():
        raise RuntimeError("IQ-TREE output not found")
    render_tree_pdf(treefile, pdf)


# ---------------------------------------------------------------------
# Step 13 - Reports
# ---------------------------------------------------------------------

def step11_outputs(sample_layout: Dict[str, Path]) -> Tuple[Path, Path]:
    return sample_layout["report"] / "summary.tsv", sample_layout["report"] / "summary.txt"


def step11_done(sample_layout: Dict[str, Path]) -> bool:
    a, b = step11_outputs(sample_layout)
    return a.exists() and b.exists()


def best_hit_from_top_hits(top_hits_tsv: Path) -> str:
    if not top_hits_tsv.exists() or top_hits_tsv.stat().st_size == 0:
        return "NA"
    with open(top_hits_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        return "NA"
    rows.sort(key=lambda r: (float(r["bitscore"]), float(r["pident"]), float(r["qcovs"])), reverse=True)
    r = rows[0]
    return f'{r["stitle"]} | pident={float(r["pident"]):.2f}% | qcov={float(r["qcovs"]):.2f}%'


def sample_summary_row(sample_layout: Dict[str, Path], shared_layout: Dict[str, Path], taxid: str,
                       reads: Path, min_q: int, sample_name: str) -> Dict[str, str]:
    renamed_refs = step2_outputs(shared_layout)[2]
    assembly_fasta = step6_outputs(sample_layout)
    target_contigs = step8_outputs(sample_layout)[0]
    preselect_stats = step5_outputs(sample_layout)[1]
    trimmed_alignment, _ = step10_outputs(shared_layout)
    identity_tsv, identity_pdf, _ = step11_identity_outputs(shared_layout)
    treefile, _, tree_pdf = step12_tree_outputs(shared_layout)

    reads_in = count_fastq_reads(reads)
    reads_qc = count_fastq_reads(step3_outputs(sample_layout, min_q)[0])
    ref_lengths = fasta_lengths(renamed_refs)
    asm_lengths = fasta_lengths(assembly_fasta)
    tgt_lengths = fasta_lengths(target_contigs) if target_contigs.exists() else {}

    selected_reads = "NA"
    selected_cov = "NA"
    if preselect_stats.exists():
        with open(preselect_stats, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        if rows:
            selected_reads = rows[0]["selected_reads"]
            selected_cov = rows[0]["estimated_selected_coverage"]

    return {
        "sample": sample_name,
        "taxid": taxid,
        "reads_input": str(reads_in),
        "reads_after_qc": str(reads_qc),
        "reads_for_assembly": selected_reads,
        "assembly_selected_coverage": selected_cov,
        "n_complete_refs": str(len(ref_lengths)),
        "n_assembly_contigs": str(len(asm_lengths)),
        "n_target_contigs": str(len(tgt_lengths)),
        "best_hit": best_hit_from_top_hits(step7_outputs(sample_layout)[1]),
        "trimmed_alignment": str(trimmed_alignment),
        "identity_matrix": str(identity_tsv),
        "identity_heatmap": str(identity_pdf),
        "treefile": str(treefile),
        "tree_pdf": str(tree_pdf),
    }


def step11_report(sample_layout: Dict[str, Path], shared_layout: Dict[str, Path], taxid: str,
                  reads: Path, min_q: int, sample_name: str) -> Dict[str, str]:
    summary_tsv, summary_txt = step11_outputs(sample_layout)
    row = sample_summary_row(sample_layout, shared_layout, taxid, reads, min_q, sample_name)

    fields = [
        "sample", "taxid", "reads_input", "reads_after_qc", "reads_for_assembly",
        "assembly_selected_coverage", "n_complete_refs", "n_assembly_contigs",
        "n_target_contigs", "best_hit", "trimmed_alignment", "identity_matrix",
        "identity_heatmap", "treefile", "tree_pdf"
    ]
    with open(summary_tsv, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    txt = textwrap.dedent(f"""\
    {PIPELINE_NAME} {PIPELINE_VERSION}

    Sample: {row['sample']}
    TaxID: {row['taxid']}
    Input reads: {row['reads_input']}
    Reads after QC: {row['reads_after_qc']}
    Reads used for assembly: {row['reads_for_assembly']}
    Selected assembly coverage: {row['assembly_selected_coverage']}x
    Complete references: {row['n_complete_refs']}
    Assembly contigs: {row['n_assembly_contigs']}
    Target contigs: {row['n_target_contigs']}
    Best BLAST hit: {row['best_hit']}

    Main outputs:
    - Assembly input subset: 03b_reads_for_assembly/reads.assembly_subset.fastq.gz
    - Target contigs: 06_taxon_filtered_contigs/contigs_target_taxon.fasta
    - Global alignment: {step9_outputs(shared_layout)[4]}
    - Trimmed alignment: {row['trimmed_alignment']}
    - Pairwise identity matrix: {row['identity_matrix']}
    - Pairwise identity heatmap: {row['identity_heatmap']}
    - Global tree: {row['treefile']}
    - Global tree PDF: {row['tree_pdf']}
    """)
    with open(summary_txt, "w", encoding="utf-8") as out:
        out.write(txt)

    return row


# ---------------------------------------------------------------------
# Step runner
# ---------------------------------------------------------------------

class Step:
    def __init__(self, label, done_check, action):
        self.label = label
        self.done_check = done_check
        self.action = action


def validate_shared_inputs(args) -> None:
    if SeqIO is None:
        die("Biopython is not installed. Please install it in your environment before running the pipeline.")
    for exe in ["datasets", "dataformat", "unzip", "fastplong", "flye", "makeblastdb", "blastn", "mafft", "trimal", "iqtree"]:
        require_executable(exe)

    if not args.refseq_virus_fasta.exists():
        die(f"RefSeq virus FASTA not found: {args.refseq_virus_fasta}")


def build_sample_jobs(args) -> List[SimpleNamespace]:
    reads_list = args.reads
    if args.sample_names and len(args.sample_names) != len(reads_list):
        die("If provided, --sample-names must contain the same number of values as --reads")

    seen = set()
    jobs = []
    for idx, reads_path in enumerate(reads_list):
        if not reads_path.exists():
            die(f"Input reads file not found: {reads_path}")
        sample_name = args.sample_names[idx] if args.sample_names else sample_name_from_reads(reads_path)
        sample_name = sanitize_field(sample_name)
        if sample_name in seen:
            die(f"Duplicated sample name detected: {sample_name}")
        seen.add(sample_name)
        jobs.append(SimpleNamespace(sample_name=sample_name, reads=reads_path))

    return jobs


def read_summary_row(summary_tsv: Path) -> Dict[str, str]:
    with open(summary_tsv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise RuntimeError(f"Summary file is empty: {summary_tsv}")
    return rows[0]


def write_batch_summary(outdir: Path, rows: List[Dict[str, str]]) -> Tuple[Path, Path]:
    report_dir = outdir / "09_report"
    mkdir(report_dir)
    summary_tsv = report_dir / "multi_sample_summary.tsv"
    summary_txt = report_dir / "multi_sample_summary.txt"

    fields = [
        "sample", "taxid", "reads_input", "reads_after_qc", "reads_for_assembly",
        "assembly_selected_coverage", "n_complete_refs", "n_assembly_contigs",
        "n_target_contigs", "best_hit", "trimmed_alignment", "identity_matrix",
        "identity_heatmap", "treefile", "tree_pdf"
    ]
    with open(summary_tsv, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    lines = [
        f"{PIPELINE_NAME} {PIPELINE_VERSION}",
        "",
        f"Processed samples: {len(rows)}",
        "",
    ]
    for row in rows:
        lines.extend([
            f"- {row['sample']}: target_contigs={row['n_target_contigs']}, assembly_contigs={row['n_assembly_contigs']}",
            f"  Best hit: {row['best_hit']}",
            f"  Trimmed alignment: {row['trimmed_alignment']}",
            f"  Identity matrix: {row['identity_matrix']}",
            f"  Identity heatmap: {row['identity_heatmap']}",
            f"  Tree: {row['treefile']}",
            f"  Tree PDF: {row['tree_pdf']}",
            "",
        ])
    with open(summary_txt, "w", encoding="utf-8") as out:
        out.write("\n".join(lines).rstrip() + "\n")

    return summary_tsv, summary_txt


def run_global_setup(shared_layout: Dict[str, Path], args) -> None:
    steps = [
        Step("Download NCBI target genomes",
             lambda: step1_done(shared_layout, args.taxid),
             lambda: step1_download_ncbi(shared_layout, args.taxid)),
        Step("Build metadata and renamed references",
             lambda: step2_done(shared_layout),
             lambda: step2_build_refs(shared_layout, args.taxid)),
        Step("Build shared BLAST database",
             lambda: step7_db_done(shared_layout),
             lambda: step7_build_blast_db(shared_layout, args.refseq_virus_fasta)),
    ]

    print_section("Shared setup")
    for idx, step in enumerate(steps, start=1):
        try:
            if step.done_check():
                draw_progress(idx, len(steps), step.label, "SKIPPED")
                continue
            draw_progress(idx - 1 if idx > 1 else 0, len(steps), step.label, "RUNNING")
            step.action()
            draw_progress(idx, len(steps), step.label, "DONE")
        except Exception as exc:
            print()
            die(f"Step failed: {step.label}\nReason: {exc}")


def run_single_sample_pipeline(args, shared_layout: Dict[str, Path], sample_name: str,
                               reads: Path, sample_index: int, total_samples: int) -> Dict[str, str]:
    layout = make_sample_layout(shared_layout, sample_name)

    steps = [
        Step("Filter long reads with fastplong",
             lambda: step3_done(layout, args.min_q),
             lambda: step3_qc_reads(layout, reads, args.min_q, args.threads)),
        Step("Shorten read headers",
             lambda: step4_done(layout),
             lambda: step4_rename_reads(layout, args.min_q)),
        Step("Preselect reads for assembly",
             lambda: step5_done(layout),
             lambda: step5_select_reads_for_assembly(
                 layout,
                 shared_layout,
                 args.assembly_min_q,
                 args.assembly_min_len,
                 args.assembly_max_len,
                 args.assembly_target_cov
             )),
        Step("Assemble with Flye",
             lambda: step6_done(layout),
             lambda: step6_assemble(layout, shared_layout, args.threads, args.flye_mode)),
        Step("BLAST contigs against RefSeq Virus",
             lambda: step7_done(layout),
             lambda: step7_blast(layout, shared_layout, args.threads)),
        Step("Select target contigs",
             lambda: step8_done(layout),
             lambda: step8_select_target_contigs(
                 layout, shared_layout, args.min_pident, args.min_qcov, args.min_contig_len_phylo, sample_name
             )),
    ]

    total = len(steps)
    print_section(f"Sample {sample_index}/{total_samples}")
    print_status_line("Sample", sample_name, "magenta")
    print_status_line("Reads", summarize_path(reads))
    print_status_line("Output", summarize_path(layout["base"]))
    print_status_line("Resume mode", "enabled", "green")
    print_status_line("Tool logs", summarize_path(layout["logs"]))
    print()

    for idx, step in enumerate(steps, start=1):
        try:
            if step.done_check():
                draw_progress(idx, total, step.label, "SKIPPED")
                continue
            draw_progress(idx - 1 if idx > 1 else 0, total, step.label, "RUNNING")
            step.action()
            draw_progress(idx, total, step.label, "DONE")
        except Exception as exc:
            print()
            die(f"Step failed: {step.label}\nReason: {exc}")

    print()
    print_status_line("Sample finished", sample_name, "green")
    print_status_line("Target contigs", summarize_path(step8_outputs(layout)[0]), "green")
    print_status_line("Logs", summarize_path(layout["logs"]))
    return {"sample": sample_name}


def run_global_phylogeny(shared_layout: Dict[str, Path], sample_layouts: List[Dict[str, Path]], args) -> None:
    steps = [
        Step("Combine target contigs and build MAFFT alignment",
             lambda: step9_done(shared_layout),
             lambda: step9_collect_and_align(
                 shared_layout, sample_layouts, args.threads, args.mafft_adjust_direction
             )),
        Step("Trim alignment with trimAl",
             lambda: step10_done(shared_layout),
             lambda: step10_trimal(shared_layout, args.trimal_gap_threshold)),
        Step("Infer batch ML tree and render PDF",
             lambda: step12_tree_done(shared_layout),
             lambda: step12_iqtree(shared_layout, args.threads)),
        Step("Build pairwise identity heatmap",
             lambda: step11_identity_done(shared_layout),
             lambda: step11_identity_plot(shared_layout, args.identity_plot_min)),
    ]

    print_section("Global phylogeny")
    for idx, step in enumerate(steps, start=1):
        try:
            if step.done_check():
                draw_progress(idx, len(steps), step.label, "SKIPPED")
                continue
            draw_progress(idx - 1 if idx > 1 else 0, len(steps), step.label, "RUNNING")
            step.action()
            draw_progress(idx, len(steps), step.label, "DONE")
        except Exception as exc:
            print()
            die(f"Step failed: {step.label}\nReason: {exc}")


def run_pipeline(args) -> None:
    validate_shared_inputs(args)
    sample_jobs = build_sample_jobs(args)
    shared_layout = make_shared_layout(args.outdir)

    print_logo()
    print_banner(f"{PIPELINE_NAME} {PIPELINE_VERSION}", "Batch viral identification and global phylogeny pipeline")
    print_status_line("TaxID", args.taxid, "magenta")
    print_status_line("Samples detected", str(len(sample_jobs)), "green")
    print_status_line("Threads", str(args.threads))
    print_status_line("RefSeq virus FASTA", summarize_path(args.refseq_virus_fasta))
    print_status_line("MAFFT strand correction", args.mafft_adjust_direction, "green")
    print_status_line("External tool output", "hidden in 00_logs/", "blue")

    run_global_setup(shared_layout, args)

    sample_layouts = []
    for idx, job in enumerate(sample_jobs, start=1):
        run_single_sample_pipeline(args, shared_layout, job.sample_name, job.reads, idx, len(sample_jobs))
        sample_layouts.append(make_sample_layout(shared_layout, job.sample_name))

    run_global_phylogeny(shared_layout, sample_layouts, args)

    batch_rows = []
    for job, sample_layout in zip(sample_jobs, sample_layouts):
        batch_rows.append(step11_report(sample_layout, shared_layout, args.taxid, job.reads, args.min_q, job.sample_name))

    batch_summary_tsv, batch_summary_txt = write_batch_summary(args.outdir, batch_rows)
    trimmed_alignment, _ = step10_outputs(shared_layout)
    identity_tsv, identity_pdf, _ = step11_identity_outputs(shared_layout)
    treefile, _, tree_pdf = step12_tree_outputs(shared_layout)

    print()
    print_banner("Batch Completed", f"{len(batch_rows)} sample(s) processed successfully")
    print_status_line("Global summary", summarize_path(batch_summary_txt), "green")
    print_status_line("Global table", summarize_path(batch_summary_tsv), "green")
    print_status_line("Trimmed alignment", summarize_path(trimmed_alignment), "green")
    print_status_line("Identity matrix", summarize_path(identity_tsv), "green")
    print_status_line("Identity heatmap", summarize_path(identity_pdf), "green")
    print_status_line("Global tree", summarize_path(treefile), "green")
    print_status_line("Global tree PDF", summarize_path(tree_pdf), "green")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()
