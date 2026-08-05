#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ViraLong-ID v5.7
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
import os
import json
import re
import shutil
import subprocess
import sys
import textwrap
from datetime import date, datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Tuple

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


PIPELINE_NAME = "ViraLong-ID"
PIPELINE_VERSION = "5.6-batch"


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


def external_tool_env() -> Dict[str, str]:
    env = os.environ.copy()
    path_parts = [part for part in env.get("PATH", "").split(":") if part]
    for required_path in ["/usr/local/bin", "/usr/bin", "/bin", "/usr/sbin", "/sbin"]:
        if required_path not in path_parts:
            path_parts.append(required_path)
    env["PATH"] = ":".join(path_parts)
    return env


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
            env=external_tool_env(),
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


def count_fasta_records(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    return sum(1 for _ in SeqIO.parse(str(path), "fasta"))


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
    def __init__(self, prog):
        super().__init__(prog, max_help_position=34, width=110)

    def _get_help_string(self, action):
        help_text = action.help or ""
        if "%(default)" in help_text:
            return help_text
        if not action.option_strings or action.required:
            return help_text
        if action.default in (None, argparse.SUPPRESS):
            return help_text
        return f"{help_text} (default: %(default)s)"


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
        "beast2": base / "10_beast2_preparation",
        "beast2_run": base / "11_beast2_run",
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
        usage=(
            "ViraLong-ID.py --taxid TAXID --reads READS [READS ...] "
            "--outdir OUTDIR --refseq-virus-fasta FASTA [options]"
        ),
        formatter_class=LogoHelpFormatter,
        description=(
            "🧬 Long-read viral identification and phylogeny pipeline for one or more samples\n"
            "🦠 Shared references and BLAST database across the batch\n"
            "🌳 One combined alignment with robust strand correction, trimAl filtering, identity heatmap generation, and one global IQ-TREE phylogeny for all retained contigs"
        ),
        epilog=textwrap.dedent("""\
            Typical BEAST 2 workflow:
              1. First run with --prepare-beast2 to create editable templates.
              2. Fill manual_dates_template.tsv and map_locations_coordinates_template.tsv.
              3. Re-run with --beast2-manual-dates and --beast2-coordinates.
              4. Optionally add --run-beast2 --beast2-xml after creating/reviewing the final XML in BEAUti.
        """)
    )

    required = p.add_argument_group("Required inputs")
    required.add_argument("--taxid", required=True, help="🧬 Target NCBI Taxonomy ID")
    required.add_argument("--reads", required=True, nargs="+", type=Path,
                          help="📥 One or more input FASTQ / FASTQ.GZ files")
    required.add_argument("--outdir", required=True, type=Path, help="📂 Output directory")
    required.add_argument("--refseq-virus-fasta", required=True, type=Path,
                          help="🗃️ Local RefSeq virus FASTA for BLAST database")

    run_options = p.add_argument_group("Run options")
    run_options.add_argument("--sample-names", nargs="*",
                             help="🏷️ Optional sample names, same order as --reads")
    run_options.add_argument("--threads", type=int, default=8, help="🧵 Threads")

    read_qc = p.add_argument_group("Read QC and assembly")
    read_qc.add_argument("--min-q", type=int, default=15,
                         help="✨ Minimum mean read quality for fastplong")
    read_qc.add_argument("--flye-mode", choices=["normal", "meta"], default="meta",
                         help="🧱 Flye mode")
    read_qc.add_argument("--flye-iterations", type=int, default=1,
                         help="🔁 Flye polishing iterations")
    read_qc.add_argument("--assembly-min-q", type=float, default=20.0,
                         help="🔬 Minimum mean Q for reads retained for assembly")
    read_qc.add_argument("--assembly-min-len", default="auto",
                         help='📐 Minimum read length for assembly preselection, or "auto"')
    read_qc.add_argument("--assembly-max-len", default="auto",
                         help='📏 Maximum read length for assembly preselection, or "auto"')
    read_qc.add_argument("--assembly-target-cov", type=int, default=300,
                         help="🧮 Target effective coverage used to cap input for Flye")
    read_qc.add_argument("--assembly-retry-all-qc", action=argparse.BooleanOptionalAction, default=True,
                         help="🛟 Retry failed assemblies with all QC-passed reads when no target contigs are detected")

    target_selection = p.add_argument_group("Target contig selection")
    target_selection.add_argument("--min-pident", type=float, default=70.0,
                                  help="🎯 Minimum BLAST identity for target contig selection")
    target_selection.add_argument("--min-qcov", type=float, default=40.0,
                                  help="📏 Minimum BLAST query coverage for target contig selection")
    target_selection.add_argument("--min-contig-len-phylo", type=int, default=300,
                                  help="🌿 Minimum contig length retained for phylogeny")

    phylogeny = p.add_argument_group("Alignment, tree, and identity plots")
    phylogeny.add_argument("--trimal-gap-threshold", type=float, default=0.8,
                           help="✂️ trimAl gap threshold used to keep well-aligned columns")
    phylogeny.add_argument("--min-alignment-occupancy", type=float, default=0.7,
                           help="🧩 Minimum fraction of sequences with a residue required to keep an alignment column")
    phylogeny.add_argument("--min-assembled-core-occupancy", type=float, default=0.7,
                           help="🧬 Minimum fraction of assembled contigs with a residue required to define the shared core block")
    phylogeny.add_argument("--mafft-adjust-direction", choices=["off", "on", "accurate"], default="on",
                           help="🔄 Automatic strand correction in MAFFT")
    phylogeny.add_argument("--identity-plot-min", type=float, default=None,
                           help="🎨 Minimum value for pairwise identity heatmap color scale")

    beast2_prepare = p.add_argument_group("Optional BEAST 2 preparation")
    beast2_prepare.add_argument("--prepare-beast2", action="store_true",
                                help="🕰️ Prepare BEAST 2 input files after the global alignment")
    beast2_prepare.add_argument("--beast2-manual-dates", type=Path, default=None,
                                help="📝 Completed/editable TSV with manual sampling dates")
    beast2_prepare.add_argument("--beast2-coordinates", type=Path, default=None,
                                help="🗺️ Completed/editable TSV with latitude/longitude for map-ready outputs")

    beast2_run = p.add_argument_group("Optional BEAST 2 execution")
    beast2_run.add_argument("--run-beast2", action="store_true",
                            help="🚀 Run BEAST 2 after validating dates, coordinates, and final XML")
    beast2_run.add_argument("--beast2-xml", type=Path, default=None,
                            help="🧾 Optional custom BEAST 2 XML; if omitted, ViraLong writes a basic XML automatically")
    beast2_run.add_argument("--beast2-burnin", type=int, default=10,
                            help="🔥 TreeAnnotator burn-in percentage")
    beast2_run.add_argument("--beast2-chain-length", type=int, default=10_000_000,
                            help="⛓️ MCMC chain length for the automatically generated BEAST 2 XML")
    beast2_run.add_argument("--beast2-log-every", type=int, default=10_000,
                            help="🧾 Logging interval for the automatically generated BEAST 2 XML")
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
        proc = subprocess.run(cmd, stdout=out, stderr=logh, text=True, check=False, env=external_tool_env())
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


def reset_sample_assembly_outputs(sample_layout: Dict[str, Path]) -> None:
    shutil.rmtree(sample_layout["assembly"], ignore_errors=True)
    shutil.rmtree(sample_layout["blast"], ignore_errors=True)
    shutil.rmtree(sample_layout["target"], ignore_errors=True)
    mkdir(sample_layout["assembly"])
    mkdir(sample_layout["blast"])
    mkdir(sample_layout["target"])


def step6_assemble(sample_layout: Dict[str, Path], shared_layout: Dict[str, Path], threads: int, flye_mode: str,
                   reads_for_flye: Path | None = None, iterations: int = 1) -> None:
    log_file = sample_layout["logs"] / "step06_flye.log"
    renamed_refs = step2_outputs(shared_layout)[2]
    subset_fastq = reads_for_flye or step5_outputs(sample_layout)[0]
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
        "--iterations", str(iterations),
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
                                sample_name: str) -> None:
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
        if rec.id in keep_ids:
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


def rescue_sample_with_all_qc_reads(sample_layout: Dict[str, Path], shared_layout: Dict[str, Path], sample_name: str,
                                    min_q: int, threads: int, min_pident: float, min_qcov: float,
                                    flye_iterations: int) -> bool:
    log_file = sample_layout["logs"] / "step06_flye.log"
    qc_reads = step3_outputs(sample_layout, min_q)[0]
    if not fastq_output_usable(qc_reads):
        return False

    reset_sample_assembly_outputs(sample_layout)
    with open(log_file, "a", encoding="utf-8") as logh:
        logh.write(
            f"[{now()}] Rescue mode: retrying assembly with all QC-passed reads using Flye meta mode.\n"
        )

    step6_assemble(
        sample_layout,
        shared_layout,
        threads,
        "meta",
        reads_for_flye=qc_reads,
        iterations=flye_iterations,
    )
    step7_blast(sample_layout, shared_layout, threads)
    step8_select_target_contigs(sample_layout, shared_layout, min_pident, min_qcov, sample_name)
    return count_fasta_records(step8_outputs(sample_layout)[0]) > 0


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


def concatenate_fastas(files: List[Path], out_fasta: Path, min_len: int = 0) -> int:
    records = []
    for f in files:
        if not f.exists() or f.stat().st_size == 0:
            continue
        for rec in SeqIO.parse(str(f), "fasta"):
            if len(rec.seq) < min_len:
                continue
            records.append(rec)
    SeqIO.write(records, str(out_fasta), "fasta")
    return len(records)


def step9_collect_and_align(shared_layout: Dict[str, Path], sample_layouts: List[Dict[str, Path]], threads: int,
                            adjust_direction: str, min_contig_len_phylo: int) -> None:
    log_file = shared_layout["logs"] / "step09_mafft.log"
    renamed_refs = step2_outputs(shared_layout)[2]
    combined_contigs, direction_guide, oriented_contigs, dataset, aln = step9_outputs(shared_layout)
    ref_aln = shared_layout["aln"] / "reference_alignment.fasta"

    target_fastas = [step8_outputs(sample_layout)[0] for sample_layout in sample_layouts]
    total_target_contigs = sum(
        1 for fasta in target_fastas if fasta.exists() for _ in SeqIO.parse(str(fasta), "fasta")
    )
    n_contigs = concatenate_fastas(target_fastas, combined_contigs, min_len=min_contig_len_phylo)
    if n_contigs == 0:
        raise RuntimeError(
            "No target contigs passed the phylogeny length filter across the batch. "
            "Try lowering --min-contig-len-phylo."
        )
    with open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"[{now()}] Target contigs detected across batch: {total_target_contigs}\n")
        logh.write(f"[{now()}] Target contigs retained for phylogeny (>= {min_contig_len_phylo} bp): {n_contigs}\n")

    cmd1 = [
        "mafft",
        "--thread", str(threads),
        *mafft_direction_flag(adjust_direction),
        "--auto",
        str(renamed_refs)
    ]
    with open(ref_aln, "w", encoding="utf-8") as out, open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"\n[{now()}] CMD: {' '.join(cmd1)}\n")
        proc = subprocess.run(cmd1, stdout=out, stderr=logh, text=True, check=False, env=external_tool_env())
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
            proc = subprocess.run(cmd_orient, stdout=out, stderr=logh, text=True, check=False, env=external_tool_env())
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
        proc = subprocess.run(cmd2, stdout=out, stderr=logh, text=True, check=False, env=external_tool_env())
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


def filter_alignment_by_occupancy(in_fasta: Path, out_fasta: Path, min_occupancy: float) -> Tuple[int, int]:
    records = list(SeqIO.parse(str(in_fasta), "fasta"))
    if not records:
        raise RuntimeError("trimAl produced an empty alignment")

    aln_len = len(records[0].seq)
    if aln_len == 0:
        raise RuntimeError("trimAl produced an alignment with zero columns")
    if any(len(rec.seq) != aln_len for rec in records):
        raise RuntimeError("Alignment records do not all have the same length")

    keep_cols = []
    n_records = len(records)
    for idx in range(aln_len):
        occupied = 0
        for rec in records:
            base = rec.seq[idx]
            if base not in "-?":
                occupied += 1
        if (occupied / n_records) >= min_occupancy:
            keep_cols.append(idx)

    if not keep_cols:
        raise RuntimeError(
            "Alignment occupancy filter removed all columns. Try lowering --min-alignment-occupancy."
        )

    filtered_records = []
    for rec in records:
        seq = "".join(rec.seq[i] for i in keep_cols)
        rec.seq = rec.seq.__class__(seq)
        filtered_records.append(rec)

    SeqIO.write(filtered_records, str(out_fasta), "fasta")
    return aln_len, len(keep_cols)


def extract_assembled_core_block(in_fasta: Path, out_fasta: Path, min_occupancy: float) -> Tuple[int, int, int]:
    records = list(SeqIO.parse(str(in_fasta), "fasta"))
    if not records:
        raise RuntimeError("trimAl produced an empty alignment")

    assembled_records = [rec for rec in records if "__" in rec.id]
    if not assembled_records:
        raise RuntimeError("No assembled contigs were found to define the core alignment block")

    aln_len = len(records[0].seq)
    if aln_len == 0:
        raise RuntimeError("trimAl produced an alignment with zero columns")
    if any(len(rec.seq) != aln_len for rec in records):
        raise RuntimeError("Alignment records do not all have the same length")

    keep_mask = []
    n_assembled = len(assembled_records)
    for idx in range(aln_len):
        occupied = 0
        for rec in assembled_records:
            base = rec.seq[idx]
            if base not in "-?":
                occupied += 1
        keep_mask.append((occupied / n_assembled) >= min_occupancy)

    best_start = None
    best_end = None
    run_start = None
    for idx, keep in enumerate(keep_mask):
        if keep and run_start is None:
            run_start = idx
        elif not keep and run_start is not None:
            run_end = idx - 1
            if best_start is None or (run_end - run_start) > (best_end - best_start):
                best_start, best_end = run_start, run_end
            run_start = None
    if run_start is not None:
        run_end = aln_len - 1
        if best_start is None or (run_end - run_start) > (best_end - best_start):
            best_start, best_end = run_start, run_end

    if best_start is None or best_end is None:
        raise RuntimeError(
            "The assembled-contig core filter removed all columns. Try lowering --min-assembled-core-occupancy."
        )

    cropped_records = []
    for rec in records:
        seq = str(rec.seq[best_start:best_end + 1])
        rec.seq = rec.seq.__class__(seq)
        cropped_records.append(rec)

    SeqIO.write(cropped_records, str(out_fasta), "fasta")
    return aln_len, best_start + 1, best_end + 1


def step10_trimal(shared_layout: Dict[str, Path], gap_threshold: float, min_occupancy: float,
                  min_assembled_core_occupancy: float) -> None:
    log_file = shared_layout["logs"] / "step10_trimal.log"
    aln = step9_outputs(shared_layout)[4]
    trimmed, html = step10_outputs(shared_layout)
    trimal_raw = shared_layout["aln"] / "alignment_mafft.trimal_raw.fasta"
    core_raw = shared_layout["aln"] / "alignment_mafft.core_raw.fasta"

    run_logged(
        [
            "trimal",
            "-in", str(aln),
            "-out", str(trimal_raw),
            "-htmlout", str(html),
            "-gt", str(gap_threshold),
        ],
        log_file
    )

    if not trimal_raw.exists() or trimal_raw.stat().st_size == 0:
        raise RuntimeError("trimAl did not produce a trimmed alignment")
    trimal_cols, core_start, core_end = extract_assembled_core_block(
        trimal_raw, core_raw, min_assembled_core_occupancy
    )
    original_cols, kept_cols = filter_alignment_by_occupancy(core_raw, trimmed, min_occupancy)
    with open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"[{now()}] Columns after trimAl: {trimal_cols}\n")
        logh.write(
            f"[{now()}] Assembled-contig core block (>= {min_assembled_core_occupancy:.2f}) "
            f"retained columns {core_start}-{core_end}\n"
        )
        logh.write(f"[{now()}] Columns after occupancy filter (>= {min_occupancy:.2f}): {kept_cols}\n")


# ---------------------------------------------------------------------
# Optional BEAST 2 preparation
# ---------------------------------------------------------------------

def beast2_outputs(shared_layout: Dict[str, Path]) -> Dict[str, Path]:
    outdir = shared_layout["beast2"]
    return {
        "outdir": outdir,
        "alignment_fasta": outdir / "alignment_beast2_safe_ids.fasta",
        "alignment_nexus": outdir / "alignment_beast2_safe_ids.nexus",
        "metadata": outdir / "metadata_beast2.tsv",
        "tip_dates": outdir / "tip_dates_beast2.tsv",
        "manual_dates_template": outdir / "manual_dates_template.tsv",
        "traits": outdir / "traits_beauti.tsv",
        "locations_template": outdir / "map_locations_coordinates_template.tsv",
        "sequence_coordinates": outdir / "sequence_coordinates_template.tsv",
        "id_map": outdir / "sequence_id_map.tsv",
        "xml_template": outdir / "CYVCV_BEAST2_template.xml",
        "readme": outdir / "README_BEAST2_preparacion.md",
    }


def beast2_done(shared_layout: Dict[str, Path]) -> bool:
    outputs = beast2_outputs(shared_layout)
    required = [
        "alignment_fasta",
        "alignment_nexus",
        "metadata",
        "tip_dates",
        "manual_dates_template",
        "traits",
        "locations_template",
        "sequence_coordinates",
        "id_map",
        "xml_template",
        "readme",
    ]
    return all(outputs[name].exists() and outputs[name].stat().st_size > 0 for name in required)


def read_ncbi_jsonl_metadata(jsonl: Path) -> Dict[str, Dict[str, str]]:
    metadata: Dict[str, Dict[str, str]] = {}
    if not jsonl.exists():
        return metadata
    with open(jsonl, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            data = json.loads(line)
            accession = data.get("accession", "")
            isolate = data.get("isolate") or {}
            host = data.get("host") or {}
            location = data.get("location") or {}
            metadata[accession] = {
                "accession": accession,
                "isolate": isolate.get("name", ""),
                "host": host.get("organismName", ""),
                "country_locality": location.get("geographicLocation", ""),
                "region": location.get("geographicRegion", ""),
                "collection_date": isolate.get("collectionDate", ""),
                "sample_source": isolate.get("source", ""),
                "completeness": data.get("completeness", ""),
                "length": str(data.get("length", "")),
                "release_date": (data.get("releaseDate", "") or "")[:10],
                "update_date": (data.get("updateDate", "") or "")[:10],
            }
    return metadata


def split_country_locality(value: str) -> Tuple[str, str]:
    text = (value or "").replace("-_", ": ").replace("_", " ").replace("-", " ")
    text = re.sub(r"\s+", " ", text).strip()
    if not text:
        return "", ""
    if ":" in text:
        country, locality = text.split(":", 1)
        locality = locality.replace(",", "_").strip().replace(" ", "_")
        return country.strip(), locality
    parts = text.split()
    two_word_countries = {"South Korea"}
    if len(parts) >= 2 and " ".join(parts[:2]) in two_word_countries:
        return " ".join(parts[:2]), "_".join(parts[2:])
    return parts[0], "_".join(parts[1:])


def decimal_sampling_date(raw: str) -> Tuple[str, str]:
    value = (raw or "").strip()
    if not value:
        return "", ""
    match = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})", value)
    if match:
        year, month, day = map(int, match.groups())
        sample_date = date(year, month, day)
        start = date(year, 1, 1)
        end = date(year + 1, 1, 1)
        return f"{year + (sample_date - start).days / (end - start).days:.6f}", "day"
    match = re.fullmatch(r"(\d{4})-(\d{2})", value)
    if match:
        year, month = map(int, match.groups())
        sample_date = date(year, month, 15)
        start = date(year, 1, 1)
        end = date(year + 1, 1, 1)
        return f"{year + (sample_date - start).days / (end - start).days:.6f}", "month_midpoint"
    if re.fullmatch(r"\d{4}", value):
        return value, "year_only"
    return "", "unparsed"


def parse_beast2_header(header: str) -> Dict[str, str]:
    if "/" in header:
        parts = header.split("/")
        return {
            "record_type": "reference",
            "accession": parts[0],
            "isolate": parts[1] if len(parts) > 1 else "",
            "host": (parts[2] if len(parts) > 2 else "").replace("_", " "),
            "country_locality": (parts[3] if len(parts) > 3 else "").replace("-_", ": ").replace("_", " "),
        }

    sample = header.split("__contig_")[0]
    contig = f"contig_{header.split('__contig_')[1]}" if "__contig_" in header else ""
    tokens = sample.split("_")
    sample_id = sample
    country = ""
    region = ""
    locality = ""
    sample_type = ""

    if "Spain" in tokens:
        idx = tokens.index("Spain")
        sample_id = "_".join(tokens[:idx])
        country = "Spain"
        rest = tokens[idx + 1:]
        region = rest[0] if rest else ""
        if rest and rest[-1] in {"Campo", "Vivero", "Germoplasma"}:
            sample_type = rest[-1]
            locality = "_".join(rest[1:-1])
        else:
            locality = "_".join(rest[1:])

    if sample_id.endswith("_PCR"):
        sample_id = sample_id[:-4]
        if not sample_type:
            sample_type = "PCR"

    return {
        "record_type": "local",
        "sample": sample,
        "sample_id": sample_id,
        "accession": "",
        "isolate": sample_id,
        "host": "Citrus",
        "country": country,
        "region": region,
        "locality": locality,
        "sample_type": sample_type,
        "contig": contig,
    }


def beast2_safe_sequence_id(header: str, ncbi_metadata: Dict[str, Dict[str, str]], seen: set[str]) -> str:
    parsed = parse_beast2_header(header)
    if parsed["record_type"] == "reference":
        accession = parsed["accession"]
        isolate = ncbi_metadata.get(accession, {}).get("isolate") or parsed.get("isolate") or accession
        sequence_id = f"{accession}_{isolate}"
    else:
        sequence_id = parsed.get("sample_id") or parsed.get("sample") or "sequence"
    sequence_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", sequence_id).strip("_") or "sequence"
    base_id = sequence_id
    suffix = 2
    while sequence_id in seen:
        sequence_id = f"{base_id}_{suffix}"
        suffix += 1
    seen.add(sequence_id)
    return sequence_id


def read_manual_beast2_dates(path: Path | None) -> Dict[str, str]:
    dates: Dict[str, str] = {}
    if path is None:
        return dates
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sequence_id = (row.get("sequence_id") or "").strip()
            if not sequence_id:
                continue
            exact_date = (row.get("collection_date_YYYY_MM_DD") or row.get("collection_date") or "").strip()
            year = (row.get("year") or "").strip()
            dates[sequence_id] = exact_date or year
    return dates


def read_beast2_coordinates(path: Path | None) -> Dict[str, Dict[str, str]]:
    coordinates: Dict[str, Dict[str, str]] = {}
    if path is None:
        return coordinates
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            location_id = (row.get("location_id") or "").strip()
            location_label = (row.get("location_label") or "").strip()
            payload = {
                "latitude": (row.get("latitude") or "").strip(),
                "longitude": (row.get("longitude") or "").strip(),
                "coordinate_precision": (row.get("coordinate_precision") or "").strip(),
                "coordinate_source": (row.get("coordinate_source") or "").strip(),
            }
            if location_id:
                coordinates[location_id] = payload
            if location_label:
                coordinates[location_label] = payload
    return coordinates


def beast2_location_label(row: Dict[str, str]) -> str:
    parts = [row.get("country", ""), row.get("region", ""), row.get("locality", "")]
    parts = [part for part in parts if part and part != "NA"]
    return " | ".join(parts) if parts else "Unknown"


def write_tsv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def step13_prepare_beast2(shared_layout: Dict[str, Path], taxid: str,
                          manual_dates_file: Path | None = None,
                          coordinates_file: Path | None = None) -> None:
    outputs = beast2_outputs(shared_layout)
    mkdir(outputs["outdir"])
    alignment = step10_outputs(shared_layout)[0]
    _, _, jsonl = step1_outputs(shared_layout, taxid)
    ncbi_metadata = read_ncbi_jsonl_metadata(jsonl)
    manual_dates = read_manual_beast2_dates(manual_dates_file)
    coordinate_overrides = read_beast2_coordinates(coordinates_file)

    records = list(SeqIO.parse(str(alignment), "fasta"))
    if not records:
        raise RuntimeError("No sequences found in trimmed alignment for BEAST 2 preparation")
    alignment_length = len(records[0].seq)
    if any(len(rec.seq) != alignment_length for rec in records):
        raise RuntimeError("Trimmed alignment is not rectangular; BEAST 2 preparation stopped")

    seen: set[str] = set()
    metadata_rows: List[Dict[str, object]] = []
    id_map_rows: List[Dict[str, object]] = []
    safe_records = []

    for rec in records:
        original_header = rec.id
        sequence_id = beast2_safe_sequence_id(original_header, ncbi_metadata, seen)
        parsed = parse_beast2_header(original_header)
        row: Dict[str, object] = {
            "sequence_id": sequence_id,
            "original_header": original_header,
            "record_type": parsed["record_type"],
            "accession": "",
            "isolate": "",
            "host": "",
            "country": "",
            "region": "",
            "locality": "",
            "sample_type": "",
            "collection_date": "",
            "decimal_date": "",
            "date_precision": "",
            "date_source": "",
            "needs_manual_date": "",
            "sequence_length": len(rec.seq),
            "ungapped_length": len(str(rec.seq).replace("-", "")),
        }

        if parsed["record_type"] == "reference":
            accession = parsed["accession"]
            ref_meta = ncbi_metadata.get(accession, {})
            country, locality = split_country_locality(ref_meta.get("country_locality") or parsed.get("country_locality", ""))
            collection_date = ref_meta.get("collection_date", "")
            date_source = "NCBI collectionDate"
            if sequence_id in manual_dates:
                collection_date = manual_dates[sequence_id]
                date_source = "manual_dates_file"
            decimal_date, precision = decimal_sampling_date(collection_date)
            row.update({
                "accession": accession,
                "isolate": ref_meta.get("isolate") or parsed.get("isolate", ""),
                "host": ref_meta.get("host") or parsed.get("host", ""),
                "country": country,
                "region": ref_meta.get("region", ""),
                "locality": locality,
                "sample_type": ref_meta.get("sample_source", ""),
                "collection_date": collection_date,
                "decimal_date": decimal_date,
                "date_precision": precision,
                "date_source": date_source if decimal_date else "missing_or_unparsed",
                "needs_manual_date": "NO" if decimal_date else "YES",
            })
        else:
            collection_date = manual_dates.get(sequence_id, "")
            decimal_date, precision = decimal_sampling_date(collection_date)
            row.update({
                "accession": parsed.get("sample_id", ""),
                "isolate": parsed.get("sample_id", ""),
                "host": parsed.get("host", ""),
                "country": parsed.get("country", ""),
                "region": parsed.get("region", ""),
                "locality": parsed.get("locality", ""),
                "sample_type": parsed.get("sample_type", ""),
                "collection_date": collection_date,
                "decimal_date": decimal_date,
                "date_precision": precision,
                "date_source": "manual_dates_file" if decimal_date else "manual_required_local_sample",
                "needs_manual_date": "NO" if decimal_date else "YES",
            })

        rec.id = sequence_id
        rec.name = sequence_id
        rec.description = ""
        safe_records.append(rec)
        metadata_rows.append(row)
        id_map_rows.append({"sequence_id": sequence_id, "original_header": original_header})

    SeqIO.write(safe_records, str(outputs["alignment_fasta"]), "fasta")
    with open(outputs["alignment_nexus"], "w", encoding="utf-8") as fh:
        fh.write("#NEXUS\n\n")
        fh.write("BEGIN DATA;\n")
        fh.write(f"  DIMENSIONS NTAX={len(safe_records)} NCHAR={alignment_length};\n")
        fh.write("  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        fh.write("  MATRIX\n")
        for rec in safe_records:
            fh.write(f"  {rec.id} {str(rec.seq)}\n")
        fh.write("  ;\nEND;\n")

    metadata_fields = [
        "sequence_id", "original_header", "record_type", "accession", "isolate", "host",
        "country", "region", "locality", "sample_type", "collection_date", "decimal_date",
        "date_precision", "date_source", "needs_manual_date", "sequence_length", "ungapped_length",
    ]
    write_tsv(outputs["metadata"], metadata_rows, metadata_fields)
    write_tsv(
        outputs["tip_dates"],
        metadata_rows,
        ["sequence_id", "collection_date", "decimal_date", "date_precision", "date_source", "needs_manual_date"],
    )
    write_tsv(outputs["id_map"], id_map_rows, ["sequence_id", "original_header"])

    manual_rows = [
        {
            "sequence_id": row["sequence_id"],
            "sample_id_or_accession": row["accession"] or row["isolate"],
            "current_collection_date": row["collection_date"],
            "year": "",
            "collection_date_YYYY_MM_DD": "",
            "notes": "",
        }
        for row in metadata_rows
        if row["needs_manual_date"] == "YES"
    ]
    manual_template_out = outputs["manual_dates_template"]
    if manual_dates_file is not None:
        try:
            same_manual_path = manual_dates_file.resolve() == manual_template_out.resolve()
        except OSError:
            same_manual_path = False
        if same_manual_path:
            manual_template_out = outputs["outdir"] / "manual_dates_missing_after_import.tsv"
    write_tsv(
        manual_template_out,
        manual_rows,
        ["sequence_id", "sample_id_or_accession", "current_collection_date", "year", "collection_date_YYYY_MM_DD", "notes"],
    )

    traits_rows = [
        {
            "sequence_id": row["sequence_id"],
            "country": row["country"] or "NA",
            "host": (row["host"] or "NA").replace(" ", "_"),
            "region": row["region"] or "NA",
            "sample_type": row["sample_type"] or "NA",
        }
        for row in metadata_rows
    ]
    write_tsv(outputs["traits"], traits_rows, ["sequence_id", "country", "host", "region", "sample_type"])

    locations: Dict[str, Dict[str, object]] = {}
    for row in metadata_rows:
        label = beast2_location_label(row)
        if label not in locations:
            location_id = f"loc_{len(locations) + 1:03d}"
            coords = coordinate_overrides.get(location_id) or coordinate_overrides.get(label) or {}
            locations[label] = {
                "location_id": location_id,
                "location_label": label,
                "country": row["country"],
                "region": row["region"],
                "locality": row["locality"],
                "latitude": coords.get("latitude", ""),
                "longitude": coords.get("longitude", ""),
                "coordinate_precision": coords.get("coordinate_precision", ""),
                "coordinate_source": coords.get("coordinate_source", ""),
                "notes": "Fill manually or geocode from country/region/locality before continuous phylogeography/map rendering.",
            }

    locations_template_out = outputs["locations_template"]
    if coordinates_file is not None:
        try:
            same_coordinates_path = coordinates_file.resolve() == locations_template_out.resolve()
        except OSError:
            same_coordinates_path = False
        if same_coordinates_path:
            locations_template_out = outputs["outdir"] / "map_locations_missing_after_import.tsv"
    write_tsv(
        locations_template_out,
        list(locations.values()),
        [
            "location_id", "location_label", "country", "region", "locality", "latitude",
            "longitude", "coordinate_precision", "coordinate_source", "notes",
        ],
    )

    sequence_coordinate_rows = []
    for row in metadata_rows:
        label = beast2_location_label(row)
        location = locations[label]
        sequence_coordinate_rows.append({
            "sequence_id": row["sequence_id"],
            "location_id": location["location_id"],
            "location_label": label,
            "latitude": location["latitude"],
            "longitude": location["longitude"],
            "needs_coordinates": "NO" if location["latitude"] and location["longitude"] else "YES",
        })
    write_tsv(
        outputs["sequence_coordinates"],
        sequence_coordinate_rows,
        ["sequence_id", "location_id", "location_label", "latitude", "longitude", "needs_coordinates"],
    )

    missing_dates = sum(1 for row in metadata_rows if row["needs_manual_date"] == "YES")
    missing_locations = sum(1 for row in locations.values() if not row["latitude"] or not row["longitude"])
    with open(outputs["xml_template"], "w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
        fh.write("<!--\n")
        fh.write("TEMPLATE ONLY - not ready to run in BEAST 2 until required dates are complete.\n")
        fh.write("Import alignment_beast2_safe_ids.fasta and tip_dates_beast2.tsv in BEAUti to create the final runnable XML.\n")
        fh.write(f"Current missing sampling dates: {missing_dates}.\n")
        fh.write(f"Current locations lacking coordinates: {missing_locations}.\n")
        fh.write("-->\n")
        fh.write('<beast beautitemplate="Standard" namespace="beast.base.core:beast.base.evolution.alignment:beast.base.evolution.tree" version="2.7">\n')
        fh.write('  <data id="CYVCV_alignment" name="alignment">\n')
        for rec in safe_records:
            fh.write(f'    <sequence id="seq_{rec.id}" taxon="{rec.id}" totalcount="4" value="{str(rec.seq)}"/>\n')
        fh.write("  </data>\n")
        fh.write("</beast>\n")

    readme = f"""# BEAST 2 preparation

Generated from ViraLong-ID optional `--prepare-beast2` mode.

## Main files

- `alignment_beast2_safe_ids.fasta`: trimmed alignment with BEAST-safe sequence IDs.
- `alignment_beast2_safe_ids.nexus`: same alignment in NEXUS format.
- `metadata_beast2.tsv`: automatically extracted metadata.
- `tip_dates_beast2.tsv`: sampling dates for BEAUti/BEAST.
- `manual_dates_template.tsv`: editable file with sequences still requiring sampling dates.
- `traits_beauti.tsv`: optional discrete traits for BEAUti.
- `map_locations_coordinates_template.tsv`: editable file with one row per unique location.
- `sequence_coordinates_template.tsv`: per-sequence coordinate table for map rendering.
- `sequence_id_map.tsv`: original headers mapped to BEAST-safe IDs.
- `CYVCV_BEAST2_template.xml`: non-runnable template until dates are complete.

## Status

- Sequences: {len(safe_records)}
- Alignment length: {alignment_length}
- Missing sampling dates: {missing_dates}
- Unique locations: {len(locations)}
- Locations without coordinates: {missing_locations}

## Required user input

1. Fill `manual_dates_template.tsv` with either `year` or `collection_date_YYYY_MM_DD`.
2. Fill `map_locations_coordinates_template.tsv` with `latitude` and `longitude` for map-ready phylogeography.

Re-run ViraLong-ID with `--prepare-beast2 --beast2-manual-dates <file> --beast2-coordinates <file>` to incorporate edited dates and coordinates.

To execute BEAST 2 from ViraLong-ID, complete the editable date and coordinate tables, then run with:

`--run-beast2 --beast2-manual-dates <file> --beast2-coordinates <file>`

If `--beast2-xml` is omitted, ViraLong-ID writes a basic automatic XML at `CYVCV_BEAST2.xml`.
Use `--beast2-xml <final.xml>` only when you want to run a custom XML generated/reviewed elsewhere.
"""
    with open(outputs["readme"], "w", encoding="utf-8") as fh:
        fh.write(readme)


def beast2_run_outputs(shared_layout: Dict[str, Path]) -> Dict[str, Path]:
    outdir = shared_layout["beast2_run"]
    return {
        "outdir": outdir,
        "log": outdir / "beast2_run.log",
        "mcc_tree": outdir / "CYVCV_BEAST2.MCC.tree",
        "map_pdf": outdir / "CYVCV_BEAST2_phylogeography_map.pdf",
        "map_png": outdir / "CYVCV_BEAST2_phylogeography_map.png",
        "map_html": outdir / "CYVCV_BEAST2_phylogeography_map.html",
        "time_tree_pdf": outdir / "CYVCV_BEAST2_time_tree.pdf",
        "time_tree_png": outdir / "CYVCV_BEAST2_time_tree.png",
        "map_edges_tsv": outdir / "CYVCV_BEAST2_phylogeography_edges.tsv",
        "visuals_version": outdir / "CYVCV_BEAST2_visuals_v6.done",
    }


def beast2_run_done(shared_layout: Dict[str, Path]) -> bool:
    outputs = beast2_run_outputs(shared_layout)
    required = ["mcc_tree", "map_pdf", "map_png", "map_html", "time_tree_pdf", "time_tree_png", "map_edges_tsv", "visuals_version"]
    return all(outputs[name].exists() and outputs[name].stat().st_size > 0 for name in required)


def validate_beast2_ready(shared_layout: Dict[str, Path]) -> None:
    outputs = beast2_outputs(shared_layout)
    missing_dates = []
    with open(outputs["tip_dates"], "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if (row.get("needs_manual_date") or "").strip().upper() == "YES" or not (row.get("decimal_date") or "").strip():
                missing_dates.append(row.get("sequence_id", "unknown"))

    missing_coordinates = []
    with open(outputs["sequence_coordinates"], "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if (row.get("needs_coordinates") or "").strip().upper() == "YES":
                missing_coordinates.append(row.get("location_label") or row.get("sequence_id") or "unknown")

    if missing_dates:
        preview = ", ".join(missing_dates[:8])
        raise RuntimeError(
            f"BEAST 2 cannot run: {len(missing_dates)} sequence(s) still lack sampling dates. "
            f"Examples: {preview}"
        )
    if missing_coordinates:
        unique_preview = ", ".join(sorted(set(missing_coordinates))[:8])
        raise RuntimeError(
            f"BEAST 2 map-ready run cannot start: {len(set(missing_coordinates))} location(s) still lack coordinates. "
            f"Examples: {unique_preview}"
        )


def xml_escape(value: str) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_auto_beast2_xml(shared_layout: Dict[str, Path], xml_path: Path,
                          chain_length: int, log_every: int) -> None:
    outputs = beast2_outputs(shared_layout)
    records = list(SeqIO.parse(str(outputs["alignment_fasta"]), "fasta"))
    if not records:
        raise RuntimeError("Cannot write automatic BEAST 2 XML: no alignment records found")

    tip_dates: Dict[str, str] = {}
    with open(outputs["tip_dates"], "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sequence_id = (row.get("sequence_id") or "").strip()
            decimal_date = (row.get("decimal_date") or "").strip()
            if sequence_id and decimal_date:
                tip_dates[sequence_id] = decimal_date

    missing = [rec.id for rec in records if rec.id not in tip_dates]
    if missing:
        preview = ", ".join(missing[:8])
        raise RuntimeError(f"Cannot write automatic BEAST 2 XML: missing dates for {preview}")

    trait_value = ",\n                ".join(
        f"{xml_escape(rec.id)} = {tip_dates[rec.id]}" for rec in records
    )

    with open(xml_path, "w", encoding="utf-8") as fh:
        fh.write("""<?xml version="1.0" encoding="UTF-8"?>
<beast version='2.0'
       namespace='beast.base.math:beast.base.evolution.alignment:beast.pkgmgmt:beast.base.core:beast.base.inference:beast.base.inference.distribution:beast.base.evolution.tree.coalescent:beast.base.inference.util:beast.evolution.nuc:beast.base.evolution.operator:beast.base.inference.operator:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood'>

    <!-- Automatic ViraLong-ID BEAST 2 XML.
         Basic first-pass configuration: HKY, strict clock, constant coalescent population.
         Coordinates are validated/exported for map post-processing but are not modelled in this basic XML. -->

    <data id="alignment" dataType="nucleotide">
""")
        for rec in records:
            fh.write(f'        <sequence taxon="{xml_escape(rec.id)}">\n')
            fh.write(f"            {str(rec.seq)}\n")
            fh.write("        </sequence>\n")
        fh.write(f"""    </data>

    <input spec='HKY' id='hky'>
        <parameter name='kappa' idref='hky.kappa'/>
        <input id='freqs' name='frequencies' spec='Frequencies'>
            <input name='data' idref='alignment'/>
        </input>
    </input>

    <input spec='SiteModel' id="siteModel">
        <input name='substModel' idref='hky'/>
    </input>

    <parameter id="hky.kappa" value="2.0" lower="0.0"/>
    <input spec='beast.base.inference.parameter.RealParameter' id="popSize">1.0</input>

    <tree spec='beast.base.evolution.tree.ClusterTree' id='tree' clusterType='upgma'>
        <trait spec='beast.base.evolution.tree.TraitSet' traitname='date-forward' units='year'
               value='
                {trait_value}
               '>
            <taxa spec='TaxonSet' alignment='@alignment'/>
        </trait>
        <input name='taxa' idref='alignment'/>
    </tree>

    <input spec='TreeLikelihood' id="treeLikelihood">
        <input name='data' idref="alignment"/>
        <input name='tree' idref="tree"/>
        <input name='siteModel' idref="siteModel"/>
        <branchRateModel id="StrictClockModel" spec="beast.base.evolution.branchratemodel.StrictClockModel">
            <parameter dimension="1" estimate="false" id="clockRate" name="clock.rate" value="1.0E-4"/>
        </branchRateModel>
    </input>

    <run spec="MCMC" id="mcmc" chainLength="{chain_length}">
        <state>
            <input name='stateNode' idref="hky.kappa"/>
            <input name='stateNode' idref="popSize"/>
            <input name='stateNode' idref="tree"/>
        </state>

        <distribution spec="CompoundDistribution" id="posterior">
            <distribution id="coalescent" spec="Coalescent">
                <treeIntervals spec='beast.base.evolution.tree.TreeIntervals' id='TreeIntervals'>
                     <tree idref="tree"/>
                </treeIntervals>
                <populationModel spec="ConstantPopulation" id='ConstantPopulation'>
                    <parameter name="popSize" idref="popSize"/>
                </populationModel>
            </distribution>
            <distribution id='likelihood' idref="treeLikelihood"/>
        </distribution>

        <operator id='kappaScaler' spec='ScaleOperator' scaleFactor="0.5" weight="1">
            <parameter idref="hky.kappa"/>
        </operator>
        <operator id='popSizeScaler' spec='ScaleOperator' scaleFactor="0.5" weight="1">
            <parameter idref="popSize"/>
        </operator>
        <operator id='treeScaler' spec='ScaleOperator' scaleFactor="0.5" weight="1">
            <tree idref="tree"/>
        </operator>
        <operator spec='beast.base.evolution.operator.Uniform' weight="10">
            <tree idref="tree"/>
        </operator>
        <operator spec='SubtreeSlide' weight="5" gaussian="true" size="1.0">
            <tree idref="tree"/>
        </operator>
        <operator id='narrow' spec='Exchange' isNarrow='true' weight="1">
            <tree idref="tree"/>
        </operator>
        <operator id='wide' spec='Exchange' isNarrow='false' weight="1">
            <tree idref="tree"/>
        </operator>
        <operator spec='WilsonBalding' weight="1">
            <tree idref="tree"/>
        </operator>

        <logger logEvery="{log_every}" fileName="CYVCV_BEAST2.log">
            <model idref='posterior'/>
            <log idref="coalescent"/>
            <log idref="likelihood"/>
            <log idref="popSize"/>
            <log idref="hky.kappa"/>
            <log idref="posterior"/>
            <log spec='beast.base.evolution.tree.TreeHeightLogger' tree='@tree'/>
        </logger>
        <logger logEvery="{log_every}" fileName="CYVCV_BEAST2.trees">
            <log idref="tree"/>
        </logger>
        <logger logEvery="{log_every}">
            <model idref='posterior'/>
            <log idref="coalescent"/>
            <log idref="likelihood"/>
            <log idref="popSize"/>
            <log idref="hky.kappa"/>
            <log idref="posterior"/>
        </logger>
    </run>

</beast>
""")


def resolve_beast2_xml(shared_layout: Dict[str, Path], beast2_xml: Path | None,
                       chain_length: int, log_every: int) -> Path:
    if beast2_xml is not None:
        return beast2_xml
    default_xml = beast2_outputs(shared_layout)["outdir"] / "CYVCV_BEAST2.xml"
    write_auto_beast2_xml(shared_layout, default_xml, chain_length, log_every)
    return default_xml


def newest_tree_file(search_dirs: List[Path]) -> Path | None:
    candidates: List[Path] = []
    for directory in search_dirs:
        if directory.exists():
            candidates.extend(directory.glob("*.trees"))
            candidates.extend(directory.glob("*.trees.gz"))
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def step14_run_beast2(shared_layout: Dict[str, Path], beast2_xml: Path | None, threads: int,
                       burnin: int, chain_length: int, log_every: int) -> None:
    validate_beast2_ready(shared_layout)
    xml = resolve_beast2_xml(shared_layout, beast2_xml, chain_length, log_every)
    outputs = beast2_run_outputs(shared_layout)
    mkdir(outputs["outdir"])

    run_xml = outputs["outdir"] / xml.name
    if xml.resolve() != run_xml.resolve():
        shutil.copy2(xml, run_xml)

    if not outputs["mcc_tree"].exists() or outputs["mcc_tree"].stat().st_size == 0:
        run_logged(
            ["beast", "-overwrite", "-threads", str(threads), str(run_xml)],
            outputs["log"],
            cwd=outputs["outdir"],
        )

        tree_file = newest_tree_file([outputs["outdir"], run_xml.parent])
        if tree_file is None:
            raise RuntimeError("BEAST 2 finished but no .trees output was found for TreeAnnotator")

        run_logged(
            ["treeannotator", "-burnin", str(burnin), str(tree_file), str(outputs["mcc_tree"])],
            outputs["log"],
            cwd=outputs["outdir"],
        )
    else:
        with open(outputs["log"], "a", encoding="utf-8") as logh:
            logh.write(f"[{now()}] Existing MCC tree found; skipping BEAST 2 and TreeAnnotator rerun.\n")

    if not outputs["mcc_tree"].exists() or outputs["mcc_tree"].stat().st_size == 0:
        raise RuntimeError("TreeAnnotator did not produce the expected MCC tree")

    render_beast2_visual_outputs(shared_layout)



def beast2_parse_float(value: str | None) -> float | None:
    text = (value or "").strip().replace(",", ".")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def read_beast2_tip_dates_table(path: Path) -> Dict[str, float]:
    dates: Dict[str, float] = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sequence_id = (row.get("sequence_id") or "").strip()
            decimal_date = beast2_parse_float(row.get("decimal_date"))
            if sequence_id and decimal_date is not None:
                dates[sequence_id] = decimal_date
    return dates


def read_beast2_metadata_table(path: Path) -> Dict[str, Dict[str, str]]:
    rows: Dict[str, Dict[str, str]] = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sequence_id = (row.get("sequence_id") or "").strip()
            if sequence_id:
                rows[sequence_id] = row
    return rows


def read_beast2_sequence_coordinates(path: Path) -> Dict[str, Dict[str, object]]:
    coordinates: Dict[str, Dict[str, object]] = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sequence_id = (row.get("sequence_id") or "").strip()
            lat = beast2_parse_float(row.get("latitude"))
            lon = beast2_parse_float(row.get("longitude"))
            if not sequence_id or lat is None or lon is None:
                continue
            coordinates[sequence_id] = {
                "lat": lat,
                "lon": lon,
                "location_id": (row.get("location_id") or "").strip(),
                "location_label": (row.get("location_label") or "").strip() or sequence_id,
            }
    return coordinates


def parse_beast2_mcc_taxlabels(path: Path) -> Dict[str, str]:
    labels: List[str] = []
    in_taxlabels = False
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            lower = line.lower()
            if lower.startswith("taxlabels"):
                in_taxlabels = True
                line = line[len("taxlabels"):].strip()
                if not line:
                    continue
            if not in_taxlabels:
                continue
            finished = line.endswith(";")
            if finished:
                line = line[:-1].strip()
            if line:
                labels.extend(part.strip("'\"") for part in line.split() if part.strip())
            if finished:
                break
    return {str(index + 1): label for index, label in enumerate(labels)}


def extract_beast2_mcc_newick(path: Path) -> str:
    collecting = False
    chunks: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not collecting and line.lower().startswith("tree "):
                collecting = True
                line = line.split("=", 1)[1].strip()
            if collecting:
                chunks.append(line)
                if line.endswith(";"):
                    break
    if not chunks:
        raise RuntimeError(f"No MCC tree line found in {path}")
    return " ".join(chunks)


def replace_beast2_numeric_tip_labels(newick: str, taxlabels: Dict[str, str]) -> str:
    out: List[str] = []
    i = 0
    n = len(newick)
    while i < n:
        ch = newick[i]
        out.append(ch)
        i += 1
        if ch not in "(,":
            continue
        while i < n and newick[i].isspace():
            out.append(newick[i])
            i += 1
        j = i
        while j < n and newick[j].isdigit():
            j += 1
        token = newick[i:j]
        if token in taxlabels and (j == n or newick[j] in "[:),"):
            out.append(taxlabels[token])
            i = j
    return "".join(out)


def read_beast2_mcc_tree(path: Path):
    from io import StringIO
    from Bio import Phylo

    taxlabels = parse_beast2_mcc_taxlabels(path)
    newick = extract_beast2_mcc_newick(path)
    if taxlabels:
        newick = replace_beast2_numeric_tip_labels(newick, taxlabels)
    return Phylo.read(StringIO(newick), "newick")


def parse_beast2_comment_number(comment: str | None, key: str) -> float | None:
    if not comment:
        return None
    match = re.search(r"(?:^|[,&])" + re.escape(key) + r"=([-+0-9.Ee]+)", comment)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def is_local_beast2_tip(sequence_id: str, metadata: Dict[str, Dict[str, str]]) -> bool:
    row = metadata.get(sequence_id, {})
    if (row.get("record_type") or "").strip().lower() == "local":
        return True
    return sequence_id.startswith(("IVIA_", "LNR_"))


def annotate_beast2_tree_for_visuals(tree, tip_dates: Dict[str, float],
                                     tip_coordinates: Dict[str, Dict[str, object]]) -> float:
    if not tip_dates:
        raise RuntimeError("Cannot render BEAST 2 visuals: no tip dates found")
    newest_date = max(tip_dates.values())
    node_counter = 0

    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            clade._beast_node_id = clade.name or "tip"
        else:
            node_counter += 1
            clade._beast_node_id = f"node_{node_counter:04d}"

    def visit(clade):
        comment = getattr(clade, "comment", None)
        parsed_height = parse_beast2_comment_number(comment, "height")
        if clade.is_terminal():
            tip_date = tip_dates.get(clade.name or "")
            fallback_height = newest_date - tip_date if tip_date is not None else 0.0
            height = parsed_height if parsed_height is not None else fallback_height
            coord_row = tip_coordinates.get(clade.name or "")
            coord = None
            if coord_row:
                coord = (float(coord_row["lon"]), float(coord_row["lat"]))
            clade._beast_tip_count = 1
            clade._beast_coord = coord
            clade._beast_height = height
            clade._beast_date = newest_date - height
            clade._beast_posterior = parse_beast2_comment_number(comment, "posterior")
            return coord, height, 1

        child_results = [visit(child) for child in clade.clades]
        weighted_lon = 0.0
        weighted_lat = 0.0
        weight = 0
        tip_count = 0
        for coord, _height, child_tips in child_results:
            tip_count += child_tips
            if coord is None:
                continue
            weighted_lon += coord[0] * child_tips
            weighted_lat += coord[1] * child_tips
            weight += child_tips
        coord = (weighted_lon / weight, weighted_lat / weight) if weight else None
        if parsed_height is not None:
            height = parsed_height
        else:
            child_heights = [h + (child.branch_length or 0.0) for child, (_c, h, _n) in zip(clade.clades, child_results)]
            height = max(child_heights) if child_heights else 0.0
        clade._beast_tip_count = max(1, tip_count)
        clade._beast_coord = coord
        clade._beast_height = height
        clade._beast_date = newest_date - height
        clade._beast_posterior = parse_beast2_comment_number(comment, "posterior")
        return coord, height, clade._beast_tip_count

    visit(tree.root)
    return newest_date


def write_beast2_phylogeography_edges(tree, path: Path) -> None:
    rows: List[Dict[str, object]] = []
    for parent in tree.find_clades(order="preorder"):
        parent_coord = getattr(parent, "_beast_coord", None)
        if parent_coord is None:
            continue
        for child in parent.clades:
            child_coord = getattr(child, "_beast_coord", None)
            if child_coord is None:
                continue
            rows.append({
                "parent_node": getattr(parent, "_beast_node_id", ""),
                "child_node": getattr(child, "_beast_node_id", ""),
                "child_is_tip": "YES" if child.is_terminal() else "NO",
                "parent_decimal_date": f"{getattr(parent, '_beast_date', 0.0):.6f}",
                "child_decimal_date": f"{getattr(child, '_beast_date', 0.0):.6f}",
                "parent_longitude": f"{parent_coord[0]:.6f}",
                "parent_latitude": f"{parent_coord[1]:.6f}",
                "child_longitude": f"{child_coord[0]:.6f}",
                "child_latitude": f"{child_coord[1]:.6f}",
                "child_descendant_tips": getattr(child, "_beast_tip_count", 1),
                "child_posterior": "" if getattr(child, "_beast_posterior", None) is None else f"{child._beast_posterior:.6f}",
            })
    write_tsv(
        path,
        rows,
        [
            "parent_node", "child_node", "child_is_tip", "parent_decimal_date", "child_decimal_date",
            "parent_longitude", "parent_latitude", "child_longitude", "child_latitude",
            "child_descendant_tips", "child_posterior",
        ],
    )


def draw_simple_world_context(ax) -> None:
    from matplotlib.patches import Polygon

    ax.set_facecolor("#eef6fb")
    land_polygons = [
        [(-168, 14), (-150, 58), (-108, 72), (-55, 52), (-50, 25), (-82, 8), (-118, 16)],
        [(-82, 12), (-74, -18), (-65, -55), (-38, -52), (-34, -10), (-52, 8)],
        [(-12, 36), (8, 58), (35, 70), (70, 60), (120, 55), (150, 45), (135, 8), (75, 5), (45, 24), (20, 34)],
        [(-18, 35), (50, 32), (48, -34), (18, -36), (-12, 5)],
        [(108, -10), (154, -12), (151, -42), (115, -36)],
        [(-10, 50), (5, 61), (31, 61), (42, 45), (18, 36), (0, 42)],
    ]
    for poly in land_polygons:
        ax.add_patch(Polygon(poly, closed=True, facecolor="#e7e2d1", edgecolor="#b8b0a0", linewidth=0.6, alpha=0.65, zorder=0))
    ax.grid(True, color="#b8ced8", linewidth=0.7, alpha=0.7)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")


def render_beast2_phylogeography_map(tree, tip_dates: Dict[str, float],
                                     tip_coordinates: Dict[str, Dict[str, object]],
                                     metadata: Dict[str, Dict[str, str]],
                                     pdf_path: Path, png_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    from matplotlib.ticker import MaxNLocator

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError as exc:
        raise RuntimeError(
            "Cannot render the BEAST 2 real-world map because cartopy is not installed. "
            "Install it with: conda install -c conda-forge cartopy"
        ) from exc

    dated_tips = [name for name in tip_coordinates if name in tip_dates]
    if not dated_tips:
        raise RuntimeError("Cannot render BEAST 2 map: no dated tips with coordinates")

    all_lons = [float(tip_coordinates[name]["lon"]) for name in dated_tips]
    all_lats = [float(tip_coordinates[name]["lat"]) for name in dated_tips]
    lon_span = max(all_lons) - min(all_lons)
    lat_span = max(all_lats) - min(all_lats)

    lon_margin = max(8.0, lon_span * 0.08)
    lon_min = max(-180.0, min(all_lons) - lon_margin)
    lon_max = min(180.0, max(all_lons) + lon_margin)

    if lon_span > 140.0:
        lat_center = (min(all_lats) + max(all_lats)) / 2.0
        min_visible_lat_span = 64.0
        lat_min = max(-20.0, min(all_lats) - 10.0, lat_center - min_visible_lat_span / 2.0)
        lat_max = min(80.0, max(all_lats) + 10.0, lat_center + min_visible_lat_span / 2.0)
        if lat_max - lat_min < min_visible_lat_span:
            extra = (min_visible_lat_span - (lat_max - lat_min)) / 2.0
            lat_min = max(-20.0, lat_min - extra)
            lat_max = min(80.0, lat_max + extra)
    else:
        lat_margin = max(5.0, lat_span * 0.22)
        lat_min = max(-70.0, min(all_lats) - lat_margin)
        lat_max = min(85.0, max(all_lats) + lat_margin)

    date_values = [tip_dates[name] for name in dated_tips]
    norm = Normalize(vmin=min(date_values), vmax=max(date_values))
    cmap = plt.get_cmap("viridis")
    data_crs = ccrs.PlateCarree()

    fig = plt.figure(figsize=(18.5, 8.6))
    ax = fig.add_axes([0.045, 0.19, 0.91, 0.73], projection=data_crs)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=data_crs)

    try:
        ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="#eaf4fb", zorder=0)
        ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="#f3efe3", edgecolor="none", zorder=0)
        ax.add_feature(cfeature.COASTLINE.with_scale("50m"), edgecolor="#5f6f73", linewidth=0.6, zorder=1)
        ax.add_feature(cfeature.BORDERS.with_scale("50m"), edgecolor="#8a8177", linewidth=0.45, zorder=1)
        ax.add_feature(cfeature.LAKES.with_scale("50m"), facecolor="#dcecf4", edgecolor="#b6cbd4", linewidth=0.28, zorder=1)
        admin1 = cfeature.NaturalEarthFeature(
            "cultural",
            "admin_1_states_provinces_lines",
            "50m",
            facecolor="none",
        )
        ax.add_feature(admin1, edgecolor="#b6aea0", linewidth=0.28, alpha=0.72, zorder=1)
    except Exception as exc:
        raise RuntimeError(
            "Cartopy is installed, but Natural Earth map data could not be loaded. "
            "Run once with internet access so Cartopy can cache the map layers."
        ) from exc

    gl = ax.gridlines(draw_labels=True, linewidth=0.45, color="#a8bfca", alpha=0.75, linestyle="-")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 10}
    gl.ylabel_style = {"size": 10}

    segments = []
    segment_dates = []
    for parent in tree.find_clades(order="preorder"):
        parent_coord = getattr(parent, "_beast_coord", None)
        if parent_coord is None:
            continue
        for child in parent.clades:
            child_coord = getattr(child, "_beast_coord", None)
            if child_coord is None:
                continue
            segments.append([parent_coord, child_coord])
            segment_dates.append(getattr(child, "_beast_date", getattr(parent, "_beast_date", min(date_values))))
    if segments:
        lc = LineCollection(
            segments,
            colors=[cmap(norm(value)) for value in segment_dates],
            linewidths=1.05,
            alpha=0.34,
            transform=data_crs,
            zorder=2,
        )
        ax.add_collection(lc)

    location_groups: Dict[str, Dict[str, object]] = {}
    for sequence_id in dated_tips:
        row = tip_coordinates[sequence_id]
        label = str(row.get("location_label") or sequence_id)
        group = location_groups.setdefault(label, {
            "lon": float(row["lon"]),
            "lat": float(row["lat"]),
            "count": 0,
            "dates": [],
            "local_count": 0,
        })
        group["count"] = int(group["count"]) + 1
        group["dates"].append(tip_dates[sequence_id])
        if is_local_beast2_tip(sequence_id, metadata):
            group["local_count"] = int(group["local_count"]) + 1

    for label, group in location_groups.items():
        count = int(group["count"])
        local_count = int(group["local_count"])
        mean_date = sum(group["dates"]) / len(group["dates"])
        size = 42 + min(430, count * 18)
        edge = "#c0392b" if local_count else "#1f2d3a"
        linewidth = 2.2 if local_count else 0.85
        ax.scatter(
            [group["lon"]],
            [group["lat"]],
            s=size,
            color=[cmap(norm(mean_date))],
            edgecolor=edge,
            linewidth=linewidth,
            alpha=0.96,
            transform=data_crs,
            zorder=4,
        )
        if local_count or count >= 3:
            short_label = label.split("|")[-1].strip() or label.split("|")[0].strip()
            short_label = short_label.replace("__", "_")
            if len(short_label) > 27:
                short_label = short_label[:24] + "..."
            ax.text(
                float(group["lon"]) + 0.75,
                float(group["lat"]) + 0.55,
                f"{short_label} ({count})",
                fontsize=8.8,
                color="#1f2d3a",
                transform=data_crs,
                zorder=5,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.58, "pad": 1.4},
            )

    ax.set_title("BEAST 2 MCC phylogeography map", fontsize=18, pad=14)
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    cax = fig.add_axes([0.22, 0.072, 0.56, 0.038])
    cbar = fig.colorbar(scalar, cax=cax, orientation="horizontal")
    cbar.locator = MaxNLocator(nbins=7)
    cbar.update_ticks()
    cbar.ax.tick_params(labelsize=11, length=4, pad=4)
    cbar.set_label("Sampling / inferred branch time (decimal year)", fontsize=12, labelpad=8)
    ax.text(
        0.01,
        0.018,
        "Map base: Natural Earth via Cartopy. Internal node coordinates are descendant-tip centroids for visualization.",
        transform=ax.transAxes,
        fontsize=9,
        color="#34495e",
        ha="left",
        va="bottom",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 3.5},
    )
    fig.savefig(str(pdf_path), format="pdf", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(str(png_path), format="png", dpi=320, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)



def render_beast2_interactive_map(tree, tip_dates: Dict[str, float],
                                  tip_coordinates: Dict[str, Dict[str, object]],
                                  metadata: Dict[str, Dict[str, str]],
                                  html_path: Path) -> None:
    import json
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    dated_tips = [name for name in tip_coordinates if name in tip_dates]
    if not dated_tips:
        raise RuntimeError("Cannot render BEAST 2 interactive map: no dated tips with coordinates")

    date_values = [tip_dates[name] for name in dated_tips]
    date_min = min(date_values)
    date_max = max(date_values)
    norm = Normalize(vmin=date_min, vmax=date_max)
    cmap = plt.get_cmap("viridis")

    def color_for(value: float) -> str:
        rgba = cmap(norm(value))
        return "#{:02x}{:02x}{:02x}".format(
            int(rgba[0] * 255),
            int(rgba[1] * 255),
            int(rgba[2] * 255),
        )

    location_groups: Dict[str, Dict[str, object]] = {}
    for sequence_id in dated_tips:
        row = tip_coordinates[sequence_id]
        label = str(row.get("location_label") or sequence_id)
        group = location_groups.setdefault(label, {
            "label": label,
            "lon": float(row["lon"]),
            "lat": float(row["lat"]),
            "count": 0,
            "dates": [],
            "local_count": 0,
            "reference_count": 0,
            "sequences": [],
        })
        group["count"] = int(group["count"]) + 1
        group["dates"].append(tip_dates[sequence_id])
        row_meta = metadata.get(sequence_id, {})
        is_local = is_local_beast2_tip(sequence_id, metadata)
        if is_local:
            group["local_count"] = int(group["local_count"]) + 1
        else:
            group["reference_count"] = int(group["reference_count"]) + 1
        group["sequences"].append({
            "id": sequence_id,
            "date": round(float(tip_dates[sequence_id]), 6),
            "record_type": "local" if is_local else "reference",
            "accession": row_meta.get("accession", ""),
            "isolate": row_meta.get("isolate", ""),
            "country": row_meta.get("country", ""),
            "region": row_meta.get("region", ""),
            "locality": row_meta.get("locality", ""),
            "sample_type": row_meta.get("sample_type", ""),
            "collection_date": row_meta.get("collection_date", ""),
        })

    locations = []
    for label, group in sorted(location_groups.items()):
        dates = [float(value) for value in group["dates"]]
        mean_date = sum(dates) / len(dates)
        sequences = sorted(group["sequences"], key=lambda item: (item["record_type"], item["id"]))
        locations.append({
            "label": label,
            "lat": round(float(group["lat"]), 8),
            "lon": round(float(group["lon"]), 8),
            "count": int(group["count"]),
            "local_count": int(group["local_count"]),
            "reference_count": int(group["reference_count"]),
            "mean_date": round(mean_date, 6),
            "min_date": round(min(dates), 6),
            "max_date": round(max(dates), 6),
            "color": color_for(mean_date),
            "sequences": sequences,
        })

    edges = []
    for parent in tree.find_clades(order="preorder"):
        parent_coord = getattr(parent, "_beast_coord", None)
        if parent_coord is None:
            continue
        for child in parent.clades:
            child_coord = getattr(child, "_beast_coord", None)
            if child_coord is None:
                continue
            child_date = float(getattr(child, "_beast_date", getattr(parent, "_beast_date", date_min)))
            edges.append({
                "parent_lat": round(float(parent_coord[1]), 8),
                "parent_lon": round(float(parent_coord[0]), 8),
                "child_lat": round(float(child_coord[1]), 8),
                "child_lon": round(float(child_coord[0]), 8),
                "date": round(child_date, 6),
                "color": color_for(child_date),
                "child_is_tip": bool(child.is_terminal()),
                "posterior": None if getattr(child, "_beast_posterior", None) is None else round(float(child._beast_posterior), 6),
            })

    center_lat = sum(item["lat"] for item in locations) / len(locations)
    center_lon = sum(item["lon"] for item in locations) / len(locations)
    max_count = max(item["count"] for item in locations)

    payload = {
        "locations": locations,
        "edges": edges,
        "center": [round(center_lat, 8), round(center_lon, 8)],
        "date_min": round(date_min, 6),
        "date_max": round(date_max, 6),
        "max_count": int(max_count),
    }
    payload_json = json.dumps(payload, ensure_ascii=False)

    html_template = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>BEAST 2 MCC interactive phylogeography map</title>
  <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
  <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
  <style>
    html, body, #map {
      height: 100%;
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }
    .panel, .legend {
      position: absolute;
      z-index: 1000;
      background: rgba(255, 255, 255, 0.94);
      border: 1px solid rgba(30, 45, 55, 0.18);
      border-radius: 8px;
      box-shadow: 0 8px 22px rgba(0, 0, 0, 0.14);
      color: #263238;
    }
    .panel {
      left: 14px;
      top: 14px;
      width: min(405px, calc(100vw - 28px));
    }
    .panel.collapsed {
      width: auto;
    }
    .panel-header, .legend-header {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 10px;
      padding: 8px 10px;
      font-weight: 700;
      font-size: 13px;
      color: #17202a;
    }
    .panel-body {
      padding: 0 10px 10px;
    }
    .panel.collapsed .panel-body,
    .legend.collapsed .legend-body {
      display: none;
    }
    .panel p, .legend p {
      font-size: 12px;
      line-height: 1.35;
      margin: 7px 0 0;
      color: #34495e;
    }
    .mini-button {
      border: 0;
      border-radius: 6px;
      padding: 5px 8px;
      background: #1f6f8b;
      color: white;
      font-size: 12px;
      cursor: pointer;
      white-space: nowrap;
    }
    .search-row {
      display: flex;
      gap: 6px;
      margin-top: 9px;
    }
    .search-row input {
      flex: 1;
      border: 1px solid #b7c3ca;
      border-radius: 6px;
      padding: 7px 8px;
      font-size: 13px;
      min-width: 0;
    }
    .legend {
      right: 14px;
      bottom: 18px;
      width: min(360px, calc(100vw - 28px));
      font-size: 12px;
    }
    .legend.collapsed {
      width: auto;
    }
    .legend-body {
      padding: 0 12px 12px;
    }
    .gradient {
      height: 13px;
      border-radius: 8px;
      background: linear-gradient(to right, #440154, #3b528b, #21918c, #5ec962, #fde725);
      margin: 7px 0 4px;
    }
    .legend-scale {
      display: flex;
      justify-content: space-between;
      font-variant-numeric: tabular-nums;
    }
    .popup-title {
      font-weight: 700;
      margin-bottom: 4px;
      color: #17202a;
    }
    .popup-meta {
      font-size: 12px;
      line-height: 1.35;
      color: #263238;
      margin-bottom: 6px;
    }
    .popup-list {
      max-height: 190px;
      overflow-y: auto;
      border-top: 1px solid #d9e1e5;
      padding-top: 5px;
      font-size: 11.5px;
      line-height: 1.35;
    }
    .local {
      color: #a93226;
      font-weight: 700;
    }
    .leaflet-control-layers {
      font-size: 12px;
    }
  </style>
</head>
<body>
  <div id="map"></div>
  <div class="panel collapsed" id="infoPanel">
    <div class="panel-header">
      <span>BEAST 2 MCC map</span>
      <button class="mini-button" id="panelToggle">Tools</button>
    </div>
    <div class="panel-body">
      <p>Zoom freely into any region. Red outlines mark local assemblies; fill color follows sampling or inferred branch time.</p>
      <div class="search-row">
        <input id="search" type="search" placeholder="Search isolate, accession, country, locality...">
        <button class="mini-button" id="searchButton">Zoom</button>
        <button class="mini-button" id="resetButton" title="Reset map view">Reset</button>
      </div>
      <p id="searchStatus"></p>
    </div>
  </div>
  <div class="legend collapsed" id="legendPanel">
    <div class="legend-header">
      <span>Time scale</span>
      <button class="mini-button" id="legendToggle">Show</button>
    </div>
    <div class="legend-body">
      <strong>Sampling / inferred branch time</strong>
      <div class="gradient"></div>
      <div class="legend-scale">
        <span>__DATE_MIN__</span>
        <span>__DATE_MAX__</span>
      </div>
      <p>Branches are drawn with a white halo plus a colored MCC line. Circle size is proportional to the number of sequences at a location.</p>
    </div>
  </div>
  <script>
    const DATA = __PAYLOAD_JSON__;

    function escapeHtml(value) {
      return String(value ?? "")
        .replaceAll("&", "&amp;")
        .replaceAll("<", "&lt;")
        .replaceAll(">", "&gt;")
        .replaceAll('"', "&quot;")
        .replaceAll("'", "&#039;");
    }

    function popupContent(loc) {
      const rows = loc.sequences.map(seq => {
        const cls = seq.record_type === "local" ? "local" : "";
        const label = seq.collection_date ? `${seq.id} · ${seq.collection_date}` : `${seq.id} · ${seq.date}`;
        const details = [seq.country, seq.region, seq.locality, seq.sample_type].filter(Boolean).join(" / ");
        return `<div class="${cls}">${escapeHtml(label)}${details ? `<br><span>${escapeHtml(details)}</span>` : ""}</div>`;
      }).join("");
      return `
        <div class="popup-title">${escapeHtml(loc.label)}</div>
        <div class="popup-meta">
          <strong>${loc.count}</strong> sequence(s): ${loc.local_count} local, ${loc.reference_count} reference<br>
          Mean date: ${loc.mean_date.toFixed(3)}<br>
          Date range: ${loc.min_date.toFixed(3)} - ${loc.max_date.toFixed(3)}<br>
          Coordinates: ${loc.lat.toFixed(5)}, ${loc.lon.toFixed(5)}
        </div>
        <div class="popup-list">${rows}</div>
      `;
    }

    const map = L.map("map", { preferCanvas: true }).setView(DATA.center, 3);
    map.createPane("branches");
    map.getPane("branches").style.zIndex = 430;
    map.createPane("samples");
    map.getPane("samples").style.zIndex = 460;

    const cleanBase = L.tileLayer("https://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}{r}.png", {
      maxZoom: 20,
      attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>'
    }).addTo(map);
    const labelOverlay = L.tileLayer("https://{s}.basemaps.cartocdn.com/light_only_labels/{z}/{x}/{y}{r}.png", {
      maxZoom: 20,
      attribution: '&copy; <a href="https://carto.com/attributions">CARTO</a>'
    });
    const osm = L.tileLayer("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png", {
      maxZoom: 19,
      attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
    });

    const mccBranchLayer = L.layerGroup().addTo(map);
    const samplingLayer = L.layerGroup().addTo(map);
    const localLayer = L.layerGroup();
    const referenceLayer = L.layerGroup();

    DATA.edges.forEach(edge => {
      const points = [[edge.parent_lat, edge.parent_lon], [edge.child_lat, edge.child_lon]];
      L.polyline(points, {
        pane: "branches",
        color: "#ffffff",
        weight: edge.child_is_tip ? 6.2 : 5.2,
        opacity: 0.78,
        smoothFactor: 0
      }).addTo(mccBranchLayer);
      L.polyline(points, {
        pane: "branches",
        color: edge.color,
        weight: edge.child_is_tip ? 3.15 : 2.45,
        opacity: 0.98,
        smoothFactor: 0
      }).addTo(mccBranchLayer);
    });

    const searchable = [];
    const bounds = [];
    DATA.locations.forEach(loc => {
      const radius = Math.max(6, Math.min(28, 5 + Math.sqrt(loc.count) * 4.0));
      const style = {
        pane: "samples",
        radius,
        fillColor: loc.color,
        color: loc.local_count > 0 ? "#c0392b" : "#1f2d3a",
        weight: loc.local_count > 0 ? 3.2 : 1.1,
        opacity: 0.98,
        fillOpacity: 0.86
      };
      const marker = L.circleMarker([loc.lat, loc.lon], style).bindPopup(popupContent(loc), { maxWidth: 410 });
      marker.addTo(samplingLayer);
      if (loc.local_count > 0) marker.addTo(localLayer);
      if (loc.reference_count > 0) marker.addTo(referenceLayer);
      const searchText = [
        loc.label,
        ...loc.sequences.flatMap(seq => [seq.id, seq.accession, seq.isolate, seq.country, seq.region, seq.locality, seq.sample_type])
      ].join(" ").toLowerCase();
      searchable.push({ loc, marker, searchText });
      bounds.push([loc.lat, loc.lon]);
    });

    if (bounds.length > 0) {
      map.fitBounds(bounds, { padding: [34, 34], maxZoom: 5 });
    }

    L.control.layers(
      {
        "Clean map": cleanBase,
        "OpenStreetMap": osm
      },
      {
        "Map labels": labelOverlay,
        "MCC branches": mccBranchLayer,
        "Sampling locations": samplingLayer,
        "Local assemblies": localLayer,
        "Reference genomes": referenceLayer
      },
      { collapsed: true }
    ).addTo(map);

    const originalBounds = bounds.length > 0 ? L.latLngBounds(bounds) : null;
    const searchInput = document.getElementById("search");
    const searchStatus = document.getElementById("searchStatus");
    const infoPanel = document.getElementById("infoPanel");
    const legendPanel = document.getElementById("legendPanel");
    const panelToggle = document.getElementById("panelToggle");
    const legendToggle = document.getElementById("legendToggle");

    panelToggle.addEventListener("click", () => {
      infoPanel.classList.toggle("collapsed");
      panelToggle.textContent = infoPanel.classList.contains("collapsed") ? "Tools" : "Hide";
    });
    legendToggle.addEventListener("click", () => {
      legendPanel.classList.toggle("collapsed");
      legendToggle.textContent = legendPanel.classList.contains("collapsed") ? "Show" : "Hide";
    });

    function runSearch() {
      const query = searchInput.value.trim().toLowerCase();
      if (!query) {
        searchStatus.textContent = "";
        return;
      }
      const matches = searchable.filter(item => item.searchText.includes(query));
      if (matches.length === 0) {
        searchStatus.textContent = "No matching isolate or location found.";
        return;
      }
      const first = matches[0];
      const zoom = first.loc.local_count > 0 ? 10 : 6;
      map.flyTo([first.loc.lat, first.loc.lon], zoom, { duration: 0.8 });
      first.marker.openPopup();
      searchStatus.textContent = `${matches.length} match(es). Showing: ${first.loc.label}`;
    }

    document.getElementById("searchButton").addEventListener("click", runSearch);
    searchInput.addEventListener("keydown", event => {
      if (event.key === "Enter") runSearch();
    });
    document.getElementById("resetButton").addEventListener("click", () => {
      if (originalBounds) map.fitBounds(originalBounds, { padding: [34, 34], maxZoom: 5 });
      searchInput.value = "";
      searchStatus.textContent = "";
    });
  </script>
</body>
</html>
"""
    html = (
        html_template
        .replace("__PAYLOAD_JSON__", payload_json)
        .replace("__DATE_MIN__", f"{date_min:.2f}")
        .replace("__DATE_MAX__", f"{date_max:.2f}")
    )
    html_path.write_text(html, encoding="utf-8")
    plt.close("all")



def assign_beast2_time_tree_y_positions(tree) -> int:
    counter = {"y": 0}

    def visit(clade):
        if clade.is_terminal():
            clade._plot_y = counter["y"]
            counter["y"] += 1
            return clade._plot_y
        child_y = [visit(child) for child in clade.clades]
        clade._plot_y = sum(child_y) / len(child_y) if child_y else counter["y"]
        return clade._plot_y

    visit(tree.root)
    return counter["y"]


def render_beast2_time_tree(tree, metadata: Dict[str, Dict[str, str]],
                            pdf_path: Path, png_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    tip_count = assign_beast2_time_tree_y_positions(tree)
    terminals = tree.get_terminals()
    tip_dates = [getattr(tip, "_beast_date", None) for tip in terminals]
    tip_dates = [float(value) for value in tip_dates if value is not None]
    if not tip_dates:
        raise RuntimeError("Cannot render BEAST 2 time tree: no terminal sampling dates found")

    min_tip_date = min(tip_dates)
    max_tip_date = max(tip_dates)
    date_span = max(1.0, max_tip_date - min_tip_date)
    margin = max(0.5, date_span * 0.035)

    def display_date(value: float | None) -> float:
        if value is None:
            return min_tip_date
        return max(float(value), min_tip_date)

    fig_width = 16
    fig_height = max(10, min(48, tip_count * 0.34))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    for parent in tree.find_clades(order="preorder"):
        if not parent.clades:
            continue
        parent_x = display_date(getattr(parent, "_beast_date", min_tip_date))
        child_ys = [getattr(child, "_plot_y", 0.0) for child in parent.clades]
        ax.plot([parent_x, parent_x], [min(child_ys), max(child_ys)], color="#5d6d7e", linewidth=0.8, alpha=0.75)
        for child in parent.clades:
            child_x = display_date(getattr(child, "_beast_date", parent_x))
            child_y = getattr(child, "_plot_y", 0.0)
            local = child.is_terminal() and is_local_beast2_tip(child.name or "", metadata)
            linewidth = 1.35 if local else 0.8
            color = "#b03a2e" if local else "#5d6d7e"
            ax.plot([parent_x, child_x], [child_y, child_y], color=color, linewidth=linewidth, alpha=0.88)

    for tip in terminals:
        x = display_date(getattr(tip, "_beast_date", max_tip_date))
        y = getattr(tip, "_plot_y", 0.0)
        local = is_local_beast2_tip(tip.name or "", metadata)
        ax.scatter([x], [y], s=18 if local else 9, color="#b03a2e" if local else "#263238", zorder=3)
        ax.text(
            max_tip_date + margin * 0.25,
            y,
            tip.name or "",
            fontsize=6.3,
            va="center",
            color="#b03a2e" if local else "#263238",
            fontweight="bold" if local else "normal",
        )

    ax.axvline(min_tip_date, color="#85929e", linewidth=0.8, linestyle="--", alpha=0.8)
    ax.set_title("BEAST 2 MCC time tree (sample-date window)", fontsize=14, pad=12)
    ax.set_xlabel("Sampling year (decimal year)")
    ax.set_yticks([])
    ax.set_ylim(-1, tip_count)
    ax.set_xlim(min_tip_date - margin, max_tip_date + margin * 2.8)
    ax.grid(axis="x", color="#d5dce2", linewidth=0.7)
    ax.text(
        min_tip_date,
        -0.8,
        f"Earliest sampled isolate: {min_tip_date:.2f}",
        fontsize=8,
        color="#5d6d7e",
        ha="left",
        va="bottom",
    )
    legend_handles = [
        Line2D([0], [0], color="#b03a2e", lw=1.6, marker="o", markersize=5, label="Local isolates"),
        Line2D([0], [0], color="#5d6d7e", lw=1.2, marker="o", markersize=4, label="Reference genomes"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(str(pdf_path), format="pdf", bbox_inches="tight")
    fig.savefig(str(png_path), format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def render_beast2_visual_outputs(shared_layout: Dict[str, Path]) -> None:
    import os

    prep_outputs = beast2_outputs(shared_layout)
    run_outputs = beast2_run_outputs(shared_layout)
    mpl_config = run_outputs["outdir"] / ".matplotlib"
    mkdir(mpl_config)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))

    tree = read_beast2_mcc_tree(run_outputs["mcc_tree"])
    tip_dates = read_beast2_tip_dates_table(prep_outputs["tip_dates"])
    tip_coordinates = read_beast2_sequence_coordinates(prep_outputs["sequence_coordinates"])
    metadata = read_beast2_metadata_table(prep_outputs["metadata"])

    annotate_beast2_tree_for_visuals(tree, tip_dates, tip_coordinates)
    write_beast2_phylogeography_edges(tree, run_outputs["map_edges_tsv"])
    render_beast2_phylogeography_map(
        tree,
        tip_dates,
        tip_coordinates,
        metadata,
        run_outputs["map_pdf"],
        run_outputs["map_png"],
    )
    render_beast2_interactive_map(
        tree,
        tip_dates,
        tip_coordinates,
        metadata,
        run_outputs["map_html"],
    )
    render_beast2_time_tree(tree, metadata, run_outputs["time_tree_pdf"], run_outputs["time_tree_png"])
    run_outputs["visuals_version"].write_text("beast2_visuals_v6\n", encoding="utf-8")

    with open(run_outputs["log"], "a", encoding="utf-8") as logh:
        logh.write(f"[{now()}] BEAST 2 phylogeography map PDF written to: {run_outputs['map_pdf']}\n")
        logh.write(f"[{now()}] BEAST 2 phylogeography map PNG written to: {run_outputs['map_png']}\n")
        logh.write(f"[{now()}] BEAST 2 interactive phylogeography HTML written to: {run_outputs['map_html']}\n")
        logh.write(f"[{now()}] BEAST 2 time tree PDF written to: {run_outputs['time_tree_pdf']}\n")
        logh.write(f"[{now()}] BEAST 2 time tree PNG written to: {run_outputs['time_tree_png']}\n")
        logh.write(f"[{now()}] BEAST 2 phylogeography edge table written to: {run_outputs['map_edges_tsv']}\n")


# ---------------------------------------------------------------------
# Step 11 - Pairwise identity heatmap
# ---------------------------------------------------------------------

def step11_identity_outputs(shared_layout: Dict[str, Path]) -> Tuple[Path, Path, Path, Path, Path, Path]:
    matrix_tsv = shared_layout["identity"] / "pairwise_identity.tsv"
    pdf = shared_layout["identity"] / "pairwise_identity_heatmap.pdf"
    png = shared_layout["identity"] / "pairwise_identity_heatmap.png"
    assembled_matrix_tsv = shared_layout["identity"] / "pairwise_identity_assembled_only.tsv"
    assembled_pdf = shared_layout["identity"] / "pairwise_identity_assembled_only_heatmap.pdf"
    assembled_png = shared_layout["identity"] / "pairwise_identity_assembled_only_heatmap.png"
    return matrix_tsv, pdf, png, assembled_matrix_tsv, assembled_pdf, assembled_png


def step11_identity_done(shared_layout: Dict[str, Path]) -> bool:
    matrix_tsv, pdf, png, assembled_matrix_tsv, assembled_pdf, assembled_png = step11_identity_outputs(shared_layout)
    return (
        matrix_tsv.exists()
        and pdf.exists()
        and png.exists()
        and assembled_matrix_tsv.exists()
        and assembled_pdf.exists()
        and assembled_png.exists()
    )


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


def subset_matrix(labels: List[str], matrix: List[List[float]], keep_labels: List[str]) -> Tuple[List[str], List[List[float]]]:
    index = {label: i for i, label in enumerate(labels)}
    keep = [label for label in keep_labels if label in index]
    subset = []
    for row_label in keep:
        i = index[row_label]
        subset.append([matrix[i][index[col_label]] for col_label in keep])
    return keep, subset


def write_identity_matrix(matrix_tsv: Path, labels: List[str], matrix: List[List[float]]) -> None:
    with open(matrix_tsv, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sequence"] + labels)
        for label, row in zip(labels, matrix):
            writer.writerow([label] + [f"{value:.4f}" for value in row])


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
    matrix_tsv, pdf, png, assembled_matrix_tsv, assembled_pdf, assembled_png = step11_identity_outputs(shared_layout)
    treefile = step12_tree_outputs(shared_layout)[0]
    records = load_alignment_records(trimmed_alignment)
    if not records:
        raise RuntimeError("Trimmed alignment is empty")
    labels, matrix = compute_identity_matrix(records)
    if treefile.exists():
        tree = Phylo.read(str(treefile), "newick")
        tree_order = [term.name for term in tree.get_terminals()]
        labels, matrix = reorder_matrix(labels, matrix, tree_order)

    write_identity_matrix(matrix_tsv, labels, matrix)
    render_identity_heatmap(labels, matrix, pdf, png, "CYVCV genome identity heatmap", plot_min)

    assembled_labels = [label for label in labels if is_assembled_tip(label)]
    if not assembled_labels:
        raise RuntimeError("No assembled isolates were found for the assembled-only identity heatmap")
    assembled_labels, assembled_matrix = subset_matrix(labels, matrix, assembled_labels)
    write_identity_matrix(assembled_matrix_tsv, assembled_labels, assembled_matrix)
    render_identity_heatmap(
        assembled_labels,
        assembled_matrix,
        assembled_pdf,
        assembled_png,
        "Assembled isolate identity heatmap",
        plot_min,
    )

    with open(log_file, "a", encoding="utf-8") as logh:
        logh.write(f"[{now()}] Pairwise identity matrix written to: {matrix_tsv}\n")
        logh.write(f"[{now()}] Heatmap PDF written to: {pdf}\n")
        logh.write(f"[{now()}] Heatmap PNG written to: {png}\n")
        logh.write(f"[{now()}] Assembled-only identity matrix written to: {assembled_matrix_tsv}\n")
        logh.write(f"[{now()}] Assembled-only heatmap PDF written to: {assembled_pdf}\n")
        logh.write(f"[{now()}] Assembled-only heatmap PNG written to: {assembled_png}\n")


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


def is_bootstrap_label(label: str) -> bool:
    clean = label.strip()
    if not clean:
        return False
    if "/" in clean:
        return all(is_bootstrap_label(part) for part in clean.split("/") if part.strip())
    try:
        value = float(clean)
    except ValueError:
        return False
    return 0.0 <= value <= 100.0


def render_tree_pdf(treefile: Path, pdf_path: Path) -> None:
    from Bio import Phylo
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    tree = Phylo.read(str(treefile), "newick")
    terminals = max(1, len(tree.get_terminals()))
    fig_width = max(14, min(42, terminals * 0.45))
    fig_height = max(9, min(52, terminals * 0.55))
    tip_font_size = max(7.5, min(10.5, 13.0 - terminals * 0.055))
    sample_font_size = min(12.0, tip_font_size + 1.3)
    bootstrap_font_size = max(8.0, min(11.0, tip_font_size + 1.0))

    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False)

    for line in ax.get_lines():
        line.set_color("#4d5656")
        line.set_linewidth(0.9)
        line.set_alpha(0.92)
    for collection in ax.collections:
        try:
            collection.set_color("#4d5656")
            collection.set_linewidth(0.9)
            collection.set_alpha(0.92)
        except Exception:
            pass

    for text in ax.texts:
        label = text.get_text().strip()
        if not label:
            continue
        if is_assembled_tip(label):
            text.set_color("#c0392b")
            text.set_fontweight("bold")
            text.set_fontsize(sample_font_size)
            text.set_zorder(5)
        elif is_bootstrap_label(label):
            text.set_color("#17202a")
            text.set_fontweight("bold")
            text.set_fontsize(bootstrap_font_size)
            text.set_zorder(6)
            text.set_bbox({"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.35})
        else:
            text.set_color("#283747")
            text.set_fontsize(tip_font_size)
            text.set_fontweight("normal")

    legend_handles = [
        Line2D([0], [0], color="#c0392b", lw=0, marker="o", markersize=7, label="Samples analyzed"),
        Line2D([0], [0], color="#283747", lw=0, marker="o", markersize=7, label="Reference genomes"),
        Line2D([0], [0], color="#17202a", lw=0, marker="$99$", markersize=10, label="Bootstrap support"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=max(9, tip_font_size))
    ax.set_title("Maximum-likelihood phylogeny", fontsize=14, fontweight="bold", pad=12)

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
    identity_tsv, identity_pdf, _, assembled_identity_tsv, assembled_identity_pdf, _ = step11_identity_outputs(shared_layout)
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
        "assembled_identity_matrix": str(assembled_identity_tsv),
        "assembled_identity_heatmap": str(assembled_identity_pdf),
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
        "identity_heatmap", "assembled_identity_matrix", "assembled_identity_heatmap",
        "treefile", "tree_pdf"
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
    - Global identity heatmap: {row['identity_heatmap']}
    - Assembled-only identity matrix: {row['assembled_identity_matrix']}
    - Assembled-only identity heatmap: {row['assembled_identity_heatmap']}
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
    if args.run_beast2:
        for exe in ["beast", "treeannotator"]:
            require_executable(exe)

    if not args.refseq_virus_fasta.exists():
        die(f"RefSeq virus FASTA not found: {args.refseq_virus_fasta}")
    if args.run_beast2:
        args.prepare_beast2 = True
        if args.beast2_manual_dates is None:
            die("--run-beast2 requires --beast2-manual-dates with completed sampling dates")
        if args.beast2_coordinates is None:
            die("--run-beast2 requires --beast2-coordinates with completed latitude/longitude")
    if args.beast2_manual_dates is not None and not args.beast2_manual_dates.exists():
        die(f"BEAST 2 manual dates file not found: {args.beast2_manual_dates}")
    if args.beast2_coordinates is not None and not args.beast2_coordinates.exists():
        die(f"BEAST 2 coordinates file not found: {args.beast2_coordinates}")
    if args.beast2_xml is not None and not args.beast2_xml.exists():
        die(f"BEAST 2 XML file not found: {args.beast2_xml}")
    if not (0 <= args.beast2_burnin < 100):
        die("--beast2-burnin must be between 0 and 99")


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
        "identity_heatmap", "assembled_identity_matrix", "assembled_identity_heatmap",
        "treefile", "tree_pdf"
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
            f"  Global identity heatmap: {row['identity_heatmap']}",
            f"  Assembled-only identity matrix: {row['assembled_identity_matrix']}",
            f"  Assembled-only identity heatmap: {row['assembled_identity_heatmap']}",
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
             lambda: step6_assemble(
                 layout, shared_layout, args.threads, args.flye_mode, iterations=args.flye_iterations
             )),
        Step("BLAST contigs against RefSeq Virus",
             lambda: step7_done(layout),
             lambda: step7_blast(layout, shared_layout, args.threads)),
        Step("Select target contigs",
             lambda: step8_done(layout),
             lambda: step8_select_target_contigs(
                 layout, shared_layout, args.min_pident, args.min_qcov, sample_name
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

    rescue_performed = False
    for idx, step in enumerate(steps, start=1):
        try:
            if step.done_check():
                draw_progress(idx, total, step.label, "SKIPPED")
                continue
            draw_progress(idx - 1 if idx > 1 else 0, total, step.label, "RUNNING")
            step.action()
            draw_progress(idx, total, step.label, "DONE")
        except Exception as exc:
            if step.label == "Assemble with Flye" and args.assembly_retry_all_qc:
                print()
                warn(
                    f"Flye failed for {sample_name} with the preselected reads. "
                    "Retrying with all QC-passed reads and Flye meta mode."
                )
                try:
                    rescued = rescue_sample_with_all_qc_reads(
                        layout,
                        shared_layout,
                        sample_name,
                        args.min_q,
                        args.threads,
                        args.min_pident,
                        args.min_qcov,
                        args.flye_iterations,
                    )
                except Exception as rescue_exc:
                    print()
                    die(f"Step failed: Rescue assembly with all QC reads\nReason: {rescue_exc}")
                draw_progress(total, total, "Rescue assembly with all QC reads", "DONE")
                print_status_line("Rescue result", "RECOVERED" if rescued else "NO TARGET CONTIGS", "yellow" if rescued else "red")
                rescue_performed = True
                break
            print()
            die(f"Step failed: {step.label}\nReason: {exc}")

    target_count = count_fasta_records(step8_outputs(layout)[0])
    if target_count == 0 and args.assembly_retry_all_qc and not rescue_performed:
        print()
        warn(
            f"No target contigs were retained for {sample_name}. "
            "Retrying assembly with all QC-passed reads and Flye meta mode."
        )
        try:
            rescued = rescue_sample_with_all_qc_reads(
                layout,
                shared_layout,
                sample_name,
                args.min_q,
                args.threads,
                args.min_pident,
                args.min_qcov,
                args.flye_iterations,
            )
        except Exception as exc:
            print()
            die(f"Step failed: Rescue assembly with all QC reads\nReason: {exc}")
        status = "DONE" if rescued else "DONE"
        detail = "RECOVERED" if rescued else "NO TARGET CONTIGS"
        draw_progress(total, total, "Rescue assembly with all QC reads", status)
        print_status_line("Rescue result", detail, "yellow" if rescued else "red")

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
                 shared_layout, sample_layouts, args.threads, args.mafft_adjust_direction, args.min_contig_len_phylo
             )),
        Step("Trim alignment with trimAl",
             lambda: step10_done(shared_layout),
             lambda: step10_trimal(
                 shared_layout,
                 args.trimal_gap_threshold,
                 args.min_alignment_occupancy,
                 args.min_assembled_core_occupancy,
             )),
    ]
    if args.prepare_beast2:
        steps.append(
            Step("Prepare optional BEAST 2 input files",
                 lambda: (
                     beast2_done(shared_layout)
                     and args.beast2_manual_dates is None
                     and args.beast2_coordinates is None
                 ),
                 lambda: step13_prepare_beast2(
                     shared_layout,
                     args.taxid,
                     args.beast2_manual_dates,
                     args.beast2_coordinates,
                 ))
        )
    if args.run_beast2:
        steps.append(
            Step("Run BEAST 2 and summarize MCC tree",
                 lambda: beast2_run_done(shared_layout),
                 lambda: step14_run_beast2(
                     shared_layout,
                     args.beast2_xml,
                     args.threads,
                     args.beast2_burnin,
                     args.beast2_chain_length,
                     args.beast2_log_every,
                 ))
        )
    steps.extend([
        Step("Infer batch ML tree and render PDF",
             lambda: step12_tree_done(shared_layout),
             lambda: step12_iqtree(shared_layout, args.threads)),
        Step("Build pairwise identity heatmap",
             lambda: step11_identity_done(shared_layout),
             lambda: step11_identity_plot(shared_layout, args.identity_plot_min)),
    ])

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
    identity_tsv, identity_pdf, _, assembled_identity_tsv, assembled_identity_pdf, _ = step11_identity_outputs(shared_layout)
    treefile, _, tree_pdf = step12_tree_outputs(shared_layout)

    print()
    print_banner("Batch Completed", f"{len(batch_rows)} sample(s) processed successfully")
    print_status_line("Global summary", summarize_path(batch_summary_txt), "green")
    print_status_line("Global table", summarize_path(batch_summary_tsv), "green")
    print_status_line("Trimmed alignment", summarize_path(trimmed_alignment), "green")
    print_status_line("Identity matrix", summarize_path(identity_tsv), "green")
    print_status_line("Global identity heatmap", summarize_path(identity_pdf), "green")
    print_status_line("Assembled-only identity matrix", summarize_path(assembled_identity_tsv), "green")
    print_status_line("Assembled-only heatmap", summarize_path(assembled_identity_pdf), "green")
    print_status_line("Global tree", summarize_path(treefile), "green")
    print_status_line("Global tree PDF", summarize_path(tree_pdf), "green")
    if args.prepare_beast2:
        beast2 = beast2_outputs(shared_layout)
        print_status_line("BEAST 2 folder", summarize_path(beast2["outdir"]), "green")
        print_status_line("BEAST 2 manual dates", summarize_path(beast2["manual_dates_template"]), "yellow")
        print_status_line("BEAST 2 map coordinates", summarize_path(beast2["locations_template"]), "yellow")
    if args.run_beast2:
        beast2_run = beast2_run_outputs(shared_layout)
        print_status_line("BEAST 2 run folder", summarize_path(beast2_run["outdir"]), "green")
        print_status_line("BEAST 2 MCC tree", summarize_path(beast2_run["mcc_tree"]), "green")
        print_status_line("BEAST 2 phylogeography map", summarize_path(beast2_run["map_pdf"]), "green")
        print_status_line("BEAST 2 time tree", summarize_path(beast2_run["time_tree_pdf"]), "green")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()
