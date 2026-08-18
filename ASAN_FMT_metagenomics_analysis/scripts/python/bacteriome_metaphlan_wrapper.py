#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
import shlex
import argparse
import subprocess
from glob import glob
from pathlib import Path
from shutil import which
from typing import Optional, List, Tuple, Dict, Set

# Accept: fastq, fq, fasta, fa, fna (+ .gz)
FASTQ_EXTS = set([".fastq", ".fq"])
FASTA_EXTS = set([".fasta", ".fa", ".fna"])
GZ_EXT = ".gz"


def die(msg, code=1):
    print(msg, file=sys.stderr)
    raise SystemExit(code)


def tool_check():
    for cmd in ("metaphlan", "merge_metaphlan_tables.py"):
        if which(cmd) is None:
            die("[x] '{}' not found in PATH. Activate the correct environment.".format(cmd))


def split_ext(p):
    """
    Returns (base_stem, main_ext, is_gz)
    main_ext excludes .gz; e.g. sample.fastq.gz -> (sample, .fastq, True)
    Python 3.7 compatible.
    """
    name = p.name
    if name.endswith(GZ_EXT):
        p2 = Path(name[:-len(GZ_EXT)])
        return p2.stem, p2.suffix.lower(), True
    return p.stem, p.suffix.lower(), False


def infer_input_type(main_ext):
    if main_ext in FASTQ_EXTS:
        return "fastq"
    if main_ext in FASTA_EXTS:
        return "fasta"
    die("[x] Unsupported extension: {}. Allowed: fastq/fq/fasta/fa/fna (+ .gz)".format(main_ext))


def list_candidate_files(input_dir):
    patterns = []
    for ext in list(FASTQ_EXTS | FASTA_EXTS):
        patterns.append(os.path.join(input_dir, "*{}".format(ext)))
        patterns.append(os.path.join(input_dir, "*{}{}".format(ext, GZ_EXT)))

    files = []
    for pat in patterns:
        files.extend(glob(pat))

    files = sorted(set(files))
    return [Path(x) for x in files]


def mate_stem(stem, mate):
    # Convert "..._R1" <-> "..._R2" or "..._1" <-> "..._2"
    if re.search(r"_R[12]$", stem):
        return re.sub(r"_R[12]$", "_R{}".format(mate), stem)
    if re.search(r"_[12]$", stem):
        return re.sub(r"_[12]$", "_{}".format(mate), stem)
    return "{}_{}".format(stem, mate)


def detect_pe_se(files):
    """
    Return:
      pe_samples: list of (sample_name, r1_path, r2_path, input_type)
      se_samples: list of (sample_name, read_path, input_type)

    Pairing rules:
      - If stem ends with _R1/_R2 or _1/_2, treat as PE candidate.
      - Only pair if both mates exist and have same main extension and gz symmetry.
      - Remaining files are SE.
    """
    by_name = dict((f.name, f) for f in files)  # type: Dict[str, Path]
    used = set()  # type: Set[str]
    pe = []  # type: List[Tuple[str, str, str, str]]
    se = []  # type: List[Tuple[str, str, str]]

    # Pairing pass
    for f in files:
        if f.name in used:
            continue

        stem, main_ext, is_gz = split_ext(f)

        m = re.match(r"^(?P<core>.+?)(?:_R)?(?P<mate>[12])$", stem)
        if not m:
            continue

        core = m.group("core")
        mate = m.group("mate")
        other_mate = "2" if mate == "1" else "1"

        other_stem = mate_stem(stem, other_mate)
        other_name = other_stem + main_ext + (GZ_EXT if is_gz else "")

        if other_name not in by_name:
            continue

        f2 = by_name[other_name]
        stem2, main_ext2, is_gz2 = split_ext(f2)

        if main_ext2 != main_ext or is_gz2 != is_gz:
            continue

        input_type = infer_input_type(main_ext)

        if mate == "1":
            r1, r2 = f, f2
        else:
            r1, r2 = f2, f

        sample_name = core
        pe.append((sample_name, str(r1), str(r2), input_type))
        used.add(r1.name)
        used.add(r2.name)

    # SE pass
    for f in files:
        if f.name in used:
            continue
        stem, main_ext, _ = split_ext(f)
        input_type = infer_input_type(main_ext)
        sample_name = stem
        se.append((sample_name, str(f), input_type))

    return pe, se


def run_metaphlan(
    sample_name,
    inputs,
    input_type,
    out_dir,
    nproc,
    bowtie2db,
    index,
    thresholds,
    extra_args,
    skip_existing,
):
    os.makedirs(out_dir, exist_ok=True)
    out_profile = os.path.join(out_dir, "{}_profile.txt".format(sample_name))
    out_bt2 = os.path.join(out_dir, "{}.bowtie2.bz2".format(sample_name))

    if skip_existing and os.path.exists(out_profile):
        print("[=] Skip {} (exists)".format(sample_name))
        return out_profile

    if len(inputs) == 1:
        input_arg = inputs[0]
    else:
        input_arg = "{},{}".format(inputs[0], inputs[1])

    cmd = [
        "metaphlan",
        input_arg,
        "--input_type", input_type,
        "--nproc", str(nproc),
        "--bowtie2out", out_bt2,
        "-o", out_profile,
    ]

    if bowtie2db:
        cmd += ["--bowtie2db", bowtie2db]
    if index:
        cmd += ["--index", index]

    # threshold-like knobs
    if thresholds:
        cmd += thresholds

    # pass-through args
    if extra_args:
        cmd += extra_args

    mode = "PE" if len(inputs) == 2 else "SE"
    print("[+] MetaPhlAn: {} ({}, {})".format(sample_name, mode, input_type))
    subprocess.run(cmd, check=True)
    return out_profile


def merge_profiles(out_dir, merged_name, overwrite):
    profs = sorted(glob(os.path.join(out_dir, "*_profile.txt")))
    if not profs:
        die("[x] No *_profile.txt found to merge.")

    merged_path = os.path.join(out_dir, merged_name)
    cmd = ["merge_metaphlan_tables.py"] + profs + ["-o", merged_path]
    if overwrite:
        cmd += ["--overwrite"]

    print("[+] Merging {} profiles -> {}".format(len(profs), merged_path))
    subprocess.run(cmd, check=True)
    return merged_path


def main():
    p = argparse.ArgumentParser(
        description="Run MetaPhlAn on SE/PE inputs with fastq/fq/fasta/fa/fna (+.gz), then merge profiles."
    )
    p.add_argument("-i", "--input_dir", required=True)
    p.add_argument("-o", "--output_dir", required=True)
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--skip-existing", action="store_true")

    # DB pinning
    p.add_argument("--bowtie2db", default=None)
    p.add_argument("--index", default=None)

    # Threshold knobs
    p.add_argument("--read_min_len", type=int, default=None, help="MetaPhlAn: minimum read length")
    p.add_argument("--stat_q", type=float, default=None, help="MetaPhlAn: stat_q threshold controlling clade reporting")
    p.add_argument("--min_ab", type=float, default=None, help="MetaPhlAn: minimum abundance threshold (analysis-type dependent)")

    # Any extra MetaPhlAn args pass-through
    p.add_argument(
        "--metaphlan-args",
        default="",
        help="Extra args passed to metaphlan as a single string. Example: '--unclassified_estimation -t rel_ab_w_read_stats'",
    )

    # Merge output
    p.add_argument("--merged_file", default="merged_metaphlan.tsv")
    p.add_argument("--merge-overwrite", action="store_true")

    args = p.parse_args()
    tool_check()

    files = list_candidate_files(args.input_dir)
    if not files:
        die("[x] No supported input files found in: {}".format(args.input_dir))

    pe, se = detect_pe_se(files)

    if not pe and not se:
        die("[x] No samples detected. Check naming and extensions.")

    thresholds = []  # type: List[str]
    if args.read_min_len is not None:
        thresholds += ["--read_min_len", str(args.read_min_len)]
    if args.stat_q is not None:
        thresholds += ["--stat_q", str(args.stat_q)]
    if args.min_ab is not None:
        thresholds += ["--min_ab", str(args.min_ab)]

    extra_args = shlex.split(args.metaphlan_args) if args.metaphlan_args.strip() else []

    print("[i] Detected PE samples: {}".format(len(pe)))
    print("[i] Detected SE samples: {}".format(len(se)))

    # Run PE then SE
    for sample, r1, r2, input_type in pe:
        run_metaphlan(
            sample_name=sample,
            inputs=[r1, r2],
            input_type=input_type,
            out_dir=args.output_dir,
            nproc=args.threads,
            bowtie2db=args.bowtie2db,
            index=args.index,
            thresholds=thresholds,
            extra_args=extra_args,
            skip_existing=args.skip_existing,
        )

    for sample, r, input_type in se:
        run_metaphlan(
            sample_name=sample,
            inputs=[r],
            input_type=input_type,
            out_dir=args.output_dir,
            nproc=args.threads,
            bowtie2db=args.bowtie2db,
            index=args.index,
            thresholds=thresholds,
            extra_args=extra_args,
            skip_existing=args.skip_existing,
        )

    merge_profiles(args.output_dir, args.merged_file, overwrite=args.merge_overwrite)
    print("[✓] Done.")


if __name__ == "__main__":
    main()

