#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Virome MEGABLAST assignment pipeline

This script performs viral read-level assignment using BLAST/MEGABLAST results
and taxonomy mapping files. It was used after shotgun metagenomic preprocessing
and host-read removal.

Workflow:
    Step 1:
        Run BLAST/MEGABLAST against a custom viral reference database.

    Step 2:
        Parse BLAST outputs, map accessions to taxids, assign taxonomy,
        and generate taxonomic aggregation tables.

    Step 0:
        Parse and aggregate existing BLAST outputs without running BLAST again.

Notes:
    - Absolute local server paths are not included in this public version.
    - Required database and taxonomy files should be provided through arguments.
    - This script requires the accompanying modules:
        blast_handler.py
        taxonomy.py
"""

import argparse
import os
from datetime import datetime

from blast_handler import BlastProcessor
from blast_handler import BlastParser
from taxonomy import Acc2TaxidLoader
from taxonomy import TaxonomyProcessor


def die(msg, code=1):
    print(msg)
    raise SystemExit(code)


def find_pairs(inputdir, findex, rindex=None):
    """
    Find paired or single-end input FASTA files.

    Parameters
    ----------
    inputdir : str
        Directory containing input FASTA files.
    findex : str
        Pattern identifying forward reads, e.g., "_1" or "_F".
    rindex : str or None
        Pattern identifying reverse reads, e.g., "_2" or "_R".

    Returns
    -------
    file_pairs : list
        List of tuples containing paired files or single files.
    unmatched_forward : list
        Forward files without matched reverse files.
    """
    files = os.listdir(inputdir)
    forward_files = [f for f in files if findex in f]
    file_pairs = []
    unmatched_forward = []

    if rindex:
        reverse_files = [f for f in files if rindex in f]

        for f_file in forward_files:
            expected_reverse_file = f_file.replace(findex, rindex)

            if expected_reverse_file in reverse_files:
                file_pairs.append(
                    (
                        os.path.abspath(os.path.join(inputdir, f_file)),
                        os.path.abspath(os.path.join(inputdir, expected_reverse_file)),
                    )
                )
            else:
                unmatched_forward.append(
                    os.path.abspath(os.path.join(inputdir, f_file))
                )

        return file_pairs, unmatched_forward

    for f_file in forward_files:
        file_pairs.append((os.path.abspath(os.path.join(inputdir, f_file)),))

    return file_pairs, []


def get_cpu_cores():
    """
    Return the number of available CPU cores.
    """
    num_cores = os.cpu_count()

    if num_cores is None:
        print("Unable to determine the number of CPU cores. Using 1 core.")
        return 1

    return num_cores


def parse_steps(step_string):
    """
    Parse step argument.

    Examples
    --------
    "1-2" -> [1, 2]
    "1,2" -> [1, 2]
    "0"   -> [0]
    """
    if "-" in step_string:
        start, end = map(int, step_string.split("-"))
        return list(range(start, end + 1))

    return list(map(int, step_string.split(",")))


def write_log(out_dir, args):
    """
    Write a run log with all parameters.
    """
    os.makedirs(out_dir, exist_ok=True)

    log_file = os.path.join(
        out_dir,
        "run_log_" + datetime.now().strftime("%Y%m%d_%H%M%S") + ".txt",
    )

    with open(log_file, "a", encoding="utf-8") as log:
        log.write(f"Run date and time: {datetime.now()}\n")
        log.write("Parameters:\n")

        for arg, value in vars(args).items():
            log.write(f"  {arg}: {value}\n")

        log.write("\n")


def require_existing_path(path, label):
    """
    Check whether a required path exists.
    """
    if path is None:
        die(f"[x] Missing required argument: {label}")

    if not os.path.exists(path):
        die(f"[x] File or directory not found for {label}: {path}")


def validate_taxonomy_arguments(args):
    """
    Validate taxonomy-related files required for BLAST parsing and aggregation.
    """
    require_existing_path(args.acc2taxid, "--acc2taxid")
    require_existing_path(args.virus_acc2taxid, "--virus-acc2taxid")
    require_existing_path(args.virus_acc2taxid_pickle, "--virus-acc2taxid-pickle")
    require_existing_path(args.taxonkit_pickle, "--taxonkit-pickle")
    require_existing_path(args.taxonomy_db, "--taxonomy-db")


def run_blast_step(args, steps):
    """
    Step 1: Run BLAST/MEGABLAST against the custom viral reference database.
    """
    if 1 not in steps:
        return

    if not args.inputdir or not args.outdir or not args.findex:
        die("[x] Step 1 requires --inputdir, --outdir, and --findex.")

    require_existing_path(args.db, "--db")

    file_pairs, unmatched_forward = find_pairs(
        inputdir=args.inputdir,
        findex=args.findex,
        rindex=args.rindex,
    )

    print("File pairs found:", file_pairs)

    if unmatched_forward:
        print("Unmatched forward files:", unmatched_forward)

    if not file_pairs:
        die("[x] No input files were detected. Check --inputdir, --findex, and --rindex.")

    num_cores = get_cpu_cores()
    blast_worker = max(1, min(max(1, args.worker // 10), max(1, num_cores // 10)))

    print("BLAST type:", args.blast_type)
    print("BLAST worker processes:", blast_worker)

    blast_processor = BlastProcessor(args.db, args.blast_type)
    blast_processor.process_fasta_multithread(
        file_pairs,
        args.outdir,
        max_workers=blast_worker,
    )

    print(f"STEP 1 DONE: {blast_processor}")


def run_parse_and_aggregate_step(args, steps, parse_step):
    """
    Step 0 or Step 2:
    Parse BLAST outputs and generate taxonomy aggregation tables.

    Step 0 parses existing BLAST outputs without best-hit or paired-read options.
    Step 2 performs the main post-BLAST parsing and aggregation.
    """
    if parse_step not in steps:
        return

    validate_taxonomy_arguments(args)

    qcoverage = args.min_query_coverage
    pidentity = args.min_percent_identity

    acc2taxid_loader = Acc2TaxidLoader(
        args.acc2taxid,
        args.virus_acc2taxid,
        args.virus_acc2taxid_pickle,
    )

    if parse_step == 0:
        blast_parser = BlastParser(
            args.outdir,
            acc2taxid_loader,
            args.taxonkit_pickle,
            qcov=qcoverage,
            pident=pidentity,
        )

        file_suffix = (
            datetime.today().strftime("%Y%m%d")
            + "_"
            + str(qcoverage)
            + "_"
            + str(pidentity)
        )

    else:
        blast_parser = BlastParser(
            args.outdir,
            acc2taxid_loader,
            args.taxonkit_pickle,
            best_hit_only=args.best_hit_only,
            qcov=qcoverage,
            pident=pidentity,
            intersect=args.paired_match_only,
        )

        file_suffix = (
            datetime.today().strftime("%Y%m%d")
            + "_"
            + str(qcoverage)
            + "_"
            + str(pidentity)
            + "_intersect_"
            + str(args.paired_match_only)
        )

    blast_parser.process_all_files(
        findex=args.findex,
        rindex=args.rindex,
        max_workers=args.worker,
    )

    taxonomy_processor = TaxonomyProcessor(
        args.taxonomy_db,
        args.taxonomy_db_name,
    )

    full_df = taxonomy_processor.process_dataframe(blast_parser.result_dict)
    taxonomy_processor.process_and_save_aggregations(
        full_df,
        args.outdir,
        file_suffix,
    )

    print(f"STEP {parse_step} DONE: taxonomy aggregation completed.")


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Run MEGABLAST/BLASTN, parse viral BLAST results, "
            "map accessions to taxids, and generate taxonomic aggregation tables."
        )
    )

    parser.add_argument(
        "-i",
        "--inputdir",
        type=str,
        required=True,
        help="Directory containing input FASTA files.",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        required=True,
        help="Directory to save output files.",
    )

    parser.add_argument(
        "-f",
        "--findex",
        type=str,
        required=True,
        help="Pattern identifying forward-read files, e.g., _1 or _F.",
    )

    parser.add_argument(
        "-r",
        "--rindex",
        type=str,
        required=False,
        default=None,
        help="Pattern identifying reverse-read files, e.g., _2 or _R.",
    )

    parser.add_argument(
        "-d",
        "--db",
        type=str,
        required=False,
        default=None,
        help="Path to the custom viral BLAST database.",
    )

    parser.add_argument(
        "-w",
        "--worker",
        type=int,
        default=10,
        help="Number of worker processes.",
    )

    parser.add_argument(
        "-m",
        "--midentity",
        "--min-query-coverage",
        dest="min_query_coverage",
        type=float,
        default=75,
        required=False,
        help=(
            "Minimum query coverage threshold for valid viral hits. "
            "The legacy option name --midentity is retained for compatibility."
        ),
    )

    parser.add_argument(
        "-c",
        "--mcoverage",
        "--min-percent-identity",
        dest="min_percent_identity",
        type=float,
        default=90,
        required=False,
        help=(
            "Minimum percent identity threshold for valid viral hits. "
            "The legacy option name --mcoverage is retained for compatibility."
        ),
    )

    parser.add_argument(
        "-s",
        "--steps",
        type=str,
        required=False,
        default="1-2",
        help=(
            "Steps to execute. Use '1-2' to run BLAST and parsing, "
            "'1' to run BLAST only, '2' to parse only, or '0' to parse existing outputs."
        ),
    )

    parser.add_argument(
        "-p",
        "--paired-match-only",
        action="store_true",
        help="Keep only matching forward-reverse read pairs during parsing.",
    )

    parser.add_argument(
        "-b",
        "--best-hit-only",
        action="store_true",
        help="Keep only best BLAST hit per read during parsing.",
    )

    parser.add_argument(
        "-t",
        "--blast-type",
        default="megablast",
        help="BLAST type to run, e.g., megablast or blastn.",
    )

    parser.add_argument(
        "--acc2taxid",
        type=str,
        required=False,
        default=None,
        help="Path to acc2taxid.tsv.",
    )

    parser.add_argument(
        "--virus-acc2taxid",
        dest="virus_acc2taxid",
        type=str,
        required=False,
        default=None,
        help="Path to virus_acc2taxid.tsv.gz.",
    )

    parser.add_argument(
        "--virus-acc2taxid-pickle",
        dest="virus_acc2taxid_pickle",
        type=str,
        required=False,
        default=None,
        help="Path to virus_acc2taxid.pkl.",
    )

    parser.add_argument(
        "--taxonkit-pickle",
        dest="taxonkit_pickle",
        type=str,
        required=False,
        default=None,
        help="Path to taxonkit.pkl.",
    )

    parser.add_argument(
        "--taxonomy-db",
        dest="taxonomy_db",
        type=str,
        required=False,
        default=None,
        help="Path to virus_taxonomy.db.",
    )

    parser.add_argument(
        "--taxonomy-db-name",
        dest="taxonomy_db_name",
        type=str,
        required=False,
        default="taxonomy",
        help="SQLite taxonomy database name.",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    steps = parse_steps(args.steps)

    os.makedirs(args.outdir, exist_ok=True)
    write_log(args.outdir, args)

    print("Steps:", steps)
    print("Paired matches only:", args.paired_match_only)
    print("Best hits only:", args.best_hit_only)
    print("Minimum query coverage:", args.min_query_coverage)
    print("Minimum percent identity:", args.min_percent_identity)

    run_parse_and_aggregate_step(args, steps, parse_step=0)
    run_blast_step(args, steps)
    run_parse_and_aggregate_step(args, steps, parse_step=2)

    print("[✓] Done.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An error occurred: {e}")
        raise
