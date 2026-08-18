from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import pickle
import subprocess
import re
import csv
import time
import pandas as pd
import pytaxonkit as ptk
from taxonomy import TaxonUtils


class BlastProcessor:

    def __init__(self, db, blast_type='megablast'):
        self.db = db
        blastdb_dir = os.path.dirname(db)
        if os.path.exists(db):
            os.environ['BLASTDB'] = blastdb_dir
        self.blast_type = blast_type

    def run_makeblastdb(self):
        if not os.path.exists(self.db + ".ndb") or not os.path.exists(self.db + ".nhr"):
            cmd = f'makeblastdb -in {self.db} -dbtype nucl'
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, encoding='utf-8')
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    print(output.strip())
            process.wait()

    def run_megablast(self, fasta_file, out_file):
        if not os.path.exists(fasta_file):
            print(f"No such file: {fasta_file}")
            return None

        cmd = (
            f'blastn -db {self.db} -num_threads 12 -query {fasta_file} '
            f'-task {self.blast_type} -out {out_file} '
            f'-qcov_hsp_perc 50 -perc_identity 50 -parse_deflines '
            f'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send '
            f'evalue bitscore sgi saccver slen qlen staxids sscinames stitle qcovs qcovhsp"'
        )
        print(cmd)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, encoding='utf-8')
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output.strip())
        return out_file

    def is_fasta(self, file_path):
        try:
            with open(file_path, 'r') as file:
                first_line = file.readline().strip()
                return first_line.startswith('>')
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            return False

    def process_fasta(self, fasta_file, output_dir):
        print(fasta_file)
        if not self.is_fasta(fasta_file):
            print(f"File {fasta_file} is not a valid FASTA file. Skipping.")
            return None

        out_file = os.path.join(output_dir, os.path.basename(fasta_file) + ".blast6.tsv")
        self.run_megablast(fasta_file, out_file)
        return out_file

    def flatten_fasta_files(self, fasta_files):
        """Flatten a list of tuples into a flat list of files."""
        flat_list = []
        for item in fasta_files:
            if isinstance(item, tuple):
                flat_list.extend(item)
            else:
                flat_list.append(item)
        return flat_list

    def process_fasta_multithread(self, fasta_files, output_dir, max_workers=4):
        results = []
        self.run_makeblastdb()

        flat_fasta_files = self.flatten_fasta_files(fasta_files)

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_fasta = {
                executor.submit(self.process_fasta, fasta, output_dir): fasta
                for fasta in flat_fasta_files
            }
            for future in as_completed(future_to_fasta):
                fasta_file = future_to_fasta[future]
                try:
                    result = future.result()
                    results.append(result)
                    print(f"Completed processing: {fasta_file}")
                except Exception as exc:
                    print(f"FASTA file {fasta_file} generated an exception: {exc}")

        return results


class BlastParser:
    EXCLUDE_PATTERN = re.compile(r'MH590414\.[0-9]|MH590388\.[0-9]|MG934417\.[0-9]|OR951374\.[0-9]')
    EXCLUDE_PATTERN2 = re.compile(r'Uncultured human fecal virus|HIV|Immunodeficiency')
    EXCLUDE_PATTERN3 = re.compile(r'NC_032111\.[0-9]')

    def __init__(
        self,
        bout_dir,
        acc2taxid_loader,
        taxon_structure_pickle='taxon_structure.pkl',
        best_hit_only=False,
        qcov=75,
        pident=90,
        intersect=True,
    ):
        self.bout_dir = bout_dir
        self.acc2taxid_loader = acc2taxid_loader
        self.taxon_structure_pickle = taxon_structure_pickle
        self.utils = self.initialize_taxon_utils()
        self.result_dict = {}
        self.best_hit_only = best_hit_only
        self.qcov = qcov
        self.pident = pident
        self.intersect = intersect

    # ------------------------------------------------------------------ #
    # Taxonomy helpers
    # ------------------------------------------------------------------ #

    def initialize_taxon_utils(self):
        taxon_structure = self.load_taxon_structure()
        taxon_utils = TaxonUtils()
        taxon_utils.build_taxon_structure(taxon_structure)
        return taxon_utils

    def load_taxon_structure(self):
        if os.path.exists(self.taxon_structure_pickle):
            print(f"Loading taxon structure from {self.taxon_structure_pickle}")
            with open(self.taxon_structure_pickle, 'rb') as f:
                taxon_structure = pickle.load(f)
        else:
            print("Generating taxon structure...")
            taxon_structure = ptk.list([10239], threads=20, raw=True)
            with open(self.taxon_structure_pickle, 'wb') as f:
                pickle.dump(taxon_structure, f)
            print(f"Taxon structure saved to {self.taxon_structure_pickle}")
        return taxon_structure

    def load_acc2taxid(self):
        return self.acc2taxid_loader.load()

    # ------------------------------------------------------------------ #
    # BLAST parsing
    # ------------------------------------------------------------------ #

    def _should_skip_row(self, row):
        try:
            sstart = int(row[8])
            send = int(row[9])
        except (IndexError, ValueError):
            return True

        in_exclude_range = 8000 <= sstart <= 8600 and 8000 <= send <= 8600
        is_excluded_seq = self.EXCLUDE_PATTERN3.search(row[1])

        return (
            float(row[19]) < self.qcov
            or float(row[2]) < self.pident
            or is_excluded_seq
            or self.EXCLUDE_PATTERN2.search(row[18])
            or (is_excluded_seq and in_exclude_range)
        )

    def _get_query_id(self, qseqid: str) -> str:
        rid = qseqid
        core = re.sub(r"/[12fr]$", "", rid)
        return core

    def parse_blast(self, blast_output_file, acc2taxid, test=False):
        print(blast_output_file)
        results = defaultdict(list)
        with open(blast_output_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                if self._should_skip_row(row):
                    continue

                match_info = {
                    'sseqid': row[1],
                    'staxids': row[16].split(';')[0],
                    'sscinames': row[17].split(';')[0],
                    'stitle': row[18].split(';')[0],
                    'qcov': float(row[19]),
                    'pident': float(row[2]),
                    'slen': int(row[14]),
                }
                query_id = self._get_query_id(row[0])

                if int(match_info["staxids"]) == 0:
                    acc = re.sub(r"\.[0-9]+", "", match_info["sseqid"])
                    match_info['staxids'] = int(acc2taxid.get(acc, 0))
                    if match_info["staxids"] == 0:
                        continue

                results[query_id].append(match_info)

                if test and len(results) >= 10:
                    print(query_id)
                    break
        return results

    def append_lineage(self, tutils, hits_list, level='species'):
        try:
            list_taxa = []
            for hit in hits_list:
                found_taxon = tutils.find_by_taxid(int(hit['staxids']))
                if found_taxon:
                    list_taxa.append(found_taxon)
            majority_taxon = tutils.find_majority_taxon(list_taxa, level)
            return majority_taxon
        except Exception as err:
            print("Unexpected error:", err, type(err))

    def update_sample_dict(self, sample_dict, taxon):
        if taxon not in sample_dict:
            sample_dict[taxon] = {
                'count': 0,
                'taxid': taxon.taxid,
                'lineage': taxon.get_all_lineage(),
                'lineage_tax': taxon.get_all_lineage_as_taxid(),
            }
        sample_dict[taxon]['count'] += 1

    # ------------------------------------------------------------------ #
    # Paired / single BLAST result handling
    # ------------------------------------------------------------------ #

    def process_paired_blast_results(self, fwd_result, rev_result):
        """
        Paired-end 모드: forward / reverse 둘 다 있을 때만 호출.
        """
        qcov = self.qcov
        pident = self.pident
        intersect = self.intersect

        uniq_read_uids = set(fwd_result.keys()).union(rev_result.keys())
        categories = {}

        for read_uid in uniq_read_uids:
            fwd_hits = fwd_result.get(read_uid, [])
            rev_hits = rev_result.get(read_uid, [])

            for_hits_subjid = {
                hit['sseqid'] for hit in fwd_hits
                if hit['pident'] >= pident and hit['qcov'] >= qcov
            }
            rev_hits_subjid = {
                hit['sseqid'] for hit in rev_hits
                if hit['pident'] >= pident and hit['qcov'] >= qcov
            }

            if intersect:
                common_sseqids = for_hits_subjid & rev_hits_subjid
                valid_hits = [
                    hit for hit in (fwd_hits + rev_hits)
                    if hit['sseqid'] in common_sseqids
                    and hit['pident'] >= pident
                    and hit['qcov'] >= qcov
                ]
            else:
                valid_hits = [
                    hit for hit in (fwd_hits + rev_hits)
                    if hit['pident'] >= pident and hit['qcov'] >= qcov
                ]

            if valid_hits:
                categories[read_uid] = {'hits': valid_hits, 'mtype': 2}
        return categories

    def process_single_blast_results(self, fwd_result):
        """
        Single-end / contig 모드: forward 결과만 사용.
        """
        qcov = self.qcov
        pident = self.pident

        categories = {}
        for read_uid, hits in fwd_result.items():
            valid_hits = [
                hit for hit in hits
                if hit['pident'] >= pident and hit['qcov'] >= qcov
            ]
            if valid_hits:
                categories[read_uid] = {'hits': valid_hits, 'mtype': 1}
        return categories

    def process_file_pair2(self, sample_id, fwd_blast_out_file, rev_blast_out_file=None):
        """
        sample_id: 파일/샘플 식별자
        fwd_blast_out_file: forward BLAST TSV
        rev_blast_out_file: reverse BLAST TSV (없을 수 있음; None 이거나 파일이 없으면 single 모드)
        """
        print(f'Processing: {fwd_blast_out_file} and {rev_blast_out_file}')

        start_time = time.time()
        acc2taxid = self.load_acc2taxid()

        # Forward parsing
        fwd_blast_result = self.parse_blast(fwd_blast_out_file, acc2taxid)

        # Paired vs single 결정
        if rev_blast_out_file is not None and os.path.exists(rev_blast_out_file):
            rev_blast_result = self.parse_blast(rev_blast_out_file, acc2taxid)
            categorized_hits = self.process_paired_blast_results(fwd_blast_result, rev_blast_result)
        else:
            categorized_hits = self.process_single_blast_results(fwd_blast_result)
            exec_time = time.time() - start_time
            print(f"Execution time (forward-only): {exec_time} seconds")

        # 유효한 히트가 하나도 없으면 에러 대신 스킵
        if not categorized_hits:
            print(f"[WARN] No valid hits for sample {sample_id}; skipping.")
            return None

        # Taxon majority / count 집계
        paired_result = {}
        for ruid, hits in categorized_hits.items():
            try:
                major_taxon = self.append_lineage(self.utils, hits['hits'], 'species')
                if major_taxon:
                    self.update_sample_dict(paired_result, major_taxon)
            except Exception as e:
                print(f"Error processing hit for {ruid}: {e}")
                continue

        return sample_id, paired_result

    def process_all_files(self, findex="_kneaddata_paired_R1", rindex="_kneaddata_paired_R2", max_workers=4):
        """
        findex: forward 파일 이름 패턴 (예: '_kneaddata_paired_R1', '.unmapped_1', '.blast6.tsv' 등)
        rindex: reverse 파일 이름 패턴 (예: '_kneaddata_paired_R2', '.unmapped_2').
                None 이면 single/contig 모드로 처리.
        """
        tasks = []

        try:
            n_files = len(
                [
                    name
                    for name in os.listdir(self.bout_dir)
                    if os.path.isfile(os.path.join(self.bout_dir, name))
                ]
            )
            print(n_files)
        except Exception:
            pass

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for parent, _, files in os.walk(self.bout_dir):
                for f in files:
                    # BLAST6 TSV만 처리
                    if ".blast6.tsv" not in f:
                        continue

                    # rindex가 주어져 있으면 reverse 파일은 여기서 건너뜀
                    if rindex and (rindex in f):
                        continue

                    idx = os.path.basename(f)
                    if findex:
                        idx = idx.replace(findex, "")

                    for_bout = os.path.join(parent, f)

                    rev_bout = None
                    if rindex:
                        rev_f = os.path.basename(f).replace(findex, rindex)
                        rev_bout_candidate = os.path.join(parent, rev_f)
                        if os.path.exists(rev_bout_candidate):
                            rev_bout = rev_bout_candidate
                        else:
                            print(
                                f"[WARN] Reverse BLAST file not found for {f} "
                                f"(expected: {rev_f}). Using forward-only."
                            )

                    tasks.append(
                        executor.submit(self.process_file_pair2, idx, for_bout, rev_bout)
                    )

            for future in as_completed(tasks):
                result = future.result()
                if result:
                    idx, pair_result = result
                    self.result_dict[idx] = pair_result

        print("All files processed.")

