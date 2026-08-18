import csv, re, os
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
import sqlite3
from lxml import etree
import pickle
import pandas as pd

class Taxon:
    # Hierarchies related to taxonomy classification
    basic_hierarchy = [
        'superkingdom', 'clade', 'kingdom', 'phylum', 'subphylum', 'class', 'subclass',
        'order', 'suborder', 'family', 'subfamily', 'genus', 'subgenus',
        'species', 'subspecies'
    ]

    major_hierarchy = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    additional_categories = ['genotype', 'isolate', 'no rank', 'serogroup', 'serotype']
    
    def __init__(self, taxid, rank, name, parent=None):
        if not isinstance(taxid, int):
            raise TypeError("taxid must be set to an integer")
        self.taxid = int(taxid)
        self.rank = rank
        self.name = name
        self.common_name = None
        self.parent = parent
        self.children = []
        if parent:
            parent.children.append(self)

    def __lt__(self, other):
        return self.name.lower() < other.name.lower() 
    
    def __repr__(self):
        return f"Taxon(taxid={self.taxid}, rank='{self.rank}', name='{self.name}')"
    
    def set_common_name(self, common_name):
        self.common_name = common_name

    def add_child(self, child):
        self.children.append(child)
        child.parent = self
        
    @lru_cache(maxsize=None)
    def get_lineage(self, include_self=True):
        lineage = [self] if include_self else []
        current = self
        while current.parent:
            lineage.insert(0, current.parent)
            current = current.parent
        return lineage
        
    def get_major_lineage_line(self):
        # Initialize a dictionary with all ranks set to an empty string
        lineage_dict = {rank: '' for rank in self.major_hierarchy}
        lineage = self.get_lineage(include_self=True)
        for ancestor in lineage:
            if ancestor.rank in lineage_dict:
                lineage_dict[ancestor.rank] = f"{ancestor.rank[0:1]}__{ancestor.name}"
        lineage_str = '; '.join([lineage_dict[rank] for rank in self.major_hierarchy])
        return lineage_str
    
    def get_all_lineage(self):
        lineage = []
        current = self
        while current:
            lineage.append(f"{current.rank[0:2]}_{current.name}")
            current = current.parent
        return '; '.join(reversed(lineage))
    
    def get_all_lineage_as_taxid(self):
        """
        Return lineage of this Taxon as a string of taxid values prefixed by rank initials.

        Example: 'su_1; ph_1234; cl_5678'
        """
        lineage = []
        current = self
        while current:
            lineage.append(f"{current.rank}_{current.taxid}")
            current = current.parent
        return '; '.join(reversed(lineage))

class NcbiTaxon(Taxon):
    """Sub Class"""

    def __init__(self, taxid, rank, name, scientific_name, parent=None, division= None):
        super().__init__(taxid, rank, name, parent)
        self.scientific_name = scientific_name
        self.division = division
        self.host = None        
        self.other_names = []  # Initialize the list of other names  
        self.func_name =  None # Initialize the list of other names        
        
    def set_other_names(self, other_names):
        self.other_names = other_names
        
    def set_host(self,host):
        self.host = host
    
    def set_functional_name(self, func_name):
        self.func_name = func_name
        
    def print(self):
        return(f'name ={self.name}, rank={self.rank}, other_names={self.other_names}, host={self.host}, div={self.division}, func_name={self.func_name}')
    

class TaxonUtils:
    def __init__(self):
        self.taxon_dict = {}    
        
    def parse_taxon_line(self, line):
        match = re.match(r'(\d+) \[([^\]]+)\] (.+)', line)
        if match:
            taxid, rank, name = match.groups()
            return Taxon(int(taxid), rank, name)
        else:
            raise ValueError(f"Line format is incorrect: {line}")
        
    def build_taxon_structure(self, d, parent=None):
        for line, children in d.items():
            current_taxon = self.parse_taxon_line(line)
            self.taxon_dict[current_taxon.taxid] = current_taxon
            if parent is None:
                parent = current_taxon  # Initialize the root on the first call
            else:
                parent.add_child(current_taxon)
            self.build_taxon_structure(children, current_taxon)
        return parent   
 
    def find_by_taxid(self, taxid):
        return self.taxon_dict.get(taxid, None)
    
    def get_taxon_dict(self):
        return self.taxon_dict
    
    def find_majority_taxon(self, list_taxa, level):
        """Find the majority taxon at a given taxonomic level considering the entire lineage."""
        count = defaultdict(int)
        found_level = False        
        for taxon in list_taxa:
            lineage = taxon.get_lineage()
            for ancestor in reversed(lineage):
                if ancestor.rank == level:
                    count[ancestor] += 1
                    found_level = True
                    break
        if not found_level:
            second_level  = lineage[-1].rank 
            for taxon in list_taxa:
                lineage = taxon.get_lineage()
                for ancestor in reversed(lineage):
                    if ancestor.rank == second_level:
                        count[ancestor] += 1
                        break
        return max(count, key=count.get) if count else None

    @staticmethod
    def collect_all_levels(taxon, tax_rank=None):
        if tax_rank is None:
            tax_rank = []
        tax_rank.append(taxon.rank)
        for child in taxon.children:
            TaxonUtils.collect_all_levels(child, tax_rank)
        return tax_rank
    
    @staticmethod
    def collect_all_taxids(taxon):
        result = [taxon.taxid]
        stack = [taxon]
        while stack:
            current = stack.pop()
            result.append(current.taxid)
            stack.extend(current.children)
        return result
        
    def process_ncbi_tax_xml_to_taxon(self, xml_input):
        context = etree.iterparse(xml_input, tag='Taxon')
        temp_parent_dict = {}  # Dictionary to store temporary parent information

        for event, taxon in context:
            if taxon.getparent().tag != 'LineageEx':
                taxid = int(taxon.findtext('TaxId'))
                name = taxon.findtext('ScientificName')
                rank = taxon.findtext('Rank')
                division = taxon.findtext('Division')
                parent_tax_id = int(taxon.findtext('ParentTaxId'))
                other_names_list = taxon.find('OtherNames')
                other_names = []
                if other_names_list is not None:
                    other_names = [child.text for child in other_names_list if child is not None]
                current_taxon = NcbiTaxon(taxid=taxid, rank=rank, name=name, scientific_name=name, parent=None, division=division)
                current_taxon.set_other_names(other_names)
                self.taxon_dict[taxid] = current_taxon
                temp_parent_dict[taxid] = parent_tax_id
                taxon.clear()
                while taxon.getprevious() is not None:
                    del taxon.getparent()[0]
                    
        for taxid, parent_tax_id in temp_parent_dict.items():
            parent = self.taxon_dict.get(parent_tax_id)
            current_taxon = self.taxon_dict[taxid]  # Access existing NcbiTaxon object
            current_taxon.parent = parent
    
        return [taxon for taxon in self.taxon_dict.values() if taxon.parent is None]
  
    def find_by_name(self, s_pattern: str, include_other_name= bool):
        matching_taxa = []
        search_pattern = re.compile(s_pattern, flags=re.IGNORECASE)  # Case-insensitive search
        for taxon in self.taxon_dict.values():
            if search_pattern.search(taxon.name) or search_pattern.search(taxon.scientific_name):
                matching_taxa.append(taxon)
            elif include_other_name and taxon.other_names:  # Check other names conditionally
                if any(search_pattern.search(other_name) for other_name in taxon.other_names):
                    matching_taxa.append(taxon)                        
        return matching_taxa



class Acc2TaxidLoader:
    # download nucl_gb.accession2taxid.gz from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
    
    def __init__(self, acc2taxid_file, virus_acc2taxid_file, pickle_file):
        self.acc2taxid_file = acc2taxid_file
        self.virus_acc2taxid_file = virus_acc2taxid_file
        self.pickle_file = pickle_file
        self.virus_acc2taxid = None

    def load(self, all_taxids=None):
        if self._is_cached():
            return self.virus_acc2taxid

        if self._pickle_exists():
            self._load_from_pickle()
        else:
            self._load_and_process_csv(all_taxids)
            self._save_to_pickle()

        return self.virus_acc2taxid

    def _is_cached(self):
        return self.virus_acc2taxid is not None

    def _pickle_exists(self):
        return os.path.exists(self.pickle_file)

    def _load_from_pickle(self):
        with open(self.pickle_file, 'rb') as pkl_file:
            self.virus_acc2taxid = pickle.load(pkl_file)
        print(f"Loaded acc2taxid dictionary from {self.pickle_file}.")

    def _load_and_process_csv(self, all_taxids):
        if os.path.exists(self.virus_acc2taxid_file):
            virus_acc2taxid = pd.read_csv(self.virus_acc2taxid_file, compression='gzip', header=0, sep='\t', index_col=0)
        else:
            acc2taxid = pd.read_csv(self.acc2taxid_file, compression='gzip', header=0, sep='\t')
            virus_acc2taxid = acc2taxid[acc2taxid['taxid'].isin(all_taxids)]
            virus_acc2taxid.to_csv(self.virus_acc2taxid_file, sep='\t', compression='gzip')
        
        self.virus_acc2taxid = dict(zip(virus_acc2taxid['accession'], virus_acc2taxid['taxid']))

    def _save_to_pickle(self):
        with open(self.pickle_file, 'wb') as pkl_file:
            pickle.dump(self.virus_acc2taxid, pkl_file)
        print(f"Created and saved acc2taxid dictionary to {self.pickle_file}.")

class TaxonomyProcessor:
    def __init__(self, sqlite3_file, db_name):
        self.sqlite3_file = sqlite3_file
        self.db_name = db_name

    def select_from_taxids(self, taxid_list, chunk_size=50):
        base_query = "SELECT taxid, func_name, original_name, division, new_division FROM {} WHERE taxid IN ({})"
        chunks = [taxid_list[i:i + chunk_size] for i in range(0, len(taxid_list), chunk_size)]

        dataframes = []
        with sqlite3.connect(self.sqlite3_file) as conn:
            for chunk in chunks:
                placeholders = ','.join('?' * len(chunk))
                query = base_query.format(self.db_name, placeholders)
                df_chunk = pd.read_sql_query(query, conn, params=chunk)
                dataframes.append(df_chunk)
        return pd.concat(dataframes, ignore_index=True)
    
    def dict_to_dataframe(self, sample_dict, sample_name):
        df = pd.DataFrame.from_dict(sample_dict, orient='index')
        df['sample'] = sample_name
        return df.reset_index().rename(columns={'index': 'taxon'})

    def process_dataframe(self, result_dict):
        dataframes = [self.dict_to_dataframe(sd, name) for name, sd in result_dict.items()]
        full_df = pd.concat(dataframes)
        unique_taxids = full_df['taxid'].unique().tolist()
        taxid_data = self.select_from_taxids(unique_taxids)
        full_df = self.merge_taxid_data(full_df, taxid_data)
        return full_df

    def merge_taxid_data(self, full_df, taxid_data):
        full_df = full_df.merge(taxid_data[['taxid', 'func_name', 'original_name', 'division']], on='taxid', how='left')
        full_df['func_name'].fillna(full_df['original_name'], inplace=True)
        full_df['simple_name'] = full_df['func_name']
         # full_df['simple_name'] = full_df.apply(
        #     lambda row: 'Streptococcus phage' if 'streptococcus phage' in row['func_name'].lower() else
        #                 ('Other phage' if 'phage' in row['func_name'].lower() or row['division'] == 'Phages' else row['simple_name']),
        #     axis=1
        # )
        return full_df

    def aggregate_and_save(self, df, group_by_column, output_path):
        aggregated_df = df.groupby(['sample', group_by_column]).agg(count_per_sample=('count', 'sum')).reset_index()
        pivot_df = aggregated_df.pivot_table(index=[group_by_column], columns='sample', values='count_per_sample', fill_value=0)
        pivot_df.to_csv(output_path, index=True, sep="\t")

    def agg_and_save(self, df, group_by_columns, output_path):
        # group_by_columns가 string이면 리스트로 변환
        if isinstance(group_by_columns, str):
            group_by_columns = [group_by_columns]
            
        aggregated_df = df.groupby(['sample'] + group_by_columns).agg(count_per_sample=('count', 'sum')).reset_index()
        
        pivot_df = aggregated_df.pivot_table(
            index=group_by_columns, columns='sample', values='count_per_sample', fill_value=0
        )
        pivot_df.to_csv(output_path, index=True, sep="\t")
        
    def to_kraken_report_from_df(self, df):
        """
        full_df만을 사용해 kraken-style report 생성
        원하는 컬럼: clade_proportion, clade_count, taxon_count, rank, taxID, name
        """
        rows = []
        total_reads = df['count'].sum()
        for idx, row in df.iterrows():
            clade_proportion = f"{100 * row['count'] / total_reads:.2f}" if total_reads > 0 else "0.00"
            clade_count = row['count']
            taxon_count = row['count']   # 실제 direct count와 clade_count가 다르다면 구분 필요
            rank = row['rank'] if 'rank' in row and pd.notnull(row['rank']) else 'U'
            taxID = row['taxid'] if 'taxid' in row else row.get('taxID', '')
            name = row.get('func_name') or row.get('original_name') or row.get('taxon') or row.get('name', '')
            rows.append([clade_proportion, clade_count, taxon_count, rank, taxID, name])
        # 정렬: clade_count (내림차순)
        rows = sorted(rows, key=lambda x: -int(x[1]))
        return rows

    def parse_lineage_tax(self, lineage_tax):
        """
        lineage_tax: 'phylum_1234; class_5678; ...; species_9999'
        return: [("phylum", 1234), ("class", 5678), ..., ("species", 9999)]
        """
        tokens = [x.strip() for x in str(lineage_tax).split(';') if x]
        return [(t.split('_',1)[0], int(t.split('_',1)[1])) for t in tokens if '_' in t]

    def build_taxid_tree(self, full_df):
        """
        full_df에는 반드시 'taxid'와 'lineage_tax' 컬럼이 있어야 함
        return: parent_map, children_map, taxid_name_map
        """
        from collections import defaultdict
        parent_map = {}     # taxid -> parent_taxid
        children_map = defaultdict(list)  # parent_taxid -> [child_taxid, ...]
        taxid_name_map = {} # taxid -> name
        for idx, row in full_df.iterrows():
            taxid = row['taxid']
            taxid_name_map[taxid] =  row.get('original_name', row.get('name', ''))
            lineage = self.parse_lineage_tax(row['lineage_tax'])
            if len(lineage) > 1:
                parent_taxid = lineage[-2][1]
            else:
                parent_taxid = None  # root
            parent_map[taxid] = parent_taxid
            if parent_taxid is not None:
                children_map[parent_taxid].append(taxid)
        return parent_map, children_map, taxid_name_map

    def kraken_report_from_lineage_tax(self, full_df, output_dir, suffix=""):
        parent_map, children_map, taxid_name_map = self.build_taxid_tree(full_df)
        sample_names = full_df['sample'].unique()
        for sample in sample_names:
            sample_df = full_df[full_df['sample'] == sample].copy()
            total_reads = sample_df['count'].sum()
            taxid2row = {row['taxid']: row for idx, row in sample_df.iterrows()}
            all_taxids = set(sample_df['taxid'])
            roots = [taxid for taxid in all_taxids if parent_map[taxid] not in all_taxids]
            columns = ["percent", "reads", "direct", "rank", "taxid", "name"]
            lines = []
            def dfs(taxid, depth=0):
                row = taxid2row[taxid]
                reads = row['count']
                direct = reads  # 실제로 자식 누적 합산이 필요하면 여기서 수정
                percent = 100 * reads / total_reads if total_reads > 0 else 0
                lineage = self.parse_lineage_tax(row['lineage_tax'])
                rank = lineage[-1][0] if lineage else "U"
                name = taxid_name_map.get(taxid, "")
                indent = "  " * depth
                lines.append([
                    f"{percent:.2f}", str(reads), str(direct),
                    rank[0].upper(),  # 한글자
                    str(taxid), indent + name
                ])
                for child_taxid in children_map.get(taxid, []):
                    dfs(child_taxid, depth + 1)
            for root in roots:
                dfs(root, 0)
            # 파일 저장
            report_file = os.path.join(output_dir, f"{sample}{'_' + suffix if suffix else ''}.report")
            with open(report_file, "w") as fout:
                fout.write('\t'.join(columns) + '\n')
                for r in lines:
                    fout.write('\t'.join(map(str, r)) + '\n')
            print(f"Saved: {report_file}")
        
    def process_and_save_aggregations(self, full_df, output_dir, suffix):
        lineage_file = os.path.join(output_dir, suffix + '_lineage.tsv')
        self.aggregate_and_save(full_df, 'lineage', lineage_file)
        
        lineage_file = os.path.join(output_dir, suffix + '_lineage_tax.tsv')
        self.aggregate_and_save(full_df, 'lineage_tax', lineage_file)

        func_name_file = os.path.join(output_dir, suffix + '_func.tsv')
        self.aggregate_and_save(full_df, 'func_name', func_name_file)

        simple_name_file = os.path.join(output_dir, suffix + '_func_simple.tsv')
        self.aggregate_and_save(full_df, 'simple_name', simple_name_file)
        
        ltax_filename = os.path.join(output_dir, suffix + '_lineage_with_tax.tsv')
        self.agg_and_save(full_df, ['lineage','lineage_tax'], ltax_filename)
        
        print("kraken")
        self.kraken_report_from_lineage_tax(full_df, output_dir, suffix)
