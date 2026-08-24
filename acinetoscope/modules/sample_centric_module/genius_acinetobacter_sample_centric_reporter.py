#!/usr/bin/env python3
"""
GENIUS ACINETOBACTER SAMPLE‑CENTRIC REPORTER
Interactive per‑isolate boxes for AMR, Virulence, BACMET, Plasmids, Mutations
Gene‑centric for MLST, Kaptive, QC, Combinations, Co‑occurrence, Patterns

Version: 1.1.0
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Date: 2026-08-22
"""

import os
import sys
import json
import re
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup


class SampleCentricParser:
    """Parse HTML reports and TSV files for AcinetoScope sample‑centric analysis."""

    def __init__(self):
        self.abricate_tsv_files = {
            'resfinder': 'acineto_resfinder_abricate_summary.tsv',
            'vfdb': 'acineto_vfdb_abricate_summary.tsv',
            'card': 'acineto_card_abricate_summary.tsv',
            'argannot': 'acineto_argannot_abricate_summary.tsv',
            'bacmet2': 'acineto_bacmet2_abricate_summary.tsv',
            'plasmidfinder': 'acineto_plasmidfinder_abricate_summary.tsv',
            'megares': 'acineto_megares_abricate_summary.tsv',
            'ncbi': 'acineto_ncbi_abricate_summary.tsv',
            'ecoh': 'acineto_ecoh_abricate_summary.tsv',
            'ecoli_vf': 'acineto_ecoli_vf_abricate_summary.tsv',
        }

    def normalize_sample_id(self, sample_id: str) -> str:
        sample = str(sample_id)
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
                break
        return sample.strip()

    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            table = tables[table_index]
            rows = table.find_all('tr')
            headers = [th.get_text().strip() for th in rows[0].find_all(['th', 'td'])]
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    while len(row_data) < len(headers):
                        row_data.append('')
                    if len(row_data) > len(headers):
                        row_data = row_data[:len(headers)]
                    data.append(row_data)
            if not data:
                return pd.DataFrame()
            df = pd.DataFrame(data, columns=headers)
            df.columns = [col.strip().replace('\n', ' ') for col in df.columns]
            return df
        except Exception as e:
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()

    # ------------------------------------------------------------------------
    # TSV loaders
    # ------------------------------------------------------------------------
    def load_amrfinder_tsv(self, input_dir: Path) -> Tuple[Dict[str, List[Dict]], Dict[str, int]]:
        path = input_dir / 'acineto_amrfinder_summary.tsv'
        if not path.exists():
            return {}, {}
        df = pd.read_csv(path, sep='\t')
        details = defaultdict(list)
        counts = Counter()
        for _, row in df.iterrows():
            sample = row['Genome']
            gene_dict = row.to_dict()
            gene_dict = {k: (None if pd.isna(v) else v) for k, v in gene_dict.items()}
            if 'Element symbol' in gene_dict:
                gene_dict['gene'] = gene_dict['Element symbol']
            elif 'Element_symbol' in gene_dict:
                gene_dict['gene'] = gene_dict['Element_symbol']
            elif 'Gene' in gene_dict:
                gene_dict['gene'] = gene_dict['Gene']
            details[sample].append(gene_dict)
            counts[gene_dict.get('gene', 'unknown')] += 1
        return dict(details), dict(counts)

    def load_abricate_tsv(self, input_dir: Path) -> Tuple[Dict[str, Dict[str, List[Dict]]], Dict[str, Dict[str, int]]]:
        details = defaultdict(lambda: defaultdict(list))
        counts = defaultdict(lambda: defaultdict(int))
        for db, fname in self.abricate_tsv_files.items():
            path = input_dir / fname
            if not path.exists():
                continue
            df = pd.read_csv(path, sep='\t')
            sample_col = 'genome' if 'genome' in df.columns else 'file'
            for _, row in df.iterrows():
                sample = str(row[sample_col])
                sample = re.sub(r'\.(fasta|fna)$', '', sample)
                gene_dict = row.to_dict()
                gene_dict = {k: (None if pd.isna(v) else v) for k, v in gene_dict.items()}
                if 'gene' not in gene_dict and 'Gene' in gene_dict:
                    gene_dict['gene'] = gene_dict['Gene']
                details[sample][db].append(gene_dict)
                gene = gene_dict.get('gene', 'unknown')
                counts[db][gene] += 1
        return dict(details), dict(counts)

    def load_apt_tsv(self, input_dir: Path) -> Dict[str, Dict]:
        path = input_dir / 'plasmid_summary.tsv'
        if not path.exists():
            return {}
        df = pd.read_csv(path, sep='\t')
        apt = {}
        for _, row in df.iterrows():
            genome = row['Genome']
            hit_count = row.get('HitCount', 0)
            rep_types = row.get('RepTypes', '')
            if pd.isna(rep_types) or rep_types == '':
                rep_types = []
            else:
                rep_types = [r.strip() for r in rep_types.split(',') if r.strip()]
            apt[genome] = {'hit_count': int(hit_count) if not pd.isna(hit_count) else 0, 'rep_types': rep_types}
        return apt

    def load_mutations_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        path = input_dir / 'mutation_summary.tsv'
        if not path.exists():
            return {}
        df = pd.read_csv(path, sep='\t')
        muts = defaultdict(list)
        for _, row in df.iterrows():
            sample = row['genome']
            def clean(v):
                return '' if pd.isna(v) else str(v)
            muts[sample].append({
                'gene': clean(row.get('gene_symbol', '')),
                'mutation': clean(row.get('element_name', '')),
                'class': clean(row.get('class', '')),
                'subclass': clean(row.get('subclass', '')),
                'contig': clean(row.get('contig_id', '')),
                'start': clean(row.get('start', '')),
                'stop': clean(row.get('stop', '')),
                'strand': clean(row.get('strand', '')),
                'coverage': clean(row.get('coverage', '')),
                'identity': clean(row.get('identity', '')),
                'accession': clean(row.get('accession', ''))
            })
        return dict(muts)

    # ------------------------------------------------------------------------
    # HTML parsers (gene‑centric)
    # ------------------------------------------------------------------------
    def parse_mlst_report(self, file_path: Path, scheme: str = "pasteur") -> Dict[str, Dict]:
        print(f"  🧬 Parsing {scheme.upper()} MLST: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                df = self.parse_html_table(html, 1)
                if df.empty:
                    return {}
            df.columns = [col.strip().replace('\n', ' ') for col in df.columns]
            sample_col = None
            for col in df.columns:
                if 'sample' in col.lower() or 'genome' in col.lower() or 'id' in col.lower() or 'strain' in col.lower():
                    if '#' not in col.lower():
                        sample_col = col
                        break
            if not sample_col and len(df.columns) > 1:
                sample_col = df.columns[0]
            if not sample_col:
                print(f"    ⚠️ Could not find sample column in {scheme} MLST")
                return {}
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            results = {}
            st_col = None
            for col in df.columns:
                if col.upper() == 'ST' or col.lower() == 'sequence type':
                    st_col = col
                    break
            if st_col is None:
                for col in df.columns:
                    if 'st' in col.lower() and 'sample' not in col.lower():
                        st_col = col
                        break
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                if pd.isna(sample) or not sample:
                    continue
                st = 'ND'
                if st_col is not None and pd.notna(row.get(st_col)):
                    st_val = str(row[st_col]).strip()
                    if st_val and st_val.upper() not in ['', 'NAN', 'NONE', 'ND', 'UNKNOWN', 'STUNKNOWN']:
                        match = re.search(r'ST(\d+)', st_val, re.IGNORECASE)
                        if match:
                            st = match.group(1)
                        else:
                            num_match = re.search(r'\d+', st_val)
                            if num_match:
                                st = num_match.group()
                else:
                    for col in df.columns:
                        if col.lower() != sample_col.lower() and pd.notna(row.get(col)):
                            cell_val = str(row[col])
                            match = re.search(r'ST(\d+)', cell_val, re.IGNORECASE)
                            if match:
                                st = match.group(1)
                                break
                ic = 'Unknown'
                if scheme.lower() == 'pasteur':
                    ic_col = None
                    for col in df.columns:
                        if 'international clone' in col.lower() or col.lower() == 'ic' or col.lower() == 'clone':
                            ic_col = col
                            break
                    if ic_col is not None and pd.notna(row.get(ic_col)):
                        ic_val = str(row[ic_col]).strip()
                        if ic_val and ic_val.lower() not in ['', 'nan', 'none', 'unknown', 'not assigned']:
                            roman_to_arabic = {'I':'1','II':'2','III':'3','IV':'4','V':'5','VI':'6','VII':'7','VIII':'8','IX':'9','X':'10'}
                            ic_match = re.search(r'IC\s*([IVXLCDM]+)', ic_val, re.I)
                            if ic_match:
                                roman = ic_match.group(1).upper()
                                arabic = roman_to_arabic.get(roman, roman)
                                ic = f"IC{arabic}"
                            else:
                                ic = ic_val.replace(' ', '')
                    if ic == 'Unknown':
                        for col in df.columns:
                            if col.lower() != sample_col.lower() and pd.notna(row.get(col)):
                                cell_val = str(row[col])
                                roman_to_arabic = {'I':'1','II':'2','III':'3','IV':'4','V':'5','VI':'6','VII':'7','VIII':'8','IX':'9','X':'10'}
                                ic_match = re.search(r'IC\s*([IVXLCDM]+)', cell_val, re.I)
                                if ic_match:
                                    roman = ic_match.group(1).upper()
                                    arabic = roman_to_arabic.get(roman, roman)
                                    ic = f"IC{arabic}"
                                    break
                allele_profile = 'ND'
                for col in df.columns:
                    if 'allele' in col.lower() and 'profile' in col.lower():
                        if pd.notna(row.get(col)):
                            allele_val = str(row[col]).strip()
                            if allele_val and allele_val.lower() not in ['', 'nan', 'none', 'nd']:
                                allele_profile = allele_val
                        break
                results[sample] = {
                    'ST': st,
                    'International_Clone': ic,
                    'Allele_Profile': allele_profile,
                    'Scheme': scheme
                }
            print(f"    ✓ Found {len(results)} samples for {scheme.upper()} scheme")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing {scheme} MLST: {e}")
            return {}

    def parse_kaptive_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing Kaptive: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                df = self.parse_html_table(html, 1)
                if df.empty:
                    return {}
            df.columns = [col.strip() for col in df.columns]
            sample_col = None
            for col in df.columns:
                if any(k in col.lower() for k in ['genome', 'sample', 'id', 'strain']):
                    sample_col = col
                    break
            if not sample_col and len(df.columns) > 0:
                sample_col = df.columns[0]
            if not sample_col:
                return {}
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                if pd.isna(sample) or not sample:
                    continue
                k_locus = 'ND'
                for k_col in ['K Locus', 'K locus', 'K-Locus', 'K-locus', 'K', 'KLocus']:
                    if k_col in df.columns and pd.notna(row.get(k_col)):
                        k_val = str(row[k_col]).strip()
                        if k_val and k_val.lower() not in ['', 'nan', 'none', 'nd', 'unknown']:
                            k_match = re.search(r'(K\d+)', k_val, re.I)
                            if k_match:
                                k_locus = k_match.group(1).upper()
                            else:
                                unknown_match = re.search(r'unknown\s*\(KL(\d+)\)', k_val, re.I)
                                if unknown_match:
                                    k_locus = f"K{unknown_match.group(1)}"
                                else:
                                    parts = k_val.split()
                                    for part in parts:
                                        if part.startswith('K') and part[1:].isdigit():
                                            k_locus = part.upper()
                                            break
                                        elif part.startswith('KL') and part[2:].isdigit():
                                            k_locus = 'K' + part[2:]
                                            break
                            if k_locus != 'ND':
                                break
                o_locus = 'ND'
                for o_col in ['O Locus', 'O locus', 'O-Locus', 'O-locus', 'O', 'OLocus']:
                    if o_col in df.columns and pd.notna(row.get(o_col)):
                        o_val = str(row[o_col]).strip()
                        if o_val and o_val.lower() not in ['', 'nan', 'none', 'nd', 'unknown']:
                            ocl_match = re.search(r'(OCL\d+)', o_val, re.I)
                            if ocl_match:
                                o_locus = ocl_match.group(1).upper()
                            else:
                                num_match = re.search(r'(\d+)', o_val)
                                if num_match:
                                    o_locus = f"OCL{num_match.group(1)}"
                                else:
                                    parts = o_val.split()
                                    for part in parts:
                                        if part.startswith('OCL') and len(part) > 3 and part[3:].isdigit():
                                            o_locus = part.upper()
                                            break
                                        elif part.startswith('OC') and len(part) > 2 and part[2:].isdigit():
                                            o_locus = f"OCL{part[2:]}"
                                            break
                                        elif part.startswith('O') and len(part) > 1 and part[1:].isdigit():
                                            o_locus = f"OCL{part[1:]}"
                                            break
                            if o_locus != 'ND':
                                break
                capsule_type = f"{k_locus}:{o_locus}"
                results[sample] = {
                    'K_Locus': k_locus,
                    'O_Locus': o_locus,
                    'Capsule_Type': capsule_type
                }
            print(f"    ✓ Found {len(results)} samples with capsule typing")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing Kaptive: {e}")
            return {}

    def parse_qc_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing FASTA QC: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}
            sample_col = None
            for col in df.columns:
                if 'filename' in col.lower() or 'sample' in col.lower() or col == df.columns[0]:
                    sample_col = col
                    break
            if not sample_col:
                return {}
            results = {}
            for _, row in df.iterrows():
                sample_raw = row[sample_col]
                if not sample_raw:
                    continue
                sample = self.normalize_sample_id(sample_raw)
                qc = {}
                for col in df.columns:
                    if col == sample_col:
                        continue
                    val = row[col]
                    if pd.isna(val) or val == '' or val == 'ND':
                        qc[col] = 'ND'
                    else:
                        cleaned = str(val).replace('%', '').replace(',', '').strip()
                        try:
                            qc[col] = float(cleaned)
                        except:
                            qc[col] = str(val)
                results[sample] = qc
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing QC: {e}")
            return {}


class SampleCentricAnalyzer:
    """Build gene‑centric frequencies, cross‑genome patterns, and typing maps."""

    def __init__(self):
        self.critical_carbapenemases = {
            'blaOXA-23', 'blaOXA-24', 'blaOXA-40', 'blaOXA-51', 'blaOXA-58',
            'blaOXA-66', 'blaOXA-69', 'blaOXA-71', 'blaOXA-143', 'blaOXA-235',
            'blaNDM', 'blaNDM-1', 'blaNDM-2', 'blaNDM-3', 'blaNDM-4', 'blaNDM-5',
            'blaVIM', 'blaVIM-1', 'blaVIM-2', 'blaVIM-3', 'blaVIM-4',
            'blaIMP', 'blaIMP-1', 'blaIMP-2', 'blaIMP-3', 'blaIMP-4', 'blaIMP-5',
            'blaKPC', 'blaKPC-2', 'blaKPC-3', 'blaKPC-4',
            'blaGES', 'blaGES-1', 'blaGES-2', 'blaGES-5', 'blaGES-14',
            'blaSIM', 'blaSPM', 'blaAIM'
        }
        self.critical_esbls = {
            'blaCTX-M', 'blaSHV', 'blaPER', 'blaVEB', 'blaBEL', 'blaGES'
        }
        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'pmrA', 'pmrB', 'pmrC', 'pmrE', 'pmrF', 'lpxA', 'lpxC', 'lpxD', 'lpxL', 'lpxM',
            'eptA', 'arnA', 'arnB', 'arnC', 'arnD', 'arnE', 'arnF'
        }
        self.critical_tigecycline = {
            'tet(X)', 'tet(X1)', 'tet(X2)', 'tet(X3)', 'tet(X4)', 'tet(X5)', 'tet(X6)',
            'tet(39)', 'tet(A)', 'tet(B)', 'tet(C)', 'tet(D)', 'tet(E)', 'tet(G)', 'tet(H)',
            'adeS', 'adeR', 'adeA', 'adeB', 'adeC', 'adeJ', 'adeK', 'adeN', 'adeT'
        }
        self.critical_biofilm = {
            'ompA', 'csuA', 'csuB', 'csuC', 'csuD', 'csuE', 'csuA/B', 'csuC/D/E',
            'bfmR', 'bfmS', 'abaI', 'abaR', 'pilA', 'pilB', 'pilC', 'pilD', 'pilE', 'pilF',
            'ptk', 'epsA', 'pgaA', 'pgaB', 'pgaC', 'pgaD', 'bap'
        }
        self.critical_efflux = {
            'adeA', 'adeB', 'adeC', 'adeF', 'adeG', 'adeH', 'adeI', 'adeJ', 'adeK', 'adeL', 'adeM', 'adeN',
            'abeM', 'abeS', 'amvA', 'adeT1', 'adeT2', 'mexJ', 'mexK', 'mexT'
        }
        self.bacmet2_markers = {
            'qacA', 'qacB', 'qacC', 'qacD', 'qacE', 'qacF', 'qacG', 'qacH', 'qacI', 'qacJ',
            'qacEA1', 'qacG2', 'qacH2', 'cepA', 'formA', 'formB', 'formC', 'oqxA', 'oqxB',
            'czcA', 'czcB', 'czcC', 'czcD', 'czcR', 'czcS',
            'merA', 'merB', 'merC', 'merD', 'merE', 'merF', 'merG', 'merH', 'merI', 'merJ', 'merP', 'merT',
            'arsA', 'arsB', 'arsC', 'arsD', 'arsE', 'arsF', 'arsG', 'arsH', 'arsI', 'arsJ', 'arsT',
            'copA', 'copB', 'copC', 'copD', 'copE', 'copF', 'copG', 'copH', 'copI', 'copJ',
            'zntA', 'zntB', 'zntC', 'zntD', 'zntE', 'zntF', 'zntG', 'zntH', 'zntI', 'zntJ',
            'chrA', 'chrB', 'chrC', 'chrD', 'chrE', 'chrF',
            'nikA', 'nikB', 'nikC', 'nikD', 'nikE', 'nikR',
            'cadA', 'cadB', 'cadC', 'cadD',
            'silA', 'silB', 'silC', 'silD', 'silE',
            'pbrA', 'pbrB', 'pbrC', 'pbrD', 'pbrR',
            'corA', 'corC', 'corR', 'zraR', 'pitA', 'nccN', 'nreB', 'fptA', 'fecE', 'fpvA', 'znuB', 'znuC', 'frnE'
        }

    def categorize_gene(self, gene: str) -> str:
        gene_lower = gene.lower()
        if any(carb in gene_lower for carb in [g.lower() for g in self.critical_carbapenemases]):
            return 'Carbapenemases'
        elif any(esbl in gene_lower for esbl in [g.lower() for g in self.critical_esbls]):
            return 'ESBLs'
        elif any(col in gene_lower for col in [g.lower() for g in self.critical_colistin]):
            return 'Colistin Resistance'
        elif any(tig in gene_lower for tig in [g.lower() for g in self.critical_tigecycline]):
            return 'Tigecycline Resistance'
        elif any(bio in gene_lower for bio in [g.lower() for g in self.critical_biofilm]):
            return 'Biofilm Formation'
        elif any(eff in gene_lower for eff in [g.lower() for g in self.critical_efflux]):
            return 'Efflux Pumps'
        elif any(bac in gene_lower for bac in [g.lower() for g in self.bacmet2_markers]):
            return 'Bacmet (Biocide/Metal)'
        else:
            return 'Other'

    def build_gene_centric_frequencies(self,
                                       amr_details: Dict[str, List[Dict]],
                                       abricate_details: Dict[str, Dict[str, List[Dict]]],
                                       total_samples: int) -> Dict[str, Any]:
        gene_centric = {
            'amr_databases': {},
            'virulence_databases': {},
            'bacmet_databases': {},
            'plasmid_databases': {},
            'gene_categories': defaultdict(list),
        }

        if amr_details:
            all_amr = []
            for genes in amr_details.values():
                all_amr.extend(genes)
            freq = {}
            for g in all_amr:
                gene = g.get('gene', 'unknown')
                freq[gene] = freq.get(gene, 0) + 1
            amr_list = []
            for gene, count in freq.items():
                category = self.categorize_gene(gene)
                percentage = (count / total_samples) * 100 if total_samples else 0
                amr_list.append({
                    'gene': gene,
                    'category': category,
                    'database': 'AMRfinder',
                    'count': count,
                    'percentage': round(percentage, 2),
                    'frequency_display': f"{count} ({percentage:.1f}%)",
                    'genomes': []
                })
            if amr_list:
                amr_list.sort(key=lambda x: x['count'], reverse=True)
                gene_centric['amr_databases']['amrfinder'] = amr_list
                for gd in amr_list:
                    gene_centric['gene_categories'][gd['category']].append(gd)

        amr_dbs = ['card', 'resfinder', 'argannot', 'megares', 'ncbi']
        vir_dbs = ['vfdb', 'ecoli_vf']
        bac_dbs = ['bacmet2']
        plas_dbs = ['plasmidfinder', 'ecoh']

        db_to_genes = defaultdict(list)
        for sample, dbs in abricate_details.items():
            for db, genes in dbs.items():
                db_to_genes[db].extend(genes)

        for db, genes in db_to_genes.items():
            if not genes:
                continue
            freq = {}
            for g in genes:
                gene = g.get('gene', 'unknown')
                freq[gene] = freq.get(gene, 0) + 1
            gene_list = []
            for gene, count in freq.items():
                category = self.categorize_gene(gene)
                percentage = (count / total_samples) * 100 if total_samples else 0
                gene_list.append({
                    'gene': gene,
                    'category': category,
                    'database': db.upper(),
                    'count': count,
                    'percentage': round(percentage, 2),
                    'frequency_display': f"{count} ({percentage:.1f}%)",
                    'genomes': []
                })
            if gene_list:
                gene_list.sort(key=lambda x: x['count'], reverse=True)
                if db in amr_dbs:
                    gene_centric['amr_databases'][db] = gene_list
                elif db in vir_dbs:
                    gene_centric['virulence_databases'][db] = gene_list
                elif db in bac_dbs:
                    gene_centric['bacmet_databases'][db] = gene_list
                elif db in plas_dbs:
                    gene_centric['plasmid_databases'][db] = gene_list
                for gd in gene_list:
                    gene_centric['gene_categories'][gd['category']].append(gd)

        return gene_centric

    def build_cross_genome_patterns(self,
                                    samples_data: Dict[str, Dict],
                                    gene_centric: Dict[str, Any]) -> Dict[str, Any]:
        patterns = {
            'pasteur_st_distribution': Counter(),
            'oxford_st_distribution': Counter(),
            'k_locus_distribution': Counter(),
            'o_locus_distribution': Counter(),
            'capsule_type_distribution': Counter(),
            'st_k_combinations': defaultdict(list),
            'st_o_combinations': defaultdict(list),
            'ko_combinations': defaultdict(list),
            'st_ko_combinations': defaultdict(list),
            'high_risk_crab': [],
            'gene_cooccurrence': defaultdict(Counter)
        }
        for sample, data in samples_data.items():
            pasteur = data.get('pasteur_mlst', {}).get('ST', 'ND')
            oxford = data.get('oxford_mlst', {}).get('ST', 'ND')
            k = data.get('kaptive', {}).get('K_Locus', 'ND')
            o = data.get('kaptive', {}).get('O_Locus', 'ND')
            capsule = data.get('kaptive', {}).get('Capsule_Type', 'ND')
            if pasteur != 'ND':
                patterns['pasteur_st_distribution'][pasteur] += 1
            if oxford != 'ND':
                patterns['oxford_st_distribution'][oxford] += 1
            if k != 'ND':
                patterns['k_locus_distribution'][k] += 1
            if o != 'ND':
                patterns['o_locus_distribution'][o] += 1
            if capsule != 'ND':
                patterns['capsule_type_distribution'][capsule] += 1
            if pasteur != 'ND' and o != 'ND':
                patterns['st_o_combinations'][f"ST{pasteur} - {o}"].append(sample)
            if pasteur != 'ND' and k != 'ND':
                patterns['st_k_combinations'][f"ST{pasteur} - {k}"].append(sample)
            if k != 'ND' and o != 'ND':
                patterns['ko_combinations'][f"{k}:{o}"].append(sample)
            if pasteur != 'ND' and k != 'ND' and o != 'ND':
                patterns['st_ko_combinations'][f"ST{pasteur} - {k}:{o}"].append(sample)
        return patterns

    def build_sample_typing_map(self, samples_data: Dict[str, Dict]) -> Dict[str, Dict]:
        typing_map = {}
        for sample, data in samples_data.items():
            typing_map[sample] = {
                'pasteur_st': data.get('pasteur_mlst', {}).get('ST', 'ND'),
                'oxford_st': data.get('oxford_mlst', {}).get('ST', 'ND'),
                'k_locus': data.get('kaptive', {}).get('K_Locus', 'ND'),
                'o_locus': data.get('kaptive', {}).get('O_Locus', 'ND'),
                'capsule': data.get('kaptive', {}).get('Capsule_Type', 'ND')
            }
        return typing_map


class SampleCentricHTMLGenerator:
    def __init__(self, analyzer: SampleCentricAnalyzer):
        self.analyzer = analyzer
        self.tab_colors = {
            'summary': '#4CAF50',
            'samples': '#2196F3',
            'qc': '#607D8B',
            'mlst': '#FF9800',
            'kaptive': '#9C27B0',
            'combinations': '#009688',
            'amr': '#F44336',
            'virulence': '#E91E63',
            'bacmet': '#795548',
            'plasmids': '#673AB7',
            'mutation': '#00BCD4',
            'cooccurrence': '#3F51B5',
            'patterns': '#FF5722',
            'aiguide': '#3F51B5',
            'citation': '#8BC34A',
            'funding': '#FFC107',
            'export': '#9E9E9E',
            'calltoaction': '#2E7D32'
        }

    def generate_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        print("\n🎨 Generating SAMPLE‑CENTRIC HTML report...")
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        qc_data = integrated_data.get('qc_data', {})
        amr_details = integrated_data.get('amr_details', {})
        abricate_details = integrated_data.get('abricate_details', {})
        apt_data = integrated_data.get('apt_data', {})
        mutation_details = integrated_data.get('mutation_details', {})
        typing_map = self.analyzer.build_sample_typing_map(samples_data)

        html = self._create_html(
            metadata=metadata,
            samples_data=samples_data,
            patterns=patterns,
            gene_centric=gene_centric,
            qc_data=qc_data,
            total_samples=len(samples_data),
            amr_details=amr_details,
            abricate_details=abricate_details,
            apt_data=apt_data,
            mutation_details=mutation_details,
            typing_map=typing_map
        )
        output_file = output_dir / "genius_acinetobacter_sample_centric_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)

    def _create_html(self, **kwargs) -> str:
        metadata = kwargs['metadata']
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        qc_data = kwargs.get('qc_data', {})
        total_samples = kwargs['total_samples']
        amr_details = kwargs.get('amr_details', {})
        abricate_details = kwargs.get('abricate_details', {})
        apt_data = kwargs.get('apt_data', {})
        mutation_details = kwargs.get('mutation_details', {})
        typing_map = kwargs.get('typing_map', {})

        total_amr = sum(len(db) for db in gene_centric.get('amr_databases', {}).values())
        total_vir = sum(len(db) for db in gene_centric.get('virulence_databases', {}).values())
        total_bacmet = sum(len(db) for db in gene_centric.get('bacmet_databases', {}).values())
        total_plasmid = sum(len(db) for db in gene_centric.get('plasmid_databases', {}).values())
        high_risk = len(patterns.get('high_risk_crab', []))
        carbapenemase_count = sum(1 for db in gene_centric.get('amr_databases', {}).values() for g in db if g['category'] == 'Carbapenemases')
        env_count = total_bacmet
        mutation_count = len([s for s, muts in mutation_details.items() if muts])

        typing_json = json.dumps(typing_map)
        css = self._get_css()
        js = self._get_js(typing_json)

        html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>GENIUS Acinetobacter baumannii Sample‑Centric Report</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
{css}{js}</head>
<body><div class="container">
<div class="main-header"><h1><i class="fas fa-bacterium"></i> GENIUS Acinetobacter baumannii Sample‑Centric Analysis Report</h1>
<p>Interactive per‑isolate boxes for AMR, Virulence, BACMET, Plasmids, Mutations – with full TSV details</p>
<div class="metadata-bar"><div class="metadata-item"><i class="fas fa-calendar"></i> {metadata.get('analysis_date','Unknown')}</div>
<div class="metadata-item"><i class="fas fa-database"></i> Samples: {total_samples}</div>
<div class="metadata-item"><i class="fas fa-vial"></i> Pathogen: Acinetobacter baumannii</div>
<div class="metadata-item"><i class="fas fa-university"></i> University of Ghana Medical School</div></div></div>
<div class="dashboard-grid">
<div class="dashboard-card card-summary" onclick="switchTab('summary')"><div class="card-number">{total_samples}</div><div class="card-label">Total Samples</div><i class="fas fa-vial fa-2x" style="color:var(--summary-color)"></i></div>
<div class="dashboard-card card-mlst" onclick="switchTab('mlst')"><div class="card-number">{len(patterns.get('pasteur_st_distribution',{}))}</div><div class="card-label">Pasteur STs</div><i class="fas fa-code-branch fa-2x" style="color:var(--mlst-color)"></i></div>
<div class="dashboard-card card-kaptive" onclick="switchTab('kaptive')"><div class="card-number">{len(patterns.get('k_locus_distribution',{}))}</div><div class="card-label">Capsule Types</div><i class="fas fa-shield-alt fa-2x" style="color:var(--kaptive-color)"></i></div>
<div class="dashboard-card card-amr" onclick="switchTab('amr')"><div class="card-number">{total_amr}</div><div class="card-label">AMR Genes</div><i class="fas fa-biohazard fa-2x" style="color:var(--amr-color)"></i></div>
<div class="dashboard-card card-virulence" onclick="switchTab('virulence')"><div class="card-number">{total_vir}</div><div class="card-label">Virulence Genes</div><i class="fas fa-virus fa-2x" style="color:var(--virulence-color)"></i></div>
<div class="dashboard-card card-bacmet" onclick="switchTab('bacmet')"><div class="card-number">{total_bacmet}</div><div class="card-label">Bacmet</div><i class="fas fa-flask fa-2x" style="color:var(--bacmet-color)"></i></div>
<div class="dashboard-card card-plasmids" onclick="switchTab('plasmids')"><div class="card-number">{total_plasmid}</div><div class="card-label">Plasmid Replicons</div><i class="fas fa-dna fa-2x" style="color:var(--plasmids-color)"></i></div>
<div class="dashboard-card card-patterns" onclick="switchTab('patterns')"><div class="card-number">{high_risk}</div><div class="card-label">CRAB High-Risk</div><i class="fas fa-exclamation-triangle fa-2x" style="color:var(--patterns-color)"></i></div>
</div>
<div class="tab-navigation">
<button class="tab-button summary active" onclick="switchTab('summary')"><i class="fas fa-chart-pie"></i> Summary</button>
<button class="tab-button samples" onclick="switchTab('samples')"><i class="fas fa-list-alt"></i> Sample Overview</button>
<button class="tab-button qc" onclick="switchTab('qc')"><i class="fas fa-chart-line"></i> FASTA QC</button>
<button class="tab-button mlst" onclick="switchTab('mlst')"><i class="fas fa-code-branch"></i> MLST</button>
<button class="tab-button kaptive" onclick="switchTab('kaptive')"><i class="fas fa-shield-alt"></i> Capsule</button>
<button class="tab-button combinations" onclick="switchTab('combinations')"><i class="fas fa-link"></i> Combinations</button>
<button class="tab-button amr" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> AMR</button>
<button class="tab-button virulence" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> Virulence</button>
<button class="tab-button bacmet" onclick="switchTab('bacmet')"><i class="fas fa-flask"></i> Bacmet</button>
<button class="tab-button plasmids" onclick="switchTab('plasmids')"><i class="fas fa-dna"></i> Plasmids</button>
<button class="tab-button mutation" onclick="switchTab('mutation')"><i class="fas fa-dna"></i> Mutations</button>
<button class="tab-button cooccurrence" onclick="switchTab('cooccurrence')"><i class="fas fa-link"></i> Co‑occurrence</button>
<button class="tab-button patterns" onclick="switchTab('patterns')"><i class="fas fa-project-diagram"></i> Patterns</button>
<button class="tab-button aiguide" onclick="switchTab('aiguide')"><i class="fas fa-robot"></i> AI Guide</button>
<button class="tab-button citation" onclick="switchTab('citation')"><i class="fas fa-book"></i> Citations</button>
<button class="tab-button funding" onclick="switchTab('funding')"><i class="fas fa-coffee"></i> Funding</button>
<button class="tab-button export" onclick="switchTab('export')"><i class="fas fa-download"></i> Export</button>
</div>
<div id="summary-tab" class="tab-content active">{self._summary_section(kwargs, carbapenemase_count, env_count, mutation_count)}</div>
<div id="samples-tab" class="tab-content">{self._samples_section(kwargs)}</div>
<div id="qc-tab" class="tab-content">{self._qc_section(kwargs)}</div>
<div id="mlst-tab" class="tab-content">{self._mlst_section(kwargs)}</div>
<div id="kaptive-tab" class="tab-content">{self._kaptive_section(kwargs)}</div>
<div id="combinations-tab" class="tab-content">{self._combinations_section(kwargs)}</div>
<div id="amr-tab" class="tab-content">{self._amr_boxes(kwargs)}</div>
<div id="virulence-tab" class="tab-content">{self._virulence_boxes(kwargs)}</div>
<div id="bacmet-tab" class="tab-content">{self._bacmet_boxes(kwargs)}</div>
<div id="plasmids-tab" class="tab-content">{self._plasmids_boxes(kwargs)}</div>
<div id="mutation-tab" class="tab-content">{self._mutation_boxes(kwargs)}</div>
<div id="cooccurrence-tab" class="tab-content">{self._cooccurrence_section(kwargs)}</div>
<div id="patterns-tab" class="tab-content">{self._patterns_section(kwargs)}</div>
<div id="aiguide-tab" class="tab-content">{self._aiguide_section()}</div>
<div id="citation-tab" class="tab-content">{self._citation_section()}</div>
<div id="funding-tab" class="tab-content">{self._funding_section()}</div>
<div id="export-tab" class="tab-content">{self._export_section()}</div>
<div class="footer"><h3>GENIUS Acinetobacter baumannii Sample‑Centric Reporter v1.1.0</h3><p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p><p>If you find this tool helpful, please <a href="https://github.com/bbeckley-hub/acinetoscope" target="_blank">⭐ star us on GitHub</a> and share with your network.</p><p>Generated on {metadata.get('analysis_date','Unknown')}</p></div>
</div></body></html>"""
        return html

    def _get_css(self) -> str:
        return """
        <style>
        :root { --summary-color:#4CAF50; --samples-color:#2196F3; --mlst-color:#FF9800; --kaptive-color:#9C27B0; --amr-color:#F44336; --virulence-color:#E91E63; --bacmet-color:#795548; --combinations-color:#009688; --patterns-color:#FF5722; --plasmids-color:#673AB7; --databases-color:#607D8B; --qc-color:#17a2b8; --aiguide-color:#3F51B5; --calltoaction-color:#2E7D32; --export-color:#3F51B5; --mutation-color:#00BCD4; --cooccurrence-color:#3F51B5; --citation-color:#8BC34A; --funding-color:#FFC107; }
        *{margin:0;padding:0;box-sizing:border-box}
        body{font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;line-height:1.6;color:#333;background:#f5f5f5;min-width:1200px}
        .container{max-width:none;margin:0 auto;padding:20px;width:100%;overflow-x:auto}
        .main-header{background:linear-gradient(135deg,#00695c 0%,#004d40 100%);color:white;padding:30px;border-radius:15px;box-shadow:0 10px 30px rgba(0,0,0,0.2);margin-bottom:30px;text-align:center}
        .main-header h1{font-size:2.8em;margin-bottom:10px;color:white}
        .metadata-bar{background:rgba(255,255,255,0.1);padding:15px;border-radius:10px;margin:20px 0;display:flex;justify-content:space-around;flex-wrap:wrap;gap:15px}
        .metadata-item{display:flex;align-items:center;gap:8px;font-size:0.95em}
        .dashboard-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:20px;margin-bottom:30px}
        .dashboard-card{background:white;padding:25px;border-radius:12px;box-shadow:0 5px 20px rgba(0,0,0,0.1);text-align:center;cursor:pointer;border-left:5px solid;transition:all 0.3s ease}
        .dashboard-card:hover{transform:translateY(-10px);box-shadow:0 15px 30px rgba(0,0,0,0.2)}
        .card-number{font-size:3em;font-weight:bold;margin:15px 0;background:linear-gradient(90deg,#00695c,#004d40);-webkit-background-clip:text;-webkit-text-fill-color:transparent}
        .tab-navigation{display:flex;gap:5px;margin-bottom:20px;flex-wrap:wrap;background:white;padding:15px;border-radius:12px;box-shadow:0 5px 20px rgba(0,0,0,0.1);position:sticky;top:10px;z-index:100}
        .tab-button{padding:12px 20px;background:#f5f5f5;border:none;border-radius:8px;cursor:pointer;font-weight:600;color:#666;transition:all 0.3s ease;display:flex;align-items:center;gap:8px;font-size:0.9em}
        .tab-button.active{color:white}
        .tab-button.summary.active{background:var(--summary-color)}
        .tab-button.samples.active{background:var(--samples-color)}
        .tab-button.qc.active{background:var(--qc-color)}
        .tab-button.mlst.active{background:var(--mlst-color)}
        .tab-button.kaptive.active{background:var(--kaptive-color)}
        .tab-button.combinations.active{background:var(--combinations-color)}
        .tab-button.amr.active{background:var(--amr-color)}
        .tab-button.virulence.active{background:var(--virulence-color)}
        .tab-button.bacmet.active{background:var(--bacmet-color)}
        .tab-button.plasmids.active{background:var(--plasmids-color)}
        .tab-button.mutation.active{background:var(--mutation-color)}
        .tab-button.cooccurrence.active{background:var(--cooccurrence-color)}
        .tab-button.patterns.active{background:var(--patterns-color)}
        .tab-button.aiguide.active{background:var(--aiguide-color)}
        .tab-button.citation.active{background:var(--citation-color)}
        .tab-button.funding.active{background:var(--funding-color)}
        .tab-button.export.active{background:var(--export-color)}
        .tab-content{display:none;background:white;padding:30px;border-radius:15px;box-shadow:0 10px 30px rgba(0,0,0,0.1);margin-bottom:30px;animation:fadeIn 0.5s ease;width:100%;overflow-x:auto}
        .tab-content.active{display:block}
        @keyframes fadeIn{from{opacity:0;transform:translateY(20px)}to{opacity:1;transform:translateY(0)}}
        .section-header{color:#2c3e50;margin-bottom:25px;padding-bottom:15px;border-bottom:3px solid;font-size:1.8em;display:flex;align-items:center;justify-content:space-between}
        .print-section-btn{background:#00695c;color:white;border:none;border-radius:5px;padding:8px 15px;cursor:pointer;display:flex;align-items:center;gap:5px;font-size:0.9em}
        .print-section-btn:hover{background:#004d40}
        .data-table{width:100%;border-collapse:collapse;margin:20px 0;font-size:0.95em;box-shadow:0 2px 10px rgba(0,0,0,0.1);border-radius:8px;overflow:hidden;table-layout:auto}
        .data-table th{background:#2c3e50;color:white;padding:15px;text-align:left;font-weight:600;white-space:nowrap;cursor:pointer}
        .data-table th:hover{background:#1a252f}
        .data-table td{padding:12px;border-bottom:1px solid #e0e0e0;vertical-align:top;white-space:nowrap}
        .data-table tr:hover{background:#f8f9fa}
        .master-scrollable-container{width:100%;overflow-x:auto;border:1px solid #e0e0e0;border-radius:8px;margin:20px 0}
        .genome-list{display:flex;flex-wrap:wrap;gap:5px;max-height:200px;overflow-y:auto;padding:5px;background:#f8f9fa;border-radius:5px}
        .genome-group{margin-bottom:10px;width:100%}
        .genome-group-header{font-weight:bold;background:#e0e0e0;padding:4px 8px;border-radius:4px;margin:5px 0;font-size:0.85em;display:inline-block}
        .genome-group-tags{display:flex;flex-wrap:wrap;gap:5px;margin-left:10px}
        .genome-tag{display:inline-block;background:#e6ffe6;color:#006400;padding:3px 10px;border-radius:12px;font-size:0.85em;border:1px solid #b3ffb3;white-space:nowrap;margin:2px}
        .genome-tag.highlight{background-color:#ffff99 !important;color:#000 !important;border:1px solid #ffc107}
        .search-box{width:100%;padding:12px;margin-bottom:15px;border:2px solid #e0e0e0;border-radius:8px;font-size:1em;transition:all 0.3s ease}
        .search-box:focus{outline:none;border-color:#00695c;box-shadow:0 0 0 3px rgba(0,105,92,0.1)}
        .badge{display:inline-block;padding:5px 15px;border-radius:20px;font-size:0.85em;font-weight:600;margin:2px}
        .badge-critical{background:#9C27B0;color:white}
        .badge-high{background:#F44336;color:white}
        .badge-medium{background:#FF9800;color:black}
        .alert-box{padding:20px;border-radius:10px;margin:20px 0;display:flex;align-items:center;gap:20px;border-left:5px solid}
        .alert-success{background:#d4edda;color:#155724;border-left-color:#28a745}
        .alert-warning{background:#fff3cd;color:#856404;border-left-color:#ffc107}
        .alert-danger{background:#f8d7da;color:#721c24;border-left-color:#dc3545}
        .alert-info{background:#d1ecf1;color:#0c5460;border-left-color:#17a2b8}
        .action-buttons{display:flex;gap:10px;margin:20px 0;flex-wrap:wrap}
        .action-btn{padding:10px 20px;border:none;border-radius:8px;cursor:pointer;font-weight:600;display:inline-flex;align-items:center;gap:8px;transition:0.3s}
        .action-btn:hover{transform:translateY(-2px);box-shadow:0 5px 15px rgba(0,0,0,0.2)}
        .btn-primary{background:#00695c;color:white}
        .btn-success{background:#28a745;color:white}
        .btn-warning{background:#ffc107;color:black}
        .btn-danger{background:#dc3545;color:white}
        .btn-info{background:#17a2b8;color:white}
        .btn-secondary{background:#6c757d;color:white}
        .btn-light{background:#f8f9fa;color:#212529;border:1px solid #dee2e6}
        .grouping-controls{background:#f0f7f0;padding:12px;border-radius:8px;margin:15px 0;display:flex;flex-wrap:wrap;gap:8px;align-items:center;border-left:4px solid #00695c}
        .grouping-controls label{font-weight:bold;margin-right:5px}
        .group-btn{background:white;border:1px solid #00695c;color:#00695c;padding:6px 12px;border-radius:20px;cursor:pointer;font-size:0.85em;transition:all 0.2s}
        .group-btn:hover{background:#00695c;color:white}
        .group-btn.active{background:#00695c;color:white}
        .footer{text-align:center;padding:30px;background:linear-gradient(135deg,#2c3e50,#34495e);color:white;border-radius:15px;margin-top:40px}
        .footer a{color:#ffc107;text-decoration:none}
        .footer a:hover{text-decoration:underline}
        .database-section{margin:30px 0;padding:25px;border-radius:12px;background:#f8f9fa;box-shadow:0 3px 15px rgba(0,0,0,0.08)}
        .database-header{font-size:1.4em;color:#2c3e50;margin-bottom:20px;padding-bottom:10px;border-bottom:2px solid #00695c;display:flex;align-items:center;justify-content:space-between}
        .citation-card{background:linear-gradient(135deg, #f8f9fa, #e9ecef);padding:15px 20px;border-radius:10px;margin-bottom:12px;border-left:6px solid;box-shadow:0 2px 8px rgba(0,0,0,0.06);transition:0.2s}
        .citation-card:hover{box-shadow:0 4px 12px rgba(0,0,0,0.12)}
        .citation-card .cite-title{font-weight:700;color:#2c3e50}
        .citation-card .cite-text{color:#333;font-size:0.95em}
        .citation-card .cite-link{color:#00695c;text-decoration:underline;font-weight:500;margin-left:8px}
        .citation-card .cite-link:hover{color:#004d40}
        .citation-card .copy-btn{margin-left:12px;background:#00695c;color:white;border:none;padding:4px 12px;border-radius:20px;cursor:pointer;font-size:0.8em}
        .citation-card .copy-btn:hover{background:#004d40}
        .citation-card.cite-main{border-left-color:#00695c;background:linear-gradient(135deg, #e0f2f1, #b2dfdb)}
        .citation-card.cite-db{border-left-color:#9C27B0;background:linear-gradient(135deg, #f3e5f5, #e1bee7)}
        .citation-card.cite-tool{border-left-color:#FF9800;background:linear-gradient(135deg, #fff3e0, #ffe0b2)}
        .citation-card.cite-other{border-left-color:#607D8B;background:linear-gradient(135deg, #eceff1, #cfd8dc)}
        .citation-card .cite-badge{display:inline-block;font-size:0.7em;font-weight:600;padding:2px 10px;border-radius:12px;background:#00695c;color:white;margin-right:10px;text-transform:uppercase}
        .citation-card.cite-main .cite-badge{background:#00695c}
        .citation-card.cite-db .cite-badge{background:#9C27B0}
        .citation-card.cite-tool .cite-badge{background:#FF9800}
        .citation-card.cite-other .cite-badge{background:#607D8B}
        .isolate-box {
            border: 1px solid #ddd;
            border-radius: 12px;
            margin-bottom: 30px;
            padding: 20px;
            background: #fafafa;
            box-shadow: 0 2px 8px rgba(0,0,0,0.06);
            border-bottom: 2px dashed #ccc;
        }
        .isolate-box .sample-header { display: flex; align-items: center; gap: 15px; flex-wrap: wrap; margin-bottom: 15px; border-bottom: 2px solid #e0e0e0; padding-bottom: 10px; }
        .isolate-box .sample-header h3 { font-size: 1.4em; margin: 0; }
        .isolate-box .sample-header .total-badge { background: #00695c; color: white; padding: 4px 16px; border-radius: 20px; font-weight: bold; font-size: 0.9em; }
        .isolate-box .sample-header .typing-info { display: flex; flex-wrap: wrap; gap: 8px; margin-left: auto; }
        .isolate-box .sample-header .typing-badge { display: inline-block; padding: 3px 10px; border-radius: 12px; font-size: 0.8em; font-weight: 600; background: #e0e0e0; color: #333; border: 1px solid #ccc; }
        .isolate-box .sample-header .typing-badge.capsule-badge { background: #9C27B0; color: white; border-color: #9C27B0; }
        .isolate-box .sample-header .typing-badge.st-badge { background: #FF9800; color: white; border-color: #FF9800; }
        .isolate-box .sample-header .typing-badge.k-badge { background: #4CAF50; color: white; border-color: #4CAF50; }
        .isolate-box .sample-header .typing-badge.ocl-badge { background: #2196F3; color: white; border-color: #2196F3; }
        .isolate-box .sample-header .typing-badge.oxford-badge { background: #795548; color: white; border-color: #795548; }
        .isolate-box .database-table-wrapper { margin: 15px 0; overflow-x: auto; }
        .isolate-box .database-table-wrapper table { width: 100%; border-collapse: collapse; font-size: 0.85em; min-width: 800px; }
        .isolate-box .database-table-wrapper table th { background: #2c3e50; color: white; padding: 8px 12px; text-align: left; white-space: nowrap; }
        .isolate-box .database-table-wrapper table td { padding: 8px 12px; border-bottom: 1px solid #e0e0e0; white-space: nowrap; }
        .isolate-box .database-table-wrapper table tr:hover { background: #f1f1f1; }
        .isolate-box .db-title { font-weight: bold; color: #00695c; margin: 10px 0 5px 0; font-size: 1.1em; border-left: 4px solid #00695c; padding-left: 10px; }
        .isolate-box .apt-rep-list { background: #f8f9fa; padding: 10px 15px; border-radius: 6px; margin: 5px 0 10px 0; font-size: 0.95em; }
        .filter-controls { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; }
        .filter-controls select { padding: 10px; border-radius: 8px; border: 2px solid #ddd; background: white; min-width: 150px; }
        @media print{body *{visibility:hidden}.tab-content.active,.tab-content.active *{visibility:visible}.tab-content.active{position:absolute;left:0;top:0;width:100%;padding:20px;box-shadow:none;border-radius:0}.print-section-btn,.tab-navigation,.dashboard-grid,.search-box,.action-buttons,.grouping-controls,.filter-controls{display:none !important}.data-table{page-break-inside:auto}.data-table tr{page-break-inside:avoid;page-break-after:auto}}
        @media (max-width:768px){.container{padding:10px}.main-header{padding:20px}.main-header h1{font-size:2em}.tab-button{padding:8px 12px;font-size:0.8em}.dashboard-grid{grid-template-columns:repeat(auto-fit,minmax(180px,1fr))}.data-table{font-size:0.8em}body{min-width:auto;overflow-x:auto}}
        </style>
        """

    def _get_js(self, typing_json: str) -> str:
        return f"""
        <script>
        var sampleTyping = {typing_json};
        var originalGenomeLists = {{}};

        function switchTab(tabName) {{
            document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(b => b.classList.remove('active'));
            document.getElementById(tabName + '-tab').classList.add('active');
            event.currentTarget.classList.add('active');
            window.location.hash = tabName;
        }}

        function searchTable(tableId, searchId) {{
            var filter = document.getElementById(searchId).value.toUpperCase();
            var rows = document.getElementById(tableId).getElementsByTagName('tr');
            for (var i = 1; i < rows.length; i++) {{
                var found = false;
                for (var cell of rows[i].getElementsByTagName('td')) {{
                    if (cell.textContent.toUpperCase().indexOf(filter) > -1) {{ found = true; break; }}
                }}
                rows[i].style.display = found ? '' : 'none';
            }}
        }}

        function highlightGenome(tableId, searchId) {{
            var filter = document.getElementById(searchId).value.toUpperCase().trim();
            var table = document.getElementById(tableId);
            var allTags = table.querySelectorAll('.genome-tag');
            allTags.forEach(tag => tag.classList.remove('highlight'));
            if (filter === '') return;
            allTags.forEach(tag => {{
                if (tag.textContent.toUpperCase().indexOf(filter) > -1) {{
                    tag.classList.add('highlight');
                }}
            }});
        }}

        function getTypingValue(genome, groupBy) {{
            var info = sampleTyping[genome];
            if (!info) return "Unknown";
            if (groupBy === "pasteur_st") {{
                var st = info.pasteur_st;
                return st !== 'ND' ? "ST " + st : "ND";
            }}
            if (groupBy === "oxford_st") {{
                var st = info.oxford_st;
                return st !== 'ND' ? "ST " + st : "ND";
            }}
            if (groupBy === "k_locus") return info.k_locus;
            if (groupBy === "o_locus") return info.o_locus;
            if (groupBy === "capsule") return info.capsule;
            if (groupBy === "st_k") {{
                var st = info.pasteur_st;
                var st_pre = st !== 'ND' ? "ST " + st : "ND";
                return st_pre + " - " + info.k_locus;
            }}
            if (groupBy === "st_o") {{
                var st = info.pasteur_st;
                var st_pre = st !== 'ND' ? "ST " + st : "ND";
                return st_pre + " - " + info.o_locus;
            }}
            if (groupBy === "k_o") return info.k_locus + ":" + info.o_locus;
            if (groupBy === "st_ko") {{
                var st = info.pasteur_st;
                var st_pre = st !== 'ND' ? "ST " + st : "ND";
                return st_pre + " - " + info.k_locus + ":" + info.o_locus;
            }}
            return "Unknown";
        }}

        function groupRowGenomes(row, groupBy, originalList) {{
            var genomesCell = null;
            for (var i = 0; i < row.cells.length; i++) {{
                if (row.cells[i].querySelector('.genome-list')) {{
                    genomesCell = row.cells[i];
                    break;
                }}
            }}
            if (!genomesCell) return;
            var genomes = originalList.slice();
            if (genomes.length === 0) {{
                genomesCell.innerHTML = '<div class="genome-list">None</div>';
                return;
            }}
            var groups = {{}};
            genomes.forEach(function(genome) {{
                var key = getTypingValue(genome, groupBy);
                if (!groups[key]) groups[key] = [];
                groups[key].push(genome);
            }});
            var html = '<div class="genome-list">';
            for (var key in groups) {{
                var tags = groups[key].map(g => `<span class="genome-tag">${{g}}</span>`).join('');
                html += `<div class="genome-group"><div class="genome-group-header">${{key}}</div><div class="genome-group-tags">${{tags}}</div></div>`;
            }}
            html += '</div>';
            genomesCell.innerHTML = html;
        }}

        function groupGenomesByTyping(tableId, groupBy) {{
            var table = document.getElementById(tableId);
            if (!table) return;
            var tbody = table.tBodies[0];
            if (!tbody) return;
            var rows = tbody.rows;
            for (var i = 0; i < rows.length; i++) {{
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                if (!originalGenomeLists[geneName]) {{
                    var genomesCell = null;
                    for (var j = 0; j < row.cells.length; j++) {{
                        if (row.cells[j].querySelector('.genome-list')) {{
                            genomesCell = row.cells[j];
                            break;
                        }}
                    }}
                    if (genomesCell) {{
                        var tags = genomesCell.querySelectorAll('.genome-tag');
                        var genomes = Array.from(tags).map(tag => tag.textContent.trim());
                        originalGenomeLists[geneName] = genomes;
                    }} else {{
                        originalGenomeLists[geneName] = [];
                    }}
                }}
            }}
            for (var i = 0; i < rows.length; i++) {{
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                var original = originalGenomeLists[geneName] || [];
                groupRowGenomes(row, groupBy, original);
            }}
            var container = table.closest('.tab-content');
            if (container) {{
                var btns = container.querySelectorAll('.group-btn');
                btns.forEach(btn => btn.classList.remove('active'));
                var activeBtn = container.querySelector(`.group-btn[data-group="${{groupBy}}"]`);
                if (activeBtn) activeBtn.classList.add('active');
            }}
        }}

        function resetGenomeList(tableId) {{
            var table = document.getElementById(tableId);
            if (!table) return;
            var tbody = table.tBodies[0];
            if (!tbody) return;
            var rows = tbody.rows;
            for (var i = 0; i < rows.length; i++) {{
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                var original = originalGenomeLists[geneName] || [];
                var genomesCell = null;
                for (var j = 0; j < row.cells.length; j++) {{
                    if (row.cells[j].querySelector('.genome-list')) {{
                        genomesCell = row.cells[j];
                        break;
                    }}
                }}
                if (genomesCell) {{
                    var tags = original.map(g => `<span class="genome-tag">${{g}}</span>`).join('');
                    genomesCell.innerHTML = `<div class="genome-list">${{tags}}</div>`;
                }}
            }}
            var container = table.closest('.tab-content');
            if (container) {{
                var btns = container.querySelectorAll('.group-btn');
                btns.forEach(btn => btn.classList.remove('active'));
            }}
        }}

        function filterBoxes(tabId) {{
            var search = document.getElementById('search-' + tabId).value.toUpperCase();
            var dbFilter = document.getElementById('dbFilter-' + tabId).value;
            var boxes = document.querySelectorAll('#' + tabId + '-tab .isolate-box');
            boxes.forEach(function(box) {{
                var sample = box.getAttribute('data-sample') || '';
                var show = true;
                if (search && sample.toUpperCase().indexOf(search) === -1) show = false;
                var dbWrappers = box.querySelectorAll('.database-table-wrapper');
                if (dbFilter !== 'all') {{
                    dbWrappers.forEach(function(wrapper) {{
                        var dbName = wrapper.getAttribute('data-db') || '';
                        if (dbName === dbFilter) {{
                            wrapper.style.display = '';
                        }} else {{
                            wrapper.style.display = 'none';
                        }}
                    }});
                    var dbTitles = box.querySelectorAll('.db-title');
                    dbTitles.forEach(function(title) {{
                        var wrapper = title.nextElementSibling;
                        if (wrapper && wrapper.style.display === 'none') {{
                            title.style.display = 'none';
                        }} else {{
                            title.style.display = '';
                        }}
                    }});
                }} else {{
                    dbWrappers.forEach(function(wrapper) {{ wrapper.style.display = ''; }});
                    box.querySelectorAll('.db-title').forEach(function(title) {{ title.style.display = ''; }});
                }}
                var aptSection = box.querySelector('.apt-section');
                if (aptSection) {{
                    if (dbFilter === 'APT' || dbFilter === 'all') {{
                        aptSection.style.display = '';
                    }} else {{
                        aptSection.style.display = 'none';
                    }}
                }}
                box.style.display = show ? '' : 'none';
            }});
        }}

        function resetBoxFilters(tabId) {{
            document.getElementById('search-' + tabId).value = '';
            document.getElementById('dbFilter-' + tabId).value = 'all';
            filterBoxes(tabId);
        }}

        function sortTable(tableId, colIndex, type) {{
            var table = document.getElementById(tableId);
            var tbody = table.tBodies[0];
            var rows = Array.from(tbody.rows);
            var isAscending = table.getAttribute('data-sort-dir') !== 'asc';
            rows.sort((a,b) => {{
                var aVal = a.cells[colIndex].innerText.trim();
                var bVal = b.cells[colIndex].innerText.trim();
                if (type === 'number') {{
                    aVal = parseFloat(aVal.replace(/,/g,'')) || 0;
                    bVal = parseFloat(bVal.replace(/,/g,'')) || 0;
                    return isAscending ? aVal - bVal : bVal - aVal;
                }} else {{
                    return isAscending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                }}
            }});
            tbody.append(...rows);
            table.setAttribute('data-sort-dir', isAscending ? 'asc' : 'desc');
            var headers = table.querySelectorAll('th');
            headers.forEach((th,idx) => {{
                var icon = th.querySelector('.sort-icon');
                if (icon) icon.innerHTML = '⇅';
            }});
            var currentHeader = headers[colIndex];
            var icon = currentHeader.querySelector('.sort-icon');
            if (icon) icon.innerHTML = isAscending ? '↑' : '↓';
        }}

        function printSection(sectionId) {{
            var content = document.getElementById(sectionId);
            var win = window.open('', '_blank');
            win.document.write('<html><head><title>Print</title><style>'+document.querySelector('style').textContent+'</style></head><body>'+content.innerHTML+'</body></html>');
            win.document.close();
            win.print();
        }}

        function exportTableToCSV(tableId, filename) {{
            var table = document.getElementById(tableId);
            var rows = table.querySelectorAll('tr');
            var csv = [];
            for (var i = 0; i < rows.length; i++) {{
                var row = [], cols = rows[i].querySelectorAll('td, th');
                for (var j = 0; j < cols.length; j++) {{
                    row.push('"' + (cols[j].innerText || '').replace(/"/g, '""') + '"');
                }}
                csv.push(row.join(','));
            }}
            var blob = new Blob([csv.join('\\n')], {{type:'text/csv'}});
            var a = document.createElement('a');
            a.href = URL.createObjectURL(blob);
            a.download = filename;
            a.click();
            URL.revokeObjectURL(a.href);
        }}

        document.addEventListener('DOMContentLoaded', function() {{
            var hash = window.location.hash.substring(1);
            if (hash) {{
                var btn = document.querySelector(`.tab-button.${{hash}}`);
                if (btn) btn.click();
            }} else {{
                document.querySelector('.tab-button').click();
            }}
            document.querySelectorAll('.data-table').forEach(table => {{
                var headers = table.querySelectorAll('th');
                headers.forEach((th, idx) => {{
                    var type = th.getAttribute('data-sort') || 'string';
                    th.style.cursor = 'pointer';
                    th.addEventListener('click', () => sortTable(table.id, idx, type));
                    var icon = document.createElement('span');
                    icon.className = 'sort-icon';
                    icon.innerHTML = '⇅';
                    th.appendChild(icon);
                }});
            }});
        }});
        </script>
        """

    def _make_genome_tags(self, genomes: List[str]) -> str:
        return ''.join(f'<span class="genome-tag">{gen}</span>' for gen in sorted(genomes))

    def _make_database_table(self, db_name: str, genes: List[Dict], db_key: str) -> str:
        if not genes:
            return ''
        keys = list(genes[0].keys()) if genes else []
        priority = ['gene', 'product', 'coverage_percent', 'identity_percent', 'accession', 'contig', 'start', 'stop', 'class', 'subclass', 'scope', 'resistance']
        ordered = []
        for p in priority:
            if p in keys:
                ordered.append(p)
        for k in keys:
            if k not in ordered:
                ordered.append(k)
        html = f"""
                <div class="database-table-wrapper" data-db="{db_key}">
                    <div class="db-title">{db_name}</div>
                    <table>
                        <thead>
                            <tr>
        """
        for col in ordered:
            display = col.replace('_', ' ').title()
            html += f'<th>{display}</th>'
        html += '</tr></thead><tbody>'
        for gene_dict in genes:
            html += '<tr>'
            for col in ordered:
                val = gene_dict.get(col, '')
                if val is None:
                    val = ''
                html += f'<td>{val}</td>'
            html += '</tr>'
        html += """
                        </tbody>
                    </table>
                </div>
        """
        return html

    def _make_mutation_table(self, mutations: List[Dict]) -> str:
        if not mutations:
            return ''
        cols = ['gene', 'mutation', 'class', 'subclass', 'contig', 'start', 'stop', 'strand', 'coverage', 'identity', 'accession']
        html = """
        <div class="database-table-wrapper" data-db="mutations">
            <div class="db-title">Mutations</div>
            <table>
                <thead>
                    <tr>
        """
        for col in cols:
            display = col.replace('_', ' ').title()
            html += f'<th>{display}</th>'
        html += '</tr></thead><tbody>'
        for mut in mutations:
            html += '<tr>'
            for col in cols:
                val = mut.get(col, '')
                if val is None or (isinstance(val, float) and val != val):
                    val = ''
                html += f'<td>{val}</td>'
            html += '</tr>'
        html += """
                </tbody>
            </table>
        </div>
        """
        return html

    def _build_isolate_box(self, sample: str, kwargs: Dict,
                           amr_details: Dict = None,
                           abricate_details: Dict = None,
                           apt_data: Dict = None,
                           db_list: List[str] = None,
                           include_apt: bool = False,
                           mutation_details: Dict = None) -> str:
        samples_data = kwargs.get('samples_data', {})
        data = samples_data.get(sample, {})
        pasteur = data.get('pasteur_mlst', {}).get('ST', 'ND')
        oxford = data.get('oxford_mlst', {}).get('ST', 'ND')
        k = data.get('kaptive', {}).get('K_Locus', 'ND')
        o = data.get('kaptive', {}).get('O_Locus', 'ND')
        capsule = data.get('kaptive', {}).get('Capsule_Type', 'ND')

        total_items = 0
        if amr_details and sample in amr_details and 'amrfinder' in db_list:
            total_items += len(amr_details[sample])
        if abricate_details and sample in abricate_details:
            for db, genes in abricate_details[sample].items():
                if db in db_list:
                    total_items += len(genes)
        if apt_data and include_apt and sample in apt_data:
            total_items += len(apt_data[sample].get('rep_types', []))
        if mutation_details and sample in mutation_details:
            total_items += len(mutation_details[sample])

        html = f'<div class="isolate-box" data-sample="{sample}">'
        html += f"""
            <div class="sample-header">
                <h3><i class="fas fa-microbe"></i> {sample}</h3>
                <span class="total-badge">Total Items: {total_items}</span>
                <div class="typing-info">
                    <span class="typing-badge st-badge">Pasteur: ST {pasteur}</span>
                    <span class="typing-badge oxford-badge">Oxford: ST {oxford}</span>
                    <span class="typing-badge k-badge">K: {k}</span>
                    <span class="typing-badge ocl-badge">OCL: {o}</span>
                    <span class="typing-badge capsule-badge">Capsule: {capsule}</span>
                </div>
            </div>
        """

        if 'amrfinder' in db_list and amr_details and sample in amr_details and amr_details[sample]:
            html += self._make_database_table('AMRfinder', amr_details[sample], 'amrfinder')

        if abricate_details and sample in abricate_details:
            for db, genes in abricate_details[sample].items():
                if db in db_list:
                    html += self._make_database_table(db.upper(), genes, db)

        if include_apt and apt_data and sample in apt_data:
            apt_info = apt_data[sample]
            rep_types = apt_info.get('rep_types', [])
            hit_count = apt_info.get('hit_count', 0)
            if rep_types:
                html += f"""
                <div class="apt-section" data-db="APT">
                    <div class="db-title">APT Rep Types</div>
                    <div class="apt-rep-list">
                        <strong>Hit Count:</strong> {hit_count} &nbsp;|&nbsp;
                        <strong>Rep Types:</strong> {', '.join(rep_types)}
                    </div>
                </div>
                """

        if mutation_details and sample in mutation_details:
            html += self._make_mutation_table(mutation_details[sample])

        html += '</div>'
        return html

    def _sample_boxes_generic(self, kwargs, tab_id: str, title: str,
                              db_list: List[str], credit_type: str = '',
                              include_apt: bool = False) -> str:
        samples_data = kwargs.get('samples_data', {})
        amr_details = kwargs.get('amr_details', {})
        abricate_details = kwargs.get('abricate_details', {})
        apt_data = kwargs.get('apt_data', {})

        relevant_samples = []
        for sample in samples_data:
            has_data = False
            if 'amrfinder' in db_list and sample in amr_details and amr_details[sample]:
                has_data = True
            for db in db_list:
                if db != 'amrfinder' and sample in abricate_details and db in abricate_details[sample] and abricate_details[sample][db]:
                    has_data = True
            if include_apt and sample in apt_data and apt_data[sample].get('rep_types'):
                has_data = True
            if has_data:
                relevant_samples.append(sample)
        relevant_samples.sort()

        if not relevant_samples:
            return f"""
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No {title} Data Available</h3><p>No samples with {title} genes were found.</p></div>
            </div>
            """

        credit_html = self._get_credit_bar(credit_type)

        html = f"""
        {credit_html}
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>🧬 {title} – Interactive Isolate Boxes</h3>
                <p>Each box represents one isolate. Inside, you will find separate tables for each database with full gene details – horizontally scrollable.</p>
                <p>Use the filters below to search by sample name or to show only a specific database.</p>
            </div>
        </div>
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-{tab_id}" onkeyup="filterBoxes('{tab_id}')" placeholder="🔍 Search sample...">
            <select id="dbFilter-{tab_id}" onchange="filterBoxes('{tab_id}')">
                <option value="all">All Databases</option>
        """
        if 'amrfinder' in db_list:
            html += '<option value="amrfinder">AMRfinder</option>'
        for db in db_list:
            if db != 'amrfinder':
                html += f'<option value="{db}">{db.upper()}</option>'
        if include_apt:
            html += '<option value="APT">APT</option>'
        html += f"""
            </select>
            <button class="action-btn btn-success" onclick="resetBoxFilters('{tab_id}')"><i class="fas fa-sync"></i> Clear Filters</button>
        </div>
        <div id="box-container-{tab_id}">
        """

        for sample in relevant_samples:
            html += self._build_isolate_box(sample, kwargs, amr_details=amr_details,
                                            abricate_details=abricate_details,
                                            apt_data=apt_data,
                                            db_list=db_list,
                                            include_apt=include_apt)
        html += """
        </div>
        """
        return html

    def _get_credit_bar(self, credit_type: str) -> str:
        if credit_type == 'amr':
            return """
            <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #f8d7da 100%); border-left: 6px solid #dc3545; margin-bottom: 20px;">
                <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                    <span style="font-size: 1.4em;">🧬</span>
                    <div>
                        <strong style="font-size: 1.2em; color: #a71d2a;">AMR Detection Tools & Databases</strong><br>
                        <span style="font-size: 0.95em; color: #333;">
                            <strong>ABRicate</strong> by <a href="https://github.com/tseemann/abricate" target="_blank" style="color: #dc3545; font-weight: bold;">Prof. Torsten Seemann</a> 
                            <i class="fas fa-arrow-right"></i> Powered by databases:
                            <a href="https://card.mcmaster.ca/" target="_blank" style="color: #28a745; font-weight: 600;">CARD</a>, 
                            <a href="https://www.mediterranee-infection.com/amr-databases/" target="_blank" style="color: #17a2b8; font-weight: 600;">ARG-ANNOT</a>, 
                            <a href="https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/" target="_blank" style="color: #007bff; font-weight: 600;">NCBI AMR</a>, 
                            <a href="https://megares.meglab.org/" target="_blank" style="color: #fd7e14; font-weight: 600;">MEGARes</a>, 
                            <a href="https://genepi.food.dtu.dk/resfinder" target="_blank" style="color: #20c997; font-weight: 600;">ResFinder</a>
                            <br>
                            <strong>AMRFinderPlus</strong> by <a href="https://github.com/ncbi/amr" target="_blank" style="color: #dc3545; font-weight: bold;">NCBI</a> 
                            <i class="fas fa-arrow-right"></i> Comprehensive resistance gene and point mutation detection.
                            <br>
                            <span style="color: #6c757d;">
                                <i class="fas fa-gratipay"></i> We are deeply grateful to all developers and curators for their open‑source contributions.
                            </span>
                        </span>
                    </div>
                </div>
            </div>
            """
        elif credit_type == 'virulence':
            return """
            <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #fce4ec 100%); border-left: 6px solid #E91E63; margin-bottom: 20px;">
                <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                    <span style="font-size: 1.4em;">🦠</span>
                    <div>
                        <strong style="font-size: 1.2em; color: #880e4f;">Virulence Detection Tools & Databases</strong><br>
                        <span style="font-size: 0.95em; color: #333;">
                            <strong>ABRicate</strong> by <a href="https://github.com/tseemann/abricate" target="_blank" style="color: #E91E63; font-weight: bold;">Prof. Torsten Seemann</a> 
                            <i class="fas fa-arrow-right"></i> Powered by the 
                            <a href="http://www.mgc.ac.cn/VFs/" target="_blank" style="color: #E91E63; font-weight: bold;">VFDB (Virulence Factor Database)</a> 
                            – an open‑source resource for bacterial virulence factors.<br>
                            <span style="color: #6c757d;">
                                <i class="fas fa-gratipay"></i> We thank the VFDB curators and Prof. Seemann for their invaluable contributions.
                            </span>
                        </span>
                    </div>
                </div>
            </div>
            """
        elif credit_type == 'bacmet':
            return """
            <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #fff3e0 100%); border-left: 6px solid #FF5722; margin-bottom: 20px;">
                <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                    <span style="font-size: 1.4em;">🧪</span>
                    <div>
                        <strong style="font-size: 1.2em; color: #bf360c;">BACMET – Biocide & Heavy Metal Resistance</strong><br>
                        <span style="font-size: 0.95em; color: #333;">
                            <strong>ABRicate</strong> by <a href="https://github.com/tseemann/abricate" target="_blank" style="color: #FF5722; font-weight: bold;">Prof. Torsten Seemann</a> 
                            <i class="fas fa-arrow-right"></i> Powered by the 
                            <a href="https://bacmet.biomedicine.gu.se/" target="_blank" style="color: #FF5722; font-weight: bold;">BacMet database</a> 
                            (Pal C, et al. Nucleic Acids Res. 2014;42:D737-43)
                            – an open‑source resource for antibacterial biocide and metal resistance genes.<br>
                            <span style="color: #6c757d;">
                                <i class="fas fa-gratipay"></i> We thank the BacMet curators and Prof. Seemann for their invaluable contributions.
                            </span>
                        </span>
                    </div>
                </div>
            </div>
            """
        elif credit_type == 'plasmids':
            return """
            <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #ede7f6 100%); border-left: 6px solid #673AB7; margin-bottom: 20px;">
                <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                    <span style="font-size: 1.4em;">🔌</span>
                    <div>
                        <strong style="font-size: 1.2em; color: #4a148c;">Plasmid Typing Tools & Databases</strong><br>
                        <span style="font-size: 0.95em; color: #333;">
                            <strong>APT (Acinetobacter Plasmid Typing)</strong> by 
                            <a href="https://github.com/MehradHamidian/AcinetobacterPlasmidTyping" target="_blank" style="color: #673AB7; font-weight: bold;">Hamidian Lab</a>
                            (Lam MMC, et al. Microbiol Spectr. 2023;11(1):e0247822.
                            <a href="https://doi.org/10.1128/spectrum.02478-22" target="_blank" style="color: #673AB7; font-weight: 600;">doi:10.1128/spectrum.02478-22</a>)
                            <i class="fas fa-arrow-right"></i> Species‑specific plasmid typing for <em>Acinetobacter</em>.
                            <br>
                            <strong>PlasmidFinder</strong> by 
                            <a href="https://genepi.food.dtu.dk/PlasmidFinder/" target="_blank" style="color: #673AB7; font-weight: bold;">Center for Genomic Epidemiology (DTU)</a>
                            (Carattoli A, et al. Antimicrob Agents Chemother. 2014;58(7):3895-903.
                            <a href="https://doi.org/10.1128/AAC.02412-14" target="_blank" style="color: #673AB7; font-weight: 600;">doi:10.1128/AAC.02412-14</a>)
                            – generalist database for plasmid replicon detection.
                            <br>
                            <span style="color: #6c757d;">
                                <i class="fas fa-gratipay"></i> We thank the developers and curators for making these resources freely available.
                            </span>
                        </span>
                    </div>
                </div>
            </div>
            """
        return ""

    # ------------------------------------------------------------------------
    # Tab content methods
    # ------------------------------------------------------------------------
    def _summary_section(self, kwargs, carb_count, env_count, mutation_count):
        samples = kwargs['samples_data']
        patterns = kwargs['patterns']
        total = len(samples)
        st_unique = len(patterns.get('pasteur_st_distribution', {}))
        k_unique = len(patterns.get('k_locus_distribution', {}))
        amr_total = sum(len(db) for db in kwargs['gene_centric'].get('amr_databases', {}).values())
        vir_total = sum(len(db) for db in kwargs['gene_centric'].get('virulence_databases', {}).values())
        bac_total = sum(len(db) for db in kwargs['gene_centric'].get('bacmet_databases', {}).values())
        high_risk = len(patterns.get('high_risk_crab', []))
        ic_unique = len(set(d.get('pasteur_mlst',{}).get('International_Clone','Unknown') for d in samples.values() if d.get('pasteur_mlst',{}).get('International_Clone') not in ['Unknown','ND']))
        return f"""
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>Sample‑Centric Analysis with Isolate Boxes</h3>
        <p>This report provides an <strong>interactive per‑isolate view</strong> of all genomic features for <strong>{total}</strong> A. baumannii genomes.</p>
        <ul><li>Each isolate is displayed in its own box with typing badges (Pasteur ST, Oxford ST, K, OCL, Capsule).</li>
        <li>Inside each box, you will find tables with full details from TSV files – horizontally scrollable.</li>
        <li>Filter boxes by sample name or by database (AMR, Virulence, BACMET, Plasmids).</li>
        <li>Gene‑centric tabs (MLST, Kaptive, QC, Combinations, Co‑occurrence, Patterns) are also available.</li></ul>
        </div></div>
        {f'<div class="alert-box alert-danger"><i class="fas fa-exclamation-triangle fa-2x"></i><div><h3>⚠️ Critical Resistance Alert</h3><p><strong>{carb_count} carbapenemase genes</strong> detected across samples. Carbapenem-resistant A. baumannii (CRAB) is a WHO Critical Priority pathogen.</p></div></div>' if carb_count>0 else ''}
        {f'<div class="alert-box alert-info"><i class="fas fa-globe-africa fa-2x"></i><div><h3>⚠️ Environmental Co-Selection Alert</h3><p><strong>{env_count} environmental resistance markers</strong> detected. These genes can co-select for antibiotic resistance in hospital environments.</p></div></div>' if env_count>0 else ''}
        <h3>Key Statistics</h3>
        <div class="master-scrollable-container"><table id="summary-table" class="data-table"><thead><tr><th data-sort="string">Metric</th><th data-sort="number">Count</th><th data-sort="string">Details</th></tr></thead><tbody>
        <tr><td>Total Samples Analyzed</td><td><strong>{total}</strong></td><td>Complete genomic analysis with all databases</td></tr>
        <tr><td>Pasteur STs</td><td><strong>{st_unique}</strong></td><td>MLST typing (Pasteur scheme)</td></tr>
        <tr><td>International Clones</td><td><strong>{ic_unique}</strong></td><td>IC classification</td></tr>
        <tr><td>Capsule (K) Loci</td><td><strong>{k_unique}</strong></td><td>Kaptive capsule typing</td></tr>
        <tr><td>Total AMR Genes</td><td><strong>{amr_total}</strong></td><td>Across all AMR databases</td></tr>
        <tr><td>Carbapenemase Genes</td><td><span class="badge {'badge-critical' if carb_count>0 else 'badge-low'}">{carb_count}</span></td><td>OXA, NDM, VIM, IMP, KPC types</td></tr>
        <tr><td>Virulence Genes</td><td><strong>{vir_total}</strong></td><td>Biofilm formation, virulence factors</td></tr>
        <tr><td>Environmental Markers</td><td><span class="badge {'badge-high' if env_count>0 else 'badge-low'}">{env_count}</span></td><td>Heavy metal, biocide, stress response</td></tr>
        <tr><td>Samples with Mutations</td><td><strong>{mutation_count}</strong></td><td>Point mutations detected (AMRFinderPlus)</td></tr>
        <tr><td>High-Risk CRAB</td><td><span class="badge {'badge-critical' if high_risk>0 else 'badge-low'}">{high_risk}</span></td><td>Carbapenemase + last-resort resistance</td></tr>
        </tbody></table></div>
        <h3>Critical Resistance Categories Tracked</h3>
        <div class="critical-cards" style="display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:20px;margin:20px 0">
            <div class="critical-card" style="background:#fff;padding:20px;border-radius:10px;box-shadow:0 2px 10px rgba(0,0,0,0.1);border-left:4px solid #F44336"><h4><i class="fas fa-skull-crossbones"></i> Carbapenemases ({carb_count})</h4><p>OXA-23, OXA-58, NDM, VIM, IMP, KPC</p></div>
            <div class="critical-card" style="background:#fff;padding:20px;border-radius:10px;box-shadow:0 2px 10px rgba(0,0,0,0.1);border-left:4px solid #9C27B0"><h4><i class="fas fa-tint"></i> Colistin Resistance</h4><p>mcr genes, pmrAB mutations, LPS modifications</p></div>
            <div class="critical-card" style="background:#fff;padding:20px;border-radius:10px;box-shadow:0 2px 10px rgba(0,0,0,0.1);border-left:4px solid #FF9800"><h4><i class="fas fa-pills"></i> Tigecycline Resistance</h4><p>tet(X) variants, efflux pumps (adeABC)</p></div>
            <div class="critical-card" style="background:#fff;padding:20px;border-radius:10px;box-shadow:0 2px 10px rgba(0,0,0,0.1);border-left:4px solid #4CAF50"><h4><i class="fas fa-layer-group"></i> Biofilm Formation</h4><p>ompA, csuABCDE, bfmRS, pili genes</p></div>
            <div class="critical-card" style="background:#fff;padding:20px;border-radius:10px;box-shadow:0 2px 10px rgba(0,0,0,0.1);border-left:4px solid #795548"><h4><i class="fas fa-globe-africa"></i> Environmental Markers ({env_count})</h4><p>Heavy metals, biocides, stress response, plasmid transfer</p></div>
        </div>
        """

    def _samples_section(self, kwargs):
        samples = kwargs['samples_data']
        abricate_details = kwargs.get('abricate_details', {})
        vir_counts = {}
        for sample, dbs in abricate_details.items():
            count = 0
            for db in ['vfdb', 'ecoli_vf']:
                if db in dbs:
                    count += len(dbs[db])
            vir_counts[sample] = count
        html = """
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><h3>Population Structure Overview</h3><p>This table summarises the key typing results for each genome. Understanding the population structure helps identify dominant clones, track outbreaks, and link genotypes to phenotypes.</p><ul><li><strong>MLST (Sequence Type)</strong>: Gold standard for global epidemiology.</li><li><strong>Capsule Typing (K:OCL)</strong>: K (capsule) and OCL (lipooligosaccharide) loci are critical for virulence and immune evasion.</li><li><strong>ND</strong>: Not Determined.</li></ul></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('samples-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById('search-samples').value=''; searchTable('samples-table','search-samples'); highlightTableCells('samples-table','highlight-samples')"><i class="fas fa-sync"></i> Clear Search</button><button class="action-btn btn-success" onclick="exportTableToCSV('samples-table','acinetobacter_samples.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        <div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-samples" onkeyup="searchTable('samples-table','search-samples')" placeholder="🔍 Filter rows..."><input type="text" class="search-box" id="highlight-samples" onkeyup="highlightTableCells('samples-table','highlight-samples')" placeholder="✨ Highlight text..."></div>
        <div class="master-scrollable-container">
            <table id="samples-table" class="data-table">
                <thead>
                    <tr><th data-sort="string">Sample</th><th data-sort="string">Pasteur ST</th><th data-sort="string">Oxford ST</th><th data-sort="string">K Locus</th><th data-sort="string">OCL Locus</th><th data-sort="number">Virulence Count</th></tr>
                </thead>
                <tbody>
        """
        for sample, data in samples.items():
            pasteur = data.get('pasteur_mlst',{}).get('ST','ND')
            oxford = data.get('oxford_mlst',{}).get('ST','ND')
            k = data.get('kaptive',{}).get('K_Locus','ND')
            o = data.get('kaptive',{}).get('O_Locus','ND')
            vcount = vir_counts.get(sample, 0)
            pasteur_disp = f"ST{pasteur}" if pasteur != 'ND' else 'ND'
            oxford_disp = f"ST{oxford}" if oxford != 'ND' else 'ND'
            html += f'<tr><td><strong>{sample}</strong></td><td>{pasteur_disp}</td><td>{oxford_disp}</td><td>{k}</td><td>{o}</td><td>{vcount}</td></tr>'
        html += """
                </tbody>
            </table>
        </div>
        """
        return html

    def _qc_section(self, kwargs):
        qc = kwargs.get('qc_data', {})
        if not qc:
            return '<div class="alert-box alert-warning"><i class="fas fa-exclamation-circle"></i><div>FASTA QC file not found or could not parse.</div></div>'

        html = """
        <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #d1ecf1 100%); border-left: 6px solid #17a2b8; margin-bottom: 20px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">📊</span>
                <div>
                    <strong style="font-size: 1.2em; color: #0c5460;">FASTA QC & Species Confirmation</strong><br>
                    <span style="font-size: 0.95em; color: #333;">
                        <strong>Biopython</strong> by 
                        <a href="https://biopython.org/" target="_blank" style="color: #17a2b8; font-weight: bold;">The Biopython Consortium</a> 
                        <i class="fas fa-arrow-right"></i> 
                        Computes assembly statistics (contigs, N50, GC%, total length) from 
                        <em>FASTA files</em> into a structured format.
                        <br>
                        <strong>FastANI</strong> by 
                        <a href="https://github.com/ParBLiSS/FastANI" target="_blank" style="color: #17a2b8; font-weight: bold;">Jain C, et al.</a> 
                        <i class="fas fa-arrow-right"></i> 
                        Species confirmation via Average Nucleotide Identity (ANI) 
                        (<a href="https://doi.org/10.1038/s41467-018-07641-9" target="_blank" style="color: #17a2b8;">Nat Commun. 2018;9:5114</a>).
                    </span>
                </div>
            </div>
        </div>
        """
        metrics = set()
        for d in qc.values():
            metrics.update(d.keys())
        metrics = sorted(metrics)
        html += '<div class="alert-box alert-info"><i class="fas fa-chart-line"></i><div><h3>FASTA Quality Control</h3><p>Assembly quality metrics for each genome: total sequences, total bases, GC%, N50, etc.</p></div></div>'
        html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'qc-tab\')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById(\'search-qc\').value=\'\'; searchTable(\'qc-table\',\'search-qc\'); highlightTableCells(\'qc-table\',\'highlight-qc\')"><i class="fas fa-sync"></i> Clear Search</button><button class="action-btn btn-success" onclick="exportTableToCSV(\'qc-table\',\'fasta_qc.csv\')"><i class="fas fa-download"></i> Export CSV</button></div>'
        html += '<div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-qc" onkeyup="searchTable(\'qc-table\',\'search-qc\')" placeholder="🔍 Filter rows..."><input type="text" class="search-box" id="highlight-qc" onkeyup="highlightTableCells(\'qc-table\',\'highlight-qc\')" placeholder="✨ Highlight text..."></div>'
        html += '<div class="master-scrollable-container"><table id="qc-table" class="data-table"><thead><tr><th data-sort="string" style="min-width:200px">Sample</th>'
        for m in metrics:
            html += f'<th data-sort="number">{m}</th>'
        html += '</tr></thead><tbody>'
        for sample, vals in sorted(qc.items()):
            html += f'<tr><td style="white-space:nowrap"><strong>{sample}</strong></td>'
            for m in metrics:
                v = vals.get(m, 'ND')
                if isinstance(v, float):
                    v = f"{v:,.0f}" if v > 1000 else f"{v:.2f}"
                html += f'<td>{v}</td>'
            html += '</tr>'
        html += '</tbody></table></div>'
        return html

    def _mlst_section(self, kwargs):
        patterns = kwargs['patterns']
        samples = kwargs['samples_data']
        pasteur = patterns.get('pasteur_st_distribution', {})
        oxford = patterns.get('oxford_st_distribution', {})
        ic_dist = Counter()
        for d in samples.values():
            ic = d.get('pasteur_mlst',{}).get('International_Clone','Unknown')
            if ic not in ['Unknown','ND']:
                ic_dist[ic] += 1
        total_pasteur = sum(pasteur.values()) if pasteur else 0
        total_oxford = sum(oxford.values()) if oxford else 0

        html = """
        <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #fff3cd 100%); border-left: 6px solid #FF9800; margin-bottom: 20px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">🧬</span>
                <div>
                    <strong style="font-size: 1.2em; color: #e65100;">MLST Typing – Acknowledgments & Licensing</strong><br>
                    <span style="font-size: 0.95em; color: #333;">
                        <strong>MLST scheme</strong> powered by 
                        <a href="https://github.com/tseemann/mlst" target="_blank" style="color: #FF9800; font-weight: bold;">Prof. Torsten Seemann's Perl scripts</a> 
                        and the 
                        <a href="https://pubmlst.org/" target="_blank" style="color: #FF9800; font-weight: bold;">PubMLST database</a> 
                        (Jolley et al., Wellcome Open Res 2018).<br>
                        <span style="color: #856404;">
                            <i class="fas fa-info-circle"></i> 
                            <strong>Note:</strong> Due to recent licensing changes and the new PubMLST data policy, 
                            the allele definitions used in this report are current as of <strong>2024</strong>. 
                            Future updates may require manual downloads from the PubMLST website.
                        </span>
                        <br>
                        <span style="font-size: 0.9em; color: #6c757d;">
                            <i class="fas fa-gratipay"></i> We are grateful to the PubMLST curators and Torsten Seemann 
                            for providing these essential open‑source tools.
                        </span>
                    </span>
                </div>
            </div>
        </div>
        """
        html += f"""
        <div class="alert-box alert-info"><i class="fas fa-code-branch"></i><div><h3>MLST Analysis - Dual Scheme (Pasteur & Oxford)</h3><p>Multi‑Locus Sequence Typing (MLST) uses internal fragments of seven housekeeping genes. It is highly reproducible and defines Sequence Types (STs). Closely related STs belong to the same Clonal Complex (CC).</p><p><strong>{len(pasteur)} unique Pasteur STs</strong> and <strong>{len(oxford)} Oxford STs</strong> identified. International Clone (IC) classification helps track global lineages.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('mlst-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById('search-mlst-pasteur').value=''; searchTable('mlst-pasteur-table','search-mlst-pasteur'); highlightGenome('mlst-pasteur-table','highlight-mlst-pasteur')"><i class="fas fa-sync"></i> Clear</button></div>
        <h3>Pasteur MLST Scheme</h3>
        <div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-mlst-pasteur" onkeyup="searchTable('mlst-pasteur-table','search-mlst-pasteur')" placeholder="🔍 Filter STs..."><input type="text" class="search-box" id="highlight-mlst-pasteur" onkeyup="highlightGenome('mlst-pasteur-table','highlight-mlst-pasteur')" placeholder="🔍 Highlight genome tags..."></div>
        <div class="master-scrollable-container"><table id="mlst-pasteur-table" class="data-table"><thead><tr><th data-sort="string">ST</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for st, cnt in sorted(pasteur.items(), key=lambda x: x[1] if isinstance(x[1], int) else x[1].get('count',0), reverse=True):
            if st=='ND': continue
            if isinstance(cnt, dict):
                count = cnt.get('count',0)
                pct = (count/total_pasteur)*100 if total_pasteur else 0
                freq = f"{count} ({pct:.1f}%)"
            else:
                count = cnt
                pct = (count/total_pasteur)*100 if total_pasteur else 0
                freq = f"{count} ({pct:.1f}%)"
            sample_list = [s for s,d in samples.items() if d.get('pasteur_mlst',{}).get('ST')==str(st)]
            tags = self._make_genome_tags(sample_list)
            html += f'<tr><td><strong>ST{st}</strong></td><td>{freq}</td><td><div class="genome-list">{tags}</div></td></tr>'
        html += '</tbody></table></div>'
        html += '<h3>Oxford MLST Scheme</h3><div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-mlst-oxford" onkeyup="searchTable(\'mlst-oxford-table\',\'search-mlst-oxford\')" placeholder="🔍 Filter STs..."><input type="text" class="search-box" id="highlight-mlst-oxford" onkeyup="highlightGenome(\'mlst-oxford-table\',\'highlight-mlst-oxford\')" placeholder="🔍 Highlight genome tags..."></div>'
        html += '<div class="master-scrollable-container"><table id="mlst-oxford-table" class="data-table"><thead><tr><th data-sort="string">ST</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>'
        for st, cnt in sorted(oxford.items(), key=lambda x: x[1] if isinstance(x[1], int) else x[1].get('count',0), reverse=True):
            if st=='ND': continue
            if isinstance(cnt, dict):
                count = cnt.get('count',0)
                pct = (count/total_oxford)*100 if total_oxford else 0
                freq = f"{count} ({pct:.1f}%)"
            else:
                count = cnt
                pct = (count/total_oxford)*100 if total_oxford else 0
                freq = f"{count} ({pct:.1f}%)"
            sample_list = [s for s,d in samples.items() if d.get('oxford_mlst',{}).get('ST')==str(st)]
            tags = self._make_genome_tags(sample_list)
            html += f'<tr><td><strong>ST{st}</strong></td><td>{freq}</td><td><div class="genome-list">{tags}</div></td></tr>'
        html += '</tbody></table></div>'
        html += '<h3>International Clone Distribution</h3><div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-ic" onkeyup="searchTable(\'ic-table\',\'search-ic\')" placeholder="🔍 Filter IC..."><input type="text" class="search-box" id="highlight-ic" onkeyup="highlightGenome(\'ic-table\',\'highlight-ic\')" placeholder="🔍 Highlight genome tags..."></div>'
        html += '<div class="master-scrollable-container"><table id="ic-table" class="data-table"><thead><tr><th data-sort="string">International Clone</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>'
        for ic, cnt in ic_dist.most_common():
            sample_list = [s for s,d in samples.items() if d.get('pasteur_mlst',{}).get('International_Clone')==ic]
            tags = self._make_genome_tags(sample_list)
            html += f'<tr><td><strong>{ic}</strong></td><td>{cnt}</td><td><div class="genome-list">{tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html

    def _kaptive_section(self, kwargs):
        patterns = kwargs['patterns']
        samples = kwargs['samples_data']
        k_dist = patterns.get('k_locus_distribution', {})
        o_dist = patterns.get('o_locus_distribution', {})

        html = """
        <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #f3e5f5 100%); border-left: 6px solid #9C27B0; margin-bottom: 20px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">🧬</span>
                <div>
                    <strong style="font-size: 1.2em; color: #4a148c;">Kaptive Capsule Typing – Acknowledgments</strong><br>
                    <span style="font-size: 0.95em; color: #333;">
                        <strong>Kaptive</strong> developed by 
                        <a href="https://github.com/klebgenomics/Kaptive" target="_blank" style="color: #9C27B0; font-weight: bold;">Stanton et al.</a> 
                        (<a href="https://doi.org/10.1099/mgen.0.001428" target="_blank" style="color: #9C27B0; font-weight: 600;">Microbial Genomics. 2025;11(6):001428</a>)
                        <i class="fas fa-arrow-right"></i> 
                        Identifies capsule (K) and lipooligosaccharide (OCL) loci in <em>Acinetobacter baumannii</em> 
                        from whole‑genome sequence data.
                        <br>
                        <span style="color: #6c757d;">
                            <i class="fas fa-gratipay"></i> We are deeply grateful to the developers for making this tool freely available 
                            for genomic surveillance of A. baumannii.
                        </span>
                    </span>
                </div>
            </div>
        </div>
        """
        html += """
        <div class="alert-box alert-info"><i class="fas fa-shield-alt"></i><div><h3>Capsule Typing (Kaptive) Analysis</h3><p>Kaptive identifies capsule (K) and lipooligosaccharide (OCL) loci in A. baumannii. The capsule is a major virulence determinant, protecting against phagocytosis and complement killing. Specific K types are associated with different clonal complexes and clinical outcomes. OCL (O-Chain Locus) influences immune evasion and serum resistance.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('kaptive-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById('search-k').value=''; searchTable('k-locus-table','search-k'); highlightGenome('k-locus-table','highlight-k')"><i class="fas fa-sync"></i> Clear</button></div>
        """
        html += '<h3>K Locus Distribution</h3><div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-k" onkeyup="searchTable(\'k-locus-table\',\'search-k\')" placeholder="🔍 Filter K locus..."><input type="text" class="search-box" id="highlight-k" onkeyup="highlightGenome(\'k-locus-table\',\'highlight-k\')" placeholder="🔍 Highlight genome tags..."></div>'
        html += '<div class="master-scrollable-container"><table id="k-locus-table" class="data-table"><thead><tr><th data-sort="string">K Locus</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>'
        totalk = sum(k_dist.values()) if k_dist else 0
        for k, cnt in sorted(k_dist.items(), key=lambda x: x[1] if isinstance(x[1], int) else x[1].get('count',0), reverse=True):
            if k=='ND': continue
            if isinstance(cnt, dict):
                count = cnt.get('count',0)
                pct = (count/totalk)*100 if totalk else 0
                freq = f"{count} ({pct:.1f}%)"
            else:
                count = cnt
                pct = (count/totalk)*100 if totalk else 0
                freq = f"{count} ({pct:.1f}%)"
            sample_list = [s for s,d in samples.items() if d.get('kaptive',{}).get('K_Locus')==k]
            tags = self._make_genome_tags(sample_list)
            html += f'<tr><td><strong>{k}</strong></td><td>{freq}</td><td><div class="genome-list">{tags}</div></td></tr>'
        html += '</tbody></table></div>'
        html += '<h3>OCL Locus Distribution</h3><div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-o" onkeyup="searchTable(\'o-locus-table\',\'search-o\')" placeholder="🔍 Filter OCL locus..."><input type="text" class="search-box" id="highlight-o" onkeyup="highlightGenome(\'o-locus-table\',\'highlight-o\')" placeholder="🔍 Highlight genome tags..."></div>'
        html += '<div class="master-scrollable-container"><table id="o-locus-table" class="data-table"><thead><tr><th data-sort="string">OCL Locus</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>'
        totalo = sum(o_dist.values()) if o_dist else 0
        for o, cnt in sorted(o_dist.items(), key=lambda x: x[1] if isinstance(x[1], int) else x[1].get('count',0), reverse=True):
            if o=='ND': continue
            if isinstance(cnt, dict):
                count = cnt.get('count',0)
                pct = (count/totalo)*100 if totalo else 0
                freq = f"{count} ({pct:.1f}%)"
            else:
                count = cnt
                pct = (count/totalo)*100 if totalo else 0
                freq = f"{count} ({pct:.1f}%)"
            sample_list = [s for s,d in samples.items() if d.get('kaptive',{}).get('O_Locus')==o]
            tags = self._make_genome_tags(sample_list)
            html += f'<tr><td><strong>{o}</strong></td><td>{freq}</td><td><div class="genome-list">{tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html

    def _combinations_section(self, kwargs):
        patterns = kwargs['patterns']
        combos = [
            ('st_o_combinations', 'ST - OCL Locus'),
            ('st_k_combinations', 'ST - K Locus'),
            ('ko_combinations', 'K:OCL (Capsule Type)'),
            ('st_ko_combinations', 'ST - K:OCL')
        ]
        html = '<div class="alert-box alert-info"><i class="fas fa-link"></i><div><h3>Combination Tables</h3><p>Associations between Sequence Types, capsule loci, and serotypes. These combinations help identify dominant clones and track epidemic strains.</p></div></div>'
        html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'combinations-tab\')"><i class="fas fa-print"></i> Print Section</button></div>'
        for key, title in combos:
            data = patterns.get(key, {})
            if not data:
                continue
            html += f'<h3>{title}</h3><div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-{key}" onkeyup="searchTable(\'{key}-table\',\'search-{key}\')" placeholder="🔍 Filter rows..."><input type="text" class="search-box" id="highlight-{key}" onkeyup="highlightGenome(\'{key}-table\',\'highlight-{key}\')" placeholder="🔍 Highlight genome tags..."></div>'
            html += f'<div class="master-scrollable-container"><table id="{key}-table" class="data-table"><thead><tr><th data-sort="string">{title}</th><th data-sort="number">Count</th><th data-sort="string">Samples</th><tr></thead><tbody>'
            for combo, samples_list in sorted(data.items(), key=lambda x: len(x[1]), reverse=True):
                tags = self._make_genome_tags(samples_list)
                html += f'<tr><td><strong>{combo}</strong></td><td>{len(samples_list)}</td><td><div class="genome-list">{tags}</div></td></tr>'
            html += '</tbody></table></div>'
        return html

    def _amr_boxes(self, kwargs):
        return self._sample_boxes_generic(kwargs, 'amr', 'AMR',
                                          db_list=['amrfinder', 'card', 'resfinder', 'argannot', 'megares', 'ncbi'],
                                          credit_type='amr')

    def _virulence_boxes(self, kwargs):
        return self._sample_boxes_generic(kwargs, 'virulence', 'Virulence',
                                          db_list=['vfdb'],
                                          credit_type='virulence')

    def _bacmet_boxes(self, kwargs):
        return self._sample_boxes_generic(kwargs, 'bacmet', 'BACMET',
                                          db_list=['bacmet2'],
                                          credit_type='bacmet')

    def _plasmids_boxes(self, kwargs):
        return self._sample_boxes_generic(kwargs, 'plasmids', 'Plasmids',
                                          db_list=['plasmidfinder'],
                                          credit_type='plasmids',
                                          include_apt=True)

    def _mutation_boxes(self, kwargs):
        samples_data = kwargs.get('samples_data', {})
        mutation_details = kwargs.get('mutation_details', {})
        relevant_samples = [s for s in samples_data if s in mutation_details and mutation_details[s]]
        relevant_samples.sort()

        if not relevant_samples:
            return """
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No Mutation Data Available</h3><p>No mutations found for any sample.</p></div>
            </div>
            """

        html = """
        <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #e0f7fa 100%); border-left: 6px solid #00BCD4; margin-bottom: 20px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">🧬</span>
                <div>
                    <strong style="font-size: 1.2em; color: #006064;">AMRFinderPlus – Point Mutation Detection</strong><br>
                    <span style="font-size: 0.95em; color: #333;">
                        <strong>AMRFinderPlus</strong> developed by 
                        <a href="https://github.com/ncbi/amr" target="_blank" style="color: #00BCD4; font-weight: bold;">NCBI</a>
                        (Feldgarden M, et al. Sci Rep. 2021;11(1):12728.
                        <a href="https://doi.org/10.1038/s41598-021-91456-0" target="_blank" style="color: #00BCD4; font-weight: 600;">doi:10.1038/s41598-021-91456-0</a>)
                        <i class="fas fa-arrow-right"></i> 
                        Comprehensive detection of antimicrobial resistance genes and point mutations.
                        <br>
                        <span style="color: #6c757d;">
                            <i class="fas fa-gratipay"></i> We thank the NCBI team for maintaining this essential resource for genomic surveillance.
                        </span>
                    </span>
                </div>
            </div>
        </div>
        """
        html += """
        <div class="alert-box alert-info">
            <i class="fas fa-dna fa-2x"></i>
            <div>
                <h3>🧬 Point Mutations – Sample‑Centric View</h3>
                <p>Each box represents one isolate. Inside, you will find a table of all detected point mutations with full details – horizontally scrollable.</p>
                <p>Use the search bar below to filter boxes by sample name.</p>
            </div>
        </div>
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-mutation" onkeyup="filterBoxes('mutation')" placeholder="🔍 Search sample...">
            <button class="action-btn btn-success" onclick="resetBoxFilters('mutation')"><i class="fas fa-sync"></i> Clear Filters</button>
        </div>
        <div id="box-container-mutation">
        """
        for sample in relevant_samples:
            html += self._build_isolate_box(sample, kwargs,
                                            db_list=['mutations'],
                                            mutation_details=mutation_details)
        html += """
        </div>
        """
        return html

    def _cooccurrence_section(self, kwargs):
        cooc = kwargs['patterns'].get('gene_cooccurrence', {})
        if not cooc:
            return '<div class="alert-box alert-warning">Co‑occurrence data not available.</div>'
        pairs = []
        for g1, partners in cooc.items():
            for g2, cnt in partners.items():
                pairs.append((g1, g2, cnt))
        pairs.sort(key=lambda x: x[2], reverse=True)
        top = pairs[:500]
        html = f"""
        <div class="alert-box alert-info"><i class="fas fa-link"></i><div><h3>Gene Co‑occurrence (Top 500)</h3><p>Pairs of genes that appear together in the same genomes. Higher counts indicate potential synergy or linked mobilisation.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('cooccurrence-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-success" onclick="exportTableToCSV('cooccurrence-table','cooccurrence.csv')"><i class="fas fa-download"></i> Export Co‑occurrence</button></div>
        <input type="text" class="search-box" id="search-cooc" onkeyup="searchTable('cooccurrence-table','search-cooc')" placeholder="🔍 Search gene pairs...">
        <div class="master-scrollable-container"><table id="cooccurrence-table" class="data-table"><thead><tr><th data-sort="string">Gene 1</th><th data-sort="string">Gene 2</th><th data-sort="number">Co‑occurrence Count</th></tr></thead><tbody>
        """
        for g1, g2, cnt in top:
            html += f'<tr><td>{g1}</td><td>{g2}</td><td>{cnt}</td></tr>'
        html += '</tbody></table></div>'
        return html

    def _patterns_section(self, kwargs):
        patterns = kwargs['patterns']
        crab = patterns.get('high_risk_crab', [])
        caveat = """
        <div class="alert-box alert-warning" style="border-left-color:#ffc107;">
            <i class="fas fa-exclamation-triangle fa-2x"></i>
            <div>
                <strong>⚠️ Important Caveat:</strong> The identification of CRAB (carbapenem‑resistant <em>A. baumannii</em>) in this section is based <strong>solely on the presence of carbapenemase genes</strong> identified in the genomic data. 
                This <strong>does not</strong> replace phenotypic antimicrobial susceptibility testing (AST). 
                Clinical decisions must always be guided by AST results and local guidelines.
            </div>
        </div>
        """
        html = f"""
        {caveat}
        <div class="alert-box alert-info"><i class="fas fa-project-diagram"></i><div><h3>Pattern Discovery – CRAB (Carbapenem‑resistant A. baumannii)</h3><p>Carbapenem‑resistant A. baumannii (CRAB) isolates with last‑resistance genes (colistin/tigecycline). These represent high‑risk, difficult‑to‑treat infections.</p></div></div>
        """
        if crab:
            html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'patterns-tab\')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById(\'search-patterns\').value=\'\'; searchTable(\'crab-table\',\'search-patterns\'); highlightTableCells(\'crab-table\',\'highlight-patterns\')"><i class="fas fa-sync"></i> Clear</button></div>'
            html += '<div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-patterns" onkeyup="searchTable(\'crab-table\',\'search-patterns\')" placeholder="🔍 Filter rows..."><input type="text" class="search-box" id="highlight-patterns" onkeyup="highlightTableCells(\'crab-table\',\'highlight-patterns\')" placeholder="✨ Highlight text..."></div>'
            html += '<div class="master-scrollable-container"><table id="crab-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">ST</th><th data-sort="string">Capsule (K:OCL)</th><th data-sort="string">Carbapenemases</th><th data-sort="string">Colistin Resistance</th><th data-sort="string">Tigecycline Resistance</th></tr></thead><tbody>'
            for c in crab:
                html += f'<tr><td><strong>{c["sample"]}</strong></td><td>{c["pasteur_st"]}</td><td>{c["capsule_type"]}</td><td>{",".join(c["carbapenemases"])}</td><td>{",".join(c["colistin_resistance"])}</td><td>{",".join(c["tigecycline_resistance"])}</td></tr>'
            html += '</tbody></table></div>'
        else:
            html += '<p>No high‑risk CRAB isolates detected.</p>'
        return html

    def _aiguide_section(self):
            return """
            <div class="alert-box alert-info"><i class="fas fa-robot fa-2x"></i><div><h3>AI Assistant Guide – Unleash the Power of AI for Genomic Epidemiology</h3><p>This guide shows you how to leverage large language models (LLMs) like <strong>ChatGPT, Claude, or Gemini</strong> to interact with your <em>A. baumannii</em> dataset – turning static reports into dynamic conversations.</p></div></div>
            <div style="margin: 20px 0;">
                <div class="database-section"><h4><i class="fas fa-brain"></i> Why Use AI for Genomic Data Analysis?</h4><p>Modern AI models excel at:</p><ul><li><strong>Pattern recognition</strong> – spotting epidemiological trends, clone associations, and co‑occurrence networks.</li><li><strong>Natural language queries</strong> – ask in plain English, get instant answers without writing code.</li><li><strong>Hypothesis generation</strong> – uncover unexpected correlations that merit experimental follow‑up.</li><li><strong>Literature synthesis</strong> – connect your findings with published resistance mechanisms and clinical guidelines.</li></ul><p><span style="background:#fff3cd; padding:2px 8px; border-radius:4px;"><i class="fas fa-lightbulb"></i> <strong>Scientific note:</strong> AI is <em>pattern‑finding</em>, not <em>causal‑inferring</em>. Use it to suggest, then verify with wet‑lab or clinical correlation.</span></p></div>
                <div class="database-section"><h4><i class="fas fa-upload"></i> How to Feed This Report to AI</h4><p>You have <strong>three powerful options</strong>:</p><ol><li><strong>Upload the JSON file</strong> – <code>genius_acinetobacter_ultimate_gene_centric_report.json</code> contains all structured data. Upload it to ChatGPT (Advanced Data Analysis), Claude, or Gemini. <em>Best for precise quantitative queries.</em></li><li><strong>Upload the HTML report</strong> – Modern AI tools can read HTML and extract tables. Upload the <code>.html</code> file directly – the AI will parse the tables and text. <em>Great for visual context.</em></li><li><strong>Copy‑paste specific tables</strong> – If you only need a quick insight, copy a table (e.g., AMR gene list) and paste it into the chat. <em>Instant, no file upload needed.</em></li></ol><p style="margin-top:10px; background:#e8f5e9; padding:10px; border-radius:5px;"><i class="fas fa-info-circle"></i> <strong>Pro tip:</strong> For best results, tell the AI: <em>"You are a bioinformatician analysing <em>A. baumannii</em> genomes. The attached data contains typing, AMR, virulence, BACMET, and mutation information. Answer my questions with references to the data."</em></p></div>
                <div class="database-section"><h4><i class="fas fa-chart-line"></i> Scientifically Relevant Questions to Ask</h4><div style="display:grid; grid-template-columns:1fr 1fr; gap:10px;"><div style="background:#f8f9fa; padding:10px; border-radius:8px;"><strong>🧬 Epidemiology &amp; Clonality</strong><ul style="margin-top:5px; font-size:0.9em;"><li>What are the most common Pasteur STs in this dataset?</li><li>Which capsule types (K:OCL) are dominant?</li><li>Are there any ST‑capsule combinations with >2 isolates?</li><li>Which clones carry the most resistance genes?</li></ul></div><div style="background:#f8f9fa; padding:10px; border-radius:8px;"><strong>💊 Antimicrobial Resistance</strong><ul style="margin-top:5px; font-size:0.9em;"><li>How many samples carry OXA-23? What are their STs?</li><li>Are there any NDM‑positive samples? What is their capsule type?</li><li>Which AMR genes co‑occur most frequently?</li><li>What is the distribution of colistin resistance genes (mcr, pmr)?</li></ul></div><div style="background:#f8f9fa; padding:10px; border-radius:8px;"><strong>🦠 Virulence &amp; Biofilm</strong><ul style="margin-top:5px; font-size:0.9em;"><li>Which samples carry biofilm genes (ompA, csu)?</li><li>List all isolates with the T6SS (tss, hcp).</li><li>Is there a correlation between capsule type and virulence gene carriage?</li></ul></div><div style="background:#f8f9fa; padding:10px; border-radius:8px;"><strong>🧪 Mutations &amp; Biocides</strong><ul style="margin-top:5px; font-size:0.9em;"><li>What are the most frequent point mutations in gyrA or parC?</li><li>Are there any 23S rRNA mutations (tigecycline resistance)?</li><li>Which samples carry qac genes (disinfectant resistance)?</li></ul></div></div><p style="margin-top:10px;"><i class="fas fa-arrow-right"></i> <strong>Beyond tables:</strong> Ask the AI to <em>“write a summary of the resistance profile for ST2”</em> or <em>“compare virulence gene carriage between K1 and K2 capsule types”</em> – the report contains all data.</p></div>
                <div class="database-section"><h4><i class="fas fa-balance-scale"></i> Scientific Rigour &amp; Ethical AI Use</h4><ul><li><strong>AI is your co‑pilot, not the pilot.</strong> Always interpret AI‑generated insights in the context of your local epidemiology, clinical guidelines, and laboratory validation.</li><li><strong>Verify, verify, verify.</strong> Cross‑check critical calls (e.g., carbapenemase presence, colistin resistance) with primary literature, genome browsers, or secondary tools.</li><li><strong>No patient‑identifiable data.</strong> Only upload aggregated, de‑identified genomic data. This report contains no clinical metadata.</li><li><strong>Transparency in publications.</strong> If you use AI for exploratory analysis, mention it in the methods (e.g., “AI‑assisted pattern discovery was performed using a large language model, followed by manual curation”).</li><li><strong>AI hallucination is real.</strong> If the AI confidently tells you that <em>blaOXA-23</em> is found in <em>E. coli</em> or that tigecycline resistance is caused by a mutation in <em>rpoB</em> – <strong>don’t believe it</strong>. Treat every AI statement as a hypothesis, not a fact.</li></ul></div>
                <div class="database-section" style="background:#fff3cd; border-left:6px solid #ffc107;"><h4><i class="fas fa-smile-wink"></i> A (Mostly Serious) AI Survival Guide</h4><ul><li><strong>If the AI says “I don’t know”</strong> – trust it. It’s being honest.</li><li><strong>If the AI says “It is widely known”</strong> – ask for a reference. It may have made it up.</li><li><strong>If the AI offers a completely novel evolutionary theory</strong> – check if your coffee is spiked.</li><li><strong>Remember:</strong> AI won’t take your job – but a microbiologist who knows how to use AI might! So learn it, use it, and always keep a healthy dose of scepticism. 😉</li></ul><p style="margin-top:10px;"><i class="fas fa-microbe"></i> <strong>Final thought:</strong> The best AI‑human partnership is one where the AI does the heavy pattern‑lifting, and you do the heavy thinking. Happy (and careful) exploring!</p></div>
            </div>
            """


    def _citation_section(self):
        citations = [
            {
                "title": "AcinetoScope",
                "text": "Beckley et al. A species‑specific computational pipeline for rapid, comprehensive Acinetobacter baumannii outbreak investigation and resistance gene tracking. GitHub 2026.",
                "url": "https://github.com/bbeckley-hub/acinetoscope",
                "category": "main"
            },
            {
                "title": "MLST",
                "text": "Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub. 2018.",
                "url": "https://github.com/tseemann/mlst",
                "category": "tool"
            },
            {
                "title": "PubMLST",
                "text": "Jolley KA, Bray JE, Maiden MCJ. Open‑access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res. 2018;3:124.",
                "url": "https://doi.org/10.12688/wellcomeopenres.14826.1",
                "category": "db"
            },
            {
                "title": "Kaptive 3",
                "text": "Stanton TD, Hetland MAK, Löhr IH, Holt KE, Wyres KL. Fast and Accurate in silico Antigen Typing with Kaptive 3. Microbial Genomics. 2025;11(6):001428.",
                "url": "https://doi.org/10.1099/mgen.0.001428",
                "category": "tool"
            },
            {
                "title": "APT (Acinetobacter Plasmid Typing)",
                "text": "Lam MMC, Koong J, Holt KE, Hall RM, Hamidian M. Detection and Typing of Plasmids in Acinetobacter baumannii Using rep Genes Encoding Replication Initiation Proteins. Microbiol Spectr. 2023;11(1):e0247822.",
                "url": "https://doi.org/10.1128/spectrum.02478-22",
                "category": "tool"
            },
            {
                "title": "AMRFinderPlus",
                "text": "Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11(1):12728.",
                "url": "https://doi.org/10.1038/s41598-021-91456-0",
                "category": "tool"
            },
            {
                "title": "ABRicate",
                "text": "Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub. 2024.",
                "url": "https://github.com/tseemann/abricate",
                "category": "tool"
            },
            {
                "title": "FastANI",
                "text": "Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat Commun. 2018;9(1):5114.",
                "url": "https://doi.org/10.1038/s41467-018-07641-9",
                "category": "tool"
            },
            {
                "title": "CARD",
                "text": "Alcock BP, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Res. 2023;51(D1):D690-D699.",
                "url": "https://doi.org/10.1093/nar/gkac920",
                "category": "db"
            },
            {
                "title": "ResFinder",
                "text": "Bortolaia V, et al. ResFinder 4.0 for predictions of phenotypes from genotypes. J Antimicrob Chemother. 2020;75(12):3491-3500.",
                "url": "https://doi.org/10.1093/jac/dkaa345",
                "category": "db"
            },
            {
                "title": "VFDB",
                "text": "Chen L, Zheng D, Liu B, Yang J, Jin Q. VFDB 2016: hierarchical and refined dataset for big data analysis – 10 years on. Nucleic Acids Res. 2016;44(D1):D694-D697.",
                "url": "https://doi.org/10.1093/nar/gkv1239",
                "category": "db"
            },
            {
                "title": "PlasmidFinder",
                "text": "Carattoli A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. Antimicrob Agents Chemother. 2014;58(7):3895-903.",
                "url": "https://doi.org/10.1128/AAC.02412-14",
                "category": "db"
            },
            {
                "title": "EcOH",
                "text": "Inouye M, et al. SRST2: Rapid genomic surveillance for public health and hospital microbiology labs. Microb Genom. 2016;2(7):e000064.",
                "url": "https://doi.org/10.1099/mgen.0.000064",
                "category": "db"
            },
            {
                "title": "MEGARes 3.0",
                "text": "Bonin N, et al. MEGARes 3.0: a curated and updated database for antimicrobial resistance, biocide, and metal resistance gene detection. Nucleic Acids Res. 2023;51(D1):D677-D684.",
                "url": "https://doi.org/10.1093/nar/gkac1047",
                "category": "db"
            },
            {
                "title": "ARG-ANNOT",
                "text": "Gupta SK, et al. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. Antimicrob Agents Chemother. 2014;58(1):212-20.",
                "url": "https://doi.org/10.1128/AAC.01310-13",
                "category": "db"
            },
            {
                "title": "BacMet",
                "text": "Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. Nucleic Acids Res. 2014;42(Database issue):D737-43.",
                "url": "https://doi.org/10.1093/nar/gkt1252",
                "category": "db"
            },
            {
                "title": "Biopython",
                "text": "Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-3.",
                "url": "https://doi.org/10.1093/bioinformatics/btp163",
                "category": "other"
            }
        ]

        html = """
        <div class="alert-box alert-info"><i class="fas fa-quote-right fa-2x"></i><div><h3>Citations & Acknowledgements</h3><p>Please cite the following tools and databases if you use this report in your research.</p></div></div>
        <div style="margin: 20px 0;">
        """
        for cite in citations:
            cat_class = cite['category']
            badge_text = {'main': 'Main Tool', 'tool': 'Tool', 'db': 'Database', 'other': 'Other'}[cat_class]
            html += f"""
            <div class="citation-card cite-{cat_class}">
                <span class="cite-badge">{badge_text}</span>
                <span class="cite-title">{cite['title']}</span>
                <div class="cite-text">{cite['text']}</div>
                <div>
                    <a href="{cite['url']}" target="_blank" class="cite-link"><i class="fas fa-external-link-alt"></i> {cite['url']}</a>
                    <button class="copy-btn" data-citation="{cite['text'].replace('"', '&quot;')}">📋 Copy</button>
                </div>
            </div>
            """
        html += """
        </div>
        <div class="alert-box alert-success" style="margin-top:20px;"><i class="fas fa-hand-peace"></i><div><strong>Suggested acknowledgement:</strong><br>
        "Genomic analysis was performed using AcinetoScope [Beckley et al., 2026], which integrates MLST [Seemann, 2018] using the PubMLST database [Jolley et al., 2018], Kaptive 3 [Stanton et al., 2025], ABRicate [Seemann, 2024], AMRFinderPlus [Feldgarden et al., 2021], APT (Acinetobacter Plasmid Typing) [Lam et al., 2023], PlasmidFinder [Carattoli et al., 2014], and FastANI [Jain et al., 2018]. Antimicrobial resistance genes were identified using CARD [Alcock et al., 2023], ResFinder [Bortolaia et al., 2020], MEGARes 3.0 [Bonin et al., 2023], and ARG-ANNOT [Gupta et al., 2014]. For biocide and heavy metal resistance genes, BacMet [Pal et al., 2014] was used. Virulence screening was performed with ABRicate using VFDB [Chen et al., 2016]. Mutation detection was performed using AMRFinderPlus. FASTA QC was performed using Biopython [Cock et al., 2009]."
        </div></div>
        <script>
        document.querySelectorAll('.copy-btn').forEach(btn => {
            btn.addEventListener('click', function() {
                const text = this.getAttribute('data-citation');
                navigator.clipboard.writeText(text).then(() => {
                    const original = this.innerHTML;
                    this.innerHTML = '✓ Copied!';
                    setTimeout(() => { this.innerHTML = original; }, 2000);
                });
            });
        });
        </script>
        """
        return html

    def _funding_section(self):
            return """
            <div class="alert-box alert-info"><i class="fas fa-coffee fa-2x"></i><div><h3>Funding & Support – Keeping the Lights On (with code and caffeine)</h3><p>AcinetoScope is an <strong>independent, unfunded project</strong> born out of passion for genomic surveillance and AMR research at the University of Ghana Medical School.</p><p>No grants, no sponsors, no institutional backing – just a laptop, a lot of coffee, and a burning desire to help researchers fight antimicrobial resistance.</p></div></div>
            <div class="alert-box alert-warning"><i class="fas fa-heart fa-2x"></i><div><h3>💡 How You Can Help (Without Opening Your Wallet)</h3><ul><li><strong>⭐ Star us on GitHub</strong> – It takes two seconds and makes us feel like rockstars.</li><li><strong>🐛 Report bugs</strong> – If something breaks, let us know. We’ll fix it with joy.</li><li><strong>💡 Suggest features</strong> – Have an idea? We’re all ears.</li><li><strong>🧬 Share your data</strong> – If you’ve used the tool and want to collaborate, we’d love to hear your story.</li><li><strong>📢 Spread the word</strong> – Tell your colleagues, tweet about it, or mention it in your next Zoom call.</li><li><strong>👋 Say hello!</strong> – Seriously, just drop an email to <strong>brownbeckley94@gmail.com</strong>. It makes our day.</li></ul><p><i class="fas fa-microbe"></i> <strong>Fun fact:</strong> This project runs on 100% volunteer tears, 0% grant money. But we’re not bitter – we’re just caffeinated.</p></div></div>
            <div class="alert-box alert-success"><i class="fas fa-hand-holding-heart"></i><div><h3>🤝 Contribute to the ESKAPE AMR Platform</h3><p>We also maintain pipelines for other ESKAPE pathogens (Kleboscope, Pseudoscope, etc.). If you’re a developer, bioinformatician, or just someone who loves clean code and bacteria, we welcome: pull requests, issues, documentation improvements, ideas for new databases. Visit our GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a> – star, fork, and let’s fight AMR together!</p><p><strong>Brown Beckley</strong> – <i class="fas fa-envelope"></i> brownbeckley94@gmail.com</p><p><i class="fas fa-laugh-beam"></i> <strong>P.S.</strong> If you ever meet Brown in person, buy him a coffee – he’ll probably talk your ear off about ST types, but it’s worth it.</p></div></div>
            """

    def _export_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-download"></i><div>Export data tables as CSV or download complete JSON.</div></div>
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','acinetobacter_samples.csv')">Sample Overview CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('amr-table','amr_genes.csv')">AMR Genes CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('vir-table','virulence_genes.csv')">Virulence Genes CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('bac-table','bacmet_genes.csv')">Bacmet CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('apt-rep-table','apt_plasmid_rep_freq.csv')">APT Plasmid Rep Frequency CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('plasmidfinder-table','plasmidfinder_markers.csv')">PlasmidFinder CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('mutation-table','mutations.csv')">Mutations CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('cooccurrence-table','cooccurrence.csv')">Co‑occurrence CSV</button>
            <button class="action-btn btn-success" onclick="location.href='genius_acinetobacter_sample_centric_report.json'">Download JSON</button>
        </div>
        """


class SampleCentricReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "GENIUS_ACINETOBACTER_SAMPLE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = SampleCentricParser()
        self.analyzer = SampleCentricAnalyzer()
        self.html_generator = SampleCentricHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "GENIUS Acinetobacter baumannii Sample‑Centric Reporter",
            "version": "1.1.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }

    def run(self):
        print("=" * 80)
        print("🧬 GENIUS ACINETOBACTER SAMPLE‑CENTRIC REPORTER v1.1.0")
        print("   Interactive per‑isolate boxes for AMR, Virulence, BACMET, Plasmids, Mutations")
        print("   Gene‑centric for MLST, Kaptive, QC, Combinations, Co‑occurrence, Patterns")
        print("=" * 80)
        print(f"📁 Input directory: {self.input_dir}")

        print("\n📂 Loading TSV files...")
        amr_details, _ = self.parser.load_amrfinder_tsv(self.input_dir)
        abricate_details, _ = self.parser.load_abricate_tsv(self.input_dir)
        apt_data = self.parser.load_apt_tsv(self.input_dir)
        mutation_details = self.parser.load_mutations_tsv(self.input_dir)

        print("\n📄 Parsing HTML reports...")
        pasteur = self.parser.parse_mlst_report(self.input_dir / 'pasteur_mlst_summary.html', 'pasteur')
        oxford = self.parser.parse_mlst_report(self.input_dir / 'oxford_mlst_summary.html', 'oxford')
        kaptive = self.parser.parse_kaptive_report(self.input_dir / 'Kaptive_summary.html')
        qc_data = self.parser.parse_qc_report(self.input_dir / 'FASTA_QC_summary.html')

        all_samples = set(pasteur.keys()) | set(oxford.keys()) | set(kaptive.keys())
        all_samples.update(amr_details.keys())
        all_samples.update(abricate_details.keys())
        all_samples.update(apt_data.keys())
        all_samples.update(mutation_details.keys())
        all_samples = sorted(list(all_samples))
        total_samples = len(all_samples)
        print(f"\n📊 Found {total_samples} unique samples")

        samples_data = {}
        for sample in all_samples:
            samples_data[sample] = {
                'pasteur_mlst': pasteur.get(sample, {'ST': 'ND', 'International_Clone': 'Unknown', 'Allele_Profile': 'ND'}),
                'oxford_mlst': oxford.get(sample, {'ST': 'ND', 'Allele_Profile': 'ND'}),
                'kaptive': kaptive.get(sample, {'K_Locus': 'ND', 'O_Locus': 'ND', 'Capsule_Type': 'ND'})
            }

        gene_centric = self.analyzer.build_gene_centric_frequencies(amr_details, abricate_details, total_samples)

        patterns = self.analyzer.build_cross_genome_patterns(samples_data, gene_centric)

        sample_genes = defaultdict(set)
        for sample, genes in amr_details.items():
            for g in genes:
                sample_genes[sample].add(g.get('gene', 'unknown'))
        for sample, dbs in abricate_details.items():
            for db, genes in dbs.items():
                for g in genes:
                    sample_genes[sample].add(g.get('gene', 'unknown'))

        cooc = defaultdict(Counter)
        for genes in sample_genes.values():
            gene_list = list(genes)
            for i, g1 in enumerate(gene_list):
                for g2 in gene_list[i+1:]:
                    cooc[g1][g2] += 1
        patterns['gene_cooccurrence'] = dict(cooc)

        high_risk_crab = []
        for sample, genes in sample_genes.items():
            carb = [g for g in genes if self.analyzer.categorize_gene(g) == 'Carbapenemases']
            colistin = [g for g in genes if self.analyzer.categorize_gene(g) == 'Colistin Resistance']
            tigecycline = [g for g in genes if self.analyzer.categorize_gene(g) == 'Tigecycline Resistance']
            if carb and (colistin or tigecycline):
                data = samples_data.get(sample, {})
                high_risk_crab.append({
                    'sample': sample,
                    'pasteur_st': data.get('pasteur_mlst', {}).get('ST', 'ND'),
                    'k_locus': data.get('kaptive', {}).get('K_Locus', 'ND'),
                    'capsule_type': data.get('kaptive', {}).get('Capsule_Type', 'ND'),
                    'carbapenemases': carb,
                    'colistin_resistance': colistin,
                    'tigecycline_resistance': tigecycline
                })
        patterns['high_risk_crab'] = high_risk_crab

        integrated_data = {
            'metadata': self.metadata,
            'samples': samples_data,
            'patterns': patterns,
            'gene_centric': gene_centric,
            'qc_data': qc_data,
            'amr_details': amr_details,
            'abricate_details': abricate_details,
            'apt_data': apt_data,
            'mutation_details': mutation_details
        }

        print("\n📊 Generating reports...")
        self._generate_json(integrated_data)
        self._generate_csv(integrated_data)
        self.html_generator.generate_report(integrated_data, self.output_dir)

        print("\n✅ Analysis complete!")
        print(f"📁 Output directory: {self.output_dir}")
        print("🌐 Open genius_acinetobacter_sample_centric_report.html in your browser.")
        return True

    def _generate_json(self, data: Dict[str, Any]) -> Path:
        print("   Generating JSON...")
        output_file = self.output_dir / "genius_acinetobacter_sample_centric_report.json"
        def convert_keys(obj):
            if isinstance(obj, dict):
                new = {}
                for k, v in obj.items():
                    if isinstance(k, tuple):
                        k = '|'.join(k)
                    new[k] = convert_keys(v)
                return new
            elif isinstance(obj, (list, tuple, set)):
                return [convert_keys(i) for i in obj]
            elif isinstance(obj, (defaultdict, Counter)):
                return convert_keys(dict(obj))
            else:
                return obj
        serializable = convert_keys(data)
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(serializable, f, indent=2, default=str)
        print(f"      ✅ JSON saved: {output_file}")
        return output_file

    def _generate_csv(self, data: Dict[str, Any]):
        print("   Generating CSVs...")
        samples_data = data['samples']
        abricate_details = data.get('abricate_details', {})
        apt_data = data.get('apt_data', {})
        mutation_details = data.get('mutation_details', {})
        patterns = data['patterns']

        vir_counts = {}
        for sample, dbs in abricate_details.items():
            count = 0
            for db in ['vfdb', 'ecoli_vf']:
                if db in dbs:
                    count += len(dbs[db])
            vir_counts[sample] = count

        samples_df = pd.DataFrame([{
            'Sample': s,
            'Pasteur_ST': d.get('pasteur_mlst', {}).get('ST', 'ND'),
            'Oxford_ST': d.get('oxford_mlst', {}).get('ST', 'ND'),
            'K_Locus': d.get('kaptive', {}).get('K_Locus', 'ND'),
            'OCL_Locus': d.get('kaptive', {}).get('O_Locus', 'ND'),
            'Virulence_Count': vir_counts.get(s, 0)
        } for s, d in samples_data.items()])
        samples_df.to_csv(self.output_dir / "sample_overview.csv", index=False)

        for cat, fname in [('amr_databases', 'amr_genes.csv'),
                           ('virulence_databases', 'virulence_genes.csv'),
                           ('bacmet_databases', 'bacmet_genes.csv'),
                           ('plasmid_databases', 'plasmid_replicons.csv')]:
            rows = []
            for db, genes in data['gene_centric'].get(cat, {}).items():
                for g in genes:
                    rows.append({'Gene': g['gene'], 'Database': g['database'], 'Count': g['count'], 'Frequency': g['frequency_display']})
            if rows:
                pd.DataFrame(rows).to_csv(self.output_dir / fname, index=False)

        if mutation_details:
            mut_rows = []
            for sample, muts in mutation_details.items():
                for m in muts:
                    mut_rows.append({'Sample': sample, **m})
            pd.DataFrame(mut_rows).to_csv(self.output_dir / "mutations.csv", index=False)

        if apt_data:
            apt_rows = [{'Genome': s, 'HitCount': info['hit_count'], 'RepTypes': ','.join(info['rep_types'])} for s, info in apt_data.items()]
            pd.DataFrame(apt_rows).to_csv(self.output_dir / "apt_plasmid_summary.csv", index=False)

        cooc = patterns.get('gene_cooccurrence', {})
        if cooc:
            cooc_rows = []
            for g1, partners in cooc.items():
                for g2, cnt in partners.items():
                    cooc_rows.append({'Gene1': g1, 'Gene2': g2, 'Count': cnt})
            cooc_rows.sort(key=lambda x: x['Count'], reverse=True)
            pd.DataFrame(cooc_rows[:500]).to_csv(self.output_dir / "cooccurrence.csv", index=False)

        high_risk = patterns.get('high_risk_crab', [])
        if high_risk:
            pd.DataFrame(high_risk).to_csv(self.output_dir / "high_risk_crab.csv", index=False)

        print("      ✅ CSVs saved.")


def main():
    parser = argparse.ArgumentParser(description='GENIUS Acinetobacter baumannii Sample‑Centric Reporter v1.1.0')
    parser.add_argument('-i', '--input-dir', required=True, help='Directory containing AcinetoScope TSV and HTML files')
    parser.add_argument('-o', '--output-dir', help='Custom output directory (optional)')
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        sys.exit(1)
    reporter = SampleCentricReporter(input_dir)
    if args.output_dir:
        reporter.output_dir = Path(args.output_dir)
        reporter.output_dir.mkdir(parents=True, exist_ok=True)
    success = reporter.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()