#!/usr/bin/env python3
"""
GENIUS ACINETOBACTER BAUMANNII ULTIMATE REPORTER
Advanced HTML Parser with Gene-Centric Cross-Genome Analysis
and Dynamic Grouping by Typing (ST, K, OCL, Capsule)

Version: 3.2.0
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Date: 2026-11-01
NEW in 3.2.0:
  - Fixed highlight search (added highlightTableCells JS + CSS)
  - Sticky headers on Sample Overview and FASTA QC tables
  - Virulence gene expander + MobileOG hits column in Sample Overview
  - AMR educational content moved into click-to-expand panels
  - Mutation tab simplified (Mutation and Subclass columns removed)
  - New MGE tab (mobileOG profile, category totals, aggregate by typing)
  - New Compare tab (pair comparison + cluster detection with heatmap)
  - New Call to Action tab
  - Updated AcinetoScope citation (Sci Rep 2026)
"""

import os
import io
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

# =============================================================================
# HTML PARSER CLASS
# =============================================================================
class UltimateHTMLParser:
    """Parse AcinetoScope HTML reports (QC, MLST, Kaptive, AMR, ABRicate, PlasmidFinder, APT, Mutations, MGE)."""

    def __init__(self):
        self.all_databases = [
            'amrfinder', 'card', 'resfinder', 'argannot', 'megares', 'bacmet2',
            'vfdb', 'ecoli_vf', 'plasmidfinder', 'ecoh', 'ncbi'
        ]
        self.db_name_mapping = {
            'acineto_amrfinder': 'amrfinder',
            'acineto_card': 'card',
            'acineto_resfinder': 'resfinder',
            'acineto_argannot': 'argannot',
            'acineto_megares': 'megares',
            'acineto_bacmet2': 'bacmet2',
            'acineto_vfdb': 'vfdb',
            'acineto_ecoli_vf': 'ecoli_vf',
            'acineto_plasmidfinder': 'plasmidfinder',
            'acineto_ecoh': 'ecoh',
            'acineto_ncbi': 'ncbi'
        }

    def normalize_sample_id(self, sample_id: str) -> str:
        """Normalise sample ID by removing path and common extensions."""
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
        """Parse an HTML table from the given content and return a pandas DataFrame."""
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

    # ---- MLST ----
    def parse_mlst_report(self, file_path: Path, scheme: str = "pasteur") -> Dict[str, Dict]:
        print(f"  🧬 Parsing {scheme.upper()} MLST: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            df = self.parse_html_table(html_content, 0)
            if df.empty:
                df = self.parse_html_table(html_content, 1)
                if df.empty:
                    return {}
            df.columns = [col.strip().replace('\n', ' ') for col in df.columns]
            sample_col = None
            for col in df.columns:
                col_lower = col.lower()
                if 'sample' in col_lower or 'genome' in col_lower or 'id' in col_lower or 'strain' in col_lower:
                    if '#' not in col_lower:
                        sample_col = col
                        break
            if not sample_col and len(df.columns) > 1:
                sample_col = df.columns[0]
            if not sample_col:
                return {}
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            results = {}
            valid_st_count = 0
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
                            valid_st_count += 1
                        else:
                            num_match = re.search(r'\d+', st_val)
                            if num_match:
                                st = num_match.group()
                                valid_st_count += 1
                else:
                    for col in df.columns:
                        if col.lower() != sample_col.lower() and pd.notna(row.get(col)):
                            cell_val = str(row[col])
                            match = re.search(r'ST(\d+)', cell_val, re.IGNORECASE)
                            if match:
                                st = match.group(1)
                                valid_st_count += 1
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
                    'ST': st, 'International_Clone': ic,
                    'Allele_Profile': allele_profile, 'Scheme': scheme
                }
            print(f"    ✓ Found {len(results)} samples for {scheme.upper()} scheme (valid STs: {valid_st_count})")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing {scheme} MLST: {e}")
            return {}

    # ---- Kaptive ----
    def parse_kaptive_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing Kaptive: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            df = self.parse_html_table(html_content, 0)
            if df.empty:
                df = self.parse_html_table(html_content, 1)
                if df.empty:
                    return {}
            df.columns = [col.strip() for col in df.columns]
            sample_col = None
            for col in df.columns:
                if any(keyword in col.lower() for keyword in ['genome', 'sample', 'id', 'strain']):
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
                k_col_names = ['K Locus', 'K locus', 'K-Locus', 'K-locus', 'K', 'KLocus']
                for k_col in k_col_names:
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
                o_col_names = ['O Locus', 'O locus', 'O-Locus', 'O-locus', 'O', 'OLocus']
                for o_col in o_col_names:
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
                    'K_Locus': k_locus, 'O_Locus': o_locus,
                    'Capsule_Type': capsule_type
                }
            print(f"    ✓ Found {len(results)} samples with capsule typing")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing Kaptive: {e}")
            return {}

    # ---- QC ----
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
                qc_data = {}
                for col in df.columns:
                    if col == sample_col:
                        continue
                    val = row[col]
                    if pd.isna(val) or val == '' or val == 'ND':
                        qc_data[col] = 'ND'
                    else:
                        cleaned = str(val).replace('%', '').replace(',', '').strip()
                        try:
                            qc_data[col] = float(cleaned)
                        except:
                            qc_data[col] = str(val)
                results[sample] = qc_data
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing QC: {e}")
            return {}

    # ---- AMRfinder ----
    def parse_amrfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        try:
            import pandas as pd
            tables = pd.read_html(file_path)
            if not tables:
                print("    ⚠️ No tables found in HTML")
                return {}, {}

            freq_df = None
            genome_df = None
            for df in tables:
                cols_lower = [str(c).lower() for c in df.columns]
                if 'gene' in cols_lower and 'genomes' in cols_lower:
                    if any('freq' in c or 'prevalence' in c or 'count' in c or 'frequency' in c for c in cols_lower):
                        freq_df = df
                        break
            if freq_df is not None:
                for df in tables:
                    if df is not freq_df:
                        genome_df = df
                        break
            if freq_df is None and len(tables) >= 2:
                genome_df = tables[0]
                freq_df = tables[1]
            elif freq_df is None and len(tables) == 1:
                genome_df = tables[0]
                freq_df = None

            gene_frequencies = {}
            if freq_df is not None and not freq_df.empty:
                cols = freq_df.columns.tolist()
                gene_col = None
                freq_col = None
                genomes_col = None
                for col in cols:
                    col_lower = str(col).lower()
                    if 'gene' in col_lower:
                        gene_col = col
                    elif 'freq' in col_lower or 'prevalence' in col_lower or 'count' in col_lower or 'frequency' in col_lower:
                        freq_col = col
                    elif 'genome' in col_lower:
                        genomes_col = col
                if gene_col is None or genomes_col is None:
                    print("    ⚠️ Could not identify required columns in frequency table")
                    return {}, {}

                if freq_col is None:
                    for col in cols:
                        if freq_df[col].dtype in ['int64', 'float64']:
                            freq_col = col
                            break
                    if freq_col is None:
                        for col in cols:
                            try:
                                pd.to_numeric(freq_df[col])
                                freq_col = col
                                break
                            except:
                                pass

                for _, row in freq_df.iterrows():
                    gene = str(row[gene_col]).strip()
                    if not gene:
                        continue
                    if freq_col is not None:
                        freq_val = str(row[freq_col]).strip()
                        match = re.search(r'(\d+)', freq_val)
                        count = int(match.group(1)) if match else 0
                    else:
                        genomes_str = str(row[genomes_col]).strip()
                        genomes = [g.strip() for g in genomes_str.split(',') if g.strip()]
                        count = len(genomes)
                    genomes_str = str(row[genomes_col]).strip()
                    genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()]
                    gene_frequencies[gene] = {
                        'count': count,
                        'percentage': round((count / total_samples) * 100, 2) if total_samples > 0 else 0,
                        'frequency_display': f"{count} ({count/total_samples*100:.1f}%)" if total_samples > 0 else "0 (0.0%)",
                        'genomes': genomes,
                        'database': 'amrfinder'
                    }

            genes_by_genome = {}
            if genome_df is not None and not genome_df.empty:
                cols = genome_df.columns.tolist()
                genome_col = None
                genes_col = None
                for col in cols:
                    col_lower = str(col).lower()
                    if 'genome' in col_lower or 'sample' in col_lower:
                        genome_col = col
                    elif 'gene' in col_lower or 'detected' in col_lower:
                        genes_col = col
                if genome_col is None or genes_col is None:
                    print("    ⚠️ Could not identify required columns in genome table")
                    return {}, {}

                for _, row in genome_df.iterrows():
                    sample = self.normalize_sample_id(str(row[genome_col]).strip())
                    if not sample:
                        continue
                    gene_str = str(row[genes_col]).strip()
                    genes = [g.strip() for g in re.split(r'[,\s]+', gene_str) if g.strip() and not g.startswith('Showing') and '(' not in g and ')' not in g]
                    genes_by_genome[sample] = genes

            print(f"    ✓ Found {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies

        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder: {e}")
            return {}, {}

    # ---- ABRicate ----
    def parse_abricate_database_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing ABRicate: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}
            db_name = 'unknown'
            filename = str(file_path.name).lower()
            for key, value in self.db_name_mapping.items():
                if key in filename:
                    db_name = value
                    break
            gene_frequencies = {}
            df_freq = self.parse_html_table(str(tables[1]), 0)
            if not df_freq.empty and 'Gene' in df_freq.columns:
                for _, row in df_freq.iterrows():
                    gene_full = str(row['Gene']).strip()
                    if not gene_full:
                        continue
                    gene = re.sub(r'^\([^)]+\)', '', gene_full).strip()
                    if not gene:
                        gene = gene_full
                    count = 0
                    frequency_str = str(row.get('Frequency', '0')).strip()
                    match = re.search(r'(\d+)', frequency_str)
                    if match:
                        count = int(match.group(1))
                    percentage = (count / total_samples) * 100 if total_samples > 0 else 0
                    genomes = []
                    if 'Genomes' in df_freq.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()]
                    gene_frequencies[gene] = {
                        'count': count,
                        'percentage': round(percentage, 2),
                        'frequency_display': f"{count} ({percentage:.1f}%)",
                        'genomes': genomes,
                        'database': db_name,
                        'full_name': gene_full
                    }
            genes_by_genome = {}
            df_genomes = self.parse_html_table(str(tables[0]), 0)
            if not df_genomes.empty:
                sample_col = None
                for col in df_genomes.columns:
                    if any(keyword in col.lower() for keyword in ['genome', 'sample', 'id']):
                        sample_col = col
                        break
                if sample_col:
                    for _, row in df_genomes.iterrows():
                        sample = self.normalize_sample_id(row[sample_col])
                        if not sample:
                            continue
                        genes = []
                        genes_col = None
                        for col in df_genomes.columns:
                            if any(keyword in col.lower() for keyword in ['genes', 'detected']):
                                genes_col = col
                                break
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            genes = [re.sub(r'^\([^)]+\)', '', g).strip() for g in gene_str.split(',') if g.strip()]
                        genes_by_genome[sample] = genes
            print(f"    ✓ {db_name.upper()}: {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
        except Exception as e:
            print(f"    ❌ Error parsing ABRicate report: {e}")
            return {}, {}

    # ---- PlasmidFinder ----
    def parse_plasmidfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing PlasmidFinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}
            plasmid_frequencies = {}
            df_freq = self.parse_html_table(str(tables[1]), 0)
            if not df_freq.empty and 'Gene' in df_freq.columns:
                for _, row in df_freq.iterrows():
                    gene_full = str(row['Gene']).strip()
                    if not gene_full:
                        continue
                    gene = self._clean_plasmid_gene_name(gene_full)
                    count = 0
                    frequency_str = str(row.get('Frequency', '0')).strip()
                    match = re.search(r'(\d+)', frequency_str)
                    if match:
                        count = int(match.group(1))
                    percentage = (count / total_samples) * 100 if total_samples > 0 else 0
                    genomes = []
                    if 'Genomes' in df_freq.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()]
                    plasmid_type = self._categorize_plasmid(gene)
                    plasmid_frequencies[gene] = {
                        'count': count,
                        'percentage': round(percentage, 2),
                        'frequency_display': f"{count} ({percentage:.1f}%)",
                        'genomes': genomes,
                        'full_name': gene_full,
                        'plasmid_type': plasmid_type,
                        'database': 'plasmidfinder'
                    }
            plasmids_by_genome = {}
            df_genomes = self.parse_html_table(str(tables[0]), 0)
            if not df_genomes.empty:
                sample_col = None
                for col in df_genomes.columns:
                    if any(keyword in col.lower() for keyword in ['genome', 'sample', 'id']):
                        sample_col = col
                        break
                if sample_col:
                    for _, row in df_genomes.iterrows():
                        sample = self.normalize_sample_id(row[sample_col])
                        if not sample:
                            continue
                        plasmids = []
                        genes_col = None
                        for col in df_genomes.columns:
                            if any(keyword in col.lower() for keyword in ['genes', 'detected']):
                                genes_col = col
                                break
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            plasmids = [self._clean_plasmid_gene_name(g.strip()) for g in gene_str.split(',') if g.strip()]
                        plasmids_by_genome[sample] = plasmids
            print(f"    ✓ PlasmidFinder: {len(plasmids_by_genome)} samples, {len(plasmid_frequencies)} plasmid markers")
            return plasmids_by_genome, plasmid_frequencies
        except Exception as e:
            print(f"    ❌ Error parsing PlasmidFinder: {e}")
            return {}, {}

    # ---- APT plasmid summary ----
    def parse_apt_plasmid_summary(self, file_path: Path) -> Dict[str, Any]:
        print(f"  🧬 Parsing APT plasmid summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            freq_table = soup.find('table', id='freqTable')
            if not freq_table:
                for table in soup.find_all('table'):
                    if 'Rep Type' in table.get_text():
                        freq_table = table
                        break
            rep_freq = []
            if freq_table:
                rows = freq_table.find_all('tr')
                for row in rows[1:]:
                    cols = row.find_all('td')
                    if len(cols) >= 3:
                        rep_type = cols[0].get_text().strip()
                        count = int(cols[1].get_text().strip())
                        genomes = [g.strip() for g in cols[2].get_text().strip().split(',')]
                        rep_freq.append({
                            'rep_type': rep_type,
                            'count': count,
                            'genomes': genomes
                        })
            genome_table = soup.find('table', id='genomeTable')
            per_genome = {}
            if genome_table:
                rows = genome_table.find_all('tr')
                for row in rows[1:]:
                    cols = row.find_all('td')
                    if len(cols) >= 2:
                        genome = cols[0].get_text().strip()
                        rep_types_raw = cols[1].get_text().strip()
                        rep_types = [rt.strip() for rt in rep_types_raw.split() if rt.strip() and rt.strip() != 'No Rep found']
                        if not rep_types and 'No Rep found' in rep_types_raw:
                            rep_types = []
                        per_genome[genome] = rep_types
            return {
                'rep_frequency': rep_freq,
                'per_genome': per_genome,
                'total_genomes': len(per_genome),
                'total_with_plasmids': sum(1 for reps in per_genome.values() if reps),
                'unique_rep_types': len(rep_freq)
            }
        except Exception as e:
            print(f"    ❌ Error parsing APT plasmid summary: {e}")
            return {}

    # ---- Mutation summary (AMRFinderPlus) ----
    def parse_mutation_summary_html(self, file_path: Path) -> Dict[str, Any]:
        print(f"  🧬 Parsing mutation summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            mutation_table = None
            for table in soup.find_all('table'):
                if table.find(string=re.compile(r'Gene', re.I)) and table.find(string=re.compile(r'Mutation', re.I)):
                    mutation_table = table
                    break
            if not mutation_table:
                print("    ⚠️ Could not find mutation table in HTML")
                return {}
            header_row = None
            thead = mutation_table.find('thead')
            if thead:
                header_row = thead.find('tr')
            if not header_row:
                first_row = mutation_table.find('tr')
                if first_row:
                    header_row = first_row
            if not header_row:
                return {}
            headers = []
            for cell in header_row.find_all(['th', 'td']):
                text = cell.get_text().strip()
                if text:
                    headers.append(text)
            col_idx = {}
            for idx, h in enumerate(headers):
                h_lower = h.lower()
                if 'gene' in h_lower:
                    col_idx['gene'] = idx
                elif 'mutation' in h_lower:
                    col_idx['mutation'] = idx
                elif 'count' in h_lower:
                    col_idx['count'] = idx
                elif 'genome' in h_lower:
                    col_idx['genomes'] = idx
                elif 'class' in h_lower:
                    col_idx['class'] = idx
                elif 'subclass' in h_lower:
                    col_idx['subclass'] = idx
            required = ['gene', 'mutation', 'count', 'genomes']
            for req in required:
                if req not in col_idx:
                    print(f"    ⚠️ Missing required column: {req}. Found headers: {headers}")
                    return {}
            tbody = mutation_table.find('tbody')
            if tbody:
                rows = tbody.find_all('tr')
            else:
                rows = mutation_table.find_all('tr')[1:]
            mutations_list = []
            genome_counts = defaultdict(int)
            for row in rows:
                cells = row.find_all('td')
                if len(cells) <= max(col_idx.values()):
                    continue
                gene = cells[col_idx['gene']].get_text().strip()
                mutation = cells[col_idx['mutation']].get_text().strip()
                count_str = cells[col_idx['count']].get_text().strip()
                count_match = re.search(r'(\d+)', count_str)
                count = int(count_match.group(1)) if count_match else 0
                genomes_str = cells[col_idx['genomes']].get_text().strip()
                genomes = [g.strip() for g in genomes_str.split(',') if g.strip()]
                if not genomes:
                    continue
                for g in genomes:
                    genome_counts[g] += 1
                class_name = cells[col_idx['class']].get_text().strip() if 'class' in col_idx else ''
                subclass = cells[col_idx['subclass']].get_text().strip() if 'subclass' in col_idx else ''
                mutations_list.append({
                    'gene': gene,
                    'mutation': mutation,
                    'class': class_name,
                    'subclass': subclass,
                    'count': count,
                    'genomes': genomes
                })
            mutations_list.sort(key=lambda x: x['count'], reverse=True)
            print(f"    ✓ Parsed {len(mutations_list)} unique mutations across {len(genome_counts)} genomes")
            return {
                'mutations': mutations_list,
                'genome_mutation_counts': dict(genome_counts)
            }
        except Exception as e:
            print(f"    ❌ Error parsing mutation summary: {e}")
            return {}

    # ---- MGE summary (NEW in v3.2.0) ----
    def parse_mge_summary_html(self, file_path: Path) -> Dict[str, Any]:
        """Parse the AcinetoScope MGE summary HTML (stat cards + per-sample table)."""
        print(f"  🧬 Parsing MGE summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            stats = {}
            for label in ('Total Samples', 'With mobileOG hits', 'Without mobileOG hits',
                          'Total mobileOG hits', 'Key MGE-associated genes', 'Runtime'):
                for tag in soup.find_all(string=re.compile(re.escape(label), re.I)):
                    parent = tag.find_parent()
                    text_block = parent.get_text(" ", strip=True) if parent else ''
                    m = re.search(rf'{re.escape(label)}\s*[:]?\s*([\d,\.]+s?)\b', text_block, re.I)
                    if m:
                        stats[label] = m.group(1)
                        break
                    if parent:
                        for sib in parent.find_all_next(limit=3):
                            s = sib.get_text(strip=True)
                            if re.fullmatch(r'[\d,\.]+s?', s):
                                stats[label] = s
                                break
                    if label in stats:
                        break
            per_sample = {}
            mge_table = None
            for table in soup.find_all('table'):
                header_text = ' '.join(th.get_text().strip() for th in table.find_all('th'))
                if 'Sample' in header_text and 'mobileOG' in header_text:
                    mge_table = table
                    break
            if mge_table:
                try:
                    df = pd.read_html(io.StringIO(str(mge_table)))[0]
                    df.columns = [str(c).strip() for c in df.columns]
                    for _, row in df.iterrows():
                        sample = self.normalize_sample_id(row.get('Sample', ''))
                        if not sample:
                            continue
                        entry = {}
                        for col in df.columns:
                            if col == 'Sample':
                                continue
                            v = row.get(col, 0)
                            try:
                                entry[col] = int(str(v).replace(',', '').strip())
                            except Exception:
                                entry[col] = v
                        per_sample[sample] = entry
                except Exception as e:
                    print(f"    Error parsing MGE table: {e}")
            note = ""
            for tag in soup.find_all(string=re.compile(r'Interpretation note', re.I)):
                parent = tag.find_parent()
                if parent:
                    raw = re.sub(r'^Interpretation note\s*:?\s*', '',
                                 parent.get_text(" ", strip=True), flags=re.I)
                    if len(raw) > 40:
                        note = raw
                        break
            print(f"    ✓ MGE: {len(stats)} stat cards, {len(per_sample)} samples")
            return {'stats': stats, 'per_sample': per_sample, 'note': note}
        except Exception as e:
            print(f"    ❌ Error parsing MGE summary: {e}")
            return {'stats': {}, 'per_sample': {}, 'note': ''}

    # ---- helpers ----
    def _clean_plasmid_gene_name(self, gene_name: str) -> str:
        gene = re.sub(r'_(\d+)$', '', gene_name)
        plasmid_match = re.search(r'\((.*?)\)', gene)
        if plasmid_match:
            plasmid_name = plasmid_match.group(1)
            base_gene = re.sub(r'\(.*?\)', '', gene).strip()
            if base_gene:
                gene = f"{base_gene}({plasmid_name})"
            else:
                gene = plasmid_name
        return gene.strip()

    def _categorize_plasmid(self, gene_name: str) -> str:
        gene_lower = gene_name.lower()
        if any(rep in gene_lower for rep in ['rep', 'inc', 'rep_']):
            if 'coli' in gene_lower or 'col' in gene_lower:
                return 'Colicin plasmid'
            elif 'broad' in gene_lower or 'multihost' in gene_lower:
                return 'Broad-host-range plasmid'
            else:
                return 'Replication protein'
        elif any(mob in gene_lower for mob in ['mob', 'tra', 'conj']):
            return 'Mobility/conjugation'
        elif 'col' in gene_lower and 'coli' not in gene_lower:
            if gene_lower.startswith('col'):
                return 'Colicin plasmid'
            else:
                return 'Other plasmid'
        elif 'inc' in gene_lower and gene_lower.startswith('inc'):
            inc_match = re.search(r'inc([a-z0-9]+)', gene_lower)
            if inc_match:
                inc_type = inc_match.group(1).upper()
                return f'Incompatibility group {inc_type}'
            return 'Incompatibility group'
        elif any(fam in gene_lower for fam in ['acinetobacter', 'baumannii']):
            return 'Acinetobacter plasmid'
        else:
            return 'Other plasmid'


# =============================================================================
# DATA ANALYZER CLASS 
# =============================================================================
class UltimateDataAnalyzer:
    """Analyze integrated data: gene categorization, cross‑genome patterns, co‑occurrence, and sample typing map."""

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
        self.critical_esbls = {'blaCTX-M', 'blaSHV', 'blaPER', 'blaVEB', 'blaBEL', 'blaGES'}
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

    def create_gene_centric_tables(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        gene_centric = {
            'amr_databases': {}, 'virulence_databases': {}, 'bacmet_databases': {},
            'plasmid_databases': {}, 'combined_gene_frequencies': [],
            'gene_categories': defaultdict(list), 'database_stats': {}
        }
        if 'amrfinder' in integrated_data.get('gene_frequencies', {}):
            amr_data = integrated_data['gene_frequencies']['amrfinder']
            gene_list = []
            for gene, data in amr_data.items():
                category = self.categorize_gene(gene)
                gene_list.append({
                    'gene': gene, 'category': category, 'database': 'AMRfinder',
                    'count': data.get('count', 0), 'percentage': data.get('percentage', 0),
                    'frequency_display': data.get('frequency_display', f"{data.get('count', 0)} ({data.get('percentage', 0):.1f}%)"),
                    'genomes': data.get('genomes', [])
                })
            if gene_list:
                gene_list.sort(key=lambda x: x['count'], reverse=True)
                gene_centric['amr_databases']['amrfinder'] = gene_list
                for gd in gene_list:
                    gene_centric['gene_categories'][gd['category']].append(gd)
        if 'abricate' in integrated_data.get('gene_frequencies', {}):
            abricate_data = integrated_data['gene_frequencies']['abricate']
            for db_name, db_genes in abricate_data.items():
                gene_list = []
                for gene, data in db_genes.items():
                    category = self.categorize_gene(gene)
                    gene_list.append({
                        'gene': gene, 'category': category, 'database': db_name.upper(),
                        'count': data.get('count', 0), 'percentage': data.get('percentage', 0),
                        'frequency_display': data.get('frequency_display', f"{data.get('count', 0)} ({data.get('percentage', 0):.1f}%)"),
                        'genomes': data.get('genomes', []),
                        'full_name': data.get('full_name', gene)
                    })
                if gene_list:
                    gene_list.sort(key=lambda x: x['count'], reverse=True)
                    db_lower = db_name.lower()
                    if 'vfdb' in db_lower or 'ecoli_vf' in db_lower:
                        gene_centric['virulence_databases'][db_name] = gene_list
                    elif 'bacmet2' in db_lower:
                        gene_centric['bacmet_databases'][db_name] = gene_list
                    elif 'plasmidfinder' in db_lower or 'ecoh' in db_lower:
                        gene_centric['plasmid_databases'][db_name] = gene_list
                    else:
                        gene_centric['amr_databases'][db_name] = gene_list
                    for gd in gene_list:
                        gene_centric['gene_categories'][gd['category']].append(gd)
        all_genes = []
        for db_type in ['amr_databases', 'virulence_databases', 'bacmet_databases', 'plasmid_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    all_genes.append(gene_data)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = all_genes
        for cat in gene_centric['gene_categories']:
            gene_centric['gene_categories'][cat].sort(key=lambda x: x['count'], reverse=True)
        for db_type in ['amr_databases', 'virulence_databases', 'bacmet_databases', 'plasmid_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                gene_centric['database_stats'][db_name] = {
                    'total_genes': len(genes),
                    'total_occurrences': sum(g['count'] for g in genes),
                    'critical_genes': sum(1 for g in genes if g['category'] in ['Carbapenemases', 'ESBLs', 'Colistin Resistance', 'Tigecycline Resistance']),
                }
        return gene_centric

    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        patterns = {
            'pasteur_st_distribution': Counter(), 'oxford_st_distribution': Counter(),
            'k_locus_distribution': Counter(), 'o_locus_distribution': Counter(),
            'capsule_type_distribution': Counter(),
            'st_k_combinations': defaultdict(list), 'st_o_combinations': defaultdict(list),
            'ko_combinations': defaultdict(list), 'st_ko_combinations': defaultdict(list),
            'carbapenemase_patterns': defaultdict(list),
            'high_risk_crab': [], 'gene_cooccurrence': defaultdict(Counter)
        }
        samples_data = integrated_data.get('samples', {})
        gene_centric = integrated_data.get('gene_centric', {})
        sample_genes = defaultdict(set)
        for db_type in ['amr_databases', 'virulence_databases', 'bacmet_databases', 'plasmid_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    for genome in gene_data['genomes']:
                        sample_genes[genome].add(gene_data['gene'])
        cooc = patterns['gene_cooccurrence']
        for genes in sample_genes.values():
            gene_list = list(genes)
            for i, g1 in enumerate(gene_list):
                for g2 in gene_list[i+1:]:
                    cooc[g1][g2] += 1

        for sample, data in samples_data.items():
            pasteur_st = data.get('pasteur_mlst', {}).get('ST', 'ND')
            oxford_st = data.get('oxford_mlst', {}).get('ST', 'ND')
            k_locus = data.get('kaptive', {}).get('K_Locus', 'ND')
            o_locus = data.get('kaptive', {}).get('O_Locus', 'ND')
            capsule = data.get('kaptive', {}).get('Capsule_Type', 'ND')
            if pasteur_st != 'ND':
                patterns['pasteur_st_distribution'][pasteur_st] += 1
            if oxford_st != 'ND':
                patterns['oxford_st_distribution'][oxford_st] += 1
            if k_locus != 'ND':
                patterns['k_locus_distribution'][k_locus] += 1
            if o_locus != 'ND':
                patterns['o_locus_distribution'][o_locus] += 1
            if capsule != 'ND':
                patterns['capsule_type_distribution'][capsule] += 1
            if pasteur_st != 'ND' and o_locus != 'ND':
                patterns['st_o_combinations'][f"ST{pasteur_st} - {o_locus}"].append(sample)
            if pasteur_st != 'ND' and k_locus != 'ND':
                patterns['st_k_combinations'][f"ST{pasteur_st} - {k_locus}"].append(sample)
            if k_locus != 'ND' and o_locus != 'ND':
                patterns['ko_combinations'][f"{k_locus}:{o_locus}"].append(sample)
            if pasteur_st != 'ND' and k_locus != 'ND' and o_locus != 'ND':
                patterns['st_ko_combinations'][f"ST{pasteur_st} - {k_locus}:{o_locus}"].append(sample)
            genes = list(sample_genes.get(sample, set()))
            carbapenemases = [g for g in genes if 'Carbapenemases' in self.categorize_gene(g)]
            colistin_res = [g for g in genes if 'Colistin Resistance' in self.categorize_gene(g)]
            tigecycline_res = [g for g in genes if 'Tigecycline Resistance' in self.categorize_gene(g)]
            if carbapenemases:
                key = '|'.join(sorted(carbapenemases))
                patterns['carbapenemase_patterns'][key].append(sample)
            if carbapenemases and (colistin_res or tigecycline_res):
                patterns['high_risk_crab'].append({
                    'sample': sample, 'pasteur_st': pasteur_st,
                    'k_locus': k_locus, 'capsule_type': capsule,
                    'carbapenemases': carbapenemases,
                    'colistin_resistance': colistin_res,
                    'tigecycline_resistance': tigecycline_res
                })
        return patterns

    def create_plasmid_analysis(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        plasmid_analysis = {'plasmid_databases': {}, 'plasmid_frequencies': [], 'sample_plasmid_profiles': defaultdict(list)}
        if 'plasmidfinder' not in integrated_data.get('gene_frequencies', {}):
            return plasmid_analysis
        plasmid_genes = integrated_data['gene_frequencies']['plasmidfinder']
        gene_list = []
        for gene, data in plasmid_genes.items():
            gene_list.append({
                'plasmid_marker': gene, 'category': data.get('plasmid_type', 'Other'),
                'database': 'PlasmidFinder', 'count': data.get('count', 0),
                'percentage': data.get('percentage', 0),
                'frequency_display': data.get('frequency_display', f"{data.get('count', 0)} ({data.get('percentage', 0):.1f}%)"),
                'genomes': data.get('genomes', [])
            })
        if gene_list:
            gene_list.sort(key=lambda x: x['count'], reverse=True)
            plasmid_analysis['plasmid_databases']['plasmidfinder'] = gene_list
            plasmid_analysis['plasmid_frequencies'] = gene_list
        samples_data = integrated_data.get('samples', {})
        for sample in samples_data:
            if 'plasmidfinder' in integrated_data.get('gene_frequencies', {}):
                db_genes = integrated_data['gene_frequencies']['plasmidfinder']
                sample_plasmids = []
                for gene, gdata in db_genes.items():
                    if sample in gdata.get('genomes', []):
                        sample_plasmids.append(gene)
                if sample_plasmids:
                    plasmid_analysis['sample_plasmid_profiles'][sample] = sample_plasmids
        return plasmid_analysis

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


# =============================================================================
# HTML GENERATOR CLASS
# =============================================================================
class UltimateHTMLGenerator:
    """Generate the interactive HTML report with all tabs and dynamic features."""

    def __init__(self, data_analyzer: UltimateDataAnalyzer):
        self.data_analyzer = data_analyzer
        self._current_samples_data = {}
        self.tab_colors = {
            'summary': '#4CAF50', 'samples': '#2196F3', 'qc': '#607D8B',
            'mlst': '#FF9800', 'kaptive': '#9C27B0', 'combinations': '#009688',
            'amr': '#F44336', 'virulence': '#E91E63', 'bacmet': '#795548',
            'plasmid': '#2196F3', 'mutation': '#00BCD4', 'mge': '#16A085',
            'cooccurrence': '#3F51B5', 'patterns': '#FF5722',
            'compare': '#0d9488',
            'aiguide': '#3F51B5', 'calltoaction': '#2E7D32',
            'citation': '#8BC34A', 'funding': '#FFC107', 'export': '#9E9E9E',
        }

    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        print("\n🎨 Generating ULTIMATE HTML report for A. baumannii...")
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        plasmid_analysis = integrated_data.get('plasmid_analysis', {})
        qc_data = integrated_data.get('qc_data', {})
        apt_plasmid = integrated_data.get('apt_plasmid', {})
        mutation_data = integrated_data.get('mutation_data', {})
        mge_data = integrated_data.get('mge_data', {})
        typing_map = self.data_analyzer.build_sample_typing_map(samples_data)
        self._current_samples_data = samples_data
        html = self._create_ultimate_html(
            metadata=metadata, samples_data=samples_data, patterns=patterns,
            gene_centric=gene_centric, plasmid_analysis=plasmid_analysis,
            qc_data=qc_data, total_samples=len(samples_data), apt_plasmid=apt_plasmid,
            mutation_data=mutation_data, mge_data=mge_data, typing_map=typing_map
        )
        output_file = output_dir / "genius_acinetobacter_ultimate_gene_centric_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)

    def _create_ultimate_html(self, **kwargs) -> str:
        metadata = kwargs['metadata']
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        plasmid_analysis = kwargs.get('plasmid_analysis', {})
        qc_data = kwargs.get('qc_data', {})
        total_samples = kwargs['total_samples']
        apt_plasmid = kwargs.get('apt_plasmid', {})
        mutation_data = kwargs.get('mutation_data', {})
        mge_data = kwargs.get('mge_data', {})
        typing_map = kwargs.get('typing_map', {})

        total_amr = sum(len(db) for db in gene_centric.get('amr_databases', {}).values())
        total_vir = sum(len(db) for db in gene_centric.get('virulence_databases', {}).values())
        total_bacmet = sum(len(db) for db in gene_centric.get('bacmet_databases', {}).values())
        if apt_plasmid and apt_plasmid.get('unique_rep_types', 0) > 0:
            total_plasmid = apt_plasmid['unique_rep_types']
        else:
            total_plasmid = sum(len(db) for db in plasmid_analysis.get('plasmid_databases', {}).values())
        high_risk = len(patterns.get('high_risk_crab', []))
        carbapenemase_count = sum(1 for db in gene_centric.get('amr_databases', {}).values()
                                  for g in db if g['category'] == 'Carbapenemases')
        env_count = total_bacmet
        mutation_count = len(mutation_data.get('mutations', []))
        mge_count = len(mge_data.get('per_sample', {}))

        typing_json = json.dumps(typing_map, default=str).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
        css = self._get_css()
        js = self._get_js(typing_json)

        # TABS — new: mge, compare, calltoaction
        TABS = [
            ('summary', 'Summary', 'chart-pie'),
            ('samples', 'Sample Overview', 'list-alt'),
            ('qc', 'FASTA QC', 'chart-line'),
            ('mlst', 'MLST', 'code-branch'),
            ('kaptive', 'Capsule', 'shield-alt'),
            ('combinations', 'Combinations', 'link'),
            ('amr', 'AMR', 'biohazard'),
            ('virulence', 'Virulence', 'virus'),
            ('bacmet', 'Bacmet', 'flask'),
            ('plasmid', 'Plasmids', 'dna'),
            ('mutation', 'Mutations', 'dna'),
            ('mge', 'MGE', 'mobile-alt'),
            ('cooccurrence', 'Co-occurrence', 'link'),
            ('patterns', 'Patterns', 'project-diagram'),
            ('compare', 'Compare', 'balance-scale'),
            ('aiguide', 'AI Guide', 'robot'),
            ('calltoaction', 'Call to Action', 'globe'),
            ('citation', 'Citations', 'book'),
            ('funding', 'Funding', 'coffee'),
            ('export', 'Export', 'download'),
        ]

        nav_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            nav_html += f'<button class="tab-button {tid}{active}" onclick="switchTab(\'{tid}\')"><i class="fas fa-{icon}"></i> {title}</button>'

        section_methods = {
            'summary': lambda: self._summary_section(kwargs, carbapenemase_count, env_count),
            'samples': lambda: self._samples_section(kwargs),
            'qc': lambda: self._qc_section(kwargs),
            'mlst': lambda: self._mlst_section(kwargs),
            'kaptive': lambda: self._kaptive_section(kwargs),
            'combinations': lambda: self._combinations_section(kwargs),
            'amr': lambda: self._amr_section(kwargs),
            'virulence': lambda: self._virulence_section(kwargs),
            'bacmet': lambda: self._bacmet_section(kwargs),
            'plasmid': lambda: self._plasmid_section(kwargs, apt_plasmid),
            'mutation': lambda: self._mutation_section(kwargs),
            'mge': lambda: self._mge_section(kwargs),
            'cooccurrence': lambda: self._cooccurrence_section(kwargs),
            'patterns': lambda: self._patterns_section(kwargs),
            'compare': lambda: self._compare_section(kwargs),
            'aiguide': lambda: self._aiguide_section(),
            'calltoaction': lambda: self._calltoaction_section(),
            'citation': lambda: self._citation_section(),
            'funding': lambda: self._funding_section(),
            'export': lambda: self._export_section(),
        }

        tabs_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            content = section_methods[tid]()
            tabs_html += f'''
            <div id="{tid}-tab" class="tab-content{active}">
                <h2 class="section-header {tid}-header">
                    <span><i class="fas fa-{icon}"></i> {title}</span>
                    <button class="print-section-btn" onclick="printSection('{tid}-tab')">
                        <i class="fas fa-print"></i> Print
                    </button>
                </h2>
                {content}
            </div>'''

        return f'''<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>GENIUS Acinetobacter baumannii Ultimate Report v3.2.0</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
{css}{js}</head>
<body><div class="container">
<div class="main-header"><h1><i class="fas fa-bacterium"></i> GENIUS Acinetobacter baumannii Ultimate Analysis Report</h1>
<p>Gene-Centric Cross-Genome Analysis with Dynamic Grouping by Typing</p>
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
<div class="dashboard-card card-plasmid" onclick="switchTab('plasmid')"><div class="card-number">{total_plasmid}</div><div class="card-label">Plasmid Rep Types</div><i class="fas fa-dna fa-2x" style="color:var(--plasmid-color)"></i></div>
<div class="dashboard-card card-mge" onclick="switchTab('mge')"><div class="card-number">{mge_count}</div><div class="card-label">MGE-profiled</div><i class="fas fa-mobile-alt fa-2x" style="color:var(--mge-color)"></i></div>
<div class="dashboard-card card-patterns" onclick="switchTab('patterns')"><div class="card-number">{high_risk}</div><div class="card-label">CRAB High-Risk</div><i class="fas fa-exclamation-triangle fa-2x" style="color:var(--patterns-color)"></i></div>
</div>
<div class="tab-navigation">{nav_html}</div>
{tabs_html}
<div class="footer"><h3>GENIUS Acinetobacter baumannii Ultimate Reporter v3.2.0</h3><p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p><p>If you find this tool helpful, please <a href="https://github.com/bbeckley-hub/acinetoscope" target="_blank">⭐ star us on GitHub</a> and share with your network.</p><p>Critical Genes Tracked: Carbapenemases • ESBLs • Colistin Resistance • Tigecycline Resistance • Biofilm Formation • Efflux Pumps • Environmental Co-Selection</p><p>Generated on {metadata.get('analysis_date','Unknown')}</p></div>
</div></body></html>'''

    # --------------------------------------------------------------------------
    # CSS — original + minimal additions for v3.2.0
    # --------------------------------------------------------------------------
    def _get_css(self) -> str:
        return """
        <style>
        :root { --summary-color:#4CAF50; --samples-color:#2196F3; --mlst-color:#FF9800; --kaptive-color:#9C27B0; --amr-color:#F44336; --virulence-color:#E91E63; --bacmet-color:#795548; --combinations-color:#009688; --patterns-color:#FF5722; --plasmid-color:#2196F3; --databases-color:#607D8B; --qc-color:#17a2b8; --aiguide-color:#3F51B5; --calltoaction-color:#2E7D32; --export-color:#3F51B5; --mutation-color:#00BCD4; --cooccurrence-color:#3F51B5; --citation-color:#8BC34A; --funding-color:#FFC107; --mge-color:#16A085; --compare-color:#0d9488; }
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
        .card-label{font-size:.9em;color:#555;font-weight:600}
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
        .tab-button.plasmid.active{background:var(--plasmid-color)}
        .tab-button.mutation.active{background:var(--mutation-color)}
        .tab-button.mge.active{background:var(--mge-color)}
        .tab-button.cooccurrence.active{background:var(--cooccurrence-color)}
        .tab-button.patterns.active{background:var(--patterns-color)}
        .tab-button.compare.active{background:var(--compare-color)}
        .tab-button.aiguide.active{background:var(--aiguide-color)}
        .tab-button.calltoaction.active{background:var(--calltoaction-color)}
        .tab-button.citation.active{background:var(--citation-color)}
        .tab-button.funding.active{background:var(--funding-color)}
        .tab-button.export.active{background:var(--export-color)}
        .tab-content{display:none;background:white;padding:30px;border-radius:15px;box-shadow:0 10px 30px rgba(0,0,0,0.1);margin-bottom:30px;animation:fadeIn 0.5s ease;width:100%;overflow-x:auto}
        .tab-content.active{display:block}
        @keyframes fadeIn{from{opacity:0;transform:translateY(20px)}to{opacity:1;transform:translateY(0)}}
        .section-header{color:#2c3e50;margin-bottom:25px;padding-bottom:15px;border-bottom:3px solid;font-size:1.8em;display:flex;align-items:center;justify-content:space-between}
        .summary-header{border-color:var(--summary-color)}
        .samples-header{border-color:var(--samples-color)}
        .qc-header{border-color:var(--qc-color)}
        .mlst-header{border-color:var(--mlst-color)}
        .kaptive-header{border-color:var(--kaptive-color)}
        .combinations-header{border-color:var(--combinations-color)}
        .amr-header{border-color:var(--amr-color)}
        .virulence-header{border-color:var(--virulence-color)}
        .bacmet-header{border-color:var(--bacmet-color)}
        .plasmid-header{border-color:var(--plasmid-color)}
        .mutation-header{border-color:var(--mutation-color)}
        .mge-header{border-color:var(--mge-color)}
        .cooccurrence-header{border-color:var(--cooccurrence-color)}
        .patterns-header{border-color:var(--patterns-color)}
        .compare-header{border-color:var(--compare-color)}
        .aiguide-header{border-color:var(--aiguide-color)}
        .calltoaction-header{border-color:var(--calltoaction-color)}
        .citation-header{border-color:var(--citation-color)}
        .funding-header{border-color:var(--funding-color)}
        .export-header{border-color:var(--export-color)}
        .print-section-btn{background:#00695c;color:white;border:none;border-radius:5px;padding:8px 15px;cursor:pointer;display:flex;align-items:center;gap:5px;font-size:0.9em}
        .print-section-btn:hover{background:#004d40}
        .data-table{width:100%;border-collapse:collapse;margin:20px 0;font-size:0.95em;box-shadow:0 2px 10px rgba(0,0,0,0.1);border-radius:8px;overflow:hidden;table-layout:auto}
        .data-table th{background:#2c3e50;color:white;padding:15px;text-align:left;font-weight:600;white-space:nowrap;cursor:pointer}
        .data-table th:hover{background:#1a252f}
        .data-table td{padding:12px;border-bottom:1px solid #e0e0e0;vertical-align:top;white-space:nowrap}
        .data-table tr:hover{background:#f8f9fa}
        .master-scrollable-container{width:100%;overflow-x:auto;border:1px solid #e0e0e0;border-radius:8px;margin:20px 0}
        /* v3.2.0: sticky headers when vertical scrolling */
        .master-scrollable-container.sticky-y{max-height:70vh;overflow-y:auto}
        .master-scrollable-container.sticky-y .data-table th{position:sticky;top:0;z-index:10}
        .genome-list{display:flex;flex-wrap:wrap;gap:5px;max-height:200px;overflow-y:auto;padding:5px;background:#f8f9fa;border-radius:5px}
        .genome-group{margin-bottom:10px;width:100%}
        .genome-group-header{font-weight:bold;background:#e0e0e0;padding:4px 8px;border-radius:4px;margin:5px 0;font-size:0.85em;display:inline-block}
        .genome-group-tags{display:flex;flex-wrap:wrap;gap:5px;margin-left:10px}
        .genome-tag{display:inline-block;background:#e6ffe6;color:#006400;padding:3px 10px;border-radius:12px;font-size:0.85em;border:1px solid #b3ffb3;white-space:nowrap;margin:2px}
        .genome-tag.highlight{background-color:#ffff99 !important;color:#000 !important;border:1px solid #ffc107}
        /* v3.2.0: cell-highlight for the highlight search box */
        .cell-highlight{background-color:#ffff99 !important;color:#000 !important}
        .search-box{width:100%;padding:12px;margin-bottom:15px;border:2px solid #e0e0e0;border-radius:8px;font-size:1em;transition:all 0.3s ease}
        .search-box:focus{outline:none;border-color:#00695c;box-shadow:0 0 0 3px rgba(0,105,92,0.1)}
        .badge{display:inline-block;padding:5px 15px;border-radius:20px;font-size:0.85em;font-weight:600;margin:2px}
        .badge-critical{background:#9C27B0;color:white}
        .badge-high{background:#F44336;color:white}
        .badge-medium{background:#FF9800;color:black}
        .badge-low{background:#9E9E9E;color:white}
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
        /* v3.2.0: expandable educational sections */
        .edu-details{margin:12px 0}
        .edu-details summary{cursor:pointer;padding:12px 16px;background:#f0f7f0;border-left:4px solid #17a2b8;border-radius:6px;font-weight:600;color:#2c3e50;list-style:none;display:flex;justify-content:space-between;align-items:center}
        .edu-details summary::-webkit-details-marker{display:none}
        .edu-details summary::before{content:'▸ ';display:inline-block;margin-right:6px;transition:transform .2s}
        .edu-details[open] summary::before{content:'▾ '}
        .edu-details > div{padding:10px 0 16px 0}
        /* v3.2.0: expandable virulence list in Sample Overview */
        .vir-details{margin-top:4px}
        .vir-details summary{cursor:pointer;color:#006400;font-weight:600;font-size:.9em;list-style:none}
        .vir-details summary::-webkit-details-marker{display:none}
        .vir-details summary::before{content:'▸ ';display:inline-block;margin-right:4px}
        .vir-details[open] summary::before{content:'▾ '}
        .vir-details .vir-gene-tags{margin-top:6px;display:flex;flex-wrap:wrap;gap:4px;max-width:600px}
        /* v3.2.0: MGE stat cards */
        .stat-card{color:white;padding:18px;border-radius:10px;text-align:center;box-shadow:0 4px 15px rgba(0,0,0,.15);transition:transform 0.2s}
        .stat-card:hover{transform:translateY(-3px)}
        .stat-card .stat-value{font-size:1.9em;font-weight:bold;margin-bottom:4px}
        .stat-card .stat-label{font-size:.85em;opacity:.95;text-transform:uppercase;letter-spacing:.5px}
        .stats-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:15px;margin:20px 0}
        /* v3.2.0: Compare tab */
        .mode-toggle{display:flex;gap:8px;margin-bottom:18px}
        .mode-btn{padding:10px 22px;border:2px solid #cbd5e1;background:white;color:#475569;border-radius:8px;font-weight:600;cursor:pointer;font-size:.95em;display:flex;align-items:center;gap:8px;transition:all .2s}
        .mode-btn:hover{border-color:#006400;color:#006400}
        .mode-btn.active{background:#006400;border-color:#006400;color:white}
        .compare-mode{display:none}.compare-mode.active{display:block}
        .compare-visuals-grid{display:grid;grid-template-columns:180px 1fr;gap:16px;margin:20px 0;align-items:center}
        .visual-card{background:white;border-radius:12px;padding:20px;box-shadow:0 2px 8px rgba(0,0,0,.06);border:1px solid #e2e8f0}
        .gauge-card{display:flex;justify-content:center;align-items:center;padding:12px}
        .similarity-gauge{width:160px;height:160px;display:block}
        .metrics-card{display:flex;flex-direction:column;gap:18px}
        .metric-block{display:flex;flex-direction:column;gap:6px}
        .metric-label{font-size:.85em;font-weight:700;color:#475569;text-transform:uppercase;letter-spacing:.5px}
        .metric-bar-track{background:#f1f5f9;height:22px;border-radius:11px;overflow:hidden;position:relative}
        .metric-bar-fill{height:100%;border-radius:11px;transition:width .8s ease;min-width:4px}
        .metric-value{font-size:.88em;color:#334155;font-weight:600}
        .typing-pill{display:inline-block;padding:5px 12px;border-radius:12px;font-size:.85em;font-weight:600;border:1.5px solid;font-family:'Consolas',monospace}
        .gene-proportion-bar{display:flex;height:12px;border-radius:6px;overflow:hidden;margin-bottom:14px;background:#f1f5f9}
        .gene-prop-segment{height:100%;transition:width .6s ease;min-width:2px}
        .gene-prop-segment.shared{background:#10b981}
        .gene-prop-segment.only-a{background:#3b82f6}
        .gene-prop-segment.only-b{background:#ec4899}
        #cluster-clusters{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:14px;margin-bottom:24px}
        .cluster-card{background:white;border-radius:10px;padding:14px 18px;border-left:5px solid;box-shadow:0 2px 8px rgba(0,0,0,.06)}
        .cluster-header{display:flex;align-items:center;gap:14px;margin-bottom:12px;flex-wrap:wrap}
        .cluster-badge{color:white;padding:4px 14px;border-radius:14px;font-size:.85em;font-weight:700}
        .cluster-stats{font-size:.9em;color:#64748b;font-weight:600}
        .cluster-samples{display:flex;flex-wrap:wrap;gap:6px}
        .cluster-sample-pill{background:#f1f5f9;color:#1e293b;padding:4px 12px;border-radius:12px;font-size:.85em;font-weight:600;font-family:'Consolas',monospace;border:1px solid #cbd5e1}
        .cluster-controls{display:flex;gap:14px;align-items:flex-end;flex-wrap:wrap;margin-bottom:20px}
        .cluster-controls label{display:flex;flex-direction:column;gap:4px}
        .cluster-controls label span{font-size:.85em;font-weight:600;color:#475569}
        .cluster-controls select{padding:10px 14px;border:2px solid #cbd5e1;border-radius:8px;background:white;font-size:.95em;min-width:220px}
        .heatmap-wrapper{background:white;border-radius:12px;padding:20px;margin-top:24px;box-shadow:0 2px 8px rgba(0,0,0,.06)}
        .heatmap-grid{display:grid;gap:2px;margin-top:16px;padding:4px}
        .heatmap-cell{width:32px;height:32px;display:flex;align-items:center;justify-content:center;font-size:.7em;font-weight:700;border-radius:3px;cursor:default;overflow:hidden}
        .heatmap-col-header,.heatmap-row-header{background:transparent;color:#64748b;font-size:.75em}
        .heatmap-row-header{justify-content:flex-start;padding-left:6px;font-family:'Consolas',monospace}
        .heatmap-col-header{writing-mode:vertical-rl;text-orientation:mixed;height:110px;width:32px;padding:4px 0;font-family:'Consolas',monospace}
        .heatmap-corner{background:transparent}
        .compare-panel{background:#f8f9fa;border-radius:12px;padding:20px;margin-bottom:24px}
        .compare-select-row{display:flex;gap:12px;align-items:flex-end;flex-wrap:wrap;margin-bottom:16px}
        .compare-select-row label{display:flex;flex-direction:column;gap:4px}
        .compare-select-row label span{font-size:.85em;font-weight:600;color:#475569}
        .compare-select-row select{padding:10px 14px;border:2px solid #cbd5e1;border-radius:8px;background:white;font-size:.95em;min-width:220px}
        .compare-verdict-box{padding:16px 20px;border-radius:10px;margin-bottom:16px;display:flex;gap:12px;align-items:flex-start;font-size:.95em}
        .compare-verdict-box.danger{background:#fef2f2;border-left:5px solid #dc2626;color:#7f1d1d}
        .compare-verdict-box.warning{background:#fffbeb;border-left:5px solid #f59e0b;color:#78350f}
        .compare-verdict-box.info{background:#eff6ff;border-left:5px solid #2563eb;color:#1e3a8a}
        .compare-warning{padding:20px;background:#fef2f2;color:#7f1d1d;border-radius:8px;border-left:4px solid #dc2626}
        .compare-section{background:white;border-radius:12px;padding:20px;margin-bottom:20px;box-shadow:0 2px 8px rgba(0,0,0,.06)}
        .compare-section h3{margin-bottom:14px;color:#1e293b;font-size:1.15em;padding-bottom:8px;border-bottom:2px solid #e2e8f0}
        .compare-table{display:flex;flex-direction:column;gap:2px}
        .compare-row{display:grid;grid-template-columns:180px 1fr 1fr;gap:12px;padding:10px 14px;border-radius:6px;align-items:center;font-size:.92em}
        .compare-row.header{background:#f1f5f9;font-weight:700;border-bottom:2px solid #cbd5e1}
        .compare-row.diff{background:#fef2f2;border-left:3px solid #dc2626}
        .compare-row.match{background:#f0fdf4;border-left:3px solid #16a34a}
        .compare-row .row-label{font-weight:600;color:#475569;font-size:.88em;text-transform:uppercase;letter-spacing:.3px}
        .compare-row .row-val{color:#1e293b;word-break:break-word}
        .compare-gene-block{margin-bottom:24px;padding-bottom:20px;border-bottom:1px dashed #e2e8f0}
        .compare-gene-block:last-child{border-bottom:none;margin-bottom:0;padding-bottom:0}
        .compare-gene-block h4{margin-bottom:12px;color:#1e293b;font-size:1.02em;display:flex;align-items:center;gap:10px;flex-wrap:wrap}
        .gene-count{font-size:.82em;font-weight:400;color:#64748b;background:#f1f5f9;padding:3px 10px;border-radius:12px}
        .gene-columns{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px}
        .gene-col{background:#f8fafc;border-radius:8px;padding:12px;border:1px solid #e2e8f0}
        .gene-col-head{font-size:.85em;font-weight:700;color:#475569;margin-bottom:8px;text-transform:uppercase;letter-spacing:.3px}
        .compare-gene{display:inline-block;font-size:.82em;padding:3px 10px;border-radius:12px;margin:2px;font-family:'Consolas',monospace}
        .compare-gene.shared{background:#dcfce7;color:#14532d;border:1px solid #86efac}
        .compare-gene.only-a{background:#dbeafe;color:#1e3a8a;border:1px solid #93c5fd}
        .compare-gene.only-b{background:#fce7f3;color:#831843;border:1px solid #f9a8d4}
        .compare-full-details{background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;padding:12px 16px;margin-bottom:10px}
        .compare-full-details summary{cursor:pointer;font-weight:700;color:#006400;list-style:none;display:flex;justify-content:space-between;align-items:center;gap:10px;flex-wrap:wrap;padding:4px 0}
        .compare-full-details summary::-webkit-details-marker{display:none}
        .compare-full-details summary::before{content:'▸';display:inline-block;margin-right:8px;color:#006400}
        .compare-full-details[open] summary::before{content:'▾'}
        .details-cat-name{font-size:1em}
        .details-counts{font-size:.85em;font-weight:500;color:#64748b;background:#f1f5f9;padding:3px 10px;border-radius:12px;margin-left:auto}
        .compare-full-profile{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-top:14px;padding-top:14px;border-top:1px dashed #e2e8f0}
        .profile-col{background:#f8fafc;border-radius:10px;padding:14px;border:1px solid #e2e8f0}
        .profile-head{display:flex;justify-content:space-between;align-items:center;gap:8px;margin-bottom:10px;padding-bottom:8px;border-bottom:2px solid #e2e8f0;flex-wrap:wrap}
        .profile-sample-label{font-weight:700;color:#1e293b;font-size:.92em;word-break:break-word}
        .profile-sample-count{font-size:.78em;font-weight:600;color:#006400;background:#dcfce7;padding:2px 10px;border-radius:12px;white-space:nowrap}
        .empty-note{font-size:.85em;color:#94a3b8;font-style:italic;padding:6px 0;display:inline-block}
        .results-counter{font-size:.9em;color:#64748b;font-weight:600;margin-left:auto}
        .no-alerts{padding:20px;background:#f0fdf4;color:#14532d;border-radius:8px;border-left:4px solid #16a34a}
        @media (max-width:900px){.compare-full-profile{grid-template-columns:1fr}.compare-row{grid-template-columns:1fr}.compare-visuals-grid{grid-template-columns:1fr}}
        @media print{body *{visibility:hidden}.tab-content.active,.tab-content.active *{visibility:visible}.tab-content.active{position:absolute;left:0;top:0;width:100%;padding:20px;box-shadow:none;border-radius:0}.print-section-btn,.tab-navigation,.dashboard-grid,.search-box,.action-buttons,.grouping-controls{display:none !important}.data-table{page-break-inside:auto}.data-table tr{page-break-inside:avoid;page-break-after:auto}}
        @media (max-width:768px){.container{padding:10px}.main-header{padding:20px}.main-header h1{font-size:2em}.tab-button{padding:8px 12px;font-size:0.8em}.dashboard-grid{grid-template-columns:repeat(auto-fit,minmax(180px,1fr))}.data-table{font-size:0.8em}body{min-width:auto;overflow-x:auto}}
        </style>
        """

    def _get_js(self, typing_json: str) -> str:
        return f"""
        <script>
        var sampleTyping = {typing_json};
        var originalGenomeLists = {{}};

        function escapeHtml(s) {{ return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }}

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

        function highlightTableCells(tableId, searchId) {{
            var el = document.getElementById(searchId);
            if (!el) return;
            var filter = el.value.toUpperCase().trim();
            var table = document.getElementById(tableId);
            if (!table) return;
            table.querySelectorAll('td').forEach(function(td) {{
                td.classList.remove('cell-highlight');
                if (filter && (td.textContent || td.innerText).toUpperCase().indexOf(filter) > -1) {{
                    td.classList.add('cell-highlight');
                }}
            }});
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
            if (groupBy === "pasteur_st") {{ var st = info.pasteur_st; return st !== 'ND' ? "ST " + st : "ND"; }}
            if (groupBy === "oxford_st") {{ var st = info.oxford_st; return st !== 'ND' ? "ST " + st : "ND"; }}
            if (groupBy === "k_locus") return info.k_locus;
            if (groupBy === "o_locus") return info.o_locus;
            if (groupBy === "capsule") return info.capsule;
            if (groupBy === "st_k") {{ var st = info.pasteur_st; var pre = st !== 'ND' ? "ST " + st : "ND"; return pre + " - " + info.k_locus; }}
            if (groupBy === "st_o") {{ var st = info.pasteur_st; var pre = st !== 'ND' ? "ST " + st : "ND"; return pre + " - " + info.o_locus; }}
            if (groupBy === "k_o") return info.k_locus + ":" + info.o_locus;
            if (groupBy === "st_ko") {{ var st = info.pasteur_st; var pre = st !== 'ND' ? "ST " + st : "ND"; return pre + " - " + info.k_locus + ":" + info.o_locus; }}
            return "Unknown";
        }}

        function groupRowGenomes(row, groupBy, originalList) {{
            var genomesCell = null;
            for (var i = 0; i < row.cells.length; i++) {{
                if (row.cells[i].querySelector('.genome-list')) {{ genomesCell = row.cells[i]; break; }}
            }}
            if (!genomesCell) return;
            var genomes = originalList.slice();
            if (genomes.length === 0) {{ genomesCell.innerHTML = '<div class="genome-list">None</div>'; return; }}
            var groups = {{}};
            genomes.forEach(function(genome) {{
                var key = getTypingValue(genome, groupBy);
                if (!groups[key]) groups[key] = [];
                groups[key].push(genome);
            }});
            var html = '<div class="genome-list">';
            for (var key in groups) {{
                var tags = groups[key].map(g => '<span class="genome-tag">' + g + '</span>').join('');
                html += '<div class="genome-group"><div class="genome-group-header">' + escapeHtml(key) + '</div><div class="genome-group-tags">' + tags + '</div></div>';
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
                        if (row.cells[j].querySelector('.genome-list')) {{ genomesCell = row.cells[j]; break; }}
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
                var activeBtn = container.querySelector('.group-btn[data-group="' + groupBy + '"]');
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
                    if (row.cells[j].querySelector('.genome-list')) {{ genomesCell = row.cells[j]; break; }}
                }}
                if (genomesCell) {{
                    var tags = original.map(g => '<span class="genome-tag">' + g + '</span>').join('');
                    genomesCell.innerHTML = '<div class="genome-list">' + tags + '</div>';
                }}
            }}
            var container = table.closest('.tab-content');
            if (container) {{
                var btns = container.querySelectorAll('.group-btn');
                btns.forEach(btn => btn.classList.remove('active'));
            }}
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
            headers.forEach((th, idx) => {{
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
            win.document.write('<html><head><title>Print</title><style>' + document.querySelector('style').textContent + '</style></head><body>' + content.innerHTML + '</body></html>');
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

        function setCompareMode(mode) {{
            document.querySelectorAll('.mode-btn').forEach(b => b.classList.toggle('active', b.dataset.mode === mode));
            document.querySelectorAll('.compare-mode').forEach(m => m.classList.remove('active'));
            var el = document.getElementById('mode-' + mode);
            if (el) el.classList.add('active');
            if (mode === 'cluster') runCluster();
        }}

        function computeSimilarity(A, B) {{
            var tk = Object.keys(A.typing);
            var matches = 0;
            tk.forEach(function(k) {{ if ((A.typing[k] || '') === (B.typing[k] || '')) matches++; }});
            var tp = tk.length ? (matches / tk.length) * 100 : 0;
            var cats = ['AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations'];
            var cs = [];
            cats.forEach(function(c) {{
                var sa = new Set(A.genes[c] || []), sb = new Set(B.genes[c] || []);
                var u = new Set([...sa, ...sb]);
                if (!u.size) return;
                var sh = 0;
                u.forEach(function(g) {{ if (sa.has(g) && sb.has(g)) sh++; }});
                cs.push((sh / u.size) * 100);
            }});
            var gp = cs.length ? cs.reduce((a, b) => a + b, 0) / cs.length : 0;
            return {{ typing: Math.round(tp), genes: Math.round(gp), overall: Math.round(tp * 0.5 + gp * 0.5), typingMatches: matches, typingTotal: tk.length }};
        }}

        function similarityColor(p) {{
            if (p >= 95) return '#059669';
            if (p >= 90) return '#10b981';
            if (p >= 80) return '#84cc16';
            if (p >= 70) return '#eab308';
            if (p >= 50) return '#f59e0b';
            if (p >= 30) return '#f97316';
            return '#dc2626';
        }}

        function drawGauge(p) {{
            var c = similarityColor(p), r = 60, cc = 2 * Math.PI * r, off = cc * (1 - p / 100);
            return '<svg viewBox="0 0 160 160" class="similarity-gauge">' +
                '<circle cx="80" cy="80" r="' + r + '" fill="none" stroke="#e2e8f0" stroke-width="14"/>' +
                '<circle cx="80" cy="80" r="' + r + '" fill="none" stroke="' + c + '" stroke-width="14" stroke-dasharray="' + cc + '" stroke-dashoffset="' + off + '" stroke-linecap="round" transform="rotate(-90 80 80)"/>' +
                '<text x="80" y="76" text-anchor="middle" font-size="28" font-weight="700" fill="' + c + '">' + p + '%</text>' +
                '<text x="80" y="98" text-anchor="middle" font-size="11" fill="#64748b">overall</text></svg>';
        }}

        function typingValuePill(va, vb) {{
            var m = va === vb, c = m ? '#10b981' : '#dc2626', bg = m ? '#dcfce7' : '#fee2e2', ic = m ? '✓' : '✗';
            return '<span class="typing-pill" style="background:' + bg + ';color:' + c + ';border-color:' + c + ';">' + ic + ' ' + escapeHtml(va) + '</span>';
        }}

        function runCompare() {{
            var data = window.STAPHSCOPE_COMPARE_DATA || {{}};
            var a = document.getElementById('compare-a').value, b = document.getElementById('compare-b').value;
            if (!a || !b || a === b) {{
                document.getElementById('compare-output').innerHTML = '<div class="compare-warning">Pick two different samples.</div>';
                document.getElementById('compare-visuals').innerHTML = '';
                document.getElementById('compare-verdict').style.display = 'none';
                return;
            }}
            var A = data[a], B = data[b];
            if (!A || !B) {{
                document.getElementById('compare-output').innerHTML = '<div class="compare-warning">Sample data not found.</div>';
                return;
            }}
            var sim = computeSimilarity(A, B);
            document.getElementById('compare-visuals').innerHTML =
                '<div class="compare-visuals-grid">' +
                '<div class="visual-card gauge-card">' + drawGauge(sim.overall) + '</div>' +
                '<div class="visual-card metrics-card">' +
                '<div class="metric-block"><div class="metric-label">Typing match</div>' +
                '<div class="metric-bar-track"><div class="metric-bar-fill" style="width:' + sim.typing + '%;background:' + similarityColor(sim.typing) + ';"></div></div>' +
                '<div class="metric-value">' + sim.typingMatches + '/' + sim.typingTotal + ' fields (' + sim.typing + '%)</div></div>' +
                '<div class="metric-block"><div class="metric-label">Gene content match</div>' +
                '<div class="metric-bar-track"><div class="metric-bar-fill" style="width:' + sim.genes + '%;background:' + similarityColor(sim.genes) + ';"></div></div>' +
                '<div class="metric-value">' + sim.genes + '% Jaccard similarity</div></div>' +
                '</div></div>';
            var tk = Object.keys(A.typing);
            var typingRows = '';
            tk.forEach(function(k) {{
                var va = A.typing[k] || 'Not Assigned', vb = B.typing[k] || 'Not Assigned';
                var diff = va !== vb;
                typingRows += '<div class="compare-row ' + (diff ? 'diff' : 'match') + '">' +
                    '<span class="row-label">' + k + '</span>' +
                    '<span class="row-val">' + typingValuePill(va, vb) + '</span>' +
                    '<span class="row-val">' + typingValuePill(vb, va) + '</span></div>';
            }});
            var cats = ['AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations'];
            var geneHtml = '';
            cats.forEach(function(cat) {{
                var ga = new Set(A.genes[cat] || []), gb = new Set(B.genes[cat] || []);
                var u = new Set([...ga, ...gb]);
                if (!u.size) return;
                var sh = [], oa = [], ob = [];
                u.forEach(function(g) {{ if (ga.has(g) && gb.has(g)) sh.push(g); else if (ga.has(g)) oa.push(g); else ob.push(g); }});
                var total = u.size;
                var sP = sh.length / total * 100, aP = oa.length / total * 100, bP = ob.length / total * 100;
                geneHtml += '<div class="compare-gene-block"><h4>' + cat + ' <span class="gene-count">' +
                    sh.length + ' shared · ' + oa.length + ' only A · ' + ob.length + ' only B (' + sP.toFixed(0) + '% match)</span></h4>' +
                    '<div class="gene-proportion-bar">' +
                    '<div class="gene-prop-segment shared" style="width:' + sP + '%;"></div>' +
                    '<div class="gene-prop-segment only-a" style="width:' + aP + '%;"></div>' +
                    '<div class="gene-prop-segment only-b" style="width:' + bP + '%;"></div>' +
                    '</div>' +
                    '<div class="gene-columns">' +
                    '<div class="gene-col"><div class="gene-col-head">Shared (' + sh.length + ')</div>' +
                    (sh.map(g => '<span class="compare-gene shared">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') + '</div>' +
                    '<div class="gene-col"><div class="gene-col-head">Only in A (' + oa.length + ')</div>' +
                    (oa.map(g => '<span class="compare-gene only-a">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') + '</div>' +
                    '<div class="gene-col"><div class="gene-col-head">Only in B (' + ob.length + ')</div>' +
                    (ob.map(g => '<span class="compare-gene only-b">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') + '</div>' +
                    '</div></div>';
            }});
            var fullProfile = '<div class="compare-section"><h3>Full Gene Profiles</h3><p style="font-size:.9em;color:#64748b;margin-bottom:12px;">Click any category to expand both samples\\' complete gene lists side by side.</p>';
            cats.forEach(function(cat) {{
                var ga = A.genes[cat] || [], gb = B.genes[cat] || [];
                if (!ga.length && !gb.length) return;
                var tagsA = ga.length ? ga.map(g => '<span class="genome-tag">' + escapeHtml(g) + '</span>').join('') : '<em class="empty-note">None detected</em>';
                var tagsB = gb.length ? gb.map(g => '<span class="genome-tag">' + escapeHtml(g) + '</span>').join('') : '<em class="empty-note">None detected</em>';
                fullProfile += '<details class="vir-details compare-full-details"><summary><span class="details-cat-name">' + cat + '</span>' +
                    '<span class="details-counts">A: ' + ga.length + ' · B: ' + gb.length + '</span></summary>' +
                    '<div class="compare-full-profile">' +
                    '<div class="profile-col"><div class="profile-head"><span class="profile-sample-label">' + escapeHtml(a) + '</span>' +
                    '<span class="profile-sample-count">' + ga.length + ' gene(s)</span></div>' +
                    '<div class="genome-list">' + tagsA + '</div></div>' +
                    '<div class="profile-col"><div class="profile-head"><span class="profile-sample-label">' + escapeHtml(b) + '</span>' +
                    '<span class="profile-sample-count">' + gb.length + ' gene(s)</span></div>' +
                    '<div class="genome-list">' + tagsB + '</div></div></div></details>';
            }});
            fullProfile += '</div>';
            var v = document.getElementById('compare-verdict');
            if (sim.overall >= 95) v.innerHTML = '<div class="compare-verdict-box danger"><i class="fas fa-exclamation-triangle"></i><strong>Near-identical profile (' + sim.overall + '%)</strong> — likely the same strain. Investigate as a possible <em>transmission pair</em> or duplicate sample.</div>';
            else if (sim.overall >= 85) v.innerHTML = '<div class="compare-verdict-box warning"><i class="fas fa-info-circle"></i><strong>High similarity (' + sim.overall + '%)</strong> — very likely the same clonal lineage.</div>';
            else if (sim.overall >= 60) v.innerHTML = '<div class="compare-verdict-box warning"><i class="fas fa-info-circle"></i><strong>Moderate similarity (' + sim.overall + '%)</strong> — related lineage, but not a recent transmission pair.</div>';
            else v.innerHTML = '<div class="compare-verdict-box info"><i class="fas fa-info-circle"></i><strong>Low similarity (' + sim.overall + '%)</strong> — distinct lineages.</div>';
            v.style.display = 'block';
            var tableHtml = '<div class="compare-section"><h3>Typing Profile</h3><div class="compare-table">' +
                '<div class="compare-row header"><span class="row-label">Field</span>' +
                '<span class="row-val"><strong>' + escapeHtml(a) + '</strong></span>' +
                '<span class="row-val"><strong>' + escapeHtml(b) + '</strong></span></div>' +
                typingRows + '</div></div>';
            document.getElementById('compare-output').innerHTML =
                tableHtml + '<div class="compare-section"><h3>Gene Content</h3>' +
                (geneHtml || '<p><em>No gene data available.</em></p>') + '</div>' + fullProfile;
        }}

        function swapCompare() {{
            var a = document.getElementById('compare-a'), b = document.getElementById('compare-b');
            var t = a.value; a.value = b.value; b.value = t; runCompare();
        }}

        function resetCompare() {{
            var sel = document.getElementById('compare-a');
            if (!sel) return;
            var opts = sel.options;
            if (opts.length >= 2) {{
                document.getElementById('compare-a').value = opts[0].value;
                document.getElementById('compare-b').value = opts[1].value;
            }}
            document.getElementById('compare-output').innerHTML = '';
            document.getElementById('compare-visuals').innerHTML = '';
            document.getElementById('compare-verdict').style.display = 'none';
        }}

        function getPrecomputedSim(i, j) {{
            if (i === j) return 100;
            if (i > j) {{ var t = i; i = j; j = t; }}
            var sim = window.STAPHSCOPE_SIM;
            if (!sim || !sim.flat) return 0;
            var n = sim.n;
            return sim.flat[i * n - i * (i + 1) / 2 + (j - i - 1)];
        }}

        function runCluster() {{
            var simData = window.STAPHSCOPE_SIM;
            if (!simData || !simData.samples || simData.samples.length < 2) {{
                var el = document.getElementById('cluster-heatmap');
                if (el) el.innerHTML = '<div class="compare-warning">Need at least two samples.</div>';
                return;
            }}
            var samples = simData.samples, n = simData.n;
            var thEl = document.getElementById('cluster-threshold');
            var threshold = thEl ? parseInt(thEl.value) : 90;
            var parent = {{}};
            for (var i = 0; i < n; i++) parent[i] = i;
            function find(x) {{ while (parent[x] !== x) {{ parent[x] = parent[parent[x]]; x = parent[x]; }} return x; }}
            function union(a, b) {{ var ra = find(a), rb = find(b); if (ra !== rb) parent[ra] = rb; }}
            for (var i = 0; i < n; i++) for (var j = i + 1; j < n; j++) if (getPrecomputedSim(i, j) >= threshold) union(i, j);
            var groups = {{}};
            for (var k = 0; k < n; k++) {{ var r = find(k); if (!groups[r]) groups[r] = []; groups[r].push(samples[k]); }}
            var clusters = Object.values(groups).filter(g => g.length >= 2).sort((a, b) => b.length - a.length);
            var cHtml = '';
            if (!clusters.length) {{
                cHtml = '<div class="no-alerts">No clusters detected at ≥' + threshold + '% similarity.</div>';
            }} else {{
                clusters.forEach(function(cl, idx) {{
                    var sum = 0, cnt = 0;
                    for (var a = 0; a < cl.length; a++) for (var b = a + 1; b < cl.length; b++) {{
                        var ia = samples.indexOf(cl[a]), ib = samples.indexOf(cl[b]);
                        sum += getPrecomputedSim(ia, ib); cnt++;
                    }}
                    var mean = cnt ? (sum / cnt).toFixed(1) : '100';
                    var color = similarityColor(parseFloat(mean));
                    cHtml += '<div class="cluster-card" style="border-left-color:' + color + ';">' +
                        '<div class="cluster-header"><span class="cluster-badge" style="background:' + color + ';">Cluster ' + (idx + 1) + '</span>' +
                        '<span class="cluster-stats">' + cl.length + ' samples · mean ' + mean + '% similar</span></div>' +
                        '<div class="cluster-samples">' + cl.map(s => '<span class="cluster-sample-pill">' + escapeHtml(s) + '</span>').join('') + '</div></div>';
                }});
            }}
            document.getElementById('cluster-clusters').innerHTML = cHtml;
            if (n <= 60) renderHeatmapDOM(samples, n); else renderHeatmapCanvas(samples, n);
            var sum2 = document.getElementById('cluster-summary');
            if (sum2) sum2.textContent = clusters.length + ' cluster(s) at ≥' + threshold + '% (' + n + ' samples)';
        }}

        function renderHeatmapDOM(samples, n) {{
            var rowLabelWidth = 220, cellSize = 28;
            var truncate = function(s, max) {{ return s.length > max ? s.substring(0, max - 1) + '…' : s; }};
            var html = '<div class="heatmap-wrapper"><h3>Pairwise Similarity Matrix</h3>' +
                '<p style="font-size:.9em;color:#64748b;margin-bottom:12px;">Hover any cell or label for the full sample name.</p>' +
                '<div style="overflow-x:auto;padding-bottom:10px;">' +
                '<div class="heatmap-grid" style="grid-template-columns: ' + rowLabelWidth + 'px repeat(' + n + ', ' + cellSize + 'px);">';
            html += '<div class="heatmap-cell heatmap-corner"></div>';
            samples.forEach(function(s) {{
                html += '<div class="heatmap-cell heatmap-col-header" title="' + escapeHtml(s) + '"><span class="heatmap-header-label">' + escapeHtml(s) + '</span></div>';
            }});
            samples.forEach(function(s1, i) {{
                html += '<div class="heatmap-cell heatmap-row-header" style="justify-content:flex-start;padding-left:8px;width:' + rowLabelWidth + 'px;" title="' + escapeHtml(s1) + '"><span class="heatmap-header-label">' + escapeHtml(truncate(s1, 26)) + '</span></div>';
                samples.forEach(function(s2, j) {{
                    var pct = (i === j) ? 100 : getPrecomputedSim(i, j);
                    var color = (i === j) ? '#1e293b' : similarityColor(pct);
                    var text = (i === j) ? '—' : Math.round(pct);
                    html += '<div class="heatmap-cell" style="background:' + color + ';color:white;" title="' + escapeHtml(s1) + ' vs ' + escapeHtml(s2) + ': ' + pct + '%">' + text + '</div>';
                }});
            }});
            html += '</div></div></div>';
            document.getElementById('cluster-heatmap').innerHTML = html;
        }}

        function renderHeatmapCanvas(samples, n) {{
            var cellSize = Math.max(3, Math.min(12, Math.floor(700 / n)));
            var rowLabelW = 160, headerH = 60;
            var canvasW = rowLabelW + n * cellSize + 20, canvasH = headerH + n * cellSize + 20;
            var canvas = document.createElement('canvas');
            canvas.width = canvasW; canvas.height = canvasH;
            canvas.style.display = 'block'; canvas.style.maxWidth = '100%'; canvas.style.cursor = 'crosshair';
            var ctx = canvas.getContext('2d');
            ctx.font = Math.max(8, Math.min(11, cellSize - 1)) + 'px sans-serif';
            ctx.textBaseline = 'middle';
            for (var i = 0; i < n; i++) for (var j = 0; j < n; j++) {{
                var pct = (i === j) ? 100 : getPrecomputedSim(i, j);
                ctx.fillStyle = (i === j) ? '#1e293b' : similarityColor(pct);
                ctx.fillRect(rowLabelW + j * cellSize, headerH + i * cellSize, cellSize, cellSize);
            }}
            ctx.fillStyle = '#333'; ctx.textAlign = 'right';
            var rowStep = Math.max(1, Math.floor(n / 40));
            for (var i = 0; i < n; i += rowStep) {{
                var lb = samples[i].length > 20 ? samples[i].substring(0, 19) + '…' : samples[i];
                ctx.fillText(lb, rowLabelW - 5, headerH + i * cellSize + cellSize / 2);
            }}
            ctx.save(); ctx.translate(0, headerH - 5); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'left';
            var colStep = Math.max(1, Math.floor(n / 40));
            for (var j = 0; j < n; j += colStep) {{
                var lb2 = samples[j].length > 15 ? samples[j].substring(0, 14) + '…' : samples[j];
                ctx.fillText(lb2, 5, rowLabelW + j * cellSize + cellSize / 2);
            }}
            ctx.restore();
            canvas.addEventListener('click', function(e) {{
                var rect = canvas.getBoundingClientRect();
                var x = (e.clientX - rect.left) * (canvas.width / rect.width) - rowLabelW;
                var y = (e.clientY - rect.top) * (canvas.height / rect.height) - headerH;
                var jj = Math.floor(x / cellSize), ii = Math.floor(y / cellSize);
                if (ii >= 0 && ii < n && jj >= 0 && jj < n && ii !== jj) {{
                    document.getElementById('compare-a').value = samples[ii];
                    document.getElementById('compare-b').value = samples[jj];
                    setCompareMode('pair');
                    runCompare();
                }}
            }});
            var wrapper = document.getElementById('cluster-heatmap');
            wrapper.innerHTML = '<div class="heatmap-wrapper"><h3>Pairwise Similarity Matrix</h3>' +
                '<p style="font-size:.9em;color:#64748b;margin-bottom:12px;">Canvas-rendered for ' + n + ' samples. Click any cell to open pair compare.</p>' +
                '<div id="canvas-heatmap-host" style="overflow-x:auto;"></div></div>';
            wrapper.querySelector('#canvas-heatmap-host').appendChild(canvas);
        }}

        document.addEventListener('DOMContentLoaded', function() {{
            var hash = window.location.hash.substring(1);
            if (hash) {{
                var btn = document.querySelector('.tab-button.' + hash);
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

    # --------------------------------------------------------------------------
    # Helper: genome tags
    # --------------------------------------------------------------------------
    def _make_genome_tags(self, genomes: List[str]) -> str:
        return ''.join(f'<span class="genome-tag">{gen}</span>' for gen in sorted(genomes))

    # --------------------------------------------------------------------------
    # SECTION: SUMMARY 
    # --------------------------------------------------------------------------
    def _summary_section(self, kwargs, carb_count, env_count):
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
        mutation_count = len(kwargs.get('mutation_data', {}).get('mutations', []))
        return f"""
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>MDR/XDR A. baumannii Analysis with Environmental Co-Selection</h3><p>This comprehensive gene-centric report analyzes <strong>{total}</strong> A. baumannii genomes with focus on carbapenem resistance, colistin resistance, biofilm formation, and environmental co-selection markers. Each gene is shown with ALL genomes that contain it for easy tracking of resistance spread.</p><p><strong>NEW:</strong> Dynamic grouping by typing (ST, K, OCL, capsule) and gene co‑occurrence analysis.</p></div></div>
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
        <tr><td>Point Mutations</td><td><strong>{mutation_count}</strong></td><td>AMRFinderPlus mutations (gyrA, parC, 23S rRNA, etc.)</td></tr>
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

    # --------------------------------------------------------------------------
    # SECTION: SAMPLES — v3.2.0 changes: sticky header, virulence expander, mobileOG column
    # --------------------------------------------------------------------------
    def _samples_section(self, kwargs):
        samples = kwargs['samples_data']
        gene_centric = kwargs['gene_centric']
        mge_per_sample = kwargs.get('mge_data', {}).get('per_sample', {})
        sample_vir = defaultdict(set)
        for db_name, genes in gene_centric.get('virulence_databases', {}).items():
            for g in genes:
                for genome in g['genomes']:
                    sample_vir[genome].add(g['gene'])

        html = """
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><h3>Population Structure Overview</h3><p>This table summarises the key typing results for each genome. Understanding the population structure helps identify dominant clones, track outbreaks, and link genotypes to phenotypes.</p><ul><li><strong>MLST (Sequence Type)</strong>: Gold standard for global epidemiology.</li><li><strong>Capsule Typing (K:OCL)</strong>: K (capsule) and OCL (lipooligosaccharide) loci are critical for virulence and immune evasion.</li><li><strong>Virulence</strong>: Expand the column to view the individual virulence genes detected in each genome.</li><li><strong>MobileOG Hits</strong>: Total number of mobileOG-db protein family matches identified in the genome (see MGE tab).</li><li><strong>ND</strong>: Not Determined.</li></ul></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('samples-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById('search-samples').value=''; searchTable('samples-table','search-samples'); document.getElementById('highlight-samples').value=''; highlightTableCells('samples-table','highlight-samples')"><i class="fas fa-sync"></i> Clear Search</button><button class="action-btn btn-success" onclick="exportTableToCSV('samples-table','acinetobacter_samples.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        <div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-samples" onkeyup="searchTable('samples-table','search-samples')" placeholder="🔍 Filter rows..."><input type="text" class="search-box" id="highlight-samples" onkeyup="highlightTableCells('samples-table','highlight-samples')" placeholder="✨ Highlight text..."></div>
        <div class="master-scrollable-container sticky-y">
            <table id="samples-table" class="data-table">
                <thead>
                    <tr><th data-sort="string">Sample</th><th data-sort="string">Pasteur ST</th><th data-sort="string">Oxford ST</th><th data-sort="string">K Locus</th><th data-sort="string">OCL Locus</th><th data-sort="number">Virulence</th><th data-sort="number">MobileOG Hits</th></tr>
                </thead>
                <tbody>
        """
        for sample, data in samples.items():
            pasteur = data.get('pasteur_mlst',{}).get('ST','ND')
            oxford = data.get('oxford_mlst',{}).get('ST','ND')
            k = data.get('kaptive',{}).get('K_Locus','ND')
            o = data.get('kaptive',{}).get('O_Locus','ND')
            vir_genes = sorted(sample_vir.get(sample, []))
            vir_count = len(vir_genes)
            mge_rec = mge_per_sample.get(sample, {})
            mge_hits = 0
            if isinstance(mge_rec, dict):
                for cand in ('mobileOG Hits', 'Key MGE-assoc.', 'mobileOG_hits', 'Total'):
                    if cand in mge_rec:
                        mge_hits = mge_rec[cand]
                        break
            if not isinstance(mge_hits, (int, float)):
                mge_hits = 0
            vir_tags = ''.join(f'<span class="genome-tag">{g}</span>' for g in vir_genes)
            vir_details = (
                f'<details class="vir-details"><summary>{vir_count} gene(s)</summary>'
                f'<div class="vir-gene-tags">{vir_tags or "<em>None</em>"}</div></details>'
            )
            pasteur_disp = f"ST{pasteur}" if pasteur != 'ND' else 'ND'
            oxford_disp = f"ST{oxford}" if oxford != 'ND' else 'ND'
            html += f'<tr><td><strong>{sample}</strong></td><td>{pasteur_disp}</td><td>{oxford_disp}</td><td>{k}</td><td>{o}</td><td>{vir_details}</td><td>{mge_hits}</td></tr>'
        html += """
                </tbody>
            </table>
        </div>
        """
        return html

    # --------------------------------------------------------------------------
    # SECTION: QC — v3.2.0 changes: sticky header
    # --------------------------------------------------------------------------
    def _qc_section(self, kwargs):
        qc = kwargs.get('qc_data', {})
        if not qc:
            return '<div class="alert-box alert-warning"><i class="fas fa-exclamation-circle"></i><div>FASTA QC file not found or could not parse.</div></div>'
        has_ani = any(
            'species_check' in vals and isinstance(vals['species_check'], dict) and 'best_match' in vals['species_check']
            for vals in qc.values()
        )
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
        html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'qc-tab\')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById(\'search-qc\').value=\'\'; searchTable(\'qc-table\',\'search-qc\'); document.getElementById(\'highlight-qc\').value=\'\'; highlightTableCells(\'qc-table\',\'highlight-qc\')"><i class="fas fa-sync"></i> Clear Search</button><button class="action-btn btn-success" onclick="exportTableToCSV(\'qc-table\',\'fasta_qc.csv\')"><i class="fas fa-download"></i> Export CSV</button></div>'
        html += '<div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-qc" onkeyup="searchTable(\'qc-table\',\'search-qc\')" placeholder="🔍 Filter rows..."><input type="text" class="search-box" id="highlight-qc" onkeyup="highlightTableCells(\'qc-table\',\'highlight-qc\')" placeholder="✨ Highlight text..."></div>'
        html += '<div class="master-scrollable-container sticky-y"><table id="qc-table" class="data-table"><thead><tr><th data-sort="string" style="min-width:200px">Sample</th>'
        for m in metrics:
            html += f'<th data-sort="number">{m}</th>'
        if has_ani:
            html += '<th data-sort="string">Best Species</th><th data-sort="number">ANI (%)</th><th data-sort="string">Confirmed</th>'
        html += '</tr></thead><tbody>'
        for sample, vals in sorted(qc.items()):
            html += f'<tr><td style="white-space:nowrap"><strong>{sample}</strong></td>'
            for m in metrics:
                v = vals.get(m, 'ND')
                if isinstance(v, float):
                    v = f"{v:,.0f}" if v > 1000 else f"{v:.2f}"
                html += f'<td>{v}</td>'
            if has_ani:
                sc = vals.get('species_check', {})
                best_species = sc.get('best_match', 'ND')
                ani = sc.get('ani_percent', 'ND')
                confirmed = '✅ Yes' if sc.get('passed') else '❌ No'
                html += f'<td>{best_species}</td><td>{ani}</td><td>{confirmed}</td>'
            html += '</tr>'
        html += '</tbody></table></div>'
        return html

    # --------------------------------------------------------------------------
    # SECTION: MLST 
    # --------------------------------------------------------------------------
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
        html += """
        <div class="card" style="margin-top:30px;">
            <h2 style="color: #333; border-bottom: 2px solid #FF9800; padding-bottom: 10px;">📚 International Clone Literature Sources</h2>
            <div class="alert-box alert-info" style="border-left-color:#FF9800; background:#fff3cd;">
                <i class="fas fa-info-circle fa-2x" style="color:#FF9800;"></i>
                <div>
                    <strong>⚠️ IC Assignment Caveat:</strong>
                    <p style="margin-top:8px; font-size:0.95em;">
                        International Clone (IC) assignments in AcinetoScope are based <strong>solely on MLST Sequence Type (ST)</strong> 
                        matching against a literature-curated lookup table. <strong>OXA-51 variants are not used</strong> for IC confirmation 
                        in this pipeline. For definitive IC assignment, users are encouraged to confirm with OXA-51-like 
                        gene analysis or whole‑genome phylogenetics.
                    </p>
                </div>
            </div>
            <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(280px, 1fr)); gap:15px; margin:20px 0;">
                <div style="background:#f8f9fa; padding:15px; border-radius:8px; border-left:4px solid #00695c;">
                    <h4 style="color:#00695c; margin:0 0 5px 0;">IC1</h4>
                    <p style="font-size:0.85em; color:#333; margin:0;"><strong>Sources:</strong><br>
                        <a href="https://doi.org/10.1016/j.ijantimicag.2012.09.008" target="_blank" style="color:#00695c;">Zarrilli R, et al. Int J Antimicrob Agents. 2013;41:11-19.</a><br>
                        <a href="https://doi.org/10.3390/microorganisms11082115" target="_blank" style="color:#00695c;">Shelenkov A, et al. Microorganisms. 2023;11(8):2115.</a>
                    </p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:8px; border-left:4px solid #00695c;">
                    <h4 style="color:#00695c; margin:0 0 5px 0;">IC2</h4>
                    <p style="font-size:0.85em; color:#333; margin:0;"><strong>Sources:</strong><br>
                        <a href="https://doi.org/10.1016/j.ijantimicag.2012.09.008" target="_blank" style="color:#00695c;">Zarrilli R, et al. Int J Antimicrob Agents. 2013;41:11-19.</a><br>
                        <a href="https://doi.org/10.3390/microorganisms11082115" target="_blank" style="color:#00695c;">Shelenkov A, et al. Microorganisms. 2023;11(8):2115.</a>
                    </p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:8px; border-left:4px solid #00695c;">
                    <h4 style="color:#00695c; margin:0 0 5px 0;">IC3 - IC9</h4>
                    <p style="font-size:0.85em; color:#333; margin:0;"><strong>Sources:</strong><br>
                        <a href="https://doi.org/10.3389/fmicb.2019.00930" target="_blank" style="color:#00695c;">Gaiarsa S, et al. Front Microbiol. 2019;10:930.</a><br>
                        <a href="https://doi.org/10.3390/microorganisms11082115" target="_blank" style="color:#00695c;">Shelenkov A, et al. Microorganisms. 2023;11(8):2115.</a>
                    </p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:8px; border-left:4px solid #FF9800;">
                    <h4 style="color:#FF9800; margin:0 0 5px 0;">IC10</h4>
                    <p style="font-size:0.85em; color:#333; margin:0;"><strong>Sources:</strong><br>
                        <a href="https://doi.org/10.3390/microorganisms11082115" target="_blank" style="color:#FF9800;">Shelenkov A, et al. Microorganisms. 2023;11(8):2115.</a>
                    </p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:8px; border-left:4px solid #9C27B0;">
                    <h4 style="color:#9C27B0; margin:0 0 5px 0;">IC11</h4>
                    <p style="font-size:0.85em; color:#333; margin:0;"><strong>Sources:</strong><br>
                        <a href="https://doi.org/10.1016/j.ijantimicag.2023.106866" target="_blank" style="color:#9C27B0;">Hansen F, et al. Int J Antimicrob Agents. 2023;62(2):106866.</a><br>
                        <a href="https://doi.org/10.1128/msphere.00276-24" target="_blank" style="color:#9C27B0;">Xu A, et al. mSphere. 2024;9(6):e0027624.</a>
                    </p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:8px; border-left:4px solid #9C27B0;">
                    <h4 style="color:#9C27B0; margin:0 0 5px 0;">IC12</h4>
                    <p style="font-size:0.85em; color:#333; margin:0;"><strong>Sources:</strong><br>
                        <a href="https://doi.org/10.1099/mgen.0.001572" target="_blank" style="color:#9C27B0;">Karah N, et al. Microb Genom. 2025;11(11):001572.</a>
                    </p>
                </div>
            </div>
        </div>
        """
        return html

    # --------------------------------------------------------------------------
    # SECTION: KAPTIVE 
    # --------------------------------------------------------------------------
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
                        <strong>Kaptive3</strong> developed by 
                        <a href="https://github.com/klebgenomics/Kaptive" target="_blank" style="color: #9C27B0; font-weight: bold;">Stanton et al.</a> 
                        (<a href="https://doi.org/10.1099/mgen.0.001428" target="_blank" style="color: #9C27B0; font-weight: 600;">Mgen. 2025;11(6):001428</a>)
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
            count = cnt if isinstance(cnt, int) else cnt.get('count',0)
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
            count = cnt if isinstance(cnt, int) else cnt.get('count',0)
            pct = (count/totalo)*100 if totalo else 0
            freq = f"{count} ({pct:.1f}%)"
            sample_list = [s for s,d in samples.items() if d.get('kaptive',{}).get('O_Locus')==o]
            tags = self._make_genome_tags(sample_list)
            html += f'<tr><td><strong>{o}</strong></td><td>{freq}</td><td><div class="genome-list">{tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html

    # --------------------------------------------------------------------------
    # SECTION: COMBINATIONS 
    # --------------------------------------------------------------------------
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

    # --------------------------------------------------------------------------
    # SECTION: AMR — v3.2.0 changes: educational content in expanders
    # --------------------------------------------------------------------------
    def _amr_section(self, kwargs):
        gene_centric = kwargs['gene_centric']
        all_genes = []
        for db in gene_centric.get('amr_databases', {}).values():
            all_genes.extend(db)
        all_genes.sort(key=lambda x: x['count'], reverse=True)

        html = """
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

        educational_box = """
        <div style="background:#e8f5e9; padding:16px; border-radius:8px;">
            <h4 style="margin:0 0 10px 0; color:#00695c; font-size:1.05em;">🔬 Understanding Database-Specific Gene Reporting</h4>
            <p style="margin:6px 0; font-size:0.95em; line-height:1.5;"><strong>Why do I see the same gene reported by multiple databases?</strong></p>
            <p style="margin:6px 0 10px 0; font-size:0.93em; line-height:1.6; color:#333;">AcinetoScope uses <strong>10 different databases</strong> because <strong>no single database is comprehensive</strong>. Each has unique strengths and limitations.</p>
            <ul style="margin:6px 0 10px 20px; font-size:0.93em; line-height:1.6; color:#333;">
                <li                    <li><strong>CARD</strong> – Strict SNP-based cutoffs; excellent for precise allele identification.</li>
                    <li><strong>ResFinder</strong> – Highly sensitive; appends <code>_1</code> to primary allele variants (e.g., <code>blaOXA-23_1</code> → same as <code>blaOXA-23</code>).</li>
                    <li><strong>MEGARes</strong> – Optimised for metagenomics; includes biocide and metal resistance.</li>
                    <li><strong>AMRFinderPlus</strong> – NCBI gold standard; includes point mutations.</li>
                    <li><strong>VFDB</strong> – Comprehensive virulence factor database.</li>
                    <li><strong>BacMet2</strong> – Biocide and heavy metal resistance genes.</li>
                </ul>
            <p style="margin:6px 0 8px 0; font-size:0.95em; line-height:1.5;"><strong>💡 How to interpret duplicate entries:</strong></p>
            <ul style="margin:6px 0 10px 20px; font-size:0.93em; line-height:1.6; color:#333;">
                <li><span style="color:#00695c; font-weight:bold;">High confidence</span> – Gene found in <strong>3+ databases</strong>.</li>
                <li><span style="color:#ff9800; font-weight:bold;">Moderate confidence</span> – Gene found in <strong>2 databases</strong>.</li>
                <li><span style="color:#dc3545; font-weight:bold;">Low confidence</span> – Gene found in <strong>only 1 database</strong>; may be a database-specific artifact.</li>
                <li>⚠️ <strong>ResFinder naming:</strong> <code>_1</code> suffix (e.g., <code>blaOXA-23_1</code>) indicates the primary allele variant — it is the same gene as <code>blaOXA-23</code>.</li>
                <li>⚠️ <strong>Intrinsic vs Acquired:</strong> Genes like <em>bla</em>OXA-66 and <em>bla</em>OXA-69 are intrinsic OXA-51-like genes in <em>A. baumannii</em> and do not indicate acquired carbapenem resistance. Users are encouraged to consult the literature for gene-specific context.</li>
            </ul>
            <p style="margin:6px 0 0 0; font-size:0.92em; line-height:1.5; background:#fff3cd; padding:8px 14px; border-radius:4px; border-left:3px solid #ffc107;">
                <i class="fas fa-lightbulb" style="color:#856404;"></i>
                <strong>Pro tip:</strong> The <strong>Database column</strong> shows exactly which database found each gene. Use this to assess confidence and identify database-specific artifacts. <strong>No results are merged or filtered</strong> — we preserve all findings for full transparency.
            </p>
        </div>
        """

        edu_db_roles = """
        <div>
            <p style="color:#666; margin-bottom:15px;">AcinetoScope uses multiple databases to maximise sensitivity and specificity. Each database is <strong>kept separate</strong> in the report to preserve provenance and allow users to assess confidence.</p>
            <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(280px, 1fr)); gap:20px; margin:20px 0;">
                <div style="background:#f8f9fa; padding:15px; border-radius:10px; border-left:4px solid #28a745;">
                    <h4 style="color:#28a745; margin:0 0 8px 0;">🟢 CARD</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Strict SNP-based cutoffs, curated protein families.<br><strong>Limitations:</strong> May miss novel variants.<br><strong>Reporting:</strong> Precise allele assignments.</p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:10px; border-left:4px solid #17a2b8;">
                    <h4 style="color:#17a2b8; margin:0 0 8px 0;">🔵 ResFinder</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Highly sensitive for known resistance genes.<br><strong>Limitations:</strong> May misclassify alleles (e.g., <em>bla</em>ADC-25 false positives).<br><strong>Reporting:</strong> Appends <code>_1</code> to primary allele variants.</p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:10px; border-left:4px solid #fd7e14;">
                    <h4 style="color:#fd7e14; margin:0 0 8px 0;">🟠 MEGARes</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Includes biocide and metal resistance genes.<br><strong>Limitations:</strong> Less focused on clinical AMR.<br><strong>Reporting:</strong> Functional metagenomic annotations.</p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:10px; border-left:4px solid #dc3545;">
                    <h4 style="color:#dc3545; margin:0 0 8px 0;">🔴 AMRFinderPlus</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> NCBI gold standard; includes point mutations.<br><strong>Limitations:</strong> Slower to update with novel variants.<br><strong>Reporting:</strong> Comprehensive AMR gene and mutation detection.</p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:10px; border-left:4px solid #E91E63;">
                    <h4 style="color:#E91E63; margin:0 0 8px 0;">🩷 VFDB</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Comprehensive virulence factor database.<br><strong>Limitations:</strong> Focused on virulence, not AMR.<br><strong>Reporting:</strong> Virulence gene detection.</p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:10px; border-left:4px solid #795548;">
                    <h4 style="color:#795548; margin:0 0 8px 0;">🤎 BacMet2</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Biocide and heavy metal resistance genes.<br><strong>Limitations:</strong> Not focused on clinical AMR.<br><strong>Reporting:</strong> Environmental co-selection markers.</p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:10px; border-left:4px solid #6f42c1;">
                    <h4 style="color:#6f42c1; margin:0 0 8px 0;">🟣 ARG-ANNOT</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Broad AMR gene coverage.<br><strong>Limitations:</strong> Older database, less frequently updated.<br><strong>Reporting:</strong> AMR gene detection.</p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:10px; border-left:4px solid #2196F3;">
                    <h4 style="color:#2196F3; margin:0 0 8px 0;">🔵 PlasmidFinder</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Generalist plasmid replicon detection.<br><strong>Limitations:</strong> Not optimised for Acinetobacter.<br><strong>Reporting:</strong> Plasmid replicon markers.</p>
                </div>
                <div style="background:#f8f9fa; padding:15px; border-radius:10px; border-left:4px solid #6c757d;">
                    <h4 style="color:#6c757d; margin:0 0 8px 0;">⚪ NCBI</h4>
                    <p style="font-size:0.9em; color:#333; margin:0;"><strong>Strengths:</strong> Broad general annotation.<br><strong>Limitations:</strong> Not specialised for AMR or virulence.<br><strong>Reporting:</strong> General gene annotation.</p>
                </div>
            </div>
            <div style="background:#e8f5e9; padding:15px; border-radius:8px; border-left:4px solid #00695c; margin-top:15px;">
                <strong>Why we keep results separate:</strong>
                <p style="margin-top:8px; font-size:0.95em;">By preserving <strong>which database</strong> found each gene, users can:</p>
                <ul style="margin:8px 0 0 20px; font-size:0.93em;">
                    <li>Assess confidence — a gene found in 3+ databases is stronger evidence than one found in only 1.</li>
                    <li>Identify database-specific artifacts (e.g., ResFinder <em>bla</em>ADC-25 false positives).</li>
                    <li>Understand differences in gene prevalence between databases.</li>
                    <li>Make informed decisions about which results to trust.</li>
                </ul>
                <p style="margin-top:10px; font-size:0.92em; background:#fff3cd; padding:8px 14px; border-radius:4px; border-left:3px solid #ffc107;">
                    <i class="fas fa-lightbulb"></i> <strong>No results are merged or filtered</strong> — we present all findings transparently.
                </p>
            </div>
        </div>
        """

        caveat = """
        <div class="alert-box alert-warning" style="border-left-color:#ffc107;">
            <i class="fas fa-exclamation-triangle fa-2x"></i>
            <div>
                <strong>⚠️ Important Caveat:</strong> Gene presence does <strong>not</strong> necessarily imply phenotypic resistance. 
                Confirmation by antimicrobial susceptibility testing (AST) is required for clinical decisions. 
                The term "CRAB" (carbapenem-resistant <em>A. baumannii</em>) used in this report is based solely on the detection of carbapenemase genes; 
                clinical correlation with MIC values is essential.
            </div>
        </div>
        """

        group_controls = self._grouping_controls('amr-table')
        search_boxes = self._search_boxes('amr', 'amr-table')
        filter_buttons = self._amr_filter_buttons()

        html += f"""
        {caveat}
        <div class="alert-box alert-info"><i class="fas fa-biohazard"></i><div><h3>AMR Genes (Carbapenemases, ESBLs, Colistin, Tigecycline)</h3><p>Each gene with frequency (count & %). Use filters below to focus on key resistance classes. Use the genome search to highlight isolates carrying specific genes. Carbapenemases are highlighted as critical.</p><p><strong>Focused on A. baumannii:</strong> OXA-type carbapenemases (OXA-23, OXA-58), NDM, VIM, IMP, KPC, and colistin resistance mechanisms (mcr, pmrAB).</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('amr-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-success" onclick="exportTableToCSV('amr-table','amr_genes.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        {group_controls}
        {search_boxes}
        <div class="action-buttons">{filter_buttons}</div>
        <div class="master-scrollable-container"><table id="amr-table" class="data-table"><thead><tr><th data-sort="string">Gene</th><th data-sort="string">Database</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead><tbody>
        """
        for g in all_genes:
            genome_tags = self._make_genome_tags(g['genomes'])
            gene_display = f"<strong>{g['gene']}</strong>" + (' 🔥' if g['category']=='Carbapenemases' else '')
            html += f'<tr><td>{gene_display}</td><td>{g["database"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'

        # v3.2.0: Educational content moved into collapsible expanders
        html += f"""
        <h3 style="margin-top:30px;">📚 Educational Context</h3>
        <p style="color:#666;font-size:.9em;margin-bottom:10px;">Click any section to expand. Collapsed by default so the gene table stays in focus.</p>
        <details class="edu-details">
            <summary>🔬 Understanding Database-Specific Gene Reporting</summary>
            <div>{educational_box}</div>
        </details>
        <details class="edu-details">
            <summary>🗃️ Database Roles & Reporting Logic</summary>
            <div>{edu_db_roles}</div>
        </details>
        """
        return html

    # --------------------------------------------------------------------------
    # SECTION: VIRULENCE 
    # --------------------------------------------------------------------------
    def _virulence_section(self, kwargs):
        gene_centric = kwargs['gene_centric']
        all_genes = []
        for db in gene_centric.get('virulence_databases', {}).values():
            all_genes.extend(db)
        all_genes.sort(key=lambda x: x['count'], reverse=True)

        html = """
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
        group_controls = self._grouping_controls('vir-table')
        search_boxes = self._search_boxes('vir', 'vir-table')
        filter_buttons = self._virulence_filter_buttons()
        html += f"""
        <div class="alert-box alert-info"><i class="fas fa-virus"></i><div><h3>Virulence Factors (Biofilm, Adhesion, Toxins)</h3><p>Key A. baumannii virulence genes: biofilm formation (ompA, csu, bfmRS), pili, and iron uptake. Biofilm formation is critical for persistence on medical devices and hospital outbreaks. Use filters to explore specific mechanisms.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('virulence-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-success" onclick="exportTableToCSV('vir-table','virulence_genes.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        {group_controls}
        {search_boxes}
        <div class="action-buttons">{filter_buttons}</div>
        <div class="master-scrollable-container"><table id="vir-table" class="data-table"><thead><tr><th data-sort="string">Gene</th><th data-sort="string">Database</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead><tbody>
        """
        for g in all_genes:
            genome_tags = self._make_genome_tags(g['genomes'])
            html += f'<tr><td><strong>{g["gene"]}</strong></td><td>{g["database"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html

    # --------------------------------------------------------------------------
    # SECTION: BACMET 
    # --------------------------------------------------------------------------
    def _bacmet_section(self, kwargs):
        gene_centric = kwargs['gene_centric']
        all_genes = []
        for db in gene_centric.get('bacmet_databases', {}).values():
            all_genes.extend(db)
        all_genes.sort(key=lambda x: x['count'], reverse=True)

        html = """
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
        group_controls = self._grouping_controls('bac-table')
        search_boxes = self._search_boxes('bac', 'bac-table')
        filter_buttons = self._bacmet_filter_buttons()
        html += f"""
        <div class="alert-box alert-info"><i class="fas fa-flask"></i><div><h3>BACMET2 Database: Biocide & Heavy Metal Resistance</h3><p>These genes confer resistance to disinfectants (quaternary ammonium, chlorhexidine) and heavy metals (mercury, copper, arsenic, etc.), which can co‑select for antibiotic resistance in hospital environments. Tracking these markers helps understand persistence and resistance spread.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('bacmet-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-success" onclick="exportTableToCSV('bac-table','bacmet_genes.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        {group_controls}
        {search_boxes}
        <div class="action-buttons">{filter_buttons}</div>
        <div class="master-scrollable-container"><table id="bac-table" class="data-table"><thead><tr><th data-sort="string">Gene</th><th data-sort="string">Database</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th><tr></thead><tbody>
        """
        for g in all_genes:
            genome_tags = self._make_genome_tags(g['genomes'])
            html += f'<tr><td><strong>{g["gene"]}</strong></td><td>{g["database"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html

    # --------------------------------------------------------------------------
    # SECTION: PLASMID 
    # --------------------------------------------------------------------------
    def _plasmid_section(self, kwargs, apt_plasmid):
        plasmidfinder = kwargs.get('plasmid_analysis', {})
        total_genomes = len(apt_plasmid.get('per_genome', {})) if apt_plasmid else 0

        html = """
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

        html += """
        <div class="alert-box alert-info"><i class="fas fa-dna"></i><div><h3>Plasmid Analysis – APT (Acinetobacter Plasmid Typing)</h3>
        <p>Species‑specific plasmid typing using rep genes (thresholds: identity ≥ 95%, subject coverage ≥ 70%). 
        Each rep type is assigned to a family (Rep1, Rep3, RepPriCT). This table shows the frequency of each rep type across the dataset, 
        with the list of genomes carrying it.</p>
        </div></div>
        """
        if apt_plasmid and apt_plasmid.get('rep_frequency'):
            group_controls_apt = self._grouping_controls('apt-rep-table')
            search_boxes_apt = self._search_boxes('apt-rep', 'apt-rep-table')
            html += f"""
            <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('plasmid-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-success" onclick="exportTableToCSV('apt-rep-table','apt_plasmid_rep_freq.csv')"><i class="fas fa-download"></i> Export APT Rep Frequency</button></div>
            {group_controls_apt}
            {search_boxes_apt}
            <div class="master-scrollable-container">
            <table id="apt-rep-table" class="data-table">
                <thead><tr><th data-sort="string">Marker</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead>
                <tbody>
            """
            for item in apt_plasmid['rep_frequency']:
                rep_type = item['rep_type']
                count = item['count']
                percentage = (count / total_genomes) * 100 if total_genomes else 0
                freq_disp = f"{count} ({percentage:.1f}%)"
                tags = self._make_genome_tags(item['genomes'])
                html += f'<tr><td><strong>{rep_type}</strong></td><td>{freq_disp}</td><td class="col-genomes"><div class="genome-list">{tags}</div></td></tr>'
            html += "</tbody></table></div>"
            html += """
            <h3>Per‑Genome Rep Types (APT)</h3>
            <div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-apt-genome" onkeyup="searchTable('apt-genome-table','search-apt-genome')" placeholder="🔍 Filter genome..."><input type="text" class="search-box" id="highlight-apt-genome" onkeyup="highlightGenome('apt-genome-table','highlight-apt-genome')" placeholder="🔍 Highlight genome or rep types..."></div>
            <div class="master-scrollable-container">
            <table id="apt-genome-table" class="data-table">
                <thead><tr><th data-sort="string">Genome</th><th data-sort="string">Rep Types</th></tr></thead>
                <tbody>
            """
            for genome, reps in sorted(apt_plasmid['per_genome'].items()):
                if reps:
                    badges = ' '.join(f'<span class="genome-tag">{r}</span>' for r in reps)
                    html += f'<tr><td><strong>{genome}</strong></td><td>{badges}</td></tr>'
                else:
                    html += f'<tr><td><strong>{genome}</strong></td><td><em>No Rep found</em></td></tr>'
            html += "</tbody></table></div>"
        else:
            html += "<p>No APT plasmid summary data found. Run the AcinetoScope plasmid typing module first.</p>"

        if plasmidfinder.get('plasmid_frequencies'):
            group_controls_pf = self._grouping_controls('plasmidfinder-table')
            search_boxes_pf = self._search_boxes('plasmidfinder', 'plasmidfinder-table')
            filter_buttons_pf = self._plasmid_filter_buttons()
            html += f"""
            <h3>Legacy PlasmidFinder Results (for reference)</h3>
            <p>Generalist database – not optimised for Acinetobacter. For species‑specific analysis, the APT results above are recommended.</p>
            {group_controls_pf}
            {search_boxes_pf}
            <div class="action-buttons">{filter_buttons_pf}</div>
            <div class="master-scrollable-container">
            <table id="plasmidfinder-table" class="data-table">
                <thead><tr><th data-sort="string">Marker</th><th data-sort="string">Category</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead>
                <tbody>
            """
            for g in plasmidfinder.get('plasmid_frequencies', []):
                count = g.get('count', 0)
                percentage = (count / total_genomes) * 100 if total_genomes else 0
                freq_disp = f"{count} ({percentage:.1f}%)"
                tags = self._make_genome_tags(g.get('genomes', []))
                html += f'<tr><td><strong>{g["plasmid_marker"]}</strong></td><td>{g.get("category", "Other")}</td><td>{freq_disp}</td><td class="col-genomes"><div class="genome-list">{tags}</div></td></tr>'
            html += "</tbody></table></div>"
        else:
            html += "<p>No PlasmidFinder data available.</p>"
        return html

    # --------------------------------------------------------------------------
    # SECTION: MUTATIONS — v3.2.0: Mutation and Subclass columns removed
    # --------------------------------------------------------------------------
    def _mutation_section(self, kwargs):
        mutation_data = kwargs.get('mutation_data', {})
        mutations = mutation_data.get('mutations', [])
        if not mutations:
            return """
            <div class="alert-box alert-warning"><i class="fas fa-dna"></i><div><h3>No Mutation Data</h3><p>mutation_summary.html not found or could not be parsed.</p></div></div>
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

        total_samples = len(kwargs['samples_data'])
        group_controls = self._grouping_controls('mutation-table')
        search_boxes = self._search_boxes('mutation', 'mutation-table')
        filter_buttons = self._mutation_filter_buttons()

        html += f"""
        <div class="alert-box alert-info"><i class="fas fa-dna"></i><div><h3>Point Mutations (AMRFinderPlus)</h3><p>Each gene/class combination is shown with all genomes that carry mutations in it. Mutations in gyrA/parC (fluoroquinolones), pmrAB (colistin), or 23S rRNA (tigecycline) are clinically relevant in A. baumannii.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('mutation-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-success" onclick="exportTableToCSV('mutation-table','mutations.csv')"><i class="fas fa-download"></i> Export Mutations</button></div>
        {group_controls}
        {search_boxes}
        <div class="action-buttons">{filter_buttons}</div>
        <div class="master-scrollable-container"><table id="mutation-table" class="data-table" style="table-layout:fixed;width:100%;">
            <colgroup>
                <col style="width:180px;">
                <col style="width:320px;">
                <col style="width:160px;">
                <col>
            </colgroup>
            <thead><tr>
                <th data-sort="string">Gene</th>
                <th data-sort="string">Class</th>
                <th data-sort="number">Count</th>
                <th data-sort="string">Genomes</th>
            </tr></thead><tbody>
        """
        for m in mutations:
            gene = m['gene']
            class_name = m['class']
            genomes = m['genomes']
            count = len(genomes)
            genome_tags = self._make_genome_tags(genomes)
            pct = (count / total_samples * 100) if total_samples else 0
            html += f'<tr><td><strong>{gene}</strong></td><td>{class_name}</td><td>{count} ({pct:.1f}%)</td><td><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html

    # --------------------------------------------------------------------------
    # SECTION: MGE — NEW in v3.2.0
    # --------------------------------------------------------------------------
    def _mge_section(self, kwargs):
        mge = kwargs.get('mge_data', {})
        stats = mge.get('stats', {})
        per_sample = mge.get('per_sample', {})
        note = mge.get('note', '')

        if not per_sample:
            return """
            <div class="alert-box alert-warning"><i class="fas fa-mobile-alt"></i><div><h3>No MGE Data</h3><p>acineto_mge_summary.html not found or could not be parsed.</p></div></div>
            """

        html = """
        <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #e0f2f1 100%); border-left: 6px solid #16A085; margin-bottom: 20px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">🧬</span>
                <div>
                    <strong style="font-size: 1.2em; color: #006064;">Mobile Genetic Elements – mobileOG-db</strong><br>
                    <span style="font-size: 0.95em; color: #333;">
                        <strong>mobileOG-db</strong> developed by 
                        <a href="https://github.com/clb21565/mobileOG-db" target="_blank" style="color: #16A085; font-weight: bold;">Brown, Mullet and colleagues</a>
                        (Brown CL, et al. Appl Environ Microbiol. 2022;88(18):e00991-22.
                        <a href="https://doi.org/10.1128/aem.00991-22" target="_blank" style="color: #16A085; font-weight: 600;">doi:10.1128/aem.00991-22</a>)
                        <i class="fas fa-arrow-right"></i> Manually curated database of protein families mediating the life cycle of bacterial mobile genetic elements.
                        <br>
                        <span style="color: #6c757d;">
                            <i class="fas fa-gratipay"></i> We are grateful to the mobileOG-db curators and contributors for making this resource freely available.
                        </span>
                    </span>
                </div>
            </div>
        </div>
        """

        cards = ''
        cards += self._stat_card(stats.get('Total Samples', len(per_sample)), 'Total Samples', '#4CAF50', 'fa-vial')
        cards += self._stat_card(stats.get('With mobileOG hits', '—'), 'With mobileOG hits', '#2196F3', 'fa-check')
        cards += self._stat_card(stats.get('Without mobileOG hits', '—'), 'Without mobileOG hits', '#9E9E9E', 'fa-times')
        cards += self._stat_card(stats.get('Total mobileOG hits', '—'), 'Total mobileOG hits', '#FF9800', 'fa-database')
        cards += self._stat_card(stats.get('Key MGE-associated genes', '—'), 'Key MGE-assoc. genes', '#8B5CF6', 'fa-star')
        cards += self._stat_card(stats.get('Runtime', '—'), 'Runtime', '#E53935', 'fa-clock')

        html += f'<div class="stats-grid">{cards}</div>'

        if note:
            html += f'<div class="alert-box alert-warning"><i class="fas fa-info-circle fa-2x"></i><div><strong>Interpretation note:</strong> {note}</div></div>'

        html += """
        <div class="alert-box alert-info"><i class="fas fa-dna"></i><div>
            <h3>Why Mobile Genetic Elements Matter</h3>
            <p>MGEs are the physical vehicles of horizontal gene transfer. In <em>A. baumannii</em>, resistance to carbapenems (OXA-23), aminoglycosides, and sulfonamides frequently travels on transposons, integrons, and conjugative plasmids. Higher mobileOG hit counts suggest a genome with more MGE machinery — potentially more permissive for acquiring new resistance determinants.</p>
        </div></div>
        """

        html += '<h3 style="margin-top:30px;">📚 Understanding the Categories</h3>'
        html += self._mge_category_cards()

        # Category totals table
        all_cols = ['Sample']
        seen = set()
        for s in per_sample.values():
            for k in s.keys():
                if k not in seen:
                    seen.add(k)
                    all_cols.append(k)
        numeric_cols = [c for c in all_cols[1:]
                        if all(isinstance(s.get(c, 0), (int, float)) for s in per_sample.values())]
        cat_rows = ''
        grand_total = 0
        cat_data = []
        for col in numeric_cols:
            total = sum(per_sample[s].get(col, 0) for s in per_sample)
            grand_total += total
            cat_data.append((col, total))
        for col, total in sorted(cat_data, key=lambda x: x[1], reverse=True):
            pct = (total / grand_total * 100) if grand_total else 0
            mean = total / len(per_sample) if per_sample else 0
            mx_sample = max(per_sample, key=lambda s: per_sample[s].get(col, 0))
            mx_val = per_sample[mx_sample].get(col, 0)
            cat_rows += (f'<tr><td><strong>{col}</strong></td><td>{total:,}</td>'
                         f'<td>{pct:.1f}%</td><td>{mean:.2f}</td>'
                         f'<td>{mx_val}</td><td>{mx_sample}</td></tr>')
        html += f"""
        <h3 style="margin-top:30px;">📊 Category Totals</h3>
        <div class="master-scrollable-container sticky-y"><table class="data-table"><thead><tr>
            <th data-sort="string">Category</th>
            <th data-sort="number">Total hits</th>
            <th data-sort="number">% of all</th>
            <th data-sort="number">Mean/sample</th>
            <th data-sort="number">Max/sample</th>
            <th data-sort="string">Sample with max</th>
        </tr></thead><tbody>{cat_rows}</tbody></table></div>
        """

        # Per-sample table
        header_cells = ''.join(f'<th data-sort="{"string" if c=="Sample" else "number"}">{c}</th>' for c in all_cols)
        sample_rows = ''
        for sample, vals in sorted(per_sample.items()):
            cells = f'<td><strong>{sample}</strong></td>'
            for c in all_cols[1:]:
                cells += f'<td>{vals.get(c, 0)}</td>'
            sample_rows += f'<tr>{cells}</tr>'

        html += f"""
        <h3 style="margin-top:30px;">🧬 Per-Sample mobileOG Profile</h3>
        <input type="text" class="search-box" id="search-mge" onkeyup="searchTable('mge-table','search-mge')" placeholder="🔍 Filter sample...">
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="exportTableToCSV('mge-table','mge_profile.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        <div class="master-scrollable-container sticky-y">
            <table id="mge-table" class="data-table">
                <thead><tr>{header_cells}</tr></thead>
                <tbody>{sample_rows}</tbody>
            </table>
        </div>
        """

        # Aggregate by typing
        agg_html = self._build_mge_aggregate(per_sample, kwargs['samples_data'])
        html += f'<h3 style="margin-top:30px;">📈 Aggregate MGE Profile by Typing</h3>{agg_html}'
        return html

    def _stat_card(self, value, label: str, color: str = '#4CAF50', icon: str = '') -> str:
        return f'''
        <div class="stat-card" style="background:linear-gradient(135deg,{color} 0%,{color}dd 100%);">
            {f'<i class="fas {icon} fa-2x" style="opacity:0.9;margin-bottom:8px;"></i>' if icon else ''}
            <div class="stat-value">{value}</div>
            <div class="stat-label">{label}</div>
        </div>'''

    def _mge_category_cards(self) -> str:
        cards = [
            ('Integr.', 'Integrase', '#ef4444', 'fa-cut',
             'Integrases catalyse integration/excision of mobile elements.',
             'High counts suggest active chromosomal integration.'),
            ('Transfer', 'Transfer / Conjugation', '#f97316', 'fa-share-alt',
             'T4SS, relaxases, coupling proteins — the machinery of HGT.',
             'The route by which resistance spreads between bacteria.'),
            ('Stab.', 'Stability / Maintenance', '#eab308', 'fa-anchor',
             'Partition systems, toxin-antitoxin cassettes, addiction modules.',
             'Plasmids with strong stability systems persist without selection pressure.'),
            ('Phage', 'Bacteriophage', '#22c55e', 'fa-virus',
             'Phage structural and lysogeny proteins.',
             'Phages can carry virulence genes via transduction.'),
            ('Repl.', 'Replication', '#14b8a6', 'fa-sync',
             'Replication initiators, primases, helicases.',
             'Elements need a functional replication module to persist extrachromosomally.'),
            ('IS-assoc.', 'Insertion Sequences', '#0891b2', 'fa-exchange-alt',
             'Transposases and IS elements — the smallest autonomous MGEs.',
             'IS elements frequently disrupt genes or carry resistance cassettes.'),
            ('ICE-assoc.', 'Integrative Conjugative Elements', '#3b82f6', 'fa-project-diagram',
             'Self-transmissible elements that integrate into the chromosome.',
             'ICEs move large cassettes (e.g., AbaR resistance islands) between lineages.'),
            ('Plasmid-assoc.', 'Plasmid', '#8b5cf6', 'fa-plug',
             'Plasmid backbone, partitioning, conjugal transfer proteins.',
             'The single most important vehicle for acquired resistance in A. baumannii.'),
            ('Phage-assoc.', 'Prophage / Lysogeny', '#d946ef', 'fa-dna',
             'Prophage integrases, repressors, lysis modules.',
             'Prophage carriage correlates with virulence gene acquisition.'),
            ('Key MGE-assoc.', 'Key MGE signatures', '#dc2626', 'fa-star',
             'High-confidence subset of mobileOG-db entries.',
             'The cleanest single indicator of MGE activity per genome.'),
        ]
        html = '<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:15px;margin:20px 0;">'
        for key, label, color, icon, definition, importance in cards:
            html += f'''
            <div style="background:#fff;padding:16px;border-radius:10px;border-left:4px solid {color};box-shadow:0 2px 8px rgba(0,0,0,.06);">
                <div style="display:flex;align-items:center;gap:10px;margin-bottom:8px;">
                    <i class="fas {icon}" style="color:{color};font-size:1.4em;"></i>
                    <div>
                        <strong style="color:{color};font-size:1.05em;">{label}</strong>
                        <div style="font-size:.8em;color:#888;font-family:monospace;">{key}</div>
                    </div>
                </div>
                <p style="font-size:.88em;color:#333;margin:0 0 8px 0;">{definition}</p>
                <p style="font-size:.85em;color:#555;margin:0;padding-top:8px;border-top:1px dashed #e0e0e0;">
                    <i class="fas fa-bullseye" style="color:{color};"></i> <em>{importance}</em>
                </p>
            </div>'''
        html += '</div>'
        return html

    def _build_mge_aggregate(self, per_sample, samples_data):
        if not samples_data or not per_sample:
            return ''
        fields = [('pasteur_st', 'Pasteur ST'), ('k_locus', 'K Locus'),
                  ('o_locus', 'OCL Locus'), ('capsule', 'Capsule')]
        sample_cols = list(next(iter(per_sample.values())).keys())
        blocks = ''
        for field_key, field_label in fields:
            groups = defaultdict(list)
            for sample, vals in per_sample.items():
                sd = samples_data.get(sample, {})
                if field_key in ('pasteur_st', 'oxford_st'):
                    typ = sd.get('pasteur_mlst', {}).get('ST', 'ND')
                elif field_key == 'capsule':
                    typ = sd.get('kaptive', {}).get('Capsule_Type', 'ND')
                elif field_key == 'k_locus':
                    typ = sd.get('kaptive', {}).get('K_Locus', 'ND')
                elif field_key == 'o_locus':
                    typ = sd.get('kaptive', {}).get('O_Locus', 'ND')
                else:
                    typ = 'ND'
                if typ in ('Not Assigned', 'ND', ''):
                    continue
                groups[typ].append(vals)
            if not groups:
                continue
            rows = ''
            for grp, recs in sorted(groups.items(), key=lambda x: -len(x[1])):
                means = {}
                for c in sample_cols:
                    vv = [r.get(c, 0) for r in recs if isinstance(r.get(c, 0), (int, float))]
                    means[c] = round(sum(vv) / len(vv), 1) if vv else 0
                cells = ''.join(f'<td>{means[c]}</td>' for c in sample_cols)
                rows += f'<tr><td><strong>{grp}</strong></td><td>{len(recs)}</td>{cells}</tr>'
            header = ''.join(f'<th>{c}</th>' for c in sample_cols)
            blocks += f'''
            <h4 style="margin-top:20px;">Mean MGE profile by {field_label}</h4>
            <div class="master-scrollable-container sticky-y"><table class="data-table"><thead><tr>
                <th>{field_label}</th><th>N samples</th>{header}
            </tr></thead><tbody>{rows}</tbody></table></div>'''
        return blocks

    # --------------------------------------------------------------------------
    # SECTION: CO-OCCURRENCE 
    # --------------------------------------------------------------------------
    def _cooccurrence_section(self, kwargs):
        cooc = kwargs['patterns'].get('gene_cooccurrence', {})
        if not cooc:
            return '<div class="alert-box alert-warning">No co‑occurrence data available.</div>'
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
        <div class="master-scrollable-container sticky-y"><table id="cooccurrence-table" class="data-table"><thead><tr><th data-sort="string">Gene 1</th><th data-sort="string">Gene 2</th><th data-sort="number">Co‑occurrence Count</th></tr></thead><tbody>
        """
        for g1, g2, cnt in top:
            html += f'<tr><td>{g1}</td><td>{g2}</td><td>{cnt}</td></tr>'
        html += '</tbody></table></div>'
        return html

    # --------------------------------------------------------------------------
    # SECTION: PATTERNS 
    # --------------------------------------------------------------------------
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
            html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'patterns-tab\')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById(\'search-patterns\').value=\'\'; searchTable(\'crab-table\',\'search-patterns\'); document.getElementById(\'highlight-patterns\').value=\'\'; highlightTableCells(\'crab-table\',\'highlight-patterns\')"><i class="fas fa-sync"></i> Clear</button></div>'
            html += '<div style="display:flex; gap:10px;"><input type="text" class="search-box" id="search-patterns" onkeyup="searchTable(\'crab-table\',\'search-patterns\')" placeholder="🔍 Filter rows..."><input type="text" class="search-box" id="highlight-patterns" onkeyup="highlightTableCells(\'crab-table\',\'highlight-patterns\')" placeholder="✨ Highlight text..."></div>'
            html += '<div class="master-scrollable-container sticky-y"><table id="crab-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">ST</th><th data-sort="string">Capsule (K:OCL)</th><th data-sort="string">Carbapenemases</th><th data-sort="string">Colistin Resistance</th><th data-sort="string">Tigecycline Resistance</th></tr></thead><tbody>'
            for c in crab:
                html += f'<tr><td><strong>{c["sample"]}</strong></td><td>{c["pasteur_st"]}</td><td>{c["capsule_type"]}</td><td>{",".join(c["carbapenemases"])}</td><td>{",".join(c["colistin_resistance"])}</td><td>{",".join(c["tigecycline_resistance"])}</td></tr>'
            html += '</tbody></table></div>'
        else:
            html += '<p>No high‑risk CRAB isolates detected.</p>'
        return html

    # --------------------------------------------------------------------------
    # SECTION: COMPARE — NEW in v3.2.0
    # --------------------------------------------------------------------------
    def _precompute_similarity(self, compare_data: Dict[str, Dict]) -> Dict[str, Any]:
        samples = sorted(compare_data.keys())
        n = len(samples)
        if n < 2:
            return {'samples': samples, 'flat': [], 'n': n}
        cats = ('AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations')
        gene_sets, typing_data = {}, {}
        for s in samples:
            d = compare_data[s]
            gene_sets[s] = {cat: set(d['genes'].get(cat, [])) for cat in cats}
            typing_data[s] = d['typing']
        typing_keys = list(typing_data[samples[0]].keys()) if samples else []
        n_keys = len(typing_keys)
        flat = []
        for i in range(n):
            s1 = samples[i]
            for j in range(i + 1, n):
                s2 = samples[j]
                matches = sum(1 for k in typing_keys if typing_data[s1].get(k) == typing_data[s2].get(k))
                tp = (matches / n_keys * 100) if n_keys else 0
                cs = []
                for cat in cats:
                    sa, sb = gene_sets[s1][cat], gene_sets[s2][cat]
                    union = sa | sb
                    if not union:
                        continue
                    cs.append(len(sa & sb) / len(union) * 100)
                gp = (sum(cs) / len(cs)) if cs else 0
                flat.append(round(tp * 0.5 + gp * 0.5, 1))
        return {'samples': samples, 'flat': flat, 'n': n}

    def _compare_section(self, kwargs):
        samples_data = kwargs['samples_data']
        mutation_data = kwargs.get('mutation_data', {})
        mutation_list = mutation_data.get('mutations', [])

        if len(samples_data) < 2:
            return self._simple_alert('warning', 'fa-balance-scale',
                f'<h3>Not Enough Samples</h3><p>Compare requires ≥2 samples (have {len(samples_data)}).</p>')

        sample_mutations = defaultdict(list)
        for m in mutation_list:
            label = f"{m.get('gene', '')} {m.get('mutation', '')}".strip()
            for genome in m.get('genomes', []):
                sample_mutations[genome].append(label)

        def db_genes(sample, db):
            return list(samples_data.get(sample, {}).get('abricate_databases', {}).get(db, []))

        compare_data = {}
        for sample, sd in samples_data.items():
            amr_genes = list(sd.get('amr_genes', []))
            for db in ('card', 'resfinder', 'argannot', 'megares', 'ncbi'):
                amr_genes += db_genes(sample, db)
            amr_genes = sorted(set(amr_genes))
            vir_genes = sorted(set(db_genes(sample, 'vfdb') + db_genes(sample, 'ecoli_vf')))
            bacmet_genes = sorted(set(db_genes(sample, 'bacmet2')))
            plasmid_genes = sorted(set(sd.get('plasmid_genes', []) + db_genes(sample, 'plasmidfinder')))
            muts = sorted(set(sample_mutations.get(sample, [])))
            p = sd.get('pasteur_mlst', {})
            o = sd.get('oxford_mlst', {})
            kd = sd.get('kaptive', {})
            compare_data[sample] = {
                'typing': {
                    'Pasteur ST': p.get('ST', 'ND'),
                    'International Clone': p.get('International_Clone', 'ND'),
                    'Oxford ST': o.get('ST', 'ND'),
                    'K Locus': kd.get('K_Locus', 'ND'),
                    'OCL Locus': kd.get('O_Locus', 'ND'),
                    'Capsule': kd.get('Capsule_Type', 'ND'),
                },
                'genes': {
                    'AMR': amr_genes, 'Virulence': vir_genes,
                    'BACMET': bacmet_genes, 'Plasmids': plasmid_genes,
                    'Mutations': muts,
                },
            }

        sorted_samples = sorted(compare_data.keys())
        options_a = ''.join(f'<option value="{s}"{" selected" if i == 0 else ""}>{s}</option>'
                            for i, s in enumerate(sorted_samples))
        options_b = ''.join(f'<option value="{s}"{" selected" if i == 1 else ""}>{s}</option>'
                            for i, s in enumerate(sorted_samples))

        def safe_json(obj):
            return (json.dumps(obj, default=str, ensure_ascii=False)
                    .replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026'))

        blob_json = safe_json(compare_data)
        precomp = self._precompute_similarity(compare_data)
        precomp_json = safe_json(precomp)

        return f"""
        <div class="alert-box alert-info"><i class="fas fa-balance-scale fa-2x"></i><div>
            <h3>⚖️ Compare &amp; Cluster Analysis</h3>
            <p>Two modes: <strong>Pair Compare</strong> for a deep side-by-side of two isolates, and <strong>Cluster Detection</strong> to find groups of highly similar samples across the whole dataset.</p>
            <ul style="margin-top:6px;">
                <li><strong>Outbreak confirmation</strong> — do two isolates share ST + K locus + OCL + capsule?</li>
                <li><strong>Discordance analysis</strong> — what changed between an ancestor and a resistant descendant?</li>
                <li><strong>Cluster detection</strong> — which samples form tight transmission clusters?</li>
            </ul>
        </div></div>

        <div class="mode-toggle">
            <button class="mode-btn active" data-mode="pair" onclick="setCompareMode('pair')">
                <i class="fas fa-balance-scale"></i> Pair Compare
            </button>
            <button class="mode-btn" data-mode="cluster" onclick="setCompareMode('cluster')">
                <i class="fas fa-project-diagram"></i> Cluster Detection
            </button>
        </div>

        <div id="mode-pair" class="compare-mode active">
            <div class="compare-panel">
                <div class="compare-select-row">
                    <label for="compare-a"><span>Sample A</span><select id="compare-a">{options_a}</select></label>
                    <button class="action-btn btn-primary" onclick="runCompare()">
                        <i class="fas fa-balance-scale"></i> Compare</button>
                    <button class="action-btn btn-light" onclick="swapCompare()">
                        <i class="fas fa-exchange-alt"></i> Swap</button>
                    <button class="action-btn btn-light" onclick="resetCompare()">
                        <i class="fas fa-sync"></i> Reset</button>
                    <label for="compare-b"><span>Sample B</span><select id="compare-b">{options_b}</select></label>
                </div>
                <div id="compare-verdict" style="display:none;"></div>
                <div id="compare-visuals"></div>
                <div id="compare-output"></div>
            </div>
        </div>

        <div id="mode-cluster" class="compare-mode">
            <div class="compare-panel">
                <div class="cluster-controls">
                    <label><span>Similarity threshold:</span>
                        <select id="cluster-threshold" onchange="runCluster()">
                            <option value="95">Tight clusters (≥95%)</option>
                            <option value="90" selected>Moderate clusters (≥90%)</option>
                            <option value="85">Loose clusters (≥85%)</option>
                            <option value="80">Broad clusters (≥80%)</option>
                        </select>
                    </label>
                    <button class="action-btn btn-primary" onclick="runCluster()">
                        <i class="fas fa-project-diagram"></i> Detect Clusters</button>
                    <span id="cluster-summary" class="results-counter"></span>
                </div>
                <div id="cluster-clusters"></div>
                <div id="cluster-heatmap"></div>
            </div>
        </div>

        <script>
        window.STAPHSCOPE_COMPARE_DATA = {blob_json};
        window.STAPHSCOPE_SIM = {precomp_json};
        </script>"""

    def _simple_alert(self, kind, icon, html_body):
        return f'<div class="alert-box alert-{kind}"><i class="fas {icon} fa-2x"></i><div>{html_body}</div></div>'

    # --------------------------------------------------------------------------
    # SECTION: AI GUIDE 
    # --------------------------------------------------------------------------
    def _aiguide_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-robot fa-2x"></i><div><h3>AI Assistant Guide – Unleash the Power of AI for Genomic Epidemiology</h3><p>This guide shows you how to leverage large language models (LLMs) like <strong>ChatGPT, Claude, or Gemini</strong> to interact with your <em>A. baumannii</em> dataset – turning static reports into dynamic conversations.</p></div></div>
        <div style="margin: 20px 0;">
            <div class="database-section"><h4><i class="fas fa-brain"></i> Why Use AI for Genomic Data Analysis?</h4><p>Modern AI models excel at:</p><ul><li><strong>Pattern recognition</strong> – spotting epidemiological trends, clone associations, and co‑occurrence networks.</li><li><strong>Natural language queries</strong> – ask in plain English, get instant answers without writing code.</li><li><strong>Hypothesis generation</strong> – uncover unexpected correlations that merit experimental follow‑up.</li><li><strong>Literature synthesis</strong> – connect your findings with published resistance mechanisms and clinical guidelines.</li></ul><p><span style="background:#fff3cd; padding:2px 8px; border-radius:4px;"><i class="fas fa-lightbulb"></i> <strong>Scientific note:</strong> AI is <em>pattern‑finding</em>, not <em>causal‑inferring</em>. Use it to suggest, then verify with wet‑lab or clinical correlation.</span></p></div>
            <div class="database-section"><h4><i class="fas fa-upload"></i> How to Feed This Report to AI</h4><p>You have <strong>three powerful options</strong>:</p><ol><li><strong>Upload the JSON file</strong> – <code>genius_acinetobacter_ultimate_gene_centric_report.json</code> contains all structured data. Upload it to ChatGPT (Advanced Data Analysis), Claude, or Gemini. <em>Best for precise quantitative queries.</em></li><li><strong>Upload the HTML report</strong> – Modern AI tools can read HTML and extract tables. Upload the <code>.html</code> file directly – the AI will parse the tables and text. <em>Great for visual context.</em></li><li><strong>Copy‑paste specific tables</strong> – If you only need a quick insight, copy a table (e.g., AMR gene list) and paste it into the chat. <em>Instant, no file upload needed.</em></li></ol><p style="margin-top:10px; background:#e8f5e9; padding:10px; border-radius:5px;"><i class="fas fa-info-circle"></i> <strong>Pro tip:</strong> For best results, tell the AI: <em>"You are a bioinformatician analysing <em>A. baumannii</em> genomes. The attached data contains typing, AMR, virulence, BACMET, mutation, and MGE information. Answer my questions with references to the data."</em></p></div>
            <div class="database-section"><h4><i class="fas fa-chart-line"></i> Scientifically Relevant Questions to Ask</h4><div style="display:grid; grid-template-columns:1fr 1fr; gap:10px;"><div style="background:#f8f9fa; padding:10px; border-radius:8px;"><strong>🧬 Epidemiology &amp; Clonality</strong><ul style="margin-top:5px; font-size:0.9em;"><li>What are the most common Pasteur STs in this dataset?</li><li>Which capsule types (K:OCL) are dominant?</li><li>Are there any ST‑capsule combinations with >2 isolates?</li><li>Which clones carry the most resistance genes?</li></ul></div><div style="background:#f8f9fa; padding:10px; border-radius:8px;"><strong>💊 Antimicrobial Resistance</strong><ul style="margin-top:5px; font-size:0.9em;"><li>How many samples carry OXA-23? What are their STs?</li><li>Are there any NDM‑positive samples? What is their capsule type?</li><li>Which AMR genes co‑occur most frequently?</li><li>What is the distribution of colistin resistance genes (mcr, pmr)?</li></ul></div><div style="background:#f8f9fa; padding:10px; border-radius:8px;"><strong>🦠 Virulence &amp; Biofilm</strong><ul style="margin-top:5px; font-size:0.9em;"><li>Which samples carry biofilm genes (ompA, csu)?</li><li>List all isolates with the T6SS (tss, hcp).</li><li>Is there a correlation between capsule type and virulence gene carriage?</li></ul></div><div style="background:#f8f9fa; padding:10px; border-radius:8px;"><strong>🧪 Mutations &amp; Biocides</strong><ul style="margin-top:5px; font-size:0.9em;"><li>What are the most frequent point mutations in gyrA or parC?</li><li>Are there any 23S rRNA mutations (tigecycline resistance)?</li><li>Which samples carry qac genes (disinfectant resistance)?</li></ul></div></div><p style="margin-top:10px;"><i class="fas fa-arrow-right"></i> <strong>Beyond tables:</strong> Ask the AI to <em>"write a summary of the resistance profile for ST2"</em> or <em>"compare virulence gene carriage between K1 and K2 capsule types"</em> – the report contains all data.</p></div>
            <div class="database-section"><h4><i class="fas fa-balance-scale"></i> Scientific Rigour &amp; Ethical AI Use</h4><ul><li><strong>AI is your co‑pilot, not the pilot.</strong> Always interpret AI‑generated insights in the context of your local epidemiology, clinical guidelines, and laboratory validation.</li><li><strong>Verify, verify, verify.</strong> Cross‑check critical calls (e.g., carbapenemase presence, colistin resistance) with primary literature, genome browsers, or secondary tools.</li><li><strong>No patient‑identifiable data.</strong> Only upload aggregated, de‑identified genomic data. This report contains no clinical metadata.</li><li><strong>Transparency in publications.</strong> If you use AI for exploratory analysis, mention it in the methods (e.g., "AI‑assisted pattern discovery was performed using a large language model, followed by manual curation").</li><li><strong>AI hallucination is real.</strong> If the AI confidently tells you that <em>blaOXA-23</em> is found in <em>E. coli</em> or that tigecycline resistance is caused by a mutation in <em>rpoB</em> – <strong>don't believe it</strong>. Treat every AI statement as a hypothesis, not a fact.</li></ul></div>
            <div class="database-section" style="background:#fff3cd; border-left:6px solid #ffc107;"><h4><i class="fas fa-smile-wink"></i> A (Mostly Serious) AI Survival Guide</h4><ul><li><strong>If the AI says "I don't know"</strong> – trust it. It's being honest.</li><li><strong>If the AI says "It is widely known"</strong> – ask for a reference. It may have made it up.</li><li><strong>If the AI offers a completely novel evolutionary theory</strong> – check if your coffee is spiked.</li><li><strong>Remember:</strong> AI won't take your job – but a microbiologist who knows how to use AI might! So learn it, use it, and always keep a healthy dose of scepticism. 😉</li></ul><p style="margin-top:10px;"><i class="fas fa-microbe"></i> <strong>Final thought:</strong> The best AI‑human partnership is one where the AI does the heavy pattern‑lifting, and you do the heavy thinking. Happy (and careful) exploring!</p></div>
        </div>
        """

    # --------------------------------------------------------------------------
    # SECTION: CALL TO ACTION — NEW in v3.2.0
    # --------------------------------------------------------------------------
    def _calltoaction_section(self):
        return """
        <div class="alert-box alert-info">
            <i class="fas fa-globe fa-2x"></i>
            <div>
                <h3>The Global Burden of AMR and Our Call to Action</h3>
                <p>Antimicrobial resistance (AMR) is one of the top global public health threats, with an estimated <strong>1.27 million direct deaths annually</strong>. <em>Acinetobacter baumannii</em> is a major contributor — skin and soft tissue infections, bloodstream infections, pneumonia, endocarditis. Tracking AMR and virulence determinants is essential to inform treatment guidelines and infection control.</p>
                <p>We developed <strong>AcinetoScope</strong> to empower researchers and clinicians — especially in low-resource settings — to analyse their own sequencing data without extensive bioinformatics expertise.</p>
            </div>
        </div>
        <div style="background:#e8f5e9;padding:20px;border-radius:12px;margin:20px 0;">
            <h3><i class="fas fa-bacterium"></i> ESCAPE AMR – Our Ongoing Project (ESKAPE Pathogens)</h3>
            <p><strong>AcinetoScope</strong> is one of the first modules of a larger initiative called <strong>ESCAPE AMR</strong> (formerly ESKAPE). We target the notorious <strong>ESKAPE pathogens</strong>:</p>
            <ul>
                <li><strong>E</strong>nterococcus faecium</li>
                <li><strong>S</strong>taphylococcus aureus</li>
                <li><strong>K</strong>lebsiella pneumoniae</li>
                <li><strong>A</strong>cinetobacter baumannii</li>
                <li><strong>P</strong>seudomonas aeruginosa</li>
                <li><strong>E</strong>nterobacter species</li>
            </ul>
            <p>These bacteria "escape" the effects of antibiotics — hence the name. But we believe the name is also a global call to action:</p>
            <div style="background:#fff3e0;padding:15px;border-radius:8px;margin:15px 0;">
                <p><strong>🔹 E</strong>veryone must join forces — researchers, clinicians, policymakers, and citizens.<br>
                <strong>🔹 S</strong>mart surveillance is our first line of defence. No more guessing — we need genomic data.<br>
                <strong>🔹 K</strong>nowledge must be shared openly. No paywalls, no closed silos.<br>
                <strong>🔹 A</strong>frica bears a heavy AMR burden, but African solutions are already emerging.<br>
                <strong>🔹 P</strong>revention is cheaper than cure. Let's stop resistant infections before they spread.<br>
                <strong>🔹 E</strong>very day we delay, more lives are at stake. The time to act is now, not tomorrow.</p>
            </div>
            <p><i class="fas fa-laugh-squint"></i> <strong>"We didn't choose the name ESKAPE because it sounds cool (though it does). We chose it because it reminds us every single day: we must ESCAPE the AMR crisis — together, urgently, and with the best science we have."</strong><br>— Brown Beckley, lead developer (who secretly hopes this pun makes you smile, not roll your eyes 😉)</p>
        </div>
        <div style="text-align:center;margin:40px 0;">
            <i class="fas fa-star" style="font-size:3em;color:#ffc107;"></i>
            <h3>🤝 We Invite You to Contribute!</h3>
            <p><strong>If you find this tool useful, please:</strong></p>
            <div class="action-buttons" style="justify-content:center;">
                <a href="https://github.com/bbeckley-hub/acinetoscope" target="_blank" class="action-btn btn-primary" style="text-decoration:none;"><i class="fab fa-github"></i> ⭐ Star us on GitHub</a>
                <a href="mailto:brownbeckley94@gmail.com" class="action-btn btn-success" style="text-decoration:none;"><i class="fas fa-envelope"></i> Contact the team</a>
                <a href="https://github.com/bbeckley-hub/acinetoscope/issues" target="_blank" class="action-btn btn-warning" style="text-decoration:none;"><i class="fas fa-bug"></i> Report issues</a>
            </div>
            <p style="margin-top:20px;"><i class="fas fa-chalkboard-user"></i> <strong>We welcome collaborations</strong> to adapt this tool for other pathogens and to improve AMR surveillance in Africa and beyond.</p>
            <p><i class="fas fa-hand-holding-heart"></i> If you are a funder or organisation interested in supporting the <strong>ESCAPE AMR</strong> project, please reach out. Together we can build a free, open-source ecosystem for genomic surveillance of all ESKAPE pathogens.</p>
        </div>
        <div style="background:#f8f9fa;padding:15px;border-radius:8px;">
            <i class="fas fa-quote-left"></i> "AMR is a silent pandemic, but we have the tools to fight it — if we share them, if we teach each other, and if we act with urgency. Let's escape the era of untreatable infections."<br>
            <strong>— The ESCAPE AMR Team, University of Ghana Medical School</strong>
        </div>
        """

    # --------------------------------------------------------------------------
    # SECTION: CITATIONS — v3.2.0: AcinetoScope Sci Rep citation updated
    # --------------------------------------------------------------------------
    def _citation_section(self):
        citations = [
            {
                "title": "AcinetoScope",
                "text": "Beckley B, Amarh V, Lopes BS, et al. AcinetoScope: a computational pipeline for rapid, comprehensive <em>Acinetobacter baumannii</em> outbreak investigation and resistance gene tracking. <em>Sci Rep</em>. 2026. <a href='https://doi.org/10.1038/s41598-026-70189-y' target='_blank'>doi:10.1038/s41598-026-70189-y</a>",
                "url": "https://doi.org/10.1038/s41598-026-70189-y",
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
                "title": "mobileOG-db",
                "text": "Brown CL, Mullet J, Hindi F, Stoll JE, Gupta S, Choi M, Keenum I, Vikesland P, Pruden A, Zhang L. mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. Appl Environ Microbiol. 2022;88(18):e00991-22.",
                "url": "https://doi.org/10.1128/aem.00991-22",
                "category": "tool"
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
            badge_text = {
                'main': 'Main Tool',
                'tool': 'Tool',
                'db': 'Database',
                'other': 'Other'
            }[cat_class]
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
        "Genomic analysis was performed using AcinetoScope (Beckley et al., <em>Sci Rep</em> 2026; doi:10.1038/s41598-026-70189-y), which integrates MLST (Seemann, 2018) using the PubMLST database (Jolley et al., 2018), Kaptive 3 (Stanton et al., 2025), ABRicate (Seemann, 2024), AMRFinderPlus (Feldgarden et al., 2021), APT (Acinetobacter Plasmid Typing) (Lam et al., 2023), PlasmidFinder (Carattoli et al., 2014), and FastANI (Jain et al., 2018). Antimicrobial resistance genes were identified using CARD (Alcock et al., 2023), ResFinder (Bortolaia et al., 2020), MEGARes 3.0 (Bonin et al., 2023), and ARG-ANNOT (Gupta et al., 2014). For biocide and heavy metal resistance genes, BacMet (Pal et al., 2014) was used. Virulence screening was performed with ABRicate using VFDB (Chen et al., 2016). Mutation detection was performed using AMRFinderPlus. Mobile genetic element profiling was performed with mobileOG-db (Brown et al., 2022). FASTA QC was performed using Biopython (Cock et al., 2009)."
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

    # --------------------------------------------------------------------------
    # SECTION: FUNDING 
    # --------------------------------------------------------------------------
    def _funding_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-coffee fa-2x"></i><div><h3>Funding & Support – Keeping the Lights On (with code and caffeine)</h3><p>AcinetoScope is an <strong>independent, unfunded project</strong> born out of passion for genomic surveillance and AMR research at the University of Ghana Medical School.</p><p>No grants, no sponsors, no institutional backing – just a laptop, a lot of coffee, and a burning desire to help researchers fight antimicrobial resistance.</p></div></div>
        <div class="alert-box alert-warning"><i class="fas fa-heart fa-2x"></i><div><h3>💡 How You Can Help (Without Opening Your Wallet)</h3><ul><li><strong>⭐ Star us on GitHub</strong> – It takes two seconds and makes us feel like rockstars.</li><li><strong>🐛 Report bugs</strong> – If something breaks, let us know. We'll fix it with joy.</li><li><strong>💡 Suggest features</strong> – Have an idea? We're all ears.</li><li><strong>🧬 Share your data</strong> – If you've used the tool and want to collaborate, we'd love to hear your story.</li><li><strong>📢 Spread the word</strong> – Tell your colleagues, tweet about it, or mention it in your next Zoom call.</li><li><strong>👋 Say hello!</strong> – Seriously, just drop an email to <strong>brownbeckley94@gmail.com</strong>. It makes our day.</li></ul><p><i class="fas fa-microbe"></i> <strong>Fun fact:</strong> This project runs on 100% volunteer tears, 0% grant money. But we're not bitter – we're just caffeinated.</p></div></div>
        <div class="alert-box alert-success"><i class="fas fa-hand-holding-heart"></i><div><h3>🤝 Contribute to the ESKAPE AMR Platform</h3><p>We also maintain pipelines for other ESKAPE pathogens (Kleboscope, Pseudoscope, etc.). If you're a developer, bioinformatician, or just someone who loves clean code and bacteria, we welcome: pull requests, issues, documentation improvements, ideas for new databases. Visit our GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a> – star, fork, and let's fight AMR together!</p><p><strong>Brown Beckley</strong> – <i class="fas fa-envelope"></i> brownbeckley94@gmail.com</p><p><i class="fas fa-laugh-beam"></i> <strong>P.S.</strong> If you ever meet Brown in person, buy him a coffee – he'll probably talk your ear off about ST types, but it's worth it.</p></div></div>
        """

    # --------------------------------------------------------------------------
    # SECTION: EXPORT — added MGE CSV in v3.2.0
    # --------------------------------------------------------------------------
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
            <button class="action-btn btn-primary" onclick="exportTableToCSV('mge-table','mge_profile.csv')">MGE Profile CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('cooccurrence-table','cooccurrence.csv')">Co‑occurrence CSV</button>
            <button class="action-btn btn-success" onclick="location.href='genius_acinetobacter_ultimate_gene_centric_report.json'">Download JSON</button>
        </div>
        """

    # --------------------------------------------------------------------------
    # Helper: Grouping controls HTML 
    # --------------------------------------------------------------------------
    def _grouping_controls(self, table_id: str) -> str:
        return f"""
        <div class="grouping-controls">
            <strong><i class="fas fa-layer-group"></i> Group genomes by:</strong>
            <button class="group-btn" data-group="pasteur_st" onclick="groupGenomesByTyping('{table_id}', 'pasteur_st')">ST</button>
            <button class="group-btn" data-group="k_locus" onclick="groupGenomesByTyping('{table_id}', 'k_locus')">K Locus</button>
            <button class="group-btn" data-group="o_locus" onclick="groupGenomesByTyping('{table_id}', 'o_locus')">OCL Locus</button>
            <button class="group-btn" data-group="capsule" onclick="groupGenomesByTyping('{table_id}', 'capsule')">Capsule</button>
            <button class="group-btn" data-group="st_k" onclick="groupGenomesByTyping('{table_id}', 'st_k')">ST‑K</button>
            <button class="group-btn" data-group="st_o" onclick="groupGenomesByTyping('{table_id}', 'st_o')">ST‑OCL</button>
            <button class="group-btn" data-group="k_o" onclick="groupGenomesByTyping('{table_id}', 'k_o')">K:OCL</button>
            <button class="group-btn" data-group="st_ko" onclick="groupGenomesByTyping('{table_id}', 'st_ko')">ST‑K:OCL</button>
            <button class="group-btn" onclick="resetGenomeList('{table_id}')">Reset</button>
        </div>
        """

    # --------------------------------------------------------------------------
    # Helper: Search boxes 
    # --------------------------------------------------------------------------
    def _search_boxes(self, prefix: str, table_id: str) -> str:
        return f"""
        <input type="text" class="search-box" id="search-{prefix}" onkeyup="searchTable('{table_id}','search-{prefix}')" placeholder="🔍 Search genes or specific database...">
        <input type="text" class="search-box" id="search-{prefix}-genome" onkeyup="highlightGenome('{table_id}','search-{prefix}-genome')" placeholder="🔍 Highlight genomes containing text...">
        """

    # --------------------------------------------------------------------------
    # Filter buttons 
    # --------------------------------------------------------------------------
    def _amr_filter_buttons(self) -> str:
        filters = [
            ('blaOXA- (Carbapenemase - OXA-type)', 'blaOXA-'),
            ('blaNDM (Carbapenemase - NDM)', 'NDM'),
            ('KPC (Carbapenemase - KPC)', 'KPC'),
            ('VIM (Carbapenemase - VIM)', 'VIM'),
            ('IMP (Carbapenemase - IMP)', 'IMP'),
            ('GES (Carbapenemase - GES)', 'GES'),
            ('CTX-M (ESBL)', 'CTX-M'),
            ('SHV (ESBL)', 'SHV'),
            ('blaTEM (ESBL / broad-spectrum)', 'blaTEM'),
            ('PER (ESBL)', 'PER'),
            ('VEB (ESBL)', 'VEB'),
            ('blaADC (AmpC cephalosporinase)', 'blaADC'),
            ('mcr (Colistin resistance - MCR)', 'mcr'),
            ('pmr (Colistin resistance - PMRAB)', 'pmr'),
            ('lpx (Colistin resistance - LPS modification)', 'lpx'),
            ('arn (Colistin resistance - L-Ara4N modification)', 'arn'),
            ('clpK (Colistin tolerance - heat shock protein)', 'clpK'),
            ('aac (Aminoglycoside acetyltransferase)', 'aac'),
            ('aph (Aminoglycoside phosphotransferase)', 'aph'),
            ('ant (Aminoglycoside nucleotidyltransferase)', 'ant'),
            ('armA (16S rRNA methylase - high-level)', 'armA'),
            ('rmtE (16S rRNA methylase - high-level)', 'rmtE'),
            ('tet (Tetracycline resistance)', 'tet'),
            ('gyrA_S (Fluoroquinolone resistance - gyrA mutation)', 'gyrA'),
            ('parC_S (Fluoroquinolone resistance - parC mutation)', 'parC'),
            ('qnr (Plasmid-mediated quinolone resistance)', 'qnr'),
            ('sul (Sulfonamide resistance)', 'sul'),
            ('dfr (Trimethoprim resistance)', 'dfr'),
            ('cat (Chloramphenicol acetyltransferase)', 'cat'),
            ('floR (Florfenicol/chloramphenicol efflux)', 'floR'),
            ('erm (Macrolide-lincosamide-streptogramin B)', 'erm'),
            ('msr (Macrolide efflux)', 'msr'),
            ('mph (Macrolide phosphotransferase)', 'mph'),
            ('ade (RND efflux pump - AdelJK / AdeFGH)', 'ade'),
            ('abeM (MFS efflux pump)', 'abeM'),
            ('abeS (MFS efflux pump)', 'abeS'),
            ('amvA (MATE efflux pump)', 'amvA'),
            ('abaF (ABC transporter)', 'abaF'),
            ('abaQ (MFS efflux pump)', 'abaQ'),
            ('mex (RND efflux pump - Pseudomonas type)', 'mex'),
            ('cxpE (Conserved hypothetical)', 'cxpE'),
            ('Clear', '')
        ]
        return ''.join(f'<button class="action-btn btn-{"danger" if "Carbapenemase" in label or "Colistin" in label else "warning" if "ESBL" in label else "info"}" onclick="document.getElementById(\'search-amr\').value=\'{term}\'; searchTable(\'amr-table\',\'search-amr\')">{label}</button>' for label, term in filters)

    def _virulence_filter_buttons(self) -> str:
        filters = [
            ('ompA (Biofilm)', 'ompA'), ('csu (Chaperone-usher pili - biofilm)', 'csu'),
            ('bfmR (Biofilm regulator)', 'bfmR'), ('bfmS (Biofilm regulator)', 'bfmS'),
            ('bap (Biofilm-associated protein)', 'bap'), ('pga (PGA biofilm synthesis)', 'pga'),
            ('pta (Biofilm polysaccharide)', 'pta'), ('ptk (Biofilm polysaccharide kinase)', 'ptk'),
            ('epsA (Exopolysaccharide)', 'epsA'), ('fim (Fimbrial adherence/assembly)', 'fim'),
            ('tsaP (Trimeric autotransporter)', 'tsaP'), ('ata (Acinetobacter trimeric autotransporter)', 'ata'),
            ('pil (Type IV pilus)', 'pil'), ('gsp (Type II secretion system)', 'gsp'),
            ('tss (Type VI secretion system)', 'tss'), ('clpV (T6SS ATPase)', 'clpV'),
            ('hcp (T6SS tube protein)', 'hcp'), ('vgrG (T6SS spike protein)', 'vgrG'),
            ('tse (T6SS effector)', 'tse'), ('bau (Acinetobactin siderophore)', 'bau'),
            ('bas (Acinetobactin biosynthesis)', 'bas'), ('entE (Enterobactin component)', 'entE'),
            ('hemO (Heme oxygenase)', 'hemO'), ('barA (Iron-regulated signalling)', 'barA'),
            ('barB (Iron-regulated signalling)', 'barB'), ('lpx (Lipid A biosynthesis)', 'lpx'),
            ('lpsB (LPS biosynthesis)', 'lpsB'), ('galE (Capsule/LPS epimerase)', 'galE'),
            ('galU (Capsule biosynthesis)', 'galU'), ('plc (Phospholipase)', 'plc'),
            ('cpaA (Complement protease)', 'cpaA'), ('abaI (Quorum sensing AHL synthase)', 'abaI'),
            ('abaR (Quorum sensing receptor)', 'abaR'), ('adeF (Efflux & virulence regulator)', 'adeF'),
            ('adeG (Efflux & virulence regulator)', 'adeG'), ('adeH (Efflux & virulence regulator)', 'adeH'),
            ('tagX (Teichoic acid biosynthesis)', 'tagX'), ('pbpG (Penicillin-binding protein)', 'pbpG'),
            ('clpK (Heat shock protein, colistin tolerance)', 'clpK'), ('Clear', '')
        ]
        return ''.join(f'<button class="action-btn btn-{"success" if "Biofilm" in label else "info"}" onclick="document.getElementById(\'search-vir\').value=\'{term}\'; searchTable(\'vir-table\',\'search-vir\')">{label}</button>' for label, term in filters)

    def _bacmet_filter_buttons(self) -> str:
        filters = [
            ('qac (Biocide - quaternary ammonium)', 'qac'), ('qacEdelta1 (Biocide - truncated qacE)', 'qacEdelta1'),
            ('cep (Biocide - chlorhexidine)', 'cep'), ('form (Biocide - formaldehyde)', 'form'),
            ('blt (Bleomycin resistance)', 'blt'), ('mer (Mercury resistance)', 'mer'),
            ('ars (Arsenic resistance)', 'ars'), ('arsT (Arsenic resistance - ArsT)', 'arsT'),
            ('cop (Copper homeostasis)', 'cop'), ('sil (Silver resistance)', 'sil'),
            ('chr (Chromate resistance)', 'chr'), ('cad (Cadmium resistance)', 'cad'),
            ('znt (Zinc efflux)', 'znt'), ('czc (Cobalt‑zinc‑cadmium efflux)', 'czc'),
            ('pbr (Lead resistance)', 'pbr'), ('hmr (Heavy metal resistance)', 'hmr'),
            ('nre (Nickel resistance regulator)', 'nre'), ('pco (Copper resistance operon)', 'pco'),
            ('nik (Nickel transport)', 'nik'), ('cor (Magnesium/cobalt transport)', 'cor'),
            ('mod (Molybdate transport)', 'mod'), ('sit (Iron/manganese transport)', 'sit'),
            ('fec (Iron‑dicitrate transport)', 'fec'), ('fpv (Pyoverdine receptor - iron)', 'fpv'),
            ('soxR (Oxidative stress regulator)', 'soxR'), ('cpxR (Envelope stress regulator)', 'cpxR'),
            ('baeR (Multidrug efflux regulator)', 'baeR'), ('ydeP (Acid resistance / efflux)', 'ydeP'),
            ('ade (RND efflux pump - Ade family)', 'ade'), ('mex (RND efflux pump - Mex family)', 'mex'),
            ('emr (MFS efflux pump - Emr family)', 'emr'), ('sme (MFS efflux pump - Sme family)', 'sme'),
            ('norA (MFS efflux pump - NorA family)', 'norA'), ('vce (Multidrug efflux - Vce family)', 'vce'),
            ('fab (Fatty acid synthesis adaptation)', 'fab'), ('oxyRkp (Oxidative stress response)', 'oxyRkp'),
            ('kla/tel/kil (Plasmid‑associated)', 'kla'), ('czrA (Cadmium‑zinc regulator)', 'czrA'),
            ('cutA (Copper tolerance)', 'cutA'), ('ter (Tellurite resistance)', 'ter'),
            ('smd (Sulfonamide/multidrug efflux)', 'smd'), ('vex (Efflux pump)', 'vex'),
            ('srp (Stress response)', 'srp'), ('dps (DNA protection under starvation)', 'dps'),
            ('Clear', '')
        ]
        return ''.join(f'<button class="action-btn btn-{"warning" if "Biocide" in label else "secondary"}" onclick="document.getElementById(\'search-bac\').value=\'{term}\'; searchTable(\'bac-table\',\'search-bac\')">{label}</button>' for label, term in filters)

    def _plasmid_filter_buttons(self) -> str:
        filters = [
            ('Rep (Rolling-circle)', 'rep'), ('Inc (Incompatibility group)', 'Inc'),
            ('Colicin plasmid', 'col'), ('Broad-host-range', 'broad'),
            ('Mobility/conjugation', 'mob'), ('Acinetobacter plasmid', 'acinetobacter'),
            ('Clear', '')
        ]
        return ''.join(f'<button class="action-btn btn-{"info" if "Rep" in label else "secondary"}" onclick="document.getElementById(\'search-plasmidfinder\').value=\'{term}\'; searchTable(\'plasmidfinder-table\',\'search-plasmidfinder\')">{label}</button>' for label, term in filters)

    def _mutation_filter_buttons(self) -> str:
        filters = [
            ('gyrA (Quinolone)', 'gyrA'), ('parC (Quinolone)', 'parC'),
            ('rpoB (Rifampin)', 'rpoB'), ('23S rRNA (Tigecycline/Linezolid)', '23S'),
            ('pmrAB (Colistin)', 'pmr'), ('Clear', '')
        ]
        return ''.join(f'<button class="action-btn btn-{"danger" if "23S" in label else "warning"}" onclick="document.getElementById(\'search-mutation\').value=\'{term}\'; searchTable(\'mutation-table\',\'search-mutation\')">{label}</button>' for label, term in filters)


# =============================================================================
# MAIN REPORTER CLASS — v3.2.0
# =============================================================================
class GeniusUltimateReporter:
    """Orchestrate parsing, analysis, and report generation."""

    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "GENIUS_ACINETOBACTER_GENE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = UltimateHTMLParser()
        self.analyzer = UltimateDataAnalyzer()
        self.html_generator = UltimateHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "GENIUS Acinetobacter baumannii Ultimate Reporter",
            "version": "3.2.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }

    def find_html_files(self) -> Dict[str, List[Path]]:
        print("🔍 Searching for AcinetoScope HTML reports...")
        html_files = {
            'pasteur_mlst': [], 'oxford_mlst': [], 'kaptive': [], 'amrfinder': [],
            'abricate': defaultdict(list), 'plasmidfinder': [], 'qc': [], 'apt_plasmid': [],
            'mutation': [], 'mge': []
        }
        all_html = list(self.input_dir.glob("**/*.html"))
        print(f"  📁 Found {len(all_html)} HTML files")
        for f in all_html:
            name = f.name.lower()
            if 'mge_summary' in name:
                html_files['mge'].append(f)
            elif name == 'plasmid_summary_report.html':
                html_files['apt_plasmid'].append(f)
            elif 'pasteur' in name and 'mlst' in name:
                html_files['pasteur_mlst'].append(f)
            elif 'oxford' in name and 'mlst' in name:
                html_files['oxford_mlst'].append(f)
            elif 'kaptive' in name:
                html_files['kaptive'].append(f)
            elif 'amrfinder' in name:
                html_files['amrfinder'].append(f)
            elif 'plasmidfinder' in name:
                html_files['plasmidfinder'].append(f)
            elif 'fasta_qc' in name or 'qc_summary' in name:
                html_files['qc'].append(f)
            elif 'mutation_summary' in name:
                html_files['mutation'].append(f)
            else:
                matched = False
                for db_key in self.parser.db_name_mapping.keys():
                    if db_key in name and 'plasmidfinder' not in db_key:
                        html_files['abricate'][db_key].append(f)
                        matched = True
                        break
                if not matched:
                    if 'card' in name:
                        html_files['abricate']['acineto_card'].append(f)
                    elif 'resfinder' in name:
                        html_files['abricate']['acineto_resfinder'].append(f)
                    elif 'argannot' in name:
                        html_files['abricate']['acineto_argannot'].append(f)
                    elif 'vfdb' in name:
                        html_files['abricate']['acineto_vfdb'].append(f)
                    elif 'bacmet2' in name:
                        html_files['abricate']['acineto_bacmet2'].append(f)
                    elif 'megares' in name:
                        html_files['abricate']['acineto_megares'].append(f)
                    elif 'ecoh' in name:
                        html_files['abricate']['acineto_ecoh'].append(f)
                    elif 'ncbi' in name:
                        html_files['abricate']['acineto_ncbi'].append(f)
        print(f"  ✅ Pasteur MLST: {len(html_files['pasteur_mlst'])} | Oxford MLST: {len(html_files['oxford_mlst'])} | Kaptive: {len(html_files['kaptive'])} | AMRfinder: {len(html_files['amrfinder'])} | PlasmidFinder: {len(html_files['plasmidfinder'])} | QC: {len(html_files['qc'])} | APT plasmid: {len(html_files['apt_plasmid'])} | Mutation: {len(html_files['mutation'])} | MGE: {len(html_files['mge'])}")
        return html_files

    def integrate_all_data(self, html_files: Dict[str, List[Path]]) -> Dict[str, Any]:
        print("\n🔗 Integrating data...")
        integrated = {'metadata': self.metadata, 'samples': {}, 'patterns': {}, 'gene_centric': {},
                      'plasmid_analysis': {}, 'qc_data': {}, 'apt_plasmid': {},
                      'mutation_data': {}, 'mge_data': {}}

        if html_files['qc']:
            integrated['qc_data'] = self.parser.parse_qc_report(html_files['qc'][0])

        if html_files['apt_plasmid']:
            integrated['apt_plasmid'] = self.parser.parse_apt_plasmid_summary(html_files['apt_plasmid'][0])

        if html_files['mge']:
            integrated['mge_data'] = self.parser.parse_mge_summary_html(html_files['mge'][0])

        pasteur = self.parser.parse_mlst_report(html_files['pasteur_mlst'][0], "pasteur") if html_files['pasteur_mlst'] else {}
        oxford = self.parser.parse_mlst_report(html_files['oxford_mlst'][0], "oxford") if html_files['oxford_mlst'] else {}
        kaptive = self.parser.parse_kaptive_report(html_files['kaptive'][0]) if html_files['kaptive'] else {}

        amr_by_sample, amr_gene_freq = {}, {}
        if html_files['amrfinder']:
            amr_by_sample, amr_gene_freq = self.parser.parse_amrfinder_report(html_files['amrfinder'][0],
                                                                              len(set(pasteur.keys()) | set(oxford.keys()) | set(kaptive.keys())))

        abricate_by_sample = defaultdict(dict)
        abricate_gene_freq = {}
        for db_key, files in html_files['abricate'].items():
            if files:
                db_name = self.parser.db_name_mapping.get(db_key, db_key)
                by_sample, freq = self.parser.parse_abricate_database_report(files[0], len(set(pasteur.keys())))
                for s, genes in by_sample.items():
                    abricate_by_sample[s][db_name] = genes
                abricate_gene_freq[db_name] = freq

        plasmid_by_sample, plasmid_gene_freq = {}, {}
        if html_files['plasmidfinder']:
            plasmid_by_sample, plasmid_gene_freq = self.parser.parse_plasmidfinder_report(html_files['plasmidfinder'][0], len(set(pasteur.keys())))

        mutation_data = {}
        if html_files['mutation']:
            mutation_data = self.parser.parse_mutation_summary_html(html_files['mutation'][0])
        integrated['mutation_data'] = mutation_data

        all_samples = (set(pasteur.keys()) | set(oxford.keys()) | set(kaptive.keys()) |
                       set(amr_by_sample.keys()) | set(abricate_by_sample.keys()) |
                       set(plasmid_by_sample.keys()) | set(integrated['qc_data'].keys()) |
                       set(integrated['mge_data'].get('per_sample', {}).keys()))
        total_samples = len(all_samples)
        print(f"📊 Found {total_samples} unique samples")

        for sample in all_samples:
            virulence = []
            for db in ['vfdb', 'ecoli_vf']:
                if db in abricate_by_sample.get(sample, {}):
                    virulence.extend(abricate_by_sample[sample][db])
            integrated['samples'][sample] = {
                'pasteur_mlst': pasteur.get(sample, {'ST': 'ND', 'International_Clone': 'Unknown', 'Allele_Profile': 'ND'}),
                'oxford_mlst': oxford.get(sample, {'ST': 'ND', 'Allele_Profile': 'ND'}),
                'kaptive': kaptive.get(sample, {'K_Locus': 'ND', 'O_Locus': 'ND', 'Capsule_Type': 'ND'}),
                'amr_genes': amr_by_sample.get(sample, []),
                'virulence_genes': list(set(virulence)),
                'environmental_genes': abricate_by_sample.get(sample, {}).get('bacmet2', []),
                'plasmid_genes': plasmid_by_sample.get(sample, []),
                'abricate_databases': dict(abricate_by_sample.get(sample, {})),
            }

        integrated['gene_frequencies'] = {'amrfinder': amr_gene_freq, 'abricate': abricate_gene_freq}
        if plasmid_gene_freq:
            integrated['gene_frequencies']['plasmidfinder'] = plasmid_gene_freq

        integrated['gene_centric'] = self.analyzer.create_gene_centric_tables(integrated, total_samples)
        integrated['patterns'] = self.analyzer.create_cross_genome_patterns(integrated, total_samples)
        integrated['plasmid_analysis'] = self.analyzer.create_plasmid_analysis(integrated, total_samples)

        return integrated

    def generate_json_report(self, data: Dict[str, Any]) -> Path:
        print("\n📝 Generating JSON report...")
        output_file = self.output_dir / "genius_acinetobacter_ultimate_gene_centric_report.json"

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
            elif isinstance(obj, defaultdict):
                return convert_keys(dict(obj))
            elif isinstance(obj, Counter):
                return dict(obj)
            else:
                return obj

        serializable = convert_keys(data)
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(serializable, f, indent=2, default=str)
        print(f"    ✅ JSON saved: {output_file}")
        return output_file

    def generate_csv_reports(self, data: Dict[str, Any]):
        print("\n📊 Generating CSV reports...")
        samples_df = pd.DataFrame([{
            'Sample': s,
            'Pasteur_ST': d['pasteur_mlst']['ST'],
            'Oxford_ST': d['oxford_mlst']['ST'],
            'K_Locus': d['kaptive']['K_Locus'],
            'O_Locus': d['kaptive']['O_Locus'],
            'Virulence_Count': len(d['virulence_genes'])
        } for s, d in data['samples'].items()])
        samples_df.to_csv(self.output_dir / "acinetobacter_samples.csv", index=False)

        amr_rows = []
        for db, genes in data['gene_centric'].get('amr_databases', {}).items():
            for g in genes:
                amr_rows.append({'Gene': g['gene'], 'Database': g['database'], 'Frequency': g['frequency_display'], 'Genomes': ';'.join(g['genomes'])})
        if amr_rows:
            pd.DataFrame(amr_rows).to_csv(self.output_dir / "amr_genes.csv", index=False)

        vir_rows = []
        for db, genes in data['gene_centric'].get('virulence_databases', {}).items():
            for g in genes:
                vir_rows.append({'Gene': g['gene'], 'Database': g['database'], 'Frequency': g['frequency_display'], 'Genomes': ';'.join(g['genomes'])})
        if vir_rows:
            pd.DataFrame(vir_rows).to_csv(self.output_dir / "virulence_genes.csv", index=False)

        bac_rows = []
        for db, genes in data['gene_centric'].get('bacmet_databases', {}).items():
            for g in genes:
                bac_rows.append({'Gene': g['gene'], 'Frequency': g['frequency_display'], 'Genomes': ';'.join(g['genomes'])})
        if bac_rows:
            pd.DataFrame(bac_rows).to_csv(self.output_dir / "bacmet_genes.csv", index=False)

        plas_rows = []
        for db, genes in data.get('plasmid_analysis', {}).get('plasmid_databases', {}).items():
            for g in genes:
                plas_rows.append({'Marker': g['plasmid_marker'], 'Frequency': g['frequency_display'], 'Genomes': ';'.join(g['genomes'])})
        if plas_rows:
            pd.DataFrame(plas_rows).to_csv(self.output_dir / "plasmid_markers.csv", index=False)

        if data.get('apt_plasmid') and data['apt_plasmid'].get('rep_frequency'):
            apt_rows = []
            for item in data['apt_plasmid']['rep_frequency']:
                apt_rows.append({'Rep_Type': item['rep_type'], 'Count': item['count'], 'Genomes': ';'.join(item['genomes'])})
            pd.DataFrame(apt_rows).to_csv(self.output_dir / "apt_plasmid_rep_frequency.csv", index=False)

        if data.get('mutation_data', {}).get('mutations'):
            mut_rows = []
            for m in data['mutation_data']['mutations']:
                mut_rows.append({'Gene': m['gene'], 'Mutation': m['mutation'], 'Class': m['class'], 'Subclass': m['subclass'], 'Count': m['count'], 'Genomes': ';'.join(m['genomes'])})
            pd.DataFrame(mut_rows).to_csv(self.output_dir / "mutations.csv", index=False)

        mge_ps = data.get('mge_data', {}).get('per_sample', {})
        if mge_ps:
            pd.DataFrame([{'Sample': s, **v} for s, v in mge_ps.items()]).to_csv(self.output_dir / "mge_profile.csv", index=False)

        cooc = data['patterns'].get('gene_cooccurrence', {})
        if cooc:
            cooc_rows = []
            for g1, partners in cooc.items():
                for g2, cnt in partners.items():
                    cooc_rows.append({'Gene1': g1, 'Gene2': g2, 'Count': cnt})
            cooc_rows.sort(key=lambda x: x['Count'], reverse=True)
            pd.DataFrame(cooc_rows[:500]).to_csv(self.output_dir / "cooccurrence.csv", index=False)

        print("    ✅ CSV reports written")

    def run(self):
        print("=" * 80)
        print("🧠 GENIUS ACINETOBACTER BAUMANNII ULTIMATE REPORTER v3.2.0")
        print("   Gene-centric with dynamic grouping, mutations, co-occurrence")
        print("   + MGE tab, Compare tab, Call to Action tab")
        print("=" * 80)
        html_files = self.find_html_files()
        if not any(html_files.values()):
            print("❌ No HTML files found!")
            return False
        data = self.integrate_all_data(html_files)
        if not data:
            return False
        self.generate_json_report(data)
        self.generate_csv_reports(data)
        self.html_generator.generate_main_report(data, self.output_dir)
        n = len(data['samples'])
        mge_n = len(data['mge_data'].get('per_sample', {}))
        print("\n" + "=" * 80)
        print("✅ REPORT COMPLETE")
        print("=" * 80)
        print(f"   Samples:          {n}")
        print(f"   MGE-profiled:     {mge_n}")
        print(f"   Output dir:       {self.output_dir}")
        print(f"   Open:             {self.output_dir / 'genius_acinetobacter_ultimate_gene_centric_report.html'}")
        print("=" * 80)
        return True


def main():
    parser = argparse.ArgumentParser(description='GENIUS Acinetobacter baumannii Ultimate Reporter v3.2.0')
    parser.add_argument('-i', '--input-dir', required=True, help='Directory with AcinetoScope HTML reports')
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        sys.exit(1)
    reporter = GeniusUltimateReporter(input_dir)
    success = reporter.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()