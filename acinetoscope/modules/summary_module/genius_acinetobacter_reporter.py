#!/usr/bin/env python3
"""
GENIUS ACINETOBACTER BAUMANNII ULTIMATE REPORTER 
Advanced HTML Parser with Gene-Centric Cross-Genome Analysis for A. baumannii
Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 2.0.0 
Date: 2026-04-28
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

# =============================================================================
# HTML PARSER CLASS 
# =============================================================================
class UltimateHTMLParser:
    """Ultimate HTML parser for AcinetoScope reports"""
    
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
        sample = str(sample_id)
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
                break
        if sample.startswith('GCF_'):
            sample = 'GCA_' + sample[4:]
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
            
            # Find sample column (usually "Sample" or "Genome")
            sample_col = None
            for col in df.columns:
                col_lower = col.lower()
                if 'sample' in col_lower or 'genome' in col_lower or 'id' in col_lower or 'strain' in col_lower:
                    if '#' not in col_lower:
                        sample_col = col
                        break
            if not sample_col and len(df.columns) > 1:
                # fallback: assume first column is sample
                sample_col = df.columns[0]
            if not sample_col:
                print(f"    ⚠️ Could not find sample column in {scheme} MLST")
                return {}
            
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            results = {}
            valid_st_count = 0
            
            # Determine the ST column: look for exact "ST" or any column containing "ST"
            st_col = None
            for col in df.columns:
                if col.upper() == 'ST':
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
                # Try to get ST from the dedicated column
                if st_col is not None and pd.notna(row.get(st_col)):
                    st_val = str(row[st_col]).strip()
                    if st_val and st_val.upper() not in ['', 'NAN', 'NONE', 'ND', 'UNKNOWN', 'STUNKNOWN']:
                        # Extract the numeric part after "ST" (case-insensitive)
                        match = re.search(r'ST(\d+)', st_val, re.IGNORECASE)
                        if match:
                            st = match.group(1)
                            valid_st_count += 1
                        else:
                            # Try to extract any number
                            num_match = re.search(r'\d+', st_val)
                            if num_match:
                                st = num_match.group()
                                valid_st_count += 1
                else:
                    # No explicit ST column, search all columns for ST pattern
                    for col in df.columns:
                        if col.lower() != sample_col.lower() and pd.notna(row.get(col)):
                            cell_val = str(row[col])
                            match = re.search(r'ST(\d+)', cell_val, re.IGNORECASE)
                            if match:
                                st = match.group(1)
                                valid_st_count += 1
                                break
                
                # International Clone (only for Pasteur scheme)
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
                            ic_match = re.search(r'(IC\s*[I|II|III|IV|V|VI|VII|VIII|IX|X]+)', ic_val, re.I)
                            if ic_match:
                                ic = ic_match.group(1).replace(' ', '')
                            else:
                                ic = ic_val.replace(' ', '')
                    if ic == 'Unknown':
                        for col in df.columns:
                            if col.lower() != sample_col.lower() and pd.notna(row.get(col)):
                                cell_val = str(row[col])
                                ic_match = re.search(r'(IC\s*[I|II|III|IV|V|VI|VII|VIII|IX|X]+)', cell_val, re.I)
                                if ic_match:
                                    ic = ic_match.group(1).replace(' ', '')
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
            
            print(f"    ✓ Found {len(results)} samples for {scheme.upper()} scheme (valid STs: {valid_st_count})")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing {scheme} MLST: {e}")
            import traceback
            traceback.print_exc()
            return {}
    
    # ---------- All other parsing methods ---------------------------

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
                            o_match = re.search(r'(OC\d+)', o_val, re.I)
                            if o_match:
                                o_locus = o_match.group(1).upper()
                            else:
                                unknown_match = re.search(r'unknown\s*\(OCL(\d+)\)', o_val, re.I)
                                if unknown_match:
                                    o_locus = f"OC{unknown_match.group(1)}"
                                else:
                                    ocl_match = re.search(r'OCL(\d+)', o_val, re.I)
                                    if ocl_match:
                                        o_locus = f"OC{ocl_match.group(1)}"
                                    else:
                                        parts = o_val.split()
                                        for part in parts:
                                            if part.startswith('OC') and part[2:].isdigit():
                                                o_locus = part.upper()
                                                break
                                            elif part.startswith('O') and part[1:].isdigit():
                                                o_locus = 'OC' + part[1:]
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
    
    def parse_amrfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}
            gene_frequencies = {}
            df_freq = self.parse_html_table(str(tables[1]), 0)
            if not df_freq.empty and 'Gene' in df_freq.columns:
                for _, row in df_freq.iterrows():
                    gene = str(row['Gene']).strip()
                    if not gene:
                        continue
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
                        'database': 'amrfinder'
                    }
            genes_by_genome = {}
            df_genomes = self.parse_html_table(str(tables[0]), 0)
            if not df_genomes.empty and 'Genome' in df_genomes.columns:
                for _, row in df_genomes.iterrows():
                    sample = self.normalize_sample_id(row['Genome'])
                    if not sample:
                        continue
                    genes = []
                    if pd.notna(row.get('Genes Detected')):
                        gene_str = str(row['Genes Detected'])
                        genes = [g.strip() for g in gene_str.split() if g.strip() and not g.startswith('Showing') and not '(' in g and not ')' in g]
                    genes_by_genome[sample] = genes
            print(f"    ✓ Found {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder: {e}")
            return {}, {}
    
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
    
    def create_gene_centric_tables(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        gene_centric = {
            'amr_databases': {},
            'virulence_databases': {},
            'bacmet_databases': {},
            'plasmid_databases': {},
            'combined_gene_frequencies': [],
            'gene_categories': defaultdict(list),
            'database_stats': {}
        }
        if 'amrfinder' in integrated_data.get('gene_frequencies', {}):
            amr_data = integrated_data['gene_frequencies']['amrfinder']
            gene_list = []
            for gene, data in amr_data.items():
                category = self.categorize_gene(gene)
                gene_list.append({
                    'gene': gene,
                    'category': category,
                    'database': 'AMRfinder',
                    'count': data.get('count', 0),
                    'percentage': data.get('percentage', 0),
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
                        'gene': gene,
                        'category': category,
                        'database': db_name.upper(),
                        'count': data.get('count', 0),
                        'percentage': data.get('percentage', 0),
                        'frequency_display': data.get('frequency_display', f"{data.get('count', 0)} ({data.get('percentage', 0):.1f}%)"),
                        'genomes': data.get('genomes', []),
                        'full_name': data.get('full_name', gene)
                    })
                if gene_list:
                    gene_list.sort(key=lambda x: x['count'], reverse=True)
                    if db_name in ['vfdb', 'ecoli_vf']:
                        gene_centric['virulence_databases'][db_name] = gene_list
                    elif db_name == 'bacmet2':
                        gene_centric['bacmet_databases'][db_name] = gene_list
                    elif db_name in ['plasmidfinder', 'ecoh']:
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
            'pasteur_st_distribution': Counter(),
            'oxford_st_distribution': Counter(),
            'k_locus_distribution': Counter(),
            'o_locus_distribution': Counter(),
            'capsule_type_distribution': Counter(),
            'st_k_combinations': defaultdict(list),
            'st_o_combinations': defaultdict(list),
            'ko_combinations': defaultdict(list),
            'st_ko_combinations': defaultdict(list),
            'carbapenemase_patterns': defaultdict(list),
            'high_risk_crab': []
        }
        samples_data = integrated_data.get('samples', {})
        gene_centric = integrated_data.get('gene_centric', {})
        sample_genes = defaultdict(set)
        for db_type in ['amr_databases', 'virulence_databases', 'bacmet_databases', 'plasmid_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    for genome in gene_data['genomes']:
                        sample_genes[genome].add(gene_data['gene'])
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
                    'sample': sample,
                    'pasteur_st': pasteur_st,
                    'k_locus': k_locus,
                    'capsule_type': capsule,
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
                'plasmid_marker': gene,
                'category': data.get('plasmid_type', 'Other'),
                'database': 'PlasmidFinder',
                'count': data.get('count', 0),
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


# =============================================================================
# HTML GENERATOR CLASS 
# =============================================================================
class UltimateHTMLGenerator:
    def __init__(self, data_analyzer: UltimateDataAnalyzer):
        self.data_analyzer = data_analyzer
    
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        print("\n🎨 Generating ULTIMATE HTML report for A. baumannii...")
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        plasmid_analysis = integrated_data.get('plasmid_analysis', {})
        qc_data = integrated_data.get('qc_data', {})
        html = self._create_ultimate_html(
            metadata=metadata, samples_data=samples_data, patterns=patterns,
            gene_centric=gene_centric, plasmid_analysis=plasmid_analysis,
            qc_data=qc_data, total_samples=len(samples_data)
        )
        output_file = output_dir / "genius_acinetobacter_ultimate_report.html"
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
        total_amr = sum(len(db) for db in gene_centric.get('amr_databases', {}).values())
        total_vir = sum(len(db) for db in gene_centric.get('virulence_databases', {}).values())
        total_bacmet = sum(len(db) for db in gene_centric.get('bacmet_databases', {}).values())
        total_plasmid = sum(len(db) for db in plasmid_analysis.get('plasmid_databases', {}).values())
        high_risk = len(patterns.get('high_risk_crab', []))
        carbapenemase_count = sum(1 for db in gene_centric.get('amr_databases', {}).values() for g in db if g['category'] == 'Carbapenemases')
        env_count = total_bacmet
        
        css = """
        <style>
        :root { --summary-color:#4CAF50; --samples-color:#2196F3; --mlst-color:#FF9800; --kaptive-color:#9C27B0; --amr-color:#F44336; --virulence-color:#E91E63; --bacmet-color:#795548; --combinations-color:#009688; --patterns-color:#FF5722; --plasmid-color:#2196F3; --databases-color:#607D8B; --qc-color:#17a2b8; --aiguide-color:#3F51B5; --calltoaction-color:#2E7D32; --export-color:#3F51B5; }
        *{margin:0;padding:0;box-sizing:border-box}
        body{font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;line-height:1.6;color:#333;background:#f5f5f5;overflow-x:auto}
        .container{width:100%;max-width:100%;margin:0 auto;padding:20px;overflow-x:hidden}
        .main-header{background:linear-gradient(135deg,#00695c 0%,#004d40 100%);color:white;padding:30px;border-radius:15px;margin-bottom:30px;text-align:center}
        .main-header h1{font-size:2.8em;margin-bottom:10px;color:white}
        .metadata-bar{background:rgba(255,255,255,0.1);padding:15px;border-radius:10px;margin:20px 0;display:flex;justify-content:space-around;flex-wrap:wrap;gap:15px}
        .dashboard-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:20px;margin-bottom:30px}
        .dashboard-card{background:white;padding:25px;border-radius:12px;box-shadow:0 5px 20px rgba(0,0,0,0.1);text-align:center;cursor:pointer;border-left:5px solid;transition:all 0.3s}
        .dashboard-card:hover{transform:translateY(-10px)}
        .card-number{font-size:3em;font-weight:bold;margin:15px 0;background:linear-gradient(90deg,#00695c,#004d40);-webkit-background-clip:text;-webkit-text-fill-color:transparent}
        .tab-navigation{display:flex;gap:5px;margin-bottom:20px;flex-wrap:wrap;background:white;padding:15px;border-radius:12px;position:sticky;top:10px;z-index:100}
        .tab-button{padding:12px 25px;background:#f5f5f5;border:none;border-radius:8px;cursor:pointer;font-weight:600;display:flex;align-items:center;gap:8px;transition:0.3s}
        .tab-button.active{color:white}
        .tab-button.summary.active{background:var(--summary-color)}
        .tab-button.samples.active{background:var(--samples-color)}
        .tab-button.mlst.active{background:var(--mlst-color)}
        .tab-button.kaptive.active{background:var(--kaptive-color)}
        .tab-button.amr.active{background:var(--amr-color)}
        .tab-button.virulence.active{background:var(--virulence-color)}
        .tab-button.bacmet.active{background:var(--bacmet-color)}
        .tab-button.combinations.active{background:var(--combinations-color)}
        .tab-button.patterns.active{background:var(--patterns-color)}
        .tab-button.plasmid.active{background:var(--plasmid-color)}
        .tab-button.databases.active{background:var(--databases-color)}
        .tab-button.qc.active{background:var(--qc-color)}
        .tab-button.aiguide.active{background:var(--aiguide-color)}
        .tab-button.calltoaction.active{background:var(--calltoaction-color)}
        .tab-button.export.active{background:var(--export-color)}
        .tab-content{display:none;background:white;padding:30px;border-radius:15px;margin-bottom:30px;animation:fadeIn 0.5s}
        .tab-content.active{display:block}
        @keyframes fadeIn{from{opacity:0;transform:translateY(20px)}to{opacity:1;transform:translateY(0)}}
        .section-header{color:#2c3e50;margin-bottom:25px;padding-bottom:15px;border-bottom:3px solid;font-size:1.8em;display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap}
        .master-scrollable-container{width:100%;overflow-x:auto;border:1px solid #e0e0e0;border-radius:8px;margin:20px 0}
        .data-table{width:100%;border-collapse:collapse;font-size:0.95em;box-shadow:0 2px 10px rgba(0,0,0,0.1);min-width:600px}
        .data-table th{background:#2c3e50;color:white;padding:15px;text-align:left;white-space:nowrap;cursor:pointer}
        .data-table th:hover{background:#1a252f}
        .data-table td{padding:12px 15px;border-bottom:1px solid #e0e0e0;vertical-align:top;word-break:break-word}
        .data-table tr:hover{background:#f8f9fa}
        .data-table td:first-child, .data-table th:first-child { min-width: 200px; white-space: nowrap; }
        .data-table td:nth-child(2), .data-table th:nth-child(2) { min-width: 100px; }
        .data-table td:nth-child(3), .data-table th:nth-child(3) { min-width: 100px; }
        .search-box{width:100%;padding:12px;margin-bottom:15px;border:2px solid #e0e0e0;border-radius:8px;font-size:1em}
        .search-box:focus{outline:none;border-color:#00695c}
        .action-buttons{display:flex;gap:10px;margin:20px 0;flex-wrap:wrap}
        .action-btn{padding:10px 20px;border:none;border-radius:8px;cursor:pointer;font-weight:600;display:inline-flex;align-items:center;gap:8px;transition:0.3s}
        .action-btn:hover{transform:translateY(-2px)}
        .btn-primary{background:#00695c;color:white}
        .btn-success{background:#28a745;color:white}
        .btn-warning{background:#ffc107;color:black}
        .btn-secondary{background:#6c757d;color:white}
        .badge{display:inline-block;padding:5px 15px;border-radius:20px;font-size:0.85em;font-weight:600;margin:2px}
        .badge-critical{background:#9C27B0;color:white}
        .badge-high{background:#F44336;color:white}
        .badge-medium{background:#FF9800;color:black}
        .alert-box{padding:20px;border-radius:10px;margin:20px 0;display:flex;align-items:center;gap:20px;border-left:5px solid}
        .alert-info{background:#d1ecf1;color:#0c5460;border-left-color:#17a2b8}
        .alert-danger{background:#f8d7da;color:#721c24;border-left-color:#dc3545}
        .genome-list{display:block;max-height:200px;overflow-y:auto;padding:5px;background:#f8f9fa;border-radius:5px}
        .genome-tag{display:inline-block;background:#e0f2f1;color:#00695c;padding:3px 10px;border-radius:12px;font-size:0.85em;margin:2px}
        .genome-tag.highlight{background-color:#ffff99 !important; color:#000 !important; border:1px solid #ffc107}
        .footer{text-align:center;padding:30px;background:linear-gradient(135deg,#2c3e50,#34495e);color:white;border-radius:15px;margin-top:40px}
        .footer a{color:#ffc107;text-decoration:none}
        .footer a:hover{text-decoration:underline}
        .info-text{background:#f8f9fa;padding:15px;border-radius:8px;margin:15px 0;border-left:4px solid #17a2b8}
        .critical-cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:20px;margin:20px 0}
        .critical-card{background:#fff;padding:20px;border-radius:10px;box-shadow:0 2px 10px rgba(0,0,0,0.1);border-left:4px solid}
        .sort-icon{margin-left:5px;font-size:0.8em;opacity:0.6}
        </style>
        """
        js = """
        <script>
        function switchTab(tabName){
            document.querySelectorAll('.tab-content').forEach(t=>t.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(b=>b.classList.remove('active'));
            document.getElementById(tabName+'-tab').classList.add('active');
            event.currentTarget.classList.add('active');
            window.location.hash=tabName;
        }
        function searchTable(tableId, searchId){
            let filter=document.getElementById(searchId).value.toUpperCase();
            let rows=document.getElementById(tableId).getElementsByTagName('tr');
            for(let i=1;i<rows.length;i++){
                let found=false;
                for(let cell of rows[i].getElementsByTagName('td')){
                    if(cell.textContent.toUpperCase().indexOf(filter)>-1){found=true;break;}
                }
                rows[i].style.display=found?'':'none';
            }
        }
        function highlightGenome(tableId, genomeSearchId){
            let filter=document.getElementById(genomeSearchId).value.toUpperCase().trim();
            let table=document.getElementById(tableId);
            let allTags=table.querySelectorAll('.genome-tag');
            allTags.forEach(tag=>tag.classList.remove('highlight'));
            if(filter==='') return;
            allTags.forEach(tag=>{
                if(tag.textContent.toUpperCase().indexOf(filter)>-1){
                    tag.classList.add('highlight');
                }
            });
        }
        function exportTableToCSV(tableId, filename){
            let rows=document.getElementById(tableId).querySelectorAll('tr');
            let csv=[];
            for(let row of rows){
                let cols=row.querySelectorAll('td,th');
                let rowData=Array.from(cols).map(c=>'"'+c.innerText.replace(/"/g,'""')+'"');
                csv.push(rowData.join(','));
            }
            let blob=new Blob([csv.join('\\n')],{type:'text/csv'});
            let a=document.createElement('a');
            a.href=URL.createObjectURL(blob);
            a.download=filename;
            a.click();
            URL.revokeObjectURL(a.href);
        }
        function printSection(sectionId){
            let content=document.getElementById(sectionId);
            let win=window.open('','_blank');
            win.document.write('<html><head><title>Print</title><style>'+document.querySelector('style').innerHTML+'</style></head><body>'+content.innerHTML+'</body></html>');
            win.document.close();
            win.print();
        }
        function sortTable(tableId, colIndex, type){
            let table=document.getElementById(tableId);
            let tbody=table.tBodies[0];
            let rows=Array.from(tbody.rows);
            let isAscending=table.getAttribute('data-sort-dir')!=='asc';
            rows.sort((a,b)=>{
                let aVal=a.cells[colIndex].innerText.trim();
                let bVal=b.cells[colIndex].innerText.trim();
                if(type==='number'){
                    aVal=parseFloat(aVal.replace(/,/g,''))||0;
                    bVal=parseFloat(bVal.replace(/,/g,''))||0;
                    return isAscending?aVal-bVal:bVal-aVal;
                }else{
                    return isAscending?aVal.localeCompare(bVal):bVal.localeCompare(aVal);
                }
            });
            tbody.append(...rows);
            table.setAttribute('data-sort-dir',isAscending?'asc':'desc');
            let headers=table.querySelectorAll('th');
            headers.forEach((th,idx)=>{
                let icon=th.querySelector('.sort-icon');
                if(icon) icon.innerHTML='⇅';
            });
            let currentHeader=headers[colIndex];
            let icon=currentHeader.querySelector('.sort-icon');
            if(icon) icon.innerHTML=isAscending?'↑':'↓';
        }
        document.addEventListener('DOMContentLoaded',function(){
            let hash=window.location.hash.substring(1);
            if(hash){
                let btn=document.querySelector(`.tab-button.${hash}`);
                if(btn) btn.click();
            } else document.querySelector('.tab-button').click();
            document.querySelectorAll('.data-table').forEach(table=>{
                let headers=table.querySelectorAll('th');
                headers.forEach((th,idx)=>{
                    let type=th.getAttribute('data-sort')||'string';
                    th.style.cursor='pointer';
                    th.addEventListener('click',()=>sortTable(table.id,idx,type));
                    let icon=document.createElement('span');
                    icon.className='sort-icon';
                    icon.innerHTML='⇅';
                    th.appendChild(icon);
                });
            });
        });
        </script>
        """
        
        html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>GENIUS Acinetobacter baumannii Ultimate Report</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
{css}{js}</head>
<body><div class="container">
<div class="main-header"><h1><i class="fas fa-bacterium"></i> GENIUS Acinetobacter baumannii Ultimate Analysis Report</h1>
<p>Gene-Centric Cross-Genome Analysis with Bacmet (Biocide/Metal) & Plasmid Tracking</p>
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
<div class="dashboard-card card-plasmid" onclick="switchTab('plasmid')"><div class="card-number">{total_plasmid}</div><div class="card-label">Plasmid Markers</div><i class="fas fa-dna fa-2x" style="color:var(--plasmid-color)"></i></div>
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
<button class="tab-button plasmid" onclick="switchTab('plasmid')"><i class="fas fa-dna"></i> Plasmids</button>
<button class="tab-button patterns" onclick="switchTab('patterns')"><i class="fas fa-project-diagram"></i> Patterns</button>
<button class="tab-button databases" onclick="switchTab('databases')"><i class="fas fa-database"></i> Database Metrics</button>
<button class="tab-button aiguide" onclick="switchTab('aiguide')"><i class="fas fa-robot"></i> AI Guide</button>
<button class="tab-button calltoaction" onclick="switchTab('calltoaction')"><i class="fas fa-globe"></i> Call to Action</button>
<button class="tab-button export" onclick="switchTab('export')"><i class="fas fa-download"></i> Export</button>
</div>
<div id="summary-tab" class="tab-content active">{self._summary_section(kwargs, carbapenemase_count, env_count)}</div>
<div id="samples-tab" class="tab-content">{self._samples_section(kwargs)}</div>
<div id="qc-tab" class="tab-content">{self._qc_section(kwargs)}</div>
<div id="mlst-tab" class="tab-content">{self._mlst_section(kwargs)}</div>
<div id="kaptive-tab" class="tab-content">{self._kaptive_section(kwargs)}</div>
<div id="combinations-tab" class="tab-content">{self._combinations_section(kwargs)}</div>
<div id="amr-tab" class="tab-content">{self._amr_section(kwargs)}</div>
<div id="virulence-tab" class="tab-content">{self._virulence_section(kwargs)}</div>
<div id="bacmet-tab" class="tab-content">{self._bacmet_section(kwargs)}</div>
<div id="plasmid-tab" class="tab-content">{self._plasmid_section(kwargs)}</div>
<div id="patterns-tab" class="tab-content">{self._patterns_section(kwargs)}</div>
<div id="databases-tab" class="tab-content">{self._databases_section(kwargs)}</div>
<div id="aiguide-tab" class="tab-content">{self._aiguide_section()}</div>
<div id="calltoaction-tab" class="tab-content">{self._calltoaction_section()}</div>
<div id="export-tab" class="tab-content">{self._export_section()}</div>
<div class="footer"><h3>GENIUS Acinetobacter baumannii Ultimate Reporter v2.0.0</h3><p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p><p>If you find this tool helpful, please <a href="https://github.com/bbeckley-hub/acinetoscope" target="_blank">⭐ star us on GitHub</a> and share with your network.</p><p>Critical Genes Tracked: Carbapenemases • ESBLs • Colistin Resistance • Tigecycline Resistance • Biofilm Formation • Efflux Pumps • Environmental Co-Selection</p><p>Generated on {metadata.get('analysis_date','Unknown')}</p></div>
</div></body></html>"""
        return html
    
    # ------------------------------ All section methods ------------------------------
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
        return f"""
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>MDR/XDR A. baumannii Analysis with Environmental Co-Selection</h3><p>This comprehensive gene-centric report analyzes <strong>{total}</strong> A. baumannii genomes with focus on carbapenem resistance, colistin resistance, biofilm formation, and environmental co-selection markers. Each gene is shown with ALL genomes that contain it for easy tracking of resistance spread.</p></div></div>
        {f'<div class="alert-box alert-danger"><i class="fas fa-exclamation-triangle fa-2x"></i><div><h3>⚠️ Critical Resistance Alert</h3><p><strong>{carb_count} carbapenemase genes</strong> detected across samples. Carbapenem-resistant A. baumannii (CRAB) is a WHO Critical Priority pathogen.</p></div></div>' if carb_count>0 else ''}
        {f'<div class="alert-box alert-info"><i class="fas fa-globe-africa fa-2x"></i><div><h3>⚠️ Environmental Co-Selection Alert</h3><p><strong>{env_count} environmental resistance markers</strong> detected. These genes can co-select for antibiotic resistance in hospital environments.</p></div></div>' if env_count>0 else ''}
        <h3>Key Statistics</h3>
        <div class="master-scrollable-container"><table id="summary-table" class="data-table"><thead><tr><th data-sort="string">Metric</th><th data-sort="number">Count</th><th data-sort="string">Details</th></tr></thead><tbody>
        <tr><td>Total Samples Analyzed</th><td><strong>{total}</strong></td><td>Complete genomic analysis with all databases</th></tr>
        <tr><td>Pasteur STs</th><td><strong>{st_unique}</strong></td><td>MLST typing (Pasteur scheme)</th></tr>
        <tr><td>International Clones</th><td><strong>{ic_unique}</strong></td><td>IC classification</th></tr>
        <tr><td>Capsule (K) Loci</th><td><strong>{k_unique}</strong></td><td>Kaptive capsule typing</th></tr>
        <tr><td>Total AMR Genes</th><td><strong>{amr_total}</strong></td><td>Across all AMR databases</th></tr>
        <tr><td>Carbapenemase Genes</th><td><span class="badge {'badge-critical' if carb_count>0 else 'badge-low'}">{carb_count}</span></td><td>OXA, NDM, VIM, IMP, KPC types</th></tr>
        <tr><td>Virulence Genes</th><td><strong>{vir_total}</strong></td><td>Biofilm formation, virulence factors</th></tr>
        <tr><td>Environmental Markers</th><td><span class="badge {'badge-high' if env_count>0 else 'badge-low'}">{env_count}</span></td><td>Heavy metal, biocide, stress response</th></tr>
        <tr><td>High-Risk Combinations</th><td><span class="badge {'badge-critical' if high_risk>0 else 'badge-low'}">{high_risk}</span></td><td>Carbapenemase + last-resort resistance</th></tr>
        </tbody></table></div>
        <h3>Critical Resistance Categories Tracked</h3>
        <div class="critical-cards">
            <div class="critical-card" style="border-left-color:#F44336"><h4><i class="fas fa-skull-crossbones"></i> Carbapenemases ({carb_count})</h4><p>OXA-23, OXA-58, NDM, VIM, IMP, KPC</p></div>
            <div class="critical-card" style="border-left-color:#9C27B0"><h4><i class="fas fa-tint"></i> Colistin Resistance</h4><p>mcr genes, pmrAB mutations, LPS modifications</p></div>
            <div class="critical-card" style="border-left-color:#FF9800"><h4><i class="fas fa-pills"></i> Tigecycline Resistance</h4><p>tet(X) variants, efflux pumps (adeABC)</p></div>
            <div class="critical-card" style="border-left-color:#4CAF50"><h4><i class="fas fa-layer-group"></i> Biofilm Formation</h4><p>ompA, csuABCDE, bfmRS, pili genes</p></div>
            <div class="critical-card" style="border-left-color:#795548"><h4><i class="fas fa-globe-africa"></i> Environmental Markers ({env_count})</h4><p>Heavy metals, biocides, stress response, plasmid transfer</p></div>
        </div>
        <div class="info-text"><i class="fas fa-chart-line"></i> <strong>About this section:</strong> The summary provides an executive overview of the analysis, highlighting key resistance threats, population structure, and environmental co-selection markers. Use the dashboard cards to navigate to specific analysis tabs.</div>
        """
    
    def _samples_section(self, kwargs):
        samples = kwargs['samples_data']
        gene_centric = kwargs['gene_centric']
        sample_vir_counts = defaultdict(int)
        for db_name, genes in gene_centric.get('virulence_databases', {}).items():
            for g in genes:
                for genome in g['genomes']:
                    sample_vir_counts[genome] += 1
        html = """
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><h3>Population Structure Overview</h3><p>This table summarises the key typing results for each genome. Understanding the population structure helps identify dominant clones, track outbreaks, and link genotypes to phenotypes.</p><ul><li><strong>MLST (Sequence Type)</strong>: Gold standard for global epidemiology. Clonal complexes indicate recent common ancestry.</li><li><strong>Capsule Typing (K:O)</strong>: K (capsule) and O (lipooligosaccharide) loci are critical for virulence and immune evasion. Specific K:O combinations are associated with different clonal complexes.</li><li><strong>ND (Not Determined)</strong>: Indicates that the typing result could not be determined from the available data.</li></ul></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('samples-tab')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById('search-samples').value=''; searchTable('samples-table','search-samples')"><i class="fas fa-sync"></i> Clear Search</button><button class="action-btn btn-success" onclick="exportTableToCSV('samples-table','acinetobacter_samples.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        <input type="text" class="search-box" id="search-samples" onkeyup="searchTable('samples-table','search-samples')" placeholder="🔍 Search samples by any field...">
        <div class="master-scrollable-container"><table id="samples-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">Pasteur ST</th><th data-sort="string">Oxford ST</th><th data-sort="string">K Locus</th><th data-sort="string">O Locus</th><th data-sort="number">Virulence Count</th></tr></thead><tbody>
        """
        for sample, data in samples.items():
            pasteur = data.get('pasteur_mlst',{}).get('ST','ND')
            oxford = data.get('oxford_mlst',{}).get('ST','ND')
            k = data.get('kaptive',{}).get('K_Locus','ND')
            o = data.get('kaptive',{}).get('O_Locus','ND')
            vcount = sample_vir_counts.get(sample,0)
            pasteur_disp = f"ST{pasteur}" if pasteur != 'ND' else 'ND'
            oxford_disp = f"ST{oxford}" if oxford != 'ND' else 'ND'
            html += f'<tr><td><strong>{sample}</strong></td><td>{pasteur_disp}</td><td>{oxford_disp}</td><td>{k}</td><td>{o}</td><td>{vcount}</td></tr>'
        html += '</tbody></table></div>'
        return html
    
    def _qc_section(self, kwargs):
        qc = kwargs.get('qc_data', {})
        if not qc:
            return '<div class="alert-box alert-warning"><i class="fas fa-exclamation-circle"></i><div>FASTA QC file not found or could not parse.</div></div>'
        metrics = set()
        for d in qc.values():
            metrics.update(d.keys())
        metrics = sorted(metrics)
        html = '<div class="alert-box alert-info"><i class="fas fa-chart-line"></i><div><h3>FASTA Quality Control</h3><p>This section provides assembly quality metrics for each genome, including total sequences (contigs), total bases, GC%, N50, and more. High N50 and low contig count indicate better assembly quality.</p></div></div>'
        html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'qc-tab\')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById(\'search-qc\').value=\'\'; searchTable(\'qc-table\',\'search-qc\')"><i class="fas fa-sync"></i> Clear Search</button><button class="action-btn btn-success" onclick="exportTableToCSV(\'qc-table\',\'fasta_qc.csv\')"><i class="fas fa-download"></i> Export CSV</button></div>'
        html += '<input type="text" class="search-box" id="search-qc" onkeyup="searchTable(\'qc-table\',\'search-qc\')" placeholder="🔍 Search sample...">'
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
    
  # ------------------------------ MLST Section ---------------------------------
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
        total_pasteur = sum(pasteur.values())
        total_oxford = sum(oxford.values())
        html = f"""
        <div class="alert-box alert-info"><i class="fas fa-code-branch"></i><div><h3>MLST Analysis - Dual Scheme (Pasteur & Oxford)</h3><p>Multi‑Locus Sequence Typing (MLST) uses internal fragments of seven housekeeping genes. It is highly reproducible and defines Sequence Types (STs). Closely related STs belong to the same Clonal Complex (CC).</p><p><strong>{len(pasteur)} unique Pasteur STs</strong> and <strong>{len(oxford)} unique Oxford STs</strong> identified. International Clone (IC) classification helps track global lineages.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('mlst-tab')"><i class="fas fa-print"></i> Print Section</button></div>
        <h3>Pasteur MLST Scheme</h3><div class="master-scrollable-container"><table id="mlst-pasteur-table" class="data-table"><thead><tr><th data-sort="string">ST</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        for st, cnt in sorted(pasteur.items(), key=lambda x: x[1] if isinstance(x[1], int) else x[1].get('count',0), reverse=True):
            if st=='ND': continue
            if isinstance(cnt, dict):
                count = cnt.get('count',0)
                pct = (count/total_pasteur)*100 if total_pasteur else 0
                freq = f"{count} ({pct:.1f}%)"
            else:
                count = cnt
                pct = (count/total_pasteur)*100
                freq = f"{count} ({pct:.1f}%)"
            sample_list = ', '.join([s for s,d in samples.items() if d.get('pasteur_mlst',{}).get('ST')==str(st)])
            html += f'<tr><td><strong>ST{st}</strong></td><td>{freq}</td><td>{sample_list}</td></tr>'
        html += '</tbody></table></div>'
        html += '<h3>Oxford MLST Scheme</h3><div class="master-scrollable-container"><table id="mlst-oxford-table" class="data-table"><thead><tr><th data-sort="string">ST</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>'
        for st, cnt in sorted(oxford.items(), key=lambda x: x[1] if isinstance(x[1], int) else x[1].get('count',0), reverse=True):
            if st=='ND': continue
            if isinstance(cnt, dict):
                count = cnt.get('count',0)
                pct = (count/total_oxford)*100 if total_oxford else 0
                freq = f"{count} ({pct:.1f}%)"
            else:
                count = cnt
                pct = (count/total_oxford)*100
                freq = f"{count} ({pct:.1f}%)"
            sample_list = ', '.join([s for s,d in samples.items() if d.get('oxford_mlst',{}).get('ST')==str(st)])
            html += f'<tr><td><strong>ST{st}</strong></td><td>{freq}</td><td>{sample_list}</td></tr>'
        html += '</tbody></table></div>'
        html += '<h3>International Clone Distribution</h3><div class="master-scrollable-container"><table id="ic-table" class="data-table"><thead><tr><th data-sort="string">International Clone</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>'
        for ic, cnt in ic_dist.most_common():
            sample_list = ', '.join([s for s,d in samples.items() if d.get('pasteur_mlst',{}).get('International_Clone')==ic])
            html += f'<tr><td><strong>{ic}</strong></td><td>{cnt}</td><td>{sample_list}</td></tr>'
        html += '</tbody></table></div>'
        return html
    
    # ------------------------------ Capsule Tab  ---------------------------------
    def _kaptive_section(self, kwargs):
        patterns = kwargs['patterns']
        samples = kwargs['samples_data']
        k_dist = patterns.get('k_locus_distribution', {})
        o_dist = patterns.get('o_locus_distribution', {})
        html = """
        <div class="alert-box alert-info"><i class="fas fa-shield-alt"></i><div><h3>Capsule Typing (Kaptive) Analysis</h3><p>Kaptive identifies capsule (K) and lipooligosaccharide (O) loci in A. baumannii. The capsule is a major virulence determinant, protecting against phagocytosis and complement killing. Specific K types are associated with different clonal complexes and clinical outcomes. O locus (OCL) influences immune evasion and serum resistance.</p><p><strong>Clinical significance:</strong> Certain K types (e.g., K1, K2, K3) are linked to hypervirulent strains. Monitoring capsule distribution helps track epidemic lineages.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('kaptive-tab')"><i class="fas fa-print"></i> Print Section</button></div>
        """
        # K Locus
        html += '<h3>K Locus Distribution</h3><div class="master-scrollable-container"><table id="k-locus-table" class="data-table"><thead><tr><th data-sort="string">K Locus</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>'
        totalk = sum(k_dist.values()) if k_dist else 0
        for k, cnt in sorted(k_dist.items(), key=lambda x: x[1] if isinstance(x[1], int) else x[1].get('count',0), reverse=True):
            if k=='ND': continue
            if isinstance(cnt, dict):
                count = cnt.get('count',0)
                pct = (count/totalk)*100 if totalk else 0
                freq = f"{count} ({pct:.1f}%)"
            else:
                count = cnt
                pct = (count/totalk)*100
                freq = f"{count} ({pct:.1f}%)"
            sample_list = ', '.join([s for s,d in samples.items() if d.get('kaptive',{}).get('K_Locus')==k])
            html += f'<tr><td><strong>{k}</strong></td><td>{freq}</td><td>{sample_list}</td></tr>'
        html += '</tbody></table></div>'
        # O Locus
        html += '<h3>O Locus Distribution</h3><div class="master-scrollable-container"><table id="o-locus-table" class="data-table"><thead><tr><th data-sort="string">O Locus</th><th data-sort="number">Frequency</th><th data-sort="string">Samples</th></tr></thead><tbody>'
        totalo = sum(o_dist.values()) if o_dist else 0
        for o, cnt in sorted(o_dist.items(), key=lambda x: x[1] if isinstance(x[1], int) else x[1].get('count',0), reverse=True):
            if o=='ND': continue
            if isinstance(cnt, dict):
                count = cnt.get('count',0)
                pct = (count/totalo)*100 if totalo else 0
                freq = f"{count} ({pct:.1f}%)"
            else:
                count = cnt
                pct = (count/totalo)*100
                freq = f"{count} ({pct:.1f}%)"
            sample_list = ', '.join([s for s,d in samples.items() if d.get('kaptive',{}).get('O_Locus')==o])
            html += f'<tr><td><strong>{o}</strong></td><td>{freq}</td><td>{sample_list}</td></tr>'
        html += '</tbody></table></div>'
        return html
    
    # ------------------------------ Combination Tables ---------------------------------
    def _combinations_section(self, kwargs):
        patterns = kwargs['patterns']
        combos = [
            ('st_o_combinations', 'ST - O Locus'),
            ('st_k_combinations', 'ST - K Locus'),
            ('ko_combinations', 'K:O (Capsule Type)'),
            ('st_ko_combinations', 'ST - K:O')
        ]
        html = '<div class="alert-box alert-info"><i class="fas fa-link"></i><div><h3>Combination Tables</h3><p>Associations between Sequence Types, capsule loci, and serotypes. These combinations help identify dominant clones and track epidemic strains.</p></div></div>'
        html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'combinations-tab\')"><i class="fas fa-print"></i> Print Section</button></div>'
        for key, title in combos:
            data = patterns.get(key, {})
            if not data:
                continue
            html += f'<h3>{title}</h3><input type="text" class="search-box" id="search-{key}" onkeyup="searchTable(\'{key}-table\',\'search-{key}\')" placeholder="🔍 Search {title}...">'
            html += f'<div class="master-scrollable-container"><table id="{key}-table" class="data-table"><thead><tr><th data-sort="string">{title}</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>'
            for combo, samples in sorted(data.items(), key=lambda x: len(x[1]), reverse=True):
                html += f'<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td>{", ".join(samples)}</td></tr>'
            html += '</tbody></table></div>'
        return html
    
    # ------------------------------ AMR Section ---------------------------------
    def _amr_section(self, kwargs):
        gene_centric = kwargs['gene_centric']
        all_genes = []
        for db in gene_centric.get('amr_databases', {}).values():
            all_genes.extend(db)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        
        # Filters: each tuple = (display_label, search_term)
        filters = [
            # Carbapenemases (Critical)
            ('blaOXA- (Carbapenemase - OXA-type)', 'blaOXA-'),
            ('blaNDM (Carbapenemase - NDM)', 'NDM'),
            ('KPC (Carbapenemase - KPC)', 'KPC'),
            ('VIM (Carbapenemase - VIM)', 'VIM'),
            ('IMP (Carbapenemase - IMP)', 'IMP'),
            ('GES (Carbapenemase - GES)', 'GES'),

            # ESBLs & other β-lactamases
            ('CTX-M (ESBL)', 'CTX-M'),
            ('SHV (ESBL)', 'SHV'),
            ('blaTEM (ESBL / broad-spectrum)', 'blaTEM'),
            ('PER (ESBL)', 'PER'),
            ('VEB (ESBL)', 'VEB'),
            ('blaADC (AmpC cephalosporinase)', 'blaADC'),

            # Colistin resistance (last-line)
            ('mcr (Colistin resistance - MCR)', 'mcr'),
            ('pmr (Colistin resistance - PMRAB)', 'pmr'),
            ('lpx (Colistin resistance - LPS modification)', 'lpx'),
            ('arn (Colistin resistance - L-Ara4N modification)', 'arn'),
            ('clpK (Colistin tolerance - heat shock protein)', 'clpK'),

            # Aminoglycoside resistance
            ('aac (Aminoglycoside acetyltransferase)', 'aac'),
            ('aph (Aminoglycoside phosphotransferase)', 'aph'),
            ('ant (Aminoglycoside nucleotidyltransferase)', 'ant'),
            ('armA (16S rRNA methylase - high-level)', 'armA'),
            ('rmtE (16S rRNA methylase - high-level)', 'rmtE'),

            # Tetracycline resistance
            ('tet (Tetracycline resistance)', 'tet'),

            # Fluoroquinolone resistance (chromosomal & plasmid)
            ('gyrA_S (Fluoroquinolone resistance - gyrA mutation)', 'gyrA'),
            ('parC_S (Fluoroquinolone resistance - parC mutation)', 'parC'),
            ('parC_E (Fluoroquinolone resistance - parC mutation)', 'parC'),
            ('qnr (Plasmid-mediated quinolone resistance)', 'qnr'),

            # Sulfonamide & trimethoprim
            ('sul (Sulfonamide resistance)', 'sul'),
            ('dfr (Trimethoprim resistance)', 'dfr'),

            # Phenicol resistance
            ('cat (Chloramphenicol acetyltransferase)', 'cat'),
            ('floR (Florfenicol/chloramphenicol efflux)', 'floR'),

            # Macrolide & lincosamide
            ('erm (Macrolide-lincosamide-streptogramin B)', 'erm'),
            ('msr (Macrolide efflux)', 'msr'),
            ('mph (Macrolide phosphotransferase)', 'mph'),

            # Efflux pumps
            ('ade (RND efflux pump - AdelJK / AdeFGH)', 'ade'),
            ('abeM (MFS efflux pump)', 'abeM'),
            ('abeS (MFS efflux pump)', 'abeS'),
            ('amvA (MATE efflux pump)', 'amvA'),
            ('abaF (ABC transporter)', 'abaF'),
            ('abaQ (MFS efflux pump)', 'abaQ'),
            ('mex (RND efflux pump - Pseudomonas type)', 'mex'),

            # Other / housekeeping
            ('cxpE (Conserved hypothetical)', 'cxpE')
        ]
        
        html = """
        <div class="alert-box alert-info"><i class="fas fa-biohazard"></i><div><h3>AMR Genes (Carbapenemases, ESBLs, Colistin, Tigecycline)</h3><p>Each gene with frequency (count & %). Use filters below to focus on key resistance classes. Use the genome search to highlight isolates carrying specific genes. Carbapenemases are highlighted as critical.</p><p><strong>Focused on A. baumannii:</strong> OXA-type carbapenemases (OXA-23, OXA-58), NDM, VIM, IMP, KPC, and colistin resistance mechanisms (mcr, pmrAB).</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('amr-tab')"><i class="fas fa-print"></i> Print Section</button></div>
        <input type="text" class="search-box" id="search-amr-gene" onkeyup="searchTable('amr-table','search-amr-gene')" placeholder="🔍 Search AMR genes...">
        <input type="text" class="search-box" id="search-amr-genome" onkeyup="highlightGenome('amr-table','search-amr-genome')" placeholder="🔍 Highlight genome tags matching text">
        <div class="action-buttons">
        """
        for display, search_term in filters:
            # Use JavaScript to set the search box value to the search_term 
            html += f'<button class="action-btn btn-warning" onclick="document.getElementById(\'search-amr-gene\').value=\'{search_term}\'; searchTable(\'amr-table\',\'search-amr-gene\')">{display}</button>'
        html += '<button class="action-btn btn-primary" onclick="document.getElementById(\'search-amr-gene\').value=\'\'; searchTable(\'amr-table\',\'search-amr-gene\')">Clear</button>'
        html += '</div><div class="master-scrollable-container"><table id="amr-table" class="data-table"><thead><tr><th data-sort="string">Gene</th><th data-sort="string">Database</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead><tbody>'
        for g in all_genes:
            genome_tags = ''.join(f'<span class="genome-tag">{gen}</span>' for gen in g['genomes'])
            gene_display = f"<strong>{g['gene']}</strong>" + (' 🔥' if g['category']=='Carbapenemases' else '')
            html += f'<tr><td>{gene_display}</td><td>{g["database"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html
    
    # ------------------------------ Virulence Section ---------------------------------
    def _virulence_section(self, kwargs):
        gene_centric = kwargs['gene_centric']
        all_genes = []
        for db in gene_centric.get('virulence_databases', {}).values():
            all_genes.extend(db)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        filters = [
            # Biofilm & Adhesion
            'ompA (Biofilm)',
            'csu (Chaperone-usher pili - biofilm)',
            'bfmR (Biofilm regulator)',
            'bfmS (Biofilm regulator)',
            'bap (Biofilm-associated protein)',
            'pga (PGA biofilm synthesis)',         
            'pta (Biofilm polysaccharide)',
            'ptk (Biofilm polysaccharide kinase)',
            'epsA (Exopolysaccharide)',
            'fim (Fimbrial adherence/assembly)',    
            'tsaP (Trimeric autotransporter)',
            'ata (Acinetobacter trimeric autotransporter)',

            # Type IV Pili (T4P)
            'pil (Type IV pilus)',                 

            # Type II Secretion System (T2SS)
            'gsp (Type II secretion system (T2SS))', 

            # Type VI Secretion System (T6SS)
            'tss (Type VI secretion system (T6SS))',
            'clpV (T6SS ATPase)',
            'hcp (T6SS tube protein)',
            'vgrG (T6SS spike protein)',
            'tse (T6SS effector)',                 

            # Iron Acquisition & Siderophores
            'bau (Acinetobactin siderophore)',
            'bas (Acinetobactin biosynthesis)',
            'entE (Enterobactin component)',
            'hemO (Heme oxygenase)',
            'barA (Iron-regulated signalling)',
            'barB (Iron-regulated signalling)',

            # Lipopolysaccharide (LPS) & Capsule
            'lpx (Lipid A biosynthesis)',           
            'lpsB (LPS biosynthesis)',
            'galE (Capsule/LPS epimerase)',
            'galU (Capsule biosynthesis)',

            # Phospholipases & Cytotoxins
            'plc (Phospholipase)',                  
            'cpaA (Complement protease)',

            # Quorum Sensing & Regulation
            'abaI (Quorum sensing AHL synthase)',
            'abaR (Quorum sensing receptor)',
            'adeF (Efflux & virulence regulator)',
            'adeG (Efflux & virulence regulator)',
            'adeH (Efflux & virulence regulator)',

            # Cell Wall & Miscellaneous
            'tagX (Teichoic acid biosynthesis)',
            'pbpG (Penicillin-binding protein)',
            'clpK (Heat shock protein, colistin tolerance)'
        ]
        html = """
        <div class="alert-box alert-info"><i class="fas fa-virus"></i><div><h3>Virulence Factors (Biofilm, Adhesion, Toxins)</h3><p>Key A. baumannii virulence genes: biofilm formation (ompA, csu, bfmRS), pili, and iron uptake. Biofilm formation is critical for persistence on medical devices and hospital outbreaks. Use filters to explore specific mechanisms.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('virulence-tab')"><i class="fas fa-print"></i> Print Section</button></div>
        <input type="text" class="search-box" id="search-vir-gene" onkeyup="searchTable('vir-table','search-vir-gene')" placeholder="🔍 Search virulence genes...">
        <input type="text" class="search-box" id="search-vir-genome" onkeyup="highlightGenome('vir-table','search-vir-genome')" placeholder="🔍 Highlight genome tags matching text">
        <div class="action-buttons">
        """
        for f in filters:
            html += f'<button class="action-btn btn-success" onclick="document.getElementById(\'search-vir-gene\').value=\'{f.split(" (")[0]}\'; searchTable(\'vir-table\',\'search-vir-gene\')">{f}</button>'
        html += '<button class="action-btn btn-primary" onclick="document.getElementById(\'search-vir-gene\').value=\'\'; searchTable(\'vir-table\',\'search-vir-gene\')">Clear</button>'
        html += '</div><div class="master-scrollable-container"><table id="vir-table" class="data-table"><thead><tr><th data-sort="string">Gene</th><th data-sort="string">Database</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead><tbody>'
        for g in all_genes:
            genome_tags = ''.join(f'<span class="genome-tag">{gen}</span>' for gen in g['genomes'])
            html += f'<tr><td><strong>{g["gene"]}</strong></td><td>{g["database"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html
    
    # ------------------------------ Bacmet Section (rich filters) ---------------------------------
    def _bacmet_section(self, kwargs):
        gene_centric = kwargs['gene_centric']
        all_genes = []
        for db in gene_centric.get('bacmet_databases', {}).values():
            all_genes.extend(db)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        
        # Filters: (display_label, search_term)
        filters = [
            # Biocides & disinfectants
            ('qac (Biocide - quaternary ammonium)', 'qac'),
            ('qacEdelta1 (Biocide - truncated qacE)', 'qacEdelta1'),
            ('cep (Biocide - chlorhexidine)', 'cep'),
            ('form (Biocide - formaldehyde)', 'form'),
            ('blt (Bleomycin resistance)', 'blt'),
            
            # Heavy metals
            ('mer (Mercury resistance)', 'mer'),
            ('ars (Arsenic resistance)', 'ars'),
            ('arsT (Arsenic resistance - ArsT)', 'arsT'),
            ('cop (Copper homeostasis)', 'cop'),
            ('sil (Silver resistance)', 'sil'),
            ('chr (Chromate resistance)', 'chr'),
            ('cad (Cadmium resistance)', 'cad'),
            ('znt (Zinc efflux)', 'znt'),
            ('czc (Cobalt‑zinc‑cadmium efflux)', 'czc'),
            ('pbr (Lead resistance)', 'pbr'),
            ('hmr (Heavy metal resistance)', 'hmr'),
            ('nre (Nickel resistance regulator)', 'nre'),
            ('pco (Copper resistance operon)', 'pco'),
            
            # Transport & stress response
            ('nik (Nickel transport)', 'nik'),
            ('cor (Magnesium/cobalt transport)', 'cor'),
            ('mod (Molybdate transport)', 'mod'),
            ('sit (Iron/manganese transport)', 'sit'),
            ('fec (Iron‑dicitrate transport)', 'fec'),
            ('fpv (Pyoverdine receptor - iron)', 'fpv'),
            ('soxR (Oxidative stress regulator)', 'soxR'),
            ('cpxR (Envelope stress regulator)', 'cpxR'),
            ('baeR (Multidrug efflux regulator)', 'baeR'),
            ('ydeP (Acid resistance / efflux)', 'ydeP'),
            
            # Multidrug efflux pumps
            ('ade (RND efflux pump - Ade family)', 'ade'),
            ('mex (RND efflux pump - Mex family)', 'mex'),
            ('emr (MFS efflux pump - Emr family)', 'emr'),
            ('sme (MFS efflux pump - Sme family)', 'sme'),
            ('norA (MFS efflux pump - NorA family)', 'norA'),
            ('vce (Multidrug efflux - Vce family)', 'vce'),
            
            # Other / general
            ('fab (Fatty acid synthesis adaptation)', 'fab'),
            ('oxyRkp (Oxidative stress response)', 'oxyRkp'),
            ('kla/tel/kil (Plasmid‑associated)', 'kla'),
            ('czrA (Cadmium‑zinc regulator)', 'czrA'),
            ('cutA (Copper tolerance)', 'cutA'),
            ('ter (Tellurite resistance)', 'ter'),
            ('smd (Sulfonamide/multidrug efflux)', 'smd'),
            ('vex (Efflux pump)', 'vex'),
            ('srp (Stress response)', 'srp'),
            ('dps (DNA protection under starvation)', 'dps')
        ]
        
        html = """
        <div class="alert-box alert-info"><i class="fas fa-flask"></i><div><h3>BACMET2 Database: Biocide & Heavy Metal Resistance</h3><p>These genes confer resistance to disinfectants (quaternary ammonium, chlorhexidine) and heavy metals (mercury, copper, arsenic, etc.), which can co‑select for antibiotic resistance in hospital environments. Tracking these markers helps understand persistence and resistance spread.</p></div></div>
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection('bacmet-tab')"><i class="fas fa-print"></i> Print Section</button></div>
        <input type="text" class="search-box" id="search-bac-gene" onkeyup="searchTable('bac-table','search-bac-gene')" placeholder="🔍 Search BACMET genes...">
        <input type="text" class="search-box" id="search-bac-genome" onkeyup="highlightGenome('bac-table','search-bac-genome')" placeholder="🔍 Highlight genome tags matching text">
        <div class="action-buttons">
        """
        for display, search_term in filters:
            html += f'<button class="action-btn btn-warning" onclick="document.getElementById(\'search-bac-gene\').value=\'{search_term}\'; searchTable(\'bac-table\',\'search-bac-gene\')">{display}</button>'
        html += '<button class="action-btn btn-primary" onclick="document.getElementById(\'search-bac-gene\').value=\'\'; searchTable(\'bac-table\',\'search-bac-gene\')">Clear</button>'
        html += '</div><div class="master-scrollable-container"><table id="bac-table" class="data-table"><thead><tr><th data-sort="string">Gene</th><th data-sort="string">Database</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead><tbody>'
        for g in all_genes:
            genome_tags = ''.join(f'<span class="genome-tag">{gen}</span>' for gen in g['genomes'])
            html += f'<tr><td><strong>{g["gene"]}</strong></td><td>{g["database"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html
    
    # ------------------------------ Plasmid Section  ---------------------------------
    def _plasmid_section(self, kwargs):
        plasmid = kwargs.get('plasmid_analysis', {})
        genes = plasmid.get('plasmid_frequencies', [])
        html = '<div class="alert-box alert-info"><i class="fas fa-dna"></i><div><h3>Plasmid Analysis</h3><p>Plasmid replicons and transfer genes track horizontal gene transfer of resistance. Plasmids are major vehicles for carbapenemase and colistin resistance genes in A. baumannii.</p></div></div>'
        if not genes:
            return html + '<p>No PlasmidFinder data available.</p>'
        html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'plasmid-tab\')"><i class="fas fa-print"></i> Print Section</button></div>'
        html += '<input type="text" class="search-box" id="search-plasmid-gene" onkeyup="searchTable(\'plasmid-table\',\'search-plasmid-gene\')" placeholder="🔍 Search plasmids...">'
        html += '<input type="text" class="search-box" id="search-plasmid-genome" onkeyup="highlightGenome(\'plasmid-table\',\'search-plasmid-genome\')" placeholder="🔍 Highlight genome tags matching text">'
        html += '<div class="master-scrollable-container"><table id="plasmid-table" class="data-table"><thead><tr><th data-sort="string">Plasmid Marker</th><th data-sort="string">Category</th><th data-sort="number">Frequency</th><th class="col-genomes" data-sort="string">Genomes</th></tr></thead><tbody>'
        for g in genes:
            genome_tags = ''.join(f'<span class="genome-tag">{gen}</span>' for gen in g['genomes'])
            html += f'<tr><td><strong>{g["plasmid_marker"]}</strong></td><td>{g["category"]}</td><td>{g["frequency_display"]}</td><td class="col-genomes"><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html
    
    # ------------------------------ Patterns Section ---------------------------------
    def _patterns_section(self, kwargs):
        patterns = kwargs['patterns']
        crab = patterns.get('high_risk_crab', [])
        html = '<div class="alert-box alert-info"><i class="fas fa-project-diagram"></i><div><h3>Pattern Discovery – CRAB (Carbapenem‑resistant A. baumannii)</h3><p>Carbapenem‑resistant A. baumannii (CRAB) isolates with last‑resistance genes (colistin/tigecycline). These represent high‑risk, difficult‑to‑treat infections.</p></div></div>'
        if crab:
            html += '<div class="action-buttons"><button class="action-btn btn-primary" onclick="printSection(\'patterns-tab\')"><i class="fas fa-print"></i> Print Section</button><button class="action-btn btn-secondary" onclick="document.getElementById(\'search-patterns\').value=\'\'; searchTable(\'crab-table\',\'search-patterns\')"><i class="fas fa-sync"></i> Clear Search</button></div>'
            html += '<input type="text" class="search-box" id="search-patterns" onkeyup="searchTable(\'crab-table\',\'search-patterns\')" placeholder="🔍 Search high‑risk CRAB samples...">'
            html += '<div class="master-scrollable-container"><table id="crab-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">ST</th><th data-sort="string">Capsule (K:O)</th><th data-sort="string">Carbapenemases</th><th data-sort="string">Colistin Resistance</th><th data-sort="string">Tigecycline Resistance</th></tr></thead><tbody>'
            for c in crab:
                html += f'<tr><td><strong>{c["sample"]}</strong></td><td>{c["pasteur_st"]}</td><td>{c["capsule_type"]}</td><td>{",".join(c["carbapenemases"])}</td><td>{",".join(c["colistin_resistance"])}</td><td>{",".join(c["tigecycline_resistance"])}</td></tr>'
            html += '</tbody></table></div>'
        else:
            html += '<p>No high‑risk CRAB isolates detected.</p>'
        return html
    
    # ------------------------------ Database Metrics ---------------------------------
    def _databases_section(self, kwargs):
        stats = kwargs['gene_centric'].get('database_stats', {})
        html = '<div class="alert-box alert-info"><i class="fas fa-database"></i><div><h3>Database Metrics</h3><p>Number of unique genes and total occurrences per database.</p></div></div>'
        html += '<div class="master-scrollable-container"><table class="data-table"><thead><tr><th data-sort="string">Database</th><th data-sort="number">Unique Genes</th><th data-sort="number">Total Occurrences</th><th data-sort="number">Critical Genes</th></tr></thead><tbody>'
        for db, d in stats.items():
            html += f'<tr><td><strong>{db.upper()}</strong></td><td>{d["total_genes"]}</td><td>{d["total_occurrences"]}</td><td>{d["critical_genes"]}</td></tr>'
        html += '</tbody></table></div>'
        return html
    
    # ------------------------------ AI Guide ---------------------------------
    def _aiguide_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-robot fa-2x"></i><div><h3>AI Assistant Guide</h3><p>Upload this HTML or the JSON file to ChatGPT, Claude, or Gemini to ask questions about A. baumannii resistance.</p></div></div>
        <div class="info-text"><h4>Example questions:</h4><ul><li>Which samples carry both carbapenemase (OXA-23) and colistin resistance (mcr) genes?</li><li>List all isolates with the K3:OC1 capsule type.</li><li>What are the most common biocide resistance genes (qac) in this collection?</li><li>Show me ST-K:O combinations associated with high-risk CRAB.</li><li>Which plasmids are most prevalent and which samples carry them?</li></ul></div>
        """
    
    # ------------------------------ Call to Action ---------------------------------
    def _calltoaction_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-globe fa-2x"></i><div><h3>The Global Burden of AMR and Our Call to Action</h3><p>Antimicrobial resistance (AMR) kills an estimated 1.27 million people annually. Acinetobacter baumannii is a WHO critical priority pathogen, especially carbapenem-resistant strains (CRAB). This tool aims to empower researchers and clinicians to track resistance locally and globally.</p></div></div>
        <div style="text-align:center; margin:40px 0;">
            <i class="fas fa-star" style="font-size:3em; color:#ffc107;"></i>
            <h3>We invite you to contribute!</h3>
            <p>If you find this report useful, please <strong>star our GitHub repository</strong> and share it with your network.<br>
            <a href="https://github.com/bbeckley-hub/acinetoscope" target="_blank" style="font-size:1.2em;">https://github.com/bbeckley-hub/acinetoscope</a></p>
            <p>For contributions, bug reports, or feature requests, please open an issue or contact <strong>brownbeckley94@gmail.com</strong>.</p>
            <p><i class="fas fa-chalkboard-user"></i> We welcome collaborations to adapt this tool for other pathogens and to improve AMR surveillance.</p>
        </div>
        """
    
    # ------------------------------ Export Section ---------------------------------
    def _export_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-download"></i><div>Export data tables as CSV or download complete JSON.</div></div>
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','acinetobacter_samples.csv')">Sample Overview CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('amr-table','amr_genes.csv')">AMR Genes CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('vir-table','virulence_genes.csv')">Virulence Genes CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('bac-table','bacmet_genes.csv')">Bacmet CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('plasmid-table','plasmid_genes.csv')">Plasmid CSV</button>
            <button class="action-btn btn-success" onclick="location.href='genius_acinetobacter_ultimate_report.json'">Download JSON</button>
        </div>
        """


# =============================================================================
# MAIN REPORTER CLASS 
# =============================================================================
class GeniusUltimateReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "GENIUS_ACINETOBACTER_ULTIMATE_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = UltimateHTMLParser()
        self.analyzer = UltimateDataAnalyzer()
        self.html_generator = UltimateHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "GENIUS Acinetobacter baumannii Ultimate Reporter",
            "version": "2.0.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }
    
    def find_html_files(self) -> Dict[str, List[Path]]:
        print("🔍 Searching for AcinetoScope HTML reports...")
        html_files = {
            'pasteur_mlst': [], 'oxford_mlst': [], 'kaptive': [], 'amrfinder': [],
            'abricate': defaultdict(list), 'plasmidfinder': [], 'qc': []
        }
        all_html = list(self.input_dir.glob("**/*.html"))
        print(f"  📁 Found {len(all_html)} HTML files")
        for f in all_html:
            name = f.name.lower()
            if 'pasteur' in name and 'mlst' in name:
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
        print(f"  ✅ Pasteur MLST: {len(html_files['pasteur_mlst'])} | Oxford MLST: {len(html_files['oxford_mlst'])} | Kaptive: {len(html_files['kaptive'])} | AMRfinder: {len(html_files['amrfinder'])} | PlasmidFinder: {len(html_files['plasmidfinder'])} | QC: {len(html_files['qc'])}")
        return html_files
    
    def integrate_all_data(self, html_files: Dict[str, List[Path]]) -> Dict[str, Any]:
        print("\n🔗 Integrating data...")
        integrated = {'metadata': self.metadata, 'samples': {}, 'patterns': {}, 'gene_centric': {}, 'plasmid_analysis': {}, 'qc_data': {}}
        if html_files['qc']:
            integrated['qc_data'] = self.parser.parse_qc_report(html_files['qc'][0])
        pasteur = self.parser.parse_mlst_report(html_files['pasteur_mlst'][0], "pasteur") if html_files['pasteur_mlst'] else {}
        oxford = self.parser.parse_mlst_report(html_files['oxford_mlst'][0], "oxford") if html_files['oxford_mlst'] else {}
        kaptive = self.parser.parse_kaptive_report(html_files['kaptive'][0]) if html_files['kaptive'] else {}
        amr_by_sample, amr_gene_freq = {}, {}
        if html_files['amrfinder']:
            amr_by_sample, amr_gene_freq = self.parser.parse_amrfinder_report(html_files['amrfinder'][0], len(set(pasteur.keys())|set(oxford.keys())|set(kaptive.keys())))
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
        all_samples = set(pasteur.keys()) | set(oxford.keys()) | set(kaptive.keys()) | set(amr_by_sample.keys()) | set(abricate_by_sample.keys()) | set(plasmid_by_sample.keys()) | set(integrated['qc_data'].keys())
        total_samples = len(all_samples)
        print(f"📊 Found {total_samples} unique samples")
        for sample in all_samples:
            virulence = []
            for db in ['vfdb', 'ecoli_vf']:
                if db in abricate_by_sample.get(sample, {}):
                    virulence.extend(abricate_by_sample[sample][db])
            integrated['samples'][sample] = {
                'pasteur_mlst': pasteur.get(sample, {'ST':'ND','International_Clone':'Unknown','Allele_Profile':'ND'}),
                'oxford_mlst': oxford.get(sample, {'ST':'ND','Allele_Profile':'ND'}),
                'kaptive': kaptive.get(sample, {'K_Locus':'ND','O_Locus':'ND','Capsule_Type':'ND'}),
                'amr_genes': amr_by_sample.get(sample, []),
                'virulence_genes': list(set(virulence)),
                'environmental_genes': abricate_by_sample.get(sample, {}).get('bacmet2', []),
                'plasmid_genes': plasmid_by_sample.get(sample, [])
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
        output_file = self.output_dir / "genius_acinetobacter_ultimate_report.json"
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
                amr_rows.append({'Gene':g['gene'], 'Database':g['database'], 'Frequency':g['frequency_display'], 'Genomes':';'.join(g['genomes'])})
        if amr_rows:
            pd.DataFrame(amr_rows).to_csv(self.output_dir / "amr_genes.csv", index=False)
        vir_rows = []
        for db, genes in data['gene_centric'].get('virulence_databases', {}).items():
            for g in genes:
                vir_rows.append({'Gene':g['gene'], 'Database':g['database'], 'Frequency':g['frequency_display'], 'Genomes':';'.join(g['genomes'])})
        if vir_rows:
            pd.DataFrame(vir_rows).to_csv(self.output_dir / "virulence_genes.csv", index=False)
        bac_rows = []
        for db, genes in data['gene_centric'].get('bacmet_databases', {}).items():
            for g in genes:
                bac_rows.append({'Gene':g['gene'], 'Frequency':g['frequency_display'], 'Genomes':';'.join(g['genomes'])})
        if bac_rows:
            pd.DataFrame(bac_rows).to_csv(self.output_dir / "bacmet_genes.csv", index=False)
        plas_rows = []
        for db, genes in data.get('plasmid_analysis', {}).get('plasmid_databases', {}).items():
            for g in genes:
                plas_rows.append({'Marker':g['plasmid_marker'], 'Frequency':g['frequency_display'], 'Genomes':';'.join(g['genomes'])})
        if plas_rows:
            pd.DataFrame(plas_rows).to_csv(self.output_dir / "plasmid_markers.csv", index=False)
    
    def run(self):
        print("="*80)
        print("🧠 GENIUS ACINETOBACTER BAUMANNII ULTIMATE REPORTER v2.0.0")
        print("="*80)
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
        print("\n✅ Analysis complete! Open genius_acinetobacter_ultimate_report.html in your browser.")
        return True

def main():
    parser = argparse.ArgumentParser(description='GENIUS Acinetobacter baumannii Ultimate Reporter v2.0.0')
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