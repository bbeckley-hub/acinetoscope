#!/usr/bin/env python3
"""
GENIUS ACINETOBACTER BAUMANNII ULTIMATE REPORTER - COMPREHENSIVE EDITION
Advanced HTML Parser with Gene-Centric Cross-Genome Analysis for A. baumannii
Author: Beckley Brown <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.0.0 - Comprehensive Gene-Centric Edition with Environmental Markers
Date: 2026-01-06
"""

import os
import sys
import json
import re
import glob
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

# HTML parsing
from bs4 import BeautifulSoup

class UltimateHTMLParser:
    """Ultimate HTML parser for AcinetoScope reports"""
    
    def __init__(self):
        # COMPREHENSIVE database list for A. baumannii
        self.all_databases = [
            # AMR databases
            'amrfinder', 'card', 'resfinder', 'argannot', 'megares', 'bacmet2',
            # Virulence databases
            'vfdb', 'victors', 'ecoli_vf',
            # Plasmid databases
            'plasmidfinder', 'ecoh',
            # Additional
            'ncbi'
        ]
        
        # Database name mapping
        self.db_name_mapping = {
            'acineto_amrfinder': 'amrfinder',
            'acineto_card': 'card',
            'acineto_resfinder': 'resfinder',
            'acineto_argannot': 'argannot',
            'acineto_megares': 'megares',
            'acineto_bacmet2': 'bacmet2',
            'acineto_vfdb': 'vfdb',
            'acineto_victors': 'victors',
            'acineto_ecoli_vf': 'ecoli_vf',
            'acineto_plasmidfinder': 'plasmidfinder',
            'acineto_ecoh': 'ecoh',
            'acineto_ncbi': 'ncbi'
        }
    
    def normalize_sample_id(self, sample_id: str) -> str:
        """Normalize sample ID - FIXED FOR CORRECT MATCHING"""
        sample = str(sample_id)
        
        # Remove path if present
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        
        # Remove ALL common extensions but keep assembly ID
        # Handle .fna, .fasta, etc. from MLST files
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
                break
        
        # Convert GCF to GCA
        if sample.startswith('GCF_'):
            sample = 'GCA_' + sample[4:]
        
        return sample.strip()
    
    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        """Parse HTML table - IMPROVED FOR BETTER PARSING"""
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            
            table = tables[table_index]
            rows = table.find_all('tr')
            
            headers = []
            for th in rows[0].find_all(['th', 'td']):
                headers.append(th.get_text().strip())
            
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    # Handle cases where row has fewer columns than headers
                    while len(row_data) < len(headers):
                        row_data.append('')
                    if len(row_data) > len(headers):
                        row_data = row_data[:len(headers)]
                    data.append(row_data)
            
            if not data:
                return pd.DataFrame()
            
            df = pd.DataFrame(data, columns=headers)
            
            # Clean column names
            df.columns = [col.strip().replace('\n', ' ') for col in df.columns]
            
            return df
            
        except Exception as e:
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()
    
    def parse_mlst_report(self, file_path: Path, scheme: str = "pasteur") -> Dict[str, Dict]:
        """Parse MLST HTML report - FIXED FOR CORRECT COLUMN SELECTION"""
        print(f"  🧬 Parsing {scheme.upper()} MLST: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            # Parse the table
            df = self.parse_html_table(html_content, 0)
            if df.empty:
                # Try different table index
                df = self.parse_html_table(html_content, 1)
                if df.empty:
                    return {}
            
            # Clean column names - handle various column name formats
            df.columns = [col.strip().replace('\n', ' ').replace('_', ' ') for col in df.columns]
            
            # FIX 1: Find sample column CORRECTLY - NOT the first column (#)
            sample_col = None
            for col in df.columns:
                col_lower = col.lower()
                # Look for "Sample" column, NOT "#" column
                if 'sample' in col_lower and '#' not in col_lower:
                    sample_col = col
                    break
            
            # If "Sample" column not found, try other variations
            if not sample_col:
                for col in df.columns:
                    col_lower = col.lower()
                    if any(keyword in col_lower for keyword in ['genome', 'id', 'strain', 'file']):
                        sample_col = col
                        break
            
            # If still not found, use second column (skip # column)
            if not sample_col and len(df.columns) > 1:
                # Skip the first column if it looks like "#" or numbers
                first_col = df.columns[0].lower()
                if first_col in ['#', 'no', 'number', 'index'] or first_col.isdigit():
                    sample_col = df.columns[1]
                else:
                    sample_col = df.columns[0]
            
            if not sample_col:
                print(f"    ⚠️ Could not find sample column in {scheme} MLST")
                return {}
            
            print(f"    Using sample column: '{sample_col}'")
            
            # Normalize ALL sample names
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                if pd.isna(sample) or not sample:
                    continue
                
                # Extract ST - FIXED EXTRACTION (remove "ST" prefix)
                st = 'ND'
                
                # Try different column names for ST
                st_col_names = ['ST', 'st', 'Sequence Type', 'Sequence type', 'Type']
                for st_col in st_col_names:
                    if st_col in df.columns and pd.notna(row.get(st_col)):
                        st_val = str(row[st_col]).strip()
                        if st_val and st_val.lower() not in ['', 'nan', 'none', 'nd', 'unknown']:
                            # Clean ST value - remove "ST" prefix
                            if st_val.startswith('ST'):
                                st = st_val[2:]  # Remove "ST" prefix
                            else:
                                st = st_val
                            
                            # Handle STUNKNOWN
                            if 'unknown' in st.lower():
                                st = 'ND'
                            break
                
                # If no ST column found, try to extract from any column
                if st == 'ND':
                    for col in df.columns:
                        if col.lower() != sample_col.lower() and pd.notna(row.get(col)):
                            cell_val = str(row[col])
                            # Look for ST pattern in any cell
                            st_match = re.search(r'ST(\d+|UNKNOWN)', cell_val, re.I)
                            if st_match:
                                st_val = st_match.group(1)
                                if st_val.upper() == 'UNKNOWN':
                                    st = 'ND'
                                else:
                                    st = st_val
                                break
                
                # Extract International Clone - FIXED EXTRACTION (clean spaces)
                ic = 'Unknown'
                
                # Try different column names for International Clone
                ic_col_names = ['International Clone', 'International clone', 'IC', 'Clone']
                for ic_col in ic_col_names:
                    if ic_col in df.columns and pd.notna(row.get(ic_col)):
                        ic_val = str(row[ic_col]).strip()
                        if ic_val and ic_val.lower() not in ['', 'nan', 'none', 'unknown', 'not assigned']:
                            # Clean IC value - remove spaces
                            ic_match = re.search(r'(IC\s*[I|II|III|IV|V|VI|VII|VIII|IX|X]+|International Clone\s*[I|II|III|IV|V|VI|VII|VIII|IX|X]+)', ic_val, re.I)
                            if ic_match:
                                ic = ic_match.group(1)
                                # Remove "International Clone" and spaces
                                ic = ic.replace('International Clone', 'IC').replace(' ', '')
                            else:
                                ic = ic_val.replace(' ', '')  # Remove spaces
                            break
                
                # If no IC column found, check for IC patterns in any column
                if ic == 'Unknown':
                    for col in df.columns:
                        if col.lower() != sample_col.lower() and pd.notna(row.get(col)):
                            cell_val = str(row[col])
                            # Look for IC pattern
                            ic_match = re.search(r'(IC\s*[I|II|III|IV|V|VI|VII|VIII|IX|X]+)', cell_val, re.I)
                            if ic_match:
                                ic = ic_match.group(1).replace(' ', '')
                                break
                
                # Extract allele profile
                allele_profile = 'ND'
                allele_col_names = ['Allele Profile', 'Allele profile', 'Profile', 'Alleles']
                for allele_col in allele_col_names:
                    if allele_col in df.columns and pd.notna(row.get(allele_col)):
                        allele_val = str(row[allele_col])
                        if allele_val and allele_val.lower() not in ['', 'nan', 'none', 'nd']:
                            allele_profile = allele_val.strip()
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
            import traceback
            traceback.print_exc()
            return {}
    
    def parse_kaptive_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse Kaptive HTML report - FIXED FOR O LOCUS EXTRACTION"""
        print(f"  🧬 Parsing Kaptive: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            df = self.parse_html_table(html_content, 0)
            if df.empty:
                # Try different table index
                df = self.parse_html_table(html_content, 1)
                if df.empty:
                    return {}
            
            # Clean column names
            df.columns = [col.strip() for col in df.columns]
            
            # Find sample column - CORRECTLY identify "Genome" column
            sample_col = None
            for col in df.columns:
                col_lower = col.lower()
                if any(keyword in col_lower for keyword in ['genome', 'sample', 'id', 'strain']):
                    sample_col = col
                    break
            
            if not sample_col and len(df.columns) > 0:
                sample_col = df.columns[0]
            
            if not sample_col:
                return {}
            
            print(f"    Using sample column: '{sample_col}'")
            
            # Normalize ALL sample names
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                if pd.isna(sample) or not sample:
                    continue
                
                # Extract K Locus - FIXED EXTRACTION
                k_locus = 'ND'
                
                # Try different column names for K Locus
                k_col_names = ['K Locus', 'K locus', 'K-Locus', 'K-locus', 'K', 'KLocus']
                for k_col in k_col_names:
                    if k_col in df.columns and pd.notna(row.get(k_col)):
                        k_val = str(row[k_col]).strip()
                        if k_val and k_val.lower() not in ['', 'nan', 'none', 'nd', 'unknown']:
                            # Extract K type from string like "K15 (99.98%, 100.00%)" or "unknown (KL18) (100.00%, 100.00%)"
                            
                            # First try to extract K with number
                            k_match = re.search(r'(K\d+)', k_val, re.I)
                            if k_match:
                                k_locus = k_match.group(1).upper()
                            else:
                                # Check for unknown with KL number
                                unknown_match = re.search(r'unknown\s*\(KL(\d+)\)', k_val, re.I)
                                if unknown_match:
                                    k_locus = f"K{unknown_match.group(1)}"
                                else:
                                    # Try to extract just the K type
                                    k_parts = k_val.split()
                                    for part in k_parts:
                                        if part.startswith('K') and part[1:].isdigit():
                                            k_locus = part.upper()
                                            break
                                        elif part.startswith('KL') and part[2:].isdigit():
                                            k_locus = 'K' + part[2:]
                                            break
                            
                            if k_locus != 'ND':
                                break
                
                # Extract O Locus - FIXED EXTRACTION for "unknown (OCL5)"
                o_locus = 'ND'
                
                # Try different column names for O Locus
                o_col_names = ['O Locus', 'O locus', 'O-Locus', 'O-locus', 'O', 'OLocus']
                for o_col in o_col_names:
                    if o_col in df.columns and pd.notna(row.get(o_col)):
                        o_val = str(row[o_col]).strip()
                        if o_val and o_val.lower() not in ['', 'nan', 'none', 'nd', 'unknown']:
                            # Extract O type from string like "OC1 (98.88%, 100.00%)" or "unknown (OCL5) (100.00%, 81.00%)"
                            
                            # First try to extract OC with number
                            o_match = re.search(r'(OC\d+)', o_val, re.I)
                            if o_match:
                                o_locus = o_match.group(1).upper()
                            else:
                                # Check for unknown with OCL number - FIX FOR "unknown (OCL5)"
                                unknown_match = re.search(r'unknown\s*\(OCL(\d+)\)', o_val, re.I)
                                if unknown_match:
                                    o_locus = f"OC{unknown_match.group(1)}"
                                else:
                                    # Try other patterns
                                    ocl_match = re.search(r'OCL(\d+)', o_val, re.I)
                                    if ocl_match:
                                        o_locus = f"OC{ocl_match.group(1)}"
                                    else:
                                        # Try to extract just the O type
                                        o_parts = o_val.split()
                                        for part in o_parts:
                                            if part.startswith('OC') and part[2:].isdigit():
                                                o_locus = part.upper()
                                                break
                                            elif part.startswith('O') and part[1:].isdigit():
                                                o_locus = 'OC' + part[1:]
                                                break
                            
                            if o_locus != 'ND':
                                break
                
                # Create capsule type
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
            import traceback
            traceback.print_exc()
            return {}    
    
    def parse_amrfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """Parse AMRfinder HTML report - WITH PERCENTAGE CALCULATION"""
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            if len(tables) < 2:
                return {}, {}
            
            # Parse Gene Frequency table (second table)
            gene_frequencies = {}
            df_freq = self.parse_html_table(str(tables[1]), 0)
            
            if not df_freq.empty and 'Gene' in df_freq.columns:
                for _, row in df_freq.iterrows():
                    gene = str(row['Gene']).strip()
                    if not gene:
                        continue
                    
                    # Get frequency count
                    count = 0
                    frequency_str = str(row.get('Frequency', '0')).strip()
                    
                    # Extract count from frequency string
                    match = re.search(r'(\d+)', frequency_str)
                    if match:
                        count = int(match.group(1))
                    
                    # Calculate percentage
                    percentage = 0
                    if total_samples > 0:
                        percentage = (count / total_samples) * 100
                    
                    # Get genomes list
                    genomes = []
                    if 'Genomes' in df_freq.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip()) 
                                  for g in genomes_str.split(',') if g.strip()]
                    
                    # Get risk level
                    risk_level = 'Standard'
                    if 'Risk Level' in df_freq.columns and pd.notna(row.get('Risk Level')):
                        risk_level = str(row['Risk Level'])
                    
                    gene_frequencies[gene] = {
                        'count': count,
                        'percentage': round(percentage, 2),
                        'frequency_display': f"{count} ({percentage:.1f}%)",  # Combined format
                        'genomes': genomes,
                        'risk_level': risk_level,
                        'database': 'amrfinder'
                    }
            
            # Parse Genes by Genome table (first table)
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
                        # Clean and split genes - handle the "Showing X genes (scroll to see all)" text
                        genes = [g.strip() for g in gene_str.split() 
                                if g.strip() and not g.startswith('Showing') and not '(' in g and not ')' in g]
                    
                    genes_by_genome[sample] = genes
            
            print(f"    ✓ Found {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
            
        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder: {e}")
            import traceback
            traceback.print_exc()
            return {}, {}
    
    def parse_abricate_database_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """Parse ABRicate database HTML report - FIXED PARSING"""
        print(f"  🧬 Parsing ABRicate: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            if len(tables) < 2:
                return {}, {}
            
            # Determine database name
            db_name = 'unknown'
            filename = str(file_path.name).lower()
            for key, value in self.db_name_mapping.items():
                if key in filename:
                    db_name = value
                    break
            
            # Parse Gene Frequency table (second table)
            gene_frequencies = {}
            df_freq = self.parse_html_table(str(tables[1]), 0)
            
            if not df_freq.empty and 'Gene' in df_freq.columns:
                for _, row in df_freq.iterrows():
                    gene_full = str(row['Gene']).strip()
                    if not gene_full:
                        continue
                    
                    # Clean gene name (remove category prefix like "(AGly)")
                    gene = re.sub(r'^\([^)]+\)', '', gene_full).strip()
                    if not gene:
                        gene = gene_full
                    
                    # Get count
                    count = 0
                    frequency_str = str(row.get('Frequency', '0')).strip()
                    
                    # Extract count
                    match = re.search(r'(\d+)', frequency_str)
                    if match:
                        count = int(match.group(1))
                    
                    # Calculate percentage
                    percentage = 0
                    if total_samples > 0:
                        percentage = (count / total_samples) * 100
                    
                    # Get genomes list
                    genomes = []
                    if 'Genomes' in df_freq.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip()) 
                                  for g in genomes_str.split(',') if g.strip()]
                    
                    gene_frequencies[gene] = {
                        'count': count,
                        'percentage': round(percentage, 2),
                        'frequency_display': f"{count} ({percentage:.1f}%)",  # Combined format
                        'genomes': genomes,
                        'database': db_name,
                        'full_name': gene_full
                    }
            
            # Parse Genes by Genome table (first table)
            genes_by_genome = {}
            df_genomes = self.parse_html_table(str(tables[0]), 0)
            
            if not df_genomes.empty:
                # Find sample column
                sample_col = None
                for col in df_genomes.columns:
                    col_lower = col.lower()
                    if any(keyword in col_lower for keyword in ['genome', 'sample', 'id']):
                        sample_col = col
                        break
                
                if sample_col:
                    for _, row in df_genomes.iterrows():
                        sample = self.normalize_sample_id(row[sample_col])
                        if not sample:
                            continue
                        
                        genes = []
                        
                        # Find genes column
                        genes_col = None
                        for col in df_genomes.columns:
                            col_lower = col.lower()
                            if any(keyword in col_lower for keyword in ['genes', 'detected']):
                                genes_col = col
                                break
                        
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            # Clean and split genes
                            genes = [re.sub(r'^\([^)]+\)', '', g).strip() 
                                    for g in gene_str.split(',') if g.strip()]
                        
                        genes_by_genome[sample] = genes
            
            print(f"    ✓ {db_name.upper()}: {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
            
        except Exception as e:
            print(f"    ❌ Error parsing ABRicate report: {e}")
            import traceback
            traceback.print_exc()
            return {}, {}

    def parse_plasmidfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """Parse PlasmidFinder HTML report - SPECIALIZED for plasmid analysis"""
        print(f"  🧬 Parsing PlasmidFinder: {file_path.name}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            
            if len(tables) < 2:
                return {}, {}
            
            # Parse Gene Frequency table (second table)
            plasmid_frequencies = {}
            df_freq = self.parse_html_table(str(tables[1]), 0)
            
            if not df_freq.empty and 'Gene' in df_freq.columns:
                for _, row in df_freq.iterrows():
                    gene_full = str(row['Gene']).strip()
                    if not gene_full:
                        continue
                    
                    # Clean gene name for PlasmidFinder specific patterns
                    gene = self._clean_plasmid_gene_name(gene_full)
                    
                    # Get count
                    count = 0
                    frequency_str = str(row.get('Frequency', '0')).strip()
                    
                    # Extract count
                    match = re.search(r'(\d+)', frequency_str)
                    if match:
                        count = int(match.group(1))
                    
                    # Calculate percentage
                    percentage = 0
                    if total_samples > 0:
                        percentage = (count / total_samples) * 100
                    
                    # Get genomes list
                    genomes = []
                    if 'Genomes' in df_freq.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip()) 
                                for g in genomes_str.split(',') if g.strip()]
                    
                    # Categorize plasmid type
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
            
            # Parse Genes by Genome table (first table)
            plasmids_by_genome = {}
            df_genomes = self.parse_html_table(str(tables[0]), 0)
            
            if not df_genomes.empty:
                # Find sample column
                sample_col = None
                for col in df_genomes.columns:
                    col_lower = col.lower()
                    if any(keyword in col_lower for keyword in ['genome', 'sample', 'id']):
                        sample_col = col
                        break
                
                if sample_col:
                    for _, row in df_genomes.iterrows():
                        sample = self.normalize_sample_id(row[sample_col])
                        if not sample:
                            continue
                        
                        plasmids = []
                        
                        # Find genes column
                        genes_col = None
                        for col in df_genomes.columns:
                            col_lower = col.lower()
                            if any(keyword in col_lower for keyword in ['genes', 'detected']):
                                genes_col = col
                                break
                        
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            # Clean and split plasmids
                            plasmids = [self._clean_plasmid_gene_name(g.strip()) 
                                    for g in gene_str.split(',') if g.strip()]
                        
                        plasmids_by_genome[sample] = plasmids
            
            print(f"    ✓ PlasmidFinder: {len(plasmids_by_genome)} samples, {len(plasmid_frequencies)} plasmid markers")
            return plasmids_by_genome, plasmid_frequencies
            
        except Exception as e:
            print(f"    ❌ Error parsing PlasmidFinder: {e}")
            import traceback
            traceback.print_exc()
            return {}, {}

    def _clean_plasmid_gene_name(self, gene_name: str) -> str:
        """Clean plasmid gene names - special handling for PlasmidFinder patterns"""
        # Remove trailing _1, _2, etc. (common in PlasmidFinder)
        gene = re.sub(r'_(\d+)$', '', gene_name)
        
        # Extract plasmid name in parentheses if present
        plasmid_match = re.search(r'\((.*?)\)', gene)
        if plasmid_match:
            plasmid_name = plasmid_match.group(1)
            # Keep the main gene name + plasmid in parentheses
            base_gene = re.sub(r'\(.*?\)', '', gene).strip()
            if base_gene:
                gene = f"{base_gene}({plasmid_name})"
            else:
                gene = plasmid_name
        
        return gene.strip()

    def _categorize_plasmid(self, gene_name: str) -> str:
        """Categorize plasmid markers"""
        gene_lower = gene_name.lower()
        
        # Plasmid replication genes
        if any(rep in gene_lower for rep in ['rep', 'inc', 'rep_']):
            if 'coli' in gene_lower or 'col' in gene_lower:
                return 'Colicin plasmid'
            elif 'broad' in gene_lower or 'multihost' in gene_lower:
                return 'Broad-host-range plasmid'
            else:
                return 'Replication protein'
        
        # Mobility genes
        elif any(mob in gene_lower for mob in ['mob', 'tra', 'conj']):
            return 'Mobility/conjugation'
        
        # Colicin plasmids
        elif 'col' in gene_lower and 'coli' not in gene_lower:
            if gene_lower.startswith('col'):
                return 'Colicin plasmid'
            else:
                return 'Other plasmid'
        
        # Incompatibility groups
        elif 'inc' in gene_lower and gene_lower.startswith('inc'):
            inc_match = re.search(r'inc([a-z0-9]+)', gene_lower)
            if inc_match:
                inc_type = inc_match.group(1).upper()
                return f'Incompatibility group {inc_type}'
            return 'Incompatibility group'
        
        # Specific plasmid families
        elif any(fam in gene_lower for fam in ['pneumoniae', 'salmonella', 'enterobacter', 'klebsiella']):
            return 'Enterobacteriales plasmid'
        
        elif any(fam in gene_lower for fam in ['pseudomonas', 'aeruginosa']):
            return 'Pseudomonas plasmid'
        
        elif any(fam in gene_lower for fam in ['acinetobacter', 'baumannii']):
            return 'Acinetobacter plasmid'
        
        elif any(fam in gene_lower for fam in ['staphylococcus', 'aureus', 'mrsa']):
            return 'Staphylococcus plasmid'
        
        else:
            return 'Other plasmid'


class UltimateDataAnalyzer:
    """Analyzes data for ultimate gene-centric reporting for A. baumannii"""
    
    def __init__(self):
        # COMPREHENSIVE CRITICAL GENE LIST FOR A. BAUMANNII
        self.critical_carbapenemases = {
            # OXA-type (most common in A. baumannii)
            'blaOXA-23', 'blaOXA-24', 'blaOXA-40', 'blaOXA-51', 'blaOXA-58',
            'blaOXA-66', 'blaOXA-69', 'blaOXA-71', 'blaOXA-143', 'blaOXA-235',
            'blaOXA-236', 'blaOXA-237', 'blaOXA-267', 'blaOXA-317', 'blaOXA-91',
            # Metallo-β-lactamases
            'blaNDM', 'blaNDM-1', 'blaNDM-2', 'blaNDM-3', 'blaNDM-4', 'blaNDM-5',
            'blaVIM', 'blaVIM-1', 'blaVIM-2', 'blaVIM-3', 'blaVIM-4',
            'blaIMP', 'blaIMP-1', 'blaIMP-2', 'blaIMP-3', 'blaIMP-4', 'blaIMP-5',
            'blaKPC', 'blaKPC-2', 'blaKPC-3', 'blaKPC-4',
            # Other carbapenemases
            'blaGES', 'blaGES-1', 'blaGES-2', 'blaGES-5', 'blaGES-14',
            'blaSIM', 'blaSPM', 'blaAIM'
        }

        self.critical_esbls = {
            # CTX-M family (ALL are ESBLs)
            'blaCTX-M', 'blaCTX-M-1', 'blaCTX-M-2', 'blaCTX-M-3', 'blaCTX-M-4', 'blaCTX-M-5',
            'blaCTX-M-6', 'blaCTX-M-7', 'blaCTX-M-8', 'blaCTX-M-9', 'blaCTX-M-10', 'blaCTX-M-11',
            'blaCTX-M-12', 'blaCTX-M-13', 'blaCTX-M-14', 'blaCTX-M-15', 'blaCTX-M-16', 'blaCTX-M-17',
            'blaCTX-M-18', 'blaCTX-M-19', 'blaCTX-M-20', 'blaCTX-M-21', 'blaCTX-M-22', 'blaCTX-M-23',
            'blaCTX-M-24', 'blaCTX-M-25', 'blaCTX-M-26', 'blaCTX-M-27', 'blaCTX-M-28', 'blaCTX-M-29',
            'blaCTX-M-30', 'blaCTX-M-31', 'blaCTX-M-32', 'blaCTX-M-33', 'blaCTX-M-34', 'blaCTX-M-35',
            'blaCTX-M-36', 'blaCTX-M-37', 'blaCTX-M-38', 'blaCTX-M-39', 'blaCTX-M-40', 'blaCTX-M-41',
            'blaCTX-M-42', 'blaCTX-M-43', 'blaCTX-M-44', 'blaCTX-M-45', 'blaCTX-M-46', 'blaCTX-M-47',
            'blaCTX-M-48', 'blaCTX-M-49', 'blaCTX-M-50', 'blaCTX-M-51', 'blaCTX-M-52', 'blaCTX-M-53',
            'blaCTX-M-54', 'blaCTX-M-55', 'blaCTX-M-56', 'blaCTX-M-57', 'blaCTX-M-58', 'blaCTX-M-59',
            'blaCTX-M-60', 'blaCTX-M-61', 'blaCTX-M-62', 'blaCTX-M-63', 'blaCTX-M-64', 'blaCTX-M-65',
            'blaCTX-M-66', 'blaCTX-M-67', 'blaCTX-M-68', 'blaCTX-M-69', 'blaCTX-M-70', 'blaCTX-M-71',
            'blaCTX-M-72', 'blaCTX-M-73', 'blaCTX-M-74', 'blaCTX-M-75', 'blaCTX-M-76', 'blaCTX-M-77',
            'blaCTX-M-78', 'blaCTX-M-79', 'blaCTX-M-80', 'blaCTX-M-81', 'blaCTX-M-82', 'blaCTX-M-83',
            'blaCTX-M-84', 'blaCTX-M-85', 'blaCTX-M-86', 'blaCTX-M-87', 'blaCTX-M-88', 'blaCTX-M-89',
            'blaCTX-M-90', 'blaCTX-M-91', 'blaCTX-M-92', 'blaCTX-M-93', 'blaCTX-M-94', 'blaCTX-M-95',
            'blaCTX-M-96', 'blaCTX-M-97', 'blaCTX-M-98', 'blaCTX-M-99', 'blaCTX-M-100', 'blaCTX-M-101',
            'blaCTX-M-102', 'blaCTX-M-103', 'blaCTX-M-104', 'blaCTX-M-105', 'blaCTX-M-106',
            
            # SHV family (ONLY ESBL variants)
            'blaSHV-2', 'blaSHV-3', 'blaSHV-4', 'blaSHV-5', 'blaSHV-7', 'blaSHV-8', 'blaSHV-9',
            'blaSHV-10', 'blaSHV-12', 'blaSHV-13', 'blaSHV-14', 'blaSHV-15', 'blaSHV-16', 'blaSHV-17',
            'blaSHV-18', 'blaSHV-19', 'blaSHV-20', 'blaSHV-21', 'blaSHV-22', 'blaSHV-23', 'blaSHV-24',
            'blaSHV-25', 'blaSHV-26', 'blaSHV-27', 'blaSHV-28', 'blaSHV-29', 'blaSHV-30', 'blaSHV-31',
            'blaSHV-32', 'blaSHV-33', 'blaSHV-34', 'blaSHV-35', 'blaSHV-36', 'blaSHV-37', 'blaSHV-38',
            'blaSHV-39', 'blaSHV-40', 'blaSHV-41', 'blaSHV-42', 'blaSHV-43', 'blaSHV-44', 'blaSHV-45',
            'blaSHV-46', 'blaSHV-47', 'blaSHV-48', 'blaSHV-49', 'blaSHV-50', 'blaSHV-51', 'blaSHV-52',
            'blaSHV-53', 'blaSHV-54', 'blaSHV-55', 'blaSHV-56', 'blaSHV-57', 'blaSHV-58', 'blaSHV-59',
            'blaSHV-60', 'blaSHV-61', 'blaSHV-62', 'blaSHV-63', 'blaSHV-64', 'blaSHV-65', 'blaSHV-66',
            'blaSHV-67', 'blaSHV-68', 'blaSHV-69', 'blaSHV-70', 'blaSHV-71', 'blaSHV-72', 'blaSHV-73',
            'blaSHV-74', 'blaSHV-75', 'blaSHV-76', 'blaSHV-77', 'blaSHV-78', 'blaSHV-79', 'blaSHV-80',
            'blaSHV-81', 'blaSHV-82', 'blaSHV-83', 'blaSHV-84', 'blaSHV-85', 'blaSHV-86', 'blaSHV-87',
            'blaSHV-88', 'blaSHV-89', 'blaSHV-90', 'blaSHV-91', 'blaSHV-92', 'blaSHV-93', 'blaSHV-94',
            'blaSHV-95', 'blaSHV-96', 'blaSHV-97', 'blaSHV-98', 'blaSHV-99', 'blaSHV-100', 'blaSHV-101',
            'blaSHV-102', 'blaSHV-103', 'blaSHV-104', 'blaSHV-105', 'blaSHV-106', 'blaSHV-107',
            
            # TEM family (ONLY ESBL variants)
            'blaTEM-3', 'blaTEM-4', 'blaTEM-5', 'blaTEM-6', 'blaTEM-7', 'blaTEM-8', 'blaTEM-9',
            'blaTEM-10', 'blaTEM-12', 'blaTEM-16', 'blaTEM-19', 'blaTEM-20', 'blaTEM-24', 'blaTEM-26',
            'blaTEM-29', 'blaTEM-30', 'blaTEM-32', 'blaTEM-34', 'blaTEM-36', 'blaTEM-39', 'blaTEM-40',
            'blaTEM-41', 'blaTEM-43', 'blaTEM-44', 'blaTEM-45', 'blaTEM-46', 'blaTEM-47', 'blaTEM-49',
            'blaTEM-50', 'blaTEM-51', 'blaTEM-52', 'blaTEM-55', 'blaTEM-56', 'blaTEM-59', 'blaTEM-63',
            'blaTEM-64', 'blaTEM-65', 'blaTEM-68', 'blaTEM-70', 'blaTEM-72', 'blaTEM-74', 'blaTEM-76',
            'blaTEM-78', 'blaTEM-79', 'blaTEM-80', 'blaTEM-81', 'blaTEM-83', 'blaTEM-84', 'blaTEM-85',
            'blaTEM-86', 'blaTEM-87', 'blaTEM-88', 'blaTEM-89', 'blaTEM-90', 'blaTEM-91', 'blaTEM-92',
            'blaTEM-93', 'blaTEM-94', 'blaTEM-95', 'blaTEM-96', 'blaTEM-97', 'blaTEM-98', 'blaTEM-99',
            'blaTEM-100', 'blaTEM-101', 'blaTEM-102', 'blaTEM-103', 'blaTEM-104', 'blaTEM-105',
            'blaTEM-106', 'blaTEM-107', 'blaTEM-108',
            
            # Other ESBL families
            'blaPER', 'blaPER-1', 'blaPER-2', 'blaPER-3', 'blaPER-4', 'blaPER-5', 'blaPER-6',
            'blaPER-7', 'blaPER-8', 'blaPER-9', 'blaPER-10', 'blaPER-11', 'blaPER-12', 'blaPER-13',
            'blaPER-14', 'blaPER-15', 'blaPER-16', 'blaPER-17', 'blaPER-18', 'blaPER-19', 'blaPER-20',
            
            'blaVEB', 'blaVEB-1', 'blaVEB-2', 'blaVEB-3', 'blaVEB-4', 'blaVEB-5', 'blaVEB-6',
            'blaVEB-7', 'blaVEB-8', 'blaVEB-9', 'blaVEB-10', 'blaVEB-11', 'blaVEB-12', 'blaVEB-13',
            
            'blaBEL', 'blaBEL-1', 'blaBEL-2', 'blaBEL-3', 'blaBEL-4', 'blaBEL-5', 'blaBEL-6',
            
            'blaGES', 'blaGES-1', 'blaGES-2', 'blaGES-3', 'blaGES-4', 'blaGES-5', 'blaGES-6',
            'blaGES-7', 'blaGES-8', 'blaGES-9', 'blaGES-10', 'blaGES-11', 'blaGES-12', 'blaGES-13',
            'blaGES-14', 'blaGES-15', 'blaGES-16', 'blaGES-17', 'blaGES-18', 'blaGES-19', 'blaGES-20',
            
            # Additional ESBL families
            'blaSFO', 'blaTLA', 'blaBES', 'blaSCO'
        }  
              
        self.critical_ampc = {
            'blaADC', 'blaADC-1', 'blaADC-2', 'blaADC-5', 'blaADC-7', 'blaADC-10',
            'blaADC-11', 'blaADC-30', 'blaADC-69', 'blaADC-75', 'blaADC-88', 'blaADC-176'
        }
        
        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'pmrA', 'pmrB', 'pmrC', 'pmrE', 'pmrF',
            'lpxA', 'lpxC', 'lpxD', 'lpxL', 'lpxM',
            'eptA', 'arnA', 'arnB', 'arnC', 'arnD', 'arnE', 'arnF',
            'pagP', 'phoP', 'phoQ'
        }
        
        self.critical_tigecycline = {
            'tet(X)', 'tet(X1)', 'tet(X2)', 'tet(X3)', 'tet(X4)', 'tet(X5)', 'tet(X6)',
            'tet(39)', 'tet(A)', 'tet(B)', 'tet(C)', 'tet(D)', 'tet(E)', 'tet(G)', 'tet(H)',
            'tet(J)', 'tet(L)', 'tet(M)', 'tet(O)', 'tet(Q)', 'tet(S)', 'tet(W)',
            'adeS', 'adeR', 'adeA', 'adeB', 'adeC', 'adeJ', 'adeK', 'adeN', 'adeT'
        }
        
        self.critical_biofilm = {
            'ompA', 'csuA', 'csuB', 'csuC', 'csuD', 'csuE', 'csuA/B', 'csuC/D/E',
            'bfmR', 'bfmS', 'abaI', 'abaR', 'pilA', 'pilB', 'pilC', 'pilD', 'pilE', 'pilF',
            'ptk', 'epsA', 'pgaA', 'pgaB', 'pgaC', 'pgaD', 'bap'
        }
        
        self.critical_efflux = {
            'adeA', 'adeB', 'adeC', 'adeF', 'adeG', 'adeH', 'adeI', 'adeJ', 'adeK', 'adeL', 'adeM', 'adeN',
            'abeM', 'abeS', 'amvA', 'craA' 'adeT1', 'adeT2', 'amvA', 'mexJ', 'mexK', 'mexT', 'mdeA', 'mdfA/cmr', 'mdtN/yjcR'
        }
        
        # ============================================================================
        # ENVIRONMENTAL RESISTANCE & CO-SELECTION MARKERS - SCIENTIFICALLY ACCURATE
        # ============================================================================
        
        # 1. TRUE BACMET2 MARKERS (Biocide & Metal Resistance only)
        # Based on actual BACMET2 database (http://bacmet.biomedicine.gu.se/)
        self.bacmet2_markers = {
            # Biocide resistance (quaternary ammonium compounds)
            'qacA', 'qacB', 'qacC', 'qacD', 'qacE', 'qacF', 'qacG', 'qacH', 'qacI', 'qacJ',
            'qacEA1', 'qacG2', 'qacH2', 'qacE', 'qacEdelta1',
            # Chlorhexidine resistance
            'cepA',
            # Formaldehyde resistance
            'formA', 'formB', 'formC',
            # Phenolic compound resistance
            'oqxA', 'oqxB',
            
            # Heavy metal resistance (from BACMET2)
            'czcA', 'czcB', 'czcC', 'czcD', 'czcR', 'czcS',  # Cadmium/Zinc/Cobalt
            'merA', 'merB', 'merC', 'merD', 'merE', 'merF', 'merG', 'merH', 'merI', 'merJ', 'merP', 'merT',  # Mercury
            'arsA', 'arsB', 'arsC', 'arsD', 'arsE', 'arsF', 'arsG', 'arsH', 'arsI', 'arsJ',  # Arsenic
            'arsT',  
            'copA', 'copB', 'copC', 'copD', 'copE', 'copF', 'copG', 'copH', 'copI', 'copJ',  # Copper
            'zntA', 'zntB', 'zntC', 'zntD', 'zntE', 'zntF', 'zntG', 'zntH', 'zntI', 'zntJ',  # Zinc
            'czcR', 'czcS', 'cadR',          
            'chrA', 'chrB', 'chrC', 'chrD', 'chrE', 'chrF',  # Chromate resistance
            'nikA', 'nikB', 'nikC', 'nikD', 'nikE', 'nikR',   # Nickel resistance
            'cadA', 'cadB', 'cadC', 'cadD',  # Cadmium resistance
            'silA', 'silB', 'silC', 'silD', 'silE',  # Silver resistance
            'pbrA', 'pbrB', 'pbrC', 'pbrD', 'pbrR',  # Lead resistance
            'corA',          # Magnesium/cobalt transporter
            'corC',          # Magnesium/cobalt transporter
            'corR',          # Magnesium/cobalt regulator
            'zraR/hydH',     # Zinc resistance
            'pitA',          # Phosphate/arsenate transporter
            'nccN',          # Nickel/cobalt/cadmium resistance
            'nreB',          # Nickel/cobalt resistance regulator
            'fptA',          # Iron transporter
            'fecE',          # Iron transport
            'fpvA',          # Iron transport
            'znuB', 'znuC',  # Zinc uptake
            'frnE',          # Iron transport
        }
        
        # 2. ENVIRONMENTAL CO-SELECTION MARKERS (Stress response, plasmid transfer, etc.)
        # These genes facilitate co-selection in hospital environments
        self.environmental_co_selection = {
            # Global regulators that induce efflux pumps
            'soxR', 'soxS', 'marA', 'marB', 'marC', 'marR', 'robA',
            'rpoS',  # Stress sigma factor
            'rpoH',  # Heat shock sigma factor
            'oxyRkp', 'cpxR', 'baeR', 'emrAsm', 'emrBsm', 'yddg/emrE', 'lmrS', 'smeF', 'sugE', 
            
            # Nutrient starvation responses
            'phoB', 'phoR', 'phoU',  # Phosphate starvation → efflux induction
            
            # Conjugation and mobilization genes (plasmid spread)
            'traA', 'traB', 'traC', 'traD', 'traE', 'traF', 'traG', 'traH', 'traI', 'traJ',
            'traK', 'traL', 'traM', 'traN', 'traO', 'traP', 'traQ', 'traR', 'traS', 'traT',
            'traU', 'traV', 'traW', 'traX', 'traY',
            'mobA', 'mobB', 'mobC', 'mobD', 'mobE', 'mobF', 'mobG', 'mobH',
            'oriT', 'oriV', 'repA', 'repB', 'repC',
            
            # Integrons and transposons
            'intI1', 'intI2', 'intI3',  # Integrase genes
            'tnpA', 'tnpB', 'tnpC', 'tnpD', 'tnpE', 'tnpF',  # Transposon genes
            'istA', 'istB',  # Insertion sequence elements
        }
        
        # 3. ENVIRONMENTAL ANTIBIOTIC RESISTANCE MARKERS (NOT BACMET2)
        # These are antibiotic resistance genes commonly found in environmental reservoirs
        self.common_antibiotic_markers = {
            # Sulfonamide resistance 
            'sul1', 'sul2', 'sul3',
            # Trimethoprim resistance
            'dfrA1', 'dfrA5', 'dfrA7', 'dfrA12', 'dfrA14', 'dfrA17', 'dfrA19', 'dfrA21',
            'dfrB1', 'dfrB2', 'dfrB3', 'dfrB4',
            # Chloramphenicol resistance
            'catA1', 'catA2', 'catB2', 'catB3', 'catB8', 'catI', 'catII', 'catIII',
            # Aminoglycoside resistance patterns in environment
            'aac', 'aad', 'ant', 'aph',
            # Tetracycline resistance in environment (non-tetX variants)
            'tet', 'tetR', 'tetA', 'tetB', 'tetC', 'tetD', 'tetE', 'tetG', 'tetH',
            'tetJ', 'tetK', 'tetL', 'tetM', 'tetO', 'tetQ', 'tetS', 'tetW', 'tetX',
            # Macrolide resistance
            'erm', 'ere', 'mef', 'msr',
            # Multidrug transporters in environment
            'mdt', 'emr', 'acr', 'tolC',
            # Beta-lactamase genes in environmental bacteria
            'blaTEM', 'blaSHV', 'blaCTX-M', 'blaOXA',
        }
        
        # 4. VICTORS DATABASE MARKERS (Virulence Factors)
        self.victors_markers = {
            # Adhesion factors
            'fimA', 'fimB', 'fimC', 'fimD', 'fimE', 'fimF', 'fimG', 'fimH',
            'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH', 'papI', 'papJ', 'papK',
            'afa', 'dra',
            
            # Toxins
            'hlyA', 'hlyB', 'hlyC', 'hlyD',  # Hemolysin
            'cnf1', 'cnf2',  # Cytotoxic necrotizing factor
            'sat',  # Secreted autotransporter toxin
            'astA',  # Heat-stable enterotoxin
            'stx1', 'stx2',  # Shiga toxin
            
            # Immune evasion
            'iss',  # Increased serum survival
            'traT',  # Serum resistance
            'kpsM', 'kpsT',  # Capsule synthesis
            'ibeA', 'ibeB',  # Invasion of brain endothelium
            
            # Iron acquisition
            'iutA',  # Aerobactin receptor
            'iroN',  # Catecholate siderophore receptor
            'fyuA',  # Yersiniabactin receptor
            'irp1', 'irp2',  # Yersiniabactin synthesis
            'chu',  # Heme uptake
            
            # Other virulence factors
            'usp',  # Uropathogenic specific protein
            'vat',  # Vacuolating autotransporter toxin
            'pic',  # Protein involved in colonization
            'sigA',  # IgA protease-like protein
            'tia',  # Invasion adhesion
        }
      
        #5. Plasmid-specific categories
        self.plasmid_categories = {
            # Common Acinetobacter plasmids
            'pAB1', 'pAB2', 'pAB3', 'pAB4', 'pAB5', 'pAB6', 'pAB7', 'pAB8',
            'pACICU', 'pACICU1', 'pACICU2',
            'pACICU1-like', 'pACICU2-like',
            'pRAY', 'pRAY-like',
            'pABVA01', 'pABVA02', 'pABVA03',
            'pAC30', 'pAC31', 'pAC32',
            'pAB3-like', 'pAB5-like',
            'pACICU-A', 'pACICU-B',
            
            # Broad-host-range plasmids in A. baumannii
            'pNDM', 'pNDM-1', 'pNDM-2',
            'pOXA', 'pOXA-23', 'pOXA-58',
            'pISAba125', 'pISAba1',
            'pACICU-1', 'pACICU-2',
            
            # Colicin plasmids (common in Gram-negatives)
            'ColE1', 'ColE2', 'ColE3',
            'ColA', 'ColB', 'ColD', 'ColE', 'ColF', 'ColH', 'ColI', 'ColK', 'ColM', 'ColN', 'ColV',
            'ColRNAI', 'ColDF13',
            'Col(MG828)', 'Col(MP18)', 'Col8282', 'Col(pHAD28)',
            
            # Replication proteins
            'rep', 'repA', 'repB', 'repC', 'repD', 'repE',
            'rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6', 'rep7', 'rep8', 'rep9', 'rep10',
            'rep_1', 'rep_2', 'rep_3', 'rep_4', 'rep_5', 'rep_6', 'rep_7', 'rep_8', 'rep_9', 'rep_10',
            
            # Incompatibility groups
            'IncA', 'IncB', 'IncC', 'IncD', 'IncF', 'IncG', 'IncH', 'IncI', 'IncJ', 'IncK', 
            'IncL', 'IncM', 'IncN', 'IncO', 'IncP', 'IncQ', 'IncR', 'IncS', 'IncT', 'IncU', 
            'IncV', 'IncW', 'IncX', 'IncY', 'IncZ',
            'IncFIA', 'IncFIB', 'IncFIC', 'IncFII', 'IncHI1', 'IncHI2', 'IncL/M', 'IncN2',
            
            # Mobility genes
            'mobA', 'mobB', 'mobC', 'mobD', 'mobE',
            'traA', 'traB', 'traC', 'traD', 'traE', 'traF', 'traG', 'traH', 'traI', 'traJ',
            'traK', 'traL', 'traM', 'traN', 'traO', 'traP', 'traQ', 'traR', 'traS', 'traT',
            'traU', 'traV', 'traW', 'traX', 'traY',
            
            # Plasmid addiction systems
            'ccdA', 'ccdB',  # CcdAB toxin-antitoxin
            'vapB', 'vapC',  # VapBC toxin-antitoxin
            'relB', 'relE',  # RelBE toxin-antitoxin
            'hok', 'sok',    # Hok/Sok system
            'pemI', 'pemK',  # PemIK system
            'maxE', 'maxF',  # MaxEF system
            'mazE', 'mazF',  # MazEF system
        }
    
        # 6. OTHER GENERAL RESISTANCE MARKERS
        self.other_resistance_markers = {
            # Fluoroquinolone resistance
            'qnrA', 'qnrB', 'qnrC', 'qnrD', 'qnrS', 'qnrVC',
            'aac(6\')-Ib-cr',  # Aminoglycoside modifying enzyme with quinolone resistance
            'oqxa', 'oqxb',  # Quinolone resistance
            # Fosfomycin resistance
            'fosA', 'fosB', 'fosC', 'fosX',
            # Rifampicin resistance
            'rpoB',  # RNA polymerase beta subunit mutations
            # Polymyxin resistance (additional to critical_colistin)
            'mgrB', 'pmrAB', 'phoPQ', 'eptB',
        }
        
        # COMBINE ALL CRITICAL GENES
        self.critical_amr_genes = (
            self.critical_carbapenemases | self.critical_esbls | self.critical_ampc |
            self.critical_colistin | self.critical_tigecycline | self.critical_efflux
        )
        
        self.critical_virulence_genes = self.critical_biofilm
        
        # COMBINE ALL ENVIRONMENTAL & CO-SELECTION MARKERS
        self.all_environmental_markers = (
            self.bacmet2_markers | 
            self.environmental_co_selection |
            self.common_antibiotic_markers |
            self.victors_markers |
            self.other_resistance_markers
        )
    
    def categorize_gene(self, gene: str) -> str:
        """Categorize a gene with SCIENTIFICALLY ACCURATE environmental markers"""
        gene_lower = gene.lower()
        
        # Check categories in order of importance
        # 1. Critical clinical resistance (highest priority)
        if any(carb in gene_lower for carb in [g.lower() for g in self.critical_carbapenemases]):
            return 'Carbapenemases'
        elif any(esbl in gene_lower for esbl in [g.lower() for g in self.critical_esbls]):
            return 'ESBLs'
        elif any(ampc in gene_lower for ampc in [g.lower() for g in self.critical_ampc]):
            return 'AmpC'
        elif any(col in gene_lower for col in [g.lower() for g in self.critical_colistin]):
            return 'Colistin Resistance'
        elif any(tig in gene_lower for tig in [g.lower() for g in self.critical_tigecycline]):
            return 'Tigecycline Resistance'
        elif any(bio in gene_lower for bio in [g.lower() for g in self.critical_biofilm]):
            return 'Biofilm Formation'
        elif any(eff in gene_lower for eff in [g.lower() for g in self.critical_efflux]):
            return 'Efflux Pumps'
        
        # 2. True BACMET2 markers (biocide & metal resistance)
        elif any(bac in gene_lower for bac in [g.lower() for g in self.bacmet2_markers]):
            # Specifically check for biocide vs metal resistance
            if any(biocide in gene_lower for biocide in ['qac', 'cep', 'form', 'oqx']):
                return 'Biocide Resistance (BACMET2)'
            elif any(metal in gene_lower for metal in ['czc', 'mer', 'ars', 'cop', 'znt', 'chr', 'nik', 'cad', 'sil', 'pbr']):
                return 'Metal Resistance (BACMET2)'
            else:
                return 'BACMET2 Other'
        
        # 3. Environmental co-selection markers
        elif any(env in gene_lower for env in [g.lower() for g in self.environmental_co_selection]):
            if any(plasmid in gene_lower for plasmid in ['tra', 'mob', 'rep', 'ori']):
                return 'Plasmid Transfer'
            elif any(stress in gene_lower for stress in ['sox', 'mar', 'rob', 'rpo']):
                return 'Stress Response'
            elif any(mobile in gene_lower for mobile in ['int', 'tnp', 'ist']):
                return 'Mobile Genetic Elements'
            else:
                return 'Environmental Co-Selection'
        
        # 4. Environmental antibiotic resistance (NOT BACMET2)
        elif any(env_abx in gene_lower for env_abx in [g.lower() for g in self.common_antibiotic_markers]):
            if 'sul' in gene_lower:
                return 'Sulfonamide Resistance'
            elif 'dfr' in gene_lower:
                return 'Trimethoprim Resistance'
            elif 'cat' in gene_lower:
                return 'Chloramphenicol Resistance'
            elif any(ag in gene_lower for ag in ['aac', 'aad', 'ant', 'aph']):
                return 'Aminoglycoside Resistance'
            elif 'tet' in gene_lower and 'tet(x)' not in gene_lower:
                return 'Tetracycline Resistance'
            elif any(ml in gene_lower for ml in ['erm', 'mef', 'msr']):
                return 'Macrolide Resistance'
            elif 'bla' in gene_lower and 'carbapenem' not in gene_lower:
                return 'Beta-lactamase'
            else:
                return 'Antibiotic Resistance'
        
        # 5. VICTORS virulence markers
        elif any(vic in gene_lower for vic in [g.lower() for g in self.victors_markers]):
            if any(adhesion in gene_lower for adhesion in ['fim', 'pap', 'afa', 'dra']):
                return 'VICTORS Adhesion'
            elif any(toxin in gene_lower for toxin in ['hly', 'cnf', 'sat', 'ast', 'stx']):
                return 'VICTORS Toxins'
            elif any(iron in gene_lower for iron in ['iut', 'iro', 'fyu', 'irp', 'chu']):
                return 'VICTORS Iron Acquisition'
            else:
                return 'VICTORS Virulence'
        
        # 6. Other resistance markers
        elif any(other in gene_lower for other in [g.lower() for g in self.other_resistance_markers]):
            if 'qnr' in gene_lower:
                return 'Quinolone Resistance'
            elif 'fos' in gene_lower:
                return 'Fosfomycin Resistance'
            else:
                return 'Other Resistance'
        
        # 7. General antibiotic resistance patterns (catch-all)
        elif any(abx in gene_lower for abx in ['sul', 'dfr', 'cat', 'aac', 'aad', 'ant', 'aph', 'tet', 'erm']):
            return 'Antibiotic Resistance'
        
        # 8. General virulence patterns
        elif any(vir in gene_lower for vir in ['tox', 'hly', 'cnf', 'iss', 'fim', 'pap']):
            return 'Virulence Factors'
        
        else:
            return 'Other'
    
    def create_gene_centric_tables(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        """Create gene-centric tables with percentages"""
        gene_centric = {
            'amr_databases': {},
            'virulence_databases': {},
            'plasmid_databases': {},
            'environmental_databases': {},
            'combined_gene_frequencies': [],
            'gene_categories': defaultdict(list),
            'database_stats': {},
            'environmental_summary': defaultdict(dict)
        }
        
        # Process AMRfinder
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
                    'risk_level': data.get('risk_level', 'Standard'),
                    'genomes': data.get('genomes', [])
                })
            
            if gene_list:
                gene_list.sort(key=lambda x: x['count'], reverse=True)
                gene_centric['amr_databases']['amrfinder'] = gene_list
                
                # Add to categories
                for gene_data in gene_list:
                    gene_centric['gene_categories'][gene_data['category']].append(gene_data)
        
        # Process ABRicate databases
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
                
                # Sort and store
                if gene_list:
                    gene_list.sort(key=lambda x: x['count'], reverse=True)
                    
                    # Categorize database type
                    if db_name in ['vfdb', 'victors', 'ecoli_vf']:
                        gene_centric['virulence_databases'][db_name] = gene_list
                    elif db_name in ['plasmidfinder', 'ecoh']:
                        gene_centric['plasmid_databases'][db_name] = gene_list
                    elif db_name in ['bacmet2']:
                        gene_centric['environmental_databases'][db_name] = gene_list
                    else:
                        gene_centric['amr_databases'][db_name] = gene_list
                    
                    # Add to categories
                    for gene_data in gene_list:
                        gene_centric['gene_categories'][gene_data['category']].append(gene_data)
        
        # Create environmental summary
        self._create_environmental_summary(gene_centric, total_samples)
        
        # Create combined gene frequencies
        all_genes = []
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'environmental_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    all_genes.append(gene_data)
        
        # Sort combined list by count
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = all_genes
        
        # Sort categories by total count
        for category in gene_centric['gene_categories']:
            gene_centric['gene_categories'][category].sort(key=lambda x: x['count'], reverse=True)
        
        # Calculate database statistics
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'environmental_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                total_genes = len(genes)
                total_occurrences = sum(g['count'] for g in genes)
                
                # Count critical genes and environmental markers
                critical_count = 0
                bacmet2_count = 0
                environmental_antibiotic_count = 0
                for gene_data in genes:
                    gene_lower = gene_data['gene'].lower()
                    if any(crit in gene_lower for crit in [g.lower() for g in self.critical_amr_genes]):
                        critical_count += 1
                    if any(bac in gene_lower for bac in [g.lower() for g in self.bacmet2_markers]):
                        bacmet2_count += 1
                    if any(env_abx in gene_lower for env_abx in [g.lower() for g in self.common_antibiotic_markers]):
                        environmental_antibiotic_count += 1
                
                gene_centric['database_stats'][db_name] = {
                    'total_genes': total_genes,
                    'total_occurrences': total_occurrences,
                    'critical_genes': critical_count,
                    'bacmet2_genes': bacmet2_count,
                    'environmental_antibiotic_genes': environmental_antibiotic_count,
                    'coverage': f"{(total_occurrences / (total_samples * max(total_genes, 1))) * 100:.1f}%" if total_genes > 0 else "0%"
                }
        
        return gene_centric
    
    def _create_environmental_summary(self, gene_centric: Dict[str, Any], total_samples: int):
        """Create summary of environmental resistance markers"""
        env_categories = {
            'BACMET2 - Biocide Resistance': [],
            'BACMET2 - Metal Resistance': [],
            'Environmental Co-Selection': [],
            'Environmental Antibiotic Resistance': [],
            'Mobile Genetic Elements': [],
            'VICTORS Virulence Factors': []
        }
        
        # Categorize environmental genes
        for db_type in ['amr_databases', 'virulence_databases', 'environmental_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    gene = gene_data['gene']
                    gene_lower = gene.lower()
                    
                    # Check each category
                    # BACMET2 categories
                    if any(biocide in gene_lower for biocide in ['qac', 'cep', 'form', 'oqx']):
                        env_categories['BACMET2 - Biocide Resistance'].append(gene_data)
                    elif any(metal in gene_lower for metal in ['czc', 'mer', 'ars', 'cop', 'znt', 'chr', 'nik', 'cad', 'sil', 'pbr']):
                        env_categories['BACMET2 - Metal Resistance'].append(gene_data)
                    
                    # Environmental co-selection
                    elif any(co_sel in gene_lower for co_sel in ['sox', 'mar', 'rob', 'rpo', 'pho']):
                        env_categories['Environmental Co-Selection'].append(gene_data)
                    
                    # Mobile genetic elements
                    elif any(mobile in gene_lower for mobile in ['tra', 'mob', 'rep', 'ori', 'int', 'tnp', 'ist']):
                        env_categories['Mobile Genetic Elements'].append(gene_data)
                    
                    # Environmental antibiotic resistance
                    elif any(env_abx in gene_lower for env_abx in 
                            ['sul', 'dfr', 'cat', 'aac', 'aad', 'ant', 'aph', 'tet', 'erm', 'mef', 'msr']):
                        env_categories['Environmental Antibiotic Resistance'].append(gene_data)
                    
                    # VICTORS virulence
                    elif any(victor in gene_lower for victor in 
                            ['fim', 'pap', 'afa', 'dra', 'hly', 'cnf', 'sat', 'ast', 'stx', 
                             'iss', 'traT', 'kps', 'ibe', 'iut', 'iro', 'fyu', 'irp', 'chu',
                             'usp', 'vat', 'pic', 'sig', 'tia']):
                        env_categories['VICTORS Virulence Factors'].append(gene_data)
        
        # Sort each category by count
        for category, genes in env_categories.items():
            if genes:
                genes.sort(key=lambda x: x['count'], reverse=True)
                gene_centric['environmental_summary'][category] = {
                    'total_genes': len(genes),
                    'total_occurrences': sum(g['count'] for g in genes),
                    'genes': genes[:10000],  # Top 10000 genes per category
                    'percentage_of_samples': (sum(g['count'] for g in genes) / (total_samples * len(genes))) * 100 if genes else 0
                }
    
    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        """Create cross-genome patterns"""
        patterns = {
            'pasteur_st_distribution': Counter(),
            'oxford_st_distribution': Counter(),
            'international_clone_distribution': Counter(),
            'k_locus_distribution': Counter(),
            'o_locus_distribution': Counter(),
            'capsule_type_distribution': Counter(),
            'st_k_locus_combinations': defaultdict(list),
            'st_capsule_combinations': defaultdict(list),
            'gene_cooccurrence': defaultdict(Counter),
            'high_risk_combinations': [],
            'carbapenemase_patterns': defaultdict(list),
            'mdr_patterns': [],
            'database_coverage': {},
            'environmental_patterns': defaultdict(list)
        }
        
        samples_data = integrated_data.get('samples', {})
        gene_centric = integrated_data.get('gene_centric', {})
        
        # Collect all genes per sample
        sample_genes = defaultdict(set)
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'environmental_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    for genome in gene_data['genomes']:
                        sample_genes[genome].add(gene_data['gene'])
        
        # Analyze each sample
        for sample, data in samples_data.items():
            # MLST data
            pasteur_st = data.get('pasteur_mlst', {}).get('ST', 'ND')
            oxford_st = data.get('oxford_mlst', {}).get('ST', 'ND')
            ic = data.get('pasteur_mlst', {}).get('International_Clone', 'Unknown')
            
            # Kaptive data
            k_locus = data.get('kaptive', {}).get('K_Locus', 'ND')
            o_locus = data.get('kaptive', {}).get('O_Locus', 'ND')
            capsule_type = data.get('kaptive', {}).get('Capsule_Type', 'ND')
            
            # Distributions
            if pasteur_st != 'ND':
                patterns['pasteur_st_distribution'][pasteur_st] += 1
            if oxford_st != 'ND':
                patterns['oxford_st_distribution'][oxford_st] += 1
            if ic != 'Unknown':
                patterns['international_clone_distribution'][ic] += 1
            if k_locus != 'ND':
                patterns['k_locus_distribution'][k_locus] += 1
            if o_locus != 'ND':
                patterns['o_locus_distribution'][o_locus] += 1
            if capsule_type != 'ND':
                patterns['capsule_type_distribution'][capsule_type] += 1
            
            # Combinations
            if pasteur_st != 'ND' and k_locus != 'ND':
                patterns['st_k_locus_combinations'][f"ST{pasteur_st}-{k_locus}"].append(sample)
            if pasteur_st != 'ND' and capsule_type != 'ND':
                patterns['st_capsule_combinations'][f"ST{pasteur_st}-{capsule_type}"].append(sample)
            
            # Gene co-occurrence
            genes = list(sample_genes.get(sample, set()))
            for i, gene1 in enumerate(genes):
                for gene2 in genes[i+1:]:
                    patterns['gene_cooccurrence'][gene1][gene2] += 1
            
            # Check for critical resistance patterns
            carbapenemases = []
            esbls = []
            colistin_res = []
            tigecycline_res = []
            environmental_markers = []
            
            for gene in genes:
                gene_lower = gene.lower()
                
                if any(carb in gene_lower for carb in [g.lower() for g in self.critical_carbapenemases]):
                    carbapenemases.append(gene)
                elif any(esbl in gene_lower for esbl in [g.lower() for g in self.critical_esbls]):
                    esbls.append(gene)
                elif any(col in gene_lower for col in [g.lower() for g in self.critical_colistin]):
                    colistin_res.append(gene)
                elif any(tig in gene_lower for tig in [g.lower() for g in self.critical_tigecycline]):
                    tigecycline_res.append(gene)
                elif any(env in gene_lower for env in [g.lower() for g in self.all_environmental_markers]):
                    environmental_markers.append(gene)
            
            # Record patterns
            if carbapenemases:
                key = tuple(sorted(set(carbapenemases)))
                patterns['carbapenemase_patterns'][key].append(sample)
            
            # Environmental patterns
            if environmental_markers:
                env_key = tuple(sorted(set(environmental_markers)))
                patterns['environmental_patterns'][env_key].append(sample)
            
            # High-risk combinations (carbapenemase + last-resort)
            if carbapenemases and (colistin_res or tigecycline_res):
                patterns['high_risk_combinations'].append({
                    'sample': sample,
                    'pasteur_st': pasteur_st,
                    'international_clone': ic,
                    'k_locus': k_locus,
                    'capsule_type': capsule_type,
                    'carbapenemases': carbapenemases,
                    'colistin_resistance': colistin_res,
                    'tigecycline_resistance': tigecycline_res,
                    'environmental_markers': environmental_markers
                })
            
            # MDR patterns (3+ critical resistance types)
            resistance_types = 0
            if carbapenemases: resistance_types += 1
            if esbls: resistance_types += 1
            if colistin_res: resistance_types += 1
            if tigecycline_res: resistance_types += 1
            
            if resistance_types >= 3:
                patterns['mdr_patterns'].append({
                    'sample': sample,
                    'pasteur_st': pasteur_st,
                    'international_clone': ic,
                    'resistance_types': resistance_types,
                    'carbapenemases': carbapenemases,
                    'esbls': esbls,
                    'colistin_resistance': colistin_res,
                    'tigecycline_resistance': tigecycline_res,
                    'environmental_markers': environmental_markers
                })
        
        # Calculate percentages for distributions
        for dist_name in ['pasteur_st_distribution', 'oxford_st_distribution', 'k_locus_distribution', 
                         'o_locus_distribution', 'capsule_type_distribution']:
            if patterns[dist_name]:
                total = sum(patterns[dist_name].values())
                for key in patterns[dist_name]:
                    patterns[dist_name][key] = {
                        'count': patterns[dist_name][key],
                        'percentage': (patterns[dist_name][key] / total) * 100,
                        'frequency_display': f"{patterns[dist_name][key]} ({(patterns[dist_name][key] / total) * 100:.1f}%)"
                    }
        
        # Calculate database coverage
        all_samples = set(samples_data.keys())
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'environmental_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                samples_with_hits = set()
                for gene_data in genes:
                    samples_with_hits.update(gene_data['genomes'])
                
                coverage = len(samples_with_hits) / len(all_samples) * 100 if all_samples else 0
                patterns['database_coverage'][db_name] = {
                    'samples_with_hits': len(samples_with_hits),
                    'total_samples': len(all_samples),
                    'coverage_percentage': round(coverage, 2),
                    'coverage_display': f"{len(samples_with_hits)} ({coverage:.1f}%)"
                }
        
        return patterns

    def create_plasmid_analysis(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        """Create comprehensive plasmid analysis tables"""
        plasmid_analysis = {
            'plasmid_databases': {},
            'plasmid_frequencies': [],
            'plasmid_categories': defaultdict(list),
            'plasmid_cooccurrence': defaultdict(Counter),
            'sample_plasmid_profiles': defaultdict(list),
            'plasmid_summary_stats': {},
            'high_frequency_plasmids': [],
            'unique_plasmid_patterns': defaultdict(list)
        }
        
        # Check if we have plasmid data - FIXED to look in the correct location
        if 'plasmidfinder' not in integrated_data.get('gene_frequencies', {}):
            return plasmid_analysis
        
        # Process PlasmidFinder data
        plasmid_genes = integrated_data['gene_frequencies']['plasmidfinder']
        gene_list = []
        
        for gene, data in plasmid_genes.items():
            # Categorize plasmid
            plasmid_type = 'Unknown'
            gene_lower = gene.lower()
            
            if any(inc in gene_lower for inc in ['inc']):
                inc_match = re.search(r'inc([a-z0-9]+)', gene_lower)
                if inc_match:
                    plasmid_type = f'Inc{inc_match.group(1).upper()}'
                else:
                    plasmid_type = 'Incompatibility group'
            elif 'col' in gene_lower:
                plasmid_type = 'Colicin plasmid'
            elif any(rep in gene_lower for rep in ['rep', 'rep_']):
                plasmid_type = 'Replication protein'
            elif any(mob in gene_lower for mob in ['mob', 'tra', 'conj']):
                plasmid_type = 'Mobility gene'
            else:
                plasmid_type = 'Other plasmid marker'
            
            gene_list.append({
                'plasmid_marker': gene,
                'full_name': data.get('full_name', gene),
                'category': plasmid_type,
                'database': 'PlasmidFinder',
                'count': data.get('count', 0),
                'percentage': data.get('percentage', 0),
                'frequency_display': data.get('frequency_display', f"{data.get('count', 0)} ({data.get('percentage', 0):.1f}%)"),
                'genomes': data.get('genomes', [])
            })
        
        if gene_list:
            # Sort by count
            gene_list.sort(key=lambda x: x['count'], reverse=True)
            plasmid_analysis['plasmid_databases']['plasmidfinder'] = gene_list
            
            # Add to plasmid categories
            for gene_data in gene_list:
                plasmid_analysis['plasmid_categories'][gene_data['category']].append(gene_data)
            
            # Create plasmid frequencies list
            plasmid_analysis['plasmid_frequencies'] = gene_list
            
            # Find high frequency plasmids (present in >30% of samples)
            high_freq_threshold = total_samples * 0.3
            plasmid_analysis['high_frequency_plasmids'] = [
                p for p in gene_list if p['count'] >= high_freq_threshold
            ]
        
        # Rest of the method continues as before...
        # Create sample plasmid profiles
        samples_data = integrated_data.get('samples', {})
        for sample, data in samples_data.items():
            # Get plasmids from all sources
            sample_plasmids = []
            
            # Check if we have plasmid data in the sample
            if 'plasmidfinder' in integrated_data.get('gene_frequencies', {}):
                db_genes = integrated_data['gene_frequencies']['plasmidfinder']
                for gene, gene_data in db_genes.items():
                    if sample in gene_data.get('genomes', []):
                        sample_plasmids.append({
                            'marker': gene,
                            'database': 'plasmidfinder',
                            'category': self._categorize_plasmid_marker(gene)
                        })
            
            if sample_plasmids:
                plasmid_analysis['sample_plasmid_profiles'][sample] = sample_plasmids
                
                # Track unique plasmid patterns
                pattern_key = tuple(sorted([p['marker'] for p in sample_plasmids]))
                plasmid_analysis['unique_plasmid_patterns'][pattern_key].append(sample)
        
        # Calculate plasmid co-occurrence
        for sample, plasmids in plasmid_analysis['sample_plasmid_profiles'].items():
            plasmid_names = [p['marker'] for p in plasmids]
            for i, p1 in enumerate(plasmid_names):
                for p2 in plasmid_names[i+1:]:
                    plasmid_analysis['plasmid_cooccurrence'][p1][p2] += 1
        
        # Calculate summary statistics
        total_plasmid_markers = sum(len(db) for db in plasmid_analysis['plasmid_databases'].values())
        total_plasmid_occurrences = sum(sum(p['count'] for p in db) for db in plasmid_analysis['plasmid_databases'].values())
        samples_with_plasmids = len(plasmid_analysis['sample_plasmid_profiles'])
        
        plasmid_analysis['plasmid_summary_stats'] = {
            'total_plasmid_markers': total_plasmid_markers,
            'total_plasmid_occurrences': total_plasmid_occurrences,
            'samples_with_plasmids': samples_with_plasmids,
            'total_samples': total_samples,
            'plasmid_prevalence': (samples_with_plasmids / total_samples * 100) if total_samples > 0 else 0,
            'unique_plasmid_patterns': len(plasmid_analysis['unique_plasmid_patterns'])
        }
        
        return plasmid_analysis

    def _categorize_plasmid_marker(self, marker: str) -> str:
        """Categorize individual plasmid marker"""
        marker_lower = marker.lower()
        
        if any(inc in marker_lower for inc in ['inc']):
            return 'Incompatibility group'
        elif 'col' in marker_lower:
            return 'Colicin plasmid'
        elif any(rep in marker_lower for rep in ['rep', 'rep_']):
            return 'Replication protein'
        elif any(mob in marker_lower for mob in ['mob', 'tra', 'conj']):
            return 'Mobility gene'
        elif 'ecoh' in marker_lower or 'e.coli' in marker_lower:
            return 'E. coli specific'
        else:
            return 'Other plasmid marker'


class UltimateHTMLGenerator:
    """Generates ultimate HTML reports for A. baumannii"""
    
    def __init__(self, data_analyzer: UltimateDataAnalyzer):
        self.data_analyzer = data_analyzer
        self.tab_colors = {
            'summary': "#5FAE62",
            'samples': '#2196F3',
            'mlst': '#FF9800',
            'kaptive': '#9C27B0',
            'amr': '#F44336',
            'virulence': '#E91E63',
            'environmental': '#795548',
            'categories': '#009688',
            'patterns': '#FF5722',
            'export': '#3F51B5',
            'databases': '#607D8B',
            'plasmid': '#2196F3'
        }
    
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        """Generate the ultimate HTML report"""
        print("\n🎨 Generating ULTIMATE HTML report for A. baumannii...")
        
        # Extract data
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        plasmid_analysis = integrated_data.get('plasmid_analysis', {})
        
        # Create HTML
        html = self._create_ultimate_html(
            metadata=metadata,
            samples_data=samples_data,
            patterns=patterns,
            gene_centric=gene_centric,
            plasmid_analysis=plasmid_analysis,
            total_samples=len(samples_data)  
    
        )
        
        # Save HTML file
        output_file = output_dir / "genius_acinetobacter_ultimate_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)
    
    def _create_ultimate_html(self, **kwargs) -> str:
        """Create ultimate HTML with all sections - FIXED FOR PLASMID TAB"""
        
        # CSS Styles - IMPROVED FOR NO TRUNCATION AND SCROLLABILITY
        css = """
        <style>
        :root {
            --summary-color: #4CAF50;
            --samples-color: #2196F3;
            --mlst-color: #FF9800;
            --kaptive-color: #9C27B0;
            --amr-color: #F44336;
            --virulence-color: #E91E63;
            --environmental-color: #795548;
            --categories-color: #009688;
            --patterns-color: #FF5722;
            --export-color: #3F51B5;
            --databases-color: #607D8B;
            --plasmid-color: #2196F3;
        }
        
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
            overflow-x: auto;
        }
        
        .container {
            width: 100%;
            max-width: 100%;
            margin: 0 auto;
            padding: 20px;
            overflow-x: hidden;
        }
        
        .main-header {
            background: linear-gradient(135deg, #00695c 0%, #004d40 100%);
            color: white;
            padding: 30px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            margin-bottom: 30px;
            text-align: center;
        }
        
        .main-header h1 {
            font-size: 2.8em;
            margin-bottom: 10px;
            color: white;
        }
        
        .metadata-bar {
            background: rgba(255,255,255,0.1);
            padding: 15px;
            border-radius: 10px;
            margin: 20px 0;
            display: flex;
            justify-content: space-around;
            flex-wrap: wrap;
            gap: 15px;
            backdrop-filter: blur(10px);
        }
        
        .metadata-item {
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 0.95em;
        }
        
        .dashboard-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        
        .dashboard-card {
            background: white;
            padding: 25px;
            border-radius: 12px;
            box-shadow: 0 5px 20px rgba(0,0,0,0.1);
            text-align: center;
            transition: all 0.3s ease;
            cursor: pointer;
            border-left: 5px solid;
            position: relative;
            overflow: hidden;
        }
        
        .dashboard-card:hover {
            transform: translateY(-10px);
            box-shadow: 0 15px 30px rgba(0,0,0,0.2);
        }
        
        .card-summary { border-left-color: var(--summary-color); }
        .card-samples { border-left-color: var(--samples-color); }
        .card-mlst { border-left-color: var(--mlst-color); }
        .card-kaptive { border-left-color: var(--kaptive-color); }
        .card-amr { border-left-color: var(--amr-color); }
        .card-virulence { border-left-color: var(--virulence-color); }
        .card-environmental { border-left-color: var(--environmental-color); }
        .card-categories { border-left-color: var(--categories-color); }
        .card-patterns { border-left-color: var(--patterns-color); }
        .card-export { border-left-color: var(--export-color); }
        .card-databases { border-left-color: var(--databases-color); }
        
        .card-number {
            font-size: 3em;
            font-weight: bold;
            margin: 15px 0;
            background: linear-gradient(90deg, #00695c, #004d40);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        
        .tab-navigation {
            display: flex;
            gap: 5px;
            margin-bottom: 20px;
            flex-wrap: wrap;
            background: white;
            padding: 15px;
            border-radius: 12px;
            box-shadow: 0 5px 20px rgba(0,0,0,0.1);
            position: sticky;
            top: 10px;
            z-index: 100;
        }
        
        .tab-button {
            padding: 12px 25px;
            background: #f5f5f5;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-weight: 600;
            color: #666;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            gap: 8px;
            position: relative;
            overflow: hidden;
        }
        
        .tab-button::after {
            content: '';
            position: absolute;
            bottom: 0;
            left: 50%;
            right: 50%;
            height: 3px;
            background: currentColor;
            transition: all 0.3s ease;
        }
        
        .tab-button:hover::after {
            left: 10%;
            right: 10%;
        }
        
        .tab-button.active {
            color: white;
        }
        
        .tab-button.active::after {
            left: 10%;
            right: 10%;
        }
        
        .tab-button.summary.active { background: var(--summary-color); }
        .tab-button.samples.active { background: var(--samples-color); }
        .tab-button.mlst.active { background: var(--mlst-color); }
        .tab-button.kaptive.active { background: var(--kaptive-color); }
        .tab-button.amr.active { background: var(--amr-color); }
        .tab-button.virulence.active { background: var(--virulence-color); }
        .tab-button.environmental.active { background: var(--environmental-color); }
        .tab-button.categories.active { background: var(--categories-color); }
        .tab-button.patterns.active { background: var(--patterns-color); }
        .tab-button.export.active { background: var(--export-color); }
        .tab-button.databases.active { background: var(--databases-color); }
        .tab-button.plasmid.active { background: #2196F3; }
        
        .tab-content {
            display: none;
            background: white;
            padding: 30px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            margin-bottom: 30px;
            animation: fadeIn 0.5s ease;
            width: 100%;
            overflow: hidden;
        }
        
        .tab-content.active {
            display: block;
        }
        
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
        }
        
        .section-header {
            color: #2c3e50;
            margin-bottom: 25px;
            padding-bottom: 15px;
            border-bottom: 3px solid;
            font-size: 1.8em;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
        .summary-header { border-color: var(--summary-color); }
        .samples-header { border-color: var(--samples-color); }
        .mlst-header { border-color: var(--mlst-color); }
        .kaptive-header { border-color: var(--kaptive-color); }
        .amr-header { border-color: var(--amr-color); }
        .virulence-header { border-color: var(--virulence-color); }
        .environmental-header { border-color: var(--environmental-color); }
        .categories-header { border-color: var(--categories-color); }
        .patterns-header { border-color: var(--patterns-color); }
        .export-header { border-color: var(--export-color); }
        .databases-header { border-color: var(--databases-color); }
        
        /* FIXED TABLE STYLES - FULL SCROLLABILITY */
        .data-table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.95em;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
        }
        
        .data-table th {
            background: #2c3e50;
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
            position: sticky;
            top: 0;
            white-space: nowrap;
            z-index: 10;
        }
        
        .data-table td {
            padding: 12px 15px;
            border-bottom: 1px solid #e0e0e0;
            vertical-align: top;
            word-wrap: break-word;
            word-break: break-word;
            white-space: normal;
        }
        
        .data-table tr:hover {
            background: #f8f9fa;
        }
        
        /* MASTER SCROLLABLE CONTAINER - NO TRUNCATION */
        .master-scrollable-container {
            width: 100%;
            max-width: 100%;
            overflow-x: auto;
            overflow-y: visible;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            margin: 20px 0;
            position: relative;
        }
        
        .master-scrollable-container table {
            min-width: 100%;
            width: auto;
        }
        
        .master-scrollable-container table.data-table {
            margin: 0;
        }
        
        .search-box {
            width: 100%;
            max-width: 100%;
            padding: 12px;
            margin-bottom: 20px;
            border: 2px solid #e0e0e0;
            border-radius: 8px;
            font-size: 1em;
            transition: all 0.3s ease;
        }
        
        .search-box:focus {
            outline: none;
            border-color: #00695c;
            box-shadow: 0 0 0 3px rgba(0, 105, 92, 0.1);
        }
        
        .badge {
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
            margin: 2px;
            white-space: nowrap;
        }
        
        .badge-low { background: #4CAF50; color: white; }
        .badge-medium { background: #FF9800; color: black; }
        .badge-high { background: #F44336; color: white; }
        .badge-critical { background: #9C27B0; color: white; }
        
        .alert-box {
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
            display: flex;
            align-items: center;
            gap: 20px;
            border-left: 5px solid;
        }
        
        .alert-success { background: #d4edda; color: #155724; border-left-color: #28a745; }
        .alert-warning { background: #fff3cd; color: #856404; border-left-color: #ffc107; }
        .alert-danger { background: #f8d7da; color: #721c24; border-left-color: #dc3545; }
        .alert-info { background: #d1ecf1; color: #0c5460; border-left-color: #17a2b8; }
        
        .action-buttons {
            display: flex;
            gap: 10px;
            margin: 20px 0;
            flex-wrap: wrap;
        }
        
        .action-btn {
            padding: 10px 20px;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-weight: 600;
            display: flex;
            align-items: center;
            gap: 8px;
            transition: all 0.3s ease;
            white-space: nowrap;
        }
        
        .action-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 5px 15px rgba(0,0,0,0.2);
        }
        
        .btn-primary { background: #00695c; color: white; }
        .btn-success { background: #28a745; color: white; }
        .btn-danger { background: #dc3545; color: white; }
        .btn-warning { background: #ffc107; color: black; }
        
        .database-section {
            margin: 30px 0;
            padding: 25px;
            border-radius: 12px;
            background: #f8f9fa;
            box-shadow: 0 3px 15px rgba(0,0,0,0.08);
            overflow: hidden;
        }
        
        .database-header {
            font-size: 1.4em;
            color: #2c3e50;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #00695c;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
        .print-section-btn {
            background: #00695c;
            color: white;
            border: none;
            border-radius: 5px;
            padding: 8px 15px;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 5px;
            font-size: 0.9em;
            white-space: nowrap;
        }
        
        .print-section-btn:hover {
            background: #004d40;
        }
        
        /* GENOME LIST - FULL DISPLAY */
        .genome-list {
            display: block;
            max-height: 200px;
            overflow-y: auto;
            padding: 5px;
            background: #f8f9fa;
            border-radius: 5px;
            border: 1px solid #e0e0e0;
        }
        
        .genome-tag {
            display: inline-block;
            background: #e0f2f1;
            color: #00695c;
            padding: 3px 10px;
            border-radius: 12px;
            font-size: 0.85em;
            border: 1px solid #b2dfdb;
            margin: 2px;
            word-break: break-all;
            white-space: normal;
        }
        
        .footer {
            text-align: center;
            padding: 30px;
            color: white;
            margin-top: 40px;
            border-radius: 15px;
            background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
        }
        
        .category-chip {
            display: inline-block;
            padding: 4px 12px;
            border-radius: 15px;
            font-size: 0.8em;
            font-weight: 600;
            margin: 2px;
            white-space: nowrap;
        }
        
        .chip-carbapenemase { background: #f8d7da; color: #721c24; border: 1px solid #f5c6cb; }
        .chip-esbl { background: #fff3cd; color: #856404; border: 1px solid #ffeaa7; }
        .chip-ampc { background: #d1ecf1; color: #0c5460; border: 1px solid #bee5eb; }
        .chip-colistin { background: #d4edda; color: #155724; border: 1px solid #c3e6cb; }
        .chip-tigecycline { background: #e2e3e5; color: #383d41; border: 1px solid #d6d8db; }
        .chip-biofilm { background: #f8f9fa; color: #343a40; border: 1px solid #e9ecef; }
        .chip-efflux { background: #e8eaf6; color: #303f9f; border: 1px solid #c5cae9; }
        .chip-environmental { background: #e8f5e8; color: #1b5e20; border: 1px solid #c8e6c9; }
        .chip-bacmet2 { background: #f3e5f5; color: #4a148c; border: 1px solid #e1bee7; }
        .chip-victors { background: #fce4ec; color: #880e4f; border: 1px solid #f8bbd0; }
        .chip-other { background: #f5f5f5; color: #212121; border: 1px solid #e0e0e0; }
        
        /* Frequency display styling */
        .frequency-display {
            font-weight: 600;
            color: #2c3e50;
        }
        
        /* Responsive design */
        @media (max-width: 768px) {
            .container {
                padding: 10px;
            }
            
            .main-header {
                padding: 20px;
            }
            
            .main-header h1 {
                font-size: 2em;
            }
            
            .tab-button {
                padding: 10px 15px;
                font-size: 0.9em;
            }
            
            .dashboard-grid {
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            }
            
            .data-table {
                font-size: 0.85em;
            }
            
            .metadata-bar {
                flex-direction: column;
                gap: 10px;
            }
        }
        
        /* Column specific styles for better readability */
        .col-gene { min-width: 200px; }
        .col-category { min-width: 180px; }
        .col-frequency { min-width: 120px; }
        .col-sample { min-width: 250px; }
        .col-genomes { min-width: 400px; }
        .col-database { min-width: 120px; }
        .col-risk { min-width: 100px; }
        .col-st { min-width: 100px; }
        .col-ic { min-width: 150px; }
        .col-k-locus { min-width: 100px; }
        .col-o-locus { min-width: 100px; }
        .col-capsule { min-width: 150px; }
        </style>
        """
        
        # JavaScript - Fixed version without backslash issues in f-string
        js = """
        <script>
        // Tab switching
        function switchTab(tabName) {
            // Hide all tabs
            document.querySelectorAll('.tab-content').forEach(tab => {
                tab.classList.remove('active');
            });
            
            // Remove active class from all buttons
            document.querySelectorAll('.tab-button').forEach(button => {
                button.classList.remove('active');
            });
            
            // Show selected tab
            document.getElementById(tabName + '-tab').classList.add('active');
            
            // Activate selected button
            event.currentTarget.classList.add('active');
            
            // Update URL hash
            window.location.hash = tabName;
        }
        
        // Search functionality
        function searchTable(tableId, searchId) {
            const input = document.getElementById(searchId);
            const filter = input.value.toUpperCase();
            const table = document.getElementById(tableId);
            const rows = table.getElementsByTagName('tr');
            
            for (let i = 1; i < rows.length; i++) {
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                
                for (let j = 0; j < cells.length; j++) {
                    const cell = cells[j];
                    if (cell) {
                        const txtValue = cell.textContent || cell.innerText;
                        if (txtValue.toUpperCase().indexOf(filter) > -1) {
                            found = true;
                            break;
                        }
                    }
                }
                
                rows[i].style.display = found ? '' : 'none';
            }
        }
        
        // Print current section
        function printSection(sectionId) {
            const content = document.getElementById(sectionId);
            const printWindow = window.open('', '_blank');
            printWindow.document.write('<html><head><title>Print Section</title>');
            printWindow.document.write('<style>' + document.querySelector('style').textContent + '</style>');
            printWindow.document.write('</head><body>');
            printWindow.document.write(content.innerHTML);
            printWindow.document.write('</body></html>');
            printWindow.document.close();
            printWindow.print();
        }
        
        // Export table to CSV
        function exportTableToCSV(tableId, filename) {
            const table = document.getElementById(tableId);
            const rows = table.querySelectorAll('tr');
            const csv = [];
            
            for (let i = 0; i < rows.length; i++) {
                const row = [], cols = rows[i].querySelectorAll('td, th');
                
                for (let j = 0; j < cols.length; j++) {
                    row.push('"' + (cols[j].innerText || '').replace(/"/g, '""') + '"');
                }
                
                csv.push(row.join(','));
            }
            
            // Use template literal with String.raw to avoid backslash issues
            const csvContent = csv.join(String.raw`\n`);
            const csvFile = new Blob([csvContent], {type: 'text/csv'});
            const downloadLink = document.createElement('a');
            downloadLink.download = filename;
            downloadLink.href = window.URL.createObjectURL(csvFile);
            downloadLink.style.display = 'none';
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
        }
        
        // Initialize from URL hash
        document.addEventListener('DOMContentLoaded', function() {
            const hash = window.location.hash.substring(1);
            if (hash) {
                const tabButton = document.querySelector(`.tab-button.${hash}`);
                if (tabButton) {
                    tabButton.click();
                }
            } else {
                // Show first tab
                document.querySelector('.tab-button').click();
            }
            
            // Add column classes for better display
            const tables = document.querySelectorAll('.data-table');
            tables.forEach(table => {
                const headers = table.querySelectorAll('th');
                headers.forEach((header, index) => {
                    const headerText = header.textContent.toLowerCase();
                    if (headerText.includes('gene')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-gene');
                        });
                    } else if (headerText.includes('category')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-category');
                        });
                    } else if (headerText.includes('frequency')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-frequency');
                        });
                    } else if (headerText.includes('sample')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-sample');
                        });
                    } else if (headerText.includes('genome')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-genomes');
                        });
                    } else if (headerText.includes('database')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-database');
                        });
                    } else if (headerText.includes('risk')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-risk');
                        });
                    } else if (headerText.includes('st') && !headerText.includes('sample')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-st');
                        });
                    } else if (headerText.includes('international') || headerText.includes('ic')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-ic');
                        });
                    } else if (headerText.includes('k locus')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-k-locus');
                        });
                    } else if (headerText.includes('o locus')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-o-locus');
                        });
                    } else if (headerText.includes('capsule')) {
                        table.querySelectorAll(`tr td:nth-child(${index + 1})`).forEach(cell => {
                            cell.classList.add('col-capsule');
                        });
                    }
                });
            });
        });
        </script>
        """
        
        # Extract data
        metadata = kwargs['metadata']
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        plasmid_analysis = kwargs.get('plasmid_analysis', {}) 
        
        total_samples = len(samples_data)
        amr_databases = gene_centric.get('amr_databases', {})
        virulence_databases = gene_centric.get('virulence_databases', {})
        environmental_databases = gene_centric.get('environmental_databases', {})
        environmental_summary = gene_centric.get('environmental_summary', {})
        
        # Calculate statistics
        total_amr_genes = sum(len(db) for db in amr_databases.values())
        total_virulence_genes = sum(len(db) for db in virulence_databases.values())
        total_environmental_genes = sum(len(db) for db in environmental_databases.values())
        high_risk_count = len(patterns.get('high_risk_combinations', []))
        mdr_count = len(patterns.get('mdr_patterns', []))
        
        # Check if we have plasmid data
        has_plasmid_data = bool(plasmid_analysis.get('plasmid_databases', {}))
        if has_plasmid_data:
            plasmid_databases = plasmid_analysis.get('plasmid_databases', {})
            total_plasmid_genes = sum(len(db) for db in plasmid_databases.values())
        else:
            total_plasmid_genes = 0
        
        # Count carbapenemases and environmental markers
        carbapenemase_count = 0
        environmental_marker_count = 0
        for db_name, genes in amr_databases.items():
            for gene_data in genes:
                if gene_data['category'] == 'Carbapenemases':
                    carbapenemase_count += 1
                if gene_data['category'] in ['Environmental Co-Selection', 'BACMET2 Markers']:
                    environmental_marker_count += 1
        
        # Environmental database specific counts
        if 'bacmet2' in environmental_databases:
            environmental_marker_count += len(environmental_databases['bacmet2'])
        
        # Build HTML - use string concatenation instead of putting js in f-string
        html_parts = []
        
        # Start building HTML
        html_parts.append(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GENIUS Acinetobacter baumannii Ultimate Report</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    {css}
""")
        
        # Add JavaScript directly (not in f-string)
        html_parts.append(js)
        
        # Continue with the rest of the HTML
        html_parts.append(f"""
</head>
<body>
    <div class="container">
        <!-- Main Header -->
        <div class="main-header">
            <h1><i class="fas fa-bacterium"></i> GENIUS Acinetobacter baumannii Ultimate Analysis Report</h1>
            <p>Comprehensive Gene-Centric Cross-Genome Analysis with Environmental Markers</p>
            
            <div class="metadata-bar">
                <div class="metadata-item">
                    <i class="fas fa-calendar"></i>
                    <span>Generated: {metadata.get('analysis_date', 'Unknown')}</span>
                </div>
                <div class="metadata-item">
                    <i class="fas fa-database"></i>
                    <span>Samples: {total_samples}</span>
                </div>
                <div class="metadata-item">
                    <i class="fas fa-vial"></i>
                    <span>Pathogen: Acinetobacter baumannii</span>
                </div>
                <div class="metadata-item">
                    <i class="fas fa-university"></i>
                    <span>University of Ghana Medical School</span>
                </div>
            </div>
        </div>
        
        <!-- Dashboard -->
        <div class="dashboard-grid">
            <div class="dashboard-card card-summary" onclick="switchTab('summary')">
                <div class="card-number">{total_samples}</div>
                <div class="card-label">Total Samples</div>
                <i class="fas fa-vial fa-2x" style="color: var(--summary-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-mlst" onclick="switchTab('mlst')">
                <div class="card-number">{len(patterns.get('pasteur_st_distribution', {}))}</div>
                <div class="card-label">Pasteur STs</div>
                <i class="fas fa-code-branch fa-2x" style="color: var(--mlst-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-kaptive" onclick="switchTab('kaptive')">
                <div class="card-number">{len(patterns.get('k_locus_distribution', {}))}</div>
                <div class="card-label">Capsule Types</div>
                <i class="fas fa-shield-alt fa-2x" style="color: var(--kaptive-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-amr" onclick="switchTab('amr')">
                <div class="card-number">{total_amr_genes}</div>
                <div class="card-label">AMR Genes</div>
                <i class="fas fa-biohazard fa-2x" style="color: var(--amr-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-virulence" onclick="switchTab('virulence')">
                <div class="card-number">{total_virulence_genes}</div>
                <div class="card-label">Virulence Genes</div>
                <i class="fas fa-virus fa-2x" style="color: var(--virulence-color); margin-top: 10px;"></i>
            </div>
            
            <div class="dashboard-card card-environmental" onclick="switchTab('environmental')">
                <div class="card-number">{total_environmental_genes}</div>
                <div class="card-label">Environmental</div>
                <i class="fas fa-globe-africa fa-2x" style="color: var(--environmental-color); margin-top: 10px;"></i>
            </div>
            """)
        
        # Add plasmid dashboard card only if we have plasmid data
        if has_plasmid_data:
            html_parts.append(f"""
            <div class="dashboard-card card-plasmid" onclick="switchTab('plasmid')" style="border-left: 5px solid #2196F3;">
                <div class="card-number">{total_plasmid_genes}</div>
                <div class="card-label">Plasmid Markers</div>
                <i class="fas fa-dna fa-2x" style="color: #2196F3; margin-top: 10px;"></i>
            </div>
            """)
        
        html_parts.append(f"""
            <div class="dashboard-card card-patterns" onclick="switchTab('patterns')">
                <div class="card-number">{high_risk_count}</div>
                <div class="card-label">High-Risk</div>
                <i class="fas fa-exclamation-triangle fa-2x" style="color: var(--patterns-color); margin-top: 10px;"></i>
            </div>
        </div>
        
        <!-- Tab Navigation -->
        <div class="tab-navigation">
            <button class="tab-button summary active" onclick="switchTab('summary')">
                <i class="fas fa-chart-pie"></i> Summary
            </button>
            <button class="tab-button samples" onclick="switchTab('samples')">
                <i class="fas fa-list-alt"></i> Sample Overview
            </button>
            <button class="tab-button mlst" onclick="switchTab('mlst')">
                <i class="fas fa-code-branch"></i> MLST Analysis
            </button>
            <button class="tab-button kaptive" onclick="switchTab('kaptive')">
                <i class="fas fa-shield-alt"></i> Capsule Typing
            </button>
            <button class="tab-button amr" onclick="switchTab('amr')">
                <i class="fas fa-biohazard"></i> AMR Genes
            </button>
            <button class="tab-button virulence" onclick="switchTab('virulence')">
                <i class="fas fa-virus"></i> Virulence Genes
            </button>
            <button class="tab-button environmental" onclick="switchTab('environmental')">
                <i class="fas fa-globe-africa"></i> Environmental
            </button>
            <button class="tab-button categories" onclick="switchTab('categories')">
                <i class="fas fa-tags"></i> Gene Categories
            </button>
            <button class="tab-button patterns" onclick="switchTab('patterns')">
                <i class="fas fa-project-diagram"></i> Pattern Discovery
            </button>
            <button class="tab-button databases" onclick="switchTab('databases')">
                <i class="fas fa-database"></i> Database Coverage
            </button>
            """)
        
        # Add plasmid button only if we have plasmid data
        if has_plasmid_data:
            html_parts.append("""<button class="tab-button plasmid" onclick="switchTab('plasmid')"><i class="fas fa-dna"></i> Plasmids</button>""")
        
        html_parts.append(f"""
            <button class="tab-button export" onclick="switchTab('export')">
                <i class="fas fa-download"></i> Export
            </button>
        </div>
        
        <!-- Summary Tab -->
        <div id="summary-tab" class="tab-content active">
            {self._generate_summary_section(kwargs)}
        </div>
        
        <!-- Sample Overview Tab -->
        <div id="samples-tab" class="tab-content">
            {self._generate_sample_overview_section(kwargs)}
        </div>
        
        <!-- MLST Analysis Tab -->
        <div id="mlst-tab" class="tab-content">
            {self._generate_mlst_section(kwargs)}
        </div>
        
        <!-- Kaptive Analysis Tab -->
        <div id="kaptive-tab" class="tab-content">
            {self._generate_kaptive_section(kwargs)}
        </div>
        
        <!-- AMR Genes Tab -->
        <div id="amr-tab" class="tab-content">
            {self._generate_amr_section(kwargs)}
        </div>
        
        <!-- Virulence Genes Tab -->
        <div id="virulence-tab" class="tab-content">
            {self._generate_virulence_section(kwargs)}
        </div>
        
        <!-- Environmental Tab -->
        <div id="environmental-tab" class="tab-content">
            {self._generate_environmental_section(kwargs)}
        </div>
        
        <!-- Gene Categories Tab -->
        <div id="categories-tab" class="tab-content">
            {self._generate_categories_section(kwargs)}
        </div>
        
        <!-- Pattern Discovery Tab -->
        <div id="patterns-tab" class="tab-content">
            {self._generate_pattern_discovery_section(kwargs)}
        </div>
        
        <!-- Database Coverage Tab -->
        <div id="databases-tab" class="tab-content">
            {self._generate_database_coverage_section(kwargs)}
        </div>
        """)
        
        # Add Plasmid Tab only if we have plasmid data
        if has_plasmid_data:
            html_parts.append(f"""
        <!-- Plasmid Analysis Tab -->
        <div id="plasmid-tab" class="tab-content">
            {self._generate_plasmid_section({
                'plasmid_analysis': plasmid_analysis,
                'total_samples': total_samples,
                'samples_data': samples_data,
                'patterns': patterns,
                'gene_centric': gene_centric,
                'metadata': metadata
            })}
        </div>
            """)
        
        html_parts.append(f"""
        <!-- Export Tab -->
        <div id="export-tab" class="tab-content">
            {self._generate_export_section(kwargs)}
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <h3>GENIUS Acinetobacter baumannii Ultimate Reporter v1.0.0</h3>
            <p>University of Ghana Medical School | Brown Beckley <brownbeckley94@gmail.com></p>
            <p>Generated on {metadata.get('analysis_date', 'Unknown')}</p>
            <p><strong>Critical Genes Tracked:</strong> Carbapenemases • ESBLs • Colistin Resistance • Tigecycline Resistance • Biofilm Formation • Efflux Pumps • Environmental Co-Selection</p>
        </div>
    </div>
</body>
</html>
        """)
        
        # Join all parts
        return ''.join(html_parts)
    
    def _generate_summary_section(self, kwargs: Dict) -> str:
        """Generate summary section"""
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        
        total_samples = len(samples_data)
        total_amr_genes = sum(len(db) for db in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(db) for db in gene_centric.get('virulence_databases', {}).values())
        total_environmental_genes = sum(len(db) for db in gene_centric.get('environmental_databases', {}).values())
        high_risk_count = len(patterns.get('high_risk_combinations', []))
        mdr_count = len(patterns.get('mdr_patterns', []))
        
        # Count carbapenemases and environmental markers
        carbapenemase_count = 0
        environmental_marker_count = 0
        for db_name, genes in gene_centric.get('amr_databases', {}).items():
            for gene_data in genes:
                if gene_data['category'] == 'Carbapenemases':
                    carbapenemase_count += 1
                if gene_data['category'] in ['Environmental Co-Selection', 'BACMET2 Markers']:
                    environmental_marker_count += 1
        
        if 'bacmet2' in gene_centric.get('environmental_databases', {}):
            environmental_marker_count += len(gene_centric['environmental_databases']['bacmet2'])
        
        html = f"""
        <h2 class="section-header summary-header">
            <i class="fas fa-chart-pie"></i> Executive Summary - A. baumannii MDR Analysis
            <button class="print-section-btn" onclick="printSection('summary-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>MDR/XDR A. baumannii Analysis with Environmental Co-Selection</h3>
                <p>This comprehensive gene-centric report analyzes <strong>{total_samples}</strong> A. baumannii genomes 
                with focus on <strong>carbapenem resistance, colistin resistance, biofilm formation, and environmental co-selection markers</strong>. 
                Each gene is shown with ALL genomes that contain it for easy tracking of resistance spread.</p>
            </div>
        </div>
        
        {f'<div class="alert-box alert-danger"><i class="fas fa-exclamation-triangle fa-2x"></i><div><h3>⚠️ Critical Resistance Alert</h3><p><strong>{carbapenemase_count} carbapenemase genes</strong> detected across samples. Carbapenem-resistant A. baumannii (CRAB) is a WHO Critical Priority pathogen.</p></div></div>' if carbapenemase_count > 0 else ''}
        
        {f'<div class="alert-box alert-warning"><i class="fas fa-globe-africa fa-2x"></i><div><h3>⚠️ Environmental Co-Selection Alert</h3><p><strong>{environmental_marker_count} environmental resistance markers</strong> detected. These genes can co-select for antibiotic resistance in hospital environments.</p></div></div>' if environmental_marker_count > 0 else ''}
        
        {f'<div class="alert-box alert-warning"><i class="fas fa-radiation fa-2x"></i><div><h3>⚠️ High-Risk Combinations Detected</h3><p><strong>{high_risk_count} samples</strong> contain dangerous combinations of carbapenemases with colistin or tigecycline resistance genes.</p></div></div>' if high_risk_count > 0 else ''}
        
        <h3><i class="fas fa-chart-bar"></i> Key Statistics</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Metric</th>
                        <th>Count</th>
                        <th>Details</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td>Total Samples Analyzed</td>
                        <td><strong>{total_samples}</strong></td>
                        <td>Complete genomic analysis with all databases</td>
                    </tr>
                    <tr>
                        <td>Pasteur STs</td>
                        <td><strong>{len(patterns.get('pasteur_st_distribution', {}))}</strong></td>
                        <td>MLST typing (Pasteur scheme)</td>
                    </tr>
                    <tr>
                        <td>International Clones</td>
                        <td><strong>{len(patterns.get('international_clone_distribution', {}))}</strong></td>
                        <td>IC classification</td>
                    </tr>
                    <tr>
                        <td>Capsule (K) Loci</td>
                        <td><strong>{len(patterns.get('k_locus_distribution', {}))}</strong></td>
                        <td>Kaptive capsule typing</td>
                    </tr>
                    <tr>
                        <td>Total AMR Genes</td>
                        <td><strong>{total_amr_genes}</strong></td>
                        <td>Across all AMR databases</td>
                    </tr>
                    <tr>
                        <td>Carbapenemase Genes</td>
                        <td><span class="badge {'badge-critical' if carbapenemase_count > 0 else 'badge-low'}">{carbapenemase_count}</span></td>
                        <td>OXA, NDM, VIM, IMP, KPC types</td>
                    </tr>
                    <tr>
                        <td>Virulence Genes</td>
                        <td><strong>{total_virulence_genes}</strong></td>
                        <td>Biofilm formation, virulence factors</td>
                    </tr>
                    <tr>
                        <td>Environmental Markers</td>
                        <td><span class="badge {'badge-high' if environmental_marker_count > 0 else 'badge-low'}">{environmental_marker_count}</span></td>
                        <td>Heavy metal, biocide, stress response</td>
                    </tr>
                    <tr>
                        <td>High-Risk Combinations</td>
                        <td><span class="badge {'badge-critical' if high_risk_count > 0 else 'badge-low'}">{high_risk_count}</span></td>
                        <td>Carbapenemase + last-resort resistance</td>
                    </tr>
                    <tr>
                        <td>MDR Patterns</td>
                        <td><span class="badge {'badge-high' if mdr_count > 0 else 'badge-low'}">{mdr_count}</span></td>
                        <td>3+ resistance classes in same isolate</td>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-dna"></i> Critical Resistance Categories Tracked</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
            <div class="database-section" style="border-left: 5px solid #F44336;">
                <h4><i class="fas fa-skull-crossbones"></i> Carbapenemases ({carbapenemase_count})</h4>
                <p>OXA-23, OXA-58, NDM, VIM, IMP, KPC</p>
            </div>
            
            <div class="database-section" style="border-left: 5px solid #9C27B0;">
                <h4><i class="fas fa-tint"></i> Colistin Resistance</h4>
                <p>mcr genes, pmrAB mutations, LPS modifications</p>
            </div>
            
            <div class="database-section" style="border-left: 5px solid #FF9800;">
                <h4><i class="fas fa-pills"></i> Tigecycline Resistance</h4>
                <p>tet(X) variants, efflux pumps (adeABC)</p>
            </div>
            
            <div class="database-section" style="border-left: 5px solid #4CAF50;">
                <h4><i class="fas fa-layer-group"></i> Biofilm Formation</h4>
                <p>ompA, csuABCDE, bfmRS, pili genes</p>
            </div>
            
            <div class="database-section" style="border-left: 5px solid #795548;">
                <h4><i class="fas fa-globe-africa"></i> Environmental Markers ({environmental_marker_count})</h4>
                <p>Heavy metals, biocides, stress response, plasmid transfer</p>
            </div>
        </div>
        """
        
        return html
    
    def _generate_sample_overview_section(self, kwargs: Dict) -> str:
        """Generate sample overview section - FULLY SCROLLABLE"""
        samples_data = kwargs['samples_data']
        gene_centric = kwargs['gene_centric']
        
        html = """
        <h2 class="section-header samples-header">
            <i class="fas fa-list-alt"></i> Complete Sample Overview
            <button class="print-section-btn" onclick="printSection('samples-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <input type="text" class="search-box" id="search-samples" 
               onkeyup="searchTable('samples-table', 'search-samples')" 
               placeholder="🔍 Search samples by any field...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table', 'acinetobacter_samples.csv')">
                <i class="fas fa-download"></i> Export to CSV
            </button>
            <button class="action-btn btn-success" onclick="document.getElementById('search-samples').value=''; searchTable('samples-table', 'search-samples')">
                <i class="fas fa-sync"></i> Clear Search
            </button>
        </div>
        
        <div class="master-scrollable-container">
            <table id="samples-table" class="data-table">
                <thead>
                    <tr>
                        <th class="col-sample">Sample</th>
                        <th class="col-st">Pasteur ST</th>
                        <th class="col-ic">International Clone</th>
                        <th class="col-st">Oxford ST</th>
                        <th class="col-k-locus">K Locus</th>
                        <th class="col-o-locus">O Locus</th>
                        <th class="col-capsule">Capsule Type</th>
                        <th class="col-frequency">AMR Genes</th>
                        <th class="col-frequency">Virulence Genes</th>
                        <th class="col-frequency">Other Genes</th>
                        <th class="col-frequency">Total Genes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Collect gene counts per sample
        sample_gene_counts = defaultdict(lambda: {'amr': 0, 'virulence': 0, 'other': 0})
        
        # Count AMR genes
        for db_name, genes in gene_centric.get('amr_databases', {}).items():
            for gene_data in genes:
                for genome in gene_data['genomes']:
                    sample_gene_counts[genome]['amr'] += 1
        
        # Count virulence genes
        for db_name, genes in gene_centric.get('virulence_databases', {}).items():
            for gene_data in genes:
                for genome in gene_data['genomes']:
                    sample_gene_counts[genome]['virulence'] += 1
        
        # Count other genes (from non-AMR, non-virulence databases)
        other_dbs = ['card', 'resfinder', 'argannot', 'megares', 'bacmet2', 'plasmidfinder', 'ecoh', 'ncbi']
        for db_name in other_dbs:
            if db_name in gene_centric.get('amr_databases', {}):
                for gene_data in gene_centric['amr_databases'][db_name]:
                    for genome in gene_data['genomes']:
                        sample_gene_counts[genome]['other'] += 1
        
        for sample, data in sorted(samples_data.items()):
            pasteur_st = data.get('pasteur_mlst', {}).get('ST', 'ND')
            ic = data.get('pasteur_mlst', {}).get('International_Clone', 'Unknown')
            oxford_st = data.get('oxford_mlst', {}).get('ST', 'ND')
            k_locus = data.get('kaptive', {}).get('K_Locus', 'ND')
            o_locus = data.get('kaptive', {}).get('O_Locus', 'ND')
            capsule_type = data.get('kaptive', {}).get('Capsule_Type', 'ND')
            
            amr_count = sample_gene_counts[sample]['amr']
            virulence_count = sample_gene_counts[sample]['virulence']
            other_count = sample_gene_counts[sample]['other']
            total_count = amr_count + virulence_count + other_count
            
            html += f"""
                    <tr>
                        <td class="col-sample"><strong>{sample}</strong></td>
                        <td class="col-st">{pasteur_st}</td>
                        <td class="col-ic">{ic}</td>
                        <td class="col-st">{oxford_st}</td>
                        <td class="col-k-locus">{k_locus}</td>
                        <td class="col-o-locus">{o_locus}</td>
                        <td class="col-capsule">{capsule_type}</td>
                        <td class="col-frequency">{amr_count}</td>
                        <td class="col-frequency">{virulence_count}</td>
                        <td class="col-frequency">{other_count}</td>
                        <td class="col-frequency">{total_count}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_mlst_section(self, kwargs: Dict) -> str:
        """Generate MLST analysis section - FIXED TO SHOW ALL SAMPLES"""
        patterns = kwargs['patterns']
        samples_data = kwargs['samples_data']
        
        pasteur_dist = patterns.get('pasteur_st_distribution', {})
        oxford_dist = patterns.get('oxford_st_distribution', {})
        ic_dist = patterns.get('international_clone_distribution', Counter())
        
        html = f"""
        <h2 class="section-header mlst-header">
            <i class="fas fa-code-branch"></i> MLST Analysis - Dual Scheme (Pasteur & Oxford)
            <button class="print-section-btn" onclick="printSection('mlst-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>MLST Analysis for A. baumannii</h3>
                <p><strong>{len(pasteur_dist)} unique Pasteur STs</strong> and <strong>{len(oxford_dist)} unique Oxford STs</strong> identified. 
                International Clone (IC) classification helps track global lineages.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-map-marker-alt"></i> Pasteur MLST Scheme</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-st">ST</th>
                        <th class="col-frequency">Frequency</th>
                        <th class="col-ic">International Clone</th>
                        <th class="col-sample">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Calculate percentages
        total_pasteur = sum(item['count'] for item in pasteur_dist.values()) if isinstance(next(iter(pasteur_dist.values()), {}), dict) else sum(pasteur_dist.values())
        
        for st, data in sorted(pasteur_dist.items(), key=lambda x: x[1]['count'] if isinstance(x[1], dict) else x[1], reverse=True):
            if st == 'ND':
                continue
                
            if isinstance(data, dict):
                count = data['count']
                frequency = data.get('frequency_display', f"{count} ({(count / total_pasteur * 100) if total_pasteur > 0 else 0:.1f}%)")
            else:
                count = data
                frequency = f"{count} ({(count / total_pasteur * 100) if total_pasteur > 0 else 0:.1f}%)"
            
            # Find IC for this ST
            ic_list = []
            for sample, sample_data in samples_data.items():
                if sample_data.get('pasteur_mlst', {}).get('ST') == st:
                    ic = sample_data.get('pasteur_mlst', {}).get('International_Clone', 'Unknown')
                    if ic != 'Unknown' and ic not in ic_list:
                        ic_list.append(ic)
            
            ic_display = ', '.join(ic_list) if ic_list else 'Unknown'
            
            # Find samples with this ST - NO TRUNCATION
            samples = []
            for sample, sample_data in samples_data.items():
                if sample_data.get('pasteur_mlst', {}).get('ST') == st:
                    samples.append(sample)
            
            # FIX: Show ALL samples, no truncation
            sample_display = ', '.join(samples)  # REMOVED truncation
            
            html += f"""
                    <tr>
                        <td class="col-st"><strong>ST{st}</strong></td>
                        <td class="col-frequency">{frequency}</td>
                        <td class="col-ic">{ic_display}</td>
                        <td class="col-sample">{sample_display}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-university"></i> Oxford MLST Scheme</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-st">ST</th>
                        <th class="col-frequency">Frequency</th>
                        <th class="col-sample">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Calculate percentages for Oxford
        total_oxford = sum(item['count'] for item in oxford_dist.values()) if isinstance(next(iter(oxford_dist.values()), {}), dict) else sum(oxford_dist.values())
        
        for st, data in sorted(oxford_dist.items(), key=lambda x: x[1]['count'] if isinstance(x[1], dict) else x[1], reverse=True):
            if st == 'ND':
                continue
                
            if isinstance(data, dict):
                count = data['count']
                frequency = data.get('frequency_display', f"{count} ({(count / total_oxford * 100) if total_oxford > 0 else 0:.1f}%)")
            else:
                count = data
                frequency = f"{count} ({(count / total_oxford * 100) if total_oxford > 0 else 0:.1f}%)"
            
            # Find samples with this ST - NO TRUNCATION
            samples = []
            for sample, sample_data in samples_data.items():
                if sample_data.get('oxford_mlst', {}).get('ST') == st:
                    samples.append(sample)
            
            # FIX: Show ALL samples, no truncation
            sample_display = ', '.join(samples)  # REMOVED truncation
            
            html += f"""
                    <tr>
                        <td class="col-st"><strong>ST{st}</strong></td>
                        <td class="col-frequency">{frequency}</td>
                        <td class="col-sample">{sample_display}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-globe"></i> International Clone Distribution</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-ic">International Clone</th>
                        <th class="col-frequency">Frequency</th>
                        <th class="col-st">Common STs</th>
                        <th class="col-sample">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total_ic = sum(ic_dist.values())
        for ic, count in ic_dist.most_common():
            if ic == 'Unknown':
                continue
                
            frequency = f"{count} ({(count / total_ic) * 100 if total_ic > 0 else 0:.1f}%)"
            
            # Find STs for this IC
            sts = []
            for sample, data in samples_data.items():
                if data.get('pasteur_mlst', {}).get('International_Clone') == ic:
                    st = data.get('pasteur_mlst', {}).get('ST', 'ND')
                    if st != 'ND' and st not in sts:
                        sts.append(st)
            
            st_list = ', '.join([f"ST{s}" for s in sts])  # REMOVED truncation
            
            # Find samples for this IC - NO TRUNCATION
            ic_samples = []
            for sample, data in samples_data.items():
                if data.get('pasteur_mlst', {}).get('International_Clone') == ic:
                    ic_samples.append(sample)
            
            sample_list = ', '.join(ic_samples)  # REMOVED truncation
            
            html += f"""
                    <tr>
                        <td class="col-ic"><strong>{ic}</strong></td>
                        <td class="col-frequency">{frequency}</td>
                        <td class="col-st">{st_list}</td>
                        <td class="col-sample">{sample_list}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_kaptive_section(self, kwargs: Dict) -> str:
        """Generate Kaptive capsule typing section - FIXED TO SHOW ALL SAMPLES"""
        patterns = kwargs['patterns']
        samples_data = kwargs['samples_data']
        
        k_dist = patterns.get('k_locus_distribution', {})
        o_dist = patterns.get('o_locus_distribution', {})
        capsule_dist = patterns.get('capsule_type_distribution', {})
        
        html = f"""
        <h2 class="section-header kaptive-header">
            <i class="fas fa-shield-alt"></i> Capsule Typing (Kaptive) Analysis
            <button class="print-section-btn" onclick="printSection('kaptive-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Capsule and Lipooligosaccharide Typing</h3>
                <p>Kaptive analysis identifies <strong>{len(k_dist)} capsule (K) loci</strong> and 
                <strong>{len(o_dist)} lipooligosaccharide (O) loci</strong> across samples. 
                Capsule type is critical for virulence and immune evasion in A. baumannii.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-capsules"></i> Capsule (K) Locus Distribution</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-k-locus">K Locus</th>
                        <th class="col-frequency">Frequency</th>
                        <th class="col-sample">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total_k = sum(item['count'] for item in k_dist.values()) if k_dist and isinstance(next(iter(k_dist.values()), {}), dict) else sum(k_dist.values())
        
        for k_locus, data in sorted(k_dist.items(), key=lambda x: x[1]['count'] if isinstance(x[1], dict) else x[1], reverse=True):
            if k_locus == 'ND':
                continue
                
            if isinstance(data, dict):
                count = data['count']
                frequency = data.get('frequency_display', f"{count} ({(count / total_k * 100) if total_k > 0 else 0:.1f}%)")
            else:
                count = data
                frequency = f"{count} ({(count / total_k * 100) if total_k > 0 else 0:.1f}%)"
            
            # Find samples with this K locus - NO TRUNCATION
            samples = []
            for sample, sample_data in samples_data.items():
                if sample_data.get('kaptive', {}).get('K_Locus') == k_locus:
                    samples.append(sample)
            
            # FIX: Show ALL samples, no truncation
            sample_display = ', '.join(samples)  # REMOVED truncation
            
            html += f"""
                    <tr>
                        <td class="col-k-locus"><strong>{k_locus}</strong></td>
                        <td class="col-frequency">{frequency}</td>
                        <td class="col-sample">{sample_display}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-layer-group"></i> Lipooligosaccharide (O) Locus Distribution</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-o-locus">O Locus</th>
                        <th class="col-frequency">Frequency</th>
                        <th class="col-sample">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total_o = sum(item['count'] for item in o_dist.values()) if o_dist and isinstance(next(iter(o_dist.values()), {}), dict) else sum(o_dist.values())
        
        for o_locus, data in sorted(o_dist.items(), key=lambda x: x[1]['count'] if isinstance(x[1], dict) else x[1], reverse=True):
            if o_locus == 'ND':
                continue
                
            if isinstance(data, dict):
                count = data['count']
                frequency = data.get('frequency_display', f"{count} ({(count / total_o * 100) if total_o > 0 else 0:.1f}%)")
            else:
                count = data
                frequency = f"{count} ({(count / total_o * 100) if total_o > 0 else 0:.1f}%)"
            
            # Find samples with this O locus - NO TRUNCATION
            samples = []
            for sample, sample_data in samples_data.items():
                if sample_data.get('kaptive', {}).get('O_Locus') == o_locus:
                    samples.append(sample)
            
            # FIX: Show ALL samples, no truncation
            sample_display = ', '.join(samples)  # REMOVED truncation
            
            html += f"""
                    <tr>
                        <td class="col-o-locus"><strong>{o_locus}</strong></td>
                        <td class="col-frequency">{frequency}</td>
                        <td class="col-sample">{sample_display}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-dna"></i> Complete Capsule Types (K:O)</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-capsule">Capsule Type</th>
                        <th class="col-frequency">Frequency</th>
                        <th class="col-sample">Samples</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        total_capsule = sum(item['count'] for item in capsule_dist.values()) if capsule_dist and isinstance(next(iter(capsule_dist.values()), {}), dict) else sum(capsule_dist.values())
        
        for capsule_type, data in sorted(capsule_dist.items(), key=lambda x: x[1]['count'] if isinstance(x[1], dict) else x[1], reverse=True):
            if capsule_type == 'ND' or 'ND:' in capsule_type or ':ND' in capsule_type:
                continue
                
            if isinstance(data, dict):
                count = data['count']
                frequency = data.get('frequency_display', f"{count} ({(count / total_capsule * 100) if total_capsule > 0 else 0:.1f}%)")
            else:
                count = data
                frequency = f"{count} ({(count / total_capsule * 100) if total_capsule > 0 else 0:.1f}%)"
            
            # Find samples with this capsule type - NO TRUNCATION
            samples = []
            for sample, sample_data in samples_data.items():
                if sample_data.get('kaptive', {}).get('Capsule_Type') == capsule_type:
                    samples.append(sample)
            
            # FIX: Show ALL samples, no truncation
            sample_display = ', '.join(samples)  # REMOVED truncation
            
            html += f"""
                    <tr>
                        <td class="col-capsule"><strong>{capsule_type}</strong></td>
                        <td class="col-frequency">{frequency}</td>
                        <td class="col-sample">{sample_display}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html

    def _generate_plasmid_section(self, kwargs: Dict) -> str:
        """Generate plasmid analysis section"""
        plasmid_analysis = kwargs['plasmid_analysis']
        
        if not plasmid_analysis.get('plasmid_databases'):
            return """
            <h2 class="section-header" style="border-color: #2196F3;">
                <i class="fas fa-dna"></i> Plasmid Analysis
                <button class="print-section-btn" onclick="printSection('plasmid-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div>
                    <h3>No Plasmid Data Available</h3>
                    <p>No PlasmidFinder or plasmid marker reports were found or successfully parsed.</p>
                </div>
            </div>
            """
        
        # Extract data
        plasmid_stats = plasmid_analysis.get('plasmid_summary_stats', {})
        plasmid_frequencies = plasmid_analysis.get('plasmid_frequencies', [])
        plasmid_categories = plasmid_analysis.get('plasmid_categories', {})
        high_freq_plasmids = plasmid_analysis.get('high_frequency_plasmids', [])
        unique_patterns = plasmid_analysis.get('unique_plasmid_patterns', {})
        sample_profiles = plasmid_analysis.get('sample_plasmid_profiles', {})
        
        html = f"""
        <h2 class="section-header" style="border-color: #2196F3;">
            <i class="fas fa-dna"></i> Plasmid Analysis - Horizontal Gene Transfer Tracking
            <button class="print-section-btn" onclick="printSection('plasmid-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Plasmid and Mobile Genetic Element Analysis</h3>
                <p>Tracking <strong>{plasmid_stats.get('total_plasmid_markers', 0)} plasmid markers</strong> across 
                <strong>{plasmid_stats.get('samples_with_plasmids', 0)} samples</strong> ({plasmid_stats.get('plasmid_prevalence', 0):.1f}% prevalence). 
                Plasmids are key vectors for antibiotic resistance gene transfer in A. baumannii.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> Plasmid Analysis Summary</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-gene">Metric</th>
                        <th class="col-frequency">Value</th>
                        <th class="col-category">Details</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td class="col-gene">Total Plasmid Markers</td>
                        <td class="col-frequency"><strong>{plasmid_stats.get('total_plasmid_markers', 0)}</strong></td>
                        <td class="col-category">Unique plasmid-associated genes</td>
                    </tr>
                    <tr>
                        <td class="col-gene">Samples with Plasmids</td>
                        <td class="col-frequency"><strong>{plasmid_stats.get('samples_with_plasmids', 0)}</strong> ({plasmid_stats.get('plasmid_prevalence', 0):.1f}%)</td>
                        <td class="col-category">Samples containing plasmid markers</td>
                    </tr>
                    <tr>
                        <td class="col-gene">Total Plasmid Occurrences</td>
                        <td class="col-frequency"><strong>{plasmid_stats.get('total_plasmid_occurrences', 0)}</strong></td>
                        <td class="col-category">Sum of all plasmid marker detections</td>
                    </tr>
                    <tr>
                        <td class="col-gene">Unique Plasmid Patterns</td>
                        <td class="col-frequency"><strong>{plasmid_stats.get('unique_plasmid_patterns', 0)}</strong></td>
                        <td class="col-category">Distinct plasmid marker combinations</td>
                    </tr>
                    <tr>
                        <td class="col-gene">High-Frequency Plasmids</td>
                        <td class="col-frequency"><span class="badge {'badge-high' if high_freq_plasmids else 'badge-low'}">{len(high_freq_plasmids)}</span></td>
                        <td class="col-category">Plasmids in >30% of samples</td>
                    </tr>
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-dna"></i> Plasmid Marker Frequency</h3>
        <div class="master-scrollable-container">
            <table id="plasmid-table" class="data-table">
                <thead>
                    <tr>
                        <th class="col-gene">Plasmid Marker</th>
                        <th class="col-category">Category</th>
                        <th class="col-frequency">Frequency</th>
                        <th class="col-database">Database</th>
                        <th class="col-genomes">Genomes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for plasmid_data in plasmid_frequencies:
            genomes = plasmid_data.get('genomes', [])
            genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
            
            # Highlight high frequency plasmids
            is_high_freq = plasmid_data['count'] >= (kwargs['total_samples'] * 0.3)
            marker_display = f"<strong>{plasmid_data['plasmid_marker']}</strong>" + (" 🔥" if is_high_freq else "")
            
            html += f"""
                    <tr>
                        <td class="col-gene">{marker_display}</td>
                        <td class="col-category"><span class="badge {'badge-info' if plasmid_data['category'] == 'Colicin plasmid' else 'badge-warning' if plasmid_data['category'] == 'Replication protein' else 'badge-success' if plasmid_data['category'] == 'Mobility gene' else 'badge-secondary'}">{plasmid_data['category']}</span></td>
                        <td class="col-frequency"><span class="frequency-display">{plasmid_data.get('frequency_display', f"{plasmid_data['count']} ({plasmid_data.get('percentage', 0):.1f}%)")}</span></td>
                        <td class="col-database">{plasmid_data['database']}</td>
                        <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-tags"></i> Plasmid Categories</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-category">Plasmid Category</th>
                        <th class="col-frequency">Unique Markers</th>
                        <th class="col-frequency">Total Occurrences</th>
                        <th class="col-gene">Top Markers</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for category, plasmids in plasmid_categories.items():
            unique_markers = len(set(p['plasmid_marker'] for p in plasmids))
            total_occurrences = sum(p['count'] for p in plasmids)
            top_markers = ', '.join([p['plasmid_marker'] for p in plasmids[:3]])
            
            html += f"""
                    <tr>
                        <td class="col-category"><strong>{category}</strong></td>
                        <td class="col-frequency">{unique_markers}</td>
                        <td class="col-frequency">{total_occurrences}</td>
                        <td class="col-gene">{top_markers}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        # High frequency plasmids section
        if high_freq_plasmids:
            html += f"""
            <h3 style="margin-top: 30px;"><i class="fas fa-exclamation-triangle"></i> High-Frequency Plasmids ({len(high_freq_plasmids)})</h3>
            <div class="alert-box alert-warning">
                <i class="fas fa-radiation fa-2x"></i>
                <div>
                    <h3>⚠️ Widespread Plasmids Detected</h3>
                    <p><strong>{len(high_freq_plasmids)} plasmid markers</strong> are present in >30% of samples, 
                    indicating potential horizontal transfer and epidemic spread.</p>
                </div>
            </div>
            
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th class="col-gene">Plasmid Marker</th>
                            <th class="col-category">Category</th>
                            <th class="col-frequency">Frequency</th>
                            <th class="col-frequency">Prevalence</th>
                            <th class="col-genomes">Sample Count</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for plasmid_data in high_freq_plasmids:
                prevalence = (plasmid_data['count'] / kwargs['total_samples']) * 100
                
                html += f"""
                        <tr>
                            <td class="col-gene"><strong>{plasmid_data['plasmid_marker']}</strong> 🔥</td>
                            <td class="col-category"><span class="badge badge-warning">{plasmid_data['category']}</span></td>
                            <td class="col-frequency">{plasmid_data.get('frequency_display', f"{plasmid_data['count']} ({plasmid_data.get('percentage', 0):.1f}%)")}</td>
                            <td class="col-frequency">{prevalence:.1f}%</td>
                            <td class="col-genomes">{plasmid_data['count']} samples</td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        
        # Unique plasmid patterns
        if unique_patterns:
            html += f"""
            <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> Unique Plasmid Patterns ({len(unique_patterns)})</h3>
            <p>Distinct combinations of plasmid markers across samples:</p>
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th class="col-gene">Plasmid Pattern</th>
                            <th class="col-frequency">Frequency</th>
                            <th class="col-sample">Samples</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for pattern, samples in sorted(unique_patterns.items(), key=lambda x: len(x[1]), reverse=True)[:20]:
                pattern_str = ', '.join(pattern)
                sample_list = ', '.join(samples)  # NO TRUNCATION
                
                html += f"""
                        <tr>
                            <td class="col-gene"><strong>{pattern_str}</strong></td>
                            <td class="col-frequency">{len(samples)}</td>
                            <td class="col-sample">{sample_list}</td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        
        # Sample plasmid profiles
        if sample_profiles:
            html += f"""
            <h3 style="margin-top: 30px;"><i class="fas fa-vial"></i> Sample Plasmid Profiles ({len(sample_profiles)} samples)</h3>
            <p>Plasmid marker composition for each sample:</p>
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th class="col-sample">Sample</th>
                            <th class="col-frequency">Plasmid Count</th>
                            <th class="col-gene">Plasmid Markers</th>
                            <th class="col-category">Categories</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for sample, plasmids in sorted(sample_profiles.items()):
                plasmid_names = [p['marker'] for p in plasmids]
                plasmid_list = ', '.join(plasmid_names)
                categories = ', '.join(sorted(set(p['category'] for p in plasmids)))
                
                html += f"""
                        <tr>
                            <td class="col-sample"><strong>{sample}</strong></td>
                            <td class="col-frequency">{len(plasmids)}</td>
                            <td class="col-gene">{plasmid_list}</td>
                            <td class="col-category">{categories}</td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        
        # Add export button
        html += """
        <div class="action-buttons" style="margin-top: 30px;">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('plasmid-table', 'acinetobacter_plasmid_analysis.csv')">
                <i class="fas fa-download"></i> Export Plasmid Analysis
            </button>
        </div>
        """
        
        return html
    
    def _generate_amr_section(self, kwargs: Dict) -> str:
        """Generate AMR genes section with combined frequency display"""
        gene_centric = kwargs['gene_centric']
        amr_databases = gene_centric.get('amr_databases', {})
        
        if not amr_databases:
            return """
            <h2 class="section-header amr-header">
                <i class="fas fa-biohazard"></i> Antimicrobial Resistance Genes
                <button class="print-section-btn" onclick="printSection('amr-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div>
                    <h3>No AMR Data Available</h3>
                    <p>No AMRfinder or ABRicate AMR database reports were found or successfully parsed.</p>
                </div>
            </div>
            """
        
        html = """
        <h2 class="section-header amr-header">
            <i class="fas fa-biohazard"></i> Antimicrobial Resistance Genes - Comprehensive Analysis
            <button class="print-section-btn" onclick="printSection('amr-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>AMR Gene Analysis for A. baumannii</h3>
                <p>Each AMR gene is shown with ALL genomes that contain it. <strong>Frequency shows count and percentage</strong> 
                in format "count (percentage%)". Carbapenemases are highlighted in red.</p>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-amr" 
               onkeyup="searchTable('amr-table', 'search-amr')" 
               placeholder="🔍 Search AMR genes by name, category, or database...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('amr-table', 'acinetobacter_amr_genes.csv')">
                <i class="fas fa-download"></i> Export All AMR Genes
            </button>
            <button class="action-btn btn-success" onclick="document.getElementById('search-amr').value='carbapenemase'; searchTable('amr-table', 'search-amr')">
                <i class="fas fa-search"></i> Show Carbapenemases
            </button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-amr').value='OXA'; searchTable('amr-table', 'search-amr')">
                <i class="fas fa-skull-crossbones"></i> Show OXA genes
            </button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-amr').value=''; searchTable('amr-table', 'search-amr')">
                <i class="fas fa-sync"></i> Clear Search
            </button>
        </div>
        
        <h3><i class="fas fa-shield-virus"></i> All AMR Genes Across Databases</h3>
        <div class="master-scrollable-container">
            <table id="amr-table" class="data-table">
                <thead>
                    <tr>
                        <th class="col-gene">Gene</th>
                        <th class="col-category">Category</th>
                        <th class="col-database">Database</th>
                        <th class="col-frequency">Frequency</th>
                        <th class="col-risk">Risk Level</th>
                        <th class="col-genomes">Genomes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Combine all AMR genes
        all_amr_genes = []
        for db_name, genes in amr_databases.items():
            for gene_data in genes:
                all_amr_genes.append(gene_data)
        
        # Sort by count
        all_amr_genes.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_amr_genes:
            gene = gene_data['gene']
            category = gene_data['category']
            database = gene_data['database']
            frequency = gene_data.get('frequency_display', f"{gene_data.get('count', 0)} ({gene_data.get('percentage', 0):.1f}%)")
            risk_level = gene_data.get('risk_level', 'Standard')
            genomes = gene_data.get('genomes', [])
            
            # Highlight carbapenemases
            is_carbapenemase = 'Carbapenemases' in category
            gene_display = f"<strong>{gene}</strong>" + (" 🔥" if is_carbapenemase else "")
            
            # Category chip class
            chip_class = f"chip-{category.lower().replace(' ', '-').replace('/', '-')}"
            
            # Create genome tags
            genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
            
            html += f"""
                    <tr>
                        <td class="col-gene">{gene_display}</td>
                        <td class="col-category"><span class="category-chip {chip_class}">{category}</span></td>
                        <td class="col-database">{database}</td>
                        <td class="col-frequency"><span class="frequency-display">{frequency}</span></td>
                        <td class="col-risk"><span class="badge {'badge-high' if risk_level == 'HIGH' else 'badge-critical' if risk_level == 'CRITICAL' else 'badge-low'}">{risk_level}</span></td>
                        <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-database"></i> AMR Databases Summary</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-database">Database</th>
                        <th class="col-frequency">Unique Genes</th>
                        <th class="col-frequency">Total Occurrences</th>
                        <th class="col-frequency">Carbapenemases</th>
                        <th class="col-gene">Top Gene</th>
                        <th class="col-frequency">Top Gene Frequency</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for db_name, genes in amr_databases.items():
            db_display = db_name.upper() if db_name != 'amrfinder' else 'AMRfinder'
            
            # Calculate statistics
            total_genes = len(genes)
            total_count = sum(g['count'] for g in genes)
            carbapenemase_count = sum(1 for g in genes if 'Carbapenemases' in g['category'])
            top_gene = genes[0]['gene'] if genes else 'None'
            top_gene_freq = genes[0]['frequency_display'] if genes else '0 (0%)'
            
            html += f"""
                    <tr>
                        <td class="col-database"><strong>{db_display}</strong></td>
                        <td class="col-frequency">{total_genes}</td>
                        <td class="col-frequency">{total_count}</td>
                        <td class="col-frequency"><span class="badge {'badge-critical' if carbapenemase_count > 0 else 'badge-low'}">{carbapenemase_count}</span></td>
                        <td class="col-gene">{top_gene}</td>
                        <td class="col-frequency">{top_gene_freq}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_virulence_section(self, kwargs: Dict) -> str:
        """Generate virulence genes section with combined frequency display"""
        gene_centric = kwargs['gene_centric']
        virulence_databases = gene_centric.get('virulence_databases', {})
        
        if not virulence_databases:
            return """
            <h2 class="section-header virulence-header">
                <i class="fas fa-virus"></i> Virulence Genes
                <button class="print-section-btn" onclick="printSection('virulence-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div>
                    <h3>No Virulence Data Available</h3>
                    <p>No VFDB, Victors, or other virulence database reports were found or successfully parsed.</p>
                </div>
            </div>
            """
        
        html = """
        <h2 class="section-header virulence-header">
            <i class="fas fa-virus"></i> Virulence Genes - Biofilm and Pathogenicity Factors
            <button class="print-section-btn" onclick="printSection('virulence-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Virulence Gene Analysis for A. baumannii</h3>
                <p>Focus on <strong>biofilm formation genes</strong> which are critical for A. baumannii persistence 
                and hospital outbreaks. Each gene shows frequency as "count (percentage%)".</p>
            </div>
        </div>
        
        <input type="text" class="search-box" id="search-virulence" 
               onkeyup="searchTable('virulence-table', 'search-virulence')" 
               placeholder="🔍 Search virulence genes by name, category, or database...">
        
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('virulence-table', 'acinetobacter_virulence_genes.csv')">
                <i class="fas fa-download"></i> Export All Virulence Genes
            </button>
            <button class="action-btn btn-success" onclick="document.getElementById('search-virulence').value='biofilm'; searchTable('virulence-table', 'search-virulence')">
                <i class="fas fa-layer-group"></i> Show Biofilm Genes
            </button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-virulence').value=''; searchTable('virulence-table', 'search-virulence')">
                <i class="fas fa-sync"></i> Clear Search
            </button>
        </div>
        
        <h3><i class="fas fa-virus"></i> All Virulence Genes Across Databases</h3>
        <div class="master-scrollable-container">
            <table id="virulence-table" class="data-table">
                <thead>
                    <tr>
                        <th class="col-gene">Gene</th>
                        <th class="col-category">Category</th>
                        <th class="col-database">Database</th>
                        <th class="col-frequency">Frequency</th>
                        <th class="col-genomes">Genomes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Combine all virulence genes
        all_virulence_genes = []
        for db_name, genes in virulence_databases.items():
            for gene_data in genes:
                all_virulence_genes.append(gene_data)
        
        # Sort by count
        all_virulence_genes.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_virulence_genes:
            gene = gene_data['gene']
            category = gene_data['category']
            database = gene_data['database']
            frequency = gene_data.get('frequency_display', f"{gene_data['count']} ({gene_data.get('percentage', 0):.1f}%)")
            genomes = gene_data.get('genomes', [])
            
            # Highlight biofilm genes
            is_biofilm = 'Biofilm' in category
            gene_display = f"<strong>{gene}</strong>" + (" 🏥" if is_biofilm else "")
            
            # Category chip class
            chip_class = f"chip-{category.lower().replace(' ', '-').replace('/', '-')}"
            
            # Create genome tags
            genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
            
            html += f"""
                    <tr>
                        <td class="col-gene">{gene_display}</td>
                        <td class="col-category"><span class="category-chip {chip_class}">{category}</span></td>
                        <td class="col-database">{database}</td>
                        <td class="col-frequency"><span class="frequency-display">{frequency}</span></td>
                        <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-database"></i> Virulence Databases Summary</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-database">Database</th>
                        <th class="col-frequency">Unique Genes</th>
                        <th class="col-frequency">Total Occurrences</th>
                        <th class="col-frequency">Biofilm Genes</th>
                        <th class="col-gene">Top Gene</th>
                        <th class="col-frequency">Top Gene Frequency</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for db_name, genes in virulence_databases.items():
            db_display = db_name.upper()
            
            # Calculate statistics
            total_genes = len(genes)
            total_count = sum(g['count'] for g in genes)
            biofilm_count = sum(1 for g in genes if 'Biofilm' in g['category'])
            top_gene = genes[0]['gene'] if genes else 'None'
            top_gene_freq = genes[0]['frequency_display'] if genes else '0 (0%)'
            
            html += f"""
                    <tr>
                        <td class="col-database"><strong>{db_display}</strong></td>
                        <td class="col-frequency">{total_genes}</td>
                        <td class="col-frequency">{total_count}</td>
                        <td class="col-frequency"><span class="badge {'badge-high' if biofilm_count > 0 else 'badge-low'}">{biofilm_count}</span></td>
                        <td class="col-gene">{top_gene}</td>
                        <td class="col-frequency">{top_gene_freq}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        """
        
        return html
    
    def _generate_environmental_section(self, kwargs: Dict) -> str:
        """Generate environmental resistance markers section"""
        gene_centric = kwargs['gene_centric']
        environmental_summary = gene_centric.get('environmental_summary', {})
        environmental_databases = gene_centric.get('environmental_databases', {})
        
        if not environmental_summary and not environmental_databases:
            return """
            <h2 class="section-header environmental-header">
                <i class="fas fa-globe-africa"></i> Environmental Resistance & Co-Selection Markers
                <button class="print-section-btn" onclick="printSection('environmental-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div>
                    <h3>No Environmental Data Available</h3>
                    <p>No BACMET2 or environmental resistance marker reports were found or successfully parsed.</p>
                </div>
            </div>
            """
        
        # Count total environmental markers
        total_env_genes = 0
        total_env_occurrences = 0
        for category, data in environmental_summary.items():
            total_env_genes += data.get('total_genes', 0)
            total_env_occurrences += data.get('total_occurrences', 0)
        
        html = f"""
        <h2 class="section-header environmental-header">
            <i class="fas fa-globe-africa"></i> Environmental Resistance & Co-Selection Markers
            <button class="print-section-btn" onclick="printSection('environmental-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Environmental Co-Selection Analysis</h3>
                <p>Tracking <strong>{total_env_genes} environmental resistance markers</strong> across {total_env_occurrences} occurrences. 
                These genes can co-select for antibiotic resistance in hospital environments through heavy metal resistance, 
                biocide resistance, and stress response mechanisms.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> Environmental Marker Categories</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-category">Category</th>
                        <th class="col-frequency">Unique Genes</th>
                        <th class="col-frequency">Total Occurrences</th>
                        <th class="col-frequency">Percentage of Samples</th>
                        <th class="col-gene">Top Genes</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for category, data in sorted(environmental_summary.items(), key=lambda x: x[1]['total_occurrences'], reverse=True):
            total_genes = data.get('total_genes', 0)
            total_occurrences = data.get('total_occurrences', 0)
            percentage = data.get('percentage_of_samples', 0)
            
            # Get top 3 genes
            top_genes = []
            if 'genes' in data and data['genes']:
                for gene_data in data['genes'][:3]:
                    top_genes.append(f"{gene_data['gene']} ({gene_data['count']})")
            
            top_genes_display = ', '.join(top_genes) if top_genes else 'None'
            
            html += f"""
                    <tr>
                        <td class="col-category"><strong>{category}</strong></td>
                        <td class="col-frequency">{total_genes}</td>
                        <td class="col-frequency">{total_occurrences}</td>
                        <td class="col-frequency">{percentage:.1f}%</td>
                        <td class="col-gene">{top_genes_display}</td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-dna"></i> Detailed Environmental Markers</h3>
        """
        
        # Show detailed lists for each category
        for category, data in sorted(environmental_summary.items(), key=lambda x: x[1]['total_occurrences'], reverse=True):
            if data.get('genes'):
                total_genes = data.get('total_genes', 0)
                total_occurrences = data.get('total_occurrences', 0)
                
                html += f"""
                <div class="database-section" style="margin-bottom: 30px;">
                    <h4>{category} ({total_genes} unique genes, {total_occurrences} total occurrences)</h4>
                    <div class="master-scrollable-container" style="max-height: 400px;">
                        <table class="data-table">
                            <thead>
                                <tr>
                                    <th class="col-gene">Gene</th>
                                    <th class="col-frequency">Frequency</th>
                                    <th class="col-database">Database</th>
                                    <th class="col-genomes">Genomes</th>
                                </tr>
                            </thead>
                            <tbody>
                """
                
                for gene_data in data['genes']:
                    genomes = gene_data.get('genomes', [])
                    genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
                    
                    html += f"""
                            <tr>
                                <td class="col-gene"><strong>{gene_data['gene']}</strong></td>
                                <td class="col-frequency"><span class="frequency-display">{gene_data.get('frequency_display', f"{gene_data['count']} ({gene_data.get('percentage', 0):.1f}%)")}</span></td>
                                <td class="col-database">{gene_data['database']}</td>
                                <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
                            </tr>
                    """
                
                html += """
                            </tbody>
                        </table>
                    </div>
                </div>
                """
        
        # Show BACMET2 database if available
        if 'bacmet2' in environmental_databases:
            html += """
            <h3 style="margin-top: 30px;"><i class="fas fa-database"></i> BACMET2 Database Analysis</h3>
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th class="col-gene">Gene</th>
                            <th class="col-category">Category</th>
                            <th class="col-frequency">Frequency</th>
                            <th class="col-genomes">Genomes</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for gene_data in environmental_databases['bacmet2']:
                genomes = gene_data.get('genomes', [])
                genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
                
                html += f"""
                        <tr>
                            <td class="col-gene"><strong>{gene_data['gene']}</strong></td>
                            <td class="col-category"><span class="category-chip chip-bacmet2">{gene_data['category']}</span></td>
                            <td class="col-frequency"><span class="frequency-display">{gene_data.get('frequency_display', f"{gene_data['count']} ({gene_data.get('percentage', 0):.1f}%)")}</span></td>
                            <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        
        return html
    
    def _generate_categories_section(self, kwargs: Dict) -> str:
        """Generate gene categories section - FULLY SCROLLABLE"""
        gene_centric = kwargs['gene_centric']
        category_data = gene_centric.get('gene_categories', {})
        
        if not category_data:
            return """
            <h2 class="section-header categories-header">
                <i class="fas fa-tags"></i> Gene Categories
                <button class="print-section-btn" onclick="printSection('categories-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div>
                    <h3>No Gene Category Data Available</h3>
                    <p>No genes were categorized. This may indicate parsing issues or no gene data.</p>
                </div>
            </div>
            """
        
        html = """
        <h2 class="section-header categories-header">
            <i class="fas fa-tags"></i> Gene Categories - Resistance Mechanism Analysis
            <button class="print-section-btn" onclick="printSection('categories-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Resistance Mechanism Categorization</h3>
                <p>Genes are categorized by resistance mechanism. This helps identify which resistance types 
                are most prevalent and track multidrug resistance patterns.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-pie"></i> Resistance Category Overview</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-category">Resistance Category</th>
                        <th class="col-frequency">Unique Genes</th>
                        <th class="col-frequency">Total Occurrences</th>
                        <th class="col-gene">Top 3 Genes</th>
                        <th class="col-risk">Critical Level</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Define critical levels
        critical_levels = {
            'Carbapenemases': 'CRITICAL',
            'ESBLs': 'HIGH',
            'AmpC': 'HIGH',
            'Colistin Resistance': 'CRITICAL',
            'Tigecycline Resistance': 'HIGH',
            'Biofilm Formation': 'MEDIUM',
            'Efflux Pumps': 'MEDIUM',
            'Environmental Co-Selection': 'MEDIUM',
            'BACMET2 Markers': 'MEDIUM',
            'VICTORS Virulence': 'MEDIUM',
            'Other Resistance': 'LOW',
            'Other': 'LOW'
        }
        
        for category, genes in sorted(category_data.items(), key=lambda x: len(x[1]), reverse=True):
            unique_genes = len(set(g['gene'] for g in genes))
            total_occurrences = sum(g['count'] for g in genes)
            top_genes = ', '.join([f"{g['gene']} ({g['count']})" for g in genes[:3]])
            
            critical_level = critical_levels.get(category, 'LOW')
            badge_class = 'badge-critical' if critical_level == 'CRITICAL' else 'badge-high' if critical_level == 'HIGH' else 'badge-medium'
            
            html += f"""
                    <tr>
                        <td class="col-category"><strong>{category}</strong></td>
                        <td class="col-frequency">{unique_genes}</td>
                        <td class="col-frequency">{total_occurrences}</td>
                        <td class="col-gene">{top_genes}</td>
                        <td class="col-risk"><span class="badge {badge_class}">{critical_level}</span></td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-dna"></i> Detailed Gene Lists by Category</h3>
        """
        
        # Show detailed lists for each category
        for category, genes in sorted(category_data.items(), key=lambda x: len(x[1]), reverse=True):
            if genes:
                unique_genes = len(set(g['gene'] for g in genes))
                total_occurrences = sum(g['count'] for g in genes)
                
                html += f"""
                <div class="database-section" style="margin-bottom: 30px;">
                    <h4>{category} ({unique_genes} unique genes, {total_occurrences} total occurrences)</h4>
                    <div class="master-scrollable-container" style="max-height: 400px;">
                        <table class="data-table">
                            <thead>
                                <tr>
                                    <th class="col-gene">Gene</th>
                                    <th class="col-frequency">Frequency</th>
                                    <th class="col-database">Database</th>
                                    <th class="col-genomes">Genomes</th>
                                </tr>
                            </thead>
                            <tbody>
                """
                
                for gene_data in genes:
                    genomes = gene_data.get('genomes', [])
                    genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in genomes])
                    
                    html += f"""
                            <tr>
                                <td class="col-gene"><strong>{gene_data['gene']}</strong></td>
                                <td class="col-frequency"><span class="frequency-display">{gene_data.get('frequency_display', f"{gene_data['count']} ({gene_data.get('percentage', 0):.1f}%)")}</span></td>
                                <td class="col-database">{gene_data['database']}</td>
                                <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
                            </tr>
                    """
                
                html += """
                            </tbody>
                        </table>
                    </div>
                </div>
                """
        
        return html
    
    def _generate_pattern_discovery_section(self, kwargs: Dict) -> str:
        """Generate pattern discovery section - NO TRUNCATION & NO ENVIRONMENTAL MARKERS"""
        patterns = kwargs['patterns']
        
        high_risk_combinations = patterns.get('high_risk_combinations', [])
        mdr_patterns = patterns.get('mdr_patterns', [])
        carbapenemase_patterns = patterns.get('carbapenemase_patterns', {})
        st_k_combinations = patterns.get('st_k_locus_combinations', {})
        st_capsule_combinations = patterns.get('st_capsule_combinations', {})
        
        html = """
        <h2 class="section-header patterns-header">
            <i class="fas fa-project-diagram"></i> Pattern Discovery - MDR/XDR Analysis
            <button class="print-section-btn" onclick="printSection('patterns-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Cross-Genome Pattern Discovery</h3>
                <p>Discover associations between resistance genes and identify 
                <strong>high-risk multidrug resistance (MDR) patterns</strong> in A. baumannii.</p>
            </div>
        </div>
        """
        
        # High-risk combinations - UPDATED: Remove Environmental Markers column
        if high_risk_combinations:
            html += f"""
            <h3><i class="fas fa-exclamation-triangle"></i> High-Risk Combinations ({len(high_risk_combinations)})</h3>
            <div class="alert-box alert-danger">
                <i class="fas fa-radiation fa-2x"></i>
                <div>
                    <h3>⚠️ CRITICAL ALERT: Carbapenemase + Last-Resort Resistance</h3>
                    <p><strong>{len(high_risk_combinations)} samples</strong> contain dangerous combinations of 
                    carbapenemases with colistin or tigecycline resistance genes.</p>
                </div>
            </div>
            
            <div class="master-scrollable-container">
                <table id="high-risk-table" class="data-table">
                    <thead>
                        <tr>
                            <th class="col-sample">Sample</th>
                            <th class="col-st">Pasteur ST</th>
                            <th class="col-ic">International Clone</th>
                            <th class="col-k-locus">K Locus</th>
                            <th class="col-capsule">Capsule Type</th>
                            <th class="col-gene">Carbapenemases</th>
                            <th class="col-gene">Colistin Resistance</th>
                            <th class="col-gene">Tigecycline Resistance</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for combo in high_risk_combinations:
                carbapenemases = ', '.join(combo['carbapenemases'])
                colistin_res = ', '.join(combo['colistin_resistance']) if combo['colistin_resistance'] else 'None'
                tigecycline_res = ', '.join(combo['tigecycline_resistance']) if combo['tigecycline_resistance'] else 'None'
                
                html += f"""
                        <tr>
                            <td class="col-sample"><strong>{combo['sample']}</strong></td>
                            <td class="col-st">{combo['pasteur_st']}</td>
                            <td class="col-ic">{combo['international_clone']}</td>
                            <td class="col-k-locus">{combo['k_locus']}</td>
                            <td class="col-capsule">{combo['capsule_type']}</td>
                            <td class="col-gene"><span class="badge badge-critical">{carbapenemases}</span></td>
                            <td class="col-gene"><span class="badge badge-high">{colistin_res}</span></td>
                            <td class="col-gene"><span class="badge badge-high">{tigecycline_res}</span></td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        else:
            html += """
            <div class="alert-box alert-success">
                <i class="fas fa-check-circle fa-2x"></i>
                <div>
                    <h3>No High-Risk Combinations Detected</h3>
                    <p>No samples were found with both carbapenemase and last-resort resistance genes.</p>
                </div>
            </div>
            """
        
        # Carbapenemase patterns - NO TRUNCATION
        if carbapenemase_patterns:
            html += f"""
            <h3 style="margin-top: 30px;"><i class="fas fa-skull-crossbones"></i> Carbapenemase Patterns ({len(carbapenemase_patterns)})</h3>
            <p>Distribution of carbapenemase gene combinations:</p>
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th class="col-gene">Carbapenemase Combination</th>
                            <th class="col-frequency">Frequency</th>
                            <th class="col-sample">Samples</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for carb_combo, samples in sorted(carbapenemase_patterns.items(), key=lambda x: len(x[1]), reverse=True):
                combo_str = ', '.join(carb_combo)
                # NO TRUNCATION - show all samples
                sample_list = ', '.join(samples)
                
                html += f"""
                        <tr>
                            <td class="col-gene"><strong>{combo_str}</strong></td>
                            <td class="col-frequency">{len(samples)}</td>
                            <td class="col-sample">{sample_list}</td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        
        # MDR patterns - UPDATED: Remove Environmental Markers
        if mdr_patterns:
            html += f"""
            <h3 style="margin-top: 30px;"><i class="fas fa-pills"></i> Multidrug Resistance (MDR) Patterns ({len(mdr_patterns)})</h3>
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div>
                    <h3>MDR/XDR A. baumannii Detected</h3>
                    <p><strong>{len(mdr_patterns)} samples</strong> show resistance to 3 or more antibiotic classes.</p>
                </div>
            </div>
            
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th class="col-sample">Sample</th>
                            <th class="col-frequency">Resistance Classes</th>
                            <th class="col-ic">International Clone</th>
                            <th class="col-gene">Carbapenemases</th>
                            <th class="col-gene">ESBLs</th>
                            <th class="col-gene">Colistin Resistance</th>
                            <th class="col-gene">Tigecycline Resistance</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for mdr in mdr_patterns:
                carbapenemases = ', '.join(mdr['carbapenemases']) if mdr['carbapenemases'] else 'None'
                esbls = ', '.join(mdr['esbls']) if mdr['esbls'] else 'None'
                colistin_res = ', '.join(mdr['colistin_resistance']) if mdr['colistin_resistance'] else 'None'
                tigecycline_res = ', '.join(mdr['tigecycline_resistance']) if mdr['tigecycline_resistance'] else 'None'
                
                html += f"""
                        <tr>
                            <td class="col-sample"><strong>{mdr['sample']}</strong></td>
                            <td class="col-frequency"><span class="badge badge-critical">{mdr['resistance_types']} classes</span></td>
                            <td class="col-ic">{mdr['international_clone']}</td>
                            <td class="col-gene">{carbapenemases}</td>
                            <td class="col-gene">{esbls}</td>
                            <td class="col-gene">{colistin_res}</td>
                            <td class="col-gene">{tigecycline_res}</td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        
        # ST-K Locus combinations - NO TRUNCATION
        if st_k_combinations:
            html += f"""
            <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> ST-K Locus Associations ({len(st_k_combinations)})</h3>
            <p>Common associations between sequence types and capsule loci:</p>
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th class="col-gene">ST-K Locus Combination</th>
                            <th class="col-frequency">Frequency</th>
                            <th class="col-sample">Samples</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for combo, samples in sorted(st_k_combinations.items(), key=lambda x: len(x[1]), reverse=True):
                # NO TRUNCATION - show all samples
                sample_list = ', '.join(samples)
                
                html += f"""
                        <tr>
                            <td class="col-gene"><strong>{combo}</strong></td>
                            <td class="col-frequency">{len(samples)}</td>
                            <td class="col-sample">{sample_list}</td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        
        # ST-Capsule combinations - NO TRUNCATION
        if st_capsule_combinations:
            html += f"""
            <h3 style="margin-top: 30px;"><i class="fas fa-project-diagram"></i> ST-Capsule Type Associations ({len(st_capsule_combinations)})</h3>
            <p>Common associations between sequence types and capsule types (K:O):</p>
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th class="col-gene">ST-Capsule Type Combination</th>
                            <th class="col-frequency">Frequency</th>
                            <th class="col-sample">Samples</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            
            for combo, samples in sorted(st_capsule_combinations.items(), key=lambda x: len(x[1]), reverse=True):
                # NO TRUNCATION - show all samples
                sample_list = ', '.join(samples)
                
                html += f"""
                        <tr>
                            <td class="col-gene"><strong>{combo}</strong></td>
                            <td class="col-frequency">{len(samples)}</td>
                            <td class="col-sample">{sample_list}</td>
                        </tr>
                """
            
            html += """
                    </tbody>
                </table>
            </div>
            """
        
        # Add export button
        html += """
        <div class="action-buttons" style="margin-top: 30px;">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('high-risk-table', 'high_risk_combinations.csv')">
                <i class="fas fa-download"></i> Export High-Risk Combinations
            </button>
        </div>
        """
        
        return html
    
    def _generate_database_coverage_section(self, kwargs: Dict) -> str:
        """Generate database coverage section - FULLY SCROLLABLE"""
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        
        database_coverage = patterns.get('database_coverage', {})
        database_stats = gene_centric.get('database_stats', {})
        
        if not database_coverage:
            return """
            <h2 class="section-header databases-header">
                <i class="fas fa-database"></i> Database Coverage
                <button class="print-section-btn" onclick="printSection('databases-tab')">
                    <i class="fas fa-print"></i> Print
                </button>
            </h2>
            
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div>
                    <h3>No Database Coverage Data Available</h3>
                    <p>No database coverage information was calculated.</p>
                </div>
            </div>
            """
        
        html = """
        <h2 class="section-header databases-header">
            <i class="fas fa-database"></i> Database Coverage and Statistics
            <button class="print-section-btn" onclick="printSection('databases-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Database Performance Analysis</h3>
                <p>This section shows which databases detected genes in which samples, helping assess 
                database sensitivity and coverage for A. baumannii analysis.</p>
            </div>
        </div>
        
        <h3><i class="fas fa-chart-bar"></i> Database Coverage Across Samples</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-database">Database</th>
                        <th class="col-frequency">Coverage</th>
                        <th class="col-frequency">Unique Genes</th>
                        <!-- REMOVED: Total Occurrences column -->
                        <th class="col-frequency">Critical Genes</th>
                        <!-- REMOVED: Environmental Genes column -->
                        <th class="col-database">Database Type</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        # Define database types for categorization
        database_types = {
            'amrfinder': 'AMR',
            'card': 'AMR',
            'resfinder': 'AMR',
            'argannot': 'AMR',
            'megares': 'AMR',
            'vfdb': 'Virulence',
            'victors': 'Virulence',
            'ecoli_vf': 'Virulence',
            'plasmidfinder': 'Plasmid',
            'ecoh': 'Plasmid',
            'bacmet2': 'Environmental',
            'ncbi': 'Reference'
        }
        
        for db_name, coverage_data in sorted(database_coverage.items(), key=lambda x: x[1]['coverage_percentage'], reverse=True):
            coverage = coverage_data.get('coverage_display', f"{coverage_data['samples_with_hits']} ({coverage_data['coverage_percentage']}%)")
            
            # Get additional stats from database_stats
            stats = database_stats.get(db_name, {})
            unique_genes = stats.get('total_genes', 0)
            critical_genes = stats.get('critical_genes', 0)
            
            # Determine coverage badge
            coverage_percentage = coverage_data['coverage_percentage']
            if coverage_percentage >= 80:
                coverage_badge = 'badge-low'
            elif coverage_percentage >= 50:
                coverage_badge = 'badge-medium'
            else:
                coverage_badge = 'badge-high'
            
            # Get database type
            db_type = database_types.get(db_name.lower(), 'Other')
            
            html += f"""
                    <tr>
                        <td class="col-database"><strong>{db_name.upper()}</strong></td>
                        <td class="col-frequency"><span class="badge {coverage_badge}">{coverage}</span></td>
                        <td class="col-frequency">{unique_genes}</td>
                        <!-- REMOVED: Total Occurrences cell -->
                        <td class="col-frequency"><span class="badge {'badge-critical' if critical_genes > 0 else 'badge-low'}">{critical_genes}</span></td>
                        <!-- REMOVED: Environmental Genes cell -->
                        <td class="col-database"><span class="badge {'badge-danger' if db_type == 'AMR' else 'badge-warning' if db_type == 'Virulence' else 'badge-success' if db_type == 'Plasmid' else 'badge-info' if db_type == 'Environmental' else 'badge-secondary'}">{db_type}</span></td>
                    </tr>
            """
        
        html += """
                </tbody>
            </table>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-lightbulb"></i> Database Recommendations</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">
            <div class="database-section">
                <h4><i class="fas fa-star"></i> Best Coverage Databases</h4>
        """
        
        # Find best databases
        best_databases = sorted(database_coverage.items(), key=lambda x: x[1]['coverage_percentage'], reverse=True)[:3]
        for db_name, data in best_databases:
            html += f"""
                <p><strong>{db_name.upper()}</strong>: {data['coverage_display']}</p>
            """
        
        html += """
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-exclamation-triangle"></i> Low Coverage Databases</h4>
        """
        
        # Find low coverage databases
        low_databases = sorted(database_coverage.items(), key=lambda x: x[1]['coverage_percentage'])[:3]
        for db_name, data in low_databases:
            if data['coverage_percentage'] < 50:
                html += f"""
                    <p><strong>{db_name.upper()}</strong>: {data['coverage_display']}</p>
                """
        
        html += """
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-dna"></i> Critical Gene Detection</h4>
        """
        
        # Find databases with most critical genes
        critical_dbs = sorted([(db, stats.get('critical_genes', 0)) for db, stats in database_stats.items()], 
                             key=lambda x: x[1], reverse=True)[:3]
        for db_name, critical_count in critical_dbs:
            if critical_count > 0:
                html += f"""
                    <p><strong>{db_name.upper()}</strong>: {critical_count} critical genes</p>
                """
        
        html += """
            </div>
        </div>
        
        <h3 style="margin-top: 30px;"><i class="fas fa-info-circle"></i> Database Type Legend</h3>
        <div style="display: flex; gap: 15px; flex-wrap: wrap; margin: 20px 0;">
            <span class="badge badge-danger">AMR</span> <span style="margin-right: 20px;">Antimicrobial Resistance</span>
            <span class="badge badge-warning">Virulence</span> <span style="margin-right: 20px;">Virulence Factors</span>
            <span class="badge badge-success">Plasmid</span> <span style="margin-right: 20px;">Plasmid Markers</span>
            <span class="badge badge-info">Environmental</span> <span style="margin-right: 20px;">Environmental Resistance</span>
            <span class="badge badge-secondary">Reference</span> <span>Reference Database</span>
        </div>
        """
        
        return html
    
    def _generate_export_section(self, kwargs: Dict) -> str:
        """Generate export section"""
        return """
        <h2 class="section-header export-header">
            <i class="fas fa-download"></i> Export Data and Reports
            <button class="print-section-btn" onclick="printSection('export-tab')">
                <i class="fas fa-print"></i> Print
            </button>
        </h2>
        
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>Export Data and Reports</h3>
                <p>Download comprehensive data in various formats for further analysis, reporting, and surveillance.</p>
            </div>
        </div>
        
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 30px 0;">
            <div class="dashboard-card card-export" onclick="exportTableToCSV('samples-table', 'acinetobacter_samples.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-file-csv"></i></div>
                <div class="card-label">Sample Overview CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All samples with MLST, capsule typing, gene counts</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="exportTableToCSV('amr-table', 'acinetobacter_amr_genes.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-biohazard"></i></div>
                <div class="card-label">AMR Genes CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All AMR genes with genomes and categories</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="exportTableToCSV('virulence-table', 'acinetobacter_virulence_genes.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-virus"></i></div>
                <div class="card-label">Virulence Genes CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All virulence genes with genomes and categories</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="exportTableToCSV('high-risk-table', 'high_risk_combinations.csv')">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-exclamation-triangle"></i></div>
                <div class="card-label">High-Risk Combos CSV</div>
                <p style="font-size: 0.9em; margin-top: 10px;">Carbapenemase + last-resort resistance combinations</p>
            </div>
            
            <div class="dashboard-card card-export" onclick="location.href='genius_acinetobacter_ultimate_report.json'">
                <div style="font-size: 2.5em; color: var(--export-color);"><i class="fas fa-file-code"></i></div>
                <div class="card-label">Complete JSON Data</div>
                <p style="font-size: 0.9em; margin-top: 10px;">All analysis data in structured JSON format</p>
            </div>
        </div>
        
        <h3><i class="fas fa-download"></i> Available Export Files</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th class="col-gene">File</th>
                        <th class="col-category">Description</th>
                        <th class="col-database">Format</th>
                        <th class="col-genomes">Contents</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td class="col-gene"><strong>genius_acinetobacter_ultimate_report.html</strong></td>
                        <td class="col-category">This interactive HTML report</td>
                        <td class="col-database">HTML</td>
                        <td class="col-genomes">Complete analysis with all sections</td>
                    </tr>
                    <tr>
                        <td class="col-gene"><strong>genius_acinetobacter_ultimate_report.json</strong></td>
                        <td class="col-category">Complete structured data</td>
                        <td class="col-database">JSON</td>
                        <td class="col-genomes">All analysis data for programmatic use</td>
                    </tr>
                    <tr>
                        <td class="col-gene"><strong>acinetobacter_samples.csv</strong></td>
                        <td class="col-category">Sample overview data</td>
                        <td class="col-database">CSV</td>
                        <td class="col-genomes">All samples with typing results</td>
                    </tr>
                    <tr>
                        <td class="col-gene"><strong>acinetobacter_amr_genes.csv</strong></td>
                        <td class="col-category">AMR gene analysis</td>
                        <td class="col-database">CSV</td>
                        <td class="col-genomes">All AMR genes with genomes and categories</td>
                    </tr>
                    <tr>
                        <td class="col-gene"><strong>acinetobacter_virulence_genes.csv</strong></td>
                        <td class="col-category">Virulence gene analysis</td>
                        <td class="col-database">CSV</td>
                        <td class="col-genomes">All virulence genes with genomes and categories</td>
                    </tr>
                    <tr>
                        <td class="col-gene"><strong>acinetobacter_environmental_markers.csv</strong></td>
                        <td class="col-category">Environmental resistance markers</td>
                        <td class="col-database">CSV</td>
                        <td class="col-genomes">Heavy metal, biocide, stress response genes</td>
                    </tr>
                    <tr>
                        <td class="col-gene"><strong>acinetobacter_gene_categories.csv</strong></td>
                        <td class="col-category">Gene category analysis</td>
                        <td class="col-database">CSV</td>
                        <td class="col-genomes">Genes by resistance mechanism</td>
                    </tr>
                    <tr>
                        <td class="col-gene"><strong>acinetobacter_patterns.csv</strong></td>
                        <td class="col-category">Pattern discovery results</td>
                        <td class="col-database">CSV</td>
                        <td class="col-genomes">Cross-genome patterns and associations</td>
                    </tr>
                    <tr>
                        <td class="col-gene"><strong>acinetobacter_database_coverage.csv</strong></td>
                        <td class="col-category">Database performance</td>
                        <td class="col-database">CSV</td>
                        <td class="col-genomes">Database coverage and statistics</td>
                    </tr>
                </tbody>
            </table>
        </div>
        """


class GeniusUltimateReporter:
    """MASTER CLASS: Generates ultimate gene-centric reports for A. baumannii"""
    
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "GENIUS_ACINETOBACTER_ULTIMATE_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize components
        self.parser = UltimateHTMLParser()
        self.analyzer = UltimateDataAnalyzer()
        self.html_generator = UltimateHTMLGenerator(self.analyzer)
        
        # Metadata
        self.metadata = {
            "tool_name": "GENIUS Acinetobacter baumannii Ultimate Reporter",
            "version": "1.0.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "pathogen": "Acinetobacter baumannii",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }
    
    def find_html_files(self) -> Dict[str, List[Path]]:
        """Find all HTML report files - SPECIFICALLY CHECK FOR PLASMIDFINDER"""
        print("🔍 Searching for AcinetoScope HTML reports...")
        
        html_files = {
            'pasteur_mlst': [],
            'oxford_mlst': [],
            'kaptive': [],
            'amrfinder': [],
            'abricate': defaultdict(list),
            'plasmidfinder': []  # NEW: Separate list for plasmidfinder
        }
        
        # First pass: collect all HTML files
        all_html_files = list(self.input_dir.glob("**/*.html"))
        
        if not all_html_files:
            print("  ⚠️ No HTML files found in the directory!")
            return html_files
        
        print(f"  📁 Found {len(all_html_files)} HTML files")
        
        # Second pass: categorize files
        for html_file in all_html_files:
            filename = html_file.name.lower()
            
            # MLST files
            if 'pasteur' in filename and 'mlst' in filename:
                html_files['pasteur_mlst'].append(html_file)
            elif 'oxford' in filename and 'mlst' in filename:
                html_files['oxford_mlst'].append(html_file)
            
            # Kaptive files
            elif 'kaptive' in filename:
                html_files['kaptive'].append(html_file)
            
            # AMRfinder files
            elif 'amrfinder' in filename:
                html_files['amrfinder'].append(html_file)
            
            # PLASMIDFINDER files - SPECIFIC CHECK
            elif 'plasmidfinder' in filename:
                html_files['plasmidfinder'].append(html_file)
                print(f"    ✅ Found PlasmidFinder file: {html_file.name}")
            
            # ABRicate database files
            else:
                # Try to match with known database names (EXCLUDING plasmidfinder)
                matched = False
                for db_key in self.parser.db_name_mapping.keys():
                    if db_key in filename and 'plasmidfinder' not in db_key:
                        html_files['abricate'][db_key].append(html_file)
                        matched = True
                        break
                
                # If no match, check for common patterns (excluding plasmidfinder)
                if not matched:
                    if 'card' in filename and 'plasmidfinder' not in filename:
                        html_files['abricate']['acineto_card'].append(html_file)
                    elif 'resfinder' in filename and 'plasmidfinder' not in filename:
                        html_files['abricate']['acineto_resfinder'].append(html_file)
                    elif 'argannot' in filename and 'plasmidfinder' not in filename:
                        html_files['abricate']['acineto_argannot'].append(html_file)
                    elif 'vfdb' in filename and 'plasmidfinder' not in filename:
                        html_files['abricate']['acineto_vfdb'].append(html_file)
                    elif 'victors' in filename and 'plasmidfinder' not in filename:
                        html_files['abricate']['acineto_victors'].append(html_file)
                    elif 'ecoli_vf' in filename and 'plasmidfinder' not in filename:
                        html_files['abricate']['acineto_ecoli_vf'].append(html_file)
                    elif 'megares' in filename and 'plasmidfinder' not in filename:
                        html_files['abricate']['acineto_megares'].append(html_file)
                    elif 'bacmet2' in filename and 'plasmidfinder' not in filename:
                        html_files['abricate']['acineto_bacmet2'].append(html_file)
                    elif 'ecoh' in filename and 'plasmidfinder' not in filename:
                        html_files['abricate']['acineto_ecoh'].append(html_file)
                    elif 'ncbi' in filename and 'plasmidfinder' not in filename:
                        html_files['abricate']['acineto_ncbi'].append(html_file)
        
        # Print findings
        print("\n📊 File Discovery Summary:")
        print(f"  ✅ Pasteur MLST: {len(html_files['pasteur_mlst'])} files")
        print(f"  ✅ Oxford MLST: {len(html_files['oxford_mlst'])} files")
        print(f"  ✅ Kaptive: {len(html_files['kaptive'])} files")
        print(f"  ✅ AMRfinder: {len(html_files['amrfinder'])} files")
        print(f"  ✅ PlasmidFinder: {len(html_files['plasmidfinder'])} files")
        
        # Print ABRicate database findings (excluding plasmidfinder)
        print(f"\n🔬 ABRicate Databases Found:")
        abricate_count = 0
        for db_key, files in html_files['abricate'].items():
            if files:
                db_name = self.parser.db_name_mapping.get(db_key, db_key)
                print(f"  ✅ {db_name.upper()}: {len(files)} files")
                abricate_count += len(files)
        
        # Check for missing databases
        expected_dbs = ['acineto_card', 'acineto_resfinder', 'acineto_argannot', 
                       'acineto_vfdb', 'acineto_victors', 'acineto_ecoli_vf',
                       'acineto_megares', 'acineto_bacmet2']
        
        missing_dbs = []
        for db in expected_dbs:
            if not html_files['abricate'][db]:
                missing_dbs.append(db.replace('acineto_', ''))
        
        if missing_dbs:
            print(f"\n⚠️  Missing Databases (OK if no hits): {', '.join(missing_dbs)}")
        
        return html_files
    
    def integrate_all_data(self, html_files: Dict[str, List[Path]]) -> Dict[str, Any]:
        """Integrate data from all HTML reports - FIXED FOR PLASMIDFINDER"""
        print("\n🔗 Integrating data from all reports...")
        
        integrated_data = {
            'metadata': self.metadata,
            'samples': {},
            'patterns': {},
            'gene_centric': {},
            'plasmid_analysis': {},  # Initialize empty plasmid analysis
            'parsing_summary': {}
        }
        
        # Parse all reports
        parsing_summary = {}
        
        # 1. Parse MLST reports
        pasteur_mlst_data = {}
        if html_files['pasteur_mlst']:
            pasteur_mlst_data = self.parser.parse_mlst_report(html_files['pasteur_mlst'][0], "pasteur")
            parsing_summary['pasteur_mlst'] = len(pasteur_mlst_data)
        else:
            print("  ⚠️ No Pasteur MLST file found")
        
        oxford_mlst_data = {}
        if html_files['oxford_mlst']:
            oxford_mlst_data = self.parser.parse_mlst_report(html_files['oxford_mlst'][0], "oxford")
            parsing_summary['oxford_mlst'] = len(oxford_mlst_data)
        else:
            print("  ⚠️ No Oxford MLST file found")
        
        # 2. Parse Kaptive report
        kaptive_data = {}
        if html_files['kaptive']:
            kaptive_data = self.parser.parse_kaptive_report(html_files['kaptive'][0])
            parsing_summary['kaptive'] = len(kaptive_data)
        else:
            print("  ⚠️ No Kaptive file found")
        
        # 3. Collect all unique samples from typing data
        all_samples = set(pasteur_mlst_data.keys())
        
        # Add samples from other sources that might have different naming
        for sample in oxford_mlst_data.keys():
            all_samples.add(sample)
        for sample in kaptive_data.keys():
            all_samples.add(sample)
        
        # Store for gene frequency calculations
        total_samples = len(all_samples)
        print(f"\n📊 Found {total_samples} unique samples from typing data")
        
        # 4. Parse AMRfinder report (needs total_samples for percentages)
        amr_by_sample, amr_gene_freq = {}, {}
        if html_files['amrfinder']:
            amr_by_sample, amr_gene_freq = self.parser.parse_amrfinder_report(
                html_files['amrfinder'][0], total_samples
            )
            parsing_summary['amrfinder'] = {
                'samples': len(amr_by_sample),
                'genes': len(amr_gene_freq)
            }
            # Add AMRfinder samples (should match with normalized names)
            for sample in amr_by_sample.keys():
                all_samples.add(sample)
        else:
            print("  ⚠️ No AMRfinder file found")
        
        # 5. Parse ALL ABRicate databases (EXCLUDING plasmidfinder)
        abricate_data = {}
        abricate_gene_freq = {}
        abricate_summary = {}
        
        for db_key, files in html_files['abricate'].items():
            if files:
                db_name = self.parser.db_name_mapping.get(db_key, db_key)
                genes_by_sample, gene_freq = self.parser.parse_abricate_database_report(
                    files[0], total_samples
                )
                
                if genes_by_sample or gene_freq:
                    abricate_data[db_name] = genes_by_sample
                    abricate_gene_freq[db_name] = gene_freq
                    abricate_summary[db_name] = {
                        'samples': len(genes_by_sample),
                        'genes': len(gene_freq)
                    }
                    # Add ABRicate samples (should match with normalized names)
                    for sample in genes_by_sample.keys():
                        all_samples.add(sample)
        
        # 6. Parse PLASMIDFINDER specifically - ONLY if file exists
        plasmid_by_sample = {}
        plasmid_gene_freq = {}
        if html_files['plasmidfinder']:
            print("\n🧬 Processing PlasmidFinder report...")
            plasmid_by_sample, plasmid_gene_freq = self.parser.parse_plasmidfinder_report(
                html_files['plasmidfinder'][0], total_samples
            )
            parsing_summary['plasmidfinder'] = {
                'samples': len(plasmid_by_sample),
                'genes': len(plasmid_gene_freq)
            }
            # Add plasmidfinder samples
            for sample in plasmid_by_sample.keys():
                all_samples.add(sample)
            print(f"    ✅ PlasmidFinder: {len(plasmid_by_sample)} samples, {len(plasmid_gene_freq)} plasmid markers")
        else:
            print("  ⚠️ No PlasmidFinder file found - plasmid section will not be created")
        
        # Update total samples with all sources
        total_samples = len(all_samples)
        all_samples = sorted(list(all_samples))
        
        if not all_samples:
            print("❌ No samples found in any report!")
            return {}
        
        print(f"\n📊 Final count: {total_samples} unique samples across all data sources")
        
        # 7. Integrate data for each sample
        samples_with_typing = 0
        samples_with_amr = 0
        samples_with_virulence = 0
        samples_with_environmental = 0
        samples_with_plasmids = 0
        
        for sample in all_samples:
            # Get typing data using normalized sample names
            pasteur_info = pasteur_mlst_data.get(sample, {
                'ST': 'ND', 
                'International_Clone': 'Unknown', 
                'Allele_Profile': 'ND', 
                'Scheme': 'pasteur'
            })
            
            oxford_info = oxford_mlst_data.get(sample, {
                'ST': 'ND', 
                'International_Clone': 'Unknown', 
                'Allele_Profile': 'ND', 
                'Scheme': 'oxford'
            })
            
            kaptive_info = kaptive_data.get(sample, {
                'K_Locus': 'ND', 
                'O_Locus': 'ND', 
                'Capsule_Type': 'ND'
            })
            
            # Get AMR genes from AMRfinder
            amr_genes = amr_by_sample.get(sample, [])
            
            # Get virulence genes from virulence databases
            virulence_genes = []
            for db_name in ['vfdb', 'victors', 'ecoli_vf']:
                if db_name in abricate_data:
                    virulence_genes.extend(abricate_data[db_name].get(sample, []))
            
            # Get environmental genes
            environmental_genes = []
            if 'bacmet2' in abricate_data:
                environmental_genes.extend(abricate_data['bacmet2'].get(sample, []))
            
            # Get plasmid genes from PlasmidFinder (ONLY if file exists)
            plasmid_genes = []
            if plasmid_by_sample:  # Only if plasmidfinder file was found and parsed
                plasmid_genes = plasmid_by_sample.get(sample, [])
            
            # Get other genes from other databases
            other_genes = []
            for db_name, db_data in abricate_data.items():
                if db_name not in ['vfdb', 'victors', 'ecoli_vf', 'bacmet2', 'plasmidfinder']:
                    other_genes.extend(db_data.get(sample, []))
            
            # Count samples with data
            if pasteur_info['ST'] != 'ND' or oxford_info['ST'] != 'ND' or kaptive_info['K_Locus'] != 'ND':
                samples_with_typing += 1
            if amr_genes:
                samples_with_amr += 1
            if virulence_genes:
                samples_with_virulence += 1
            if environmental_genes:
                samples_with_environmental += 1
            if plasmid_genes:  # Only count if plasmidfinder file exists
                samples_with_plasmids += 1
            
            # Create sample data
            sample_data = {
                'pasteur_mlst': pasteur_info,
                'oxford_mlst': oxford_info,
                'kaptive': kaptive_info,
                'amr_genes': list(set(amr_genes)),
                'virulence_genes': list(set(virulence_genes)),
                'environmental_genes': list(set(environmental_genes)),
                'plasmid_genes': list(set(plasmid_genes)),  # Store plasmid genes
                'other_genes': list(set(other_genes))
            }
            
            integrated_data['samples'][sample] = sample_data
        
        # 8. Store gene frequencies 
            integrated_data['gene_frequencies'] = {
                'amrfinder': amr_gene_freq,
                'abricate': abricate_gene_freq
        }
        
        # IMPORTANT: Store plasmidfinder data separately so create_plasmid_analysis can find it
        if plasmid_gene_freq:
            integrated_data['gene_frequencies']['plasmidfinder'] = plasmid_gene_freq
            integrated_data['plasmidfinder_data'] = {
                'by_sample': plasmid_by_sample,
                'gene_freq': plasmid_gene_freq
            }
        
        # 9. Process gene-centric data and patterns
        print("\n🧠 Processing gene-centric analysis...")
        integrated_data['gene_centric'] = self.analyzer.create_gene_centric_tables(
            integrated_data, total_samples
        )
        integrated_data['patterns'] = self.analyzer.create_cross_genome_patterns(
            integrated_data, total_samples
        )
        
        # 10. Process plasmid analysis ONLY if plasmidfinder file exists AND has data
        if plasmid_gene_freq:  # Check if we have plasmid gene frequency data
            print("\n🧬 Processing plasmid analysis...")
            integrated_data['plasmid_analysis'] = self.analyzer.create_plasmid_analysis(integrated_data, total_samples)
            print(f"    ✅ Plasmid analysis created with {len(integrated_data['plasmid_analysis'].get('plasmid_databases', {}).get('plasmidfinder', []))} markers")
        else:
            print("\n⚠️ Skipping plasmid analysis - no PlasmidFinder data found")
            integrated_data['plasmid_analysis'] = {}  # Empty dict

        # 11. Store parsing summary
        integrated_data['parsing_summary'] = {
            'total_samples': total_samples,
            'samples_with_typing': samples_with_typing,
            'samples_with_amr': samples_with_amr,
            'samples_with_virulence': samples_with_virulence,
            'samples_with_environmental': samples_with_environmental,
            'samples_with_plasmids': samples_with_plasmids,
            'databases_parsed': len(abricate_summary) + (1 if amr_by_sample else 0) + (1 if plasmid_by_sample else 0),
            'abricate_summary': abricate_summary,
            **parsing_summary
        }
        
        # Print integration summary
        print("\n📈 Integration Summary:")
        print(f"  ✅ Total samples integrated: {total_samples}")
        print(f"  ✅ Samples with typing data: {samples_with_typing} ({samples_with_typing/total_samples*100:.1f}%)")
        print(f"  ✅ Samples with AMR genes: {samples_with_amr} ({samples_with_amr/total_samples*100:.1f}%)")
        print(f"  ✅ Samples with virulence genes: {samples_with_virulence} ({samples_with_virulence/total_samples*100:.1f}%)")
        print(f"  ✅ Samples with environmental markers: {samples_with_environmental} ({samples_with_environmental/total_samples*100:.1f}%)")
        if plasmid_by_sample:
            print(f"  ✅ Samples with plasmid markers: {samples_with_plasmids} ({samples_with_plasmids/total_samples*100:.1f}%)")
        print(f"  ✅ Databases successfully parsed: {len(abricate_summary) + (1 if amr_by_sample else 0) + (1 if plasmid_by_sample else 0)}")
        
        return integrated_data 
    
    def generate_json_report(self, integrated_data: Dict[str, Any]) -> Path:
        """Generate comprehensive JSON report"""
        print("\n📝 Generating JSON report...")
        
        output_file = self.output_dir / "genius_acinetobacter_ultimate_report.json"
        
        # Create serializable copy
        def make_serializable(obj):
            if obj is None:
                return None
            elif isinstance(obj, (str, int, float, bool)):
                return obj
            elif isinstance(obj, (list, tuple)):
                return [make_serializable(item) for item in obj]
            elif isinstance(obj, dict):
                return {str(k): make_serializable(v) for k, v in obj.items()}
            elif isinstance(obj, set):
                return [make_serializable(item) for item in obj]
            elif isinstance(obj, (Counter, defaultdict)):
                return {str(k): make_serializable(v) for k, v in dict(obj).items()}
            elif isinstance(obj, Path):
                return str(obj)
            elif hasattr(obj, 'isoformat'):
                return obj.isoformat()
            else:
                try:
                    return str(obj)
                except:
                    return None
        
        # Create serializable data
        serializable_data = make_serializable(integrated_data)
        
        # Write to file
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(serializable_data, f, indent=2, ensure_ascii=False)
        
        print(f"    ✅ JSON report saved: {output_file}")
        return output_file
    
    def generate_csv_reports(self, integrated_data: Dict[str, Any]):
        """Generate multiple CSV reports for A. baumannii - WITH COMBINED FREQUENCY"""
        print("\n📊 Generating CSV reports...")
        
        total_samples = len(integrated_data['samples'])
        
        # 1. Sample summary
        samples_data = []
        for sample, data in integrated_data['samples'].items():
            row = {
                'Sample': sample,
                'Pasteur_ST': data['pasteur_mlst']['ST'],
                'International_Clone': data['pasteur_mlst']['International_Clone'],
                'Oxford_ST': data['oxford_mlst']['ST'],
                'K_Locus': data['kaptive']['K_Locus'],
                'O_Locus': data['kaptive']['O_Locus'],
                'Capsule_Type': data['kaptive']['Capsule_Type'],
                'AMR_Gene_Count': len(data['amr_genes']),
                'Virulence_Gene_Count': len(data['virulence_genes']),
                'Environmental_Gene_Count': len(data['environmental_genes']),
                'Plasmid_Gene_Count': len(data.get('plasmid_genes', [])),  
                'Other_Gene_Count': len(data['other_genes']),
                'Total_Gene_Count': len(data['amr_genes']) + len(data['virulence_genes']) + len(data['environmental_genes']) + len(data.get('plasmid_genes', [])) + len(data['other_genes'])
            }
            samples_data.append(row)
        
        df_samples = pd.DataFrame(samples_data)
        samples_file = self.output_dir / "acinetobacter_samples.csv"
        df_samples.to_csv(samples_file, index=False)
        print(f"    ✅ Sample overview: {len(samples_data)} samples")
        
        # 2. AMR genes (gene-centric)
        amr_data = []
        gene_centric = integrated_data.get('gene_centric', {})
        
        for db_name, genes in gene_centric.get('amr_databases', {}).items():
            for gene_info in genes:
                amr_data.append({
                    'Gene': gene_info['gene'],
                    'Category': gene_info['category'],
                    'Database': gene_info['database'],
                    'Frequency': gene_info.get('frequency_display', f"{gene_info['count']} ({gene_info.get('percentage', 0):.1f}%)"),
                    'Count': gene_info['count'],
                    'Percentage': round(gene_info.get('percentage', 0), 2),
                    'Risk_Level': gene_info.get('risk_level', 'Standard'),
                    'Genomes': ';'.join(gene_info.get('genomes', []))
                })
        
        if amr_data:
            df_amr = pd.DataFrame(amr_data)
            amr_file = self.output_dir / "acinetobacter_amr_genes.csv"
            df_amr.to_csv(amr_file, index=False)
            print(f"    ✅ AMR genes: {len(amr_data)} genes")
        
        # 3. Virulence genes (gene-centric)
        virulence_data = []
        for db_name, genes in gene_centric.get('virulence_databases', {}).items():
            for gene_info in genes:
                virulence_data.append({
                    'Gene': gene_info['gene'],
                    'Category': gene_info['category'],
                    'Database': gene_info['database'],
                    'Frequency': gene_info.get('frequency_display', f"{gene_info['count']} ({gene_info.get('percentage', 0):.1f}%)"),
                    'Count': gene_info['count'],
                    'Percentage': round(gene_info.get('percentage', 0), 2),
                    'Genomes': ';'.join(gene_info.get('genomes', []))
                })
        
        if virulence_data:
            df_virulence = pd.DataFrame(virulence_data)
            virulence_file = self.output_dir / "acinetobacter_virulence_genes.csv"
            df_virulence.to_csv(virulence_file, index=False)
            print(f"    ✅ Virulence genes: {len(virulence_data)} genes")
        
        # 4. Environmental markers
        environmental_data = []
        environmental_summary = gene_centric.get('environmental_summary', {})
        
        for category, data in environmental_summary.items():
            if 'genes' in data:
                for gene_info in data['genes']:
                    environmental_data.append({
                        'Category': category,
                        'Gene': gene_info['gene'],
                        'Database': gene_info['database'],
                        'Frequency': gene_info.get('frequency_display', f"{gene_info['count']} ({gene_info.get('percentage', 0):.1f}%)"),
                        'Count': gene_info['count'],
                        'Percentage': round(gene_info.get('percentage', 0), 2),
                        'Genomes': ';'.join(gene_info.get('genomes', []))
                    })
        
        if environmental_data:
            df_environmental = pd.DataFrame(environmental_data)
            environmental_file = self.output_dir / "acinetobacter_environmental_markers.csv"
            df_environmental.to_csv(environmental_file, index=False)
            print(f"    ✅ Environmental markers: {len(environmental_data)} genes")
        
        # 5. Gene categories
        category_data = []
        for category, genes in gene_centric.get('gene_categories', {}).items():
            for gene_info in genes:
                category_data.append({
                    'Category': category,
                    'Gene': gene_info['gene'],
                    'Database': gene_info['database'],
                    'Frequency': gene_info.get('frequency_display', f"{gene_info['count']} ({gene_info.get('percentage', 0):.1f}%)"),
                    'Count': gene_info['count'],
                    'Percentage': round(gene_info.get('percentage', 0), 2),
                    'Genomes': ';'.join(gene_info.get('genomes', []))
                })
        
        if category_data:
            df_categories = pd.DataFrame(category_data)
            categories_file = self.output_dir / "acinetobacter_gene_categories.csv"
            df_categories.to_csv(categories_file, index=False)
            print(f"    ✅ Gene categories: {len(category_data)} entries")
        
        # 6. Pattern discovery
        pattern_data = []
        patterns = integrated_data['patterns']
        
        # High-risk combinations
        for combo in patterns.get('high_risk_combinations', []):
            pattern_data.append({
                'Pattern_Type': 'High_Risk_Combination',
                'Sample': combo['sample'],
                'Pasteur_ST': combo['pasteur_st'],
                'International_Clone': combo['international_clone'],
                'K_Locus': combo['k_locus'],
                'Capsule_Type': combo['capsule_type'],
                'Carbapenemases': ';'.join(combo['carbapenemases']),
                'Colistin_Resistance_Genes': ';'.join(combo['colistin_resistance']) if combo['colistin_resistance'] else '',
                'Tigecycline_Resistance_Genes': ';'.join(combo['tigecycline_resistance']) if combo['tigecycline_resistance'] else '',
                'Environmental_Markers': ';'.join(combo['environmental_markers']) if combo['environmental_markers'] else ''
            })
        
        # MDR patterns
        for mdr in patterns.get('mdr_patterns', []):
            pattern_data.append({
                'Pattern_Type': 'MDR_Pattern',
                'Sample': mdr['sample'],
                'Pasteur_ST': mdr['pasteur_st'],
                'International_Clone': mdr['international_clone'],
                'Resistance_Classes': mdr['resistance_types'],
                'Carbapenemases': ';'.join(mdr['carbapenemases']) if mdr['carbapenemases'] else '',
                'ESBLs': ';'.join(mdr['esbls']) if mdr['esbls'] else '',
                'Colistin_Resistance': ';'.join(mdr['colistin_resistance']) if mdr['colistin_resistance'] else '',
                'Tigecycline_Resistance': ';'.join(mdr['tigecycline_resistance']) if mdr['tigecycline_resistance'] else '',
                'Environmental_Markers': ';'.join(mdr['environmental_markers']) if mdr['environmental_markers'] else ''
            })
        
        # ST-K Locus combinations
        for combo, samples in patterns.get('st_k_locus_combinations', {}).items():
            pattern_data.append({
                'Pattern_Type': 'ST_K_Locus_Association',
                'Combination': combo,
                'Samples': ';'.join(samples),
                'Count': len(samples)
            })
        
        # Environmental patterns
        for env_combo, samples in patterns.get('environmental_patterns', {}).items():
            pattern_data.append({
                'Pattern_Type': 'Environmental_Pattern',
                'Combination': ';'.join(env_combo),
                'Samples': ';'.join(samples),
                'Count': len(samples)
            })
        
        if pattern_data:
            df_patterns = pd.DataFrame(pattern_data)
            patterns_file = self.output_dir / "acinetobacter_patterns.csv"
            df_patterns.to_csv(patterns_file, index=False)
            print(f"    ✅ Patterns: {len(pattern_data)} patterns")
        
        # 7. Database coverage
        coverage_data = []
        database_coverage = patterns.get('database_coverage', {})
        database_stats = gene_centric.get('database_stats', {})
        
        for db_name, coverage in database_coverage.items():
            stats = database_stats.get(db_name, {})
            
            coverage_data.append({
                'Database': db_name.upper(),
                'Coverage': coverage.get('coverage_display', f"{coverage['samples_with_hits']} ({coverage['coverage_percentage']}%)"),
                'Samples_With_Hits': coverage['samples_with_hits'],
                'Total_Samples': coverage['total_samples'],
                'Coverage_Percentage': coverage['coverage_percentage'],
                'Unique_Genes': stats.get('total_genes', 0),
                'Total_Occurrences': stats.get('total_occurrences', 0),
                'Critical_Genes': stats.get('critical_genes', 0),
                'Environmental_Genes': stats.get('environmental_genes', 0)
            })
        
        if coverage_data:
            df_coverage = pd.DataFrame(coverage_data)
            coverage_file = self.output_dir / "acinetobacter_database_coverage.csv"
            df_coverage.to_csv(coverage_file, index=False)
            print(f"    ✅ Database coverage: {len(coverage_data)} databases")

        # 8. Plasmid analysis - FIXED: Check for plasmid data correctly
        plasmid_data = []
        plasmid_analysis = integrated_data.get('plasmid_analysis', {})
        
        # Check if we have plasmid data - fixed condition
        if plasmid_analysis.get('plasmid_databases'):
            for db_name, plasmids in plasmid_analysis.get('plasmid_databases', {}).items():
                for plasmid_info in plasmids:
                    plasmid_data.append({
                        'Plasmid_Marker': plasmid_info['plasmid_marker'],
                        'Full_Name': plasmid_info.get('full_name', plasmid_info['plasmid_marker']),
                        'Category': plasmid_info['category'],
                        'Database': plasmid_info['database'],
                        'Frequency': plasmid_info.get('frequency_display', f"{plasmid_info['count']} ({plasmid_info.get('percentage', 0):.1f}%)"),
                        'Count': plasmid_info['count'],
                        'Percentage': round(plasmid_info.get('percentage', 0), 2),
                        'Genomes': ';'.join(plasmid_info.get('genomes', []))
                    })

            if plasmid_data:
                df_plasmid = pd.DataFrame(plasmid_data)
                plasmid_file = self.output_dir / "acinetobacter_plasmid_analysis.csv"
                df_plasmid.to_csv(plasmid_file, index=False)
                print(f"    ✅ Plasmid analysis: {len(plasmid_data)} markers")
            else:
                print(f"    ⚠️ No plasmid markers found in plasmid analysis")
        else:
            print(f"    ⚠️ No plasmid analysis data available - skipping plasmid CSV")     
    
    def run(self):
        """Run the complete analysis for A. baumannii"""
        print("=" * 80)
        print("🧠 GENIUS ACINETOBACTER BAUMANNII ULTIMATE REPORTER v1.0.0")
        print("=" * 80)
        print(f"📁 Input directory: {self.input_dir}")
        print(f"🦠 Pathogen: Acinetobacter baumannii (MDR/XDR Focus)")
        print(f"🎯 Analysis: Gene-centric with environmental co-selection markers")
        print("=" * 80)
        
        # Find HTML files
        html_files = self.find_html_files()
        
        if not any(html_files.values()):
            print("\n❌ No HTML report files found in the directory!")
            print("   Please ensure HTML files are in the correct location.")
            return False
        
        # Check for critical files
        critical_files_found = False
        if html_files['pasteur_mlst'] or html_files['oxford_mlst']:
            critical_files_found = True
        
        if not critical_files_found:
            print("\n⚠️  WARNING: No MLST files found. Analysis will proceed without typing data.")
        
        # Integrate data
        print("\n" + "=" * 80)
        print("🔗 INTEGRATING DATA FROM ALL REPORTS")
        print("=" * 80)
        
        integrated_data = self.integrate_all_data(html_files)
        if not integrated_data:
            print("\n❌ Data integration failed!")
            return False
        
        # Generate reports
        print("\n" + "=" * 80)
        print("📊 GENERATING ULTIMATE REPORTS FOR A. BAUMANNII")
        print("=" * 80)
        
        # Generate JSON
        json_file = self.generate_json_report(integrated_data)
        
        # Generate CSV
        self.generate_csv_reports(integrated_data)
        
        # Generate HTML
        html_file = self.html_generator.generate_main_report(integrated_data, self.output_dir)
        
        # Print comprehensive summary
        total_samples = len(integrated_data['samples'])
        patterns = integrated_data['patterns']
        high_risk = len(patterns.get('high_risk_combinations', []))
        mdr_patterns = len(patterns.get('mdr_patterns', []))
        gene_centric = integrated_data['gene_centric']
        parsing_summary = integrated_data.get('parsing_summary', {})
        
        # Count genes by category
        total_amr_genes = sum(len(db) for db in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(db) for db in gene_centric.get('virulence_databases', {}).values())
        total_environmental_genes = sum(len(db) for db in gene_centric.get('environmental_databases', {}).values())
        
        # Check if we have plasmid data - FIXED LOGIC
        plasmid_analysis = integrated_data.get('plasmid_analysis', {})
        has_plasmid_data = False
        total_plasmid_genes = 0
        samples_with_plasmids = 0
        
        # Check for plasmid data more thoroughly
        if plasmid_analysis:
            plasmid_databases = plasmid_analysis.get('plasmid_databases', {})
            if plasmid_databases:
                has_plasmid_data = True
                total_plasmid_genes = sum(len(db) for db in plasmid_databases.values())
                plasmid_stats = plasmid_analysis.get('plasmid_summary_stats', {})
                samples_with_plasmids = plasmid_stats.get('samples_with_plasmids', 0)
        
        # Count carbapenemases and environmental markers
        carbapenemase_count = 0
        environmental_marker_count = 0
        for db_name, genes in gene_centric.get('amr_databases', {}).items():
            for gene_data in genes:
                if gene_data['category'] == 'Carbapenemases':
                    carbapenemase_count += 1
                if gene_data['category'] in ['Environmental Co-Selection', 'BACMET2 Markers']:
                    environmental_marker_count += 1
        
        if 'bacmet2' in gene_centric.get('environmental_databases', {}):
            environmental_marker_count += len(gene_centric['environmental_databases']['bacmet2'])
        
        print("\n" + "=" * 80)
        print("✅ ULTIMATE ANALYSIS COMPLETE FOR A. BAUMANNII!")
        print("=" * 80)
        print(f"📁 Output directory: {self.output_dir}")
        print(f"\n📄 FILES GENERATED:")
        print(f"   • genius_acinetobacter_ultimate_report.html (Interactive report)")
        print(f"   • genius_acinetobacter_ultimate_report.json (Complete data)")
        print(f"   • acinetobacter_samples.csv (Sample overview with counts)")
        print(f"   • acinetobacter_amr_genes.csv (AMR genes with frequency 'count (percentage%)')")
        print(f"   • acinetobacter_virulence_genes.csv (Virulence genes with frequency 'count (percentage%)')")
        print(f"   • acinetobacter_environmental_markers.csv (Environmental resistance markers)")
        print(f"   • acinetobacter_gene_categories.csv (Resistance mechanism analysis)")
        print(f"   • acinetobacter_patterns.csv (MDR/XDR patterns)")
        print(f"   • acinetobacter_database_coverage.csv (Database performance)")
        if has_plasmid_data:
            print(f"   • acinetobacter_plasmid_analysis.csv (Plasmid marker analysis)")
        
        print(f"\n🔬 KEY IMPROVEMENTS:")
        print(f"   • ✅ Environmental resistance & co-selection markers added")
        print(f"   • ✅ Frequency displayed as 'count (percentage%)' in single column")
        print(f"   • ✅ All tables now fully scrollable (no truncation)")
        print(f"   • ✅ BACMET2 and VICTORS database markers categorized")
        print(f"   • ✅ Heavy metal, biocide, stress response tracking")
        print(f"   • ✅ Environmental co-selection pattern discovery")
        if has_plasmid_data:
            print(f"   • ✅ Plasmid analysis from PlasmidFinder")
        
        print(f"\n📈 ANALYSIS SUMMARY:")
        print(f"   • Total samples analyzed: {total_samples}")
        print(f"   • Samples with typing data: {parsing_summary.get('samples_with_typing', 0)} ({parsing_summary.get('samples_with_typing', 0)/total_samples*100:.1f}%)")
        print(f"   • Samples with AMR genes: {parsing_summary.get('samples_with_amr', 0)} ({parsing_summary.get('samples_with_amr', 0)/total_samples*100:.1f}%)")
        print(f"   • Samples with virulence genes: {parsing_summary.get('samples_with_virulence', 0)} ({parsing_summary.get('samples_with_virulence', 0)/total_samples*100:.1f}%)")
        print(f"   • Samples with environmental markers: {parsing_summary.get('samples_with_environmental', 0)} ({parsing_summary.get('samples_with_environmental', 0)/total_samples*100:.1f}%)")
        if has_plasmid_data:
            print(f"   • Samples with plasmid markers: {samples_with_plasmids} ({samples_with_plasmids/total_samples*100:.1f}%)")
            print(f"   • Plasmid markers detected: {total_plasmid_genes}")
        print(f"   • Carbapenemase genes detected: {carbapenemase_count}")
        print(f"   • Environmental resistance markers: {environmental_marker_count}")
        print(f"   • Total AMR genes: {total_amr_genes}")
        print(f"   • Total virulence genes: {total_virulence_genes}")
        print(f"   • Total environmental genes: {total_environmental_genes}")
        print(f"   • High-risk combinations: {high_risk}")
        print(f"   • MDR/XDR patterns: {mdr_patterns}")
        print(f"   • Databases successfully parsed: {parsing_summary.get('databases_parsed', 0)}")
        
        print(f"\n🎯 NEXT STEPS:")
        print(f"   1. Open genius_acinetobacter_ultimate_report.html in your browser")
        print(f"   2. Check 'Environmental' tab for co-selection markers")
        print(f"   3. Review 'AMR Genes' tab - frequency shown as 'count (percentage%)'")
        if has_plasmid_data:
            print(f"   4. Check 'Plasmids' tab for plasmid marker analysis")
            step_offset = 1
        else:
            step_offset = 0
        print(f"   {5 + step_offset}. Check 'Pattern Discovery' for environmental co-selection patterns")
        print(f"   {6 + step_offset}. Use 'Gene Categories' to analyze resistance mechanisms")
        print(f"   {7 + step_offset}. Export data for surveillance reporting")
        
        print(f"\n⚠️  CRITICAL FINDINGS ALERT:")
        if carbapenemase_count > 0:
            print(f"   • 🔴 {carbapenemase_count} CARBAPENEMASE genes detected - WHO CRITICAL PRIORITY")
        if environmental_marker_count > 0:
            print(f"   • 🟡 {environmental_marker_count} ENVIRONMENTAL resistance markers - risk of co-selection")
        if high_risk > 0:
            print(f"   • 🔴 {high_risk} samples with carbapenemase + last-resort resistance")
        if mdr_patterns > 0:
            print(f"   • 🟡 {mdr_patterns} MDR/XDR A. baumannii isolates identified")
        if has_plasmid_data and samples_with_plasmids > 0:
            print(f"   • 🔵 {samples_with_plasmids} samples with plasmid markers - potential for horizontal gene transfer")
        
        if carbapenemase_count == 0 and high_risk == 0 and mdr_patterns == 0:
            print(f"   • ✅ No critical resistance patterns detected")
        
        print(f"\n📧 Report issues: brownbeckley94@gmail.com")
        print(f"🏫 University of Ghana Medical School - Department of Medical Biochemistry")
        print("\n" + "=" * 80)
        return True


def main():
    """Main function for A. baumannii reporter"""
    parser = argparse.ArgumentParser(
        description='GENIUS Acinetobacter baumannii Ultimate Reporter - Comprehensive Edition',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  python genius_acinetobacter_reporter.py -i /path/to/acinetoscope/html/reports
  
  python genius_acinetobacter_reporter.py -i ./summary_module -o ./my_reports

REQUIREMENTS:
  - AcinetoScope HTML report files in the input directory
  - Python 3.7+ with pandas and beautifulsoup4
  - Files should have standardized names (e.g., acineto_card_summary_report.html)

NEW FEATURES v1.0.0:
  ✅ Environmental resistance & co-selection markers (BACMET2/VICTORS)
  ✅ Frequency displayed as "count (percentage%)" in single column
  ✅ All tables fully scrollable (no truncation)
  ✅ Heavy metal, biocide, stress response gene tracking
  ✅ Environmental co-selection pattern discovery
  ✅ Improved categorization with environmental markers

AUTHOR: Brown Beckley <brownbeckley94@gmail.com>
AFFILIATION: University of Ghana Medical School
PATHOGEN: Acinetobacter baumannii (MDR/XDR focus with environmental co-selection)
        """
    )
    
    parser.add_argument('-i', '--input-dir', required=True,
                       help='Directory containing AcinetoScope HTML report files')
    parser.add_argument('-o', '--output-dir',
                       help='Custom output directory (default: input_dir/GENIUS_ACINETOBACTER_ULTIMATE_REPORTS)')
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    
    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        sys.exit(1)
    
    if not any(input_dir.glob("*.html")):
        print(f"❌ No HTML files found in: {input_dir}")
        print(f"   Please ensure HTML report files are in the directory.")
        sys.exit(1)
    
    # Create and run reporter
    try:
        reporter = GeniusUltimateReporter(input_dir)
        
        if args.output_dir:
            reporter.output_dir = Path(args.output_dir)
            reporter.output_dir.mkdir(parents=True, exist_ok=True)
        
        success = reporter.run()
        
        if not success:
            print("❌ Report generation failed!")
            sys.exit(1)
            
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        print(f"\n💡 TROUBLESHOOTING:")
        print(f"   1. Ensure all HTML files are in the input directory")
        print(f"   2. Check file permissions")
        print(f"   3. Verify Python packages are installed: pandas, beautifulsoup4")
        print(f"   4. Contact: brownbeckley94@gmail.com")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()