#!/usr/bin/env python3
"""
AcinetoScope - MLST Module for Acinetobacter baumannii
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2025-12-29
Need help? Reach out by mail!!!
"""

import os
import sys
import json
import glob
import argparse
import subprocess
import random
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import pandas as pd
from datetime import datetime

class AcinetoMLSTAnalyzer:
    def __init__(self, database_dir: Path, script_dir: Path):
        self.database_dir = database_dir
        self.script_dir = script_dir
        self.mlst_bin = script_dir / "mlst"
        
        # Verify mlst binary exists
        if not self.mlst_bin.exists():
            # Try to find mlst in PATH
            mlst_path = shutil.which("mlst")
            if mlst_path:
                self.mlst_bin = Path(mlst_path)
            else:
                raise FileNotFoundError(f"MLST binary not found at: {self.mlst_bin}")
        
        # Check for Excel support
        self.has_excel_support = self.check_excel_support()
        
        # Science quotes
        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
            {"text": "The universe is not required to be in harmony with human ambition.", "author": "Carl Sagan"},
        ]
        
        # A. baumannii MLST database names
        self.scheme_databases = {
            "oxford": "abaumannii",      # Oxford scheme
            "pasteur": "abaumannii_2"    # Pasteur scheme
        }
        
        # Scheme display names
        self.scheme_display_names = {
            "abaumannii": "OXFORD",
            "abaumannii_2": "PASTEUR"
        }
        
        # International Clone mapping
        self.ic_mapping = {
            "IC I": {"abaumannii_2": "1", "abaumannii": "231"},
            "IC II": {"abaumannii_2": "2", "abaumannii": "208"},
            "IC III": {"abaumannii_2": "3", "abaumannii": "452"},
            "IC IV": {"abaumannii_2": "4", "abaumannii": "195"},
        }

    def check_excel_support(self):
        """Check if Excel export is available"""
        try:
            import openpyxl
            return True
        except ImportError:
            return False

    def get_random_quote(self):
        """Get a random science quote"""
        return random.choice(self.science_quotes)
    
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns with proper handling"""
        # If it's a file, return it
        if os.path.isfile(input_path):
            return [Path(input_path)]
        
        # If it's a directory, search recursively
        if os.path.isdir(input_path):
            patterns = [
                "**/*.fna",
                "**/*.fasta", 
                "**/*.fa",
                "**/*.fn",
                "**/*.fna.gz",
                "**/*.fasta.gz",
                "**/*.fa.gz",
                "**/*.gb",
                "**/*.gbk",
                "**/*.gbff"
            ]
        else:
            # It's a pattern
            patterns = [input_path]
        
        fasta_files = []
        for pattern in patterns:
            try:
                # Handle both directory search and direct pattern
                if os.path.isdir(input_path):
                    matched_files = glob.glob(os.path.join(input_path, pattern), recursive=True)
                else:
                    matched_files = glob.glob(pattern, recursive=True)
                
                for file_path in matched_files:
                    path = Path(file_path)
                    if path.is_file():
                        # Check if it looks like a FASTA file
                        if self.is_fasta_file(path):
                            fasta_files.append(path)
            except Exception as e:
                print(f"⚠️ Warning: Pattern {pattern} failed: {e}")
                continue
        
        # Remove duplicates and sort
        fasta_files = sorted(list(set(fasta_files)))
        return fasta_files
    
    def is_fasta_file(self, file_path: Path) -> bool:
        """Check if file is likely a FASTA file"""
        # Check extension
        valid_extensions = {'.fna', '.fasta', '.fa', '.fn', '.gb', '.gbk', '.gbff'}
        if file_path.suffix.lower() in valid_extensions:
            return True
        
        # Check first line for FASTA header
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                return first_line.startswith('>')
        except:
            return False

    def run_mlst_single(self, input_file: Path, output_dir: Path, scheme: str = "abaumannii") -> Dict:
        """Run MLST analysis for a single file"""
        print(f"🔬 Processing: {input_file.name} with {scheme} scheme")
        
        # Create output directory
        sample_output_dir = output_dir / input_file.stem
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save raw MLST output
        raw_output_file = sample_output_dir / "mlst_raw_output.txt"
        
        # Run MLST command
        mlst_cmd = [
            str(self.mlst_bin),
            str(input_file),
            "--scheme", scheme,
            "--csv",
            "--nopath"
        ]
        
        try:
            print(f"   Running: {' '.join(mlst_cmd)}")
            result = subprocess.run(mlst_cmd, capture_output=True, text=True, check=True)
            
            # Save raw output
            with open(raw_output_file, 'w') as f:
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
            
            # Parse the CSV output - FIXED PARSING
            mlst_results = self.parse_mlst_csv_fixed(result.stdout, input_file.name, scheme)
            
            # Add lineage information
            lineage_info = self.get_lineage_info(mlst_results.get('st', 'UNKNOWN'), scheme)
            mlst_results.update(lineage_info)
            
            # Generate output files
            self.generate_output_files(mlst_results, sample_output_dir)
            
            st_display = mlst_results.get('st', 'UNKNOWN')
            print(f"✅ Completed: {input_file.name} -> ST{st_display}")
            return mlst_results
            
        except subprocess.CalledProcessError as e:
            print(f"❌ MLST failed for {input_file.name}: {e}")
            # Save error output
            with open(raw_output_file, 'w') as f:
                f.write(f"ERROR: {e}\n")
                f.write(f"STDOUT: {e.stdout}\n")
                f.write(f"STDERR: {e.stderr}\n")
            
            error_result = self.get_fallback_results(input_file.name, scheme)
            self.generate_output_files(error_result, sample_output_dir)
            return error_result

    def parse_mlst_csv_fixed(self, stdout: str, sample_name: str, scheme: str) -> Dict:
        """Parse MLST CSV output - FIXED VERSION"""
        lines = stdout.strip().split('\n')
        if not lines:
            return self.get_empty_results(sample_name, scheme)
        
        # Find the result line (usually the last non-empty line)
        result_line = None
        for line in reversed(lines):
            line = line.strip()
            if line and ',' in line:
                # Check if it's a result line (not a header or message)
                if not line.startswith('[') and not line.startswith('file'):
                    result_line = line
                    break
        
        if not result_line:
            return self.get_empty_results(sample_name, scheme)
        
        # Split by comma - handle quoted fields
        parts = []
        current_part = []
        in_quotes = False
        
        for char in result_line:
            if char == '"':
                in_quotes = not in_quotes
            elif char == ',' and not in_quotes:
                parts.append(''.join(current_part))
                current_part = []
            else:
                current_part.append(char)
        
        if current_part:
            parts.append(''.join(current_part))
        
        # Clean parts
        parts = [p.strip() for p in parts]
        
        # Expected format: filename, scheme, ST, allele1(allele_val), allele2(allele_val), ...
        if len(parts) < 3:
            return self.get_empty_results(sample_name, scheme)
        
        filename = parts[0]
        detected_scheme = parts[1]
        st_raw = parts[2]
        
        # Handle ST values: could be "-", "STxxx", or just "xxx"
        if st_raw == '-' or st_raw == '' or st_raw == '0':
            st = "UNKNOWN"
        elif st_raw.startswith('ST'):
            st = st_raw[2:]  # Remove "ST" prefix
        else:
            st = st_raw
        
        # Extract alleles
        alleles = {}
        allele_parts = []
        
        for i in range(3, len(parts)):
            allele_str = parts[i].strip()
            if not allele_str:
                continue
                
            # Parse allele in format "gene(allele)"
            if '(' in allele_str and ')' in allele_str:
                gene = allele_str.split('(')[0].strip()
                allele = allele_str.split('(')[1].rstrip(')').strip()
                alleles[gene] = allele
                allele_parts.append(f"{gene}({allele})")
            elif allele_str:  # Just gene name without allele
                alleles[allele_str] = "?"
                allele_parts.append(f"{allele_str}(?)")
        
        allele_profile = '-'.join(allele_parts) if allele_parts else ""
        
        # Determine confidence
        if st != "UNKNOWN" and st != "-" and st != "":
            confidence = "HIGH"
            mlst_assigned = True
        else:
            confidence = "LOW"
            mlst_assigned = False
        
        return {
            "sample": sample_name,
            "original_filename": filename,
            "st": st,
            "scheme": scheme,
            "scheme_display": self.scheme_display_names.get(scheme, scheme.upper()),
            "alleles": alleles,
            "allele_profile": allele_profile,
            "confidence": confidence,
            "mlst_assigned": mlst_assigned,
            "detected_scheme": detected_scheme
        }

    def get_lineage_info(self, st: str, scheme: str) -> Dict:
        """Get lineage information based on ST"""
        # Comprehensive database for A. baumannii
        lineage_db = {
            "abaumannii": {  # Oxford scheme
                '231': {
                    "international_clone": "IC I",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=231"
                },
                '208': {
                    "international_clone": "IC II",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=208"
                },
                '452': {
                    "international_clone": "IC III",
                    "clonal_complex": "CC3",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly South America",
                    "clinical_significance": "Third major global clone, associated with hospital-acquired infections",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-72"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=452"
                },
                '195': {
                    "international_clone": "IC IV",
                    "clonal_complex": "CC4",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Emerging global clone, increasing in prevalence",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=195"
                },
                '81': {
                    "international_clone": "IC V",
                    "clonal_complex": "CC10",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, Europe",
                    "clinical_significance": "Fifth international clone, often OXA-23 producer",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=81"
                },
                '619': {
                    "international_clone": "IC VI",
                    "clonal_complex": "CC109",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "South America, Europe",
                    "clinical_significance": "Sixth international clone, often associated with OXA-58",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-58", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=619"
                },
                '114': {
                    "international_clone": "IC VII",
                    "clonal_complex": "CC25",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Middle East, Asia",
                    "clinical_significance": "Seventh international clone, emerging threat",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-23", "NDM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=114"
                },
                '449': {
                    "international_clone": "IC VIII",
                    "clonal_complex": "CC164",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Europe, Middle East",
                    "clinical_significance": "Eighth international clone, often colistin-resistant",
                    "common_resistance": ["Colistin-resistant", "Carbapenem-resistant", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like", "mcr-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=449"
                }
            },
            "abaumannii_2": {  # Pasteur scheme - COMPREHENSIVE IC DATABASE
                '1': {
                    "international_clone": "IC I",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-66", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1"
                },
                '2': {
                    "international_clone": "IC II",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=2"
                },
                '3': {
                    "international_clone": "IC III",
                    "clonal_complex": "CC3",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly South America and Europe",
                    "clinical_significance": "Third major global clone, often associated with OXA-72 carbapenemase",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-72", "OXA-51-like", "ADC-73"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=3"
                },
                '4': {
                    "international_clone": "IC IV",
                    "clonal_complex": "CC4",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Europe",
                    "clinical_significance": "Emerging global clone, increasing in prevalence. Often associated with OXA-58",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-58", "OXA-51-like", "ADC-96"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=4"
                },
                '5': {
                    "international_clone": "IC V",
                    "clonal_complex": "CC5",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, Europe, Middle East",
                    "clinical_significance": "Fifth international clone, often OXA-23 producer. Emerging in Asian hospitals",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=5"
                },
                '6': {
                    "international_clone": "IC VI",
                    "clonal_complex": "CC6",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "South America, Europe",
                    "clinical_significance": "Sixth international clone, often associated with OXA-58 and colistin resistance",
                    "common_resistance": ["Carbapenem-resistant", "Colistin-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-58", "OXA-51-like", "mcr-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=6"
                },
                '7': {
                    "international_clone": "IC VII",
                    "clonal_complex": "CC7",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Middle East, Asia",
                    "clinical_significance": "Seventh international clone, often NDM-1 metallo-β-lactamase producer",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["NDM-1", "OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=7"
                },
                '8': {
                    "international_clone": "IC VIII",
                    "clonal_complex": "CC8",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Europe, Middle East",
                    "clinical_significance": "Eighth international clone, often pan-drug resistant with colistin resistance",
                    "common_resistance": ["Colistin-resistant", "Carbapenem-resistant", "XDR", "PDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like", "mcr-1", "NDM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=8"
                },
                '10': {
                    "international_clone": "IC IX",
                    "clonal_complex": "CC10",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, particularly China and India",
                    "clinical_significance": "Ninth international clone, emerging threat with OXA-23 and NDM co-production",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "NDM-1", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=10"
                },
                '15': {
                    "international_clone": "IC X",
                    "clonal_complex": "CC15",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Europe, North America",
                    "clinical_significance": "Tenth international clone, often associated with OXA-40/24 carbapenemase",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-40", "OXA-51-like", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=15"
                },
                '20': {
                    "international_clone": "IC XI",
                    "clonal_complex": "CC20",
                    "classification": "Regional Epidemic Clone",
                    "geographic_distribution": "Middle East, Mediterranean region",
                    "clinical_significance": "Eleventh international clone, often OXA-143 producer",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-143", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=20"
                },
                '25': {
                    "international_clone": "IC XII",
                    "clonal_complex": "CC25",
                    "classification": "Regional Epidemic Clone",
                    "geographic_distribution": "South America, Caribbean",
                    "clinical_significance": "Twelfth international clone, often associated with OXA-72 and OXA-58",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-72", "OXA-58", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=25"
                },
                '30': {
                    "international_clone": "IC XIII",
                    "clonal_complex": "CC30",
                    "classification": "Regional Epidemic Clone",
                    "geographic_distribution": "Africa, Middle East",
                    "clinical_significance": "Thirteenth international clone, emerging in African hospitals",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=30"
                },
                '40': {
                    "international_clone": "IC XIV",
                    "clonal_complex": "CC40",
                    "classification": "Regional Epidemic Clone",
                    "geographic_distribution": "Asia-Pacific region",
                    "clinical_significance": "Fourteenth international clone, often IMP-type metallo-β-lactamase producer",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["IMP-type", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=40"
                },
                '49': {
                    "international_clone": "IC XV",
                    "clonal_complex": "CC49",
                    "classification": "Regional Epidemic Clone",
                    "geographic_distribution": "Europe, particularly Mediterranean",
                    "clinical_significance": "Fifteenth international clone, often OXA-58 and PER-1 producer",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-58", "PER-1", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=49"
                },
                '60': {
                    "international_clone": "IC XVI",
                    "clonal_complex": "CC60",
                    "classification": "Regional Epidemic Clone",
                    "geographic_distribution": "Latin America",
                    "clinical_significance": "Sixteenth international clone, often KPC producer",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["KPC", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=60"
                },
                '78': {
                    "international_clone": "IC XVII",
                    "clonal_complex": "CC78",
                    "classification": "Regional Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly in ICU settings",
                    "clinical_significance": "Seventeenth international clone, often associated with ventilator-associated pneumonia",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "Colistin-resistant"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like", "mcr-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=78"
                },
                '85': {
                    "international_clone": "IC XVIII",
                    "clonal_complex": "CC85",
                    "classification": "Regional Epidemic Clone",
                    "geographic_distribution": "Middle East, North Africa",
                    "clinical_significance": "Eighteenth international clone, often VIM producer",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["VIM", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=85"
                },
                '92': {
                    "international_clone": "IC XIX",
                    "clonal_complex": "CC92",
                    "classification": "Regional Epidemic Clone",
                    "geographic_distribution": "Asia, particularly India and Pakistan",
                    "clinical_significance": "Nineteenth international clone, often NDM and OXA-23 co-producer",
                    "common_resistance": ["Carbapenem-resistant", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["NDM-1", "OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=92"
                },
                '103': {
                    "international_clone": "IC XX",
                    "clonal_complex": "CC103",
                    "classification": "Regional Epidemic Clone",
                    "geographic_distribution": "Europe, North America",
                    "clinical_significance": "Twentieth international clone, often associated with OXA-48-like enzymes",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-48-like", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=103"
                },
                # Additional notable STs in Pasteur scheme
                '79': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC79",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Widely distributed clone, often carbapenem-resistant",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=79"
                },
                '84': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC84",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Middle East",
                    "clinical_significance": "Emerging clone in Middle Eastern hospitals",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=84"
                },
                '94': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC94",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Europe",
                    "clinical_significance": "Clone associated with hospital outbreaks in Europe",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-58", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=94"
                },
                '98': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC98",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Asia",
                    "clinical_significance": "Clone prevalent in Asian hospital settings",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=98"
                },
                '106': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC106",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "South America",
                    "clinical_significance": "Clone associated with outbreaks in South America",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-72", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=106"
                },
                '111': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC111",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Africa",
                    "clinical_significance": "Emerging clone in African healthcare settings",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=111"
                },
                '136': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC136",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Europe, Mediterranean",
                    "clinical_significance": "Clone associated with OXA-58 production",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-58", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=136"
                },
                '158': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC158",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Middle East",
                    "clinical_significance": "Clone associated with NDM-1 production",
                    "common_resistance": ["Carbapenem-resistant", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["NDM-1", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=158"
                },
                '176': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC176",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Asia",
                    "clinical_significance": "Clone associated with OXA-23 and IMP co-production",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "IMP-type", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=176"
                },
                '195': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC195",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Widely distributed clone with diverse resistance patterns",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=195"
                },
                '208': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC208",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Common clone in hospital environments",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=208"
                },
                '229': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC229",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Europe",
                    "clinical_significance": "Clone associated with nosocomial infections in European hospitals",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-58", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=229"
                },
                '254': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC254",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Asia-Pacific",
                    "clinical_significance": "Emerging clone in Asia-Pacific region",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=254"
                },
                '281': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC281",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Middle East, North Africa",
                    "clinical_significance": "Clone associated with OXA-72 production",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-72", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=281"
                },
                '369': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC369",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "South America",
                    "clinical_significance": "Clone prevalent in Brazilian hospitals",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=369"
                },
                '450': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC450",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Widely distributed clone with high resistance potential",
                    "common_resistance": ["Carbapenem-resistant", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "NDM-1", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=450"
                },
                '499': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC499",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Europe",
                    "clinical_significance": "Clone associated with OXA-48-like enzymes in Europe",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-48-like", "OXA-51-like"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=499"
                },
                '513': {
                    "international_clone": "Regional Clone",
                    "clonal_complex": "CC513",
                    "classification": "Epidemic Clone",
                    "geographic_distribution": "Asia",
                    "clinical_significance": "Emerging clone in Asian ICUs",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "Colistin-resistant"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-51-like", "mcr-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=513"
                }
            }
        }
        
        # Get the appropriate database
        db = lineage_db.get(scheme, {})
        
        # Check if ST is in database
        if st in db:
            return db[st]
        else:
            # For unknown STs
            if st.isdigit():
                # Try to find IC mapping
                ic_status = "Unknown"
                for ic, mapping in self.ic_mapping.items():
                    if mapping.get(scheme) == st:
                        ic_status = ic
                        break
                
                # Generate PubMLST link
                scheme_id = "1" if scheme == "abaumannii" else "2"
                pubmlst_link = f"https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id={scheme_id}&ST={st}"
                
                return {
                    "international_clone": ic_status,
                    "clonal_complex": f"Unknown (ST{st})",
                    "classification": "Not in database - novel or uncommon",
                    "geographic_distribution": "Unknown",
                    "clinical_significance": f"ST{st} is not currently in the AcinetoScope {scheme} MLST database.",
                    "common_resistance": ["Unknown - requires further analysis"],
                    "outbreak_potential": "UNKNOWN",
                    "typical_resistance_genes": ["Unknown"],
                    "pubmlst_link": pubmlst_link
                }
            else:
                # For non-numeric STs (UNKNOWN, -, etc.)
                return {
                    "international_clone": "Not Assigned",
                    "clonal_complex": "Not Assigned",
                    "classification": "MLST typing failed",
                    "geographic_distribution": "N/A",
                    "clinical_significance": "Could not determine sequence type.",
                    "common_resistance": ["Cannot determine"],
                    "outbreak_potential": "UNKNOWN",
                    "typical_resistance_genes": ["Cannot determine"],
                    "pubmlst_link": "https://pubmlst.org/organisms/acinetobacter-baumannii"
                }    

    def get_empty_results(self, sample_name: str, scheme: str) -> Dict:
        """Return empty results structure"""
        return {
            "sample": sample_name,
            "st": "UNKNOWN",
            "scheme": scheme,
            "scheme_display": self.scheme_display_names.get(scheme, scheme.upper()),
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "mlst_assigned": False,
            "detected_scheme": "Unknown"
        }

    def get_fallback_results(self, sample_name: str, scheme: str) -> Dict:
        """Fallback when MLST fails"""
        return {
            "sample": sample_name,
            "st": "UNKNOWN",
            "scheme": scheme,
            "scheme_display": self.scheme_display_names.get(scheme, scheme.upper()),
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "mlst_assigned": False,
            "detected_scheme": "Unknown",
            "error": "MLST analysis failed"
        }

    def generate_output_files(self, mlst_results: Dict, output_dir: Path):
        """Generate output files: HTML, TXT, and TSV"""
        # 1. HTML Report
        self.generate_html_report(mlst_results, output_dir)
        
        # 2. Text Report
        self.generate_text_report(mlst_results, output_dir)
        
        # 3. TSV Report
        self.generate_tsv_report(mlst_results, output_dir)
        
        # 4. JSON Report
        self.generate_json_report(mlst_results, output_dir)

    def generate_text_report(self, mlst_results: Dict, output_dir: Path):
        """Generate detailed text report"""
        # NO TRUNCATION: Show full sample name
        sample_display = mlst_results['sample']
        
        report = f"""ACINETOSCOPE - MLST Analysis Report
======================================

Sample: {sample_display}
Original File: {mlst_results.get('original_filename', mlst_results['sample'])}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
MLST Scheme: {mlst_results['scheme_display']}
MLST Database: {mlst_results['scheme']}
Detected Scheme: {mlst_results.get('detected_scheme', 'Unknown')}

MLST TYPING RESULTS:
-------------------
Sequence Type (ST): {mlst_results['st']}
Confidence: {mlst_results['confidence']}
MLST Status: {'Assigned' if mlst_results['mlst_assigned'] else 'Not Assigned'}

Allele Profile ({mlst_results['scheme_display']} scheme):
{mlst_results['allele_profile']}

Detailed Alleles:
"""
        for gene, allele in mlst_results['alleles'].items():
            report += f"- {gene}: {allele}\n"
        
        if not mlst_results['alleles']:
            report += "- No alleles detected\n"
        
        # Lineage information
        report += f"""
LINEAGE ANALYSIS:
-----------------
International Clone: {mlst_results.get('international_clone', 'Unknown')}
Clonal Complex: {mlst_results.get('clonal_complex', 'Unknown')}
Classification: {mlst_results.get('classification', 'Unknown')}
Geographic Distribution: {mlst_results.get('geographic_distribution', 'Unknown')}
Clinical Significance: {mlst_results.get('clinical_significance', 'Unknown')}
Outbreak Potential: {mlst_results.get('outbreak_potential', 'UNKNOWN')}

Common Resistance: {', '.join(mlst_results.get('common_resistance', ['Unknown']))}
Typical Resistance Genes: {', '.join(mlst_results.get('typical_resistance_genes', ['Unknown']))}

PubMLST Link: {mlst_results.get('pubmlst_link', 'https://pubmlst.org/organisms/acinetobacter-baumannii')}
"""
        
        if 'error' in mlst_results:
            report += f"\nERROR: {mlst_results['error']}\n"
        
        with open(output_dir / "mlst_report.txt", 'w') as f:
            f.write(report)

    def generate_tsv_report(self, mlst_results: Dict, output_dir: Path):
        """Generate simple TSV report"""
        sample_display = mlst_results['sample']
        
        tsv_content = f"""Sample\tOriginal_File\tScheme\tMLST_Database\tST\tMLST_Status\tInternational_Clone\tClonal_Complex\tClassification\tOutbreak_Potential\tResistance_Profile\tBiofilm_Formation\tSurvival_on_Surfaces\tConfidence\tAllele_Profile
{sample_display}\t{mlst_results.get('original_filename', mlst_results['sample'])}\t{mlst_results['scheme_display'].lower()}\t{mlst_results['scheme']}\t{mlst_results['st']}\t{'Assigned' if mlst_results['mlst_assigned'] else 'Not Assigned'}\t{mlst_results.get('international_clone', 'Unknown')}\t{mlst_results.get('clonal_complex', 'Unknown')}\t{mlst_results.get('classification', 'Unknown')}\t{mlst_results.get('outbreak_potential', 'UNKNOWN')}\t{', '.join(mlst_results.get('common_resistance', ['Unknown']))}\tUnknown\tUnknown\t{mlst_results['confidence']}\t{mlst_results['allele_profile']}
"""
        
        with open(output_dir / "mlst_report.tsv", 'w') as f:
            f.write(tsv_content)

    def generate_json_report(self, mlst_results: Dict, output_dir: Path):
        """Generate JSON report"""
        json_report = {
            "metadata": {
                "sample": mlst_results['sample'],
                "original_filename": mlst_results.get('original_filename', mlst_results['sample']),
                "analysis_date": datetime.now().isoformat(),
                "scheme": mlst_results['scheme'],
                "scheme_display": mlst_results['scheme_display'],
                "version": "1.0",
                "tool": "AcinetoScope"
            },
            "mlst_results": {
                "sequence_type": mlst_results['st'],
                "confidence": mlst_results['confidence'],
                "mlst_assigned": mlst_results['mlst_assigned'],
                "allele_profile": mlst_results['allele_profile'],
                "alleles": mlst_results['alleles'],
                "detected_scheme": mlst_results.get('detected_scheme', 'Unknown')
            },
            "lineage_information": {
                "international_clone": mlst_results.get('international_clone', 'Unknown'),
                "clonal_complex": mlst_results.get('clonal_complex', 'Unknown'),
                "classification": mlst_results.get('classification', 'Unknown'),
                "geographic_distribution": mlst_results.get('geographic_distribution', 'Unknown'),
                "clinical_significance": mlst_results.get('clinical_significance', 'Unknown'),
                "outbreak_potential": mlst_results.get('outbreak_potential', 'UNKNOWN'),
                "common_resistance": mlst_results.get('common_resistance', ['Unknown']),
                "typical_resistance_genes": mlst_results.get('typical_resistance_genes', ['Unknown']),
                "pubmlst_link": mlst_results.get('pubmlst_link', 'https://pubmlst.org/organisms/acinetobacter-baumannii')
            }
        }
        
        if 'error' in mlst_results:
            json_report["error"] = mlst_results['error']
        
        with open(output_dir / "mlst_report.json", 'w') as f:
            json.dump(json_report, f, indent=2)

    def generate_html_report(self, mlst_results: Dict, output_dir: Path):
        """Generate beautiful HTML report with NO TRUNCATION of sample names or allele profiles"""
        random_quote = self.get_random_quote()
        
        # NO TRUNCATION: Show full names
        sample = mlst_results['sample']
        sample_display = sample  # No truncation
        
        original_filename = mlst_results.get('original_filename', sample)
        original_display = original_filename  # No truncation
        
        st = mlst_results['st']
        scheme = mlst_results['scheme']
        scheme_display = mlst_results['scheme_display']
        confidence = mlst_results['confidence']
        international_clone = mlst_results.get('international_clone', 'Unknown')
        clonal_complex = mlst_results.get('clonal_complex', 'Unknown')
        allele_profile = mlst_results['allele_profile'] or "No allele profile detected"
        classification = mlst_results.get('classification', 'Unknown')
        geographic_distribution = mlst_results.get('geographic_distribution', 'Unknown')
        outbreak_potential = mlst_results.get('outbreak_potential', 'UNKNOWN')
        clinical_significance = mlst_results.get('clinical_significance', 'Unknown')
        common_resistance = mlst_results.get('common_resistance', ['Unknown'])
        typical_resistance_genes = mlst_results.get('typical_resistance_genes', ['Unknown'])
        pubmlst_link = mlst_results.get('pubmlst_link', 'https://pubmlst.org/organisms/acinetobacter-baumannii')
        
        # Build alleles HTML - ensure gene names are fully visible
        alleles_html = ''
        for gene, allele in mlst_results['alleles'].items():
            alleles_html += f'''                <div class="allele-card">
                    <div style="font-size: 14px; opacity: 0.9; word-wrap: break-word;">{gene}</div>
                    <div style="font-size: 18px; word-wrap: break-word;">{allele}</div>
                </div>
'''
        
        if not mlst_results['alleles']:
            alleles_html = '''                <div class="allele-card" style="grid-column: 1 / -1; background: linear-gradient(135deg, #fca5a5 0%, #ef4444 100%);">
                    <div style="font-size: 14px; word-wrap: break-word;">No alleles detected</div>
                    <div style="font-size: 12px; word-wrap: break-word;">MLST typing may have failed</div>
                </div>
'''
        
        # Build resistance genes HTML - no truncation
        resistance_html = ''
        for gene in typical_resistance_genes:
            resistance_html += f'''                        <div class="resistance-card" style="word-wrap: break-word;">{gene}</div>
'''
        
        # Determine confidence class
        confidence_class = confidence.lower()
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ACINETOSCOPE - MLST Analysis Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 1800px; /* Increased max width for better display */
            margin: 0 auto;
        }}
        
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }}
        
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
        }}
        
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
        }}
        
        .quote-text {{
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
            color: #ffffff;
            word-wrap: break-word;
        }}
        
        .quote-author {{
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
            word-wrap: break-word;
        }}
        
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            overflow-wrap: break-word;
        }}
        
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
            word-wrap: break-word;
        }}
        
        .report-section h3 {{
            color: #1e40af;
            margin-top: 20px;
            margin-bottom: 10px;
            font-size: 18px;
            word-wrap: break-word;
        }}
        
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); /* Increased min width */
            gap: 20px;
            margin-top: 15px;
        }}
        
        .metric-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
            overflow-wrap: break-word;
        }}
        
        .metric-label {{
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 5px;
            word-wrap: break-word;
        }}
        
        .metric-value {{
            font-size: 24px;
            font-weight: bold;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }}
        
        .allele-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); /* Increased min width */
            gap: 15px;
            margin-top: 15px;
        }}
        
        .allele-card {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
            text-align: center;
            font-weight: bold;
            overflow-wrap: break-word;
            word-wrap: break-word;
        }}
        
        .confidence-high {{
            color: #16a34a;
            font-weight: bold;
            word-wrap: break-word;
        }}
        
        .confidence-medium {{
            color: #f59e0b;
            font-weight: bold;
            word-wrap: break-word;
        }}
        
        .confidence-low {{
            color: #dc2626;
            font-weight: bold;
            word-wrap: break-word;
        }}
        
        .resistance-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(250px, 1fr)); /* Increased min width */
            gap: 10px;
            margin-top: 15px;
        }}
        
        .resistance-card {{
            background: #fee2e2;
            color: #991b1b;
            padding: 12px;
            border-radius: 6px;
            text-align: center;
            font-weight: bold;
            border-left: 4px solid #ef4444;
            overflow-wrap: break-word;
            word-wrap: break-word;
        }}
        
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
            overflow-wrap: break-word;
        }}
        
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
            word-wrap: break-word;
        }}
        
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
            overflow-wrap: break-word;
        }}
        
        .profile-box {{
            background: #f8fafc;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border-left: 4px solid #3b82f6;
            overflow-x: auto; /* Allow horizontal scrolling for long profiles */
            white-space: pre-wrap; /* Preserve formatting but allow wrapping */
            word-wrap: break-word;
            overflow-wrap: break-word;
        }}
        
        .clinical-box {{
            background: #f0f9ff;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border-left: 4px solid #0ea5e9;
            overflow-wrap: break-word;
            word-wrap: break-word;
        }}
        
        .filename-info {{
            background: #f3f4f6;
            padding: 10px;
            border-radius: 6px;
            margin: 10px 0;
            font-size: 12px;
            color: #6b7280;
            border-left: 3px solid #9ca3af;
            overflow-wrap: break-word;
            word-wrap: break-word;
        }}
        
        .full-width {{
            grid-column: 1 / -1;
        }}
        
        .no-truncate {{
            white-space: normal !important;
            word-wrap: break-word !important;
            overflow-wrap: break-word !important;
        }}
        
        @media (max-width: 768px) {{
            .ascii-art {{
                font-size: 6px;
            }}
            .allele-grid {{
                grid-template-columns: 1fr;
            }}
            .metric-card {{
                padding: 15px;
            }}
            .metrics-grid {{
                grid-template-columns: 1fr;
            }}
            .container {{
                max-width: 100%;
                overflow-x: auto;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art"> █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝</div>
            </div>
            
            <div class="quote-container">
                <div class="quote-text">"{random_quote['text']}"</div>
                <div class="quote-author">— {random_quote['author']}</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 Sample Information</h2>
            <div class="filename-info">
                <strong>Original File:</strong> {original_display}<br>
                <strong>Display Name:</strong> {sample_display}
            </div>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sample Name</div>
                    <div class="metric-value no-truncate">{sample_display}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">MLST Scheme</div>
                    <div class="metric-value">{scheme_display}</div>
                </div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>🧬 MLST Typing Results</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sequence Type</div>
                    <div class="metric-value">ST{st}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">International Clone</div>
                    <div class="metric-value no-truncate">{international_clone}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Confidence</div>
                    <div class="metric-value confidence-{confidence_class}">{confidence}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Clonal Complex</div>
                    <div class="metric-value no-truncate">{clonal_complex}</div>
                </div>
            </div>
            
            <h3>Allele Profile</h3>
            <div class="profile-box">
                <code style="font-size: 16px; color: #1e40af; font-weight: bold; white-space: pre-wrap; word-wrap: break-word;">{allele_profile}</code>
            </div>
            
            <h3>Individual Alleles</h3>
            <div class="allele-grid">
{alleles_html}            </div>
        </div>
        
        <div class="report-section">
            <h2>🌍 Lineage Information</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Classification</div>
                    <div class="metric-value no-truncate">{classification}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Geographic Distribution</div>
                    <div class="metric-value no-truncate">{geographic_distribution}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Outbreak Potential</div>
                    <div class="metric-value no-truncate">{outbreak_potential}</div>
                </div>
            </div>
            
            <h3>Clinical Significance</h3>
            <div class="clinical-box">
                <p style="margin: 0; line-height: 1.6; color: #374151; word-wrap: break-word;">{clinical_significance}</p>
            </div>
            
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin-top: 20px;">
                <div>
                    <h3>Resistance Profile</h3>
                    <p><strong>Common Resistance:</strong> {', '.join(common_resistance)}</p>
                </div>
                <div>
                    <h3>Typical Resistance Genes</h3>
                    <div class="resistance-grid">
{resistance_html}                    </div>
                </div>
            </div>
            
            <div style="margin-top: 20px; text-align: center;">
                <a href="{pubmlst_link}" target="_blank" style="background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 10px 20px; border-radius: 6px; text-decoration: none; font-weight: bold; display: inline-block;">
                    🔗 View on PubMLST
                </a>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>ACINETOSCOPE</strong> - A. baumannii MLST Analysis Module</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <div class="authorship">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>
</body>
</html>'''
        
        with open(output_dir / "mlst_report.html", 'w', encoding='utf-8') as f:
            f.write(html_content)

    def run_mlst_dual_scheme(self, input_file: Path, output_dir: Path) -> Dict[str, Dict]:
        """Run MLST analysis with both schemes"""
        print(f"\n{'='*80}")
        print(f"🧬 ACINETOSCOPE - Dual-Scheme MLST Analysis")
        print(f"{'='*80}")
        print(f"Sample: {input_file.name}")
        print(f"{'='*80}")
        
        # Create main output directory
        main_output_dir = output_dir / "mlst_results"
        main_output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n📁 Output directory: {main_output_dir}")
        
        # Create scheme-specific directories
        oxford_dir = main_output_dir / "Oxford"
        pasteur_dir = main_output_dir / "Pasteur"
        oxford_dir.mkdir(parents=True, exist_ok=True)
        pasteur_dir.mkdir(parents=True, exist_ok=True)
        
        print("\n🔍 Running Oxford MLST Scheme (abaumannii)...")
        oxford_results = self.run_mlst_single(input_file, oxford_dir, "abaumannii")
        
        print("\n🔍 Running Pasteur MLST Scheme (abaumannii_2)...")
        pasteur_results = self.run_mlst_single(input_file, pasteur_dir, "abaumannii_2")
        
        # Create combined summary
        self.create_combined_summary(oxford_results, pasteur_results, main_output_dir, input_file.name)
        
        return {
            "abaumannii": oxford_results,
            "abaumannii_2": pasteur_results
        }

    def create_combined_summary(self, oxford_results: Dict, pasteur_results: Dict, 
                               output_dir: Path, sample_name: str):
        """Create combined summary report"""
        print("\n📊 Creating Combined Summary...")
        
        combined_dir = output_dir / "Combined"
        combined_dir.mkdir(parents=True, exist_ok=True)
        
        # Extract key information
        oxford_st = oxford_results.get('st', 'UNKNOWN')
        pasteur_st = pasteur_results.get('st', 'UNKNOWN')
        oxford_ic = oxford_results.get('international_clone', 'Unknown')
        pasteur_ic = pasteur_results.get('international_clone', 'Unknown')
        
        # Check for International Clone consistency
        ic_match = "✅ MATCH" if oxford_ic == pasteur_ic else "⚠️ MISMATCH"
        
        # Create combined report
        combined_report = f"""ACINETOSCOPE - Dual-Scheme MLST Summary
==========================================

Sample: {sample_name}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

SCHEME COMPARISON:
==================

Oxford Scheme (abaumannii):
--------------------------------------------------------
Sequence Type: ST{oxford_st}
International Clone: {oxford_ic}
Confidence: {oxford_results.get('confidence', 'LOW')}
Allele Profile: {oxford_results.get('allele_profile', 'No profile')}

Pasteur Scheme (abaumannii_2):
----------------------------------------------------------
Sequence Type: ST{pasteur_st}
International Clone: {pasteur_ic}
Confidence: {pasteur_results.get('confidence', 'LOW')}
Allele Profile: {pasteur_results.get('allele_profile', 'No profile')}

COMPARISON RESULTS:
------------------
International Clone Match: {ic_match}

RECOMMENDATIONS:
---------------
"""
        
        # Add recommendations
        if oxford_st == 'UNKNOWN' and pasteur_st == 'UNKNOWN':
            combined_report += "1. ⚠️ MLST typing failed for both schemes.\n"
            combined_report += "2. Check input FASTA file quality and completeness.\n"
        elif oxford_st != 'UNKNOWN' and pasteur_st != 'UNKNOWN':
            if oxford_ic == pasteur_ic:
                combined_report += f"1. ✅ Consistent typing: Strain belongs to {oxford_ic}.\n"
                combined_report += f"2. Use ST{pasteur_st} (Pasteur) as primary identifier.\n"
                combined_report += f"3. Corresponds to ST{oxford_st} in Oxford scheme.\n"
            else:
                combined_report += f"1. ⚠️ Scheme discrepancy detected.\n"
                combined_report += f"2. Verify with additional typing methods.\n"
                combined_report += f"3. Oxford: ST{oxford_st} ({oxford_ic})\n"
                combined_report += f"4. Pasteur: ST{pasteur_st} ({pasteur_ic})\n"
        elif oxford_st != 'UNKNOWN':
            combined_report += f"1. Use Oxford ST{oxford_st} as primary result.\n"
            combined_report += f"2. Pasteur scheme failed or returned unknown.\n"
            combined_report += f"3. International Clone: {oxford_ic}\n"
        else:
            combined_report += f"1. Use Pasteur ST{pasteur_st} as primary result.\n"
            combined_report += f"2. Oxford scheme failed or returned unknown.\n"
            combined_report += f"3. International Clone: {pasteur_ic}\n"
        
        combined_report += f"""
FILES GENERATED:
----------------
• Oxford Scheme: {output_dir}/Oxford/{sample_name.split('.')[0]}/
• Pasteur Scheme: {output_dir}/Pasteur/{sample_name.split('.')[0]}/
• Combined Summary: {combined_dir}/

For detailed analysis, refer to individual scheme reports.
"""
        
        # Write combined report
        with open(combined_dir / "combined_summary.txt", 'w') as f:
            f.write(combined_report)
        
        # Also create JSON summary
        json_summary = {
            "sample": sample_name,
            "analysis_date": datetime.now().isoformat(),
            "oxford": {
                "st": oxford_st,
                "international_clone": oxford_ic,
                "confidence": oxford_results.get('confidence', 'LOW'),
                "allele_profile": oxford_results.get('allele_profile', '')
            },
            "pasteur": {
                "st": pasteur_st,
                "international_clone": pasteur_ic,
                "confidence": pasteur_results.get('confidence', 'LOW'),
                "allele_profile": pasteur_results.get('allele_profile', '')
            },
            "comparison": {
                "ic_match": oxford_ic == pasteur_ic,
                "recommendation": "See combined_summary.txt for details"
            }
        }
        
        with open(combined_dir / "combined_summary.json", 'w') as f:
            json.dump(json_summary, f, indent=2)
        
        print(f"✅ Combined summary created: {combined_dir}/")

    def run_mlst_dual_batch(self, input_path: str, output_dir: Path) -> Dict[str, Dict]:
        """Run dual-scheme MLST analysis for multiple files with glob pattern support"""
        print("🔍 Searching for FASTA files...")
        fasta_files = self.find_fasta_files(input_path)
        
        if not fasta_files:
            print("❌ No FASTA files found!")
            return {}
        
        print(f"📁 Found {len(fasta_files)} FASTA files")
        for i, fasta_file in enumerate(fasta_files, 1):
            print(f"  {i}. {fasta_file}")
        
        # Create main output directory
        main_output_dir = output_dir / "Dual_Scheme_MLST"
        main_output_dir.mkdir(parents=True, exist_ok=True)
        
        results = {}
        for fasta_file in fasta_files:
            # Create individual sample directory
            sample_base = fasta_file.stem
            sample_output_dir = main_output_dir / sample_base
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            
            # Run dual-scheme MLST
            print(f"\n{'='*80}")
            print(f"🧬 Processing: {fasta_file.name} with dual schemes")
            print(f"{'='*80}")
            
            # Create scheme-specific directories
            oxford_dir = sample_output_dir / "Oxford"
            pasteur_dir = sample_output_dir / "Pasteur"
            oxford_dir.mkdir(parents=True, exist_ok=True)
            pasteur_dir.mkdir(parents=True, exist_ok=True)
            
            # Run both schemes
            print("\n🔍 Running Oxford MLST Scheme (abaumannii)...")
            oxford_results = self.run_mlst_single(fasta_file, oxford_dir, "abaumannii")
            
            print("\n🔍 Running Pasteur MLST Scheme (abaumannii_2)...")
            pasteur_results = self.run_mlst_single(fasta_file, pasteur_dir, "abaumannii_2")
            
            # Create combined summary for this sample
            combined_dir = sample_output_dir / "Combined"
            combined_dir.mkdir(parents=True, exist_ok=True)
            self.create_combined_summary(oxford_results, pasteur_results, combined_dir, fasta_file.name)
            
            results[fasta_file.name] = {
                "abaumannii": oxford_results,
                "abaumannii_2": pasteur_results
            }
            
            st_display_oxford = oxford_results.get('st', 'UNKNOWN')
            st_display_pasteur = pasteur_results.get('st', 'UNKNOWN')
            print(f"✅ Completed: {fasta_file.name} -> Oxford:ST{st_display_oxford}, Pasteur:ST{st_display_pasteur}")
        
        # Create batch summary
        self.create_dual_batch_summary(results, main_output_dir)
        
        return results

    def create_dual_batch_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create comprehensive summary for dual-scheme batch analysis"""
        print(f"\n📊 Creating Dual-Scheme Batch Summary...")
        
        summary_data = []
        detailed_results = []
        
        for sample_name, schemes in all_results.items():
            oxford_results = schemes.get("abaumannii", {})
            pasteur_results = schemes.get("abaumannii_2", {})
            
            # For TSV/Excel
            summary_data.append({
                "Sample": sample_name,
                "Original_Filename": oxford_results.get('original_filename', sample_name),
                "Oxford_ST": oxford_results.get('st', 'UNKNOWN'),
                "Pasteur_ST": pasteur_results.get('st', 'UNKNOWN'),
                "Oxford_IC": oxford_results.get('international_clone', 'Unknown'),
                "Pasteur_IC": pasteur_results.get('international_clone', 'Unknown'),
                "Oxford_Confidence": oxford_results.get('confidence', 'LOW'),
                "Pasteur_Confidence": pasteur_results.get('confidence', 'LOW'),
                "Oxford_Allele_Profile": oxford_results.get('allele_profile', ''),
                "Pasteur_Allele_Profile": pasteur_results.get('allele_profile', ''),
                "IC_Match": "MATCH" if oxford_results.get('international_clone') == pasteur_results.get('international_clone') else "MISMATCH"
            })
            
            # For detailed JSON
            detailed_results.append({
                "sample": sample_name,
                "original_filename": oxford_results.get('original_filename', sample_name),
                "oxford_results": {
                    "sequence_type": oxford_results.get('st', 'UNKNOWN'),
                    "confidence": oxford_results.get('confidence', 'LOW'),
                    "mlst_assigned": oxford_results.get('mlst_assigned', False),
                    "allele_profile": oxford_results.get('allele_profile', ''),
                    "alleles": oxford_results.get('alleles', {})
                },
                "pasteur_results": {
                    "sequence_type": pasteur_results.get('st', 'UNKNOWN'),
                    "confidence": pasteur_results.get('confidence', 'LOW'),
                    "mlst_assigned": pasteur_results.get('mlst_assigned', False),
                    "allele_profile": pasteur_results.get('allele_profile', ''),
                    "alleles": pasteur_results.get('alleles', {})
                }
            })
        
        # Create DataFrame
        df = pd.DataFrame(summary_data)
        
        # Save to CSV (always works)
        csv_file = output_dir / "dual_scheme_mlst_summary.csv"
        df.to_csv(csv_file, index=False)
        print(f"  ✅ CSV summary: {csv_file}")
        
        # Save to Excel if available
        if self.has_excel_support:
            try:
                excel_file = output_dir / "dual_scheme_mlst_summary.xlsx"
                df.to_excel(excel_file, index=False)
                print(f"  ✅ Excel summary: {excel_file}")
            except Exception as e:
                print(f"  ⚠️ Excel export failed: {e}")
        else:
            print(f"  ⚠️ Excel export not available (install openpyxl: pip install openpyxl)")
        
        # Save detailed JSON
        json_summary = {
            "metadata": {
                "analysis_date": datetime.now().isoformat(),
                "schemes": ["Oxford (abaumannii)", "Pasteur (abaumannii_2)"],
                "samples_analyzed": len(all_results),
                "version": "1.0",
                "tool": "AcinetoScope"
            },
            "results": detailed_results,
            "summary_statistics": {
                "total_samples": len(all_results),
                "oxford_assigned": sum(1 for r in all_results.values() if r.get("abaumannii", {}).get('mlst_assigned', False)),
                "pasteur_assigned": sum(1 for r in all_results.values() if r.get("abaumannii_2", {}).get('mlst_assigned', False)),
                "ic_matches": sum(1 for r in summary_data if r['IC_Match'] == 'MATCH'),
                "ic_mismatches": sum(1 for r in summary_data if r['IC_Match'] == 'MISMATCH')
            }
        }
        
        json_file = output_dir / "dual_scheme_mlst_summary.json"
        with open(json_file, 'w') as f:
            json.dump(json_summary, f, indent=2)
        print(f"  ✅ JSON summary: {json_file}")
        
        # Create HTML summary
        self.create_dual_batch_html_summary(summary_data, output_dir)
        
        print(f"✅ Dual-scheme batch summary created in: {output_dir}/")

    def create_dual_batch_html_summary(self, summary_data: List[Dict], output_dir: Path):
        """Create HTML batch summary report for dual-scheme analysis"""
        total_samples = len(summary_data)
        ic_matches = sum(1 for row in summary_data if row['IC_Match'] == 'MATCH')
        ic_mismatches = total_samples - ic_matches
        
        # Build table rows with NO TRUNCATION
        table_rows = ""
        for i, row in enumerate(summary_data, 1):
            sample_name = row['Sample']
            sample_display = sample_name  # NO TRUNCATION
            
            ic_match_class = "match" if row['IC_Match'] == 'MATCH' else "mismatch"
            
            table_rows += f"""
            <tr>
                <td>{i}</td>
                <td><span class="sample-name no-truncate" title="{sample_name}">{sample_display}</span></td>
                <td>ST{row['Oxford_ST']}</td>
                <td>ST{row['Pasteur_ST']}</td>
                <td class="no-truncate">{row['Oxford_IC']}</td>
                <td class="no-truncate">{row['Pasteur_IC']}</td>
                <td class="{ic_match_class} no-truncate">{row['IC_Match']}</td>
                <td><code class="no-truncate" style="white-space: pre-wrap; word-wrap: break-word;">{row['Oxford_Allele_Profile']}</code></td>
                <td><code class="no-truncate" style="white-space: pre-wrap; word-wrap: break-word;">{row['Pasteur_Allele_Profile']}</code></td>
            </tr>
            """
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ACINETOSCOPE - Dual-Scheme MLST Batch Summary</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 2000px; /* Very wide for dual-scheme data */
            margin: 0 auto;
        }}
        
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 15px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }}
        
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 8px;
            line-height: 1.1;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 25px;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        
        .stat-value {{
            font-size: 28px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .stat-label {{
            font-size: 12px;
            opacity: 0.9;
        }}
        
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            overflow-x: auto; /* Allow horizontal scrolling */
        }}
        
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 15px;
            font-size: 22px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
            table-layout: auto; /* Allow columns to expand */
        }}
        
        th {{
            background: #3b82f6;
            color: white;
            padding: 12px 15px;
            text-align: left;
            position: sticky;
            top: 0;
            white-space: nowrap;
        }}
        
        td {{
            padding: 10px 15px;
            border-bottom: 1px solid #e5e7eb;
            vertical-align: top;
        }}
        
        tr:hover {{
            background-color: #f3f4f6;
        }}
        
        .match {{
            color: #16a34a;
            font-weight: bold;
        }}
        
        .mismatch {{
            color: #dc2626;
            font-weight: bold;
        }}
        
        .sample-name {{
            cursor: help;
        }}
        
        .no-truncate {{
            white-space: normal !important;
            word-wrap: break-word !important;
            overflow-wrap: break-word !important;
        }}
        
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
        }}
        
        .download-links {{
            display: flex;
            justify-content: center;
            gap: 15px;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        
        .download-btn {{
            background: linear-gradient(135deg, #10b981 0%, #059669 100%);
            color: white;
            padding: 10px 20px;
            border-radius: 6px;
            text-decoration: none;
            font-weight: bold;
            display: inline-flex;
            align-items: center;
            gap: 8px;
        }}
        
        .download-btn:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
        }}
        
        @media (max-width: 768px) {{
            .ascii-art {{
                font-size: 5px;
            }}
            table {{
                display: block;
                overflow-x: auto;
            }}
            th, td {{
                padding: 8px 10px;
                font-size: 14px;
            }}
            .download-links {{
                flex-direction: column;
                align-items: center;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art"> █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝</div>
            </div>
            
            <h1 style="color: white; margin-bottom: 10px;">Dual-Scheme MLST Batch Analysis</h1>
            <p style="color: #d1d5db; margin-bottom: 20px;">Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{total_samples}</div>
                <div class="stat-label">SAMPLES ANALYZED</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{ic_matches}</div>
                <div class="stat-label">IC MATCHES</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{ic_mismatches}</div>
                <div class="stat-label">IC MISMATCHES</div>
            </div>
        </div>
        
        <div class="download-links">
            <a href="dual_scheme_mlst_summary.csv" class="download-btn">
                📄 Download CSV Summary
            </a>
            <a href="dual_scheme_mlst_summary.json" class="download-btn">
                🔧 Download JSON Summary
            </a>
'''
        
        # Add Excel download link only if available
        if self.has_excel_support:
            html_content += f'''            <a href="dual_scheme_mlst_summary.xlsx" class="download-btn">
                📊 Download Excel Summary
            </a>
'''
        
        html_content += f'''        </div>
        
        <div class="report-section">
            <h2>📋 Dual-Scheme Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>#</th>
                        <th>Sample</th>
                        <th>Oxford ST</th>
                        <th>Pasteur ST</th>
                        <th>Oxford IC</th>
                        <th>Pasteur IC</th>
                        <th>IC Match</th>
                        <th>Oxford Allele Profile</th>
                        <th>Pasteur Allele Profile</th>
                    </tr>
                </thead>
                <tbody>
                    {table_rows}
                </tbody>
            </table>
        </div>
        
        <div class="footer">
            <p><strong>ACINETOSCOPE</strong> - A. baumannii MLST Analysis Module</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <p>For individual sample reports, navigate to respective sample directories.</p>
            <div style="margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1); border-radius: 8px; font-size: 12px;">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>
</body>
</html>'''
        
        html_file = output_dir / "dual_scheme_mlst_summary.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"  ✅ HTML summary: {html_file}")

    def run_mlst_batch(self, input_path: str, output_dir: Path, scheme: str = "abaumannii") -> Dict[str, Dict]:
        """Run MLST analysis for multiple files with glob pattern support"""
        print("🔍 Searching for FASTA files...")
        fasta_files = self.find_fasta_files(input_path)
        
        if not fasta_files:
            print("❌ No FASTA files found!")
            return {}
        
        print(f"📁 Found {len(fasta_files)} FASTA files")
        for i, fasta_file in enumerate(fasta_files, 1):
            print(f"  {i}. {fasta_file}")
        
        # Create scheme-specific output directory
        scheme_display = self.scheme_display_names.get(scheme, scheme.upper())
        scheme_output_dir = output_dir / f"{scheme_display}_MLST"
        scheme_output_dir.mkdir(parents=True, exist_ok=True)
        
        results = {}
        for fasta_file in fasta_files:
            # Create individual sample directory
            sample_base = fasta_file.stem
            sample_output_dir = scheme_output_dir / sample_base
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            
            # Run MLST
            result = self.run_mlst_single(fasta_file, sample_output_dir, scheme)
            results[fasta_file.name] = result
        
        # Create scheme-specific summaries
        self.create_scheme_summary(results, scheme_output_dir, scheme)
        
        return results

    def create_scheme_summary(self, all_results: Dict[str, Dict], output_dir: Path, scheme: str):
        """Create comprehensive summary for a specific scheme"""
        print(f"\n📊 Creating {scheme} Summary...")
        
        summary_data = []
        detailed_results = []
        
        for sample_name, result in all_results.items():
            # For TSV/Excel
            summary_data.append({
                "Sample": sample_name,
                "Original_Filename": result.get('original_filename', sample_name),
                "ST": result.get('st', 'UNKNOWN'),
                "Scheme": result.get('scheme_display', 'Unknown'),
                "International_Clone": result.get('international_clone', 'Unknown'),
                "Clonal_Complex": result.get('clonal_complex', 'Unknown'),
                "Classification": result.get('classification', 'Unknown'),
                "Confidence": result.get('confidence', 'LOW'),
                "Allele_Profile": result.get('allele_profile', ''),
                "MLST_Assigned": result.get('mlst_assigned', False)
            })
            
            # For detailed JSON
            detailed_results.append({
                "sample": sample_name,
                "original_filename": result.get('original_filename', sample_name),
                "mlst_results": {
                    "sequence_type": result.get('st', 'UNKNOWN'),
                    "confidence": result.get('confidence', 'LOW'),
                    "mlst_assigned": result.get('mlst_assigned', False),
                    "allele_profile": result.get('allele_profile', ''),
                    "alleles": result.get('alleles', {})
                },
                "lineage_information": {
                    "international_clone": result.get('international_clone', 'Unknown'),
                    "clonal_complex": result.get('clonal_complex', 'Unknown'),
                    "classification": result.get('classification', 'Unknown'),
                    "geographic_distribution": result.get('geographic_distribution', 'Unknown'),
                    "outbreak_potential": result.get('outbreak_potential', 'UNKNOWN')
                }
            })
        
        # Create DataFrame
        df = pd.DataFrame(summary_data)
        
        # Save to CSV (always works)
        csv_file = output_dir / f"{self.scheme_display_names.get(scheme, scheme).lower()}_mlst_summary.csv"
        df.to_csv(csv_file, index=False)
        print(f"  ✅ CSV summary: {csv_file}")
        
        # Save to Excel if available
        if self.has_excel_support:
            try:
                excel_file = output_dir / f"{self.scheme_display_names.get(scheme, scheme).lower()}_mlst_summary.xlsx"
                df.to_excel(excel_file, index=False)
                print(f"  ✅ Excel summary: {excel_file}")
            except Exception as e:
                print(f"  ⚠️ Excel export failed: {e}")
        else:
            print(f"  ⚠️ Excel export not available (install openpyxl: pip install openpyxl)")
        
        # Save detailed JSON
        json_summary = {
            "metadata": {
                "analysis_date": datetime.now().isoformat(),
                "scheme": scheme,
                "scheme_display": self.scheme_display_names.get(scheme, scheme.upper()),
                "samples_analyzed": len(all_results),
                "version": "1.0",
                "tool": "AcinetoScope"
            },
            "results": detailed_results,
            "summary_statistics": {
                "total_samples": len(all_results),
                "assigned_st": sum(1 for r in all_results.values() if r.get('mlst_assigned', False)),
                "unknown_st": sum(1 for r in all_results.values() if not r.get('mlst_assigned', False)),
                "high_confidence": sum(1 for r in all_results.values() if r.get('confidence') == 'HIGH'),
                "low_confidence": sum(1 for r in all_results.values() if r.get('confidence') == 'LOW')
            }
        }
        
        json_file = output_dir / f"{self.scheme_display_names.get(scheme, scheme).lower()}_mlst_summary.json"
        with open(json_file, 'w') as f:
            json.dump(json_summary, f, indent=2)
        print(f"  ✅ JSON summary: {json_file}")
        
        # Create HTML summary
        self.create_html_batch_summary(summary_data, output_dir, scheme)
        
        print(f"✅ {self.scheme_display_names.get(scheme, scheme)} summary created in: {output_dir}/")

    def create_html_batch_summary(self, summary_data: List[Dict], output_dir: Path, scheme: str):
        """Create HTML batch summary report"""
        scheme_display = self.scheme_display_names.get(scheme, scheme.upper())
        
        # Statistics
        total_samples = len(summary_data)
        assigned_st = sum(1 for r in summary_data if r['ST'] != 'UNKNOWN')
        unknown_st = total_samples - assigned_st
        high_confidence = sum(1 for r in summary_data if r['Confidence'] == 'HIGH')
        
        # Build table rows with NO TRUNCATION
        table_rows = ""
        for i, row in enumerate(summary_data, 1):
            sample_name = row['Sample']
            sample_display = sample_name  # NO TRUNCATION
            
            st_class = "st-assigned" if row['ST'] != 'UNKNOWN' else "st-unknown"
            confidence_class = row['Confidence'].lower()
            
            table_rows += f"""
            <tr>
                <td>{i}</td>
                <td><span class="sample-name no-truncate" title="{sample_name}">{sample_display}</span></td>
                <td class="{st_class}">ST{row['ST']}</td>
                <td class="no-truncate">{row['International_Clone']}</td>
                <td class="confidence-{confidence_class}">{row['Confidence']}</td>
                <td><code class="no-truncate" style="white-space: pre-wrap; word-wrap: break-word;">{row['Allele_Profile']}</code></td>
            </tr>
            """
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ACINETOSCOPE - {scheme_display} MLST Batch Summary</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 2000px;
            margin: 0 auto;
        }}
        
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 15px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }}
        
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 8px;
            line-height: 1.1;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 25px;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        
        .stat-value {{
            font-size: 28px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .stat-label {{
            font-size: 12px;
            opacity: 0.9;
        }}
        
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            overflow-x: auto;
        }}
        
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 15px;
            font-size: 22px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
            table-layout: auto;
        }}
        
        th {{
            background: #3b82f6;
            color: white;
            padding: 12px 15px;
            text-align: left;
            position: sticky;
            top: 0;
            white-space: nowrap;
        }}
        
        td {{
            padding: 10px 15px;
            border-bottom: 1px solid #e5e7eb;
            vertical-align: top;
        }}
        
        tr:hover {{
            background-color: #f3f4f6;
        }}
        
        .st-assigned {{
            color: #16a34a;
            font-weight: bold;
        }}
        
        .st-unknown {{
            color: #dc2626;
            font-weight: bold;
        }}
        
        .confidence-high {{
            color: #16a34a;
            font-weight: bold;
        }}
        
        .confidence-low {{
            color: #dc2626;
            font-weight: bold;
        }}
        
        .sample-name {{
            cursor: help;
        }}
        
        .no-truncate {{
            white-space: normal !important;
            word-wrap: break-word !important;
            overflow-wrap: break-word !important;
        }}
        
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
        }}
        
        .download-links {{
            display: flex;
            justify-content: center;
            gap: 15px;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        
        .download-btn {{
            background: linear-gradient(135deg, #10b981 0%, #059669 100%);
            color: white;
            padding: 10px 20px;
            border-radius: 6px;
            text-decoration: none;
            font-weight: bold;
            display: inline-flex;
            align-items: center;
            gap: 8px;
        }}
        
        .download-btn:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
        }}
        
        @media (max-width: 768px) {{
            .ascii-art {{
                font-size: 5px;
            }}
            table {{
                display: block;
                overflow-x: auto;
            }}
            th, td {{
                padding: 8px 10px;
                font-size: 14px;
            }}
            .download-links {{
                flex-direction: column;
                align-items: center;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art"> █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝</div>
            </div>
            
            <h1 style="color: white; margin-bottom: 10px;">{scheme_display} MLST Batch Analysis</h1>
            <p style="color: #d1d5db; margin-bottom: 20px;">Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{total_samples}</div>
                <div class="stat-label">SAMPLES ANALYZED</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{assigned_st}</div>
                <div class="stat-label">STs ASSIGNED</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{unknown_st}</div>
                <div class="stat-label">UNKNOWN STs</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{high_confidence}</div>
                <div class="stat-label">HIGH CONFIDENCE</div>
            </div>
        </div>
        
        <div class="download-links">
            <a href="{scheme_display.lower()}_mlst_summary.csv" class="download-btn">
                📄 Download CSV Summary
            </a>
            <a href="{scheme_display.lower()}_mlst_summary.json" class="download-btn">
                🔧 Download JSON Summary
            </a>
'''
        
        # Add Excel download link only if available
        if self.has_excel_support:
            html_content += f'''            <a href="{scheme_display.lower()}_mlst_summary.xlsx" class="download-btn">
                📊 Download Excel Summary
            </a>
'''
        
        html_content += f'''        </div>
        
        <div class="report-section">
            <h2>📋 Detailed Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>#</th>
                        <th>Sample</th>
                        <th>ST</th>
                        <th>International Clone</th>
                        <th>Confidence</th>
                        <th>Allele Profile (Full)</th>
                    </tr>
                </thead>
                <tbody>
                    {table_rows}
                </tbody>
            </table>
        </div>
        
        <div class="footer">
            <p><strong>ACINETOSCOPE</strong> - A. baumannii MLST Analysis Module</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <p>For individual sample reports, navigate to respective sample directories.</p>
            <div style="margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1); border-radius: 8px; font-size: 12px;">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>
</body>
</html>'''
        
        html_file = output_dir / f"{self.scheme_display_names.get(scheme, scheme).lower()}_mlst_summary.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"  ✅ HTML summary: {html_file}")

def main():
    parser = argparse.ArgumentParser(
        description='AcinetoScope MLST Analyzer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single sample with Oxford scheme
  python mlst_module.py -i sample.fasta -o results/ -db db/ -sc bin/
  
  # Single sample with dual schemes (Oxford + Pasteur)
  python mlst_module.py -i sample.fasta -o results/ -db db/ -sc bin/ --dual
  
  # Batch processing with glob pattern
  python mlst_module.py -i "*.fasta" -o results/ -db db/ -sc bin/ --batch
  
  # Batch processing with directory
  python mlst_module.py -i /path/to/fasta_files/ -o results/ -db db/ -sc bin/ --batch
  
  # Dual-scheme with glob pattern (NEW FEATURE)
  python mlst_module.py -i "*.fasta" -o results/ -db db/ -sc bin/ --dual
  
  # Dual-scheme with directory
  python mlst_module.py -i /path/to/fasta_files/ -o results/ -db db/ -sc bin/ --dual
  
  # Specify Pasteur scheme
  python mlst_module.py -i sample.fasta -o results/ -db db/ -sc bin/ -s pasteur
        """
    )
    
    parser.add_argument('-i', '--input', required=True, 
                       help='Input FASTA file, directory, or glob pattern (e.g., "*.fasta", "data/*.fna")')
    parser.add_argument('-o', '--output-dir', required=True, 
                       help='Output directory for results')
    parser.add_argument('-db', '--database-dir', required=True,
                       help='Database directory')
    parser.add_argument('-sc', '--script-dir', required=True,
                       help='Script directory (contains mlst binary)')
    parser.add_argument('-s', '--scheme', default='oxford',
                       choices=['oxford', 'pasteur'],
                       help='MLST scheme: oxford or pasteur (default: oxford)')
    parser.add_argument('--batch', action='store_true',
                       help='Process multiple files with single scheme (supports glob patterns)')
    parser.add_argument('--dual', action='store_true',
                       help='Run both Oxford and Pasteur schemes (supports glob patterns and multiple files)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.batch and args.dual:
        print("❌ Error: --batch and --dual are mutually exclusive")
        print("   Use --batch for single scheme batch processing")
        print("   Use --dual for dual-scheme processing (single or multiple files)")
        sys.exit(1)
    
    # Map scheme name to database name
    scheme_map = {
        'oxford': 'abaumannii',
        'pasteur': 'abaumannii_2'
    }
    database_scheme = scheme_map[args.scheme]
    
    try:
        analyzer = AcinetoMLSTAnalyzer(
            database_dir=Path(args.database_dir),
            script_dir=Path(args.script_dir)
        )
    except FileNotFoundError as e:
        print(f"❌ {e}")
        print("   Please ensure mlst is installed and in PATH or script directory")
        sys.exit(1)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"""
    {'='*80}
    🧬 ACINETOSCOPE - A. baumannii MLST Analyzer
    {'='*80}
    Author: Brown Beckley
    Email: brownbeckley94@gmail.com
    GitHub: bbeckley-hub
    {'='*80}
    Input: {args.input}
    Output: {args.output_dir}
    Scheme: {args.scheme.upper()} ({database_scheme})
    Mode: {'Dual (with glob support)' if args.dual else 'Batch' if args.batch else 'Single'}
    {'='*80}
    """)
    
    if args.dual:
        # Find all FASTA files using glob pattern support
        fasta_files = analyzer.find_fasta_files(args.input)
        
        if not fasta_files:
            print(f"❌ No FASTA files found for input: {args.input}")
            sys.exit(1)
        
        if len(fasta_files) == 1:
            # Single file mode - use original dual-scheme method
            input_file = fasta_files[0]
            results = analyzer.run_mlst_dual_scheme(input_file, output_dir)
            print(f"\n{'='*80}")
            print(f"🎉 Dual-scheme MLST completed for {input_file.name}")
            print(f"📁 Results saved in: {output_dir}/mlst_results/")
            print(f"{'='*80}")
        else:
            # Multiple files - use new dual-scheme batch method
            results = analyzer.run_mlst_dual_batch(args.input, output_dir)
            print(f"\n{'='*80}")
            print(f"🎉 Dual-scheme MLST completed for {len(fasta_files)} samples")
            print(f"📁 Results saved in: {output_dir}/Dual_Scheme_MLST/")
            print(f"{'='*80}")
            
    elif args.batch:
        try:
            results = analyzer.run_mlst_batch(args.input, output_dir, database_scheme)
            print(f"\n{'='*80}")
            print(f"🎉 Batch MLST completed! Processed {len(results)} samples")
            print(f"📁 Results saved in: {output_dir}/")
            print(f"{'='*80}")
        except Exception as e:
            print(f"\n❌ Batch processing failed: {e}")
            sys.exit(1)
    else:
        # Single file mode
        input_file = Path(args.input)
        if input_file.exists() and input_file.is_file():
            # Create output directory for single sample
            sample_output_dir = output_dir / input_file.stem
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            
            result = analyzer.run_mlst_single(input_file, sample_output_dir, database_scheme)
            
            print(f"\n{'='*80}")
            print(f"🎉 MLST completed for {input_file.name}: ST{result.get('st', 'UNKNOWN')}")
            print(f"📁 Results saved in: {sample_output_dir}/")
            print(f"{'='*80}")
        else:
            # Try to find the file using glob patterns
            fasta_files = analyzer.find_fasta_files(args.input)
            if len(fasta_files) == 1:
                input_file = fasta_files[0]
                sample_output_dir = output_dir / input_file.stem
                sample_output_dir.mkdir(parents=True, exist_ok=True)
                
                result = analyzer.run_mlst_single(input_file, sample_output_dir, database_scheme)
                
                print(f"\n{'='*80}")
                print(f"🎉 MLST completed for {input_file.name}: ST{result.get('st', 'UNKNOWN')}")
                print(f"📁 Results saved in: {sample_output_dir}/")
                print(f"{'='*80}")
            elif len(fasta_files) > 1:
                print(f"❌ Found {len(fasta_files)} files for input: {args.input}")
                print(f"   Use --batch for multiple files or --dual for dual-scheme analysis")
                for i, f in enumerate(fasta_files[:5], 1):
                    print(f"   {i}. {f}")
                if len(fasta_files) > 5:
                    print(f"   ... and {len(fasta_files) - 5} more")
                sys.exit(1)
            else:
                print(f"❌ Input file not found: {args.input}")
                print(f"   Did you mean to use --batch for multiple files or --dual for dual-scheme?")
                sys.exit(1)

if __name__ == "__main__":
    main()