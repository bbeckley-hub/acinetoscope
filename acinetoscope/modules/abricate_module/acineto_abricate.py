#!/usr/bin/env python3
"""
AcinetoScope ABRicate Module - MAXIMUM SPEED VERSION
Comprehensive ABRicate analysis for Acinetobacter baumannii with HTML reporting
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2025-12-25
Send a quick mail for any issues or further explanations.
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any
import argparse
import re
from datetime import datetime
import psutil
import json
import pandas as pd
from collections import defaultdict, Counter

class AcinetoAbricateExecutor:
    """ABRicate executor for A. baumannii with comprehensive HTML reporting - MAXIMUM SPEED"""
    
    def __init__(self, cpus: int = None):
        # Setup logging FIRST
        self.logger = self._setup_logging()
        
        # Initialize available_ram before calculating cpus
        self.available_ram = self._get_available_ram()
        
        # Then calculate resources - MAXIMUM SPEED MODE
        self.cpus = self._calculate_optimal_cpus(cpus)
        
        # A. baumannii specific databases
        self.required_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot', 
            'plasmidfinder', 'megares', 'ecoh', 'bacmet2', 'ecoli_vf', 'victors'
        ]
        
        # CRITICAL CARBAPENEMASE genes for A. baumannii (🔴 HIGHEST PRIORITY)
        self.critical_carbapenemases = {
            'blaOXA-23', 'blaOXA-24', 'blaOXA-40', 'blaOXA-51', 'blaOXA-58', 'blaOXA-66', 'blaOXA-72', 
            'blaOXA-143', 'blaOXA-235', 'blaOXA-181', 'blaOXA-232', 'blaOXA-438', 'blaOXA-48',
            'blaKPC', 'blaNDM', 'blaVIM', 'blaIMP', 'blaSIM', 'blaGIM', 'blaSPM',
            'OXA-23', 'OXA-24', 'OXA-40', 'OXA-51', 'OXA-58', 'OXA-66', 'OXA-72', 'OXA-143', 'OXA-235',
            'KPC', 'NDM', 'VIM', 'IMP', 'SIM', 'GIM', 'SPM'
        }
        
        # CRITICAL ESBL genes for A. baumannii (🔴 CRITICAL)
        self.critical_esbls = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'blaPER', 'blaVEB', 'blaGES', 'blaBEL', 'blaBES',
            'CTX-M', 'SHV', 'TEM', 'PER', 'VEB', 'GES', 'BEL', 'BES',
            'blaADC', 'blaOXA-10', 'blaOXA-2', 'blaOXA-9',
            'ADC', 'OXA-10', 'OXA-2', 'OXA-9'
        }
        
        # CRITICAL Colistin resistance genes (🔴 CRITICAL)
        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'pmrA', 'pmrB', 'pmrC', 'lpxA', 'lpxC', 'lpxD', 'arnA', 'arnB', 'arnC', 'arnD', 'eptA'
        }
        
        # CRITICAL Aminoglycoside resistance (16S rRNA methyltransferases) (🔴 CRITICAL)
        self.critical_aminoglycoside = {
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH',
            'APH(3\')', 'APH(3\')-VI', 'APH(6)', 'AAC(3)', 'AAC(6\')', 'ANT(2")', 'ANT(3")', 'ANT(4")',
            'aacC1', 'aacC2', 'aacC4', 'aacA4', 'aacA7', 'aadA1', 'aadA2', 'aadA5', 'aadA7',
            'strA', 'strB', 'aphA1', 'aphA2', 'aphA3', 'aphA6', 'aac3', 'aac6', 'aadA', 'aadB'
        }
        
        # HIGH RISK resistance genes for A. baumannii (🟡 HIGH RISK)
        self.high_risk_resistance = {
            # Tetracycline resistance
            'tetA', 'tetB', 'tetC', 'tetD', 'tetE', 'tetG', 'tetH', 'tetK', 'tetL', 'tetM', 'tetO', 'tetQ', 'tetS', 'tetX',
            
            # Sulfonamide resistance
            'sul1', 'sul2', 'sul3', 'sul4',
            
            # Trimethoprim resistance
            'dfrA1', 'dfrA5', 'dfrA7', 'dfrA8', 'dfrA12', 'dfrA14', 'dfrA17', 'dfrA19', 'dfrA20', 'dfrA21',
            'dfrB1', 'dfrB2', 'dfrB3', 'dfrB4', 'dfrB5', 'dfrB6', 'dfrB7',
            
            # Chloramphenicol resistance
            'catA1', 'catA2', 'catB2', 'catB3', 'catB8', 'catI', 'catII', 'catIII',
            'cmlA', 'cmlA1', 'cmlA5', 'cmlA6', 'cmlA7', 'floR',
            
            # Macrolide resistance
            'ermA', 'ermB', 'ermC', 'ermF', 'ermG', 'ermX', 'ermY',
            'mphA', 'mphB', 'mphC', 'mphD', 'mphE', 'msrA', 'msrB', 'msrC', 'msrD',
            
            # Quinolone resistance
            'qnrA1', 'qnrB1', 'qnrB2', 'qnrB4', 'qnrB6', 'qnrB10', 'qnrB19', 'qnrS1', 'qnrS2', 'qnrVC1', 'qnrVC4',
            'aac(6\')-Ib-cr', 'qepA', 'qepA1', 'qepA2', 'qepA3', 'qepA4',
            
            # Fosfomycin resistance
            'fosA', 'fosB', 'fosC', 'fosX',
            
            # Rifampicin resistance
            'arr-2', 'arr-3', 'arr-4', 'arr-5', 'arr-6', 'arr-7',
            
            # Multidrug efflux pumps (CRITICAL for A. baumannii)
            'adeA', 'adeB', 'adeC', 'adeG', 'adeH', 'adeI', 'adeJ', 'adeK', 'adeL', 'adeM', 'adeN',
            'abeM', 'abeS', 'acrA', 'acrB', 'tolC', 'mexA', 'mexB', 'mexC', 'mexD', 'mexE', 'mexF',
            'oprM', 'oprN', 'oprJ', 'mdtA', 'mdtB', 'mdtC', 'emrA', 'emrB', 'emrD', 'emrE', 'emrK', 'emrY'
        }
        
        # VIRULENCE FACTORS for A. baumannii (🟢 VIRULENCE)
        self.virulence_genes = {
            # Biofilm formation
            'csuA', 'csuB', 'csuC', 'csuD', 'csuE', 'csuA/B', 'csuC/D/E',
            'bfmR', 'bfmS', 'ompA', 'bap', 'blaPER-1', 'epsA', 'ptk',
            
            # Siderophores and iron acquisition
            'basA', 'basB', 'basC', 'basD', 'basE', 'basF', 'basG', 'basH', 'basI', 'basJ',
            'barA', 'barB', 'bauA', 'bauB', 'bauC', 'bauD', 'bauE', 'bauF',
            'entA', 'entB', 'entC', 'entD', 'entE', 'entF', 'entS',
            'fhuA', 'fhuB', 'fhuC', 'fhuD', 'fhuE', 'fhuF',
            'fepA', 'fepB', 'fepC', 'fepD', 'fepE', 'fepG',
            
            # Quorum sensing
            'abaI', 'abaR', 'luxI', 'luxR',
            
            # Toxins and secretion systems
            'hlyA', 'hlyB', 'hlyC', 'hlyD', 'rtxA', 'rtxB', 'rtxC', 'rtxD', 'rtxE',
            'tssA', 'tssB', 'tssC', 'tssD', 'tssE', 'tssF', 'tssG', 'tssH', 'tssI', 'tssJ', 'tssK', 'tssL', 'tssM',
            'vgrG', 'vgrG1', 'vgrG2', 'vgrG3', 'vgrG4', 'vgrG5',
            
            # Phospholipases and proteases
            'plc1', 'plc2', 'plcD', 'lipA', 'lipB', 'lipC', 'lipH',
            'aprA', 'aprB', 'aprC', 'aprD', 'aprE', 'aprF',
            
            # Capsule and surface antigens
            'cpsA', 'cpsB', 'cpsC', 'cpsD', 'cpsE', 'cpsF', 'cpsG', 'cpsH', 'cpsI', 'cpsJ',
            'kpsM', 'kpsS', 'kpsT', 'kpsU', 'kpsV', 'kpsW',
            'wza', 'wzb', 'wzc', 'wzd', 'wze', 'wzf', 'wzg',
            
            # Adhesins
            'fimA', 'fimB', 'fimC', 'fimD', 'fimE', 'fimF', 'fimG', 'fimH',
            'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH',
            'sfmA', 'sfmB', 'sfmC', 'sfmD', 'sfmE', 'sfmF', 'sfmG', 'sfmH',
            
            # Type IV pili
            'pilA', 'pilB', 'pilC', 'pilD', 'pilE', 'pilF', 'pilG', 'pilH', 'pilI', 'pilJ', 'pilK', 'pilL', 'pilM', 'pilN', 'pilO', 'pilP', 'pilQ'
        }

        # Combined risk categories for better reporting
        self.critical_resistance_genes = self.critical_carbapenemases.union(
            self.critical_esbls
        ).union(self.critical_colistin).union(self.critical_aminoglycoside)
        
        self.critical_virulence_genes = {
            'hlyA', 'hlyB', 'hlyC', 'hlyD', 'rtxA', 'rtxB', 'rtxC', 'rtxD', 'rtxE'
        }
        
        # Science quotes
        self.science_quotes = [
            "The important thing is not to stop questioning. Curiosity has its own reason for existing. - Albert Einstein",
            "Science is not only a disciple of reason but also one of romance and passion. - Stephen Hawking", 
            "Somewhere, something incredible is waiting to be known. - Carl Sagan",
            "The good thing about science is that it's true whether or not you believe in it. - Neil deGrasse Tyson",
            "In science, there are no shortcuts to truth. - Karl Popper",
            "Science knows no country, because knowledge belongs to humanity. - Louis Pasteur",
            "The science of today is the technology of tomorrow. - Edward Teller",
            "Nothing in life is to be feared, it is only to be understood. - Marie Curie",
            "Research is what I'm doing when I don't know what I'm doing. - Wernher von Braun",
            "The universe is not required to be in harmony with human ambition. - Carl Sagan"
        ]
        
        # ASCII Art for AcinetoScope
        self.ascii_art = r"""
 █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
"""
        
        self.metadata = {
            "tool_name": "AcinetoScope ABRicate",
            "version": "1.0.0", 
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
    
    def _setup_logging(self):
        """Setup logging - must be called first in __init__"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def _get_available_ram(self) -> int:
        """Get available RAM in GB"""
        try:
            ram_gb = psutil.virtual_memory().available / (1024 ** 3)
            return ram_gb
        except Exception as e:
            self.logger.warning(f"Could not detect RAM: {e}")
            return 8  # Assume 8GB as fallback
    
    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        """Calculate optimal number of CPU cores for MAXIMUM SPEED"""
        if user_cpus is not None:
            self._log_resource_info(user_cpus)
            return user_cpus
            
        try:
            # Get total PHYSICAL CPU cores (not logical threads)
            total_physical_cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            
            # MAXIMUM SPEED RULES - AGGRESSIVE CPU USAGE
            if total_physical_cores <= 4:
                optimal_cpus = total_physical_cores  # Use ALL cores on small systems
            elif total_physical_cores <= 8:
                optimal_cpus = total_physical_cores - 1  # Use 7/8, 6/7, etc.
            elif total_physical_cores <= 16:
                optimal_cpus = max(8, total_physical_cores - 1)  # Use 15/16, 14/15, etc.
            elif total_physical_cores <= 32:
                optimal_cpus = max(16, total_physical_cores - 2)  # Use 30/32, 29/31, etc.
            else:
                optimal_cpus = min(32, int(total_physical_cores * 0.95))  # Use 95% on huge systems
            
            # Ensure at least 1 CPU and not more than available cores
            optimal_cpus = max(1, min(optimal_cpus, total_physical_cores))
            
            self._log_resource_info(optimal_cpus, total_physical_cores)
            return optimal_cpus
            
        except Exception as e:
            # Fallback to using all available cores for maximum speed
            self.logger.warning(f"Could not detect CPU cores, using maximum available: {e}")
            return os.cpu_count() or 4
    
    def _log_resource_info(self, cpus: int, total_cores: int = None):
        """Log resource allocation information"""
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            utilization = (cpus / total_cores) * 100
            self.logger.info(f"Using CPU cores: {cpus} ({utilization:.1f}% of available cores)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        
        # Performance recommendations - MAXIMUM SPEED FOCUS
        if cpus == 1:
            self.logger.info("💡 Performance: Single-core (max speed for 1-core systems)")
        elif cpus <= 4:
            self.logger.info("💡 Performance: Multi-core (max speed for small systems)")
        elif cpus <= 8:
            self.logger.info("💡 Performance: High-speed mode")
        else:
            self.logger.info("💡 Performance: MAXIMUM SPEED MODE 🚀")

    def check_abricate_installed(self) -> bool:
        """Check if ABRicate is installed and meets version requirements"""
        try:
            result = subprocess.run(['abricate', '--version'], 
                                  capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            self.logger.info("ABRicate version: %s", version_line)
            
            version_match = re.search(r'(\d+\.\d+\.\d+)', version_line)
            if version_match:
                version_str = version_match.group(1)
                if version_str >= "1.0.1":
                    self.logger.info("✓ ABRicate version meets requirement (>=1.2.0)")
                    return True
                else:
                    self.logger.error("ABRicate version too old: %s. Required >=1.2.0", version_str)
                    return False
            self.logger.info("✓ ABRicate installed (version check skipped)")
            return True
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("ABRicate not found. Please install with: conda install -c bioconda abricate")
            return False
    
    def setup_abricate_databases(self):
        """Setup ABRicate databases if they don't exist"""
        self.logger.info("Setting up ABRicate databases for A. baumannii analysis...")
        
        available_dbs = []
        missing_dbs = []
        
        try:
            # Check which databases exist
            check_result = subprocess.run(['abricate', '--list'], 
                                        capture_output=True, text=True, check=True)
            
            for db in self.required_databases:
                if db in check_result.stdout:
                    self.logger.info("✓ Database available: %s", db)
                    available_dbs.append(db)
                else:
                    self.logger.warning("Database not available: %s", db)
                    missing_dbs.append(db)
            
            # Setup missing databases
            for db in missing_dbs:
                self.logger.info("Attempting to setup database: %s", db)
                try:
                    result = subprocess.run(
                        ['abricate', '--setupdb', '--db', db],
                        capture_output=True, text=True, check=True
                    )
                    self.logger.info("✓ Database setup completed: %s", db)
                    available_dbs.append(db)
                except subprocess.CalledProcessError as e:
                    self.logger.error("Failed to setup database %s: %s", db, e.stderr)
            
            # Update required databases to only include available ones
            self.required_databases = available_dbs
            self.logger.info("Using databases: %s", ", ".join(self.required_databases))
            
        except subprocess.CalledProcessError as e:
            self.logger.error("Error checking ABRicate databases: %s", e.stderr)
        except Exception as e:
            self.logger.error("Unexpected error setting up databases: %s", e)
    
    def run_abricate_single_db(self, genome_file: str, database: str, output_dir: str) -> Dict[str, Any]:
        """Run ABRicate on a single genome with specific database"""
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"abricate_{database}.txt")
        
        cmd = [
            'abricate', 
            genome_file, 
            '--db', database,
            '--minid', '80',
            '--mincov', '80'
        ]
        
        self.logger.info("Running ABRicate: %s --db %s", genome_name, database)
        
        try:
            with open(output_file, 'w') as outfile:
                result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
            
            # Parse results for reporting
            hits = self._parse_abricate_output(output_file)
            
            # Create individual database HTML report
            self._create_database_html_report(genome_name, database, hits, output_dir)
            
            return {
                'database': database,
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'status': 'success'
            }
            
        except subprocess.CalledProcessError as e:
            self.logger.error("ABRicate failed for %s on %s: %s", database, genome_name, e.stderr)
            return {
                'database': database,
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed'
            }
    
    def _parse_abricate_output(self, abricate_file: str) -> List[Dict]:
        """Parse ABRicate output file - ROBUST VERSION that handles tabs in fields"""
        hits = []
        try:
            with open(abricate_file, 'r') as f:
                lines = f.readlines()
                
            if not lines:
                return hits
                
            # Find header line
            headers = []
            data_lines = []
            
            for line in lines:
                if line.startswith('#FILE') and not headers:
                    # This is the header line
                    headers = line.strip().replace('#', '').split('\t')
                elif line.strip() and not line.startswith('#'):
                    data_lines.append(line.strip())
            
            if not headers:
                self.logger.warning("No headers found in %s", abricate_file)
                return hits
            
            # Expected column count based on standard ABRicate output
            expected_columns = len(headers)
                
            # Parse data lines with robust tab handling
            for line_num, line in enumerate(data_lines, 1):
                # Split by tab but be careful about fields that contain tabs
                parts = line.split('\t')
                
                # Handle cases where there are more parts than headers due to tabs in fields
                if len(parts) > expected_columns:
                    # Combine extra fields into the last column (usually PRODUCT)
                    combined_parts = parts[:expected_columns-1]  # Take all but the last expected column
                    combined_parts.append('\t'.join(parts[expected_columns-1:]))  # Combine the rest into PRODUCT
                    parts = combined_parts
                elif len(parts) < expected_columns:
                    # Pad with empty strings if fewer columns
                    parts.extend([''] * (expected_columns - len(parts)))
                
                if len(parts) == expected_columns:
                    hit = {}
                    for i, header in enumerate(headers):
                        hit[header] = parts[i] if i < len(parts) else ''
                    
                    # Map to consistent field names
                    processed_hit = {
                        'file': hit.get('FILE', ''),
                        'sequence': hit.get('SEQUENCE', ''),
                        'start': hit.get('START', ''),
                        'end': hit.get('END', ''),
                        'strand': hit.get('STRAND', ''),
                        'gene': hit.get('GENE', ''),
                        'coverage': hit.get('COVERAGE', ''),
                        'coverage_map': hit.get('COVERAGE_MAP', ''),
                        'gaps': hit.get('GAPS', ''),
                        'coverage_percent': hit.get('%COVERAGE', ''),
                        'identity_percent': hit.get('%IDENTITY', ''),
                        'database': hit.get('DATABASE', ''),
                        'accession': hit.get('ACCESSION', ''),
                        'product': hit.get('PRODUCT', ''),
                        'resistance': hit.get('RESISTANCE', '')
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning("Line %d has %d parts, expected %d: %s", 
                                      line_num, len(parts), expected_columns, line[:10000] + "...")
                    
        except Exception as e:
            self.logger.error("Error parsing %s: %s", abricate_file, e)
            
        self.logger.info("Parsed %d hits from %s", len(hits), abricate_file)
        return hits
    
    def _create_database_html_report(self, genome_name: str, database: str, hits: List[Dict], output_dir: str):
        """Create individual HTML report for each database with AcinetoScope styling"""
        
        # JavaScript for rotating quotes
        quotes_js = f"""
        <script>
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            
            // Rotate every 10 seconds
            setInterval(rotateQuote, 10000);
            
            // Initial display
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
            }});
        </script>
        """
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope ABRicate - {database.upper()} Database</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
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
            max-width: 1200px;
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
        
        .card {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        
        .gene-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }}
        
        .gene-table th, .gene-table td {{
            padding: 15px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            font-weight: 600;
        }}
        
        tr:hover {{ background-color: #f8f9fa; }}
        
        .summary-stats {{
            display: flex;
            justify-content: space-around;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #fff;
        }}
        
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
        }}
        
        .footer a {{
            color: #667eea;
            text-decoration: none;
        }}
        
        .footer a:hover {{
            text-decoration: underline;
        }}
        
        .present {{ background-color: #d4edda; }}
        .high-risk {{ background-color: #fff3cd; }}
        .critical {{ background-color: #f8d7da; font-weight: bold; }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 AcinetoScope ABRicate - {database.upper()} Database</h1>
                <p style="color: #666; font-size: 1.2em;">Genome: {genome_name} | Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
            </div>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📊 Database Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{len(hits)}</p>
                </div>
                <div class="stat-card">
                    <h3>Database</h3>
                    <p style="font-size: 1.5em; margin: 0;">{database.upper()}</p>
                </div>
            </div>
        </div>
"""
        
        if hits:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🔍 Genes Detected</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Product</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Accession</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for hit in hits:
                # Determine row class based on gene risk
                row_class = "present"
                gene_base = hit['gene'].split('-')[0]  # Get base gene name (handle blaOXA-23)
                
                if any(crit_gene in gene_base for crit_gene in self.critical_resistance_genes):
                    row_class = "critical"
                elif any(vf_gene in gene_base for vf_gene in self.virulence_genes):
                    row_class = "high-risk"
                
                # Truncate very long product descriptions for display
                product_display = hit['product']
                if len(product_display) > 500:
                    product_display = product_display[:1000] + "..."
                
                html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{hit['gene']}</strong></td>
                        <td title="{hit['product']}">{product_display}</td>
                        <td>{hit['coverage_percent']}%</td>
                        <td>{hit['identity_percent']}%</td>
                        <td>{hit['accession']}</td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        else:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">✅ No Genes Detected</h2>
            <p>No significant hits found in the {database.upper()} database.</p>
        </div>
"""
        
        html_content += f"""
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using AcinetoScope ABRicate v1.2.0
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write individual database HTML report
        html_file = os.path.join(output_dir, f"abricate_{database}_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("Individual database report: %s", html_file)
    
    def analyze_acineto_resistance(self, all_hits: List[Dict]) -> Dict[str, Any]:
        """Enhanced A. baumannii resistance analysis with comprehensive risk assessment"""
        analysis = {
            'carbapenemase_status': 'negative',
            'esbl_status': 'negative',
            'colistin_resistance': 'negative',
            'critical_carbapenemase_genes': [],
            'critical_esbl_genes': [],
            'critical_colistin_genes': [],
            'critical_aminoglycoside_genes': [],
            'high_risk_resistance_genes': [],
            'critical_virulence_genes': [],
            'high_risk_virulence_genes': [],
            'moderate_risk_genes': [],
            'resistance_classes': {},
            'total_critical_resistance': 0,
            'total_high_risk_resistance': 0,
            'total_critical_virulence': 0,
            'total_high_risk_virulence': 0,
            'total_hits': len(all_hits)
        }
        
        for hit in all_hits:
            gene = hit['gene']
            gene_base = gene.split('-')[0] if '-' in gene else gene
            gene_lower = gene.lower()
            
            # Check for CRITICAL CARBAPENEMASE patterns (🔴 HIGHEST PRIORITY for A. baumannii)
            if any(carba in gene_lower for carba in ['oxa', 'kpc', 'ndm', 'vim', 'imp', 'sim']):
                if any(crit_gene in gene_lower for crit_gene in [g.lower() for g in self.critical_carbapenemases]):
                    analysis['carbapenemase_status'] = 'positive'
                    risk_level = 'CARBAPENEMASE'
                    
                    analysis['critical_carbapenemase_genes'].append({
                        'gene': gene,
                        'product': hit['product'],
                        'database': hit['database'],
                        'coverage': hit['coverage_percent'],
                        'identity': hit['identity_percent'],
                        'risk_level': risk_level
                    })
            
            # Check for CRITICAL ESBL patterns
            elif any(esbl in gene_lower for esbl in ['ctx', 'shv', 'tem', 'per', 'veb', 'ges']):
                if any(crit_gene in gene_lower for crit_gene in [g.lower() for g in self.critical_esbls]):
                    analysis['esbl_status'] = 'positive'
                    risk_level = 'ESBL'
                    
                    analysis['critical_esbl_genes'].append({
                        'gene': gene,
                        'product': hit['product'],
                        'database': hit['database'],
                        'coverage': hit['coverage_percent'],
                        'identity': hit['identity_percent'],
                        'risk_level': risk_level
                    })
            
            # Check for CRITICAL COLISTIN resistance
            elif any(mcr in gene_lower for mcr in ['mcr']):
                if any(crit_gene in gene_lower for crit_gene in [g.lower() for g in self.critical_colistin]):
                    analysis['colistin_resistance'] = 'positive'
                    risk_level = 'COLISTIN-RES'
                    
                    analysis['critical_colistin_genes'].append({
                        'gene': gene,
                        'product': hit['product'],
                        'database': hit['database'],
                        'coverage': hit['coverage_percent'],
                        'identity': hit['identity_percent'],
                        'risk_level': risk_level
                    })
            
            # Check for CRITICAL AMINOGLYCOSIDE resistance
            elif any(ag in gene_lower for ag in ['arm', 'rmt', 'aph', 'aac', 'aad', 'str', 'ant']):
                if any(crit_gene in gene_lower for crit_gene in [g.lower() for g in self.critical_aminoglycoside]):
                    analysis['critical_aminoglycoside_genes'].append({
                        'gene': gene,
                        'product': hit['product'],
                        'database': hit['database'],
                        'coverage': hit['coverage_percent'],
                        'identity': hit['identity_percent'],
                        'risk_level': 'CRITICAL_AMINOGLYCOSIDE'
                    })
            
            # Check for HIGH RISK resistance genes
            elif any(hr_gene in gene_lower for hr_gene in [g.lower() for g in self.high_risk_resistance]):
                analysis['high_risk_resistance_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': 'HIGH-RISK'
                })
            
            # Check for CRITICAL virulence genes
            if any(crit_vf in gene_base for crit_vf in self.critical_virulence_genes):
                analysis['critical_virulence_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': 'CRITICAL-VIRULENCE'
                })
            
            # Check for HIGH RISK virulence genes
            elif any(vf_gene in gene_lower for vf_gene in [g.lower() for g in self.virulence_genes]):
                # Check if it's high risk virulence (not moderate)
                is_high_risk_virulence = any(hr_vf in gene_lower for hr_vf in [
                    'hly', 'rtx', 'tss', 'vgr', 'plc', 'lip', 'apr'
                ])
                
                if is_high_risk_virulence:
                    analysis['high_risk_virulence_genes'].append({
                        'gene': gene,
                        'product': hit['product'],
                        'database': hit['database'],
                        'coverage': hit['coverage_percent'],
                        'identity': hit['identity_percent'],
                        'risk_level': 'HIGH-VIRULENCE'
                    })
                else:
                    analysis['moderate_risk_genes'].append({
                        'gene': gene,
                        'product': hit['product'],
                        'database': hit['database'],
                        'coverage': hit['coverage_percent'],
                        'identity': hit['identity_percent'],
                        'risk_level': 'MODERATE'
                    })
            
            # Track resistance classes
            resistance_class = self._classify_resistance(hit['product'], gene)
            if resistance_class:
                if resistance_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][resistance_class] = []
                if gene not in [g['gene'] for g in analysis['resistance_classes'][resistance_class]]:
                    analysis['resistance_classes'][resistance_class].append({
                        'gene': gene,
                        'product': hit['product']
                    })
        
        # Calculate totals
        analysis['total_critical_resistance'] = (
            len(analysis['critical_carbapenemase_genes']) +
            len(analysis['critical_esbl_genes']) +
            len(analysis['critical_colistin_genes']) +
            len(analysis['critical_aminoglycoside_genes'])
        )
        analysis['total_high_risk_resistance'] = len(analysis['high_risk_resistance_genes'])
        analysis['total_critical_virulence'] = len(analysis['critical_virulence_genes'])
        analysis['total_high_risk_virulence'] = len(analysis['high_risk_virulence_genes'])
        
        return analysis

    def _classify_resistance(self, product: str, gene: str) -> str:
        """Enhanced resistance classification for A. baumannii"""
        product_lower = product.lower()
        gene_lower = gene.lower()
        
        if any(term in product_lower or term in gene_lower for term in ['carbapenem', 'oxa', 'kpc', 'ndm', 'vim', 'imp', 'sim']):
            return 'Carbapenem resistance'
        elif any(term in product_lower or term in gene_lower for term in ['beta-lactam', 'esbl', 'ctx', 'shv', 'tem', 'per', 'veb', 'ges']):
            return 'Beta-lactam resistance'
        elif any(term in product_lower or term in gene_lower for term in ['aminoglycoside', 'aac', 'aad', 'aph', 'str', 'arm', 'rmt']):
            return 'Aminoglycoside resistance'
        elif any(term in product_lower or term in gene_lower for term in ['tetracycline', 'tet']):
            return 'Tetracycline resistance'
        elif any(term in product_lower or term in gene_lower for term in ['sulfonamide', 'sul']):
            return 'Sulfonamide resistance'
        elif any(term in product_lower or term in gene_lower for term in ['trimethoprim', 'dfr']):
            return 'Trimethoprim resistance'
        elif any(term in product_lower or term in gene_lower for term in ['chloramphenicol', 'cat', 'cml', 'flor']):
            return 'Chloramphenicol resistance'
        elif any(term in product_lower or term in gene_lower for term in ['macrolide', 'erm', 'mph', 'msr']):
            return 'Macrolide resistance'
        elif any(term in product_lower or term in gene_lower for term in ['quinolone', 'qnr', 'qep']):
            return 'Quinolone resistance'
        elif any(term in product_lower or term in gene_lower for term in ['colistin', 'polymyxin', 'mcr']):
            return 'Polymyxin resistance'
        elif any(term in product_lower or term in gene_lower for term in ['fosfomycin', 'fos']):
            return 'Fosfomycin resistance'
        elif any(term in product_lower or term in gene_lower for term in ['rifampicin', 'arr']):
            return 'Rifampicin resistance'
        elif any(term in product_lower or term in gene_lower for term in ['efflux', 'ade', 'abe', 'mex', 'acr', 'emr']):
            return 'Efflux pumps'
        elif any(term in product_lower or term in gene_lower for term in ['virulence', 'toxin', 'hemolysin', 'siderophore']):
            return 'Virulence factors'
        else:
            return 'Other resistance'
    
    def create_comprehensive_html_report(self, genome_name: str, results: Dict, output_dir: str):
        """Create comprehensive HTML report for A. baumannii with AcinetoScope styling"""
        
        # Collect all hits
        all_hits = []
        for db_result in results.values():
            all_hits.extend(db_result['hits'])
        
        # Analyze A. baumannii resistance
        analysis = self.analyze_acineto_resistance(all_hits)
        
        # JavaScript for rotating quotes
        quotes_js = f"""
        <script>
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            
            // Rotate every 10 seconds
            setInterval(rotateQuote, 10000);
            
            // Initial display
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
            }});
        </script>
        """
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope ABRicate Analysis Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
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
            max-width: 1400px;
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
        
        .card {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        
        .gene-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }}
        
        .gene-table th, .gene-table td {{
            padding: 15px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            font-weight: 600;
        }}
        
        tr:hover {{ background-color: #f8f9fa; }}
        
        .success {{ color: #28a745; font-weight: 600; }}
        .warning {{ color: #ffc107; font-weight: 600; }}
        .error {{ color: #dc3545; font-weight: 600; }}
        .critical {{ background-color: #f8d7da; font-weight: bold; }}
        .high-risk {{ background-color: #fff3cd; }}
        
        .summary-stats {{
            display: flex;
            justify-content: space-around;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #fff;
        }}
        
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
        }}
        
        .footer a {{
            color: #667eea;
            text-decoration: none;
        }}
        
        .footer a:hover {{
            text-decoration: underline;
        }}
        
        .risk-badge {{
            display: inline-block;
            background: #dc3545;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        
        .warning-badge {{
            display: inline-block;
            background: #ffc107;
            color: black;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        
        .safe-badge {{
            display: inline-block;
            background: #28a745;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 AcinetoScope ABRicate Analysis Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive A. baumannii Antimicrobial Resistance & Virulence Analysis</p>
            </div>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📊 A. baumannii AMR Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_hits']}</p>
                </div>
                <div class="stat-card">
                    <h3>Critical Resistance</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_critical_resistance']}</p>
                </div>
                <div class="stat-card">
                    <h3>Virulence Factors</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_high_risk_virulence'] + analysis['total_critical_virulence']}</p>
                </div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
        </div>
"""
        
        # Critical resistance alerts
        critical_alerts = []
        if analysis['carbapenemase_status'] == 'positive':
            critical_alerts.append("🔴 CARBAPENEMASE DETECTED")
        if analysis['esbl_status'] == 'positive':
            critical_alerts.append("🟡 ESBL DETECTED")
        if analysis['colistin_resistance'] == 'positive':
            critical_alerts.append("🔴 COLISTIN RESISTANCE DETECTED")
        
        if critical_alerts:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #dc3545;">
            <h2 style="color: #dc3545;">⚠️ CRITICAL RESISTANCE ALERTS</h2>
            <div style="margin: 10px 0;">
"""
            for alert in critical_alerts:
                html_content += f'<span class="risk-badge">{alert}</span>'
            html_content += """
            </div>
        </div>
"""
        
        # Resistance classes summary
        if analysis['resistance_classes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🧪 Resistance Classes Detected</h2>
            <div style="margin: 20px 0;">
"""
            
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join([g['gene'] for g in genes])
                html_content += f"""
                <div style="margin: 10px 0; padding: 10px; background: #f8f9fa; border-radius: 8px;">
                    <strong style="color: #3b82f6;">{class_name}</strong> ({len(genes)} genes)
                    <br><span style="color: #666; font-size: 0.9em;">{gene_list}</span>
                </div>
"""
            
            html_content += "</div></div>"
        
        # Critical carbapenemase genes table
        if analysis['critical_carbapenemase_genes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🔴 CRITICAL Carbapenemase Genes</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Product</th>
                        <th>Database</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Risk Level</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for gene_info in analysis['critical_carbapenemase_genes']:
                product_display = gene_info['product']
                if len(product_display) > 100:
                    product_display = gene_info['product'][:97] + "..."
                
                html_content += f"""
                    <tr class="critical">
                        <td><strong>{gene_info['gene']}</strong></td>
                        <td title="{gene_info['product']}">{product_display}</td>
                        <td>{gene_info['database']}</td>
                        <td>{gene_info['coverage']}%</td>
                        <td>{gene_info['identity']}%</td>
                        <td><span class="risk-badge">{gene_info['risk_level']}</span></td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        
        # Critical ESBL genes table
        if analysis['critical_esbl_genes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🟡 Critical ESBL Genes</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Product</th>
                        <th>Database</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Risk Level</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for gene_info in analysis['critical_esbl_genes']:
                product_display = gene_info['product']
                if len(product_display) > 100:
                    product_display = gene_info['product'][:97] + "..."
                
                html_content += f"""
                    <tr class="critical">
                        <td><strong>{gene_info['gene']}</strong></td>
                        <td title="{gene_info['product']}">{product_display}</td>
                        <td>{gene_info['database']}</td>
                        <td>{gene_info['coverage']}%</td>
                        <td>{gene_info['identity']}%</td>
                        <td><span class="warning-badge">{gene_info['risk_level']}</span></td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        
        # Critical colistin resistance genes table
        if analysis['critical_colistin_genes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🔴 Colistin Resistance Genes</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Product</th>
                        <th>Database</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Risk Level</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for gene_info in analysis['critical_colistin_genes']:
                product_display = gene_info['product']
                if len(product_display) > 100:
                    product_display = gene_info['product'][:97] + "..."
                
                html_content += f"""
                    <tr class="critical">
                        <td><strong>{gene_info['gene']}</strong></td>
                        <td title="{gene_info['product']}">{product_display}</td>
                        <td>{gene_info['database']}</td>
                        <td>{gene_info['coverage']}%</td>
                        <td>{gene_info['identity']}%</td>
                        <td><span class="risk-badge">{gene_info['risk_level']}</span></td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        
        # High-risk resistance genes table
        if analysis['high_risk_resistance_genes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🟡 High-Risk Resistance Genes</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Product</th>
                        <th>Database</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Risk Level</th>
                    </tr>
                </thead>
                <tbody>
"""
            
            for gene_info in analysis['high_risk_resistance_genes']:
                product_display = gene_info['product']
                if len(product_display) > 1000:
                    product_display = gene_info['product'][:500] + "..."
                
                html_content += f"""
                    <tr class="high-risk">
                        <td><strong>{gene_info['gene']}</strong></td>
                        <td title="{gene_info['product']}">{product_display}</td>
                        <td>{gene_info['database']}</td>
                        <td>{gene_info['coverage']}%</td>
                        <td>{gene_info['identity']}%</td>
                        <td><span class="warning-badge">{gene_info['risk_level']}</span></td>
                    </tr>
"""
            
            html_content += """
                </tbody>
            </table>
        </div>
"""
        
        # Database summary
        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🗃️ Database Results Summary</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Database</th>
                        <th>Hits</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>
"""
        
        for db, result in results.items():
            status_icon = "✅" if result['status'] == 'success' else "❌"
            html_content += f"""
                    <tr>
                        <td>{db}</td>
                        <td>{result['hit_count']}</td>
                        <td>{status_icon} {result['status']}</td>
                    </tr>
"""
        
        html_content += """
                </tbody>
            </table>
        </div>
        
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using AcinetoScope (ABRicate v1.2.0)
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write comprehensive HTML report
        html_file = os.path.join(output_dir, f"{genome_name}_comprehensive_abricate_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("Comprehensive A. baumannii HTML report generated: %s", html_file)
    
    def create_database_summaries(self, all_results: Dict[str, Any], output_base: str):
        """Create ABRicate summary files and HTML reports for each database across all genomes"""
        self.logger.info("Creating A. baumannii database summary files and HTML reports...")
        
        # Group results by database
        db_results = {}
        for genome_name, genome_result in all_results.items():
            for db, db_result in genome_result['results'].items():
                if db not in db_results:
                    db_results[db] = []
                
                # Add hits with genome name
                for hit in db_result['hits']:
                    hit_with_genome = hit.copy()
                    hit_with_genome['genome'] = genome_name
                    db_results[db].append(hit_with_genome)
        
        # Create summary files for each database
        for db, hits in db_results.items():
            if hits:
                # Create TSV summary
                summary_file = os.path.join(output_base, f"acineto_{db}_abricate_summary.tsv")
                
                # Get headers from first hit
                headers = list(hits[0].keys())
                
                with open(summary_file, 'w') as f:
                    # Write header
                    f.write('\t'.join(headers) + '\n')
                    
                    # Write data
                    for hit in hits:
                        row = [str(hit.get(header, '')) for header in headers]
                        f.write('\t'.join(row) + '\n')
                
                self.logger.info("✓ Created %s summary: %s (%d hits)", db, summary_file, len(hits))
                
                # Create JSON summary for this database
                self._create_database_json_summary(db, hits, output_base)
                
                # Create HTML summary report for this database
                self._create_database_summary_html(db, hits, output_base)
            else:
                self.logger.info("No hits for database %s, skipping summary", db)
    
    def _create_database_json_summary(self, database: str, hits: List[Dict], output_base: str):
        """Create JSON summary report for a specific database across all genomes"""
        
        # Calculate statistics
        unique_genomes = list(set(hit['genome'] for hit in hits))
        unique_genes = list(set(hit['gene'] for hit in hits))
        
        # Calculate gene frequency
        gene_frequency = {}
        gene_details = {}
        for hit in hits:
            gene = hit['gene']
            if gene not in gene_frequency:
                gene_frequency[gene] = set()
                gene_details[gene] = {
                    'product': hit['product'],
                    'databases': set(),
                    'avg_coverage': 0,
                    'avg_identity': 0
                }
            gene_frequency[gene].add(hit['genome'])
            gene_details[gene]['databases'].add(hit['database'])
            
            # Track coverage and identity for averaging
            if 'coverage_sum' not in gene_details[gene]:
                gene_details[gene]['coverage_sum'] = 0
                gene_details[gene]['identity_sum'] = 0
                gene_details[gene]['count'] = 0
            
            try:
                coverage = float(hit['coverage_percent'].replace('%', ''))
                identity = float(hit['identity_percent'].replace('%', ''))
                gene_details[gene]['coverage_sum'] += coverage
                gene_details[gene]['identity_sum'] += identity
                gene_details[gene]['count'] += 1
            except:
                pass
        
        # Calculate averages
        for gene in gene_details:
            if gene_details[gene]['count'] > 0:
                gene_details[gene]['avg_coverage'] = round(gene_details[gene]['coverage_sum'] / gene_details[gene]['count'], 2)
                gene_details[gene]['avg_identity'] = round(gene_details[gene]['identity_sum'] / gene_details[gene]['count'], 2)
            del gene_details[gene]['coverage_sum']
            del gene_details[gene]['identity_sum']
            del gene_details[gene]['count']
            gene_details[gene]['databases'] = list(gene_details[gene]['databases'])
        
        # Calculate per-genome statistics
        genomes_summary = {}
        for genome in unique_genomes:
            genome_hits = [h for h in hits if h['genome'] == genome]
            genomes_summary[genome] = {
                'total_hits': len(genome_hits),
                'unique_genes': len(set(h['gene'] for h in genome_hits)),
                'genes': list(set(h['gene'] for h in genome_hits))
            }
        
        # Risk analysis for this database
        risk_genes = {
            'critical_carbapenemase': [],
            'critical_esbl': [],
            'critical_colistin': [],
            'critical_aminoglycoside': [],
            'high_risk_resistance': [],
            'virulence_factors': []
        }
        
        for gene in unique_genes:
            gene_base = gene.split('-')[0] if '-' in gene else gene
            gene_lower = gene.lower()
            
            # Check critical carbapenemases
            if any(carba in gene_lower for carba in ['oxa', 'kpc', 'ndm', 'vim', 'imp', 'sim']):
                if any(crit in gene_lower for crit in [c.lower() for c in self.critical_carbapenemases]):
                    risk_genes['critical_carbapenemase'].append(gene)
            
            # Check critical ESBLs
            elif any(esbl in gene_lower for esbl in ['ctx', 'shv', 'tem', 'per', 'veb', 'ges']):
                if any(crit in gene_lower for crit in [c.lower() for c in self.critical_esbls]):
                    risk_genes['critical_esbl'].append(gene)
            
            # Check critical colistin
            elif any(mcr in gene_lower for mcr in ['mcr']):
                if any(crit in gene_lower for crit in [c.lower() for c in self.critical_colistin]):
                    risk_genes['critical_colistin'].append(gene)
            
            # Check critical aminoglycoside
            elif any(ag in gene_lower for ag in ['arm', 'rmt', 'aph', 'aac', 'aad', 'str', 'ant']):
                if any(crit in gene_lower for crit in [c.lower() for c in self.critical_aminoglycoside]):
                    risk_genes['critical_aminoglycoside'].append(gene)
            
            # Check high risk resistance
            elif any(hr in gene_lower for hr in [g.lower() for g in self.high_risk_resistance]):
                risk_genes['high_risk_resistance'].append(gene)
            
            # Check virulence factors
            elif any(vf in gene_lower for vf in [g.lower() for g in self.virulence_genes]):
                risk_genes['virulence_factors'].append(gene)
        
        # Build JSON report
        json_report = {
            "metadata": {
                "database_name": database,
                "analysis_type": "ABRicate Antimicrobial Resistance & Virulence",
                "total_genomes": len(unique_genomes),
                "total_hits": len(hits),
                "analysis_date": self.metadata['analysis_date'],
                "tool_version": self.metadata['version']
            },
            "summary_statistics": {
                "unique_genes": len(unique_genes),
                "unique_genomes": len(unique_genomes),
                "hits_per_genome": round(len(hits) / len(unique_genomes), 2) if unique_genomes else 0,
                "risk_analysis": {
                    "critical_carbapenemase_genes": len(risk_genes['critical_carbapenemase']),
                    "critical_esbl_genes": len(risk_genes['critical_esbl']),
                    "critical_colistin_genes": len(risk_genes['critical_colistin']),
                    "critical_aminoglycoside_genes": len(risk_genes['critical_aminoglycoside']),
                    "high_risk_resistance_genes": len(risk_genes['high_risk_resistance']),
                    "virulence_factors": len(risk_genes['virulence_factors'])
                }
            },
            "gene_analysis": {
                "most_frequent_genes": dict(sorted(
                    {gene: len(genomes) for gene, genomes in gene_frequency.items()}.items(),
                    key=lambda x: x[1],
                    reverse=True
                )),
                "gene_details": gene_details
            },
            "genome_analysis": genomes_summary,
            "risk_genes": risk_genes,
            "raw_hits_summary": {
                "total_hits": len(hits),
                "hits_by_genome": {genome: len([h for h in hits if h['genome'] == genome]) for genome in unique_genomes}
            }
        }
        
        # Write JSON report
        json_file = os.path.join(output_base, f"acineto_{database}_summary.json")
        with open(json_file, 'w') as f:
            json.dump(json_report, f, indent=4, default=str)
        
        self.logger.info("✓ Created JSON summary: %s", json_file)
    
    def create_master_json_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create master JSON summary containing everything from all databases and all genomes"""
        self.logger.info("Creating master JSON summary for all databases and genomes...")
        
        master_report = {
            "metadata": {
                **self.metadata,
                "analysis_type": "A. baumannii Comprehensive Antimicrobial Resistance & Virulence Analysis",
                "databases_analyzed": self.required_databases,
                "total_genomes": len(all_results),
                "analysis_parameters": {
                    "minimum_identity": 80,
                    "minimum_coverage": 80
                }
            },
            "overall_summary": {
                "total_genomes": len(all_results),
                "total_hits": 0,
                "critical_carbapenemase_genes": 0,
                "critical_esbl_genes": 0,
                "critical_colistin_genes": 0,
                "critical_aminoglycoside_genes": 0,
                "carbapenemase_positive_genomes": 0,
                "esbl_positive_genomes": 0,
                "colistin_resistant_genomes": 0
            },
            "database_summaries": {},
            "genome_summaries": {},
            "risk_assessment": {
                "critical_carbapenemase_genes_found": [],
                "critical_esbl_genes_found": [],
                "critical_colistin_genes_found": [],
                "critical_aminoglycoside_genes_found": [],
                "high_risk_genes_found": [],
                "virulence_genes_found": [],
                "resistance_classes_found": defaultdict(list)
            },
            "detailed_results": {}
        }
        
        # Process each genome
        for genome_name, genome_result in all_results.items():
            # Collect all hits for this genome
            all_genome_hits = []
            for db_result in genome_result['results'].values():
                all_genome_hits.extend(db_result['hits'])
            
            # Analyze resistance for this genome
            analysis = self.analyze_acineto_resistance(all_genome_hits)
            
            # Update overall summary
            master_report["overall_summary"]["total_hits"] += analysis['total_hits']
            master_report["overall_summary"]["critical_carbapenemase_genes"] += len(analysis['critical_carbapenemase_genes'])
            master_report["overall_summary"]["critical_esbl_genes"] += len(analysis['critical_esbl_genes'])
            master_report["overall_summary"]["critical_colistin_genes"] += len(analysis['critical_colistin_genes'])
            master_report["overall_summary"]["critical_aminoglycoside_genes"] += len(analysis['critical_aminoglycoside_genes'])
            
            if analysis['carbapenemase_status'] == 'positive':
                master_report["overall_summary"]["carbapenemase_positive_genomes"] += 1
            if analysis['esbl_status'] == 'positive':
                master_report["overall_summary"]["esbl_positive_genomes"] += 1
            if analysis['colistin_resistance'] == 'positive':
                master_report["overall_summary"]["colistin_resistant_genomes"] += 1
            
            # Add to genome summaries
            master_report["genome_summaries"][genome_name] = {
                "total_hits": analysis['total_hits'],
                "critical_carbapenemase": len(analysis['critical_carbapenemase_genes']),
                "critical_esbl": len(analysis['critical_esbl_genes']),
                "critical_colistin": len(analysis['critical_colistin_genes']),
                "critical_aminoglycoside": len(analysis['critical_aminoglycoside_genes']),
                "high_risk_genes": analysis['total_high_risk_resistance'],
                "virulence_factors": analysis['total_high_risk_virulence'] + analysis['total_critical_virulence'],
                "carbapenemase_status": analysis['carbapenemase_status'],
                "esbl_status": analysis['esbl_status'],
                "colistin_resistance": analysis['colistin_resistance'],
                "database_hits": {db: result['hit_count'] for db, result in genome_result['results'].items()}
            }
            
            # Update risk assessment
            for gene_info in analysis['critical_carbapenemase_genes']:
                if gene_info['gene'] not in master_report["risk_assessment"]["critical_carbapenemase_genes_found"]:
                    master_report["risk_assessment"]["critical_carbapenemase_genes_found"].append(gene_info['gene'])
            
            for gene_info in analysis['critical_esbl_genes']:
                if gene_info['gene'] not in master_report["risk_assessment"]["critical_esbl_genes_found"]:
                    master_report["risk_assessment"]["critical_esbl_genes_found"].append(gene_info['gene'])
            
            for gene_info in analysis['critical_colistin_genes']:
                if gene_info['gene'] not in master_report["risk_assessment"]["critical_colistin_genes_found"]:
                    master_report["risk_assessment"]["critical_colistin_genes_found"].append(gene_info['gene'])
            
            for gene_info in analysis['critical_aminoglycoside_genes']:
                if gene_info['gene'] not in master_report["risk_assessment"]["critical_aminoglycoside_genes_found"]:
                    master_report["risk_assessment"]["critical_aminoglycoside_genes_found"].append(gene_info['gene'])
            
            for gene_info in analysis['high_risk_resistance_genes']:
                if gene_info['gene'] not in master_report["risk_assessment"]["high_risk_genes_found"]:
                    master_report["risk_assessment"]["high_risk_genes_found"].append(gene_info['gene'])
            
            for gene_info in analysis['high_risk_virulence_genes'] + analysis['critical_virulence_genes']:
                if gene_info['gene'] not in master_report["risk_assessment"]["virulence_genes_found"]:
                    master_report["risk_assessment"]["virulence_genes_found"].append(gene_info['gene'])
            
            # Update resistance classes
            for class_name, genes in analysis['resistance_classes'].items():
                for gene_info in genes:
                    if gene_info['gene'] not in [g['gene'] for g in master_report["risk_assessment"]["resistance_classes_found"][class_name]]:
                        master_report["risk_assessment"]["resistance_classes_found"][class_name].append(gene_info)
            
            # Add detailed results
            master_report["detailed_results"][genome_name] = {
                "total_hits": genome_result['total_hits'],
                "database_results": genome_result['results']
            }
        
        # Create database summaries
        for db in self.required_databases:
            # Collect all hits for this database across all genomes
            db_hits = []
            for genome_name, genome_result in all_results.items():
                if db in genome_result['results']:
                    for hit in genome_result['results'][db]['hits']:
                        hit_with_genome = hit.copy()
                        hit_with_genome['genome'] = genome_name
                        db_hits.append(hit_with_genome)
            
            if db_hits:
                unique_genes = len(set(h['gene'] for h in db_hits))
                unique_genomes = len(set(h['genome'] for h in db_hits))
                
                master_report["database_summaries"][db] = {
                    "total_hits": len(db_hits),
                    "unique_genes": unique_genes,
                    "unique_genomes": unique_genomes,
                    "hits_per_genome": round(len(db_hits) / unique_genomes, 2) if unique_genomes else 0,
                    "most_frequent_genes": dict(Counter([h['gene'] for h in db_hits]).most_common(10))
                }
            else:
                master_report["database_summaries"][db] = {
                    "total_hits": 0,
                    "unique_genes": 0,
                    "unique_genomes": 0,
                    "hits_per_genome": 0,
                    "most_frequent_genes": {}
                }
        
        # Sort and clean up
        for key in ["critical_carbapenemase_genes_found", "critical_esbl_genes_found", 
                   "critical_colistin_genes_found", "critical_aminoglycoside_genes_found",
                   "high_risk_genes_found", "virulence_genes_found"]:
            master_report["risk_assessment"][key].sort()
        
        master_report["risk_assessment"]["resistance_classes_found"] = dict(master_report["risk_assessment"]["resistance_classes_found"])
        
        # Add performance metrics
        master_report["performance_metrics"] = {
            "cpu_cores_used": self.cpus,
            "available_ram_gb": round(self.available_ram, 2),
            "processing_mode": "MAXIMUM SPEED" if self.cpus > 8 else "STANDARD"
        }
        
        # Add Excel summary
        try:
            # Create DataFrame for Excel export
            excel_data = []
            for genome_name, genome_result in all_results.items():
                all_genome_hits = []
                for db_result in genome_result['results'].values():
                    all_genome_hits.extend(db_result['hits'])
                
                analysis = self.analyze_acineto_resistance(all_genome_hits)
                
                excel_data.append({
                    'Genome': genome_name,
                    'Total_Hits': genome_result['total_hits'],
                    'Critical_Carbapenemase': len(analysis['critical_carbapenemase_genes']),
                    'Critical_ESBL': len(analysis['critical_esbl_genes']),
                    'Critical_Colistin': len(analysis['critical_colistin_genes']),
                    'Critical_Aminoglycoside': len(analysis['critical_aminoglycoside_genes']),
                    'High_Risk_Resistance': len(analysis['high_risk_resistance_genes']),
                    'Virulence_Factors': len(analysis['high_risk_virulence_genes']) + len(analysis['critical_virulence_genes']),
                    'Carbapenemase_Status': analysis['carbapenemase_status'],
                    'ESBL_Status': analysis['esbl_status'],
                    'Colistin_Resistance': analysis['colistin_resistance']
                })
            
            df = pd.DataFrame(excel_data)
            excel_file = os.path.join(output_base, "acineto_abricate_master_summary.xlsx")
            df.to_excel(excel_file, index=False)
            master_report["excel_summary_file"] = excel_file
            self.logger.info("✓ Created Excel summary: %s", excel_file)
        except Exception as e:
            self.logger.warning("Could not create Excel summary: %s", e)
        
        # Write master JSON report
        master_file = os.path.join(output_base, "acineto_abricate_master_summary.json")
        with open(master_file, 'w') as f:
            json.dump(master_report, f, indent=4, default=str)
        
        self.logger.info("✓ Created master JSON summary: %s", master_file)
        return master_file
    
    def _create_database_summary_html(self, database: str, hits: List[Dict], output_base: str):
        """Create HTML summary report for a specific database across all genomes"""
        
        # JavaScript for rotating quotes
        quotes_js = f"""
        <script>
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            
            // Rotate every 10 seconds
            setInterval(rotateQuote, 10000);
            
            // Initial display
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
            }});
        </script>
        """
        
        # Count unique genomes
        unique_genomes = list(set(hit['genome'] for hit in hits))
        
        # Count genes per genome
        genes_per_genome = {}
        for hit in hits:
            genome = hit['genome']
            if genome not in genes_per_genome:
                genes_per_genome[genome] = set()
            genes_per_genome[genome].add(hit['gene'])
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope ABRicate - {database.upper()} Database Summary</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
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
            max-width: 1400px;
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
        
        .card {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        
        .gene-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }}
        
        .gene-table th, .gene-table td {{
            padding: 15px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            font-weight: 600;
        }}
        
        tr:hover {{ background-color: #f8f9fa; }}
        
        .summary-stats {{
            display: flex;
            justify-content: space-around;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #fff;
        }}
        
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
        }}
        
        .footer a {{
            color: #667eea;
            text-decoration: none;
        }}
        
        .footer a:hover {{
            text-decoration: underline;
        }}
        
        .present {{ background-color: #d4edda; }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 AcinetoScope ABRicate - {database.upper()} Database Summary</h1>
                <p style="color: #666; font-size: 1.2em;">Cross-genome analysis of {database.upper()} database results</p>
            </div>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📊 Database Overview</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Hits</h3>
                    <p style="font-size: 2em; margin: 0;">{len(hits)}</p>
                </div>
                <div class="stat-card">
                    <h3>Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{len(unique_genomes)}</p>
                </div>
                <div class="stat-card">
                    <h3>Unique Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{len(set(hit['gene'] for hit in hits))}</p>
                </div>
            </div>
            <p><strong>Database:</strong> {database.upper()}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🔍 Genes by Genome</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Genome</th>
                        <th>Gene Count</th>
                        <th>Genes Detected</th>
                    </tr>
                </thead>
                <tbody>
"""
        
        for genome in sorted(unique_genomes):
            genes = genes_per_genome.get(genome, set())
            gene_list = ", ".join(sorted(genes))
            html_content += f"""
                    <tr class="present">
                        <td><strong>{genome}</strong></td>
                        <td>{len(genes)}</td>
                        <td>{gene_list}</td>
                    </tr>
"""
        
        html_content += """
                </tbody>
            </table>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📈 Gene Frequency</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Frequency</th>
                        <th>Genomes</th>
                    </tr>
                </thead>
                <tbody>
"""
        
        # Calculate gene frequency
        gene_frequency = {}
        for hit in hits:
            gene = hit['gene']
            if gene not in gene_frequency:
                gene_frequency[gene] = set()
            gene_frequency[gene].add(hit['genome'])
        
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            genome_list = ", ".join(sorted(genomes))
            html_content += f"""
                    <tr>
                        <td><strong>{gene}</strong></td>
                        <td>{len(genomes)}</td>
                        <td>{genome_list}</td>
                    </tr>
"""
        
        html_content += """
                </tbody>
            </table>
        </div>
        
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using AcinetoScope (ABRicate v1.2.0)
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write database summary HTML report
        html_file = os.path.join(output_base, f"acineto_{database}_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info("Database summary HTML report: %s", html_file)
    
    def process_single_genome(self, genome_file: str, output_base: str = "acineto_abricate_results") -> Dict[str, Any]:
        """Process a single A. baumannii genome with all databases and HTML reporting"""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        
        self.logger.info("=== PROCESSING A. BAUMANNII GENOME: %s ===", genome_name)
        
        # Create output directory
        os.makedirs(results_dir, exist_ok=True)
        
        databases = self.required_databases
        
        # Run ABRicate on all databases
        results = {}
        for db in databases:
            result = self.run_abricate_single_db(genome_file, db, results_dir)
            results[db] = result
            status_icon = "✓" if result['status'] == 'success' else "✗"
            self.logger.info("%s %s: %d hits", status_icon, db, result['hit_count'])
        
        # Create comprehensive HTML report
        self.create_comprehensive_html_report(genome_name, results, results_dir)
        
        return {
            'genome': genome_name,
            'results': results,
            'total_hits': sum(r['hit_count'] for r in results.values())
        }
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "acineto_abricate_results") -> Dict[str, Any]:
        """Process multiple A. baumannii genomes using wildcard pattern - MAXIMUM SPEED"""
        
        # Check ABRicate installation
        if not self.check_abricate_installed():
            raise RuntimeError("ABRicate not properly installed")
        
        # Setup databases
        self.setup_abricate_databases()
        
        # Find genome files (support all FASTA extensions)
        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa", 
                         f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        
        # Remove duplicates
        genome_files = list(set(genome_files))
        
        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")
        
        self.logger.info("Found %d A. baumannii genomes: %s", len(genome_files), [Path(f).name for f in genome_files])
        
        # Create output directory
        os.makedirs(output_base, exist_ok=True)
        
        # Process genomes with parallel execution - MAXIMUM SPEED CONFIGURATION
        all_results = {}
        
        if len(genome_files) > 1 and self.cpus > 1:
            # Use ThreadPoolExecutor for parallel processing of multiple genomes
            self.logger.info("Using parallel processing with %d CPU cores (MAXIMUM SPEED)", self.cpus)
            
            with ThreadPoolExecutor(max_workers=self.cpus) as executor:
                # Submit all genomes for processing
                future_to_genome = {
                    executor.submit(self.process_single_genome, genome, output_base): genome 
                    for genome in genome_files
                }
                
                # Collect results as they complete
                for future in as_completed(future_to_genome):
                    genome = future_to_genome[future]
                    try:
                        result = future.result()
                        all_results[Path(genome).stem] = result
                        self.logger.info("✓ Completed: %s (%d total hits)", result['genome'], result['total_hits'])
                    except Exception as e:
                        self.logger.error("✗ Failed: %s - %s", genome, e)
        else:
            # Process genomes sequentially
            for genome in genome_files:
                try:
                    result = self.process_single_genome(genome, output_base)
                    all_results[Path(genome).stem] = result
                    self.logger.info("✓ Completed: %s (%d total hits)", result['genome'], result['total_hits'])
                except Exception as e:
                    self.logger.error("✗ Failed: %s - %s", genome, e)
        
        # Create database summary files and HTML reports after processing all genomes
        self.create_database_summaries(all_results, output_base)
        
        # Create master JSON summary with Excel export
        master_json_file = self.create_master_json_summary(all_results, output_base)
        
        self.logger.info("=== A. BAUMANNII ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        self.logger.info("Master JSON summary: %s", master_json_file)
        
        return all_results


def main():
    """Command line interface for A. baumannii ABRicate analysis"""
    parser = argparse.ArgumentParser(
        description='AcinetoScope ABRicate Analysis - MAXIMUM SPEED VERSION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all A. baumannii FASTA files (auto-detect optimal CPU cores - MAXIMUM SPEED)
  python acineto_abricate.py "*.fna"
  
  # Run on specific pattern with auto CPU detection
  python acineto_abricate.py "AB_*.fasta"
  
  # Run with custom output directory and auto CPUs
  python acineto_abricate.py "*.fa" --output my_acineto_results

  # Force specific number of CPU cores
  python acineto_abricate.py "*.fna" --cpus 4

MAXIMUM SPEED RESOURCE MANAGEMENT:
  • 1-4 cores: Uses ALL CPU cores (100% utilization)
  • 5-8 cores: Uses (cores-1) for optimal performance  
  • 9-16 cores: Uses (cores-2) for high performance
  • 17-32 cores: Uses (cores-4) for maximum throughput
  • 32+ cores: Uses 85% of cores (capped at 32)

Supported FASTA extensions: .fasta, .fa, .fna, .faa
        """
    )
    
    parser.add_argument('pattern', help='File pattern for A. baumannii genomes (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--cpus', '-c', type=int, default=None, 
                       help='Number of CPU cores to use (default: auto-detect optimal for MAXIMUM SPEED)')
    parser.add_argument('--output', '-o', default='acineto_abricate_results', 
                       help='Output directory (default: acineto_abricate_results)')
    
    args = parser.parse_args()
    
    # Print AcinetoScope banner
    print("\n" + "="*80)
    print("""
 █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
""")
    print("A. baumannii ABRicate Analysis - MAXIMUM SPEED VERSION")
    print("="*80)
    print(f"Author: Brown Beckley | Email: brownbeckley94@gmail.com")
    print(f"Affiliation: University of Ghana Medical School - Department of Medical Biochemistry")
    print("="*80)
    
    executor = AcinetoAbricateExecutor(cpus=args.cpus)
    
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        
        # Print summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🧬 AcinetoScope ABRicate FINAL SUMMARY")
        executor.logger.info("="*50)
        
        total_critical_carbapenemase = 0
        total_critical_esbl = 0
        total_critical_colistin = 0
        
        for genome_name, result in results.items():
            # Collect all hits for this genome
            all_genome_hits = []
            for db_result in result['results'].values():
                all_genome_hits.extend(db_result['hits'])
            
            # Analyze for this genome
            analysis = executor.analyze_acineto_resistance(all_genome_hits)
            
            executor.logger.info("✓ %s: %d total hits, %d critical carbapenemase, %d critical ESBL, %d critical colistin", 
                               genome_name, result['total_hits'], 
                               len(analysis['critical_carbapenemase_genes']),
                               len(analysis['critical_esbl_genes']),
                               len(analysis['critical_colistin_genes']))
            
            total_critical_carbapenemase += len(analysis['critical_carbapenemase_genes'])
            total_critical_esbl += len(analysis['critical_esbl_genes'])
            total_critical_colistin += len(analysis['critical_colistin_genes'])
        
        # Database usage summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🗃️  DATABASE USAGE SUMMARY")
        executor.logger.info("="*50)
        executor.logger.info("Used databases: %s", ", ".join(executor.required_databases))
        
        # Summary files
        executor.logger.info("\n" + "="*50)
        executor.logger.info("📊 SUMMARY FILES")
        executor.logger.info("="*50)
        executor.logger.info("Created per-database summaries:")
        for db in executor.required_databases:
            executor.logger.info("  • acineto_%s_summary.[tsv|json|html]", db)
        executor.logger.info("Created master summaries:")
        executor.logger.info("  • acineto_abricate_master_summary.json")
        executor.logger.info("  • acineto_abricate_master_summary.xlsx (Excel)")
        
        # Performance summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("⚡ MAXIMUM SPEED PERFORMANCE SUMMARY")
        executor.logger.info("="*50)
        executor.logger.info("CPU cores utilized: %d cores", executor.cpus)
        executor.logger.info("Available RAM: %.1f GB", executor.available_ram)
        executor.logger.info("Total A. baumannii genomes processed: %d", len(results))
        executor.logger.info("Total critical carbapenemase genes found: %d", total_critical_carbapenemase)
        executor.logger.info("Total critical ESBL genes found: %d", total_critical_esbl)
        executor.logger.info("Total critical colistin resistance genes found: %d", total_critical_colistin)
        executor.logger.info("Processing mode: MAXIMUM SPEED 🚀")
        
        import random
        executor.logger.info("\n💡 %s", random.choice(executor.science_quotes))
        
        # Output directory structure
        executor.logger.info("\n" + "="*50)
        executor.logger.info("📁 OUTPUT STRUCTURE")
        executor.logger.info("="*50)
        executor.logger.info("%s/", args.output)
        executor.logger.info("├── [SAMPLE_DIRECTORIES]/")
        executor.logger.info("│   ├── abricate_*.txt           # Raw ABRicate outputs")
        executor.logger.info("│   ├── abricate_*_report.html   # Individual database reports")
        executor.logger.info("│   └── [SAMPLE]_comprehensive_abricate_report.html")
        executor.logger.info("├── acineto_*_abricate_summary.tsv   # Per-database TSV summaries")
        executor.logger.info("├── acineto_*_summary.json           # Per-database JSON summaries")
        executor.logger.info("├── acineto_*_summary_report.html    # Per-database HTML summaries")
        executor.logger.info("├── acineto_abricate_master_summary.json")
        executor.logger.info("└── acineto_abricate_master_summary.xlsx")
        
    except Exception as e:
        executor.logger.error("A. baumannii analysis failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()