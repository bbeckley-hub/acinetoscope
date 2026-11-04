#!/usr/bin/env python3
"""
AcinetoScope ABRicate Module - MAXIMUM SPEED VERSION with configurable thresholds
Comprehensive ABRicate analysis for Acinetobacter baumannii with HTML reporting
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2025-12-25 (updated 2026-11-01)

Outputs per genome (2 files):
  - {genome}_comprehensive_abricate_report.html
  - abricate_raw.tsv                    (all databases merged)

Batch-level outputs:
  - acineto_{db}_abricate_summary.tsv   (per database)
  - acineto_{db}_summary_report.html    (per database)
  - acineto_abricate_master_summary.json  (single batch JSON)
  - acineto_abricate_master_summary.xlsx  (optional Excel summary)

Send a quick mail for any issues or further explanations.
"""

import subprocess
import sys
import os
import glob
import logging
import tempfile
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

    def __init__(self, cpus: int = None, min_identity: float = 80, min_coverage: float = 80):
        self.logger = self._setup_logging()
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)
        self.min_identity = min_identity
        self.min_coverage = min_coverage

        self.required_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoh', 'bacmet2', 'ecoli_vf'
        ]

        # CRITICAL CARBAPENEMASE genes (highest priority)
        self.critical_carbapenemases = {
            'blaOXA-23', 'blaOXA-24', 'blaOXA-40', 'blaOXA-51', 'blaOXA-58', 'blaOXA-66', 'blaOXA-72',
            'blaOXA-143', 'blaOXA-235', 'blaOXA-181', 'blaOXA-232', 'blaOXA-438', 'blaOXA-48',
            'blaKPC', 'blaNDM', 'blaVIM', 'blaIMP', 'blaSIM', 'blaGIM', 'blaSPM',
            'OXA-23', 'OXA-24', 'OXA-40', 'OXA-51', 'OXA-58', 'OXA-66', 'OXA-72', 'OXA-143', 'OXA-235',
            'KPC', 'NDM', 'VIM', 'IMP', 'SIM', 'GIM', 'SPM'
        }

        self.critical_esbls = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'blaPER', 'blaVEB', 'blaGES', 'blaBEL', 'blaBES',
            'CTX-M', 'SHV', 'TEM', 'PER', 'VEB', 'GES', 'BEL', 'BES',
            'blaADC', 'blaOXA-10', 'blaOXA-2', 'blaOXA-9',
            'ADC', 'OXA-10', 'OXA-2', 'OXA-9'
        }

        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'pmrA', 'pmrB', 'pmrC', 'lpxA', 'lpxC', 'lpxD', 'arnA', 'arnB', 'arnC', 'arnD', 'eptA'
        }

        self.critical_aminoglycoside = {
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH',
            'APH(3\')', 'APH(3\')-VI', 'APH(6)', 'AAC(3)', 'AAC(6\')', 'ANT(2")', 'ANT(3")', 'ANT(4")',
            'aacC1', 'aacC2', 'aacC4', 'aacA4', 'aacA7', 'aadA1', 'aadA2', 'aadA5', 'aadA7',
            'strA', 'strB', 'aphA1', 'aphA2', 'aphA3', 'aphA6', 'aac3', 'aac6', 'aadA', 'aadB'
        }

        self.high_risk_resistance = {
            'tetA', 'tetB', 'tetC', 'tetD', 'tetE', 'tetG', 'tetH', 'tetK', 'tetL', 'tetM', 'tetO', 'tetQ', 'tetS', 'tetX',
            'sul1', 'sul2', 'sul3', 'sul4',
            'dfrA1', 'dfrA5', 'dfrA7', 'dfrA8', 'dfrA12', 'dfrA14', 'dfrA17', 'dfrA19', 'dfrA20', 'dfrA21',
            'dfrB1', 'dfrB2', 'dfrB3', 'dfrB4', 'dfrB5', 'dfrB6', 'dfrB7',
            'catA1', 'catA2', 'catB2', 'catB3', 'catB8', 'catI', 'catII', 'catIII',
            'cmlA', 'cmlA1', 'cmlA5', 'cmlA6', 'cmlA7', 'floR',
            'ermA', 'ermB', 'ermC', 'ermF', 'ermG', 'ermX', 'ermY',
            'mphA', 'mphB', 'mphC', 'mphD', 'mphE', 'msrA', 'msrB', 'msrC', 'msrD',
            'qnrA1', 'qnrB1', 'qnrB2', 'qnrB4', 'qnrB6', 'qnrB10', 'qnrB19', 'qnrS1', 'qnrS2', 'qnrVC1', 'qnrVC4',
            'aac(6\')-Ib-cr', 'qepA', 'qepA1', 'qepA2', 'qepA3', 'qepA4',
            'fosA', 'fosB', 'fosC', 'fosX',
            'arr-2', 'arr-3', 'arr-4', 'arr-5', 'arr-6', 'arr-7',
            'adeA', 'adeB', 'adeC', 'adeG', 'adeH', 'adeI', 'adeJ', 'adeK', 'adeL', 'adeM', 'adeN',
            'abeM', 'abeS', 'acrA', 'acrB', 'tolC', 'mexA', 'mexB', 'mexC', 'mexD', 'mexE', 'mexF',
            'oprM', 'oprN', 'oprJ', 'mdtA', 'mdtB', 'mdtC', 'emrA', 'emrB', 'emrD', 'emrE', 'emrK', 'emrY'
        }

        self.virulence_genes = {
            'csuA', 'csuB', 'csuC', 'csuD', 'csuE', 'csuA/B', 'csuC/D/E',
            'bfmR', 'bfmS', 'ompA', 'bap', 'blaPER-1', 'epsA', 'ptk',
            'basA', 'basB', 'basC', 'basD', 'basE', 'basF', 'basG', 'basH', 'basI', 'basJ',
            'barA', 'barB', 'bauA', 'bauB', 'bauC', 'bauD', 'bauE', 'bauF',
            'entA', 'entB', 'entC', 'entD', 'entE', 'entF', 'entS',
            'fhuA', 'fhuB', 'fhuC', 'fhuD', 'fhuE', 'fhuF',
            'fepA', 'fepB', 'fepC', 'fepD', 'fepE', 'fepG',
            'abaI', 'abaR', 'luxI', 'luxR',
            'hlyA', 'hlyB', 'hlyC', 'hlyD', 'rtxA', 'rtxB', 'rtxC', 'rtxD', 'rtxE',
            'tssA', 'tssB', 'tssC', 'tssD', 'tssE', 'tssF', 'tssG', 'tssH', 'tssI', 'tssJ', 'tssK', 'tssL', 'tssM',
            'vgrG', 'vgrG1', 'vgrG2', 'vgrG3', 'vgrG4', 'vgrG5',
            'plc1', 'plc2', 'plcD', 'lipA', 'lipB', 'lipC', 'lipH',
            'aprA', 'aprB', 'aprC', 'aprD', 'aprE', 'aprF',
            'cpsA', 'cpsB', 'cpsC', 'cpsD', 'cpsE', 'cpsF', 'cpsG', 'cpsH', 'cpsI', 'cpsJ',
            'kpsM', 'kpsS', 'kpsT', 'kpsU', 'kpsV', 'kpsW',
            'wza', 'wzb', 'wzc', 'wzd', 'wze', 'wzf', 'wzg',
            'fimA', 'fimB', 'fimC', 'fimD', 'fimE', 'fimF', 'fimG', 'fimH',
            'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH',
            'sfmA', 'sfmB', 'sfmC', 'sfmD', 'sfmE', 'sfmF', 'sfmG', 'sfmH',
            'pilA', 'pilB', 'pilC', 'pilD', 'pilE', 'pilF', 'pilG', 'pilH', 'pilI', 'pilJ',
            'pilK', 'pilL', 'pilM', 'pilN', 'pilO', 'pilP', 'pilQ'
        }

        self.critical_resistance_genes = (
            self.critical_carbapenemases
            .union(self.critical_esbls)
            .union(self.critical_colistin)
            .union(self.critical_aminoglycoside)
        )

        self.critical_virulence_genes = {
            'hlyA', 'hlyB', 'hlyC', 'hlyD', 'rtxA', 'rtxB', 'rtxC', 'rtxD', 'rtxE'
        }

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
        self.ascii_art = r""" █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝""".strip("\n")

        self.metadata = {
            "tool_name": "AcinetoScope ABRicate",
            "version": "1.4.1",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "min_identity": self.min_identity,
            "min_coverage": self.min_coverage,
        }

    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)

    def _get_available_ram(self) -> float:
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except Exception as e:
            self.logger.warning(f"Could not detect RAM: {e}")
            return 8

    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        if user_cpus is not None:
            self._log_resource_info(user_cpus)
            return user_cpus
        try:
            total_physical_cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            if total_physical_cores <= 4:
                optimal_cpus = total_physical_cores
            elif total_physical_cores <= 8:
                optimal_cpus = total_physical_cores - 1
            elif total_physical_cores <= 16:
                optimal_cpus = max(8, total_physical_cores - 1)
            elif total_physical_cores <= 32:
                optimal_cpus = max(16, total_physical_cores - 2)
            else:
                optimal_cpus = min(32, int(total_physical_cores * 0.95))
            optimal_cpus = max(1, min(optimal_cpus, total_physical_cores))
            self._log_resource_info(optimal_cpus, total_physical_cores)
            return optimal_cpus
        except Exception as e:
            self.logger.warning(f"Could not detect CPU cores, using maximum available: {e}")
            return os.cpu_count() or 4

    def _log_resource_info(self, cpus: int, total_cores: int = None):
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            utilization = (cpus / total_cores) * 100
            self.logger.info(f"Using CPU cores: {cpus} ({utilization:.1f}% of available cores)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        if cpus == 1:
            self.logger.info("Performance: single-core")
        elif cpus <= 4:
            self.logger.info("Performance: multi-core (small system)")
        elif cpus <= 8:
            self.logger.info("Performance: high-speed mode")
        else:
            self.logger.info("Performance: MAXIMUM SPEED MODE")

    def check_abricate_installed(self) -> bool:
        try:
            result = subprocess.run(['abricate', '--version'],
                                    capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            self.logger.info("ABRicate version: %s", version_line)
            version_match = re.search(r'(\d+\.\d+\.\d+)', version_line)
            if version_match:
                version_str = version_match.group(1)
                if version_str >= "1.0.1":
                    self.logger.info("ABRicate version meets requirement (>=1.2.0)")
                    return True
                self.logger.error("ABRicate version too old: %s. Required >=1.2.0", version_str)
                return False
            self.logger.info("ABRicate installed (version check skipped)")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("ABRicate not found. Install with: conda install -c bioconda abricate")
            return False

    def setup_abricate_databases(self):
        self.logger.info("Setting up ABRicate databases for A. baumannii analysis...")
        available_dbs = []
        missing_dbs = []
        try:
            check_result = subprocess.run(['abricate', '--list'],
                                          capture_output=True, text=True, check=True)
            for db in self.required_databases:
                if db in check_result.stdout:
                    self.logger.info("Database available: %s", db)
                    available_dbs.append(db)
                else:
                    self.logger.warning("Database not available: %s", db)
                    missing_dbs.append(db)
            for db in missing_dbs:
                self.logger.info("Attempting to setup database: %s", db)
                try:
                    subprocess.run(['abricate', '--setupdb', '--db', db],
                                   capture_output=True, text=True, check=True)
                    self.logger.info("Database setup completed: %s", db)
                    available_dbs.append(db)
                except subprocess.CalledProcessError as e:
                    self.logger.error("Failed to setup database %s: %s", db, e.stderr)
            self.required_databases = available_dbs
            self.logger.info("Using databases: %s", ", ".join(self.required_databases))
        except subprocess.CalledProcessError as e:
            self.logger.error("Error checking ABRicate databases: %s", e.stderr)
        except Exception as e:
            self.logger.error("Unexpected error setting up databases: %s", e)

    def run_abricate_single_db(self, genome_file: str, database: str, temp_dir: str,
                               min_identity: float = None, min_coverage: float = None) -> Dict[str, Any]:
        """Run ABRicate on one genome/one DB; write to temp_dir; return hits + raw content."""
        genome_name = Path(genome_file).stem
        raw_file = os.path.join(temp_dir, f"abricate_{database}.txt")
        min_id = min_identity if min_identity is not None else self.min_identity
        min_cov = min_coverage if min_coverage is not None else self.min_coverage

        cmd = ['abricate', genome_file, '--db', database,
               '--minid', str(min_id), '--mincov', str(min_cov)]
        self.logger.info("Running ABRicate: %s --db %s (minid=%.1f, mincov=%.1f)",
                         genome_name, database, min_id, min_cov)
        try:
            with open(raw_file, 'w') as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE,
                               text=True, check=True)
            hits = self._parse_abricate_output(raw_file)
            for hit in hits:
                if not hit.get('database'):
                    hit['database'] = database
            with open(raw_file, 'r') as f:
                raw_content = f.read()
            return {
                'database': database,
                'genome': genome_name,
                'hits': hits,
                'hit_count': len(hits),
                'raw_content': raw_content,
                'status': 'success',
            }
        except subprocess.CalledProcessError as e:
            self.logger.error("ABRicate failed for %s on %s: %s", database, genome_name, e.stderr)
            return {
                'database': database,
                'genome': genome_name,
                'hits': [],
                'hit_count': 0,
                'raw_content': '',
                'status': 'failed',
            }

    def _parse_abricate_output(self, abricate_file: str) -> List[Dict]:
        """Parse ABRicate output file (robust to tabs inside fields)."""
        hits = []
        try:
            with open(abricate_file, 'r') as f:
                lines = f.readlines()
            if not lines:
                return hits
            headers = []
            data_lines = []
            for line in lines:
                if line.startswith('#FILE') and not headers:
                    headers = line.strip().replace('#', '').split('\t')
                elif line.strip() and not line.startswith('#'):
                    data_lines.append(line.strip())
            if not headers:
                self.logger.warning("No headers found in %s", abricate_file)
                return hits
            expected_columns = len(headers)
            for line_num, line in enumerate(data_lines, 1):
                parts = line.split('\t')
                if len(parts) > expected_columns:
                    combined = parts[:expected_columns - 1]
                    combined.append('\t'.join(parts[expected_columns - 1:]))
                    parts = combined
                elif len(parts) < expected_columns:
                    parts.extend([''] * (expected_columns - len(parts)))
                if len(parts) == expected_columns:
                    hit = dict(zip(headers, parts))
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
                        'resistance': hit.get('RESISTANCE', ''),
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning(
                        "Line %d has %d parts, expected %d: %s",
                        line_num, len(parts), expected_columns, line[:200] + "..."
                    )
        except Exception as e:
            self.logger.error("Error parsing %s: %s", abricate_file, e)
        self.logger.info("Parsed %d hits from %s", len(hits), abricate_file)
        return hits

    def _write_merged_raw_tsv(self, results: Dict[str, Dict], output_dir: str):
        """Merge all per-database ABRicate outputs into a single abricate_raw.tsv."""
        output_file = os.path.join(output_dir, "abricate_raw.tsv")
        header_written = False
        data_lines = 0
        with open(output_file, 'w') as f:
            for db, result in results.items():
                raw = result.get('raw_content', '') or ''
                if not raw:
                    continue
                for line in raw.splitlines():
                    if not line.strip():
                        continue
                    if line.startswith('#'):
                        if not header_written:
                            f.write(line + '\n')
                            header_written = True
                    else:
                        f.write(line + '\n')
                        data_lines += 1
        self.logger.info("Merged raw TSV: %s (%d hits across %d databases)",
                         output_file, data_lines, len(results))

    # ------------------------------------------------------------------
    # Shared footer (used by every HTML report for consistency)
    # ------------------------------------------------------------------
    def _get_footer_html(self, title: str) -> str:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return f"""
        <div class="footer">
            <p><strong>ACINETOSCOPE</strong> &mdash; ABRICATE Analysis Module</p>
            <p class="timestamp">Generated: <span id="acinetoscope-timestamp"></span></p>
            <div class="authorship">
                <p><strong>Technical Support &amp; Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
            <script>document.addEventListener("DOMContentLoaded", () => document.getElementById("acinetoscope-timestamp").textContent = new Date().toLocaleString());</script>
        </div>
        """

    # ------------------------------------------------------------------
    # Shared ASCII styling (used by every HTML report)
    # ------------------------------------------------------------------
    def _get_ascii_css(self) -> str:
        """CSS block for reliably rendering box-drawing ASCII art in HTML."""
        return """
        .ascii-container {
            background: rgba(0, 0, 0, 0.9);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
            overflow-x: auto;
        }
        pre.ascii-art {
            font-family: 'Consolas', 'Courier New', 'DejaVu Sans Mono',
                         'Liberation Mono', 'Menlo', 'Monaco', monospace;
            font-size: 12px;
            line-height: 1.05;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            font-variant-ligatures: none;
            font-kerning: none;
            letter-spacing: 0;
            word-spacing: 0;
            margin: 0;
            padding: 0;
            display: inline-block;
            min-width: max-content;
            tab-size: 4;
        }
        @media (max-width: 768px) {
            pre.ascii-art { font-size: 6px; }
        }
        """

    # ------------------------------------------------------------------
    # Resistance analysis
    # ------------------------------------------------------------------
    def analyze_acineto_resistance(self, all_hits: List[Dict]) -> Dict[str, Any]:
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

            if any(carba in gene_lower for carba in ['oxa', 'kpc', 'ndm', 'vim', 'imp', 'sim']):
                if any(crit in gene_lower for crit in [g.lower() for g in self.critical_carbapenemases]):
                    analysis['carbapenemase_status'] = 'positive'
                    analysis['critical_carbapenemase_genes'].append({
                        'gene': gene, 'product': hit['product'], 'database': hit['database'],
                        'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                    })
            elif any(esbl in gene_lower for esbl in ['ctx', 'shv', 'tem', 'per', 'veb', 'ges']):
                if any(crit in gene_lower for crit in [g.lower() for g in self.critical_esbls]):
                    analysis['esbl_status'] = 'positive'
                    analysis['critical_esbl_genes'].append({
                        'gene': gene, 'product': hit['product'], 'database': hit['database'],
                        'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                    })
            elif any(mcr in gene_lower for mcr in ['mcr']):
                if any(crit in gene_lower for crit in [g.lower() for g in self.critical_colistin]):
                    analysis['colistin_resistance'] = 'positive'
                    analysis['critical_colistin_genes'].append({
                        'gene': gene, 'product': hit['product'], 'database': hit['database'],
                        'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                    })
            elif any(ag in gene_lower for ag in ['arm', 'rmt', 'aph', 'aac', 'aad', 'str', 'ant']):
                if any(crit in gene_lower for crit in [g.lower() for g in self.critical_aminoglycoside]):
                    analysis['critical_aminoglycoside_genes'].append({
                        'gene': gene, 'product': hit['product'], 'database': hit['database'],
                        'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                    })
            elif any(hr in gene_lower for hr in [g.lower() for g in self.high_risk_resistance]):
                analysis['high_risk_resistance_genes'].append({
                    'gene': gene, 'product': hit['product'], 'database': hit['database'],
                    'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                })

            if any(crit_vf in gene_base for crit_vf in self.critical_virulence_genes):
                analysis['critical_virulence_genes'].append({
                    'gene': gene, 'product': hit['product'], 'database': hit['database'],
                    'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                })
            elif any(vf in gene_lower for vf in [g.lower() for g in self.virulence_genes]):
                is_high = any(hr in gene_lower for hr in ['hly', 'rtx', 'tss', 'vgr', 'plc', 'lip', 'apr'])
                if is_high:
                    analysis['high_risk_virulence_genes'].append({
                        'gene': gene, 'product': hit['product'], 'database': hit['database'],
                        'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                    })
                else:
                    analysis['moderate_risk_genes'].append({
                        'gene': gene, 'product': hit['product'], 'database': hit['database'],
                        'coverage': hit['coverage_percent'], 'identity': hit['identity_percent'],
                    })

            resistance_class = self._classify_resistance(hit['product'], gene)
            if resistance_class:
                if resistance_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][resistance_class] = []
                if gene not in [g['gene'] for g in analysis['resistance_classes'][resistance_class]]:
                    analysis['resistance_classes'][resistance_class].append({
                        'gene': gene, 'product': hit['product'],
                    })

        analysis['total_critical_resistance'] = (
            len(analysis['critical_carbapenemase_genes'])
            + len(analysis['critical_esbl_genes'])
            + len(analysis['critical_colistin_genes'])
            + len(analysis['critical_aminoglycoside_genes'])
        )
        analysis['total_high_risk_resistance'] = len(analysis['high_risk_resistance_genes'])
        analysis['total_critical_virulence'] = len(analysis['critical_virulence_genes'])
        analysis['total_high_risk_virulence'] = len(analysis['high_risk_virulence_genes'])
        return analysis

    def _classify_resistance(self, product: str, gene: str) -> str:
        product_lower = product.lower()
        gene_lower = gene.lower()

        def matches(*terms):
            return any(t in product_lower or t in gene_lower for t in terms)

        if matches('carbapenem', 'oxa', 'kpc', 'ndm', 'vim', 'imp', 'sim'):
            return 'Carbapenem resistance'
        elif matches('beta-lactam', 'esbl', 'ctx', 'shv', 'tem', 'per', 'veb', 'ges'):
            return 'Beta-lactam resistance'
        elif matches('aminoglycoside', 'aac', 'aad', 'aph', 'str', 'arm', 'rmt'):
            return 'Aminoglycoside resistance'
        elif matches('tetracycline', 'tet'):
            return 'Tetracycline resistance'
        elif matches('sulfonamide', 'sul'):
            return 'Sulfonamide resistance'
        elif matches('trimethoprim', 'dfr'):
            return 'Trimethoprim resistance'
        elif matches('chloramphenicol', 'cat', 'cml', 'flor'):
            return 'Chloramphenicol resistance'
        elif matches('macrolide', 'erm', 'mph', 'msr'):
            return 'Macrolide resistance'
        elif matches('quinolone', 'qnr', 'qep'):
            return 'Quinolone resistance'
        elif matches('colistin', 'polymyxin', 'mcr'):
            return 'Polymyxin resistance'
        elif matches('fosfomycin', 'fos'):
            return 'Fosfomycin resistance'
        elif matches('rifampicin', 'arr'):
            return 'Rifampicin resistance'
        elif matches('efflux', 'ade', 'abe', 'mex', 'acr', 'emr'):
            return 'Efflux pumps'
        elif matches('virulence', 'toxin', 'hemolysin', 'siderophore'):
            return 'Virulence factors'
        return 'Other resistance'

    # ------------------------------------------------------------------
    # Per-genome comprehensive report
    # ------------------------------------------------------------------
    def create_comprehensive_html_report(self, genome_name: str, results: Dict, output_dir: str):
        all_hits = []
        for db_result in results.values():
            all_hits.extend(db_result['hits'])
        analysis = self.analyze_acineto_resistance(all_hits)

        quotes_js = f"""
        <script>
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            setInterval(rotateQuote, 10000);
            document.addEventListener('DOMContentLoaded', rotateQuote);
        </script>
        """

        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope ABRicate Analysis Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff; padding: 20px; min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        {self._get_ascii_css()}
        .card {{
            background: rgba(255, 255, 255, 0.95); color: #1f2937;
            padding: 25px; margin: 20px 0; border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        .gene-table {{
            width: 100%; border-collapse: collapse; margin: 20px 0;
            background: white; border-radius: 8px; overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }}
        .gene-table th, .gene-table td {{
            padding: 15px; text-align: left; border-bottom: 1px solid #e0e0e0;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white; font-weight: 600;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .critical {{ background-color: #f8d7da; font-weight: bold; }}
        .high-risk {{ background-color: #fff3cd; }}
        .summary-stats {{
            display: flex; justify-content: space-around; margin: 20px 0; flex-wrap: wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white; padding: 20px; border-radius: 12px; text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2); margin: 10px; flex: 1; min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1); color: white; padding: 20px;
            border-radius: 12px; margin: 20px 0; text-align: center;
            font-style: italic; border-left: 4px solid #fff;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8); color: white; padding: 30px;
            border-radius: 12px; margin-top: 40px;
        }}
        .footer a {{ color: #667eea; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .risk-badge {{
            display: inline-block; background: #dc3545; color: white;
            padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em;
        }}
        .warning-badge {{
            display: inline-block; background: #ffc107; color: black;
            padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em;
        }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <pre class="ascii-art">{self.ascii_art}</pre>
            </div>
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">AcinetoScope ABRicate Analysis Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive A. baumannii Antimicrobial Resistance &amp; Virulence Analysis</p>
            </div>
        </div>

        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">A. baumannii AMR Summary</h2>
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
            <p><strong>Min Identity:</strong> {self.metadata['min_identity']}%</p>
            <p><strong>Min Coverage:</strong> {self.metadata['min_coverage']}%</p>
        </div>
"""

        critical_alerts = []
        if analysis['carbapenemase_status'] == 'positive':
            critical_alerts.append("CARBAPENEMASE DETECTED")
        if analysis['esbl_status'] == 'positive':
            critical_alerts.append("ESBL DETECTED")
        if analysis['colistin_resistance'] == 'positive':
            critical_alerts.append("COLISTIN RESISTANCE DETECTED")
        if critical_alerts:
            html_content += """
        <div class="card" style="border-left: 4px solid #dc3545;">
            <h2 style="color: #dc3545;">CRITICAL RESISTANCE ALERTS</h2>
            <div style="margin: 10px 0;">
"""
            for alert in critical_alerts:
                html_content += f'<span class="risk-badge">{alert}</span>'
            html_content += "</div></div>"

        if analysis['resistance_classes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Resistance Classes Detected</h2>
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

        # ---- Render risk group table WITHOUT Risk Level column ----
        def render_group(title, genes, row_class):
            if not genes:
                return ""
            block = f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">{title}</h2>
            <table class="gene-table">
                <thead><tr>
                    <th>Gene</th><th>Product</th><th>Database</th>
                    <th>Coverage</th><th>Identity</th>
                </tr></thead><tbody>
"""
            for g in genes:
                prod = g['product']
                disp = prod[:97] + "..." if len(prod) > 100 else prod
                block += f"""
                    <tr class="{row_class}">
                        <td><strong>{g['gene']}</strong></td>
                        <td title="{prod}">{disp}</td>
                        <td>{g['database']}</td>
                        <td>{g['coverage']}%</td>
                        <td>{g['identity']}%</td>
                    </tr>
"""
            block += "</tbody></table></div>"
            return block

        html_content += render_group("CRITICAL Carbapenemase Genes",
                                     analysis['critical_carbapenemase_genes'], "critical")
        html_content += render_group("Critical ESBL Genes",
                                     analysis['critical_esbl_genes'], "critical")
        html_content += render_group("Colistin Resistance Genes",
                                     analysis['critical_colistin_genes'], "critical")
        html_content += render_group("Critical Aminoglycoside Resistance Genes",
                                     analysis['critical_aminoglycoside_genes'], "critical")
        html_content += render_group("High-Risk Resistance Genes",
                                     analysis['high_risk_resistance_genes'], "high-risk")
        html_content += render_group("Critical Virulence Genes",
                                     analysis['critical_virulence_genes'], "critical")
        html_content += render_group("High-Risk Virulence Genes",
                                     analysis['high_risk_virulence_genes'], "high-risk")

        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Database Results Summary</h2>
            <table class="gene-table">
                <thead><tr><th>Database</th><th>Hits</th><th>Status</th></tr></thead><tbody>
"""
        for db, result in results.items():
            status_icon = "OK" if result['status'] == 'success' else "FAIL"
            html_content += f"""
                    <tr>
                        <td>{db}</td>
                        <td>{result['hit_count']}</td>
                        <td>{status_icon} {result['status']}</td>
                    </tr>
"""
        html_content += "</tbody></table></div>"

        html_content += self._get_footer_html("AcinetoScope - Comprehensive ABRicate Report")
        html_content += "</div></body></html>"

        html_file = os.path.join(output_dir, f"{genome_name}_comprehensive_abricate_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Comprehensive HTML report: %s", html_file)

    # ------------------------------------------------------------------
    # Batch-level summaries
    # ------------------------------------------------------------------
    def create_database_summaries(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating A. baumannii database summary files and HTML reports...")
        db_results = {}
        for genome_name, genome_result in all_results.items():
            for db, db_result in genome_result['results'].items():
                db_results.setdefault(db, [])
                for hit in db_result['hits']:
                    hit_with_genome = hit.copy()
                    hit_with_genome['genome'] = genome_name
                    db_results[db].append(hit_with_genome)

        for db, hits in db_results.items():
            if hits:
                summary_file = os.path.join(output_base, f"acineto_{db}_abricate_summary.tsv")
                headers = list(hits[0].keys())
                with open(summary_file, 'w') as f:
                    f.write('\t'.join(headers) + '\n')
                    for hit in hits:
                        f.write('\t'.join(str(hit.get(h, '')) for h in headers) + '\n')
                self.logger.info("Created %s summary: %s (%d hits)",
                                 db, summary_file, len(hits))
                self._create_database_summary_html(db, hits, output_base)
            else:
                self.logger.info("No hits for database %s, skipping summary", db)

    def _create_database_summary_html(self, database: str, hits: List[Dict], output_base: str):
        quotes_js = f"""
        <script>
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            setInterval(rotateQuote, 10000);
            document.addEventListener('DOMContentLoaded', rotateQuote);
        </script>
        """

        unique_genomes = list(set(hit['genome'] for hit in hits))
        genes_per_genome = {}
        for hit in hits:
            genes_per_genome.setdefault(hit['genome'], set()).add(hit['gene'])

        genes_by_genome_rows = ""
        for genome in sorted(unique_genomes):
            genes = genes_per_genome.get(genome, set())
            genes_by_genome_rows += f"""
                    <tr><td><strong>{genome}</strong></td>
                        <td>{len(genes)}</td>
                        <td>{", ".join(sorted(genes))}</td></tr>
"""

        gene_frequency = {}
        for hit in hits:
            gene_frequency.setdefault(hit['gene'], set()).add(hit['genome'])
        frequency_rows = ""
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            frequency_rows += f"""
                    <tr><td><strong>{gene}</strong></td>
                        <td>{len(genomes)}</td>
                        <td>{", ".join(sorted(genomes))}</td></tr>
"""

        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope - {database.upper()} Batch Summary</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff; padding: 20px; min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        {self._get_ascii_css()}
        .card {{
            background: rgba(255, 255, 255, 0.95); color: #1f2937;
            padding: 25px; margin: 20px 0; border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        .gene-table {{
            width: 100%; border-collapse: collapse; margin: 20px 0;
            background: white; border-radius: 8px; overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }}
        .gene-table th, .gene-table td {{
            padding: 15px; text-align: left; border-bottom: 1px solid #e0e0e0;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white; font-weight: 600;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .summary-stats {{
            display: flex; justify-content: space-around; margin: 20px 0; flex-wrap: wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white; padding: 20px; border-radius: 12px; text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2); margin: 10px; flex: 1; min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1); color: white; padding: 20px;
            border-radius: 12px; margin: 20px 0; text-align: center;
            font-style: italic; border-left: 4px solid #fff;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8); color: white; padding: 30px;
            border-radius: 12px; margin-top: 40px;
        }}
        .footer a {{ color: #667eea; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .present {{ background-color: #d4edda; }}
    
        /* ════════ AcinetoScope uniformity override (auto-added) ════════ */
        .ascii-container {{ text-align: center !important; }}
        .ascii-art {{ display: inline-block !important; text-align: left !important; }}
        .footer {{ text-align: center !important; margin-top: 30px !important; padding: 20px !important; background: rgba(0, 0, 0, 0.3) !important; border-radius: 10px !important; font-size: 14px !important; overflow-wrap: break-word !important; }}
        .footer a {{ color: #fbbf24 !important; text-decoration: none !important; }}
        .footer a:hover {{ text-decoration: underline !important; }}
        .footer .timestamp {{ color: #fbbf24 !important; font-weight: bold !important; }}
        .footer .authorship {{ margin-top: 15px !important; padding: 15px !important; background: rgba(255, 255, 255, 0.1) !important; border-radius: 8px !important; font-size: 12px !important; line-height: 1.7 !important; }}
    </style>
    {quotes_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <pre class="ascii-art">{self.ascii_art}</pre>
            </div>
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">AcinetoScope ABRicate - {database.upper()} Database Summary</h1>
                <p style="color: #666; font-size: 1.2em;">Cross-genome analysis of {database.upper()} database results</p>
            </div>
        </div>

        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Database Overview</h2>
            <div class="summary-stats">
                <div class="stat-card"><h3>Total Hits</h3>
                    <p style="font-size: 2em; margin: 0;">{len(hits)}</p></div>
                <div class="stat-card"><h3>Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{len(unique_genomes)}</p></div>
                <div class="stat-card"><h3>Unique Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{len(set(hit['gene'] for hit in hits))}</p></div>
            </div>
            <p><strong>Database:</strong> {database.upper()}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Min Identity:</strong> {self.metadata['min_identity']}%</p>
            <p><strong>Min Coverage:</strong> {self.metadata['min_coverage']}%</p>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Genes by Genome</h2>
            <table class="gene-table">
                <thead><tr><th>Genome</th><th>Gene Count</th><th>Genes Detected</th></tr></thead>
                <tbody>{genes_by_genome_rows}</tbody>
            </table>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Gene Frequency</h2>
            <table class="gene-table">
                <thead><tr><th>Gene</th><th>Frequency</th><th>Genomes</th></tr></thead>
                <tbody>{frequency_rows}</tbody>
            </table>
        </div>

        {self._get_footer_html(f"AcinetoScope - {database.upper()} Batch Summary")}
    </div>
</body>
</html>"""

        html_file = os.path.join(output_base, f"acineto_{database}_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Database summary HTML report: %s", html_file)

    # ------------------------------------------------------------------
    # Single master JSON (+ optional Excel)
    # ------------------------------------------------------------------
    def create_master_json_summary(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating master JSON summary for all databases and genomes...")

        master_report = {
            "metadata": {
                **self.metadata,
                "analysis_type": "A. baumannii Comprehensive Antimicrobial Resistance & Virulence Analysis",
                "databases_analyzed": self.required_databases,
                "total_genomes": len(all_results),
                "analysis_parameters": {
                    "minimum_identity": self.metadata['min_identity'],
                    "minimum_coverage": self.metadata['min_coverage'],
                },
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
                "colistin_resistant_genomes": 0,
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
                "resistance_classes_found": defaultdict(list),
            },
            "detailed_results": {},
        }

        for genome_name, genome_result in all_results.items():
            all_genome_hits = []
            for db_result in genome_result['results'].values():
                all_genome_hits.extend(db_result['hits'])
            analysis = self.analyze_acineto_resistance(all_genome_hits)

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
                "database_hits": {db: r['hit_count'] for db, r in genome_result['results'].items()},
            }

            for key, bucket in [
                ("critical_carbapenemase_genes_found", analysis['critical_carbapenemase_genes']),
                ("critical_esbl_genes_found", analysis['critical_esbl_genes']),
                ("critical_colistin_genes_found", analysis['critical_colistin_genes']),
                ("critical_aminoglycoside_genes_found", analysis['critical_aminoglycoside_genes']),
                ("high_risk_genes_found", analysis['high_risk_resistance_genes']),
                ("virulence_genes_found", analysis['high_risk_virulence_genes'] + analysis['critical_virulence_genes']),
            ]:
                for g in bucket:
                    if g['gene'] not in master_report["risk_assessment"][key]:
                        master_report["risk_assessment"][key].append(g['gene'])

            for class_name, genes in analysis['resistance_classes'].items():
                for g in genes:
                    if g['gene'] not in [x['gene'] for x in master_report["risk_assessment"]["resistance_classes_found"][class_name]]:
                        master_report["risk_assessment"]["resistance_classes_found"][class_name].append(g)

            master_report["detailed_results"][genome_name] = {
                "total_hits": genome_result['total_hits'],
                "database_results": genome_result['results'],
            }

        for db in self.required_databases:
            db_hits = []
            for genome_name, genome_result in all_results.items():
                if db in genome_result['results']:
                    for hit in genome_result['results'][db]['hits']:
                        h = hit.copy()
                        h['genome'] = genome_name
                        db_hits.append(h)
            if db_hits:
                unique_genomes = len(set(h['genome'] for h in db_hits))
                master_report["database_summaries"][db] = {
                    "total_hits": len(db_hits),
                    "unique_genes": len(set(h['gene'] for h in db_hits)),
                    "unique_genomes": unique_genomes,
                    "hits_per_genome": round(len(db_hits) / unique_genomes, 2) if unique_genomes else 0,
                    "most_frequent_genes": dict(Counter([h['gene'] for h in db_hits]).most_common(10)),
                }
            else:
                master_report["database_summaries"][db] = {
                    "total_hits": 0, "unique_genes": 0, "unique_genomes": 0,
                    "hits_per_genome": 0, "most_frequent_genes": {},
                }

        for key in ["critical_carbapenemase_genes_found", "critical_esbl_genes_found",
                    "critical_colistin_genes_found", "critical_aminoglycoside_genes_found",
                    "high_risk_genes_found", "virulence_genes_found"]:
            master_report["risk_assessment"][key].sort()
        master_report["risk_assessment"]["resistance_classes_found"] = dict(
            master_report["risk_assessment"]["resistance_classes_found"]
        )

        master_report["performance_metrics"] = {
            "cpu_cores_used": self.cpus,
            "available_ram_gb": round(self.available_ram, 2),
            "processing_mode": "MAXIMUM SPEED" if self.cpus > 8 else "STANDARD",
        }

        # Optional Excel export
        try:
            excel_data = []
            for genome_name, genome_result in all_results.items():
                all_hits = []
                for r in genome_result['results'].values():
                    all_hits.extend(r['hits'])
                a = self.analyze_acineto_resistance(all_hits)
                excel_data.append({
                    'Genome': genome_name,
                    'Total_Hits': genome_result['total_hits'],
                    'Critical_Carbapenemase': len(a['critical_carbapenemase_genes']),
                    'Critical_ESBL': len(a['critical_esbl_genes']),
                    'Critical_Colistin': len(a['critical_colistin_genes']),
                    'Critical_Aminoglycoside': len(a['critical_aminoglycoside_genes']),
                    'High_Risk_Resistance': len(a['high_risk_resistance_genes']),
                    'Virulence_Factors': len(a['high_risk_virulence_genes']) + len(a['critical_virulence_genes']),
                    'Carbapenemase_Status': a['carbapenemase_status'],
                    'ESBL_Status': a['esbl_status'],
                    'Colistin_Resistance': a['colistin_resistance'],
                })
            df = pd.DataFrame(excel_data)
            excel_file = os.path.join(output_base, "acineto_abricate_master_summary.xlsx")
            df.to_excel(excel_file, index=False)
            master_report["excel_summary_file"] = excel_file
            self.logger.info("Created Excel summary: %s", excel_file)
        except Exception as e:
            self.logger.warning("Could not create Excel summary: %s", e)

        master_file = os.path.join(output_base, "acineto_abricate_master_summary.json")
        with open(master_file, 'w') as f:
            json.dump(master_report, f, indent=4, default=str)
        self.logger.info("Created master JSON summary: %s", master_file)
        return master_file

    # ------------------------------------------------------------------
    # Per-genome orchestration
    # ------------------------------------------------------------------
    def process_single_genome(self, genome_file: str,
                              output_base: str = "acineto_abricate_results") -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        self.logger.info("=== PROCESSING A. BAUMANNII GENOME: %s ===", genome_name)
        os.makedirs(results_dir, exist_ok=True)

        results = {}
        with tempfile.TemporaryDirectory() as temp_dir:
            for db in self.required_databases:
                result = self.run_abricate_single_db(genome_file, db, temp_dir)
                results[db] = result
                status_icon = "OK" if result['status'] == 'success' else "FAIL"
                self.logger.info("%s %s: %d hits", status_icon, db, result['hit_count'])
            self._write_merged_raw_tsv(results, results_dir)

        for r in results.values():
            r.pop('raw_content', None)

        self.create_comprehensive_html_report(genome_name, results, results_dir)

        return {
            'genome': genome_name,
            'results': results,
            'total_hits': sum(r['hit_count'] for r in results.values()),
        }

    def process_multiple_genomes(self, genome_pattern: str,
                                 output_base: str = "acineto_abricate_results") -> Dict[str, Any]:
        # NOTE: ASCII art is printed once from main(); do not duplicate here.
        print("\n" + "=" * 80)
        print("AcinetoScope ABRicate - Analysis Started")
        print(f"   Thresholds: minID={self.min_identity}%, minCov={self.min_coverage}%")
        print("=" * 80)

        if not self.check_abricate_installed():
            raise RuntimeError("ABRicate not properly installed")
        self.setup_abricate_databases()

        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta",
                          f"{genome_pattern}.fa", f"{genome_pattern}.fna",
                          f"{genome_pattern}.faa"]
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        genome_files = list(set(genome_files))

        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")

        self.logger.info("Found %d A. baumannii genomes: %s",
                         len(genome_files), [Path(f).name for f in genome_files])

        os.makedirs(output_base, exist_ok=True)
        all_results = {}

        if len(genome_files) > 1 and self.cpus > 1:
            self.logger.info("Using parallel processing with %d CPU cores", self.cpus)
            with ThreadPoolExecutor(max_workers=self.cpus) as executor:
                future_to_genome = {
                    executor.submit(self.process_single_genome, g, output_base): g
                    for g in genome_files
                }
                for future in as_completed(future_to_genome):
                    genome = future_to_genome[future]
                    try:
                        result = future.result()
                        all_results[Path(genome).stem] = result
                        self.logger.info("Completed: %s (%d total hits)",
                                         result['genome'], result['total_hits'])
                    except Exception as e:
                        self.logger.error("Failed: %s - %s", genome, e)
        else:
            for genome in genome_files:
                try:
                    result = self.process_single_genome(genome, output_base)
                    all_results[Path(genome).stem] = result
                    self.logger.info("Completed: %s (%d total hits)",
                                     result['genome'], result['total_hits'])
                except Exception as e:
                    self.logger.error("Failed: %s - %s", genome, e)

        self.create_database_summaries(all_results, output_base)
        master_json_file = self.create_master_json_summary(all_results, output_base)

        self.logger.info("=== A. BAUMANNII ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        self.logger.info("Master JSON summary: %s", master_json_file)
        return all_results


def main():
    parser = argparse.ArgumentParser(
        description='AcinetoScope ABRicate Analysis - MAXIMUM SPEED VERSION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python acineto_abricate.py "*.fna"
  python acineto_abricate.py "*.fna" --min-identity 90 --min-coverage 85
  python acineto_abricate.py "AB_*.fasta" --cpus 4
  python acineto_abricate.py "*.fa" --output my_acineto_results

OUTPUTS:
  Per genome (2 files each):
    - {genome}_comprehensive_abricate_report.html
    - abricate_raw.tsv                    (all databases merged)

  Batch-level:
    - acineto_{db}_abricate_summary.tsv   (per database)
    - acineto_{db}_summary_report.html    (per database)
    - acineto_abricate_master_summary.json  (single batch JSON)
    - acineto_abricate_master_summary.xlsx  (optional Excel)
        """
    )
    parser.add_argument('pattern', help='File pattern for A. baumannii genomes')
    parser.add_argument('--cpus', '-c', type=int, default=None,
                        help='Number of CPU cores to use (default: auto-detect)')
    parser.add_argument('--output', '-o', default='acineto_abricate_results',
                        help='Output directory (default: acineto_abricate_results)')
    parser.add_argument('--min-identity', type=float, default=80,
                        help='Minimum identity percentage (default: 80)')
    parser.add_argument('--min-coverage', type=float, default=80,
                        help='Minimum coverage percentage (default: 80)')
    args = parser.parse_args()

    banner_art = r"""
 █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
""".strip("\n")
    print("\n" + "=" * 80)
    print(banner_art)
    print("=" * 80)
    print("A. baumannii ABRicate Analysis")
    print("=" * 80)
    print("Author: Brown Beckley | Email: brownbeckley94@gmail.com")
    print("Affiliation: University of Ghana Medical School - Department of Medical Biochemistry")
    print("=" * 80)

    executor = AcinetoAbricateExecutor(cpus=args.cpus,
                                       min_identity=args.min_identity,
                                       min_coverage=args.min_coverage)

    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("AcinetoScope FINAL SUMMARY")
        executor.logger.info("=" * 50)

        total_carba = total_esbl = total_colistin = 0
        for genome_name, result in results.items():
            all_hits = []
            for r in result['results'].values():
                all_hits.extend(r['hits'])
            a = executor.analyze_acineto_resistance(all_hits)
            executor.logger.info(
                "%s: %d total hits, %d carbapenemase, %d ESBL, %d colistin",
                genome_name, result['total_hits'],
                len(a['critical_carbapenemase_genes']),
                len(a['critical_esbl_genes']),
                len(a['critical_colistin_genes']),
            )
            total_carba += len(a['critical_carbapenemase_genes'])
            total_esbl += len(a['critical_esbl_genes'])
            total_colistin += len(a['critical_colistin_genes'])

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("OUTPUT FILES SUMMARY")
        executor.logger.info("=" * 50)
        executor.logger.info("Per genome (2 files each):")
        executor.logger.info("  - {genome}_comprehensive_abricate_report.html")
        executor.logger.info("  - abricate_raw.tsv  (all databases merged)")
        executor.logger.info("Batch-level:")
        executor.logger.info("  - acineto_{db}_abricate_summary.tsv")
        executor.logger.info("  - acineto_{db}_summary_report.html")
        executor.logger.info("  - acineto_abricate_master_summary.json  (single batch JSON)")
        executor.logger.info("  - acineto_abricate_master_summary.xlsx  (Excel)")
        executor.logger.info("Location: %s", args.output)

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("PERFORMANCE SUMMARY")
        executor.logger.info("=" * 50)
        executor.logger.info("CPU cores utilized: %d", executor.cpus)
        executor.logger.info("Available RAM: %.1f GB", executor.available_ram)
        executor.logger.info("Genomes processed: %d", len(results))
        executor.logger.info("Carbapenemase genes: %d", total_carba)
        executor.logger.info("ESBL genes: %d", total_esbl)
        executor.logger.info("Colistin resistance genes: %d", total_colistin)

        import random
        executor.logger.info("\n%s", random.choice(executor.science_quotes))

    except Exception as e:
        executor.logger.error("A. baumannii analysis failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()