#!/usr/bin/env python3
"""
AcinetoScope AMRfinderPlus - A. bamannii AMR Analysis (BUNDLED VERSION with DYNAMIC DATABASE)
Comprehensive AMR analysis with HTML, TSV, and JSON reporting - MAXIMUM SPEED VERSION
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-08-22 (updated 2026-11-01)
Version: 1.4.2

Per-genome outputs (2 files each):
  - {genome}_amrfinder_report.html   (AMR genes + point mutations combined)
  - {genome}_amrfinder_raw.tsv       (merged: Source = AMR / MUTATION)

Batch-level outputs:
  - acineto_amrfinder_summary.tsv
  - acineto_amrfinder_statistics_summary.tsv
  - acineto_amrfinder_summary_report.html
  - acineto_amrfinder_master_summary.json   (SINGLE batch JSON: AMR + mutations)
  - mutation_summary.tsv
  - mutation_summary.html
"""

import subprocess
import sys
import os
import glob
import logging
import shutil
import tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import psutil
import json
import random
from collections import defaultdict, Counter


class AcinetoAMRfinderPlus:
    """
    AMRfinderPlus executor for A. baumannii with dynamic database detection,
    maximum-speed concurrent processing, and comprehensive reporting.
    """

    def __init__(self, cpus: int = None):
        self.logger = self._setup_logging()
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)

        self.bundled_amrfinder = os.path.join(self.module_dir, "bin", "amrfinder")
        self.bundled_update = os.path.join(self.module_dir, "bin", "amrfinder_update")
        self.bundled_database = self._get_latest_database()

        if self.bundled_database is None:
            self.logger.warning("No AMRfinderPlus database found. Please run with --update-db or --force-update.")

        self.metadata = {
            "tool_name": "AcinetoScope AMRfinderPlus",
            "version": "1.4.2",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "amrfinder_version": "4.2.7",
            "database_version": self._get_database_version() if self.bundled_database else "Unknown",
        }

        self.critical_carbapenemases = {
            'blaOXA-23', 'blaOXA-24', 'blaOXA-40', 'blaOXA-51', 'blaOXA-58',
            'blaOXA-66', 'blaOXA-72', 'blaOXA-143', 'blaOXA-235', 'blaOXA-181',
            'blaOXA-232', 'blaOXA-438', 'blaOXA-48',
            'blaKPC', 'blaNDM', 'blaVIM', 'blaIMP', 'blaSIM', 'blaGIM', 'blaSPM'
        }

        self.critical_esbls = {
            'blaCTX-M', 'blaSHV', 'blaTEM', 'blaPER', 'blaVEB', 'blaGES',
            'blaBEL', 'blaBES', 'blaADC', 'blaOXA-10', 'blaOXA-2', 'blaOXA-9'
        }

        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7',
            'mcr-8', 'mcr-9', 'mcr-10',
            'pmrA', 'pmrB', 'pmrC', 'lpxA', 'lpxC', 'lpxD',
            'arnA', 'arnB', 'arnC', 'arnD', 'eptA'
        }

        self.critical_aminoglycoside = {
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG',
            'APH(3\')', 'APH(3\')-VI', 'APH(6)', 'AAC(3)', 'AAC(6\')', 'ANT(2")',
            'ANT(3")', 'ANT(4")', 'aacC1', 'aacC2', 'aacC4', 'aacA4', 'aacA7',
            'aadA1', 'aadA2', 'aadA5', 'aadA7', 'strA', 'strB', 'aphA1', 'aphA2',
            'aphA3', 'aphA6', 'aac3', 'aac6', 'aadA', 'aadB'
        }

        self.high_risk_genes = (
            self.critical_carbapenemases
            | self.critical_esbls
            | self.critical_colistin
            | self.critical_aminoglycoside
            | {
                'tetA', 'tetB', 'tetC', 'tetD', 'tetE', 'tetG', 'tetH', 'tetK', 'tetL',
                'tetM', 'tetO', 'tetQ', 'tetS', 'tetX',
                'sul1', 'sul2', 'sul3', 'sul4',
                'dfrA1', 'dfrA5', 'dfrA7', 'dfrA8', 'dfrA12', 'dfrA14', 'dfrA17',
                'dfrA19', 'dfrA20', 'dfrA21', 'dfrB1', 'dfrB2', 'dfrB3', 'dfrB4',
                'dfrB5', 'dfrB6', 'dfrB7',
                'catA1', 'catA2', 'catB2', 'catB3', 'catB8', 'catI', 'catII', 'catIII',
                'cmlA', 'cmlA1', 'cmlA5', 'cmlA6', 'cmlA7', 'floR',
                'ermA', 'ermB', 'ermC', 'ermF', 'ermG', 'ermX', 'ermY',
                'mphA', 'mphB', 'mphC', 'mphD', 'mphE', 'msrA', 'msrB', 'msrC', 'msrD',
                'qnrA1', 'qnrB1', 'qnrB2', 'qnrB4', 'qnrB6', 'qnrB10', 'qnrB19',
                'qnrS1', 'qnrS2', 'qnrVC1', 'qnrVC4',
                'aac(6\')-Ib-cr', 'qepA', 'qepA1', 'qepA2', 'qepA3', 'qepA4',
                'fosA', 'fosB', 'fosC', 'fosX',
                'arr-2', 'arr-3', 'arr-4', 'arr-5', 'arr-6', 'arr-7',
                'adeA', 'adeB', 'adeC', 'adeG', 'adeH', 'adeI', 'adeJ', 'adeK',
                'adeL', 'adeM', 'adeN',
                'abeM', 'abeS', 'acrA', 'acrB', 'tolC', 'mexA', 'mexB', 'mexC',
                'mexD', 'mexE', 'mexF', 'oprM', 'oprN', 'oprJ',
                'mdtA', 'mdtB', 'mdtC', 'emrA', 'emrB', 'emrD', 'emrE', 'emrK', 'emrY'
            }
        )

        self.critical_risk_genes = (
            self.critical_carbapenemases
            | self.critical_colistin
            | {'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5'}
        )

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
            {"text": "AcinetoScope represents the convergence of genomic surveillance and clinical diagnostics.", "author": "Brown Beckley"},
        ]

        self.ascii_art = r"""█████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝""".strip("\n")

    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)

    def _get_available_ram(self) -> float:
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except Exception:
            return 8

    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        if user_cpus is not None:
            self._log_resource_info(user_cpus)
            return user_cpus
        try:
            total_physical = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            if total_physical <= 4:
                optimal = total_physical
            elif total_physical <= 8:
                optimal = total_physical - 1
            elif total_physical <= 16:
                optimal = max(8, total_physical - 1)
            elif total_physical <= 32:
                optimal = max(16, total_physical - 3)
            else:
                optimal = min(32, int(total_physical * 0.95))
            optimal = max(1, min(optimal, total_physical))
            self._log_resource_info(optimal, total_physical)
            return optimal
        except Exception:
            return os.cpu_count() or 4

    def _log_resource_info(self, cpus: int, total_cores: int = None):
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            self.logger.info(f"Using CPU cores: {cpus} ({cpus/total_cores*100:.1f}% of available cores)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        if cpus <= 4:
            self.logger.info("Performance: multi-core (max speed for small systems)")
        elif cpus <= 8:
            self.logger.info("Performance: high-speed mode")
        else:
            self.logger.info("Performance: MAXIMUM SPEED MODE")

    def _get_latest_database(self) -> Optional[str]:
        db_root = os.path.join(self.module_dir, "data", "amrfinder_db")
        if not os.path.exists(db_root):
            return None
        candidates = []
        for item in os.listdir(db_root):
            full_path = os.path.join(db_root, item)
            if os.path.isdir(full_path) and item.startswith('20'):
                candidates.append(item)
        if not candidates:
            return None
        latest = sorted(candidates)[-1]
        return os.path.join(db_root, latest)

    def _get_database_version(self) -> str:
        if not self.bundled_database:
            return "Unknown"
        version_file = os.path.join(self.bundled_database, "version.txt")
        if os.path.exists(version_file):
            with open(version_file, 'r') as f:
                return f.read().strip()
        return os.path.basename(self.bundled_database)

    def update_database(self, force: bool = False) -> bool:
        if not os.path.exists(self.bundled_update):
            self.logger.error(f"amrfinder_update not found at {self.bundled_update}")
            return False
        if not os.access(self.bundled_update, os.X_OK):
            os.chmod(self.bundled_update, 0o755)

        db_dir = os.path.join(self.module_dir, "data", "amrfinder_db")
        os.makedirs(db_dir, exist_ok=True)

        if force:
            self.logger.info("Forcing full database update - removing old folders...")
            for item in os.listdir(db_dir):
                full_path = os.path.join(db_dir, item)
                if os.path.isdir(full_path) and item.startswith('20'):
                    self.logger.info(f"  Removing {item}")
                    shutil.rmtree(full_path)

        self.logger.info("Updating AMRfinderPlus database...")
        try:
            subprocess.run([self.bundled_update, "--database", db_dir],
                           capture_output=True, text=True, check=True)
            self.logger.info("Database update completed successfully.")
            self.bundled_database = self._get_latest_database()
            if self.bundled_database:
                self.metadata['database_version'] = self._get_database_version()
                self.logger.info(f"New database version: {self.metadata['database_version']}")
                return True
            self.logger.error("Database update succeeded but no database folder found.")
            return False
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Database update failed: {e.stderr}")
            return False

    def check_amrfinder_installed(self) -> bool:
        try:
            if not os.path.exists(self.bundled_amrfinder):
                self.logger.error(f"Bundled AMRfinderPlus not found at {self.bundled_amrfinder}")
                return False
            if not os.access(self.bundled_amrfinder, os.X_OK):
                os.chmod(self.bundled_amrfinder, 0o755)
            subprocess.run([self.bundled_amrfinder, '--version'],
                           capture_output=True, text=True, check=True)
            self.logger.info("Bundled AMRfinderPlus: OK")
            if self.bundled_database and os.path.exists(self.bundled_database):
                self.logger.info(f"Bundled database found: {self.bundled_database}")
                return True
            self.logger.warning("Bundled database not found. Please run with --update-db or --force-update.")
            return False
        except Exception as e:
            self.logger.error(f"AMRfinderPlus check failed: {e}")
            return False

    # ------------------------------------------------------------------
    # Shared footer (ABRicate-style, used by every HTML report)
    # ------------------------------------------------------------------
    def _get_footer_html(self, title: str) -> str:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return f"""
    <div class="footer">
            <p><strong>ACINETOSCOPE</strong> &mdash; AMR Analysis Module</p>
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
    # Run + merged raw TSV
    # ------------------------------------------------------------------
    def _write_merged_raw_tsv(self, amr_file: str, mut_file: Optional[str], output_file: str):
        """Merge AMRFinderPlus main output and --mutation_all output into one TSV.

        A leading 'Source' column distinguishes AMR gene hits from MUTATION rows.
        Lines appearing identically in both files are written only once (as AMR).
        """
        amr_lines_seen = set()
        with open(output_file, 'w') as out:
            if amr_file and os.path.exists(amr_file):
                with open(amr_file, 'r') as f:
                    amr_lines = f.readlines()
                if amr_lines:
                    header = amr_lines[0].rstrip('\n')
                    out.write('Source\t' + header + '\n')
                    for line in amr_lines[1:]:
                        if not line.strip():
                            continue
                        amr_lines_seen.add(line.rstrip('\n'))
                        out.write('AMR\t' + line)

            if mut_file and os.path.exists(mut_file):
                with open(mut_file, 'r') as f:
                    mut_lines = f.readlines()
                if mut_lines:
                    for line in mut_lines[1:]:
                        if not line.strip():
                            continue
                        if line.rstrip('\n') in amr_lines_seen:
                            continue
                        out.write('MUTATION\t' + line)

    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str,
                                     min_identity: float = None, min_coverage: float = None,
                                     report_mutations: bool = True) -> Dict[str, Any]:
        """Run AMRfinderPlus on one genome. Emits exactly 2 files per sample:
             - {genome}_amrfinder_raw.tsv  (AMR + MUTATION merged)
             - {genome}_amrfinder_report.html  (comprehensive AMR + mutations)
        """
        genome_name = Path(genome_file).stem
        merged_raw_file = os.path.join(output_dir, f"{genome_name}_amrfinder_raw.tsv")

        if not os.path.exists(self.bundled_amrfinder):
            self.logger.error(f"Bundled AMRfinderPlus not found at {self.bundled_amrfinder}")
            return {
                'genome': genome_name, 'output_file': merged_raw_file,
                'hits': [], 'hit_count': 0,
                'mutations': [], 'mutation_count': 0,
                'status': 'failed', 'error': 'AMRfinder binary not found',
            }

        with tempfile.TemporaryDirectory() as temp_dir:
            amr_output_file = os.path.join(temp_dir, f"{genome_name}_amrfinder.txt")
            mut_file = os.path.join(temp_dir, f"{genome_name}_mutations.tsv") if report_mutations else None

            cmd = [
                self.bundled_amrfinder,
                '-n', genome_file,
                '--output', amr_output_file,
                '--threads', str(self.cpus),
                '--organism', 'Acinetobacter_baumannii',
            ]

            if self.bundled_database and os.path.exists(self.bundled_database):
                cmd.extend(['--database', self.bundled_database])
                self.logger.info(f"Using bundled database: {self.bundled_database}")
            else:
                self.logger.warning("Using default AMRfinderPlus database location")

            if min_identity is not None:
                cmd.extend(['--ident_min', str(min_identity)])
            if min_coverage is not None:
                cmd.extend(['--coverage_min', str(min_coverage)])
            if report_mutations and mut_file:
                cmd.extend(['--mutation_all', mut_file])

            self.logger.info(f"Running BUNDLED AMRfinderPlus: {genome_name} (using {self.cpus} CPU cores)")
            try:
                subprocess.run(cmd, capture_output=True, text=True, check=True)

                hits = self._parse_amrfinder_output(amr_output_file)
                mutations = (self._parse_mutations_file(mut_file)
                             if mut_file and os.path.exists(mut_file) else [])

                self._write_merged_raw_tsv(amr_output_file, mut_file, merged_raw_file)
                self._create_amrfinder_html_report(genome_name, hits, mutations, output_dir)

                return {
                    'genome': genome_name,
                    'output_file': merged_raw_file,
                    'hits': hits,
                    'hit_count': len(hits),
                    'mutations': mutations,
                    'mutation_count': len(mutations),
                    'status': 'success',
                }
            except subprocess.CalledProcessError as e:
                self.logger.error(f"AMRfinderPlus failed for {genome_name}")
                self.logger.error(f"STDERR: {e.stderr}")
                return {
                    'genome': genome_name,
                    'output_file': merged_raw_file,
                    'hits': [], 'hit_count': 0,
                    'mutations': [], 'mutation_count': 0,
                    'status': 'failed', 'error': e.stderr,
                }

    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        hits = []
        try:
            with open(amrfinder_file, 'r') as f:
                lines = f.readlines()
            if len(lines) < 2:
                return hits
            headers = lines[0].strip().split('\t')
            for line in lines[1:]:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    hit = dict(zip(headers, parts))
                    processed = {
                        'protein_id': hit.get('Protein id', ''),
                        'contig_id': hit.get('Contig id', ''),
                        'start': hit.get('Start', ''),
                        'stop': hit.get('Stop', ''),
                        'strand': hit.get('Strand', ''),
                        'gene_symbol': hit.get('Element symbol', ''),
                        'sequence_name': hit.get('Element name', ''),
                        'scope': hit.get('Scope', ''),
                        'element_type': hit.get('Type', ''),
                        'element_subtype': hit.get('Subtype', ''),
                        'class': hit.get('Class', ''),
                        'subclass': hit.get('Subclass', ''),
                        'method': hit.get('Method', ''),
                        'target_length': hit.get('Target length', ''),
                        'ref_length': hit.get('Reference sequence length', ''),
                        'coverage': hit.get('% Coverage of reference', '').replace('%', ''),
                        'identity': hit.get('% Identity to reference', '').replace('%', ''),
                        'alignment_length': hit.get('Alignment length', ''),
                        'accession': hit.get('Closest reference accession', ''),
                        'closest_name': hit.get('Closest reference name', ''),
                        'hmm_id': hit.get('HMM accession', ''),
                        'hmm_description': hit.get('HMM description', ''),
                    }
                    hits.append(processed)
        except Exception as e:
            self.logger.error(f"Error parsing {amrfinder_file}: {e}")
        self.logger.info(f"Parsed {len(hits)} AMR hits from {amrfinder_file}")
        return hits

    def _parse_mutations_file(self, mut_file: str) -> List[Dict]:
        return self._parse_amrfinder_output(mut_file)

    # ------------------------------------------------------------------
    # HTML report (per sample, combined AMR + mutations)
    # ------------------------------------------------------------------
    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict],
                                      mutations: List[Dict], output_dir: str):
        analysis = self._analyze_acineto_amr_results(hits)
        random_quote = random.choice(self.science_quotes)

        # --- Build mutation rows ---
        mutation_section = ""
        if mutations:
            mut_rows = ""
            for m in mutations:
                mut_rows += f"""
                    <tr>
                        <td>{m.get('gene_symbol', '')}</td>
                        <td>{m.get('sequence_name', '')}</td>
                        <td>{m.get('class', '')}</td>
                        <td>{m.get('subclass', '')}</td>
                        <td>{m.get('contig_id', '')}</td>
                        <td>{m.get('start', '')}</td>
                        <td>{m.get('stop', '')}</td>
                        <td>{m.get('strand', '')}</td>
                        <td>{m.get('coverage', '')}</td>
                        <td>{m.get('identity', '')}</td>
                        <td>{m.get('accession', '')}</td>
                    </tr>
                """
            mutation_section = f"""
    <div class="card">
        <h2>Point Mutations Detected</h2>
        <p>All point mutations detected by AMRFinderPlus (including synonymous variants).</p>
        <div class="table-responsive">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Gene Symbol</th><th>Mutation</th><th>Class</th><th>Subclass</th>
                        <th>Contig</th><th>Start</th><th>Stop</th><th>Strand</th>
                        <th>Coverage (%)</th><th>Identity (%)</th><th>Accession</th>
                    </tr>
                </thead>
                <tbody>
                    {mut_rows}
                </tbody>
            </table>
        </div>
    </div>
    """

        # --- Build AMR section ---
        amr_section = ""
        if hits:
            amr_rows = ""
            for hit in hits:
                if hit['gene_symbol'] in analysis['critical_risk_list']:
                    row_class = "critical-row"
                elif hit['gene_symbol'] in analysis['high_risk_list']:
                    row_class = "high-risk-row"
                else:
                    row_class = "present"
                amr_rows += f"""
                    <tr class="{row_class}">
                        <td><strong>{hit['gene_symbol']}</strong></td>
                        <td>{hit['sequence_name']}</td>
                        <td>{hit['class']}</td>
                        <td>{hit['subclass']}</td>
                        <td>{hit['coverage']}%</td>
                        <td>{hit['identity']}%</td>
                        <td>{hit['scope']}</td>
                    </tr>
                """
            amr_section = f"""
    <div class="card">
        <h2>Detailed AMR Genes</h2>
        <div class="table-responsive">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Gene</th><th>Sequence Name</th><th>Class</th>
                        <th>Subclass</th><th>Coverage</th><th>Identity</th><th>Scope</th>
                    </tr>
                </thead>
                <tbody>{amr_rows}</tbody>
            </table>
        </div>
    </div>
    """
        else:
            amr_section = """
    <div class="card">
        <h2>No AMR Genes Detected</h2>
        <p>No antimicrobial resistance genes found in this A. baumannii genome.</p>
    </div>
    """

        # --- Alerts ---
        critical_alert = ""
        if analysis['critical_risk_genes'] > 0:
            badges = "".join(f'<span class="risk-badge">{g}</span>'
                             for g in analysis['critical_risk_list'])
            critical_alert = f"""
    <div class="card" style="border-left: 4px solid #dc2626;">
        <h2 style="color: #dc2626;">CRITICAL RISK AMR GENES DETECTED</h2>
        <p><strong>{analysis['critical_risk_genes']} critical risk genes:</strong></p>
        <div style="margin: 10px 0;">{badges}</div>
    </div>
    """

        high_risk_alert = ""
        if analysis['high_risk_genes'] > 0 and analysis['critical_risk_genes'] == 0:
            badges = "".join(f'<span class="warning-badge">{g}</span>'
                             for g in analysis['high_risk_list'])
            high_risk_alert = f"""
    <div class="card" style="border-left: 4px solid #f59e0b;">
        <h2 style="color: #f59e0b;">High-Risk AMR Genes</h2>
        <p><strong>{analysis['high_risk_genes']} high-risk genes:</strong></p>
        <div style="margin: 10px 0;">{badges}</div>
    </div>
    """

        # --- Mechanisms ---
        mech_section = ""
        if any(analysis['resistance_mechanisms'].values()):
            blocks = ""
            for mechanism, genes in analysis['resistance_mechanisms'].items():
                if genes:
                    label = mechanism.replace('_', ' ').title()
                    color = "#fee2e2" if mechanism in ('carbapenemase', 'colistin_resistance') else "#fef3c7"
                    blocks += f"""
        <div style="margin: 10px 0; padding: 10px; background: {color}; border-radius: 5px;">
            <strong>{label}:</strong> {', '.join(genes)}
        </div>
        """
            mech_section = f"""
    <div class="card">
        <h2>Resistance Mechanism Breakdown</h2>
        {blocks}
    </div>
    """

        # --- Resistance classes ---
        classes_section = ""
        if analysis['resistance_classes']:
            class_rows = ""
            for cls, genes in analysis['resistance_classes'].items():
                class_rows += f"""
                    <tr>
                        <td><strong>{cls}</strong></td>
                        <td>{len(genes)}</td>
                        <td>{', '.join(genes)}</td>
                    </tr>
                """
            classes_section = f"""
    <div class="card">
        <h2>Resistance Classes</h2>
        <div class="table-responsive">
            <table class="data-table">
                <thead><tr><th>Class</th><th>Count</th><th>Genes</th></tr></thead>
                <tbody>{class_rows}</tbody>
            </table>
        </div>
    </div>
    """

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AcinetoScope - AMRfinderPlus Report: {genome_name}</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{ background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%); font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; color: #ffffff; padding: 20px; min-height: 100vh; }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{ background: rgba(0,0,0,0.7); padding: 20px; border-radius: 15px; margin-bottom: 20px; box-shadow: 0 8px 32px rgba(0,0,0,0.4); border: 2px solid rgba(0,255,0,0.3); }}
        .ascii-art {{ font-family: 'Courier New', monospace; font-size: 10px; line-height: 1.1; white-space: pre; color: #00ff00; text-shadow: 0 0 10px rgba(0,255,0,0.5); overflow-x: auto; }}
        .quote-container {{ background: rgba(255,255,255,0.1); backdrop-filter: blur(10px); padding: 20px; border-radius: 10px; margin-bottom: 30px; text-align: center; min-height: 100px; display: flex; flex-direction: column; justify-content: center; box-shadow: 0 4px 20px rgba(0,0,0,0.3); border: 1px solid rgba(255,255,255,0.2); }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; color: #ffffff; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .card {{ background: rgba(255,255,255,0.95); color: #1f2937; padding: 25px; border-radius: 10px; margin-bottom: 20px; box-shadow: 0 4px 15px rgba(0,0,0,0.2); }}
        .card h2 {{ color: #1e3a8a; border-bottom: 3px solid #3b82f6; padding-bottom: 10px; margin-bottom: 20px; font-size: 24px; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-bottom: 20px; }}
        .stat-card {{ background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%); color: white; padding: 15px; border-radius: 8px; text-align: center; box-shadow: 0 4px 12px rgba(0,0,0,0.15); }}
        .stat-card.critical-card {{ background: linear-gradient(135deg, #dc2626 0%, #b91c1c 100%); }}
        .stat-value {{ font-size: 24px; font-weight: bold; margin-bottom: 5px; }}
        .stat-label {{ font-size: 12px; opacity: 0.9; }}
        .risk-badge {{ display: inline-block; background: #dc2626; color: white; padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em; }}
        .warning-badge {{ display: inline-block; background: #f59e0b; color: black; padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em; }}
        .table-responsive {{ width: 100%; overflow-x: auto; margin: 20px 0; }}
        .data-table {{ width: 100%; border-collapse: collapse; font-size: 13px; min-width: 900px; }}
        .data-table th {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 12px; text-align: left; font-weight: bold; }}
        .data-table td {{ padding: 10px; border-bottom: 1px solid #e5e7eb; }}
        .data-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .critical-row {{ background-color: #fee2e2; font-weight: bold; border-left: 4px solid #dc2626; }}
        .high-risk-row {{ background-color: #fef3c7; border-left: 4px solid #f59e0b; }}
        .present {{ background-color: #d1fae5; }}
        .footer {{ text-align: center; margin-top: 30px; padding: 20px; background: rgba(0,0,0,0.3); border-radius: 10px; font-size: 14px; }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .authorship {{ margin-top: 15px; padding: 15px; background: rgba(255,255,255,0.1); border-radius: 8px; font-size: 12px; }}
        .authorship a {{ color: #fbbf24; text-decoration: none; }}
        .authorship a:hover {{ text-decoration: underline; }}
        @media (max-width: 768px) {{ .ascii-art {{ font-size: 6px; }} .stats-grid {{ grid-template-columns: 1fr; }} }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="quote-container">
            <div class="quote-text">"{random_quote['text']}"</div>
            <div class="quote-author">&mdash; {random_quote['author']}</div>
        </div>
    </div>

    <div class="card">
        <h2>AMRfinderPlus Analysis: {genome_name}</h2>
        <div class="stats-grid">
            <div class="stat-card"><div class="stat-value">{analysis['total_genes']}</div><div class="stat-label">Total AMR Genes</div></div>
            <div class="stat-card"><div class="stat-value">{analysis['high_risk_genes']}</div><div class="stat-label">High Risk Genes</div></div>
            <div class="stat-card critical-card"><div class="stat-value">{analysis['critical_risk_genes']}</div><div class="stat-label">Critical Risk</div></div>
            <div class="stat-card"><div class="stat-value">{len(mutations)}</div><div class="stat-label">Point Mutations</div></div>
        </div>
        <p><strong>Genome:</strong> {genome_name}</p>
        <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
        <p><strong>AMRfinderPlus:</strong> {self.metadata['amrfinder_version']}</p>
        <p><strong>Database:</strong> {self.metadata['database_version']}</p>
        <p><strong>Analysis Date:</strong> {datetime.now().strftime('%Y-%m-%d')}</p>
    </div>

    {critical_alert}
    {high_risk_alert}
    {mech_section}
    {classes_section}
    {amr_section}
    {mutation_section}

    {self._get_footer_html("AcinetoScope - AMRfinderPlus Per-Sample Report")}
</div>
</body>
</html>"""

        out_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"AMR HTML report: {out_file}")

    def _analyze_acineto_amr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        analysis = {
            'total_genes': len(hits),
            'resistance_classes': {},
            'high_risk_genes': 0,
            'critical_risk_genes': 0,
            'high_risk_list': [],
            'critical_risk_list': [],
            'resistance_mechanisms': {
                'carbapenemase': [],
                'esbl': [],
                'colistin_resistance': [],
                'fluoroquinolone_resistance': [],
                'aminoglycoside_resistance': [],
                'efflux_pumps': [],
                'other_amr': [],
            },
        }
        for hit in hits:
            gene = hit.get('gene_symbol', '')
            class_name = hit.get('class', '')
            if not gene:
                continue
            self._categorize_acineto_mechanism(gene, class_name, analysis)
            if gene in self.critical_risk_genes:
                analysis['critical_risk_genes'] += 1
                if gene not in analysis['critical_risk_list']:
                    analysis['critical_risk_list'].append(gene)
            if gene in self.high_risk_genes:
                analysis['high_risk_genes'] += 1
                if gene not in analysis['high_risk_list']:
                    analysis['high_risk_list'].append(gene)
            if class_name:
                analysis['resistance_classes'].setdefault(class_name, []).append(gene)
        return analysis

    def _categorize_acineto_mechanism(self, gene_symbol: str, resistance_class: str, analysis: Dict):
        gene_lower = gene_symbol.lower()
        if any(c in gene_lower for c in ['oxa', 'kpc', 'ndm', 'vim', 'imp', 'sim']):
            if gene_symbol in self.critical_carbapenemases:
                analysis['resistance_mechanisms']['carbapenemase'].append(gene_symbol)
        elif any(c in gene_lower for c in ['ctx', 'shv', 'tem', 'per', 'veb', 'ges']):
            if gene_symbol in self.critical_esbls:
                analysis['resistance_mechanisms']['esbl'].append(gene_symbol)
        elif 'mcr' in gene_lower:
            if gene_symbol in self.critical_colistin:
                analysis['resistance_mechanisms']['colistin_resistance'].append(gene_symbol)
        elif any(c in gene_lower for c in ['arm', 'rmt', 'aph', 'aac', 'aad', 'str', 'ant']):
            if gene_symbol in self.critical_aminoglycoside:
                analysis['resistance_mechanisms']['aminoglycoside_resistance'].append(gene_symbol)
        elif any(c in gene_lower for c in ['qnr', 'qep', 'gyra', 'parc']):
            analysis['resistance_mechanisms']['fluoroquinolone_resistance'].append(gene_symbol)
        elif any(c in gene_lower for c in ['ade', 'abe', 'mex', 'acr', 'emr', 'mdt']):
            analysis['resistance_mechanisms']['efflux_pumps'].append(gene_symbol)
        else:
            analysis['resistance_mechanisms']['other_amr'].append(gene_symbol)

    # ------------------------------------------------------------------
    # Orchestration
    # ------------------------------------------------------------------
    def process_single_genome(self, genome_file: str, output_base: str = "acineto_amrfinder_results",
                              min_identity: float = None, min_coverage: float = None,
                              report_mutations: bool = True) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        os.makedirs(results_dir, exist_ok=True)
        self.logger.info(f"=== PROCESSING GENOME: {genome_name} ===")
        result = self.run_amrfinder_single_genome(genome_file, results_dir,
                                                  min_identity, min_coverage, report_mutations)
        status = "OK" if result['status'] == 'success' else "FAIL"
        self.logger.info(f"{status} {genome_name}: {result['hit_count']} AMR hits, "
                         f"{result.get('mutation_count', 0)} mutations")
        return result

    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "acineto_amrfinder_results",
                                 min_identity: float = None, min_coverage: float = None,
                                 report_mutations: bool = True) -> Dict[str, Any]:
        if not self.check_amrfinder_installed():
            raise RuntimeError("BUNDLED AMRfinderPlus not properly installed")

        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa",
                          f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        genome_files = list(set(genome_files))

        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")

        self.logger.info(f"Found {len(genome_files)} genomes: {[Path(f).name for f in genome_files]}")
        os.makedirs(output_base, exist_ok=True)

        all_results = {}
        max_concurrent = max(1, min(self.cpus, len(genome_files), int(self.available_ram / 2.5)))
        self.logger.info(f"Using {max_concurrent} concurrent jobs "
                         f"(each AMRfinderPlus instance uses {self.cpus} threads internally)")

        with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            futures = {
                executor.submit(self.process_single_genome, genome, output_base,
                                min_identity, min_coverage, report_mutations): genome
                for genome in genome_files
            }
            for future in as_completed(futures):
                genome = futures[future]
                try:
                    result = future.result()
                    all_results[result['genome']] = result
                    self.logger.info(f"COMPLETED: {result['genome']} "
                                     f"({result['hit_count']} AMR hits, "
                                     f"{result.get('mutation_count', 0)} mutations)")
                except Exception as e:
                    self.logger.error(f"FAILED: {genome} - {e}")
                    all_results[Path(genome).stem] = {
                        'genome': Path(genome).stem,
                        'hits': [], 'hit_count': 0,
                        'mutations': [], 'mutation_count': 0,
                        'status': 'failed',
                    }

        self.create_amr_summary(all_results, output_base)
        self.create_mutation_summary(all_results, output_base)

        self.logger.info("=== A. BAUMANNII AMR ANALYSIS COMPLETE ===")
        self.logger.info(f"Processed {len(all_results)} genomes")
        self.logger.info(f"Results saved to: {output_base}")
        return all_results

    # ------------------------------------------------------------------
    # Batch summaries
    # ------------------------------------------------------------------
    def create_amr_summary(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating A. baumannii AMR summary files...")
        summary_file = os.path.join(output_base, "acineto_amrfinder_summary.tsv")
        with open(summary_file, 'w') as f:
            f.write("Genome\tProtein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\t"
                    "Scope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t"
                    "% Coverage of reference\t% Identity to reference\tAlignment length\t"
                    "Closest reference accession\tClosest reference name\tHMM accession\tHMM description\n")
            for genome_name, result in all_results.items():
                for hit in result['hits']:
                    row = [
                        genome_name,
                        hit.get('protein_id', ''),
                        hit.get('contig_id', ''),
                        hit.get('start', ''),
                        hit.get('stop', ''),
                        hit.get('strand', ''),
                        hit.get('gene_symbol', ''),
                        hit.get('sequence_name', ''),
                        hit.get('scope', ''),
                        hit.get('element_type', ''),
                        hit.get('element_subtype', ''),
                        hit.get('class', ''),
                        hit.get('subclass', ''),
                        hit.get('method', ''),
                        hit.get('target_length', ''),
                        hit.get('ref_length', ''),
                        hit.get('coverage', ''),
                        hit.get('identity', ''),
                        hit.get('alignment_length', ''),
                        hit.get('accession', ''),
                        hit.get('closest_name', ''),
                        hit.get('hmm_id', ''),
                        hit.get('hmm_description', ''),
                    ]
                    f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info(f"AMR summary TSV: {summary_file}")

        stats_file = os.path.join(output_base, "acineto_amrfinder_statistics_summary.tsv")
        with open(stats_file, 'w') as f:
            f.write("Genome\tTotal_AMR_Genes\tHigh_Risk_Genes\tCritical_Risk_Genes\t"
                    "Resistance_Classes\tTotal_Mutations\tGene_List\n")
            for genome_name, result in all_results.items():
                genes = list(set(hit.get('gene_symbol', '')
                                 for hit in result['hits'] if hit.get('gene_symbol')))
                gene_list = ",".join(genes)
                high_risk_count = sum(1 for g in genes if g in self.high_risk_genes)
                critical_risk_count = sum(1 for g in genes if g in self.critical_risk_genes)
                classes = list(set(hit.get('class', '')
                                   for hit in result['hits'] if hit.get('class')))
                class_list = ",".join(classes)
                mutation_count = len(result.get('mutations', []))
                f.write(f"{genome_name}\t{result['hit_count']}\t{high_risk_count}\t"
                        f"{critical_risk_count}\t{class_list}\t{mutation_count}\t{gene_list}\n")
        self.logger.info(f"Statistics summary: {stats_file}")

        self._create_amr_summary_html(all_results, output_base)
        self._create_amr_master_json(all_results, output_base)

    def _create_amr_summary_html(self, all_results: Dict[str, Any], output_base: str):
        all_hits = []
        for genome_name, result in all_results.items():
            for hit in result['hits']:
                hit_with_genome = hit.copy()
                hit_with_genome['genome'] = genome_name
                all_hits.append(hit_with_genome)

        total_genomes = len(all_results)
        total_hits = len(all_hits)
        total_mutations = sum(len(r.get('mutations', [])) for r in all_results.values())

        critical_genes_found = set()
        high_risk_genes_found = set()
        genomes_with_critical = 0
        genomes_with_high_risk = 0

        gene_frequency = {}
        for hit in all_hits:
            gene = hit.get('gene_symbol', '')
            if not gene:
                continue
            gene_frequency.setdefault(gene, set()).add(hit['genome'])

        for genome_name, result in all_results.items():
            genome_genes = {hit['gene_symbol'] for hit in result['hits'] if hit['gene_symbol']}
            if any(g in genome_genes for g in self.critical_risk_genes):
                genomes_with_critical += 1
                critical_genes_found.update(genome_genes.intersection(self.critical_risk_genes))
            if any(g in genome_genes for g in self.high_risk_genes):
                genomes_with_high_risk += 1
                high_risk_genes_found.update(genome_genes.intersection(self.high_risk_genes))

        random_quote = random.choice(self.science_quotes)

        critical_alert = ""
        if critical_genes_found:
            badges = "".join(f'<span class="risk-badge">{g}</span>'
                             for g in sorted(critical_genes_found))
            critical_alert = f"""
    <div class="card" style="border-left: 4px solid #dc2626;">
        <h2 style="color: #dc2626;">CRITICAL RISK GENES ACROSS ALL GENOMES</h2>
        <p><strong>{len(critical_genes_found)} unique critical genes in
        {genomes_with_critical} genomes:</strong></p>
        <div style="margin: 10px 0;">{badges}</div>
    </div>
    """

        high_risk_alert = ""
        if high_risk_genes_found and not critical_genes_found:
            badges = "".join(f'<span class="warning-badge">{g}</span>'
                             for g in sorted(high_risk_genes_found))
            high_risk_alert = f"""
    <div class="card" style="border-left: 4px solid #f59e0b;">
        <h2 style="color: #f59e0b;">High-Risk Genes</h2>
        <p><strong>{len(high_risk_genes_found)} unique high-risk genes:</strong></p>
        <div style="margin: 10px 0;">{badges}</div>
    </div>
    """

        freq_rows = ""
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            freq = len(genomes)
            pct = freq / total_genomes * 100 if total_genomes else 0
            genome_list = ', '.join(sorted(genomes))
            freq_rows += f"""
                    <tr>
                        <td><strong>{gene}</strong></td>
                        <td>{freq} ({pct:.1f}%)</td>
                        <td>{genome_list}</td>
                    </tr>
            """

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AcinetoScope - AMRfinderPlus Batch Summary</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{ background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%); font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; color: #ffffff; padding: 20px; min-height: 100vh; }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{ background: rgba(0,0,0,0.7); padding: 20px; border-radius: 15px; margin-bottom: 20px; box-shadow: 0 8px 32px rgba(0,0,0,0.4); border: 2px solid rgba(0,255,0,0.3); }}
        .ascii-art {{ font-family: 'Courier New', monospace; font-size: 10px; line-height: 1.1; white-space: pre; color: #00ff00; text-shadow: 0 0 10px rgba(0,255,0,0.5); overflow-x: auto; }}
        .quote-container {{ background: rgba(255,255,255,0.1); backdrop-filter: blur(10px); padding: 20px; border-radius: 10px; margin-bottom: 30px; text-align: center; min-height: 100px; display: flex; flex-direction: column; justify-content: center; box-shadow: 0 4px 20px rgba(0,0,0,0.3); border: 1px solid rgba(255,255,255,0.2); }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; color: #ffffff; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .card {{ background: rgba(255,255,255,0.95); color: #1f2937; padding: 25px; border-radius: 10px; margin-bottom: 20px; box-shadow: 0 4px 15px rgba(0,0,0,0.2); }}
        .card h2 {{ color: #1e3a8a; border-bottom: 3px solid #3b82f6; padding-bottom: 10px; margin-bottom: 20px; font-size: 24px; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-bottom: 20px; }}
        .stat-card {{ background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%); color: white; padding: 15px; border-radius: 8px; text-align: center; }}
        .stat-card.critical-card {{ background: linear-gradient(135deg, #dc2626 0%, #b91c1c 100%); }}
        .stat-value {{ font-size: 24px; font-weight: bold; }}
        .stat-label {{ font-size: 12px; opacity: 0.9; }}
        .risk-badge {{ display: inline-block; background: #dc2626; color: white; padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em; }}
        .warning-badge {{ display: inline-block; background: #f59e0b; color: black; padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em; }}
        .table-responsive {{ width: 100%; overflow-x: auto; margin: 20px 0; }}
        .data-table {{ width: 100%; border-collapse: collapse; font-size: 13px; min-width: 900px; }}
        .data-table th {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 12px; text-align: left; font-weight: bold; }}
        .data-table td {{ padding: 10px; border-bottom: 1px solid #e5e7eb; }}
        .data-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .footer {{ text-align: center; margin-top: 30px; padding: 20px; background: rgba(0,0,0,0.3); border-radius: 10px; font-size: 14px; }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .authorship {{ margin-top: 15px; padding: 15px; background: rgba(255,255,255,0.1); border-radius: 8px; font-size: 12px; }}
        .authorship a {{ color: #fbbf24; text-decoration: none; }}
        .authorship a:hover {{ text-decoration: underline; }}
        @media (max-width: 768px) {{ .ascii-art {{ font-size: 6px; }} .stats-grid {{ grid-template-columns: 1fr; }} }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="quote-container">
            <div class="quote-text">"{random_quote['text']}"</div>
            <div class="quote-author">&mdash; {random_quote['author']}</div>
        </div>
    </div>

    <div class="card">
        <h2>Batch AMR Summary</h2>
        <div class="stats-grid">
            <div class="stat-card"><div class="stat-value">{total_genomes}</div><div class="stat-label">Total Genomes</div></div>
            <div class="stat-card"><div class="stat-value">{total_hits}</div><div class="stat-label">Total AMR Genes</div></div>
            <div class="stat-card critical-card"><div class="stat-value">{genomes_with_critical}</div><div class="stat-label">Critical Risk Genomes</div></div>
            <div class="stat-card"><div class="stat-value">{total_mutations}</div><div class="stat-label">Total Mutations</div></div>
        </div>
        <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
        <p><strong>AMRfinderPlus:</strong> {self.metadata['amrfinder_version']}</p>
        <p><strong>Database:</strong> {self.metadata['database_version']}</p>
        <p><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
    </div>

    {critical_alert}
    {high_risk_alert}

    <div class="card">
        <h2>Gene Frequency</h2>
        <div class="table-responsive">
            <table class="data-table">
                <thead><tr><th>Gene</th><th>Frequency</th><th>Genomes</th></tr></thead>
                <tbody>{freq_rows}</tbody>
            </table>
        </div>
    </div>

    {self._get_footer_html("AcinetoScope - AMRfinderPlus Batch Summary")}
</div>
</body>
</html>"""

        out_file = os.path.join(output_base, "acineto_amrfinder_summary_report.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"AMR summary HTML: {out_file}")

    def _create_amr_master_json(self, all_results: Dict[str, Any], output_base: str):
        """Write the SINGLE batch JSON — AMR + mutation data combined."""
        gene_freq = defaultdict(lambda: {'count': 0, 'genomes': set()})
        for genome_name, result in all_results.items():
            for hit in result['hits']:
                gene = hit.get('gene_symbol', '')
                if gene:
                    gene_freq[gene]['count'] += 1
                    gene_freq[gene]['genomes'].add(genome_name)
        for gene in gene_freq:
            gene_freq[gene]['genomes'] = list(gene_freq[gene]['genomes'])

        # --- Mutation aggregation (merged from prior _create_mutation_master_json) ---
        all_mutations = []
        genome_mutation_counts = {}
        for genome_name, result in all_results.items():
            muts = result.get('mutations', []) or []
            genome_mutation_counts[genome_name] = len(muts)
            for m in muts:
                m_copy = m.copy()
                m_copy['genome'] = genome_name
                all_mutations.append(m_copy)

        gene_mut_freq = defaultdict(lambda: {'count': 0, 'genomes': set()})
        for m in all_mutations:
            gene = m.get('gene_symbol', 'unknown')
            mut = m.get('sequence_name', '')
            key = f"{gene}_{mut}"
            gene_mut_freq[key]['count'] += 1
            gene_mut_freq[key]['genomes'].add(m.get('genome', ''))
        for k in gene_mut_freq:
            gene_mut_freq[k]['genomes'] = list(gene_mut_freq[k]['genomes'])

        master = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results),
                'database_version': self.metadata['database_version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
            },
            'amr': {
                'gene_frequency': dict(gene_freq),
                'genome_summaries': {
                    genome: {
                        'total_hits': result['hit_count'],
                        'genes': list({hit['gene_symbol']
                                       for hit in result['hits'] if hit['gene_symbol']}),
                    } for genome, result in all_results.items()
                },
            },
            'mutations': {
                'genome_mutation_counts': genome_mutation_counts,
                'mutation_frequency': dict(gene_mut_freq),
                'all_mutations': all_mutations,
            },
        }
        json_file = os.path.join(output_base, "acineto_amrfinder_master_summary.json")
        with open(json_file, 'w') as f:
            json.dump(master, f, indent=2, default=str)
        self.logger.info(f"Master JSON (single): {json_file}")

    def create_mutation_summary(self, all_results: Dict[str, Any], output_base: str):
        """Write mutation batch TSV + HTML (JSON is merged into the master JSON)."""
        self.logger.info("Creating mutation batch summaries...")
        all_mutations = []
        genome_mutation_counts = {}
        for genome_name, result in all_results.items():
            muts = result.get('mutations', []) or []
            genome_mutation_counts[genome_name] = len(muts)
            for m in muts:
                m_copy = m.copy()
                m_copy['genome'] = genome_name
                all_mutations.append(m_copy)

        if not all_mutations:
            self.logger.info("No mutations found in any genome; skipping mutation summaries.")
            return

        tsv_file = os.path.join(output_base, "mutation_summary.tsv")
        with open(tsv_file, 'w') as f:
            f.write("genome\tgene_symbol\tsequence_name\tclass\tsubclass\tcontig_id\tstart\tstop\t"
                    "strand\tcoverage\tidentity\taccession\n")
            for m in all_mutations:
                row = [
                    m.get('genome', ''),
                    m.get('gene_symbol', ''),
                    m.get('sequence_name', ''),
                    m.get('class', ''),
                    m.get('subclass', ''),
                    m.get('contig_id', ''),
                    m.get('start', ''),
                    m.get('stop', ''),
                    m.get('strand', ''),
                    m.get('coverage', ''),
                    m.get('identity', ''),
                    m.get('accession', ''),
                ]
                f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info(f"Mutation TSV: {tsv_file}")

        self._create_mutation_summary_html(all_mutations, genome_mutation_counts, output_base)

    def _create_mutation_summary_html(self, all_mutations: List[Dict],
                                      genome_counts: Dict[str, int], output_base: str):
        gene_freq = {}
        for m in all_mutations:
            gene = m.get('gene_symbol', 'unknown')
            mut = m.get('sequence_name', '')
            key = f"{gene}_{mut}"
            if key not in gene_freq:
                gene_freq[key] = {
                    'count': 0, 'genomes': set(),
                    'gene': gene, 'mutation': mut,
                    'class': m.get('class', ''),
                    'subclass': m.get('subclass', ''),
                }
            gene_freq[key]['count'] += 1
            gene_freq[key]['genomes'].add(m.get('genome', ''))
        for k in gene_freq:
            gene_freq[k]['genomes'] = ', '.join(sorted(gene_freq[k]['genomes']))
        sorted_freq = sorted(gene_freq.values(), key=lambda x: x['count'], reverse=True)

        random_quote = random.choice(self.science_quotes)
        total_genomes_with_mutations = sum(1 for c in genome_counts.values() if c > 0)

        rows = ""
        for item in sorted_freq:
            rows += f"""
                    <tr>
                        <td><strong>{item['gene']}</strong></td>
                        <td>{item['mutation']}</td>
                        <td>{item['count']}</td>
                        <td>{item['genomes']}</td>
                        <td>{item['class']}</td>
                        <td>{item['subclass']}</td>
                    </tr>
            """

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AcinetoScope - Mutation Batch Summary</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{ background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%); font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; color: #ffffff; padding: 20px; min-height: 100vh; }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{ background: rgba(0,0,0,0.7); padding: 20px; border-radius: 15px; margin-bottom: 20px; box-shadow: 0 8px 32px rgba(0,0,0,0.4); border: 2px solid rgba(0,255,0,0.3); }}
        .ascii-art {{ font-family: 'Courier New', monospace; font-size: 10px; line-height: 1.1; white-space: pre; color: #00ff00; text-shadow: 0 0 10px rgba(0,255,0,0.5); overflow-x: auto; }}
        .quote-container {{ background: rgba(255,255,255,0.1); backdrop-filter: blur(10px); padding: 20px; border-radius: 10px; margin-bottom: 30px; text-align: center; min-height: 100px; display: flex; flex-direction: column; justify-content: center; box-shadow: 0 4px 20px rgba(0,0,0,0.3); border: 1px solid rgba(255,255,255,0.2); }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; color: #ffffff; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .card {{ background: rgba(255,255,255,0.95); color: #1f2937; padding: 25px; border-radius: 10px; margin-bottom: 20px; box-shadow: 0 4px 15px rgba(0,0,0,0.2); }}
        .card h2 {{ color: #1e3a8a; border-bottom: 3px solid #3b82f6; padding-bottom: 10px; margin-bottom: 20px; font-size: 24px; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-bottom: 20px; }}
        .stat-card {{ background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%); color: white; padding: 15px; border-radius: 8px; text-align: center; }}
        .stat-value {{ font-size: 24px; font-weight: bold; }}
        .stat-label {{ font-size: 12px; opacity: 0.9; }}
        .table-responsive {{ width: 100%; overflow-x: auto; margin: 20px 0; }}
        .data-table {{ width: 100%; border-collapse: collapse; font-size: 13px; min-width: 900px; }}
        .data-table th {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 12px; text-align: left; font-weight: bold; }}
        .data-table td {{ padding: 10px; border-bottom: 1px solid #e5e7eb; }}
        .data-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .footer {{ text-align: center; margin-top: 30px; padding: 20px; background: rgba(0,0,0,0.3); border-radius: 10px; font-size: 14px; }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .authorship {{ margin-top: 15px; padding: 15px; background: rgba(255,255,255,0.1); border-radius: 8px; font-size: 12px; }}
        .authorship a {{ color: #fbbf24; text-decoration: none; }}
        .authorship a:hover {{ text-decoration: underline; }}
        @media (max-width: 768px) {{ .ascii-art {{ font-size: 6px; }} .stats-grid {{ grid-template-columns: 1fr; }} }}
    
        /* ════════ AcinetoScope uniformity override (auto-added) ════════ */
        .ascii-container {{ text-align: center !important; }}
        .ascii-art {{ display: inline-block !important; text-align: left !important; }}
        .footer {{ text-align: center !important; margin-top: 30px !important; padding: 20px !important; background: rgba(0, 0, 0, 0.3) !important; border-radius: 10px !important; font-size: 14px !important; overflow-wrap: break-word !important; }}
        .footer a {{ color: #fbbf24 !important; text-decoration: none !important; }}
        .footer a:hover {{ text-decoration: underline !important; }}
        .footer .timestamp {{ color: #fbbf24 !important; font-weight: bold !important; }}
        .footer .authorship {{ margin-top: 15px !important; padding: 15px !important; background: rgba(255, 255, 255, 0.1) !important; border-radius: 8px !important; font-size: 12px !important; line-height: 1.7 !important; }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="quote-container">
            <div class="quote-text">"{random_quote['text']}"</div>
            <div class="quote-author">&mdash; {random_quote['author']}</div>
        </div>
    </div>

    <div class="card">
        <h2>Mutation Summary Across All Genomes</h2>
        <div class="stats-grid">
            <div class="stat-card"><div class="stat-value">{len(all_mutations)}</div><div class="stat-label">Total Mutation Events</div></div>
            <div class="stat-card"><div class="stat-value">{total_genomes_with_mutations}</div><div class="stat-label">Genomes With Mutations</div></div>
            <div class="stat-card"><div class="stat-value">{len(sorted_freq)}</div><div class="stat-label">Unique Mutations</div></div>
            <div class="stat-card"><div class="stat-value">{datetime.now().strftime('%Y-%m-%d')}</div><div class="stat-label">Analysis Date</div></div>
        </div>
    </div>

    <div class="card">
        <h2>Mutation Frequency by Gene / Mutation</h2>
        <div class="table-responsive">
            <table class="data-table">
                <thead>
                    <tr><th>Gene</th><th>Mutation</th><th>Count</th><th>Genomes</th><th>Class</th><th>Subclass</th></tr>
                </thead>
                <tbody>{rows}</tbody>
            </table>
        </div>
    </div>

    {self._get_footer_html("AcinetoScope - Mutation Batch Summary")}
</div>
</body>
</html>"""

        out_file = os.path.join(output_base, "mutation_summary.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"Mutation summary HTML: {out_file}")


def main():
    parser = argparse.ArgumentParser(
        description='AcinetoScope AMRfinderPlus - A. baumannii AMR Analysis (BUNDLED VERSION)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python acineto_amrfinder.py "*.fna"
  python acineto_amrfinder.py "*.fna" --min-identity 0.95 --min-coverage 0.9
  python acineto_amrfinder.py "*.fna" --skip-mutations
  python acineto_amrfinder.py --force-update
  python acineto_amrfinder.py --db-version

PER-GENOME OUTPUTS (2 files each):
  - {genome}_amrfinder_report.html   (AMR genes + point mutations combined)
  - {genome}_amrfinder_raw.tsv       (Source column: AMR / MUTATION)

BATCH-LEVEL OUTPUTS:
  - acineto_amrfinder_summary.tsv
  - acineto_amrfinder_statistics_summary.tsv
  - acineto_amrfinder_summary_report.html
  - acineto_amrfinder_master_summary.json   (SINGLE batch JSON: AMR + mutations)
  - mutation_summary.tsv
  - mutation_summary.html
        """
    )

    parser.add_argument('pattern', nargs='?',
                        help='File pattern for A. baumannii genomes (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--cpus', '-c', type=int, default=None,
                        help='Number of CPU cores to use (default: auto-detect optimal)')
    parser.add_argument('--output', '-o', default='acineto_amrfinder_results',
                        help='Output directory (default: acineto_amrfinder_results)')
    parser.add_argument('--min-identity', type=float, default=None,
                        help='Minimum identity (0..1) for hits. Default: AMRfinder auto threshold')
    parser.add_argument('--min-coverage', type=float, default=None,
                        help='Minimum coverage of reference (0..1). Default: 0.5')
    parser.add_argument('--skip-mutations', action='store_true',
                        help='Skip point mutation reporting (mutations are reported by default)')
    parser.add_argument('--update-db', action='store_true',
                        help='Update AMRfinderPlus database to latest version (incremental) and exit')
    parser.add_argument('--force-update', action='store_true',
                        help='Force complete database update (overwrites existing folders) and exit')
    parser.add_argument('--db-version', action='store_true',
                        help='Show current database version and exit')

    args = parser.parse_args()

    executor = AcinetoAMRfinderPlus(cpus=args.cpus)

    if args.update_db or args.force_update or args.db_version:
        if args.force_update:
            print("Forcing full database update (this will overwrite existing database)...")
            success = executor.update_database(force=True)
            sys.exit(0 if success else 1)
        if args.update_db:
            print("Updating AMRfinderPlus database (incremental)...")
            success = executor.update_database(force=False)
            sys.exit(0 if success else 1)
        if args.db_version:
            print(f"Database version: {executor.metadata['database_version']}")
            print(f"Database path: {executor.bundled_database or 'Not found'}")
            sys.exit(0)

    if not args.pattern:
        parser.error("Please provide a file pattern for genomes (or use --update-db / --force-update / --db-version)")

    try:
        results = executor.process_multiple_genomes(
            args.pattern,
            args.output,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            report_mutations=not args.skip_mutations,
        )

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("AcinetoScope AMRfinderPlus FINAL SUMMARY")
        executor.logger.info("=" * 50)

        total_hits = sum(r['hit_count'] for r in results.values())
        total_mutations = sum(len(r.get('mutations', [])) for r in results.values())
        high_risk_count = 0
        critical_risk_count = 0
        for genome_name, result in results.items():
            genes = [hit['gene_symbol'] for hit in result['hits'] if hit['gene_symbol']]
            high_risk_count += sum(1 for g in genes if g in executor.high_risk_genes)
            critical_risk_count += sum(1 for g in genes if g in executor.critical_risk_genes)
            executor.logger.info(f"{genome_name}: {result['hit_count']} AMR hits, "
                                 f"{len(result.get('mutations', []))} mutations")

        executor.logger.info("\nA. BAUMANNII SUMMARY STATISTICS:")
        executor.logger.info(f"   Total genomes processed: {len(results)}")
        executor.logger.info(f"   Total AMR hits: {total_hits}")
        executor.logger.info(f"   Total point mutations: {total_mutations}")
        executor.logger.info(f"   High-risk genes detected: {high_risk_count}")
        executor.logger.info(f"   CRITICAL RISK genes detected: {critical_risk_count}")
        executor.logger.info(f"   Average AMR hits per genome: {total_hits / len(results) if results else 0:.1f}")

        executor.logger.info("\nPER-GENOME OUTPUTS (2 files each):")
        executor.logger.info("   - {genome}_amrfinder_report.html   (AMR + mutations combined)")
        executor.logger.info("   - {genome}_amrfinder_raw.tsv       (Source column: AMR / MUTATION)")

        executor.logger.info("\nBATCH-LEVEL OUTPUTS:")
        executor.logger.info(f"   AMR summary TSV:     {args.output}/acineto_amrfinder_summary.tsv")
        executor.logger.info(f"   Statistics TSV:      {args.output}/acineto_amrfinder_statistics_summary.tsv")
        executor.logger.info(f"   Summary HTML:        {args.output}/acineto_amrfinder_summary_report.html")
        executor.logger.info(f"   Master JSON (1):     {args.output}/acineto_amrfinder_master_summary.json")
        executor.logger.info(f"   Mutation TSV:        {args.output}/mutation_summary.tsv")
        executor.logger.info(f"   Mutation HTML:       {args.output}/mutation_summary.html")

        executor.logger.info("\nPERFORMANCE SUMMARY:")
        executor.logger.info(f"   CPU cores utilized: {executor.cpus}")
        executor.logger.info(f"   Available RAM: {executor.available_ram:.1f} GB")
        executor.logger.info(f"   Processing mode: MAXIMUM SPEED CONCURRENT MODE")
        executor.logger.info(f"   Bundled AMRfinderPlus: {executor.metadata['amrfinder_version']}")
        executor.logger.info(f"   Bundled database: {executor.metadata['database_version']}")

        if critical_risk_count > 0:
            executor.logger.info("\nCRITICAL RISK ALERT: last-resort antibiotic resistance genes detected.")
            executor.logger.info("   Immediate clinical attention and infection control measures required.")

        quote = random.choice(executor.science_quotes)
        executor.logger.info(f"\n{quote['text']} - {quote['author']}")

    except Exception as e:
        executor.logger.error(f"A. baumannii AMR analysis failed: {e}")
        import traceback
        executor.logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()