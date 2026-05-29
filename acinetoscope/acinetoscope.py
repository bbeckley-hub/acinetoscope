#!/usr/bin/env python3
"""
AcinetoScope Main Orchestrator - Colored Sequential Execution with Scientific Quotes
Complete A. baumannii typing pipeline - MLST, K/O, AMR, Plasmid, Virulence, QC, Summary
Author: Brown Beckley <brownbeckley94@gmail.com>
Date: 2026-05-29
Version: 1.2.0 - Added APT plasmid typing module and improved error handling
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import random
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set


# ----------------------------------------------------------------------
# Version and color definitions
# ----------------------------------------------------------------------
__version__ = "1.2.0"

class Color:
    """ANSI color codes for colored output"""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    
    BRIGHT_BLACK = '\033[90m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'

class AcinetoScopeOrchestrator:
    """AcinetoScope orchestrator with colored sequential execution and scientific quotes"""
    
    def __init__(self):
        self.base_dir = Path(__file__).parent
        self.setup_colors()
        self.quotes = self._get_scientific_quotes()
        self.quote_colors = [
            Color.BRIGHT_CYAN, Color.BRIGHT_GREEN, Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA, Color.BRIGHT_BLUE, Color.BRIGHT_RED,
            Color.CYAN, Color.GREEN, Color.YELLOW, Color.MAGENTA
        ]
        
        # Output directory names (must match what modules produce)
        self.output_dirs = {
            'qc': 'fasta_qc_results',
            'abricate': 'acineto_abricate_results',
            'amr': 'acineto_amrfinder_results',
            'kaptive': 'kaptive_results',
            'mlst_pasteur': 'PASTEUR_MLST',
            'mlst_oxford': 'OXFORD_MLST',
            'plasmid': 'acineto_plasmid_results'        # new
        }
        
        # HTML files required for the summary reporter (including plasmid)
        self.summary_html_files = {
            # FASTA QC
            'FASTA_QC_summary.html': 'fasta_qc_results/FASTA_QC_summary.html',
            # MLST
            'pasteur_mlst_summary.html': 'pasteur_mlst_summary.html',
            'oxford_mlst_summary.html': 'oxford_mlst_summary.html',
            # Kaptive
            'Kaptive_summary.html': 'kaptive_results/Kaptive_summary.html',
            # AMR
            'acineto_amrfinder_summary_report.html': 'acineto_amrfinder_results/acineto_amrfinder_summary_report.html',
            # ABRicate databases
            'acineto_card_summary_report.html': 'acineto_abricate_results/acineto_card_summary_report.html',
            'acineto_ncbi_summary_report.html': 'acineto_abricate_results/acineto_ncbi_summary_report.html',
            'acineto_resfinder_summary_report.html': 'acineto_abricate_results/acineto_resfinder_summary_report.html',
            'acineto_vfdb_summary_report.html': 'acineto_abricate_results/acineto_vfdb_summary_report.html',
            'acineto_argannot_summary_report.html': 'acineto_abricate_results/acineto_argannot_summary_report.html',
            'acineto_megares_summary_report.html': 'acineto_abricate_results/acineto_megares_summary_report.html',
            'acineto_ecoli_vf_summary_report.html': 'acineto_abricate_results/acineto_ecoli_vf_summary_report.html',
            'acineto_bacmet2_summary_report.html': 'acineto_abricate_results/acineto_bacmet2_summary_report.html',
            'acineto_ecoh_summary_report.html': 'acineto_abricate_results/acineto_ecoh_summary_report.html',
            'acineto_plasmidfinder_summary_report.html': 'acineto_abricate_results/acineto_plasmidfinder_summary_report.html',
            # New plasmid typing (APT)
            'plasmid_summary_report.html': 'acineto_plasmid_results/plasmid_summary_report.html'
        }
    
    def _get_scientific_quotes(self):
        """Curated scientific quotes (same as original)"""
        return [
            {"quote": "Science is organized knowledge.", "author": "Herbert Spencer", "theme": "knowledge"},
            {"quote": "The science of today is the technology of tomorrow.", "author": "Edward Teller", "theme": "technology"},
            {"quote": "Nature is the source of all true knowledge.", "author": "Leonardo da Vinci", "theme": "nature"},
            {"quote": "Biology is the most powerful technology ever created.", "author": "Freeman Dyson", "theme": "biology"},
            {"quote": "Genomics is a lens on biology.", "author": "Eric Lander", "theme": "genomics"},
            {"quote": "Every microbe has its own story.", "author": "Anonymous", "theme": "microbiology"},
            {"quote": "Data beats emotions.", "author": "Sean Rad", "theme": "data"},
            {"quote": "In every walk with nature, one receives far more than he seeks.", "author": "John Muir", "theme": "discovery"},
            {"quote": "The microbe is nothing; the terrain is everything.", "author": "Louis Pasteur", "theme": "microbiology"},
            {"quote": "What we know is a drop, what we don't know is an ocean.", "author": "Isaac Newton", "theme": "knowledge"},
            {"quote": "The good physician treats the disease; the great physician treats the patient who has the disease.", "author": "William Osler", "theme": "medicine"},
            {"quote": "In science, there are no shortcuts to truth.", "author": "Karl Popper", "theme": "science"},
            {"quote": "The art of research is the art of making difficult problems soluble.", "author": "Peter Medawar", "theme": "research"},
            {"quote": "Equipped with his five senses, man explores the universe around him and calls the adventure Science.", "author": "Edwin Hubble", "theme": "exploration"},
            {"quote": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein", "theme": "curiosity"},
            {"quote": "One must learn by doing the thing; though you think you know it, you have no certainty until you try.", "author": "Sophocles", "theme": "practice"},
            {"quote": "The secret of getting ahead is getting started.", "author": "Mark Twain", "theme": "motivation"},
            {"quote": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie", "theme": "understanding"},
            {"quote": "The greatest enemy of knowledge is not ignorance, it is the illusion of knowledge.", "author": "Stephen Hawking", "theme": "knowledge"},
            {"quote": "We are just an advanced breed of monkeys on a minor planet of a very average star. But we can understand the Universe.", "author": "Stephen Hawking", "theme": "perspective"},
            {"quote": "To raise new questions, new possibilities, to regard old problems from a new angle, requires creative imagination.", "author": "Albert Einstein", "theme": "innovation"},
            {"quote": "The most exciting phrase to hear in science, the one that heralds new discoveries, is not 'Eureka!' but 'That's funny...'", "author": "Isaac Asimov", "theme": "discovery"},
            {"quote": "In science the credit goes to the man who convinces the world, not to the man to whom the idea first occurs.", "author": "Francis Darwin", "theme": "recognition"},
            {"quote": "The aim of science is not to open the door to infinite wisdom, but to set a limit to infinite error.", "author": "Bertolt Brecht", "theme": "purpose"},
            {"quote": "Science knows no country, because knowledge belongs to humanity, and is the torch which illuminates the world.", "author": "Louis Pasteur", "theme": "global"},
            {"quote": "The first rule of intelligent tinkering is to save all the parts.", "author": "Paul Ehrlich", "theme": "conservation"},
            {"quote": "DNA is like a computer program but far, far more advanced than any software ever created.", "author": "Bill Gates", "theme": "genomics"},
            {"quote": "The beauty of a living thing is not the atoms that go into it, but the way those atoms are put together.", "author": "Carl Sagan", "theme": "biology"},
            {"quote": "If I have seen further it is by standing on the shoulders of Giants.", "author": "Isaac Newton", "theme": "collaboration"}
        ]
    
    def display_random_quote(self):
        if not self.quotes:
            return
        quote_data = random.choice(self.quotes)
        quote = quote_data["quote"]
        author = quote_data["author"]
        theme = quote_data.get("theme", "science")
        quote_color = random.choice(self.quote_colors)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        theme_icons = {"microbiology": "🦠", "discovery": "🔬", "knowledge": "📚", "medicine": "⚕️",
                       "science": "🧪", "research": "🔍", "exploration": "🚀", "curiosity": "🤔",
                       "practice": "🛠️", "motivation": "💪", "nature": "🌿", "inquiry": "❓",
                       "technology": "💻", "understanding": "🧠", "perspective": "👁️", "innovation": "💡",
                       "recognition": "🏆", "purpose": "🎯", "biology": "🧬", "genomics": "🧬",
                       "data": "📊", "global": "🌍", "conservation": "🌱", "collaboration": "🤝"}
        icon = theme_icons.get(theme, "💭")
        print()
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}[{current_time}] {icon} SCIENTIFIC INSIGHT: {Color.RESET}")
        print()
        print(f"{quote_color}   \"{quote}\"{Color.RESET}")
        print(f"{Color.BOLD}{Color.WHITE}   — {author}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}   Theme: {theme.capitalize()}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print()
    
    def setup_colors(self):
        self.color_info = Color.CYAN
        self.color_success = Color.BRIGHT_GREEN
        self.color_warning = Color.BRIGHT_YELLOW
        self.color_error = Color.BRIGHT_RED
        self.color_highlight = Color.BRIGHT_CYAN
        self.color_banner = Color.BRIGHT_MAGENTA
        self.color_module = Color.BRIGHT_BLUE
        self.color_sample = Color.GREEN
        self.color_file = Color.YELLOW
        self.color_reset = Color.RESET
    
    def print_color(self, message: str, color: str = Color.RESET, bold: bool = False):
        style = Color.BOLD if bold else ''
        print(f"{style}{color}{message}{Color.RESET}")
    
    def print_header(self, title: str, subtitle: str = ""):
        print()
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_CYAN}{' ' * 20}{title}{Color.RESET}")
        if subtitle:
            print(f"{Color.DIM}{Color.WHITE}{' ' * 22}{subtitle}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print()
    
    def print_info(self, message: str):
        print(f"{self.color_info}[INFO]{Color.RESET} {message}")
    
    def print_success(self, message: str):
        print(f"{self.color_success}✓{Color.RESET} {message}")
    
    def print_warning(self, message: str):
        print(f"{self.color_warning}⚠️{Color.RESET} {message}")
    
    def print_error(self, message: str):
        print(f"{self.color_error}✗{Color.RESET} {message}")
    
    def print_command(self, command: str):
        print(f"{Color.DIM}{Color.WHITE}  $ {command}{Color.RESET}")
    
    def display_banner(self):
        banner = f"""{Color.BOLD}{Color.BRIGHT_MAGENTA}
{'='*80}
{' '*20}🦠 ACINETOSCOPE v{__version__} - A. baumannii Genomic Analysis Pipeline
{'='*80}
{Color.RESET}{Color.BRIGHT_CYAN}
Complete A. baumannii genomic analysis pipeline
MLST | K/O Locus | AMR | Virulence | Plasmid (APT) | Quality Control | Critical Genes Flagging | Summary Reports
{Color.RESET}{Color.DIM}
Author: Brown Beckley | Email: brownbeckley94@gmail.com
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
{'='*80}{Color.RESET}
"""
        print(banner)
    
    def find_fasta_files(self, input_path: str) -> List[Path]:
        self.print_info(f"Searching for files with pattern: {input_path}")
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files if Path(f).is_file() and 
                          f.lower().endswith(('.fna', '.fasta', '.fa', '.fn')) and
                          not Path(f).name.startswith('.')]
            self.print_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)
        input_path_obj = Path(input_path)
        if input_path_obj.is_file() and input_path_obj.suffix.lower() in ['.fna', '.fasta', '.fa', '.fn']:
            self.print_success(f"Found single FASTA file: {input_path_obj.name}")
            return [input_path_obj]
        if input_path_obj.is_dir():
            patterns = [f"{input_path}/*.fna", f"{input_path}/*.fasta", f"{input_path}/*.fa", f"{input_path}/*.fn"]
            fasta_files = []
            for pattern in patterns:
                for file_path in glob.glob(pattern):
                    if Path(file_path).is_file() and not Path(file_path).name.startswith('.'):
                        fasta_files.append(Path(file_path))
            fasta_files = sorted(list(set(fasta_files)))
            if fasta_files:
                self.print_success(f"Found {len(fasta_files)} FASTA files in directory")
            else:
                self.print_warning(f"No FASTA files found in directory: {input_path}")
            return fasta_files
        self.print_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        if not fasta_files:
            return '"*.fna"'
        extensions = set(f.suffix.lower() for f in fasta_files)
        if len(extensions) == 1:
            ext = list(extensions)[0]
            return f'"*{ext}"'
        return '"*"'

    def cleanup_module_directory(self, module_path: Path, fasta_files: List[Path]):
        try:
            self.print_info(f"Cleaning up {module_path.name}...")
            for fasta_file in fasta_files:
                temp_file = module_path / fasta_file.name
                if temp_file.exists():
                    temp_file.unlink()
            exact_dirs_to_remove = ["fasta_qc_results", "acineto_abricate_results", "acineto_amrfinder_results",
                                    "kaptive_results", "PASTEUR_MLST", "OXFORD_MLST", "mlst_pasteur_results",
                                    "mlst_oxford_results", "results", "acineto_plasmid_results"]
            for dir_name in exact_dirs_to_remove:
                dir_path = module_path / dir_name
                if dir_path.exists():
                    shutil.rmtree(dir_path)
            for html_file in self.summary_html_files.keys():
                html_path = module_path / html_file
                if html_path.exists():
                    html_path.unlink()
            for html_file in module_path.glob("*.html"):
                if html_file.is_file():
                    html_file.unlink()
            self.print_success(f"✅ {module_path.name} cleaned up successfully")
        except Exception as e:
            self.print_warning(f"⚠️  Partial cleanup issue in {module_path.name}: {str(e)}")

    # ----------------------------------------------------------------------
    # Module execution methods (each returns True on completion, False on fatal error)
    # ----------------------------------------------------------------------
    def run_qc_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        qc_module_path = self.base_dir / "modules" / "qc_module"
        try:
            self.print_header("FASTA QC ANALYSIS", "Comprehensive Quality Control")
            qc_script = qc_module_path / "acineto_fasta_qc.py"
            if not qc_script.exists():
                self.print_error(f"QC script not found at: {qc_script}")
                return False
            for fasta_file in fasta_files:
                shutil.copy2(fasta_file, qc_module_path / fasta_file.name)
            self.print_info(f"Copied {len(fasta_files)} files to QC module")
            file_pattern = self.get_file_pattern(fasta_files)
            pattern_clean = file_pattern.strip('"')
            cmd = [sys.executable, str(qc_script), pattern_clean]
            self.print_info(f"Running QC analysis with pattern: {file_pattern}")
            self.print_command(f"python3 {qc_script.name} {pattern_clean}")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=qc_module_path)
            if result.returncode == 0:
                self.print_success("QC analysis completed!")
                qc_source = qc_module_path / "fasta_qc_results"
                qc_target = output_dir / "fasta_qc_results"
                if qc_source.exists():
                    if qc_target.exists():
                        shutil.rmtree(qc_target)
                    shutil.copytree(qc_source, qc_target)
                    self.print_success(f"QC results copied to: {qc_target}")
                else:
                    self.print_warning("QC results directory not found: fasta_qc_results")
                self.display_random_quote()
                return True
            else:
                self.print_warning("QC analysis had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:5]:
                        if "error" in line.lower() or "failed" in line.lower():
                            print(f"{self.color_warning}  {line}{Color.RESET}")
                return True
        except Exception as e:
            self.print_error(f"QC analysis failed: {str(e)}")
            return False
        finally:
            self.cleanup_module_directory(qc_module_path, fasta_files)

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int, scheme: str) -> bool:
        mlst_module_path = self.base_dir / "modules" / "mlst_module"
        try:
            scheme_name = "PASTEUR" if scheme == "pasteur" else "OXFORD"
            self.print_header(f"MLST ANALYSIS - {scheme_name}", "Multi-Locus Sequence Typing")
            mlst_script = mlst_module_path / "mlst_module.py"
            if not mlst_script.exists():
                self.print_error(f"MLST script not found at: {mlst_script}")
                return False
            for fasta_file in fasta_files:
                shutil.copy2(fasta_file, mlst_module_path / fasta_file.name)
            self.print_info(f"Copied {len(fasta_files)} files to MLST module")
            file_pattern = self.get_file_pattern(fasta_files)
            pattern_clean = file_pattern.strip('"')
            output_subdir = f"mlst_{scheme}_results"
            cmd = [sys.executable, str(mlst_script), "-i", pattern_clean, "-o", output_subdir, "-db", "db", "-sc", "bin", "--batch", "-s", scheme]
            self.print_info(f"Running MLST analysis with scheme: {scheme_name}")
            self.print_command(f"python {mlst_script.name} -i {pattern_clean} -o {output_subdir} -db db -sc bin --batch -s {scheme}")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=mlst_module_path)
            if result.returncode == 0:
                self.print_success(f"MLST analysis completed for {scheme_name} scheme!")
                scheme_dir = "PASTEUR_MLST" if scheme == "pasteur" else "OXFORD_MLST"
                mlst_source = mlst_module_path / output_subdir / scheme_dir
                mlst_target = output_dir / scheme_dir
                if mlst_source.exists():
                    if mlst_target.exists():
                        shutil.rmtree(mlst_target)
                    shutil.copytree(mlst_source, mlst_target)
                    self.print_success(f"MLST results copied to: {mlst_target}")
                else:
                    self.print_warning(f"MLST results directory not found: {mlst_source}")
                html_filename = f"{scheme}_mlst_summary.html"
                html_source = mlst_source / html_filename
                if html_source.exists():
                    html_target = output_dir / html_filename
                    shutil.copy2(html_source, html_target)
                    self.print_success(f"MLST HTML report copied to: {html_target}")
                else:
                    self.print_warning(f"MLST HTML report not found: {html_filename}")
                self.display_random_quote()
                return True
            else:
                self.print_warning(f"MLST analysis for {scheme_name} had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:5]:
                        if "error" in line.lower() or "failed" in line.lower():
                            print(f"{self.color_warning}  {line}{Color.RESET}")
                return True
        except Exception as e:
            self.print_error(f"MLST analysis failed for {scheme}: {str(e)}")
            return False
        finally:
            self.cleanup_module_directory(mlst_module_path, fasta_files)

    def run_kaptive_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        kaptive_module_path = self.base_dir / "modules" / "k_o_module"
        try:
            self.print_header("KAPTIVE ANALYSIS", "K and O Locus Typing")
            kaptive_script = kaptive_module_path / "acineto_kaptive.py"
            if not kaptive_script.exists():
                self.print_error(f"Kaptive script not found at: {kaptive_script}")
                return False
            for fasta_file in fasta_files:
                shutil.copy2(fasta_file, kaptive_module_path / fasta_file.name)
            self.print_info(f"Copied {len(fasta_files)} files to Kaptive module")
            file_pattern = self.get_file_pattern(fasta_files)
            pattern_clean = file_pattern.strip('"')
            cmd = [sys.executable, str(kaptive_script), pattern_clean]
            self.print_info(f"Running Kaptive analysis with pattern: {file_pattern}")
            self.print_command(f"python3 {kaptive_script.name} {pattern_clean}")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=kaptive_module_path)
            if result.returncode == 0:
                self.print_success("Kaptive analysis completed!")
                kaptive_source = kaptive_module_path / "kaptive_results"
                kaptive_target = output_dir / "kaptive_results"
                if kaptive_source.exists():
                    if kaptive_target.exists():
                        shutil.rmtree(kaptive_target)
                    shutil.copytree(kaptive_source, kaptive_target)
                    self.print_success(f"Kaptive results copied to: {kaptive_target}")
                else:
                    self.print_warning("Kaptive results directory not found: kaptive_results")
                html_source = kaptive_module_path / "kaptive_results" / "Kaptive_summary.html"
                if html_source.exists():
                    html_target = output_dir / "Kaptive_summary.html"
                    shutil.copy2(html_source, html_target)
                    self.print_success(f"Kaptive HTML report copied to: {html_target}")
                else:
                    self.print_warning("Kaptive HTML report not found: Kaptive_summary.html")
                self.display_random_quote()
                return True
            else:
                self.print_warning("Kaptive analysis had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:5]:
                        if "error" in line.lower() or "failed" in line.lower():
                            print(f"{self.color_warning}  {line}{Color.RESET}")
                return True
        except Exception as e:
            self.print_error(f"Kaptive analysis failed: {str(e)}")
            return False
        finally:
            self.cleanup_module_directory(kaptive_module_path, fasta_files)

    def run_amr_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        amr_module_path = self.base_dir / "modules" / "amr_module"
        try:
            self.print_header("AMR ANALYSIS", "Antimicrobial Resistance Gene Detection")
            amr_script = amr_module_path / "acineto_amrfinder.py"
            if not amr_script.exists():
                self.print_error(f"AMR script not found at: {amr_script}")
                return False
            for fasta_file in fasta_files:
                shutil.copy2(fasta_file, amr_module_path / fasta_file.name)
            self.print_info(f"Copied {len(fasta_files)} files to AMR module")
            file_pattern = self.get_file_pattern(fasta_files)
            pattern_clean = file_pattern.strip('"')
            cmd = [sys.executable, str(amr_script), pattern_clean]
            self.print_info(f"Running AMR analysis with pattern: {file_pattern}")
            self.print_command(f"python3 {amr_script.name} {pattern_clean}")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
            if result.returncode == 0:
                self.print_success("AMR analysis completed!")
                amr_source = amr_module_path / "acineto_amrfinder_results"
                amr_target = output_dir / "acineto_amrfinder_results"
                if amr_source.exists():
                    if amr_target.exists():
                        shutil.rmtree(amr_target)
                    shutil.copytree(amr_source, amr_target)
                    self.print_success(f"AMR results copied to: {amr_target}")
                else:
                    self.print_warning("AMR results directory not found: acineto_amrfinder_results")
                html_source = amr_module_path / "acineto_amrfinder_results" / "acineto_amrfinder_summary_report.html"
                if html_source.exists():
                    html_target = output_dir / "acineto_amrfinder_summary_report.html"
                    shutil.copy2(html_source, html_target)
                    self.print_success(f"AMR HTML report copied to: {html_target}")
                else:
                    self.print_warning("AMR HTML report not found: acineto_amrfinder_summary_report.html")
                self.display_random_quote()
                return True
            else:
                self.print_warning("AMR analysis had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:5]:
                        if "error" in line.lower() or "failed" in line.lower():
                            print(f"{self.color_warning}  {line}{Color.RESET}")
                return True
        except Exception as e:
            self.print_error(f"AMR analysis failed: {str(e)}")
            return False
        finally:
            self.cleanup_module_directory(amr_module_path, fasta_files)

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        abricate_module_path = self.base_dir / "modules" / "abricate_module"
        try:
            self.print_header("ABRICATE ANALYSIS", "Comprehensive Resistance & Virulence Gene Screening & Plasmid Profiling")
            abricate_script = abricate_module_path / "acineto_abricate.py"
            if not abricate_script.exists():
                self.print_error(f"ABRicate script not found at: {abricate_script}")
                return False
            for fasta_file in fasta_files:
                shutil.copy2(fasta_file, abricate_module_path / fasta_file.name)
            self.print_info(f"Copied {len(fasta_files)} files to ABRicate module")
            file_pattern = self.get_file_pattern(fasta_files)
            pattern_clean = file_pattern.strip('"')
            cmd = [sys.executable, str(abricate_script), pattern_clean]
            self.print_info(f"Running ABRicate analysis with pattern: {file_pattern}")
            self.print_command(f"python3 {abricate_script.name} {pattern_clean}")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=abricate_module_path)
            if result.returncode == 0:
                self.print_success("ABRicate analysis completed!")
                abricate_source = abricate_module_path / "acineto_abricate_results"
                abricate_target = output_dir / "acineto_abricate_results"
                if abricate_source.exists():
                    if abricate_target.exists():
                        shutil.rmtree(abricate_target)
                    shutil.copytree(abricate_source, abricate_target)
                    self.print_success(f"ABRicate results copied to: {abricate_target}")
                else:
                    self.print_warning("ABRicate results directory not found: acineto_abricate_results")
                html_files_copied = 0
                for html_filename, relative_path in self.summary_html_files.items():
                    if html_filename.startswith('acineto_') and html_filename != 'acineto_amrfinder_summary_report.html':
                        html_source = abricate_module_path / "acineto_abricate_results" / html_filename
                        if html_source.exists():
                            html_target = output_dir / html_filename
                            shutil.copy2(html_source, html_target)
                            html_files_copied += 1
                            self.print_info(f"  ✓ Copied: {html_filename}")
                if html_files_copied > 0:
                    self.print_success(f"Copied {html_files_copied} ABRicate HTML reports")
                else:
                    self.print_warning("No ABRicate HTML reports found")
                self.display_random_quote()
                return True
            else:
                self.print_warning("ABRicate analysis had warnings")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:5]:
                        if "error" in line.lower() or "failed" in line.lower():
                            print(f"{self.color_warning}  {line}{Color.RESET}")
                return True
        except Exception as e:
            self.print_error(f"ABRicate analysis failed: {str(e)}")
            return False
        finally:
            self.cleanup_module_directory(abricate_module_path, fasta_files)

    def run_plasmid_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run APT‑based plasmid typing module."""
        plasmid_module_path = self.base_dir / "modules" / "plasmid_typing_module"
        try:
            self.print_header("PLASMID TYPING ANALYSIS", "Acinetobacter Plasmid Typing (APT)")
            plasmid_script = plasmid_module_path / "plasmid_typer.py"
            if not plasmid_script.exists():
                self.print_error(f"Plasmid script not found at: {plasmid_script}")
                return False
            
            # Copy FASTA files to the module directory
            for fasta_file in fasta_files:
                shutil.copy2(fasta_file, plasmid_module_path / fasta_file.name)
            self.print_info(f"Copied {len(fasta_files)} files to plasmid module")
            
            file_pattern = self.get_file_pattern(fasta_files)
            pattern_clean = file_pattern.strip('"')
            cmd = [sys.executable, str(plasmid_script), pattern_clean]
            self.print_info(f"Running APT plasmid typing with pattern: {file_pattern}")
            self.print_command(f"python3 {plasmid_script.name} {pattern_clean}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=plasmid_module_path)
            if result.returncode == 0:
                self.print_success("Plasmid typing completed!")
                plasmid_source = plasmid_module_path / "acineto_plasmid_results"
                plasmid_target = output_dir / "acineto_plasmid_results"
                if plasmid_source.exists():
                    if plasmid_target.exists():
                        shutil.rmtree(plasmid_target)
                    shutil.copytree(plasmid_source, plasmid_target)
                    self.print_success(f"Plasmid results copied to: {plasmid_target}")
                else:
                    self.print_warning("Plasmid results directory not found: acineto_plasmid_results")
                
                # Copy the plasmid summary HTML to the output root (for later summary module)
                html_source = plasmid_module_path / "acineto_plasmid_results" / "plasmid_summary_report.html"
                if html_source.exists():
                    html_target = output_dir / "plasmid_summary_report.html"
                    shutil.copy2(html_source, html_target)
                    self.print_success(f"Plasmid HTML report copied to: {html_target}")
                else:
                    self.print_warning("Plasmid HTML report not found: plasmid_summary_report.html")
                self.display_random_quote()
                return True
            else:
                self.print_warning("Plasmid typing had warnings (non‑fatal, continuing)")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines[:5]:
                        if "error" in line.lower() or "failed" in line.lower():
                            print(f"{self.color_warning}  {line}{Color.RESET}")
                return True   # non‑fatal, continue pipeline
        except Exception as e:
            self.print_error(f"Plasmid analysis failed: {str(e)}")
            return False
        finally:
            self.cleanup_module_directory(plasmid_module_path, fasta_files)

    def copy_files_to_summary_module(self, output_dir: Path) -> Dict[str, bool]:
        try:
            self.print_header("PREPARING SUMMARY MODULE", "Copying required HTML files")
            summary_module_path = self.base_dir / "modules" / "summary_module"
            required_files = {}
            copied_count = 0
            missing_count = 0
            for target_filename, relative_path in self.summary_html_files.items():
                source_path = output_dir / relative_path
                if source_path.exists():
                    target_path = summary_module_path / target_filename
                    shutil.copy2(source_path, target_path)
                    copied_count += 1
                    required_files[target_filename] = True
                    self.print_success(f"  ✓ {target_filename}")
                else:
                    required_files[target_filename] = False
                    missing_count += 1
                    self.print_warning(f"  ✗ {target_filename} (not found at: {relative_path})")
            self.print_info(f"Copied {copied_count} files, {missing_count} files missing")
            # Critical files that must exist for the reporter to run
            critical_files = ["pasteur_mlst_summary.html", "oxford_mlst_summary.html", "Kaptive_summary.html",
                              "acineto_amrfinder_summary_report.html", "acineto_card_summary_report.html"]
            missing_critical = [f for f in critical_files if not required_files.get(f, False)]
            if missing_critical:
                self.print_warning(f"Missing critical files: {', '.join(missing_critical)}")
                return {"success": False, "files": required_files}
            return {"success": True, "files": required_files}
        except Exception as e:
            self.print_error(f"Error copying files to summary module: {str(e)}")
            return {"success": False, "files": {}}

    def run_summary_analysis(self, output_dir: Path) -> bool:
        try:
            self.print_header("ULTIMATE REPORTER", "Gene-centric Integrated Analysis")
            summary_module_path = self.base_dir / "modules" / "summary_module"
            summary_script = summary_module_path / "genius_acinetobacter_reporter.py"
            if not summary_script.exists():
                self.print_error(f"Summary script not found at: {summary_script}")
                return False
            cmd = [sys.executable, str(summary_script), "-i", "."]
            self.print_info("Running ultimate reporter...")
            self.print_command(f"python3 {summary_script.name} -i .")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=summary_module_path)
            if result.returncode == 0:
                self.print_success("Ultimate reporter completed successfully!")
                ultimate_source = summary_module_path / "GENIUS_ACINETOBACTER_ULTIMATE_REPORTS"
                ultimate_target = output_dir / "GENIUS_ACINETOBACTER_ULTIMATE_REPORTS"
                if ultimate_source.exists():
                    if ultimate_target.exists():
                        shutil.rmtree(ultimate_target)
                    shutil.copytree(ultimate_source, ultimate_target)
                    self.print_success(f"Ultimate reports copied to: {ultimate_target}")
                else:
                    self.print_warning("Ultimate reports directory not found: GENIUS_ACINETOBACTER_ULTIMATE_REPORTS")
                self.display_random_quote()
                return True
            else:
                self.print_warning("Ultimate reporter had issues")
                if result.stderr:
                    error_lines = result.stderr.strip().split('\n')
                    for line in error_lines:
                        if "error" in line.lower() or "failed" in line.lower():
                            print(f"{self.color_warning}  {line}{Color.RESET}")
                return True
        except Exception as e:
            self.print_error(f"Ultimate reporter failed: {str(e)}")
            return False

    def run_sequential_analyses(self, fasta_files: List[Path], output_dir: Path, threads: int, 
                               skip_modules: Dict[str, bool], mlst_scheme: str = "both") -> Dict[str, bool]:
        analysis_functions = []
        if not skip_modules.get('qc', False):
            analysis_functions.append(("QC Analysis", self.run_qc_analysis, True))
        if not skip_modules.get('mlst', False):
            if mlst_scheme in ["pasteur", "both"]:
                analysis_functions.append(("MLST Pasteur", lambda f, o, t: self.run_mlst_analysis(f, o, t, "pasteur"), True))
            if mlst_scheme in ["oxford", "both"]:
                analysis_functions.append(("MLST Oxford", lambda f, o, t: self.run_mlst_analysis(f, o, t, "oxford"), True))
        if not skip_modules.get('kaptive', False):
            analysis_functions.append(("Kaptive Analysis", self.run_kaptive_analysis, True))
        if not skip_modules.get('amr', False):
            analysis_functions.append(("AMR Analysis", self.run_amr_analysis, True))
        if not skip_modules.get('abricate', False):
            analysis_functions.append(("ABRicate Analysis", self.run_abricate_analysis, True))
        if not skip_modules.get('plasmid', False):
            analysis_functions.append(("Plasmid Typing", self.run_plasmid_analysis, True))
        
        if not analysis_functions:
            self.print_warning("All analyses were skipped! Nothing to run.")
            return {}
        
        self.print_info(f"Running {len(analysis_functions)} analyses sequentially")
        results = {}
        for analysis_name, analysis_func, _ in analysis_functions:
            try:
                success = analysis_func(fasta_files, output_dir, max(1, threads))
                results[analysis_name] = success
                if success:
                    self.print_success(f"✅ {analysis_name} completed")
                else:
                    self.print_error(f"❌ {analysis_name} failed")
            except Exception as e:
                self.print_error(f"❌ {analysis_name} failed with exception: {str(e)}")
                results[analysis_name] = False
        return results

    def run_force_update_amr(self) -> bool:
        amr_module_path = self.base_dir / "modules" / "amr_module"
        amr_script = amr_module_path / "acineto_amrfinder.py"
        if not amr_script.exists():
            self.print_error(f"AMR script not found at {amr_script}")
            return False
        self.print_header("FORCE UPDATE - AMR DATABASE", "Overwriting local AMR database with latest version")
        self.print_info("Running AMR database force update (this may take several minutes)...")
        cmd = [sys.executable, str(amr_script), "--force-update"]
        self.print_command(f"python3 {amr_script.name} --force-update")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
            if result.returncode == 0:
                self.print_success("✅ AMR database successfully updated (force).")
                return True
            else:
                self.print_error("AMR database update failed.")
                if result.stderr:
                    self.print_error(result.stderr.strip())
                return False
        except Exception as e:
            self.print_error(f"Error during AMR update: {str(e)}")
            return False

    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1, 
                             skip_modules: Dict[str, bool] = None, mlst_scheme: str = "both",
                             skip_summary: bool = False):
        if skip_modules is None:
            skip_modules = {}
        start_time = datetime.now()
        try:
            self.display_banner()
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            fasta_files = self.find_fasta_files(input_path)
            if not fasta_files:
                self.print_error("No FASTA files found! Analysis stopped.")
                return
            extensions = set(f.suffix.lower() for f in fasta_files)
            self.print_success(f"Starting analysis of {len(fasta_files)} A. baumannii samples")
            self.print_info(f"File formats detected: {', '.join(extensions)}")
            subdirs = ["fasta_qc_results", "PASTEUR_MLST", "OXFORD_MLST", "kaptive_results",
                       "acineto_amrfinder_results", "acineto_abricate_results", "acineto_plasmid_results",
                       "GENIUS_ACINETOBACTER_ULTIMATE_REPORTS"]
            for subdir in subdirs:
                (output_path / subdir).mkdir(exist_ok=True)
            self.print_header("ANALYSIS PLAN", "Modules to be executed")
            analyses_to_run = [
                ("QC Analysis", not skip_modules.get('qc', False)),
                ("MLST Pasteur", not skip_modules.get('mlst', False) and mlst_scheme in ["pasteur", "both"]),
                ("MLST Oxford", not skip_modules.get('mlst', False) and mlst_scheme in ["oxford", "both"]),
                ("Kaptive Analysis", not skip_modules.get('kaptive', False)),
                ("AMR Analysis", not skip_modules.get('amr', False)),
                ("ABRicate Analysis", not skip_modules.get('abricate', False)),
                ("Plasmid Typing", not skip_modules.get('plasmid', False)),
                ("Ultimate Reporter", not skip_summary),
            ]
            for analysis, enabled in analyses_to_run:
                if enabled:
                    print(f"   {Color.BRIGHT_GREEN}✅ ENABLED{Color.RESET} - {analysis}")
                else:
                    print(f"   {Color.YELLOW}⏸️  SKIPPED{Color.RESET} - {analysis}")
            print()
            analysis_results = self.run_sequential_analyses(fasta_files, output_path, threads, skip_modules, mlst_scheme)
            if not skip_summary:
                copy_result = self.copy_files_to_summary_module(output_path)
                if copy_result["success"]:
                    summary_success = self.run_summary_analysis(output_path)
                    analysis_results["Ultimate Reporter"] = summary_success
                    if summary_success:
                        ultimate_dir = output_path / "GENIUS_ACINETOBACTER_ULTIMATE_REPORTS"
                        if ultimate_dir.exists():
                            self.print_header("ANALYSIS COMPLETE", "All reports generated")
                            self.print_success(f"🎉 Ultimate reports available in: {ultimate_dir}")
                            html_files = list(ultimate_dir.glob("*.html"))
                            if html_files:
                                self.print_info("Main HTML reports:")
                                for html_file in sorted(html_files):
                                    self.print_info(f"  📄 {html_file.name}")
                else:
                    self.print_warning("Skipping ultimate reporter due to missing required files")
            analysis_time = datetime.now() - start_time
            analysis_time_str = str(analysis_time).split('.')[0]
            successful_count = sum(analysis_results.values())
            total_count = len(analysis_results)
            print()
            print(f"{Color.BOLD}{Color.BRIGHT_GREEN}{'='*80}{Color.RESET}")
            print(f"{Color.BOLD}{Color.BRIGHT_CYAN}{' '*25}ANALYSIS COMPLETE{Color.RESET}")
            print(f"{Color.BOLD}{Color.BRIGHT_GREEN}{'='*80}{Color.RESET}")
            print()
            print(f"{Color.BOLD}📊 Summary:{Color.RESET}")
            print(f"  ⏱️  Time elapsed: {analysis_time_str}")
            print(f"  🧫 Samples processed: {len(fasta_files)}")
            print(f"  ✅ Successful analyses: {successful_count}/{total_count}")
            print()
            print(f"{Color.BOLD}📁 Output directory:{Color.RESET} {output_path}")
            self.print_info("Generated directories:")
            for subdir in sorted(output_path.iterdir()):
                if subdir.is_dir():
                    file_count = len(list(subdir.glob("*")))
                    self.print_info(f"  📁 {subdir.name} ({file_count} files)")
            if successful_count == total_count:
                self.print_success(f"🎉 All analyses completed successfully!")
            else:
                self.print_warning(f"⚠️  {successful_count}/{total_count} analyses completed successfully.")
            print()
            self.print_header("INSPIRATION FOR THE JOURNEY", "Scientific Wisdom")
            self.display_random_quote()
        except KeyboardInterrupt:
            self.print_error("Analysis interrupted by user")
        except Exception as e:
            self.print_error(f"Critical error in analysis pipeline: {str(e)}")
            import traceback
            traceback.print_exc()

def main():
    # Handle --version first
    if '--version' in sys.argv or '-V' in sys.argv:
        print(f"AcinetoScope v{__version__}")
        sys.exit(0)

    parser = argparse.ArgumentParser(
        description="AcinetoScope: Complete A. baumannii Typing Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False
    )
    parser.add_argument('-i', '--input', help='Input FASTA file(s) - can use glob patterns like "*.fna" or "*.fasta"')
    parser.add_argument('-o', '--output', help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2, help='Number of threads (default: 2)')
    parser.add_argument('--mlst-scheme', choices=['pasteur', 'oxford', 'both'], default='both',
                        help='MLST scheme to use: pasteur, oxford, or both (default: both)')
    parser.add_argument('--skip-qc', action='store_true', help='Skip QC analysis')
    parser.add_argument('--skip-mlst', action='store_true', help='Skip MLST analysis')
    parser.add_argument('--skip-kaptive', action='store_true', help='Skip Kaptive analysis')
    parser.add_argument('--skip-amr', action='store_true', help='Skip AMR analysis')
    parser.add_argument('--skip-abricate', action='store_true', help='Skip ABRicate analysis')
    parser.add_argument('--skip-plasmid', action='store_true', help='Skip plasmid typing (APT) analysis')
    parser.add_argument('--skip-summary', action='store_true', help='Skip ultimate reporter generation')
    parser.add_argument('--force-update', action='store_true', help='Force update of AMR database (overwrites existing)')
    parser.add_argument('-h', '--help', action='store_true', help='Show this help message and exit')
    parser.add_argument('--version', '-V', action='store_true', help='Show version and exit')

    args = parser.parse_args()

    if args.help:
        orchestrator = AcinetoScopeOrchestrator()
        orchestrator.display_banner()
        print(f"\n{Color.BRIGHT_YELLOW}{Color.BOLD}USAGE:{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope{Color.RESET} {Color.CYAN}-i INPUT -o OUTPUT{Color.RESET} [OPTIONS]")
        print(f"\n{Color.BRIGHT_YELLOW}{Color.BOLD}REQUIRED ARGUMENTS:{Color.RESET}")
        print(f"  {Color.GREEN}-i, --input{Color.RESET} INPUT    Input FASTA file(s)")
        print(f"  {Color.GREEN}-o, --output{Color.RESET} OUTPUT  Output directory for results")
        print(f"\n{Color.BRIGHT_YELLOW}{Color.BOLD}OPTIONAL ARGUMENTS:{Color.RESET}")
        print(f"  {Color.GREEN}-t, --threads{Color.RESET} THREADS Number of threads (default: 2)")
        print(f"  {Color.GREEN}--mlst-scheme{Color.RESET} SCHEME  MLST scheme (pasteur/oxford/both)")
        print(f"  {Color.GREEN}--skip-qc{Color.RESET}           Skip QC analysis")
        print(f"  {Color.GREEN}--skip-mlst{Color.RESET}         Skip MLST analysis")
        print(f"  {Color.GREEN}--skip-kaptive{Color.RESET}      Skip Kaptive analysis")
        print(f"  {Color.GREEN}--skip-amr{Color.RESET}          Skip AMR analysis")
        print(f"  {Color.GREEN}--skip-abricate{Color.RESET}     Skip ABRicate analysis")
        print(f"  {Color.GREEN}--skip-plasmid{Color.RESET}      Skip plasmid typing (APT) analysis")
        print(f"  {Color.GREEN}--skip-summary{Color.RESET}      Skip ultimate reporter")
        print(f"  {Color.GREEN}--force-update{Color.RESET}      Force update of AMR database (overwrites existing)")
        print(f"  {Color.GREEN}-h, --help{Color.RESET}           Show this help message")
        print(f"  {Color.GREEN}--version, -V{Color.RESET}        Show version and exit")
        print(f"\n{Color.BRIGHT_YELLOW}Examples:{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope -i genome.fna -o results/{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope -i \"*.fna\" -o batch_results --threads 4{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope -i \"*.fasta\" -o analysis --threads 8 --skip-qc{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope -i \"genome*.fa\" -o results/ --threads 2 --skip-summary{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope -i \"*.fna\" -o results/ --mlst-scheme pasteur{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope --force-update{Color.RESET}")
        print(f"\n{Color.BRIGHT_YELLOW}🔧 DATABASE UPDATE NOTES:{Color.RESET}")
        print(f"  • {Color.CYAN}Conda users{Color.RESET}: Run {Color.GREEN}acinetoscope --force-update{Color.RESET} after installation to download the latest AMR database.")
        print(f"  • {Color.CYAN}Docker users{Color.RESET}: The database is updated automatically when pulling the latest image.")
        print(f"\n{Color.BRIGHT_YELLOW}Supported FASTA formats:{Color.RESET} {Color.CYAN}.fna, .fasta, .fa, .fn{Color.RESET}")
        print(f"\n{Color.BRIGHT_YELLOW}Critical Genes Tracked{Color.RESET}: Carbapenemases • ESBLs • Colistin Resistance • Tigecycline Resistance • Biofilm Formation • Efflux Pumps • Environmental Co-Selection")
        print(f"\n{Color.CYAN}⭐ Star us on GitHub if you find this tool useful!{Color.RESET}")
        sys.exit(0)

    if args.version:
        print(f"AcinetoScope v{__version__}")
        sys.exit(0)

    if args.force_update:
        orchestrator = AcinetoScopeOrchestrator()
        orchestrator.display_banner()
        success = orchestrator.run_force_update_amr()
        sys.exit(0 if success else 1)

    if not args.input or not args.output:
        print(f"{Color.BRIGHT_RED}Error: Both -i/--input and -o/--output are required for analysis.{Color.RESET}")
        print(f"Use {Color.GREEN}acinetoscope --help{Color.RESET} for usage information.")
        sys.exit(1)

    skip_modules = {
        'qc': args.skip_qc,
        'mlst': args.skip_mlst,
        'kaptive': args.skip_kaptive,
        'amr': args.skip_amr,
        'abricate': args.skip_abricate,
        'plasmid': args.skip_plasmid,
    }
    orchestrator = AcinetoScopeOrchestrator()
    orchestrator.run_complete_analysis(
        input_path=args.input,
        output_dir=args.output,
        threads=args.threads,
        skip_modules=skip_modules,
        mlst_scheme=args.mlst_scheme,
        skip_summary=args.skip_summary
    )

if __name__ == "__main__":
    main()