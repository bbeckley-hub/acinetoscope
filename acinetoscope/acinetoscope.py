#!/usr/bin/env python3
"""
AcinetoScope Main Orchestrator - Colored Sequential Execution with Scientific Quotes
Complete A. baumannii typing pipeline - MLST, K/O, AMR, Plasmid, Virulence, QC, Summary
Author: Brown Beckley <brownbeckley94@gmail.com>
Date: 2026-01-06
Send a quick mail for any issues or further explanations.
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
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

class Color:
    """ANSI color codes for colored output"""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    
    # Regular colors
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    
    # Bright colors
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
            Color.BRIGHT_CYAN,
            Color.BRIGHT_GREEN,
            Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA,
            Color.BRIGHT_BLUE,
            Color.BRIGHT_RED,
            Color.CYAN,
            Color.GREEN,
            Color.YELLOW,
            Color.MAGENTA
        ]
        
        # EXACT OUTPUT DIRECTORY NAMES
        self.output_dirs = {
            'qc': 'fasta_qc_results',
            'abricate': 'acineto_abricate_results',
            'amr': 'acineto_amrfinder_results',
            'kaptive': 'kaptive_results',
            'mlst_pasteur': 'PASTEUR_MLST',
            'mlst_oxford': 'OXFORD_MLST'
        }
        
        # EXACT HTML FILES REQUIRED FOR SUMMARY MODULE
        self.summary_html_files = {
            # MLST files
            'pasteur_mlst_summary.html': 'pasteur_mlst_summary.html',
            'oxford_mlst_summary.html': 'oxford_mlst_summary.html',
            # Kaptive file
            'Kaptive_summary.html': 'kaptive_results/Kaptive_summary.html',
            # AMR file
            'acineto_amrfinder_summary_report.html': 'acineto_amrfinder_results/acineto_amrfinder_summary_report.html',
            # ABRicate files
            'acineto_card_summary_report.html': 'acineto_abricate_results/acineto_card_summary_report.html',
            'acineto_ncbi_summary_report.html': 'acineto_abricate_results/acineto_ncbi_summary_report.html',
            'acineto_resfinder_summary_report.html': 'acineto_abricate_results/acineto_resfinder_summary_report.html',
            'acineto_vfdb_summary_report.html': 'acineto_abricate_results/acineto_vfdb_summary_report.html',
            'acineto_argannot_summary_report.html': 'acineto_abricate_results/acineto_argannot_summary_report.html',
            'acineto_megares_summary_report.html': 'acineto_abricate_results/acineto_megares_summary_report.html',
            'acineto_ecoli_vf_summary_report.html': 'acineto_abricate_results/acineto_ecoli_vf_summary_report.html',
            'acineto_bacmet2_summary_report.html': 'acineto_abricate_results/acineto_bacmet2_summary_report.html',
            'acineto_plasmidfinder_summary_report.html': 'acineto_abricate_results/acineto_bacmet2_summary_report.html',
            'acineto_ecoh_summary_report.html': 'acineto_abricate_results/acineto_bacmet2_summary_report.html',
            'acineto_victors_summary_report.html': 'acineto_abricate_results/acineto_victors_summary_report.html',
            'acineto_plasmidfinder_summary_report.html': 'acineto_abricate_results/acineto_plasmidfinder_summary_report.html'
        }
    
    def _get_scientific_quotes(self):
        """Curated scientific quotes about microbiology, genomics, and discovery"""
        return [
            # Short quotes for quick inspiration
            {
                "quote": "Science is organized knowledge.",
                "author": "Herbert Spencer",
                "theme": "knowledge"
            },
            {
                "quote": "The science of today is the technology of tomorrow.",
                "author": "Edward Teller",
                "theme": "technology"
            },
            {
                "quote": "Nature is the source of all true knowledge.",
                "author": "Leonardo da Vinci",
                "theme": "nature"
            },
            {
                "quote": "Biology is the most powerful technology ever created.",
                "author": "Freeman Dyson",
                "theme": "biology"
            },
            {
                "quote": "Genomics is a lens on biology.",
                "author": "Eric Lander",
                "theme": "genomics"
            },
            {
                "quote": "Every microbe has its own story.",
                "author": "Anonymous",
                "theme": "microbiology"
            },
            {
                "quote": "Data beats emotions.",
                "author": "Sean Rad",
                "theme": "data"
            },
            {
                "quote": "Code is poetry.",
                "author": "WordPress",
                "theme": "programming"
            },
            {
                "quote": "Sequence today, understand tomorrow.",
                "author": "Anonymous",
                "theme": "sequencing"
            },
            {
                "quote": "Microbes rule the world.",
                "author": "Paul Stamets",
                "theme": "microbiology"
            },
            {
                "quote": "In every drop, a universe.",
                "author": "Antonie van Leeuwenhoek",
                "theme": "microscopy"
            },
            {
                "quote": "Genes are the language of life.",
                "author": "Francis Collins",
                "theme": "genetics"
            },
            {
                "quote": "Resistance is not futile.",
                "author": "Antibiotic Researcher",
                "theme": "resistance"
            },
            {
                "quote": "Evolution in a petri dish.",
                "author": "Richard Lenski",
                "theme": "evolution"
            },
            {
                "quote": "Small things, big impact.",
                "author": "Microbiologist",
                "theme": "microbes"
            },
            # Original quotes
            {
                "quote": "In every walk with nature, one receives far more than he seeks.",
                "author": "John Muir",
                "theme": "discovery"
            },
            {
                "quote": "The microbe is nothing; the terrain is everything.",
                "author": "Louis Pasteur",
                "theme": "microbiology"
            },
            {
                "quote": "What we know is a drop, what we don't know is an ocean.",
                "author": "Isaac Newton",
                "theme": "knowledge"
            },
            {
                "quote": "The good physician treats the disease; the great physician treats the patient who has the disease.",
                "author": "William Osler",
                "theme": "medicine"
            },
            {
                "quote": "In science, there are no shortcuts to truth.",
                "author": "Karl Popper",
                "theme": "science"
            },
            {
                "quote": "The art of research is the art of making difficult problems soluble by devising means of getting at them.",
                "author": "Peter Medawar",
                "theme": "research"
            },
            {
                "quote": "Equipped with his five senses, man explores the universe around him and calls the adventure Science.",
                "author": "Edwin Hubble",
                "theme": "exploration"
            },
            {
                "quote": "The important thing is not to stop questioning. Curiosity has its own reason for existing.",
                "author": "Albert Einstein",
                "theme": "curiosity"
            },
            {
                "quote": "One must learn by doing the thing; though you think you know it, you have no certainty until you try.",
                "author": "Sophocles",
                "theme": "practice"
            },
            {
                "quote": "The secret of getting ahead is getting started.",
                "author": "Mark Twain",
                "theme": "motivation"
            },
            {
                "quote": "Nature is not a place to visit. It is home.",
                "author": "Gary Snyder",
                "theme": "nature"
            },
            {
                "quote": "Every scientific advance begins with the asking of a question.",
                "author": "Anonymous",
                "theme": "inquiry"
            },
            {
                "quote": "Nothing in life is to be feared, it is only to be understood. Now is the time to understand more, so that we may fear less.",
                "author": "Marie Curie",
                "theme": "understanding"
            },
            {
                "quote": "The greatest enemy of knowledge is not ignorance, it is the illusion of knowledge.",
                "author": "Stephen Hawking",
                "theme": "knowledge"
            },
            {
                "quote": "We are just an advanced breed of monkeys on a minor planet of a very average star. But we can understand the Universe. That makes us something very special.",
                "author": "Stephen Hawking",
                "theme": "perspective"
            },
            {
                "quote": "To raise new questions, new possibilities, to regard old problems from a new angle, requires creative imagination and marks real advance in science.",
                "author": "Albert Einstein",
                "theme": "innovation"
            },
            {
                "quote": "The most exciting phrase to hear in science, the one that heralds new discoveries, is not 'Eureka!' but 'That's funny...'",
                "author": "Isaac Asimov",
                "theme": "discovery"
            },
            {
                "quote": "In science the credit goes to the man who convinces the world, not to the man to whom the idea first occurs.",
                "author": "Francis Darwin",
                "theme": "recognition"
            },
            {
                "quote": "The aim of science is not to open the door to infinite wisdom, but to set a limit to infinite error.",
                "author": "Bertolt Brecht",
                "theme": "purpose"
            },
            {
                "quote": "Science knows no country, because knowledge belongs to humanity, and is the torch which illuminates the world.",
                "author": "Louis Pasteur",
                "theme": "global"
            },
            {
                "quote": "The first rule of intelligent tinkering is to save all the parts.",
                "author": "Paul Ehrlich",
                "theme": "conservation"
            },
            {
                "quote": "DNA is like a computer program but far, far more advanced than any software ever created.",
                "author": "Bill Gates",
                "theme": "genomics"
            },
            {
                "quote": "The human body is the best picture of the human soul.",
                "author": "Ludwig Wittgenstein",
                "theme": "medicine"
            },
            {
                "quote": "The beauty of a living thing is not the atoms that go into it, but the way those atoms are put together.",
                "author": "Carl Sagan",
                "theme": "biology"
            },
            {
                "quote": "A man who dares to waste one hour of time has not discovered the value of life.",
                "author": "Charles Darwin",
                "theme": "time"
            },
            {
                "quote": "If I have seen further it is by standing on the shoulders of Giants.",
                "author": "Isaac Newton",
                "theme": "collaboration"
            }
        ]
    
    def display_random_quote(self):
        """Display a random scientific quote with timestamp and colored formatting"""
        if not self.quotes:
            return
        
        quote_data = random.choice(self.quotes)
        quote = quote_data["quote"]
        author = quote_data["author"]
        theme = quote_data.get("theme", "science")
        
        # Choose a random color for this quote
        quote_color = random.choice(self.quote_colors)
        
        # Current date and time
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Quote icon based on theme
        theme_icons = {
            "microbiology": "🦠",
            "discovery": "🔬",
            "knowledge": "📚",
            "medicine": "⚕️",
            "science": "🧪",
            "research": "🔍",
            "exploration": "🚀",
            "curiosity": "🤔",
            "practice": "🛠️",
            "motivation": "💪",
            "nature": "🌿",
            "inquiry": "❓",
            "technology": "💻",
            "understanding": "🧠",
            "perspective": "👁️",
            "innovation": "💡",
            "recognition": "🏆",
            "purpose": "🎯",
            "biology": "🧬",
            "genomics": "🧬",
            "data": "📊",
            "programming": "💻",
            "sequencing": "🧬",
            "microscopy": "🔬",
            "genetics": "🧬",
            "resistance": "🛡️",
            "evolution": "🔄",
            "microbes": "🦠",
            "global": "🌍",
            "conservation": "🌱",
            "time": "⏳",
            "collaboration": "🤝"
        }
        
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
        """Setup color aliases for different message types"""
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
        """Print colored message"""
        style = Color.BOLD if bold else ''
        print(f"{style}{color}{message}{Color.RESET}")
    
    def print_header(self, title: str, subtitle: str = ""):
        """Print module header"""
        print()
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_CYAN}{' ' * 20}{title}{Color.RESET}")
        if subtitle:
            print(f"{Color.DIM}{Color.WHITE}{' ' * 22}{subtitle}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print()
    
    def print_info(self, message: str):
        """Print info message"""
        print(f"{self.color_info}[INFO]{Color.RESET} {message}")
    
    def print_success(self, message: str):
        """Print success message"""
        print(f"{self.color_success}✓{Color.RESET} {message}")
    
    def print_warning(self, message: str):
        """Print warning message"""
        print(f"{self.color_warning}⚠️{Color.RESET} {message}")
    
    def print_error(self, message: str):
        """Print error message"""
        print(f"{self.color_error}✗{Color.RESET} {message}")
    
    def print_command(self, command: str):
        """Print command being executed"""
        print(f"{Color.DIM}{Color.WHITE}  $ {command}{Color.RESET}")
    
    def display_banner(self):
        """Display AcinetoScope banner"""
        banner = f"""{Color.BOLD}{Color.BRIGHT_MAGENTA}
{'='*80}
{' '*20}🦠 ACINETOSCOPE - A. baumannii Genomic Analysis Pipeline
{'='*80}
{Color.RESET}{Color.BRIGHT_CYAN}
Complete A. baumannii genomic analysis pipeline
MLST | K/O Locus | AMR | Virulence | Plasmid | Quality Control | Critical Genes Flagging | Summary Reports
{Color.RESET}{Color.DIM}
Author: Brown Beckley | Email: brownbeckley94@gmail.com
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
{'='*80}{Color.RESET}
"""
        print(banner)
    
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files using glob patterns"""
        self.print_info(f"Searching for files with pattern: {input_path}")
        
        # Handle quoted wildcards properly
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files if Path(f).is_file() and 
                          f.lower().endswith(('.fna', '.fasta', '.fa', '.fn')) and
                          not Path(f).name.startswith('.')]
            self.print_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)
        
        # Handle direct file path
        input_path_obj = Path(input_path)
        if input_path_obj.is_file() and input_path_obj.suffix.lower() in ['.fna', '.fasta', '.fa', '.fn']:
            self.print_success(f"Found single FASTA file: {input_path_obj.name}")
            return [input_path_obj]
        
        # Handle directory
        if input_path_obj.is_dir():
            patterns = [
                f"{input_path}/*.fna", f"{input_path}/*.fasta",
                f"{input_path}/*.fa", f"{input_path}/*.fn"
            ]
            fasta_files = []
            for pattern in patterns:
                matched_files = glob.glob(pattern)
                for file_path in matched_files:
                    path = Path(file_path)
                    if path.is_file() and not path.name.startswith('.'):
                        fasta_files.append(path)
            fasta_files = sorted(list(set(fasta_files)))
            
            if fasta_files:
                self.print_success(f"Found {len(fasta_files)} FASTA files in directory")
            else:
                self.print_warning(f"No FASTA files found in directory: {input_path}")
            return fasta_files
        
        self.print_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        """Get the correct file pattern based on actual file extensions"""
        if not fasta_files:
            return '"*.fna"'  # Default fallback
        
        # Get all unique extensions from the input files
        extensions = set(f.suffix.lower() for f in fasta_files)
        
        # If all files have the same extension, use that
        if len(extensions) == 1:
            ext = list(extensions)[0]
            return f'"*{ext}"'
        
        # If mixed extensions, use a pattern that matches all FASTA files
        return '"*"'

    def cleanup_module_directory(self, module_path: Path, fasta_files: List[Path]):
        """Cleanup module directory after analysis - EXACT NAMES ONLY"""
        try:
            self.print_info(f"Cleaning up {module_path.name}...")
            
            # Remove copied input files
            for fasta_file in fasta_files:
                temp_file = module_path / fasta_file.name
                if temp_file.exists():
                    temp_file.unlink()
            
            # Remove exact output directories that might exist
            exact_dirs_to_remove = [
                "fasta_qc_results",          # QC module
                "acineto_abricate_results",  # ABRicate module
                "acineto_amrfinder_results", # AMR module
                "kaptive_results",           # Kaptive module
                "PASTEUR_MLST",              # MLST Pasteur
                "OXFORD_MLST",               # MLST Oxford
                "mlst_pasteur_results",      # MLST Pasteur wrapper
                "mlst_oxford_results",       # MLST Oxford wrapper
                "results"                    # MLST general results
            ]
            
            for dir_name in exact_dirs_to_remove:
                dir_path = module_path / dir_name
                if dir_path.exists():
                    shutil.rmtree(dir_path)
            
            # Remove specific HTML files (EXACT NAMES from summary_html_files dictionary)
            for html_file in self.summary_html_files.keys():
                html_path = module_path / html_file
                if html_path.exists():
                    html_path.unlink()
            
            # Also remove any other HTML files that might have been created
            for html_file in module_path.glob("*.html"):
                if html_file.is_file():
                    html_file.unlink()
            
            self.print_success(f"✅ {module_path.name} cleaned up successfully")
            
        except Exception as e:
            self.print_warning(f"⚠️  Partial cleanup issue in {module_path.name}: {str(e)}")

    def run_qc_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run QC analysis - EXACT OUTPUT: fasta_qc_results"""
        qc_module_path = self.base_dir / "modules" / "qc_module"
        
        try:
            self.print_header("FASTA QC ANALYSIS", "Comprehensive Quality Control")
            
            qc_script = qc_module_path / "acineto_fasta_qc.py"
            
            if not qc_script.exists():
                self.print_error(f"QC script not found at: {qc_script}")
                return False
            
            # Copy files to QC module directory
            for fasta_file in fasta_files:
                target_file = qc_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.print_info(f"Copied {len(fasta_files)} files to QC module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            file_pattern_clean = file_pattern.strip('"')
            
            # Build command - QC module uses pattern directly
            cmd = [sys.executable, str(qc_script), file_pattern_clean]
            
            self.print_info(f"Running QC analysis with pattern: {file_pattern}")
            self.print_command(f"python3 {qc_script.name} {file_pattern_clean}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=qc_module_path)
            
            if result.returncode == 0:
                self.print_success("QC analysis completed!")
                
                # Copy results to output directory - EXACT NAME: fasta_qc_results
                qc_source = qc_module_path / "fasta_qc_results"
                qc_target = output_dir / "fasta_qc_results"
                
                if qc_source.exists():
                    if qc_target.exists():
                        shutil.rmtree(qc_target)
                    shutil.copytree(qc_source, qc_target)
                    self.print_success(f"QC results copied to: {qc_target}")
                else:
                    self.print_warning("QC results directory not found: fasta_qc_results")
                
                # Display random scientific quote
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
            self.print_info("Cleaning up qc_module...")
            self.cleanup_module_directory(qc_module_path, fasta_files)

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int, scheme: str) -> bool:
        """Run MLST analysis for specific scheme - SIMPLIFIED OUTPUT"""
        mlst_module_path = self.base_dir / "modules" / "mlst_module"
        
        try:
            scheme_name = "PASTEUR" if scheme == "pasteur" else "OXFORD"
            self.print_header(f"MLST ANALYSIS - {scheme_name}", "Multi-Locus Sequence Typing")
            
            mlst_script = mlst_module_path / "mlst_module.py"
            
            if not mlst_script.exists():
                self.print_error(f"MLST script not found at: {mlst_script}")
                return False
            
            # Copy files to MLST module directory
            for fasta_file in fasta_files:
                target_file = mlst_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.print_info(f"Copied {len(fasta_files)} files to MLST module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            file_pattern_clean = file_pattern.strip('"')
            
            # Create output directory name based on scheme
            output_subdir = f"mlst_{scheme}_results"
            
            # Build command - MLST module uses specific arguments
            cmd = [
                sys.executable, str(mlst_script),
                "-i", file_pattern_clean,
                "-o", output_subdir,
                "-db", "db",
                "-sc", "bin",
                "--batch",
                "-s", scheme
            ]
            
            self.print_info(f"Running MLST analysis with scheme: {scheme_name}")
            self.print_command(f"python {mlst_script.name} -i {file_pattern_clean} -o {output_subdir} -db db -sc bin --batch -s {scheme}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=mlst_module_path)
            
            if result.returncode == 0:
                self.print_success(f"MLST analysis completed for {scheme_name} scheme!")
                
                # MLST creates: mlst_pasteur_results/PASTEUR_MLST/ OR mlst_oxford_results/OXFORD_MLST/
                scheme_dir = "PASTEUR_MLST" if scheme == "pasteur" else "OXFORD_MLST"
                mlst_source = mlst_module_path / output_subdir / scheme_dir
                mlst_target = output_dir / scheme_dir
                
                if mlst_source.exists():
                    if mlst_target.exists():
                        shutil.rmtree(mlst_target)
                    shutil.copytree(mlst_source, mlst_target)
                    self.print_success(f"MLST results copied to: {mlst_target}")
                    
                    # Count files in the directory
                    file_count = len(list(mlst_target.glob("*")))
                    self.print_info(f"  Contains {file_count} files")
                else:
                    self.print_warning(f"MLST results directory not found: {mlst_source}")
                
                # Copy the HTML summary file
                html_filename = f"{scheme}_mlst_summary.html"
                html_source = mlst_source / html_filename
                if html_source.exists():
                    html_target = output_dir / html_filename
                    shutil.copy2(html_source, html_target)
                    self.print_success(f"MLST HTML report copied to: {html_target}")
                else:
                    self.print_warning(f"MLST HTML report not found: {html_filename}")
                
                # Display random scientific quote
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
            self.print_info("Cleaning up mlst_module...")
            self.cleanup_module_directory(mlst_module_path, fasta_files)

    def run_kaptive_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run Kaptive K/O locus analysis - EXACT OUTPUT: kaptive_results"""
        kaptive_module_path = self.base_dir / "modules" / "k_o_module"
        
        try:
            self.print_header("KAPTIVE ANALYSIS", "K and O Locus Typing")
            
            kaptive_script = kaptive_module_path / "acineto_kaptive.py"
            
            if not kaptive_script.exists():
                self.print_error(f"Kaptive script not found at: {kaptive_script}")
                return False
            
            # Copy files to Kaptive module directory
            for fasta_file in fasta_files:
                target_file = kaptive_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.print_info(f"Copied {len(fasta_files)} files to Kaptive module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            file_pattern_clean = file_pattern.strip('"')
            
            # Build command - Kaptive uses pattern directly
            cmd = [sys.executable, str(kaptive_script), file_pattern_clean]
            
            self.print_info(f"Running Kaptive analysis with pattern: {file_pattern}")
            self.print_command(f"python3 {kaptive_script.name} {file_pattern_clean}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=kaptive_module_path)
            
            if result.returncode == 0:
                self.print_success("Kaptive analysis completed!")
                
                # Copy results to output directory - EXACT NAME: kaptive_results
                kaptive_source = kaptive_module_path / "kaptive_results"
                kaptive_target = output_dir / "kaptive_results"
                
                if kaptive_source.exists():
                    if kaptive_target.exists():
                        shutil.rmtree(kaptive_target)
                    shutil.copytree(kaptive_source, kaptive_target)
                    self.print_success(f"Kaptive results copied to: {kaptive_target}")
                else:
                    self.print_warning("Kaptive results directory not found: kaptive_results")
                
                # Copy specific HTML summary file
                html_source = kaptive_module_path / "kaptive_results" / "Kaptive_summary.html"
                if html_source.exists():
                    html_target = output_dir / "Kaptive_summary.html"
                    shutil.copy2(html_source, html_target)
                    self.print_success(f"Kaptive HTML report copied to: {html_target}")
                else:
                    self.print_warning("Kaptive HTML report not found: Kaptive_summary.html")
                
                # Display random scientific quote
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
            self.print_info("Cleaning up k_o_module...")
            self.cleanup_module_directory(kaptive_module_path, fasta_files)

    def run_amr_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run AMR analysis - EXACT OUTPUT: acineto_amrfinder_results"""
        amr_module_path = self.base_dir / "modules" / "amr_module"
        
        try:
            self.print_header("AMR ANALYSIS", "Antimicrobial Resistance Gene Detection")
            
            amr_script = amr_module_path / "acineto_amrfinder.py"
            
            if not amr_script.exists():
                self.print_error(f"AMR script not found at: {amr_script}")
                return False
            
            # Copy files to AMR module directory
            for fasta_file in fasta_files:
                target_file = amr_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.print_info(f"Copied {len(fasta_files)} files to AMR module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            file_pattern_clean = file_pattern.strip('"')
            
            # Build command - AMR uses pattern directly
            cmd = [sys.executable, str(amr_script), file_pattern_clean]
            
            self.print_info(f"Running AMR analysis with pattern: {file_pattern}")
            self.print_command(f"python3 {amr_script.name} {file_pattern_clean}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
            
            if result.returncode == 0:
                self.print_success("AMR analysis completed!")
                
                # Copy results to output directory - EXACT NAME: acineto_amrfinder_results
                amr_source = amr_module_path / "acineto_amrfinder_results"
                amr_target = output_dir / "acineto_amrfinder_results"
                
                if amr_source.exists():
                    if amr_target.exists():
                        shutil.rmtree(amr_target)
                    shutil.copytree(amr_source, amr_target)
                    self.print_success(f"AMR results copied to: {amr_target}")
                else:
                    self.print_warning("AMR results directory not found: acineto_amrfinder_results")
                
                # Copy specific HTML summary file
                html_source = amr_module_path / "acineto_amrfinder_results" / "acineto_amrfinder_summary_report.html"
                if html_source.exists():
                    html_target = output_dir / "acineto_amrfinder_summary_report.html"
                    shutil.copy2(html_source, html_target)
                    self.print_success(f"AMR HTML report copied to: {html_target}")
                else:
                    self.print_warning("AMR HTML report not found: acineto_amrfinder_summary_report.html")
                
                # Display random scientific quote
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
            self.print_info("Cleaning up amr_module...")
            self.cleanup_module_directory(amr_module_path, fasta_files)

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        """Run ABRicate analysis - EXACT OUTPUT: acineto_abricate_results"""
        abricate_module_path = self.base_dir / "modules" / "abricate_module"
        
        try:
            self.print_header("ABRICATE ANALYSIS", "Comprehensive Resistance & Virulence Gene Screening & Plasmid Profiling")
            
            abricate_script = abricate_module_path / "acineto_abricate.py"
            
            if not abricate_script.exists():
                self.print_error(f"ABRicate script not found at: {abricate_script}")
                return False
            
            # Copy files to ABRicate module directory
            for fasta_file in fasta_files:
                target_file = abricate_module_path / fasta_file.name
                shutil.copy2(fasta_file, target_file)
            
            self.print_info(f"Copied {len(fasta_files)} files to ABRicate module")
            
            # Get correct file pattern based on actual files
            file_pattern = self.get_file_pattern(fasta_files)
            file_pattern_clean = file_pattern.strip('"')
            
            # Build command - ABRicate uses pattern directly
            cmd = [sys.executable, str(abricate_script), file_pattern_clean]
            
            self.print_info(f"Running ABRicate analysis with pattern: {file_pattern}")
            self.print_command(f"python3 {abricate_script.name} {file_pattern_clean}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=abricate_module_path)
            
            if result.returncode == 0:
                self.print_success("ABRicate analysis completed!")
                
                # Copy results to output directory - EXACT NAME: acineto_abricate_results
                abricate_source = abricate_module_path / "acineto_abricate_results"
                abricate_target = output_dir / "acineto_abricate_results"
                
                if abricate_source.exists():
                    if abricate_target.exists():
                        shutil.rmtree(abricate_target)
                    shutil.copytree(abricate_source, abricate_target)
                    self.print_success(f"ABRicate results copied to: {abricate_target}")
                else:
                    self.print_warning("ABRicate results directory not found: acineto_abricate_results")
                
                # Copy all HTML summary files (EXACT NAMES from summary_html_files dictionary)
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
                
                # Display random scientific quote
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
            self.print_info("Cleaning up abricate_module...")
            self.cleanup_module_directory(abricate_module_path, fasta_files)

    def copy_files_to_summary_module(self, output_dir: Path) -> Dict[str, bool]:
        """Copy required files to summary_module for ultimate reporter - EXACT NAMES ONLY"""
        try:
            self.print_header("PREPARING SUMMARY MODULE", "Copying required HTML files")
            
            summary_module_path = self.base_dir / "modules" / "summary_module"
            required_files = {}
            
            copied_count = 0
            missing_count = 0
            
            # Copy each required HTML file using exact names from dictionary
            for target_filename, relative_path in self.summary_html_files.items():
                # Build source path
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
            
            # Check critical files
            critical_files = [
                "pasteur_mlst_summary.html",
                "oxford_mlst_summary.html", 
                "Kaptive_summary.html",
                "acineto_amrfinder_summary_report.html",
                "acineto_card_summary_report.html"
            ]
            
            missing_critical = [f for f in critical_files if not required_files.get(f, False)]
            
            if missing_critical:
                self.print_warning(f"Missing critical files: {', '.join(missing_critical)}")
                return {"success": False, "files": required_files}
            
            return {"success": True, "files": required_files}
            
        except Exception as e:
            self.print_error(f"Error copying files to summary module: {str(e)}")
            return {"success": False, "files": {}}

    def run_summary_analysis(self, output_dir: Path) -> bool:
        """Run ultimate reporter - EXACT OUTPUT: GENIUS_ACINETOBACTER_ULTIMATE_REPORTS"""
        try:
            self.print_header("ULTIMATE REPORTER", "Gene-centric Integrated Analysis")
            
            summary_module_path = self.base_dir / "modules" / "summary_module"
            summary_script = summary_module_path / "genius_acinetobacter_reporter.py"
            
            if not summary_script.exists():
                self.print_error(f"Summary script not found at: {summary_script}")
                return False
            
            # Build command
            cmd = [sys.executable, str(summary_script), "-i", "."]
            
            self.print_info("Running ultimate reporter...")
            self.print_command(f"python3 {summary_script.name} -i .")
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=summary_module_path)
            
            if result.returncode == 0:
                self.print_success("Ultimate reporter completed successfully!")
                
                # Copy ultimate reports to output directory - EXACT NAME: GENIUS_ACINETOBACTER_ULTIMATE_REPORTS
                ultimate_source = summary_module_path / "GENIUS_ACINETOBACTER_ULTIMATE_REPORTS"
                ultimate_target = output_dir / "GENIUS_ACINETOBACTER_ULTIMATE_REPORTS"
                
                if ultimate_source.exists():
                    if ultimate_target.exists():
                        shutil.rmtree(ultimate_target)
                    shutil.copytree(ultimate_source, ultimate_target)
                    
                    # Count files
                    report_files = list(ultimate_target.glob("*"))
                    html_count = len([f for f in report_files if f.suffix == '.html'])
                    json_count = len([f for f in report_files if f.suffix == '.json'])
                    csv_count = len([f for f in report_files if f.suffix == '.csv'])
                    
                    self.print_success(f"Ultimate reports copied to: {ultimate_target}")
                    self.print_info(f"  📊 Reports: {html_count} HTML, {json_count} JSON, {csv_count} CSV files")
                else:
                    self.print_warning("Ultimate reports directory not found: GENIUS_ACINETOBACTER_ULTIMATE_REPORTS")
                
                # Display random scientific quote
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
        """Run analyses sequentially with colored output"""
        analysis_functions = []
        
        # Define analysis functions and their skip flags
        if not skip_modules.get('qc', False):
            analysis_functions.append(("QC Analysis", self.run_qc_analysis, True))
        
        # MLST analysis
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
        
        if not analysis_functions:
            self.print_warning("All analyses were skipped! Nothing to run.")
            return {}
        
        self.print_info(f"Running {len(analysis_functions)} analyses sequentially")
        
        results = {}
        
        # Run analyses SEQUENTIALLY
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

    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1, 
                             skip_modules: Dict[str, bool] = None, mlst_scheme: str = "both",
                             skip_summary: bool = False):
        """Run complete AcinetoScope analysis pipeline"""
        if skip_modules is None:
            skip_modules = {}
        
        start_time = datetime.now()
        
        try:
            # Display banner
            self.display_banner()
            
            # Create output directory
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            
            # Find input files
            fasta_files = self.find_fasta_files(input_path)
            
            if not fasta_files:
                self.print_error("No FASTA files found! Analysis stopped.")
                return
            
            # Show file formats detected
            extensions = set(f.suffix.lower() for f in fasta_files)
            self.print_success(f"Starting analysis of {len(fasta_files)} A. baumannii samples")
            self.print_info(f"File formats detected: {', '.join(extensions)}")
            
            # Create output subdirectories with EXACT NAMES
            subdirs = [
                "fasta_qc_results",                    # QC module
                "PASTEUR_MLST",                        # MLST Pasteur
                "OXFORD_MLST",                         # MLST Oxford
                "kaptive_results",                     # Kaptive module
                "acineto_amrfinder_results",           # AMR module
                "acineto_abricate_results",            # ABRicate module
                "GENIUS_ACINETOBACTER_ULTIMATE_REPORTS" # Summary module output
            ]
            for subdir in subdirs:
                (output_path / subdir).mkdir(exist_ok=True)
            
            # Display analysis plan
            self.print_header("ANALYSIS PLAN", "Modules to be executed")
            
            analyses_to_run = [
                ("QC Analysis", not skip_modules.get('qc', False)),
                ("MLST Pasteur", not skip_modules.get('mlst', False) and mlst_scheme in ["pasteur", "both"]),
                ("MLST Oxford", not skip_modules.get('mlst', False) and mlst_scheme in ["oxford", "both"]),
                ("Kaptive Analysis", not skip_modules.get('kaptive', False)),
                ("AMR Analysis", not skip_modules.get('amr', False)),
                ("ABRicate Analysis", not skip_modules.get('abricate', False)),
                ("Ultimate Reporter", not skip_summary),
            ]
            
            for analysis, enabled in analyses_to_run:
                if enabled:
                    print(f"   {Color.BRIGHT_GREEN}✅ ENABLED{Color.RESET} - {analysis}")
                else:
                    print(f"   {Color.YELLOW}⏸️  SKIPPED{Color.RESET} - {analysis}")
            
            print()
            
            # Run main analyses SEQUENTIALLY
            analysis_results = self.run_sequential_analyses(
                fasta_files, output_path, threads, skip_modules, mlst_scheme
            )
            
            # Run summary analysis if not skipped
            if not skip_summary:
                # Copy required files to summary module
                copy_result = self.copy_files_to_summary_module(output_path)
                
                if copy_result["success"]:
                    summary_success = self.run_summary_analysis(output_path)
                    analysis_results["Ultimate Reporter"] = summary_success
                    
                    if summary_success:
                        # Display final report location
                        ultimate_dir = output_path / "GENIUS_ACINETOBACTER_ULTIMATE_REPORTS"
                        if ultimate_dir.exists():
                            self.print_header("ANALYSIS COMPLETE", "All reports generated")
                            self.print_success(f"🎉 Ultimate reports available in: {ultimate_dir}")
                            
                            # List generated files
                            html_files = list(ultimate_dir.glob("*.html"))
                            if html_files:
                                self.print_info("Main HTML reports:")
                                for html_file in sorted(html_files):
                                    self.print_info(f"  📄 {html_file.name}")
                else:
                    self.print_warning("Skipping ultimate reporter due to missing required files")
            
            # Calculate analysis time
            analysis_time = datetime.now() - start_time
            analysis_time_str = str(analysis_time).split('.')[0]
            
            # Display completion summary
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
            
            # List output directories with exact names
            self.print_info("Generated directories:")
            for subdir in sorted(output_path.iterdir()):
                if subdir.is_dir():
                    file_count = len(list(subdir.glob("*")))
                    self.print_info(f"  📁 {subdir.name} ({file_count} files)")
            
            if successful_count == total_count:
                self.print_success(f"🎉 All analyses completed successfully!")
            else:
                self.print_warning(f"⚠️  {successful_count}/{total_count} analyses completed successfully.")
            
            # Display a final inspirational quote
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
    """Main entry point for AcinetoScope"""
    
    # Check for help flag FIRST before any argparse processing
    if '-h' in sys.argv or '--help' in sys.argv:
        # Create orchestrator and display banner
        orchestrator = AcinetoScopeOrchestrator()
        orchestrator.display_banner()
        
        # Display custom colored help
        print(f"\n{Color.BRIGHT_YELLOW}{Color.BOLD}USAGE:{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope{Color.RESET} {Color.CYAN}-i INPUT -o OUTPUT{Color.RESET} [OPTIONS]")
        
        print(f"\n{Color.BRIGHT_YELLOW}{Color.BOLD}REQUIRED ARGUMENTS:{Color.RESET}")
        print(f"  {Color.GREEN}-i, --input{Color.RESET} INPUT    Input FASTA file(s)")
        print(f"  {Color.GREEN}-o, --output{Color.RESET} OUTPUT  Output directory for results")
        
        print(f"\n{Color.BRIGHT_YELLOW}{Color.BOLD}OPTIONAL ARGUMENTS:{Color.RESET}")
        print(f"  {Color.GREEN}-h, --help{Color.RESET}           Show this help message")
        print(f"  {Color.GREEN}-t, --threads{Color.RESET} THREADS Number of threads (default: 2)")
        
        print(f"\n{Color.BRIGHT_YELLOW}{Color.BOLD}ANALYSIS OPTIONS:{Color.RESET}")
        print(f"  {Color.GREEN}--mlst-scheme{Color.RESET} SCHEME  MLST scheme (pasteur/oxford/both)")
        print(f"  {Color.GREEN}--skip-qc{Color.RESET}           Skip QC analysis")
        print(f"  {Color.GREEN}--skip-mlst{Color.RESET}         Skip MLST analysis")
        print(f"  {Color.GREEN}--skip-kaptive{Color.RESET}      Skip Kaptive analysis")
        print(f"  {Color.GREEN}--skip-amr{Color.RESET}          Skip AMR analysis")
        print(f"  {Color.GREEN}--skip-abricate{Color.RESET}     Skip ABRicate analysis")
        print(f"  {Color.GREEN}--skip-summary{Color.RESET}      Skip ultimate reporter")
        
        # Display examples and other info
        print(f"\n{Color.BRIGHT_YELLOW}Examples:{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope -i genome.fna -o results/{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope -i \"*.fna\" -o batch_results --threads 4{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope -i \"*.fasta\" -o analysis --threads 8 --skip-qc{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope -i \"genome*.fa\" -o results/ --threads 2 --skip-summary{Color.RESET}")
        print(f"  {Color.GREEN}acinetoscope -i \"*.fna\" -o results/ --mlst-scheme pasteur{Color.RESET}")
        
        print(f"\n{Color.BRIGHT_YELLOW}Supported FASTA formats:{Color.RESET} {Color.CYAN}.fna, .fasta, .fa, .fn{Color.RESET}")
        
        print(f"\n{Color.BRIGHT_YELLOW}Analysis Modules:{Color.RESET}")
        print(f"  {Color.GREEN}•{Color.RESET} QC Analysis (Quality control)")
        print(f"  {Color.GREEN}•{Color.RESET} MLST Analysis (Pasteur and/or Oxford schemes)")
        print(f"  {Color.GREEN}•{Color.RESET} Kaptive Analysis (K and O locus typing)")
        print(f"  {Color.GREEN}•{Color.RESET} AMR Analysis (Antimicrobial resistance)")
        print(f"  {Color.GREEN}•{Color.RESET} ABRicate Analysis (Comprehensive resistance/virulence screening/Plasmid screening)")
        print(f"  {Color.GREEN}•{Color.RESET} Ultimate Reporter (Gene-centric integrated analysis)")
        print(f"  {Color.GREEN}•{Color.RESET} Critical Genes Tracked: Carbapenemases • ESBLs • Colistin Resistance • Tigecycline Resistance • Biofilm Formation • Efflux Pumps • Environmental Co-Selection")
        
        print(f"\n{Color.BRIGHT_YELLOW}Run this command prior analysis:{Color.RESET} {Color.CYAN} abricate --setupdb{Color.RESET}")

        print(f"\n{Color.BRIGHT_YELLOW}Output:{Color.RESET} Comprehensive results for all analyses in organized directories")
        print(f"{Color.CYAN}⭐ Star us on GitHub if you find this tool useful!{Color.RESET}")
        sys.exit(0)
    
    # Now parse arguments normally
    parser = argparse.ArgumentParser(
        description="AcinetoScope: Complete A. baumannii Typing Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False  
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input FASTA file(s) - can use glob patterns like "*.fna" or "*.fasta"')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2,
                       help='Number of threads (default: 2)')
    
    # MLST scheme option
    parser.add_argument('--mlst-scheme', choices=['pasteur', 'oxford', 'both'], default='both',
                       help='MLST scheme to use: pasteur, oxford, or both (default: both)')
    
    # Skip options
    parser.add_argument('--skip-qc', action='store_true',
                       help='Skip QC analysis')
    parser.add_argument('--skip-mlst', action='store_true',
                       help='Skip MLST analysis')
    parser.add_argument('--skip-kaptive', action='store_true',
                       help='Skip Kaptive analysis')
    parser.add_argument('--skip-amr', action='store_true',
                       help='Skip AMR analysis')
    parser.add_argument('--skip-abricate', action='store_true',
                       help='Skip ABRicate analysis')
    parser.add_argument('--skip-summary', action='store_true',
                       help='Skip ultimate reporter generation')
    
    args = parser.parse_args()
    
    # Create skip modules dictionary
    skip_modules = {
        'qc': args.skip_qc,
        'mlst': args.skip_mlst,
        'kaptive': args.skip_kaptive,
        'amr': args.skip_amr,
        'abricate': args.skip_abricate
    }
    
    # Create and run AcinetoScope
    acinetoscope = AcinetoScopeOrchestrator()
    
    try:
        acinetoscope.run_complete_analysis(
            input_path=args.input,
            output_dir=args.output,
            threads=args.threads,
            skip_modules=skip_modules,
            mlst_scheme=args.mlst_scheme,
            skip_summary=args.skip_summary
        )
    except KeyboardInterrupt:
        print(f"\n{Color.BRIGHT_RED}❌ Analysis interrupted by user{Color.RESET}")
    except Exception as e:
        print(f"\n{Color.BRIGHT_RED}💥 Critical error: {e}{Color.RESET}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()