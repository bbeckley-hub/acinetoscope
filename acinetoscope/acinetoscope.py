#!/usr/bin/env python3
"""
AcinetoScope Main Orchestrator - v1.3.2
All module writes happen in /tmp, final results are copied to user output.
HPC/ Docker-friendly with full flag support.

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Date: 2026-08-22
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import tempfile
import logging
import traceback
import signal
import random
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

__version__ = "1.3.2"


# ----------------------------------------------------------------------
# Color definitions
# ----------------------------------------------------------------------
class Color:
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


# ----------------------------------------------------------------------
# Custom Help Formatter with colors
# ----------------------------------------------------------------------
class ColoredHelpFormatter(argparse.HelpFormatter):
    HEADER = Color.BRIGHT_MAGENTA
    BLUE = Color.BRIGHT_BLUE
    GREEN = Color.BRIGHT_GREEN
    YELLOW = Color.BRIGHT_YELLOW
    RED = Color.BRIGHT_RED
    CYAN = Color.BRIGHT_CYAN
    BOLD = Color.BOLD
    RESET = Color.RESET

    def _format_usage(self, usage, actions, groups, prefix):
        usage_str = super()._format_usage(usage, actions, groups, prefix)
        if usage_str:
            usage_str = self.GREEN + self.BOLD + "usage: " + self.RESET + usage_str
        return usage_str

    def _format_action(self, action):
        action_str = super()._format_action(action)
        if not action_str:
            return action_str
        lines = action_str.split('\n')
        colored_lines = []
        for line in lines:
            if line.strip():
                if line.lstrip().startswith('-'):
                    parts = line.split('  ', 1)
                    if len(parts) == 2:
                        options = parts[0].strip()
                        help_text = parts[1]
                        colored_line = f"  {self.CYAN}{options}{self.RESET}  {help_text}"
                    else:
                        colored_line = f"  {self.CYAN}{line.strip()}{self.RESET}"
                else:
                    colored_line = f"  {self.YELLOW}{line}{self.RESET}"
                colored_lines.append(colored_line)
            else:
                colored_lines.append(line)
        return '\n'.join(colored_lines)

    def _format_text(self, text):
        if not text:
            return text
        return f"{self.BLUE}{text}{self.RESET}"

    def start_section(self, heading):
        heading = f"{self.BOLD}{self.GREEN}{heading}{self.RESET}"
        super().start_section(heading)


# ----------------------------------------------------------------------
# Main Orchestrator
# ----------------------------------------------------------------------
class AcinetoScopeOrchestrator:
    """
    HPC/Docker‑friendly orchestrator for AcinetoScope.
    Each module runs in a temporary directory; final reports are collected.
    """

    def __init__(self):
        self.base_dir = Path(__file__).parent
        self.user_output_dir = None
        self.logger = None
        self.keep_temp = False
        self.temp_dirs = []

        self.quotes = self._get_scientific_quotes()
        self.quote_colors = [
            Color.BRIGHT_CYAN, Color.BRIGHT_GREEN, Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA, Color.BRIGHT_BLUE, Color.BRIGHT_RED,
            Color.CYAN, Color.GREEN, Color.YELLOW, Color.MAGENTA
        ]

        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

    # ------------------------------------------------------------------
    # Signal handling and logging
    # ------------------------------------------------------------------
    def _signal_handler(self, signum, frame):
        self.print_color("Interrupted by user. Cleaning up temporary directories...", Color.BRIGHT_RED)
        for temp_dir in self.temp_dirs:
            if Path(temp_dir).exists():
                shutil.rmtree(temp_dir, ignore_errors=True)
                if self.logger:
                    self.logger.info(f"Removed temporary directory: {temp_dir}")
        sys.exit(1)

    def _register_temp_dir(self, path: Path):
        self.temp_dirs.append(str(path))

    def setup_logging(self, output_dir: Path):
        log_file = output_dir / "acinetoscope_run.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler(log_file, mode='w')]
        )
        self.logger = logging.getLogger("AcinetoScope")
        self.logger.info(f"Logging to {log_file}")
        self.user_output_dir = output_dir

    # ------------------------------------------------------------------
    # Colored printing (console with timestamps and quotes)
    # ------------------------------------------------------------------
    def print_color(self, message: str, color: str = Color.RESET, bold: bool = False):
        style = Color.BOLD if bold else ''
        print(f"{style}{color}{message}{Color.RESET}")

    def print_header(self, title: str, subtitle: str = ""):
        print()
        self.print_color("=" * 80, Color.BRIGHT_BLUE, bold=True)
        self.print_color(f"{' ' * 20}{title}", Color.BRIGHT_CYAN, bold=True)
        if subtitle:
            self.print_color(f"{' ' * 22}{subtitle}", Color.DIM + Color.WHITE)
        self.print_color("=" * 80, Color.BRIGHT_BLUE, bold=True)
        print()

    def print_info(self, message: str):
        self.print_color(f"[INFO] {message}", Color.CYAN)

    def print_success(self, message: str):
        self.print_color(f"✓ {message}", Color.BRIGHT_GREEN)

    def print_warning(self, message: str):
        self.print_color(f"⚠️  {message}", Color.BRIGHT_YELLOW)

    def print_error(self, message: str):
        self.print_color(f"✗ {message}", Color.BRIGHT_RED)

    def print_command(self, command: str):
        # Shorten command for display
        parts = command.split()
        if len(parts) > 1 and (parts[0].endswith('python') or parts[0].endswith('python3') or 'python' in parts[0]):
            script = Path(parts[1]).name if len(parts) > 1 else "script"
            args = " ".join(parts[2:]) if len(parts) > 2 else ""
            display = f"$ {script} {args}".strip()
        else:
            display = command
        self.print_color(f"  {display}", Color.DIM + Color.WHITE)

    def display_random_quote(self):
        if not self.quotes:
            return
        quote_data = random.choice(self.quotes)
        quote = quote_data["quote"]
        author = quote_data["author"]
        theme = quote_data.get("theme", "science")
        quote_color = random.choice(self.quote_colors)
        theme_icons = {
            "microbiology": "🦠", "discovery": "🔬", "knowledge": "📚", "medicine": "⚕️",
            "science": "🧪", "research": "🔍", "exploration": "🚀", "curiosity": "🤔",
            "practice": "🛠️", "motivation": "💪", "nature": "🌿", "inquiry": "❓",
            "technology": "💻", "understanding": "🧠", "perspective": "👁️", "innovation": "💡",
            "recognition": "🏆", "purpose": "🎯", "biology": "🧬", "genomics": "🧬",
            "data": "📊", "global": "🌍", "conservation": "🌱", "collaboration": "🤝"
        }
        icon = theme_icons.get(theme, "💭")
        print()
        self.print_color("─" * 80, Color.DIM + Color.WHITE)
        self.print_color(f"{icon} SCIENTIFIC INSIGHT:", Color.DIM + Color.WHITE)
        print()
        self.print_color(f"   \"{quote}\"", quote_color)
        self.print_color(f"   — {author}", Color.BOLD + Color.WHITE)
        self.print_color(f"   Theme: {theme.capitalize()}", Color.DIM + Color.WHITE)
        self.print_color("─" * 80, Color.DIM + Color.WHITE)
        print()

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

    # ------------------------------------------------------------------
    # Scientific quotes
    # ------------------------------------------------------------------
    def _get_scientific_quotes(self) -> List[Dict[str, str]]:
        return [
            {"quote": "Science is organized knowledge.", "author": "Herbert Spencer", "theme": "knowledge"},
            {"quote": "The science of today is the technology of tomorrow.", "author": "Edward Teller", "theme": "technology"},
            {"quote": "Nature is the source of all true knowledge.", "author": "Leonardo da Vinci", "theme": "nature"},
            {"quote": "Biology is the most powerful technology ever created.", "author": "Freeman Dyson", "theme": "biology"},
            {"quote": "Genomics is a lens on biology.", "author": "Eric Lander", "theme": "genomics"},
            {"quote": "Every microbe has its own story.", "author": "Anonymous", "theme": "microbiology"},
            {"quote": "Data beats emotions.", "author": "Sean Rad", "theme": "data"},
            {"quote": "In every walk with nature, one receives far more than one seeks.", "author": "John Muir", "theme": "discovery"},
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

    # ------------------------------------------------------------------
    # Module execution in temporary directory
    # ------------------------------------------------------------------
    def get_module_path(self, module_name: str) -> Path:
        if hasattr(sys, 'prefix'):
            share_path = Path(sys.prefix) / "share" / "acinetoscope" / "modules" / module_name
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        if not fasta_files:
            return '"*.fna"'
        extensions = set(f.suffix.lower() for f in fasta_files)
        if len(extensions) == 1:
            ext = list(extensions)[0]
            return f'"*{ext}"'
        return '"*"'

    # ================================================================
    # MODIFIED: run_module_in_temp now copies ref_db/ to temp root
    # ================================================================
    def run_module_in_temp(self, module_name: str, fasta_files: List[Path],
                           cmd_str: str, result_subdir: str = None,
                           extra_files: List[str] = None) -> bool:
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.print_error(f"Module directory not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"acinetoscope_{module_name}_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for {module_name}: {temp_dir}")

        try:
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True)

            # ------------------------------------------------------------
            # NEW: copy ref_db/ to temp_dir root so QC script finds it
            # ------------------------------------------------------------
            ref_db_src = module_orig / "ref_db"
            if ref_db_src.exists():
                shutil.copytree(ref_db_src, temp_dir / "ref_db")
                self.logger.info(f"Copied ref_db to temporary directory root: {temp_dir / 'ref_db'}")
            else:
                self.logger.warning(f"ref_db not found in {module_orig}, species confirmation will be disabled.")

            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.print_info(f"Running {module_name}...")
            self.print_command(cmd_str)
            result = subprocess.run(cmd_str, shell=True, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"{module_name} failed with return code {result.returncode}")
                return False

            if result_subdir:
                src = temp_dir / result_subdir
                if src.exists():
                    dst = self.user_output_dir / result_subdir
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.logger.info(f"Results copied to {dst}")

            if extra_files:
                for fname in extra_files:
                    src = temp_dir / fname
                    if src.exists():
                        shutil.copy2(src, self.user_output_dir / fname)
                        self.logger.info(f"Copied {fname} to output directory")

            return True

        except Exception as e:
            self.logger.error(f"Exception in {module_name}: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    # ------------------------------------------------------------------
    # Module‑specific run methods 
    # ------------------------------------------------------------------
    def run_qc_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"{sys.executable} qc_module/acineto_fasta_qc.py {pattern}"
        return self.run_module_in_temp("qc_module", fasta_files, cmd, "fasta_qc_results")

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path,
                          threads: int, scheme: str) -> bool:
        pattern_unquoted = self.get_file_pattern(fasta_files).strip('"')
        temp_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_mlst_"))
        self._register_temp_dir(temp_dir)
        try:
            module_orig = self.get_module_path("mlst_module")
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True)
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            script = temp_dir / "mlst_module.py"
            output_subdir = f"mlst_{scheme}_results"
            cmd = [sys.executable, str(script), "-i", pattern_unquoted,
                   "-o", output_subdir, "-db", "db", "-sc", "bin", "--batch", "-s", scheme]
            self.print_command(" ".join(cmd))
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"MLST {scheme} failed: {result.stderr}")
                return False

            scheme_dir = "PASTEUR_MLST" if scheme == "pasteur" else "OXFORD_MLST"
            src = temp_dir / output_subdir / scheme_dir
            if src.exists():
                dst = self.user_output_dir / scheme_dir
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
                self.logger.info(f"MLST {scheme} results copied to {dst}")
                # Do not copy the summary HTML to top-level – keep it inside the subdirectory
            else:
                self.logger.warning(f"MLST {scheme} results directory not found: {src}")
                return False
            return True
        except Exception as e:
            self.logger.error(f"MLST {scheme} exception: {e}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_kaptive_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"{sys.executable} k_o_module/acineto_kaptive.py {pattern}"
        return self.run_module_in_temp("k_o_module", fasta_files, cmd, "kaptive_results")

    def run_amr_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int,
                         min_identity: float = None, min_coverage: float = None,
                         skip_mutations: bool = False, force_update: bool = False) -> bool:
        cmd = f"{sys.executable} amr_module/acineto_amrfinder.py " + self.get_file_pattern(fasta_files)
        if min_identity is not None:
            cmd += f" --min-identity {min_identity}"
        if min_coverage is not None:
            cmd += f" --min-coverage {min_coverage}"
        if skip_mutations:
            cmd += " --skip-mutations"
        if force_update:
            self.print_info("Forcing AMR database update before analysis...")
            self.update_amr_database(force=True)

        return self.run_module_in_temp(
            "amr_module", fasta_files, cmd,
            "acineto_amrfinder_results",
            extra_files=["mutation_summary.tsv", "mutation_summary.html"]
        )

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int,
                              min_identity: float = None, min_coverage: float = None) -> bool:
        cmd = f"{sys.executable} abricate_module/acineto_abricate.py " + self.get_file_pattern(fasta_files)
        if min_identity is not None:
            cmd += f" --min-identity {min_identity}"
        if min_coverage is not None:
            cmd += f" --min-coverage {min_coverage}"
        return self.run_module_in_temp("abricate_module", fasta_files, cmd, "acineto_abricate_results")

    def run_plasmid_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"{sys.executable} plasmid_typing_module/plasmid_typer.py {pattern}"
        return self.run_module_in_temp("plasmid_typing_module", fasta_files, cmd, "acineto_plasmid_results")

    # ------------------------------------------------------------------
    # AMR Database update
    # ------------------------------------------------------------------
    def update_amr_database(self, force: bool = False) -> bool:
        amr_module_path = self.get_module_path("amr_module")
        amr_script = amr_module_path / "acineto_amrfinder.py"
        if not amr_script.exists():
            self.print_error(f"AMR script not found at: {amr_script}")
            return False

        # Ensure logger exists for standalone commands
        if self.logger is None:
            import logging
            logging.basicConfig(level=logging.INFO, format='%(message)s')
            self.logger = logging.getLogger("AcinetoScope")

        self.print_info("Updating AMRfinderPlus database...")
        flag = "--force-update" if force else "--update-db"
        cmd = [sys.executable, str(amr_script), flag]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.stdout:
            self.logger.info(result.stdout)
        if result.stderr:
            self.logger.warning(result.stderr)
        if result.returncode == 0:
            self.print_success("AMR database updated successfully.")
            version_cmd = [sys.executable, str(amr_script), "--db-version"]
            version_result = subprocess.run(version_cmd, capture_output=True, text=True, cwd=amr_module_path)
            if version_result.returncode == 0:
                self.print_info(f"New database version: {version_result.stdout.strip()}")
            return True
        else:
            self.print_error("AMR database update failed.")
            if result.stderr:
                self.logger.warning(result.stderr)
            return False

    def ensure_amr_database(self) -> bool:
        amr_module_path = self.get_module_path("amr_module")
        amr_script = amr_module_path / "acineto_amrfinder.py"
        if not amr_script.exists():
            self.print_error("AMR script not found, cannot check database.")
            return False
        cmd = [sys.executable, str(amr_script), "--db-version"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.stdout and "Unknown" not in result.stdout and "No database" not in result.stdout:
            self.print_success(f"AMR database already present: {result.stdout.strip()}")
            return True
        else:
            self.print_warning("AMR database not found or outdated. Attempting automatic update...")
            return self.update_amr_database(force=False)

    # ------------------------------------------------------------------
    # Gene‑centric and sample‑centric reporters
    # ------------------------------------------------------------------
    def run_gene_centric_reporter(self, output_dir: Path) -> bool:
        module_name = "gene_centric_module"
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.print_error(f"Gene‑centric module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_gene_centric_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for gene‑centric: {temp_dir}")

        try:
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True)

            html_files = {
                "FASTA_QC_summary.html": "fasta_qc_results",
                "pasteur_mlst_summary.html": "PASTEUR_MLST",
                "oxford_mlst_summary.html": "OXFORD_MLST",
                "Kaptive_summary.html": "kaptive_results",
                "acineto_amrfinder_summary_report.html": "acineto_amrfinder_results",
                "acineto_card_summary_report.html": "acineto_abricate_results",
                "acineto_ncbi_summary_report.html": "acineto_abricate_results",
                "acineto_resfinder_summary_report.html": "acineto_abricate_results",
                "acineto_vfdb_summary_report.html": "acineto_abricate_results",
                "acineto_argannot_summary_report.html": "acineto_abricate_results",
                "acineto_megares_summary_report.html": "acineto_abricate_results",
                "acineto_ecoli_vf_summary_report.html": "acineto_abricate_results",
                "acineto_bacmet2_summary_report.html": "acineto_abricate_results",
                "acineto_ecoh_summary_report.html": "acineto_abricate_results",
                "acineto_plasmidfinder_summary_report.html": "acineto_abricate_results",
                "plasmid_summary_report.html": "acineto_plasmid_results",
                "mutation_summary.html": "acineto_amrfinder_results",
            }

            for fname, subdir in html_files.items():
                src = output_dir / subdir / fname
                if src.exists():
                    shutil.copy2(src, temp_dir / fname)
                    self.logger.info(f"Copied {fname} to gene‑centric temp dir")
                else:
                    self.logger.warning(f"Required HTML not found: {src}")

            self.print_info("Running gene‑centric reporter...")
            script = temp_dir / "genius_acinetobacter_gene_centric_reporter.py"
            if not script.exists():
                self.logger.error("genius_acinetobacter_gene_centric_reporter.py not found")
                return False

            cmd = [sys.executable, str(script), "-i", "."]
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"Gene‑centric reporter failed with return code {result.returncode}")
                return False

            src_dir = temp_dir / "GENIUS_ACINETOBACTER_GENE_CENTRIC_REPORTS"
            if src_dir.exists():
                dst_dir = self.user_output_dir / "GENIUS_ACINETOBACTER_GENE_CENTRIC_REPORTS"
                if dst_dir.exists():
                    shutil.rmtree(dst_dir)
                shutil.copytree(src_dir, dst_dir)
                self.logger.info(f"Gene‑centric reports copied to {dst_dir}")
                self.print_success("Gene‑centric reporter completed successfully.")
                return True
            else:
                self.logger.error("GENIUS_ACINETOBACTER_GENE_CENTRIC_REPORTS not found")
                return False
        except Exception as e:
            self.logger.error(f"Gene‑centric reporter exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info("Removed temporary gene‑centric directory")

    def run_sample_centric_reporter(self, output_dir: Path) -> bool:
        module_name = "sample_centric_module"
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.print_error(f"Sample‑centric module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_sample_centric_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for sample‑centric: {temp_dir}")

        try:
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True)

            sample_files = {
                "acineto_amrfinder_summary.tsv": "acineto_amrfinder_results",
                "acineto_card_abricate_summary.tsv": "acineto_abricate_results",
                "acineto_ncbi_abricate_summary.tsv": "acineto_abricate_results",
                "acineto_resfinder_abricate_summary.tsv": "acineto_abricate_results",
                "acineto_vfdb_abricate_summary.tsv": "acineto_abricate_results",
                "acineto_argannot_abricate_summary.tsv": "acineto_abricate_results",
                "acineto_megares_abricate_summary.tsv": "acineto_abricate_results",
                "acineto_ecoli_vf_abricate_summary.tsv": "acineto_abricate_results",
                "acineto_bacmet2_abricate_summary.tsv": "acineto_abricate_results",
                "acineto_ecoh_abricate_summary.tsv": "acineto_abricate_results",
                "acineto_plasmidfinder_abricate_summary.tsv": "acineto_abricate_results",
                "plasmid_summary.tsv": "acineto_plasmid_results",
                "mutation_summary.tsv": "acineto_amrfinder_results",
                "FASTA_QC_summary.html": "fasta_qc_results",
                "pasteur_mlst_summary.html": "PASTEUR_MLST",
                "oxford_mlst_summary.html": "OXFORD_MLST",
                "Kaptive_summary.html": "kaptive_results",
            }

            for fname, subdir in sample_files.items():
                src = output_dir / subdir / fname
                if src.exists():
                    shutil.copy2(src, temp_dir / fname)
                    self.logger.info(f"Copied {fname} to sample‑centric temp dir")
                else:
                    self.logger.warning(f"Required file not found: {src}")

            self.print_info("Running sample‑centric reporter...")
            script = temp_dir / "genius_acinetobacter_sample_centric_reporter.py"
            if not script.exists():
                self.logger.error("genius_acinetobacter_sample_centric_reporter.py not found")
                return False

            cmd = [sys.executable, str(script), "-i", "."]
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"Sample‑centric reporter failed with return code {result.returncode}")
                return False

            src_dir = temp_dir / "GENIUS_ACINETOBACTER_SAMPLE_CENTRIC_REPORTS"
            if src_dir.exists():
                dst_dir = self.user_output_dir / "GENIUS_ACINETOBACTER_SAMPLE_CENTRIC_REPORTS"
                if dst_dir.exists():
                    shutil.rmtree(dst_dir)
                shutil.copytree(src_dir, dst_dir)
                self.logger.info(f"Sample‑centric reports copied to {dst_dir}")
                self.print_success("Sample‑centric reporter completed successfully.")
                return True
            else:
                self.logger.error("GENIUS_ACINETOBACTER_SAMPLE_CENTRIC_REPORTS not found")
                return False
        except Exception as e:
            self.logger.error(f"Sample‑centric reporter exception: {e}\n{traceback.format_exc()}")
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info("Removed temporary sample‑centric directory")

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------
    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1,
                              skip_modules: Dict[str, bool] = None,
                              amr_min_identity: float = None,
                              amr_min_coverage: float = None,
                              amr_skip_mutations: bool = False,
                              amr_force_update: bool = False,
                              abricate_min_identity: float = None,
                              abricate_min_coverage: float = None,
                              mlst_scheme: str = "both",
                              skip_sample_centric: bool = False,
                              clean_output: bool = False):
        if skip_modules is None:
            skip_modules = {}

        start_time = datetime.now()
        try:
            self.display_banner()
            output_path = Path(output_dir)
            if clean_output and output_path.exists():
                shutil.rmtree(output_path)
            output_path.mkdir(parents=True, exist_ok=True)
            self.setup_logging(output_path)

            fasta_files = self.find_fasta_files(input_path)
            if not fasta_files:
                self.print_error("No FASTA files found! Analysis stopped.")
                return

            extensions = set(f.suffix.lower() for f in fasta_files)
            self.print_success(f"Starting analysis of {len(fasta_files)} samples")
            self.print_info(f"File formats: {', '.join(extensions)}")

            # Create output subdirectories
            subdirs = ["fasta_qc_results", "PASTEUR_MLST", "OXFORD_MLST", "kaptive_results",
                       "acineto_amrfinder_results", "acineto_abricate_results", "acineto_plasmid_results"]
            for subdir in subdirs:
                (output_path / subdir).mkdir(exist_ok=True)

            # Display analysis plan
            self.print_header("Analysis Plan", "Modules to be executed")
            analyses = [
                ("QC", not skip_modules.get('qc', False)),
                ("MLST Pasteur", not skip_modules.get('mlst', False) and mlst_scheme in ["pasteur", "both"]),
                ("MLST Oxford", not skip_modules.get('mlst', False) and mlst_scheme in ["oxford", "both"]),
                ("Kaptive", not skip_modules.get('kaptive', False)),
                ("AMRfinderPlus", not skip_modules.get('amr', False)),
                ("ABRicate", not skip_modules.get('abricate', False)),
                ("Plasmid Typing (APT)", not skip_modules.get('plasmid', False)),
                ("Gene‑centric Reporter", not skip_modules.get('gene_centric', False)),
                ("Sample‑centric Reporter", not skip_sample_centric and not skip_modules.get('sample_centric', False)),
            ]
            for name, enabled in analyses:
                status = "✅ ENABLED" if enabled else "⏸️  SKIPPED"
                print(f"   {status} - {name}")
            sys.stdout.flush()

            # Ensure AMR database if not skipped
            if not skip_modules.get('amr', False):
                if not self.ensure_amr_database():
                    self.print_warning("AMR database check failed; AMR analysis may fail.")

            # Run each module sequentially with full headers
            module_configs = [
                ("QC Analysis", "Quality Control and Statistics", not skip_modules.get('qc', False),
                 lambda: self.run_qc_analysis(fasta_files, output_path, threads)),
                ("MLST Analysis - PASTEUR", "Multi‑Locus Sequence Typing",
                 not skip_modules.get('mlst', False) and mlst_scheme in ["pasteur", "both"],
                 lambda: self.run_mlst_analysis(fasta_files, output_path, threads, "pasteur")),
                ("MLST Analysis - OXFORD", "Multi‑Locus Sequence Typing",
                 not skip_modules.get('mlst', False) and mlst_scheme in ["oxford", "both"],
                 lambda: self.run_mlst_analysis(fasta_files, output_path, threads, "oxford")),
                ("Kaptive Analysis", "K and O Locus Typing", not skip_modules.get('kaptive', False),
                 lambda: self.run_kaptive_analysis(fasta_files, output_path, threads)),
                ("AMR Analysis", "Antimicrobial Resistance Gene Detection", not skip_modules.get('amr', False),
                 lambda: self.run_amr_analysis(fasta_files, output_path, threads,
                                                amr_min_identity, amr_min_coverage,
                                                amr_skip_mutations, amr_force_update)),
                ("ABRicate Analysis", "Comprehensive Resistance & Virulence Screening", not skip_modules.get('abricate', False),
                 lambda: self.run_abricate_analysis(fasta_files, output_path, threads,
                                                     abricate_min_identity, abricate_min_coverage)),
                ("Plasmid Typing (APT)", "Acinetobacter Plasmid Typing", not skip_modules.get('plasmid', False),
                 lambda: self.run_plasmid_analysis(fasta_files, output_path, threads)),
            ]

            for title, subtitle, enabled, func in module_configs:
                if enabled:
                    self.print_header(title, subtitle)
                    success = func()
                    if success:
                        self.print_success(f"{title} completed successfully.")
                    else:
                        self.print_warning(f"{title} had issues (check log).")
                    # Display a random quote after each module
                    self.display_random_quote()

            # Now run reporters
            self.print_header("Generating Reports", "Ultimate Reporters")
            if not skip_modules.get('gene_centric', False):
                self.print_header("Gene‑centric Reporter", "Ultimate gene‑centric integrated report")
                gene_success = self.run_gene_centric_reporter(output_path)
                if gene_success:
                    self.print_success("Gene‑centric reporter completed successfully.")
                else:
                    self.print_warning("Gene‑centric reporter had issues (check log).")
                self.display_random_quote()

            if not skip_sample_centric and not skip_modules.get('sample_centric', False):
                self.print_header("Sample‑centric Reporter", "Interactive isolate‑centric report")
                sample_success = self.run_sample_centric_reporter(output_path)
                if sample_success:
                    self.print_success("Sample‑centric reporter completed successfully.")
                else:
                    self.print_warning("Sample‑centric reporter had issues (check log).")
                self.display_random_quote()

            # Finalise output: copy ultimate reports to AcinetoScope_final_reports, then remove duplicates
            final_dir = output_path / "AcinetoScope_final_reports"
            final_dir.mkdir(parents=True, exist_ok=True)

            # Copy ultimate report directories and remove top-level duplicates
            for src_name, dst_name in [
                ("GENIUS_ACINETOBACTER_GENE_CENTRIC_REPORTS", "GENIUS_ACINETOBACTER_GENE_CENTRIC_REPORTS"),
                ("GENIUS_ACINETOBACTER_SAMPLE_CENTRIC_REPORTS", "GENIUS_ACINETOBACTER_SAMPLE_CENTRIC_REPORTS")
            ]:
                src = output_path / src_name
                if src.exists():
                    dst = final_dir / dst_name
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.logger.info(f"Copied {src_name} to {dst}")
                    # Remove top-level duplicate
                    shutil.rmtree(src)
                    self.logger.info(f"Removed top-level duplicate: {src_name}")

            # Copy log file to final reports, then remove top-level log
            log_src = output_path / "acinetoscope_run.log"
            if log_src.exists():
                shutil.copy2(log_src, final_dir / "acinetoscope_run.log")
                log_src.unlink()
                self.logger.info("Copied log to final reports and removed top-level copy")

            analysis_time = datetime.now() - start_time
            analysis_time_str = str(analysis_time).split('.')[0]

            self.print_header("Analysis Complete", f"Time elapsed: {analysis_time_str}")
            self.print_success(f"🎉 Output ready in: {output_path}")
            self.print_info(f"Final reports in: {final_dir}")
            self.print_info("All intermediate module results are preserved.")
            self.display_random_quote()

        except KeyboardInterrupt:
            self.print_error("Analysis interrupted by user")
        except Exception as e:
            self.print_error(f"Critical error: {str(e)}")
            import traceback
            traceback.print_exc()

    # ------------------------------------------------------------------
    # File finding helper
    # ------------------------------------------------------------------
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


# ----------------------------------------------------------------------
# Main entry point
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="AcinetoScope: Complete A. baumannii Typing Pipeline",
        formatter_class=ColoredHelpFormatter,
        epilog=f"""
{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Examples:{ColoredHelpFormatter.RESET}
  # Basic analysis
  acinetoscope -i genome.fna -o results/

  # Batch with glob pattern
  acinetoscope -i "*.fna" -o batch_results --threads 8

  # Skip some modules
  acinetoscope -i "*.fasta" -o analysis --threads 16 --skip-abricate --skip-plasmid

  # AMR with custom thresholds and no mutation reporting
  acinetoscope -i "*.fna" -o results --amr-min-identity 0.95 --amr-min-coverage 0.9 --skip-amr-mutations

  # Force update AMR database before analysis
  acinetoscope -i "*.fna" -o results --amr-force-update

  # Update AMR database only (standalone)
  acinetoscope --update-amr-db

  # Force update AMR database only
  acinetoscope --force-update-amr-db

  # ABRicate with custom thresholds
  acinetoscope -i "*.fna" -o results --abricate-min-identity 90 --abricate-min-coverage 85

  # Skip sample‑centric reporter
  acinetoscope -i "*.fna" -o results --skip-sample-centric

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Supported FASTA formats:{ColoredHelpFormatter.RESET} .fna, .fasta, .fa, .fn

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Analysis Modules:{ColoredHelpFormatter.RESET}
  • FASTA QC (Quality Control & Statistics)
  • MLST (Pasteur & Oxford schemes)
  • Kaptive (K and O locus typing)
  • AMR profiling (AMRFinderPlus) – point mutations reported by default
  • ABRicate (Comprehensive resistance/plasmid/virulence)
  • Plasmid typing (APT – Acinetobacter Plasmid Typing)
  • Gene‑centric ultimate reporter
  • Sample‑centric ultimate reporter (isolate‑centric)

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Output:{ColoredHelpFormatter.RESET}
  All intermediate results are preserved in the output directory.
  Final reports are placed in AcinetoScope_final_reports/ containing:
    • GENIUS_ACINETOBACTER_GENE_CENTRIC_REPORTS/
    • GENIUS_ACINETOBACTER_SAMPLE_CENTRIC_REPORTS/
  A detailed log file is included in the final reports directory.

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Support & Contributions:{ColoredHelpFormatter.RESET}
  • Issues & feature requests: https://github.com/bbeckley-hub/acinetoscope/issues
  • Email: brownbeckley94@gmail.com

{ColoredHelpFormatter.YELLOW}⭐ Star us on GitHub if you find this tool useful! ⭐{ColoredHelpFormatter.RESET}
        """
    )

    # Required args
    parser.add_argument('-i', '--input', help='Input FASTA file(s) - can use glob patterns')
    parser.add_argument('-o', '--output', help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2, help='Number of threads (default: 2)')

    # Skip flags
    parser.add_argument('--skip-qc', action='store_true', help='Skip QC analysis')
    parser.add_argument('--skip-mlst', action='store_true', help='Skip MLST analysis')
    parser.add_argument('--skip-kaptive', action='store_true', help='Skip Kaptive analysis')
    parser.add_argument('--skip-amr', action='store_true', help='Skip AMR analysis (AMRfinderPlus)')
    parser.add_argument('--skip-abricate', action='store_true', help='Skip ABRicate analysis')
    parser.add_argument('--skip-plasmid', action='store_true', help='Skip plasmid typing (APT) analysis')
    parser.add_argument('--skip-gene-centric', action='store_true', help='Skip gene‑centric ultimate reporter')
    parser.add_argument('--skip-sample-centric', action='store_true', help='Skip sample‑centric ultimate reporter')

    # MLST scheme
    parser.add_argument('--mlst-scheme', choices=['pasteur', 'oxford', 'both'], default='both',
                        help='MLST scheme to use (default: both)')

    # AMR flags
    parser.add_argument('--amr-min-identity', type=float, help='Minimum identity for AMR hits (0..1)')
    parser.add_argument('--amr-min-coverage', type=float, help='Minimum coverage for AMR hits (0..1)')
    parser.add_argument('--skip-amr-mutations', action='store_true', help='Disable point mutation reporting in AMR')
    parser.add_argument('--amr-force-update', action='store_true', help='Force update AMR database before analysis')
    parser.add_argument('--update-amr-db', action='store_true', help='Update AMR database (incremental) and exit')
    parser.add_argument('--force-update-amr-db', action='store_true', help='Force complete AMR database update and exit')

    # ABRicate flags
    parser.add_argument('--abricate-min-identity', type=float, help='Minimum identity for ABRicate hits (default: 80)')
    parser.add_argument('--abricate-min-coverage', type=float, help='Minimum coverage for ABRicate hits (default: 80)')

    # Misc flags
    parser.add_argument('--keep-temp', action='store_true', help='Do not delete temporary directories (debugging)')
    parser.add_argument('--clean-output', action='store_true', help='Delete output directory before analysis')

    args = parser.parse_args()

    # Handle AMR database update commands
    if args.update_amr_db or args.force_update_amr_db:
        orch = AcinetoScopeOrchestrator()
        if args.force_update_amr_db:
            orch.update_amr_database(force=True)
        else:
            orch.update_amr_database(force=False)
        sys.exit(0)

    if not args.input or not args.output:
        parser.error("Both -i/--input and -o/--output are required (unless using database update flags).")

    skip_modules = {
        'qc': args.skip_qc,
        'mlst': args.skip_mlst,
        'kaptive': args.skip_kaptive,
        'amr': args.skip_amr,
        'abricate': args.skip_abricate,
        'plasmid': args.skip_plasmid,
        'gene_centric': args.skip_gene_centric,
        'sample_centric': args.skip_sample_centric,
    }

    orch = AcinetoScopeOrchestrator()
    orch.keep_temp = args.keep_temp
    orch.run_complete_analysis(
        input_path=args.input,
        output_dir=args.output,
        threads=args.threads,
        skip_modules=skip_modules,
        amr_min_identity=args.amr_min_identity,
        amr_min_coverage=args.amr_min_coverage,
        amr_skip_mutations=args.skip_amr_mutations,
        amr_force_update=args.amr_force_update,
        abricate_min_identity=args.abricate_min_identity,
        abricate_min_coverage=args.abricate_min_coverage,
        mlst_scheme=args.mlst_scheme,
        skip_sample_centric=args.skip_sample_centric,
        clean_output=args.clean_output
    )


if __name__ == "__main__":
    main()