#!/usr/bin/env python3
"""
AcinetoScope Main Orchestrator - v1.4.0
All module writes happen in /tmp, final results are copied to user output.
HPC / Docker-friendly with full flag support.

Changelog (1.4.0):
  • NEW: MGE module (mobileOG-db) integrated into the pipeline
  • NEW: --pull-mge-db downloads mobileOG-db to a writable location
    (conda users only; Docker/Singularity already bundle the data)
  • NEW: Species gate (FASTA_QC_species_passed.txt) filters downstream modules
  • NEW: --skip-mge flag
  • NEW: Citations printed at end of every run + in --help epilog
  • NEW: Simple ASCII banner printed at startup
  • IMPROVED: --pull-mlst-db prunes to A. baumannii schemes only
  • IMPROVED: MLST output path auto-detects both old and new folder names
  • IMPROVED: Gene-centric + sample-centric reporters now receive MGE summaries
  • Folder names consistent: abricate_results, amrfinder_results, plasmid_results,
    mge_results, pasteur_mlst_results, oxford_mlst_results, kaptive_results

Changelog (1.3.3):
  • Add --pull-mlst-db / --update-mlst-db (PubMLST + GitHub)
  • Add --report-level {final,summary,full} (disk footprint control)
  • Add AcinetoScope_final_report/ consolidated folder
  • Clean folder names: no more acineto_* prefixes
  • HPC / Docker env-var overrides (ACINETOSCOPE_*) for tool binaries + DBs
  • Per-sample copy churn eliminated

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-08-22 / 2026-10-24 / 2026-10-29 / 2026-11-04
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
import urllib.request
import zipfile
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Set

__version__ = "1.4.0"


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
# Canonical output directory names (user-facing)
# ----------------------------------------------------------------------
DIR_FASTA_QC = "fasta_qc_results"
DIR_MLST_PASTEUR = "pasteur_mlst_results"
DIR_MLST_OXFORD = "oxford_mlst_results"
DIR_KAPTIVE = "kaptive_results"
DIR_AMRFINDER = "amrfinder_results"
DIR_ABRICATE = "abricate_results"
DIR_PLASMID = "plasmid_results"
DIR_MGE = "mge_results"
DIR_FINAL = "AcinetoScope_final_report"
DIR_GENE_CENTRIC = "GENIUS_ACINETOBACTER_GENE_CENTRIC_REPORTS"
DIR_SAMPLE_CENTRIC = "GENIUS_ACINETOBACTER_SAMPLE_CENTRIC_REPORTS"

ALL_MODULE_DIRS = [
    DIR_FASTA_QC, DIR_MLST_PASTEUR, DIR_MLST_OXFORD, DIR_KAPTIVE,
    DIR_AMRFINDER, DIR_ABRICATE, DIR_PLASMID, DIR_MGE,
]


CITATION_LINES = [
    "Beckley, B., Amarh, V., Lopes, B.S. et al. AcinetoScope: a computational pipeline for rapid, comprehensive Acinetobacter baumannii outbreak investigation and resistance gene tracking. Sci Rep (2026). https://doi.org/10.1038/s41598-026-70189-y",
]


# Simple ASCII banner printed once at the top of every pipeline run.
ASCII_ART_TOP = r"""
################################################################################
#                                                                              #
#             A   C   I   N   E   T   O   S   C   O   P   E                    #
#                                                                              #
#                  A. baumannii  Genomic  Analysis  Pipeline                   #
#                                                                              #
################################################################################
"""


# ----------------------------------------------------------------------
# Main Orchestrator
# ----------------------------------------------------------------------
class AcinetoScopeOrchestrator:
    """HPC/Docker-friendly orchestrator for AcinetoScope."""

    def __init__(self):
        self.base_dir = Path(__file__).parent
        self.user_output_dir = None
        self.logger = None
        self.keep_temp = False
        self.report_level = 'summary'
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
    # Colored printing
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
        parts = command.split()
        if len(parts) > 1 and 'python' in parts[0]:
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
        print(ASCII_ART_TOP)
        banner = f"""{Color.BOLD}{Color.BRIGHT_MAGENTA}
{'='*80}
{' '*20}🦠 ACINETOSCOPE v{__version__} - A. baumannii Genomic Analysis Pipeline
{'='*80}
{Color.RESET}{Color.BRIGHT_CYAN}
Complete A. baumannii genomic analysis pipeline
QC | MLST | K/O Locus | AMR | Virulence | Plasmid (APT) | MGE | Summary Reports
{Color.RESET}{Color.DIM}
Author: Brown Beckley | Email: brownbeckley94@gmail.com
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
{'='*80}{Color.RESET}
"""
        print(banner)

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
    # Path resolution (HPC / Docker friendly)
    # ------------------------------------------------------------------
    def get_module_path(self, module_name: str) -> Path:
        """Return module dir with env override + packaged + dev fallback."""
        env_override = os.environ.get(f'ACINETOSCOPE_{module_name.upper()}_MODULE_PATH')
        if env_override:
            return Path(env_override)
        if hasattr(sys, 'prefix'):
            share_path = Path(sys.prefix) / "share" / "acinetoscope" / "modules" / module_name
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def get_shared_data_root(self) -> Path:
        """Writable user-local data root for read-only module installs."""
        root = Path(os.environ.get(
            "ACINETOSCOPE_SHARED_DATA",
            Path.home() / ".local" / "share" / "acinetoscope",
        ))
        root.mkdir(parents=True, exist_ok=True)
        return root

    def get_mlst_db_dir(self, for_update: bool = False) -> Path:
        """Return the MLST database directory.

        1. $ACINETOSCOPE_MLST_DB env override
        2. User-local ~/.local/share/acinetoscope/mlst_module/db/ if populated
        3. Packaged <module>/db/
        """
        env_db = os.environ.get("ACINETOSCOPE_MLST_DB")
        if env_db:
            env_path = Path(env_db)
            if for_update:
                env_path.mkdir(parents=True, exist_ok=True)
                return env_path
            if self._mlst_db_has_scheme(env_path):
                return env_path

        user_module = Path(os.environ.get(
            "ACINETOSCOPE_MLST_MODULE_PATH",
            self.get_shared_data_root() / "mlst_module",
        ))
        user_db = user_module / "db"

        if for_update:
            user_db.mkdir(parents=True, exist_ok=True)
            return user_db

        if self._mlst_db_has_scheme(user_db):
            if self.logger:
                self.logger.info(f"Using user-local MLST database: {user_db}")
            return user_db

        module_db = self.get_module_path("mlst_module") / "db"
        if module_db.exists():
            return module_db

        self.print_warning(
            "No populated MLST database found. Run `acinetoscope --pull-mlst-db`."
        )
        user_db.mkdir(parents=True, exist_ok=True)
        return user_db

    @staticmethod
    def _mlst_db_has_scheme(db_dir: Path) -> bool:
        """Check whether MLST db dir has a usable A. baumannii scheme."""
        for scheme in ("abaumannii", "abaumannii_2", "pasteur", "oxford"):
            if (db_dir / "pubmlst" / scheme).exists():
                return True
        return False

    def get_amrfinder_paths(self) -> tuple:
        """Resolve AMRFinderPlus bin + db directories (env override or bundled)."""
        env_bin = os.environ.get("ACINETOSCOPE_AMRFINDER_BIN_DIR")
        if env_bin and Path(env_bin).is_dir():
            bin_dir = Path(env_bin)
        else:
            bundled_bin = self.get_module_path("amr_module") / "bin"
            bin_dir = bundled_bin if bundled_bin.is_dir() else None

        env_db = os.environ.get("ACINETOSCOPE_AMRFINDER_DB")
        if env_db and Path(env_db).is_dir():
            db_dir = Path(env_db)
        else:
            bundled_db_root = self.get_module_path("amr_module") / "data" / "amrfinder_db"
            db_dir = None
            if bundled_db_root.is_dir():
                candidates = sorted([
                    d for d in bundled_db_root.iterdir()
                    if d.is_dir() and d.name.startswith("20")
                ])
                if candidates:
                    db_dir = candidates[-1]

        return bin_dir, db_dir

    # ------------------------------------------------------------------
    # MGE database resolution + management
    # ------------------------------------------------------------------
    def get_mge_data_dir(self, for_update: bool = False) -> Optional[Path]:
        """Return the mobileOG-db data directory.

        Resolution order:
          1. $ACINETOSCOPE_MGE_DATA_DIR env override
          2. User-local ~/.local/share/acinetoscope/mge_module/mobileOG-db/
             (populated by `acinetoscope --pull-mge-db`)
          3. Packaged <module>/mobileOG-db/ (Docker / Singularity)
        """
        env_dir = os.environ.get("ACINETOSCOPE_MGE_DATA_DIR")
        if env_dir:
            env_path = Path(env_dir)
            if for_update:
                env_path.mkdir(parents=True, exist_ok=True)
                return env_path
            if self._mge_db_has_data(env_path):
                return env_path

        user_dir = self.get_shared_data_root() / "mge_module" / "mobileOG-db"
        if for_update:
            user_dir.mkdir(parents=True, exist_ok=True)
            return user_dir
        if self._mge_db_has_data(user_dir):
            if self.logger:
                self.logger.info(f"Using user-local mobileOG-db: {user_dir}")
            return user_dir

        module_db = self.get_module_path("mge_module") / "mobileOG-db"
        if self._mge_db_has_data(module_db):
            return module_db

        return None

    @staticmethod
    def _mge_db_has_data(db_dir: Optional[Path]) -> bool:
        """Return True if the directory contains the beatrix DIAMOND index."""
        if db_dir is None or not db_dir.exists():
            return False
        candidates = [
            db_dir / "beatrix-1-6_v1_all" / "mobileOG-db-beatrix-1.6.dmnd",
            db_dir / "mobileOG-db-beatrix-1.6.dmnd",
            db_dir / "beatrix-1-6_v1_all.dmnd",
        ]
        return any(c.exists() for c in candidates)

    def pull_mge_database(self) -> bool:
        """Download the mobileOG-db payload from the GitHub release.

        Docker and Singularity users do not need this — the data is baked
        into the image and get_mge_data_dir() finds it there automatically.
        Conda users run `acinetoscope --pull-mge-db` once after install.
        """
        module_path = self.get_module_path("mge_module")

        if module_path.exists() and os.access(module_path, os.W_OK):
            target_root = module_path
            self.print_info(f"Using writable module directory: {target_root}")
        else:
            target_root = self.get_shared_data_root() / "mge_module"
            target_root.mkdir(parents=True, exist_ok=True)
            self.print_info(
                f"Module dir is read-only; using user-local path: {target_root}"
            )

        repo_url = os.environ.get(
            "ACINETOSCOPE_MGE_PULL_URL",
            "https://github.com/bbeckley-hub/acinetoscope/releases/download/"
            "v1.4.0/acinetoscope-mge-data-v1.4.0.zip",
        )
        temp_zip = Path(tempfile.mktemp(suffix=".zip"))
        extract_dir = None

        try:
            self.print_info(f"Downloading mobileOG-db from {repo_url}")
            self.print_info(
                "Payload is approximately 300 MB, this may take a few minutes..."
            )
            urllib.request.urlretrieve(repo_url, temp_zip)
            self.print_success("Download complete")

            extract_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_mge_extract_"))
            with zipfile.ZipFile(temp_zip, "r") as zip_ref:
                zip_ref.extractall(extract_dir)
            self.print_info(f"Extracted to {extract_dir}")

            beatrix_dir = None
            for root, dirs, _ in os.walk(extract_dir):
                if "beatrix-1-6_v1_all" in dirs:
                    beatrix_dir = Path(root) / "beatrix-1-6_v1_all"
                    break

            if beatrix_dir is None:
                self.print_error(
                    "Could not find 'beatrix-1-6_v1_all' directory "
                    "in downloaded zip"
                )
                return False

            dmnd = beatrix_dir / "mobileOG-db-beatrix-1.6.dmnd"
            if not dmnd.exists():
                self.print_error(
                    f"DIAMOND index missing from archive: {dmnd.name}"
                )
                return False

            target_beatrix = target_root / "mobileOG-db" / "beatrix-1-6_v1_all"
            if target_beatrix.exists():
                shutil.rmtree(target_beatrix)
            target_beatrix.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(beatrix_dir, target_beatrix)

            dmnd_target = target_beatrix / "mobileOG-db-beatrix-1.6.dmnd"
            if not dmnd_target.exists():
                self.print_error(
                    f"Install verification failed: {dmnd_target} not present"
                )
                return False

            size_mb = dmnd_target.stat().st_size / (1024 * 1024)
            self.print_success(
                f"mobileOG-db installed to {target_root / 'mobileOG-db'} "
                f"(DIAMOND index: {size_mb:.0f} MB)"
            )

            os.environ["ACINETOSCOPE_MGE_DATA_DIR"] = str(
                target_root / "mobileOG-db"
            )
            os.environ["ACINETOSCOPE_MGE_MODULE_PATH"] = str(target_root)
            return True

        except Exception as e:
            self.print_error(f"Failed to download/extract mobileOG-db: {e}")
            return False
        finally:
            if temp_zip.exists():
                temp_zip.unlink()
            if extract_dir and extract_dir.exists():
                shutil.rmtree(extract_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # File helpers
    # ------------------------------------------------------------------
    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        if not fasta_files:
            return '"*.fna"'
        extensions = set(f.suffix.lower() for f in fasta_files)
        if len(extensions) == 1:
            ext = list(extensions)[0]
            return f'"*{ext}"'
        return '"*"'

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
            patterns = [f"{input_path}/*.fna", f"{input_path}/*.fasta",
                        f"{input_path}/*.fa", f"{input_path}/*.fn"]
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

    # ------------------------------------------------------------------
    # Species gate (mirrors EcoliTyper / FASTA_QC_species_passed.txt)
    # ------------------------------------------------------------------
    def _load_species_gate(self, output_dir: Path) -> Optional[Set[str]]:
        """Read the species-gate pass list produced by FASTA QC.

        Returns None when the gate file does not exist (fail-open: downstream
        modules will run on every sample).
        """
        gate_file = output_dir / DIR_FASTA_QC / "FASTA_QC_species_passed.txt"
        if not gate_file.exists():
            return None
        with open(gate_file, "r", encoding="utf-8") as fh:
            names = {line.strip() for line in fh if line.strip()}
        if self.logger:
            self.logger.info(
                f"Species gate loaded from {gate_file}: {len(names)} sample(s) passed"
            )
        return names

    def _filter_by_species_gate(self, fasta_files: List[Path],
                                gate: Optional[Set[str]]) -> List[Path]:
        """Keep only FASTA files that passed the species gate (A. baumannii)."""
        if gate is None:
            self.print_warning(
                "Species gate not available (FASTA QC failed or was skipped). "
                "Running ALL samples through downstream modules."
            )
            return fasta_files

        passed = [f for f in fasta_files if f.name in gate]
        excluded = [f for f in fasta_files if f.name not in gate]

        if excluded:
            self.print_info(
                f"Species gate: {len(passed)} of {len(fasta_files)} samples passed"
            )
            for f in excluded:
                self.print_warning(
                    f"  Excluded (not A. baumannii): {f.name}"
                )
        else:
            self.print_success(
                f"Species gate: all {len(passed)} sample(s) passed species check"
            )
        return passed

    # ------------------------------------------------------------------
    # Copy helper
    # ------------------------------------------------------------------
    def _copy_module_results(self, src: Path, dst: Path,
                             skip_dirs: Optional[Set[str]] = None) -> None:
        """Copy module results to user output, respecting report_level.

        Batch-level files (top-level files in src) are always copied.
        Per-sample folders (subdirectories) are only copied when
        report_level == 'full'. This avoids copy-then-delete churn.
        """
        if not src.exists():
            return
        if dst.exists():
            shutil.rmtree(dst)
        dst.mkdir(parents=True, exist_ok=True)

        skip_dirs = skip_dirs or set()
        for item in src.iterdir():
            if item.name in skip_dirs:
                continue
            if item.is_file():
                shutil.copy2(item, dst / item.name)
            elif item.is_dir() and self.report_level == 'full':
                shutil.copytree(item, dst / item.name)

    # ------------------------------------------------------------------
    # Generic module runner
    # ------------------------------------------------------------------
    def run_module_in_temp(self, module_name: str, fasta_files: List[Path],
                           cmd_str: str, result_subdir: str = None,
                           target_dir: str = None,
                           extra_files: List[str] = None,
                           env: Dict[str, str] = None,
                           ignore_dirs: Optional[List[str]] = None) -> bool:
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.print_error(f"Module directory not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix=f"acinetoscope_{module_name}_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for {module_name}: {temp_dir}")

        try:
            ignore_patterns = ['__pycache__', '*.pyc'] + (ignore_dirs or [])
            shutil.copytree(module_orig, temp_dir / module_name, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns(*ignore_patterns))

            # QC module: copy ref_db to temp root so its script finds it
            ref_db_src = module_orig / "ref_db"
            if ref_db_src.exists():
                shutil.copytree(ref_db_src, temp_dir / "ref_db")
                self.logger.info(f"Copied ref_db to temp root: {temp_dir / 'ref_db'}")

            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.print_info(f"Running {module_name}...")
            self.print_command(cmd_str)
            result = subprocess.run(cmd_str, shell=True, cwd=temp_dir,
                                    capture_output=True, text=True, env=env)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(
                    f"{module_name} failed with return code {result.returncode}"
                )
                return False

            if result_subdir:
                src = temp_dir / result_subdir
                if src.exists():
                    dst_name = target_dir if target_dir else result_subdir
                    dst = self.user_output_dir / dst_name
                    self._copy_module_results(src, dst)
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
    # Module-specific runners
    # ------------------------------------------------------------------
    def run_qc_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"{sys.executable} qc_module/acineto_fasta_qc.py {pattern}"
        return self.run_module_in_temp(
            "qc_module", fasta_files, cmd,
            "fasta_qc_results", DIR_FASTA_QC,
        )

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path,
                          threads: int, scheme: str) -> bool:
        pattern_unquoted = self.get_file_pattern(fasta_files).strip('"')
        temp_dir = Path(tempfile.mkdtemp(prefix=f"acinetoscope_mlst_{scheme}_"))
        self._register_temp_dir(temp_dir)

        env = os.environ.copy()
        resolved_db = self.get_mlst_db_dir()
        if resolved_db.exists() and self._mlst_db_has_scheme(resolved_db):
            env['MLST_DATADIR'] = str(resolved_db)
            env['ACINETOSCOPE_MLST_DB'] = str(resolved_db)

        try:
            module_orig = self.get_module_path("mlst_module")
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            script = temp_dir / "mlst_module.py"
            if not script.exists():
                self.print_error("mlst_module.py not found")
                return False

            output_subdir = f"mlst_{scheme}_results"
            cmd = [sys.executable, str(script), "-i", pattern_unquoted,
                   "-o", output_subdir, "-db", "db", "-sc", "bin", "--batch", "-s", scheme]
            self.print_command(" ".join(cmd))
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True,
                                    text=True, env=env)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.logger.error(f"MLST {scheme} failed: {result.stderr}")
                return False

            # Try both old (PASTEUR_MLST) and new (pasteur_results) folder names.
            inner_candidates = [
                f"{scheme}_results",          # new (post-patch)
                f"{scheme.upper()}_MLST",     # legacy
            ]
            inner_scheme_dir = None
            for cand in inner_candidates:
                cand_path = temp_dir / output_subdir / cand
                if cand_path.is_dir():
                    inner_scheme_dir = cand
                    self.logger.info(f"MLST {scheme} wrote to {cand_path}")
                    break

            if inner_scheme_dir is None:
                self.logger.error(
                    f"MLST {scheme} output not found. Tried: {inner_candidates} "
                    f"under {temp_dir / output_subdir}"
                )
                return False

            src = temp_dir / output_subdir / inner_scheme_dir
            dst_name = DIR_MLST_PASTEUR if scheme == "pasteur" else DIR_MLST_OXFORD
            dst = self.user_output_dir / dst_name
            self._copy_module_results(src, dst)
            self.logger.info(f"MLST {scheme} results copied to {dst}")
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
        return self.run_module_in_temp(
            "k_o_module", fasta_files, cmd,
            "kaptive_results", DIR_KAPTIVE,
        )

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

        env = os.environ.copy()
        bin_dir, db_dir = self.get_amrfinder_paths()
        if db_dir:
            env['AMRFINDER_DB'] = str(db_dir)
            env['ACINETOSCOPE_AMRFINDER_DB'] = str(db_dir)
        if bin_dir:
            env['ACINETOSCOPE_AMRFINDER_BIN_DIR'] = str(bin_dir)

        return self.run_module_in_temp(
            "amr_module", fasta_files, cmd,
            "acineto_amrfinder_results", DIR_AMRFINDER,
            extra_files=["mutation_summary.tsv", "mutation_summary.html"],
            env=env,
        )

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int,
                              min_identity: float = None, min_coverage: float = None) -> bool:
        cmd = f"{sys.executable} abricate_module/acineto_abricate.py " + self.get_file_pattern(fasta_files)
        if min_identity is not None:
            cmd += f" --min-identity {min_identity}"
        if min_coverage is not None:
            cmd += f" --min-coverage {min_coverage}"
        return self.run_module_in_temp(
            "abricate_module", fasta_files, cmd,
            "acineto_abricate_results", DIR_ABRICATE,
        )

    def run_plasmid_analysis(self, fasta_files: List[Path], output_dir: Path, threads: int) -> bool:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"{sys.executable} plasmid_typing_module/plasmid_typer.py {pattern}"
        return self.run_module_in_temp(
            "plasmid_typing_module", fasta_files, cmd,
            "acineto_plasmid_results", DIR_PLASMID,
        )

    def run_mge_analysis(self, fasta_files: List[Path], output_dir: Path,
                         threads: int) -> bool:
        """Run the MGE profiler (mobileOG-db).

        Resolves the mobileOG-db payload via get_mge_data_dir():
          - Docker/Singularity: found inside the packaged mge_module/
          - Conda + `--pull-mge-db`: found at ~/.local/share/acinetoscope/
        Either way, the beatrix directory is symlinked (or copied as a
        fallback) into the temp working directory so MGEFinder.py sees the
        layout it expects:
            mobileOG-db/beatrix-1-6_v1_all/mobileOG-db-beatrix-1.6.dmnd
        """
        module_orig = self.get_module_path("mge_module")
        if not module_orig.exists():
            self.print_error(f"mge_module not found: {module_orig}")
            return False

        mge_data = self.get_mge_data_dir()
        if mge_data is None:
            self.print_error(
                "mobileOG-db not found. Run `acinetoscope --pull-mge-db` once, "
                "or use the Docker/Singularity image (which has it bundled)."
            )
            return False
        self.logger.info(f"Resolved mobileOG-db: {mge_data}")

        temp_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_mge_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for mge: {temp_dir}")

        try:
            # Copy the module tree, excluding any existing mobileOG-db —
            # the resolved copy below is authoritative.
            shutil.copytree(
                module_orig, temp_dir, dirs_exist_ok=True,
                ignore=shutil.ignore_patterns(
                    "mge_results", "annotations", "*.fna", "*.fasta",
                    "mobileOG-db", "__pycache__", "*.pyc",
                ),
            )

            # Link or copy the resolved mobileOG-db into the temp tree
            temp_mge = temp_dir / "mobileOG-db"
            temp_mge.mkdir(parents=True, exist_ok=True)

            beatrix_src = mge_data / "beatrix-1-6_v1_all"
            if not beatrix_src.exists() and mge_data.name == "beatrix-1-6_v1_all":
                beatrix_src = mge_data
            if not beatrix_src.exists():
                self.print_error(
                    f"beatrix-1-6_v1_all not found under {mge_data}"
                )
                return False

            beatrix_dst = temp_mge / "beatrix-1-6_v1_all"
            if beatrix_dst.exists() or beatrix_dst.is_symlink():
                if beatrix_dst.is_symlink() or beatrix_dst.is_file():
                    beatrix_dst.unlink()
                else:
                    shutil.rmtree(beatrix_dst)

            try:
                os.symlink(beatrix_src.resolve(), beatrix_dst)
                self.logger.info(f"Symlinked {beatrix_src} -> {beatrix_dst}")
            except OSError as exc:
                self.logger.warning(
                    f"Symlink failed ({exc}); falling back to copy"
                )
                shutil.copytree(beatrix_src, beatrix_dst)
                self.logger.info(f"Copied {beatrix_src} -> {beatrix_dst}")

            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.print_info(
                f"Running MGE analysis on {len(fasta_files)} genome(s)..."
            )
            script = temp_dir / "MGEFinder.py"
            if not script.exists():
                self.print_error("MGEFinder.py not found in mge_module")
                return False

            mge_inputs = [f.name for f in fasta_files]
            cmd = [sys.executable, str(script), "-i"] + mge_inputs + [
                "-d", "mge_results",
                "-db_dir", "mobileOG-db/beatrix-1-6_v1_all",
            ]
            display_pattern = self.get_file_pattern(fasta_files).strip('"')
            display_cmd = (
                f"{sys.executable} {script} -i {display_pattern} "
                f"-d mge_results -db_dir mobileOG-db/beatrix-1-6_v1_all"
            )
            self.print_command(display_cmd)

            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True,
                                    text=True)
            if result.stdout:
                self.logger.info(result.stdout[:2000])
            if result.stderr:
                self.logger.warning(result.stderr[:1000])

            if result.returncode != 0:
                self.print_error(
                    f"MGE analysis failed with return code {result.returncode}"
                )
                return False

            src = temp_dir / "mge_results"
            if not src.exists():
                self.print_error(f"MGE results directory not found: {src}")
                return False

            dst = self.user_output_dir / DIR_MGE
            self._copy_module_results(src, dst, skip_dirs={"annotations"})
            expected_files = [
                "acineto_mge_summary.tsv",
                "acineto_mge_all_hits.tsv",
                "acineto_mge_summary.json",
                "acineto_mge_summary.html",
            ]
            missing = [f for f in expected_files if not (dst / f).exists()]
            if missing:
                self.print_warning(
                    f"MGE ran but these summary files are missing: "
                    f"{', '.join(missing)}. Reporters may not include MGE data."
                )

            self.print_success(
                f"MGE results copied to {dst} (annotations skipped)"
            )
            return True
        except Exception as e:
            self.print_error(f"MGE exception: {e}")
            self.logger.error(traceback.format_exc())
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # AMR database management
    # ------------------------------------------------------------------
    def update_amr_database(self, force: bool = False) -> bool:
        amr_module_path = self.get_module_path("amr_module")
        amr_script = amr_module_path / "acineto_amrfinder.py"
        if not amr_script.exists():
            self.print_error(f"AMR script not found at: {amr_script}")
            return False

        if self.logger is None:
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
            version_result = subprocess.run(version_cmd, capture_output=True,
                                            text=True, cwd=amr_module_path)
            if version_result.returncode == 0:
                self.print_info(f"New database version: {version_result.stdout.strip()}")
            return True
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
        self.print_warning("AMR database not found or outdated. Attempting automatic update...")
        return self.update_amr_database(force=False)

    # ------------------------------------------------------------------
    # MLST database management
    # ------------------------------------------------------------------
    def pull_mlst_database(self) -> bool:
        """Download prebuilt MLST database from GitHub, then prune to
        A. baumannii schemes (Oxford + Pasteur) only."""
        module_path = self.get_module_path("mlst_module")
        system_path = module_path

        if module_path.exists() and os.access(module_path, os.W_OK):
            target_path = module_path
            self.print_info(f"Using writable module directory: {target_path}")
        else:
            target_path = self.get_shared_data_root() / "mlst_module"
            target_path.mkdir(parents=True, exist_ok=True)
            self.print_info(
                f"System MLST module is read-only; using user-local path: {target_path}"
            )

        repo_url = os.environ.get(
            "ACINETOSCOPE_MLST_PULL_URL",
            "https://github.com/bbeckley-hub/mlst/archive/refs/heads/master.zip",
        )
        temp_zip = Path(tempfile.mktemp(suffix=".zip"))
        extract_dir = None

        try:
            self.print_info(f"Downloading MLST data from {repo_url}...")
            urllib.request.urlretrieve(repo_url, temp_zip)

            extract_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_mlst_extract_"))
            with zipfile.ZipFile(temp_zip, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)

            src_root = extract_dir / "mlst-master"
            if not src_root.exists():
                self.print_error("Extracted repository root not found.")
                return False

            for dirname in ['bin', 'db', 'perl5', 'scripts']:
                src = src_root / dirname
                if src.exists():
                    dst = target_path / dirname
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.logger.info(f"Copied {dirname} to {dst}")
                else:
                    self.print_warning(f"'{dirname}' not found in repository, skipping.")

            if target_path != system_path:
                for script in ['mlst_module.py', 'mlst_database.sh']:
                    src_script = system_path / script
                    if src_script.exists():
                        shutil.copy2(src_script, target_path / script)
                        self.logger.info(f"Copied {script} to {target_path / script}")

            target_db = target_path / "db"

            # Prune to A. baumannii schemes only (Oxford + Pasteur).
            # AcinetoScope --mlst-scheme both needs both; everything else
            # (E. coli, Salmonella, etc.) is unnecessary bulk.
            target_pubmlst = target_db / "pubmlst"
            keep_schemes = {"abaumannii", "abaumannii_2"}
            if target_pubmlst.exists():
                pruned = 0
                for scheme_dir in list(target_pubmlst.iterdir()):
                    if scheme_dir.is_dir() and scheme_dir.name not in keep_schemes:
                        shutil.rmtree(scheme_dir, ignore_errors=True)
                        pruned += 1
                if pruned:
                    self.print_info(
                        f"Pruned {pruned} non-A. baumannii MLST scheme(s) "
                        f"from {target_pubmlst}"
                    )

            if not self._mlst_db_has_scheme(target_db):
                self.print_error(
                    f"A. baumannii MLST scheme missing after pull. "
                    f"Expected pubmlst/abaumannii or pubmlst/abaumannii_2 "
                    f"under {target_db}"
                )
                return False

            self.print_success(f"MLST module pulled to {target_path}")
            os.environ['ACINETOSCOPE_MLST_MODULE_PATH'] = str(target_path)
            os.environ['ACINETOSCOPE_MLST_DB'] = str(target_db)
            return True

        except Exception as e:
            self.print_error(f"Failed to download/extract MLST data: {e}")
            return False
        finally:
            if temp_zip.exists():
                temp_zip.unlink()
            if extract_dir and extract_dir.exists():
                shutil.rmtree(extract_dir, ignore_errors=True)

    def update_mlst_database(self) -> bool:
        """Refresh MLST scheme from PubMLST (API key required)."""
        module_path = self.get_module_path("mlst_module")
        if module_path.exists() and os.access(module_path, os.W_OK):
            target_module = module_path
        else:
            target_module = self.get_shared_data_root() / "mlst_module"

        db_dir = target_module / "db"
        script_dir = target_module / "bin"
        scheme_map = db_dir / "scheme_species_map.tab"

        if not scheme_map.exists():
            self.print_error(
                f"MLST module tree is not populated at {target_module}"
            )
            self.print_info(
                "Run `acinetoscope --pull-mlst-db` first to install the module, "
                "then run `acinetoscope --update-mlst-db` to refresh the scheme."
            )
            return False

        if not script_dir.exists():
            self.print_error(f"MLST bin/ not found at {script_dir}")
            return False

        os.environ['ACINETOSCOPE_MLST_MODULE_PATH'] = str(target_module)
        os.environ['ACINETOSCOPE_MLST_DB'] = str(db_dir)

        # Try to import the module's update function
        try:
            try:
                from acinetoScope.modules.mlst_module.mlst_module import update_mlst_database as _update
            except ImportError:
                sys.path.insert(0, str(self.base_dir / "modules" / "mlst_module"))
                from mlst_module import update_mlst_database as _update
        except (ImportError, AttributeError):
            self.print_warning(
                "The MLST module does not expose an `update_mlst_database()` helper."
            )
            self.print_info(
                "Option 1: refresh manually with `mlstdb update` "
                "(requires: mlstdb connect --db pubmlst --api-key)"
            )
            self.print_info(
                "Option 2: delete mlst_module/db/pubmlst and run "
                "`acinetoscope --pull-mlst-db` for a fresh prebuilt copy"
            )
            return False

        self.print_info(f"Using MLST database directory: {db_dir}")
        try:
            return _update(db_dir, script_dir, scheme="abaumannii")
        except TypeError:
            # Module's signature may be (db_dir, script_dir) without scheme
            return _update(db_dir, script_dir)

    # ------------------------------------------------------------------
    # Gene-centric reporter
    # ------------------------------------------------------------------
    def run_gene_centric_reporter(self, output_dir: Path) -> bool:
        module_name = "gene_centric_module"
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.print_error(f"Gene-centric module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_gene_centric_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for gene-centric: {temp_dir}")

        try:
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))

            # Copy required HTMLs from module output dirs to reporter work dir
            html_files = {
                "FASTA_QC_summary.html": DIR_FASTA_QC,
                "pasteur_mlst_summary.html": DIR_MLST_PASTEUR,
                "oxford_mlst_summary.html": DIR_MLST_OXFORD,
                "Kaptive_summary.html": DIR_KAPTIVE,
                "acineto_amrfinder_summary_report.html": DIR_AMRFINDER,
                "acineto_card_summary_report.html": DIR_ABRICATE,
                "acineto_ncbi_summary_report.html": DIR_ABRICATE,
                "acineto_resfinder_summary_report.html": DIR_ABRICATE,
                "acineto_vfdb_summary_report.html": DIR_ABRICATE,
                "acineto_argannot_summary_report.html": DIR_ABRICATE,
                "acineto_megares_summary_report.html": DIR_ABRICATE,
                "acineto_ecoli_vf_summary_report.html": DIR_ABRICATE,
                "acineto_bacmet2_summary_report.html": DIR_ABRICATE,
                "acineto_ecoh_summary_report.html": DIR_ABRICATE,
                "acineto_plasmidfinder_summary_report.html": DIR_ABRICATE,
                "plasmid_summary_report.html": DIR_PLASMID,
                "mutation_summary.html": DIR_AMRFINDER,
                "acineto_mge_summary.html": DIR_MGE,
            }

            for fname, subdir in html_files.items():
                src = output_dir / subdir / fname
                if src.exists():
                    shutil.copy2(src, temp_dir / fname)
                    self.logger.info(f"Copied {fname} to gene-centric temp dir")
                else:
                    self.logger.warning(f"Required HTML not found: {src}")

            self.print_info("Running gene-centric reporter...")
            script = temp_dir / "genius_acinetobacter_gene_centric_reporter.py"
            if not script.exists():
                self.print_error("genius_acinetobacter_gene_centric_reporter.py not found")
                return False

            cmd = [sys.executable, str(script), "-i", "."]
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.warning(result.stderr)

            if result.returncode != 0:
                self.print_error(f"Gene-centric reporter failed (rc={result.returncode})")
                return False

            src_dir = temp_dir / DIR_GENE_CENTRIC
            if src_dir.exists():
                dst_dir = self.user_output_dir / DIR_GENE_CENTRIC
                if dst_dir.exists():
                    shutil.rmtree(dst_dir)
                shutil.copytree(src_dir, dst_dir)
                self.print_success("Gene-centric reporter completed successfully.")
                return True
            else:
                self.print_error(f"{DIR_GENE_CENTRIC} not produced by reporter")
                return False

        except Exception as e:
            self.print_error(f"Gene-centric reporter exception: {e}")
            self.logger.error(traceback.format_exc())
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Sample-centric reporter
    # ------------------------------------------------------------------
    def run_sample_centric_reporter(self, output_dir: Path) -> bool:
        module_name = "sample_centric_module"
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            self.print_error(f"Sample-centric module not found: {module_orig}")
            return False

        temp_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_sample_centric_"))
        self._register_temp_dir(temp_dir)
        self.logger.info(f"Temporary directory for sample-centric: {temp_dir}")

        try:
            shutil.copytree(module_orig, temp_dir, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))

            sample_files = {
                "acineto_amrfinder_summary.tsv": DIR_AMRFINDER,
                "acineto_card_abricate_summary.tsv": DIR_ABRICATE,
                "acineto_ncbi_abricate_summary.tsv": DIR_ABRICATE,
                "acineto_resfinder_abricate_summary.tsv": DIR_ABRICATE,
                "acineto_vfdb_abricate_summary.tsv": DIR_ABRICATE,
                "acineto_argannot_abricate_summary.tsv": DIR_ABRICATE,
                "acineto_megares_abricate_summary.tsv": DIR_ABRICATE,
                "acineto_ecoli_vf_abricate_summary.tsv": DIR_ABRICATE,
                "acineto_bacmet2_abricate_summary.tsv": DIR_ABRICATE,
                "acineto_ecoh_abricate_summary.tsv": DIR_ABRICATE,
                "acineto_plasmidfinder_abricate_summary.tsv": DIR_ABRICATE,
                "plasmid_summary.tsv": DIR_PLASMID,
                "mutation_summary.tsv": DIR_AMRFINDER,
                "acineto_mge_summary.tsv": DIR_MGE,
                "FASTA_QC_summary.html": DIR_FASTA_QC,
                "pasteur_mlst_summary.html": DIR_MLST_PASTEUR,
                "oxford_mlst_summary.html": DIR_MLST_OXFORD,
                "Kaptive_summary.html": DIR_KAPTIVE,
            }

            for fname, subdir in sample_files.items():
                src = output_dir / subdir / fname
                if src.exists():
                    shutil.copy2(src, temp_dir / fname)
                    self.logger.info(f"Copied {fname} to sample-centric temp dir")
                else:
                    self.logger.warning(f"Required file not found: {src}")

            self.print_info("Running sample-centric reporter...")
            script = temp_dir / "genius_acinetobacter_sample_centric_reporter.py"
            if not script.exists():
                self.print_error("genius_acinetobacter_sample_centric_reporter.py not found")
                return False

            cmd = [sys.executable, str(script), "-i", "."]
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True,
                                    text=True, input="")

            if result.stdout:
                self.logger.info(f"Sample-centric stdout:\n{result.stdout}")
            if result.stderr:
                self.logger.warning(f"Sample-centric stderr:\n{result.stderr}")

            if result.returncode != 0:
                self.print_error(f"Sample-centric reporter failed (rc={result.returncode})")
                return False

            src_dir = temp_dir / DIR_SAMPLE_CENTRIC
            if not src_dir.exists():
                # Fallback: any new directory with HTMLs
                for alt in ["SAMPLE_CENTRIC_REPORTS", "GENIUS_ACINETOBACTER_REPORTS",
                            "output", "reports"]:
                    alt_dir = temp_dir / alt
                    if alt_dir.exists() and any(alt_dir.glob("*.html")):
                        src_dir = alt_dir
                        self.print_info(f"Reporter used alternative output dir: {alt_dir.name}")
                        break

            if src_dir.exists():
                dst_dir = self.user_output_dir / DIR_SAMPLE_CENTRIC
                if dst_dir.exists():
                    shutil.rmtree(dst_dir)
                shutil.copytree(src_dir, dst_dir)
                self.print_success("Sample-centric reporter completed successfully.")
                return True
            else:
                self.print_error(f"{DIR_SAMPLE_CENTRIC} not produced by reporter")
                return False

        except Exception as e:
            self.print_error(f"Sample-centric reporter exception: {e}")
            self.logger.error(traceback.format_exc())
            return False
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Final report + de-duplication
    # ------------------------------------------------------------------
    def copy_summary_results_to_final_directory(self, output_dir: Path) -> None:
        """Copy deliverable folders into AcinetoScope_final_report/."""
        try:
            self.print_info("Copying summary results to final directory...")
            final_dir = output_dir / DIR_FINAL
            final_dir.mkdir(parents=True, exist_ok=True)

            for src_name in [DIR_GENE_CENTRIC, DIR_SAMPLE_CENTRIC]:
                src = output_dir / src_name
                if src.exists() and src.is_dir():
                    dst = final_dir / src_name
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.print_info(f"  * Copied: {src_name}/")

            self.print_success(f"All results copied to: {final_dir}")
        except Exception as e:
            self.print_error(f"Error copying summary results: {e}")

    def _cleanup_top_level_duplicates(self, output_path: Path) -> None:
        """Remove top-level dirs copied into AcinetoScope_final_report/."""
        removed = 0
        for dirname in [DIR_GENE_CENTRIC, DIR_SAMPLE_CENTRIC]:
            d = output_path / dirname
            if d.exists():
                shutil.rmtree(d, ignore_errors=True)
                removed += 1
        if removed:
            self.logger.info(
                f"Cleaned up {removed} duplicate dir(s) from top level"
            )

    def _finalize_output_by_report_level(self, output_path: Path) -> None:
        """Prune intermediate files based on report level.

        final   : Only AcinetoScope_final_report/ is kept.
        summary : AcinetoScope_final_report/ + batch files in module dirs.
        full    : No pruning.
        """
        if self.report_level == 'full':
            self.print_info("Report level: full (all intermediate files kept)")
            return

        final_dir = output_path / DIR_FINAL
        final_exists = final_dir.exists() and any(final_dir.iterdir())
        if not final_exists:
            self.print_warning(
                "Final report directory is missing or empty - skipping prune "
                "to avoid losing all output."
            )
            return

        if self.report_level == 'final':
            removed = 0
            for module_dir_name in ALL_MODULE_DIRS:
                module_dir = output_path / module_dir_name
                if not module_dir.exists():
                    continue
                shutil.rmtree(module_dir, ignore_errors=True)
                removed += 1
            self.print_success(
                f"Report level: final - removed {removed} module directories"
            )
            return

        if self.report_level == 'summary':
            removed = 0
            for module_dir_name in ALL_MODULE_DIRS:
                module_dir = output_path / module_dir_name
                if not module_dir.exists():
                    continue
                for item in list(module_dir.iterdir()):
                    if item.is_dir():
                        shutil.rmtree(item, ignore_errors=True)
                        removed += 1
            if removed:
                self.print_success(
                    f"Report level: summary - removed {removed} per-sample folders"
                )
            else:
                self.print_success(
                    "Report level: summary - only batch summaries retained"
                )

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
                              clean_output: bool = False,
                              report_level: str = 'summary'):
        if skip_modules is None:
            skip_modules = {}
        self.report_level = report_level

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

            # Pre-create output dirs
            for subdir in [DIR_FASTA_QC, DIR_MLST_PASTEUR, DIR_MLST_OXFORD,
                           DIR_KAPTIVE, DIR_AMRFINDER, DIR_ABRICATE,
                           DIR_PLASMID, DIR_MGE]:
                (output_path / subdir).mkdir(exist_ok=True)

            # Analysis plan
            self.print_header("Analysis Plan", "Modules to be executed")

            # Warn early if MGE is enabled but the DB is missing
            if not skip_modules.get('mge', False) and self.get_mge_data_dir() is None:
                self.print_warning(
                    "MGE profiling is enabled but mobileOG-db is not installed. "
                    "Run `acinetoscope --pull-mge-db` first, or add --skip-mge."
                )

            analyses = [
                ("QC", not skip_modules.get('qc', False)),
                ("MLST Pasteur", not skip_modules.get('mlst', False) and mlst_scheme in ["pasteur", "both"]),
                ("MLST Oxford", not skip_modules.get('mlst', False) and mlst_scheme in ["oxford", "both"]),
                ("Kaptive", not skip_modules.get('kaptive', False)),
                ("AMRfinderPlus", not skip_modules.get('amr', False)),
                ("ABRicate", not skip_modules.get('abricate', False)),
                ("Plasmid Typing (APT)", not skip_modules.get('plasmid', False)),
                ("MGE Profiling", not skip_modules.get('mge', False)),
                ("Gene-centric Reporter", not skip_modules.get('gene_centric', False)),
                ("Sample-centric Reporter", not skip_sample_centric and not skip_modules.get('sample_centric', False)),
            ]
            for name, enabled in analyses:
                status = "✅ ENABLED" if enabled else "⏸️  SKIPPED"
                print(f"   {status} - {name}")
            print(f"\n   Report level: {self.report_level}")
            sys.stdout.flush()

            # Ensure AMR database
            if not skip_modules.get('amr', False):
                if not self.ensure_amr_database():
                    self.print_warning("AMR database check failed; AMR analysis may fail.")

            # ── Step 1: QC runs first — its species gate filters downstream ──
            qc_success = False
            if not skip_modules.get('qc', False):
                self.print_header("QC Analysis", "Quality Control and Statistics")
                qc_success = self.run_qc_analysis(fasta_files, output_path, threads)
                if qc_success:
                    self.print_success("QC Analysis completed successfully.")
                else:
                    self.print_warning("QC Analysis had issues (check log).")
                self.display_random_quote()

            # ── Step 2: Species gate (fail-open if QC failed or gate missing) ──
            downstream_files = fasta_files
            if qc_success:
                gate = self._load_species_gate(output_path)
                downstream_files = self._filter_by_species_gate(fasta_files, gate)
                if not downstream_files:
                    self.print_error(
                        "No samples passed the species gate - "
                        "downstream modules will be skipped."
                    )

            # ── Step 3: Downstream modules (run on gated file set) ──
            if downstream_files:
                module_configs = [
                    ("MLST Analysis - PASTEUR", "Multi-Locus Sequence Typing",
                     not skip_modules.get('mlst', False) and mlst_scheme in ["pasteur", "both"],
                     lambda: self.run_mlst_analysis(downstream_files, output_path, threads, "pasteur")),
                    ("MLST Analysis - OXFORD", "Multi-Locus Sequence Typing",
                     not skip_modules.get('mlst', False) and mlst_scheme in ["oxford", "both"],
                     lambda: self.run_mlst_analysis(downstream_files, output_path, threads, "oxford")),
                    ("Kaptive Analysis", "K and O Locus Typing",
                     not skip_modules.get('kaptive', False),
                     lambda: self.run_kaptive_analysis(downstream_files, output_path, threads)),
                    ("AMR Analysis", "Antimicrobial Resistance Gene Detection",
                     not skip_modules.get('amr', False),
                     lambda: self.run_amr_analysis(downstream_files, output_path, threads,
                                                    amr_min_identity, amr_min_coverage,
                                                    amr_skip_mutations, amr_force_update)),
                    ("ABRicate Analysis", "Comprehensive Resistance & Virulence Screening",
                     not skip_modules.get('abricate', False),
                     lambda: self.run_abricate_analysis(downstream_files, output_path, threads,
                                                         abricate_min_identity, abricate_min_coverage)),
                    ("Plasmid Typing (APT)", "Acinetobacter Plasmid Typing",
                     not skip_modules.get('plasmid', False),
                     lambda: self.run_plasmid_analysis(downstream_files, output_path, threads)),
                    ("MGE Profiling", "Mobile Genetic Element Detection (mobileOG-db)",
                     not skip_modules.get('mge', False),
                     lambda: self.run_mge_analysis(downstream_files, output_path, threads)),
                ]

                for title, subtitle, enabled, func in module_configs:
                    if enabled:
                        self.print_header(title, subtitle)
                        success = func()
                        if success:
                            self.print_success(f"{title} completed successfully.")
                        else:
                            self.print_warning(f"{title} had issues (check log).")
                        self.display_random_quote()

            # Reporters
            self.print_header("Generating Reports", "Ultimate Reporters")
            if not skip_modules.get('gene_centric', False):
                self.print_header("Gene-centric Reporter", "Ultimate gene-centric integrated report")
                gene_success = self.run_gene_centric_reporter(output_path)
                if gene_success:
                    self.print_success("Gene-centric reporter completed successfully.")
                else:
                    self.print_warning("Gene-centric reporter had issues (check log).")
                self.display_random_quote()

            if not skip_sample_centric and not skip_modules.get('sample_centric', False):
                self.print_header("Sample-centric Reporter", "Interactive isolate-centric report")
                sample_success = self.run_sample_centric_reporter(output_path)
                if sample_success:
                    self.print_success("Sample-centric reporter completed successfully.")
                else:
                    self.print_warning("Sample-centric reporter had issues (check log).")
                self.display_random_quote()

            # Consolidate into final report folder + clean top-level duplicates
            self.copy_summary_results_to_final_directory(output_path)
            self._cleanup_top_level_duplicates(output_path)

            # Prune per report level
            self._finalize_output_by_report_level(output_path)

            # Footer
            analysis_time = datetime.now() - start_time
            analysis_time_str = str(analysis_time).split('.')[0]

            self.print_header("Analysis Complete", f"Time elapsed: {analysis_time_str}")
            self.print_success(f"🎉 Output ready in: {output_path}")
            self.print_info(f"Final reports in: {output_path / DIR_FINAL}")
            self.display_random_quote()

            # Print final structure
            print("\nFinal output structure:")
            for item in sorted(output_path.glob("*")):
                if item.is_dir():
                    count = len(list(item.rglob("*")))
                    print(f"   {item.name}/ ({count} items)")
                else:
                    print(f"   {item.name}")

            # Print citations
            self.print_color(
                "\n📚 Please cite AcinetoScope:",
                Color.BOLD + Color.BRIGHT_CYAN,
            )
            for line in CITATION_LINES:
                self.print_color(f"   {line}", Color.CYAN)
            print()

        except KeyboardInterrupt:
            self.print_error("Analysis interrupted by user")
        except Exception as e:
            self.print_error(f"Critical error: {str(e)}")
            traceback.print_exc()


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

  # Report level (disk footprint control)
  acinetoscope -i "*.fna" -o results --report-level final     # minimal disk usage
  acinetoscope -i "*.fna" -o results --report-level summary   # default
  acinetoscope -i "*.fna" -o results --report-level full      # everything kept

  # Update AMR database only
  acinetoscope --update-amr-db
  acinetoscope --force-update-amr-db

  # Update MLST database
  acinetoscope --pull-mlst-db                # clone from GitHub (no API key)
  acinetoscope --update-mlst-db              # refresh from PubMLST (API key required)

  # Pull the MGE profiling database (mobileOG-db, ~300 MB)
  acinetoscope --pull-mge-db

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}MLST Database Management:{ColoredHelpFormatter.RESET}
  --pull-mlst-db     : Download a prebuilt A. baumannii MLST database from GitHub.
                       Prunes to Oxford + Pasteur schemes only. No API key needed.
  --update-mlst-db   : Update the MLST database directly from PubMLST.
                       Requires one-time API key setup (see below).
                       Download may take 10-20 minutes depending on internet speed.

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}MGE Database Management:{ColoredHelpFormatter.RESET}
  --pull-mge-db      : Download the mobileOG-db payload (~300 MB) from GitHub.
                       Required once after conda install if you want MGE profiling.
                       Docker and Singularity users do NOT need this —
                       the data is already baked into the image.

                       Location: ~/.local/share/acinetoscope/mge_module/mobileOG-db/
                       Env override: ACINETOSCOPE_MGE_DATA_DIR
                       Override URL: ACINETOSCOPE_MGE_PULL_URL

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}PubMLST API Key Setup (for --update-mlst-db):{ColoredHelpFormatter.RESET}
  1. Register/login at https://pubmlst.org/site-accounts
  2. Generate an API key from your profile (under "API keys").
  3. Run: mlstdb connect --db pubmlst --api-key
  4. Paste your API key when prompted.
  5. After setup, run acinetoscope --update-mlst-db whenever you want to refresh the database.

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}HPC / Docker Environment Variables:{ColoredHelpFormatter.RESET}
  All overrides are optional. Useful when the module directory is read-only
  (Docker images, system-wide conda, shared HPC installs) or when databases
  live on a shared filesystem (Lustre, NFS, BeeGFS).

  ACINETOSCOPE_SHARED_DATA         Writable root for user-local data
                                    (default: ~/.local/share/acinetoscope/)

  ACINETOSCOPE_MLST_DB             MLST database directory
  ACINETOSCOPE_MLST_MODULE_PATH    MLST module root (bin/, db/, perl5/, scripts/)
  ACINETOSCOPE_MLST_PULL_URL       Override default GitHub pull URL

  ACINETOSCOPE_MGE_DATA_DIR        mobileOG-db data directory
  ACINETOSCOPE_MGE_MODULE_PATH     MGE module root
  ACINETOSCOPE_MGE_PULL_URL        Override default GitHub pull URL

  ACINETOSCOPE_AMRFINDER_BIN_DIR   Directory containing amrfinder + amrfinder_update
  ACINETOSCOPE_AMRFINDER_DB        AMRFinderPlus database versioned directory

  Example (HPC with shared DBs on Lustre):
    export ACINETOSCOPE_MLST_DB=/lustre/scratch/dbs/mlst/pubmlst
    export ACINETOSCOPE_MGE_DATA_DIR=/lustre/scratch/dbs/mge/mobileOG-db
    export ACINETOSCOPE_AMRFINDER_DB=/lustre/scratch/dbs/amrfinder/2026-08-07.1
    acinetoscope -i "*.fna" -o results --report-level final

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Report Levels:{ColoredHelpFormatter.RESET}
  final    : Only AcinetoScope_final_report/ is kept. All module dirs are deleted.
             Best for archival and very large batches.
  summary  : (default) AcinetoScope_final_report/ + batch TSVs/HTMLs/JSONs in each
             module directory. Per-sample folders are removed.
  full     : Everything is kept, including all per-sample HTMLs.
             Best for small batches or when debugging.

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Supported FASTA formats:{ColoredHelpFormatter.RESET} .fna, .fasta, .fa, .fn

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Analysis Modules:{ColoredHelpFormatter.RESET}
  • FASTA QC (Quality Control & Species Confirmation via fastANI)
  • MLST (Pasteur & Oxford schemes)
  • Kaptive (K and O locus typing)
  • AMR profiling (AMRFinderPlus) – point mutations reported by default
  • ABRicate (Comprehensive resistance/plasmid/virulence)
  • Plasmid typing (APT – Acinetobacter Plasmid Typing)
  • MGE profiling (mobileOG-db) – Mobile Genetic Element detection
  • Gene-centric ultimate reporter
  • Sample-centric ultimate reporter (isolate-centric)

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Species Gate:{ColoredHelpFormatter.RESET}
  After FASTA QC completes, only samples that pass fastANI species confirmation
  (A. baumannii, ANI ≥ 95%) proceed through downstream modules. The gate is
  read from fasta_qc_results/FASTA_QC_species_passed.txt. If the gate file is
  missing (QC skipped or failed), ALL samples run through downstream modules.

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Output:{ColoredHelpFormatter.RESET}
  All intermediate results are preserved in the output directory.
  Final reports are placed in AcinetoScope_final_report/ containing:
    • GENIUS_ACINETOBACTER_GENE_CENTRIC_REPORTS/
    • GENIUS_ACINETOBACTER_SAMPLE_CENTRIC_REPORTS/
  A detailed log file is written to the output directory root.

{ColoredHelpFormatter.BOLD}{ColoredHelpFormatter.GREEN}Please cite AcinetoScope:{ColoredHelpFormatter.RESET}
  Beckley, B., Amarh, V., Lopes, B.S. et al. AcinetoScope: a computational pipeline for rapid, comprehensive Acinetobacter baumannii outbreak investigation and resistance gene tracking. Sci Rep (2026). https://doi.org/10.1038/s41598-026-70189-y

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
    parser.add_argument('--skip-mge', action='store_true', help='Skip MGE profiling (mobileOG) analysis')
    parser.add_argument('--skip-gene-centric', action='store_true', help='Skip gene-centric ultimate reporter')
    parser.add_argument('--skip-sample-centric', action='store_true', help='Skip sample-centric ultimate reporter')

    # MLST scheme
    parser.add_argument('--mlst-scheme', choices=['pasteur', 'oxford', 'both'], default='both',
                        help='MLST scheme to use (default: both)')

    # MLST DB management
    parser.add_argument('--pull-mlst-db', action='store_true',
                        help='Pull A. baumannii MLST database from GitHub (recommended for first-time setup)')
    parser.add_argument('--update-mlst-db', action='store_true',
                        help='Update A. baumannii MLST database from PubMLST (API key required)')

    # MGE DB management
    parser.add_argument('--pull-mge-db', action='store_true',
                        help='Pull the mobileOG-db payload (MGE profiling data) '
                             'from GitHub (~300 MB, no API key)')

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

    # Report level
    parser.add_argument('--report-level', choices=['final', 'summary', 'full'],
                        default='summary',
                        help='Controls which files are kept after the run. '
                             'final = only AcinetoScope_final_report/ ; '
                             'summary = (default) final report + batch files, per-sample folders removed ; '
                             'full = everything including per-sample folders.')

    # Misc flags
    parser.add_argument('--keep-temp', action='store_true', help='Do not delete temporary directories (debugging)')
    parser.add_argument('--clean-output', action='store_true', help='Delete output directory before analysis')

    args = parser.parse_args()

    # DB management early exits
    if args.update_amr_db or args.force_update_amr_db:
        orch = AcinetoScopeOrchestrator()
        if args.force_update_amr_db:
            orch.update_amr_database(force=True)
        else:
            orch.update_amr_database(force=False)
        sys.exit(0)

    if args.pull_mlst_db or args.update_mlst_db:
        orch = AcinetoScopeOrchestrator()
        temp_log_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_mlst_db_"))
        orch.setup_logging(temp_log_dir)
        if args.pull_mlst_db:
            success = orch.pull_mlst_database()
        else:
            success = orch.update_mlst_database()
        if not orch.keep_temp:
            shutil.rmtree(temp_log_dir, ignore_errors=True)
        sys.exit(0 if success else 1)

    if args.pull_mge_db:
        orch = AcinetoScopeOrchestrator()
        temp_log_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_mge_db_"))
        orch.setup_logging(temp_log_dir)
        success = orch.pull_mge_database()
        if not orch.keep_temp:
            shutil.rmtree(temp_log_dir, ignore_errors=True)
        sys.exit(0 if success else 1)

    if not args.input or not args.output:
        parser.error("Both -i/--input and -o/--output are required "
                     "(unless using --update-amr-db / --pull-mlst-db / --update-mlst-db / --pull-mge-db).")

    # Guard: report_level=final requires at least one reporter
    if args.report_level == 'final' and args.skip_gene_centric and args.skip_sample_centric:
        parser.error(
            "--report-level final requires at least one reporter to run "
            "(--skip-gene-centric and --skip-sample-centric both set would leave "
            "the final report empty and pruning would delete all output)."
        )

    skip_modules = {
        'qc': args.skip_qc,
        'mlst': args.skip_mlst,
        'kaptive': args.skip_kaptive,
        'amr': args.skip_amr,
        'abricate': args.skip_abricate,
        'plasmid': args.skip_plasmid,
        'mge': args.skip_mge,
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
        clean_output=args.clean_output,
        report_level=args.report_level,
    )


if __name__ == "__main__":
    main()