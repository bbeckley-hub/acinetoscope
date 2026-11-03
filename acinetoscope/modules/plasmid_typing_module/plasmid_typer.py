#!/usr/bin/env python3
"""
AcinetoScope Plasmid Typing Module - AcinetobacterPlasmidTyping (APT)
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-07-17 (updated 2026-10-29) (updated 2026-11-01)
Version: 1.4.0

Uses APT database v3.0 (Feb 2025) for Acinetobacter plasmid rep typing.

Per-sample outputs (1 file each, inside ONE level of subfolder):
    <output>/<sample>/<sample>_plasmid_report.html

Batch-level outputs (at top of -o):
    plasmid_summary.tsv
    plasmid_summary.json    (single batch JSON)
    plasmid_summary_report.html

All intermediate BLAST output lives inside a TemporaryDirectory that
self-destructs on exit — no intermediates left on disk.
"""

import subprocess
import sys
import os
import glob
import logging
import tempfile
import random
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any
import argparse
import re
from datetime import datetime
import psutil
import json
from collections import defaultdict


class AcinetoPlasmidTyper:
    def __init__(self, cpus: int = None):
        self.logger = self._setup_logging()
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        self.db_path = os.path.join(self.module_dir, "data", "apt_reps_v3.fasta")
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)
        self.min_identity = 95.0
        self.min_subject_coverage = 70.0
        self.rep_pattern = re.compile(r'^>?([rR][1-9P]+[-_][Tt][0-9]+)', re.IGNORECASE)

        self.metadata = {
            "tool_name": "AcinetoScope Plasmid Typing",
            "version": "1.4.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "database": "APT v3.0 (Feb 2025)",
            "blast_version": self._get_blast_version(),
            "thresholds": f"identity >= {self.min_identity}%, subject coverage >= {self.min_subject_coverage}%",
        }

        self.ascii_art = r""" █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝""".strip("\n")

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
            "The universe is not required to be in harmony with human ambition. - Carl Sagan",
        ]

    def _setup_logging(self):
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s')
        return logging.getLogger(__name__)

    def _get_available_ram(self) -> float:
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except Exception:
            return 8.0

    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        if user_cpus is not None:
            return user_cpus
        try:
            total_physical = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            if total_physical <= 4:
                return total_physical
            elif total_physical <= 8:
                return total_physical - 1
            else:
                return max(4, total_physical - 2)
        except Exception:
            return 4

    def _get_blast_version(self) -> str:
        try:
            res = subprocess.run(['blastn', '-version'], capture_output=True, text=True)
            return res.stdout.splitlines()[0] if res.stdout else "unknown"
        except Exception:
            return "not found"

    # ------------------------------------------------------------------
    # Shared CSS + ASCII CSS + footer (ABRicate-style consistency)
    # ------------------------------------------------------------------
    def _get_shared_css(self) -> str:
        return """
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff; padding: 20px; min-height: 100vh;
        }
        .container { max-width: 1400px; margin: 0 auto; }
        .card {
            background: rgba(255, 255, 255, 0.95); color: #1f2937;
            padding: 25px; margin: 20px 0; border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }
        .stats {
            display: flex; justify-content: space-around;
            margin: 20px 0; flex-wrap: wrap;
        }
        .stat-card {
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white; padding: 20px; border-radius: 12px;
            text-align: center; margin: 10px; flex: 1; min-width: 200px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }
        table {
            width: 100%; border-collapse: collapse; margin: 20px 0;
            background: white; border-radius: 8px; overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #e0e0e0; vertical-align: top; }
        th {
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white; font-weight: 600;
        }
        tr:hover { background-color: #f8f9fa; }
        a { color: #2563eb; text-decoration: none; }
        a:hover { text-decoration: underline; }
        .rep-badge {
            display: inline-block; background: #e74c3c; color: white;
            padding: 5px 12px; border-radius: 20px; margin: 3px; font-size: 0.9em;
        }
        .quote-container {
            background: rgba(255, 255, 255, 0.1); color: white;
            padding: 20px; border-radius: 12px; margin: 20px 0;
            text-align: center; font-style: italic; border-left: 4px solid #fff;
        }
        .interactive-controls {
            background: #f8f9fa; padding: 15px; border-radius: 8px;
            margin: 15px 0; display: flex; flex-wrap: wrap;
            gap: 10px; align-items: center;
        }
        .search-box { flex: 1; min-width: 200px; }
        .search-box input {
            width: 100%; padding: 8px 12px; border: 1px solid #ddd;
            border-radius: 4px; font-size: 14px;
        }
        .export-buttons { display: flex; gap: 8px; flex-wrap: wrap; }
        .export-buttons button {
            padding: 8px 16px; background: #3b82f6; color: white;
            border: none; border-radius: 4px; cursor: pointer; transition: background 0.3s;
        }
        .export-buttons button:hover { background: #2563eb; }
        .export-buttons button.print { background: #10b981; }
        .export-buttons button.print:hover { background: #059669; }
        .result-counter { font-size: 0.9em; color: #666; font-style: italic; }
        .footer {
            background: rgba(0, 0, 0, 0.8); color: white; padding: 30px;
            border-radius: 12px; margin-top: 40px;
        }
        .footer a { color: #667eea; text-decoration: none; }
        .footer a:hover { text-decoration: underline; }
        .footer .timestamp { color: #fbbf24; font-weight: bold; margin-top: 10px; }
        @media print {
            .no-print, .interactive-controls, .quote-container, .ascii-container { display: none !important; }
            body { background: white !important; color: black !important; }
            .card { background: white !important; color: black !important; box-shadow: none !important; }
        }
        """

    def _get_ascii_css(self) -> str:
        return """
        .ascii-container {
            background: rgba(0, 0, 0, 0.9); padding: 20px; border-radius: 15px;
            margin-bottom: 20px; box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3); overflow-x: auto;
        }
        pre.ascii-art {
            font-family: 'Consolas', 'Courier New', 'DejaVu Sans Mono',
                         'Liberation Mono', 'Menlo', 'Monaco', monospace;
            font-size: 12px; line-height: 1.05; white-space: pre;
            color: #00ff00; text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            font-variant-ligatures: none; font-kerning: none;
            letter-spacing: 0; word-spacing: 0; margin: 0; padding: 0;
            display: inline-block; min-width: max-content; tab-size: 4;
        }
        @media (max-width: 768px) {
            pre.ascii-art { font-size: 6px; }
        }
        """

    def _get_footer_html(self, title: str) -> str:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return f"""
        <div class="footer">
            <p><strong>ACINETOSCOPE</strong> &mdash; Plasmid Typing Analysis Module</p>
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
    # Environment / database checks
    # ------------------------------------------------------------------
    def check_blast_installed(self) -> bool:
        try:
            subprocess.run(['blastn', '-version'], capture_output=True, check=True)
            return True
        except Exception:
            self.logger.error("blastn not found. Install BLAST+ (conda install -c bioconda blast)")
            return False

    def check_database(self) -> bool:
        if not os.path.exists(self.db_path):
            self.logger.error(f"APT database missing: {self.db_path}")
            self.logger.info("Download from https://github.com/MehradHamidian/AcinetobacterPlasmidTyping and place as data/apt_reps_v3.fasta")
            return False
        self.logger.info(f"Using APT database: {self.db_path}")
        return True

    def parse_rep_type(self, header: str) -> str:
        match = self.rep_pattern.search(header)
        if match:
            return match.group(1).upper().replace('_', '-')
        return "unknown"

    # ------------------------------------------------------------------
    # Per-genome BLAST + parse (all intermediates in temp dir)
    # ------------------------------------------------------------------
    def run_blast_on_genome(self, genome_file: str, out_dir: str) -> Dict:
        genome = Path(genome_file).stem
        os.makedirs(out_dir, exist_ok=True)

        try:
            with tempfile.TemporaryDirectory(prefix=f"apt_{genome}_") as temp_dir:
                blast_out = os.path.join(temp_dir, f"{genome}_blast.tab")
                cmd = [
                    'blastn', '-query', genome_file, '-subject', self.db_path,
                    '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
                    '-out', blast_out,
                    '-perc_identity', str(self.min_identity),
                    '-evalue', '1e-5',
                    '-num_threads', '1',
                ]
                try:
                    subprocess.run(cmd, check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e:
                    self.logger.error(f"BLAST failed {genome}: {e.stderr}")
                    return {'genome': genome, 'hits': [], 'hit_count': 0, 'status': 'failed'}

                raw_hits = self._parse_blast_output(blast_out)
                filtered_hits = [h for h in raw_hits
                                 if h['subject_coverage'] >= self.min_subject_coverage]

                result = {
                    'genome': genome,
                    'hits': filtered_hits,
                    'hit_count': len(filtered_hits),
                    'status': 'success',
                }

                # ONLY per-sample output: HTML report
                self._create_per_sample_html(genome, result, out_dir)
                return result

        except Exception as e:
            self.logger.error(f"Error processing {genome}: {e}")
            return {'genome': genome, 'hits': [], 'hit_count': 0, 'status': 'failed'}

    def _parse_blast_output(self, blast_file: str) -> List[Dict]:
        hits = []
        try:
            with open(blast_file) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 14:
                        continue
                    (qseqid, sseqid, pident, length, mismatch, gapopen,
                     qstart, qend, sstart, send, evalue, bitscore, qlen, slen) = parts
                    pident = float(pident)
                    align_len = int(length)
                    subj_len = int(slen)
                    subject_cov = (align_len / subj_len) * 100
                    rep_type = self.parse_rep_type(sseqid)
                    hits.append({
                        'query_sequence': qseqid,
                        'subject_header': sseqid,
                        'rep_type': rep_type,
                        'identity': pident,
                        'alignment_length': align_len,
                        'subject_coverage': round(subject_cov, 2),
                        'evalue': float(evalue),
                        'bitscore': float(bitscore),
                    })
        except Exception as e:
            self.logger.error(f"Parse error {blast_file}: {e}")
        return hits

    # ------------------------------------------------------------------
    # Per-sample HTML report (the ONLY per-sample output)
    # ------------------------------------------------------------------
    def _create_per_sample_html(self, genome: str, result: Dict, out_dir: str):
        hits = result['hits']
        sample_name = Path(genome).name
        unique_reps = sorted(set(h['rep_type'] for h in hits))
        random_quote = random.choice(self.science_quotes)

        # Build rows outside the f-string to keep it readable
        hit_rows = ""
        for h in hits:
            hit_rows += f"""
                <tr>
                    <td>{h['query_sequence']}</td>
                    <td><strong>{h['rep_type']}</strong></td>
                    <td>{h['identity']:.2f}</td>
                    <td>{h['subject_coverage']:.1f}</td>
                    <td>{h['evalue']:.1e}</td>
                    <td>{h['bitscore']:.1f}</td>
                </tr>
            """

        rep_badges = "".join(f'<span class="rep-badge">{r}</span> ' for r in unique_reps)

        if hits:
            hits_table_section = f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">
                All Plasmid Hits (Filtered by Coverage &ge; {self.min_subject_coverage}%)
            </h2>
            <div class="interactive-controls no-print">
                <div class="search-box">
                    <input type="text" id="search-hits" placeholder="Search rep type, contig..."
                           onkeyup="searchTable('hits-table', this.value)">
                </div>
                <div class="export-buttons">
                    <button onclick="exportToCSV('hits-table', '{sample_name}_plasmid_hits.csv')">Export CSV</button>
                    <button onclick="printReport()" class="print">Print</button>
                </div>
                <div class="result-counter" id="result-counter-hits-table">{len(hits)} results found</div>
            </div>
            <table id="hits-table">
                <thead>
                    <tr>
                        <th>Query Contig</th>
                        <th>Rep Type</th>
                        <th>Identity (%)</th>
                        <th>Subject Coverage (%)</th>
                        <th>E-value</th>
                        <th>Bitscore</th>
                    </tr>
                </thead>
                <tbody>{hit_rows}</tbody>
            </table>
        </div>
            """
        else:
            hits_table_section = """
        <div class="card">
            <h2>No Plasmid Rep Types Detected</h2>
            <p>No hits passed the identity and coverage thresholds.</p>
        </div>
            """

        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AcinetoScope Plasmid Typing Report - {sample_name}</title>
    <style>{self._get_shared_css()}{self._get_ascii_css()}</style>
    <script>
        function searchTable(tableId, term) {{
            const table = document.getElementById(tableId);
            const rows = table.getElementsByTagName('tr');
            let visible = 0;
            for (let i = 1; i < rows.length; i++) {{
                const text = rows[i].textContent.toLowerCase();
                if (text.includes(term.toLowerCase())) {{
                    rows[i].style.display = ''; visible++;
                }} else {{ rows[i].style.display = 'none'; }}
            }}
            const counter = document.getElementById('result-counter-' + tableId);
            if (counter) counter.innerText = visible + ' results found';
        }}
        function exportToCSV(tableId, filename) {{
            const table = document.getElementById(tableId);
            let csv = [];
            const headers = table.querySelectorAll('th');
            const headerRow = [];
            headers.forEach(h => headerRow.push(h.textContent));
            csv.push(headerRow.join(','));
            const rows = table.querySelectorAll('tbody tr');
            rows.forEach(row => {{
                if (row.style.display !== 'none') {{
                    const cells = row.querySelectorAll('td');
                    const rowData = [];
                    cells.forEach(c => rowData.push(c.textContent.trim().replace(/,/g,';')));
                    csv.push(rowData.join(','));
                }}
            }});
            const blob = new Blob([csv.join('\\n')], {{type: 'text/csv'}});
            const a = document.createElement('a');
            a.href = URL.createObjectURL(blob);
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(a.href);
        }}
        function printReport() {{ window.print(); }}
        let quotes = {json.dumps(self.science_quotes)};
        let currentQuote = 0;
        function rotateQuote() {{
            document.getElementById('science-quote').innerHTML = quotes[currentQuote];
            currentQuote = (currentQuote + 1) % quotes.length;
        }}
        document.addEventListener('DOMContentLoaded', function() {{
            rotateQuote();
            setInterval(rotateQuote, 10000);
        }});
    </script>
</head>
<body>
<div class="container">
    <div class="ascii-container">
        <pre class="ascii-art">{self.ascii_art}</pre>
    </div>

    <div class="quote-container">
        <div id="science-quote" style="font-size: 1.1em;">"{random_quote}"</div>
    </div>

    <div class="card">
        <h1 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">AcinetoScope Plasmid Typing Report</h1>
        <p><strong>Sample Name:</strong> {sample_name}</p>
        <p><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
        <p><strong>Database:</strong> {self.metadata['database']}</p>
        <p><strong>Thresholds:</strong> {self.metadata['thresholds']}</p>
        <div class="stats">
            <div class="stat-card">Total Hits: {len(hits)}</div>
            <div class="stat-card">Unique Rep Types: {len(unique_reps)}</div>
        </div>
        <h2>Detected Rep Types</h2>
        <div>{rep_badges if rep_badges else '<em>None detected</em>'}</div>
    </div>

    {hits_table_section}

    {self._get_footer_html("AcinetoScope - Plasmid Typing Per-Sample Report")}
</div>
</body>
</html>
"""

        html_file = os.path.join(out_dir, f"{sample_name}_plasmid_report.html")
        with open(html_file, 'w') as f:
            f.write(html)
        self.logger.info(f"HTML report: {html_file}")

    # ------------------------------------------------------------------
    # Orchestration
    # ------------------------------------------------------------------
    def process_single_genome(self, genome_file: str, out_base: str) -> Dict:
        """Creates the ONE per-sample subfolder and passes it to run_blast_on_genome.
        Only ONE place appends genome_name -> one level of nesting.
        """
        genome = Path(genome_file).stem
        genome_dir = os.path.join(out_base, genome)
        os.makedirs(genome_dir, exist_ok=True)
        return self.run_blast_on_genome(genome_file, genome_dir)

    def process_multiple_genomes(self, pattern: str,
                                 out_base: str = "acineto_plasmid_results") -> Dict:
        exts = [pattern, f"{pattern}.fasta", f"{pattern}.fa", f"{pattern}.fna"]
        files = []
        for e in exts:
            files.extend(glob.glob(e))
        files = list(set(files))
        if not files:
            raise FileNotFoundError(f"No FASTA files matching {pattern}")
        self.logger.info(f"Found {len(files)} genomes")

        os.makedirs(out_base, exist_ok=True)
        max_concurrent = max(1, min(self.cpus, len(files), int(self.available_ram)))
        self.logger.info(f"Using {max_concurrent} concurrent BLAST jobs")

        results = {}
        with ThreadPoolExecutor(max_workers=max_concurrent) as ex:
            futures = {ex.submit(self.process_single_genome, f, out_base): f for f in files}
            for fut in as_completed(futures):
                f = futures[fut]
                try:
                    res = fut.result()
                    results[res['genome']] = res
                    self.logger.info(f"{res['genome']}: {res['hit_count']} rep hits")
                except Exception as e:
                    self.logger.error(f"FAILED: {f} - {e}")
                    results[Path(f).stem] = {
                        'genome': Path(f).stem,
                        'hits': [], 'hit_count': 0, 'status': 'failed',
                    }

        self._create_summary_tsv(results, out_base)
        self._create_summary_json(results, out_base)
        self._create_summary_html(results, out_base)
        return results

    # ------------------------------------------------------------------
    # Batch summaries: TSV, single JSON, HTML
    # ------------------------------------------------------------------
    def _create_summary_tsv(self, results: Dict, out_base: str):
        tsv_file = os.path.join(out_base, "plasmid_summary.tsv")
        with open(tsv_file, 'w') as f:
            f.write("Genome\tHitCount\tRepTypes\n")
            for name, res in results.items():
                reps = ",".join(sorted(set(h['rep_type'] for h in res['hits'])))
                f.write(f"{name}\t{res['hit_count']}\t{reps}\n")
        self.logger.info(f"Summary TSV: {tsv_file}")

    def _create_summary_json(self, results: Dict, out_base: str):
        """Write the SINGLE batch JSON."""
        rep_freq = defaultdict(int)
        rep_genomes = defaultdict(list)
        total_hits = 0

        genome_summaries = {}
        for name, res in results.items():
            reps = sorted(set(h['rep_type'] for h in res['hits']))
            for r in reps:
                rep_freq[r] += 1
                rep_genomes[r].append(name)
            total_hits += res['hit_count']
            genome_summaries[name] = {
                'status': res.get('status', 'unknown'),
                'hit_count': res['hit_count'],
                'rep_types': reps,
                'hits': res['hits'],
            }

        master = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'database': self.metadata['database'],
                'blast_version': self.metadata['blast_version'],
                'thresholds': self.metadata['thresholds'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(results),
            },
            'statistics': {
                'total_genomes': len(results),
                'genomes_with_plasmids': sum(1 for r in results.values() if r['hit_count'] > 0),
                'genomes_without_plasmids': sum(1 for r in results.values() if r['hit_count'] == 0),
                'total_hits': total_hits,
                'unique_rep_types': len(rep_freq),
                'rep_type_distribution': dict(sorted(rep_freq.items(),
                                                     key=lambda x: x[1], reverse=True)),
                'rep_type_genomes': {k: sorted(v) for k, v in rep_genomes.items()},
            },
            'genome_summaries': genome_summaries,
        }

        json_file = os.path.join(out_base, "plasmid_summary.json")
        with open(json_file, 'w') as f:
            json.dump(master, f, indent=2, default=str)
        self.logger.info(f"Master JSON (single): {json_file}")

    def _create_summary_html(self, results: Dict, out_base: str):
        rep_freq = defaultdict(int)
        genome_reps = {}
        for name, res in results.items():
            reps = sorted(set(h['rep_type'] for h in res['hits']))
            genome_reps[name] = reps
            for r in reps:
                rep_freq[r] += 1

        random_quote = random.choice(self.science_quotes)
        genomes_with_plasmids = sum(1 for r in results.values() if r['hit_count'] > 0)

        # Frequency rows
        freq_rows = ""
        for rt, cnt in sorted(rep_freq.items(), key=lambda x: x[1], reverse=True):
            glist = [name for name, reps in genome_reps.items() if rt in reps]
            freq_rows += f"<tr><td><strong>{rt}</strong></td><td>{cnt}</td><td>{', '.join(glist)}</td></tr>"

        # Per-genome rows (genome name links to its HTML report)
        genome_rows = ""
        for name, reps in sorted(genome_reps.items()):
            link = f"{name}/{name}_plasmid_report.html"
            if reps:
                badges = ' '.join(f'<span class="rep-badge">{r}</span>' for r in reps)
            else:
                badges = '<em>No Rep found</em>'
            genome_rows += f'<tr><td><a href="{link}"><strong>{name}</strong></a></td><td>{badges}</td></tr>'

        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AcinetoScope Plasmid Typing Summary</title>
    <style>{self._get_shared_css()}{self._get_ascii_css()}
        /* ════════ AcinetoScope uniformity override (auto-added) ════════ */
        .ascii-container {{ text-align: center !important; }}
        .ascii-art {{ display: inline-block !important; text-align: left !important; }}
        .footer {{ text-align: center !important; margin-top: 30px !important; padding: 20px !important; background: rgba(0, 0, 0, 0.3) !important; border-radius: 10px !important; font-size: 14px !important; overflow-wrap: break-word !important; }}
        .footer a {{ color: #fbbf24 !important; text-decoration: none !important; }}
        .footer a:hover {{ text-decoration: underline !important; }}
        .footer .timestamp {{ color: #fbbf24 !important; font-weight: bold !important; }}
        .footer .authorship {{ margin-top: 15px !important; padding: 15px !important; background: rgba(255, 255, 255, 0.1) !important; border-radius: 8px !important; font-size: 12px !important; line-height: 1.7 !important; }}
    </style>
    <script>
        function searchTable(tid, term) {{
            const t = document.getElementById(tid);
            const rows = t.getElementsByTagName('tr');
            let vis = 0;
            for (let i = 1; i < rows.length; i++) {{
                const txt = rows[i].textContent.toLowerCase();
                if (txt.includes(term.toLowerCase())) {{
                    rows[i].style.display = ''; vis++;
                }} else {{ rows[i].style.display = 'none'; }}
            }}
            const cnt = document.getElementById('result-counter-' + tid);
            if (cnt) cnt.innerText = vis + ' results found';
        }}
        function exportToCSV(tid, fname) {{
            const t = document.getElementById(tid); let csv = [];
            const headers = t.querySelectorAll('th');
            const headerRow = []; headers.forEach(h => headerRow.push(h.textContent));
            csv.push(headerRow.join(','));
            const rows = t.querySelectorAll('tbody tr');
            rows.forEach(row => {{
                if (row.style.display !== 'none') {{
                    const cells = row.querySelectorAll('td');
                    const rowData = [];
                    cells.forEach(c => rowData.push(c.textContent.trim().replace(/,/g,';')));
                    csv.push(rowData.join(','));
                }}
            }});
            const blob = new Blob([csv.join('\\n')], {{type: 'text/csv'}});
            const a = document.createElement('a');
            a.href = URL.createObjectURL(blob); a.download = fname;
            document.body.appendChild(a); a.click(); document.body.removeChild(a);
            URL.revokeObjectURL(a.href);
        }}
        function printReport() {{ window.print(); }}
        let quotes = {json.dumps(self.science_quotes)};
        let currentQuote = 0;
        function rotateQuote() {{
            document.getElementById('science-quote').innerHTML = quotes[currentQuote];
            currentQuote = (currentQuote + 1) % quotes.length;
        }}
        document.addEventListener('DOMContentLoaded', function() {{
            rotateQuote();
            setInterval(rotateQuote, 10000);
        }});
    </script>
</head>
<body>
<div class="container">
    <div class="ascii-container">
        <pre class="ascii-art">{self.ascii_art}</pre>
    </div>

    <div class="quote-container">
        <div id="science-quote" style="font-size: 1.1em;">"{random_quote}"</div>
    </div>

    <div class="card">
        <h1 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">AcinetoScope Plasmid Typing Summary</h1>
        <p><strong>Database:</strong> {self.metadata['database']} | <strong>Thresholds:</strong> {self.metadata['thresholds']}</p>
        <p><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
        <div class="stats">
            <div class="stat-card">Genomes: {len(results)}</div>
            <div class="stat-card">With Plasmids: {genomes_with_plasmids}</div>
            <div class="stat-card">Unique Rep Types: {len(rep_freq)}</div>
        </div>
    </div>

    <div class="card">
        <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Rep Type Frequency</h2>
        <div class="interactive-controls no-print">
            <div class="search-box">
                <input type="text" id="freqSearch" placeholder="Search rep type..."
                       onkeyup="searchTable('freqTable', this.value)">
            </div>
            <div class="export-buttons">
                <button onclick="exportToCSV('freqTable', 'rep_frequency.csv')">Export CSV</button>
                <button onclick="printReport()" class="print">Print</button>
            </div>
            <div class="result-counter" id="result-counter-freqTable">{len(freq_rows.split('</tr>'))-1 if freq_rows else 0} results found</div>
        </div>
        <table id="freqTable">
            <thead><tr><th>Rep Type</th><th>Count</th><th>Genomes</th></tr></thead>
            <tbody>{freq_rows}</tbody>
        </table>
    </div>

    <div class="card">
        <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Per-Genome Rep Types</h2>
        <div class="interactive-controls no-print">
            <div class="search-box">
                <input type="text" id="genomeSearch" placeholder="Search genome..."
                       onkeyup="searchTable('genomeTable', this.value)">
            </div>
            <div class="export-buttons">
                <button onclick="exportToCSV('genomeTable', 'genome_reps.csv')">Export CSV</button>
            </div>
        </div>
        <table id="genomeTable">
            <thead><tr><th>Genome</th><th>Rep Types</th></tr></thead>
            <tbody>{genome_rows}</tbody>
        </table>
    </div>

    {self._get_footer_html("AcinetoScope - Plasmid Typing Batch Summary")}
</div>
</body>
</html>
"""
        html_file = os.path.join(out_base, "plasmid_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html)
        self.logger.info(f"Summary HTML: {html_file}")


def main():
    parser = argparse.ArgumentParser(
        description="AcinetoScope Plasmid Typing (APT)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
PER-SAMPLE OUTPUTS (1 file each, inside ONE level of subfolder):
    <output>/<sample>/<sample>_plasmid_report.html

BATCH-LEVEL OUTPUTS:
    plasmid_summary.tsv
    plasmid_summary.json    (single batch JSON)
    plasmid_summary_report.html
        """
    )
    parser.add_argument('pattern', help='FASTA file pattern (e.g., "*.fna")')
    parser.add_argument('--cpus', '-c', type=int, help='Number of CPU cores')
    parser.add_argument('--output', '-o', default='acineto_plasmid_results',
                        help='Output directory')
    args = parser.parse_args()

    typer = AcinetoPlasmidTyper(cpus=args.cpus)

    # Print banner ONCE
    print("\n" + "=" * 80)
    print(typer.ascii_art)
    print("=" * 80)
    print("AcinetoScope Plasmid Typing (APT)")
    print("=" * 80)
    print("Author: Brown Beckley | Email: brownbeckley94@gmail.com")
    print("Affiliation: University of Ghana Medical School - Department of Medical Biochemistry")
    print("=" * 80)

    if not typer.check_blast_installed():
        sys.exit(1)
    if not typer.check_database():
        sys.exit(1)

    results = typer.process_multiple_genomes(args.pattern, args.output)

    total_hits = sum(r['hit_count'] for r in results.values())
    with_plasmids = sum(1 for r in results.values() if r['hit_count'] > 0)

    print(f"\nAnalysis Complete. Results in {args.output}")
    print(f"  Genomes processed: {len(results)}")
    print(f"  Genomes with plasmids: {with_plasmids}")
    print(f"  Total rep hits: {total_hits}")

    print("\nPER-SAMPLE OUTPUTS (1 file each):")
    print("  <sample>/<sample>_plasmid_report.html")

    print("\nBATCH-LEVEL OUTPUTS:")
    print(f"  {args.output}/plasmid_summary.tsv")
    print(f"  {args.output}/plasmid_summary.json   (single batch JSON)")
    print(f"  {args.output}/plasmid_summary_report.html")

    print(f"\n{random.choice(typer.science_quotes)}")


if __name__ == "__main__":
    main()