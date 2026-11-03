#!/usr/bin/env python3
"""
AcinetoScope Kaptive K/O Locus Analysis - A. baumannii Capsule and Lipooligosaccharide Typing
Comprehensive Kaptive analysis for A. baumannii with beautiful HTML reporting
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-07-18 (updated 2026-08-23) final (updated 2026-11-01)
Version: 1.4.0

Per-sample outputs (1 file each, inside ONE level of subfolder):
    <output>/<sample>/<sample>_kaptive_report.html

Batch-level outputs (at top of -o):
    Kaptive_summary.tsv
    Kaptive_summary.json    (single batch JSON)
    Kaptive_summary.html

All intermediate files (raw K/O TSVs, combined TSV) live inside a
TemporaryDirectory that self-destructs on exit — no intermediates on disk.
"""

import subprocess
import sys
import os
import glob
import logging
import tempfile
import random
from pathlib import Path
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import json
from collections import defaultdict


class AcinetoKaptive:
    """Kaptive executor for A. baumannii K/O locus typing with comprehensive HTML reporting."""

    def __init__(self):
        self.logger = self._setup_logging()

        self.k_db_keyword = "ab_k"
        self.o_db_keyword = "ab_o"

        self.ascii_art = r""" █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝""".strip("\n")

        self.metadata = {
            "tool_name": "AcinetoScope Kaptive K/O Analysis",
            "version": "1.4.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "kaptive_version": self._get_kaptive_version(),
            "database_version": "2025-12-26 (ab_k/ab_o)",
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
            "The universe is not required to be in harmony with human ambition. - Carl Sagan",
        ]

    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)

    def _get_kaptive_version(self) -> str:
        try:
            result = subprocess.run(["kaptive", "--version"],
                                    capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            match = re.search(r'v?(\d+\.\d+\.\d+)', version_line)
            if match:
                return match.group(1)
            return version_line
        except Exception as e:
            self.logger.warning(f"Could not determine Kaptive version: {e}")
            return "Unknown"

    def check_kaptive_installed(self) -> bool:
        try:
            result = subprocess.run(["kaptive", "--version"],
                                    capture_output=True, text=True, check=True)
            self.logger.info(f"Kaptive version: {result.stdout.strip()}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"Kaptive check failed: {e}")
            self.logger.info("Please install Kaptive: conda install -c bioconda kaptive")
            return False

    def _get_footer_html(self, title: str) -> str:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return f"""
        <div class="footer">
            <p><strong>ACINETOSCOPE</strong> &mdash; K/O Typing Analysis Module</p>
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

    def _combine_k_o_results(self, k_tsv: str, o_tsv: str, combined_tsv: str):
        try:
            k_lines, o_lines = [], []
            if os.path.exists(k_tsv):
                with open(k_tsv, 'r') as f:
                    k_lines = f.readlines()
            if os.path.exists(o_tsv):
                with open(o_tsv, 'r') as f:
                    o_lines = f.readlines()
            with open(combined_tsv, 'w') as f:
                if k_lines:
                    f.write(k_lines[0])
                    for line in k_lines[1:]:
                        f.write(line)
                if o_lines:
                    for line in o_lines[1:]:
                        f.write(line)
            self.logger.info(f"Combined K and O results to {combined_tsv}")
        except Exception as e:
            self.logger.error(f"Error combining K/O results: {e}")

    def _parse_kaptive_output(self, kaptive_file: str) -> List[Dict]:
        hits = []
        try:
            if not os.path.exists(kaptive_file):
                return hits
            with open(kaptive_file, 'r') as f:
                lines = f.readlines()
            if not lines or len(lines) < 2:
                return hits
            headers = lines[0].strip().split('\t')
            num_headers = len(headers)
            for line_num, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < num_headers:
                    parts.extend([''] * (num_headers - len(parts)))
                elif len(parts) > num_headers:
                    parts = parts[:num_headers]
                hit = {h: parts[i] if i < len(parts) else ''
                       for i, h in enumerate(headers)}
                locus = hit.get('Best match locus', '')
                if 'KL' in locus:
                    hit['Locus Type'] = 'K'
                elif 'OCL' in locus:
                    hit['Locus Type'] = 'O'
                else:
                    hit['Locus Type'] = 'Unknown'
                gene_details = hit.get('Expected genes in locus, details', '')
                hit['Gene Details Parsed'] = self._parse_gene_details(gene_details)
                hits.append(hit)
        except Exception as e:
            self.logger.error(f"Error parsing {kaptive_file}: {e}")
        self.logger.info(f"Parsed {len(hits)} Kaptive hits from {kaptive_file}")
        return hits

    def _parse_gene_details(self, gene_details_str: str) -> List[Dict]:
        genes = []
        try:
            if not gene_details_str:
                return genes
            for gene_part in gene_details_str.split(';'):
                if not gene_part:
                    continue
                subparts = gene_part.split(',')
                if len(subparts) >= 3:
                    genes.append({
                        'name': subparts[0],
                        'identity': subparts[1].replace('%', ''),
                        'coverage': subparts[2].replace('%', ''),
                    })
        except Exception as e:
            self.logger.warning(f"Could not parse gene details: {e}")
        return genes

    # ------------------------------------------------------------------
    # ★ FIXED: output_dir is used DIRECTLY — NO genome_name appended.
    # ------------------------------------------------------------------
    def run_kaptive_single_genome(self, genome_file: str, output_dir: str) -> Dict[str, Any]:
        """Run Kaptive on one genome. Emits exactly 1 file:
             - <output_dir>/<genome>_kaptive_report.html
        output_dir is ALREADY the per-sample folder — do NOT append genome_name.
        All raw K/O/combined TSVs live inside a TemporaryDirectory.
        """
        genome_name = Path(genome_file).stem
        os.makedirs(output_dir, exist_ok=True)

        self.logger.info(f"Running Kaptive K locus analysis on {genome_name}")

        try:
            with tempfile.TemporaryDirectory(prefix=f"kaptive_{genome_name}_") as temp_dir:
                k_tsv = os.path.join(temp_dir, f"{genome_name}_K.tsv")
                o_tsv = os.path.join(temp_dir, f"{genome_name}_O.tsv")
                combined_tsv = os.path.join(temp_dir, f"{genome_name}_combined.tsv")

                subprocess.run(
                    ["kaptive", "assembly", self.k_db_keyword,
                     genome_file, "-o", k_tsv, "--verbose"],
                    capture_output=True, text=True, check=True
                )
                self.logger.info(f"Running Kaptive O locus analysis on {genome_name}")
                subprocess.run(
                    ["kaptive", "assembly", self.o_db_keyword,
                     genome_file, "-o", o_tsv, "--verbose"],
                    capture_output=True, text=True, check=True
                )

                self._combine_k_o_results(k_tsv, o_tsv, combined_tsv)
                hits = self._parse_kaptive_output(combined_tsv)

                # Write HTML DIRECTLY into output_dir (already the per-sample folder)
                self._create_kaptive_html_report(genome_name, hits, output_dir)

                return {
                    'genome': genome_name,
                    'hits': hits,
                    'hit_count': len(hits),
                    'status': 'success',
                }

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Kaptive failed for {genome_name}: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            return {
                'genome': genome_name,
                'hits': [], 'hit_count': 0,
                'status': 'failed',
                'error': str(e),
            }

    def _create_kaptive_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        basic_js = f"""
        <script>
            function sortTable(tableId, columnIndex) {{
                const table = document.getElementById(tableId);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                rows.sort((a, b) => {{
                    const aText = a.children[columnIndex].textContent.trim();
                    const bText = b.children[columnIndex].textContent.trim();
                    const aNum = parseFloat(aText.replace('%', ''));
                    const bNum = parseFloat(bText.replace('%', ''));
                    if (!isNaN(aNum) && !isNaN(bNum)) return aNum - bNum;
                    return aText.localeCompare(bText);
                }});
                rows.forEach(row => tbody.appendChild(row));
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
        """

        export_data_js = f"""
        <script>
            window.reportData = {{
                metadata: {{
                    genome: '{genome_name}',
                    date: '{self.metadata['analysis_date']}',
                    tool: '{self.metadata['tool_name']}',
                    version: '{self.metadata['version']}'
                }},
                hits: {json.dumps(hits, indent=2)}
            }};
            function exportToJSON() {{
                const dataStr = JSON.stringify(window.reportData, null, 2);
                const dataBlob = new Blob([dataStr], {{ type: 'application/json' }});
                const url = URL.createObjectURL(dataBlob);
                const a = document.createElement('a');
                a.href = url;
                a.download = '{genome_name}_kaptive_report.json';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                URL.revokeObjectURL(url);
            }}
        </script>
        """

        k_count = sum(1 for hit in hits if hit.get('Locus Type') == 'K')
        o_count = sum(1 for hit in hits if hit.get('Locus Type') == 'O')

        results_section = ""
        if hits:
            headers = list(hits[0].keys())
            headers_to_show = [h for h in headers if h not in ['Gene Details Parsed']]

            header_cells = "".join(
                f'<th class="sortable" onclick="sortTable(\'kaptive-results-table\', {i})">{h}</th>'
                for i, h in enumerate(headers_to_show)
            )

            body_rows = ""
            for hit in hits:
                locus_type = hit.get('Locus Type', '')
                row_class = ("k-locus-row" if locus_type == 'K'
                             else "o-locus-row" if locus_type == 'O' else "")
                cells = ""
                for header in headers_to_show:
                    value = hit.get(header, '')
                    if len(str(value)) > 100:
                        value = str(value)[:100] + '...'
                    cells += f'<td>{value}</td>'
                body_rows += f'<tr class="{row_class}">{cells}</tr>'

            gene_details_html = ""
            for hit in hits:
                locus = hit.get('Best match locus', '')
                locus_type = hit.get('Locus Type', '')
                gene_details = hit.get('Gene Details Parsed', [])
                if not gene_details:
                    continue
                color = '#28a745' if locus_type == 'K' else '#17a2b8'
                gene_items = "".join(
                    f'<span class="gene-item">{g["name"]}: {g["identity"]}% identity, '
                    f'{g["coverage"]}% coverage</span>'
                    for g in gene_details
                )
                gene_details_html += f"""
            <div class="gene-details-container">
                <h4 style="margin-top: 0; color: {color};">{locus} ({locus_type} Locus) - Gene Details:</h4>
                <div style="margin-top: 10px;">{gene_items}</div>
            </div>
                """

            results_section = f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Complete Kaptive Analysis Results</h2>
            <div class="table-container">
                <table class="gene-table" id="kaptive-results-table">
                    <thead><tr>{header_cells}</tr></thead>
                    <tbody>{body_rows}</tbody>
                </table>
            </div>
            <h3 style="margin-top: 30px; color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Detailed Gene Information</h3>
            {gene_details_html}
        </div>
            """
        else:
            results_section = """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">No K/O Loci Detected</h2>
            <p>No K or O loci found in this A. baumannii genome.</p>
        </div>
            """

        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope Kaptive K/O Analysis Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff; padding: 20px; min-height: 100vh;
        }}
        .container {{ max-width: 1800px; margin: 0 auto; }}
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
            box-shadow: 0 2px 10px rgba(0,0,0,0.1); font-size: 12px;
        }}
        .gene-table th, .gene-table td {{
            padding: 12px 8px; text-align: left;
            border-bottom: 1px solid #e0e0e0;
            max-width: 300px; overflow-wrap: break-word;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white; font-weight: 600; cursor: pointer; position: relative;
        }}
        .gene-table th:hover {{
            background: linear-gradient(135deg, #2563eb 0%, #1e3a8a 100%);
        }}
        .gene-table th.sortable::after {{
            content: "\\2195"; position: absolute; right: 8px; opacity: 0.6;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .summary-stats {{
            display: flex; justify-content: space-around;
            margin: 20px 0; flex-wrap: wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white; padding: 20px; border-radius: 12px;
            text-align: center; box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px; flex: 1; min-width: 200px;
        }}
        .k-stat-card {{
            background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%);
            color: white; padding: 20px; border-radius: 12px;
            text-align: center; box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px; flex: 1; min-width: 200px;
        }}
        .o-stat-card {{
            background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%);
            color: white; padding: 20px; border-radius: 12px;
            text-align: center; box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px; flex: 1; min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1); color: white;
            padding: 20px; border-radius: 12px; margin: 20px 0;
            text-align: center; font-style: italic; border-left: 4px solid #fff;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8); color: white; padding: 30px;
            border-radius: 12px; margin-top: 40px;
        }}
        .footer a {{ color: #667eea; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .footer .timestamp {{ color: #fbbf24; font-weight: bold; margin-top: 10px; }}
        .k-locus-row {{ background-color: #d4edda; border-left: 4px solid #28a745; }}
        .o-locus-row {{ background-color: #d1ecf1; border-left: 4px solid #17a2b8; }}
        .controls {{
            background: #f8f9fa; padding: 15px; border-radius: 8px;
            margin: 15px 0; display: flex; flex-wrap: wrap;
            gap: 10px; align-items: center;
        }}
        .export-buttons {{ display: flex; gap: 8px; flex-wrap: wrap; }}
        .export-buttons button {{
            padding: 8px 16px; background: #3b82f6; color: white;
            border: none; border-radius: 4px; cursor: pointer;
            font-size: 14px; transition: background 0.3s;
        }}
        .export-buttons button:hover {{ background: #2563eb; }}
        .export-buttons button.print {{ background: #10b981; }}
        .export-buttons button.print:hover {{ background: #059669; }}
        .export-buttons button.json {{ background: #8b5cf6; }}
        .export-buttons button.json:hover {{ background: #7c3aed; }}
        .gene-details-container {{
            max-height: 300px; overflow-y: auto; padding: 10px;
            background: #f8f9fa; border-radius: 5px;
            border: 1px solid #e0e0e0; margin-top: 10px; font-size: 11px;
        }}
        .gene-item {{
            display: inline-block; margin: 2px; padding: 3px 8px;
            background: #e9ecef; border-radius: 3px;
            font-size: 0.85em; white-space: nowrap;
        }}
        .table-container {{ overflow-x: auto; margin: 20px 0; }}
        @media (max-width: 1200px) {{
            .gene-table {{ font-size: 11px; }}
            .gene-table th, .gene-table td {{ padding: 8px 6px; }}
        }}
        @media print {{
            .controls {{ display: none !important; }}
            body {{ background: white !important; color: black !important; }}
            .card {{ background: white !important; color: black !important; box-shadow: none !important; }}
            .gene-table {{ box-shadow: none !important; }}
            .gene-details-container {{ max-height: none !important; }}
        }}
    </style>
    {basic_js}
    {export_data_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <pre class="ascii-art">{self.ascii_art}</pre>
            </div>
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">AcinetoScope Kaptive K/O Analysis Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive A. baumannii Capsule (K) and Lipooligosaccharide (O) Locus Analysis</p>
                <div class="controls">
                    <div class="export-buttons">
                        <button onclick="exportToJSON()" class="json">Export JSON</button>
                        <button onclick="printReport()" class="print">Print Report</button>
                    </div>
                </div>
            </div>
        </div>

        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">A. baumannii K/O Locus Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Loci Detected</h3>
                    <p style="font-size: 2em; margin: 0;">{len(hits)}</p>
                </div>
                <div class="k-stat-card">
                    <h3>K Locus (Capsule)</h3>
                    <p style="font-size: 2em; margin: 0;">{k_count}</p>
                </div>
                <div class="o-stat-card">
                    <h3>O Locus (Lipooligosaccharide)</h3>
                    <p style="font-size: 2em; margin: 0;">{o_count}</p>
                </div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>Kaptive Version:</strong> {self.metadata['kaptive_version']}</p>
            <p><strong>Database Version:</strong> {self.metadata['database_version']}</p>
        </div>

        {results_section}

        {self._get_footer_html("AcinetoScope - Kaptive K/O Per-Sample Report")}
    </div>
</body>
</html>
"""

        html_file = os.path.join(output_dir, f"{genome_name}_kaptive_report.html")
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        self.logger.info(f"Kaptive HTML report: {html_file}")

    # ------------------------------------------------------------------
    # ★ FIXED: creates the ONE per-sample folder and passes it directly
    # ------------------------------------------------------------------
    def process_single_genome(self, genome_file: str,
                              output_base: str = "kaptive_results") -> Dict[str, Any]:
        """Process one genome. Creates ONE per-sample subfolder and puts the
        HTML report inside it:
            <output_base>/<genome>/<genome>_kaptive_report.html
        """
        genome_name = Path(genome_file).stem
        self.logger.info(f"=== PROCESSING A. BAUMANNII GENOME: {genome_name} ===")

        if not self.check_kaptive_installed():
            self.logger.error("Kaptive not available!")
            return {
                'genome': genome_name,
                'hits': [], 'hit_count': 0,
                'status': 'failed',
                'error': 'Kaptive not available',
            }

        # Create the per-sample folder here — ONLY place that appends genome_name
        results_dir = os.path.join(output_base, genome_name)
        os.makedirs(results_dir, exist_ok=True)

        # Pass results_dir as-is; run_kaptive_single_genome will NOT append genome_name
        result = self.run_kaptive_single_genome(genome_file, results_dir)
        status_icon = "OK" if result['status'] == 'success' else "FAIL"
        self.logger.info(f"{status_icon} {genome_name}: {result['hit_count']} K/O loci found")
        return result

    def process_multiple_genomes(self, genome_pattern: str,
                                 output_base: str = "kaptive_results") -> Dict[str, Any]:
        fasta_extensions = ['.fna', '.fasta', '.fa', '.faa']
        genome_files = []
        genome_files.extend(glob.glob(genome_pattern))
        if not genome_files:
            for ext in fasta_extensions:
                genome_files.extend(glob.glob(f"{genome_pattern}{ext}"))
        genome_files = list(set([f for f in genome_files if os.path.exists(f)]))

        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")

        self.logger.info(f"Found {len(genome_files)} A. baumannii genomes: "
                         f"{[Path(f).name for f in genome_files]}")

        os.makedirs(output_base, exist_ok=True)
        all_results = {}

        for genome_file in genome_files:
            genome_name = Path(genome_file).stem
            try:
                result = self.process_single_genome(genome_file, output_base)
                all_results[genome_name] = result
                self.logger.info(f"COMPLETED: {genome_name} ({result['hit_count']} K/O loci)")
            except Exception as e:
                self.logger.error(f"FAILED: {genome_name} - {e}")
                all_results[genome_name] = {
                    'genome': genome_name,
                    'hits': [], 'hit_count': 0,
                    'status': 'failed',
                    'error': str(e),
                }

        self.create_kaptive_summary(all_results, output_base)
        self.logger.info("=== A. BAUMANNII K/O ANALYSIS COMPLETE ===")
        self.logger.info(f"Processed {len(all_results)} genomes")
        self.logger.info(f"Results saved to: {output_base}")
        return all_results

    def create_kaptive_summary(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating A. baumannii Kaptive summary files and HTML reports...")
        summary_file = os.path.join(output_base, "Kaptive_summary.tsv")
        all_headers = set()
        successful_hits = []

        for genome_name, result in all_results.items():
            if result.get('status') == 'success' and result.get('hits'):
                for hit in result['hits']:
                    hit_with_genome = hit.copy()
                    hit_with_genome['Genome'] = genome_name
                    successful_hits.append(hit_with_genome)
                    all_headers.update(hit.keys())

        if successful_hits:
            headers = ['Genome'] + [h for h in all_headers
                                    if h not in ['Genome', 'Gene Details Parsed']]
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write('\t'.join(headers) + '\n')
                for hit in successful_hits:
                    row = []
                    for header in headers:
                        value = hit.get(header, '')
                        if isinstance(value, str):
                            value = value.replace('\t', ' ').replace('\n', ' ')
                        row.append(str(value))
                    f.write('\t'.join(row) + '\n')

        self.logger.info(f"Kaptive summary TSV: {summary_file}")
        self.create_json_summary(all_results, output_base)
        self._create_summary_html_report(all_results, output_base)

    def create_json_summary(self, all_results: Dict[str, Any], output_base: str):
        master_summary = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'kaptive_version': self.metadata['kaptive_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results),
            },
            'genome_summaries': {},
            'statistics': self._create_summary_statistics(all_results),
        }

        for genome_name, result in all_results.items():
            master_summary['genome_summaries'][genome_name] = {
                'status': result.get('status', 'unknown'),
                'hit_count': result.get('hit_count', 0),
                'k_loci': sum(1 for hit in result.get('hits', [])
                              if hit.get('Locus Type') == 'K'),
                'o_loci': sum(1 for hit in result.get('hits', [])
                              if hit.get('Locus Type') == 'O'),
                'k_types': list(set(hit.get('Best match type', '')
                                    for hit in result.get('hits', [])
                                    if hit.get('Locus Type') == 'K')),
                'o_types': list(set(hit.get('Best match type', '')
                                    for hit in result.get('hits', [])
                                    if hit.get('Locus Type') == 'O')),
                'error': result.get('error', '') if result.get('status') == 'failed' else '',
            }

        json_file = os.path.join(output_base, "Kaptive_summary.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(master_summary, f, indent=2, ensure_ascii=False)
        self.logger.info(f"JSON summary (single): {json_file}")

    def _create_summary_statistics(self, all_results: Dict[str, Any]) -> Dict[str, Any]:
        stats = {
            'total_genomes': len(all_results),
            'successful_genomes': sum(1 for r in all_results.values()
                                      if r.get('status') == 'success'),
            'failed_genomes': sum(1 for r in all_results.values()
                                  if r.get('status') == 'failed'),
            'total_k_loci': 0,
            'total_o_loci': 0,
            'k_type_distribution': defaultdict(int),
            'o_type_distribution': defaultdict(int),
            'genomes_with_k': 0,
            'genomes_with_o': 0,
            'genomes_with_both': 0,
            'common_k_types': [],
            'common_o_types': [],
        }

        for result in all_results.values():
            if result.get('status') != 'success':
                continue
            hits = result.get('hits', [])
            k_hits = [h for h in hits if h.get('Locus Type') == 'K']
            o_hits = [h for h in hits if h.get('Locus Type') == 'O']

            stats['total_k_loci'] += len(k_hits)
            stats['total_o_loci'] += len(o_hits)

            if k_hits:
                stats['genomes_with_k'] += 1
                for hit in k_hits:
                    stats['k_type_distribution'][hit.get('Best match type', 'Unknown')] += 1
            if o_hits:
                stats['genomes_with_o'] += 1
                for hit in o_hits:
                    stats['o_type_distribution'][hit.get('Best match type', 'Unknown')] += 1
            if k_hits and o_hits:
                stats['genomes_with_both'] += 1

        if stats['k_type_distribution']:
            common_k = sorted(stats['k_type_distribution'].items(),
                              key=lambda x: x[1], reverse=True)[:5]
            stats['common_k_types'] = [{'type': k, 'count': v} for k, v in common_k]
        if stats['o_type_distribution']:
            common_o = sorted(stats['o_type_distribution'].items(),
                              key=lambda x: x[1], reverse=True)[:5]
            stats['common_o_types'] = [{'type': o, 'count': v} for o, v in common_o]

        stats['k_type_distribution'] = dict(stats['k_type_distribution'])
        stats['o_type_distribution'] = dict(stats['o_type_distribution'])
        return stats

    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        stats = self._create_summary_statistics(all_results)

        summary_rows = []
        for genome_name, result in all_results.items():
            if result.get('status') == 'success':
                hits = result.get('hits', [])
                k_hits = [h for h in hits if h.get('Locus Type') == 'K']
                o_hits = [h for h in hits if h.get('Locus Type') == 'O']

                k_info = "None"
                if k_hits:
                    k_hit = k_hits[0]
                    k_info = (f"{k_hit.get('Best match type', 'Unknown')} "
                              f"({k_hit.get('Identity', 'N/A')}, "
                              f"{k_hit.get('Coverage', 'N/A')})")
                o_info = "None"
                if o_hits:
                    o_hit = o_hits[0]
                    o_info = (f"{o_hit.get('Best match type', 'Unknown')} "
                              f"({o_hit.get('Identity', 'N/A')}, "
                              f"{o_hit.get('Coverage', 'N/A')})")

                summary_rows.append({
                    'genome': genome_name,
                    'status': 'Success',
                    'k_locus': k_info,
                    'o_locus': o_info,
                    'total_loci': len(hits),
                    'k_details': k_info if k_hits else "No K locus",
                    'o_details': o_info if o_hits else "No O locus",
                })
            else:
                summary_rows.append({
                    'genome': genome_name,
                    'status': 'Failed',
                    'k_locus': 'N/A',
                    'o_locus': 'N/A',
                    'total_loci': 0,
                    'error': result.get('error', 'Unknown error'),
                })

        interactive_js = f"""
        <script>
            function searchTable(tableId, searchTerm) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let visibleCount = 0;
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const text = row.textContent.toLowerCase();
                    if (text.includes(searchTerm.toLowerCase())) {{
                        row.style.display = ''; visibleCount++;
                    }} else {{ row.style.display = 'none'; }}
                }}
                const resultCounter = document.getElementById('result-counter-' + tableId);
                if (resultCounter) resultCounter.textContent = visibleCount + ' genomes found';
            }}
            function exportToCSV(tableId, filename) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let csv = [];
                const headerCells = rows[0].getElementsByTagName('th');
                const headerRow = [];
                for (let cell of headerCells) headerRow.push(cell.textContent);
                csv.push(headerRow.join(','));
                for (let i = 1; i < rows.length; i++) {{
                    if (rows[i].style.display !== 'none') {{
                        const cells = rows[i].getElementsByTagName('td');
                        const row = [];
                        for (let cell of cells) {{
                            let text = cell.textContent.trim();
                            text = text.replace(/\\n/g, ' ').replace(/,/g, ';');
                            row.push(text);
                        }}
                        csv.push(row.join(','));
                    }}
                }}
                const blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url; a.download = filename;
                document.body.appendChild(a); a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            function exportToJSON() {{
                const data = {{
                    metadata: {{
                        total_genomes: {stats['total_genomes']},
                        successful_genomes: {stats['successful_genomes']},
                        date: '{self.metadata['analysis_date']}',
                        tool: '{self.metadata['tool_name']}',
                        version: '{self.metadata['version']}'
                    }},
                    statistics: {json.dumps(stats, indent=2)},
                    summary_rows: {json.dumps(summary_rows, indent=2)}
                }};
                const jsonStr = JSON.stringify(data, null, 2);
                const blob = new Blob([jsonStr], {{ type: 'application/json' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url; a.download = 'Kaptive_summary_data.json';
                document.body.appendChild(a); a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            function printReport() {{ window.print(); }}
            function sortTable(tableId, columnIndex) {{
                const table = document.getElementById(tableId);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                rows.sort((a, b) => {{
                    const aText = a.children[columnIndex].textContent.trim();
                    const bText = b.children[columnIndex].textContent.trim();
                    const aNum = parseFloat(aText);
                    const bNum = parseFloat(bText);
                    if (!isNaN(aNum) && !isNaN(bNum)) return aNum - bNum;
                    return aText.localeCompare(bText);
                }});
                rows.forEach(row => tbody.appendChild(row));
            }}
            function filterByStatus(status) {{
                const table = document.getElementById('detailed-results-table');
                const rows = table.getElementsByTagName('tr');
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const statusCell = row.children[1].textContent.trim();
                    if (status === 'all' || statusCell === status) {{
                        row.style.display = '';
                    }} else {{ row.style.display = 'none'; }}
                }}
            }}
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
        """

        type_distribution_html = ""
        if stats['common_k_types'] or stats['common_o_types']:
            blocks = ""
            if stats['common_k_types']:
                items = ""
                for info in stats['common_k_types']:
                    pct = (info['count'] / stats['genomes_with_k'] * 100
                           if stats['genomes_with_k'] > 0 else 0)
                    items += f"""
                        <li>
                            <span class="type-name">{info['type']}</span>
                            <span class="type-count">{info['count']} genomes ({pct:.1f}%)</span>
                        </li>
                    """
                blocks += f"""
                <div class="type-list">
                    <h3 style="color: #28a745;">Top K Locus Types</h3>
                    <ul>{items}</ul>
                </div>
                """
            if stats['common_o_types']:
                items = ""
                for info in stats['common_o_types']:
                    pct = (info['count'] / stats['genomes_with_o'] * 100
                           if stats['genomes_with_o'] > 0 else 0)
                    items += f"""
                        <li>
                            <span class="type-name">{info['type']}</span>
                            <span class="type-count">{info['count']} genomes ({pct:.1f}%)</span>
                        </li>
                    """
                blocks += f"""
                <div class="type-list">
                    <h3 style="color: #17a2b8;">Top O Locus Types</h3>
                    <ul>{items}</ul>
                </div>
                """
            type_distribution_html = f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Most Common Locus Types</h2>
            <div class="type-distribution">{blocks}</div>
        </div>
            """

        detail_rows = ""
        for row in summary_rows:
            row_class = "success-row" if row['status'] == 'Success' else "failed-row"
            error_note = f"Error: {row.get('error', '')}" if row.get('error') else ""
            details = error_note if error_note else (
                f"K: {row.get('k_details', 'N/A')} | O: {row.get('o_details', 'N/A')}"
            )
            status_color = '#28a745' if row['status'] == 'Success' else '#dc3545'
            detail_rows += f"""
                    <tr class="{row_class}">
                        <td><strong>{row['genome']}</strong></td>
                        <td><span style="color: {status_color}; font-weight: bold;">{row['status']}</span></td>
                        <td>{row['k_locus']}</td>
                        <td>{row['o_locus']}</td>
                        <td>{row['total_loci']}</td>
                        <td>{details}</td>
                    </tr>
            """

        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope Kaptive K/O Analysis - Summary Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff; padding: 20px; min-height: 100vh;
        }}
        .container {{ max-width: 1600px; margin: 0 auto; }}
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
            padding: 12px; text-align: left; border-bottom: 1px solid #e0e0e0;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white; font-weight: 600; cursor: pointer; position: relative;
        }}
        .gene-table th:hover {{
            background: linear-gradient(135deg, #2563eb 0%, #1e3a8a 100%);
        }}
        .gene-table th.sortable::after {{
            content: "\\2195"; position: absolute; right: 8px; opacity: 0.6;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .summary-stats {{
            display: flex; justify-content: space-around;
            margin: 20px 0; flex-wrap: wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white; padding: 20px; border-radius: 12px;
            text-align: center; box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px; flex: 1; min-width: 200px;
        }}
        .k-stat-card {{
            background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%);
            color: white; padding: 20px; border-radius: 12px;
            text-align: center; box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px; flex: 1; min-width: 200px;
        }}
        .o-stat-card {{
            background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%);
            color: white; padding: 20px; border-radius: 12px;
            text-align: center; box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px; flex: 1; min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1); color: white;
            padding: 20px; border-radius: 12px; margin: 20px 0;
            text-align: center; font-style: italic; border-left: 4px solid #fff;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8); color: white; padding: 30px;
            border-radius: 12px; margin-top: 40px;
        }}
        .footer a {{ color: #667eea; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .footer .timestamp {{ color: #fbbf24; font-weight: bold; margin-top: 10px; }}
        .interactive-controls {{
            background: #f8f9fa; padding: 15px; border-radius: 8px;
            margin: 15px 0; display: flex; flex-wrap: wrap;
            gap: 10px; align-items: center;
        }}
        .search-box {{ flex: 1; min-width: 200px; }}
        .search-box input {{
            width: 100%; padding: 8px 12px; border: 1px solid #ddd;
            border-radius: 4px; font-size: 14px;
        }}
        .filter-buttons {{ display: flex; gap: 5px; flex-wrap: wrap; }}
        .filter-buttons button {{
            padding: 6px 12px; background: #e9ecef; color: #495057;
            border: 1px solid #dee2e6; border-radius: 4px;
            cursor: pointer; font-size: 13px;
        }}
        .filter-buttons button.active {{
            background: #3b82f6; color: white; border-color: #3b82f6;
        }}
        .filter-buttons button:hover {{ background: #dee2e6; }}
        .export-buttons {{ display: flex; gap: 8px; flex-wrap: wrap; }}
        .export-buttons button {{
            padding: 8px 16px; background: #3b82f6; color: white;
            border: none; border-radius: 4px; cursor: pointer;
            font-size: 14px; transition: background 0.3s;
        }}
        .export-buttons button:hover {{ background: #2563eb; }}
        .export-buttons button.print {{ background: #10b981; }}
        .export-buttons button.print:hover {{ background: #059669; }}
        .export-buttons button.json {{ background: #8b5cf6; }}
        .export-buttons button.json:hover {{ background: #7c3aed; }}
        .result-counter {{
            font-size: 0.9em; color: #666; font-style: italic; margin-left: auto;
        }}
        .success-row {{ border-left: 4px solid #28a745; }}
        .failed-row {{ border-left: 4px solid #dc3545; background-color: #f8d7da; }}
        .type-distribution {{
            display: flex; flex-wrap: wrap; gap: 20px; margin-top: 20px;
        }}
        .type-list {{ flex: 1; min-width: 300px; }}
        .type-list ul {{ list-style-type: none; padding: 0; }}
        .type-list li {{
            padding: 8px 0; border-bottom: 1px solid #eee;
            display: flex; justify-content: space-between;
        }}
        .type-list .type-name {{ font-weight: 600; }}
        .type-list .type-count {{
            color: #666; background: #f8f9fa;
            padding: 2px 8px; border-radius: 12px; font-size: 0.9em;
        }}
        @media print {{
            .interactive-controls {{ display: none !important; }}
            body {{ background: white !important; color: black !important; }}
            .card {{ background: white !important; color: black !important; box-shadow: none !important; }}
            .gene-table {{ box-shadow: none !important; }}
        }}
    
        /* ════════ AcinetoScope uniformity override (auto-added) ════════ */
        .ascii-container {{ text-align: center !important; }}
        .ascii-art {{ display: inline-block !important; text-align: left !important; }}
        .footer {{ text-align: center !important; margin-top: 30px !important; padding: 20px !important; background: rgba(0, 0, 0, 0.3) !important; border-radius: 10px !important; font-size: 14px !important; overflow-wrap: break-word !important; }}
        .footer a {{ color: #fbbf24 !important; text-decoration: none !important; }}
        .footer a:hover {{ text-decoration: underline !important; }}
        .footer .timestamp {{ color: #fbbf24 !important; font-weight: bold !important; }}
        .footer .authorship {{ margin-top: 15px !important; padding: 15px !important; background: rgba(255, 255, 255, 0.1) !important; border-radius: 8px !important; font-size: 12px !important; line-height: 1.7 !important; }}
    </style>
    {interactive_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <pre class="ascii-art">{self.ascii_art}</pre>
            </div>
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">AcinetoScope Kaptive K/O Analysis - Summary Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive A. baumannii Capsule (K) and Lipooligosaccharide (O) Locus Analysis Across All Genomes</p>
                <p style="color: #666; font-size: 1.1em;">Kaptive {self.metadata['kaptive_version']} | Database: {self.metadata['database_version']}</p>
                <div class="interactive-controls">
                    <div class="search-box">
                        <input type="text" id="search-detailed-results"
                               placeholder="Search genomes, loci, or types..."
                               onkeyup="searchTable('detailed-results-table', this.value)">
                    </div>
                    <div class="filter-buttons">
                        <button class="active" onclick="filterByStatus('all')">All</button>
                        <button onclick="filterByStatus('Success')">Success</button>
                        <button onclick="filterByStatus('Failed')">Failed</button>
                    </div>
                    <div class="result-counter" id="result-counter-detailed-results-table">
                        {len(summary_rows)} genomes found
                    </div>
                    <div class="export-buttons">
                        <button onclick="exportToCSV('detailed-results-table', 'Kaptive_summary.csv')">Export CSV</button>
                        <button onclick="exportToJSON()" class="json">Export JSON</button>
                        <button onclick="printReport()" class="print">Print</button>
                    </div>
                </div>
            </div>
        </div>

        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Overall Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{stats['total_genomes']}</p>
                </div>
                <div class="stat-card">
                    <h3>Successful Analyses</h3>
                    <p style="font-size: 2em; margin: 0;">{stats['successful_genomes']}</p>
                    <p style="font-size: 0.9em; margin-top: 5px;">({(stats['successful_genomes']/stats['total_genomes']*100) if stats['total_genomes'] else 0:.1f}%)</p>
                </div>
                <div class="k-stat-card">
                    <h3>K Loci Found</h3>
                    <p style="font-size: 2em; margin: 0;">{stats['total_k_loci']}</p>
                    <p style="font-size: 0.9em; margin-top: 5px;">{stats['genomes_with_k']} genomes</p>
                </div>
                <div class="o-stat-card">
                    <h3>O Loci Found</h3>
                    <p style="font-size: 2em; margin: 0;">{stats['total_o_loci']}</p>
                    <p style="font-size: 0.9em; margin-top: 5px;">{stats['genomes_with_o']} genomes</p>
                </div>
            </div>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>Kaptive:</strong> {self.metadata['kaptive_version']}</p>
            <p><strong>Database:</strong> {self.metadata['database_version']}</p>
        </div>

        {type_distribution_html}

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Detailed Results by Genome</h2>
            <table class="gene-table" id="detailed-results-table">
                <thead>
                    <tr>
                        <th class="sortable" onclick="sortTable('detailed-results-table', 0)">Genome</th>
                        <th class="sortable" onclick="sortTable('detailed-results-table', 1)">Status</th>
                        <th class="sortable" onclick="sortTable('detailed-results-table', 2)">K Locus</th>
                        <th class="sortable" onclick="sortTable('detailed-results-table', 3)">O Locus</th>
                        <th class="sortable" onclick="sortTable('detailed-results-table', 4)">Total Loci</th>
                        <th>Details</th>
                    </tr>
                </thead>
                <tbody>
                    {detail_rows}
                </tbody>
            </table>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">Generated Files</h2>
            <ul style="color: #666; font-size: 1.1em;">
                <li><strong>Kaptive_summary.tsv</strong> - Complete Kaptive data for all genomes</li>
                <li><strong>Kaptive_summary.json</strong> - Single batch JSON with statistics and per-genome summaries</li>
                <li><strong>Per-sample HTML reports</strong> - One report per genome in <code>&lt;genome&gt;/&lt;genome&gt;_kaptive_report.html</code></li>
                <li><strong>This interactive summary report</strong> - Cross-genome analysis with search/filter</li>
            </ul>
        </div>

        {self._get_footer_html("AcinetoScope - Kaptive K/O Batch Summary")}
    </div>
</body>
</html>
"""

        html_file = os.path.join(output_base, "Kaptive_summary.html")
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        self.logger.info(f"Kaptive summary HTML report: {html_file}")


def main():
    parser = argparse.ArgumentParser(
        description='AcinetoScope Kaptive K/O Analysis - A. baumannii Capsule and Lipooligosaccharide Typing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python acineto_kaptive.py "*.fna"
  python acineto_kaptive.py "AB_*.fasta"
  python acineto_kaptive.py "*.fa" --output my_kaptive_results

Supported FASTA extensions: .fasta, .fa, .fna, .faa

PER-SAMPLE OUTPUTS (1 file each):
    <sample>/<sample>_kaptive_report.html

BATCH-LEVEL OUTPUTS:
    Kaptive_summary.tsv
    Kaptive_summary.json    (single batch JSON)
    Kaptive_summary.html
        """
    )
    parser.add_argument('pattern', help='File pattern for A. baumannii genomes (e.g., "*.fasta")')
    parser.add_argument('--output', '-o', default='kaptive_results',
                        help='Output directory (default: kaptive_results)')
    args = parser.parse_args()

    print("\n" + "=" * 80)
    executor = AcinetoKaptive()
    print(executor.ascii_art)
    print("=" * 80)
    print("A. baumannii Kaptive K/O Locus Analysis")
    print("=" * 80)
    print("Author: Brown Beckley | Email: brownbeckley94@gmail.com")
    print("Affiliation: University of Ghana Medical School - Department of Medical Biochemistry")
    print("=" * 80)

    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("AcinetoScope Kaptive K/O Analysis FINAL SUMMARY")
        executor.logger.info("=" * 50)

        total_hits = 0
        k_hits = 0
        o_hits = 0
        successful_genomes = 0

        for genome_name, result in results.items():
            if result['status'] == 'success':
                successful_genomes += 1
                total_hits += result['hit_count']
                k_hits += sum(1 for hit in result['hits'] if hit.get('Locus Type') == 'K')
                o_hits += sum(1 for hit in result['hits'] if hit.get('Locus Type') == 'O')
                executor.logger.info(f"{genome_name}: {result['hit_count']} K/O loci")

        executor.logger.info("\nA. BAUMANNII K/O SUMMARY STATISTICS:")
        executor.logger.info(f"   Total genomes processed: {len(results)}")
        executor.logger.info(f"   Successfully analyzed: {successful_genomes}")
        executor.logger.info(f"   Total K/O loci found: {total_hits}")
        executor.logger.info(f"   K loci (capsule): {k_hits}")
        executor.logger.info(f"   O loci (lipooligosaccharide): {o_hits}")

        executor.logger.info("\nPER-SAMPLE OUTPUTS (1 file each):")
        executor.logger.info("   <sample>/<sample>_kaptive_report.html")

        executor.logger.info("\nBATCH-LEVEL OUTPUTS:")
        executor.logger.info(f"   Kaptive summary TSV:  {args.output}/Kaptive_summary.tsv")
        executor.logger.info(f"   Master JSON (1):      {args.output}/Kaptive_summary.json")
        executor.logger.info(f"   Summary HTML:         {args.output}/Kaptive_summary.html")

        executor.logger.info(f"\n{random.choice(executor.science_quotes)}")

    except Exception as e:
        executor.logger.error(f"A. baumannii K/O analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()