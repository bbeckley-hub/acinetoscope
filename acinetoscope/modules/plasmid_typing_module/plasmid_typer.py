#!/usr/bin/env python3
"""
AcinetoScope Plasmid Typing Module - AcinetobacterPlasmidTyping (APT)
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-07-17

Uses APT database v3.0 (Feb 2025) for Acinetobacter plasmid rep typing.
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
from collections import defaultdict
import random

class AcinetoPlasmidTyper:
    def __init__(self, cpus: int = None):
        self.logger = self._setup_logging()
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        self.db_path = os.path.join(self.module_dir, "data", "apt_reps_v3.fasta")
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)
        self.min_identity = 95.0
        self.min_subject_coverage = 70.0   # coverage of the rep gene itself
        self.rep_pattern = re.compile(r'^>?([rR][1-9P]+[-_][Tt][0-9]+)', re.IGNORECASE)
        
        self.metadata = {
            "tool_name": "AcinetoScope Plasmid Typing",
            "version": "1.3.2",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "database": "APT v3.0 (Feb 2025)",
            "blast_version": self._get_blast_version(),
            "thresholds": f"identity ≥ {self.min_identity}%, subject coverage ≥ {self.min_subject_coverage}%"
        }
        
        # Exact same ASCII art as AMR module
        self.ascii_art = r"""
 █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
"""
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
    
    def _setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        return logging.getLogger(__name__)
    
    def _get_available_ram(self) -> float:
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except:
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
        except:
            return 4
    
    def _get_blast_version(self) -> str:
        try:
            res = subprocess.run(['blastn', '-version'], capture_output=True, text=True)
            return res.stdout.splitlines()[0] if res.stdout else "unknown"
        except:
            return "not found"
    
    def check_blast_installed(self) -> bool:
        try:
            subprocess.run(['blastn', '-version'], capture_output=True, check=True)
            return True
        except:
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
            rep = match.group(1).upper().replace('_', '-')
            return rep
        return "unknown"
    
    def run_blast_on_genome(self, genome_file: str, out_dir: str) -> Dict:
        genome = Path(genome_file).stem
        blast_out = os.path.join(out_dir, f"{genome}_blast.tab")
        cmd = [
            'blastn', '-query', genome_file, '-subject', self.db_path,
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
            '-out', blast_out,
            '-perc_identity', str(self.min_identity),
            '-evalue', '1e-5',
            '-num_threads', '1'
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BLAST failed {genome}: {e.stderr}")
            return {'genome': genome, 'hits': [], 'hit_count': 0, 'status': 'failed'}
        
        # Parse all hits (no deduplication)
        raw_hits = self._parse_blast_output(blast_out)
        # Filter by subject coverage (≥70% of rep gene length)
        filtered_hits = [h for h in raw_hits if h['subject_coverage'] >= self.min_subject_coverage]
        result = {
            'genome': genome,
            'hits': filtered_hits,
            'hit_count': len(filtered_hits),
            'status': 'success'
        }
        # Per‑sample TSV and HTML
        self._create_per_sample_tsv(genome, result, out_dir)
        self._create_per_sample_html(genome, result, out_dir)
        return result
    
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
                        'bitscore': float(bitscore)
                    })
        except Exception as e:
            self.logger.error(f"Parse error {blast_file}: {e}")
        return hits
    
    def _create_per_sample_tsv(self, genome: str, result: Dict, out_dir: str):
        tsv_file = os.path.join(out_dir, f"{genome}_plasmid_hits.tsv")
        with open(tsv_file, 'w') as f:
            f.write("Query_Contig\tRep_Type\tIdentity(%)\tSubject_Coverage(%)\tE-value\tBitscore\tSubject_Header\n")
            for h in result['hits']:
                f.write(f"{h['query_sequence']}\t{h['rep_type']}\t{h['identity']:.2f}\t{h['subject_coverage']:.1f}\t{h['evalue']:.1e}\t{h['bitscore']:.1f}\t{h['subject_header']}\n")
        self.logger.info(f"TSV for {genome}: {tsv_file}")
    
    def _create_per_sample_html(self, genome: str, result: Dict, out_dir: str):
        hits = result['hits']
        unique_reps = sorted(set(h['rep_type'] for h in hits))
        # Extract sample name from genome (remove path)
        sample_name = Path(genome).name
        
        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>AcinetoScope Plasmid Typing Report - {sample_name}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
            text-align: center;
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
            display: inline-block;
            text-align: left;
        }}
        .card {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        .stats {{
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
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .rep-badge {{
            display: inline-block;
            background: #e74c3c;
            color: white;
            padding: 5px 12px;
            border-radius: 20px;
            margin: 3px;
            font-size: 0.9em;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
            text-align: center;
        }}
        .footer a {{ color: #667eea; text-decoration: none; }}
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
        .interactive-controls {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            align-items: center;
        }}
        .search-box {{
            flex: 1;
            min-width: 200px;
        }}
        .search-box input {{
            width: 100%;
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        .export-buttons {{
            display: flex;
            gap: 8px;
            flex-wrap: wrap;
        }}
        .export-buttons button {{
            padding: 8px 16px;
            background: #3b82f6;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            transition: background 0.3s;
        }}
        .export-buttons button:hover {{ background: #2563eb; }}
        .export-buttons button.print {{ background: #10b981; }}
        .export-buttons button.print:hover {{ background: #059669; }}
        .result-counter {{
            font-size: 0.9em;
            color: #666;
            font-style: italic;
        }}
        .no-print {{ }}
        @media print {{
            .no-print, .interactive-controls, .quick-search-bar, .quote-container, .ascii-container {{ display: none !important; }}
            body {{ background: white !important; color: black !important; }}
            .card {{ background: white !important; color: black !important; box-shadow: none !important; }}
        }}
    </style>
    <script>
        function searchTable(tableId, term) {{
            const table = document.getElementById(tableId);
            const rows = table.getElementsByTagName('tr');
            let visible = 0;
            for (let i = 1; i < rows.length; i++) {{
                const text = rows[i].textContent.toLowerCase();
                if (text.includes(term.toLowerCase())) {{
                    rows[i].style.display = '';
                    visible++;
                }} else {{
                    rows[i].style.display = 'none';
                }}
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
        function printReport() {{
            const originalTitle = document.title;
            const printWindow = window.open('', '_blank');
            const content = document.querySelector('.container').cloneNode(true);
            const noPrint = content.querySelectorAll('.no-print, .interactive-controls, .quick-search-bar, .quote-container, .ascii-container');
            noPrint.forEach(el => el.remove());
            printWindow.document.write('<html><head><title>' + originalTitle + '</title><style>' +
                'body{{font-family:Arial;margin:20px;}} table{{border-collapse:collapse;width:100%;}} ' +
                'th,td{{border:1px solid #ddd;padding:8px;text-align:left;}} th{{background:#f2f2f2;}}' +
                '</style></head><body>' + content.innerHTML + '</body></html>');
            printWindow.document.close();
            printWindow.print();
        }}
        // Rotating quotes
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
        <div class="ascii-art">{self.ascii_art}</div>
    </div>
    
    <div class="quote-container">
        <div id="science-quote" style="font-size: 1.1em;"></div>
    </div>
    
    <div class="card">
        <h1 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🧬 AcinetoScope Plasmid Typing Report</h1>
        <p><strong>Sample Name:</strong> {sample_name}</p>
        <p><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
        <p><strong>Database:</strong> {self.metadata['database']}</p>
        <p><strong>Thresholds:</strong> {self.metadata['thresholds']}</p>
        <div class="stats">
            <div class="stat-card">Total Hits: {len(hits)}</div>
            <div class="stat-card">Unique Rep Types: {len(unique_reps)}</div>
        </div>
        <h2>Detected Rep Types</h2>
        <div>
"""
        for rt in unique_reps:
            html += f'<span class="rep-badge">{rt}</span> '
        html += """
        </div>
    </div>
"""
        if hits:
            html += f"""
    <div class="card">
        <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📋 All Plasmid Hits (Filtered by Coverage ≥ {self.min_subject_coverage}%)</h2>
        <div class="interactive-controls no-print">
            <div class="search-box">
                <input type="text" id="search-hits" placeholder="Search rep type, contig..." onkeyup="searchTable('hits-table', this.value)">
            </div>
            <div class="export-buttons">
                <button onclick="exportToCSV('hits-table', '{sample_name}_plasmid_hits.csv')">📥 Export CSV</button>
                <button onclick="printReport()" class="print">🖨️ Print</button>
            </div>
            <div class="result-counter" id="result-counter-hits-table">{len(hits)} results found</div>
        </div>
        <table id="hits-table">
            <thead>
                <tr><th>Query Contig</th><th>Rep Type</th><th>Identity (%)</th><th>Subject Coverage (%)</th><th>E-value</th><th>Bitscore</th></tr>
            </thead>
            <tbody>
"""
            for h in hits:
                html += f"""
                <tr>
                    <td>{h['query_sequence']}</td>
                    <td><strong>{h['rep_type']}</strong></td>
                    <td>{h['identity']:.2f}</td>
                    <td>{h['subject_coverage']:.1f}</td>
                    <td>{h['evalue']:.1e}</td>
                    <td>{h['bitscore']:.1f}</td>
                </tr>
"""
            html += """
            </tbody>
        </table>
    </div>
"""
        else:
            html += """
    <div class="card">
        <h2>✅ No Plasmid Rep Types Detected</h2>
        <p>No hits passed the identity and coverage thresholds.</p>
    </div>
"""
        html += f"""
    <div class="footer">
        <h3 style="color: #fff; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">👥 Contact Information</h3>
        <p><strong>Author:</strong> Brown Beckley</p>
        <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
        <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
        <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
        <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;"></p>
    </div>
</div>
</body>
</html>
"""
        html_file = os.path.join(out_dir, f"{sample_name}_plasmid_report.html")
        with open(html_file, 'w') as f:
            f.write(html)
    
    def process_single_genome(self, genome_file: str, out_base: str) -> Dict:
        genome = Path(genome_file).stem
        genome_dir = os.path.join(out_base, genome)
        os.makedirs(genome_dir, exist_ok=True)
        return self.run_blast_on_genome(genome_file, genome_dir)
    
    def process_multiple_genomes(self, pattern: str, out_base: str = "acineto_plasmid_results") -> Dict:
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
                    self.logger.info(f"✓ {res['genome']}: {res['hit_count']} rep hits")
                except Exception as e:
                    self.logger.error(f"✗ {f}: {e}")
                    results[Path(f).stem] = {'genome': Path(f).stem, 'hits': [], 'hit_count': 0, 'status': 'failed'}
        self._create_summary_tsv(results, out_base)
        self._create_summary_html(results, out_base)
        return results
    
    def _create_summary_tsv(self, results: Dict, out_base: str):
        tsv_file = os.path.join(out_base, "plasmid_summary.tsv")
        with open(tsv_file, 'w') as f:
            f.write("Genome\tHitCount\tRepTypes\n")
            for name, res in results.items():
                reps = ",".join(sorted(set(h['rep_type'] for h in res['hits'])))
                f.write(f"{name}\t{res['hit_count']}\t{reps}\n")
        self.logger.info(f"Summary TSV: {tsv_file}")
    
    def _create_summary_html(self, results: Dict, out_base: str):
        rep_freq = defaultdict(int)
        genome_reps = {}
        for name, res in results.items():
            reps = sorted(set(h['rep_type'] for h in res['hits']))
            genome_reps[name] = reps
            for r in reps:
                rep_freq[r] += 1
        
        html = f"""<!DOCTYPE html>
    <html>
    <head><meta charset="UTF-8"><title>AcinetoScope Plasmid Typing Summary</title>
    <style>
    * {{ margin:0; padding:0; box-sizing:border-box; }}
    body {{
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        color: #ffffff; padding: 20px; min-height: 100vh;
    }}
    .container {{ max-width: 1400px; margin: 0 auto; }}
    .ascii-container {{
        background: rgba(0, 0, 0, 0.7); padding: 20px; border-radius: 15px;
        border: 2px solid rgba(0, 255, 0, 0.3); margin-bottom: 20px; text-align: center;
    }}
    .ascii-art {{
        font-family: 'Courier New', monospace; font-size: 10px; color: #00ff00;
        white-space: pre; overflow-x: auto; display: inline-block; text-align: left;
    }}
    .card {{ background: rgba(255, 255, 255, 0.95); color: #1f2937; padding: 25px; border-radius: 12px; margin: 20px 0; }}
    .stats {{ display: flex; justify-content: space-around; margin: 20px 0; flex-wrap: wrap; }}
    .stat-card {{
        background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
        color: white; padding: 20px; border-radius: 12px; text-align: center;
        margin: 10px; flex: 1; min-width: 200px;
    }}
    table {{ width: 100%; border-collapse: collapse; margin: 20px 0; background: white; border-radius: 8px; overflow: hidden; }}
    th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #e0e0e0; vertical-align: top; }}
    th {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; }}
    .interactive-controls {{ background: #f8f9fa; padding: 15px; border-radius: 8px; display: flex; gap: 10px; align-items: center; flex-wrap: wrap; margin: 10px 0; }}
    .search-box input {{ padding: 8px; width: 250px; border: 1px solid #ddd; border-radius: 4px; }}
    button {{ padding: 8px 16px; background: #3b82f6; color: white; border: none; border-radius: 4px; cursor: pointer; }}
    .rep-badge {{ display: inline-block; background: #e74c3c; color: white; padding: 4px 10px; border-radius: 20px; margin: 2px; }}
    .footer {{ background: rgba(0, 0, 0, 0.8); padding: 30px; border-radius: 12px; margin-top: 40px; text-align: center; }}
    .quote-container {{ background: rgba(255, 255, 255, 0.1); padding: 20px; border-radius: 12px; margin: 20px 0; text-align: center; font-style: italic; }}
    </style>
    <script>
    function searchTable(tid, term) {{
        const t = document.getElementById(tid); const rows = t.getElementsByTagName('tr'); let vis = 0;
        for(let i=1;i<rows.length;i++) {{
            const txt = rows[i].textContent.toLowerCase();
            if(txt.includes(term.toLowerCase())) {{ rows[i].style.display = ''; vis++; }}
            else {{ rows[i].style.display = 'none'; }}
        }}
        const cnt = document.getElementById('result-counter');
        if(cnt) cnt.innerText = vis + ' results found';
    }}
    function exportToCSV(tid, fname) {{
        const t = document.getElementById(tid); let csv = [];
        const headers = t.querySelectorAll('th');
        let headerRow = []; headers.forEach(h => headerRow.push(h.textContent));
        csv.push(headerRow.join(','));
        const rows = t.querySelectorAll('tbody tr');
        rows.forEach(row => {{
            if(row.style.display !== 'none') {{
                const cells = row.querySelectorAll('td');
                let rowData = [];
                cells.forEach(c => rowData.push(c.textContent.trim().replace(/,/g,';')));
                csv.push(rowData.join(','));
            }}
        }});
        const blob = new Blob([csv.join('\\n')], {{type:'text/csv'}});
        const a = document.createElement('a');
        a.href = URL.createObjectURL(blob); a.download = fname;
        document.body.appendChild(a); a.click(); document.body.removeChild(a);
        URL.revokeObjectURL(a.href);
    }}
    function printReport() {{
        const win = window.open('', '_blank');
        const content = document.querySelector('.container').cloneNode(true);
        content.querySelectorAll('.no-print, .interactive-controls, .quote-container, .ascii-container').forEach(el => el.remove());
        win.document.write('<html><head><title>Plasmid Summary</title><style>body{{font-family:Arial;margin:20px;}}table{{border-collapse:collapse;width:100%;}}th,td{{border:1px solid #ddd;padding:8px;}}</style></head><body>' + content.innerHTML + '</body></html>');
        win.document.close(); win.print();
    }}
    let quotes = {json.dumps(self.science_quotes)};
    let currentQuote = 0;
    function rotateQuote() {{
        document.getElementById('science-quote').innerHTML = quotes[currentQuote];
        currentQuote = (currentQuote + 1) % quotes.length;
    }}
    document.addEventListener('DOMContentLoaded', function() {{
        rotateQuote(); setInterval(rotateQuote, 10000);
    }});
    </script>
    </head>
    <body>
    <div class="container">
    <div class="ascii-container"><div class="ascii-art">{self.ascii_art}</div></div>
    <div class="quote-container"><div id="science-quote"></div></div>
    <div class="card">
    <h1>📊 AcinetoScope Plasmid Typing Summary</h1>
    <p>Database: {self.metadata['database']} | Thresholds: {self.metadata['thresholds']}</p>
    <div class="stats">
    <div class="stat-card">Genomes: {len(results)}</div>
    <div class="stat-card">With Plasmids: {sum(1 for r in results.values() if r['hit_count']>0)}</div>
    <div class="stat-card">Unique Rep Types: {len(rep_freq)}</div>
    </div>
    </div>
    <div class="card">
    <h2>📈 Rep Type Frequency</h2>
    <div class="interactive-controls no-print">
    <div class="search-box"><input type="text" id="freqSearch" placeholder="Search rep type..." onkeyup="searchTable('freqTable',this.value)"></div>
    <button onclick="exportToCSV('freqTable','rep_frequency.csv')">📥 CSV</button>
    <button onclick="printReport()">🖨️ Print</button>
    <span id="result-counter"></span>
    </div>
    <table id="freqTable"><thead><tr><th>Rep Type</th><th>Count</th><th>Genomes</th></tr></thead><tbody>
    """
        for rt, cnt in sorted(rep_freq.items(), key=lambda x: x[1], reverse=True):
            glist = [name for name, reps in genome_reps.items() if rt in reps]
            # No truncation – show all genome names
            html += f"<tr><td><strong>{rt}</strong></td><td>{cnt}</td><td>{', '.join(glist)}</td></tr>"
        html += """
    </tbody></table></div>
    <div class="card"><h2>🧬 Per‑Genome Rep Types</h2>
    <div class="interactive-controls no-print">
    <input type="text" id="genomeSearch" placeholder="Search genome..." onkeyup="searchTable('genomeTable',this.value)">
    <button onclick="exportToCSV('genomeTable','genome_reps.csv')">📥 CSV</button>
    </div>
    <table id="genomeTable"><thead><tr><th>Genome</th><th>Rep Types</th></tr></thead><tbody>
    """
        for name, reps in sorted(genome_reps.items()):
            if reps:
                badges = ' '.join(f'<span class="rep-badge">{r}</span>' for r in reps)
                html += f"<tr><td>{name}</td><td>{badges}</td></tr>"
            else:
                html += f"<tr><td>{name}</td><td><em>No Rep found</em></td></tr>"
        html += f"""
    </tbody></table></div>
    <div class="footer">
    <h3 style="color:#fff; border-bottom:2px solid #3b82f6; padding-bottom:10px;">👥 Contact Information</h3>
    <p><strong>Author:</strong> Brown Beckley</p>
    <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
    <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">bbeckley-hub</a></p>
    <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
    <p>Analysis date: {self.metadata['analysis_date']}</p>
    </div>
    </div>
    <script>document.getElementById('result-counter').innerText = document.querySelectorAll('#freqTable tbody tr').length + ' results found';</script>
    </body>
    </html>
    """
        html_file = os.path.join(out_base, "plasmid_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html)
        self.logger.info(f"Summary HTML: {html_file}")

def main():
    parser = argparse.ArgumentParser(description="AcinetoScope Plasmid Typing (APT)")
    parser.add_argument('pattern', help='FASTA file pattern (e.g., "*.fna")')
    parser.add_argument('--cpus', '-c', type=int, help='Number of CPU cores')
    parser.add_argument('--output', '-o', default='acineto_plasmid_results', help='Output directory')
    args = parser.parse_args()
    
    typer = AcinetoPlasmidTyper(cpus=args.cpus)
    if not typer.check_blast_installed():
        sys.exit(1)
    if not typer.check_database():
        sys.exit(1)
    
    results = typer.process_multiple_genomes(args.pattern, args.output)
    print(f"\n✅ Completed. Results in {args.output}")
    print(f"   Summary: {args.output}/plasmid_summary_report.html")
    print(f"   Per‑sample TSV files are in each genome subdirectory.")
    print(f"   Per‑sample HTML reports are also in each genome subdirectory.")

if __name__ == "__main__":
    main()