#!/usr/bin/env python3
"""
AcinetoScope Kaptive K/O Locus Analysis - A. baumannii Capsule and Lipooligosaccharide Typing
Comprehensive Kaptive analysis for A. baumannii with beautiful HTML reporting
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2025-12-28
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import json
import pandas as pd
import csv
from collections import defaultdict

class AcinetoKaptive:
    """Kaptive executor for A. baumannii K/O locus typing with comprehensive HTML reporting"""
    
    def __init__(self):
        # Setup logging
        self.logger = self._setup_logging()
        
        # Kaptive database keywords for A. baumannii
        self.k_db_keyword = "ab_k"  # K locus database
        self.o_db_keyword = "ab_o"  # O locus database
        
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
            "tool_name": "AcinetoScope Kaptive K/O Analysis",
            "version": "1.0.0", 
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "kaptive_version": self._get_kaptive_version(),
            "database_version": "2025-12-26 (ab_k/ab_o)"
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
    
    def _setup_logging(self):
        """Setup logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def _get_kaptive_version(self) -> str:
        """Get Kaptive version"""
        try:
            result = subprocess.run(
                ["kaptive", "--version"], 
                capture_output=True, 
                text=True, 
                check=True
            )
            version_line = result.stdout.strip()
            # Extract version number
            match = re.search(r'v?(\d+\.\d+\.\d+)', version_line)
            if match:
                return match.group(1)
            return version_line
        except Exception as e:
            self.logger.warning(f"Could not determine Kaptive version: {e}")
            return "Unknown"

    def check_kaptive_installed(self) -> bool:
        """Check if Kaptive is available"""
        try:
            # Test Kaptive version
            result = subprocess.run(
                ["kaptive", "--version"], 
                capture_output=True, 
                text=True, 
                check=True
            )
            
            version_line = result.stdout.strip()
            self.logger.info(f"Kaptive version: {version_line}")
            
            return True
            
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"Kaptive check failed: {e}")
            self.logger.info("Please install Kaptive: conda install -c bioconda kaptive")
            self.logger.info("Or: pip install kaptive")
            return False

    def run_kaptive_single_genome(self, genome_file: str, output_dir: str) -> Dict[str, Any]:
        """Run Kaptive on a single A. baumannii genome for both K and O loci"""
        genome_name = Path(genome_file).stem
        output_base = os.path.join(output_dir, genome_name)
        os.makedirs(output_base, exist_ok=True)
        
        k_tsv = os.path.join(output_base, f"{genome_name}_K.tsv")
        o_tsv = os.path.join(output_base, f"{genome_name}_O.tsv")
        combined_tsv = os.path.join(output_base, f"{genome_name}_combined.tsv")
        
        self.logger.info(f"Running Kaptive K locus analysis on {genome_name}")
        
        try:
            # Run K locus analysis
            cmd_k = [
                "kaptive", "assembly",
                self.k_db_keyword,
                genome_file,
                "-o", k_tsv,
                "--verbose"
            ]
            
            result_k = subprocess.run(
                cmd_k,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Run O locus analysis
            self.logger.info(f"Running Kaptive O locus analysis on {genome_name}")
            
            cmd_o = [
                "kaptive", "assembly",
                self.o_db_keyword,
                genome_file,
                "-o", o_tsv,
                "--verbose"
            ]
            
            result_o = subprocess.run(
                cmd_o,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Combine K and O results
            self._combine_k_o_results(k_tsv, o_tsv, combined_tsv)
            
            # Parse the combined results
            hits = self._parse_kaptive_output(combined_tsv)
            
            # Create reports
            self._create_kaptive_html_report(genome_name, hits, output_base)
            self._create_kaptive_json_report(genome_name, hits, output_base)
            
            # Try to create Excel report if pandas is available
            try:
                import pandas as pd
                self._create_kaptive_excel_report(genome_name, hits, output_base)
            except ImportError:
                self.logger.warning("pandas not installed, skipping Excel report")
            
            return {
                'genome': genome_name,
                'output_file': combined_tsv,
                'hits': hits,
                'hit_count': len(hits),
                'status': 'success',
                'logs': {
                    'k_stdout': result_k.stdout[:1000] if result_k.stdout else "",
                    'k_stderr': result_k.stderr[:1000] if result_k.stderr else "",
                    'o_stdout': result_o.stdout[:1000] if result_o.stdout else "",
                    'o_stderr': result_o.stderr[:1000] if result_o.stderr else ""
                }
            }
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Kaptive failed for {genome_name}: {e}")
            return {
                'genome': genome_name,
                'output_file': combined_tsv,
                'hits': [],
                'hit_count': 0,
                'status': 'failed',
                'error': str(e)
            }
    
    def _combine_k_o_results(self, k_tsv: str, o_tsv: str, combined_tsv: str):
        """Combine K and O results into a single TSV file"""
        try:
            # Read K results
            k_lines = []
            if os.path.exists(k_tsv):
                with open(k_tsv, 'r') as f:
                    k_lines = f.readlines()
            
            # Read O results
            o_lines = []
            if os.path.exists(o_tsv):
                with open(o_tsv, 'r') as f:
                    o_lines = f.readlines()
            
            # Write combined results
            with open(combined_tsv, 'w') as f:
                # Write header from K results (they should be the same)
                if k_lines:
                    f.write(k_lines[0])
                    # Write K results (skip header)
                    for line in k_lines[1:]:
                        f.write(line)
                # Write O results (skip header)
                if o_lines:
                    for line in o_lines[1:]:
                        f.write(line)
            
            self.logger.info(f"Combined K and O results to {combined_tsv}")
            
        except Exception as e:
            self.logger.error(f"Error combining K/O results: {e}")
    
    def _parse_kaptive_output(self, kaptive_file: str) -> List[Dict]:
        """Parse Kaptive output file into structured data"""
        hits = []
        try:
            if not os.path.exists(kaptive_file):
                return hits
                
            with open(kaptive_file, 'r') as f:
                lines = f.readlines()
                
            if not lines or len(lines) < 2:
                return hits
                
            # Parse header
            headers = lines[0].strip().split('\t')
            num_headers = len(headers)
            
            # Parse data lines
            for line_num, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                
                # Handle missing columns by padding with empty strings
                if len(parts) < num_headers:
                    parts.extend([''] * (num_headers - len(parts)))
                elif len(parts) > num_headers:
                    # If more parts than headers, truncate (shouldn't happen)
                    parts = parts[:num_headers]
                
                # Create hit with all headers
                hit = {}
                for i, header in enumerate(headers):
                    hit[header] = parts[i] if i < len(parts) else ''
                
                # Add locus type (K or O)
                locus = hit.get('Best match locus', '')
                if 'KL' in locus:
                    hit['Locus Type'] = 'K'
                elif 'OCL' in locus:
                    hit['Locus Type'] = 'O'
                else:
                    hit['Locus Type'] = 'Unknown'
                
                # Parse gene details if present
                gene_details = hit.get('Expected genes in locus, details', '')
                hit['Gene Details Parsed'] = self._parse_gene_details(gene_details)
                
                hits.append(hit)
                    
        except Exception as e:
            self.logger.error(f"Error parsing {kaptive_file}: {e}")
            
        self.logger.info(f"Parsed {len(hits)} Kaptive hits from {kaptive_file}")
        return hits
    
    def _parse_gene_details(self, gene_details_str: str) -> List[Dict]:
        """Parse gene details string into structured list"""
        genes = []
        try:
            if not gene_details_str:
                return genes
                
            # Split by semicolon
            gene_parts = gene_details_str.split(';')
            for gene_part in gene_parts:
                if gene_part:
                    # Format: gene_name,identity%,coverage%
                    subparts = gene_part.split(',')
                    if len(subparts) >= 3:
                        gene = {
                            'name': subparts[0],
                            'identity': subparts[1].replace('%', '') if len(subparts) > 1 else '',
                            'coverage': subparts[2].replace('%', '') if len(subparts) > 2 else ''
                        }
                        genes.append(gene)
        except Exception as e:
            self.logger.warning(f"Could not parse gene details: {e}")
        
        return genes
    
    def _create_kaptive_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Create comprehensive HTML report for Kaptive results"""
        
        # JavaScript for basic features (only for individual reports)
        basic_js = f"""
        <script>
            // Simple table sorting
            function sortTable(tableId, columnIndex) {{
                const table = document.getElementById(tableId);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                
                rows.sort((a, b) => {{
                    const aText = a.children[columnIndex].textContent.trim();
                    const bText = b.children[columnIndex].textContent.trim();
                    
                    // Try to convert to number if possible
                    const aNum = parseFloat(aText.replace('%', ''));
                    const bNum = parseFloat(bText.replace('%', ''));
                    
                    if (!isNaN(aNum) && !isNaN(bNum)) {{
                        return aNum - bNum;
                    }}
                    
                    return aText.localeCompare(bText);
                }});
                
                // Clear and re-add rows
                rows.forEach(row => tbody.appendChild(row));
            }}
            
            // Print report
            function printReport() {{
                window.print();
            }}
            
            // Rotating quotes
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            
            // Initialize
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
                setInterval(rotateQuote, 10000);
            }});
        </script>
        """
        
        # Store data for JSON export
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
            
            // Simple export to JSON
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
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope Kaptive K/O Analysis Report</title>
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
            max-width: 1800px;
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
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            font-size: 12px;
        }}
        
        .gene-table th, .gene-table td {{ 
            padding: 12px 8px; 
            text-align: left; 
            border-bottom: 1px solid #e0e0e0; 
            max-width: 300px;
            overflow-wrap: break-word;
        }}
        
        .gene-table th {{ 
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            font-weight: 600;
            cursor: pointer;
            position: relative;
        }}
        
        .gene-table th:hover {{ 
            background: linear-gradient(135deg, #2563eb 0%, #1e3a8a 100%);
        }}
        
        .gene-table th.sortable::after {{
            content: "↕";
            position: absolute;
            right: 8px;
            opacity: 0.6;
        }}
        
        tr:hover {{ background-color: #f8f9fa; }}
        .success {{ color: #28a745; font-weight: 600; }}
        .warning {{ color: #ffc107; font-weight: 600; }}
        .error {{ color: #dc3545; font-weight: 600; }}
        
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
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        
        .k-stat-card {{
            background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        
        .o-stat-card {{
            background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
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
        
        .k-locus-row {{ background-color: #d4edda; border-left: 4px solid #28a745; }}
        .o-locus-row {{ background-color: #d1ecf1; border-left: 4px solid #17a2b8; }}
        
        /* Controls */
        .controls {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            align-items: center;
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
            font-size: 14px;
            transition: background 0.3s;
        }}
        
        .export-buttons button:hover {{
            background: #2563eb;
        }}
        
        .export-buttons button.print {{
            background: #10b981;
        }}
        
        .export-buttons button.print:hover {{
            background: #059669;
        }}
        
        .export-buttons button.json {{
            background: #8b5cf6;
        }}
        
        .export-buttons button.json:hover {{
            background: #7c3aed;
        }}
        
        .gene-details-container {{
            max-height: 300px;
            overflow-y: auto;
            padding: 10px;
            background: #f8f9fa;
            border-radius: 5px;
            border: 1px solid #e0e0e0;
            margin-top: 10px;
            font-size: 11px;
        }}
        
        .gene-item {{
            display: inline-block;
            margin: 2px;
            padding: 3px 8px;
            background: #e9ecef;
            border-radius: 3px;
            font-size: 0.85em;
            white-space: nowrap;
        }}
        
        /* Responsive table */
        .table-container {{
            overflow-x: auto;
            margin: 20px 0;
        }}
        
        @media (max-width: 1200px) {{
            .gene-table {{
                font-size: 11px;
            }}
            
            .gene-table th, .gene-table td {{
                padding: 8px 6px;
            }}
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
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 AcinetoScope Kaptive K/O Analysis Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive A. baumannii Capsule (K) and Lipooligosaccharide (O) Locus Analysis</p>
                
                <div class="controls">
                    <div class="export-buttons">
                        <button onclick="exportToJSON()" class="json">📥 Export JSON</button>
                        <button onclick="printReport()" class="print">🖨️ Print Report</button>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""
        
        # Count K and O loci
        k_count = sum(1 for hit in hits if hit.get('Locus Type') == 'K')
        o_count = sum(1 for hit in hits if hit.get('Locus Type') == 'O')
        
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📊 A. baumannii K/O Locus Summary</h2>
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
"""
        
        # Detailed Kaptive results table
        if hits:
            html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🔬 Complete Kaptive Analysis Results</h2>
            
            <div class="table-container">
                <table class="gene-table" id="kaptive-results-table">
                    <thead>
                        <tr>
"""
            
            # Get all column headers from the first hit
            if hits:
                headers = list(hits[0].keys())
                # Remove parsed gene details from display
                headers_to_show = [h for h in headers if h not in ['Gene Details Parsed']]
                
                for i, header in enumerate(headers_to_show):
                    html_content += f'<th class="sortable" onclick="sortTable(\'kaptive-results-table\', {i})">{header}</th>'
                
                html_content += """
                        </tr>
                    </thead>
                    <tbody>
"""
            
                for hit in hits:
                    locus_type = hit.get('Locus Type', '')
                    row_class = "k-locus-row" if locus_type == 'K' else "o-locus-row" if locus_type == 'O' else ""
                    
                    html_content += f'<tr class="{row_class}">'
                    for header in headers_to_show:
                        value = hit.get(header, '')
                        # Truncate very long values
                        if len(str(value)) > 100:
                            value = str(value)[:100] + '...'
                        html_content += f'<td>{value}</td>'
                    html_content += '</tr>'
                
                html_content += """
                    </tbody>
                </table>
            </div>
            
            <h3 style="margin-top: 30px; color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🧬 Detailed Gene Information</h3>
"""
            
            # Show gene details for each hit
            for i, hit in enumerate(hits):
                locus = hit.get('Best match locus', '')
                locus_type = hit.get('Locus Type', '')
                gene_details = hit.get('Gene Details Parsed', [])
                
                if gene_details:
                    html_content += f"""
            <div class="gene-details-container">
                <h4 style="margin-top: 0; color: {'#28a745' if locus_type == 'K' else '#17a2b8'};">{locus} ({locus_type} Locus) - Gene Details:</h4>
                <div style="margin-top: 10px;">
"""
                    
                    for gene in gene_details:
                        html_content += f'<span class="gene-item">{gene["name"]}: {gene["identity"]}% identity, {gene["coverage"]}% coverage</span>'
                    
                    html_content += """
                </div>
            </div>
"""
            
            html_content += """
        </div>
"""
        else:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">⚠️ No K/O Loci Detected</h2>
            <p>No K or O loci found in this A. baumannii genome.</p>
        </div>
"""
        
        # Footer
        html_content += f"""
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using AcinetoScope Kaptive K/O Analysis v{self.metadata['version']}
                with Kaptive {self.metadata['kaptive_version']} and A. baumannii databases
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write HTML report
        html_file = os.path.join(output_dir, f"{genome_name}_kaptive_report.html")
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        self.logger.info(f"A. baumannii Kaptive HTML report generated: {html_file}")
    
    def _create_kaptive_json_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Create JSON report for Kaptive results"""
        json_data = {
            'metadata': {
                'genome': genome_name,
                'analysis_date': self.metadata['analysis_date'],
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'kaptive_version': self.metadata['kaptive_version'],
                'database_version': self.metadata['database_version']
            },
            'summary': {
                'total_loci': len(hits),
                'k_loci': sum(1 for hit in hits if hit.get('Locus Type') == 'K'),
                'o_loci': sum(1 for hit in hits if hit.get('Locus Type') == 'O'),
                'k_types': list(set(hit.get('Best match type', '') for hit in hits if hit.get('Locus Type') == 'K')),
                'o_types': list(set(hit.get('Best match type', '') for hit in hits if hit.get('Locus Type') == 'O'))
            },
            'hits': hits
        }
        
        json_file = os.path.join(output_dir, f"{genome_name}_kaptive_report.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, ensure_ascii=False)
        
        self.logger.info(f"A. baumannii Kaptive JSON report generated: {json_file}")
    
    def _create_kaptive_excel_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Create Excel report for Kaptive results"""
        try:
            # Create DataFrame from hits
            df = pd.DataFrame(hits)
            
            # Create Excel writer
            excel_file = os.path.join(output_dir, f"{genome_name}_kaptive_report.xlsx")
            
            # Remove 'Gene Details Parsed' column if it exists
            if 'Gene Details Parsed' in df.columns:
                df = df.drop(columns=['Gene Details Parsed'])
            
            df.to_excel(excel_file, index=False, engine='openpyxl')
            
            self.logger.info(f"A. baumannii Kaptive Excel report generated: {excel_file}")
            
        except Exception as e:
            self.logger.error(f"Could not create Excel report: {e}")
            # Create a simple CSV as fallback
            if hits:
                csv_file = os.path.join(output_dir, f"{genome_name}_kaptive_report.csv")
                df = pd.DataFrame(hits)
                if 'Gene Details Parsed' in df.columns:
                    df = df.drop(columns=['Gene Details Parsed'])
                df.to_csv(csv_file, index=False)
                self.logger.info(f"Created CSV report instead: {csv_file}")
    
    def process_single_genome(self, genome_file: str, output_base: str = "kaptive_results") -> Dict[str, Any]:
        """Process a single A. baumannii genome with Kaptive"""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        
        self.logger.info(f"=== PROCESSING A. BAUMANNII GENOME: {genome_name} ===")
        
        # Create output directory
        os.makedirs(results_dir, exist_ok=True)
        
        # Check Kaptive before running
        if not self.check_kaptive_installed():
            self.logger.error("Kaptive not available!")
            return {
                'genome': genome_name,
                'hits': [],
                'hit_count': 0,
                'status': 'failed',
                'error': 'Kaptive not available'
            }
        
        # Run Kaptive
        result = self.run_kaptive_single_genome(genome_file, results_dir)
        
        status_icon = "✓" if result['status'] == 'success' else "✗"
        self.logger.info(f"{status_icon} {genome_name}: {result['hit_count']} K/O loci found")
        
        return result
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "kaptive_results") -> Dict[str, Any]:
        """Process multiple A. baumannii genomes using wildcard pattern"""
        
        # Find genome files (support all FASTA extensions)
        fasta_extensions = ['.fna', '.fasta', '.fa', '.faa']
        genome_files = []
        
        # Try the pattern as-is first
        genome_files.extend(glob.glob(genome_pattern))
        
        # If no files found, try with extensions
        if not genome_files:
            for ext in fasta_extensions:
                genome_files.extend(glob.glob(f"{genome_pattern}{ext}"))
        
        # Remove duplicates and non-existent files
        genome_files = list(set([f for f in genome_files if os.path.exists(f)]))
        
        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")
        
        self.logger.info(f"Found {len(genome_files)} A. baumannii genomes: {[Path(f).name for f in genome_files]}")
        
        # Create output directory
        os.makedirs(output_base, exist_ok=True)
        
        # Process genomes sequentially (no threading)
        all_results = {}
        
        for genome_file in genome_files:
            genome_name = Path(genome_file).stem
            try:
                result = self.process_single_genome(genome_file, output_base)
                all_results[genome_name] = result
                self.logger.info(f"✓ COMPLETED: {genome_name} ({result['hit_count']} K/O loci)")
            except Exception as e:
                self.logger.error(f"✗ FAILED: {genome_name} - {e}")
                all_results[genome_name] = {
                    'genome': genome_name,
                    'hits': [],
                    'hit_count': 0,
                    'status': 'failed',
                    'error': str(e)
                }
        
        # Create Kaptive summary files after processing all genomes
        self.create_kaptive_summary(all_results, output_base)
        
        self.logger.info("=== A. BAUMANNII K/O ANALYSIS COMPLETE ===")
        self.logger.info(f"Processed {len(all_results)} genomes")
        self.logger.info(f"Results saved to: {output_base}")
        
        return all_results
    
    def create_kaptive_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create comprehensive Kaptive summary files and HTML reports for all A. baumannii samples"""
        self.logger.info("Creating A. baumannii Kaptive summary files and HTML reports...")
        
        # Create TSV summary file
        summary_file = os.path.join(output_base, "Kaptive_summary.tsv")
        
        # Get all headers from successful results
        all_headers = set()
        successful_hits = []
        
        for genome_name, result in all_results.items():
            if result.get('status') == 'success' and result.get('hits'):
                for hit in result['hits']:
                    # Add genome name to hit
                    hit_with_genome = hit.copy()
                    hit_with_genome['Genome'] = genome_name
                    successful_hits.append(hit_with_genome)
                    all_headers.update(hit.keys())
        
        # Write TSV summary
        if successful_hits:
            headers = ['Genome'] + [h for h in all_headers if h not in ['Genome', 'Gene Details Parsed']]
            
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write('\t'.join(headers) + '\n')
                
                for hit in successful_hits:
                    row = []
                    for header in headers:
                        value = hit.get(header, '')
                        # Replace tabs in values to avoid breaking TSV format
                        if isinstance(value, str):
                            value = value.replace('\t', ' ').replace('\n', ' ')
                        row.append(str(value))
                    f.write('\t'.join(row) + '\n')
        
        self.logger.info(f"✓ A. baumannii Kaptive summary file created: {summary_file}")
        
        # Create JSON summary
        self.create_json_summary(all_results, output_base)
        
        # Create comprehensive HTML summary report with interactive features
        self._create_summary_html_report(all_results, output_base)
    
    def create_json_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create JSON summary file"""
        master_summary = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'kaptive_version': self.metadata['kaptive_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results)
            },
            'genome_summaries': {},
            'statistics': self._create_summary_statistics(all_results)
        }
        
        for genome_name, result in all_results.items():
            master_summary['genome_summaries'][genome_name] = {
                'status': result.get('status', 'unknown'),
                'hit_count': result.get('hit_count', 0),
                'k_loci': sum(1 for hit in result.get('hits', []) if hit.get('Locus Type') == 'K'),
                'o_loci': sum(1 for hit in result.get('hits', []) if hit.get('Locus Type') == 'O'),
                'k_types': list(set(hit.get('Best match type', '') for hit in result.get('hits', []) if hit.get('Locus Type') == 'K')),
                'o_types': list(set(hit.get('Best match type', '') for hit in result.get('hits', []) if hit.get('Locus Type') == 'O')),
                'error': result.get('error', '') if result.get('status') == 'failed' else ''
            }
        
        json_file = os.path.join(output_base, "Kaptive_summary.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(master_summary, f, indent=2, ensure_ascii=False)
        
        self.logger.info(f"✓ JSON summary created: {json_file}")
    
    def _create_summary_statistics(self, all_results: Dict[str, Any]) -> Dict[str, Any]:
        """Create summary statistics"""
        stats = {
            'total_genomes': len(all_results),
            'successful_genomes': sum(1 for r in all_results.values() if r.get('status') == 'success'),
            'failed_genomes': sum(1 for r in all_results.values() if r.get('status') == 'failed'),
            'total_k_loci': 0,
            'total_o_loci': 0,
            'k_type_distribution': defaultdict(int),
            'o_type_distribution': defaultdict(int),
            'genomes_with_k': 0,
            'genomes_with_o': 0,
            'genomes_with_both': 0,
            'common_k_types': [],
            'common_o_types': []
        }
        
        for genome_name, result in all_results.items():
            if result.get('status') == 'success':
                hits = result.get('hits', [])
                k_hits = [h for h in hits if h.get('Locus Type') == 'K']
                o_hits = [h for h in hits if h.get('Locus Type') == 'O']
                
                stats['total_k_loci'] += len(k_hits)
                stats['total_o_loci'] += len(o_hits)
                
                if k_hits:
                    stats['genomes_with_k'] += 1
                    for hit in k_hits:
                        k_type = hit.get('Best match type', 'Unknown')
                        stats['k_type_distribution'][k_type] += 1
                
                if o_hits:
                    stats['genomes_with_o'] += 1
                    for hit in o_hits:
                        o_type = hit.get('Best match type', 'Unknown')
                        stats['o_type_distribution'][o_type] += 1
                
                if k_hits and o_hits:
                    stats['genomes_with_both'] += 1
        
        # Get top 5 most common K and O types
        if stats['k_type_distribution']:
            common_k = sorted(stats['k_type_distribution'].items(), key=lambda x: x[1], reverse=True)[:5]
            stats['common_k_types'] = [{'type': k, 'count': v} for k, v in common_k]
        
        if stats['o_type_distribution']:
            common_o = sorted(stats['o_type_distribution'].items(), key=lambda x: x[1], reverse=True)[:5]
            stats['common_o_types'] = [{'type': o, 'count': v} for o, v in common_o]
        
        # Convert defaultdict to regular dict for JSON serialization
        stats['k_type_distribution'] = dict(stats['k_type_distribution'])
        stats['o_type_distribution'] = dict(stats['o_type_distribution'])
        
        return stats
    
    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        """Create comprehensive HTML summary report with interactive features"""
        
        # Calculate statistics
        stats = self._create_summary_statistics(all_results)
        
        # Prepare data for tables
        summary_rows = []
        
        for genome_name, result in all_results.items():
            if result.get('status') == 'success':
                hits = result.get('hits', [])
                k_hits = [h for h in hits if h.get('Locus Type') == 'K']
                o_hits = [h for h in hits if h.get('Locus Type') == 'O']
                
                k_info = "None"
                if k_hits:
                    k_hit = k_hits[0]  # Take first K hit
                    k_type = k_hit.get('Best match type', 'Unknown')
                    identity = k_hit.get('Identity', 'N/A')
                    coverage = k_hit.get('Coverage', 'N/A')
                    k_info = f"{k_type} ({identity}, {coverage})"
                
                o_info = "None"
                if o_hits:
                    o_hit = o_hits[0]  # Take first O hit
                    o_type = o_hit.get('Best match type', 'Unknown')
                    identity = o_hit.get('Identity', 'N/A')
                    coverage = o_hit.get('Coverage', 'N/A')
                    o_info = f"{o_type} ({identity}, {coverage})"
                
                summary_rows.append({
                    'genome': genome_name,
                    'status': 'Success',
                    'k_locus': k_info,
                    'o_locus': o_info,
                    'total_loci': len(hits),
                    'k_details': k_info if k_hits else "No K locus",
                    'o_details': o_info if o_hits else "No O locus"
                })
            else:
                summary_rows.append({
                    'genome': genome_name,
                    'status': 'Failed',
                    'k_locus': 'N/A',
                    'o_locus': 'N/A',
                    'total_loci': 0,
                    'error': result.get('error', 'Unknown error')
                })
        
        # JavaScript for interactive features (only in summary report)
        interactive_js = f"""
        <script>
            // Search functionality
            function searchTable(tableId, searchTerm) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let visibleCount = 0;
                
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const text = row.textContent.toLowerCase();
                    if (text.includes(searchTerm.toLowerCase())) {{
                        row.style.display = '';
                        visibleCount++;
                    }} else {{
                        row.style.display = 'none';
                    }}
                }}
                
                // Update result count
                const resultCounter = document.getElementById('result-counter-' + tableId);
                if (resultCounter) {{
                    resultCounter.textContent = visibleCount + ' genomes found';
                }}
            }}
            
            // Export to CSV
            function exportToCSV(tableId, filename) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let csv = [];
                
                // Add headers
                const headerCells = rows[0].getElementsByTagName('th');
                const headerRow = [];
                for (let cell of headerCells) {{
                    headerRow.push(cell.textContent);
                }}
                csv.push(headerRow.join(','));
                
                // Add data
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
                
                // Create download
                const blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            
            // Export to JSON
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
                a.href = url;
                a.download = 'Kaptive_summary_data.json';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            
            // Print report
            function printReport() {{
                window.print();
            }}
            
            // Sort table
            function sortTable(tableId, columnIndex) {{
                const table = document.getElementById(tableId);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                
                rows.sort((a, b) => {{
                    const aText = a.children[columnIndex].textContent.trim();
                    const bText = b.children[columnIndex].textContent.trim();
                    
                    // Try to convert to number if possible
                    const aNum = parseFloat(aText);
                    const bNum = parseFloat(bText);
                    
                    if (!isNaN(aNum) && !isNaN(bNum)) {{
                        return aNum - bNum;
                    }}
                    
                    return aText.localeCompare(bText);
                }});
                
                // Clear and re-add rows
                rows.forEach(row => tbody.appendChild(row));
            }}
            
            // Filter by status
            function filterByStatus(status) {{
                const table = document.getElementById('detailed-results-table');
                const rows = table.getElementsByTagName('tr');
                
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const statusCell = row.children[1].textContent.trim();
                    
                    if (status === 'all' || statusCell === status) {{
                        row.style.display = '';
                    }} else {{
                        row.style.display = 'none';
                    }}
                }}
            }}
            
            // Rotating quotes
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            
            // Initialize
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
                setInterval(rotateQuote, 10000);
            }});
        </script>
        """
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope Kaptive K/O Analysis - Summary Report</title>
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
            max-width: 1600px;
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
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            font-weight: 600;
            cursor: pointer;
            position: relative;
        }}
        
        .gene-table th:hover {{
            background: linear-gradient(135deg, #2563eb 0%, #1e3a8a 100%);
        }}
        
        .gene-table th.sortable::after {{
            content: "↕";
            position: absolute;
            right: 8px;
            opacity: 0.6;
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
        
        .k-stat-card {{
            background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        
        .o-stat-card {{
            background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%);
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
        
        /* Interactive controls */
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
        
        .filter-buttons {{
            display: flex;
            gap: 5px;
            flex-wrap: wrap;
        }}
        
        .filter-buttons button {{
            padding: 6px 12px;
            background: #e9ecef;
            color: #495057;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            cursor: pointer;
            font-size: 13px;
        }}
        
        .filter-buttons button.active {{
            background: #3b82f6;
            color: white;
            border-color: #3b82f6;
        }}
        
        .filter-buttons button:hover {{
            background: #dee2e6;
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
            font-size: 14px;
            transition: background 0.3s;
        }}
        
        .export-buttons button:hover {{
            background: #2563eb;
        }}
        
        .export-buttons button.print {{
            background: #10b981;
        }}
        
        .export-buttons button.print:hover {{
            background: #059669;
        }}
        
        .export-buttons button.json {{
            background: #8b5cf6;
        }}
        
        .export-buttons button.json:hover {{
            background: #7c3aed;
        }}
        
        .result-counter {{
            font-size: 0.9em;
            color: #666;
            font-style: italic;
            margin-left: auto;
        }}
        
        .success-row {{ border-left: 4px solid #28a745; }}
        .failed-row {{ border-left: 4px solid #dc3545; background-color: #f8d7da; }}
        
        /* Type distribution */
        .type-distribution {{
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
            margin-top: 20px;
        }}
        
        .type-list {{
            flex: 1;
            min-width: 300px;
        }}
        
        .type-list ul {{
            list-style-type: none;
            padding: 0;
        }}
        
        .type-list li {{
            padding: 8px 0;
            border-bottom: 1px solid #eee;
            display: flex;
            justify-content: space-between;
        }}
        
        .type-list .type-name {{
            font-weight: 600;
        }}
        
        .type-list .type-count {{
            color: #666;
            background: #f8f9fa;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 0.9em;
        }}
        
        @media print {{
            .interactive-controls {{ display: none !important; }}
            body {{ background: white !important; color: black !important; }}
            .card {{ background: white !important; color: black !important; box-shadow: none !important; }}
            .gene-table {{ box-shadow: none !important; }}
        }}
    </style>
    {interactive_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 AcinetoScope Kaptive K/O Analysis - Summary Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive A. baumannii Capsule (K) and Lipooligosaccharide (O) Locus Analysis Across All Genomes</p>
                <p style="color: #666; font-size: 1.1em;">Kaptive {self.metadata['kaptive_version']} | Database: {self.metadata['database_version']}</p>
                
                <div class="interactive-controls">
                    <div class="search-box">
                        <input type="text" id="search-detailed-results" 
                               placeholder="🔍 Search genomes, loci, or types..." 
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
                        <button onclick="exportToCSV('detailed-results-table', 'Kaptive_summary.csv')">📥 Export CSV</button>
                        <button onclick="exportToJSON()" class="json">📥 Export JSON</button>
                        <button onclick="printReport()" class="print">🖨️ Print</button>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""
        
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📊 Overall Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{stats['total_genomes']}</p>
                </div>
                <div class="stat-card">
                    <h3>Successful Analyses</h3>
                    <p style="font-size: 2em; margin: 0;">{stats['successful_genomes']}</p>
                    <p style="font-size: 0.9em; margin-top: 5px;">({(stats['successful_genomes']/stats['total_genomes']*100):.1f}%)</p>
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
"""
        
        # Type distribution
        if stats['common_k_types'] or stats['common_o_types']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📈 Most Common Locus Types</h2>
            <div class="type-distribution">
"""
            
            if stats['common_k_types']:
                html_content += """
                <div class="type-list">
                    <h3 style="color: #28a745;">Top K Locus Types</h3>
                    <ul>
"""
                for type_info in stats['common_k_types']:
                    percentage = (type_info['count'] / stats['genomes_with_k'] * 100) if stats['genomes_with_k'] > 0 else 0
                    html_content += f"""
                        <li>
                            <span class="type-name">{type_info['type']}</span>
                            <span class="type-count">{type_info['count']} genomes ({percentage:.1f}%)</span>
                        </li>
"""
                html_content += """
                    </ul>
                </div>
"""
            
            if stats['common_o_types']:
                html_content += """
                <div class="type-list">
                    <h3 style="color: #17a2b8;">Top O Locus Types</h3>
                    <ul>
"""
                for type_info in stats['common_o_types']:
                    percentage = (type_info['count'] / stats['genomes_with_o'] * 100) if stats['genomes_with_o'] > 0 else 0
                    html_content += f"""
                        <li>
                            <span class="type-name">{type_info['type']}</span>
                            <span class="type-count">{type_info['count']} genomes ({percentage:.1f}%)</span>
                        </li>
"""
                html_content += """
                    </ul>
                </div>
"""
            
            html_content += """
            </div>
        </div>
"""
        
        # Detailed results table
        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🔍 Detailed Results by Genome</h2>
            
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
"""
        
        for row in summary_rows:
            row_class = "success-row" if row['status'] == 'Success' else "failed-row"
            error_note = f"Error: {row.get('error', '')}" if row.get('error') else ""
            details = error_note if error_note else f"K: {row.get('k_details', 'N/A')} | O: {row.get('o_details', 'N/A')}"
            
            html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{row['genome']}</strong></td>
                        <td><span style="color: {'#28a745' if row['status'] == 'Success' else '#dc3545'}; font-weight: bold;">{row['status']}</span></td>
                        <td>{row['k_locus']}</td>
                        <td>{row['o_locus']}</td>
                        <td>{row['total_loci']}</td>
                        <td>{details}</td>
                    </tr>
"""
        
        html_content += """
                </tbody>
            </table>
        </div>
        
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📁 Generated Files</h2>
            <ul style="color: #666; font-size: 1.1em;">
                <li><strong>Kaptive_summary.tsv</strong> - Complete Kaptive data for all genomes</li>
                <li><strong>Kaptive_summary.json</strong> - JSON summary with all data</li>
                <li><strong>Individual genome HTML reports</strong> - Detailed analysis per genome (in genome folders)</li>
                <li><strong>Individual genome JSON reports</strong> - JSON data per genome (in genome folders)</li>
                <li><strong>This interactive summary report</strong> - Cross-genome analysis with search/filter</li>
            </ul>
            
            <div class="export-buttons" style="margin-top: 20px;">
                <button onclick="window.open('Kaptive_summary.tsv')">📄 View TSV Summary</button>
                <button onclick="window.open('Kaptive_summary.json')">📄 View JSON Summary</button>
            </div>
        </div>
        
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
                Analysis performed using AcinetoScope Kaptive K/O Analysis v1.0.0
                with Kaptive v3.1.0 and A. baumannii databases
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write summary HTML report
        html_file = os.path.join(output_base, "Kaptive_summary.html")
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        self.logger.info(f"✓ A. baumannii Kaptive summary HTML report created: {html_file}")


def main():
    """Command line interface for A. baumannii K/O analysis"""
    parser = argparse.ArgumentParser(
        description='AcinetoScope Kaptive K/O Analysis - A. baumannii Capsule and Lipooligosaccharide Typing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all A. baumannii FASTA files
  python acineto_kaptive.py "*.fna"
  
  # Run on specific pattern
  python acineto_kaptive.py "AB_*.fasta"
  
  # Run with custom output directory
  python acineto_kaptive.py "*.fa" --output my_kaptive_results

Supported FASTA extensions: .fasta, .fa, .fna, .faa
        """
    )
    
    parser.add_argument('pattern', help='File pattern for A. baumannii genomes (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--output', '-o', default='kaptive_results', 
                       help='Output directory (default: kaptive_results)')
    
    args = parser.parse_args()
    
    # Print AcinetoScope banner
    print("\n" + "="*80)
    print(r"""
 █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
""")
    print("A. baumannii Kaptive K/O Locus Analysis")
    print("="*80)
    print(f"Author: Brown Beckley | Email: brownbeckley94@gmail.com")
    print(f"Affiliation: University of Ghana Medical School - Department of Medical Biochemistry")
    print("="*80)
    
    executor = AcinetoKaptive()
    
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        
        # Print summary
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🧬 AcinetoScope Kaptive K/O Analysis FINAL SUMMARY")
        executor.logger.info("="*50)
        
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
                
                executor.logger.info(f"✓ {genome_name}: {result['hit_count']} K/O loci")
        
        executor.logger.info("\n📊 A. BAUMANNII K/O SUMMARY STATISTICS:")
        executor.logger.info(f"   Total genomes processed: {len(results)}")
        executor.logger.info(f"   Successfully analyzed: {successful_genomes}")
        executor.logger.info(f"   Total K/O loci found: {total_hits}")
        executor.logger.info(f"   K loci (capsule): {k_hits}")
        executor.logger.info(f"   O loci (lipooligosaccharide): {o_hits}")
        
        # Show summary file locations
        executor.logger.info("\n📁 SUMMARY FILES CREATED:")
        executor.logger.info(f"   Comprehensive Kaptive data: {args.output}/Kaptive_summary.tsv")
        executor.logger.info(f"   JSON summary: {args.output}/Kaptive_summary.json")
        executor.logger.info(f"   Summary HTML report: {args.output}/Kaptive_summary.html")
        executor.logger.info(f"   Individual genome reports in: {args.output}/*/")
        
        import random
        executor.logger.info(f"\n💡 {random.choice(executor.science_quotes)}")
        
    except Exception as e:
        executor.logger.error(f"A. baumannii K/O analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()