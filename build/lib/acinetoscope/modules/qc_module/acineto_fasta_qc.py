#!/usr/bin/env python3
"""
AcinetoScope FASTA QC - Comprehensive Quality Control with Beautiful HTML Reports
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2025-12-28
"""

import os
import sys
import glob
import json
import math
import statistics
from pathlib import Path
from datetime import datetime
from collections import Counter, defaultdict
import argparse
import logging
import subprocess
from typing import List, Dict, Any, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

# BioPython imports
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

class AcinetoFASTAQC:
    """Comprehensive FASTA QC with beautiful HTML"""
    
    def __init__(self, cpus: int = None):
        # Setup logging
        self.logger = self._setup_logging()
        
        # ASCII Art for AcinetoScope
        self.ascii_art = r"""
 █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
"""
        
        # Metadata
        self.metadata = {
            "tool_name": "AcinetoScope FASTA QC Analysis",
            "version": "1.0.0", 
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "biopython_version": self._get_biopython_version()
        }
        
        # Colors matching
        self.colors = {
            'primary': '#3b82f6',
            'success': '#28a745',
            'warning': '#ffc107',
            'danger': '#dc3545',
            'info': '#17a2b8',
            'dark': '#1f2937',
            'light': '#f8f9fa'
        }
        
        # Science quotes
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
        
        # Thresholds
        self.thresholds = {
            'gc_normal': (35, 65),
            'short_seq': 100,
            'long_seq': 1000000,
            'max_n_run': 100,
            'max_homopolymer': 20,
            'ambiguous_critical': 5.0,
            'ambiguous_warning': 1.0,
        }
        
        self.cpus = cpus or os.cpu_count() or 1
    
    def _setup_logging(self):
        """Setup logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def _get_biopython_version(self) -> str:
        """Get BioPython version"""
        try:
            import Bio
            return Bio.__version__
        except:
            return "Unknown"
    
    def analyze_file(self, fasta_file: str) -> Dict[str, Any]:
        """Comprehensive analysis of a FASTA file"""
        try:
            self.logger.info(f"Analyzing {os.path.basename(fasta_file)}...")
            
            # Read sequences
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            if not sequences:
                return {
                    'filename': os.path.basename(fasta_file),
                    'filepath': fasta_file,
                    'status': 'error',
                    'error': 'No sequences found'
                }
            
            # Basic stats
            seq_lengths = [len(seq) for seq in sequences]
            total_length = sum(seq_lengths)
            
            # Sort for N statistics
            sorted_lengths = sorted(seq_lengths, reverse=True)
            
            # N statistics
            n50 = self._calculate_nx(sorted_lengths, total_length, 50)
            n75 = self._calculate_nx(sorted_lengths, total_length, 75)
            n90 = self._calculate_nx(sorted_lengths, total_length, 90)
            
            # L statistics
            l50 = self._calculate_lx(sorted_lengths, total_length, 50)
            l75 = self._calculate_lx(sorted_lengths, total_length, 75)
            l90 = self._calculate_lx(sorted_lengths, total_length, 90)
            
            # Nucleotide composition
            base_counts = Counter()
            gc_contents = []
            
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                base_counts.update(seq_str)
                gc_contents.append(gc_fraction(seq_str) * 100)
            
            total_bases = sum(base_counts.values())
            
            # Calculate percentages
            a_count = base_counts.get('A', 0)
            t_count = base_counts.get('T', 0)
            g_count = base_counts.get('G', 0)
            c_count = base_counts.get('C', 0)
            
            a_percent = (a_count / total_bases) * 100 if total_bases > 0 else 0
            t_percent = (t_count / total_bases) * 100 if total_bases > 0 else 0
            g_percent = (g_count / total_bases) * 100 if total_bases > 0 else 0
            c_percent = (c_count / total_bases) * 100 if total_bases > 0 else 0
            gc_percent = g_percent + c_percent
            at_percent = a_percent + t_percent
            
            # Ambiguous bases
            ambiguous_bases = sum(base_counts.get(b, 0) for b in ['N', 'Y', 'R', 'W', 'S', 'K', 'M', 'B', 'D', 'H', 'V'])
            ambiguous_percent = (ambiguous_bases / total_bases) * 100 if total_bases > 0 else 0
            
            # N statistics
            sequences_with_n = sum(1 for seq in sequences if 'N' in str(seq.seq).upper())
            total_n_bases = base_counts.get('N', 0)
            
            # Find longest N-run
            max_n_run = 0
            n_runs = []
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                import re
                for match in re.finditer(r'N+', seq_str):
                    run_len = len(match.group())
                    n_runs.append(run_len)
                    if run_len > max_n_run:
                        max_n_run = run_len
            
            # Homopolymers
            homopolymers = []
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                for base in ['A', 'T', 'G', 'C']:
                    for match in re.finditer(f'{base}+', seq_str):
                        if len(match.group()) > 4:
                            homopolymers.append({
                                'base': base,
                                'length': len(match.group()),
                                'position': match.start()
                            })
            
            max_homopolymer = max([h['length'] for h in homopolymers]) if homopolymers else 0
            
            # Duplicate sequences
            seq_hashes = set()
            duplicate_sequences = 0
            for seq in sequences:
                seq_hash = hash(str(seq.seq))
                if seq_hash in seq_hashes:
                    duplicate_sequences += 1
                else:
                    seq_hashes.add(seq_hash)
            
            # Short and long sequences
            short_sequences = sum(1 for length in seq_lengths if length < self.thresholds['short_seq'])
            long_sequences = sum(1 for length in seq_lengths if length > self.thresholds['long_seq'])
            
            # Length distribution
            length_distribution = self._create_length_bins(seq_lengths)
            
            # Compile results
            results = {
                'filename': os.path.basename(fasta_file),
                'filepath': fasta_file,
                'file_size_mb': os.path.getsize(fasta_file) / (1024 ** 2),
                'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'status': 'success',
                'total_sequences': len(sequences),
                'total_length': total_length,
                'total_bases': total_bases,
                'longest_sequence': max(seq_lengths),
                'shortest_sequence': min(seq_lengths),
                'mean_length': statistics.mean(seq_lengths) if seq_lengths else 0,
                'median_length': statistics.median(seq_lengths) if seq_lengths else 0,
                'n50': n50,
                'n75': n75,
                'n90': n90,
                'l50': l50,
                'l75': l75,
                'l90': l90,
                'gc_percent': gc_percent,
                'at_percent': at_percent,
                'a_percent': a_percent,
                't_percent': t_percent,
                'g_percent': g_percent,
                'c_percent': c_percent,
                'ambiguous_percent': ambiguous_percent,
                'sequences_with_n': sequences_with_n,
                'total_n_bases': total_n_bases,
                'max_n_run': max_n_run,
                'homopolymers_count': len(homopolymers),
                'max_homopolymer': max_homopolymer,
                'duplicate_sequences': duplicate_sequences,
                'short_sequences': short_sequences,
                'long_sequences': long_sequences,
                'gc_distribution': {
                    'mean': statistics.mean(gc_contents) if gc_contents else 0,
                    'median': statistics.median(gc_contents) if gc_contents else 0,
                    'min': min(gc_contents) if gc_contents else 0,
                    'max': max(gc_contents) if gc_contents else 0,
                },
                'length_distribution': length_distribution,
                'base_counts': dict(base_counts),
                'n_runs': n_runs,
                'warnings': self._generate_warnings(
                    gc_percent, ambiguous_percent, max_n_run, max_homopolymer,
                    short_sequences, long_sequences, duplicate_sequences, len(sequences)
                )
            }
            
            self.logger.info(f"✓ {os.path.basename(fasta_file)}: {len(sequences)} sequences, {total_length:,} bp, N50: {n50:,}, GC: {gc_percent:.1f}%")
            return results
            
        except Exception as e:
            self.logger.error(f"✗ Error analyzing {fasta_file}: {e}")
            return {
                'filename': os.path.basename(fasta_file),
                'status': 'error',
                'error': str(e)
            }
    
    def _calculate_nx(self, sorted_lengths: List[int], total_length: int, x: int) -> int:
        """Calculate Nx (N50, N75, N90)"""
        if not sorted_lengths:
            return 0
        
        target = total_length * (x / 100)
        cumulative = 0
        
        for length in sorted_lengths:
            cumulative += length
            if cumulative >= target:
                return length
        
        return sorted_lengths[-1]
    
    def _calculate_lx(self, sorted_lengths: List[int], total_length: int, x: int) -> int:
        """Calculate Lx (L50, L75, L90)"""
        if not sorted_lengths:
            return 0
        
        target = total_length * (x / 100)
        cumulative = 0
        
        for i, length in enumerate(sorted_lengths, 1):
            cumulative += length
            if cumulative >= target:
                return i
        
        return len(sorted_lengths)
    
    def _create_length_bins(self, lengths: List[int]) -> Dict[str, int]:
        """Create length distribution bins"""
        bins = {
            '< 100 bp': 0,
            '100-500 bp': 0,
            '500-1k bp': 0,
            '1k-5k bp': 0,
            '5k-10k bp': 0,
            '10k-50k bp': 0,
            '50k-100k bp': 0,
            '100k-500k bp': 0,
            '500k-1M bp': 0,
            '> 1M bp': 0
        }
        
        for length in lengths:
            if length < 100:
                bins['< 100 bp'] += 1
            elif length < 500:
                bins['100-500 bp'] += 1
            elif length < 1000:
                bins['500-1k bp'] += 1
            elif length < 5000:
                bins['1k-5k bp'] += 1
            elif length < 10000:
                bins['5k-10k bp'] += 1
            elif length < 50000:
                bins['10k-50k bp'] += 1
            elif length < 100000:
                bins['50k-100k bp'] += 1
            elif length < 500000:
                bins['100k-500k bp'] += 1
            elif length < 1000000:
                bins['500k-1M bp'] += 1
            else:
                bins['> 1M bp'] += 1
        
        return bins
    
    def _generate_warnings(self, gc_percent: float, ambiguous_percent: float, 
                          max_n_run: int, max_homopolymer: int,
                          short_sequences: int, long_sequences: int, 
                          duplicate_sequences: int, total_sequences: int) -> List[Dict]:
        """Generate warning messages"""
        warnings = []
        
        # GC content
        low, high = self.thresholds['gc_normal']
        if gc_percent < low:
            warnings.append({
                'level': 'warning',
                'message': f'Low GC content ({gc_percent:.1f}%) - below normal range ({low}-{high}%)'
            })
        elif gc_percent > high:
            warnings.append({
                'level': 'warning',
                'message': f'High GC content ({gc_percent:.1f}%) - above normal range ({low}-{high}%)'
            })
        
        # Ambiguous bases
        if ambiguous_percent > self.thresholds['ambiguous_critical']:
            warnings.append({
                'level': 'danger',
                'message': f'High ambiguous bases ({ambiguous_percent:.2f}%) - may indicate poor quality'
            })
        elif ambiguous_percent > self.thresholds['ambiguous_warning']:
            warnings.append({
                'level': 'warning',
                'message': f'Elevated ambiguous bases ({ambiguous_percent:.2f}%)'
            })
        
        # N-runs
        if max_n_run > 100:
            warnings.append({
                'level': 'danger',
                'message': f'Very long N-run detected ({max_n_run} bases) - may indicate assembly gaps'
            })
        elif max_n_run > 10:
            warnings.append({
                'level': 'warning',
                'message': f'Long N-run detected ({max_n_run} bases)'
            })
        
        # Homopolymers
        if max_homopolymer > 20:
            warnings.append({
                'level': 'danger',
                'message': f'Very long homopolymer ({max_homopolymer} bases) - may cause sequencing errors'
            })
        elif max_homopolymer > 10:
            warnings.append({
                'level': 'warning',
                'message': f'Long homopolymer ({max_homopolymer} bases)'
            })
        
        # Short sequences
        if short_sequences > total_sequences * 0.5:
            warnings.append({
                'level': 'danger',
                'message': f'Many short sequences ({short_sequences}) - may indicate poor assembly'
            })
        elif short_sequences > total_sequences * 0.1:
            warnings.append({
                'level': 'warning',
                'message': f'Some short sequences ({short_sequences})'
            })
        
        # Long sequences
        if long_sequences > 0:
            warnings.append({
                'level': 'warning',
                'message': f'Very long sequences detected ({long_sequences}) - may indicate contamination'
            })
        
        # Duplicate sequences
        duplicate_percent = (duplicate_sequences / total_sequences) * 100 if total_sequences > 0 else 0
        if duplicate_percent > 10:
            warnings.append({
                'level': 'warning',
                'message': f'Duplicate sequences detected ({duplicate_sequences}, {duplicate_percent:.1f}%)'
            })
        
        return warnings
    
    def create_individual_html_report(self, results: Dict[str, Any], output_dir: str) -> str:
        """Create comprehensive HTML report for a single FASTA file"""
        # Create individual folder for this FASTA file
        filename_no_ext = Path(results['filename']).stem
        sample_dir = os.path.join(output_dir, filename_no_ext)
        os.makedirs(sample_dir, exist_ok=True)
        
        html_file = os.path.join(sample_dir, f"{filename_no_ext}_fasta_qc_report.html")
        
        # JavaScript for interactive features
        js_content = f"""
        <script>
            // Print report
            function printReport() {{
                window.print();
            }}
            
            // Export to JSON
            function exportToJSON() {{
                const reportData = {json.dumps(results, indent=2)};
                const dataStr = JSON.stringify(reportData, null, 2);
                const dataBlob = new Blob([dataStr], {{ type: 'application/json' }});
                const url = URL.createObjectURL(dataBlob);
                const a = document.createElement('a');
                a.href = url;
                a.download = '{filename_no_ext}_fasta_qc_report.json';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                URL.revokeObjectURL(url);
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
                
                // Add tooltips
                const tooltips = document.querySelectorAll('[data-toggle="tooltip"]');
                tooltips.forEach(el => {{
                    el.addEventListener('mouseenter', function(e) {{
                        const tooltip = document.createElement('div');
                        tooltip.className = 'tooltip-custom';
                        tooltip.textContent = this.getAttribute('title');
                        tooltip.style.position = 'absolute';
                        tooltip.style.background = 'rgba(0,0,0,0.8)';
                        tooltip.style.color = 'white';
                        tooltip.style.padding = '5px 10px';
                        tooltip.style.borderRadius = '4px';
                        tooltip.style.zIndex = '1000';
                        tooltip.style.left = e.pageX + 'px';
                        tooltip.style.top = (e.pageY - 30) + 'px';
                        document.body.appendChild(tooltip);
                        this._tooltip = tooltip;
                    }});
                    
                    el.addEventListener('mouseleave', function() {{
                        if (this._tooltip) {{
                            document.body.removeChild(this._tooltip);
                            this._tooltip = null;
                        }}
                    }});
                }});
            }});
            
            // Sort table
            function sortTable(tableId, columnIndex) {{
                const table = document.getElementById(tableId);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                
                rows.sort((a, b) => {{
                    const aText = a.children[columnIndex].textContent.trim();
                    const bText = b.children[columnIndex].textContent.trim();
                    
                    // Try to convert to number if possible
                    const aNum = parseFloat(aText.replace(/,/g, '').replace('%', ''));
                    const bNum = parseFloat(bText.replace(/,/g, '').replace('%', ''));
                    
                    if (!isNaN(aNum) && !isNaN(bNum)) {{
                        return aNum - bNum;
                    }}
                    
                    return aText.localeCompare(bText);
                }});
                
                // Clear and re-add rows
                rows.forEach(row => tbody.appendChild(row));
            }}
        </script>
        """
        
        # HTML Content
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope FASTA QC Report - {results['filename']}</title>
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
        
        .qc-table {{ 
            width: 100%; 
            border-collapse: collapse; 
            margin: 20px 0; 
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            font-size: 13px;
        }}
        
        .qc-table th, .qc-table td {{ 
            padding: 12px 8px; 
            text-align: left; 
            border-bottom: 1px solid #e0e0e0; 
            word-wrap: break-word;
        }}
        
        .qc-table th {{ 
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            font-weight: 600;
            cursor: pointer;
            position: relative;
        }}
        
        .qc-table th:hover {{ 
            background: linear-gradient(135deg, #2563eb 0%, #1e3a8a 100%);
        }}
        
        .qc-table th.sortable::after {{
            content: "↕";
            position: absolute;
            right: 8px;
            opacity: 0.6;
        }}
        
        tr:hover {{ background-color: #f8f9fa; }}
        .success {{ color: #28a745; font-weight: 600; }}
        .warning {{ color: #ffc107; font-weight: 600; }}
        .danger {{ color: #dc3545; font-weight: 600; }}
        
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
            min-width: 180px;
        }}
        
        .gc-stat-card {{
            background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
        }}
        
        .at-stat-card {{
            background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%);
            color: white;
            padding: 20px; 
            border-radius: 12px; 
            text-align: center; 
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
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
        
        /* Warning boxes */
        .warning-box {{
            background: #fff3cd;
            border: 1px solid #ffc107;
            border-radius: 5px;
            padding: 15px;
            margin: 10px 0;
            color: #856404;
        }}
        
        .danger-box {{
            background: #f8d7da;
            border: 1px solid #dc3545;
            border-radius: 5px;
            padding: 15px;
            margin: 10px 0;
            color: #721c24;
        }}
        
        /* Responsive table */
        .table-container {{
            overflow-x: auto;
            margin: 20px 0;
        }}
        
        /* Nucleotide composition */
        .composition-grid {{
            display: flex;
            justify-content: space-around;
            flex-wrap: wrap;
            margin: 20px 0;
        }}
        
        .composition-item {{
            text-align: center;
            margin: 10px;
            flex: 1;
            min-width: 120px;
        }}
        
        .composition-value {{
            font-size: 24px;
            font-weight: bold;
        }}
        
        .composition-label {{
            font-size: 12px;
            color: #666;
            margin-top: 5px;
        }}
        
        @media (max-width: 1200px) {{
            .qc-table {{
                font-size: 11px;
            }}
            
            .qc-table th, .qc-table td {{
                padding: 8px 6px;
            }}
        }}
        
        @media print {{
            .controls {{ display: none !important; }}
            body {{ background: white !important; color: black !important; }}
            .card {{ background: white !important; color: black !important; box-shadow: none !important; }}
            .qc-table {{ box-shadow: none !important; }}
        }}
    </style>
    {js_content}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 AcinetoScope FASTA QC Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive Quality Control for FASTA Files</p>
                <p style="color: #666; font-size: 1.1em;"><strong>Sample:</strong> {results['filename']}</p>
                
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
        
        # Summary Statistics
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📊 FASTA QC Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Sequences</h3>
                    <p style="font-size: 2em; margin: 0;">{results['total_sequences']:,}</p>
                </div>
                <div class="stat-card">
                    <h3>Total Length</h3>
                    <p style="font-size: 2em; margin: 0;">{results['total_length']:,}</p>
                    <p style="font-size: 0.9em; margin-top: 5px;">bp</p>
                </div>
                <div class="stat-card">
                    <h3>N50</h3>
                    <p style="font-size: 2em; margin: 0;">{results['n50']:,}</p>
                    <p style="font-size: 0.9em; margin-top: 5px;">bp</p>
                </div>
                <div class="gc-stat-card">
                    <h3>GC Content</h3>
                    <p style="font-size: 2em; margin: 0;">{results['gc_percent']:.1f}%</p>
                    <p style="font-size: 0.9em; margin-top: 5px;">Normal: {self.thresholds['gc_normal'][0]}-{self.thresholds['gc_normal'][1]}%</p>
                </div>
                <div class="at-stat-card">
                    <h3>AT Content</h3>
                    <p style="font-size: 2em; margin: 0;">{results['at_percent']:.1f}%</p>
                </div>
            </div>
            <p><strong>File:</strong> {results['filename']}</p>
            <p><strong>File Size:</strong> {results['file_size_mb']:.2f} MB</p>
            <p><strong>Analysis Date:</strong> {results['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
        </div>
"""
        
        # Basic Statistics Table
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📈 Basic Statistics</h2>
            
            <div class="table-container">
                <table class="qc-table" id="basic-stats-table">
                    <thead>
                        <tr>
                            <th class="sortable" onclick="sortTable('basic-stats-table', 0)">Metric</th>
                            <th class="sortable" onclick="sortTable('basic-stats-table', 1)">Value</th>
                            <th>Description</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td><strong>Total Sequences</strong></td>
                            <td>{results['total_sequences']:,}</td>
                            <td>Number of sequences in the file</td>
                        </tr>
                        <tr>
                            <td><strong>Total Length</strong></td>
                            <td>{results['total_length']:,} bp</td>
                            <td>Total number of bases (including Ns)</td>
                        </tr>
                        <tr>
                            <td><strong>Total Bases</strong></td>
                            <td>{results['total_bases']:,} bp</td>
                            <td>Total bases excluding ambiguous characters</td>
                        </tr>
                        <tr>
                            <td><strong>Longest Sequence</strong></td>
                            <td>{results['longest_sequence']:,} bp</td>
                            <td>Length of the longest sequence</td>
                        </tr>
                        <tr>
                            <td><strong>Shortest Sequence</strong></td>
                            <td>{results['shortest_sequence']:,} bp</td>
                            <td>Length of the shortest sequence</td>
                        </tr>
                        <tr>
                            <td><strong>Mean Length</strong></td>
                            <td>{results['mean_length']:,.0f} bp</td>
                            <td>Average sequence length</td>
                        </tr>
                        <tr>
                            <td><strong>Median Length</strong></td>
                            <td>{results['median_length']:,} bp</td>
                            <td>Median sequence length</td>
                        </tr>
                        <tr>
                            <td><strong>N50</strong></td>
                            <td>{results['n50']:,} bp</td>
                            <td>Length for which 50% of total bases are in longer sequences</td>
                        </tr>
                        <tr>
                            <td><strong>N75</strong></td>
                            <td>{results['n75']:,} bp</td>
                            <td>Length for which 75% of total bases are in longer sequences</td>
                        </tr>
                        <tr>
                            <td><strong>N90</strong></td>
                            <td>{results['n90']:,} bp</td>
                            <td>Length for which 90% of total bases are in longer sequences</td>
                        </tr>
                        <tr>
                            <td><strong>L50</strong></td>
                            <td>{results['l50']:,}</td>
                            <td>Number of sequences that make up 50% of total length</td>
                        </tr>
                        <tr>
                            <td><strong>L75</strong></td>
                            <td>{results['l75']:,}</td>
                            <td>Number of sequences that make up 75% of total length</td>
                        </tr>
                        <tr>
                            <td><strong>L90</strong></td>
                            <td>{results['l90']:,}</td>
                            <td>Number of sequences that make up 90% of total length</td>
                        </tr>
                        <tr>
                            <td><strong>GC Content</strong></td>
                            <td>{results['gc_percent']:.1f}%</td>
                            <td>Percentage of G and C bases</td>
                        </tr>
                        <tr>
                            <td><strong>AT Content</strong></td>
                            <td>{results['at_percent']:.1f}%</td>
                            <td>Percentage of A and T bases</td>
                        </tr>
                        <tr>
                            <td><strong>Ambiguous Bases</strong></td>
                            <td>{results['ambiguous_percent']:.2f}%</td>
                            <td>Percentage of N and other ambiguous bases</td>
                        </tr>
                        <tr>
                            <td><strong>Sequences with Ns</strong></td>
                            <td>{results['sequences_with_n']:,}</td>
                            <td>Number of sequences containing N bases</td>
                        </tr>
                        <tr>
                            <td><strong>Total N Bases</strong></td>
                            <td>{results['total_n_bases']:,}</td>
                            <td>Total count of N bases</td>
                        </tr>
                        <tr>
                            <td><strong>Max N-run</strong></td>
                            <td>{results['max_n_run']:,}</td>
                            <td>Longest consecutive run of N bases</td>
                        </tr>
                        <tr>
                            <td><strong>Homopolymers</strong></td>
                            <td>{results['homopolymers_count']:,}</td>
                            <td>Number of homopolymers longer than 4 bases</td>
                        </tr>
                        <tr>
                            <td><strong>Max Homopolymer</strong></td>
                            <td>{results['max_homopolymer']:,}</td>
                            <td>Longest homopolymer run</td>
                        </tr>
                        <tr>
                            <td><strong>Duplicate Sequences</strong></td>
                            <td>{results['duplicate_sequences']:,}</td>
                            <td>Number of identical sequences</td>
                        </tr>
                        <tr>
                            <td><strong>Short Sequences (&lt;{self.thresholds['short_seq']} bp)</strong></td>
                            <td>{results['short_sequences']:,}</td>
                            <td>Sequences shorter than {self.thresholds['short_seq']} bp</td>
                        </tr>
                        <tr>
                            <td><strong>Long Sequences (&gt;{self.thresholds['long_seq']:,} bp)</strong></td>
                            <td>{results['long_sequences']:,}</td>
                            <td>Sequences longer than {self.thresholds['long_seq']:,} bp</td>
                        </tr>
                    </tbody>
                </table>
            </div>
        </div>
"""
        
        # Nucleotide Composition
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">🧬 Nucleotide Composition</h2>
            
            <div class="composition-grid">
                <div class="composition-item">
                    <div class="composition-value" style="color: #28a745;">{results['gc_percent']:.1f}%</div>
                    <div class="composition-label">GC Content</div>
                </div>
                <div class="composition-item">
                    <div class="composition-value" style="color: #17a2b8;">{results['at_percent']:.1f}%</div>
                    <div class="composition-label">AT Content</div>
                </div>
                <div class="composition-item">
                    <div class="composition-value">{results['a_percent']:.1f}%</div>
                    <div class="composition-label">Adenine (A)</div>
                </div>
                <div class="composition-item">
                    <div class="composition-value">{results['t_percent']:.1f}%</div>
                    <div class="composition-label">Thymine (T)</div>
                </div>
                <div class="composition-item">
                    <div class="composition-value">{results['g_percent']:.1f}%</div>
                    <div class="composition-label">Guanine (G)</div>
                </div>
                <div class="composition-item">
                    <div class="composition-value">{results['c_percent']:.1f}%</div>
                    <div class="composition-label">Cytosine (C)</div>
                </div>
                <div class="composition-item">
                    <div class="composition-value" style="color: #dc3545;">{results['ambiguous_percent']:.2f}%</div>
                    <div class="composition-label">Ambiguous Bases</div>
                </div>
            </div>
        </div>
"""
        
        # Length Distribution
        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📊 Length Distribution</h2>
            
            <div class="table-container">
                <table class="qc-table" id="length-dist-table">
                    <thead>
                        <tr>
                            <th class="sortable" onclick="sortTable('length-dist-table', 0)">Length Range</th>
                            <th class="sortable" onclick="sortTable('length-dist-table', 1)">Count</th>
                            <th class="sortable" onclick="sortTable('length-dist-table', 2)">Percentage</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        
        total_seqs = results['total_sequences']
        for length_range, count in results['length_distribution'].items():
            percentage = (count / total_seqs * 100) if total_seqs > 0 else 0
            html_content += f"""
                        <tr>
                            <td>{length_range}</td>
                            <td>{count:,}</td>
                            <td>{percentage:.1f}%</td>
                        </tr>
"""
        
        html_content += """
                    </tbody>
                </table>
            </div>
        </div>
"""
        
        # Warnings section
        if results.get('warnings'):
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">⚠️ Quality Warnings</h2>
"""
            
            for warning in results['warnings']:
                box_class = 'danger-box' if warning['level'] == 'danger' else 'warning-box'
                html_content += f"""
            <div class="{box_class}">
                <strong>{warning['level'].upper()}:</strong> {warning['message']}
            </div>
"""
            
            html_content += """
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
                Analysis performed using AcinetoScope FASTA QC v{self.metadata['version']}
                with BioPython {self.metadata['biopython_version']}
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write HTML file
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        # Create JSON report
        json_file = os.path.join(sample_dir, f"{filename_no_ext}_fasta_qc_report.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, default=str)
        
        self.logger.info(f"✓ HTML report generated: {html_file}")
        self.logger.info(f"✓ JSON report generated: {json_file}")
        
        return html_file
    
    def create_summary_report(self, all_results: List[Dict[str, Any]], output_dir: str):
        """Create summary report for multiple FASTA files"""
        successful_results = [r for r in all_results if r.get('status') == 'success']
        
        if not successful_results:
            self.logger.warning("No successful analyses to create summary report")
            return
        
        # Create TSV summary
        tsv_file = os.path.join(output_dir, "FASTA_QC_summary.tsv")
        with open(tsv_file, 'w', encoding='utf-8') as f:
            # Header
            headers = [
                'Filename', 'Total Sequences', 'Total Length', 'Total Bases',
                'GC Content (%)', 'AT Content (%)', 'N50', 'N75', 'N90',
                'Median Length', 'Mean Length', 'Longest Sequence', 'Shortest Sequence',
                'Ambiguous Bases (%)', 'Sequences with Ns', 'Max N-run',
                'Homopolymers', 'Max Homopolymer', 'Duplicate Sequences',
                'Short Sequences (<100 bp)', 'Long Sequences (>1M bp)',
                'File Size (MB)', 'Warnings'
            ]
            f.write('\t'.join(headers) + '\n')
            
            for result in successful_results:
                row = [
                    result['filename'],
                    str(result['total_sequences']),
                    str(result['total_length']),
                    str(result['total_bases']),
                    f"{result['gc_percent']:.2f}",
                    f"{result['at_percent']:.2f}",
                    str(result['n50']),
                    str(result['n75']),
                    str(result['n90']),
                    str(result['median_length']),
                    f"{result['mean_length']:.0f}",
                    str(result['longest_sequence']),
                    str(result['shortest_sequence']),
                    f"{result['ambiguous_percent']:.2f}",
                    str(result['sequences_with_n']),
                    str(result['max_n_run']),
                    str(result['homopolymers_count']),
                    str(result['max_homopolymer']),
                    str(result['duplicate_sequences']),
                    str(result['short_sequences']),
                    str(result['long_sequences']),
                    f"{result['file_size_mb']:.2f}",
                    str(len(result.get('warnings', [])))
                ]
                f.write('\t'.join(row) + '\n')
        
        # Create JSON summary
        json_summary = self._create_json_summary(successful_results)
        json_file = os.path.join(output_dir, "FASTA_QC_summary.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(json_summary, f, indent=2, default=str)
        
        # Create HTML summary report
        html_file = os.path.join(output_dir, "FASTA_QC_summary.html")
        self._create_summary_html_report(successful_results, html_file, json_summary)
        
        self.logger.info(f"✓ TSV summary created: {tsv_file}")
        self.logger.info(f"✓ JSON summary created: {json_file}")
        self.logger.info(f"✓ HTML summary created: {html_file}")
    
    def _create_json_summary(self, successful_results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Create JSON summary structure"""
        summary_data = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'biopython_version': self.metadata['biopython_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_files': len(successful_results),
                'total_sequences': sum(r['total_sequences'] for r in successful_results),
                'total_length': sum(r['total_length'] for r in successful_results)
            },
            'statistics': {
                'gc_content_range': {
                    'min': min(r['gc_percent'] for r in successful_results) if successful_results else 0,
                    'max': max(r['gc_percent'] for r in successful_results) if successful_results else 0,
                    'mean': statistics.mean(r['gc_percent'] for r in successful_results) if successful_results else 0,
                    'median': statistics.median(r['gc_percent'] for r in successful_results) if successful_results else 0,
                },
                'n50_range': {
                    'min': min(r['n50'] for r in successful_results) if successful_results else 0,
                    'max': max(r['n50'] for r in successful_results) if successful_results else 0,
                    'mean': statistics.mean(r['n50'] for r in successful_results) if successful_results else 0,
                    'median': statistics.median(r['n50'] for r in successful_results) if successful_results else 0,
                },
                'files_with_warnings': sum(1 for r in successful_results if r.get('warnings')),
                'total_warnings': sum(len(r.get('warnings', [])) for r in successful_results),
                'files_with_high_gc': sum(1 for r in successful_results if r['gc_percent'] > self.thresholds['gc_normal'][1]),
                'files_with_low_gc': sum(1 for r in successful_results if r['gc_percent'] < self.thresholds['gc_normal'][0]),
            },
            'files': successful_results
        }
        
        return summary_data
    
    def _create_summary_html_report(self, successful_results: List[Dict[str, Any]], 
                                  html_file: str, json_summary: Dict[str, Any]):
        """Create HTML summary report with scrollable table"""
        
        # JavaScript for interactive features
        js_content = f"""
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
                    resultCounter.textContent = visibleCount + ' files found';
                }}
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
                    const aNum = parseFloat(aText.replace(/,/g, '').replace('%', ''));
                    const bNum = parseFloat(bText.replace(/,/g, '').replace('%', ''));
                    
                    if (!isNaN(aNum) && !isNaN(bNum)) {{
                        return aNum - bNum;
                    }}
                    
                    return aText.localeCompare(bText);
                }});
                
                // Clear and re-add rows
                rows.forEach(row => tbody.appendChild(row));
            }}
            
            // Export to CSV
            function exportToCSV() {{
                const rows = document.querySelectorAll('#qc-summary-table tr');
                let csv = [];
                
                // Add headers
                const headerCells = rows[0].querySelectorAll('th');
                const headerRow = [];
                for (let cell of headerCells) {{
                    headerRow.push(cell.textContent);
                }}
                csv.push(headerRow.join(','));
                
                // Add data
                for (let i = 1; i < rows.length; i++) {{
                    if (rows[i].style.display !== 'none') {{
                        const cells = rows[i].querySelectorAll('td');
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
                a.download = 'FASTA_QC_summary.csv';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            
            // Export to JSON
            function exportToJSON() {{
                const data = {json.dumps(json_summary, indent=2)};
                const jsonStr = JSON.stringify(data, null, 2);
                const blob = new Blob([jsonStr], {{ type: 'application/json' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = 'FASTA_QC_summary.json';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            
            // Print report
            function printReport() {{
                window.print();
            }}
            
            // Filter by warnings
            function filterByWarnings(level) {{
                const table = document.getElementById('qc-summary-table');
                const rows = table.getElementsByTagName('tr');
                
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const warningsCell = row.children[row.children.length - 1];
                    const warnings = parseInt(warningsCell.textContent);
                    
                    if (level === 'all') {{
                        row.style.display = '';
                    }} else if (level === 'with_warnings' && warnings > 0) {{
                        row.style.display = '';
                    }} else if (level === 'no_warnings' && warnings === 0) {{
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
        
        # HTML Content
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AcinetoScope FASTA QC - Summary Report</title>
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
        
        .qc-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
            font-size: 12px;
        }}
        
        .qc-table th, .qc-table td {{
            padding: 10px 8px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
            white-space: nowrap;
        }}
        
        .qc-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            font-weight: 600;
            cursor: pointer;
            position: sticky;
            top: 0;
            z-index: 10;
        }}
        
        .qc-table th:hover {{
            background: linear-gradient(135deg, #2563eb 0%, #1e3a8a 100%);
        }}
        
        .qc-table th.sortable::after {{
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
            min-width: 180px;
        }}
        
        .seq-stat-card {{
            background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
        }}
        
        .warning-stat-card {{
            background: linear-gradient(135deg, #ffc107 0%, #e0a800 100%);
            color: black;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
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
        
        /* Table container with horizontal scroll */
        .table-container {{
            overflow-x: auto;
            margin: 20px 0;
            max-height: 600px;
            overflow-y: auto;
        }}
        
        .table-container::-webkit-scrollbar {{
            width: 8px;
            height: 8px;
        }}
        
        .table-container::-webkit-scrollbar-track {{
            background: #f1f1f1;
            border-radius: 4px;
        }}
        
        .table-container::-webkit-scrollbar-thumb {{
            background: #888;
            border-radius: 4px;
        }}
        
        .table-container::-webkit-scrollbar-thumb:hover {{
            background: #555;
        }}
        
        .warning-count {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 11px;
            font-weight: bold;
        }}
        
        .warning-count.none {{
            background: #d4edda;
            color: #155724;
        }}
        
        .warning-count.low {{
            background: #fff3cd;
            color: #856404;
        }}
        
        .warning-count.high {{
            background: #f8d7da;
            color: #721c24;
        }}
        
        @media print {{
            .interactive-controls {{ display: none !important; }}
            body {{ background: white !important; color: black !important; }}
            .card {{ background: white !important; color: black !important; box-shadow: none !important; }}
            .qc-table {{ box-shadow: none !important; }}
            .table-container {{ max-height: none !important; overflow: visible !important; }}
        }}
    </style>
    {js_content}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 AcinetoScope FASTA QC - Summary Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive Quality Control Analysis Across All FASTA Files</p>
                
                <div class="interactive-controls">
                    <div class="search-box">
                        <input type="text" id="search-summary" 
                               placeholder="🔍 Search filenames or values..." 
                               onkeyup="searchTable('qc-summary-table', this.value)">
                    </div>
                    
                    <div class="filter-buttons">
                        <button class="active" onclick="filterByWarnings('all')">All Files</button>
                        <button onclick="filterByWarnings('with_warnings')">With Warnings</button>
                        <button onclick="filterByWarnings('no_warnings')">No Warnings</button>
                    </div>
                    
                    <div class="result-counter" id="result-counter-qc-summary-table">
                        {len(successful_results)} files found
                    </div>
                    
                    <div class="export-buttons">
                        <button onclick="exportToCSV()">📥 Export CSV</button>
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
        
        # Overall statistics
        total_sequences = sum(r['total_sequences'] for r in successful_results)
        total_length = sum(r['total_length'] for r in successful_results)
        total_warnings = sum(len(r.get('warnings', [])) for r in successful_results)
        files_with_warnings = sum(1 for r in successful_results if r.get('warnings'))
        
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📊 Overall Statistics</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Files</h3>
                    <p style="font-size: 2em; margin: 0;">{len(successful_results)}</p>
                </div>
                <div class="seq-stat-card">
                    <h3>Total Sequences</h3>
                    <p style="font-size: 2em; margin: 0;">{total_sequences:,}</p>
                </div>
                <div class="seq-stat-card">
                    <h3>Total Bases</h3>
                    <p style="font-size: 2em; margin: 0;">{total_length:,}</p>
                    <p style="font-size: 0.9em; margin-top: 5px;">bp</p>
                </div>
                <div class="warning-stat-card">
                    <h3>Quality Warnings</h3>
                    <p style="font-size: 2em; margin: 0;">{total_warnings}</p>
                    <p style="font-size: 0.9em; margin-top: 5px;">in {files_with_warnings} files</p>
                </div>
            </div>
            <p><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>BioPython Version:</strong> {self.metadata['biopython_version']}</p>
        </div>
"""
        
        # Detailed table with scroll
        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #3b82f6; padding-bottom: 10px;">📈 Detailed FASTA QC Results</h2>
            <p style="color: #666; margin-bottom: 15px;">Scroll horizontally to view all columns. Click headers to sort.</p>
            
            <div class="table-container">
                <table class="qc-table" id="qc-summary-table">
                    <thead>
                        <tr>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 0)">Filename</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 1)">Total Sequences</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 2)">Total Length</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 3)">Total Bases</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 4)">GC Content (%)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 5)">AT Content (%)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 6)">N50</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 7)">N75</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 8)">N90</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 9)">Median Length</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 10)">Mean Length</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 11)">Longest Sequence</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 12)">Shortest Sequence</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 13)">Ambiguous Bases (%)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 14)">Sequences with Ns</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 15)">Max N-run</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 16)">Homopolymers</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 17)">Max Homopolymer</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 18)">Duplicate Sequences</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 19)">Short Sequences (&lt;100 bp)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 20)">Long Sequences (&gt;1M bp)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 21)">File Size (MB)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 22)">Warnings</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        
        for result in successful_results:
            warning_count = len(result.get('warnings', []))
            warning_class = 'none' if warning_count == 0 else 'low' if warning_count < 3 else 'high'
            
            html_content += f"""
                        <tr>
                            <td><strong>{result['filename']}</strong></td>
                            <td>{result['total_sequences']:,}</td>
                            <td>{result['total_length']:,}</td>
                            <td>{result['total_bases']:,}</td>
                            <td>{result['gc_percent']:.1f}</td>
                            <td>{result['at_percent']:.1f}</td>
                            <td>{result['n50']:,}</td>
                            <td>{result['n75']:,}</td>
                            <td>{result['n90']:,}</td>
                            <td>{result['median_length']:,}</td>
                            <td>{result['mean_length']:.0f}</td>
                            <td>{result['longest_sequence']:,}</td>
                            <td>{result['shortest_sequence']:,}</td>
                            <td>{result['ambiguous_percent']:.2f}</td>
                            <td>{result['sequences_with_n']:,}</td>
                            <td>{result['max_n_run']:,}</td>
                            <td>{result['homopolymers_count']:,}</td>
                            <td>{result['max_homopolymer']:,}</td>
                            <td>{result['duplicate_sequences']:,}</td>
                            <td>{result['short_sequences']:,}</td>
                            <td>{result['long_sequences']:,}</td>
                            <td>{result['file_size_mb']:.2f}</td>
                            <td><span class="warning-count {warning_class}">{warning_count}</span></td>
                        </tr>
"""
        
        html_content += """
                    </tbody>
                </table>
            </div>
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
                Analysis performed using AcinetoScope FASTA QC v{self.metadata['version']}
                with BioPython {self.metadata['biopython_version']}
            </p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write HTML file
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
    
    def process_files(self, pattern: str, output_dir: str = "fasta_qc_results") -> List[Dict[str, Any]]:
        """Process multiple FASTA files with parallel execution and comprehensive reporting"""
        print("\n" + "="*80)
        print(self.ascii_art)
        print("AcinetoScope FASTA QC - Comprehensive Quality Control")
        print("="*80)
        print(f"Author: Brown Beckley | Email: brownbeckley94@gmail.com")
        print(f"Affiliation: University of Ghana Medical School - Department of Medical Biochemistry")
        print("="*80)
        print(f"Output directory: {output_dir}")
        print(f"Using {self.cpus} CPU cores")
        print("="*80)
        
        # Find FASTA files
        fasta_extensions = ['.fna', '.fasta', '.fa', '.faa']
        files = []
        
        # Try the pattern as-is first
        files.extend(glob.glob(pattern))
        
        # If no files found, try with extensions
        if not files:
            for ext in fasta_extensions:
                files.extend(glob.glob(f"{pattern}{ext}"))
        
        # Remove duplicates and non-existent files
        files = list(set([f for f in files if os.path.exists(f)]))
        
        if not files:
            self.logger.error(f"No FASTA files found matching pattern: {pattern}")
            self.logger.info(f"Supported extensions: {', '.join(fasta_extensions)}")
            return []
        
        self.logger.info(f"Found {len(files)} FASTA files: {[Path(f).name for f in files]}")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Process files in parallel
        all_results = []
        successful_files = 0
        
        with ThreadPoolExecutor(max_workers=self.cpus) as executor:
            future_to_file = {executor.submit(self.analyze_file, f): f for f in files}
            
            for future in as_completed(future_to_file):
                file = future_to_file[future]
                try:
                    result = future.result()
                    all_results.append(result)
                    
                    if result.get('status') == 'success':
                        successful_files += 1
                        # Create individual HTML report
                        self.create_individual_html_report(result, output_dir)
                    else:
                        self.logger.error(f"✗ {result['filename']}: {result.get('error', 'Unknown error')}")
                        
                except Exception as e:
                    self.logger.error(f"✗ {os.path.basename(file)}: ERROR - {str(e)}")
                    all_results.append({
                        'filename': os.path.basename(file),
                        'status': 'error',
                        'error': str(e)
                    })
        
        # Create summary reports if we have successful analyses
        if successful_files > 0:
            self.create_summary_report(all_results, output_dir)
        
        # Print final summary
        print("\n" + "="*80)
        print("ANALYSIS COMPLETE")
        print("="*80)
        print(f"Total files: {len(files)}")
        print(f"Successful analyses: {successful_files}")
        print(f"Failed analyses: {len(files) - successful_files}")
        print("\n📁 OUTPUT FILES:")
        print(f"   Individual HTML reports: {output_dir}/*/*_fasta_qc_report.html")
        print(f"   Individual JSON reports: {output_dir}/*/*_fasta_qc_report.json")
        if successful_files > 1:
            print(f"   Summary TSV: {output_dir}/FASTA_QC_summary.tsv")
            print(f"   Summary JSON: {output_dir}/FASTA_QC_summary.json")
            print(f"   Summary HTML: {output_dir}/FASTA_QC_summary.html")
        
        print("\n" + "="*80)
        
        return all_results

def main():
    """Command line interface for FASTA QC analysis"""
    parser = argparse.ArgumentParser(
        description='AcinetoScope FASTA QC - Comprehensive Quality Control with HTML Reports',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all FASTA files
  python acineto_fasta_qc.py "*.fna"
  
  # Run on specific pattern
  python acineto_fasta_qc.py "genomes/*.fasta"
  
  # Run with custom output directory
  python acineto_fasta_qc.py "*.fa" --output my_qc_results

Supported FASTA extensions: .fasta, .fa, .fna, .faa
        """
    )
    
    parser.add_argument('pattern', help='File pattern for FASTA files (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--output', '-o', default='fasta_qc_results', 
                       help='Output directory (default: fasta_qc_results)')
    parser.add_argument('--cpus', '-c', type=int, default=None,
                       help='Number of CPU cores to use (default: all available)')
    
    args = parser.parse_args()
    
    try:
        qc = AcinetoFASTAQC(cpus=args.cpus)
        results = qc.process_files(args.pattern, args.output)
        
    except Exception as e:
        print(f"FASTA QC analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()