#!/usr/bin/env python3
"""
AcinetoScope - MLST Module for Acinetobacter baumannii (UPDATED WITH NEW ICS AND STs)
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Open an issue if you found any miss matches!!
Date: 2026-08-19 (updated 2026-11-01)
"""

import os
import sys
import json
import glob
import tempfile
import argparse
import subprocess
import random
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import pandas as pd
from datetime import datetime
from html import escape

class AcinetoMLSTAnalyzer:
    def __init__(self, database_dir: Path, script_dir: Path):
        self.database_dir = database_dir
        self.script_dir = script_dir
        self.mlst_bin = script_dir / "mlst"
        
        if not self.mlst_bin.exists():
            mlst_path = shutil.which("mlst")
            if mlst_path:
                self.mlst_bin = Path(mlst_path)
            else:
                raise FileNotFoundError(f"MLST binary not found at: {self.mlst_bin}")
        
        self.has_excel_support = self.check_excel_support()
        
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
            {"text": "The universe is not required to be in harmony with human ambition.", "author": "Carl Sagan"},
        ]
        
        self.scheme_databases = {
            "oxford": "abaumannii",
            "pasteur": "abaumannii_2"
        }
        
        self.scheme_display_names = {
            "abaumannii": "OXFORD",
            "abaumannii_2": "PASTEUR"
        }
        
        # ============================================================
        # IC REFERENCE DICTIONARY - Literature sources for each IC
        # ============================================================
        self.ic_references = {
            "IC1": [
                "Zarrilli R, et al. Int J Antimicrob Agents. 2013",
                "Shelenkov A, et al. Microorganisms. 2023"
            ],
            "IC2": [
                "Zarrilli R, et al. Int J Antimicrob Agents. 2013",
                "Shelenkov A, et al. Microorganisms. 2023"
            ],
            "IC3": [
                "Gaiarsa S, et al. Front Microbiol. 2019",
                "Shelenkov A, et al. Microorganisms. 2023"
            ],
            "IC4": [
                "Gaiarsa S, et al. Front Microbiol. 2019",
                "Shelenkov A, et al. Microorganisms. 2023"
            ],
            "IC5": [
                "Gaiarsa S, et al. Front Microbiol. 2019",
                "Shelenkov A, et al. Microorganisms. 2023"
            ],
            "IC6": [
                "Gaiarsa S, et al. Front Microbiol. 2019",
                "Shelenkov A, et al. Microorganisms. 2023"
            ],
            "IC7": [
                "Gaiarsa S, et al. Front Microbiol. 2019",
                "Shelenkov A, et al. Microorganisms. 2023"
            ],
            "IC8": [
                "Gaiarsa S, et al. Front Microbiol. 2019",
                "Shelenkov A, et al. Microorganisms. 2023"
            ],
            "IC9": [
                "Gaiarsa S, et al. Front Microbiol. 2019",
                "Shelenkov A, et al. Microorganisms. 2023"
            ],
            "IC10": [
                "Shelenkov A, et al. Microorganisms. 2023"
            ],
            "IC11": [
                "Hansen F, et al. Int J Antimicrob Agents. 2023",
                "Xu A, et al. mSphere. 2024"
            ],
            "IC12": [
                "Karah N, et al. Microb Genom. 2025"
            ],
            "NOIC": [
                "No International Clone has been assigned to this ST."
            ]
        }

        # ============================================================
        # IC MAPPING - ST to IC for both schemes
        # ============================================================
        self.ic_mapping = {
            "IC1": {
                "abaumannii_2": ["1", "19", "20", "81", "94", "315", "460", "623", "717", "734", "986", "1090", "1106"],
                "abaumannii": ["81", "231"]
            },
            "IC2": {
                "abaumannii_2": ["2", "45", "187", "195", "414", "492", "570", "571", "577", "600", "604", "664", "745", "823", "1537", "1550", "1555", "1579"],
                "abaumannii": ["195", "208"]
            },
            "IC3": {
                "abaumannii_2": ["3", "124", "229", "500", "1822"],
                "abaumannii": ["452"]
            },
            "IC4": {
                "abaumannii_2": ["15", "238"],
                "abaumannii": []
            },
            "IC5": {
                "abaumannii_2": ["79", "156", "175", "422", "730", "1163", "1196"],
                "abaumannii": []
            },
            "IC6": {
                "abaumannii_2": ["78"],
                "abaumannii": []
            },
            "IC7": {
                "abaumannii_2": ["25", "113", "619", "945", "1487"],
                "abaumannii": ["619"]
            },
            "IC8": {
                "abaumannii_2": ["10", "23", "82", "575", "613", "1512"],
                "abaumannii": []
            },
            "IC9": {
                "abaumannii_2": ["6", "85", "464"],
                "abaumannii": []
            },
            "IC10": {
                "abaumannii_2": ["329", "336", "339", "340", "341", "348", "349", "350", "355", "356", "357", "358", "359", "360", "361", "362", "363", "364", "365", "366", "367", "368", "369", "370", "371", "372", "373", "374", "375", "376", "377", "378", "379", "380", "381", "382", "383", "384", "385", "386", "387", "388", "389", "390", "391", "392", "393", "394", "395", "396", "397", "398", "399", "400", "401", "402", "403", "404", "405", "406", "407", "408", "409", "410", "411", "412", "413", "414", "415", "416", "417", "418", "419", "420", "421", "422", "423", "424", "425", "426", "427", "428", "429", "430", "431", "432", "433", "434", "435", "436", "437", "438", "439", "440", "441", "442", "443", "444", "445", "446", "447", "448", "449"],
                "abaumannii": []
            },
            "IC11": {
                "abaumannii_2": ["164"],
                "abaumannii": ["1418"]
            },
            "IC12": {
                "abaumannii_2": ["158"],
                "abaumannii": []
            },
        }

        # ============================================================
        # LINEAGE DATABASE - ST → IC and lineage information
        # Based on literature: Zarrilli et al. 2013, Gaiarsa et al. 2019,
        # Shelenkov et al. 2023, Hansen et al. 2023, Karah et al. 2025
        # ============================================================
        self.lineage_db = {
            # ========== OXFORD SCHEME (abaumannii) ==========
            "abaumannii": {
                '81': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=81"
                },
                '231': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=231"
                },
                '195': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=195"
                },
                '208': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=208"
                },
                '452': {
                    "international_clone": "IC3",
                    "clonal_complex": "CC3",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly South America",
                    "clinical_significance": "Third major global clone, associated with hospital-acquired infections",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-72"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=452"
                },
                '619': {
                    "international_clone": "IC7",
                    "clonal_complex": "CC25",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Middle East, Asia",
                    "clinical_significance": "Seventh international clone, emerging threat",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["NDM-1", "OXA-23", "OXA-64"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=619"
                },
                '1418': {
                    "international_clone": "IC11",
                    "clonal_complex": "CC164",
                    "classification": "Emerging International Clone",
                    "geographic_distribution": "Denmark, China, Middle East",
                    "clinical_significance": "Emerging international clone associated with NDM-1 and OXA-23 co-production. First described in Danish patients (travel-related) and subsequently in China.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "NDM-1", "OXA-91", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=1&ST=1418"
                },
            },

            # ========== PASTEUR SCHEME (abaumannii_2) ==========
            "abaumannii_2": {
                # ---------- IC1 (CC1, OXA-69 lineage) ----------
                '1': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1"
                },
                '19': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=19"
                },
                '20': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=20"
                },
                '81': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=81"
                },
                '94': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=94"
                },
                '315': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=315"
                },
                '460': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=460"
                },
                '623': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=623"
                },
                '717': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=717"
                },
                '734': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=734"
                },
                '986': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=986"
                },
                '1090': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1090"
                },
                '1106': {
                    "international_clone": "IC1",
                    "clonal_complex": "CC1",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide",
                    "clinical_significance": "Most prevalent global clone, responsible for majority of hospital outbreaks worldwide. Often OXA-23 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-69", "ADC-30", "TEM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1106"
                },

                # ---------- IC2 (CC2, OXA-66 lineage) ----------
                '2': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=2"
                },
                '45': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=45"
                },
                '187': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=187"
                },
                '195': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=195"
                },
                '414': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=414"
                },
                '492': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=492"
                },
                '570': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=570"
                },
                '571': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=571"
                },
                '577': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=577"
                },
                '600': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=600"
                },
                '604': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=604"
                },
                '664': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=664"
                },
                '745': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=745"
                },
                '823': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=823"
                },
                '1537': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1537"
                },
                '1550': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1550"
                },
                '1555': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1555"
                },
                '1579': {
                    "international_clone": "IC2",
                    "clonal_complex": "CC2",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Asia and Middle East",
                    "clinical_significance": "Second most prevalent global clone, associated with nosocomial outbreaks. Frequently OXA-23 and OXA-58 producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "VERY HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-58", "OXA-66", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1579"
                },

                # ---------- IC3 (CC3, OXA-71 lineage) ----------
                '3': {
                    "international_clone": "IC3",
                    "clonal_complex": "CC3",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly South America and Europe",
                    "clinical_significance": "Third major global clone, often associated with OXA-72 carbapenemase.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-72", "OXA-71", "ADC-73"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=3"
                },
                '124': {
                    "international_clone": "IC3",
                    "clonal_complex": "CC3",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly South America and Europe",
                    "clinical_significance": "Third major global clone, often associated with OXA-72 carbapenemase.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-72", "OXA-71", "ADC-73"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=124"
                },
                '229': {
                    "international_clone": "IC3",
                    "clonal_complex": "CC3",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly South America and Europe",
                    "clinical_significance": "Third major global clone, often associated with OXA-72 carbapenemase.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-72", "OXA-71", "ADC-73"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=229"
                },
                '500': {
                    "international_clone": "IC3",
                    "clonal_complex": "CC3",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly South America and Europe",
                    "clinical_significance": "Third major global clone, often associated with OXA-72 carbapenemase.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-72", "OXA-71", "ADC-73"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=500"
                },
                '1822': {
                    "international_clone": "IC3",
                    "clonal_complex": "CC3",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly South America and Europe",
                    "clinical_significance": "Third major global clone, often associated with OXA-72 carbapenemase.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-72", "OXA-71", "ADC-73"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1822"
                },

                # ---------- IC4 (CC15, OXA-51 lineage) ----------
                '15': {
                    "international_clone": "IC4",
                    "clonal_complex": "CC15",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Europe",
                    "clinical_significance": "Emerging global clone, increasing in prevalence. Often associated with OXA-58.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-58", "OXA-51", "ADC-96"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=15"
                },
                '238': {
                    "international_clone": "IC4",
                    "clonal_complex": "CC15",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Worldwide, particularly Europe",
                    "clinical_significance": "Emerging global clone, increasing in prevalence. Often associated with OXA-58.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-58", "OXA-51", "ADC-96"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=238"
                },

                # ---------- IC5 (CC79, OXA-65 lineage) ----------
                '79': {
                    "international_clone": "IC5",
                    "clonal_complex": "CC79",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, Europe, Middle East",
                    "clinical_significance": "Fifth international clone, often OXA-23 producer. Emerging in Asian hospitals.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-65", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=79"
                },
                '156': {
                    "international_clone": "IC5",
                    "clonal_complex": "CC79",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, Europe, Middle East",
                    "clinical_significance": "Fifth international clone, often OXA-23 producer. Emerging in Asian hospitals.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-65", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=156"
                },
                '175': {
                    "international_clone": "IC5",
                    "clonal_complex": "CC79",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, Europe, Middle East",
                    "clinical_significance": "Fifth international clone, often OXA-23 producer. Emerging in Asian hospitals.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-65", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=175"
                },
                '422': {
                    "international_clone": "IC5",
                    "clonal_complex": "CC79",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, Europe, Middle East",
                    "clinical_significance": "Fifth international clone, often OXA-23 producer. Emerging in Asian hospitals.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-65", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=422"
                },
                '730': {
                    "international_clone": "IC5",
                    "clonal_complex": "CC79",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, Europe, Middle East",
                    "clinical_significance": "Fifth international clone, often OXA-23 producer. Emerging in Asian hospitals.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-65", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=730"
                },
                '1163': {
                    "international_clone": "IC5",
                    "clonal_complex": "CC79",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, Europe, Middle East",
                    "clinical_significance": "Fifth international clone, often OXA-23 producer. Emerging in Asian hospitals.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-65", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1163"
                },
                '1196': {
                    "international_clone": "IC5",
                    "clonal_complex": "CC79",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, Europe, Middle East",
                    "clinical_significance": "Fifth international clone, often OXA-23 producer. Emerging in Asian hospitals.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-65", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1196"
                },

                # ---------- IC6 (CC78, OXA-90 lineage) ----------
                '78': {
                    "international_clone": "IC6",
                    "clonal_complex": "CC78",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "South America, Europe",
                    "clinical_significance": "Sixth international clone, often associated with OXA-58 and colistin resistance.",
                    "common_resistance": ["Carbapenem-resistant", "Colistin-resistant", "MDR"],
                    "outbreak_potential": "MODERATE",
                    "typical_resistance_genes": ["OXA-58", "OXA-90", "mcr-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=78"
                },

                # ---------- IC7 (CC25, OXA-64 lineage) ----------
                '25': {
                    "international_clone": "IC7",
                    "clonal_complex": "CC25",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Middle East, Asia",
                    "clinical_significance": "Seventh international clone, often NDM-1 metallo-β-lactamase producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["NDM-1", "OXA-23", "OXA-64"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=25"
                },
                '113': {
                    "international_clone": "IC7",
                    "clonal_complex": "CC25",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Middle East, Asia",
                    "clinical_significance": "Seventh international clone, often NDM-1 metallo-β-lactamase producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["NDM-1", "OXA-23", "OXA-64"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=113"
                },
                '619': {
                    "international_clone": "IC7",
                    "clonal_complex": "CC25",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Middle East, Asia",
                    "clinical_significance": "Seventh international clone, often NDM-1 metallo-β-lactamase producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["NDM-1", "OXA-23", "OXA-64"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=619"
                },
                '945': {
                    "international_clone": "IC7",
                    "clonal_complex": "CC25",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Middle East, Asia",
                    "clinical_significance": "Seventh international clone, often NDM-1 metallo-β-lactamase producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["NDM-1", "OXA-23", "OXA-64"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=945"
                },
                '1487': {
                    "international_clone": "IC7",
                    "clonal_complex": "CC25",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Middle East, Asia",
                    "clinical_significance": "Seventh international clone, often NDM-1 metallo-β-lactamase producer.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["NDM-1", "OXA-23", "OXA-64"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1487"
                },

                # ---------- IC8 (CC10, OXA-68 lineage) ----------
                '10': {
                    "international_clone": "IC8",
                    "clonal_complex": "CC10",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Europe, Middle East",
                    "clinical_significance": "Eighth international clone, often pan-drug resistant with colistin resistance.",
                    "common_resistance": ["Colistin-resistant", "Carbapenem-resistant", "XDR", "PDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-68", "mcr-1", "NDM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=10"
                },
                '23': {
                    "international_clone": "IC8",
                    "clonal_complex": "CC10",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Europe, Middle East",
                    "clinical_significance": "Eighth international clone, often pan-drug resistant with colistin resistance.",
                    "common_resistance": ["Colistin-resistant", "Carbapenem-resistant", "XDR", "PDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-68", "mcr-1", "NDM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=23"
                },
                '82': {
                    "international_clone": "IC8",
                    "clonal_complex": "CC10",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Europe, Middle East",
                    "clinical_significance": "Eighth international clone, often pan-drug resistant with colistin resistance.",
                    "common_resistance": ["Colistin-resistant", "Carbapenem-resistant", "XDR", "PDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-68", "mcr-1", "NDM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=82"
                },
                '575': {
                    "international_clone": "IC8",
                    "clonal_complex": "CC10",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Europe, Middle East",
                    "clinical_significance": "Eighth international clone, often pan-drug resistant with colistin resistance.",
                    "common_resistance": ["Colistin-resistant", "Carbapenem-resistant", "XDR", "PDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-68", "mcr-1", "NDM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=575"
                },
                '613': {
                    "international_clone": "IC8",
                    "clonal_complex": "CC10",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Europe, Middle East",
                    "clinical_significance": "Eighth international clone, often pan-drug resistant with colistin resistance.",
                    "common_resistance": ["Colistin-resistant", "Carbapenem-resistant", "XDR", "PDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-68", "mcr-1", "NDM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=613"
                },
                '1512': {
                    "international_clone": "IC8",
                    "clonal_complex": "CC10",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Europe, Middle East",
                    "clinical_significance": "Eighth international clone, often pan-drug resistant with colistin resistance.",
                    "common_resistance": ["Colistin-resistant", "Carbapenem-resistant", "XDR", "PDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-68", "mcr-1", "NDM-1"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=1512"
                },

                # ---------- IC9 (CC464, OXA-94 lineage) ----------
                '6': {
                    "international_clone": "IC9",
                    "clonal_complex": "CC464",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, particularly China and India",
                    "clinical_significance": "Ninth international clone, emerging threat with OXA-23 and NDM co-production.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "NDM-1", "OXA-94"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=6"
                },
                '85': {
                    "international_clone": "IC9",
                    "clonal_complex": "CC464",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, particularly China and India",
                    "clinical_significance": "Ninth international clone, emerging threat with OXA-23 and NDM co-production.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "NDM-1", "OXA-94"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=85"
                },
                '464': {
                    "international_clone": "IC9",
                    "clonal_complex": "CC464",
                    "classification": "Global Epidemic Clone",
                    "geographic_distribution": "Asia, particularly China and India",
                    "clinical_significance": "Ninth international clone, emerging threat with OXA-23 and NDM co-production.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "NDM-1", "OXA-94"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=464"
                },

                # ---------- IC10 (ST329-449, OXA-329 lineage) ----------
                '329': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=329"
                },
                '336': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=336"
                },
                '339': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=339"
                },
                '340': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=340"
                },
                '341': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=341"
                },
                '348': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=348"
                },
                '349': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=349"
                },
                '350': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=350"
                },
                '355': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=355"
                },
                '356': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=356"
                },
                '357': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=357"
                },
                '358': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=358"
                },
                '359': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=359"
                },
                '360': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=360"
                },
                '361': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=361"
                },
                '362': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=362"
                },
                '363': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=363"
                },
                '364': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=364"
                },
                '365': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=365"
                },
                '366': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=366"
                },
                '367': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=367"
                },
                '368': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=368"
                },
                '369': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=369"
                },
                '370': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=370"
                },
                '371': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=371"
                },
                '372': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=372"
                },
                '373': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=373"
                },
                '374': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=374"
                },
                '375': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=375"
                },
                '376': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=376"
                },
                '377': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=377"
                },
                '378': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=378"
                },
                '379': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=379"
                },
                '380': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=380"
                },
                '381': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=381"
                },
                '382': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=382"
                },
                '383': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=383"
                },
                '384': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=384"
                },
                '385': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=385"
                },
                '386': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=386"
                },
                '387': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=387"
                },
                '388': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=388"
                },
                '389': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=389"
                },
                '390': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=390"
                },
                '391': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=391"
                },
                '392': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=392"
                },
                '393': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=393"
                },
                '394': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=394"
                },
                '395': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=395"
                },
                '396': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=396"
                },
                '397': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=397"
                },
                '398': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=398"
                },
                '399': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=399"
                },
                '400': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=400"
                },
                '401': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=401"
                },
                '402': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=402"
                },
                '403': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=403"
                },
                '404': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=404"
                },
                '405': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=405"
                },
                '406': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=406"
                },
                '407': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=407"
                },
                '408': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=408"
                },
                '409': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=409"
                },
                '410': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=410"
                },
                '411': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=411"
                },
                '412': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=412"
                },
                '413': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=413"
                },
                '415': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=415"
                },
                '416': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=416"
                },
                '417': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=417"
                },
                '418': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=418"
                },
                '419': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=419"
                },
                '420': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=420"
                },
                '421': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=421"
                },
                '423': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=423"
                },
                '424': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=424"
                },
                '425': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=425"
                },
                '426': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=426"
                },
                '427': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=427"
                },
                '428': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=428"
                },
                '429': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=429"
                },
                '430': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=430"
                },
                '431': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=431"
                },
                '432': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=432"
                },
                '433': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=433"
                },
                '434': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=434"
                },
                '435': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=435"
                },
                '436': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=436"
                },
                '437': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=437"
                },
                '438': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=438"
                },
                '439': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=439"
                },
                '440': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=440"
                },
                '441': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=441"
                },
                '442': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=442"
                },
                '443': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=443"
                },
                '444': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=444"
                },
                '445': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=445"
                },
                '446': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=446"
                },
                '447': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=447"
                },
                '448': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=448"
                },
                '449': {
                    "international_clone": "IC10",
                    "clonal_complex": "CC329",
                    "classification": "Emerging International Clone (Proposed)",
                    "geographic_distribution": "Currently under investigation; emerging",
                    "clinical_significance": "Proposed as the tenth international clone. Still being characterised; further surveillance needed.",
                    "common_resistance": ["Carbapenem-resistant", "MDR"],
                    "outbreak_potential": "MODERATE-HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-329", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=449"
                },

                # ---------- IC11 (ST164, OXA-91 lineage) ----------
                '164': {
                    "international_clone": "IC11",
                    "clonal_complex": "CC164",
                    "classification": "Emerging International Clone",
                    "geographic_distribution": "Denmark, China, Middle East",
                    "clinical_significance": "Emerging international clone associated with NDM-1 and OXA-23 co-production. First described in Danish patients (travel-related) and subsequently in China.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "NDM-1", "OXA-91", "ADC-30"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=164"
                },

                # ---------- IC12 (ST158, OXA-65 lineage) ----------
                '158': {
                    "international_clone": "IC12",
                    "clonal_complex": "CC158",
                    "classification": "Emerging International Clone",
                    "geographic_distribution": "Middle East (since 2007), Europe, Asia, South America",
                    "clinical_significance": "Newly described international clone predominantly found in the Middle East. Characterised by OXA-65, ADC-117, and a unique 80 kb plasmid (pIC12-2) with RP-T1 replicon. Some isolates carry tet(X3)-mediated tigecycline resistance.",
                    "common_resistance": ["Carbapenem-resistant", "MDR", "XDR"],
                    "outbreak_potential": "HIGH",
                    "typical_resistance_genes": ["OXA-23", "OXA-65", "ADC-117", "tet(X3)"],
                    "pubmlst_link": "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates&page=query&scheme_id=2&ST=158"
                },
            }
        }


        # Canonical Pasteur ST -> IC lookup derived from the supplied sample-level lookup table.
        # The source table contains repeated observations; conflicts are resolved by majority vote.
        # This compact map is loaded from a companion TSV when available, otherwise embedded below.
        self.pasteur_ic_lookup = {'437': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 43}, '2': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 30434}, '1': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 1397}, '52': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 102}, '81': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 32}, '49': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 33}, '20': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 25}, '16': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-402', 'n': 264}, '25': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-64', 'n': 561}, '404': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-89', 'n': 37}, '267': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 61}, '78': {'international_clone': 'IC6', 'oxa51_variant': 'OXA-90', 'n': 304}, '3': {'international_clone': 'IC3', 'oxa51_variant': 'OXA-71', 'n': 468}, '32': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-100', 'n': 197}, '415': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 16}, '218': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 27}, '414': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-260', 'n': 5}, '136': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-317', 'n': 13}, '513': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-90', 'n': 3}, '419': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 3}, '417': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 19}, '113': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 103}, '241': {'international_clone': 'IC11', 'oxa51_variant': 'OXA-91', 'n': 26}, '416': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 15}, '406': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-173', 'n': 521}, '413': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 2}, '422': {'international_clone': 'IC5', 'oxa51_variant': 'OXA-65', 'n': 26}, '429': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 3}, '438': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 4}, '405': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 1}, '19': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 203}, '10': {'international_clone': 'IC8', 'oxa51_variant': 'OXA-68', 'n': 280}, '515': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-90', 'n': 5}, '412': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 19}, '403': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-263', 'n': 8}, '428': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 7}, '155': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 35}, '158': {'international_clone': 'IC12', 'oxa51_variant': 'OXA-343', 'n': 39}, '431': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-117', 'n': 1}, '665': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 2}, '255': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-99', 'n': 4}, '516': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 4}, '47': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-260', 'n': 3}, '12': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 2}, '15': {'international_clone': 'IC4', 'oxa51_variant': 'OXA-51', 'n': 227}, '9': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '34': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 21}, '38': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-261', 'n': 3}, '11': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-263', 'n': 3}, '40': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 195}, '35': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-262', 'n': 2}, '36': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 14}, '37': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-78', 'n': 5}, '239': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 25}, '39': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 3}, '215': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 62}, '23': {'international_clone': 'IC8', 'oxa51_variant': 'OXA-68', 'n': 68}, '425': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '6': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-94', 'n': 11}, '427': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 4}, '33': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 54}, '79': {'international_clone': 'IC5', 'oxa51_variant': 'OXA-65', 'n': 289}, '411': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-383', 'n': 3}, '221': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 56}, '639': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 3}, '524': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '669': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 2}, '514': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 21}, '523': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '537': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-387', 'n': 3}, '500': {'international_clone': 'IC3', 'oxa51_variant': 'OXA-71', 'n': 12}, '527': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '499': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-95', 'n': 3563}, '529': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '203': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 49}, '512': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 5}, '503': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 2}, '528': {'international_clone': 'IC9', 'oxa51_variant': 'OXA-94', 'n': 1}, '536': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 9}, '240': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-374', 'n': 31}, '508': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-100', 'n': 1}, '103': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 71}, '509': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '506': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '522': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 1}, '650': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '154': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 17}, '108': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-132', 'n': 116}, '526': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '525': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '195': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 45}, '535': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 1}, '651': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 3}, '510': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 1}, '54': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-365', 'n': 13}, '671': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '532': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '652': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 7}, '521': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 1}, '285': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-337', 'n': 5}, '661': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '662': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '656': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '657': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '658': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '659': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '660': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '673': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-82', 'n': 1}, '674': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '336': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-441', 'n': 4}, '1143': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-82', 'n': 12}, '56': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 2}, '126': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 56}, '668': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '664': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 10}, '621': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '94': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 24}, '638': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 6}, '271': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 8}, '1554': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-424', 'n': 5}, '105': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 4}, '653': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 2}, '1088': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-172', 'n': 40}, '250': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-407', 'n': 70}, '464': {'international_clone': 'IC9', 'oxa51_variant': 'OXA-94', 'n': 10}, '622': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 26}, '654': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 2}, '402': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-64', 'n': 1}, '1426': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 2}, '1140': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-111', 'n': 1}, '1545': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '149': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 81}, '2480': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-144', 'n': 1}, '1549': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '492': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 88}, '498': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 5}, '1337': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 5}, '2248': {'international_clone': 'IC5', 'oxa51_variant': 'OXA-65', 'n': 15}, '2481': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-340', 'n': 1}, '85': {'international_clone': 'IC9', 'oxa51_variant': 'OXA-94', 'n': 108}, '372': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 10}, '1022': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-241', 'n': 76}, '742': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-497', 'n': 1}, '187': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 96}, '745': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 19}, '48': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 18}, '717': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 8}, '604': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 40}, '1035': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '317': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 2}, '107': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 39}, '632': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 83}, '1079': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 6}, '111': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 4}, '388': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 24}, '943': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-381', 'n': 1}, '821': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 28}, '309': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 11}, '80': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 3}, '1546': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 1}, '46': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-531', 'n': 97}, '880': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '338': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 17}, '875': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 1}, '878': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-95', 'n': 1}, '879': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '920': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '1580': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 50}, '156': {'international_clone': 'IC5', 'oxa51_variant': 'OXA-65', 'n': 111}, '738': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-65', 'n': 3}, '1168': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 2}, '1195': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 1}, '941': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 5}, '150': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-121', 'n': 96}, '810': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 5}, '1048': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 4}, '142': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 19}, '494': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 14}, '848': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 6}, '973': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-94', 'n': 4}, '216': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 29}, '2482': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '1165': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-100', 'n': 1}, '571': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 126}, '730': {'international_clone': 'IC5', 'oxa51_variant': 'OXA-65', 'n': 7}, '307': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-64', 'n': 5}, '294': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 10}, '448': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 3}, '164': {'international_clone': 'IC11', 'oxa51_variant': 'OXA-91', 'n': 460}, '138': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 13}, '923': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 2}, '927': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 12}, '1179': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-374', 'n': 1}, '690': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-374', 'n': 16}, '938': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 1}, '3289': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 1}, '575': {'international_clone': 'IC8', 'oxa51_variant': 'OXA-144', 'n': 95}, '286': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-338', 'n': 3}, '1543': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '106': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 35}, '77': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 21}, '229': {'international_clone': 'IC3', 'oxa51_variant': 'OXA-71', 'n': 25}, '626': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-123', 'n': 3}, '636': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 232}, '2167': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-132', 'n': 5}, '490': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 8}, '193': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 39}, '329': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 2}, '1156': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '1157': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 3}, '341': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-344', 'n': 2}, '1106': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-371', 'n': 5}, '1544': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '911': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 3}, '1158': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 1}, '991': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-64', 'n': 17}, '922': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 6}, '1542': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '57': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 7}, '132': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 41}, '2146': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-555', 'n': 2}, '145': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 7}, '1436': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-432', 'n': 1}, '1438': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 2}, '1437': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '152': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 20}, '570': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 119}, '1303': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '1197': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-78', 'n': 1}, '1039': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 3}, '2250': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-345', 'n': 4}, '2251': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-123', 'n': 4}, '345': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 7}, '2484': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '1512': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-68', 'n': 13}, '2485': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '400': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-100', 'n': 47}, '315': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 4}, '213': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 21}, '1130': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 1}, '1093': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-429', 'n': 12}, '963': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 1}, '1095': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 2}, '1108': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '1126': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-424', 'n': 1}, '1118': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-90', 'n': 1}, '374': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 60}, '972': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-94', 'n': 2}, '109': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 39}, '1302': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-383', 'n': 3}, '237': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-64', 'n': 6}, '1547': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2486': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '151': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 7}, '2487': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 1}, '130': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 16}, '235': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 4}, '2488': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '2489': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '1136': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 2}, '727': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 6}, '456': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 5}, '266': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 2}, '1548': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 5}, '2490': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 5}, '1133': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-402', 'n': 7}, '2491': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 2}, '2492': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '1646': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-130', 'n': 6}, '2494': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '647': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-261', 'n': 11}, '334': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 8}, '2495': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 2}, '1878': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 9}, '236': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '2034': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 16}, '331': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-441', 'n': 11}, '424': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 4}, '2238': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 11}, '2264': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 6}, '1464': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 3}, '924': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-78', 'n': 2}, '2111': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-555', 'n': 2}, '118': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 5}, '2558': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '1433': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 2}, '134': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '42': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 4}, '1102': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 7}, '328': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 9}, '623': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-371', 'n': 28}, '162': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 27}, '1448': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-64', 'n': 3}, '825': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-408', 'n': 3}, '2478': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 6}, '2497': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 2}, '1278': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-337', 'n': 1}, '172': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-431', 'n': 8}, '1472': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '3312': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 2}, '1479': {'international_clone': 'IC11', 'oxa51_variant': 'OXA-91', 'n': 12}, '1535': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-345', 'n': 1}, '600': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 57}, '2424': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '1904': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 13}, '763': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 1}, '1137': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 5}, '2499': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 2}, '1562': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-223', 'n': 2}, '1091': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-223', 'n': 203}, '1573': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 3}, '238': {'international_clone': 'IC4', 'oxa51_variant': 'OXA-51', 'n': 5}, '1447': {'international_clone': 'IC4', 'oxa51_variant': 'OXA-51', 'n': 1}, '163': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '2252': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-78', 'n': 4}, '375': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 5}, '2501': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 1}, '1555': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 19}, '2502': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 3}, '1142': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 139}, '2504': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-402', 'n': 7}, '2511': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 1}, '45': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 16}, '2512': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-223', 'n': 37}, '734': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 27}, '2513': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '724': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 15}, '577': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 11}, '2514': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 2}, '1537': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 2}, '1822': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 5}, '129': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 18}, '243': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 17}, '318': {'international_clone': 'IC4', 'oxa51_variant': 'OXA-219', 'n': 20}, '837': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 5}, '2517': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 4}, '578': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 12}, '1280': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-377', 'n': 3}, '2518': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-429', 'n': 1}, '1191': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 2}, '348': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-432', 'n': 6}, '342': {'international_clone': 'IC12', 'oxa51_variant': 'OXA-98', 'n': 3}, '1279': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-431', 'n': 1}, '740': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 4}, '1252': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '1190': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 2}, '1550': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 34}, '2522': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-555', 'n': 6}, '1111': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 10}, '280': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-377', 'n': 3}, '1230': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 13}, '1229': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '141': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 14}, '1228': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-89', 'n': 4}, '1224': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '1223': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 9}, '1225': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-346', 'n': 2}, '1227': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '462': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-262', 'n': 22}, '139': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-117', 'n': 3}, '1632': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 4}, '2524': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 1}, '2525': {'international_clone': 'IC8', 'oxa51_variant': 'OXA-68', 'n': 1}, '624': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-411', 'n': 8}, '188': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-407', 'n': 5}, '625': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-90', 'n': 9}, '966': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 1}, '968': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-531', 'n': 1}, '1115': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 2}, '970': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 1}, '962': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 2}, '1132': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-65', 'n': 1}, '53': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-100', 'n': 3}, '1104': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-528', 'n': 1}, '967': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '971': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 1}, '976': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 5}, '2194': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-341', 'n': 3}, '977': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-532', 'n': 1}, '1090': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 16}, '960': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 4}, '1487': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 11}, '2527': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 7}, '1077': {'international_clone': 'IC6', 'oxa51_variant': 'OXA-90', 'n': 19}, '1475': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 1}, '265': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 2}, '721': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-132', 'n': 7}, '2542': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '1083': {'international_clone': 'IC4', 'oxa51_variant': 'OXA-51', 'n': 2}, '1853': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-241', 'n': 3}, '613': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-68', 'n': 4}, '876': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 3}, '1357': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 3}, '764': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-429', 'n': 10}, '370': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 2}, '2559': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-94', 'n': 1}, '1756': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 3}, '2561': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-555', 'n': 2}, '2562': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '2563': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-408', 'n': 1}, '2119': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 2}, '1435': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-68', 'n': 5}, '2207': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '2523': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 5}, '2564': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 2}, '1863': {'international_clone': 'IC11', 'oxa51_variant': 'OXA-91', 'n': 1}, '648': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '649': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-562', 'n': 1}, '584': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-125', 'n': 12}, '602': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 45}, '98': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 8}, '1196': {'international_clone': 'IC5', 'oxa51_variant': 'OXA-65', 'n': 1}, '2253': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-554', 'n': 2}, '2139': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-95', 'n': 120}, '1163': {'international_clone': 'IC5', 'oxa51_variant': 'OXA-65', 'n': 4}, '975': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 3}, '1574': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 3}, '2566': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-383', 'n': 1}, '2565': {'international_clone': 'IC5', 'oxa51_variant': 'OXA-65', 'n': 4}, '729': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 2}, '2571': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '619': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-64', 'n': 11}, '823': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 2}, '459': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 12}, '2540': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-223', 'n': 2}, '157': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-68', 'n': 2}, '1579': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 104}, '2217': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-317', 'n': 4}, '2089': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 1}, '697': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 6}, '2573': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '866': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-385', 'n': 25}, '190': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '2569': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-431', 'n': 1}, '380': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-82', 'n': 1}, '367': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 4}, '1293': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '1296': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-89', 'n': 1}, '1564': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 2}, '1566': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-410', 'n': 3}, '1295': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '1298': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 1}, '1299': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '1212': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 1}, '1301': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 3}, '1297': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-431', 'n': 1}, '1294': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-383', 'n': 1}, '1211': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '351': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-99', 'n': 7}, '1577': {'international_clone': 'IC11', 'oxa51_variant': 'OXA-91', 'n': 1}, '1290': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-342', 'n': 4}, '1300': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 1}, '1578': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-408', 'n': 1}, '1567': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 2}, '444': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-121', 'n': 5}, '30': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '1117': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 2}, '2135': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 2}, '2136': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-402', 'n': 1}, '311': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 1}, '1014': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-150', 'n': 4}, '1520': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-376', 'n': 1}, '160': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 3}, '144': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 2}, '1057': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 7}, '1131': {'international_clone': 'IC11', 'oxa51_variant': 'OXA-91', 'n': 6}, '1124': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-64', 'n': 1}, '1109': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 3}, '256': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-68', 'n': 5}, '1620': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 3}, '2108': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '1411': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 3}, '284': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 2}, '3354': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-408', 'n': 4}, '762': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 2}, '2334': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 27}, '2336': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2337': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 1}, '58': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-78', 'n': 1}, '867': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '268': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-95', 'n': 5}, '2694': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 3}, '757': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2266': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 8}, '2306': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 2}, '504': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-429', 'n': 11}, '1730': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 1}, '3151': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '353': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 11}, '765': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-99', 'n': 2}, '2836': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 2}, '3336': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-68', 'n': 1}, '1828': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-441', 'n': 9}, '440': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 2}, '104': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 4}, '357': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-78', 'n': 4}, '641': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 144}, '275': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 2}, '231': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 6}, '1902': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 2}, '739': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '1112': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-89', 'n': 11}, '360': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-117', 'n': 26}, '585': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 1}, '1013': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 23}, '1074': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-95', 'n': 2}, '1742': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '1740': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '1744': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-339', 'n': 2}, '1635': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '1514': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 1}, '634': {'international_clone': 'IC4', 'oxa51_variant': 'OXA-51', 'n': 9}, '3345': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 45}, '2595': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 1}, '2711': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 2}, '2541': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-223', 'n': 13}, '1344': {'international_clone': 'IC6', 'oxa51_variant': 'OXA-90', 'n': 1}, '1627': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 2}, '1159': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 7}, '723': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '131': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 7}, '2267': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 1}, '773': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 6}, '1641': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 10}, '2265': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '439': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 1}, '2263': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2160': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '333': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-111', 'n': 3}, '2158': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-68', 'n': 1}, '2159': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-64', 'n': 1}, '858': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 3}, '2157': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 3}, '282': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 4}, '2683': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '2236': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '2177': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 1}, '2555': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 3}, '1405': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 6}, '2575': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2688': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '1908': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 1}, '1907': {'international_clone': 'IC6', 'oxa51_variant': 'OXA-90', 'n': 2}, '1017': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 5}, '2423': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '881': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 2}, '2397': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 8}, '2422': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 1}, '766': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 5}, '2332': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-400', 'n': 1}, '2213': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 4}, '2421': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 1}, '753': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-342', 'n': 1}, '466': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-409', 'n': 1}, '2420': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 1}, '2414': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 1}, '2404': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 2}, '2419': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-411', 'n': 1}, '2330': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 2}, '2416': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-78', 'n': 1}, '1059': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 1}, '2415': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-346', 'n': 1}, '2176': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-556', 'n': 5}, '935': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-261', 'n': 22}, '2407': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-385', 'n': 2}, '2413': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-390', 'n': 1}, '2276': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 4}, '932': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 2}, '2409': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-431', 'n': 1}, '2408': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '2410': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-344', 'n': 1}, '2405': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-130', 'n': 1}, '2406': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-223', 'n': 1}, '933': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-93', 'n': 3}, '22': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 3}, '2403': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-376', 'n': 1}, '1829': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '2401': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-383', 'n': 1}, '2400': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-555', 'n': 2}, '2402': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '465': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-408', 'n': 2}, '820': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '2398': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 1}, '2399': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-345', 'n': 1}, '2396': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 1}, '2395': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 2}, '2392': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-89', 'n': 1}, '2391': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2394': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '2393': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 1}, '2389': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-516', 'n': 1}, '2388': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 1}, '2390': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 1}, '1733': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-376', 'n': 1}, '2387': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-95', 'n': 1}, '2385': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-411', 'n': 1}, '2381': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-100', 'n': 1}, '2383': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-375', 'n': 1}, '2380': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '2382': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 1}, '2384': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 1}, '2378': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 1}, '1027': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 22}, '2379': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 1}, '928': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-432', 'n': 2}, '2377': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 1}, '2375': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-385', 'n': 3}, '2374': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2376': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '2373': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 2}, '2372': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-341', 'n': 2}, '2371': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-344', 'n': 1}, '2370': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-341', 'n': 1}, '2368': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '1019': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-65', 'n': 5}, '2367': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-383', 'n': 2}, '2366': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '926': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-379', 'n': 2}, '2365': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 2}, '1060': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '2364': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '2363': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-407', 'n': 1}, '461': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 3}, '2362': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-381', 'n': 1}, '2361': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 1}, '1018': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 5}, '2360': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 1}, '2359': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-130', 'n': 1}, '1526': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 3}, '1020': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-383', 'n': 1}, '2358': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-337', 'n': 1}, '1367': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-380', 'n': 2}, '1147': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-179', 'n': 3}, '2357': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2356': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 1}, '1042': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 4}, '2354': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 3}, '2355': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-179', 'n': 1}, '2353': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 3}, '2352': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-344', 'n': 1}, '814': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 2}, '2351': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '283': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 12}, '961': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '808': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-344', 'n': 6}, '2350': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-497', 'n': 1}, '2349': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2348': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 1}, '481': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-391', 'n': 1}, '2347': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-380', 'n': 1}, '812': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 3}, '2346': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '1858': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 3}, '2338': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '2333': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2328': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 1}, '2325': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-95', 'n': 5}, '934': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 8}, '2326': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 6}, '2335': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 6}, '2293': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 3}, '2324': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 1}, '2323': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-344', 'n': 1}, '2321': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-508', 'n': 2}, '2322': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2320': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '293': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-342', 'n': 1}, '2319': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '1037': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 7}, '2318': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '2188': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '1741': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '44': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 16}, '2315': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-95', 'n': 1}, '2313': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2316': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2312': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-508', 'n': 1}, '2314': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2311': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '332': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 1}, '2310': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 1}, '2309': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 1}, '2308': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '228': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-64', 'n': 1}, '2307': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-429', 'n': 1}, '936': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '1016': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 7}, '940': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '939': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-346', 'n': 4}, '473': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-412', 'n': 1}, '277': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-431', 'n': 2}, '2305': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 2}, '2304': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-383', 'n': 2}, '1058': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 2}, '2303': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '2302': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 1}, '2301': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 1}, '1511': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 2}, '2300': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-132', 'n': 1}, '2299': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-75', 'n': 1}, '2297': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 2}, '2298': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 1}, '2296': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-383', 'n': 1}, '2331': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 1}, '1069': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '2295': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-432', 'n': 1}, '274': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 3}, '2290': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 4}, '2291': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-89', 'n': 1}, '2289': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-431', 'n': 1}, '2329': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-93', 'n': 1}, '1061': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 3}, '1046': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '2292': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-123', 'n': 2}, '2288': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 1}, '2287': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-530', 'n': 1}, '1906': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2286': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '931': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 1}, '937': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-317', 'n': 1}, '2283': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-429', 'n': 1}, '2285': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-93', 'n': 1}, '2284': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '2282': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2281': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '2280': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-65', 'n': 1}, '2279': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-389', 'n': 1}, '1029': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-75', 'n': 3}, '2278': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-408', 'n': 1}, '2277': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 1}, '2275': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2274': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '2273': {'international_clone': 'IC12', 'oxa51_variant': 'OXA-98', 'n': 1}, '620': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-65', 'n': 7}, '2272': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '1036': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '1368': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 1}, '929': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-379', 'n': 1}, '2271': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-433', 'n': 2}, '930': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-390', 'n': 1}, '2270': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '2269': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-382', 'n': 1}, '2268': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-388', 'n': 1}, '1416': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-386', 'n': 1}, '2435': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 2}, '2434': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-374', 'n': 1}, '2430': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-89', 'n': 1}, '2426': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 3}, '2425': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2244': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 3}, '409': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 1}, '1422': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-411', 'n': 2}, '2090': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 1}, '1572': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 1}, '257': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '2585': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '1338': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-407', 'n': 1}, '1127': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-68', 'n': 1}, '806': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-408', 'n': 7}, '2180': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 3}, '2583': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 2}, '857': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 1}, '1100': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-317', 'n': 1}, '979': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 2}, '1866': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 5}, '2681': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 2}, '2327': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 14}, '2118': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 10}, '2317': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 3}, '1153': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-68', 'n': 4}, '1336': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 11}, '2187': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 1}, '582': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-317', 'n': 12}, '2183': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 3}, '687': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-554', 'n': 1}, '1639': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 9}, '1623': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 1}, '2186': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 1}, '1470': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 5}, '2185': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 1}, '2110': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-338', 'n': 1}, '1328': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-408', 'n': 1}, '2184': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-556', 'n': 1}, '2182': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-130', 'n': 1}, '2191': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-402', 'n': 3}, '349': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 6}, '3341': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '2620': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 2}, '2757': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-377', 'n': 4}, '3317': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '335': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-100', 'n': 1}, '245': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 4}, '1480': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-402', 'n': 1}, '1233': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 2}, '2545': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-377', 'n': 2}, '2676': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 2}, '2114': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 4}, '452': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-402', 'n': 1}, '2549': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 4}, '592': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 1}, '1450': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-89', 'n': 1}, '470': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 3}, '350': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-374', 'n': 4}, '2635': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 3}, '2643': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-508', 'n': 3}, '2633': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '2642': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 2}, '2197': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 8}, '2640': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 1}, '2576': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-389', 'n': 4}, '1463': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 2}, '946': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 5}, '2552': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-407', 'n': 1}, '2732': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 5}, '447': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 1}, '2822': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 3}, '2874': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 2}, '2150': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 3}, '1314': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-408', 'n': 1}, '1596': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 2}, '1378': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 1}, '686': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-508', 'n': 2}, '3349': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 3}, '583': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-376', 'n': 1}, '1094': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 3}, '1418': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 5}, '2703': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-424', 'n': 1}, '153': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '2493': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 2}, '2247': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 2}, '1384': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '2532': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 6}, '291': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-341', 'n': 1}, '2526': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 2}, '2651': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-431', 'n': 1}, '877': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-508', 'n': 2}, '2655': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-130', 'n': 1}, '3210': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2604': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-150', 'n': 1}, '916': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 3}, '2520': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 2}, '2418': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '2847': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '1642': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-410', 'n': 3}, '1851': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '2716': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '1509': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-69', 'n': 3}, '467': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-410', 'n': 1}, '2695': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 3}, '2199': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-337', 'n': 1}, '1043': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-241', 'n': 2}, '2828': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '2702': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 2}, '2729': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 2}, '2733': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 2}, '2543': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 3}, '2728': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 2}, '2727': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 2}, '2560': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 2}, '2140': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-80', 'n': 7}, '2726': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 2}, '1356': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 4}, '2731': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 2}, '2725': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-555', 'n': 2}, '2030': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '2031': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 3}, '2035': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 1}, '2036': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '397': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2032': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2029': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 4}, '2033': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 2}, '2037': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 3}, '2689': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-317', 'n': 2}, '2106': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '1425': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-377', 'n': 1}, '1653': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-383', 'n': 1}, '1643': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-554', 'n': 3}, '1647': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '1650': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 1}, '1645': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '1648': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-89', 'n': 2}, '1651': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-429', 'n': 1}, '1644': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 2}, '455': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-374', 'n': 1}, '1654': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 1}, '1652': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 1}, '3327': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-132', 'n': 1}, '3332': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-95', 'n': 2}, '298': {'international_clone': 'IC5', 'oxa51_variant': 'OXA-65', 'n': 1}, '2622': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 1}, '1160': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 41}, '719': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 1}, '990': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 1}, '2826': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-90', 'n': 2}, '2966': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 4}, '185': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 4}, '1377': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-68', 'n': 1}, '3299': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-89', 'n': 1}, '2132': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2116': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-339', 'n': 2}, '242': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 3}, '2115': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-124', 'n': 2}, '354': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 15}, '2708': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-71', 'n': 1}, '28': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-100', 'n': 2}, '2631': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '29': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 1}, '458': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-414', 'n': 2}, '756': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '3331': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '2145': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-402', 'n': 1}, '2143': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '2149': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 1}, '273': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 4}, '358': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 2}, '2148': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '352': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-78', 'n': 1}, '2141': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-508', 'n': 1}, '2142': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 1}, '364': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 2}, '2577': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-117', 'n': 1}, '2147': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '356': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 1}, '2578': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-432', 'n': 1}, '2823': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-88', 'n': 1}, '2819': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-376', 'n': 2}, '2820': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '2816': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 18}, '2817': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '2818': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-374', 'n': 1}, '2815': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-100', 'n': 1}, '2810': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 1}, '2812': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 1}, '2811': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 1}, '2813': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '1910': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-121', 'n': 1}, '2113': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2475': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '2474': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 8}, '2473': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-117', 'n': 1}, '2109': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2243': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '2102': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-432', 'n': 1}, '252': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '137': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 1}, '755': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 1}, '3240': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 1}, '1097': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 1}, '2519': {'international_clone': 'IC8', 'oxa51_variant': 'OXA-259', 'n': 1}, '3333': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-554', 'n': 5}, '3329': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-391', 'n': 1}, '3330': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-95', 'n': 1}, '534': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 2}, '3337': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-65', 'n': 2}, '292': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '3335': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-528', 'n': 1}, '1138': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 3}, '1412': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 2}, '1625': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-441', 'n': 1}, '454': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-259', 'n': 1}, '2698': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-179', 'n': 1}, '2648': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 1}, '2696': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-441', 'n': 5}, '770': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '2682': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2615': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2692': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '2690': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '1459': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-217', 'n': 2}, '2553': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-111', 'n': 1}, '2701': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '1362': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 3}, '760': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 2}, '2685': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2700': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 2}, '2112': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-78', 'n': 1}, '2677': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 2}, '992': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-98', 'n': 1}, '2684': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2679': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-562', 'n': 1}, '2697': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '2588': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 1}, '2693': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-555', 'n': 1}, '2680': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '2686': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 1}, '479': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-126', 'n': 1}, '3351': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 2}, '3352': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 2}, '3346': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-407', 'n': 11}, '471': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-411', 'n': 3}, '3149': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 23}, '3350': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 2}, '2339': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-144', 'n': 8}, '2833': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-424', 'n': 1}, '2832': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-377', 'n': 1}, '392': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-378', 'n': 1}, '1649': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 1}, '1365': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 1}, '2858': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 1}, '3303': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '125': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 1}, '2122': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-242', 'n': 2}, '2799': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '3348': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '3347': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '1560': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2724': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '474': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 3}, '3313': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-375', 'n': 1}, '1070': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-203', 'n': 4}, '3310': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 2}, '3309': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '2744': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 6}, '3308': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '3307': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 2}, '3306': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '2105': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '912': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 2}, '1044': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-429', 'n': 1}, '3300': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-123', 'n': 1}, '3298': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-380', 'n': 2}, '3297': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-375', 'n': 1}, '3281': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 12}, '1056': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-408', 'n': 2}, '3283': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-314', 'n': 1}, '3282': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '2603': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-343', 'n': 1}, '21': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-389', 'n': 1}, '1363': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '3270': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-391', 'n': 4}, '3302': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 1}, '3296': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-180', 'n': 2}, '3284': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-380', 'n': 1}, '3278': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-411', 'n': 2}, '3273': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-70', 'n': 4}, '1731': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-65', 'n': 1}, '3272': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-413', 'n': 5}, '1049': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-375', 'n': 1}, '3305': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '3304': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-208', 'n': 1}, '2582': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 2}, '3277': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-51', 'n': 1}, '3276': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 1}, '2580': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-91', 'n': 1}, '3322': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-241', 'n': 1}, '1359': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-67', 'n': 7}, '2260': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '3058': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-510', 'n': 1}, '368': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-558', 'n': 1}, '2743': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-430', 'n': 1}, '1581': {'international_clone': 'IC7', 'oxa51_variant': 'OXA-64', 'n': 2}, '1575': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '986': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 1}, '2803': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-64', 'n': 1}, '143': {'international_clone': 'IC2', 'oxa51_variant': 'OXA-66', 'n': 1}, '3202': {'international_clone': 'NOIC', 'oxa51_variant': '', 'n': 2}, '2806': {'international_clone': 'NOIC', 'oxa51_variant': '', 'n': 7}, '1232': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-338', 'n': 1}, '2544': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-66', 'n': 1}, '460': {'international_clone': 'IC1', 'oxa51_variant': 'OXA-69', 'n': 1}, '1432': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 3}, '1430': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-407', 'n': 3}, '2500': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '1431': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-120', 'n': 1}, '1622': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-130', 'n': 2}, '2546': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-106', 'n': 1}, '445': {'international_clone': 'IC10', 'oxa51_variant': 'OXA-120', 'n': 2}, '1028': {'international_clone': 'NOIC', 'oxa51_variant': 'OXA-346', 'n': 1}}

    def check_excel_support(self):
        try:
            import openpyxl
            return True
        except ImportError:
            return False

    def get_random_quote(self):
        return random.choice(self.science_quotes)

    def find_fasta_files(self, input_path: str) -> List[Path]:
        if os.path.isfile(input_path):
            return [Path(input_path)]

        if os.path.isdir(input_path):
            patterns = [
                "**/*.fna", "**/*.fasta", "**/*.fa", "**/*.fn",
                "**/*.fna.gz", "**/*.fasta.gz", "**/*.fa.gz",
                "**/*.gb", "**/*.gbk", "**/*.gbff"
            ]
        else:
            patterns = [input_path]

        fasta_files = []
        for pattern in patterns:
            try:
                if os.path.isdir(input_path):
                    matched_files = glob.glob(os.path.join(input_path, pattern), recursive=True)
                else:
                    matched_files = glob.glob(pattern, recursive=True)

                for file_path in matched_files:
                    path = Path(file_path)
                    if path.is_file() and self.is_fasta_file(path):
                        fasta_files.append(path)
            except Exception as e:
                print(f"⚠️ Warning: Pattern {pattern} failed: {e}")
                continue

        return sorted(list(set(fasta_files)))

    def is_fasta_file(self, file_path: Path) -> bool:
        valid_extensions = {'.fna', '.fasta', '.fa', '.fn', '.gb', '.gbk', '.gbff'}
        if file_path.suffix.lower() in valid_extensions:
            return True
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                return first_line.startswith('>')
        except:
            return False

    def run_mlst_single(self, input_file: Path, output_dir: Path, scheme: str = "abaumannii") -> Dict:
        print(f"🔬 Processing: {input_file.name} with {scheme} scheme")

        sample_output_dir = output_dir / input_file.stem
        sample_output_dir.mkdir(parents=True, exist_ok=True)

        raw_output_file = sample_output_dir / "mlst_raw_output.txt"

        mlst_cmd = [
            str(self.mlst_bin),
            str(input_file),
            "--scheme", scheme,
            "--csv",
            "--nopath"
        ]

        try:
            print(f"   Running: {' '.join(mlst_cmd)}")
            result = subprocess.run(mlst_cmd, capture_output=True, text=True, check=True)

            # Raw output is written only on failure (see except block below).

            mlst_results = self.parse_mlst_csv_fixed(result.stdout, input_file.name, scheme)

            lineage_info = self.get_lineage_info(mlst_results.get('st', 'UNKNOWN'), scheme)
            mlst_results.update(lineage_info)

            self.generate_output_files(mlst_results, sample_output_dir)

            st_display = mlst_results.get('st', 'UNKNOWN')
            ic_display = mlst_results.get('international_clone', 'NOIC')
            print(f"✅ Completed: {input_file.name} -> ST{st_display} ({ic_display})")
            return mlst_results

        except subprocess.CalledProcessError as e:
            print(f"❌ MLST failed for {input_file.name}: {e}")
            with open(raw_output_file, 'w') as f:
                f.write(f"ERROR: {e}\n")
                f.write(f"STDOUT: {e.stdout}\n")
                f.write(f"STDERR: {e.stderr}\n")

            error_result = self.get_fallback_results(input_file.name, scheme)
            self.generate_output_files(error_result, sample_output_dir)
            return error_result

    def parse_mlst_csv_fixed(self, stdout: str, sample_name: str, scheme: str) -> Dict:
        lines = stdout.strip().split('\n')
        if not lines:
            return self.get_empty_results(sample_name, scheme)

        result_line = None
        for line in reversed(lines):
            line = line.strip()
            if line and ',' in line:
                if not line.startswith('[') and not line.startswith('file'):
                    result_line = line
                    break

        if not result_line:
            return self.get_empty_results(sample_name, scheme)

        parts = []
        current_part = []
        in_quotes = False

        for char in result_line:
            if char == '"':
                in_quotes = not in_quotes
            elif char == ',' and not in_quotes:
                parts.append(''.join(current_part))
                current_part = []
            else:
                current_part.append(char)

        if current_part:
            parts.append(''.join(current_part))

        parts = [p.strip() for p in parts]

        if len(parts) < 3:
            return self.get_empty_results(sample_name, scheme)

        filename = parts[0]
        detected_scheme = parts[1]
        st_raw = parts[2]

        if st_raw == '-' or st_raw == '' or st_raw == '0':
            st = "UNKNOWN"
        elif st_raw.startswith('ST'):
            st = st_raw[2:]
        else:
            st = st_raw

        alleles = {}
        allele_parts = []

        for i in range(3, len(parts)):
            allele_str = parts[i].strip()
            if not allele_str:
                continue

            if '(' in allele_str and ')' in allele_str:
                gene = allele_str.split('(')[0].strip()
                allele = allele_str.split('(')[1].rstrip(')').strip()
                alleles[gene] = allele
                allele_parts.append(f"{gene}({allele})")
            elif allele_str:
                alleles[allele_str] = "?"
                allele_parts.append(f"{allele_str}(?)")

        allele_profile = '-'.join(allele_parts) if allele_parts else ""

        if st != "UNKNOWN" and st != "-" and st != "":
            confidence = "HIGH"
            mlst_assigned = True
        else:
            confidence = "LOW"
            mlst_assigned = False

        return {
            "sample": sample_name,
            "original_filename": filename,
            "st": st,
            "scheme": scheme,
            "scheme_display": self.scheme_display_names.get(scheme, scheme.upper()),
            "alleles": alleles,
            "allele_profile": allele_profile,
            "confidence": confidence,
            "mlst_assigned": mlst_assigned,
            "detected_scheme": detected_scheme
        }

    def get_lineage_info(self, st: str, scheme: str) -> Dict:
        """Return lineage metadata, applying the corrected Pasteur lookup consensus."""
        st_key = str(st).removeprefix("ST")
        db = self.lineage_db.get(scheme, {})
        info = db.get(st_key, {}).copy()

        if not info:
            if st_key.isdigit():
                scheme_id = "1" if scheme == "abaumannii" else "2"
                pubmlst_link = (
                    "https://pubmlst.org/bigsdb?db=pubmlst_abaumannii_isolates"
                    f"&page=query&scheme_id={scheme_id}&ST={st_key}"
                )
                info = {
                    "international_clone": "NOIC",
                    "clonal_complex": f"Unknown (ST{st_key})",
                    "classification": "Not in database - novel or uncommon",
                    "geographic_distribution": "Unknown",
                    "clinical_significance": (
                        f"ST{st_key} is not currently described in the AcinetoScope lineage "
                        "metadata. International Clone assignment is based on the supplied "
                        "Pasteur lookup table where available."
                    ),
                    "common_resistance": ["Unknown - requires further analysis"],
                    "outbreak_potential": "UNKNOWN",
                    "typical_resistance_genes": ["Unknown"],
                    "pubmlst_link": pubmlst_link,
                }
            else:
                info = {
                    "international_clone": "NOIC",
                    "clonal_complex": "Not Assigned",
                    "classification": "MLST typing failed",
                    "geographic_distribution": "N/A",
                    "clinical_significance": "Could not determine sequence type.",
                    "common_resistance": ["Cannot determine"],
                    "outbreak_potential": "UNKNOWN",
                    "typical_resistance_genes": ["Cannot determine"],
                    "pubmlst_link": "https://pubmlst.org/organisms/acinetobacter-baumannii",
                }

        # Track whether this ST was already resolved from the curated lineage_db
        from_curated_db = bool(db.get(st_key, {}))

        # Pasteur lookup: fills IC/OXA-51 for STs missing from the curated DB.
        # When both exist, the curated (literature-backed) lineage_db wins.
        if scheme == "abaumannii_2" and st_key in self.pasteur_ic_lookup:
            entry = self.pasteur_ic_lookup[st_key]
            if not from_curated_db:
                info["international_clone"] = entry["international_clone"]
                info["ic_assignment_source"] = "Supplied Pasteur ST/IC lookup table (majority consensus)"
            else:
                info["ic_assignment_source"] = "AcinetoScope curated lineage metadata (literature-backed)"
            info["oxa51_variant"] = entry["oxa51_variant"]
            info["lookup_observations"] = entry["n"]
        else:
            info.setdefault("oxa51_variant", "Not assessed")
            info.setdefault("ic_assignment_source", "AcinetoScope lineage metadata")

        # Final fallback: resolve IC via ic_mapping when absent from both sources
        if info.get("international_clone", "NOIC") == "NOIC" and st_key.isdigit():
            for ic_name, scheme_map in self.ic_mapping.items():
                if st_key in scheme_map.get(scheme, []):
                    info["international_clone"] = ic_name
                    info["ic_assignment_source"] = "AcinetoScope ic_mapping fallback"
                    break

        return info

    def get_empty_results(self, sample_name: str, scheme: str) -> Dict:
        return {
            "sample": sample_name,
            "st": "UNKNOWN",
            "scheme": scheme,
            "scheme_display": self.scheme_display_names.get(scheme, scheme.upper()),
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "mlst_assigned": False,
            "detected_scheme": "Unknown"
        }

    def get_fallback_results(self, sample_name: str, scheme: str) -> Dict:
        return {
            "sample": sample_name,
            "st": "UNKNOWN",
            "scheme": scheme,
            "scheme_display": self.scheme_display_names.get(scheme, scheme.upper()),
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "mlst_assigned": False,
            "detected_scheme": "Unknown",
            "error": "MLST analysis failed"
        }

    def generate_output_files(self, mlst_results: Dict, output_dir: Path):
        # One self-contained HTML report per sample/scheme. CSS and report data are embedded.
        self.generate_html_report(mlst_results, output_dir)

    def generate_text_report(self, mlst_results: Dict, output_dir: Path):
        sample_display = mlst_results['sample']

        report = f"""ACINETOSCOPE - MLST Analysis Report
======================================

Sample: {sample_display}
Original File: {mlst_results.get('original_filename', mlst_results['sample'])}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
MLST Scheme: {mlst_results['scheme_display']}
MLST Database: {mlst_results['scheme']}
Detected Scheme: {mlst_results.get('detected_scheme', 'Unknown')}

MLST TYPING RESULTS:
-------------------
Sequence Type (ST): {mlst_results['st']}
Confidence: {mlst_results['confidence']}
MLST Status: {'Assigned' if mlst_results['mlst_assigned'] else 'Not Assigned'}

Allele Profile ({mlst_results['scheme_display']} scheme):
{mlst_results['allele_profile']}

Detailed Alleles:
"""
        for gene, allele in mlst_results['alleles'].items():
            report += f"- {gene}: {allele}\n"

        if not mlst_results['alleles']:
            report += "- No alleles detected\n"

        ic = mlst_results.get('international_clone', 'NOIC')
        report += f"""
LINEAGE ANALYSIS:
-----------------
International Clone: {ic}
Clonal Complex: {mlst_results.get('clonal_complex', 'Unknown')}
Classification: {mlst_results.get('classification', 'Unknown')}
Geographic Distribution: {mlst_results.get('geographic_distribution', 'Unknown')}
Clinical Significance: {mlst_results.get('clinical_significance', 'Unknown')}
Outbreak Potential: {mlst_results.get('outbreak_potential', 'UNKNOWN')}

Common Resistance: {', '.join(mlst_results.get('common_resistance', ['Unknown']))}
Typical Resistance Genes: {', '.join(mlst_results.get('typical_resistance_genes', ['Unknown']))}

Literature Sources:
"""
        sources = self.ic_references.get(ic, ["No specific references available"])
        for source in sources:
            report += f"- {source}\n"

        report += f"\nPubMLST Link: {mlst_results.get('pubmlst_link', 'https://pubmlst.org/organisms/acinetobacter-baumannii')}\n"

        if 'error' in mlst_results:
            report += f"\nERROR: {mlst_results['error']}\n"

        with open(output_dir / "mlst_report.txt", 'w') as f:
            f.write(report)

    def generate_tsv_report(self, mlst_results: Dict, output_dir: Path):
        sample_display = mlst_results['sample']

        tsv_content = f"""Sample\tOriginal_File\tScheme\tMLST_Database\tST\tMLST_Status\tInternational_Clone\tClonal_Complex\tClassification\tOutbreak_Potential\tResistance_Profile\tConfidence\tAllele_Profile
{sample_display}\t{mlst_results.get('original_filename', mlst_results['sample'])}\t{mlst_results['scheme_display'].lower()}\t{mlst_results['scheme']}\t{mlst_results['st']}\t{'Assigned' if mlst_results['mlst_assigned'] else 'Not Assigned'}\t{mlst_results.get('international_clone', 'NOIC')}\t{mlst_results.get('clonal_complex', 'Unknown')}\t{mlst_results.get('classification', 'Unknown')}\t{mlst_results.get('outbreak_potential', 'UNKNOWN')}\t{', '.join(mlst_results.get('common_resistance', ['Unknown']))}\t{mlst_results['confidence']}\t{mlst_results['allele_profile']}
"""

        with open(output_dir / "mlst_report.tsv", 'w') as f:
            f.write(tsv_content)

    def generate_json_report(self, mlst_results: Dict, output_dir: Path):
        ic = mlst_results.get('international_clone', 'NOIC')
        json_report = {
            "metadata": {
                "sample": mlst_results['sample'],
                "original_filename": mlst_results.get('original_filename', mlst_results['sample']),
                "analysis_date": datetime.now().isoformat(),
                "scheme": mlst_results['scheme'],
                "scheme_display": mlst_results['scheme_display'],
                "version": "1.2.0",
                "tool": "AcinetoScope"
            },
            "mlst_results": {
                "sequence_type": mlst_results['st'],
                "confidence": mlst_results['confidence'],
                "mlst_assigned": mlst_results['mlst_assigned'],
                "allele_profile": mlst_results['allele_profile'],
                "alleles": mlst_results['alleles'],
                "detected_scheme": mlst_results.get('detected_scheme', 'Unknown')
            },
            "lineage_information": {
                "international_clone": ic,
                "clonal_complex": mlst_results.get('clonal_complex', 'Unknown'),
                "classification": mlst_results.get('classification', 'Unknown'),
                "geographic_distribution": mlst_results.get('geographic_distribution', 'Unknown'),
                "clinical_significance": mlst_results.get('clinical_significance', 'Unknown'),
                "outbreak_potential": mlst_results.get('outbreak_potential', 'UNKNOWN'),
                "common_resistance": mlst_results.get('common_resistance', ['Unknown']),
                "typical_resistance_genes": mlst_results.get('typical_resistance_genes', ['Unknown']),
                "sources": self.ic_references.get(ic, []),
                "pubmlst_link": mlst_results.get('pubmlst_link', 'https://pubmlst.org/organisms/acinetobacter-baumannii')
            }
        }

        if 'error' in mlst_results:
            json_report["error"] = mlst_results['error']

        with open(output_dir / "mlst_report.json", 'w') as f:
            json.dump(json_report, f, indent=2)

    def generate_html_report(self, mlst_results: Dict, output_dir: Path):
        random_quote = self.get_random_quote()

        sample = mlst_results['sample']
        sample_display = sample
        original_filename = mlst_results.get('original_filename', sample)
        original_display = original_filename
        st = mlst_results['st']
        scheme = mlst_results['scheme']
        scheme_display = mlst_results['scheme_display']
        confidence = mlst_results['confidence']
        international_clone = mlst_results.get('international_clone', 'NOIC')
        oxa51_variant = mlst_results.get('oxa51_variant', 'Not assessed')
        ic_assignment_source = mlst_results.get('ic_assignment_source', 'AcinetoScope lineage metadata')
        clonal_complex = mlst_results.get('clonal_complex', 'Unknown')
        allele_profile = mlst_results['allele_profile'] or "No allele profile detected"
        classification = mlst_results.get('classification', 'Unknown')
        geographic_distribution = mlst_results.get('geographic_distribution', 'Unknown')
        outbreak_potential = mlst_results.get('outbreak_potential', 'UNKNOWN')
        clinical_significance = mlst_results.get('clinical_significance', 'Unknown')
        common_resistance = mlst_results.get('common_resistance', ['Unknown'])
        typical_resistance_genes = mlst_results.get('typical_resistance_genes', ['Unknown'])
        pubmlst_link = mlst_results.get('pubmlst_link', 'https://pubmlst.org/organisms/acinetobacter-baumannii')

        ic_sources = self.ic_references.get(international_clone, ["No specific references available"])
        sources_html = ''.join(f'<li style="font-size: 0.9em; color: #4b5563; margin: 4px 0; word-wrap: break-word;">{source}</li>' for source in ic_sources)

        if international_clone in ["IC10", "IC11", "IC12"]:
            ic_note = f"🌟 {international_clone} is a newly described or emerging international clone. See literature sources below."
        elif international_clone == "NOIC":
            ic_note = "This ST has not been assigned to any known International Clone."
        else:
            ic_note = f"This strain belongs to {international_clone}, one of the globally distributed epidemic clones."

        if international_clone in ["IC10"]:
            ic_status_badge = "Proposed/Emerging"
            ic_badge_color = "#f59e0b"
        elif international_clone in ["IC11", "IC12"]:
            ic_status_badge = "Recently Described"
            ic_badge_color = "#8b5cf6"
        elif international_clone != "NOIC":
            ic_status_badge = "Established"
            ic_badge_color = "#16a34a"
        else:
            ic_status_badge = "Not Assigned"
            ic_badge_color = "#6b7280"

        alleles_html = ''
        for gene, allele in mlst_results['alleles'].items():
            alleles_html += f'''                <div class="allele-card">
                    <div style="font-size: 14px; opacity: 0.9; word-wrap: break-word;">{gene}</div>
                    <div style="font-size: 18px; word-wrap: break-word;">{allele}</div>
                </div>
'''

        if not mlst_results['alleles']:
            alleles_html = '''                <div class="allele-card" style="grid-column: 1 / -1; background: linear-gradient(135deg, #fca5a5 0%, #ef4444 100%);">
                    <div style="font-size: 14px; word-wrap: break-word;">No alleles detected</div>
                    <div style="font-size: 12px; word-wrap: break-word;">MLST typing may have failed</div>
                </div>
'''

        resistance_html = ''
        for gene in typical_resistance_genes:
            resistance_html += f'''                        <div class="resistance-card" style="word-wrap: break-word;">{gene}</div>
'''

        confidence_class = confidence.lower()

        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ACINETOSCOPE - MLST Analysis Report</title>
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

        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
        }}

        .quote-text {{
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
            color: #ffffff;
            word-wrap: break-word;
        }}

        .quote-author {{
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
            word-wrap: break-word;
        }}

        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            overflow-wrap: break-word;
        }}

        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
            word-wrap: break-word;
        }}

        .report-section h3 {{
            color: #1e40af;
            margin-top: 20px;
            margin-bottom: 10px;
            font-size: 18px;
            word-wrap: break-word;
        }}

        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}

        .metric-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
            overflow-wrap: break-word;
        }}

        .metric-label {{
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 5px;
            word-wrap: break-word;
        }}

        .metric-value {{
            font-size: 24px;
            font-weight: bold;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }}

        .allele-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}

        .allele-card {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
            text-align: center;
            font-weight: bold;
            overflow-wrap: break-word;
            word-wrap: break-word;
        }}

        .confidence-high {{ color: #16a34a; font-weight: bold; }}
        .confidence-medium {{ color: #f59e0b; font-weight: bold; }}
        .confidence-low {{ color: #dc2626; font-weight: bold; }}

        .resistance-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
            gap: 10px;
            margin-top: 15px;
        }}

        .resistance-card {{
            background: #fee2e2;
            color: #991b1b;
            padding: 12px;
            border-radius: 6px;
            text-align: center;
            font-weight: bold;
            border-left: 4px solid #ef4444;
            overflow-wrap: break-word;
            word-wrap: break-word;
        }}

        .alert-box {{
            display: flex;
            align-items: flex-start;
            gap: 15px;
            padding: 15px 20px;
            border-radius: 8px;
            margin: 15px 0;
            border-left: 6px solid;
        }}

        .alert-box .alert-icon {{
            font-size: 24px;
            flex-shrink: 0;
            margin-top: 2px;
        }}

        .alert-box .alert-content {{
            flex: 1;
        }}

        .alert-box .alert-content strong {{
            display: block;
            margin-bottom: 4px;
        }}

        .alert-box.alert-info {{
            background: #eff6ff;
            border-left-color: #3b82f6;
            color: #1e3a8a;
        }}

        .alert-box.alert-purple {{
            background: #f5f3ff;
            border-left-color: #8b5cf6;
            color: #5b21b6;
        }}

        .ic-badge {{
            display: inline-block;
            padding: 3px 12px;
            border-radius: 12px;
            font-size: 12px;
            font-weight: bold;
            margin-left: 8px;
        }}

        .ic-badge.established {{ background: #d4edda; color: #155724; }}
        .ic-badge.emerging {{ background: #fff3cd; color: #856404; }}
        .ic-badge.new {{ background: #e8d5f5; color: #4a148c; }}
        .ic-badge.noic {{ background: #e5e7eb; color: #4b5563; }}

        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
            overflow-wrap: break-word;
        }}

        .timestamp {{ color: #fbbf24; font-weight: bold; }}

        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
            overflow-wrap: break-word;
        }}

        .profile-box {{
            background: #f8fafc;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border-left: 4px solid #3b82f6;
            overflow-x: auto;
            white-space: pre-wrap;
            word-wrap: break-word;
            overflow-wrap: break-word;
        }}

        .clinical-box {{
            background: #f0f9ff;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border-left: 4px solid #0ea5e9;
            overflow-wrap: break-word;
            word-wrap: break-word;
        }}

        .sources-box {{
            background: #f8fafc;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border-left: 4px solid #8b5cf6;
        }}

        .filename-info {{
            background: #f3f4f6;
            padding: 10px;
            border-radius: 6px;
            margin: 10px 0;
            font-size: 12px;
            color: #6b7280;
            border-left: 3px solid #9ca3af;
            overflow-wrap: break-word;
            word-wrap: break-word;
        }}

        .no-truncate {{
            white-space: normal !important;
            word-wrap: break-word !important;
            overflow-wrap: break-word !important;
        }}

        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 6px; }}
            .allele-grid {{ grid-template-columns: 1fr; }}
            .metric-card {{ padding: 15px; }}
            .metrics-grid {{ grid-template-columns: 1fr; }}
            .container {{ max-width: 100%; overflow-x: auto; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art"> █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝</div>
            </div>

            <div class="quote-container">
                <div class="quote-text">"{random_quote['text']}"</div>
                <div class="quote-author">— {random_quote['author']}</div>
            </div>
        </div>

        <div class="report-section">
            <h2>📊 Sample Information</h2>
            <div class="filename-info">
                <strong>Original File:</strong> {original_display}<br>
                <strong>Display Name:</strong> {sample_display}
            </div>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sample Name</div>
                    <div class="metric-value no-truncate">{sample_display}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">MLST Scheme</div>
                    <div class="metric-value">{scheme_display}</div>
                </div>
            </div>
        </div>

        <div class="report-section">
            <h2>🧬 MLST Typing Results</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Sequence Type</div>
                    <div class="metric-value">ST{st}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">International Clone</div>
                    <div class="metric-value no-truncate">{international_clone}
                        <span class="ic-badge {{
                            'established' if '{international_clone}' in ['IC1','IC2','IC3','IC4','IC5','IC6','IC7','IC8','IC9']
                            else 'emerging' if '{international_clone}' == 'IC10'
                            else 'new' if '{international_clone}' in ['IC11','IC12']
                            else 'noic'
                        }}">{ic_status_badge}</span>
                    </div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Confidence</div>
                    <div class="metric-value confidence-{confidence_class}">{confidence}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Clonal Complex</div>
                    <div class="metric-value no-truncate">{clonal_complex}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Intrinsic OXA-51-like Variant</div>
                    <div class="metric-value no-truncate">{oxa51_variant}</div>
                </div>
            </div>
            <p class="small-note">IC assignment source: {ic_assignment_source}. For Pasteur STs represented in the supplied lookup table, the IC label and OXA-51-like variant use the table's majority consensus across repeated rows.</p>

            <div class="alert-box alert-info">
                <span class="alert-icon">💡</span>
                <div class="alert-content">
                    <strong>IC Assignment Note</strong>
                    <span style="font-size: 0.95em;">{ic_note}</span>
                </div>
            </div>

            <h3>Allele Profile</h3>
            <div class="profile-box">
                <code style="font-size: 16px; color: #1e40af; font-weight: bold; white-space: pre-wrap; word-wrap: break-word;">{allele_profile}</code>
            </div>

            <h3>Individual Alleles</h3>
            <div class="allele-grid">
{alleles_html}            </div>
        </div>

        <div class="report-section">
            <h2>🌍 Lineage Information</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Classification</div>
                    <div class="metric-value no-truncate">{classification}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Geographic Distribution</div>
                    <div class="metric-value no-truncate">{geographic_distribution}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Outbreak Potential</div>
                    <div class="metric-value no-truncate">{outbreak_potential}</div>
                </div>
            </div>

            <h3>Clinical Significance</h3>
            <div class="clinical-box">
                <p style="margin: 0; line-height: 1.6; color: #374151; word-wrap: break-word;">{clinical_significance}</p>
            </div>

            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin-top: 20px;">
                <div>
                    <h3>Resistance Profile</h3>
                    <p><strong>Common Resistance:</strong> {', '.join(common_resistance)}</p>
                </div>
                <div>
                    <h3>Typical Resistance Genes</h3>
                    <div class="resistance-grid">
{resistance_html}                    </div>
                </div>
            </div>

            <h3>Literature Sources</h3>
            <div class="sources-box">
                <ul style="margin: 0; padding-left: 20px; list-style-type: disc;">
                    {sources_html}
                </ul>
            </div>

            <div style="margin-top: 20px; text-align: center;">
                <a href="{pubmlst_link}" target="_blank" style="background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 10px 20px; border-radius: 6px; text-decoration: none; font-weight: bold; display: inline-block;">
                    🔗 View on PubMLST
                </a>
            </div>
        </div>

        <div class="footer">
            <p><strong>ACINETOSCOPE</strong> &mdash; MLST Analysis Module</p>
            <p class="timestamp">Generated: <span id="acinetoscope-timestamp"></span></p>
            <div class="authorship">
                <p><strong>Technical Support &amp; Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
            <script>document.addEventListener("DOMContentLoaded", () => document.getElementById("acinetoscope-timestamp").textContent = new Date().toLocaleString());</script>
        </div>
    </div>
</body>
</html>'''

        # Strip FASTA-style extension so we write
        #   "<accession>_mlst_report.html"  (not "<accession>.fasta_mlst_report.html")
        display_name = sample
        for _ in range(2):  # handle .fna.gz → .fna → ""
            if Path(display_name).suffix.lower() in (
                ".fna", ".fasta", ".fa", ".fn", ".gz",
                ".gb", ".gbk", ".gbff",
            ):
                display_name = Path(display_name).stem
            else:
                break
        safe_sample = "".join(c if c.isalnum() or c in "._-" else "_" for c in display_name)
        out_file = output_dir / f"{safe_sample}_mlst_report.html"
        with open(out_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

    def run_mlst_dual_scheme(self, input_file: Path, output_dir: Path) -> Dict[str, Dict]:
        print(f"\n{'='*80}")
        print(f"🧬 ACINETOSCOPE - Dual-Scheme MLST Analysis")
        print(f"{'='*80}")
        print(f"Sample: {input_file.name}")
        print(f"{'='*80}")

        main_output_dir = output_dir / "mlst_results"
        main_output_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n📁 Output directory: {main_output_dir}")

        oxford_dir = main_output_dir / "oxford_results"
        pasteur_dir = main_output_dir / "pasteur_results"
        oxford_dir.mkdir(parents=True, exist_ok=True)
        pasteur_dir.mkdir(parents=True, exist_ok=True)

        print("\n🔍 Running Oxford MLST Scheme (abaumannii)...")
        oxford_results = self.run_mlst_single(input_file, oxford_dir, "abaumannii")

        print("\n🔍 Running Pasteur MLST Scheme (abaumannii_2)...")
        pasteur_results = self.run_mlst_single(input_file, pasteur_dir, "abaumannii_2")

        self.create_combined_summary(oxford_results, pasteur_results, main_output_dir, input_file.name)

        return {
            "abaumannii": oxford_results,
            "abaumannii_2": pasteur_results
        }

    def create_combined_summary(self, oxford_results: Dict, pasteur_results: Dict, output_dir: Path, sample_name: str):
        print("\n📊 Creating Combined Summary...")

        combined_dir = output_dir / "combined_results"
        combined_dir.mkdir(parents=True, exist_ok=True)

        oxford_st = oxford_results.get('st', 'UNKNOWN')
        pasteur_st = pasteur_results.get('st', 'UNKNOWN')
        oxford_ic = oxford_results.get('international_clone', 'NOIC')
        pasteur_ic = pasteur_results.get('international_clone', 'NOIC')

        ic_match = "✅ MATCH" if oxford_ic == pasteur_ic else "⚠️ MISMATCH"

        combined_report = f"""ACINETOSCOPE - Dual-Scheme MLST Summary
==========================================

Sample: {sample_name}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

SCHEME COMPARISON:
==================

Oxford Scheme (abaumannii):
--------------------------------------------------------
Sequence Type: ST{oxford_st}
International Clone: {oxford_ic}
Confidence: {oxford_results.get('confidence', 'LOW')}
Allele Profile: {oxford_results.get('allele_profile', 'No profile')}

Pasteur Scheme (abaumannii_2):
----------------------------------------------------------
Sequence Type: ST{pasteur_st}
International Clone: {pasteur_ic}
Confidence: {pasteur_results.get('confidence', 'LOW')}
Allele Profile: {pasteur_results.get('allele_profile', 'No profile')}

COMPARISON RESULTS:
------------------
International Clone Match: {ic_match}

RECOMMENDATIONS:
---------------
"""

        if oxford_st == 'UNKNOWN' and pasteur_st == 'UNKNOWN':
            combined_report += "1. ⚠️ MLST typing failed for both schemes.\n"
            combined_report += "2. Check input FASTA file quality and completeness.\n"
        elif oxford_st != 'UNKNOWN' and pasteur_st != 'UNKNOWN':
            if oxford_ic == pasteur_ic:
                combined_report += f"1. ✅ Consistent typing: Strain belongs to {oxford_ic}.\n"
                combined_report += f"2. Use ST{pasteur_st} (Pasteur) as primary identifier.\n"
                combined_report += f"3. Corresponds to ST{oxford_st} in Oxford scheme.\n"
            else:
                combined_report += f"1. ⚠️ Scheme discrepancy detected.\n"
                combined_report += f"2. Verify with additional typing methods.\n"
                combined_report += f"3. Oxford: ST{oxford_st} ({oxford_ic})\n"
                combined_report += f"4. Pasteur: ST{pasteur_st} ({pasteur_ic})\n"
        elif oxford_st != 'UNKNOWN':
            combined_report += f"1. Use Oxford ST{oxford_st} as primary result.\n"
            combined_report += f"2. Pasteur scheme failed or returned unknown.\n"
            combined_report += f"3. International Clone: {oxford_ic}\n"
        else:
            combined_report += f"1. Use Pasteur ST{pasteur_st} as primary result.\n"
            combined_report += f"2. Oxford scheme failed or returned unknown.\n"
            combined_report += f"3. International Clone: {pasteur_ic}\n"

        safe_sample = "".join(c if c.isalnum() or c in "._-" else "_" for c in sample_name)
        match_label = "MATCH" if oxford_ic == pasteur_ic else "MISMATCH"
        match_class = "match" if oxford_ic == pasteur_ic else "mismatch"

        combined_html = (
            "<!DOCTYPE html><html lang='en'><head><meta charset='UTF-8'>"
            f"<title>ACINETOSCOPE - Combined MLST Summary - {sample_name}</title>"
            "<style>"
            "body{font-family:'Segoe UI',sans-serif;background:linear-gradient(135deg,#1a1a2e 0%,#16213e 50%,#0f3460 100%);color:#fff;padding:20px;min-height:100vh;}"
            ".container{max-width:1200px;margin:0 auto;}"
            ".report-section{background:rgba(255,255,255,.95);color:#1f2937;padding:25px;border-radius:10px;margin-bottom:20px;box-shadow:0 4px 15px rgba(0,0,0,.2);}"
            "h1,h2{color:#1e3a8a;}h2{border-bottom:3px solid #3b82f6;padding-bottom:10px;}"
            ".match{color:#16a34a;font-weight:bold;}.mismatch{color:#dc2626;font-weight:bold;}"
            "pre{background:#f8fafc;padding:15px;border-radius:6px;border-left:4px solid #3b82f6;overflow-x:auto;white-space:pre-wrap;word-wrap:break-word;font-size:13px;color:#1f2937;}"
            ".footer{text-align:center;margin-top:30px;padding:20px;background:rgba(0,0,0,.3);border-radius:10px;font-size:14px;}"
            ".timestamp{color:#fbbf24;font-weight:bold;}"
            "</style></head><body><div class='container'>"
            "<h1>ACINETOSCOPE - Combined Dual-Scheme MLST Summary</h1>"
            f"<p style='color:#d1d5db;'>Sample: <strong>{sample_name}</strong> | Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>"
            "<div class='report-section'>"
            f"<h2>IC Match: <span class='{match_class}'>{match_label}</span></h2>"
            f"<pre>{escape(combined_report)}</pre>"
            "</div>"
            "<div class='footer'>"
            "<p><strong>ACINETOSCOPE</strong> - A. baumannii MLST Analysis Module</p>"
            f"<p class='timestamp'>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>"
            "<div style='margin-top:15px;padding:15px;background:rgba(255,255,255,.1);border-radius:8px;font-size:12px;'>"
            "<p><strong>Technical Support & Inquiries:</strong></p>"
            "<p>Author: Brown Beckley | GitHub: bbeckley-hub</p>"
            "<p>Email: brownbeckley94@gmail.com</p>"
            "<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>"
            "</div></div></div></body></html>"
        )
        with open(combined_dir / f"{safe_sample}_combined_summary.html", 'w', encoding='utf-8') as f:
            f.write(combined_html)

        print(f"✅ Combined HTML summary created: {combined_dir}/")

    def run_mlst_dual_batch(self, input_path: str, output_dir: Path) -> Dict[str, Dict]:
        print("🔍 Searching for FASTA files...")
        fasta_files = self.find_fasta_files(input_path)

        if not fasta_files:
            print("❌ No FASTA files found!")
            return {}

        print(f"📁 Found {len(fasta_files)} FASTA files")
        for i, fasta_file in enumerate(fasta_files, 1):
            print(f"  {i}. {fasta_file}")

        main_output_dir = output_dir / "dual_scheme_mlst"
        main_output_dir.mkdir(parents=True, exist_ok=True)

        results = {}
        for fasta_file in fasta_files:
            sample_base = fasta_file.stem
            sample_output_dir = main_output_dir / sample_base
            sample_output_dir.mkdir(parents=True, exist_ok=True)

            print(f"\n{'='*80}")
            print(f"🧬 Processing: {fasta_file.name} with dual schemes")
            print(f"{'='*80}")

            oxford_dir = sample_output_dir / "Oxford"
            pasteur_dir = sample_output_dir / "Pasteur"
            oxford_dir.mkdir(parents=True, exist_ok=True)
            pasteur_dir.mkdir(parents=True, exist_ok=True)

            print("\n🔍 Running Oxford MLST Scheme (abaumannii)...")
            oxford_results = self.run_mlst_single(fasta_file, oxford_dir, "abaumannii")

            print("\n🔍 Running Pasteur MLST Scheme (abaumannii_2)...")
            pasteur_results = self.run_mlst_single(fasta_file, pasteur_dir, "abaumannii_2")

            self.create_combined_summary(oxford_results, pasteur_results, sample_output_dir, fasta_file.name)

            results[fasta_file.name] = {
                "abaumannii": oxford_results,
                "abaumannii_2": pasteur_results
            }

            st_display_oxford = oxford_results.get('st', 'UNKNOWN')
            st_display_pasteur = pasteur_results.get('st', 'UNKNOWN')
            print(f"✅ Completed: {fasta_file.name} -> Oxford:ST{st_display_oxford}, Pasteur:ST{st_display_pasteur}")

        self.create_dual_batch_summary(results, main_output_dir)

        return results

    def create_dual_batch_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        print(f"\n📊 Creating Dual-Scheme Batch Summary...")

        summary_data = []
        detailed_results = []

        for sample_name, schemes in all_results.items():
            oxford_results = schemes.get("abaumannii", {})
            pasteur_results = schemes.get("abaumannii_2", {})

            summary_data.append({
                "Sample": sample_name,
                "Original_Filename": oxford_results.get('original_filename', sample_name),
                "Oxford_ST": oxford_results.get('st', 'UNKNOWN'),
                "Pasteur_ST": pasteur_results.get('st', 'UNKNOWN'),
                "Oxford_IC": oxford_results.get('international_clone', 'NOIC'),
                "Pasteur_IC": pasteur_results.get('international_clone', 'NOIC'),
                "Oxford_Confidence": oxford_results.get('confidence', 'LOW'),
                "Pasteur_Confidence": pasteur_results.get('confidence', 'LOW'),
                "Oxford_Allele_Profile": oxford_results.get('allele_profile', ''),
                "Pasteur_Allele_Profile": pasteur_results.get('allele_profile', ''),
                "IC_Match": "MATCH" if oxford_results.get('international_clone') == pasteur_results.get('international_clone') else "MISMATCH"
            })

            detailed_results.append({
                "sample": sample_name,
                "original_filename": oxford_results.get('original_filename', sample_name),
                "oxford_results": {
                    "sequence_type": oxford_results.get('st', 'UNKNOWN'),
                    "confidence": oxford_results.get('confidence', 'LOW'),
                    "mlst_assigned": oxford_results.get('mlst_assigned', False),
                    "allele_profile": oxford_results.get('allele_profile', ''),
                    "alleles": oxford_results.get('alleles', {}),
                    "international_clone": oxford_results.get('international_clone', 'NOIC')
                },
                "pasteur_results": {
                    "sequence_type": pasteur_results.get('st', 'UNKNOWN'),
                    "confidence": pasteur_results.get('confidence', 'LOW'),
                    "mlst_assigned": pasteur_results.get('mlst_assigned', False),
                    "allele_profile": pasteur_results.get('allele_profile', ''),
                    "alleles": pasteur_results.get('alleles', {}),
                    "international_clone": pasteur_results.get('international_clone', 'NOIC')
                }
            })

        df = pd.DataFrame(summary_data)
        csv_file = output_dir / "dual_scheme_mlst_summary.csv"
        df.to_csv(csv_file, index=False)
        print(f"  ✅ CSV summary: {csv_file}")

        if self.has_excel_support:
            try:
                excel_file = output_dir / "dual_scheme_mlst_summary.xlsx"
                df.to_excel(excel_file, index=False)
                print(f"  ✅ Excel summary: {excel_file}")
            except Exception as e:
                print(f"  ⚠️ Excel export failed: {e}")
        else:
            print(f"  ⚠️ Excel export not available (install openpyxl: pip install openpyxl)")

        json_summary = {
            "metadata": {
                "analysis_date": datetime.now().isoformat(),
                "schemes": ["Oxford (abaumannii)", "Pasteur (abaumannii_2)"],
                "samples_analyzed": len(all_results),
                "version": "1.1.0",
                "tool": "AcinetoScope"
            },
            "results": detailed_results,
            "summary_statistics": {
                "total_samples": len(all_results),
                "oxford_assigned": sum(1 for r in all_results.values() if r.get("abaumannii", {}).get('mlst_assigned', False)),
                "pasteur_assigned": sum(1 for r in all_results.values() if r.get("abaumannii_2", {}).get('mlst_assigned', False)),
                "ic_matches": sum(1 for r in summary_data if r['IC_Match'] == 'MATCH'),
                "ic_mismatches": sum(1 for r in summary_data if r['IC_Match'] == 'MISMATCH')
            }
        }

        json_file = output_dir / "dual_scheme_mlst_summary.json"
        with open(json_file, 'w') as f:
            json.dump(json_summary, f, indent=2)
        print(f"  ✅ JSON summary: {json_file}")

        self.create_dual_batch_html_summary(summary_data, output_dir)

        print(f"✅ Dual-scheme batch summary created in: {output_dir}/")

    def create_dual_batch_html_summary(self, summary_data: List[Dict], output_dir: Path):
        total_samples = len(summary_data)
        ic_matches = sum(1 for row in summary_data if row['IC_Match'] == 'MATCH')
        ic_mismatches = total_samples - ic_matches

        table_rows = ""
        for i, row in enumerate(summary_data, 1):
            sample_name = row['Sample']
            sample_display = sample_name
            ic_match_class = "match" if row['IC_Match'] == 'MATCH' else "mismatch"

            table_rows += f"""
            <tr>
                <td>{i}</td>
                <td><span class="sample-name no-truncate" title="{sample_name}">{sample_display}</span></td>
                <td>ST{row['Oxford_ST']}</td>
                <td>ST{row['Pasteur_ST']}</td>
                <td class="no-truncate">{row['Oxford_IC']}</td>
                <td class="no-truncate">{row['Pasteur_IC']}</td>
                <td class="{ic_match_class} no-truncate">{row['IC_Match']}</td>
                <td><code class="no-truncate" style="white-space: pre-wrap; word-wrap: break-word;">{row['Oxford_Allele_Profile']}</code></td>
                <td><code class="no-truncate" style="white-space: pre-wrap; word-wrap: break-word;">{row['Pasteur_Allele_Profile']}</code></td>
            </tr>
            """

        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ACINETOSCOPE - Dual-Scheme MLST Batch Summary</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 2200px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 15px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 8px;
            line-height: 1.1;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 25px;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        .stat-value {{ font-size: 28px; font-weight: bold; margin-bottom: 5px; }}
        .stat-label {{ font-size: 12px; opacity: 0.9; }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            overflow-x: auto;
        }}
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 15px;
            font-size: 22px;
        }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 10px; table-layout: auto; }}
        th {{
            background: #3b82f6;
            color: white;
            padding: 12px 15px;
            text-align: left;
            position: sticky;
            top: 0;
            white-space: nowrap;
        }}
        td {{ padding: 10px 15px; border-bottom: 1px solid #e5e7eb; vertical-align: top; }}
        tr:hover {{ background-color: #f3f4f6; }}
        .match {{ color: #16a34a; font-weight: bold; }}
        .mismatch {{ color: #dc2626; font-weight: bold; }}
        .sample-name {{ cursor: help; }}
        .no-truncate {{ white-space: normal !important; word-wrap: break-word !important; overflow-wrap: break-word !important; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .download-links {{ display: flex; justify-content: center; gap: 15px; margin: 20px 0; flex-wrap: wrap; }}
        .download-btn {{
            background: linear-gradient(135deg, #10b981 0%, #059669 100%);
            color: white;
            padding: 10px 20px;
            border-radius: 6px;
            text-decoration: none;
            font-weight: bold;
            display: inline-flex;
            align-items: center;
            gap: 8px;
        }}
        .download-btn:hover {{ transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2); }}
        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 5px; }}
            table {{ display: block; overflow-x: auto; }}
            th, td {{ padding: 8px 10px; font-size: 14px; }}
            .download-links {{ flex-direction: column; align-items: center; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art"> █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝</div>
            </div>
            <h1 style="color: white; margin-bottom: 10px;">Dual-Scheme MLST Batch Analysis</h1>
            <p style="color: #d1d5db; margin-bottom: 20px;">Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>

        <div class="stats-grid">
            <div class="stat-card"><div class="stat-value">{total_samples}</div><div class="stat-label">SAMPLES ANALYZED</div></div>
            <div class="stat-card"><div class="stat-value">{ic_matches}</div><div class="stat-label">IC MATCHES</div></div>
            <div class="stat-card"><div class="stat-value">{ic_mismatches}</div><div class="stat-label">IC MISMATCHES</div></div>
        </div>

        <div class="download-links">
            <a href="dual_scheme_mlst_summary.csv" class="download-btn">📄 Download CSV Summary</a>
            <a href="dual_scheme_mlst_summary.json" class="download-btn">🔧 Download JSON Summary</a>
'''

        if self.has_excel_support:
            html_content += f'''            <a href="dual_scheme_mlst_summary.xlsx" class="download-btn">📊 Download Excel Summary</a>
'''

        html_content += f'''        </div>

        <div class="report-section">
            <h2>📋 Dual-Scheme Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>#</th>
                        <th>Sample</th>
                        <th>Oxford ST</th>
                        <th>Pasteur ST</th>
                        <th>Oxford IC</th>
                        <th>Pasteur IC</th>
                        <th>IC Match</th>
                        <th>Oxford Allele Profile</th>
                        <th>Pasteur Allele Profile</th>
                    </tr>
                </thead>
                <tbody>
                    {table_rows}
                </tbody>
            </table>
        </div>

        <div class="footer">
            <p><strong>ACINETOSCOPE</strong> - A. baumannii MLST Analysis Module</p>
            <p class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>For individual sample reports, navigate to respective sample directories.</p>
            <div style="margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1); border-radius: 8px; font-size: 12px;">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>
</body>
</html>'''

        html_file = output_dir / "dual_scheme_mlst_summary.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"  ✅ HTML summary: {html_file}")

    def run_mlst_batch(self, input_path: str, output_dir: Path, scheme: str = "abaumannii") -> Dict[str, Dict]:
        print("🔍 Searching for FASTA files...")
        fasta_files = self.find_fasta_files(input_path)

        if not fasta_files:
            print("❌ No FASTA files found!")
            return {}

        print(f"📁 Found {len(fasta_files)} FASTA files")
        for i, fasta_file in enumerate(fasta_files, 1):
            print(f"  {i}. {fasta_file}")

        scheme_display = self.scheme_display_names.get(scheme, scheme.upper())
        scheme_output_dir = output_dir / f"{scheme_display.lower()}_results"
        scheme_output_dir.mkdir(parents=True, exist_ok=True)

        results = {}
        for fasta_file in fasta_files:
            sample_base = fasta_file.stem
            sample_output_dir = scheme_output_dir / sample_base
            sample_output_dir.mkdir(parents=True, exist_ok=True)

            result = self.run_mlst_single(fasta_file, sample_output_dir, scheme)
            results[fasta_file.name] = result

        self.create_scheme_summary(results, scheme_output_dir, scheme)

        return results

    def create_scheme_summary(self, all_results: Dict[str, Dict], output_dir: Path, scheme: str):
        print(f"\n📊 Creating {scheme} Summary...")

        summary_data = []
        detailed_results = []

        for sample_name, result in all_results.items():
            summary_data.append({
                "Sample": sample_name,
                "Original_Filename": result.get('original_filename', sample_name),
                "ST": result.get('st', 'UNKNOWN'),
                "Scheme": result.get('scheme_display', 'Unknown'),
                "International_Clone": result.get('international_clone', 'NOIC'),
                "Clonal_Complex": result.get('clonal_complex', 'Unknown'),
                "Classification": result.get('classification', 'Unknown'),
                "Confidence": result.get('confidence', 'LOW'),
                "Allele_Profile": result.get('allele_profile', ''),
                "MLST_Assigned": result.get('mlst_assigned', False)
            })

            detailed_results.append({
                "sample": sample_name,
                "original_filename": result.get('original_filename', sample_name),
                "mlst_results": {
                    "sequence_type": result.get('st', 'UNKNOWN'),
                    "confidence": result.get('confidence', 'LOW'),
                    "mlst_assigned": result.get('mlst_assigned', False),
                    "allele_profile": result.get('allele_profile', ''),
                    "alleles": result.get('alleles', {})
                },
                "lineage_information": {
                    "international_clone": result.get('international_clone', 'NOIC'),
                    "clonal_complex": result.get('clonal_complex', 'Unknown'),
                    "classification": result.get('classification', 'Unknown'),
                    "geographic_distribution": result.get('geographic_distribution', 'Unknown'),
                    "outbreak_potential": result.get('outbreak_potential', 'UNKNOWN')
                }
            })

        df = pd.DataFrame(summary_data)
        csv_file = output_dir / f"{self.scheme_display_names.get(scheme, scheme).lower()}_mlst_summary.csv"
        df.to_csv(csv_file, index=False)
        print(f"  ✅ CSV summary: {csv_file}")

        if self.has_excel_support:
            try:
                excel_file = output_dir / f"{self.scheme_display_names.get(scheme, scheme).lower()}_mlst_summary.xlsx"
                df.to_excel(excel_file, index=False)
                print(f"  ✅ Excel summary: {excel_file}")
            except Exception as e:
                print(f"  ⚠️ Excel export failed: {e}")
        else:
            print(f"  ⚠️ Excel export not available (install openpyxl: pip install openpyxl)")

        json_summary = {
            "metadata": {
                "analysis_date": datetime.now().isoformat(),
                "scheme": scheme,
                "scheme_display": self.scheme_display_names.get(scheme, scheme.upper()),
                "samples_analyzed": len(all_results),
                "version": "1.1.0",
                "tool": "AcinetoScope"
            },
            "results": detailed_results,
            "summary_statistics": {
                "total_samples": len(all_results),
                "assigned_st": sum(1 for r in all_results.values() if r.get('mlst_assigned', False)),
                "unknown_st": sum(1 for r in all_results.values() if not r.get('mlst_assigned', False)),
                "high_confidence": sum(1 for r in all_results.values() if r.get('confidence') == 'HIGH'),
                "low_confidence": sum(1 for r in all_results.values() if r.get('confidence') == 'LOW')
            }
        }

        json_file = output_dir / f"{self.scheme_display_names.get(scheme, scheme).lower()}_mlst_summary.json"
        with open(json_file, 'w') as f:
            json.dump(json_summary, f, indent=2)
        print(f"  ✅ JSON summary: {json_file}")

        self.create_html_batch_summary(summary_data, output_dir, scheme)

        print(f"✅ {self.scheme_display_names.get(scheme, scheme)} summary created in: {output_dir}/")

    def create_html_batch_summary(self, summary_data: List[Dict], output_dir: Path, scheme: str):
        scheme_display = self.scheme_display_names.get(scheme, scheme.upper())

        total_samples = len(summary_data)
        assigned_st = sum(1 for r in summary_data if r['ST'] != 'UNKNOWN')
        unknown_st = total_samples - assigned_st
        high_confidence = sum(1 for r in summary_data if r['Confidence'] == 'HIGH')

        table_rows = ""
        for i, row in enumerate(summary_data, 1):
            sample_name = row['Sample']
            sample_display = sample_name
            st_class = "st-assigned" if row['ST'] != 'UNKNOWN' else "st-unknown"
            confidence_class = row['Confidence'].lower()

            table_rows += f"""
            <tr>
                <td>{i}</td>
                <td><span class="sample-name no-truncate" title="{sample_name}">{sample_display}</span></td>
                <td class="{st_class}">ST{row['ST']}</td>
                <td class="no-truncate">{row['International_Clone']}</td>
                <td class="confidence-{confidence_class}">{row['Confidence']}</td>
                <td><code class="no-truncate" style="white-space: pre-wrap; word-wrap: break-word;">{row['Allele_Profile']}</code></td>
            </tr>
            """

        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ACINETOSCOPE - {scheme_display} MLST Batch Summary</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 2000px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 15px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 8px;
            line-height: 1.1;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            overflow-x: auto;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 25px;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        .stat-value {{ font-size: 28px; font-weight: bold; margin-bottom: 5px; }}
        .stat-label {{ font-size: 12px; opacity: 0.9; }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            overflow-x: auto;
        }}
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 15px;
            font-size: 22px;
        }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 10px; table-layout: auto; }}
        th {{
            background: #3b82f6;
            color: white;
            padding: 12px 15px;
            text-align: left;
            position: sticky;
            top: 0;
            white-space: nowrap;
        }}
        td {{ padding: 10px 15px; border-bottom: 1px solid #e5e7eb; vertical-align: top; }}
        tr:hover {{ background-color: #f3f4f6; }}
        .st-assigned {{ color: #16a34a; font-weight: bold; }}
        .st-unknown {{ color: #dc2626; font-weight: bold; }}
        .confidence-high {{ color: #16a34a; font-weight: bold; }}
        .confidence-low {{ color: #dc2626; font-weight: bold; }}
        .sample-name {{ cursor: help; }}
        .no-truncate {{ white-space: normal !important; word-wrap: break-word !important; overflow-wrap: break-word !important; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .download-links {{ display: flex; justify-content: center; gap: 15px; margin: 20px 0; flex-wrap: wrap; }}
        .download-btn {{
            background: linear-gradient(135deg, #10b981 0%, #059669 100%);
            color: white;
            padding: 10px 20px;
            border-radius: 6px;
            text-decoration: none;
            font-weight: bold;
            display: inline-flex;
            align-items: center;
            gap: 8px;
        }}
        .download-btn:hover {{ transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2); }}
        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 5px; }}
            table {{ display: block; overflow-x: auto; }}
            th, td {{ padding: 8px 10px; font-size: 14px; }}
            .download-links {{ flex-direction: column; align-items: center; }}
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
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art"> █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝</div>
            </div>
            <h1 style="color: white; margin-bottom: 10px;">{scheme_display} MLST Batch Analysis</h1>
            <p style="color: #d1d5db; margin-bottom: 20px;">Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>

        <div class="stats-grid">
            <div class="stat-card"><div class="stat-value">{total_samples}</div><div class="stat-label">SAMPLES ANALYZED</div></div>
            <div class="stat-card"><div class="stat-value">{assigned_st}</div><div class="stat-label">STs ASSIGNED</div></div>
            <div class="stat-card"><div class="stat-value">{unknown_st}</div><div class="stat-label">UNKNOWN STs</div></div>
            <div class="stat-card"><div class="stat-value">{high_confidence}</div><div class="stat-label">HIGH CONFIDENCE</div></div>
        </div>

        <div class="download-links">
            <a href="{scheme_display.lower()}_mlst_summary.csv" class="download-btn">📄 Download CSV Summary</a>
            <a href="{scheme_display.lower()}_mlst_summary.json" class="download-btn">🔧 Download JSON Summary</a>
'''

        if self.has_excel_support:
            html_content += f'''            <a href="{scheme_display.lower()}_mlst_summary.xlsx" class="download-btn">📊 Download Excel Summary</a>
'''

        html_content += f'''        </div>

        <div class="report-section">
            <h2>📋 Detailed Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>#</th>
                        <th>Sample</th>
                        <th>ST</th>
                        <th>International Clone</th>
                        <th>Confidence</th>
                        <th>Allele Profile</th>
                    </tr>
                </thead>
                <tbody>
                    {table_rows}
                </tbody>
            </table>
        </div>

        <div class="footer">
            <p><strong>ACINETOSCOPE</strong> - A. baumannii MLST Analysis Module</p>
            <p class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>For individual sample reports, navigate to respective sample directories.</p>
            <div style="margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1); border-radius: 8px; font-size: 12px;">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>
</body>
</html>'''

        html_file = output_dir / f"{self.scheme_display_names.get(scheme, scheme).lower()}_mlst_summary.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"  ✅ HTML summary: {html_file}")


def update_mlst_database(db_dir: Path, script_dir: Path) -> bool:
    """
    Download / refresh PubMLST Acinetobacter baumannii MLST schemes using mlstdb.

    Both schemes are preserved because AcinetoScope dual-scheme analysis needs
    Oxford (``abaumannii``) AND Pasteur (``abaumannii_2``) simultaneously.

    Returns True on success, False on any failure.
    """
    if not shutil.which('mlstdb'):
        print("❌ 'mlstdb' not found. Install with: pip install mlstdb")
        print("   Then configure PubMLST access: mlstdb connect --db pubmlst --api-key")
        return False

    db_dir = Path(db_dir).resolve()
    db_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix="acinetoscope_mlst_update_"))
    print(f"📁 Temporary download directory: {temp_dir}")

    try:
        print("⬇️  Running mlstdb update (full download, 10-20 min)...")
        result = subprocess.run(
            ['mlstdb', 'update'],
            input='y\n',
            text=True,
            cwd=temp_dir,
            env=os.environ.copy(),
        )
        if result.returncode != 0:
            print(f"❌ mlstdb update failed with exit code {result.returncode}")
            return False
        print("\n✅ Download completed.")

        temp_pubmlst = temp_dir / "pubmlst"
        if not temp_pubmlst.exists():
            print("❌ pubmlst directory not found in downloaded data.")
            return False

        target_pubmlst = db_dir / "pubmlst"
        target_pubmlst.mkdir(parents=True, exist_ok=True)

        schemes_needed = ['abaumannii', 'abaumannii_2']
        aliases = {
            'abaumannii':   ['abaumannii', 'acinetobacter_baumannii', 'abaumannii_oxford'],
            'abaumannii_2': ['abaumannii_2', 'acinetobacter_baumannii_2', 'abaumannii_pasteur'],
        }
        installed = []
        for scheme_name in schemes_needed:
            src_dir = None
            for alias in aliases[scheme_name]:
                candidate = temp_pubmlst / alias
                if candidate.exists() and candidate.is_dir():
                    src_dir = candidate
                    if alias != scheme_name:
                        print(f"  Source folder '{alias}' → installed as '{scheme_name}'")
                    break
            if src_dir is None:
                print(f"  ⚠️  Scheme '{scheme_name}' not found. Available:")
                for d in sorted(p.name for p in temp_pubmlst.iterdir() if p.is_dir()):
                    print(f"     - {d}")
                continue
            target = target_pubmlst / scheme_name
            if target.exists():
                shutil.rmtree(target)
            shutil.copytree(src_dir, target)
            installed.append(scheme_name)
            print(f"  ✓ Installed: {scheme_name}")

        if not installed:
            print("❌ No A. baumannii schemes could be installed.")
            return False

        # BLAST db from mlstdb download
        blast_src = temp_dir / "blast"
        if blast_src.exists():
            target_blast = db_dir / "blast"
            if target_blast.exists():
                shutil.rmtree(target_blast)
            shutil.copytree(blast_src, target_blast)
            print(f"  ✓ BLAST database copied to {target_blast}")

        # Prune any other schemes so the DB stays lean
        removed = 0
        for item in target_pubmlst.iterdir():
            if item.is_dir() and item.name not in schemes_needed:
                shutil.rmtree(item)
                removed += 1
                print(f"  🗑  Removed scheme: {item.name}")
        if removed:
            print(f"  Removed {removed} other scheme(s)")

        # Optional BLAST rebuild via mlst-make_blast_db
        rebuild_script = Path(script_dir).resolve().parent / "scripts" / "mlst-make_blast_db"
        if not rebuild_script.exists():
            alt = Path(script_dir).resolve() / "scripts" / "mlst-make_blast_db"
            if alt.exists():
                rebuild_script = alt
        if rebuild_script.exists():
            print("🔨 Rebuilding BLAST database...")
            rebuild_env = os.environ.copy()
            rebuild_env["MLST_DATADIR"] = str(db_dir)
            try:
                subprocess.run(
                    ['perl', str(rebuild_script)],
                    check=True,
                    cwd=str(Path(script_dir).resolve().parent),
                    env=rebuild_env,
                )
                print("  ✓ BLAST database rebuilt")
            except subprocess.CalledProcessError as e:
                print(f"  ⚠️  BLAST rebuild failed ({e}); mlstdb output likely sufficient")
        else:
            print(f"  ⚠️  mlst-make_blast_db not found at {rebuild_script}")
            print("      mlstdb download already produced blast/, so analysis can proceed")

        print(f"\n✅ MLST database updated. Schemes: {', '.join(installed)} in {db_dir}")
        return True

    except subprocess.CalledProcessError as e:
        print(f"❌ Update failed: {e}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)
        print(f"🧹 Removed temporary directory: {temp_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='AcinetoScope MLST Analyzer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single sample with Oxford scheme
  python mlst_module.py -i sample.fasta -o results/ -db db/ -sc bin/

  # Single sample with dual schemes (Oxford + Pasteur)
  python mlst_module.py -i sample.fasta -o results/ -db db/ -sc bin/ --dual

  # Batch processing with glob pattern
  python mlst_module.py -i "*.fasta" -o results/ -db db/ -sc bin/ --batch

  # Dual-scheme with glob pattern
  python mlst_module.py -i "*.fasta" -o results/ -db db/ -sc bin/ --dual

  # Specify Pasteur scheme
  python mlst_module.py -i sample.fasta -o results/ -db db/ -sc bin/ -s pasteur
        """
    )

    parser.add_argument('-i', '--input', required=False, default=None,
                       help='Input FASTA file, directory, or glob pattern '
                            '(required unless using --pull-db / --update-db)')
    parser.add_argument('-o', '--output-dir', required=False, default=None,
                       help='Output directory for results '
                            '(required unless using --pull-db / --update-db)')
    parser.add_argument('-db', '--database-dir', required=True,
                       help='Database directory')
    parser.add_argument('-sc', '--script-dir', required=True,
                       help='Script directory (contains mlst binary)')
    parser.add_argument('-s', '--scheme', default='oxford',
                       choices=['oxford', 'pasteur'],
                       help='MLST scheme: oxford or pasteur (default: oxford)')
    parser.add_argument('--batch', action='store_true',
                       help='Process multiple files with single scheme')
    parser.add_argument('--dual', action='store_true',
                       help='Run both Oxford and Pasteur schemes')
    parser.add_argument('--pull-db', action='store_true',
                       help='Download a fresh PubMLST A. baumannii database '
                            '(Oxford + Pasteur) and exit')
    parser.add_argument('--update-db', action='store_true',
                       help='Update the PubMLST A. baumannii database in place and exit')

    args = parser.parse_args()

    # ── Database operations (do not require -i / -o) ────────────
    if args.pull_db or args.update_db:
        try:
            success = update_mlst_database(
                db_dir=Path(args.database_dir),
                script_dir=Path(args.script_dir),
            )
        except Exception as e:
            print(f"❌ Database update failed: {e}")
            success = False
        sys.exit(0 if success else 1)

    if not args.input or not args.output_dir:
        parser.error(
            "When not using --pull-db or --update-db, -i/--input and "
            "-o/--output-dir are required."
        )

    if args.batch and args.dual:
        print("❌ Error: --batch and --dual are mutually exclusive")
        sys.exit(1)

    scheme_map = {
        'oxford': 'abaumannii',
        'pasteur': 'abaumannii_2'
    }
    database_scheme = scheme_map[args.scheme]

    try:
        analyzer = AcinetoMLSTAnalyzer(
            database_dir=Path(args.database_dir),
            script_dir=Path(args.script_dir)
        )
    except FileNotFoundError as e:
        print(f"❌ {e}")
        print("   Please ensure mlst is installed and in PATH or script directory")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"""
    {'='*80}
    🧬 ACINETOSCOPE - A. baumannii MLST Analyzer
    {'='*80}
    Author: Brown Beckley
    Email: brownbeckley94@gmail.com
    GitHub: bbeckley-hub
    {'='*80}
    Input: {args.input}
    Output: {args.output_dir}
    Scheme: {args.scheme.upper()} ({database_scheme})
    Mode: {'Dual' if args.dual else 'Batch' if args.batch else 'Single'}
    IC Database: IC1-IC12 including newly described IC10, IC11, IC12
    {'='*80}
    """)

    if args.dual:
        fasta_files = analyzer.find_fasta_files(args.input)

        if not fasta_files:
            print(f"❌ No FASTA files found for input: {args.input}")
            sys.exit(1)

        if len(fasta_files) == 1:
            input_file = fasta_files[0]
            results = analyzer.run_mlst_dual_scheme(input_file, output_dir)
            print(f"\n{'='*80}")
            print(f"🎉 Dual-scheme MLST completed for {input_file.name}")
            print(f"📁 Results saved in: {output_dir}/mlst_results/")
            print(f"{'='*80}")
        else:
            results = analyzer.run_mlst_dual_batch(args.input, output_dir)
            print(f"\n{'='*80}")
            print(f"🎉 Dual-scheme MLST completed for {len(fasta_files)} samples")
            print(f"📁 Results saved in: {output_dir}/Dual_Scheme_MLST/")
            print(f"{'='*80}")

    elif args.batch:
        try:
            results = analyzer.run_mlst_batch(args.input, output_dir, database_scheme)
            print(f"\n{'='*80}")
            print(f"🎉 Batch MLST completed! Processed {len(results)} samples")
            print(f"📁 Results saved in: {output_dir}/")
            print(f"{'='*80}")
        except Exception as e:
            print(f"\n❌ Batch processing failed: {e}")
            sys.exit(1)
    else:
        input_file = Path(args.input)
        if input_file.exists() and input_file.is_file():
            # run_mlst_single() itself appends input_file.stem — do NOT pre-nest here.
            result = analyzer.run_mlst_single(input_file, output_dir, database_scheme)

            print(f"\n{'='*80}")
            print(f"🎉 MLST completed for {input_file.name}: ST{result.get('st', 'UNKNOWN')}")
            print(f"📁 Results saved in: {sample_output_dir}/")
            print(f"{'='*80}")
        else:
            fasta_files = analyzer.find_fasta_files(args.input)
            if len(fasta_files) == 1:
                input_file = fasta_files[0]
                result = analyzer.run_mlst_single(input_file, output_dir, database_scheme)

                print(f"\n{'='*80}")
                print(f"🎉 MLST completed for {input_file.name}: ST{result.get('st', 'UNKNOWN')}")
                print(f"📁 Results saved in: {sample_output_dir}/")
                print(f"{'='*80}")
            elif len(fasta_files) > 1:
                print(f"❌ Found {len(fasta_files)} files for input: {args.input}")
                print(f"   Use --batch for multiple files or --dual for dual-scheme analysis")
                for i, f in enumerate(fasta_files[:5], 1):
                    print(f"   {i}. {f}")
                if len(fasta_files) > 5:
                    print(f"   ... and {len(fasta_files) - 5} more")
                sys.exit(1)
            else:
                print(f"❌ Input file not found: {args.input}")
                print(f"   Did you mean to use --batch for multiple files or --dual for dual-scheme?")
                sys.exit(1)


if __name__ == "__main__":
    main()