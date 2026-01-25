
# 🦠 AcinetoScope - The Ultimate A. baumannii Genomic Analysis Pipeline

```bash
 █████╗  ██████╗██╗███╗   ██╗███████╗████████╗ ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██║████╗  ██║██╔════╝╚══██╔══╝██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████║██║     ██║██╔██╗ ██║█████╗     ██║   ██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔══██║██║     ██║██║╚██╗██║██╔══╝     ██║   ██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║  ██║╚██████╗██║██║ ╚████║███████╗   ██║   ╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝  ╚═╝ ╚═════╝╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
                                                                                               
 A C I N E T O B A C T E R   B A U M A N N I I   G E N O M I C   S U R V E I L L A N C E   T O O L
                              S C I E N C E · P R E C I S I O N · I M P A C T
```
### **A species-specific computational pipeline for rapid, comprehensive *Acinetobacter baumannii* outbreak investigation and resistance gene tracking**

**Complete genomic surveillance in a single automated workflow — from FASTA to actionable insights.**

📖 Documentation • ⚡ Quick Start • ✨ Features • 🔧 Installation • 🚀 Usage • 📊 Output • 🤖 AI-Enhanced Analysis • 📈 Performance

* * *
<div align="center">

# 🦠 AcinetoScope
![Conda Version](https://anaconda.org/bbeckley-hub/acinetoscope/badges/version.svg),
![Conda Downloads](https://anaconda.org/bbeckley-hub/acinetoscope/badges/downloads.svg)
![Platform](https://anaconda.org/bbeckley-hub/acinetoscope/badges/platforms.svg)
![License](https://img.shields.io/github/license/bbeckley-hub/acinetoscope)
![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)
![Bioconda](https://img.shields.io/badge/bioconda-✓-44A833.svg?logo=conda)
![Last updated](https://anaconda.org/bbeckley-hub/acinetoscope/badges/latest_release_date.svg)

![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/acinetoscope)
![GitHub Issues](https://img.shields.io/github/issues/bbeckley-hub/acinetoscope)
![Research Use](https://img.shields.io/badge/research%20use-✓-purple)

</div>



## 📋 **Table of Contents**
*   [🎯 Overview](#-overview)
*   [✨ Key Features](#-key-features)
*   [⚡ Quick Start](#-quick-start)
*   [🔧 Installation](#-installation)
*   [🚀 Usage Guide](#-usage-guide)
*   [📊 Output Structure](#-output-structure)
*   [🔍 Analytical Modules](#-analytical-modules)
*   [🤖 AI-Enhanced Analysis](#-ai-enhanced-analysis)
*   [📈 Performance & Validation](#-performance--validation)
*   [📚 Citation](#-citation)
*   [👥 Authors & Contact](#-authors--contact)
*   [📄 License](#-license)

* * *

## 🎯 **Overview**

**AcinetoScope** is an automated, comprehensive bioinformatics pipeline designed specifically for the genomic analysis of *Acinetobacter baumannii*, a critical multidrug-resistant nosocomial pathogen. It integrates fragmented analysis steps—quality control, dual-scheme MLST typing, capsule (K/O) typing, antimicrobial resistance (AMR) detection, virulence profiling, and environmental co-selection marker screening—into a single, cohesive workflow.

### 🌍 **The Problem**
*   **Fragmented Workflows**: Analyzing *A. baumannii* requires manual chaining of 6+ separate tools (MLST, Kaptive, AMRFinder, ABRicate, etc.).
*   **Interpretation Barrier**: Raw outputs from multiple tools need manual integration to form an epidemiological narrative.
*   **Time-Consuming Process**: Generalist pipelines like Bactopia perform unnecessary steps, slowing down outbreak response.

### 💡 **Our Solution**
AcinetoScope delivers:
*   **✅ End-to-End Automation**: One command runs the entire analysis from raw FASTA to a consolidated report.
*   **✅ *A. baumannii*-Optimized**: Pre-configured with species-specific databases and typing schemes (Pasteur & Oxford MLST, Kaptive K/O loci).
*   **✅ Actionable Intelligence**: Features a four-tier risk flagging system (CRITICAL, HIGH, MEDIUM, LOW) and gene-centric tracking to highlight high-threat resistance determinants.
*   **✅ Speed & Efficiency**: Benchmarked **40-75% faster** than generalist pipelines by eliminating redundant processing.
*   **✅ AI-Ready Outputs**: Generates interactive HTML reports designed for seamless exploration with modern AI browser extensions.

**Perfect for**: Hospital outbreak investigation, public health surveillance, antimicrobial resistance (AMR) research, and clinical microbiology.

* * *

## ✨ **Key Features**

### 🔬 **Core Analytical Modules**
| Module | 🎯 Purpose | 📊 Key Outputs | ⚡ Speed |
| :--- | :--- | :--- | :--- |
| **Quality Control** | Assembly metric assessment & integrity checking | N50/N75, GC%, ambiguous bases, homopolymers | <1 min |
| **Dual MLST Typing** | Phylogenetic classification via Pasteur & Oxford schemes | Sequence Type (ST), International Clone (IC), novel alleles | <1 min |
| **Capsule (K/O) Typing** | Polysaccharide capsule & lipooligosaccharide typing via Kaptive | K type, O type, locus coverage/identity | 1-2 min |
| **AMR Detection** | Comprehensive resistance gene detection with AMRFinderPlus | Carbapenemases, ESBLs, colistin/tigecycline resistance, 4-tier risk flags | 2-3 min |
| **Multi-DB Screening** | Screening across 11 curated databases via ABRicate | Virulence factors, plasmid replicons, metal/biocide resistance, stress regulators | 3-4 min |
| **Integrated Reporting** | Synthesizes all results into gene-centric, interactive reports | HTML dashboard, JSON/CSV/TSV exports, pattern discovery | Instant |

### 🚨 **Innovations for *A. baumannii* Surveillance**
*   **Four-Tier Risk Flagging**: Automatically categorizes resistance genes (e.g., OXA-23 → CRITICAL; *qacE* → ENVIRONMENTAL).
*   **Environmental Co-Selection Tracking**: Uniquely screens for heavy metal (*czc, mer, ars*) and biocide (*qac*) resistance genes.
*   **Gene-Centric Analysis Framework**: Tracks each resistance gene across all samples for clear visualization of dissemination patterns.
*   **Cross-Genome Pattern Discovery**: Automatically identifies high-risk combinations (e.g., carbapenemase + last-resort resistance).
*   **Dynamic Resource Allocation**: Uses Python's `psutil` to optimize parallel processing for any system (from laptops to HPC clusters).

* * *

## ⚡ **Quick Start**

### **Install in 60 Seconds**
```bash
# Method 1: Conda (Recommended - handles all dependencies)
conda create -n acinetoscope -c conda-forge -c bioconda acinetoscope -y
conda activate acinetoscope

# Method 2: From source
git clone https://github.com/bbeckley-hub/acinetoscope.git
cd acinetoscope
pip install -e .
```

### **Run Your First Analysis**
```bash
# Analyze a single genome
acinetoscope -i sample.fasta -o results/

# Batch process multiple genomes
acinetoscope -i "*.fasta" -o batch_results --threads 8

# Analysis complete! Explore the interactive report.
# The main report is in: batch_results/GENIUS_ACINETOBACTER_ULTIMATE_REPORTS/
```

* * *

## 🔧 **Installation**

### **System Requirements**
| Resource | Minimum | Recommended |
| :--- | :--- | :--- |
| **CPU Cores** | 2 | 4+ |
| **RAM** | 4 GB | 8 GB |
| **Storage** | 2 GB | 10 GB+ |
| **OS** | Linux, macOS, WSL2 | Linux |

### **Step-by-Step Installation**
1.  **Install Miniconda** (if needed):
    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    # Follow prompts, then:
    source ~/.bashrc
    ```
2.  **Install AcinetoScope**:
    ```bash
    conda create -n acinetoscope -c conda-forge -c bioconda acinetoscope
    conda activate acinetoscope
    ```
3.  **(Recommended) Update ABRicate Databases**:
    ```bash
    abricate --setupdb
    ```

* * *

## 🚀 **Usage Guide**

### **Basic Command**
```bash
acinetoscope -i <INPUT_PATTERN> -o <OUTPUT_DIR> [OPTIONS]
```
**Example**: `acinetoscope -i "genomes/*.fna" -o my_analysis -t 4`

### **Command Line Options**
| Flag | Description | Default |
| :--- | :--- | :--- |
| `-i, --input` | Input FASTA file(s). Supports wildcards (`*.fna`). | **Required** |
| `-o, --output` | Directory for all results. | **Required** |
| `-t, --threads` | Number of CPU threads to use. | Auto-detected |
| `--skip-qc` | Skip the quality control module. | False |
| `--skip-summary` | Skip the final integrated report generation. | False |
| `--mlst-scheme` | Specify scheme: `pasteur`, `oxford`, or `both`. | `both` |
| `--verbose` | Print detailed progress messages. | False |

### **Input Requirements**
*   **Format**: Assembled genomes in FASTA format (`.fna`, `.fasta`, `.fa`, `.fn`).
*   **Content**: The pipeline is designed exclusively for *Acinetobacter baumannii* genomes.

* * *

## 📊 **Output Structure**
AcinetoScope generates a well-organized directory structure. A typical run on 10 genomes produces over 500 files, systematically organized as shown in your `tree` output. The key directories are:

```
analysis/
├── fasta_qc_results/                 # Quality control reports per sample
├── PASTEUR_MLST/                     # MLST results (Pasteur scheme)
├── OXFORD_MLST/                      # MLST results (Oxford scheme)
├── kaptive_results/                  # Capsule (K) and lipooligosaccharide (O) typing
├── acineto_amrfinder_results/        # AMR gene detection with risk stratification
├── acineto_abricate_results/         # Multi-database screening (11 DBs)
└── GENIUS_ACINETOBACTER_ULTIMATE_REPORTS/  # 🎯 FINAL INTEGRATED REPORT
    ├── genius_acinetobacter_ultimate_report.html  # Interactive HTML Dashboard
    ├── genius_acinetobacter_ultimate_report.json  # Complete data (machine-readable)
    └── *.csv files for easy import into spreadsheets
```

* * *

## 🔍 **Analytical Modules**

### **1. Quality Control Module**
Validates input FASTA files using *A. baumannii*-specific thresholds (GC% 35-65%, ambiguous bases <5%). Generates per-sample and summary reports.

### **2. Dual MLST Typing Module**
Runs the `mlst` tool with both **Pasteur** and **Oxford** schemes. Identifies International Clones (IC1-IC10) and novel alleles.

### **3. Capsule Typing Module**
Utilizes **Kaptive** with *A. baumannii*-specific databases (`ab_k`, `ab_o`) for high-resolution capsule (K) and lipooligosaccharide (O) typing.

### **4. AMR Detection Module**
Leverages **NCBI's AMRFinderPlus** to identify and categorize resistance determinants, applying the four-tier risk flagging system (CRITICAL, HIGH, MEDIUM, LOW).

### **5. Comprehensive Screening Module**
Executes **ABRicate** across 11 databases (CARD, ResFinder, VFDB, PlasmidFinder, BacMet2, etc.) to detect virulence factors, plasmid replicons, and environmental co-selection markers.

### **6. Integrated Reporting Module**
The **GENIUS Acinetobacter Reporter** synthesizes all results into a **gene-centric interactive HTML report**, enabling cross-sample pattern discovery and high-risk combination identification.

* * *

## 🤖 **AI-Enhanced Analysis**

AcinetoScope's interactive HTML reports are **designed for AI augmentation**. You can use browser AI extensions (ChatGPT, Claude, Gemini, Copilot) to interrogate your genomic data conversationally.

### **🎯 Guiding Principle: AI Assists, Experts Decide!**
Use AI as a collaborative tool to explore data and generate hypotheses, but always apply your domain expertise for final interpretation and clinical decisions.

### **🚀 Step-by-Step AI Integration**

1.  **Generate Your Report**: Run AcinetoScope to create the `genius_acinetobacter_ultimate_report.html`.
2.  **Install an AI Assistant**: Add a browser extension like **ChatGPT for Chrome**, **Claude**, or **Microsoft Copilot** to your browser.
3.  **Open and Explore**:
    *   Open the HTML report in your browser.
    *   Use the AI extension's "ask about this page" feature or simply copy-paste findings into the chat.
4.  **Ask Powerful Questions**:
    *   **"Summarize the primary resistance threat in these isolates."**
    *   **"Is there evidence of an outbreak cluster based on the ST and capsule types?"**
    *   **"Which samples carry both a carbapenemase and a colistin resistance mechanism? List them."**
    *   **"Generate a concise clinical risk assessment for infection control."**
    *   **"Suggest antibiotic treatment options based on this resistance profile."**

### **💡 Example AI Interaction**
**Your Prompt**: "I'm looking at the AcinetoScope report for 10 ICU isolates. The summary says 8 are ST2 and carry OXA-23. What's the immediate implication?"

**AI Assistant Response**: "This suggests a likely clonal outbreak of a high-risk carbapenem-resistant *A. baumannii* (CRAB) strain in your ICU. Immediate actions should include: 1) Reviewing infection control practices, 2) Patient cohorting, 3) Environmental decontamination focus. The co-presence of [other genes from report] indicates limited treatment options, necessitating an infectious disease consult."

* * *

## 📈 **Performance & Validation**

### **⚡ Benchmarking vs. Bactopia**
AcinetoScope is purpose-built for *A. baumannii*, making it significantly faster than generalist pipelines.

| System Config | Pipeline | Time (50 genomes) | Speed Gain |
| :--- | :--- | :--- | :--- |
| 2 CPU, 8 GB RAM | **AcinetoScope** | **~2.5 hours** | **≈40% faster** |
| | Bactopia | ~4 hours | |
| 16 CPU, 16 GB RAM | **AcinetoScope** | **~35 minutes** | **≈75% faster** |
| | Bactopia | ~2.5 hours | |

### **✅ Validation Results**
Tested on 10 well-characterized reference genomes, AcinetoScope achieved **100% accuracy** in:
*   MLST typing (Pasteur & Oxford schemes)
*   Capsule (K/O) type determination
*   Identification of known antimicrobial resistance genes

### **🔬 Key Findings from Clinical Genomes**
Analysis of 50 clinical genomes revealed:
*   **High-Risk Clones**: Dominance of International Clone II (ST2, 46%) and I (ST1, 26%).
*   **Critical Resistance**: 56% harbored carbapenemase genes (*bla_OXA-23/66*); 96% co-harbored carbapenemase + last-resort resistance genes.
*   **Environmental Persistence**: 100% contained heavy metal resistance genes; 58% had the biocide resistance gene *qacEdelta1*.

* * *

## 📚 **Citation**

If you use AcinetoScope in your research, please cite:

**Beckley, B. et al.** (2026). AcinetoScope: A Tool for Enhanced Outbreak Investigation and Resistance Gene Tracking in *Acinetobacter baumannii*. *[Journal Name],* [Volume], [Pages]. DOI: [To be assigned].

```bibtex
@software{acinetoscope2026,
  title = {AcinetoScope: A Tool for Enhanced Outbreak Investigation and Resistance Gene Tracking in Acinetobacter baumannii},
  author = {Beckley, Brown and Amarh, Vincent and Lopes, Bruno Silvester and Kakah, John and Kwarteng, Alexander and Olalekan, Adesola and Afeke, Innocent},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/bbeckley-hub/acinetoscope}
}
```

* * *

## 👥 **Authors & Contact**

*   **Brown Beckley** (Corresponding Author) – Department of Medical Biochemistry, University of Ghana Medical School. Email: `brownbeckley94@gmail.com`
*   Vincent Amarh – University of Ghana Medical School
*   Bruno Silvester Lopes – Teesside University, UK
*   John Kakah – University of Ghana Medical School
*   Alexander Kwarteng – Kwame Nkrumah University of Science and Technology (KNUST)
*   Adesola Olalekan – University of Lagos
*   Innocent Afeke – University of Health and Allied Sciences

**GitHub Repository**: [https://github.com/bbeckley-hub/acinetoscope](https://github.com/bbeckley-hub/acinetoscope)

* * *

## 📄 **License**

This project is licensed under the **MIT License**. See the LICENSE file in the repository for details.

---
**⭐ Star us on GitHub if you find AcinetoScope useful!**

**Transforming complex genomic data into clear, actionable insights for tackling AMR.** 🧬✨
<div>
