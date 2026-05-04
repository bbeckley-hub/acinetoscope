<p align="center">
  <img src="https://raw.githubusercontent.com/bbeckley-hub/acinetoscope/main/acinetoscope_banner.png" alt="AcinetoScope Banner" width="100%">
</p>

<div align="center">

### **A species-specific computational pipeline for rapid, comprehensive _Acinetobacter baumannii_ outbreak investigation and resistance gene tracking**

**Complete genomic surveillance in a single automated workflow — from FASTA to actionable insights.**

[![Conda Version](https://anaconda.org/bbeckley-hub/acinetoscope/badges/version.svg)](https://anaconda.org/bbeckley-hub/acinetoscope)
[![Platform](https://anaconda.org/bbeckley-hub/acinetoscope/badges/platforms.svg)](https://anaconda.org/bbeckley-hub/acinetoscope)
[![License](https://img.shields.io/github/license/bbeckley-hub/acinetoscope)](LICENSE)
[![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![Bioconda](https://img.shields.io/badge/bioconda-✓-44A833.svg?logo=conda)](https://anaconda.org/bioconda/acinetoscope)
[![Last updated](https://anaconda.org/bbeckley-hub/acinetoscope/badges/latest_release_date.svg)](https://anaconda.org/bbeckley-hub/acinetoscope)
[![Docker Pulls](https://img.shields.io/docker/pulls/bbeckleyhub/acinetoscope)](https://hub.docker.com/r/bbeckleyhub/acinetoscope)

[![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/acinetoscope)](https://github.com/bbeckley-hub/acinetoscope/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/bbeckley-hub/acinetoscope)](https://github.com/bbeckley-hub/acinetoscope/issues)
[![Research Use](https://img.shields.io/badge/research%20use-✓-purple)]()
![Latest Release Date](https://anaconda.org/bbeckley-hub/acinetoscope/badges/latest_release_date.svg)
![License](https://anaconda.org/bbeckley-hub/acinetoscope/badges/license.svg)


[![Docker Image Size](https://img.shields.io/docker/image-size/bbeckleyhub/acinetoscope/latest)](https://hub.docker.com/r/bbeckleyhub/acinetoscope)
[![Docker Version](https://img.shields.io/docker/v/bbeckleyhub/acinetoscope?sort=semver)](https://hub.docker.com/r/bbeckleyhub/acinetoscope)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](#)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Profile-0A66C2?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/brown-beckley-190315319)
[![Stage](https://img.shields.io/badge/status-active-brightgreen)](#)


[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://docs.conda.io/en/latest/)
[![GitHub Issues](https://img.shields.io/github/issues/bbeckley-hub/acinetoscope)](https://github.com/bbeckley-hub/acinetoscope/issues)
[![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/acinetoscope)](https://github.com/bbeckley-hub/acinetoscope/stargazers)
[![Sample Report](https://img.shields.io/badge/📊-View_Sample_Report-blue)](https://htmlpreview.github.io/?https://bbeckley-hub.github.io/acinetoscope/#summary)
![Profile Views](https://komarev.com/ghpvc/?username=bbeckley-hub&label=Profile%20Views&color=0e75b6&style=flat)
[![Google Scholar](https://img.shields.io/badge/Google%20Scholar-Profile-4285F4?style=flat&logo=googlescholar&logoColor=white)](https://scholar.google.com/citations?user=CYNOsqIAAAAJ&hl=en)


![GitHub stats](https://github-readme-stats.vercel.app/api?username=bbeckley-hub&show_icons=true&theme=radical)
![Top Langs](https://github-readme-stats.vercel.app/api/top-langs/?username=bbeckley-hub&layout=compact&theme=radical)
[![GitHub Streak](https://streak-stats.demolab.com?user=bbeckley-hub&theme=radical&date_format=j%20M%5B%20Y%5D)](https://git.io/streak-stats)

</div>

---

## 📋 **Table of Contents**

- [🎯 Overview](#-overview)
- [✨ Key Features](#-key-features)
- [⚡ Quick Start](#-quick-start)
- [🔧 Installation](#-installation)
- [🐳 Docker Usage](#-docker-usage)
- [🚀 Usage Guide](#-usage-guide)
- [📊 Output Structure](#-output-structure)
- [🔍 Analytical Modules](#-analytical-modules)
- [🔗 Integrated External Tools & Dependencies](#-integrated-external-tools--dependencies)
- [🤖 AI-Enhanced Analysis](#-ai-enhanced-analysis)
- [📈 Performance & Validation](#-performance--validation)
- [📚 Citation](#-citation)
- [📚 Third-Party Tool Citations](#-third-party-tool-citations)
- [👥 Authors & Contact](#-authors--contact)
- [📄 License](#-license)

---

## 🎯 **Overview**

**AcinetoScope** is an automated, comprehensive bioinformatics pipeline designed specifically for the genomic analysis of *Acinetobacter baumannii*, a critical multidrug-resistant nosocomial pathogen. It integrates fragmented analysis steps—quality control, dual-scheme MLST typing, capsule (K/O) typing, antimicrobial resistance (AMR) detection, virulence profiling, and environmental co-selection marker screening—into a single, cohesive workflow.

### 🌍 **The Problem**
- **Fragmented Workflows**: Analyzing *A. baumannii* requires manual chaining of 6+ separate tools (MLST, Kaptive, AMRFinder, ABRicate, etc.).
- **Interpretation Barrier**: Raw outputs from multiple tools need manual integration to form an epidemiological narrative.
- **Time-Consuming Process**: Generalist pipelines like Bactopia perform unnecessary steps, slowing down outbreak response.

### 💡 **Our Solution**
AcinetoScope delivers:
- **✅ End-to-End Automation**: One command runs the entire analysis from raw FASTA to a consolidated report.
- **✅ *A. baumannii*-Optimized**: Pre-configured with species-specific databases and typing schemes (Pasteur & Oxford MLST, Kaptive K/O loci).
- **✅ Actionable Intelligence**: Features a four-tier risk flagging system (CRITICAL, HIGH, MEDIUM, LOW) and gene-centric tracking to highlight high-threat resistance determinants.
- **✅ Speed & Efficiency**: Benchmarked **40-75% faster** than generalist pipelines by eliminating redundant processing.
- **✅ AI-Ready Outputs**: Generates interactive HTML reports designed for seamless exploration with modern AI browser extensions.

**Perfect for**: Hospital outbreak investigation, public health surveillance, antimicrobial resistance (AMR) research, and clinical microbiology.

---

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
- **Four-Tier Risk Flagging**: Automatically categorizes resistance genes (e.g., OXA-23 → CRITICAL; *qacE* → ENVIRONMENTAL).
- **Environmental Co-Selection Tracking**: Uniquely screens for heavy metal (*czc, mer, ars*) and biocide (*qac*) resistance genes.
- **Gene-Centric Analysis Framework**: Tracks each resistance gene across all samples for clear visualization of dissemination patterns.
- **Cross-Genome Pattern Discovery**: Automatically identifies high-risk combinations (e.g., carbapenemase + last-resort resistance).
- **Dynamic Resource Allocation**: Uses Python's `psutil` to optimize parallel processing for any system (from laptops to HPC clusters).


## 📊 Sample Output

See a complete interactive report generated by AcinetoScope:

[![Sample Report](https://img.shields.io/badge/📊-View_Sample_Report-blue)](https://htmlpreview.github.io/?https://bbeckley-hub.github.io/acinetoscope/#summary)

*The report includes AMR and virulence gene tables, filter buttons, combination tables, and FASTA QC metrics.*

---

## ⚡ **Quick Start**

### **Install in 60 Seconds**
```bash
# Method 1: Conda (Recommended)
conda create -n acinetoscope -c conda-forge -c bioconda  -c bbeckley-hub acinetoscope -y
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
# Main report: batch_results/GENIUS_ACINETOBACTER_ULTIMATE_REPORTS/
```

---

## 🔧 **Installation**

### **System Requirements**
| Resource | Minimum | Recommended |
| :--- | :--- | :--- |
| **CPU Cores** | 2 | 4+ |
| **RAM** | 4 GB | 8 GB |
| **Storage** | 2 GB | 10 GB+ |
| **OS** | Linux, macOS, WSL2 | Linux |

### **Step-by-Step Installation**
1. **Install Miniconda** (if needed):
    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    source ~/.bashrc
    ```
2. **Install AcinetoScope**:
    ```bash
    conda create -n acinetoscope -c conda-forge -c bioconda -c bbeckley-hub acinetoscope -y
    conda activate acinetoscope
    ```
3. **(Recommended) Update ABRicate Databases**:
    ```bash
    abricate --setupdb
    ```

---

## 🐳 **Docker Usage**

For users who prefer a containerized environment or cannot install Conda, we provide a Docker image with all dependencies pre‑installed and ABRicate databases pre‑configured. Run the complete *Acinetobacter baumannii* typing pipeline with zero installation – just Docker.

---

## 📦 What’s inside this Docker image

- Full **AcinetoScope** pipeline (all modules)
- All dependencies pre‑installed (Conda environment, Perl, BLAST, ABRicate, Kaptive, etc.)
- **ABRicate databases** pre‑configured (`abricate --setupdb` already run)
- **jq** installed for reliable JSON parsing
- No need for Conda, no manual setup, no “read‑only filesystem” errors

---

## 🚀 Quick Start

### Pull the image

```bash
docker pull bbeckleyhub/acinetoscope:latest
```

### Run on a single FASTA file

```bash
docker run --rm -v $(pwd):/data bbeckleyhub/acinetoscope:latest -i "/data/genome.fna" -o /data/output
```

After the run, output files are owned by `root` on your host. To reclaim ownership:

```bash
sudo chown -R $USER:$USER ./output
```

### Run on all FASTA files in the current directory

```bash
docker run --rm -v $(pwd):/data bbeckleyhub/acinetoscope:latest -i "/data/*.fna" -o /data/output
```

---

## 📖 Detailed Usage

### Basic syntax

```bash
docker run --rm -v $(pwd):/data bbeckleyhub/acinetoscope:latest [OPTIONS]
```

- `--rm` : remove container after exit
- `-v $(pwd):/data` : mount current directory to `/data` inside container
- Input files must be under `/data` (e.g., `/data/*.fna`)
- Output directory must also be under `/data` (e.g., `/data/output`)

### All AcinetoScope options work

```bash
docker run --rm -v $(pwd):/data bbeckleyhub/acinetoscope:latest \
  -i "/data/*.fna" -o /data/output \
  --threads 8 --skip-qc --skip-amr
```

See `docker run --rm bbeckleyhub/acinetoscope:latest -h` for all options.

### Using custom threads

```bash
docker run --rm -v $(pwd):/data bbeckleyhub/acinetoscope:latest \
  -i "/data/*.fna" -o /data/output -t 16
```

---

## 🔧 Handling File Permissions (The “Padlock” Issue)

By default, Docker runs as `root` inside the container. Any files written to your mounted directory will be owned by `root:root`.  
You have three options:

### 1. Change ownership after the run (easiest)

```bash
sudo chown -R $USER:$USER ./output
```

### 2. Run with your host user ID (requires a small code fix – coming soon)

Currently not fully supported because AcinetoScope needs to write to its own installation directory. A future update will fix this.

### 3. Use Singularity (recommended for HPC, no `sudo` needed)

See the [Singularity section](#singularity-for-hpc-no-sudo) below.

---

## 🧪 Testing Your Docker Setup

### Check help message

```bash
docker run --rm bbeckleyhub/acinetoscope:latest -h
```

### Verify ABRicate databases are installed

```bash
docker run --rm --entrypoint /bin/bash bbeckleyhub/acinetoscope:latest -c "abricate --list | head -5"
```

Expected output: list of databases (ncbi, card, vfdb, etc.)

### Verify jq is installed (important for correct summaries)

```bash
docker run --rm --entrypoint /bin/bash bbeckleyhub/acinetoscope:latest -c "jq --version"
```

Should output `jq-1.6` or similar.

---

## 🖥️ Singularity for HPC (no `sudo`, correct ownership)

On HPC clusters that support [Singularity/Apptainer](https://sylabs.io/singularity/), you can run AcinetoScope **without `sudo`** and output files will be owned by your user automatically.

> **Important:** AcinetoScope writes temporary files inside its own installation directory (`/opt/acinetoscope/...`). Singularity mounts containers as read‑only by default, so you **must** add the `--writable-tmpfs` flag to allow these writes. The flag creates an ephemeral, writable overlay in memory – no permanent changes are made to the container.

### Option A: Direct pull (if network allows)

```bash
singularity pull acinetoscope.sif docker://bbeckleyhub/acinetoscope:latest
singularity run --writable-tmpfs -B $(pwd):/data acinetoscope.sif -i "/data/*.fna" -o /data/output
```

### Option B: Convert from a local Docker image (when `singularity pull` fails)

If you encounter TLS timeouts or other network errors (common on some HPCs), convert an existing Docker image to a Singularity SIF file on a machine with Docker, then transfer the `.sif` file to the HPC.

**Step 1 – on a machine with Docker (e.g., your laptop):**

```bash
docker pull bbeckleyhub/acinetoscope:latest
docker save bbeckleyhub/acinetoscope:latest -o acinetoscope.tar
singularity build acinetoscope.sif docker-archive://acinetoscope.tar
```

Now copy `acinetoscope.sif` to your HPC home or project directory (e.g., using `scp`).

**Step 2 – on the HPC (no sudo needed):**

```bash
singularity run --writable-tmpfs -B $(pwd):/data acinetoscope.sif -i "/data/*.fna" -o /data/output
```

### Explanation of flags

| Flag | Purpose |
|------|---------|
| `--writable-tmpfs` | Creates a temporary writable overlay – **required** for AcinetoScope to write intermediate files to `/opt/...` |
| `-B $(pwd):/data` | Binds your current directory to `/data` inside the container (input files are read from here, output is written here) |
| `-i "/data/*.fna"` | Input pattern – use quotes to prevent shell expansion on the host |
| `-o /data/output` | Output directory (will appear as `./output` on your host) |

### Additional options

You can use any AcinetoScope flag, e.g.:

```bash
singularity run --writable-tmpfs -B $(pwd):/data acinetoscope.sif \
    -i "/data/*.fna" -o /data/output --threads 8 --skip-qc
```

### Verify it works

After a successful run, you will see output like:

```
✓ QC analysis completed!
✓ MLST Pasteur completed
...
✓ 🎉 All analyses completed successfully!
```

All result files in `./output` will be owned by **your HPC user** – no `sudo chown` needed.

---

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

### **Input Requirements**
- **Format**: Assembled genomes in FASTA format (`.fna`, `.fasta`, `.fa`, `.fn`).
- **Content**: Designed exclusively for *Acinetobacter baumannii* genomes.

---

## 📊 **Output Structure**

AcinetoScope generates a well-organized directory. Key directories include:

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

---

## 🔍 **Analytical Modules**

1. **Quality Control Module**: Validates input using *A. baumannii*-specific thresholds (GC% 35-65%, ambiguous bases <5%).
2. **Dual MLST Typing Module**: Runs `mlst` with both **Pasteur** and **Oxford** schemes. Identifies International Clones (IC1-IC10).
3. **Capsule Typing Module**: Uses **Kaptive** with *A. baumannii*-specific databases (`ab_k`, `ab_o`).
4. **AMR Detection Module**: Leverages **NCBI's AMRFinderPlus** with a four-tier risk flagging system.
5. **Comprehensive Screening Module**: Executes **ABRicate** across 11 databases (CARD, ResFinder, VFDB, PlasmidFinder, BacMet2, etc.).
6. **Integrated Reporting Module**: The **GENIUS Acinetobacter Reporter** synthesizes all results into a gene-centric interactive HTML report.

---

## 🔗 **Integrated External Tools & Dependencies**

AcinetoScope integrates several powerful open-source tools and databases. These are **not bundled directly in this repository**. Instead, they are automatically installed as **dependencies via Conda** (as defined in the Conda recipe). The MIT license applies only to the AcinetoScope pipeline code (the workflow engine, report generation, and Python modules written by the authors). Each tool is used under the terms of its own license, and we gratefully acknowledge their authors.

| Tool/Database | Purpose | Source | License |
|---------------|---------|--------|---------|
| **MLST** | Multi-locus sequence typing | [tseemann/mlst](https://github.com/tseemann/mlst) | GPL v2 |
| **ABRicate** | Mass screening for resistance/virulence | [tseemann/abricate](https://github.com/tseemann/abricate) | GPL v2 |
| **AMRFinderPlus** | AMR gene detection | [ncbi/amr](https://github.com/ncbi/amr) | Public Domain |
| **Kaptive** | Capsule (K/O) typing | [katholt/Kaptive](https://github.com/katholt/Kaptive) | GPL v3 |
| **Pasteur MLST DB** | *A. baumannii* MLST scheme | [PubMLST](https://pubmlst.org/) | Free for research |
| **Oxford MLST DB** | *A. baumannii* MLST scheme | [PubMLST](https://pubmlst.org/) | Free for research |
| **Kaptive DB** | K/O locus databases | [katholt/Kaptive](https://github.com/katholt/Kaptive) | GPL v3 |
| **CARD** | AMR database | [card.mcmaster.ca](https://card.mcmaster.ca/) | ODbL |
| **ResFinder** | Acquired AMR genes | [cge.cbs.dtu.dk](https://cge.cbs.dtu.dk/services/ResFinder/) | Free for research |
| **VFDB** | Virulence factors | [mgc.ac.cn/VFs](http://www.mgc.ac.cn/VFs/) | Free for research |
| **PlasmidFinder** | Plasmid replicons | [cge.cbs.dtu.dk](https://cge.cbs.dtu.dk/services/PlasmidFinder/) | Free for research |
| **BacMet2** | Biocide/metal resistance | [bacmet.biomedicine.gu.se](http://bacmet.biomedicine.gu.se/) | Free for research |

By using AcinetoScope, you agree to comply with the licenses of these third-party tools and databases.

---

## 🤖 **AI-Enhanced Analysis**

AcinetoScope's interactive HTML reports are **designed for AI augmentation**. You can use browser AI extensions (ChatGPT, Claude, Gemini, Copilot) to interrogate your genomic data conversationally.

### 🎯 **Guiding Principle: AI Assists, Experts Decide!**
Use AI as a collaborative tool to explore data and generate hypotheses, but always apply your domain expertise for final interpretation and clinical decisions.

### 🚀 **Step-by-Step AI Integration**
1. **Generate Your Report**: Run AcinetoScope to create `genius_acinetobacter_ultimate_report.html`.
2. **Install an AI Assistant**: Add a browser extension like **ChatGPT for Chrome**, **Claude**, or **Microsoft Copilot**.
3. **Open and Explore**:
   - Open the HTML report in your browser.
   - Use the AI extension's "ask about this page" feature or copy-paste findings into the chat.
4. **Ask Powerful Questions**:
   - *"Summarize the primary resistance threat in these isolates."*
   - *"Is there evidence of an outbreak cluster based on the ST and capsule types?"*
   - *"Which samples carry both a carbapenemase and a colistin resistance mechanism? List them."*
   - *"Generate a concise clinical risk assessment for infection control."*
   - *"Suggest antibiotic treatment options based on this resistance profile."*

### 💡 **Example AI Interaction**
**Your Prompt**: "I'm looking at the AcinetoScope report for 10 ICU isolates. The summary says 8 are ST2 and carry OXA-23. What's the immediate implication?"

**AI Assistant Response**: "This suggests a likely clonal outbreak of a high-risk carbapenem-resistant *A. baumannii* (CRAB) strain in your ICU. Immediate actions should include: 1) Reviewing infection control practices, 2) Patient cohorting, 3) Environmental decontamination focus. The co-presence of [other genes from report] indicates limited treatment options, necessitating an infectious disease consult."

---

## 📈 **Performance & Validation**

### ⚡ **Benchmarking vs. Bactopia**
AcinetoScope is purpose-built for *A. baumannii*, making it significantly faster than generalist pipelines.

| System Config | Pipeline | Time (50 genomes) | Speed Gain |
| :--- | :--- | :--- | :--- |
| 2 CPU, 8 GB RAM | **AcinetoScope** | **~2.5 hours** | **≈40% faster** |
| | Bactopia | ~4 hours | |
| 16 CPU, 16 GB RAM | **AcinetoScope** | **~35 minutes** | **≈75% faster** |
| | Bactopia | ~2.5 hours | |

### ✅ **Validation Results**
Tested on 10 well-characterized reference genomes, AcinetoScope achieved **100% accuracy** in:
- MLST typing (Pasteur & Oxford schemes)
- Capsule (K/O) type determination
- Identification of known antimicrobial resistance genes

### 🔬 **Key Findings from Clinical Genomes**
Analysis of 50 clinical genomes revealed:
- **High-Risk Clones**: Dominance of International Clone II (ST2, 46%) and I (ST1, 26%).
- **Critical Resistance**: 56% harbored carbapenemase genes (*bla_OXA-23/66*); 96% co-harbored carbapenemase + last-resort resistance genes.
- **Environmental Persistence**: 100% contained heavy metal resistance genes; 58% had the biocide resistance gene *qacEdelta1*.

---

## 📚 **Citation**

If you use AcinetoScope in your research, please cite:

```bibtex
@software{acinetoscope2026,
  title = {AcinetoScope: A Tool for Enhanced Outbreak Investigation and Resistance Gene Tracking in Acinetobacter baumannii},
  author = {Beckley, B. and Amarh, V. and Lopes, B. S. and Kakah, J. and Kwarteng, A. and Olalekan, A. and Afeke, I.},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/bbeckley-hub/acinetoscope}
}
```

---

## 📚 **Third-Party Tool Citations**

AcinetoScope integrates several essential third-party tools and databases. If you use AcinetoScope in your research, please also cite the following:

#### **MLST (Torsten Seemann)**
```bibtex
@software{seemann_mlst_2018,
  author = {Seemann, T.},
  title = {MLST: Scan contig files against traditional PubMLST typing schemes},
  year = {2018},
  publisher = {GitHub},
  url = {https://github.com/tseemann/mlst}
}
```

#### **ABRicate (Torsten Seemann)**
```bibtex
@software{seemann_abricate_2018,
  author = {Seemann, T.},
  title = {ABRicate: Mass screening of contigs for antimicrobial resistance and virulence genes},
  year = {2018},
  publisher = {GitHub},
  url = {https://github.com/tseemann/abricate}
}
```

#### **AMRFinderPlus (NCBI)**
```bibtex
@article{feldgarden_amrfinderplus_2021,
  author = {Feldgarden, M. et al.},
  title = {AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence},
  journal = {Scientific Reports},
  volume = {11},
  pages = {12728},
  year = {2021},
  doi = {10.1038/s41598-021-91456-0}
}
```

#### **Kaptive (Kath Holt Lab)**
```bibtex
@article{wyres_kaptive_2025,
  author = {Wyres, K. L. et al.},
  title = {Kaptive: a tool for identification of Klebsiella pneumoniae and Acinetobacter baumannii capsule loci},
  journal = {Microbial Genomics},
  volume = {6},
  number = {3},
  year = {2025},
  doi = {10.1099/mgen.0.000334}
}
```

#### **CARD Database**
```bibtex
@article{alcock_card_2023,
  author = {Alcock, B. P. et al.},
  title = {CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database},
  journal = {Nucleic Acids Research},
  volume = {51},
  number = {D1},
  pages = {D690-D699},
  year = {2023},
  doi = {10.1093/nar/gkac920}
}
```

#### **ResFinder**
```bibtex
@article{bortolaia_resfinder_2020,
  author = {Bortolaia, V. et al.},
  title = {ResFinder 4.0 for predictions of phenotypes from genotypes},
  journal = {Journal of Antimicrobial Chemotherapy},
  volume = {75},
  number = {12},
  pages = {3491-3500},
  year = {2020},
  doi = {10.1093/jac/dkaa345}
}
```

#### **VFDB**
```bibtex
@article{chen_vfdb_2016,
  author = {Chen, L. et al.},
  title = {VFDB 2016: hierarchical and refined dataset for big data analysis—10 years on},
  journal = {Nucleic Acids Research},
  volume = {44},
  number = {D1},
  pages = {D694-D697},
  year = {2016},
  doi = {10.1093/nar/gkv1239}
}
```

#### **PlasmidFinder**
```bibtex
@article{carattoli_plasmidfinder_2014,
  author = {Carattoli, A. et al.},
  title = {In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing},
  journal = {Antimicrobial Agents and Chemotherapy},
  volume = {58},
  number = {7},
  pages = {3895-3903},
  year = {2014},
  doi = {10.1128/AAC.02412-14}
}
```

#### **BacMet2**
```bibtex
@article{Pal_bacmet_2014,
  author = {Pal, C. et al.},
  title = {BacMet: antibacterial biocide and metal resistance genes database},
  journal = {Nucleic Acids Research},
  volume = {42},
  pages = {D737-D743},
  year = {2014},
  doi = {10.1093/nar/gkt1252}
}
```

---

## 👥 **Authors & Contact**

- **Brown Beckley** (Corresponding Author) – Department of Medical Biochemistry, University of Ghana Medical School. Email: `brownbeckley94@gmail.com`
- **Vincent Amarh** – University of Ghana Medical School
- **Bruno Silvester Lopes** – Teesside University, UK
- **John Kakah** – University of Ghana Medical School
- **Alexander Kwarteng** – Kwame Nkrumah University of Science and Technology (KNUST)
- **Adesola Olalekan** – University of Lagos
- **Innocent Afeke** – University of Health and Allied Sciences

**GitHub Repository**: [https://github.com/bbeckley-hub/acinetoscope](https://github.com/bbeckley-hub/acinetoscope)

---

## 📄 **License**

### Core AcinetoScope Code
The AcinetoScope pipeline code (the workflow engine, report generation, HTML templates, and Python modules written by the authors) is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.

### Third-Party Tools
AcinetoScope executes several external bioinformatics tools, which are installed as Conda dependencies. Each tool is the property of its respective developers and is used under its own license. Key dependencies include:

| Tool | License |
|------|---------|
| `mlst` (Torsten Seemann) | GPL v2 |
| `ABRicate` (Torsten Seemann) | GPL v2 |
| `AMRFinderPlus` (NCBI) | Public Domain |
| `Kaptive` (Kath Holt) | GPL v3 |
| Various ABRicate databases | Various (see above) |

By using AcinetoScope, you agree to comply with the licenses of these third-party tools and databases.

---

<div align="center">

## **⭐ Star us on GitHub if you find AcinetoScope useful!**

**Transforming complex genomic data into clear, actionable insights for tackling AMR.** 🧬✨

</div>
