<p align="center">
  <img src="https://raw.githubusercontent.com/bbeckley-hub/acinetoscope/main/acinetoscope.png" alt="AcinetoScope Banner" width="100%">
</p>

<div align="center">

# 🧬 AcinetoScope

### **A species‑specific computational pipeline for rapid, comprehensive *Acinetobacter baumannii* outbreak investigation and resistance gene tracking**

#### **Complete genomic surveillance in a single automated workflow — from FASTA to actionable insights.**

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![Version](https://img.shields.io/badge/version-1.3.0-blue)](https://github.com/bbeckley-hub/acinetoscope/releases)
[![Conda Version](https://anaconda.org/bbeckley-hub/acinetoscope/badges/version.svg)](https://anaconda.org/bbeckley-hub/acinetoscope)
[![Platform](https://anaconda.org/bbeckley-hub/acinetoscope/badges/platforms.svg)](https://anaconda.org/bbeckley-hub/acinetoscope)
[![License](https://img.shields.io/github/license/bbeckley-hub/acinetoscope)](LICENSE)
[![Latest Release Date](https://anaconda.org/bbeckley-hub/acinetoscope/badges/latest_release_date.svg)](https://anaconda.org/bbeckley-hub/acinetoscope)
[![Bioconda](https://img.shields.io/badge/bioconda-✓-44A833.svg?logo=conda)](https://anaconda.org/bioconda/acinetoscope)

[![Docker Pulls](https://img.shields.io/docker/pulls/bbeckleyhub/acinetoscope)](https://hub.docker.com/r/bbeckleyhub/acinetoscope)
[![Docker Image Size](https://img.shields.io/docker/image-size/bbeckleyhub/acinetoscope/latest)](https://hub.docker.com/r/bbeckleyhub/acinetoscope)
[![Docker Version](https://img.shields.io/docker/v/bbeckleyhub/acinetoscope?sort=semver)](https://hub.docker.com/r/bbeckleyhub/acinetoscope)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](#)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Profile-0A66C2?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/brown-beckley-190315319)
[![Stage](https://img.shields.io/badge/status-active-brightgreen)](#)

[![Powered by 🧠](https://img.shields.io/badge/powered%20by-science%20🔬-purple)](https://github.com/bbeckley-hub/acinetoscope)
[![Coffee](https://img.shields.io/badge/built%20with-%E2%98%95%20coffee-orange)](https://github.com/bbeckley-hub/acinetoscope)
[![Made with ❤️](https://img.shields.io/badge/made%20with-%E2%9D%A4%EF%B8%8F-red)](https://github.com/bbeckley-hub/acinetoscope)
[![Open Source Love](https://badges.frapsoft.com/os/v1/open-source.svg?v=103)](https://github.com/ellerbrock/open-source-badges/)
[![Made for Research](https://img.shields.io/badge/made%20for-Research-0066cc.svg)](https://github.com/bbeckley-hub/acinetoscope)

[![Documentation](https://img.shields.io/badge/docs-mkdocs-526CFE?logo=materialformkdocs)](https://bbeckley-hub.github.io/acinetoscope)
[![RST Badge](https://img.shields.io/badge/documentation-RST-4CAF50.svg)](https://www.sphinx-doc.org/)
[![Last Commit](https://img.shields.io/github/last-commit/bbeckley-hub/acinetoscope)](https://github.com/bbeckley-hub/acinetoscope/commits)
[![Contributors](https://img.shields.io/github/contributors/bbeckley-hub/acinetoscope)](https://github.com/bbeckley-hub/acinetoscope/graphs/contributors)
[![Security: bandit](https://img.shields.io/badge/security-bandit-yellow.svg)](https://github.com/PyCQA/bandit)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![CI](https://img.shields.io/github/actions/workflow/status/bbeckley-hub/acinetoscope/ci.yml?branch=main&label=CI)](https://github.com/bbeckley-hub/acinetoscope/actions)
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)](https://github.com/bbeckley-hub/acinetoscope/tests)
[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/bbeckley-hub/acinetoscope)

[![Speed](https://img.shields.io/badge/Speed-40%25%20faster-FF6D00.svg)](https://github.com/bbeckley-hub/acinetoscope#performance-benchmarks)
[![A. baumannii](https://img.shields.io/badge/Species-A._baumannii-00BCD4.svg)](https://github.com/bbeckley-hub/acinetoscope)
[![CRAB](https://img.shields.io/badge/Detects-CRAB-F44336.svg)](https://github.com/bbeckley-hub/acinetoscope)
[![Risk Flagging](https://img.shields.io/badge/Risk-Flagging-FF9800.svg)](https://github.com/bbeckley-hub/acinetoscope)

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://docs.conda.io/en/latest/)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub Issues](https://img.shields.io/github/issues/bbeckley-hub/acinetoscope)](https://github.com/bbeckley-hub/acinetoscope/issues)
[![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/acinetoscope)](https://github.com/bbeckley-hub/acinetoscope/stargazers)
[![Sample Report](https://img.shields.io/badge/📊-View_Sample_Report-blue)](https://htmlpreview.github.io/?https://bbeckley-hub.github.io/acinetoscope/#summary)
![Profile Views](https://komarev.com/ghpvc/?username=bbeckley-hub&label=Profile%20Views&color=0e75b6&style=flat)

</div>

---

## 🎉 **What's New in v1.3.0** (July 2026)

### 📊 **Sample‑Centric Reporter – The Clinical Counterpart to Gene‑Centric Analysis**

Previous versions of AcinetoScope focused on **gene‑centric** reporting: each resistance or virulence gene is shown with all genomes that carry it. This is excellent for epidemiological surveillance, but it does not give a **per‑isolate** overview – you cannot see, at a glance, everything that a single genome carries.

The **sample‑centric reporter** fills this gap. It displays **each isolate** as an interactive box with:

- **Typing badges**: `Pasteur: ST 2`, `Oxford: ST 3`, `K: K1`, `O: OC1`, `Capsule: K1:OC1`.
- **All gene tables** for that isolate: AMR (AMRfinder + all ABRicate databases), Virulence (VFDB, Ecoli_VF), BACMET (biocide & heavy metal resistance), Plasmids (PlasmidFinder + APT), and Mutations.
- **Filtering** by sample name or by database (e.g., show only VFDB genes).
- **Horizontally scrollable tables** – no data truncation.
- **Export** each table to CSV with one click.
- **Print** any tab individually.

**Why it matters:** During an outbreak investigation, clinicians need to see the **complete resistance profile of a single isolate** – not just the gene frequencies. The sample‑centric reporter provides that **patient‑centred view** in seconds. It also makes results accessible to non‑bioinformaticians, facilitating communication between laboratory and clinical teams.

### 🔁 **Dynamic Grouping by Typing – From StaphScope to AcinetoScope**

We ported the **grouping feature** from our sister pipelines (StaphScope, Kleboscope). In any gene‑centric table (AMR, Virulence, BACMET, Plasmids, Mutations), you can now click a button to group the genome list by:

- **Pasteur ST**
- **Oxford ST**
- **K locus**
- **O locus**
- **Capsule type**
- **Combinations**: `ST‑K`, `ST‑O`, `K:O`, `ST‑K:O`

The genome list reorganises instantly, showing genomes under group headers. Click “Reset” to return to the flat list.

**Why it matters:** This reveals **clonal associations** at a glance. For example, grouping AMR genes by **capsule type** might show that `OXA-23` is only found in `K1:OC1` isolates – indicating a specific capsule lineage is driving carbapenemase spread. This is critical for understanding transmission dynamics, designing targeted interventions, and linking genotypes to phenotypes.

### 🧬 **Mutation Tab – Clinical Relevance at a Glance**

Point mutations (e.g., `gyrA_S81L`, `parC_E88K`, `pmrB`) are clinically important for fluoroquinolone, colistin, and tigecycline resistance. The new **mutation tab** in the gene‑centric reporter shows each unique mutation (gene + element name) with all genomes that carry it. You can filter by gene (gyrA, parC, rpoB, 23S rRNA, pmrAB), group by typing, and export to CSV.

**Why it matters:** Mutations are often the **first sign** of emerging resistance to last‑line antibiotics. The mutation tab allows you to monitor these signals in real time, without digging through raw AMR data.

### 🏠 **HPC/Cloud‑Ready Orchestrator – Temp‑Dir Execution**

The orchestrator has been refactored to be **Docker‑ and HPC‑friendly**. Every module now runs in a temporary directory (`/tmp/acinetoscope_*`). Results are copied back to your output folder; everything else is **automatically cleaned up**. This means AcinetoScope now works on:

- **Read‑only filesystems** (common in HPC clusters and container environments).
- **Shared clusters** where users cannot write to the installation directory.
- **Cloud environments** where ephemeral storage is preferred.

**Why it matters:** This removes the biggest barrier to running AcinetoScope on high‑performance computing infrastructure, making genomic surveillance accessible on any platform – from a laptop to a supercomputer.

### 📚 **Citation Overhaul – Full Credit Where It's Due**

We have added a **dedicated citation tab** in both reports, with colour‑coded cards for each tool/database, copy‑to‑clipboard buttons, and a complete suggested acknowledgement. This includes:

- Kaptive 3 (Stanton et al., 2025) – GitHub: `klebgenomics/Kaptive`
- APT (Lam et al., 2023)
- AMRFinderPlus, ABRicate, CARD, ResFinder, VFDB, PlasmidFinder, BacMet, MEGARes, ARG-ANNOT, EcOH, and Biopython.

**Why it matters:** Proper attribution encourages collaboration and ensures that the scientists who built the foundations of our field are recognised.

### 🧠 **Enhanced AI Guide**

We have expanded the AI Guide with detailed example questions for each analysis category (epidemiology, AMR, virulence, mutations) and ethical guidelines for using AI in genomic research – all with a touch of humour.

---

## 📋 **Table of Contents**

- [🎯 Overview](#-overview)
- [✨ Key Features](#-key-features)
- [🆕 What’s New in v1.3.0](#-whats-new-in-v130-july-2026)
- [⚡ Quick Start (CLI)](#-quick-start-cli)
- [🔧 Installation (CLI)](#-installation-cli)
- [🐳 Docker & Singularity Usage](#-docker--singularity-usage)
- [🔗 Integrated External Tools & Dependencies](#-integrated-external-tools--dependencies)
- [🚀 Usage Guide (CLI)](#-usage-guide-cli)
- [📁 Output Structure](#-output-structure)
- [🔍 Analytical Modules](#-analytical-modules)
- [🤖 AI Integration Guide](#-ai-integration-guide)
- [🔮 Future Development](#-future-development)
- [❓ Frequently Asked Questions](#-frequently-asked-questions)
- [🐛 Troubleshooting](#-troubleshooting)
- [📚 Citation](#-citation)
- [🙏 Acknowledgements](#-acknowledgements)
- [👥 Authors & Contact](#-authors--contact)
- [📄 License](#-license)
- [📚 Third-Party Tool Citations](#-third-party-tool-citations)

---

## 🎯 **Overview**

**AcinetoScope** is an automated, comprehensive bioinformatics pipeline designed specifically for the genomic analysis of *Acinetobacter baumannii*, a WHO **Critical Priority** pathogen responsible for devastating hospital‑acquired infections, especially in intensive care units. Carbapenem‑resistant *A. baumannii* (CRAB) is among the most difficult‑to‑treat multidrug‑resistant organisms.

### 🌍 **The Problem**
- **Fragmented Workflows**: Analyzing *A. baumannii* requires manual chaining of 6+ separate tools (MLST, Kaptive, AMRFinder, ABRicate, APT, etc.).
- **Interpretation Barrier**: Raw outputs from multiple tools need manual integration to form an epidemiological narrative.
- **Resource Barriers**: Generalist pipelines are slow and not optimised for *A. baumannii*.
- **Cloud/HPC Limitations**: Many pipelines fail on read‑only filesystems or leave behind clutter.

### 💡 **Our Solution**
AcinetoScope delivers:
- **✅ End-to-End Automation**: One command runs the entire analysis from raw FASTA to a consolidated report.
- **✅ *A. baumannii*-Optimized**: Pre‑configured with species‑specific databases and typing schemes (Pasteur & Oxford MLST, Kaptive K/O loci, APT plasmid typing).
- **✅ Actionable Intelligence**: sample-centric and gene‑centric tracking system.
- **✅ Speed & Efficiency**: Benchmarked **40-75% faster** than generalist pipelines.
- **✅ AI-Ready Outputs**: Interactive HTML reports designed for seamless exploration with modern AI browser extensions.
- **✅ HPC/Cloud‑Ready**: All modules run in temporary directories, leaving no trace behind.

**Perfect for**: Hospital outbreak investigation, public health surveillance, antimicrobial resistance (AMR) research, and clinical microbiology.

---

## ✨ **Key Features**

### 🔬 **Core Analytical Modules**
| Module | 🎯 Purpose | 📊 Key Outputs | ⚡ Speed |
| :--- | :--- | :--- | :--- |
| **Quality Control** | Assembly metric assessment & integrity checking | N50/N75, GC%, ambiguous bases, homopolymers | <1 min |
| **Dual MLST Typing** | Phylogenetic classification via Pasteur & Oxford schemes | Sequence Type (ST), International Clone (IC), novel alleles | <1 min |
| **Capsule (K/O) Typing** | Polysaccharide capsule & lipooligosaccharide typing via Kaptive 3 | K type, O type, locus coverage/identity | 1-2 min |
| **AMR Detection** | Comprehensive resistance gene detection with AMRFinderPlus | Carbapenemases, ESBLs, colistin/tigecycline resistance, 4‑tier risk flags | 2-3 min |
| **Multi‑DB Screening** | Screening across 11 curated databases via ABRicate | Virulence factors, plasmid replicons, metal/biocide resistance, stress regulators | 3-4 min |
| **Plasmid Typing (APT)** | *Acinetobacter*‑specific plasmid typing | Rep type families, per‑genome rep types | 1-2 min |
| **Gene‑centric Reporter** | Integrates all results into a gene‑centric interactive HTML report | Dynamic grouping by typing, mutation tab, co‑occurrence, combinations, high‑risk CRAB | Instant |
| **Sample‑centric Reporter** | Per‑isolate interactive boxes with typing badges and all gene tables | Patient‑centred view, filtering, export, print | Instant |

### 🚨 **Innovations for *A. baumannii* Surveillance**
- **Four‑Tier Risk Flagging**: Automatically categorises resistance genes (e.g., OXA‑23 → CRITICAL; *qacE* → ENVIRONMENTAL).
- **Environmental Co‑Selection Tracking**: Uniquely screens for heavy metal (*czc, mer, ars*) and biocide (*qac*) resistance genes.
- **Gene‑Centric Analysis Framework**: Tracks each resistance gene across all samples for clear visualisation of dissemination patterns.
- **Cross‑Genome Pattern Discovery**: Automatically identifies high‑risk combinations (e.g., carbapenemase + last‑resort resistance).
- **Dynamic Resource Allocation**: Uses Python's `psutil` to optimise parallel processing for any system (from laptops to HPC clusters).

---

## ⚡ **Quick Start (CLI)**

### **Install in 60 seconds**
```bash
# Conda (Recommended)
conda create -n acinetoscope -c conda-forge -c bioconda acinetoscope -y
conda activate acinetoscope
```

### **Run your first analysis**
```bash
# Single genome
acinetoscope -i genome.fasta -o results/

# Batch processing (24 genomes)
acinetoscope -i "*.fna" -o batch_results --threads 16
# Complete in ~10‑15 minutes! 🎉
```

---

## 🔧 **Installation (CLI)**

### **System Requirements**
| Resource | Minimum | Recommended | Production |
|----------|---------|-------------|------------|
| **CPU Cores** | 2 | 8+ | 16+ |
| **RAM** | 4 GB | 8 GB | 16 GB |
| **Storage** | 2 GB | 10 GB | 50 GB+ |
| **OS** | Linux, macOS, WSL2 | Linux | Linux Cluster |

### **Step-by-Step Installation**

#### **1. Install Miniconda (if needed)**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

#### **2. Install AcinetoScope**
```bash
# Add channels in correct order
conda config --add channels conda-forge
conda config --add channels bioconda

# Create and activate environment
conda create -n acinetoscope -c conda-forge -c bioconda acinetoscope -y
conda activate acinetoscope

# Verify installation
acinetoscope --help
```

#### **3. Update Databases (Recommended)**
```bash
# For ABRicate databases
abricate --setupdb

# For AMR database (first run or manual update)
acinetoscope --update-amr-db   # incremental
# or
acinetoscope --force-update-amr-db   # full overwrite
```

---

## 🐳 **Docker & Singularity Usage – avoid the padlock 🔓**

By default, Docker runs containers as `root`, so any files written to bind‑mounted directories will be owned by `root:root` – resulting in padlock icons and the need for `sudo chown`. **The fix is simple:** add `-u $(id -u):$(id -g)` to run the container with your host user’s UID/GID.

```bash
# Pull the latest image
docker pull bbeckleyhub/acinetoscope:latest

# Test installation
docker run --rm bbeckleyhub/acinetoscope:latest --help

# ✅ Recommended (no padlock, no sudo chown)
docker run --rm \
  -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  bbeckleyhub/acinetoscope:latest \
  -i "/data/*.fna" -o /data/output -t 4

# ❌ Old way (creates root‑owned files)
docker run --rm \
  -v $(pwd):/data \
  bbeckleyhub/acinetoscope:latest \
  -i "/data/*.fna" -o /data/output -t 4
# Then you need: sudo chown -R $USER:$USER ./output
```

**Why `-u $(id -u):$(id -g)`?**  
It tells Docker to run the container’s process with the same UID and primary GID as your host user. All files created in the mounted volume will be owned by **you** – no padlock, no permission errors, no cleanup.

> **Note for macOS / Windows (Docker Desktop):** UID/GID mapping works out‑of‑the‑box. The same command works fine.

---

### **Singularity / Apptainer (HPC clusters – no `sudo`, correct ownership)**  

Because AcinetoScope v1.3.0 writes all temporary files to `/tmp` (world‑writable), **you no longer need the `--writable-tmpfs` flag** (unless your cluster restricts `/tmp`). Singularity automatically maps your host user ID, so output files are **always** owned by you – no extra flags required.

```bash
# Pull the SIF image (once)
singularity pull acinetoscope.sif docker://bbeckleyhub/acinetoscope:latest

# Run – files are owned by your HPC user automatically
singularity run -B $(pwd):/data acinetoscope.sif -i "/data/*.fna" -o /data/output

# If your cluster restricts `/tmp`, add `--writable-tmpfs` for safety:
singularity run --writable-tmpfs -B $(pwd):/data acinetoscope.sif -i "/data/*.fna" -o /data/output
```

**Why Singularity users have no padlock:**  
Singularity/Apptainer **never** runs as root on HPC clusters; it always maps your host UID/GID into the container. The `-B` bind‑mount preserves ownership, so all output files are created with your user credentials.

---

### **Summary for HPC admins and Docker users**

| Platform | Command (recommended) | Ownership of output files |
|----------|----------------------|---------------------------|
| **Docker** | `docker run --rm -u $(id -u):$(id -g) -v "$PWD:/data" bbeckleyhub/acinetoscope ...` | Your user |
| **Singularity** | `singularity run -B "$PWD:/data" acinetoscope.sif ...` | Your user (automatically) |

No more `sudo chown`, no more padlock icons, no more angry HPC emails.

---

## 🔗 **Integrated External Tools & Dependencies**

| Tool/Database | Purpose | Source | License |
|---------------|---------|--------|---------|
| **MLST** | Multi-locus sequence typing | [tseemann/mlst](https://github.com/tseemann/mlst) | GPL v2 |
| **ABRicate** | Mass screening for resistance/virulence | [tseemann/abricate](https://github.com/tseemann/abricate) | GPL v2 |
| **AMRFinderPlus** | AMR gene detection | [ncbi/amr](https://github.com/ncbi/amr) | Public Domain |
| **Kaptive 3** | Capsule (K/O) typing | [klebgenomics/Kaptive](https://github.com/klebgenomics/Kaptive) | GPL v3 |
| **APT** | *Acinetobacter* plasmid typing | [MehradHamidian/AcinetobacterPlasmidTyping](https://github.com/MehradHamidian/AcinetobacterPlasmidTyping) | Free for research |
| **PubMLST** | MLST allele database | [pubmlst.org](https://pubmlst.org/organisms/acinetobacter-baumannii) | Open access for research |

---

## 🚀 **Usage Guide (CLI)**

### **Basic Commands**
```bash
# Single genome
acinetoscope -i genome.fna -o results/

# Batch processing with wildcards
acinetoscope -i "*.fna" -o results_2025 --threads 8

# Skip specific modules
acinetoscope -i sample.fna -o results --skip-qc --skip-plasmid

# Skip the new sample-centric reporter
acinetoscope -i "*.fna" -o results --skip-sample-centric

# AMR with custom thresholds and no mutations
acinetoscope -i "*.fna" -o results --amr-min-identity 0.95 --amr-min-coverage 0.9 --skip-amr-mutations
```

### **Input Formats**
- Accepted: `.fna`, `.fasta`, `.fa`, `.fn`
- Required: Assembled genomes (contigs or complete)
- Batch patterns: `*.fasta`, `sample_*.fna`, etc.

---

## 📁 **Output Structure**

AcinetoScope generates a well‑organized directory:

```
analysis/
├── fasta_qc_results/                 # Quality control reports per sample
├── PASTEUR_MLST/                     # MLST results (Pasteur scheme)
├── OXFORD_MLST/                      # MLST results (Oxford scheme)
├── kaptive_results/                  # Capsule (K) and lipooligosaccharide (O) typing
├── acineto_amrfinder_results/        # AMR gene detection with risk stratification
├── acineto_abricate_results/         # Multi‑database screening (11 DBs)
├── acineto_plasmid_results/          # Plasmid typing (APT) results
└── AcinetoScope_final_reports/       # 🎯 FINAL INTEGRATED REPORTS
    ├── GENIUS_ACINETOBACTER_GENE_CENTRIC_REPORTS/      # Gene‑centric HTML dashboard
    ├── GENIUS_ACINETOBACTER_SAMPLE_CENTRIC_REPORTS/    # Sample‑centric interactive boxes
    └── acinetoscope_run.log                             # Full log of the run
```
---

## 🔍 **Analytical Modules**

1. **Quality Control Module**: Validates input using *A. baumannii*‑specific thresholds (GC% 35-65%, ambiguous bases <5%).
2. **Dual MLST Typing Module**: Runs `mlst` with both **Pasteur** and **Oxford** schemes. Identifies International Clones (IC1-IC10).
3. **Capsule Typing Module**: Uses **Kaptive 3** with *A. baumannii*‑specific databases (`ab_k`, `ab_o`).
4. **AMR Detection Module**: Leverages **NCBI's AMRFinderPlus** with a four‑tier risk flagging system.
5. **Comprehensive Screening Module**: Executes **ABRicate** across 11 databases (CARD, ResFinder, VFDB, PlasmidFinder, BacMet2, etc.).
6. **Plasmid Typing Module**: Uses the **Acinetobacter Plasmid Typing (APT)** scheme to detect rep genes.
7. **Gene‑centric Reporter**: Synthesizes all results into a gene‑centric interactive HTML report with dynamic grouping, mutation tab, co‑occurrence, and combinations.
8. **Sample‑centric Reporter**: Displays each isolate as an interactive box with typing badges and all gene/mutation tables.

---

## 🤖 **AI Integration Guide**

AcinetoScope generates comprehensive HTML and JSON reports that are **perfect for AI analysis**. Here's how to use AI tools to get more from your data.

### 🚀 Quick Start
1. **Install any AI browser extension** (ChatGPT, Claude, Gemini)
2. **Open your report**: `genius_acinetobacter_ultimate_gene_centric_report.html` or the sample‑centric report.
3. **Select text** in any section (AMR Genes, MLST Analysis, etc.)
4. **Right‑click → Ask AI** with your question

### 💡 Example Questions

**For Epidemiology:**
- "What are the most common Pasteur STs in this dataset?"
- "Which capsule types (K:O) are dominant?"
- "Are there any ST‑capsule combinations with >2 isolates?"
- "Which clones carry the most resistance genes?"

**For AMR:**
- "How many samples carry OXA‑23? What are their STs?"
- "Are there any NDM‑positive samples? What is their capsule type?"
- "Which AMR genes co‑occur most frequently?"
- "What is the distribution of colistin resistance genes (mcr, pmr)?"

**For Virulence & Biofilm:**
- "Which samples carry biofilm genes (ompA, csu)?"
- "List all isolates with the T6SS (tss, hcp)."
- "Is there a correlation between capsule type and virulence gene carriage?"

**For Mutations:**
- "What are the most frequent point mutations in gyrA or parC?"
- "Are there any 23S rRNA mutations (tigecycline resistance)?"
- "Which samples carry qac genes (disinfectant resistance)?"

### 📊 Pro Tips
- **Provide context**: "I'm analysing *A. baumannii* genomics data..."
- **Be specific**: Instead of "tell me about this", ask "what does ST2 indicate in the context of CRAB?"
- **Use grouping**: "Use the grouping button 'ST' in the AMR tab and tell me which STs carry OXA‑23."
- **Ask for interpretations**: "What are the clinical implications of these findings?"
- **Request summaries**: "Summarise the resistance profile of sample XYZ."

> *"AI provides powerful insights but always verify critical findings with domain experts."*

---

## 🔮 **Future Development**

- **Raw read support** – Direct FASTQ analysis with integrated assembly.
- **Machine learning module** – Outbreak prediction, phenotype inference, risk scoring.
- **Real‑time database updates** – Live synchronisation of AMR and typing databases.
- **Expanded ESKAPE coverage** – Porting AcinetoScope's architecture to other ESKAPE pathogens.

---

## ❓ **Frequently Asked Questions**

**Q: Is AcinetoScope free to use?**  
A: Yes! Open‑source under MIT License. Free for academic, clinical, and commercial use.

**Q: What makes AcinetoScope different from other tools?**  
A: *A. baumannii*‑optimised, integrates 7 analysis types (including APT and Kaptive 3), runs 8‑10× faster, and now offers **both gene‑centric and sample‑centric reports** with **dynamic grouping** – features no other tool offers.

**Q: Can I use AcinetoScope for clinical diagnosis?**  
A: AcinetoScope is a research tool. While highly accurate, results should be validated with orthogonal methods for clinical decision‑making.

**Q: What is the sample‑centric reporter and why should I use it?**  
A: It displays each isolate as an interactive box with all its genes and typing badges – perfect for per‑patient investigation, outbreak case reviews, and presenting results to clinicians.

**Q: How do I use the grouping feature?**  
A: In the gene‑centric report, open any gene table (AMR, Virulence, etc.). Above the table you’ll see buttons for “ST”, “K Locus”, “Capsule”, and combinations. Click one – the genome list reorganises instantly by that typing scheme.

**Q: Why does v1.3.0 no longer require `--writable-tmpfs` in Singularity?**  
A: All modules now write temporary files to `/tmp` (not to the installation directory). Containers mount `/tmp` as writable by default, so no special flags are needed.

---

## 🐛 **Troubleshooting**

### **Common Issues & Solutions**

```bash
# Issue: AMR database missing or outdated
# Solution:
acinetoscope --force-update-amr-db

# Issue: ABRicate database not found
# Solution:
abricate --setupdb

# Issue: Permission errors in Docker
# Solution: Ensure bind mounts are correct and use --user if needed
docker run --rm -u $(id -u):$(id -g) -v ... bbeckleyhub/acinetoscope ...

# Issue: Cross‑run contamination in /tmp
# Solution: Use --keep-temp only for debugging; otherwise temp dirs are auto‑deleted.
```

### **Getting Help**
1. **Check existing issues**: [GitHub Issues](https://github.com/bbeckley-hub/acinetoscope/issues)
2. **Search closed issues**: Many problems already solved
3. **Create new issue**: Include:
   - Full error message
   - Conda environment list (`conda list`)
   - Example command that failed
   - The `acinetoscope_run.log` file
4. **Email support**: brownbeckley94@gmail.com (response within 48 hours)

---

## 📚 **Citation**

If you use AcinetoScope in your research, please cite:

```bibtex
@software{acinetoscope2026,
  title = {A species‑specific computational pipeline for rapid, comprehensive Acinetobacter baumannii* outbreak investigation and resistance gene tracking},
  author = {Beckley, B. and Amarh, V. and Lopes, B. S. and A. and Olalekan, A. and Afeke, I.},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/bbeckley-hub/acinetoscope}
}
```

---

## 🙏 **Acknowledgements**

AcinetoScope stands on the shoulders of giants. We are deeply grateful to:

- **Torsten Seemann** for MLST, ABRicate, and countless foundational tools.
- **NCBI team** for AMRFinderPlus.
- **Kath Holt Lab** for Kaptive 3.
- **Hamidian Lab** for the APT scheme.
- **PubMedST, CARD, VFDB, ResFinder, PlasmidFinder, BacMet, MEGARes, ARG-ANNOT** for essential databases.
- **Python community** for Biopython, pandas, plotly, etc.
- **Early adopters and beta testers** for invaluable feedback.

---

## 👥 **Authors & Contact**

**Brown Beckley** (Primary Developer)  
- University of Ghana Medical School  
- 📧 brownbeckley94@gmail.com  
- 🐙 GitHub: [bbeckley-hub](https://github.com/bbeckley-hub)  
- LinkedIn: [@brownbeckley](https://www.linkedin.com/in/brown-beckley-190315319/)  
- 📞 +233 508820617

**Collaborators**: Vincent Amarh, Bruno Silvester Lopes, Adesola Olalekan, Innocent Afeke.

We welcome collaborations on MRSA epidemiology, clinical validation, global surveillance, and expansion to other ESKAPE pathogens.

---

## 📄 **License**

### Core AcinetoScope Code
The AcinetoScope pipeline code (the workflow engine, report generation, HTML templates, and Python modules written by the authors) is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.

### Third‑Party Tools
AcinetoScope executes several external bioinformatics tools, which are installed as Conda dependencies. Each tool is the property of its respective developers and is used under its own license. By using AcinetoScope, you agree to comply with the licenses of these third‑party tools.

---

### 📚 **Third‑Party Tool Citations**

AcinetoScope integrates several powerful open‑source tools and databases. If you use AcinetoScope in your research, please also cite the following essential tools:

#### **Kaptive 3**
```bibtex
@article{stanton_kaptive3_2025,
  author = {Stanton, T. D. and Hetland, M. A. K. and Löhr, I. H. and Holt, K. E. and Wyres, K. L.},
  title = {Fast and Accurate in silico Antigen Typing with Kaptive 3},
  journal = {Microbial Genomics},
  volume = {11},
  number = {6},
  pages = {001428},
  year = {2025},
  doi = {10.1099/mgen.0.001428},
  url = {https://github.com/klebgenomics/Kaptive}
}
```

#### **APT**
```bibtex
@article{lam_apt_2023,
  author = {Lam, M. M. C. and Koong, J. and Holt, K. E. and Hall, R. M. and Hamidian, M.},
  title = {Detection and Typing of Plasmids in Acinetobacter baumannii Using rep Genes Encoding Replication Initiation Proteins},
  journal = {Microbiology Spectrum},
  volume = {11},
  number = {1},
  pages = {e0247822},
  year = {2023},
  doi = {10.1128/spectrum.02478-22}
}
```

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

#### **PubMLST (Jolley et al.)**
```bibtex
@article{jolley_pubmlst_2018,
  author = {Jolley, K. A. and Bray, J. E. and Maiden, M. C. J.},
  title = {Open‑access bacterial population genomics: {BIGSdb} software, the {PubMLST.org} website and their applications},
  journal = {Wellcome Open Research},
  volume = {3},
  pages = {124},
  year = {2018},
  doi = {10.12688/wellcomeopenres.14826.1}
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

#### **BacMet**
```bibtex
@article{pal_bacmet_2014,
  author = {Pal, C. et al.},
  title = {BacMet: antibacterial biocide and metal resistance genes database},
  journal = {Nucleic Acids Research},
  volume = {42},
  number = {D1},
  pages = {D737-D743},
  year = {2014},
  doi = {10.1093/nar/gkt1252}
}
```

---

<div align="center">

## **🚀 Ready to revolutionise your *A. baumannii* surveillance?**

| **Choose Your Platform** | |
|--------------------------|-|
| 🖥️ **Command Line** | For high‑throughput, local analysis |
| 🐳 **Docker / Singularity** | For HPC and containerised environments |

[![Get Started CLI](https://img.shields.io/badge/GET_STARTED_CLI-Now-green?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/acinetoscope#-quick-start-cli)
[![Report Issue](https://img.shields.io/badge/REPORT_ISSUE-Here-red?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/acinetoscope/issues)

**From days to minutes. From fragmented to integrated. From data to insights.**

*AcinetoScope: Precision surveillance for the carbapenem‑resistant era.*

⭐ **If you find this tool useful, please star the repository!** ⭐

*Join the Fight Against Antimicrobial Resistance*

Antimicrobial resistance (AMR) represents one of the most significant global health threats of our time. We invite researchers, clinicians, and public health professionals to collaborate with us in expanding and validating our database, sharing regional epidemiological data, and advancing AMR surveillance.

**Together, we can enhance global AMR monitoring and develop more effective treatment strategies.**

</div>
