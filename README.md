
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
<div align="center">

# 🦠 AcinetoScope
![Conda Version](https://anaconda.org/bbeckley-hub/acinetoscope/badges/version.svg)
![Conda Downloads](https://anaconda.org/bbeckley-hub/acinetoscope/badges/downloads.svg)
![Platform](https://anaconda.org/bbeckley-hub/acinetoscope/badges/platforms.svg)
![License](https://img.shields.io/github/license/bbeckley-hub/acinetoscope)
![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)
![Bioconda](https://img.shields.io/badge/bioconda-✓-44A833.svg?logo=conda)
![Last updated](https://anaconda.org/bbeckley-hub/acinetoscope/badges/latest_release_date.svg)

![Docker Pulls](https://img.shields.io/docker/pulls/bbeckleyhub/acinetoscope)
![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/acinetoscope)
![GitHub Issues](https://img.shields.io/github/issues/bbeckley-hub/acinetoscope)
![Research Use](https://img.shields.io/badge/research%20use-✓-purple)

</div>

## 🚀 Revolutionizing A. baumannii Genomic Analysis

**AcinetoScope** is a comprehensive, parallelized bioinformatics pipeline that solves the fragmentation problem in *Acinetobacter baumannii* genomic analysis. By integrating multiple analysis modules into a unified, automated workflow, it provides clinical microbiologists and researchers with a complete solution from raw sequencing data to meaningful insights.

```
📊 BEFORE: 6 separate tools, manual integration, hours of work
╔══════════════════════════════════════════════════════════════════════════════╗
║  FASTQC → MLST → Kaptive → AMRFinder → ABRicate → Summary Reports            ║
║  │         │        │          │           │              │                  ║
║  │         │        │          │           │              │                  ║
║  │         │        │          │           │              │                  ║
║  └─────┬───┴───┬────┴────┬─────┴─────┬─────┴──────┬───────┘                  ║
║        │       │         │           │            │                          ║
║  6 HOURS     MANUAL    ERROR     INCONSISTENT  FRAGMENTED                    ║
║              MERGE     PRONE       OUTPUTS       RESULTS                     ║
╚══════════════════════════════════════════════════════════════════════════════╝

📊 AFTER: 1 command, automated integration, parallel execution
╔══════════════════════════════════════════════════════════════════════════════╗
║                    acinetoscope -i "*.fna" -o results --threads 8            ║
║                                                                              ║
║                    🚀 PARALLEL EXECUTION (8 workers)                         ║
║                    ┌─────────┬─────────┬─────────┬─────────┐                 ║
║                    │  QC     │ MLST    │ Kaptive │ AMR     │                 ║
║                    │         │         │         │         │                 ║
║                    │  MODULE │ MODULE  │ MODULE  │ MODULE  │                 ║
║                    └─────────┴─────────┴─────────┴─────────┘                 ║
║                    ┌─────────┬─────────┬───────────────────┐                 ║
║                    │ABRicate │Summary  │    45 MINUTES     │                 ║
║                    │ MODULE  │MODULE   │  (not 6 hours!)   │                 ║
║                    └─────────┴─────────┴───────────────────┘                 ║
║                                                                              ║
║                    ✅ COMPLETE  ✅ INTEGRATED  ✅ CONSISTENT                 ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

## 🌟 Why AcinetoScope?

### The Problem: Bioinformatics Fragmentation
Current *A. baumannii* analysis requires running 6+ separate tools, manually integrating results, dealing with inconsistent outputs, and spending hours on data processing. This fragmentation leads to:
- **Time-consuming manual integration**
- **Inconsistent result formats**
- **Missed critical genetic markers**
- **Difficulty in outbreak tracking**
- **Reproducibility challenges**

### The Solution: Unified Parallel Pipeline
AcinetoScope eliminates fragmentation by providing:
- **🚀 Parallel Execution**: All modules run simultaneously
- **🎯 Critical Gene Tracking**: 300+ scientifically curated markers
- **📊 Interactive Reports**: Gene-centric HTML summaries
- **⚡ One-Command Analysis**: From FASTA to final reports
- **🔬 Scientific Accuracy**: Research-grade analysis pipelines

## 🧬 Comprehensive Critical Gene Tracking

AcinetoScope tracks **300+ scientifically curated critical genes** across multiple categories:

### 🚨 Priority 1: Life-Threatening Resistance
```
╔══════════════════════════════════════════════════════════════════════════════╗
║                     CRITICAL CARBAPENEMASE SURVEILLANCE                      ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  OXA-type (23, 24, 40, 51, 58, 66, 69, 71, 143, 235, 236, 237, 267, 317)     ║
║  Metallo-β-lactamases (NDM-1-5, VIM-1-4, IMP-1-5, KPC-2-4)                   ║
║  Other carbapenemases (GES, SIM, SPM, AIM)                                   ║ 
╚══════════════════════════════════════════════════════════════════════════════╝
```

### 💊 ESBL & AmpC Detection
- **ESBLs**: CTX-M (1,2,3,9,14,15,27,55), SHV (1,2,5,11,12), TEM (1,2,10,52), PER (1-3), VEB (1-2)
- **AmpC**: ADC (1,2,5,7,10,11,30,69,75,88,176)

### 🛡️ Last-Resort Antibiotic Resistance
```
╔══════════════════════════════════════════════════════════════════════════════╗
║                   LAST LINE OF DEFENSE - TRACKING GENES                      ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  🚨 COLISTIN: mcr-1 to mcr-10, pmrABC, lpxACD, eptA, arnA-F, phoPQ           ║
║  🚨 TIGECYCLINE: tet(X1-6), tet(39), adeRSABC, adeJKNT                       ║
║  🏥 BIOFILM: ompA, csuABCDE, bfmRS, pgaABCD, bap, pilABCDEF                  ║
║  🧬 EFFLUX PUMPS: adeA-N, abeMS, amvA, craA, mexJKT, mdeA, mdfA/cmr          ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

### 🏥 Hospital Environment Survival Genes
```
╔══════════════════════════════════════════════════════════════════════════════╗
║            ENVIRONMENTAL CO-SELECTION & HOSPITAL PERSISTENCE                 ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  🧪 Biocide Resistance: qacA-J, qacEA1, cepA, formABC, oqxAB                 ║
║  ⚙️ Metal Resistance: czcABCD, merA-J, arsA-J, copA-J, zntA-J                ║
║  🧬 Stress Response: soxRS, marABC, robA, rpoS, phoBR                        ║
║  🔄 Mobile Elements: traA-X, mobA-H, intI1-3, tnpA-F, istAB                  ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

### 🔬 Multi-Database Integration
```
╔══════════════════════════════════════════════════════════════════════════════╗
║                    COMPREHENSIVE DATABASE COVERAGE                           ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  🎯 AMR: AMRFinder, CARD, ResFinder, NCBI, MEGARes                           ║
║  🦠 Virulence: VFDB, Victors, BacMet2, EcoH_VF                               ║
║  🧬 Typing: MLST (Oxford & Pasteur), Kaptive K/O Locus                       ║
║  📊 Quality: Comprehensive FASTA QC with multi-metric scoring                ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

## 🚀 Key Features

### 📊 Interactive HTML Reports
```
╔══════════════════════════════════════════════════════════════════════════════╗
║                    GENE-CENTRIC INTERACTIVE DASHBOARD                        ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  📈 Interactive HTML reports   🔬 Cross-genome pattern discovery             ║
║  🎯 Risk Stratification        📊 Multi-Sample Comparisons                   ║
║  🧬 Phylogenetic Context       🏥 Clinical Relevance Scoring                 ║
║  📋 Export to CSV/JSON         🖨️ Statistical analysis                       ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

### 🎯 Scientific Accuracy & Curation
- **300+ scientifically validated gene markers**
- **BACMET2-compliant environmental resistance tracking**
- **Clinically relevant gene categorization**
- **Research-grade analysis pipelines**

## 📦 Installation

### Quick Install with Conda
```bash
conda create -n acinetoscope -c bioconda -c conda-forge acinetoscope
conda activate acinetoscope
```

### From Source
```bash
git clone https://github.com/bbeckley-hub/acinetoscope.git
cd acinetoscope
pip install -e .
```

### Docker
```bash
docker pull bbeckleyhub/acinetoscope:latest
docker run -v $(pwd)/data:/data bbeckley/acinetoscope \
  -i /data/*.fna -o /data/results --threads 8
```

## 🚦 Quick Start

### Single Genome Analysis
```bash
acinetoscope -i sample.fasta -o results/ --threads 4
```

### Batch Processing
```bash
acinetoscope -i "*.fna" -o batch_results/ --threads 8
```

### Specific Module Selection
```bash
# Skip QC and summary for faster analysis
acinetoscope -i "*.fasta" -o results/ --threads 6 --skip-qc --skip-summary

# Only MLST Pasteur scheme
acinetoscope -i genomes/ -o analysis/ --threads 4 --mlst-scheme pasteur
```

## 📁 Output Structure

```
results/
├── 📊 fasta_qc_results/                    # Quality control reports
├── 🧬 PASTEUR_MLST/                        # MLST Pasteur scheme results
├── 🧬 OXFORD_MLST/                         # MLST Oxford scheme results  
├── 🔬 kaptive_results/                     # K/O locus typing
├── 💊 acineto_amrfinder_results/          # AMR gene detection
├── 🦠 acineto_abricate_results/           # Multi-database screening
└── 📈 GENIUS_ACINETOBACTER_ULTIMATE_REPORTS/
    ├── 📄 gene_centric_summary.html       # Interactive dashboard
    ├── 📋 complete_analysis.json          # Machine-readable output
    
```

## 🔬 Analysis Modules

### 1. 🧬 MLST Typing Module
- **Dual scheme support**: Oxford & Pasteur
- **Batch processing**: Multiple genomes simultaneously with glob patterns "*.fasta"
- **HTML summaries**: Interactive typing reports
- **Sequence type tracking**: Outbreak cluster detection

### 2. 🔬 Kaptive K/O Locus Module
- **Capsule typing**: Complete K locus analysis
- **Lipooligosaccharide**: O locus determination
- **Confidence scoring**: Quality metrics for typing

### 3. 💊 AMR Detection Module
- **AMRFinderPlus**: Comprehensive resistance gene detection
- **Critical gene highlighting**: Priority resistance markers
- **Percentage-based reporting**: Sample frequency analysis
- **Risk stratification**: Clinical relevance scoring

### 4. 🦠 Multi-Database Screening Module
- **8 databases**: CARD, ResFinder, NCBI, VFDB, ARG-ANNOT, MEGARES, VICTORS, BacMet2
- **Environmental markers**: Hospital persistence genes
- **Plasmid screening**: Mobile genetic element detection
- **Virulence factors**: Pathogenicity assessment

### 5. 📊 Quality Control Module
- **FASTA validation**: Sequence integrity checking
- **Contamination screening**: Multi-metric quality scoring
- **Completeness assessment**: Assembly quality metrics
- **Visual reports**: Interactive QC dashboards

### 6. 📈 Ultimate Reporter Module
- **Gene-centric integration**: All results in one view
- **Interactive filtering**: Dynamic data exploration
- **Export functionality**: CSV, JSON, HTML outputs

## 🏥 Clinical & Research Applications

### Hospital Infection Control
```
╔══════════════════════════════════════════════════════════════════════════════╗
║                    OUTBREAK DETECTION WORKFLOW                               ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  1. 📥 Upload patient isolates (n=50)                                        ║
║  2. 🚀 Run: acinetoscope -i "*.fna" -o outbreak/ --threads 8                 ║
║  3. 🔬 Detect: ST2/OXA-23 cluster in 12/50 samples                           ║
║  4. 🎯 Identify: mcr-1 in 3/50 → immediate infection control measures        ║
║  5. 📊 Report: Interactive dashboard for hospital epidemiology team          ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

### Research Studies
- **Longitudinal surveillance**: Track resistance evolution
- **Geographic distribution**: Map resistance gene spread  
- **Treatment efficacy**: Correlate genotypes with outcomes
- **Novel gene discovery**: Identify emerging resistance markers

### Public Health
- **Regional surveillance**: Monitor resistance trends
- **Early warning systems**: Detect emerging threats
- **Policy guidance**: Data-driven antibiotic stewardship
- **Global tracking**: International resistance monitoring

## 🧪 Scientific Validation

### Peer-Reviewed Methodology
- **MLST**: Based on Institut Pasteur & Oxford University schemes
- **AMR detection**: Validated against AMRFinderPlus gold standard
- **Kaptive**: Published typing algorithm for *A. baumannii*
- **Gene categorization**: Clinically validated marker lists

## 📚 Citation & Acknowledgements

If you use AcinetoScope in your research, please cite:

```bibtex
@software{acinetoscope2026,
  title = {AcinetoScope: A Unified Parallelized Pipeline for Comprehensive Genomic Surveillance of Acinetobacter baumannii and Tracking of Critical Resistance Determinants},
  author = {Beckley, Brown},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/bbeckley-hub/acinetoscope},
  note = {Complete A. baumannii genomic analysis pipeline with parallel execution}
}
```

### 🙏 Acknowledgments
- **University of Ghana Medical School** - Department of Medical Biochemistry
- **Global AMR Surveillance Network** - Methodology validation
- **Open Source Bioinformatics Community** - Tool integration
- **Clinical Microbiology Partners** - Real-world testing and feedback

## 🐛 Bug Reports & Feature Requests

Found a bug or have a feature request? Please open an issue on our [GitHub repository](https://github.com/bbeckley-hub/acinetoscope/issues).

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

## 📄 License

AcinetoScope is released under the MIT License. See [LICENSE](LICENSE) for details.

## 📞 Contact & Support

- **Author**: Brown Beckley
- **Email**: brownbeckley94@gmail.com  
- **Affiliation**: University of Ghana Medical School - Department of Medical Biochemistry
- **GitHub**: [bbeckley-hub](https://github.com/bbeckley-hub)
- **Documentation**: [GitHub Wiki](https://github.com/bbeckley-hub/acinetoscope/wiki)

---

</div>
"Every microbe has its story. We provide the translation."
