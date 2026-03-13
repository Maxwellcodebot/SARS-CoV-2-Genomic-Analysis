

# 🦠 SARS-CoV-2 Genomic Analysis Pipeline

A complete bioinformatics pipeline for analyzing SARS-CoV-2 genome sequences, comparing variants (Alpha, Delta, Omicron BA.1, Omicron XBB), and visualizing evolutionary patterns.

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![BioPython](https://img.shields.io/badge/BioPython-1.81-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

---

## 📋 Table of Contents
1. [Project Overview](#-project-overview)
2. [Key Findings](#-key-findings)
3. [Requirements](#-requirements)
4. [How to Run](#-how-to-run)
5. [Repository Structure](#-repository-structure)
6. [Output Files Explained](#-output-files-explained)
7. [Biological Interpretation](#-biological-interpretation)
8. [For Beginners: Understanding the Code](#-for-beginners-understanding-the-code)
9. [Future Improvements](#-future-improvements)
10. [References](#-references)
11. [Author](#-author)
12. [License](#-license)
13. [Acknowledgments](#-acknowledgments)
14. [Contact](#-contact)
15. [Quick Start](#-quick-start-one-liner)

---

## 📊 Project Overview

This project analyzes **5 SARS-CoV-2 genomes** (Wuhan-Hu-1 reference, Alpha, Delta, Omicron BA.1, and Omicron XBB) to identify:

- 🧬 **GC Content variation** across variants  
- 🔬 **Point mutations** relative to the reference  
- 🗺️ **Mutation hotspots** in the genome  
- 📈 **Evolutionary patterns** among variants  

The pipeline generates **4 publication-quality figures** and comprehensive CSV reports automatically.

---

## 🔬 Key Findings

| Variant | GC Content | Mutations | Notable Finding |
|---------|------------|-----------|----------------|
| **Reference (Wuhan-Hu-1)** | 37.973% | - | Original strain |
| **Alpha Variant** | 37.982% | 84 | Lowest mutation count |
| **Delta Variant** | 37.964% | 102 | Contains 40 ambiguous bases (N's) |
| **Omicron BA.1** | **37.989%** | **167** | **Highest GC% and mutations** |
| **Omicron XBB** | 37.953% | 126 | Lowest GC% among Omicron lineage |

### 🎯 Major Discovery
A significant **mutation hotspot** was identified at **29500-30000 bp**, corresponding to the **spike protein gene and 3'UTR**. Omicron variants show the highest mutation density in this region, correlating with increased transmissibility and immune evasion.

---

## 🛠️ Requirements

```bash
biopython==1.81
pandas==2.0.3
numpy==1.24.3
matplotlib==3.7.2
seaborn==0.12.2
tabulate==0.9.0

Install with:

pip install -r requirements.txt


---

🚀 How to Run

1. Clone the repository

git clone https://github.com/Maxwellcodebot/SARS-CoV-2-Genomic-Analysis.git
cd SARS-CoV-2-Genomic-Analysis

2. Install dependencies

pip install -r requirements.txt

3. Run the analysis

python sars_cov2_analysis.py

The script will automatically:

✅ Read all FASTA files

✅ Calculate GC content for each variant

✅ Identify mutations compared to reference

✅ Generate 4 figures

✅ Save all results to CSV files



---

📁 Repository Structure

SARS-CoV-2-Genomic-Analysis/
│
├── sars_cov2_analysis.py               # Main analysis script (beginner-friendly)
├── README.md                           # This file
├── requirements.txt                    # Python package dependencies
│
├── COVID_Sequence_Reference.fasta      # Wuhan-Hu-1 reference genome
├── COVID_Sequence_Alpha_Variant.fasta  # Alpha variant
├── COVID_Sequence_Delta_Variant.fasta  # Delta variant
├── Omicron_BA.1_COVID_Variant.fasta    # Omicron BA.1
├── Omicron_XBB_COVID_Variant.fasta     # Omicron XBB
├── sequence.gb                          # GenBank annotation file
│
├── Complete_GC_Report.csv               # GC content summary table
├── gc_content_bar_chart.png             # Figure 1: GC comparison
├── mutation_count_chart.png             # Figure 2: Mutation counts
├── mutation_distribution.png            # Figure 3: Mutation positions
├── mutation_heatmap.png                 # Figure 4: Mutation density heatmap
├── figure_captions.txt                  # All figure captions
│
└── mutation_results/                     # Folder with mutation data
    ├── Alpha_Variant_mutations.csv
    ├── Delta_Variant_mutations.csv
    ├── Omicron_BA.1_mutations.csv
    ├── Omicron_XBB_mutations.csv
    ├── mutation_summary.csv
    └── *_annotated_mutations.csv        # Mutations with gene info (if available)


---

📊 Output Files Explained

CSV Files

File	Description

Complete_GC_Report.csv	GC content for all variants with statistics
mutation_results/mutation_summary.csv	Mutation counts per variant
mutation_results/*_mutations.csv	Individual mutation lists (position, reference, variant)
mutation_results/*_annotated_mutations.csv	Mutations with gene info (if available)
gene_annotations.csv	Gene positions from GenBank file


Figures

Figure	File	Description

Figure 1	gc_content_bar_chart.png	Bar chart comparing GC% across variants
Figure 2	mutation_count_chart.png	Bar chart showing mutation counts
Figure 3	mutation_distribution.png	Scatter plot of mutation positions
Figure 4	mutation_heatmap.png	Heatmap of mutation density (500bp bins)



---

🧪 Biological Interpretation

1. GC Content (37.95% - 37.99%)

Omicron BA.1 has the highest GC% (37.989%) → potentially more stable genome

Omicron XBB has the lowest GC% (37.953%) → divergent evolution within Omicron lineage

Very narrow range shows evolutionary constraint on base composition



2. Mutation Burden

Alpha: 84 mutations (lowest)

Delta: 102 mutations

Omicron XBB: 126 mutations

Omicron BA.1: 167 mutations (highest)

Pattern: Mutation count increased over time → accelerated evolution in later variants



3. Mutation Hotspot (29500-30000 bp)

Located in spike protein gene / 3'UTR

Omicron variants show highest density here → important for transmission and immune evasion



4. Quality Control

Delta variant contains 40 'N' bases → possible sequencing quality issue





---

👨‍💻 For Beginners: Understanding the Code

# GC content = how many G and C letters in the DNA
# Higher GC = more stable genome

def GC_content_calculator(sequence):
    # Count G's and C's, divide by total, multiply by 100
    return (sequence.count("G") + sequence.count("C")) / len(sequence) * 100


---

🔄 Future Improvements

[ ] Add more variants (Beta, Gamma, other Omicron sublineages)

[ ] Analyze multiple sequences per variant for statistics

[ ] Add transition/transversion ratio analysis

[ ] Create interactive dashboard with Plotly

[ ] Add protein-level mutation analysis

[ ] Compare with published literature automatically



---

📚 References

NCBI GenBank: https://www.ncbi.nlm.nih.gov/

SARS-CoV-2 Reference Genome: NC_045512.2

WHO SARS-CoV-2 Variants: https://www.who.int/activities/tracking-SARS-CoV-2-variants



---

👤 Author

Your Name

GitHub: @Maxwellcodebot

Project Link: GitHub Repo



---

📄 License

This project is licensed under the MIT License - see the LICENSE file for details.


---

⭐ Acknowledgments

Thanks to NCBI for providing genomic data

BioPython documentation and community

All scientists working on SARS-CoV-2 genomics



---

📬 Contact

Questions or suggestions? Reach out via:

GitHub Issues

Email: maxwellwisaac@gmail.com



---

⚡ Quick Start (One-Liner)

pip install biopython pandas numpy matplotlib seaborn tabulate && python sars_cov2_analysis.py


---

If you find this project useful, please give it a star! ⭐

---

