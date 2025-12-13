# insilicopcr
InSilicoPCR – Fast Parallel In-Silico PCR Tool: Python tool for scanning multiple genomes with multiple primer pairs and producing multi-FASTA amplicons.

# 🔬 InSilicoPCR – Parallel Multi-Genome In-Silico PCR Tool  
### Developed by **Dr. Vishal**  
**Assistant Professor, School of Forensic Science  
National Forensic Sciences University (NFSU), Gandhinagar, India**

---

## 📘 Overview

**InSilicoPCR** is a high-performance Python tool designed for fast, parallel *in-silico PCR* across multiple genomes.  
It accepts:

- A directory containing genome FASTA files  
- A CSV file containing primer pairs  
- Optional mismatch, size, and CPU parameters  

It outputs:

- **Amplicon multi-FASTA files per primer pair**  
- **A comprehensive CSV summary** with amplicon metadata  
- IUPAC-aware matching (degenerate primers supported)  
- Full parallelization across all available CPU cores  

This tool is especially useful in forensic genomics, microbial identification, DNA barcoding, wildlife forensics, assay design, and molecular biology research.

---

## ✨ Features

- 🔥 **High-speed parallel processing**
- 🧬 **IUPAC degenerate primer support** (R, Y, N, etc.)
- 📁 **Multi-FASTA output per primer pair**
- 📊 **Amplicon summary table**
- 🎯 **Configurable mismatches and amplicon size**
- ⚙️ **Zero external dependencies** (pure Python)
- 🧪 Ideal for forensic biology, teaching, and research workflows

---

## 📂 Directory Structure
InSilicoPCR/
│
├── insilico_pcr_parallel_multifasta.py # Main tool
├── primers.csv # Example primer file
├── example_genomes/ # Example FASTA files
├── LICENSE # MIT License
└── README.md
---

## 🧬 Primer CSV Format

Primer file must contain at least:
```bash
pair_id,forward,reverse
16S,AGAGTTTGATCCTGGCTCAG,GGTTACCTTGTTACGACTT
gyrB,CTTCGACATCGACGACGA,ACGACGACGACTTCCAG
```
Column names can include variations like:
```bash
`forward, fwd, left`, and `reverse, rev, right`.  
`pair_id` is optional — will be auto-generated if missing.
```
---

## 📁 FASTA Input Directory

Place any number of FASTA files in a folder:
```bash
genomes/
├── strain1.fasta
├── strain2.fa
└── environmental_sample.fna
```
Supported extensions: `.fa`, `.fasta`, `.fna`, `.ffn`

---

## ⚡ Usage

### **Basic command**
```bash
python insilico_pcr_parallel_multifasta.py \
    --fasta_dir genomes \
    --primers primers.csv \
    --out_dir output

python insilico_pcr_parallel_multifasta.py \
    --fasta_dir genomes \
    --primers primers.csv \
    --out_dir output \
    --max_mismatch 2 \
    --min_len 50 \
    --max_len 2000 \
    --workers 0
```
--workers 0 automatically uses all CPU cores.


The tool produces two main types of output:

1. **Multi-FASTA files** (one per primer pair)  
2. **A combined summary CSV** containing metadata for every detected amplicon  

All output files are saved inside the directory specified using `--out_dir`.

---

### **1. Multi-FASTA Files**

For each primer pair, the tool generates a separate multi-FASTA file containing all amplicons detected from all genome FASTA files.

For each primer pair, the tool generates a separate multi-FASTA file containing all amplicons detected across all input genome sequences.

Example output files:
16S_amplicons.fasta
gyrB_amplicons.fasta
toxA_amplicons.fasta

Each FASTA file contains one or more amplicons.  
The header format is structured to include primer ID, sample name, hit number, amplicon length, and primer sequences:

H16S_Strain1_hit1|len=452|fwd=AGAGTTTGATCCTGGCTCAG|rev=GGTTACCTTGTTACGACTT
ATGCGTACGTT...


**Header Components Explained:**

| Component | Meaning |
|----------|---------|
| `16S` | Primer pair ID |
| `Strain1` | Genome/sample name extracted from FASTA filename |
| `hit1` | First amplicon detected for this primer pair in this genome |
| `len=452` | Amplicon length in base pairs |
| `fwd=...` | Forward primer sequence |
| `rev=...` | Reverse primer sequence |

Amplicon sequences are automatically wrapped at 80 characters per line for FASTA compatibility.
ATGCGTACGTTCAGAGTTGATCCTG...
CGTACGTAGCAGTAGCTAGCTGATGA...

If no amplicons are found for a primer pair, its corresponding file will be created but remain empty.

### **2. Summary CSV File**

A single file named:
amplicons_summary.csv
contains metadata for **every amplicon discovered** across all primer pairs and genomes.

**Example rows:**
**Column Description Table:**

| Column | Description |
|--------|-------------|
| `pair_id` | Primer pair name or ID |
| `fasta_file` | Full path of genome FASTA file |
| `sample_name` | Sample name extracted from filename |
| `seq_id` | Sequence identifier inside FASTA |
| `fwd_start`, `fwd_end` | Coordinates of forward primer match |
| `rev_start`, `rev_end` | Coordinates of reverse primer match |
| `fwd_mm`, `rev_mm` | Number of mismatches for each primer |
| `amp_len` | Length of amplicon (bp) |
| `forward_primer` | Exact forward primer sequence |
| `reverse_primer` | Exact reverse primer sequence |

This file is ideal for downstream analysis, reporting, forensics documentation, and pipeline integration.

---

### **3. Output Directory Structure**

After running the tool, the output directory may look like this:
```bash
output/
├── 16S_amplicons.fasta
├── gyrB_amplicons.fasta
├── toxA_amplicons.fasta
└── amplicons_summary.csv
```
---
---
If a primer pair produces no amplicons, the corresponding FASTA file will be created but left empty.

---

### **4. Example FASTA Output Snippet**
gyrB_Strain2_hit3|len=510|fwd=CTTCGACATCGACGACGA|rev=ACGACGACGACTTCCAG
ATGGCTGATCGTACGATCGAGATCGACAGT...
CGTAGCTAGCTAGCTAGCTACGATCAGGAT...

### **Example CSV Output Snippet**
pair_id,fasta_file,sample_name,seq_id,fwd_start,fwd_end,rev_start,rev_end,fwd_mm,rev_mm,amp_len,forward_primer,reverse_primer
gyrB,genomes/Strain2.fna,Strain2,chromosome,8421,8440,8910,8931,1,2,510,CTTCGACATCGACGACGA,ACGACGACGACTTCCAG

---




# 🔬 InSilicoPCR – Parallel Multi-Genome In-Silico PCR Tool
### Developed by **Dr. Vishal**  
**Assistant Professor, School of Forensic Science  
National Forensic Sciences University (NFSU), Gandhinagar, India**

---

<p align="center">

<img src="https://img.shields.io/badge/Python-3.7+-blue.svg">
<img src="https://img.shields.io/badge/License-MIT-green.svg">
<img src="https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg">
<img src="https://img.shields.io/badge/Status-Active-success.svg">
<img src="https://img.shields.io/badge/Parallel-CPU%20Optimized-orange.svg">

</p>

---

# 📘 Overview

**InSilicoPCR** is a fast, parallelized tool for performing *in-silico PCR* across large collections of genomes using multiple primer pairs.  
It is built for forensic genomics, microbial identification, wildlife forensics, DNA barcoding, environmental DNA, and molecular diagnostics.

✔ Multi-FASTA output  
✔ Supports mismatches  
✔ Supports degenerate (IUPAC) primers  
✔ Auto-parallelization using all CPU cores  
✔ Zero external dependencies (pure Python)

---

# 📥 Installation

<details>
<summary><strong>📌 Click to expand Installation Instructions</strong></summary>

### **1️⃣ Clone the Repository**

```bash
git clone https://github.com/<your-username>/InSilicoPCR.git
cd InSilicoPCR
```

### **2️⃣ Make the tool executable**

```bash
chmod +x insilico_pcr_parallel_multifasta.py
```

### **3️⃣ Confirm Python version**

```bash
python3 --version
```

Must be **Python 3.7+**.

</details>

---

# 🚀 Usage Guide

<details open>
<summary><strong>📌 Basic Usage</strong></summary>

```bash
python insilico_pcr_parallel_multifasta.py     --fasta_dir genomes     --primers primers.csv     --out_dir output
```

</details>

<details>
<summary><strong>⚙️ Advanced Usage (Recommended)</strong></summary>

```bash
python insilico_pcr_parallel_multifasta.py     --fasta_dir genomes     --primers primers.csv     --out_dir output     --max_mismatch 2     --min_len 50     --max_len 2000     --workers 0
```

`workers=0` = use all available CPU cores  
`max_mismatch` = allowed mismatches per primer  

</details>

---

# 📂 Primer CSV Format

<details>
<summary><strong>Show primer format</strong></summary>

```
pair_id,forward,reverse
16S,AGAGTTTGATCCTGGCTCAG,GGTTACCTTGTTACGACTT
gyrB,CTTCGACATCGACGACGA,ACGACGACGACTTCCAG
```

✔ Column names can be flexible (`forward, fwd, left`...) – the tool auto-detects them.  
</details>

---

# 📁 FASTA Input Directory

<details>
<summary><strong>Show FASTA format</strong></summary>

```
genomes/
├── StrainA.fasta
├── Sample12.fa
└── environmental_01.fna
```

Allowed extensions: `.fa`, `.fasta`, `.fna`, `.ffn`

</details>

---

# 📦 Output Files

<details open>
<summary><strong>1️⃣ Multi-FASTA Files</strong></summary>

Each primer pair generates a multi-FASTA file:

```
16S_amplicons.fasta
gyrB_amplicons.fasta
```

FASTA header format:

```
>16S_Strain1_hit1|len=452|fwd=AGAGTTTGATCCTGGCTCAG|rev=GGTTACCTTGTTACGACTT
ATGCGTACGTT...
```

**Header fields:**

| Field | Meaning |
|-------|---------|
| `16S` | Primer pair ID |
| `Strain1` | Genome/sample name |
| `hit1` | First hit in this genome |
| `len=452` | Amplicon length |
| `fwd=` / `rev=` | Primer sequences |

</details>

---

<details open>
<summary><strong>2️⃣ Summary CSV (Metadata)</strong></summary>

A file named:

```
amplicons_summary.csv
```

Sample rows:

```
pair_id,fasta_file,sample_name,seq_id,fwd_start,fwd_end,rev_start,rev_end,fwd_mm,rev_mm,amp_len,forward_primer,reverse_primer
16S,genomes/Strain1.fasta,Strain1,contig0001,120,140,560,580,0,1,460,AGAGTTTGATCCTGGCTCAG,GGTTACCTTGTTACGACTT
```

✔ Ideal for forensic reports  
✔ Suitable for downstream bioinformatics pipelines  
✔ Contains primer mismatches & binding coordinates  

</details>

---

# 🧬 Output Workflow Diagram

```plaintext
        ┌───────────────────┐
        │   primers.csv     │
        └────────┬──────────┘
                 │
                 ▼
        ┌───────────────────┐
        │  InSilicoPCR Tool │
        └────────┬──────────┘
                 │
        ┌────────┴───────────────┐
        ▼                        ▼
┌─────────────────┐      ┌──────────────────────────┐
│ Multi-FASTA per │      │  amplicons_summary.csv    │
│   primer pair   │      │  (all detected amplicons) │
└─────────────────┘      └──────────────────────────┘
```

---

# 🧬 Forensic Applications

- Microbial forensic screening  
- Wildlife identification from trace DNA  
- Environmental DNA (eDNA) monitoring  
- Species differentiation using barcoding markers  
- STR/marker validation in forensic labs  
- PCR assay development & validation  
- Academic teaching and training modules  

---

# ⚙️ Performance Tips

- Place genomes on SSD storage  
- Set `--workers 0` to use all CPUs  
- Lower `--max_mismatch` for faster scanning  
- Increase `--min_len` to filter out noisy hits  
- Avoid extremely short or overly degenerate primers  

---

# 🧾 Citation

```
Vishal, V. (2025). InSilicoPCR: A Parallel Python Tool for Multi-Genome In-Silico PCR Analysis. 
School of Forensic Science, National Forensic Sciences University (NFSU), Gandhinagar, India.
GitHub: https://github.com/<your-username>/InSilicoPCR
```

---

# 📜 License

This project is released under the **MIT License**, the most appropriate for academic, forensic, and research tools.  
It permits free use, modification, and distribution with attribution.

---

# 🤝 Contributing

Pull requests are welcome.  
Feature ideas include:

- BLAST-gated validation  
- Web UI (Flask/Streamlit)  
- Nextflow or Snakemake automation  
- CRISPR-based primer support  

---


