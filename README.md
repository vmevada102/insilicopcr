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
git clone https://github.com/vmevada102/InSilicoPCR.git
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
will be added later
```

---



---

# 🤝 Contributing

Pull requests are welcome.  
Feature ideas include:

- BLAST-gated validation  
- Web UI (Flask/Streamlit)  
- Nextflow or Snakemake automation  
- CRISPR-based primer support  

---


