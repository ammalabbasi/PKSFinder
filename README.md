# 🧬 PKSFinder

**PKSFinder** detects and profiles **polyketide synthase (PKS)** gene signatures from **unmapped sequencing reads**. It identifies PKS presence, estimates abundance, assigns a microbial species origin, and provides a confidence score for each hit.

---

## 🔧 Features

- ✅ Accepts **unmapped FASTQ/BAM** reads
- 🔍 Detects **PKS domains** using HMMs or alignments
- 🧬 Assigns **species origin** for each PKS hit
- 📊 Reports **abundance** and **confidence scores**
- 📁 Outputs include tables, sequences, and taxonomic reports

---

## 🚀 Quick Start

```bash
git clone https://github.com/yourusername/PKSFinder.git
cd PKSFinder
pip install -r requirements.txt

python pksfinder.py \
    --input unmapped_reads.fastq \
    --output results/ \
    --threads 8
