# HGVS Variant Normalizer Pipeline (v0.5.0)

A Python-based bioinformatic pipeline designed for normalization, coordinate validation, and standardized annotation of biological genomic variants utilizing `biocommons/hgvs` and UTA.

## Features

- **Multi-variant Support:** Handles Single Nucleotide Variants (SNVs), Small Insertions/Deletions (Indels), Structural Variants (SVs), and Copy Number Variants (CNVs).
- **Strict Coordinate Validation:** Cross-validates variant representation and boundaries against biological databases.
- **Clinical Logic Integration:** Utilizes robust variant resolution tools.

## Installation & Setup

### Prerequisites
- Python 3.9+
- Docker (optional)

### Local Environment Setup
```bash
git clone [https://github.com/edademiral/hgvs-normalizer.git](https://github.com/edademiral/hgvs-normalizer.git)
cd hgvs-normalizer
pip install -r requirements.txt
