# CRISPR Sanger Sequencing Analysis Toolkit

A comprehensive toolkit for analyzing CRISPR-induced mutations from Sanger sequencing data (.ab1 files). Supports detection of insertions, deletions, and complex compound heterozygous mutations.

## Overview

This toolkit provides automated analysis of Sanger sequencing traces from CRISPR-edited samples. It uses local alignment algorithms and chromatogram peak pattern analysis to:

1. **Local Alignment**: Find the best matching region between sequencing and reference sequences
2. **Indel Detection**: Identify insertions and deletions at target sites
3. **Heterozygosity Verification**: Distinguish homozygous vs heterozygous mutations using peak patterns
4. **Compound Heterozygote Detection**: Identify samples with different mutations on each allele

## Table of Contents

1. [Installation](#installation)
2. [Input Files](#input-files)
3. [Basic Usage](#basic-usage)
4. [Algorithm Details](#algorithm-details)
5. [Output Interpretation](#output-interpretation)
6. [Troubleshooting](#troubleshooting)

---

## 1. Installation

### 1.1 Requirements

| Software | Version | Required |
|----------|---------|----------|
| Python | 3.8+ | Yes |
| Biopython | 1.78+ | Yes |
| NumPy | 1.20+ | Yes |
| Matplotlib | 3.3+ | Yes |
| Biopython (for ABI files) | - | Yes |

### 1.2 Install Dependencies

```bash
# Using pip
pip install numpy biopython matplotlib

# Or using conda
conda install numpy biopython matplotlib
```

### 1.3 Download the Tool

```bash
# Clone the repository
git clone https://github.com/IceBear321/CRISPR_Sanger_Analysis.git
cd CRISPR_Sanger_Analysis
```

---

## 2. Input Files

### 2.1 Reference Sequence (FASTA)

Prepare a reference sequence file in FASTA format:

```
>demo
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
```

**Example Input Files:**
```
example/
├── demo.fasta              # Reference sequence
├── test1.ab1              # Sample 1
├── test2.ab1              # Sample 2
└── test3.ab1              # Sample 3
```

Sanger sequencing files in ABI format (.ab1). These contain:
- Base calls
- Chromatogram data (G/A/T/C channels)
- Quality scores
- Peak location information

---

## 3. Basic Usage

### 3.1 Single File Analysis

```bash
python3 crispr_final_analyzer.py \
    -i example/test1.ab1 \
    -r example/demo.fasta \
    -p 250 \
    -o my_analysis
```

Parameters:
- `-i`: Input .ab1 file
- `-r`: Reference sequence (FASTA)
- `-p`: Target position in reference (1-based)
- `-o`: Output prefix

### 3.2 Batch Analysis

Process multiple files in a directory:

```bash
python3 crispr_final_analyzer.py \
    -i example/ \
    -r example/demo.fasta \
    -p 250 \
    -o batch_results
```

### 3.3 With Visualization

Generate chromatogram visualization:

```bash
python3 crispr_final_analyzer.py \
    -i example/test1.ab1 \
    -r example/demo.fasta \
    -p 250 \
    -o my_analysis \
    --visualize
```

---

## 4. Algorithm Details

### 4.1 Local Alignment Strategy

The tool uses **local alignment** (Smith-Waterman algorithm) to find the best matching region between the sequencing read and reference:

```python
aligner = Align.PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -2
```

**Why Local Alignment?**
- Handles poor quality at sequencing primers
- Ignores low-quality start/end regions
- Focuses on the core matching region

### 4.2 Coordinate Mapping

After alignment, the tool builds a coordinate mapping:
```
Reference Position → Sequencing Position
```

This allows precise identification of where mutations occur relative to the reference.

### 4.3 Peak Pattern Analysis

For each position, the tool calculates the **secondary/primary peak ratio**:

```
Ratio = Secondary Peak Height / Primary Peak Height
```

**Interpretation:**
- **Ratio < 0.3**: Homozygous (single allele)
- **Ratio > 0.3**: Heterozygous (two alleles)

### 4.4 Frameshift Detection

After detecting a heterozygous position, the tool analyzes downstream positions:

- If subsequent positions show high ratios (>0.5), it indicates a **frameshift**
- Frameshift suggests different mutations on each allele (compound heterozygote)

### 4.5 Mutation Classification

| Pattern | Genotype | Description |
|---------|-----------|-------------|
| Low ratio throughout | Wild Type / Homozygous | Single allele |
| High ratio at target | Heterozygous | Single substitution |
| High ratio with downstream effect | Compound Heterozygote | Different indels on each allele |

---

## 5. Output Interpretation

### 5.1 Text Report

The tool generates a detailed report:

```
File: sample.ab1
========================================
Genotype: Compound Heterozygous
Mutation site: Reference 250 / Sequencing 186
Reference base: C
Peak analysis: Primary peak C (height: 2450), Secondary peak A (height: 1580), Ratio 0.645
Frameshift: Yes (average ratio for next 20 positions: 0.520)

--- Inferred Results ---
Allele 1: 1bp C insertion at reference position 248 (C->CC)
Allele 2: C to A mutation at reference position 248 (C->A)
```

### 5.2 CSV Summary

For batch analysis, a CSV file is generated with:
- Sample name
- Genotype
- Mutation type
- Position
- Peak ratio

### 5.3 Visualization (Optional)

When `--visualize` is enabled, generates:
- Chromatogram plot with base calls
- Secondary peak ratio analysis
- Mutation site highlighting

---

## 6. Troubleshooting

### 6.1 No Alignment Found

**Problem**: "No alignment found"

**Solutions**:
1. Check that reference and sequencing are in the same orientation
2. Verify the target position is within the sequencing range
3. Ensure reference sequence is correct

### 6.2 Position Not in Alignment Range

**Problem**: "Reference position not in alignment range"

**Solutions**:
1. The sequencing may not cover the target region
2. Check sequencing quality at that position
3. Try a different sequencing direction

### 6.3 Incorrect Genotype

**Problem**: Wild type called as mutant or vice versa

**Solutions**:
1. Adjust the heterozygous threshold (default: 0.3)
2. Check chromatogram quality
3. Verify reference sequence is correct
4. Try visualization to inspect the data manually

### 6.4 Low Quality Sequencing

**Problem**: High background noise in chromatogram

**Solutions**:
1. Re-sequence with better quality
2. Use local alignment (already default)
3. Check for contamination

---

## Example Analysis Result

For a CRISPR-edited sample with compound heterozygous mutations:

| Metric | Value |
|--------|-------|
| Sample | test1.ab1 |
| Genotype | Compound Heterozygous |
| Indel 1 | +1bp (C insertion) at position 260 |
| Indel 2 | -1bp (G deletion) at position 402 |
| Secondary Peak Ratio | 0.64 (first indel), 0.60 (second) |

This indicates two different mutations occurred on each chromosome - a complex but valid CRISPR outcome.

---

## Citation

If you use this toolkit, please cite:
- Biopython: Cock PJA et al. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics
- Original CRISPR analysis methodology

---

## License

MIT License

## Author

Zhihe Cai

## Contact

zh.cai@pku.edu.cn
