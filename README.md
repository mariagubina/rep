# Nextflow pipeline for detecting changes in allele-specific binding & expression events

This repository contains a reproducible pipeline for detecting allele-specific SNPs in RNA-seq, ChIP-seq, or ATAC-seq data. It includes two workflows:

1. **asymmetry.nf** — Identifies SNPs with significant allele-specific binding or expression.
2. **delta_asymmetry.nf** — Detects SNPs with changes in allele-specific signal between conditions.

---

## Requirements

- [Nextflow](https://www.nextflow.io/)
- Conda (https://anaconda.org/anaconda/conda)
- Reference genome (FASTA)
- HISAT2 genome index (GRCh38, genome_snp_tran)
- Metadata in CSV format

---

## Installation

```bash
git clone https://github.com/mariagubina/rsnps
cd rsnps
```

---

## Usage

### 1. Detect allele-specific SNPs

```bash
nextflow run pipeline.nf \
  --mode AS_SNPs \
  --metadata <metadata.csv> \
  --fasta <genome.fa> \
  -profile conda
  
```

### 2. Detect changes in allele-specific signal

```bash
nextflow run pipeline.nf \
  --mode asymmetry \
  --fasta <genome.fa> \
  -profile conda
```

