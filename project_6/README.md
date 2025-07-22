# Transcript Quantification Snakemake Workflow

This Snakemake workflow performs transcript-level expression analysis using multiple quantification strategies including alignment-based and alignment-free tools. It supports STAR, Kallisto, RSEM, featureCounts, and Trinity, and outputs a unified MultiQC report.

---

## Quick Start

```bash
conda config --set channel_priority flexible
snakemake --use-conda --cores all
```

---

## 📁 Project Structure

```
.
├── config/                      # Configuration files
│   └── config.yaml              # Sample info, parameters, toggles
├── data/                        # Input FASTQ and annotations
│   ├── reads/                   # FASTQ files
│   └── annotation/              # GTF/GFF and reference FASTA
├── workflow/
│   ├── Snakefile                # Main workflow entrypoint
│   └── rules/                   # Individual modules
│       ├── preprocess_qc.smk    # Trimming, QC (e.g., FastQC)
│       ├── star_index.smk       # STAR genome index
│       ├── star_align.smk       # STAR alignment
│       ├── kallisto.smk         # Kallisto quantification
│       ├── rsem.smk             # RSEM alignment & quantification
│       ├── featurecounts.smk    # featureCounts gene quantification
│       ├── merge_samples.smk    # Merging outputs
│       ├── multiqc.smk          # Summary report
│       └── trinity.smk          # De novo transcriptome assembly
├── envs/                        # Conda environment YAMLs
├── logs/                        # Logs for each rule
├── results/                     # Outputs
│   ├── qc/
│   ├── star/
│   ├── kallisto/
│   ├── rsem/
│   ├── featurecounts/
│   ├── trinity/
│   ├── merged/
│   └── multiqc/
└── README.md
```

---

## ⚙️ Configuration

Edit `config/config.yaml` with paths and options:

| Key                 | Description                          |
| ------------------- | ------------------------------------ |
| `samples`           | Sample sheet: sample ID, FASTQ files |
| `threads`           | Number of threads per rule           |
| `gtf`               | Path to GTF file                     |
| `genome_fasta`      | Reference genome FASTA               |
| `run_star`          | Enable STAR alignment                |
| `run_kallisto`      | Enable Kallisto quantification       |
| `run_rsem`          | Enable RSEM analysis                 |
| `run_featurecounts` | Enable featureCounts step            |
| `run_trinity`       | Enable Trinity assembly (optional)   |

---

## ▶️ Running the Workflow

```bash
snakemake \
  --use-conda \
  -s workflow/Snakefile \
  --cores 8
```

To run a specific rule:

```bash
snakemake results/kallisto/sample1/abundance.tsv
```

---

## Clean Run

```bash
rm -rf results/ logs/ .snakemake/
```

---
