# ADVSeq: Bowtie2 + Samtools RNA-Seq Mapping Workflow

This Snakemake workflow performs basic RNA-Seq mapping using **Bowtie2** and **Samtools**. It takes in paired-end FASTQ files, maps them to a reference genome, converts and sorts the alignments, and computes per-sample statistics. It is intended as a minimal, reproducible pipeline for quick quality assessments and summaries of RNA-Seq alignment results.

---

## 📁 Directory Layout

```
avdseq/
├── data/                         # Input files
│   ├── reads/                    # Raw FASTQ files
│   └── reference/                # Reference FASTA
├── envs/
│   └── environment.yaml          # Conda spec (e.g., samtools, bowtie2)
├── index/                        # Bowtie2 index files
├── logs/                         # Per-rule log files
├── results/                      # Output directory
│   ├── sam/                      # Uncompressed SAM files
│   ├── bam/                      # Converted BAM files
│   ├── bam_sorted/               # Sorted BAM + BAI index
│   ├── stats/                    # Per-sample idxstats
│   └── idxstats_summary.tsv      # Aggregated idxstats summary
├── samples.tsv                   # Sample table (sample, fq1, fq2)
├── workflow/
│   ├── Snakefile                 # Main Snakemake entrypoint
│   ├── config/
│   │   └── config.yaml           # Pipeline parameters
│   └── rules/
│       ├── bowtie.smk            # Bowtie2 alignment rules
│       └── samtools.smk          # Samtools conversion and stats
└── README.md                     # This file
```

---

## 🛠️ Requirements

* **Conda** (with strict channel priority)
* **Snakemake**

---

## Quick Start

1. Place FASTQ files into `data/reads/`, and the reference genome (`reference.fa`) into `data/reference/`.
2. Edit configuration in `workflow/config/config.yaml` as needed.
3. Create a `samples.tsv` file in the project root. Example:

   ```tsv
   sample    fq1                          fq2
   ERR024604 data/reads/ERR024604_1.fq.gz data/reads/ERR024604_2.fq.gz
   ```

---

## ▶️ Run the Workflow

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --use-conda \
  --cores 4
```

---

## Output Summary

All pipeline results will be saved in the `results/` directory, including:

* `sam/`: Raw alignments (SAM format)
* `bam/`: Converted BAM files
* `bam_sorted/`: Sorted BAM files with BAI index
* `stats/`: `samtools idxstats` outputs per sample
* `idxstats_summary.tsv`: Merged summary table across all samples
