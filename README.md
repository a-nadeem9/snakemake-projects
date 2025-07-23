# 🧬 Snakemake Projects

A collection of modular, reproducible Snakemake workflows for viral genomics, transcriptomics, phylogenomics, and metagenome analysis.


## 📂 Projects

| Folder       | Project Title                        | Description                                                                                                                     |
| ------------ | ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------- |
| `project_1/` | **SAMtools Alignment Processing**    | Converts SAM to BAM, sorts and indexes files, and outputs mapping statistics using `samtools idxstats`.                         |
| `project_2/` | **ADVSeq: RNA-Seq Mapping Workflow** | Minimal RNA-Seq alignment pipeline using Bowtie2 and Samtools with per-sample QC stats.                                         |
| `project_3/` | **RNA-Seq QC & Mapping**             | End-to-end RNA-Seq workflow including trimming, alignment, and MultiQC summarization.                                           |
| `project_4/` | **Phylogenomics Workflow**           | Genome-to-tree pipeline with MSA, tree inference, and SNP-based variability stats.                                              |
| `project_5/` | **Metagenome Assembly & Analysis**   | Full metagenome pipeline: QC, assembly (MEGAHIT), polishing (Racon), classification (Kraken2), and reference-guided refinement. |
| `project_6/` | **Transcript Quantification**        | Expression quantification using STAR, Kallisto, RSEM, featureCounts, and Trinity, with MultiQC reporting.                       |

---
