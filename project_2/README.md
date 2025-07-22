## Directory Layout

advseq/
├── data/
│ ├── reads/ Raw FASTQ files
│ └── reference/ Reference FASTA
├── envs/
│ └── environment.yaml Conda spec (samtools, bowtie2)
├── index/ Bowtie2 index files
├── logs/ Per-rule log files
├── results/
│ ├── sam/ SAM outputs
│ ├── bam/ BAM files
│ ├── bam_sorted/ Sorted BAM + BAI
│ ├── stats/ Per-sample idxstats
│ └── idxstats_summary.tsv Aggregated stats
├── samples.tsv Sample table (sample, fq1, fq2)
├── workflow/
│ ├── Snakefile Main workflow
│ ├── config/
│ │ └── config.yaml Pipeline parameters
│ └── rules/
│ ├── bowtie.smk
│ └── samtools.smk
└── README.md This file


## Requirements

- **Conda** with strict channel priority  
- **Snakemake**

## Quick Start

1. Place FASTQ files in `data/reads/` and `reference.fa` in `data/reference/`.  
2. Edit `workflow/config/config.yaml` if needed.  
3. Populate `samples.tsv` in the project root:
   ```tsv
   sample    fq1                         fq2
   ERR024604 data/reads/ERR024604_1.fq.gz data/reads/ERR024604_2.fq.gz

## Run

snakemake \
  --snakefile workflow/Snakefile \
  --use-conda \
  --cores 4

## Outputs

All outputs will appear under results/
