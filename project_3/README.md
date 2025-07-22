# RNA-Seq QC & Mapping Snakemake Workflow

This Snakemake pipeline performs paired-end RNA-Seq quality control, trimming, alignment, and basic statistics. All results are summarized in a final **MultiQC** report. Each step can be toggled on or off via configuration flags.

---

## 📁 Project Structure

```
.
├── config/                  
│   └── config.yaml            # Paths, threads, adapters, skip flags
├── data/
│   ├── reads/                 # Raw FASTQ files
│   └── reference/
│       └── reference.fa       # Reference genome (FASTA)
├── resources/
│   └── adapters/
│       └── TruSeq3-PE.fa      # Trimmomatic adapter file
├── samples.tsv                # Sample table: sample, fq1, fq2
├── workflow/
│   ├── Snakefile              # Main Snakemake entrypoint
│   ├── envs/
│   │   └── environment.yaml   # Conda environment specification
│   └── rules/                 # Modular Snakemake rules
│       ├── bowtie.smk
│       ├── samtools.smk
│       └── qc.smk
├── logs/                      # Per-rule logs
├── results/                   # Output files
│   ├── sam/
│   ├── bam/
│   ├── bam_sorted/
│   ├── stats/
│   ├── trimmed/
│   └── qc/
│       ├── fastqc/{raw,trimmed}/
│       ├── qualimap/
│       └── multiqc_report.html
└── README.Rmd
```

---

## ⚙️ Configuration

Edit `config/config.yaml` to set up paths and pipeline options:

| Key             | Description                      |
| --------------- | -------------------------------- |
| `samples`       | Path to `samples.tsv`            |
| `reference`     | Path to reference genome (FASTA) |
| `threads`       | Number of cores to use per step  |
| `trim.adapter`  | Trimmomatic adapter file path    |
| `trim.options`  | Trimming options                 |
| `bowtie.minins` | Minimum insert size              |
| `bowtie.maxins` | Maximum insert size              |
| `skip_trimming` | `true` to skip trimming step     |
| `skip_fastqc`   | `true` to skip FastQC            |
| `skip_qualimap` | `true` to skip Qualimap QC       |

> **Note:**
> `samples.tsv` must contain the following header:
> `sample	fq1	fq2`

---

## ▶️ Running the Workflow

From the project root, run:

```bash
snakemake \
  --use-conda \
  -j 4 \
  -s workflow/Snakefile
```

---

## Clean Start (Optional)

To clear previous outputs and start fresh:

```bash
rm -rf results/ logs/ index/ .snakemake/
mkdir -p results/{sam,bam,bam_sorted,stats,trimmed,qc/fastqc/{raw,trimmed},qc/qualimap} logs index
snakemake --use-conda -j 4 -s workflow/Snakefile
```
