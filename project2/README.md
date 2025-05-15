# Project 1: SAM → BAM → Sorted → Indexed → Stats

This Snakemake workflow takes `.sam` files from a user-specified directory, converts them to `.bam`, sorts and indexes them, and finally runs `samtools idxstats` to produce mapping statistics.

``` bash
project2/
├── config
│   ├── config.yaml
│   └── samples.tsv
├── envs
│   └── env.yaml
├── logs
├── results
├── rules
│   ├── bowtie.smk
│   └── samtools.smk
├── data
│   ├── reads    ← (user-provided FASTQ files)
│   └── reference.fa    ← (user-provided reference)
├── run.sh
└── Snakefile

```

## Prerequisites

- WSL with Ubuntu (or similar)  
- Miniconda3 installed  
- Dev environment:

  \`\`\`bash
  cd ~/snakemake-projects
  conda env create -f environment.yml
  conda activate snakemake-dev
  \`\`\`

## Configuration

Edit \`config.yaml\` to point at your input folder:

\`\`\`yaml
sam_dir: "/home/loki/snakemake-projects/project1/data"
\`\`\`

## Usage

1. **Copy your SAMs** into \`project1/data/\`:

   \`\`\`bash
   cp /mnt/c/Users/Lenovo/Desktop/*.sam project1/data/
   \`\`\`

2. **Run the pipeline**:

   \`\`\`bash
   cd project1
   rm -rf .snakemake      # clean previous runs
   ./run.sh
   \`\`\`

   This will:
   - Build (via Conda) a project-local environment for \`samtools\`  
   - Convert, sort, index, and generate stats for all \`.sam\` files  
   - Output into \`results/bam/\`, \`results/bam_sorted/\`, and \`results/stats/\`

3. **Inspect outputs**:

   \`\`\`bash
   ls results/bam
   ls results/bam_sorted
   ls results/stats
   \`\`\`

4. **Target a single sample**:

   \`\`\`bash
   ./run.sh results/stats/ERR024604_tiny.txt
   \`\`\`

## Reproducibility

- Tool versions pinned in \`envs/samtools.yaml\`  
- Environments built locally under \`.conda_envs/\` via Conda  
- Single-script invocation ensures consistent behavior

## Cleaning Up

To remove all generated files:

\`\`\`bash
rm -rf .snakemake .conda_envs results
\`\`\`

---
