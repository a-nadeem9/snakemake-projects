# SAMtools Alignment Processing Workflow

This lightweight Snakemake pipeline automates downstream processing of SAM alignment files. It converts them to BAM, sorts and indexes them, and produces per-sample mapping statistics using `samtools idxstats`. Simple, fast, and ideal for quality control or downstream quantification.

---

## 📁 Directory Layout

```
project1/
├── Snakefile               # Main Snakemake workflow
├── envs/
│   └── samtools.yaml       # Conda environment for samtools
├── data/                   # Input SAM files
├── results/                # Pipeline outputs
│   ├── bam/                # BAM files
│   ├── bam_sorted/         # Sorted BAM + BAI index
│   └── stats/              # idxstats per sample
├── .conda_envs/            # Auto-generated Conda environments
└── run.sh                  # Run wrapper script
```

---

## Quick Start

1. Copy your `.sam` files into the `data/` folder:

   ```bash
   cp /mnt/c/Users/Lenovo/Desktop/*.sam project1/data/
   ```

2. Create and activate the environment:

   ```bash
   conda env create -f environment.yml
   conda activate snakemake-dev
   ```

3. Run the pipeline:

   ```bash
   cd project1
   ./run.sh
   ```

4. (Optional) Run a specific target:

   ```bash
   ./run.sh results/stats/ERR024604_tiny.txt
   ```

---

## Output Summary

All outputs are saved in the `results/` directory:

* `bam/`: BAM-converted alignments
* `bam_sorted/`: Sorted BAMs with `.bai` index
* `stats/`: Mapping statistics from `samtools idxstats`
