# Metagenome Assembly & Analysis Snakemake Workflow

This Snakemake pipeline performs end-to-end metagenome analysis—from read QC, assembly, and polishing, to taxonomic classification, contaminant removal, phylogenetic analysis, and final reporting. Modular rule files support flexible, reproducible processing across diverse metagenomic datasets.

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
├── config/                      # Configuration directory
│   └── config.yaml              # Paths, toggles, and parameters
├── data/                        # Input raw reads or references
│   ├── reads/                   # Input FASTQ files
│   └── references/              # Optional reference genomes
├── workflow/
│   ├── Snakefile                # Main Snakemake entrypoint
│   └── rules/                   # Workflow rules
│       ├── qc.smk               # Read QC (fastqc, trimmomatic)
│       ├── assembly.smk         # Assembly (megahit)
│       ├── polish.smk           # Polishing (racon)
│       ├── kraken.smk           # Taxonomic classification (kraken2)
│       ├── decontam.smk         # Decontamination
│       ├── candidates.smk       # Genome bin candidate filtering
│       ├── msa.smk              # Multiple sequence alignment
│       ├── phylogeny.smk        # Phylogenetic inference
│       ├── variability.smk      # Diversity and SNP metrics
│       ├── multiqc.smk          # Aggregated reporting
│       └── reference_pipeline.smk # For reference-based workflows
├── envs/                        # Conda environments
│   └── *.yaml                   # Rule-specific envs
├── logs/                        # Log files by rule
├── results/                     # All output files
│   ├── qc/
│   ├── assembly/
│   ├── polished/
│   ├── kraken/
│   ├── decontam/
│   ├── msa/
│   ├── phylogeny/
│   ├── variability/
│   └── multiqc/
└── README.md
```

---

## ⚙️ Configuration

Edit `config/config.yaml` to define input paths, tools, and control flags. Common keys:

| Key             | Description                           |
| --------------- | ------------------------------------- |
| `samples`       | Sample sheet with IDs and FASTQ paths |
| `threads`       | Threads to use per rule               |
| `run_decontam`  | Enable/disable contaminant filtering  |
| `run_phylogeny` | Enable/disable phylogenetic analysis  |
| `qc_tools`      | Tools used: `fastqc`, `trimmomatic`   |
| `assembly_tool` | Set to: `megahit`                     |
| `polish_tool`   | Set to: `racon`                       |
| `classifier`    | Set to: `kraken2`                     |

---

## 📌 Reference Genome Selection (`candidates.fasta`)

For reference-guided assembly and consensus calling, the workflow uses a `candidates.fasta` file containing full genomes from eight Human Adenovirus strains. These genomes were downloaded from NCBI RefSeq and selected to represent diversity across species:

* Human Adenovirus A
* Human Adenovirus B
* Human Adenovirus C
* Human Adenovirus D
* Human Adenovirus E
* Human Adenovirus F
* Human Adenovirus 1
* Human Adenovirus 2

The reference selection pipeline works as follows:

1. Each assembled sample is compared to all references using **Mash**.
2. The most similar reference (based on Mash distance) is selected.
3. **RagTag** patches the assembly using this selected reference.
4. Reads are aligned to the patched assembly and a consensus is built.

---

## ▶️ Running the Workflow

```bash
snakemake \
  --use-conda \
  -s workflow/Snakefile \
  --cores 8
```

To run a specific target:

```bash
snakemake results/phylogeny/tree.nwk
```

---

## Clean Slate

Remove intermediate files before rerunning:

```bash
rm -rf results/ logs/ .snakemake/
```

---
