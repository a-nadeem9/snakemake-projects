# Phylogenomics Snakemake Workflow

This Snakemake pipeline performs a complete phylogenomics analysis: from raw assembly and quality control, through multiple sequence alignment (MSA), to phylogenetic inference and variability statistics.

---

## 📁 Project Structure

```
.
├── config/
│   └── config.yaml              # Config: input paths, parameters, toggles
├── data/                        # Input raw and intermediate data
│   ├── reads/                   # Input genome assemblies
│   └── reference/               
├── workflow/
│   ├── Snakefile                # Main Snakemake entrypoint
│   └── rules/                   # Workflow rules
│       ├── assembly.smk         # Assembly rules (if applicable)
│       ├── qc.smk               # Quality control (e.g., QUAST, BUSCO)
│       ├── msa.smk              # Multiple sequence alignment (e.g., MAFFT)
│       ├── phylogeny.smk        # Tree inference (e.g., FastTree, RAxML)
│       └── variability.smk      # Variability stats (e.g., SNP density)
│   └──envs/ 
│       └── *.yaml
├── logs/                        # Rule-specific log files
├── results/                     # All pipeline outputs
│   ├── qc/
│   ├── msa/
│   ├── phylogeny/
│   └── variability/
└── README.md
```

---

## ⚙️ Configuration

Customize pipeline behavior in `config/config.yaml`. Example keys:

| Key             | Description                                    |
| --------------- | ---------------------------------------------- |
| `genome_dir`    | Path to input genomes                          |
| `metadata`      | Metadata file with sample IDs and descriptions |
| `threads`       | Number of threads per rule                     |
| `run_qc`        | `true/false` to enable/disable QC steps        |
| `run_phylogeny` | `true/false` to enable/disable tree inference  |
| `msa_tool`      | e.g., `mafft`, `clustalo`                      |
| `tree_tool`     | e.g., `fasttree`, `raxml`                      |

---

## ▶️ Running the Workflow

Activate the environment and execute:

```bash
snakemake \
  --use-conda \
  --cores 4 \
  -s workflow/Snakefile
```

To run specific steps:

```bash
snakemake results/phylogeny/tree.nwk
```

---

## Clean Start

To wipe all results and logs:

```bash
rm -rf results/ logs/ .snakemake/
```

---