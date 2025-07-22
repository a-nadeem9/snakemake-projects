# workflow/rules/msa.smk

import pandas as pd

# Load sample names from the configuration
SAMPLES = pd.read_csv(config["samples_tsv"], sep=r"\s+", header=0)["sample"].tolist()

# Directories from config
ASSEMBLY_DIR = config["results"]["assembly_dir"]
MSA_DIR      = config["results"]["msa_dir"]

rule concat_scaffolds:
    input:
        expand(
            f"{ASSEMBLY_DIR}/{{sample}}/merged/merged.fasta",
            sample=SAMPLES
        )
    output:
        combined = f"{MSA_DIR}/scaffolds_all.fasta"
    run:
        shell("mkdir -p {MSA_DIR}")
        with open(output.combined, "w") as out_fh:
            for sample in SAMPLES:
                path = f"{ASSEMBLY_DIR}/{sample}/merged/merged.fasta"
                out_fh.write(open(path).read())

rule run_mafft:
    input:
        combined = f"{MSA_DIR}/scaffolds_all.fasta"
    output:
        aligned = f"{MSA_DIR}/aligned_scaffolds.fasta"
    params:
        threads = config["parameters"]["mafft"]["threads"]
    log:
        f"{config['logs']['mafft']}/mafft.log"
    conda:
        "../envs/mafft.yaml"
    shell:
        """
        mafft --thread {params.threads} {input.combined} \
            > {output.aligned} 2> {log}
        """
