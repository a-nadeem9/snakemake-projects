import pandas as pd
from pathlib import Path

# Load sample names from the configuration
SAMPLES = pd.read_csv(config["samples_tsv"], sep=r"\s+", header=0)["sample"].tolist()

# Define directories from config
ASSEMBLY_DIR = config["results"]["assembly_dir"]
REF = config["reference"]

rule spades_assemble:
    input:
        r1 = lambda wc: f"{config['results']['qc_trimmed']}/{wc.sample}_1.trimmed.fastq.gz",
        r2 = lambda wc: f"{config['results']['qc_trimmed']}/{wc.sample}_2.trimmed.fastq.gz"
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/contigs/contigs.fasta"
    params:
        threads = config["parameters"]["spades"]["threads"],
        memory = config["parameters"]["spades"]["memory"]
    conda:
        "../envs/spades.yaml"
    log:
        f"{config['logs']['spades']}/{{sample}}.log"
    shell:
        """
        mkdir -p {ASSEMBLY_DIR}/{wildcards.sample}/contigs
        spades.py \
          --phred-offset 33 \
          -1 {input.r1} \
          -2 {input.r2} \
          -o {ASSEMBLY_DIR}/{wildcards.sample}/contigs \
          -t {params.threads} \
          -m {params.memory} \
          > {log} 2>&1
        """

rule ragtag_scaffold:
    input:
        contigs = lambda wc: f"{ASSEMBLY_DIR}/{wc.sample}/contigs/contigs.fasta",
        ref = REF
    output:
        scaffolds = f"{ASSEMBLY_DIR}/{{sample}}/scaffolds.fasta"
    conda:
        "../envs/ragtag.yaml"
    log:
        f"{config['logs']['ragtag']}/{{sample}}.log"
    shell:
        """
        mkdir -p {ASSEMBLY_DIR}/{wildcards.sample}/ragtag_output
        ragtag.py scaffold {input.ref} {input.contigs} \
            -o {ASSEMBLY_DIR}/{wildcards.sample}/ragtag_output \
            > {log} 2>&1
        mv {ASSEMBLY_DIR}/{wildcards.sample}/ragtag_output/ragtag.scaffold.fasta {output.scaffolds}
        """

rule pick_longest_scaffold:
    input:
        scaff = lambda wc: f"{config['results']['assembly_dir']}/{wc.sample}/scaffolds.fasta"
    output:
        merged = f"{config['results']['assembly_dir']}/{{sample}}/merged/merged.fasta"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        mkdir -p $(dirname {output.merged})
        seqkit fx2tab {input.scaff} \
          | awk -F'\t' 'BEGIN{{max=0}}{{if(length($2)>max){{max=length($2); header=$1; seq=$2}}}}END{{print ">{wildcards.sample}\\n"seq}}' \
          > {output.merged}
        """
