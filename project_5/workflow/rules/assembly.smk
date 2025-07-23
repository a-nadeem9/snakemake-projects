import os

# pick correct read set depending on de-contam flag
USE_DECONT = bool(config["contaminants_fasta"])
READ_DIR   = QC_DECONT_DIR if USE_DECONT else QC_TRIM_DIR
READ_TAG   = "clean"       if USE_DECONT else "trimmed"

# ---------------- SPAdes ----------------
rule spades_assemble:
    threads: config["parameters"]["spades"]["threads"]
    input:
        r1 = lambda wc: f"{READ_DIR}/{wc.sample}_1.{READ_TAG}.fastq.gz",
        r2 = lambda wc: f"{READ_DIR}/{wc.sample}_2.{READ_TAG}.fastq.gz",
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/contigs/contigs.fasta",
    params:
        mem = config["parameters"]["spades"]["memory"]
    conda: "../envs/spades.yaml"
    log:   "logs/spades/{sample}.log"
    shell: r"""
        spades.py \
          --phred-offset 33 \
          -1 {input.r1} -2 {input.r2} \
          -o {ASSEMBLY_DIR}/{wildcards.sample}/contigs \
          -t {threads} -m {params.mem} \
          2> {log}
    """

# ---------------- RagTag scaffolding ----------------
rule ragtag_scaffold:
    threads: config["parameters"]["ragtag"]["threads"]
    input:
        contigs = rules.spades_assemble.output.contigs,
        ref     = config["reference"],
    output:
        scaffolds = f"{ASSEMBLY_DIR}/{{sample}}/scaffolds.fasta",
    conda: "../envs/ragtag.yaml"
    log:   "logs/ragtag/{sample}.log"
    shell: r"""
        ragtag.py scaffold {input.ref} {input.contigs} \
               -o {ASSEMBLY_DIR}/{wildcards.sample}/ragtag_out 2> {log}
        mv {ASSEMBLY_DIR}/{wildcards.sample}/ragtag_out/ragtag.scaffold.fasta {output.scaffolds}
    """

# ---------------- extract main scaffold & rename ----------------
rule extract_ragtag_scaffold:
    input:  scaff = rules.ragtag_scaffold.output.scaffolds
    output: merged = f"{ASSEMBLY_DIR}/{{sample}}/merged/merged.fasta"
    conda:  "../envs/seqkit.yaml"
    shell:  r"""
        seqkit grep -r -p "_RagTag$" {input.scaff} \
        | sed 's/^>.*/>{wildcards.sample}/' > {output.merged}
    """
