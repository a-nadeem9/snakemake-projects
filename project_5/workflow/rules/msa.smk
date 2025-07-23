# workflow/rules/msa.smk
# ------------------------------------------------------------
# Build a multi-sequence alignment of per-sample FINAL assemblies
# ------------------------------------------------------------
# Globals exported from Snakefile:
#   SAMPLES, ASSEMBLY_DIR, MSA_DIR

SCAFF_PATTERN = f"{ASSEMBLY_DIR}/{{sample}}/final.fasta"

# ------------------------------------------------------------
rule concat_scaffolds:
    input:
        expand(SCAFF_PATTERN, sample=SAMPLES)
    output:
        combined = f"{MSA_DIR}/scaffolds_all.fasta"
    shell:
        """
        cat {input} > {output.combined}
        """

# ------------------------------------------------------------
rule run_mafft:
    threads: config["parameters"]["mafft"]["threads"]
    input:
        combined = rules.concat_scaffolds.output.combined
    output:
        aligned = f"{MSA_DIR}/aligned_scaffolds.fasta"
    log:
        "logs/mafft/mafft.log"
    conda: "../envs/mafft.yaml"
    shell:
        """
        mafft --thread {threads} {input.combined} \
              > {output.aligned} 2> {log}
        """
