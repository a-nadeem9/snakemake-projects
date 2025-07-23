# ------------------------------------------------------------
# candidates.smk – minimal, assumes candidates.fasta exists
# ------------------------------------------------------------
CAND_FASTA = "resources/reference/candidates.fasta"
MASH_DIR   = "resources/mash"

rule build_mash_sketch:
    input:
        CAND_FASTA
    output:
        f"{MASH_DIR}/candidates.msh"
    conda:
        "envs/mash.yaml"
    shell:
        """
        mkdir -p {MASH_DIR}
        mash sketch -o {output.rsplit('.',1)[0]} {input}
        """
