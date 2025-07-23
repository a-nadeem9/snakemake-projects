# ------------------------------------------------------------
#   OPTIONAL contaminant removal
#   (enabled when config["contaminants_fasta"] is *not* empty)
# ------------------------------------------------------------

import os

CONTAM_FASTA = config["contaminants_fasta"]

if CONTAM_FASTA:          # ---------- real de-contamination ----------

    IDX_DIR    = "resources/contaminants/bt2_idx"
    IDX_PREFIX = f"{IDX_DIR}/contaminants"

    rule build_contaminant_index:
        input:  fasta = CONTAM_FASTA
        output: f"{IDX_DIR}/index.done"
        params: prefix = IDX_PREFIX
        conda:  "../envs/bowtie2.yaml"
        log:    "logs/decont/bowtie2_build.log"
        shell:  r"""
            bowtie2-build {input.fasta} {params.prefix} 2> {log}
            touch {output}
        """

    rule remove_contaminants:
        threads: 8
        input:
            fq1 = f"{QC_TRIM_DIR}/{{sample}}_1.trimmed.fastq.gz",
            fq2 = f"{QC_TRIM_DIR}/{{sample}}_2.trimmed.fastq.gz",
            idx = rules.build_contaminant_index.output,
        output:
            clean1 = f"{QC_DECONT_DIR}/{{sample}}_1.clean.fastq.gz",
            clean2 = f"{QC_DECONT_DIR}/{{sample}}_2.clean.fastq.gz",
        params:
            prefix = IDX_PREFIX
        conda: "../envs/bowtie2.yaml"
        log:   "logs/decont/{sample}.log"
        shell: r"""
            bowtie2 --very-sensitive --threads {threads} \
                    -x {params.prefix} -1 {input.fq1} -2 {input.fq2} \
                    2> {log} \
            | samtools fastq -@ {threads} -f 12 -F 256 \
                  -1 {output.clean1} -2 {output.clean2} -
        """

else:                      # ---------- de-contamination skipped ----------
    # Nothing to do – downstream rules will just use the trimmed reads.
    pass
