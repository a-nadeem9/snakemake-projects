# --- at top of workflow/rules/polish.smk ---
def _as_bool(val):
    return str(val).lower() in {"true", "yes", "1"}

ENABLE_POLISH = _as_bool(config.get("enable_polishing", "true"))
# workflow/rules/polish.smk
# ------------------------------------------------------------
# Optional short-read assembly polishing with Pilon
# ------------------------------------------------------------
import os

ENABLE_POLISH = bool(config.get("enable_polishing", True))

# ---------- constants ------------------------------------------------
MERGED_SCAF      = f"{ASSEMBLY_DIR}" + "/{sample}/merged/merged.fasta"
PILON_FOLDER     = f"{ASSEMBLY_DIR}" + "/{sample}/pilon"
BOWTIE2_IDX_PREF = PILON_FOLDER + "/idx"          # *.1.bt2 etc.
IDX_FLAG         = PILON_FOLDER + "/index.done"   # sentinel file

# choose read set depending on de-contam flag
READ_DIR = QC_DECONT_DIR if config["contaminants_fasta"] else QC_TRIM_DIR
READ_TAG = "clean"       if config["contaminants_fasta"] else "trimmed"

# ==================== 1. build Bowtie2 index =========================
if ENABLE_POLISH:

    rule build_bowtie2_index_pilon:
        input:
            fasta = MERGED_SCAF
        output:
            flag  = temp(IDX_FLAG)
        params:
            prefix = BOWTIE2_IDX_PREF
        conda: "../envs/bowtie2.yaml"
        shell: r"""
            bowtie2-build {input.fasta} {params.prefix} 2> /dev/null
            touch {output.flag}
        """

    # ================= 2. map reads back to assembly =================
    rule map_reads_for_pilon:
        threads: 8
        input:
            r1  = lambda wc: f"{READ_DIR}/{wc.sample}_1.{READ_TAG}.fastq.gz",
            r2  = lambda wc: f"{READ_DIR}/{wc.sample}_2.{READ_TAG}.fastq.gz",
            idx = IDX_FLAG,
        output:
            bam = temp(PILON_FOLDER + "/aln.sorted.bam")
        params:
            prefix = BOWTIE2_IDX_PREF
        conda: "../envs/bowtie2.yaml"
        shell: r"""
            bowtie2 --very-sensitive --threads {threads} \
                    -x {params.prefix} \
                    -1 {input.r1} -2 {input.r2} 2> /dev/null \
            | samtools sort -@ {threads} -o {output.bam}
            samtools index {output.bam}
        """

    # ===================== 3. run Pilon ==============================
    rule run_pilon:
        threads: config["parameters"]["pilon"]["threads"]
        input:
            bam = rules.map_reads_for_pilon.output.bam,
            ref = MERGED_SCAF
        output:
            polished = PILON_FOLDER + "/polished.fasta"
        conda: "../envs/pilon.yaml"
        log:   "logs/pilon/{sample}.log"
        shell: r"""
            pilon --genome {input.ref} --frags {input.bam} \
                  --threads {threads} --output pilon_tmp 2> {log}
            mv pilon_tmp.fasta {output.polished}
        """

    FINAL_FASTA = PILON_FOLDER + "/polished.fasta"

else:
    # polishing disabled → just use merged scaffold
    FINAL_FASTA = MERGED_SCAF

# ===================== 4. final symlink ==============================
rule finalize_assembly:
    input:
        fasta = FINAL_FASTA
    output:
        final = f"{ASSEMBLY_DIR}" + "/{sample}/final.fasta"
    shell:
        "ln -sf $(realpath {input.fasta}) {output.final}"
