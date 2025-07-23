# workflow/rules/reference_pipeline.smk
# ------------------------------------------------------------
 
# ------------------------------------------------------------
# reference selection -> patch -> mapping consensus
# ------------------------------------------------------------
import os, json, subprocess

# --- read config & decide if the reference pipeline should run ----------
CAND_FASTA = config.get("candidates_fasta", "")           
PIPE_ON = (
    str(config.get("enable_reference_pipeline", "true")).lower() in {"true", "1", "yes"}
    and CAND_FASTA != ""                                  
)

# ---------- Mash sketch of candidate references ----------
if PIPE_ON:

    rule mash_sketch_candidates:
        input:
            fasta=CAND_FASTA
        output:
            sketch="resources/mash/candidates.msh"
        conda: "../envs/mash.yaml"
        shell: """
            mkdir -p resources/mash
            mash sketch -o {output.sketch} {input.fasta}
        """

    # ---------- Mash sketch of each polished assembly -----------
    rule mash_sketch_sample:
        input:
            fasta=f"{ASSEMBLY_DIR}/{{sample}}/final.fasta"
        output:
            sketch=temp(f"{ASSEMBLY_DIR}/{{sample}}/mash.msh")
        conda: "../envs/mash.yaml"
        shell: "mash sketch -o {output.sketch} {input.fasta}"

    # ---------- Pick best reference for each sample -------------
    rule select_best_reference:
        input:
            sample_sk=rules.mash_sketch_sample.output.sketch,
            cand_sk=rules.mash_sketch_candidates.output.sketch
        output:
            best_ref=f"{REFSEL_DIR}/{{sample}}_best_ref.txt"
        conda: "../envs/mash.yaml"
        shell: r"""
            mkdir -p {REFSEL_DIR}
            mash dist {input.sample_sk} {input.cand_sk} \
              | sort -k3,3n | head -n1 | cut -f2 | cut -d':' -f1 > {output.best_ref}
        """

    # ---------- RagTag patch with chosen reference --------------
    rule ragtag_patch:
        threads: config["parameters"]["ragtag"]["threads"]
        input:
            contigs = f"{ASSEMBLY_DIR}/{{sample}}/final.fasta",
            ref_txt = rules.select_best_reference.output.best_ref
        output:
            patched = f"{PATCH_DIR}/{{sample}}/patched.fasta"
        conda: "../envs/ragtag.yaml"
        log: "logs/ragtag/{sample}_patch.log"
        shell: r"""
            set -euo pipefail
            REF=$(cat {input.ref_txt})
            outdir={PATCH_DIR}/{wildcards.sample}
            mkdir -p "$outdir"

            # 1) run RagTag - minimap2 aligner
            ragtag.py patch --aligner mm2 -x asm10 "$REF" {input.contigs} -o "$outdir" 2> {log} || true

            # 2) find whatever file RagTag produced (version-dependent names)
            if   [[ -f "$outdir/ragtag.patch.fasta.gz" ]]; then gunzip -c "$outdir/ragtag.patch.fasta.gz" > {output.patched}
            elif [[ -f "$outdir/ragtag.patch.fa.gz"   ]]; then gunzip -c "$outdir/ragtag.patch.fa.gz"    > {output.patched}
            elif [[ -f "$outdir/ragtag.patch.fasta"   ]]; then cp       "$outdir/ragtag.patch.fasta"     {output.patched}
            elif [[ -f "$outdir/ragtag.patch.fa"      ]]; then cp       "$outdir/ragtag.patch.fa"        {output.patched}
            else
                echo "RagTag found no alignments – using un-patched assembly." >&2
                cp {input.contigs} {output.patched}
            fi
        """
    # ---------- BWA-MEM2 index of patched assembly --------------
    rule bwa_index_patched:
        input:
            fasta=rules.ragtag_patch.output.patched
        output:
            done=temp(f"{PATCH_DIR}/{{sample}}/index.done")
        conda: "../envs/bwa.yaml"
        shell: """
            bwa-mem2 index {input.fasta}
            touch {output.done}
        """

    # ---------- Map reads, call variants, create consensus ------
    READ_DIR = QC_DECONT_DIR if config["contaminants_fasta"] else QC_TRIM_DIR
    READ_TAG = "clean" if config["contaminants_fasta"] else "trimmed"

    rule consensus_from_mapping:
        threads: 8
        input:
            r1  = lambda wc: f"{READ_DIR}/{wc.sample}_1.{READ_TAG}.fastq.gz",
            r2  = lambda wc: f"{READ_DIR}/{wc.sample}_2.{READ_TAG}.fastq.gz",
            ref = rules.ragtag_patch.output.patched,
            idx = rules.bwa_index_patched.output.done
        output:
            vcf      = f"{CONS_DIR}/{{sample}}/{{sample}}.vcf",
            cons     = f"{CONS_DIR}/{{sample}}/consensus.fasta",
            vcf_gz   = temp(f"{CONS_DIR}/{{sample}}/tmp.vcf.gz")  # ← NEW: treat tmp VCF as temp()
        conda: "../envs/bwa.yaml"
        log:   "logs/consensus/{sample}.log"
        shell: r"""
            set -euo pipefail
            mkdir -p {CONS_DIR}/{wildcards.sample}

            ## 1. map reads -------------------------------------------------------
            bwa-mem2 mem -t {threads} {input.ref} {input.r1} {input.r2} \
                | samtools sort -@ {threads} -o aln.bam
            samtools index aln.bam

            ## 2. rebuild a guaranteed-fresh FASTA index -------------------------
            rm -f {input.ref}.fai          # remove any stale index
            samtools faidx {input.ref}

            ## 3. call variants ---------------------------------------------------
            bcftools mpileup -Ou -f {input.ref} aln.bam \
                | bcftools call -mv -Oz -o {output.vcf_gz}
            bcftools index {output.vcf_gz}
            bcftools view {output.vcf_gz} > {output.vcf}

            ## 4. build consensus -------------------------------------------------
            bcftools consensus -f {input.ref} {output.vcf_gz} > {output.cons} 2>>{log}

            ## 5. sanity-check ----------------------------------------------------
            if grep -q 'not found' {log}; then
                echo "ERROR: contig mismatch between VCF and reference – check RagTag/VCF steps." >&2
                exit 1
            fi
            
        """
else:
    # Pipeline disabled -> dummy rule so DAG resolves
    rule consensus_from_mapping:
        shell: "echo 'Reference pipeline disabled.'"
