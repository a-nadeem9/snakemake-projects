# workflow/rules/qc.smk

# FastQC on raw reads using Snakemake wrapper
rule fastqc_raw_1:
    input:
        fq=lambda wc: samples.at[wc.sample, "fq1"]
    output:
        html="results/qc/fastqc/raw/{sample}_1_fastqc.html",
        zip="results/qc/fastqc/raw/{sample}_1_fastqc.zip"
    log:
        "logs/fastqc/raw/{sample}_1.log"
    threads: config["threads"]["qc"]
    wrapper:
        "0.72.0/bio/fastqc"

rule fastqc_raw_2:
    input:
        fq=lambda wc: samples.at[wc.sample, "fq2"]
    output:
        html="results/qc/fastqc/raw/{sample}_2_fastqc.html",
        zip="results/qc/fastqc/raw/{sample}_2_fastqc.zip"
    log:
        "logs/fastqc/raw/{sample}_2.log"
    threads: config["threads"]["qc"]
    wrapper:
        "0.72.0/bio/fastqc"

# Trimmomatic paired-end trimming
rule trim_reads:
    input:
        fq1=lambda wc: samples.at[wc.sample, "fq1"],
        fq2=lambda wc: samples.at[wc.sample, "fq2"]
    output:
        trimmed_r1  = "results/trimmed/{sample}_1.trimmed.fastq.gz",
        unpaired_r1 = "results/trimmed/{sample}_1.unpaired.fastq.gz",
        trimmed_r2  = "results/trimmed/{sample}_2.trimmed.fastq.gz",
        unpaired_r2 = "results/trimmed/{sample}_2.unpaired.fastq.gz"
    params:
        adapter = config["trim"]["adapter"],
        options = lambda wc: " ".join(
            opt.format(adapter=config["trim"]["adapter"])
            for opt in config["trim"]["options"]
        )
    threads: 1
    log:
        "logs/trimmomatic/{sample}.log"
    conda: "../envs/environment.yaml"
    shell:
        """
        trimmomatic PE \
          -threads {threads} \
          {input.fq1} {input.fq2} \
          {output.trimmed_r1} {output.unpaired_r1} \
          {output.trimmed_r2} {output.unpaired_r2} \
          {params.options}
        """

# FastQC on trimmed reads using Snakemake wrapper
rule fastqc_trimmed_1:
    input:
        fq=lambda wc: f"results/trimmed/{wc.sample}_1.trimmed.fastq.gz"
    output:
        html="results/qc/fastqc/trimmed/{sample}_1.trimmed_fastqc.html",
        zip="results/qc/fastqc/trimmed/{sample}_1.trimmed_fastqc.zip"
    log:
        "logs/fastqc/trimmed/{sample}_1.trimmed.log"
    threads: config["threads"]["qc"]
    wrapper:
        "0.72.0/bio/fastqc"

rule fastqc_trimmed_2:
    input:
        fq=lambda wc: f"results/trimmed/{wc.sample}_2.trimmed.fastq.gz"
    output:
        html="results/qc/fastqc/trimmed/{sample}_2.trimmed_fastqc.html",
        zip="results/qc/fastqc/trimmed/{sample}_2.trimmed_fastqc.zip"
    log:
        "logs/fastqc/trimmed/{sample}_2.trimmed.log"
    threads: config["threads"]["qc"]
    wrapper:
        "0.72.0/bio/fastqc"

# Qualimap2 BAM QC (skip via skip_qualimap flag)
rule qualimap_bamqc:
    input:
        bam="results/bam_sorted/{sample}_sorted.bam"
    output:
        html="results/qc/qualimap/{sample}/qualimapReport.html"
    log:
        "logs/qualimap/{sample}.log"
    threads: config["threads"]["qc"]
    conda: "../envs/environment.yaml"
    shell:
        """
        qualimap bamqc \
          -bam {input.bam} \
          -outdir results/qc/qualimap/{wildcards.sample} \
          &> {log}
        """

# MultiQC aggregation
rule multiqc:
    input:
        # Raw FastQC (if not skipped)
        *([] if config["skip_fastqc"] else expand(
            "results/qc/fastqc/raw/{sample}_1_fastqc.html", sample=SAMPLES)),
        *([] if config["skip_fastqc"] else expand(
            "results/qc/fastqc/raw/{sample}_2_fastqc.html", sample=SAMPLES)),
        # Trimmed FastQC (if not skipped and trimming performed)
        *([] if config["skip_fastqc"] or config["skip_trimming"] else expand(
            "results/qc/fastqc/trimmed/{sample}_1.trimmed_fastqc.html", sample=SAMPLES)),
        *([] if config["skip_fastqc"] or config["skip_trimming"] else expand(
            "results/qc/fastqc/trimmed/{sample}_2.trimmed_fastqc.html", sample=SAMPLES)),
        # Qualimap reports (if not skipped)
        *([] if config["skip_qualimap"] else expand(
            "results/qc/qualimap/{sample}/qualimapReport.html", sample=SAMPLES))
    output:
        html="results/qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    threads: 1
    conda: "../envs/environment.yaml"
    shell:
        "multiqc results/qc -o results/qc &> {log}"
