############################
# Fastp trimming & QC
############################
rule fastp:
    input:
        r1=fq1,
        r2=fq2
    output:
        r1="trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="trimmed/{sample}_R2.trimmed.fastq.gz",
        html="qc/fastp/{sample}.fastp.html",
        json="qc/fastp/{sample}.fastp.json"
    threads: DEFAULT_THREADS
    log:      "logs/fastp/{sample}.log"
    conda:    "envs/fastp.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              --thread {threads} \
              --html {output.html} --json {output.json} \
              > /dev/null 2> {log}
        """

rule fastqc:
    input:
        r1="trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        "qc/fastqc/{sample}_R1.trimmed_fastqc.html",
        "qc/fastqc/{sample}_R2.trimmed_fastqc.html",
        "qc/fastqc/{sample}_R1.trimmed_fastqc.zip",
        "qc/fastqc/{sample}_R2.trimmed_fastqc.zip"
    threads: 2
    log:   "logs/fastqc/{sample}.log"
    conda: "envs/fastqc.yaml"
    shell:
        """
        fastqc -t {threads} -o qc/fastqc {input.r1} {input.r2} > /dev/null 2> {log}
        """
