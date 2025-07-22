rule merge_samples:
    input:
        r1=fq1,
        r2=fq2
    output:
        merged_r1="merged/{sample}_R1.merged.fastq",
        merged_r2="merged/{sample}_R2.merged.fastq"
    shell:
        "zcat {input.r1} > {output.merged_r1} && zcat {input.r2} > {output.merged_r2}"
