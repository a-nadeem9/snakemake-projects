rule star_align:
    input:
        index = rules.star_index.output,
        r1    = "trimmed/{sample}_R1.trimmed.fastq.gz",
        r2    = "trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        bam = temp("align/{sample}.sorted.bam")
    threads: DEFAULT_THREADS
    log:
        "logs/star/{sample}.log"
    conda:
        "envs/star.yaml"
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix align/{wildcards.sample}. \
             > /dev/null 2> {log}

        mv align/{wildcards.sample}.Aligned.sortedByCoord.out.bam {output.bam}
        """
