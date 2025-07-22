rule featurecounts:
    input:
        bam = rules.star_align.output.bam,
        gtf = config["annotation_gtf"]
    output:
        "counts/{sample}.counts.txt"
    threads: DEFAULT_THREADS
    log:
        "logs/featurecounts/{sample}.log"
    conda: "envs/subread.yaml"
    shell:
        """
        featureCounts -T {threads} \
                      -p -B -C \
                      -a {input.gtf} \
                      -o {output} {input.bam} \
                      > /dev/null 2> {log}
        """
