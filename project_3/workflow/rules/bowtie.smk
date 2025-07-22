# workflow/rules/bowtie.smk

INDEX_PREFIX = "index/genome"
INDEX_SUFFIXES = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]

rule build_bowtie_index:
    input:
        ref = config["reference"]
    output:
        expand("{prefix}.{suffix}", prefix=INDEX_PREFIX, suffix=INDEX_SUFFIXES)
    log:
        "logs/build_index.log"
    threads: config["threads"]["build_index"]
    conda:
        "../envs/environment.yaml"
    shell:
        """
        bowtie2-build --threads {threads} {input.ref} {INDEX_PREFIX} > {log} 2>&1
        """

rule map_reads:
    input:
        index = expand("{prefix}.{suffix}", prefix=INDEX_PREFIX, suffix=INDEX_SUFFIXES),
        r1 = lambda wc: samples.at[wc.sample, "fq1"] if config["skip_trimming"]
                        else f"results/trimmed/{wc.sample}_1.trimmed.fastq.gz",
        r2 = lambda wc: samples.at[wc.sample, "fq2"] if config["skip_trimming"]
                        else f"results/trimmed/{wc.sample}_2.trimmed.fastq.gz"
    output:
        sam = "results/sam/{sample}.sam"
    params:
        insert_params = lambda wc: f"-I {config['bowtie']['minins']} -X {config['bowtie']['maxins']}"
    log:
        "logs/{sample}.bowtie2.log"
    threads: config["threads"]["map_reads"]
    conda:
        "../envs/environment.yaml"
    shell:
        """
        bowtie2 \
          --threads {threads} \
          -x {INDEX_PREFIX} \
          -1 {input.r1} \
          -2 {input.r2} \
          {params.insert_params} \
          -S {output.sam} 2> {log}
        """

