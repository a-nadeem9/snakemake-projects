# workflow/rules/bowtie.smk
import os

# Build the Bowtie2 index and touch a marker file so Snakemake sees it
rule build_bowtie_index:
    input:
        ref = config["ref"]
    output:
        marker = "index/genome"
    threads: config["threads"]["bowtie"]
    log: "logs/build_index.log"
    conda: "../../envs/environment.yaml"
    shell:
        """
        bowtie2-build {input.ref} index/genome > {log} 2>&1
        touch {output.marker}
        """

# Map paired reads to the reference
rule map_reads:
    input:
        index = "index/genome",
        r1    = lambda wc: samples_df.at[wc.sample, "fq1"],
        r2    = lambda wc: samples_df.at[wc.sample, "fq2"]
    output:
        sam = "results/sam/{sample}.sam"
    params:
        insert_params = lambda wc: f"-I {config['bowtie']['minins']} -X {config['bowtie']['maxins']}"
    threads: config["threads"]["bowtie"]
    log: "logs/{sample}.bowtie2.log"
    conda: "../../envs/environment.yaml"
    shell:
        """
        bowtie2 \
          --threads {threads} \
          -x {input.index} \
          -1 {input.r1} \
          -2 {input.r2} \
          {params.insert_params} \
          -S {output.sam} \
          2> {log}
        """
