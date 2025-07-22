# workflow/rules/samtools.smk

# Convert SAM to BAM
rule sam_to_bam:
    input:
        sam="results/sam/{sample}.sam"
    output:
        bam="results/bam/{sample}.bam"
    threads: config["threads"]["samtools"]
    log: "logs/{sample}.sam2bam.log"
    conda: "../../envs/environment.yaml"
    shell:
        """
        samtools view -@ {threads} -bS {input.sam} > {output.bam} 2> {log}
        """

# Sort BAM
rule sort_bam:
    input:
        bam="results/bam/{sample}.bam"
    output:
        sorted="results/bam_sorted/{sample}.sorted.bam"
    threads: config["threads"]["samtools"]
    log: "logs/{sample}.sort.log"
    conda: "../../envs/environment.yaml"
    shell:
        """
        samtools sort --threads {threads} {input.bam} -o {output.sorted} 2> {log}
        """

# Index sorted BAM
rule index_bam:
    input:
        sorted="results/bam_sorted/{sample}.sorted.bam"
    output:
        bai="results/bam_sorted/{sample}.sorted.bam.bai"
    threads: config["threads"]["samtools"]
    log: "logs/{sample}.index.log"
    conda: "../../envs/environment.yaml"
    shell:
        """
        samtools index {input.sorted} 2> {log}
        """

# Compute idxstats
rule idxstats:
    input:
        sorted="results/bam_sorted/{sample}.sorted.bam",
        bai="results/bam_sorted/{sample}.sorted.bam.bai"
    output:
        stats="results/stats/{sample}.stats"
    threads: config["threads"]["samtools"]
    log: "logs/{sample}.idxstats.log"
    conda: "../../envs/environment.yaml"
    shell:
        """
        samtools idxstats {input.sorted} > {output.stats} 2> {log}
        """

# Aggregate all idxstats
rule aggregate_idxstats:
    input:
        stats=expand("results/stats/{sample}.stats", sample=SAMPLES)
    output:
        summary="results/idxstats_summary.tsv"
    log: "logs/aggregate.log"
    conda: "../../envs/environment.yaml"
    run:
        import pandas as pd
        dfs = []
        for s in SAMPLES:
            df = pd.read_csv(f"results/stats/{s}.stats", sep="\t",
                             names=["ref","length",s,"unmapped"])
            dfs.append(df.set_index("ref")[[ "length", s ]])
        merged = pd.concat(dfs, axis=1)
        merged.to_csv(output.summary, sep="\t")
        with open(log[0], "w") as f:
            f.write("Aggregated idxstats into " + output.summary)
