# workflow/rules/samtools.smk

rule convert_bam:
    input:
        sam="results/sam/{sample}.sam"
    output:
        bam="results/bam/{sample}.bam"
    threads: config["threads"]["map_reads"]
    log:
        "logs/samtools/convert_{sample}.log"
    conda: "../envs/environment.yaml"
    shell:
        "samtools view -@ {threads} -bS {input.sam} > {output.bam} 2> {log} && touch {output.bam}"

rule sort_bam:
    input:
        bam="results/bam/{sample}.bam"
    output:
        sorted="results/bam_sorted/{sample}_sorted.bam"
    threads: config["threads"]["map_reads"]
    log:
        "logs/samtools/sort_{sample}.log"
    conda: "../envs/environment.yaml"
    shell:
        "samtools sort --threads {threads} -o {output.sorted} {input.bam} &> {log}"

rule index_bam:
    input:
        sorted="results/bam_sorted/{sample}_sorted.bam"
    output:
        bai="results/bam_sorted/{sample}_sorted.bam.bai"
    log:
        "logs/samtools/index_{sample}.log"
    conda: "../envs/environment.yaml"
    shell:
        "samtools index {input.sorted} &> {log}"

rule idxstats:
    input:
        sorted="results/bam_sorted/{sample}_sorted.bam"
    output:
        stats="results/stats/{sample}.stats.txt"
    log:
        "logs/samtools/idxstats_{sample}.log"
    conda: "../envs/environment.yaml"
    shell:
        "samtools idxstats {input.sorted} > {output.stats} 2> {log}"

rule aggregate_idxstats:
    input:
        expand("results/stats/{sample}.stats.txt", sample=SAMPLES)
    output:
        "results/idxstats_summary.tsv"
    conda: "../envs/environment.yaml"
    run:
        import pandas as pd
        df_list = []
        for sample_name, fn in zip(SAMPLES, input):
            df = pd.read_csv(fn, sep="\t", header=None,
                             names=["ref","length","mapped","unmapped"])
            df = df[["ref","mapped"]].set_index("ref")
            df.columns = [sample_name]
            df_list.append(df)
        summary = pd.concat(df_list, axis=1)
        summary.to_csv(output[0], sep="\t")
