rule rsem:
    input:
        assembly="de_novo/Trinity/{sample}.fasta",
        r1="merged/{sample}_R1.merged.fastq",
        r2="merged/{sample}_R2.merged.fastq"
    output:
        rsem="rsem/{sample}.rsem.txt"
    threads: 4
    conda: "../envs/rsem.yaml"
    shell:
        """
        mkdir -p rsem_ref rsem_out rsem

        # 1) transcript→gene map
        get_Trinity_gene_to_trans_map.pl {input.assembly} \
            > de_novo/Trinity/{wildcards.sample}.map

        # 2) build Bowtie2 index & prepare RSEM reference
        rsem-prepare-reference \
            --transcript-to-gene-map de_novo/Trinity/{wildcards.sample}.map \
            {input.assembly} \
            rsem_ref/{wildcards.sample} \
            --bowtie2

        # 3) quantify expression
        rsem-calculate-expression \
            --bowtie2 \
            --paired-end \
            --num-threads {threads} \
            {input.r1} {input.r2} \
            rsem_ref/{wildcards.sample} \
            rsem_out/{wildcards.sample}

        # 4) copy gene‐level results
        cp rsem_out/{wildcards.sample}.genes.results {output.rsem}
        """
