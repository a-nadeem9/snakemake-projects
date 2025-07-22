rule trinity:
    input:
        merged_r1="merged/{sample}_R1.merged.fastq",
        merged_r2="merged/{sample}_R2.merged.fastq"
    output:
        assembly="de_novo/Trinity/{sample}.fasta",
        map="de_novo/Trinity/{sample}.map"
    threads: 8
    conda: "envs/trinity.yaml"  # now using the environment file with trinity=2.14.0
    shell:
        """
        Trinity --seqType fq --left {input.merged_r1} --right {input.merged_r2} \
                --max_memory 10G --CPU {threads} --output trinity_out/Trinity_{wildcards.sample}
        cp trinity_out/Trinity_{wildcards.sample}/Trinity.fasta {output.assembly}
        get_Trinity_gene_to_trans_map.pl {output.assembly} > {output.map}
        if [ ! -s {output.map} ]; then
            echo "Error: Mapping file is empty. Ensure get_Trinity_gene_to_trans_map.pl is installed and the assembly headers are correct."
            exit 1;
        fi
        """
