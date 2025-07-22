rule kallisto:
    input:
        assembly="de_novo/Trinity/{sample}.fasta",
        r1="merged/{sample}_R1.merged.fastq",
        r2="merged/{sample}_R2.merged.fastq"
    output:
        quant="kallisto/{sample}.kallisto.txt"
    threads: 4
    conda: "envs/kallisto.yaml"
    shell:
        r"""
        # make sure all output dirs exist
        mkdir -p kallisto_index kallisto_out/{wildcards.sample} kallisto

        # only build the index once
        if [ ! -f kallisto_index/{wildcards.sample}.idx ]; then
            kallisto index \
                -i kallisto_index/{wildcards.sample}.idx \
                {input.assembly}
        fi

        # run quantification
        kallisto quant \
            -i kallisto_index/{wildcards.sample}.idx \
            -t {threads} \
            -o kallisto_out/{wildcards.sample} \
            {input.r1} {input.r2}

        # copy out the abundance table
        cp kallisto_out/{wildcards.sample}/abundance.tsv {output.quant}
        """
