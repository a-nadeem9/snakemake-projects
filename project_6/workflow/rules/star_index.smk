rule star_index:
    input:
        fasta=config["genome_fasta"],
        gtf=config["annotation_gtf"]
    output:
        directory(config["star_index"])
    threads: DEFAULT_THREADS
    params:
        overhang = int(config["read_length"]) - 1
    log:
        "logs/star/index.log"
    conda:
        "envs/star.yaml"
    shell:
        """
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.overhang} \
             > /dev/null 2> {log}
        """
