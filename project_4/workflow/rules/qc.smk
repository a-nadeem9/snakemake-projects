rule fastp_trim:
    input:
        fq1 = lambda wc: os.path.join(config["input_dir"], f"{wc.sample}_1.fastq.gz"),
        fq2 = lambda wc: os.path.join(config["input_dir"], f"{wc.sample}_2.fastq.gz")
    output:
        trimmed1 = os.path.join(config["results"]["qc_trimmed"], "{sample}_1.trimmed.fastq.gz"),
        trimmed2 = os.path.join(config["results"]["qc_trimmed"], "{sample}_2.trimmed.fastq.gz"),
        html     = os.path.join(config["logs"]["fastp"], "{sample}_fastp.html"),
        json     = os.path.join(config["logs"]["fastp"], "{sample}_fastp.json")
    params:
        adapter = config["adapter_path"],
        quality = config["parameters"]["fastp"]["qualified_quality_phred"],
        unqual  = config["parameters"]["fastp"]["unqualified_percent_limit"],
        nlimit  = config["parameters"]["fastp"]["n_base_limit"],
        length  = config["parameters"]["fastp"]["length_required"],
        detect  = "--detect_adapter_for_pe" if config["parameters"]["fastp"]["detect_adapter_for_pe"] else ""
    threads: 4
    conda:
        "../envs/fastp.yaml"
    log:
        os.path.join(config["logs"]["fastp"], "{sample}.log")
    shell:
        """
        fastp \
            --in1 {input.fq1} \
            --in2 {input.fq2} \
            --out1 {output.trimmed1} \
            --out2 {output.trimmed2} \
            --adapter_fasta {params.adapter} \
            --qualified_quality_phred {params.quality} \
            --unqualified_percent_limit {params.unqual} \
            --n_base_limit {params.nlimit} \
            --length_required {params.length} \
            {params.detect} \
            --thread {threads} \
            --html {output.html} \
            --json {output.json} \
            > {log} 2>&1
        """
