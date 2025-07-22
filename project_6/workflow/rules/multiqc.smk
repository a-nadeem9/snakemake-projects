rule multiqc:
    input:
        expand("qc/fastqc/{sample}_R1.trimmed_fastqc.zip", sample=SAMPLES),
        expand("qc/fastp/{sample}.fastp.json", sample=SAMPLES)
    output:
        "qc/multiqc/multiqc_report.html"
    conda: "envs/multiqc.yaml"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        """
        mkdir -p qc/multiqc logs/multiqc
        
        # Run MultiQC without redirect to capture errors in console
        multiqc qc/fastqc qc/fastp \
            --name multiqc_report.html \
            --outdir qc/multiqc 2>&1 | tee {log}

        # Verify report
        test -f {output} || (echo "Error: {output} not found" >&2; exit 1)
        """
