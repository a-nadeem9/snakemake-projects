rule compute_variability:
    input:
        msa = config["results"]["msa_dir"] + "/aligned_scaffolds.fasta"
    output:
        txt = config["results"]["variability_dir"] + "/variability.txt",
        plot = config["results"]["variability_dir"] + "/variability_plot.png",
        windowed = config["results"]["variability_dir"] + "/variability_windowed.txt"
    params:
        script = "scripts/compute_variability.py"
    conda:
        "../envs/variability.yaml"
    shell:
        """
        mkdir -p {config[results][variability_dir]}
        python {params.script}
        """
