# workflow/rules/phylogeny.smk
# ------------------------------------------------------------
# IQ-TREE + Toytree plotting, using the fixed folder layout
# ------------------------------------------------------------
# Globals:  MSA_DIR, PLOTS_DIR

rule run_iqtree:
    threads: config["parameters"]["iqtree"]["threads"]
    input:
        aln = f"{MSA_DIR}/aligned_scaffolds.fasta"
    output:
        tree = f"{MSA_DIR}/tree.nwk"
    params:
        out_prefix = f"{MSA_DIR}/tree"
    log:
        "logs/iqtree/iqtree.log"
    conda: "../envs/iqtree.yaml"
    shell:
        """
        iqtree2 -s {input.aln} -nt {threads} -m GTR+G -bb 1000 \
                -pre {params.out_prefix} -redo 2> {log}
        mv {params.out_prefix}.treefile {output.tree}
        """

rule plot_tree_toytree:
    input:
        tree = rules.run_iqtree.output.tree
    output:
        svg = f"{PLOTS_DIR}/tree.svg",
        pdf = f"{PLOTS_DIR}/tree.pdf",
        png = f"{PLOTS_DIR}/tree.png"
    conda: "../envs/toytree.yaml"
    script: "../../scripts/plot_tree_toytree.py"
