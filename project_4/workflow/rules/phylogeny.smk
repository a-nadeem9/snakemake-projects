# workflow/rules/phylogeny.smk

rule run_iqtree:
    input:
        aln = "results/msa/aligned_scaffolds.fasta"
    output:
        tree = "results/msa/tree.nwk"
    log:
        "logs/iqtree/iqtree.log"
    conda:
        "../envs/iqtree.yaml"
    shell:
        """
        mkdir -p results/msa
        iqtree2 -s {input.aln} -nt AUTO -m GTR+G -bb 1000 -pre results/msa/tree -redo > {log} 2>&1
        mv results/msa/tree.treefile {output.tree}
        """

rule plot_tree_toytree:
    input:
        tree="results/msa/tree.nwk"
    output:
        svg="results/plots/tree.svg"
    conda:
        "../envs/toytree.yaml"
    script:
        "../../scripts/plot_tree_toytree.py"
