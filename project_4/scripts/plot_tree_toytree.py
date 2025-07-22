import toytree
import toyplot.svg

tree_file = snakemake.input.tree
output_svg = snakemake.output.svg

tree = toytree.tree(tree_file)
draw_output = tree.draw()
canvas = draw_output[0] if isinstance(draw_output, tuple) else draw_output
toyplot.svg.render(canvas, output_svg)
