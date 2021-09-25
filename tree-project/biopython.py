from Bio import Phylo
from Bio.Phylo import PhyloXML

import matplotlib

"""
parse(), read(), write() and convert()
--Each function accepts either a file name or an open file handle, 
    so data can be also loaded from compressed files, StringIO objects, and so on.
-- The second argument to each function is the target format. Currently, the following formats are supported:
    newick, nexus, nexml, phyloxml, cdao
"""

tree = Phylo.read('iqtree_newick.txt', 'newick')
# print(tree.root)
# print(tree)
# tree.to_alignment()
# print(tree.distance('Node3', 'Zea'))
# print(tree.clade[1][1])
# print(tree.get_terminals())
# print(tree.get_nonterminals())
# print(tree.trace('Oryza','Triticum'))
for clade in tree.find_clades():
    print('Clade:', clade)
    if len(clade.clades) > 0:
        print(clade.clades)
        print(clade.clades[0])
        print('====')

# draw tree
# Phylo.draw_ascii(tree)
# tree.ladderize()
# Phylo.draw(tree)










