from Bio import Phylo

"""
parse(), read(), write() and convert()
--Each function accepts either a file name or an open file handle, 
    so data can be also loaded from compressed files, StringIO objects, and so on.
-- The second argument to each function is the target format. Currently, the following formats are supported:
    newick, nexus, nexml, phyloxml, cdao
"""

tree = Phylo.parse()