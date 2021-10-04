"""
This file iterates through the Newick tree parsed data and applies a function on each branch.
A class that creates this tree dictionary and then functions can be called.
"""
import numpy as np
from Bio import Phylo


class TreeTracer:
    """
    newick_file = 'iqtree_newick.txt'

    """

    def __init__(self, newick_file: str, file_name):
        self.tree = Phylo.read(newick_file, 'newick')
        self.file_name = file_name
        self.seq_dict = self.read_seqs()

    def read_seqs(self):
        all_sequences = {}
        infile = open(self.file_name, 'r')
        lines = infile.readlines()
        infile.close()
        for line in range(0, len(lines), 2):
            gene_name = lines[line].strip()[1:]
            gene_name = gene_name.strip(' ').split(" ")[0]
            gene_name = gene_name.split(" ")[0]
            all_sequences[gene_name] = lines[line + 1].strip()
        return all_sequences


    def draw_tree(self, tree_type:str = 'ascii'):
        """
        :param tree_type: ascii or normal
        :return:
        """
        if tree_type == 'ascii':
            Phylo.draw_ascii(self.tree)
        elif tree_type == 'normal':
            self.tree.ladderize()
            Phylo.draw(self.tree)

    def print_tree(self):
        print(self.tree)
        return






tree_obj = TreeTracer('iqtree_newick.txt', 'grass_rbcl_nodes_seq_fasta.txt')
print(tree_obj.seq_dict.keys())
#tree_obj.draw_tree('ascii')
#tree_obj.print_tree()