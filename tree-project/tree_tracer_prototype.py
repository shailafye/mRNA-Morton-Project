"""
This file iterates through the Newick tree parsed data and applies a function on each branch.
A class that creates this tree dictionary and then functions can be called.
"""
import itertools
import numpy as np
from Bio import Phylo
from tree_functions import n0_context, n1_context
import collections, functools, operator



class TreeTracer:
    """
    newick_file = 'iqtree_newick.txt'
    Create helper functions that take in two sequences to compute along branch
    The tree tracer method takes in a function as a parameters and inputs two sequences at a time
    """

    def __init__(self, newick_file: str, file_name):
        self.tree = Phylo.read(newick_file, 'newick')
        self.file_name = file_name
        self.seq_dict = self.read_seqs()
        self.final_matrix_dict = {}
        self.sum_matrix_dict = {}

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

    def draw_tree(self, tree_type: str = 'ascii'):
        """
        :param tree_type: ascii or normal
        :return: Bool
            True if tree is drawn
        """
        if tree_type == 'ascii':
            Phylo.draw_ascii(self.tree)
            return True
        elif tree_type == 'normal':
            self.tree.ladderize()
            Phylo.draw(self.tree)
            return True
        return False

    def print_tree(self):
        print(self.tree)
        return

    def trace_tree_function(self, function_called):
        # iterate through all nodes
        list_dicts = []
        for clade in self.tree.find_clades():
            #print('Clade:', clade)
            if len(clade.clades) > 0:
                #print('children:', clade.clades)
                for child in clade.clades:
                    #print('child: ', child)
                    # add if statement to check if in dictionary
                    if str(clade) in self.seq_dict and str(child) in self.seq_dict:
                        # print('pair:',clade, ' and ',child)
                        # print(type(clade.name))
                        seq1 = self.seq_dict[clade.name]
                        seq2 = self.seq_dict[child.name]
                        key = str(clade.name + ', ' + child.name)
                        self.final_matrix_dict[key] = function_called(seq1, seq2)
                        list_dicts.append(function_called(seq1, seq2))
                    else:
                        print('not in dictionary with sequences')
            #print('====')
        self.sum_matrix_dict = dict(functools.reduce(operator.add, map(collections.Counter, list_dicts)))
        return True

    def individual_matrix_conversion(self) -> np.array:
        for key in self.final_matrix_dict:
            arr_zero = np.full((4, 4), 0, dtype=int)
            matrix_dict = self.final_matrix_dict[key]  # 'Node1, Lilium'
            matrix_arr = arr_zero
            conversion = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
            for k in matrix_dict:
                val = matrix_dict[k]
                i = k[0]
                j = k[1]
                matrix_arr[conversion[i]][conversion[j]] = int(val)
            print(key, '\n',matrix_arr)
        return True

    def cumulative_matrix_conversion(self):
        arr_zero = np.full((4, 4), 0, dtype=int)
        matrix_dict = self.sum_matrix_dict
        matrix_arr = arr_zero
        conversion = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
        for k in matrix_dict:
            val = matrix_dict[k]
            i = k[0]
            j = k[1]
            matrix_arr[conversion[i]][conversion[j]] = int(val)
        print(matrix_arr)
        return True

tree_obj = TreeTracer('iqtree_newick.txt', 'grass_rbcl_nodes_seq_fasta.txt')
tree_obj.trace_tree_function(n1_context)
print(tree_obj.final_matrix_dict)
tree_obj.cumulative_matrix_conversion()
#tree_obj.matrix_conversion()
# print(tree_obj.seq_dict.keys())
# tree_obj.draw_tree('normal')
# tree_obj.draw_tree('ascii')
# tree_obj.print_tree()
