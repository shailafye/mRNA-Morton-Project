"""
This file iterates through the Newick tree parsed data and applies a function on each branch.
A class that creates this tree dictionary and then functions can be called.
"""
import itertools
import numpy as np
from Bio import Phylo
from tree_functions import n0_context, n1_context, globdict
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

    def trace_tree_function(self, function_called, branch_length = True):
        # iterate through all nodes
        list_dicts = []
        for clade in self.tree.find_clades():
            #print('Clade:', clade)
            if len(clade.clades) > 0:
                #print('children:', clade.clades)
                for child in clade.clades:
                    #print('child: ', child.branch_length)
                    # add if statement to check if in dictionary
                    if str(clade) in self.seq_dict and str(child) in self.seq_dict:
                        # print('pair:',clade, ' and ',child)
                        # print(type(clade.name))
                        seq1 = self.seq_dict[clade.name]
                        seq2 = self.seq_dict[child.name]
                        key = str(clade.name + ', ' + child.name)
                        self.final_matrix_dict[key] = function_called(seq1, seq2)
                        if branch_length:
                            list_dicts.append(function_called(seq1, seq2, increment=child.branch_length))
                        else:
                            list_dicts.append(function_called(seq1, seq2, increment=1))
                    else:
                        print('not in dictionary with sequences')
            # print('====')
        try:
            self.sum_matrix_dict = dict(functools.reduce(operator.add, map(collections.Counter, list_dicts)))
        except:
            pass
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
            print(key, '\n', matrix_arr)
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


class Matrix:
    def __init__(self, total_matrix_dict, sep_matrix_dict, sum_matrix_dict):
        self.total_matrix_dict = total_matrix_dict
        self.sep_matrix_dict = sep_matrix_dict
        self.comb_matrix_dict = sum_matrix_dict

    def n1_matrix(self):
        for key in self.total_matrix_dict:
            arr_zero = np.full((4, 4), 0, dtype=int)
            matrix_dict = self.total_matrix_dict[key]
            matrix_arr = arr_zero
            conversion = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
            for k in matrix_dict:
                val = matrix_dict[k]
                i = k[0]
                j = k[1]
                matrix_arr[conversion[i]][conversion[j]] = float(val)
            print(key, '\n', matrix_arr)
        return True


tree_obj = TreeTracer('iqtree_newick.txt', 'grass_rbcl_nodes_seq_fasta.txt')
tree_obj.trace_tree_function(n1_context, branch_length=False)
sep_matrix_dict = tree_obj.final_matrix_dict
sum_matrix_dict = tree_obj.sum_matrix_dict
#print(tree_obj.final_matrix_dict)
#print(tree_obj.sum_matrix_dict)
#print('global',globdict)
tree_matrix = Matrix(globdict, sep_matrix_dict, sum_matrix_dict)
tree_matrix.n1_matrix()




