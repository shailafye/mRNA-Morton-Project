"""
The TreeTracer class takes a given Newick file for a tree representation and a corresponding list of
FASTA sequences for the tree.
This class will parse the information and create a tree object using Phylo package from Biopython
You can apply specific functions/analysis by specifying when creating an instance of TreeTracer class
"""

from Bio import Phylo
from treetracer_functions import n0_context, n1_context, n2_context, fourfold_n0_context, fourfold_n1_context, fourfold_n2_context
import collections, functools, operator
from collections import Counter
from typing import Callable
from matrix import Matrix


class TreeTracer:

    def __init__(self, newick_file: str, sequences_file: str):
        """
        Create a new instance of TreeTracer with the Newick file and Phylo reads it in and stores it in self.tree
        :param newick_file: file corresponding to Newick file with tree representation in Newick format
        :param sequences_file: file with sequences of species and nodes for tree **Make sure names are the same
        """
        self.tree = Phylo.read(newick_file, 'newick')
        self.file_name = sequences_file
        self.seq_dict = self.read_seqs()
        self.final_matrix_dict = {}
        self.cumulative_mat = {}
        self.function_called = Callable
        self.third_codon_sites = self.get_3rd_codon_sites()

    # def trace_tree_function(self, function_called: Callable, branch_length=True):
    #     # iterate through all nodes and call function on parent node and child
    #     self.function_called = function_called
    #     list_dicts_n0 = []
    #     for clade in self.tree.find_clades():
    #         # print('Clade:', clade)
    #         if len(clade.clades) > 0:
    #             # print('children:', clade.clades)
    #             for child in clade.clades:
    #                 # print('child: ', child.branch_length)
    #                 # add if statement to check if in dictionary
    #                 if str(clade) in self.seq_dict and str(child) in self.seq_dict:
    #                     #print('pair:',clade, ' and ',child)
    #                     # print(type(clade.name))
    #                     seq1 = self.seq_dict[clade.name]
    #                     seq2 = self.seq_dict[child.name]
    #                     key = str(clade.name + ', ' + child.name)
    #                     if branch_length:
    #                         output_dict = function_called(seq1, seq2, increment=child.branch_length)
    #                     else:
    #                         output_dict = function_called(seq1, seq2, increment=1.0)
    #                     self.final_matrix_dict[key] = output_dict
    #                     if function_called != n0_context:
    #                         # update cumulative dictionary of total changes at neighboring sites
    #                         for neighbor in output_dict:
    #                             if neighbor in self.cumulative_mat:
    #                                 one = Counter(self.cumulative_mat[neighbor])
    #                                 two = Counter(output_dict[neighbor])
    #                                 self.cumulative_mat[neighbor] = dict(one + two)
    #                             else:
    #                                 self.cumulative_mat[neighbor] = output_dict[neighbor]
    #                     elif function_called == n0_context:
    #                         list_dicts_n0.append(output_dict)
    #             if len(list_dicts_n0) > 0 and function_called == n0_context:
    #                 self.cumulative_mat = dict(functools.reduce(operator.add, map(collections.Counter, list_dicts_n0)))
    #     return True

    def trace_tree_function(self, function_called: Callable, branch_length=True):
        # iterate through all nodes and call function on parent node and child
        self.function_called = function_called
        list_dicts_n0 = []
        for clade in self.tree.find_clades():
            # print('Clade:', clade)
            if len(clade.clades) > 0:
                # print('children:', clade.clades)
                for child in clade.clades:
                    # print('child: ', child.branch_length)
                    # add if statement to check if in dictionary
                    if str(clade) in self.seq_dict and str(child) in self.seq_dict:
                        #print('pair:',clade, ' and ',child)
                        # print(type(clade.name))
                        seq1 = self.seq_dict[clade.name]
                        seq2 = self.seq_dict[child.name]
                        key = str(clade.name + ', ' + child.name)
                        if branch_length:
                            output_dict = function_called(seq1, seq2, increment=child.branch_length, codon_sites=self.third_codon_sites)
                        else:
                            output_dict = function_called(seq1, seq2, increment=1.0, codon_sites=self.third_codon_sites)
                        self.final_matrix_dict[key] = output_dict
                        if function_called == n1_context or function_called == n2_context or \
                                function_called == fourfold_n1_context or function_called == fourfold_n2_context:
                            # update cumulative dictionary of total changes at neighboring sites
                            for neighbor in output_dict:
                                if neighbor in self.cumulative_mat:
                                    one = Counter(self.cumulative_mat[neighbor])
                                    two = Counter(output_dict[neighbor])
                                    self.cumulative_mat[neighbor] = dict(one + two)
                                else:
                                    self.cumulative_mat[neighbor] = output_dict[neighbor]
                        elif function_called == n0_context or function_called == fourfold_n0_context:
                            list_dicts_n0.append(output_dict)
                if len(list_dicts_n0) > 0 and function_called == n0_context or function_called == fourfold_n0_context:
                    self.cumulative_mat = dict(functools.reduce(operator.add, map(collections.Counter, list_dicts_n0)))
        return True

    def print_cumulative_matrices(self):
        matrix = Matrix(self.cumulative_mat)
        if self.function_called == n0_context or self.function_called == fourfold_n0_context:
            matrix.n0_matrix()
        else:
            matrix.ngt0_matrix()
        return True

    def get_3rd_codon_sites(self):
        """
        find all third codon sites
        in the functions, check if those sites are 2 or 4 fold
        Assign:
        self.third_codon_sites
        :return: True if successfully assign self.third_codon_sites list []
        """
        sites = []
        min_val = max([len(self.seq_dict[ele]) for ele in self.seq_dict])
        for nt_idx in range(2, min_val, 3):
            sites.append(nt_idx)
        return sites

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

    def read_seqs(self) -> dict:
        """
        Method to create a dictionary of the sequences
        Key is the name
        Value is the sequence
        :return: dictionary of sequences for nodes and species
        """
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


tree_obj = TreeTracer('../iqtree_newick.txt', '../grass_rbcl_nodes_seq_fasta.txt')
print('fourfold')
tree_obj.trace_tree_function(fourfold_n0_context, branch_length=False)
tree_obj.print_cumulative_matrices()
print('normal')
tree_obj.trace_tree_function(n0_context, branch_length=False)
tree_obj.print_cumulative_matrices()
#print(tree_obj.final_matrix_dict)
#print("\ncumulative:")
#print(tree_obj.cumulative_mat)
#tree_obj.print_cumulative_matrices()
