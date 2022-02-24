"""
The TreeTracer class takes a given Newick file for a tree representation and a corresponding list of
FASTA sequences for the tree.
This class will parse the information and create a tree object using Phylo package from Biopython
You can apply specific functions/analysis by specifying when creating an instance of TreeTracer class
"""

"""
TODO:
- outgroup specification for loop, don't include the branch comparisons - DONE
- store the dictionary of all of the species to species relationships - ?DONE
- function to only include same codons throughout the tree
- sum of the off diagonals  

"""


from Bio import Phylo
from treetracer_functions import n0_context, n1_context, n2_context, fourfold_n0_context, fourfold_n1_context, \
    fourfold_n2_context, site_changes
import collections, functools, operator
from collections import Counter
from typing import Callable, List
from matrix import Matrix
import pandas as pd


class TreeTracer:

    def __init__(self, newick_file: str, sequences_file: str, outgroups: List):
        """
        Create a new instance of TreeTracer with the Newick file and Phylo reads it in and stores it in self.tree
        :param newick_file: file corresponding to Newick file with tree representation in Newick format
        :param sequences_file: file with sequences of species and nodes for tree **Make sure names are the same
        :param outgroups: which outgroups correspond to the newick tree --> list of strings. May only have one
        """
        self.tree = Phylo.read(newick_file, 'newick')
        self.file_name = sequences_file
        self.seq_dict = self.read_seqs()
        self.outgroups = outgroups
        self.final_matrix_dict = {}
        self.cumulative_mat = {}
        self.total_branch_length = self.cumulative_branch_length()
        self.function_called = Callable
        self.third_codon_sites = self.get_3rd_codon_sites()
        # changes dictionary of each site --> key is node2:Zea and value is [(CTCG,CTGG),nuc_change,codon_dif=T or F]
        self.site_changes_dict = {}
        # dataframe of all sites and where there are changes
        self.change_site_df = pd.DataFrame
        # call the tree_trace function and get all the possible matrix dictionaries and store as instance variables

    def cumulative_branch_length(self):
        total_length = 0
        for clade in self.tree.find_clades():
            if len(clade.clades) > 0:
                for child in clade.clades:
                    if str(child) in self.outgroups:
                        continue
                    total_length += child.branch_length
        return total_length

    def trace_tree_function(self, function_called: Callable, branch_length=True):
        # iterate through all nodes and call function on parent node and child
        self.function_called = function_called
        list_dicts_n0 = []
        for clade in self.tree.find_clades():
            # print('Clade:', clade)
            if len(clade.clades) > 0:
                # print('children:', clade.clades)
                for child in clade.clades:
                    # skip comparison of outgroup --> go to next for loop iteration line above
                    if str(child) in self.outgroups:
                        continue
                    # print('child: ', child.branch_length)
                    # add if statement to check if in dictionary
                    if str(clade) in self.seq_dict and str(child) in self.seq_dict:
                        # print('pair:',clade, ' and ',child)
                        # print(type(clade.name))
                        seq1 = self.seq_dict[clade.name]
                        seq2 = self.seq_dict[child.name]
                        key = str(clade.name + ', ' + child.name)
                        # call the function and return that dictionary
                        if branch_length:
                            output_dict = function_called(seq1, seq2, increment=child.branch_length, codon_sites=self.third_codon_sites)
                        else:
                            output_dict = function_called(seq1, seq2, increment=1.0, codon_sites=self.third_codon_sites)
                        self.final_matrix_dict[key] = output_dict
                        # update the cumulative matrix --> different if its n0 context because no neighboring pairs
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
                if len(list_dicts_n0) > 0 and (function_called == n0_context or function_called == fourfold_n0_context):
                    self.cumulative_mat = dict(functools.reduce(operator.add, map(collections.Counter, list_dicts_n0)))
        return self.cumulative_mat

    def print_cumulative_matrices(self):
        # create an instance of Matrix class
        matrix = Matrix(self.cumulative_mat)
        if self.function_called == n0_context or self.function_called == fourfold_n0_context:
            return matrix.n0_matrix()
        elif self.function_called == n1_context or self.function_called == n2_context or \
                self.function_called == fourfold_n1_context or self.function_called == fourfold_n2_context:
            return matrix.ngt0_matrix()
        return None

    def get_3rd_codon_sites(self):
        """
        #consider adding in degeneracy here!!
        find all third codon sites
        in the functions, check if those sites are 2 or 4 fold
        Assign:
        self.third_codon_sites
        :return: sites if successfully assign self.third_codon_sites list []
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

    def special_trace_tree_function(self):
        # iterate through all nodes and call function on parent node and child
        list_dicts_n0 = []
        for clade in self.tree.find_clades():
            # print('Clade:', clade)
            if len(clade.clades) > 0:
                # print('children:', clade.clades)
                for child in clade.clades:
                    # skip comparison of outgroup --> go to next for loop iteration line above
                    if str(child) in self.outgroups:
                        #print('outgroup:', child)
                        continue
                    # print('child: ', child.branch_length)
                    # add if statement to check if in dictionary
                    if str(clade) in self.seq_dict and str(child) in self.seq_dict:
                        #print('pair:',clade, ' and ',child)
                        # print(type(clade.name))
                        seq1 = self.seq_dict[clade.name]
                        seq2 = self.seq_dict[child.name]
                        key = str(clade.name + ', ' + child.name)
                        # call the function and return that dictionary
                        rd = site_changes(seq1, seq2, child.branch_length, self.third_codon_sites)
                        self.site_changes_dict[str(clade.name+':'+child.name)] = rd
        return self.site_changes_dict

    def get_site_matrix(self):
        matrix = Matrix(all_site_dict=self.site_changes_dict)
        print(matrix.site_matrix())
        print("off diagonal")
        matrix.sum_off_diagonal()

    def get_site_dataframe(self):
        matrix = Matrix(all_site_dict=self.site_changes_dict)
        #matrix.site_changes_table()
        # add
        self.change_site_df = matrix.site_changes_dictionary_to_df(round(self.total_branch_length, 4))
        #print(self.change_site_df)
        cdf = matrix.condensed_change_site_dict(min_tree_prop=0.6)
        #print(cdf)
        #matrix.graph_freq_distribution_matplotlib()
        matrix.graph_freq_distribution_seaborn()
        return self.change_site_df


if __name__ == '__main__':
    # newick_path = '/Users/shailafye/Documents/Morton-Research/2021-research/mRNA-Morton-Project/tree-project/iqtree_newick.txt'
    # seq_path = '/Users/shailafye/Documents/Morton-Research/2021-research/mRNA-Morton-Project/tree-project/grass_rbcl_nodes_seq_fasta.txt'
    # tree_obj = TreeTracer(newick_path, seq_path, outgroups=['Lilium'])

    newick_path = '/Users/shailafye/Documents/Morton-Research/2021-research/all_rbcl_seqs_Newick.txt'
    seq_path = '/Users/shailafye/Documents/Morton-Research/2021-research/all_rbcL_seqs.txt'
    tree_obj = TreeTracer(newick_path, seq_path, outgroups=['Pinus', 'Ginkgo', 'Zamia', 'Welwitschi'])

    tree_obj.special_trace_tree_function()
    tree_obj.get_site_dataframe()
    #tree_obj.draw_tree(tree_type="normal")





