"""
The TreeTracer class takes a given Newick file for a tree representation and a corresponding list of
FASTA sequences for the tree.
This class will parse the information and create a tree object using Phylo package from Biopython
You can apply specific functions/analysis by specifying when creating an instance of TreeTracer class
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
        # dataframe of sites with same context and corresponding changes/counts
        self.condensed_final_site_df = pd.DataFrame
        # call the tree_trace function and get all the possible matrix dictionaries and store as instance variables

    def cumulative_branch_length(self):
        """
        sum the total length of branches without the outgroup branches
        """
        total_length = 0
        for clade in self.tree.find_clades():
            if len(clade.clades) > 0:
                for child in clade.clades:
                    if str(child) in self.outgroups:
                        continue
                    total_length += child.branch_length
        return total_length

    def get_3rd_codon_sites(self):
        """
        Find all third codon sites indices/positions of the sequence
        Assume that the sequences start inframe
        Assign: self.third_codon_sites
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

    def trace_tree_function(self, function_called: Callable, branch_length=True):
        """

        :param function_called:
        :param branch_length:
        :return:
        """
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
                        # store the output in a dictionary with the key being the two sequences
                        # and the output of the function called
                        self.final_matrix_dict[key] = output_dict
                        # update the cumulative matrix --> different if its n0 context because no neighboring pairs
                        # add the output all the branches and store in self.cumulative_mat
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
        return True

    def site_trace_tree_function(self):
        # iterate through all nodes and call function on parent node and child
        for clade in self.tree.find_clades():
            # print('Clade:', clade)
            if len(clade.clades) > 0:
                # print('children:', clade.clades)
                for child in clade.clades:
                    # skip comparison of outgroup --> go to next for loop iteration line above
                    if str(child) in self.outgroups:
                        continue
                    # print('child: ', child.branch_length)
                    # add if statement to check if in dictionary of sequences
                    if str(clade) in self.seq_dict and str(child) in self.seq_dict:
                        # print('pair:',clade, ' and ',child)
                        # print(type(clade.name))
                        seq1 = self.seq_dict[clade.name]
                        seq2 = self.seq_dict[child.name]
                        key = str(clade.name + ', ' + child.name)
                        # call the function and return that dictionary
                        rd = site_changes(seq1, seq2, child.branch_length, self.third_codon_sites)
                        self.site_changes_dict[str(clade.name+':'+child.name)] = rd
        return self.site_changes_dict

    def site_change_analysis(self, site_matrices=False, to_csv=False, show_graphs=False, save_graphs=False, run_stats=False):
        """
        Function to analyze the site by site changes with the same context
        Create an instance of matrix with self.site_changes_dict
        site_matrices: if you want to print out matrices for each site, the
        :return:
        """
        # create a matrix object with the self.site_changes_dict
        matrix = Matrix(all_site_dict=self.site_changes_dict)

        # print out site matrices
        if site_matrices:
            matrix.site_matrix()

        # generate a dataframe of the different sites and context
        # this is an expanded dataframe view and includes ALL branches and contexts
        self.change_site_df = matrix.site_changes_dictionary_to_df(round(self.total_branch_length, 4))

        # need to reduce the dataframe to only sites of interest
        # sites that have the same context and are at minimum a certain percentage of the tree
        self.condensed_final_site_df = matrix.condensed_change_site_dict(min_tree_prop=0.6, to_csv=to_csv)

        # show frequency distribution graphs of sites with different GC context
        matrix.graph_freq_distribution_seaborn(show_graphs=show_graphs, save_graphs=save_graphs, run_stats=run_stats)

        return True


if __name__ == '__main__':
    # newick_path = '/Users/shailafye/Documents/Morton-Research/2021-research/mRNA-Morton-Project/tree-project/iqtree_newick.txt'
    # seq_path = '/Users/shailafye/Documents/Morton-Research/2021-research/mRNA-Morton-Project/tree-project/grass_rbcl_nodes_seq_fasta.txt'
    # tree_obj = TreeTracer(newick_path, seq_path, outgroups=['Lilium'])

    newick_path = '/Users/shailafye/Documents/Morton-Research/2021-research/all_rbcl_seqs_Newick.txt'
    seq_path = '/Users/shailafye/Documents/Morton-Research/2021-research/all_rbcL_seqs.txt'
    tree_obj = TreeTracer(newick_path, seq_path, outgroups=['Pinus', 'Ginkgo', 'Zamia', 'Welwitschi'])

    tree_obj.trace_tree_function(n2_context, branch_length=False)
    print(tree_obj.cumulative_mat)
    tree_obj.print_cumulative_matrices()


    # tree_obj.site_trace_tree_function()
    # tree_obj.site_change_analysis(to_csv=False, show_graphs=False, save_graphs=False, run_stats=False)
    #
    # print(tree_obj.condensed_final_site_df)
    #tree_obj.draw_tree(tree_type="normal")
