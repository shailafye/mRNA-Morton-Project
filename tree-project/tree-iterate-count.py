"""
This is the final file to iterate through the tree and apply functions to get matrix for each node

"""
import numpy as np
from Bio import Phylo
from seq_dict import open_seq
from tree_functions import n0_context

#get sequence dictionary
ancestral = open_seq('./ancestralseqs.fasta.nodes.txt')
rbcl = open_seq('./grass_rbcl.txt')
ancestral.update(rbcl)
all_seqs_dict = ancestral
#print(all_seqs_dict)

# Read in the tree
tree = Phylo.read('iqtree_newick.txt', 'newick')

final_matrix = {} # list of dictionaries

# draw tree
# Phylo.draw_ascii(tree)
# tree.ladderize()
# Phylo.draw(tree)
list_dicts = []
#iterate through all nodes
for clade in tree.find_clades():
    print('Clade:', clade)
    if len(clade.clades) > 0:
        print('children:', clade.clades)
        for child in clade.clades:
            print('child: ', child)
            #add if statement to check if in dictionary
            if str(clade) in all_seqs_dict and str(child) in all_seqs_dict:
                #print('pair:',clade, ' and ',child)
                #print(type(clade.name))
                seq1 = all_seqs_dict[clade.name]
                seq2 = all_seqs_dict[child.name]
                key = str(clade.name+', '+ child.name)
                final_matrix[key] = n0_context(seq1,seq2)
                list_dicts.append(n0_context(seq1,seq2))
            else:
                print('not in dictionary with sequences')
    print('====')


print(final_matrix)

def matrix_conversion(key, final_matrix) -> np.array:
    arr_zero = np.full((4, 4), 0, dtype=int)
    if key:
        matrix_dict = final_matrix[key] #'Node1, Lilium'
    else:
        matrix_dict = final_matrix
    matrix_arr = arr_zero
    conversion = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    for k in matrix_dict:
        val = matrix_dict[k]
        i = k[0]
        j = k[1]
        matrix_arr[conversion[i]][conversion[j]] = int(val)
    return matrix_arr

import collections, functools, operator
print(final_matrix)
# sum the values with same keys
result = dict(functools.reduce(operator.add, map(collections.Counter, list_dicts)))
print("resultant dictionary : ", str(result))
print(matrix_conversion(None, result))

final_matrix_all_pairs = {}
for i in final_matrix:
    #print(i)
    #print(matrix_conversion(i, final_matrix))
    final_matrix_all_pairs[i] = matrix_conversion(i, final_matrix)

#print(final_matrix_all_pairs)
