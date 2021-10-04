import numpy as np
from seq_dict import open_seq
import itertools

def check_seq(seq1, seq2):
    # find shorter and do that range
    change_dict = {}
    nucleotides = ['A', 'T', 'G', 'C']
    combos = list(itertools.product(nucleotides, nucleotides))
    # initialize dictionary
    for i in combos:
        key = str(i[0]+i[1])
        change_dict[key] = 0
    #print(change_dict)
    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] == seq2[i]:
            continue
        elif seq1[i] != seq2[i]:
            # want GT or GC or CG
            nuc_change = str(seq1[i]) + str(seq2[i])
            change_dict[nuc_change] += 1
    return change_dict



#print(check_seq('ATCGATCGATCG', 'ATCGATGGATCC'))


# ancestral = open_seq('./ancestralseqs.fasta.nodes.txt')
# rbcl = open_seq('./grass_rbcl.txt')
#
# ancestral.update(rbcl)
# all_seqs_dict = ancestral
#print(all_seqs_dict.keys())
#print(ancestral.keys())
#print(ancestral['Node1'])

#print(check_seq(all_seqs_dict['Node1'],all_seqs_dict['Node2']))


#
#
# def two_d_arr() -> np.array:
#     arr_zero = np.full((5, 5), 0, dtype=int)
#     return arr_zero
#
# print(two_d_arr())
