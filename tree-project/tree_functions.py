import numpy as np
from seq_dict import *

def check_seq(seq1, seq2):
    # find shorter and do that range
    change_dict = {'AA': 0, 'TT': 0, 'CC': 0, 'GG': 0, 'AT': 0, 'TA': 0, 'CG': 0, 'GC': 0,
                   'TC': 0, 'TG': 0, 'CT': 0, 'GT': 0, 'AC': 0, 'AG': 0, 'CA': 0, 'GA': 0}
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            continue
        elif seq1[i] != seq2[i]:
            # want GT or GC or CG
            nuc_change = str(seq1[i]) + str(seq2[i])
            change_dict[nuc_change] += 1
    return change_dict


for k,v in ancestral.items():
    print(k)
    #seqs = list(ancestral.values())

#print(check_seq(seqs[0], seqs[1]))


'''def two_d_arr() -> np.array:
     arr_zero = np.full((5, 5), 0, dtype=int)
     return arr_zero

print(two_d_arr())'''
