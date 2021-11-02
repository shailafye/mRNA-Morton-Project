"""
This file: treetracer_functions.py is a file with helper functions for TreeTracer class
"""

import itertools


def n0_context(seq1, seq2, increment=1.0):
    change_dict = {}
    nucleotides = ['A', 'T', 'G', 'C']
    combos = list(itertools.product(nucleotides, nucleotides))
    # initialize dictionary
    for i in combos:
        key = str(i[0] + i[1])
        change_dict[key] = 0
    for i in range(min(len(seq1), len(seq2))):
        nuc_change = str(seq1[i]) + str(seq2[i])
        change_dict[nuc_change] += increment
    return change_dict


def n1_context(seq1, seq2, increment=1.0):
    increment = float(increment)
    neighboring_nuc_dict = {}
    for i in range(1, min(len(seq1), len(seq2))-1):
        nuc_change_key = str(seq1[i]) + str(seq2[i])  # AT, TG
        # change_dict[nuc_change] += 1
        # check to see if neighboring nt exists in dictionary
        neighbor_key1 = str(seq1[i - 1]) + '_' + str(seq1[i + 1])  # key = X_X of sequence 1
        neighbor_key2 = str(seq2[i - 1]) + '_' + str(seq2[i + 1])  # key = X_X of sequence 2
        if neighbor_key1 != neighbor_key2:
            continue
        if neighbor_key1 not in neighboring_nuc_dict:
            # initialize a new dictionary inside of neighboring_nuc_dict
            change_dict = {}
            nucleotides = ['A', 'T', 'G', 'C']
            combos = list(itertools.product(nucleotides, nucleotides))
            # initialize dictionary
            for k in combos:
                key = str(k[0] + k[1])
                change_dict[key] = 0
            neighboring_nuc_dict[neighbor_key1] = change_dict
        neighboring_nuc_dict[neighbor_key1][nuc_change_key] += increment
    return neighboring_nuc_dict


def n2_context(seq1, seq2, increment=1.0):
    increment = float(increment)
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    # loop through sequence, if one is different, make the key the one before and one after
    neighboring_nuc_dict = {}
    # change_dict = {}
    for i in range(2, min(len(seq1), len(seq2))-2):
        nuc_change_key = str(seq1[i]) + str(seq2[i])  # AT, TG
        # change_dict[nuc_change] += 1
        # check to see if neighboring nt exists in dictionary
        neighbor_key1 = str(seq1[i - 2:i]) + '_' + str(seq1[i+1:i + 3])  # key = X_X of sequence 1
        neighbor_key2 = str(seq2[i - 2:i]) + '_' + str(seq2[i+1:i + 3])  # key = X_X of sequence 2
        if neighbor_key1 != neighbor_key2:
            continue
        if neighbor_key1 not in neighboring_nuc_dict:
            # initialize a new dictionary inside of neighboring_nuc_dict
            change_dict = {}
            nucleotides = ['A', 'T', 'G', 'C']
            combos = list(itertools.product(nucleotides, nucleotides))
            # initialize dictionary
            for k in combos:
                key = str(k[0] + k[1])
                change_dict[key] = 0
            neighboring_nuc_dict[neighbor_key1] = change_dict
        neighboring_nuc_dict[neighbor_key1][nuc_change_key] += increment
    return neighboring_nuc_dict

#
# seq1 = 'ATGACACAGATTTCGAGACGGAGGGCCCATATATAGCGATCTAGCGAC'
# seq2 = 'ATGcCACAGATTaCGAGACGGAGGGCCCATATATAGCGATCTAGCGAC'
# # TG_CA and TT_CG
# dicts = n2_context(seq1, seq2)
# print(dicts['TT_CG'])