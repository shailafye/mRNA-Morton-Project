"""
This file: treetracer_functions.py is a file with helper functions for TreeTracer class
"""

import itertools

FOURFOLD = ['CTT', 'CTC', 'CTA', 'CTG', 'GTT', 'GTC', 'GTA', 'GTG', 'ACT', 'ACC', 'ACA', 'ACG',
            'TCT', 'TCC', 'TCA', 'TCG', 'CCT', 'CCC', 'CCA', 'CCG', 'GCT', 'GCC', 'GCA', 'GCG',
            'GGT', 'GGC', 'GGA', 'GGG', 'CGT', 'CGC', 'CGA', 'CGG']
TWOFOLD = ['TTT', 'TTC', 'TTA', 'TTG', 'TAT', 'TAC', 'TGT', 'TGC', 'CAT', 'CAC', 'CAA', 'CAG',
           'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'AGT', 'AGC', 'AGA', 'AGG']


def n0_context(seq1, seq2, increment=1.0, codon_sites=[]):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
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


def n1_context(seq1, seq2, increment=1.0, codon_sites=[]):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    increment = float(increment)
    neighboring_nuc_dict = {}
    for i in range(1, min(len(seq1), len(seq2)) - 1):
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


def n2_context(seq1, seq2, increment=1.0, codon_sites=[]):
    increment = float(increment)
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    # loop through sequence, if one is different, make the key the one before and one after
    neighboring_nuc_dict = {}
    # change_dict = {}
    for i in range(2, min(len(seq1), len(seq2)) - 2):
        nuc_change_key = str(seq1[i]) + str(seq2[i])  # AT, TG
        # change_dict[nuc_change] += 1
        # check to see if neighboring nt exists in dictionary
        neighbor_key1 = str(seq1[i - 2:i]) + '_' + str(seq1[i + 1:i + 3])  # key = X_X of sequence 1
        neighbor_key2 = str(seq2[i - 2:i]) + '_' + str(seq2[i + 1:i + 3])  # key = X_X of sequence 2
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


def fourfold_n0_context(seq1, seq2, increment=1.0, codon_sites=[]):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    change_dict = {}
    nucleotides = ['A', 'T', 'G', 'C']
    combos = list(itertools.product(nucleotides, nucleotides))
    # initialize dictionary
    for i in combos:
        key = str(i[0] + i[1])
        change_dict[key] = 0
    for i in range(min(len(seq1), len(seq2))):
        if i not in codon_sites:
            continue
        if not check_4fold(seq1[i - 2:i+1]):
            continue
        nuc_change = str(seq1[i]) + str(seq2[i])
        change_dict[nuc_change] += increment
    return change_dict


def fourfold_n1_context(seq1, seq2, increment=1.0, codon_sites=[]):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    increment = float(increment)
    neighboring_nuc_dict = {}
    for i in range(2, min(len(seq1), len(seq2)) - 1):
        if i not in codon_sites:
            continue
        if not check_4fold(seq1[i - 2:i+1]):
            continue
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


def fourfold_n2_context(seq1, seq2, increment=1.0, codon_sites=[]):
    increment = float(increment)
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    # loop through sequence, if one is different, make the key the one before and one after
    neighboring_nuc_dict = {}
    # change_dict = {}
    for i in range(2, min(len(seq1), len(seq2)) - 2):
        if i not in codon_sites:
            continue
        if not check_4fold(seq1[i - 2:i+1]):
            continue
        nuc_change_key = str(seq1[i]) + str(seq2[i])  # AT, TG
        # change_dict[nuc_change] += 1
        # check to see if neighboring nt exists in dictionary
        neighbor_key1 = str(seq1[i - 2:i]) + '_' + str(seq1[i + 1:i + 3])  # key = X_X of sequence 1
        neighbor_key2 = str(seq2[i - 2:i]) + '_' + str(seq2[i + 1:i + 3])  # key = X_X of sequence 2
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


def check_4fold(codon):
    """
    input is a string codon --> length 3
    checks if site is a 4fold degenerate
    :return: true or false if not a 4 fold deg site
    """
    if len(codon) != 3:
        return False
    return codon in FOURFOLD


def check_2fold(codon):
    """
    input is a string codon --> length 3
    checks if site is a 2fold degenerate
    :return: true or false if not a 2 fold deg site
    """
    if len(codon) != 3:
        return False
    return codon in TWOFOLD




# #
# seq1 = 'ATGACACAGATTTCGAGACGGAGGGCCCATATATAGCGATCTAGCGAGGCGACTCAT'
# seq2 = 'ATGACcCAGATaTCGAGACGGAGGGCCCATATATAGCGATCTAGCGAGGCGACTCAT'
# # AC_CA and AT_TC
#
# list_sites = []
# for i in range(2,len(seq2),3):
#     list_sites.append(i)
#
# dicts = fourfold_n1_context(seq1, seq2, codon_sites = list_sites)
# #print(dicts)
# print(dicts['C_C'])
