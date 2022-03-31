"""
This file: treetracer_functions.py is a file with helper functions for TreeTracer class
"""
import pandas as pd
import itertools

FOURFOLD = ['CTT', 'CTC', 'CTA', 'CTG', 'GTT', 'GTC', 'GTA', 'GTG', 'ACT', 'ACC', 'ACA', 'ACG',
            'TCT', 'TCC', 'TCA', 'TCG', 'CCT', 'CCC', 'CCA', 'CCG', 'GCT', 'GCC', 'GCA', 'GCG',
            'GGT', 'GGC', 'GGA', 'GGG', 'CGT', 'CGC', 'CGA', 'CGG']
TWOFOLD = ['TTT', 'TTC', 'TTA', 'TTG', 'TAT', 'TAC', 'TGT', 'TGC', 'CAT', 'CAC', 'CAA', 'CAG',
           'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'AGT', 'AGC', 'AGA', 'AGG']


def check_nt(nt1, nt2):
    """
    Function that takes two nucleotides and confirms they are ATGC and not - or X
    """
    nucleotides = ['A', 'T', 'G', 'C']
    if nt1.upper() in nucleotides and nt2.upper() in nucleotides:
        return True
    else:
        return False


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
        if not check_nt(str(seq1[i]), str(seq2[i])):
            continue
        nuc_change = str(seq1[i]) + str(seq2[i])
        change_dict[nuc_change] += increment
    return change_dict


def n1_context(seq1, seq2, increment=1.0, codon_sites=[]):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    increment = float(increment)
    neighboring_nuc_dict = {}
    for i in range(1, min(len(seq1), len(seq2)) - 1):
        if not check_nt(str(seq1[i]), str(seq2[i])):
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


def n2_context(seq1, seq2, increment=1.0, codon_sites=[]):
    increment = float(increment)
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    # loop through sequence, if one is different, make the key the one before and one after
    neighboring_nuc_dict = {}
    # change_dict = {}
    for i in range(2, min(len(seq1), len(seq2)) - 2):
        if not check_nt(str(seq1[i]), str(seq2[i])):
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
        if not check_nt(str(seq1[i]), str(seq2[i])):
            continue
        if i not in codon_sites:
            continue
        if not check_4fold(seq1[i - 2:i + 1]):
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
        if not check_nt(str(seq1[i]), str(seq2[i])):
            continue
        if i not in codon_sites:
            continue
        if not check_4fold(seq1[i - 2:i + 1]):
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
        if not check_nt(str(seq1[i]), str(seq2[i])):  # makes sure its ATGC not X or -
            continue
        if i not in codon_sites:
            continue
        if not check_4fold(seq1[i - 2:i + 1]):
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


"""
Create function to iterate and store each codon site as a dictionary with 
    key being codon site and value is a tuple
    tuple is (seq1-name_seq2-name, codon_site CTCG) --> C is the third codon position 
    and we are keeping track of the next codon
"""


def site_changes(seq1, seq2, branch_length, codon_sites=[]):
    """
    :param seq1: sequence of node
    :param seq2: sequence of child
    :param branch_length: length of child to node
    :param codon_sites: list of codon sites
    :return: a dictionary with 1178: [('TCTT', 'ATCT'), 'TC', False, branch_length] -
            site, context, change, if site context same, branch length
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    # site[2] = [(CTCG,CTGG),nuc_change,codon_dif=T or F]
    site_dict = {}
    for i in range(2, min(len(seq1), len(seq2)) - 1):
        # if you have something like CTCG CATG --> FALSE because there was a change in sequence/codon context
        codon_good = True
        if not check_nt(str(seq1[i]), str(seq2[i])):
            continue
        if i not in codon_sites:
            continue
        if not check_4fold(seq1[i - 2:i + 1]):
            continue
        # check this
        neighbor_key1 = str(seq1[i - 1]) + '_' + str(seq1[i + 1])  # key = X_X of sequence 1
        neighbor_key2 = str(seq2[i - 1]) + '_' + str(seq2[i + 1])  # key = X_X of sequence 2
        if neighbor_key1 != neighbor_key2:
            codon_good = False
            # can add continue back in if you want to ignore the ones that aren't same codon position
            # continue
        nuc_change = str(seq1[i]) + str(seq2[i])
        site_dict[i] = [(seq1[i - 2:i + 2], seq2[i - 2:i + 2]), nuc_change, codon_good, round(branch_length, 4)]
    #print(site_dict)
    return site_dict
