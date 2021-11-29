import itertools


def n0_context(seq1, seq2):
    change_dict = {}
    nucleotides = ['A', 'T', 'G', 'C']
    combos = list(itertools.product(nucleotides, nucleotides))
    # initialize dictionary
    for i in combos:
        key = str(i[0] + i[1])
        change_dict[key] = 0
    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] == seq2[i]:
            continue
        elif seq1[i] != seq2[i]:
            nuc_change = str(seq1[i]) + str(seq2[i])
            change_dict[nuc_change] += 1
    return change_dict


def n1_context(seq1, seq2):
    """
    A_A, A_T, A_G, A_C,
    T_A, T_T, T_G, T_C,
    G_A, G_T, G_G, G_C,
    C_A, C_T, C_G, C_C

    neighboring_nuc_dict = {A_A: {'AA': 0, 'AT': 7, 'AG': 22, 'AC': 5...}, A_T: {'AA': 9, 'AT': 2, 'AG': 12...},

    ATGCAGACgCTCGACaC
    ATGCAGACtCTCGACgC

    C_C: {'AG': 1, 'AT': 0, 'AG': 0, 'GT': 1...}

    Edge case:
    - two or more mutations in a row -- what is neighbor
    - mutation at beginning or end
    """
    # loop through sequence, if one is different, make the key the one before and one after
    neighboring_nuc_dict = {}
    #change_dict = {}
    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] == seq2[i]:
            continue
        elif seq1[i] != seq2[i] and i > 0:
            nuc_change_key = str(seq1[i]) + str(seq2[i]) #AT, TG
            #change_dict[nuc_change] += 1
            # check to see if neighboring nt exists in dictionary
            neighbor_key = ''
            if seq1[i-1] is None:
               neighbor_key = '_' + str(seq1[i+1])
            elif seq1[i+1] is None:
               neighbor_key = str(seq1[i-1]) + '_'
            else:
                neighbor_key = str(seq1[i-1]) + '_' + str(seq1[i+1]) # key = X_X of sequence 1

            other_neighbor_key = str(seq2[i - 1]) + '_' + str(seq2[i + 1])
            if neighbor_key != other_neighbor_key:
                continue
            if neighbor_key not in neighboring_nuc_dict:
                # initialize a new dictionary inside of neighboring_nuc_dict
                change_dict = {}
                nucleotides = ['A', 'T', 'G', 'C']
                combos = list(itertools.product(nucleotides, nucleotides))
                # initialize dictionary
                for k in combos:
                    key = str(k[0] + k[1])
                    change_dict[key] = 0
                neighboring_nuc_dict[neighbor_key] = change_dict
            neighboring_nuc_dict[neighbor_key][nuc_change_key] += 1
    return neighboring_nuc_dict

print(n1_context('ATGCAGACGCTCAATGTGCAGATATACCA', 'CTGCACACGCTCAATGTACAGATTTACCA'))
"""
ATGCAGACGCTCAATGTGCAGATATACCA
cTGCAcACGCTCAATGTaCAGATtTACCA
AA --> gc
TC --> ga
TT --> at
"""

def n2_context(seq1, seq2):
    neighboring_nuc_dict = {}
    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] == seq2[i]:
            continue
        elif seq1[i] != seq2[i] and i > 0:
            nuc_change_key = str(seq1[i]) + str(seq2[i]) #AT, TG
            #change_dict[nuc_change] += 1
            # check to see if neighboring nt exists in dictionary
            '''if seq1[i-1] is None:
                neighbor_key = '_' + str(seq1[i+1])
            elif seq1[i+1] is None:
                neighbor_key = str(seq1[i-1]) + '_' '''
            neighbor_key = str(seq1[i-2] + str(seq1[i-1]) + '_' + str(seq1[i+1]) + str(seq1[i+2])) # key = X_X of sequence 1
            other_neighbor_key = str(seq2[i - 2] + str(seq2[i - 1]) + '_' + str(seq2[i + 1]) + str(seq2[i + 2]))
            if neighbor_key != other_neighbor_key:
                continue
            if neighbor_key not in neighboring_nuc_dict:
                # initialize a new dictionary inside of neighboring_nuc_dict
                change_dict = {}
                nucleotides = ['A', 'T', 'G', 'C']
                combos = list(itertools.product(nucleotides, nucleotides))
                # initialize dictionary
                for k in combos:
                    key = str(k[0] + k[1])
                    change_dict[key] = 0
                neighboring_nuc_dict[neighbor_key] = change_dict
            neighboring_nuc_dict[neighbor_key][nuc_change_key] += 1
    return neighboring_nuc_dict

#print(n2_context('ATGCAGACGCTCAATGTGCAGATATACCA', 'ATGCACTCGCTCAATGTACAGATTTACCA'))
"""
ATGCAGACGCTCAATGTGCAGATATACCA
ATGCActCGCTCAATGTaCAGATtTACCA
AA --> gc
TC --> ga
TT --> at
"""