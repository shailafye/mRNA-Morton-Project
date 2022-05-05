giimport itertools


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
            nuc_change = str(seq1[i]) + str(seq2[i])
            change_dict[nuc_change] += 1
        elif seq1[i] != seq2[i]:
            nuc_change = str(seq1[i]) + str(seq2[i])
            change_dict[nuc_change] += 1
    return change_dict

globdict = {}
def n1_context(seq1, seq2,increment=1.0):
    global globdict
    increment = float(increment)
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
            neighbor_key = str(seq1[i-1]) + '_' + str(seq1[i+1]) # key = X_X of sequence 1
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
                if neighbor_key not in globdict:
                    globdict[neighbor_key] = change_dict
            neighboring_nuc_dict[neighbor_key][nuc_change_key] += increment
            globdict[neighbor_key][nuc_change_key] += increment

    return neighboring_nuc_dict

# print(globdict)
# n1_context('ATGCAGACGCTCAATGTGCAGATATACCA', 'ATGCACACGCTCAATGTACAGATTTACCA')
# #print(n1_context('ATGCAGACGCTCAATGTGCAGATATACCA', 'ATGCACACGCTCAATGTACAGATTTACCA'))
# print('\nglobe_dict:',globdict,'\n')
# n1_context('ATGCAGACGCTCAATGTGCAGATATACCA', 'ATGGAGACGCTCAATGTCCAGATCTAACA')
# print('\nglobe_dict:',globdict,'\n')
"""
ATGCAGACGCTCAATGTGCAGATATACCA
ATGCAcACGCTCAATGTaCAGATtTACCA
AA --> gc
TC --> ga
TT --> at
"""


def n2_context(seq1, seq2):
    """


    """


>>>>>>> 542c96290ecdd8a31f279baf4aadc84a2c60eaf8
