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
    """