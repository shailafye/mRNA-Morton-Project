# import numpy as np

# def check_seq(seq1, seq2):
#     # find shorter and do that range
#     change_dict = {'AA': 0, 'CG': 0, 'GC': 0}
#     for i in range(len(seq1)):
#         if seq1[i] == seq2[i]:
#             continue
#         elif seq1[i] != seq2[i]:
#             # want GT or GC or CG
#             nuc_change = str(seq1[i]) + str(seq2[i])
#             change_dict[nuc_change] += 1
#     return change_dict
#
#
# print(check_seq('ATCGATCGATCG', 'ATCGATGGATCC'))


# def two_d_arr() -> np.array:
#     arr_zero = np.full((5, 5), 0, dtype=int)
#     return arr_zero
#
# print(two_d_arr())
