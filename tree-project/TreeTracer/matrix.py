"""
Class
"""


class Matrix:

    def __init__(self, total_matrix_dict, sep_matrix_dict, sum_matrix_dict):
        self.total_matrix_dict = total_matrix_dict
        self.sep_matrix_dict = sep_matrix_dict
        self.comb_matrix_dict = sum_matrix_dict

    def n1_matrix(self):
        for key in self.total_matrix_dict:
            arr_zero = np.full((4, 4), 0, dtype=int)
            matrix_dict = self.total_matrix_dict[key]
            matrix_arr = arr_zero
            conversion = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
            for k in matrix_dict:
                val = matrix_dict[k]
                i = k[0]
                j = k[1]
                matrix_arr[conversion[i]][conversion[j]] = float(val)
            print(key, '\n', matrix_arr)
        return True
