"""
Class
"""
import numpy as np


class Matrix:

    def __init__(self, cumulative_matrix_dict: dict):
        self.cumulative_matrix_dict = cumulative_matrix_dict

    def ngt0_matrix(self):
        all_matrices = {}
        for key in self.cumulative_matrix_dict:
            arr_zero = np.full((4, 4), 0, dtype=int)
            matrix_dict = self.cumulative_matrix_dict[key]
            matrix_arr = arr_zero
            conversion = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
            for k in matrix_dict:
                val = matrix_dict[k]
                i = k[0]
                j = k[1]
                matrix_arr[conversion[i]][conversion[j]] = float(val)
            print(key, '\n', matrix_arr)
            all_matrices[key] = matrix_arr
        return all_matrices

    def n0_matrix(self):
        arr_zero = np.full((4, 4), 0, dtype=int)
        matrix_dict = self.cumulative_matrix_dict
        matrix_arr = arr_zero
        conversion = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
        for k in matrix_dict:
            val = matrix_dict[k]
            i = k[0]
            j = k[1]
            matrix_arr[conversion[i]][conversion[j]] = float(val)
        print(matrix_arr)
        return matrix_arr
