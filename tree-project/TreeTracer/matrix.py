"""
Class
"""
import numpy as np
import itertools

def create_nt_dict():
    change_dict = {}
    nucleotides = ['A', 'T', 'G', 'C']
    combos = list(itertools.product(nucleotides, nucleotides))
    # initialize dictionary
    for i in combos:
        key = str(i[0] + i[1])
        change_dict[key] = 0
    return change_dict


class Matrix:

    def __init__(self, cumulative_matrix_dict: dict = {}, all_site_dict: dict = {}):
        self.cumulative_matrix_dict = cumulative_matrix_dict
        self.all_site_dict = all_site_dict
        self.site_dict = {}

    def ngt0_matrix(self, input_dict = None):
        all_matrices = {}
        print(self.cumulative_matrix_dict)
        mat_dict = input_dict or self.cumulative_matrix_dict
        for key in mat_dict:
            arr_zero = np.full((4, 4), 0, dtype=int)
            matrix_dict = mat_dict[key]
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

    """
    Create a function that makes a matrix for each site
    Input is the dictionary of self.site_changes_dict
        changes dictionary of each site --> key is node2:Zea and value is [(CTCG,CTGG),nuc_change,codon_same=T or F]
    """

    def site_matrix(self, specific_site: int = None):
        # print(self.all_site_dict)
        site_dict = {}  # site_dict[5] = {'AA':3, 'TT':5....}
        if specific_site:
            for key in self.all_site_dict.keys():
                print(self.all_site_dict[key][5])
        else:
            for key in self.all_site_dict.keys():
                for site in self.all_site_dict[key].keys():
                    # initialize dictionary inside site_dict
                    if site not in site_dict.keys():
                        site_dict[site] = create_nt_dict()
                    change_key = self.all_site_dict[key][site][1]
                    # make sure its true or don't count
                    if self.all_site_dict[key][site][2]:
                        site_dict[site][change_key] += 1
        print(site_dict)
        self.site_dict = site_dict
        #self.ngt0_matrix(site_dict)
        return #self.ngt0_matrix(site_dict)


    """
    Create a function that calculates the off diagonals
    Take into account transitions and transversions
    return total sum of off diagonals and sum of transitions and transversions
    """
    def sum_off_diagonal(self):
        off_diagonals_sum = {}
        self.site_dict
        for site in self.site_dict.keys():
            total = 0
            for k in self.site_dict[site]:
                if k[0]!=k[1]:
                    total += self.site_dict[site][k]
            if total >0:
                off_diagonals_sum[site] = total
        print(off_diagonals_sum)
        return off_diagonals_sum
