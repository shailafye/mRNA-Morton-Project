"""
Class
"""
import numpy as np
import seaborn as sns
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from statistics import mean, median, mode, stdev, variance


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
        self.all_change_site_dict = {}
        self.final_changes_df = pd.DataFrame

    def ngt0_matrix(self, input_dict=None):
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
        # self.ngt0_matrix(site_dict)
        return  # self.ngt0_matrix(site_dict)

    """
    Create a function that calculates the off diagonals
    Take into account transitions and transversions
    return total sum of off diagonals and sum of transitions and transversions
    """

    def sum_off_diagonal(self):
        off_diagonals_sum = {}
        # self.site_dict
        for site in self.site_dict.keys():
            total = 0
            for k in self.site_dict[site]:
                if k[0] != k[1]:
                    total += self.site_dict[site][k]
            if total > 0:
                off_diagonals_sum[site] = total
        print(off_diagonals_sum)
        return off_diagonals_sum

    def site_changes_dictionary_to_df(self, total_branch_length):
        """
        Attempting to do this but create a dictionary
        Key = site_and_context
            Key = site and context ==> 5_('TCAC', 'TCAC') string
        Values = total_branch_length, change_nt, branch_context
            total_branch_length = sum of all branch lengths of site occurrence corresponding to the key
            proportion_branch_length = what proportion of branches have that context
            change_nt = a list of all the nt changes ==> [AA, AT, AA, AG] eg
            branch_context = a list of all the branches with this context ==> [Node1:Oryza, Node1:Node2]
        """
        print("matrix site changes dictionary function")

        site_change_dict = {}
        condense_site_changes = {}
        for key in self.all_site_dict:
            cur_dict = self.all_site_dict[key]
            # need to now calculate for each site, IF same codon context then add to site matrix?
            for site in cur_dict:
                # skip if a site doesn't have same codon context --> CHECK THIS...
                if not cur_dict[site][2]:
                    continue
                # if the site already has a corresponding codon context that is same then add to that,
                # if different, then create a new row
                cod1 = cur_dict[site][0][0]
                cod2 = cur_dict[site][0][1]
                site_and_context = str(site) + "_" + str(cod1) + "_" + str(cod2)
                # if row already exists then update!
                b_len = round(cur_dict[site][3], 4)
                change_nt = [cur_dict[site][1]]
                if site_and_context in site_change_dict.keys():
                    # update length --> add together
                    site_change_dict[site_and_context][0] += round(b_len, 4)
                    # update branch length proportion to total tree
                    site_change_dict[site_and_context][1] \
                        = round(site_change_dict[site_and_context][0] / total_branch_length, 3)
                    # append the change_nt
                    site_change_dict[site_and_context][2] += change_nt
                    # append to list of branches
                    site_change_dict[site_and_context][3] += [key]
                else:
                    # add dictionary element
                    b_len_prop = round(b_len / total_branch_length, 3)
                    site_change_dict[site_and_context] = [b_len, b_len_prop, change_nt, [key]]

        df = pd.DataFrame.from_dict(site_change_dict, orient='index')
        # print(df.sort_index(axis = 0))
        self.all_change_site_dict = df.sort_index(axis=0)
        return df.sort_index(axis=0)

    def condensed_change_site_dict(self, min_tree_prop=0.6):
        # print(self.all_change_site_dict)
        df = self.all_change_site_dict
        condense_site_changes = {}
        # check if site and cod1 and cod2 equal each other then add to dictionary
        # 1130_GT_G: NumChanges, {AA:2,GA:1,GG:3}, {AA:0.103,GA:0.236,GG:0.659}, total proportion, number of flanking GC
        for i, row in df.iterrows():
            # print('index: ', i, row)
            i = i.split("_")
            site = i[0]
            cod1 = i[1]
            cod2 = i[2]
            # print(site, cod1,cod2)
            context1 = cod1[0:2] + "_" + cod1[3]
            context2 = cod2[0:2] + "_" + cod2[3]
            if context1 == context2:
                key = str(site) + "_" + context1
                # print(key)
                # create new instance of that key in dictionary - intialize
                if key not in condense_site_changes:
                    condense_site_changes[key] = [0, {}, {}, 0, 0]
                nuc = row[2][0]
                # print(nuc)
                # now update the dictionary
                # update changes -- if like AT or not same
                if nuc[0] != nuc[1]:
                    condense_site_changes[key][0] += len(row[2])
                    # print(nuc, len(row[2]))
                # update first inside dictionary with count  --> {AA:2,GA:1,GG:3}
                condense_site_changes[key][1][nuc] = len(row[2])
                # update second inside dictionary with proportions --> {AA:0.103,GA:0.236,GG:0.659}
                condense_site_changes[key][2][nuc] = row[1]
                # update total branch length
                condense_site_changes[key][3] += row[1]
                # num of surround G or C
                # count occurrences of G and C in flank context
                flanked_bases = context1[1:]
                # print(flanked_bases)
                counter = flanked_bases.count('G')
                counter += flanked_bases.count('C')
                condense_site_changes[key][4] = counter
        # print(condense_site_changes)
        # remove sites with a minimum proportion --> parameter
        # go through condense_site_changes and make new dictionary with tree proportions greater than whats specified
        site_changes_proportion = {}
        for key in condense_site_changes.keys():
            if condense_site_changes[key][3] > min_tree_prop:
                site_changes_proportion[key] = condense_site_changes[key]

        # condensed_df = pd.DataFrame.from_dict(condense_site_changes, orient='index')
        condensed_df = pd.DataFrame.from_dict(site_changes_proportion, orient='index')
        condensed_df.columns = ['total_changes', 'change_type_count', 'change_tree_proportion',
                                'tree_branch_proportion', 'GC_flanking_count']
        condensed_df.to_csv('/Users/shailafye/Documents/Morton-Research/2021-research/condensed-df-169taxa.csv',
                            index=True)
        self.final_changes_df = condensed_df
        return condensed_df

    def graph_freq_distribution_matplotlib(self):
        # fix the
        plt.style.use("seaborn")
        # Later in the code
        # self.final_changes_df
        # use matplotlib to create graphs of frequency distribution of total changes
        # will need to separate on GC flanking count

        fig, ax = plt.subplots(figsize=(6, 4))
        self.final_changes_df['weighted_changes'] = self.final_changes_df['total_changes'] / self.final_changes_df[
            'tree_branch_proportion']
        changes = self.final_changes_df['weighted_changes']
        # plt.hist(x, bins=17)
        # Plot histogram
        changes.plot(kind="hist", density=True, bins=18)  # change density to true, because KDE uses density
        # Plot KDE
        changes.plot(kind="kde")
        # x.plot(kind="hist", density=False, bins=17)  # change density to true, because KDE uses density
        # X #
        ax.set_xlabel("Changes")
        ax.set_xlim(0)
        # Y #
        # Overall #
        ax.set_title("Frequency Distribution of Changes per Site", size=15, pad=10)
        ax.grid(False)
        plt.show()
        return

    def graph_freq_distribution_seaborn(self):
        self.final_changes_df['weighted_changes'] = self.final_changes_df['total_changes'] / self.final_changes_df[
            'tree_branch_proportion']

        def plot_changes(gc=None, df=pd.DataFrame):
            # df.loc[df['favorite_color'] == 'yellow']
            gc_context = "All"
            if gc is not None:
                print(df.loc[df['GC_flanking_count'] == gc])
                df = df.loc[df['GC_flanking_count'] == gc]
                gc_context = str(gc)
            changes = df['weighted_changes']
            sns.set_style("white")
            plt.figure(figsize=(10, 7), dpi=80)
            p = sns.histplot(changes, color="dodgerblue", label="Compact", kde=True, bins=20)
            p.set_xlabel("Number of Changes", fontsize=13)
            p.set_ylabel("Sites", fontsize=13)
            title_gc = 'Frequency of Changes per Sites ' + gc_context + ' GC Context'
            p.set_title(label=title_gc, fontsize=19)
            plt.show()
            statistical_test(changes)

        # gc content = all
        plot_changes(gc=None, df=self.final_changes_df)
        # gc content = 0
        plot_changes(gc=0, df=self.final_changes_df)
        # gc content = 1
        plot_changes(gc=1, df=self.final_changes_df)
        # gc content = 2
        plot_changes(gc=2, df=self.final_changes_df)


def statistical_test(all_changes=pd.DataFrame):
    # skew is variance over mean
    print(all_changes)
    var = variance(all_changes)
    print('variance:', var)
    me = mean(all_changes)
    print('mean:', me)
    print('skew:', var / me)
