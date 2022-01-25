from treetracer import TreeTracer
from treetracer_functions import n0_context, n1_context, n2_context, fourfold_n0_context, fourfold_n1_context, fourfold_n2_context


tree_obj = TreeTracer('/Users/cindyruan/Documents/Morton-Research/2021-research/edited_poaceae_rbcl_alignment_txt_phyml_tree.txt', '/Users/cindyruan/Documents/Morton-Research/2021-research/rbcl.fasta')
print('fourfold')

tree_obj.trace_tree_function(fourfold_n1_context, branch_length=False)
tree_obj.print_cumulative_matrices()
'''tree_obj.trace_tree_function(fourfold_n0_context, branch_length=False)
tree_obj.print_cumulative_matrices()'''
'''print('normal')
tree_obj.trace_tree_function(n1_context, branch_length=False)
tree_obj.print_cumulative_matrices()
print('\nn0')
tree_obj.trace_tree_function(n0_context, branch_length=False)
tree_obj.print_cumulative_matrices()'''