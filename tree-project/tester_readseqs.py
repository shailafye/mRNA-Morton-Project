# test file
# read fasta into dictionary

seq_file = "/Users/shailafye/Documents/Morton-Research/2021-research/data/atpA/all_atpA_seqs.txt"



from Bio import SeqIO
seq_dict = {rec.id : rec.seq for rec in SeqIO.parse("/Users/shailafye/Documents/Morton-Research/2021-research/data/atpA/all_atpA_seqs.txt", "fasta")}
seq_dict_final = {}
for key in seq_dict:
    seq_temp = str(seq_dict[key])
    seq_dict_final[key] = seq_temp

print(seq_dict_final)
