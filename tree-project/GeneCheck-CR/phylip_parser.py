gene_sequences = {}
infile = open('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_rbcl_alignment.txt', 'r')
file = infile.read()
infile.close()

segments = file.split('\n\n')
species_seg = segments[0].split('\n')[1:]
sequence_segs = segments[1:]
count = 0
for line in species_seg:
    species_name = line[:10]
    full_seq = line[31:]
    for seq in sequence_segs:
        seq_lines = seq.split('\n')
        full_seq += seq_lines[count]
    count += 1
    gene_sequences[species_name] = full_seq

print(gene_sequences)

