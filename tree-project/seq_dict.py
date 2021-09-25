def open_seq(file_name):
    gene_sequences = {}
    infile = open(file_name, 'r')
    lines = infile.readlines()
    infile.close()

    for line in range(0, len(lines), 2):
        gene_name = lines[line].strip()[1:]
        gene_sequences[gene_name] = lines[line+1].strip()

    print(gene_sequences)

open_seq('./ancestralseqs.fasta.nodes.txt')
open_seq('./grass_rbcl.txt')