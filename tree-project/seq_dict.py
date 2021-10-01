def open_seq(file_name):
    gene_sequences = {}
    infile = open(file_name, 'r')
    lines = infile.readlines()
    infile.close()

    for line in range(0, len(lines), 2):
        gene_name = lines[line].strip()[1:]
        gene_sequences[gene_name] = lines[line+1].strip()

    return gene_sequences

ancestral = open_seq('./ancestralseqs.fasta.nodes.txt')
rbcl = open_seq('./grass_rbcl.txt')

ancestral.update(rbcl)
print(ancestral.keys())