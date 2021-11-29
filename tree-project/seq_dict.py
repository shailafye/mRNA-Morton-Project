'''from Bio import SeqIO
ancestral = {rec.id : rec.seq for rec in SeqIO.parse("./ancestralseqs.fasta.nodes.txt", "fasta")}
rbcl = {rec.id : rec.seq for rec in SeqIO.parse("./grass_rbcl.txt", "fasta")}
ancestral.update(rbcl)'''


def open_seq(file_name):
    import dna_sequence as dna
    gene_sequences = {}
    infile = open(file_name, 'r')
    lines = infile.readlines()
    infile.close()

    '''for line in range(0, len(lines), 2):
        gene_name = lines[line].strip(' ')[1:].split(" ")[0]
        gene_sequences[gene_name] = lines[line+1].strip()
        #gene_sequences[gene_name] = dna.DNASequence(lines[line+1].strip(), True)'''

    '''{'Node': 'AGCTAGTCGATGC', 'Node2': 'TGACTAGCTA'} '''
    for line in range(len(lines)):
        if lines[line].startswith('>'):
            gene_name = lines[line].strip()[1:].split()[0]
            #gene_name = gene_name.strip(' ')
            #gene_name = gene_name.split(" ")[0]
            sequence = ''

            for num in range(line + 1, len(lines)):
                if lines[num].startswith('>'):
                    line = num - 1
                    break
                sequence += lines[num].strip()
            gene_sequences[gene_name] = sequence


    return gene_sequences

ancestral = open_seq('./ancestralseqs.fasta.nodes.txt')
rbcl = open_seq('./grass_rbcl.txt')
ancestral.update(rbcl)
print(ancestral)
#names = open_seq('./MultipleLineTester.txt')
'''for s in names:
    print(s)
    print(ancestral[s].sequence[:20])
    print(ancestral[s].percentGC())'''

'''seq = open_seq('./MultipleLineTester.txt')
print(seq)'''


