import os

def parse_genome(path, file_name):
    #import dna_sequence as dna
    f1 = file_name.split('.')[0]
    species = f1.split('_')
    #if len(species) == 0:
        #continue
    genus = species[0]
    specific = species[1]

    infile = open(path + file_name, 'r')
    lines = infile.readlines()
    infile.close()

    gene_sequences = {}

    for line in range(0, len(lines), 2):
        gene_name = lines[line].strip(' ')[1:].strip()
        gene_sequences[gene_name] = lines[line+1].strip()
        #gene_sequences[gene_name] = dna.DNASequence(lines[line+1].strip(), True)

    #if os.path.exists(file_name):
        outfile = open('./Sorting/' + gene_name + '_nc.txt', 'a+')
        outfile.write('>' + genus + " " + specific + '\n' + gene_sequences[gene_name] + '\n')


files = os.listdir('./parsed_files/')
for f in files:
    if f.endswith('nc.txt'):

        parse_genome('./parsed_files/', f)