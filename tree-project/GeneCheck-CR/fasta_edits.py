import os

def species_edit(path, file_name):
    f1 = file_name.split('.')[0]
    #gene = f1.split('_')[0]
    print(file_name)
    infile = open(path + file_name, 'r')
    lines = infile.readlines()
    infile.close()

    #if f1.endswith('cds'):
    outfile = open(
        '/Users/cindyruan/Documents/Morton-Research/2021-research/edited_fasta/edited_' + file_name, 'w')
    #elif f1.endswith('nc'):
        #outfile = open(
            #'/Users/cindyruan/Documents/Morton-Research/2021-research/edited_fasta/edited_' + gene + '_nc.txt', 'w')


    for line in range(0, len(lines), 2):
        species_name = lines[line].strip(' ')[1:].strip()
        genus = species_name.split(' ')[0][:2]
        specific = species_name.split(' ')[1]
        sequence = lines[line + 1].strip()
        #outfile = open('/Users/cindyruan/Documents/Morton-Research/2021-research/edited_fasta/edited_' + gene + '_cds.txt', 'a+')
        outfile.write('>' + genus + "." + specific + '\n' + sequence + '\n')

        #outfile = open('/Users/cindyruan/Documents/Morton-Research/2021-research/edited_fasta/edited_' + gene + '_nc.txt', 'a+')
        #outfile.write('>' + genus + " " + specific + '\n' + sequence + '\n')
    outfile.close()
    #return gene_name, gene_sequence, gene

files = os.listdir('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_fasta_files')
for f in files:
    if f.endswith('.txt'):
        species_edit('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_fasta_files/', f)


