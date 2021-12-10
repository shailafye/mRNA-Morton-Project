import os
'''problem_species is a list of species we don't want to include in
the analysis'''
problem_species = ['bambusa oldhamii']
'''gene_check function takes in a file that should contain all the sequences of a given gene
in fasta format
a dictionary is created where the key is the name of the gene and the value is the list
of the species for which a sequence exists in the fasta file

gene length must be greater than the value of min_length.
species_lst will exclude any species within the problem_species list.'''
def gene_check(path, file_name):
    min_length = 300
    f1 = file_name.split('.')[0]
    gene = f1.split('_')[0]
    infile = open(path + file_name, 'r')
    lines = infile.readlines()
    infile.close()
    '''gene_sequences dictionary: keys are gene names, values are the list of species 
    that have that gene sequence.
    for each gene, the list of species is contained in species_lst before it gets 
    stored as the value.
    species_lst gets returned separately to be utilized in other applications'''
    gene_sequences = {}
    species_lst = []

    for line in range(0, len(lines), 2):
        gene_name = lines[line].strip(' ')[1:].strip()
        gene_sequence = lines[line+1].strip()
        if len(gene_sequence) > min_length and gene_name not in problem_species:
            species_lst.append(gene_name)
    gene_sequences[gene] = species_lst

    return(gene, species_lst)

'''species_check is sent a file that contains all the genes for a 
particular species in fasta format
it returns a list of all the genes that are in that file for which
a sequence exists'''
def species_check(path, file_name):
    f1 = file_name.split('.')[0]
    species = f1.split('_')
    genus = species[0]
    specific = species[1]
    infile = open(path + file_name, 'r')
    lines = infile.readlines()
    infile.close()
    gene_lst = []

    for line in range(0, len(lines), 2):
        gene_name = lines[line].strip(' ')[1:].strip()
        gene_lst.append(gene_name)
        #gene_sequence = lines[line+1].strip()
    return(genus, specific, gene_lst)


files = os.listdir('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_fasta_files')
files2 = os.listdir('/Users/cindyruan/Documents/Morton-Research/2021-research/parsed_files')
outfile1 = open('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_analysis/list_of_genes.txt', 'w')
outfile2 = open('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_analysis/list_of_species.txt', 'w')

cds_list = [] #list of all genes in the provided folder
nc_list = [] #list of all nc regions in the provided folder
for g in files:
    if g.endswith('cds.txt'):
        gene, species_cds = gene_check('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_fasta_files/', g)
        cds_list.append(gene)
    if g.endswith('nc.txt'):
        nc_region, species_nc = gene_check('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_fasta_files/', g)
        nc_list.append(nc_region)
print(cds_list)

'''the genes_missing_in_species and ncr_missing_in_species dictionaries
are used for finding all the genes and nc regions that are missing from 
each species in the files given'''
genes_missing_in_species = {}
ncr_missing_in_species = {}
cds_species = []
nc_species = []
lack_nc_regions = []
lack_genes = []
print('\n\n\n')
for f in files2:
    if f.endswith('cds.txt'):
        genus, specific, species_gene_lst = species_check('/Users/cindyruan/Documents/Morton-Research/2021-research/parsed_files/', f)
        species = genus + ' ' + specific
        outfile1.write('\n' + species + '\n')
        for element in species_gene_lst:
            outfile1.write(element + '\n')
        if species not in problem_species:
            cds_species.append(species) # check genus + specific isn't in problem species
            outfile2.write(species + '\n')
        lack_genes = []
        for gene in cds_list:
            if gene not in species_gene_lst:
                lack_genes.append(gene)
        genes_missing_in_species[species] = lack_genes
    if f.endswith('nc.txt'):
        genus, specific, species_nc_lst = species_check('/Users/cindyruan/Documents/Morton-Research/2021-research/parsed_files/', f)
        species = genus + ' ' + specific
        if species not in problem_species:
            nc_species.append(species)  # check genus + specific isn't in problem species
        lack_nc_regions = []
        for nc in nc_list:
            if nc not in species_nc_lst:
                lack_nc_regions.append(nc)
        ncr_missing_in_species[species] = lack_nc_regions

'''Does the opposite of above: for each gene or nc region, 
this for loop generates a list of species that do not have sequence data
for that gene.'''
species_without_gene_seqs = {}
species_without_ncr = {}
lack_nc_species = []
lack_species = []
for h in files:
    if h.endswith('cds.txt'):
        gene, species_cds = gene_check('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_fasta_files/', h)
        lack_species = []
        for species in cds_species:
            if species not in species_cds and species not in problem_species:
                lack_species.append(species)
        species_without_gene_seqs[gene] = lack_species
    if h.endswith('nc.txt'):
        gene, species_nc = gene_check('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_fasta_files/', h)
        lack_nc_species = []
        for species in nc_species:
            if species not in species_nc and species not in problem_species:
                lack_nc_species.append(species)
        species_without_ncr[gene] = lack_nc_species

genes_in_all_species = []
ncr_in_all_species = []
for gene in species_without_gene_seqs:
    if len(species_without_gene_seqs[gene]) < 10 and len(species_without_gene_seqs[gene]) != 0:
        pass
        #print(gene, len(species_without_gene_seqs[gene]), species_without_gene_seqs[gene])
    elif len(species_without_gene_seqs[gene]) == 0:
        print('gene in all species: ', gene)
        genes_in_all_species.append(gene)
for ncr in species_without_ncr:
    if len(species_without_ncr[ncr]) < 10 and len(species_without_ncr[ncr]) != 0:
        pass
        #print(ncr, len(species_without_ncr[ncr]), species_without_ncr[ncr])
    elif len(species_without_ncr[ncr]) == 0:
        print('nc region in all species: ', ncr)
        ncr_in_all_species.append(ncr)


'''print('following are the genes that exist in all species: ')
print(genes_in_all_species)
print('following are the ncrs that exist in all species: ')
print(ncr_in_all_species)'''

'''create a folder: poaceae_analysis
put two text files in folder: list of genes present in all species line by line
    and list of species used in poaceae: species_cds not including problem species line by line'''

outfile1.close()
outfile2.close()
