import os
'''eliminate problem species '''
problem_species = ['bambusa oldhamii']
def gene_check(path, file_name):
    f1 = file_name.split('.')[0]
    gene = f1.split('_')[0]
    infile = open(path + file_name, 'r')
    lines = infile.readlines()
    infile.close()
    gene_sequences = {}
    species_lst = []

    for line in range(0, len(lines), 2):
        gene_name = lines[line].strip(' ')[1:].strip()
        gene_sequence = lines[line+1].strip()
        if len(gene_sequence) > 300 and gene_name not in problem_species:
            species_lst.append(gene_name)
    gene_sequences[gene] = species_lst
    '''if len(species_lst) > num/2:
        print(gene)
        print(len(species_lst), num)'''

    return(gene, species_lst)


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


files = os.listdir('./Sorting/')
files2 = os.listdir('./parsed_files/')
'''check for problem species within file name
find number of species you want to use for num
ex: number of species - number of problem species'''

cds_list = []
nc_list = []
for g in files:
    if g.endswith('cds.txt'):
        gene1, species_cds = gene_check('./Sorting/', g)
        cds_list.append(gene1)
    if g.endswith('nc.txt'):
        gene2, species_nc = gene_check('./Sorting/', g)
        nc_list.append(gene2)

lack_genes = []
lack_nc_genes = []
lack_dict = {}
lack_nc_dict = {}
cds_species = []
nc_species = []
print('\n\n\n')
for f in files2:
    if f.endswith('cds.txt'):
        genus, specific, species_gene_lst = species_check('./parsed_files/', f)
        species = genus + ' ' + specific
        if species not in problem_species:
            cds_species.append(species) # check genus + specific isn't in problem species
        lack_genes = []
        for gene in cds_list:
            if gene not in species_gene_lst:
                lack_genes.append(gene)
    lack_dict[species] = lack_genes
    if f.endswith('nc.txt'):
        genus, specific, species_nc_gene_lst = species_check('./parsed_files/', f)
        species = genus + ' ' + specific
        if species not in problem_species:
            nc_species.append(species)  # check genus + specific isn't in problem species
        lack_nc_genes = []
        for gene in nc_list:
            if gene not in species_nc_gene_lst:
                lack_nc_genes.append(gene)
    lack_nc_dict[genus + ' ' + specific] = lack_nc_genes
#print(lack_dict)
#print(lack_nc_dict)


species_lack_dict = {}
lack_species = []
species_nc_lack_dict = {}
lack_nc_species = []
for h in files:
    if h.endswith('cds.txt'):
        gene, species_cds = gene_check('./Sorting/', h)
        lack_species = []
        for species in cds_species:
            if species not in species_cds and species not in problem_species:
                lack_species.append(species)
    species_lack_dict[gene] = lack_species
    if h.endswith('nc.txt'):
        gene, species_nc = gene_check('./Sorting/', h)
        lack_nc_species = []
        for species in nc_species:
            if species not in species_nc and species not in problem_species:
                lack_nc_species.append(species)
    species_nc_lack_dict[gene] = lack_nc_species

zero_gene = []
zero_nc_gene = []
for gene in species_lack_dict:
    if len(species_lack_dict[gene]) < 10 and len(species_lack_dict[gene]) != 0:
        print(gene, len(species_lack_dict[gene]), species_lack_dict[gene])
    elif len(species_lack_dict[gene]) == 0:
        zero_gene.append(gene)
for gene in species_nc_lack_dict:
    if len(species_nc_lack_dict[gene]) < 10 and len(species_nc_lack_dict[gene]) != 0:
        print(gene, len(species_nc_lack_dict[gene]), species_nc_lack_dict[gene])
    elif len(species_nc_lack_dict[gene]) == 0:
        zero_nc_gene.append(gene)
#print(zero_gene)

'''for each gene (key), make list of species that doesn't have it (value)'''
