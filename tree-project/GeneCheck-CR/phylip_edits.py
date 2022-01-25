infile = open('/Users/cindyruan/Documents/Morton-Research/2021-research/poaceae_rbcl_alignment.txt', 'r')
file = infile.read()
infile.close()

outfile = open('/Users/cindyruan/Documents/Morton-Research/2021-research/edited_poaceae_rbcl_alignment.txt', 'w')

species_segment = file.split('\n\n')[0]
species_list = species_segment.split('\n')[1:]
outfile.write(species_segment.split('\n')[0] + '\n')

for species in species_list:
        species = species[:10] + "  " + species[10:]
        outfile.write(species + '\n')

outfile.write('\n')
segments = file.split('\n\n')[1:]
for line in segments:
        outfile.write(line)
        outfile.write('\n\n')


