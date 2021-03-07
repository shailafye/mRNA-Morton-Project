#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 10:32:10 2021

@author: shailafye
"""

'''
    IMPORTANT: RUN TWO SECTIONS SEPARATELY TO CREATE A FILE codon_change.txt 
    first
    --  F9 to run highlighted code

    IMPORTANT: repalce any blank codons in Excel with XXX

    This file takes in the excel file as with the corresponding gene 
    with the change/variation per codon data and extracts and outputs the
    resulting codon sequence
    Make sure to adjust after looking at an alignement with this sequence and
    the species you are interested in looking at
    Edit the excel file and add or remove any extra lines
    
    New file: codon_change.txt is created in corresponding folder of gene
    New file: codon_change_sequence.txt in FASTA format
    
    *no need to run 'get_codon_change'
    
'''


import pandas as pd

gene = 'psbd'

gene_change_data = '/Users/shailafye/Documents/RNA Secondary Structures/' +gene+ '/' +gene+ '_change_data.xlsx'

change_data = pd.read_excel(gene_change_data,
                            sheet_name='Sheet1',
                            header=0,
                            index_col=False,
                            keep_default_na=True
                            )

change_data = change_data[['Anc. Codon']]
#print(change_data)

#seq = []

f8 = open( '/Users/shailafye/Documents/RNA Secondary Structures/' +gene+ '/codon_change.txt', 'w')
for index, row in change_data.iterrows():
    j = (row['Anc. Codon'])
    print(j)
    #seq.append(j.strip())
    f8.write(str(j)+"\n")
    



'''


RUN SEPARATELY ABOVE AS ONE RUN AND BELOW AS ANOTHER RUN!!!!



'''



f = open('/Users/shailafye/Documents/RNA Secondary Structures/' +gene+ '/codon_change.txt', 'r')

seq = []

for line in f:
    #print(line)
    #seq += line
    seq.append(line.strip())
    
#print(seq)

#print(len(seq[0]))

print((seq[3][0:1]))

for i in range(len(seq)):
    if(len(seq[i])>3):
        print(len(seq[i]))
        print(seq[i])
        #seq[i] = seq[i][2:5] #this is choosing the last three
        o = seq[i][0:1]
        l = seq[i][3:5]
        seq[i] = o+l
        print(seq[i])
        
        
print(seq)

sequence = ''.join(seq)

#print(sequence)
sequence = ">codon_seq\n" + sequence
print(sequence)

f8 = open( '/Users/shailafye/Documents/RNA Secondary Structures/' +gene+ '/codon_change_sequence.txt', 'w')
f8.write(sequence)
