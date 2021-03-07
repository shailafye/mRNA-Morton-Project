#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 16:11:26 2021

@author: shailafye
"""

sequences = {}

                                #CHANGE FILE NAME DOWN BELOW ALSO!!!!
f = open("/Users/shailafye/Documents/RNA Secondary Structures/psaa/psaa alignment good.txt")
seq_length = 0
for line in f:
    line_data = line.split()
    if(len(line_data) == 3 ):
        #check to make sure doesn't start with an *
        try:
            sequences[line_data[0]]+=line_data[1]
        except:
            sequences[line_data[0]] = line_data[1]
 
seq_list = []
for s in sequences:
    #print(s) #key
    #print(sequences[s])
    seq_length = len(sequences[s])
    seq_list.append(sequences[s])
    
print(seq_list)



def result_nuc(nuc_dict, num, param):
    cutoff = int(num*param)
    for k in nuc_dict:
        if(nuc_dict.get(k)>cutoff):
            return str(k)
    return "N"


consensus = ""
nucleotides = {"A":0, "T":0, "G":0, "C":0}
param = .5

for site in range(seq_length):
    nucleotides = {"A":0, "T":0, "G":0, "C":0}
    base_count = 0 
    for j in range(len(seq_list)):
        cur_nuc = seq_list[j][site]
        try:
            nucleotides[cur_nuc]+=1
            base_count += 1
        except:
            pass
    consensus = consensus + result_nuc(nucleotides,base_count,param)
    #print(nucleotides)
    #nucleotides.clear()

    
print(consensus)
f1 = open('/Users/shailafye/Documents/RNA Secondary Structures/psaa/psaa consenus (.5) sequence.txt', 'w')
f1.write(consensus)

#write a second function that returns the degeneracy code 
# =============================================================================
# def deg_nuc(nuc_dict, num, param):
#     iupac ={"GA":"R", "CT":"Y"}#... continue
#     #make sure in the right order
#     #AG not GA
#     cutoff = int(num*param)
#     
#     for k in nuc_dict:
#         if(nuc_dict.get(k)>cutoff):
#             yield str(k)
# 
# =============================================================================

#yield sends back values --> at the recieving end, collect all and then find common
#single degenerate code








