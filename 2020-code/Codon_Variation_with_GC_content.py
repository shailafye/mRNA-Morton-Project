#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 14:59:43 2021

@author: shailafye
"""

'''

    This file goes through the excel with the change data and
    finds the average across the window --> make sure the window size
    is the same as the window size set for the MFE generating program
    
    This program also accounts for surrounding GC content within the 
    degenerate site.
    It saves 8 files total
    
    #0 GC, 1 GC, 2 GC, aggregate
    
    *changes written as variable name or column name can mean the same to GC
    surrounding content
    
    Double check that the excel document has the correct column names as the
    program assumes
    
    IMPORTANT
    

'''



import pandas as pd
import math
import numpy as np
                            
gene = 'psaa'
#rbcl_change_data = '/Users/shailafye/Documents/RNA Secondary Structures/psab/psab_change_data.xlsx'
gene_change_data = '/Users/shailafye/Documents/RNA Secondary Structures/' +gene+ '/' +gene+ '_change_data.xlsx'

change_data = pd.read_excel(gene_change_data,
                            sheet_name='Sheet1',
                            header=0,
                            index_col=False,
                            keep_default_na=True
                            )

#adding two new columns
#codon and next codon
change_data['codon'] = change_data['Anc. Codon']
change_data['next codon'] = change_data['Anc. Codon'].shift(-1)

print(change_data)
#print(change_data['Anc. Codon'][3])
#print(change_data['Anc. Codon'][3][1])
length = len(change_data)
#print(change_data['next codon'].str[0:1])

change_data['first'] = change_data['codon'].str[1:2]
#change c and g to 1 and then will have last column be 0 1 or 2 --> sum of G or C
#print(change_data)
#change_data['one'] = np.where(change_data['first'] == 'G', 1, 0)
change_data['one'] = [1 if x == 'C'or x == 'G' else 0 for x in change_data["first"]]
#df['color'] = ['red' if x == 'Z' else 'green' for x in df['Set']]

change_data['third'] = change_data['next codon'].str[0:1]
change_data['three'] = np.where(((change_data['third'] == 'G') | (change_data['first'] == 'C')), 1, 0)
change_data["total_GC"] = change_data['three'] + change_data['one']



change_data = change_data[0:length] #NEED TO CHANGE THIS LENGTH!!

#change_data = change_data[['Codon #', "Anc. Codon", "Changes", "Deg"]]
change_data = change_data[["Changes", "Deg", "total_GC"]]
change_data['Changes'] = change_data['Changes'].fillna(0)
change_data['Deg'] = change_data['Deg'].fillna(0)
#print(change_data)
#print(change_data[0:50])

window_avg_4= [[],[],[],[]] #0 GC, 1 GC, 2 GC, aggregate
window_avg_2= [[],[],[],[]] #0 GC, 1 GC, 2 GC, aggregate
deg_value = 4.0
offset = 5
i = 0

print(change_data)

windowsize = 40 
length = 3*len(change_data)
while i<(length-windowsize):
    end = i + windowsize
    #print(end - i) --> this should be window size
    if(i%3 != 0):
        start_codon_num = math.ceil(i/3) 
    else:
        start_codon_num = math.ceil(i/3)
    if(end%3 != 0):
        end_codon_num = round(end/3)-1    
    else:
        end_codon_num = (end/3)-1 
    #get average for span of rows based on condition if deg = 2 or 4 
    #variable above to choose
    subset = change_data[int(start_codon_num):int(end_codon_num)+1]
    #print(subset)

    subset_deg = subset.loc[subset['Deg'] == deg_value]

    #zero changes
    new_1 = subset_deg.loc[subset['total_GC'] == 0]
    deg_avg = new_1['Changes'].mean() 
    window_avg_4[0].append(deg_avg)
    
    #one change
    new_1 = subset_deg.loc[subset['total_GC'] == 1]
    deg_avg = new_1['Changes'].mean() 
    window_avg_4[1].append(deg_avg)

    #two changes
    new_1 = subset_deg.loc[subset['total_GC'] == 2]
    deg_avg = new_1['Changes'].mean() 
    window_avg_4[2].append(deg_avg)
    
    #aggregate
    deg_avg = subset_deg['Changes'].mean() 
    window_avg_4[3].append(deg_avg)
    
    #for 2

    subset_deg2 = subset.loc[subset['Deg'] == 2.0]
    #print(subset_deg2)
    #zero changes
    new_2 = subset_deg2.loc[subset['total_GC'] == 0]
    deg_avg_2 = new_2['Changes'].mean() 
    window_avg_2[0].append(deg_avg_2)
    
    #one change

    new_2 = subset_deg2.loc[subset['total_GC'] == 1]
    deg_avg_2 = new_2['Changes'].mean() 
    window_avg_2[1].append(deg_avg_2)
    
    #two changes

    new_2 = subset_deg2.loc[subset['total_GC'] == 2]
    deg_avg_2 = new_2['Changes'].mean() 
    window_avg_2[2].append(deg_avg_2)
    
    #aggregate

    deg_avg_2 = subset_deg2['Changes'].mean() 
    window_avg_2[3].append(deg_avg_2)
    
    i = i+offset

for i in range(len(window_avg_4)):
    window_avg_4[i] = [0 if math.isnan(x) else x for x in window_avg_4[i]]
    window_avg_2[i] = [0 if math.isnan(x) else x for x in window_avg_2[i]]
    
    
#print(len(window_avg_4[0]))

#print(window_avg_2[3]) --> can manually check with excell file



'''
    SAVING ALL OF THE FILES OF 0,1,2 and aggregate GC content on each side of the third codon

'''

f4 = open("/Users/shailafye/Documents/RNA Secondary Structures/" +gene+  "/window_avg_4_aggregate.txt", "w")
f2 = open("/Users/shailafye/Documents/RNA Secondary Structures/" +gene+  "/window_avg_2_aggregate.txt", "w")
for i in window_avg_2[3]:
    f2.write(str(i)+"\n")
for i in window_avg_4[3]:
    f4.write(str(i)+"\n")
    

f4 = open("/Users/shailafye/Documents/RNA Secondary Structures/" +gene+  "/window_avg_4_0GC.txt", "w")
f2 = open("/Users/shailafye/Documents/RNA Secondary Structures/" +gene+  "/window_avg_2_0GC.txt", "w")
for i in window_avg_2[0]:
    f2.write(str(i)+"\n")
for i in window_avg_4[0]:
    f4.write(str(i)+"\n")


f4 = open("/Users/shailafye/Documents/RNA Secondary Structures/" +gene+  "/window_avg_4_1GC.txt", "w")
f2 = open("/Users/shailafye/Documents/RNA Secondary Structures/" +gene+  "/window_avg_2_1GC.txt", "w")
for i in window_avg_2[1]:
    f2.write(str(i)+"\n")
for i in window_avg_4[1]:
    f4.write(str(i)+"\n")


f4 = open("/Users/shailafye/Documents/RNA Secondary Structures/" +gene+  "/window_avg_4_2GC.txt", "w")
f2 = open("/Users/shailafye/Documents/RNA Secondary Structures/" +gene+  "/window_avg_2_2GC.txt", "w")
for i in window_avg_2[2]:
    f2.write(str(i)+"\n")
for i in window_avg_4[2]:
    f4.write(str(i)+"\n")

