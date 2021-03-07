#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 15:49:33 2021

@author: shailafye
"""


"""
    IMPORTANT
    *may need to hit run in the beggining 2 times to clear and write files
    Two separate chunks, can hit run once then go to terminal and run:
        FIRST:
        RNAfold < sequence_windows.txt > output_original.txt
            File name: output_original.txt
        SECOND:
        RNAfold < all_random_seq.txt > random_output.txt
            File name: random_output.
    Then hit run on this file again


"""




#IF THIS ISN'T FIRST RUN, THEN RUN THIS COMMAND in terminal/director
    # rm all_random_seq.txt sequence_windows.txt output_original.txt random_output.txt

import random
from collections import defaultdict

path = '/Users/shailafye/Documents/RNA Secondary Structures/'
file_name = '/psbc/Nicotiana tabacum_dna psbc.txt' #change for each gene #check correct folder
#file_name = 'rbcl/rbcl consenus (.5) sequence.txt'

input_file = open(path+file_name, 'r')

sequence = input_file.read()
input_file.close()

sequence = sequence.strip()
print(len(sequence))

# Below is code to generate the windows --> window size is 40 and offset by 5
windows = []
offset = 5
i = 0
windowsize = 40 
while i<(len(sequence)-windowsize):
    windows.append(sequence[i:windowsize+i]) #0-41
    #print(i)
    #print(windowsize+i)
    #print(len(sequence[i:windowsize+i]))
    i = i+offset
    
print(len(sequence))

print(len(windows))


# This is code for putting the original 40 window sequences into a file
f1 = open('/Users/shailafye/Tutorial/Progs/sequence_windows.txt', 'w')
count = 0
for i in windows:
    f1.write(i + "\n")
    count = count +1


##after
print(count)
print(len(windows))


'''
    Randomize the sequence within each window
    
'''


cur = windows[0]
#print(''.join(random.sample(cur, len(cur))))
cur = cur.replace("T", "U")
cur = list(cur)
random.shuffle(cur)
r = ''.join(cur)
print(r)
print(sorted(r.replace("T", "U")) == sorted(windows[0].replace("T", "U")))
all_seq_random = defaultdict(list)

'''
Now go through those sliding windows and randomize the sequence
Create 100 random sequences
This is stored in a dictionary all_seq_random with the key being the original 
40bp window and the value being a list of 100 randomized variations of that sequence
'''


for cur_seq in windows:
    j = 0
    random_sequence = []
    seq = ''
    cur = windows[j]
    for i in range(100):
        cur = list(cur)
        random.shuffle(cur)
        seq = ''.join(cur)
        #seq = ''.join(random.sample(cur, len(cur)))
        random_sequence.append(seq)
    all_seq_random[cur_seq].append(random_sequence)
    j = j+1
    
    
# This is to check to make sure the length is the same as windows length
print(len(windows))
print(len(all_seq_random))


#Save randomized sequences to a file now
only_random = []
j=0
cur = ""
for cur_seq in windows:
    seq = ''
    cur = windows[j]
    for i in range(100):
        cur = list(cur)
        random.shuffle(cur)
        seq = ''.join(cur)
        #seq = ''.join(random.sample(cur, len(cur)))
        only_random.append(seq)
    j = j + 1
    
 

len(only_random)

#check to make sure this is true
#print(only_random[0])
sorted(only_random[0].replace("T", "U")) == sorted(windows[0].replace("T", "U"))


f8 = open('/Users/shailafye/Tutorial/Progs/all_random_seq.txt','w')
c = 0
for i in only_random:
    f8.write(i + "\n")
    c= c+1
print("Lines:", c)
#make sure this number matches

    
'''
GO TO TERMINAL

    After running to get MFE, do this to get output
    This is just for the first 40 window one for the original sequence, 
        not the randomized one.
    Run in C and write to output_original.txt
    
    FIRST:
        RNAfold < sequence_windows.txt > output_original.txt
            File name: output_original.txt
    SECOND:
        RNAfold < all_random_seq.txt > random_output.txt
            File name: random_output.txt

Now get back the 100*windows of randomized sequences MFE
- Next parse it
     - get just the values of MFE
     - split by every 100
     - check if that is the same sequence though
     

'''

#NONRANDOMIZED

f2 = open('/Users/shailafye/Tutorial/Progs/output_original.txt', 'r')
#extract the MFE from the file and put it into mfe list
count = 0
linecount = 0
mfe = []
value = ''
for line in f2:
    linecount = linecount +1
    #print(line)
    if ' (' in line:
        value = float(line[-8:-2])
        mfe.append(value)
        count = count + 1

f = open('/Users/shailafye/Documents/RNA Secondary Structures/original_mfe.txt', 'w')
for i in mfe:
    f.write(str(i) + "\n")



#RANDOMIZED

f1 = open('/Users/shailafye/Tutorial/Progs/random_output.txt', 'r')

#saving just the sequences to check the randomize

mfe_random_seq = []
for line in f1:
    #print(line)
    if ' (' not in line:
        mfe_random_seq.append(line.replace('\n', ""))

print(len(mfe_random_seq))
print(only_random[100])
print(mfe_random_seq[100])
sorted(mfe_random_seq[2030]) == sorted(windows[20].replace("T", "U"))
#make sure this is true



#now getting MFE of randomized sequences
    # Get average of each 100 MFE
    # Below just sends all of the mfe values for the all 100*windows
f1 = open('/Users/shailafye/Tutorial/Progs/random_output.txt', 'r')

count = 0
linecount = 0
mfe_random = []
value = 0.0
for line in f1:
    linecount = linecount +1
    if ' (' in line:
        value = float(line[-8:-2])
        mfe_random.append(value)
        count = count + 1


# This is the code to get the averages of the 100 MFE of the random sequences
i = 0
total = 0;
count = 0
avg = []
while i <= (len(mfe_random)-100):
    avg.append(sum(mfe_random[i:i+100])/100)
    i = i + 100
    count = count +1;

print(len(avg))

f1 = open('/Users/shailafye/Documents/RNA Secondary Structures/random_avg.txt', 'w')
for i in avg:
    f1.write(str(i) + "\n")













