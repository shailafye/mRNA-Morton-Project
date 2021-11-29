import os

files = os.listdir('./Sorted/')
for f in files:
    if f.endswith('.txt'):
        f1 = f.split('.')[0]
        gene = f1.split('_')[0]
        print(gene)