#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 13:22:46 2020

@author: bmorton
"""
    
class DNASequence:
    
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'R':'Y', 'Y':'R', 'r':'y', 'y':'r', 'N':'N', 'n':'n'}
    nucleotides = ['A','C','G','T', 'N', 'Y', 'R']
    bases = ['A','C','G','T']

    codon_translation = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
                           "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
                           "TAT": "Y", "TAC": "Y", "TGT": "C", "TGC": "C",
                           "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L",
                           "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P",
                           "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q",
                           "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R",
                           "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I",
                           "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T",
                           "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K",
                           "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R",
                           "AGG": "R", "GTT": "V", "GTC": "V", "GTA": "V",
                           "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A",
                           "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E",
                           "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
                           "GGG": "G", "TAG": "*", "TAA": "*", "TGA": "*"}

    genetic_code = {'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 
                'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Y': ['TAT', 'TAC'], 
                'C': ['TGT', 'TGC'], 'W': ['TGG'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 
                'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 
                     'CGG', 'AGA', 'AGG'], 'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 
                     'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'],
                     'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
                     'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
                     '*': ['TAG', 'TAA', 'TGA']}
    
    def __init__(self, s = 'N', isCoding = False):
        self.sequence = s
        self.is_Coding = isCoding
        
        
    def countBasesInSequence(self):
        count_num = {}
        for n in DNASequence.bases:
           count_num[n] = self.sequence.upper().count(n)
        return count_num
       
    def countBasesInSequenceWithTotal(self):
        count_num = {}
        total = 0
        for n in DNASequence.bases:
            count_num[n] = self.sequence.upper().count(n)
            total += self.sequence.upper().count(n)
        return count_num, total
    
    def percentGC(self):
        n, t = self.countBasesInSequenceWithTotal()
        try:
            return (n['G']+n['C'])/t
        except:
            return -1.0
    
    
    
    