from __future__ import division
import os, re, math
import pandas as pd
#import  matplotlib.pyplot as plt

mydir = os.path.expanduser('~/github/Evol_16S/')

def Fig1():
    L = 1600



class classFASTA:

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print("Not in FASTA format.")

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list

def numberVarSites():
    read_FASTA = classFASTA(mydir + 'data/phylogeny/alignment/nmicrobiol201648-s7_clean_noEuks.txt')
    seqs = [x[1] for x in read_FASTA.readFASTA()]
    N = len(seqs)
    sites = zip(*seqs)
    sites_filtered = [x for x in sites if (x.count('-') / N) < 0.9 ]
    sites_filtered_same = [len(set(x)) for x in sites if (x.count('-') / N) < 0.9 ]
    print sites_filtered_same

def testGaussian():
    # in units kcal per mol^-1
    L_phys = 1600
    L_eff = 100
    rescale = math.sqrt(L_phys/ L_eff)
    sample = np.random.normal(0,1 * rescale)

numberVarSites()
