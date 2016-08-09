from __future__ import division
import os, argparse, random, math
from collections import Counter
from itertools import product
import scipy
from scipy.stats import chisquare
from Bio import Phylo

mydir = os.path.expanduser("~/github/Evol_16S/")


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


def rename_fasta():
    read_FASTA = classFASTA(mydir + 'data/nmicrobiol201648-s7.txt')
    OUT = open(mydir + 'data/nmicrobiol201648-s7_clean.txt','w+')
    for x in read_FASTA.readFASTA():
        header = x[0].translate(None, '-.')
        header = header.replace(',', '_')
        header = header.replace(':', '_')
        header = header.replace('/', '_')
        header = header.replace("'", '_')
        if header == 'Archaea_Diapherotrites_candidate_division_DUSEL3_archaeon_SCGC_AAA011E11_DUSEL_001_206':
            header = header[:8] + 'DPANN_' + header[8:]
        if header == 'Bacteria_bacterium_JKG1':
            header = 'Bacteria_Chloroflexi_bacterium_JKG1'
        if header == 'Bacteria_unclassified_Bacteria_Thermobaculum_Thermobaculum_terrenum_ATCC_BAA798':
            header = 'Bacteria_Chloroflexi_Thermobaculum_Thermobaculum_terrenum_ATCC_BAA798'
        if 'Eukaryota' in header:
            continue
        print>> OUT, '>' +header, '\n', x[1]

    tree = Phylo.read(mydir + 'data/nmicrobiol201648-s8.txt', 'newick')
    for node in tree.find_clades():
        #Don't print anything if there is no name
        if node.name:
            node.name = node.name.translate(None, '-.')
            node.name  = node.name.replace(',', '_')
            if 'CPpSL4' in node.name:
                node.name = node.name.replace('CPpSL4','CP_pSL4')
            if 'CPEM' in node.name:
                node.name = node.name.replace('CPEM','CP_EM')
            if 'CPDadabacteria' in node.name:
                node.name = node.name.replace('CPDadabacteria','CP_Dadabacteria')
            if 'CPOP9' in node.name:
                node.name = node.name.replace('CPOP9','CP_OP9')
            if 'CPRokubacteria' in node.name:
                node.name = node.name.replace('CPRokubacteria','CP_Rokubacteria')
            if 'CPOP8' in node.name:
                node.name = node.name.replace('CPOP8','CP_OP8')
            if 'CPCD12' in node.name:
                node.name = node.name.replace('CPCD12','CP_CD12')
            if 'CP' in node.name[:12] and 'CP_' not in node.name[:12]  :
                node.name = node.name[:12].replace('CP', 'CP_') + node.name[12:]

    #net = Phylo.to_networkx(tree)
    Phylo.write(tree, mydir + 'data/nmicrobiol201648-s8_clean.txt', 'newick')

import re
def unAlignFasta():
    read_FASTA = classFASTA(mydir + 'data/nmicrobiol201648-s7.txt')
    OUT = open(mydir + 'data/nmicrobiol201648-s7_unAligned.txt','w+')
    headers = []
    for x in read_FASTA.readFASTA():
        header = x[0].translate(None, '-.')
        seq = x[1]
        if 'Eukaryota' in header:
            continue
        seq = seq.replace('-', '')
        header = header.replace(',', '_')
        header = header.replace(':', '_')
        header = header.replace('/', '_')
        header = header.replace("'", '_')
        header = header.translate(None, '!')
        headers.append(header)
        print>> OUT, '>' +header, '\n', seq


unAlignFasta()
