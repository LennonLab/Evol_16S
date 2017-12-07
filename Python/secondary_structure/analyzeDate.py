from __future__ import division
import os, re
import pandas as pd

mydir = os.path.expanduser('~/github/Evol_16S/')

def read_bpseq(file_path):
    name = '>' + file_path.split('/')[-1]
    IN = pd.read_csv(file_path, skiprows=4, sep = ' ', header = None)
    IN.columns = ['Position', 'Nucleotide', 'base_pair']
    seq_list = IN['Nucleotide'].values.tolist()
    seq_string = ''.join(seq_list)
    return (name, seq_string)


def bpseq_to_fasta():
    '''This function takes all the secondary structure files and
    outputs the sequences as a single fasta.'''
    OUT = open(mydir + 'data/secondary_structure/structure.fasta','w+')
    data_dir = mydir + 'data/secondary_structure/structure.bpseq/'
    for filename in os.listdir(data_dir):
        if filename.endswith(".bpseq"):
            file_path = os.path.join(data_dir, filename)
            data = read_bpseq(file_path)
        # [-1] because the last string is empty
        seq_string_split = re.findall('.{0,80}', data[1])[:-1]
        print>> OUT, data[0]
        for x in seq_string_split:
            print>> OUT, x
    OUT.close()

def clean_align():
    IN = mydir + 'data/secondary_structure/structure.align'
    OUT = open(mydir + 'data/secondary_structure/structure_clean.align','w+')
    with open(IN) as f:
        for line in f:
            if line[0] == '>':
                line = line.replace(".", "_")
                line = line.replace("-", "_")
            else:
                line = line.replace(".", "-")
                line = line.replace("-", "-")
            line = line.strip('\n')
            print>> OUT, line
    OUT.close()


clean_align()
