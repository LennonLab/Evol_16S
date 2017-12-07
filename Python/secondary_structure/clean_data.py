from __future__ import division
import os, re, math
import pandas as pd

mydir = os.path.expanduser('~/github/Evol_16S/')


'''
data acquired from the Comparatice RNA Web Site and Project

http://www.rna.ccbb.utexas.edu/DAT/3C/SBPI/index.php
'''

def merge_files(bracket, nopct, sample_name):
    br_clean = ''
    for i, line in enumerate(open(bracket, 'r')):
        if i < 4:
            continue
        br_clean += line.strip()
    ct_clean = pd.read_csv(nopct, skiprows=[0,1,2,3,4], header=None, sep = r"\s*", engine='python')
    ct_clean.columns = ['Index', 'Base', 'Index_minus_one', 'Index_plus_one', 'Paired_base' , 'Nat_num']
    paired_base = ct_clean.Paired_base.values
    br_clean_test = ''.join([ '.' if x == '.' else '|' for x in br_clean])
    ct_clean_test = ''.join([ '.' if x == 0 else '|' for x in paired_base])
    if len(br_clean) != len(paired_base):
        print('Files have unequal length:' + sample_name)
        print(len(br_clean), len(paired_base))
    else:
        ct_clean['Bracket'] = pd.Series(list(br_clean), index=ct_clean.index)
        ct_clean.to_csv(mydir + 'data/secondary_structure/no_pseudoknot/SB.bacilli.merged/' + sample_name+ '.txt' , sep = '\t')


def iterate_files():
    directory = mydir + 'data/secondary_structure/no_pseudoknot/'
    if not os.path.exists(directory +  'SB.bacilli.merged/'):
        os.makedirs(directory +  'SB.bacilli.merged/')
    for br in os.listdir(directory + 'SB.bacilli.bracket/'):
        sample_name = br.rsplit('.', 1)[0]
        ct_path = directory + 'SB.bacilli.nopct/' + sample_name + '.nopct'
        if br.endswith(".bracket") and os.path.exists(ct_path):
            merge_files(os.path.join(directory + 'SB.bacilli.bracket/', br), ct_path, sample_name)


#br = mydir + 'data/secondary_structure/no_pseudoknot/SB.bacilli.bracket/AB021192.bracket'

#no = mydir + 'data/secondary_structure/no_pseudoknot/SB.bacilli.nopct/AB021192.nopct'

#merge_files(br, no, 'AB021192')

iterate_files()
