from __future__ import division
import pandas as pd
import os

mydir = os.path.expanduser("/N/dc2/projects/Lennon_Sequences/2016_Evol_16S/")

def parseMothur(M = 65, Ns = [10, 100, 1000, 10000]):
    df = pd.DataFrame()
    for N in Ns:
        directory = mydir + 'Combinatorics/data/M' + str(M) + '/N' + str(N) +  '/'
        list_S = []
        for filename in os.listdir(directory):
            if 'phylip' in filename:
                IN = pd.read_csv(directory + filename, sep ='\t')
                S = IN.iloc[0,1]
                list_S.append(S)
        data = pd.DataFrame({'N' + str(N): list_S})
        df = pd.concat([df,data], ignore_index=True, axis=1)
    df.to_csv(path_or_buf = mydir + 'Combinatorics/data/M' + str(M) + '/' + 'M' + str(M) +'_SbyS.txt', \
        sep=' ')

parseMothur()
