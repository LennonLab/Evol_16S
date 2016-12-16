from __future__ import division
import numpy as np
import pandas as pd
import os

mydir = os.path.expanduser('~/github/Evol_16S/')

def parseMothur(Ns = [10, 100, 1000, 10000]):
    df = pd.DataFrame()
    for N in Ns:
        directory = mydir + 'Simulation/data/N' + str(N) +  '/'
        list_S = []
        for filename in os.listdir(directory):
            print N, filename
            if 'phylip' in filename:
                IN = pd.read_csv(directory + filename, sep ='\t')
                S = IN.iloc[0,1]
                list_S.append(S)
        data = pd.DataFrame({'N' + str(N): list_S})
        df = pd.concat([df,data], ignore_index=True, axis=1)
    df.to_csv(path_or_buf = mydir + 'Simulation/data/test_SbyS.txt', sep=' ')

parseMothur()
