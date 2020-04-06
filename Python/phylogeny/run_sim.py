from __future__ import division
import os
import numpy as np
import pandas as pd

mydir = os.path.expanduser("~/GitHub/Evol_16S/")


class generate_pop:

    def __init__(self, N, L, t):
        self.L = L
        self.N = N
        self.t = t

    def get_seq(self, sample):
        seq = list('A' * self.L)
        rates = sample.Rate.values
        rates = rates / sum(rates)
        index = sample.index.values
        #sites = sample.Site.values
        sub_sites = np.random.choice(index, p = rates, size = self.t)
        sub_sites = list(set(sub_sites))
        for sub_site in  sub_sites:
            site_row = sample.loc[[sub_site]]
            nuc_rates = site_row[['A', 'C', 'G', 'T']].values[0]
            new_nuc = np.random.choice(['A', 'C', 'G', 'T'], p=nuc_rates)
            seq[sub_site] = new_nuc
        seq = "".join(seq)
        return seq


    def get_pop(self):
        df_path = mydir + 'data/phylogeny/HYPHYMP/rates_freqs.txt'
        df = pd.read_csv(df_path, sep = '\t')
        OUT_path = mydir + 'data/phylogeny/simulation/sim_align/test.txt'
        OUT = open(OUT_path, 'w')
        for n in range(self.N):
            print(n)
            sample = df.sample(n = self.L)
            n_seq = self.get_seq(sample.reset_index(drop=True))
            OUT.write('>Sequence_' + str(n) + '\n')
            #[line[i:i+70] for i in range(0, len(line), 70)]
            for i in range(0, self.L, 70):
                OUT.write(n_seq[i:i+70] + '\n')

            #rates = sample.Rate.values
            #rates = rates / sum(rates)
            #print(rates)
        OUT.close()


generate_pop(N = 10, L = 1425, t = 4907).get_pop()
