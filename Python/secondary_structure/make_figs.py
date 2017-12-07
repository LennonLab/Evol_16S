from __future__ import division
import os
import pandas as pd
import  matplotlib.pyplot as plt

mydir = os.path.expanduser('~/github/Evol_16S/')

def hist():
    df = pd.read_csv(mydir + 'data/secondary_structure/no_pseudoknot/mfe.txt', sep = '\t')

    mfe_clean = [x for x in df.MFE.values if x < 0]
    fig = plt.figure()

    plt.hist(mfe_clean, normed=True, bins=30)
    plt.ylabel('Frequency')
    plt.xlabel('Free energy, ' + r'$\Delta G$')

    out_path = mydir + 'figs/secondary_structure/delta_g.png'
    plt.savefig(out_path, bbox_inches='tight',  dpi = 600)
    plt.close()

def delta_g_vs_p():
    df = pd.read_csv(mydir + 'data/secondary_structure/no_pseudoknot/mfe.txt', sep = '\t')
    df_clean = df.loc[df['MFE'] < 0]
    fig, ax = plt.subplots()
    ax.plot(df_clean.MFE.values, df_clean.Percent_bad.values, marker='o', alpha = 0.8, \
            linestyle='', ms=12)
    fig.savefig(mydir + 'figs/test_bad.png', bbox_inches='tight',  dpi = 600)
    plt.close()


delta_g_vs_p()
