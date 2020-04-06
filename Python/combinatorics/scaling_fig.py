from __future__ import division
import os
import numpy as np
import  matplotlib.pyplot as plt


mydir = os.path.expanduser('~/github/Evol_16S/')

def scaling_fig(k = 4):
    x1 = np.linspace(1, 30, num=30)
    y_drift = 4 ** x1
    y_drift_std= np.sqrt(4 ** (2*x1))
    print y_drift_std
    x2 = np.linspace(1, 1600, num=1600)
    y_select = x2 * (np.log(k) + 0.5772)
    y_select_var = np.sqrt(y_select)

    fig, ax = plt.subplots()

    ax.plot(x1,y_drift, ls="-", c="#FF6347", lw = 3, label='Neutral limit')
    ax.plot(x2,y_select, ls="-", c="#87CEEB", lw  = 3, label='Strong-selection limit')
    plt.axhline(y = 13440, c = 'grey', linestyle = '--', lw = 3, label='Average number of \n16S rRNA substitutions')
    #plt.fill_between(x1, y_drift + y_drift_std, y_drift - y_drift_std , facecolor="#FF6347", alpha=0.5)
    ax.set_yscale('log', basey=10)
    ax.set_xscale('log', basex=10)
    ax.set_ylim([1, max(y_drift)])
    ax.set_xlim([1, 1600])
    ax.legend(loc='upper right', fontsize = 12)
    ax.set_xlabel('Sequence length', fontsize=16)
    ax.set_ylabel('Average number of substitutions \nto explore sequence space', fontsize=14)
    fig.savefig(mydir + 'figs/pop_gen_scale.png', bbox_inches='tight',  dpi = 600)
    plt.close()




scaling_fig()
