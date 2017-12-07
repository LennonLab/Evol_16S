from __future__ import division
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
import  matplotlib.pyplot as plt

mydir = os.path.expanduser("~/github/Evol_16S/")

def CV_KDE(oneD_array, expand = 1000):
    # use model_selection instead, grid search is deprecated
    # remove +/- inf
    oneD_array = oneD_array[np.logical_not(np.isnan(oneD_array))]
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.logspace(0.1, 5.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(oneD_array[:, None])
    x_grid = np.linspace(np.amin(oneD_array), np.amax(oneD_array), 10000)
    # add nothing to the end of grid and pdf so you can get a nice looking kde
    kde = grid.best_estimator_
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple


def makePlot():
    df = pd.read_csv(mydir + 'data/phylogeny/HYPHYMP/out_clean.txt', sep = '\t')
    values = np.sort(df.Rate.values)
    getKDE = CV_KDE(values)
    fig = plt.figure()
    #fig.suptitle('Distribution of site-specific subsitution rates \n in the 16S ribosomal RNA gene', fontsize=20)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('Distribution of site-specific subsitution rates \n in the 16S ribosomal RNA gene', fontsize=20)

    plt.plot(getKDE[0], getKDE[1],color = 'b', linestyle = '-')
    ax.fill_between(getKDE[0], 0, getKDE[1], color='b', alpha = 0.60)
    ax.set_xlim([0, 50])
    ax.set_ylim([0, 0.07])
    plt.xlabel('Substitution rate', fontsize = 18)
    plt.ylabel('Frequency', fontsize = 18)
    #ax.text(1.8, 0.075, 'Distribution of site-specific subsitution \n rates in the 16S ribosomal RNA gene', fontsize=20)
    output =  mydir + 'figs/rateKDE.png'
    fig.tight_layout()
    plt.savefig(output, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

makePlot()
