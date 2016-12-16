from __future__ import division
import numpy as np
import pandas as pd
import os, argparse
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
import  matplotlib.pyplot as plt
from scipy.misc import comb

mydir = os.path.expanduser('~/github/Evol_16S/')

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
            print "Not in FASTA format."

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


def CV_KDE(oneD_array):
    # remove +/- inf
    oneD_array = oneD_array[np.logical_not(np.isnan(oneD_array))]
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.logspace(0.1, 5.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(oneD_array[:, None])
    x_grid = np.linspace(np.amin(oneD_array), np.amax(oneD_array), 10000)
    kde = grid.best_estimator_
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple

def plotSeqSubRate():
    siteRates =  np.loadtxt(open(mydir + 'Tree/data/rates/siteRate.txt',"rb"),delimiter=",",skiprows=1)
    rates = siteRates[:,0]
    print min(rates)
    x = range(1, len(rates) + 1)
    fig, ax = plt.subplots()
    plt.plot(x, rates, marker='o')
    ax.set_xlabel('16S rRNA position', fontsize = 16  )
    ax.set_ylabel('Substitutions (rate * tree length)', fontsize = 16)
    plt.savefig(mydir + 'Simulation/figs/SeqSubRates1.png', dpi=600,)
    plt.close()

    fig, ax = plt.subplots()
    plt.plot(x, rates, marker='o')
    #fig, ax = plt.subplots()
    plt.ylim(0,6)
    ax.set_xlabel('16S rRNA position', fontsize = 16  )
    ax.set_ylabel('Substitutions (rate * tree length)', fontsize = 16)
    plt.savefig(mydir + 'Simulation/figs/SeqSubRates2.png', dpi=600,)
    plt.close()




def plotSubRate():
    siteRates =  np.loadtxt(open(mydir + 'Tree/data/rates/siteRate.txt',"rb"),delimiter=",",skiprows=1)
    rates = siteRates[:,0]
    rates = [x/np.mean(rates) for x in rates ]
    #print len(rates)
    #KDE = CV_KDE(rates)
    fig, ax = plt.subplots()
    #ax.plot(KDE[0], KDE[1], linewidth=3, alpha=0.5, label='bw=%.2f' % KDE[2])
    #ax.hist(rates, 100, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
    ax.hist(rates, 30, fc='gray', histtype='stepfilled', alpha=0.5, normed=False)
    ax.set_xlabel('Substitutions (rate * tree length)', fontsize = 16  )
    ax.set_ylabel('Number of sites', fontsize = 16)
    ax.text(2, 1650, 'Site-specific substitutions in the 16S rRNA gene ', fontsize=15)
    #ax.legend(loc='upper left')
    plt.savefig(mydir + 'Simulation/figs/SubRates.png', dpi=600,)
    plt.close()


def binomFig():
    fig, ax = plt.subplots()
    N = 1000
    x = []
    y = []
    for j in range(1, N+1):
        binSample = comb(N,j)*(0.25)**j*(0.75)**(N-j)
        x.append(j)
        y.append(binSample)
        #print binSample
    plt.plot(x, y, marker='o')
    plt.savefig(mydir + 'Simulation/figs/Binom.png', dpi=600,)
    plt.close()

def plotDivergence(obligate = False):
    # 3.8 billion years / divergence per 1 million years
    time_steps = 3800
    X = np.logspace(0, np.log10(time_steps), num = 1000, base=10.0)
    # percent per 50 million years
    fig, ax = plt.subplots()
    if obligate == True:
        minDiv = (0.025 / 0.03 )
        maxDiv = (0.091 / 0.03)
        meanDiv = (minDiv + maxDiv) / 2
        y_min = np.asarray([(2**(minDiv * x)) for x in X ])[:-300]
        y_max = np.asarray([(2**(maxDiv * x)) for x in X ])[:-300]
        y_mean = np.asarray([(2**(meanDiv * x)) for x in X ])[:-300]
        X = np.logspace(0, np.log10(time_steps), num = 1000, base=10.0)[:-300]
        name = mydir + '/Simulation/figs/DivergenceObligate.png'

    else:
        minDiv = (0.025 / 0.03 )  / 25
        maxDiv = (0.091 / 0.03) / 25
        meanDiv = (minDiv + maxDiv) / 2
        y_min = np.asarray([(2**(minDiv * x)) for x in X ])
        y_max = np.asarray([(2**(maxDiv * x)) for x in X ])
        y_mean = np.asarray([(2**(meanDiv * x)) for x in X ])
        name = mydir + '/Simulation/figs/DivergenceFree.png'

    y_min = np.asarray([(2**(minDiv * x)) for x in X ])
    y_max = np.asarray([(2**(maxDiv * x)) for x in X ])
    y_mean = np.asarray([(2**(meanDiv * x)) for x in X ])

    ax.plot(X, y_mean, lw=2, color='black', alpha = 0.9)
    #for k in zip(X, y_mean):
    #    print k
    ax.fill_between(X, y_min, y_max, facecolor='#87CEEB', alpha=0.5)
    plt.axhline(y=1000000000000,  linestyle='dashed', color = 'darkgrey', linewidth=4)
    ax.grid()
    ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)
    ax.set_xlabel('Time (millions of years)', fontsize=20)
    if obligate == True:
        ax.set_xlim([0,500])
    else:
        ax.set_xlim([0,3800])
    ax.set_ylabel('Species', fontsize=20)
    fig.savefig(name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def gcContent():
    IN = mydir + 'Tree/data/alignment/nmicrobiol201648-s7_clean_noEuks.txt'
    classIN = classFASTA(IN)
    seqs = [ x[1] for x in classIN.readFASTA() ]
    sites = zip(*seqs)
    GCs = []
    for site in sites:
        if site.count('-') / len(site) > 0.9:
            continue
        else:
            site1 = [s.replace('-', '') for s in site]
            GC = (site1.count('G') + site1.count('C')) / len(site1)r
            GCs.append(GC)
    GCs = np.asarray(GCs)
    print "Mean = "  + str(np.mean(GCs))
    print "Standard deviation = " + str( np.std(GCs))
    x = range(1, len(GCs) + 1)
    fig, ax = plt.subplots()
    plt.plot(x, GCs, marker='o', alpha = 0.5)
    plt.ylim(0,1)
    ax.set_xlabel('16S rRNA position', fontsize = 16  )
    ax.set_ylabel('GC content', fontsize = 16)
    plt.savefig(mydir + 'Simulation/figs/GC.png', dpi=600,)
    plt.close()


gcContent()
#plotDivergence(obligate = False)
