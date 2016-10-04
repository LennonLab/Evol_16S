from __future__ import division
import generatePop as gp
#import quantifyPop as qp
import numpy as np
import os, argparse

#mydir = os.path.expanduser('/N/dc2/projects/Lennon_Sequences/2016_Evol_16S/')
mydir = os.path.expanduser('~/github/Evol_16S/')

#pop_size = 100
mutation_rate = 0.01
#generations = 10

siteRates =  np.loadtxt(open(mydir + 'Tree/data/rates/siteRate.txt',"rb"),delimiter=",",skiprows=1)
rates = siteRates[:,0]
seq_length = len(rates)
print seq_length

def runSimulation(pop_size):
    OUT = open(mydir + 'Simulation/data/simulation_out.fna', 'w')
    generations = 1
    pop = gp.generate_pop(pop_size, seq_length)
    pop_out = gp.simulate(pop, mutation_rate, generations, seq_length, pop_size, rates)
    pop_list = []

    for key, value in pop_out.iteritems():
        pop_list += value * [key]
    for i, sequence in enumerate(pop_list):
        header = '>'+str(i)
        print>> OUT, header
        print>> OUT, sequence

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Run analyses and/or generate figures.")
    parser.add_argument('-n', type = int, default = 10, help = "How many sequences")

    args = parser.parse_args()
    N = args.n

    runSimulation(N)
