from __future__ import division
import randomDNA as rd
import numpy as np
import os

mydir = os.path.expanduser('~/github/Evol_16S/')


def SimN():
    siteRates =  np.loadtxt(open(mydir + \
    'Tree/data/rates/siteRate.txt','rb'),delimiter=',',skiprows=1)
    rates = siteRates[:,0]
    probRates = [x/sum(rates) for x in rates]
    L = len(probRates)
    setSub = set(rd.subSeq(10))
    L_set = set(range(L))
    print len(L_set)
    print len(L_set - setSub)

def get_distance(seq_a, seq_b):
    diffs = 0
    length = len(seq_a)
    assert len(seq_a) == len(seq_b)
    for chr_a, chr_b in zip(seq_a, seq_b):
        if chr_a != chr_b:
            diffs += 1
    #return diffs / float(length)
    return diffs

def tajimas_theta(population):
    '''Accepts a single dictionary of a population where keys are haplotypes \
        and values are the counts'''
    # AKA Pi
    pop_size = sum(population.values())
    haplotypes = population.keys()
    haplotype_count = len(haplotypes)
    diversity = 0
    for i in range(haplotype_count):
        haplotype_a = haplotypes[i]
        frequency_a = population[haplotype_a] / float(pop_size)
        for j in range(0, i):
            haplotype_b = haplotypes[j]
            frequency_b = population[haplotype_b] / float(pop_size)
            frequency_pair = frequency_a * frequency_b
            diversity += frequency_pair * get_distance(haplotype_a, haplotype_b)
    return diversity * 2


def testString(L = 10):
    Ns = [10, 100, 1000, 10000]
    for N in Ns:
        test_dict = {}
        for x in range(N):
            string = rd.String(L)
            if string not in  test_dict:
                test_dict[string] = 1
            else:
                test_dict[string] += 1
        #print test_dict
        print tajimas_theta(test_dict) / L


def testSub(iterations = 10):
    IN = siteRates =  np.loadtxt(open(mydir + \
    'Tree/data/rates/siteRate.txt','rb'),delimiter=',',skiprows=1)
    rates = siteRates[:,0]
    L = len(rates)
    ks = [ 1000]
    #Ns = [10, 100, 1000]
    Ns = [10]
    for k in ks:
        OUT = open(mydir + 'Combinatorics/data/SubRates/testSub_' + str(k) + '.txt', 'w+')
        for N in Ns:
            for i in range(iterations):

                test_dict = {}
                for x in range(N):
                    string = rd.subSeq(k)
                    if string not in  test_dict:
                        test_dict[string] = 1
                    else:
                        test_dict[string] += 1

                pi = tajimas_theta(test_dict) / L
                print  i,k, N, pi
                print>> OUT,i, k, N, pi
        OUT.close()


testSub()
