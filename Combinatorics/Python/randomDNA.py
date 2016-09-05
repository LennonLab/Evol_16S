from __future__ import division
from random import choice
import argparse


def String(length):

       DNA=""
       for count in range(length):
          DNA+=choice("CGTA")
       return DNA

def randomFASTA(n, path, length):
    OUT = open(path + '/randomDNA.fna', 'w+')
    for x in range(n):
        DNA = String(length)
        print>> OUT, '>' + str(x+1)
        print>> OUT, DNA
    OUT.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Run analyses and/or generate figures.")
    parser.add_argument('-f', default = '/N/dc2/projects/Lennon_Sequences/2016_Evol_16S/data/',\
        help = 'Whats the file path?')
    parser.add_argument('-n', type = int, default = 10, help = "How many sequences")
    parser.add_argument('-m', type = int, default = 500, help = "How many bases")

    args = parser.parse_args()
    path = args.f
    N = args.n
    M = args.m

    randomFASTA(N, path, M)
