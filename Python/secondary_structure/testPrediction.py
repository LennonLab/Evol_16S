from __future__ import division
import os
import pandas as pd
import RNA
import sys
sys.path.append("/usr/local/ViennaRNA")
import RNA

mydir = os.path.expanduser('~/github/Evol_16S/')

energy = { ('A', 'U'): -4.6,  ('G', 'C'): -6.9, \
        ('U', 'A'): -4.6, ('C', 'G'): -6.9}


def calculateFE():
    path = mydir + 'data/secondary_structure/structure.bpseq/d.16.b.B.subtilis.bpseq'
    IN = pd.read_csv(path, skiprows=4, header=None, sep = ' ', names = ["position", "nuc", "pair"])
    print(IN)


calculateFE()

#./configure --enable-universal-binary --disable-openmp --without-perl
