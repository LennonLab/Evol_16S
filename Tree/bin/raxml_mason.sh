#!/bin/bash
#PBS -k o
#PBS -l nodes=2:ppn=8,vmem=100gb,walltime=96:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load raxml/8.0.26

cd /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Evol_16S/data

raxmlHPC-PTHREADS -T 4 -m GTRCAT -p 12345 -x 12345 -# 300 \
  -s nmicrobiol201648-s7_clean.txt \
  -n T14 \
  -w /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Evol_16S/data/tree

#raxmlHPC-PTHREADS -T 4 -m GTRCAT -J MRE -# 100 \
#  -z /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Evol_16S/data/tree/RAxML_bootstrap.T14 \
#  -n T18 -w /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Evol_16S/data/tree
