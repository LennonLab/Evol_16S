#!/bin/bash
#PBS -k o
#PBS -l nodes=2:ppn=8,vmem=100gb,walltime=120:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load raxml/8.0.26

cd /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Evol_16S/data

#raxmlHPC -m GTRGAMMA -p 12345 -# 20 -s nmicrobiol201648-s7_clean.txt -n T13

#raxmlHPC-PTHREADS -T 4 -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE \
#  -s nmicrobiol201648-s7_clean.txt \
#  -n T14 \
#  -w /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Evol_16S/data/tree

#raxmlHPC-PTHREADS -T 4 -m PROTGAMMALG -J MRE -# autoMRE -f b \
#  -z /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Evol_16S/data/tree/RAxML_bootstrap.T14 \
#  -n T18 -w /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Evol_16S/data/tree


# Rapid Bootstrapping, does ML search + Bootstrapping in one step
#raxmlHPC-PTHREADS -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -# autoMRE \
#       -s /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Tree/data/alignment/nmicrobiol201648-s7_clean.txt -n T19 \
#       -w /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Tree/data/tree/
raxmlHPC-PTHREADS -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -# autoMRE \
        -s /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Tree/data/alignment/nmicrobiol201648-s7_clean_noEuks.txt -n T20 \
        -w /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Tree/data/tree/
