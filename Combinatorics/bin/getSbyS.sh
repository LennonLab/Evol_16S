#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=0:30:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load python

python /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/Python/readList.py
