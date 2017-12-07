#!/bin/bash
#PBS -k o
#PBS -l nodes=2:ppn=8,vmem=100gb,walltime=48:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load raxml/8.0.26

# delete second copy of n-16-b-T-thermophilus-2-num-bpseq before running this

cd /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/secondary_structure

raxmlHPC-PTHREADS -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -# autoMRE \
        -s ./structure_clean.align -n T1 -w /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/secondary_structure/tree/
