#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=12:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc/4.9.2
module load boost/1.52.0
module load openmpi
module load mothur/1.38.1

# align seqs
cd /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/secondary_structure
# wget https://www.mothur.org/w/images/1/15/Silva.seed_v123.tgz
# tar -xvzf Silva.seed_v123.tgz
mothur "#align.seqs(candidate=structure.fasta, template=silva.seed_v123.align)"
