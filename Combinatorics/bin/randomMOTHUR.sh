#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=72:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc/4.9.2
module load boost/1.52.0
module load openmpi
module load mothur/1.38.1

## declare an array variable
declare -a arr=(10 100 1000 10000 100000)
M=576
mkdir "/N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/data/M${M}" -p
#declare -a arr=(10 100)
## now loop through the above array
for i in "${arr[@]}"
do
    OUT="/N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/data/M${M}/N${i}"
    mkdir $OUT -p
    for j in $(seq 1000 $END);
    do
      python /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/Python/randomDNA.py \
          -n $i -m $M -f /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/data
      mothur /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/bin/randomMOTHUR.batch
      mv /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/data/randomDNA.phylip.an.list \
          "${OUT}/randomDNA.phylip.an.list.${j}"
    done
    rm /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/bin/mothur.*.logfile
done
