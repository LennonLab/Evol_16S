#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=24:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load python
module load openmpi
module load mothur

## declare an array variable
declare -a arr=(10 100 1000 10000 100000)
#declare -a arr=(10)
## now loop through the above array
for i in "${arr[@]}"
do
    OUT="/N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/data/N${i}"
    mkdir $OUT -p
    for j in $(seq 1000 $END);
    do
      python /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/Python/randomDNA.py \
          -n $i -m 10 -f /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/data
      mothur /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/bin/randomMOTHUR.batch
      mv /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/data/randomDNA.phylip.an.list \
          "${OUT}/randomDNA.phylip.an.list.${j}"
    done
    rm /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Combinatorics/bin/mothur.*.logfile
done
