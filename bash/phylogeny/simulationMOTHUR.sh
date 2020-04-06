#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=72:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load python/2.7.13
module load gcc/6.3.0
module load boost/1.52.0
module load openmpi/gnu/2.1.0
module load mothur/1.39.5

## declare an array variable
declare -a arr=(10 100 1000)

## now loop through the above array
for i in "${arr[@]}"
do
    OUT="/N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Simulation/data/N${i}"
    mkdir $OUT -p
    for j in $(seq 100 $END);
    do
      python /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Simulation/Python/test.py \
          -n $i
      mothur /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Simulation/bin/simulationMOTHUR.batch
      mv /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Simulation/data/simulation_out.phylip.an.list \
          "${OUT}/simulationMOTHUR.phylip.an.list.${j}"
    done
    rm /N/dc2/projects/Lennon_Sequences/2016_Evol_16S/Simulation/bin/mothur.*.logfile
done
