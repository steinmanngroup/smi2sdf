#!/usr/bin/env bash
#SBATCH --job-name sdf2mae
#SBATCH --partition kemi1
#SBATCH --time 48:00:00
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
export SCHRODINGER=/groups/kemi/cstein/programs/schrodinger2021-1

rm -f temp.sdf
for i in `seq -w 1 33708`
do
   if [ -e input/0${i}.sdf ]
   then
       cat input/0${i}.sdf >> temp.sdf
   fi
done

$SCHRODINGER/utilities/sdconvert -isd temp.sdf -omae rdkit.maegz
rm -f temp.sdf
