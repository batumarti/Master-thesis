#!/bin/bash

#TRAJECTORY PROCESSING
#make the protein whole
echo "Protein" | gmx trjconv -f traj.xtc -s sim.tpr -pbc whole -o traj_whole.xtc

#center the protein in a compact box
echo "Protein Protein" | gmx trjconv -f traj_whole.xtc -s sim.tpr -n index.ndx -pbc cluster -center -ur compact -o traj_center.xtc

#fit rotation and translation to the reference structure
echo "Protein Protein" | gmx trjconv -f traj_center.xtc -s sim.tpr -n index.ndx -fit rot+trans -o traj_fit.xtc

rm traj_whole.xtc traj_center.xtc

#CLUSTERING
echo "Protein Protein" | gmx cluster -f traj_fit.xtc -s sim.tpr -n index.ndx -method gromos -cutoff 0.3 -cl clusters.pdb -skip 500 -g cluster.log
