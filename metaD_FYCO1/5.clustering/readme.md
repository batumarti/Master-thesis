---------------------------------------------
TRAJECTORY PROCESSING & CLUSTERING
---------------------------------------------
This procedure utilizes GROMACS tools to center and fit the trajectory, followed by an RMSD-based clustering using the GROMOS method. Finally, a Python script evaluates the thermodynamic population and relative Free Energy Surface (FES) of each cluster by reweighting the frames according to the previously calculated Metadynamics Boltzmann weights.

## REQUIRED FILES:
Ensure the following files are available in the working directory:

- traj_comp.xtc --> The raw/compressed trajectory file to analyze.
- sim.tpr --> GROMACS run input file containing the system topology.
- index.ndx --> GROMACS index file defining system groups (e.g., Protein, Complex, LC3B).
- COLVAR_WEIGHTED --> The colvar file with appended weights generated from step 2.
- weight_clusters.py --> Python script utilized for reweighting the clusters.

## EXECUTION

1. Make the protein whole to remove periodic boundary condition artifacts:
`gmx trjconv -f traj_comp.xtc -s sim.tpr -pbc whole -o traj_whole.xtc`
(Select the "Protein" group when prompted).

2. Center the trajectory with respect to a specific reference (e.g., LC3B) and compact the unit cell:
`gmx trjconv -s sim.tpr -f traj_whole.xtc -o traj_center.xtc -n index.ndx -center -pbc mol -ur compact`
(Select the reference group to center on, e.g., group 17, and then select the group to output, e.g., group 1 for the whole complex).

3. Fit rotation and translation to the reference structure:
`gmx trjconv -f traj_center.xtc -s sim.tpr -n index.ndx -fit rot+trans -o traj_fit.xtc`
(Select the "Protein" group when prompted).

4. Perform RMSD-based clustering using the GROMOS method:
`gmx cluster -f traj_fit.xtc -s sim.tpr -n index.ndx -method gromos -cutoff 0.3 -cl clusters.pdb -b 400000 -skip 500 -g cluster.log`

5. Execute the reweighting script to calculate relative FES and populations:
```
module load python
python weight_clusters.py -f COLVAR_WEIGHTED -c cluster.log -o clusters_reweighted.txt
```

## OUTPUT

The execution of the procedure will generate the following files:
- traj_whole.xtc --> Trajectory with whole molecules.
- traj_center.xtc --> Centered trajectory.
- traj_fit.xtc --> Final trajectory fitted to the reference structure.
- clusters.pdb --> PDB file containing the representative structures for each cluster.
- cluster.log --> Raw GROMACS log file detailing the frames belonging to each cluster.
- clusters_reweighted.txt --> A text file containing the sorted Cluster IDs, their reweighted Relative Free Energy (in kcal/mol), and their actual Thermodynamic Population (%).

## DATA
Original trajectory available at:
```
/data/raw_data/computational_data/simulations_data/lir_atg8/lc3b/lir_complexes/fyco1/fyco1_5WRDac_1-120_1270-1294/mod_maxg_219/CHARMM22star/metad/metad_pb/pb_2CVs
/data/raw_data/computational_data/simulations_data/lir_atg8/lc3b/lir_complexes/fyco1/fyco1_5WRDbd_1-120_1270-1294/mod_maxg_980/CHARMM22star/metad/metad_pb/pb_2CV_from_gefion/
```
