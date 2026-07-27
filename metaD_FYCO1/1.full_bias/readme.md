---------------------------------------------
TOTAL BIAS RECALCULATION - PLUMED DRIVER
---------------------------------------------
This procedure uses the PLUMED driver to read an existing trajectory and recalculate the exact potential applied to the system. The calculation includes both the Metadynamics bias (pb.bias) and the harmonic walls bias (walls.bias), yielding the total bias required for a correct Boltzmann reweighting.

## REQUIRED FILES:
Ensure the following files are linked or available in the working directory:

- traj_comp.xtc --> Trajectory to analyze.
- reference.pdb --> Reference topology/structure.
- plumed.dat --> PLUMED input file modified for passive reading (must include: RESTART, HEIGHT=0.0, PACE=100000000, walls grouped under a single label, and a customized colvar file name).
- HILLS_{CVs} --> Files containing the history of the Gaussians for all Collective Variables.

## EXECUTION

Run the following command to execute the PLUMED driver:
'plumed driver --mf_xtc traj_comp.xtc --plumed plumed.dat'

## OUTPUT

The execution of the driver will generate the following file:
- <customized_colvar_name> --> A colvar file where the last two columns contain `pb.bias` (the potential derived from Metadynamics in that specific frame) and `walls.bias` (the potential derived from the harmonic walls).
