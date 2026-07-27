---------------------------------------------
BLOCK ANALYSIS FOR ERROR ESTIMATION
---------------------------------------------
This procedure utilizes scripts to perform block analysis on Metadynamics data, computing the statistical error as a function of block size to help select the optimal block dimension.

## REQUIRED FILES:
Ensure the following files are linked or available in the working directory:

- colvar_weighted --> The weighted colvar file generated from step 2.
- find_block_dim.sh --> Bash script to automate the block dimension screening.
- do_block_fes_v2.py --> Python script utilized for discrete binning and error computation.
- plot_block_err.py --> Python script to plot the resulting block error.

## EXECUTION

1. Inspect the "hills height vs. time" plot from the preliminary analysis to properly select the equilibration time.

2. Customize the parameters inside find_block_dim.sh if necessary. Default variables are:
```
        - MIN_BSIZE=10
        - NUM_POINTS=50
        - BIN_RESOLUTION=0.05
        - ORIGINAL_MASTER="colvar_weighted"
        - W_NAME="weight"  
        - KBT=2.477710
```
3. Run the bash script:
```
module load python
./find_block_dim.sh <CV_NAME> <EQUIL_TIME>
```
*Note: Replace <CV_NAME> with the exact name of the collective variable column in your colvar file (e.g., armsdC), and <EQUIL_TIME> with your chosen equilibration time.*

## OUTPUT

The execution of the procedure will compute the error as a function of block size, outputting data files and plots (via plot_block_err.py) that allow to evaluate convergence and select the appropriate block size for final error propagation.

---------------------------------------------
REFERENCE:
This script was developed following the PLUMED MASTERCLASS and the do_block_fes_v2.py script is a modified version of do_block_fes.py from M. Bonomi:
https://github.com/plumed/masterclass-21-4
