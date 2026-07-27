---------------------------------------------
BATCH FREE ENERGY SURFACE (FES) PIPELINE
---------------------------------------------
This procedure provides an automated pipeline to compute and plot 1D and 2D Free Energy Surfaces with statistical errors using block analysis. Calculations are performed in native PLUMED units (kJ/mol and nm) and are automatically converted to kcal/mol and Å during the plotting phase.

## REQUIRED FILES:
Ensure the following files are linked or available in the working directory:

- colvar --> The colvar file (e.g., COLVAR_WEIGHTED) containing the Collective Variables (CVs) and appended weights.
- config_fes.yaml --> Configuration file defining the CVs, block size, and grid bins.
- run_fes_pipeline.sh --> Main bash execution script.
- do_block_fes_v2.py --> Python script utilized for block analysis.
- plot_fes_1D.py --> Python script for 1D FES plotting.
- plot_fes_2D.py --> Python script for 2D FES plotting.

## EXECUTION

1. Edit the configuration file (config_fes.yaml) to set your specific tasks and general parameters (e.g., KBT, BSIZE).

2. Load the Python environment and execute the pipeline by running the main bash script:
'''
module load python
bash run_fes_pipeline.sh
'''
## CONFIGURATION FORMAT (config_fes.yaml)

Tasks are dynamically read from the TASKS array in the configuration file. The pipeline automatically detects whether to run a 1D or 2D analysis based on the string format:

- 1D Analysis format: "CV_NAME BINS" (e.g., "armsdC 80")
- 2D Analysis format: "CV1_NAME CV1_BINS CV2_NAME CV2_BINS" (e.g., "armsdC 80 d1 50")

## OUTPUT

The pipeline generates the following files for each requested task:

- fes_<bsize>_<cv(s)>_main.dat --> Raw data grid containing CVs, free energy, and statistical error.
- fes_<bsize>_<cv(s)>_walls.dat --> Raw data grid computed with the secondary walls bias (if PLOT_WALLS="true").
- plot_1D_<cv>.png --> Standard 1D plot with statistical error bands.
- plot_1D_<cv>_comp.png --> 1D comparison plot showing Main Bias vs +Walls Bias (no error bands).
- plot_2D_<cv1>_<cv2>.png --> 2D side-by-side contour maps showing the FES and its Statistical Error.
- plot_2D_<cv1>_<cv2>_walls.png --> 2D maps computed with the walls bias (if PLOT_WALLS="true").

---------------------------------------------
REFERENCE:
This script was developed following the PLUMED MASTERCLASS and the do_block_fes_v2.py script is a modified version of do_block_fes.py from M. Bonomi:
https://github.com/plumed/masterclass-21-4
