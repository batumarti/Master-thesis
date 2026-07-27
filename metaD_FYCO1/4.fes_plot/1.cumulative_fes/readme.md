---------------------------------------------
CUMULATIVE FES PLOTTING
---------------------------------------------
This procedure plots cumulative 1D Free Energy Surfaces over time to evaluate the convergence of the simulation. It calculates histograms using native PLUMED units (kJ/mol and nm) and converts them to kcal/mol and Å for plotting.

## REQUIRED FILES:
Ensure the following file is linked or available in the working directory:

- COLVAR_WEIGHTED --> The colvar file with appended weights generated from step 2.

## EXECUTION

Load the Python environment and execute the script:
'''
module load python
python plot_cumulative_fes.py -f COLVAR_WEIGHTED --cv armsdC --min 0.0 --max 1.5 --bins 100 --nm-to-a
'''
*Additional useful flags:*
        --wcol      --> Specify the weight column (default is "weight", use "weight_walls" if needed).
        --eq-time   --> Time to discard as equilibration in ps (default: 400000).
        --step-time --> Time step for cumulative blocks in ps (default: 200000).
        -o          --> Specify the output file name (default: cumulative_fes_<cv>.png).

## OUTPUT

The execution of the script will generate the following file:
- cumulative_fes_<cv>.png --> A plot showing the FES profile at different time slices.
