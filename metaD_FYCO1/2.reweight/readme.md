---------------------------------------------
REWEIGHTING
---------------------------------------------
This procedure calculates weights considering both the Metadynamics bias and the walls bias by including the --walls flag.

## REQUIRED FILES:
Ensure the following files are linked or available in the working directory:

- colvar --> The colvar file generated from step 1.
- append_weights.py --> The Python script required to append the Boltzmann weights.

## EXECUTION

1. Execute the reweighting script:
module load python
python append_weights.py -f colvar --walls

## OUTPUT

The execution of the script will generate the following file:

        - colvar_weighted --> A new colvar file (default name) containing two additional columns:
          * weight --> calculated using only the primary bias (pb.bias).
          * weight_walls --> calculated using the sum of the primary bias and walls.bias (triggered by the --walls option).

---------------------------------------------
REFERENCE:
This script was developed following the PLUMED MASTERCLASS:
https://www.plumed.org/doc-v2.7/user-doc/html/masterclass-21-1.html
