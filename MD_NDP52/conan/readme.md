---------------------------------------------
BAYESIAN ANALYSIS OF INTERMOLECULAR CONTACTS
---------------------------------------------
This procedure utilizes Python scripts to estimate the effective sample size (N_eff) from CONAN timeline data and subsequently calculate the posterior probability of specific intermolecular contacts. The Bayesian approach combines contact frequencies from two replicates, scoring and ranking interaction hubs based on region-specific structural thresholds.

## REQUIRED FILES:
Ensure the following files are available in the working directory:

- timeline1.dat --> CONAN timeline output file containing contact encounter history.
- rep1.csv --> Contact frequencies for replicate 1 (extracted from CONAN).
- rep2.csv --> Contact frequencies for replicate 2 (extracted from CONAN).
- calcu_neff.py --> Script to calculate the effective number of observations (N_eff).
- plot_conan_bayesian.py --> Main script for Bayesian calculation and plotting.

## SCRIPT EXECUTION

1. To estimate the Effective Sample Size (ESS) based on independent encounters:
```
module load pyhon
python calcu_neff.py`
```
2. Update the ESS_PER_REP variable in the plot_conan_bayesian.py script using the results from the previous step.

3. To run the Bayesian ranking and generate the plots:
`python plot_conan_bayesian.py`

## OUTPUT

The execution of the main script will generate the following files:
- Persistence_BayesianRank.png --> Dot plot of the top 20 actual contact frequencies ranked by their Bayesian reliability.
- Hub_Scores_Bayesian.png --> Dot plot of the normalized Bayesian Hub Scores for NDP52 LIR residues.
- Bayesian_Results.csv --> A dataset containing combined frequencies, alpha/beta parameters, structural thresholds, and posterior probabilities for all analyzed contacts.

## DATA
Original data available at:
```
/data/user/shared_projects/mavisp/CALCOCO2/simulations_analysis/complexes/lc3c/AF_125-147/replicate1/CHARMM36/model1/md/9.conan
/data/user/shared_projects/mavisp/CALCOCO2/simulations_analysis/complexes/lc3c/AF_125-147/replicate1/CHARMM36/model2/md/9.conan
```
