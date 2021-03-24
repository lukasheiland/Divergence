# Divergent occurrences of juvenile and adult trees are explained by both environmental change and ontogenetic niche shifts

This repository provides code and aggregated data to reproduce the simulations and analyses in "Divergent occurrences of juvenile and adult trees are explained by both environmental change and ontogenetic niche shifts". Session info is provided in 'DEF_sessionInfo__2020-11-02.*'.

The simulations are contained in
- 'Sim ontogenetic stages.R'

The main script for the empirical analysis
- 'Fit beta.R' performs the complete analysis with already aggregated data that is provided with 'Data/taxtables_pres_thresholdsubset.rds'.

In addition, several other scripts are provided for reference, which can be explored interactively, but would need additional unanomyzed data to be run completely:
- 'Prepare data.R' had been run for wrangling the unanonymized data prior to the analysis.
- 'Publishing/Publish fit.R' creates maps, plots, and tables. (For maps, original coordinates would be necessary.)
