# Multivariate_resilience_indicators

Repositories with the codes to reproduces the analyses of the manuscript "Multivariate resilience indicators to anticipate vector-borne disease outbreaks: a West Nile virus case-study".

## Model
The model is implemented with the package SimInf, in R. It is provided in siminf4_extra.reservoir_wfeedingpref.R

## Perturbation recovery experiments
Perturbation recovery experiments are performed using the model in SimInf. They can be reproduced using perturbation_recovery_siminf2.R

## Resilience indicators analyses
The analyses to assess the performance of resilience indicators are performed in Matlab, using the generic_ews package.
They can be reproduced using:
- auc_results_new_scenario.m for the analyses of the performance of all indicators and all scenarios
- downsampling_new_scenarios.m for the analyses of the performance in data-poor scenarios, when reducing the sampling resolution
- obserror_new_scenarios.m for the analyses of the performance in data-poor scenarios, when decreasing the reporting probability


