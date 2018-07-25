Combintation of photon v2 results from PCM and PHOS and significance calculation

Macros:
read_pcm.C: read data from Mike and write output data_pcm.root in the form needed for combination macro
read_phos.C: read data from Dmitri and write output data_phos.root in the form needed for combination macro
v2dir_pcm_phos_comb.C: combination macro to calculate the direct photon v2 (PHOS, PCM, combined)
v2dir_significance.C: calculate significance of the measured v2inc w.r.t. to a null hypothesis in a Bayesian-inspired approach

Further macros (less important):
test_v2dir_significance_different_methods.C: compare different methods to calculate the significance
v2dir_vs_Rgamma.C: create plot v2dir vs. Rgamma as, e.g., shown in the paper

Both centralities in one go (output on terminal and also written to log file):
./run_v2dir_pcm_phos_comb.sh 2>&1 | tee logs/run_all_centralities.log

Significance calculation (log files created inside the script)
./run_v2dir_significance.sh

Significance for individual pt bins (log files created inside the script)
./run_v2dir_significance_individual_bins.sh

Extract significances from log file output:
./parse_significance_logfile.py <logfile.log>

Step-by-step (output directory must exist):
> root read_pcm.C
> root read_phos.C
> root
>> .x v2dir_pcm_phos_comb.C++
> root -b -q 'v2dir_pcm_phos_comb.C++("20-40", 10000), "output"'

Some macros for testing purposes:
test_fit_v2inc.C, test_fit_v2dec.C, test_fit_Rgamma.C --> test whether statistical (or more previsely, pt uncorrelated errors, make sense)



Notes:
- Results directory for the paper: output_rgam_corr_coeff_025

- Results directory for the paper with correct covariance matrix for cocktail v2: output_rgam_corr_coeff_025_2018_07_04
-- correct errors: significances: 0.9 sigma for 0-20% class, 0.5 sigma for 20-40% class)
-- artificially reduced statistical PHOS err (v2inc_stat_err * 0.25: significances: 0.6 sigma for 0-20% class, 0.8 sigma for 20-40% class)

- Test directory output_test: changes: corr. coeff. = 1 also for Rgamma, "new" method in v2dir_significance.C
-- correct errors: significances: 0.7 sigma for 0-20% class, 0.7 sigma for 20-40% class)
-- artificially reduced statistical PHOS err (v2inc_stat_err * 0.25: significances: 0.7 sigma for 0-20% class, 2.5 sigma for 20-40% class)


