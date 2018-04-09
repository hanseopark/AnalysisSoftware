Combintation of photon v2 results from PCM and PHOS

Macros:
read_pcm.C: read data from Mike and write output data_pcm.root in the form needed for combination macro
read_phos.C: read data from Dmitri and write output data_phos.root in the form needed for combination macro
v2dir_pcm_phos_comb.C: combination macro to calculate the direct photon v2 (PHOS, PCM, combined)
v2dir_significance.C: calculate significance of the measured v2inc w.r.t. to a null hypothesis in a Bayesian-inspired approach

Further macros (less important):
test_v2dir_significance_different_methods.C: compare different methods to calculate the significance
v2dir_vs_Rgamma.C: create plot v2dir vs. Rgamma as, e.g., shown in the paper

Both centralities in one go (output on terminal and also written to log file):
./run_all_centralities.sh 2>&1 | tee run_all_centralities.log

Significance calculation
./run_significance.sh 2>&1 | tee run_significance.log

Step-by-step (output directory must exist):
> root read_pcm.C
> root read_phos.C
> root
>> .x v2dir_pcm_phos_comb.C++
> root -b -q 'v2dir_pcm_phos_comb.C++("20-40", 10000), "output"'



