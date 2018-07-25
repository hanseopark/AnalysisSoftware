echo "run_all_centralities.sh: starting at `date`"
root -b -q read_pcm.C
root -b -q read_phos.C
root -b -q read_cocktail.C
# root -b -q 'v2dir_pcm_phos_comb.C("00-20", 10000, "output")'
# root -b -q 'v2dir_pcm_phos_comb_test.C("00-20", 100000, "output")'
# root -b -q 'v2dir_pcm_phos_comb_test.C("20-40", 100, "output")'

# root -b -q 'v2dir_pcm_phos_comb.C("00-20", 100000, "output_rgam_corr_coeff_025_2018_07_04")'
# root -b -q 'v2dir_pcm_phos_comb.C("20-40", 100000, "output_rgam_corr_coeff_025_2018_07_04")'
# root -b -q 'v2dir_pcm_phos_comb.C("00-20", 10000, "output_rgam_corr_coeff_025_2018_07_04")'
# root -b -q 'v2dir_pcm_phos_comb.C("20-40", 10000, "output_rgam_corr_coeff_025_2018_07_04")'
# root -b -q 'v2dir_pcm_phos_comb.C("00-20", 10000, "output_test")' 
# root -b -q 'v2dir_pcm_phos_comb.C("20-40", 10000, "output_test")'
# root -b -q 'v2dir_pcm_phos_comb.C("00-20", 10000, "output_rgam_corr_coeff_025_2018_07_20")'
# root -b -q 'v2dir_pcm_phos_comb.C("20-40", 10000, "output_rgam_corr_coeff_025_2018_07_20")'
root -b -q 'v2dir_pcm_phos_comb.C("00-20", 10000, "output_projection_for_yellow_report")' 
root -b -q 'v2dir_pcm_phos_comb.C("20-40", 10000, "output_projection_for_yellow_report")'
