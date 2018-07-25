echo "run_v2dir_significance_individual_bins.sh: starting at `date`"

#for i in {0..15}
#do
#    root -l -b -q 'v2dir_significance.C(0, "00-20",'"$i"','"$i"', 1000, "output_rgam_corr_coeff_025")' >> run_significance_individual_bins_00_20_output_rgam_corr_coeff_025.log
#done

#for i in {0..15}
#do
#    root -l -b -q 'v2dir_significance.C(0, "20-40",'"$i"','"$i"', 1000, "output_rgam_corr_coeff_025")' >> run_significance_individual_bins_20_40_output_rgam_corr_coeff_025.log
#done


for i in {0..15}
do
    root -l -b -q 'v2dir_significance.C(0, "00-20",'"$i"','"$i"', 1000, "output_rgam_corr_coeff_025_2018_07_04")' >> logs/run_significance_individual_bins_00_20_output_rgam_corr_coeff_025_2018_07_04.log
    echo "Hallo"
done

for i in {0..15}
do
    root -l -b -q 'v2dir_significance.C(0, "20-40",'"$i"','"$i"', 1000, "output_rgam_corr_coeff_025_2018_07_04")' >> logs/run_significance_individual_bins_20_40_output_rgam_corr_coeff_025_2018_07_04.log
done
