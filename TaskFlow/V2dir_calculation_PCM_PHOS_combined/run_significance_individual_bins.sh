echo "run_significance_individual_bins.sh: starting at `date`"

for i in {0..15}
do
    root -l -b -q 'v2dir_significance.C(0, "00-20",'"$i"','"$i"', 1000, "output_corr_coeff_1")' >> run_significance_individual_bins_00_20_corr_coeff_1.log
done

for i in {0..15}
do
    root -l -b -q 'v2dir_significance.C(0, "20-40",'"$i"','"$i"', 1000, "output_corr_coeff_1")' >> run_significance_individual_bins_20_40_corr_coeff_1.log
done

for i in {0..15}
do
    root -l -b -q 'v2dir_significance.C(0, "00-20",'"$i"','"$i"', 1000, "output_corr_coeff_05")' >> run_significance_individual_bins_00_20_corr_coeff_05.log
done

for i in {0..15}
do
    root -l -b -q 'v2dir_significance.C(0, "20-40",'"$i"','"$i"', 1000, "output_corr_coeff_05")' >> run_significance_individual_bins_20_40_corr_coeff_05.log
done

for i in {0..15}
do
    root -l -b -q 'v2dir_significance.C(0, "00-20",'"$i"','"$i"', 1000, "output_corr_coeff_0")' >> run_significance_individual_bins_00_20_corr_coeff_0.log
done

for i in {0..15}
do
    root -l -b -q 'v2dir_significance.C(0, "20-40",'"$i"','"$i"', 1000, "output_corr_coeff_0")' >> run_significance_individual_bins_20_40_corr_coeff_0.log
done
