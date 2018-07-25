echo "run_v2dir_significance.sh: starting at `date`"
# root -l -b -q 'v2dir_significance.C(0, "00-20", 0, 5, 10000, "output_test")' > output_test/v2dir_significance_00_20.log
# root -l -b -q 'v2dir_significance.C(0, "20-40", 0, 5, 10000, "output_test")' > output_test/v2dir_significance_20_40.log
root -l -b -q 'v2dir_significance.C(0, "00-20", 0, 5, 10000, "output_projection_for_yellow_report")' > output_projection_for_yellow_report/v2dir_significance_00_20.log
root -l -b -q 'v2dir_significance.C(0, "20-40", 0, 5, 10000, "output_projection_for_yellow_report")' > output_projection_for_yellow_report/v2dir_significance_20_40.log
