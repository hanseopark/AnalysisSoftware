echo "run_significance.sh: starting at `date`"

# one pt bin: around pT = 1.6 GeV/c
# root -l -b -q 'v2dir_significance.C(0, "00-20", 3, 3, 1000, "output")'

# pt range 0.9 - 2.1 GeV/c
# root -l -b -q 'v2dir_significance.C(0, "00-20", 0, 5, 1000, "output")'

# pt range 0.9 - 3.0 GeV/c
# root -l -b -q 'v2dir_significance.C(0, "00-20", 0, 9, 1000, "output")'

# full pt range
# root -l -b -q 'v2dir_significance.C(0, "00-20", 0, 15, 1000, "output")'

# one pt bin: around pT = 1.6 GeV/c
# root -l -b -q 'v2dir_significance.C(0, "20-40", 3, 3, 1000, "output")'

# pt range 0.9 - 2.1 GeV/c
# root -l -b -q 'v2dir_significance.C(0, "20-40", 0, 5, 1000, "output")'

# pt range 0.9 - 3.0 GeV/c
# root -l -b -q 'v2dir_significance.C(0, "20-40", 0, 9, 1000, "output")'

# full pt range
# root -l -b -q 'v2dir_significance.C(0, "20-40", 0, 15, 1000, "output")'

root -l -b -q 'v2dir_significance.C(0, "00-20", 0, 5, 10000, "output_corr_coeff_1")'
root -l -b -q 'v2dir_significance.C(0, "20-40", 0, 5, 10000, "output_corr_coeff_1")'

root -l -b -q 'v2dir_significance.C(0, "00-20", 0, 5, 10000, "output_corr_coeff_05")'
root -l -b -q 'v2dir_significance.C(0, "20-40", 0, 5, 10000, "output_corr_coeff_05")'

root -l -b -q 'v2dir_significance.C(0, "00-20", 0, 5, 10000, "output_corr_coeff_0")'
root -l -b -q 'v2dir_significance.C(0, "20-40", 0, 5, 10000, "output_corr_coeff_0")'


