root -b -q read_pcm.C
root -b -q read_phos.C
root -b -q read_cocktail.C
root -b -q 'v2dir_pcm_phos_comb.C++("00-20", 10000, "output")'
root -b -q 'v2dir_pcm_phos_comb.C++("20-40", 10000, "output")'
# root -b -q 'v2dir_pcm_phos_comb.C++("00-20", 100000, "output")'
# root -b -q 'v2dir_pcm_phos_comb.C++("20-40", 100000, "output")'
