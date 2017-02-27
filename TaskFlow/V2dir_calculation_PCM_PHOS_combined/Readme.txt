Combintation of photon v2 results from PCM and PHOS

Macros:
read_pcm.C: read data from Mike and write output data_pcm.root in the form needed for combination macro
read_phos.C: read data from Dmitri and write output data_phos.root in the form needed for combination macro
v2dir_pcm_phos_comb.C: combination macro to calculate the direct photon v2 (PHOS, PCM, combined)

Step-by-step:
> root read_pcm.C
> root read_phos.C
> root
>> .x v2dir_pcm_phos_comb.C++
