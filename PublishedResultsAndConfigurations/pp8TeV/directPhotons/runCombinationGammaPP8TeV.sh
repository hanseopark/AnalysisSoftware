cd ../../..
root -x -l -b -q 'CombineGammaResultsPP8TeV.C+("ExternalInput/PCM/directPhotons/data_PCMResultsFullCorrection_PP8TeV_BinShifted_20180724.root", kFALSE, "", kTRUE, "ExternalInput/EMCAL/directPhotons/data_EMCResultsFullCorrection_PP8TeV_BinShifted_20180724.root", kTRUE, "ExternalInput/PCMEMC/directPhotons/data_PCM-EMCResultsFullCorrection_PP8TeV_BinShifted_20180724.root", "ExternalInput/CombDirectPhotons/GammaCocktailRatios_PP2760GeV_PP8TeV_20180724.root", "eps", "ExternalInput/CombDirectPhotons/dirGammaCorrelationInput_PP8TeV_20180724.root",kFALSE,1.28,kTRUE,"")'