cd ../../..
filePCM=PreliminaryResultsAndConfigurations/pPb5TeV/directPhotons/data_PCMResultsFullCorrection_pPb_CentandMB_Run1_20180920.root
filePHOSCent=PreliminaryResultsAndConfigurations/pPb5TeV/directPhotons/PHOS_gamma_05092018_v2.root
doPHOSfix=1
filePHOSMB=PreliminaryResultsAndConfigurations/pPb5TeV/directPhotons/PHOS_pPb_final_19042018_directPhotonpPbMB_prelimInputQM18.root
fileEMC=PreliminaryResultsAndConfigurations/pPb5TeV/directPhotons/data_EMCResultsFullCorrection_pPb_Cent_20180926.root
filePCMEMC=PreliminaryResultsAndConfigurations/pPb5TeV/directPhotons/data_PCM-EMCResultsFullCorrection_pPb_Cent_20180926.root
fileCorrelations=PreliminaryResultsAndConfigurations/pPb5TeV/directPhotons/pPb5TeV_CF_20180926_MBandCent.root

root -x -l -b -q 'CombineGammaResultspPb_V2.C+("'$filePCM'","'$filePHOSCent'","'$filePHOSMB'","'$fileEMC'","'$filePCMEMC'","eps",kFALSE,kFALSE,'$doPHOSfix',"'$fileCorrelations'")'
