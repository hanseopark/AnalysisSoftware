#void ExtractSignalPiPlPiMiPiZero(TString meson="", TString subDir="", TString file="", TString cutSelection="", TString Suffix="", TString optionMC="", TString optionEnergy="", TString optionCrystalBall="", TString directphotonPlots="", TString optionUseMinBiasEff="", TString optionPeriod="", TString optionAdvancedMesonQA="",Int_t numberOfBins = 30, Bool_t addSig = kFALSE, Bool_t makeQAplots = kTRUE)


root -l -b -q -x 'TaskV1/ExtractSignalPiPlPiMiPiZero.C++("Omega","/home/jens/Cloud/Sciebo/Linux_Arbeitsbereich/Data/7TeV_trainRun/LHC14j4_GammaConvNeutralMesonPiPlPiMiPiZero_2_22.root","2_00000113_1111111067032230000_0103503400000000_002010706_0103503000000000","eps","kFALSE","7TeV","Gaussian","","","","",8,kFALSE,kTRUE)'


