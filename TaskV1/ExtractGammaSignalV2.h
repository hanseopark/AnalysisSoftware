// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

#ifndef GAMMACONV_ExtractGammaSignal
#define GAMMACONV_ExtractGammaSignal

    //!! -> what is commented here is defined in ExtractSignal.h
    Bool_t      fEnablePCM                                              = 0;
    Bool_t      fEnableCalo                                             = 0;
    Bool_t      fUseCocktail                                            = 0;
    Bool_t      fUseDataDrivenPurity                                    = 0;

    TString     fOutputDir                                              = "";
    TString     fSuffix                                                 = "";
    TString     fMeson                                                  = "";

    TString     fCutSelectionRead                                       = "";
    TString     fEventCutSelectionRead                                  = "";
    TString     fGammaCutSelectionRead                                  = "";
    TString     fMesonCutSelectionRead                                  = "";
    TString     fClusterCutSelectionRead                                = "";

    TString     fHistogramDimension                                     = "";
    Int_t       nHistogramDimension                                     = 1;
    Int_t       nPileupMethodUsed                                       = 0;

    Int_t       fNBinsPtDummy                                           = 0;
    Double_t*   fBinsPtDummy                                            = NULL;

    Bool_t      fEnableDCConv                                           = kFALSE;

    TString     fDecays[7]                                              = { "Pi0",  "Eta",   "Etap", "Omega",   "Rho",
                                                                            "Phi",  "Sigma"
                                                                        };
    Int_t       nCombinatorics                                          = 16;
    TString     fCombinatorics[17]                                      = { "Elec+Elec",     "Elec+Pion",   "Elec+Kaon",        "Elec+Proton",   "Elec+Muon",
                                                                            "Pion+Pion",     "Pion+Kaon",   "Pion+Proton",      "Pion+Muon",     "Kaon+Kaon",
                                                                            "Kaon+Proton",   "Kaon+Muon",   "Proton+Proton",    "Proton+Muon",   "Muon+Muon",
                                                                            "Rest",          "All"
                                                                        };
    Int_t       nContamination                                          = 10;
    TString     fContamination[11]                                      = { "Electron", "Pion", "Proton",   "Kaon", "Neutron", "K0s",
                                                                            "Lambda",   "Muon", "K0l",      "Rest", "All"
                                                                        };

    TString     fSecondaries[4]                                         = {"K0s", "K0l", "Lambda", "Rest"};

    // exemplary pt bins for DCAz distributions
    Int_t       exemplaryLowPtBin                                       = 1;
    Int_t       exemplaryHighPtBin                                      = 15;

    // ShowBackground arguments
    Int_t       nIterationsShowBackground[4]                            = {0};
    TString     optionShowBackground[5]                                 = {""};
    Int_t       fNOOBEstMethods                                         = 5;
    Int_t       fNOOBCat                                                = 3;

    // binning
    TH1D*       fDeltaPtDummy                                           = NULL;

    //*****************************************************************************************************
    //******************** Histograms for conversion analysis *********************************************
    //*****************************************************************************************************
    TH1D*       fHistoPurityKappaTemplates                                                          = NULL;

    TH1D*       fHistoGammaConvPt                                                                   = NULL;
    TH1D*       fHistoGammaConvPtOrBin                                                              = NULL;
    TH2D*       f2DHistoSecondaryGammaMCConvPt                                                      = NULL;
    TH2D*       f2DHistoAllSecondaryGammaMCPt                                                       = NULL;
    TH1D*       fHistoGammaMCrecConvPt                                                              = NULL;
    TH1D*       fHistoGammaMCrecConvPtOrBin                                                         = NULL;
    TH1D*       fHistoGammaMCConvPt                                                                 = NULL;
    TH1D*       fHistoGammaMCConvPtOrBin                                                            = NULL;
    TH1D*       fHistoGammaMCConvRSPt                                                               = NULL;
    TH1D*       fHistoGammaMCAllPt                                                                  = NULL;
    TH1D*       fHistoGammaMCAllPtOrBin                                                             = NULL;

    TH1D*       fHistoSecondaryGammaConvFromXPt[4]                                                  = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoSecondaryGammaConvFromXPtOrBin[4]                                             = { NULL, NULL, NULL, NULL };

    TH1D*       fHistoSecondaryGammaCocktailFromXPt[4]                                              = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoSecondaryGammaCocktailFromXPtOrBin[4]                                         = { NULL, NULL, NULL, NULL };

    TH1D*       fHistoAllSecondaryGammaFromXPt[4]                                                   = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoAllSecondaryGammaFromXPtOrBin[4]                                              = { NULL, NULL, NULL, NULL };

    TH1D*       fHistoGammaTrueConvPt                                                               = NULL;
    TH1D*       fHistoGammaTrueConvPtOrBin                                                          = NULL;
    TH1D*       fHistoGammaTruePrimaryConvPt                                                        = NULL;
    TH1D*       fHistoGammaTruePrimaryConvPtOrBin                                                   = NULL;
    TH1D*       fHistoGammaTrueSecondaryConvPt                                                      = NULL;
    TH2D*       f2DHistoGammaTrueSecondaryConvPt                                                    = NULL;
    TH2D*       f2DHistoGammaTrueSecondaryConvMCPt                                                  = NULL;
    TH1D*       fHistoGammaTrueSecondaryConvPtOrBin                                                 = NULL;
    TH1D*       fHistoGammaTrueSecondaryConvGammaFromXMCPt[4]                                       = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaTrueSecondaryConvGammaFromXMCPtOrBin[4]                                  = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaTrueSecondaryConvGammaFromXPt[4]                                         = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaTrueSecondaryConvGammaFromXPtOrBin[4]                                    = { NULL, NULL, NULL, NULL };
    TH2D*       fHistoGammaTrueSecondaryFromXConv_MCPt_recPt_MC[3]                                  = { NULL, NULL, NULL };

    TH1D**      fHistoGammaMCDecayPt                                                                = NULL;
    TH1D*       fHistoGammaTruePrimaryConvMCPt                                                      = NULL;
    TH1D*       fHistoGammaTruePrimaryConvMCPtOrBin                                                 = NULL;

    TH1D*       fHistoFracAllGammaToSecOrBin                                                        = NULL;
    TH1D*       fHistoFracAllGammaToSec                                                             = NULL;
    TH1D*       fHistoFracAllGammaToSecFromX[4]                                                     = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoFracAllGammaToSecFromXOrBin[4]                                                = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaMCConvProb                                                               = NULL;
    TH1D*       fHistoGammaMCConvProbOrBin                                                          = NULL;
    TH1D*       fHistoSecondaryGammaFromXMCConvProb[3]                                              = { NULL, NULL, NULL };
    TH1D*       fHistoSecondaryGammaFromXMCConvProbOrBin[3]                                         = { NULL, NULL, NULL };
    TH1D*       fHistoGammaMCPurity                                                                 = NULL;
    TH1D*       fHistoGammaMCTruePurityOrBin                                                        = NULL;
    TH1D*       fHistoGammaMCrecPrimaryConvPt                                                       = NULL;
    TH1D*       fHistoGammaMCTruePurity                                                             = NULL;
    TH1D*       fHistoGammaMCrecPrimaryConvPtOrBin                                                  = NULL;
    TH1D*       fHistoGammaMCRecoEff                                                                = NULL;
    TH1D*       fHistoGammaMCPrimaryRecoEff                                                         = NULL;
    TH1D*       fHistoGammaMCPrimaryRecoEffOrBin                                                    = NULL;
    TH1D*       fHistoGammaMCPrimaryRecoEffMCPt                                                     = NULL;
    TH1D*       fHistoGammaMCPrimaryRecoEffMCPtOrBin                                                = NULL;
    TH1D*       fHistoGammaMCBackground                                                             = NULL;

    TH1D*       fHistoSecondaryGammaFromXMCRecoEffMCPt[3]                                           = { NULL, NULL, NULL };
    TH1D*       fHistoSecondaryGammaFromXMCRecoEffMCPtOrBin[3]                                      = { NULL, NULL, NULL };

    TH1D*       fHistoSecondaryGammaFromXMCRecoEffPt[3]                                             = { NULL, NULL, NULL };
    TH1D*       fHistoSecondaryGammaFromXMCRecoEffPtOrBin[3]                                        = { NULL, NULL, NULL };

    TH2D*       fHistoGammaTrueSecondaryFromXConv_MCPt_recPt_MC_Rebin[3]                            = { NULL, NULL, NULL };

    TH2D*       fHistoCombinatorialBackground                                                       = NULL;
    TH2D*       fHistoGammaTruePrimaryConv_recPt_MCPt_MC                                            = NULL;
    TH2D*       fHistoGammaTruePrimaryConv_recPt_MCPt_MC_Rebin                                      = NULL;

    TH2F*       fMCrecGammaPtDCAz                                                                   = NULL;
    TH2F**      fESDGammaPtDCAz                                                                     = NULL;
    TH1D***     fESDGammaPtDCAzBins                                                                 = NULL;
    TH1D****    fESDGammaPtDCAzBinsBack                                                             = NULL;
    TH1D****    fESDSubGammaPtDCAzBins                                                              = NULL;

    TH1D***     fMCrecGammaPtDCAzBins                                                               = NULL;
    TH1D***     fMCrecGammaPtDCAzBinsBack                                                           = NULL;
    TH1D***     fMCrecSubGammaPtDCAzBins                                                            = NULL;

    TH2F**      fTruePrimaryPhotonPtDCAz                                                            = NULL;
    TH1D***     fTruePrimaryGammaPtDCAzBins                                                         = NULL;
    TH1D***     fTruePrimarySubGammaPtDCAzBins                                                      = NULL;
    TH2F**      fTrueSecondaryPhotonPtDCAz                                                          = NULL;
    TH1D***     fTrueSecondaryGammaPtDCAzBins                                                       = NULL;
    TH1D***     fTrueSecondarySubGammaPtDCAzBins                                                    = NULL;
    TH2F**      fTrueSecondaryPhotonFromXPtDCAz[4]                                                  = {NULL, NULL, NULL, NULL};
    TH1D***     fTrueSecondaryGammaFromXPtDCAzBins[4]                                               = {NULL, NULL, NULL, NULL};
    TH1D***     fTrueSecondarySubGammaFromXPtDCAzBins[4]                                            = {NULL, NULL, NULL, NULL};
    TH1D***     fTrueBackgroundPtDCAzBins                                                           = NULL;
    TH1D***     fTrueGammaPtDCAzBins                                                                = NULL;
    TH1D***     fTrueSubGammaPtDCAzBins                                                             = NULL;
    TH1D**      fHistoCombinatorialSpecies                                                          = NULL;

    TH1D**      fESDGammaPerCatPtDCAzBins                                                           = NULL;
    TH1D**      fESDGammaRatioCatToCombinedPtDCAzBins                                               = NULL;

    TH1D*       fESDGammaPtPileUpAllCat                                                             = NULL;
    TH1D*       fESDGammaPtPileUp                                                                   = NULL;
    TH1D**      fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat                              = NULL;
    TH1D**      fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning                                    = NULL;
    TF1**       fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat                           = NULL;
    TF1**       fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning                                 = NULL;

    TH1D*       fMCrecGammaPtPileUpAllCat                                                           = NULL;
    TH1D*       fMCrecGammaPtPileUp                                                                 = NULL;
    TH1D*       fTruePrimaryConvGammaPtPileUpAllCat                                                 = NULL;
    TH1D*       fTruePrimaryConvGammaPtPileUp                                                       = NULL;
    TH1D*       fTrueSecondaryConvGammaPtPileUpAllCat                                               = NULL;
    TH1D*       fTrueSecondaryConvGammaPtPileUp                                                     = NULL;
    TH1D*       fHistoGammaTrueSecondaryConvGammaFromXPtPileUpAllCat[4]                             = {NULL, NULL, NULL, NULL};
    TH1D*       fHistoGammaTrueSecondaryConvGammaFromXPtPileUp[4]                                   = {NULL, NULL, NULL, NULL};

    TH1D*       fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat                            = NULL;
    TH1D*       fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning                                  = NULL;
    TF1*        fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat                         = NULL;
    TF1*        fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinning                               = NULL;

    TH1D*       fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat                  = NULL;
    TH1D*       fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning                        = NULL;
    TF1*        fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat               = NULL;
    TF1*        fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning                     = NULL;

    TH1D*       fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat                = NULL;
    TH1D*       fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning                      = NULL;
    TF1*        fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat             = NULL;
    TF1*        fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning                   = NULL;

    TH1D*       fTrueSecondaryFromXConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[4]        = {NULL, NULL, NULL, NULL};
    TH1D*       fTrueSecondaryFromXConvGammaPtRatioWithWithoutPileUpDCAzDistBinning[4]              = {NULL, NULL, NULL, NULL};
    TF1*        fTrueSecondaryFromXConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat[4]     = {NULL, NULL, NULL, NULL};
    TF1*        fTrueSecondaryFromXConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[4]           = {NULL, NULL, NULL, NULL};

    TH1D**      fESDGammaPileUpCorrFactorAllCat                                                     = NULL;
    TH1D**      fESDGammaPileUpCorrFactor                                                           = NULL;
    TGraphAsymmErrors* sysErrOOBPileupDown                                                          = NULL;
    TGraphAsymmErrors* sysErrOOBPileupUp                                                            = NULL;

    TH1D*       fMCrecGammaPileUpCorrFactorAllCat                                                   = NULL;
    TH1D*       fMCrecGammaPileUpCorrFactor                                                         = NULL;
    TH1D*       fTruePrimaryConvGammaPileUpCorrFactorAllCat                                         = NULL;
    TH1D*       fTruePrimaryConvGammaPileUpCorrFactor                                               = NULL;
    TH1D*       fTrueSecondaryConvGammaPileUpCorrFactorAllCat                                       = NULL;
    TH1D*       fTrueSecondaryConvGammaPileUpCorrFactor                                             = NULL;
    TH1D*       fTrueSecondaryFromXConvGammaPileUpCorrFactorAllCat[4]                               = {NULL, NULL, NULL, NULL};
    TH1D*       fTrueSecondaryFromXConvGammaPileUpCorrFactor[4]                                     = {NULL, NULL, NULL, NULL};

    TH1D*       fHistoFracAllGammaToSecPileUp                                                       = NULL;
    TH1D*       fHistoFracAllGammaToSecFromXPileUp[4]                                               = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaMCPurityPileUp                                                           = NULL;
    TH1D*       fHistoGammaMCrecPrimaryConvPtPileUp                                                 = NULL;
    TH1D*       fHistoGammaMCTruePurityPileUp                                                       = NULL;
    TH1D*       fHistoGammaMCRecoEffPileUp                                                          = NULL;
    TH1D*       fHistoGammaMCPrimaryRecoEffPileUp                                                   = NULL;

    TH1D*       fHistoPhotonIsSelected                                                              = NULL;

    //******************* tagging outputHistograms ********************************************************
    TH1D**      fHistoGGInvMassPtGConvBin                                                           = NULL;
    TH1D**      fHistoSignalInvMassPtGConvBin                                                       = NULL;
    TH1D**      fHistoSignalInvMassLeftPtGConvBin                                                   = NULL;
    TH1D**      fHistoBackInvMassPtGconvBin                                                         = NULL;
    TH1D**      fHistoBackNormInvMassPtGconvBin                                                     = NULL;
    TH1D**      fHistoBackNormInvMassLeftPtGconvBin                                                 = NULL;
    TH1D**      fHistoTruePrimMesonInvMassPtBins                                                    = NULL;
    TH1D**      fHistoTrueFullMesonInvMassPtBins                                                    = NULL;
    TH1D**      fHistoTrueSecMesonInvMassPtBins[4]                                                  = {NULL, NULL, NULL, NULL};
    TH2D*       fHistoTruePrimMesonInvMassVSPt                                                      = NULL;
    TF1**       fFitTrueFullSignalInvMassPtBin                                                      = NULL;
    TH2D*       fHistoMCPrimPi0PtGammaLeg                                                           = NULL;
    TH2D*       fHistoMCPrimPi0InAccPtGammaLeg                                                      = NULL;
    TH2D*       fHistoMCSecPi0PtGammaLeg[2]                                                         = {NULL, NULL };
    TH2D*       fHistoMCSecPi0InAccPtGammaLeg[2]                                                    = {NULL, NULL };
    TH1D*       fHistoMCPrimPi0PtGamma[3]                                                           = {NULL, NULL, NULL };
    TH1D*       fHistoMCPrimPi0InAccPtGamma[3]                                                      = {NULL, NULL, NULL };
    TH1D*       fHistoMCSecPi0PtGamma[4][3]                                                         = {{NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }};
    TH1D*       fHistoMCSecPi0InAccPtGamma[4][3]                                                    = {{NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }};
    TH1D*       fHistoMCPrimPi0PtGammaReb[3]                                                        = {NULL, NULL, NULL };
    TH1D*       fHistoMCPrimPi0InAccPtGammaReb[3]                                                   = {NULL, NULL, NULL };
    TH1D*       fHistoMCSecPi0PtGammaReb[4][3]                                                      = {{NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }};
    TH1D*       fHistoMCSecPi0InAccPtGammaReb[4][3]                                                 = {{NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }};
    TH1D*       fHistoAcceptancePi0PtGammaReb[3]                                                    = {NULL, NULL, NULL };
    TH1D*       fHistoAcceptanceSecPi0PtGammaReb[4][3]                                              = {{NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }};
    TH1D*       fHistoTruePi0EffiGammaPt[3]                                                         = {NULL, NULL, NULL };
    TH1D*       fHistoTrueSecPi0EffiGammaPt[4][3]                                                   = {{NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }, {NULL, NULL, NULL }};

    //*****************************************************************************************************
    //*********************** Pi0Tagging histograms *******************************************************
    //*****************************************************************************************************

    TH1D*       fHistoTruePrimaryPi0MissingPtGconv                                                  = NULL;
    TH1D*       fHistoTruePrimaryPi0MissingPtGconv_Rebin                                            = NULL;
    TH1D*       fHistoTruePrimaryPi0Sum                                                             = NULL;
    TH1D*       fHistoTruePi0TaggingEfficiency                                                      = NULL;
    TH1D*       fHistoPi0TaggingEfficiency                                                          = NULL;
    TH2D*       fHistoTrueSecondaryMesonMissingPtGconv[4]                                           = {NULL, NULL, NULL, NULL};

    //*****************************************************************************************************
    //******************** Check multiple counts of photons ***********************************************
    //*****************************************************************************************************
    TH2D*       fHistoTrueGammaConvDCRVSPt                                                          = NULL;
    TH1D*       fHistoTrueGammaConvDCR                                                              = NULL;
    TH1D*       fHistoTrueGammaConvDCPt                                                             = NULL;
    TH1F*       fHistoTrueGammaConvMultipleCount                                                    = NULL;
    TString     ObjectNameDCGammaConvRPt                                                            = "";
    TString     ObjectNameGammaConvMultipleCount                                                    = "";

    TTree*      dcaTree                                                                             = NULL;
    TList*      DCAContainer                                                                        = NULL;
    Bool_t      pileUpCorrection                                                                    = kFALSE;

    //*****************************************************************************************************
    //******************** Histograms for calo analysis ***************************************************
    //*****************************************************************************************************
    TList*      CaloContainer                                                                       = NULL;
    TH1D*       fHistoGammaCaloPt                                                                   = NULL;
    TH1D*       fHistoGammaCaloPtOrBin                                                              = NULL;
    TH1D*       fHistoGammaMCrecCaloPt                                                              = NULL;
    TH1D*       fHistoGammaMCrecCaloPtOrBin                                                         = NULL;
    TH1D*       fHistoGammaMCAllInEMCAccPt                                                          = NULL;
    TH1D*       fHistoGammaMCAllInEMCAccPtOrBin                                                     = NULL;
    TH1D*       fHistoGammaTrueCaloPt                                                               = NULL;
    TH1D*       fHistoGammaTrueCaloPtOrBin                                                          = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloPt                                                        = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloPtOrBin                                                   = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloPt                                                      = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloMCPt                                                    = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloPtOrBin                                                 = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloMCPtOrBin                                               = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloFromXPt[4]                                              = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaTrueSecondaryCaloFromXMCPt[4]                                            = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaTrueSecondaryCaloFromXPtOrBin[4]                                         = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaTrueSecondaryCaloFromXMCPtOrBin[4]                                       = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaTrueSecondaryCaloFromXFromEtaPt                                          = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloFromXFromEtaMCPt                                        = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloFromXFromEtaPtOrBin                                     = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloFromXFromEtaMCPtOrBin                                   = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloRestNoEtaPt                                             = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloRestNoEtaMCPt                                           = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloRestNoEtaPtOrBin                                        = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloRestNoEtaMCPtOrBin                                      = NULL;

    TH2D*       fHistoGammaTruePrimaryCalo_recPt_MCPt_MC                                            = NULL;
    TH2D*       fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin                                      = NULL;
    TH2D*       fHistoGammaTrueSecondaryFromXCaloUnConv_MCPt_recPt_MC[4]                            = { NULL, NULL, NULL, NULL };
    TH2D*       fHistoGammaTrueSecondaryFromXCaloConv_MCPt_recPt_MC[4]                              = { NULL, NULL, NULL, NULL };
    TH2D*       fHistoGammaTrueSecondaryFromXCalo_MCPt_recPt_MC[4]                                  = { NULL, NULL, NULL, NULL };
    TH2D*       fHistoGammaTrueSecondaryFromXCalo_MCPt_recPt_MC_Rebin[4]                            = { NULL, NULL, NULL, NULL };

    TH1D*       fHistoGammaTruePrimaryCaloMCPt                                                      = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloMCPtOrBin                                                 = NULL;
    TH1D*       fHistoFracAllGammaCaloToSec                                                         = NULL;
    TH1D*       fHistoFracAllGammaCaloToSecOrBin                                                    = NULL;
    TH1D*       fHistoFracAllGammaCaloToSecFromX[4]                                                 = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoFracAllGammaCaloToSecFromXOrBin[4]                                            = { NULL, NULL, NULL, NULL };

    TH1D*       fHistoGammaCaloMCPurity                                                             = NULL;
    TH1D*       fHistoGammaMCrecPrimaryCaloPt                                                       = NULL;
    TH1D*       fHistoGammaMCrecPrimaryCaloPtOrBin                                                  = NULL;
    TH1D*       fHistoGammaCaloMCTruePurity                                                         = NULL;
    TH1D*       fHistoGammaCaloMCTruePurityOrBin                                                    = NULL;
    TH1D*       fHistoGammaCaloMCRecoEff                                                            = NULL;
    TH1D*       fHistoGammaCaloMCPrimaryRecoEff                                                     = NULL;
    TH1D*       fHistoGammaCaloMCPrimaryRecoEffMCPt                                                 = NULL;
    TH1D*       fHistoGammaCaloMCBackground                                                         = NULL;

    TH1D*       fHistoGammaTrueCaloConvPt                                                           = NULL;
    TH1D*       fHistoGammaTrueCaloConvPtOrBin                                                      = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloConvPt                                                    = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloConvPtOrBin                                               = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloConvPt                                                  = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloConvPtOrBin                                             = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloConvFromXPt[4]                                          = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaTrueSecondaryCaloConvFromXPtOrBin[4]                                     = { NULL, NULL, NULL, NULL };
    TH2D*       fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC                                        = NULL;
    TH2D*       fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC_Rebin                                  = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloConvMCPt                                                  = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloConvMCPtOrBin                                             = NULL;

    TH1D*       fHistoGammaTrueCaloUnConvPt                                                         = NULL;
    TH1D*       fHistoGammaTrueCaloUnConvPtOrBin                                                    = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloUnConvPt                                                  = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloUnConvPtOrBin                                             = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloUnConvPt                                                = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloUnConvPtOrBin                                           = NULL;
    TH1D*       fHistoGammaTrueSecondaryCaloUnConvFromXPt[4]                                        = { NULL, NULL, NULL, NULL };
    TH1D*       fHistoGammaTrueSecondaryCaloUnConvFromXPtOrBin[4]                                   = { NULL, NULL, NULL, NULL };
    TH2D*       fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC                                      = NULL;
    TH2D*       fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC_Rebin                                = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloUnConvMCPt                                                = NULL;
    TH1D*       fHistoGammaTruePrimaryCaloUnConvMCPtOrBin                                           = NULL;

    TH2D*       f2DHistoGammaTrueSecondaryCaloUnConvPt                                              = NULL;
    TH2D*       f2DHistoGammaTrueSecondaryCaloUnConvMCPt                                            = NULL;
    TH2D*       f2DHistoGammaTrueSecondaryCaloConvPt                                                = NULL;
    TH2D*       f2DHistoGammaTrueSecondaryCaloConvMCPt                                              = NULL;
    TH2D*       f2DHistoGammaTrueSecondaryCaloPt                                                    = NULL;
    TH2D*       f2DHistoGammaTrueSecondaryCaloMCPt                                                  = NULL;

    //******************** Definition specific for pileup loading in PbPb 2.76Tev *************************
    TFile* fForPileUp                       = NULL;
    TString autoDetectedMainDirForPileUp    = "";
    TList *TopDirForPileUp                  = NULL;
    TList* HistosGammaConversionForPileUp   = NULL;


    //******************** Definition of functions ****************************
    void     RebinSpectrum                      (   TH1D*       Spectrum,
                                                    TString     NewName = ""                        );
    void     RebinSpectrumToDCAzDistBinning     (   TH1D*       Spectrum,
                                                    TString     NewName = ""                        );
    void     CalculatePileUpBackground          (   Bool_t      doMC                                );
    void     CalculateGammaCorrection           (                                                   );
    void     CalculatePileUpGammaCorrection     (                                                   );
    void     Initialize                         (   TString     setPi0,
                                                    TString     Energy,
                                                    Int_t       numberOfBins,
                                                    Int_t       mode,
                                                    Bool_t      addSig                              );
    void     SaveHistos                         (   Int_t       isMC,
                                                    TString     fCutID,
                                                    TString     fPrefix3,
                                                    Bool_t      PileUpCorrection                    );
    void     SaveCorrectionHistos               (   TString     fCutID,
                                                    TString     fPrefix3,
                                                    Bool_t      PileUpCorrection                    );
    void     SaveDCAHistos                      (   Int_t       isMC,
                                                    TString     fCutID,
                                                    TString     fPrefix3                            );
    void     FillDCAHistogramsFromTree          (   TTree*      dcaTree,
                                                    Bool_t      isMC                                );
    void     FillMassHistosArray                (   TH2D*       fGammaGammaInvMassVSPtDummy         );
    void     ProduceBckProperWeighting          (   TH2D*       fHistoBckZM,
                                                    TH2D*       fHistoMotherZM                      );
    Bool_t   CompareArrays                      (   Int_t       nEntriesA,
                                                    Double_t*   arrayA,
                                                    Int_t       nEntriesB,
                                                    Double_t*   arrayB                              );
    Bool_t   CalculatePileUpSubtractedDCAz      (   TH1D*       trueGamma,
                                                    TH1D*       trueSubGamma,
                                                    TH1D*       trueGammaX,
                                                    TH1D*       &trueSubGammaX                      );
    Bool_t   CalculateDCAzDistributionRatio     (   TH1D***     numerator,
                                                    TH1D***     denominator,
                                                    Int_t       categoryFirst,
                                                    Int_t       categoryLast,
                                                    TH1D*       &ratio                              );
    Bool_t   CalculateDCAzDistributionRatio     (   TH1D***     numerator,
                                                    TH1D****     denominator,
                                                    Int_t       categoryFirst,
                                                    Int_t       backgroundExtractionMethod,
                                                    Int_t       categoryLast,
                                                    TH1D*       &ratio                              );
    Bool_t    CalculatePileUpCorrectionFactor   (   TH1D*       ratioWithWithoutPileUp,
                                                    TH1D*       &pileupCorrectionFactor,
                                                    TF1*        &fitToRatio                         );
    Bool_t   LoadSecondariesFromCocktailFile    (   TString,
                                                    TString                                         );
    void     Delete                             (                                                   );
    void     FillMassMCTrueMesonHistosArrays    (   TH2D* fHistoTrueMesonPrimInvMassVSPtFill,
                                                    TH2D** fHistoTrueMesonSecInvMassVSPtFill        );
    void     CheckForNULLForPointer             (   TH1D* fDummy1                                   );
    void     PlotAdditionalDCAz                 (   Int_t       isMC,
                                                    TString     fCutID                              );
#endif
