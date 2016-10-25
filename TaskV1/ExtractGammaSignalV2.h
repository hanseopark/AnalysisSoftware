// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion


//!! -> what is commented here is defined in ExtractSignal.h

Double_t    fNEvnt                                                  = 0;
Bool_t      fEnablePCM                                              = 0;
Bool_t      fEnableCalo                                             = 0;

TString     fDate                                                   = "";
TString     fDirectPhoton                                           = "";

TString     fOutputDir                                              = "";
TString     fSuffix                                                 = "";
TString     fMeson                                                  = "";

TString     fEventCutNumber                                         = "";
TString     fGammaCutNumber                                         = "";
TString     fClusterCutNumber                                       = "";
TString     fElectronCutNumber                                      = "";
TString     fMesonCutNumber                                         = "";

TString     fCutSelectionRead                                       = ""; 
TString     fEventCutSelectionRead                                  = "";
TString     fGammaCutSelectionRead                                  = "";
TString     fMesonCutSelectionRead                                  = "";

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

// exemplary pt bins for DCAz distributions
Int_t       exemplaryLowPtBin                                       = 1;
Int_t       exemplaryHighPtBin                                      = 15;

// photon categories
TString     categoryName[4]                                         = {"all", "cat1", "cat2", "cat3"};
TString     backgroundExtractionMethod[3]                           = {"std", "var1", "var2"};
Color_t     backgroundColor[3]                                      = {kBlue+2, kGreen+3, kGray+2};

// ShowBackground arguments
Int_t nIterationsShowBackground[4]                                  = {0};
TString optionShowBackground[3]                                     = {""};

// binning
TH1D*       fDeltaPtDummy                                           = NULL;
                                                                      
// Histograms for conversion analysis
TH1D*       fHistoGammaConvPt                                                                   = NULL;
TH1D*       fHistoGammaConvPtOrBin                                                              = NULL;
TH1D*       fHistoGammaMCrecConvPt                                                              = NULL;
TH1D*       fHistoGammaMCrecConvPtOrBin                                                         = NULL;
TH1D*       fHistoGammaMCConvPt                                                                 = NULL;
TH1D*       fHistoGammaMCConvRSPt                                                               = NULL;
TH1D*       fHistoGammaMCAllPt                                                                  = NULL;
TH1D*       fHistoGammaMCAllPtOrBin                                                             = NULL;
TH1D*       fHistoGammaTrueConvPt                                                               = NULL;
TH1D*       fHistoGammaTrueConvPtOrBin                                                          = NULL;
TH1D*       fHistoGammaTruePrimaryConvPt                                                        = NULL;
TH1D*       fHistoGammaTruePrimaryConvPtOrBin                                                   = NULL;
TH1D*       fHistoGammaTrueSecondaryConvPt                                                      = NULL;
TH1D*       fHistoGammaTrueSecondaryConvPtOrBin                                                 = NULL;
TH1D*       fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt                                     = NULL;
TH1D*       fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin                                = NULL;
TH1D*       fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt                                  = NULL;
TH1D*       fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin                             = NULL;
TH1D**      fHistoGammaMCDecayPt                                                                = NULL;
TH1D*       fHistoGammaTruePrimaryConvMCPt                                                      = NULL;

TH1D*       fHistoFracAllGammaToSecOrBin                                                        = NULL;
TH1D*       fHistoFracAllGammaToSec                                                             = NULL;
TH1D*       fHistoFracAllGammaToSecFromXFromK0s                                                 = NULL;
TH1D*       fHistoFracAllGammaToSecFromXFromK0sOrBin                                            = NULL;
TH1D*       fHistoFracAllGammaToSecFromXFromLambda                                              = NULL;
TH1D*       fHistoFracAllGammaToSecFromXFromLambdaOrBin                                         = NULL;
TH1D*       fHistoGammaMCConvProb                                                               = NULL;
TH1D*       fHistoGammaMCPurity                                                                 = NULL;
TH1D*       fHistoGammaMCTruePurityOrBin                                                        = NULL;
TH1D*       fHistoGammaMCrecPrimaryConvPt                                                       = NULL;
TH1D*       fHistoGammaMCTruePurity                                                             = NULL;
TH1D*       fHistoGammaMCrecPrimaryConvPtOrBin                                                  = NULL;
TH1D*       fHistoGammaMCRecoEff                                                                = NULL;
TH1D*       fHistoGammaMCPrimaryRecoEff                                                         = NULL;
TH1D*       fHistoGammaMCPrimaryRecoEffMCPt                                                     = NULL;
TH1D*       fHistoGammaMCPrimaryRecoEffRSMCPt                                                   = NULL;
TH1D*       fHistoGammaMCBackground                                                             = NULL;

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
TH2F**      fTrueSecondaryPhotonFromXFromK0sPtDCAz                                              = NULL;
TH1D***     fTrueSecondaryGammaFromXFromK0sPtDCAzBins                                           = NULL;
TH1D***     fTrueSecondarySubGammaFromXFromK0sPtDCAzBins                                        = NULL;
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
TH1D*       fTrueSecondaryFromXFromK0sConvGammaPtPileUpAllCat                                   = NULL;
TH1D*       fTrueSecondaryFromXFromK0sConvGammaPtPileUp                                         = NULL;

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

TH1D*       fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat    = NULL;
TH1D*       fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning          = NULL;
TF1*        fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat = NULL;
TF1*        fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning       = NULL;

TH1D**      fESDGammaPileUpCorrFactorAllCat                                                     = NULL;
TH1D**      fESDGammaPileUpCorrFactor                                                           = NULL;
TH1D*       fMCrecGammaPileUpCorrFactorAllCat                                                   = NULL;
TH1D*       fMCrecGammaPileUpCorrFactor                                                         = NULL;
TH1D*       fTruePrimaryConvGammaPileUpCorrFactorAllCat                                         = NULL;
TH1D*       fTruePrimaryConvGammaPileUpCorrFactor                                               = NULL;
TH1D*       fTrueSecondaryConvGammaPileUpCorrFactorAllCat                                       = NULL;
TH1D*       fTrueSecondaryConvGammaPileUpCorrFactor                                             = NULL;
TH1D*       fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactorAllCat                           = NULL;
TH1D*       fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactor                                 = NULL;

TH1D*       fHistoFracAllGammaToSecPileUp                                                       = NULL;
TH1D*       fHistoFracAllGammaToSecFromXFromK0sPileUp                                           = NULL;
TH1D*       fHistoGammaMCPurityPileUp                                                           = NULL;
TH1D*       fHistoGammaMCrecPrimaryConvPtPileUp                                                 = NULL;
TH1D*       fHistoGammaMCTruePurityPileUp                                                       = NULL;
TH1D*       fHistoGammaMCRecoEffPileUp                                                          = NULL;
TH1D*       fHistoGammaMCPrimaryRecoEffPileUp                                                   = NULL;

TH1D*       fHistoPhotonIsSelected                                                              = NULL;

//******************* tagging outputHistograms ********************************************************
TH1D**      fHistoGconvGInvMassPtGConvBin                                                       = NULL;

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


// Histograms for calo analysis
TList*      CaloContainer                                           = NULL;
TH1D*       fHistoGammaCaloPt                                       = NULL;
TH1D*       fHistoGammaCaloPtOrBin                                  = NULL;
TH1D*       fHistoGammaMCrecCaloPt                                  = NULL;
TH1D*       fHistoGammaMCrecCaloPtOrBin                             = NULL;
TH1D*       fHistoGammaMCAllInEMCAccPt                              = NULL;
TH1D*       fHistoGammaMCAllInEMCAccPtOrBin                         = NULL;
TH1D*       fHistoGammaTrueCaloPt                                   = NULL;
TH1D*       fHistoGammaTrueCaloPtOrBin                              = NULL;
TH1D*       fHistoGammaTruePrimaryCaloPt                            = NULL;
TH1D*       fHistoGammaTruePrimaryCaloPtOrBin                       = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloPt                          = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloPtOrBin                     = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloFromK0sPt                   = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloFromK0sPtOrBin              = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloFromLambdaPt                = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloFromLambdaPtOrBin           = NULL;
TH2D*       fHistoGammaTruePrimaryCalo_recPt_MCPt_MC                = NULL;
TH2D*       fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin          = NULL;
TH1D*       fHistoGammaTruePrimaryCaloMCPt                          = NULL;
TH1D*       fHistoFracAllGammaCaloToSec                             = NULL;
TH1D*       fHistoFracAllGammaCaloToSecOrBin                        = NULL;
TH1D*       fHistoFracAllGammaCaloToSecFromK0s                      = NULL;
TH1D*       fHistoFracAllGammaCaloToSecFromK0sOrBin                 = NULL;
TH1D*       fHistoFracAllGammaCaloToSecFromLambda                   = NULL;
TH1D*       fHistoFracAllGammaCaloToSecFromLambdaOrBin              = NULL;

TH1D*       fHistoGammaCaloMCPurity                                 = NULL;
TH1D*       fHistoGammaMCrecPrimaryCaloPt                           = NULL;
TH1D*       fHistoGammaMCrecPrimaryCaloPtOrBin                      = NULL;
TH1D*       fHistoGammaCaloMCTruePurity                             = NULL;
TH1D*       fHistoGammaCaloMCTruePurityOrBin                        = NULL;
TH1D*       fHistoGammaCaloMCRecoEff                                = NULL;
TH1D*       fHistoGammaCaloMCPrimaryRecoEff                         = NULL;
TH1D*       fHistoGammaCaloMCPrimaryRecoEffMCPt                     = NULL;
TH1D*       fHistoGammaCaloMCBackground                             = NULL;

TH1D*       fHistoGammaTrueCaloConvPt                               = NULL;
TH1D*       fHistoGammaTrueCaloConvPtOrBin                          = NULL;
TH1D*       fHistoGammaTruePrimaryCaloConvPt                        = NULL;
TH1D*       fHistoGammaTruePrimaryCaloConvPtOrBin                   = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloConvPt                      = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloConvPtOrBin                 = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloConvFromK0sPt               = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloConvFromK0sPtOrBin          = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloConvFromLambdaPt            = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloConvFromLambdaPtOrBin       = NULL;
TH2D*       fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC            = NULL;
TH2D*       fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC_Rebin      = NULL;
TH1D*       fHistoGammaTruePrimaryCaloConvMCPt                      = NULL;

TH1D*       fHistoGammaTrueCaloUnConvPt                             = NULL;
TH1D*       fHistoGammaTrueCaloUnConvPtOrBin                        = NULL;
TH1D*       fHistoGammaTruePrimaryCaloUnConvPt                      = NULL;
TH1D*       fHistoGammaTruePrimaryCaloUnConvPtOrBin                 = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloUnConvPt                    = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloUnConvPtOrBin               = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloUnConvFromK0sPt             = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloUnConvFromK0sPtOrBin        = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloUnConvFromLambdaPt          = NULL;
TH1D*       fHistoGammaTrueSecondaryCaloUnConvFromLambdaPtOrBin     = NULL;
TH2D*       fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC          = NULL;
TH2D*       fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC_Rebin    = NULL;
TH1D*       fHistoGammaTruePrimaryCaloUnConvMCPt                    = NULL;



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
void     PlotAdditionalDCAz                 (   Int_t       isMC, 
                                                TString     fCutID                              );
void     PlotDCAzInPtBinsWithBack           (   TH1D**      ESDGammaPtDCAzBins, 
                                                TH1D**      ESDGammaPtDCAzBinsBack,
                                                TH1D**      ESDGammaPtDCAzBinsBackB, 
                                                TString     namePlot, 
                                                TString     nameCanvas, 
                                                TString     namePad, 
                                                TString     dateDummy, 
                                                TString     fMesonType,  
                                                Int_t       fRowPlot, 
                                                Int_t       fColumnPlot, 
                                                Int_t       fStartBinPtRange, 
                                                Int_t       fNumberPtBins, 
                                                Double_t*   fRangeBinsPt, 
                                                TString     fDecayChannel, 
                                                Bool_t      fMonteCarloInfo, 
                                                TString     textCent = "MinBias"                );
void    PlotDCAzInPtBinsWithBack            (   TH1D**      ESDGammaPtDCAzBins,
                                                TH1D**      ESDGammaPtDCAzBinsBack,
                                                TH1D**      ESDGammaPtDCAzBinsBackB,
                                                TString     namePlot,
                                                TString     nameCanvas,
                                                TString     namePad,
                                                TString     dateDummy,
                                                TString     fMesonType,
                                                Int_t       fStartBinPtRange,
                                                Int_t       fNumberPtBins,
                                                Double_t*   fRangeBinsPt,
                                                TString     fDecayChannel,
                                                Bool_t      fMonteCarloInfo,
                                                TString     textCent = "MinBias"                );
void    PlotDCAzInPtBinsWithBack            (   TH1D**      ESDGammaPtDCAzBins,
                                                TH1D***     ESDGammaPtDCAzBinsBack,
                                                TH1D**      ESDGammaPtDCAzBinsBackB,
                                                TString     namePlot,
                                                TString     nameCanvas,
                                                TString     namePad,
                                                TString     dateDummy,
                                                TString     fMesonType,
                                                Int_t       fStartBinPtRange,
                                                Int_t       fNumberPtBins,
                                                Double_t*   fRangeBinsPt,
                                                TString     fDecayChannel,
                                                Bool_t      fMonteCarloInfo,
                                                TString     textCent = "MinBias"                );
Int_t   CalculateNumberOfRowsForDCAzPlots   (   Int_t       numberOfPads,
                                                Int_t       numberOfColumns                     );
void    DrawDCAzHisto                       (   TH1*        histo1,
                                                TString     Title,
                                                TString     XTitle,
                                                TString     YTitle,
                                                Float_t     xMin,
                                                Float_t     xMax,
                                                Int_t       bck,
                                                Color_t     color = kBlue                       );
void     DrawFractionPerCat                 (   TH1D**      frac,
                                                TString     fOutputDir,
                                                TString     fPrefix,
                                                TString     fPrefix2,
                                                TString     fCutSelection,
                                                TString     fSuffix                             );
void     FillMassHistosArray                (   TH2D*       fGammaGammaInvMassVSPtDummy         );
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
