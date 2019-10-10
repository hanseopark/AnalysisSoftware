//****************************************************************************
//********* provided by Gamma Conversion Group *******************************
//********* $ALICE_ROOT/PWGGA/GammaConv **************************************
//***** https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion *******
//****************************************************************************

#ifndef GAMMACONV_ExtractSignalV2
#define GAMMACONV_ExtractSignalV2

    #include <fstream>
    #include <memory>
    #include <vector>
    #include "TF1.h"
    #include "TH1D.h"
    #include "TString.h"
    #include "TObjectTable.h"

    //****************************************************************************
    //******************** Global variable for setup of macro ********************
    //****************************************************************************
    TDatime     now;
    std::fstream fFileErrLog;
    std::fstream fFileDataLog;
    TString     fTextCent                                                   = "";
    TString     fEnergyFlag                                                 = "";
    TString     fPrefix                                                     = "";
    TString     fPrefix2                                                    = "";
    TString     fPeriodFlag                                                 = "";
    TString     ftextDayth                                                  = "";
    TString     fDate                                                       = "";
    TString     fCollisionSystem                                            = "";
    TString     fTextMeasurement                                            = "";
    TString     fDetectionProcess                                           = "";
    TString     fCutSelection                                               = "";
    TString     fBackgroundMultCutNumber                                    = "";
    Int_t       fMode                                                       = -1;
    Int_t       fModeHeavy                                                  = -1; // used to retain mode >100 info
    TString     fDirectPhoton                                               = "";
    TString     fThesis                                                     = "";
    Bool_t      fAdvancedMesonQA                                            = kFALSE;
    Bool_t      fEstimateTrainPileUp                                        = kFALSE;
    Bool_t      fEnableDCMeson                                              = kFALSE;
    Bool_t      fEnableDCCluster                                            = kFALSE;
    Bool_t      fUseRPBackground                                            = kFALSE;
    Bool_t      fNewMCOutput                                                = kFALSE;
    Bool_t      fEnableNormBckHistoComparisonToTrueBck                      = kFALSE;

    TString     fEventCutSelection                                          = "";
    TString     fGammaCutSelection                                          = "";
    TString     fClusterCutSelection                                        = "";
    TString     fElectronCutSelection                                       = "";
    TString     fMesonCutSelection                                          = "";

    Int_t       fCrysFitting                                                = 0;
    Int_t       fIsMC                                                       = 0;
    Int_t       fMesonId                                                    = 0;
    Float_t     fBackgroundMultNumber                                       = 0.;
    Double_t    fYMaxMeson                                                  = 0.9;
    Double_t    fMesonMassExpect                                            = 0;   // expected meson mass
    Double_t    fNEvents                                                    = 0;
    Double_t*   fPeakRange                                                  = nullptr;
    Double_t*   fFitRange                                                   = nullptr;
    Double_t*   fBGFitRange                                                 = nullptr;
    Double_t*   fMesonIntDeltaRange                                         = nullptr;
    TString     fcMonth[12]                                                 = {"Jan","Feb","Mar","Apr","May","Jun",
                                                                            "Jul","Aug","Sep","Oct","Nov","Dec"};
    TString     centralityString                                            = "";
    TH1D*       fCopySignal                                                 = nullptr;
    TH1D*       fCopyOnlyBG                                                 = nullptr;
    Double_t    fScalingFactorBck[7][4];
    Double_t    fScalingFactorBckFullPt;
    Double_t    fScalingFactorBckMidPt;
    Double_t    fCBAlpha                                                    = 0;
    Double_t    fCBn                                                        = 0;
    Float_t     pictDrawingCoordinatesFWHM[9]                               = {0.6, 0.8, 0.30, 0.04, 0.15,0.7, 0.1, 0.035,0};
    Int_t       fNRebinGlobal                                               = 2;
    Int_t       maxNSec                                                     = 4;
    TString     nameSecondaries[4]                                          = {"K0S", "Lambda", "K0L", "Rest"};
    TString     nameSecondariesCocktail[4]                                  = {"K0s", "Lambda", "K0l", "Rest"};
    TString     nameResonanceFeedDownContributions[15]                      = {"Eta", "Rho0", "EtaPrime", "Omega", "Rho+", "Rho-", "Phi", "JPsi", "Delta0",
                                                                                "Delta+", "K+", "K-", "Omega+", "Omega-", "K*(892)0"};
    TString     nameResonanceFeedDownContributionsCocktail[15]              = {"Eta", "rho0", "EtaPrime", "omega", "rho+", "rho-", "phi", "J/psi", "Delta0",
                                                                                "Delta+", "K+", "K-", "Omega+", "Omega-", "K*(892)0"};
    TString     nameIntRange[6]                                             = {"", "Wide", "Narrow", "Left", "LeftWide", "LeftNarrow"};
    TString     nameIntBck[6]                                               = {"Minpol2","Pluspol2","Minexp","Plusexp","Minexp2","Plusexp2"};
    TString     nameIntBckRatios[6]                                         = {"RatioMinpol2","RatioPluspol2","RatioMinexp","RatioPlusexp","RatioMinexp2","RatioPlusexp2"};
    TString     nameIntBckResult[3]                                         = {"pol2_normal","exp_normal","exp2_normal"};

    Bool_t      fDoJetAnalysis                                              = kFALSE;
    Bool_t      fUsingUnfolding_AsData                                      = kFALSE;
    Bool_t      fUsingUnfolding_Missed                                      = kFALSE;
    Bool_t      fUsingUnfolding_Reject                                      = kFALSE;
    Int_t       fNJetEvents                                                 = 0;

    //****************************************************************************
    //******************************** Output files ******************************
    //****************************************************************************
    TFile*      fOutput                                                     = nullptr;
    TFile*      fOutput1                                                    = nullptr;
    TFile*      fOutput2                                                    = nullptr;

    //****************************************************************************
    //******************************** Object names ******************************
    //****************************************************************************
    TString     ObjectNameTrue                                              = "";
    TString     ObjectNameTrueFull                                          = "";
    TString     ObjectNameTrueWOWeights                                     = "";
    TString     ObjectNameProfileWeights                                    = "";
    TString     ObjectNameTrueSec                                           = "";
    TString     ObjectNameTrueSecFromK0S                                    = "";
    TString     ObjectNameTrueSecFromK0L                                    = "";
    TString     ObjectNameTrueSecFromLambda                                 = "";
    TString     ObjectNameMCPi0Acc                                          = "";
    TString     ObjectNameMCPi0AccWOWeights                                 = "";
    TString     ObjectNameMCPi0AccWOEvtWeights                              = "";
    TString     ObjectNameMCEtaAcc                                          = "";
    TString     ObjectNameMCEtaAccWOWeights                                 = "";
    TString     ObjectNameMCEtaAccWOEvtWeights                              = "";
    TString     ObjectNameMCEtaPrimeAcc                                     = "";
    TString     ObjectNameMCEtaPrimeAccWOWeights                            = "";
    TString     ObjectNameMCEtaPrimeAccWOEvtWeights                         = "";
    TString     ObjectNameMCPi0                                             = "";
    TString     ObjectNameMCPi0WOWeights                                    = "";
    TString     ObjectNameMCPi0WOEvtWeights                                 = "";
    TString     ObjectNameMCEta                                             = "";
    TString     ObjectNameMCEtaWOWeights                                    = "";
    TString     ObjectNameMCEtaWOEvtWeights                                 = "";
    TString     ObjectNameMCEtaPrime                                        = "";
    TString     ObjectNameMCEtaPrimeWOWeights                               = "";
    TString     ObjectNameMCEtaPrimeWOEvtWeights                            = "";
    TString     ObjectNameTrueGGBck                                         = "";
    TString     ObjectNameTrueContBck                                       = "";
    TString     ObjectNameTrueBckFullMesonContained                         = "";
    TString     ObjectNameTrueBckAsymEClus                                  = "";
    TString     ObjectNameTrueAllBck                                        = "";
    TString     ObjectNameTrueCaloPhoton                                    = "";
    TString     ObjectNameTrueCaloConvPhoton                                = "";
    TString     ObjectNameTrueCaloMerged                                    = "";
    TString     ObjectNameTrueMixedCaloConvPhoton                           = "";
    TString     ObjectNameTrueCaloMergedPartConv                            = "";
    TString     ObjectNameK0sRecPi0                                         = "";
    TString     ObjectNameLambdaRecPi0                                      = "";
    TString     ObjectNameDCMesonInvMassPt                                  = "";
    TString     ObjectNameDCGammaClusPt                                     = "";
    TString     ObjectNameMesonMultipleCount                                = "";
    TString     ObjectNameGammaClusMultipleCount                            = "";
    // object names for secondary MC histos
    TString     ObjectNameMCSecPi0Acc                                       = "";
    TString     ObjectNameMCSecPi0                                          = "";

    TString     fNameHistoGG                                                = "";
    TString     fNameHistoBack                                              = "";
    TString     fNameHistoPP                                                = "";
    TString     fNameHistoBackNormOut                                       = "";
    TString     fNameHistoBackNormLeftOut                                   = "";
    TString     fNameHistoTrue                                              = "";
    TString     fNameHistoTrueGGBck                                         = "";
    TString     fNameHistoTrueContBck                                       = "";
    TString     fNameHistoTrueAllBck                                        = "";
    TString     fNameHistoTrueMesonContained                                = "";
    TString     fNameHistoTrueAsymEClus                                     = "";
    TString     fNameHistoTrueSec                                           = "";
    TString     fNameHistoEffi                                              = "";
    TString     fNameHistoFrac                                              = "";
    TString     fNameHistoMotherZM                                          = "";
    TString     fNameHistoBckZM                                             = "";
    TString     fNameFitSignalPos                                           = "";

    //****************************************************************************
    //************************ global variable ranges ****************************
    //****************************************************************************
    Double_t*    fBGFitRangeLeft                                            = nullptr;
    Double_t*    fMesonPlotRange                                            = nullptr;
    Double_t*    fIntFixedRange                                             = nullptr;
    Double_t*    fMesonIntDeltaRangeWide                                    = nullptr;
    Double_t*    fMesonIntDeltaRangeNarrow                                  = nullptr;
    Double_t*    fMesonCurIntRange[6]                                       = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*    fMesonTrueIntRange[3]                                      = {nullptr, nullptr, nullptr};
    Double_t*    fMesonTrueIntReweightedRange[3]                            = {nullptr, nullptr, nullptr};
    Double_t*    fMesonTrueIntUnweightedRange[3]                            = {nullptr, nullptr, nullptr};
    Double_t*    fMesonMassRange                                            = nullptr;
    Double_t*    fMesonMassPlotRange                                        = nullptr;
    Double_t*    fMesonFitRange                                             = nullptr;
    Double_t*    fMesonFitRangeWithoutPeak                                  = nullptr;
    Double_t     fMesonWidthExpect                                          = 0;
    Double_t     fMesonWidthExpectMC                                        = 0;
    Double_t     fMesonLambdaTail                                           = 0;
    Double_t     fMesonLambdaTailMC                                         = 0;
    Double_t*    fMesonWidthRange                                           = nullptr;
    Double_t*    fMesonWidthRangeMC                                         = nullptr;
    Double_t*    fMesonLambdaTailRange                                      = nullptr;
    Double_t*    fMesonLambdaTailRangeNominal                               = nullptr;
    Double_t*    fMesonLambdaTailRangeMC                                    = nullptr;
    Double_t*    fMidPt                                                     = nullptr;
    Double_t*    fFullPt                                                    = nullptr;

    //****************************************************************************
    //************************ correction histograms *****************************
    //****************************************************************************
    TH1D*     fDeltaPt                                                        = nullptr;

    //****************************************************************************
    //******************************** Functions *********************************
    //****************************************************************************
    void ProcessEM_switch(TH1D*,TH1D*,Double_t *);                                                              // Function that calls different Normalization of Signal and BG Methods
    void ProcessEM(TH1D*,TH1D*,Double_t *);                                                                     // Normalization of Signal and BG
    Double_t FitFunctionPHOSBck(Double_t *, Double_t *);
    Double_t FitFunctionPHOSBckPol1(Double_t *, Double_t *);
    Double_t FitFunctionPHOSBckPol2(Double_t *, Double_t *);
    void ProcessEM_FitBins(TH1D*,TH1D*,Double_t *);                                                             // Normalization of Signal and BG For PHOS
    void ProcessRatioSignalBackground(TH1D* , TH1D* );                                                          // Calculate Ratio of Signal and BG in momentum slices
    void FillMassHistosArray(TH2D*);                                                                            // Fill invariant mass histograms for Signal and Backgrounf
    TH1D* FillProjectionX (TH2*, TString, Double_t, Double_t, Int_t);                                           // Fill Projection in according to Y bins
    void FillMassMCTrueMesonHistosArray(TH2D*);                                                                 // Fill invariant mass histograms for validated mesons
    void FillMassMCTrueFullMesonHistosArray(TH2D* );                                                            // Fill invariant mass histograms for validated mesons full momentum range
    void FillMassMCTrueMesonDCHistosArray(TH2D*);                                                               // Fill invariant mass histograms for double counted validated mesons
    void FillMassMCTrueMesonCaloPhotonHistosArray(TH2D*);                                                       // Fill invariant mass histograms for validated mesons true calo photons
    void FillMassMCTrueMesonCaloConvPhotonHistosArray(TH2D*);                                                   // Fill invariant mass histograms for validated mesons true calo conv photons
    void FillMassMCTrueMesonCaloElectronHistosArray(TH2D*);                                                     // Fill invariant mass histograms for validated mesons true calo electrons
    void FillMassMCTrueMesonMixedCaloConvPhotonHistosArray(TH2D*);                                              // Fill invariant mass histograms for validated mesons true calo photons, conv photons
    void FillMassMCTrueMesonCaloMergedClusterHistosArray(TH2D*);                                                // Fill invariant mass histograms for validated mesons true calo merged
    void FillMassMCTrueMesonCaloMergedClusterPartConvHistosArray(TH2D*);                                        // Fill invariant mass histograms for validated mesons true calo merged part conv
    void FillMassMCTrueReweightedMesonHistosArray(TH2D*);                                                       // Fill invariant mass histograms for validated mesons reweighted
    void FillMassMCTrueUnweightedMesonHistosArray(TH2D*);                                                       // Fill invariant mass histograms for validated mesons unweighted
    void FillMassMCTrueGGBckHistosArray(TH2D*);                                                                 // Fill invariant mass histograms for validated gamma gamma BG
    void FillMassMCTrueContBckHistosArray(TH2D*);                                                               // Fill invariant mass histograms for validated contamination
    void FillMassMCTrueAllBckHistosArray(TH2D*);                                                                // Fill invariant mass histograms for validated BG
    void FillMassMCTrueFullMesonContainedHistosArray(TH2D*);                                                    // Fill invariant mass histograms for validated fully contained mesons in one cluster
    void FillMassMCTrueAsymEClusHistosArray(TH2D*);                                                             // Fill invariant mass histograms for validated asym Eclusters
    void FillMassMCTrueSecMesonHistosArray(TH2D**);                                                             // Fill invariant mass histograms for validated secondaries
    TH1D* CalculateSecondaryFractions(TH1D* histoRawYield, TH1D* histoRawYieldSec, TString nameHistoFrac);      // Calculate fraction of secondaries
    void CreatePtHistos();                                                                                      // Creat pt dependent histograms
    void FillPtHistos();                                                                                        // Fill pt dependent histograms
    void FitSubtractedInvMassInPtBins(TH1D * ,Double_t *, Int_t, Bool_t );                                      // Fits the invariant mass histos with a gaussian plus exponential plus lin BG
    void FitSubtractedPol2InvMassInPtBins(TH1D*, Double_t*, Int_t, Bool_t );                                    // Fits the invariant mass histos with a gaussian plus exponential plus pol2 BG
    void FitSubtractedExp1InvMassInPtBins(TH1D* , Double_t* , Int_t , Bool_t );                                 // Fits the invariant mass histos with a gaussian plus exponential plus exp BG
    void FitSubtractedExp2InvMassInPtBins(TH1D* , Double_t* , Int_t , Bool_t );                                 // Fits the invariant mass histos with a gaussian plus exponential plus exp BG
    void FitSubtractedPureGaussianInvMassInPtBins(TH1D*, Int_t);                                                // Fits the invariant mass histos with a pure gaussian plus lin BG
    void FitTrueInvMassInPtBins(TH1D * ,Double_t *, Int_t, Bool_t);                                             // Fits the true invariant mass histos with a gaussian plus exponential plus lin BG
    void FitTrueInvMassPureGaussianInPtBins(TH1D * , Int_t);                                                    // Fits the true invariant mass histos with a gaussian plus lin BG
    void FitCBSubtractedInvMassInPtBins(TH1D* ,Double_t * , Int_t ,Bool_t,TString, Bool_t );                    // Fits the invariant mass histos with a CB function
    void ProduceBckProperWeighting(TList*, TList*, TList*, TList* ,Bool_t);                                     // Create BG with proper weighting
    void ProduceBckWithoutWeighting(TH2D *);                                                                    // Create BG without proper weighting
    void IntegrateHistoInvMassStream(TH1D * , Double_t *);                                                      // Integrate invariant mass histogram with output to ifstream
    void IntegrateHistoInvMass(TH1D * , Double_t *);                                                            // Integrate invariant mass histogram
    void IntegrateFitFunc(TF1 * , TH1D *, Double_t *);                                                          // Integrate fit function
    void IntegrateFitFuncAndError(TF1 * , TH1D *, Double_t *);                                                          // Integrate fit function
    void FillHistosArrayMC(TH1D* , TH1D*, TH1D*);                                                               // Fill MC input histograms
    void FillHistosArrayMCWOWeights(TH1D* , TH1D*, TH1D*);                                                      // Fill MC input histograms
    void FillHistosArrayMCWOEvtWeights(TH1D* , TH1D*, TH1D*);                                                   // Fill MC input histograms
    void CalculateMesonAcceptance();                                                                            // Calculation of meson acceptance
    void CalculateMesonAcceptanceWOWeights();                                                                   // Calculation of meson acceptance
    void CalculateMesonAcceptanceWOEvtWeights();                                                                // Calculation of meson acceptance
    TH1D* CalculateMesonEfficiency(TH1D*, TH1D**, TH1D*, TString);                                              // Calculation of meson efficiencies
    void SaveHistos(Int_t, TString, TString, Bool_t);                                                           // Saving standard histograms to a file
    void SaveCorrectionHistos(TString , TString);                                                               // Saving correction histograms to a file
    void Initialize(TString setPi0, Int_t, Int_t);                                                              // Initialization of global variables depending on meson analysed
    void CalculateFWHM(TF1 *);                                                                                  // Calculation of FWHM
    Double_t CrystalBallBck(Double_t *,Double_t *);                                                             // Definition of CrystalBall with linear BG
    Double_t LinearBackground(Double_t *,Double_t *);                                                           // Definition of linear BG
    Double_t LinearBGExclusion(Double_t *,Double_t *);                                                          // Definition of linear BG with excluded region
    Double_t LinearBGExclusionnew(Double_t *,Double_t *);                                                       // Definition of linear BG with excluded region
    Double_t CrystalBall(Double_t *,Double_t *);                                                                // Definition of CrystalBall
    void Delete();                                                                                              // Deleting all pointers
    void SetCorrectMCHistogrammNames(TString);                                                                  // Setting correct histogram names
    void FillMCSecondaryHistAndCalculateAcceptance(TH2D*, TH2D*);                                               // Fill secondary MC input histograms and calculate respective acceptance
    Bool_t LoadSecondaryPionsFromExternalFile();                                                                // Loads secondary neutral pion input graphs from file
    Bool_t LoadSecondaryPionsFromCocktailFile(TString, TString);                                                // Loads secondary neutral pion input graphs from file

    Double_t fitGaussianPileUp(Double_t *x, Double_t *par);
    Double_t fitGaussianPileUp2(Double_t *x, Double_t *par);
    Int_t GetHeavyMesonDigit(TString);
    void PlotJetPlots(TH1D*, TString, TString, TString, TString, TString, TString, Bool_t, Bool_t, Bool_t);
    void PlotJetPlots(TH1D*, TH1D*, TString, TString, TString, TString, TString, TString, Bool_t, Bool_t, Bool_t);
    void PlotJetPlots(TH1F*, TString, TString, TString, TString, TString, TString, Bool_t, Bool_t, Bool_t);
    void PlotJetPlots(TH2D*, TString, TString, TString, TString, TString, Bool_t);

    //****************************************************************************
    //************************** input histograms ********************************
    //****************************************************************************
    TH2D*       fHistoTrueMesonInvMassVSPt                                  = nullptr;
    TH2D*       fHistoTrueFullMesonInvMassVSPt                              = nullptr;
    TH2D*       fHistoTrueMesonInvMassVSPtWOWeights                         = nullptr;
    TH2D*       fHistoTrueMesonInvMassVSPtReweighted                        = nullptr;
    TH2D*       fHistoTrueMesonInvMassVSPtUnweighted                        = nullptr;
    TProfile2D* fProfileTrueMesonInvMassVSPtWeights                         = nullptr;
    TH2D*       fHistoTrueMesonInvMassVSPtSec                               = nullptr;
    TH2D*       fHistoTrueGGBckInvMassVSPt                                  = nullptr;
    TH2D*       fHistoTrueFullMesonContainedInvMassVSPt                     = nullptr;
    TH2D*       fHistoTrueAsymEClusInvMassVSPt                              = nullptr;
    TH2D*       fHistoTrueContBckInvMassVSPt                                = nullptr;
    TH2D*       fHistoTrueAllBckInvMassVSPt                                 = nullptr;
    TH2D*       fHistoTrueSecMesonInvMassVSPt[4]                            = { nullptr, nullptr, nullptr, nullptr };
    TH2D*       fHistoTrueMesonCaloPhotonInvMassVSPt                        = nullptr;
    TH2D*       fHistoTrueMesonCaloConvPhotonInvMassVSPt                    = nullptr;
    TH2D*       fHistoTrueMesonCaloElectronInvMassVSPt                      = nullptr;
    TH2D*       fHistoTrueMesonMergedClusterInvMassVSPt                     = nullptr;
    TH2D*       fHistoTrueMesonMergedClusterPartConvInvMassVSPt             = nullptr;
    TH2D*       fHistoTrueMesonMixedCaloConvPhotonInvMassVSPt               = nullptr;
    TH2D*       fGammaGammaInvMassVSPt                                      = nullptr;
    TH2D*       fBckInvMassVSPt                                             = nullptr;
    TH1D*       fNumberOfGoodESDTracks                                      = nullptr;
    TH1D*       fEventQuality                                               = nullptr;
    TH1D*       fHistoYieldK0sWithPi0DaughterRec                            = nullptr;
    TH1D*       fHistoYieldLambdaWithPi0DaughterRec                         = nullptr;
    TH2D*       fPi0ResponseMatrix                                          = nullptr;
    TH2D*       fPi0ResponseMatrixRebin                                     = nullptr;
    TH2D*       fHistoMotherZM                                              = nullptr;
    TH2D*       fHistoBckZM                                                 = nullptr;
    TH1D*       fHistJetPt                                                  = nullptr;
    TH1F*       fHistJetEta                                                 = nullptr;
    TH1F*       fHistJetPhi                                                 = nullptr;
    TH1D*       fHistJetArea                                                = nullptr;
    TH1D*       fHistNJetsEvents                                            = nullptr;
    TH1D*       fHistNEventswithJets                                        = nullptr;
    TH1D*       fHistRatioPtPi0Jet                                          = nullptr;
    TH1D*       fHistDoubleCounting                                         = nullptr;
    TH2D*       fHistGammaGammaPi0JetEvent                                  = nullptr;
    TH2D*       fHistRPi0Jet                                                = nullptr;
    TH2D*       fHistEtaPhiPi0Jet                                           = nullptr;
    TH2D*       fHistEtaPhiPi0inJet                                         = nullptr;
    TH1D*       fHistoDoubleCountTruePi0                                    = nullptr;
    TH1D*       fHistoDoubleCountTrueEta                                    = nullptr;
    TH2D*       fHistoJetUnfold                                             = nullptr;
    TH2D*       fHistFragmFunct                                             = nullptr;
    TH2D*       fHistFragmFunctChargedPart                                  = nullptr;
    TH2D*       fHistoTruePi0FragmFunc                                      = nullptr;
    TH2D*       fHistoTruePi0FragmFuncChargPart                             = nullptr;
    TH2D*       fHistoTrueEtaFragmFunc                                      = nullptr;
    TH2D*       fHistoTrueEtaFragmFuncChargPart                             = nullptr;
    TH2D*       fHistFragmZInvMass                                          = nullptr;
    TH2D*       fHistoTruePi0FramZInvMass                                   = nullptr;
    TH2D*       fHistoTrueEtaFramZInvMass                                   = nullptr;


    //****************************************************************************
    //************************** background histograms ***************************
    //****************************************************************************
    TH2F**      fHistoWeightsBGZbinVsMbin                                   = 0x0;
    TH2F**      fHistoFillPerEventBGZbinVsMbin                              = 0x0;
    TH2F**      fHistoWeightsBGZbinVsPsibin                                 = 0x0;
    TH2F**      fHistoFillPerEventBGZbinVsPsibin                            = 0x0;

    //****************************************************************************
    //************************ sample histograms for inv Mass ********************
    //****************************************************************************
    Int_t       iBckSwitch                                                  = 0;
    Int_t       iNumberOfOtherSigToBckRatioFits                             = 1; //If u change this value remember to also change: fHistoChi2SigToBckFit, fHistoChi2SigToBckFit, colorFitSigToBckFit, styleFitSigToBckFit, fSigToBckFitChi2, fFitPHOSAllOtherSigToBckFits
    TH1D*       fBckNorm                                                    = nullptr;
    TH1D*       fSignal                                                     = nullptr;
    TH1D*       fRatioSB                                                    = nullptr;
    TH1D*       fFittingHistMidPtSignal                                     = nullptr;
    TH1D*       fFittingHistMidPtBackground                                 = nullptr;
    TH1D*       fFittingHistMidPtSignalSub                                  = nullptr;
    TH1D*       fMesonFullPtSignal                                          = nullptr;
    TH1D*       fMesonFullPtBackground                                      = nullptr;
    TH1D*       fMesonFullPtBackNorm                                        = nullptr;
    TF1 *       fFitSignalInvMassMidPt                                      = nullptr;
    TF1 *       fFitSignalInvMassMidPt2                                     = nullptr;

    //****************************************************************************
    //************************ inv mass histograms in pT bins ********************
    //****************************************************************************
    TH1D*       fHistoMotherZMProj                                          = nullptr;
    TH1D*       fHistoBckZMProj                                             = nullptr;
    TH1D*       fHistoMotherZMProjFullPt                                    = nullptr;
    TH1D*       fHistoBckZMProjFullPt                                       = nullptr;
    TH1D*       fHistoMotherZMProjMidPt                                     = nullptr;
    TH1D*       fHistoBckZMProjMidPt                                        = nullptr;
    TH1D*       fHistoMotherZPsiProj                                        = nullptr;
    TH1D*       fHistoBckZPsiProj                                           = nullptr;
    TH1D*       fHistoMotherZPsiProjFullPt                                  = nullptr;
    TH1D*       fHistoBckZPsiProjFullPt                                     = nullptr;
    TH1D*       fHistoMotherZPsiProjMidPt                                   = nullptr;
    TH1D*       fHistoBckZPsiProjMidPt                                      = nullptr;
    TH1D**      fHistoMappingTrueMesonInvMassPtBins                         = nullptr;
    TH1D**      fHistoMappingTrueFullMesonInvMassPtBins                     = nullptr;
    TH1D**      fHistoMappingTrueMesonInvMassPtReweightedBins               = nullptr;
    TH1D**      fHistoMappingTrueMesonInvMassPtUnweightedBins               = nullptr;
    TH1D**      fHistoMappingTrueGGBckInvMassPtBins                         = nullptr;
    TH1D**      fHistoMappingTrueContBckInvMassPtBins                       = nullptr;
    TH1D**      fHistoMappingTrueAllBckInvMassPtBins                        = nullptr;
    TH1D**      fHistoMappingTrueMesonContainedInvMassPtBins                = nullptr;
    TH1D**      fHistoMappingTrueAsymEClusInvMassPtBins                     = nullptr;
    TH1D**      fHistoMappingTrueSecMesonInvMassPtBins[4]                   = {nullptr, nullptr, nullptr, nullptr};
    TH1D**      fHistoMappingGGInvMassPtBin                                 = nullptr;
    TH1D**      fHistoMappingBackInvMassPtBin                               = nullptr;
    TH1D**      fHistoMappingBackNormInvMassPtBin                           = nullptr;
    TH1D**      fHistoMappingBackNormAndRemainingBGInvMassPtBin             = nullptr;
    TH1D**      fHistoMappingRatioSBInvMassPtBin                            = nullptr;
    TH1D**      fHistoMappingSignalInvMassPtBin                             = nullptr;
    TH1D**      fHistoMappingRemainingBGInvMassPtBin                        = nullptr;
    TH1D**      fHistoMappingSignalRemainingBGSubInvMassPtBin               = nullptr;
    TH1D**      fHistoMappingBackNormInvMassLeftPtBin                       = nullptr;
    TH1D**      fHistoMappingSignalInvMassLeftPtBin                         = nullptr;
    TH1D**      fHistoMappingSignalRemainingBGSubInvMassLeftPtBin           = nullptr;
    TH1D**      fHistoMappingRemainingBGInvMassLeftPtBin                    = nullptr;
    TH1D**      fHistoMappingTrueMesonCaloPhotonInvMassPtBins               = nullptr;
    TH1D**      fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins           = nullptr;
    TH1D**      fHistoMappingTrueMesonCaloElectronInvMassPtBins             = nullptr;
    TH1D**      fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins        = nullptr;
    TH1D**      fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins= nullptr;
    TH1D**      fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins      = nullptr;
    TH1D**      fHistoMapping_Fragm_ZInvMassZBin                            = nullptr;

    //****************************************************************************
    //**************************** global fit functions **************************
    //****************************************************************************
    TF1 *       fFitReco                                                    = nullptr;
    TF1 *       fFitGausExp                                                 = nullptr;
    TF1 *       fFitLinearBck                                               = nullptr;
    TF1 *       fFitLinearBckExcl                                           = nullptr;
    TF1 *       fFitLinearBckOut                                            = nullptr;
    TF1 *       fFitPHOSPol1                                                = nullptr;
    TF1 *       fFitPHOSPol2                                                = nullptr;
    Double_t    fYields;
    Double_t    fYieldsError;
    Double_t    fFWHMFunc;
    Double_t    fFWHMFuncError;
    Double_t    fYieldsFunc;
    Double_t    fYieldsFuncError;
    Double_t    fIntLinearBck;
    Double_t    fIntLinearBckError;
    Double_t    fIntLinearBckOut;
    Double_t    fIntLinearBckErrorOut;

    //****************************************************************************
    //******************** yield extraction standard window **********************
    //****************************************************************************
    Double_t*   fGGYields[6]                                                = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fBckYields[6]                                               = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fTotalBckYields[6]                                          = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYields[6]                                             = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fGGYieldsError[6]                                           = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fBckYieldsError[6]                                          = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fTotalBckYieldsError[6]                                     = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsError[6]                                        = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldMeson[6]                                         = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldDiffBck[6]                                       = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldDiffBckRatios[6]                                 = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldDiffBckResult[3]                                 = {nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldMesonPerEvent[6]                                 = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldTrueMeson[6]                                     = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldTrueMesonFromFit[6]                              = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsFunc[6]                                         = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsResidualBckFunc[6]                              = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsCorResidualBckFunc[6]                           = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsPerEvent[6]                                     = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsFuncError[6]                                    = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsResidualBckFuncError[6]                         = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsCorResidualBckFuncError[6]                      = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsPerEventError[6]                                = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

    Double_t*   fMesonTrueYields[3]                                         = {nullptr, nullptr, nullptr};
    Double_t*   fMesonTrueYieldsError[3]                                    = {nullptr, nullptr, nullptr};
    Double_t*   fMesonTrueYieldsFromFit[3]                                  = {nullptr, nullptr, nullptr};
    Double_t*   fMesonTrueYieldsFromFitError[3]                             = {nullptr, nullptr, nullptr};
    Double_t*   fMesonTrueYieldsReweighted[3]                               = {nullptr, nullptr, nullptr};
    Double_t*   fMesonTrueYieldsReweightedError[3]                          = {nullptr, nullptr, nullptr};
    Double_t*   fMesonTrueYieldsUnweighted[3]                               = {nullptr, nullptr, nullptr};
    Double_t*   fMesonTrueYieldsUnweightedError[3]                          = {nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldTrueMesonReweighted[3]                           = {nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldTrueMesonUnweighted[3]                           = {nullptr, nullptr, nullptr};

    Double_t*   fMesonYieldsResBckOtherFunc[3]                              = { nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsResBckOtherFuncError[3]                         = { nullptr, nullptr, nullptr};

    TH1D*       fHistoYieldMesonPerJetEvent[6]                              = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsPerJetEvent[6]                                  = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    Double_t*   fMesonYieldsPerJetEventError[6]                             = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

    //****************************************************************************
    //************* histos, fits, doubles for pure Standard fitting **************
    //****************************************************************************
    TF1**       fFitSignalInvMassPtBin                                      = nullptr;
    TF1**       fFitSignalWithOtherBGInvMassPtBin[3]                        = { nullptr, nullptr, nullptr };
    TF1**       fFitRemainingBGInvMassPtBin                                 = nullptr;
    TF1**       fFitBckInvMassPtBin                                         = nullptr;
    TF1**       fFitBckOtherInvMassPtBin[3]                                 = { nullptr, nullptr, nullptr };
    TF1**       fFitTrueSignalInvMassPtBin                                  = nullptr;
    TF1**       fFitTrueSignalInvMassPtReweightedBin                        = nullptr;
    TF1**       fFitTrueSignalInvMassPtUnweightedBin                        = nullptr;
    TF1**       fFitPHOSAllOtherSigToBckFits[1]                             = { nullptr };
    TF1**       fFitPHOSPol2PtBin                                           = nullptr;
    Double_t*   fMesonMass                                                  = nullptr;
    Double_t*   fMesonMassError                                             = nullptr;
    Double_t*   fMesonTrueMass                                              = nullptr;
    Double_t*   fMesonTrueMassError                                         = nullptr;
    Double_t*   fMesonFWHM                                                  = nullptr;
    Double_t*   fMesonFWHMError                                             = nullptr;
    Double_t*   fMesonTrueFWHM                                              = nullptr;
    Double_t*   fMesonTrueFWHMError                                         = nullptr;
    Double_t*   fMesonLambdaTailpar                                         = nullptr;
    Double_t*   fMesonLambdaTailparError                                    = nullptr;
    Double_t*   fMesonLambdaTailMCpar                                       = nullptr;
    Double_t*   fMesonLambdaTailMCparError                                  = nullptr;
    Double_t*   fMesonAmplitudepar                                          = nullptr;
    Double_t*   fMesonAmplitudeparError                                     = nullptr;
    Double_t*   fMesonSigmapar                                              = nullptr;
    Double_t*   fMesonSigmaparError                                         = nullptr;
    Double_t*   fMesonTrueSigmapar                                          = nullptr;
    Double_t*   fMesonTrueSigmaparError                                     = nullptr;
    Double_t*   fMesonResidualBGlin                                         = nullptr;
    Double_t*   fMesonResidualBGlinError                                    = nullptr;
    Double_t*   fMesonResidualBGcon                                         = nullptr;
    Double_t*   fMesonResidualBGconError                                    = nullptr;
    Double_t*   fMesonChi2[4]                                               = {nullptr, nullptr, nullptr, nullptr};
    Double_t*   fSigToBckFitChi2[2]                                         = {nullptr, nullptr};


    TH1D*       fHistoMassMeson                                             = nullptr;
    TH1D*       fHistoFWHMMeson                                             = nullptr;
    TH1D*       fHistoTrueMassMeson                                         = nullptr;
    TH1D*       fHistoTrueFWHMMeson                                         = nullptr;
    Double_t*   fMesonTrueMassReweighted                                    = nullptr;
    Double_t*   fMesonTrueFWHMReweighted                                    = nullptr;
    Double_t*   fMesonTrueMassReweightedError                               = nullptr;
    Double_t*   fMesonTrueFWHMReweightedError                               = nullptr;
    TH1D*       fHistoTrueMassMesonReweighted                               = nullptr;
    TH1D*       fHistoTrueFWHMMesonReweighted                               = nullptr;
    Double_t*   fMesonTrueMassUnweighted                                    = nullptr;
    Double_t*   fMesonTrueFWHMUnweighted                                    = nullptr;
    Double_t*   fMesonTrueMassUnweightedError                               = nullptr;
    Double_t*   fMesonTrueFWHMUnweightedError                               = nullptr;
    TH1D*       fHistoTrueMassMesonUnweighted                               = nullptr;
    TH1D*       fHistoTrueFWHMMesonUnweighted                               = nullptr;

    //****************************************************************************
    //********* histos, fits, doubles for pure Standard fitting left *************
    //****************************************************************************
    TF1**       fFitInvMassLeftPtBin                                        = nullptr;
    TF1**       fFitInvMassLeftPtBin2                                       = nullptr;
    TF1**       fFitRemainingBGInvMassLeftPtBin                             = nullptr;
    TF1**       fFitRemainingBGInvMassLeftPtBin2                            = nullptr;
    TF1**       fFitBckInvMassLeftPtBin                                     = nullptr;
    TF1**       fFitBckInvMassLeftPtBin2                                    = nullptr;
    Double_t*   fMesonMassLeft                                              = nullptr;
    Double_t*   fMesonFWHMLeft                                              = nullptr;
    Double_t*   fMesonMassLeftError                                         = nullptr;
    Double_t*   fMesonFWHMLeftError                                         = nullptr;
    TH1D*       fHistoMassMesonLeft                                         = nullptr;
    TH1D*       fHistoFWHMMesonLeft                                         = nullptr;

    //****************************************************************************
    //************** histos, fits, doubles for pure Gaussian fitting *************
    //****************************************************************************
    TF1**       fFitSignalGaussianInvMassPtBin                              = nullptr;
    TF1**       fFitTrueSignalGaussianInvMassPtBin                          = nullptr;
    Double_t*   fMesonMassGaussian                                          = nullptr;
    Double_t*   fMesonMassGaussianError                                     = nullptr;
    Double_t*   fMesonTrueMassGaussian                                      = nullptr;
    Double_t*   fMesonTrueMassGaussianError                                 = nullptr;
    Double_t*   fMesonWidthGaussian                                         = nullptr;
    Double_t*   fMesonWidthGaussianError                                    = nullptr;
    Double_t*   fMesonTrueWidthGaussian                                     = nullptr;
    Double_t*   fMesonTrueWidthGaussianError                                = nullptr;
    TH1D*       fHistoMassGaussianMeson                                     = nullptr;
    TH1D*       fHistoTrueMassGaussianMeson                                 = nullptr;
    TH1D*       fHistoWidthGaussianMeson                                    = nullptr;
    TH1D*       fHistoTrueWidthGaussianMeson                                = nullptr;

    //****************************************************************************
    //******************** Peak position and calibration *************************
    //****************************************************************************
    TF1**       fFitSignalPeakPosInvMassLeftPtBin                           = nullptr;
    TF1**       fFitSignalPeakPosInvMassPtBin                               = nullptr;

    //****************************************************************************
    //**************** Decomposition of signal for calo clusters *****************
    //****************************************************************************
    TF1**       fFitTrueSignalCaloPhotonInvMassPtBin                        = nullptr;
    TF1**       fFitTrueSignalCaloConvPhotonInvMassPtBin                    = nullptr;
    TF1**       fFitTrueSignalCaloElectronInvMassPtBin                      = nullptr;
    TF1**       fFitTrueSignalMixedCaloConvPhotonInvMassPtBin               = nullptr;
    TF1**       fFitTrueSignalCaloMergedClusterInvMassPtBin                 = nullptr;
    TF1**       fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin         = nullptr;
    Double_t*   fMesonTrueMassCaloPhoton                                    = nullptr;
    Double_t*   fMesonTrueMassCaloElectron                                  = nullptr;
    Double_t*   fMesonTrueMassCaloConvPhoton                                = nullptr;
    Double_t*   fMesonTrueMassCaloMergedCluster                             = nullptr;
    Double_t*   fMesonTrueMassCaloMergedClusterPartConv                     = nullptr;
    Double_t*   fMesonTrueMassMixedCaloConvPhoton                           = nullptr;
    Double_t*   fMesonTrueFWHMCaloPhoton                                    = nullptr;
    Double_t*   fMesonTrueFWHMCaloElectron                                  = nullptr;
    Double_t*   fMesonTrueFWHMCaloConvPhoton                                = nullptr;
    Double_t*   fMesonTrueFWHMCaloMergedCluster                             = nullptr;
    Double_t*   fMesonTrueFWHMCaloMergedClusterPartConv                     = nullptr;
    Double_t*   fMesonTrueFWHMMixedCaloConvPhoton                           = nullptr;
    Double_t*   fMesonTrueMassErrorCaloPhoton                               = nullptr;
    Double_t*   fMesonTrueMassErrorCaloElectron                             = nullptr;
    Double_t*   fMesonTrueMassErrorCaloConvPhoton                           = nullptr;
    Double_t*   fMesonTrueMassErrorCaloMergedCluster                        = nullptr;
    Double_t*   fMesonTrueMassErrorCaloMergedClusterPartConv                = nullptr;
    Double_t*   fMesonTrueMassErrorMixedCaloConvPhoton                      = nullptr;
    Double_t*   fMesonTrueFWHMErrorCaloPhoton                               = nullptr;
    Double_t*   fMesonTrueFWHMErrorCaloElectron                             = nullptr;
    Double_t*   fMesonTrueFWHMErrorCaloConvPhoton                           = nullptr;
    Double_t*   fMesonTrueFWHMErrorCaloMergedCluster                        = nullptr;
    Double_t*   fMesonTrueFWHMErrorCaloMergedClusterPartConv                = nullptr;
    Double_t*   fMesonTrueFWHMErrorMixedCaloConvPhoton                      = nullptr;
    TH1D*       fHistoTrueMassMesonCaloPhoton                               = nullptr;
    TH1D*       fHistoTrueMassMesonCaloElectron                             = nullptr;
    TH1D*       fHistoTrueMassMesonCaloConvPhoton                           = nullptr;
    TH1D*       fHistoTrueMassMesonCaloMergedCluster                        = nullptr;
    TH1D*       fHistoTrueMassMesonCaloMergedPartConvCluster                = nullptr;
    TH1D*       fHistoTrueMassMesonMixedCaloConvPhoton                      = nullptr;
    TH1D*       fHistoTrueFWHMMesonCaloPhoton                               = nullptr;
    TH1D*       fHistoTrueFWHMMesonCaloElectron                             = nullptr;
    TH1D*       fHistoTrueFWHMMesonCaloConvPhoton                           = nullptr;
    TH1D*       fHistoTrueFWHMMesonCaloMergedCluster                        = nullptr;
    TH1D*       fHistoTrueFWHMMesonCaloMergedPartConvCluster                = nullptr;
    TH1D*       fHistoTrueFWHMMesonMixedCaloConvPhoton                      = nullptr;

    //****************************************************************************
    //****************** significance & S/B doubles, histograms ******************
    //****************************************************************************
    Double_t*   fMesonSBdefault[3]                                          = { nullptr, nullptr, nullptr};
    Double_t*   fMesonSigndefault[3]                                        = { nullptr, nullptr, nullptr};
    Double_t*   fMesonSBdefaultError[3]                                     = { nullptr, nullptr, nullptr};
    Double_t*   fMesonSigndefaultError[3]                                   = { nullptr, nullptr, nullptr};
    Double_t*   fMesonTrueSB                                                = nullptr;
    Double_t*   fMesonTrueSign                                              = nullptr;
    Double_t*   fMesonTrueSBError                                           = nullptr;
    Double_t*   fMesonTrueSignError                                         = nullptr;
    TH1D*       fHistoSBdefaultMeson[3]                                     = { nullptr, nullptr, nullptr};
    TH1D*       fHistoSigndefaultMeson[3]                                   = { nullptr, nullptr, nullptr};
    TH1D*       fHistoTrueSignMeson                                         = nullptr;
    TH1D*       fHistoTrueSBMeson                                           = nullptr;
    TH1D*       fHistoLambdaTail                                            = nullptr;
    TH1D*       fHistoTrueLambdaTail                                        = nullptr;
    TH1D*       fHistoAmplitude                                             = nullptr;
    TH1D*       fHistoSigma                                                 = nullptr;
    TH1D*       fHistoTrueSigma                                             = nullptr;
    TH1D*       fHistoResidualBGlin                                         = nullptr;
    TH1D*       fHistoResidualBGcon                                         = nullptr;
    TH1D*       fHistoChi2[4]                                               = { nullptr, nullptr, nullptr, nullptr };
    TH1D*       fHistoChi2SigToBckFit[2]                                    = { nullptr, nullptr };
    TH1D*       fHistoRatioResBGYield                                       = nullptr;
    TH1D*       fHistoRatioResBGYieldToSPlusResBG                           = nullptr;
    TH1D*       fHistoResBGYield[4]                                         = { nullptr, nullptr, nullptr, nullptr };

    Double_t*   fMassWindowHigh[3]                                          = { nullptr, nullptr, nullptr};
    Double_t*   fMassWindowLow[3]                                           = { nullptr, nullptr, nullptr};
    TH1D*       fHistoMassWindowHigh[3]                                     = { nullptr, nullptr, nullptr};
    TH1D*       fHistoMassWindowLow[3]                                      = { nullptr, nullptr, nullptr};

    //****************************************************************************
    //**************************** MC input histograms ***************************
    //****************************************************************************
    TH1D*       fHistoMCMesonPtWithinAcceptance                             = nullptr;
    TH1D*       fHistoMCMesonPtWithinAcceptanceWOWeights                    = nullptr;
    TH1D*       fHistoMCMesonPtWithinAcceptanceWOEvtWeights                 = nullptr;
    TH1D*       fHistoMCMesonPt                                             = nullptr;
    TH1D*       fHistoMCMesonPtWOWeights                                    = nullptr;
    TH1D*       fHistoMCMesonPtWOEvtWeights                                 = nullptr;
    TH1D*       fHistoMCMesonPtWeights                                      = nullptr;
    TH1D*       fHistoMCMesonWithinAccepPt                                  = nullptr;
    TH1D*       fHistoMCMesonWithinAccepPtWOWeights                         = nullptr;
    TH1D*       fHistoMCMesonWithinAccepPtWOEvtWeights                      = nullptr;
    TH1D*       fHistoMCMesonPt1                                            = nullptr;
    TH1D*       fHistoMCMesonPt1WOWeights                                   = nullptr;
    TH1D*       fHistoMCMesonPt1WOEvtWeights                                = nullptr;
    //****************************************************************************
    //************************ sec MC input histograms ***************************
    //****************************************************************************
    TH2D*       fHistoMCSecPi0SourcePt                                      = nullptr;
    TH2D*       fHistoMCSecPi0WAccSourcePt                                  = nullptr;
    TH1D*       fHistoMCSecPi0Pt[4]                                         = { nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoMCSecPi0PtWAcc[4]                                     = { nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoMCSecPi0PtReb[4]                                      = { nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoMCSecPi0PtWAccReb[4]                                  = { nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoMCSecPi0AcceptPt[4]                                   = { nullptr, nullptr, nullptr, nullptr};

    //****************************************************************************
    //*********************** MC efficiency histograms ***************************
    //****************************************************************************
    TH1D*       fHistoMCMesonAcceptPt                                       = nullptr;
    TH1D*       fHistoMCMesonAcceptPtWOWeights                              = nullptr;
    TH1D*       fHistoMCMesonAcceptPtWOEvtWeights                           = nullptr;

    TH1D*       fHistoMonteMesonEffiPt[6]                                   = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoMCTrueMesonEffiPt[3]                                  = {nullptr, nullptr, nullptr};
    TH1D*       fHistoMCTrueMesonEffiPtReweighted[3]                        = {nullptr, nullptr, nullptr};
    TH1D*       fHistoMCTrueMesonEffiPtUnweighted[3]                        = {nullptr, nullptr, nullptr};
    TH1D*       fHistoMCTrueSecMesonEffiPt[3][4]                            = {{nullptr, nullptr, nullptr, nullptr},{nullptr, nullptr, nullptr, nullptr},{nullptr, nullptr, nullptr, nullptr}};

    //****************************************************************************
    //**************************** MC rec sec mesons  ****************************
    //****************************************************************************
    Double_t*   fMesonTrueSecYields[3][4]                                   = {{nullptr, nullptr, nullptr, nullptr},{nullptr, nullptr, nullptr, nullptr},{nullptr, nullptr, nullptr, nullptr}};
    Double_t*   fMesonTrueSecYieldsError[3][4]                              = {{nullptr, nullptr, nullptr, nullptr},{nullptr, nullptr, nullptr, nullptr},{nullptr, nullptr, nullptr, nullptr}};
    TH1D*       fHistoYieldTrueSecMeson[3][4]                               = {{nullptr, nullptr, nullptr, nullptr},{nullptr, nullptr, nullptr, nullptr},{nullptr, nullptr, nullptr, nullptr}};
    TH1D*       fHistoYieldTrueSecFracMeson[3][4]                           = {{nullptr, nullptr, nullptr, nullptr},{nullptr, nullptr, nullptr, nullptr},{nullptr, nullptr, nullptr, nullptr}};

    //****************************************************************************
    //******************* yield extraction fixed windows val. MC *****************
    //****************************************************************************
    Double_t*   fMesonTrueYieldFixedWindow                                  = nullptr;
    Double_t*   fMesonTrueYieldGammaFixedWindow                             = nullptr;
    Double_t*   fMesonTrueYieldGammaConvGammaFixedWindow                    = nullptr;
    Double_t*   fMesonTrueYieldConvGammaConvGammaFixedWindow                = nullptr;
    Double_t*   fMesonTrueYieldErrorFixedWindow                             = nullptr;
    Double_t*   fMesonTrueYieldGammaErrorFixedWindow                        = nullptr;
    Double_t*   fMesonTrueYieldGammaConvGammaErrorFixedWindow               = nullptr;
    Double_t*   fMesonTrueYieldConvGammaConvGammaErrorFixedWindow           = nullptr;
    TH1D*       fHistoYieldTrueMesonFixedWindow                             = nullptr;
    TH1D*       fHistoYieldTrueMesonGammaFixedWindow                        = nullptr;
    TH1D*       fHistoYieldTrueMesonGammaConvGammaFixedWindow               = nullptr;
    TH1D*       fHistoYieldTrueMesonConvGammaConvGammaFixedWindow           = nullptr;

    //*****************************************************************************
    //*********************** dedicated cluster histograms ************************
    //*****************************************************************************
    TH1D*       fHistoClustersPt                                            = nullptr;
    TH1D*       fHistoClustersE                                             = nullptr;
    TH1D*       fHistoClustersOverlapHeadersPt                              = nullptr;

    //*****************************************************************************
    //******** Check multiple counts of neutral mesons, clusters ******************
    //*****************************************************************************
    TH2D*       fHistoTrueMesonDCInvMassVSPt                                = nullptr;
    TH1D**      fHistoMappingTrueMesonDCInvMassPtBins                       = nullptr;
    TH1F*       fHistoTrueMesonMultipleCount                                = nullptr;
    Double_t*   fMesonTrueYieldsDC                                          = nullptr;
    Double_t*   fMesonTrueYieldsDCError                                     = nullptr;
    TH1D*       fHistoYieldTrueMesonDC                                      = nullptr;
    TH1F*       fHistoTrueGammaClusPt                                       = nullptr;
    TH1D*       fHistoTrueGammaDCClusPt                                     = nullptr;
    TH1F*       fHistoTrueGammaClusMultipleCount                            = nullptr;

    //*****************************************************************************
    //************ Load secondary pion histograms from external file **************
    //*****************************************************************************
    Bool_t      fHaveToyMCInputForSec                                       = kFALSE;
    Bool_t      fHaveCocktailInputForSec                                    = kFALSE;
    TFile*      fFileToyMCInput[3]                                          = {nullptr, nullptr, nullptr};
    TFile*      fFileCocktailInput                                          = nullptr;
    TH1D*       fHistoYieldExternSecInput[3]                                = {nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldExternSecInputReb[3]                             = {nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldExternResonanceFeedDownInput[15]                 = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                                                                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    TH1D*       fHistoYieldExternResonanceFeedDownInputReb[15]              = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                                                                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

    //*****************************************************************************
    //******** monitoring of SPD pileup histograms ********************************
    //*****************************************************************************
    TH1F*       fHistoPileUpVertexDistance                                  = nullptr;
    TH1F*       fHistoPileUpVertexDistance_SPDPileup                        = nullptr;
    TH1F*       fHistoPileUpVertexDistance_TrackletHits                     = nullptr;

    //****************************************************************************************************
    //*************** Auxiliary function for setting a range *********************************************
    //****************************************************************************************************
    template<class TYPE_t>
    void SetRangeArray(TYPE_t* array, TYPE_t from, TYPE_t to){
        array[0] = from;
        array[1] = to;
    }


    //****************************************************************************************************
    //*************** Function to initalize different fitting, plotting and integration windows***********
    //****************************************************************************************************
    void InitializeWindows(TString setPi0, Int_t mode, TString trigger, Int_t triggerSet = -1){

        // Heavy meson analysis
        if(mode>=100) mode -= 100;

        fPeakRange                  = new Double_t[2];
        fIntFixedRange              = new Double_t[2];
        fFitRange                   = new Double_t[2];
        fBGFitRange                 = new Double_t[2];
        fBGFitRangeLeft             = new Double_t[2];
        fMesonPlotRange             = new Double_t[2];
        fMesonIntDeltaRange         = new Double_t[2];
        fMesonIntDeltaRangeWide     = new Double_t[2];
        fMesonIntDeltaRangeNarrow   = new Double_t[2];
        fMesonMassRange             = new Double_t[2];
        fMesonMassPlotRange         = new Double_t[2];
        fMesonFitRange              = new Double_t[2];
        fMesonWidthRange            = new Double_t[2];
        fMesonLambdaTailRange       = new Double_t[2];
        fMesonLambdaTailRangeNominal= new Double_t[2];
        fMesonLambdaTailRangeMC     = new Double_t[2];
        fMesonWidthRangeMC          = new Double_t[2];
        fMesonLambdaTailRangeMC[0]  = -10;
        fMesonLambdaTailRangeMC[1]  = -10;
        fMesonWidthRangeMC[0]       = -10;
        fMesonWidthRangeMC[1]       = -10;
        fMesonLambdaTailMC          = -10;
        fMidPt                      = new Double_t[2];

        fFullPt                     = new Double_t[2];
        fFullPt[0]                  = 0.4;
        fFullPt[1]                  = 15;

        //****************************************************************************************************
        // Initialization for pi0 meson
        //****************************************************************************************************
        if (setPi0.Contains("Pi0")){

            // set meson ID according to PDG
            fMesonId                    = 111;

            // set medium pt range
            fMidPt[0]                   = 0.8;
            fMidPt[1]                   = 2.5;
            if(mode == 0 && fEnergyFlag.CompareTo("PbPb_5.02TeV") == 0){  // just temporary (18.1.2017)
            fMidPt[0]                   = 4.0;                           // due to strange BG shape and wrong fits at lower pT
            fMidPt[1]                   = 6.0;
            }

            // Initialize peak range
            if (mode == 2 || mode == 13){                           // PCM-EMC
                fPeakRange[0]               = 0.05;
                fPeakRange[1]               = 0.145;
            } else if ( mode == 4 || mode == 12) {                  // EMC
                fPeakRange[0]               = 0.05;
                fPeakRange[1]               = 0.17;
                if ( fEnergyFlag.CompareTo("7TeV") == 0 && ( trigger.CompareTo("52") == 0 || triggerSet == 2))
                    fPeakRange[1] = 0.19;
                if ( fEnergyFlag.BeginsWith("8TeV") && (trigger.CompareTo("52") == 0 || triggerSet == 1 || trigger.CompareTo("81") == 0 || triggerSet == 2))
                    fPeakRange[1] = 0.19;
            } else if ( mode == 5) {                                // PHOS
                fPeakRange[0]               = 0.05;
                fPeakRange[1]               = 0.145;
            } else if ( fMode == 0 && fEnergyFlag.CompareTo("13TeVLowB") == 0) {
                fPeakRange[0]               = 0.05;
                fPeakRange[1]               = 0.145;
            } else {                                                // defaults
                fPeakRange[0]               = 0.1;
                fPeakRange[1]               = 0.145;
            }
            // Initialize fit range
            fIntFixedRange[0]           = 0.08;
            fIntFixedRange[1]           = 0.2;
            if (mode == 2 || mode == 13){                           // PCM-EMC
                fFitRange[0]                = 0.02;
                fFitRange[1]                = 0.25;
            } else if (mode == 4 || mode == 12 ) {                  // EMC
                fFitRange[0]                = 0.07;
                fFitRange[1]                = 0.25;
            } else if ( mode == 5) {                                // PHOS
                fFitRange[0]                = 0.07;
                fFitRange[1]                = 0.25;
            } else {                                                // defaults
                fFitRange[0]                = 0.05;
                fFitRange[1]                = 0.25;
            }

            // Initialize default BG fit range right & left
            if (mode == 2 || mode == 13){                           // PCM-EMC
                fBGFitRange[0]              = 0.19;
                fBGFitRange[1]              = 0.3;
                fBGFitRangeLeft[0]          = 0.03;
                fBGFitRangeLeft[1]          = 0.05;
                if ( fEnergyFlag.CompareTo("13TeV") == 0 ){
                    if ( trigger.CompareTo("8e")==0 ){
                        fBGFitRange[0] = 0.25;
                        fBGFitRange[1] = 0.30;
                    } else if ( trigger.CompareTo("8d")==0 ){
                        fBGFitRange[0] = 0.25;
                        fBGFitRange[1] = 0.30;
                    }
                }
            } else if (mode == 4 || mode == 12 ) {                  // EMC
                fBGFitRange[0]              = 0.19;
                fBGFitRange[1]              = 0.3;
                fBGFitRangeLeft[0]          = 0.05;
                fBGFitRangeLeft[1]          = 0.08;
                if ( fEnergyFlag.CompareTo("7TeV") == 0 ){
                    if( trigger.CompareTo("52") == 0 || triggerSet == 1 ){
                        fBGFitRange[0] = 0.26;
                    }
                }
                if ( fEnergyFlag.BeginsWith("8TeV") ){
                    if( trigger.CompareTo("52") == 0 || triggerSet == 1 ){
                        fBGFitRange[0] = 0.25;
                    } else if ( trigger.CompareTo("81")==0 || triggerSet == 2 ){
                        fBGFitRange[0] = 0.26;
                    }
                }
                if ( fEnergyFlag.Contains("pPb_8TeV") ){
                    if( trigger.CompareTo("8e") == 0 || triggerSet == 1 ){
                        fBGFitRange[0] = 0.25;
                    } else if ( trigger.CompareTo("8d")==0 || triggerSet == 2 ){
                        fBGFitRange[0] = 0.26;
                    }
                }
                if ( fEnergyFlag.CompareTo("13TeV") == 0 ){
                    if( trigger.CompareTo("52") == 0 || triggerSet == 1 ){
                        fBGFitRange[0] = 0.25;
                    } else if ( trigger.CompareTo("85")==0 || triggerSet == 3 ){
                        fBGFitRange[0] = 0.26;
                        fBGFitRange[1] = 0.31;
                    } else if ( trigger.CompareTo("83")==0 || triggerSet == 2 ){
                        fBGFitRange[0] = 0.28;
                        fBGFitRange[1] = 0.34;
                    }
                }
                if ( fEnergyFlag.CompareTo("5TeV2017") == 0 ){
                  if( trigger.CompareTo("a1") == 0 || trigger.CompareTo("a2") == 0){
                    fBGFitRange[0]              = 0.23;
                    fBGFitRange[1]              = 0.3;
                  }
                }
                if ( fEnergyFlag.CompareTo("PbPb_5.02TeV") == 0 ){
                  fBGFitRange[0]              = 0.24;
                  fBGFitRange[1]              = 0.3;
                }
            } else if ( mode == 5) {                                // PHOS
                fBGFitRange[0]              = 0.19;
                fBGFitRange[1]              = 0.3;
                fBGFitRangeLeft[0]          = 0.04;
                fBGFitRangeLeft[1]          = 0.08;
            } else {                                                // defaults
                fBGFitRange[0]              = 0.17;
                fBGFitRange[1]              = 0.3;
                fBGFitRangeLeft[0]          = 0.05;
                fBGFitRangeLeft[1]          = 0.08;
            }

            // Initialize default Plot range for meson
            fMesonPlotRange[0]          = 0.13;
            fMesonPlotRange[1]          = 0.138;

            // Initialize default Plot default integration ranges
            if (mode == 0){                                                         // PCM
                if ( fEnergyFlag.BeginsWith("8TeV")){
                    fMesonIntDeltaRange[0]      = -0.035;
                    fMesonIntDeltaRange[1]      =  0.010;
                } else if ( !fEnergyFlag.CompareTo("pPb_8TeV")){
                    fMesonIntDeltaRange[0]      = -0.030;
                    fMesonIntDeltaRange[1]      =  0.015;
                } else {
                    fMesonIntDeltaRange[0]      = -0.035;
                    fMesonIntDeltaRange[1]      = 0.012;
                }
                fMesonIntDeltaRangeWide[0]      = -0.055;
                fMesonIntDeltaRangeWide[1]      = 0.025;
                fMesonIntDeltaRangeNarrow[0]    = -0.015;
                fMesonIntDeltaRangeNarrow[1]    = 0.005;
            } else if (mode == 2 || mode == 13){                                    // PCM-EMC, PCM-DMC
                fMesonIntDeltaRange[0]      = -0.032;
                fMesonIntDeltaRange[1]      = 0.022;
                fMesonIntDeltaRangeWide[0]  = -0.048;
                fMesonIntDeltaRangeWide[1]  = 0.028;
                fMesonIntDeltaRangeNarrow[0]= -0.016;
                fMesonIntDeltaRangeNarrow[1]= 0.016;
                if( fEnergyFlag.Contains("5TeV2017") ){
                    fMesonIntDeltaRange[0]      = fMesonIntDeltaRange[0]*1.3;
                    fMesonIntDeltaRange[1]      = fMesonIntDeltaRange[1]*1.3;
                    fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRangeWide[0]*1.3;
                    fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRangeWide[1]*1.3;
                    fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRangeNarrow[1]*1.3;
                    fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRangeNarrow[0]*1.3;
                }
                if( fEnergyFlag.CompareTo("7TeV") == 0 ){
                    if( trigger.CompareTo("52") == 0 || triggerSet == 1){
                        fMesonIntDeltaRange[0]      = fMesonIntDeltaRange[0]*1.3;
                        fMesonIntDeltaRange[1]      = fMesonIntDeltaRange[1]*1.3;
                        fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRangeWide[0]*1.3;
                        fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRangeWide[1]*1.3;
                        fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRangeNarrow[1]*1.3;
                        fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRangeNarrow[0]*1.3;
                    }
                }
                if( fEnergyFlag.BeginsWith("8TeV") ){
                    if( trigger.CompareTo("81") == 0 || triggerSet == 2){
                        fMesonIntDeltaRange[0]      = fMesonIntDeltaRange[0]*1.3;
                        fMesonIntDeltaRange[1]      = fMesonIntDeltaRange[1]*1.3;
                        fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRangeWide[0]*1.3;
                        fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRangeWide[1]*1.3;
                        fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRangeNarrow[1]*1.3;
                        fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRangeNarrow[0]*1.3;
                    }else if( trigger.CompareTo("52") == 0 || triggerSet == 1){
                        fMesonIntDeltaRange[0]      = fMesonIntDeltaRange[0]*1.1;
                        fMesonIntDeltaRange[1]      = fMesonIntDeltaRange[1]*1.1;
                        fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRangeWide[0]*1.1;
                        fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRangeWide[1]*1.1;
                        fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRangeNarrow[1]*1.1;
                        fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRangeNarrow[0]*1.1;
                    }
                }
                if( fEnergyFlag.BeginsWith("pPb_8TeV") ){
                    fMesonIntDeltaRange[0]      = -0.032;
                    fMesonIntDeltaRange[1]      = 0.022;
                    fMesonIntDeltaRangeWide[0]  = -0.048;
                    fMesonIntDeltaRangeWide[1]  = 0.028;
                    fMesonIntDeltaRangeNarrow[0]= -0.026;
                    fMesonIntDeltaRangeNarrow[1]= 0.016;
                    if( trigger.CompareTo("8d") == 0 || triggerSet == 2){
                        fMesonIntDeltaRange[0]      = fMesonIntDeltaRange[0]*1.3;
                        fMesonIntDeltaRange[1]      = fMesonIntDeltaRange[1]*1.5;
                        fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRangeWide[0]*1.3;
                        fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRangeWide[1]*1.5;
                        fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRangeNarrow[0]*1.3;
                        fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRangeNarrow[1]*1.5;
                    }else if( trigger.CompareTo("8e") == 0 || triggerSet == 1){
                        fMesonIntDeltaRange[0]      = fMesonIntDeltaRange[0]*1.1;
                        fMesonIntDeltaRange[1]      = fMesonIntDeltaRange[1]*1.3;
                        fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRangeWide[0]*1.1;
                        fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRangeWide[1]*1.3;
                        fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRangeNarrow[0]*1.1;
                        fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRangeNarrow[1]*1.3;
                    }
                }
            } else if (mode == 3){                                                  // PCM-PHOS
                fMesonIntDeltaRange[0]      = -0.038;
                fMesonIntDeltaRange[1]      =  0.018;
                fMesonIntDeltaRangeWide[0]  = -0.055;
                fMesonIntDeltaRangeWide[1]  =  0.028;
                fMesonIntDeltaRangeNarrow[0]= -0.015;
                fMesonIntDeltaRangeNarrow[1]=  0.008;
                if( fEnergyFlag.CompareTo("PbPb_5.02TeV") == 0 ){
                  fMesonIntDeltaRange[0]      = -0.030;
                  fMesonIntDeltaRange[1]      =  0.025;
                  fMesonIntDeltaRangeWide[0]  = -0.040;
                  fMesonIntDeltaRangeWide[1]  =  0.035;
                  fMesonIntDeltaRangeNarrow[0]= -0.020;
                  fMesonIntDeltaRangeNarrow[1]=  0.015;
                }
            } else if (mode == 4 || mode == 12 ) {                                  // EMC, DMC
                fMesonIntDeltaRange[0]      = -0.05;
                fMesonIntDeltaRange[1]      = 0.04;
                fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRange[0]*1.4;
                fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRange[1]*1.4;
                fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRange[0]*0.6;
                fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRange[1]*0.6;
                if(fEnergyFlag.CompareTo("7TeV") == 0){
                    if(trigger.CompareTo("52")==0 || triggerSet == 1){
                        fMesonIntDeltaRange[0]      = -0.06;
                        fMesonIntDeltaRange[1]      = 0.08;
                        fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRange[0]*1.2;
                        fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRange[1]*1.2;
                        fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRange[0]*0.8;
                        fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRange[1]*0.8;
                    }
                }
                if(fEnergyFlag.BeginsWith("8TeV")){
                    if(trigger.CompareTo("52")==0 || triggerSet == 1){
                        fMesonIntDeltaRange[0]      = -0.05;
                        fMesonIntDeltaRange[1]      = 0.06;
                        fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRange[0]*1.4;
                        fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRange[1]*1.4;
                        fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRange[0]*0.6;
                        fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRange[1]*0.6;
                    }else if(trigger.CompareTo("81")==0 || triggerSet == 2 ){
                        fMesonIntDeltaRange[0]      = -0.06;
                        fMesonIntDeltaRange[1]      = 0.08;
                        fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRange[0]*1.2;
                        fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRange[1]*1.2;
                        fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRange[0]*0.8;
                        fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRange[1]*0.8;
                    }
                }
                if(fEnergyFlag.Contains("pPb_8TeV")){
                    if(trigger.CompareTo("8e")==0 || triggerSet == 2){
                        fMesonIntDeltaRange[0]      = -0.05;
                        fMesonIntDeltaRange[1]      = 0.06;
                        fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRange[0]*1.4;
                        fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRange[1]*1.4;
                        fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRange[0]*0.6;
                        fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRange[1]*0.6;
                    }else if(trigger.CompareTo("8d")==0 || triggerSet == 3 ){
                        fMesonIntDeltaRange[0]      = -0.06;
                        fMesonIntDeltaRange[1]      = 0.08;
                        fMesonIntDeltaRangeWide[0]  = fMesonIntDeltaRange[0]*1.2;
                        fMesonIntDeltaRangeWide[1]  = fMesonIntDeltaRange[1]*1.2;
                        fMesonIntDeltaRangeNarrow[0]= fMesonIntDeltaRange[0]*0.8;
                        fMesonIntDeltaRangeNarrow[1]= fMesonIntDeltaRange[1]*0.8;
                    }
                }
            } else if (mode == 5) {                                                 // PHOS
                fMesonIntDeltaRange[0]      = -0.035;
                fMesonIntDeltaRange[1]      = 0.023;
                fMesonIntDeltaRangeWide[0]  = -0.055;
                fMesonIntDeltaRangeWide[1]  = 0.028;
                fMesonIntDeltaRangeNarrow[0]= -0.018;
                fMesonIntDeltaRangeNarrow[1]= 0.015;
            } else {                                                                // defaults
                fMesonIntDeltaRange[0]      = -0.035;
                fMesonIntDeltaRange[1]      = 0.010;
                fMesonIntDeltaRangeWide[0]  = -0.055;
                fMesonIntDeltaRangeWide[1]  = 0.025;
                fMesonIntDeltaRangeNarrow[0]= -0.015;
                fMesonIntDeltaRangeNarrow[1]= 0.005;
            }
            // Set meson mass ranges for fitting and plotting
            fMesonMassPlotRange[0]      = 0.;
            fMesonMassPlotRange[1]      = 0.3;
            fMesonMassRange[0]          = 0.;
            fMesonMassRange[1]          = 0.3;
            if (mode == 0 && fEnergyFlag.CompareTo("pPb_5.023TeV") == 0 )
                fMesonMassRange[1]      = 0.14;

            // Set Meson fit range
            if (mode == 0){
                if( fEnergyFlag.Contains("PbPb")){                      // PCM
                    fMesonFitRange[0]       = 0.07; //0.07 -> 0.05
                    fMesonFitRange[1]       = 0.25; //0.22 -> 0.25
                } else if( fEnergyFlag.Contains("pPb_5.023TeV")){
                    fMesonFitRange[0]       = 0.03;
                    fMesonFitRange[1]       = 0.25;
               } else if( fEnergyFlag.Contains("13TeV")){
                    fMesonFitRange[0]       = 0.05;
                    fMesonFitRange[1]       = 0.25;
                } else {
                    fMesonFitRange[0]       = 0.03;
                    fMesonFitRange[1]       = 0.25;
                }
            } else if (mode == 3){                        // PCM-PHOS
                fMesonFitRange[0]       = 0.01;
                fMesonFitRange[1]       = 0.3;
                if (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0){
                    fMesonFitRange[0]           = 0.05;
                    fMesonFitRange[1]           = 0.25;
                }
            } else if (mode == 2 || mode == 13){                        // PCM-EMC
                fMesonFitRange[0]       = 0.05;
                fMesonFitRange[1]       = 0.25;
                if( fEnergyFlag.CompareTo("7TeV") == 0){
                    if( trigger.CompareTo("52") == 0 ){
                        fMesonFitRange[0]       = 0.05;
                        fMesonFitRange[1]       = 0.27;
                        if (setPi0.CompareTo("Pi0EtaBinning") == 0){
                            fMesonFitRange[1]   = 0.29;
                        }
                    }
                }else if( fEnergyFlag.BeginsWith("8TeV")){
                    if( trigger.CompareTo("81") == 0 ){
                        fMesonFitRange[0]       = 0.05;
                        fMesonFitRange[1]       = 0.27;
                    } else if( trigger.CompareTo("52") == 0 ){
                        fMesonFitRange[0]       = 0.03;
                        fMesonFitRange[1]       = 0.25;
                        if (setPi0.CompareTo("Pi0EtaBinning") == 0){
                            fMesonFitRange[1]   = 0.29;
                        }
                    }
                } else if (fEnergyFlag.Contains("pPb_5.023TeV") ){
                    fMesonFitRange[0]           = 0.03;
                    fMesonFitRange[1]           = 0.25;
                } else if (fEnergyFlag.CompareTo("XeXe_5.44TeV") == 0 ){
                    fMesonFitRange[0]           = 0.03;
                    fMesonFitRange[1]           = 0.25;
                }
            } else if (mode == 4 || mode == 12 ) {                      // EMC
                fMesonFitRange[0]               = 0.06;
                fMesonFitRange[1]               = 0.25;
                if( fEnergyFlag.CompareTo("7TeV") == 0 ){
                    if (setPi0.CompareTo("Pi0EtaBinning") == 0){
                        fMesonFitRange[0]               = 0.06;
                        fMesonFitRange[1]               = 0.30;
                        if( trigger.CompareTo("52") == 0){
                            fMesonFitRange[0] = 0.08;
                            fMesonFitRange[1] = 0.28;
                        }
                    } else if( trigger.CompareTo("52") == 0 ){
                    fMesonFitRange[0] = 0.08;
                    fMesonFitRange[1] = 0.29;
                    } else {
                    fMesonFitRange[0] = 0.08;
                    fMesonFitRange[1] = 0.25;
                    }
                } else if( fEnergyFlag.BeginsWith("8TeV") ){
                    if (setPi0.CompareTo("Pi0EtaBinning") == 0){
                        fMesonFitRange[0]               = 0.06;
                        fMesonFitRange[1]               = 0.30;
                        if( trigger.CompareTo("81") == 0 || triggerSet == 2){
                            fMesonFitRange[0] = 0.08;
                            fMesonFitRange[1] = 0.28;
                        } else if( trigger.CompareTo("52") == 0 ){
                            fMesonFitRange[0] = 0.06;
                            fMesonFitRange[1] = 0.27;
                        }
                    } else if ( trigger.CompareTo("81") == 0 || triggerSet == 2){
                    fMesonFitRange[0] = 0.08;
                    fMesonFitRange[1] = 0.29;
                    } else if( trigger.CompareTo("52") == 0 ){
                    fMesonFitRange[0] = 0.04;
                    fMesonFitRange[1] = 0.26;
                    } else {
                    fMesonFitRange[0] = 0.08;
                    fMesonFitRange[1] = 0.25;
                    }
                } else if( fEnergyFlag.Contains("pPb_8TeV") ){
                    if (setPi0.CompareTo("Pi0EtaBinning") == 0){
                        fMesonFitRange[0]               = 0.06;
                        fMesonFitRange[1]               = 0.30;
                        if( trigger.CompareTo("8d") == 0 || triggerSet == 3){
                            fMesonFitRange[0] = 0.08;
                            fMesonFitRange[1] = 0.28;
                        } else if( trigger.CompareTo("8e") == 0 ){
                            fMesonFitRange[0] = 0.06;
                            fMesonFitRange[1] = 0.27;
                        }
                    } else if ( trigger.CompareTo("8d") == 0 || triggerSet == 2){
                        fMesonFitRange[0] = 0.08;
                        fMesonFitRange[1] = 0.30;
                    } else if( trigger.CompareTo("8e") == 0 ){
                        fMesonFitRange[0] = 0.065;
                        fMesonFitRange[1] = 0.30;
                    } else {
                        fMesonFitRange[0] = 0.055;
                        fMesonFitRange[1] = 0.28;
                    }
                } else if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 ){
                    fMesonFitRange[0]       = 0.06;
                    fMesonFitRange[1]       = 0.30;
                } else if( fEnergyFlag.CompareTo("PbPb_5.02TeV") == 0 ){
                    fMesonFitRange[0]       = 0.075;
                    fMesonFitRange[1]       = 0.25;
                }
            } else if (mode == 5){                                      // PHOS
                if (fEnergyFlag.CompareTo("XeXe_5.44TeV") == 0 ){
                    fMesonFitRange[0]       = 0.03;
                    fMesonFitRange[1]       = 0.25;
                } else {
                    fMesonFitRange[0]       = 0.08;
                    fMesonFitRange[1]       = 0.25;
                }
            } else {                                                    // default
                fMesonFitRange[0]       = 0.05;
                fMesonFitRange[1]       = 0.25;
            }

            // Set remaining parameters for fitting
            if (mode == 2 || mode == 13){                               // PCM-EMC
                fMesonWidthExpect               = 0.005;
                fMesonWidthRange[0]             = 0.001;
                fMesonWidthRange[1]             = 0.025;
                fMesonLambdaTail                = 0.012;
                fMesonLambdaTailRange[0]        = 0.001;
                fMesonLambdaTailRange[1]        = 0.09;
                if( fEnergyFlag.CompareTo("7TeV") == 0){
                    if( trigger.CompareTo("52") == 0){
                        fMesonLambdaTail                = 0.017;
                        fMesonLambdaTailRange[0]        = 0.01;
                        fMesonWidthExpect               = 0.010;
                        fMesonWidthRange[1]             = 0.050;
                    }
                }else if( fEnergyFlag.BeginsWith("8TeV")){
                    if( trigger.CompareTo("81") == 0){
                        fMesonLambdaTail                = 0.017;
                        fMesonLambdaTailRange[0]        = 0.01;
                        fMesonWidthExpect               = 0.010;
                        fMesonWidthRange[1]             = 0.050;
                    } else if( trigger.CompareTo("52") == 0){
                        fMesonWidthExpect               = 0.010;
                        fMesonWidthRange[1]             = 0.040;
                    }
                }else if( fEnergyFlag.CompareTo("pPb_8TeV")==0){
                    if( trigger.CompareTo("8d") == 0){
                        fMesonLambdaTail                = 0.017;
                        fMesonLambdaTailRange[0]        = 0.01;
                        fMesonWidthExpect               = 0.020;
                        fMesonWidthRange[0]             = 0.0035;
                        fMesonWidthRange[1]             = 0.025;
                    } else if( trigger.CompareTo("8e") == 0){
                        fMesonWidthExpect               = 0.010;
                        fMesonWidthRange[1]             = 0.040;
                    }
                } else if( fEnergyFlag.Contains("pPb_5.023TeV") ){
                    fMesonLambdaTailRange[1]        = 0.025;
                    fMesonWidthRange[1]             = 0.032;
                } else if( fEnergyFlag.CompareTo("900GeV") == 0 ){
                    fMesonLambdaTail                = 0.01;
                    fMesonLambdaTailRange[0]        = 0.01;
                    fMesonLambdaTailRange[1]        = 0.01;
                } else if( fEnergyFlag.CompareTo("XeXe_5.44TeV") == 0 ){
                    fMesonLambdaTail                = 0.012;
                    fMesonLambdaTailRange[0]        = 0.012;
                    fMesonLambdaTailRange[1]        = 0.012;
                }
            } else if (mode == 3){
                fMesonWidthExpect           = 0.006;
                fMesonLambdaTail            = 0.012;
                fMesonWidthRange[0]         = 0.001;
                fMesonWidthRange[1]         = 0.015;
                fMesonLambdaTailRange[0]    = 0.001;
                fMesonLambdaTailRange[1]    = 0.02;
                if( fEnergyFlag.CompareTo("XeXe_5.44TeV") == 0 ){
                    fMesonLambdaTail            = 0.007;
                    fMesonLambdaTailRange[0]    = 0.007;
                    fMesonLambdaTailRange[1]    = 0.007;
                }
            } else if (mode == 4 || mode == 12 ) {                      // EMC
                fMesonWidthExpect               = 0.01;
                fMesonWidthRange[0]             = 0.006;
                fMesonWidthRange[1]             = 0.028;
                fMesonLambdaTail                = 0.02;
                fMesonLambdaTailRange[0]        = 0.01;
                fMesonLambdaTailRange[1]        = 0.03;
                if (fEnergyFlag.CompareTo("7TeV") == 0){
                    fMesonLambdaTailRange[0]    = 0.005;
                    if ( trigger.CompareTo("52") == 0 ){
                        fMesonWidthExpect       = 0.035;
                        fMesonWidthRange[0]     = 0.001;
                        fMesonWidthRange[1]     = 0.040;
                    }
                }
                if (fEnergyFlag.BeginsWith("8TeV")){
                    fMesonLambdaTailRange[0]    = 0.005;
                    if ( trigger.CompareTo("81") == 0 || triggerSet == 2){
                        fMesonWidthExpect       = 0.035;
                        fMesonWidthRange[0]     = 0.001;
                        fMesonWidthRange[1]     = 0.040;
                    } else if ( trigger.CompareTo("52") == 0 ){
                        fMesonLambdaTailRange[0]= 0.003;
                        fMesonWidthExpect       = 0.012;
                        fMesonWidthRange[0]     = 0.001;
                        fMesonWidthRange[1]     = 0.030;
                    }
                }
                if (fEnergyFlag.Contains("pPb_8TeV") ){
                    fMesonLambdaTailRange[0]    = 0.005;
                    if ( trigger.CompareTo("8d") == 0 || triggerSet == 3){
                        fMesonWidthExpect       = 0.035;
                        fMesonWidthRange[0]     = 0.001;
                        fMesonWidthRange[1]     = 0.040;
                    } else if ( trigger.CompareTo("8e") == 0 ){
                        fMesonLambdaTailRange[0]= 0.003;
                        fMesonWidthExpect       = 0.012;
                        fMesonWidthRange[0]     = 0.001;
                        fMesonWidthRange[1]     = 0.030;
                    }
                } else if (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0){
                    fMesonLambdaTail            = 0.011;
                    fMesonLambdaTailRange[0]    = 0.005;
                    fMesonLambdaTailRange[1]    = 0.03;
                    fMesonWidthRange[0]         = 0.003;
                    fMesonWidthRange[1]         = 0.056;
                } else if (fEnergyFlag.CompareTo("pPb_5.023TeVRun2") == 0){
                    fMesonLambdaTail            = 0.011;
                    fMesonLambdaTailRange[0]    = 0.005;
                    fMesonLambdaTailRange[1]    = 0.04;
                    fMesonWidthRange[0]         = 0.003;
                    fMesonWidthRange[1]         = 0.056;
                }else if( fEnergyFlag.CompareTo("900GeV") == 0 ){
                    fMesonLambdaTail            = 0.012;
                    fMesonLambdaTailRange[0]    = 0.012;
                    fMesonLambdaTailRange[1]    = 0.012;
                } else if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0){
                    fMesonLambdaTail            = 0.011;
                    fMesonLambdaTailRange[0]    = 0.011;
                    fMesonLambdaTailRange[1]    = 0.011;
                    fMesonWidthExpect           = 0.015;
                    fMesonWidthRange[0]         = 0.01;
                    fMesonWidthRange[1]         = 0.04;
                } else if( fEnergyFlag.CompareTo("XeXe_5.44TeV") == 0 ){
                    fMesonLambdaTail            = 0.012;
                    fMesonLambdaTailRange[0]    = 0.012;
                    fMesonLambdaTailRange[1]    = 0.012;
                    fMesonWidthExpect           = 0.015;
                    fMesonWidthRange[0]         = 0.01;
                    fMesonWidthRange[1]         = 0.04;
                }
            } else if (mode == 5){                                      // PHOS
                fMesonWidthExpect               = 0.005;
                fMesonLambdaTail                = 0.012;
                fMesonWidthRange[0]             = 0.001;
                fMesonWidthRange[1]             = 0.040;
                fMesonLambdaTailRange[0]        = 0.001;
                fMesonLambdaTailRange[1]        = 0.03;
                if( fEnergyFlag.CompareTo("XeXe_5.44TeV") == 0 ){
                    fMesonLambdaTail            = 0.007;
                    fMesonLambdaTailRange[0]    = 0.007;
                    fMesonLambdaTailRange[1]    = 0.007;
                }
            } else {                                                    // default
                fMesonWidthExpect           = 0.003;
                fMesonLambdaTail            = 0.012;
                fMesonWidthRange[0]         = 0.001;
                fMesonLambdaTailRange[0]    = 0.001;
                if  ( fEnergyFlag.Contains("LowB") ){
                    fMesonLambdaTailRange[1]     = 0.03;
                    fMesonWidthRange[1]         = 0.012;
                } else {
                    fMesonLambdaTailRange[1]    = 0.02;
                    fMesonWidthRange[1]         = 0.009;
                }
            }

        //****************************************************************************************************
        // Initialization for eta meson
        //****************************************************************************************************
        } else if (setPi0.CompareTo("Eta") == 0){

            // set meson ID according to PDG
            fMesonId                    = 221;

            // set medium pt range
            fMidPt[0]                   = 1.5;
            fMidPt[1]                   = 2.5;

            // Initialize peak range
            if(mode ==0){
                fPeakRange[0]           = 0.48;
                fPeakRange[1]           = 0.58;
            if (fEnergyFlag.Contains("13TeV")){
                fPeakRange[0]           = 0.5;
                fPeakRange[1]           = 0.57;
            }
            }else if (mode == 2 || mode == 13){
                fPeakRange[0]           = 0.48;
                fPeakRange[1]           = 0.58;
            } else if (mode == 3){
                fPeakRange[0]           = 0.48;
                fPeakRange[1]           = 0.58;
            } else if ( mode == 4 || mode == 12) {
                fPeakRange[0]           = 0.51;
                fPeakRange[1]           = 0.59;
            } else if ( mode == 5 ) {
                fPeakRange[0]           = 0.51;
                fPeakRange[1]           = 0.59;
            } else {
                fPeakRange[0]           = 0.48;
                fPeakRange[1]           = 0.58;
            }
            // Initialize fit range
            fIntFixedRange[0]           = 0.48;
            fIntFixedRange[1]           = 0.58;
            if (mode == 2 || mode == 13){
                fFitRange[0]                = 0.4;
                fFitRange[1]                = 0.7;
                if( fEnergyFlag.CompareTo("7TeV") == 0 )
                    fFitRange[1]                = 0.72;
                if( fEnergyFlag.BeginsWith("8TeV") )
                    fFitRange[1]                = 0.72;
            } else if (mode == 3){
                fFitRange[0]                = 0.38;
                fFitRange[1]                = 0.7;
            } else if (mode == 4 || mode == 12 ) {
                if( fEnergyFlag.CompareTo("7TeV") == 0 ){
                    fFitRange[0]                = 0.38;
                    fFitRange[1]                = 0.70;
                    if( trigger.CompareTo("52") == 0 || triggerSet == 1)
                        fFitRange[1]            = 0.74;
                } else if( fEnergyFlag.BeginsWith("8TeV") ){
                    fFitRange[0]                = 0.34;
                    if( trigger.CompareTo("52") == 0 || triggerSet == 1 || trigger.CompareTo("81") == 0  || triggerSet == 2)
                        fFitRange[1]                = 0.74;
                    else
                        fFitRange[1]                = 0.7;
                } else {
                    fFitRange[0]                = 0.38;
                    fFitRange[1]                = 0.73;
                }
            } else if (mode == 5) {
                fFitRange[0]                = 0.45;
                fFitRange[1]                = 0.65;
            } else if (mode == 0) {
                if( fEnergyFlag.Contains("13TeV")  ){
                    fFitRange[0]                = 0.44;
                    fFitRange[1]                = 0.65;
                }

             } else {
                fFitRange[0]            = 0.4;
                fFitRange[1]            = 0.65;
            }

            // Initialize default BG fit range right & left
            if (mode == 2 || mode == 13){
                fBGFitRange[0]                  = 0.65;
                fBGFitRange[1]                  = 0.75;
                fBGFitRangeLeft[0]              = 0.35;
                fBGFitRangeLeft[1]              = 0.42;
                if( fEnergyFlag.Contains("5TeV2017") || fEnergyFlag.CompareTo("7TeV") == 0 )
                    fBGFitRangeLeft[1] = 0.46;
                if( fEnergyFlag.BeginsWith("8TeV") )
                    fBGFitRangeLeft[1] = 0.46;
            } else if (mode == 4 || mode == 12 ) {
                if( fEnergyFlag.BeginsWith("8TeV") || fEnergyFlag.CompareTo("7TeV") == 0 ){
                    fBGFitRangeLeft[0]          = 0.34;
                    fBGFitRangeLeft[1]          = 0.44;
                    fBGFitRange[0]              = 0.67;
                    fBGFitRange[1]              = 0.795;
                } else if( fEnergyFlag.Contains("pPb_5.023TeV")){
                    fBGFitRangeLeft[0]          = 0.35;
                    fBGFitRangeLeft[1]          = 0.42;
                    fBGFitRange[0]              = 0.67;
                    fBGFitRange[1]              = 0.79;
                } else {
                    fBGFitRangeLeft[0]          = 0.34;
                    fBGFitRangeLeft[1]          = 0.44;
                    fBGFitRange[0]              = 0.67;
                    fBGFitRange[1]              = 0.79;
                }
            } else if ( mode == 5) {
                fBGFitRange[0]                  = 0.62;
                fBGFitRange[1]                  = 0.79;
                fBGFitRangeLeft[0]              = 0.35;
                fBGFitRangeLeft[1]              = 0.48;
            } else if ( mode == 0 ) {
                if (fEnergyFlag.BeginsWith("8TeV")) {
                    fBGFitRange[0]              = 0.60;
                    fBGFitRange[1]              = 0.79;
                    fBGFitRangeLeft[0]          = 0.35;
                    fBGFitRangeLeft[1]          = 0.48;
                } else {
                    fBGFitRange[0]              = 0.58;
                    fBGFitRange[1]              = 0.79;
                    fBGFitRangeLeft[0]          = 0.35;
                    fBGFitRangeLeft[1]          = 0.48;
                }
            } else {
                fBGFitRange[0]          = 0.58;
                fBGFitRange[1]          = 0.79;
                fBGFitRangeLeft[0]      = 0.35;
                fBGFitRangeLeft[1]      = 0.48;
            }

            // Initialize default Plot range for meson
            fMesonPlotRange[0]          = 0.53;
            fMesonPlotRange[1]          = 0.560;

            // Initialize default Plot default integration ranges
            if (mode == 0){
                if(fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0){
                    fMesonIntDeltaRange[0]          = -0.047;
                    fMesonIntDeltaRange[1]          = 0.022;
                } else {
                    fMesonIntDeltaRange[0]          = -0.036;
                    fMesonIntDeltaRange[1]          = 0.018;
                }
                fMesonIntDeltaRangeWide[0]      = -0.068;
                fMesonIntDeltaRangeWide[1]      = 0.032;
                fMesonIntDeltaRangeNarrow[0]    = -0.033;
                fMesonIntDeltaRangeNarrow[1]    = 0.012;
                if(fEnergyFlag.BeginsWith("8TeV")){
                    fMesonIntDeltaRangeNarrow[0]    = -0.030;
                    fMesonIntDeltaRangeNarrow[1]    = 0.010;
                    fMesonIntDeltaRange[0]          = -0.045;
                    fMesonIntDeltaRange[1]          = 0.025;
                }
            } else if (mode == 3) {
                fMesonIntDeltaRange[0]          = -0.080;
                fMesonIntDeltaRange[1]          =  0.040;
                fMesonIntDeltaRangeWide[0]      = -0.100;
                fMesonIntDeltaRangeWide[1]      =  0.060;
                fMesonIntDeltaRangeNarrow[0]    = -0.033;
                fMesonIntDeltaRangeNarrow[1]    =  0.012;
            } else if (mode == 2 || mode == 13){
                fMesonIntDeltaRange[0]          = -0.060;
                fMesonIntDeltaRange[1]          = 0.055;
                fMesonIntDeltaRangeWide[0]      = -0.080;
                fMesonIntDeltaRangeWide[1]      = 0.065;
                fMesonIntDeltaRangeNarrow[0]    = -0.040;
                fMesonIntDeltaRangeNarrow[1]    = 0.045;
                if( fEnergyFlag.Contains("5TeV2017") ){
                    fMesonIntDeltaRange[0]          *= 1.2;
                    fMesonIntDeltaRange[1]          *= 1.2;
                    fMesonIntDeltaRangeWide[0]      *= 1.2;
                    fMesonIntDeltaRangeWide[1]      *= 1.2;
                    fMesonIntDeltaRangeNarrow[0]    *= 1.2;
                    fMesonIntDeltaRangeNarrow[1]    *= 1.2;
                }
                if( fEnergyFlag.CompareTo("7TeV") == 0 ){
                    if( trigger.CompareTo("52") == 0 || triggerSet == 1  || triggerSet == 2 ){
                        fMesonIntDeltaRange[0]          *= 1.2;
                        fMesonIntDeltaRange[1]          *= 1.2;
                        fMesonIntDeltaRangeWide[0]      *= 1.2;
                        fMesonIntDeltaRangeWide[1]      *= 1.2;
                        fMesonIntDeltaRangeNarrow[0]    *= 1.2;
                        fMesonIntDeltaRangeNarrow[1]    *= 1.2;
                    }
                }
                if( fEnergyFlag.BeginsWith("8TeV") ){
                    if( trigger.CompareTo("81") == 0 ||  trigger.CompareTo("52") == 0 || triggerSet == 1  || triggerSet == 2 ){
                        fMesonIntDeltaRange[0]          *= 1.2;
                        fMesonIntDeltaRange[1]          *= 1.2;
                        fMesonIntDeltaRangeWide[0]      *= 1.2;
                        fMesonIntDeltaRangeWide[1]      *= 1.2;
                        fMesonIntDeltaRangeNarrow[0]    *= 1.2;
                        fMesonIntDeltaRangeNarrow[1]    *= 1.2;
                    }
                }
                if( fEnergyFlag.BeginsWith("pPb_8TeV") ){
                    if( trigger.CompareTo("8d") == 0 ||  trigger.CompareTo("8e") == 0 || triggerSet == 1  || triggerSet == 2 ){
                        fMesonIntDeltaRange[0]          *= 1.2;
                        fMesonIntDeltaRange[1]          *= 1.2;
                        fMesonIntDeltaRangeWide[0]      *= 1.2;
                        fMesonIntDeltaRangeWide[1]      *= 1.2;
                        fMesonIntDeltaRangeNarrow[0]    *= 1.2;
                        fMesonIntDeltaRangeNarrow[1]    *= 1.2;
                    } else {
                        fMesonIntDeltaRange[0]          *= 1.1;
                        fMesonIntDeltaRange[1]          *= 1.1;
                        fMesonIntDeltaRangeWide[0]      *= 1.1;
                        fMesonIntDeltaRangeWide[1]      *= 1.1;
                        fMesonIntDeltaRangeNarrow[0]    *= 1.1;
                        fMesonIntDeltaRangeNarrow[1]    *= 1.1;
                    }
                }
            } else if (mode == 4 || mode == 12 ) {
                fMesonIntDeltaRange[0]          = -0.080;
                fMesonIntDeltaRange[1]          = 0.080;
                fMesonIntDeltaRangeWide[0]      = -0.10;
                fMesonIntDeltaRangeWide[1]      = 0.10;
                fMesonIntDeltaRangeNarrow[0]    = -0.060;
                fMesonIntDeltaRangeNarrow[1]    = 0.060;
                if( fEnergyFlag.CompareTo("7TeV") == 0 || fEnergyFlag.BeginsWith("8TeV") ){
                    if( trigger.CompareTo("81") == 0 ||  trigger.CompareTo("52") == 0 || triggerSet == 1  || triggerSet == 2 ){
                        fMesonIntDeltaRange[0]          *= 0.8;
                        fMesonIntDeltaRange[1]          *= 0.8;
                        fMesonIntDeltaRangeWide[0]      *= 0.8;
                        fMesonIntDeltaRangeWide[1]      *= 0.8;
                        fMesonIntDeltaRangeNarrow[0]    *= 0.8;
                        fMesonIntDeltaRangeNarrow[1]    *= 0.8;
                    }
                }
            } else if (mode == 5){
                fMesonIntDeltaRange[0]          = -0.05;
                fMesonIntDeltaRange[1]          = 0.035;
                fMesonIntDeltaRangeWide[0]      = -0.068;
                fMesonIntDeltaRangeWide[1]      = 0.040;
                fMesonIntDeltaRangeNarrow[0]    = -0.04;
                fMesonIntDeltaRangeNarrow[1]    = 0.02;
            } else {
                fMesonIntDeltaRange[0]          = -0.048;
                fMesonIntDeltaRange[1]          = 0.022;
                fMesonIntDeltaRangeWide[0]      = -0.068;
                fMesonIntDeltaRangeWide[1]      = 0.032;
                fMesonIntDeltaRangeNarrow[0]    = -0.033;
                fMesonIntDeltaRangeNarrow[1]    = 0.012;
            }

            // Set meson mass ranges for fitting and plotting
            fMesonMassPlotRange[0]      = 0.35;
            fMesonMassPlotRange[1]      = 0.79;
            fMesonMassRange[0]          = 0.35;
            fMesonMassRange[1]          = 0.79;

            // Set Meson fit range
            if (mode == 2 || mode == 13){
                fMesonFitRange[0]               = 0.35;
                fMesonFitRange[1]               = 0.7;
                if( fEnergyFlag.CompareTo("7TeV") == 0 ){
                    fMesonFitRange[1] = 0.72;
                    if(  trigger.CompareTo("52") == 0  || triggerSet == 1)
                    fMesonFitRange[0] = 0.4;
                }
                if( fEnergyFlag.BeginsWith("8TeV") ){
                    fMesonFitRange[1] = 0.72;
                    if(  trigger.CompareTo("52") == 0  || triggerSet == 1)
                    fMesonFitRange[0] = 0.4;
                    if(  trigger.CompareTo("81") == 0  || triggerSet == 2){
                        fMesonFitRange[0] = 0.35;
                        fMesonFitRange[1] = 0.73;
                    }
                }
            } else if (mode == 4 || mode == 12 ) {
                if( fEnergyFlag.CompareTo("7TeV") == 0 ){
                    fMesonFitRange[0]           = 0.37;
                    fMesonFitRange[1]           = 0.72;
                    if( trigger.CompareTo("52") == 0  || triggerSet == 1){
                    fMesonFitRange[0]         = 0.4;
                    fMesonFitRange[1]         = 0.7;
                    }
                } else if( fEnergyFlag.BeginsWith("8TeV") ){
                    fMesonFitRange[0]           = 0.34;
                    fMesonFitRange[1]           = 0.7;
                    if( trigger.CompareTo("52") == 0 || triggerSet == 1 ){
                    fMesonFitRange[0]         = 0.4;
                    fMesonFitRange[1]         = 0.74;
                    }else if( trigger.CompareTo("81") == 0  || triggerSet == 2){
                    fMesonFitRange[0]         = 0.4;
                    fMesonFitRange[1]         = 0.7;
                    }
                } else {
                    fMesonFitRange[0]           = 0.38;
                    fMesonFitRange[1]           = 0.73;
                }
                if( fEnergyFlag.Contains("13TeV")  ){
                    fMesonFitRange[0]                = 0.44;
                    fMesonFitRange[1]                = 0.7;
                }else if( fEnergyFlag.Contains("PbPb") ){
                    fMesonFitRange[0]                = 0.44;
                    fMesonFitRange[1]                = 0.66;
                }else{
                    fMesonFitRange[0]                = 0.4;
                    fMesonFitRange[1]                = 0.7;
                }
            } else {
                fMesonFitRange[0]           = 0.4;
                fMesonFitRange[1]           = 0.7;
            }

            // Set remaining parameters for fitting
            if (mode == 0){
                fMesonWidthExpect           = 0.005;
                if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0)
                    fMesonWidthExpect = 0.010;
                fMesonLambdaTail            = 0.007;
                fMesonWidthRange[0]         = 0.002;
                fMesonWidthRange[1]         = 0.020;
                fMesonLambdaTailRange[0]    = 0.004;
                fMesonLambdaTailRange[1]    = 0.03;
            } else if (mode == 2 || mode == 13){
                fMesonWidthRange[0]             = 0.010;
                fMesonWidthRange[1]             = 0.050;
                if (fEnergyFlag.Contains("5TeV2017") || fEnergyFlag.CompareTo("7TeV") == 0 || fEnergyFlag.BeginsWith("8TeV") || fEnergyFlag.BeginsWith("pPb_8TeV")){
                    fMesonWidthExpect           = 0.020;
                    fMesonLambdaTail            = 0.025;
                    fMesonLambdaTailRange[0]    = 0.025;
                    fMesonLambdaTailRange[1]    = 0.025;
                    if( trigger.CompareTo("81") == 0  || triggerSet == 2){
                        fMesonLambdaTail            = 0.02;
                        fMesonLambdaTailRange[0]    = 0.02;
                        fMesonLambdaTailRange[1]    = 0.02;
                        fMesonWidthExpect           = 0.025;
                    } else if ( trigger.CompareTo("52") == 0  || triggerSet == 1){
                        fMesonLambdaTail            = 0.02;
                        fMesonLambdaTailRange[0]    = 0.02;
                        fMesonLambdaTailRange[1]    = 0.02;
                        fMesonWidthExpect           = 0.025;
                    }
                } else if (fEnergyFlag.Contains("13TeV") ){
                    fMesonWidthExpect           = 0.03;
                    fMesonLambdaTail            = 0.04;
                    fMesonLambdaTailRange[0]    = 0.03;
                    fMesonLambdaTailRange[1]    = 0.05;
                } else if (fEnergyFlag.Contains("pPb_5.023TeV") ){
                    fMesonWidthExpect               = 0.020;
                    fMesonWidthRange[0]             = 0.010;
                    fMesonWidthRange[1]             = 0.060;
                    if (fIsMC > 0){
                        fMesonLambdaTail                = 0.032;
                        fMesonLambdaTailRange[0]        = 0.032;
                        fMesonLambdaTailRange[1]        = 0.032;
                        fMesonLambdaTailMC              = 0.032;
                        fMesonLambdaTailRangeMC[0]      = 0.032;
                        fMesonLambdaTailRangeMC[1]      = 0.032;
                        fMesonWidthRangeMC[0]           = 0.005;
                    } else {
                        fMesonLambdaTail                = 0.030;
                        fMesonLambdaTailRange[0]        = 0.030;
                        fMesonLambdaTailRange[1]        = 0.030;
                    }
                } else {
                    fMesonWidthExpect               = 0.020;
                    fMesonLambdaTail                = 0.020;
                    fMesonLambdaTailRange[0]        = 0.018;
                    fMesonLambdaTailRange[1]        = 0.020;
                }
            } else if (mode == 3) {
                fMesonWidthExpect               = 0.018;
                fMesonWidthRange[0]             = 0.012;
                fMesonWidthRange[1]             = 0.025;
                fMesonLambdaTail                = 0.007;
                fMesonLambdaTailRange[0]        = 0.007;
                fMesonLambdaTailRange[1]        = 0.007;
                if (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) {
                fMesonWidthExpect           = 0.017;
                fMesonWidthRange[0]         = 0.004;
                fMesonWidthRange[1]         = 0.032;
                fMesonLambdaTail                = 0.007;
                fMesonLambdaTailRange[0]        = 0.006;
                fMesonLambdaTailRange[1]        = 0.009;
                }
            } else if (mode == 4 || mode == 12 ) {
                if (fEnergyFlag.CompareTo("2.76TeV") == 0 || fEnergyFlag.CompareTo("7TeV") == 0 || fEnergyFlag.BeginsWith("8TeV") ){
                    fMesonWidthExpect           = 0.025;
                    fMesonWidthRange[0]         = 0.018;
                    fMesonWidthRange[1]         = 0.070;
                    fMesonLambdaTail            = 0.025;
                    fMesonLambdaTailRange[0]    = 0.025;
                    fMesonLambdaTailRange[1]    = 0.025;
                } else if (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0){
                    fMesonWidthExpect           = 0.036;
                    fMesonWidthRange[0]         = 0.010;
                    fMesonWidthRange[1]         = 0.140;
                    fMesonLambdaTail            = 0.033;
                    fMesonLambdaTailRange[0]    = 0.033;
                    fMesonLambdaTailRange[1]    = 0.033;
                    fMesonLambdaTailMC          = 0.033;
                    fMesonLambdaTailRangeMC[0]  = 0.010;
                    fMesonLambdaTailRangeMC[1]  = 0.080;
                } else if (mode == 12 && fEnergyFlag.CompareTo("5TeV2017") == 0){
                    fMesonWidthExpect           = 0.025;
                    fMesonWidthRange[0]         = 0.018;
                    fMesonWidthRange[1]         = 0.050;
                    fMesonLambdaTail            = 0.02;
                    fMesonLambdaTailRange[0]    = 0.001;
                    fMesonLambdaTailRange[1]    = 0.1;
                    fMesonLambdaTailMC          = 0.02;
                    fMesonLambdaTailRangeMC[0]  = 0.001;
                    fMesonLambdaTailRangeMC[1]  = 0.1;
                } else {
                    fMesonWidthExpect           = 0.025;
                    fMesonLambdaTail            = 0.012;
                    fMesonWidthRange[0]         = 0.018;
                    fMesonWidthRange[1]         = 0.070;
                    fMesonLambdaTailRange[0]    = 0.001;
                    fMesonLambdaTailRange[1]    = 0.025;
                }
            } else {
                fMesonWidthExpect           = 0.005;
                fMesonLambdaTail            = 0.007;
                fMesonWidthRange[0]         = 0.002;
                fMesonWidthRange[1]         = 0.020;
                fMesonLambdaTailRange[0]    = 0.004;
                fMesonLambdaTailRange[1]    = 0.03;
            }


        //****************************************************************************************************
        // Initialization for eta' meson
        //****************************************************************************************************
        } else if (setPi0.CompareTo("EtaPrime") == 0){

            // set meson ID according to PDG
            fMesonId                    = 331;

            switch(mode) {
                case 0 : // PCM-PCM
                    fMesonWidthExpect = 0.05;
                    fMesonLambdaTail  = 0.007;
                    SetRangeArray(fPeakRange               ,  0.92  , 1.00);
                    SetRangeArray(fIntFixedRange           ,  0.75  , 1.10); // not optimized
                    SetRangeArray(fFitRange                ,  0.75  , 1.10);
                    SetRangeArray(fBGFitRange              ,  1.07  , 1.17);
                    SetRangeArray(fBGFitRangeLeft          ,  0.75  , 0.80);
                    SetRangeArray(fMesonPlotRange          ,  0.71  , 1.10);
                    SetRangeArray(fMesonIntDeltaRange      , -0.10  , 0.10);
                    SetRangeArray(fMesonIntDeltaRangeWide  , -0.15  , 0.15);
                    SetRangeArray(fMesonIntDeltaRangeNarrow, -0.08  , 0.08);
                    SetRangeArray(fMesonMassRange          ,  0.71  , 1.10);
                    SetRangeArray(fMesonMassPlotRange      ,  0.71  , 1.17);
                    SetRangeArray(fMesonFitRange           ,  0.75  , 1.10);
                    SetRangeArray(fMesonWidthRange         ,  0.002 , 0.10);
                    SetRangeArray(fMesonLambdaTailRange    ,  0.0005, 0.03);
                    SetRangeArray(fMidPt                   ,  1.20  , 2.50);
                    break;
                case 2 : // PCM-EMC
                    fMesonWidthExpect = 0.05;
                    fMesonLambdaTail  = 0.007;
                    SetRangeArray(fPeakRange               ,  0.82  , 1.02);
                    SetRangeArray(fIntFixedRange           ,  0.75  , 1.10); // not optimized
                    SetRangeArray(fFitRange                ,  0.75  , 1.10);
                    SetRangeArray(fBGFitRange              ,  1.07  , 1.17);
                    SetRangeArray(fBGFitRangeLeft          ,  0.75  , 0.80);
                    SetRangeArray(fMesonPlotRange          ,  0.71  , 1.10);
                    SetRangeArray(fMesonIntDeltaRange      , -0.10  , 0.10);
                    SetRangeArray(fMesonIntDeltaRangeWide  , -0.15  , 0.15);
                    SetRangeArray(fMesonIntDeltaRangeNarrow, -0.08  , 0.08);
                    SetRangeArray(fMesonMassRange          ,  0.71  , 1.10);
                    SetRangeArray(fMesonMassPlotRange      ,  0.71  , 1.17);
                    SetRangeArray(fMesonFitRange           ,  0.75  , 1.10);
                    SetRangeArray(fMesonWidthRange         ,  0.002 , 0.10);
                    SetRangeArray(fMesonLambdaTailRange    ,  0.0005, 0.03);
                    SetRangeArray(fMidPt                   ,  1.20  , 2.50);
                    break;
                case 3 : // PCM-PHOS
                    fMesonWidthExpect = 0.05;
                    fMesonLambdaTail  = 0.007;
                    SetRangeArray(fPeakRange               ,  0.91  , 1.03);
                    SetRangeArray(fIntFixedRange           ,  0.75  , 1.10); // not optimized
                    SetRangeArray(fFitRange                ,  0.75  , 1.10);
                    SetRangeArray(fBGFitRange              ,  1.07  , 1.17);
                    SetRangeArray(fBGFitRangeLeft          ,  0.75  , 0.80);
                    SetRangeArray(fMesonPlotRange          ,  0.71  , 1.10);
                    SetRangeArray(fMesonIntDeltaRange      , -0.10  , 0.10);
                    SetRangeArray(fMesonIntDeltaRangeWide  , -0.15  , 0.15);
                    SetRangeArray(fMesonIntDeltaRangeNarrow, -0.08  , 0.08);
                    SetRangeArray(fMesonMassRange          ,  0.71  , 1.10);
                    SetRangeArray(fMesonMassPlotRange      ,  0.71  , 1.17);
                    SetRangeArray(fMesonFitRange           ,  0.75  , 1.10);
                    SetRangeArray(fMesonWidthRange         ,  0.002 , 0.10);
                    SetRangeArray(fMesonLambdaTailRange    ,  0.0005, 0.03);
                    SetRangeArray(fMidPt                   ,  1.20  , 2.50);
                    break;
                case 4 : // EMC-EMC
                    fMesonWidthExpect = 0.05;
                    fMesonLambdaTail  = 0.007;
                    SetRangeArray(fPeakRange               ,  0.78  , 1.00);
                    SetRangeArray(fIntFixedRange           ,  0.75  , 1.10); // not optimized
                    SetRangeArray(fFitRange                ,  0.75  , 1.10);
                    SetRangeArray(fBGFitRange              ,  1.07  , 1.17);
                    SetRangeArray(fBGFitRangeLeft          ,  0.75  , 0.80);
                    SetRangeArray(fMesonPlotRange          ,  0.71  , 1.10);
                    SetRangeArray(fMesonIntDeltaRange      , -0.10  , 0.10);
                    SetRangeArray(fMesonIntDeltaRangeWide  , -0.15  , 0.15);
                    SetRangeArray(fMesonIntDeltaRangeNarrow, -0.08  , 0.08);
                    SetRangeArray(fMesonMassRange          ,  0.71  , 1.10);
                    SetRangeArray(fMesonMassPlotRange      ,  0.71  , 1.17);
                    SetRangeArray(fMesonFitRange           ,  0.75  , 1.10);
                    SetRangeArray(fMesonWidthRange         ,  0.002 , 0.10);
                    SetRangeArray(fMesonLambdaTailRange    ,  0.0005, 0.03);
                    SetRangeArray(fMidPt                   ,  1.20  , 2.50);
                    break;
                case 5 : // PHOS-PHOS
                    fMesonWidthExpect = 0.05;
                    fMesonLambdaTail  = 0.007;
                    SetRangeArray(fPeakRange               ,  0.90  , 1.04);
                    SetRangeArray(fIntFixedRange           ,  0.75  , 1.10); // not optimized
                    SetRangeArray(fFitRange                ,  0.75  , 1.10);
                    SetRangeArray(fBGFitRange              ,  1.07  , 1.17);
                    SetRangeArray(fBGFitRangeLeft          ,  0.75  , 0.80);
                    SetRangeArray(fMesonPlotRange          ,  0.71  , 1.10);
                    SetRangeArray(fMesonIntDeltaRange      , -0.10  , 0.10);
                    SetRangeArray(fMesonIntDeltaRangeWide  , -0.15  , 0.15);
                    SetRangeArray(fMesonIntDeltaRangeNarrow, -0.08  , 0.08);
                    SetRangeArray(fMesonMassRange          ,  0.71  , 1.10);
                    SetRangeArray(fMesonMassPlotRange      ,  0.71  , 1.17);
                    SetRangeArray(fMesonFitRange           ,  0.75  , 1.10);
                    SetRangeArray(fMesonWidthRange         ,  0.002 , 0.10);
                    SetRangeArray(fMesonLambdaTailRange    ,  0.0005, 0.03);
                    SetRangeArray(fMidPt                   ,  1.20  , 2.50);
                    break;
                default :
                    std::cout << "FATAL ERROR: fit ranges not defined for  \"" << setPi0 << "\", mode " << mode << std::endl;
                    std::cout << "--> ExtractSignalV2.h: line " << __LINE__ << std::endl;
                    std::terminate();
            }
        }
        if (fMesonLambdaTailMC == -10)
            fMesonLambdaTailMC          = fMesonLambdaTail;
        fMesonLambdaTailRangeNominal[0] = fMesonLambdaTailRange[0];
        fMesonLambdaTailRangeNominal[1] = fMesonLambdaTailRange[1];
        if (fMesonLambdaTailRangeMC[0] == -10)
            fMesonLambdaTailRangeMC[0]  = fMesonLambdaTailRange[0];
        if (fMesonLambdaTailRangeMC[1] == -10)
            fMesonLambdaTailRangeMC[1]  = fMesonLambdaTailRange[1];

        fMesonWidthExpectMC     = fMesonWidthExpect;
        if (fMesonWidthRangeMC[0] == -10)
            fMesonWidthRangeMC[0]   = fMesonWidthRange[0];
        if (fMesonWidthRangeMC[1] == -10)
            fMesonWidthRangeMC[1]   = fMesonWidthRange[1];
    }
#endif
