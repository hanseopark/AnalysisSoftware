//****************************************************************************
//********* provided by Gamma Conversion Group *******************************
//********* $ALICE_ROOT/PWGGA/GammaConv **************************************
//***** https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion *******
//****************************************************************************

//****************************************************************************
//******************** Global variable for setup of macro ********************
//****************************************************************************
fstream     fFileErrLog;
fstream     fFileDataLog;
TString     fTextCent                                                   = "";
TString     fEnergyFlag                                                 = "";
TString     fPrefix                                                     = "";
TString     fPrefix2                                                    = "";
TString     fPeriodFlag                                                 = "";
TString     fdate                                                       = "";
TString     fCollisionSystem                                            = "";
TString     fTextMeasurement                                            = "";
TString     fDetectionProcess                                           = "";
TString     fDetectionProcessPtBins                                     = "";
TString     fDecayChannel                                               = "";
TString     fCutSelection                                               = "";
Int_t       fMode                                                       = -1;
Bool_t      fAdvancedMesonQA                                            = kFALSE;
Bool_t      fAdvancedClusterQA                                          = kFALSE;
Bool_t      fAdditionalLabels                                           = kFALSE;
Bool_t      fNewMCOutput                                                = kFALSE;
TString     fEventCutSelection                                          = "";
TString     fClusterCutSelection                                        = "";
TString     fClusterMergedCutSelection                                  = "";
TString     fMesonCutSelection                                          = "";
TString     fAdditionalName                                             = "";

Int_t       fIsMC                                                       = 0;
Int_t       fIsMCGammaTrig                                              = 0;
Int_t       fMesonId                                                    = 0;
Double_t    fYMaxMeson                                                  = 0.9;
Double_t    fMesonMassExpect                                            = 0;   // expected meson mass
Double_t    fNEvents                                                    = 0;
TString     fcMonth[12]                                                 = {"Jan","Feb","Mar","Apr","May","Jun",
                                                                           "Jul","Aug","Sep","Oct","Nov","Dec"};
TString     fCentralityString                                           = "";
TString     fNLMString                                                  = "";
Int_t       fNLMmin                                                     = 0;
TFile*      fOutput1                                                    = NULL;
TFile*      fOutput2                                                    = NULL;
Double_t    fYields                                                     = 0.;
Double_t    fYieldsError                                                = 0.;
TString     labelsBG[10]                                                = { "PiCh", "P", "KCh", "N", "K0s", 
                                                                            "Lambda", "Mu", "K0l", "Rest", "AllNonEM"};
TString     labelsGamma[8]                                              = { "DirGamma", "PiN", "Eta", "EtaPrim", "Omega", 
                                                                            "Phi", "Lambda", "Rest"};
TString     labelsElectron[9]                                           = { "DirE", "Gamma", "PiN", "Eta", "C", 
                                                                            "B", "WOrZ", "Tau", "Rest"};
TString     nameSecondaries[4]                                          = {"K0S", "Lambda", "K0L", "Rest"};

//****************************************************************************
//************************* MC object names **********************************
//****************************************************************************
TString     fObjectNameMCMesonAcc                                       = "";
TString     fObjectNameMCMeson                                          = "";
TString     fObjectNameMCMesonAccSecPi0                                 = "";
TString     fObjectNameMCMesonSecPi0                                    = "";
TString     fObjectNameMCMesonWOWeights                                 = "";
TString     fObjectNameMCMesonDalitzAcc                                 = "";
TString     fObjectNameMCMesonDalitz                                    = "";
TString     fObjectNameMCMesonDalitzWOWeights                           = "";
TString     fObjectNameTrueMergedM02                                    = "";
TString     fObjectNameTrueFromPi0M02                                   = "";
TString     fObjectNameTrueFromPi0DCM02                                 = "";
TString     fObjectNameTrueFromPi0GGM02                                 = "";
TString     fObjectNameTrueFromPi0DalitzM02                             = "";
TString     fObjectNameTrueFromEtaM02                                   = "";
TString     fObjectNameTrueFromEtaDCM02                                 = "";
TString     fObjectNameTrueFromEtaGGM02                                 = "";
TString     fObjectNameTrueFromEtaDalitzM02                             = "";
TString     fObjectNameTrueMergedPartConvM02                            = "";
TString     fObjectNameTrueMergedPartConvFromPi0M02                     = "";
TString     fObjectNameTrueMergedPartConvFromEtaM02                     = "";
TString     fObjectNameTrueMergedPureM02                                = "";
TString     fObjectNameTrueMergedPureFromPi0M02                         = "";
TString     fObjectNameTrueMergedPureFromEtaM02                         = "";
TString     fObjectNameTrueMergedPartConvLeadEM02                       = "";
TString     fObjectNameTrueMergedOneGammaM02                            = "";
TString     fObjectNameTrueMergedOneGammaFromPi0M02                     = "";
TString     fObjectNameTrueMergedOneGammaFromEtaM02                     = "";
TString     fObjectNameTrueMergedOneElectronM02                         = "";
TString     fObjectNameTrueMergedOneElectronFromPi0M02                  = "";
TString     fObjectNameTrueMergedOneElectronFromEtaM02                  = "";
TString     fObjectNameTrueClusBGM02                                    = "";
TString     fObjectNameTrueClusGammaM02                                 = "";
TString     fObjectNameTrueClusElectronM02                              = "";
TString     fObjectNameTrueClusBG_Source                                = "";
TString     fObjectNameTrueClusGamma_Source                             = "";
TString     fObjectNameTrueClusElectron_Source                          = "";
TString     fObjectNameTrueClusPrimMesonM02                             = "";
TString     fObjectNameTrueClusSecMesonM02                              = "";
TString     fObjectNameTrueClusSecMesonFromK0sM02                       = "";
TString     fObjectNameTrueClusSecMesonFromK0lM02                       = "";
TString     fObjectNameTrueClusSecMesonFromLambdaM02                    = "";


//****************************************************************************
//******************************** Functions *********************************
//****************************************************************************
void CreatePtHistos();                                                                                      // Creat pt dependent histograms
void FillPtHistos();                                                                                        // Fill pt dependent histograms
void Initialize(TString setPi0, Int_t, Int_t);                                                              // Initialization of global variables depending on meson analysed 
Double_t CrystalBallBck(Double_t *,Double_t *);                                                             // Definition of CrystalBall with linear BG
Double_t LinearBackground(Double_t *,Double_t *);                                                           // Definition of linear BG
Double_t CrystalBall(Double_t *,Double_t *);                                                                // Definition of CrystalBall
void DrawAutoGammaHistoPaper2D( TH2* histo1,
                                TString Title, TString XTitle, TString YTitle,
                                Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
                                Bool_t YRange, Double_t YMin ,Double_t YMax,
                                Bool_t XRange, Double_t XMin, Double_t XMax, Double_t xOffset=1, Double_t yOffset=1.);
Double_t FindSmallestEntryIn2D(TH2F* histo);
void FillDataHistosArray(TH2F*,TH2F*);                                                                      // Fill invariant mass histograms for Signal 
TH1D* FillProjectionX (TH2F*, TString, Double_t, Double_t, Int_t);                                          // Fill Projection in according to Y bins
TH1D* FillProjectionY (TH2F*, TString, Double_t, Double_t, Int_t);                                          // Fill Projection in according to X bins
TH1D* AddPossiblyNotExistingHists(TH1D*, TH1D*, TString );                                                  // Add possibly none existing histos
void CheckForNULLForPointer(TH1D* );                                                                        // Check if histo is NULL
void FillMCM02HistosArray( TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*);                                       // Fill M02 histograms for MC 
void FillMCM02AdditionHistosArray( TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*,    // Fill additional M02 histograms for MC     
                                   TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F* );                                                   
void FillMCPrimSecM02HistosArray( TH2F*, TH2F** );                                                          // Fill M02 histograms for MC for Pi0 Prim Sec 
void FillMCBGSeparated(TH2F* );                                                                             // Fill MC BG separated 
void FillMCElectronSeparated(TH2F* );                                                                       // Fill MC Electron separated 
void FillMCGammaSeparated(TH2F* );                                                                          // Fill MC Gamma separated 
void FillHistosArrayMC(TH1D*, TH1D * , TH1D* );                                                             // Rebin MC input histo's
void FillHistosArrayMCSecAndCalcAcceptance(TH2D*, TH2D * );                                                 // Rebin Sec MC input histo's and calc acceptance
void IntegrateHistoInvMass(TH1D*, Double_t* );                                                              // Integrate invariant mass histogram
void IntegrateHistoM02(TH1D*, Double_t* );                                                                  // Integrate M02 histogram
void CalculateMesonAcceptance();                                                                            // Calculation of meson acceptance
TH1D* CalculateMesonEfficiency(TH1D*, TH1D*,TString);                                                       // Calculation of meson efficiencies 
TH1D* CalculateMesonEfficiencySec(TH1D*, TH1D*,TString);                                                    // Calculation of meson efficiencies 
TH1D* CalculatePurity(TH1D*, TH1D*,TString);                                                                // Calculation of purity
TH1D* CalculateSecondaryFractions(TH1D*, TH1D*, TString );                                                  // Calculate fraction of secondaries
void SaveHistos(Int_t, TString, TString);                                                                   // Saving standard histograms to a file
void SaveCorrectionHistos(TString , TString);                                                               // Saving correction histograms to a file
void Delete();                                                                                              // Deleting all pointers
void SetCorrectMCHistogrammNames(TString);                                                                  // Setting correct histogram names
Bool_t LoadSecondaryPionsFromExternalFile();                                                                // Loads secondary neutral pion input graphs from file
void CreateRatioHistos();                                                                                   // creates ratio histos for MC inputs

//****************************************************************************
//************************** input histograms ********************************
//****************************************************************************
TH2F*       fHistoInvMassVsPt                                           = NULL;
TH1D*       fNumberOfGoodESDTracks                                      = NULL;
TH1D*       fEventQuality                                               = NULL;

//****************************************************************************
//*********************** variables for Plotting header **********************
//****************************************************************************
Double_t    fIntLinearBck;
Double_t    fIntLinearBckError;
Double_t*   fBGFitRangeLeft                                             = NULL;
Double_t*   fBGFitRange                                                 = NULL;
Double_t*   fMesonIntDeltaRange                                         = NULL;
Double_t*   fMesonIntDeltaRangeWide                                     = NULL;
Double_t*   fMesonIntDeltaRangeNarrow                                   = NULL;
TH1D*       fCopySignal                                                 = NULL;
TH1D*       fCopyOnlyBG                                                 = NULL;
Double_t*   fMesonMassPlotRange                                         = NULL;
Double_t*   fMesonM02PlotRange                                          = NULL;
Double_t*   fMesonFitRange                                              = NULL;
TF1 *       fFitReco                                                    = NULL;
TF1 *       fFitGausExp                                                 = NULL;
TF1 *       fFitLinearBck                                               = NULL;
TF1 *       fFitLinearBckExcl                                           = NULL;
Double_t    fMesonWidthExpect                                           = 0;
Double_t    fMesonLambdaTail                                            = 0;
Double_t*   fMesonWidthRange                                            = NULL;
Double_t*   fMesonLambdaTailRange                                       = NULL;
Double_t*   fMesonMassRange                                             = NULL;
Double_t*   fMesonMassIntRange                                          = NULL;
Double_t*   fMesonM02IntRange                                           = NULL;
Double_t*   fMesonMass                                                  = NULL;
Double_t*   fMesonMassLeft                                              = NULL;
Double_t*   fMesonTrueMass                                              = NULL;

//****************************************************************************
//***************************** histos in Pt Bins ****************************
//****************************************************************************
TH1D**      fHistoInvMassPtBin                                          = NULL;    
TH1D**      fHistoM02PtBin                                              = NULL;    

//****************************************************************************
//**************************** global fit functions **************************
//****************************************************************************
Double_t    fFWHMFunc;
Double_t    fFWHMFuncError;

//*****************************************************************************
//*********************** dedicated cluster histograms ************************
//*****************************************************************************
TH1D*       fHistoClustersPt                                            = NULL;
TH1D*       fHistoClustersOverlapHeadersPt                              = NULL;
TH2F*       fHistoClustersMergedPtM02                                   = NULL;
TH2F*       fHistoClustersMergedPtM02AccMeson                           = NULL;
TH2F*       fHistoClustersMergedEM02AccMeson                            = NULL;
TH1D*       fHistoClusterCandidates                                     = NULL;
TH1D*       fHistoClusterMergedCandidates                               = NULL;
TH2F*       fHistoTrackVsClusterCandidates                              = NULL;
TH2F*       fHistoSPDtrackletvsSPDclusters                              = NULL;

//*****************************************************************************
//************************* QA histos for clusters ****************************
//*****************************************************************************
TH2F*       fHistoClusterNCellsPt                                       = NULL;
TH2F*       fHistoClusterMergedNCellsPt                                 = NULL;
TH2F*       fHistoClusterMergedNCellsArClPt                             = NULL;
TH2F*       fHistoClusterMergedNCellsArAInclPt                          = NULL;
TH2F*       fHistoClusterMergedEAroundClE                               = NULL;
TH2F*       fHistoClusterMergedNPartPt                                     = NULL;

//*****************************************************************************
//************************* true MC histograms ********************************
//*****************************************************************************
TH2F*       fHistoTrueClustersMergedPtM02                               = NULL;
TH2F*       fHistoTrueClustersPi0PtM02                                  = NULL;
TH2F*       fHistoTrueClustersPi0DCPtM02                                = NULL;
TH2F*       fHistoTrueClustersPi0GGPtM02                                = NULL;
TH2F*       fHistoTrueClustersPi0DalitzPtM02                            = NULL;
TH2F*       fHistoTrueClustersEtaPtM02                                  = NULL;
TH2F*       fHistoTrueClustersEtaDCPtM02                                = NULL;
TH2F*       fHistoTrueClustersEtaGGPtM02                                = NULL;
TH2F*       fHistoTrueClustersEtaDalitzPtM02                            = NULL;
TH2F*       fHistoTrueClustersBGPtM02                                   = NULL;
TH2F*       fHistoTrueClustersBGPtSource                                = NULL;
TH1D*       fHistoTrueClustersBGPt[10]                                   = { NULL, NULL, NULL, NULL, NULL,
                                                                            NULL, NULL, NULL, NULL, NULL};
TH1D*       fHistoRatioTrueClustersBGPt[10]                             = { NULL, NULL, NULL, NULL, NULL,
                                                                            NULL, NULL, NULL, NULL, NULL};

TH2F*       fHistoTrueClustersGammaPtM02                                = NULL;
TH2F*       fHistoTrueClustersGammaPtSource                             = NULL;
TH1D*       fHistoTrueClustersGammaPt[8]                                = { NULL, NULL, NULL, NULL, NULL,
                                                                            NULL, NULL, NULL};
TH2F*       fHistoTrueClustersElectronPtM02                             = NULL;
TH2F*       fHistoTrueClustersElectronPtSource                          = NULL;
TH1D*       fHistoTrueClustersElectronPt[9]                             = { NULL, NULL, NULL, NULL, NULL,
                                                                            NULL, NULL, NULL, NULL};
TH2F*       fHistoTrueClusPartConvMergedPtM02                           = NULL;
TH2F*       fHistoTrueClusPartConvMergedFromPi0PtM02                    = NULL;
TH2F*       fHistoTrueClusPartConvMergedFromEtaPtM02                    = NULL;
TH2F*       fHistoTrueClusPureMergedPtM02                               = NULL;
TH2F*       fHistoTrueClusPureMergedFromPi0PtM02                        = NULL;
TH2F*       fHistoTrueClusPureMergedFromEtaPtM02                        = NULL;
TH2F*       fHistoTrueClusOneGammaPtM02                                 = NULL;
TH2F*       fHistoTrueClusOneGammaFromPi0PtM02                          = NULL;
TH2F*       fHistoTrueClusOneGammaFromEtaPtM02                          = NULL;
TH2F*       fHistoTrueClusOneElectronPtM02                              = NULL;
TH2F*       fHistoTrueClusOneElectronFromPi0PtM02                       = NULL;
TH2F*       fHistoTrueClusOneElectronFromEtaPtM02                       = NULL;

TH2F*       fHistoTrueClustersPrimPi0PtM02                              = NULL;

//****************************************************************************
//***************************** histos in Pt Bins ****************************
//****************************************************************************
TH1D**      fHistoTrueClusMergedM02PtBin                                = NULL;
TH1D**      fHistoTrueClusPartConvMergedM02PtBin                        = NULL;
TH1D**      fHistoTrueClusPartConvMergedFromPi0M02PtBin                 = NULL;
TH1D**      fHistoTrueClusPartConvMergedFromEtaM02PtBin                 = NULL;
TH1D**      fHistoTrueClusPureMergedM02PtBin                            = NULL;
TH1D**      fHistoTrueClusPureMergedFromPi0M02PtBin                     = NULL;
TH1D**      fHistoTrueClusPureMergedFromEtaM02PtBin                     = NULL;
TH1D**      fHistoTrueClusOneGammaM02PtBin                              = NULL;
TH1D**      fHistoTrueClusOneGammaFromPi0M02PtBin                       = NULL;
TH1D**      fHistoTrueClusOneGammaFromEtaM02PtBin                       = NULL;
TH1D**      fHistoTrueClusOneElectronM02PtBin                           = NULL;
TH1D**      fHistoTrueClusOneElectronFromPi0M02PtBin                    = NULL;
TH1D**      fHistoTrueClusOneElectronFromEtaM02PtBin                    = NULL;
TH1D**      fHistoTrueClusPi0M02PtBin                                   = NULL;
TH1D**      fHistoTrueClusPi0DCM02PtBin                                 = NULL;
TH1D**      fHistoTrueClusPi0GGM02PtBin                                 = NULL;
TH1D**      fHistoTrueClusPi0DalitzM02PtBin                             = NULL;
TH1D**      fHistoTrueClusEtaM02PtBin                                   = NULL;
TH1D**      fHistoTrueClusEtaDCM02PtBin                                 = NULL;
TH1D**      fHistoTrueClusEtaGGM02PtBin                                 = NULL;
TH1D**      fHistoTrueClusEtaDalitzM02PtBin                             = NULL;
TH1D**      fHistoTrueClusGammaM02PtBin                                 = NULL;
TH1D**      fHistoTrueClusElectronM02PtBin                              = NULL;
TH1D**      fHistoTrueClusBGM02PtBin                                    = NULL;

TH1D**      fHistoTrueClusPrimPi0M02PtBin                               = NULL;

//****************************************************************************
//************************** Yields vs Pt ************************************
//****************************************************************************
Double_t*   fMesonM02Yields                                             = NULL;
Double_t*   fMesonM02TrueMergedYields                                   = NULL;
Double_t*   fMesonM02TruePi0Yields                                      = NULL;
Double_t*   fMesonM02TruePi0DCYields                                    = NULL;
Double_t*   fMesonM02TruePi0GGYields                                    = NULL;
Double_t*   fMesonM02TruePi0DalitzYields                                = NULL;
Double_t*   fMesonM02TrueEtaYields                                      = NULL;
Double_t*   fMesonM02TrueEtaDCYields                                    = NULL;
Double_t*   fMesonM02TrueEtaGGYields                                    = NULL;
Double_t*   fMesonM02TrueEtaDalitzYields                                = NULL;
Double_t*   fMesonM02TrueGammaYields                                    = NULL;
Double_t*   fMesonM02TrueElectronYields                                 = NULL;
Double_t*   fMesonM02TrueBGYields                                       = NULL;
Double_t*   fMesonM02YieldsError                                        = NULL;
Double_t*   fMesonM02TrueMergedYieldsError                              = NULL;
Double_t*   fMesonM02TruePi0YieldsError                                 = NULL;
Double_t*   fMesonM02TruePi0DCYieldsError                               = NULL;
Double_t*   fMesonM02TruePi0GGYieldsError                               = NULL;
Double_t*   fMesonM02TruePi0DalitzYieldsError                           = NULL;
Double_t*   fMesonM02TrueEtaYieldsError                                 = NULL;
Double_t*   fMesonM02TrueEtaDCYieldsError                               = NULL;
Double_t*   fMesonM02TrueEtaGGYieldsError                               = NULL;
Double_t*   fMesonM02TrueEtaDalitzYieldsError                           = NULL;
Double_t*   fMesonM02TrueGammaYieldsError                               = NULL;
Double_t*   fMesonM02TrueElectronYieldsError                            = NULL;
Double_t*   fMesonM02TrueBGYieldsError                                  = NULL;
Double_t*   fMesonM02TrueMergedPartConvYields                           = NULL;
Double_t*   fMesonM02TrueMergedPartConvYieldsError                      = NULL;
Double_t*   fMesonM02TrueMergedPartConvFromPi0Yields                    = NULL;
Double_t*   fMesonM02TrueMergedPartConvFromPi0YieldsError               = NULL;
Double_t*   fMesonM02TrueMergedPartConvFromEtaYields                    = NULL;
Double_t*   fMesonM02TrueMergedPartConvFromEtaYieldsError               = NULL;
Double_t*   fMesonM02TrueMergedPureYields                               = NULL;
Double_t*   fMesonM02TrueMergedPureYieldsError                          = NULL;
Double_t*   fMesonM02TrueMergedPureFromPi0Yields                        = NULL;
Double_t*   fMesonM02TrueMergedPureFromPi0YieldsError                   = NULL;
Double_t*   fMesonM02TrueMergedPureFromEtaYields                        = NULL;
Double_t*   fMesonM02TrueMergedPureFromEtaYieldsError                   = NULL;
Double_t*   fMesonM02TrueMergedOneGammaYields                           = NULL;
Double_t*   fMesonM02TrueMergedOneGammaYieldsError                      = NULL;
Double_t*   fMesonM02TrueMergedOneGammaFromPi0Yields                    = NULL;
Double_t*   fMesonM02TrueMergedOneGammaFromPi0YieldsError               = NULL;
Double_t*   fMesonM02TrueMergedOneGammaFromEtaYields                    = NULL;
Double_t*   fMesonM02TrueMergedOneGammaFromEtaYieldsError               = NULL;
Double_t*   fMesonM02TrueMergedOneElectronYields                        = NULL;
Double_t*   fMesonM02TrueMergedOneElectronYieldsError                   = NULL;
Double_t*   fMesonM02TrueMergedOneElectronFromPi0Yields                 = NULL;
Double_t*   fMesonM02TrueMergedOneElectronFromPi0YieldsError            = NULL;
Double_t*   fMesonM02TrueMergedOneElectronFromEtaYields                 = NULL;
Double_t*   fMesonM02TrueMergedOneElectronFromEtaYieldsError            = NULL;

Double_t*   fMesonM02TruePrimPi0Yields                                  = NULL;
Double_t*   fMesonM02TruePrimPi0YieldsError                             = NULL;

//****************************************************************************
//************************ correction histograms *****************************
//****************************************************************************
TH1D*       fDeltaPt                                                    = NULL;
TH1D*       fHistoMCAcceptancePt                                        = NULL;
TH1D*       fHistoTrueEffiMerged                                        = NULL;
TH1D*       fHistoTrueEffiPrimMeson                                     = NULL;
TH1D*       fHistoTruePurityMerged                                      = NULL;
TH1D*       fHistoTruePi0PurityMerged                                   = NULL;
TH1D*       fHistoTrueEtaPurityMerged                                   = NULL;

//****************************************************************************
//****************************************************************************
//****************************************************************************
TH1D*       fHistoYieldMesonM02                                         = NULL;
TH1D*       fHistoTrueYieldMergedM02                                    = NULL;
TH1D*       fHistoTrueYieldMergedPartConvM02                            = NULL;
TH1D*       fHistoTrueYieldMergedPartConvFromPi0M02                     = NULL;
TH1D*       fHistoTrueYieldMergedPartConvFromEtaM02                     = NULL;
TH1D*       fHistoTrueYieldMergedPureM02                                = NULL;
TH1D*       fHistoTrueYieldMergedPureFromPi0M02                         = NULL;
TH1D*       fHistoTrueYieldMergedPureFromEtaM02                         = NULL;
TH1D*       fHistoTrueYieldMergedOneGammaM02                            = NULL;
TH1D*       fHistoTrueYieldMergedOneGammaFromPi0M02                     = NULL;
TH1D*       fHistoTrueYieldMergedOneGammaFromEtaM02                     = NULL;
TH1D*       fHistoTrueYieldMergedOneElectronM02                         = NULL;
TH1D*       fHistoTrueYieldMergedOneElectronFromPi0M02                  = NULL;
TH1D*       fHistoTrueYieldMergedOneElectronFromEtaM02                  = NULL;
TH1D*       fHistoTrueYieldPi0M02                                       = NULL;
TH1D*       fHistoTrueYieldPi0DCM02                                     = NULL;
TH1D*       fHistoTrueYieldPi0GGM02                                     = NULL;
TH1D*       fHistoTrueYieldPi0DalitzM02                                 = NULL;
TH1D*       fHistoTrueYieldEtaM02                                       = NULL;
TH1D*       fHistoTrueYieldEtaDCM02                                     = NULL;
TH1D*       fHistoTrueYieldEtaGGM02                                     = NULL;
TH1D*       fHistoTrueYieldEtaDalitzM02                                 = NULL;
TH1D*       fHistoTrueYieldGammaM02                                     = NULL;
TH1D*       fHistoTrueYieldElectronM02                                  = NULL;
TH1D*       fHistoTrueYieldBGM02                                        = NULL;
TH1D*       fHistoTrueYieldPrimPi0M02                                   = NULL;

TH1D*       fHistoRatioTrueYieldEtaM02                                  = NULL;
TH1D*       fHistoRatioTrueYieldPi0M02                                  = NULL;
TH1D*       fHistoRatioTrueYieldGammaM02                                = NULL;
TH1D*       fHistoRatioTrueYieldElectronM02                             = NULL;
TH1D*       fHistoRatioPi0DCFrac                                        = NULL;
TH1D*       fHistoRatioPi0GGFrac                                        = NULL;
TH1D*       fHistoRatioPi0DalitzFrac                                    = NULL;
TH1D*       fHistoRatioEtaDCFrac                                        = NULL;
TH1D*       fHistoRatioEtaGGFrac                                        = NULL;
TH1D*       fHistoRatioEtaDalitzFrac                                    = NULL;
TH1D*       fHistoRatioMergedPureFracPi0                                = NULL;
TH1D*       fHistoRatioMergedPartConvFracPi0                            = NULL;
TH1D*       fHistoRatioMergedOneGammaFracPi0                            = NULL;
TH1D*       fHistoRatioMergedOneElectronFracPi0                         = NULL;
TH1D*       fHistoRatioMergedPureFracEta                                = NULL;
TH1D*       fHistoRatioMergedPartConvFracEta                            = NULL;
TH1D*       fHistoRatioMergedOneGammaFracEta                            = NULL;
TH1D*       fHistoRatioMergedOneElectronFracEta                         = NULL;
    


//****************************************************************************
//******************* MC input histograms ************************************
//****************************************************************************
TH1D*       fHistoMCMesonPt                                             = NULL;
TH1D*       fHistoMCMesonGGPt                                           = NULL;
TH1D*       fHistoMCMesonDalitzPt                                       = NULL;
TH1D*       fHistoMCMesonPtRebin                                        = NULL;
TH1D*       fHistoMCMesonPtWOWeights                                    = NULL;
TH1D*       fHistoMCMesonGGPtWOWeights                                  = NULL;
TH1D*       fHistoMCMesonDalitzPtWOWeights                              = NULL;
TH1D*       fHistoMCMesonWithinAccepPt                                  = NULL;
TH1D*       fHistoMCMesonGGWithinAccepPt                                = NULL;
TH1D*       fHistoMCMesonDalitzWithinAccepPt                            = NULL;
TH1D*       fHistoMCMesonWithinAccepPtRebin                             = NULL;

//****************************************************************************
//**************************** MC rec sec mesons  ****************************
//****************************************************************************
TH2F*       fHistoTrueClustersSecPi0PtM02[4]                            = { NULL, NULL, NULL, NULL };
TH1D**      fHistoTrueClusSecPi0M02PtBin[4]                             = { NULL, NULL, NULL, NULL};    

Double_t*   fMesonM02TrueSecPi0Yields[4]                                = { NULL, NULL, NULL, NULL};    
Double_t*   fMesonM02TrueSecPi0YieldsError[4]                           = { NULL, NULL, NULL, NULL};    
TH1D*       fHistoTrueYieldSecPi0M02[4]                                 = { NULL, NULL, NULL, NULL};    
TH1D*       fHistoTruePi0SecFrac[4]                                     = { NULL, NULL, NULL, NULL};
TH1D*       fHistoTrueEffiSecPi0[4]                                     = { NULL, NULL, NULL, NULL};

//****************************************************************************
//******************* Secondary correction histograms ************************
//****************************************************************************
TH2D*       fHistoMCSecPi0PtSource                                      = NULL;
TH2D*       fHistoMCSecPi0WithinAccepPtSource                           = NULL;
TH1D*       fHistoMCSecPi0Pt[4]                                         = { NULL, NULL, NULL, NULL};
TH1D*       fHistoMCSecPi0PtWAcc[4]                                     = { NULL, NULL, NULL, NULL};
TH1D*       fHistoMCSecPi0PtReb[4]                                      = { NULL, NULL, NULL, NULL};
TH1D*       fHistoMCSecPi0PtWAccReb[4]                                  = { NULL, NULL, NULL, NULL};
TH1D*       fHistoMCSecPi0AcceptPt[4]                                   = { NULL, NULL, NULL, NULL};

//*****************************************************************************
//************ Load secondary pion histograms from external file **************
//*****************************************************************************
Bool_t      fHaveToyMCInputForSec                                       = kFALSE;
TFile*      fFileToyMCInput[3]                                          = {NULL, NULL, NULL};
TH1D*       fHistoYieldToyMCSecInput[3]                                 = {NULL, NULL, NULL};
TH1D*       fHistoYieldToyMCSecInputReb[3]                              = {NULL, NULL, NULL};