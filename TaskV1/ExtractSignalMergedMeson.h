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
TString     fEventCutSelection                                          = "";
TString     fClusterCutSelection                                        = "";
TString     fClusterMergedCutSelection                                  = "";
TString     fMesonCutSelection                                          = "";
TString     fAdditionalName                                             = "";

Int_t       fCrysFitting                                                = 0;
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
TString     labelsBG[9]                                                 = { "PiCh", "P", "KCh", "N", "K0s", 
                                                                            "Lambda", "Mu", "K0l", "Rest"};

//****************************************************************************
//************************* MC object names **********************************
//****************************************************************************
TString     fObjectNameMCMesonAcc                                       = "";
TString     fObjectNameMCMeson                                          = "";
TString     fObjectNameMCMesonWOWeights                                 = "";
TString     fObjectNameTrueMergedM02                                    = "";
TString     fObjectNameTrueMergedInvMass                                = "";
TString     fObjectNameTrueFromPi0M02                                   = "";
TString     fObjectNameTrueFromPi0InvMass                               = "";
TString     fObjectNameTrueFromEtaM02                                   = "";
TString     fObjectNameTrueFromEtaInvMass                               = "";
TString     fObjectNameTrueMergedPartConvM02                            = "";
TString     fObjectNameTrueMergedPartConvLeadEM02                       = "";
TString     fObjectNameTrueMergedPartConvLeadEInvMass                   = "";
TString     fObjectNameTrueMergedPartConvInvMass                        = "";
TString     fObjectNameTrueClusPartConvFromPi0M02                       = "";
TString     fObjectNameTrueClusPartConvFromPi0InvMass                   = "";
TString     fObjectNameTrueClusPartConvFromEtaM02                       = "";
TString     fObjectNameTrueClusPartConvFromEtaInvMass                   = "";
TString     fObjectNameTrueClusBGM02                                    = "";
TString     fObjectNameTrueClusBGInvMass                                = "";
TString     fObjectNameTrueClusGammaM02                                 = "";
TString     fObjectNameTrueClusGammaInvMass                             = "";
TString     fObjectNameTrueClusElectronM02                              = "";
TString     fObjectNameTrueClusElectronInvMass                          = "";
TString     fObjectNameTrueClusBG_Source                                = "";
TString     fObjectNameTrueClusPrimMesonM02                             = "";
TString     fObjectNameTrueClusPrimMesonInvMass                         = "";
TString     fObjectNameTrueClusSecMesonM02                              = "";
TString     fObjectNameTrueClusSecMesonInvMass                          = "";
TString     fObjectNameTrueClusSecMesonFromK0sM02                       = "";
TString     fObjectNameTrueClusSecMesonFromK0sInvMass                   = "";
TString     fObjectNameTrueClusSecMesonFromLambdaM02                    = "";
TString     fObjectNameTrueClusSecMesonFromLambdaInvMass                = "";
TString     fObjectNameTrueClusPartConvPrimMesonM02                     = "";
TString     fObjectNameTrueClusPartConvPrimMesonInvMass                 = "";
TString     fObjectNameTrueClusPartConvSecMesonM02                      = "";
TString     fObjectNameTrueClusPartConvSecMesonInvMass                  = "";
TString     fObjectNameTrueClusPartConvSecMesonFromK0sM02               = "";
TString     fObjectNameTrueClusPartConvSecMesonFromK0sInvMass           = "";
TString     fObjectNameTrueClusPartConvSecMesonFromLambdaM02            = "";
TString     fObjectNameTrueClusPartConvSecMesonFromLambdaInvMass        = "";

//****************************************************************************
//******************************** Functions *********************************
//****************************************************************************
void CreatePtHistos();                                                                                      // Creat pt dependent histograms
void FillPtHistos();                                                                                        // Fill pt dependent histograms
void Initialize(TString setPi0, Int_t);                                                                     // Initialization of global variables depending on meson analysed 
void CalculateFWHM(TF1 *);                                                                                  // Calculation of FWHM
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
void FillMCInvMassHistosArray( TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F* );             // Fill invariant mass histograms for MC 
void FillMCPrimSecInvMassHistosArray( TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F* );             // Fill invariant mass histograms for MC for Pi0 Prim Sec 
void FillMCM02HistosArray( TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F* );                 // Fill M02 histograms for MC 
void FillMCPrimSecM02HistosArray( TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F*, TH2F* );                 // Fill M02 histograms for MC for Pi0 Prim Sec 
void FillMCBGSeparated(TH2F* );                                                                             // Fill MC BG separated 
void FillHistosArrayMC(TH1D*, TH1D * , TH1D * );                                                            // Rebin MC input histo's
void IntegrateHistoInvMass(TH1D*, Double_t* );                                                              // Integrate invariant mass histogram
void IntegrateHistoM02(TH1D*, Double_t* );                                                                  // Integrate M02 histogram
void CalculateMesonAcceptance();                                                                            // Calculation of meson acceptance
TH1D* CalculateMesonEfficiency(TH1D*, TH1D*,TString);                                                       // Calculation of meson efficiencies 
TH1D* CalculatePurity(TH1D*, TH1D*,TString);                                                                // Calculation of purity
TH1D* CalculateSecondaryFractions(TH1D*, TH1D*, TString );                                                  // Calculate fraction of secondaries
void SaveHistos(Int_t, TString, TString);                                                                   // Saving standard histograms to a file
void SaveCorrectionHistos(TString , TString);                                                               // Saving correction histograms to a file
void Delete();                                                                                              // Deleting all pointers
void SetCorrectMCHistogrammNames(TString);                                                                  // Setting correct histogram names

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
TH2F*       fHistoTrueClustersMergedInvMassPt                           = NULL;
TH2F*       fHistoTrueClustersPi0PtM02                                  = NULL;
TH2F*       fHistoTrueClustersPi0InvMassPt                              = NULL;
TH2F*       fHistoTrueClustersEtaPtM02                                  = NULL;
TH2F*       fHistoTrueClustersEtaInvMassPt                              = NULL;
TH2F*       fHistoTrueClustersBGPtM02                                   = NULL;
TH2F*       fHistoTrueClustersBGPtSource                                = NULL;
TH1D*       fHistoTrueClustersBGPt[9]                                   = { NULL, NULL, NULL, NULL, NULL,
                                                                            NULL, NULL, NULL, NULL};
TH2F*       fHistoTrueClustersBGInvMassPt                               = NULL;
TH2F*       fHistoTrueClustersGammaPtM02                                = NULL;
TH2F*       fHistoTrueClustersGammaInvMassPt                            = NULL;
TH2F*       fHistoTrueClustersElectronPtM02                             = NULL;
TH2F*       fHistoTrueClustersElectronInvMassPt                         = NULL;
TH2F*       fHistoTrueClusPartConvMergedPtM02                           = NULL;
TH2F*       fHistoTrueClusPartConvMergedInvMassPt                       = NULL;
TH2F*       fHistoTrueClusPartConvPi0PtM02                              = NULL;
TH2F*       fHistoTrueClusPartConvPi0InvMassPt                          = NULL;
TH2F*       fHistoTrueClusPartConvEtaPtM02                              = NULL;
TH2F*       fHistoTrueClusPartConvEtaInvMassPt                          = NULL;

TH2F*       fHistoTrueClustersPrimPi0PtM02                              = NULL;
TH2F*       fHistoTrueClustersPrimPi0InvMassPt                          = NULL;
TH2F*       fHistoTrueClustersSecPi0PtM02                               = NULL;
TH2F*       fHistoTrueClustersSecPi0InvMassPt                           = NULL;
TH2F*       fHistoTrueClustersSecPi0FK0sPtM02                           = NULL;
TH2F*       fHistoTrueClustersSecPi0FK0sInvMassPt                       = NULL;
TH2F*       fHistoTrueClustersSecPi0FLambdaPtM02                        = NULL;
TH2F*       fHistoTrueClustersSecPi0FLambdaInvMassPt                    = NULL;
TH2F*       fHistoTrueClusPartConvPrimPi0PtM02                          = NULL;
TH2F*       fHistoTrueClusPartConvPrimPi0InvMassPt                      = NULL;
TH2F*       fHistoTrueClusPartConvSecPi0PtM02                           = NULL;
TH2F*       fHistoTrueClusPartConvSecPi0InvMassPt                       = NULL;
TH2F*       fHistoTrueClusPartConvSecPi0FK0sPtM02                       = NULL;
TH2F*       fHistoTrueClusPartConvSecPi0FK0sInvMassPt                   = NULL;
TH2F*       fHistoTrueClusPartConvSecPi0FLambdaPtM02                    = NULL;
TH2F*       fHistoTrueClusPartConvSecPi0FLambdaInvMassPt                = NULL;

//****************************************************************************
//***************************** histos in Pt Bins ****************************
//****************************************************************************
TH1D**      fHistoTrueClusMergedInvMassPtBin                            = NULL;
TH1D**      fHistoTrueClusPi0InvMassPtBin                               = NULL;
TH1D**      fHistoTrueClusEtaInvMassPtBin                               = NULL;
TH1D**      fHistoTrueClusGammaInvMassPtBin                             = NULL;
TH1D**      fHistoTrueClusElectronInvMassPtBin                          = NULL;
TH1D**      fHistoTrueClusBGInvMassPtBin                                = NULL;
TH1D**      fHistoTrueClusMergedM02PtBin                                = NULL;
TH1D**      fHistoTrueClusPi0M02PtBin                                   = NULL;
TH1D**      fHistoTrueClusEtaM02PtBin                                   = NULL;
TH1D**      fHistoTrueClusGammaM02PtBin                                 = NULL;
TH1D**      fHistoTrueClusElectronM02PtBin                              = NULL;
TH1D**      fHistoTrueClusBGM02PtBin                                    = NULL;
TH1D**      fHistoTrueClusPartConvMergedInvMassPtBin                    = NULL;
TH1D**      fHistoTrueClusPartConvPi0InvMassPtBin                       = NULL;
TH1D**      fHistoTrueClusPartConvEtaInvMassPtBin                       = NULL;
TH1D**      fHistoTrueClusPartConvMergedM02PtBin                        = NULL;
TH1D**      fHistoTrueClusPartConvPi0M02PtBin                           = NULL;
TH1D**      fHistoTrueClusPartConvEtaM02PtBin                           = NULL;
TH1D**      fHistoTrueClusFullMergedInvMassPtBin                        = NULL;
TH1D**      fHistoTrueClusFullPi0InvMassPtBin                           = NULL;
TH1D**      fHistoTrueClusFullEtaInvMassPtBin                           = NULL;
TH1D**      fHistoTrueClusFullMergedM02PtBin                            = NULL;
TH1D**      fHistoTrueClusFullPi0M02PtBin                               = NULL;
TH1D**      fHistoTrueClusFullEtaM02PtBin                               = NULL;

TH1D**      fHistoTrueClusPrimPi0InvMassPtBin                           = NULL;
TH1D**      fHistoTrueClusSecPi0InvMassPtBin                            = NULL;
TH1D**      fHistoTrueClusSecPi0FK0sInvMassPtBin                        = NULL;
TH1D**      fHistoTrueClusSecPi0FLambdaInvMassPtBin                     = NULL;
TH1D**      fHistoTrueClusPartConvPrimPi0InvMassPtBin                   = NULL;
TH1D**      fHistoTrueClusPartConvSecPi0InvMassPtBin                    = NULL;
TH1D**      fHistoTrueClusPartConvSecPi0FK0sInvMassPtBin                = NULL;
TH1D**      fHistoTrueClusPartConvSecPi0FLambdaInvMassPtBin             = NULL;
TH1D**      fHistoTrueClusFullPrimPi0InvMassPtBin                       = NULL;
TH1D**      fHistoTrueClusFullSecPi0InvMassPtBin                        = NULL;
TH1D**      fHistoTrueClusFullSecPi0FK0sInvMassPtBin                    = NULL;
TH1D**      fHistoTrueClusFullSecPi0FLambdaInvMassPtBin                 = NULL;
TH1D**      fHistoTrueClusPrimPi0M02PtBin                               = NULL;
TH1D**      fHistoTrueClusSecPi0M02PtBin                                = NULL;
TH1D**      fHistoTrueClusSecPi0FK0sM02PtBin                            = NULL;
TH1D**      fHistoTrueClusSecPi0FLambdaM02PtBin                         = NULL;
TH1D**      fHistoTrueClusPartConvPrimPi0M02PtBin                       = NULL;
TH1D**      fHistoTrueClusPartConvSecPi0M02PtBin                        = NULL;
TH1D**      fHistoTrueClusPartConvSecPi0FK0sM02PtBin                    = NULL;
TH1D**      fHistoTrueClusPartConvSecPi0FLambdaM02PtBin                 = NULL;
TH1D**      fHistoTrueClusFullPrimPi0M02PtBin                           = NULL;
TH1D**      fHistoTrueClusFullSecPi0M02PtBin                            = NULL;
TH1D**      fHistoTrueClusFullSecPi0FK0sM02PtBin                        = NULL;
TH1D**      fHistoTrueClusFullSecPi0FLambdaM02PtBin                     = NULL;

//****************************************************************************
//************************** Yields vs Pt ************************************
//****************************************************************************
Double_t*   fMesonM02Yields                                             = NULL;
Double_t*   fMesonM02TrueMergedYields                                   = NULL;
Double_t*   fMesonM02TruePi0Yields                                      = NULL;
Double_t*   fMesonM02TrueEtaYields                                      = NULL;
Double_t*   fMesonM02TrueGammaYields                                    = NULL;
Double_t*   fMesonM02TrueElectronYields                                 = NULL;
Double_t*   fMesonM02TrueBGYields                                       = NULL;
Double_t*   fMesonM02YieldsError                                        = NULL;
Double_t*   fMesonM02TrueMergedYieldsError                              = NULL;
Double_t*   fMesonM02TruePi0YieldsError                                 = NULL;
Double_t*   fMesonM02TrueEtaYieldsError                                 = NULL;
Double_t*   fMesonM02TrueGammaYieldsError                               = NULL;
Double_t*   fMesonM02TrueElectronYieldsError                            = NULL;
Double_t*   fMesonM02TrueBGYieldsError                                  = NULL;
Double_t*   fMesonM02TrueMergedPartConvYields                           = NULL;
Double_t*   fMesonM02TruePi0PartConvYields                              = NULL;
Double_t*   fMesonM02TrueEtaPartConvYields                              = NULL;
Double_t*   fMesonM02TrueMergedPartConvYieldsError                      = NULL;
Double_t*   fMesonM02TruePi0PartConvYieldsError                         = NULL;
Double_t*   fMesonM02TrueEtaPartConvYieldsError                         = NULL;
Double_t*   fMesonM02TrueMergedFullYields                               = NULL;
Double_t*   fMesonM02TruePi0FullYields                                  = NULL;
Double_t*   fMesonM02TrueEtaFullYields                                  = NULL;
Double_t*   fMesonM02TrueMergedFullYieldsError                          = NULL;
Double_t*   fMesonM02TruePi0FullYieldsError                             = NULL;
Double_t*   fMesonM02TrueEtaFullYieldsError                             = NULL;

Double_t*   fMesonM02TruePrimPi0Yields                                  = NULL;
Double_t*   fMesonM02TrueSecPi0Yields                                   = NULL;
Double_t*   fMesonM02TrueSecPi0FK0sYields                               = NULL;
Double_t*   fMesonM02TrueSecPi0FLambdaYields                            = NULL;
Double_t*   fMesonM02TruePrimPi0YieldsError                             = NULL;
Double_t*   fMesonM02TrueSecPi0YieldsError                              = NULL;
Double_t*   fMesonM02TrueSecPi0FK0sYieldsError                          = NULL;
Double_t*   fMesonM02TrueSecPi0FLambdaYieldsError                       = NULL;
Double_t*   fMesonM02TruePrimPi0PartConvYields                          = NULL;
Double_t*   fMesonM02TrueSecPi0PartConvYields                           = NULL;
Double_t*   fMesonM02TrueSecPi0FK0sPartConvYields                       = NULL;
Double_t*   fMesonM02TrueSecPi0FLambdaPartConvYields                    = NULL;
Double_t*   fMesonM02TruePrimPi0PartConvYieldsError                     = NULL;
Double_t*   fMesonM02TrueSecPi0PartConvYieldsError                      = NULL;
Double_t*   fMesonM02TrueSecPi0FK0sPartConvYieldsError                  = NULL;
Double_t*   fMesonM02TrueSecPi0FLambdaPartConvYieldsError               = NULL;
Double_t*   fMesonM02TruePrimPi0FullYields                              = NULL;
Double_t*   fMesonM02TrueSecPi0FullYields                               = NULL;
Double_t*   fMesonM02TrueSecPi0FK0sFullYields                           = NULL;
Double_t*   fMesonM02TrueSecPi0FLambdaFullYields                        = NULL;
Double_t*   fMesonM02TruePrimPi0FullYieldsError                         = NULL;
Double_t*   fMesonM02TrueSecPi0FullYieldsError                          = NULL;
Double_t*   fMesonM02TrueSecPi0FK0sFullYieldsError                      = NULL;
Double_t*   fMesonM02TrueSecPi0FLambdaFullYieldsError                   = NULL;

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
TH1D*       fHistoTruePi0SecFrac                                        = NULL;
TH1D*       fHistoTruePi0SecFracFK0S                                    = NULL;
TH1D*       fHistoTruePi0SecFracFLambda                                 = NULL;

//****************************************************************************
//****************************************************************************
//****************************************************************************
TH1D*       fHistoYieldMesonM02                                         = NULL;
TH1D*       fHistoTrueYieldMergedM02                                    = NULL;
TH1D*       fHistoTrueYieldMergedPartConvM02                            = NULL;
TH1D*       fHistoTrueYieldMergedFullM02                                = NULL;
TH1D*       fHistoTrueYieldPi0M02                                       = NULL;
TH1D*       fHistoTrueYieldPi0PartConvM02                               = NULL;
TH1D*       fHistoTrueYieldPi0FullM02                                   = NULL;
TH1D*       fHistoTrueYieldEtaM02                                       = NULL;
TH1D*       fHistoTrueYieldEtaPartConvM02                               = NULL;
TH1D*       fHistoTrueYieldEtaFullM02                                   = NULL;
TH1D*       fHistoTrueYieldGammaM02                                     = NULL;
TH1D*       fHistoTrueYieldElectronM02                                  = NULL;
TH1D*       fHistoTrueYieldBGM02                                        = NULL;
TH1D*       fHistoTrueYieldPrimPi0M02                                   = NULL;
TH1D*       fHistoTrueYieldPrimPi0PartConvM02                           = NULL;
TH1D*       fHistoTrueYieldPrimPi0FullM02                               = NULL;
TH1D*       fHistoTrueYieldSecPi0M02                                    = NULL;
TH1D*       fHistoTrueYieldSecPi0PartConvM02                            = NULL;
TH1D*       fHistoTrueYieldSecPi0FullM02                                = NULL;
TH1D*       fHistoTrueYieldSecPi0FK0sM02                                = NULL;
TH1D*       fHistoTrueYieldSecPi0FK0sPartConvM02                        = NULL;
TH1D*       fHistoTrueYieldSecPi0FK0sFullM02                            = NULL;
TH1D*       fHistoTrueYieldSecPi0FLambdaM02                             = NULL;
TH1D*       fHistoTrueYieldSecPi0FLambdaPartConvM02                     = NULL;
TH1D*       fHistoTrueYieldSecPi0FLambdaFullM02                         = NULL;

//****************************************************************************
//******************* MC input histograms ************************************
//****************************************************************************
TH1D*       fHistoMCMesonPt                                             = NULL;
TH1D*       fHistoMCMesonPtRebin                                        = NULL;
TH1D*       fHistoMCMesonPtWOWeights                                    = NULL;
TH1D*       fHistoMCMesonWithinAccepPt                                  = NULL;
TH1D*       fHistoMCMesonWithinAccepPtRebin                             = NULL;

