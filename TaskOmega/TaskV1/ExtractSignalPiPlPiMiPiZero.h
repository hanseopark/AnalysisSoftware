// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

// Double, Int, etc.
TDatime	now;
fstream fFileErrLog;
fstream	fFileDataLog;

TString fTypeCutSelection                                   = "";
TString fEventCutSelection                                  = "";
TString fGammaCutSelection                                  = "";
TString fClusterCutSelection                                = "";
TString fPionCutSelection                                   = "";
TString fNeutralPionCutSelection                            = "";
TString fMesonCutSelection                                  = "";

Int_t fCrysFitting                                          = 0;
Int_t fIsMC                                                 = 0;
Int_t fMesonId                                              = 0;
Float_t fBackgroundMultNumber;

Double_t fYMaxMeson                                         = 0.9;
Double_t fMesonMassExpect                                   = 0;   // expected meson mass
Double_t fNEvents                                           = 0;

Double_t *fPeakRange;
Double_t *fIntFixedRange;
Double_t *fFitRange;
Double_t *fBGFitRange                                       = NULL;
Double_t *fMesonIntRange                                    = NULL;

// TSring
TString fTextCent;
TString fEnergyFlag;
TString fPrefix;
TString fPrefix2;
TString fPeriodFlag;
TString ftextDayth;
TString fdate;
TString fCollisionSystem;
TString fTextMeasurement;
TString fDetectionProcess;
TString fCutSelection;
TString fBackgroundMultCutNumber;
Int_t   fMode;


TString fdirectphoton;
TString fThesis                                             = "";

// Output Files
TFile*  fOutput                                             = NULL;
TFile*  fOutput1                                            = NULL;
TFile*  fOutput2                                            = NULL;

//TH1D
TH1D*   fBckNorm                                            = NULL;
TH1D*   fSignal                                             = NULL;
TH1D*   fRatioSB                                            = NULL;
TH1D*   fFittingHistMidPtSignal                             = NULL;
TH1D*   fFittingHistMidPtBackground                         = NULL;
TH1D*   fFittingHistMidPtSignalSub                          = NULL;
TH1D*   fMesonFullPtSignal                                  = NULL;
TH1D*   fMesonFullPtBackground                              = NULL;
TH1D*   fOmegaFullPtSignal                                  = NULL;
TH2D*   fOmegaFullPtSignalR2                                = NULL;
TH1D*   fOmegaFullPtSignalR1                                = NULL;
TH1D*   fOmegaFullPtBackNorm                                = NULL;
TH1D*   fOmegaFullPtBack2                                   = NULL;
TH1D*   fOmegaFullPtBack3                                   = NULL;
TH1D*   fOmegaFullPtBack4                                   = NULL;
TH1D*   fRelation                                           = NULL;
TH1D*   fMesonFullPtBackNorm                                = NULL;
TH1D*   fHistoMotherZMProj                                  = NULL;
TH1D*   fHistoBckZMProj                                     = NULL;
TH1D*   fNumberOfGoodESDTracks                              = NULL;
TH1D*   fEventQuality                                       = NULL;
TH1D*   fHistoMappingBackNormInvMass                        = NULL;
TH1D*   fHistoMappingSignalInvMass                          = NULL;
TH1D*   fHistoYieldMeson[6]                                 = {NULL, NULL, NULL, NULL, NULL, NULL};
TH1D*   fHistoYieldMesonBackFit                             = NULL;
TH1D*   fHistoYieldMesonPerEvent[6]                         = {NULL, NULL, NULL, NULL, NULL, NULL};
TH1D*   fHistoYieldMesonPerEventBackFit                     = NULL;
TH1D*   fHistoSBdefaultMeson[6]                             = { NULL, NULL, NULL, NULL, NULL, NULL};
TH1D*   fHistoSigndefaultMeson[6]                           = { NULL, NULL, NULL, NULL, NULL, NULL};
TH1D*   fHistoTrueSignMeson                                 = NULL;
TH1D*   fHistoTrueSBMeson                                   = NULL;
TH1D*   fHistoMassPosition=									NULL;
TH1D*   fHistoMassMeson=									NULL;
TH1D*   fHistoWidthMeson=									NULL;
TH1D*   fHistoFWHMMeson=									NULL;
TH1D*   fHistoTrueMassMeson=								NULL;
TH1D*   fHistoTrueMassMesonReweighted=						NULL;
TH1D*   fHistoTrueFWHMMeson=								NULL;
TH1D*   fHistoTrueFWHMMesonReweighted=						NULL;
TH1D*   fDeltaPt=											NULL;
TH1D*   fHistoMassMesonLeft=								NULL;
TH1D*   fHistoWidthMesonLeft=								NULL;
TH1D*   fHistoFWHMMesonLeft=								NULL;
TH1D*   fHistoMCMesonPtWithinAcceptance=					NULL;
TH1D*   fHistoMCMesonPt=									NULL;
TH1D*   fHistoMCMesonPtWOWeights=							NULL;
TH1D*   fHistoMCMesonPtWeights=								NULL;
TH1D*   fHistoMCMesonWithinAccepPt=							NULL; // Proper bins in Pt
TH1D*   fHistoMCMesonPt1=									NULL; // Proper bins in Pt
TH1D*   fHistoYieldTrueMeson[3]                             = {NULL, NULL, NULL};
TH1D*   fHistoYieldTrueMesonReweighted[3]                   = {NULL, NULL, NULL};

TH1D*   fHistoMCMesonAcceptPt                               = NULL;
TH1D*   fHistoMCMesonEffiPt                                 = NULL;
TH1D*   fHistoMCMesonEffiFitPt                              = NULL;
TH1D*   fHistoMonteMesonEffiPt[6]                           = {NULL, NULL, NULL, NULL, NULL, NULL};
TH1D*   fHistoMonteMesonEffiBackFitPt                       = NULL;
TH1D*   fHistoMonteMesonEffiMCAll                           = NULL;
TH1D*   fHistoMCTrueMesonEffiPt[3]                          = {NULL, NULL, NULL};
TH1D*   fHistoMCTrueMesonEffiPtReweighted[3]                = {NULL, NULL, NULL};
TH1D*   fHistoMCTrueMesonEffiFitPt[3]                       =  {NULL, NULL, NULL};

TH1D** 	fHistoMappingTrueMesonInvMassPtBins=				NULL;
TH1D**	fHistoMappingTrueMesonInvMassPtReweightedBins=		NULL;
TH1D** 	fHistoMappingTrueGGBckInvMassPtBins=				NULL;
TH1D** 	fHistoMappingTrueContBckInvMassPtBins=				NULL;
TH1D** 	fHistoMappingTrueAllBckInvMassPtBins=				NULL;
TH1D** 	fHistoMappingGGInvMassPtBin=						NULL;
TH1D** 	fHistoMappingGGInvMassBackFitPtBin=					NULL;
TH1D** 	fHistoMappingGGInvMassBackFitWithoutSignalPtBin=	NULL;
TH1D** 	fHistoMappingBackInvMassPtBin=						NULL;
TH1D** 	fHistoMappingBackNormInvMassPtBin=					NULL;
//TH1D** 	fRatio=                                             NULL;
TH1D** 	fHistoMappingRatioSBInvMassPtBin=					NULL;
TH1D** 	fHistoMappingSignalInvMassPtBin=					NULL;
TH1D** 	fHistoMappingPeakPosInvMassPtBin=					NULL;
// Histograms for normalization on the left of the peak
TH1D** 	fHistoMappingBackNormInvMassLeftPtBin=				NULL;
TH1D** 	fHistoMappingSignalInvMassLeftPtBin=				NULL;
TH1D** 	fHistoMappingTrueMesonCaloPhotonInvMassPtBins =					NULL;
TH1D** 	fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins =				NULL;
TH1D** 	fHistoMappingTrueMesonCaloElectronInvMassPtBins =				NULL;
TH1D** 	fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins =			NULL;
TH1D** 	fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins =	NULL;
TH1D** 	fHistoMappingTrueMesonCaloEMNonLeadingInvMassPtBins =			NULL;

TF1**  fBackgroundFitPol=NULL;
TF1 * 	fFitReco=											NULL;
TF1 * 	fFitGausExp=										NULL;
TF1 * 	fFitLinearBck=										NULL;
TF1 * 	fFitLinearBckOut=									NULL;
TF1 * 	fFitSignalInvMassMidPt=								NULL;

Double_t    fYields                                         = 0;
Double_t    fYieldsError                                    = 0;
Double_t    fYieldsFunc                                     = 0;
Double_t 	fFWHMFunc;
Double_t 	fFWHMFuncError;
Double_t 	fIntLinearBck;
Double_t 	fIntLinearBckError;
Double_t 	fIntLinearBckOut;
Double_t 	fIntLinearBckErrorOut;

Float_t 	pictDrawingCoordinatesFWHM[9] = 				{0.6, 0.8, 0.30, 0.04, 0.15,0.7, 0.1, 0.035,0};
Int_t 	fNRebinGlobal = 									2;

// common meson analysis variables (So it is easy to switch between the pi0 and the eta without changing the code)
Double_t 	*fBGFitRangeLeft = 								NULL;
Double_t 	*fMesonPlotRange = 								NULL;

// Names
TString     nameIntRange[6]                                             = {"", "Wide", "Narrow", "Left", "LeftWide", "LeftNarrow"};
TString     nameIntBck[6]                                               = {"Minpol2","Pluspol2","Minexp","Plusexp","Minexp2","Plusexp2"};
TString     nameIntBckRatios[6]                                         = {"RatioMinpol2","RatioPluspol2","RatioMinexp","RatioPlusexp","RatioMinexp2","RatioPlusexp2"};
TString     nameIntBckResult[3]                                         = {"pol2_normal","exp_normal","exp2_normal"};


// Added because of change of how integration ranges are stored
Double_t    *fMesonIntDeltaRange =                          NULL;
Double_t    *fMesonIntDeltaRangeWide =                      NULL;
Double_t    *fMesonIntDeltaRangeNarrow =                    NULL;
//
Double_t 	*fMesonIntRangeWide = 							NULL;
Double_t 	*fMesonIntRangeNarrow = 						NULL;
Double_t    *fMesonCurIntRange[6] =                         {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t 	*fMesonCurIntRangeBackFit = 					NULL;
Double_t 	*fMesonTrueIntRange[3] =						{NULL, NULL, NULL};
Double_t    *fMesonTrueIntReweightedRange[3] =				{NULL, NULL, NULL};

Double_t 	*fMesonMassRange = 								NULL;
Double_t    *fMesonMassPlotRange =                          NULL; // not implemented yet
Double_t 	*fMesonFitRange = 								NULL;
Double_t 	*fMesonFitRangeWithoutPeak = 					NULL;

Double_t 	fMesonWidthExpect=								0;
Double_t 	fMesonWidthExpectMC=							0;
Double_t 	fMesonLambdaTail=								0;
Double_t 	fMesonLambdaTailMC=								0;
Double_t 	*fMesonWidthRange = 							NULL;
Double_t 	*fMesonWidthRangeMC = 							NULL;
Double_t 	*fMesonLambdaTailRange = 						NULL;
Double_t 	*fMidPt = 										NULL;
Double_t 	*fFullPt = 										NULL;
// end common meson analysis variables

//Background histograms in different M and Z bins

TH2D* fHistoMotherZM;
TH2D* fHistoBckZM;
TH2D* fPi0ResponseMatrix;
TH2D* fPi0ResponseMatrixRebin;
TH2D* fPeakPosAlpha01;

Double_t 	fScalingFactorBck[7][4];
Double_t 	fCBAlpha;
Double_t 	fCBn;

///////////////////////////
TString 	fNameHistoGG;
TString 	fNameHistoBack;
TString 	fNameHistoPP;
TString 	fNameHistoBackNormOut;
TString 	fNameHistoBackNormLeftOut;
TString 	fNameHistoTrue;
TString 	fNameHistoTrueGGBck;
TString 	fNameHistoTrueContBck;
TString 	fNameHistoTrueAllBck;
TString 	fNameHistoEffi;
TString 	fNameHistoFrac;
TString 	fNameHistoMotherZM;
TString 	fNameHistoBckZM;
TString 	fNameFitSignalPos;

Double_t* fGGYields[6] =                                    {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fBckYields[6] =                                   {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fTotalBckYields[6]= 								{NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYields[6]=                                  {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsBackFit= 							     NULL;
Double_t* fMesonTrueYields[3]= 								{NULL, NULL, NULL};
Double_t* fMesonTrueYieldsReweighted[3]=                    {NULL, NULL, NULL};
Double_t* fMesonYieldsFunc[6] =                             {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsResidualBckFunc[6]= 					{NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsResidualBckFuncBackFit= 				NULL;
Double_t* fMesonYieldsCorResidualBckFunc[6]=				{NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsCorResidualBckFuncBackFit=			NULL;
Double_t* fMesonYieldsPerEvent[6]= 							{NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsPerEventBackFit= 						NULL;
Double_t* fMesonMass= 										NULL;
Double_t* fMesonMassBackFit= 								NULL;
Double_t* fMesonWidth= 										NULL;
Double_t* fMesonWidthBackFit= 								NULL;
//Double_t* fMesonSB= 										NULL;
Double_t* fMesonSBdefault[3]                                = { NULL, NULL, NULL};
Double_t* fMesonSigndefault[3]                              = { NULL, NULL, NULL};
Double_t* fMesonSBdefaultError[3]                           = { NULL, NULL, NULL};
Double_t* fMesonSigndefaultError[3]                         = { NULL, NULL, NULL};
Double_t* fMesonTrueSB= 									NULL;
Double_t* fMesonTrueSign= 									NULL;
//Double_t* fMesonSign= 										NULL;
Double_t* fMesonFWHM= 										NULL;
Double_t* fMesonTrueMass= 									NULL;
Double_t* fMesonTrueMassCaloPhoton= 						NULL;
Double_t* fMesonTrueMassCaloElectron= 						NULL;
Double_t* fMesonTrueMassCaloConvPhoton=						NULL;
Double_t* fMesonTrueMassCaloMergedCluster= 					NULL;
Double_t* fMesonTrueMassCaloMergedClusterPartConv= 			NULL;
Double_t* fMesonTrueMassCaloEMNonLeading= 					NULL;
Double_t* fMesonTrueMassReweighted=							NULL;
Double_t* fMesonTrueFWHM= 									NULL;
Double_t* fMesonTrueFWHMCaloPhoton= 						NULL;
Double_t* fMesonTrueFWHMCaloElectron= 						NULL;
Double_t* fMesonTrueFWHMCaloConvPhoton=						NULL;
Double_t* fMesonTrueFWHMCaloMergedCluster= 					NULL;
Double_t* fMesonTrueFWHMCaloMergedClusterPartConv= 			NULL;
Double_t* fMesonTrueFWHMCaloEMNonLeading= 					NULL;
Double_t* fMesonTrueFWHMReweighted=							NULL;
Double_t* fMesonFWHMAlpha01= 								NULL;

// Normalization at the left of the peak
Double_t* fMesonMassLeft=									NULL;
Double_t* fMesonWidthLeft=									NULL;
Double_t* fMesonFWHMLeft=									NULL;
Double_t* fMesonMassLeftError=                              NULL;
Double_t* fMesonFWHMLeftError=                              NULL;
Double_t* fMesonWidthLeftError=                              NULL;

Double_t fScaleFac  =                                         1.;

Double_t* fGGYieldsError[6]                             = {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fBckYieldsError[6] =                           {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsError[6]=							   {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsBackFitError=						    NULL;
Double_t* fMesonTrueYieldsError[3]=                        {NULL, NULL, NULL};
Double_t* fMesonTrueYieldsReweightedError[3]=              {NULL, NULL, NULL};
Double_t* fMesonYieldsFuncError[6]=							{NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsResidualBckFuncError[6]=              {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsResidualBckFuncBackFitError=			NULL;
Double_t* fMesonYieldsCorResidualBckFuncError[6]=			{NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsCorResidualBckFuncBackFitError=		NULL;
Double_t* fMesonYieldsPerEventError[6]=						{NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsPerEventBackFitError=					NULL;
Double_t* fMesonMassError=									NULL;
Double_t* fMesonMassBackFitError=							NULL;
Double_t* fMesonWidthError=									NULL;
Double_t* fMesonWidthBackFitError=							NULL;
Double_t* fMesonTrueMassError=								NULL;
Double_t* fMesonTrueMassReweightedError=					NULL;
Double_t* fMesonTrueFWHMReweightedError=					NULL;
Double_t* fMesonTrueFWHMError=								NULL;

Double_t* fMesonSBError=									NULL;
Double_t* fMesonSignError=									NULL;
Double_t* fMesonTrueSBError= 								NULL;
Double_t* fMesonTrueSignError=								NULL;
Double_t* fMesonFWHMError=									NULL;

Double_t* fTotalBckYieldsError[6]=							{NULL, NULL, NULL, NULL, NULL, NULL};

TF1** 	fFitSignalInvMassPtBin=								NULL;
TF1** 	fFitSignalInvMassBackFitPtBin=						NULL;
TF1** 	fFitTrueSignalInvMassPtBin=							NULL;
TF1**	fFitTrueSignalInvMassPtReweightedBin=				NULL;
TF1** 	fFitBckInvMassPtBin=								NULL;
TF1** 	fFitBckInvMassBackFitPtBin=							NULL;
TF1** 	fFitRatioInvMassPtBin=								NULL;
TF1**	fFitWithPol2ForBG = 								NULL;

TF1** 	fFitInvMassLeftPtBin=								NULL;
TF1** 	fFitSignalPeakPosInvMassLeftPtBin=					NULL;
TF1** 	fFitBckInvMassLeftPtBin=							NULL;

TH2D* 	fHistoTrueMesonInvMassVSPt=							NULL;
TH2D*	fHistoTrueMesonInvMassVSPtWOWeights=				NULL;
TH2D*	fHistoTrueMesonInvMassVSPtReweighted=				NULL;
TProfile2D*	fProfileTrueMesonInvMassVSPtWeights=		 	NULL;

TH2D*	fGammaGammaInvMassVSPt=								NULL;
TH2D*   hist_bck2=                                          NULL;
TH2D*   hist_bck3=                                          NULL;
TH2D*   hist_bck4=                                          NULL;
TH2D*	fBckInvMassVSPt=									NULL;
TH2D*	fBckInvMassVSPt2=									NULL;
TH2D*	fBckInvMassVSPt3=									NULL;
TH2D*	fBckInvMassVSPt4=									NULL;

TH2F**	fHistoWeightsBGZbinVsMbin = 						0x0;
TH2F**	fHistoFillPerEventBGZbinVsMbin = 					0x0;


TString ObjectNameTrue;
TString ObjectNameTrueWOWeights;
TString ObjectNameProfileWeights;
TString ObjectNameMCPi0Acc;
TString ObjectNameMCEtaAcc;
TString ObjectNameMCOmegaAcc;
TString ObjectNameMCPi0;
TString ObjectNameMCPi0WOWeights;
TString ObjectNameMCEta;
TString ObjectNameMCEtaWOWeights;
TString ObjectNameMCOmega;
TString ObjectNameMCOmegaWOWeights;


Bool_t fAdvancedMesonQA = kFALSE;
Bool_t fEstimateTrainPileUp = kFALSE;

void ProcessEM(TH1D*,TH1D*,Double_t *);
void ProcessEMLeftRight(    TH1D* , TH1D*, Double_t*, Double_t* );

void ProcessBckFitSubtraction(TH1D *fGammaGamma, Int_t i, Double_t * fPeakRange, Double_t *fFitRange, TString, TString, TString, TString);
void ProcessRatioSignalBackground(TH1D* , TH1D* );
void FillMassHistosArray(TH2D*);
void FillMassMCTrueMesonHistosArray(TH2D*);
void FillMassMCTrueReweightedMesonHistosArray(TH2D*);
void CreatePtHistos();
void FillPtHistos();
void FitSubtractedInvMassInPtBins(TH1D * , Double_t *fMesonIntDeltaRangeFit, Int_t, Bool_t );  // Fits the Invariant Mass histos with a given function
void FitTrueInvMassInPtBins(TH1D * , Double_t *fMesonIntDeltaRangeFit, Int_t, Bool_t);  // Fits the Invariant Mass histos with a given function
void FitCBSubtractedInvMassInPtBins(TH1D * , Double_t *fMesonIntDeltaRangeFit, Int_t, Bool_t , TString );  // Fits the Invariant Mass histos with a given function;
void FitWithPol2ForBG(TH1D*, Double_t*fMesonFitRangeCur, Int_t, Bool_t);
void ProduceBckProperWeighting(TList*, TList* );
void ProduceBckWithoutWeighting(TH2D *);
void IntegrateHistoInvMassStream(TH1D * , Double_t *);
void IntegrateHistoInvMass(TH1D * , Double_t *);
void IntegrateFitFunc(TF1 * , TH1D *, Double_t *);
void FillHistosArrayMC(TH1D* , TH1D*, TH1D* ,TString);
void CalculateMesonAcceptance();
TH1D* CalculateMesonEfficiency(TH1D*, TH1D*, TString);
void SaveHistos(Int_t, TString, TString);
void SaveCorrectionHistos(TString , TString);
void Initialize(TString setPi0, Int_t);
void CalculateFWHM(TF1 *);
Double_t CrystalBallBck(Double_t *,Double_t *);
Double_t CrystalBall(Double_t *,Double_t *);
void Delete();
void SetCorrectMCHistogrammNames();

TString centralityString = 										"";

void InitializeWindows(TString setPi0, Int_t mode, TString trigger, Int_t triggerSet = -1){

    fPeakRange                  = new Double_t[2];
    fIntFixedRange              = new Double_t[2]; // not yet implemented
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
    fMidPt                      = new Double_t[2];
    fMesonWidthRangeMC          = new Double_t[2];
    fFullPt                     = new Double_t[2];
    //****************************************************************************************************
    // Initialization for Omega meson
    //****************************************************************************************************
    if (setPi0.CompareTo("Omega") == 0){

        // set meson ID according to PDG
        fMesonId                  = 223;

        // set medium pt range (currently for all modes the same)
        if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45){
           fMidPt[0]                   = 1.5;
           fMidPt[1]                   = 2.5;
         }

        // Initialize peak range
        if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45){
             fPeakRange[0]             = 0.75;
             fPeakRange[1]             = 0.81;
         }

        // Initialze fit range
        if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45){
            fFitRange[0]             = 0.615;
            fFitRange[1]             = 0.89;
            fIntFixedRange[0]        = 0.615; // not yet implemented
            fIntFixedRange[1]        = 0.89;  // not yet implemented
        }

        // Initialize default BG fit range right & left
        if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45){
            fBGFitRange[0]             = 0.82;
            fBGFitRange[1]             = 0.89;
            fBGFitRangeLeft[0]         = 0.615;
            fBGFitRangeLeft[1]         = 0.69;
        }

        // Initialize default Plot range for meson
        if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45){
            fMesonPlotRange[0]         = 0.75;
            fMesonPlotRange[1]         = 0.79;
        }

        // Initialize default Plot default integration ranges
        if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45){
            fMesonIntDeltaRange[0]      = -0.04;
            fMesonIntDeltaRange[1]      =  0.05;
            fMesonIntDeltaRangeWide[0]  = -0.06;
            fMesonIntDeltaRangeWide[1]  = 0.05;
            fMesonIntDeltaRangeNarrow[0]= -0.03;
            fMesonIntDeltaRangeNarrow[1]= 0.03;
        }

        // Set meson mass ranges (here same for fitting and plotting)
         if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45){
             fMesonMassPlotRange[0]      = 0.65;
             fMesonMassPlotRange[1]      = 0.89;
             fMesonMassRange[0]          = 0.65;
             fMesonMassRange[1]          = 0.89;
         }

         // Set meson fit range
         if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45){
             fMesonFitRange[0]           = 0.65;
             fMesonFitRange[1]           = 0.85;
         }

         // Set remaining parameters for fitting
         if(mode == 40 || mode == 41 || mode == 42 || mode == 44 || mode == 45){
             fMesonWidthExpect           = 0.01;
             fMesonWidthRange[0]         = 0.005;
             fMesonWidthRange[1]         = 0.01;
             fMesonLambdaTail            = 0.007;
             fMesonLambdaTailRange[0]    = 0.004;
             fMesonLambdaTailRange[1]    = 0.01;

             fFullPt[0]                  = 0.4;
             fFullPt[1]                  = 15;

             // Settings for MC
             fMesonLambdaTailMC    = fMesonLambdaTail;
             fMesonWidthExpectMC   = fMesonWidthExpect;
             fMesonWidthRangeMC[0] = fMesonWidthRange[0];
             fMesonWidthRangeMC[1] = fMesonWidthRange[1];

         }

    }
}
