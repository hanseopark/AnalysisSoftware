// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

// Double, Int, etc.
TDatime	now;
fstream fFileErrLog;
fstream	fFileDataLog;

Int_t fCrysFitting = 0;
Int_t fIsMC = 0;
Int_t fMesonId = 0;
Float_t fBackgroundMultNumber;

Double_t fYMaxMeson = 0.9;
Double_t fMesonMassExpect = 0;   // expected meson mass
Double_t fNEvents = 0;

Double_t *fPeakRange;
Double_t *fFitRange;
Double_t *fBGFitRange =	NULL;
Double_t *fMesonIntRange = NULL;

// TSring 
TString fTextCent;
TString fEnergyFlag;
TString fPrefix;
TString fPrefix2;
TString fPeriodFlag;
TString ftextDayth;
TString fdate;
TString	fCollisionSystem;
TString	fTextMeasurement;
TString fCutSelection;
TString fBackgroundMultCutNumber;

TString fcMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
                    "Jul","Aug","Sep","Oct","Nov","Dec"};

TString 	fdirectphoton;
TString 	fThesis ="";

// Output Files
TFile* 	fOutput = NULL;
TFile* 	fOutput1 = NULL;
TFile* 	fOutput2 = NULL;

//TH1D
TH1D*  	fBckNorm=							NULL;
TH1D*  	fSignal=							NULL;
TH1D*  	fRatioSB=							NULL;
TH1D*  	fFittingHistMidPtSignal=				NULL;
TH1D*  	fFittingHistMidPtBackground=			NULL;
TH1D*  	fFittingHistMidPtSignalSub=			NULL;
TH1D*  	fMesonFullPtSignal=					NULL;
TH1D*  	fMesonFullPtBackground=				NULL;
TH1D*  	fMesonFullPtBackNorm=				NULL;
TH1D* 	fHistoMotherZMProj;
TH1D* 	fHistoBckZMProj;
TH1D* 	fNumberOfGoodESDTracks;
TH1D* 	fEventQuality;
TH1D* 	fHistoMappingBackNormInvMass=			NULL;
TH1D* 	fHistoMappingSignalInvMass=			NULL;
TH1D*	fHistoYieldMeson=					NULL;
TH1D*	fHistoYieldMesonBackFit=					NULL;
TH1D*	fHistoYieldMesonPerEvent=			NULL;
TH1D*	fHistoYieldMesonPerEventBackFit=			NULL;
TH1D*	fHistoSignMeson=					NULL;
TH1D*	fHistoSBMeson=						NULL;
TH1D*	fHistoTrueSignMeson=				NULL;
TH1D*	fHistoTrueSBMeson=					NULL;
TH1D*	fHistoYieldMesonNarrow=				NULL;
TH1D*	fHistoYieldMesonPerEventNarrow=		NULL;
TH1D*	fHistoSignMesonNarrow=				NULL;
TH1D*	fHistoSBMesonNarrow=				NULL;
TH1D*	fHistoYieldMesonWide=				NULL;
TH1D*	fHistoYieldMesonPerEventWide=			NULL;
TH1D*	fHistoSignMesonWide=				NULL;
TH1D*	fHistoSBMesonWide=					NULL;
TH1D*	fHistoMassPosition=					NULL;
TH1D*	fHistoMassMeson=					NULL;
TH1D*	fHistoWidthMeson=					NULL;
TH1D*	fHistoFWHMMeson=					NULL;
TH1D*	fHistoTrueMassMeson=				NULL;
TH1D* fHistoTrueMassMesonReweighted=          NULL;
TH1D*	fHistoTrueFWHMMeson=				NULL;
TH1D* fHistoTrueFWHMMesonReweighted=          NULL;
TH1D*	fHistoFWHMMesonAlpha01=				NULL;
TH1D* 	fDeltaPt=							NULL;
TH1D* 	fHistoYieldMesonLeft=				NULL;
TH1D* 	fHistoYieldMesonLeftPerEvent=			NULL;
TH1D* 	fHistoSignMesonLeft=				NULL;
TH1D* 	fHistoSBMesonLeft=					NULL;
TH1D* 	fHistoYieldMesonLeftNarrow=			NULL;
TH1D* 	fHistoYieldMesonLeftPerEventNarrow=	NULL;
TH1D* 	fHistoSignMesonLeftNarrow=			NULL;
TH1D* 	fHistoSBMesonLeftNarrow=				NULL;
TH1D* 	fHistoYieldMesonLeftWide=			NULL;
TH1D* 	fHistoYieldMesonLeftPerEventWide=		NULL;
TH1D* 	fHistoSignMesonLeftWide=				NULL;
TH1D* 	fHistoSBMesonLeftWide;
TH1D* 	fHistoMassMesonLeft=				NULL;
TH1D* 	fHistoWidthMesonLeft=				NULL;
TH1D* 	fHistoFWHMMesonLeft=				NULL;
TH1D* 	fHistoMCMesonPtWithinAcceptance=	NULL;
TH1D*	fHistoMCMesonPt=					NULL;
TH1D* fHistoMCMesonPtWOWeights=              NULL;
TH1D* fHistoMCMesonPtWeights=              NULL;
TH1D* 	fHistoMCMesonWithinAccepPt=			NULL; // Proper bins in Pt
TH1D* 	fHistoMCMesonPt1=					NULL; // Proper bins in Pt
TH1D* 	fHistoYieldTrueMeson=				NULL;
TH1D*  	fHistoYieldTrueMesonWide=			NULL;
TH1D* 	fHistoYieldTrueMesonNarrow=			NULL;
TH1D*    fHistoYieldTrueMesonReweighted=            NULL;
TH1D*    fHistoYieldTrueMesonReweightedWide=        NULL;
TH1D*    fHistoYieldTrueMesonReweightedNarrow=         NULL;

TH1D* 	fHistoYieldTrueGGMeson=				NULL;
TH1D* 	fHistoYieldTrueGGCont=				NULL;
TH1D* 	fHistoYieldTrueDalitzMeson=				NULL;
TH1D* 	fHistoYieldTrueDalitzCont=				NULL;
TH1D* 	fHistoYieldTrueGGMesonWide=				NULL;
TH1D* 	fHistoYieldTrueGGContWide=				NULL;
TH1D* 	fHistoYieldTrueDalitzMesonWide=				NULL;
TH1D* 	fHistoYieldTrueDalitzContWide=				NULL;
TH1D* 	fHistoYieldTrueGGMesonNarrow=				NULL;
TH1D* 	fHistoYieldTrueGGContNarrow=				NULL;
TH1D* 	fHistoYieldTrueDalitzMesonNarrow=				NULL;
TH1D* 	fHistoYieldTrueDalitzContNarrow=				NULL;
TH1D*	fHistoMCMesonAcceptPt=				NULL;
TH1D* 	fHistoMCMesonEffiPt=				NULL;
TH1D* 	fHistoMCMesonEffiFitPt=				NULL;
TH1D*  	fHistoTrueMesonEffiPt=				NULL;
TH1D* 	fHistoMonteMesonEffiPt=				NULL;
TH1D* 	fHistoMonteMesonEffiBackFitPt=				NULL;
TH1D* 	fHistoMonteMesonNarrowEffiPt=			NULL;
TH1D*	fHistoMonteMesonWideEffiPt=			NULL;
TH1D* 	fHistoMonteMesonLeftEffiPt=			NULL;
TH1D*	fHistoMonteMesonLeftNarrowEffiPt=		NULL;
TH1D* 	fHistoMonteMesonLeftWideEffiPt=		NULL;
TH1D*   fHistoMonteMesonEffiMCAll=		NULL;
TH1D * 	fHistoMCTrueMesonEffiPt=				NULL;
TH1D *   fHistoMCTrueMesonEffiPtReweighted=            NULL;
TH1D * 	fHistoMCTrueMesonNarrowEffiPt=		NULL;
TH1D *   fHistoMCTrueMesonNarrowEffiPtReweighted=      NULL;
TH1D *	fHistoMCTrueMesonWideEffiPt=			NULL;
TH1D *   fHistoMCTrueMesonWideEffiPtReweighted=        NULL;
TH1D*  	fHistoTrueMesonEffiFitPt=				NULL;
TH1D* 	fHistoMonteMesonEffiFitPt=				NULL;
TH1D* 	fHistoMonteMesonNarrowEffiFitPt=			NULL;
TH1D*	fHistoMonteMesonWideEffiFitPt=			NULL;
TH1D* 	fHistoMonteMesonLeftEffiFitPt=			NULL;
TH1D*	fHistoMonteMesonLeftNarrowEffiFitPt=		NULL;
TH1D* 	fHistoMonteMesonLeftWideEffiFitPt=		NULL;
TH1D*   fHistoMonteMesonEffiFitMCAll=		NULL;
TH1D * 	fHistoMCTrueMesonEffiFitPt=				NULL;
TH1D * 	fHistoMCTrueMesonNarrowEffiFitPt=		NULL;
TH1D *	fHistoMCTrueMesonWideEffiFitPt=			NULL;

TH1D** 	fHistoMappingTrueMesonInvMassPtBins=	NULL;    
TH1D**   fHistoMappingTrueMesonInvMassPtReweightedBins=   NULL;    
TH1D** 	fHistoMappingTrueGGMesonInvMassPtBins=	NULL;    
TH1D** 	fHistoMappingTrueDalitzInvMassPtBins=	NULL;     
TH1D** 	fHistoMappingPiPlPiMiGammaInvMassPtBin=			NULL;    
TH1D** 	fHistoMappingGGInvMassBackFitPtBin=		NULL;    
TH1D** 	fHistoMappingGGInvMassBackFitWithoutSignalPtBin=		NULL;    
TH1D** 	fHistoMappingBackInvMassPtBin=		NULL;
TH1D** 	fHistoMappingBackNormInvMassPtBin=		NULL;
TH1D** 	fHistoMappingRatioSBInvMassPtBin=		NULL;
TH1D** 	fHistoMappingSignalInvMassPtBin=		NULL;
TH1D** 	fHistoMappingPeakPosInvMassPtBin=		NULL;
// Histograms for normalization on the left of the peak
TH1D** 	fHistoMappingBackNormInvMassLeftPtBin=	NULL;
TH1D** 	fHistoMappingSignalInvMassLeftPtBin=	NULL;


TF1**  fBackgroundFitPol=NULL;
TF1 * 	fFitReco=							NULL;
TF1 * 	fFitGausExp=						NULL;
TF1 * 	fFitLinearBck=						NULL;
TF1 * 	fFitLinearBckOut=					NULL;
TF1 * 	fFitSignalInvMassMidPt=				NULL;

Double_t 	fYields;
Double_t 	fYieldsError;

Double_t 	fFWHMFunc;
Double_t 	fFWHMFuncError;
Double_t 	fYieldsFunc;
Double_t 	fYieldsFuncError;
Double_t 	fIntLinearBck;
Double_t 	fIntLinearBckError;
Double_t 	fIntLinearBckOut;
Double_t 	fIntLinearBckErrorOut;

Float_t 	pictDrawingCoordinatesFWHM[9] = 		{0.6, 0.8, 0.30, 0.04, 0.15,0.7, 0.1, 0.035,0};
Int_t 	fNRebinGlobal = 					2;

// common meson analysis variables (So it is easy to switch between the pi0 and the eta without changing the code)
Double_t 	*fBGFitRangeLeft = 					NULL;
Double_t 	*fMesonPlotRange = 					NULL;

Double_t 	*fMesonIntRangeWide = 				NULL;
Double_t 	*fMesonIntRangeNarrow = 				NULL;
Double_t 	*fMesonCurIntRange = 				NULL;
Double_t 	*fMesonCurIntRangeBackFit = 				NULL;
Double_t 	*fMesonCurIntRangeWide = 			NULL;
Double_t 	*fMesonCurIntRangeNarrow = 			NULL;
Double_t 	*fMesonCurLeftIntRange = 			NULL;
Double_t 	*fMesonCurLeftIntRangeWide = 			NULL;
Double_t 	*fMesonCurLeftIntRangeNarrow = 		NULL;
Double_t 	*fMesonTrueIntRange =				NULL;
Double_t 	*fMesonTrueIntRangeWide = 			NULL;
Double_t 	*fMesonTrueIntRangeNarrow = 			NULL;
Double_t    *fMesonTrueIntReweightedRange =            NULL;
Double_t    *fMesonTrueIntReweightedRangeWide =        NULL;
Double_t    *fMesonTrueIntReweightedRangeNarrow =         NULL;

Double_t 	*fMesonMassRange = 					NULL;
Double_t 	*fMesonFitRange = 					NULL;
Double_t 	*fMesonFitRangeWithoutPeak = 			NULL;

Double_t 	fMesonWidthExpect=					0;
Double_t 	fMesonLambdaTail=					0;
Double_t 	*fMesonWidthRange = 				NULL;
Double_t 	*fMesonLambdaTailRange = 			NULL;
Double_t 	*fMidPt = 						NULL;
Double_t 	*fFullPt = 						NULL;
// end common meson analysis variables

//Background histograms in different M and Z bins

TH2D* fHistoMotherZM;
TH2D* fHistoBckZM;

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
TString 	fNameHistoTrueGG;
TString 	fNameHistoTrueDalitz;
TString 	fNameHistoEffi;
TString 	fNameHistoFrac;
TString 	fNameHistoMotherZM;
TString 	fNameHistoBckZM;
TString 	fNameFitSignalPos;

Double_t* fPiPlPiMiGammaYields = 						NULL;
Double_t* fBckYields = 						NULL;
Double_t* fTotalBckYields= 					NULL;
Double_t* fMesonYields= 						NULL;
Double_t* fMesonYieldsBackFit= 						NULL;
Double_t* fMesonTrueYields= 					NULL;
Double_t* fMesonTrueYieldsReweighted=               NULL;
Double_t* fMesonTrueGGContYields= 					NULL;
Double_t* fMesonTrueDalitzContYields= 					NULL;
Double_t* fMesonTrueGGContYieldsWide= 					NULL;
Double_t* fMesonTrueDalitzContYieldsWide= 					NULL;
Double_t* fMesonTrueGGContYieldsNarrow= 					NULL;
Double_t* fMesonTrueDalitzContYieldsNarrow= 					NULL;
Double_t* fMesonYieldsFunc= 					NULL;
Double_t* fMesonYieldsResidualBckFunc= 			NULL;
Double_t* fMesonYieldsResidualBckFuncBackFit= 			NULL;
Double_t* fMesonYieldsCorResidualBckFunc=		NULL;
Double_t* fMesonYieldsCorResidualBckFuncBackFit=		NULL;
Double_t* fMesonYieldsPerEvent= 				NULL;
Double_t* fMesonYieldsPerEventBackFit= 				NULL;
Double_t* fMesonMass= 						NULL;
Double_t* fMesonMassBackFit= 						NULL;
Double_t* fMesonWidth= 						NULL;
Double_t* fMesonWidthBackFit= 						NULL;
Double_t* fMesonSB= 						NULL;
Double_t* fMesonTrueSB= 						NULL;
Double_t* fMesonTrueSign= 					NULL;
Double_t* fMesonSign= 						NULL;
Double_t* fMesonFWHM= 						NULL;
Double_t* fMesonTrueMass= 					NULL;
Double_t* fMesonTrueMassReweighted=              NULL;
Double_t* fMesonTrueFWHM= 					NULL;
Double_t* fMesonTrueFWHMReweighted=              NULL;
Double_t* fMesonFWHMAlpha01= 					NULL;

// Normalization at the left of the peak
Double_t* fPiPlPiMiGammaYieldsLeft=						NULL;
Double_t* fBckYieldsLeft=					NULL;
Double_t* fTotalBckYieldsLeft=				NULL;
Double_t* fMesonYieldsLeft=					NULL;
Double_t* fMesonYieldsFuncLeft=				NULL;
Double_t* fMesonYieldsResidualBckFuncLeft=		NULL;
Double_t* fMesonYieldsCorResidualBckFuncLeft=	NULL;
Double_t* fMesonYieldsLeftPerEvent=			NULL;
Double_t* fMesonMassLeft=					NULL;
Double_t* fMesonWidthLeft=					NULL;
Double_t* fMesonSBLeft=						NULL;
Double_t* fMesonSignLeft=					NULL;
Double_t* fMesonFWHMLeft=					NULL;

// Narrow Integration Window
Double_t* fPiPlPiMiGammaYieldsNarrow=					NULL;
Double_t* fBckYieldsNarrow=					NULL;
Double_t* fTotalBckYieldsNarrow=				NULL;
Double_t* fMesonYieldsNarrow=					NULL;
Double_t* fMesonTrueYieldsNarrow=				NULL;
Double_t* fMesonTrueYieldsReweightedNarrow=            NULL;
Double_t* fMesonYieldsFuncNarrow=				NULL;
Double_t* fMesonYieldsResidualBckFuncNarrow=		NULL;
Double_t* fMesonYieldsCorResidualBckFuncNarrow=	NULL;
Double_t* fMesonYieldsPerEventNarrow=			NULL;
Double_t* fMesonSBNarrow=					NULL;
Double_t* fMesonSignNarrow=					NULL;

Double_t* fPiPlPiMiGammaYieldsLeftNarrow=				NULL;
Double_t* fBckYieldsLeftNarrow=				NULL;
Double_t* fTotalBckYieldsLeftNarrow=			NULL;
Double_t* fMesonYieldsLeftNarrow=				NULL;
Double_t* fMesonYieldsFuncLeftNarrow=			NULL;
Double_t* fMesonYieldsResidualBckFuncLeftNarrow=	NULL;
Double_t* fMesonYieldsCorResidualBckFuncLeftNarrow=NULL;
Double_t* fMesonYieldsLeftPerEventNarrow=		NULL;
Double_t* fMesonSBLeftNarrow=					NULL;
Double_t* fMesonSignLeftNarrow=				NULL;

// Wide Integration Window
Double_t* fPiPlPiMiGammaYieldsWide=						NULL;
Double_t* fBckYieldsWide=					NULL;
Double_t* fTotalBckYieldsWide=				NULL;
Double_t* fMesonYieldsWide=					NULL;
Double_t* fMesonTrueYieldsWide=				NULL;
Double_t* fMesonTrueYieldsReweightedWide=           NULL;
Double_t* fMesonYieldsFuncWide=				NULL;
Double_t* fMesonYieldsResidualBckFuncWide=		NULL;
Double_t* fMesonYieldsCorResidualBckFuncWide=	NULL;
Double_t* fMesonYieldsPerEventWide=			NULL;
Double_t* fMesonSBWide=						NULL;
Double_t* fMesonSignWide=					NULL;

Double_t* fPiPlPiMiGammaYieldsLeftWide=					NULL;
Double_t* fBckYieldsLeftWide=					NULL;
Double_t* fTotalBckYieldsLeftWide=				NULL;
Double_t* fMesonYieldsLeftWide=				NULL;
Double_t* fMesonYieldsFuncLeftWide=			NULL;
Double_t* fMesonYieldsResidualBckFuncLeftWide=	NULL;
Double_t* fMesonYieldsCorResidualBckFuncLeftWide=	NULL;
Double_t* fMesonYieldsLeftPerEventWide=			NULL;
Double_t* fMesonSBLeftWide=					NULL;
Double_t* fMesonSignLeftWide=					NULL;

Double_t* fPiPlPiMiGammaYieldsError=					NULL;
Double_t* fBckYieldsError=					NULL;
Double_t* fTotalBckYieldsError=				NULL;
Double_t* fMesonYieldsError=					NULL;
Double_t* fMesonYieldsBackFitError=					NULL;
Double_t* fMesonTrueYieldsError=				NULL;
Double_t* fMesonTrueYieldsReweightedError=          NULL;
Double_t* fMesonTrueGGContYieldsError=				NULL;
Double_t* fMesonTrueDalitzContYieldsError=				NULL;
Double_t* fMesonTrueGGContYieldsWideError=				NULL;
Double_t* fMesonTrueDalitzContYieldsWideError=				NULL;
Double_t* fMesonTrueGGContYieldsNarrowError=				NULL;
Double_t* fMesonTrueDalitzContYieldsNarrowError=				NULL;
Double_t* fMesonYieldsFuncError=				NULL;
Double_t* fMesonYieldsResidualBckFuncError=		NULL;
Double_t* fMesonYieldsResidualBckFuncBackFitError=		NULL;
Double_t* fMesonYieldsCorResidualBckFuncError=	NULL;
Double_t* fMesonYieldsCorResidualBckFuncBackFitError=	NULL;
Double_t* fMesonYieldsPerEventError=			NULL;
Double_t* fMesonYieldsPerEventBackFitError=			NULL;
Double_t* fMesonMassError=					NULL;
Double_t* fMesonMassBackFitError=					NULL;
Double_t* fMesonWidthError=					NULL;
Double_t* fMesonWidthBackFitError=					NULL;
Double_t* fMesonTrueMassError=				NULL;
Double_t* fMesonTrueMassReweightedError=            NULL;
Double_t* fMesonTrueFWHMReweightedError=				NULL;
Double_t* fMesonTrueFWHMError=           NULL;

Double_t* fMesonSBError=						NULL;
Double_t* fMesonSignError=					NULL;
Double_t* fMesonTrueSBError= 					NULL;
Double_t* fMesonTrueSignError=				NULL;
Double_t* fMesonFWHMError=					NULL;

Double_t* fPiPlPiMiGammaYieldsLeftError=					NULL;
Double_t* fBckYieldsLeftError=				NULL;
Double_t* fTotalBckYieldsLeftError=			NULL;
Double_t* fMesonYieldsLeftError=				NULL;
Double_t* fMesonYieldsFuncLeftError=			NULL;
Double_t* fMesonYieldsResidualBckFuncLeftError=	NULL;
Double_t* fMesonYieldsCorResidualBckFuncLeftError=NULL;
Double_t* fMesonYieldsLeftPerEventError=		NULL;
Double_t* fMesonMassLeftError=				NULL;
Double_t* fMesonWidthLeftError=				NULL;
Double_t* fMesonSBLeftError=					NULL;
Double_t* fMesonSignLeftError=				NULL;
Double_t* fMesonFWHMLeftError=				NULL;
Double_t* fMesonFWHMAlpha01Error=				NULL;

// Narrow integration Window
Double_t* fPiPlPiMiGammaYieldsNarrowError=				NULL;
Double_t* fBckYieldsNarrowError=				NULL;
Double_t* fTotalBckYieldsNarrowError=			NULL;
Double_t* fMesonYieldsNarrowError=				NULL;
Double_t* fMesonTrueYieldsNarrowError=			NULL;
Double_t* fMesonTrueYieldsReweightedNarrowError=       NULL;
Double_t* fMesonYieldsFuncNarrowError=			NULL;
Double_t* fMesonYieldsResidualBckFuncNarrowError=	NULL;
Double_t* fMesonYieldsCorResidualBckFuncNarrowError=NULL;
Double_t* fMesonYieldsPerEventNarrowError=		NULL;
Double_t* fMesonSBNarrowError=				NULL;
Double_t* fMesonSignNarrowError=				NULL;

Double_t* fPiPlPiMiGammaYieldsLeftNarrowError=			NULL;
Double_t* fBckYieldsLeftNarrowError=			NULL;
Double_t* fTotalBckYieldsLeftNarrowError=		NULL;
Double_t* fMesonYieldsLeftNarrowError=			NULL;
Double_t* fMesonYieldsFuncLeftNarrowError=		NULL;
Double_t* fMesonYieldsResidualBckFuncLeftNarrowError=NULL;
Double_t* fMesonYieldsCorResidualBckFuncLeftNarrowError=NULL;
Double_t* fMesonYieldsLeftPerEventNarrowError=	NULL;
Double_t* fMesonSBLeftNarrowError=				NULL;
Double_t* fMesonSignLeftNarrowError=			NULL;

// Wide integration Window
Double_t* fPiPlPiMiGammaYieldsWideError=					NULL;
Double_t* fBckYieldsWideError=				NULL;
Double_t* fTotalBckYieldsWideError=			NULL;
Double_t* fMesonYieldsWideError=				NULL;
Double_t* fMesonTrueYieldsWideError=			NULL;
Double_t* fMesonTrueYieldsReweightedWideError=         NULL;
Double_t* fMesonYieldsFuncWideError=			NULL;
Double_t* fMesonYieldsResidualBckFuncWideError=	NULL;
Double_t* fMesonYieldsCorResidualBckFuncWideError=NULL;
Double_t* fMesonYieldsPerEventWideError=		NULL;
Double_t* fMesonSBWideError=					NULL;
Double_t* fMesonSignWideError=				NULL;

Double_t* fPiPlPiMiGammaYieldsLeftWideError=				NULL;
Double_t* fBckYieldsLeftWideError=				NULL;
Double_t* fTotalBckYieldsLeftWideError=			NULL;
Double_t* fMesonYieldsLeftWideError=			NULL;
Double_t* fMesonYieldsFuncLeftWideError=		NULL;
Double_t* fMesonYieldsResidualBckFuncLeftWideError=NULL;
Double_t* fMesonYieldsCorResidualBckFuncLeftWideError=NULL;
Double_t* fMesonYieldsLeftPerEventWideError=		NULL;
Double_t* fMesonSBLeftWideError=				NULL;
Double_t* fMesonSignLeftWideError=				NULL;

TF1** 	fFitSignalInvMassPtBin=				NULL;
TF1** 	fFitSignalInvMassBackFitPtBin=			NULL;
TF1** 	fFitTrueSignalInvMassPtBin=			NULL;
TF1**    fFitTrueSignalInvMassPtReweightedBin=         NULL;
TF1** 	fFitSignalPeakPosInvMassPtBin=		NULL;
TF1** 	fFitSignalPeakPosInvMassBackFitPtBin=		NULL;
TF1** 	fFitBckInvMassPtBin=				NULL;
TF1** 	fFitBckInvMassBackFitPtBin=				NULL;
TF1** 	fFitRatioInvMassPtBin=				NULL;
TF1**	fFitWithPol2ForBG = 				NULL;
TF1*     fFitDCAZPhotonTrainBGH = NULL;
TF1*     fFitDCAZPhotonTrainBG = NULL;
TH1F*    fHistDCAZPhotonTrainBG = NULL;
TF1*     fFitDCAZPhotonInterTrainBG = NULL;
TF1*     fFitDCAZPhotonInterTrainBGH = NULL;
TF1*     fFitDCAZPhotonBGH = NULL;
TF1*     fFitDCAZPhotonBG = NULL;
TF1** 	fFitPeakPosPtBin=					NULL;
Double_t* fMesonMassPeakPos=					NULL;
Double_t* fMesonMassPeakPosError=				NULL;

TF1** 	fFitInvMassLeftPtBin=				NULL;
TF1** 	fFitSignalPeakPosInvMassLeftPtBin=		NULL;
TF1** 	fFitBckInvMassLeftPtBin=				NULL;

TH2D* 	fHistoTrueMesonInvMassVSPt=			NULL;
TH2D*    fHistoTrueMesonInvMassVSPtWOWeights=         NULL;
TH2D*    fHistoTrueMesonInvMassVSPtReweighted=         NULL;
TProfile2D*    fProfileTrueMesonInvMassVSPtWeights=         NULL;
TH2D* 	fHistoTrueMesonInvMassVSPtSec=			NULL;
TH2D* 	fHistoTrueGGMesonInvMassVSPt=			NULL;
TH2D* 	fHistoTrueDalitzMesonInvMassVSPt=		NULL;
TH2D*   fGammaGammaInvMassVSPt=				NULL;
TH2D*   fBckInvMassVSPt=				NULL;

TH2F**    fHistoWeightsBGZbinVsMbin = 0x0;
TH2F**    fHistoFillPerEventBGZbinVsMbin = 0x0;


TString ObjectNameTrue;
TString ObjectNameTrueWOWeights;
TString ObjectNameProfileWeights;
TString ObjectNameContaminationGG;
TString ObjectNameContaminationDalitz;
TString ObjectNameMCEtaAcc;
TString ObjectNameMCEta;
TString ObjectNameMCEtaWOWeights;
   
Bool_t fAdvancedMesonQA = kFALSE;
Bool_t fEstimateTrainPileUp = kFALSE;

void ProcessEM(TH1D*,TH1D*,Double_t *);
void ProcessBckFitSubtraction(TH1D *fGammaGamma, Int_t i, Double_t * fPeakRange, Double_t *fFitRange);
void ProcessRatioSignalBackground(TH1D* , TH1D* );
void FillMassHistosArray(TH2D*);
void FillMassMCTrueMesonHistosArray(TH2D*);
void FillMassMCTrueGGMesonHistosArray(TH2D*);
void FillMassMCTrueDalitzMesonHistosArray(TH2D*);
TH1D* CalulateContaminationFraction(TH1D* histoRawYield, TH1D* histoRawYieldSec, TString nameHistoFrac);
void CreatePtHistos();
void FillPtHistos();
void FitSubtractedInvMassInPtBins(TH1D * ,Double_t *, Int_t, Bool_t );  // Fits the Invariant Mass histos with a given function
void FitTrueInvMassInPtBins(TH1D * ,Double_t *, Int_t);  // Fits the Invariant Mass histos with a given function
void FitWithPol2ForBG(TH1D*,Double_t*,Int_t,Bool_t);
void ProduceBckProperWeighting(TList*, TList* );
void ProduceBckWithoutWeighting(TH2D *);
void IntegrateHistoInvMassStream(TH1D * , Double_t *);
void IntegrateHistoInvMass(TH1D * , Double_t *);
void IntegrateFitFunc(TF1 * , TH1D *, Double_t *);
void FillHistosArrayMC(TH1D* , TH1D*, TH1D* ,TString);
void CalculateMesonAcceptance();
void CalculateMesonEfficiency(TH1D*, TH1D*,TString,TString);
void SaveHistos(Int_t, TString, TString);
void SaveCorrectionHistos(TString , TString);
void Initialize(TString setPi0, Int_t);
void CalculateFWHM(TF1 *);
void Delete();
void SetCorrectMCHistogrammNames(TString, TString);

TString centralityString = "";