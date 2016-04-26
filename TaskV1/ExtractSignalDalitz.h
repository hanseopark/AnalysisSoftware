// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion


Int_t           fMode                            = -1;
TString 	fPrefix;
TString 	fPrefix2;
TDatime 	now;

TString 	cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
			"Jul","Aug","Sep","Oct","Nov","Dec"};
TString 	textDayth;
TString 	date;
/*TString 	fSuffix;*/
TString 	fConference;
TString 	fThesis;
TString 	fTextCent;
TString 	fEnergyFlag;
TString 	fPeriodFlag;
TString 	fBackgroundMultCutNumber;
Int_t 	fBackgroundMultNumber;
TString 	fTextMeasurement;
TString 	fCollisionSystem;
TString         fDetectionProcess;
TFile* 	fOutput=					NULL;
TFile* 	fOutput1=					NULL;
TFile* 	fOutput2=					NULL;

TH1D*  	fBckNorm=					NULL;
TH1D*  	fSignal=					NULL;
TH1D*  	fRatioSB=					NULL;
TH1D*  	fFittingHistMidPtSignal=			NULL;
TH1D*  	fFittingHistMidPtBackground=			NULL;
TH1D*  	fFittingHistMidPtSignalSub=			NULL;
TH1D*  	fMesonFullPtSignal=				NULL;
TH1D*  	fMesonFullPtBackground=				NULL;
TH1D*  	fMesonFullPtBackNorm=				NULL;


TF1 * 	fFitReco=					NULL;
TF1 * 	fFitGausExp=					NULL;
TF1 * 	fFitLinearBck=					NULL;
TF1 * 	fFitLinearBckOut=				NULL;
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

const Double_t  factorGGCont = (98.826*10)/(1.174*90);


const Double_t  fPi0GGBRDPG  	= 0.98823; //DPG
const Double_t  fPi0DalitzBRDPG = 0.01174; //DPG 0.01174 
		


Float_t 	pictDrawingCoordinatesFWHM[9] = 		{0.6, 0.8, 0.30, 0.04, 0.15,0.7, 0.15, 0.035,0};

Int_t 	fCrysFitting=						0;
Int_t 	fIsMC=							0;
Double_t 	fYMaxMeson=						0.9;
TString 	fCutSelection;
Int_t 	fNRebinGlobal = 					2;
fstream 	fFileErrLog;
fstream 	fFileDataLog;
// common meson analysis variables (So it is easy to switch between the pi0 and the eta without changing the code)

TArrayD         *fArrayDBinsPt		=			NULL;
TArrayD         *fArrayDMesonMassRange 	=			NULL;
TArrayI         *fArrayIParametersBins	=			NULL;
TArrayD         *fArrayBRPi0Meson       =                       NULL;


Double_t 	*fBGFitRange = 					NULL;
Double_t 	*fBGFitRangeLeft = 				NULL;
Double_t 	*fMesonPlotRange = 				NULL;
Double_t 	*fMesonIntRange = 				NULL;
Double_t 	*fMesonIntRangeWide = 				NULL;
Double_t 	*fMesonIntRangeNarrow = 			NULL;
Double_t 	*fMesonCurIntRange = 				NULL;
Double_t 	*fMesonCurIntRangeWide = 			NULL;
Double_t 	*fMesonCurIntRangeNarrow = 			NULL;
Double_t 	*fMesonCurLeftIntRange = 			NULL;
Double_t 	*fMesonCurLeftIntRangeWide = 			NULL;
Double_t 	*fMesonCurLeftIntRangeNarrow = 			NULL;
Double_t 	*fMesonTrueIntRange =				NULL;
Double_t 	*fMesonTrueIntRangeWide = 			NULL;
Double_t 	*fMesonTrueIntRangeNarrow = 			NULL;

Double_t 	*fMesonMassRange = 				NULL;
Double_t 	*fMesonFitRange = 				NULL;
Double_t 	*fMesonFitRangeWithoutPeak = 			NULL;
Double_t 	fMesonMassExpect=					0;   // expected meson mass
Int_t 	fMesonId=							0;
Double_t 	fMesonWidthExpect=					0;
Double_t 	fMesonLambdaTail=					0;
Double_t 	*fMesonWidthRange = 				NULL;
Double_t 	*fMesonLambdaTailRange = 			NULL;
Double_t 	*fMidPt = 					NULL;
Double_t 	*fFullPt = 					NULL;
// end common meson analysis variables

//Background histograms in different M and Z bins
TH2D* 	 fHistoMotherZM;
TH2D* 	 fHistoBckZM;
TH2F* 	 fPeakPosAlpha01;
TH1D* 	 fHistoMotherZMProj;
TH1D* 	 fHistoBckZMProj;
TH1D* 	 fHistoMotherZMProjFullPt;
TH1D* 	 fHistoBckZMProjFullPt;
TH1D*    fHistoMotherZMProjMidPt;
TH1D*    fHistoBckZMProjMidPt;
Double_t fScalingFactorBckFullPt;
Double_t fScalingFactorBckMidPt;



Double_t 	fScalingFactorBck[7][4];


Double_t 	fCBAlpha;
Double_t 	fCBn;

///////////////////////////

Double_t fNEvents;
TH1F* 	 fNumberOfGoodESDTracksVtx;
TH1F* 	 fEventQuality;
TH2F*    fESDEposEnegInvMassPt;

TString 	fNameHistoGG;
TString 	fNameHistoBack;
TString 	fNameHistoPP;
TString 	fNameHistoBackNormOut;
TString 	fNameHistoBackNormLeftOut;
TString 	fNameHistoTrue;
TString         fNameHistoTrueGG;
TString 	fNameHistoTrueGGBck;
TString 	fNameHistoTrueContBck;
TString         fNameHistoSignalInvMassW0TruePi0; 
TString 	fNameHistoGGInvMassW0TruePi0;
TString 	fNameHistoTrueAllBck;
TString 	fNameHistoTrueSec;
TString 	fNameHistoTrueSecFromK0S;
TString 	fNameHistoEffi;
TString 	fNameHistoFrac;
TString 	fNameHistoMotherZM;
TString 	fNameHistoBckZM;
TString 	fNameFitSignalPos;

Double_t* fGGYields = 						NULL;
Double_t* fBckYields = 						NULL;
Double_t* fTotalBckYields= 					NULL;
Double_t* fMesonYields= 						NULL;
Double_t* fMesonTrueYields= 					NULL;
Double_t* fMesonTrueSecYields= 					NULL;
Double_t* fMesonTrueSecFromK0SYields= 					NULL;
Double_t* fMesonTrueSecYieldsWide= 					NULL;
Double_t* fMesonTrueSecFromK0SYieldsWide= 					NULL;
Double_t* fMesonTrueSecYieldsNarrow= 					NULL;
Double_t* fMesonTrueSecFromK0SYieldsNarrow= 					NULL;
Double_t* fMesonYieldsFunc= 					NULL;
Double_t* fMesonYieldsResidualBckFunc= 			NULL;
Double_t* fMesonYieldsCorResidualBckFunc=		NULL;
Double_t* fMesonYieldsPerEvent= 				NULL;
Double_t* fMesonMass= 						NULL;
Double_t* fMesonWidth= 						NULL;
Double_t* fMesonSB= 						NULL;
Double_t* fMesonTrueSB= 						NULL;
Double_t* fMesonTrueSign= 					NULL;
Double_t* fMesonSign= 						NULL;
Double_t* fMesonFWHM= 						NULL;
Double_t* fMesonTrueMass= 					NULL;
Double_t* fMesonTrueFWHM= 					NULL;
Double_t* fMesonFWHMAlpha01= 					NULL;

Double_t* fMesonTrueGGYields= 					NULL;
Double_t* fMesonTrueGGYieldsError= 					NULL;
Double_t* fMesonTrueGGYieldsNarrow= 					NULL;
Double_t* fMesonTrueGGYieldsNarrowError= 					NULL;
Double_t* fMesonTrueGGYieldsWide= 					NULL;
Double_t* fMesonTrueGGYieldsWideError= 					NULL;


Double_t* 	fMesonIntDeltaRange 				= NULL;
Double_t* 	fMesonIntDeltaRangeWide 			= NULL;
Double_t* 	fMesonIntDeltaRangeNarrow 			= NULL;





// Normalization at the left of the peak
Double_t* fGGYieldsLeft=						NULL;
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
Double_t* fGGYieldsNarrow=					NULL;
Double_t* fBckYieldsNarrow=					NULL;
Double_t* fTotalBckYieldsNarrow=				NULL;
Double_t* fMesonYieldsNarrow=					NULL;
Double_t* fMesonTrueYieldsNarrow=				NULL;
Double_t* fMesonYieldsFuncNarrow=				NULL;
Double_t* fMesonYieldsResidualBckFuncNarrow=		NULL;
Double_t* fMesonYieldsCorResidualBckFuncNarrow=	NULL;
Double_t* fMesonYieldsPerEventNarrow=			NULL;
Double_t* fMesonSBNarrow=					NULL;
Double_t* fMesonSignNarrow=					NULL;

Double_t* fGGYieldsLeftNarrow=				NULL;
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
Double_t* fGGYieldsWide=						NULL;
Double_t* fBckYieldsWide=					NULL;
Double_t* fTotalBckYieldsWide=				NULL;
Double_t* fMesonYieldsWide=					NULL;
Double_t* fMesonTrueYieldsWide=				NULL;
Double_t* fMesonYieldsFuncWide=				NULL;
Double_t* fMesonYieldsResidualBckFuncWide=		NULL;
Double_t* fMesonYieldsCorResidualBckFuncWide=	NULL;
Double_t* fMesonYieldsPerEventWide=			NULL;
Double_t* fMesonSBWide=						NULL;
Double_t* fMesonSignWide=					NULL;

Double_t* fGGYieldsLeftWide=					NULL;
Double_t* fBckYieldsLeftWide=					NULL;
Double_t* fTotalBckYieldsLeftWide=				NULL;
Double_t* fMesonYieldsLeftWide=				NULL;
Double_t* fMesonYieldsFuncLeftWide=			NULL;
Double_t* fMesonYieldsResidualBckFuncLeftWide=	NULL;
Double_t* fMesonYieldsCorResidualBckFuncLeftWide=	NULL;
Double_t* fMesonYieldsLeftPerEventWide=			NULL;
Double_t* fMesonSBLeftWide=					NULL;
Double_t* fMesonSignLeftWide=					NULL;

Double_t* fGGYieldsError=					NULL;
Double_t* fBckYieldsError=					NULL;
Double_t* fTotalBckYieldsError=				NULL;
Double_t* fMesonYieldsError=					NULL;
Double_t* fMesonTrueYieldsError=				NULL;
Double_t* fMesonTrueSecYieldsError=				NULL;
Double_t* fMesonTrueSecFromK0SYieldsError=				NULL;
Double_t* fMesonTrueSecYieldsWideError=				NULL;
Double_t* fMesonTrueSecFromK0SYieldsWideError=				NULL;
Double_t* fMesonTrueSecYieldsNarrowError=				NULL;
Double_t* fMesonTrueSecFromK0SYieldsNarrowError=				NULL;
Double_t* fMesonYieldsFuncError=				NULL;
Double_t* fMesonYieldsResidualBckFuncError=		NULL;
Double_t* fMesonYieldsCorResidualBckFuncError=	NULL;
Double_t* fMesonYieldsPerEventError=			NULL;
Double_t* fMesonMassError=					NULL;
Double_t* fMesonWidthError=					NULL;
Double_t* fMesonTrueMassError=				NULL;
Double_t* fMesonTrueFWHMError=				NULL;

Double_t* fMesonSBError=						NULL;
Double_t* fMesonSignError=					NULL;
Double_t* fMesonTrueSBError= 					NULL;
Double_t* fMesonTrueSignError=				NULL;
Double_t* fMesonFWHMError=					NULL;

Double_t* fGGYieldsLeftError=					NULL;
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
Double_t* fGGYieldsNarrowError=				NULL;
Double_t* fBckYieldsNarrowError=				NULL;
Double_t* fTotalBckYieldsNarrowError=			NULL;
Double_t* fMesonYieldsNarrowError=				NULL;
Double_t* fMesonTrueYieldsNarrowError=			NULL;
Double_t* fMesonYieldsFuncNarrowError=			NULL;
Double_t* fMesonYieldsResidualBckFuncNarrowError=	NULL;
Double_t* fMesonYieldsCorResidualBckFuncNarrowError=NULL;
Double_t* fMesonYieldsPerEventNarrowError=		NULL;
Double_t* fMesonSBNarrowError=				NULL;
Double_t* fMesonSignNarrowError=				NULL;

Double_t* fGGYieldsLeftNarrowError=			NULL;
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
Double_t* fGGYieldsWideError=					NULL;
Double_t* fBckYieldsWideError=				NULL;
Double_t* fTotalBckYieldsWideError=			NULL;
Double_t* fMesonYieldsWideError=				NULL;
Double_t* fMesonTrueYieldsWideError=			NULL;
Double_t* fMesonYieldsFuncWideError=			NULL;
Double_t* fMesonYieldsResidualBckFuncWideError=	NULL;
Double_t* fMesonYieldsCorResidualBckFuncWideError=NULL;
Double_t* fMesonYieldsPerEventWideError=		NULL;
Double_t* fMesonSBWideError=					NULL;
Double_t* fMesonSignWideError=				NULL;

Double_t* fGGYieldsLeftWideError=				NULL;
Double_t* fBckYieldsLeftWideError=				NULL;
Double_t* fTotalBckYieldsLeftWideError=			NULL;
Double_t* fMesonYieldsLeftWideError=			NULL;
Double_t* fMesonYieldsFuncLeftWideError=		NULL;
Double_t* fMesonYieldsResidualBckFuncLeftWideError=NULL;
Double_t* fMesonYieldsCorResidualBckFuncLeftWideError=NULL;
Double_t* fMesonYieldsLeftPerEventWideError=		NULL;
Double_t* fMesonSBLeftWideError=				NULL;
Double_t* fMesonSignLeftWideError=				NULL;

TH1D* 	fHistoMappingBackNormInvMass=			NULL;
TH1D* 	fHistoMappingSignalInvMass=			NULL;

TH1D** 	fHistoMappingTrueMesonInvMassPtBins=	NULL;    
TH1D** 	fHistoMappingTrueGGMesonInvMassPtBins=	NULL;
//TH1D**  fHistoMappingTrueBckGGInvMassPtBins=    NULL;
//TH1D**  fHistoMappingTrueBckContInvMassPtBins = NULL;
TH1D** 	fHistoMappingTrueGGBckInvMassPtBins       = NULL;   
TH1D** 	fHistoMappingTrueContBckInvMassPtBins     = NULL;    
TH1D**  fHistoMappingSignalInvMassW0TruePi0PtBins = NULL;
TH1D**  fHistoMappingGGInvMassW0TruePi0PtBins = NULL;
TH1D** 	fHistoMappingTrueAllBckInvMassPtBins=	NULL;    
TH1D** 	fHistoMappingTrueSecMesonInvMassPtBins=	NULL;    
TH1D** 	fHistoMappingTrueSecFromK0SMesonInvMassPtBins=	NULL;     
TH1D** 	fHistoMappingGGInvMassPtBin=			NULL;  
TH1D**  fHistoMappingEEInvMassPtBin =                   NULL;
TH1D** 	fHistoMappingBackInvMassPtBin=		NULL;
TH1D** 	fHistoMappingBackNormInvMassPtBin=		NULL;
TH1D** 	fHistoMappingRatioSBInvMassPtBin=		NULL;
TH1D** 	fHistoMappingSignalInvMassPtBin=		NULL;
TH1D** 	fHistoMappingPeakPosInvMassPtBin=		NULL;
TF1** 	fFitSignalInvMassPtBin=				NULL;
TF1** 	fFitTrueSignalInvMassPtBin=			NULL;
TF1** 	fFitSignalPeakPosInvMassPtBin=		NULL;
TF1** 	fFitBckInvMassPtBin=				NULL;
TF1** 	fFitRatioInvMassPtBin=				NULL;
TF1**	fFitWithPol2ForBG = 				NULL;
TH1D*	fHistoYieldMeson=					NULL;
TH1D*	fHistoYieldMesonPerEvent=			NULL;
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
TH1D*	fHistoTrueFWHMMeson=				NULL;
TH1D*	fHistoFWHMMesonAlpha01=				NULL;

TH1D* 	fDeltaPt=							NULL;

TF1** 	fFitPeakPosPtBin=					NULL;
Double_t* fMesonMassPeakPos=					NULL;
Double_t* fMesonMassPeakPosError=				NULL;

// Histograms for normalization on the left of the peak
TH1D** 	fHistoMappingBackNormInvMassLeftPtBin=	NULL;
TH1D** 	fHistoMappingSignalInvMassLeftPtBin=	NULL;

TF1** 	fFitInvMassLeftPtBin=				NULL;
TF1** 	fFitSignalPeakPosInvMassLeftPtBin=		NULL;
TF1** 	fFitBckInvMassLeftPtBin=				NULL;

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

TH2F* 	fHistoMCMesonPtEtaWithinAcceptance=		NULL;
TH1D*	fHistoMCMesonPt=				NULL;
TH1D*   fHistoMCMesonDalitzPt =				NULL;
TH1D*   fHistoMCMesonGGPt     =          		NULL;

TH1D*  	fHistoMCMesonPtWithinAcceptance=		NULL;
TH1D*   fHistoMCMesonPtWithinAcceptance0=               NULL;
TH1D* 	fHistoMCMesonWithinAccepPt=			NULL; // Proper bins in Pt
TH1D* 	fHistoMCMesonPt1=					NULL; // Proper bins in Pt

TH1D* 	fHistoYieldTrueMeson=				NULL;
TH1D*  	fHistoYieldTrueMesonWide=			NULL;
TH1D* 	fHistoYieldTrueMesonNarrow=			NULL;

//TH1D*   fHistoRAWYieldTrueGGFracMeson=                           NULL;
//TH1D*   fHistoRAWYieldTrueGGFracMesonWide=                       NULL;
//TH1D*   fHistoRAWYieldTrueGGFracMesonNarrow=                     NULL;


TH1D* 	fHistoYieldTrueSecMeson=				NULL;
TH1D* 	fHistoYieldTrueSecFracMeson=				NULL;
TH1D* 	fHistoYieldTrueSecFromK0SMeson=				NULL;
TH1D* 	fHistoYieldTrueSecFracFromK0SMeson=				NULL;
TH1D* 	fHistoYieldTrueSecMesonWide=				NULL;
TH1D* 	fHistoYieldTrueSecFracMesonWide=				NULL;
TH1D* 	fHistoYieldTrueSecFromK0SMesonWide=				NULL;
TH1D* 	fHistoYieldTrueSecFracFromK0SMesonWide=				NULL;
TH1D* 	fHistoYieldTrueSecMesonNarrow=				NULL;
TH1D* 	fHistoYieldTrueSecFracMesonNarrow=				NULL;
TH1D* 	fHistoYieldTrueSecFromK0SMesonNarrow=				NULL;
TH1D* 	fHistoYieldTrueSecFracFromK0SMesonNarrow=				NULL;

TH1D*	fHistoMCMesonAcceptPt=				NULL;
TH2F* 	fHistoTrueMesonInvMassVSPt=			NULL;
TH2F* 	fHistoTrueGGMesonInvMassVSPt=			NULL;
TH2F*   fHistoTrueMesonBckGGInvMassVSPt =               NULL;
TH2F*   fHistoTrueMesonBckContInvMassVSPt =             NULL;
TH2F* 	fHistoTrueMesonInvMassVSPtSec=			NULL;
TH2F* 	fHistoTrueGGBckInvMassVSPt=			NULL;
TH2F* 	fHistoTrueContBckInvMassVSPt=			NULL;
TH2F* 	fHistoTrueAllBckInvMassVSPt=			NULL;
TH2F* 	fHistoTrueSecMesonInvMassVSPt=			NULL;
TH2F* 	fHistoTrueSecFromK0SMesonInvMassVSPt=			NULL;

TH2F**    fHistoWeightsBGZbinVsMbin 		= 0x0;
TH2F**    fHistoFillPerEventBGZbinVsMbin 	= 0x0;



TH1D* 	fHistoMCMesonEffiPt=				NULL;
TH1D*  	fHistoTrueMesonEffiPt=				NULL;
TH1D*   fHistoMCRAWYieldTrueGGMeson =                   NULL;

TH1D* 	fHistoMonteMesonEffiPt=				NULL;
TH1D* 	fHistoMonteMesonNarrowEffiPt=			NULL;
TH1D*	fHistoMonteMesonWideEffiPt=			NULL;
TH1D* 	fHistoMonteMesonLeftEffiPt=			NULL;
TH1D*	fHistoMonteMesonLeftNarrowEffiPt=		NULL;
TH1D* 	fHistoMonteMesonLeftWideEffiPt=		NULL;

TH1D * 	fHistoMCTrueMesonEffiPt=				NULL;
TH1D * 	fHistoMCTrueMesonNarrowEffiPt=		NULL;
TH1D *	fHistoMCTrueMesonWideEffiPt=			NULL;

TH1D * 	fHistoYieldTrueGGMeson=				NULL;
TH1D * 	fHistoYieldTrueGGMesonWide=				NULL;
TH1D * 	fHistoYieldTrueGGMesonNarrow=				NULL;
TH1D * 	fHistoYieldTrueGGFracMeson=				NULL;
TH1D * 	fHistoYieldTrueGGFracMesonWide=				NULL;
TH1D * 	fHistoYieldTrueGGFracMesonNarrow=			NULL;
TH1D * 	fHistoYieldTrueGGFracMesonForData=			NULL;
TH1D * 	fHistoYieldTrueGGFracMesonWideForData=			NULL;
TH1D * 	fHistoYieldTrueGGFracMesonNarrowForData=		NULL;


void ProcessEM(TH1D*,TH1D*,Double_t *);
void FillMassHistosArray(TH2F*);
void FillEposEnegHistosArray(TH2F*);
void FillMassMCTrueMesonHistosArray(TH2F*);
void FillMassMCTrueGGMesonHistosArray(TH2F*);
void FillMassMCTrueMesonGGBckHistosArray(TH2F*);
void FillMassMCTrueMesonBckContHistosArray(TH2F*);
void FillSignalInvMassW0TruePi0HistosArray();
void FillGGInvMassW0TruePi0HistosArray();
void CreatePtHistos();
void FillPtHistos();
void FitSubtractedInvMassInPtBins(TH1D * ,Double_t *, Int_t, Bool_t );  // Fits the Invariant Mass histos with a given function
void FitTrueInvMassInPtBins(TH1D * ,Double_t *, Int_t);  // Fits the Invariant Mass histos with a given function
void FitCBSubtractedInvMassInPtBins(TH1D * ,Double_t *, Int_t, Bool_t , TString );  // Fits the Invariant Mass histos with a given function;
void ProduceBckProperWeighting(TList*,TList*, TList* );
void ProduceBckWithoutWeighting(TH2F *);
void IntegrateHistoInvMassStream(TH1D * , Double_t *);
void IntegrateHistoInvMass(TH1D * , Double_t *);
void IntegrateFitFunc(TF1 * , TH1D *, Double_t *);
void FillHistosArrayMC(TH1D* , TH1D * , TH1D * );
void CalculateMesonAcceptance();
void CalculateMesonEfficiency(TH1D*, TString);
void CalculateMesonEfficiencyWithoutGGCont(TH1D*, TH1D*, TString, TString);
TH1D* CalculateSecondaryFractions(TH1D* , TH1D* , TString /*, Bool_t */);
void SaveHistos(Int_t, TString, TString);
void SaveCorrectionHistos(TString , TString);
void Initialize(TString setPi0, Int_t);
void CalculateFWHM(TF1 *);
Double_t CrystalBallBck(Double_t *,Double_t *);
Double_t CrystalBall(Double_t *,Double_t *);
void Delete();


