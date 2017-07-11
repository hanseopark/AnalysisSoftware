// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

// Double, Int, etc.
TDatime	now;
fstream fFileErrLog;
fstream	fFileDataLog;

TString fTypeCutSelection                                   ="";
TString fEventCutSelection                                  ="";
TString fGammaCutSelection                                  ="";
TString fClusterCutSelection                                ="";
TString fPionCutSelection                                   ="";
TString fNeutralPionCutSelection                            ="";
TString fMesonCutSelection                                  ="";

Int_t fCrysFitting = 										0;
Int_t fIsMC = 												0;
Int_t fMesonId = 											0;
Float_t fBackgroundMultNumber;

Double_t fYMaxMeson = 										0.9;
Double_t fMesonMassExpect = 								0;   // expected meson mass
Double_t fNEvents = 										0;

Double_t *fPeakRange;
Double_t *fIntFixedRange;
Double_t *fFitRange;
Double_t *fBGFitRange =										NULL;
Double_t *fMesonIntRange = 									NULL;

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
TString fDetectionProcess;
TString fCutSelection;
TString fBackgroundMultCutNumber;
Int_t 	fMode;

TString fcMonth[12]=										{"Jan","Feb","Mar","Apr","May","Jun",
															"Jul","Aug","Sep","Oct","Nov","Dec"};

TString 	fdirectphoton;
TString 	fThesis ="";

// Output Files
TFile* 	fOutput =							 				NULL;
TFile* 	fOutput1 = 											NULL;
TFile* 	fOutput2 = 											NULL;

//TH1D
TH1D*  	fBckNorm=											NULL;
TH1D*  	fSignal=											NULL;
TH1D*  	fRatioSB=											NULL;
TH1D*  	fFittingHistMidPtSignal=							NULL;
TH1D*  	fFittingHistMidPtBackground=						NULL;
TH1D*  	fFittingHistMidPtSignalSub=							NULL;
TH1D*  	fMesonFullPtSignal=									NULL;
TH1D*  	fMesonFullPtBackground=								NULL;
TH1D*  	fOmegaFullPtSignal=									NULL;
TH2D*  	fOmegaFullPtSignalR2=								NULL;
TH1D*  	fOmegaFullPtSignalR1=								NULL;
TH1D*  	fOmegaFullPtBackNorm=								NULL;
TH1D*  	fOmegaFullPtBack2=                                  NULL;
TH1D*  	fOmegaFullPtBack3=                                  NULL;
TH1D*  	fOmegaFullPtBack4=                                  NULL;
TH1D*   fRelation=                                          NULL;
TH1D*  	fMesonFullPtBackNorm=								NULL;
TH1D* 	fHistoMotherZMProj=									NULL;
TH1D* 	fHistoBckZMProj=									NULL;
TH1D* 	fNumberOfGoodESDTracks=								NULL;
TH1D* 	fEventQuality=										NULL;
TH1D* 	fHistoMappingBackNormInvMass=						NULL;
TH1D* 	fHistoMappingSignalInvMass=							NULL;
TH1D*   fHistoYieldMeson[6]                                         = {NULL, NULL, NULL, NULL, NULL, NULL};
TH1D*	fHistoYieldMesonBackFit=							NULL;
TH1D*   fHistoYieldMesonPerEvent[6]                                 = {NULL, NULL, NULL, NULL, NULL, NULL};
TH1D*	fHistoYieldMesonPerEventBackFit=					NULL;
TH1D*   fHistoSBdefaultMeson[3]=                            { NULL, NULL, NULL};
TH1D*   fHistoSigndefaultMeson[3]=                          { NULL, NULL, NULL};
TH1D*	fHistoTrueSignMeson=								NULL;
TH1D*	fHistoTrueSBMeson=									NULL;
TH1D*	fHistoYieldMesonNarrow=								NULL;
TH1D*	fHistoYieldMesonPerEventNarrow=						NULL;
TH1D*	fHistoSignMesonNarrow=								NULL;
TH1D*	fHistoSBMesonNarrow=								NULL;
TH1D*	fHistoYieldMesonWide=								NULL;
TH1D*	fHistoYieldMesonPerEventWide=						NULL;
TH1D*	fHistoSignMesonWide=								NULL;
TH1D*	fHistoSBMesonWide=									NULL;
TH1D*	fHistoMassPosition=									NULL;
TH1D*	fHistoMassMeson=									NULL;
TH1D*	fHistoWidthMeson=									NULL;
TH1D*	fHistoFWHMMeson=									NULL;
TH1D*	fHistoTrueMassMeson=								NULL;
TH1D*	fHistoTrueMassMesonCaloPhoton=						NULL;
TH1D*	fHistoTrueMassMesonCaloElectron=					NULL;
TH1D*	fHistoTrueMassMesonCaloConvPhoton=					NULL;
TH1D*	fHistoTrueMassMesonCaloMergedCluster=				NULL;
TH1D*	fHistoTrueMassMesonCaloMergedPartConvCluster=		NULL;
TH1D*	fHistoTrueMassMesonCaloEMNonLeading=				NULL;
TH1D* 	fHistoTrueMassMesonReweighted=						NULL;
TH1D*	fHistoTrueFWHMMeson=								NULL;
TH1D*	fHistoTrueFWHMMesonCaloPhoton=						NULL;
TH1D*	fHistoTrueFWHMMesonCaloElectron=					NULL;
TH1D*	fHistoTrueFWHMMesonCaloConvPhoton=					NULL;
TH1D*	fHistoTrueFWHMMesonCaloMergedCluster=				NULL;
TH1D*	fHistoTrueFWHMMesonCaloMergedPartConvCluster=		NULL;
TH1D*	fHistoTrueFWHMMesonCaloEMNonLeading=				NULL;
TH1D*	fHistoTrueFWHMMesonReweighted=						NULL;
TH1D*	fHistoFWHMMesonAlpha01=								NULL;
TH1D* 	fDeltaPt=											NULL;
TH1D* 	fHistoYieldMesonLeft=								NULL;
TH1D* 	fHistoYieldMesonLeftPerEvent=						NULL;
TH1D* 	fHistoSignMesonLeft=								NULL;
TH1D* 	fHistoSBMesonLeft=									NULL;
TH1D* 	fHistoYieldMesonLeftNarrow=							NULL;
TH1D* 	fHistoYieldMesonLeftPerEventNarrow=					NULL;
TH1D* 	fHistoSignMesonLeftNarrow=							NULL;
TH1D* 	fHistoSBMesonLeftNarrow=							NULL;
TH1D* 	fHistoYieldMesonLeftWide=							NULL;
TH1D* 	fHistoYieldMesonLeftPerEventWide=					NULL;
TH1D* 	fHistoSignMesonLeftWide=							NULL;
TH1D* 	fHistoSBMesonLeftWide= 								NULL;
TH1D* 	fHistoMassMesonLeft=								NULL;
TH1D* 	fHistoWidthMesonLeft=								NULL;
TH1D* 	fHistoFWHMMesonLeft=								NULL;
TH1D* 	fHistoMCMesonPtWithinAcceptance=					NULL;
TH1D*	fHistoMCMesonPt=									NULL;
TH1D* 	fHistoMCMesonPtWOWeights=							NULL;
TH1D* 	fHistoMCMesonPtWeights=								NULL;
TH1D* 	fHistoMCMesonWithinAccepPt=							NULL; // Proper bins in Pt
TH1D* 	fHistoMCMesonPt1=									NULL; // Proper bins in Pt
TH1D*   fHistoYieldTrueMeson[6]                           = {NULL, NULL, NULL, NULL, NULL, NULL};
TH1D*  	fHistoYieldTrueMesonWide=							NULL;
TH1D* 	fHistoYieldTrueMesonNarrow=							NULL;
TH1D*   fHistoYieldTrueMesonReweighted[3]                           = {NULL, NULL, NULL};
TH1D*	fHistoYieldTrueMesonReweightedWide=					NULL;
TH1D*	fHistoYieldTrueMesonReweightedNarrow=				NULL;

TH1D*   fHistoYieldTrueSecMeson[3][4]                   = {{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL}};
TH1D*   fHistoYieldTrueSecFracMeson[3][4]               = {{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL}};
TH1D* 	fHistoYieldTrueSecFromK0SMeson=						NULL;
TH1D* 	fHistoYieldTrueSecFracFromK0SMeson=					NULL;
TH1D* 	fHistoYieldTrueSecMesonWide=						NULL;
TH1D* 	fHistoYieldTrueSecFracMesonWide=					NULL;
TH1D* 	fHistoYieldTrueSecFromK0SMesonWide=					NULL;
TH1D* 	fHistoYieldTrueSecFracFromK0SMesonWide=				NULL;
TH1D* 	fHistoYieldTrueSecMesonNarrow=						NULL;
TH1D* 	fHistoYieldTrueSecFracMesonNarrow=					NULL;
TH1D* 	fHistoYieldTrueSecFromK0SMesonNarrow=				NULL;
TH1D* 	fHistoYieldTrueSecFracFromK0SMesonNarrow=			NULL;
TH1D*	fHistoMCMesonAcceptPt=								NULL;
TH1D* 	fHistoMCMesonEffiPt=								NULL;
TH1D* 	fHistoMCMesonEffiFitPt=								NULL;
TH1D*  	fHistoTrueMesonEffiPt=								NULL;
TH1D*   fHistoMonteMesonEffiPt[6]                         = {NULL, NULL, NULL, NULL, NULL, NULL};
TH1D* 	fHistoMonteMesonEffiBackFitPt=						NULL;
TH1D* 	fHistoMonteMesonNarrowEffiPt=						NULL;
TH1D*	fHistoMonteMesonWideEffiPt=							NULL;
TH1D* 	fHistoMonteMesonLeftEffiPt=							NULL;
TH1D*	fHistoMonteMesonLeftNarrowEffiPt=					NULL;
TH1D* 	fHistoMonteMesonLeftWideEffiPt=						NULL;
TH1D*   fHistoMonteMesonEffiMCAll=							NULL;
TH1D*   fHistoMCTrueMesonEffiPt[3]                        = {NULL, NULL, NULL};
TH1D*	fHistoMCTrueMesonEffiPtReweighted[3]			  = {NULL, NULL, NULL};
TH1D* 	fHistoMCTrueMesonNarrowEffiPt=						NULL;
TH1D*	fHistoMCTrueMesonNarrowEffiPtReweighted=			NULL;
TH1D*	fHistoMCTrueMesonWideEffiPt=						NULL;
TH1D*	fHistoMCTrueMesonWideEffiPtReweighted=				NULL;
TH1D*  	fHistoTrueMesonEffiFitPt=							NULL;
TH1D* 	fHistoMonteMesonEffiFitPt=							NULL;
TH1D* 	fHistoMonteMesonNarrowEffiFitPt=					NULL;
TH1D*	fHistoMonteMesonWideEffiFitPt=						NULL;
TH1D* 	fHistoMonteMesonLeftEffiFitPt=						NULL;
TH1D*	fHistoMonteMesonLeftNarrowEffiFitPt=				NULL;
TH1D* 	fHistoMonteMesonLeftWideEffiFitPt=					NULL;
TH1D*   fHistoMonteMesonEffiFitMCAll=						NULL;
TH1D* 	fHistoMCTrueMesonEffiFitPt[3]				    =  {NULL, NULL, NULL};
TH1D* 	fHistoMCTrueMesonNarrowEffiFitPt=					NULL;
TH1D*	fHistoMCTrueMesonWideEffiFitPt=						NULL;

TH1D** 	fHistoMappingTrueMesonInvMassPtBins=				NULL;    
TH1D**	fHistoMappingTrueMesonInvMassPtReweightedBins=		NULL;    
TH1D** 	fHistoMappingTrueGGBckInvMassPtBins=				NULL;   
TH1D** 	fHistoMappingTrueContBckInvMassPtBins=				NULL;    
TH1D** 	fHistoMappingTrueAllBckInvMassPtBins=				NULL;    
TH1D** 	fHistoMappingTrueSecMesonInvMassPtBins=				NULL;    
TH1D** 	fHistoMappingTrueSecFromK0SMesonInvMassPtBins=		NULL;     
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

Double_t 	fYields;
Double_t 	fFWHMFunc;
Double_t 	fFWHMFuncError;
Double_t 	fYieldsFunc;
Double_t 	fYieldsFuncError;
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
TString 	fNameHistoTrueSec;
TString 	fNameHistoTrueSecFromK0S;
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
Double_t* fMesonTrueSecYields[3][4]                                   = {{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL}};
Double_t* fMesonTrueSecFromK0SYields= 						NULL;
Double_t* fMesonTrueSecYieldsWide= 							NULL;
Double_t* fMesonTrueSecFromK0SYieldsWide= 					NULL;
Double_t* fMesonTrueSecYieldsNarrow= 						NULL;
Double_t* fMesonTrueSecFromK0SYieldsNarrow= 				NULL;
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
Double_t*   fMesonSBdefault[3]                                          = { NULL, NULL, NULL};
Double_t*   fMesonSigndefault[3]                                        = { NULL, NULL, NULL};
Double_t*   fMesonSBdefaultError[3]                                     = { NULL, NULL, NULL};
Double_t*   fMesonSigndefaultError[3]                                   = { NULL, NULL, NULL};
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
Double_t* fMesonSBLeft=										NULL;
Double_t* fMesonSignLeft=									NULL;
Double_t* fMesonFWHMLeft=									NULL;
Double_t* fMesonMassLeftError=                              NULL;
Double_t* fMesonFWHMLeftError=                              NULL;
Double_t* fMesonWidthLeftError=                              NULL;

// Narrow Integration Window
Double_t* fMesonSBNarrow=									NULL;
Double_t* fMesonSignNarrow=									NULL;

Double_t* fMesonSBLeftNarrow=								NULL;
Double_t* fMesonSignLeftNarrow=								NULL;

Double_t fScaleFac  =                                         1.;

// Wide Integration Window
Double_t* fMesonSBWide=										NULL;
Double_t* fMesonSignWide=									NULL;

Double_t* fMesonSBLeftWide=									NULL;
Double_t* fMesonSignLeftWide=								NULL;

Double_t*   fGGYieldsError[6] =                            {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t*   fBckYieldsError[6] =                           {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t*   fYieldsError[6] =                              {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsError[6]=							   {NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonYieldsBackFitError=						    NULL;
Double_t* fMesonTrueYieldsError[3]=                        {NULL, NULL, NULL};
Double_t* fMesonTrueYieldsReweightedError[3]=              {NULL, NULL, NULL};
Double_t* fMesonTrueSecYieldsError[3][4] =                  {{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL}};
Double_t* fMesonTrueSecFromK0SYieldsError=					NULL;
Double_t* fMesonTrueSecYieldsWideError=						NULL;
Double_t* fMesonTrueSecFromK0SYieldsWideError=				NULL;
Double_t* fMesonTrueSecYieldsNarrowError=					NULL;
Double_t* fMesonTrueSecFromK0SYieldsNarrowError=			NULL;
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
Double_t* fMesonTrueMassErrorCaloPhoton= 					NULL;
Double_t* fMesonTrueMassErrorCaloElectron= 					NULL;
Double_t* fMesonTrueMassErrorCaloConvPhoton=				NULL;
Double_t* fMesonTrueMassErrorCaloMergedCluster= 			NULL;
Double_t* fMesonTrueMassErrorCaloMergedClusterPartConv= 	NULL;
Double_t* fMesonTrueMassErrorCaloEMNonLeading= 				NULL;
Double_t* fMesonTrueFWHMErrorCaloPhoton= 					NULL;
Double_t* fMesonTrueFWHMErrorCaloElectron= 					NULL;
Double_t* fMesonTrueFWHMErrorCaloConvPhoton=				NULL;
Double_t* fMesonTrueFWHMErrorCaloMergedCluster= 			NULL;
Double_t* fMesonTrueFWHMErrorCaloMergedClusterPartConv= 	NULL;
Double_t* fMesonTrueFWHMErrorCaloEMNonLeading= 				NULL;


Double_t* fMesonSBError=									NULL;
Double_t* fMesonSignError=									NULL;
Double_t* fMesonTrueSBError= 								NULL;
Double_t* fMesonTrueSignError=								NULL;
Double_t* fMesonFWHMError=									NULL;

Double_t* fTotalBckYieldsError[6]=							{NULL, NULL, NULL, NULL, NULL, NULL};
Double_t* fMesonFWHMAlpha01Error=							NULL;

// Narrow integration Window
Double_t* fMesonSBNarrowError=								NULL;
Double_t* fMesonSignNarrowError=							NULL;

Double_t* fMesonSBLeftNarrowError=							NULL;
Double_t* fMesonSignLeftNarrowError=						NULL;

// Wide integration Window
Double_t* fMesonSBWideError=								NULL;
Double_t* fMesonSignWideError=								NULL;

Double_t* fMesonSBLeftWideError=							NULL;
Double_t* fMesonSignLeftWideError=							NULL;

TF1** 	fFitSignalInvMassPtBin=								NULL;
TF1** 	fFitSignalInvMassBackFitPtBin=						NULL;
TF1** 	fFitTrueSignalInvMassPtBin=							NULL;
TF1**	fFitTrueSignalInvMassPtReweightedBin=				NULL;
TF1** 	fFitTrueSignalCaloPhotonInvMassPtBin=				NULL;
TF1** 	fFitTrueSignalCaloConvPhotonInvMassPtBin=			NULL;
TF1** 	fFitTrueSignalCaloElectronInvMassPtBin=				NULL;
TF1** 	fFitTrueSignalCaloEMNonLeadingInvMassPtBin=			NULL;
TF1** 	fFitTrueSignalCaloMergedClusterInvMassPtBin=		NULL;
TF1** 	fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin=NULL;
TF1** 	fFitSignalPeakPosInvMassPtBin=						NULL;
TF1** 	fFitSignalPeakPosInvMassBackFitPtBin=				NULL;
TF1** 	fFitBckInvMassPtBin=								NULL;
TF1** 	fFitBckInvMassBackFitPtBin=							NULL;
TF1** 	fFitRatioInvMassPtBin=								NULL;
TF1**	fFitWithPol2ForBG = 								NULL;
TF1*	fFitDCAZPhotonTrainBGH = 							NULL;
TF1*	fFitDCAZPhotonTrainBG = 							NULL;
TH1F*	fHistDCAZPhotonTrainBG = 							NULL;
TF1*	fFitDCAZPhotonInterTrainBG = 						NULL;
TF1*	fFitDCAZPhotonInterTrainBGH = 						NULL;
TF1*	fFitDCAZPhotonBGH = 								NULL;
TF1*	fFitDCAZPhotonBG = 									NULL;
TF1** 	fFitPeakPosPtBin=									NULL;
Double_t* fMesonMassPeakPos=								NULL;
Double_t* fMesonMassPeakPosError=							NULL;

TF1** 	fFitInvMassLeftPtBin=								NULL;
TF1** 	fFitSignalPeakPosInvMassLeftPtBin=					NULL;
TF1** 	fFitBckInvMassLeftPtBin=							NULL;

TH2D* 	fHistoTrueMesonInvMassVSPt=							NULL;
TH2D*	fHistoTrueMesonInvMassVSPtWOWeights=				NULL;
TH2D*	fHistoTrueMesonInvMassVSPtReweighted=				NULL;
TProfile2D*	fProfileTrueMesonInvMassVSPtWeights=		 	NULL;
TH2D* 	fHistoTrueMesonInvMassVSPtSec=						NULL;
TH2D* 	fHistoTrueGGBckInvMassVSPt=							NULL;
TH2D* 	fHistoTrueContBckInvMassVSPt=						NULL;
TH2D* 	fHistoTrueAllBckInvMassVSPt=						NULL;
TH2D* 	fHistoTrueSecMesonInvMassVSPt=						NULL;
TH2D* 	fHistoTrueSecFromK0SMesonInvMassVSPt=				NULL;
TH2D* 	fHistoTrueMesonCaloPhotonInvMassVSPt=				NULL;
TH2D* 	fHistoTrueMesonCaloConvPhotonInvMassVSPt=			NULL;
TH2D* 	fHistoTrueMesonCaloElectronInvMassVSPt=				NULL;
TH2D* 	fHistoTrueMesonMergedClusterInvMassVSPt=			NULL;
TH2D* 	fHistoTrueMesonMergedClusterPartConvInvMassVSPt=	NULL;
TH2D* 	fHistoTrueMesonCaloEMNonLeading=					NULL;

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
TString ObjectNameTrueSec;
TString ObjectNameTrueSecFromK0S;
TString ObjectNameMCPi0Acc;
TString ObjectNameMCEtaAcc;
TString ObjectNameMCOmegaAcc;
TString ObjectNameMCPi0;
TString ObjectNameMCPi0WOWeights;
TString ObjectNameMCEta;
TString ObjectNameMCEtaWOWeights;
TString ObjectNameMCOmega;
TString ObjectNameMCOmegaWOWeights;
TString ObjectNameTrueGGBck;
TString ObjectNameTrueContBck;
TString ObjectNameTrueAllBck;
TString ObjectNameTrueCaloPhoton;
TString ObjectNameTrueCaloConvPhoton;
TString ObjectNameTrueCaloElectron;
TString ObjectNameTrueCaloMerged;
TString ObjectNameTrueCaloEMNonLeading;
TString ObjectNameTrueCaloMergedPartConv;


Bool_t fAdvancedMesonQA = kFALSE;
Bool_t fEstimateTrainPileUp = kFALSE;

void ProcessEM(TH1D*,TH1D*,Double_t *);
void ProcessBckFitSubtraction(TH1D *fGammaGamma, Int_t i, Double_t * fPeakRange, Double_t *fFitRange, TString, TString, TString, TString);
void ProcessRatioSignalBackground(TH1D* , TH1D* );
void FillMassHistosArray(TH2D*,TH2D*);
void FillMassMCTrueMesonHistosArray(TH2D*);
void FillMassMCTrueMesonCaloPhotonHistosArray(TH2D*);
void FillMassMCTrueMesonCaloConvPhotonHistosArray(TH2D*);
void FillMassMCTrueMesonCaloElectronHistosArray(TH2D*);
void FillMassMCTrueMesonCaloEMNonLeadingHistosArray(TH2D*);
void FillMassMCTrueMesonCaloMergedClusterHistosArray(TH2D*);
void FillMassMCTrueMesonCaloMergedClusterPartConvHistosArray(TH2D*);
void FillMassMCTrueReweightedMesonHistosArray(TH2D*);
void FillMassMCTrueGGBckHistosArray(TH2D*);
void FillMassMCTrueContBckHistosArray(TH2D*);
void FillMassMCTrueAllBckHistosArray(TH2D*);
void FillMassMCTrueSecMesonHistosArray(TH2D*);
void FillMassMCTrueSecFromK0SMesonHistosArray(TH2D*);
TH1D* CalculateSecondaryFractions(TH1D* histoRawYield, TH1D* histoRawYieldSec, TString nameHistoFrac);
void CreatePtHistos();
void FillPtHistos();
void FitSubtractedInvMassInPtBins(TH1D * , Double_t *fMesonIntDeltaRangeFit, Int_t, Bool_t );  // Fits the Invariant Mass histos with a given function
void FitTrueInvMassInPtBins(TH1D * , Double_t *fMesonIntDeltaRangeFit, Int_t, Bool_t);  // Fits the Invariant Mass histos with a given function
void FitPeakPosInvMassInPtBins(TH1D * , Int_t, Bool_t );  // Fits the Invariant Mass histos with a given function
void FitCBSubtractedInvMassInPtBins(TH1D * , Double_t *fMesonIntDeltaRangeFit, Int_t, Bool_t , TString );  // Fits the Invariant Mass histos with a given function;
void FitWithPol2ForBG(TH1D*, Double_t*fMesonFitRangeCur, Int_t, Bool_t);
void ProduceBckProperWeighting(TList*, TList* );
void ProduceBckWithoutWeighting(TH2D *);
void IntegrateHistoInvMassStream(TH1D * , Double_t *);
void IntegrateHistoInvMass(TH1D * , Double_t *);
void IntegrateFitFunc(TF1 * , TH1D *, Double_t *);
void FillHistosArrayMC(TH1D* , TH1D*, TH1D* ,TString);
void CalculateMesonAcceptance();
TH1D* CalculateMesonEfficiency(TH1D*, TH1D**,TString,TString);
void SaveHistos(Int_t, TString, TString);
void SaveCorrectionHistos(TString , TString);
void Initialize(TString setPi0, Int_t);
void CalculateFWHM(TF1 *);
Double_t CrystalBallBck(Double_t *,Double_t *);
Double_t CrystalBall(Double_t *,Double_t *);
void Delete();
void SetCorrectMCHistogrammNames(TString, TString);

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
        if(mode == 40 || mode == 41 || mode == 44){
           fMidPt[0]                   = 1.5;
           fMidPt[1]                   = 2.5;
         }

        // Initialize peak range
        if(mode == 40 || mode == 41 || mode == 44){
             fPeakRange[0]             = 0.75;
             fPeakRange[1]             = 0.81;
         }

        // Initialze fit range
        if(mode == 40 || mode == 41 || mode == 44){
            fFitRange[0]             = 0.615;
            fFitRange[1]             = 0.89;
            fIntFixedRange[0]        = 0.615; // not yet implemented
            fIntFixedRange[1]        = 0.89;  // not yet implemented
        }

        // Initialize default BG fit range right & left
        if(mode == 40 || mode == 41 || mode == 44){
            fBGFitRange[0]             = 0.82;
            fBGFitRange[1]             = 0.89;
            fBGFitRangeLeft[0]         = 0.615;
            fBGFitRangeLeft[1]         = 0.69;
        }

        // Initialize default Plot range for meson
        if(mode == 40 || mode == 41 || mode == 44){
            fMesonPlotRange[0]         = 0.75;
            fMesonPlotRange[1]         = 0.79;
        }

        // Initialize default Plot default integration ranges
        if(mode == 40 || mode == 41 || mode == 44){
            fMesonIntDeltaRange[0]      = -0.04;
            fMesonIntDeltaRange[1]      =  0.05;
            fMesonIntDeltaRangeWide[0]  = -0.06;
            fMesonIntDeltaRangeWide[1]  = 0.05;
            fMesonIntDeltaRangeNarrow[0]= -0.03;
            fMesonIntDeltaRangeNarrow[1]= 0.03;
        }

        // Set meson mass ranges (here same for fitting and plotting)
         if(mode == 40 || mode == 41 || mode == 44){
             fMesonMassPlotRange[0]      = 0.65;
             fMesonMassPlotRange[1]      = 0.89;
             fMesonMassRange[0]          = 0.65;
             fMesonMassRange[1]          = 0.89;
         }

         // Set meson fit range
         if(mode == 40 || mode == 41 || mode == 44){
             fMesonFitRange[0]           = 0.65;
             fMesonFitRange[1]           = 0.85;
         }

         // Set remaining parameters for fitting
         if(mode == 40 || mode == 41 || mode == 44){
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
