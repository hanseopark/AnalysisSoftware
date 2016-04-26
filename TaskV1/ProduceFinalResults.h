// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

// void CorrectConvPHOS(TH1D*, Double_t);

TString 		collisionSystem;
const char* 	fileNameCaloPhos;
const char* 	fileNameChargedSpectra;
const char* 	fileNameChargedExpectation;
const char* 	fileNameSysErrEta;

Bool_t 		thesis =		 	kFALSE;

Bool_t 		conference =		 	kFALSE;
Double_t		minPtForFits = 		0.2;
Double_t		minPtForFitsEta = 		0.4;
TString 		rapitdityCutNumber;
Double_t 		deltaEta;
TString 		rapidityRange;

Float_t 		pictDrawingCoordinatesAcc[9] = 	{0.63, 0.2, 0.30, 0.04, 0.7,0.33, 0.12, 0.035,0};
Float_t 		pictDrawingCoordinatesPi0Eta[9] = 	{0.64, 0.82, 0.3, 0.04, 0.13,0.77, 0.08, 0.04,0};
Float_t 		pictDrawingCoordinatesAccEff[9] = 	{0.65, 0.25, 0.4, 0.04, 0.18,0.15, 0.1, 0.035,0};
Float_t 		pictDrawingCoordinatesInv[9] = 	{0.43, 0.825, 0.15, 0.04, 0.825, 0.78, 0.08, 0.035,0};
Float_t 		pictDrawingCoordinatesMassFWHM[9] = {0.3, 0.8, 0.55, 0.04, 0.16, 0.75, 0.08, 0.035,0};
Float_t 		pictDrawingCoordinatesMassFWHM2[9] = {0.3, 0.83, 0.55, 0.04, 0.12, 0.78, 0.08, 0.035,0};
Float_t 		pictDrawingCoordinatesMassFWHM3[9] = {0.3, 0.83, 0.55, 0.04, 0.15, 0.8, 0.08, 0.03,0};
Float_t 		pictDrawingCoordinatesMass[9] = 	{0.4, 0.8, 0.75, 0.04, 0.2,0.68, 0.18, 0.035,0};
Float_t 		pictDrawingCoordinatesFWHM[9] = 	{0.4, 0.8, 0.70, 0.04, 0.2,0.68, 0.12, 0.035,0};
Float_t 		pictDrawingCoordinatesSign[9] = 	{0.63, 0.8, 0.30, 0.04, 0.7, 0.47, 0.12, 0.035,0};
Float_t 		pictDrawingCoordinatesSB[9] = 	{0.4, 0.8, 0.70, 0.04, 0.2,0.68, 0.12, 0.035,0};
Float_t 		pictDrawingCoordinatesCombineMeasCross[9] = 	{0.64, 0.86, 0.27, 0.02, 0.71,0.73, 0.09, 0.026,0};
//Float_t 		pictDrawingCoordinatesCombineMeasCross[9] = 	{0.67, 0.89, 0.25, 0.02, 0.72,0.75, 0.15, 0.025,0};
//Float_t 		pictDrawingCoordinatesOnlyCTSMeasCross[9] = 	{0.65, 0.86, 0.3, 0.02, 0.72,0.70, 0.15, 0.025,0};
Float_t 		pictDrawingCoordinatesOnlyCTSMeasCross[9] = 	{0.64, 0.86, 0.3, 0.02, 0.72,0.70, 0.1, 0.025,0};
Float_t 		pictDrawingCoordinatesCombineMeas[9] = 	{0.68, 0.8, 0.2, 0.04, 0.75,0.65, 0.1, 0.03,0};
Float_t 		pictDrawingCoordinatesOnlySpectrum[9] = 	{0.7, 0.86, 0.25, 0.02, 0.28,0.2, 0.1, 0.028,0};
Float_t 		pictDrawingCoordinatesOnlyRatio[9] = 		{0.7, 0.2, 0.25, 0.02, 0.85,0.55, 0.1, 0.05,0};
Float_t 		pictDrawingCoordinatesOnlySpectrumOnlyCTS[9] = 	{0.7, 0.84, 0.35, 0.02, 0.2,0.15, 0.18, 0.028,0};
Float_t 		pictDrawingCoordinates[9] = 		{0.55, 0.8, 0.25, 0.04, 0.7,0.5, 0.12, 0.035,0};
Float_t 		pictDrawingCoordinatesInvX[9] = 		{0.7, 0.8, 0.35, 0.04, 0.75,0.6, 0.15, 0.035,0};
Float_t 		pictDrawingCoordinatesRat[9] = 	{0.55, 0.8, 0.20, 0.04, 0.7,0.5, 0.18, 0.035,0};
Bool_t 		pictDrawingOptions[4] = 			{kFALSE, kFALSE, kFALSE, kTRUE};

TString 		prefix2;
Bool_t	 	kLevy;
Bool_t 		kHag;
Double_t 		scaling = 					1./(2.*TMath::Pi());
Double_t 		mesonMassExpectPi0 = 			0;
Double_t 		mesonMassExpectEta = 			0;

TString 		fileNamePi0ch;
TFile* 		filePi0;
TH1D*		histoCorrectedYieldPi0;
TH1D*		histoUncorrectedYieldPi0;
TH1D*		histoFWHMMesonPi0;
TH1D*		histoMassMesonPi0;
TH1D*		histoSBMesonPi0;
TH1D*		histoSignMesonPi0;
TH1D*		histoAccPi0;
TH1D*		histoTrueEffPtPi0;
TH1D*		histoTrueFWHMMesonPi0;
TH1D*		histoTrueMassMesonPi0;
TH1D*    histoMCInputPi0;
TH1D*		histoCorrectedYieldMtPi0;
TH1D*		histoCorrectedYieldXtPi0;
TH1F*		histoEventQualtityPi0;
TH1D* 		histoMassMesonPi0MinusExp;
TH1D* 		histoTrueMassMesonPi0MinusExp;
TH1D* 		histoMesonSignalFullPtInvMass;
TH1D* 		histoFWHMMesonPi0MeV;
TH1D*		histoTrueFWHMMesonPi0MeV;

TLatex * 		textInvMass;
Double_t 		widthMeV;
TString 		date;

Float_t 		nEvt;
TFile* 		fileEta;
TH1D *		histoCorrectedYieldEta;
TH1D *		histoCorrectedYieldMtEta;
TH1D *		histoCorrectedYieldXtEta;
TH1D *		histoUnCorrectedYieldEta;
TH1D *		histoFWHMMesonEta;
TH1D *		histoMassMesonEta;
TH1D *		histoSBMesonEta;
TH1D *		histoSignMesonEta;
TH1D *		histoAccEta;
TH1D *		histoTrueEffPtEta;
TH1D *		histoTrueFWHMMesonEta;
TH1D *		histoTrueMassMesonEta;
TH1D*       histoMCInputEta;
TH1D * 		histoMassMesonEtaMinusExp;
TH1D * 		histoTrueMassMesonEtaMinusExp;
TH1D* 		histoFWHMMesonEtaMeV;
TH1D*		histoTrueFWHMMesonEtaMeV;
TH1D* 		histoDeltaPtEta;

TFile* 		filePhos;
TH1D *		histoPi0Phos = 		0x00;
TH1D *		histoPi0PhosSys;

Double_t maxPtPi0;
Double_t maxPtEta;
Double_t maxPtEtaPHOS;

TH1D* 		histoEtaToPi0Pythia7TeV;
TH1D* 		histoEtaToPi0Phojet7TeV;

TH1D* 		histoEtaToPi0Pythia900GeV;
TH1D* 		histoEtaToPi0Phojet900GeV;

TH1D *		histoEtaToPi0Phojet;
TH1D *		histoEtaToPi0Pythia;
TH1D * 		histoRatioEtaPi0;
TH1D * 		histoRatioEtaPi0BinShifted;
//Ratio as TGraphAsymError with +20 - 40
Double_t 		ratioXValue[50];
Double_t 		ratioYValue[50];
Double_t 		ratioYValueBinShifted[50];
Double_t 		ratioXError[50];
Double_t 		ratioSysUpError[50];
Double_t 		ratioSysDownError[50];
Double_t 		ratioSysUpErrorBinShifted[50];
Double_t 		ratioSysDownErrorBinShifted[50];


const char *	nameFinalResDat;
fstream 		fileFinalResults;

//************ Fitting variables Pi0 *****************
TH1D *		histoFittingPi0;
Bool_t	 	kLevySuccPi0 = 		kFALSE;
Bool_t 		kHagSuccPi0 = 			kFALSE;
Bool_t 		kExpSuccPi0 = 			kFALSE;
Bool_t 		kPowSuccPi0 = 			kFALSE;
Bool_t 		kBoltzSuccPi0 = 		kFALSE;
Bool_t 		kModPowSuccPi0 = 		kFALSE;
Bool_t          kTwoCompModelSuccPi0 =          kFALSE;
TF1 *		fitPtLevyPi0;
TF1 *		fitPtHagedornPi0;
TF1 *		fitPtHagedornPi02;
TF1 *		fitPtBoltzmannPi0;
TF1 *		fitPtExpPi0;
TF1 *		fitPtPowerlawPi0;
TF1 *		fitPtModPowerlawPi0;
TF1 *           fitPtTwoCompModelPi0;
TH1D *		histoRatioFitLevyPi0;
TH1D *		histoRatioFitHagPi0;
TH1D *		histoRatioFitBoltzPi0;
TH1D *		histoRatioFitExpPi0;
TH1D *		histoRatioFitPowPi0;
TH1D *		histoRatioFitModPowPi0;
TH1D *          histoRatioFitTwoCompModelPi0;

//************ Bin shift variables pi0 +***********************
TH1D *		histoRecBinShiftPi0;
TH1D *		histoPi0CorrYieldBinShifted;
Int_t 		nb;
TF1 *		fitBinShifting;
const char* 	legendEntryFunction;
Double_t 		globalRatio = 			0.;
Double_t 		testGlobalRatio = 		1.;
Double_t 		color = 				1;
TH1D *		histoRatioBinShiftedPi0;

//************* NLO reading variables Pi0 ********************
Double_t 		ptNLOPi0[100];
Double_t 		muHalf[100];
Double_t 		muOne[100];
Double_t 		muTwo[100];
Double_t 		ptNLOBKKPi0[100];
Double_t 		muTwoBKK[100];

Double_t 		ptNLODSSPi0[100];
Double_t		energyDSS[100];
Double_t 		muTwoDSS[100];

Double_t 		muXError[100];
Double_t 		muHalfError[100];
Double_t 		muOneError[100]; 
Double_t 		muTwoError[100];

ifstream 		in;
ifstream 		inBKK;
ifstream 		inDSS;

Int_t       nlinesNLO               = 0;
Int_t       nlinesNLOBKK            = 0;
Int_t       nlinesNLODSS            = 0;
Double_t    xSection                = 0; 


//*************** Fitting Eta variables **********************
TH1D *		histoFittingEta;
Bool_t 		kLevySuccEta = 		kFALSE;
Bool_t 		kHagSuccEta = 			kFALSE;
Bool_t 		kPowSuccEta = 			kFALSE;
Bool_t 		kModPowSuccEta = 		kFALSE;
Bool_t          kTwoCompSuccEta =               kFALSE;
TF1 *		fitPtLevyEta;
TF1 *		fitPtHagedornEta;
TF1 *		fitPtPowerlawEta;
TF1 *		fitPtModPowerlawEta;
TF1 *           fitPtTwoCompModelEta;
TH1D *		histoRatioFitLevyEta;
TH1D *		histoRatioFitHagEta;
TH1D *		histoRatioFitPowEta;
TH1D *		histoRatioFitModPowEta;
TH1D*           histoRatioFitTwoCompEta;

//*************** Binshift Eta variables **********************
TH1D *		histoRecBinShiftEta;
Int_t 		nbEta; 
TH1D *		histoEtaCorrYieldBinShifted;
TF1 *		fitBinShiftingEta;
const char* legendEntryFunctionEta;
Double_t 	globalRatioEta = 		0.;
Double_t 	testGlobalRatioEta = 	1.;
TH1D *		histoRatioBinShiftedEta;


//*************** Variables Systematic errors *******************
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErr;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErrA;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErrBinShifted;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErrABinShifted;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysRatioFitConv;
TGraphAsymmErrors* 	graphCorrectedYieldPi0StatPlusSys;
TGraphAsymmErrors* 	graphCorrectedYieldPi0StatPlusSysBinShifted;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErrForFit;
TH1D*			histoInvCrossSectionPi0;
TH1D*			histoInvCrossSectionPi0Xt;
TGraphAsymmErrors* 	graphInvCrossSectionSysPi0;
TGraphAsymmErrors* 	graphInvCrossSectionSysAPi0;
TGraphErrors* 	        graphInvCrossSectionStaPi0;

TH1D* 			histoPi0CorrYieldBinShiftedMt;
TGraphAsymmErrors*	graphCorrectedYieldPi0SysErrMt;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErrBinShiftedMt;
TGraphAsymmErrors*	graphCorrectedYieldPi0StatPlusSysMt;
TGraphAsymmErrors* 	graphCorrectedYieldPi0StatPlusSysBinShiftedMt;

TH1D* 			histoPi0CorrYieldBinShiftedXt;
TGraphAsymmErrors*	graphCorrectedYieldPi0SysErrXt;
TGraphAsymmErrors*	graphCorrectedYieldPi0SysErrBinShiftedXt;
TGraphAsymmErrors*	graphCorrectedYieldPi0StatPlusSysXt;
TGraphAsymmErrors*	graphCorrectedYieldPi0StatPlusSysBinShiftedXt;

TGraphAsymmErrors* 	graphCorrectedYieldEtaSysErr;
TGraphAsymmErrors* 	graphCorrectedYieldEtaSysErrA;
TGraphAsymmErrors* 	graphCorrectedYieldEtaSysErrBinShifted;
TGraphAsymmErrors* 	graphCorrectedYieldEtaSysErrABinShifted;
TGraphAsymmErrors* 	graphCorrectedYieldEtaSysRatioFitConv;
TGraphAsymmErrors* 	graphCorrectedYieldEtaStatPlusSys;
TGraphAsymmErrors* 	graphCorrectedYieldEtaStatPlusSysBinShifted;
TGraphAsymmErrors* 	graphCorrectedYieldEtaSysErrForFit;
TH1D	*			histoInvCrossSectionEta;
TH1D	*			histoInvCrossSectionEtaXt;
TGraphAsymmErrors* 	graphInvCrossSectionSysEta;
TGraphAsymmErrors* 	graphInvCrossSectionSysAEta;
TGraphErrors* 	        graphInvCrossSectionStaEta;

TH1D* 			histoEtaCorrYieldBinShiftedMt;
TGraphAsymmErrors*	graphCorrectedYieldEtaSysErrMt;
TGraphAsymmErrors* 	graphCorrectedYieldEtaSysErrBinShiftedMt;
TGraphAsymmErrors*	graphCorrectedYieldEtaStatPlusSysMt;
TGraphAsymmErrors* 	graphCorrectedYieldEtaStatPlusSysBinShiftedMt;

TH1D* 			histoEtaCorrYieldBinShiftedXt;
TGraphAsymmErrors*	graphCorrectedYieldEtaSysErrXt;
TGraphAsymmErrors*	graphCorrectedYieldEtaSysErrBinShiftedXt;
TGraphAsymmErrors*	graphCorrectedYieldEtaStatPlusSysXt;
TGraphAsymmErrors*	graphCorrectedYieldEtaStatPlusSysBinShiftedXt;

ifstream 		fileSysErrPi0;
Int_t 		nPointsPi0 = 			0;
Double_t 		relSystErrorPi0Up[50];
Double_t 		relSystErrorPi0Down[50];
Double_t 		relSystErrorWOMaterialPi0Up[50];
Double_t 		relSystErrorWOMaterialPi0Down[50];
Double_t 		systErrorPi0Up[50];
Double_t 		systErrorPi0Down[50];
Double_t 		systErrorPi0BinShiftedUp[50];
Double_t 		systErrorPi0BinShiftedDown[50];
Double_t 		xValueCorrYieldPi0[50];
Double_t 		yValueCorrYieldPi0[50];
Double_t 		yErrorCorrYieldPi0[50];
Double_t 		xErrorCorrYieldPi0[50];
Double_t 		yErrorCombinedPi0Up[50];
Double_t 		yErrorCombinedPi0Down[50];
Double_t 		xValueCorrYieldPi0BinShifted[50];
Double_t 		yValueCorrYieldPi0BinShifted[50];
Double_t 		yErrorCorrYieldPi0BinShifted[50];
Double_t 		xErrorCorrYieldPi0BinShifted[50];
Double_t 		yErrorBinShiftedCombinedPi0Up[50];
Double_t 		yErrorBinShiftedCombinedPi0Down[50];
Double_t 		sysErrPi0RatioFitConvUp[50];
Double_t 		sysErrPi0RatioFitConvDown[50];
Double_t 		yValuePi0RatioFitConv[50];
ifstream 		fileSysErrEta;
Int_t 		nPointsEta = 0;
Double_t 		relSystErrorEtaUp[50];
Double_t 		relSystErrorEtaDown[50];
Double_t 		relSystErrorWOMaterialEtaUp[50];
Double_t 		relSystErrorWOMaterialEtaDown[50];
Double_t 		systErrorEtaUp[50];
Double_t 		systErrorEtaDown[50];
Double_t 		systErrorEtaBinShiftedUp[50];
Double_t 		systErrorEtaBinShiftedDown[50];
Double_t 		xValueCorrYieldEta[50];
Double_t 		yValueCorrYieldEta[50];
Double_t 		yErrorCorrYieldEta[50];
Double_t 		xErrorCorrYieldEta[50];
Double_t 		yErrorCombinedEtaUp[50];
Double_t 		yErrorCombinedEtaDown[50];
Double_t 		xValueCorrYieldEtaBinShifted[50];
Double_t 		yValueCorrYieldEtaBinShifted[50];
Double_t 		yErrorCorrYieldEtaBinShifted[50];
Double_t 		xErrorCorrYieldEtaBinShifted[50];
Double_t 		yErrorBinShiftedCombinedEtaUp[50];
Double_t 		yErrorBinShiftedCombinedEtaDown[50];
Double_t 		sysErrEtaRatioFitConvUp[50];
Double_t 		sysErrEtaRatioFitConvDown[50];
Double_t 		yValueEtaRatioFitConv[50];

TF1 *		fitPtLevyPi0SysErr;
TF1 *		fitPtHagedornPi0SysErr;
TF1 *		fitPtBoltzmannPi0SysErr;
TF1 *		fitPtExpPi0SysErr;
TF1 *		fitPtPowerlawPi0SysErr;
TF1 *		fitPtModPowerlawPi0SysErr;

TF1 *		fitPtLevyEtaSysErr;
TF1 *		fitPtHagedornEtaSysErr;
TF1 *		fitPtBoltzmannEtaSysErr;
TF1 *		fitPtExpEtaSysErr;
TF1 *		fitPtPowerlawEtaSysErr;
TF1 *		fitPtModPowerlawEtaSysErr;


TString 		forOutput;
TH1D* 		histoNumberOfEvents;

TH1D*		histoPi0ToChargedPhojet;
TH1D*		histoPi0ToChargedPythia;
TGraphAsymmErrors* graphSystErrRatio;
TGraphAsymmErrors* graphSystErrRatioBinShifted;

Double_t fPileUpCorrectionConv7TeV = 1-0.0105;
