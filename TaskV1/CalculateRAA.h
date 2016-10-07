// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

// void CorrectConvPHOS(TH1D*, Double_t);

TString 	collisionSystem;
TString 	fileNameSysErr;
TFile* 	fileConversionsPi0PP;
TFile* 	fileConversionsEtaPP;
TF1* 		fitInvCrossSectionPi0Comb2760GeV;
TF1* 		fitInvCrossSectionEtaComb2760GeV;

Double_t		minPtForFits = 		0.2;
Double_t		minPtForFitsEta = 0.4;
Double_t 	deltaEta;
TString 		rapidityRange;

Float_t 		pictDrawingCoordinatesMass[9] = 	{0.4, 0.8, 0.75, 0.04, 0.2,0.68, 0.12, 0.035,0};
Float_t 		pictDrawingCoordinates[9] = 		{0.55, 0.8, 0.25, 0.04, 0.7,0.5, 0.12, 0.035,0};
Bool_t 			pictDrawingOptions[4] = 			{kFALSE, kFALSE, kFALSE, kTRUE};

TString 		prefix2;
Bool_t	 	kLevy;
Bool_t 		kHag;
Double_t 	mesonMassExpect = 			0;

TFile* 	file;
TH1D*		histoCorrectedYield;;
TH1D*    histoMCInput;
TH1D*    histoMCInputWOWeight = NULL;
TH1D*    histoMCInputWeights = NULL;
TH1D*    histoAccept;
TH1D*    histoEffi;
TH1D*    histoEffiFitted;
TH1D*    histoRawYield; 
TH1D*		histoFWHMMeson;
TH1D*		histoMassMeson;
TH1D*		histoTrueFWHMMeson;
TH1D*		histoTrueMassMeson;
TH1F*		histoEventQualtity;
TH1D* 	histoMassMesonMinusExp;
TH1D* 	histoTrueMassMesonMinusExp;
TH1D* 	histoMesonSignalFullPtInvMass;
TH1D* 	histoFWHMMesonMeV;
TH1D*		histoTrueFWHMMesonMeV;

TString 		date;
Float_t 		nEvt;
Double_t 	maxPt;

TString	nameFinalResDat;
fstream 	fileFinalResults;

//************ Fitting variables Pi0 *****************
TH1D *		histoFitting;
Bool_t	 	kLevySucc = 		kFALSE;
Bool_t	 	kRaduSucc = 		kFALSE;
Bool_t 		kHagSucc = 			kFALSE;
Bool_t 		kExpSucc = 			kFALSE;
Bool_t 		kPowSucc = 			kFALSE;
Bool_t 		kBoltzSucc = 		kFALSE;
Bool_t 		kModPowSucc = 		kFALSE;
TF1 *		fitPtLevy;
TF1 *		fitPtRadu;
TF1 *		fitPtHagedorn;
TF1 *		fitPtHagedorn2;
TF1 *		fitPtBoltzmann;
TF1 *		fitPtExp;
TF1 *		fitPtPowerlaw;
TF1 *		fitPtModPowerlaw;
TH1D *		histoRatioFitLevy;
TH1D *		histoRatioFitRadu;
TH1D *		histoRatioFitHag;
TH1D *		histoRatioFitBoltz;
TH1D *		histoRatioFitExp;
TH1D *		histoRatioFitPow;
TH1D *		histoRatioFitModPow;

//************ Bin shift variables pi0 +***********************
TH1D *		histoRecBinShift;
TH1D *		histoCorrYieldBinShifted;
Int_t 		nb;
TF1 *		fitBinShifting;
const char* 	legendEntryFunction;
Double_t 		globalRatio = 			0.;
Double_t 		testGlobalRatio = 		1.;
Double_t 		color = 				1;
TH1D *			histoRatioBinShifted;

//*************** Variables Systematic errors *******************
TGraphAsymmErrors* 	graphCorrectedYieldSysErr;
TGraphAsymmErrors* 	graphCorrectedYieldSysErrA;
TGraphAsymmErrors* 	graphCorrectedYieldSysErrBinShifted;
TGraphAsymmErrors* 	graphCorrectedYieldSysErrABinShifted;
TGraphAsymmErrors* 	graphCorrectedYieldSysRatioFitConv;
TGraphAsymmErrors* 	graphCorrectedYieldStatPlusSys;
TGraphAsymmErrors* 	graphCorrectedYieldStatPlusSysBinShifted;
TGraphAsymmErrors* 	graphCorrectedYieldSysErrForFit;

ifstream 		fileSysErr;
Int_t 		nPoints = 			0;
Double_t 		relSystErrorUp[50];
Double_t 		relSystErrorDown[50];
Double_t 		relSystErrorWOMaterialUp[50];
Double_t 		relSystErrorWOMaterialDown[50];

TF1 *		fitPtLevySysErr;
TF1 *		fitPtHagedornSysErr;
TF1 *		fitPtBoltzmannSysErr;
TF1 *		fitPtExpSysErr;
TF1 *		fitPtPowerlawSysErr;
TF1 *		fitPtModPowerlawSysErr;

TString 		forOutput;
TH1D* 		histoNumberOfEvents;
TString centralityString;

// Double_t factorToInel = 1/1.12;		// this factor is multiplied with Raa and comes from trigger inelastic effiency for pp
Double_t	xSection2760GeVpp = 		55.43*1e-3; //sigmaOR
Double_t	xSection2760GeVErrpp = 	3.9;
// Double_t recalcBarn = 			1e12; //NLO in pbarn!!!!

Double_t xPtLimitsPbPbPHOS[21] =  {	0.4,0.6,0.8,1.0,1.2,
												1.4,1.6,1.8,2.0,2.2,
												2.4,2.6,3.0,3.5,4.0,
												5.0,6.0,8.0,10.,12.,
												15.};
Double_t nColl;
Double_t nCollErr;
TGraphAsymmErrors* graphInvSectionCombPi02760GeV;
TGraphAsymmErrors* graphInvSectionPCMPi02760GeV;
TGraphAsymmErrors* graphInvSectionPCMSysPi02760GeV;

TGraphAsymmErrors* graphInvSectionCombEta2760GeV;
TGraphAsymmErrors* graphInvSectionPCMEta2760GeV;
TGraphAsymmErrors* graphInvSectionPCMSysEta2760GeV;