
/****************************************************************************************************************************
****** 	provided by Gamma Conversion Group, PWG4, 														*****
******		Ana Marin, marin@physi.uni-heidelberg.de													*****
******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
*****************************************************************************************************************************/
/* Color_t	colorComb0005				= kRed+1; */
/* Color_t	colorComb0010				= kRed+1; */
/* Color_t	colorComb0510				= 807; */
/* Color_t	colorComb1020				= 800; */
/* Color_t	colorComb2040				= kGreen+2; */
/* Color_t	colorComb4060				= kCyan+2; */
/* Color_t	colorComb6080				= kBlue+1; */
//************************* Variables for styling combined plots inv Yield + Mass *******************************************
Color_t	colorComb0020				= kRed+1;
Color_t	colorComb2040				= 807;
Color_t	colorComb4060				=800;
Color_t	colorComb6080				= kGreen+2;
Color_t  colorComb80100         =kCyan+2;
Color_t  colorComb60100         =kBlue+1;

Color_t	colorComb0020Box			= kRed-6;
Color_t	colorComb2040Box			= 807-6;
Color_t	colorComb4060Box			=800-6;
Color_t	colorComb6080Box			= kGreen-6;
Color_t  colorComb80100Box      = kCyan-6;
Color_t  colorComb60100Box      = kBlue-6;

Style_t 	markerStylePCM		= 20 ;
Style_t 	markerStyleDalitz		= 21 ;
Style_t  markerStylePCMMC   = 24 ;
Style_t  markerStyleDalitzMC   = 25 ;
Color_t 	colorPCM 			   = kBlack;
Color_t 	colorPCMMC 			= kGray+1;
Color_t 	colorDalitz				= kRed+1;
Color_t 	colorDalitzMC				= kRed-7;
Color_t 	colorPCMSyst			= kGray;
Color_t 	colorDalitzSyst			= kRed-8;
Style_t 	fillStylePCM			= 0 ;
Style_t 	fillStyleDalitz			= 0;
Size_t 	markerSizeInvYield	= 1.5;
Size_t 	markerSizeMass			= 3.2;

//********************* Variables Combmon Spectrum ************************************************************************
Style_t 	markerStyleCommmonSpectrumpp 	= 20 ;
Style_t 	markerStyleCommmonSpectrum0020 	= 20 ;
Style_t 	markerStyleCommmonSpectrum2040 	= 21 ;
Style_t 	markerStyleCommmonSpectrum4060 	= 29 ;
Style_t 	markerStyleCommmonSpectrum6080 	= 33 ;
Style_t  markerStyleCommmonSpectrum80100 = 34 ;
Style_t  markerStyleCommmonSpectrum60100 = 34 ;

Style_t 	markerStyleCommmonSpectrumComp0020 	= 24 ;
Style_t 	markerStyleCommmonSpectrumComp2040 	= 25 ;
Style_t 	markerStyleCommmonSpectrumComp4060 	= 30 ;
Style_t 	markerStyleCommmonSpectrumComp6080 	= 27 ;
Style_t  markerStyleCommmonSpectrumComp80100 = 28 ;
Style_t  markerStyleCommmonSpectrumComp60100 = 28 ;

Size_t 	markerSizeCommonSpectrum0020 	= 2.;
Size_t 	markerSizeCommonSpectrum2040 	= 2.;
Size_t 	markerSizeCommonSpectrum4060 	= 2.5;
Size_t 	markerSizeCommonSpectrum6080 	= 2.5;
Size_t   markerSizeCommonSpectrum80100  = 2.;
Size_t   markerSizeCommonSpectrum60100  = 2.;
Size_t 	markerSizeSpectrum 			= 2.;

Style_t 	styleFitCommonSpectrum 		= 1;
Width_t 	widthLinesBoxes;
Width_t 	widthStatErrBars;
Width_t 	widthCommonFit;
Width_t 	widthCommonSpectrumBoxes;
Width_t 	widthCommonErrors;

Bool_t 		pictDrawingOptions[4] = 			{kFALSE, kFALSE, kFALSE, kTRUE};
TString 		date;
TString 		collisionSystemPP;		
TString 		collisionSystem0020;		
TString 		collisionSystem2040;		
TString 		collisionSystem4060;		
TString 		collisionSystem6080;
TString     collisionSystem80100;
TString     collisionSystem60100;
TString     collisionSystempPb;
TString 		forOutput;

TString 		nameFinalResDat;
TString 		prefix2;
Double_t 	minPtForFits;
Double_t 	mesonMassExpectPi0;

Double_t maxPtPi0pPb0020;
Double_t maxPtPi0pPb2040;
Double_t maxPtPi0pPb4060;
Double_t maxPtPi0pPb6080;
Double_t maxPtPi0pPb80100;
Double_t maxPtPi0pPb60100;
Double_t maxPtPi0pPb;

Double_t nColl0020;
Double_t nColl2040;
Double_t nColl4060;
Double_t nColl6080;
Double_t nColl80100;
Double_t nColl60100;
Double_t nCollpPb;


TFile* filePCMpPb;
TDirectory* directoryPCMPi0pPb;
TH1D* histoPCMNumberOfEventspPb;
TH1D* histoPCMYieldPi0pPb;
TH1D* histoPCMRAWYieldPi0pPb;
TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb;
TH1D* histoPCMMassPi0DatapPb;
TH1D* histoPCMMassPi0MCpPb;
TH1D* histoPCMWidthPi0DatapPb;
TH1D* histoPCMWidthPi0MCpPb;
Double_t nPCMEventpPb;
TDirectory* directoryPCMEtapPb;
TH1D* histoPCMYieldEtapPb;
TH1D* histoPCMRAWYieldEtapPb;
TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb;
TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb;
TH1D* histoPCMEtaPi0RatiopPb;
TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb;
TH1D* histoPCMMassEtaDatapPb;
TH1D* histoPCMMassEtaMCpPb;
TH1D* histoPCMWidthEtaDatapPb;
TH1D* histoPCMWidthEtaMCpPb;

TDirectory* directoryPCMPi0pPb0020;
TH1D* histoPCMNumberOfEvents0020;
TH1D* histoPCMYieldPi0pPb0020;
TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb0020;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb0020;
Double_t nPCMEventpPb0020;
TH1D* histoPCMMassPi0DatapPb0020;
TH1D* histoPCMMassPi0MCpPb0020;
TH1D* histoPCMWidthPi0DatapPb0020;
TH1D* histoPCMWidthPi0MCpPb0020;
TDirectory* directoryPCMEtapPb0020;
TH1D* histoPCMYieldEtapPb0020;
TH1D* histoPCMRAWYieldEtapPb0020;
TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb0020;
TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb0020;
TH1D* histoPCMEtaPi0RatiopPb0020;
TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb0020;
TH1D* histoPCMMassEtaDatapPb0020;
TH1D* histoPCMMassEtaMCpPb0020;
TH1D* histoPCMWidthEtaDatapPb0020;
TH1D* histoPCMWidthEtaMCpPb0020;

TDirectory* directoryPCMPi0pPb2040;
TH1D* histoPCMNumberOfEvents2040;
TH1D* histoPCMYieldPi0pPb2040;
TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb2040;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb2040;
Double_t nPCMEventpPb2040;
TH1D* histoPCMMassPi0DatapPb2040;
TH1D* histoPCMMassPi0MCpPb2040;
TH1D* histoPCMWidthPi0DatapPb2040;
TH1D* histoPCMWidthPi0MCpPb2040;
TDirectory* directoryPCMEtapPb2040;
TH1D* histoPCMYieldEtapPb2040;
TH1D* histoPCMRAWYieldEtapPb2040;
TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb2040;
TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb2040;
TH1D* histoPCMEtaPi0RatiopPb2040;
TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb2040;
TH1D* histoPCMMassEtaDatapPb2040;
TH1D* histoPCMMassEtaMCpPb2040;
TH1D* histoPCMWidthEtaDatapPb2040;
TH1D* histoPCMWidthEtaMCpPb2040;

TDirectory* directoryPCMPi0pPb4060;
TH1D* histoPCMNumberOfEvents4060;
TH1D* histoPCMYieldPi0pPb4060;
TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb4060;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb4060;
Double_t nPCMEventpPb4060;
TH1D* histoPCMMassPi0DatapPb4060;
TH1D* histoPCMMassPi0MCpPb4060;
TH1D* histoPCMWidthPi0DatapPb4060;
TH1D* histoPCMWidthPi0MCpPb4060;
TDirectory* directoryPCMEtapPb4060;
TH1D* histoPCMYieldEtapPb4060;
TH1D* histoPCMRAWYieldEtapPb4060;
TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb4060;
TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb4060;
TH1D* histoPCMEtaPi0RatiopPb4060;
TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb4060;
TH1D* histoPCMMassEtaDatapPb4060;
TH1D* histoPCMMassEtaMCpPb4060;
TH1D* histoPCMWidthEtaDatapPb4060;
TH1D* histoPCMWidthEtaMCpPb4060;

TDirectory* directoryPCMPi0pPb6080;
TH1D* histoPCMNumberOfEvents6080;
TH1D* histoPCMYieldPi0pPb6080;
TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb6080;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb6080;
TH1D* histoPCMMassPi0DatapPb6080;
TH1D* histoPCMMassPi0MCpPb6080;
TH1D* histoPCMWidthPi0DatapPb6080;
TH1D* histoPCMWidthPi0MCpPb6080;
Double_t nPCMEventpPb6080;
TDirectory* directoryPCMEtapPb6080;
TH1D* histoPCMYieldEtapPb6080;
TH1D* histoPCMRAWYieldEtapPb6080;
TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb6080;
TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb6080;
TH1D* histoPCMEtaPi0RatiopPb6080;
TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb6080;
TH1D* histoPCMMassEtaDatapPb6080;
TH1D* histoPCMMassEtaMCpPb6080;
TH1D* histoPCMWidthEtaDatapPb6080;
TH1D* histoPCMWidthEtaMCpPb6080;

TDirectory* directoryPCMPi0pPb80100;
TH1D* histoPCMNumberOfEvents80100;
TH1D* histoPCMYieldPi0pPb80100;
TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb80100;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb80100;
TH1D* histoPCMMassPi0DatapPb80100;
TH1D* histoPCMMassPi0MCpPb80100;
TH1D* histoPCMWidthPi0DatapPb80100;
TH1D* histoPCMWidthPi0MCpPb80100;
Double_t nPCMEventpPb80100;
TDirectory* directoryPCMEtapPb80100;
TH1D* histoPCMYieldEtapPb80100;
TH1D* histoPCMRAWYieldEtapPb80100;
TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb80100;
TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb80100;
TH1D* histoPCMEtaPi0RatiopPb80100;
TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb80100;
TH1D* histoPCMMassEtaDatapPb80100;
TH1D* histoPCMMassEtaMCpPb80100;
TH1D* histoPCMWidthEtaDatapPb80100;
TH1D* histoPCMWidthEtaMCpPb80100;

TDirectory* directoryPCMPi0pPb60100;
TH1D* histoPCMNumberOfEvents60100;
TH1D* histoPCMYieldPi0pPb60100;
TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb60100;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb60100;
TH1D* histoPCMMassPi0DatapPb60100;
TH1D* histoPCMMassPi0MCpPb60100;
TH1D* histoPCMWidthPi0DatapPb60100;
TH1D* histoPCMWidthPi0MCpPb60100;
Double_t nPCMEventpPb60100;
TDirectory* directoryPCMEtapPb60100;
TH1D* histoPCMYieldEtapPb60100;
TH1D* histoPCMRAWYieldEtapPb60100;
TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb60100;
TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb60100;
TH1D* histoPCMEtaPi0RatiopPb60100;
TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb60100;
TH1D* histoPCMMassEtaDatapPb60100;
TH1D* histoPCMMassEtaMCpPb60100;
TH1D* histoPCMWidthEtaDatapPb60100;
TH1D* histoPCMWidthEtaMCpPb60100;

TFile* fileDalitzpPb;
TDirectory* directoryDalitzPi0pPb;
TH1D* histoDalitzNumberOfEventspPb;
TH1D* histoDalitzYieldPi0pPb;
TH1D* histoDalitzRAWYieldPi0pPb;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPb;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrRAApPb;
TH1D* histoDalitzMassPi0DatapPb;
TH1D* histoDalitzMassPi0MCpPb;
TH1D* histoDalitzWidthPi0DatapPb;
TH1D* histoDalitzWidthPi0MCpPb;
Double_t nDalitzEventpPb;
TDirectory* directoryDalitzEtapPb;
TH1D* histoDalitzYieldEtapPb;
TH1D* histoDalitzRAWYieldEtapPb;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrpPb;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrRAApPb;
TH1D* histoDalitzEtaPi0RatiopPb;
TGraphAsymmErrors* graphDalitzEtaPi0RatioSysErrpPb;
TH1D* histoDalitzMassEtaDatapPb;
TH1D* histoDalitzMassEtaMCpPb;
TH1D* histoDalitzWidthEtaDatapPb;
TH1D* histoDalitzWidthEtaMCpPb;

TDirectory* directoryDalitzPi0pPb0020;
TH1D* histoDalitzNumberOfEvents0020;
TH1D* histoDalitzYieldPi0pPb0020;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPb0020;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrRAApPb0020;
Double_t nDalitzEventpPb0020;
TH1D* histoDalitzMassPi0DatapPb0020;
TH1D* histoDalitzMassPi0MCpPb0020;
TH1D* histoDalitzWidthPi0DatapPb0020;
TH1D* histoDalitzWidthPi0MCpPb0020;
TDirectory* directoryDalitzEtapPb0020;
TH1D* histoDalitzYieldEtapPb0020;
TH1D* histoDalitzRAWYieldEtapPb0020;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrpPb0020;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrRAApPb0020;
TH1D* histoDalitzEtaPi0RatiopPb0020;
TGraphAsymmErrors* graphDalitzEtaPi0RatioSysErrpPb0020;
TH1D* histoDalitzMassEtaDatapPb0020;
TH1D* histoDalitzMassEtaMCpPb0020;
TH1D* histoDalitzWidthEtaDatapPb0020;
TH1D* histoDalitzWidthEtaMCpPb0020;

TDirectory* directoryDalitzPi0pPb2040;
TH1D* histoDalitzNumberOfEvents2040;
TH1D* histoDalitzYieldPi0pPb2040;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPb2040;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrRAApPb2040;
Double_t nDalitzEventpPb2040;
TH1D* histoDalitzMassPi0DatapPb2040;
TH1D* histoDalitzMassPi0MCpPb2040;
TH1D* histoDalitzWidthPi0DatapPb2040;
TH1D* histoDalitzWidthPi0MCpPb2040;
TDirectory* directoryDalitzEtapPb2040;
TH1D* histoDalitzYieldEtapPb2040;
TH1D* histoDalitzRAWYieldEtapPb2040;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrpPb2040;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrRAApPb2040;
TH1D* histoDalitzEtaPi0RatiopPb2040;
TGraphAsymmErrors* graphDalitzEtaPi0RatioSysErrpPb2040;
TH1D* histoDalitzMassEtaDatapPb2040;
TH1D* histoDalitzMassEtaMCpPb2040;
TH1D* histoDalitzWidthEtaDatapPb2040;
TH1D* histoDalitzWidthEtaMCpPb2040;

TDirectory* directoryDalitzPi0pPb4060;
TH1D* histoDalitzNumberOfEvents4060;
TH1D* histoDalitzYieldPi0pPb4060;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPb4060;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrRAApPb4060;
Double_t nDalitzEventpPb4060;
TH1D* histoDalitzMassPi0DatapPb4060;
TH1D* histoDalitzMassPi0MCpPb4060;
TH1D* histoDalitzWidthPi0DatapPb4060;
TH1D* histoDalitzWidthPi0MCpPb4060;
TDirectory* directoryDalitzEtapPb4060;
TH1D* histoDalitzYieldEtapPb4060;
TH1D* histoDalitzRAWYieldEtapPb4060;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrpPb4060;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrRAApPb4060;
TH1D* histoDalitzEtaPi0RatiopPb4060;
TGraphAsymmErrors* graphDalitzEtaPi0RatioSysErrpPb4060;
TH1D* histoDalitzMassEtaDatapPb4060;
TH1D* histoDalitzMassEtaMCpPb4060;
TH1D* histoDalitzWidthEtaDatapPb4060;
TH1D* histoDalitzWidthEtaMCpPb4060;

TDirectory* directoryDalitzPi0pPb6080;
TH1D* histoDalitzNumberOfEvents6080;
TH1D* histoDalitzYieldPi0pPb6080;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPb6080;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrRAApPb6080;
TH1D* histoDalitzMassPi0DatapPb6080;
TH1D* histoDalitzMassPi0MCpPb6080;
TH1D* histoDalitzWidthPi0DatapPb6080;
TH1D* histoDalitzWidthPi0MCpPb6080;
Double_t nDalitzEventpPb6080;
TDirectory* directoryDalitzEtapPb6080;
TH1D* histoDalitzYieldEtapPb6080;
TH1D* histoDalitzRAWYieldEtapPb6080;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrpPb6080;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrRAApPb6080;
TH1D* histoDalitzEtaPi0RatiopPb6080;
TGraphAsymmErrors* graphDalitzEtaPi0RatioSysErrpPb6080;
TH1D* histoDalitzMassEtaDatapPb6080;
TH1D* histoDalitzMassEtaMCpPb6080;
TH1D* histoDalitzWidthEtaDatapPb6080;
TH1D* histoDalitzWidthEtaMCpPb6080;

TDirectory* directoryDalitzPi0pPb80100;
TH1D* histoDalitzNumberOfEvents80100;
TH1D* histoDalitzYieldPi0pPb80100;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPb80100;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrRAApPb80100;
TH1D* histoDalitzMassPi0DatapPb80100;
TH1D* histoDalitzMassPi0MCpPb80100;
TH1D* histoDalitzWidthPi0DatapPb80100;
TH1D* histoDalitzWidthPi0MCpPb80100;
Double_t nDalitzEventpPb80100;
TDirectory* directoryDalitzEtapPb80100;
TH1D* histoDalitzYieldEtapPb80100;
TH1D* histoDalitzRAWYieldEtapPb80100;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrpPb80100;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrRAApPb80100;
TH1D* histoDalitzEtaPi0RatiopPb80100;
TGraphAsymmErrors* graphDalitzEtaPi0RatioSysErrpPb80100;
TH1D* histoDalitzMassEtaDatapPb80100;
TH1D* histoDalitzMassEtaMCpPb80100;
TH1D* histoDalitzWidthEtaDatapPb80100;
TH1D* histoDalitzWidthEtaMCpPb80100;

TDirectory* directoryDalitzPi0pPb60100;
TH1D* histoDalitzNumberOfEvents60100;
TH1D* histoDalitzYieldPi0pPb60100;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPb60100;
TGraphAsymmErrors* graphDalitzYieldPi0SysErrRAApPb60100;
TH1D* histoDalitzMassPi0DatapPb60100;
TH1D* histoDalitzMassPi0MCpPb60100;
TH1D* histoDalitzWidthPi0DatapPb60100;
TH1D* histoDalitzWidthPi0MCpPb60100;
Double_t nDalitzEventpPb60100;
TDirectory* directoryDalitzEtapPb60100;
TH1D* histoDalitzYieldEtapPb60100;
TH1D* histoDalitzRAWYieldEtapPb60100;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrpPb60100;
TGraphAsymmErrors* graphDalitzYieldEtaSysErrRAApPb60100;
TH1D* histoDalitzEtaPi0RatiopPb60100;
TGraphAsymmErrors* graphDalitzEtaPi0RatioSysErrpPb60100;
TH1D* histoDalitzMassEtaDatapPb60100;
TH1D* histoDalitzMassEtaMCpPb60100;
TH1D* histoDalitzWidthEtaDatapPb60100;
TH1D* histoDalitzWidthEtaMCpPb60100;

Double_t ptBinningPCMpPb[31] =               {0.0,0.4,0.5,0.6,0.7,
                                            0.8,1.0,1.2,1.4,1.6,
                                            1.8,2.0,2.2,2.4,2.6,
                                            2.8,3.0,3.2,3.4,3.6,
                                            3.8,4.0,4.5,5.0,5.5,
                                            6.0,7.0,8.0,10.0,12.0, 14.0};
/* Double_t ptBinningPCMpPb [30] = {0.0,0.4,0.5,0.6,0.7, */
/*                               0.8,1.0,1.2,1.4,1.6, */
/*                               1.8,2.0,2.2,2.4,2.6, */
/*                               2.8,3.0,3.2,3.4,3.6, */
/*                               3.8,4.0,4.5,5.0,5.5, */
/*                               6.0,7.0,8.0,10.0,15.0}; */
/* Double_t ptBinningDalitzpPb [30] = {0.0,0.4,0.5,0.6,0.7, */
/*                               0.8,1.0,1.2,1.4,1.6, */
/*                               1.8,2.0,2.2,2.4,2.6, */
/*                               2.8,3.0,3.2,3.4,3.6, */
/*                               3.8,4.0,4.5,5.0,5.5, */
/*                               6.0,7.0,8.0,10.0,15.0}; */

Double_t ptBinningDalitzpPb[19] = {0.0,0.5,0.6,0.8,1.0,
				   1.2,1.4,1.6,1.8,2.0,
				   2.2,2.4,2.8,3.2,3.6,
				   4.0,6.0,10.,15.};
