/****************************************************************************************************************************
****** 	provided by Gamma Conversion Group, PWG4, 														*****
******		Ana Marin, marin@physi.uni-heidelberg.de													*****
******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
*****************************************************************************************************************************/

//************************* Variables for styling combined plots inv Yield + Mass *******************************************
Color_t	colorComb0005				= kRed+1;
Color_t	colorComb0010				= kRed+1;
Color_t	colorComb0510				= 807;
Color_t	colorComb1020				= 800;
Color_t	colorComb2040				= kGreen+2;
Color_t	colorComb4060				= kCyan+2;
Color_t	colorComb6080				= kBlue+1;
Color_t	colorWHDG0005				= kRed-4;
Color_t	colorWHDG0510				= 807+1;
Color_t	colorWHDG1020				= kYellow-6;
Color_t	colorWHDG2040				= kGreen-3;
Color_t	colorWHDG4060				= kCyan-3;
Color_t	colorWHDG6080				= kBlue-3;
Color_t	colorXiao0005				= kRed+3;
Color_t	colorXiao0510				= 807+2;
Color_t	colorXiao1020				= 800+1;
Color_t	colorXiao2040				= kGreen+3;
Color_t	colorXiao4060				= kCyan+3;
Color_t	colorXiao6080				= kBlue+3;
Color_t  colorEPOS               = kBlack;
Color_t  colorKopeliovichHydro   = kBlue+2;
Color_t  colorKopeliovichELoss   = kGreen+4;
Color_t  colorKopeliovichComb    = kCyan+2;
Color_t	colorEPOS0005				= kRed+4;
Color_t	colorEPOS0510				= 807-2;
Color_t	colorEPOS1020				= 800-3;
Color_t	colorEPOS2040				= kGreen+4;
Color_t	colorEPOS4060				= kCyan+4;
Color_t	colorEPOS6080				= kBlue+4;
Color_t  colorEPOSFill0005       = kRed-6;
Color_t  colorEPOSFill0510       = kOrange;
Color_t  colorEPOSFill1020       = kYellow-6;
Color_t  colorEPOSFill2040       = kGreen-5;
Color_t  colorEPOSFill4060       = kCyan-5;
Color_t  colorEPOSFill6080       = kBlue-6;

Style_t  styleLineNLOMuHalf         = 8;
Style_t  styleLineNLOMuOne       = 7;
Style_t  styleLineNLOMuTwo       = 4;

Style_t  styleKopeliovichHydro   = 8;
Style_t  styleKopeliovichELoss   = 7;
Style_t  styleKopeliovichComb    = 4;


Color_t	colorVitevBas0005				= kRed-6;
Color_t	colorVitevBas0510				= kOrange+1;
Color_t	colorVitevBas1020				= kOrange-5;
Color_t	colorVitevBas2040				= kGreen-6;
Color_t	colorVitevBas4060				= kCyan-6;
Color_t	colorVitevBas6080				= kBlue-6;

Style_t fillStyleVitev = 3766;
Style_t fillStyleEPOS = 0;//3454;
Style_t fillStyleEPOSRatio = 1001;//3454;
Style_t fillStyleXiao = 3002;
Style_t fillStyleWHDG = 3545;


Color_t	colorComb0005Box				= kRed-6;
Color_t	colorComb0510Box				= 807-6;
Color_t	colorComb1020Box				= 800-6;
Color_t	colorComb2040Box				= kGreen-6;
Color_t	colorComb4060Box				= kCyan-6;
Color_t	colorComb6080Box				= kBlue-6;

Style_t 	markerStyleConv		 	= 20 ;
Style_t 	markerStylePHOS		 	= 21 ;
Color_t 	colorConv 				= kBlack;
Color_t 	colorConvMC 				= kGray+1;
Color_t 	colorPHOS					= kRed+1;
Color_t 	colorPHOSMC				= kRed-7;
Color_t 	colorConvSyst				= kGray;
Color_t 	colorPHOSSyst				= kBlue-8;
Style_t 	fillStyleConv			 	= 0 ;
Style_t 	fillStylePHOS			 	= 0;
Style_t 	fillStyleEMCAL			 	= 0;
Size_t 	markerSizeInvYield			= 1.5;
Size_t 	markerSizeMass				= 3.2;
Style_t 	markerStyleConvMC			= 24 ;
Color_t 	colorPHOSMass					= kRed+1;
Color_t 	colorPHOSMCMass				= kRed-7;
Style_t 	markerStylePHOSMC		 	= 25 ;
Style_t  markerStyleEMCAL        = 33 ;
Style_t  markerStyleEMCALMC         = 27 ;
Color_t  colorEMCALMass              = kGreen+2;
Color_t  colorEMCALMCMass            = kGreen-6;


//********************* Variables Combmon Spectrum ************************************************************************
Style_t 	markerStyleCommmonSpectrumpp 	= 20 ;
Style_t 	markerStyleCommmonSpectrum0005 	= 20 ;
Style_t 	markerStyleCommmonSpectrum0010 	= 20 ;
Style_t 	markerStyleCommmonSpectrum0510 	= 21 ;
Style_t 	markerStyleCommmonSpectrum1020 	= 29 ;
Style_t 	markerStyleCommmonSpectrum2040 	= 33 ;
Style_t 	markerStyleCommmonSpectrum4060 	= 20 ;
Style_t 	markerStyleCommmonSpectrum6080 	= 21 ;

Style_t 	markerStylePHENIX200GeV 	= 25 ;
Style_t 	markerStylePHENIX62GeV 	= 27 ;
Style_t 	markerStylePHENIX39GeV 	= 24 ;
Style_t 	markerStyleWA98 	= 28 ;
Size_t 	markerSizePHENIX200GeV 	= 1.95;
Size_t 	markerSizePHENIX62GeV 	= 3;
Size_t 	markerSizePHENIX39GeV 	= 1.95;
Size_t 	markerSizeWA98 	= 1.95;

Size_t 	markerSizeCommonSpectrum0005 	= 2.;
Size_t 	markerSizeCommonSpectrum0010 	= 2.;
Size_t 	markerSizeCommonSpectrum0510 	= 2.;
Size_t 	markerSizeCommonSpectrum1020 	= 2.5;
Size_t 	markerSizeCommonSpectrum2040 	= 2.5;
Size_t 	markerSizeCommonSpectrum4060 	= 2.;
Size_t 	markerSizeCommonSpectrum6080 	= 2.;
Size_t 	markerSizeCommonSpectrumPi07TeV 	= 1.8;
Size_t 	markerSizeCommonSpectrumPi0900GeV = 1.8;
Size_t 	markerSizeCommonSpectrumEta7TeV 	= 2.2;
Size_t 	markerSizeSpectrum 			= 2.;
Size_t 	markerSizeChargedHadronSpectrum = 1.5;
Style_t 	styleFitCommonSpectrum 		= 1;
Width_t 	widthLinesBoxes;
Width_t 	widthStatErrBars;
Width_t 	widthCommonFit;
Width_t 	widthCommonSpectrumBoxes;
Width_t 	widthCommonErrors;
Color_t 	colorFitIndivid			= kRed+2;
Style_t	styleFitIndivid			= 1;

Bool_t 		pictDrawingOptions[4] = 			{kFALSE, kFALSE, kFALSE, kTRUE};
Bool_t 		conference;
TString 		date;
TString 		collisionSystemPP;		
TString 		collisionSystemCent0005;		
TString 		collisionSystemCent0510;		
TString 		collisionSystemCent1020;		
TString     collisionSystemPbPb;      
TString 		collisionSystemCent0020;		
TString 		collisionSystemCent2040;		
TString 		collisionSystemCent4060;		
TString 		collisionSystemCent6080;
TString 		forOutput;

TString     fileNameCaloEMCALPbPb;
TString 		fileNameCaloPHOSPbPb;
TString 		fileNameCaloPHOSPP;
TString 		fileNamePCMPP2760GeV;
TString		fileNamePreliminaryPbPb;
TString 		fileNameTheoryInput;
TString 		fileNameDataOtherEnergyInput;
TString 		nameFinalResDat;
TString 		prefix2;
Double_t 	minPtForFits;

Double_t 	mesonMassExpectPi0;

TFile*		filePCMPbPb;

// all given in %
Double_t commonCentralityErr0005 = 0.2;
Double_t commonCentralityErr0510 = 0.3;
Double_t commonCentralityErr0010 = 0.25;
Double_t commonCentralityErr1020 = 0.7;
Double_t commonCentralityErr2040 = 1.5;
Double_t commonCentralityErr4060 = 3;
Double_t commonCentralityErr6080 = 6.1;

TDirectory* directoryPi0PbPb0005;
TH1D*	histoNumberOfEvents0005;
TH1D*	histoYieldPbPb0005;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb0005;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrRAAPbPb0005;
TH1D*	histoPCMMassData0005;
TH1D*	histoPCMMassMC0005;
TH1D*	histoPCMWidthData0005;
TH1D*	histoPCMWidthMC0005;
Int_t	nEvent0005;
TGraphAsymmErrors* graphPCMYieldPi0SysErrPbPb0005Red;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAAPbPb0005Red;
TH1D*	histoYieldPbPb0005Red;
Int_t	nPointsPi00005;

TDirectory* directoryPi0PbPb0010;
TH1D*	histoNumberOfEvents0010;
TH1D*	histoYieldPbPb0010;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb0010;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrRAAPbPb0010;
TH1D* histoPCMMassData0010;
TH1D* histoPCMMassMC0010;
TH1D* histoPCMWidthData0010;
TH1D* histoPCMWidthMC0010;
Int_t	nEvent0010;
TGraphAsymmErrors* graphPCMYieldPi0SysErrPbPb0010Red;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAAPbPb0010Red;
TH1D*	histoYieldPbPb0010Red;
Int_t	nPointsPi00010;


TDirectory* directoryPi0PbPb0510;
TH1D*	histoNumberOfEvents0510;
TH1D*	histoYieldPbPb0510;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb0510;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrRAAPbPb0510;
TH1D*	histoPCMMassData0510;
TH1D*	histoPCMMassMC0510;
TH1D*	histoPCMWidthData0510;
TH1D*	histoPCMWidthMC0510;
Int_t	nEvent0510;
TGraphAsymmErrors* graphPCMYieldPi0SysErrPbPb0510Red;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAAPbPb0510Red;
TH1D*	histoYieldPbPb0510Red;
Int_t	nPointsPi00510;

TDirectory* directoryPi0PbPb1020;
TH1D*	histoNumberOfEvents1020;
TH1D*	histoYieldPbPb1020;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb1020;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrRAAPbPb1020;
TH1D*	histoPCMMassData1020;
TH1D*	histoPCMMassMC1020;
TH1D*	histoPCMWidthData1020;
TH1D*	histoPCMWidthMC1020;
Int_t	nEvent1020;
TGraphAsymmErrors* graphPCMYieldPi0SysErrPbPb1020Red;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAAPbPb1020Red;
TH1D*	histoYieldPbPb1020Red;
Int_t	nPointsPi01020;

TDirectory* directoryPi0PbPb2040;
TH1D*	histoNumberOfEvents2040;
TH1D*	histoYieldPbPb2040;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb2040;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrRAAPbPb2040;
TH1D*	histoPCMMassData2040;
TH1D*	histoPCMMassMC2040;
TH1D*	histoPCMWidthData2040;
TH1D*	histoPCMWidthMC2040;
Int_t	nEvent2040;
TGraphAsymmErrors* graphPCMYieldPi0SysErrPbPb2040Red;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAAPbPb2040Red;
TH1D*	histoYieldPbPb2040Red;
Int_t	nPointsPi02040;

TDirectory* directoryPi0PbPb4060;
TH1D*	histoNumberOfEvents4060;
TH1D*	histoYieldPbPb4060;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb4060;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrRAAPbPb4060;
TH1D*	histoPCMMassData4060;
TH1D*	histoPCMMassMC4060;
TH1D*	histoPCMWidthData4060;
TH1D*	histoPCMWidthMC4060;
Int_t	nEvent4060;
TGraphAsymmErrors* graphPCMYieldPi0SysErrPbPb4060Red;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAAPbPb4060Red;
TH1D*	histoYieldPbPb4060Red;
Int_t	nPointsPi04060;

TDirectory* directoryPi0PbPb6080;
TH1D*	histoNumberOfEvents6080;
TH1D*	histoYieldPbPb6080;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb6080;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrRAAPbPb6080;
TH1D*	histoPCMMassData6080;
TH1D*	histoPCMMassMC6080;
TH1D*	histoPCMWidthData6080;
TH1D*	histoPCMWidthMC6080;
Int_t	nEvent6080;
TGraphAsymmErrors* graphPCMYieldPi0SysErrPbPb6080Red;
TGraphAsymmErrors* graphPCMYieldPi0SysErrRAAPbPb6080Red;
TH1D*	histoYieldPbPb6080Red;
Int_t	nPointsPi06080;

TFile*	filePHOSPbPb;
TDirectory* directoryPHOSPi0PbPb0005;
TDirectory* directoryPHOSPi0PbPb0510;
TDirectory* directoryPHOSPi0PbPb0010;
TDirectory* directoryPHOSPi0PbPb1020;
TDirectory* directoryPHOSPi0PbPb2040;
TDirectory* directoryPHOSPi0PbPb4060;
TDirectory* directoryPHOSPi0PbPb6080;

TH1D* histoPi0PHOSPbPb0005;
TH1D*	histoPi0PHOSSysPbPb0005;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb0005;
TH1D*	histoPi0PHOSSysRAAPbPb0005;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb0005;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb0005Red;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb0005Red;
TH1D*	histoPi0PHOSPbPb0005Red;
TH1D* histoPHOSMassData0005;
TH1D* histoPHOSMassMC0005;
TH1D* histoPHOSWidthData0005;
TH1D* histoPHOSWidthMC0005;


TH1D* histoPi0PHOSPbPb0010;
TH1D*	histoPi0PHOSSysPbPb0010;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb0010;
TH1D*	histoPi0PHOSSysRAAPbPb0010;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb0010;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb0010Red;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb0010Red;
TH1D*	histoPi0PHOSPbPb0010Red;
TH1D* histoPHOSMassData0010;
TH1D* histoPHOSMassMC0010;
TH1D* histoPHOSWidthData0010;
TH1D* histoPHOSWidthMC0010;

TH1D* histoPi0PHOSPbPb0510;
TH1D*	histoPi0PHOSSysPbPb0510;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb0510;
TH1D*	histoPi0PHOSSysRAAPbPb0510;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb0510;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb0510Red;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb0510Red;
TH1D*	histoPi0PHOSPbPb0510Red;
TH1D* histoPHOSMassData0510;
TH1D* histoPHOSMassMC0510;
TH1D* histoPHOSWidthData0510;
TH1D* histoPHOSWidthMC0510;


TH1D* histoPi0PHOSPbPb1020;
TH1D*	histoPi0PHOSSysPbPb1020;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb1020;
TH1D*	histoPi0PHOSSysRAAPbPb1020;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb1020;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb1020Red;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb1020Red;
TH1D*	histoPi0PHOSPbPb1020Red;
TH1D* histoPHOSMassData1020;
TH1D* histoPHOSMassMC1020;
TH1D* histoPHOSWidthData1020;
TH1D* histoPHOSWidthMC1020;

TH1D* histoPi0PHOSPbPb2040;
TH1D*	histoPi0PHOSSysPbPb2040;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb2040;
TH1D*	histoPi0PHOSSysRAAPbPb2040;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb2040;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb2040Red;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb2040Red;
TH1D*	histoPi0PHOSPbPb2040Red;
TH1D*	histoPHOSMassData2040;
TH1D*	histoPHOSMassMC2040;
TH1D*	histoPHOSWidthData2040;
TH1D*	histoPHOSWidthMC2040;

TH1D* histoPi0PHOSPbPb4060;
TH1D*	histoPi0PHOSSysPbPb4060;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb4060;
TH1D*	histoPi0PHOSSysRAAPbPb4060;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb4060;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb4060Red;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb4060Red;
TH1D*	histoPi0PHOSPbPb4060Red;
TH1D*	histoPHOSMassData4060;
TH1D*	histoPHOSMassMC4060;
TH1D*	histoPHOSWidthData4060;
TH1D*	histoPHOSWidthMC4060;

TH1D* histoPi0PHOSPbPb6080;
TH1D*	histoPi0PHOSSysPbPb6080;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb6080;
TH1D*	histoPi0PHOSSysRAAPbPb6080;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb6080;
TGraphAsymmErrors* graphPHOSYieldPi0SysErrPbPb6080Red;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb6080Red;
TH1D*	histoPi0PHOSPbPb6080Red;
TH1D*	histoPHOSMassData6080;
TH1D*	histoPHOSMassMC6080;
TH1D*	histoPHOSWidthData6080;
TH1D*	histoPHOSWidthMC6080;

Double_t maxPtPi0PbPb0005;
Double_t maxPtPi0PbPb0010;
Double_t maxPtPi0PbPb0510;
Double_t maxPtPi0PbPb1020;
Double_t maxPtPi0PbPb2040;
Double_t maxPtPi0PbPb4060;
Double_t maxPtPi0PbPb6080;

Double_t xPtLimitsPbPbPHOS[21] =  {	0.4,0.6,0.8,1.0,1.2,
												1.4,1.6,1.8,2.0,2.2,
												2.4,2.6,3.0,3.5,4.0,
												5.0,6.0,8.0,10.,12.,
												15.};
Double_t xPtLimitsPbPbPHOS2[20] =  {	0.6,0.8,1.0,1.2,
												1.4,1.6,1.8,2.0,2.2,
												2.4,2.6,3.0,3.5,4.0,
												5.0,6.0,8.0,10.,12.,
												15.};

Double_t xPtLimitsPbPbPHOSCent[20] =  {	0.6,0.8,1.0,1.2,
												1.4,1.6,1.8,2.0,2.2,
												2.4,2.6,3.0,3.5,4.0,
												5.0,6.0,8.0,10.,12.,
												15.};

Double_t	fBinsPi0HIPtPCM[22] = 	{0.0,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,3.0,3.5,4.0,5.0,6.0,8.0,10.0,12.0,14.0};
Double_t	fBinsPi0HIPtPHOS[19] = 	{0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,3.0,3.5,4.0,5.0,6.0,8.0,10.0,12.0,15.0};

TGraphAsymmErrors*	graphYieldPi0CombPbPb0005;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0005StatErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0005SysErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0005Unshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0005StatErrUnshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0005SysErrUnshifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb0005;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb0005;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb0005Unshifted;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb0005Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb0005Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0SysErrPbPb0005Unshifted;
TF1*	fitModTsallisPi0PbPb0005;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0005Unscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0005StatErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0005SysErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0005SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0005StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb0005SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb0005SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb0005StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0005SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0005SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0005StatErrYShifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb0005Unscaled;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb0005;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb0005Sys;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb0005;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb0005Sys;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb0005;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb0005Stat;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb0005Sys;
TF1* fitYieldPi0PbPb0005;
TGraphAsymmErrors* 	graphYieldPi0CombPbPb0005YShifted;

TGraphAsymmErrors*	graphYieldPi0CombPbPb0010;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0010StatErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0010SysErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0010Unshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0010StatErrUnshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0010SysErrUnshifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb0010;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb0010;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb0010Unshifted;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb0010Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb0010Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0SysErrPbPb0010Unshifted;
TF1*	fitModTsallisPi0PbPb0010;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0010Unscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0010StatErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0010SysErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0010SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0010StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb0010SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb0010SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb0010StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0010SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0010SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0010StatErrYShifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb0010Unscaled;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb0010;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb0010Sys;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb0010;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb0010Sys;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb0010;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb0010Stat;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb0010Sys;
TF1* fitYieldPi0PbPb0010;
TGraphAsymmErrors* 	graphYieldPi0CombPbPb0010YShifted;



TGraphAsymmErrors*	graphYieldPi0CombPbPb0510;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0510StatErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0510SysErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0510Unshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0510StatErrUnshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0510SysErrUnshifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb0510;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb0510;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb0510Unshifted;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb0510Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb0510Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0SysErrPbPb0510Unshifted;
TF1*	fitModTsallisPi0PbPb0510;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0510Unscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0510StatErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0510SysErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0510SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb0510StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb0510SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb0510SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb0510StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0510SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0510SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0510StatErrYShifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb0510Unscaled;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb0510;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb0510Sys;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb0510;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb0510Sys;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb0510;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb0510Stat;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb0510Sys;
TF1* fitYieldPi0PbPb0510;
TGraphAsymmErrors* 	graphYieldPi0CombPbPb0510YShifted;

TGraphAsymmErrors*	graphYieldPi0CombPbPb1020;
TGraphAsymmErrors*	graphYieldPi0CombPbPb1020StatErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb1020SysErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb1020Unshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb1020StatErrUnshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb1020SysErrUnshifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb1020;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb1020;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb1020Unshifted;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb1020Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb1020Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0SysErrPbPb1020Unshifted;
TF1*	fitModTsallisPi0PbPb1020;
TGraphAsymmErrors*	graphYieldPi0CombPbPb1020Unscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb1020StatErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb1020SysErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb1020SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb1020StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb1020SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb1020SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb1020StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb1020SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb1020SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb1020StatErrYShifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb1020Unscaled;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb1020;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb1020Sys;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb1020;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb1020Sys;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb1020;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb1020Stat;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb1020Sys;
TF1* fitYieldPi0PbPb1020;
TGraphAsymmErrors* 	graphYieldPi0CombPbPb1020YShifted;

TGraphAsymmErrors*	graphYieldPi0CombPbPb2040;
TGraphAsymmErrors*	graphYieldPi0CombPbPb2040StatErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb2040SysErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb2040Unshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb2040StatErrUnshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb2040SysErrUnshifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb2040;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb2040;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb2040Unshifted;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb2040Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb2040Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0SysErrPbPb2040Unshifted;
TF1*	fitModTsallisPi0PbPb2040;
TGraphAsymmErrors*	graphYieldPi0CombPbPb2040Unscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb2040StatErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb2040SysErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb2040SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb2040StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb2040SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb2040SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb2040StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb2040SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb2040SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb2040StatErrYShifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb2040Unscaled;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb2040;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb2040Sys;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb2040;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb2040Sys;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb2040;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb2040Stat;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb2040Sys;
TF1* fitYieldPi0PbPb2040;
TGraphAsymmErrors* 	graphYieldPi0CombPbPb2040YShifted;

TGraphAsymmErrors*	graphYieldPi0CombPbPb4060;
TGraphAsymmErrors*	graphYieldPi0CombPbPb4060StatErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb4060SysErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb4060Unshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb4060StatErrUnshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb4060SysErrUnshifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb4060;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb4060;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb4060Unshifted;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb4060Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb4060Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0SysErrPbPb4060Unshifted;
TF1*	fitModTsallisPi0PbPb4060;
TGraphAsymmErrors*	graphYieldPi0CombPbPb4060Unscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb4060StatErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb4060SysErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb4060SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb4060StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb4060SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb4060SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb4060StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb4060SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb4060SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb4060StatErrYShifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb4060Unscaled;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb4060;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb4060Sys;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb4060;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb4060Sys;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb4060;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb4060Stat;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb4060Sys;
TF1* fitYieldPi0PbPb4060;
TGraphAsymmErrors* 	graphYieldPi0CombPbPb4060YShifted;

TGraphAsymmErrors*	graphYieldPi0CombPbPb6080;
TGraphAsymmErrors*	graphYieldPi0CombPbPb6080StatErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb6080SysErr;
TGraphAsymmErrors*	graphYieldPi0CombPbPb6080Unshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb6080StatErrUnshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb6080SysErrUnshifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb6080;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb6080;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb6080Unshifted;
TGraphAsymmErrors*	graphPCMYieldPi0SysErrPbPb6080Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0PbPb6080Unshifted;
TGraphAsymmErrors*	graphPHOSYieldPi0SysErrPbPb6080Unshifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb6080Unscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb6080StatErrUnscaled;
TGraphAsymmErrors*	graphYieldPi0CombPbPb6080SysErrUnscaled;
TF1* fitModTsallisPi0PbPb6080;
TF1* fitModTsallisPi0PbPb6080YShift;
TGraphAsymmErrors*	graphYieldPi0CombPbPb6080SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0CombPbPb6080StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb6080SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb6080SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PCMPbPb6080StatErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb6080SysErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb6080SysRAAErrYShifted;
TGraphAsymmErrors*	graphYieldPi0PHOSPbPb6080StatErrYShifted;
TGraphAsymmErrors*	graphPCMYieldPi0PbPb6080Unscaled;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb6080;
TGraphAsymmErrors*	graphRatioCombPHOSPi0PbPb6080Sys;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb6080;
TGraphAsymmErrors*	graphRatioCombConvPi0PbPb6080Sys;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb6080;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb6080Stat;
TGraphAsymmErrors*	graphRatioCombCombFitPbPb6080Sys;
TF1* fitYieldPi0PbPb6080;
TGraphAsymmErrors* 	graphYieldPi0CombPbPb6080YShifted;


Double_t nColl0005;
Double_t nColl0010;
Double_t nColl0510;
Double_t nColl1020;
Double_t nColl2040;
Double_t nColl4060;
Double_t nColl6080;

// Double_t nCollErr0020;
// Double_t nCollErr2040;
// Double_t nCollErr4060;
// Double_t nCollErr6080;

Double_t factorToInel = 1/1.117; //1; //							// this factor is multiplied with Raa and comes from trigger inelastic effiency for pp
Double_t xSection2760GeVpp = 		55.416*1e-3;
Double_t xSection2760GeVErrpp = 	3.9;
Double_t xSection2769GeVppINEL = 62.8*1e9;
Double_t recalcBarn = 			1e12; //NLO in pbarn!!!!

TFile*	fileFinalResultsPP;
TGraphAsymmErrors* graphInvSectionCombPi02760GeV;
TGraphAsymmErrors* graphInvSectionCombSysPi02760GeV;
TGraphAsymmErrors* graphInvSectionCombSysPi02760GeVPlot;
TGraphAsymmErrors* graphInvSectionPCMSysPi02760GeVPlot;
TGraphAsymmErrors* graphInvSectionPHOSSysPi02760GeVPlot;

TGraphErrors* Xiao_Raa_0020;
TGraphErrors* Xiao_Raa_0005;
TGraphErrors* Xiao_Raa_0510;
TGraphErrors* Xiao_Raa_1020;
TGraphErrors* 	Xiao_Raa_2040;
TGraphErrors* Xiao_Raa_4060;
TGraphErrors* Xiao_Raa_6080;
TGraphAsymmErrors*	gWHDG_Raa_0020;
TGraphAsymmErrors*	gWHDG_Raa_0005;
TGraphAsymmErrors*	gWHDG_Raa_0510;
TGraphAsymmErrors*	gWHDG_Raa_1020;
TGraphAsymmErrors*	gWHDG_Raa_2040;
TGraphAsymmErrors*	gWHDG_Raa_4060;
TGraphAsymmErrors*	gWHDG_Raa_6080;
TGraph* gEPOS_Spec_0005;
TGraph* gEPOS_Spec_0510;
TGraph* gEPOS_Spec_1020;
TGraph* gEPOS_Spec_0020;
TGraph* gEPOS_Spec_2040;
TGraph* gEPOS_Spec_4060;
TGraph* gEPOS_Spec_6080;
TGraph* gKopeliovichELoss_Spec_0005;
TGraph* gKopeliovichELoss_Spec_0510;
TGraph* gKopeliovichELoss_Spec_1020;
TGraph* gKopeliovichELoss_Spec_0020;
TGraph* gKopeliovichELoss_Spec_2040;
TGraph* gKopeliovichELoss_Spec_4060;
TGraph* gKopeliovichELoss_Spec_6080;
TGraph* gKopeliovichTotal_Spec_0005;
TGraph* gKopeliovichTotal_Spec_0510;
TGraph* gKopeliovichTotal_Spec_1020;
TGraph* gKopeliovichTotal_Spec_0020;
TGraph* gKopeliovichTotal_Spec_2040;
TGraph* gKopeliovichTotal_Spec_4060;
TGraph* gKopeliovichTotal_Spec_6080;

TGraph* gKopeliovichHydro_Spec_0005;
TGraph* gKopeliovichHydro_Spec_0510;
TGraph* gKopeliovichHydro_Spec_1020;
TGraph* gKopeliovichHydro_Spec_0020;
TGraph* gKopeliovichHydro_Spec_2040;
TGraph* gKopeliovichHydro_Spec_4060;
TGraph* gKopeliovichHydro_Spec_6080;

TGraph* gEPOS_RAA_0005;
TGraph* gEPOS_RAA_0510;
TGraph* gEPOS_RAA_1020;
TGraph* gEPOS_RAA_0020;
TGraph* gEPOS_RAA_2040;
TGraph* gEPOS_RAA_4060;
TGraph* gEPOS_RAA_6080;
TGraph* gRatioEPOSToFit0005;
TGraph* gRatioEPOSToFit0510;
TGraph* gRatioEPOSToFit1020;
TGraph* gRatioEPOSToFit2040;
TGraph* gRatioEPOSToFit4060;
TGraph* gRatioEPOSToFit6080;
TGraph* gRatioKopeliovichELossToFit0005;
TGraph* gRatioKopeliovichELossToFit0510;
TGraph* gRatioKopeliovichELossToFit1020;
TGraph* gRatioKopeliovichELossToFit2040;
TGraph* gRatioKopeliovichELossToFit4060;
TGraph* gRatioKopeliovichELossToFit6080;
TGraph* gRatioKopeliovichTotalToFit0005;
TGraph* gRatioKopeliovichTotalToFit0510;
TGraph* gRatioKopeliovichTotalToFit1020;
TGraph* gRatioKopeliovichTotalToFit2040;
TGraph* gRatioKopeliovichTotalToFit4060;
TGraph* gRatioKopeliovichTotalToFit6080;
TGraph* gRatioKopeliovichHydroToFit0005;
TGraph* gRatioKopeliovichHydroToFit0510;
TGraph* gRatioKopeliovichHydroToFit1020;
TGraph* gRatioKopeliovichHydroToFit2040;
TGraph* gRatioKopeliovichHydroToFit4060;
TGraph* gRatioKopeliovichHydroToFit6080;


TGraphErrors*	Vitev_Bas_Raa_0020;
TGraphErrors*	Vitev_Bas_Raa_0005;
TGraphErrors*	Vitev_Bas_Raa_0510;
TGraphErrors*	Vitev_Bas_Raa_1020;
TGraphErrors*  Vitev_Bas_Raa_2040;
TGraphErrors*  Vitev_Bas_Raa_4060;
TGraphErrors*  Vitev_Bas_Raa_6080;
TGraphErrors*	Vitev_ShlSel_Raa_0020;
TGraphErrors*	Vitev_ShlSel_Raa_0005;
TGraphErrors*	Vitev_ShlSel_Raa_0510;
TGraphErrors*	Vitev_ShlSel_Raa_1020;

TGraphErrors* graphPHENIX200GeVRAA_0010;
TGraphErrors* graphPHENIX200GeVRAA_1020;
TGraphErrors* graphPHENIX200GeVRAA_0020;
TGraphErrors* 	graphPHENIX200GeVRAA_2040;
TGraphErrors* graphPHENIX200GeVRAA_4060;
TGraphErrors* 	graphPHENIX200GeVRAA_6080;
TGraphErrors* 	graphPHENIX39GeVRAA_0010;
TGraphErrors* 	graphPHENIX39GeVRAA_1020;
TGraphErrors* 	graphPHENIX39GeVRAA_0020;
TGraphErrors* 	graphPHENIX39GeVRAA_2040;
TGraphErrors* 	graphPHENIX39GeVRAA_4060;
TGraphErrors* 	graphPHENIX39GeVRAA_6080;

TGraphErrors* 	graphPHENIX62GeVRAA_0010;
TGraphErrors* 	graphPHENIX62GeVRAA_1020;
TGraphErrors* 	graphPHENIX62GeVRAA_0020;
TGraphErrors* 	graphPHENIX62GeVRAA_2040;
TGraphErrors* 	graphPHENIX62GeVRAA_4060;
TGraphErrors* 	graphPHENIX62GeVRAA_6080;

TGraphErrors* 	graphWA98_17_3GeVRAA_0013;

TGraphAsymmErrors* 	graphInvSectionCombStatPi02760GeV;
TGraphAsymmErrors* 	graphInvSectionCombStatPi02760GeVPlot;
TGraphAsymmErrors* 	graphInvSectionPCMStatPi02760GeVPlot;
TGraphAsymmErrors* 	graphInvSectionPHOSStatPi02760GeVPlot;
TGraphAsymmErrors* 	graphInvSectionPHOSSysPi02760GeV;
TGraphAsymmErrors* 	graphInvSectionPHOSStatPi02760GeV;
TGraphAsymmErrors* 	graphInvSectionPCMSysPi02760GeV;
TGraphAsymmErrors* 	graphInvSectionPCMStatPi02760GeV;

TGraphAsymmErrors* 	graphInvSectionPCMStatPi02760GeVRed;
TGraphAsymmErrors* 	graphInvSectionPCMSysPi02760GeVRed;
TGraphAsymmErrors* 	graphInvSectionPHOSStatPi02760GeVRed;
TGraphAsymmErrors* 	graphInvSectionPHOSSysPi02760GeVRed;


TFile*	filePCMPP;
TDirectory* directoryPi0PP;
TH1D*	histoPCMMassDataPP;
TH1D*	histoPCMWidthDataPP;
TH1D*	histoPCMMassMCPP;
TH1D*	histoPCMWidthMCPP;

TFile*	filePHOSPP;
TDirectory* directoryPHOSPi0PP;
TH1D* 	histoPHOSMassDataPP;
TH1D* 	histoPHOSWidthDataPP;
TH1D* 	histoPHOSMassMCPP;
TH1D* 	histoPHOSWidthMCPP;

TFile*   fileEMCALPbPb;
TDirectory* directoryEMCALPi0PbPb0010;
   
TGraphErrors*  graphEMCALMassData0010;
TGraphErrors*  graphEMCALMassMC0010;
TGraphErrors*  graphEMCALWidthData0010;
TGraphErrors*  graphEMCALWidthMC0010;
TH1D*  histoEMCALMassData0010;
TH1D*  histoEMCALMassMC0010;
TH1D*  histoEMCALWidthData0010;
TH1D*  histoEMCALWidthMC0010;
