/****************************************************************************************************************************
****** 	provided by Gamma Conversion Group, PWG4, 														*****
******		Ana Marin, marin@physi.uni-heidelberg.de													*****
******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
*****************************************************************************************************************************/
	Style_t 			markerStylePCM		 		= 20;
	Style_t 			markerStylePHOS		 		= 21;
	Style_t 			markerStyleEMCAL		 	= 33;
	Style_t 			markerStylePCMMC			= 24;
	Style_t 			markerStylePHOSMC		 	= 25;
	Style_t 			markerStyleEMCALMC		 	= 27;
	Style_t				markerStylePCMEMCAL			= 34;
	Style_t				markerStylePCMEMCALMC		= 28;
	Style_t				markerStyleDalitz			= 29;
	Style_t				markerStyleDalitzMC			= 30;
	Style_t				markerStylePCMPHOS			= 34;
	Style_t				markerStylePCMPHOSMC		= 28;
	
	
	Size_t 				markerSizeMass				= 2.2;
	Size_t 				markerSizeMassMC			= 1.5;
	Size_t 				markerSizeInvYield			= 1.5;
	Size_t	 			markerSizeComparison 		= 1.2;
	
	Color_t 			colorPCM 					= kBlack;
	Color_t 			colorPCMMC 					= kGray+1;
	Color_t				colorEMCAL					= kGreen+2;
	Color_t				colorEMCALMC				= kGreen-6;
	Color_t 			colorPHOS					= kRed+1;
	Color_t 			colorPHOSMC					= kRed-7;
	Color_t				colorPCMEMCAL				= kCyan+2;
	Color_t				colorPCMEMCALMC				= kCyan-6;
	Color_t				colorPCMPHOS				= 807;
	Color_t				colorPCMPHOSMC				= kOrange+5;
	Color_t 			colorDalitz 				= kBlue+1;
	Color_t 			colorDalitzMC 				= kBlue-6;
	
	Style_t 			styleFitCommonSpectrum 		= 1;
	Width_t 			widthLinesBoxes;
	Width_t 			widthStatErrBars;
	Width_t 			widthCommonFit;
	Width_t 			widthCommonSpectrumBoxes;
	Width_t 			widthCommonErrors;

	Bool_t 				pictDrawingOptions[4] 		= {kFALSE, kFALSE, kFALSE, kTRUE};
	TString 			date;
	TString 		    collisionSystempPb;
	TString 			forOutput;

	TString 			nameFinalResDat;
	TString 			prefix2;
	Double_t		 	minPtForFits;
	Double_t 			mesonMassExpectPi0;
	Double_t 			mesonMassExpectEta;
	
	Double_t			maxPtPi0pPb;
	Double_t 			nCollpPb;

	
	TFile* 				filePCMpPb;
	TDirectory* 		directoryPCMPi0pPb;
	TH1D* 				histoPCMNumberOfEventspPb;
	TH1D* 				histoPCMYieldPi0pPb;
	TH1D* 				histoPCMRAWYieldPi0pPb;
	TGraphAsymmErrors* 	graphPCMYieldPi0SysErrpPb;
	TGraphAsymmErrors* 	graphPCMYieldPi0SysErrRAApPb;
	TH1D* 				histoPCMMassPi0DatapPb;
	TH1D* 				histoPCMMassPi0MCpPb;
	TH1D*			 	histoPCMWidthPi0DatapPb;
	TH1D* 				histoPCMWidthPi0MCpPb;
	TH1D* 				histoRatioPCMMassPi0DiffDataMC;
	TH1D* 				histoRatioPCMWidthPi0DiffDataMC;
	Double_t 			nPCMEventpPb;
	TDirectory* 		directoryPCMEtapPb;
	TH1D* 				histoPCMYieldEtapPb;
	TH1D* 				histoPCMRAWYieldEtapPb;
	TGraphAsymmErrors* 	graphPCMYieldEtaSysErrpPb;
	TGraphAsymmErrors* 	graphPCMYieldEtaSysErrRAApPb;
	TH1D* 				histoPCMEtaPi0RatiopPb;
	TGraphAsymmErrors* 	graphPCMEtaPi0RatioSysErrpPb;
	TH1D* 				histoPCMMassEtaDatapPb;
	TH1D* 				histoPCMMassEtaMCpPb;
	TH1D* 				histoPCMWidthEtaDatapPb;
	TH1D* 				histoPCMWidthEtaMCpPb;
	TH1D* 				histoRatioPCMMassEtaDiffDataMC;

	
	TString 			fileNameCaloEMCALpPb;
	TFile* 				fileEMCALpPb;
	TDirectory* 		directoryPi0EMCALpPb;
	TH1D* 				histoEMCALNumberOfEventspPb;
	TH1D* 				histoEMCALYieldPi0pPb;
	TH1D*				histoEMCALYieldPi0pPbSys;
	TH1D* 				histoEMCALRawYieldPi0pPb;
	TH1D* 				histoEMCALMassPi0DatapPb;
	TH1D* 				histoEMCALMassPi0MCpPb;
	TH1D* 				histoEMCALWidthPi0DatapPb;
	TH1D* 				histoEMCALWidthPi0MCpPb;
	TH1D*				histoRatioEMCALMassPi0DiffDataMC;
	TH1D*				histoRatioEMCALWidthPi0DiffDataMC;
	Double_t 			nEMCALEventpPb;

	TString 			fileNamePCMEMCALpPb;
	TFile* 				filePCMEMCALpPb;
	TDirectory* 		directoryPCMEMCALPi0pPb;
	TH1D* 				histoPCMEMCALYieldPi0pPb;
	TH1D* 				histoPCMEMCALRAWYieldPi0pPb;
	TGraphAsymmErrors* 	graphPCMEMCALYieldPi0SysErrpPb;
	TGraphAsymmErrors* 	graphPCMEMCALYieldPi0SysErrRAApPb;
	TH1D* 				histoPCMEMCALMassPi0DatapPb;
	TH1D* 				histoPCMEMCALMassPi0MCpPb;
	TH1D*			 	histoPCMEMCALWidthPi0DatapPb;
	TH1D* 				histoPCMEMCALWidthPi0MCpPb;
	TH1D*				histoRatioPCMEMCALMassPi0DiffDataMC;
	TH1D*				histoRatioPCMEMCALWidthPi0DiffDataMC;
	Double_t 			nPCMEMCALEventpPb;

	TDirectory* 		directoryPCMEMCALEtapPb;
	TH1D* 				histoPCMEMCALYieldEtapPb;
	TH1D* 				histoPCMEMCALRAWYieldEtapPb;
	TGraphAsymmErrors* 	graphPCMEMCALYieldEtaSysErrpPb;
	TGraphAsymmErrors* 	graphPCMEMCALYieldEtaSysErrRAApPb;
	TH1D* 				histoPCMEMCALEtaPi0RatiopPb;
	TGraphAsymmErrors* 	graphPCMEMCALEtaPi0RatioSysErrpPb;
	TH1D* 				histoPCMEMCALMassEtaDatapPb;
	TH1D* 				histoPCMEMCALMassEtaMCpPb;
	TH1D* 				histoPCMEMCALWidthEtaDatapPb;
	TH1D* 				histoPCMEMCALWidthEtaMCpPb;
	TH1D*				histoRatioPCMEMCALMassEtaDiffDataMC;

	TString 			fileNamePCMPHOSpPb;
	TFile* 				filePCMPHOSpPb;
	TDirectory* 		directoryPCMPHOSPi0pPb;
	TH1D* 				histoPCMPHOSYieldPi0pPb;
	TH1D* 				histoPCMPHOSRAWYieldPi0pPb;
	TGraphAsymmErrors* 	graphPCMPHOSYieldPi0SysErrpPb;
	TGraphAsymmErrors* 	graphPCMPHOSYieldPi0SysErrRAApPb;
	TH1D* 				histoPCMPHOSMassPi0DatapPb;
	TH1D* 				histoPCMPHOSMassPi0MCpPb;
	TH1D*			 	histoPCMPHOSWidthPi0DatapPb;
	TH1D* 				histoPCMPHOSWidthPi0MCpPb;
	TH1D*				histoRatioPCMPHOSMassPi0DiffDataMC;
	TH1D*				histoRatioPCMPHOSWidthPi0DiffDataMC;
	Double_t 			nPCMPHOSEventpPb;

	TDirectory* 		directoryPCMPHOSEtapPb;
	TH1D* 				histoPCMPHOSYieldEtapPb;
	TH1D* 				histoPCMPHOSRAWYieldEtapPb;
	TGraphAsymmErrors* 	graphPCMPHOSYieldEtaSysErrpPb;
	TGraphAsymmErrors* 	graphPCMPHOSYieldEtaSysErrRAApPb;
	TH1D* 				histoPCMPHOSEtaPi0RatiopPb;
	TGraphAsymmErrors* 	graphPCMPHOSEtaPi0RatioSysErrpPb;
	TH1D* 				histoPCMPHOSMassEtaDatapPb;
	TH1D* 				histoPCMPHOSMassEtaMCpPb;
	TH1D* 				histoPCMPHOSWidthEtaDatapPb;
	TH1D* 				histoPCMPHOSWidthEtaMCpPb;
	TH1D*				histoRatioPCMPHOSMassEtaDiffDataMC;
	
	TString 			fileNameEMCALEMCALpPb;
	TFile* 				fileEMCALEMCALpPb;
	TDirectory* 		directoryEMCALEMCALPi0pPb;
	TH1D* 				histoEMCALEMCALYieldPi0pPb;
	TH1D* 				histoEMCALEMCALRAWYieldPi0pPb;
	TGraphAsymmErrors* 	graphEMCALEMCALYieldPi0SysErrpPb;
	TGraphAsymmErrors* 	graphEMCALEMCALYieldPi0SysErrRAApPb;
	TH1D* 				histoEMCALEMCALMassPi0DatapPb;
	TH1D* 				histoEMCALEMCALMassPi0MCpPb;
	TH1D*			 	histoEMCALEMCALWidthPi0DatapPb;
	TH1D* 				histoEMCALEMCALWidthPi0MCpPb;
	TH1D*				histoRatioEMCALEMCALMassPi0DiffDataMC;
	TH1D*				histoRatioEMCALEMCALWidthPi0DiffDataMC;
	Double_t 			nEMCALEMCALEventpPb;

	TDirectory* 		directoryEMCALEMCALEtapPb;
	TH1D* 				histoEMCALEMCALYieldEtapPb;
	TH1D* 				histoEMCALEMCALRAWYieldEtapPb;
	TGraphAsymmErrors* 	graphEMCALEMCALYieldEtaSysErrpPb;
	TGraphAsymmErrors* 	graphEMCALEMCALYieldEtaSysErrRAApPb;
	TH1D* 				histoEMCALEMCALEtaPi0RatiopPb;
	TGraphAsymmErrors* 	graphEMCALEMCALEtaPi0RatioSysErrpPb;
	TH1D* 				histoEMCALEMCALMassEtaDatapPb;
	TH1D* 				histoEMCALEMCALMassEtaMCpPb;
	TH1D* 				histoEMCALEMCALWidthEtaDatapPb;
	TH1D* 				histoEMCALEMCALWidthEtaMCpPb;
	TH1D*				histoRatioEMCALEMCALMassEtaDiffDataMC;

	TString 			fileNamePHOSPHOSpPb;
	TFile* 				filePHOSPHOSpPb;
	TDirectory* 		directoryPHOSPHOSPi0pPb;
	TH1D* 				histoPHOSPHOSYieldPi0pPb;
	TH1D* 				histoPHOSPHOSRAWYieldPi0pPb;
	TGraphAsymmErrors* 	graphPHOSPHOSYieldPi0SysErrpPb;
	TGraphAsymmErrors* 	graphPHOSPHOSYieldPi0SysErrRAApPb;
	TH1D* 				histoPHOSPHOSMassPi0DatapPb;
	TH1D* 				histoPHOSPHOSMassPi0MCpPb;
	TH1D*			 	histoPHOSPHOSWidthPi0DatapPb;
	TH1D* 				histoPHOSPHOSWidthPi0MCpPb;
	TH1D*				histoRatioPHOSPHOSMassPi0DiffDataMC;
	TH1D*				histoRatioPHOSPHOSWidthPi0DiffDataMC;
	Double_t 			nPHOSPHOSEventpPb;
	
	TFile* 				fileChargedPionspPb;
	TH1D* 				histoChargedPionSyspPb;
	TH1D* 				histoChargedPionStatpPb;

	TString 			fileNameDalitzpPb;
	TFile* 				fileDalitzpPb;
	TDirectory* 		directoryDalitzPi0pPb;
	TGraphAsymmErrors* 	graphDalitzYieldPi0pPb;
	TGraphAsymmErrors* 	graphDalitzYieldPi0SysErrpPb;
	TH1D* 				histoDalitzMassPi0DatapPb;
	TH1D* 				histoDalitzMassPi0MCpPb;
	TH1D*			 	histoDalitzWidthPi0DatapPb;
	TH1D* 				histoDalitzWidthPi0MCpPb;
	TH1D*				histoRatioDalitzMassPi0DiffDataMC;
	TH1D*				histoRatioDalitzWidthPi0DiffDataMC;
	
	TString				fileNameCaloPHOSpPb;
	TFile*   			filePHOSpPb;
	TDirectory*			directoryPHOSPi0pPb;
	TH1D* 				histoPHOSYieldPi0pPb;
	TH1D* 				histoPHOSYieldPi0pPbSys;
	TH1D* 				histoPHOSMassPi0DatapPb;
	TH1D* 				histoPHOSMassPi0MCpPb;
	TH1D*			 	histoPHOSWidthPi0DatapPb;
	TH1D* 				histoPHOSWidthPi0MCpPb;
	TH1D*				histoRatioPHOSMassPi0DiffDataMC;
	TH1D*				histoRatioPHOSWidthPi0DiffDataMC;