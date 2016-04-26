


Color_t 	colorPi02760GeV			= kMagenta+1;
Color_t		colorPi07TeV			= kBlue;
Color_t		colorPi0pPb5023GeV		= kGreen;


Style_t 	styleMarkerNLOMuHalf		= 24;
Style_t 	styleMarkerNLOMuOne		= 27;
Style_t 	styleMarkerNLOMuTwo		= 30;
Style_t         styleMarkerNLODSS14MuOne	= 28; //Tobe complete


Style_t 	styleLineNLOMuHalf		= 8;
Style_t 	styleLineNLODSS14MuHalf		= 8;
Style_t 	styleLineNLOMuOne		= 7;
Style_t         styleLineNLODSS14MuOne          = 7;
Style_t 	styleLineNLOMuTwo		= 4;
Style_t		styleLineNLODSS14MuTwo		= 4;

Style_t 	styleLineNLOMuTwoBKK		= 3;
Style_t 	styleLineNLOMuTwoDSS		= 6;


Size_t		sizeMarkerNLO			= 5.0;
Width_t		widthLineNLO			= 6.;
	 
Color_t 	colorNLOPi07TeVMuHalf 		= colorPi07TeV +2;
Color_t 	colorNLOPi07TeVMuOne		= colorPi07TeV;
Color_t 	colorNLOPi07TeVMuTwo		= colorPi07TeV -6 ;
Color_t 	colorNLOBKKPi07TeVMuTwo		= kCyan+3;
Color_t 	colorNLODSSPi07TeVMuTwo		= kBlack;
Color_t         colorNLODSS14Pi07TeVMuOne       = kMagenta + 3;//Tob complete
Color_t 	colorNLOPi02760GeVMuHalf	= colorPi02760GeV + 2;
Color_t         colorNLODSS14Pi02760GeVMuHalf	= colorPi02760GeV + 2;
Color_t 	colorNLOPi02760GeVMuOne		= colorPi02760GeV;
Color_t  	colorNLODSS14Pi02760GeVMuOne	= colorPi02760GeV;
Color_t 	colorNLOPi02760GeVMuTwo		= colorPi02760GeV - 6;
Color_t		colorNLODSS14Pi02760GeVMuTwo	= colorPi02760GeV - 6;
Color_t 	colorNLODSSPi02760GeVMuTwo	= kBlack;


//********************* Variables *******************



//********************* Variables Combmon Spectrum ************************************************************************
Style_t 	markerStyleCommmonSpectrumPi0pPb5023GeV =20;

Style_t         markerStyleTsallisFitRatioPi0pPb5023GeV  = 20;
Style_t         markerStyleBylinkinFitRatioPi0pPb5023GeV = 24;


Style_t 	markerStyleCommmonSpectrum 		= 20 ;
Style_t 	markerStyleCommmonSpectrum7TeV 		= 20 ;
Style_t 	markerStyleCommmonSpectrumPi07TeV 	= 20 ;
Style_t         markerStyleTsallisFitRatioPi07TeV 	= 20;
Style_t		markerStyleBylinkinFitRatioPi07TeV 	= 24;


Style_t 	markerStyleCommmonSpectrum2760GeV 	= 29;
Style_t		markerStyleTsallisFitRatioPi02760GeV	= 29;
Style_t		markerStyleBylinkinFitRatioPi02760GeV	= 30;


Style_t 	markerStyleCommmonSpectrum900GeV = 21 ;

Size_t 		markerSizeCommonSpectrum 	= 1.;
Size_t 		markerSizeCommonSpectrumPi07TeV 	= 1.0;
Color_t 	colorCommonSpectrumPi02760GeV = colorPi02760GeV+2;
Color_t 	colorCommonSpectrumPi02760GeVBox = colorPi02760GeV-8;
Color_t 	colorCommonSpectrumPi07TeV 	= colorPi07TeV+3;
Color_t 	colorCommonSpectrumPi07TeVBox = colorPi07TeV-8;
Color_t 	colorCommonSpectrumPi0pPb5023GeV  = colorPi0pPb5023GeV+3;
Color_t         colorCommonSpectrumPi0pPb5023GeVBox  = colorPi0pPb5023GeV-8;

Size_t 		markerSizeSpectrum 			= 1.;
Color_t 	colorConvPi07TeV			= colorPi07TeV+2;
Style_t 	styleFitCommonSpectrum 		= 1;
Style_t         styleBylinkinFitCommonSpectrum	= 2;
Width_t 	widthLinesBoxes;
Width_t 	widthStatErrBars;
Width_t 	widthCommonFit			= 4;
Width_t 	widthCommonSpectrumBoxes;
Width_t 	widthCommonErrors;
Color_t 	colorFitIndivid			= kRed+2;
Style_t		styleFitIndivid			= 1;


Double_t maxPtPi0 = 10.0;
Double_t minPtPi0 =  0.6;
TF1* 			fitInvCrossSectionPi02760GeV;







Double_t 		xSection;
// Double_t		xSection8TeVV0AND =		55.74*1e-3;		// from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
// Double_t		xSection8TeVErrUp =		0.46;			// from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
// Double_t		xSection8TeVErrDown =	0.46;			// from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
// Double_t		xSection8TeVT0AND =		55.74*1e-3;		// from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
// Double_t		xSection8TeVT0ErrUp =	0.46;			// from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
// Double_t		xSection8TeVT0ErrDown =	0.46;			// from https://indico.cern.ch/event/276276/contribution/3/material/slides/0.pdf
// Double_t		xSection7TeV =			62.22*1e-3;
// Double_t		xSection7TeVV0AND =		54.31*1e-3;
// Double_t		xSection7TeVErrUp =		2.18;
// Double_t		xSection7TeVErrDown =	2.18;
// Double_t		xSection900GeV =		47.78*1e-3;
// Double_t		xSection900GeVV0AND =	40.06*1e-3;
// Double_t		xSection900GeVErrUp =	2.39;
// Double_t		xSection900GeVErrDown =	1.86;
// Double_t		xSection2760GeV       = 	55.416*1e-3;
// Double_t		xSection2760GeVV0AND = 	47.73*1e-3;
// Double_t		xSection2760GeVErr = 	3.9;
// Double_t 		recalcBarn = 			1e12; //NLO in pbarn!!!!


TGraphAsymmErrors* 	graphDalitzPi0InvYieldStat7TeV;
TGraphAsymmErrors* 	graphDalitzPi0InvYieldSys7TeV;
TGraphAsymmErrors* 	graphDalitzPi0InvYieldSysA7TeV;
TGraphAsymmErrors* 	graphDalitzPi0InvYieldCompl7TeV;

TGraphAsymmErrors* 	graphDalitzPi0InvYieldStat7TeVUnShifted;  
TGraphAsymmErrors* 	graphDalitzPi0InvYieldSys7TeVUnShifted; 
TGraphAsymmErrors* 	graphDalitzPi0InvYieldSysA7TeVUnShifted;  
TGraphAsymmErrors* 	graphDalitzPi0InvYieldCompl7TeVUnShifted; 

TGraphAsymmErrors*	graphDalitzPi0InvCrossSectionStat7TeV;  		
TGraphAsymmErrors*      graphDalitzPi0InvCrossSectionSys7TeV;   		
TGraphAsymmErrors*	graphDalitzPi0InvCrossSectionSysA7TeV;  		
TGraphAsymmErrors* 	graphDalitzPi0InvCrossSectionCompl7TeV; 

TGraphAsymmErrors* 	graphDalitzPi0InvCrossSectionStat7TeVUnShifted; 
TGraphAsymmErrors* 	graphDalitzPi0InvCrossSectionSys7TeVUnShifted;   
TGraphAsymmErrors*      graphDalitzPi0InvCrossSectionSysA7TeVUnShifted;
TGraphAsymmErrors* 	graphDalitzPi0InvCrossSectionCompl7TeVUnShifted; 





TGraphAsymmErrors* 	graphDalitzPi0InvYieldStat2760GeV;
TGraphAsymmErrors* 	graphDalitzPi0InvYieldSys2760GeV;
TGraphAsymmErrors* 	graphDalitzPi0InvYieldSysA2760GeV;
TGraphAsymmErrors* 	graphDalitzPi0InvYieldCompl2760GeV;

TGraphAsymmErrors* 	graphDalitzPi0InvYieldStat2760GeVUnShifted;  
TGraphAsymmErrors* 	graphDalitzPi0InvYieldSys2760GeVUnShifted; 
TGraphAsymmErrors* 	graphDalitzPi0InvYieldSysA2760GeVUnShifted;  
TGraphAsymmErrors* 	graphDalitzPi0InvYieldCompl2760GeVUnShifted; 

TGraphAsymmErrors*	graphDalitzPi0InvCrossSectionStat2760GeV;  		
TGraphAsymmErrors*      graphDalitzPi0InvCrossSectionSys2760GeV;   		
TGraphAsymmErrors*	graphDalitzPi0InvCrossSectionSysA2760GeV;  		
TGraphAsymmErrors* 	graphDalitzPi0InvCrossSectionCompl2760GeV; 

TGraphAsymmErrors* 	graphDalitzPi0InvCrossSectionStat2760GeVUnShifted; 
TGraphAsymmErrors* 	graphDalitzPi0InvCrossSectionSys2760GeVUnShifted;   
TGraphAsymmErrors*      graphDalitzPi0InvCrossSectionSysA2760GeVUnShifted;
TGraphAsymmErrors* 	graphDalitzPi0InvCrossSectionCompl2760GeVUnShifted; 





TGraph*			graphRatioNLODalitzPi02760GeVMuHalf;
TGraph*			graphRatioNLODalitzPi02760GeVMuOne;
TGraph*			graphRatioNLODalitzPi02760GeVMuTwo;

//TGraphAsymmErrors*      histoRatioPythia8ToFit;
TGraphAsymmErrors*	graphRatioNLODSS14DalitzPi02760GeVMuHalf;
TGraphAsymmErrors*	graphRatioNLODSS14DalitzPi02760GeVMuOne;
TGraphAsymmErrors*	graphRatioNLODSS14DalitzPi02760GeVMuTwo;

TGraphAsymmErrors*	graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuHalf;
TGraphAsymmErrors*	graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuOne;
TGraphAsymmErrors*	graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuTwo;



TGraph*			graphRatioNLODalitzPi07TeVMuHalf;
TGraph*			graphRatioNLODalitzPi07TeVMuOne;
TGraph*			graphRatioNLODalitzPi07TeVMuTwo;
TGraph*			graphRatioNLODSS14DalitzPi07TeVMuOne;

TGraphAsymmErrors*	graphRatioInvCrossSectionPi07TeVTsallisFit;
TGraphAsymmErrors* 	graphRatioInvCrossSectionPi02760GeVTsallisFit;

TGraphAsymmErrors*	graphRatioInvCrossSectionPi02760GeVTBylinkinFit;
TGraphAsymmErrors*      graphRatioInvCrossSectionPi07TeVBylinkinFit;



	













