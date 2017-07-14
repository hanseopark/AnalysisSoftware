
TString cent[5] = {"0-5%","5-10%","0-10%","20-40%","20-50%"};
TString labelcent[5] = {"0#font[122]{-}5%","5#font[122]{-}10%","0#font[122]{-}10%","20#font[122]{-}40%","20#font[122]{-}50%"};    
TString collisionSystPbPb[5] = { "0#font[122]{-}5% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","5#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","0#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV", "20#font[122]{-}40% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","20#font[122]{-}50% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV"};


//___________________________________ Labels definition _____________________________________________

TString collisionSystem2760GeV      = "Pb#font[122]{-}Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
TString collisionSystemPbPb0010     = "0#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
TString collisionSystemPbPb2040     = "20#font[122]{-}40% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
TString collisionSystemPbPb2050     = "20#font[122]{-}50% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";

TString collisionSystemPP2760GeV    = "pp, #sqrt{#it{s}} = 2.76 TeV";
TString collisionSystemPP7TeV       = "pp, #sqrt{#it{s}} = 7 TeV";

TString cent0010 = "0#font[122]{-}10%";
TString cent2050 = "20#font[122]{-}50%";
TString cent2040 = "20#font[122]{-}40%";



//___________________________________ Colors and Markers definition _____________________________________________

Style_t     markerStyle0010     = 20;
Style_t     markerStyle2050     = 21;
Style_t     markerStylepp       = 29;

Size_t      markerSize0010      = 2;
Size_t      markerSize2050      = 2;
Size_t      markerSizepp      = 3;

//plot range for yields
Double_t minYaxisYields = 1e-9;
Double_t maxYaxisYields = 1e3;
Double_t minPtYields = 0.25;
Double_t maxPtYields = 70.;
Double_t FontSize = 0.035;
Double_t maxYEtatoPi0 = 1.05;
Double_t minPtRange = 0.4;
Double_t maxPtRange = 30.;
Double_t minfitPt;

Color_t colorPCM0010        = kRed+1;
Color_t colorEMCal0010      = kRed+2;

Color_t colorPCM2050        = kBlue-7;
Color_t colorEMCal2050      = kBlue+2;

Color_t colorCombo0010      = kRed+1;
Color_t colorCombo2050      = kAzure+1;
Color_t colorComboRAA2050   = kAzure+2;

Color_t colorCharged        = kGray+1;

//charged pions: kOrange+1, kOrange-2    kBlue-3, kBlue-10
//charged kaons: kPink+5, kPink+1        kTeal-6, kCyan-8

Color_t colorPhenix         = kBlack;

//eta raa in phenix central kOrange+1, semicentral kGreen+2 nad kGreen-6

Color_t colorCracow0010     = kOrange+6;
Color_t colorCracow2050     = kAzure-4;//kCyan-9;

Color_t colorEPOS0010       = kPink+5;
Color_t colorEPOS2050       = kBlue;

Color_t colorCracowRatio    = kSpring+9;
Color_t colorEPOSRatio      = kViolet+6;

Color_t colorNLO            = kGray;
Color_t colorNLO0010        = kOrange-4;
Color_t colorNLO2050        = kCyan+1;

Width_t     widthStatErrBars;
Width_t     widthCommonSpectrumBoxes;
Width_t     widthCommonErrors;
Width_t     widthLinesBoxes     = 1.4;
Width_t     widthCommonFit      = 2;
Style_t     styleFitIndivid         = 1;


Color_t     colorComb           = kMagenta+2;
Color_t     colorFitIndivid     = kRed+2;
Style_t     markerStyleComb     = 20;
Size_t      markerSizeComb      = 2;

Color_t     colorCombLowPt          = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kFALSE);
Color_t     colorCombHighPt         = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kTRUE);
Style_t     markerStyleCombLowPt    = 20;
Style_t     markerStyleCombHighPt   = 20;
Size_t      markerSizeComparison = 2;

TString     nameMeasGlobal[11]  = {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM - EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCM 2010"};
Color_t     colorDet[11];
Color_t     colorDetMC[11];
Style_t     markerStyleDet[11];
Style_t     markerStyleDetMC[11];
Size_t      markerSizeDet[11];
Size_t      markerSizeDetMC[11];

Style_t     styleMarkerNLOMuHalf    = 24;
Style_t     styleMarkerNLOMuOne     = 27;
Style_t     styleMarkerNLOMuTwo     = 30;
Style_t     styleLineNLOMuHalf      = 8;
Style_t     styleLineNLOMuOne       = 7;
Style_t     styleLineNLOMuTwo       = 4;
Style_t     styleLineNLOMuTwoBKK    = 3;
Style_t     styleLineNLOMuTwoDSS    = 6;
Size_t      sizeMarkerNLO           = 1;
Width_t     widthLineNLO            = 2.;

Color_t colorComb0005               = kRed+1;
Color_t colorComb0010               = kRed+1;
Color_t colorComb0510               = 807;
Color_t colorComb1020               = 800;
Color_t colorComb2040               = kGreen+2;
Color_t colorComb4060               = kCyan+2;
Color_t colorComb6080               = kBlue+1;
Color_t colorWHDG0005               = kRed-4;
Color_t colorWHDG0510               = 807+1;
Color_t colorWHDG1020               = kYellow-6;
Color_t colorWHDG2040               = kGreen-3;
Color_t colorWHDG4060               = kCyan-3;
Color_t colorWHDG6080               = kBlue-3;
Color_t colorXiao0005               = kRed+3;
Color_t colorXiao0510               = 807+2;
Color_t colorXiao1020               = 800+1;
Color_t colorXiao2040               = kGreen+3;
Color_t colorXiao4060               = kCyan+3;
Color_t colorXiao6080               = kBlue+3;
Color_t  colorEPOS               = kBlack;
Color_t  colorKopeliovichHydro   = kBlue+2;
Color_t  colorKopeliovichELoss   = kGreen+4;
Color_t  colorKopeliovichComb    = kCyan+2;
Color_t colorEPOS0005               = kRed+4;
Color_t colorEPOS0510               = 807-2;
Color_t colorEPOS1020               = 800-3;
Color_t colorEPOS2040               = kGreen+4;
Color_t colorEPOS4060               = kCyan+4;
Color_t colorEPOS6080               = kBlue+4;
Color_t  colorEPOSFill0005       = kRed-6;
Color_t  colorEPOSFill0510       = kOrange;
Color_t  colorEPOSFill1020       = kYellow-6;
Color_t  colorEPOSFill2040       = kGreen-5;
Color_t  colorEPOSFill4060       = kCyan-5;
Color_t  colorEPOSFill6080       = kBlue-6;

Style_t  styleKopeliovichHydro   = 8;
Style_t  styleKopeliovichELoss   = 7;
Style_t  styleKopeliovichComb    = 4;

Color_t colorVitevBas0005               = kRed-6;
Color_t colorVitevBas0510               = kOrange+1;
Color_t colorVitevBas1020               = kOrange-5;
Color_t colorVitevBas2040               = kGreen-6;
Color_t colorVitevBas4060               = kCyan-6;
Color_t colorVitevBas6080               = kBlue-6;

Style_t fillStyleVitev = 3766;
Style_t fillStyleEPOS = 0;//3454;
Style_t fillStyleEPOSRatio = 1001;//3454;
Style_t fillStyleXiao = 3002;
Style_t fillStyleWHDG = 3545;

Color_t colorComb0005Box                = kRed-6;
Color_t colorComb0510Box                = 807-6;
Color_t colorComb1020Box                = 800-6;
Color_t colorComb2040Box                = kGreen-6;
Color_t colorComb2050Box                = kBlue-8;
Color_t colorComb4060Box                = kCyan-6;
Color_t colorComb6080Box                = kBlue-6;

Style_t     markerStyleConv         = 20 ;
Style_t     markerStylePHOS         = 21 ;
Color_t     colorConv               = kBlack;
Color_t     colorConvMC                 = kGray+1;
Color_t     colorPHOS                   = kRed+1;
Color_t     colorPHOSMC             = kRed-7;
Color_t     colorConvSyst               = kGray;
Color_t     colorPHOSSyst               = kBlue-8;
Style_t     fillStyleConv               = 0 ;
Style_t     fillStylePHOS               = 0;
Style_t     fillStyleEMCAL              = 0;
Size_t  markerSizeInvYield          = 1.5;
Size_t  markerSizeMass              = 3.2;
Style_t     markerStyleConvMC           = 24 ;
Color_t     colorPHOSMass                   = kRed+1;
Color_t     colorPHOSMCMass             = kRed-7;
Style_t     markerStylePHOSMC           = 25 ;
Style_t  markerStyleEMCAL        = 33 ;
Style_t  markerStyleEMCALMC         = 27 ;
Color_t  colorEMCALMass              = kGreen+2;
Color_t  colorEMCALMCMass            = kGreen-6;



//___________________________________ Bins/Offsets definition _____________________________________________


// Definition of final pt binning (has to be set manually)
Double_t xPtLimitsPi0[24] =  {0.0, 0.4, 0.6, 0.8, 1.0,
                              1.2, 1.4, 1.6, 1.8, 2.0,
                              2.2, 2.4, 2.6, 3.0, 3.5,
                              4.0, 5.0, 6.0, 8.0, 10.0,
                              12.0, 14.,17.,20.};
Double_t xPtLimitsEta[14] =  {0.0, 0.5, 1.0, 1.5, 2,
                              3.0, 4.0, 6.0, 8.0, 10.,
                              12., 14., 17., 20.};

// matrix order:    = {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM - EMCal",
//                    "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged"};
// Definition of offsets for stat & sys see output of function in shell, make sure pt bins match
Int_t offSetsPi0[11]    =   { 0, 2, 15, 0, 0,
                              0, 0, 0, 0, 0, 0};
Int_t offSetsPi0Sys[11]=    { 1, 4, 15, 0, 0,
                              0, 0, 0, 0, 0, 1};
Int_t offSetsPi0RAA[11]    =   { 4, 4, 15, 0, 0,
                              0, 0, 0, 0, 0, 0};
Int_t offSetsPi0RAASys[11]=    { 4, 4, 15, 0, 0,
                              0, 0, 0, 0, 0, 1};

Int_t offSetsEta[11]    =   { 0, 0, 6, 0, 0,
                              0, 0, 0, 0, 0, 0};//qui
Int_t offSetsEtaSys[11]=    { 2, 0, 6, 0, 0,
                              0, 0, 0, 0, 0, 0};
Int_t offSetsEtaRAA[11]    =   { 2, 0, 6, 0, 0,
                              0, 0, 0, 0, 0, 0};//qui
Int_t offSetsEtaRAASys[11]=    { 2, 0, 6, 0, 0,
                              0, 0, 0, 0, 0, 0};
//  Int_t offSetsEta[11]    =   { 0, 0, 10, 0, 0,
//                                0, 0, 0, 0, 0, 0};//qui
//  Int_t offSetsEtaSys[11]=    { 1, 0, 10, 0, 0,
//                                0, 0, 0, 0, 0, 0};
Bool_t      pictDrawingOptions[4] =             {kFALSE, kFALSE, kFALSE, kTRUE};


Style_t     markerStyleCommmonSpectrumpp    = 20;
Style_t     markerStyleCommmonSpectrum0005  = 20;
Style_t     markerStyleCommmonSpectrum0010  = 20;
Style_t     markerStyleCommmonSpectrum0510  = 21;
Style_t     markerStyleCommmonSpectrum1020  = 29;
Style_t     markerStyleCommmonSpectrum2040  = 33;
Style_t     markerStyleCommmonSpectrum4060  = 20;
Style_t     markerStyleCommmonSpectrum6080  = 21;
Style_t     markerStylePHENIX200GeV         = 25;
Style_t     markerStylePHENIX62GeV          = 27;
Style_t     markerStylePHENIX39GeV          = 24;
Style_t     markerStyleWA98                 = 28;
Style_t     styleFitCommonSpectrum          = 1;

Size_t  markerSizePHENIX200GeV  = 1.95;
Size_t  markerSizePHENIX62GeV   = 3;
Size_t  markerSizePHENIX39GeV   = 1.95;
Size_t  markerSizeWA98  = 1.95;
Size_t  markerSizeCommonSpectrum0005    = 2.;
Size_t  markerSizeCommonSpectrum0010    = 2.;
Size_t  markerSizeCommonSpectrum0510    = 2.;
Size_t  markerSizeCommonSpectrum1020    = 2.5;
Size_t  markerSizeCommonSpectrum2040    = 2.5;
Size_t  markerSizeCommonSpectrum4060    = 2.;
Size_t  markerSizeCommonSpectrum6080    = 2.;
Size_t  markerSizeCommonSpectrumPi07TeV = 1.8;
Size_t  markerSizeCommonSpectrumPi0900GeV = 1.8;
Size_t  markerSizeCommonSpectrumEta7TeV  = 2.2;
Size_t  markerSizeSpectrum               = 2.;
Size_t  markerSizeChargedHadronSpectrum  = 1.5;

// all given in %
Double_t commonCentralityErr0005 = 0.2;
Double_t commonCentralityErr0510 = 0.3;
Double_t commonCentralityErr0010 = 0.25;
Double_t commonCentralityErr1020 = 0.7;
Double_t commonCentralityErr2040 = 1.5;
Double_t commonCentralityErr4060 = 3;
Double_t commonCentralityErr6080 = 6.1;
Double_t nColl0010 = 1500;
Double_t nColl2040 = 438.4;
Double_t nColl2050 = 349.1;
Double_t nColl[5] = {1684.4, 1316, 1500,438.4,349.1};
Double_t nCollErr[5] = {190, 140, 165,42.,51.};


//Theory graph
TGraphErrors *graphEPOSpred_0010;
TGraphErrors *graphEPOSpred_2050;
TGraphAsymmErrors*  gWHDG_Pi0_Raa_0020;
TGraphAsymmErrors*  gWHDG_Pi0_Raa_0005;
TGraphAsymmErrors*  gWHDG_Pi0_Raa_0510;
TGraphAsymmErrors*  gWHDG_Pi0_Raa_0010;
TGraphAsymmErrors*  gWHDG_Pi0_Raa_1020;
TGraphAsymmErrors*  gWHDG_Pi0_Raa_2040;
TGraphAsymmErrors*  gWHDG_Pi0_Raa_2050;
TGraphAsymmErrors*  gWHDG_Pi0_Raa_4060;
TGraphAsymmErrors*  gWHDG_Pi0_Raa_6080;
TGraphAsymmErrors*  gWHDG_Eta_Raa_0010;
TGraphAsymmErrors*  gWHDG_Eta_Raa_2050;

//Charged particles graphs
TH1D*  histoChargedPionSpectraStat0010;
TH1D*  histoChargedPionSpectraSyst0010;
TH1D*  histoChargedKaonSpectraStat0010;
TH1D*  histoChargedKaonSpectraSyst0010;
TH1D*  histoChargedPionSpectraStat2040;
TH1D*  histoChargedPionSpectraSyst2040;
TH1D*  histoChargedKaonSpectraStat2040;
TH1D*  histoChargedKaonSpectraSyst2040;
TGraphAsymmErrors*  graphChargedRatioKaonToPion0010;
TGraphAsymmErrors*  graphChargedRatioKaonToPionSys0010;
TGraphAsymmErrors*  graphChargedRatioKaonToPion2040;
TGraphAsymmErrors*  graphChargedRatioKaonToPionSys2040;
TGraphAsymmErrors*  graphChargedPionRAA0010;
TGraphAsymmErrors*  graphChargedPionRAASys0010;
TGraphAsymmErrors*  graphChargedKaonRAA0010;
TGraphAsymmErrors*  graphChargedKaonRAASys0010;
TGraphAsymmErrors*  graphChargedPionRAA2040;
TGraphAsymmErrors*  graphChargedPionRAASys2040;
TGraphAsymmErrors*  graphChargedKaonRAA2040;
TGraphAsymmErrors*  graphChargedKaonRAASys2040;

//PHENIX graphs
TGraphErrors*   graphPHENIX200GeVPi0RAA_0010;
TGraphErrors*   graphPHENIX200GeVPi0RAA_1020;
TGraphErrors*   graphPHENIX200GeVPi0RAA_0020;
TGraphErrors*   graphPHENIX200GeVPi0RAA_2040;
TGraphErrors*   graphPHENIX200GeVPi0RAA_4060;
TGraphErrors*   graphPHENIX200GeVPi0RAA_6080;
TGraphErrors*   graphPHENIX39GeVPi0RAA_0010;
TGraphErrors*   graphPHENIX39GeVPi0RAA_1020;
TGraphErrors*   graphPHENIX39GeVPi0RAA_0020;
TGraphErrors*   graphPHENIX39GeVPi0RAA_2040;
TGraphErrors*   graphPHENIX39GeVPi0RAA_4060;
TGraphErrors*   graphPHENIX39GeVPi0RAA_6080;
TGraphErrors*   graphPHENIX62GeVPi0RAA_0010;
TGraphErrors*   graphPHENIX62GeVPi0RAA_1020;
TGraphErrors*   graphPHENIX62GeVPi0RAA_0020;
TGraphErrors*   graphPHENIX62GeVPi0RAA_2040;
TGraphErrors*   graphPHENIX62GeVPi0RAA_4060;
TGraphErrors*   graphPHENIX62GeVPi0RAA_6080;
TGraphErrors*   graphPHENIX200GeVEtaToPi0Ratio_0020;
TGraphErrors*   graphPHENIX200GeVEtaToPi0Ratio_2060;
TGraphErrors*   graphPHENIX200GeVEtaRAA_0020;
TGraphErrors*   graphPHENIX200GeVEtaRAA_2060;
TGraphErrors*   graphPHENIX200GeVEtaRAA_0010;
TGraphErrors*   graphPHENIX200GeVEtaRAA_2040;
TGraphErrors*   graphPHENIX200GeVEtaHighPtRAA_2060;
TGraphErrors*   graphWA98_17_3GeVPi0RAA_0013;

TH1D* histoChargedPionSpectraStat0005;
TH1D* histoChargedPionSpectraSyst0005;
TH1D* histoChargedPionSpectraStat0510;
TH1D* histoChargedPionSpectraSyst0510;

TH1D* histoChargedKaonSpectraStat0005;
TH1D* histoChargedKaonSpectraSyst0005;
TH1D* histoChargedKaonSpectraStat0510;
TH1D* histoChargedKaonSpectraSyst0510;

TGraphAsymmErrors* graphChargedRatioKaonToPion0005;
TGraphAsymmErrors* graphChargedRatioKaonToPionSys0005;
TGraphAsymmErrors* graphChargedRatioKaonToPion0510;
TGraphAsymmErrors* graphChargedRatioKaonToPionSys0510;

TGraphAsymmErrors* graphChargedPionRAA0005;
TGraphAsymmErrors* graphChargedPionRAASys0005;
TGraphAsymmErrors* graphChargedPionRAA0510;
TGraphAsymmErrors* graphChargedPionRAASys0510;

TGraphAsymmErrors* graphChargedKaonRAA0005;
TGraphAsymmErrors* graphChargedKaonRAASys0005;
TGraphAsymmErrors* graphChargedKaonRAA0510;
TGraphAsymmErrors* graphChargedKaonRAASys0510;

TGraphAsymmErrors* graphChargedPionRAA[5];
TGraphAsymmErrors* graphChargedPionRAASys[5];
TGraphAsymmErrors* graphChargedKaonRAA[5];
TGraphAsymmErrors* graphChargedKaonRAASys[5];
TGraphAsymmErrors* graphChargedKaonToPion[5];
TGraphAsymmErrors* graphChargedKaonToPionSys[5];


//*********************************************************************************************************//
//*****************************************  Neutral mesons   *********************************************//
//*********************************************************************************************************//

TGraphAsymmErrors* graphCombEtaToPi0Ratiopp7TeVNoXErrors;
TGraphAsymmErrors* graphCombEtaToPi0RatioSysErrpp7TeV;

TGraphAsymmErrors* graphInvSectionCombStatPi02760GeV;
TGraphAsymmErrors* graphInvSectionCombSysPi02760GeV;
TGraphAsymmErrors* graphInvSectionCombStatPi02760GeVPlot;
TGraphAsymmErrors* graphInvSectionCombSysPi02760GeVPlot;

TGraphAsymmErrors* graphInvSectionCombStatPi02760GeVforRAA;
TGraphAsymmErrors* graphInvSectionCombSysPi02760GeVforRAA;

TGraphAsymmErrors* graphInvSectionPCMStatPi02760GeV;
TGraphAsymmErrors* graphInvSectionPCMSysPi02760GeV;
TGraphAsymmErrors* graphInvSectionPCMStatPi02760GeVforRAA;
TGraphAsymmErrors* graphInvSectionPCMSysPi02760GeV_yShifted;
TGraphAsymmErrors* graphInvSectionPCMSysPi02760GeVforRAA;

TGraphAsymmErrors* graphInvSectionPHOSStatPi02760GeV;
TGraphAsymmErrors* graphInvSectionPHOSSysPi02760GeV;
TGraphAsymmErrors* graphInvSectionPHOSStatPi02760GeVforRAA;
TGraphAsymmErrors* graphInvSectionPHOSSysPi02760GeV_yShifted;
TGraphAsymmErrors* graphInvSectionPHOSSysPi02760GeVforRAA;

TGraphAsymmErrors* graphInvSectionEMCalStatPi02760GeV;
TGraphAsymmErrors* graphInvSectionEMCalSysPi02760GeV;
TGraphAsymmErrors* graphInvSectionEMCalStatPi02760GeVforRAA;
TGraphAsymmErrors* graphInvSectionEMCalSysPi02760GeV_yShifted;
TGraphAsymmErrors* graphInvSectionEMCalSysPi02760GeVforRAA;

TGraphAsymmErrors* graphInvSectionCombStatEta2760GeV;
TGraphAsymmErrors* graphInvSectionCombSysEta2760GeV;
TGraphAsymmErrors* graphInvSectionCombStatEta2760GeVPlot;
TGraphAsymmErrors* graphInvSectionCombSysEta2760GeVPlot;

TGraphAsymmErrors* graphInvSectionCombStatEta2760GeVforRAA;
TGraphAsymmErrors* graphInvSectionCombSysEta2760GeVforRAA;

TGraphAsymmErrors* graphInvSectionPCMStatEta2760GeV;
TGraphAsymmErrors* graphInvSectionPCMSysEta2760GeV;
TGraphAsymmErrors* graphInvSectionPCMStatEta2760GeVforRAA;
TGraphAsymmErrors* graphInvSectionPCMSysEta2760GeV_yShifted;
TGraphAsymmErrors* graphInvSectionPCMSysEta2760GeVforRAA;

TGraphAsymmErrors* graphInvSectionEMCalStatEta2760GeV;
TGraphAsymmErrors* graphInvSectionEMCalSysEta2760GeV;
TGraphAsymmErrors* graphInvSectionEMCalStatEta2760GeVforRAA;
TGraphAsymmErrors* graphInvSectionEMCalSysEta2760GeV_yShifted;
TGraphAsymmErrors* graphInvSectionEMCalSysEta2760GeVforRAA;

TGraphAsymmErrors* graphRatioEtaToPi0Comb2760GeVStatErr;
TGraphAsymmErrors* graphRatioEtaToPi0Comb2760GeVSysErr;
    
TGraphAsymmErrors* graphPCMPubPi0InvYieldStatPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMPubPi0InvYieldSysPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMPubPi0InvYieldStatPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMPubPi0InvYieldSysPbPb2760GeV_2040;

TH1D* histoPi0PHOSPbPb0010;
TGraphAsymmErrors* graphPHOSPi0InvYieldStatPbPb2760GeV_0010;
TGraphAsymmErrors* graphPHOSPi0InvYieldSysPbPb2760GeV_0010;
TGraphAsymmErrors* graphSysErrRAAYieldPi0PHOSPbPb0010;


TH1D* histoPCMPi0InvYieldPbPb2760GeV_0005;
TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_0005;
TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_0005;
TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0005;
TGraphAsymmErrors* graphPCMPi0InvYieldSysWOMat2760GeV_0005;
TH1D* histoPCMPi0RAAStatPbPb2760GeV_0005;
TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_0005;
TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_0005;

TH1D* histoPCMPi0InvYieldPbPb2760GeV_0510;
TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_0510;
TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_0510;
TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0510;
TGraphAsymmErrors* graphPCMPi0InvYieldSysWOMat2760GeV_0510;
TH1D* histoPCMPi0RAAStatPbPb2760GeV_0510;
TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_0510;
TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_0510;

TH1D* histoPCMPi0InvYieldPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010;
TGraphAsymmErrors* graphPCMPi0InvYieldSysWOMat2760GeV_0010;
TH1D* histoPCMPi0RAAStatPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_0010;

TH1D* histoPCMPi0InvYieldPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2040;
TGraphAsymmErrors* graphPCMPi0InvYieldSysWOMat2760GeV_2040;
TH1D* histoPCMPi0RAAStatPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_2040;

TH1D* histoPCMPi0InvYieldPbPb2760GeV_2050;
TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_2050;
TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_2050;
TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050;
TGraphAsymmErrors* graphPCMPi0InvYieldSysWOMat2760GeV_2050;
TH1D* histoPCMPi0RAAStatPbPb2760GeV_2050;
TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_2050;
TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_2050;

TH1D* histoPCMEtaInvYieldPbPb2760GeV_0005;
TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_0005;
TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_0005;
TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0005;
TGraphAsymmErrors* graphPCMEtaInvYieldSysWOMat2760GeV_0005;
TH1D* histoPCMEtaRAAStatPbPb2760GeV_0005;
TGraphAsymmErrors* graphPCMEtaRAAStatPbPb2760GeV_0005;
TGraphAsymmErrors* graphPCMEtaRAASysPbPb2760GeV_0005;
TH1D* histoPCMEtatoPi0Stat2760GeV_0005;
TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_0005;
TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_0005;

TH1D* histoPCMEtaInvYieldPbPb2760GeV_0510;
TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_0510;
TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_0510;
TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0510;
TGraphAsymmErrors* graphPCMEtaInvYieldSysWOMat2760GeV_0510;
TH1D* histoPCMEtaRAAStatPbPb2760GeV_0510;
TGraphAsymmErrors* graphPCMEtaRAAStatPbPb2760GeV_0510;
TGraphAsymmErrors* graphPCMEtaRAASysPbPb2760GeV_0510;
TH1D* histoPCMEtatoPi0Stat2760GeV_0510;
TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_0510;
TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_0510;

TH1D* histoPCMEtaInvYieldPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010;
TGraphAsymmErrors* graphPCMEtaInvYieldSysWOMat2760GeV_0010;
TH1D* histoPCMEtaRAAStatPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMEtaRAAStatPbPb2760GeV_0010;
TGraphAsymmErrors* graphPCMEtaRAASysPbPb2760GeV_0010;
TH1D* histoPCMEtatoPi0Stat2760GeV_0010;
TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_0010;
TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_0010;

TH1D* histoPCMEtaInvYieldPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2040;
TGraphAsymmErrors* graphPCMEtaInvYieldSysWOMat2760GeV_2040;
TH1D* histoPCMEtaRAAStatPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMEtaRAAStatPbPb2760GeV_2040;
TGraphAsymmErrors* graphPCMEtaRAASysPbPb2760GeV_2040;
TH1D* histoPCMEtatoPi0Stat2760GeV_2040;
TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_2040;
TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_2040;

TH1D* histoPCMEtaInvYieldPbPb2760GeV_2050;
TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_2050;
TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_2050;
TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050;
TGraphAsymmErrors* graphPCMEtaInvYieldSysWOMat2760GeV_2050;
TH1D* histoPCMEtaRAAStatPbPb2760GeV_2050;
TGraphAsymmErrors* graphPCMEtaRAAStatPbPb2760GeV_2050;
TGraphAsymmErrors* graphPCMEtaRAASysPbPb2760GeV_2050;
TH1D* histoPCMEtatoPi0Stat2760GeV_2050;
TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_2050;
TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_2050;

// histoEMCalPi0InvYieldStatPbPb2760GeV_0010;
// graphEMCalPi0InvYieldStatPbPb2760GeV_0010;
// graphEMCalPi0InvYieldSysPbPb2760GeV_0010;
// graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_0010;
// histoEMCalPi0RAAStatPbPb2760GeV_0010;
// graphEMCalPi0RAAStatPbPb2760GeV_0010;
// graphEMCalPi0RAASysPbPb2760GeV_0010;
// 
// histoEMCalPi0InvYieldStatPbPb2760GeV_2050;
// graphEMCalPi0InvYieldStatPbPb2760GeV_2050;
// graphEMCalPi0InvYieldSysPbPb2760GeV_2050;
// graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_2050;
// histoEMCalPi0RAAStatPbPb2760GeV_2050;
// graphEMCalPi0RAAStatPbPb2760GeV_2050;
// graphEMCalPi0RAASysPbPb2760GeV_2050;
// 
// histoEMCalEtaInvYieldStatPbPb2760GeV_0010;
// graphEMCalEtaInvYieldStatPbPb2760GeV_0010;
// graphEMCalEtaInvYieldSysPbPb2760GeV_0010;
// graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_0010;
// histoEMCalEtaRAAStatPbPb2760GeV_0010;
// graphEMCalEtaRAAStatPbPb2760GeV_0010;
// graphEMCalEtaRAASysPbPb2760GeV_0010;
// histoEMCalEtatoPi0Stat2760GeV_0010;
// graphEMCalEtatoPi0Stat2760GeV_0010;
// graphEMCalEtatoPi0Sys2760GeV_0010;
// 
// histoEMCalEtaInvYieldStatPbPb2760GeV_2050;
// graphEMCalEtaInvYieldStatPbPb2760GeV_2050;
// graphEMCalEtaInvYieldSysPbPb2760GeV_2050;
// graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_2050;
// histoEMCalEtaRAAStatPbPb2760GeV_2050;
// graphEMCalEtaRAAStatPbPb2760GeV_2050;
// graphEMCalEtaRAASysPbPb2760GeV_2050;
// histoEMCalEtatoPi0Stat2760GeV_2050;
// graphEMCalEtatoPi0Stat2760GeV_2050;
// graphEMCalEtatoPi0Sys2760GeV_2050;


//*********************************************************************************************************//
//***************************************     generalise arrays     ***************************************//
//*********************************************************************************************************//

TH1D *histoarrayPCMPi0InvYieldPbPb2760GeV[5];
TH1D* histoarrayPCMPi0InvYieldPbPb2760GeVYshifted[5];
TH1D *histoarrayPCMEtaInvYieldPbPb2760GeV[5];
TH1D* histoarrayPCMEtaInvYieldPbPb2760GeVYshifted[5];
TH1D* histoarrayPCMEtatoPi0Stat2760GeV[5];
TGraphAsymmErrors* grapharrayPCMPi0InvYieldStatPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMPi0InvYieldSysPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMPi0InvYieldSysWOMat2760GeV[5];
TGraphAsymmErrors* grapharrayPCMEtaInvYieldStatPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMEtaInvYieldSysPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMEtaInvYieldSysWOMat2760GeV[5];
TGraphAsymmErrors* grapharrayPCMEtatoPi0Stat2760GeV[5];
TGraphAsymmErrors* grapharrayPCMEtatoPi0Sys2760GeV[5];
TGraphAsymmErrors* grapharrayPCMInvYieldTotPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMInvYieldStatPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMInvYieldSysPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMInvYieldSysWOMat2760GeV[5];
TGraphAsymmErrors* grapharrayPCMPi0RAAStatPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMPi0RAASysPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMEtaRAAStatPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMEtaRAASysPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMRAAStatPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayPCMRAASysPbPb2760GeV[5];

//syst and stat relative errors
TH1D *statErrorCollectionLHC11h[5][11];
TH1D *statErrorCollectionEtatoPi0LHC11h[5][11];
TH1D *statErrorCollectionRaaLHC11h[5][11];
TGraphAsymmErrors *sysErrorCollectionLHC11h[5][11];
TGraphAsymmErrors *sysErrorCollectionEtatoPi0LHC11h[5][11];
TGraphAsymmErrors *sysErrorCollectionRaaLHC11h[5][11];


//**********************************************************************************************************************//
//************************************* Calculating bin shifted spectra & fitting **************************************//
//**********************************************************************************************************************//
//cloning spectra for shifting
TGraphAsymmErrors* grapharrayPCMInvYieldTotPbPb2760GeVUnshifted[5];
TGraphAsymmErrors* grapharrayPCMInvYieldStatPbPb2760GeVUnshifted[5];
TGraphAsymmErrors* grapharrayPCMInvYieldSysPbPb2760GeVUnshifted[5];
TGraphAsymmErrors* grapharrayPCMInvYieldSysWOMat2760GeVUnshifted[5];
TGraphAsymmErrors* grapharrayPCMInvYieldStatPbPb2760GeVYshifted[5];
TGraphAsymmErrors* grapharrayPCMInvYieldSysPbPb2760GeVYshifted[5];
TGraphAsymmErrors* grapharrayRatioToFitStatPbPb2760GeV[5];
TGraphAsymmErrors* grapharrayRatioToFitSysPbPb2760GeV[5];

