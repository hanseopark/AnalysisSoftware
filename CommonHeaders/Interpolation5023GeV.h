extern TRandom*         gRandom;	
extern TBenchmark*      gBenchmark;
extern TSystem*         gSystem;
extern TMinuit*         gMinuit;

void CalcRpPb(  TGraphAsymmErrors* PPSpectrumSystErr, TGraphAsymmErrors*  PPSpectrumStatErr, TGraphAsymmErrors* pPbSpectrumSystErr, TGraphAsymmErrors* pPbSpectrumStatErr,
                TGraphAsymmErrors** graphRpPbSystErr, TGraphAsymmErrors** graphRpPbStatErr, Bool_t Individual=kTRUE);

void CalcRpPbInvYield(  TGraphAsymmErrors* PPSpectrumSystErr, TGraphAsymmErrors*  PPSpectrumStatErr, TGraphAsymmErrors* pPbSpectrumSystErr, TGraphAsymmErrors* pPbSpectrumStatErr,
                        TGraphAsymmErrors** graphRpPbSystErr, TGraphAsymmErrors** graphRpPbStatErr);
void GetTGraphErrorsUpDown(TGraphAsymmErrors* graphYieldError, TGraphErrors** graphRelErrUpDown);
Double_t* ExtractRelErrUpAsymmGraphClone(TGraphAsymmErrors* graph);
Double_t* ExtractRelErrDownAsymmGraphClone(TGraphAsymmErrors* graph);
TGraphAsymmErrors* CalculateSystErrors(TGraphErrors* spectrum,TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);
TGraphAsymmErrors* CalculateSystErrors(TGraphAsymmErrors* spectrum,TGraphErrors* g1,TGraphErrors* g2);
//TF1* RebinWithFitToTGraph(TGraphAsymmErrors *spectrum, TGraphAsymmErrors** newSpectrum, TGraphAsymmErrors *newBins, TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters,Double_t probability);
TF1* RebinWithFitToTGraph(TGraphAsymmErrors *spectrum, TGraphAsymmErrors** newSpectrum, TGraphAsymmErrors *newBins,TF1* CurrentFit, TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters,Int_t fixParNumber);
TF1* FillTGraphEYWithFitErr(TGraphAsymmErrors *spectrum, TGraphAsymmErrors** newSpectrum, TGraphAsymmErrors *newBins, TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters, TString meson);
TF1* RebinWithFitToTGraphWithMeanErr(TGraphAsymmErrors *spectrum, TGraphAsymmErrors** newSpectrum, TGraphAsymmErrors *newBins,TF1* CurrentFit, TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters, Int_t fixParNumber, TString meson);
TGraphErrors *GetInterpolSpectrum2D(TGraphErrors *g1, TGraphErrors *g2,Double_t d1, Double_t d2,Double_t dSqrts);
TGraphAsymmErrors* GetChargeParticlesRpPb2013(TString typeErr);
TGraphAsymmErrors* GetChargeParticlesRpPb2012(TString typeErr);
TGraphAsymmErrors* BinYShiftwithErrhigh(TGraphAsymmErrors* spectrum,TGraphAsymmErrors* spectrumErrhigh);
TGraphAsymmErrors* BinYShiftwithErrlow(TGraphAsymmErrors* spectrum,TGraphAsymmErrors* spectrumErrlow);
TF1* RebinWithFitToTGraphWithUpDownYShifted(TGraphAsymmErrors *spectrumStatErr, TGraphAsymmErrors *spectrumSystErr,TGraphAsymmErrors *newBins,TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters,  TGraphAsymmErrors** newSpectrumStatErr,
					     TGraphAsymmErrors** newSpectrumSystErr,TGraphAsymmErrors** spectrumYShiftedDownStatErr,TGraphAsymmErrors** spectrumYShiftedUpStatErr,
					     TF1** FitToSpectrumYShiftedDown, TF1** FitToSpectrumYShiftedUp, TString meson="",TString energy="",TString system="");

void ExtrapolateSpectrum(TF1* FitToSpectrum, TGraphAsymmErrors *spectrumStatErr, TGraphAsymmErrors *spectrumSystErr, TGraphAsymmErrors** newSpectrumStatErr, TGraphAsymmErrors** newSpectrumSystErr, TString system, TString energy, TString meson);

TGraphErrors *ConvertTGraphAsymmErrorstoTGraphErrors(TGraphAsymmErrors* g1);
TGraphAsymmErrors* ConvertTGraphErrorstoTGraphAsymmErrors(TGraphErrors* g1);
TH1F* GetPPReferenceFromPythia(TString fileName);
void SavingFiles(TString outputDir,TString meson, TString System);

Double_t xSection2760GeVppINEL  = 62.8*1e-3;
Double_t xSectionpPb5023GeVINEL = 70*1e-3;
//Double_t	xSection7TeVINEL         =  73.2*1e-3;
//Double_t        recalcBarn               =  1e12; 
Double_t fNcoll                 = 6.9;
Double_t fTpPb                  = 0.0983e3*(1/recalcBarn);
Double_t fTpPbErr               = 0.0035e3*(1/recalcBarn);

TGraphAsymmErrors* graphAErrosInterPolation5023GeVBinStatErr;
TGraphAsymmErrors* graphAErrosInterPolation5023GeVBinSystErr;	
TGraphAsymmErrors* graphErrosInterPolation5023GeVBinSystWOMatErr = NULL;
TGraphAsymmErrors* graphInvYieldPi0pPb5023GeVXShiftedSystWOMatErr = NULL;
TGraphAsymmErrors* graphInvYieldPi0pPb5023GeVYShiftedSystWOMatErr = NULL;

TGraphAsymmErrors* graphRpPbStatErr;
TGraphAsymmErrors* graphRpPbSystErr;
TGraphAsymmErrors* graphRpPbPHOSStatErr;
TGraphAsymmErrors* graphRpPbPHOSSystErr;

TGraphErrors* graphAlpha[3];
TGraphErrors** graphPtvsSqrts[3]; 
TGraphErrors** gPtvsEnergiesSystem[3];
TF1** fPowerlawSystem[3];

TGraphAsymmErrors* graphAsymmChargedParticlesRpPbSystErr;
TGraphAsymmErrors* graphAsymmChargedParticlesRpPbStatErr;
TGraphAsymmErrors* graphChargedPionRpPbSystErr;
TGraphAsymmErrors* graphChargedPionRpPbStatErr;
TGraph* graphPi0ESP09sPi0KKP;
TGraph* graphPi0ESP09sPi0AKK;
TGraph* graphPi0DSS5000;
TGraphAsymmErrors* graphAsymmErrorsPi0DSS5000;
TGraph* graphPi0CGC;

/*Labels*/

TString invYieldLabel       = "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})";
TString pTLabel             = "#it{p}_{T} (GeV/#it{c})";
TString RpPbLabel           = "#it{R}_{pPb}";
TString energyLabel         = "p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV";


TString DalitzLabel         = "#pi^{0} #rightarrow e^{+}e^{-}#gamma";
TString PCMLabel            = "#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}";
TString CaloLabel           = "#pi^{0} #rightarrow #gamma #gamma";
TString PCMEMCalLabel       = "#pi^{0} #rightarrow #gamma_{conv}#gamma_{calo}";

TString DalitzLabelEta         = "#eta #rightarrow e^{+}e^{-}#gamma";
TString PCMLabelEta            = "#eta #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}";
TString CaloLabelEta           = "#eta #rightarrow #gamma #gamma";
TString PCMEMCalLabelEta       = "#eta #rightarrow #gamma_{conv}#gamma_{calo}";


TString ChargedPionsLabel   = "#pi^{+} + #pi^{-}";
TString thesisPlotLabel     = "";

enum{kDalitz,kPCM,kPHOS,kEMCal,kPCMEMCal,kPion,kKaon,kProton,kNAna};
Int_t colorsArray[kNAna];

TGraphAsymmErrors* GetChargeParticlesRpPb2012(TString typeErr){
  // Plot: p8424_d3x1y1
    double p8424_d3x1y1_xval[]      = { 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 
                                        0.975, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 
                                        1.95, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 
                                        3.9, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 
                                        10.5, 11.5, 12.5, 13.5, 15.0, 18.0 };
    double p8424_d3x1y1_xerrminus[] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025,
                                        0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050,
                                        0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100,
                                        0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5,
                                        0.5, 0.5, 0.5, 1.0, 2.0  };
    double p8424_d3x1y1_xerrplus[]  = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025,
                                        0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050,
                                        0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100,
                                        0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5,
                                        0.5, 0.5, 0.5, 1.0, 2.0  };
    double p8424_d3x1y1_yval[]      = { 0.5824, 0.5987, 0.6139, 0.6337, 0.6514, 0.6716, 0.6878, 0.6989, 0.7118, 
                                        0.722, 0.7485, 0.7822, 0.8046, 0.8306, 0.8539, 0.8811, 0.8947, 0.9145, 0.9431, 
                                        0.9584, 0.9844, 1.01, 1.034, 1.051, 1.051, 1.08, 1.083, 1.104, 1.102, 
                                        1.111, 1.1, 1.091, 1.111, 1.058, 1.035, 1.087, 1.064, 1.064, 1.11, 
                                        0.9183, 1.144, 1.061, 1.17, 0.9558, 1.121 };
    double p8424_d3x1y1_yerrminus[] = { 0.05600723167591842, 0.05670881765651617, 0.056910631695668255, 0.058712264476853564, 0.060314011639087645, 
                                        0.06221575363201831, 0.06371765846294103, 0.06472232999514155, 0.0659245781177248, 0.06692697512961422, 
                                        0.06931623186527092, 0.07251992829560713, 0.07452684080249208, 0.07693438762997987, 0.07914271918502674, 
                                        0.08165512843661445, 0.08296969326205805, 0.08478519918004557, 0.08740583504549339, 0.08892429364352579, 
                                        0.09128335006998813, 0.09413288479590966, 0.09618731725128839, 0.09725224933131367, 0.09732933781753578, 
                                        0.10031948963187563, 0.1004987562112089, 0.10259142264341595, 0.10270345661174213, 0.1039471019317037, 
                                        0.10270345661174213, 0.10196568050084304, 0.10932977636490435, 0.10413932974625868, 0.10235721762533408, 
                                        0.1078563859954523, 0.10412012293500235, 0.1077032961426901, 0.11764777940955791, 0.10736125930707034, 
                                        0.13917614738165443, 0.14509651959988565, 0.1776091213873882, 0.14606854555310667, 0.17462244987400674 };
    double p8424_d3x1y1_yerrplus[]  = { 0.05600723167591842, 0.05670881765651617, 0.056910631695668255, 0.058712264476853564, 0.060314011639087645,
                                        0.06221575363201831, 0.06371765846294103, 0.06472232999514155, 0.0659245781177248, 0.06692697512961422, 
                                        0.06931623186527092, 0.07251992829560713, 0.07452684080249208, 0.07693438762997987, 0.07914271918502674, 
                                        0.08165512843661445, 0.08296969326205805, 0.08478519918004557, 0.08740583504549339, 0.08892429364352579, 
                                        0.09128335006998813, 0.09413288479590966, 0.09618731725128839, 0.09725224933131367, 0.09732933781753578,
                                        0.10031948963187563, 0.1004987562112089, 0.10259142264341595, 0.10270345661174213, 0.1039471019317037,  
                                        0.10270345661174213, 0.10196568050084304, 0.10932977636490435, 0.10413932974625868, 0.10235721762533408,
                                        0.1078563859954523, 0.10412012293500235, 0.1077032961426901, 0.11764777940955791, 0.10736125930707034, 
                                        0.13917614738165443, 0.14509651959988565, 0.1776091213873882, 0.14606854555310667, 0.17462244987400674 };
    double p8424_d3x1y1_ystatminus[]= { 9.0E-4, 0.001, 0.0011, 0.0012, 0.0013, 0.0014, 0.0015, 0.0017, 0.0018, 
                                        0.0019, 0.0015, 0.0017, 0.002, 0.0023, 0.0026, 0.003, 0.0034, 0.0038, 0.0043, 
                                        0.0047, 0.0039, 0.005, 0.006, 0.007, 0.008, 0.008, 0.01, 0.011, 0.012, 
                                        0.014, 0.012, 0.014, 0.017, 0.021, 0.026, 0.032, 0.029, 0.04, 0.055, 
                                        0.064, 0.089, 0.107, 0.141, 0.1159, 0.138 };
    double p8424_d3x1y1_ystatplus[] = { 9.0E-4, 0.001, 0.0011, 0.0012, 0.0013, 0.0014, 0.0015, 0.0017, 0.0018, 
                                        0.0019, 0.0015, 0.0017, 0.002, 0.0023, 0.0026, 0.003, 0.0034, 0.0038, 0.0043, 
                                        0.0047, 0.0039, 0.005, 0.006, 0.007, 0.008, 0.008, 0.01, 0.011, 0.012, 
                                        0.014, 0.012, 0.014, 0.017, 0.021, 0.026, 0.032, 0.029, 0.04, 0.055, 
                                        0.064, 0.089, 0.107, 0.141, 0.1159, 0.138 };
    int p8424_d3x1y1_numpoints      = 45;
    TGraphAsymmErrors* p8424_d3x1y1;
    if( typeErr.CompareTo("Syst") == 0 ) {
        p8424_d3x1y1 = new TGraphAsymmErrors(p8424_d3x1y1_numpoints, p8424_d3x1y1_xval, p8424_d3x1y1_yval, p8424_d3x1y1_xerrminus, p8424_d3x1y1_xerrplus, p8424_d3x1y1_yerrminus, p8424_d3x1y1_yerrplus);
    } else if ( typeErr.CompareTo("Stat") == 0 ) {
        p8424_d3x1y1 = new TGraphAsymmErrors(p8424_d3x1y1_numpoints, p8424_d3x1y1_xval, p8424_d3x1y1_yval, p8424_d3x1y1_xerrminus, p8424_d3x1y1_xerrplus, p8424_d3x1y1_ystatminus, p8424_d3x1y1_ystatplus);
    }
    p8424_d3x1y1->SetName("/HepData/8424/d3x1y1");
    p8424_d3x1y1->SetTitle("/HepData/8424/d3x1y1");
    return p8424_d3x1y1;
}

TGraphAsymmErrors* GetChargeParticlesRpPb2013(TString typeErr){
    // Plot: p8550_d4x1y2
    double p8550_d4x1y1_xval[]      = { 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 
                                        0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.05, 1.15, 
                                        1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 
                                        2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.25, 4.75, 
                                        5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 
                                        13.5, 14.5, 15.5, 17.0, 19.0, 21.0, 23.0, 26.0, 30.0, 36.0, 
                                        45.0 };
    double p8550_d4x1y1_xerrminus[] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                        0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 
                                        0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.100, 0.100, 
                                        0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.25, 0.25, 
                                        0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                                        0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 4.0, 
                                        5.0 };
    double p8550_d4x1y1_xerrplus[]  = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                        0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 
                                        0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.100, 0.100, 
                                        0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.25, 0.25, 
                                        0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                                        0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 4.0, 
                                        5.0 };
    double p8550_d4x1y1_yval[]      = { 0.5346, 0.5265, 0.5276, 0.5324, 0.5399, 0.5505, 0.5646, 0.5811, 0.5975, 
                                        0.6143, 0.6316, 0.6496, 0.667, 0.6845, 0.6986, 0.7119, 0.7252, 0.7499, 0.7806, 
                                        0.8074, 0.8354, 0.86, 0.8836, 0.9024, 0.9284, 0.9492, 0.9732, 0.9932, 1.021, 
                                        1.049, 1.065, 1.075, 1.083, 1.087, 1.082, 1.082, 1.093, 1.099, 1.084, 
                                        1.083, 1.07, 1.065, 1.035, 1.027, 1.027, 0.9953, 1.001, 0.9924, 0.9966, 
                                        0.9812, 1.001, 0.9878, 0.9951, 0.9866, 0.9764, 1.048, 0.9631, 0.9364, 1.108, 
                                        0.8506 };
    double p8550_d4x1y1_yerrminus[] = { 0.05200104710484203, 0.042000524996718795, 0.04300046511376359, 0.04500040110932346, 0.048000416664858235, 
                                        0.05100039215535504, 0.052000465382532876, 0.05300054339344079, 0.053000637732012246, 0.052000753840689654, 
                                        0.05300090659601965, 0.05500105089905101, 0.05600128927087304, 0.05700154734741856, 0.058001904279083805, 
                                        0.059002204196114565, 0.060002613276423214, 0.06200129030915405, 0.06500155767364349, 0.070002160680939, 
                                        0.0730025485856487, 0.07100385623330609, 0.07300449369730606, 0.07600593726808452, 0.07700649323271383, 
                                        0.07800923022309604, 0.07901069547852367, 0.08100888840121187, 0.08301018009858792, 0.08501699830033993, 
                                        0.08602098581160297, 0.08704188646852733, 0.08803272118933958, 0.08804141071109663, 0.08805112151472007, 
                                        0.08806185326235191, 0.08807774974418908, 0.0890649201425567, 0.08808637806153684, 0.0941029223775755, 
                                        0.09216425554410994, 0.0912463149940862, 0.08834500551813894, 0.08728923186739587, 0.08755021416307329, 
                                        0.08485281374238571, 0.085, 0.08551023330572781, 0.08514693182963201, 0.0848999411071645, 
                                        0.0891852005660132, 0.08953770155638349, 0.0904433524367601, 0.09546203433826454, 0.09888377015466188,
                                        0.11399122773266371, 0.10430723848324237, 0.11745637488020819, 0.14142135623730953, 0.14477914214416385 };
    double p8550_d4x1y1_yerrplus[]  = { 0.05200104710484203, 0.042000524996718795, 0.04300046511376359, 0.04500040110932346, 0.048000416664858235, 
                                        0.05100039215535504, 0.052000465382532876, 0.05300054339344079, 0.053000637732012246, 0.052000753840689654, 
                                        0.05300090659601965, 0.05500105089905101, 0.05600128927087304, 0.05700154734741856, 0.058001904279083805, 
                                        0.059002204196114565, 0.060002613276423214, 0.06200129030915405, 0.06500155767364349, 0.070002160680939, 
                                        0.0730025485856487, 0.07100385623330609, 0.07300449369730606, 0.07600593726808452, 0.07700649323271383, 
                                        0.07800923022309604, 0.07901069547852367, 0.08100888840121187, 0.08301018009858792, 0.08501699830033993, 
                                        0.08602098581160297, 0.08704188646852733, 0.08803272118933958, 0.08804141071109663, 0.08805112151472007, 
                                        0.08806185326235191, 0.08807774974418908, 0.0890649201425567, 0.08808637806153684, 0.0941029223775755, 
                                        0.09216425554410994, 0.0912463149940862, 0.08834500551813894, 0.08728923186739587, 0.08755021416307329, 
                                        0.08485281374238571, 0.085, 0.08551023330572781, 0.08514693182963201, 0.0848999411071645, 
                                        0.0891852005660132, 0.08953770155638349, 0.0904433524367601, 0.09546203433826454, 0.09888377015466188, 
                                        0.11399122773266371, 0.10430723848324237, 0.11745637488020819, 0.14142135623730953, 0.14477914214416385 };
    double p8550_d4x1y1_ystatminus[]= { 3.3E-4, 2.1E-4, 2.0E-4, 1.9E-4, 2.0E-4, 2.0E-4, 2.2E-4, 2.4E-4, 2.6E-4, 
                                        2.8E-4, 3.1E-4, 3.4E-4, 3.8E-4, 4.2E-4, 4.7E-4, 5.1E-4, 5.6E-4, 4.0E-4, 4.5E-4, 
                                        5.5E-4, 6.1E-4, 7.4E-4, 8.1E-4, 9.5E-4, 0.001, 0.0012, 0.0013, 0.0012, 0.0013, 
                                        0.0017, 0.0019, 0.0027, 0.0024, 0.0027, 0.003, 0.0033, 0.0037, 0.0034, 0.0039, 
                                        0.0044, 0.0055, 0.0067, 0.0078, 0.0071, 0.0098, 0.012, 0.013, 0.016, 0.019, 
                                        0.022, 0.027, 0.031, 0.028, 0.037, 0.047, 0.063, 0.056, 0.08, 0.1, 
                                        0.12 };
    double p8550_d4x1y1_ystatplus[] = { 3.3E-4, 2.1E-4, 2.0E-4, 1.9E-4, 2.0E-4, 2.0E-4, 2.2E-4, 2.4E-4, 2.6E-4, 
                                        2.8E-4, 3.1E-4, 3.4E-4, 3.8E-4, 4.2E-4, 4.7E-4, 5.1E-4, 5.6E-4, 4.0E-4, 4.5E-4, 
                                        5.5E-4, 6.1E-4, 7.4E-4, 8.1E-4, 9.5E-4, 0.001, 0.0012, 0.0013, 0.0012, 0.0013, 
                                        0.0017, 0.0019, 0.0027, 0.0024, 0.0027, 0.003, 0.0033, 0.0037, 0.0034, 0.0039, 
                                        0.0044, 0.0055, 0.0067, 0.0078, 0.0071, 0.0098, 0.012, 0.013, 0.016, 0.019, 
                                        0.022, 0.027, 0.031, 0.028, 0.037, 0.047, 0.063, 0.056, 0.08, 0.1, 
                                        0.12 };
    int p8550_d4x1y1_numpoints      = 60;
    TGraphAsymmErrors* p8550_d4x1y1;
    
    if( typeErr.CompareTo("Syst") == 0 ) {
        p8550_d4x1y1 = new TGraphAsymmErrors(p8550_d4x1y1_numpoints, p8550_d4x1y1_xval, p8550_d4x1y1_yval, p8550_d4x1y1_xerrminus, p8550_d4x1y1_xerrplus, p8550_d4x1y1_yerrminus, p8550_d4x1y1_yerrplus);
    } else if ( typeErr.CompareTo("Stat") == 0 ) {
        p8550_d4x1y1  = new TGraphAsymmErrors(p8550_d4x1y1_numpoints, p8550_d4x1y1_xval, p8550_d4x1y1_yval, p8550_d4x1y1_xerrminus, p8550_d4x1y1_xerrplus, p8550_d4x1y1_ystatminus, p8550_d4x1y1_ystatplus);
    }
    p8550_d4x1y1->SetName("/HepData/8550/d4x1y1");
    p8550_d4x1y1->SetTitle("/HepData/8550/d4x1y1");
    
    return p8550_d4x1y1;
}

TGraphErrors* RemovePointsFromGraph(TGraphErrors *graph, Int_t NbPoints){
    
    TGraphErrors* dummyGraph= (TGraphErrors*) graph->Clone();
    Double_t * xValue       = dummyGraph->GetX();
    Double_t * yValue       = dummyGraph->GetY();
    Double_t* xError        = dummyGraph->GetEX();
    Double_t* yError        = dummyGraph->GetEY();
    Int_t nPoints           = dummyGraph->GetN();
    Int_t nPoints_new       = nPoints  - NbPoints;
    
    for (Int_t i = NbPoints; i < nPoints; i++){
        yValue[i-NbPoints]  = yValue[i];
        xValue[i-NbPoints]  = xValue[i];
        xError[i-NbPoints]  = xError[i];
        yError[i-NbPoints]  = yError[i];
    }
    TGraphErrors* returnGraph = new TGraphErrors(nPoints_new,xValue,yValue, xError,yError);
    return returnGraph;
}

TF1* RebinWithFitToTGraph(TGraphAsymmErrors *spectrum, TGraphAsymmErrors** newSpectrum, TGraphAsymmErrors *newBins,TF1* CurrentFit, TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters, Int_t fixParNumber){
    
    Int_t     nOldPoints        = spectrum->GetN();
    Double_t *xOldValue         = spectrum->GetX();
    Double_t *yOldValueErrlow   = spectrum->GetEYlow();
    Double_t *yOldValueErrhigh  = spectrum->GetEYhigh();
    Double_t *yOldValue         = spectrum->GetY();
    
    
    ///////////////////////////////////////////////////////////
    Int_t nNewPoints            = newBins->GetN();
    Double_t *xNewValue         = newBins->GetX();
    Double_t *xNewValueErrlow   = newBins->GetEXlow();
    Double_t *xNewValueErrhigh  = newBins->GetEXhigh();
    //////////////////////////////////////////////////////////////    
    
 
    
    
    Double_t  yNewValueRelErrlow[nNewPoints];
    Double_t  yNewValueRelErrhigh[nNewPoints];
    (*newSpectrum)              = new TGraphAsymmErrors (nNewPoints);
    
    if( CurrentFit == 0) {
        CurrentFit = new TF1();
        if( FitType.BeginsWith("l") || FitType.BeginsWith("L") ) {
            CurrentFit = FitObject("l","fitInvCrossSectionPi0","Pi0");
            CurrentFit->SetRange(minPt,maxPt);
            CurrentFit->SetParameters(parameters[0],parameters[1],parameters[2]); // standard
    
            if( fixParNumber > -1 && fixParNumber < 3 ){
                CurrentFit->FixParameter(fixParNumber,parameters[fixParNumber]);
            } else if ( fixParNumber > -1 ){
                cout<<"WARNING: the number of parameter is wrong"<<endl;
            }
            spectrum->Fit(CurrentFit,"SQNRME+","",minPt,maxPt);  //One time
        } else if (FitType.BeginsWith("tcm") || FitType.BeginsWith("TCM") ){
            CurrentFit = FitObject("tcm","fitInvCrossSectionPi0","Pi0");
            CurrentFit->SetRange(minPt,maxPt);
            CurrentFit->SetParameters(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]); // standard
            spectrum->Fit(CurrentFit,"SQNRME+","",minPt,maxPt);  //One time
        }
    } 
    cout<<"iPoint"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"xElow"<<"\t"<<"xEhigh"<<"\t"<<"yElow"<<"\t"<<"yEhigh"<<endl;
    
    Double_t decisionBoundary           = 0.0000001;
        
    for( Int_t iPoint = 0; iPoint < nNewPoints; iPoint++){
        Double_t yEval;
        Int_t index                     = 0;
        
        while(  ( xNewValue[iPoint] - xOldValue[index] ) > decisionBoundary  && index < nOldPoints-1 ) {index++;}
        
        cout<<"xNewValue "<<xNewValue[iPoint]<<" xOldValue "<<xOldValue[index]<<"  "<<(xNewValue[iPoint] - xOldValue[index])<<" nOldPoints "<<nOldPoints-1<<" index "<<index<<endl;
        
        if( index == 0 ){
            yNewValueRelErrlow[iPoint]  = yOldValueErrlow[index]  / yOldValue[index];
            yNewValueRelErrhigh[iPoint] = yOldValueErrhigh[index] / yOldValue[index];
        } else if (  TMath::Abs(  xNewValue[iPoint] - xOldValue[index-1] ) < TMath::Abs(  xNewValue[iPoint] - xOldValue[index] ) ) {
            yNewValueRelErrlow[iPoint]  = yOldValueErrlow[index-1] / yOldValue[index-1];
            yNewValueRelErrhigh[iPoint] = yOldValueErrhigh[index-1] / yOldValue[index-1];
        } else {
            yNewValueRelErrlow[iPoint]  = yOldValueErrlow[index]  / yOldValue[index]; 
            yNewValueRelErrhigh[iPoint] = yOldValueErrhigh[index] / yOldValue[index];
        }
        yEval                           = CurrentFit->Eval(xNewValue[iPoint]);
        Double_t yNewValueErrlow        = yEval * yNewValueRelErrlow[iPoint];
        Double_t yNewValueErrhigh       = yEval * yNewValueRelErrhigh[iPoint];
        (*newSpectrum)->SetPoint(iPoint,xNewValue[iPoint],yEval);
        (*newSpectrum)->SetPointError(iPoint,xNewValueErrlow[iPoint],xNewValueErrhigh[iPoint],yNewValueErrlow,yNewValueErrhigh);
        
        cout<<iPoint<<"\t"<<xNewValue[iPoint]<<"\t"<<yEval<<"\t"<<xNewValueErrlow[iPoint]<<"\t"<<xNewValueErrhigh[iPoint]<<"\t"<<yNewValueErrlow<<"\t"<<yNewValueErrhigh<<"\t"<<yNewValueRelErrlow[iPoint]<<endl;    
    }
    return CurrentFit;
}


void ExtrapolateSpectrum(TF1* FitToSpectrum, TGraphAsymmErrors *spectrumStatErr, TGraphAsymmErrors *spectrumSystErr, TGraphAsymmErrors** newSpectrumStatErr, TGraphAsymmErrors** newSpectrumSystErr, TString system, TString energy, TString meson){

    
    
    (*newSpectrumStatErr) = (TGraphAsymmErrors*) spectrumStatErr->Clone();
    (*newSpectrumSystErr) = (TGraphAsymmErrors*) spectrumSystErr->Clone();
    
    Double_t newXBins[20];
    Double_t newXerrlow[20];
    Double_t newXerrhigh[20];
    Double_t newYerrlowStat[20];
    Double_t newYerrhighStat[20];
    Double_t newYerrlowSyst[20];
    Double_t newYerrhighSyst[20];
    
    Int_t nNewPoints = -1;
    
    //Initialation
    for(Int_t i=0; i<20; i++){
        newXBins[i]=0.;
        newXerrhigh[i]=0.;
        newXerrlow[i]=0.;
    }
    
   
    

    
    if( system.CompareTo("PCM") == 0  && energy.CompareTo("2.76TeV") == 0 && meson.CompareTo("Eta") == 0 ) {
        
        nNewPoints = 3;    
        newXBins[0]      = 5.5; newXBins[1] = 7.0;newXBins[2] = 9.0;
        newXerrhigh[0]   = 0.5; newXerrlow[0] = 0.5;
        newXerrhigh[1]   = 1.0; newXerrlow[1] = 1.0;
        newXerrhigh[2]   = 1.0; newXerrhigh[2] = 1.0;
        
    
        
    } else if ( system.CompareTo("PCM") == 0 && ( energy.CompareTo("7TeV") == 0 || energy.CompareTo("8TeV") == 0 ) && meson.CompareTo("Eta") == 0 ) {
        
        nNewPoints       = 1;    
        newXBins[0]      = 9.0; 
        newXerrhigh[0]   = 1.; newXerrlow[0] = 1.;
        
        
    } else if ( system.CompareTo("PCM") == 0  && energy.CompareTo("2.76TeV") == 0 && meson.CompareTo("Pi0") == 0 ) {
        
    

        nNewPoints = 4; 
        newXBins[0]      = 7.5; newXBins[1] = 9;newXBins[2] = 11; newXBins[3] = 14;
        newXerrhigh[0]   = 0.5; newXerrlow[0] = 0.5;
        newXerrhigh[1]   = 1.0; newXerrlow[1] = 1.0;
        newXerrhigh[2]   = 1.0; newXerrlow[2] = 1.0;
        newXerrhigh[3]   = 2.0; newXerrlow[3] = 2.0;
        
        
    } else if ( system.CompareTo("PHOS") ==0 && energy.CompareTo("2.76TeV") == 0 && meson.CompareTo("Pi0") == 0 ) {
        
        nNewPoints       = 2;    
        newXBins[0]      = 14.0; 
        newXBins[1]      = 18.0;
        newXerrhigh[0]   = 2.; newXerrlow[0] = 2.;
        newXerrhigh[1]   = 2.; newXerrlow[1] = 2.;
        
        
    }

   
    
    
    for(Int_t i=0; i<nNewPoints; i++) {
        newYerrhighStat[i]   = spectrumStatErr->GetEYhigh()[spectrumStatErr->GetN()-1] / spectrumStatErr->GetY()[spectrumStatErr->GetN()-1]; 
        newYerrlowStat[i]    = spectrumStatErr->GetEYlow()[spectrumStatErr->GetN()-1] / spectrumStatErr->GetY()[spectrumStatErr->GetN()-1]; 
        newYerrhighSyst[i]   = spectrumSystErr->GetEYhigh()[spectrumSystErr->GetN()-1] / spectrumSystErr->GetY()[spectrumSystErr->GetN()-1]; ; 
        newYerrlowSyst[i]    = spectrumSystErr->GetEYlow()[spectrumSystErr->GetN()-1] / spectrumSystErr->GetY()[spectrumSystErr->GetN()-1]; ; 
    }
    
    
    
    
    
    Int_t currentPoint = (*newSpectrumSystErr)->GetN();
    
    
    for(Int_t i=0; i<nNewPoints; i++){
       
        Double_t yEval      = FitToSpectrum->Eval(newXBins[i]);
        newYerrlowStat[i]   = yEval* newYerrlowStat[i];
        newYerrhighStat[i]  = yEval* newYerrhighStat[i];
        newYerrlowSyst[i]   = yEval* newYerrlowSyst[i];
        newYerrhighSyst[i]  = yEval* newYerrhighSyst[i];
        
        (*newSpectrumStatErr)->SetPoint(currentPoint,newXBins[i],yEval);
        (*newSpectrumStatErr)->SetPointError(currentPoint,newXerrlow[i],newXerrhigh[i],newYerrlowStat[i],newYerrhighStat[i]);
        (*newSpectrumSystErr)->SetPoint(currentPoint,newXBins[i],yEval);
        (*newSpectrumSystErr)->SetPointError(currentPoint,newXerrlow[i],newXerrhigh[i],newYerrlowStat[i],newYerrhighStat[i]);
        
        currentPoint++;
        
    }
    

    
}



TF1* RebinWithFitToTGraphWithUpDownYShifted(TGraphAsymmErrors *spectrumStatErr, TGraphAsymmErrors *spectrumSystErr,TGraphAsymmErrors *newBins,TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters,  TGraphAsymmErrors** newSpectrumStatErr,
					     TGraphAsymmErrors** newSpectrumSystErr,TGraphAsymmErrors** spectrumYShiftedDownStatErr,TGraphAsymmErrors** spectrumYShiftedUpStatErr,
					     TF1** FitToSpectrumYShiftedDown, TF1** FitToSpectrumYShiftedUp,TString meson,TString energy,TString system){

    //Shift the spectra up and down acording to the systematic errors
   
    (*spectrumYShiftedUpStatErr)    = BinYShiftwithErrhigh( spectrumStatErr,spectrumSystErr);
    (*spectrumYShiftedDownStatErr)  = BinYShiftwithErrlow(  spectrumStatErr,spectrumSystErr);
    
    
    
    const Int_t nNewPoints        = newBins->GetN();
    Double_t *xNewValue           = newBins->GetX();
    Double_t *xNewValueErrlow     = newBins->GetEXlow();
    Double_t *xNewValueErrhigh    = newBins->GetEXhigh();
  
   
    //Preparing fit functions
   TF1* FitToSpectrum 		  = 0;
   
   TGraphErrors *grint   = new TGraphErrors(nNewPoints);
       
   for (Int_t i=0; i<nNewPoints; i++)
        grint->SetPoint(i, xNewValue[i], 0);
   
   
       
   if( FitType.BeginsWith("l") || FitType.BeginsWith("L") ) {
     
            FitToSpectrum = FitObject("l",Form("fitInvCrossSection%s",meson.Data()),meson.Data());
            FitToSpectrum->SetRange(minPt,maxPt);
            FitToSpectrum->SetParameters(parameters[0],parameters[1],parameters[2]); // standard
	    spectrumStatErr->Fit(FitToSpectrum,"SNRME+","",minPt,maxPt);  //One time
            (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
            
            //parameters[0]= FitToSpectrum->GetParameter(0);
	    //parameters[1]= FitToSpectrum->GetParameter(1);
	    //parameters[2]= FitToSpectrum->GetParameter(2);
	    
        
	    //To compute systematic errors
	    (*FitToSpectrumYShiftedUp) = FitObject("l",Form("fitInvCrossSection%sYShiftedUp",meson.Data()),meson.Data());
            (*FitToSpectrumYShiftedUp)->SetRange(minPt,maxPt);
            (*FitToSpectrumYShiftedUp)->SetParameters(parameters[0],parameters[1],parameters[2]); // standard
            (*spectrumYShiftedUpStatErr)->Fit((*FitToSpectrumYShiftedUp),"SNRME+","",minPt,maxPt);  //One time
	    
	    
	    
	    (*FitToSpectrumYShiftedDown) = FitObject("l",Form("fitInvCrossSection%sYShiftedDown",meson.Data()),meson.Data());
            (*FitToSpectrumYShiftedDown)->SetRange(minPt,maxPt);
            (*FitToSpectrumYShiftedDown)->SetParameters(parameters[0] ,parameters[1],parameters[2]); // standard 
            (*FitToSpectrumYShiftedDown)->SetParLimits(0,0.8*parameters[0],1.2*parameters[0]);
	    (*FitToSpectrumYShiftedDown)->SetParLimits(1,0.8*parameters[1],1.2*parameters[1]);
            (*spectrumYShiftedDownStatErr)->Fit((*FitToSpectrumYShiftedDown),"SNRME+","",minPt,maxPt);  //One time
            
             cout<<(*FitToSpectrumYShiftedDown)->GetParameter(0)<<endl;
             cout<<(*FitToSpectrumYShiftedDown)->GetParameter(1)<<endl;
             cout<<(*FitToSpectrumYShiftedDown)->GetParameter(2)<<endl;
             cout<<(*FitToSpectrumYShiftedDown)->GetParameter(1)<<endl;
            
	   
   } else if ( FitType.BeginsWith("tcm") || FitType.BeginsWith("tcm") ){
            
            FitToSpectrum = FitObject("tcm",Form("fitInvCrossSection%s",meson.Data()),meson.Data());
            FitToSpectrum->SetRange(minPt,maxPt);
            FitToSpectrum->SetParameters(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]); // standard
            
            if( meson.CompareTo("Pi0") == 0 && system.CompareTo("PCM-EMCal") == 0 ) {
            FitToSpectrum->SetParLimits(0,0.6*parameters[0],1.2*parameters[0]);
            } else {
            FitToSpectrum->SetParLimits(0,0.7*parameters[0],1.2*parameters[0]);   
            }
	    FitToSpectrum->SetParLimits(1,0.8*parameters[1],1.2*parameters[1]);
            FitToSpectrum->SetParLimits(2,0.8*parameters[2],1.2*parameters[2]);	    
            FitToSpectrum->SetParLimits(3,0.8*parameters[3],1.2*parameters[3]);	    
            FitToSpectrum->SetParLimits(4,0.8*parameters[4],1.2*parameters[4]);
            
            
            spectrumStatErr->Fit(FitToSpectrum,"SNRME+","",minPt,maxPt);  //One time
            //spectrumStatErr->Fit(FitToSpectrum,"SNRME+","",minPt,maxPt);  //One time
            (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
  
            
            parameters[0]= FitToSpectrum->GetParameter(0);
	    parameters[1]= FitToSpectrum->GetParameter(1);
	    parameters[2]= FitToSpectrum->GetParameter(2);
	    parameters[3]= FitToSpectrum->GetParameter(3);
	    parameters[4]= FitToSpectrum->GetParameter(4);
	    
	    //To compute systematic errors
	    (*FitToSpectrumYShiftedUp) = FitObject("tcm",Form("fitInvCrossSection%sYShiftedUp",meson.Data()),meson.Data());
            (*FitToSpectrumYShiftedUp)->SetRange(minPt,maxPt);
            (*FitToSpectrumYShiftedUp)->SetParameters(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]); // standard
            (*FitToSpectrumYShiftedUp)->SetParLimits(0,0.7*parameters[0],1.3*parameters[0]);
	    (*FitToSpectrumYShiftedUp)->SetParLimits(1,0.7*parameters[1],1.3*parameters[1]);
           
            (*spectrumYShiftedUpStatErr)->Fit((*FitToSpectrumYShiftedUp),"SNRME+","",minPt,maxPt);  //One time
            //(*spectrumYShiftedUpStatErr)->Fit((*FitToSpectrumYShiftedUp),"SNRME+","",minPt,maxPt);  //One time
	   
	    
	    //To compute systematic errors
	    
	    (*FitToSpectrumYShiftedDown) = FitObject("tcm",Form("fitInvCrossSection%sYShiftedDown",meson.Data()),meson.Data());
            (*FitToSpectrumYShiftedDown)->SetRange(minPt,maxPt);
            (*FitToSpectrumYShiftedDown)->SetParameters(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]); //standard
            (*FitToSpectrumYShiftedDown)->SetParLimits(0,0.7*parameters[0],1.3*parameters[0]);
	    (*FitToSpectrumYShiftedDown)->SetParLimits(1,0.7*parameters[1],1.3*parameters[1]);
            (*spectrumYShiftedDownStatErr)->Fit((*FitToSpectrumYShiftedDown),"SNRME+","",minPt,maxPt);  //One time
            //(*spectrumYShiftedDownStatErr)->Fit((*FitToSpectrumYShiftedDown),"SNRME+","",minPt,maxPt);  //One time
       
       
   }
   
    (*newSpectrumStatErr) = new TGraphAsymmErrors(nNewPoints);
    
    (*newSpectrumSystErr) = new TGraphAsymmErrors(nNewPoints);
    
    for( Int_t iPoint = 0; iPoint < nNewPoints; iPoint++){
        
        Double_t yEval      = FitToSpectrum->Eval(xNewValue[iPoint]);
	
	Double_t yNewValueStatErrlow    = grint->GetEY()[iPoint]; //Getting the confidencial error from the fit
	Double_t yNewValueStatErrhigh   = grint->GetEY()[iPoint]; //Getting the confidencial error from the fit
  
        (*newSpectrumStatErr)->SetPoint(iPoint,xNewValue[iPoint],yEval);
        (*newSpectrumStatErr)->SetPointError(iPoint,xNewValueErrlow[iPoint],xNewValueErrhigh[iPoint],yNewValueStatErrlow,yNewValueStatErrhigh);
        cout<<"Stat "<<iPoint<<"\t"<<xNewValue[iPoint]<<"\t"<<yEval<<"\t"<<xNewValueErrlow[iPoint]<<"\t"<<xNewValueErrhigh[iPoint]<<"\t"<<yNewValueStatErrlow<<"\t"<<yNewValueStatErrhigh<<endl;
	
	
	Double_t yEvalUp    = (*FitToSpectrumYShiftedUp)->Eval(xNewValue[iPoint]);
	Double_t yEvalDown  = (*FitToSpectrumYShiftedDown)->Eval(xNewValue[iPoint]);
	
	Double_t yNewValueSystErrlow	= 0;
	Double_t yNewValueSystErrhigh   = 0;
	
	Double_t yErrUp   = 0;
	Double_t yErrDown = 0;
	
	if( yEval != 0 ){
	
	yErrUp   = TMath::Abs( (yEval-yEvalUp)/yEval );
	yErrDown = TMath::Abs( (yEval-yEvalDown)/yEval);
	
	}
	
	cout<<"Syst "<<iPoint<<"  ErrUp "<<yErrUp<<"  ErrDown "<<yErrDown<<endl;
	
	
	if( yErrUp >= yErrDown){
	  
	yNewValueSystErrlow  = yErrUp*yEval;
	yNewValueSystErrhigh = yErrUp*yEval;
	  
	} else {
	  
	yNewValueSystErrlow  = yErrDown*yEval;
	yNewValueSystErrhigh = yErrDown*yEval;
	  
	}
	
	(*newSpectrumSystErr)->SetPoint(iPoint,xNewValue[iPoint],yEval);
        (*newSpectrumSystErr)->SetPointError(iPoint,xNewValueErrlow[iPoint],xNewValueErrhigh[iPoint],yNewValueSystErrlow,yNewValueSystErrhigh);
        cout<<"Syst "<<iPoint<<"\t"<<xNewValue[iPoint]<<"\t"<<yEval<<"\t"<<xNewValueErrlow[iPoint]<<"\t"<<xNewValueErrhigh[iPoint]<<"\t"<<yNewValueSystErrlow<<"\t"<<yNewValueSystErrhigh<<endl;
	
	
    }
    
    return FitToSpectrum;
}






TF1* RebinWithFitToTGraphWithMeanErr(TGraphAsymmErrors *spectrum, TGraphAsymmErrors** newSpectrum, TGraphAsymmErrors *newBins,TF1* CurrentFit, TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters,Int_t fixParNumber, TString meson){
    
    Int_t     nOldPoints       = spectrum->GetN();
    Double_t *xOldValue        = spectrum->GetX();
    Double_t *yOldValueErrlow  = spectrum->GetEYlow();
    Double_t *yOldValueErrhigh = spectrum->GetEYhigh();
    Double_t *yOldValue        = spectrum->GetY();    
    ///////////////////////////////////////////////////////////
    Int_t nNewPoints            = newBins->GetN();
    Double_t *xNewValue         = newBins->GetX();
    Double_t *xNewValueErrlow   = newBins->GetEXlow();
    Double_t *xNewValueErrhigh  = newBins->GetEXhigh();
    //////////////////////////////////////////////////////////////    
    Double_t  yNewValueRelErrlow[nNewPoints];
    Double_t  yNewValueRelErrhigh[nNewPoints];
    
    (*newSpectrum)        = new TGraphAsymmErrors (nNewPoints);
    
    if( CurrentFit == 0) {
        CurrentFit = new TF1();
        if( FitType.BeginsWith("l") || FitType.BeginsWith("L") ) {
            CurrentFit = FitObject("l","fitInvCrossSectionPi0",meson.Data());
            CurrentFit->SetRange(minPt,maxPt);
            CurrentFit->SetParameters(parameters[0],parameters[1],parameters[2]); // standard
            
            if( fixParNumber > -1 && fixParNumber < 3 ){
                CurrentFit->FixParameter(fixParNumber,parameters[fixParNumber]);
            } else if ( fixParNumber > -1 ){
            cout<<"WARNING: the number of parameter is wrong"<<endl;
            }
            spectrum->Fit(CurrentFit,"SQNRME+","",minPt,maxPt);  //One time
        } else if (FitType.BeginsWith("tcm") || FitType.BeginsWith("TCM") ){
            CurrentFit = FitObject("tcm","fitInvCrossSectionPi0",meson.Data());
            CurrentFit->SetRange(minPt,maxPt);
            CurrentFit->SetParameters(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]); // standard
            spectrum->Fit(CurrentFit,"SQNRME+","",minPt,maxPt);  //One time
        }
    } 
    cout<<"iPoint"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"xElow"<<"\t"<<"xEhigh"<<"\t"<<"yElow"<<"\t"<<"yEhigh"<<endl;
    Double_t decisionBoundary = 0.0000001;
        
    for( Int_t iPoint = 0; iPoint < nNewPoints; iPoint++){
        Double_t yEval;
        Int_t index         = 0;

        while(  ( xNewValue[iPoint] - xOldValue[index] ) > decisionBoundary  && index < nOldPoints-1 ) {index++;}
    
        cout<<"xNewValue "<<xNewValue[iPoint]<<" xOldValue "<<xOldValue[index]<<"  "<<(xNewValue[iPoint] - xOldValue[index])<<" nOldPoints "<<nOldPoints-1<<" index "<<index<<endl;
        
        if( index == 0 ){
            yNewValueRelErrlow[iPoint]  = yOldValueErrlow[index]  / yOldValue[index];
            yNewValueRelErrhigh[iPoint] = yOldValueErrhigh[index] / yOldValue[index];
        } else {
            Double_t meanErrlow         = ((yOldValueErrlow[index]  / yOldValue[index])  +  ( yOldValueErrlow[index-1]   / yOldValue[index-1] ) ) / 2 ;
            Double_t meanErrhigh        = ((yOldValueErrhigh[index] / yOldValue[index])  +  ( yOldValueErrhigh[index-1]  / yOldValue[index-1] ) ) / 2 ;
            yNewValueRelErrlow[iPoint]  = meanErrlow;
            yNewValueRelErrhigh[iPoint] = meanErrlow;
        }
        
        yEval                           = CurrentFit->Eval(xNewValue[iPoint]);
        Double_t yNewValueErrlow        = yEval * yNewValueRelErrlow[iPoint];
        Double_t yNewValueErrhigh       = yEval * yNewValueRelErrhigh[iPoint];
        
        (*newSpectrum)->SetPoint(iPoint,xNewValue[iPoint],yEval);
        (*newSpectrum)->SetPointError(iPoint,xNewValueErrlow[iPoint],xNewValueErrhigh[iPoint],yNewValueErrlow,yNewValueErrhigh);
        cout<<iPoint<<"\t"<<xNewValue[iPoint]<<"\t"<<yEval<<"\t"<<xNewValueErrlow[iPoint]<<"\t"<<xNewValueErrhigh[iPoint]<<"\t"<<yNewValueErrlow<<"\t"<<yNewValueErrhigh<<"\t"<<yNewValueRelErrlow[iPoint]<<endl;
    }
    return CurrentFit;
}


TF1* FillTGraphEYWithFitErr(TGraphAsymmErrors *spectrum, TGraphAsymmErrors** newSpectrum, TGraphAsymmErrors *newBins,TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters, TString meson){
    //This function compute the errors of Y from the fit function

    ///////////////////////////////////////////////////////////    
    Int_t nNewPoints            = newBins->GetN();
    Double_t *xNewValue         = newBins->GetX();
    Double_t *yNewValue         = newBins->GetY();
    Double_t *xNewValueErrlow   = newBins->GetEXlow();
    Double_t *xNewValueErrhigh  = newBins->GetEXhigh();
    
    //////////////////////////////////////////////////////////////    
    //We are going to obtain the values mentioned above from the fit
    
    Double_t  yNewValueRelErrlow[nNewPoints];
    Double_t  yNewValueRelErrhigh[nNewPoints];
    
    (*newSpectrum)              = new TGraphAsymmErrors (nNewPoints);
    TGraphErrors *grint         = new TGraphErrors(nNewPoints);
    
    for (Int_t i=0; i<nNewPoints; i++)
        grint->SetPoint(i, xNewValue[i], 0);
    
    TF1* CurrentFit             = new TF1();
    
    if( FitType.BeginsWith("l") || FitType.BeginsWith("L") ) {
        CurrentFit = FitObject("l","fitInvCrossSectionPi0",meson.Data());
        CurrentFit->SetRange(minPt,maxPt);
        CurrentFit->SetParameters(parameters[0],parameters[1],parameters[2]); // standard
        spectrum->Fit(CurrentFit,"SQNRME+","",minPt,maxPt);  //One time

        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
    } else if (FitType.BeginsWith("tcm") || FitType.BeginsWith("TCM") ){
        CurrentFit = FitObject("tcm","fitInvCrossSectionPi0",meson.Data());
        CurrentFit->SetRange(minPt,maxPt);
        CurrentFit->SetParameters(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]); // standard
        spectrum->Fit(CurrentFit,"SQNRME+","",minPt,maxPt);  //One time
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
    }
    cout<<"iPoint"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"xElow"<<"\t"<<"xEhigh"<<"\t"<<"yElow"<<"\t"<<"yEhigh"<<endl;
       
    for( Int_t iPoint = 0; iPoint < nNewPoints; iPoint++){
        Double_t yNewValueErrlow    = grint->GetEY()[iPoint]; //Getting the confidencial error from the fit
        Double_t yNewValueErrhigh   = grint->GetEY()[iPoint]; //Getting the confidencial error from the fit
        (*newSpectrum)->SetPoint(iPoint,xNewValue[iPoint],yNewValue[iPoint]);
        (*newSpectrum)->SetPointError(iPoint,xNewValueErrlow[iPoint],xNewValueErrhigh[iPoint],yNewValueErrlow,yNewValueErrhigh);
        cout<<iPoint<<"\t"<<xNewValue[iPoint]<<"\t"<<yNewValue[iPoint]<<"\t"<<xNewValueErrlow[iPoint]<<"\t"<<xNewValueErrhigh[iPoint]<<"\t"<<yNewValueErrlow<<"\t"<<yNewValueErrhigh<<"\t"<<yNewValueRelErrlow[iPoint]<<endl;
    }
    return CurrentFit;
    
}

TGraphAsymmErrors* BinYShiftwithErrhigh(TGraphAsymmErrors* spectrum,TGraphAsymmErrors* spectrumErrhigh){
  
    TGraphAsymmErrors* spectrumClone = (TGraphAsymmErrors*)spectrum->Clone();
     
    const Int_t nNewPoints         = spectrumClone->GetN();
    Double_t *xValue         = spectrumClone->GetX();
    Double_t *yValue	     = spectrumClone->GetY();
    Double_t *xValueErrlow   = spectrumClone->GetEXlow();
    Double_t *xValueErrhigh  = spectrumClone->GetEXhigh();
    Double_t *yValueErrlow   = spectrumClone->GetEYlow();
    Double_t *yValueErrhigh  = spectrumClone->GetEYhigh();
    Double_t yNewValue[nNewPoints];
    
    //Get errors Y high to shift the spectrum
    Double_t *yValueErrShift = spectrumErrhigh->GetEYhigh();
    Double_t *yValueShift    = spectrumErrhigh->GetY();
    
    for(Int_t iPoint = 0; iPoint < nNewPoints; iPoint++){
        
      yNewValue[iPoint] = yValue[iPoint] + TMath::Abs((yValueErrShift[iPoint]/yValueShift[iPoint])*yValue[iPoint]);
     
      yValueErrlow[iPoint]  = (yValueErrlow[iPoint]/yValue[iPoint])*yNewValue[iPoint];
      yValueErrhigh[iPoint] = (yValueErrhigh[iPoint]/yValue[iPoint])*yNewValue[iPoint];
  
    
    }
    
    TGraphAsymmErrors* newSpectrumYShifted = new TGraphAsymmErrors(nNewPoints,xValue,yNewValue,xValueErrlow,xValueErrhigh,yValueErrlow,yValueErrhigh);
    
    return newSpectrumYShifted;
  
    
}


TGraphAsymmErrors* BinYShiftwithErrlow(TGraphAsymmErrors* spectrum,TGraphAsymmErrors* spectrumErrlow){
  
    cout<<"Cross checking Bin Shift"<<endl;
    
    spectrumErrlow->Print();
    
    
  
  
    TGraphAsymmErrors* spectrumClone = (TGraphAsymmErrors*)spectrum->Clone();
     
    const Int_t nNewPoints         = spectrumClone->GetN();
    Double_t *xValue         = spectrumClone->GetX();
    Double_t *yValue	     = spectrumClone->GetY();
    Double_t *xValueErrlow   = spectrumClone->GetEXlow();
    Double_t *xValueErrhigh  = spectrumClone->GetEXhigh();
    Double_t *yValueErrlow   = spectrumClone->GetEYlow();
    Double_t *yValueErrhigh  = spectrumClone->GetEYhigh();
    Double_t *yValueErrShift = spectrumErrlow->GetEYlow();
    Double_t *yValueShift    = spectrumErrlow->GetY();
    
    //Get errors Y high to shift the spectrum
    Double_t yNewValue[nNewPoints];
    
    for(Int_t iPoint = 0; iPoint < nNewPoints; iPoint++){
        
      yNewValue[iPoint] = yValue[iPoint] - TMath::Abs((yValueErrShift[iPoint]/yValueShift[iPoint])*yValue[iPoint]);
      yValueErrlow[iPoint]  = (yValueErrlow[iPoint]/yValue[iPoint])*yNewValue[iPoint];
      yValueErrhigh[iPoint] = (yValueErrhigh[iPoint]/yValue[iPoint])*yNewValue[iPoint];
  
    
    }
    
    TGraphAsymmErrors* newSpectrumYShifted = new TGraphAsymmErrors(nNewPoints,xValue,yNewValue,xValueErrlow,xValueErrhigh,yValueErrlow,yValueErrhigh);
    
    cout<<"///"<<endl;
    newSpectrumYShifted->Print();
    
    return newSpectrumYShifted;
  
    
}

TGraphAsymmErrors* CalculateSystErrors(TGraphErrors* spectrum,TGraphAsymmErrors* g1,TGraphAsymmErrors* g2){

    Int_t nPoints               = spectrum->GetN();
    Double_t* xBins             = spectrum->GetX();
    Double_t* yBins             = spectrum->GetY();
    Double_t *g1yErrlow         = g1->GetEYlow();
    Double_t *g1xErrlow         = g1->GetEXlow();
    Double_t *g1xErrhigh        = g1->GetEXhigh();
    Double_t *g1yVal            = g1->GetY();
    Double_t *g2yErrlow         = g2->GetEYlow();
    Double_t *g2yVal            = g2->GetY();
    TGraphAsymmErrors* graphErr = new TGraphAsymmErrors(nPoints);
    cout<<"iPoin\tx\ty\tXerrLow\tXerrHigh\tYerrLow\tYerrHigh"<<endl;;
    
    for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
        Double_t g1yrelErr  = g1yErrlow[iPoint]/g1yVal[iPoint];
        Double_t g2yrelErr  = g2yErrlow[iPoint]/g2yVal[iPoint];
        Double_t ySysErr    = 0.0;
        
        if(  g1yrelErr >= g2yrelErr ){
            ySysErr         = yBins[iPoint]*g1yrelErr;
        } else {
            ySysErr         = yBins[iPoint]*g2yrelErr;
        }
        graphErr->SetPoint(iPoint,xBins[iPoint],yBins[iPoint]);
        graphErr->SetPointError(iPoint,g1xErrlow[iPoint],g1xErrhigh[iPoint],ySysErr,ySysErr);
        cout<<iPoint<<"\t"<<xBins[iPoint]<<"\t"<<yBins[iPoint]<<"\t"<<g1xErrlow[iPoint]<<"\t"<<g1xErrhigh[iPoint]<<"\t"<<ySysErr<<"\t"<<ySysErr<<endl;
    }
    return graphErr;
}


TGraphAsymmErrors* CalculateSystErrors(TGraphErrors* spectrum,TGraphAsymmErrors* g1,TGraphAsymmErrors* g2,TGraphAsymmErrors* g3){

    Int_t nPoints               = spectrum->GetN();
    Double_t* xBins             = spectrum->GetX();
    Double_t* yBins             = spectrum->GetY();
    Double_t *g1yErrlow         = g1->GetEYlow();
    Double_t *g1xErrlow         = g1->GetEXlow();
    Double_t *g1xErrhigh        = g1->GetEXhigh();
    Double_t *g1yVal            = g1->GetY();
    Double_t *g2yErrlow         = g2->GetEYlow();
    Double_t *g2yVal            = g2->GetY();
    
    Double_t *g3yErrlow         = g3->GetEYlow();
    Double_t *g3yVal            = g3->GetY();
    
    
    
    
    
    TGraphAsymmErrors* graphErr = new TGraphAsymmErrors(nPoints);
    cout<<"iPoin\tx\ty\tXerrLow\tXerrHigh\tYerrLow\tYerrHigh"<<endl;;
    
    for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
        Double_t g1yrelErr  = g1yErrlow[iPoint]/g1yVal[iPoint];
        Double_t g2yrelErr  = g2yErrlow[iPoint]/g2yVal[iPoint];
        Double_t g3yrelErr  = g3yErrlow[iPoint]/g3yVal[iPoint];
        Double_t ySysErr    = 0.0;
        
        if(  g1yrelErr >= g2yrelErr ){
            
            if( g1yrelErr >= g3yrelErr ) {
            ySysErr         = yBins[iPoint]*g1yrelErr;
            } else {
                
            ySysErr         = yBins[iPoint]*g3yrelErr;
            
            }
        } else if( g2yrelErr >= g3yrelErr ) {
            ySysErr         = yBins[iPoint]*g2yrelErr;
        } else {
            
             ySysErr         = yBins[iPoint]*g3yrelErr;
        }
            
        graphErr->SetPoint(iPoint,xBins[iPoint],yBins[iPoint]);
        graphErr->SetPointError(iPoint,g1xErrlow[iPoint],g1xErrhigh[iPoint],ySysErr,ySysErr);
        cout<<iPoint<<"\t"<<xBins[iPoint]<<"\t"<<yBins[iPoint]<<"\t"<<g1xErrlow[iPoint]<<"\t"<<g1xErrhigh[iPoint]<<"\t"<<ySysErr<<"\t"<<ySysErr<<endl;
    }
    return graphErr;
}




TGraphAsymmErrors* CalculateSystErrors(TGraphAsymmErrors* spectrum,TGraphErrors* g1,TGraphErrors* g2){

    Int_t nPoints               = spectrum->GetN();
    Double_t* xBins             = spectrum->GetX();
    Double_t* yBins             = spectrum->GetY();
    Double_t *g1xErrlow         = spectrum->GetEXlow();
    Double_t *g1xErrhigh        = spectrum->GetEXhigh();
    
    
  //  Double_t *g1yErr         = g1->GetEY();
    Double_t *g1yVal         = g1->GetY();
    //Double_t *g2yErr         = g2->GetEY();
    Double_t *g2yVal         = g2->GetY();
    
    TGraphAsymmErrors* graphErr = new TGraphAsymmErrors(nPoints);
    cout<<"iPoin\tx\ty\tXerrLow\tXerrHigh\tYerrLow\tYerrHigh"<<endl;;
    
    for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
        
        Double_t g1yrelErr  = 0;
	Double_t g2yrelErr  = 0;
        
	g1yrelErr  = TMath::Abs( (yBins[iPoint]-g1yVal[iPoint])/yBins[iPoint] );
	g2yrelErr  = TMath::Abs( (yBins[iPoint]-g2yVal[iPoint])/yBins[iPoint] );
	
        Double_t ySysErr    = 0.0;
        
        if(  g1yrelErr >= g2yrelErr ){
             ySysErr         = yBins[iPoint]*g1yrelErr;
        } else {
             ySysErr         = yBins[iPoint]*g2yrelErr;
        }
        graphErr->SetPoint(iPoint,xBins[iPoint],yBins[iPoint]);
        graphErr->SetPointError(iPoint,g1xErrlow[iPoint],g1xErrhigh[iPoint],ySysErr,ySysErr);
        cout<<iPoint<<"\t"<<xBins[iPoint]<<"\t"<<yBins[iPoint]<<"\t"<<g1xErrlow[iPoint]<<"\t"<<g1xErrhigh[iPoint]<<"\t"<<ySysErr<<"\t"<<ySysErr<<endl;
    }
    return graphErr;
}



void CalcRpPb(  TGraphAsymmErrors* PPSpectrumSystErr, 
                TGraphAsymmErrors* PPSpectrumStatErr, 
                TGraphAsymmErrors* pPbSpectrumSystErr, 
                TGraphAsymmErrors* pPbSpectrumStatErr,
                TGraphAsymmErrors** graphRpPbSystErr, 
                TGraphAsymmErrors** graphRpPbStatErr, 
                Bool_t Individual
             ){
    
    //Computing RpPb using the interpolated pp@5.023 TeV and PCM pp@7 TeV spectra
    //Double_t        xSectionpPb5023GeVINEL   =  70*1e-3;  //Dangerous to should be put in a library
        
    Int_t nPoints               = pPbSpectrumStatErr->GetN();

    (*graphRpPbStatErr)         = new TGraphAsymmErrors( nPoints);
    (*graphRpPbSystErr)         = new TGraphAsymmErrors( nPoints);
        
    Double_t *xBins             =  pPbSpectrumStatErr->GetX();
    Double_t *xErrlow           =  pPbSpectrumStatErr->GetEXlow();
    Double_t *xErrhigh          =  pPbSpectrumStatErr->GetEXhigh();
    Double_t *ypPbBins          =  pPbSpectrumStatErr->GetY();
    Double_t *yPPBins           =  PPSpectrumStatErr->GetY();
    Double_t *ypPbStatErrlow    = pPbSpectrumStatErr->GetEYlow();
    //Double_t *ypPbStatErrhigh = pPbSpectrumStatErr->GetEYlow();
    Double_t *ypPbSystErrlow    = pPbSpectrumSystErr->GetEYlow();
    Double_t *ypPbSystErrhigh   = pPbSpectrumSystErr->GetEYlow();
    Double_t *yPPStatErrlow     = PPSpectrumStatErr->GetEYlow();
    //Double_t *yPPStatErrhigh = PPSpectrumStatErr->GetEYlow();
    Double_t *yPPSystErrlow     = PPSpectrumSystErr->GetEYlow();
    Double_t *yPPSystErrhigh    = PPSpectrumSystErr->GetEYlow();
    Double_t *RpPb              = new Double_t[nPoints];
    Double_t fNcoll             = 6.9;
    Double_t fTpPb              = 0.0983e3*(1/recalcBarn);
    Double_t fTpPbErr           = 0.0035e3*(1/recalcBarn);
    
    for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
        RpPb[iPoint]            = ypPbBins[iPoint] / (fTpPb* yPPBins[iPoint]);
        
        (*graphRpPbSystErr)->SetPoint( iPoint, xBins[iPoint], RpPb[iPoint]);
        (*graphRpPbStatErr)->SetPoint( iPoint, xBins[iPoint], RpPb[iPoint] );
    
        Double_t errYStat       = pow( pow( ypPbStatErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( yPPStatErrlow[iPoint]/yPPBins[iPoint],  2.), 0.5)* RpPb[iPoint];
        Double_t errYSystlow;
        Double_t errYSysthigh; 
        if (Individual) {//assigne TpPbErr only to combined result
            cout<<"Test Individual!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<< endl;
            errYSystlow     = pow( pow( ypPbSystErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( yPPSystErrlow[iPoint]/yPPBins[iPoint]  ,2. ) ,   0.5) * RpPb[iPoint];
            errYSysthigh    = pow( pow( ypPbSystErrhigh[iPoint]/ypPbBins[iPoint],2. ) + pow( yPPSystErrhigh[iPoint]/yPPBins[iPoint], 2. ) ,   0.5) * RpPb[iPoint];    
        } else{	   cout<<"Test Comb!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<< endl;
            errYSystlow     = pow( pow( ypPbSystErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( yPPSystErrlow[iPoint]/yPPBins[iPoint]  ,2. ) + pow ((fTpPbErr/fTpPb) , 2),   0.5) * RpPb[iPoint];
            errYSysthigh    = pow( pow( ypPbSystErrhigh[iPoint]/ypPbBins[iPoint],2. ) + pow( yPPSystErrhigh[iPoint]/yPPBins[iPoint], 2. ) + pow ((fTpPbErr/fTpPb) , 2),   0.5) * RpPb[iPoint];
        }
        
        (*graphRpPbStatErr)->SetPointError(iPoint, xErrlow[iPoint],xErrhigh[iPoint],errYStat,errYStat);
        (*graphRpPbSystErr)->SetPointError(iPoint, xErrlow[iPoint],xErrhigh[iPoint],errYSystlow,errYSysthigh);
    }
}

void CalcRpPbInvYield(  TGraphAsymmErrors* PPSpectrumSystErr, 
                        TGraphAsymmErrors*  PPSpectrumStatErr, 
                        TGraphAsymmErrors* pPbSpectrumSystErr, 
                        TGraphAsymmErrors* pPbSpectrumStatErr,
                        TGraphAsymmErrors** graphRpPbSystErr, 
                        TGraphAsymmErrors** graphRpPbStatErr
                     ){
    //Computing RpPb using the interpolated pp@5.023 TeV and PCM pp@7 TeV spectra
    //Double_t        xSectionpPb5023GeVINEL   =  70*1e-3;  //Dangerous to should be put in a library
    
    Int_t nPoints               = pPbSpectrumStatErr->GetN();    
    (*graphRpPbStatErr)         = new TGraphAsymmErrors( nPoints);
    (*graphRpPbSystErr)         = new TGraphAsymmErrors( nPoints);
    
    Double_t *xBins             = pPbSpectrumStatErr->GetX();
    Double_t *xErrlow           = pPbSpectrumStatErr->GetEXlow();
    Double_t *xErrhigh          = pPbSpectrumStatErr->GetEXhigh();
    Double_t *ypPbBins          = pPbSpectrumStatErr->GetY();
    Double_t *yPPBins           = PPSpectrumStatErr->GetY();
    Double_t *ypPbStatErrlow    = pPbSpectrumStatErr->GetEYlow();
    //Double_t *ypPbStatErrhigh = pPbSpectrumStatErr->GetEYlow();
    Double_t *ypPbSystErrlow    = pPbSpectrumSystErr->GetEYlow();
    Double_t *ypPbSystErrhigh   = pPbSpectrumSystErr->GetEYlow();
    Double_t *yPPStatErrlow     = PPSpectrumStatErr->GetEYlow();
    //Double_t *yPPStatErrhigh  = PPSpectrumStatErr->GetEYlow();
    Double_t *yPPSystErrlow     = PPSpectrumSystErr->GetEYlow();
    Double_t *yPPSystErrhigh    = PPSpectrumSystErr->GetEYlow();
    Double_t *RpPb   = new Double_t[nPoints];
    Double_t fNcoll             = 6.9;
    Double_t fNcollErr          = 0.7;
    
    for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
        RpPb[iPoint]            = ypPbBins[iPoint] / (fNcoll* yPPBins[iPoint]);        
        (*graphRpPbSystErr)->SetPoint( iPoint, xBins[iPoint], RpPb[iPoint]);
        (*graphRpPbStatErr)->SetPoint( iPoint, xBins[iPoint], RpPb[iPoint] );
        
        Double_t errYStat       = pow( pow( ypPbStatErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( yPPStatErrlow[iPoint]/yPPBins[iPoint],  2.), 0.5)* RpPb[iPoint];
        Double_t errYSystlow    = pow( pow( ypPbSystErrlow[iPoint]/ypPbBins[iPoint], 2. ) + pow( yPPSystErrlow[iPoint]/yPPBins[iPoint]  ,2. ) + pow ((fNcollErr/fNcoll) , 2),   0.5) * RpPb[iPoint];
        Double_t errYSysthigh   = pow( pow( ypPbSystErrhigh[iPoint]/ypPbBins[iPoint],2. ) + pow( yPPSystErrhigh[iPoint]/yPPBins[iPoint], 2. ) + pow ((fNcollErr/fNcoll) , 2),   0.5) * RpPb[iPoint];
        
        (*graphRpPbStatErr)->SetPointError(iPoint, xErrlow[iPoint],xErrhigh[iPoint],errYStat,errYStat);
        (*graphRpPbSystErr)->SetPointError(iPoint, xErrlow[iPoint],xErrhigh[iPoint],errYSystlow,errYSysthigh);
    }
}

TGraphAsymmErrors* ConvertTGraphErrorstoTGraphAsymmErrors(TGraphErrors* g1){
    
    Int_t nPoints = g1->GetN();
    Double_t *xValue = g1->GetX();
    Double_t *yValue = g1->GetY();
    TGraphAsymmErrors* gNew = new TGraphAsymmErrors(nPoints);
  
    for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
        Double_t xError,yError;
        xError = g1->GetErrorX(iPoint);
        yError = g1->GetErrorY(iPoint);
        gNew->SetPoint(iPoint,xValue[iPoint],yValue[iPoint]);
        gNew->SetPointError(iPoint,xError,xError,yError,yError);
    }
    return gNew;
}


TGraphErrors *ConvertTGraphAsymmErrorstoTGraphErrors(TGraphAsymmErrors* g1){

    Int_t nPoints = g1->GetN();
    Double_t *xValue = g1->GetX();
    Double_t *yValue = g1->GetY();
    TGraphErrors* gNew = new TGraphErrors(nPoints);

    for(Int_t iPoint = 0; iPoint < nPoints; iPoint++){
        Double_t xError,yError;
        xError = g1->GetErrorX(iPoint);
        yError = g1->GetErrorY(iPoint);
        gNew->SetPoint(iPoint,xValue[iPoint],yValue[iPoint]);
        gNew->SetPointError(iPoint,xError,yError);
    }
    return gNew;
}

TGraphErrors *GetInterpolSpectrum2D(TGraphErrors *g1, TGraphErrors *g2, Double_t d1, Double_t d2,Double_t dSqrts,TString opc){
         
    if(!g1) return 0x0;
    if(!g2) return 0x0;

    TGraphErrors  *gInterpol     = new TGraphErrors(g1->GetN());
    TGraphErrors  *gAlpha        = new TGraphErrors(g1->GetN());
    TGraphErrors** gPtvsSqrts    = new TGraphErrors*[g1->GetN()];
    TGraphErrors** gPtvsEnergies = new TGraphErrors*[g1->GetN()];
    TF1**          fPowerlawFits = new TF1*[g1->GetN()];

    for(Int_t i = 0; i < g1->GetN(); i++){
        TGraphErrors *grint = new TGraphErrors(1);
        grint->SetPoint(0, dSqrts, 0);
        TGraphErrors *gToFit = new TGraphErrors(2);

        gToFit->SetPoint(0, d1, g1->GetY()[i]);
        gToFit->SetPointError(0, 0, g1->GetEY()[i]);
        gToFit->SetPoint(1, d2, g2->GetY()[i]);
        gToFit->SetPointError(1, 0, g2->GetEY()[i]);

        TF1 *fPowerlaw = new TF1("fPowerlaw","[0]*x^([1])", 0,10000);
        fPowerlaw->SetParameters(0, 0.1);
        fPowerlaw->SetParameters(1, 2.0);

        for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlaw,"Q");
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
        
        Double_t alpha  = fPowerlaw->GetParameter(1);
        Double_t alphaE = fPowerlaw->GetParError(1);
        
        gInterpol->SetPoint(i, g1->GetX()[i],fPowerlaw->Eval(dSqrts));
        gInterpol->SetPointError(i, 0, grint->GetEY()[0]);
        
        gAlpha->SetPoint(i, g1->GetX()[i],alpha);
        //gAlpha->SetPointError(i, 0,alphaE);
        
        gPtvsSqrts[i]= new TGraphErrors(1);
        gPtvsSqrts[i]->SetPoint(0,dSqrts,gInterpol->GetY()[i]);
        //gPtvsSqrts[i]->Sort();
        
        fPowerlawFits[i] = fPowerlaw;
        gPtvsEnergies[i] = gToFit;
        delete grint;
        //delete fPowerlaw;
    }
    
    if( opc.EqualTo("SystUp") || opc.EqualTo("SystDown") || opc.EqualTo("Stat") ) {
      
        Int_t index = 0;
        
        if( opc.EqualTo("Stat") ) index = 0;
	else if (opc.EqualTo("SystDown") ) index = 1;
	else if (opc.EqualTo("SystUp") ) index = 2;
       
        graphAlpha[index]    = gAlpha;
        graphPtvsSqrts[index]      = new TGraphErrors*[g1->GetN()];
        fPowerlawSystem[index]     = new TF1*[g1->GetN()];
        gPtvsEnergiesSystem[index] = new TGraphErrors*[g1->GetN()];
        for ( Int_t i = 0; i < g1->GetN(); i++ ){
            graphPtvsSqrts[index][i]       = gPtvsSqrts[i];
            fPowerlawSystem[index][i]      = fPowerlawFits[i];
            gPtvsEnergiesSystem[index][i]  = gPtvsEnergies[i];
        }
    }
    return gInterpol;
}

TGraphErrors *GetInterpolSpectrum3D(TGraphErrors *g1, TGraphErrors *g2,TGraphErrors *g3, Double_t d1, Double_t d2, Double_t d3, Double_t dSqrts,TString opc){
    if(!g1) return 0x0;
    if(!g2) return 0x0;
    if(!g3) return 0x0;
    
    TGraphErrors  *gInterpol     = new TGraphErrors(g1->GetN());
    TGraphErrors  *gAlpha        = new TGraphErrors(g1->GetN());
    TGraphErrors** gPtvsSqrts    = new TGraphErrors*[g1->GetN()];
    TGraphErrors** gPtvsEnergies = new TGraphErrors*[g1->GetN()];
    TF1**          fPowerlawFits = new TF1*[g1->GetN()];
    

   // TGraphErrors *gInterpol = new TGraphErrors(g1->GetN());

    for(Int_t i = 0; i < g1->GetN(); i++){
        TGraphErrors *grint = new TGraphErrors(1);
        grint->SetPoint(0, dSqrts, 0);
        TGraphErrors *gToFit = new TGraphErrors(3);
        gToFit->SetPoint(0, d1, g1->GetY()[i]);
        gToFit->SetPointError(0, 0, g1->GetEY()[i]);
        gToFit->SetPoint(1, d2, g2->GetY()[i]);
        gToFit->SetPointError(1, 0, g2->GetEY()[i]);
        gToFit->SetPoint(2, d3, g3->GetY()[i]);
        gToFit->SetPointError(2, 0, g3->GetEY()[i]);

        TF1 *fPowerlaw = new TF1("fPowerlaw","[0]*x^([1])", 0,10000);
        fPowerlaw->SetParameters(0, 0.1);
        fPowerlaw->SetParameters(1, 2.0);

        for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlaw,"Q");
        
         
        Double_t alpha  = fPowerlaw->GetParameter(1);
        Double_t alphaE = fPowerlaw->GetParError(1);

        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.48);
        gInterpol->SetPoint(i, g1->GetX()[i],
        fPowerlaw->Eval(dSqrts));
        gInterpol->SetPointError(i, 0, grint->GetEY()[0]);
       
                
        gAlpha->SetPoint(i, g1->GetX()[i],alpha);
        //gAlpha->SetPointError(i, 0,alphaE);
        
        gPtvsSqrts[i]= new TGraphErrors(1);
        gPtvsSqrts[i]->SetPoint(0,dSqrts,gInterpol->GetY()[i]);
        //gPtvsSqrts[i]->Sort();
        
        fPowerlawFits[i] = fPowerlaw;
        gPtvsEnergies[i] = gToFit;
        
        
        

        delete grint;
        delete fPowerlaw;
    }
    
     if( opc.EqualTo("SystUp") || opc.EqualTo("SystDown") || opc.EqualTo("Stat") ) {
      
        Int_t index = 0;
        
        if( opc.EqualTo("Stat") ) index = 0;
	else if (opc.EqualTo("SystDown") ) index = 1;
	else if (opc.EqualTo("SystUp") ) index = 2;
       
        graphAlpha[index]    = gAlpha;
        graphPtvsSqrts[index]      = new TGraphErrors*[g1->GetN()];
        fPowerlawSystem[index]     = new TF1*[g1->GetN()];
        gPtvsEnergiesSystem[index] = new TGraphErrors*[g1->GetN()];
        for ( Int_t i = 0; i < g1->GetN(); i++ ){
            graphPtvsSqrts[index][i]       = gPtvsSqrts[i];
            fPowerlawSystem[index][i]      = fPowerlawFits[i];
            gPtvsEnergiesSystem[index][i]  = gPtvsEnergies[i];
        }
    }
    
    
    

    return gInterpol;
}


TGraphAsymmErrors* CancelOutMaterialError(TGraphAsymmErrors* graphYieldPi0, TString method, TGraphAsymmErrors* GraphErrCancellation){

    TGraphAsymmErrors* graphYieldPi0Clone   = (TGraphAsymmErrors*) graphYieldPi0->Clone();
    Double_t* valueX                        = graphYieldPi0Clone->GetX();
    Double_t* valueY                        = graphYieldPi0Clone->GetY();
    Double_t* errorXlow                     = graphYieldPi0Clone->GetEXlow();
    Double_t* errorXhigh                    = graphYieldPi0Clone->GetEXhigh();
    Double_t* errorYlow                     = graphYieldPi0Clone->GetEYlow();
    Double_t* errorYhigh                    = graphYieldPi0Clone->GetEYhigh();
    Int_t nPoints                           = graphYieldPi0Clone->GetN();

    if( method.CompareTo("PCMPCM") == 0 ){  //For RpPb PCM/PCM cancelout the error of the material    
        Double_t  errorMaterial = 4.5;
        for(Int_t i = 0; i < nPoints; i++){
            Double_t materialTwoE   =  ( 2*errorMaterial * valueY[i] ) / 100.0;
            //cout<<"Before "<<i<<" "<<valueY[i]<<" "<<errorYlow[i]<<" errorYLow/valueY "<<errorYlow[i]/valueY[i]<<endl;
            errorYlow[i]  = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( materialTwoE * materialTwoE ) );
            //cout<<"After  "<<i<<" "<<valueY[i]<<" "<<errorYlow[i]<<" errorYLow/valueY "<<errorYlow[i]/valueY[i]<<endl;
            //cout<<"Before "<<i<<" "<<valueY[i]<<" "<<errorYhigh[i]<<" errorYhigh/valueY "<<errorYhigh[i]/valueY[i]<<endl;
            errorYhigh[i] = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( materialTwoE * materialTwoE ) );
            //cout<<"After  "<<i<<" "<<valueY[i]<<" "<<errorYhigh[i]<<" errorYhigh/valueY "<<errorYhigh[i]/valueY[i]<<endl;
        }
    } else if( method.CompareTo("PCMDalitz") == 0) {
        Double_t  errorMaterial = 4.5;
        for(Int_t i = 0; i < nPoints; i++){
            Double_t materialE      = (  errorMaterial   * valueY[i] ) / 100.0;
            Double_t materialTwoE   = (  2*errorMaterial * valueY[i] ) / 100.0;
            errorYlow[i]            = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( materialTwoE * materialTwoE ) + ( materialE*materialE ) );
            errorYhigh[i]           = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( materialTwoE * materialTwoE ) + ( materialE*materialE ) );
        }
    } else if ( method.CompareTo("Dalitz") == 0 ) {
        Double_t  errorMaterial = 4.5;
        for(Int_t i = 0; i < nPoints; i++){
            Double_t materialE 	=  (  errorMaterial   * valueY[i] ) / 100.0;
            //cout<<valueX[i]<<" "<<"errorLow:  "<<errorYlow[i] * errorYlow[i]<<"  material:    "<<materialE*materialE<<endl;
            //cout<<valueX[i]<<" "<<"errorHigh: "<<errorYhigh[i] * errorYhigh[i]<<"  material:  "<<materialE*materialE<<endl;
            errorYlow[i]  = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( materialE*materialE ) );
            errorYhigh[i] = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( materialE*materialE ) );
        }
    } else if ( method.CompareTo("PHOSPHOS") == 0 ) {
        cout<<"Entro PHOS"<<endl;
        TGraphAsymmErrors* graphErrCancellationClone = (TGraphAsymmErrors*) GraphErrCancellation->Clone();	
        Double_t* valueSystErrCancellation     = graphErrCancellationClone->GetY();	
        Double_t  errorMaterial = 0.5;
        for(Int_t i = 0; i < nPoints; i++){
            errorMaterial=valueSystErrCancellation[i]*100.;
            cout<<"Cancel out: "<< errorMaterial << endl;  
            Double_t materialE 	=  (  errorMaterial   * valueY[i] ) / 100.0;
            cout<<valueX[i]<<" "<<"errorLow:  "<< errorYlow[i]/ valueY[i]*100<<"  material:    "<<errorMaterial<<endl;
            cout<<valueX[i]<<" "<<"errorHigh: "<< errorYhigh[i]/ valueY[i]*100<<"  material:  "<<errorMaterial<<endl;
            errorYlow[i]  = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( materialE*materialE ) );
            errorYhigh[i] = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( materialE*materialE ) );
            cout<<valueX[i]<<" "<<"errorLow:  "<< errorYlow[i]/ valueY[i]*100<<endl;
            cout<<valueX[i]<<" "<<"errorHigh: "<< errorYhigh[i]/ valueY[i]*100<<endl;
        }
    } else if (method.CompareTo("PHOSPHOSv2") == 0 ){  
        
        cout<<"Entro PHOSv2"<<endl;
        
        Double_t  errorMaterial = 3.5;
        for(Int_t i = 0; i < nPoints; i++){
            Double_t materialE 	=  (  errorMaterial   * valueY[i] ) / 100.0;
            errorYlow[i]  = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( materialE*materialE ) );
            errorYhigh[i] = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( materialE*materialE ) );
        }
    } else if ( method.CompareTo("Comb") == 0 ) {
      cout<<"Entro Comb"<<endl;
      Double_t  errorMaterial = 4.5;
        for(Int_t i = 0; i < nPoints; i++){
            Double_t materialE      = (  errorMaterial   * valueY[i] ) / 100.0;
            Double_t materialTwoE   = (  2*errorMaterial * valueY[i] ) / 100.0;
            errorYlow[i]            = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( materialTwoE * materialTwoE ) + ( materialE*materialE ) );
            errorYhigh[i]           = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( materialTwoE * materialTwoE ) + ( materialE*materialE ) );
        }
      
    } else if ( method.CompareTo("PCM-EMCal") == 0 || method.CompareTo("PCM-EMCAL") == 0){
      
        Double_t  errorMaterial = 4.5;
      
        for(Int_t i = 0; i < nPoints; i++){
            Double_t materialE 	=  ( errorMaterial   * valueY[i] ) / 100.0;
            errorYlow[i]  = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( materialE*materialE ) );
            errorYhigh[i] = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( materialE*materialE ) );
        }
      
    } else if ( method.CompareTo("EMCEMC") == 0 ){
    
        Double_t  errorClusterMaterialTRD = 4.24;
      
        for(Int_t i = 0; i < nPoints; i++){
            
            Double_t errorClsMatTRD	=  ( errorClusterMaterialTRD   * valueY[i] ) / 100.0;
            errorYlow[i]  = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( errorClsMatTRD*errorClsMatTRD ) );
            errorYhigh[i] = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( errorClsMatTRD*errorClsMatTRD ) );
        }
        
    } else if ( method.CompareTo("PCMEMC") == 0 ) {
        
        Double_t  errorClusterMaterialTRD = 4.24;
        Double_t  errorInnerMaterial      = 4.5;
        
        for(Int_t i = 0; i < nPoints; i++){
            
            Double_t errorClsMatTRD	=  ( errorClusterMaterialTRD   * valueY[i] ) / 100.0;
            Double_t errorInnerMat      =  ( errorInnerMaterial        * valueY[i] ) / 100.0;
            
            errorYlow[i]  = TMath::Sqrt(   ( errorYlow[i] * errorYlow[i]  )  - ( errorClsMatTRD*errorClsMatTRD ) - (errorInnerMat*errorInnerMat) );
            errorYhigh[i] = TMath::Sqrt(   ( errorYhigh[i]* errorYhigh[i] )  - ( errorClsMatTRD*errorClsMatTRD ) - (errorInnerMat*errorInnerMat) );
        }
        
        
    }    
    TGraphAsymmErrors* graphInvYieldPi0CancelOutMaterial = new TGraphAsymmErrors(nPoints,valueX,valueY,errorXlow,errorXhigh,errorYlow,errorYhigh);
    return graphInvYieldPi0CancelOutMaterial;
}

void GetTGraphErrorsUpDown(TGraphAsymmErrors* graphYieldError, TGraphErrors** graphRelErrUpDown){
    
    Double_t* RelErrUp      = ExtractRelErrUpAsymmGraphClone(graphYieldError);
    Int_t nPoints           = graphYieldError->GetN();
    Double_t *xValue        = graphYieldError->GetX();
    (*graphRelErrUpDown)    = new TGraphErrors(nPoints);

    for(Int_t iPoint = 0; iPoint < nPoints; iPoint++ ){
        (*graphRelErrUpDown)->SetPoint( iPoint, xValue[iPoint], RelErrUp[iPoint] );
        (*graphRelErrUpDown)->SetPointError(iPoint, 0, 0);
    }
}

Double_t* ExtractRelErrUpAsymmGraphClone(TGraphAsymmErrors* graph){
    TGraphAsymmErrors* graphClone   = (TGraphAsymmErrors*)graph->Clone();
    Double_t* yValue                = graphClone->GetY();
    Double_t* yErrorHigh            = graphClone->GetEYhigh();
    Double_t* yErrorHigh2           = yErrorHigh;
    Int_t nPoints                   = graphClone->GetN();
    for (Int_t i = 0; i < nPoints; i++){
        yErrorHigh2[i]              = yErrorHigh2[i]/yValue[i]*100;
    }
    return yErrorHigh2;
}

Double_t* ExtractRelErrDownAsymmGraphClone(TGraphAsymmErrors* graph){
    TGraphAsymmErrors* graphClone   = (TGraphAsymmErrors*)graph->Clone();
    Double_t* yValue                = graphClone->GetY();
    Double_t* yErrorLow             = graphClone->GetEYlow();
    Double_t* yErrorLow2            = yErrorLow;
    Int_t nPoints                   = graphClone->GetN();
    for (Int_t i = 0; i < nPoints; i++){
        yErrorLow2[i]               = -yErrorLow2[i]/yValue[i]*100;
    }
    return yErrorLow2;
}

TGraphAsymmErrors* ProduceTGraphAsymmToPlotErrors(TGraphAsymmErrors *graph){

    TGraphAsymmErrors* newGraph     = (TGraphAsymmErrors*)graph->Clone();
    Int_t nNewPoints                = newGraph->GetN();
    Double_t *xNewValue             = newGraph->GetX();
    Double_t *yNewValue             = newGraph->GetY();
    Double_t *yErrlow               = newGraph->GetEYlow();
    Double_t *yErrhigh              = newGraph->GetEYhigh();
    Double_t *xErrlow               = newGraph->GetEXlow();	
    Double_t *xErrhigh              = newGraph->GetEXhigh();	
    
    for(Int_t iNewPoint = 0 ; iNewPoint < nNewPoints; iNewPoint++){
        yNewValue[iNewPoint] = yErrlow[iNewPoint]/yNewValue[iNewPoint]*100;
        yErrlow[iNewPoint]      = 0;
        yErrhigh[iNewPoint]     = 0;
        xErrlow[iNewPoint]      = 0;
        xErrhigh[iNewPoint]     = 0;
    }
    return newGraph;    
}
