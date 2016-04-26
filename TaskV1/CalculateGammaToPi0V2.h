#include "../CommonHeaders/ExtractSignalBinning.h"

TString outputDir;

TObjArray *array1;
TLatex *gammaL;
TLatex *nlo;
TLatex *cocktail;
TLatex *tpi;
TLatex *teta;
TLatex *tomega;
TLatex *tetaprime;
TLatex *tphi;
TLatex *trho;

//declaration for printing logo 
Float_t purityDrawingCoordinates[5] = {0.25, 0.76, 0.2, 0.05, 0.0};
Float_t convProbDrawingCoordinates[5] = {0.25, 0.76, 0.2, 0.05, 0.0};
Float_t recoEffDrawingCoordinates[5] =  {0.25, 0.76, 0.2, 0.05, 0.0};
Float_t incRatioDrawingCoordinates[5] =  {0.65, 0.75, 0.2, 0.05, 0.0};
Float_t incRatioAllDrawingCoordinates[5] =  {0.45, 0.75, 0.2, 0.05, 0.0};
Float_t decayRatioDrawingCoordinates[5] =  {0.2, 0.75, 0.2, 0.05, 0.0};


Double_t binningRAA[12] = {0.4,0.8,1.0,1.4,1.8,2.2,2.6,3.0,4.0,5.0,6.0,8.0};
Double_t multiplicity[10] = {1502.7,923.26,558.68,321.20,171.67,85.13,38.51,15.78,6.32,2.63}; 
//                           0-1,   1-2,   2-3,   3-4,   4-5,   5-6,  6-7,  7-8,  8-9, 9-10
Double_t multiplicityCent[6] = {1686.87,1319.89, 1031.9,807.90,627.99,483.95};
//                               0-05,   05-1  

TH1D *histoBinningDoubleRatioRaa;

Bool_t pictDrawingOptions[4] = {kFALSE, kFALSE, kFALSE, kTRUE};

TString textFinder;
TString textPi0New;

TString textPrefix2;
TString Meson_text;

TFile *fileGamma;
TFile *filePi0;
TFile *fileEta;
TFile *filePhi;
TFile *fileQuarkMatter;
TFile *fileSystErrorPbPb;

TString fcocktailFunc = "";
TString fcocktailFuncModA = "";
TString fcocktailFuncModB = "";
TString fcocktailFuncMtScaledEta = "";
TString fcocktailFuncChargedPions = "";
TString fcocktailFuncNarrow = "";
TString fcocktailFuncWide = "";
TString fcocktailFuncSyst = "";
Double_t fcmult;
Bool_t ffive = kFALSE;

TGraphErrors *graphNLOCalcMuHalfcopyA;
TGraphErrors *graphNLOCalcMuHalfcopyB;
TGraphErrors *graphNLOCalcMuHalfcopyC;

TGraphErrors *half;
TGraphErrors *one;
TGraphErrors *two;


Double_t *xHalf;
Double_t *yHalf;
Double_t *eyHalf;

Double_t *xOne;
Double_t *yOne;
Double_t *eyOne;

Double_t *xTwo;
Double_t *yTwo;
Double_t *eyTwo;


TFile *combinedFile;
TFile *combinedFileOld;
TFile *cocktailFile;
TFile *cocktailFileModA;
TFile *cocktailFileModB;
TFile *cocktailFileMtScaledEta;
TFile *cocktailFileHadrons;
TFile *cocktailFileChargedPions;
TFile *cocktailFileWide;
TFile *cocktailFileNarrow;
TFile *cocktailFileSyst;
TDirectoryFile *cocktailDirModA;
TDirectoryFile *cocktailDirModB;
TDirectoryFile *cocktailDir;
TDirectoryFile *cocktailDirMtScaledEta;
TDirectoryFile *cocktailDirHadrons;
TDirectoryFile *cocktailDirChargedPions;
TDirectoryFile *cocktailDirWide;
TDirectoryFile *cocktailDirNarrow;
TDirectoryFile *cocktailDirSyst;

TH1D *histoCorrectedPi0YieldFit;
TH1D *histoCorrectedPi0YieldFitWide;
TH1D *histoCorrectedPi0YieldFitNarrow;


TH1D *histoCorrectedEtaYieldNormalEff;
TH1D *histoCorrectedEtaYield;
TH1D* histoCorrectedEtaYieldPatched;

TH1D *histoCorrectedPhiYieldStat;
TH1D *histoCorrectedPhiYieldSyst;
TF1* fitPhiYield;

const char* fileNameNLOPhotonTwo;
const char* fileNameNLOPhotonOne;
const char* fileNameNLOPhotonHalf;

ifstream inHalf;
ifstream inOne;
ifstream inTwo;
Double_t xSection;
//Double_t recalcBarn = 1e12; //NLO in pbarn!!!!


//Double_t	xSection7TeV =			62.22*1e-3;
//Double_t	xSection7TeVV0AND =			54.31*1e-3;
//Double_t	xSection7TeVErrUp =		2.18;
//Double_t	xSection7TeVErrDown =	2.18;
//Double_t	xSection900GeV =		47.78*1e-3;
//Double_t	xSection900GeVV0AND =		40.06*1e-3;
//Double_t	xSection900GeVErrUp =		 2.39;
//Double_t	xSection900GeVErrDown =		 1.86;
//Double_t	xSection2760GeV = 		55.416*1e-3;
//Double_t	xSection2760GeVV0AND = 		47.73*1e-3;
//Double_t	xSection2760GeVErr = 	3.9;
Double_t xSection2760GeVpp = 		55.416*1e-3;
//Double_t xSection2760GeVErrpp = 	3.9;
Double_t xSection2760GeVppINEL = 62.8*1e9;
Double_t xSection5023GeVppINEL = 		70*1e-3;
Double_t xSection900GeVppINEL = 52.5*1e9;
Double_t xSection7TeVppINEL = 73.2*1e9;	


Int_t nlinesNLOTwo = 0;
Int_t nlinesNLOOne = 0;
Int_t nlinesNLOHalf = 0;
Double_t muXU[100];
Double_t muXD[100];
Double_t ptNLOPhotonTwo[100];
Double_t ptNLOPhotonOne[100];
Double_t ptNLOPhotonHalf[100];
Double_t muHalfDF[100];
Double_t muHalfD[100];
Double_t muHalfF[100];
Double_t muOneDF[100];
Double_t muOneD[100];
Double_t muOneF[100];
Double_t muTwoDF[100];
Double_t muTwoD[100];
Double_t muTwoF[100];

TGraphErrors *graphNLOCalcMuTwo;
TGraphErrors *graphNLOCalcMuOne;
TGraphErrors *graphNLOCalcMuHalf;
TGraphErrors *graphNLOCalcMuTwoRebin;
TGraphErrors *graphNLOCalcMuOneRebin;
TGraphErrors *graphNLOCalcMuHalfRebin;
TH1D *histNLOCalcMuTwoRebin;
TH1D *histNLOCalcMuOneRebin;
TH1D *histNLOCalcMuHalfRebin;

TGraphAsymmErrors* combined7TeV;
TGraphAsymmErrors* combinedEta7TeV;
TGraphAsymmErrors* combined7TeVRebin;

TGraphAsymmErrors* combined7TeVOld;
TGraphAsymmErrors* combinedEta7TeVOld;
TGraphAsymmErrors* combined7TeVRebinOld;

TF1 *One;

TF1 *fitNLOMuHalf;
TF1 *fitNLOMuOne;
TF1 *fitNLOMuTwo;

TF1 *CombinedFitA;
TF1 *CombinedFitB;
TF1 *CombinedFitEtaA;
TF1 *CombinedFitEtaB;

TF1 *CombinedFitHighA;
TF1 *CombinedFitHighEtaA;
TF1 *CombinedFitLowA;
TF1 *CombinedFitLowEtaA;

TF1 *fitPi0YieldA;
TF1 *fitPi0YieldAHighpt;
TF1 *fitPi0YieldAPions;
TF1 *fitPi0YieldAPionsHighpt;
TF1 *fitPi0YieldB;
TF1 *fitEtaYieldA;
TF1 *fitEtaYieldB;
TF1 *fitEtaYieldC;
TF1 *fitPi0YieldC;
TF1 *fitPi0YieldCWide;
TF1 *fitPi0YieldCNarrow;


TF1 *ConversionGammaFitA;
TF1 *ConversionGammaFitB;

TH1D *histRatioNLOMuHalf;
TH1D *histRatioNLOMuOne;
TH1D *histRatioNLOMuTwo;

TH1D *histRatioConversionEtaA;
TH1D *histRatioConversionEtaB;
TH1D* histRatioConversionEtaC;
TH1D *histRatioConversionPi0A;
TH1D *histRatioConversionPi0AHighpt;
TH1D *histRatioConversionPi0B;
TH1D *histRatioConversionPi0C;
TH1D *histRatioConversionGammaA;
TH1D *histRatioConversionGammaB;

TGraphAsymmErrors* graphRatioCombinedEtaA;
TGraphAsymmErrors* graphRatioCombinedEtaB;
TGraphAsymmErrors* graphRatioCombinedPi0A;
TGraphAsymmErrors* graphRatioCombinedPi0B;
TGraphAsymmErrors* histoDirectPhotonNLORAA;

TH1D *histoTruePi0MCData;
TH1D *histoNormalPi0MCData;
TH1D *histoPurityGammaMCData;

TH1D *histoGammaSpecCorrPurity;
TH1D *histoGammaSpecLateCorr;
TH1D *histoGammaSpecMCAll;
TH1D *histoGammaSpecCorrESDMC;
TH1D *histoMCDecaySumGammaPt;
TH1D *histoMCDecayPi0GammaPt;
TH1D *histoMCDecayEtaGammaPt;
TH1D *histoMCDecayEtapGammaPt;
TH1D *histoMCDecayOmegaGammaPt;
TH1D *histoMCDecayRho0GammaPt;
TH1D *histoMCDecayK0sGammaPt;
TH1D *histoMCDecayPhiGammaPt;
TH1D *histMCAllMinusDecay;
TH1D *histoCombinedPi0;
TH1D *histoConversionCombinedRatio;
TH1D *histoMCGammaConvProb;
TH1D *histoGammaPrimaryRecoEffMCPt;
TH1D *histoEventQuality;


TH1D *histoErrorGamma;
TH1D *histoErrorPi0;
TH1D *histoErrorIncRatio;
TH1D *histoErrorCocktailRatio;
TH1D *histoErrorDoubleRatio;

TH1D *histoCorrectedPi0YieldNormalEff;
TH1D *histoCorrectedPi0Yield;
TH1D *histoCorrectedPi0YieldWide;
TH1D *histoCorrectedPi0YieldNarrow;
TH1D *histoMCYieldMeson;
TH1D *histoMCYieldMesonOldBin;

TH1D *histoIncRatioFitPurity;
TH1D *histoIncRatioFitPurityWide;
TH1D *histoIncRatioFitPurityNarrow;
TH1D *histoIncRatioHighFitPurity;
TH1D *histoIncRatioLowFitPurity;
TH1D *histoIncRatioFitFit;
TH1D *histoIncRatioCombinedPurity;
TH1D *histoIncRatioPurity;
TH1D *histoIncRatioPurityTrueEff;
TH1D *histoIncRatioPurityTrueEffWide;
TH1D *histoIncRatioPurityTrueEffNarrow;

TH1D *histoMCIncRatio;
TH1D *histoMCesdIncRatio;
TH1D *histoIncRatioGammaMC;
TH1D *histoIncRatioPurityGammaData;

TH1D *histoDecayRatioSumGamma;
TH1D *histoDecayRatioPi0Gamma;
TH1D *histoDecayRatioEtaGamma;
TH1D *histoDecayRatioEtapGamma;
TH1D *histoDecayRatioOmegaGamma;
TH1D *histoDecayRatioRho0Gamma;
TH1D *histoDecayRatioK0sGamma;
TH1D *histoDecayRatioPhiGamma;

TFile *fileCorrectedOutput;

TH1D *cocktailAllGamma;
TH1D *cocktailPi0Gamma;
TH1D *cocktailEtaGamma;
TH1D *cocktailEtapGamma;
TH1D *cocktailOmegaGamma;
TH1D *cocktailPhiGamma;
TH1D *cocktailRhoGamma;
TH1D *cocktailSigmaGamma;

TH1D *cocktailAllGammaMtScaledEta;
TH1D *cocktailPi0GammaMtScaledEta;
TH1D *cocktailEtaGammaMtScaledEta;
TH1D *cocktailEtapGammaMtScaledEta;
TH1D *cocktailOmegaGammaMtScaledEta;
TH1D *cocktailPhiGammaMtScaledEta;
TH1D *cocktailRhoGammaMtScaledEta;
TH1D *cocktailSigmaGammaMtScaledEta;

TH1D *cocktailAllGammaChargedPions;
TH1D *cocktailPi0GammaChargedPions;
TH1D *cocktailEtaGammaChargedPions;
TH1D *cocktailEtapGammaChargedPions;
TH1D *cocktailOmegaGammaChargedPions;
TH1D *cocktailPhiGammaChargedPions;
TH1D *cocktailRhoGammaChargedPions;
TH1D *cocktailSigmaGammaChargedPions;


TH1D *cocktailPi0;
TH1D *cocktailEta;
TH1D *cocktailPhi;

TH1D *cocktailAllGammaPi0Hadrons;
TH1D *cocktailAllGammaPi0Pions;
TH1D *cocktailAllGammaPi0Wide;
TH1D *cocktailAllGammaPi0Narrow;
TH1D *cocktailAllGammaPi0Syst;
TH1D *cocktailAllGammaPi0;
TH1D *cocktailPi0GammaPi0;
TH1D *cocktailEtaGammaPi0;
TH1D *cocktailEtapGammaPi0;
TH1D *cocktailOmegaGammaPi0;
TH1D *cocktailPhiGammaPi0;
TH1D *cocktailRhoGammaPi0;
TH1D *cocktailSigmaGammaPi0;
TH1D *cocktailEta1;
TH1D *cocktailEta2;

TH1D *cocktailAllGammaPi0ModA;
TH1D *cocktailAllGammaPi0ModB;
TH1D *cocktailAllGammaModA;
TH1D *cocktailAllGammaModB;
TH1D *cocktailPi0ModA;
TH1D *cocktailPi0ModB;

TH1D *cocktailPi0MtScaledEta;
TH1D *cocktailEtaMtScaledEta;
TH1D *cocktailAllGammaPi0MtScaledEta;
TH1D *cocktailPi0GammaPi0MtScaledEta;
TH1D *cocktailEtaGammaPi0MtScaledEta;
TH1D *cocktailEtapGammaPi0MtScaledEta;
TH1D *cocktailOmegaGammaPi0MtScaledEta;
TH1D *cocktailPhiGammaPi0MtScaledEta;
TH1D *cocktailRhoGammaPi0MtScaledEta;
TH1D *cocktailSigmaGammaPi0MtScaledEta;

TH1D *cocktailPi0ChargedPions;
TH1D *cocktailEtaChargedPions;
TH1D *cocktailAllGammaPi0ChargedPions;
TH1D *cocktailPi0GammaPi0ChargedPions;
TH1D *cocktailEtaGammaPi0ChargedPions;
TH1D *cocktailEtapGammaPi0ChargedPions;
TH1D *cocktailOmegaGammaPi0ChargedPions;
TH1D *cocktailPhiGammaPi0ChargedPions;
TH1D *cocktailRhoGammaPi0ChargedPions;
TH1D *cocktailSigmaGammaPi0ChargedPions;

TH1D *cocktailPi0Rebin;
TH1D *cocktailGammaPi0Rebin;

TH1D *histoDoubleRatioConversionNormalPurity;
TH1D *histoDoubleRatioConversionTrueEffPurity;
TH1D *histoDoubleRatioConversionTrueEffPurityModA;
TH1D *histoDoubleRatioConversionTrueEffPurityModB;
TH1D *histoDoubleRatioConversionTrueEffPurityMtScaledEta;
TH1D *histoDoubleRatioConversionTrueEffPurityChargedPions;
TH1D *histoDoubleRatioConversionOnlyGamma;
TH1D *histoDoubleRatioConversionTrueEffPurityHadrons;
TH1D *histoDoubleRatioConversionTrueEffPurityPions;
TH1D *histoDoubleRatioConversionTrueEffPurityWide;
TH1D *histoDoubleRatioConversionTrueEffPurityNarrow;
TH1D *histoDoubleRatioConversionTrueEffPuritySyst;

TH1D *histoDoubleRatioFitPi0YieldPurity;
TH1D *histoDoubleRatioFitPi0YieldPurityModA;
TH1D *histoDoubleRatioFitPi0YieldPurityModB;
TH1D *histoDoubleRatioFitPi0YieldPurityMtScaledEta;
TH1D *histoDoubleRatioFitPi0YieldPurityChargedPions;
TH1D *histoDoubleRatioFitPi0YieldPurityPions;
TH1D *histoDoubleRatioConversionLowFitPurity;
TH1D *histoDoubleRatioConversionHighFitPurity;
TH1D *histoDoubleRatioFitPi0YieldPurityWide;
TH1D *histoDoubleRatioFitPi0YieldPurityNarrow;
TH1D *histoDoubleRatioFitPi0YieldPuritySyst;
TH1D *histoDoubleRatiofitPi0YieldFit;

TH1D *histoDoubleRatioCombinedPurity;
TH1D *histoMCDoubleRatioSum;
TH1D *histoMCesdDoubleRatioSum;

// Systematic Error Graphs
ifstream fileSysErrPi0;

TGraphAsymmErrors *graphCorrectedYieldPi0SysErr;
TGraphAsymmErrors *graphCorrectedYieldPi0SysErrA;

TDirectoryFile *cocktailHighError;
TDirectoryFile *cocktailLowError;

TF1 *FitRandomHigh;
TF1 *FitRandomLow;
TGraphErrors* combined7TeVRandomHighError;
TGraphErrors* combined7TeVRandomLowError;

TF1 *CombinedFitEtaRandomHighA;
TF1 *CombinedFitEtaRandomLowA;
TGraphErrors* combinedEta7TeVRandomHighError;
TGraphErrors* combinedEta7TeVRandomLowError;

TGraphErrors* combined7TeVHighError;
TGraphErrors* combined7TeVLowError;
TGraphErrors* combinedEta7TeVHighError;
TGraphErrors* combinedEta7TeVLowError;

TGraphAsymmErrors *AsymmCombined7TeVHighError;
TGraphAsymmErrors *AsymmCombined7TeVLowError;

TH1D *histcombined7TeVHighError;
TH1D *histcombined7TeVLowError;
TH1D *histcombined7TeVIncRatioHighError;
TH1D *histcombined7TeVIncRatioLowError;


TGraphAsymmErrors* graphAllGammaPi0SysError;
TH1D *HighErrorAllGammaPi0;
TH1D *LowErrorAllGammaPi0;


void FindCocktailFile(TString eventCutSelection, TString gammaCutSelection, TString centrality, TString option){
   
   fcmult = 1;
   
   ffive = kFALSE;
   
   TString finderType(gammaCutSelection(0,1));
   TString extraSignals(eventCutSelection(GetEventRejectExtraSignalsCutPosition(),1));
   
   
   Int_t FinderType = 0;      
   FinderType = finderType.Atoi();
   
   TString centralityCutNumberLow = eventCutSelection(GetEventCentralityMinCutPosition(),1);
   TString centralityCutNumberHigh = eventCutSelection(GetEventCentralityMaxCutPosition(),1);
   Int_t centCutNumberI  = centralityCutNumberLow.Atoi();
   Int_t centCutNumberBI = centralityCutNumberHigh.Atoi();

   if(!option.CompareTo("PbPb_2.76TeV")){    
      if(!centrality.CompareTo("0-5%")){
         fcocktailFunc = "cocktail_HI0005_oHagBW";//cocktail_HI_0005_QCD";
         fcocktailFuncSyst = "cocktail_HI_0005_QCD";
         fcocktailFuncNarrow = "cocktail_HI_0005_QCD";
         fcocktailFuncWide = "cocktail_HI_0005_QCD";
         fcocktailFuncChargedPions =  "cocktail_HI_0005_QCD";
         ffive = kTRUE;
      }
      if(!centrality.CompareTo("5-10%")){
         fcocktailFunc = "cocktail_HI0510_qcdBW";//cocktail_HI_0510_QCD";
         fcocktailFuncSyst = "cocktail_HI_0510_QCD";
         fcocktailFuncNarrow = "cocktail_HI_0510_QCD";
         fcocktailFuncWide = "cocktail_HI_0510_QCD";
         fcocktailFuncChargedPions =  "cocktail_HI_0510_QCD";
         ffive = kTRUE;
      }
      if(!centrality.CompareTo("0-10%")){
         fcocktailFunc = "cocktail_PbPb0010_WeightedQGP";
         if(extraSignals.CompareTo("3")) 
            fcocktailFunc = "cocktail_PbPb0010_WeightedQGP";//cocktail_HI_0010_qcdBWEta";
         fcocktailFuncSyst = "cocktail_HI_0010_qcdBWEta";
         fcocktailFuncNarrow = "cocktail_HI_0010_qcdBWEta";
         fcocktailFuncWide = "cocktail_HI_0010_qcdBWEta";
         fcocktailFuncChargedPions =  "cocktail_HI_0010_qcdBWEta";
      }
      if(!centrality.CompareTo("10-20%")){
         fcocktailFunc = "cocktail_HI1020_oHagBW";//cocktail_HI_1020_oHagBWEta_CorrCent";
         if(extraSignals.CompareTo("3")) 
            fcocktailFunc = "cocktail_PbPb1020_WeightedQGP";
         fcocktailFuncSyst = "cocktail_HI_1020_oHagBWEta_CorrCent";
         fcocktailFuncNarrow = "cocktail_HI_1020_oHagBWEta_CorrCent";
         fcocktailFuncWide = "cocktail_HI_1020_oHagBWEta_CorrCent";
         fcocktailFuncChargedPions =  "cocktail_HI_1020_oHagBWEta_CorrCent";
      }
      if(!centrality.CompareTo("0-40%")){
         fcocktailFunc = "cocktail_HI_00-40_QCD_mtFixBW_9002972094503042212000010420020";
         fcocktailFuncChargedPions =  "cocktail_HI_00-40_QCD_mtFixBW_9002972094503042212000010420020";
         fcocktailFuncSyst = "cocktail_HI_00-40_QCD_mtFixBW_9002972094503042212000010420020";
         fcocktailFuncNarrow = "cocktail_HI_00-40_QCD_mtFixBW_9002972094503042212000010420020";
         fcocktailFuncWide = "cocktail_HI_00-40_QCD_mtFixBW_9002972094503042212000010420020";
      }
      if(!centrality.CompareTo("0-20%")){
         if(FinderType) fcocktailFunc = "cocktail_HI_0020_qcdOffline";
         //fcocktailFunc = "cocktail_HI_0020_modHagBWEta_HighR"; //9
         //if(extraSignals.CompareTo("3"))
         //fcocktailFunc = "cocktail_PbPb0020NewMC_MCoHag";//cocktail_PbPb0020NewMC_oHag";
         //fcocktailFunc = "cocktail_PbPb0020NewMC_MCoHag";//cocktail_PbPb0020NewMC_oHag";
         fcocktailFunc = "cocktail_PbPb0020NewMC_oHag";
         fcocktailFuncSyst = "cocktail_HI_0020_modHagBWEta";
         fcocktailFuncNarrow = "cocktail_HI_0020_modHagBWEta";
         fcocktailFuncWide = "cocktail_HI_0020_modHagBWEta";
         fcocktailFuncChargedPions =  "cocktail_HI_0020_modHagBWEta";
      }
      if(!centrality.CompareTo("20-40%")){
         fcocktailFuncNarrow = "cocktail_HI_2040_qcdBWEta";
         fcocktailFuncWide =   "cocktail_HI_2040_qcdBWEta";
         fcocktailFuncSyst = "cocktail_HI_2040_qcdBWEta";
         fcocktailFunc = "cocktail_PbPb0020newMC_oHag";//cocktail_HI_2040_qcdBWEta";//cocktail_HI_oHag_12400010420927700237000000_010220450000";//cocktail_HI_2040_qcdBWEta";
         /* if(extraSignals.CompareTo("3")) */
         /*    fcocktailFunc = "cocktail_PbPb2040_WeightedQGP"; */
         if(FinderType) fcocktailFunc = "cocktail_HI_2040_qcdOffline";
         fcocktailFuncChargedPions =  "cocktail_HI_2040_qcdBWEta";
      }
      if(!centrality.CompareTo("40-80%")){
         fcocktailFuncNarrow = "cocktail_HI_40-80_oHag";
         fcocktailFuncWide =   "cocktail_HI_40-80_oHag";
         fcocktailFuncSyst = "cocktail_HI_40-80_oHag";
         //fcocktailFunc = "cocktail_HI_4080_qcdpp";//cocktail_PbPb4080NewMC_oHag";//";
         fcocktailFunc = "cocktail_PbPb4080NewMC_oHag";
         //fcocktailFunc = "cocktail_PbPb4080NewMC_qcd";
         if(FinderType) fcocktailFunc = "cocktail_HI_4080_qcdOffline";
         fcocktailFuncChargedPions =  "cocktail_HI_40-80_oHag";
      }
      if(!centrality.CompareTo("60-80%")){
         fcocktailFuncNarrow = "cocktail_HI_40-80_oHag";
         fcocktailFuncWide =   "cocktail_HI_40-80_oHag";
         fcocktailFuncSyst = "cocktail_HI_40-80_oHag";
         fcocktailFunc = "cocktail_HI_4080_qcdpp";
         if(FinderType) fcocktailFunc = "cocktail_HI_4080_qcdOffline";
         fcocktailFuncChargedPions =  "cocktail_HI_40-80_oHag";
      }
      if(!centrality.CompareTo("0-80%")){
         fcocktailFuncNarrow = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncWide =   "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncSyst = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFunc = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncChargedPions =  "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
      }
      if(!centrality.CompareTo("0-90%")){
         fcocktailFuncNarrow = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncWide =   "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncSyst = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFunc = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncChargedPions =  "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
      }
   
      fcmult = 0;
      for(Int_t i = centCutNumberI; i<centCutNumberBI;i++){
         if(ffive){
            fcmult += multiplicityCent[i];
            cout<<multiplicityCent[i]<<endl;
         }
         else{
            fcmult += multiplicity[i];
            cout<<multiplicity[i]<<endl;
         }
         cout<<" Sum:"<<fcmult<<endl;
      }
      cout<<" / "<<centCutNumberBI-centCutNumberI<<"  =  "<<endl;
      fcmult = fcmult/(centCutNumberBI-centCutNumberI);
      cout<<fcmult<<endl;
   }
   else if(!option.CompareTo("7TeV")){    
      //fcocktailFuncNarrow = "cocktail_7TeV_pi0LevyetaLevy";
      //fcocktailFuncWide =   "cocktail_7TeV_pi0LevyetaLevy";
      //fcocktailFuncSyst = "cocktail_7TeV_pi0LevyetaLevy";
       fcocktailFuncNarrow = "cocktail_pp_7TeV";
       fcocktailFuncWide =   "cocktail_pp_7TeV";
       fcocktailFuncSyst = "cocktail_pp_7TeV";
      //fcocktailFunc = "cocktail_7TeV_pi0LevyetaLevy";
      //fcocktailFunc = "cocktail_pp7TeV_HagLevy";
      //fcocktailFunc = "cocktail_pp7TeV_HagHag";;
      //fcocktailFunc = "cocktail_pp7TeV_oHagHag";
      //fcocktailFunc = "cocktail_pp7TeV_TsallisTsallis";
       fcocktailFunc = "cocktail_pp_7TeV";
      //fcocktailFunc = "cocktail_pp7TeV_QCDTsallis";
      //fcocktailFunc = "cocktail_7TeV_pi0LevyetaLevy";
      //fcocktailFunc = "cocktail_pp7TeV_qcdqcd";
      //fcocktailFuncChargedPions = "cocktail_276TeV";
       fcocktailFuncChargedPions = "";
   }
   else if(!option.CompareTo("2.76TeV")){
      fcocktailFuncNarrow = "cocktail_pp2760GeV_oHagmt";
      fcocktailFuncWide =   "cocktail_pp2760GeV_oHagmt";
      fcocktailFuncSyst = "cocktail_pp2760GeV_oHagmt";
      fcocktailFunc = "cocktail_pp2760GeVEta_Hag";//cocktail_pp2760GeV_oHagmt";//cocktail_pp2760GeV_qcd";//cocktail_pp2760GeV_qcd";//cocktail_pp2760GeV_oHagmt";//cocktail_pp2760GeVmtunfold_oHag_0000013002093663003800000_016310310090";//ocktail_pp2760GeV_oHag_0000013002093663003800000_016310310090";//cocktail_pp2760GeV_oHagmt";//cocktail_pp2760GeV_TsallisTsallis";////";
      fcocktailFuncChargedPions = "cocktail_pp2760GeV_oHagmt";
   }
   else if(!option.CompareTo("900GeV")){
      fcocktailFuncNarrow = "cocktail_pp2760GeV_oHagmt";
      fcocktailFuncWide =   "cocktail_pp2760GeV_oHagmt";
      fcocktailFuncSyst = "cocktail_pp2760GeV_oHagmt";
      fcocktailFunc = "cocktail_pp900GeV_oHag";//cocktail_pp900GeV_Iqcd";
      fcocktailFuncChargedPions = "cocktail_pp2760GeV_oHagmt";
   }
   else if(!option.CompareTo("pPb_5.023TeV")){
	  fcmult = 6.9;
      fcocktailFuncNarrow = "cocktail_pPb_DPMJET_MB";
      fcocktailFuncWide =   "cocktail_pPb_DPMJET_MB";
      fcocktailFuncSyst = "cocktail_pPb_DPMJET_MB";
      fcocktailFunc = "cocktail_pPb_HIJING_MB_pileupCorr_20042014_EtaAndPhimeasured";
	  fcocktailFuncModA = "cocktail_pPb_HIJING_MB_pileupCorr_05042014_etameasured";
	  fcocktailFuncModB = "cocktail_pPb_HIJING_MB_pileupCorr_05042014_etameasured";
	  fcocktailFuncMtScaledEta = "cocktail_pPb_HIJING_MB_pileupCorr";
      fcocktailFuncChargedPions = "cocktail_pPb_HIJING_MB_pileupCorr_chargedPionsScaledEta";
   }
   else {
      fcocktailFuncNarrow = "cocktail_Test";
      fcocktailFuncWide =   "cocktail_Test";
      fcocktailFuncSyst = "cocktail_Test";
      fcocktailFunc = "cocktail_Test";
      fcocktailFuncChargedPions =  "cocktail_Test";
   }

   //return 0;

}

