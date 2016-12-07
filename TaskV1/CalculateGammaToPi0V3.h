#include "../CommonHeaders/ExtractSignalBinning.h"

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

TString textPrefix2;
TString Meson_text;
TString textPi0New;

//declaration for printing logo 
Float_t purityDrawingCoordinates[5]                         = {0.25, 0.76, 0.2, 0.05, 0.0};
Float_t convProbDrawingCoordinates[5]                       = {0.25, 0.76, 0.2, 0.05, 0.0};
Float_t recoEffDrawingCoordinates[5]                        = {0.25, 0.76, 0.2, 0.05, 0.0};
Float_t incRatioDrawingCoordinates[5]                       = {0.65, 0.75, 0.2, 0.05, 0.0};
Float_t incRatioAllDrawingCoordinates[5]                    = {0.45, 0.75, 0.2, 0.05, 0.0};
Float_t decayRatioDrawingCoordinates[5]                     = {0.2, 0.75, 0.2, 0.05, 0.0};

Double_t binningRAA[12]                                     = {0.4,0.8,1.0,1.4,1.8,2.2,2.6,3.0,4.0,5.0,6.0,8.0};
Double_t multiplicity[10]                                   = {1502.7,923.26,558.68,321.20,171.67,85.13,38.51,15.78,6.32,2.63};
//                                                              0-1,   1-2,   2-3,   3-4,   4-5,   5-6,  6-7,  7-8,  8-9, 9-10
Double_t multiplicityCent[6]                                = {1686.87,1319.89, 1031.9,807.90,627.99,483.95};
//                                                                      0-05,   05-1

Bool_t pictDrawingOptions[4]                                = {kFALSE, kFALSE, kFALSE, kTRUE};

TFile *fileGamma                                            = NULL;
TFile *filePi0                                              = NULL;
TFile *fileEta                                              = NULL;
TFile *filePhi                                              = NULL;
TFile *fileQuarkMatter                                      = NULL;
TFile *fileSystErrorPbPb                                    = NULL;
TFile* fileTheoryCompilation                                = NULL;
TFile *cocktailFile                                         = NULL;
TFile *cocktailFileModA                                     = NULL;
TFile *cocktailFileModB                                     = NULL;
TFile *cocktailFileMtScaledEta                              = NULL;
TFile *cocktailFileChargedPions                             = NULL;
TFile *fileCorrectedOutput                                  = NULL;
TDirectoryFile* directoryGamma                              = NULL;
TDirectoryFile *cocktailDirModA                             = NULL;
TDirectoryFile *cocktailDirModB                             = NULL;
TDirectoryFile *cocktailDir                                 = NULL;
TDirectoryFile *cocktailDirMtScaledEta                      = NULL;
TDirectoryFile *cocktailDirChargedPions                     = NULL;

TString fcocktailFunc                                       = "";
TString fcocktailFuncModA                                   = "";
TString fcocktailFuncModB                                   = "";
TString fcocktailFuncMtScaledEta                            = "";
TString fcocktailFuncChargedPions                           = "";
TString fcocktailFuncNarrow                                 = "";
TString fcocktailFuncWide                                   = "";
TString fcocktailFuncSyst                                   = "";
Double_t fcmult                                             = 0;
Bool_t ffive                                                = kFALSE;

TF1 *One                                                    = NULL;

TH1D *histoBinningDoubleRatioRaa                            = NULL;

TH1D *histoCorrectedPi0YieldFit                             = NULL;
TH1D *histoCorrectedPi0YieldFitWide                         = NULL;
TH1D *histoCorrectedPi0YieldFitNarrow                       = NULL;

TH1D *histoCorrectedEtaYieldNormalEff                       = NULL;
TH1D *histoCorrectedEtaYield                                = NULL;
TH1D* histoCorrectedEtaYieldPatched                         = NULL;

TH1D *histoCorrectedPhiYieldStat                            = NULL;
TH1D *histoCorrectedPhiYieldSyst                            = NULL;
TF1* fitPhiYield                                            = NULL;

TH1D *cocktailAllGammaNLO                                   = NULL;
TGraphAsymmErrors *graphDirectPhotonNLO                     = NULL;
TGraphAsymmErrors *graphPromptPhotonNLO                     = NULL;
TGraphAsymmErrors *graphFragmentationPhotonNLO              = NULL;
TGraphAsymmErrors* graphDirectPhotonNLOCopy                 = NULL;
TGraphAsymmErrors* NLODoubleRatio                           = NULL;
TGraphAsymmErrors* NLO                                      = NULL;

Double_t* xVal                                              = NULL;
Double_t* xErr                                              = NULL;
Double_t* yVal                                              = NULL;
Double_t* yErrUp                                            = NULL;
Double_t* yErrDown                                          = NULL;

TF1* fitNLODirectPhoton                                     = NULL;
TF1* fitNLOPromptPhoton                                     = NULL;
TF1* fitNLOFragmentationPhoton                              = NULL;

TF1 *fitPi0YieldA                                           = NULL;
TF1 *fitPi0YieldAHighpt                                     = NULL;
TF1 *fitPi0YieldB                                           = NULL;
TF1 *fitEtaYieldA                                           = NULL;
TF1 *fitEtaYieldB                                           = NULL;
TF1 *fitEtaYieldC                                           = NULL;
TF1 *fitPi0YieldC                                           = NULL;
TF1 *fitPi0YieldCWide                                       = NULL;
TF1 *fitPi0YieldCNarrow                                     = NULL;

TF1 *ConversionGammaFitA                                    = NULL;
TF1 *ConversionGammaFitB                                    = NULL;

TH1D *histRatioNLODirectPhoton                              = NULL;

TH1D *histRatioConversionEtaA                               = NULL;
TH1D *histRatioConversionEtaB                               = NULL;
TH1D* histRatioConversionEtaC                               = NULL;
TH1D *histRatioConversionPi0A                               = NULL;
TH1D *histRatioConversionPi0AHighpt                         = NULL;
TH1D *histRatioConversionPi0B                               = NULL;
TH1D *histRatioConversionPi0C                               = NULL;
TH1D *histRatioConversionGammaA                             = NULL;
TH1D *histRatioConversionGammaB                             = NULL;

TH1D *histoTruePi0MCData                                    = NULL;
TH1D *histoNormalPi0MCData                                  = NULL;
TH1D *histoPurityGammaMCData                                = NULL;

TH1D *histoGammaSpecCorrPurity                              = NULL;
TH1D *histoGammaSpecMCAll                                   = NULL;
TH1D *histoGammaSpecCorrESDMC                               = NULL;
TH1D *histoMCDecaySumGammaPt                                = NULL;
TH1D *histoMCDecayPi0GammaPt                                = NULL;
TH1D *histoMCDecayEtaGammaPt                                = NULL;
TH1D *histoMCDecayEtapGammaPt                               = NULL;
TH1D *histoMCDecayOmegaGammaPt                              = NULL;
TH1D *histoMCDecayRho0GammaPt                               = NULL;
TH1D *histoMCDecayK0sGammaPt                                = NULL;
TH1D *histoMCDecayPhiGammaPt                                = NULL;
TH1D *histMCAllMinusDecay                                   = NULL;

TH1D *histoCorrectedPi0YieldNormalEff                       = NULL;
TH1D *histoCorrectedPi0Yield                                = NULL;
TH1D *histoCorrectedPi0YieldWide                            = NULL;
TH1D *histoCorrectedPi0YieldNarrow                          = NULL;
TH1D *histoMCYieldMeson                                     = NULL;
TH1D *histoMCYieldMesonOldBin                               = NULL;

TH1D *histoIncRatioFitPurity                                = NULL;
TH1D *histoIncRatioFitPurityWide                            = NULL;
TH1D *histoIncRatioFitPurityNarrow                          = NULL;
TH1D *histoIncRatioHighFitPurity                            = NULL;
TH1D *histoIncRatioLowFitPurity                             = NULL;
TH1D *histoIncRatioCombinedPurity                           = NULL;
TH1D *histoIncRatioPurity                                   = NULL;
TH1D *histoIncRatioPurityTrueEff                            = NULL;
TH1D *histoIncRatioPurityTrueEffWide                        = NULL;
TH1D *histoIncRatioPurityTrueEffNarrow                      = NULL;

TH1D *histoMCIncRatio                                       = NULL;
TH1D *histoMCesdIncRatio                                    = NULL;
TH1D *histoIncRatioGammaMC                                  = NULL;

TH1D *histoDecayRatioSumGamma                               = NULL;
TH1D *histoDecayRatioPi0Gamma                               = NULL;
TH1D *histoDecayRatioEtaGamma                               = NULL;
TH1D *histoDecayRatioEtapGamma                              = NULL;
TH1D *histoDecayRatioOmegaGamma                             = NULL;
TH1D *histoDecayRatioRho0Gamma                              = NULL;
TH1D *histoDecayRatioK0sGamma                               = NULL;
TH1D *histoDecayRatioPhiGamma                               = NULL;

TH1D *cocktailAllGamma                                      = NULL;
TH1D *cocktailPi0Gamma                                      = NULL;
TH1D *cocktailEtaGamma                                      = NULL;
TH1D *cocktailEtapGamma                                     = NULL;
TH1D *cocktailOmegaGamma                                    = NULL;
TH1D *cocktailPhiGamma                                      = NULL;
TH1D *cocktailRhoGamma                                      = NULL;
TH1D *cocktailSigmaGamma                                    = NULL;

TH1D *cocktailAllGammaMtScaledEta                           = NULL;
TH1D *cocktailPi0GammaMtScaledEta                           = NULL;
TH1D *cocktailEtaGammaMtScaledEta                           = NULL;
TH1D *cocktailEtapGammaMtScaledEta                          = NULL;
TH1D *cocktailOmegaGammaMtScaledEta                         = NULL;
TH1D *cocktailPhiGammaMtScaledEta                           = NULL;
TH1D *cocktailRhoGammaMtScaledEta                           = NULL;
TH1D *cocktailSigmaGammaMtScaledEta                         = NULL;

TH1D *cocktailAllGammaChargedPions                          = NULL;
TH1D *cocktailPi0GammaChargedPions                          = NULL;
TH1D *cocktailEtaGammaChargedPions                          = NULL;
TH1D *cocktailEtapGammaChargedPions                         = NULL;
TH1D *cocktailOmegaGammaChargedPions                        = NULL;
TH1D *cocktailPhiGammaChargedPions                          = NULL;
TH1D *cocktailRhoGammaChargedPions                          = NULL;
TH1D *cocktailSigmaGammaChargedPions                        = NULL;

TH1D *cocktailPi0                                           = NULL;
TH1D *cocktailEta                                           = NULL;
TH1D *cocktailPhi                                           = NULL;

TH1D *cocktailAllGammaPi0                                   = NULL;
TH1D *cocktailPi0GammaPi0                                   = NULL;
TH1D *cocktailEtaGammaPi0                                   = NULL;
TH1D *cocktailEtapGammaPi0                                  = NULL;
TH1D *cocktailOmegaGammaPi0                                 = NULL;
TH1D *cocktailPhiGammaPi0                                   = NULL;
TH1D *cocktailRhoGammaPi0                                   = NULL;
TH1D *cocktailSigmaGammaPi0                                 = NULL;

TH1D *cocktailAllGammaPi0ModA                               = NULL;
TH1D *cocktailAllGammaPi0ModB                               = NULL;
TH1D *cocktailAllGammaModA                                  = NULL;
TH1D *cocktailAllGammaModB                                  = NULL;
TH1D *cocktailPi0ModA                                       = NULL;
TH1D *cocktailPi0ModB                                       = NULL;

TH1D *cocktailPi0MtScaledEta                                = NULL;
TH1D *cocktailEtaMtScaledEta                                = NULL;
TH1D *cocktailAllGammaPi0MtScaledEta                        = NULL;
TH1D *cocktailPi0GammaPi0MtScaledEta                        = NULL;
TH1D *cocktailEtaGammaPi0MtScaledEta                        = NULL;
TH1D *cocktailEtapGammaPi0MtScaledEta                       = NULL;
TH1D *cocktailOmegaGammaPi0MtScaledEta                      = NULL;
TH1D *cocktailPhiGammaPi0MtScaledEta                        = NULL;
TH1D *cocktailRhoGammaPi0MtScaledEta                        = NULL;
TH1D *cocktailSigmaGammaPi0MtScaledEta                      = NULL;

TH1D *cocktailPi0ChargedPions                               = NULL;
TH1D *cocktailEtaChargedPions                               = NULL;
TH1D *cocktailAllGammaPi0ChargedPions                       = NULL;
TH1D *cocktailPi0GammaPi0ChargedPions                       = NULL;
TH1D *cocktailEtaGammaPi0ChargedPions                       = NULL;
TH1D *cocktailEtapGammaPi0ChargedPions                      = NULL;
TH1D *cocktailOmegaGammaPi0ChargedPions                     = NULL;
TH1D *cocktailPhiGammaPi0ChargedPions                       = NULL;
TH1D *cocktailRhoGammaPi0ChargedPions                       = NULL;
TH1D *cocktailSigmaGammaPi0ChargedPions                     = NULL;

TH1D *cocktailPi0Rebinned                                   = NULL;
TH1D *cocktailAllGammaRebinned                              = NULL;

TH1D *histoDoubleRatioConversionNormalPurity                = NULL;
TH1D *histoDoubleRatioConversionTrueEffPurity               = NULL;
TH1D *histoDoubleRatioConversionTrueEffPurityModA           = NULL;
TH1D *histoDoubleRatioConversionTrueEffPurityModB           = NULL;
TH1D *histoDoubleRatioConversionTrueEffPurityMtScaledEta    = NULL;
TH1D *histoDoubleRatioConversionTrueEffPurityChargedPions   = NULL;
TH1D *histoDoubleRatioConversionOnlyGamma                   = NULL;
TH1D *histoDoubleRatioConversionTrueEffPurityWide           = NULL;
TH1D *histoDoubleRatioConversionTrueEffPurityNarrow         = NULL;

TH1D *histoDoubleRatioFitPi0YieldPurity                     = NULL;
TH1D *histoDoubleRatioFitPi0YieldPurityModA                 = NULL;
TH1D *histoDoubleRatioFitPi0YieldPurityModB                 = NULL;
TH1D *histoDoubleRatioFitPi0YieldPurityMtScaledEta          = NULL;
TH1D *histoDoubleRatioFitPi0YieldPurityChargedPions         = NULL;
TH1D *histoDoubleRatioConversionLowFitPurity                = NULL;
TH1D *histoDoubleRatioConversionHighFitPurity               = NULL;
TH1D *histoDoubleRatioFitPi0YieldPurityWide                 = NULL;
TH1D *histoDoubleRatioFitPi0YieldPurityNarrow               = NULL;

TH1D *histoDoubleRatioCombinedPurity                        = NULL;
TH1D *histoMCDoubleRatioSum                                 = NULL;
TH1D *histoMCesdDoubleRatioSum                              = NULL;

TF1 *FitRandomHigh                                          = NULL;
TF1 *FitRandomLow                                           = NULL;
TF1 *NLOdoubleRatioFit                                      = NULL;
TF1 *cocktailFitAllGammaForNLO                              = NULL;

TH1D *histoDirectPhotonSpectrum                             = NULL;
TH1D *histoThermalPhotonSpectrum                            = NULL;
TH1D *histoPromptPhotonSpectrum                             = NULL;

TH1D *cocktailAllGammaRebinnedModA                          = NULL;
TH1D *cocktailPi0RebinnedModA                               = NULL;
TH1D *cocktailAllGammaRebinnedModB                          = NULL;
TH1D *cocktailPi0RebinnedModB                               = NULL;

Color_t colorCocktailAllDecay               = kBlack;
Color_t colorCocktailPi0                    = kRed+2;
Color_t colorCocktailEta                    = kBlue+1;
Color_t colorCocktailEtaP                   = kOrange+1;
Color_t colorCocktailOmega                  = kYellow+2;
Color_t colorCocktailPhi                    = kViolet;
Color_t colorCocktailRho0                   = kAzure-2;

TString fEventCutSelection                  = "";
TString fGammaCutSelection                  = "";
TString fClusterCutSelection                = "";
TString fElectronCutSelection               = "";
TString fMesonCutSelection                  = "";
TString centralityCutNumberLow              = "";
Int_t centCutNumberI                        = 0;
TString centrality                          = "";
TString centralityAdd                       = "";
Double_t fNcoll                             = 0;
Double_t fNcollErr                          = 0;
TString triggerCutNumber                    = "";
TString subTriggerCutNumber                 = "";

Bool_t kDoPileup                            = kFALSE;
Double_t textSizeSpectra                    = 0.035;

Double_t fitMinPt                           = 0.4;
Double_t fitMaxPt                           = 0;
TString forOutput                           = "";
TH2F * histo2DSpectraPi0                    = NULL;

ifstream fileSysErrGamma;
Int_t nPointsGamma                              = 0;
Double_t relSystErrorGammaUp[50];
Double_t relSystErrorGammaDown[50];
Double_t relSystErrorWOMaterialGammaUp[50];
Double_t relSystErrorWOMaterialGammaDown[50];
Double_t systErrorGammaUp[50];
Double_t systErrorGammaDown[50];

ifstream fileSysErrInclRatio;
Int_t nPointsInclRatio                          = 0;
Double_t relSystErrorInclRatioUp[50];
Double_t relSystErrorInclRatioDown[50];
Double_t relSystErrorWOMaterialInclRatioUp[50];
Double_t relSystErrorWOMaterialInclRatioDown[50];
Double_t systErrorInclRatioUp[50];
Double_t systErrorInclRatioDown[50];

ifstream fileSysErrDoubleRatio;
Int_t nPointsDoubleRatio                        = 0;
Double_t relSystErrorDoubleRatioUp[50];
Double_t relSystErrorDoubleRatioDown[50];
Double_t relSystErrorWOMaterialDoubleRatioUp[50];
Double_t relSystErrorWOMaterialDoubleRatioDown[50];
Double_t systErrorDoubleRatioUp[50];
Double_t systErrorDoubleRatioDown[50];

TGraphAsymmErrors*  graphGammaYieldSysErr       = NULL;
TGraphAsymmErrors*  graphInclRatioSysErr        = NULL;
TGraphAsymmErrors*  graphDoubleRatioSysErr      = NULL;
TGraphAsymmErrors*  graphDoubleRatioFitSysErr   = NULL;

void SeparateCutnumberString(TString cutSelStr, Int_t mode)
{
   TString fCutSelection                       = cutSelStr;
   if (mode == 9){
      ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
      fEventCutSelection                      = fGammaCutSelection(0,7);
      fGammaCutSelection                      = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
      cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
   } else {
      ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
   }
   centralityCutNumberLow             = fEventCutSelection(GetEventCentralityMinCutPosition(),1);
   centCutNumberI                     = centralityCutNumberLow.Atoi();
   centrality                         = GetCentralityString(fEventCutSelection);
   centralityAdd                      = GetCentralityStringWoPer(fEventCutSelection);
   triggerCutNumber                   = fEventCutSelection(3,1);
   subTriggerCutNumber                = fEventCutSelection(4,1);
   fNcoll                             = GetNCollFromCutNumber(fEventCutSelection);
   fNcollErr                          = GetNCollErrFromCutNumber(fEventCutSelection);
}

void CreateNamingStrings(TString nameMeson, TString isMC)
{
   textPi0New                          = Form("Gamma_%s",nameMeson.Data());
   if (isMC.CompareTo("kTRUE") == 0){
      textPrefix2                             = "recMC";
      pictDrawingOptions[1]                   = kTRUE;
   } else {
      textPrefix2                             = "data";
      pictDrawingOptions[1]                   = kFALSE;
   }

   TString mesonType;
   if ((nameMeson.CompareTo("Pi0") == 0) || (nameMeson.CompareTo("Pi0EtaBinning") == 0)){
      pictDrawingOptions[3]                   = kTRUE;
      Meson_text                              = "#pi^{0}";
      mesonType                               = "Pi0";
   } else {
      pictDrawingOptions[3]                   = kFALSE;
      Meson_text                              = "#eta";
      mesonType                               = "Eta";
   }
}

void PlotLatexLegend(Double_t xPosition=0.5,Double_t yPosition=0.5,Double_t textSizeSpectra=0.04,TString collisionSystem ="",TString detectionProcess ="", Int_t labelNums = 3)
{
   TLatex* labelPi0Plot; TLatex* labelEnergyPlot;TLatex* labelDetProcPlot; 
   labelEnergyPlot = new TLatex(xPosition, yPosition+3*0.85*textSizeSpectra,collisionSystem.Data());
   if(labelNums<3)
   {
      labelPi0Plot    = new TLatex(xPosition, yPosition+200*0.85*textSizeSpectra,"#pi^{0} #rightarrow #gamma#gamma");
      labelDetProcPlot= new TLatex(xPosition, yPosition+2*0.85*textSizeSpectra,detectionProcess.Data());
   }else{
      labelPi0Plot    = new TLatex(xPosition, yPosition+2*0.85*textSizeSpectra,"#pi^{0} #rightarrow #gamma#gamma");
      labelDetProcPlot= new TLatex(xPosition, yPosition+0.85*textSizeSpectra,detectionProcess.Data());
   }
   SetStyleTLatex( labelEnergyPlot, 0.85*textSizeSpectra,4);
   SetStyleTLatex( labelPi0Plot, 0.85*textSizeSpectra,4);
   SetStyleTLatex( labelDetProcPlot, 0.85*textSizeSpectra,4);

   labelEnergyPlot->Draw();
   if(labelNums==3)
      labelPi0Plot->Draw();
   labelDetProcPlot->Draw();
}

void drawLatex(TString textString, Double_t xCoord, Double_t yCoord, Color_t txtclr = kBlack,Size_t txtsize = 0.1){
   TLatex *TextLatex = new TLatex(xCoord,yCoord,Form("#font[62]{%s}",textString.Data()));
   TextLatex->SetTextColor(txtclr);
   TextLatex->SetTextSize(txtsize);
   TextLatex->Draw();
}

Double_t GetUpperLimit(Double_t mean, Double_t statErr, Double_t sysErr, Double_t confidenceLevel, Double_t& confidenceLevelReached, Double_t accuracy = 1e-9, Int_t maxNIterations = 1e8) {

    // function to return upper limit on photon excess, using a Bayesian approach
    // with the heaviside function used as prior (excluding R_gamma < 1)

    // set range in R_gamma
    Double_t    minRGamma           = 0.;
    Double_t    maxRGamma           = 4.;

    // conditional probability to measure R_gamma above R_gamma_true = 1
    TF1         condProb("condProb", Form("[0]*(TMath::Gaus(x, [1], [2])*((x>=1)*1 + (x<1)*0))"),minRGamma, maxRGamma);
    condProb.SetParameter(0, 1.);
    condProb.SetParameter(1, mean);
    condProb.SetParameter(2, TMath::Sqrt(statErr*statErr + sysErr*sysErr));

    // normalize conditional probability to one
    Int_t       np                  = 10000;
    Double_t*   x                   = new Double_t[np];
    Double_t*   w                   = new Double_t[np];
    condProb.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
    Double_t    norm                = condProb.IntegralFast(np,x,w,minRGamma,maxRGamma);
    condProb.SetParameter(0, 1/norm);

    // iteratively find upper limit (interval bisection)
    Double_t    upperLimit          = (maxRGamma-1)/2;
    Double_t    upperLimitPrev      = upperLimit;
    Double_t    step                = 0.;
    Int_t       nIterations         = 0;
    while (((condProb.IntegralFast(np,x,w,1,upperLimit) < (confidenceLevel-accuracy)) || (condProb.IntegralFast(np,x,w,1,upperLimit) > (confidenceLevel+accuracy))) && nIterations < maxNIterations) {
        
        if (condProb.IntegralFast(np,x,w,1,upperLimit) > confidenceLevel)
            step                    = - TMath::Abs(upperLimit-1)/2;
        else
            step                    = TMath::Abs(upperLimitPrev-upperLimit)/2;
        upperLimitPrev              = upperLimit;
        upperLimit                  = upperLimit + step;

        //if ( !(nIterations%10) ) cout << "   condProb.IntegralFast(1, " << upperLimit << ") = " << condProb.IntegralFast(np,x,w,1, upperLimit) << endl;

        nIterations++;
    }

    confidenceLevelReached          = condProb.IntegralFast(np,x,w,1, upperLimit);
    return upperLimit;
}

TH1D* GetUpperLimitsHisto(TH1D* histo, TGraphAsymmErrors* sysErrGraph, Double_t confidenceLevel = 0.95, Double_t accuracy = 0.004, Int_t maxNIterations = 1e3) {

    cout << endl;
    cout << "*************************************************************" << endl;
    cout << "**                                                         **" << endl;
    cout << "**      STARTING UPPER LIMIT CALCULATION                   **" << endl;
    cout << "**                                                         **" << endl;
    cout << "*************************************************************" << endl;
    cout << endl;

    // upper limits histo
    TH1D*       upperLimits             = (TH1D*)histo->Clone("upperLimits");

    // get graph quantities
    Int_t       nBinsGraph              = sysErrGraph->GetN();
    Double_t*   xValueGraph             = sysErrGraph->GetX();
    Double_t*   xErrorLowGraph          = sysErrGraph->GetEXlow();
    Double_t*   xErrorHighGraph         = sysErrGraph->GetEXhigh();

    // fill upper limits histo
    for (Int_t i=1; i<histo->GetNbinsX()+1; i++) {
        if (!histo->GetBinContent(i)) {
            upperLimits->SetBinContent( i, 0);
            upperLimits->SetBinError(   i, 0);
        } else {
            for (Int_t j=0; j<sysErrGraph->GetN(); j++) {
                if (xValueGraph[j] == histo->GetBinCenter(i)) {

                    Double_t reached    = 0.;

                    upperLimits->SetBinContent( i, GetUpperLimit(histo->GetBinContent(i),histo->GetBinError(i),(xErrorLowGraph[j]+xErrorHighGraph[j])/2,confidenceLevel,reached,accuracy,maxNIterations));
                    upperLimits->SetBinError(   i, 0);

                    cout << "p_T = " << histo->GetBinCenter(i) << ": " << histo->GetBinContent(i) << " -> " << upperLimits->GetBinContent(i) << " at CL = " << reached << endl;
                }
            }
        }
    }

    cout << endl;
    cout << "*************************************************************" << endl;
    cout << "**                                                         **" << endl;
    cout << "**      DONE WITH UPPER LIMIT CALCULATION                  **" << endl;
    cout << "**                                                         **" << endl;
    cout << "*************************************************************" << endl;
    cout << endl;
    
    return upperLimits;
}