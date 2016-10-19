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