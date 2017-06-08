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

void FindCocktailFile(TString eventCutSelection, TString gammaCutSelection, TString centrality, TString option){
   
   fcmult                           = 1;
   ffive                            = kFALSE;
   
   TString finderType(gammaCutSelection(0,1));
   TString extraSignals(eventCutSelection(GetEventRejectExtraSignalsCutPosition(),1));
   
   Int_t FinderType                 = 0;
   FinderType                       = CutNumberToInteger(finderType);
   
   TString centralityCutNumberLow   = eventCutSelection(GetEventCentralityMinCutPosition(),1);
   TString centralityCutNumberHigh  = eventCutSelection(GetEventCentralityMaxCutPosition(),1);
   Int_t centCutNumberI             = CutNumberToInteger(centralityCutNumberLow);
   Int_t centCutNumberBI            = CutNumberToInteger(centralityCutNumberHigh);

   if(!option.CompareTo("PbPb_2.76TeV")){    
      if(!centrality.CompareTo("0-5%")){
         fcocktailFunc              = "cocktail_HI0005_oHagBW";//cocktail_HI_0005_QCD";
         fcocktailFuncSyst          = "cocktail_HI_0005_QCD";
         fcocktailFuncNarrow        = "cocktail_HI_0005_QCD";
         fcocktailFuncWide          = "cocktail_HI_0005_QCD";
         fcocktailFuncChargedPions  =  "cocktail_HI_0005_QCD";
         ffive                      = kTRUE;
      }
      if(!centrality.CompareTo("5-10%")){
         fcocktailFunc              = "cocktail_HI0510_qcdBW";//cocktail_HI_0510_QCD";
         fcocktailFuncSyst          = "cocktail_HI_0510_QCD";
         fcocktailFuncNarrow        = "cocktail_HI_0510_QCD";
         fcocktailFuncWide          = "cocktail_HI_0510_QCD";
         fcocktailFuncChargedPions  = "cocktail_HI_0510_QCD";
         ffive                      = kTRUE;
      }
      if(!centrality.CompareTo("0-10%")){
         fcocktailFunc              = "cocktail_PbPb0010_WeightedQGP";
         if(extraSignals.CompareTo("3")) 
            fcocktailFunc           = "cocktail_PbPb0010_WeightedQGP";//cocktail_HI_0010_qcdBWEta";
         fcocktailFuncSyst          = "cocktail_HI_0010_qcdBWEta";
         fcocktailFuncNarrow        = "cocktail_HI_0010_qcdBWEta";
         fcocktailFuncWide          = "cocktail_HI_0010_qcdBWEta";
         fcocktailFuncChargedPions  = "cocktail_HI_0010_qcdBWEta";
      }
      if(!centrality.CompareTo("10-20%")){
         fcocktailFunc              = "cocktail_HI1020_oHagBW";//cocktail_HI_1020_oHagBWEta_CorrCent";
         if(extraSignals.CompareTo("3")) 
            fcocktailFunc           = "cocktail_PbPb1020_WeightedQGP";
         fcocktailFuncSyst          = "cocktail_HI_1020_oHagBWEta_CorrCent";
         fcocktailFuncNarrow        = "cocktail_HI_1020_oHagBWEta_CorrCent";
         fcocktailFuncWide          = "cocktail_HI_1020_oHagBWEta_CorrCent";
         fcocktailFuncChargedPions  = "cocktail_HI_1020_oHagBWEta_CorrCent";
      }
      if(!centrality.CompareTo("0-40%")){
         fcocktailFunc              = "cocktail_HI_00-40_QCD_mtFixBW_9002972094503042212000010420020";
         fcocktailFuncChargedPions  = "cocktail_HI_00-40_QCD_mtFixBW_9002972094503042212000010420020";
         fcocktailFuncSyst          = "cocktail_HI_00-40_QCD_mtFixBW_9002972094503042212000010420020";
         fcocktailFuncNarrow        = "cocktail_HI_00-40_QCD_mtFixBW_9002972094503042212000010420020";
         fcocktailFuncWide          = "cocktail_HI_00-40_QCD_mtFixBW_9002972094503042212000010420020";
      }
      if(!centrality.CompareTo("0-20%")){
         if(FinderType)
             fcocktailFunc          = "cocktail_HI_0020_qcdOffline";
         //fcocktailFunc            = "cocktail_HI_0020_modHagBWEta_HighR"; //9
         //if(extraSignals.CompareTo("3"))
         // fcocktailFunc           = "cocktail_PbPb0020NewMC_MCoHag";//cocktail_PbPb0020NewMC_oHag";
         //fcocktailFunc            = "cocktail_PbPb0020NewMC_MCoHag";//cocktail_PbPb0020NewMC_oHag";
         fcocktailFunc              = "cocktail_PbPb0020NewMC_oHag";
         fcocktailFuncSyst          = "cocktail_HI_0020_modHagBWEta";
         fcocktailFuncNarrow        = "cocktail_HI_0020_modHagBWEta";
         fcocktailFuncWide          = "cocktail_HI_0020_modHagBWEta";
         fcocktailFuncChargedPions  = "cocktail_HI_0020_modHagBWEta";
      }
      if(!centrality.CompareTo("20-40%")){
         fcocktailFuncNarrow        = "cocktail_HI_2040_qcdBWEta";
         fcocktailFuncWide          = "cocktail_HI_2040_qcdBWEta";
         fcocktailFuncSyst          = "cocktail_HI_2040_qcdBWEta";
         fcocktailFunc              = "cocktail_PbPb0020newMC_oHag";//cocktail_HI_2040_qcdBWEta";//cocktail_HI_oHag_12400010420927700237000000_010220450000";//cocktail_HI_2040_qcdBWEta";
         //if(extraSignals.CompareTo("3"))
         //    fcocktailFunc        = "cocktail_PbPb2040_WeightedQGP";
         if(FinderType)
             fcocktailFunc          = "cocktail_HI_2040_qcdOffline";
         fcocktailFuncChargedPions  = "cocktail_HI_2040_qcdBWEta";
      }
      if(!centrality.CompareTo("40-80%")){
         fcocktailFuncNarrow        = "cocktail_HI_40-80_oHag";
         fcocktailFuncWide          = "cocktail_HI_40-80_oHag";
         fcocktailFuncSyst          = "cocktail_HI_40-80_oHag";
         //fcocktailFunc            = "cocktail_HI_4080_qcdpp";//cocktail_PbPb4080NewMC_oHag";//";
         fcocktailFunc              = "cocktail_PbPb4080NewMC_oHag";
         //fcocktailFunc            = "cocktail_PbPb4080NewMC_qcd";
         if(FinderType)
             fcocktailFunc          = "cocktail_HI_4080_qcdOffline";
         fcocktailFuncChargedPions  = "cocktail_HI_40-80_oHag";
      }
      if(!centrality.CompareTo("60-80%")){
         fcocktailFuncNarrow        = "cocktail_HI_40-80_oHag";
         fcocktailFuncWide          = "cocktail_HI_40-80_oHag";
         fcocktailFuncSyst          = "cocktail_HI_40-80_oHag";
         fcocktailFunc              = "cocktail_HI_4080_qcdpp";
         if(FinderType)
             fcocktailFunc          = "cocktail_HI_4080_qcdOffline";
         fcocktailFuncChargedPions  = "cocktail_HI_40-80_oHag";
      }
      if(!centrality.CompareTo("0-80%")){
         fcocktailFuncNarrow        = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncWide          = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncSyst          = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFunc              = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncChargedPions  = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
      }
      if(!centrality.CompareTo("0-90%")){
         fcocktailFuncNarrow        = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncWide          = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncSyst          = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFunc              = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
         fcocktailFuncChargedPions  = "cocktail_HI_0080_qcd";//cocktail_40-80_QCD_Neu_pp";
      }
   
      fcmult                        = 0;
      for(Int_t i = centCutNumberI; i<centCutNumberBI;i++){
         if(ffive){
            fcmult                  += multiplicityCent[i];
            cout<<multiplicityCent[i]<<endl;
         }
         else{
            fcmult                  += multiplicity[i];
            cout<<multiplicity[i]<<endl;
         }
         cout<<" Sum:"<<fcmult<<endl;
      }
      cout<<" / "<<centCutNumberBI-centCutNumberI<<"  =  "<<endl;
      fcmult                        = fcmult/(centCutNumberBI-centCutNumberI);
      cout<<fcmult<<endl;
   }
   else if(!option.CompareTo("pPb_5.023TeV")){
       fcmult = 6.9;
       fcocktailFuncNarrow          = "cocktail_pPb_DPMJET_MB";
       fcocktailFuncWide            = "cocktail_pPb_DPMJET_MB";
       fcocktailFuncSyst            = "cocktail_pPb_DPMJET_MB";
       fcocktailFunc                = "cocktail_pPb_HIJING_MB_pileupCorr_20042014_EtaAndPhimeasured";
       fcocktailFuncModA            = "cocktail_pPb_HIJING_MB_pileupCorr_05042014_etameasured";
       fcocktailFuncModB            = "cocktail_pPb_HIJING_MB_pileupCorr_05042014_etameasured";
       fcocktailFuncMtScaledEta     = "cocktail_pPb_HIJING_MB_pileupCorr";
       fcocktailFuncChargedPions    = "cocktail_pPb_HIJING_MB_pileupCorr_chargedPionsScaledEta";
   }
   else if(!option.CompareTo("900GeV")){
       fcocktailFuncNarrow          = "cocktail_pp2760GeV_oHagmt";
       fcocktailFuncWide            = "cocktail_pp2760GeV_oHagmt";
       fcocktailFuncSyst            = "cocktail_pp2760GeV_oHagmt";
       fcocktailFunc                = "cocktail_pp900GeV_oHag";//cocktail_pp900GeV_Iqcd";
       fcocktailFuncChargedPions    = "cocktail_pp2760GeV_oHagmt";
   }
   else if(!option.CompareTo("2.76TeV")){
       fcocktailFuncNarrow          = "cocktail_pp2760GeV_oHagmt";
       fcocktailFuncWide            = "cocktail_pp2760GeV_oHagmt";
       fcocktailFuncSyst            = "cocktail_pp2760GeV_oHagmt";
       fcocktailFunc                = "cocktail_pp2760GeVEta_Hag";//cocktail_pp2760GeV_oHagmt";//cocktail_pp2760GeV_qcd";//cocktail_pp2760GeV_qcd";//cocktail_pp2760GeV_oHagmt";//cocktail_pp2760GeVmtunfold_oHag_0000013002093663003800000_016310310090";//ocktail_pp2760GeV_oHag_0000013002093663003800000_016310310090";//cocktail_pp2760GeV_oHagmt";//cocktail_pp2760GeV_TsallisTsallis";////";
       fcocktailFuncChargedPions    = "cocktail_pp2760GeV_oHagmt";
   }
   else if(!option.CompareTo("7TeV")){    
      //fcocktailFuncNarrow         = "cocktail_7TeV_pi0LevyetaLevy";
      //fcocktailFuncWide           = "cocktail_7TeV_pi0LevyetaLevy";
      //fcocktailFuncSyst           = "cocktail_7TeV_pi0LevyetaLevy";
       fcocktailFuncNarrow          = "cocktail_pp_7TeV";
       fcocktailFuncWide            = "cocktail_pp_7TeV";
       fcocktailFuncSyst            = "cocktail_pp_7TeV";
      //fcocktailFunc               = "cocktail_7TeV_pi0LevyetaLevy";
      //fcocktailFunc               = "cocktail_pp7TeV_HagLevy";
      //fcocktailFunc               = "cocktail_pp7TeV_HagHag";;
      //fcocktailFunc               = "cocktail_pp7TeV_oHagHag";
      //fcocktailFunc               = "cocktail_pp7TeV_TsallisTsallis";
       fcocktailFunc                = "cocktail_pp_7TeV";
      //fcocktailFunc               = "cocktail_pp7TeV_QCDTsallis";
      //fcocktailFunc               = "cocktail_7TeV_pi0LevyetaLevy";
      //fcocktailFunc               = "cocktail_pp7TeV_qcdqcd";
      //fcocktailFuncChargedPions   = "cocktail_276TeV";
       fcocktailFuncChargedPions    = "";
   }
   else if(!option.CompareTo("13TeV")){             // for testing purposes only
       fcocktailFuncNarrow          = "cocktail_pp_7TeV";
       fcocktailFuncWide            = "cocktail_pp_7TeV";
       fcocktailFuncSyst            = "cocktail_pp_7TeV";
       fcocktailFunc                = "cocktail_pp_7TeV";
       fcocktailFuncChargedPions    = "";
   }
   else {
      fcocktailFuncNarrow           = "cocktail_Test";
      fcocktailFuncWide             = "cocktail_Test";
      fcocktailFuncSyst             = "cocktail_Test";
      fcocktailFunc                 = "cocktail_Test";
      fcocktailFuncChargedPions     =  "cocktail_Test";
   }
}

