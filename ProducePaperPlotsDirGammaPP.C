#include <Riostream.h>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include "TGaxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"
#include "CombineGammaResultsPP.h"


struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //    TString name;
};

void drawhMarker(TH1D* histoDummy, Double_t column, Double_t row, Double_t markerScale){
    TMarker* markerdummy;
    markerdummy= CreateMarkerFromHisto(histoDummy,column ,row ,markerScale);
    markerdummy->DrawMarker(column ,row);
}


void ProducePaperPlotsDirGammaPP(
    TString inputFileName2760           = "CombinationInputGammaPP/CombinedGammaResultPP2760GeV_2017_09_14.root",
    TString inputFileName2760PCM        = "CombinationInputGammaPP/InputPCMGammapp2760GeV.root",
    TString inputFileName2760EMC        = "CombinationInputGammaPP/InputEMCGammapp2760GeV.root",
    TString inputFileName2760PCMEMC     = "CombinationInputGammaPP/InputPCMEMCGammapp2760GeV.root",
    TString inputFileName8000           = "CombinationInputGammaPP/CombinedGammaResultPP8TeV_2017_11_06.root",
    TString inputFileName8000PCM        = "CombinationInputGammaPP/InputPCMGammapp8TeV.root",
    TString inputFileName8000EMC        = "CombinationInputGammaPP/InputEMCGammapp8TeV.root",
    TString inputFileName8000PCMEMC     = "CombinationInputGammaPP/InputPCMEMCGammapp8TeV.root",
    Bool_t plotEMC                      = kTRUE,
    Bool_t plotPCMEMC                   = kTRUE,
    TString suffix                      = "pdf"
){
  //*******************************************************************************************************************************************
  //*********************************************************** Set main style choices ********************************************************
  //*******************************************************************************************************************************************
  StyleSettingsThesis();
  SetPlotStyle();
  
  //*******************************************************************************************************************************************
  //********************************************* Create output directory and copy input files there ******************************************
  //*******************************************************************************************************************************************
  TString dateForOutput                                       = ReturnDateStringForOutput();
  TString outputDir                                           = Form("%s/%s/ProducePaperPlotsDirGammaPP",suffix.Data(),dateForOutput.Data());
  TString fileNameTheorypp                                    = "ExternalInput/Theory/TheoryCompilationPP.root";

  gSystem->Exec("mkdir -p "+outputDir);
  gSystem->Exec(Form("cp %s %s/InputPCMGammapp2760GeV.root",      inputFileName2760.Data(),       outputDir.Data()));
  gSystem->Exec(Form("cp %s %s/InputPCMGammapp2760GeVPCM.root",   inputFileName2760PCM.Data(),    outputDir.Data()));
  gSystem->Exec(Form("cp %s %s/InputPCMGammapp2760GeVEMC.root",   inputFileName2760EMC.Data(),    outputDir.Data()));
  gSystem->Exec(Form("cp %s %s/InputPCMGammapp2760GeVPCMEMC.root",inputFileName2760PCMEMC.Data(), outputDir.Data()));
  gSystem->Exec(Form("cp %s %s/InputPCMGammapp8000GeV.root",      inputFileName8000.Data(),       outputDir.Data()));
  gSystem->Exec(Form("cp %s %s/InputPCMGammapp8000GeVPCM.root",   inputFileName8000PCM.Data(),    outputDir.Data()));
  gSystem->Exec(Form("cp %s %s/InputPCMGammapp8000GeVEMC.root",   inputFileName8000EMC.Data(),    outputDir.Data()));
  gSystem->Exec(Form("cp %s %s/InputPCMGammapp8000GeVPCMEMC.root",inputFileName8000PCMEMC.Data(), outputDir.Data()));
  
  gSystem->Exec(Form("cp %s %s/Theorypp.root",                    fileNameTheorypp.Data(),        outputDir.Data()));
  
  //*******************************************************************************************************************************************
  //******************************************************* set ranges for plotting ***********************************************************
  //*******************************************************************************************************************************************
  Double_t doubleRatio[4];
  Double_t indMeasRatio[4];
  Double_t doubleRatioX[4];
  Double_t doubleRatioXpp[4];
  doubleRatio[0]      = 0.75;     doubleRatio[1]      = 1.65;     doubleRatio[2]      = 1.65;     doubleRatio[3]      = 1.65;
  indMeasRatio[0]     = 0.65;     indMeasRatio[1]     = 1.45;     indMeasRatio[2]     = 1.45;     indMeasRatio[3]     = 1.45;
  doubleRatioX[0]     = 0.7;      doubleRatioX[1]     = 14.5;     doubleRatioX[2]     = 14.5;     doubleRatioX[3]     = 14.5;
  doubleRatioXpp[0]   = 0.23;     doubleRatioXpp[1]   = 14.5;     doubleRatioXpp[2]   = 24.5;     doubleRatioXpp[3]   = 34.5;

  TString  nameMeasGlobal[11]                     = { "PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal",
                                                      "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged",
                                                      "PCMOtherDataset"};
  TString  nameMeasGlobalLabel[11]                = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC",
                                                      "PCM-Dalitz", "PHOS-Dalitz", "EMC-Dalitz", "EMC high pT", "mEMC",
                                                      "PCMOtherDataset"};
  TString  nameMeasGlobalLabelGamma[11]           = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM*",
                                                      "PCM-Dalitz", "PHOS-Dalitz", "EMC-Dalitz", "EMC high pT", "mEMC",
                                                      "PCMOtherDataset"};


  Color_t  colorDet[11];
  Color_t  colorDetMC[11];
  Style_t  markerStyleDet[2][11];
  Style_t  markerStyleDetMC[11];
  Size_t   markerSizeDet[11];
  Size_t   markerSizeDetMC[11];


  for (Int_t i = 0; i < 11; i++){
      colorDet[i]                                     = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
      colorDetMC[i]                                   = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
      markerStyleDet[0][i]                               = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
      markerStyleDet[1][i]                               = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
      markerStyleDetMC[i]                             = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
      markerSizeDet[i]                                = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
      markerSizeDetMC[i]                              = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
  }
  markerStyleDet[1][0]                                = 24;
  markerStyleDet[1][2]                                = 27;
  markerStyleDet[1][4]                                = 28;
  Color_t colorComb2760                               = GetColorDefaultColor("2.76TeV", "", "");
  Style_t markerStyleComb2760                         = GetDefaultMarkerStyle("2.76TeV", "", "");
  Size_t markerSizeComb2760                           = GetDefaultMarkerSize("2.76TeV", "", "");
  Color_t colorComb8000                               = GetColorDefaultColor("8TeV", "", "");
  Style_t markerStyleComb8000                         = GetDefaultMarkerStyle("8TeV", "", "");
  Size_t markerSizeComb8000                           = GetDefaultMarkerSize("8TeV", "", "");
  Color_t defaultColor                                = kBlack;
  Color_t defaultColorBox                             = kGray+2;
  Style_t defaultMarker                               = 20;
  Size_t defaultMarkerSize                            = 1.8;

  Width_t widthLinesBoxes                             = 1.4;
  Width_t widthCommonFit                              = 2.4;

  TString collisionSystempp2760GeV                    = "pp #sqrt{#it{s}} = 2.76 TeV";
  TString collisionSystempp8000GeV                    = "pp #sqrt{#it{s}} = 8 TeV";
  TString textALICE                                   = "ALICE";
  
  TGraphAsymmErrors*  graphCombDRSys[2]               = {NULL, NULL};
  TGraphAsymmErrors*  graphCombDRStatPlot[2]          = {NULL, NULL};
  TGraphAsymmErrors*  graphCombIncGammaTot[2]         = {NULL, NULL};
  TGraphAsymmErrors*  graphCombIncGammaStat[2]        = {NULL, NULL};
  TGraphAsymmErrors*  graphCombIncGammaSys[2]         = {NULL, NULL};
  TGraphAsymmErrors*  graphCombDirGammaSpectrumSystErr[2]         = {NULL, NULL};
  TGraphAsymmErrors*  graphCombDirGammaSpectrumStatErr[2]        = {NULL, NULL};
  TGraphAsymmErrors*  graphCombDirGammaSpectrumSumErrAr[2]         = {NULL, NULL};
  
  TH1D* histoDRPi0FitStatErr[2][11]                  = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TGraphAsymmErrors*  graphDRPi0FitSysErr[2][11]     = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TH1D* histoIncGammaRatioStatErr[2][11]             = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TGraphAsymmErrors* graphIncGammaRatioSysErr[2][11] = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TH1D* histoIncGammaStatErr[2][11]                  = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TGraphAsymmErrors* graphIncGammaStatErr[2][11]     = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TGraphAsymmErrors* graphIncGammaSysErr[2][11]      = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TH1D* histoConvProb[2][11]                         = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TH1D* histoEffi[2][11]                             = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TH1D* histoEffiMCPt[2][11]                         = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TH1D* histoPurity[2][11]                           = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TH1D* histoResolCorr[2][11]                        = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  TH1D* histoEffSecCorr[2][4][11]                    = {{ {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                      {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                      {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                      {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL} },
                                                      { {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                      {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                      {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                      {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL} }};
  TH1D* histoPileupCorr[2][11]                       = {{NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
  //*******************************************************************************************************************************************
  //*********************************************** Load Comb histograms from 2.76 file *******************************************************
  //*******************************************************************************************************************************************
  TFile* fileCombGammapp2760GeV                         = new TFile( inputFileName2760.Data());
  graphCombDRStatPlot[0]                                = (TGraphAsymmErrors*) fileCombGammapp2760GeV->Get("Gamma2.76TeV/graphRGammaCombStatErr");
  graphCombDRSys[0]                                     = (TGraphAsymmErrors*) fileCombGammapp2760GeV->Get("Gamma2.76TeV/graphRGammaCombSysErr");
  graphCombIncGammaTot[0]                               = (TGraphAsymmErrors*) fileCombGammapp2760GeV->Get("Gamma2.76TeV/graphInvYieldIncGammaTotErr");
  graphCombIncGammaStat[0]                              = (TGraphAsymmErrors*) fileCombGammapp2760GeV->Get("Gamma2.76TeV/graphInvYieldIncGammaStatErr");
  graphCombIncGammaSys[0]                               = (TGraphAsymmErrors*) fileCombGammapp2760GeV->Get("Gamma2.76TeV/graphInvYieldIncGammaSysErr");
  graphCombDirGammaSpectrumSystErr[0]                   = (TGraphAsymmErrors*) fileCombGammapp2760GeV->Get("Gamma2.76TeV/graphInvYieldDirGammaSysErr");
  graphCombDirGammaSpectrumStatErr[0]                   = (TGraphAsymmErrors*) fileCombGammapp2760GeV->Get("Gamma2.76TeV/graphInvYieldDirGammaStatErr");
  graphCombDirGammaSpectrumSumErrAr[0]                  = (TGraphAsymmErrors*) fileCombGammapp2760GeV->Get("Gamma2.76TeV/graphInvYieldDirGammaSumErrAr");
  
  TFile* fileCombGammapp8TeV                            = new TFile( inputFileName8000.Data());
  graphCombDRStatPlot[1]                                = (TGraphAsymmErrors*) fileCombGammapp8TeV->Get("Gamma8TeV/graphRGammaCombStatErr");
  graphCombDRSys[1]                                     = (TGraphAsymmErrors*) fileCombGammapp8TeV->Get("Gamma8TeV/graphRGammaCombSysErr");
  graphCombIncGammaTot[1]                               = (TGraphAsymmErrors*) fileCombGammapp8TeV->Get("Gamma8TeV/graphInvYieldIncGammaTotErr");
  graphCombIncGammaStat[1]                              = (TGraphAsymmErrors*) fileCombGammapp8TeV->Get("Gamma8TeV/graphInvYieldIncGammaStatErr");
  graphCombIncGammaSys[1]                               = (TGraphAsymmErrors*) fileCombGammapp8TeV->Get("Gamma8TeV/graphInvYieldIncGammaSysErr");
  graphCombDirGammaSpectrumSystErr[1]                   = (TGraphAsymmErrors*) fileCombGammapp8TeV->Get("Gamma8TeV/graphInvYieldDirGammaSysErr");
  graphCombDirGammaSpectrumStatErr[1]                   = (TGraphAsymmErrors*) fileCombGammapp8TeV->Get("Gamma8TeV/graphInvYieldDirGammaStatErr");
  graphCombDirGammaSpectrumSumErrAr[1]                  = (TGraphAsymmErrors*) fileCombGammapp8TeV->Get("Gamma8TeV/graphInvYieldDirGammaSumErrAr");
  

  //*******************************************************************************************************************************************
  //*********************************************** Load PCM histograms from PCM file ********************************************************
  //*******************************************************************************************************************************************
  TFile* filePCMGammapp2760GeV                    = new TFile( inputFileName2760PCM.Data());
  //________________________________________________ Load PCM 2.76TeV _________________________________________________________________________
  TDirectory* directoryPCMGammapp2760GeV          = (TDirectory*)filePCMGammapp2760GeV->Get("Gamma_pp2760GeV");
      histoDRPi0FitStatErr[0][0]                         = (TH1D*) directoryPCMGammapp2760GeV->Get("DoubleRatioPi0FitStatError");
      graphDRPi0FitSysErr[0][0]                          = (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("DoubleRatioPi0FitSystError");
      histoIncGammaRatioStatErr[0][0]                    = (TH1D*) directoryPCMGammapp2760GeV->Get("IncRatioStatError");
      graphIncGammaRatioSysErr[0][0]                     = (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("IncRatioSystError");
      histoConvProb[0][0]                                = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaConversionProbability");
      histoEffi[0][0]                                    = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaRecoEfficiency");
      histoEffiMCPt[0][0]                                = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaRecoEfficiencyMCPt");
      histoPurity[0][0]                                  = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaTruePurity");
      histoResolCorr[0][0]                               = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaResolCorr");
      histoIncGammaStatErr[0][0]                         = (TH1D*) directoryPCMGammapp2760GeV->Get("IncGammaStatError");
      graphIncGammaStatErr[0][0]                         = new TGraphAsymmErrors(histoIncGammaStatErr[0][0]);
      graphIncGammaSysErr[0][0]                          = (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("IncGammaSystError");
      histoPileupCorr[0][0]                              = (TH1D*) directoryPCMGammapp2760GeV->Get("PileUpCorrectionFactor");
      histoEffSecCorr[0][0][0]                           = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0s");
      histoEffSecCorr[0][1][0]                           = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0l");
      histoEffSecCorr[0][2][0]                           = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Lambda");
      histoEffSecCorr[0][3][0]                           = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Rest");
      
  TFile* filePCMGammapp8TeV                    = new TFile( inputFileName8000PCM.Data());
  //________________________________________________ Load PCM 2.76TeV _________________________________________________________________________
  TDirectory* directoryPCMGammapp8TeV          = (TDirectory*)filePCMGammapp8TeV->Get("Gamma_pp8TeV");
      histoDRPi0FitStatErr[1][0]                         = (TH1D*) directoryPCMGammapp8TeV->Get("DoubleRatioPi0FitStatError");
      graphDRPi0FitSysErr[1][0]                          = (TGraphAsymmErrors*) directoryPCMGammapp8TeV->Get("DoubleRatioPi0FitSystError");
      histoIncGammaRatioStatErr[1][0]                    = (TH1D*) directoryPCMGammapp8TeV->Get("IncRatioStatError");
      graphIncGammaRatioSysErr[1][0]                     = (TGraphAsymmErrors*) directoryPCMGammapp8TeV->Get("IncRatioSystError");
      histoConvProb[1][0]                                = (TH1D*) directoryPCMGammapp8TeV->Get("GammaConversionProbability");
      histoEffi[1][0]                                    = (TH1D*) directoryPCMGammapp8TeV->Get("GammaRecoEfficiency");
      histoEffiMCPt[1][0]                                = (TH1D*) directoryPCMGammapp8TeV->Get("GammaRecoEfficiencyMCPt");
      histoPurity[1][0]                                  = (TH1D*) directoryPCMGammapp8TeV->Get("GammaTruePurity");
      histoResolCorr[1][0]                               = (TH1D*) directoryPCMGammapp8TeV->Get("GammaResolCorr");
      histoIncGammaStatErr[1][0]                         = (TH1D*) directoryPCMGammapp8TeV->Get("IncGammaStatError");
      graphIncGammaStatErr[1][0]                         = new TGraphAsymmErrors(histoIncGammaStatErr[1][0]);
      graphIncGammaSysErr[1][0]                          = (TGraphAsymmErrors*) directoryPCMGammapp8TeV->Get("IncGammaSystError");
      histoPileupCorr[1][0]                              = (TH1D*) directoryPCMGammapp8TeV->Get("PileUpCorrectionFactor");
      histoEffSecCorr[1][0][0]                           = (TH1D*) directoryPCMGammapp8TeV->Get("GammaEffectiveSecondaryCorr_K0s");
      histoEffSecCorr[1][1][0]                           = (TH1D*) directoryPCMGammapp8TeV->Get("GammaEffectiveSecondaryCorr_K0l");
      histoEffSecCorr[1][2][0]                           = (TH1D*) directoryPCMGammapp8TeV->Get("GammaEffectiveSecondaryCorr_Lambda");
      histoEffSecCorr[1][3][0]                           = (TH1D*) directoryPCMGammapp8TeV->Get("GammaEffectiveSecondaryCorr_Rest");
      
      //*******************************************************************************************************************************************
      //*********************************************** Load PCMEMC histograms from 2.76TeV PCM file **********************************************
      //*******************************************************************************************************************************************
      if (plotPCMEMC){
          TFile* filePCMEMCGammapp2760GeV                 = new TFile( inputFileName2760PCMEMC.Data());
          //________________________________________________ Load PCM-EMC 2.76TeV _________________________________________________________________________
          TDirectory* directoryPCMEMCGammapp2760GeV       = (TDirectory*)filePCMEMCGammapp2760GeV->Get("Gamma_pp2760GeV");
              histoDRPi0FitStatErr[0][4]                         = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("DoubleRatioPi0FitStatError");
              graphDRPi0FitSysErr[0][4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapp2760GeV->Get("DoubleRatioPi0FitSystError");
              histoIncGammaRatioStatErr[0][4]                    = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("IncRatioStatError");
              graphIncGammaRatioSysErr[0][4]                     = (TGraphAsymmErrors*) directoryPCMEMCGammapp2760GeV->Get("IncRatioSystError");
              histoConvProb[0][4]                                = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaConversionProbability");
              histoEffi[0][4]                                    = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaRecoEfficiency");
              histoEffiMCPt[0][4]                                = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaRecoEfficiencyMCPt");
              histoPurity[0][4]                                  = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaTruePurity");
              histoResolCorr[0][4]                               = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaResolCorr");
              histoIncGammaStatErr[0][4]                         = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("IncGammaStatError");
              graphIncGammaStatErr[0][4]                         = new TGraphAsymmErrors(histoIncGammaStatErr[0][4]);
              graphIncGammaSysErr[0][4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapp2760GeV->Get("IncGammaSystError");
              histoPileupCorr[0][4]                              = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("PileUpCorrectionFactor");
              if (histoPileupCorr[0][4]) histoPileupCorr[0][4]->GetXaxis()->SetRangeUser(0.8,8);
              histoEffSecCorr[0][0][4]                           = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0s");
              histoEffSecCorr[0][1][4]                           = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0l");
              histoEffSecCorr[0][2][4]                           = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Lambda");
              histoEffSecCorr[0][3][4]                           = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Rest");
              
          TFile* filePCMEMCGammapp8000GeV                 = new TFile( inputFileName8000PCMEMC.Data());
          //________________________________________________ Load PCM-EMC 2.76TeV _________________________________________________________________________
          TDirectory* directoryPCMEMCGammapp8000GeV       = (TDirectory*)filePCMEMCGammapp8000GeV->Get("Gamma_pp8TeV");
              histoDRPi0FitStatErr[1][4]                         = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("DoubleRatioPi0FitStatError");
              graphDRPi0FitSysErr[1][4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapp8000GeV->Get("DoubleRatioPi0FitSystError");
              histoIncGammaRatioStatErr[1][4]                    = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("IncRatioStatError");
              graphIncGammaRatioSysErr[1][4]                     = (TGraphAsymmErrors*) directoryPCMEMCGammapp8000GeV->Get("IncRatioSystError");
              histoConvProb[1][4]                                = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("GammaConversionProbability");
              histoEffi[1][4]                                    = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("GammaRecoEfficiency");
              histoEffiMCPt[1][4]                                = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("GammaRecoEfficiencyMCPt");
              histoPurity[1][4]                                  = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("GammaTruePurity");
              histoResolCorr[1][4]                               = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("GammaResolCorr");
              histoIncGammaStatErr[1][4]                         = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("IncGammaStatError");
              graphIncGammaStatErr[1][4]                         = new TGraphAsymmErrors(histoIncGammaStatErr[1][4]);
              graphIncGammaSysErr[1][4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapp8000GeV->Get("IncGammaSystError");
              histoPileupCorr[1][4]                              = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("PileUpCorrectionFactor");
              if (histoPileupCorr[1][4]) histoPileupCorr[1][4]->GetXaxis()->SetRangeUser(0.8,8);
              histoEffSecCorr[1][0][4]                           = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("GammaEffectiveSecondaryCorr_K0s");
              histoEffSecCorr[1][1][4]                           = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("GammaEffectiveSecondaryCorr_K0l");
              histoEffSecCorr[1][2][4]                           = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("GammaEffectiveSecondaryCorr_Lambda");
              histoEffSecCorr[1][3][4]                           = (TH1D*) directoryPCMEMCGammapp8000GeV->Get("GammaEffectiveSecondaryCorr_Rest");
      }
      //*******************************************************************************************************************************************
      //*********************************************** Load PCM histograms from 2.76TeV EMC file *************************************************
      //*******************************************************************************************************************************************
      if (plotEMC){
          TFile* fileEMCGammapp2760GeV                    = new TFile( inputFileName2760EMC.Data());
          //________________________________________________ Load EMC 2.76TeV _________________________________________________________________________
          TDirectory* directoryEMCGammapp2760GeV          = (TDirectory*)fileEMCGammapp2760GeV->Get("Gamma_pp2760GeV");
              histoDRPi0FitStatErr[0][2]                         = (TH1D*) directoryEMCGammapp2760GeV->Get("DoubleRatioPi0FitStatError");
              graphDRPi0FitSysErr[0][2]                          = (TGraphAsymmErrors*) directoryEMCGammapp2760GeV->Get("DoubleRatioPi0FitSystError");
              histoIncGammaRatioStatErr[0][2]                    = (TH1D*) directoryEMCGammapp2760GeV->Get("IncRatioStatError");
              graphIncGammaRatioSysErr[0][2]                     = (TGraphAsymmErrors*) directoryEMCGammapp2760GeV->Get("IncRatioSystError");
              histoEffi[0][2]                                    = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaRecoEfficiency");
              histoEffiMCPt[0][2]                                = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaRecoEfficiencyMCPt");
              histoPurity[0][2]                                  = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaTruePurity");
              histoResolCorr[0][2]                               = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaResolCorr");
              histoIncGammaStatErr[0][2]                         = (TH1D*) directoryEMCGammapp2760GeV->Get("IncGammaStatError");
              graphIncGammaStatErr[0][2]                         = new TGraphAsymmErrors(histoIncGammaStatErr[0][2]);
              graphIncGammaSysErr[0][2]                          = (TGraphAsymmErrors*) directoryEMCGammapp2760GeV->Get("IncGammaSystError");
              histoEffSecCorr[0][0][2]                           = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0s");
              histoEffSecCorr[0][1][2]                           = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0l");
              histoEffSecCorr[0][2][2]                           = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Lambda");
              histoEffSecCorr[0][3][2]                           = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Rest");
              
          TFile* fileEMCGammapp8000GeV                    = new TFile( inputFileName8000EMC.Data());
          //________________________________________________ Load EMC 2.76TeV _________________________________________________________________________
          TDirectory* directoryEMCGammapp8000GeV          = (TDirectory*)fileEMCGammapp8000GeV->Get("Gamma_pp8TeV");
              histoDRPi0FitStatErr[1][2]                         = (TH1D*) directoryEMCGammapp8000GeV->Get("DoubleRatioPi0FitStatError");
              graphDRPi0FitSysErr[1][2]                          = (TGraphAsymmErrors*) directoryEMCGammapp8000GeV->Get("DoubleRatioPi0FitSystError");
              histoIncGammaRatioStatErr[1][2]                    = (TH1D*) directoryEMCGammapp8000GeV->Get("IncRatioStatError");
              graphIncGammaRatioSysErr[1][2]                     = (TGraphAsymmErrors*) directoryEMCGammapp8000GeV->Get("IncRatioSystError");
              histoEffi[1][2]                                    = (TH1D*) directoryEMCGammapp8000GeV->Get("GammaRecoEfficiency");
              histoEffiMCPt[1][2]                                = (TH1D*) directoryEMCGammapp8000GeV->Get("GammaRecoEfficiencyMCPt");
              histoPurity[1][2]                                  = (TH1D*) directoryEMCGammapp8000GeV->Get("GammaTruePurity");
              histoResolCorr[1][2]                               = (TH1D*) directoryEMCGammapp8000GeV->Get("GammaResolCorr");
              histoIncGammaStatErr[1][2]                         = (TH1D*) directoryEMCGammapp8000GeV->Get("IncGammaStatError");
              graphIncGammaStatErr[1][2]                         = new TGraphAsymmErrors(histoIncGammaStatErr[1][2]);
              graphIncGammaSysErr[1][2]                          = (TGraphAsymmErrors*) directoryEMCGammapp8000GeV->Get("IncGammaSystError");
              histoEffSecCorr[1][0][2]                           = (TH1D*) directoryEMCGammapp8000GeV->Get("GammaEffectiveSecondaryCorr_K0s");
              histoEffSecCorr[1][1][2]                           = (TH1D*) directoryEMCGammapp8000GeV->Get("GammaEffectiveSecondaryCorr_K0l");
              histoEffSecCorr[1][2][2]                           = (TH1D*) directoryEMCGammapp8000GeV->Get("GammaEffectiveSecondaryCorr_Lambda");
              histoEffSecCorr[1][3][2]                           = (TH1D*) directoryEMCGammapp8000GeV->Get("GammaEffectiveSecondaryCorr_Rest");

      }
      //*******************************************************************************************************************************************
      //************************************************ Load theory curves from external input ***************************************************
      //*******************************************************************************************************************************************
      TFile* fileTheory                               = new TFile( fileNameTheorypp.Data());
          TGraphAsymmErrors* graphTheoryNLODRpp[2];
          TGraph* graphTheoryNLODRppCenter[2];
          TGraphAsymmErrors* graphTheoryNLOpp[2];
          TGraph* graphTheoryNLODRppPaquettCenter[2];
          TGraph* graphTheoryNLOppPaquettCenter[2];
          
          graphTheoryNLODRpp[0]               = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonNLOVogelsangInvYieldINT1_pp2760GeV_CT10_ALICECocktail");
          graphTheoryNLODRppCenter[0]         = (TGraph*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonNLOVogelsangInvYieldINT1_pp2760GeV_CT10_ALICECocktail_Center");
          graphTheoryNLOpp[0]                 = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphDirectPhotonNLOVogelsangInvYieldINT1_2760GeV");
          graphTheoryNLODRppPaquettCenter[0]  = (TGraph*) fileTheory->Get("DirectPhoton/graphDRNLOPaquett_2760GeV_ALICECocktail");
          graphTheoryNLOppPaquettCenter[0]    = (TGraph*) fileTheory->Get("DirectPhoton/graphNLOPaquett_2760GeV");
          
          graphTheoryNLODRpp[1]               = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonNLOVogelsangInvYieldINT1_pp2760GeV_CT10_ALICECocktail");
          graphTheoryNLODRppCenter[1]         = (TGraph*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonNLOVogelsangInvYieldINT1_pp2760GeV_CT10_ALICECocktail_Center");
          graphTheoryNLOpp[1]                 = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphDirectPhotonNLOVogelsangInvYieldINT1_2760GeV");
          graphTheoryNLODRppPaquettCenter[1]  = (TGraph*) fileTheory->Get("DirectPhoton/graphDRNLOPaquett_2760GeV_ALICECocktail");
          graphTheoryNLOppPaquettCenter[1]    = (TGraph*) fileTheory->Get("DirectPhoton/graphNLOPaquett_2760GeV");

          TF1* fitHagGammaComb[2];
          fitHagGammaComb[0]                  = FitObject("h","fitHagGammaComb2760","Gamma0",graphCombIncGammaTot[0],graphCombIncGammaTot[0]->GetX()[0], graphCombIncGammaTot[0]->GetX()[graphCombIncGammaTot[0]->GetN()],NULL,"QNRME+");
          // fitHagGammaComb[1]               = FitObject("h","fitHagGammaComb8000","Gamma1",graphCombIncGammaStat[1],graphCombIncGammaStat[1]->GetX()[0], graphCombIncGammaStat[1]->GetX()[graphCombIncGammaStat[1]->GetN()],NULL,"QNRME+");
          fitHagGammaComb[1]                  = FitObject("oHag","fitHagGammaComb8000","Gamma1",histoIncGammaStatErr[1][0],histoIncGammaStatErr[1][0]->GetXaxis()->GetBinLowEdge(1), histoIncGammaStatErr[1][4]->GetXaxis()->GetBinUpEdge(histoIncGammaStatErr[1][0]->GetNbinsX()),NULL,"QNRME+");
          
          TGraphAsymmErrors* graphRatioGammaCombCombFitTot[2]     = { NULL, NULL};
          TGraphAsymmErrors* graphRatioGammaCombCombFitStat[2]    = { NULL, NULL};
          TGraphAsymmErrors* graphRatioGammaCombCombFitSys[2]     = { NULL, NULL};
          TGraphAsymmErrors* graphRatioGammaIndCombFitStat[2][11] = {{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                     { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};
          TGraphAsymmErrors* graphRatioGammaIndCombFitSys[2][11]  = {{ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                     { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

          for (Int_t j= 0; j< 2; j++){
            graphRatioGammaCombCombFitTot[j]     = (TGraphAsymmErrors*)graphCombIncGammaTot[j]->Clone();
            graphRatioGammaCombCombFitTot[j]     = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitTot[j], fitHagGammaComb[j]);
            graphRatioGammaCombCombFitStat[j]    = (TGraphAsymmErrors*)graphCombIncGammaStat[j]->Clone();
            graphRatioGammaCombCombFitStat[j]    = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitStat[j], fitHagGammaComb[j]);
            graphRatioGammaCombCombFitSys[j]     = (TGraphAsymmErrors*)graphCombIncGammaSys[j]->Clone();
            graphRatioGammaCombCombFitSys[j]     = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitSys[j], fitHagGammaComb[j]);

            for (Int_t i= 0; i< 11; i++){
                if (graphIncGammaStatErr[j][i]){
                  graphRatioGammaIndCombFitStat[j][i]             = (TGraphAsymmErrors*)graphIncGammaStatErr[j][i]->Clone(Form("RatioGamma%sToCombFitStat", nameMeasGlobalLabel[i].Data()));
                  graphRatioGammaIndCombFitStat[j][i]             = CalculateGraphErrRatioToFit(graphRatioGammaIndCombFitStat[j][i], fitHagGammaComb[j]);
                }
              if (graphIncGammaSysErr[j][i]){
                  graphRatioGammaIndCombFitSys[j][i]              = (TGraphAsymmErrors*)graphIncGammaSysErr[j][i]->Clone(Form("RatioGamma%sToCombFitSyst", nameMeasGlobalLabel[i].Data()));
                  graphRatioGammaIndCombFitSys[j][i]              = CalculateGraphErrRatioToFit(graphRatioGammaIndCombFitSys[j][i], fitHagGammaComb[j]);
              }
            }
          }











  
          Double_t textSizeLabelsPixel                    = 50;
          Double_t textSizeLabelsRel                      = 50./1200;
          // **********************************************************************************************************************
          // ******************************** InvYield for indiv. inc gamma measurement *****************************************
          // **********************************************************************************************************************
            
            TCanvas* canvasInvYieldGamma          = new TCanvas("canvasInvYieldGamma","",200,10,1350,1350*1.15);  // gives the page size
            DrawGammaCanvasSettings( canvasInvYieldGamma, 0.16, 0.02, 0.02, 0.08);
            canvasInvYieldGamma->SetLogx();
            canvasInvYieldGamma->SetLogy();


            TH1F * histo2DYieldGamma              = new TH1F("histo2DYieldGamma","histo2DYieldGamma",11000,doubleRatioXpp[0], doubleRatioXpp[2]);
            SetStyleHistoTH1ForGraphs(histo2DYieldGamma, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.7);
            histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-9,3.5e1);
            histo2DYieldGamma->GetXaxis()->SetMoreLogLabels();
            histo2DYieldGamma->GetXaxis()->SetLabelOffset(-0.01);
            histo2DYieldGamma->Draw("copy");
            
            TLegend* legendYieldIncGammaInd2760       = GetAndSetLegend2(0.20, 0.12, 0.5, 0.12+(3*textSizeLabelsRel*0.85),textSizeLabelsPixel);
            // TLegend* legendYieldIncGammaInd8000       = GetAndSetLegend2(0.62, 0.83-(3*textSizeLabelsRel*0.85), 0.95, 0.83,textSizeLabelsPixel);
            TLegend* legendYieldIncGammaInd8000       = GetAndSetLegend2(0.72, 0.83-(3*textSizeLabelsRel*0.85), 0.95, 0.83,textSizeLabelsPixel);
            for (Int_t i = 0; i < 2; i++){
              for (Int_t j = 0; j < 11; j++){
                if(graphIncGammaSysErr[i][j]){
                  DrawGammaSetMarkerTGraphAsym(graphIncGammaSysErr[i][j], markerStyleDet[i][j], markerSizeDet[j], colorDet[j] , colorDet[j],widthLinesBoxes, kTRUE);
                  graphIncGammaSysErr[i][j]->Draw("E2same");
                  DrawGammaSetMarkerTGraphAsym(graphIncGammaStatErr[i][j],  markerStyleDet[i][j], markerSizeDet[j], colorDet[j] , colorDet[j]);
                  graphIncGammaStatErr[i][j]->Draw("Epsame,x0");
                  if(i==0)
                    legendYieldIncGammaInd2760->AddEntry(graphIncGammaSysErr[i][j], Form("#gamma_{inc} %s", nameMeasGlobalLabelGamma[j].Data()),"pf");
                  else if(i==1)
                    legendYieldIncGammaInd8000->AddEntry(graphIncGammaSysErr[i][j], Form("#gamma_{inc} %s", nameMeasGlobalLabelGamma[j].Data()),"pf");
                }
              }
            }
            legendYieldIncGammaInd2760->Draw();
            legendYieldIncGammaInd8000->Draw();
            
            TLatex *labelEnergyInvYieldPaperAll8000 = new TLatex(0.72, 0.92, "pp #sqrt{#it{s}} = 8 TeV");
            SetStyleTLatex( labelEnergyInvYieldPaperAll8000, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
            labelEnergyInvYieldPaperAll8000->Draw();
            TLatex *labelALICEInvYieldPaperAll8000  = new TLatex(0.72,0.92-0.04,"ALICE");
            SetStyleTLatex( labelALICEInvYieldPaperAll8000, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
            labelALICEInvYieldPaperAll8000->Draw();
            TLatex *labelALICENormUnPaperAll8000    = new TLatex(0.72,0.92-0.07,"Norm. unc. 2.6%");
            SetStyleTLatex( labelALICENormUnPaperAll8000, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
            labelALICENormUnPaperAll8000->Draw();
            
            TLatex *labelEnergyInvYieldPaperAll = new TLatex(0.20, 0.20+0.04*3, "pp #sqrt{#it{s}} = 2.76 TeV");
            SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
            labelEnergyInvYieldPaperAll->Draw();
            TLatex *labelALICEInvYieldPaperAll  = new TLatex(0.20,0.20+0.04*2,"ALICE");
            SetStyleTLatex( labelALICEInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
            labelALICEInvYieldPaperAll->Draw();
            TLatex *labelALICENormUnPaperAll    = new TLatex(0.20,0.20+0.05*1,"Norm. unc. 2.5%");
            SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
            labelALICENormUnPaperAll->Draw();
            

            histo2DYieldGamma->Draw("same,axis");
            canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_IndMeas.%s",outputDir.Data(),suffix.Data()));
            
            
          
          
            histo2DYieldGamma->Draw("copy");
            DrawGammaSetMarkerTF1( fitHagGammaComb[0], 3, 2.5, kGray+3);
            DrawGammaSetMarkerTF1( fitHagGammaComb[1], 2, 2.5, kGray+1);
            fitHagGammaComb[0]->SetRange(0.38,10.5);
            fitHagGammaComb[1]->SetRange(0.28,16.5);
            TLegend* legendYieldIncGammaInd2760_comb       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.12+(2*textSizeLabelsRel*0.85),textSizeLabelsPixel);
            TLegend* legendYieldIncGammaInd8000_comb       = GetAndSetLegend2(0.64, 0.86-(2*textSizeLabelsRel*0.85), 0.95, 0.87,textSizeLabelsPixel);
            for (Int_t i = 0; i < 2; i++){
                if(graphIncGammaSysErr[i]){
                  DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys[i], markerStyleDet[i][0], markerSizeDet[i], defaultColor , defaultColor,widthLinesBoxes, kTRUE);
                  graphCombIncGammaSys[i]->Draw("E2same");
                  DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStat[i],  markerStyleDet[i][0], markerSizeDet[i], defaultColor , defaultColor);
                  graphCombIncGammaStat[i]->Draw("Epsame,x0");
                  if(i==0){
                    legendYieldIncGammaInd2760_comb->AddEntry(graphCombIncGammaSys[i], Form("#gamma_{inc} %s", "ALICE"),"pf");
                    legendYieldIncGammaInd2760_comb->AddEntry(fitHagGammaComb[i], "Hagedorn fit","l");
                  }
                  else if(i==1){
                    legendYieldIncGammaInd8000_comb->AddEntry(graphCombIncGammaSys[i], Form("#gamma_{inc} %s", "ALICE"),"pf");
                    legendYieldIncGammaInd8000_comb->AddEntry(fitHagGammaComb[i], "mod. Hag. fit","l");
                  }
                  fitHagGammaComb[i]->Draw("same");
                }
                
            }
            
            legendYieldIncGammaInd2760_comb->Draw();
            legendYieldIncGammaInd8000_comb->Draw();
            
            labelEnergyInvYieldPaperAll8000 = new TLatex(0.64, 0.92, "pp #sqrt{#it{s}} = 8 TeV");
            SetStyleTLatex( labelEnergyInvYieldPaperAll8000, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
            labelEnergyInvYieldPaperAll8000->Draw();
            // labelALICEInvYieldPaperAll8000  = new TLatex(0.62,0.92-0.04,"ALICE");
            // SetStyleTLatex( labelALICEInvYieldPaperAll8000, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
            // labelALICEInvYieldPaperAll8000->Draw();
            labelALICENormUnPaperAll8000    = new TLatex(0.64,0.92-0.03,"Norm. unc. 2.6%");
            SetStyleTLatex( labelALICENormUnPaperAll8000, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
            labelALICENormUnPaperAll8000->Draw();
            
            labelEnergyInvYieldPaperAll = new TLatex(0.20, 0.12+0.04*3, "pp #sqrt{#it{s}} = 2.76 TeV");
            SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
            labelEnergyInvYieldPaperAll->Draw();
            // labelALICEInvYieldPaperAll  = new TLatex(0.20,0.12+0.04*2,"ALICE");
            // SetStyleTLatex( labelALICEInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
            // labelALICEInvYieldPaperAll->Draw();
            labelALICENormUnPaperAll    = new TLatex(0.20,0.12+0.05+0.04,"Norm. unc. 2.5%");
            SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
            labelALICENormUnPaperAll->Draw();

            histo2DYieldGamma->Draw("same,axis");
            canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_Comb.%s",outputDir.Data(),suffix.Data()));
            
            
  // **********************************************************************************************************************
  // ******************************** InvYield for combined inc gamma measurement *****************************************
  // **********************************************************************************************************************

        TCanvas *canvasDoubleRatio = new TCanvas("canvasDoubleRatio","",0.095,0.09,1000,815);
        DrawGammaCanvasSettings( canvasDoubleRatio, 0.086, 0.01, 0.01, 0.105);
        canvasDoubleRatio->cd();
        canvasDoubleRatio->SetLogx();

        Double_t minY                            = 0.85;
        Double_t maxY                            = 1.65;
        Double_t textSizeSinglePad               = 0.05;
        TH2F * hist2DDRDummySingle       = new TH2F("hist2DDRDummySingle","hist2DDRDummySingle",1000,0.01, 50,1000,doubleRatio[0], doubleRatio[1]);
        SetStyleHistoTH2ForGraphs(hist2DDRDummySingle, "#it{p}_{T} (GeV/#it{c})","#it{R}_{#gamma}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
        hist2DDRDummySingle->GetXaxis()->SetRangeUser(doubleRatioXpp[0], doubleRatioXpp[1]);
        hist2DDRDummySingle->GetXaxis()->SetLabelOffset(-0.01);
        hist2DDRDummySingle->GetXaxis()->SetMoreLogLabels(kTRUE);
        hist2DDRDummySingle->DrawCopy();
        TLegend* legendDRSingle = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*2,0.5,0.953, textSizeSinglePad, 2, "", 42, 0.3);
        legendDRSingle->SetTextAlign(11);

        TLatex *labelDRSingle[2];
        labelDRSingle[0] = new TLatex(0.95,0.92,Form("%s",collisionSystempp2760GeV.Data()));
        SetStyleTLatex( labelDRSingle[0], textSizeSinglePad,4, 1, 42, kTRUE, 31);
        labelDRSingle[1] = new TLatex(0.95,0.92,Form("%s",collisionSystempp8000GeV.Data()));
        SetStyleTLatex( labelDRSingle[1], textSizeSinglePad,4, 1, 42, kTRUE, 31);
        TLatex *labelALICEDRSingle = new TLatex(0.95,0.87,"ALICE");
        SetStyleTLatex( labelALICEDRSingle, textSizeSinglePad,4, 1, 42, kTRUE, 31);

    for (Int_t j = 0; j < 2; j++){
      // if(j==1)
        // legendDRSingle = GetAndSetLegend2(0.12,0.953-textSizeSinglePad,0.5,0.953, textSizeSinglePad, 2, "", 42, 0.3);
      hist2DDRDummySingle->GetXaxis()->SetRangeUser(doubleRatioXpp[0], doubleRatioXpp[j+1]);
      hist2DDRDummySingle->DrawCopy();
      // if(j==1)
      //   legendDRSingle = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*1,0.5,0.953, textSizeSinglePad, 2, "", 42, 0.3);
      legendDRSingle->Clear();
      DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[j+1], 1., 1., 1.2, kGray+2, 7);
      for (Int_t i = 10; i > -1; i--){
        // if(!(j==1&&i==2)){
          if (graphDRPi0FitSysErr[j][i]){
              DrawGammaSetMarkerTGraphAsym(graphDRPi0FitSysErr[j][i], markerStyleDet[0][i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
              graphDRPi0FitSysErr[j][i]->Draw("E2same");
              legendDRSingle->AddEntry(graphDRPi0FitSysErr[j][i],nameMeasGlobalLabel[i],"pf");
    
          }
          if (histoDRPi0FitStatErr[j][i]){
              DrawGammaSetMarker(histoDRPi0FitStatErr[j][i],  markerStyleDet[0][i], markerSizeDet[i], colorDet[i] , colorDet[i]);
              histoDRPi0FitStatErr[j][i]->Draw("p,same,e0,X0");
              if (!graphDRPi0FitSysErr[j][i])legendDRSingle->AddEntry(histoDRPi0FitStatErr[j][i],nameMeasGlobalLabel[i],"p");
          }
        // }
      }
    
      // if (j!=1&&histoDRPi0FitStatErr[j][2]) histoDRPi0FitStatErr[j][2]->Draw("p,same,e0,X0");
      if (histoDRPi0FitStatErr[j][2]) histoDRPi0FitStatErr[j][2]->Draw("p,same,e0,X0");
      if (histoDRPi0FitStatErr[j][0]) histoDRPi0FitStatErr[j][0]->Draw("p,same,e0,X0");
      if (histoDRPi0FitStatErr[j][4]) histoDRPi0FitStatErr[j][4]->Draw("p,same,e0,X0");
      legendDRSingle->Draw();
      labelDRSingle[j]->Draw();
      labelALICEDRSingle->Draw();
      hist2DDRDummySingle->Draw("same,axis");
    
      canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_%d.pdf", outputDir.Data(),j));
    }
    
    
    TLegend* legendDRComb = GetAndSetLegend2(0.12,0.95-textSizeSinglePad*1,0.5,0.95, textSizeSinglePad, 1, "", 42, 0.15);
    for (Int_t j = 0; j < 2; j++){
      hist2DDRDummySingle->DrawCopy();
      legendDRComb->Clear();
      DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
      DrawGammaSetMarkerTGraphAsym(graphCombDRSys[j], defaultMarker, defaultMarkerSize, defaultColor , defaultColor,widthLinesBoxes, kTRUE);
      DrawGammaSetMarkerTGraphAsym(graphCombDRStatPlot[j], defaultMarker, defaultMarkerSize, defaultColor , defaultColor, widthLinesBoxes);
      legendDRComb->AddEntry(graphCombDRSys[j],"ALICE","pf");
      graphCombDRSys[j]->Draw("E2same");
      graphCombDRStatPlot[j]->Draw("Epsame");
      legendDRComb->Draw();
      labelDRSingle[j]->Draw();
      labelALICEDRSingle->Draw();
      hist2DDRDummySingle->Draw("same,axis");
    
      canvasDoubleRatio->Print(Form("%s/DR_CombMeasurements_%d.pdf", outputDir.Data(),j));
    }
    
    
    Color_t  colorNLOWerner                         = kAzure+2;
    Color_t  colorNLOWernerBand                     = kAzure-9;
    Color_t  colorNLOMcGill                         = kGreen+2;
    Color_t  colorNLOMcGillBand                     = kGreen-6;
    Style_t  styleMarkerNLOWerner                   = 24;
    Style_t  styleLineNLOWerner                     = 5;
    Style_t  styleLineMcGill                        = 7;
    Width_t  widthLineNLO                           = 2.;
    
    hist2DDRDummySingle->DrawCopy();

    TLegend* legendDRTheoryComb = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*1,0.5,0.96, textSizeSinglePad, 1, "", 42, 0.15);
    legendDRTheoryComb->SetTextAlign(11);
    TLegend* legendDRTheoryComb2 = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*4,0.5,0.96-textSizeSinglePad*1, textSizeSinglePad, 1, "NLO pQCD:", 42, 0.15);
    legendDRTheoryComb2->SetTextAlign(11);

    ProduceGraphAsymmWithoutXErrors(graphCombDRStatPlot[0]);
    DrawGammaSetMarkerTGraphAsym(graphCombDRSys[0],defaultMarker, defaultMarkerSize, defaultColor , defaultColor,widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphCombDRStatPlot[0], defaultMarker, defaultMarkerSize, defaultColor , defaultColor, widthLinesBoxes);
    legendDRTheoryComb->AddEntry(graphCombDRSys[0],"ALICE","pf");

    if (graphTheoryNLODRpp[0]) {
        DrawGammaSetMarkerTGraphAsym(graphTheoryNLODRpp[0], 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
        graphTheoryNLODRpp[0]->Draw("3,same");
        legendDRTheoryComb2->AddEntry(graphTheoryNLODRpp[0],"PDF: CT10, FF: GRV","f");
    }
    if (graphTheoryNLODRppCenter[0]){
        DrawGammaNLOTGraph( graphTheoryNLODRppCenter[0], 2, styleLineNLOWerner, colorNLOWerner);
        graphTheoryNLODRppCenter[0]->Draw("lc,same");
    }
    if (graphTheoryNLODRppPaquettCenter[0]){
        DrawGammaNLOTGraph( graphTheoryNLODRppPaquettCenter[0], 2, styleLineMcGill, colorNLOMcGill );
        graphTheoryNLODRppPaquettCenter[0]->RemovePoint(0);
        graphTheoryNLODRppPaquettCenter[0]->Draw("lc,same");
        legendDRTheoryComb2->AddEntry(graphTheoryNLODRppPaquettCenter[0],"PDF: CTEQ6.1M, FF: BFG2","l");
    }
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

    graphCombDRSys[0]->Draw("E2same");
    graphCombDRStatPlot[0]->Draw("Epsame");
    legendDRTheoryComb->Draw();
    legendDRTheoryComb2->Draw();
    DrawGammaLines(0.25, 0.31, 1.492, 1.492, 2, kAzure+2, 7);

    labelALICEDRSingle->Draw();

    labelDRSingle[0]->Draw();
    hist2DDRDummySingle->Draw("same,axis");

canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pp2760GeV.pdf", outputDir.Data()));

    hist2DDRDummySingle->DrawCopy();

    legendDRTheoryComb = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*1,0.5,0.96, textSizeSinglePad, 1, "", 42, 0.15);
    legendDRTheoryComb->SetTextAlign(11);
    legendDRTheoryComb2 = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*4,0.5,0.96-textSizeSinglePad*1, textSizeSinglePad, 1, "NLO pQCD:", 42, 0.15);
    legendDRTheoryComb2->SetTextAlign(11);
    
    ProduceGraphAsymmWithoutXErrors(graphCombDRStatPlot[1]);
    DrawGammaSetMarkerTGraphAsym(graphCombDRSys[1],defaultMarker, defaultMarkerSize, defaultColor , defaultColor,widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphCombDRStatPlot[1], defaultMarker, defaultMarkerSize, defaultColor , defaultColor, widthLinesBoxes);
    legendDRTheoryComb->AddEntry(graphCombDRSys[1],"ALICE","pf");

    if (graphTheoryNLODRpp[1]) {
        DrawGammaSetMarkerTGraphAsym(graphTheoryNLODRpp[1], 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
        graphTheoryNLODRpp[1]->Draw("3,same");
        legendDRTheoryComb2->AddEntry(graphTheoryNLODRpp[1],"PDF: CT10, FF: GRV","f");
    }
    if (graphTheoryNLODRppCenter[1]){
        DrawGammaNLOTGraph( graphTheoryNLODRppCenter[1], 2, styleLineNLOWerner, colorNLOWerner);
        graphTheoryNLODRppCenter[1]->Draw("lc,same");
    }
    if (graphTheoryNLODRppPaquettCenter[1]){
        DrawGammaNLOTGraph( graphTheoryNLODRppPaquettCenter[1], 2, styleLineMcGill, colorNLOMcGill );
        graphTheoryNLODRppPaquettCenter[1]->RemovePoint(0);
        graphTheoryNLODRppPaquettCenter[1]->Draw("lc,same");
        legendDRTheoryComb2->AddEntry(graphTheoryNLODRppPaquettCenter[1],"PDF: CTEQ6.1M, FF: BFG2","l");
    }
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

    graphCombDRSys[1]->Draw("E2same");
    graphCombDRStatPlot[1]->Draw("Epsame");
    legendDRTheoryComb->Draw();
    legendDRTheoryComb2->Draw();
    DrawGammaLines(0.25, 0.31, 1.492, 1.492, 2, kAzure+2, 7);

    labelALICEDRSingle->Draw();

    labelDRSingle[1]->Draw();
    hist2DDRDummySingle->Draw("same,axis");

canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pp8TeV.pdf", outputDir.Data()));
    
    
    
    // **********************************************************************************************************************
    // ******************************************* Plot Ratio of Comb to Fit ************************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 54;
    TCanvas* canvasRatioToCombFit       = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToCombFit, 0.08, 0.01, 0.01, 0.125);
    canvasRatioToCombFit->SetLogx();

    Double_t textsizeLabelsPPb      = 0;
    if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
        textsizeLabelsPPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
    } else {
        textsizeLabelsPPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
    }
    cout << textsizeLabelsPPb << endl;

    TH2F * histo2DGammaRatioToCombFit;
    histo2DGammaRatioToCombFit               = new TH2F("histo2DGammaRatioToCombFit","histo2DGammaRatioToCombFit",1000,0.23, 25.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DGammaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPPb, textsizeLabelsPPb,
                              0.85*textsizeLabelsPPb,textsizeLabelsPPb, 0.9, 0.65, 510, 505);
    histo2DGammaRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DGammaRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
    //  histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(0.59,1.42);
    histo2DGammaRatioToCombFit->Draw("copy");
    
    ProduceGraphAsymmWithoutXErrors(graphRatioGammaCombCombFitStat[0]);
    
    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitSys[0], defaultMarker, defaultMarkerSize, defaultColor , defaultColor, widthLinesBoxes, kTRUE);
    graphRatioGammaCombCombFitSys[0]->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitStat[0],defaultMarker, defaultMarkerSize, defaultColor , defaultColor);
    graphRatioGammaCombCombFitStat[0]->Draw("p,same,z");
    
    DrawGammaLines(0.23, 25. , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0.23, 25. , 1.05, 1.05,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.95, 0.95,0.1, kGray, 7);
    
    TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, collisionSystempp2760GeV.Data());
    SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergy->Draw();
    TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICE.Data());
    SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitALICE->Draw();
    TLatex *labelRatioToFitGamma      = new TLatex(0.12, 0.92, "#gamma_{inc}");
    SetStyleTLatex( labelRatioToFitGamma, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
    labelRatioToFitGamma->Draw();
    
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_0.pdf",outputDir.Data()));
    
    
    histo2DGammaRatioToCombFit->Draw("copy");
    
    ProduceGraphAsymmWithoutXErrors(graphRatioGammaCombCombFitStat[1]);
    
    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitSys[1], defaultMarker, defaultMarkerSize, defaultColor , defaultColor, widthLinesBoxes, kTRUE);
    graphRatioGammaCombCombFitSys[1]->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitStat[1],defaultMarker, defaultMarkerSize, defaultColor , defaultColor);
    graphRatioGammaCombCombFitStat[1]->Draw("p,same,z");
    
    DrawGammaLines(0.23, 25. , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0.23, 25. , 1.05, 1.05,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.95, 0.95,0.1, kGray, 7);
    
    labelRatioToFitEnergy   = new TLatex(0.95, 0.92, collisionSystempp8000GeV.Data());
    SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergy->Draw();
    labelRatioToFitALICE->Draw();
    labelRatioToFitGamma->Draw();
    
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_1.pdf",outputDir.Data()));
    
    
    
    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DGammaRatioToCombFit->Draw("copy");

    for (Int_t i = 10; i > -1 ; i--){
        if (graphRatioGammaIndCombFitSys[0][i]){
            DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitSys[0][i], markerStyleDet[0][i] ,markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
            graphRatioGammaIndCombFitSys[0][i]->Draw("E2same");
        }
        if (graphRatioGammaIndCombFitStat[0][i]){
            ProduceGraphAsymmWithoutXErrors(graphRatioGammaIndCombFitStat[0][i]);
            DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitStat[0][i], markerStyleDet[0][i] ,markerSizeDet[i], colorDet[i], colorDet[i]);
            graphRatioGammaIndCombFitStat[0][i]->Draw("p,same,z");
        }
    }
    graphRatioGammaIndCombFitStat[0][4]->Draw("p,same,z");

    DrawGammaLines(0.23, 25. , 1., 1.,0.5, kGray+2);
    DrawGammaLines(0.23, 25. , 1.05, 1.05,0.5, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.95, 0.95,0.5, kGray, 7);
    DrawGammaLines(0.23, 25. , 1.1, 1.1,0.5, kGray, 9);
    DrawGammaLines(0.23, 25. , 0.9, 0.9,0.5, kGray, 9);

    labelRatioToFitEnergy   = new TLatex(0.95, 0.92, collisionSystempp2760GeV.Data());
    SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergy->Draw();
    labelRatioToFitALICE->Draw();
    labelRatioToFitGamma->Draw();
    histo2DGammaRatioToCombFit->Draw("same,axis");


    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    Double_t rowsLegendGammaRatio[5]          = {0.31, 0.26, 0.21, 0.16, 0.16};
    Double_t rowsLegendGammaRatioAbs[5]       = {0.74, 0.69, 0.64, 0.59 };
    Double_t columnsLegendGammaRatio[7]       = {0.14, 0.33, 0.43, 0.48, 0.7, 0.8};
    Double_t columnsLegendGammaRatioAbs[7]    = {0.215, 0.99, 1.62, 2, 6.6, 9.6};
    Double_t lengthBox                          = 0.22;
    Double_t heightBox                          = 0.03/2;
    //****************** first Column **************************************************
    TLatex *textPCMRatioGamma                 = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[1],nameMeasGlobalLabelGamma[0]);
    SetStyleTLatex( textPCMRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textPCMRatioGamma->Draw();
    TLatex *textEMCALRatioGamma               = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[2],nameMeasGlobalLabelGamma[2]);
    SetStyleTLatex( textEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textEMCALRatioGamma->Draw();
    // TLatex *textPCMEMCALRatioGamma            = new TLatex(columnsLegendGammaRatio[3],rowsLegendGammaRatio[1],nameMeasGlobalLabelGamma[4]);
    // SetStyleTLatex( textPCMEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    // textPCMEMCALRatioGamma->Draw();
    TLatex *textPCMEMCALRatioGamma            = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[3],nameMeasGlobalLabelGamma[4]);
    SetStyleTLatex( textPCMEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textPCMEMCALRatioGamma->Draw();

    //****************** second Column *************************************************
    TLatex *textStatRatioGamma                = new TLatex(columnsLegendGammaRatio[1],rowsLegendGammaRatio[0] ,"stat");
    SetStyleTLatex( textStatRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textStatRatioGamma->Draw();
    TLatex *textSysRatioGamma                 = new TLatex(columnsLegendGammaRatio[2] ,rowsLegendGammaRatio[0],"syst");
    SetStyleTLatex( textSysRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textSysRatioGamma->Draw();

    TMarker* markerPCMGammaRatio           = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[0][0],columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[1],1);
    markerPCMGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[0]);
    TMarker* markerEMCALGammaRatio         = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[0][2], columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[2],1);
    markerEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[1]);
    TMarker* markerPCMEMCALGammaRatio      = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[0][4], columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[3],1);
    markerPCMEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[2]);

    TBox* boxPCMGammaRatio                 = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[0][0], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[0]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[0]+ heightBox);
    boxPCMGammaRatio->Draw("l");
    TBox* boxEMCALGammaRatio               = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[0][2], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[1]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[1]+ heightBox);
    boxEMCALGammaRatio->Draw("l");
    TBox* boxPCMEMCALGammaRatio            = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[0][4], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[2]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[2]+ heightBox);
    boxPCMEMCALGammaRatio->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit_0.pdf",outputDir.Data()));
    
    
    canvasRatioToCombFit->cd();
    histo2DGammaRatioToCombFit->Draw("copy");

    for (Int_t i = 10; i > -1 ; i--){
      // if(i!=2){
        if (graphRatioGammaIndCombFitSys[1][i]){
            DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitSys[1][i], markerStyleDet[0][i] ,markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
            graphRatioGammaIndCombFitSys[1][i]->Draw("E2same");
        }
        if (graphRatioGammaIndCombFitStat[1][i]){
            ProduceGraphAsymmWithoutXErrors(graphRatioGammaIndCombFitStat[1][i]);
            DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitStat[1][i], markerStyleDet[0][i] ,markerSizeDet[i], colorDet[i], colorDet[i]);
            graphRatioGammaIndCombFitStat[1][i]->Draw("p,same,z");
        }
      // }
    }

    DrawGammaLines(0.23, 25. , 1., 1.,0.5, kGray+2);
    DrawGammaLines(0.23, 25. , 1.05, 1.05,0.5, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.95, 0.95,0.5, kGray, 7);
    DrawGammaLines(0.23, 25. , 1.1, 1.1,0.5, kGray, 9);
    DrawGammaLines(0.23, 25. , 0.9, 0.9,0.5, kGray, 9);

    labelRatioToFitEnergy   = new TLatex(0.95, 0.92, collisionSystempp8000GeV.Data());
    SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergy->Draw();
    labelRatioToFitALICE->Draw();
    labelRatioToFitGamma->Draw();
    histo2DGammaRatioToCombFit->Draw("same,axis");
    
    textPCMRatioGamma->Draw();
    textEMCALRatioGamma->Draw();
    textPCMEMCALRatioGamma->Draw();
    // textPCMEMCALRatioGamma            = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[2],nameMeasGlobalLabel[4]);
    // SetStyleTLatex( textPCMEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    // textPCMEMCALRatioGamma->Draw();

    //****************** second Column *************************************************
    textStatRatioGamma->Draw();
    textSysRatioGamma->Draw();

    
    markerPCMGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[0]);
    markerEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[1]);
    markerPCMEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[2]);
    // markerPCMEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[1]);

    boxPCMGammaRatio->Draw("l");
    boxEMCALGammaRatio->Draw("l");
    boxPCMEMCALGammaRatio->Draw("l");
    // boxPCMEMCALGammaRatio            = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[0][4], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[1]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[1]+ heightBox);
    // boxPCMEMCALGammaRatio->Draw("l");
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit_1.pdf",outputDir.Data()));
    
    
    
    
    
    
    TCanvas *canvasDirGamma = new TCanvas("canvasDirGamma","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasDirGamma, 0.175, 0.01, 0.01, 0.07);
    canvasDirGamma->SetLogy();
    canvasDirGamma->SetLogx();

    Int_t textSizeLabelsPixelDirGam = 48;
    Double_t textsizeLabelsDirGamma = 0;
    if (canvasDirGamma->XtoPixel(canvasDirGamma->GetX2()) < canvasDirGamma->YtoPixel(canvasDirGamma->GetY1())){
        textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGamma->XtoPixel(canvasDirGamma->GetX2()) ;
    } else {
        textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGamma->YtoPixel(canvasDirGamma->GetY1());
    }

    TH1D* dummyDirGamma = new TH1D("histo2DYieldGamma","histo2DYieldGamma",11000,doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs( dummyDirGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.8);
    dummyDirGamma->GetYaxis()->SetRangeUser(7e-11,3.5e0);
    dummyDirGamma->GetXaxis()->SetLabelOffset(-0.015);
    dummyDirGamma->GetXaxis()->SetTickLength(0.025);
    dummyDirGamma->GetYaxis()->SetTickLength(0.025);
    dummyDirGamma->GetXaxis()->SetMoreLogLabels();

    dummyDirGamma->GetXaxis()->SetRangeUser(doubleRatioXpp[0],doubleRatioXpp[1]);
    //     dummyDirGamma->GetXaxis()->SetRangeUser(0,16);
    dummyDirGamma->DrawCopy();
    
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErrPlot[2] = {NULL,NULL};
    TGraphAsymmErrors* graphCombIncGammaStatPlot[2] ={NULL,NULL};
    

    
    for(Int_t i=0;i<2;i++){
      if (graphCombDirGammaSpectrumStatErr[i]) graphCombDirGammaSpectrumStatErrPlot[i]       = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr[i]->Clone(Form("graphCombDirGammaSpectrumStatErrPlot_%d",i));
      if (graphCombDirGammaSpectrumStatErrPlot[i]) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErrPlot[i]);
      graphCombIncGammaStatPlot[i]  = (TGraphAsymmErrors*)graphCombIncGammaStat[i]->Clone("graphCombIncGammaStatPlot");
      ProduceGraphAsymmWithoutXErrors(graphCombIncGammaStatPlot[i]);
      DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys[i], defaultMarker, defaultMarkerSize, defaultColor , defaultColor,widthLinesBoxes, kTRUE);
      DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot[i], defaultMarker, defaultMarkerSize, defaultColor , defaultColor,widthLinesBoxes);
    }
    
    TLegend* legendYieldDirGamma        = GetAndSetLegend2(0.70, 0.84-(2*textSizeLabelsRel*0.85), 0.93, 0.84,textSizeLabelsPixel, 1, "", 43, 0.3);
    
    DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys[0], defaultMarker+4, defaultMarkerSize+0.2, defaultColor , defaultColor,widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot[0], defaultMarker+4, defaultMarkerSize+0.2, defaultColor , defaultColor, widthLinesBoxes);
    legendYieldDirGamma->AddEntry(graphCombIncGammaSys[0], "#gamma_{inc} ALICE","pf");
    graphCombIncGammaSys[0]->Draw("E2same");
    graphCombIncGammaStatPlot[0]->Draw("Epsame");
        
        if (graphCombDirGammaSpectrumSystErr[0]){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr[0], defaultMarker, defaultMarkerSize, defaultColor , defaultColor, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr[0]->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErrPlot[0]){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlot[0], defaultMarker, defaultMarkerSize, defaultColor , defaultColor);
            graphCombDirGammaSpectrumStatErrPlot[0]->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr[0]){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrAr[0] , 1, 3, defaultColor, defaultColor, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErrAr[0]->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr[0]);
            legendYieldDirGamma->AddEntry((TObject*)0, "#gamma_{dir} ALICE","");
        }
        TLatex *labelEnergyDGInvYieldPaperAll = new TLatex(0.94, 0.965-0.04*1, collisionSystempp2760GeV.Data());
        SetStyleTLatex( labelEnergyDGInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
        labelEnergyDGInvYieldPaperAll->Draw();
        TLatex *labelALICEDGInvYieldPaperAll  = new TLatex(0.94,0.965-0.04*2,textALICE.Data());
        SetStyleTLatex( labelALICEDGInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
        labelALICEDGInvYieldPaperAll->Draw();
        TLatex *labelALICEDGNormUnPaperAll    = new TLatex(0.94,0.965-(0.04*2+0.04*0.8),"Norm. unc. 2.5%");
        SetStyleTLatex( labelALICEDGNormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 31);
        labelALICEDGNormUnPaperAll->Draw();
        
        legendYieldDirGamma->Draw();


      TGraphAsymmErrors* dummyForLegend    = new TGraphAsymmErrors(1);
        dummyForLegend->SetPoint(0,4.2,0.02);
        dummyForLegend->SetPointError(0,0,0,0.009,0);
        DrawGammaSetMarkerTGraphAsym(dummyForLegend , 1, 3, defaultColor, defaultColor, 1.8, kTRUE);
        if (graphCombDirGammaSpectrumSumErrAr[0]){
            dummyForLegend->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend);
        }


    canvasDirGamma->Print(Form("%s/DirGammaSpectra_Unscaled_0.pdf",outputDir.Data()));
    
    dummyDirGamma->GetYaxis()->SetRangeUser(7e-10,3.5e1);
    dummyDirGamma->DrawCopy();
    legendYieldDirGamma->Clear();
    DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys[1], defaultMarker+4, defaultMarkerSize+0.2, defaultColor , defaultColor,widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot[1], defaultMarker+4, defaultMarkerSize+0.2, defaultColor , defaultColor, widthLinesBoxes);
    legendYieldDirGamma->AddEntry(graphCombIncGammaSys[1], "#gamma_{inc} ALICE","pf");
    graphCombIncGammaSys[1]->Draw("E2same");
    graphCombIncGammaStatPlot[1]->Draw("Epsame");
    
        if (graphCombDirGammaSpectrumSystErr[1]){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr[1], defaultMarker, defaultMarkerSize, defaultColor , defaultColor, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr[1]->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErrPlot[1]){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlot[1], defaultMarker, defaultMarkerSize, defaultColor , defaultColor);
            graphCombDirGammaSpectrumStatErrPlot[1]->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr[1]){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrAr[1] , 1, 3, defaultColor, defaultColor, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErrAr[1]->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr[1]);
            legendYieldDirGamma->AddEntry((TObject*)0, "#gamma_{dir} ALICE","");
        }
        labelEnergyDGInvYieldPaperAll = new TLatex(0.94, 0.965-0.04*1, collisionSystempp8000GeV.Data());
        SetStyleTLatex( labelEnergyDGInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll    = new TLatex(0.94,0.965-(0.04*2+0.04*0.8),"Norm. unc. 2.6%");
        SetStyleTLatex( labelALICEDGNormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 31);
        labelALICEDGNormUnPaperAll->Draw();
        
        legendYieldDirGamma->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectra_Unscaled_1.pdf",outputDir.Data()));


}