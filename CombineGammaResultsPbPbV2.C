#include <Riostream.h>
#include <fstream>
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
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
//#include "AliHEPDataParser.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"


void CombineGammaResultsPbPbV2( TString inputFileNamePCM    = "",
                                TString suffix              = "eps",
                                Bool_t combinePHOS          = kTRUE,
                                TString inputFileNamePHOS   = "ExternalInputPbPb/PHOS/PHOS_GammaDir_final_09022015.root",
                                Bool_t enablepValueCalcPCM  = kFALSE,
                                Bool_t enablepValueCalc     = kFALSE
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
    TString outputDir                                           = Form("%s/%s/CombineGammaMeasurementsPbPb",suffix.Data(),dateForOutput.Data());
    TString fileNameTheoryPbPb                                  = "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";
    TString fileNameExperimentPbPb                              = "ExternalInputPbPb/OtherExperiments/phenix_200.root";
    TString fileNameMesonPbPb = Form("../FinalMeson25July2017/pdf/2018_02_23/CombineMesonMeasurementsPbPb2760GeVX/CombinedResultsPaperPbPb2760GeV_2018_02_23.root",dateForOutput.Data(),dateForOutput.Data());
    TString fileNameMesonpPb = "ResultsRpPbpPb_2016_06_16_Preliminary.root";
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCMGammaPbPb.root", inputFileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOSGammaPbPb.root", inputFileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/TheoryPbPb.root", fileNameTheoryPbPb.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/OtherExperimentsPbPb.root", fileNameExperimentPbPb.Data(), outputDir.Data()));

    TString nameFinalResDat                                     = Form("%s/CombinedResultsGamma_FitResults_%s.dat",outputDir.Data(),dateForOutput.Data());
    fstream fileFinalResults;
    fileFinalResults.open(nameFinalResDat.Data(), ios::out);

    //*******************************************************************************************************************************************
    //******************************************************* set ranges for plotting ***********************************************************
    //*******************************************************************************************************************************************
    Double_t doubleRatio[2];
    Double_t indMeasRatio[2];
//     Double_t incRatio[2];
    Double_t incSpectraX[2];
    Double_t doubleRatioX[2];
    Double_t doubleRatioXpp[2];
    doubleRatio[0]      = 0.75;     doubleRatio[1]      = 1.95;
    indMeasRatio[0]     = 0.65;     indMeasRatio[1]     = 1.45;
//  incRatio[0]         = 0.0;      incRatio[1]         = 1.7;
    incSpectraX[0]      = 0.5;      incSpectraX[1]      = 20.;
    doubleRatioX[0]     = 0.7;      doubleRatioX[1]     = 16.;
    doubleRatioXpp[0]   = 0.23;     doubleRatioXpp[1]   = 23;

    Double_t nCollErr0020                                       = GetNCollErrFromName("0020");
    Double_t nCollErr0010                                       = GetNCollErrFromName("0010");
    Double_t nCollErr2040                                       = GetNCollErrFromName("2040");
    Double_t nCollErr2050                                       = GetNCollErrFromName("2050");

    Double_t nColl0020                                          = GetNCollFromName("0020");
    Double_t nColl0010                                          = GetNCollFromName("0010");
    Double_t nColl2040                                          = GetNCollFromName("2040");
    Double_t nColl2050                                          = GetNCollFromName("2050");

    Double_t normErr0020                                        = nCollErr0020/nColl0020;
    Double_t normErr0010                                        = nCollErr0010/nColl0010;
    Double_t normErr2040                                        = nCollErr2040/nColl2040;
    Double_t normErr2050                                        = nCollErr2050/nColl2050;
    Double_t normErrpPb = TMath::Sqrt(pow(0.031,2)+pow(0.036,2)+pow(0.036,2));//pPb normalization,TpPb and pp normalization errors

    Color_t colorCocktailPi0                                    = kRed+2;
    Color_t colorCocktailEta                                    = kBlue+1;
    Color_t colorCocktailEtaP                                   = kOrange+1;
    Color_t colorCocktailOmega                                  = kYellow+2;
    Color_t colorCocktailPhi                                    = kViolet;
    Color_t colorCocktailRho0                                   = kAzure-2;
    Color_t colorCocktailSigma0                                 = kGray+1;

    Color_t colorNLOcalc                                        = kBlack;//kBlue-7;
    Style_t fillStyleNLO                                        = 1001;
    Color_t colorEPS09calc                                      = kGray;
    Color_t colorEPS09calc2                                     = kGray+1;
    Style_t fillStyleEPS09                                      = 3008;
    Color_t colorCT10calc                                       = kAzure-4;
    Style_t fillStyleCT10                                       = 3001;
    Color_t colorNLOMcGill                                      = kOrange-3;
    Style_t styleNLOMcGill                                      = 7;
    Color_t colorPHSD                                           = kGray+2;
    Style_t stylePHSD                                           = 2;
    Color_t colorChatterjeeThermal                              = kTeal-6;
    Color_t colorChatterjeePrompt                               = kCyan+1;
    Color_t colorChatterjee                                     = kBlue+1;
    Style_t styleChatterjee                                     = 4;
    Color_t colorHees                                           = kBlack;
    Style_t styleHees                                           = 3;
    Color_t colorHe                                             = kBlue-7;
    Style_t styleHe                                             = 8;
    Style_t styleFit                                            = 1;
    Color_t colorRapp                                           = kMagenta-8;

    Color_t colorPHENIX                                         = kGray+2;
    Style_t markerStylePHENIX                                   = 24;
    Size_t markerSizePHENIX                                     = 2;

    Color_t colorPCM                                            = GetDefaultColorDiffDetectors("PCM", kFALSE, kFALSE, kTRUE);
    Color_t colorPHOS                                           = GetDefaultColorDiffDetectors("PHOS", kFALSE, kFALSE, kTRUE);
    Color_t colorPCMBox                                         = GetDefaultColorDiffDetectors("PCM", kFALSE, kTRUE, kTRUE);
    Color_t colorPHOSBox                                        = GetDefaultColorDiffDetectors("PHOS", kFALSE, kTRUE, kTRUE);
    Style_t markerStylePCM                                      = GetDefaultMarkerStyleDiffDetectors("PCM", kFALSE);
    Style_t markerStylePHOS                                     = GetDefaultMarkerStyleDiffDetectors("PHOS", kFALSE);
    Size_t markerSizePCM                                        = GetDefaultMarkerSizeDiffDetectors("PCM", kFALSE);
    Size_t markerSizePHOS                                       = GetDefaultMarkerSizeDiffDetectors("PHOS", kFALSE);

    Color_t colorComb0010                                       = GetColorDefaultColor("PbPb_2.76TeV", "", "0-10%");
    Color_t colorComb2040                                       = GetColorDefaultColor("PbPb_2.76TeV", "", "20-40%");
    Color_t colorComb2050                                       = GetColorDefaultColor("PbPb_2.76TeV", "", "20-50%");
    Style_t colorCombNotFit0010                                 = GetColorDefaultColor("PbPb_2.76TeV", "", "0-10%")-7;
    Style_t colorCombNotFit2040                                 = GetColorDefaultColor("PbPb_2.76TeV", "", "20-40%")-7;
    Style_t colorCombNotFit2050                                 = GetColorDefaultColor("PbPb_2.76TeV", "", "20-50%")-7;

    Color_t colorComb0010Box                                    = GetColorDefaultColor("PbPb_2.76TeV", "", "0-10%", kTRUE);
    Color_t colorComb2040Box                                    = GetColorDefaultColor("PbPb_2.76TeV", "", "20-40%", kTRUE);
    Color_t colorComb2050Box                                    = GetColorDefaultColor("PbPb_2.76TeV", "", "20-50%", kTRUE);
    Color_t colorCombpPb                                        = GetColorDefaultColor("pPb_5.023TeV", "", "");
    Color_t colorCombpPbBox                                     = GetColorDefaultColor("pPb_5.023TeV", "", "", kTRUE);
    Color_t colorCombpp7TeV                                     = GetColorDefaultColor("7TeV", "", "");
    Color_t colorCombpp7TeVBox                                  = GetColorDefaultColor("7TeV", "", "", kTRUE);
    Color_t colorCombpp2760GeV                                  = GetColorDefaultColor("2.76TeV", "", "");
    Color_t colorCombpp2760GeVBox                               = GetColorDefaultColor("2.76TeV", "", "", kTRUE);

    Style_t markerStyleComb0010                                 = GetDefaultMarkerStyle("PbPb_2.76TeV", "", "0-10%");
    Style_t markerStyleComb2040                                 = GetDefaultMarkerStyle("PbPb_2.76TeV", "", "20-40%");
    Style_t markerStyleComb2050                                 = GetDefaultMarkerStyle("PbPb_2.76TeV", "", "20-50%");

    Style_t markerStyleCombpPb                                  = GetDefaultMarkerStyle("pPb_5.023TeV", "", "");
    Style_t markerStyleCombpp7TeV                               = GetDefaultMarkerStyle("7TeV", "", "");
    Style_t markerStyleCombpp2760GeV                            = GetDefaultMarkerStyle("2.76TeV", "", "");
    Size_t markerSizeComb0010                                   = GetDefaultMarkerSize("PbPb_2.76TeV", "", "0-10%");
    Size_t markerSizeComb2040                                   = GetDefaultMarkerSize("PbPb_2.76TeV", "", "20-40%");
    Size_t markerSizeComb2050                                   = GetDefaultMarkerSize("PbPb_2.76TeV", "", "20-50%");
    Size_t markerSizeCombpPb                                    = GetDefaultMarkerSize("pPb_5.023TeV", "", "");
    Size_t markerSizeCombpp7TeV                                 = GetDefaultMarkerSize("7TeV", "", "");
    Size_t markerSizeCombpp2760GeV                              = GetDefaultMarkerSize("2.76TeV", "", "");

    Width_t widthLinesBoxes                                     = 1.4;
    Width_t widthCommonFit                                      = 2.4;

    TString collisionSystem                                     = "Pb#font[122]{-}Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString collisionSystemCent0010                             = "0#font[122]{-}10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString collisionSystemCent2040                             = "20#font[122]{-}40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString collisionSystemCent2050                             = "20#font[122]{-}50% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString collisionSystempPb                                  = "0#font[122]{-}100% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    TString collisionSystempp7TeV                               = "pp, #sqrt{#it{s}} = 7 TeV";
    TString collisionSystempp2760GeV                            = "pp, #sqrt{#it{s}} = 2.76 TeV";
    TString collisionSystemRHIC0010                             = "0#font[122]{-}10% Au-Au #sqrt{#it{s}_{_{NN}}} = 0.2 TeV";
    TString collisionSystemRHIC2040                             = "20#font[122]{-}40% Au-Au #sqrt{#it{s}_{_{NN}}} = 0.2 TeV";

    TString cent0010 = "0#font[122]{-}10%";
    TString cent2050 = "20#font[122]{-}50%";
    TString cent2040 = "20#font[122]{-}40%";

    TLatex *labelScalingDirGamma0010 = new TLatex(11.2,/*1.3*/3.5E-3/*3*/,"5 x 10^{2}");
    TLatex *labelScalingDirGamma2040 = new TLatex(11.2,/*6.0*/2.5E-4/*5*/,"5 x 10^{1}");
    TLatex *labelScalingDirGamma2050 = new TLatex(11.2,7.5E-7,"x 10^{0}");
    TLatex *labelDirGammaColl = new TLatex(0.25,0.94,Form("%s",collisionSystem.Data()));

    TLatex *labelThesisHighLeft = new TLatex(0.15,0.9,"This thesis");
    SetStyleTLatex( labelThesisHighLeft, 0.04,4);
    TLatex *labelThesisHighRight = new TLatex(0.82,0.83,"This thesis");
    SetStyleTLatex( labelThesisHighRight, 0.04,4);

    Int_t scaleFactor1 = 100;
    Int_t scaleFactor2 = 10;
    if(!combinePHOS){ scaleFactor1 = 500; scaleFactor2 = 50;}

    TFile* fileMesonPbPb                                      = new TFile( fileNameMesonPbPb.Data());
    TGraphAsymmErrors* graphRAAPi0CombPbPb2760GeVSysErr_0010  = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphRAAPi0CombPbPb2760GeVSysErr_0010");
    TGraphAsymmErrors* graphRAAPi0CombPbPb2760GeVStatErr_0010 = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphRAAPi0CombPbPb2760GeVStatErr_0010");
    TGraphAsymmErrors* graphRAAEtaCombPbPb2760GeVSysErr_0010  = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphRAAEtaCombPbPb2760GeVSysErr_0010");
    TGraphAsymmErrors* graphRAAEtaCombPbPb2760GeVStatErr_0010 = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphRAAEtaCombPbPb2760GeVStatErr_0010");
    TGraphAsymmErrors* graphRAAPi0CombPbPb2760GeVSysErr_2040  = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphPCMPi0RAASysPbPb2760GeV_2040");
    TGraphAsymmErrors* graphRAAPi0CombPbPb2760GeVStatErr_2040 = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphPCMPi0RAAStatPbPb2760GeV_2040");
    TGraphAsymmErrors* graphRAAEtaCombPbPb2760GeVSysErr_2040  = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphPCMEtaRAASysPbPb2760GeV_2040");
    TGraphAsymmErrors* graphRAAEtaCombPbPb2760GeVStatErr_2040 = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphPCMEtaRAAStatPbPb2760GeV_2040");
    TGraphAsymmErrors* graphRAAPi0CombPbPb2760GeVSysErr_2050  = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphRAAPi0CombPbPb2760GeVSysErr_2050");
    TGraphAsymmErrors* graphRAAPi0CombPbPb2760GeVStatErr_2050 = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphRAAPi0CombPbPb2760GeVStatErr_2050");
    TGraphAsymmErrors* graphRAAEtaCombPbPb2760GeVSysErr_2050  = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphRAAEtaCombPbPb2760GeVSysErr_2050");
    TGraphAsymmErrors* graphRAAEtaCombPbPb2760GeVStatErr_2050 = (TGraphAsymmErrors*)fileMesonPbPb->Get("graphRAAEtaCombPbPb2760GeVStatErr_2050");

    TFile* fileMesonpPb                             = new TFile( fileNameMesonpPb.Data());
    TGraphAsymmErrors*	CombinedPi0RpPbSystErr      = (TGraphAsymmErrors*)fileMesonpPb->Get("CombinedPi0RpPbSystErr");
    CombinedPi0RpPbSystErr->RemovePoint(CombinedPi0RpPbSystErr->GetN()-1);
    TGraphAsymmErrors*	CombinedPi0RpPbStatErr      = (TGraphAsymmErrors*)fileMesonpPb->Get("CombinedPi0RpPbStatErr");
    CombinedPi0RpPbStatErr->RemovePoint(CombinedPi0RpPbStatErr->GetN()-1);

    //*******************************************************************************************************************************************
    //***************************************************** Load theory predictions *************************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMGamma                             = new TFile( inputFileNamePCM.Data());
    TDirectory* directoryTheory                     = (TDirectory*)filePCMGamma->Get("Theory");
//         TGraphAsymmErrors* graphTheoryNLODR0010                 = (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatio_0-10%");
//         TGraphAsymmErrors* graphTheoryNLODR2040                 = (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatio_20-40%");
//         TGraphAsymmErrors* graphTheoryNLODR2050                 = (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatio_20-50%");
        TGraphAsymmErrors* graphTheoryNLODR0010                 = (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatiowithErr_0-10%");
        while(graphTheoryNLODR0010->GetX()[graphTheoryNLODR0010->GetN()-1]>14.)graphTheoryNLODR0010->RemovePoint(graphTheoryNLODR0010->GetN()-1);
        TGraphAsymmErrors* graphTheoryNLODR2040                 = (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatiowithErr_20-40%");
        while(graphTheoryNLODR2040->GetX()[graphTheoryNLODR2040->GetN()-1]>14.)graphTheoryNLODR2040->RemovePoint(graphTheoryNLODR2040->GetN()-1);
        TGraphAsymmErrors* graphTheoryNLODR2050                 = (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatiowithErr_20-50%");
        while(graphTheoryNLODR2050->GetX()[graphTheoryNLODR2050->GetN()-1]>14.)graphTheoryNLODR2050->RemovePoint(graphTheoryNLODR2050->GetN()-1);
        TGraphAsymmErrors* graphTheoryNLO0010                   = (TGraphAsymmErrors*) directoryTheory->Get("NLO_0-10%");
        while(graphTheoryNLO0010->GetX()[graphTheoryNLO0010->GetN()-1]>14.)graphTheoryNLO0010->RemovePoint(graphTheoryNLO0010->GetN()-1);
        TGraphAsymmErrors* graphTheoryNLO2040                   = (TGraphAsymmErrors*) directoryTheory->Get("NLO_20-40%");
        while(graphTheoryNLO2040->GetX()[graphTheoryNLO2040->GetN()-1]>14.)graphTheoryNLO2040->RemovePoint(graphTheoryNLO2040->GetN()-1);
        TGraphAsymmErrors* graphTheoryNLO2050                   = (TGraphAsymmErrors*) directoryTheory->Get("NLO_20-50%");
        while(graphTheoryNLO2050->GetX()[graphTheoryNLO2050->GetN()-1]>14.)graphTheoryNLO2050->RemovePoint(graphTheoryNLO2050->GetN()-1);
        TGraphAsymmErrors* graphTheoryEPS090010                 = (TGraphAsymmErrors*) directoryTheory->Get("EPS09_0-10%");
        TGraphAsymmErrors* graphTheoryRelErrEPS090010           = CalculateRelErrAsymmGraphAround1(graphTheoryEPS090010, "graphTheoryRelErrEPS090010");
        TGraphAsymmErrors* graphTheoryEPS092040                 = (TGraphAsymmErrors*) directoryTheory->Get("EPS09_20-40%");
        TGraphAsymmErrors* graphTheoryRelErrEPS092040           = CalculateRelErrAsymmGraphAround1(graphTheoryEPS092040, "graphTheoryRelErrEPS092040");
        TGraphAsymmErrors* graphTheoryEPS092050                 = (TGraphAsymmErrors*) directoryTheory->Get("EPS09_20-50%");
        TGraphAsymmErrors* graphTheoryRelErrEPS092050           = CalculateRelErrAsymmGraphAround1(graphTheoryEPS092050, "graphTheoryRelErrEPS092050");
        TGraphAsymmErrors* graphTheoryEPS09DR0010               = (TGraphAsymmErrors*) directoryTheory->Get("EPS09DoubleRatio_0-10%");
        TGraphAsymmErrors* graphTheoryEPS09DR2040               = (TGraphAsymmErrors*) directoryTheory->Get("EPS09DoubleRatio_20-40%");
        TGraphAsymmErrors* graphTheoryEPS09DR2050               = (TGraphAsymmErrors*) directoryTheory->Get("EPS09DoubleRatio_20-50%");
        TGraphAsymmErrors* graphTheoryCT100010                  = (TGraphAsymmErrors*) directoryTheory->Get("CT10BF_0-10%");
        TGraphAsymmErrors* graphTheoryCT102040                  = (TGraphAsymmErrors*) directoryTheory->Get("CT10BF_20-40%");
        TGraphAsymmErrors* graphTheoryCT102050                  = (TGraphAsymmErrors*) directoryTheory->Get("CT10BF_20-50%");
        TGraphAsymmErrors* graphTheoryCT10DR0010                = (TGraphAsymmErrors*) directoryTheory->Get("CT10BFDoubleRatio_0-10%");
        TGraphAsymmErrors* graphTheoryCT10DR2040                = (TGraphAsymmErrors*) directoryTheory->Get("CT10BFDoubleRatio_20-40%");
        TGraphAsymmErrors* graphTheoryCT10DR2050                = (TGraphAsymmErrors*) directoryTheory->Get("CT10BFDoubleRatio_20-50%");
        TGraphErrors* graphTheoryDRMcGillfromTF10010                     = (TGraphErrors*) directoryTheory->Get("McGillfromTF1DoubleRatio_0-10%");
        TGraphErrors* graphTheoryDRMcGillfromTF12040                     = (TGraphErrors*) directoryTheory->Get("McGillfromTF1DoubleRatio_20-40%");
        TGraphErrors* graphTheoryDRMcGillfromTF12050                     = (TGraphErrors*) directoryTheory->Get("McGillfromTF1DoubleRatio_20-50%");
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRMcGillfromTF10010, 3.0, styleNLOMcGill, colorNLOMcGill, 3015, colorNLOMcGill, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRMcGillfromTF12040, 3.0, styleNLOMcGill, colorNLOMcGill, 3015, colorNLOMcGill, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRMcGillfromTF12050, 3.0, styleNLOMcGill, colorNLOMcGill, 3015, colorNLOMcGill, 0);
        TGraphErrors* graphTheoryDRMcGill0010                     = (TGraphErrors*) directoryTheory->Get("McGillDoubleRatio_0-10%");
        TGraphErrors* graphTheoryDRMcGill2040                     = (TGraphErrors*) directoryTheory->Get("McGillDoubleRatio_20-40%");
        TGraphErrors* graphTheoryDRMcGill2050                     = (TGraphErrors*) directoryTheory->Get("McGillDoubleRatio_20-50%");
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRMcGill0010, 3.0, styleNLOMcGill, colorNLOMcGill, 3015, colorNLOMcGill, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRMcGill2040, 3.0, styleNLOMcGill, colorNLOMcGill, 3015, colorNLOMcGill, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRMcGill2050, 3.0, styleNLOMcGill, colorNLOMcGill, 3015, colorNLOMcGill, 0);
        TGraphErrors* graphTheoryDRPHSD0010                     = (TGraphErrors*) directoryTheory->Get("PHSDDoubleRatio_0-10%");
        TGraphErrors* graphTheoryDRPHSD2040                     = (TGraphErrors*) directoryTheory->Get("PHSDDoubleRatio_20-40%");
        TGraphErrors* graphTheoryDRPHSD2050                     = (TGraphErrors*) directoryTheory->Get("PHSDDoubleRatio_20-50%");
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRPHSD0010, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRPHSD2040, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRPHSD2050, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);

        TGraphErrors* graphTheoryDRChatterjee0010               = (TGraphErrors*) directoryTheory->Get("ChatterjeeDoubleRatio_0-10%");
        TGraphErrors* graphTheoryDRChatterjeeSummed0010         = (TGraphErrors*) directoryTheory->Get("ChatterjeeSummedDoubleRatio_0-10%");
        TGraphErrors* graphTheoryDRChatterjee2040               = (TGraphErrors*) directoryTheory->Get("ChatterjeeDoubleRatio_20-40%");
        TGraphErrors* graphTheoryDRChatterjee2050               = (TGraphErrors*) directoryTheory->Get("ChatterjeeDoubleRatio_20-50%");
        TGraphErrors* graphTheoryDRChatterjeeSummed2050         = (TGraphErrors*) directoryTheory->Get("ChatterjeeSummedDoubleRatio_20-50%");
        TGraphErrors* graphTheoryDRChatterjeePrompt0010         = (TGraphErrors*) directoryTheory->Get("ChatterjeePromptDoubleRatio_0-10%");
        TGraphErrors* graphTheoryDRChatterjeePrompt2050         = (TGraphErrors*) directoryTheory->Get("ChatterjeePromptDoubleRatio_20-50%");
        TGraphErrors* graphTheoryDRChatterjeeThermal0010        = (TGraphErrors*) directoryTheory->Get("ChatterjeeThermalDoubleRatio_0-10%");
        TGraphErrors* graphTheoryDRChatterjeeThermal2050        = (TGraphErrors*) directoryTheory->Get("ChatterjeeThermalDoubleRatio_20-50%");
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRChatterjee0010, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRChatterjee2040, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRChatterjee2050, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRChatterjeeSummed0010, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRChatterjeeSummed2050, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRChatterjeePrompt0010, 3.0, styleChatterjee, colorChatterjeePrompt, 3015, colorChatterjeePrompt, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRChatterjeePrompt2050, 3.0, styleChatterjee, colorChatterjeePrompt, 3015, colorChatterjeePrompt, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRChatterjeeThermal0010, 3.0, styleChatterjee, colorChatterjeeThermal, 3015, colorChatterjeeThermal, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRChatterjeeThermal2050, 3.0, styleChatterjee, colorChatterjeeThermal, 3015, colorChatterjeeThermal, 0);
        TGraphErrors* graphTheoryDRRapp0010                     = (TGraphErrors*) directoryTheory->Get("RappDoubleRatio_0-10%");
        TGraphErrors* graphTheoryDRRapp2040                     = (TGraphErrors*) directoryTheory->Get("RappDoubleRatio_20-40%");
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRRapp0010, 3.0, styleHe, colorHe, 3015, colorHe, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryDRRapp2040, 3.0, styleHe, colorHe, 3015, colorHe, 0);

        TGraphAsymmErrors* graphTheoryNLO0010Plot                                                = ScaleGraph(graphTheoryNLO0010,scaleFactor1);
        TGraphAsymmErrors* graphTheoryEPS090010Plot                                              = ScaleGraph(graphTheoryEPS090010,scaleFactor1);
        TGraphAsymmErrors* graphTheoryCT100010Plot                                               = ScaleGraph(graphTheoryCT100010,scaleFactor1);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLO0010Plot, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryEPS090010Plot, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryCT100010Plot, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);

        TGraphAsymmErrors* graphTheoryNLO2040Plot                                                = ScaleGraph(graphTheoryNLO2040,scaleFactor2);
        TGraphAsymmErrors* graphTheoryEPS092040Plot                                              = ScaleGraph(graphTheoryEPS092040,scaleFactor2);
        TGraphAsymmErrors* graphTheoryCT102040Plot                                               = ScaleGraph(graphTheoryCT102040,scaleFactor2);
        SetStyleGammaNLOTGraphWithBand( graphTheoryEPS092040Plot, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryCT102040Plot, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLO2040Plot, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);

        TGraphAsymmErrors* graphTheoryNLO2050Plot                                                = ScaleGraph(graphTheoryNLO2050,1);
        TGraphAsymmErrors* graphTheoryEPS092050Plot                                              = ScaleGraph(graphTheoryEPS092050,1);
        TGraphAsymmErrors* graphTheoryCT102050Plot                                               = ScaleGraph(graphTheoryCT102050,1);
        SetStyleGammaNLOTGraphWithBand( graphTheoryEPS092050Plot, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryCT102050Plot, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLO2050Plot, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);

        TF1* fitTheoryGammaNLO0010  = FitObject("l","fitTheoryGammaNLO0010","Photon",graphTheoryNLO0010,1.6,12.5,NULL,"QNRMEX0+");
        cout << WriteParameterToFile(fitTheoryGammaNLO0010)<< endl;
        fitTheoryGammaNLO0010->SetRange(0.9,14);
        TF1* fitTheoryGammaNLO2050  = FitObject("l","fitTheoryGammaNLO2050","Photon",graphTheoryNLO2050,1.6,12.5,NULL,"QNRMEX0+");
        cout << WriteParameterToFile(fitTheoryGammaNLO2050)<< endl;
        fitTheoryGammaNLO2050->SetRange(0.9,14);

        TF1* fitTheoryGammaEPSO90010  = FitObject("l","fitTheoryGammaEPSO90010","Photon",graphTheoryEPS090010,1.6,12.5,NULL,"QNRMEX0+");
        cout << WriteParameterToFile(fitTheoryGammaEPSO90010)<< endl;
        fitTheoryGammaEPSO90010->SetRange(0.9,14);
        TF1* fitTheoryGammaEPSO92040  = FitObject("l","fitTheoryGammaEPSO90010","Photon",graphTheoryEPS092040,1.6,12.5,NULL,"QNRMEX0+");
        cout << WriteParameterToFile(fitTheoryGammaEPSO92040)<< endl;
        fitTheoryGammaEPSO92040->SetRange(0.9,14);
        TF1* fitTheoryGammaEPSO92050  = FitObject("l","fitTheoryGammaEPSO90010","Photon",graphTheoryEPS092050,1.6,12.5,NULL,"QNRMEX0+");
        cout << WriteParameterToFile(fitTheoryGammaEPSO92050)<< endl;
        fitTheoryGammaEPSO92050->SetRange(0.9,14);

    TFile* fileTheoryPbPb                           = new TFile( fileNameTheoryPbPb.Data());
    TDirectory* directoryTheoryGamma                = (TDirectory*)fileTheoryPbPb->Get("DirectPhoton");

    //______________________________________________________McGill_________________________________________
        TGraphErrors* graphTheoryMcGill0010                     = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_McGill_0010");
        TGraphErrors* graphTheoryMcGill2050                     = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_McGill_2050");
        TGraphErrors* graphTheoryMcGill2040                     = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_McGill_2040");

        TGraphErrors* graphTheoryMcGillforFit0010 = (TGraphErrors*)graphTheoryMcGill0010->Clone("graphTheoryMcGillforFit0010");
        TGraphErrors* graphTheoryMcGillforFit2040 = (TGraphErrors*)graphTheoryMcGill2040->Clone("graphTheoryMcGillforFit2040");
        for(Int_t i=0; i<graphTheoryMcGillforFit0010->GetN(); i++){
            graphTheoryMcGillforFit0010->SetPointError(i,0.,0.05);
            graphTheoryMcGillforFit2040->SetPointError(i,0.,0.05);
        }
//         TF1* fitExpDirGammaTheoryMcGill0010  = FitObject("e","fitExpDirGammaTheoryMcGill0010","Photon",graphTheoryMcGillforFit0010,0.97,1.6,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryMcGill0010)<< endl;
//         TF1* fitExpDirGammaTheoryMcGill0010_22  = FitObject("e","fitExpDirGammaTheoryMcGill0010_22","Photon",graphTheoryMcGillforFit0010,0.97,2.3,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryMcGill0010_22)<< endl;
//         TF1* fitExpDirGammaTheoryMcGill2040  = FitObject("e","fitExpDirGammaTheoryMcGill2040","Photon",graphTheoryMcGillforFit2040,0.97,1.6,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryMcGill2040)<< endl;
//         TF1* fitExpDirGammaTheoryMcGill2040_22  = FitObject("e","fitExpDirGammaTheoryMcGill2040_22","Photon",graphTheoryMcGillforFit2040,0.97,2.3,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryMcGill2040_22)<< endl;

        TGraphErrors* graphTheoryPromptMcGill2040               = (TGraphErrors*) directoryTheoryGamma->Get("graphPromptPhotonYield_McGill_2040");

        TF1* fitTheoryPromptMcGill0010  = new TF1 ("fitTheoryPromptMcGill0010","[0]/64.*1e-9*TMath::Exp(16.20-3.94*TMath::Log(x)-0.269*TMath::Log(x)**2)*1./(0.865779/TMath::Power(4,0.0694875))",0.1, 20);
        fitTheoryPromptMcGill0010->SetParameter(0,nColl0010);
        TF1* fitTheoryPromptMcGill2040  = new TF1 ("fitTheoryPromptMcGill2040","[0]/64.*1e-9*TMath::Exp(16.20-3.94*TMath::Log(x)-0.269*TMath::Log(x)**2)*1./(0.865779/TMath::Power(4,0.0694875))",0.1, 20);
        fitTheoryPromptMcGill2040->SetParameter(0,nColl2040);
        TF1* fitTheoryPromptMcGill2050  = new TF1 ("fitTheoryPromptMcGill2050","[0]/64.*1e-9*TMath::Exp(16.20-3.94*TMath::Log(x)-0.269*TMath::Log(x)**2)*1./(0.865779/TMath::Power(4,0.0694875))",0.1, 20);
        fitTheoryPromptMcGill2050->SetParameter(0,nColl2050);

        TF1* fitTheoryPromptforNLOMcGill0010 = (TF1*) fitTheoryPromptMcGill0010->Clone();
        fitTheoryPromptforNLOMcGill0010->SetParameter(0,fitTheoryPromptforNLOMcGill0010->GetParameter(0)*scaleFactor1);
        fitTheoryPromptforNLOMcGill0010->SetLineColor(colorNLOMcGill);
        fitTheoryPromptforNLOMcGill0010->SetLineStyle(styleNLOMcGill);
        fitTheoryPromptforNLOMcGill0010->SetRange(1.,14.);
        TF1* fitTheoryPromptforNLOMcGill2040 = (TF1*) fitTheoryPromptMcGill2040->Clone();
        fitTheoryPromptforNLOMcGill2040->SetParameter(0,fitTheoryPromptforNLOMcGill2040->GetParameter(0)*scaleFactor2);
        fitTheoryPromptforNLOMcGill2040->SetLineColor(colorNLOMcGill);
        fitTheoryPromptforNLOMcGill2040->SetLineStyle(styleNLOMcGill);
        fitTheoryPromptforNLOMcGill2040->SetRange(1.,14.);
        TF1* fitTheoryPromptforNLOMcGill2050 = (TF1*) fitTheoryPromptMcGill2050->Clone();
        fitTheoryPromptforNLOMcGill2050->SetLineColor(colorNLOMcGill);
        fitTheoryPromptforNLOMcGill2050->SetLineStyle(styleNLOMcGill);
        fitTheoryPromptforNLOMcGill2050->SetRange(1.,14.);

    //______________________________________________________PHSD_________________________________________
        TGraphErrors* graphTheoryPHSD0010                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_PHSD_0010");
        graphTheoryPHSD0010->RemovePoint(0);
        TGraphErrors* graphTheoryPHSD2050                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_PHSD_2050");
        graphTheoryPHSD2050->RemovePoint(0);
        TGraphErrors* graphTheoryPHSDplusPrompt0010             = (TGraphErrors*) graphTheoryPHSD0010->Clone("graphDirectPhotonYieldplusPrompt_PHSD_0010");
        Double_t* xValue0010 = graphTheoryPHSDplusPrompt0010->GetX();
        Double_t* yValue0010 = graphTheoryPHSDplusPrompt0010->GetY();
        for(Int_t i = 0; i < graphTheoryPHSDplusPrompt0010->GetN(); i++){
                Double_t newyValue0010 = yValue0010[i] + fitTheoryPromptMcGill0010->Eval(xValue0010[i]);
                graphTheoryPHSDplusPrompt0010->SetPoint(i,xValue0010[i],newyValue0010);
        }

        TGraphErrors* graphTheoryPHSDplusPrompt2050             = (TGraphErrors*) graphTheoryPHSD2050->Clone("graphDirectPhotonYieldplusPrompt_PHSD_2050");
        Double_t* xValue2050 = graphTheoryPHSDplusPrompt2050->GetX();
        Double_t* yValue2050 = graphTheoryPHSDplusPrompt2050->GetY();
        for(Int_t i = 0; i < graphTheoryPHSDplusPrompt2050->GetN(); i++){
                Double_t newyValue2050 = yValue2050[i] + fitTheoryPromptMcGill2050->Eval(xValue2050[i]);
                graphTheoryPHSDplusPrompt2050->SetPoint(i,xValue2050[i],newyValue2050);
        }
        TGraphErrors* graphTheoryPHSD2040                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_PHSD_2040");


        TGraphErrors* graphTheoryPHSDforFit0010 = (TGraphErrors*)graphTheoryPHSDplusPrompt0010->Clone("graphTheoryPHSDforFit0010");
        TGraphErrors* graphTheoryPHSDforFit2040 = (TGraphErrors*)graphTheoryPHSD2040->Clone("graphTheoryPHSDforFit2040");
        for(Int_t i=0; i<graphTheoryPHSDforFit0010->GetN(); i++){
            graphTheoryPHSDforFit0010->SetPointError(i,0.,0.05);
            graphTheoryPHSDforFit2040->SetPointError(i,0.,0.05);
        }
//         TF1* fitExpDirGammaTheoryPHSDplusPrompt0010  = FitObject("e","fitExpDirGammaTheoryPHSDplusPrompt0010","Photon",graphTheoryPHSDforFit0010,0.97,1.6,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryPHSDplusPrompt0010)<< endl;
//         TF1* fitExpDirGammaTheoryPHSDplusPrompt0010_22  = FitObject("e","fitExpDirGammaTheoryPHSDplusPrompt0010_22","Photon",graphTheoryPHSDforFit0010,0.97,2.3,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryPHSDplusPrompt0010_22)<< endl;
//         TF1* fitExpDirGammaTheoryPHSDplusPrompt2040  = FitObject("e","fitExpDirGammaTheoryPHSDplusPrompt2040","Photon",graphTheoryPHSDforFit2040,0.97,1.6,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryPHSDplusPrompt2040)<< endl;
//         TF1* fitExpDirGammaTheoryPHSDplusPrompt2040_22  = FitObject("e","fitExpDirGammaTheoryPHSDplusPrompt2040_22","Photon",graphTheoryPHSDforFit2040,0.97,2.3,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryPHSDplusPrompt2040_22)<< endl;


    //______________________________________________________Chatterjee_________________________________________
        TGraphErrors* graphTheoryChatterjee0010                 = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_0010");
        TGraphErrors* graphTheoryChatterjeeSummed0010           = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonSummedYield_Chatterjee_0010");
        while(graphTheoryChatterjeeSummed0010->GetX()[graphTheoryChatterjeeSummed0010->GetN()-1]>14.) graphTheoryChatterjeeSummed0010->RemovePoint(graphTheoryChatterjeeSummed0010->GetN()-1);

        TGraphErrors* graphTheoryChatterjeePrompt0010           = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonPromptYield_Chatterjee_0010");
        TGraphErrors* graphTheoryChatterjeeThermal0010          = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonThermalYield_Chatterjee_0010");
        TGraphErrors* graphTheoryChatterjee2050                 = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_2050");
        TGraphErrors* graphTheoryChatterjeeSummed2050           = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonSummedYield_Chatterjee_2050");
        while(graphTheoryChatterjeeSummed2050->GetX()[graphTheoryChatterjeeSummed2050->GetN()-1]>14.) graphTheoryChatterjeeSummed2050->RemovePoint(graphTheoryChatterjeeSummed2050->GetN()-1);
        TGraphErrors* graphTheoryChatterjeePrompt2050           = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonPromptYield_Chatterjee_2050");
        TGraphErrors* graphTheoryChatterjeeThermal2050          = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonThermalYield_Chatterjee_2050");
        TGraphErrors* graphTheoryChatterjee2040                 = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_2040_2");

        TGraphErrors* graphTheoryChatterjeeforFit0010 = (TGraphErrors*)graphTheoryChatterjeeSummed0010->Clone("graphTheoryChatterjeeforFit0010");
        TGraphErrors* graphTheoryChatterjeeforFit2040 = (TGraphErrors*)graphTheoryChatterjee2040->Clone("graphTheoryChatterjeeforFit2040");
        for(Int_t i=0; i<graphTheoryChatterjeeforFit0010->GetN(); i++){
            graphTheoryChatterjeeforFit0010->SetPointError(i,0.,0.05);
            graphTheoryChatterjeeforFit2040->SetPointError(i,0.,0.05);
        }
//         TF1* fitExpDirGammaTheoryChatterjee0010  = FitObject("e","fitExpDirGammaTheoryChatterjee0010","Photon",graphTheoryChatterjeeforFit0010,1.,1.6,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryChatterjee0010)<< endl;
//         TF1* fitExpDirGammaTheoryChatterjee0010_22  = FitObject("e","fitExpDirGammaTheoryChatterjee0010_22","Photon",graphTheoryChatterjeeforFit0010,1.,2.3,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryChatterjee0010_22)<< endl;
//         TF1* fitExpDirGammaTheoryChatterjee2040  = FitObject("e","fitExpDirGammaTheoryChatterjee2040","Photon",graphTheoryChatterjeeforFit2040,1.,1.6,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryChatterjee2040)<< endl;
//         TF1* fitExpDirGammaTheoryChatterjee2040_22  = FitObject("e","fitExpDirGammaTheoryChatterjee2040_22","Photon",graphTheoryChatterjeeforFit2040,1.,2.3,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryChatterjee2040_22)<< endl;


    //______________________________________________________new from Rapp_________________________________________
        TGraphErrors* graphTheoryRapp0010                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonRapp276GeV_0010");
        TGraphErrors* graphTheoryRapp2040                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonRapp276GeV_2040");

        TGraphErrors* graphTheoryRappforFit0010 = (TGraphErrors*)graphTheoryRapp0010->Clone("graphTheoryRappforFit0010");
        TGraphErrors* graphTheoryRappforFit2040 = (TGraphErrors*)graphTheoryRapp2040->Clone("graphTheoryRappforFit2040");
        for(Int_t i=0; i<graphTheoryRappforFit0010->GetN(); i++){
            graphTheoryRappforFit0010->SetPointError(i,0.,0.05);
            graphTheoryRappforFit2040->SetPointError(i,0.,0.05);
        }
//         TF1* fitExpDirGammaTheoryRapp0010  = FitObject("e","fitExpDirGammaTheoryRapp0010","Photon",graphTheoryRappforFit0010,0.97,1.6,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryRapp0010)<< endl;
//         TF1* fitExpDirGammaTheoryRapp0010_22  = FitObject("e","fitExpDirGammaTheoryRapp0010_22","Photon",graphTheoryRappforFit0010,0.97,2.3,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryRapp0010_22)<< endl;
//         TF1* fitExpDirGammaTheoryRapp2040  = FitObject("e","fitExpDirGammaTheoryRapp2040","Photon",graphTheoryRappforFit2040,0.97,1.6,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryRapp2040)<< endl;
//         TF1* fitExpDirGammaTheoryRapp2040_22  = FitObject("e","fitExpDirGammaTheoryRapp2040_22","Photon",graphTheoryRappforFit2040,0.97,2.3,NULL,"NRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitExpDirGammaTheoryRapp2040_22)<< endl;


    //______________________________________________________VanHees_________________________________________
        TGraphErrors* graphTheoryHees2040                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_VanHees_2040");
        TGraphErrors* graphTheoryHe2040                         = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_He_2040");

        TGraphErrors* graphTheoryMcGill0010Plot    = ScaleGraph(graphTheoryMcGill0010,scaleFactor1);
        TGraph* graphTheoryPHSD0010Plot            = ScaleGraph(graphTheoryPHSD0010,scaleFactor1);
        TGraph* graphTheoryPHSDplusPrompt0010Plot            = ScaleGraph(graphTheoryPHSDplusPrompt0010,scaleFactor1);
        TGraph* graphTheoryChatterjee0010Plot      = ScaleGraph(graphTheoryChatterjee0010,scaleFactor1);
        TGraph* graphTheoryChatterjeeSummed0010Plot      = ScaleGraph(graphTheoryChatterjeeSummed0010,scaleFactor1);
        TGraph* graphTheoryChatterjeePrompt0010Plot      = ScaleGraph(graphTheoryChatterjeePrompt0010,scaleFactor1);
        TGraph* graphTheoryChatterjeeThermal0010Plot      = ScaleGraph(graphTheoryChatterjeeThermal0010,scaleFactor1);
        TGraph* graphTheoryRapp0010Plot            = ScaleGraph(graphTheoryRapp0010,scaleFactor1);
        SetStyleGammaNLOTGraphWithBand( graphTheoryMcGill0010Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryPHSDplusPrompt0010Plot, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryChatterjee0010Plot, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryChatterjeeSummed0010Plot, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRapp0010Plot, 3.0, styleHe, colorHe, 3015, colorHe, 0);

        TGraphErrors* graphTheoryMcGill2040Plot    = ScaleGraph(graphTheoryMcGill2040,scaleFactor2);
        TGraphErrors* graphTheoryPromptMcGill2040Plot           = ScaleGraph(graphTheoryPromptMcGill2040,scaleFactor2);
        SetStyleGammaNLOTGraphWithBand( graphTheoryPromptMcGill2040Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
        TGraph* graphTheoryPHSD2040Plot            = ScaleGraph(graphTheoryPHSD2040,scaleFactor2);
        graphTheoryPHSD2040Plot->RemovePoint(0);
        TGraph* graphTheoryChatterjee2040Plot      = ScaleGraph(graphTheoryChatterjee2040,scaleFactor2);
        TGraph* graphTheoryRapp2040Plot            = ScaleGraph(graphTheoryRapp2040,scaleFactor2);
        TGraph* graphTheoryHees2040Plot            = ScaleGraph(graphTheoryHees2040,scaleFactor2);
        TGraph* graphTheoryHe2040Plot              = ScaleGraph(graphTheoryHe2040,scaleFactor2);
        SetStyleGammaNLOTGraphWithBand( graphTheoryMcGill2040Plot, 3.0, styleChatterjee, colorNLOMcGill+11, 3001, colorNLOMcGill+11, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryPHSD2040Plot, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryChatterjee2040Plot, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryHees2040Plot, 3.0, styleHees, colorHees, 3015, colorHees, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryHe2040Plot, 3.0, styleHe, colorHe, 3015, colorHe, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRapp2040Plot, 3.0, styleHe, colorHe, 3015, colorHe, 0);

        TGraphErrors* graphTheoryMcGill2050Plot    = ScaleGraph(graphTheoryMcGill2050,1);
        TGraph* graphTheoryPHSD2050Plot                      = ScaleGraph(graphTheoryPHSD2050,1);
        TGraph* graphTheoryPHSDplusPrompt2050Plot            = ScaleGraph(graphTheoryPHSDplusPrompt2050,1);
        TGraph* graphTheoryChatterjee2050Plot                = ScaleGraph(graphTheoryChatterjee2050,1);
        TGraph* graphTheoryChatterjeeSummed2050Plot          = ScaleGraph(graphTheoryChatterjeeSummed2050,1);
        TGraph* graphTheoryChatterjeePrompt2050Plot          = ScaleGraph(graphTheoryChatterjeePrompt2050,1);
        TGraph* graphTheoryChatterjeeThermal2050Plot         = ScaleGraph(graphTheoryChatterjeeThermal2050,1);
        SetStyleGammaNLOTGraphWithBand( graphTheoryMcGill2050Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryPHSDplusPrompt2050Plot, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryChatterjee2050Plot, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryChatterjeeSummed2050Plot, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);


    //*******************************************************************************************************************************************
    //***************************************** Loading previous combined and PCM measurement ***************************************************
    //*******************************************************************************************************************************************
    TString PubCombGammaFileName = "ExternalInputPbPb/CombDirGamma/Gamma_CombResults_PbPb_2.76TeV_20150729_Pub2015.root";
    TString PubPCMGammaFileName  = "ExternalInputPbPb/PCM/Gamma_PCMResults_PbPb_2.76TeV_20150729_Pub2015.root";
    TFile *filePubGammaPCM       = new TFile(PubPCMGammaFileName.Data());
    TFile *filePubGammaComb      = new TFile(PubCombGammaFileName.Data());
    if (filePubGammaPCM->IsZombie() || filePubGammaComb->IsZombie()) {
        cout << "published gamma file couldn't be read, aborting....";
        return;
    }
    TDirectoryFile* directoryCombGamma_0020                         = (TDirectoryFile*)filePubGammaComb->Get("Gamma_PbPb_2.76TeV_0-20%");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumStat_0020    = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DirGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumSyst_0020    = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DirGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumStat_0020   = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("IncGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumSyst_0020   = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("IncGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioStat_0020         = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DR_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioSyst_0020         = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DR_comb_SysErr");

    TDirectoryFile* directoryCombGamma_2040                         = (TDirectoryFile*)filePubGammaComb->Get("Gamma_PbPb_2.76TeV_20-40%");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumStat_2040    = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DirGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumSyst_2040    = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DirGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumStat_2040   = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("IncGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumSyst_2040   = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("IncGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioStat_2040         = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DR_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioSyst_2040         = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DR_comb_SysErr");

    DrawGammaSetMarkerTGraphAsym(graphPubCombDirGammaSpectrumStat_0020 , 20,2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombDirGammaSpectrumSyst_0020 , 20, 2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombInclGammaSpectrumStat_0020 , 20,2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombInclGammaSpectrumSyst_0020 , 20, 2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombDoubleRatioStat_0020 , 20,2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombDoubleRatioSyst_0020 , 20, 2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombDirGammaSpectrumStat_2040 , 20,2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombDirGammaSpectrumSyst_2040 , 20, 2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombInclGammaSpectrumStat_2040 , 20,2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombInclGammaSpectrumSyst_2040 , 20, 2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombDoubleRatioStat_2040 , 20,2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPubCombDoubleRatioSyst_2040 , 20, 2, kGray+1, kGray+1, 1, kTRUE);



    TDirectoryFile* directoryPCMGamma_0010                         = (TDirectoryFile*)filePubGammaPCM->Get("Gamma_PbPb_2.76TeV_0-10%");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumStat_0010    = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("graphDirGammaSpectrumStat");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumSyst_0010   = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("graphDirGammaSpectrumSyst");
        TH1D *histoPubPCMInclGammaSpectrumStat_0010                 = (TH1D*)directoryPCMGamma_0010->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumStat_0010    = new TGraphAsymmErrors(histoPubPCMInclGammaSpectrumStat_0010);
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumSyst_0010    = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("IncGammaSystError");
        TH1D *histoPubPCMInclRatioSpectrumStat_0010                 = (TH1D*)directoryPCMGamma_0010->Get("IncRatioPi0FitStatError");
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumStat_0010    = new TGraphAsymmErrors(histoPubPCMInclRatioSpectrumStat_0010);
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumSyst_0010    = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("IncRatioPi0FitSystError");
        TH1D *histoPubPCMDoubleRatioStat_0010                       = (TH1D*)directoryPCMGamma_0010->Get("DoubleRatioPi0FitStatError");
        TGraphAsymmErrors* graphPubPCMDoubleRatioStat_0010          = new TGraphAsymmErrors(histoPubPCMDoubleRatioStat_0010);
        TGraphAsymmErrors* graphPubPCMDoubleRatioSyst_0010          = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("DoubleRatioPi0FitSystError");
    TDirectoryFile* directoryPCMGamma_2040                         = (TDirectoryFile*)filePubGammaPCM->Get("Gamma_PbPb_2.76TeV_20-40%");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumStat_2040    = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("graphDirGammaSpectrumStat");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumSyst_2040   = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("graphDirGammaSpectrumSyst");
        TH1D *histoPubPCMInclGammaSpectrumStat_2040                 = (TH1D*)directoryPCMGamma_2040->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumStat_2040    = new TGraphAsymmErrors(histoPubPCMInclGammaSpectrumStat_2040);
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumSyst_2040    = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("IncGammaSystError");
        TH1D *histoPubPCMInclRatioSpectrumStat_2040                 = (TH1D*)directoryPCMGamma_2040->Get("IncRatioPi0FitStatError");
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumStat_2040    = new TGraphAsymmErrors(histoPubPCMInclRatioSpectrumStat_2040);
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumSyst_2040    = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("IncRatioPi0FitSystError");
        TH1D *histoPubPCMDoubleRatioStat_2040                       = (TH1D*)directoryPCMGamma_2040->Get("DoubleRatioPi0FitStatError");
        TGraphAsymmErrors* graphPubPCMDoubleRatioStat_2040          = new TGraphAsymmErrors(histoPubPCMDoubleRatioStat_2040);
        TGraphAsymmErrors* graphPubPCMDoubleRatioSyst_2040          = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("DoubleRatioPi0FitSystError");


    //*******************************************************************************************************************************************
    //**************************************************** Loading PHOS measurement *************************************************************
    //*******************************************************************************************************************************************
    TFile* filePHOSGamma                            = new TFile( inputFileNamePHOS.Data());
    //________________________________________________ Load PHOS 0-10% __________________________________________________________________________
    TDirectory* directoryPHOSGamma0010              = (TDirectory*)filePHOSGamma->Get("PHOS_PbPb_2760_Centrality_00-10");
        TH1D* histoPHOSDRPi0FitStatErr0010                      = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_Stat");
        TGraphAsymmErrors* graphPHOSDRPi0FitStatErr0010          = new TGraphAsymmErrors(histoPHOSDRPi0FitStatErr0010);
        TH1D* histoPHOSDRPi0FitSysErr0010                       = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_Syst");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysErr0010          = new TGraphAsymmErrors(histoPHOSDRPi0FitSysErr0010);
        TH1D* histoPHOSDRPi0FitSysAErr0010                      = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_SystA");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysAErr0010         = new TGraphAsymmErrors(histoPHOSDRPi0FitSysAErr0010);
        TH1D* histoPHOSDRPi0FitSysBcontErr0010                  = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_SystBcont");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBcontErr0010     = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcontErr0010);
        TH1D* histoPHOSDRPi0FitSysBpi0Err0010                   = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_SystBpi0");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBpi0Err0010      = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBpi0Err0010);
        TH1D* histoPHOSDRPi0FitSysBcocktErr0010                 = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_SystBcockt");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBcocktErr0010    = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcocktErr0010);
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBErr0010         = Add3ErrorsOfGraphsQuadratically (graphPHOSDRPi0FitSysBcontErr0010, graphPHOSDRPi0FitSysBpi0Err0010, graphPHOSDRPi0FitSysBcocktErr0010);
        graphPHOSDRPi0FitSysBErr0010->SetName("hPHOS_DoubleRatio_PbPb_cen00-10_SystBtot");

        TH1D* histoPHOSDRPi0FitSysCErr0010                      = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_SystC");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysCErr0010         = new TGraphAsymmErrors(histoPHOSDRPi0FitSysCErr0010);
//         TH1D* histoPHOSPi0StatErr0010                        = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_pi0_PbPb_cen00-10_Stat");
//         TH1D* histoPHOSPi0SysErr0010                         = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_pi0_PbPb_cen00-10_Syst");
//         TH1D* histoPHOSPi0SysAErr0010                        = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_pi0_PbPb_cen00-10_SystA");
//         TH1D* histoPHOSPi0SysBErr0010                        = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_pi0_PbPb_cen00-10_SystB");
//         TH1D* histoPHOSPi0SysCErr0010                        = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_pi0_PbPb_cen00-10_SystC");

        TH1D* histoPHOSIncGammaStatErr0010                      = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_gammaIncl_PbPb_cen00-10_Stat");
        TGraphAsymmErrors* graphPHOSIncGammaStatErr0010          = new TGraphAsymmErrors(histoPHOSIncGammaStatErr0010);
        TH1D* histoPHOSIncGammaSysErr0010                       = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_gammaIncl_PbPb_cen00-10_Syst");
        TGraphAsymmErrors* graphPHOSIncGammaSysErr0010          = new TGraphAsymmErrors(histoPHOSIncGammaSysErr0010);
        TH1D* histoPHOSIncGammaSysAErr0010                      = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_gammaIncl_PbPb_cen00-10_SystA");
        TGraphAsymmErrors* graphPHOSIncGammaSysAErr0010         = new TGraphAsymmErrors(histoPHOSIncGammaSysAErr0010);
        TH1D* histoPHOSIncGammaSysBErr0010                      = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_gammaIncl_PbPb_cen00-10_SystBtot");
        TGraphAsymmErrors* graphPHOSIncGammaSysBErr0010         = new TGraphAsymmErrors(histoPHOSIncGammaSysBErr0010);
        TH1D* histoPHOSIncGammaSysCErr0010                      = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_gammaIncl_PbPb_cen00-10_SystC");
        TGraphAsymmErrors* graphPHOSIncGammaSysCErr0010         = new TGraphAsymmErrors(histoPHOSIncGammaSysCErr0010);

        TH1D*  histoPHOSDirGammaStatErr0010                     = NULL;
        TH1D*  histoPHOSDirGammaSysErr0010                      = NULL;
//         TH1D*  histoPHOSDirGammaSysBErrnl0010                = NULL;
//         TH1D*  histoPHOSDirGammaSysBErrglobalE0010           = NULL;
//         TH1D*  histoPHOSDirGammaSysCErr0010                  = NULL;
        histoPHOSDirGammaStatErr0010                            = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_gammaDir_PbPb_cen00-10_Stat");
        histoPHOSDirGammaSysErr0010                             = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_gammaDir_PbPb_cen00-10_Syst");
//         histoPHOSDirGammaSysBErrnl0010                       = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_gammaDir_PbPb_cen00-10_SystBnl");
//         histoPHOSDirGammaSysBErrglobalE0010                  = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_gammaDir_PbPb_cen00-10_SystBglobalE");
//         histoPHOSDirGammaSysCErr0010                         = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_gammaDir_PbPb_cen00-10_SystC");
    //________________________________________________ Load PHOS 20-40% _________________________________________________________________________
    TDirectory* directoryPHOSGamma2040                  = (TDirectory*)filePHOSGamma->Get("PHOS_PbPb_2760_Centrality_20-40");
        TH1D* histoPHOSDRPi0FitStatErr2040                      = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_Stat");
        TGraphAsymmErrors* graphPHOSDRPi0FitStatErr2040          = new TGraphAsymmErrors(histoPHOSDRPi0FitStatErr2040);
        TH1D* histoPHOSDRPi0FitSysErr2040                       = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_Syst");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysErr2040          = new TGraphAsymmErrors(histoPHOSDRPi0FitSysErr2040);
        TH1D* histoPHOSDRPi0FitSysAErr2040                      = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystA");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysAErr2040         = new TGraphAsymmErrors(histoPHOSDRPi0FitSysAErr2040);
        TH1D* histoPHOSDRPi0FitSysBcontErr2040                  = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystBcont");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBcontErr2040     = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcontErr2040);
        TH1D* histoPHOSDRPi0FitSysBpi0Err2040                   = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystBpi0");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBpi0Err2040      = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBpi0Err2040);
        TH1D* histoPHOSDRPi0FitSysBcocktErr2040                 = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystBcockt");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBcocktErr2040    = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcocktErr2040);
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBErr2040         = Add3ErrorsOfGraphsQuadratically (graphPHOSDRPi0FitSysBcontErr2040, graphPHOSDRPi0FitSysBpi0Err2040, graphPHOSDRPi0FitSysBcocktErr2040);
        graphPHOSDRPi0FitSysBErr2040->SetName("hPHOS_DoubleRatio_PbPb_cen20-40_SystBtot");
        TH1D* histoPHOSDRPi0FitSysCErr2040                      = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystC");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysCErr2040         = new TGraphAsymmErrors(histoPHOSDRPi0FitSysCErr2040);

//         TH1D* histoPHOSPi0StatErr2040                        = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_pi0_PbPb_cen20-40_Stat");
//         TH1D* histoPHOSPi0SysErr2040                         = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_pi0_PbPb_cen20-40_Syst");
//         TH1D* histoPHOSPi0SysAErr2040                        = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_pi0_PbPb_cen20-40_SystA");
//         TH1D* histoPHOSPi0SysBErr2040                        = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_pi0_PbPb_cen20-40_SystB");
//         TH1D* histoPHOSPi0SysCErr2040                        = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_pi0_PbPb_cen20-40_SystC");

        TH1D* histoPHOSIncGammaStatErr2040                      = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaIncl_PbPb_cen20-40_Stat");
        TGraphAsymmErrors* graphPHOSIncGammaStatErr2040             = new TGraphAsymmErrors(histoPHOSIncGammaStatErr2040);
        TH1D* histoPHOSIncGammaSysErr2040                       = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaIncl_PbPb_cen20-40_Syst");
        TGraphAsymmErrors* graphPHOSIncGammaSysErr2040          = new TGraphAsymmErrors(histoPHOSIncGammaSysErr2040);
        TH1D* histoPHOSIncGammaSysAErr2040                      = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaIncl_PbPb_cen20-40_SystA");
        TGraphAsymmErrors* graphPHOSIncGammaSysAErr2040         = new TGraphAsymmErrors(histoPHOSIncGammaSysAErr2040);
        TH1D* histoPHOSIncGammaSysBErr2040                      = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaIncl_PbPb_cen20-40_SystBtot");
        TGraphAsymmErrors* graphPHOSIncGammaSysBErr2040         = new TGraphAsymmErrors(histoPHOSIncGammaSysBErr2040);
        TH1D* histoPHOSIncGammaSysCErr2040                      = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaIncl_PbPb_cen20-40_SystC");
        TGraphAsymmErrors* graphPHOSIncGammaSysCErr2040         = new TGraphAsymmErrors(histoPHOSIncGammaSysCErr2040);

        TH1D*  histoPHOSDirGammaStatErr2040                     = NULL;
        TH1D*  histoPHOSDirGammaSysErr2040                      = NULL;
//         TH1D*  histoPHOSDirGammaSysBErrnl2040                = NULL;
//         TH1D*  histoPHOSDirGammaSysBErrglobalE2040           = NULL;
//         TH1D*  histoPHOSDirGammaSysCErr2040                  = NULL;
        histoPHOSDirGammaStatErr2040                            = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaDir_PbPb_cen20-40_Stat");
        histoPHOSDirGammaSysErr2040                             = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaDir_PbPb_cen20-40_Syst");
//         histoPHOSDirGammaSysBErrnl2040                       = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaDir_PbPb_cen20-40_SystBnl");
//         histoPHOSDirGammaSysBErrglobalE2040                  = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaDir_PbPb_cen20-40_SystBglobalE");
//         histoPHOSDirGammaSysCErr2040                         = (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaDir_PbPb_cen20-40_SystC");

    //________________________________________________ Load PHOS 20-50% _________________________________________________________________________
    TDirectory* directoryPHOSGamma2050                  = (TDirectory*)filePHOSGamma->Get("PHOS_PbPb_2760_Centrality_20-40");
        TH1D* histoPHOSDRPi0FitStatErr2050                      = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_DoubleRatio_PbPb_cen20-40_Stat");
        TH1D* histoPHOSDRPi0FitSysErr2050                       = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_DoubleRatio_PbPb_cen20-40_Syst");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysErr2050          = new TGraphAsymmErrors(histoPHOSDRPi0FitSysErr2050);
        TH1D* histoPHOSDRPi0FitSysAErr2050                      = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystA");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysAErr2050         = new TGraphAsymmErrors(histoPHOSDRPi0FitSysAErr2050);
        TH1D* histoPHOSDRPi0FitSysBcontErr2050                  = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystBcont");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBcontErr2050     = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcontErr2050);
        TH1D* histoPHOSDRPi0FitSysBpi0Err2050                   = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystBpi0");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBpi0Err2050      = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBpi0Err2050);
        TH1D* histoPHOSDRPi0FitSysBcocktErr2050                 = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystBcockt");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBcocktErr2050    = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcocktErr2050);
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBErr2050         = Add3ErrorsOfGraphsQuadratically (graphPHOSDRPi0FitSysBcontErr2050, graphPHOSDRPi0FitSysBpi0Err2050, graphPHOSDRPi0FitSysBcocktErr2050);
        graphPHOSDRPi0FitSysBErr2050->SetName("hPHOS_DoubleRatio_PbPb_cen20-40_SystBtot");
        TH1D* histoPHOSDRPi0FitSysCErr2050                      = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystC");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysCErr2050         = new TGraphAsymmErrors(histoPHOSDRPi0FitSysCErr2050);

//         TH1D* histoPHOSPi0StatErr2050                        = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_pi0_PbPb_cen20-50_Stat");
//         TH1D* histoPHOSPi0SysErr2050                         = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_pi0_PbPb_cen20-50_Syst");
//         TH1D* histoPHOSPi0SysAErr2050                        = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_pi0_PbPb_cen20-50_SystA");
//         TH1D* histoPHOSPi0SysBErr2050                        = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_pi0_PbPb_cen20-50_SystB");
//         TH1D* histoPHOSPi0SysCErr2050                        = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_pi0_PbPb_cen20-50_SystC");

        TH1D* histoPHOSIncGammaStatErr2050                      = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_gammaIncl_PbPb_cen20-40_Stat");
        TH1D* histoPHOSIncGammaSysErr2050                       = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_gammaIncl_PbPb_cen20-40_Syst");
        TGraphAsymmErrors* graphPHOSIncGammaSysErr2050          = new TGraphAsymmErrors(histoPHOSIncGammaSysErr2050);
        TH1D* histoPHOSIncGammaSysAErr2050                      = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_gammaIncl_PbPb_cen20-40_SystA");
        TGraphAsymmErrors* graphPHOSIncGammaSysAErr2050         = new TGraphAsymmErrors(histoPHOSIncGammaSysAErr2050);
        TH1D* histoPHOSIncGammaSysBErr2050                      = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_gammaIncl_PbPb_cen20-40_SystBtot");
        TGraphAsymmErrors* graphPHOSIncGammaSysBErr2050         = new TGraphAsymmErrors(histoPHOSIncGammaSysBErr2050);
        TH1D* histoPHOSIncGammaSysCErr2050                      = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_gammaIncl_PbPb_cen20-40_SystC");
        TGraphAsymmErrors* graphPHOSIncGammaSysCErr2050         = new TGraphAsymmErrors(histoPHOSIncGammaSysCErr2050);
//
        TH1D*  histoPHOSDirGammaStatErr2050                     = NULL;
        TH1D*  histoPHOSDirGammaSysErr2050                      = NULL;
//         TH1D*  histoPHOSDirGammaSysBErrnl2050                = NULL;
//         TH1D*  histoPHOSDirGammaSysBErrglobalE2050           = NULL;
//         TH1D*  histoPHOSDirGammaSysCErr2050                  = NULL;
        histoPHOSDirGammaStatErr2050                            = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_gammaDir_PbPb_cen20-40_Stat");
        histoPHOSDirGammaSysErr2050                             = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_gammaDir_PbPb_cen20-40_Syst");
//         histoPHOSDirGammaSysBErrnl2050                       = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_gammaDir_PbPb_cen20-50_SystBnl");
//         histoPHOSDirGammaSysBErrglobalE2050                  = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_gammaDir_PbPb_cen20-50_SystBglobalE");
//         histoPHOSDirGammaSysCErr2050                         = (TH1D*) directoryPHOSGamma2050->Get("hPHOS_gammaDir_PbPb_cen20-50_SystC");
//

    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from PCM file *********************************************************
    //*******************************************************************************************************************************************
    //________________________________________________ Load PCM 0-10% ___________________________________________________________________________
    TDirectory* directoryPCMGamma0010               = (TDirectory*)filePCMGamma->Get("Gamma_PbPb_2.76TeV_0-10%");
        TH1D* histoPCMDRStatErr0010                             = (TH1D*) directoryPCMGamma0010->Get("DoubleRatioStatError");
        TGraphAsymmErrors* graphPCMDRStatErr0010                = new TGraphAsymmErrors(histoPCMDRStatErr0010);
        TGraphAsymmErrors* graphPCMDRSysErr0010                 = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("DoubleRatioSystError");
        TH1D* histoPCMDRPi0FitStatErr0010                       = (TH1D*) directoryPCMGamma0010->Get("DoubleRatioPi0FitStatError");
        TGraphAsymmErrors* graphPCMDRPi0FitStatErr0010        = new TGraphAsymmErrors(histoPCMDRPi0FitStatErr0010);
        graphPCMDRPi0FitStatErr0010->RemovePoint(0);
        graphPCMDRPi0FitStatErr0010->RemovePoint(0);
        graphPCMDRPi0FitStatErr0010->RemovePoint(0);
        TGraphAsymmErrors* graphPCMDRPi0FitSysErr0010           = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("DoubleRatioPi0FitSystError");
        TGraphAsymmErrors* graphPCMDRPi0FitSysAErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("DoubleRatioPi0FitSystErrorA");
        TGraphAsymmErrors* graphPCMDRPi0FitSysBErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("DoubleRatioPi0FitSystErrorB");
        TGraphAsymmErrors* graphPCMDRPi0FitSysCErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("DoubleRatioPi0FitSystErrorC");
        TH1D* histoPCMIncRStatErr0010                        = (TH1D*) directoryPCMGamma0010->Get("IncRatioStatError");
        TGraphAsymmErrors* graphPCMIncRStatErr0010        = new TGraphAsymmErrors(histoPCMIncRStatErr0010);
        TGraphAsymmErrors* graphPCMIncRSysErr0010            = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioSystError");
        while(graphPCMIncRStatErr0010->GetY()[0]==0) graphPCMIncRStatErr0010->RemovePoint(0);
        while(graphPCMIncRSysErr0010->GetY()[0]==0) graphPCMIncRSysErr0010->RemovePoint(0);
//         TGraphAsymmErrors* graphPCMIncRSysAErr0010           = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRSysBErr0010           = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRSysCErr0010           = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioSystErrorC");
        TH1D* histoPCMIncRPi0FitStatErr0010                     = (TH1D*) directoryPCMGamma0010->Get("IncRatioPi0FitStatError");
        TGraphAsymmErrors* graphPCMIncRFitPi0StatErr0010        = new TGraphAsymmErrors(histoPCMIncRPi0FitStatErr0010);
        TGraphAsymmErrors* graphPCMIncRFitPi0SysErr0010         = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioPi0FitSystError");
        while(graphPCMIncRFitPi0StatErr0010->GetY()[0]==0) graphPCMIncRFitPi0StatErr0010->RemovePoint(0);
        while(graphPCMIncRFitPi0SysErr0010->GetY()[0]==0) graphPCMIncRFitPi0SysErr0010->RemovePoint(0);
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysAErr0010     = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioPi0FitSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysBErr0010     = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioPi0FitSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysCErr0010     = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioPi0FitSystErrorC");
        TH1D* histoPCMIncGammaStatErr0010                       = (TH1D*) directoryPCMGamma0010->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPCMIncGammaStatErr0010              = new TGraphAsymmErrors(histoPCMIncGammaStatErr0010);
        while(graphPCMIncGammaStatErr0010->GetX()[0]<=1.) graphPCMIncGammaStatErr0010->RemovePoint(0);
        TGraphAsymmErrors* graphPCMIncGammaSysErr0010           = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncGammaSystError");
        TGraphAsymmErrors* graphPCMIncGammaSysAErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncGammaSystErrorA");
        TGraphAsymmErrors* graphPCMIncGammaSysBErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncGammaSystErrorB");
        TGraphAsymmErrors* graphPCMIncGammaSysCErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncGammaSystErrorC");
//         TH1D* histoPCMPi0StatErr0010                         = (TH1D*) directoryPCMGamma0010->Get("Pi0StatError");
//         TGraphAsymmErrors* graphPCMPi0SysErr0010             = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("Pi0SystError");
//         TH1D* histoPCMPi0FitStatErr0010                      = (TH1D*) directoryPCMGamma0010->Get("Pi0FitStatError");
//         TGraphAsymmErrors* graphPCMPi0FitSysErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("Pi0FitSystError");
//         TH1D* histoIncGammaCocktail0010                         = (TH1D*) directoryPCMGamma0010->Get("IncRatioCocktail");
//         TGraphAsymmErrors *graphIncGammaCocktail0010            = new TGraphAsymmErrors(histoIncGammaCocktail0010);
//         while(graphIncGammaCocktail0010->GetX()[0]<1.) graphIncGammaCocktail0010->RemovePoint(0);

        TGraphAsymmErrors* graphPCMDirGammaStatErr0010          = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSysErr0010           = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSumErrAr0010         = NULL;
        graphPCMDirGammaStatErr0010                             = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("graphDirGammaSpectrumStat");
        graphPCMDirGammaSysErr0010                              = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("graphDirGammaSpectrumSyst");
        graphPCMDirGammaSumErrAr0010                            = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("graphDirGammaSpectrumSummedAr");

        TH1D* histoDirGammaSpectrumErrorSummed_0010                  = (TH1D*)directoryPCMGamma0010->Get("histoDirGammaSpectrumErrorSummed");


    //________________________________________________ Load PCM 20-40% __________________________________________________________________________
    TDirectory* directoryPCMGamma2040               = (TDirectory*)filePCMGamma->Get("Gamma_PbPb_2.76TeV_20-40%");
        TH1D* histoPCMDRStatErr2040                          = (TH1D*) directoryPCMGamma2040->Get("DoubleRatioStatError");
        TGraphAsymmErrors* graphPCMDRSysErr2040              = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystError");
//         TGraphAsymmErrors* graphPCMDRSysAErr2040             = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystErrorA");
//         TGraphAsymmErrors* graphPCMDRSysBErr2040             = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystErrorB");
//         TGraphAsymmErrors* graphPCMDRSysCErr2040             = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystErrorC");
        TH1D* histoPCMDRPi0FitStatErr2040                       = (TH1D*) directoryPCMGamma2040->Get("DoubleRatioPi0FitStatError");
        TGraphAsymmErrors* graphPCMDRPi0FitStatErr2040        = new TGraphAsymmErrors(histoPCMDRPi0FitStatErr2040);
        graphPCMDRPi0FitStatErr2040->RemovePoint(0);
        graphPCMDRPi0FitStatErr2040->RemovePoint(0);
        graphPCMDRPi0FitStatErr2040->RemovePoint(0);
        TGraphAsymmErrors* graphPCMDRPi0FitSysErr2040           = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystError");
        TGraphAsymmErrors* graphPCMDRPi0FitSysAErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystErrorA");
        TGraphAsymmErrors* graphPCMDRPi0FitSysBErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystErrorB");
        TGraphAsymmErrors* graphPCMDRPi0FitSysCErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystErrorC");
        TH1D* histoPCMIncRStatErr2040                        = (TH1D*) directoryPCMGamma2040->Get("IncRatioStatError");
        TGraphAsymmErrors* graphPCMIncRStatErr2040        = new TGraphAsymmErrors(histoPCMIncRStatErr2040);
        TGraphAsymmErrors* graphPCMIncRSysErr2040            = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystError");
        while(graphPCMIncRStatErr2040->GetY()[0]==0) graphPCMIncRStatErr2040->RemovePoint(0);
        while(graphPCMIncRSysErr2040->GetY()[0]==0) graphPCMIncRSysErr2040->RemovePoint(0);
//         TGraphAsymmErrors* graphPCMIncRSysAErr2040           = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRSysBErr2040           = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRSysCErr2040           = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystErrorC");
        TH1D* histoPCMIncRPi0FitStatErr2040                     = (TH1D*) directoryPCMGamma2040->Get("IncRatioPi0FitStatError");
        TGraphAsymmErrors* graphPCMIncRFitPi0StatErr2040        = new TGraphAsymmErrors(histoPCMIncRPi0FitStatErr2040);
        TGraphAsymmErrors* graphPCMIncRFitPi0SysErr2040         = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystError");
        while(graphPCMIncRFitPi0StatErr2040->GetY()[0]==0) graphPCMIncRFitPi0StatErr2040->RemovePoint(0);
        while(graphPCMIncRFitPi0SysErr2040->GetY()[0]==0) graphPCMIncRFitPi0SysErr2040->RemovePoint(0);
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysAErr2040     = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysBErr2040     = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysCErr2040     = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystErrorC");
        TH1D* histoPCMIncGammaStatErr2040                       = (TH1D*) directoryPCMGamma2040->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPCMIncGammaStatErr2040              = new TGraphAsymmErrors(histoPCMIncGammaStatErr2040);
        while(graphPCMIncGammaStatErr2040->GetX()[0]<=1.) graphPCMIncGammaStatErr2040->RemovePoint(0);
        TGraphAsymmErrors* graphPCMIncGammaSysErr2040           = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystError");
        TGraphAsymmErrors* graphPCMIncGammaSysAErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystErrorA");
        TGraphAsymmErrors* graphPCMIncGammaSysBErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystErrorB");
        TGraphAsymmErrors* graphPCMIncGammaSysCErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystErrorC");
//         TH1D* histoPCMPi0StatErr2040                         = (TH1D*) directoryPCMGamma2040->Get("Pi0StatError");
//         TGraphAsymmErrors* graphPCMPi0SysErr2040             = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("Pi0SystError");
//         TH1D* histoPCMPi0FitStatErr2040                      = (TH1D*) directoryPCMGamma2040->Get("Pi0FitStatError");
//         TGraphAsymmErrors* graphPCMPi0FitSysErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("Pi0FitSystError");
//         TH1D* histoIncGammaCocktail2040                         = (TH1D*) directoryPCMGamma2040->Get("IncRatioCocktail");
//         TGraphAsymmErrors *graphIncGammaCocktail2040            = new TGraphAsymmErrors(histoIncGammaCocktail2040);

        TGraphAsymmErrors* graphPCMDirGammaStatErr2040          = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSysErr2040           = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSumErrAr2040         = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSumErrUL2040         = NULL;
        graphPCMDirGammaStatErr2040                             = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("graphDirGammaSpectrumStat");
        graphPCMDirGammaSysErr2040                              = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("graphDirGammaSpectrumSyst");
        graphPCMDirGammaSumErrAr2040                            = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("graphDirGammaSpectrumSummedAr");
        graphPCMDirGammaSumErrUL2040                            = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("graphDirGammaSpectrumSummedUL");

        TH1D* histoDirGammaSpectrumErrorSummed_2040             = (TH1D*)directoryPCMGamma2040->Get("histoDirGammaSpectrumErrorSummed");


    //________________________________________________ Load PCM 20-50% __________________________________________________________________________
    TDirectory* directoryPCMGamma2050               = (TDirectory*)filePCMGamma->Get("Gamma_PbPb_2.76TeV_20-50%");
        TH1D* histoPCMDRStatErr2050                          = (TH1D*) directoryPCMGamma2050->Get("DoubleRatioStatError");
        TGraphAsymmErrors* graphPCMDRSysErr2050              = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioSystError");
//         TGraphAsymmErrors* graphPCMDRSysAErr2050             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioSystErrorA");
//         TGraphAsymmErrors* graphPCMDRSysBErr2050             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioSystErrorB");
//         TGraphAsymmErrors* graphPCMDRSysCErr2050             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioSystErrorC");
        TH1D* histoPCMDRPi0FitStatErr2050                       = (TH1D*) directoryPCMGamma2050->Get("DoubleRatioPi0FitStatError");
        TGraphAsymmErrors* graphPCMDRPi0FitStatErr2050        = new TGraphAsymmErrors(histoPCMDRPi0FitStatErr2050);
        graphPCMDRPi0FitStatErr2050->RemovePoint(0);
        graphPCMDRPi0FitStatErr2050->RemovePoint(0);
        graphPCMDRPi0FitStatErr2050->RemovePoint(0);
        TGraphAsymmErrors* graphPCMDRPi0FitSysErr2050         = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioPi0FitSystError");
        TGraphAsymmErrors* graphPCMDRPi0FitSysAErr2050        = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioPi0FitSystErrorA");
        TGraphAsymmErrors* graphPCMDRPi0FitSysBErr2050        = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioPi0FitSystErrorB");
        TGraphAsymmErrors* graphPCMDRPi0FitSysCErr2050        = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioPi0FitSystErrorC");
        TH1D* histoPCMIncRStatErr2050                         = (TH1D*) directoryPCMGamma2050->Get("IncRatioStatError");
        TGraphAsymmErrors* graphPCMIncRStatErr2050            = new TGraphAsymmErrors(histoPCMIncRStatErr2050);
        TGraphAsymmErrors* graphPCMIncRSysErr2050             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioSystError");
        while(graphPCMIncRStatErr2050->GetY()[0]==0) graphPCMIncRStatErr2050->RemovePoint(0);
        while(graphPCMIncRSysErr2050->GetY()[0]==0) graphPCMIncRSysErr2050->RemovePoint(0);
//         TGraphAsymmErrors* graphPCMIncRSysAErr2050           = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRSysBErr2050           = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRSysCErr2050           = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioSystErrorC");
        TH1D* histoPCMIncRPi0FitStatErr2050                     = (TH1D*) directoryPCMGamma2050->Get("IncRatioPi0FitStatError");
        TGraphAsymmErrors* graphPCMIncRFitPi0StatErr2050        = new TGraphAsymmErrors(histoPCMIncRPi0FitStatErr2050);
        TGraphAsymmErrors* graphPCMIncRFitPi0SysErr2050         = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioPi0FitSystError");
        while(graphPCMIncRFitPi0StatErr2050->GetY()[0]==0) graphPCMIncRFitPi0StatErr2050->RemovePoint(0);
        while(graphPCMIncRFitPi0SysErr2050->GetY()[0]==0) graphPCMIncRFitPi0SysErr2050->RemovePoint(0);
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysAErr2050     = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioPi0FitSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysBErr2050     = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioPi0FitSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysCErr2050     = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioPi0FitSystErrorC");
        TH1D* histoPCMIncGammaStatErr2050                       = (TH1D*) directoryPCMGamma2050->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPCMIncGammaStatErr2050              = new TGraphAsymmErrors(histoPCMIncGammaStatErr2050);
        while(graphPCMIncGammaStatErr2050->GetX()[0]<=1.) graphPCMIncGammaStatErr2050->RemovePoint(0);
        TGraphAsymmErrors* graphPCMIncGammaSysErr2050           = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncGammaSystError");
        TGraphAsymmErrors* graphPCMIncGammaSysAErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncGammaSystErrorA");
        TGraphAsymmErrors* graphPCMIncGammaSysBErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncGammaSystErrorB");
        TGraphAsymmErrors* graphPCMIncGammaSysCErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncGammaSystErrorC");
//         TH1D* histoPCMPi0StatErr2050                         = (TH1D*) directoryPCMGamma2050->Get("Pi0StatError");
//         TGraphAsymmErrors* graphPCMPi0SysErr2050             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("Pi0SystError");
//         TH1D* histoPCMPi0FitStatErr2050                      = (TH1D*) directoryPCMGamma2050->Get("Pi0FitStatError");
//         TGraphAsymmErrors* graphPCMPi0FitSysErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("Pi0FitSystError");
//         TH1D* histoIncGammaCocktail2050                         = (TH1D*) directoryPCMGamma2050->Get("IncRatioCocktail");
//         TGraphAsymmErrors *graphIncGammaCocktail2050            = new TGraphAsymmErrors(histoIncGammaCocktail2050);
//         while(graphIncGammaCocktail2050->GetX()[0]<1.) graphIncGammaCocktail2050->RemovePoint(0);

        TGraphAsymmErrors* graphPCMDirGammaStatErr2050          = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSysErr2050           = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSumErrAr2050         = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSumErrUL2050         = NULL;
        graphPCMDirGammaStatErr2050                             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("graphDirGammaSpectrumStat");
        graphPCMDirGammaSysErr2050                              = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("graphDirGammaSpectrumSyst");
        graphPCMDirGammaSumErrAr2050                            = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("graphDirGammaSpectrumSummedAr");
        graphPCMDirGammaSumErrUL2050                            = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("graphDirGammaSpectrumSummedUL");

        TH1D* histoDirGammaSpectrumErrorSummed_2050             = (TH1D*)directoryPCMGamma2050->Get("histoDirGammaSpectrumErrorSummed");


        //********************************************************************************************************************************
        //********************************************************* Thermal fits *********************************************************
        //********************************************************************************************************************************
        TF1* fitThermalGammaPCMOnly0010Sum    = FitObject("e","fitThermalGammaPCMOnly0010Sum","Photon",histoDirGammaSpectrumErrorSummed_0010,0.97,1.8,NULL,"QNRMEX0+");
            //fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly0010Sum)<< endl;
        TF1* fitThermalGammaPCMOnly0010Stat   = FitObject("e","fitThermalGammaPCMOnly0010Stat","Photon",graphPCMDirGammaStatErr0010,0.97,2.3/*1.6*/,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly0010Stat)<< endl;
        TF1* fitThermalGammaPCMOnly0010Sys    = FitObject("e","fitThermalGammaPCMOnly0010Sys","Photon",graphPCMDirGammaSysErr0010,0.97,2.3/*1.6*/,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly0010Sys)<< endl;

        TF1* fitThermalGammaPCMOnly0010Stat18  = FitObject("e","fitThermalGammaPCMOnly0010Stat18","Photon",graphPCMDirGammaStatErr0010,0.97,1.8,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly0010Stat18)<< endl;
        TF1* fitThermalGammaPCMOnly0010Sys18    = FitObject("e","fitThermalGammaPCMOnly0010Sys18","Photon",graphPCMDirGammaSysErr0010,0.97,1.8,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly0010Sys18)<< endl;

        TF1* fitThermalGammaPCMOnly0010Stat20  = FitObject("e","fitThermalGammaPCMOnly0010Stat20","Photon",graphPCMDirGammaStatErr0010,0.97,2.,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly0010Stat20)<< endl;
        TF1* fitThermalGammaPCMOnly0010Sys20    = FitObject("e","fitThermalGammaPCMOnly0010Sys20","Photon",graphPCMDirGammaSysErr0010,0.97,2.,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly0010Sys20)<< endl;

        TF1* fitThermalGammaPCMOnly0010Stat22  = FitObject("e","fitThermalGammaPCMOnly0010Stat22","Photon",graphPCMDirGammaStatErr0010,0.97,2.3,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly0010Stat22)<< endl;
        TF1* fitThermalGammaPCMOnly0010Sys22    = FitObject("e","fitThermalGammaPCMOnly0010Sys22","Photon",graphPCMDirGammaSysErr0010,0.97,2.3,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly0010Sys22)<< endl;

        TGraphAsymmErrors* graphCombDirGammaSpectrumErrSumPCMOnly0010      = new TGraphAsymmErrors(histoDirGammaSpectrumErrorSummed_0010);
        TF1* fitFullDirGammaPCMOnly0010Sys            = FitObject("qcd","fitFullDirGammaPCMOnly0010Sys","Photon",graphPCMDirGammaSysErr0010,0.97,14,NULL,"QNRMEX0+");
            //fileFinalResults << WriteParameterToFile(fitFullDirGammaPCMOnly0010Sys)<< endl;
        TF1* fitFullDirGammaPCMOnly0010Stat           = FitObject("qcd","fitFullDirGammaPCMOnly0010Stat","Photon",graphPCMDirGammaStatErr0010,0.97,14,NULL,"QNRMEX0+");
            //fileFinalResults << WriteParameterToFile(fitFullDirGammaPCMOnly0010Stat)<< endl;

        // Calculate thermal spectrum
        TF1* fitPureThermalGammaPCMOnly0010Stat                     = NULL;
        TF1* fitPureThermalGammaPCMOnly0010Stat18                   = NULL;
        TF1* fitPureThermalGammaPCMOnly0010Stat20                   = NULL;
        TF1* fitPureThermalGammaPCMOnly0010Stat22                   = NULL;
        TF1* fitPureThermalGammaPCMOnly0010Sys                      = NULL;
        TF1* fitPureThermalGammaPCMOnly0010Sys18                    = NULL;
        TF1* fitPureThermalGammaPCMOnly0010Sys20                    = NULL;
        TF1* fitPureThermalGammaPCMOnly0010Sys22                    = NULL;
        TF1* fitFullThermalGammaPCMOnly0010Stat                     = NULL;
        TF1* fitFullThermalGammaPCMOnly0010Sys                      = NULL;
        TGraphAsymmErrors* graphPCMThermalGammaSpectrumStatErr0010  = NULL;
        TGraphAsymmErrors* graphRatioPCMFitThermalGammaStatErr0010  = NULL;
        TGraphAsymmErrors* graphPCMThermalGammaSpectrumSysErr0010   = NULL;
        TGraphAsymmErrors* graphRatioPCMFitThermalGammaSysErr0010   = NULL;
        TGraphAsymmErrors* graphPCMThermalGammaSpectrumSumErr0010Ar = NULL;
        if (graphPCMDirGammaStatErr0010){
            graphPCMThermalGammaSpectrumStatErr0010 = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0010, graphPCMDirGammaStatErr0010);
            graphPCMThermalGammaSpectrumStatErr0010->Print();
            fitPureThermalGammaPCMOnly0010Stat = FitObject("e","fitPureDirGamma0010Stat","Photon",graphPCMThermalGammaSpectrumStatErr0010,0.97,2.3/*1.6*/,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly0010Stat)<< endl;
            fitPureThermalGammaPCMOnly0010Stat18 = FitObject("e","fitPureDirGamma0010Stat18","Photon",graphPCMThermalGammaSpectrumStatErr0010,0.97,1.8,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly0010Stat18)<< endl;
            fitPureThermalGammaPCMOnly0010Stat20 = FitObject("e","fitPureDirGamma0010Stat20","Photon",graphPCMThermalGammaSpectrumStatErr0010,0.97,2.0,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly0010Stat20)<< endl;
            fitPureThermalGammaPCMOnly0010Stat22 = FitObject("e","fitPureDirGamma0010Stat22","Photon",graphPCMThermalGammaSpectrumStatErr0010,0.97,2.3,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly0010Stat22)<< endl;

            fitFullThermalGammaPCMOnly0010Stat = FitObject("qcd","fitFullDirGamma0010Stat","Photon",graphPCMThermalGammaSpectrumStatErr0010,0.97,14,NULL,"QNRMEX0+");
                //fileFinalResults << WriteParameterToFile(fitFullThermalGammaPCMOnly0010Stat)<< endl;
            graphRatioPCMFitThermalGammaStatErr0010 = (TGraphAsymmErrors*)graphPCMThermalGammaSpectrumStatErr0010->Clone();
            graphRatioPCMFitThermalGammaStatErr0010 = CalculateGraphErrRatioToFit(graphRatioPCMFitThermalGammaStatErr0010, fitFullThermalGammaPCMOnly0010Stat);
        }
        if (graphPCMDirGammaSysErr0010){
            graphPCMThermalGammaSpectrumSysErr0010  = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0010, graphPCMDirGammaSysErr0010);
            graphPCMThermalGammaSpectrumSysErr0010->Print();
            fitPureThermalGammaPCMOnly0010Sys = FitObject("e","fitPureDirGamma0010Sys","Photon",graphPCMThermalGammaSpectrumSysErr0010,0.97,2.3/*1.6*/,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly0010Sys)<< endl;
            fitPureThermalGammaPCMOnly0010Sys18 = FitObject("e","fitPureDirGamma0010Sys18","Photon",graphPCMThermalGammaSpectrumSysErr0010,0.97,1.8,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly0010Sys18)<< endl;
            fitPureThermalGammaPCMOnly0010Sys20 = FitObject("e","fitPureDirGamma0010Sys20","Photon",graphPCMThermalGammaSpectrumSysErr0010,0.97,2.0,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly0010Sys20)<< endl;
            fitPureThermalGammaPCMOnly0010Sys22 = FitObject("e","fitPureDirGamma0010Sys22","Photon",graphPCMThermalGammaSpectrumSysErr0010,0.97,2.3,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly0010Sys22)<< endl;
            fitFullThermalGammaPCMOnly0010Sys = FitObject("qcd","fitFullDirGamma0010Sys","Photon",graphPCMThermalGammaSpectrumSysErr0010,0.97,14,NULL,"QNRMEX0+");
                //fileFinalResults << WriteParameterToFile(fitFullThermalGammaPCMOnly0010Sys)<< endl;
            graphRatioPCMFitThermalGammaSysErr0010 = (TGraphAsymmErrors*)graphPCMThermalGammaSpectrumSysErr0010->Clone();
            graphRatioPCMFitThermalGammaSysErr0010 = CalculateGraphErrRatioToFit(graphRatioPCMFitThermalGammaSysErr0010, fitFullThermalGammaPCMOnly0010Sys);

        }
//         if (graphPCMDirGammaSumErrAr0010){
//             graphPCMThermalGammaSpectrumSumErr0010Ar  = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0010, graphPCMDirGammaSumErrAr0010, newBinsPCM, 20);
//         }
        cout << "Thermal (sum) fit for 0-10%:";
        cout << WriteParameterToFile(fitThermalGammaPCMOnly0010Stat) << endl;
        cout << WriteParameterToFile(fitThermalGammaPCMOnly0010Sys) << endl;
        TH1D* histoFitThermalGammaPCMOnly0010Stat                                                       = (TH1D*)fitThermalGammaPCMOnly0010Stat->GetHistogram();
        histoFitThermalGammaPCMOnly0010Stat->Scale(scaleFactor1);
        TH1D* histoFitPureThermalGammaPCMOnly0010Stat                                                   = NULL;
        if (graphPCMThermalGammaSpectrumStatErr0010){
            cout << "Pure thermal fit for 0-10%:";
            cout << WriteParameterToFile(fitPureThermalGammaPCMOnly0010Stat) << endl;
            cout << WriteParameterToFile(fitPureThermalGammaPCMOnly0010Sys) << endl;
            histoFitPureThermalGammaPCMOnly0010Stat                                                     = (TH1D*)fitPureThermalGammaPCMOnly0010Stat->GetHistogram();
            histoFitPureThermalGammaPCMOnly0010Stat->Scale(scaleFactor1);
        }
        SetStyleHisto(histoFitThermalGammaPCMOnly0010Stat, 3, styleFit, colorComb0010+1 );
        SetStyleHisto(histoFitPureThermalGammaPCMOnly0010Stat, 3, styleFit, colorComb0010+1 );


        TF1* fitThermalGammaPCMOnly2040Sum    = FitObject("e","fitThermalGammaPCMOnly2040Sum","Photon",histoDirGammaSpectrumErrorSummed_2040,0.97,1.8,NULL,"QNRMEX0+");
            //fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2040Sum)<< endl;
        TF1* fitThermalGammaPCMOnly2040Stat   = FitObject("e","fitThermalGammaPCMOnly2040Stat","Photon",graphPCMDirGammaStatErr2040,0.97,2.3/*1.6*/,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2040Stat)<< endl;
        TF1* fitThermalGammaPCMOnly2040Sys    = FitObject("e","fitThermalGammaPCMOnly2040Sys","Photon",graphPCMDirGammaSysErr2040,0.97,2.3/*1.6*/,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2040Sys)<< endl;

        TF1* fitThermalGammaPCMOnly2040Stat18  = FitObject("e","fitThermalGammaPCMOnly2040Stat18","Photon",graphPCMDirGammaStatErr2040,0.97,1.8,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2040Stat18)<< endl;
        TF1* fitThermalGammaPCMOnly2040Sys18    = FitObject("e","fitThermalGammaPCMOnly2040Sys18","Photon",graphPCMDirGammaSysErr2040,0.97,1.8,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2040Sys18)<< endl;

        TF1* fitThermalGammaPCMOnly2040Stat20  = FitObject("e","fitThermalGammaPCMOnly2040Stat20","Photon",graphPCMDirGammaStatErr2040,0.97,2.,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2040Stat20)<< endl;
        TF1* fitThermalGammaPCMOnly2040Sys20    = FitObject("e","fitThermalGammaPCMOnly2040Sys20","Photon",graphPCMDirGammaSysErr2040,0.97,2.,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2040Sys20)<< endl;

        TF1* fitThermalGammaPCMOnly2040Stat22  = FitObject("e","fitThermalGammaPCMOnly2040Stat22","Photon",graphPCMDirGammaStatErr2040,0.97,2.3,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2040Stat22)<< endl;
        TF1* fitThermalGammaPCMOnly2040Sys22    = FitObject("e","fitThermalGammaPCMOnly2040Sys22","Photon",graphPCMDirGammaSysErr2040,0.97,2.3,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2040Sys22)<< endl;

        TGraphAsymmErrors* graphCombDirGammaSpectrumErrSumPCMOnly2040      = new TGraphAsymmErrors(histoDirGammaSpectrumErrorSummed_2040);
        TF1* fitFullDirGammaPCMOnly2040Sys            = FitObject("qcd","fitFullDirGammaPCMOnly2040Sys","Photon",graphPCMDirGammaSysErr2040,0.97,14,NULL,"QNRMEX0+");
            //fileFinalResults << WriteParameterToFile(fitFullDirGammaPCMOnly2040Sys)<< endl;
        TF1* fitFullDirGammaPCMOnly2040Stat           = FitObject("qcd","fitFullDirGammaPCMOnly2040Stat","Photon",graphPCMDirGammaStatErr2040,0.97,14,NULL,"QNRMEX0+");
            //fileFinalResults << WriteParameterToFile(fitFullDirGammaPCMOnly2040Stat)<< endl;

        // Calculate thermal spectrum
        TF1* fitPureThermalGammaPCMOnly2040Stat                     = NULL;
        TF1* fitPureThermalGammaPCMOnly2040Stat18                   = NULL;
        TF1* fitPureThermalGammaPCMOnly2040Stat20                   = NULL;
        TF1* fitPureThermalGammaPCMOnly2040Stat22                   = NULL;
        TF1* fitPureThermalGammaPCMOnly2040Sys                      = NULL;
        TF1* fitPureThermalGammaPCMOnly2040Sys18                    = NULL;
        TF1* fitPureThermalGammaPCMOnly2040Sys20                    = NULL;
        TF1* fitPureThermalGammaPCMOnly2040Sys22                    = NULL;
        TF1* fitFullThermalGammaPCMOnly2040Stat                     = NULL;
        TF1* fitFullThermalGammaPCMOnly2040Sys                      = NULL;
        TGraphAsymmErrors* graphPCMThermalGammaSpectrumStatErr2040  = NULL;
        TGraphAsymmErrors* graphRatioPCMFitThermalGammaStatErr2040  = NULL;
        TGraphAsymmErrors* graphPCMThermalGammaSpectrumSysErr2040   = NULL;
        TGraphAsymmErrors* graphRatioPCMFitThermalGammaSysErr2040   = NULL;
        TGraphAsymmErrors* graphPCMThermalGammaSpectrumSumErr2040Ar = NULL;
        if (graphPCMDirGammaStatErr2040){
            graphPCMThermalGammaSpectrumStatErr2040 = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill2040, graphPCMDirGammaStatErr2040);
            graphPCMThermalGammaSpectrumStatErr2040->Print();
            fitPureThermalGammaPCMOnly2040Stat = FitObject("e","fitPureDirGamma2040Stat","Photon",graphPCMThermalGammaSpectrumStatErr2040,0.97,2.3/*1.6*/,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2040Stat)<< endl;
            fitPureThermalGammaPCMOnly2040Stat18 = FitObject("e","fitPureDirGamma2040Stat18","Photon",graphPCMThermalGammaSpectrumStatErr2040,0.97,1.8,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2040Stat18)<< endl;
            fitPureThermalGammaPCMOnly2040Stat20 = FitObject("e","fitPureDirGamma2040Stat20","Photon",graphPCMThermalGammaSpectrumStatErr2040,0.97,2.0,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2040Stat20)<< endl;
            fitPureThermalGammaPCMOnly2040Stat22 = FitObject("e","fitPureDirGamma2040Stat22","Photon",graphPCMThermalGammaSpectrumStatErr2040,0.97,2.3,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2040Stat22)<< endl;

            fitFullThermalGammaPCMOnly2040Stat = FitObject("qcd","fitFullDirGamma2040Stat","Photon",graphPCMThermalGammaSpectrumStatErr2040,0.97,14,NULL,"QNRMEX0+");
                //fileFinalResults << WriteParameterToFile(fitFullThermalGammaPCMOnly2040Stat)<< endl;
            graphRatioPCMFitThermalGammaStatErr2040 = (TGraphAsymmErrors*)graphPCMThermalGammaSpectrumStatErr2040->Clone();
            graphRatioPCMFitThermalGammaStatErr2040 = CalculateGraphErrRatioToFit(graphRatioPCMFitThermalGammaStatErr2040, fitFullThermalGammaPCMOnly2040Stat);
        }
        if (graphPCMDirGammaSysErr2040){
            graphPCMThermalGammaSpectrumSysErr2040  = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill2040, graphPCMDirGammaSysErr2040);
            graphPCMThermalGammaSpectrumSysErr2040->Print();
            fitPureThermalGammaPCMOnly2040Sys = FitObject("e","fitPureDirGamma2040Sys","Photon",graphPCMThermalGammaSpectrumSysErr2040,0.97,2.3/*1.6*/,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2040Sys)<< endl;
            fitPureThermalGammaPCMOnly2040Sys18 = FitObject("e","fitPureDirGamma2040Sys18","Photon",graphPCMThermalGammaSpectrumSysErr2040,0.97,1.8,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2040Sys18)<< endl;
            fitPureThermalGammaPCMOnly2040Sys20 = FitObject("e","fitPureDirGamma2040Sys20","Photon",graphPCMThermalGammaSpectrumSysErr2040,0.97,2.0,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2040Sys20)<< endl;
            fitPureThermalGammaPCMOnly2040Sys22 = FitObject("e","fitPureDirGamma2040Sys22","Photon",graphPCMThermalGammaSpectrumSysErr2040,0.97,2.3,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2040Sys22)<< endl;
            fitFullThermalGammaPCMOnly2040Sys = FitObject("qcd","fitFullDirGamma2040Sys","Photon",graphPCMThermalGammaSpectrumSysErr2040,0.97,14,NULL,"QNRMEX0+");
                //fileFinalResults << WriteParameterToFile(fitFullThermalGammaPCMOnly2040Sys)<< endl;
            graphRatioPCMFitThermalGammaSysErr2040 = (TGraphAsymmErrors*)graphPCMThermalGammaSpectrumSysErr2040->Clone();
            graphRatioPCMFitThermalGammaSysErr2040 = CalculateGraphErrRatioToFit(graphRatioPCMFitThermalGammaSysErr2040, fitFullThermalGammaPCMOnly2040Sys);

        }

//         if (graphPCMDirGammaSumErrAr2040){
//             graphPCMThermalGammaSpectrumSumErr2040Ar  = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill2040, graphPCMDirGammaSumErrAr2040, newBinsPCM, 20);
//         }
        cout << "Thermal (sum) fit for 20-40%:";
        cout << WriteParameterToFile(fitThermalGammaPCMOnly2040Stat) << endl;
        cout << WriteParameterToFile(fitThermalGammaPCMOnly2040Sys) << endl;
        TH1D* histoFitThermalGammaPCMOnly2040Stat                                                       = (TH1D*)fitThermalGammaPCMOnly2040Stat->GetHistogram();
        histoFitThermalGammaPCMOnly2040Stat->Scale(scaleFactor2);
        TH1D* histoFitPureThermalGammaPCMOnly2040Stat                                                   = NULL;
        if (graphPCMThermalGammaSpectrumStatErr2040){
            cout << "Pure thermal fit for 20-40%:";
            cout << WriteParameterToFile(fitPureThermalGammaPCMOnly2040Stat) << endl;
            cout << WriteParameterToFile(fitPureThermalGammaPCMOnly2040Sys) << endl;
            histoFitPureThermalGammaPCMOnly2040Stat                                                     = (TH1D*)fitPureThermalGammaPCMOnly2040Stat->GetHistogram();
            histoFitPureThermalGammaPCMOnly2040Stat->Scale(scaleFactor2);
        }
        SetStyleHisto(histoFitThermalGammaPCMOnly2040Stat, 3, styleFit, colorComb2040+1 );
        SetStyleHisto(histoFitPureThermalGammaPCMOnly2040Stat, 3, styleFit, colorComb2040+1 );



        TF1* fitThermalGammaPCMOnly2050Sum    = FitObject("e","fitThermalGammaPCMOnly2050Sum","Photon",histoDirGammaSpectrumErrorSummed_2050,0.97,1.8,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2050Sum)<< endl;
        TF1* fitThermalGammaPCMOnly2050Stat16  = FitObject("e","fitThermalGammaPCMOnly2050Stat16","Photon",graphPCMDirGammaStatErr2050,0.97,1.6,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2050Stat16)<< endl;
        TF1* fitThermalGammaPCMOnly2050Stat22  = FitObject("e","fitThermalGammaPCMOnly2050Sum22","Photon",graphPCMDirGammaStatErr2050,0.97,2.3,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2050Stat22)<< endl;
        TF1* fitThermalGammaPCMOnly2050Stat   = FitObject("e","fitThermalGammaPCMOnly2050Stat","Photon",graphPCMDirGammaStatErr2050,0.97,1.8,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2050Stat)<< endl;
        TF1* fitThermalGammaPCMOnly2050Sys    = FitObject("e","fitThermalGammaPCMOnly2050Sys","Photon",graphPCMDirGammaSysErr2050,0.97,1.8,NULL,"QNRMEX0+");
            fileFinalResults << WriteParameterToFile(fitThermalGammaPCMOnly2050Sys)<< endl;

        TGraphAsymmErrors* graphCombDirGammaSpectrumErrSumPCMOnly2050      = new TGraphAsymmErrors(histoDirGammaSpectrumErrorSummed_2050);
        TF1* fitFullDirGammaPCMOnly2050Sys            = FitObject("qcd","fitFullDirGammaPCMOnly2050Sys","Photon",graphPCMDirGammaSysErr2050,0.97,14,NULL,"QNRMEX0+");
            //fileFinalResults << WriteParameterToFile(fitFullDirGammaPCMOnly2050Sys)<< endl;
        TF1* fitFullDirGammaPCMOnly2050Stat           = FitObject("qcd","fitFullDirGammaPCMOnly2050Stat","Photon",graphPCMDirGammaStatErr2050,0.97,14,NULL,"QNRMEX0+");
            //fileFinalResults << WriteParameterToFile(fitFullDirGammaPCMOnly2050Stat)<< endl;

        // Calculate thermal spectrum
        TF1* fitPureThermalGammaPCMOnly2050Stat                     = NULL;
        TF1* fitPureThermalGammaPCMOnly2050Stat16                   = NULL;
        TF1* fitPureThermalGammaPCMOnly2050Stat22                   = NULL;
        TF1* fitPureThermalGammaPCMOnly2050Sys                      = NULL;
        TF1* fitFullThermalGammaPCMOnly2050Stat                     = NULL;
        TF1* fitFullThermalGammaPCMOnly2050Sys                      = NULL;
        TGraphAsymmErrors* graphPCMThermalGammaSpectrumStatErr2050  = NULL;
        TGraphAsymmErrors* graphRatioPCMFitThermalGammaStatErr2050  = NULL;
        TGraphAsymmErrors* graphPCMThermalGammaSpectrumSysErr2050   = NULL;
        TGraphAsymmErrors* graphRatioPCMFitThermalGammaSysErr2050   = NULL;
        TGraphAsymmErrors* graphPCMThermalGammaSpectrumSumErr2050Ar = NULL;
        if (graphPCMDirGammaStatErr2050){
            graphPCMThermalGammaSpectrumStatErr2050 = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill2050, graphPCMDirGammaStatErr2050);
            graphPCMThermalGammaSpectrumStatErr2050->Print();
            fitPureThermalGammaPCMOnly2050Stat = FitObject("e","fitPureDirGamma2050Stat","Photon",graphPCMThermalGammaSpectrumStatErr2050,0.97,1.8,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2050Stat)<< endl;
            fitPureThermalGammaPCMOnly2050Stat16 = FitObject("e","fitPureDirGamma2050Stat16","Photon",graphPCMThermalGammaSpectrumStatErr2050,0.97,1.6,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2050Stat16)<< endl;
            fitPureThermalGammaPCMOnly2050Stat22 = FitObject("e","fitPureDirGamma2050Stat22","Photon",graphPCMThermalGammaSpectrumStatErr2050,0.97,2.3,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2050Stat22)<< endl;

            fitFullThermalGammaPCMOnly2050Stat = FitObject("qcd","fitFullDirGamma2050Stat","Photon",graphPCMThermalGammaSpectrumStatErr2050,0.97,14,NULL,"QNRMEX0+");
                //fileFinalResults << WriteParameterToFile(fitFullThermalGammaPCMOnly2050Stat)<< endl;
            graphRatioPCMFitThermalGammaStatErr2050 = (TGraphAsymmErrors*)graphPCMThermalGammaSpectrumStatErr2050->Clone();
            graphRatioPCMFitThermalGammaStatErr2050 = CalculateGraphErrRatioToFit(graphRatioPCMFitThermalGammaStatErr2050, fitFullThermalGammaPCMOnly2050Stat);
        }
        if (graphPCMDirGammaSysErr2050){
            graphPCMThermalGammaSpectrumSysErr2050  = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill2050, graphPCMDirGammaSysErr2050);
            graphPCMThermalGammaSpectrumSysErr2050->Print();
            fitPureThermalGammaPCMOnly2050Sys = FitObject("e","fitPureDirGamma2050Sys","Photon",graphPCMThermalGammaSpectrumSysErr2050,0.97,1.8,NULL,"NRMEX0+");
                fileFinalResults << WriteParameterToFile(fitPureThermalGammaPCMOnly2050Sys)<< endl;
            fitFullThermalGammaPCMOnly2050Sys = FitObject("qcd","fitFullDirGamma2050Sys","Photon",graphPCMThermalGammaSpectrumSysErr2050,0.97,14,NULL,"QNRMEX0+");
                //fileFinalResults << WriteParameterToFile(fitFullThermalGammaPCMOnly2050Sys)<< endl;
            graphRatioPCMFitThermalGammaSysErr2050 = (TGraphAsymmErrors*)graphPCMThermalGammaSpectrumSysErr2050->Clone();
            graphRatioPCMFitThermalGammaSysErr2050 = CalculateGraphErrRatioToFit(graphRatioPCMFitThermalGammaSysErr2050, fitFullThermalGammaPCMOnly2050Sys);

        }
//         if (graphPCMDirGammaSumErrAr2050){
//             graphPCMThermalGammaSpectrumSumErr2050Ar  = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill2050, graphPCMDirGammaSumErrAr2050, newBinsPCM, 20);
//         }
        cout << "Thermal (sum) fit for 20-50%:";
        cout << WriteParameterToFile(fitThermalGammaPCMOnly2050Stat) << endl;
        cout << WriteParameterToFile(fitThermalGammaPCMOnly2050Sys) << endl;
        TH1D* histoFitThermalGammaPCMOnly2050Stat                                                       = (TH1D*)fitThermalGammaPCMOnly2050Stat->GetHistogram();
        histoFitThermalGammaPCMOnly2050Stat->Scale(1);
        TH1D* histoFitPureThermalGammaPCMOnly2050Stat                                                   = NULL;
        if (graphPCMThermalGammaSpectrumStatErr2050){
            cout << "Pure thermal fit for 20-50%:";
            cout << WriteParameterToFile(fitPureThermalGammaPCMOnly2050Stat) << endl;
            cout << WriteParameterToFile(fitPureThermalGammaPCMOnly2050Sys) << endl;
            histoFitPureThermalGammaPCMOnly2050Stat                                                     = (TH1D*)fitPureThermalGammaPCMOnly2050Stat->GetHistogram();
            histoFitPureThermalGammaPCMOnly2050Stat->Scale(1);
        }
        SetStyleHisto(histoFitThermalGammaPCMOnly2050Stat, 3, styleFit, colorComb2050+1 );
        SetStyleHisto(histoFitPureThermalGammaPCMOnly2050Stat, 3, styleFit, colorComb2050+1 );



        //********************************************************************************************************************************
        //************************************************ Compare errors PCM measurements ***********************************************
        //********************************************************************************************************************************
        //____________________________________________ Inlcusive Gamma spectrum __________________________________________________________
        TGraphAsymmErrors *statErrorIncGammaPCMPublished_0010  = NULL;
        TGraphAsymmErrors *systErrorIncGammaPCMPublished_0010  = NULL;
        TGraphAsymmErrors *statErrorIncGammaPCM_0010           = NULL;
        TGraphAsymmErrors *systErrorIncGammaPCM_0010           = NULL;
        TGraphAsymmErrors *statErrorIncGammaPHOSPublished_0010 = NULL;
        TGraphAsymmErrors *systErrorIncGammaPHOSPublished_0010 = NULL;
        statErrorIncGammaPCM_0010 = CalculateRelErrUpAsymmGraph( graphPCMIncGammaStatErr0010,"statErrorIncGammaPCM_0010");
        statErrorIncGammaPCMPublished_0010 = CalculateRelErrUpAsymmGraph( graphPubPCMInclGammaSpectrumStat_0010,"statErrorIncGammaPCMPublished_0010");
        statErrorIncGammaPHOSPublished_0010 = CalculateRelErrUpAsymmGraph( graphPHOSIncGammaStatErr0010,"statErrorIncGammaPHOSPublished_0010");
        systErrorIncGammaPCM_0010 = CalculateRelErrUpAsymmGraph( graphPCMIncGammaSysErr0010,"systErrorIncGammaPCM_0010");
        systErrorIncGammaPCMPublished_0010 = CalculateRelErrUpAsymmGraph( graphPubPCMInclGammaSpectrumSyst_0010,"systErrorIncGammaPCMPublished_0010");
        systErrorIncGammaPHOSPublished_0010 = CalculateRelErrUpAsymmGraph( graphPHOSIncGammaSysErr0010,"systErrorIncGammaPHOSPublished_0010");
        TGraphAsymmErrors *statErrorIncGammaPCMPublished_2040  = NULL;
        TGraphAsymmErrors *systErrorIncGammaPCMPublished_2040  = NULL;
        TGraphAsymmErrors *statErrorIncGammaPCM_2040           = NULL;
        TGraphAsymmErrors *systErrorIncGammaPCM_2040           = NULL;
        TGraphAsymmErrors *statErrorIncGammaPHOSPublished_2040 = NULL;
        TGraphAsymmErrors *systErrorIncGammaPHOSPublished_2040 = NULL;
        statErrorIncGammaPCM_2040 = CalculateRelErrUpAsymmGraph( graphPCMIncGammaStatErr2040,"statErrorIncGammaPCM_2040");
        statErrorIncGammaPCMPublished_2040 = CalculateRelErrUpAsymmGraph( graphPubPCMInclGammaSpectrumStat_2040,"statErrorIncGammaPCMPublished_2040");
        statErrorIncGammaPHOSPublished_2040 = CalculateRelErrUpAsymmGraph( graphPHOSIncGammaStatErr2040,"statErrorIncGammaPHOSPublished_2040");
        systErrorIncGammaPCM_2040 = CalculateRelErrUpAsymmGraph( graphPCMIncGammaSysErr2040,"systErrorIncGammaPCM_2040");
        systErrorIncGammaPCMPublished_2040 = CalculateRelErrUpAsymmGraph( graphPubPCMInclGammaSpectrumSyst_2040,"systErrorIncGammaPCMPublished_2040");
        systErrorIncGammaPHOSPublished_2040 = CalculateRelErrUpAsymmGraph( graphPHOSIncGammaSysErr2040,"systErrorIncGammaPHOSPublished_2040");


        TCanvas* canvasRelStatErr = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
        canvasRelStatErr->SetLogx();
        TH2F * histo2IncGammaelStatErrLHC11h = new TH2F("histo2IncGammaelStatErrLHC11h","histo2IncGammaelStatErrLHC11h",11000,0.23,70.,1000,0,80.5);
        SetStyleHistoTH2ForGraphs(histo2IncGammaelStatErrLHC11h, "#it{p}_{T} (GeV/#it{c})","statistical error (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2IncGammaelStatErrLHC11h->GetXaxis()->SetMoreLogLabels();
        histo2IncGammaelStatErrLHC11h->GetXaxis()->SetLabelOffset(-0.005);
        histo2IncGammaelStatErrLHC11h->GetYaxis()->SetRangeUser(0.,20.);
        histo2IncGammaelStatErrLHC11h->GetXaxis()->SetRangeUser(0.7,20.);
        histo2IncGammaelStatErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(statErrorIncGammaPCM_0010, 20, 2, kBlack, kBlack);
            statErrorIncGammaPCM_0010->Draw("p,same");
            DrawGammaSetMarkerTGraph(statErrorIncGammaPCMPublished_0010, 20, 2,  kRed, kRed);
            statErrorIncGammaPCMPublished_0010->Draw("p,same");
            DrawGammaSetMarkerTGraph(statErrorIncGammaPHOSPublished_0010, 24, 2,kGray+1, kGray+1);
            statErrorIncGammaPHOSPublished_0010->Draw("p,same");

            TLegend* legendRelStatErr = GetAndSetLegend2(0.12, 0.92-0.04*4.2*1.35, 0.5, 0.92,32);
            legendRelStatErr->SetHeader(collisionSystemCent0010.Data());
            legendRelStatErr->AddEntry(statErrorIncGammaPCM_0010,"PCM (this thesis)","p");
            legendRelStatErr->AddEntry(statErrorIncGammaPCMPublished_0010,"PCM   published","p");
            legendRelStatErr->AddEntry(statErrorIncGammaPHOSPublished_0010,"PHOS published","p");
            legendRelStatErr->AddEntry((TObject*)0,"published from PLB 754 (2016)",""); //Phys.Lett. B754 (2016) 235-248
            legendRelStatErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/IncGamma_RelStatErr2010and2011_0010.%s",outputDir.Data(),suffix.Data()));

        canvasRelStatErr->cd();
        TH2F * histo2IncGammaelSystErrLHC11h = new TH2F("histo2IncGammaelSystErrLHC11h","histo2IncGammaelSystErrLHC11h",11000,0.23,70.,1000,0,80.5);
        SetStyleHistoTH2ForGraphs(histo2IncGammaelSystErrLHC11h, "#it{p}_{T} (GeV/#it{c})","systematic error (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2IncGammaelSystErrLHC11h->GetXaxis()->SetMoreLogLabels();
        histo2IncGammaelSystErrLHC11h->GetXaxis()->SetLabelOffset(-0.005);
        histo2IncGammaelSystErrLHC11h->GetYaxis()->SetRangeUser(0.,20.);
        histo2IncGammaelSystErrLHC11h->GetXaxis()->SetRangeUser(0.7,20.);
        histo2IncGammaelSystErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(systErrorIncGammaPCM_0010, 20, 2, kBlack, kBlack);
            systErrorIncGammaPCM_0010->Draw("p,same");
            DrawGammaSetMarkerTGraph(systErrorIncGammaPCMPublished_0010, 20, 2,kRed, kRed);
            systErrorIncGammaPCMPublished_0010->Draw("p,same");
            DrawGammaSetMarkerTGraph(systErrorIncGammaPHOSPublished_0010, 24, 2,  kGray+1, kGray+1);
            systErrorIncGammaPHOSPublished_0010->Draw("p,same");

            legendRelStatErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/IncGamma_RelSystErr2010and2011_0010.%s",outputDir.Data(),suffix.Data()));

        canvasRelStatErr->SetLogx();
        histo2IncGammaelStatErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(statErrorIncGammaPCM_2040, 20, 2, kBlack, kBlack);
            statErrorIncGammaPCM_2040->Draw("p,same");
            DrawGammaSetMarkerTGraph(statErrorIncGammaPCMPublished_2040, 20, 2, kRed, kRed);
            statErrorIncGammaPCMPublished_2040->Draw("p,same");
            DrawGammaSetMarkerTGraph(statErrorIncGammaPHOSPublished_2040, 24, 2, kGray+1, kGray+1);
            statErrorIncGammaPHOSPublished_2040->Draw("p,same");

            TLegend* legendRelStatErr2 = GetAndSetLegend2(0.12, 0.92-0.04*4.2*1.35, 0.5, 0.92,32);
            legendRelStatErr2->SetHeader(collisionSystemCent2040.Data());
            legendRelStatErr2->AddEntry(statErrorIncGammaPCM_2040,"PCM (this thesis)","p");
            legendRelStatErr2->AddEntry(statErrorIncGammaPCMPublished_2040,"PCM   published","p");
            legendRelStatErr2->AddEntry(statErrorIncGammaPHOSPublished_2040,"PHOS published","p");
            legendRelStatErr2->AddEntry((TObject*)0,"published from PLB 754 (2016)",""); //Phys.Lett. B754 (2016) 235-248
            legendRelStatErr2->Draw();

        canvasRelStatErr->SaveAs(Form("%s/IncGamma_RelStatErr2010and2011_2040.%s",outputDir.Data(),suffix.Data()));

        canvasRelStatErr->cd();
        histo2IncGammaelSystErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(systErrorIncGammaPCM_2040, 20, 2, kBlack, kBlack);
            systErrorIncGammaPCM_2040->Draw("p,same");
            DrawGammaSetMarkerTGraph(systErrorIncGammaPCMPublished_2040, 20, 2,  kRed, kRed);
            systErrorIncGammaPCMPublished_2040->Draw("p,same");
            DrawGammaSetMarkerTGraph(systErrorIncGammaPHOSPublished_2040, 24, 2,kGray+1, kGray+1);
            systErrorIncGammaPHOSPublished_2040->Draw("p,same");

            legendRelStatErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/IncGamma_RelSystErr2010and2011_2040.%s",outputDir.Data(),suffix.Data()));


        //********************************************************************************************************************************
        //_______________________________________________ Inlcusive Gamma ratio __________________________________________________________
        TGraphAsymmErrors *statErrorIncRatioPCMPublished_0010  = NULL;
        TGraphAsymmErrors *systErrorIncRatioPCMPublished_0010  = NULL;
        TGraphAsymmErrors *statErrorIncRatioPCM_0010           = NULL;
        TGraphAsymmErrors *systErrorIncRatioPCM_0010           = NULL;
        statErrorIncRatioPCM_0010 = CalculateRelErrUpAsymmGraph( graphPCMIncRFitPi0StatErr0010,"statErrorIncRatioPCM_0010");
        statErrorIncRatioPCMPublished_0010 = CalculateRelErrUpAsymmGraph( graphPubPCMInclRatioSpectrumStat_0010,"statErrorIncRatioPCMPublished_0010");
        systErrorIncRatioPCM_0010 = CalculateRelErrUpAsymmGraph( graphPCMIncRFitPi0SysErr0010,"systErrorIncRatioPCM_0010");
        systErrorIncRatioPCMPublished_0010 = CalculateRelErrUpAsymmGraph( graphPubPCMInclRatioSpectrumSyst_0010,"systErrorIncRatioPCMPublished_0010");
        TGraphAsymmErrors *statErrorIncRatioPCMPublished_2040  = NULL;
        TGraphAsymmErrors *systErrorIncRatioPCMPublished_2040  = NULL;
        TGraphAsymmErrors *statErrorIncRatioPCM_2040           = NULL;
        TGraphAsymmErrors *systErrorIncRatioPCM_2040           = NULL;
        statErrorIncRatioPCM_2040 = CalculateRelErrUpAsymmGraph( graphPCMIncRFitPi0StatErr2040,"statErrorIncRatioPCM_2040");
        statErrorIncRatioPCMPublished_2040 = CalculateRelErrUpAsymmGraph( graphPubPCMInclRatioSpectrumStat_2040,"statErrorIncRatioPCMPublished_2040");
        systErrorIncRatioPCM_2040 = CalculateRelErrUpAsymmGraph( graphPCMIncRFitPi0SysErr2040,"systErrorIncRatioPCM_2040");
        systErrorIncRatioPCMPublished_2040 = CalculateRelErrUpAsymmGraph( graphPubPCMInclRatioSpectrumSyst_2040,"systErrorIncRatioPCMPublished_2040");


        canvasRelStatErr->cd();
        histo2IncGammaelStatErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(statErrorIncRatioPCM_0010, 20, 2, kBlack, kBlack);
            statErrorIncRatioPCM_0010->Draw("p,same");
            DrawGammaSetMarkerTGraph(statErrorIncRatioPCMPublished_0010, 20, 2,  kRed, kRed);
            statErrorIncRatioPCMPublished_0010->Draw("p,same");

            legendRelStatErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/IncRatio_RelStatErr2010and2011_0010.%s",outputDir.Data(),suffix.Data()));

        canvasRelStatErr->cd();
        histo2IncGammaelSystErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(systErrorIncRatioPCM_0010, 20, 2, kBlack, kBlack);
            systErrorIncRatioPCM_0010->Draw("p,same");
            DrawGammaSetMarkerTGraph(systErrorIncRatioPCMPublished_0010, 20, 2,kRed, kRed);
            systErrorIncRatioPCMPublished_0010->Draw("p,same");

            legendRelStatErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/IncRatio_RelSystErr2010and2011_0010.%s",outputDir.Data(),suffix.Data()));

        canvasRelStatErr->SetLogx();
        histo2IncGammaelStatErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(statErrorIncRatioPCM_2040, 20, 2, kBlack, kBlack);
            statErrorIncRatioPCM_2040->Draw("p,same");
            DrawGammaSetMarkerTGraph(statErrorIncRatioPCMPublished_2040, 20, 2, kRed, kRed);
            statErrorIncRatioPCMPublished_2040->Draw("p,same");

            legendRelStatErr2->Draw();

        canvasRelStatErr->SaveAs(Form("%s/IncRatio_RelStatErr2010and2011_2040.%s",outputDir.Data(),suffix.Data()));

        canvasRelStatErr->cd();
        histo2IncGammaelSystErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(systErrorIncRatioPCM_2040, 20, 2, kBlack, kBlack);
            systErrorIncRatioPCM_2040->Draw("p,same");
            DrawGammaSetMarkerTGraph(systErrorIncRatioPCMPublished_2040, 20, 2,  kRed, kRed);
            systErrorIncRatioPCMPublished_2040->Draw("p,same");

            legendRelStatErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/IncRatio_RelSystErr2010and2011_2040.%s",outputDir.Data(),suffix.Data()));


        //********************************************************************************************************************************
        TGraphAsymmErrors *statErrorDRPCMPublished_0010  = NULL;
        TGraphAsymmErrors *systErrorDRPCMPublished_0010  = NULL;
        TGraphAsymmErrors *statErrorDRPCM_0010           = NULL;
        TGraphAsymmErrors *systErrorDRPCM_0010           = NULL;
        TGraphAsymmErrors *statErrorDRPHOSPublished_0010 = NULL;
        TGraphAsymmErrors *systErrorDRPHOSPublished_0010 = NULL;
        statErrorDRPCM_0010 = CalculateRelErrUpAsymmGraph( graphPCMDRPi0FitStatErr0010,"statErrorDRPCM_0010");
        statErrorDRPCMPublished_0010 = CalculateRelErrUpAsymmGraph( graphPubPCMDoubleRatioStat_0010,"statErrorDRPCMPublished_0010");
        statErrorDRPHOSPublished_0010 = CalculateRelErrUpAsymmGraph( graphPHOSDRPi0FitStatErr0010,"statErrorDRPHOSPublished_0010");
        systErrorDRPCM_0010 = CalculateRelErrUpAsymmGraph( graphPCMDRPi0FitSysErr0010,"systErrorDRPCM_0010");
        systErrorDRPCMPublished_0010 = CalculateRelErrUpAsymmGraph( graphPubPCMDoubleRatioStat_0010,"systErrorDRPCMPublished_0010");
        systErrorDRPHOSPublished_0010 = CalculateRelErrUpAsymmGraph( graphPHOSDRPi0FitSysErr0010,"systErrorDRPHOSPublished_0010");
        TGraphAsymmErrors *statErrorDRPCMPublished_2040  = NULL;
        TGraphAsymmErrors *systErrorDRPCMPublished_2040  = NULL;
        TGraphAsymmErrors *statErrorDRPCM_2040           = NULL;
        TGraphAsymmErrors *systErrorDRPCM_2040           = NULL;
        TGraphAsymmErrors *statErrorDRPHOSPublished_2040 = NULL;
        TGraphAsymmErrors *systErrorDRPHOSPublished_2040 = NULL;
        statErrorDRPCM_2040 = CalculateRelErrUpAsymmGraph( graphPCMDRPi0FitStatErr2040,"statErrorDRPCM_2040");
        statErrorDRPCMPublished_2040 = CalculateRelErrUpAsymmGraph( graphPubPCMDoubleRatioStat_2040,"statErrorDRPCMPublished_2040");
        statErrorDRPHOSPublished_2040 = CalculateRelErrUpAsymmGraph( graphPHOSDRPi0FitStatErr2040,"statErrorDRPHOSPublished_2040");
        systErrorDRPCM_2040 = CalculateRelErrUpAsymmGraph( graphPCMDRPi0FitSysErr2040,"systErrorDRPCM_2040");
        systErrorDRPCMPublished_2040 = CalculateRelErrUpAsymmGraph( graphPubPCMDoubleRatioStat_2040,"systErrorDRPCMPublished_2040");
        systErrorDRPHOSPublished_2040 = CalculateRelErrUpAsymmGraph( graphPHOSDRPi0FitSysErr2040,"systErrorDRPHOSPublished_2040");

        canvasRelStatErr->cd();
        histo2IncGammaelStatErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(statErrorDRPCM_0010, 20, 2, kBlack, kBlack);
            statErrorDRPCM_0010->Draw("p,same");
            DrawGammaSetMarkerTGraph(statErrorDRPCMPublished_0010, 20, 2,  kRed, kRed);
            statErrorDRPCMPublished_0010->Draw("p,same");
            DrawGammaSetMarkerTGraph(statErrorDRPHOSPublished_0010, 24, 2,kGray+1, kGray+1);
            statErrorDRPHOSPublished_0010->Draw("p,same");

            legendRelStatErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr2010and2011_0010.%s",outputDir.Data(),suffix.Data()));

        canvasRelStatErr->cd();
        histo2IncGammaelSystErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(systErrorDRPCM_0010, 20, 2, kBlack, kBlack);
            systErrorDRPCM_0010->Draw("p,same");
            DrawGammaSetMarkerTGraph(systErrorDRPCMPublished_0010, 20, 2,kRed, kRed);
            systErrorDRPCMPublished_0010->Draw("p,same");
            DrawGammaSetMarkerTGraph(systErrorDRPHOSPublished_0010, 24, 2,  kGray+1, kGray+1);
            systErrorDRPHOSPublished_0010->Draw("p,same");

            legendRelStatErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/DR_RelSystErr2010and2011_0010.%s",outputDir.Data(),suffix.Data()));

        canvasRelStatErr->SetLogx();
        histo2IncGammaelStatErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(statErrorDRPCM_2040, 20, 2, kBlack, kBlack);
            statErrorDRPCM_2040->Draw("p,same");
            DrawGammaSetMarkerTGraph(statErrorDRPCMPublished_2040, 20, 2, kRed, kRed);
            statErrorDRPCMPublished_2040->Draw("p,same");
            DrawGammaSetMarkerTGraph(statErrorDRPHOSPublished_2040, 24, 2, kGray+1, kGray+1);
            statErrorDRPHOSPublished_2040->Draw("p,same");

            legendRelStatErr2->Draw();

        canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr2010and2011_2040.%s",outputDir.Data(),suffix.Data()));

        canvasRelStatErr->cd();
        histo2IncGammaelSystErrLHC11h->Draw("copy");

            DrawGammaSetMarkerTGraph(systErrorDRPCM_2040, 20, 2, kBlack, kBlack);
            systErrorDRPCM_2040->Draw("p,same");
            DrawGammaSetMarkerTGraph(systErrorDRPCMPublished_2040, 20, 2,  kRed, kRed);
            systErrorDRPCMPublished_2040->Draw("p,same");
            DrawGammaSetMarkerTGraph(systErrorDRPHOSPublished_2040, 24, 2,kGray+1, kGray+1);
            systErrorDRPHOSPublished_2040->Draw("p,same");

            legendRelStatErr->Draw();

        canvasRelStatErr->SaveAs(Form("%s/DR_RelSystErr2010and2011_2040.%s",outputDir.Data(),suffix.Data()));


    //*******************************************************************************************************************************************
    //**************************************************** Load results from PHENIX *************************************************************
    //*******************************************************************************************************************************************
    TFile* fileExperimentPbPb                       = new TFile( fileNameExperimentPbPb.Data() );
        TGraphAsymmErrors* graphPHENIXAuAuStat0020              = (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_0_20_StatErr");
        TGraphAsymmErrors* graphPHENIXAuAuSys0020               = (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_0_20_SysErr");
        TGraphAsymmErrors* graphPHENIXAuAuStat2040              = (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_20_40_StatErr");
        TGraphAsymmErrors* graphPHENIXAuAuSys2040               = (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_20_40_SysErr");


    //*******************************************************************************************************************************************
    //******************************************* DR plot with individual measurements **********************************************************
    //*******************************************************************************************************************************************
    cout << "Plotting DR individual measurements at " << __LINE__ << endl;
    Double_t arrayBoundsXIndMeasRatio[2];
    Double_t arrayBoundsYIndMeasRatio[4];
    Double_t relativeMarginsIndMeasRatioX[3];
    Double_t relativeMarginsIndMeasRatioY[3];
    ReturnCorrectValuesForCanvasScaling(1200, 1400, 1, 3, 0.09, 0.01, 0.01, 0.065, arrayBoundsXIndMeasRatio, arrayBoundsYIndMeasRatio, relativeMarginsIndMeasRatioX, relativeMarginsIndMeasRatioY);

    TCanvas * canvasRatioIndDR = new TCanvas("canvasRatioIndDR","",10,10,1200,1400);  // gives the page size
    canvasRatioIndDR->cd();

    TPad* padPartRatioInDR1 = new TPad("padPartRatioInDR1", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
    DrawGammaPadSettings( padPartRatioInDR1, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
    padPartRatioInDR1->Draw();
    TPad* padPartRatioInDR2 = new TPad("padPartRatioInDR2", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
    DrawGammaPadSettings( padPartRatioInDR2, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[1]);
    padPartRatioInDR2->Draw();
    TPad* padPartRatioInDR3 = new TPad("padPartRatioInDR3", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[3], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[2],-1, -1, -2);
    DrawGammaPadSettings( padPartRatioInDR3, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
    padPartRatioInDR3->Draw();

    //_______________________________________________________________ define text sizes _________________________________________________________
    Int_t textSizeLabelsPixel = 43;
    Double_t margin = relativeMarginsIndMeasRatioX[0]*1200;
    Double_t textsizeLabelsPad1 = 0;
    Double_t textsizeFacPad1 = 0;
    if (padPartRatioInDR1->XtoPixel(padPartRatioInDR1->GetX2()) < padPartRatioInDR1->YtoPixel(padPartRatioInDR1->GetY1())){
        textsizeLabelsPad1 = (Double_t)textSizeLabelsPixel/padPartRatioInDR1->XtoPixel(padPartRatioInDR1->GetX2()) ;
        textsizeFacPad1 = (Double_t)1./padPartRatioInDR1->XtoPixel(padPartRatioInDR1->GetX2()) ;
    } else {
        textsizeLabelsPad1 = (Double_t)textSizeLabelsPixel/padPartRatioInDR1->YtoPixel(padPartRatioInDR1->GetY1());
        textsizeFacPad1 = (Double_t)1./padPartRatioInDR1->YtoPixel(padPartRatioInDR1->GetY1());
    }
    Double_t textsizeLabelsPad2 = 0;
    Double_t textsizeFacPad2 = 0;
    if (padPartRatioInDR2->XtoPixel(padPartRatioInDR2->GetX2()) <padPartRatioInDR2->YtoPixel(padPartRatioInDR2->GetY1()) ){
        textsizeLabelsPad2 = (Double_t)textSizeLabelsPixel/padPartRatioInDR2->XtoPixel(padPartRatioInDR2->GetX2()) ;
        textsizeFacPad2 = (Double_t)1./padPartRatioInDR2->XtoPixel(padPartRatioInDR2->GetX2()) ;
    } else {
        textsizeLabelsPad2 = (Double_t)textSizeLabelsPixel/padPartRatioInDR2->YtoPixel(padPartRatioInDR2->GetY1());
        textsizeFacPad2 = (Double_t)1./padPartRatioInDR2->YtoPixel(padPartRatioInDR2->GetY1());
    }
    Double_t textsizeLabelsPad3 = 0;
    Double_t textsizeFacPad3= 0;
    if (padPartRatioInDR3->XtoPixel(padPartRatioInDR3->GetX2()) <padPartRatioInDR3->YtoPixel(padPartRatioInDR3->GetY1()) ){
        textsizeLabelsPad3 = (Double_t)textSizeLabelsPixel/padPartRatioInDR3->XtoPixel(padPartRatioInDR3->GetX2()) ;
        textsizeFacPad3 = (Double_t)1./padPartRatioInDR3->XtoPixel(padPartRatioInDR3->GetX2()) ;
    } else {
        textsizeLabelsPad3 = (Double_t)textSizeLabelsPixel/padPartRatioInDR3->YtoPixel(padPartRatioInDR3->GetY1());
        textsizeFacPad3 = (Double_t)1./padPartRatioInDR3->YtoPixel(padPartRatioInDR3->GetY1());
    }

    //_______________________________________________________________ 0-10% dummy upper panel ___________________________________________________
    TH2D *dummyDR1 ;
    dummyDR1 = new TH2D("dummyDR1", "dummyDR1", 1000, 0., 22, 1000., doubleRatio[0]+0.05, doubleRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyDR1, "#it{p}_{T} (GeV/#it{c})", "",
                            0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.95,0.10/(textsizeFacPad1*margin), 510, 505);
    dummyDR1->GetXaxis()->SetLabelOffset(-0.015);
    dummyDR1->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDR1->GetXaxis()->SetTickLength(0.06);
    dummyDR1->GetYaxis()->SetTickLength(0.028);

    //_______________________________________________________________ 20-40% dummy middle panel _________________________________________________
    TH2D *dummyDR2 ;
    dummyDR2 = new TH2D("dummyDR2", "dummyDR2", 1000, 0., 22, 1000., doubleRatio[0]+0.05, doubleRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyDR2, "#it{p}_{T} (GeV/#it{c})", "#it{R}_{#gamma}", // = (#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})
                            0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.95,0.10/(textsizeFacPad2*margin), 510, 505);
    dummyDR2->GetXaxis()->SetLabelOffset(-0.015);
    dummyDR2->GetYaxis()->CenterTitle(kTRUE);
    dummyDR2->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDR2->GetXaxis()->SetTickLength(0.06);
    dummyDR2->GetYaxis()->SetTickLength(0.028);

    //_______________________________________________________________ 20-50% dummy lower panel __________________________________________________
    TH2D *dummyDR3 ;
    dummyDR3 = new TH2D("dummyDR3", "dummyDR3", 1000, 0., 22, 1000., doubleRatio[0]+0.05, doubleRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyDR3, "#it{p}_{T} (GeV/#it{c})", "",
                            0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.92,0.10/(textsizeFacPad3*margin), 510, 505);
    dummyDR3->GetXaxis()->SetLabelOffset(-0.015);
    dummyDR3->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDR3->GetXaxis()->SetTickLength(0.055);
    dummyDR3->GetYaxis()->SetTickLength(0.035);

    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr0010,  markerStylePCM, markerSizePCM, colorPCM , colorPCM);//markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);//
        graphPCMDRPi0FitSysErr0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr0010, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        graphPHOSDRPi0FitSysErr0010->Draw("E2same");


        DrawGammaSetMarker(histoPCMDRPi0FitStatErr0010,  markerStylePCM, markerSizePCM, colorPCM , colorPCM); //markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);//
        histoPCMDRPi0FitStatErr0010->Draw("p,same,e0,X0");
        DrawGammaSetMarker(histoPHOSDRPi0FitStatErr0010, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        histoPHOSDRPi0FitStatErr0010->Draw("p,same,e0,X0");

        TLatex *labelDRCent0010 = new TLatex(0.12,0.85,collisionSystemCent0010.Data());
        SetStyleTLatex( labelDRCent0010, 0.85*textsizeLabelsPad1,4);
        labelDRCent0010->Draw();

        TLatex *labelALICEInd = new TLatex(0.82,0.85,"ALICE");
        SetStyleTLatex( labelALICEInd, 0.85*textsizeLabelsPad1,4);
        labelALICEInd->Draw();

        TLegend* legendDRIndMeas = new TLegend(0.15,0.6,0.4,0.78);
        legendDRIndMeas->SetFillStyle(0);
        legendDRIndMeas->SetFillColor(0);
        legendDRIndMeas->SetLineColor(0);
        legendDRIndMeas->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRIndMeas->SetMargin(0.2);
        legendDRIndMeas->SetTextFont(42);
        legendDRIndMeas->AddEntry(graphPCMDRPi0FitSysErr0010, "PCM","pf");
        legendDRIndMeas->AddEntry(graphPHOSDRPi0FitSysErr0010,"PHOS","pf");
        legendDRIndMeas->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2040, markerStylePCM, markerSizePCM, colorPCM , colorPCM);//markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);//
        graphPCMDRPi0FitSysErr2040->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        graphPHOSDRPi0FitSysErr2040->Draw("E2same");

        DrawGammaSetMarker(histoPCMDRPi0FitStatErr2040, markerStylePCM, markerSizePCM, colorPCM , colorPCM); // markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);//
        histoPCMDRPi0FitStatErr2040->Draw("p,same,e0,X0");
        DrawGammaSetMarker(histoPHOSDRPi0FitStatErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        histoPHOSDRPi0FitStatErr2040->Draw("p,same,e0,X0");

        TLatex *labelDRCent2040 = new TLatex(0.12,0.88,collisionSystemCent2040.Data());
        SetStyleTLatex( labelDRCent2040, 0.85*textsizeLabelsPad2,4);
        labelDRCent2040->Draw();

    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2050, markerStylePCM, markerSizePCM, colorPCM , colorPCM); //markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);//
        graphPCMDRPi0FitSysErr2050->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr2050, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        graphPHOSDRPi0FitSysErr2050->Draw("E2same");

        DrawGammaSetMarker(histoPCMDRPi0FitStatErr2050,  markerStylePCM,markerSizePCM, colorPCM , colorPCM);//markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);//
        histoPCMDRPi0FitStatErr2050->Draw("p,same,e0,X0");
        DrawGammaSetMarker(histoPHOSDRPi0FitStatErr2050, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        histoPHOSDRPi0FitStatErr2050->Draw("p,same,e0,X0");

        TLatex *labelDRCent2050 = new TLatex(0.12,0.9,collisionSystemCent2050.Data());
        SetStyleTLatex( labelDRCent2050, 0.85*textsizeLabelsPad3,4);
        labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_individualMeasurements.%s", outputDir.Data(), suffix.Data()));


    cout << "Plotting PCM DR with Theory at " << __LINE__ << endl;
    //*******************************************************************************************************************************************
    //********************************************** DR plot with PCM measurement + NLO  ********************************************************
    //*******************************************************************************************************************************************
   canvasRatioIndDR->cd();
    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();
    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
        graphPCMDRPi0FitSysErr0010->Draw("E2same");
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        histoPCMDRPi0FitStatErr0010->Draw("p,same,e0,X0");

        TLegend* legendDRPCM0010 = new TLegend(0.12,0.82-1.1*0.85*textsizeLabelsPad1*1.,0.5,0.82);
        legendDRPCM0010->SetFillStyle(0);
        legendDRPCM0010->SetFillColor(0);
        legendDRPCM0010->SetLineColor(0);
        legendDRPCM0010->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRPCM0010->SetMargin(0.2);
        legendDRPCM0010->SetTextFont(42);
        legendDRPCM0010->AddEntry(graphPCMDRPi0FitSysErr0010,"PCM (this thesis)","pf");
        legendDRPCM0010->Draw();

        labelDRCent0010->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);


        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphPCMDRPi0FitSysErr2040->Draw("E2same");
        histoPCMDRPi0FitStatErr2040->Draw("p,same,e0,X0");

        TLegend* legendDRPCM2040 = new TLegend(0.12,0.85-1.1*0.85*textsizeLabelsPad2*1.,0.5,0.85);
        legendDRPCM2040->SetFillStyle(0);
        legendDRPCM2040->SetFillColor(0);
        legendDRPCM2040->SetLineColor(0);
        legendDRPCM2040->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRPCM2040->SetMargin(0.2);
        legendDRPCM2040->SetTextFont(42);
        legendDRPCM2040->AddEntry(graphPCMDRPi0FitSysErr2040,"PCM (this thesis)","pf");
        legendDRPCM2040->Draw();
        labelDRCent2040->Draw();

    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        graphPCMDRPi0FitSysErr2050->Draw("E2same");
        histoPCMDRPi0FitStatErr2050->Draw("p,same,e0,X0");


        TLegend* legendDRPCM2050 = new TLegend(0.12,0.87-1.1*0.85*textsizeLabelsPad3*1.,0.5,0.87);
        legendDRPCM2050->SetFillStyle(0);
        legendDRPCM2050->SetFillColor(0);
        legendDRPCM2050->SetLineColor(0);
        legendDRPCM2050->SetTextSize(0.85*textsizeLabelsPad3);
        legendDRPCM2050->SetMargin(0.2);
        legendDRPCM2050->SetTextFont(42);
        legendDRPCM2050->AddEntry(graphPCMDRPi0FitSysErr2050,"PCM (this thesis)","pf");
        legendDRPCM2050->Draw();

        labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_PCMMeasurement.%s", outputDir.Data(), suffix.Data()));

    TGraphAsymmErrors* graphPCMDRPi0FitStatSysAErr0010 = AddErrorsOfGraphsQuadratically (graphPCMDRPi0FitStatErr0010, graphPCMDRPi0FitSysAErr0010);
    TGraphAsymmErrors* graphPCMDRPi0FitStatSysAErr2040 = AddErrorsOfGraphsQuadratically (graphPCMDRPi0FitStatErr2040, graphPCMDRPi0FitSysAErr2040);
    TGraphAsymmErrors* graphPCMDRPi0FitStatSysAErr2050 = AddErrorsOfGraphsQuadratically (graphPCMDRPi0FitStatErr2050, graphPCMDRPi0FitSysAErr2050);
    ProduceGraphAsymmWithoutXErrors(graphPCMDRPi0FitStatSysAErr0010);
    ProduceGraphAsymmWithoutXErrors(graphPCMDRPi0FitStatSysAErr2040);
    ProduceGraphAsymmWithoutXErrors(graphPCMDRPi0FitStatSysAErr2050);
    Double_t SysCPCMDRPi0Fit0010 = graphPCMDRPi0FitSysCErr0010->GetErrorYlow(4)/graphPCMDRPi0FitSysCErr0010->GetY()[4];
    Double_t SysCPCMDRPi0Fit2040 = graphPCMDRPi0FitSysCErr2040->GetErrorYlow(4)/graphPCMDRPi0FitSysCErr2040->GetY()[4];
    Double_t SysCPCMDRPi0Fit2050 = graphPCMDRPi0FitSysCErr2050->GetErrorYlow(4)/graphPCMDRPi0FitSysCErr2050->GetY()[4];

    canvasRatioIndDR->cd();
    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();
    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");

        TBox* boxSysCErrPCM0010 = CreateBoxConv(colorComb0010Box, 0.74, 1.-(SysCPCMDRPi0Fit0010) , 0.83, 1.+(SysCPCMDRPi0Fit0010));
        boxSysCErrPCM0010->Draw();
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysCErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010Box , colorComb0010Box, widthLinesBoxes,  kTRUE, colorComb0010Box);
        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysBErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitStatSysAErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        graphPCMDRPi0FitSysBErr0010->Draw("E2same");
        graphPCMDRPi0FitStatSysAErr0010->Draw("p,E1Z,same");

        TLegend* legendDRPCMDiffErrRepPad1 = new TLegend(0.15,0.85-1.1*0.85*textsizeLabelsPad1*4,0.35,0.85);//0.15,0.82-1.1*0.85*textsizeLabelsPad1*4,0.35,0.82);
        legendDRPCMDiffErrRepPad1->SetFillStyle(0);
        legendDRPCMDiffErrRepPad1->SetFillColor(0);
        legendDRPCMDiffErrRepPad1->SetLineColor(0);
        legendDRPCMDiffErrRepPad1->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRPCMDiffErrRepPad1->SetMargin(0.25);
        legendDRPCMDiffErrRepPad1->SetTextFont(42);
        legendDRPCMDiffErrRepPad1->SetHeader(collisionSystemCent0010.Data());//"This thesis");
        legendDRPCMDiffErrRepPad1->AddEntry(graphPCMDRPi0FitStatSysAErr0010,"Stat. #oplus Syst. A","pe");
        legendDRPCMDiffErrRepPad1->AddEntry(graphPCMDRPi0FitSysBErr0010,"Syst. B","f");
        legendDRPCMDiffErrRepPad1->AddEntry(graphPCMDRPi0FitSysCErr0010,"Syst. C","f");
        legendDRPCMDiffErrRepPad1->Draw();

//         labelDRCent0010->Draw();
    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");

        TBox* boxSysCErrPCM2040 = CreateBoxConv(colorComb2040Box, 0.74, 1.-(SysCPCMDRPi0Fit2040) , 0.83, 1.+(SysCPCMDRPi0Fit2040));
        boxSysCErrPCM2040->Draw();
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysCErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040Box , colorComb2040Box, widthLinesBoxes,  kTRUE, colorComb2040Box);
        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysBErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitStatSysAErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphPCMDRPi0FitSysBErr2040->Draw("E2same");
        graphPCMDRPi0FitStatSysAErr2040->Draw("p,E1Z,same");

        TLegend* legendDRPCMDiffErrRepPad2 = new TLegend(0.15,0.88-1.1*0.85*textsizeLabelsPad2*4,0.35,0.88);//0.15,0.85-1.1*0.85*textsizeLabelsPad2*4,0.35,0.85);
        legendDRPCMDiffErrRepPad2->SetFillStyle(0);
        legendDRPCMDiffErrRepPad2->SetFillColor(0);
        legendDRPCMDiffErrRepPad2->SetLineColor(0);
        legendDRPCMDiffErrRepPad2->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRPCMDiffErrRepPad2->SetMargin(0.25);
        legendDRPCMDiffErrRepPad2->SetTextFont(42);
        legendDRPCMDiffErrRepPad2->SetHeader(collisionSystemCent2040.Data()); //"This thesis");
        legendDRPCMDiffErrRepPad2->AddEntry(graphPCMDRPi0FitStatSysAErr2040,"Stat. #oplus Syst. A","pe");
        legendDRPCMDiffErrRepPad2->AddEntry(graphPCMDRPi0FitSysBErr2040,"Syst. B","f");
        legendDRPCMDiffErrRepPad2->AddEntry(graphPCMDRPi0FitSysCErr2040,"Syst. C","f");
        legendDRPCMDiffErrRepPad2->Draw();

//         labelDRCent2040->Draw();
    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");

        TBox* boxSysCErrPCM2050 = CreateBoxConv(colorComb2050Box, 0.74, 1.-(SysCPCMDRPi0Fit2050) , 0.83, 1.+(SysCPCMDRPi0Fit2050));
        boxSysCErrPCM2050->Draw();
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysCErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050Box , colorComb2050Box, widthLinesBoxes,  kTRUE, colorComb2050Box);
        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysBErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitStatSysAErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        graphPCMDRPi0FitSysBErr2050->Draw("E2same");
        graphPCMDRPi0FitStatSysAErr2050->Draw("p,E1Z,same");

        TLegend* legendDRPCMDiffErrRepPad3 = new TLegend(0.15,0.9-1.1*0.85*textsizeLabelsPad3*4,0.35,0.9);//0.15,0.87-1.1*0.85*textsizeLabelsPad3*4,0.35,0.87);
        legendDRPCMDiffErrRepPad3->SetFillStyle(0);
        legendDRPCMDiffErrRepPad3->SetFillColor(0);
        legendDRPCMDiffErrRepPad3->SetLineColor(0);
        legendDRPCMDiffErrRepPad3->SetTextSize(0.85*textsizeLabelsPad3);
        legendDRPCMDiffErrRepPad3->SetMargin(0.25);
        legendDRPCMDiffErrRepPad3->SetTextFont(42);
        legendDRPCMDiffErrRepPad3->SetHeader(collisionSystemCent2050.Data());//"This thesis");
        legendDRPCMDiffErrRepPad3->AddEntry(graphPCMDRPi0FitStatSysAErr2050,"Stat. #oplus Syst. A","pe");
        legendDRPCMDiffErrRepPad3->AddEntry(graphPCMDRPi0FitSysBErr2050,"Syst. B","f");
        legendDRPCMDiffErrRepPad3->AddEntry(graphPCMDRPi0FitSysCErr2050,"Syst. C","f");
        legendDRPCMDiffErrRepPad3->Draw();

//         labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_PCMMeasurement_diffErrRep.%s", outputDir.Data(), suffix.Data()));


   canvasRatioIndDR->cd();
    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();
    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        graphPCMDRPi0FitSysErr0010->Draw("E2same");
        histoPCMDRPi0FitStatErr0010->Draw("p,same,e0,X0");
        DrawGammaSetMarkerTGraphAsym(graphPCMDRSysErr0010, markerStyleComb0010+4, markerSizeComb0010, colorCombNotFit0010 , colorCombNotFit0010,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRStatErr0010, markerStyleComb0010+4, markerSizeComb0010, colorCombNotFit0010 , colorCombNotFit0010);
        graphPCMDRSysErr0010->Draw("E2same");
        histoPCMDRStatErr0010->Draw("p,same,e0,X0");

        TLatex *labelThesisPad1 = new TLatex(0.82,0.07,"This thesis");
        SetStyleTLatex( labelThesisPad1, 0.85*textsizeLabelsPad1*1.05,4);
        labelThesisPad1->Draw();

        TLegend* legendDRPCMCompToFit0010 = new TLegend(0.12,0.82-1.1*0.85*textsizeLabelsPad1*2.7,0.5,0.82);
        legendDRPCMCompToFit0010->SetFillStyle(0);
        legendDRPCMCompToFit0010->SetFillColor(0);
        legendDRPCMCompToFit0010->SetLineColor(0);
        legendDRPCMCompToFit0010->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRPCMCompToFit0010->SetMargin(0.2);
        legendDRPCMCompToFit0010->SetTextFont(42);
        legendDRPCMCompToFit0010->AddEntry(graphPCMDRPi0FitSysErr0010,"#pi^{0} fitted","pf");
        legendDRPCMCompToFit0010->AddEntry(graphPCMDRSysErr0010,"#pi^{0} measured","pf");
        legendDRPCMCompToFit0010->Draw();
        labelDRCent0010->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);


        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphPCMDRPi0FitSysErr2040->Draw("E2same");
        histoPCMDRPi0FitStatErr2040->Draw("p,same,e0,X0");
        DrawGammaSetMarkerTGraphAsym(graphPCMDRSysErr2040, markerStyleComb2040-6, markerSizeComb2040, colorCombNotFit2040 , colorCombNotFit2040,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRStatErr2040, markerStyleComb2040-6, markerSizeComb2040, colorCombNotFit2040 , colorCombNotFit2040);
        graphPCMDRSysErr2040->Draw("E2same");
        histoPCMDRStatErr2040->Draw("p,same,e0,X0");

        TLatex *labelThesisPad2 = new TLatex(0.82,0.07,"This thesis");
        SetStyleTLatex( labelThesisPad2, 0.85*textsizeLabelsPad2*1.05,4);
        labelThesisPad2->Draw();

        TLegend* legendDRPCMCompToFit2040 = new TLegend(0.12,0.85-1.1*0.85*textsizeLabelsPad2*2.7,0.5,0.85);
        legendDRPCMCompToFit2040->SetFillStyle(0);
        legendDRPCMCompToFit2040->SetFillColor(0);
        legendDRPCMCompToFit2040->SetLineColor(0);
        legendDRPCMCompToFit2040->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRPCMCompToFit2040->SetMargin(0.2);
        legendDRPCMCompToFit2040->SetTextFont(42);
        legendDRPCMCompToFit2040->AddEntry(graphPCMDRPi0FitSysErr2040,"#pi^{0} fitted","pf");
        legendDRPCMCompToFit2040->AddEntry(graphPCMDRSysErr2040,"#pi^{0} measured","pf");
        legendDRPCMCompToFit2040->Draw();
        labelDRCent2040->Draw();

    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        graphPCMDRPi0FitSysErr2050->Draw("E2same");
        histoPCMDRPi0FitStatErr2050->Draw("p,same,e0,X0");
        DrawGammaSetMarkerTGraphAsym(graphPCMDRSysErr2050, markerStyleComb2050-6, markerSizeComb2050, colorCombNotFit2050 , colorCombNotFit2050,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRStatErr2050, markerStyleComb2050-6, markerSizeComb2050, colorCombNotFit2050 , colorCombNotFit2050);
        graphPCMDRSysErr2050->Draw("E2same");
        histoPCMDRStatErr2050->Draw("p,same,e0,X0");

        TLatex *labelThesisPad3 = new TLatex(0.82,0.23,"This thesis");
        SetStyleTLatex( labelThesisPad3, 0.85*textsizeLabelsPad3*1.05,4);
        labelThesisPad3->Draw();

        TLegend* legendDRPCMCompToFit2050 = new TLegend(0.12,0.87-1.1*0.85*textsizeLabelsPad3*2.7,0.5,0.87);
        legendDRPCMCompToFit2050->SetFillStyle(0);
        legendDRPCMCompToFit2050->SetFillColor(0);
        legendDRPCMCompToFit2050->SetLineColor(0);
        legendDRPCMCompToFit2050->SetTextSize(0.85*textsizeLabelsPad3);
        legendDRPCMCompToFit2050->SetMargin(0.2);
        legendDRPCMCompToFit2050->SetTextFont(42);
        legendDRPCMCompToFit2050->AddEntry(graphPCMDRPi0FitSysErr2050,"#pi^{0} fitted","pf");
        legendDRPCMCompToFit2050->AddEntry(graphPCMDRSysErr2050,"#pi^{0} measured","pf");
        legendDRPCMCompToFit2050->Draw();
        labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_PCMMeasurement_Pi0MeasAndFit.%s", outputDir.Data(), suffix.Data()));

    canvasRatioIndDR->cd();
    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();
    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR0010, 3, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryEPS09DR0010, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryCT10DR0010, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
        graphTheoryEPS09DR0010->Draw("p3lsame");
        graphTheoryCT10DR0010->Draw("p3lsame");
        graphTheoryNLODR0010->Draw("p3lsame");

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
        graphPCMDRPi0FitSysErr0010->Draw("E2same");
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        histoPCMDRPi0FitStatErr0010->Draw("p,same,e0,X0");

        TLegend* legendDRPCMNLO0010 = new TLegend(0.12,0.82-1.1*0.85*textsizeLabelsPad1*4.7,0.5,0.82);
        legendDRPCMNLO0010->SetFillStyle(0);
        legendDRPCMNLO0010->SetFillColor(0);
        legendDRPCMNLO0010->SetLineColor(0);
        legendDRPCMNLO0010->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRPCMNLO0010->SetMargin(0.2);
        legendDRPCMNLO0010->SetTextFont(42);
        legendDRPCMNLO0010->AddEntry(graphPCMDRPi0FitSysErr0010,"PCM","pf"); // (this thesis)
        legendDRPCMNLO0010->AddEntry(graphTheoryNLODR0010,"NLO pQCD PDF: CTEQ6M5 FF: GRV ","l");
        legendDRPCMNLO0010->AddEntry(graphTheoryCT10DR0010,"JETPHOX PDF: CT10, FF: BFG2","f");
        legendDRPCMNLO0010->AddEntry(graphTheoryEPS09DR0010,"JETPHOX nPDF: EPS09, FF: BFG2","f");
        legendDRPCMNLO0010->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
        legendDRPCMNLO0010->Draw();

        labelDRCent0010->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR2040, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryCT10DR2040, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryEPS09DR2040, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryEPS09DR2040->Draw("p3,l,same");
        graphTheoryCT10DR2040->Draw("p3,l,same");
        graphTheoryNLODR2040->Draw("p3lsame");

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphPCMDRPi0FitSysErr2040->Draw("E2same");
        histoPCMDRPi0FitStatErr2040->Draw("p,same,e0,X0");

        TLegend* legendDRPCMNLO2040 = new TLegend(0.12,0.85-1.1*0.85*textsizeLabelsPad2*4.7,0.5,0.85);
        legendDRPCMNLO2040->SetFillStyle(0);
        legendDRPCMNLO2040->SetFillColor(0);
        legendDRPCMNLO2040->SetLineColor(0);
        legendDRPCMNLO2040->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRPCMNLO2040->SetMargin(0.2);
        legendDRPCMNLO2040->SetTextFont(42);
        legendDRPCMNLO2040->AddEntry(graphPCMDRPi0FitSysErr2040,"PCM","pf"); // (this thesis)
        legendDRPCMNLO2040->AddEntry(graphTheoryNLODR2040,"NLO pQCD PDF: CTEQ6M5 FF: GRV ","l");
        legendDRPCMNLO2040->AddEntry(graphTheoryCT10DR2040,"JETPHOX PDF: CT10, FF: BFG2","f");
        legendDRPCMNLO2040->AddEntry(graphTheoryEPS09DR2040,"JETPHOX nPDF: EPS09, FF: BFG2","f");
        legendDRPCMNLO2040->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
        legendDRPCMNLO2040->Draw();

        labelDRCent2040->Draw();

    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR2050, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryCT10DR2050, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryEPS09DR2050, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryEPS09DR2050->Draw("p3lsame");
        graphTheoryCT10DR2050->Draw("p3lsame");
        graphTheoryNLODR2050->Draw("p3lsame");

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        graphPCMDRPi0FitSysErr2050->Draw("E2same");
        histoPCMDRPi0FitStatErr2050->Draw("p,same,e0,X0");


        TLegend* legendDRPCMNLO2050 = new TLegend(0.12,0.87-1.1*0.85*textsizeLabelsPad3*4.7,0.5,0.87);
        legendDRPCMNLO2050->SetFillStyle(0);
        legendDRPCMNLO2050->SetFillColor(0);
        legendDRPCMNLO2050->SetLineColor(0);
        legendDRPCMNLO2050->SetTextSize(0.85*textsizeLabelsPad3);
        legendDRPCMNLO2050->SetMargin(0.2);
        legendDRPCMNLO2050->SetTextFont(42);
        legendDRPCMNLO2050->AddEntry(graphPCMDRPi0FitSysErr2050,"PCM","pf"); // (this thesis)
        legendDRPCMNLO2050->AddEntry(graphTheoryNLODR2050,"NLO pQCD PDF: CTEQ6M5 FF: GRV ","l");
        legendDRPCMNLO2050->AddEntry(graphTheoryCT10DR2050,"JETPHOX PDF: CT10, FF: BFG2","f");
        legendDRPCMNLO2050->AddEntry(graphTheoryEPS09DR2050,"JETPHOX nPDF: EPS09, FF: BFG2","f");
        legendDRPCMNLO2050->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
        legendDRPCMNLO2050->Draw();

        labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_PCMMeasurementTheory.%s", outputDir.Data(), suffix.Data()));

    //_______________________________________________________________ 0-10% dummy upper panel ___________________________________________________
    TH2D *dummyDRwithModels1 = new TH2D("dummyDRwithModels1", "dummyDRwithModels1", 1000, 0., 22, 1000., doubleRatio[0]+0.1, doubleRatio[1]+0.25);
    SetStyleHistoTH2ForGraphs( dummyDRwithModels1, "#it{p}_{T} (GeV/#it{c})", "",0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.95,0.10/(textsizeFacPad1*margin), 510, 505);
    dummyDRwithModels1->GetXaxis()->SetLabelOffset(-0.015);
    dummyDRwithModels1->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDRwithModels1->GetXaxis()->SetTickLength(0.06);
    dummyDRwithModels1->GetYaxis()->SetTickLength(0.028);

    //_______________________________________________________________ 20-40% dummy middle panel _________________________________________________
    TH2D *dummyDRwithModels2 = new TH2D("dummyDRwithModels2", "dummyDRwithModels2", 1000, 0., 22, 1000., doubleRatio[0]+0.1, doubleRatio[1]+0.25);
    SetStyleHistoTH2ForGraphs( dummyDRwithModels2, "#it{p}_{T} (GeV/#it{c})", "#it{R}_{#gamma}",0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.95,0.10/(textsizeFacPad2*margin), 510, 505);
    dummyDRwithModels2->GetXaxis()->SetLabelOffset(-0.015);
    dummyDRwithModels2->GetYaxis()->CenterTitle(kTRUE);
    dummyDRwithModels2->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDRwithModels2->GetXaxis()->SetTickLength(0.06);
    dummyDRwithModels2->GetYaxis()->SetTickLength(0.028);

    //_______________________________________________________________ 20-50% dummy lower panel __________________________________________________
    TH2D *dummyDRwithModels3 = new TH2D("dummyDRwithModels3", "dummyDRwithModels3", 1000, 0., 22, 1000., doubleRatio[0]+0.1, doubleRatio[1]+0.25);
    SetStyleHistoTH2ForGraphs( dummyDRwithModels3, "#it{p}_{T} (GeV/#it{c})", "",0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.92,0.10/(textsizeFacPad3*margin), 510, 505);
    dummyDRwithModels3->GetXaxis()->SetLabelOffset(-0.015);
    dummyDRwithModels3->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDRwithModels3->GetXaxis()->SetTickLength(0.055);
    dummyDRwithModels3->GetYaxis()->SetTickLength(0.035);

    canvasRatioIndDR->cd();
    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();
    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDRwithModels1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

//         graphTheoryEPS09DR0010->Draw("p3lsame");
//         graphTheoryCT10DR0010->Draw("p3lsame");
//         graphTheoryNLODR0010->Draw("p3lsame");

        TLegend* legendDRPCMOnlyTheory1 = new TLegend(0.12,0.82-1.1*0.85*textsizeLabelsPad1*6,0.5,0.82);
        legendDRPCMOnlyTheory1->SetFillStyle(0);
        legendDRPCMOnlyTheory1->SetFillColor(0);
        legendDRPCMOnlyTheory1->SetLineColor(0);
        //legendDRPCMOnlyTheory1->SetNColumns(2);
        legendDRPCMOnlyTheory1->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRPCMOnlyTheory1->SetMargin(0.15);
        legendDRPCMOnlyTheory1->SetTextFont(42);
        legendDRPCMOnlyTheory1->AddEntry(graphPCMDRPi0FitSysErr0010,"PCM","pf"); // (this thesis)
        legendDRPCMOnlyTheory1->AddEntry(graphTheoryMcGill0010Plot,"McGill group,      in preparation","l"); //"Paquet et al.,      PRC 93 (2016) 044906","l");
        legendDRPCMOnlyTheory1->AddEntry(graphTheoryPHSDplusPrompt0010Plot,"Linnyk et al.,       PRC 92 (2015) 054914","l");
        legendDRPCMOnlyTheory1->AddEntry(graphTheoryRapp0010Plot,"v. Hees et al.,     NPA 933 (2015) 256","l");
        legendDRPCMOnlyTheory1->AddEntry(graphTheoryChatterjee0010Plot,"Chatterjee et al.,","l");
        legendDRPCMOnlyTheory1->AddEntry((TObject*)0,"PRC 85 (2012) 064910","");
        legendDRPCMOnlyTheory1->AddEntry((TObject*)0,"+ arXiv: 1305.0624","");
        legendDRPCMOnlyTheory1->Draw();
        labelDRCent0010->Draw();

        graphTheoryDRChatterjeeSummed0010->Draw("p3lsame");
        graphTheoryDRPHSD0010->Draw("p3lsame");
        graphTheoryDRMcGill0010->Draw("p3lsame");
        //graphTheoryDRMcGillfromTF10010->Draw("p3lsame");
        graphTheoryDRRapp0010->Draw("p3lsame");
//         graphTheoryDRChatterjee0010->Draw("p3lsame");
//         graphTheoryDRChatterjeePrompt0010->Draw("p3lsame");
//         graphTheoryDRChatterjeeThermal0010->Draw("p3lsame");

        graphPCMDRPi0FitSysErr0010->Draw("E2same");
        histoPCMDRPi0FitStatErr0010->Draw("p,same,e0,X0");

        dummyDRwithModels1->Draw("axis,same");
    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDRwithModels2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

//         graphTheoryEPS09DR2040->Draw("p3,l,same");
//         graphTheoryCT10DR2040->Draw("p3,l,same");
//         graphTheoryNLODR2040->Draw("p3lsame");
        graphTheoryDRChatterjee2040->Draw("p3lsame");
        graphTheoryDRRapp2040->Draw("p3lsame");
        graphTheoryDRPHSD2040->Draw("p3lsame");
        graphTheoryDRMcGill2040->Draw("p3lsame");
        //graphTheoryDRMcGillfromTF12040->Draw("p3lsame");

        graphPCMDRPi0FitSysErr2040->Draw("E2same");
        histoPCMDRPi0FitStatErr2040->Draw("p,same,e0,X0");

        TLegend* legendDRPCMOnlyTheory2 = new TLegend(0.12,0.85-1.1*0.85*textsizeLabelsPad2*5.5,0.5,0.85);
        legendDRPCMOnlyTheory2->SetFillStyle(0);
        legendDRPCMOnlyTheory2->SetFillColor(0);
        legendDRPCMOnlyTheory2->SetLineColor(0);
        //legendDRPCMOnlyTheory2->SetNColumns(2);
        legendDRPCMOnlyTheory2->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRPCMOnlyTheory2->SetMargin(0.15);
        legendDRPCMOnlyTheory2->SetTextFont(42);
        legendDRPCMOnlyTheory2->AddEntry(graphPCMDRPi0FitSysErr2040,"PCM","pf"); // (this thesis)
        legendDRPCMOnlyTheory2->AddEntry(graphTheoryMcGill0010Plot,"Paquet et al.,      PRC 93 (2016) 044906","l");
        legendDRPCMOnlyTheory2->AddEntry(graphTheoryPHSDplusPrompt0010Plot,"Linnyk et al.,       PRC 92 (2015) 054914","l");
        legendDRPCMOnlyTheory2->AddEntry(graphTheoryRapp0010Plot,"v. Hees et al.,     NPA 933 (2015) 256","l");
        legendDRPCMOnlyTheory2->AddEntry(graphTheoryChatterjee0010Plot,"Chatterjee et al., PRC 85 (2012) 064910","l");
        legendDRPCMOnlyTheory2->AddEntry((TObject*)0,"                              + arXiv: 1305.0624","");
        legendDRPCMOnlyTheory2->Draw();
        labelDRCent2040->Draw();

        dummyDRwithModels2->Draw("axis,same");
    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDRwithModels3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

//         graphTheoryEPS09DR2050->Draw("p3lsame");
//         graphTheoryCT10DR2050->Draw("p3lsame");
//         graphTheoryNLODR2050->Draw("p3lsame");
        graphTheoryDRChatterjeeSummed2050->Draw("p3lsame");
        graphTheoryDRPHSD2050->Draw("p3lsame");
        graphTheoryDRMcGill2050->Draw("p3lsame");
        //graphTheoryDRMcGillfromTF12050->Draw("p3lsame");
//         graphTheoryDRChatterjee2050->Draw("p3lsame");
//         graphTheoryDRChatterjeePrompt2050->Draw("p3lsame");
//         graphTheoryDRChatterjeeThermal2050->Draw("p3lsame");

        graphPCMDRPi0FitSysErr2050->Draw("E2same");
        histoPCMDRPi0FitStatErr2050->Draw("p,same,e0,X0");

        TLegend* legendDRPCMOnlyTheory3 = new TLegend(0.12,0.87-1.1*0.85*textsizeLabelsPad3*4.5,0.5,0.87);
        legendDRPCMOnlyTheory3->SetFillStyle(0);
        legendDRPCMOnlyTheory3->SetFillColor(0);
        legendDRPCMOnlyTheory3->SetLineColor(0);
        //legendDRPCMOnlyTheory3->SetNColumns(2);
        legendDRPCMOnlyTheory3->SetTextSize(0.85*textsizeLabelsPad3);
        legendDRPCMOnlyTheory3->SetMargin(0.15);
        legendDRPCMOnlyTheory3->SetTextFont(42);
        legendDRPCMOnlyTheory3->AddEntry(graphPCMDRPi0FitSysErr2050,"PCM","pf"); // (this thesis)
        legendDRPCMOnlyTheory3->AddEntry(graphTheoryMcGill0010Plot,"McGill group,      in preparation","l"); //""Paquet et al.,      PRC 93 (2016) 044906","l");
        legendDRPCMOnlyTheory3->AddEntry(graphTheoryPHSDplusPrompt0010Plot,"Linnyk et al.,       PRC 92 (2015) 054914","l");
        //legendDRPCMOnlyTheory3->AddEntry(graphTheoryRapp0010Plot,"v. Hees et al.,     NPA 933 (2015) 256","l");
        legendDRPCMOnlyTheory3->AddEntry(graphTheoryChatterjee0010Plot,"Chatterjee et al., PRC 85 (2012) 064910","l");
        legendDRPCMOnlyTheory3->AddEntry((TObject*)0,"                              + arXiv: 1305.0624","");
        legendDRPCMOnlyTheory3->Draw();
        labelDRCent2050->Draw();

        dummyDRwithModels3->Draw("axis,same");
    canvasRatioIndDR->SaveAs(Form("%s/DR_PCMMeasurementNLOplusModels.%s", outputDir.Data(), suffix.Data()));


    cout << "Plotting inclusive gamma PCM only " << __LINE__ << endl;
    //*******************************************************************************************************************************************
    //*************************************************** Plotting inclusive Gamma Spectrum PCM only ***********************************************
    //*******************************************************************************************************************************************
    TGraphAsymmErrors* graphPCMIncGammaSysErr0010Plot = ScaleGraph(graphPCMIncGammaSysErr0010,100);
    TGraphAsymmErrors* graphPCMIncGammaStatErr0010Plot = ScaleGraph(graphPCMIncGammaStatErr0010,100);
    TGraphAsymmErrors* graphPCMIncGammaSysErr2040Plot = ScaleGraph(graphPCMIncGammaSysErr2040,10);
    TGraphAsymmErrors* graphPCMIncGammaStatErr2040Plot = ScaleGraph(graphPCMIncGammaStatErr2040,10);
    TGraphAsymmErrors* graphPCMIncGammaSysErr2050Plot = ScaleGraph(graphPCMIncGammaSysErr2050,1);
    TGraphAsymmErrors* graphPCMIncGammaStatErr2050Plot = ScaleGraph(graphPCMIncGammaStatErr2050,1);

    TCanvas *canvasIncGammaIndMeas = new TCanvas("canvasIncGammaIndMeas","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasIncGammaIndMeas, 0.16, 0.01, 0.01, 0.07);
    canvasIncGammaIndMeas->SetLogy();
    canvasIncGammaIndMeas->SetLogx();

    Int_t textSizeLabelsPixelIncGam = 48;
    Double_t textsizeLabelsIncGamma = 0;
    if (canvasIncGammaIndMeas->XtoPixel(canvasIncGammaIndMeas->GetX2()) < canvasIncGammaIndMeas->YtoPixel(canvasIncGammaIndMeas->GetY1())){
        textsizeLabelsIncGamma = (Double_t)textSizeLabelsPixelIncGam/canvasIncGammaIndMeas->XtoPixel(canvasIncGammaIndMeas->GetX2()) ;
    } else {
        textsizeLabelsIncGamma = (Double_t)textSizeLabelsPixelIncGam/canvasIncGammaIndMeas->YtoPixel(canvasIncGammaIndMeas->GetY1());
    }

    TH2D *dummyGammaIndMeas = new TH2D("dummyGammaIndMeas", "dummyGammaIndMeas", 1000, 0., 22, 1000., 6e-8,8e3);
    SetStyleHistoTH2ForGraphs( dummyGammaIndMeas, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{inc}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})", 0.85*textsizeLabelsIncGamma, textsizeLabelsIncGamma, 0.85*textsizeLabelsIncGamma, textsizeLabelsIncGamma, 0.75, 1.6);
    dummyGammaIndMeas->GetXaxis()->SetLabelOffset(-0.015);
    dummyGammaIndMeas->GetXaxis()->SetRangeUser(incSpectraX[0],incSpectraX[1]);
    dummyGammaIndMeas->DrawCopy();

    DrawGammaSetMarkerTGraphAsym(graphPCMIncGammaSysErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPCMIncGammaStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
    graphPCMIncGammaSysErr0010Plot->Draw("E2same");
    graphPCMIncGammaStatErr0010Plot->Draw("p,E1Z,same");

    DrawGammaSetMarkerTGraphAsym(graphPCMIncGammaSysErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPCMIncGammaStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
    graphPCMIncGammaSysErr2040Plot->Draw("E2same");
    graphPCMIncGammaStatErr2040Plot->Draw("p,E1Z,same");

    DrawGammaSetMarkerTGraphAsym(graphPCMIncGammaSysErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPCMIncGammaStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
    graphPCMIncGammaSysErr2050Plot->Draw("E2same");
    graphPCMIncGammaStatErr2050Plot->Draw("p,E1Z,same");

    TLatex *labelScalingIncGamma0010 = new TLatex(11.25,1.8E-3,"x 10^{2}");
    SetStyleTLatex( labelScalingIncGamma0010, 0.85*textsizeLabelsIncGamma,4,colorComb0010,42,kFALSE);
    labelScalingIncGamma0010->Draw();

    TLatex *labelScalingIncGamma2040 = new TLatex(11.25,1E-4,"x 10^{1}");
    SetStyleTLatex( labelScalingIncGamma2040, 0.85*textsizeLabelsIncGamma,4,colorComb2040,42,kFALSE);
    labelScalingIncGamma2040->Draw();

    TLatex *labelIncGammaThesis = new TLatex(0.6,0.93,"This thesis");
    SetStyleTLatex( labelIncGammaThesis, 0.85*textsizeLabelsIncGamma,4);
//     labelIncGammaThesis->Draw();

    TLegend* legendIncGammaPCMOnly = new TLegend(0.6,0.95-1.*0.85*textsizeLabelsIncGamma*4,0.9,0.95);//0.6,0.92-1.*0.85*textsizeLabelsIncGamma*4,0.9,0.92);
    legendIncGammaPCMOnly->SetFillStyle(0);
    legendIncGammaPCMOnly->SetFillColor(0);
    legendIncGammaPCMOnly->SetLineColor(0);
    legendIncGammaPCMOnly->SetTextSize(0.85*textsizeLabelsIncGamma);
    legendIncGammaPCMOnly->SetMargin(0.2);
    legendIncGammaPCMOnly->SetTextFont(42);
    legendIncGammaPCMOnly->SetHeader(collisionSystem.Data());
    legendIncGammaPCMOnly->AddEntry(graphPCMIncGammaSysErr0010Plot,"  0-10%","pf");
    legendIncGammaPCMOnly->AddEntry(graphPCMIncGammaSysErr2040Plot,cent2040.Data(),"pf");
    legendIncGammaPCMOnly->AddEntry(graphPCMIncGammaSysErr2050Plot,cent2050.Data(),"pf");
    legendIncGammaPCMOnly->Draw();

    dummyGammaIndMeas->Draw("axis,same");
    canvasIncGammaIndMeas->Print(Form("%s/IncGammaSpectrumPCMOnly.%s",outputDir.Data(),suffix.Data()));


    cout << "Plotting inc gamma ratio PCM only " << __LINE__ << endl;
    // **************************************************************************
    // ****************** plotting inclusive ratio with pi0 fitted **************
    // **************************************************************************
    TCanvas *canvasIncRatio = GetAndSetCanvas("canvasIncRatioFinal");
    TH2D *dummyIncR = new TH2D("dummyIncR", "dummyIncR", 120, 0., 15.5, 1000., 0., 1.6);
    SetStyleHistoTH2ForGraphs( dummyIncR, "#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}",0.045, 0.05, 0.045, 0.05, 0.85, 0.85);
    dummyIncR->GetXaxis()->SetRangeUser(0., 15.5);
    dummyIncR->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0SysErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0StatErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        graphPCMIncRFitPi0SysErr0010->Draw("E2same");
        graphPCMIncRFitPi0StatErr0010->Draw("p,E1Z,same");

        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0SysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0StatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphPCMIncRFitPi0SysErr2040->Draw("E2same");
        graphPCMIncRFitPi0StatErr2040->Draw("p,E1Z,same");

        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0SysErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0StatErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        graphPCMIncRFitPi0SysErr2050->Draw("E2same");
        graphPCMIncRFitPi0StatErr2050->Draw("p,E1Z,same");

        TLatex *labelIncGammaRatioThesis = new TLatex(0.7,0.9,"This thesis");
        SetStyleTLatex( labelIncGammaRatioThesis, 0.9*textsizeLabelsIncGamma,4);
        labelIncGammaRatioThesis->Draw();

        TLegend* legendIncGammaRatioPCMOnly = new TLegend(0.7,0.89-1.*0.9*textsizeLabelsIncGamma*4,0.9,0.89);
        legendIncGammaRatioPCMOnly->SetFillStyle(0);
        legendIncGammaRatioPCMOnly->SetFillColor(0);
        legendIncGammaRatioPCMOnly->SetLineColor(0);
        legendIncGammaRatioPCMOnly->SetTextSize(0.9*textsizeLabelsIncGamma);
        legendIncGammaRatioPCMOnly->SetMargin(0.2);
        legendIncGammaRatioPCMOnly->SetTextFont(42);
        legendIncGammaRatioPCMOnly->SetHeader(collisionSystem.Data());
        legendIncGammaRatioPCMOnly->AddEntry(graphPCMIncRFitPi0SysErr0010,cent0010.Data(),"pf");
        legendIncGammaRatioPCMOnly->AddEntry(graphPCMIncRFitPi0SysErr2040,cent2040.Data(),"pf");
        legendIncGammaRatioPCMOnly->AddEntry(graphPCMIncRFitPi0SysErr2050,cent2050.Data(),"pf");
        legendIncGammaRatioPCMOnly->Draw();

    canvasIncRatio->Print(Form("%s/IncGammaRatioPCMOnly_Pi0Fit_AllCent.%s",outputDir.Data(),suffix.Data()));

    canvasIncRatio->cd();
    dummyIncR->DrawCopy();

//         SetStyleHisto(histoIncGammaCocktail0010, 1, styleFit, kBlue+2);
//         histoIncGammaCocktail0010->Draw("l,x0,same");

        graphPCMIncRFitPi0SysErr0010->Draw("E2same");
        graphPCMIncRFitPi0StatErr0010->Draw("p,E1Z,same");

//         labelThesisHighLeft->Draw();
        TLegend* legendIncRatioFit_0010 = GetAndSetLegend(0.65,0.785,3.2);
        legendIncRatioFit_0010->SetMargin(0.2);
        legendIncRatioFit_0010->SetHeader(collisionSystem.Data());
        legendIncRatioFit_0010->AddEntry(graphPCMIncRFitPi0SysErr0010,Form("%s (#pi^{0} fitted)",cent0010.Data()),"pf");
//         legendIncRatioFit_0010->AddEntry(histoIncGammaCocktail0010,"cocktail","lp");
        legendIncRatioFit_0010->Draw();

    canvasIncRatio->Print(Form("%s/IncGammaRatioPCMOnly_Pi0Fit_0010.%s",outputDir.Data(),suffix.Data()));

    canvasIncRatio->cd();
    dummyIncR->DrawCopy();

//         SetStyleHisto(histoIncGammaCocktail2040, 1, styleFit, kBlue+2);
//         histoIncGammaCocktail2040->Draw("l,same");

        graphPCMIncRFitPi0SysErr2040->Draw("E2same");
        graphPCMIncRFitPi0StatErr2040->Draw("p,E1Z,same");
        graphPCMIncRFitPi0SysErr2050->Draw("E2same");
        graphPCMIncRFitPi0StatErr2050->Draw("p,E1Z,same");

//         labelThesisHighLeft->Draw();
        TLegend* legendIncRatioFit_SC = GetAndSetLegend(0.65,0.742,4.2);
        legendIncRatioFit_SC->SetMargin(0.2);
        legendIncRatioFit_SC->SetHeader(collisionSystem.Data());
        legendIncRatioFit_SC->AddEntry(graphPCMIncRFitPi0SysErr2040,Form("%s (#pi^{0} fitted)",cent2040.Data()),"pf");
        legendIncRatioFit_SC->AddEntry(graphPCMIncRFitPi0SysErr2050,Form("%s (#pi^{0} fitted)",cent2050.Data()),"pf");
//         legendIncRatioFit_SC->AddEntry(histoIncGammaCocktail2040,"cocktail","pl");
        legendIncRatioFit_SC->Draw();

    canvasIncRatio->Print(Form("%s/IncGammaRatioPCMOnly_Pi0Fit_SC.%s",outputDir.Data(),suffix.Data()));

    canvasIncRatio->cd();
    dummyIncR->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphPCMIncRSysErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRStatErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        graphPCMIncRSysErr0010->Draw("E2same");
        graphPCMIncRStatErr0010->Draw("p,E1Z,same");

        DrawGammaSetMarkerTGraphAsym(graphPCMIncRSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRStatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphPCMIncRSysErr2040->Draw("E2same");
        graphPCMIncRStatErr2040->Draw("p,E1Z,same");

        DrawGammaSetMarkerTGraphAsym(graphPCMIncRSysErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRStatErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        graphPCMIncRSysErr2050->Draw("E2same");
        graphPCMIncRStatErr2050->Draw("p,E1Z,same");

        labelIncGammaRatioThesis->Draw();
        legendIncGammaRatioPCMOnly->Draw();

    canvasIncRatio->Print(Form("%s/IncGammaRatioPCMOnly_AllCent.%s",outputDir.Data(),suffix.Data()));

    canvasIncRatio->cd();
    dummyIncR->DrawCopy();

        graphPCMIncRSysErr0010->Draw("E2same");
        graphPCMIncRStatErr0010->Draw("p,E1Z,same");

        labelThesisHighLeft->Draw();
        TLegend* legendIncRatio_0010 = GetAndSetLegend(0.15,0.15,2.2);
        legendIncRatio_0010->SetMargin(0.2);
        legendIncRatio_0010->SetHeader(collisionSystem.Data());
        legendIncRatio_0010->AddEntry(graphPCMIncRSysErr0010,cent0010.Data(),"pf");
        legendIncRatio_0010->Draw();

    canvasIncRatio->Print(Form("%s/IncGammaRatioPCMOnly_0010.%s",outputDir.Data(),suffix.Data()));

    canvasIncRatio->cd();
    dummyIncR->DrawCopy();

        graphPCMIncRSysErr2040->Draw("E2same");
        graphPCMIncRStatErr2040->Draw("p,E1Z,same");

        labelThesisHighLeft->Draw();
        TLegend* legendIncRatio_2040 = GetAndSetLegend(0.15,0.15,2.2);
        legendIncRatio_2040->SetMargin(0.2);
        legendIncRatio_2040->SetHeader(collisionSystem.Data());
        legendIncRatio_2040->AddEntry(graphPCMIncRSysErr2040,cent2040.Data(),"pf");
        legendIncRatio_2040->Draw();

    canvasIncRatio->Print(Form("%s/IncGammaRatioPCMOnly_2040.%s",outputDir.Data(),suffix.Data()));

    canvasIncRatio->cd();
    dummyIncR->DrawCopy();

        graphPCMIncRSysErr2050->Draw("E2same");
        graphPCMIncRStatErr2050->Draw("p,E1Z,same");

        labelThesisHighLeft->Draw();
        TLegend* legendIncRatioFit_2050 = GetAndSetLegend(0.15,0.15,2.2);
        legendIncRatioFit_2050->SetMargin(0.2);
        legendIncRatioFit_2050->SetHeader(collisionSystem.Data());
        legendIncRatioFit_2050->AddEntry(graphPCMIncRSysErr2050,cent2050.Data(),"pf");
        legendIncRatioFit_2050->Draw();

    canvasIncRatio->Print(Form("%s/IncGammaRatioPCMOnly_2050.%s",outputDir.Data(),suffix.Data()));

    canvasIncRatio->cd();
    dummyIncR->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0SysErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0StatErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        graphPCMIncRFitPi0SysErr0010->Draw("E2same");
        graphPCMIncRFitPi0StatErr0010->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRSysErr0010, markerStyleComb0010+4, markerSizeComb0010, colorCombNotFit0010 , colorCombNotFit0010,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRStatErr0010, markerStyleComb0010+4, markerSizeComb0010, colorCombNotFit0010 , colorCombNotFit0010);
        graphPCMIncRSysErr0010->Draw("E2same");
        graphPCMIncRStatErr0010->Draw("p,E1Z,same");

        labelThesisHighLeft->Draw();
        TLegend* legendIncRatioCompFit_0010 = GetAndSetLegend(0.15,0.75,3.2);
        legendIncRatioCompFit_0010->SetMargin(0.2);
        legendIncRatioCompFit_0010->SetHeader(collisionSystem.Data());
        legendIncRatioCompFit_0010->AddEntry(graphPCMIncRSysErr0010,cent0010.Data(),"pf");
        legendIncRatioCompFit_0010->AddEntry(graphPCMIncRFitPi0SysErr0010,Form("%s (#pi^{0} fitted)",cent0010.Data()),"pf");
        legendIncRatioCompFit_0010->Draw();

    canvasIncRatio->Print(Form("%s/IncGammaRatioPCMOnly_compPi0MeasAndFit_0010.%s",outputDir.Data(),suffix.Data()));

    canvasIncRatio->cd();
    dummyIncR->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0SysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0StatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphPCMIncRFitPi0SysErr2040->Draw("E2same");
        graphPCMIncRFitPi0StatErr2040->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRSysErr2040, markerStyleComb2040-6, markerSizeComb2040, colorCombNotFit2040 , colorCombNotFit2040, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRStatErr2040, markerStyleComb2040-6, markerSizeComb2040, colorCombNotFit2040 , colorCombNotFit2040);
        graphPCMIncRSysErr2040->Draw("E2same");
        graphPCMIncRStatErr2040->Draw("p,E1Z,same");

        labelThesisHighLeft->Draw();
        TLegend* legendIncRatioCompFit_2040 = GetAndSetLegend(0.15,0.75,3.2);
        legendIncRatioCompFit_2040->SetMargin(0.2);
        legendIncRatioCompFit_2040->SetHeader(collisionSystem.Data());
        legendIncRatioCompFit_2040->AddEntry(graphPCMIncRSysErr2040,cent2040.Data(),"pf");
        legendIncRatioCompFit_2040->AddEntry(graphPCMIncRFitPi0SysErr2040,Form("%s (#pi^{0} fitted)",cent2040.Data()),"pf");
        legendIncRatioCompFit_2040->Draw();

    canvasIncRatio->Print(Form("%s/IncGammaRatioPCMOnly_compPi0MeasAndFit_2040.%s",outputDir.Data(),suffix.Data()));

    canvasIncRatio->cd();
    dummyIncR->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0SysErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRFitPi0StatErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        graphPCMIncRFitPi0SysErr2050->Draw("E2same");
        graphPCMIncRFitPi0StatErr2050->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRSysErr2050, markerStyleComb2050-6, markerSizeComb2050, colorCombNotFit2050 , colorCombNotFit2050, widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphPCMIncRStatErr2050, markerStyleComb2050-6, markerSizeComb2050, colorCombNotFit2050 , colorCombNotFit2050);
        graphPCMIncRSysErr2050->Draw("E2same");
        graphPCMIncRStatErr2050->Draw("p,E1Z,same");

        labelThesisHighLeft->Draw();
        TLegend* legendIncRatioCompFit_2050 = GetAndSetLegend(0.15,0.75,3.2);
        legendIncRatioCompFit_2050->SetMargin(0.2);
        legendIncRatioCompFit_2050->SetHeader(collisionSystem.Data());
        legendIncRatioCompFit_2050->AddEntry(graphPCMIncRSysErr2050,cent2050.Data(),"pf");
        legendIncRatioCompFit_2050->AddEntry(graphPCMIncRFitPi0SysErr2050,Form("%s (#pi^{0} fitted)",cent2050.Data()),"pf");
        legendIncRatioCompFit_2050->Draw();

    canvasIncRatio->Print(Form("%s/IncGammaRatioPCMOnly_compPi0MeasAndFit_2050.%s",outputDir.Data(),suffix.Data()));


    cout << "Plotting direct gamma PCM only " << __LINE__ << endl;
    //*******************************************************************************************************************************************
    //*************************************************** Plotting direct Gamma Spectrum PCM only ***********************************************
    //*******************************************************************************************************************************************
    TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr0010Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr0010Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr0010ArPlot = NULL;
    if (graphPCMDirGammaStatErr0010) graphPCMDirGammaSpectrumStatErr0010Plot         = ScaleGraph(graphPCMDirGammaStatErr0010,scaleFactor1);
    if (graphPCMDirGammaSysErr0010) graphPCMDirGammaSpectrumSystErr0010Plot         = ScaleGraph(graphPCMDirGammaSysErr0010,scaleFactor1);
    if (graphPCMDirGammaSumErrAr0010) graphPCMDirGammaSpectrumSumErr0010ArPlot     = ScaleGraph(graphPCMDirGammaSumErrAr0010,scaleFactor1);
    TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr2040Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr2040Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr2040ArPlot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr2040ULPlot = NULL;
    if (graphPCMDirGammaStatErr2040) graphPCMDirGammaSpectrumStatErr2040Plot         = ScaleGraph(graphPCMDirGammaStatErr2040,scaleFactor2);
    if (graphPCMDirGammaSysErr2040) graphPCMDirGammaSpectrumSystErr2040Plot         = ScaleGraph(graphPCMDirGammaSysErr2040,scaleFactor2);
    if (graphPCMDirGammaSumErrAr2040) graphPCMDirGammaSpectrumSumErr2040ArPlot     = ScaleGraph(graphPCMDirGammaSumErrAr2040,scaleFactor2);
    if (graphPCMDirGammaSumErrUL2040) graphPCMDirGammaSpectrumSumErr2040ULPlot     = ScaleGraph(graphPCMDirGammaSumErrUL2040,scaleFactor2);
    TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr2050Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr2050Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr2050ArPlot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr2050UPPlot = NULL;
    if (graphPCMDirGammaStatErr2050) graphPCMDirGammaSpectrumStatErr2050Plot         = ScaleGraph(graphPCMDirGammaStatErr2050,1);
    if (graphPCMDirGammaSysErr2050) graphPCMDirGammaSpectrumSystErr2050Plot         = ScaleGraph(graphPCMDirGammaSysErr2050,1);
    if (graphPCMDirGammaSumErrAr2050) graphPCMDirGammaSpectrumSumErr2050ArPlot     = ScaleGraph(graphPCMDirGammaSumErrAr2050,1);
    if (graphPCMDirGammaSumErrUL2050) graphPCMDirGammaSpectrumSumErr2050UPPlot     = ScaleGraph(graphPCMDirGammaSumErrUL2050,1);


    TCanvas *canvasDirGammaIndMeas = new TCanvas("canvasDirGammaIndMeas","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasDirGammaIndMeas, 0.165, 0.01, 0.01, 0.07);
    canvasDirGammaIndMeas->SetLogy();
    canvasDirGammaIndMeas->SetLogx();

    Int_t textSizeLabelsPixelDirGam = 48;
    Double_t textsizeLabelsDirGamma = 0;
    if (canvasDirGammaIndMeas->XtoPixel(canvasDirGammaIndMeas->GetX2()) < canvasDirGammaIndMeas->YtoPixel(canvasDirGammaIndMeas->GetY1())){
        textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGammaIndMeas->XtoPixel(canvasDirGammaIndMeas->GetX2()) ;
    } else {
        textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGammaIndMeas->YtoPixel(canvasDirGammaIndMeas->GetY1());
    }

        TH2D *dummyDirGammaIndMeas = new TH2D("dummyDirGammaIndMeas", "dummyDirGammaIndMeas", 1000, 0., 30, 1000., 1.e-7,1e7);
        SetStyleHistoTH2ForGraphs( dummyDirGammaIndMeas, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
        dummyDirGammaIndMeas->GetXaxis()->SetLabelOffset(-0.015);
        dummyDirGammaIndMeas->GetXaxis()->SetTickLength(0.025);
        dummyDirGammaIndMeas->GetYaxis()->SetTickLength(0.025);
        dummyDirGammaIndMeas->GetXaxis()->SetRangeUser(0.5, 30.);//doubleRatioX[0],doubleRatioX[1]);
        dummyDirGammaIndMeas->GetYaxis()->SetRangeUser(1e-7,5e4);

   canvasDirGammaIndMeas->cd();
    dummyDirGammaIndMeas->DrawCopy();

        if (graphPCMDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr0010Plot->Draw("Z2same");
        }
        if (graphPCMDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPCMDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }

        if (graphPCMDirGammaSpectrumSystErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2040Plot->Draw("Z2,same");
        }
        if (graphPCMDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPCMDirGammaSpectrumStatErr2040Plot->Draw("zpE1,same");
        }

        if (graphPCMDirGammaSpectrumSystErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphPCMDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSumErr2050ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr2050ArPlot , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphPCMDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2050ArPlot);
        }

        SetStyleTLatex( labelScalingDirGamma0010, 0.85*textsizeLabelsDirGamma,4,colorComb0010,42,kFALSE);
        SetStyleTLatex( labelScalingDirGamma2040, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
        SetStyleTLatex( labelScalingDirGamma2050, 0.85*textsizeLabelsDirGamma,4,colorComb2050,42,kFALSE);
        SetStyleTLatex( labelDirGammaColl, 0.85*textsizeLabelsDirGamma,4);
        labelScalingDirGamma0010->Draw();
        labelScalingDirGamma2040->Draw();

        TLatex *labelThesisDG = new TLatex(0.63,0.94,"This thesis");
        SetStyleTLatex( labelThesisDG, 1.*0.85*textsizeLabelsDirGamma,4);
//         labelThesisDG->Draw();
        TLegend* legendDirGammaPCMOnlyForThesis = new TLegend(0.63,0.96-1.*0.85*textsizeLabelsDirGamma*4,0.9,0.96);//0.63,0.93-1.*0.85*textsizeLabelsDirGamma*4,0.9,0.93);
        legendDirGammaPCMOnlyForThesis->SetFillStyle(0);
        legendDirGammaPCMOnlyForThesis->SetFillColor(0);
        legendDirGammaPCMOnlyForThesis->SetLineColor(0);
        legendDirGammaPCMOnlyForThesis->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaPCMOnlyForThesis->SetMargin(0.2);
        legendDirGammaPCMOnlyForThesis->SetTextFont(42);
        legendDirGammaPCMOnlyForThesis->SetHeader(collisionSystem.Data());
        legendDirGammaPCMOnlyForThesis->AddEntry(graphPCMDirGammaSpectrumSystErr0010Plot,"  0-10%","pf");
        legendDirGammaPCMOnlyForThesis->AddEntry(graphPCMDirGammaSpectrumSystErr2040Plot,"20-40%","pf");
        legendDirGammaPCMOnlyForThesis->AddEntry(graphPCMDirGammaSpectrumSystErr2050Plot,"20-50%","pf");
        legendDirGammaPCMOnlyForThesis->Draw();

    canvasDirGammaIndMeas->Print(Form("%s/DirGammaSpectrumPCMOnly.%s",outputDir.Data(),suffix.Data()));

    canvasDirGammaIndMeas->cd();
    dummyDirGammaIndMeas->GetYaxis()->SetRangeUser(1e-7,2e6);
    dummyDirGammaIndMeas->DrawCopy();

        graphTheoryEPS090010Plot->Draw("p3lsame");
        graphTheoryCT100010Plot->Draw("p3lsame");
        graphTheoryNLO0010Plot->Draw("p3lsame");
        fitTheoryPromptforNLOMcGill0010->Draw("p3lsame");
        if (graphPCMDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr0010Plot->Draw("Z2same");
        }
        if (graphPCMDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPCMDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }

        graphTheoryEPS092040Plot->Draw("p3lsame");
        graphTheoryCT102040Plot->Draw("p3lsame");
        graphTheoryNLO2040Plot->Draw("p3lsame");
        fitTheoryPromptforNLOMcGill2040->Draw("p3lsame");
        if (graphPCMDirGammaSpectrumSystErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2040Plot->Draw("Z2,same");
        }
        if (graphPCMDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPCMDirGammaSpectrumStatErr2040Plot->Draw("zpE1,same");
        }

        graphTheoryEPS092050Plot->Draw("p3lsame");
        graphTheoryCT102050Plot->Draw("p3lsame");
        graphTheoryNLO2050Plot->Draw("p3lsame");
        fitTheoryPromptforNLOMcGill2050->Draw("p3lsame");
        if (graphPCMDirGammaSpectrumSystErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphPCMDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSumErr2050ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr2050ArPlot , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphPCMDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2050ArPlot);
        }

        labelScalingDirGamma0010->Draw();
        labelScalingDirGamma2040->Draw();

        //TLatex *labelThesisDGwithNLO = new TLatex(0.205,0.265,"This thesis");
        TLatex *labelThesisDGwithNLO = new TLatex(0.2,0.94,"This thesis");
        SetStyleTLatex( labelThesisDGwithNLO, 1.*0.85*textsizeLabelsDirGamma,4);
        labelThesisDGwithNLO->Draw();

        //TLegend* legendDirGammaPCMOnlyNLO = new TLegend(0.2,0.12,0.3+0.2,0.12+1.*0.85*textsizeLabelsDirGamma*4);
        TLegend* legendDirGammaPCMOnlyNLO = new TLegend(0.2,0.93-1.*0.85*textsizeLabelsDirGamma*4,0.3+0.2,0.93);
        legendDirGammaPCMOnlyNLO->SetFillStyle(0);
        legendDirGammaPCMOnlyNLO->SetFillColor(0);
        legendDirGammaPCMOnlyNLO->SetLineColor(0);
        legendDirGammaPCMOnlyNLO->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaPCMOnlyNLO->SetMargin(0.2);
        legendDirGammaPCMOnlyNLO->SetTextFont(42);
        legendDirGammaPCMOnlyNLO->SetHeader(collisionSystem.Data());
        legendDirGammaPCMOnlyNLO->AddEntry(graphPCMDirGammaSpectrumSystErr0010Plot,"  0-10%","pf");
        legendDirGammaPCMOnlyNLO->AddEntry(graphPCMDirGammaSpectrumSystErr2040Plot,"20-40%","pf");
        legendDirGammaPCMOnlyNLO->AddEntry(graphPCMDirGammaSpectrumSystErr2050Plot,"20-50%","pf");
        legendDirGammaPCMOnlyNLO->Draw();

        TLegend* legendDirGammaNLO = new TLegend(0.6,0.96-1.*0.85*textsizeLabelsDirGamma*5.5,0.6+0.21,0.96);
        legendDirGammaNLO->SetFillStyle(0);
        legendDirGammaNLO->SetFillColor(0);
        legendDirGammaNLO->SetLineColor(0);
        legendDirGammaNLO->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaNLO->SetMargin(0.25);
        legendDirGammaNLO->SetTextFont(42);
        legendDirGammaNLO->SetHeader("NLO pQCD");
        legendDirGammaNLO->AddEntry(graphTheoryNLO0010Plot,"PDF: CTEQ6M5","l"); //"PDF: CTEQ6M5, FF: GRV ","l");
        legendDirGammaNLO->AddEntry((TObject*)0,"FF: GRV","");
        legendDirGammaNLO->AddEntry(graphTheoryPromptMcGill2040Plot,"(n)PDF: ","l");
        legendDirGammaNLO->AddEntry((TObject*)0,"    CTEQ6.1M/EPS09","");
        legendDirGammaNLO->AddEntry((TObject*)0,"FF: BFG2","");
        legendDirGammaNLO->Draw();

        TLegend* legendDirGammaNLO2 = new TLegend(0.6,0.89-1.*0.85*textsizeLabelsDirGamma*7.5-0.02,0.6+0.21,0.96-1.*0.85*textsizeLabelsDirGamma*5.5-0.01);
        legendDirGammaNLO2->SetFillStyle(0);
        legendDirGammaNLO2->SetFillColor(0);
        legendDirGammaNLO2->SetLineColor(0);
        legendDirGammaNLO2->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaNLO2->SetMargin(0.25);
        legendDirGammaNLO2->SetTextFont(42);
        legendDirGammaNLO2->SetHeader("JETPHOX");
        legendDirGammaNLO2->AddEntry(graphTheoryCT100010Plot,"PDF: CT10","f");
        legendDirGammaNLO2->AddEntry(graphTheoryEPS090010Plot,"nPDF: EPS09","f");
        legendDirGammaNLO2->AddEntry((TObject*)0,"FF: BFG2","");
        legendDirGammaNLO2->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
        legendDirGammaNLO2->Draw();

    canvasDirGammaIndMeas->Print(Form("%s/DirGammaSpectrumPCMOnly_withNLO.%s",outputDir.Data(),suffix.Data()));

//     canvasDirGammaIndMeas->cd();
//     dummyDirGammaIndMeas->GetYaxis()->SetRangeUser(1e-7,1e2);
//     dummyDirGammaIndMeas->DrawCopy();
//
//         TLegend* legendDirGammaPHSD = new TLegend(0.2,0.97-1.*0.85*textsizeLabelsDirGamma*3,0.3+0.2,0.97);
//         legendDirGammaPHSD->SetFillStyle(0);
//         legendDirGammaPHSD->SetFillColor(0);
//         legendDirGammaPHSD->SetLineColor(0);
//         legendDirGammaPHSD->SetTextSize(0.85*textsizeLabelsDirGamma);
//         legendDirGammaPHSD->SetMargin(0.2);
//         legendDirGammaPHSD->SetTextFont(42);
//         legendDirGammaPHSD->SetHeader(collisionSystemCent2050.Data());
//         legendDirGammaPHSD->AddEntry(graphTheoryPHSD2050Plot,"PHSD","l");
//         legendDirGammaPHSD->AddEntry(graphTheoryPHSDplusPrompt2050Plot,"PHSD + prompt","l");
//         legendDirGammaPHSD->Draw();
//
//         graphTheoryEPS092050Plot->Draw("p3lsame");
//         graphTheoryCT102050Plot->Draw("p3lsame");
//         graphTheoryNLO2050Plot->Draw("p3lsame");
//         SetStyleGammaNLOTGraphWithBand( graphTheoryPHSD2050Plot, 3.0, stylePHSD, kRed+1, 3015, kRed+1, 0);
//         graphTheoryPHSD2050Plot->Draw("p3lsame");
//         SetStyleGammaNLOTGraphWithBand( graphTheoryPHSDplusPrompt2050Plot, 3.0, stylePHSD, kBlack, 3015, kBlack, 0);
//         graphTheoryPHSDplusPrompt2050Plot->Draw("p3lsame");
//         fitTheoryPromptMcGill2050->SetLineColor(kOrange+1);
//         fitTheoryPromptMcGill2050->Draw("l,same");
//
//         legendDirGammaNLO->Draw();
//         legendDirGammaNLO2->Draw();
//
//     canvasDirGammaIndMeas->Print(Form("%s/DirGammaPHSD.%s",outputDir.Data(),suffix.Data()));

    canvasDirGammaIndMeas->cd();
        TH2D *dummyDirGammaPCMOnly_withTheory = new TH2D("dummyDirGammaIndMeas", "dummyDirGammaIndMeas", 1000, 0., 30, 1000., 1.e-7,1e5);
        SetStyleHistoTH2ForGraphs( dummyDirGammaPCMOnly_withTheory, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
        dummyDirGammaPCMOnly_withTheory->GetXaxis()->SetLabelOffset(-0.015);
        dummyDirGammaPCMOnly_withTheory->GetXaxis()->SetTickLength(0.025);
        dummyDirGammaPCMOnly_withTheory->GetYaxis()->SetTickLength(0.025);
        dummyDirGammaPCMOnly_withTheory->GetXaxis()->SetRangeUser(0.5, 30.);//doubleRatioX[0],doubleRatioX[1]);
        dummyDirGammaPCMOnly_withTheory->DrawCopy();

        TLatex *labelScalingDirGamma0010_withTheory = new TLatex(11.5,3.5E-3,"5 x 10^{2}");
        SetStyleTLatex( labelScalingDirGamma0010_withTheory, 0.85*textsizeLabelsDirGamma,4,colorComb0010,42,kFALSE);
        TLatex *labelScalingDirGamma2040_withTheory = new TLatex(11.5,3.E-4,"5 x 10^{1}");
        SetStyleTLatex( labelScalingDirGamma2040_withTheory, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
        labelScalingDirGamma0010_withTheory->Draw();
        labelScalingDirGamma2040_withTheory->Draw();

        TLatex *labelThesisDGwithTheory = new TLatex(/*0.62,0.94*/0.205,0.265,"This thesis");
        SetStyleTLatex( labelThesisDGwithTheory, 1.*0.85*textsizeLabelsDirGamma,4);
//         labelThesisDGwithTheory->Draw();
        TLegend* legendDirGammaPCMOnlywithTheory = new TLegend(0.2,0.12,0.3+0.2,0.12+1.*0.85*textsizeLabelsDirGamma*4);//0.62,0.93-1.*0.85*textsizeLabelsDirGamma*4,0.62+0.21,0.93);
        legendDirGammaPCMOnlywithTheory->SetFillStyle(0);
        legendDirGammaPCMOnlywithTheory->SetFillColor(0);
        legendDirGammaPCMOnlywithTheory->SetLineColor(0);
        legendDirGammaPCMOnlywithTheory->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaPCMOnlywithTheory->SetMargin(0.2);
        legendDirGammaPCMOnlywithTheory->SetTextFont(42);
        legendDirGammaPCMOnlywithTheory->SetHeader(collisionSystem.Data());
        legendDirGammaPCMOnlywithTheory->AddEntry(graphPCMDirGammaSpectrumSystErr0010Plot,cent0010,"pf");
        legendDirGammaPCMOnlywithTheory->AddEntry(graphPCMDirGammaSpectrumSystErr2040Plot,cent2040,"pf");
        legendDirGammaPCMOnlywithTheory->AddEntry(graphPCMDirGammaSpectrumSystErr2050Plot,cent2050,"pf");
        legendDirGammaPCMOnlywithTheory->Draw();

//         TLegend* legendDirGammaPCMOnlyTheory = new TLegend(0.56,0.96-0.95*0.85*textsizeLabelsDirGamma*9,0.92,0.96);
//         legendDirGammaPCMOnlyTheory->SetFillStyle(0);
//         legendDirGammaPCMOnlyTheory->SetFillColor(0);
//         legendDirGammaPCMOnlyTheory->SetLineColor(0);
//         //legendDirGammaPCMOnlyTheory->SetNColumns(2);
//         legendDirGammaPCMOnlyTheory->SetTextSize(0.85*textsizeLabelsDirGamma);
//         legendDirGammaPCMOnlyTheory->SetMargin(0.15);
//         legendDirGammaPCMOnlyTheory->SetTextFont(42);
//         legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryMcGill0010Plot,"Paquet et al.","l");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"PRC 93 (2016) 044906","");
//         legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryPHSDplusPrompt0010Plot,"Linnyk et al.","l");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"PRC 92 (2015) 054914 ","");
//         legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryRapp0010Plot,"v. Hees et al.","l");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"NPA 933 (2015) 256","");
//         legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryChatterjee0010Plot,"Chatterjee et al.","l");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"PRC 85 (2012) 064910","");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"+ arXiv: 1305.0624","");
        TLegend* legendDirGammaPCMOnlyTheory = new TLegend(0.575,0.965-0.95*0.85*textsizeLabelsDirGamma*10.8,0.925,0.965);
        legendDirGammaPCMOnlyTheory->SetFillStyle(0);
        legendDirGammaPCMOnlyTheory->SetFillColor(0);
        legendDirGammaPCMOnlyTheory->SetLineColor(0);
        //legendDirGammaPCMOnlyTheory->SetNColumns(2);
        legendDirGammaPCMOnlyTheory->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaPCMOnlyTheory->SetMargin(0.15);
        legendDirGammaPCMOnlyTheory->SetTextFont(42);
        legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryMcGill0010Plot,"McGill group","l");
        legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"(0#font[122]{-}10% and 20#font[122]{-}50%)","");
        legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"in preparation","");
        legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryMcGill2040Plot,"Paquet et al. (20#font[122]{-}40%)","l");
        legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"PRC 93 (2016) 044906","");
        legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryPHSDplusPrompt0010Plot,"Linnyk et al.","l");
        legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"PRC 92 (2015) 054914 ","");
        legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryRapp0010Plot,"v. Hees et al.","l");
        legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"NPA 933 (2015) 256","");
        legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryChatterjee0010Plot,"Chatterjee et al.","l");
        legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"PRC 85 (2012) 064910","");
        legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"+ arXiv: 1305.0624","");
        legendDirGammaPCMOnlyTheory->Draw();

//         TLegend* legendDirGammaPCMOnlyTheory = new TLegend(0.19,0.35-1.*0.85*textsizeLabelsDirGamma*7.5,0.94,0.35);
//         legendDirGammaPCMOnlyTheory->SetFillStyle(0);
//         legendDirGammaPCMOnlyTheory->SetFillColor(0);
//         legendDirGammaPCMOnlyTheory->SetLineColor(0);
//         legendDirGammaPCMOnlyTheory->SetNColumns(2);
//         legendDirGammaPCMOnlyTheory->SetTextSize(0.85*textsizeLabelsDirGamma);
//         legendDirGammaPCMOnlyTheory->SetMargin(0.12);
//         legendDirGammaPCMOnlyTheory->SetTextFont(42);
//         legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryMcGill0010Plot,"Paquet et al.","l");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"","");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"PRC 93 (2016) 044906","");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"","");
//         legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryPHSDplusPrompt0010Plot,"Linnyk et al.","l");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"","");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"PRC 92 (2015) 054914 ","");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"","");
//         legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryRapp0010Plot,"v.Hees et al.","l");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"","");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"NPA 933 (2015) 256","");
//         legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryChatterjee0010Plot,"Chatterjee et al.","l");
//         legendDirGammaPCMOnlyTheory->AddEntry(graphTheoryRapp0010Plot,"Rapp et al.","l");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"PRC 85 (2012) 064910","");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"arXiv:","");
//         legendDirGammaPCMOnlyTheory->AddEntry((TObject*)0,"+ arXiv: 1305.0624","");
//         legendDirGammaPCMOnlyTheory->Draw();

        //_________________________________________________ 0-10% ____________________________________________________
        //graphTheoryChatterjee0010Plot->Draw("p3lsame");
        graphTheoryChatterjeeSummed0010Plot->Draw("p3lsame");
        graphTheoryPHSDplusPrompt0010Plot->Draw("p3lsame");
        graphTheoryRapp0010Plot->Draw("p3lsame");
        graphTheoryMcGill0010Plot->Draw("p3lsame");
        if (graphPCMDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr0010Plot->Draw("Z2same");
        }
        if (graphPCMDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPCMDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }

        //_________________________________________________ 20-40% ____________________________________________________
        //graphTheoryHe2040Plot->Draw("p3lsame"); //->same as Rapp
        graphTheoryPHSD2040Plot->Draw("p3lsame");
        graphTheoryRapp2040Plot->Draw("p3lsame");
        graphTheoryChatterjee2040Plot->Draw("p3lsame");
        graphTheoryMcGill2040Plot->Draw("p3lsame");

        if (graphPCMDirGammaSpectrumSystErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2040Plot->Draw("Z2,same");
        }
        if (graphPCMDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPCMDirGammaSpectrumStatErr2040Plot->Draw("zpE1,same");
        }

        //_________________________________________________ 20-50% ____________________________________________________
        graphTheoryChatterjeeSummed2050Plot->Draw("p3lsame");
        graphTheoryPHSDplusPrompt2050Plot->Draw("p3lsame");
        //graphTheoryChatterjee2050Plot->Draw("p3lsame");
        graphTheoryMcGill2050Plot->Draw("p3lsame");
        if (graphPCMDirGammaSpectrumSystErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphPCMDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSumErr2050ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr2050ArPlot , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphPCMDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2050ArPlot);
        }

    canvasDirGammaIndMeas->Print(Form("%s/DirGammaSpectrumPCMOnly_withModels.%s",outputDir.Data(),suffix.Data()));


    // Cut plotting graphs
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr0010Plot2 = NULL;
    if (graphPCMDirGammaSpectrumSystErr0010Plot) {
        graphPCMDirGammaSpectrumSystErr0010Plot2 = (TGraphAsymmErrors*)graphPCMDirGammaSpectrumSystErr0010Plot->Clone("graphPCMDirGammaSpectrumSystErr0010Plot2");
        while (graphPCMDirGammaSpectrumSystErr0010Plot2->GetX()[graphPCMDirGammaSpectrumSystErr0010Plot2->GetN()-1] > 6. && graphPCMDirGammaSpectrumSystErr0010Plot2->GetN()>0)
            graphPCMDirGammaSpectrumSystErr0010Plot2->RemovePoint(graphPCMDirGammaSpectrumSystErr0010Plot2->GetN()-1);
        if (graphPCMDirGammaSpectrumSystErr0010Plot2->GetN()==0) graphPCMDirGammaSpectrumSystErr0010Plot2 = NULL;
    }
    TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr0010Plot2 = NULL;
    if (graphPCMDirGammaSpectrumStatErr0010Plot) {
        graphPCMDirGammaSpectrumStatErr0010Plot2 = (TGraphAsymmErrors*)graphPCMDirGammaSpectrumStatErr0010Plot->Clone("graphPCMDirGammaSpectrumStatErr0010Plot2");
        while (graphPCMDirGammaSpectrumStatErr0010Plot2->GetX()[graphPCMDirGammaSpectrumStatErr0010Plot2->GetN()-1] > 6. && graphPCMDirGammaSpectrumStatErr0010Plot2->GetN()>0)
            graphPCMDirGammaSpectrumStatErr0010Plot2->RemovePoint(graphPCMDirGammaSpectrumStatErr0010Plot2->GetN()-1);
        if (graphPCMDirGammaSpectrumStatErr0010Plot2->GetN()==0) graphPCMDirGammaSpectrumStatErr0010Plot2 = NULL;
    }
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr2040Plot2 = NULL;
    if (graphPCMDirGammaSpectrumSystErr2040Plot) {
        graphPCMDirGammaSpectrumSystErr2040Plot2 = (TGraphAsymmErrors*)graphPCMDirGammaSpectrumSystErr2040Plot->Clone("graphPCMDirGammaSpectrumSystErr2040Plot2");
        while (graphPCMDirGammaSpectrumSystErr2040Plot2->GetX()[graphPCMDirGammaSpectrumSystErr2040Plot2->GetN()-1] > 6. && graphPCMDirGammaSpectrumSystErr2040Plot2->GetN()>0)
            graphPCMDirGammaSpectrumSystErr2040Plot2->RemovePoint(graphPCMDirGammaSpectrumSystErr2040Plot2->GetN()-1);
        if (graphPCMDirGammaSpectrumSystErr2040Plot2->GetN()==0) graphPCMDirGammaSpectrumSystErr2040Plot2 = NULL;
    }
    TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr2040Plot2 = NULL;
    if (graphPCMDirGammaSpectrumStatErr2040Plot) {
        graphPCMDirGammaSpectrumStatErr2040Plot2 = (TGraphAsymmErrors*)graphPCMDirGammaSpectrumStatErr2040Plot->Clone("graphPCMDirGammaSpectrumStatErr2040Plot2");
        while (graphPCMDirGammaSpectrumStatErr2040Plot2->GetX()[graphPCMDirGammaSpectrumStatErr2040Plot2->GetN()-1] > 6. && graphPCMDirGammaSpectrumStatErr2040Plot2->GetN()>0)
            graphPCMDirGammaSpectrumStatErr2040Plot2->RemovePoint(graphPCMDirGammaSpectrumStatErr2040Plot2->GetN()-1);
        if (graphPCMDirGammaSpectrumStatErr2040Plot2->GetN()==0) graphPCMDirGammaSpectrumStatErr2040Plot2 = NULL;
    }

    TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr2040ArPlot2 = NULL;
    if (graphPCMDirGammaSpectrumSumErr2040ArPlot) {
        graphPCMDirGammaSpectrumSumErr2040ArPlot2 = (TGraphAsymmErrors*)graphPCMDirGammaSpectrumSumErr2040ArPlot->Clone("graphPCMDirGammaSpectrumSumErr2040ArPlot2");
        while (graphPCMDirGammaSpectrumSumErr2040ArPlot2->GetX()[graphPCMDirGammaSpectrumSumErr2040ArPlot2->GetN()-1] > 6. && graphPCMDirGammaSpectrumSumErr2040ArPlot2->GetN()>0)
            graphPCMDirGammaSpectrumSumErr2040ArPlot2->RemovePoint(graphPCMDirGammaSpectrumSumErr2040ArPlot2->GetN()-1);
        if (graphPCMDirGammaSpectrumSumErr2040ArPlot2->GetN()==0) graphPCMDirGammaSpectrumSumErr2040ArPlot2 = NULL;
    }

    TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr2050Plot2 = NULL;
    if (graphPCMDirGammaSpectrumSystErr2050Plot) {
        graphPCMDirGammaSpectrumSystErr2050Plot2 = (TGraphAsymmErrors*)graphPCMDirGammaSpectrumSystErr2050Plot->Clone("graphPCMDirGammaSpectrumSystErr2050Plot2");
        while (graphPCMDirGammaSpectrumSystErr2050Plot2->GetX()[graphPCMDirGammaSpectrumSystErr2050Plot2->GetN()-1] > 6. && graphPCMDirGammaSpectrumSystErr2050Plot2->GetN()>0)
            graphPCMDirGammaSpectrumSystErr2050Plot2->RemovePoint(graphPCMDirGammaSpectrumSystErr2050Plot2->GetN()-1);
        if (graphPCMDirGammaSpectrumSystErr2050Plot2->GetN()==0) graphPCMDirGammaSpectrumSystErr2050Plot2 = NULL;
    }
    TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr2050Plot2 = NULL;
    if (graphPCMDirGammaSpectrumStatErr2050Plot) {
        graphPCMDirGammaSpectrumStatErr2050Plot2 = (TGraphAsymmErrors*)graphPCMDirGammaSpectrumStatErr2050Plot->Clone("graphPCMDirGammaSpectrumStatErr2050Plot2");
        while (graphPCMDirGammaSpectrumStatErr2050Plot2->GetX()[graphPCMDirGammaSpectrumStatErr2050Plot2->GetN()-1] > 6. && graphPCMDirGammaSpectrumStatErr2050Plot2->GetN()>0)
            graphPCMDirGammaSpectrumStatErr2050Plot2->RemovePoint(graphPCMDirGammaSpectrumStatErr2050Plot2->GetN()-1);
        if (graphPCMDirGammaSpectrumStatErr2050Plot2->GetN()==0) graphPCMDirGammaSpectrumStatErr2050Plot2 = NULL;
    }

    TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr2050ArPlot2 = NULL;
    if (graphPCMDirGammaSpectrumSumErr2050ArPlot) {
        graphPCMDirGammaSpectrumSumErr2050ArPlot2 = (TGraphAsymmErrors*)graphPCMDirGammaSpectrumSumErr2050ArPlot->Clone("graphPCMDirGammaSpectrumSumErr2050ArPlot2");
        while (graphPCMDirGammaSpectrumSumErr2050ArPlot2->GetX()[graphPCMDirGammaSpectrumSumErr2050ArPlot2->GetN()-1] > 6. && graphPCMDirGammaSpectrumSumErr2050ArPlot2->GetN()>0)
            graphPCMDirGammaSpectrumSumErr2050ArPlot2->RemovePoint(graphPCMDirGammaSpectrumSumErr2050ArPlot2->GetN()-1);
        if (graphPCMDirGammaSpectrumSumErr2050ArPlot2->GetN()==0) graphPCMDirGammaSpectrumSumErr2050ArPlot2 = NULL;
    }


    cout << "Plotting direct gamma spectra, PCM only, reduced X at " << __LINE__ << endl;
    canvasDirGammaIndMeas->cd();
    TH2D *dummDirGammaIndMeasRedX = new TH2D("dummDirGammaIndMeasRedX", "dummDirGammaIndMeasRedX", 100000, 0., 8., 1000., 1e-6,3e5);
    SetStyleHistoTH2ForGraphs( dummDirGammaIndMeasRedX, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})", 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummDirGammaIndMeasRedX->GetXaxis()->SetLabelOffset(-0.002);
    dummDirGammaIndMeasRedX->GetXaxis()->SetTickLength(0.025);
    dummDirGammaIndMeasRedX->GetYaxis()->SetTickLength(0.025);
    dummDirGammaIndMeasRedX->GetXaxis()->SetRangeUser(0.7,incSpectraX[1]);
    dummDirGammaIndMeasRedX->GetYaxis()->SetRangeUser(1e-5,2e3);
    dummDirGammaIndMeasRedX->DrawCopy();

        TLatex *labelScalingDirGamma0010_3 = new TLatex(5.5,5E-2,"5 x 10^{2}");
        SetStyleTLatex( labelScalingDirGamma0010_3, 0.85*textsizeLabelsDirGamma,4,colorComb0010,42,kFALSE);
        labelScalingDirGamma0010_3->Draw();

        TLatex *labelScalingDirGamma2040_3 = new TLatex(5.5,1E-3,"5 x 10^{1}");
        SetStyleTLatex( labelScalingDirGamma2040_3, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
        labelScalingDirGamma2040_3->Draw();

        TLatex *labelDirGammaThesis2 = new TLatex(0.205,0.265,"This thesis");
        SetStyleTLatex( labelDirGammaThesis2, 0.85*textsizeLabelsDirGamma,4);
        labelDirGammaThesis2->Draw();
        TLegend* legendDirGammaPCMOnlyRedX = new TLegend(0.2,0.12,0.3+0.2,0.12+1.*0.85*textsizeLabelsDirGamma*4);
        legendDirGammaPCMOnlyRedX->SetFillStyle(0);
        legendDirGammaPCMOnlyRedX->SetFillColor(0);
        legendDirGammaPCMOnlyRedX->SetLineColor(0);
        legendDirGammaPCMOnlyRedX->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaPCMOnlyRedX->SetMargin(0.20);
        legendDirGammaPCMOnlyRedX->SetTextFont(42);
        legendDirGammaPCMOnlyRedX->SetHeader(collisionSystem.Data());
        legendDirGammaPCMOnlyRedX->AddEntry(graphPCMDirGammaSpectrumSystErr0010Plot,"  0-10%","pf");
        legendDirGammaPCMOnlyRedX->AddEntry(graphPCMDirGammaSpectrumSystErr2040Plot,"20-40%","pf");
        legendDirGammaPCMOnlyRedX->AddEntry(graphPCMDirGammaSpectrumSystErr2050Plot,"20-50%","pf");
        legendDirGammaPCMOnlyRedX->Draw();

        if (graphPCMDirGammaSpectrumSystErr0010Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr0010Plot2->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr0010Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPCMDirGammaSpectrumStatErr0010Plot2->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSystErr2040Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2040Plot2, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2040Plot2->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2040Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPCMDirGammaSpectrumStatErr2040Plot2->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSystErr2050Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2050Plot2, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2050Plot2->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2050Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2050Plot2, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphPCMDirGammaSpectrumStatErr2050Plot2->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSumErr2050ArPlot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr2050ArPlot2 , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphPCMDirGammaSpectrumSumErr2050ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2050ArPlot2);
        }

    canvasDirGammaIndMeas->Print(Form("%s/DirGammaSpectrumPCMOnly_ReducedX.%s",outputDir.Data(),suffix.Data()));

    canvasDirGammaIndMeas->cd();
    canvasDirGammaIndMeas->SetLogx(0);
    dummDirGammaIndMeasRedX->GetXaxis()->SetRangeUser(0,6.5);
    dummDirGammaIndMeasRedX->DrawCopy();

        labelScalingDirGamma0010_3->Draw();
        labelScalingDirGamma2040_3->Draw();

        TLatex *labelDirGammaThesisWithFit = new TLatex(0.20,0.93,"ALICE work in progress");
        SetStyleTLatex( labelDirGammaThesisWithFit, 0.85*textsizeLabelsDirGamma,4);
        labelDirGammaThesisWithFit->Draw();
        TLegend* legendDirGammaPCMOnlyRedXWithFit = new TLegend(0.62,0.97-1*0.85*textsizeLabelsDirGamma*4,0.62+0.3,0.97);//0.62,0.94-1*0.85*textsizeLabelsDirGamma*4,0.62+0.3,0.94);
        legendDirGammaPCMOnlyRedXWithFit->SetFillStyle(0);
        legendDirGammaPCMOnlyRedXWithFit->SetFillColor(0);
        legendDirGammaPCMOnlyRedXWithFit->SetLineColor(0);
        legendDirGammaPCMOnlyRedXWithFit->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaPCMOnlyRedXWithFit->SetMargin(0.2);
        legendDirGammaPCMOnlyRedXWithFit->SetTextFont(42);
        legendDirGammaPCMOnlyRedXWithFit->SetHeader(collisionSystem.Data());
        legendDirGammaPCMOnlyRedXWithFit->AddEntry(graphPCMDirGammaSpectrumSystErr0010Plot,cent0010.Data(),"pf");
        legendDirGammaPCMOnlyRedXWithFit->AddEntry(graphPCMDirGammaSpectrumSystErr2040Plot,cent2040.Data(),"pf");
        legendDirGammaPCMOnlyRedXWithFit->AddEntry(graphPCMDirGammaSpectrumSystErr2050Plot,cent2050.Data(),"pf");
        legendDirGammaPCMOnlyRedXWithFit->Draw();

        if (graphPCMDirGammaSpectrumSystErr0010Plot2){
            graphPCMDirGammaSpectrumSystErr0010Plot2->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr0010Plot2){
            graphPCMDirGammaSpectrumStatErr0010Plot2->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSystErr2040Plot2){
            graphPCMDirGammaSpectrumSystErr2040Plot2->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2040Plot2){
            graphPCMDirGammaSpectrumStatErr2040Plot2->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSystErr2050Plot2){
            graphPCMDirGammaSpectrumSystErr2050Plot2->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2050Plot2){
            graphPCMDirGammaSpectrumStatErr2050Plot2->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSumErr2050ArPlot2){
            graphPCMDirGammaSpectrumSumErr2050ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2050ArPlot2);
        }

        TLegend* legendDirGammaPCMOnlyExpFit = new TLegend(0.2,0.1,0.3+0.2,0.1+1.*0.83*textsizeLabelsDirGamma*5);
        legendDirGammaPCMOnlyExpFit->SetFillStyle(0);
        legendDirGammaPCMOnlyExpFit->SetFillColor(0);
        legendDirGammaPCMOnlyExpFit->SetLineColor(0);
        legendDirGammaPCMOnlyExpFit->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaPCMOnlyExpFit->SetMargin(0.2);
        legendDirGammaPCMOnlyExpFit->SetTextFont(42);
        legendDirGammaPCMOnlyExpFit->SetHeader("#it{A} exp(-#it{p}_{T}/#it{T}_{eff})");
        legendDirGammaPCMOnlyExpFit->AddEntry(histoFitThermalGammaPCMOnly0010Stat,cent0010.Data(),"l");
        legendDirGammaPCMOnlyExpFit->AddEntry((TObject*)0,"#it{T}_{eff} = 270 #pm 30^{#it{stat}} #pm 65^{#it{syst}} MeV","");
        legendDirGammaPCMOnlyExpFit->AddEntry(histoFitThermalGammaPCMOnly2040Stat,cent2040.Data(),"l");
        legendDirGammaPCMOnlyExpFit->AddEntry((TObject*)0,"#it{T}_{eff} = 294 #pm 45^{#it{stat}} #pm 78^{#it{syst}} MeV","");
//         legendDirGammaPCMOnlyExpFit->AddEntry((TObject*)0,Form("#it{T}_{eff} = %3.f #pm %3.f^{#it{stat}} #pm %3.f^{#it{syst}} MeV",1e3*fitThermalGammaPCMOnly0010Stat->GetParameter(1),1e3*fitThermalGammaPCMOnly0010Stat->GetParError(1),1e3*fitThermalGammaPCMOnly0010Sys->GetParError(1)),"");
//         legendDirGammaPCMOnlyExpFit->AddEntry(histoFitThermalGammaPCMOnly2040Stat,cent2040.Data(),"l");
//         legendDirGammaPCMOnlyExpFit->AddEntry((TObject*)0,Form("#it{T}_{eff} = %3.f #pm %3.f^{#it{stat}} #pm %3.f^{#it{syst}} MeV",1e3*fitThermalGammaPCMOnly2040Stat->GetParameter(1),1e3*fitThermalGammaPCMOnly2040Stat->GetParError(1),1e3*fitThermalGammaPCMOnly2040Sys->GetParError(1)),"");
//         legendDirGammaPCMOnlyExpFit->AddEntry(histoFitThermalGammaPCMOnly2050Stat,cent2050.Data(),"l");
//         legendDirGammaPCMOnlyExpFit->AddEntry((TObject*)0,Form("#it{T}_{eff} = %3.3f #pm %3.3f^{#it{stat}} #pm %3.3f^{#it{syst}} GeV",fitThermalGammaPCMOnly2050Stat->GetParameter(1),fitThermalGammaPCMOnly2040Stat->GetParError(1),fitThermalGammaPCMOnly2050Sys->GetParError(1)),"");
        legendDirGammaPCMOnlyExpFit->Draw();

        histoFitThermalGammaPCMOnly0010Stat->Draw("same,l");
        histoFitThermalGammaPCMOnly2040Stat->Draw("same,l");
//         histoFitThermalGammaPCMOnly2050Stat->Draw("same,l");

    canvasDirGammaIndMeas->Print(Form("%s/DirGammaSpectrumPCMOnlyPlusFit_ReducedX.%s",outputDir.Data(),suffix.Data()));

    canvasDirGammaIndMeas->cd();
    canvasDirGammaIndMeas->SetLogx(0);
    TH2D *dummDirGammaIndMeasRedX_withTheory = new TH2D("dummDirGammaIndMeasRedX_withTheory", "dummDirGammaIndMeasRedX_withTheory", 100000, 0., 8., 1000., 2e-6,6e5);
    SetStyleHistoTH2ForGraphs( dummDirGammaIndMeasRedX_withTheory, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})", 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummDirGammaIndMeasRedX_withTheory->GetXaxis()->SetLabelOffset(-0.002);
    dummDirGammaIndMeasRedX_withTheory->GetXaxis()->SetTickLength(0.025);
    dummDirGammaIndMeasRedX_withTheory->GetYaxis()->SetTickLength(0.025);
    dummDirGammaIndMeasRedX_withTheory->GetXaxis()->SetRangeUser(0.7,incSpectraX[1]);
    //dummDirGammaIndMeasRedX_withTheory->GetYaxis()->SetRangeUser(1e-5,1.3e4);
    dummDirGammaIndMeasRedX_withTheory->GetXaxis()->SetRangeUser(0,6.5);
    dummDirGammaIndMeasRedX_withTheory->DrawCopy();

        TLegend* legendDirGammaTheoryPlusMcGillRedX = new TLegend(0.19,0.295-1.*0.85*textsizeLabelsDirGamma*6,0.94,0.295);
        legendDirGammaTheoryPlusMcGillRedX->SetFillStyle(0);
        legendDirGammaTheoryPlusMcGillRedX->SetFillColor(0);
        legendDirGammaTheoryPlusMcGillRedX->SetLineColor(0);
        legendDirGammaTheoryPlusMcGillRedX->SetNColumns(2);
        legendDirGammaTheoryPlusMcGillRedX->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaTheoryPlusMcGillRedX->SetMargin(0.15);
        legendDirGammaTheoryPlusMcGillRedX->SetTextFont(42);
        legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryMcGill0010Plot,"Paquet et al.","l");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"PRC 93 (2016) 044906","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryPHSDplusPrompt0010Plot,"Linnyk et al.","l");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"PRC 92 (2015) 054914 ","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryChatterjee0010Plot,"Chatterjee et al.","l");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryRapp0010Plot,"v.Hees et al.","l");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"PRC 85 (2012) 064910","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"NPA 933 (2015) 256","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"+ arXiv: 1305.0624","");
        //legendDirGammaTheoryPlusMcGillRedX->Draw();
        legendDirGammaPCMOnlyTheory->Draw();

        labelScalingDirGamma2040_3->Draw();
        labelScalingDirGamma0010_3->Draw();

        graphTheoryChatterjeeSummed0010Plot->Draw("p3lsame");
        graphTheoryPHSDplusPrompt0010Plot->Draw("p3lsame");
        graphTheoryRapp0010Plot->Draw("p3lsame");
        //graphTheoryChatterjee0010Plot->Draw("p3lsame");
        graphTheoryMcGill0010Plot->Draw("p3lsame");
        if (graphPCMDirGammaSpectrumSystErr0010Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr0010Plot2->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr0010Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPCMDirGammaSpectrumStatErr0010Plot2->Draw("p,E1Z,same");
        }

        graphTheoryChatterjee2040Plot->Draw("p3lsame");
        //graphTheoryHees2040Plot->Draw("p3lsame");
        //graphTheoryHe2040Plot->Draw("p3lsame");
        graphTheoryPHSD2040Plot->Draw("p3lsame");
        graphTheoryRapp2040Plot->Draw("p3lsame");
        graphTheoryMcGill2040Plot->Draw("p3lsame");

        if (graphPCMDirGammaSpectrumSystErr2040Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2040Plot2, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2040Plot2->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2040Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPCMDirGammaSpectrumStatErr2040Plot2->Draw("p,E1Z,same");
        }

        graphTheoryChatterjeeSummed2050Plot->Draw("p3lsame");
        graphTheoryPHSDplusPrompt2050Plot->Draw("p3lsame");
        //graphTheoryChatterjee2050Plot->Draw("p3lsame");
        graphTheoryMcGill2050Plot->Draw("p3lsame");
        if (graphPCMDirGammaSpectrumSystErr2050Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2050Plot2, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2050Plot2->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2050Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2050Plot2, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphPCMDirGammaSpectrumStatErr2050Plot2->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSumErr2050ArPlot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr2050ArPlot2 , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphPCMDirGammaSpectrumSumErr2050ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2050ArPlot2);
        }

        labelDirGammaThesisWithFit->Draw();
        legendDirGammaPCMOnlywithTheory->Draw();

    canvasDirGammaIndMeas->Print(Form("%s/DirGammaSpectrumPCMOnly_withModels_ReducedX.%s",outputDir.Data(),suffix.Data()));

   canvasDirGammaIndMeas->cd();
    canvasDirGammaIndMeas->SetLogx(0);
    TH2D *dummyDirGammaPCMOnly_linX = new TH2D("dummyDirGammaPCMOnly_linX", "dummyDirGammaPCMOnly_linX", 1000, 0., 30, 1000., 2.e-8,1.5e4);
    SetStyleHistoTH2ForGraphs( dummyDirGammaPCMOnly_linX, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummyDirGammaPCMOnly_linX->GetXaxis()->SetLabelOffset(-0.002);
    dummyDirGammaPCMOnly_linX->GetXaxis()->SetTickLength(0.025);
    dummyDirGammaPCMOnly_linX->GetYaxis()->SetTickLength(0.025);
    dummyDirGammaPCMOnly_linX->GetXaxis()->SetRangeUser(0,15);
    dummyDirGammaPCMOnly_linX->DrawCopy();

        if (graphPCMDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr0010Plot->Draw("Z2same");
        }
        if (graphPCMDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPCMDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSystErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2040Plot->Draw("Z2,same");
        }
        if (graphPCMDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPCMDirGammaSpectrumStatErr2040Plot->Draw("zpE1,same");
        }
        if (graphPCMDirGammaSpectrumSystErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphPCMDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSumErr2050ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr2050ArPlot , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphPCMDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2050ArPlot);
        }

        labelScalingDirGamma0010->Draw();
        labelScalingDirGamma2040->Draw();
        legendDirGammaPCMOnlyForThesis->Draw();
        labelThesisDG->Draw();

    canvasDirGammaIndMeas->Print(Form("%s/DirGammaSpectrumPCMOnly_linX.%s",outputDir.Data(),suffix.Data()));

    canvasDirGammaIndMeas->cd();
    canvasDirGammaIndMeas->SetLogx(0);
    dummyDirGammaPCMOnly_linX->GetYaxis()->SetRangeUser(2e-8,8e3);
    dummyDirGammaPCMOnly_linX->DrawCopy();

        graphTheoryEPS090010Plot->Draw("p3lsame");
        graphTheoryCT100010Plot->Draw("p3lsame");
        graphTheoryNLO0010Plot->Draw("p3lsame");
        fitTheoryPromptforNLOMcGill0010->Draw("p3lsame");
        if (graphPCMDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr0010Plot->Draw("Z2same");
        }
        if (graphPCMDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPCMDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }

        graphTheoryEPS092040Plot->Draw("p3lsame");
        graphTheoryCT102040Plot->Draw("p3lsame");
        graphTheoryNLO2040Plot->Draw("p3lsame");
        fitTheoryPromptforNLOMcGill2040->Draw("p3lsame");
        //graphTheoryPromptMcGill2040Plot->Draw("p3lsame");
        if (graphPCMDirGammaSpectrumSystErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2040Plot->Draw("Z2,same");
        }
        if (graphPCMDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPCMDirGammaSpectrumStatErr2040Plot->Draw("zpE1,same");
        }

        graphTheoryEPS092050Plot->Draw("p3lsame");
        graphTheoryCT102050Plot->Draw("p3lsame");
        graphTheoryNLO2050Plot->Draw("p3lsame");
        fitTheoryPromptforNLOMcGill2050->Draw("p3lsame");
        if (graphPCMDirGammaSpectrumSystErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphPCMDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSumErr2050ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr2050ArPlot , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphPCMDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2050ArPlot);
        }

        labelScalingDirGamma0010->Draw();
        labelScalingDirGamma2040->Draw();
        legendDirGammaPCMOnlyRedX->Draw();
        labelDirGammaThesis2->Draw();
        legendDirGammaNLO->Draw();
        legendDirGammaNLO2->Draw();

    canvasDirGammaIndMeas->Print(Form("%s/DirGammaSpectrumPCMOnly_linX_withNLO.%s",outputDir.Data(),suffix.Data()));

    //*******************************************************************************************************************************************
    //********************************************** DR plot with PCM measurements all cents in one plot  ***************************************
    //*******************************************************************************************************************************************
    TCanvas *canvasDoubleRatio = GetAndSetCanvas("canvasDoubleRatioFinal");
    canvasDoubleRatio->SetLogx();

    TH2D *dummyDR ;
    dummyDR = new TH2D("dummyDR", "dummyDR", 1000, 0., 16, 1000., doubleRatio[0]+0.1, doubleRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyDR, "#it{p}_{T} (GeV/#it{c})", "#it{R}_{#gamma} = (#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})",
                            0.045, 0.05, 0.045, 0.05, 0.85, 0.85);
    dummyDR->GetXaxis()->SetLabelOffset(-0.015);
    dummyDR->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDR->DrawCopy();

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        graphPCMDRPi0FitSysErr0010->Draw("E2same");
        histoPCMDRPi0FitStatErr0010->Draw("p,same,e0,X0");
        graphPCMDRPi0FitSysErr2040->Draw("E2same");
        histoPCMDRPi0FitStatErr2040->Draw("p,same,e0,X0");
        graphPCMDRPi0FitSysErr2050->Draw("E2same");
        histoPCMDRPi0FitStatErr2050->Draw("p,same,e0,X0");

        labelThesisHighLeft->Draw();

        TLegend* legendDRAllCentPCM = GetAndSetLegend(0.15,0.67,5);
        legendDRAllCentPCM->SetHeader(collisionSystem.Data());
        legendDRAllCentPCM->AddEntry(graphPCMDRPi0FitSysErr0010,"  0-10%","pf");
        legendDRAllCentPCM->AddEntry(graphPCMDRPi0FitSysErr2040,"20-40%","pf");
        legendDRAllCentPCM->AddEntry(graphPCMDRPi0FitSysErr2050,"20-50%","pf");

        legendDRAllCentPCM->Draw();

    canvasDoubleRatio->Print(Form("%s/DR_PCMMeasurementAllCentsInOne.%s", outputDir.Data(), suffix.Data()));


    // Calculate RAA for PCM only measurement
    Double_t newBinsPCM[19]                                     = {0.0, 0.4, 0.8, 1.,  1.2,
                                                                   1.4, 1.6, 1.8, 2.,  2.3,
                                                                   2.6, 3.0, 3.5, 4.,  5.0,
                                                                   6.0, 8.0, 10., 14.};
    cout << __LINE__ << endl;
    cout << endl << "Calculating RAA for PCM only" << endl;
    TGraphAsymmErrors* graphPCMRAADirGammaStat0010                 = NULL;
    TGraphAsymmErrors* graphPCMRAADirGammaSys0010                  = NULL;
    TGraphAsymmErrors* graphPCMRAADirGammaSum0010Ar                = NULL;
    if (graphPCMDirGammaStatErr0010)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0010, graphPCMDirGammaStatErr0010, &graphPCMRAADirGammaStat0010);
    if (graphPCMDirGammaSysErr0010)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0010, graphPCMDirGammaSysErr0010, &graphPCMRAADirGammaSys0010);
    if (graphPCMDirGammaSumErrAr0010)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0010, graphPCMDirGammaSumErrAr0010, &graphPCMRAADirGammaSum0010Ar, newBinsPCM, 18);
    //if (graphPCMRAADirGammaStat0010) graphPCMRAADirGammaStat0010->Print();
    //if (graphPCMRAADirGammaSys0010) graphPCMRAADirGammaSys0010->Print();
    //if (graphPCMRAADirGammaSum0010Ar) graphPCMRAADirGammaSum0010Ar->Print();
    TGraphAsymmErrors* graphPCMRAADirGammaStat2040                 = NULL;
    TGraphAsymmErrors* graphPCMRAADirGammaSys2040                  = NULL;
    TGraphAsymmErrors* graphPCMRAADirGammaSum2040Ar                = NULL;
    if (graphPCMDirGammaStatErr2040)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2040, graphPCMDirGammaStatErr2040, &graphPCMRAADirGammaStat2040);
    if (graphPCMDirGammaSysErr2040)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2040, graphPCMDirGammaSysErr2040, &graphPCMRAADirGammaSys2040);
    if (graphPCMDirGammaSumErrAr2040)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2040, graphPCMDirGammaSumErrAr2040, &graphPCMRAADirGammaSum2040Ar, newBinsPCM, 18);
    //if (graphPCMRAADirGammaStat2040) graphPCMRAADirGammaStat2040->Print();
    //if (graphPCMRAADirGammaSys2040) graphPCMRAADirGammaSys2040->Print();
    //if (graphPCMRAADirGammaSum2040Ar) graphPCMRAADirGammaSum2040Ar->Print();
    TGraphAsymmErrors* graphPCMRAADirGammaStat2050                 = NULL;
    TGraphAsymmErrors* graphPCMRAADirGammaSys2050                  = NULL;
    TGraphAsymmErrors* graphPCMRAADirGammaSum2050Ar                = NULL;
    if (graphPCMDirGammaStatErr2050)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2050, graphPCMDirGammaStatErr2050, &graphPCMRAADirGammaStat2050);
    if (graphPCMDirGammaSysErr2050)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2050, graphPCMDirGammaSysErr2050, &graphPCMRAADirGammaSys2050);
    if (graphPCMDirGammaSumErrAr2050)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2050, graphPCMDirGammaSumErrAr2050, &graphPCMRAADirGammaSum2050Ar, newBinsPCM, 18);
    //if (graphPCMRAADirGammaStat2050) graphPCMRAADirGammaStat2050->Print();
    //if (graphPCMRAADirGammaSys2050) graphPCMRAADirGammaSys2050->Print();
    //if (graphPCMRAADirGammaSum2050Ar) graphPCMRAADirGammaSum2050Ar->Print();

    cout << endl << "Calculating RAA for models only" << endl;
    TGraphErrors* graphTheoryRAADirGammaPHSD0010       = NULL;
    TGraphErrors* graphTheoryRAADirGammaRapp0010       = NULL;
    TGraphErrors* graphTheoryRAADirGammaMcGill0010     = NULL;
    TGraphErrors* graphTheoryRAADirGammaChatterjee0010 = NULL;
    if (graphTheoryPHSD0010)                CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill0010, graphTheoryPHSD0010, &graphTheoryRAADirGammaPHSD0010);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaPHSD0010, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
    if (graphTheoryRapp0010)                CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill0010, graphTheoryRapp0010, &graphTheoryRAADirGammaRapp0010);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaRapp0010, 3.0, styleHe, colorHe, 3015, colorHe, 0);
    if (graphTheoryMcGill0010)              CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill0010, graphTheoryMcGill0010, &graphTheoryRAADirGammaMcGill0010);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaMcGill0010, 3.0, styleNLOMcGill, colorNLOMcGill, 3015, colorNLOMcGill, 0);
    if (graphTheoryChatterjeeSummed0010)    CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill0010, graphTheoryChatterjeeSummed0010, &graphTheoryRAADirGammaChatterjee0010);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaChatterjee0010, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
    TGraphErrors* graphTheoryRAADirGammaPHSD2040       = NULL;
    TGraphErrors* graphTheoryRAADirGammaRapp2040       = NULL;
    TGraphErrors* graphTheoryRAADirGammaMcGill2040     = NULL;
    TGraphErrors* graphTheoryRAADirGammaChatterjee2040 = NULL;
    if (graphTheoryPHSD2040)                CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill2040, graphTheoryPHSD2040, &graphTheoryRAADirGammaPHSD2040);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaPHSD2040, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
    if (graphTheoryRapp2040)                CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill2040, graphTheoryRapp2040, &graphTheoryRAADirGammaRapp2040);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaRapp2040, 3.0,  styleHe, colorHe, 3015, colorHe, 0);
    if (graphTheoryMcGill2040)              CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill2040, graphTheoryMcGill2040, &graphTheoryRAADirGammaMcGill2040);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaMcGill2040, 3.0, styleNLOMcGill, colorNLOMcGill, 3015, colorNLOMcGill, 0);
    if (graphTheoryChatterjee2040)          CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill2040, graphTheoryChatterjee2040, &graphTheoryRAADirGammaChatterjee2040);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaChatterjee2040, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
    TGraphErrors* graphTheoryRAADirGammaPHSD2050       = NULL;
    TGraphErrors* graphTheoryRAADirGammaMcGill2050     = NULL;
    TGraphErrors* graphTheoryRAADirGammaChatterjee2050 = NULL;
    if (graphTheoryPHSD2050)                CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill2050, graphTheoryPHSD2050, &graphTheoryRAADirGammaPHSD2050);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaPHSD2050, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
    if (graphTheoryMcGill2050)              CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill2050, graphTheoryMcGill2050, &graphTheoryRAADirGammaMcGill2050);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaMcGill2050, 3.0, styleNLOMcGill, colorNLOMcGill, 3015, colorNLOMcGill, 0);
    if (graphTheoryChatterjeeSummed2050)    CalcRaaTheoryWithTheoryFit( fitTheoryPromptMcGill2050, graphTheoryChatterjeeSummed2050, &graphTheoryRAADirGammaChatterjee2050);
        SetStyleGammaNLOTGraphWithBand( graphTheoryRAADirGammaChatterjee2050, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);

    // *****************************************************************************************************
    // ***************************************** Plotting PCM RAA 0-10% ****************************************
    // *****************************************************************************************************
    cout << "Plotting PCM direct gamma RAA " << __LINE__ << endl;
    TGraphAsymmErrors* graphPCMRAADirGammaSys0010Plot = NULL;
    if (graphPCMRAADirGammaSys0010){
        graphPCMRAADirGammaSys0010Plot = (TGraphAsymmErrors*)graphPCMRAADirGammaSys0010->Clone("graphPCMRAADirGammaSys0010Plot");
    }
    TGraphAsymmErrors* graphPCMRAADirGammaStat0010Plot = NULL;
    if (graphPCMRAADirGammaStat0010){
        graphPCMRAADirGammaStat0010Plot = (TGraphAsymmErrors*)graphPCMRAADirGammaStat0010->Clone("graphPCMRAADirGammaStat0010Plot");
        ProduceGraphAsymmWithoutXErrors(graphPCMRAADirGammaStat0010Plot);
    }

    TCanvas* canvasRAA_0010 = new TCanvas("canvasRAA_0010","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRAA_0010,  0.1, 0.01, 0.015, 0.1);
    canvasRAA_0010->SetLogx(1);

    Int_t textSizeLabelsPixelRAA = 48;
    Double_t textsizeLabelsRAA = 0;
    if (canvasRAA_0010->XtoPixel(canvasRAA_0010->GetX2()) < canvasRAA_0010->YtoPixel(canvasRAA_0010->GetY1())){
        textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAA_0010->XtoPixel(canvasRAA_0010->GetX2()) ;
    } else {
        textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAA_0010->YtoPixel(canvasRAA_0010->GetY1());
    }

    TH2F * histo2DRAADummy;
    histo2DRAADummy = new TH2F("histo2DRAADummy","histo2DRAADummy",1000,0.,20.,1000,0.0,11.5);
    SetStyleHistoTH2ForGraphs(histo2DRAADummy, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA,0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.95,1., 510, 510); //#frac{#frac{1
    //     histo2DRAAAll3->GetYaxis()->SetRangeUser(0.05,8.);
    histo2DRAADummy->GetYaxis()->SetLabelOffset(0.005);
    histo2DRAADummy->GetXaxis()->SetLabelOffset(-0.005);
    histo2DRAADummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
//     histo2DRAADummy->GetXaxis()->SetMoreLogLabels();
    histo2DRAADummy->DrawCopy("");

        TLatex *labelRAAEnergy0010 = new TLatex(0.48,0.92,collisionSystemCent0010.Data());
        SetStyleTLatex( labelRAAEnergy0010, 0.85*textsizeLabelsRAA,4);
        //labelRAAEnergy0010->Draw();

        TBox* boxErrorNorm0010_Single = CreateBoxConv(colorComb0010Box, 0.75, 1.-normErr0010 , 0.8, 1.+normErr0010);
        boxErrorNorm0010_Single->Draw();

        SetStyleGammaNLOTGraphWithBand( graphTheoryRelErrEPS090010, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryRelErrEPS090010->Draw("p3lsame");

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1] , 1, 1 ,1.5,kGray,7);

        if (graphPCMRAADirGammaSys0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSys0010Plot, markerStyleComb0010,markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPCMRAADirGammaSys0010Plot->Draw("E2same");
        }
        if (graphPCMRAADirGammaStat0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaStat0010Plot, markerStyleComb0010,markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPCMRAADirGammaStat0010Plot->Draw("p,same,e1Z");
        }

        TLegend* legendRAAPCMOnly0010 = new TLegend(0.48,0.91-1.25*0.85*textsizeLabelsRAA*6,0.48+0.21,0.91);
        legendRAAPCMOnly0010->SetFillStyle(0);
        legendRAAPCMOnly0010->SetFillColor(0);
        legendRAAPCMOnly0010->SetLineColor(0);
        legendRAAPCMOnly0010->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAAPCMOnly0010->SetMargin(0.3);
        legendRAAPCMOnly0010->SetTextFont(42);
        legendRAAPCMOnly0010->SetHeader(collisionSystemCent0010.Data());
        if (graphPCMRAADirGammaSys0010Plot)legendRAAPCMOnly0010->AddEntry(graphPCMRAADirGammaSys0010Plot,"#it{#gamma}_{dir} (this thesis)","fp");
        legendRAAPCMOnly0010->AddEntry((TObject*)0,"pp reference: ","");
        legendRAAPCMOnly0010->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
        legendRAAPCMOnly0010->AddEntry((TObject*)0,"FF: BFG2","");
        legendRAAPCMOnly0010->AddEntry(graphTheoryRelErrEPS090010,"Rel. error JETPHOX","fp");
        legendRAAPCMOnly0010->AddEntry((TObject*)0,"nPDF: EPS09, FF: BFG2","");
        legendRAAPCMOnly0010->Draw();

        histo2DRAADummy->Draw("axis,same");
    canvasRAA_0010->Update();
    canvasRAA_0010->Print(Form("%s/DirGammaPCMOnlyRAA_0010.%s",outputDir.Data(),suffix.Data()));


    canvasRAA_0010->cd();
    histo2DRAADummy->DrawCopy("");

        boxErrorNorm0010_Single->Draw();
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1] , 1, 1 ,1.5,kGray+1,7);

        TLegend* legendRAATheory = new TLegend(0.48,0.685-1.25*0.85*textsizeLabelsRAA*7.5,0.48+0.21,0.685);
        legendRAATheory->SetFillStyle(0);
        legendRAATheory->SetFillColor(0);
        legendRAATheory->SetLineColor(0);
        legendRAATheory->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAATheory->SetMargin(0.3);
        legendRAATheory->SetTextFont(42);
        legendRAATheory->AddEntry(graphTheoryRAADirGammaMcGill0010,"McGill group, in preparation","l");//"Paquet et al.,      PRC 93 (2016) 044906","l");
        //legendRAATheory->AddEntry((TObject*)0,"in preparation","");
        legendRAATheory->AddEntry(graphTheoryRAADirGammaPHSD0010,"Linnyk et al.","l");
        legendRAATheory->AddEntry((TObject*)0,"PRC 92 (2015) 054914","");
        legendRAATheory->AddEntry(graphTheoryRAADirGammaRapp0010,"v. Hees et al.","l");
        legendRAATheory->AddEntry((TObject*)0,"NPA 933 (2015) 256","");
        legendRAATheory->AddEntry(graphTheoryRAADirGammaChatterjee0010,"Chatterjee et al.","l");
        legendRAATheory->AddEntry((TObject*)0,"PRC 85 (2012) 064910","");
        legendRAATheory->AddEntry((TObject*)0,"+ arXiv: 1305.0624","");
        legendRAATheory->Draw();

        if (graphPCMRAADirGammaSys0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSys0010Plot, markerStyleComb0010,markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPCMRAADirGammaSys0010Plot->Draw("E2same");
        }
        if (graphPCMRAADirGammaStat0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaStat0010Plot, markerStyleComb0010,markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPCMRAADirGammaStat0010Plot->Draw("p,same,e1Z");
        }

        TLegend* legendRAAPCMOnly0010withModels = new TLegend(0.48,0.95-1.25*0.85*textsizeLabelsRAA*5,0.48+0.21,0.95);
        legendRAAPCMOnly0010withModels->SetFillStyle(0);
        legendRAAPCMOnly0010withModels->SetFillColor(0);
        legendRAAPCMOnly0010withModels->SetLineColor(0);
        legendRAAPCMOnly0010withModels->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAAPCMOnly0010withModels->SetMargin(0.3);
        legendRAAPCMOnly0010withModels->SetTextFont(42);
        legendRAAPCMOnly0010withModels->SetHeader(collisionSystemCent0010.Data());
        if (graphPCMRAADirGammaSys0010Plot)legendRAAPCMOnly0010withModels->AddEntry(graphPCMRAADirGammaSys0010Plot,"#it{#gamma}_{dir}","fp"); // (this thesis)
        legendRAAPCMOnly0010withModels->AddEntry((TObject*)0,"pp reference: ","");
        legendRAAPCMOnly0010withModels->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
        legendRAAPCMOnly0010withModels->AddEntry((TObject*)0,"FF: BFG2","");
        legendRAAPCMOnly0010withModels->Draw();

        graphTheoryRAADirGammaPHSD0010->Draw("p3lsame");
        graphTheoryRAADirGammaRapp0010->Draw("p3lsame");
        graphTheoryRAADirGammaMcGill0010->Draw("p3lsame");
        graphTheoryRAADirGammaChatterjee0010->Draw("p3lsame");

        histo2DRAADummy->Draw("axis,same");
    canvasRAA_0010->Update();
    canvasRAA_0010->Print(Form("%s/DirGammaPCMOnlyRAAwithModels_0010.%s",outputDir.Data(),suffix.Data()));


    // *****************************************************************************************************
    // ***************************************** Plotting RAA 20-40% ****************************************
    // *****************************************************************************************************
    TGraphAsymmErrors* graphPCMRAADirGammaSys2040Plot = NULL;
    if (graphPCMRAADirGammaSys2040){
        graphPCMRAADirGammaSys2040Plot = (TGraphAsymmErrors*)graphPCMRAADirGammaSys2040->Clone("graphPCMRAADirGammaSys2040Plot");
    }
    TGraphAsymmErrors* graphPCMRAADirGammaStat2040Plot = NULL;
    if (graphPCMRAADirGammaStat2040){
        graphPCMRAADirGammaStat2040Plot = (TGraphAsymmErrors*)graphPCMRAADirGammaStat2040->Clone("graphPCMRAADirGammaStat2040Plot");
        ProduceGraphAsymmWithoutXErrors(graphPCMRAADirGammaStat2040Plot);
    }

    TCanvas* canvasRAA_2040 = new TCanvas("canvasRAA_2040","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRAA_2040,  0.1, 0.01, 0.015, 0.1);
    canvasRAA_2040->SetLogx(1);

    histo2DRAADummy->DrawCopy("");

        TBox* boxErrorNorm2040_Single = CreateBoxConv(colorComb2040Box, 0.75, 1.-normErr2040 , 0.8, 1.+normErr2040);
        boxErrorNorm2040_Single->Draw();

        TLatex *labelRAAEnergy2040 = new TLatex(0.48,0.92,collisionSystemCent2040.Data());
        SetStyleTLatex( labelRAAEnergy2040, 0.85*textsizeLabelsRAA,4);
        labelRAAEnergy2040->Draw();

        SetStyleGammaNLOTGraphWithBand( graphTheoryRelErrEPS092040, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryRelErrEPS092040->Draw("p3lsame");

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1] ,1, 1 ,1.5,kGray+1,7);

        if (graphPCMRAADirGammaSys2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSys2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPCMRAADirGammaSys2040Plot->Draw("E2same");
        }
        if (graphPCMRAADirGammaStat2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaStat2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPCMRAADirGammaStat2040Plot->Draw("p,E1Z,same");
        }
        if (graphPCMRAADirGammaSum2040Ar){
//             graphPCMRAADirGammaSum2040Ar->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSum2040Ar , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
            for (Int_t i = 0; i < graphPCMRAADirGammaSum2040Ar->GetN(); i++){
                graphPCMRAADirGammaSum2040Ar->SetPointEYhigh(i,0);
            }
            graphPCMRAADirGammaSum2040Ar->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMRAADirGammaSum2040Ar);
            graphPCMRAADirGammaSum2040Ar->Print();
        }
        histo2DRAADummy->Draw("axis,same");

        TLegend* legendRAAPCMOnly2040 = new TLegend(0.48,0.91-1.25*0.85*textsizeLabelsRAA*6,0.48+0.21,0.91);
        legendRAAPCMOnly2040->SetFillStyle(0);
        legendRAAPCMOnly2040->SetFillColor(0);
        legendRAAPCMOnly2040->SetLineColor(0);
        legendRAAPCMOnly2040->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAAPCMOnly2040->SetMargin(0.3);
        legendRAAPCMOnly2040->SetTextFont(42);
        if (graphPCMRAADirGammaSys2040Plot)legendRAAPCMOnly2040->AddEntry(graphPCMRAADirGammaSys2040Plot,"#it{#gamma}_{dir} (this thesis)","fp");
        legendRAAPCMOnly2040->AddEntry((TObject*)0,"pp reference: ","");
        legendRAAPCMOnly2040->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
        legendRAAPCMOnly2040->AddEntry((TObject*)0,"FF: BFG2","");
        legendRAAPCMOnly2040->AddEntry(graphTheoryRelErrEPS092040,"Rel. error JETPHOX","fp");
        legendRAAPCMOnly2040->AddEntry((TObject*)0,"nPDF: EPS09, FF: BFG2","");
        legendRAAPCMOnly2040->Draw();

        histo2DRAADummy->Draw("axis,same");
    canvasRAA_2040->Update();
    canvasRAA_2040->Print(Form("%s/DirGammaPCMOnlyRAA_2040.%s",outputDir.Data(),suffix.Data()));

    canvasRAA_2040->cd();
    histo2DRAADummy->DrawCopy("");

        boxErrorNorm2040_Single->Draw();
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1] ,1, 1 ,1.5,kGray+1,7);

        TLegend* legendRAATheory2 = new TLegend(0.48,0.688-1.25*0.85*textsizeLabelsRAA*7.5,0.48+0.21,0.688);
        legendRAATheory2->SetFillStyle(0);
        legendRAATheory2->SetFillColor(0);
        legendRAATheory2->SetLineColor(0);
        legendRAATheory2->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAATheory2->SetMargin(0.3);
        legendRAATheory2->SetTextFont(42);
        legendRAATheory2->AddEntry(graphTheoryRAADirGammaMcGill2040,"McGill group,in preparation","l");//"Paquet et al.,      PRC 93 (2016) 044906","l");
        //legendRAATheory2->AddEntry((TObject*)0,"in preparation","");
        legendRAATheory2->AddEntry(graphTheoryRAADirGammaPHSD2040,"Linnyk et al.","l");
        legendRAATheory2->AddEntry((TObject*)0,"PRC 92 (2015) 054914","");
        legendRAATheory2->AddEntry(graphTheoryRAADirGammaRapp2040,"v. Hees et al.","l");
        legendRAATheory2->AddEntry((TObject*)0,"NPA 933 (2015) 256","");
        legendRAATheory2->AddEntry(graphTheoryRAADirGammaChatterjee2040,"Chatterjee et al.","l");
        legendRAATheory2->AddEntry((TObject*)0,"PRC 85 (2012) 064910","");
        legendRAATheory2->AddEntry((TObject*)0,"       + arXiv: 1305.0624","");
        legendRAATheory2->Draw();

        if (graphPCMRAADirGammaSys2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSys2040Plot, markerStyleComb2040,markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPCMRAADirGammaSys2040Plot->Draw("E2same");
        }
        if (graphPCMRAADirGammaStat2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaStat2040Plot, markerStyleComb2040,markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPCMRAADirGammaStat2040Plot->Draw("p,same,e1Z");
        }

        TLegend* legendRAAPCMOnly2040withModels = new TLegend(0.48,0.95-1.25*0.85*textsizeLabelsRAA*5,0.48+0.21,0.95);
        legendRAAPCMOnly2040withModels->SetFillStyle(0);
        legendRAAPCMOnly2040withModels->SetFillColor(0);
        legendRAAPCMOnly2040withModels->SetLineColor(0);
        legendRAAPCMOnly2040withModels->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAAPCMOnly2040withModels->SetMargin(0.3);
        legendRAAPCMOnly2040withModels->SetTextFont(42);
        legendRAAPCMOnly2040withModels->SetHeader(collisionSystemCent2040.Data());
        if (graphPCMRAADirGammaSys2040Plot)legendRAAPCMOnly2040withModels->AddEntry(graphPCMRAADirGammaSys2040Plot,"#it{#gamma}_{dir} (this thesis)","fp");
        legendRAAPCMOnly2040withModels->AddEntry((TObject*)0,"pp reference: ","");
        legendRAAPCMOnly2040withModels->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
        legendRAAPCMOnly2040withModels->AddEntry((TObject*)0,"FF: BFG2","");
        legendRAAPCMOnly2040withModels->Draw();

        graphTheoryRAADirGammaPHSD2040->Draw("p3lsame");
        graphTheoryRAADirGammaRapp2040->Draw("p3lsame");
        graphTheoryRAADirGammaMcGill2040->Draw("p3lsame");
        graphTheoryRAADirGammaChatterjee2040->Draw("p3lsame");

        histo2DRAADummy->Draw("axis,same");
    canvasRAA_2040->Update();
    canvasRAA_2040->Print(Form("%s/DirGammaPCMOnlyRAAwithModels_2040.%s",outputDir.Data(),suffix.Data()));

    // *****************************************************************************************************
    // ***************************************** Plotting RAA 20-50% ****************************************
    // *****************************************************************************************************
    TGraphAsymmErrors* graphPCMRAADirGammaSys2050Plot = NULL;
    if (graphPCMRAADirGammaSys2050){
        graphPCMRAADirGammaSys2050Plot = (TGraphAsymmErrors*)graphPCMRAADirGammaSys2050->Clone("graphPCMRAADirGammaSys2050Plot");
    }
    TGraphAsymmErrors* graphPCMRAADirGammaStat2050Plot = NULL;
    if (graphPCMRAADirGammaStat2050){
        graphPCMRAADirGammaStat2050Plot = (TGraphAsymmErrors*)graphPCMRAADirGammaStat2050->Clone("graphPCMRAADirGammaStat2050Plot");
        ProduceGraphAsymmWithoutXErrors(graphPCMRAADirGammaStat2050Plot);
    }

    TCanvas* canvasRAA_2050 = new TCanvas("canvasRAA_2050","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRAA_2050,  0.1, 0.01, 0.015, 0.1);
    canvasRAA_2050->SetLogx(1);

    histo2DRAADummy->GetYaxis()->SetRangeUser(-0.5,14.);
    histo2DRAADummy->DrawCopy("");

        TBox* boxErrorNorm2050_Single = CreateBoxConv(colorComb2050Box, 0.75, 1.-normErr2050 , 0.8, 1.+normErr2050);
        boxErrorNorm2050_Single->Draw();

        TLatex *labelRAAEnergy2050 = new TLatex(0.48,0.92,collisionSystemCent2050.Data());
        SetStyleTLatex( labelRAAEnergy2050, 0.85*textsizeLabelsRAA,4);
        labelRAAEnergy2050->Draw();

        SetStyleGammaNLOTGraphWithBand( graphTheoryRelErrEPS092050, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryRelErrEPS092050->Draw("p3lsame");

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1] ,1, 1 ,1.5,kGray+1,7);

        if (graphPCMRAADirGammaSys2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSys2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphPCMRAADirGammaSys2050Plot->Draw("E2same");
        }
        if (graphPCMRAADirGammaStat2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaStat2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphPCMRAADirGammaStat2050Plot->Draw("p,E1Z,same");
        }
        if (graphPCMRAADirGammaSum2050Ar){
//             graphPCMRAADirGammaSum2050Ar->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSum2050Ar , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            for (Int_t i = 0; i < graphPCMRAADirGammaSum2050Ar->GetN(); i++){
                graphPCMRAADirGammaSum2050Ar->SetPointEYhigh(i,0);
            }
            graphPCMRAADirGammaSum2050Ar->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMRAADirGammaSum2050Ar);
            graphPCMRAADirGammaSum2050Ar->Print();
        }

        TLegend* legendRAAPCMOnly2050 = new TLegend(0.48,0.91-1.25*0.85*textsizeLabelsRAA*6,0.48+0.21,0.91);
        legendRAAPCMOnly2050->SetFillStyle(0);
        legendRAAPCMOnly2050->SetFillColor(0);
        legendRAAPCMOnly2050->SetLineColor(0);
        legendRAAPCMOnly2050->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAAPCMOnly2050->SetMargin(0.3);
        legendRAAPCMOnly2050->SetTextFont(42);
        if (graphPCMRAADirGammaSys2050Plot)legendRAAPCMOnly2050->AddEntry(graphPCMRAADirGammaSys2050Plot,"#it{#gamma}_{dir} (this thesis)","fp");
        legendRAAPCMOnly2050->AddEntry((TObject*)0,"pp reference: ","");
        legendRAAPCMOnly2050->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
        legendRAAPCMOnly2050->AddEntry((TObject*)0,"FF: BFG2","");
        legendRAAPCMOnly2050->AddEntry(graphTheoryRelErrEPS092050,"Rel. error JETPHOX","fp");
        legendRAAPCMOnly2050->AddEntry((TObject*)0,"nPDF: EPS09, FF: BFG2","");
        legendRAAPCMOnly2050->Draw();


        histo2DRAADummy->Draw("axis,same");
    canvasRAA_2050->Update();
    canvasRAA_2050->Print(Form("%s/DirGammaPCMOnlyRAA_2050.%s",outputDir.Data(),suffix.Data()));

    canvasRAA_2050->cd();
    histo2DRAADummy->DrawCopy("");

        boxErrorNorm2050_Single->Draw();
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1] ,1, 1 ,1.5,kGray+1,7);

        TLegend* legendRAAPCMOnly2050withModels = new TLegend(0.48,0.95-1.25*0.85*textsizeLabelsRAA*5,0.48+0.21,0.95);
        legendRAAPCMOnly2050withModels->SetFillStyle(0);
        legendRAAPCMOnly2050withModels->SetFillColor(0);
        legendRAAPCMOnly2050withModels->SetLineColor(0);
        legendRAAPCMOnly2050withModels->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAAPCMOnly2050withModels->SetMargin(0.3);
        legendRAAPCMOnly2050withModels->SetTextFont(42);
        legendRAAPCMOnly2050withModels->SetHeader(collisionSystemCent2050.Data());
        if (graphPCMRAADirGammaSys2050Plot)legendRAAPCMOnly2050withModels->AddEntry(graphPCMRAADirGammaSys2050Plot,"#it{#gamma}_{dir} (this thesis)","fp");
        legendRAAPCMOnly2050withModels->AddEntry((TObject*)0,"pp reference: ","");
        legendRAAPCMOnly2050withModels->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
        legendRAAPCMOnly2050withModels->AddEntry((TObject*)0,"FF: BFG2","");
        legendRAAPCMOnly2050withModels->Draw();

        TLegend* legendRAATheory3 = new TLegend(0.48,0.685-1.25*0.85*textsizeLabelsRAA*5.5,0.48+0.21,0.685);
        legendRAATheory3->SetFillStyle(0);
        legendRAATheory3->SetFillColor(0);
        legendRAATheory3->SetLineColor(0);
        legendRAATheory3->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAATheory3->SetMargin(0.3);
        legendRAATheory3->SetTextFont(42);
        legendRAATheory3->AddEntry(graphTheoryRAADirGammaMcGill2050,"McGill group, in preparation","l");//"Paquet et al.,      PRC 93 (2016) 044906","l");
        //legendRAATheory3->AddEntry((TObject*)0,"in preparation","");
        legendRAATheory3->AddEntry(graphTheoryRAADirGammaPHSD2050,"Linnyk et al.","l");
        legendRAATheory3->AddEntry((TObject*)0,"PRC 92 (2015) 054914","");
        legendRAATheory3->AddEntry(graphTheoryRAADirGammaChatterjee2050,"Chatterjee et al.","l");
        legendRAATheory3->AddEntry((TObject*)0,"PRC 85 (2012) 064910","");
        legendRAATheory3->AddEntry((TObject*)0,"   + arXiv: 1305.0624","");
        legendRAATheory3->Draw();

        if (graphPCMRAADirGammaSys2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSys2050Plot, markerStyleComb2050,markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphPCMRAADirGammaSys2050Plot->Draw("E2same");
        }
        if (graphPCMRAADirGammaStat2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaStat2050Plot, markerStyleComb2050,markerSizeComb2050, colorComb2050 , colorComb2050);
            graphPCMRAADirGammaStat2050Plot->Draw("p,same,e1Z");
        }
        if (graphPCMRAADirGammaSum2050Ar){
//             graphPCMRAADirGammaSum2050Ar->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSum2050Ar , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            for (Int_t i = 0; i < graphPCMRAADirGammaSum2050Ar->GetN(); i++){
                graphPCMRAADirGammaSum2050Ar->SetPointEYhigh(i,0);
            }
            graphPCMRAADirGammaSum2050Ar->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMRAADirGammaSum2050Ar);
            graphPCMRAADirGammaSum2050Ar->Print();
        }

        graphTheoryRAADirGammaPHSD2050->Draw("p3lsame");
        graphTheoryRAADirGammaMcGill2050->Draw("p3lsame");
        graphTheoryRAADirGammaChatterjee2050->Draw("p3lsame");

        histo2DRAADummy->Draw("axis,same");
    canvasRAA_2050->Update();
    canvasRAA_2050->Print(Form("%s/DirGammaPCMOnlyRAAwithModels_2050.%s",outputDir.Data(),suffix.Data()));




    //*******************************************************************************************************************************************
    //************************************************* Plotting direct Gamma RAA with other meas ***********************************************
    //*******************************************************************************************************************************************
    cout << "Plotting PCM direct gamma RAA with Pi0 RAA for PbPb and pPb" << __LINE__ << endl;
    TGraphAsymmErrors *graphPCMRAADirGammaSys0010Plot2 = (TGraphAsymmErrors*)graphPCMRAADirGammaSys0010Plot->Clone("graphPCMRAADirGammaSys0010Plot2");
    TGraphAsymmErrors *graphPCMRAADirGammaStat0010Plot2 = (TGraphAsymmErrors*)graphPCMRAADirGammaStat0010Plot->Clone("graphPCMRAADirGammaStat0010Plot2");
    TGraphAsymmErrors *graphPCMRAADirGammaSys2040Plot2 = (TGraphAsymmErrors*)graphPCMRAADirGammaSys2040Plot->Clone("graphPCMRAADirGammaSys2040Plot2");
    TGraphAsymmErrors *graphPCMRAADirGammaStat2040Plot2 = (TGraphAsymmErrors*)graphPCMRAADirGammaStat2040Plot->Clone("graphPCMRAADirGammaStat2040Plot2");
    TGraphAsymmErrors *graphPCMRAADirGammaSys2050Plot2 = (TGraphAsymmErrors*)graphPCMRAADirGammaSys2050Plot->Clone("graphPCMRAADirGammaSys2050Plot2");
    TGraphAsymmErrors *graphPCMRAADirGammaStat2050Plot2 = (TGraphAsymmErrors*)graphPCMRAADirGammaStat2050Plot->Clone("graphPCMRAADirGammaStat2050Plot2");

    TCanvas* canvasRAADirGammaPi0_0010 = new TCanvas("canvasRAADirGammaPi0_0010","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRAADirGammaPi0_0010,  0.1, 0.01, 0.015, 0.1);
    canvasRAADirGammaPi0_0010->SetLogx(1);

    textSizeLabelsPixelRAA = 48;
    textsizeLabelsRAA = 0;
    if (canvasRAADirGammaPi0_0010->XtoPixel(canvasRAADirGammaPi0_0010->GetX2()) < canvasRAADirGammaPi0_0010->YtoPixel(canvasRAADirGammaPi0_0010->GetY1())){
        textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAADirGammaPi0_0010->XtoPixel(canvasRAADirGammaPi0_0010->GetX2()) ;
    } else {
        textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAADirGammaPi0_0010->YtoPixel(canvasRAADirGammaPi0_0010->GetY1());
    }

    TH2F * histo2DRAADummyFullRange2 = new TH2F("histo2DRAADummyFullRange2","histo2DRAADummyFullRange2",1000,0.,40.,1000,0.0,11.5);
    SetStyleHistoTH2ForGraphs(histo2DRAADummyFullRange2, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA,0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.95,1., 510, 510);
    histo2DRAADummyFullRange2->GetYaxis()->SetLabelOffset(0.005);
    histo2DRAADummyFullRange2->GetXaxis()->SetLabelOffset(-0.005);
    histo2DRAADummyFullRange2->GetXaxis()->SetRangeUser(0.8,25.);
    histo2DRAADummyFullRange2->GetYaxis()->SetRangeUser(0.,7.5);
    histo2DRAADummyFullRange2->GetXaxis()->SetMoreLogLabels();
    histo2DRAADummyFullRange2->DrawCopy("");

        DrawGammaLines(0.8, 25. , 1, 1 ,1.5,kGray+1,7);

        TLatex *labelThesisRAARight = new TLatex(0.52,0.91,"This thesis");
        SetStyleTLatex( labelThesisRAARight, 0.85*textsizeLabelsRAA,4);
//         labelThesisRAARight->Draw();

        TBox* boxErrorNorm0010_Single2 = CreateBoxConv(colorNLOMcGill, 0.8, 1.-normErr0010 , 0.85, 1.+normErr0010);
        boxErrorNorm0010_Single2->Draw();
        TBox* boxErrorNorm0010_Single3 = CreateBoxConv(colorComb0010Box, 0.85, 1.-normErr0010 , 0.9, 1.+normErr0010);
        boxErrorNorm0010_Single3->Draw();
        TBox* boxErrorNormpPb_Single = CreateBoxConv(colorPHSD, 0.9, 1.-normErrpPb , 0.95, 1.+normErrpPb);
        boxErrorNormpPb_Single->Draw();

        TLegend* legendRAADirGammaPi00010V2 = new TLegend(0.52,0.92-1.25*0.85*textsizeLabelsRAA*6.5,0.52+0.21,0.92);//0.52,0.9-1.25*0.85*textsizeLabelsRAA*6.5,0.52+0.21,0.9);
        legendRAADirGammaPi00010V2->SetFillStyle(0);
        legendRAADirGammaPi00010V2->SetFillColor(0);
        legendRAADirGammaPi00010V2->SetLineColor(0);
        legendRAADirGammaPi00010V2->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAADirGammaPi00010V2->SetMargin(0.3);
        legendRAADirGammaPi00010V2->SetTextFont(42);
        legendRAADirGammaPi00010V2->SetHeader(collisionSystemCent0010.Data());
        if (graphPCMRAADirGammaSys0010Plot)legendRAADirGammaPi00010V2->AddEntry(graphPCMRAADirGammaSys0010Plot2,"#it{#gamma}_{dir} ","fp");
        legendRAADirGammaPi00010V2->AddEntry((TObject*)0,"pp reference: ","");
        legendRAADirGammaPi00010V2->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09","");
        legendRAADirGammaPi00010V2->AddEntry((TObject*)0,"    FF: BFG2","");
        legendRAADirGammaPi00010V2->AddEntry(graphRAAPi0CombPbPb2760GeVSysErr_0010,"#pi^{0} combined","fp");
        legendRAADirGammaPi00010V2->AddEntry(graphRAAEtaCombPbPb2760GeVSysErr_0010,"#eta combined","fp");
        legendRAADirGammaPi00010V2->Draw();

        TLegend* legendRAADirGammaPi0pPb0100V2 = new TLegend(0.52,0.59-1.25*0.85*textsizeLabelsRAA*2,0.52+0.21,0.59);//0.52,0.57-1.25*0.85*textsizeLabelsRAA*2,0.52+0.21,0.57);
        legendRAADirGammaPi0pPb0100V2->SetFillStyle(0);
        legendRAADirGammaPi0pPb0100V2->SetFillColor(0);
        legendRAADirGammaPi0pPb0100V2->SetLineColor(0);
        legendRAADirGammaPi0pPb0100V2->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAADirGammaPi0pPb0100V2->SetMargin(0.3);
        legendRAADirGammaPi0pPb0100V2->SetTextFont(42);
        legendRAADirGammaPi0pPb0100V2->SetHeader(collisionSystempPb.Data());
        legendRAADirGammaPi0pPb0100V2->AddEntry(CombinedPi0RpPbSystErr,"#pi^{0} ALICE preliminary","fp");
        legendRAADirGammaPi0pPb0100V2->Draw();

        if (graphPCMRAADirGammaSys0010Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSys0010Plot2, markerStyleComb0010,markerSizeComb0010, colorNLOMcGill , colorNLOMcGill, widthLinesBoxes, kTRUE);
            graphPCMRAADirGammaSys0010Plot2->Draw("E2same");
        }
        if (graphPCMRAADirGammaStat0010Plot2){
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaStat0010Plot2, markerStyleComb0010,markerSizeComb0010, colorNLOMcGill , colorNLOMcGill);
            graphPCMRAADirGammaStat0010Plot2->Draw("p,same,e1Z");
        }

        DrawGammaSetMarkerTGraphAsym(graphRAAPi0CombPbPb2760GeVSysErr_0010, markerStyleComb0010+13,markerSizeComb0010+1, colorComb0010-7 , colorComb0010-7, widthLinesBoxes, kTRUE);
        graphRAAPi0CombPbPb2760GeVSysErr_0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRAAPi0CombPbPb2760GeVStatErr_0010, markerStyleComb0010+13,markerSizeComb0010+1, colorComb0010-7, colorComb0010-7);
        graphRAAPi0CombPbPb2760GeVStatErr_0010->Draw("p,same,e1Z");

        DrawGammaSetMarkerTGraphAsym(graphRAAEtaCombPbPb2760GeVSysErr_0010, markerStyleComb0010,markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
        graphRAAEtaCombPbPb2760GeVSysErr_0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRAAEtaCombPbPb2760GeVStatErr_0010, markerStyleComb0010,markerSizeComb0010, colorComb0010 , colorComb0010);
        graphRAAEtaCombPbPb2760GeVStatErr_0010->Draw("p,same,e1Z");

        DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbSystErr, markerStyleComb0010,markerSizeComb0010, colorPHSD , colorPHSD, widthLinesBoxes, kTRUE);
        CombinedPi0RpPbSystErr->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbStatErr, markerStyleComb0010,markerSizeComb0010, colorPHSD , colorPHSD);
        CombinedPi0RpPbStatErr->Draw("p,same,e1Z");

        histo2DRAADummyFullRange2->Draw("axis,same");
    canvasRAADirGammaPi0_0010->Update();
    canvasRAADirGammaPi0_0010->Print(Form("%s/DirGammaPi0PbPbandpPbRAAfullrange_0010.%s",outputDir.Data(),suffix.Data()));


    TH2F * histo2DRAADummy2 = new TH2F("histo2DRAADummy2","histo2DRAADummy2",1000,0.,40.,1000,0.0,11.5);
    SetStyleHistoTH2ForGraphs(histo2DRAADummy2, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA,0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.95,1., 510, 510);
    histo2DRAADummy2->GetYaxis()->SetLabelOffset(0.005);
    histo2DRAADummy2->GetXaxis()->SetLabelOffset(-0.005);
    histo2DRAADummy2->GetXaxis()->SetRangeUser(0.8,25.);
    histo2DRAADummy2->GetYaxis()->SetRangeUser(0.,4.);
//     histo2DRAADummy2->GetXaxis()->SetMoreLogLabels();
    histo2DRAADummy2->DrawCopy("");

        DrawGammaLines(0.8, 25. , 1, 1 ,1.5,kGray+1,7);

        TLatex *labelThesis = new TLatex(0.15,0.91,"This thesis");
        SetStyleTLatex( labelThesis, 0.85*textsizeLabelsRAA,4);
//         labelThesis->Draw();

//         TBox* boxErrorNorm0010_Single2 = CreateBoxConv(colorNLOMcGill, 0.8, 1.-normErr0010 , 0.85, 1.+normErr0010);
        boxErrorNorm0010_Single2->Draw();
//         TBox* boxErrorNorm0010_Single3 = CreateBoxConv(colorComb0010Box, 0.85, 1.-normErr0010 , 0.9, 1.+normErr0010);
        boxErrorNorm0010_Single3->Draw();
//         TBox* boxErrorNormpPb_Single = CreateBoxConv(colorPHSD, 0.9, 1.-normErrpPb , 0.95, 1.+normErrpPb);
        boxErrorNormpPb_Single->Draw();

        TLegend* legendRAADirGammaPi00010 = new TLegend(0.52,0.95-1.25*0.85*textsizeLabelsRAA*6.5,0.52+0.21,0.95);
        legendRAADirGammaPi00010->SetFillStyle(0);
        legendRAADirGammaPi00010->SetFillColor(0);
        legendRAADirGammaPi00010->SetLineColor(0);
        legendRAADirGammaPi00010->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAADirGammaPi00010->SetMargin(0.3);
        legendRAADirGammaPi00010->SetTextFont(42);
        legendRAADirGammaPi00010->SetHeader(collisionSystemCent0010.Data());
        if (graphPCMRAADirGammaSys0010Plot)legendRAADirGammaPi00010->AddEntry(graphPCMRAADirGammaSys0010Plot2,"#it{#gamma}_{dir} ","fp");
        legendRAADirGammaPi00010->AddEntry((TObject*)0,"pp reference: ","");
        legendRAADirGammaPi00010->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09","");
        legendRAADirGammaPi00010->AddEntry((TObject*)0,"    FF: BFG2","");
        legendRAADirGammaPi00010->AddEntry(graphRAAPi0CombPbPb2760GeVSysErr_0010,"#pi^{0} combined","fp");
        legendRAADirGammaPi00010->AddEntry(graphRAAEtaCombPbPb2760GeVSysErr_0010,"#eta combined","fp");
        legendRAADirGammaPi00010->Draw();

        TLegend* legendRAADirGammaPi0pPb0100 = new TLegend(0.52,0.62-1.25*0.85*textsizeLabelsRAA*2,0.52+0.21,0.62);
        legendRAADirGammaPi0pPb0100->SetFillStyle(0);
        legendRAADirGammaPi0pPb0100->SetFillColor(0);
        legendRAADirGammaPi0pPb0100->SetLineColor(0);
        legendRAADirGammaPi0pPb0100->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAADirGammaPi0pPb0100->SetMargin(0.3);
        legendRAADirGammaPi0pPb0100->SetTextFont(42);
        legendRAADirGammaPi0pPb0100->SetHeader(collisionSystempPb.Data());
        legendRAADirGammaPi0pPb0100->AddEntry(CombinedPi0RpPbSystErr,"#pi^{0} ALICE preliminary","fp");
        legendRAADirGammaPi0pPb0100->Draw();

        if (graphPCMRAADirGammaSys0010Plot2){
            while(graphPCMRAADirGammaSys0010Plot2->GetX()[0]<2.3) graphPCMRAADirGammaSys0010Plot2->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSys0010Plot2, markerStyleComb0010,markerSizeComb0010, colorNLOMcGill , colorNLOMcGill, widthLinesBoxes, kTRUE);
            graphPCMRAADirGammaSys0010Plot2->Draw("E2same");
        }
        if (graphPCMRAADirGammaStat0010Plot2){
            while(graphPCMRAADirGammaStat0010Plot2->GetX()[0]<2.3) graphPCMRAADirGammaStat0010Plot2->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaStat0010Plot2, markerStyleComb0010,markerSizeComb0010, colorNLOMcGill , colorNLOMcGill);
            graphPCMRAADirGammaStat0010Plot2->Draw("p,same,e1Z");
        }

        DrawGammaSetMarkerTGraphAsym(graphRAAPi0CombPbPb2760GeVSysErr_0010, markerStyleComb0010+13,markerSizeComb0010+1, colorComb0010-7 , colorComb0010-7, widthLinesBoxes, kTRUE);
        graphRAAPi0CombPbPb2760GeVSysErr_0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRAAPi0CombPbPb2760GeVStatErr_0010, markerStyleComb0010+13,markerSizeComb0010+1, colorComb0010-7, colorComb0010-7);
        graphRAAPi0CombPbPb2760GeVStatErr_0010->Draw("p,same,e1Z");

        DrawGammaSetMarkerTGraphAsym(graphRAAEtaCombPbPb2760GeVSysErr_0010, markerStyleComb0010,markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
        graphRAAEtaCombPbPb2760GeVSysErr_0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRAAEtaCombPbPb2760GeVStatErr_0010, markerStyleComb0010,markerSizeComb0010, colorComb0010 , colorComb0010);
        graphRAAEtaCombPbPb2760GeVStatErr_0010->Draw("p,same,e1Z");

        DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbSystErr, markerStyleComb0010,markerSizeComb0010, colorPHSD , colorPHSD, widthLinesBoxes, kTRUE);
        CombinedPi0RpPbSystErr->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbStatErr, markerStyleComb0010,markerSizeComb0010, colorPHSD , colorPHSD);
        CombinedPi0RpPbStatErr->Draw("p,same,e1Z");

        histo2DRAADummy2->Draw("axis,same");
    canvasRAADirGammaPi0_0010->Update();
    canvasRAADirGammaPi0_0010->Print(Form("%s/DirGammaPi0PbPbandpPbRAA_0010.%s",outputDir.Data(),suffix.Data()));

        return;


    canvasRAADirGammaPi0_0010->cd();
    histo2DRAADummy2->DrawCopy("");

        DrawGammaLines(0.8, 25. , 1, 1 ,1.5,kGray+1,7);

        labelThesis->Draw();

        TBox* boxErrorNorm2040_Single2 = CreateBoxConv(colorNLOMcGill, 0.8, 1.-normErr2040 , 0.85, 1.+normErr2040);
        boxErrorNorm2040_Single2->Draw();
        TBox* boxErrorNorm2040_Single3 = CreateBoxConv(colorComb2040Box, 0.85, 1.-normErr2040 , 0.9, 1.+normErr2040);
        boxErrorNorm2040_Single3->Draw();
        boxErrorNormpPb_Single->Draw();

        TLegend* legendRAADirGammaPi02040 = new TLegend(0.52,0.95-1.25*0.85*textsizeLabelsRAA*6.5,0.52+0.21,0.95);
        legendRAADirGammaPi02040->SetFillStyle(0);
        legendRAADirGammaPi02040->SetFillColor(0);
        legendRAADirGammaPi02040->SetLineColor(0);
        legendRAADirGammaPi02040->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAADirGammaPi02040->SetMargin(0.3);
        legendRAADirGammaPi02040->SetTextFont(42);
        legendRAADirGammaPi02040->SetHeader(collisionSystemCent2040.Data());
        if (graphPCMRAADirGammaSys2040Plot)legendRAADirGammaPi02040->AddEntry(graphPCMRAADirGammaSys2040Plot2,"#it{#gamma}_{dir} ","fp");
        legendRAADirGammaPi02040->AddEntry((TObject*)0,"pp reference: ","");
        legendRAADirGammaPi02040->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09","");
        legendRAADirGammaPi02040->AddEntry((TObject*)0,"    FF: BFG2","");
        legendRAADirGammaPi02040->AddEntry(graphRAAPi0CombPbPb2760GeVSysErr_2040,"#pi^{0} combined","fp");
        legendRAADirGammaPi02040->AddEntry(graphRAAEtaCombPbPb2760GeVSysErr_2040,"#eta combined","fp");
        legendRAADirGammaPi02040->Draw();

        legendRAADirGammaPi0pPb0100->Draw();

        if (graphPCMRAADirGammaSys2040Plot2){
            while(graphPCMRAADirGammaSys2040Plot2->GetX()[0]<2.3) graphPCMRAADirGammaSys2040Plot2->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSys2040Plot2, markerStyleComb2040,markerSizeComb2040, colorNLOMcGill , colorNLOMcGill, widthLinesBoxes, kTRUE);
            graphPCMRAADirGammaSys2040Plot2->Draw("E2same");
        }
        if (graphPCMRAADirGammaStat2040Plot2){
            while(graphPCMRAADirGammaStat2040Plot2->GetX()[0]<2.3) graphPCMRAADirGammaStat2040Plot2->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaStat2040Plot2, markerStyleComb2040,markerSizeComb2040, colorNLOMcGill , colorNLOMcGill);
            graphPCMRAADirGammaStat2040Plot2->Draw("p,same,e1Z");
        }


        DrawGammaSetMarkerTGraphAsym(graphRAAPi0CombPbPb2760GeVSysErr_2040, markerStyleComb2040,markerSizeComb2040, colorComb2040-7 , colorComb2040-7, widthLinesBoxes, kTRUE);
        graphRAAPi0CombPbPb2760GeVSysErr_2040->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRAAPi0CombPbPb2760GeVStatErr_2040, markerStyleComb2040,markerSizeComb2040, colorComb2040-7 , colorComb2040-7);
        graphRAAPi0CombPbPb2760GeVStatErr_2040->Draw("p,same,e1Z");

        DrawGammaSetMarkerTGraphAsym(graphRAAEtaCombPbPb2760GeVSysErr_2040, markerStyleComb2040,markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
        graphRAAEtaCombPbPb2760GeVSysErr_2040->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRAAEtaCombPbPb2760GeVStatErr_2040, markerStyleComb2040,markerSizeComb2040, colorComb2040 , colorComb2040);
        graphRAAEtaCombPbPb2760GeVStatErr_2040->Draw("p,same,e1Z");

        CombinedPi0RpPbSystErr->Draw("E2same");
        CombinedPi0RpPbStatErr->Draw("p,same,e1Z");

        histo2DRAADummy2->Draw("axis,same");
    canvasRAADirGammaPi0_0010->Update();
    canvasRAADirGammaPi0_0010->Print(Form("%s/DirGammaPi0PbPbandpPbRAA_2040.%s",outputDir.Data(),suffix.Data()));

    canvasRAADirGammaPi0_0010->cd();
    histo2DRAADummy2->DrawCopy("");

        DrawGammaLines(0.8, 25. , 1, 1 ,1.5,kGray+1,7);

        labelThesis->Draw();

        TBox* boxErrorNorm2050_Single2 = CreateBoxConv(colorNLOMcGill, 0.8, 1.-normErr2050 , 0.85, 1.+normErr2050);
        boxErrorNorm2050_Single2->Draw();
        TBox* boxErrorNorm2050_Single3 = CreateBoxConv(colorComb2050Box, 0.85, 1.-normErr2050 , 0.9, 1.+normErr2050);
        boxErrorNorm2050_Single3->Draw();
        boxErrorNormpPb_Single->Draw();

        TLegend* legendRAADirGammaPi02050 = new TLegend(0.52,0.95-1.25*0.85*textsizeLabelsRAA*6.5,0.52+0.21,0.95);
        legendRAADirGammaPi02050->SetFillStyle(0);
        legendRAADirGammaPi02050->SetFillColor(0);
        legendRAADirGammaPi02050->SetLineColor(0);
        legendRAADirGammaPi02050->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAADirGammaPi02050->SetMargin(0.3);
        legendRAADirGammaPi02050->SetTextFont(42);
        legendRAADirGammaPi02050->SetHeader(collisionSystemCent2050.Data());
        if (graphPCMRAADirGammaSys2050Plot)legendRAADirGammaPi02050->AddEntry(graphPCMRAADirGammaSys2050Plot2,"#it{#gamma}_{dir} ","fp");
        legendRAADirGammaPi02050->AddEntry((TObject*)0,"pp reference: ","");
        legendRAADirGammaPi02050->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09","");
        legendRAADirGammaPi02050->AddEntry((TObject*)0,"    FF: BFG2","");
        legendRAADirGammaPi02050->AddEntry(graphRAAPi0CombPbPb2760GeVSysErr_2050,"#pi^{0} combined","fp");
        legendRAADirGammaPi02050->AddEntry(graphRAAEtaCombPbPb2760GeVSysErr_2050,"#eta combined","fp");
        legendRAADirGammaPi02050->Draw();

        legendRAADirGammaPi0pPb0100->Draw();

        if (graphPCMRAADirGammaSys2050Plot2){
            while(graphPCMRAADirGammaSys2050Plot2->GetX()[0]<2.3) graphPCMRAADirGammaSys2050Plot2->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaSys2050Plot2, markerStyleComb2050,markerSizeComb2050, colorNLOMcGill , colorNLOMcGill, widthLinesBoxes, kTRUE);
            graphPCMRAADirGammaSys2050Plot2->Draw("E2same");
        }
        if (graphPCMRAADirGammaStat2050Plot2){
            while(graphPCMRAADirGammaStat2050Plot2->GetX()[0]<2.3) graphPCMRAADirGammaStat2050Plot2->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphPCMRAADirGammaStat2050Plot2, markerStyleComb2050,markerSizeComb2050, colorNLOMcGill , colorNLOMcGill);
            graphPCMRAADirGammaStat2050Plot2->Draw("p,same,e1Z");
        }

        DrawGammaSetMarkerTGraphAsym(graphRAAPi0CombPbPb2760GeVSysErr_2050, markerStyleComb2050,markerSizeComb2050, colorComb2050-7 , colorComb2050-7, widthLinesBoxes, kTRUE);
        graphRAAPi0CombPbPb2760GeVSysErr_2050->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRAAPi0CombPbPb2760GeVStatErr_2050, markerStyleComb2050,markerSizeComb2050, colorComb2050-7 , colorComb2050-7);
        graphRAAPi0CombPbPb2760GeVStatErr_2050->Draw("p,same,e1Z");

        DrawGammaSetMarkerTGraphAsym(graphRAAEtaCombPbPb2760GeVSysErr_2050, markerStyleComb2050,markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
        graphRAAEtaCombPbPb2760GeVSysErr_2050->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRAAEtaCombPbPb2760GeVStatErr_2050, markerStyleComb2050,markerSizeComb2050, colorComb2050 , colorComb2050);
        graphRAAEtaCombPbPb2760GeVStatErr_2050->Draw("p,same,e1Z");

        CombinedPi0RpPbSystErr->Draw("E2same");
        CombinedPi0RpPbStatErr->Draw("p,same,e1Z");

        histo2DRAADummy2->Draw("axis,same");
    canvasRAADirGammaPi0_0010->Update();
    canvasRAADirGammaPi0_0010->Print(Form("%s/DirGammaPi0PbPbandpPbRAA_2050.%s",outputDir.Data(),suffix.Data()));

    //*******************************************************************************************************************************************
    //*************************************************** Plotting direct Gamma Spectrum PCM only ***********************************************
    //*******************************************************************************************************************************************
    cout << "Plotting direct gamma PHOS only " << __LINE__ << endl;
    TGraphAsymmErrors* graphPHOSDirGammaSpectrumStatErr0010Plot = NULL;
    TGraphAsymmErrors* graphPHOSDirGammaSpectrumSystErr0010Plot = NULL;
    if (histoPHOSDirGammaStatErr0010){
        graphPHOSDirGammaSpectrumStatErr0010Plot         = new TGraphAsymmErrors(histoPHOSDirGammaStatErr0010);
        graphPHOSDirGammaSpectrumStatErr0010Plot         = ScaleGraph(graphPHOSDirGammaSpectrumStatErr0010Plot,100);
    }
    if (histoPHOSDirGammaSysErr0010){
        graphPHOSDirGammaSpectrumSystErr0010Plot         = new TGraphAsymmErrors(histoPHOSDirGammaSysErr0010);
        graphPHOSDirGammaSpectrumSystErr0010Plot         = ScaleGraph(graphPHOSDirGammaSpectrumSystErr0010Plot,100);
    }
    TGraphAsymmErrors* graphPHOSDirGammaSpectrumStatErr2040Plot = NULL;
    TGraphAsymmErrors* graphPHOSDirGammaSpectrumSystErr2040Plot = NULL;
    if (histoPHOSDirGammaStatErr2040){
        graphPHOSDirGammaSpectrumStatErr2040Plot         = new TGraphAsymmErrors(histoPHOSDirGammaStatErr2040);
        graphPHOSDirGammaSpectrumStatErr2040Plot         = ScaleGraph(graphPHOSDirGammaSpectrumStatErr2040Plot,10);
    }
    if (histoPHOSDirGammaSysErr2040){
        graphPHOSDirGammaSpectrumSystErr2040Plot         = new TGraphAsymmErrors(histoPHOSDirGammaSysErr2040);
        graphPHOSDirGammaSpectrumSystErr2040Plot         = ScaleGraph(graphPHOSDirGammaSpectrumSystErr2040Plot,10);
    }

    canvasDirGammaIndMeas->SetLogx(1);
    canvasDirGammaIndMeas->cd();
    dummyDirGammaIndMeas->DrawCopy();

        graphTheoryEPS090010Plot->Draw("p3lsame");
        graphTheoryCT100010Plot->Draw("p3lsame");
        graphTheoryNLO0010Plot->Draw("p3lsame");
        //graphTheoryPromptMcGill0020Plot->Draw("p3lsame");

        if (graphPHOSDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPHOSDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPHOSDirGammaSpectrumSystErr0010Plot->Draw("E2same");
        }
        if (graphPHOSDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPHOSDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPHOSDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }

        graphTheoryEPS092040Plot->Draw("p3lsame");
        graphTheoryCT102040Plot->Draw("p3lsame");
        graphTheoryNLO2040Plot->Draw("p3lsame");
        //graphTheoryPromptMcGill2040Plot->Draw("p3lsame");

        if (graphPHOSDirGammaSpectrumSystErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPHOSDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPHOSDirGammaSpectrumSystErr2040Plot->Draw("E2same");
        }
        if (graphPHOSDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPHOSDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPHOSDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
        }
        graphTheoryEPS092050Plot->Draw("p3lsame");
        graphTheoryCT102050Plot->Draw("p3lsame");
        graphTheoryNLO2050Plot->Draw("p3lsame");


        labelScalingDirGamma0010->Draw();
        labelScalingDirGamma2040->Draw();
        labelScalingDirGamma2050->Draw();
        labelDirGammaColl->Draw();

        TLegend* legendDirGammaPHOSOnly = new TLegend(0.24,0.93-1.*0.85*textsizeLabelsDirGamma*3,0.24+0.21,0.93);
        legendDirGammaPHOSOnly->SetFillStyle(0);
        legendDirGammaPHOSOnly->SetFillColor(0);
        legendDirGammaPHOSOnly->SetLineColor(0);
        legendDirGammaPHOSOnly->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaPHOSOnly->SetMargin(0.2);
        legendDirGammaPHOSOnly->SetTextFont(42);
        legendDirGammaPHOSOnly->AddEntry(graphPCMDirGammaSpectrumSystErr0010Plot,"  0-10% PHOS","pf");
        legendDirGammaPHOSOnly->AddEntry(graphPCMDirGammaSpectrumSystErr2040Plot,"20-40% PHOS","pf");
        legendDirGammaPHOSOnly->AddEntry(graphPCMDirGammaSpectrumSystErr2050Plot,"20-50% PHOS","pf");
        legendDirGammaPHOSOnly->Draw();
        legendDirGammaNLO->Draw();
        legendDirGammaNLO2->Draw();

//     canvasDirGammaIndMeas->Print(Form("%s/DirGammaSpectrumPHOSonly_wihtTheory.%s",outputDir.Data(),suffix.Data()));


    cout << "Plotting PHOS DR with NLO at " << __LINE__ << endl;
    //*******************************************************************************************************************************************
    //********************************************** DR plot with PHOS measurement + NLO  ********************************************************
    //*******************************************************************************************************************************************
    canvasRatioIndDR->cd();
    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();
    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPHOSDRPi0FitStatErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR0010, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        graphTheoryNLODR0010->Draw("pE3lsame");
        graphPHOSDRPi0FitSysErr0010->Draw("E2same");
        histoPHOSDRPi0FitStatErr0010->Draw("p,same,e0,X0");

        TLegend* legendDRPHOSNLO0010 = new TLegend(0.15,0.82-1.1*0.85*textsizeLabelsPad1*3,0.5,0.82);
        legendDRPHOSNLO0010->SetFillStyle(0);
        legendDRPHOSNLO0010->SetFillColor(0);
        legendDRPHOSNLO0010->SetLineColor(0);
        legendDRPHOSNLO0010->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRPHOSNLO0010->SetMargin(0.2);
        legendDRPHOSNLO0010->SetTextFont(42);
        legendDRPHOSNLO0010->AddEntry(graphPHOSDRPi0FitSysErr0010,"PHOS","pf");
        legendDRPHOSNLO0010->AddEntry(graphTheoryNLODR0010,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
        legendDRPHOSNLO0010->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
        legendDRPHOSNLO0010->Draw();

        labelDRCent0010->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPHOSDRPi0FitStatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR2040, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        graphTheoryNLODR2040->Draw("p3lsame");
        graphPHOSDRPi0FitSysErr2040->Draw("E2same");
        histoPHOSDRPi0FitStatErr2040->Draw("p,same,e0,X0");

        TLegend* legendDRPHOSNLO2040 = new TLegend(0.15,0.85-1.1*0.85*textsizeLabelsPad2*3,0.5,0.85);
        legendDRPHOSNLO2040->SetFillStyle(0);
        legendDRPHOSNLO2040->SetFillColor(0);
        legendDRPHOSNLO2040->SetLineColor(0);
        legendDRPHOSNLO2040->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRPHOSNLO2040->SetMargin(0.2);
        legendDRPHOSNLO2040->SetTextFont(42);
        legendDRPHOSNLO2040->AddEntry(graphPHOSDRPi0FitSysErr2040,"PHOS","pf");
        legendDRPHOSNLO2040->AddEntry(graphTheoryNLODR2040,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
        legendDRPHOSNLO2040->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
        legendDRPHOSNLO2040->Draw();

        labelDRCent2040->Draw();

    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPHOSDRPi0FitStatErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR2050, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
//         graphTheoryNLODR2050->Draw("p3lsame");
//         graphPHOSDRPi0FitSysErr2050->Draw("E2same");
//         histoPHOSDRPi0FitStatErr2050->Draw("p,same,e0,X0");

//         TLegend* legendDRPHOSNLO2050 = new TLegend(0.15,0.87-1.1*0.85*textsizeLabelsPad3*3,0.5,0.87);
//         legendDRPHOSNLO2050->SetFillStyle(0);
//         legendDRPHOSNLO2050->SetFillColor(0);
//         legendDRPHOSNLO2050->SetLineColor(0);
//         legendDRPHOSNLO2050->SetTextSize(0.85*textsizeLabelsPad3);
//         legendDRPHOSNLO2050->SetMargin(0.2);
//         legendDRPHOSNLO2050->SetTextFont(42);
//         legendDRPHOSNLO2050->AddEntry(graphPHOSDRPi0FitSysErr2050,"PHOS","pf");
//         legendDRPHOSNLO2050->AddEntry(graphTheoryNLODR2050,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
//         legendDRPHOSNLO2050->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
//         legendDRPHOSNLO2050->Draw();

        labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_PHOSMeasurementTheory.%s", outputDir.Data(), suffix.Data()));

    //*******************************************************************************************************************************************
    //********************************************** DR plot with PHOS measurements all cents in one plot  ***************************************
    //*******************************************************************************************************************************************
    canvasDoubleRatio->cd();
    dummyDR->DrawCopy();

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        graphPHOSDRPi0FitSysErr0010->Draw("E2same");
        histoPHOSDRPi0FitStatErr0010->Draw("p,same,e0,X0");
        graphPHOSDRPi0FitSysErr2040->Draw("E2same");
        histoPHOSDRPi0FitStatErr2040->Draw("p,same,e0,X0");
        graphPHOSDRPi0FitSysErr2050->Draw("E2same");
        histoPHOSDRPi0FitStatErr2050->Draw("p,same,e0,X0");

        TLegend* legendDRAllCentPHOS = GetAndSetLegend(0.15,0.75,4);
        legendDRAllCentPHOS->AddEntry(graphPHOSDRPi0FitSysErr0010,"  0-10% PHOS","pf");
        legendDRAllCentPHOS->AddEntry(graphPHOSDRPi0FitSysErr2040,"20-40% PHOS","pf");
        legendDRAllCentPHOS->AddEntry(graphPHOSDRPi0FitSysErr2050,"20-50% PHOS","pf");
        legendDRAllCentPHOS->Draw();

    canvasDoubleRatio->Print(Form("%s/DR_PHOSMeasurementAllCentsInOne.%s", outputDir.Data(), suffix.Data()));

    if (enablepValueCalcPCM){

        //Int_t nbins = 2;
        for(Int_t nbins=2; nbins<6; nbins++){
            //*********************************************************************************
            // Calculating Significance of combined double ratio for 0-10%
            //*********************************************************************************
            Int_t minBinSig0010 = 0;
            Int_t maxBinSig0010 = nbins;
            const Int_t nPointSig0010 = maxBinSig0010-minBinSig0010 +1;

            // filling of graphs
            TGraph* graphAbsTypeAPlusStatErr0010                     = new TGraph(nPointSig0010);
            TGraph* graphRelTypeAPlusStatErr0010                     = new TGraph(nPointSig0010);
            TGraph* graphRelTypeBErr0010                             = new TGraph(nPointSig0010);
            TGraphErrors* graphDoubleRatioTypeAPlusStatErr0010       = new TGraphErrors(nPointSig0010);
            for (Int_t i = 0; i < nPointSig0010; i++){
                Double_t pt             = graphPCMDRPi0FitSysBErr0010->GetX()[minBinSig0010+i];
                Double_t R              = graphPCMDRPi0FitSysBErr0010->GetY()[minBinSig0010+i];
                Double_t statErr        = graphPCMDRPi0FitStatErr0010->GetEYlow()[minBinSig0010+i];
                Double_t sysAErr        = graphPCMDRPi0FitSysAErr0010->GetEYlow()[minBinSig0010+i];
                Double_t relErrSysA     = sysAErr/ R;
                Double_t statSysAErr    = TMath::Sqrt(statErr*statErr+sysAErr*sysAErr);
                Double_t relErrStatSysA = statSysAErr/ R;
                Double_t relErrSysB     = graphPCMDRPi0FitSysBErr0010->GetEYlow()[minBinSig0010+i]/ graphPCMDRPi0FitSysBErr0010->GetY()[minBinSig0010+i];

                graphAbsTypeAPlusStatErr0010->SetPoint(i, pt, statSysAErr);
                graphRelTypeAPlusStatErr0010->SetPoint(i, pt, relErrStatSysA);
                graphRelTypeBErr0010->SetPoint(i, pt, relErrSysB);
                graphDoubleRatioTypeAPlusStatErr0010->SetPoint(i, pt, R);
                graphDoubleRatioTypeAPlusStatErr0010->SetPointError(i,0,statSysAErr);
                cout << "pt: " << pt << endl;
            }
            Double_t relErrSysC0010 = graphPCMDRPi0FitSysCErr0010->GetEYlow()[minBinSig0010]/ graphPCMDRPi0FitSysCErr0010->GetY()[minBinSig0010];

            // Generating pseudo data
            TH1F* histoChi2NullHypo0010 = new TH1F("histoChi2NullHypo0010", "histoChi2NullHypo0010", 3000, 0., 1500);
            FillChi2HistForNullHypo(100000000, histoChi2NullHypo0010, graphRelTypeAPlusStatErr0010, graphRelTypeBErr0010, relErrSysC0010, kFALSE);

            // calculating significance
            Double_t chi2Data0010     = Chi2ForNullHypo(graphDoubleRatioTypeAPlusStatErr0010);
            Int_t binChi2Data0010     = histoChi2NullHypo0010->FindBin(chi2Data0010);
            Int_t binFirst0010         = 1;
            Int_t binLast0010        = histoChi2NullHypo0010->GetNbinsX();
            Double_t intTot0010        = histoChi2NullHypo0010->Integral(binFirst0010, binLast0010);
            Double_t intData0010    = histoChi2NullHypo0010->Integral(binChi2Data0010, binLast0010);
            Double_t pValue0010        = intData0010/ intTot0010;
            Double_t nSigma0010        = PValueToNSigma(pValue0010);

            cout << "0-10%: chi2data =" << chi2Data0010 << "\t, pVal = " << pValue0010 << "\t, nSigma = " << nSigma0010 << endl;

            // preparing histo for drawing
            TH1F* histoChi2NullData0010 = (TH1F*)histoChi2NullHypo0010->Clone("histoChi2NullData0010");
            for (Int_t i = 1; i < binChi2Data0010; i++ ){
                histoChi2NullData0010->SetBinContent(i,-1);
            }
            for (Int_t i = 1; i < binLast0010; i++ ){
                histoChi2NullData0010->SetBinError(i,0);
            }

            //*********************************************************************************
            // Calculating Significance of combined double ratio for 20-40%
            //*********************************************************************************
            Int_t minBinSig2040 = 0;
            Int_t maxBinSig2040 = nbins;
            const Int_t nPointSig2040 = maxBinSig2040- minBinSig2040 +1;

            // filling of graphs
            TGraph* graphAbsTypeAPlusStatErr2040                     = new TGraph(nPointSig2040);
            TGraph* graphRelTypeAPlusStatErr2040                     = new TGraph(nPointSig2040);
            TGraph* graphRelTypeBErr2040                             = new TGraph(nPointSig2040);
            TGraphErrors* graphDoubleRatioTypeAPlusStatErr2040         = new TGraphErrors(nPointSig2040);

            for (Int_t i = 0; i < nPointSig2040; i++){
                Double_t pt             = graphPCMDRPi0FitSysBErr2040->GetX()[minBinSig2040+i];
                Double_t R                = graphPCMDRPi0FitSysBErr2040->GetY()[minBinSig2040+i];
                Double_t statErr        = graphPCMDRPi0FitStatErr2040->GetEYlow()[minBinSig2040+i];
                Double_t sysAErr        = graphPCMDRPi0FitSysAErr2040->GetEYlow()[minBinSig2040+i];
                Double_t relErrSysA     = sysAErr/ R;
                Double_t statSysAErr    = TMath::Sqrt(statErr*statErr+sysAErr*sysAErr);
                Double_t relErrStatSysA = statSysAErr/ R;
                Double_t relErrSysB        = graphPCMDRPi0FitSysBErr2040->GetEYlow()[minBinSig2040+i]/ graphPCMDRPi0FitSysBErr2040->GetY()[minBinSig2040+i];

                graphAbsTypeAPlusStatErr2040->SetPoint(i, pt, statSysAErr);
                graphRelTypeAPlusStatErr2040->SetPoint(i, pt, relErrStatSysA);
                graphRelTypeBErr2040->SetPoint(i, pt, relErrSysB);
                graphDoubleRatioTypeAPlusStatErr2040->SetPoint(i, pt, R);
                graphDoubleRatioTypeAPlusStatErr2040->SetPointError(i,0,statSysAErr);
                cout << "pt: " << pt << endl;
            }
            Double_t relErrSysC2040 = graphPCMDRPi0FitSysCErr2040->GetEYlow()[minBinSig2040]/ graphPCMDRPi0FitSysCErr2040->GetY()[minBinSig2040];

            // Generating pseudo data
            TH1F* histoChi2NullHypo2040 = new TH1F("histoChi2NullHypo2040", "histoChi2NullHypo2040", 3000, 0., 1500);
            FillChi2HistForNullHypo(100000000, histoChi2NullHypo2040, graphRelTypeAPlusStatErr2040, graphRelTypeBErr2040, relErrSysC2040, kFALSE);

            // calculating significance
            Double_t chi2Data2040     = Chi2ForNullHypo(graphDoubleRatioTypeAPlusStatErr2040);
            Int_t binChi2Data2040     = histoChi2NullHypo2040->FindBin(chi2Data2040);
            Int_t binFirst2040         = 1;
            Int_t binLast2040        = histoChi2NullHypo2040->GetNbinsX();
            Double_t intTot2040        = histoChi2NullHypo2040->Integral(binFirst2040, binLast2040);
            Double_t intData2040    = histoChi2NullHypo2040->Integral(binChi2Data2040, binLast2040);
            Double_t pValue2040        = intData2040/ intTot2040;
            Double_t nSigma2040        = PValueToNSigma(pValue2040);

            cout << "20-40%:  chi2data =" << chi2Data2040 << "\t, pVal = " << pValue2040 << "\t, nSigma = " << nSigma2040 << endl;

            // preparing histo for drawing
            TH1F* histoChi2NullData2040 = (TH1F*)histoChi2NullHypo2040->Clone("histoChi2NullData2040");
            for (Int_t i = 1; i < binChi2Data2040; i++ ){
                histoChi2NullData2040->SetBinContent(i,-1);
            }
            for (Int_t i = 1; i < binLast2040; i++ ){
                histoChi2NullData2040->SetBinError(i,0);
            }

            //*********************************************************************************
            // Calculating Significance of combined double ratio for 20-50%
            //*********************************************************************************

            Int_t minBinSig2050 = 0;
            Int_t maxBinSig2050 = nbins;
            const Int_t nPointSig2050 = maxBinSig2050- minBinSig2050 +1;

            // filling of graphs
            TGraph* graphAbsTypeAPlusStatErr2050                     = new TGraph(nPointSig2050);
            TGraph* graphRelTypeAPlusStatErr2050                     = new TGraph(nPointSig2050);
            TGraph* graphRelTypeBErr2050                             = new TGraph(nPointSig2050);
            TGraphErrors* graphDoubleRatioTypeAPlusStatErr2050         = new TGraphErrors(nPointSig2050);

            for (Int_t i = 0; i < nPointSig2050; i++){
                Double_t pt             = graphPCMDRPi0FitSysBErr2050->GetX()[minBinSig2050+i];
                Double_t R                = graphPCMDRPi0FitSysBErr2050->GetY()[minBinSig2050+i];
                Double_t statErr        = graphPCMDRPi0FitStatErr2050->GetEYlow()[minBinSig2050+i];
                Double_t sysAErr        = graphPCMDRPi0FitSysAErr2050->GetEYlow()[minBinSig2050+i];
                Double_t relErrSysA     = sysAErr/ R;
                Double_t statSysAErr    = TMath::Sqrt(statErr*statErr+sysAErr*sysAErr);
                Double_t relErrStatSysA = statSysAErr/ R;
                Double_t relErrSysB        = graphPCMDRPi0FitSysBErr2050->GetEYlow()[minBinSig2050+i]/ graphPCMDRPi0FitSysBErr2050->GetY()[minBinSig2050+i];

                graphAbsTypeAPlusStatErr2050->SetPoint(i, pt, statSysAErr);
                graphRelTypeAPlusStatErr2050->SetPoint(i, pt, relErrStatSysA);
                graphRelTypeBErr2050->SetPoint(i, pt, relErrSysB);
                graphDoubleRatioTypeAPlusStatErr2050->SetPoint(i, pt, R);
                graphDoubleRatioTypeAPlusStatErr2050->SetPointError(i,0,statSysAErr);
                cout << "pt: " << pt << endl;

            }
            Double_t relErrSysC2050 = graphPCMDRPi0FitSysCErr2050->GetEYlow()[minBinSig2050]/ graphPCMDRPi0FitSysCErr2050->GetY()[minBinSig2050];

            // Generating pseudo data
            TH1F* histoChi2NullHypo2050 = new TH1F("histoChi2NullHypo2050", "histoChi2NullHypo2050", 3000, 0., 1500);
            FillChi2HistForNullHypo(100000000, histoChi2NullHypo2050, graphRelTypeAPlusStatErr2050, graphRelTypeBErr2050, relErrSysC2050, kFALSE);

            // calculating significance
            Double_t chi2Data2050     = Chi2ForNullHypo(graphDoubleRatioTypeAPlusStatErr2050);
            Int_t binChi2Data2050     = histoChi2NullHypo2050->FindBin(chi2Data2050);
            Int_t binFirst2050         = 1;
            Int_t binLast2050        = histoChi2NullHypo2050->GetNbinsX();
            Double_t intTot2050        = histoChi2NullHypo2050->Integral(binFirst2050, binLast2050);
            Double_t intData2050    = histoChi2NullHypo2050->Integral(binChi2Data2050, binLast2050);
            Double_t pValue2050        = intData2050/ intTot2050;
            Double_t nSigma2050        = PValueToNSigma(pValue2050);
            cout << "20-50%:  chi2data =" << chi2Data2050 << "\t, pVal = " << pValue2050 << "\t, nSigma = " << nSigma2050 << endl;
            cout << "**********************************************************************************" << endl;

            // preparing histo for drawing
//             TH1F* histoChi2NullData2050 = (TH1F*)histoChi2NullHypo2050->Clone("histoChi2NullData2050");
//             for (Int_t i = 1; i < binChi2Data2050; i++ ){
//                 histoChi2NullData2050->SetBinContent(i,-1);
//             }
//             for (Int_t i = 1; i < binLast2050; i++ ){
//                 histoChi2NullData2050->SetBinError(i,0);
//             }

            //*********************************************************************************
            // Plotting Significance of combined double ratio for 0-10%
            //*********************************************************************************

//             TCanvas* canvasSignificance = new TCanvas("canvasSignificance","",200,10,1400,1100);  // gives the page size
//             DrawGammaCanvasSettings( canvasSignificance,  0.09, 0.01, 0.015, 0.1);
//             canvasSignificance->SetLogy(1);
//
//             Int_t textSizeLabelsPixelSignificance = 48;
//             Double_t textsizeLabelsSignificance = 0;
//             if (canvasSignificance->XtoPixel(canvasSignificance->GetX2()) < canvasSignificance->YtoPixel(canvasSignificance->GetY1())){
//                 textsizeLabelsSignificance = (Double_t)textSizeLabelsPixelSignificance/canvasSignificance->XtoPixel(canvasSignificance->GetX2()) ;
//             } else {
//                 textsizeLabelsSignificance = (Double_t)textSizeLabelsPixelSignificance/canvasSignificance->YtoPixel(canvasSignificance->GetY1());
//             }
//
//             SetStyleHistoTH1ForGraphs(histoChi2NullHypo0010, "Test statistics #it{t}","Counts",
//                                     0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
//                                     0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
//                                     0.95, 1., 510, 510);
//             histoChi2NullHypo0010->GetYaxis()->SetLabelOffset(0.005);
//
//             Int_t firstAbove0010 = histoChi2NullHypo0010->FindFirstBinAbove();
//             Int_t lastAbove0010 = histoChi2NullHypo0010->FindLastBinAbove();
//
//             histoChi2NullHypo0010->GetXaxis()->SetRange(1,lastAbove0010+1);
//             DrawGammaSetMarker(histoChi2NullHypo0010, 1, 0.1, kGray+2 , kGray+2);
//             histoChi2NullHypo0010->DrawCopy("");
//
//             DrawGammaSetMarker(histoChi2NullData0010, 1, 0.1, colorComb0010Box , colorComb0010Box);
//             histoChi2NullData0010->SetFillColor(colorComb0010Box);
//             histoChi2NullData0010->SetFillStyle(3356);
//             histoChi2NullData0010->Draw("same,lf");
//             histoChi2NullHypo0010->DrawCopy("same");
//             histoChi2NullHypo0010->DrawCopy("same,axis");
//
//                 TLatex *labelSignificanceEnergy0010 = new TLatex(0.58,0.92,collisionSystemCent0010.Data());
//                 SetStyleTLatex( labelSignificanceEnergy0010, 0.85*textsizeLabelsSignificance,4);
//                 labelSignificanceEnergy0010->Draw();
//
//                 TLatex *labelSignificanceALICE = new TLatex(0.58,0.87,"ALICE pseudo data");
//                 SetStyleTLatex( labelSignificanceALICE, 0.85*textsizeLabelsSignificance,4);
//                 labelSignificanceALICE->Draw();
//
//                 DrawGammaLines(chi2Data0010, chi2Data0010, 0, histoChi2NullHypo0010->GetMaximum()*0.1, 3, colorComb0010, 7);
//
//                 TLatex *labelTData0010 = new TLatex(chi2Data0010+10,histoChi2NullHypo0010->GetMaximum()*0.1,"#it{t}_{data}");
//                 SetStyleTLatex( labelTData0010, 0.85*textsizeLabelsSignificance,4,colorComb0010,42,kFALSE);
//                 labelTData0010->Draw();
//
//                 TLatex *labelSignificancePValue0010 = new TLatex(0.48,0.57,Form("#it{p}-value = %1.4f",pValue0010));
//                 SetStyleTLatex( labelSignificancePValue0010, 0.85*textsizeLabelsSignificance,4);
//                 labelSignificancePValue0010->Draw();
//                 TLatex *labelSignificanceNSigma0010 = new TLatex(0.48,0.53,Form("%1.2f #sigma",nSigma0010));
//                 SetStyleTLatex( labelSignificanceNSigma0010, 0.85*textsizeLabelsSignificance,4);
//                 labelSignificanceNSigma0010->Draw();
//
//
//             canvasSignificance->Update();
//             canvasSignificance->Print(Form("%s/DirGammaSignificanceTest_0010.%s",outputDir.Data(),suffix.Data()));
//
//             //*********************************************************************************
//             // Plotting Significance of combined double ratio for 20-40%
//             //*********************************************************************************
//
//             SetStyleHistoTH1ForGraphs(histoChi2NullHypo2040, "Test statistics #it{t}","Counts",
//                                     0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
//                                     0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
//                                     0.95, 1., 510, 510);
//             histoChi2NullHypo2040->GetYaxis()->SetLabelOffset(0.005);
//
//             Int_t firstAbove2040 = histoChi2NullHypo2040->FindFirstBinAbove();
//             Int_t lastAbove2040 = histoChi2NullHypo2040->FindLastBinAbove();
//
//             histoChi2NullHypo2040->GetXaxis()->SetRange(1,lastAbove2040+1);
//             DrawGammaSetMarker(histoChi2NullHypo2040, 1, 0.1, kGray+2 , kGray+2);
//             histoChi2NullHypo2040->DrawCopy("");
//
//             DrawGammaSetMarker(histoChi2NullData2040, 1, 0.1, colorComb2040Box , colorComb2040Box);
//             histoChi2NullData2040->SetFillColor(colorComb2040Box);
//             histoChi2NullData2040->SetFillStyle(3356);
//             histoChi2NullData2040->Draw("same,lf");
//             histoChi2NullHypo2040->DrawCopy("same");
//             histoChi2NullHypo2040->DrawCopy("same,axis");
//
//                 TLatex *labelSignificanceEnergy2040 = new TLatex(0.58,0.92,collisionSystemCent2040.Data());
//                 SetStyleTLatex( labelSignificanceEnergy2040, 0.85*textsizeLabelsSignificance,4);
//                 labelSignificanceEnergy2040->Draw();
//
//                 labelSignificanceALICE->Draw();
//
//                 DrawGammaLines(chi2Data2040, chi2Data2040, 0, histoChi2NullHypo2040->GetMaximum()*0.5, 3, colorComb2040, 7);
//
//                 TLatex *labelTData2040 = new TLatex(chi2Data2040+10,histoChi2NullHypo2040->GetMaximum()*0.5,"#it{t}_{data}");
//                 SetStyleTLatex( labelTData2040, 0.85*textsizeLabelsSignificance,4,colorComb2040,42,kFALSE);
//                 labelTData2040->Draw();
//
//                 TLatex *labelSignificancePValue2040 = new TLatex(0.48,0.57,Form("#it{p}-value = %1.4f",pValue2040));
//                 SetStyleTLatex( labelSignificancePValue2040, 0.85*textsizeLabelsSignificance,4);
//                 labelSignificancePValue2040->Draw();
//                 TLatex *labelSignificanceNSigma2040 = new TLatex(0.48,0.53,Form("%1.2f #sigma",nSigma2040));
//                 SetStyleTLatex( labelSignificanceNSigma2040, 0.85*textsizeLabelsSignificance,4);
//                 labelSignificanceNSigma2040->Draw();
//
//
//             canvasSignificance->Update();
//             canvasSignificance->Print(Form("%s/DirGammaSignificanceTest_2040.%s",outputDir.Data(),suffix.Data()));
//
//             //*********************************************************************************
//             // Plotting Significance of combined double ratio for 20-50%
//             //*********************************************************************************
//
//             SetStyleHistoTH1ForGraphs(histoChi2NullHypo2050, "Test statistics #it{t}","Counts",
//                                     0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
//                                     0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
//                                     0.95, 1., 510, 510);
//             histoChi2NullHypo2050->GetYaxis()->SetLabelOffset(0.005);
//
//             Int_t firstAbove2050 = histoChi2NullHypo2050->FindFirstBinAbove();
//             Int_t lastAbove2050 = histoChi2NullHypo2050->FindLastBinAbove();
//
//             histoChi2NullHypo2050->GetXaxis()->SetRange(1,lastAbove2050+1);
//             histoChi2NullHypo2050->GetYaxis()->SetRangeUser(1, 10*histoChi2NullHypo2050->GetMaximum());
//             DrawGammaSetMarker(histoChi2NullHypo2050, 1, 0.1, kGray+2 , kGray+2);
//             histoChi2NullHypo2050->DrawCopy("");
//
//             DrawGammaSetMarker(histoChi2NullData2050, 1, 0.1, colorComb2050Box , colorComb2050Box);
//             histoChi2NullData2050->SetFillColor(colorComb2050Box);
//             histoChi2NullData2050->SetFillStyle(3356);
//             histoChi2NullData2050->Draw("same,lf");
//             histoChi2NullHypo2050->DrawCopy("same");
//             histoChi2NullHypo2050->DrawCopy("same,axis");
//
//                 TLatex *labelSignificanceEnergy2050 = new TLatex(0.58,0.92,collisionSystemCent2050.Data());
//                 SetStyleTLatex( labelSignificanceEnergy2050, 0.85*textsizeLabelsSignificance,4);
//                 labelSignificanceEnergy2050->Draw();
//
//                 labelSignificanceALICE->Draw();
//
//                 DrawGammaLines(chi2Data2050, chi2Data2050, 0, histoChi2NullHypo2050->GetMaximum()*0.2, 3, colorComb2050, 7);
//
//                 TLatex *labelTData2050 = new TLatex(chi2Data2050+10,histoChi2NullHypo2050->GetMaximum()*0.2,"#it{t}_{data}");
//                 SetStyleTLatex( labelTData2050, 0.85*textsizeLabelsSignificance,4,colorComb2050,42,kFALSE);
//                 labelTData2050->Draw();
//
//                 TLatex *labelSignificancePValue2050 = new TLatex(0.48,0.57,Form("#it{p}-value = %1.4f",pValue2050));
//                 SetStyleTLatex( labelSignificancePValue2050, 0.85*textsizeLabelsSignificance,4);
//                 labelSignificancePValue2050->Draw();
//                 TLatex *labelSignificanceNSigma2050 = new TLatex(0.48,0.53,Form("%1.2f #sigma",nSigma2050));
//                 SetStyleTLatex( labelSignificanceNSigma2050, 0.85*textsizeLabelsSignificance,4);
//                 labelSignificanceNSigma2050->Draw();
//
//
//             canvasSignificance->Update();
//             canvasSignificance->Print(Form("%s/DirGammaSignificanceTest_2050.%s",outputDir.Data(),suffix.Data()));
        }
    }
    if(!combinePHOS){
        const char* fileNameOutputCompPCM = "Gamma_PCMOnlyResults_PbPb_2.76TeV.root";
        TFile* fileGammaSpectrumPCMOnly = new TFile(fileNameOutputCompPCM,"RECREATE");
            fileGammaSpectrumPCMOnly->mkdir("Gamma_PbPb_2.76TeV_0-10%");
            TDirectoryFile* directoryGammaPCMOnly_0010 = (TDirectoryFile*)fileGammaSpectrumPCMOnly->Get("Gamma_PbPb_2.76TeV_0-10%");
            fileGammaSpectrumPCMOnly->cd("Gamma_PbPb_2.76TeV_0-10%");

                if (histoPCMIncGammaStatErr0010) histoPCMIncGammaStatErr0010->Write("hIncGammaSpec_PCM_StatErr",TObject::kOverwrite);
                if (graphPCMIncGammaSysAErr0010) graphPCMIncGammaSysAErr0010->Write("IncGammaSpec_PCM_SysAErr",TObject::kOverwrite);
                if (graphPCMIncGammaSysBErr0010) graphPCMIncGammaSysBErr0010->Write("IncGammaSpec_PCM_SysBErr",TObject::kOverwrite);
                if (graphPCMIncGammaSysCErr0010) graphPCMIncGammaSysCErr0010->Write("IncGammaSpec_PCM_SysCErr",TObject::kOverwrite);

                if (histoPCMDRPi0FitStatErr0010) histoPCMDRPi0FitStatErr0010->Write("hDR_PCM_StatErr",TObject::kOverwrite);
                if (graphPCMDRPi0FitSysErr0010) graphPCMDRPi0FitSysErr0010->Write("DR_PCM_SysTotErr",TObject::kOverwrite);
                if (graphPCMDRPi0FitSysAErr0010) graphPCMDRPi0FitSysAErr0010->Write("DR_PCM_SysAErr",TObject::kOverwrite);
                if (graphPCMDRPi0FitSysBErr0010) graphPCMDRPi0FitSysBErr0010->Write("DR_PCM_SysBErr",TObject::kOverwrite);
                if (graphPCMDRPi0FitSysCErr0010) graphPCMDRPi0FitSysCErr0010->Write("DR_PCM_SysCErr",TObject::kOverwrite);

            fileGammaSpectrumPCMOnly->mkdir("Gamma_PbPb_2.76TeV_20-40%");
            TDirectoryFile* directoryGammaPCMOnly_2040 = (TDirectoryFile*)fileGammaSpectrumPCMOnly->Get("Gamma_PbPb_2.76TeV_20-40%");
            fileGammaSpectrumPCMOnly->cd("Gamma_PbPb_2.76TeV_20-40%");

                if (histoPCMIncGammaStatErr2040) histoPCMIncGammaStatErr2040->Write("hIncGammaSpec_PCM_StatErr",TObject::kOverwrite);
                if (graphPCMIncGammaSysAErr2040) graphPCMIncGammaSysAErr2040->Write("IncGammaSpec_PCM_SysAErr",TObject::kOverwrite);
                if (graphPCMIncGammaSysBErr2040) graphPCMIncGammaSysBErr2040->Write("IncGammaSpec_PCM_SysBErr",TObject::kOverwrite);
                if (graphPCMIncGammaSysCErr2040) graphPCMIncGammaSysCErr2040->Write("IncGammaSpec_PCM_SysCErr",TObject::kOverwrite);

                if (histoPCMDRPi0FitStatErr2040) histoPCMDRPi0FitStatErr2040->Write("hDR_PCM_StatErr",TObject::kOverwrite);
                if (graphPCMDRPi0FitSysAErr2040) graphPCMDRPi0FitSysAErr2040->Write("DR_PCM_SysAErr",TObject::kOverwrite);
                if (graphPCMDRPi0FitSysBErr2040) graphPCMDRPi0FitSysBErr2040->Write("DR_PCM_SysBErr",TObject::kOverwrite);
                if (graphPCMDRPi0FitSysCErr2040) graphPCMDRPi0FitSysCErr2040->Write("DR_PCM_SysCErr",TObject::kOverwrite);

            fileGammaSpectrumPCMOnly->mkdir("Gamma_PbPb_2.76TeV_20-50%");
            TDirectoryFile* directoryGammaPCMOnly_2050 = (TDirectoryFile*)fileGammaSpectrumPCMOnly->Get("Gamma_PbPb_2.76TeV_20-50%");
            fileGammaSpectrumPCMOnly->cd("Gamma_PbPb_2.76TeV_20-50%");

                if (histoPCMIncGammaStatErr2050) histoPCMIncGammaStatErr2050->Write("hIncGammaSpec_PCM_StatErr",TObject::kOverwrite);
                if (graphPCMIncGammaSysAErr2050) graphPCMIncGammaSysAErr2050->Write("IncGammaSpec_PCM_SysAErr",TObject::kOverwrite);
                if (graphPCMIncGammaSysBErr2050) graphPCMIncGammaSysBErr2050->Write("IncGammaSpec_PCM_SysBErr",TObject::kOverwrite);
                if (graphPCMIncGammaSysCErr2050) graphPCMIncGammaSysCErr2050->Write("IncGammaSpec_PCM_SysCErr",TObject::kOverwrite);

                if (histoPCMDRPi0FitStatErr2050) histoPCMDRPi0FitStatErr2050->Write("hDR_PCM_StatErr",TObject::kOverwrite);
                if (graphPCMDRPi0FitSysAErr2050) graphPCMDRPi0FitSysAErr2050->Write("DR_PCM_SysAErr",TObject::kOverwrite);
                if (graphPCMDRPi0FitSysBErr2050) graphPCMDRPi0FitSysBErr2050->Write("DR_PCM_SysBErr",TObject::kOverwrite);
                if (graphPCMDRPi0FitSysCErr2050) graphPCMDRPi0FitSysCErr2050->Write("DR_PCM_SysCErr",TObject::kOverwrite);

        fileGammaSpectrumPCMOnly->Write();
        fileGammaSpectrumPCMOnly->Close();
    }
    if(!combinePHOS) return;



    //*******************************************************************************************************************************************
    //***************************************************** Combine DR of PCM and PHOS **********************************************************
    //*******************************************************************************************************************************************
    Double_t newBinsComb[21]                                    = { 0.9, 1.1, 1.3,
                                                                    1.5, 1.7, 1.9, 2.1, 2.3,
                                                                    2.5, 2.7, 3.,  3.3, 3.7,
                                                                    4.1, 4.6, 5.4, 6.2, 7.,
                                                                    8., 11., 14.};
    TGraphAsymmErrors *graphCombDRPi0FitSysErr0010;
    TGraphAsymmErrors *graphCombDRPi0FitSysAErr0010;
    TGraphAsymmErrors *graphCombDRPi0FitSysBErr0010;
    TGraphAsymmErrors *graphCombDRPi0FitSysCErr0010;
    TGraphAsymmErrors *graphCombDRPi0FitStatErr0010;
    TGraphAsymmErrors *graphCombDRPi0FitSumErr0010;
    graphCombDRPi0FitSumErr0010     = CombinePtPointsSpectraAdv(    histoPCMDRPi0FitStatErr0010, graphPCMDRPi0FitSysErr0010,
                                                                    graphPCMDRPi0FitSysAErr0010, graphPCMDRPi0FitSysBErr0010, graphPCMDRPi0FitSysCErr0010,
                                                                    histoPHOSDRPi0FitStatErr0010, graphPHOSDRPi0FitSysErr0010 ,
                                                                    graphPHOSDRPi0FitSysAErr0010, graphPHOSDRPi0FitSysBErr0010, graphPHOSDRPi0FitSysCErr0010,
                                                                    graphCombDRPi0FitStatErr0010, graphCombDRPi0FitSysErr0010,
                                                                    graphCombDRPi0FitSysAErr0010, graphCombDRPi0FitSysBErr0010, graphCombDRPi0FitSysCErr0010,
                                                                    newBinsComb, 21, 3, 1, -1);
//     graphCombDRPi0FitSysErr0010->RemovePoint(0);
//     graphCombDRPi0FitStatErr0010->RemovePoint(0);
//     graphCombDRPi0FitSumErr0010->RemovePoint(0);
    cout << __LINE__ << endl;
//     graphCombDRPi0FitSumErr0010->Print();
    TGraphAsymmErrors* graphCombDRPi0FitStatSysAErr0010         = AddErrorsOfGraphsQuadratically (graphCombDRPi0FitStatErr0010, graphCombDRPi0FitSysAErr0010);
    Double_t SysCCombDRPi0Fit0010                               = graphCombDRPi0FitSysCErr0010->GetErrorYlow(4)/graphCombDRPi0FitSysCErr0010->GetY()[4];



    TGraphAsymmErrors *graphCombDRPi0FitSysErr2040;
    TGraphAsymmErrors *graphCombDRPi0FitSysAErr2040;
    TGraphAsymmErrors *graphCombDRPi0FitSysBErr2040;
    TGraphAsymmErrors *graphCombDRPi0FitSysCErr2040;
    TGraphAsymmErrors *graphCombDRPi0FitStatErr2040;
    TGraphAsymmErrors *graphCombDRPi0FitSumErr2040;
    graphCombDRPi0FitSumErr2040     = CombinePtPointsSpectraAdv(    histoPCMDRPi0FitStatErr2040, graphPCMDRPi0FitSysErr2040,
                                                                    graphPCMDRPi0FitSysAErr2040, graphPCMDRPi0FitSysBErr2040, graphPCMDRPi0FitSysCErr2040,
                                                                    histoPHOSDRPi0FitStatErr2040, graphPHOSDRPi0FitSysErr2040 ,
                                                                    graphPHOSDRPi0FitSysAErr2040, graphPHOSDRPi0FitSysBErr2040, graphPHOSDRPi0FitSysCErr2040,
                                                                    graphCombDRPi0FitStatErr2040, graphCombDRPi0FitSysErr2040,
                                                                    graphCombDRPi0FitSysAErr2040, graphCombDRPi0FitSysBErr2040, graphCombDRPi0FitSysCErr2040,
                                                                    newBinsComb, 21, 3, 1, -1);
//     graphCombDRPi0FitSysErr2040->RemovePoint(0);
//     graphCombDRPi0FitStatErr2040->RemovePoint(0);
//     graphCombDRPi0FitSumErr2040->RemovePoint(0);
    cout << __LINE__ << endl;
//     graphCombDRPi0FitSumErr2040->Print();
    TGraphAsymmErrors* graphCombDRPi0FitStatSysAErr2040         = AddErrorsOfGraphsQuadratically (graphCombDRPi0FitStatErr2040, graphCombDRPi0FitSysAErr2040);
    Double_t SysCCombDRPi0Fit2040                               = graphCombDRPi0FitSysCErr2040->GetErrorYlow(4)/graphCombDRPi0FitSysCErr2040->GetY()[4];

    TGraphAsymmErrors *graphCombDRPi0FitSysErr2050;
    TGraphAsymmErrors *graphCombDRPi0FitSysAErr2050;
    TGraphAsymmErrors *graphCombDRPi0FitSysBErr2050;
    TGraphAsymmErrors *graphCombDRPi0FitSysCErr2050;
    TGraphAsymmErrors *graphCombDRPi0FitStatErr2050;
    TGraphAsymmErrors *graphCombDRPi0FitSumErr2050;
    graphCombDRPi0FitSumErr2050     = CombinePtPointsSpectraAdv(    histoPCMDRPi0FitStatErr2050, graphPCMDRPi0FitSysErr2050,
                                                                    graphPCMDRPi0FitSysAErr2050, graphPCMDRPi0FitSysBErr2050, graphPCMDRPi0FitSysCErr2050,
                                                                    histoPHOSDRPi0FitStatErr2040, graphPHOSDRPi0FitSysErr2040 ,
                                                                    graphPHOSDRPi0FitSysAErr2040, graphPHOSDRPi0FitSysBErr2040, graphPHOSDRPi0FitSysCErr2040,
                                                                    graphCombDRPi0FitStatErr2050, graphCombDRPi0FitSysErr2050,
                                                                    graphCombDRPi0FitSysAErr2050, graphCombDRPi0FitSysBErr2050, graphCombDRPi0FitSysCErr2050,
                                                                    newBinsComb,  21, 3, 1, -1);
//     graphCombDRPi0FitSysErr2050->RemovePoint(0);
//     graphCombDRPi0FitStatErr2050->RemovePoint(0);
//     graphCombDRPi0FitSumErr2050->RemovePoint(0);
    cout << __LINE__ << endl;
//     graphCombDRPi0FitSumErr2050->Print();
    TGraphAsymmErrors* graphCombDRPi0FitStatSysAErr2050         = AddErrorsOfGraphsQuadratically (graphCombDRPi0FitStatErr2050, graphCombDRPi0FitSysAErr2050);
    Double_t SysCCombDRPi0Fit2050                               = graphCombDRPi0FitSysCErr2050->GetErrorYlow(4)/graphCombDRPi0FitSysCErr2050->GetY()[4];

    //*******************************************************************************************************************************************
    //******************************************** Combine Inclusive gamma of PCM and PHOS ******************************************************
    //*******************************************************************************************************************************************
    //__________________________________________________ 0-10% combine spectra __________________________________________________________________
    TGraphAsymmErrors *graphCombIncGammaSysErr0010;
    TGraphAsymmErrors *graphCombIncGammaSysAErr0010;
    TGraphAsymmErrors *graphCombIncGammaSysBErr0010;
    TGraphAsymmErrors *graphCombIncGammaSysCErr0010;
    TGraphAsymmErrors *graphCombIncGammaStatErr0010;
    TGraphAsymmErrors *graphCombIncGammaSumErr0010;
    graphCombIncGammaSumErr0010     = CombinePtPointsSpectraAdv(    histoPCMIncGammaStatErr0010, graphPCMIncGammaSysErr0010,
                                                                    graphPCMIncGammaSysAErr0010, graphPCMIncGammaSysBErr0010, graphPCMIncGammaSysCErr0010,
                                                                    histoPHOSIncGammaStatErr0010, graphPHOSIncGammaSysErr0010 ,
                                                                    graphPHOSIncGammaSysAErr0010, graphPHOSIncGammaSysBErr0010, graphPHOSIncGammaSysCErr0010,
                                                                    graphCombIncGammaStatErr0010, graphCombIncGammaSysErr0010,
                                                                    graphCombIncGammaSysAErr0010, graphCombIncGammaSysBErr0010, graphCombIncGammaSysCErr0010,
                                                                    newBinsComb, 21, 3, 1, -1);
    cout << __LINE__ << endl;
//     graphCombIncGammaSumErr0010->Print();

    //__________________________________________ 0-10% fit combined and build ratio of individual to fit ________________________________________
    TF1* fitIncGammaCombQCD0010                                 = FitObject("qcd","fitIncGammaCombQCD0010","Photon",graphCombIncGammaSumErr0010,0.9,14,NULL,"QNRMEX0+");
    cout << WriteParameterToFile(fitIncGammaCombQCD0010)<< endl;
    TH1D* histoFitQCDIncGammaComb0010                           = (TH1D*)fitIncGammaCombQCD0010->GetHistogram();

    TGraphAsymmErrors* graphPCMIncGammaStatSysAErr0010          = AddErrorsOfGraphsQuadratically (graphPCMIncGammaStatErr0010, graphPCMIncGammaSysAErr0010);

    graphPHOSIncGammaStatErr0010->RemovePoint(0);
    graphPHOSIncGammaSysErr0010->RemovePoint(0);
    graphPHOSIncGammaSysAErr0010->RemovePoint(0);
    graphPHOSIncGammaSysBErr0010->RemovePoint(0);
    graphPHOSIncGammaSysCErr0010->RemovePoint(0);
    TGraphAsymmErrors* graphPHOSIncGammaStatSysAErr0010         = AddErrorsOfGraphsQuadratically (graphPHOSIncGammaStatErr0010, graphPHOSIncGammaSysAErr0010);

    TGraphAsymmErrors* graphRatioCombPCMIncGammaStatErr0010     = (TGraphAsymmErrors*) graphPCMIncGammaStatErr0010->Clone();
    TGraphAsymmErrors* graphRatioCombPCMIncGammaSysErr0010      = (TGraphAsymmErrors*) graphPCMIncGammaSysErr0010->Clone();
    graphRatioCombPCMIncGammaStatErr0010                        = CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaStatErr0010, fitIncGammaCombQCD0010);
    graphRatioCombPCMIncGammaSysErr0010                         = CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaSysErr0010, fitIncGammaCombQCD0010);
    TGraphAsymmErrors* graphRatioCombPHOSIncGammaStatErr0010    = (TGraphAsymmErrors*) graphPHOSIncGammaStatErr0010->Clone();
    TGraphAsymmErrors* graphRatioCombPHOSIncGammaSysErr0010     = (TGraphAsymmErrors*) graphPHOSIncGammaSysErr0010->Clone();
    graphRatioCombPHOSIncGammaStatErr0010                       = CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaStatErr0010, fitIncGammaCombQCD0010);
    graphRatioCombPHOSIncGammaSysErr0010                        = CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaSysErr0010, fitIncGammaCombQCD0010);

    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatErr0010   = (TGraphAsymmErrors*) graphPHOSIncGammaStatErr0010->Clone();
    graphRatioCombPPHOSIncGammaStatErr0010                      = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatErr0010, graphCombIncGammaSumErr0010);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysErr0010    = (TGraphAsymmErrors*) graphPHOSIncGammaSysErr0010->Clone();
    graphRatioCombPPHOSIncGammaSysErr0010                       = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysErr0010, graphCombIncGammaSumErr0010);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysBErr0010   = (TGraphAsymmErrors*) graphPHOSIncGammaSysBErr0010->Clone();
    graphRatioCombPPHOSIncGammaSysBErr0010                      = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysBErr0010, graphCombIncGammaSumErr0010);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysCErr0010   = (TGraphAsymmErrors*) graphPHOSIncGammaSysCErr0010->Clone();
    graphRatioCombPPHOSIncGammaSysCErr0010                      = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysCErr0010, graphCombIncGammaSumErr0010);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr0010   = (TGraphAsymmErrors*) graphPHOSIncGammaStatSysAErr0010->Clone();
    graphRatioCombPPHOSIncGammaStatSysAErr0010                  = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatSysAErr0010, graphCombIncGammaSumErr0010);

    TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatErr0010    = (TGraphAsymmErrors*) graphPCMIncGammaStatErr0010->Clone();
    graphRatioCombPPCMIncGammaStatErr0010                       = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatErr0010, graphCombIncGammaSumErr0010);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysErr0010     = (TGraphAsymmErrors*) graphPCMIncGammaSysErr0010->Clone();
    graphRatioCombPPCMIncGammaSysErr0010                        = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysErr0010, graphCombIncGammaSumErr0010);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysBErr0010    = (TGraphAsymmErrors*) graphPCMIncGammaSysBErr0010->Clone();
    graphRatioCombPPCMIncGammaSysBErr0010                       = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysBErr0010, graphCombIncGammaSumErr0010);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysCErr0010    = (TGraphAsymmErrors*) graphPCMIncGammaSysCErr0010->Clone();
    graphRatioCombPPCMIncGammaSysCErr0010                       = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysCErr0010, graphCombIncGammaSumErr0010);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr0010= (TGraphAsymmErrors*) graphPCMIncGammaStatSysAErr0010->Clone();
    graphRatioCombPPCMIncGammaStatSysAErr0010                   = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatSysAErr0010, graphCombIncGammaSumErr0010);

    Double_t SysCPHOSIncGamma0010                               = graphPHOSIncGammaSysCErr0010->GetErrorYlow(4)/graphPHOSIncGammaSysCErr0010->GetY()[4];
    Double_t SysCPCMIncGamma0010                                = graphPCMIncGammaSysCErr0010->GetErrorYlow(4)/graphPCMIncGammaSysCErr0010->GetY()[4];


    TGraphAsymmErrors* graphRatioCombCombFitIncGammaStatErr0010 = (TGraphAsymmErrors*)graphCombIncGammaStatErr0010->Clone();
    TGraphAsymmErrors* graphRatioCombCombFitIncGammaSysErr0010  = (TGraphAsymmErrors*)graphCombIncGammaSysErr0010->Clone();
    graphRatioCombCombFitIncGammaStatErr0010                    = CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaStatErr0010, fitIncGammaCombQCD0010);
    graphRatioCombCombFitIncGammaSysErr0010                     = CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaSysErr0010, fitIncGammaCombQCD0010);

    //__________________________________________________ 20-40% combine spectra _________________________________________________________________
    TGraphAsymmErrors *graphCombIncGammaSysErr2040;
    TGraphAsymmErrors *graphCombIncGammaSysAErr2040;
    TGraphAsymmErrors *graphCombIncGammaSysBErr2040;
    TGraphAsymmErrors *graphCombIncGammaSysCErr2040;
    TGraphAsymmErrors *graphCombIncGammaStatErr2040;
    TGraphAsymmErrors *graphCombIncGammaSumErr2040;
    graphCombIncGammaSumErr2040     = CombinePtPointsSpectraAdv(    histoPCMIncGammaStatErr2040, graphPCMIncGammaSysErr2040,
                                                                    graphPCMIncGammaSysAErr2040, graphPCMIncGammaSysBErr2040, graphPCMIncGammaSysCErr2040,
                                                                    histoPHOSIncGammaStatErr2040, graphPHOSIncGammaSysErr2040 ,
                                                                    graphPHOSIncGammaSysAErr2040, graphPHOSIncGammaSysBErr2040, graphPHOSIncGammaSysCErr2040,
                                                                    graphCombIncGammaStatErr2040, graphCombIncGammaSysErr2040,
                                                                    graphCombIncGammaSysAErr2040, graphCombIncGammaSysBErr2040, graphCombIncGammaSysCErr2040,
                                                                    newBinsComb, 21, 3, 1, -1);
    cout << __LINE__ << endl;
//     graphCombIncGammaSumErr2040->Print();

    //__________________________________________ 20-40% fit combined and build ratio of individual to fit _______________________________________
    TF1* fitIncGammaCombQCD2040                                 = FitObject("qcd","fitIncGammaCombQCD2040","Photon",graphCombIncGammaSumErr2040,0.9,14,NULL,"QNRMEX0+");
    cout << WriteParameterToFile(fitIncGammaCombQCD2040)<< endl;
    TH1D* histoFitQCDIncGammaComb2040                           = (TH1D*)fitIncGammaCombQCD2040->GetHistogram();

    TGraphAsymmErrors* graphPCMIncGammaStatSysAErr2040          = AddErrorsOfGraphsQuadratically (graphPCMIncGammaStatErr2040, graphPCMIncGammaSysAErr2040);

    graphPHOSIncGammaStatErr2040->RemovePoint(0);
    graphPHOSIncGammaSysErr2040->RemovePoint(0);
    graphPHOSIncGammaSysAErr2040->RemovePoint(0);
    graphPHOSIncGammaSysBErr2040->RemovePoint(0);
    graphPHOSIncGammaSysCErr2040->RemovePoint(0);
    TGraphAsymmErrors* graphPHOSIncGammaStatSysAErr2040         = AddErrorsOfGraphsQuadratically (graphPHOSIncGammaStatErr2040, graphPHOSIncGammaSysAErr2040);

    TGraphAsymmErrors* graphRatioCombPHOSIncGammaStatErr2040    = (TGraphAsymmErrors*) graphPHOSIncGammaStatErr2040->Clone();
    TGraphAsymmErrors* graphRatioCombPHOSIncGammaSysErr2040     = (TGraphAsymmErrors*) graphPHOSIncGammaSysErr2040->Clone();
    TGraphAsymmErrors* graphRatioCombPCMIncGammaStatErr2040     = (TGraphAsymmErrors*) graphPCMIncGammaStatErr2040->Clone();
    TGraphAsymmErrors* graphRatioCombPCMIncGammaSysErr2040      = (TGraphAsymmErrors*) graphPCMIncGammaSysErr2040->Clone();
    graphRatioCombPHOSIncGammaStatErr2040                       = CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaStatErr2040, fitIncGammaCombQCD2040);
    graphRatioCombPHOSIncGammaSysErr2040                        = CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaSysErr2040, fitIncGammaCombQCD2040);
    graphRatioCombPCMIncGammaStatErr2040                        = CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaStatErr2040, fitIncGammaCombQCD2040);
    graphRatioCombPCMIncGammaSysErr2040                         = CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaSysErr2040, fitIncGammaCombQCD2040);
    TGraphAsymmErrors* graphRatioCombCombFitIncGammaStatErr2040 = (TGraphAsymmErrors*)graphCombIncGammaStatErr2040->Clone();
    TGraphAsymmErrors* graphRatioCombCombFitIncGammaSysErr2040  = (TGraphAsymmErrors*)graphCombIncGammaSysErr2040->Clone();
    graphRatioCombCombFitIncGammaStatErr2040                    = CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaStatErr2040, fitIncGammaCombQCD2040);
    graphRatioCombCombFitIncGammaSysErr2040                     = CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaSysErr2040, fitIncGammaCombQCD2040);

    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatErr2040   = (TGraphAsymmErrors*) graphPHOSIncGammaStatErr2040->Clone();
    graphRatioCombPPHOSIncGammaStatErr2040                      = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatErr2040, graphCombIncGammaSumErr2040);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysErr2040    = (TGraphAsymmErrors*) graphPHOSIncGammaSysErr2040->Clone();
    graphRatioCombPPHOSIncGammaSysErr2040                       = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysErr2040, graphCombIncGammaSumErr2040);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysBErr2040   = (TGraphAsymmErrors*) graphPHOSIncGammaSysBErr2040->Clone();
    graphRatioCombPPHOSIncGammaSysBErr2040                      = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysBErr2040, graphCombIncGammaSumErr2040);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr2040   = (TGraphAsymmErrors*) graphPHOSIncGammaStatSysAErr2040->Clone();
    graphRatioCombPPHOSIncGammaStatSysAErr2040                  = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatSysAErr2040, graphCombIncGammaSumErr2040);

    TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatErr2040    = (TGraphAsymmErrors*) graphPCMIncGammaStatErr2040->Clone();
    graphRatioCombPPCMIncGammaStatErr2040                       = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatErr2040, graphCombIncGammaSumErr2040);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysErr2040     = (TGraphAsymmErrors*) graphPCMIncGammaSysErr2040->Clone();
    graphRatioCombPPCMIncGammaSysErr2040                        = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysErr2040, graphCombIncGammaSumErr2040);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysBErr2040    = (TGraphAsymmErrors*) graphPCMIncGammaSysBErr2040->Clone();
    graphRatioCombPPCMIncGammaSysBErr2040                       = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysBErr2040, graphCombIncGammaSumErr2040);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr2040= (TGraphAsymmErrors*) graphPCMIncGammaStatSysAErr2040->Clone();
    graphRatioCombPPCMIncGammaStatSysAErr2040                   = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatSysAErr2040, graphCombIncGammaSumErr2040);

    Double_t SysCPHOSIncGamma2040                               = graphPHOSIncGammaSysCErr2040->GetErrorYlow(4)/graphPHOSIncGammaSysCErr2040->GetY()[4];
    Double_t SysCPCMIncGamma2040                                = graphPCMIncGammaSysCErr2040->GetErrorYlow(4)/graphPCMIncGammaSysCErr2040->GetY()[4];


    //__________________________________________________ 20-50% combine spectra _________________________________________________________________
    TGraphAsymmErrors *graphCombIncGammaSysErr2050;
    TGraphAsymmErrors *graphCombIncGammaSysAErr2050;
    TGraphAsymmErrors *graphCombIncGammaSysBErr2050;
    TGraphAsymmErrors *graphCombIncGammaSysCErr2050;
    TGraphAsymmErrors *graphCombIncGammaStatErr2050;
    TGraphAsymmErrors *graphCombIncGammaSumErr2050;
    graphCombIncGammaSumErr2050     = CombinePtPointsSpectraAdv(    histoPCMIncGammaStatErr2050, graphPCMIncGammaSysErr2050,
                                                                    graphPCMIncGammaSysAErr2050, graphPCMIncGammaSysBErr2050, graphPCMIncGammaSysCErr2050,
                                                                    histoPHOSIncGammaStatErr2040, graphPHOSIncGammaSysErr2040 ,
                                                                    graphPHOSIncGammaSysAErr2040, graphPHOSIncGammaSysBErr2040, graphPHOSIncGammaSysCErr2040,
                                                                    graphCombIncGammaStatErr2050, graphCombIncGammaSysErr2050,
                                                                    graphCombIncGammaSysAErr2050, graphCombIncGammaSysBErr2050, graphCombIncGammaSysCErr2050,
                                                                    newBinsComb, 21, 3, 1, -1);
    cout << __LINE__ << endl;
//     graphCombIncGammaSumErr2050->Print();

    //__________________________________________ 20-50% fit combined and build ratio of individual to fit _______________________________________
    TF1* fitIncGammaCombQCD2050                                 = FitObject("qcd","fitIncGammaCombQCD2050","Photon",graphCombIncGammaSumErr2050,0.9,14,NULL,"QNRMEX0+");
    cout << WriteParameterToFile(fitIncGammaCombQCD2050)<< endl;
    TH1D* histoFitQCDIncGammaComb2050                           = (TH1D*)fitIncGammaCombQCD2050->GetHistogram();

    TGraphAsymmErrors* graphPCMIncGammaStatSysAErr2050          = AddErrorsOfGraphsQuadratically (graphPCMIncGammaStatErr2050, graphPCMIncGammaSysAErr2050);

    TGraphAsymmErrors* graphPHOSIncGammaStatErr2050             = new TGraphAsymmErrors(histoPHOSIncGammaStatErr2050);
    graphPHOSIncGammaStatErr2050->RemovePoint(0);
    graphPHOSIncGammaSysErr2050->RemovePoint(0);
    graphPHOSIncGammaSysAErr2050->RemovePoint(0);
    graphPHOSIncGammaSysBErr2050->RemovePoint(0);
    graphPHOSIncGammaSysCErr2050->RemovePoint(0);
    TGraphAsymmErrors* graphPHOSIncGammaStatSysAErr2050         = AddErrorsOfGraphsQuadratically (graphPHOSIncGammaStatErr2050, graphPHOSIncGammaSysAErr2050);

    TGraphAsymmErrors* graphRatioCombPHOSIncGammaStatErr2050    = (TGraphAsymmErrors*) graphPHOSIncGammaStatErr2050->Clone();
    TGraphAsymmErrors* graphRatioCombPHOSIncGammaSysErr2050     = (TGraphAsymmErrors*) graphPHOSIncGammaSysErr2050->Clone();
    TGraphAsymmErrors* graphRatioCombPCMIncGammaStatErr2050     = (TGraphAsymmErrors*) graphPCMIncGammaStatErr2050->Clone();
    TGraphAsymmErrors* graphRatioCombPCMIncGammaSysErr2050      = (TGraphAsymmErrors*) graphPCMIncGammaSysErr2050->Clone();
    graphRatioCombPHOSIncGammaStatErr2050                       = CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaStatErr2050, fitIncGammaCombQCD2050);
    graphRatioCombPHOSIncGammaSysErr2050                        = CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaSysErr2050, fitIncGammaCombQCD2050);
    graphRatioCombPCMIncGammaStatErr2050                        = CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaStatErr2050, fitIncGammaCombQCD2050);
    graphRatioCombPCMIncGammaSysErr2050                         = CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaSysErr2050, fitIncGammaCombQCD2050);
    TGraphAsymmErrors* graphRatioCombCombFitIncGammaStatErr2050 = (TGraphAsymmErrors*)graphCombIncGammaStatErr2050->Clone();
    TGraphAsymmErrors* graphRatioCombCombFitIncGammaSysErr2050  = (TGraphAsymmErrors*)graphCombIncGammaSysErr2050->Clone();
    graphRatioCombCombFitIncGammaStatErr2050                    = CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaStatErr2050, fitIncGammaCombQCD2050);
    graphRatioCombCombFitIncGammaSysErr2050                     = CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaSysErr2050, fitIncGammaCombQCD2050);

    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatErr2050   = (TGraphAsymmErrors*) graphPHOSIncGammaStatErr2050->Clone();
    graphRatioCombPPHOSIncGammaStatErr2050                      = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatErr2050, graphCombIncGammaSumErr2050);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysErr2050    = (TGraphAsymmErrors*) graphPHOSIncGammaSysErr2050->Clone();
    graphRatioCombPPHOSIncGammaSysErr2050                       = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysErr2050, graphCombIncGammaSumErr2050);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysBErr2050   = (TGraphAsymmErrors*) graphPHOSIncGammaSysBErr2050->Clone();
    graphRatioCombPPHOSIncGammaSysBErr2050                      = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysBErr2050, graphCombIncGammaSumErr2050);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr2050   = (TGraphAsymmErrors*) graphPHOSIncGammaStatSysAErr2050->Clone();
    graphRatioCombPPHOSIncGammaStatSysAErr2050                  = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatSysAErr2050, graphCombIncGammaSumErr2050);

    TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatErr2050    = (TGraphAsymmErrors*) graphPCMIncGammaStatErr2050->Clone();
    graphRatioCombPPCMIncGammaStatErr2050                       = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatErr2050, graphCombIncGammaSumErr2050);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysErr2050     = (TGraphAsymmErrors*) graphPCMIncGammaSysErr2050->Clone();
    graphRatioCombPPCMIncGammaSysErr2050                        = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysErr2050, graphCombIncGammaSumErr2050);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysBErr2050    = (TGraphAsymmErrors*) graphPCMIncGammaSysBErr2050->Clone();
    graphRatioCombPPCMIncGammaSysBErr2050                       = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysBErr2050, graphCombIncGammaSumErr2050);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr2050= (TGraphAsymmErrors*) graphPCMIncGammaStatSysAErr2050->Clone();
    graphRatioCombPPCMIncGammaStatSysAErr2050                   = CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatSysAErr2050, graphCombIncGammaSumErr2050);

    Double_t SysCPHOSIncGamma2050                               = graphPHOSIncGammaSysCErr2050->GetErrorYlow(4)/graphPHOSIncGammaSysCErr2050->GetY()[4];
    Double_t SysCPCMIncGamma2050                                = graphPCMIncGammaSysCErr2050->GetErrorYlow(4)/graphPCMIncGammaSysCErr2050->GetY()[4];

    //*******************************************************************************************************************************************
    //**************************************************** Calculate direct photon spectrum 0010 ************************************************
    //*******************************************************************************************************************************************
    cout << endl;
    cout << __LINE__ << endl;
    //graphCombIncGammaStatErr0010->Print();
    Double_t xArrayCombined[graphCombDRPi0FitStatErr0010->GetN()+1];
    xArrayCombined[0] = graphCombDRPi0FitStatErr0010->GetX()[0] - graphCombDRPi0FitStatErr0010->GetEXhigh()[0];
    //cout << "Binning \n" << xArrayCombined[0] << endl;
    for (Int_t i = 1; i<graphCombDRPi0FitStatErr0010->GetN()+1;i++){
        xArrayCombined[i] = graphCombDRPi0FitStatErr0010->GetX()[i-1] + graphCombDRPi0FitStatErr0010->GetEXhigh()[i-1];
        //cout << xArrayCombined[i] << endl;
    }

    //_______________________ copy inclusive photon spectra _____________________
    TH1D *histoCombDirGammaSpectrumErrSum0010                   = new TH1D("histoCombDirGammaSpectrumErrSum0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSys0010                   = new TH1D("histoCombDirGammaSpectrumErrSys0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSysA0010                  = new TH1D("histoCombDirGammaSpectrumErrSysA0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSysB0010                  = new TH1D("histoCombDirGammaSpectrumErrSysB0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSysC0010                  = new TH1D("histoCombDirGammaSpectrumErrSysC0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrStat0010                  = new TH1D("histoCombDirGammaSpectrumErrStat0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);

    //_______________________ get arrays of double ratio errors __________________
    Double_t *SystErrorsCombDR0010                              = new Double_t[graphCombIncGammaStatErr0010->GetN()];
    Double_t *SystAErrorsCombDR0010                             = new Double_t[graphCombIncGammaStatErr0010->GetN()];
    Double_t *SystBErrorsCombDR0010                             = new Double_t[graphCombIncGammaStatErr0010->GetN()];
    Double_t *SystCErrorsCombDR0010                             = new Double_t[graphCombIncGammaStatErr0010->GetN()];
    Double_t *sumErrorsCombDR0010                               = new Double_t[graphCombIncGammaStatErr0010->GetN()];
    Double_t *StatErrorsCombDR0010                              = new Double_t[graphCombIncGammaStatErr0010->GetN()];
    Double_t *xErrorsDR0010                                     = new Double_t[graphCombIncGammaStatErr0010->GetN()];
    for (Int_t i = 0; i< graphCombDRPi0FitStatErr0010->GetN(); i++){
        SystErrorsCombDR0010[i]                                 = graphCombDRPi0FitSysErr0010->GetEYhigh()[i]/graphCombDRPi0FitSysErr0010->GetY()[i] *100;
        SystAErrorsCombDR0010[i]                                = graphCombDRPi0FitSysAErr0010->GetEYhigh()[i]/graphCombDRPi0FitSysAErr0010->GetY()[i] *100;
        SystBErrorsCombDR0010[i]                                = graphCombDRPi0FitSysBErr0010->GetEYhigh()[i]/graphCombDRPi0FitSysBErr0010->GetY()[i] *100;
        SystCErrorsCombDR0010[i]                                = graphCombDRPi0FitSysCErr0010->GetEYhigh()[i]/graphCombDRPi0FitSysCErr0010->GetY()[i] *100;
        StatErrorsCombDR0010[i]                                 = graphCombDRPi0FitStatErr0010->GetEYhigh()[i]/graphCombDRPi0FitStatErr0010->GetY()[i] *100;
        sumErrorsCombDR0010[i]                                  = graphCombDRPi0FitSumErr0010->GetEYhigh()[i]/graphCombDRPi0FitSumErr0010->GetY()[i] *100;
        //cout << i << "\t" << graphCombDRPi0FitSysErr0010->GetY()[i] << "\t" << graphCombDRPi0FitSysErr0010->GetEYhigh()[i] << "\t" <<SystErrorsCombDR0010[i] << endl;
    }
    xErrorsDR0010                                               = graphCombDRPi0FitStatErr0010->GetX();

    cout << __LINE__ << endl;
    //graphCombDRPi0FitSumErr0010->Print();

    //_______________________ copy inclusive photon spectra _____________________
    TH1D* histoCombErrorsForDRSum0010                           = new TH1D("histoCombErrorsForDRSum0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRStat0010                          = new TH1D("histoCombErrorsForDRStat0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSys0010                           = new TH1D("histoCombErrorsForDRSys0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSysA0010                          = new TH1D("histoCombErrorsForDRSysA0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSysB0010                          = new TH1D("histoCombErrorsForDRSysB0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSysC0010                          = new TH1D("histoCombErrorsForDRSysC0010","",graphCombDRPi0FitStatErr0010->GetN(),xArrayCombined);

    for(Int_t i = 1; i<graphCombDRPi0FitStatErr0010->GetN()+1;i++){
        //cout<< i << "\t"<<xErrorsDR0010[i-1]<<"  "<<histoCombErrorsForDRSum0010->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum0010->GetBinWidth(i) <<endl;
        Double_t binErrorSummed                                 = sumErrorsCombDR0010[i-1];
        Double_t binErrorSyst                                   = SystErrorsCombDR0010[i-1];
        Double_t binErrorSystA                                  = SystAErrorsCombDR0010[i-1];
        Double_t binErrorSystB                                  = SystBErrorsCombDR0010[i-1];
        Double_t binErrorSystC                                  = SystCErrorsCombDR0010[i-1];
        Double_t binErrorStat                                   = StatErrorsCombDR0010[i-1];
        Double_t DR                                             = graphCombDRPi0FitStatErr0010->GetY()[i-1];

        //cout << DR << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
        histoCombErrorsForDRSum0010->SetBinContent(i,DR);
        histoCombErrorsForDRSys0010->SetBinContent(i,DR);
        histoCombErrorsForDRSysA0010->SetBinContent(i,DR);
        histoCombErrorsForDRSysB0010->SetBinContent(i,DR);
        histoCombErrorsForDRSysC0010->SetBinContent(i,DR);
        histoCombErrorsForDRStat0010->SetBinContent(i,DR);
        histoCombErrorsForDRSum0010->SetBinError(i,(binErrorSummed/100)*DR);
        histoCombErrorsForDRSys0010->SetBinError(i,(binErrorSyst/100)*DR);
        histoCombErrorsForDRSysA0010->SetBinError(i,(binErrorSystA/100)*DR);
        histoCombErrorsForDRSysB0010->SetBinError(i,(binErrorSystB/100)*DR);
        histoCombErrorsForDRSysC0010->SetBinError(i,(binErrorSystC/100)*DR);
        histoCombErrorsForDRStat0010->SetBinError(i,(binErrorStat/100)*DR);
    }

    for(Int_t i = 1; i<histoCombErrorsForDRSum0010->GetNbinsX()+1;i++){
        histoCombDirGammaSpectrumErrSum0010->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSys0010->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSysA0010->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSysB0010->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSysC0010->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrStat0010->SetBinContent(i+1,-1);

        histoCombDirGammaSpectrumErrSum0010->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSys0010->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSysA0010->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSysB0010->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSysC0010->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrStat0010->SetBinError(i+1,0);
    }

    // get the binning of the direct photons from the DR
    TH1D *histoCombDirGammaSpecSysErr0010                        = new TH1D(*histoCombErrorsForDRSys0010);
    TH1D *histoCombDirGammaSpecSysAErr0010                       = new TH1D(*histoCombErrorsForDRSysA0010);
    TH1D *histoCombDirGammaSpecSysBErr0010                       = new TH1D(*histoCombErrorsForDRSysB0010);
    TH1D *histoCombDirGammaSpecSysCErr0010                       = new TH1D(*histoCombErrorsForDRSysC0010);
    TH1D *histoCombDirGammaSpecStatErr0010                       = new TH1D(*histoCombErrorsForDRStat0010);
    TH1D *histoCombDirGammaSpecSumErr0010                        = new TH1D(*histoCombErrorsForDRSum0010);

    cout << __LINE__ << endl;
    //graphCombDRPi0FitStatErr0010->Print();
    for(Int_t i = 1; i<graphCombDRPi0FitStatErr0010->GetN()+1; i++){
        // obtain common quantities
        Double_t Rgamma                 = histoCombErrorsForDRSys0010->GetBinContent(i);
        Double_t nIncGamma              = graphCombIncGammaStatErr0010->GetY()[i-1];

        // calculating Systematics graph
        Double_t errRgamma              = histoCombErrorsForDRSys0010->GetBinError(i);
        Double_t errNIncGam             = graphCombIncGammaSysErr0010->GetEYhigh()[i-1];
        Double_t q1                     = 1 - 1/ Rgamma;

        Double_t q1Error                = errRgamma/(Rgamma*Rgamma);
        Double_t content                = nIncGamma * ( 1 - 1/ Rgamma);
        Double_t error                  = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        Double_t errDR                  = content - error;
        histoCombDirGammaSpecSysErr0010->SetBinError(i, error);
        histoCombDirGammaSpecSysErr0010->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSys0010->SetBinContent(i, errDR);

        // calculating Systematics A graph
        errRgamma                       = histoCombErrorsForDRSysA0010->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSysAErr0010->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSysAErr0010->SetBinError(i, error);
        histoCombDirGammaSpecSysAErr0010->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSysA0010->SetBinContent(i, errDR);

        // calculating Systematics B graph
        errRgamma                       = histoCombErrorsForDRSysB0010->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSysBErr0010->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSysBErr0010->SetBinError(i, error);
        histoCombDirGammaSpecSysBErr0010->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSysB0010->SetBinContent(i, errDR);

        // calculating Systematics C graph
        errRgamma                       = histoCombErrorsForDRSysC0010->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSysCErr0010->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSysCErr0010->SetBinError(i, error);
        histoCombDirGammaSpecSysCErr0010->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSysC0010->SetBinContent(i, errDR);

        // calculating Stat graphs
        errRgamma                       = histoCombErrorsForDRStat0010->GetBinError(i);
        errNIncGam                      = graphCombIncGammaStatErr0010->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecStatErr0010->SetBinError(i, error);
        histoCombDirGammaSpecStatErr0010->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrStat0010->SetBinContent(i, errDR);

        // calculating summed error graphs
        errRgamma                       = histoCombErrorsForDRSum0010->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSumErr0010->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSumErr0010->SetBinError(i, error);
        histoCombDirGammaSpecSumErr0010->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSum0010->SetBinContent(i, errDR);
    }

    // purely calculating points based on all Systematic errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystErr0010 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSys0010,histoCombDirGammaSpecStatErr0010,0,0.5);
    if(graphCombDirGammaSpectrumSystErr0010)graphCombDirGammaSpectrumSystErr0010->SetName("graphCombDirGammaSpectrumSystErr0010");
    //if(graphCombDirGammaSpectrumSystErr0010)graphCombDirGammaSpectrumSystErr0010->Print();
    //if(graphCombDirGammaSpectrumSystErr0010) cout << "graph has been found" << endl;
    // purely calculating points based on all Systematic errors A
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystAErr0010 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysA0010,histoCombDirGammaSpecStatErr0010,0,0.5);
    if(graphCombDirGammaSpectrumSystAErr0010)graphCombDirGammaSpectrumSystAErr0010->SetName("graphCombDirGammaSpectrumSystAErr0010");
    //if(graphCombDirGammaSpectrumSystAErr0010)graphCombDirGammaSpectrumSystAErr0010->Print();
    // purely calculating points based on all Systematic errors B
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystBErr0010 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysB0010,histoCombDirGammaSpecStatErr0010,0,0.5);
    if(graphCombDirGammaSpectrumSystBErr0010)graphCombDirGammaSpectrumSystBErr0010->SetName("graphCombDirGammaSpectrumSystBErr0010");
    //if(graphCombDirGammaSpectrumSystBErr0010)graphCombDirGammaSpectrumSystBErr0010->Print();
    // purely calculating points based on all Systematic errors C
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystCErr0010 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysC0010,histoCombDirGammaSpecStatErr0010,0,0.5);
    if(graphCombDirGammaSpectrumSystCErr0010)graphCombDirGammaSpectrumSystCErr0010->SetName("graphCombDirGammaSpectrumSystCErr0010");
    //if(graphCombDirGammaSpectrumSystCErr0010)graphCombDirGammaSpectrumSystCErr0010->Print();

    // purely calculating points based on Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumStatErr0010 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrStat0010,histoCombDirGammaSpecStatErr0010,0,0.5);
    if(graphCombDirGammaSpectrumStatErr0010)graphCombDirGammaSpectrumStatErr0010->SetName("graphCombDirGammaSpectrumStatErr0010");
    //if(graphCombDirGammaSpectrumStatErr0010)graphCombDirGammaSpectrumStatErr0010->Print();
    // purely calculating points based on all Systematic + Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr0010 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum0010,histoCombDirGammaSpecStatErr0010,0,0.5);
    if(graphCombDirGammaSpectrumSumErr0010)graphCombDirGammaSpectrumSumErr0010->SetName("graphCombDirGammaSpectrumSumErr0010");
    //if(graphCombDirGammaSpectrumSumErr0010)graphCombDirGammaSpectrumSumErr0010->Print();
    // calculate points above confidence level summed errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr0010Confi = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum0010,histoCombDirGammaSpecStatErr0010,2,0.5);
    if(graphCombDirGammaSpectrumSumErr0010Confi)graphCombDirGammaSpectrumSumErr0010Confi->SetName("graphCombDirGammaSpectrumSumErr0010Confi");
    //if(graphCombDirGammaSpectrumSumErr0010Confi)graphCombDirGammaSpectrumSumErr0010Confi->Print();
    // calculate arrows for points with 0, error summed
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr0010Ar = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum0010,histoCombDirGammaSpecStatErr0010,5,0.5);
    if(graphCombDirGammaSpectrumSumErr0010Ar)graphCombDirGammaSpectrumSumErr0010Ar->SetName("graphCombDirGammaSpectrumSumErr0010Ar");
    //if(graphCombDirGammaSpectrumSumErr0010Ar)graphCombDirGammaSpectrumSumErr0010Ar->Print();

    cout << __LINE__ << endl;

    TF1* fitThermalGamma0010Sum                                 = FitObject("e","fitThermalGamma0010Sum","Photon",histoCombDirGammaSpectrumErrSum0010,0.9,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma0010Sum)<< endl;
    TF1* fitThermalGamma0010Sum23                               = FitObject("e","fitThermalGamma0010Sum23","Photon",histoCombDirGammaSpectrumErrSum0010,0.9,2.3,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma0010Sum23)<< endl;
    TF1* fitThermalGamma0010Stat                                = FitObject("e","fitThermalGamma0010Stat","Photon",graphCombDirGammaSpectrumStatErr0010,0.9,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma0010Stat)<< endl;

    TF1* fitThermalGamma0010Sys                                 = FitObject("e","fitThermalGamma0010Sys","Photon",graphCombDirGammaSpectrumSystErr0010,0.9,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma0010Sys)<< endl;
    TF1* fitThermalGamma0010SysA                                = FitObject("e","fitThermalGamma0010SysA","Photon",graphCombDirGammaSpectrumSystAErr0010,0.9,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma0010SysA)<< endl;
    TF1* fitThermalGamma0010SysB                                = FitObject("e","fitThermalGamma0010SysB","Photon",graphCombDirGammaSpectrumSystBErr0010,0.9,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma0010SysB)<< endl;
    TF1* fitThermalGamma0010SysC                                = FitObject("e","fitThermalGamma0010SysC","Photon",graphCombDirGammaSpectrumSystCErr0010,0.9,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma0010SysC)<< endl;

    TGraphAsymmErrors* graphCombDirGammaSpectrumErrSum0010      = new TGraphAsymmErrors(histoCombDirGammaSpectrumErrSum0010);
    TF1* fitFullDirGamma0010Stat                                = FitObject("qcd","fitFullDirGamma0010Stat","Photon",graphCombDirGammaSpectrumStatErr0010,1.,14,NULL,"NRMEX0+");
    fileFinalResults << WriteParameterToFile(fitFullDirGamma0010Stat)<< endl;
    TF1* fitFullDirGamma0010Sys                                 = FitObject("qcd","fitFullDirGamma0010Sys","Photon",graphCombDirGammaSpectrumSystErr0010,1.,14,NULL,"NRMEX0+");
    fileFinalResults << WriteParameterToFile(fitFullDirGamma0010Sys)<< endl;
    fitFullDirGamma0010Sys->SetParameter(0,fitFullDirGamma0010Stat->GetParameter(0));
    fitFullDirGamma0010Sys->SetParameter(1,fitFullDirGamma0010Stat->GetParameter(1));
    fitFullDirGamma0010Sys->SetParameter(2,fitFullDirGamma0010Stat->GetParameter(2));
    fitFullDirGamma0010Sys->SetParameter(3,fitFullDirGamma0010Stat->GetParameter(3));
    fitFullDirGamma0010Sys->SetParameter(4,fitFullDirGamma0010Stat->GetParameter(4));

    //return;
    TGraphAsymmErrors* graphRatioCombFitDirGammaStatErr0010         = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr0010->Clone();
    TGraphAsymmErrors* graphRatioCombFitDirGammaSysErr0010          = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr0010->Clone();
    graphRatioCombFitDirGammaStatErr0010                            = CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaStatErr0010, fitFullDirGamma0010Stat);
    graphRatioCombFitDirGammaSysErr0010                             = CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaSysErr0010, fitFullDirGamma0010Sys);

    // Calculate thermal spectrum
    TGraphAsymmErrors* graphCombThermalGammaSpectrumStatErr0010     = NULL;
    TF1* fitPureThermalGamma0010Stat                                = NULL;
    TGraphAsymmErrors* graphCombThermalGammaSpectrumSysErr0010      = NULL;
    TGraphAsymmErrors* graphCombThermalGammaSpectrumSumErr0010Ar    = NULL;
    TF1* fitFullThermalGamma0010Sys                                 = NULL;
    TF1* fitFullThermalGamma0010Stat                                = NULL;
    TGraphAsymmErrors* graphRatioCombFitThermalGammaStatErr0010     = NULL;
    TGraphAsymmErrors* graphRatioCombFitThermalGammaSysErr0010      = NULL;
    if (graphCombDirGammaSpectrumStatErr0010){
        graphCombThermalGammaSpectrumStatErr0010                    = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0010, graphCombDirGammaSpectrumStatErr0010);
        graphCombThermalGammaSpectrumStatErr0010->Print();
        fitPureThermalGamma0010Stat                                 = FitObject("e","fitPureThermalGamma0010Stat","Photon",graphCombThermalGammaSpectrumStatErr0010,0.9,2.1,NULL,"QNRMEX0+");
        fileFinalResults << WriteParameterToFile(fitPureThermalGamma0010Stat)<< endl;
        fitFullThermalGamma0010Stat                                 = FitObject("qcd","fitFullDirGamma0010Stat","Photon",graphCombThermalGammaSpectrumStatErr0010,0.9,14,NULL,"QNRMEX0+");
        fileFinalResults << WriteParameterToFile(fitFullThermalGamma0010Stat)<< endl;
        graphRatioCombFitThermalGammaStatErr0010                    = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumStatErr0010->Clone();
        graphRatioCombFitThermalGammaStatErr0010                    = CalculateGraphErrRatioToFit(graphRatioCombFitThermalGammaStatErr0010, fitFullThermalGamma0010Stat);
    }
    if (graphCombDirGammaSpectrumSystErr0010){
        graphCombThermalGammaSpectrumSysErr0010                     = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0010, graphCombDirGammaSpectrumSystErr0010);
        graphCombThermalGammaSpectrumSysErr0010->Print();
        fitFullThermalGamma0010Sys                                  = FitObject("qcd","fitFullDirGamma0010Sys","Photon",graphCombThermalGammaSpectrumSysErr0010,0.9,14,NULL,"QNRMEX0+");
        fileFinalResults << WriteParameterToFile(fitFullThermalGamma0010Sys)<< endl;
        graphRatioCombFitThermalGammaSysErr0010                     = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumSysErr0010->Clone();
        graphRatioCombFitThermalGammaSysErr0010                     = CalculateGraphErrRatioToFit(graphRatioCombFitThermalGammaSysErr0010, fitFullThermalGamma0010Sys);

    }
    if (graphCombDirGammaSpectrumSumErr0010Ar){
        graphCombThermalGammaSpectrumSumErr0010Ar                   = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0010, graphCombDirGammaSpectrumSumErr0010Ar, newBinsComb, 20);
    }

    // Calculate RAA
    cout << __LINE__ << endl;
    cout << endl << "Calculating RAA" << endl;
    TGraphAsymmErrors* graphCombRAADirGammaStat0010                 = NULL;
    TGraphAsymmErrors* graphCombRAADirGammaSys0010                  = NULL;
    TGraphAsymmErrors* graphCombRAADirGammaSum0010Ar                = NULL;
    if (graphCombDirGammaSpectrumStatErr0010)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0010, graphCombDirGammaSpectrumStatErr0010, &graphCombRAADirGammaStat0010);
    if (graphCombDirGammaSpectrumSystErr0010)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0010, graphCombDirGammaSpectrumSystErr0010, &graphCombRAADirGammaSys0010);
    if (graphCombDirGammaSpectrumSumErr0010Ar)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0010, graphCombDirGammaSpectrumSumErr0010Ar, &graphCombRAADirGammaSum0010Ar, newBinsComb, 20);
    //if (graphCombRAADirGammaStat0010) graphCombRAADirGammaStat0010->Print();
    //if (graphCombRAADirGammaSys0010) graphCombRAADirGammaSys0010->Print();
    //if (graphCombRAADirGammaSum0010Ar) graphCombRAADirGammaSum0010Ar->Print();

    //*******************************************************************************************************************************************
    //*********************************************** Calculate direct photon spectrum 2040 *****************************************************
    //*******************************************************************************************************************************************
    //_______________________ copy inclusive photon spectra _____________________
    TH1D *histoCombDirGammaSpectrumErrSum2040                       = new TH1D("histoCombDirGammaSpectrumErrSum2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSys2040                       = new TH1D("histoCombDirGammaSpectrumErrSys2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSysA2040                      = new TH1D("histoCombDirGammaSpectrumErrSysA2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSysB2040                      = new TH1D("histoCombDirGammaSpectrumErrSysB2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSysC2040                      = new TH1D("histoCombDirGammaSpectrumErrSysC2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrStat2040                      = new TH1D("histoCombDirGammaSpectrumErrStat2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);

    //_______________________ get arrays of double ratio errors __________________
    Double_t *SystErrorsCombDR2040                                  = new Double_t[graphCombIncGammaStatErr2040->GetN()];
    Double_t *SystAErrorsCombDR2040                                 = new Double_t[graphCombIncGammaStatErr2040->GetN()];
    Double_t *SystBErrorsCombDR2040                                 = new Double_t[graphCombIncGammaStatErr2040->GetN()];
    Double_t *SystCErrorsCombDR2040                                 = new Double_t[graphCombIncGammaStatErr2040->GetN()];
    Double_t *sumErrorsCombDR2040                                   = new Double_t[graphCombIncGammaStatErr2040->GetN()];
    Double_t *StatErrorsCombDR2040                                  = new Double_t[graphCombIncGammaStatErr2040->GetN()];
    Double_t *xErrorsDR2040                                         = new Double_t[graphCombIncGammaStatErr2040->GetN()];
    graphCombIncGammaSysErr2040->Print();
    for (Int_t i = 0; i< graphCombDRPi0FitStatErr2040->GetN(); i++){
        SystErrorsCombDR2040[i]                                     = graphCombDRPi0FitSysErr2040->GetEYhigh()[i]/graphCombDRPi0FitSysErr2040->GetY()[i] *100;
        SystAErrorsCombDR2040[i]                                    = graphCombDRPi0FitSysAErr2040->GetEYhigh()[i]/graphCombDRPi0FitSysAErr2040->GetY()[i] *100;
        SystBErrorsCombDR2040[i]                                    = graphCombDRPi0FitSysBErr2040->GetEYhigh()[i]/graphCombDRPi0FitSysBErr2040->GetY()[i] *100;
        SystCErrorsCombDR2040[i]                                    = graphCombDRPi0FitSysCErr2040->GetEYhigh()[i]/graphCombDRPi0FitSysCErr2040->GetY()[i] *100;
        StatErrorsCombDR2040[i]                                     = graphCombDRPi0FitStatErr2040->GetEYhigh()[i]/graphCombDRPi0FitStatErr2040->GetY()[i] *100;
        sumErrorsCombDR2040[i]                                      = graphCombDRPi0FitSumErr2040->GetEYhigh()[i]/graphCombDRPi0FitSumErr2040->GetY()[i] *100;
        //cout << i << "\t" << graphCombDRPi0FitSysErr2040->GetY()[i] << "\t" << graphCombDRPi0FitSysErr2040->GetEYhigh()[i] << "\t" <<SystErrorsCombDR2040[i] << endl;
    }
    xErrorsDR2040                                                   = graphCombDRPi0FitStatErr2040->GetX();

    cout << __LINE__ << endl;
    //graphCombDRPi0FitSumErr2040->Print();


    //_______________________ copy inclusive photon spectra _____________________
    TH1D* histoCombErrorsForDRSum2040                               = new TH1D("histoCombErrorsForDRSum2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRStat2040                              = new TH1D("histoCombErrorsForDRStat2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSys2040                               = new TH1D("histoCombErrorsForDRSys2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSysA2040                              = new TH1D("histoCombErrorsForDRSysA2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSysB2040                              = new TH1D("histoCombErrorsForDRSysB2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSysC2040                              = new TH1D("histoCombErrorsForDRSysC2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);

    for(Int_t i = 1; i<graphCombDRPi0FitStatErr2040->GetN()+1;i++){
        //cout<< i << "\t"<<xErrorsDR2040[i-1]<<"  "<<histoCombErrorsForDRSum2040->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum2040->GetBinWidth(i) <<endl;
        Double_t binErrorSummed                                     = sumErrorsCombDR2040[i-1];
        Double_t binErrorSyst                                       = SystErrorsCombDR2040[i-1];
        Double_t binErrorSystA                                      = SystAErrorsCombDR2040[i-1];
        Double_t binErrorSystB                                      = SystBErrorsCombDR2040[i-1];
        Double_t binErrorSystC                                      = SystCErrorsCombDR2040[i-1];
        Double_t binErrorStat                                       = StatErrorsCombDR2040[i-1];
        Double_t DR                                                 = graphCombDRPi0FitStatErr2040->GetY()[i-1];

        //cout << DR << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
        histoCombErrorsForDRSum2040->SetBinContent(i,DR);
        histoCombErrorsForDRSys2040->SetBinContent(i,DR);
        histoCombErrorsForDRSysA2040->SetBinContent(i,DR);
        histoCombErrorsForDRSysB2040->SetBinContent(i,DR);
        histoCombErrorsForDRSysC2040->SetBinContent(i,DR);
        histoCombErrorsForDRStat2040->SetBinContent(i,DR);
        histoCombErrorsForDRSum2040->SetBinError(i,(binErrorSummed/100)*DR);
        histoCombErrorsForDRSys2040->SetBinError(i,(binErrorSyst/100)*DR);
        histoCombErrorsForDRSysA2040->SetBinError(i,(binErrorSystA/100)*DR);
        histoCombErrorsForDRSysB2040->SetBinError(i,(binErrorSystB/100)*DR);
        histoCombErrorsForDRSysC2040->SetBinError(i,(binErrorSystC/100)*DR);
        histoCombErrorsForDRStat2040->SetBinError(i,(binErrorStat/100)*DR);
    }

    for(Int_t i = 1; i<histoCombErrorsForDRSum2040->GetNbinsX()+1;i++){
        histoCombDirGammaSpectrumErrSum2040->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSys2040->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSysA2040->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSysB2040->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSysC2040->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrStat2040->SetBinContent(i+1,-1);

        histoCombDirGammaSpectrumErrSum2040->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSys2040->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSysA2040->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSysB2040->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSysC2040->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrStat2040->SetBinError(i+1,0);
    }

    // get the binning of the direct photons from the DR
    TH1D *histoCombDirGammaSpecSysErr2040                            = new TH1D(*histoCombErrorsForDRSys2040);
    TH1D *histoCombDirGammaSpecSysAErr2040                           = new TH1D(*histoCombErrorsForDRSysA2040);
    TH1D *histoCombDirGammaSpecSysBErr2040                           = new TH1D(*histoCombErrorsForDRSysB2040);
    TH1D *histoCombDirGammaSpecSysCErr2040                           = new TH1D(*histoCombErrorsForDRSysC2040);
    TH1D *histoCombDirGammaSpecStatErr2040                           = new TH1D(*histoCombErrorsForDRStat2040);
    TH1D *histoCombDirGammaSpecSumErr2040                            = new TH1D(*histoCombErrorsForDRSum2040);

    for(Int_t i = 1; i<graphCombDRPi0FitStatErr2040->GetN()+1; i++){
        // obtain common quantities
        Double_t Rgamma                 = histoCombErrorsForDRSys2040->GetBinContent(i);
        Double_t nIncGamma              = graphCombIncGammaStatErr2040->GetY()[i-1];

        // calculating Systematics graph
        Double_t errRgamma              = histoCombErrorsForDRSys2040->GetBinError(i);
        Double_t errNIncGam             = graphCombIncGammaSysErr2040->GetEYhigh()[i-1];
        Double_t q1                     = 1 - 1/ Rgamma;

        Double_t q1Error                = errRgamma/(Rgamma*Rgamma);
        Double_t content                = nIncGamma * ( 1 - 1/ Rgamma);
        Double_t error                  = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        Double_t errDR                  = content - error;
        histoCombDirGammaSpecSysErr2040->SetBinError(i, error);
        histoCombDirGammaSpecSysErr2040->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSys2040->SetBinContent(i, errDR);

        // calculating Systematics A graph
        errRgamma                       = histoCombErrorsForDRSysA2040->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSysAErr2040->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSysAErr2040->SetBinError(i, error);
        histoCombDirGammaSpecSysAErr2040->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSysA2040->SetBinContent(i, errDR);

        // calculating Systematics B graph
        errRgamma                       = histoCombErrorsForDRSysB2040->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSysBErr2040->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSysBErr2040->SetBinError(i, error);
        histoCombDirGammaSpecSysBErr2040->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSysB2040->SetBinContent(i, errDR);

        // calculating Systematics C graph
        errRgamma                       = histoCombErrorsForDRSysC2040->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSysCErr2040->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSysCErr2040->SetBinError(i, error);
        histoCombDirGammaSpecSysCErr2040->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSysC2040->SetBinContent(i, errDR);

        // calculating Stat graphs
        errRgamma                       = histoCombErrorsForDRStat2040->GetBinError(i);
        errNIncGam                      = graphCombIncGammaStatErr2040->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecStatErr2040->SetBinError(i, error);
        histoCombDirGammaSpecStatErr2040->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrStat2040->SetBinContent(i, errDR);

        // calculating summed error graphs
        errRgamma                       = histoCombErrorsForDRSum2040->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSumErr2040->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSumErr2040->SetBinError(i, error);
        histoCombDirGammaSpecSumErr2040->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSum2040->SetBinContent(i, errDR);
    }
    // purely calculating points based on all Systematic errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystErr2040         = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSys2040,histoCombDirGammaSpecStatErr2040,0,0.5);
    if(graphCombDirGammaSpectrumSystErr2040)graphCombDirGammaSpectrumSystErr2040->SetName("graphCombDirGammaSpectrumSystErr2040");
    //if(graphCombDirGammaSpectrumSystErr2040)graphCombDirGammaSpectrumSystErr2040->Print();
    //if(graphCombDirGammaSpectrumSystErr2040)cout << "graph has been found" << endl;
    // purely calculating points based on all Systematic errors A
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystAErr2040        = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysA2040,histoCombDirGammaSpecStatErr2040,0,0.5);
    if(graphCombDirGammaSpectrumSystAErr2040)graphCombDirGammaSpectrumSystAErr2040->SetName("graphCombDirGammaSpectrumSystAErr2040");
    //if(graphCombDirGammaSpectrumSystAErr2040)graphCombDirGammaSpectrumSystAErr2040->Print();
    // purely calculating points based on all Systematic errors B
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystBErr2040        = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysB2040,histoCombDirGammaSpecStatErr2040,0,0.5);
    if(graphCombDirGammaSpectrumSystBErr2040)graphCombDirGammaSpectrumSystBErr2040->SetName("graphCombDirGammaSpectrumSystBErr2040");
    //if(graphCombDirGammaSpectrumSystBErr2040)graphCombDirGammaSpectrumSystBErr2040->Print();
    // purely calculating points based on all Systematic errors C
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystCErr2040        = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysC2040,histoCombDirGammaSpecStatErr2040,0,0.5);
    if(graphCombDirGammaSpectrumSystCErr2040)graphCombDirGammaSpectrumSystCErr2040->SetName("graphCombDirGammaSpectrumSystCErr2040");
    //if(graphCombDirGammaSpectrumSystCErr2040)graphCombDirGammaSpectrumSystCErr2040->Print();

    // purely calculating points based on Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumStatErr2040         = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrStat2040,histoCombDirGammaSpecStatErr2040,0,0.5);
    if(graphCombDirGammaSpectrumStatErr2040)graphCombDirGammaSpectrumStatErr2040->SetName("graphCombDirGammaSpectrumStatErr2040");
    //if(graphCombDirGammaSpectrumStatErr2040)graphCombDirGammaSpectrumStatErr2040->Print();
    // purely calculating points based on all Systematic + Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr2040          = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum2040,histoCombDirGammaSpecStatErr2040,0,0.5);
    if(graphCombDirGammaSpectrumSumErr2040)graphCombDirGammaSpectrumSumErr2040->SetName("graphCombDirGammaSpectrumSumErr2040");
    //if(graphCombDirGammaSpectrumSumErr2040)graphCombDirGammaSpectrumSumErr2040->Print();
    // calculate arrows for points with 0, error summed
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr2040Ar        = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum2040,histoCombDirGammaSpecStatErr2040,5,0.5);
    if(graphCombDirGammaSpectrumSumErr2040Ar)graphCombDirGammaSpectrumSumErr2040Ar->SetName("graphCombDirGammaSpectrumSumErr2040Ar");
    //if(graphCombDirGammaSpectrumSumErr2040Ar)graphCombDirGammaSpectrumSumErr2040Ar->Print();
    // calculate points below confidence level summed errors with arrows
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr2040ArConfi   = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum2040,histoCombDirGammaSpecStatErr2040,7,0.5);
    if(graphCombDirGammaSpectrumSumErr2040ArConfi)graphCombDirGammaSpectrumSumErr2040ArConfi->SetName("graphCombDirGammaSpectrumSumErr2040Ar");
    //if(graphCombDirGammaSpectrumSumErr2040ArConfi)graphCombDirGammaSpectrumSumErr2040ArConfi->Print();

    TF1* fitThermalGamma2040Sum                                     = FitObject("e","fitThermalGamma2040Sum","Photon",histoCombDirGammaSpectrumErrSum2040,1.1,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma2040Sum)<< endl;
    TF1* fitThermalGamma2040Sum23                                   = FitObject("e","fitThermalGamma2040Sum23","Photon",histoCombDirGammaSpectrumErrSum2040,1.1,2.3,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma2040Sum23)<< endl;
    TF1* fitThermalGamma2040Stat                                    = FitObject("e","fitThermalGamma2040Stat","Photon",graphCombDirGammaSpectrumStatErr2040,1.1,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma2040Stat)<< endl;
    TF1* fitThermalGamma2040Sys                                     = FitObject("e","fitThermalGamma2040Sys","Photon",graphCombDirGammaSpectrumSystErr2040,1.1,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma2040Sys)<< endl;
    TF1* fitThermalGamma2040SysA                                    = FitObject("e","fitThermalGamma2040SysA","Photon",graphCombDirGammaSpectrumSystAErr2040,1.1,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma2040SysA)<< endl;
    TF1* fitThermalGamma2040SysB                                    = FitObject("e","fitThermalGamma2040SysB","Photon",graphCombDirGammaSpectrumSystBErr2040,1.1,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma2040SysB)<< endl;
    TF1* fitThermalGamma2040SysC                                    = FitObject("e","fitThermalGamma2040SysC","Photon",graphCombDirGammaSpectrumSystCErr2040,1.1,2.1,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitThermalGamma2040SysC)<< endl;

    TF1* fitFullDirGamma2040Sys                                     = FitObject("qcd","fitFullDirGamma2040Sys","Photon",graphCombDirGammaSpectrumSystErr2040,0.9,14,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitFullDirGamma2040Sys)<< endl;
    TF1* fitFullDirGamma2040Stat                                    = FitObject("qcd","fitFullDirGamma2040Stat","Photon",graphCombDirGammaSpectrumStatErr2040,0.9,14,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitFullDirGamma2040Stat)<< endl;

    TGraphAsymmErrors* graphRatioCombFitDirGammaStatErr2040         = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr2040->Clone();
    TGraphAsymmErrors* graphRatioCombFitDirGammaSysErr2040          = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr2040->Clone();
    graphRatioCombFitDirGammaStatErr2040                            = CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaStatErr2040, fitFullDirGamma2040Stat);
    graphRatioCombFitDirGammaSysErr2040                             = CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaSysErr2040, fitFullDirGamma2040Sys);

    TGraphAsymmErrors* graphCombThermalGammaSpectrumStatErr2040     = SubtractPromptPhotonsViaFit(     fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumStatErr2040);
    TGraphAsymmErrors* graphCombThermalGammaSpectrumSysErr2040      = SubtractPromptPhotonsViaFit(     fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumSystErr2040);
    TGraphAsymmErrors* graphCombThermalGammaSpectrumSumErr2040Ar    = SubtractPromptPhotonsViaFit(     fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumSumErr2040Ar, newBinsComb, 20);

    cout << endl << "Calculating RAA" << endl;
    TGraphAsymmErrors* graphCombRAADirGammaStat2040                 = NULL;
    TGraphAsymmErrors* graphCombRAADirGammaSys2040                  = NULL;
    TGraphAsymmErrors* graphCombRAADirGammaSum2040Ar                = NULL;
    if (graphCombDirGammaSpectrumStatErr2040)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumStatErr2040, &graphCombRAADirGammaStat2040);
    if (graphCombDirGammaSpectrumSystErr2040)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumSystErr2040, &graphCombRAADirGammaSys2040);
    if (graphCombDirGammaSpectrumSumErr2040Ar)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumSumErr2040Ar, &graphCombRAADirGammaSum2040Ar, newBinsComb, 20);
    //if (graphCombRAADirGammaStat2040) graphCombRAADirGammaStat2040->Print();
    //if (graphCombRAADirGammaSys2040) graphCombRAADirGammaSys2040->Print();
    //if (graphCombRAADirGammaSum2040Ar) graphCombRAADirGammaSum2040Ar->Print();

    //*******************************************************************************************************************************************
    //*********************************************** Calculate direct photon spectrum 2050 *****************************************************
    //*******************************************************************************************************************************************
    //_______________________ copy inclusive photon spectra _____________________
    TH1D *histoCombDirGammaSpectrumErrSum2050                       = new TH1D("histoCombDirGammaSpectrumErrSum2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSys2050                       = new TH1D("histoCombDirGammaSpectrumErrSys2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSysA2050                      = new TH1D("histoCombDirGammaSpectrumErrSysA2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSysB2050                      = new TH1D("histoCombDirGammaSpectrumErrSysB2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrSysC2050                      = new TH1D("histoCombDirGammaSpectrumErrSysC2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);
    TH1D *histoCombDirGammaSpectrumErrStat2050                      = new TH1D("histoCombDirGammaSpectrumErrStat2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);

    //_______________________ get arrays of double ratio errors __________________
    Double_t *SystErrorsCombDR2050                                  = new Double_t[graphCombIncGammaStatErr2050->GetN()];
    Double_t *SystAErrorsCombDR2050                                 = new Double_t[graphCombIncGammaStatErr2050->GetN()];
    Double_t *SystBErrorsCombDR2050                                 = new Double_t[graphCombIncGammaStatErr2050->GetN()];
    Double_t *SystCErrorsCombDR2050                                 = new Double_t[graphCombIncGammaStatErr2050->GetN()];
    Double_t *sumErrorsCombDR2050                                   = new Double_t[graphCombIncGammaStatErr2050->GetN()];
    Double_t *StatErrorsCombDR2050                                  = new Double_t[graphCombIncGammaStatErr2050->GetN()];
    Double_t *xErrorsDR2050                                         = new Double_t[graphCombIncGammaStatErr2050->GetN()];
    graphCombIncGammaSysErr2050->Print();
    for (Int_t i = 0; i< graphCombDRPi0FitStatErr2050->GetN(); i++){
        SystErrorsCombDR2050[i]                                     = graphCombDRPi0FitSysErr2050->GetEYhigh()[i]/graphCombDRPi0FitSysErr2050->GetY()[i] *100;
        SystAErrorsCombDR2050[i]                                    = graphCombDRPi0FitSysAErr2050->GetEYhigh()[i]/graphCombDRPi0FitSysAErr2050->GetY()[i] *100;
        SystBErrorsCombDR2050[i]                                    = graphCombDRPi0FitSysBErr2050->GetEYhigh()[i]/graphCombDRPi0FitSysBErr2050->GetY()[i] *100;
        SystCErrorsCombDR2050[i]                                    = graphCombDRPi0FitSysCErr2050->GetEYhigh()[i]/graphCombDRPi0FitSysCErr2050->GetY()[i] *100;
        StatErrorsCombDR2050[i]                                     = graphCombDRPi0FitStatErr2050->GetEYhigh()[i]/graphCombDRPi0FitStatErr2050->GetY()[i] *100;
        sumErrorsCombDR2050[i]                                      = graphCombDRPi0FitSumErr2050->GetEYhigh()[i]/graphCombDRPi0FitSumErr2050->GetY()[i] *100;
        //cout << i << "\t" << graphCombDRPi0FitSysErr2050->GetY()[i] << "\t" << graphCombDRPi0FitSysErr2050->GetEYhigh()[i] << "\t" <<SystErrorsCombDR2050[i] << endl;
    }
    xErrorsDR2050                                                   = graphCombDRPi0FitStatErr2050->GetX();

    cout << __LINE__ << endl;
    //graphCombDRPi0FitSumErr2050->Print();


    //_______________________ copy inclusive photon spectra _____________________
    TH1D* histoCombErrorsForDRSum2050                               = new TH1D("histoCombErrorsForDRSum2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRStat2050                              = new TH1D("histoCombErrorsForDRStat2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSys2050                               = new TH1D("histoCombErrorsForDRSys2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSysA2050                              = new TH1D("histoCombErrorsForDRSysA2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSysB2050                              = new TH1D("histoCombErrorsForDRSysB2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);
    TH1D* histoCombErrorsForDRSysC2050                              = new TH1D("histoCombErrorsForDRSysC2050","",graphCombDRPi0FitStatErr2050->GetN(),xArrayCombined);

    for(Int_t i = 1; i<graphCombDRPi0FitStatErr2050->GetN()+1;i++){
        //cout<< i << "\t"<<xErrorsDR2050[i-1]<<"  "<<histoCombErrorsForDRSum2050->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum2050->GetBinWidth(i) <<endl;
        Double_t binErrorSummed                                     = sumErrorsCombDR2050[i-1];
        Double_t binErrorSyst                                       = SystErrorsCombDR2050[i-1];
        Double_t binErrorSystA                                      = SystAErrorsCombDR2050[i-1];
        Double_t binErrorSystB                                      = SystBErrorsCombDR2050[i-1];
        Double_t binErrorSystC                                      = SystCErrorsCombDR2050[i-1];
        Double_t binErrorStat                                       = StatErrorsCombDR2050[i-1];
        Double_t DR                                                 = graphCombDRPi0FitStatErr2050->GetY()[i-1];

        //cout << DR << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
        histoCombErrorsForDRSum2050->SetBinContent(i,DR);
        histoCombErrorsForDRSys2050->SetBinContent(i,DR);
        histoCombErrorsForDRSysA2050->SetBinContent(i,DR);
        histoCombErrorsForDRSysB2050->SetBinContent(i,DR);
        histoCombErrorsForDRSysC2050->SetBinContent(i,DR);
        histoCombErrorsForDRStat2050->SetBinContent(i,DR);
        histoCombErrorsForDRSum2050->SetBinError(i,(binErrorSummed/100)*DR);
        histoCombErrorsForDRSys2050->SetBinError(i,(binErrorSyst/100)*DR);
        histoCombErrorsForDRSysA2050->SetBinError(i,(binErrorSystA/100)*DR);
        histoCombErrorsForDRSysB2050->SetBinError(i,(binErrorSystB/100)*DR);
        histoCombErrorsForDRSysC2050->SetBinError(i,(binErrorSystC/100)*DR);
        histoCombErrorsForDRStat2050->SetBinError(i,(binErrorStat/100)*DR);
    }

    for(Int_t i = 1; i<histoCombErrorsForDRSum2050->GetNbinsX()+1;i++){
        histoCombDirGammaSpectrumErrSum2050->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSys2050->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSysA2050->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSysB2050->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSysC2050->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrStat2050->SetBinContent(i+1,-1);

        histoCombDirGammaSpectrumErrSum2050->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSys2050->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSysA2050->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSysB2050->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSysC2050->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrStat2050->SetBinError(i+1,0);
    }

    // get the binning of the direct photons from the DR
    TH1D *histoCombDirGammaSpecSysErr2050                            = new TH1D(*histoCombErrorsForDRSys2050);
    TH1D *histoCombDirGammaSpecSysAErr2050                           = new TH1D(*histoCombErrorsForDRSysA2050);
    TH1D *histoCombDirGammaSpecSysBErr2050                           = new TH1D(*histoCombErrorsForDRSysB2050);
    TH1D *histoCombDirGammaSpecSysCErr2050                           = new TH1D(*histoCombErrorsForDRSysC2050);
    TH1D *histoCombDirGammaSpecStatErr2050                           = new TH1D(*histoCombErrorsForDRStat2050);
    TH1D *histoCombDirGammaSpecSumErr2050                            = new TH1D(*histoCombErrorsForDRSum2050);

    for(Int_t i = 1; i<graphCombDRPi0FitStatErr2050->GetN()+1; i++){
        // obtain common quantities
        Double_t Rgamma                 = histoCombErrorsForDRSys2050->GetBinContent(i);
        Double_t nIncGamma              = graphCombIncGammaStatErr2050->GetY()[i-1];

        // calculating Systematics graph
        Double_t errRgamma              = histoCombErrorsForDRSys2050->GetBinError(i);
        Double_t errNIncGam             = graphCombIncGammaSysErr2050->GetEYhigh()[i-1];
        Double_t q1                     = 1 - 1/ Rgamma;

        Double_t q1Error                = errRgamma/(Rgamma*Rgamma);
        Double_t content                = nIncGamma * ( 1 - 1/ Rgamma);
        Double_t error                  = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        Double_t errDR                  = content - error;
        histoCombDirGammaSpecSysErr2050->SetBinError(i, error);
        histoCombDirGammaSpecSysErr2050->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSys2050->SetBinContent(i, errDR);

        // calculating Systematics A graph
        errRgamma                       = histoCombErrorsForDRSysA2050->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSysAErr2050->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSysAErr2050->SetBinError(i, error);
        histoCombDirGammaSpecSysAErr2050->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSysA2050->SetBinContent(i, errDR);

        // calculating Systematics B graph
        errRgamma                       = histoCombErrorsForDRSysB2050->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSysBErr2050->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSysBErr2050->SetBinError(i, error);
        histoCombDirGammaSpecSysBErr2050->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSysB2050->SetBinContent(i, errDR);

        // calculating Systematics C graph
        errRgamma                       = histoCombErrorsForDRSysC2050->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSysCErr2050->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSysCErr2050->SetBinError(i, error);
        histoCombDirGammaSpecSysCErr2050->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSysC2050->SetBinContent(i, errDR);

        // calculating Stat graphs
        errRgamma                       = histoCombErrorsForDRStat2050->GetBinError(i);
        errNIncGam                      = graphCombIncGammaStatErr2050->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecStatErr2050->SetBinError(i, error);
        histoCombDirGammaSpecStatErr2050->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrStat2050->SetBinContent(i, errDR);

        // calculating summed error graphs
        errRgamma                       = histoCombErrorsForDRSum2050->GetBinError(i);
        errNIncGam                      = graphCombIncGammaSumErr2050->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSumErr2050->SetBinError(i, error);
        histoCombDirGammaSpecSumErr2050->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSum2050->SetBinContent(i, errDR);
    }
    // purely calculating points based on all Systematic errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystErr2050     = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSys2050,histoCombDirGammaSpecStatErr2050,0,0.5);
    if(graphCombDirGammaSpectrumSystErr2050)graphCombDirGammaSpectrumSystErr2050->SetName("graphCombDirGammaSpectrumSystErr2050");
    //if(graphCombDirGammaSpectrumSystErr2050)graphCombDirGammaSpectrumSystErr2050->Print();
    //if(graphCombDirGammaSpectrumSystErr2050)cout << "graph has been found" << endl;
    // purely calculating points based on all Systematic errors A
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystAErr2050    = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysA2050,histoCombDirGammaSpecStatErr2050,0,0.5);
    if(graphCombDirGammaSpectrumSystAErr2050)graphCombDirGammaSpectrumSystAErr2050->SetName("graphCombDirGammaSpectrumSystAErr2050");
    //if(graphCombDirGammaSpectrumSystAErr2050)graphCombDirGammaSpectrumSystAErr2050->Print();
    // purely calculating points based on all Systematic errors B
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystBErr2050    = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysB2050,histoCombDirGammaSpecStatErr2050,0,0.5);
    if(graphCombDirGammaSpectrumSystBErr2050)graphCombDirGammaSpectrumSystBErr2050->SetName("graphCombDirGammaSpectrumSystBErr2050");
    //if(graphCombDirGammaSpectrumSystBErr2050)graphCombDirGammaSpectrumSystBErr2050->Print();
    // purely calculating points based on all Systematic errors C
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystCErr2050    = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysC2050,histoCombDirGammaSpecStatErr2050,0,0.5);
    if(graphCombDirGammaSpectrumSystCErr2050)graphCombDirGammaSpectrumSystCErr2050->SetName("graphCombDirGammaSpectrumSystCErr2050");
    //if(graphCombDirGammaSpectrumSystCErr2050)graphCombDirGammaSpectrumSystCErr2050->Print();

    // purely calculating points based on Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumStatErr2050     = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrStat2050,histoCombDirGammaSpecStatErr2050,0,0.5);
    if(graphCombDirGammaSpectrumStatErr2050)graphCombDirGammaSpectrumStatErr2050->SetName("graphCombDirGammaSpectrumStatErr2050");
    //if(graphCombDirGammaSpectrumStatErr2050)graphCombDirGammaSpectrumStatErr2050->Print();
    // purely calculating points based on all Systematic + Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr2050      = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum2050,histoCombDirGammaSpecStatErr2050,0,0.5);
    if(graphCombDirGammaSpectrumSumErr2050)graphCombDirGammaSpectrumSumErr2050->SetName("graphCombDirGammaSpectrumSumErr2050");
    //if(graphCombDirGammaSpectrumSumErr2050)graphCombDirGammaSpectrumSumErr2050->Print();
    // calculate arrows for points with 0, error summed
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr2050Ar    = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum2050,histoCombDirGammaSpecStatErr2050,5,0.5);
    if(graphCombDirGammaSpectrumSumErr2050Ar)graphCombDirGammaSpectrumSumErr2050Ar->SetName("graphCombDirGammaSpectrumSumErr2050Ar");
    if(graphCombDirGammaSpectrumSumErr2050Ar)graphCombDirGammaSpectrumSumErr2050Ar->Print();

    TF1* fitFullDirGamma2050Sys                                 = FitObject("qcd","fitFullDirGamma2050Sys","Photon",graphCombDirGammaSpectrumSystErr2050,0.9,14,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitFullDirGamma2050Sys)<< endl;
    TF1* fitFullDirGamma2050Stat                                = FitObject("qcd","fitFullDirGamma2050Stat","Photon",graphCombDirGammaSpectrumStatErr2050,0.9,14,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitFullDirGamma2050Stat)<< endl;

    TGraphAsymmErrors* graphRatioCombFitDirGammaStatErr2050     = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr2050->Clone();
    TGraphAsymmErrors* graphRatioCombFitDirGammaSysErr2050      = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr2050->Clone();
    graphRatioCombFitDirGammaStatErr2050                        = CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaStatErr2050, fitFullDirGamma2050Stat);
    graphRatioCombFitDirGammaSysErr2050                         = CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaSysErr2050, fitFullDirGamma2050Sys);


    cout << endl << "Calculating RAA" << endl;
    TGraphAsymmErrors* graphCombRAADirGammaStat2050             = NULL;
    TGraphAsymmErrors* graphCombRAADirGammaSys2050              = NULL;
    TGraphAsymmErrors* graphCombRAADirGammaSum2050Ar            = NULL;
//     if (graphCombDirGammaSpectrumStatErr2050)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2050, graphCombDirGammaSpectrumStatErr2050, &graphCombRAADirGammaStat2050);
//     if (graphCombDirGammaSpectrumSystErr2050)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2050, graphCombDirGammaSpectrumSystErr2050, &graphCombRAADirGammaSys2050);
//     if (graphCombDirGammaSpectrumSumErr2050Ar)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2050, graphCombDirGammaSpectrumSumErr2050Ar, &graphCombRAADirGammaSum2050Ar, newBinsComb, 20);
    //if (graphCombRAADirGammaStat2050) graphCombRAADirGammaStat2050->Print();
    //if (graphCombRAADirGammaSys2050) graphCombRAADirGammaSys2050->Print();
    //if (graphCombRAADirGammaSum2050Ar) graphCombRAADirGammaSum2050Ar->Print();


    //*******************************************************************************************************************************************
    //********************************************** DR plot with combined measurement **********************************************************
    //*******************************************************************************************************************************************
    cout << "Plotting DR combined with theory at " << __LINE__ << endl;
    canvasRatioIndDR->cd();

    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();

    TGraphAsymmErrors* graphCombDRPi0FitStatErr0010Plot = (TGraphAsymmErrors*)graphCombDRPi0FitStatErr0010->Clone("graphCombDRPi0FitStatErr0010Plot");
    ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatErr0010Plot);
    TGraphAsymmErrors* graphCombDRPi0FitStatErr2040Plot = (TGraphAsymmErrors*)graphCombDRPi0FitStatErr2040->Clone("graphCombDRPi0FitStatErr2040Plot");
    ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatErr2040Plot);
    TGraphAsymmErrors* graphCombDRPi0FitStatErr2050Plot = (TGraphAsymmErrors*)graphCombDRPi0FitStatErr2050->Clone("graphCombDRPi0FitStatErr2050Plot");
    ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatErr2050Plot);
    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        SetStyleGammaNLOTGraphWithBand( graphTheoryEPS09DR0010, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryEPS09DR0010->Draw("p3lsame");
        SetStyleGammaNLOTGraphWithBand( graphTheoryCT10DR0010, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
        graphTheoryCT10DR0010->Draw("p3lsame");
        graphTheoryNLODR0010->Draw("p3lsame");

        graphCombDRPi0FitSysErr0010->Draw("E2same");
        graphCombDRPi0FitStatErr0010Plot->Draw("p,E1Z,same");

        TLegend* legendDRCombNLO0010 = new TLegend(0.12,0.82-1.02*0.85*textsizeLabelsPad1*5,0.42,0.82);
        legendDRCombNLO0010->SetFillStyle(0);
        legendDRCombNLO0010->SetFillColor(0);
        legendDRCombNLO0010->SetLineColor(0);
        legendDRCombNLO0010->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRCombNLO0010->SetMargin(0.23);
        legendDRCombNLO0010->SetTextFont(42);
        legendDRCombNLO0010->AddEntry(graphCombDRPi0FitSysErr0010,"ALICE","pf");
        legendDRCombNLO0010->AddEntry(graphTheoryNLODR0010,"NLO pQCD PDF: CTEQ6M5 FF: GRV ","l");
        legendDRCombNLO0010->AddEntry(graphTheoryCT10DR0010,"JETPHOX PDF: CT10, FF: BFG2","f");
        legendDRCombNLO0010->AddEntry(graphTheoryEPS09DR0010,"JETPHOX nPDF: EPS09, FF: BFG2","f");
        legendDRCombNLO0010->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
        legendDRCombNLO0010->Draw();

        labelDRCent0010->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->DrawCopy();
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);

        SetStyleGammaNLOTGraphWithBand( graphTheoryEPS09DR2040, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryEPS09DR2040->Draw("p3,l,same");
        SetStyleGammaNLOTGraphWithBand( graphTheoryCT10DR2040, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
        graphTheoryCT10DR2040->Draw("p3,l,same");
        graphTheoryNLODR2040->Draw("p3,l,same");

        graphCombDRPi0FitSysErr2040->Draw("E2same");
        graphCombDRPi0FitStatErr2040Plot->Draw("p,E1Z,same");

        TLegend* legendDRCombNLO2040 = new TLegend(0.12,0.85-1.02*0.85*textsizeLabelsPad2*5,0.42,0.85);
        legendDRCombNLO2040->SetFillStyle(0);
        legendDRCombNLO2040->SetFillColor(0);
        legendDRCombNLO2040->SetLineColor(0);
        legendDRCombNLO2040->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRCombNLO2040->SetMargin(0.23);
        legendDRCombNLO2040->SetTextFont(42);
        legendDRCombNLO2040->AddEntry(graphCombDRPi0FitSysErr2040,"ALICE","pf");
        legendDRCombNLO2040->AddEntry(graphTheoryNLODR0010,"NLO pQCD PDF: CTEQ6M5 FF: GRV ","l");
        legendDRCombNLO2040->AddEntry(graphTheoryCT10DR0010,"JETPHOX PDF: CT10, FF: BFG2","f");
        legendDRCombNLO2040->AddEntry(graphTheoryEPS09DR0010,"JETPHOX nPDF: EPS09, FF: BFG2","f");
        legendDRCombNLO2040->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
        legendDRCombNLO2040->Draw();

        labelDRCent2040->Draw();

    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        SetStyleGammaNLOTGraphWithBand( graphTheoryEPS09DR2050, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryEPS09DR2050->Draw("p3lsame");
        SetStyleGammaNLOTGraphWithBand( graphTheoryCT10DR2050, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
        graphTheoryCT10DR2050->Draw("p3lsame");
        graphTheoryNLODR2050->Draw("p3lsame");

        graphCombDRPi0FitSysErr2050->Draw("E2same");
        graphCombDRPi0FitStatErr2050Plot->Draw("p,E1Z,same");

        TLegend* legendDRCombNLO2050 = new TLegend(0.12,0.87-1.02*0.85*textsizeLabelsPad3*5,0.42,0.87);
        legendDRCombNLO2050->SetFillStyle(0);
        legendDRCombNLO2050->SetFillColor(0);
        legendDRCombNLO2050->SetLineColor(0);
        legendDRCombNLO2050->SetTextSize(0.85*textsizeLabelsPad3);
        legendDRCombNLO2050->SetMargin(0.23);
        legendDRCombNLO2050->SetTextFont(42);
        legendDRCombNLO2050->AddEntry(graphCombDRPi0FitSysErr2050,"ALICE","pf");
        legendDRCombNLO2050->AddEntry(graphTheoryNLODR0010,"NLO pQCD PDF: CTEQ6M5 FF: GRV ","l");
        legendDRCombNLO2050->AddEntry(graphTheoryCT10DR0010,"JETPHOX PDF: CT10, FF: BFG2","f");
        legendDRCombNLO2050->AddEntry(graphTheoryEPS09DR0010,"JETPHOX nPDF: EPS09, FF: BFG2","f");
        legendDRCombNLO2050->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
        legendDRCombNLO2050->Draw();

        labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_combMeasurement_incNLO.%s", outputDir.Data(), suffix.Data()));


    //*******************************************************************************************************************************************
    //********************************************** DR plot with combined measurement inc NLO **************************************************
    //*******************************************************************************************************************************************
    canvasRatioIndDR->cd();
    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();
    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRPi0FitSysErr0010->Draw("E2same");
        graphCombDRPi0FitStatErr0010Plot->Draw("p,E1Z,same");

        labelDRCent0010->Draw();

        TLatex *labelALICECent0010 = new TLatex(0.87,0.85,"ALICE");
        SetStyleTLatex( labelALICECent0010, 0.85*textsizeLabelsPad1,4);
        labelALICECent0010->Draw();
    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRPi0FitSysErr2040->Draw("E2same");
        graphCombDRPi0FitStatErr2040Plot->Draw("p,E1Z,same");

        labelDRCent2040->Draw();

        TLatex *labelALICECent2040 = new TLatex(0.87,0.88,"ALICE");
        SetStyleTLatex( labelALICECent2040, 0.85*textsizeLabelsPad2,4);
        labelALICECent2040->Draw();
    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRPi0FitSysErr2050->Draw("E2same");
        graphCombDRPi0FitStatErr2050Plot->Draw("p,E1Z,same");

        labelDRCent2050->Draw();

        TLatex *labelALICECent2050 = new TLatex(0.87,0.9,"ALICE");
        SetStyleTLatex( labelALICECent2050, 0.85*textsizeLabelsPad3,4);
        labelALICECent2050->Draw();


    canvasRatioIndDR->SaveAs(Form("%s/DR_combMeasurement.%s", outputDir.Data(), suffix.Data()));

    canvasRatioIndDR->cd();
    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();
    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRPi0FitSysErr0010->Draw("E2same");
        graphCombDRPi0FitStatErr0010Plot->Draw("p,E1Z,same");

        graphPubCombDoubleRatioStat_0020->Draw("p,E1Z,same");
        graphPubCombDoubleRatioSyst_0020->Draw("E2same");

        TLegend* legendDRCombWithPublished_0010 = new TLegend(0.12,0.85-1.02*0.85*textsizeLabelsPad2*2,0.42,0.85);
        legendDRCombWithPublished_0010->SetFillStyle(0);
        legendDRCombWithPublished_0010->SetFillColor(0);
        legendDRCombWithPublished_0010->SetLineColor(0);
        legendDRCombWithPublished_0010->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRCombWithPublished_0010->SetMargin(0.23);
        legendDRCombWithPublished_0010->SetTextFont(42);
        legendDRCombWithPublished_0010->AddEntry(graphCombDRPi0FitSysErr0010,"0-10% combined PCM (2011) + PHOS (2010)","pf");
        legendDRCombWithPublished_0010->AddEntry(graphPubCombDoubleRatioSyst_0020,"0-20% published combined","pf");
        legendDRCombWithPublished_0010->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRPi0FitSysErr2040->Draw("E2same");
        graphCombDRPi0FitStatErr2040Plot->Draw("p,E1Z,same");

        graphPubCombDoubleRatioStat_2040->Draw("p,E1Z,same");
        graphPubCombDoubleRatioSyst_2040->Draw("E2same");

        TLegend* legendDRCombWithPublished_2040 = new TLegend(0.12,0.85-1.02*0.85*textsizeLabelsPad2*2,0.42,0.85);
        legendDRCombWithPublished_2040->SetFillStyle(0);
        legendDRCombWithPublished_2040->SetFillColor(0);
        legendDRCombWithPublished_2040->SetLineColor(0);
        legendDRCombWithPublished_2040->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRCombWithPublished_2040->SetMargin(0.23);
        legendDRCombWithPublished_2040->SetTextFont(42);
        legendDRCombWithPublished_2040->AddEntry(graphCombDRPi0FitSysErr2040,"0-10% combined PCM (2011) + PHOS (2010)","pf");
        legendDRCombWithPublished_2040->AddEntry(graphPubCombDoubleRatioSyst_2040,"0-20% published combined","pf");
        legendDRCombWithPublished_2040->Draw();
    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRPi0FitSysErr2050->Draw("E2same");
        graphCombDRPi0FitStatErr2050Plot->Draw("p,E1Z,same");

        labelDRCent2050->Draw();


    canvasRatioIndDR->SaveAs(Form("%s/DR_combMeasurementWithPublished.%s", outputDir.Data(), suffix.Data()));



    //*******************************************************************************************************************************************
    //********************************************** DR plot with combined measurement diff Syst Err ********************************************
    //*******************************************************************************************************************************************
    ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatSysAErr0010);
    ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatSysAErr2040);
    ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatSysAErr2050);

    canvasRatioIndDR->cd();
    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();
    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");

        TBox* boxSysCErrComb0010 = CreateBoxConv(colorComb0010Box, 0.74, 1.-(SysCCombDRPi0Fit0010) , 0.83, 1.+(SysCCombDRPi0Fit0010));
        boxSysCErrComb0010->Draw();
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysCErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010Box , colorComb0010Box, widthLinesBoxes,  kTRUE, colorComb0010Box);
        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysBErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatSysAErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        graphCombDRPi0FitSysBErr0010->Draw("E2same");
        graphCombDRPi0FitStatSysAErr0010->Draw("p,E1Z,same");

        TLegend* legendDRCombDiffErrRepPad1 = new TLegend(0.15,0.82-1.1*0.85*textsizeLabelsPad1*3,0.35,0.82);
        legendDRCombDiffErrRepPad1->SetFillStyle(0);
        legendDRCombDiffErrRepPad1->SetFillColor(0);
        legendDRCombDiffErrRepPad1->SetLineColor(0);
        legendDRCombDiffErrRepPad1->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRCombDiffErrRepPad1->SetMargin(0.25);
        legendDRCombDiffErrRepPad1->SetTextFont(42);
        legendDRCombDiffErrRepPad1->AddEntry(graphCombDRPi0FitStatSysAErr0010,"Stat. #oplus Syst. A","pe");
        legendDRCombDiffErrRepPad1->AddEntry(graphCombDRPi0FitSysBErr0010,"Syst. B","f");
        legendDRCombDiffErrRepPad1->AddEntry(graphCombDRPi0FitSysCErr0010,"Syst. C","f");
        legendDRCombDiffErrRepPad1->Draw();

        labelDRCent0010->Draw();
        labelALICECent0010->Draw();
    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");

        TBox* boxSysCErrComb2040 = CreateBoxConv(colorComb2040Box, 0.74, 1.-(SysCCombDRPi0Fit2040) , 0.83, 1.+(SysCCombDRPi0Fit2040));
        boxSysCErrComb2040->Draw();
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysCErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040Box , colorComb2040Box, widthLinesBoxes,  kTRUE, colorComb2040Box);
        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysBErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatSysAErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphCombDRPi0FitSysBErr2040->Draw("E2same");
        graphCombDRPi0FitStatSysAErr2040->Draw("p,E1Z,same");

        TLegend* legendDRCombDiffErrRepPad2 = new TLegend(0.15,0.86-1.1*0.85*textsizeLabelsPad2*3,0.35,0.86);
        legendDRCombDiffErrRepPad2->SetFillStyle(0);
        legendDRCombDiffErrRepPad2->SetFillColor(0);
        legendDRCombDiffErrRepPad2->SetLineColor(0);
        legendDRCombDiffErrRepPad2->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRCombDiffErrRepPad2->SetMargin(0.25);
        legendDRCombDiffErrRepPad2->SetTextFont(42);
        legendDRCombDiffErrRepPad2->AddEntry(graphCombDRPi0FitStatSysAErr2040,"Stat. #oplus Syst. A","pe");
        legendDRCombDiffErrRepPad2->AddEntry(graphCombDRPi0FitSysBErr2040,"Syst. B","f");
        legendDRCombDiffErrRepPad2->AddEntry(graphCombDRPi0FitSysCErr2040,"Syst. C","f");
        legendDRCombDiffErrRepPad2->Draw();

        labelDRCent2040->Draw();
        labelALICECent2040->Draw();
    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");
        TBox* boxSysCErrComb2050 = CreateBoxConv(colorComb2050Box, 0.74, 1.-(SysCCombDRPi0Fit2050) , 0.83, 1.+(SysCCombDRPi0Fit2050));
        boxSysCErrComb2050->Draw();
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysCErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050Box , colorComb2050Box, widthLinesBoxes,  kTRUE, colorComb2050Box);
        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysBErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatSysAErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        graphCombDRPi0FitSysBErr2050->Draw("E2same");
        graphCombDRPi0FitStatSysAErr2050->Draw("p,E1Z,same");

        TLegend* legendDRCombDiffErrRepPad3 = new TLegend(0.15,0.88-1.1*0.85*textsizeLabelsPad3*3,0.35,0.88);
        legendDRCombDiffErrRepPad3->SetFillStyle(0);
        legendDRCombDiffErrRepPad3->SetFillColor(0);
        legendDRCombDiffErrRepPad3->SetLineColor(0);
        legendDRCombDiffErrRepPad3->SetTextSize(0.85*textsizeLabelsPad3);
        legendDRCombDiffErrRepPad3->SetMargin(0.25);
        legendDRCombDiffErrRepPad3->SetTextFont(42);
        legendDRCombDiffErrRepPad3->AddEntry(graphCombDRPi0FitStatSysAErr2050,"Stat. #oplus Syst. A","pe");
        legendDRCombDiffErrRepPad3->AddEntry(graphCombDRPi0FitSysBErr2050,"Syst. B","f");
        legendDRCombDiffErrRepPad3->AddEntry(graphCombDRPi0FitSysCErr2050,"Syst. C","f");
        legendDRCombDiffErrRepPad3->Draw();

        labelDRCent2050->Draw();
        labelALICECent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_combMeasurement_diffErrRep.%s", outputDir.Data(), suffix.Data()));



    //*******************************************************************************************************************************************
    //************************************ plot ratio of individual measurements to combined fit for inc gamma **********************************
    //*******************************************************************************************************************************************
    canvasRatioIndDR->cd();

    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();

    //_______________________________________________________________ 0-10% dummy upper panel ___________________________________________________
    TH2D *dummyIncGammaIndMeas1 ;
    dummyIncGammaIndMeas1 = new TH2D("dummyIncGammaIndMeas1", "dummyIncGammaIndMeas1", 1000, 0., 22, 1000., indMeasRatio[0], indMeasRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyIncGammaIndMeas1, "#it{p}_{T} (GeV/#it{c})", "",
                            0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.95,0.105/(textsizeFacPad1*margin), 510, 505);
    dummyIncGammaIndMeas1->GetXaxis()->SetLabelOffset(-0.015);
    dummyIncGammaIndMeas1->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyIncGammaIndMeas1->GetXaxis()->SetTickLength(0.06);
    dummyIncGammaIndMeas1->GetYaxis()->SetTickLength(0.028);

    //_______________________________________________________________ 20-40% dummy middle panel _________________________________________________
    TH2D *dummyIncGammaIndMeas2 ;
    dummyIncGammaIndMeas2 = new TH2D("dummyIncGammaIndMeas2", "dummyIncGammaIndMeas2", 1000, 0., 22, 1000., indMeasRatio[0], indMeasRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyIncGammaIndMeas2, "#it{p}_{T} (GeV/#it{c})", "#gamma_{inc} (Data / Combined)",
                            0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.95,0.105/(textsizeFacPad2*margin), 510, 505);
    dummyIncGammaIndMeas2->GetXaxis()->SetLabelOffset(-0.015);
    dummyIncGammaIndMeas2->GetYaxis()->CenterTitle(kTRUE);
    dummyIncGammaIndMeas2->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyIncGammaIndMeas2->GetXaxis()->SetTickLength(0.06);
    dummyIncGammaIndMeas2->GetYaxis()->SetTickLength(0.028);

    //_______________________________________________________________ 20-50% dummy lower panel __________________________________________________
    TH2D *dummyIncGammaIndMeas3 ;
    dummyIncGammaIndMeas3 = new TH2D("dummyIncGammaIndMeas3", "dummyIncGammaIndMeas3", 1000, 0., 22, 1000., indMeasRatio[0], indMeasRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyIncGammaIndMeas3, "#it{p}_{T} (GeV/#it{c})", "",
                            0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.95,0.105/(textsizeFacPad3*margin), 510, 505);
    dummyIncGammaIndMeas3->GetXaxis()->SetLabelOffset(-0.015);
    dummyIncGammaIndMeas3->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyIncGammaIndMeas3->GetXaxis()->SetTickLength(0.055);
    dummyIncGammaIndMeas3->GetYaxis()->SetTickLength(0.035);


    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyIncGammaIndMeas1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysErr0010, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
        ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatErr0010);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatErr0010, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        graphRatioCombPPHOSIncGammaSysErr0010->Draw("E2same");
        graphRatioCombPPHOSIncGammaStatErr0010->Draw("p,E1Z,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysErr0010, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
        ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatErr0010);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatErr0010, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
        graphRatioCombPPCMIncGammaSysErr0010->Draw("E2same");
        graphRatioCombPPCMIncGammaStatErr0010->Draw("p,E1Z,same");

        TLatex *labelALICEIncGamma = new TLatex(0.12,0.76,"ALICE");
        SetStyleTLatex( labelALICEIncGamma, 0.85*textsizeLabelsPad1,4);
        labelALICEIncGamma->Draw();

        TLegend* legendIncGammaRIndMeas = new TLegend(0.8,0.92-1.1*0.85*textsizeLabelsPad1*2,0.95,0.92);
        legendIncGammaRIndMeas->SetFillStyle(0);
        legendIncGammaRIndMeas->SetFillColor(0);
        legendIncGammaRIndMeas->SetLineColor(0);
        legendIncGammaRIndMeas->SetTextSize(0.85*textsizeLabelsPad1);
        legendIncGammaRIndMeas->SetMargin(0.4);
        legendIncGammaRIndMeas->SetTextFont(42);
        legendIncGammaRIndMeas->AddEntry(graphRatioCombPPCMIncGammaSysErr0010,"PCM","pf");
        legendIncGammaRIndMeas->AddEntry(graphRatioCombPPHOSIncGammaSysErr0010,"PHOS","pf");
        legendIncGammaRIndMeas->Draw();

        labelDRCent0010->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyIncGammaIndMeas2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
        ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatErr2040);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        graphRatioCombPPHOSIncGammaSysErr2040->Draw("E2same");
        graphRatioCombPPHOSIncGammaStatErr2040->Draw("p,E1Z,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysErr2040, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
        ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatErr2040);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatErr2040, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
        graphRatioCombPPCMIncGammaSysErr2040->Draw("E2same");
        graphRatioCombPPCMIncGammaStatErr2040->Draw("p,E1Z,same");

        labelDRCent2040->Draw();

    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyIncGammaIndMeas3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
        ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatErr2040);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        graphRatioCombPPHOSIncGammaSysErr2040->Draw("E2same");
        graphRatioCombPPHOSIncGammaStatErr2040->Draw("p,E1Z,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysErr2050, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
        ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatErr2050);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatErr2050, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
        graphRatioCombPPCMIncGammaSysErr2050->Draw("E2same");
        graphRatioCombPPCMIncGammaStatErr2050->Draw("p,E1Z,same");

        labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/Ratio_indMeasurement_incGamma.%s", outputDir.Data(), suffix.Data()));

    //*******************************************************************************************************************************************
    //************************************ plot ratio of individual measurements to combined fit for inc gamma **********************************
    //*******************************************************************************************************************************************
    canvasRatioIndDR->cd();

    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();

    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr0010Plot = (TGraphAsymmErrors*)graphRatioCombPPHOSIncGammaStatSysAErr0010->Clone("graphRatioCombPPHOSIncGammaStatSysAErr0010Plot");
    ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatSysAErr0010Plot);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr2040Plot = (TGraphAsymmErrors*)graphRatioCombPPHOSIncGammaStatSysAErr2040->Clone("graphRatioCombPPHOSIncGammaStatSysAErr2040Plot");
    ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatSysAErr2040Plot);
    TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr2050Plot = (TGraphAsymmErrors*)graphRatioCombPPHOSIncGammaStatSysAErr2050->Clone("graphRatioCombPPHOSIncGammaStatSysAErr2050Plot");
    ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatSysAErr2050Plot);

    TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr0010Plot = (TGraphAsymmErrors*)graphRatioCombPPCMIncGammaStatSysAErr0010->Clone("graphRatioCombPPCMIncGammaStatSysAErr0010Plot");
    ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatSysAErr0010Plot);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr2040Plot = (TGraphAsymmErrors*)graphRatioCombPPCMIncGammaStatSysAErr2040->Clone("graphRatioCombPPCMIncGammaStatSysAErr2040Plot");
    ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatSysAErr2040Plot);
    TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr2050Plot = (TGraphAsymmErrors*)graphRatioCombPPCMIncGammaStatSysAErr2050->Clone("graphRatioCombPPCMIncGammaStatSysAErr2050Plot");
    ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatSysAErr2050Plot);

    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyIncGammaIndMeas1->Draw("");

        boxSysCErrPCM0010->Draw();
        TBox* boxSysCErrPHOS0010 = CreateBoxConv(colorPHOSBox, 0.815, 1.-(SysCPHOSIncGamma0010) , 0.88, 1.+(SysCPHOSIncGamma0010));
        boxSysCErrPHOS0010->Draw();

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysBErr0010, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatSysAErr0010Plot, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        graphRatioCombPPHOSIncGammaSysBErr0010->Draw("E2same");
        graphRatioCombPPHOSIncGammaStatSysAErr0010Plot->Draw("p,E1Z,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysCErr0010, markerStylePCM, markerSizePCM, colorPCMBox , colorPCMBox,widthLinesBoxes, kTRUE, colorPCMBox);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysBErr0010, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatSysAErr0010Plot, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
        graphRatioCombPPCMIncGammaSysBErr0010->Draw("E2same");
        graphRatioCombPPCMIncGammaStatSysAErr0010Plot->Draw("p,E1Z,same");

        legendIncGammaRIndMeas->Draw();
        dummyIncGammaIndMeas1->Draw("axis,same");
        labelDRCent0010->Draw();
        labelALICEIncGamma->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyIncGammaIndMeas2->Draw("");

        TLegend* legendIncGammaRIndMeas2 = new TLegend(0.65,0.95-1.1*0.85*textsizeLabelsPad2*3,0.9,0.95);
        legendIncGammaRIndMeas2->SetFillStyle(0);
        legendIncGammaRIndMeas2->SetFillColor(0);
        legendIncGammaRIndMeas2->SetLineColor(0);
        legendIncGammaRIndMeas2->SetTextSize(0.85*textsizeLabelsPad2);
        legendIncGammaRIndMeas2->SetMargin(0.45);
        legendIncGammaRIndMeas2->SetTextFont(42);
        legendIncGammaRIndMeas2->SetNColumns(2);
        legendIncGammaRIndMeas2->SetHeader("For both methods:");
        legendIncGammaRIndMeas2->AddEntry(graphRatioCombPPCMIncGammaStatSysAErr0010Plot,"Stat. #oplus Syst. A","pel");
        legendIncGammaRIndMeas2->AddEntry((TObject*)0,"","");
        legendIncGammaRIndMeas2->AddEntry(graphRatioCombPPCMIncGammaSysBErr0010,"Syst. B","f");
        legendIncGammaRIndMeas2->AddEntry(graphRatioCombPPCMIncGammaSysCErr0010,"Syst. C","f");
        legendIncGammaRIndMeas2->Draw();

        boxSysCErrPCM2040->Draw();
        TBox* boxSysCErrPHOS2040 = CreateBoxConv(colorPHOSBox, 0.815, 1.-(SysCPHOSIncGamma2040) , 0.88, 1.+(SysCPHOSIncGamma2040));
        boxSysCErrPHOS2040->Draw();

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysBErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatSysAErr2040Plot, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        graphRatioCombPPHOSIncGammaSysBErr2040->Draw("E2same");
        graphRatioCombPPHOSIncGammaStatSysAErr2040Plot->Draw("p,E1Z,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysBErr2040, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatSysAErr2040Plot, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
        graphRatioCombPPCMIncGammaSysBErr2040->Draw("E2same");
        graphRatioCombPPCMIncGammaStatSysAErr2040Plot->Draw("p,E1Z,same");

        dummyIncGammaIndMeas2->Draw("axis,same");
        labelDRCent2040->Draw();

    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyIncGammaIndMeas3->Draw("");
        boxSysCErrPCM2050->Draw();
        TBox* boxSysCErrPHOS2050 = CreateBoxConv(colorPHOSBox, 0.815, 1.-(SysCPHOSIncGamma2050) , 0.88, 1.+(SysCPHOSIncGamma2050));
        boxSysCErrPHOS2050->Draw();

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysBErr2050, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatSysAErr2050Plot, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        graphRatioCombPPHOSIncGammaSysBErr2050->Draw("E2same");
        graphRatioCombPPHOSIncGammaStatSysAErr2050Plot->Draw("p,E1Z,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysBErr2050, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatSysAErr2050Plot, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
        graphRatioCombPPCMIncGammaSysBErr2050->Draw("E2same");
        graphRatioCombPPCMIncGammaStatSysAErr2050Plot->Draw("p,E1Z,same");

        dummyIncGammaIndMeas3->Draw("axis,same");

        labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/Ratio_indMeasurementDiffErrorRep_incGamma.%s", outputDir.Data(), suffix.Data()));

    //*******************************************************************************************************************************************
    //************************************ plot ratio of dirGamma measurement to combined fit for  **********************************************
    //*******************************************************************************************************************************************
    canvasRatioIndDR->cd();

    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();

    //_______________________________________________________________ 0-10% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyIncGammaIndMeas1->Draw("");

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
        if (graphRatioCombFitDirGammaSysErr0010){
            DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaSysErr0010, markerStyleComb0010, markerSizeComb2040, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
            graphRatioCombFitDirGammaSysErr0010->Draw("E2same");
        }
        if (graphRatioCombFitDirGammaStatErr0010){
            DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaStatErr0010, markerStyleComb0010, markerSizeComb2040, colorComb0010 , colorComb0010);
            graphRatioCombFitDirGammaStatErr0010->Draw("p,E1Z,same");
        }
        legendIncGammaRIndMeas->Draw();
        dummyIncGammaIndMeas1->Draw("axis,same");
        labelDRCent0010->Draw();
        labelALICEIncGamma->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyIncGammaIndMeas2->Draw("");

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaStatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphRatioCombFitDirGammaSysErr2040->Draw("E2same");
        graphRatioCombFitDirGammaStatErr2040->Draw("p,E1Z,same");

        dummyIncGammaIndMeas2->Draw("axis,same");
        labelDRCent2040->Draw();

    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyIncGammaIndMeas3->Draw("");

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaSysErr2050, markerStyleComb2050, markerSizeComb2040, colorComb2050 , colorComb2050,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaStatErr2050, markerStyleComb2050, markerSizeComb2040, colorComb2050 , colorComb2050);
        graphRatioCombFitDirGammaSysErr2050->Draw("E2same");
        graphRatioCombFitDirGammaStatErr2050->Draw("p,E1Z,same");

        dummyIncGammaIndMeas3->Draw("axis,same");

        labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/Ratio_DirGammaCombMeasurementToFit.%s", outputDir.Data(), suffix.Data()));


    //*******************************************************************************************************************************************
    //*************************************************** Plotting inclusive Gamma Spectrum *****************************************************
    //*******************************************************************************************************************************************

    TGraphAsymmErrors* graphCombIncGammaSysErr0010Plot = ScaleGraph(graphCombIncGammaSysErr0010,100);
    TGraphAsymmErrors* graphCombIncGammaStatErr0010Plot = ScaleGraph(graphCombIncGammaStatErr0010,100);

    TGraphAsymmErrors* graphCombIncGammaSysErr2040Plot = ScaleGraph(graphCombIncGammaSysErr2040,10);
    graphCombIncGammaSysErr2040Plot->RemovePoint(0);
    TGraphAsymmErrors* graphCombIncGammaStatErr2040Plot = ScaleGraph(graphCombIncGammaStatErr2040,10);
    graphCombIncGammaStatErr2040Plot->RemovePoint(0);

    TGraphAsymmErrors* graphCombIncGammaSysErr2050Plot = ScaleGraph(graphCombIncGammaSysErr2050,1);
    TGraphAsymmErrors* graphCombIncGammaStatErr2050Plot = ScaleGraph(graphCombIncGammaStatErr2050,1);

    TH1D* histoFitDummyPlotting = (TH1D*) histoFitQCDIncGammaComb2040->Clone("histoFitDummyPlotting");

    TCanvas *canvasIncGamma = new TCanvas("canvasIncGamma","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasIncGamma, 0.16, 0.01, 0.01, 0.07);
    canvasIncGamma->SetLogy();
    canvasIncGamma->SetLogx();

    textSizeLabelsPixelIncGam = 48;
    textsizeLabelsIncGamma = 0;
    if (canvasIncGamma->XtoPixel(canvasIncGamma->GetX2()) < canvasIncGamma->YtoPixel(canvasIncGamma->GetY1())){
        textsizeLabelsIncGamma = (Double_t)textSizeLabelsPixelIncGam/canvasIncGamma->XtoPixel(canvasIncGamma->GetX2()) ;
    } else {
        textsizeLabelsIncGamma = (Double_t)textSizeLabelsPixelIncGam/canvasIncGamma->YtoPixel(canvasIncGamma->GetY1());
    }

    TH2D *dummyGamma ;
    dummyGamma = new TH2D("dummyGamma", "dummyGamma", 1000, 0., 22, 1000., 6e-7,1.2e4);
    SetStyleHistoTH2ForGraphs( dummyGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{inc}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                            0.85*textsizeLabelsIncGamma, textsizeLabelsIncGamma, 0.85*textsizeLabelsIncGamma, textsizeLabelsIncGamma, 0.75, 1.6);
    dummyGamma->GetXaxis()->SetLabelOffset(-0.015);
    dummyGamma->GetXaxis()->SetRangeUser(0.5,30);
//     dummyGamma->GetXaxis()->SetRangeUser(0,16);
    dummyGamma->DrawCopy();

    DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSysErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010,widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
    graphCombIncGammaSysErr0010Plot->Draw("E2same");
    graphCombIncGammaStatErr0010Plot->Draw("p,E1Z,same");

    DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSysErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
    graphCombIncGammaSysErr2040Plot->Draw("E2same");
    graphCombIncGammaStatErr2040Plot->Draw("p,E1Z,same");

    DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSysErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
    graphCombIncGammaSysErr2050Plot->Draw("E2same");
    graphCombIncGammaStatErr2050Plot->Draw("p,E1Z,same");

    labelScalingIncGamma0010->Draw();
    labelScalingIncGamma2040->Draw();

    TLegend* legendIncGamma = new TLegend(0.6,0.9-1.*0.85*textsizeLabelsIncGamma*3,0.9,0.9);
    legendIncGamma->SetFillStyle(0);
    legendIncGamma->SetFillColor(0);
    legendIncGamma->SetLineColor(0);
    legendIncGamma->SetTextSize(0.85*textsizeLabelsIncGamma);
    legendIncGamma->SetMargin(0.2);
    legendIncGamma->SetTextFont(42);
    legendIncGamma->AddEntry(graphCombIncGammaSysErr0010Plot,"  0-10% ALICE","pf");
    legendIncGamma->AddEntry(graphCombIncGammaSysErr2040Plot,"20-40% ALICE","pf");
    legendIncGamma->AddEntry(graphCombIncGammaSysErr2050Plot,"20-50% ALICE","pf");
//     legendIncGamma->AddEntry(histoFitDummyPlotting,"Fits to #gamma_{inc}","l");
    legendIncGamma->Draw();
    dummyGamma->Draw("axis,same");
    canvasIncGamma->Print(Form("%s/IncGammaSpectrum.%s",outputDir.Data(),suffix.Data()));


    cout << "Plotting direct gamma spectra at " << __LINE__ << endl;
    //*******************************************************************************************************************************************
    //*************************************************** Plotting direct Gamma Spectrum ********************************************************
    //*******************************************************************************************************************************************

    TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr0010Plot = NULL;
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr0010Plot = NULL;
    TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr0010Plot = NULL;
    TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr0010ArPlot = NULL;
    if (graphCombDirGammaSpectrumSumErr0010) graphCombDirGammaSpectrumSumErr0010Plot         = ScaleGraph(graphCombDirGammaSpectrumSumErr0010,100);
    if (graphCombDirGammaSpectrumStatErr0010) graphCombDirGammaSpectrumStatErr0010Plot       = ScaleGraph(graphCombDirGammaSpectrumStatErr0010,100);
    if (graphCombDirGammaSpectrumStatErr0010Plot) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErr0010Plot);
    if (graphCombDirGammaSpectrumSystErr0010) graphCombDirGammaSpectrumSystErr0010Plot       = ScaleGraph(graphCombDirGammaSpectrumSystErr0010,100);
    if (graphCombDirGammaSpectrumSumErr0010Ar) graphCombDirGammaSpectrumSumErr0010ArPlot     = ScaleGraph(graphCombDirGammaSpectrumSumErr0010Ar,100);
    TH1D* histoFitThermalGamma0010Stat                                                       = (TH1D*)fitThermalGamma0010Stat->GetHistogram();
    histoFitThermalGamma0010Stat->Scale(100);
    TH1D* histoFitPureThermalGamma0010Stat                                                   = NULL;
    if (graphCombThermalGammaSpectrumStatErr0010){
        histoFitPureThermalGamma0010Stat                                                     = (TH1D*)fitPureThermalGamma0010Stat->GetHistogram();
        histoFitPureThermalGamma0010Stat->Scale(100);
    }

    TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2040Plot = NULL;
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr2040Plot = NULL;
    TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr2040Plot = NULL;
    TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2040ArPlot = NULL;
    if (graphCombDirGammaSpectrumSumErr2040) graphCombDirGammaSpectrumSumErr2040Plot         = ScaleGraph(graphCombDirGammaSpectrumSumErr2040,10);
    if (graphCombDirGammaSpectrumStatErr2040) graphCombDirGammaSpectrumStatErr2040Plot       = ScaleGraph(graphCombDirGammaSpectrumStatErr2040,10);
    if (graphCombDirGammaSpectrumStatErr2040Plot) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErr2040Plot);
    if (graphCombDirGammaSpectrumSystErr2040) graphCombDirGammaSpectrumSystErr2040Plot       = ScaleGraph(graphCombDirGammaSpectrumSystErr2040,10);
    if (graphCombDirGammaSpectrumSumErr2040Ar) graphCombDirGammaSpectrumSumErr2040ArPlot     = ScaleGraph(graphCombDirGammaSpectrumSumErr2040Ar,10);
    TH1D* histoFitThermalGamma2040Stat                                                       = (TH1D*)fitThermalGamma2040Stat->GetHistogram();
    histoFitThermalGamma2040Stat->Scale(10);

    TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2050Plot = NULL;
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr2050Plot = NULL;
    TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr2050Plot = NULL;
    TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2050ArPlot = NULL;
    if (graphCombDirGammaSpectrumSumErr2050) graphCombDirGammaSpectrumSumErr2050Plot         = ScaleGraph(graphCombDirGammaSpectrumSumErr2050,1);
    if (graphCombDirGammaSpectrumStatErr2050) graphCombDirGammaSpectrumStatErr2050Plot       = ScaleGraph(graphCombDirGammaSpectrumStatErr2050,1);
    if (graphCombDirGammaSpectrumStatErr2050Plot) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErr2050Plot);
    if (graphCombDirGammaSpectrumSystErr2050) graphCombDirGammaSpectrumSystErr2050Plot       = ScaleGraph(graphCombDirGammaSpectrumSystErr2050,1);
    if (graphCombDirGammaSpectrumSumErr2050Ar) graphCombDirGammaSpectrumSumErr2050ArPlot     = ScaleGraph(graphCombDirGammaSpectrumSumErr2050Ar,1);


    TCanvas *canvasDirGamma = new TCanvas("canvasDirGamma","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasDirGamma, 0.165, 0.01, 0.01, 0.07);
    canvasDirGamma->SetLogy();
    canvasDirGamma->SetLogx();

    textSizeLabelsPixelDirGam = 48;
    textsizeLabelsDirGamma = 0;
    if (canvasDirGamma->XtoPixel(canvasDirGamma->GetX2()) < canvasDirGamma->YtoPixel(canvasDirGamma->GetY1())){
        textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGamma->XtoPixel(canvasDirGamma->GetX2()) ;
    } else {
        textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGamma->YtoPixel(canvasDirGamma->GetY1());
    }

        TH2D *dummyDirGamma ;
        dummyDirGamma = new TH2D("dummyDirGamma", "dummyDirGamma", 1000, 0., 22, 1000., 1.2e-7,1.5e5);
        SetStyleHistoTH2ForGraphs( dummyDirGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                                0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
        dummyDirGamma->GetXaxis()->SetLabelOffset(-0.015);
        dummyDirGamma->GetXaxis()->SetTickLength(0.025);
        dummyDirGamma->GetYaxis()->SetTickLength(0.025);
        dummyDirGamma->GetXaxis()->SetRangeUser(0.5,30.);
    //     dummyDirGamma->GetXaxis()->SetRangeUser(0,16);
        dummyDirGamma->DrawCopy();

        TLatex *labelScalingDirGammaComb0010 = new TLatex(11.2,7.E-4,"x 10^{2}");
        SetStyleTLatex( labelScalingDirGammaComb0010, 0.85*textsizeLabelsDirGamma,4,colorComb0010,42,kFALSE);
        TLatex *labelScalingDirGammaComb2040 = new TLatex(11.2,3.5E-5,"x 10^{1}");
        SetStyleTLatex( labelScalingDirGammaComb2040, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
        //labelScalingDirGamma2050->Draw();
        //labelDirGammaColl->Draw();
        SetStyleHisto(histoFitDummyPlotting, widthCommonFit, 5, kGray+1);
        labelScalingDirGammaComb0010->Draw();
        labelScalingDirGammaComb2040->Draw();
    //
        TLegend* legendDirGamma = new TLegend(0.24,0.96-1.*0.85*textsizeLabelsDirGamma*3,0.24+0.21,0.96);//0.24,0.93-1.*0.85*textsizeLabelsDirGamma*3,0.24+0.21,0.93);
        legendDirGamma->SetFillStyle(0);
        legendDirGamma->SetFillColor(0);
        legendDirGamma->SetLineColor(0);
        legendDirGamma->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGamma->SetMargin(0.25);
        legendDirGamma->SetTextFont(42);
        legendDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot,"  0-10% ALICE","pf");
        legendDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
        legendDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr2050Plot,"20-50% ALICE","pf");
//         legendDirGamma->Draw();

        if (graphCombDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr0010Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphCombDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr0010ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr0010ArPlot , 1, 3, colorComb0010, colorComb0010, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErr0010ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr0010ArPlot);
        }
        if (graphCombDirGammaSpectrumSystErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
        }
        if (graphCombDirGammaSpectrumSystErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphCombDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2050ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2050ArPlot , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2050ArPlot);
        }


        TLegend* legendDirGammaForThesis = new TLegend(0.2,0.96-1.*0.85*textsizeLabelsDirGamma*4,0.24+0.21,0.96);//0.2,0.93-1.*0.85*textsizeLabelsDirGamma*4,0.24+0.21,0.93);
        legendDirGammaForThesis->SetFillStyle(0);
        legendDirGammaForThesis->SetFillColor(0);
        legendDirGammaForThesis->SetLineColor(0);
        legendDirGammaForThesis->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaForThesis->SetMargin(0.2);
        legendDirGammaForThesis->SetTextFont(42);
        legendDirGammaForThesis->SetHeader(collisionSystem.Data());
        legendDirGammaForThesis->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot,"  0-10%","pf");
        legendDirGammaForThesis->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40%","pf");
        legendDirGammaForThesis->AddEntry(graphCombDirGammaSpectrumSystErr2050Plot,"20-50%","pf");
        legendDirGammaForThesis->Draw();
//         labelThesisDGwithNLO->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum.%s",outputDir.Data(),suffix.Data()));


    cout << "Plotting direct gamma spectra with models at " << __LINE__ << endl;
        graphTheoryEPS090010Plot->Draw("p3lsame");
        graphTheoryCT100010Plot->Draw("p3lsame");
        graphTheoryNLO0010Plot->Draw("p3lsame");
    //     graphTheoryPromptMcGill0020Plot->Draw("lsame");

        if (graphCombDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr0010Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphCombDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }

        graphTheoryEPS092040Plot->Draw("p3lsame");
        graphTheoryCT102040Plot->Draw("p3lsame");
        graphTheoryNLO2040Plot->Draw("p3lsame");
//         graphTheoryPromptMcGill2040Plot->Draw("lsame");

        if (graphCombDirGammaSpectrumSystErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
        }

        graphTheoryEPS092050Plot->Draw("p3lsame");
        graphTheoryCT102050Plot->Draw("p3lsame");
        graphTheoryNLO2050Plot->Draw("p3lsame");

        if (graphCombDirGammaSpectrumSystErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphCombDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2050ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2050ArPlot , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2050ArPlot);
        }

        legendDirGammaNLO->Draw();
        legendDirGammaNLO2->Draw();
        labelScalingDirGammaComb0010->Draw();
        labelScalingDirGammaComb2040->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum_withNLO.%s",outputDir.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    //************************************************************************************************************************************

    TH2D *dummyDirGammaTheory = new TH2D("dummyDirGammaTheory", "dummyDirGamma", 1000, 0., 22, 1000., 4e-8,9e2);
    SetStyleHistoTH2ForGraphs( dummyDirGammaTheory, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})", 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummyDirGammaTheory->GetXaxis()->SetLabelOffset(-0.015);
    dummyDirGammaTheory->GetXaxis()->SetTickLength(0.025);
    dummyDirGammaTheory->GetYaxis()->SetTickLength(0.025);
    dummyDirGammaTheory->GetXaxis()->SetRangeUser(doubleRatioX[0],incSpectraX[1]);
    dummyDirGammaTheory->DrawCopy();

        if (graphCombDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        }
        if (graphCombDirGammaSpectrumSystErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , styleNLOMcGill, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
        }

        TLatex *labelDirGammaCollRedX = new TLatex(0.6,0.93,"This thesis");
        SetStyleTLatex( labelDirGammaCollRedX, 0.85*textsizeLabelsDirGamma,4);
//         labelDirGammaCollRedX->Draw();
        TLatex *labelDirGammaLowLeft = new TLatex(0.21,0.27,"This thesis");
        SetStyleTLatex( labelDirGammaLowLeft, 0.85*textsizeLabelsDirGamma,4);
        labelDirGammaLowLeft->Draw();

        TLegend* legendDirGammaRedX = new TLegend(0.6,0.92-1.1*0.85*textsizeLabelsDirGamma*4,0.6+0.21,0.92);
        legendDirGammaRedX->SetFillStyle(0);
        legendDirGammaRedX->SetFillColor(0);
        legendDirGammaRedX->SetLineColor(0);
        legendDirGammaRedX->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaRedX->SetMargin(0.3);
        legendDirGammaRedX->SetTextFont(42);
        legendDirGammaRedX->SetHeader(collisionSystem.Data());
        legendDirGammaRedX->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot,"  0-10% ALICE","pf");
        legendDirGammaRedX->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
        legendDirGammaRedX->AddEntry(graphCombDirGammaSpectrumSystErr2050Plot,"20-50% ALICE","pf");
//         legendDirGammaRedX->Draw();
        TLegend* legendDirGammaLowLeft = new TLegend(0.2,0.12,0.555+0.21,0.12+1.*0.85*textsizeLabelsDirGamma*4);
        legendDirGammaLowLeft->SetFillStyle(0);
        legendDirGammaLowLeft->SetFillColor(0);
        legendDirGammaLowLeft->SetLineColor(0);
        legendDirGammaLowLeft->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaLowLeft->SetMargin(0.15);
        legendDirGammaLowLeft->SetTextFont(42);
        legendDirGammaLowLeft->SetHeader(collisionSystem.Data());
        legendDirGammaLowLeft->AddEntry(graphPCMDirGammaSpectrumSystErr0010Plot,"  0-10%","pf");
        legendDirGammaLowLeft->AddEntry(graphPCMDirGammaSpectrumSystErr2040Plot,"20-40%","pf");
        legendDirGammaLowLeft->AddEntry(graphPCMDirGammaSpectrumSystErr2050Plot,"20-50%","pf");
        legendDirGammaLowLeft->Draw();

        if (graphCombDirGammaSpectrumSystErr0010Plot){
            graphCombDirGammaSpectrumSystErr0010Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot){
            graphCombDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }

        if (graphCombDirGammaSpectrumSystErr2040Plot){
            graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot){
            graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot){
            graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
        }
        if (graphCombDirGammaSpectrumSystErr2050Plot){
            graphCombDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2050Plot){
            graphCombDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2050ArPlot){
            graphCombDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2050ArPlot);
        }
        labelScalingDirGammaComb0010->Draw();
        labelScalingDirGammaComb2040->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum_LogX.%s",outputDir.Data(),suffix.Data()));

    legendDirGammaPCMOnlyTheory->Draw();

        graphTheoryChatterjeeSummed0010Plot->Draw("p3lsame");
        graphTheoryPHSDplusPrompt0010Plot->Draw("p3lsame");
        //graphTheoryChatterjee0010Plot->Draw("p3lsame");
        graphTheoryRapp0010Plot->Draw("p3lsame");
        graphTheoryMcGill0010Plot->Draw("p3lsame");

        if (graphCombDirGammaSpectrumSystErr0010Plot){
            graphCombDirGammaSpectrumSystErr0010Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot){
            graphCombDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }
        graphTheoryMcGill2040Plot->Draw("p3lsame");
        graphTheoryPHSD2040Plot->Draw("p3lsame");
        graphTheoryChatterjee2040Plot->Draw("p3lsame");
//         graphTheoryHees2040Plot->Draw("p3lsame");
        graphTheoryHe2040Plot->Draw("p3lsame");
        graphTheoryRapp2040Plot->Draw("p3lsame");

        if (graphCombDirGammaSpectrumSystErr2040Plot){
            graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot){
            graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot){
            graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
        }

        //graphTheoryChatterjee2050Plot->Draw("p3lsame");
        graphTheoryChatterjeeSummed2050Plot->Draw("p3lsame");
        graphTheoryPHSDplusPrompt2050Plot->Draw("p3lsame");
        if (graphCombDirGammaSpectrumSystErr2050Plot){
            graphCombDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2050Plot){
            graphCombDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2050ArPlot){
            graphCombDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2050ArPlot);
        }
        TLegend* legendDirGammaTheoryRapp_2 = new TLegend(0.6,0.12,0.94,0.12+1.*0.85*textsizeLabelsDirGamma*2);
        legendDirGammaTheoryRapp_2->SetFillStyle(0);
        legendDirGammaTheoryRapp_2->SetFillColor(0);
        legendDirGammaTheoryRapp_2->SetLineColor(0);
//         legendDirGammaTheoryRapp_2->SetNColumns(2);
        legendDirGammaTheoryRapp_2->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaTheoryRapp_2->SetMargin(0.15);
        legendDirGammaTheoryRapp_2->SetTextFont(42);
        legendDirGammaTheoryRapp_2->AddEntry(graphTheoryRapp0010Plot,"Rapp et al.","l");
        legendDirGammaTheoryRapp_2->AddEntry((TObject*)0,"arXiv:","");
//         legendDirGammaTheoryRapp_2->AddEntry((TObject*)0,"","");
//         legendDirGammaTheoryRapp_2->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryRapp_2->Draw();

        legendDirGammaLowLeft->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusMcGill_LogX.%s",outputDir.Data(),suffix.Data()));

    //*****************************************************************************************************************
    //******************************* Plotting spectrum with reduced x range ******************************************
    //*****************************************************************************************************************
    // Cut plotting graphs
    TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr0010Plot2 = NULL;
    if (graphCombDirGammaSpectrumSystErr0010Plot) {
        graphCombDirGammaSpectrumSystErr0010Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr0010Plot->Clone("graphCombDirGammaSpectrumSystErr0010Plot2");
        while (graphCombDirGammaSpectrumSystErr0010Plot2->GetX()[graphCombDirGammaSpectrumSystErr0010Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumSystErr0010Plot2->GetN()>0)
            graphCombDirGammaSpectrumSystErr0010Plot2->RemovePoint(graphCombDirGammaSpectrumSystErr0010Plot2->GetN()-1);
        if (graphCombDirGammaSpectrumSystErr0010Plot2->GetN()==0) graphCombDirGammaSpectrumSystErr0010Plot2 = NULL;
    }
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr0010Plot2 = NULL;
    if (graphCombDirGammaSpectrumStatErr0010Plot) {
        graphCombDirGammaSpectrumStatErr0010Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr0010Plot->Clone("graphCombDirGammaSpectrumStatErr0010Plot2");
        while (graphCombDirGammaSpectrumStatErr0010Plot2->GetX()[graphCombDirGammaSpectrumStatErr0010Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumStatErr0010Plot2->GetN()>0)
            graphCombDirGammaSpectrumStatErr0010Plot2->RemovePoint(graphCombDirGammaSpectrumStatErr0010Plot2->GetN()-1);
        if (graphCombDirGammaSpectrumStatErr0010Plot2->GetN()==0) graphCombDirGammaSpectrumStatErr0010Plot2 = NULL;
    }
    TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr2040Plot2 = NULL;
    if (graphCombDirGammaSpectrumSystErr2040Plot) {
        graphCombDirGammaSpectrumSystErr2040Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr2040Plot->Clone("graphCombDirGammaSpectrumSystErr2040Plot2");
        while (graphCombDirGammaSpectrumSystErr2040Plot2->GetX()[graphCombDirGammaSpectrumSystErr2040Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumSystErr2040Plot2->GetN()>0)
            graphCombDirGammaSpectrumSystErr2040Plot2->RemovePoint(graphCombDirGammaSpectrumSystErr2040Plot2->GetN()-1);
        if (graphCombDirGammaSpectrumSystErr2040Plot2->GetN()==0) graphCombDirGammaSpectrumSystErr2040Plot2 = NULL;
    }
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr2040Plot2 = NULL;
    if (graphCombDirGammaSpectrumStatErr2040Plot) {
        graphCombDirGammaSpectrumStatErr2040Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr2040Plot->Clone("graphCombDirGammaSpectrumStatErr2040Plot2");
        while (graphCombDirGammaSpectrumStatErr2040Plot2->GetX()[graphCombDirGammaSpectrumStatErr2040Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumStatErr2040Plot2->GetN()>0)
            graphCombDirGammaSpectrumStatErr2040Plot2->RemovePoint(graphCombDirGammaSpectrumStatErr2040Plot2->GetN()-1);
        if (graphCombDirGammaSpectrumStatErr2040Plot2->GetN()==0) graphCombDirGammaSpectrumStatErr2040Plot2 = NULL;
    }

    TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2040ArPlot2 = NULL;
    if (graphCombDirGammaSpectrumSumErr2040ArPlot) {
        graphCombDirGammaSpectrumSumErr2040ArPlot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSumErr2040ArPlot->Clone("graphCombDirGammaSpectrumSumErr2040ArPlot2");
        while (graphCombDirGammaSpectrumSumErr2040ArPlot2->GetX()[graphCombDirGammaSpectrumSumErr2040ArPlot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumSumErr2040ArPlot2->GetN()>0)
            graphCombDirGammaSpectrumSumErr2040ArPlot2->RemovePoint(graphCombDirGammaSpectrumSumErr2040ArPlot2->GetN()-1);
        if (graphCombDirGammaSpectrumSumErr2040ArPlot2->GetN()==0) graphCombDirGammaSpectrumSumErr2040ArPlot2 = NULL;
    }

    TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr2050Plot2 = NULL;
    if (graphCombDirGammaSpectrumSystErr2050Plot) {
        graphCombDirGammaSpectrumSystErr2050Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr2050Plot->Clone("graphCombDirGammaSpectrumSystErr2050Plot2");
        while (graphCombDirGammaSpectrumSystErr2050Plot2->GetX()[graphCombDirGammaSpectrumSystErr2050Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumSystErr2050Plot2->GetN()>0)
            graphCombDirGammaSpectrumSystErr2050Plot2->RemovePoint(graphCombDirGammaSpectrumSystErr2050Plot2->GetN()-1);
        if (graphCombDirGammaSpectrumSystErr2050Plot2->GetN()==0) graphCombDirGammaSpectrumSystErr2050Plot2 = NULL;
    }
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr2050Plot2 = NULL;
    if (graphCombDirGammaSpectrumStatErr2050Plot) {
        graphCombDirGammaSpectrumStatErr2050Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr2050Plot->Clone("graphCombDirGammaSpectrumStatErr2050Plot2");
        while (graphCombDirGammaSpectrumStatErr2050Plot2->GetX()[graphCombDirGammaSpectrumStatErr2050Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumStatErr2050Plot2->GetN()>0)
            graphCombDirGammaSpectrumStatErr2050Plot2->RemovePoint(graphCombDirGammaSpectrumStatErr2050Plot2->GetN()-1);
        if (graphCombDirGammaSpectrumStatErr2050Plot2->GetN()==0) graphCombDirGammaSpectrumStatErr2050Plot2 = NULL;
    }

    TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2050ArPlot2 = NULL;
    if (graphCombDirGammaSpectrumSumErr2050ArPlot) {
        graphCombDirGammaSpectrumSumErr2050ArPlot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSumErr2050ArPlot->Clone("graphCombDirGammaSpectrumSumErr2050ArPlot2");
        while (graphCombDirGammaSpectrumSumErr2050ArPlot2->GetX()[graphCombDirGammaSpectrumSumErr2050ArPlot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumSumErr2050ArPlot2->GetN()>0)
            graphCombDirGammaSpectrumSumErr2050ArPlot2->RemovePoint(graphCombDirGammaSpectrumSumErr2050ArPlot2->GetN()-1);
        if (graphCombDirGammaSpectrumSumErr2050ArPlot2->GetN()==0) graphCombDirGammaSpectrumSumErr2050ArPlot2 = NULL;
    }

    cout << "Plotting direct gamma spectra with thermal fit at " << __LINE__ << endl;
    canvasDirGamma->SetLogx(0);
    TH2D *dummDirGammaRedX = new TH2D("dummDirGammaRedX", "dummDirGammaRedX", 100000, 0., 7., 1000., 9e-6,9e2);
    SetStyleHistoTH2ForGraphs( dummDirGammaRedX, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                            0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummDirGammaRedX->GetXaxis()->SetLabelOffset(-0.002);
    dummDirGammaRedX->GetXaxis()->SetTickLength(0.025);
    dummDirGammaRedX->GetYaxis()->SetTickLength(0.025);
    dummDirGammaRedX->GetXaxis()->SetRangeUser(0.5,5.);
    dummDirGammaRedX->DrawCopy();

        SetStyleHisto(histoFitThermalGamma0010Stat, 3, styleFit, colorComb0010+1 );
        SetStyleHisto(histoFitThermalGamma2040Stat, 3, styleFit, colorComb2040+1 );

        TLatex *labelScalingDirGammaComb0010_2 = new TLatex(4.,2.2E-1,"x 10^{2}");
        SetStyleTLatex( labelScalingDirGammaComb0010_2, 0.85*textsizeLabelsDirGamma,4,colorComb0010,42,kFALSE);
        TLatex *labelScalingDirGammaComb2040_2 = new TLatex(4.,1E-2,"x 10^{1}");
        SetStyleTLatex( labelScalingDirGammaComb2040_2, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
        labelScalingDirGammaComb0010_2->Draw();
        labelScalingDirGammaComb2040_2->Draw();

        labelDirGammaLowLeft ->Draw();
        legendDirGammaLowLeft->Draw();

        TLegend* legendDirGammaRedXWithFit = new TLegend(0.455,0.915-0.9*0.85*textsizeLabelsDirGamma*3,0.555+0.21,0.915);
        legendDirGammaRedXWithFit->SetFillStyle(0);
        legendDirGammaRedXWithFit->SetFillColor(0);
        legendDirGammaRedXWithFit->SetLineColor(0);
        legendDirGammaRedXWithFit->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaRedXWithFit->SetMargin(0.2);
        legendDirGammaRedXWithFit->SetTextFont(42);
        legendDirGammaRedXWithFit->SetHeader(collisionSystem.Data());
        legendDirGammaRedXWithFit->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot,"  0-10% ALICE","pf");
        legendDirGammaRedXWithFit->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
        legendDirGammaRedXWithFit->AddEntry(graphCombDirGammaSpectrumSystErr2050Plot,"20-50% ALICE","pf");
//         legendDirGammaRedXWithFit->Draw();


        if (graphCombDirGammaSpectrumSystErr0010Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr0010Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphCombDirGammaSpectrumStatErr0010Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSystErr2040Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot2, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr2040Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphCombDirGammaSpectrumStatErr2040Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , styleNLOMcGill, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErr2040ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot2);
        }
        if (graphCombDirGammaSpectrumSystErr2050Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2050Plot2, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr2050Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2050Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2050Plot2, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphCombDirGammaSpectrumStatErr2050Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2050ArPlot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2050ArPlot2 , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErr2050ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2050ArPlot2);
        }

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum_ReducedX.%s",outputDir.Data(),suffix.Data()));

    canvasDirGamma->cd();
    dummDirGammaRedX->DrawCopy();

        labelScalingDirGammaComb0010_2->Draw();
        labelScalingDirGammaComb2040_2->Draw();

        if (graphCombDirGammaSpectrumSystErr0010Plot2){
            graphCombDirGammaSpectrumSystErr0010Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot2){
            graphCombDirGammaSpectrumStatErr0010Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSystErr2040Plot2){
            graphCombDirGammaSpectrumSystErr2040Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot2){
            graphCombDirGammaSpectrumStatErr2040Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot2){
            graphCombDirGammaSpectrumSumErr2040ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot2);
        }
        if (graphCombDirGammaSpectrumSystErr2050Plot2){
            graphCombDirGammaSpectrumSystErr2050Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2050Plot2){
            graphCombDirGammaSpectrumStatErr2050Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2050ArPlot2){
            graphCombDirGammaSpectrumSumErr2050ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2050ArPlot2);
        }

        TLatex *labelDirGammaCombWithFit = new TLatex(0.62,0.93,"This thesis");
        SetStyleTLatex( labelDirGammaCombWithFit, 0.85*textsizeLabelsDirGamma,4);
        labelDirGammaCombWithFit->Draw();
        TLegend* legendDirGammaCombRedXWithFit = new TLegend(0.62,0.92-1*0.85*textsizeLabelsDirGamma*4,0.96,0.92);
        legendDirGammaCombRedXWithFit->SetFillStyle(0);
        legendDirGammaCombRedXWithFit->SetFillColor(0);
        legendDirGammaCombRedXWithFit->SetLineColor(0);
        legendDirGammaCombRedXWithFit->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaCombRedXWithFit->SetMargin(0.2);
        legendDirGammaCombRedXWithFit->SetTextFont(42);
        legendDirGammaCombRedXWithFit->SetHeader(collisionSystem.Data());
        legendDirGammaCombRedXWithFit->AddEntry(graphPCMDirGammaSpectrumSystErr0010Plot,"  0-10%","pf");
        legendDirGammaCombRedXWithFit->AddEntry(graphPCMDirGammaSpectrumSystErr2040Plot,"20-40%","pf");
        legendDirGammaCombRedXWithFit->AddEntry(graphPCMDirGammaSpectrumSystErr2050Plot,"20-50%","pf");
        legendDirGammaCombRedXWithFit->Draw();

        TLegend* legendDirGammaExpFit = new TLegend(0.2,0.1,0.6,0.1+1.*0.83*textsizeLabelsDirGamma*5);
        legendDirGammaExpFit->SetFillStyle(0);
        legendDirGammaExpFit->SetFillColor(0);
        legendDirGammaExpFit->SetLineColor(0);
        legendDirGammaExpFit->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaExpFit->SetMargin(0.15);
        legendDirGammaExpFit->SetTextFont(42);
        legendDirGammaExpFit->SetHeader("#it{A} exp(-#it{p}_{T}/#it{T}_{eff})");
        legendDirGammaExpFit->AddEntry(histoFitThermalGamma0010Stat,cent0010.Data(),"l");
        legendDirGammaExpFit->AddEntry((TObject*)0,Form("#it{T}_{eff} = %3.f #pm %3.f^{#it{stat}} #pm %3.f^{#it{syst}} MeV",1e3*fitThermalGamma0010Stat->GetParameter(1),1e3*fitThermalGamma0010Stat->GetParError(1),1e3*fitThermalGamma0010Sys->GetParError(1)),"");
        legendDirGammaExpFit->AddEntry(histoFitThermalGamma2040Stat,cent2040.Data(),"l");
        legendDirGammaExpFit->AddEntry((TObject*)0,Form("#it{T}_{eff} = %3.f #pm %3.f^{#it{stat}} #pm %3.f^{#it{syst}} MeV",1e3*fitThermalGamma2040Stat->GetParameter(1),1e3*fitThermalGamma2040Stat->GetParError(1),1e3*fitThermalGamma2040Sys->GetParError(1)),"");
        legendDirGammaExpFit->Draw();

        histoFitThermalGamma0010Stat->Draw("same,l");
        histoFitThermalGamma2040Stat->Draw("same,l");

    canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusFit_ReducedX.%s",outputDir.Data(),suffix.Data()));

    canvasDirGamma->cd();
    dummDirGammaRedX->DrawCopy();

        labelDirGammaCombWithFit->Draw();
        legendDirGammaCombRedXWithFit->Draw();
        labelScalingDirGammaComb0010_2->Draw();
        labelScalingDirGammaComb2040_2->Draw();

        if (graphCombDirGammaSpectrumSystErr0010Plot2){
            graphCombDirGammaSpectrumSystErr0010Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot2){
            graphCombDirGammaSpectrumStatErr0010Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSystErr2040Plot2){
            graphCombDirGammaSpectrumSystErr2040Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot2){
            graphCombDirGammaSpectrumStatErr2040Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot2){
            graphCombDirGammaSpectrumSumErr2040ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot2);
        }
        if (graphCombDirGammaSpectrumSystErr2050Plot2){
            graphCombDirGammaSpectrumSystErr2050Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2050Plot2){
            graphCombDirGammaSpectrumStatErr2050Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2050ArPlot2){
            graphCombDirGammaSpectrumSumErr2050ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2050ArPlot2);
        }
        TLegend* legendDirGammaTheoryPlusMcGillRedX_2 = new TLegend(0.19,0.35-1.*0.85*textsizeLabelsDirGamma*7.5,0.94,0.35);
        legendDirGammaTheoryPlusMcGillRedX_2->SetFillStyle(0);
        legendDirGammaTheoryPlusMcGillRedX_2->SetFillColor(0);
        legendDirGammaTheoryPlusMcGillRedX_2->SetLineColor(0);
        legendDirGammaTheoryPlusMcGillRedX_2->SetNColumns(2);
        legendDirGammaTheoryPlusMcGillRedX_2->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaTheoryPlusMcGillRedX_2->SetMargin(0.12);
        legendDirGammaTheoryPlusMcGillRedX_2->SetTextFont(42);
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry(graphTheoryMcGill0010Plot,"Paquet et al.","l");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"PRC 93 (2016) 044906","");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry(graphTheoryPHSDplusPrompt0010Plot,"Linnyk et al.","l");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"PRC 92 (2015) 054914 ","");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry(graphTheoryRapp0010Plot,"v. Hees et al.","l");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"NPA 933 (2015) 256","");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry(graphTheoryChatterjee0010Plot,"Chatterjee et al.","l");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"PRC 85 (2012) 064910","");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX_2->AddEntry((TObject*)0,"+ arXiv: 1305.0624","");
        legendDirGammaTheoryPlusMcGillRedX_2->Draw();

        graphTheoryMcGill0010Plot->Draw("p3lsame");
        graphTheoryPHSDplusPrompt0010Plot->Draw("p3lsame");
        graphTheoryRapp0010Plot->Draw("p3lsame");
        //graphTheoryChatterjee0010Plot->Draw("p3lsame");
        graphTheoryChatterjeeSummed0010Plot->Draw("p3lsame");
        if (graphCombDirGammaSpectrumSystErr0010Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr0010Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphCombDirGammaSpectrumStatErr0010Plot2->Draw("p,E1Z,same");
        }
        histoFitThermalGamma0010Stat->Draw("same,l");

        graphTheoryMcGill2040Plot->Draw("p3lsame");
        graphTheoryPHSD2040Plot->Draw("p3lsame");
        graphTheoryChatterjee2040Plot->Draw("p3lsame");
//         graphTheoryHees2040Plot->Draw("p3lsame");
        graphTheoryHe2040Plot->Draw("p3lsame");
        graphTheoryRapp2040Plot->Draw("p3lsame");

        if (graphCombDirGammaSpectrumSystErr2040Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot2, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr2040Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphCombDirGammaSpectrumStatErr2040Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , styleNLOMcGill, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErr2040ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot2);
        }
        histoFitThermalGamma2040Stat->Draw("same,l");

        graphTheoryChatterjee2050Plot->Draw("p3lsame");
        graphTheoryPHSDplusPrompt2050Plot->Draw("p3lsame");
        if (graphCombDirGammaSpectrumSystErr2050Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2050Plot2, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr2050Plot2->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2050Plot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2050Plot2, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphCombDirGammaSpectrumStatErr2050Plot2->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2050ArPlot2){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2050ArPlot2 , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErr2050ArPlot2->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2050ArPlot2);
        }

    canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusMcGill_ReducedX.%s",outputDir.Data(),suffix.Data()));



    cout << "Plotting direct gamma ALICE vs PHENIX at " << __LINE__ << endl;
    //*******************************************************************************************************
    //********************************PHENIX and ALICE together for 0-10% ***********************************
    //*******************************************************************************************************
//     graphPHENIXAuAuStat0020->RemovePoint(graphPHENIXAuAuStat0020->GetN()-1);
//     graphPHENIXAuAuSys0020->RemovePoint(graphPHENIXAuAuSys0020->GetN()-1);
    TF1* fitThermalGamma0020PHENIXStat                = FitObject("e","fitThermalGamma0020PHENIXStat","Photon",graphPHENIXAuAuStat0020,0.6,2,NULL,"QNRMEX0+");
    fitThermalGamma0020PHENIXStat->FixParameter(1,0.239);
    graphPHENIXAuAuStat0020->Fit(fitThermalGamma0020PHENIXStat,"QNRMEX0+","",0.6,2.);
    cout << WriteParameterToFile(fitThermalGamma0020PHENIXStat)<< endl;
    TF1* fitThermalGamma0020PHENIXSyst                = FitObject("e","fitThermalGamma0020PHENIXSyst","Photon",graphPHENIXAuAuSys0020,0.6,2,NULL,"QNRMEX0+");
    fitThermalGamma0020PHENIXSyst->FixParameter(1,0.239);
    graphPHENIXAuAuSys0020->Fit(fitThermalGamma0020PHENIXSyst,"QNRMEX0+","",0.6,2.);
    cout << WriteParameterToFile(fitThermalGamma0020PHENIXSyst)<< endl;

    TH1D* histoFitThermalGamma0020PHENIXStat        = (TH1D*)fitThermalGamma0020PHENIXStat->GetHistogram();
    SetStyleHisto(histoFitThermalGamma0020PHENIXStat, 3, styleFit, colorPHENIX );

    TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr0010Plot3 = NULL;
    if (graphCombDirGammaSpectrumSystErr0010Plot) {
        graphCombDirGammaSpectrumSystErr0010Plot3 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr0010Plot->Clone("graphCombDirGammaSpectrumSystErr0010Plot3");
        while (graphCombDirGammaSpectrumSystErr0010Plot3->GetX()[graphCombDirGammaSpectrumSystErr0010Plot3->GetN()-1] > 5.1 && graphCombDirGammaSpectrumSystErr0010Plot3->GetN()>0)
            graphCombDirGammaSpectrumSystErr0010Plot3->RemovePoint(graphCombDirGammaSpectrumSystErr0010Plot3->GetN()-1);
        if (graphCombDirGammaSpectrumSystErr0010Plot3->GetN()==0) graphCombDirGammaSpectrumSystErr0010Plot3 = NULL;
    }
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr0010Plot3 = NULL;
    if (graphCombDirGammaSpectrumStatErr0010Plot) {
        graphCombDirGammaSpectrumStatErr0010Plot3 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr0010Plot->Clone("graphCombDirGammaSpectrumStatErr0010Plot3");
        while (graphCombDirGammaSpectrumStatErr0010Plot3->GetX()[graphCombDirGammaSpectrumStatErr0010Plot3->GetN()-1] > 5.1 && graphCombDirGammaSpectrumStatErr0010Plot3->GetN()>0)
            graphCombDirGammaSpectrumStatErr0010Plot3->RemovePoint(graphCombDirGammaSpectrumStatErr0010Plot3->GetN()-1);
        if (graphCombDirGammaSpectrumStatErr0010Plot3->GetN()==0) graphCombDirGammaSpectrumStatErr0010Plot3 = NULL;
    }

//     cout << "here" << endl;
//     TGraphAsymmErrors* graphCombThermalGammaSpectrumSystErr0010Plot2 = NULL;
//     if (graphCombThermalGammaSpectrumSysErr0010) {
//         graphCombThermalGammaSpectrumSystErr0010Plot2 = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumSysErr0010->Clone("graphCombThermalGammaSpectrumSystErr0010Plot2");
//         while (graphCombThermalGammaSpectrumSystErr0010Plot2->GetX()[graphCombThermalGammaSpectrumSystErr0010Plot2->GetN()-1] > 4.1 && graphCombThermalGammaSpectrumSystErr0010Plot2->GetN()>0)
//             graphCombThermalGammaSpectrumSystErr0010Plot2->RemovePoint(graphCombThermalGammaSpectrumSystErr0010Plot2->GetN()-1);
//         if (graphCombThermalGammaSpectrumSystErr0010Plot2->GetN()==0) graphCombThermalGammaSpectrumSystErr0010Plot2 = NULL;
//     }
//     cout << "here" << endl;
//     TGraphAsymmErrors* graphCombThermalGammaSpectrumStatErr0010Plot2 = NULL;
//     if (graphCombThermalGammaSpectrumStatErr0010) {
//         graphCombThermalGammaSpectrumStatErr0010Plot2 = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumStatErr0010->Clone("graphCombThermalGammaSpectrumStatErr0010Plot2");
//         ProduceGraphAsymmWithoutXErrors(graphCombThermalGammaSpectrumStatErr0010Plot2);
//         while (graphCombThermalGammaSpectrumStatErr0010Plot2->GetX()[graphCombThermalGammaSpectrumStatErr0010Plot2->GetN()-1] > 4.1 && graphCombThermalGammaSpectrumStatErr0010Plot2->GetN()>0)
//             graphCombThermalGammaSpectrumStatErr0010Plot2->RemovePoint(graphCombThermalGammaSpectrumStatErr0010Plot2->GetN()-1);
//         if (graphCombThermalGammaSpectrumStatErr0010Plot2->GetN()==0) graphCombThermalGammaSpectrumStatErr0010Plot2 = NULL;
//     }
//     cout << "here" << endl;

    TH1D *dummDirGammaPHENIX0020 = new TH1D("dummDirGammaPHENIX0020", "dummDirGammaPHENIX0020", 1000, 0., 5.6);
    SetStyleHistoTH1ForGraphs( dummDirGammaPHENIX0020, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})", 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummDirGammaPHENIX0020->GetXaxis()->SetLabelOffset(-0.002);
    dummDirGammaPHENIX0020->GetXaxis()->SetTickLength(0.025);
    dummDirGammaPHENIX0020->GetYaxis()->SetTickLength(0.025);
    dummDirGammaPHENIX0020->GetYaxis()->SetRangeUser(0.9e-5,50);
    dummDirGammaPHENIX0020->GetXaxis()->SetRangeUser(0.,incSpectraX[1]);
    dummDirGammaPHENIX0020->DrawCopy();
    dummDirGammaPHENIX0020->DrawCopy("same,axis");

        DrawGammaSetMarkerTGraphAsym(graphPHENIXAuAuSys0020, markerStylePHENIX, markerSizePHENIX, colorPHENIX , colorPHENIX, widthLinesBoxes,  kTRUE);
        graphPHENIXAuAuSys0020->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPHENIXAuAuStat0020, markerStylePHENIX, markerSizePHENIX, colorPHENIX , colorPHENIX);
        graphPHENIXAuAuStat0020->Draw("p,E1Z,same");
        histoFitThermalGamma0020PHENIXStat->Draw("same,l");

        if (histoFitThermalGamma0010Stat){
            SetStyleHisto(histoFitThermalGamma0010Stat, 3, styleFit, colorComb0010+1 );
            histoFitThermalGamma0010Stat->Scale(1e-2);
            histoFitThermalGamma0010Stat->Draw("same,l");
        }
        if (graphCombDirGammaSpectrumSystErr0010Plot3){
            graphCombDirGammaSpectrumSystErr0010Plot3= ScaleGraph(graphCombDirGammaSpectrumSystErr0010Plot3,1e-2);
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0010Plot3, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr0010Plot3->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot3){
            graphCombDirGammaSpectrumStatErr0010Plot3= ScaleGraph(graphCombDirGammaSpectrumStatErr0010Plot3,1e-2);
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0010Plot3, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphCombDirGammaSpectrumStatErr0010Plot3->Draw("p,E1Z,same");
        }

        TLegend* legendDirGammaWithPHENIX0020 = new TLegend(0.455,0.935-1.15*0.85*textsizeLabelsDirGamma*8,0.555+0.21,0.935);
        legendDirGammaWithPHENIX0020->SetFillStyle(0);
        legendDirGammaWithPHENIX0020->SetFillColor(0);
        legendDirGammaWithPHENIX0020->SetLineColor(0);
        legendDirGammaWithPHENIX0020->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaWithPHENIX0020->SetMargin(0.2);
        legendDirGammaWithPHENIX0020->SetTextFont(42);
        legendDirGammaWithPHENIX0020->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot3,"ALICE","pf");
        legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,Form("%s",collisionSystemCent0010.Data()),"");
        legendDirGammaWithPHENIX0020->AddEntry(histoFitThermalGamma0010Stat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
        legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,Form("#it{T}_{eff} = %3.f #pm %3.f^{#it{stat}} #pm %3.f^{#it{syst}} MeV",1e3*fitThermalGamma0010Stat->GetParameter(1),1e3*fitThermalGamma0010Stat->GetParError(1),1e3*fitThermalGamma0010Sys->GetParError(1)),"");
        legendDirGammaWithPHENIX0020->AddEntry(graphPHENIXAuAuSys0020,"PHENIX","pf");
        legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,Form("%s",collisionSystemRHIC0010.Data()),"");
        legendDirGammaWithPHENIX0020->AddEntry(histoFitThermalGamma0020PHENIXStat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
//         legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,Form("#it{T}_{eff} = %3.f #pm %3.f^{#it{stat}} #pm %3.f^{#it{syst}} MeV",1e3*fitThermalGamma0020PHENIXStat->GetParameter(1),1e3*fitThermalGamma0020PHENIXStat->GetParError(1),1e3*fitThermalGamma0020PHENIXSyst->GetParError(1)),"");
//         TLegend* legendDirGammaWithPHENIX0020 = new TLegend(0.455,0.935-1.15*0.85*textsizeLabelsDirGamma*8,0.555+0.21,0.935);
//         legendDirGammaWithPHENIX0020->SetFillStyle(0);
//         legendDirGammaWithPHENIX0020->SetFillColor(0);
//         legendDirGammaWithPHENIX0020->SetLineColor(0);
//         legendDirGammaWithPHENIX0020->SetTextSize(0.85*textsizeLabelsDirGamma);
//         legendDirGammaWithPHENIX0020->SetMargin(0.2);
//         legendDirGammaWithPHENIX0020->SetTextFont(42);
//         legendDirGammaWithPHENIX0020->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot3,"ALICE","pf");
//         legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,Form("%s",collisionSystemCent0010.Data()),"");
//         legendDirGammaWithPHENIX0020->AddEntry(histoFitThermalGamma0010Stat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
//     //     legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,"#it{T}_{eff} = 297 #pm 12^{#it{stat}} #pm 41^{#it{sys}}",""); // numbers for thermal
//         legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,"#it{T}_{eff} = 304 #pm 11^{ #it{stat}} #pm 40^{#it{sys}} MeV","");
//         legendDirGammaWithPHENIX0020->AddEntry(graphPHENIXAuAuSys0020,"PHENIX","pf");
//         legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,Form("%s",collisionSystemRHIC0010.Data()),"");
//         legendDirGammaWithPHENIX0020->AddEntry(histoFitThermalGamma0020PHENIXStat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
        legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,"#it{T}_{eff} = 239 #pm 25^{#it{stat}} #pm 7^{ #it{sys}} MeV","");
        legendDirGammaWithPHENIX0020->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusPHENIX0020_ReducedX.%s",outputDir.Data(),suffix.Data()));

    //*******************************************************************************************************
    //********************************PHENIX and ALICE together for 0-10% ***********************************
    //*******************************************************************************************************
//     graphPHENIXAuAuStat2040->RemovePoint(4);
//     graphPHENIXAuAuStat2040->RemovePoint(graphPHENIXAuAuStat2040->GetN()-1);
//     graphPHENIXAuAuSys2040->RemovePoint(graphPHENIXAuAuSys2040->GetN()-1);
        TF1* fitThermalGamma2040PHENIXStat                = FitObject("e","fitThermalGamma2040PHENIXStat","Photon",graphPHENIXAuAuStat2040,0.6,2,NULL,"QNRMEX0+");
        fitThermalGamma2040PHENIXStat->FixParameter(1,0.260);
        fitThermalGamma2040PHENIXStat->SetParLimits(0,0,50);
        graphPHENIXAuAuStat2040->Fit(fitThermalGamma2040PHENIXStat,"QNRMEX0+","",0.6,2.);
        cout << WriteParameterToFile(fitThermalGamma2040PHENIXStat)<< endl;
        TF1* fitThermalGamma2040PHENIXSyst                = FitObject("e","fitThermalGamma2040PHENIXSyst","Photon",graphPHENIXAuAuSys2040,0.6,2,NULL,"QNRMEX0+");
        fitThermalGamma2040PHENIXSyst->FixParameter(1,0.260);
        fitThermalGamma2040PHENIXSyst->SetParLimits(0,0,50);
        graphPHENIXAuAuSys2040->Fit(fitThermalGamma2040PHENIXSyst,"QNRMEX0+","",0.6,2.);
        cout << WriteParameterToFile(fitThermalGamma2040PHENIXSyst)<< endl;

        TH1D* histoFitThermalGamma2040PHENIXStat        = (TH1D*)fitThermalGamma2040PHENIXStat->GetHistogram();
        SetStyleHisto(histoFitThermalGamma2040PHENIXStat, 3, styleFit, colorPHENIX );

        TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr2040Plot3 = NULL;
        if (graphCombDirGammaSpectrumSystErr2040Plot) {
            graphCombDirGammaSpectrumSystErr2040Plot3 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr2040Plot->Clone("graphCombDirGammaSpectrumSystErr2040Plot3");
            while (graphCombDirGammaSpectrumSystErr2040Plot3->GetX()[graphCombDirGammaSpectrumSystErr2040Plot3->GetN()-1] > 5.1 && graphCombDirGammaSpectrumSystErr2040Plot3->GetN()>0)
                graphCombDirGammaSpectrumSystErr2040Plot3->RemovePoint(graphCombDirGammaSpectrumSystErr2040Plot3->GetN()-1);
            if (graphCombDirGammaSpectrumSystErr2040Plot3->GetN()==0) graphCombDirGammaSpectrumSystErr2040Plot3 = NULL;
        }
        TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr2040Plot3 = NULL;
        if (graphCombDirGammaSpectrumStatErr2040Plot) {
            graphCombDirGammaSpectrumStatErr2040Plot3 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr2040Plot->Clone("graphCombDirGammaSpectrumStatErr2040Plot3");
            while (graphCombDirGammaSpectrumStatErr2040Plot3->GetX()[graphCombDirGammaSpectrumStatErr2040Plot3->GetN()-1] > 5.1 && graphCombDirGammaSpectrumStatErr2040Plot3->GetN()>0)
                graphCombDirGammaSpectrumStatErr2040Plot3->RemovePoint(graphCombDirGammaSpectrumStatErr2040Plot3->GetN()-1);
            if (graphCombDirGammaSpectrumStatErr2040Plot3->GetN()==0) graphCombDirGammaSpectrumStatErr2040Plot3 = NULL;
        }

        TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2040ArPlot3 = NULL;
        if (graphCombDirGammaSpectrumSumErr2040ArPlot) {
            graphCombDirGammaSpectrumSumErr2040ArPlot3 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSumErr2040ArPlot->Clone("graphCombDirGammaSpectrumSumErr2040ArPlot3");
            while (graphCombDirGammaSpectrumSumErr2040ArPlot3->GetX()[graphCombDirGammaSpectrumSumErr2040ArPlot3->GetN()-1] > 5.1 && graphCombDirGammaSpectrumSumErr2040ArPlot3->GetN()>0)
                graphCombDirGammaSpectrumSumErr2040ArPlot3->RemovePoint(graphCombDirGammaSpectrumSumErr2040ArPlot3->GetN()-1);
            if (graphCombDirGammaSpectrumSumErr2040ArPlot3->GetN()==0) graphCombDirGammaSpectrumSumErr2040ArPlot3 = NULL;
        }

//     cout << "here" << endl;
//     TGraphAsymmErrors* graphCombThermalGammaSpectrumSystErr2040Plot2 = NULL;
//     if (graphCombThermalGammaSpectrumSysErr2040) {
//         graphCombThermalGammaSpectrumSystErr2040Plot2 = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumSysErr2040->Clone("graphCombThermalGammaSpectrumSystErr2040Plot2");
//         while (graphCombThermalGammaSpectrumSystErr2040Plot2->GetX()[graphCombThermalGammaSpectrumSystErr2040Plot2->GetN()-1] > 4.1 && graphCombThermalGammaSpectrumSystErr2040Plot2->GetN()>0)
//             graphCombThermalGammaSpectrumSystErr2040Plot2->RemovePoint(graphCombThermalGammaSpectrumSystErr2040Plot2->GetN()-1);
//         if (graphCombThermalGammaSpectrumSystErr2040Plot2->GetN()==0) graphCombThermalGammaSpectrumSystErr2040Plot2 = NULL;
//     }
//     cout << "here" << endl;
//     TGraphAsymmErrors* graphCombThermalGammaSpectrumStatErr2040Plot2 = NULL;
//     if (graphCombThermalGammaSpectrumStatErr2040) {
//         graphCombThermalGammaSpectrumStatErr2040Plot2 = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumStatErr2040->Clone("graphCombThermalGammaSpectrumStatErr2040Plot2");
//         ProduceGraphAsymmWithoutXErrors(graphCombThermalGammaSpectrumStatErr2040Plot2);
//         while (graphCombThermalGammaSpectrumStatErr2040Plot2->GetX()[graphCombThermalGammaSpectrumStatErr2040Plot2->GetN()-1] > 4.1 && graphCombThermalGammaSpectrumStatErr2040Plot2->GetN()>0)
//             graphCombThermalGammaSpectrumStatErr2040Plot2->RemovePoint(graphCombThermalGammaSpectrumStatErr2040Plot2->GetN()-1);
//         if (graphCombThermalGammaSpectrumStatErr2040Plot2->GetN()==0) graphCombThermalGammaSpectrumStatErr2040Plot2 = NULL;
//     }
//     cout << "here" << endl;
    canvasDirGamma->cd();
    TH2D *dummDirGammaPHENIX2040 ;
    dummDirGammaPHENIX2040 = new TH2D("dummDirGammaPHENIX2040", "dummDirGammaPHENIX2040", 100, 0., 5.6,100000,2e-5,20);
    SetStyleHistoTH2ForGraphs( dummDirGammaPHENIX2040, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                            0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummDirGammaPHENIX2040->GetXaxis()->SetLabelOffset(-0.002);
    dummDirGammaPHENIX2040->GetXaxis()->SetTickLength(0.025);
    dummDirGammaPHENIX2040->GetYaxis()->SetTickLength(0.025);
    dummDirGammaPHENIX2040->GetYaxis()->SetRangeUser(1.1e-5,20);
    dummDirGammaPHENIX2040->GetXaxis()->SetRangeUser(0.,incSpectraX[1]);
    dummDirGammaPHENIX2040->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphPHENIXAuAuSys2040, markerStylePHENIX, markerSizePHENIX, colorPHENIX , colorPHENIX, widthLinesBoxes,  kTRUE);
        graphPHENIXAuAuSys2040->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPHENIXAuAuStat2040, markerStylePHENIX, markerSizePHENIX, colorPHENIX , colorPHENIX);
        graphPHENIXAuAuStat2040->Draw("p,E1Z,same");
        histoFitThermalGamma2040PHENIXStat->Draw("same,l");

        if (histoFitThermalGamma2040Stat){
            SetStyleHisto(histoFitThermalGamma2040Stat, 3, styleFit, colorComb2040+1 );
            histoFitThermalGamma2040Stat->Scale(1e-1);
            histoFitThermalGamma2040Stat->Draw("same,l");
        }

        if (graphCombDirGammaSpectrumSystErr2040Plot3){
            graphCombDirGammaSpectrumSystErr2040Plot3= ScaleGraph(graphCombDirGammaSpectrumSystErr2040Plot3,1e-1);
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot3, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr2040Plot3->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot3){
            graphCombDirGammaSpectrumStatErr2040Plot3= ScaleGraph(graphCombDirGammaSpectrumStatErr2040Plot3,1e-1);
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot3, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphCombDirGammaSpectrumStatErr2040Plot3->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot3){
            graphCombDirGammaSpectrumSumErr2040ArPlot3= ScaleGraph(graphCombDirGammaSpectrumSumErr2040ArPlot3,1e-1);
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot3, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphCombDirGammaSpectrumSumErr2040ArPlot3->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot3);
        }

        TLegend* legendDirGammaWithPHENIX2040 = new TLegend(0.455,0.935-1.15*0.85*textsizeLabelsDirGamma*8,0.555+0.21,0.935);
        legendDirGammaWithPHENIX2040->SetFillStyle(0);
        legendDirGammaWithPHENIX2040->SetFillColor(0);
        legendDirGammaWithPHENIX2040->SetLineColor(0);
        legendDirGammaWithPHENIX2040->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaWithPHENIX2040->SetMargin(0.2);
        legendDirGammaWithPHENIX2040->SetTextFont(42);
        legendDirGammaWithPHENIX2040->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot3,"ALICE","pf");
        legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,Form("%s",collisionSystemCent2040.Data()),"");
        legendDirGammaWithPHENIX2040->AddEntry(histoFitThermalGamma2040Stat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
        legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,Form("#it{T}_{eff} = %3.f #pm %3.f^{#it{stat}} #pm %3.f^{#it{syst}} MeV",1e3*fitThermalGamma2040Stat->GetParameter(1),1e3*fitThermalGamma2040Stat->GetParError(1),1e3*fitThermalGamma2040Sys->GetParError(1)),"");
        legendDirGammaWithPHENIX2040->AddEntry(graphPHENIXAuAuSys2040,"PHENIX","pf");
        legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,Form("%s",collisionSystemRHIC2040.Data()),"");
        legendDirGammaWithPHENIX2040->AddEntry(histoFitThermalGamma2040PHENIXStat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
//         legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,Form("#it{T}_{eff} = %3.f #pm %3.f^{#it{stat}} #pm %3.f^{#it{syst}} MeV",1e3*fitThermalGamma2040PHENIXStat->GetParameter(1),1e3*fitThermalGamma2040PHENIXStat->GetParError(1),1e3*fitThermalGamma2040PHENIXSyst->GetParError(1)),"");
//         TLegend* legendDirGammaWithPHENIX2040 = new TLegend(0.455,0.935-1.15*0.85*textsizeLabelsDirGamma*8,0.555+0.21,0.935);
//         legendDirGammaWithPHENIX2040->SetFillStyle(0);
//         legendDirGammaWithPHENIX2040->SetFillColor(0);
//         legendDirGammaWithPHENIX2040->SetLineColor(0);
//         legendDirGammaWithPHENIX2040->SetTextSize(0.85*textsizeLabelsDirGamma);
//         legendDirGammaWithPHENIX2040->SetMargin(0.2);
//         legendDirGammaWithPHENIX2040->SetTextFont(42);
//         legendDirGammaWithPHENIX2040->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot3,"ALICE","pf");
//         legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,Form("%s",collisionSystemCent2040.Data()),"");
//         legendDirGammaWithPHENIX2040->AddEntry(histoFitThermalGamma2040Stat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
//         legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,"#it{T}_{eff} = 407 #pm 61^{ #it{stat}} #pm 96^{#it{sys}} MeV","");
//         //     legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,"#it{T}_{eff} = 410 #pm 84^{#it{stat}} #pm 140^{#it{sys}}",""); // numbers for thermal
//         legendDirGammaWithPHENIX2040->AddEntry(graphPHENIXAuAuSys2040,"PHENIX","pf");
//         legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,Form("%s",collisionSystemRHIC2040.Data()),"");
//         legendDirGammaWithPHENIX2040->AddEntry(histoFitThermalGamma2040PHENIXStat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
        legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,"#it{T}_{eff} = 260 #pm 33^{#it{ stat}} #pm 8^{#it{sys}} MeV","");
        legendDirGammaWithPHENIX2040->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusPHENIX2040_ReducedX.%s",outputDir.Data(),suffix.Data()));

    //*******************************************************************************************************************************************
    //********************************************** Plotting direct Gamma Spectrum LinX ********************************************************
    //*******************************************************************************************************************************************
    TCanvas *canvasDirGammaLinX = new TCanvas("canvasDirGammaLinX","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasDirGammaLinX, 0.165, 0.01, 0.01, 0.07);
    canvasDirGammaLinX->SetLogy();

    TH2D *dummyDirGammaLinX ;
    dummyDirGammaLinX = new TH2D("dummyDirGammaLinX", "dummyDirGammaLinX", 1000, 0., 22, 1000., 6e-9,9e3);
    SetStyleHistoTH2ForGraphs( dummyDirGammaLinX, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                            0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummyDirGammaLinX->GetXaxis()->SetLabelOffset(-0.002);
    dummyDirGammaLinX->GetXaxis()->SetRangeUser(0,14.5);
    dummyDirGammaLinX->DrawCopy();

        if (graphCombDirGammaSpectrumSystErr0010Plot){
            graphCombDirGammaSpectrumSystErr0010Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot){
            graphCombDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }

        if (graphCombDirGammaSpectrumSystErr2040Plot){
            graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot){
            graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot){
            graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
        }

        if (graphCombDirGammaSpectrumSystErr2050Plot){
            graphCombDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2050Plot){
            graphCombDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2050ArPlot){
            graphCombDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2050ArPlot);
        }

        TLatex *labelScalingDirGamma0010_2 = new TLatex(12.5,6.E-4,"x 10^{2}");
        SetStyleTLatex( labelScalingDirGamma0010_2, 0.85*textsizeLabelsDirGamma,4,colorComb0010,42,kFALSE);
        labelScalingDirGamma0010_2->Draw();

        TLatex *labelScalingDirGamma2040_2 = new TLatex(12.5,3.5E-5,"x 10^{1}");
        SetStyleTLatex( labelScalingDirGamma2040_2, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
        labelScalingDirGamma2040_2->Draw();

        TLatex *labelScalingDirGamma2050_2 = new TLatex(12.5,6.5E-7,"x 10^{0}");
        SetStyleTLatex( labelScalingDirGamma2050_2, 0.85*textsizeLabelsDirGamma,4,colorComb2050,42,kFALSE);
//         labelScalingDirGamma2050_2->Draw();

        TLatex *labelDirGammaColl_2 = new TLatex(0.26,0.94,Form("%s",collisionSystem.Data()));
        SetStyleTLatex( labelDirGammaColl_2, 0.85*textsizeLabelsDirGamma,4);
        TLatex *labelDirGammaCollThesis_2 = new TLatex(0.22,0.26,"This thesis");
        SetStyleTLatex( labelDirGammaCollThesis_2, 0.85*textsizeLabelsDirGamma,4);
        labelDirGammaCollThesis_2->Draw();
        TLegend* legendDirGammaLinX = new TLegend(0.22,0.25-1.1*0.83*textsizeLabelsDirGamma*4,0.22+0.21,0.25);
        legendDirGammaLinX->SetFillStyle(0);
        legendDirGammaLinX->SetFillColor(0);
        legendDirGammaLinX->SetLineColor(0);
        legendDirGammaLinX->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaLinX->SetMargin(0.25);
        legendDirGammaLinX->SetTextFont(42);
        legendDirGammaLinX->SetHeader(collisionSystem.Data());
        legendDirGammaLinX->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot,"  0-10% ALICE","pf");
        legendDirGammaLinX->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
        legendDirGammaLinX->AddEntry(graphCombDirGammaSpectrumSystErr2050Plot,"20-50% ALICE","pf");
        legendDirGammaLinX->Draw();

    canvasDirGammaLinX->Print(Form("%s/DirGammaSpectrum_LinX_withoutFit.%s",outputDir.Data(),suffix.Data()));

        graphTheoryEPS090010Plot->Draw("p3lsame");
        graphTheoryCT100010Plot->Draw("p3lsame");
        graphTheoryNLO0010Plot->Draw("p3lsame");
        //graphTheoryPromptMcGill0020Plot->Draw("p3lsame");
        if (graphCombDirGammaSpectrumSystErr0010Plot){
            graphCombDirGammaSpectrumSystErr0010Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot){
            graphCombDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }

        graphTheoryEPS092040Plot->Draw("p3lsame");
        graphTheoryCT102040Plot->Draw("p3lsame");
        graphTheoryNLO2040Plot->Draw("p3lsame");
        graphTheoryPromptMcGill2040Plot->Draw("p3lsame");
        if (graphCombDirGammaSpectrumSystErr2040Plot){
            graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot){
            graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot){
            graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
        }

        graphTheoryEPS092050Plot->Draw("p3lsame");
        graphTheoryCT102050Plot->Draw("p3lsame");
        graphTheoryNLO2050Plot->Draw("p3lsame");

        if (graphCombDirGammaSpectrumSystErr2050Plot){
            graphCombDirGammaSpectrumSystErr2050Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr2050Plot){
            graphCombDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErr2050ArPlot){
            graphCombDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2050ArPlot);
        }

        TLegend* legendDirGammaTheoryLinX = new TLegend(0.5,0.93-1.*0.83*textsizeLabelsDirGamma*3,0.5+0.21,0.93);
        legendDirGammaTheoryLinX->SetFillStyle(0);
        legendDirGammaTheoryLinX->SetFillColor(0);
        legendDirGammaTheoryLinX->SetLineColor(0);
        legendDirGammaTheoryLinX->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaTheoryLinX->SetMargin(0.23);
        legendDirGammaTheoryLinX->SetTextFont(42);
        legendDirGammaTheoryLinX->AddEntry(graphTheoryNLO0010Plot,"PDF: CTEQ6M5, FF: GRV ","l");
        legendDirGammaTheoryLinX->AddEntry(graphTheoryPromptMcGill2040Plot,"(n)PDF: CTEQ6.1M/EPS09,","l");
        legendDirGammaTheoryLinX->AddEntry((TObject*)0,"FF: BFG2","");
        legendDirGammaTheoryLinX->Draw();

        TLegend* legendDirGammaTheoryLinX2 = new TLegend(0.5,0.93-1.*0.83*textsizeLabelsDirGamma*7-0.02,0.5+0.21,0.93-1.*0.85*textsizeLabelsDirGamma*3-0.02);
        legendDirGammaTheoryLinX2->SetFillStyle(0);
        legendDirGammaTheoryLinX2->SetFillColor(0);
        legendDirGammaTheoryLinX2->SetLineColor(0);
        legendDirGammaTheoryLinX2->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaTheoryLinX2->SetMargin(0.23);
        legendDirGammaTheoryLinX2->SetTextFont(42);
        legendDirGammaTheoryLinX2->AddEntry((TObject*)0,"#it{JETPHOX}","");
        legendDirGammaTheoryLinX2->AddEntry(graphTheoryCT100010Plot,"PDF: CT10, FF: BFG2","f");
        legendDirGammaTheoryLinX2->AddEntry(graphTheoryEPS090010Plot,"nPDF: EPS09, FF: BFG2","f");
        legendDirGammaTheoryLinX2->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
        legendDirGammaTheoryLinX2->Draw();

    canvasDirGammaLinX->Print(Form("%s/DirGammaSpectrum_LinX_withNLO_withoutFit.%s",outputDir.Data(),suffix.Data()));


    // *****************************************************************************************************
    // ***************************************** Plotting Comb RAA 0-10% ****************************************
    // *****************************************************************************************************
    cout << "Plotting Comb direct gamma RAA " << __LINE__ << endl;
    TGraphAsymmErrors* graphCombRAADirGammaSys0010Plot = NULL;
    if (graphCombRAADirGammaSys0010){
        graphCombRAADirGammaSys0010Plot = (TGraphAsymmErrors*)graphCombRAADirGammaSys0010->Clone("graphCombRAADirGammaSys0010Plot");
    }
    TGraphAsymmErrors* graphCombRAADirGammaStat0010Plot = NULL;
    if (graphCombRAADirGammaStat0010){
        graphCombRAADirGammaStat0010Plot = (TGraphAsymmErrors*)graphCombRAADirGammaStat0010->Clone("graphCombRAADirGammaStat0010Plot");
        ProduceGraphAsymmWithoutXErrors(graphCombRAADirGammaStat0010Plot);
    }

//     TCanvas* canvasRAA_0010 = new TCanvas("canvasRAA_0010","",200,10,1200,1100);  // gives the page size
//     DrawGammaCanvasSettings( canvasRAA_0010,  0.1, 0.01, 0.015, 0.1);
//     canvasRAA_0010->SetLogx(1);
//     Int_t textSizeLabelsPixelRAA = 48;
//     Double_t textsizeLabelsRAA = 0;
//     if (canvasRAA_0010->XtoPixel(canvasRAA_0010->GetX2()) < canvasRAA_0010->YtoPixel(canvasRAA_0010->GetY1())){
//         textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAA_0010->XtoPixel(canvasRAA_0010->GetX2()) ;
//     } else {
//         textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAA_0010->YtoPixel(canvasRAA_0010->GetY1());
//     }

    canvasRAA_0010->cd();

//     TH2F * histo2DRAADummy;
//     histo2DRAADummy = new TH2F("histo2DRAADummy","histo2DRAADummy",1000,0.,20.,1000,0.0,11.5);
//     SetStyleHistoTH2ForGraphs(histo2DRAADummy, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA,0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.95,1., 510, 510); //#frac{#frac{1
//     //     histo2DRAAAll3->GetYaxis()->SetRangeUser(0.05,8.);
//     histo2DRAADummy->GetYaxis()->SetLabelOffset(0.005);
//     histo2DRAADummy->GetXaxis()->SetLabelOffset(-0.005);
//     histo2DRAADummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
//     histo2DRAADummy->GetXaxis()->SetMoreLogLabels();
    histo2DRAADummy->DrawCopy("");

//         TLatex *labelRAAEnergy0010 = new TLatex(0.48,0.92,collisionSystemCent0010.Data());
//         SetStyleTLatex( labelRAAEnergy0010, 0.85*textsizeLabelsRAA,4);
        labelRAAEnergy0010->Draw();

//         TBox* boxErrorNorm0010_Single = CreateBoxConv(colorComb0010Box, 0.75, 1.-normErr0010 , 0.8, 1.+normErr0010);
        boxErrorNorm0010_Single->Draw();

        SetStyleGammaNLOTGraphWithBand( graphTheoryRelErrEPS090010, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryRelErrEPS090010->Draw("p3lsame");

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1] , 1, 1 ,3,kGray+1,7);

        if (graphCombRAADirGammaSys0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaSys0010Plot, markerStyleComb0010,markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphCombRAADirGammaSys0010Plot->Draw("E2same");
        }
        if (graphCombRAADirGammaStat0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaStat0010Plot, markerStyleComb0010,markerSizeComb0010, colorComb0010 , colorComb0010);
            graphCombRAADirGammaStat0010Plot->Draw("p,same,e1Z");
        }

        TLegend* legendRAA0010 = new TLegend(0.48,0.91-1.25*0.85*textsizeLabelsRAA*6,0.48+0.21,0.91);
        legendRAA0010->SetFillStyle(0);
        legendRAA0010->SetFillColor(0);
        legendRAA0010->SetLineColor(0);
        legendRAA0010->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAA0010->SetMargin(0.3);
        legendRAA0010->SetTextFont(42);
        if (graphCombRAADirGammaSys0010Plot)legendRAA0010->AddEntry(graphCombRAADirGammaSys0010Plot,"ALICE #it{#gamma}_{dir}","fp");
        legendRAA0010->AddEntry((TObject*)0,"pp reference: ","");
        legendRAA0010->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
        legendRAA0010->AddEntry((TObject*)0,"FF: BFG2","");
        legendRAA0010->AddEntry(graphTheoryRelErrEPS090010,"Rel. error #it{JETPHOX}","fp");
        legendRAA0010->AddEntry((TObject*)0,"nPDF: EPS09, FF: BFG2","");
        legendRAA0010->Draw();

        histo2DRAADummy->Draw("axis,same");

    canvasRAA_0010->Update();
    canvasRAA_0010->Print(Form("%s/DirGammaRAA_0010.%s",outputDir.Data(),suffix.Data()));

    // *****************************************************************************************************
    // ***************************************** Plotting RAA 20-40% ****************************************
    // *****************************************************************************************************
    TGraphAsymmErrors* graphCombRAADirGammaSys2040Plot = NULL;
    if (graphCombRAADirGammaSys2040){
        graphCombRAADirGammaSys2040Plot = (TGraphAsymmErrors*)graphCombRAADirGammaSys2040->Clone("graphCombRAADirGammaSys2040Plot");
    }
    TGraphAsymmErrors* graphCombRAADirGammaStat2040Plot = NULL;
    if (graphCombRAADirGammaStat2040){
        graphCombRAADirGammaStat2040Plot = (TGraphAsymmErrors*)graphCombRAADirGammaStat2040->Clone("graphCombRAADirGammaStat2040Plot");
        ProduceGraphAsymmWithoutXErrors(graphCombRAADirGammaStat2040Plot);
    }

//     TCanvas* canvasRAA_2040 = new TCanvas("canvasRAA_2040","",200,10,1200,1100);  // gives the page size
//     DrawGammaCanvasSettings( canvasRAA_2040,  0.1, 0.01, 0.015, 0.1);
//     canvasRAA_2040->SetLogx(1);
    canvasRAA_2040->cd();

    histo2DRAADummy->DrawCopy("");

//         TBox* boxErrorNorm2040_Single = CreateBoxConv(colorComb2040Box, 0.75, 1.-normErr2040 , 0.8, 1.+normErr2040);
        boxErrorNorm2040_Single->Draw();

//         TLatex *labelRAAEnergy2040 = new TLatex(0.48,0.92,collisionSystemCent2040.Data());
//         SetStyleTLatex( labelRAAEnergy2040, 0.85*textsizeLabelsRAA,4);
        labelRAAEnergy2040->Draw();

        SetStyleGammaNLOTGraphWithBand( graphTheoryRelErrEPS092040, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryRelErrEPS092040->Draw("p3lsame");

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1] , 1, 1 ,3,kGray+1,7);

        if (graphCombRAADirGammaSys2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaSys2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphCombRAADirGammaSys2040Plot->Draw("E2same");
        }
        if (graphCombRAADirGammaStat2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaStat2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphCombRAADirGammaStat2040Plot->Draw("p,E1Z,same");
        }
        if (graphCombRAADirGammaSum2040Ar){
//             graphCombRAADirGammaSum2040Ar->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaSum2040Ar , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
            for (Int_t i = 0; i < graphCombRAADirGammaSum2040Ar->GetN(); i++){
                graphCombRAADirGammaSum2040Ar->SetPointEYhigh(i,0);
            }
            graphCombRAADirGammaSum2040Ar->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombRAADirGammaSum2040Ar);
            graphCombRAADirGammaSum2040Ar->Print();
        }
        histo2DRAADummy->Draw("axis,same");

        TLegend* legendRAA2040 = new TLegend(0.48,0.91-1.25*0.85*textsizeLabelsRAA*6,0.48+0.21,0.91);
        legendRAA2040->SetFillStyle(0);
        legendRAA2040->SetFillColor(0);
        legendRAA2040->SetLineColor(0);
        legendRAA2040->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAA2040->SetMargin(0.3);
        legendRAA2040->SetTextFont(42);
        if (graphCombRAADirGammaSys2040Plot)legendRAA2040->AddEntry(graphCombRAADirGammaSys2040Plot,"ALICE #it{#gamma}_{dir}","fp");
        legendRAA2040->AddEntry((TObject*)0,"pp reference: ","");
        legendRAA2040->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
        legendRAA2040->AddEntry((TObject*)0,"FF: BFG2","");
        legendRAA2040->AddEntry(graphTheoryRelErrEPS092040,"Rel. error #it{JETPHOX}","fp");
        legendRAA2040->AddEntry((TObject*)0,"nPDF: EPS09, FF: BFG2","");
        legendRAA2040->Draw();


    canvasRAA_2040->Update();
    canvasRAA_2040->Print(Form("%s/DirGammaRAA_2040.%s",outputDir.Data(),suffix.Data()));

    // *****************************************************************************************************
    // ***************************************** Plotting RAA 20-50% ****************************************
    // *****************************************************************************************************
    TGraphAsymmErrors* graphCombRAADirGammaSys2050Plot = NULL;
    if (graphCombRAADirGammaSys2050){
        graphCombRAADirGammaSys2050Plot = (TGraphAsymmErrors*)graphCombRAADirGammaSys2050->Clone("graphCombRAADirGammaSys2050Plot");
    }
    TGraphAsymmErrors* graphCombRAADirGammaStat2050Plot = NULL;
    if (graphCombRAADirGammaStat2050){
        graphCombRAADirGammaStat2050Plot = (TGraphAsymmErrors*)graphCombRAADirGammaStat2050->Clone("graphCombRAADirGammaStat2050Plot");
        ProduceGraphAsymmWithoutXErrors(graphCombRAADirGammaStat2050Plot);
    }

//     TCanvas* canvasRAA_2050 = new TCanvas("canvasRAA_2050","",200,10,1200,1100);  // gives the page size
//     DrawGammaCanvasSettings( canvasRAA_2050,  0.1, 0.01, 0.015, 0.1);
//     canvasRAA_2050->SetLogx(1);
    canvasRAA_2050->cd();
    histo2DRAADummy->DrawCopy("");

//         TBox* boxErrorNorm2050_Single = CreateBoxConv(colorComb2050Box, 0.75, 1.-normErr2050 , 0.8, 1.+normErr2050);
        boxErrorNorm2050_Single->Draw();

//         TLatex *labelRAAEnergy2050 = new TLatex(0.48,0.92,collisionSystemCent2050.Data());
//         SetStyleTLatex( labelRAAEnergy2050, 0.85*textsizeLabelsRAA,4);
        labelRAAEnergy2050->Draw();

        SetStyleGammaNLOTGraphWithBand( graphTheoryRelErrEPS092050, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryRelErrEPS092050->Draw("p3lsame");

        DrawGammaLines(doubleRatioX[0], doubleRatioX[1] , 1, 1 ,3,kGray+1,7);

        if (graphCombRAADirGammaSys2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaSys2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
            graphCombRAADirGammaSys2050Plot->Draw("E2same");
        }
        if (graphCombRAADirGammaStat2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaStat2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphCombRAADirGammaStat2050Plot->Draw("p,E1Z,same");
        }
        if (graphCombRAADirGammaSum2050Ar){
//             graphCombRAADirGammaSum2050Ar->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaSum2050Ar , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            for (Int_t i = 0; i < graphCombRAADirGammaSum2050Ar->GetN(); i++){
                graphCombRAADirGammaSum2050Ar->SetPointEYhigh(i,0);
            }
            graphCombRAADirGammaSum2050Ar->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombRAADirGammaSum2050Ar);
            graphCombRAADirGammaSum2050Ar->Print();
        }
        histo2DRAADummy->Draw("axis,same");

        TLegend* legendRAA2050 = new TLegend(0.48,0.91-1.25*0.85*textsizeLabelsRAA*6,0.48+0.21,0.91);
        legendRAA2050->SetFillStyle(0);
        legendRAA2050->SetFillColor(0);
        legendRAA2050->SetLineColor(0);
        legendRAA2050->SetTextSize(0.85*textsizeLabelsRAA);
        legendRAA2050->SetMargin(0.3);
        legendRAA2050->SetTextFont(42);
        if (graphCombRAADirGammaSys2050Plot)legendRAA2050->AddEntry(graphCombRAADirGammaSys2050Plot,"ALICE #it{#gamma}_{dir}","fp");
        legendRAA2050->AddEntry((TObject*)0,"pp reference: ","");
        legendRAA2050->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
        legendRAA2050->AddEntry((TObject*)0,"FF: BFG2","");
        legendRAA2050->AddEntry(graphTheoryRelErrEPS092050,"Rel. error #it{JETPHOX}","fp");
        legendRAA2050->AddEntry((TObject*)0,"nPDF: EPS09, FF: BFG2","");
        legendRAA2050->Draw();


    canvasRAA_2050->Update();
    canvasRAA_2050->Print(Form("%s/DirGammaRAA_2050.%s",outputDir.Data(),suffix.Data()));

    //*******************************************************************************************************************
    //******************************************** Significance calculations ********************************************
    //*******************************************************************************************************************
    // originally written by Klaus Reygers
    cout << "Significance test at " << __LINE__ << endl;
    if (enablepValueCalc){

        //*********************************************************************************
        // Calculating Significance of combined double ratio for 0-10%
        //*********************************************************************************
        Int_t minBinSig0010 = 0;
        Int_t maxBinSig0010 = 3;
        const Int_t nPointSig0010 = maxBinSig0010- minBinSig0010 +1;

        // filling of graphs
        TGraph* graphAbsTypeAPlusStatErr0010                     = new TGraph(nPointSig0010);
        TGraph* graphRelTypeAPlusStatErr0010                     = new TGraph(nPointSig0010);
        TGraph* graphRelTypeBErr0010                             = new TGraph(nPointSig0010);
        TGraphErrors* graphDoubleRatioTypeAPlusStatErr0010         = new TGraphErrors(nPointSig0010);
        for (Int_t i = 0; i < nPointSig0010; i++){
            Double_t pt             = graphCombDRPi0FitSysBErr0010->GetX()[minBinSig0010+i];
            Double_t R                = graphCombDRPi0FitSysBErr0010->GetY()[minBinSig0010+i];
            Double_t statErr        = graphCombDRPi0FitStatErr0010->GetEYlow()[minBinSig0010+i];
            Double_t sysAErr        = graphCombDRPi0FitSysAErr0010->GetEYlow()[minBinSig0010+i];
            Double_t relErrSysA     = sysAErr/ R;
            Double_t statSysAErr    = TMath::Sqrt(statErr*statErr+sysAErr*sysAErr);
            Double_t relErrStatSysA = statSysAErr/ R;
            Double_t relErrSysB        = graphCombDRPi0FitSysBErr0010->GetEYlow()[minBinSig0010+i]/ graphCombDRPi0FitSysBErr0010->GetY()[minBinSig0010+i];

            graphAbsTypeAPlusStatErr0010->SetPoint(i, pt, statSysAErr);
            graphRelTypeAPlusStatErr0010->SetPoint(i, pt, relErrStatSysA);
            graphRelTypeBErr0010->SetPoint(i, pt, relErrSysB);
            graphDoubleRatioTypeAPlusStatErr0010->SetPoint(i, pt, R);
            graphDoubleRatioTypeAPlusStatErr0010->SetPointError(i,0,statSysAErr);
            cout << "pt: " << pt << endl;
        }
        Double_t relErrSysC0010 = graphCombDRPi0FitSysCErr0010->GetEYlow()[minBinSig0010]/ graphCombDRPi0FitSysCErr0010->GetY()[minBinSig0010];

        // Generating pseudo data
        TH1F* histoChi2NullHypo0010 = new TH1F("histoChi2NullHypo0010", "histoChi2NullHypo0010", 3000, 0., 1500);
        FillChi2HistForNullHypo(100000000, histoChi2NullHypo0010, graphRelTypeAPlusStatErr0010, graphRelTypeBErr0010, relErrSysC0010, kFALSE);

        // calculating significance
        Double_t chi2Data0010     = Chi2ForNullHypo(graphDoubleRatioTypeAPlusStatErr0010);
        Int_t binChi2Data0010     = histoChi2NullHypo0010->FindBin(chi2Data0010);
        Int_t binFirst0010         = 1;
        Int_t binLast0010        = histoChi2NullHypo0010->GetNbinsX();
        Double_t intTot0010        = histoChi2NullHypo0010->Integral(binFirst0010, binLast0010);
        Double_t intData0010    = histoChi2NullHypo0010->Integral(binChi2Data0010, binLast0010);
        Double_t pValue0010        = intData0010/ intTot0010;
        Double_t nSigma0010        = PValueToNSigma(pValue0010);

        cout << "0-10%: chi2data =" << chi2Data0010 << "\t, pVal = " << pValue0010 << "\t, nSigma = " << nSigma0010 << endl;

        // preparing histo for drawing
        TH1F* histoChi2NullData0010 = (TH1F*)histoChi2NullHypo0010->Clone("histoChi2NullData0010");
        for (Int_t i = 1; i < binChi2Data0010; i++ ){
            histoChi2NullData0010->SetBinContent(i,-1);
        }
        for (Int_t i = 1; i < binLast0010; i++ ){
            histoChi2NullData0010->SetBinError(i,0);
        }

        //*********************************************************************************
        // Calculating Significance of combined double ratio for 20-40%
        //*********************************************************************************
        Int_t minBinSig2040 = 0;
        Int_t maxBinSig2040 = 5;
        const Int_t nPointSig2040 = maxBinSig2040- minBinSig2040 +1;

        // filling of graphs
        TGraph* graphAbsTypeAPlusStatErr2040                     = new TGraph(nPointSig2040);
        TGraph* graphRelTypeAPlusStatErr2040                     = new TGraph(nPointSig2040);
        TGraph* graphRelTypeBErr2040                             = new TGraph(nPointSig2040);
        TGraphErrors* graphDoubleRatioTypeAPlusStatErr2040         = new TGraphErrors(nPointSig2040);

        for (Int_t i = 0; i < nPointSig2040; i++){
            Double_t pt             = graphCombDRPi0FitSysBErr2040->GetX()[minBinSig2040+i];
            Double_t R                = graphCombDRPi0FitSysBErr2040->GetY()[minBinSig2040+i];
            Double_t statErr        = graphCombDRPi0FitStatErr2040->GetEYlow()[minBinSig2040+i];
            Double_t sysAErr        = graphCombDRPi0FitSysAErr2040->GetEYlow()[minBinSig2040+i];
            Double_t relErrSysA     = sysAErr/ R;
            Double_t statSysAErr    = TMath::Sqrt(statErr*statErr+sysAErr*sysAErr);
            Double_t relErrStatSysA = statSysAErr/ R;
            Double_t relErrSysB        = graphCombDRPi0FitSysBErr2040->GetEYlow()[minBinSig2040+i]/ graphCombDRPi0FitSysBErr2040->GetY()[minBinSig2040+i];

            graphAbsTypeAPlusStatErr2040->SetPoint(i, pt, statSysAErr);
            graphRelTypeAPlusStatErr2040->SetPoint(i, pt, relErrStatSysA);
            graphRelTypeBErr2040->SetPoint(i, pt, relErrSysB);
            graphDoubleRatioTypeAPlusStatErr2040->SetPoint(i, pt, R);
            graphDoubleRatioTypeAPlusStatErr2040->SetPointError(i,0,statSysAErr);
        }
        Double_t relErrSysC2040 = graphCombDRPi0FitSysCErr2040->GetEYlow()[minBinSig2040]/ graphCombDRPi0FitSysCErr2040->GetY()[minBinSig2040];

        // Generating pseudo data
        TH1F* histoChi2NullHypo2040 = new TH1F("histoChi2NullHypo2040", "histoChi2NullHypo2040", 3000, 0., 1500);
        FillChi2HistForNullHypo(100000000, histoChi2NullHypo2040, graphRelTypeAPlusStatErr2040, graphRelTypeBErr2040, relErrSysC2040, kFALSE);

        // calculating significance
        Double_t chi2Data2040     = Chi2ForNullHypo(graphDoubleRatioTypeAPlusStatErr2040);
        Int_t binChi2Data2040     = histoChi2NullHypo2040->FindBin(chi2Data2040);
        Int_t binFirst2040         = 1;
        Int_t binLast2040        = histoChi2NullHypo2040->GetNbinsX();
        Double_t intTot2040        = histoChi2NullHypo2040->Integral(binFirst2040, binLast2040);
        Double_t intData2040    = histoChi2NullHypo2040->Integral(binChi2Data2040, binLast2040);
        Double_t pValue2040        = intData2040/ intTot2040;
        Double_t nSigma2040        = PValueToNSigma(pValue2040);

        cout << "20-40%:  chi2data =" << chi2Data2040 << "\t, pVal = " << pValue2040 << "\t, nSigma = " << nSigma2040 << endl;

        // preparing histo for drawing
        TH1F* histoChi2NullData2040 = (TH1F*)histoChi2NullHypo2040->Clone("histoChi2NullData2040");
        for (Int_t i = 1; i < binChi2Data2040; i++ ){
            histoChi2NullData2040->SetBinContent(i,-1);
        }
        for (Int_t i = 1; i < binLast2040; i++ ){
            histoChi2NullData2040->SetBinError(i,0);
        }

        //*********************************************************************************
        // Calculating Significance of combined double ratio for 20-50%
        //*********************************************************************************

        Int_t minBinSig2050 = 0;
        Int_t maxBinSig2050 = 5;
        const Int_t nPointSig2050 = maxBinSig2050- minBinSig2050 +1;

        // filling of graphs
        TGraph* graphAbsTypeAPlusStatErr2050                     = new TGraph(nPointSig2050);
        TGraph* graphRelTypeAPlusStatErr2050                     = new TGraph(nPointSig2050);
        TGraph* graphRelTypeBErr2050                             = new TGraph(nPointSig2050);
        TGraphErrors* graphDoubleRatioTypeAPlusStatErr2050         = new TGraphErrors(nPointSig2050);

        for (Int_t i = 0; i < nPointSig2050; i++){
            Double_t pt             = graphCombDRPi0FitSysBErr2050->GetX()[minBinSig2050+i];
            Double_t R                = graphCombDRPi0FitSysBErr2050->GetY()[minBinSig2050+i];
            Double_t statErr        = graphCombDRPi0FitStatErr2050->GetEYlow()[minBinSig2050+i];
            Double_t sysAErr        = graphCombDRPi0FitSysAErr2050->GetEYlow()[minBinSig2050+i];
            Double_t relErrSysA     = sysAErr/ R;
            Double_t statSysAErr    = TMath::Sqrt(statErr*statErr+sysAErr*sysAErr);
            Double_t relErrStatSysA = statSysAErr/ R;
            Double_t relErrSysB        = graphCombDRPi0FitSysBErr2050->GetEYlow()[minBinSig2050+i]/ graphCombDRPi0FitSysBErr2050->GetY()[minBinSig2050+i];

            graphAbsTypeAPlusStatErr2050->SetPoint(i, pt, statSysAErr);
            graphRelTypeAPlusStatErr2050->SetPoint(i, pt, relErrStatSysA);
            graphRelTypeBErr2050->SetPoint(i, pt, relErrSysB);
            graphDoubleRatioTypeAPlusStatErr2050->SetPoint(i, pt, R);
            graphDoubleRatioTypeAPlusStatErr2050->SetPointError(i,0,statSysAErr);
        }
        Double_t relErrSysC2050 = graphCombDRPi0FitSysCErr2050->GetEYlow()[minBinSig2050]/ graphCombDRPi0FitSysCErr2050->GetY()[minBinSig2050];

        // Generating pseudo data
        TH1F* histoChi2NullHypo2050 = new TH1F("histoChi2NullHypo2050", "histoChi2NullHypo2050", 3000, 0., 1500);
        FillChi2HistForNullHypo(100000000, histoChi2NullHypo2050, graphRelTypeAPlusStatErr2050, graphRelTypeBErr2050, relErrSysC2050, kFALSE);

        // calculating significance
        Double_t chi2Data2050     = Chi2ForNullHypo(graphDoubleRatioTypeAPlusStatErr2050);
        Int_t binChi2Data2050     = histoChi2NullHypo2050->FindBin(chi2Data2050);
        Int_t binFirst2050         = 1;
        Int_t binLast2050        = histoChi2NullHypo2050->GetNbinsX();
        Double_t intTot2050        = histoChi2NullHypo2050->Integral(binFirst2050, binLast2050);
        Double_t intData2050    = histoChi2NullHypo2050->Integral(binChi2Data2050, binLast2050);
        Double_t pValue2050        = intData2050/ intTot2050;
        Double_t nSigma2050        = PValueToNSigma(pValue2050);
        cout << "20-50%:  chi2data =" << chi2Data2050 << "\t, pVal = " << pValue2050 << "\t, nSigma = " << nSigma2050 << endl;

        // preparing histo for drawing
        TH1F* histoChi2NullData2050 = (TH1F*)histoChi2NullHypo2050->Clone("histoChi2NullData2050");
        for (Int_t i = 1; i < binChi2Data2050; i++ ){
            histoChi2NullData2050->SetBinContent(i,-1);
        }
        for (Int_t i = 1; i < binLast2050; i++ ){
            histoChi2NullData2050->SetBinError(i,0);
        }

        //*********************************************************************************
        // Plotting Significance of combined double ratio for 0-10%
        //*********************************************************************************

        TCanvas* canvasSignificance = new TCanvas("canvasSignificance","",200,10,1400,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasSignificance,  0.09, 0.01, 0.015, 0.1);
        canvasSignificance->SetLogy(1);

        Int_t textSizeLabelsPixelSignificance = 48;
        Double_t textsizeLabelsSignificance = 0;
        if (canvasSignificance->XtoPixel(canvasSignificance->GetX2()) < canvasSignificance->YtoPixel(canvasSignificance->GetY1())){
            textsizeLabelsSignificance = (Double_t)textSizeLabelsPixelSignificance/canvasSignificance->XtoPixel(canvasSignificance->GetX2()) ;
        } else {
            textsizeLabelsSignificance = (Double_t)textSizeLabelsPixelSignificance/canvasSignificance->YtoPixel(canvasSignificance->GetY1());
        }

        SetStyleHistoTH1ForGraphs(histoChi2NullHypo0010, "Test statistics #it{t}","Counts",
                                0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
                                0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
                                0.95, 1., 510, 510);
        histoChi2NullHypo0010->GetYaxis()->SetLabelOffset(0.005);

        Int_t firstAbove0010 = histoChi2NullHypo0010->FindFirstBinAbove();
        Int_t lastAbove0010 = histoChi2NullHypo0010->FindLastBinAbove();

        histoChi2NullHypo0010->GetXaxis()->SetRange(1,lastAbove0010+1);
        DrawGammaSetMarker(histoChi2NullHypo0010, 1, 0.1, kGray+2 , kGray+2);
        histoChi2NullHypo0010->DrawCopy("");

        DrawGammaSetMarker(histoChi2NullData0010, 1, 0.1, colorComb0010Box , colorComb0010Box);
        histoChi2NullData0010->SetFillColor(colorComb0010Box);
        histoChi2NullData0010->SetFillStyle(3356);
        histoChi2NullData0010->Draw("same,lf");
        histoChi2NullHypo0010->DrawCopy("same");
        histoChi2NullHypo0010->DrawCopy("same,axis");

            TLatex *labelSignificanceEnergy0010 = new TLatex(0.58,0.92,collisionSystemCent0010.Data());
            SetStyleTLatex( labelSignificanceEnergy0010, 0.85*textsizeLabelsSignificance,4);
            labelSignificanceEnergy0010->Draw();

            TLatex *labelSignificanceALICE = new TLatex(0.58,0.87,"ALICE pseudo data");
            SetStyleTLatex( labelSignificanceALICE, 0.85*textsizeLabelsSignificance,4);
            labelSignificanceALICE->Draw();

            DrawGammaLines(chi2Data0010, chi2Data0010, 0, histoChi2NullHypo0010->GetMaximum()*0.1, 3, colorComb0010, 7);

            TLatex *labelTData0010 = new TLatex(chi2Data0010+10,histoChi2NullHypo0010->GetMaximum()*0.1,"#it{t}_{data}");
            SetStyleTLatex( labelTData0010, 0.85*textsizeLabelsSignificance,4,colorComb0010,42,kFALSE);
            labelTData0010->Draw();

            TLatex *labelSignificancePValue0010 = new TLatex(0.48,0.57,Form("#it{p}-value = %1.4f",pValue0010));
            SetStyleTLatex( labelSignificancePValue0010, 0.85*textsizeLabelsSignificance,4);
            labelSignificancePValue0010->Draw();
            TLatex *labelSignificanceNSigma0010 = new TLatex(0.48,0.53,Form("%1.2f #sigma",nSigma0010));
            SetStyleTLatex( labelSignificanceNSigma0010, 0.85*textsizeLabelsSignificance,4);
            labelSignificanceNSigma0010->Draw();


        canvasSignificance->Update();
        canvasSignificance->Print(Form("%s/DirGammaSignificanceTest_0010.%s",outputDir.Data(),suffix.Data()));

        //*********************************************************************************
        // Plotting Significance of combined double ratio for 20-40%
        //*********************************************************************************

        SetStyleHistoTH1ForGraphs(histoChi2NullHypo2040, "Test statistics #it{t}","Counts",
                                0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
                                0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
                                0.95, 1., 510, 510);
        histoChi2NullHypo2040->GetYaxis()->SetLabelOffset(0.005);

        Int_t firstAbove2040 = histoChi2NullHypo2040->FindFirstBinAbove();
        Int_t lastAbove2040 = histoChi2NullHypo2040->FindLastBinAbove();

        histoChi2NullHypo2040->GetXaxis()->SetRange(1,lastAbove2040+1);
        DrawGammaSetMarker(histoChi2NullHypo2040, 1, 0.1, kGray+2 , kGray+2);
        histoChi2NullHypo2040->DrawCopy("");

        DrawGammaSetMarker(histoChi2NullData2040, 1, 0.1, colorComb2040Box , colorComb2040Box);
        histoChi2NullData2040->SetFillColor(colorComb2040Box);
        histoChi2NullData2040->SetFillStyle(3356);
        histoChi2NullData2040->Draw("same,lf");
        histoChi2NullHypo2040->DrawCopy("same");
        histoChi2NullHypo2040->DrawCopy("same,axis");

            TLatex *labelSignificanceEnergy2040 = new TLatex(0.58,0.92,collisionSystemCent2040.Data());
            SetStyleTLatex( labelSignificanceEnergy2040, 0.85*textsizeLabelsSignificance,4);
            labelSignificanceEnergy2040->Draw();

            labelSignificanceALICE->Draw();

            DrawGammaLines(chi2Data2040, chi2Data2040, 0, histoChi2NullHypo2040->GetMaximum()*0.5, 3, colorComb2040, 7);

            TLatex *labelTData2040 = new TLatex(chi2Data2040+10,histoChi2NullHypo2040->GetMaximum()*0.5,"#it{t}_{data}");
            SetStyleTLatex( labelTData2040, 0.85*textsizeLabelsSignificance,4,colorComb2040,42,kFALSE);
            labelTData2040->Draw();

            TLatex *labelSignificancePValue2040 = new TLatex(0.48,0.57,Form("#it{p}-value = %1.4f",pValue2040));
            SetStyleTLatex( labelSignificancePValue2040, 0.85*textsizeLabelsSignificance,4);
            labelSignificancePValue2040->Draw();
            TLatex *labelSignificanceNSigma2040 = new TLatex(0.48,0.53,Form("%1.2f #sigma",nSigma2040));
            SetStyleTLatex( labelSignificanceNSigma2040, 0.85*textsizeLabelsSignificance,4);
            labelSignificanceNSigma2040->Draw();


        canvasSignificance->Update();
        canvasSignificance->Print(Form("%s/DirGammaSignificanceTest_2040.%s",outputDir.Data(),suffix.Data()));

        //*********************************************************************************
        // Plotting Significance of combined double ratio for 20-50%
        //*********************************************************************************

        SetStyleHistoTH1ForGraphs(histoChi2NullHypo2050, "Test statistics #it{t}","Counts",
                                0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
                                0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
                                0.95, 1., 510, 510);
        histoChi2NullHypo2050->GetYaxis()->SetLabelOffset(0.005);

        Int_t firstAbove2050 = histoChi2NullHypo2050->FindFirstBinAbove();
        Int_t lastAbove2050 = histoChi2NullHypo2050->FindLastBinAbove();

        histoChi2NullHypo2050->GetXaxis()->SetRange(1,lastAbove2050+1);
        histoChi2NullHypo2050->GetYaxis()->SetRangeUser(1, 10*histoChi2NullHypo2050->GetMaximum());
        DrawGammaSetMarker(histoChi2NullHypo2050, 1, 0.1, kGray+2 , kGray+2);
        histoChi2NullHypo2050->DrawCopy("");

        DrawGammaSetMarker(histoChi2NullData2050, 1, 0.1, colorComb2050Box , colorComb2050Box);
        histoChi2NullData2050->SetFillColor(colorComb2050Box);
        histoChi2NullData2050->SetFillStyle(3356);
        histoChi2NullData2050->Draw("same,lf");
        histoChi2NullHypo2050->DrawCopy("same");
        histoChi2NullHypo2050->DrawCopy("same,axis");

            TLatex *labelSignificanceEnergy2050 = new TLatex(0.58,0.92,collisionSystemCent2050.Data());
            SetStyleTLatex( labelSignificanceEnergy2050, 0.85*textsizeLabelsSignificance,4);
            labelSignificanceEnergy2050->Draw();

            labelSignificanceALICE->Draw();

            DrawGammaLines(chi2Data2050, chi2Data2050, 0, histoChi2NullHypo2050->GetMaximum()*0.2, 3, colorComb2050, 7);

            TLatex *labelTData2050 = new TLatex(chi2Data2050+10,histoChi2NullHypo2050->GetMaximum()*0.2,"#it{t}_{data}");
            SetStyleTLatex( labelTData2050, 0.85*textsizeLabelsSignificance,4,colorComb2050,42,kFALSE);
            labelTData2050->Draw();

            TLatex *labelSignificancePValue2050 = new TLatex(0.48,0.57,Form("#it{p}-value = %1.4f",pValue2050));
            SetStyleTLatex( labelSignificancePValue2050, 0.85*textsizeLabelsSignificance,4);
            labelSignificancePValue2050->Draw();
            TLatex *labelSignificanceNSigma2050 = new TLatex(0.48,0.53,Form("%1.2f #sigma",nSigma2050));
            SetStyleTLatex( labelSignificanceNSigma2050, 0.85*textsizeLabelsSignificance,4);
            labelSignificanceNSigma2050->Draw();


        canvasSignificance->Update();
        canvasSignificance->Print(Form("%s/DirGammaSignificanceTest_2050.%s",outputDir.Data(),suffix.Data()));

    }

    // ************************************************************************************************************************************************
    // ************************************************************ Write data to output file *********************************************************
    // ************************************************************************************************************************************************
    const char* fileNameOutputComp = "Gamma_CombResults_PbPb_2.76TeV.root";
    TFile* fileGammaSpectrum = new TFile(fileNameOutputComp,"UPDATE");
        fileGammaSpectrum->mkdir("Gamma_PbPb_2.76TeV_0-10%");
        fileGammaSpectrum->cd("Gamma_PbPb_2.76TeV_0-10%");
            if (graphCombDRPi0FitStatErr0010) graphCombDRPi0FitStatErr0010->Write("DR_comb_StatErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysErr0010) graphCombDRPi0FitSysErr0010->Write("DR_comb_SysErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSumErr0010) graphCombDRPi0FitSumErr0010->Write("DR_comb_totErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysAErr0010) graphCombDRPi0FitSysAErr0010->Write("DR_comb_SysAErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysBErr0010) graphCombDRPi0FitSysBErr0010->Write("DR_comb_SysBErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysCErr0010) graphCombDRPi0FitSysCErr0010->Write("DR_comb_SysCErr",TObject::kOverwrite);

            if (graphCombIncGammaStatErr0010) graphCombIncGammaStatErr0010->Write("IncGammaSpec_comb_StatErr",TObject::kOverwrite);
            if (graphCombIncGammaSysErr0010) graphCombIncGammaSysErr0010->Write("IncGammaSpec_comb_SysErr",TObject::kOverwrite);
            if (graphCombIncGammaSumErr0010) graphCombIncGammaSumErr0010->Write("IncGammaSpec_comb_totErr",TObject::kOverwrite);
            if (graphCombIncGammaSysAErr0010) graphCombIncGammaSysAErr0010->Write("IncGammaSpec_comb_SysAErr",TObject::kOverwrite);
            if (graphCombIncGammaSysBErr0010) graphCombIncGammaSysBErr0010->Write("IncGammaSpec_comb_SysBErr",TObject::kOverwrite);
            if (graphCombIncGammaSysCErr0010) graphCombIncGammaSysCErr0010->Write("IncGammaSpec_comb_SysCErr",TObject::kOverwrite);

            if (graphCombDirGammaSpectrumStatErr0010) graphCombDirGammaSpectrumStatErr0010->Write("DirGammaSpec_comb_StatErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystErr0010) graphCombDirGammaSpectrumSystErr0010->Write("DirGammaSpec_comb_SysErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystAErr0010) graphCombDirGammaSpectrumSystAErr0010->Write("DirGammaSpec_comb_SysAErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystBErr0010) graphCombDirGammaSpectrumSystBErr0010->Write("DirGammaSpec_comb_SysBErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystCErr0010) graphCombDirGammaSpectrumSystCErr0010->Write("DirGammaSpec_comb_SysCErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSumErr0010) graphCombDirGammaSpectrumSumErr0010->Write("DirGammaSpec_comb_totErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSumErr0010Ar) graphCombDirGammaSpectrumSumErr0010Ar->Write("DirGammaSpec_comb_upperLimits",TObject::kOverwrite);

            if (graphCombThermalGammaSpectrumStatErr0010) graphCombThermalGammaSpectrumStatErr0010->Write("ThermalGammaSpec_comb_StatErr",TObject::kOverwrite);
            if (graphCombThermalGammaSpectrumSysErr0010) graphCombThermalGammaSpectrumSysErr0010->Write("ThermalGammaSpec_comb_SysErr",TObject::kOverwrite);
            if (fitPureThermalGamma0010Stat) fitPureThermalGamma0010Stat->Write("PureThermal_ExpFit_Stat",TObject::kOverwrite);

            if (histoPCMIncGammaStatErr0010) histoPCMIncGammaStatErr0010->Write("hIncGammaSpec_PCM_StatErr",TObject::kOverwrite);
            if (graphPCMIncGammaSysAErr0010) graphPCMIncGammaSysAErr0010->Write("IncGammaSpec_PCM_SysAErr",TObject::kOverwrite);
            if (graphPCMIncGammaSysBErr0010) graphPCMIncGammaSysBErr0010->Write("IncGammaSpec_PCM_SysBErr",TObject::kOverwrite);
            if (graphPCMIncGammaSysCErr0010) graphPCMIncGammaSysCErr0010->Write("IncGammaSpec_PCM_SysCErr",TObject::kOverwrite);

            if (histoPHOSIncGammaStatErr0010) histoPHOSIncGammaStatErr0010->Write("hIncGammaSpec_PHOS_StatErr",TObject::kOverwrite);
            if (graphPHOSIncGammaSysAErr0010) graphPHOSIncGammaSysAErr0010->Write("IncGammaSpec_PHOS_SysAErr",TObject::kOverwrite);
            if (graphPHOSIncGammaSysBErr0010) graphPHOSIncGammaSysBErr0010->Write("IncGammaSpec_PHOS_SysBErr",TObject::kOverwrite);
            if (graphPHOSIncGammaSysCErr0010) graphPHOSIncGammaSysCErr0010->Write("IncGammaSpec_PHOS_SysCErr",TObject::kOverwrite);

            if (histoPCMDRPi0FitStatErr0010) histoPCMDRPi0FitStatErr0010->Write("hDR_PCM_StatErr",TObject::kOverwrite);
            if (graphPCMDRPi0FitSysAErr0010) graphPCMDRPi0FitSysAErr0010->Write("DR_PCM_SysAErr",TObject::kOverwrite);
            if (graphPCMDRPi0FitSysBErr0010) graphPCMDRPi0FitSysBErr0010->Write("DR_PCM_SysBErr",TObject::kOverwrite);
            if (graphPCMDRPi0FitSysCErr0010) graphPCMDRPi0FitSysCErr0010->Write("DR_PCM_SysCErr",TObject::kOverwrite);

            if (histoPHOSDRPi0FitStatErr0010) histoPHOSDRPi0FitStatErr0010->Write("hDR_PHOS_StatErr",TObject::kOverwrite);
            if (graphPHOSDRPi0FitSysAErr0010) graphPHOSDRPi0FitSysAErr0010->Write("DR_PHOS_SysAErr",TObject::kOverwrite);
            if (graphPHOSDRPi0FitSysBErr0010) graphPHOSDRPi0FitSysBErr0010->Write("DR_PHOS_SysBErr",TObject::kOverwrite);
            if (graphPHOSDRPi0FitSysCErr0010) graphPHOSDRPi0FitSysCErr0010->Write("DR_PHOS_SysCErr",TObject::kOverwrite);

            if (graphCombRAADirGammaStat0010) graphCombRAADirGammaStat0010->Write("RAA_comb_StatErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSys0010) graphCombRAADirGammaSys0010->Write("RAA_comb_SysErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSum0010Ar) graphCombRAADirGammaSum0010Ar->Write("RAA_comb_upperLimits",TObject::kOverwrite);

        fileGammaSpectrum->mkdir("Gamma_PbPb_2.76TeV_20-40%");
        fileGammaSpectrum->cd("Gamma_PbPb_2.76TeV_20-40%");
            if (graphCombDRPi0FitStatErr2040) graphCombDRPi0FitStatErr2040->Write("DR_comb_StatErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysErr2040) graphCombDRPi0FitSysErr2040->Write("DR_comb_SysErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSumErr2040) graphCombDRPi0FitSumErr2040->Write("DR_comb_totErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysAErr2040) graphCombDRPi0FitSysAErr2040->Write("DR_comb_SysAErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysBErr2040) graphCombDRPi0FitSysBErr2040->Write("DR_comb_SysBErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysCErr2040) graphCombDRPi0FitSysCErr2040->Write("DR_comb_SysCErr",TObject::kOverwrite);

            if (graphCombIncGammaStatErr2040) graphCombIncGammaStatErr2040->Write("IncGammaSpec_comb_StatErr",TObject::kOverwrite);
            if (graphCombIncGammaSysErr2040) graphCombIncGammaSysErr2040->Write("IncGammaSpec_comb_SysErr",TObject::kOverwrite);
            if (graphCombIncGammaSumErr2040) graphCombIncGammaSumErr2040->Write("IncGammaSpec_comb_totErr",TObject::kOverwrite);
            if (graphCombIncGammaSysAErr2040) graphCombIncGammaSysAErr2040->Write("IncGammaSpec_comb_SysAErr",TObject::kOverwrite);
            if (graphCombIncGammaSysBErr2040) graphCombIncGammaSysBErr2040->Write("IncGammaSpec_comb_SysBErr",TObject::kOverwrite);
            if (graphCombIncGammaSysCErr2040) graphCombIncGammaSysCErr2040->Write("IncGammaSpec_comb_SysCErr",TObject::kOverwrite);

            if (graphCombDirGammaSpectrumStatErr2040) graphCombDirGammaSpectrumStatErr2040->Write("DirGammaSpec_comb_StatErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystErr2040) graphCombDirGammaSpectrumSystErr2040->Write("DirGammaSpec_comb_SysErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystAErr2040) graphCombDirGammaSpectrumSystAErr2040->Write("DirGammaSpec_comb_SysAErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystBErr2040) graphCombDirGammaSpectrumSystBErr2040->Write("DirGammaSpec_comb_SysBErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystCErr2040) graphCombDirGammaSpectrumSystCErr2040->Write("DirGammaSpec_comb_SysCErr",TObject::kOverwrite);

            if (graphCombDirGammaSpectrumSumErr2040) graphCombDirGammaSpectrumSumErr2040->Write("DirGammaSpec_comb_totErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSumErr2040Ar) graphCombDirGammaSpectrumSumErr2040Ar->Write("DirGammaSpec_comb_upperLimits",TObject::kOverwrite);

            if (histoPCMIncGammaStatErr2040) histoPCMIncGammaStatErr2040->Write("hIncGammaSpec_PCM_StatErr",TObject::kOverwrite);
            if (graphPCMIncGammaSysAErr2040) graphPCMIncGammaSysAErr2040->Write("IncGammaSpec_PCM_SysAErr",TObject::kOverwrite);
            if (graphPCMIncGammaSysBErr2040) graphPCMIncGammaSysBErr2040->Write("IncGammaSpec_PCM_SysBErr",TObject::kOverwrite);
            if (graphPCMIncGammaSysCErr2040) graphPCMIncGammaSysCErr2040->Write("IncGammaSpec_PCM_SysCErr",TObject::kOverwrite);

            if (histoPHOSIncGammaStatErr2040) histoPHOSIncGammaStatErr2040->Write("hIncGammaSpec_PHOS_StatErr",TObject::kOverwrite);
            if (graphPHOSIncGammaSysAErr2040) graphPHOSIncGammaSysAErr2040->Write("IncGammaSpec_PHOS_SysAErr",TObject::kOverwrite);
            if (graphPHOSIncGammaSysBErr2040) graphPHOSIncGammaSysBErr2040->Write("IncGammaSpec_PHOS_SysBErr",TObject::kOverwrite);
            if (graphPHOSIncGammaSysCErr2040) graphPHOSIncGammaSysCErr2040->Write("IncGammaSpec_PHOS_SysCErr",TObject::kOverwrite);

            if (histoPCMDRPi0FitStatErr2040) histoPCMDRPi0FitStatErr2040->Write("hDR_PCM_StatErr",TObject::kOverwrite);
            if (graphPCMDRPi0FitSysAErr2040) graphPCMDRPi0FitSysAErr2040->Write("DR_PCM_SysAErr",TObject::kOverwrite);
            if (graphPCMDRPi0FitSysBErr2040) graphPCMDRPi0FitSysBErr2040->Write("DR_PCM_SysBErr",TObject::kOverwrite);
            if (graphPCMDRPi0FitSysCErr2040) graphPCMDRPi0FitSysCErr2040->Write("DR_PCM_SysCErr",TObject::kOverwrite);

            if (histoPHOSDRPi0FitStatErr2040) histoPHOSDRPi0FitStatErr2040->Write("hDR_PHOS_StatErr",TObject::kOverwrite);
            if (graphPHOSDRPi0FitSysAErr2040) graphPHOSDRPi0FitSysAErr2040->Write("DR_PHOS_SysAErr",TObject::kOverwrite);
            if (graphPHOSDRPi0FitSysBErr2040) graphPHOSDRPi0FitSysBErr2040->Write("DR_PHOS_SysBErr",TObject::kOverwrite);
            if (graphPHOSDRPi0FitSysCErr2040) graphPHOSDRPi0FitSysCErr2040->Write("DR_PHOS_SysCErr",TObject::kOverwrite);

            if (graphCombRAADirGammaStat2040) graphCombRAADirGammaStat2040->Write("RAA_comb_StatErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSys2040) graphCombRAADirGammaSys2040->Write("RAA_comb_SysErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSum2040Ar) graphCombRAADirGammaSum2040Ar->Write("RAA_comb_upperLimits",TObject::kOverwrite);

        fileGammaSpectrum->mkdir("Gamma_PbPb_2.76TeV_20-50%");
        fileGammaSpectrum->cd("Gamma_PbPb_2.76TeV_20-50%");
            if (graphCombDRPi0FitStatErr2050) graphCombDRPi0FitStatErr2050->Write("DR_comb_StatErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysErr2050) graphCombDRPi0FitSysErr2050->Write("DR_comb_SysErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSumErr2050) graphCombDRPi0FitSumErr2050->Write("DR_comb_totErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysAErr2050) graphCombDRPi0FitSysAErr2050->Write("DR_comb_SysAErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysBErr2050) graphCombDRPi0FitSysBErr2050->Write("DR_comb_SysBErr",TObject::kOverwrite);
            if (graphCombDRPi0FitSysCErr2050) graphCombDRPi0FitSysCErr2050->Write("DR_comb_SysCErr",TObject::kOverwrite);

            if (graphCombIncGammaStatErr2050) graphCombIncGammaStatErr2050->Write("IncGammaSpec_comb_StatErr",TObject::kOverwrite);
            if (graphCombIncGammaSysErr2050) graphCombIncGammaSysErr2050->Write("IncGammaSpec_comb_SysErr",TObject::kOverwrite);
            if (graphCombIncGammaSumErr2050) graphCombIncGammaSumErr2050->Write("IncGammaSpec_comb_totErr",TObject::kOverwrite);
            if (graphCombIncGammaSysAErr2050) graphCombIncGammaSysAErr2050->Write("IncGammaSpec_comb_SysAErr",TObject::kOverwrite);
            if (graphCombIncGammaSysBErr2050) graphCombIncGammaSysBErr2050->Write("IncGammaSpec_comb_SysBErr",TObject::kOverwrite);
            if (graphCombIncGammaSysCErr2050) graphCombIncGammaSysCErr2050->Write("IncGammaSpec_comb_SysCErr",TObject::kOverwrite);

            if (graphCombDirGammaSpectrumStatErr2050) graphCombDirGammaSpectrumStatErr2050->Write("DirGammaSpec_comb_StatErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystErr2050) graphCombDirGammaSpectrumSystErr2050->Write("DirGammaSpec_comb_SysErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystAErr2050) graphCombDirGammaSpectrumSystAErr2050->Write("DirGammaSpec_comb_SysAErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystBErr2050) graphCombDirGammaSpectrumSystBErr2050->Write("DirGammaSpec_comb_SysBErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSystCErr2050) graphCombDirGammaSpectrumSystCErr2050->Write("DirGammaSpec_comb_SysCErr",TObject::kOverwrite);

            if (graphCombDirGammaSpectrumSumErr2050) graphCombDirGammaSpectrumSumErr2050->Write("DirGammaSpec_comb_totErr",TObject::kOverwrite);
            if (graphCombDirGammaSpectrumSumErr2050Ar) graphCombDirGammaSpectrumSumErr2050Ar->Write("DirGammaSpec_comb_upperLimits",TObject::kOverwrite);

            if (histoPCMIncGammaStatErr2050) histoPCMIncGammaStatErr2050->Write("hIncGammaSpec_PCM_StatErr",TObject::kOverwrite);
            if (graphPCMIncGammaSysAErr2050) graphPCMIncGammaSysAErr2050->Write("IncGammaSpec_PCM_SysAErr",TObject::kOverwrite);
            if (graphPCMIncGammaSysBErr2050) graphPCMIncGammaSysBErr2050->Write("IncGammaSpec_PCM_SysBErr",TObject::kOverwrite);
            if (graphPCMIncGammaSysCErr2050) graphPCMIncGammaSysCErr2050->Write("IncGammaSpec_PCM_SysCErr",TObject::kOverwrite);

            if (histoPHOSIncGammaStatErr2050) histoPHOSIncGammaStatErr2050->Write("hIncGammaSpec_PHOS_StatErr",TObject::kOverwrite);
            if (graphPHOSIncGammaSysAErr2050) graphPHOSIncGammaSysAErr2050->Write("IncGammaSpec_PHOS_SysAErr",TObject::kOverwrite);
            if (graphPHOSIncGammaSysBErr2050) graphPHOSIncGammaSysBErr2050->Write("IncGammaSpec_PHOS_SysBErr",TObject::kOverwrite);
            if (graphPHOSIncGammaSysCErr2050) graphPHOSIncGammaSysCErr2050->Write("IncGammaSpec_PHOS_SysCErr",TObject::kOverwrite);

            if (histoPCMDRPi0FitStatErr2050) histoPCMDRPi0FitStatErr2050->Write("hDR_PCM_StatErr",TObject::kOverwrite);
            if (graphPCMDRPi0FitSysAErr2050) graphPCMDRPi0FitSysAErr2050->Write("DR_PCM_SysAErr",TObject::kOverwrite);
            if (graphPCMDRPi0FitSysBErr2050) graphPCMDRPi0FitSysBErr2050->Write("DR_PCM_SysBErr",TObject::kOverwrite);
            if (graphPCMDRPi0FitSysCErr2050) graphPCMDRPi0FitSysCErr2050->Write("DR_PCM_SysCErr",TObject::kOverwrite);

            if (histoPHOSDRPi0FitStatErr2050) histoPHOSDRPi0FitStatErr2050->Write("hDR_PHOS_StatErr",TObject::kOverwrite);
            if (graphPHOSDRPi0FitSysAErr2050) graphPHOSDRPi0FitSysAErr2050->Write("DR_PHOS_SysAErr",TObject::kOverwrite);
            if (graphPHOSDRPi0FitSysBErr2050) graphPHOSDRPi0FitSysBErr2050->Write("DR_PHOS_SysBErr",TObject::kOverwrite);
            if (graphPHOSDRPi0FitSysCErr2050) graphPHOSDRPi0FitSysCErr2050->Write("DR_PHOS_SysCErr",TObject::kOverwrite);

            if (graphCombRAADirGammaStat2050) graphCombRAADirGammaStat2050->Write("RAA_comb_StatErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSys2050) graphCombRAADirGammaSys2050->Write("RAA_comb_SysErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSum2050Ar) graphCombRAADirGammaSum2050Ar->Write("RAA_comb_upperLimits",TObject::kOverwrite);

    fileGammaSpectrum->Write();
    fileGammaSpectrum->Close();

}