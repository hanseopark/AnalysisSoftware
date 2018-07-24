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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/ExtractSignalBinning.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"



void CompareGammaResultsDiffSystems(    TString inputFileNamePP2760GeV      = "",
                                        TString inputFileNamePP8TeV         = "",
                                        TString inputFileNamePPb5TeV        = "",
                                        TString inputFileNamePbPb2760GeV    = "",
                                        TString suffix                      = "eps"
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
    TString outputDir                                           = Form("%s/%s/CompareDiffSystems",suffix.Data(),dateForOutput.Data());

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputGammapp2760GeV.root", inputFileNamePP2760GeV.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputGammapp8TeV.root", inputFileNamePP8TeV.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputGammapPb5TeV.root", inputFileNamePPb5TeV.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputGammaPbPb2760GeV.root", inputFileNamePbPb2760GeV.Data(), outputDir.Data()));

    TString fileNameTheoryPP                                    = "ExternalInput/Theory/TheoryCompilationPP.root";
    TString fileNameTheoryPPb                                   = "ExternalInputPPb/Theory/TheoryCompilationPPb.root";
    TString fileNameTheoryPbPb                                  = "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";


    //*******************************************************************************************************************************************
    //******************************************************* set ranges for plotting ***********************************************************
    //*******************************************************************************************************************************************
    Double_t doubleRatio[2];
    Double_t doubleRatioX[2];
    doubleRatio[0]      = 0.75;     doubleRatio[1]      = 1.65;
    doubleRatioX[0]     = 0.23;     doubleRatioX[1]     = 25.5;

    Color_t colorCocktailPi0                        = kRed+2;
    Color_t colorCocktailEta                        = kBlue+1;
    Color_t colorCocktailEtaP                       = kOrange+1;
    Color_t colorCocktailOmega                      = kYellow+2;
    Color_t colorCocktailPhi                        = kViolet;
    Color_t colorCocktailRho0                       = kAzure-2;
    Color_t colorCocktailSigma0                     = kGray+1;

    Color_t  colorNLOWerner                         = kAzure+2;
    Color_t  colorNLOWernerBand                     = kAzure-9;
    Color_t  colorNLOMcGill                         = kGreen+2;
    Color_t  colorNLOMcGillBand                     = kGreen-6;
    Style_t  styleMarkerNLOWerner                   = 24;
    Style_t  styleLineNLOWerner                     = 5;
    Style_t  styleLineMcGill                        = 7;
    Width_t  widthLineNLO                           = 2.;

    Color_t colorCombpPb                            = GetColorDefaultColor("pPb_5.023TeV", "", "");
    Color_t colorCombpPbBox                         = GetColorDefaultColor("pPb_5.023TeV", "", "", kTRUE);
    Style_t markerStyleCombpPb                      = GetDefaultMarkerStyle("pPb_5.023TeV", "", "");
    Size_t markerSizeCombpPb                        = GetDefaultMarkerSize("pPb_5.023TeV", "", "");

    Color_t colorComb0020                           = GetColorDefaultColor("PbPb_2.76TeV", "", "0-20%");
    Color_t colorComb2040                           = GetColorDefaultColor("PbPb_2.76TeV", "", "20-40%");
    Color_t colorComb4080                           = GetColorDefaultColor("PbPb_2.76TeV", "", "40-80%");
    Color_t colorComb0020Ar                         = kRed-5;
    Color_t colorComb2040Ar                         = kYellow-5;
    Color_t colorComb4080Ar                         = kAzure-7;
    Style_t markerStyleComb0020                     = GetDefaultMarkerStyle("PbPb_2.76TeV", "", "0-20%");
    Style_t markerStyleComb2040                     = GetDefaultMarkerStyle("PbPb_2.76TeV", "", "20-40%");
    Style_t markerStyleComb4080                     = GetDefaultMarkerStyle("PbPb_2.76TeV", "", "40-80%");
    Style_t markerStyleComb0020MC                   = GetDefaultMarkerStyle("PbPb_2.76TeV", "MC", "0-20%");
    Style_t markerStyleComb2040MC                   = GetDefaultMarkerStyle("PbPb_2.76TeV", "MC", "20-40%");
    Style_t markerStyleComb4080MC                   = GetDefaultMarkerStyle("PbPb_2.76TeV", "MC", "40-80%");
    Size_t markerSizeComb0020                       = GetDefaultMarkerSize("PbPb_2.76TeV", "", "0-20%");
    Size_t markerSizeComb2040                       = GetDefaultMarkerSize("PbPb_2.76TeV", "", "20-40%");
    Size_t markerSizeComb4080                       = GetDefaultMarkerSize("PbPb_2.76TeV", "", "40-80%");

    TString collisionSystempp2760GeV                = "pp, #sqrt{#it{s}} = 2.76 TeV";
    TString collisionSystempp8TeV                   = "pp, #sqrt{#it{s}} = 8 TeV";
    TString collisionSystemPbPb760GeV               = "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString collisionSystempPb5TeV                  = "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    TString textALICE                               = "ALICE";

    Color_t colorCombpp2760GeV                      = kBlack;
    Style_t markerStyleCombpp2760GeV                = 20;
    Size_t markerSizeCombpp2760GeV                  = 1.8;
    Color_t colorCombpp8TeV                         = kGray;
    Style_t markerStyleCombpp8TeV                   = 24;
    Size_t markerSizeCombpp8TeV                     = 1.8;
    Width_t widthLinesBoxes                         = 1.4;
    Width_t widthCommonFit                          = 2.4;

    Double_t nColl0020                              = GetNCollFromName("0020");
    Double_t nColl2040                              = GetNCollFromName("2040");
    Double_t nColl4080                              = GetNCollFromName("4080");
    
    Color_t colorEpp2760GeV                         = GetColorDefaultColor("2.76TeV", "", "");
    Style_t markerStyleEpp2760GeV                   = GetDefaultMarkerStyle("2.76TeV", "", "");
    Size_t markerSizeEpp2760GeV                     = GetDefaultMarkerSize("2.76TeV", "", "");
    
    Color_t colorEpp8TeV                            = GetColorDefaultColor("8TeV", "", "");
    Style_t markerStyleEpp8TeV                      = GetDefaultMarkerStyle("8TeV", "", "");
    Size_t markerSizeEpp8TeV                        = GetDefaultMarkerSize("8TeV", "", "");


    //*************************************************************************************************************************************************
    //*************************************** Read in data ********************************************************************************************
    //*************************************************************************************************************************************************
    //--------------------------------------- pp 2.76TeV --------------------------------------------
    TFile* fileCombPP2760GeV                                = new TFile( inputFileNamePP2760GeV.Data());
    TGraphAsymmErrors* graphDRStatpp2760GeV                 = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphRGammaNonFitCombStatErr");
    TGraphAsymmErrors* graphDRSyspp2760GeV                  = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphRGammaNonFitCombSysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaStatpp2760GeV   = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvYieldDirGammaNonFitStatErr");
    TGraphAsymmErrors* graphInvYieldDirGammaSyspp2760GeV    = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvYieldDirGammaNonFitSysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaTotArpp2760GeV  = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvYieldDirGammaNonFitSumErrAr");
    TGraphAsymmErrors* graphInvYieldIncGammaStatpp2760GeV   = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvYieldIncGammaStatErr");
    TGraphAsymmErrors* graphInvYieldIncGammaSyspp2760GeV    = (TGraphAsymmErrors*) fileCombPP2760GeV->Get("Gamma2.76TeV/graphInvYieldIncGammaSysErr");
    TF1* fitInvYieldIncGammapp2760GeV                       = (TF1*) fileCombPP2760GeV->Get("Gamma2.76TeV/Fits/TwoComponentModelFitGamma");

    //--------------------------------------- pp 8TeV --------------------------------------------
    TFile* fileCombPP8TeV                                   = new TFile( inputFileNamePP8TeV.Data());
    TGraphAsymmErrors* graphDRStatpp8TeV                    = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphRGammaCombNonFitStatErr");
    TGraphAsymmErrors* graphDRSyspp8TeV                     = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphRGammaCombNonFitSysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaStatpp8TeV      = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvYieldDirGammaNonFitStatErr");
    TGraphAsymmErrors* graphInvYieldDirGammaSyspp8TeV       = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvYieldDirGammaNonFitSysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaTotArpp8TeV     = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvYieldDirGammaNonFitSumErrAr");
    TGraphAsymmErrors* graphInvYieldIncGammaStatpp8TeV      = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvYieldIncGammaStatErr");
    TGraphAsymmErrors* graphInvYieldIncGammaSyspp8TeV       = (TGraphAsymmErrors*) fileCombPP8TeV->Get("Gamma8TeV/graphInvYieldIncGammaSysErr");
    TF1* fitInvYieldIncGammapp8TeV                          = (TF1*) fileCombPP8TeV->Get("Gamma8TeV/Fits/TwoComponentModelFitGamma");
    //--------------------------------------- pPb 5TeV --------------------------------------------
    TFile* fileCombPPb5TeV                                  = new TFile( inputFileNamePPb5TeV.Data());
    TGraphAsymmErrors* graphDRStatpPb5TeV                   = (TGraphAsymmErrors*) fileCombPPb5TeV->Get("Gamma_pPb5TeV/graphRGammaCombStatErr");
    TGraphAsymmErrors* graphDRSyspPb5TeV                    = (TGraphAsymmErrors*) fileCombPPb5TeV->Get("Gamma_pPb5TeV/graphRGammaCombSysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaStatpPb5TeV     = (TGraphAsymmErrors*) fileCombPPb5TeV->Get("Gamma_pPb5TeV/graphInvYieldDirGammaStatErr");
    TGraphAsymmErrors* graphInvYieldDirGammaSyspPb5TeV      = (TGraphAsymmErrors*) fileCombPPb5TeV->Get("Gamma_pPb5TeV/graphInvYieldDirGammaSysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaTotArpPb5TeV    = (TGraphAsymmErrors*) fileCombPPb5TeV->Get("Gamma_pPb5TeV/graphInvYieldDirGammaSumErrAr");
    //--------------------------------------- PbPb 2.76TeV --------------------------------------------
    TFile* fileCombPbPb2760GeV                              = new TFile( inputFileNamePbPb2760GeV.Data());
    TDirectory* dirGammaPbPb0020                            = (TDirectory*)fileCombPbPb2760GeV->Get("Gamma_PbPb_2.76TeV_0-20%");
    TGraphAsymmErrors* graphDRStatPbPb0020                  = (TGraphAsymmErrors*) dirGammaPbPb0020->Get("DR_comb_StatErr");
    TGraphAsymmErrors* graphDRSysPbPb0020                   = (TGraphAsymmErrors*) dirGammaPbPb0020->Get("DR_comb_SysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaStatPbPb0020    = (TGraphAsymmErrors*) dirGammaPbPb0020->Get("DirGammaSpec_comb_StatErr");
    TGraphAsymmErrors* graphInvYieldDirGammaSysPbPb0020     = (TGraphAsymmErrors*) dirGammaPbPb0020->Get("DirGammaSpec_comb_SysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaTotArPbPb0020   = (TGraphAsymmErrors*) dirGammaPbPb0020->Get("DirGammaSpec_comb_upperLimits");
    TDirectory* dirGammaPbPb2040                            = (TDirectory*)fileCombPbPb2760GeV->Get("Gamma_PbPb_2.76TeV_20-40%");
    TGraphAsymmErrors* graphDRStatPbPb2040                  = (TGraphAsymmErrors*) dirGammaPbPb2040->Get("DR_comb_StatErr");
    TGraphAsymmErrors* graphDRSysPbPb2040                   = (TGraphAsymmErrors*) dirGammaPbPb2040->Get("DR_comb_SysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaStatPbPb2040    = (TGraphAsymmErrors*) dirGammaPbPb2040->Get("DirGammaSpec_comb_StatErr");
    TGraphAsymmErrors* graphInvYieldDirGammaSysPbPb2040     = (TGraphAsymmErrors*) dirGammaPbPb2040->Get("DirGammaSpec_comb_SysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaTotArPbPb2040   = (TGraphAsymmErrors*) dirGammaPbPb2040->Get("DirGammaSpec_comb_upperLimits");
    TDirectory* dirGammaPbPb4080                            = (TDirectory*)fileCombPbPb2760GeV->Get("Gamma_PbPb_2.76TeV_40-80%");
    TGraphAsymmErrors* graphDRStatPbPb4080                  = (TGraphAsymmErrors*) dirGammaPbPb4080->Get("DR_comb_StatErr");
    TGraphAsymmErrors* graphDRSysPbPb4080                   = (TGraphAsymmErrors*) dirGammaPbPb4080->Get("DR_comb_SysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaStatPbPb4080    = (TGraphAsymmErrors*) dirGammaPbPb4080->Get("DirGammaSpec_comb_StatErr");
    TGraphAsymmErrors* graphInvYieldDirGammaSysPbPb4080     = (TGraphAsymmErrors*) dirGammaPbPb4080->Get("DirGammaSpec_comb_SysErr");
    TGraphAsymmErrors* graphInvYieldDirGammaTotArPbPb4080   = (TGraphAsymmErrors*) dirGammaPbPb4080->Get("DirGammaSpec_comb_upperLimits");


    TFile* fileTheoryPP                                     = new TFile( fileNameTheoryPP.Data());
    TGraph* graphTheoryNLODRpp2760GeVPaquettCenter          = (TGraph*) fileTheoryPP->Get("DirectPhoton/graphDRNLOPaquett_2760GeV_ALICECocktail");

//     TFile* fileTheoryPPb                                    = new TFile( fileNameTheoryPPb.Data());
//     TGraph* graphTheoryNLODRpPbCenter                       = (TGraph*) fileTheoryPPb->Get("pPb_5.023TeV/graphRGammaDirectPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10_ALICECocktail_Center");
//     while (graphTheoryNLODRpPbCenter->GetX()[0] < 3)    graphTheoryNLODRpPbCenter->RemovePoint(0);
//     TGraph* graphTheoryMCGillDRpPbCenter                    = (TGraph*) fileTheoryPPb->Get("pPb_5.023TeV/graphRGammaDirectPhotonSpecMcGill5023GeV_ALICECocktail_Center");

//     TFile* fileTheoryPbPb                                   = new TFile( fileNameTheorypp2760GeV.Data());

    //*************************************************************************************************************************************************
    //*************************************** Prepare for plotting ************************************************************************************
    //*************************************************************************************************************************************************

    TGraphAsymmErrors* graphInvYieldDirGammaStatPbPb0020Plot;
    TGraphAsymmErrors* graphInvYieldDirGammaSysPbPb0020Plot;
    TGraphAsymmErrors* graphInvYieldDirGammaTotArPbPb0020Plot;
    if (graphInvYieldDirGammaStatPbPb0020) graphInvYieldDirGammaStatPbPb0020Plot        = ScaleGraph(graphInvYieldDirGammaStatPbPb0020,100);
    if (graphInvYieldDirGammaStatPbPb0020Plot) ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPb0020Plot);
    if (graphInvYieldDirGammaSysPbPb0020) graphInvYieldDirGammaSysPbPb0020Plot          = ScaleGraph(graphInvYieldDirGammaSysPbPb0020,100);
    if (graphInvYieldDirGammaTotArPbPb0020) graphInvYieldDirGammaTotArPbPb0020Plot      = ScaleGraph(graphInvYieldDirGammaTotArPbPb0020,100);

    TGraphAsymmErrors* graphInvYieldDirGammaStatPbPb2040Plot;
    TGraphAsymmErrors* graphInvYieldDirGammaSysPbPb2040Plot;
    TGraphAsymmErrors* graphInvYieldDirGammaTotArPbPb2040Plot;
    if (graphInvYieldDirGammaStatPbPb2040) graphInvYieldDirGammaStatPbPb2040Plot        = ScaleGraph(graphInvYieldDirGammaStatPbPb2040,1);
    if (graphInvYieldDirGammaStatPbPb2040Plot) ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPb2040Plot);
    if (graphInvYieldDirGammaSysPbPb2040) graphInvYieldDirGammaSysPbPb2040Plot          = ScaleGraph(graphInvYieldDirGammaSysPbPb2040,1);
    if (graphInvYieldDirGammaTotArPbPb2040) graphInvYieldDirGammaTotArPbPb2040Plot      = ScaleGraph(graphInvYieldDirGammaTotArPbPb2040,1);

    TGraphAsymmErrors* graphInvYieldDirGammaStatPbPb4080Plot;
    TGraphAsymmErrors* graphInvYieldDirGammaSysPbPb4080Plot;
    TGraphAsymmErrors* graphInvYieldDirGammaTotArPbPb4080Plot;
    if (graphInvYieldDirGammaStatPbPb4080) graphInvYieldDirGammaStatPbPb4080Plot        = ScaleGraph(graphInvYieldDirGammaStatPbPb4080,1e-2);
    if (graphInvYieldDirGammaStatPbPb4080Plot) ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPb4080Plot);
    if (graphInvYieldDirGammaSysPbPb4080) graphInvYieldDirGammaSysPbPb4080Plot          = ScaleGraph(graphInvYieldDirGammaSysPbPb4080,1e-2);
    if (graphInvYieldDirGammaTotArPbPb4080) graphInvYieldDirGammaTotArPbPb4080Plot      = ScaleGraph(graphInvYieldDirGammaTotArPbPb4080,1e-2);


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

    TH1D* dummyDirGamma = new TH1D("dummyDirGamma", "dummyDirGamma", 1000, 0., 22.);
    SetStyleHistoTH1ForGraphs( dummyDirGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{inel}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.8);
    dummyDirGamma->GetYaxis()->SetRangeUser( 1.2e-9,1.5e1);
    dummyDirGamma->GetXaxis()->SetLabelOffset(-0.015);
    dummyDirGamma->GetXaxis()->SetTickLength(0.025);
    dummyDirGamma->GetYaxis()->SetTickLength(0.025);
    dummyDirGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    //     dummyDirGamma->GetXaxis()->SetRangeUser(0,16);
    dummyDirGamma->DrawCopy();

//         TLatex *labelScalingDirGamma0020 = new TLatex(11.2,1.3E-3,"x 10^{2}");
//         SetStyleTLatex( labelScalingDirGamma0020, 0.85*textsizeLabelsDirGamma,4,colorComb0020,42,kFALSE);
//         labelScalingDirGamma0020->Draw();
//
//         TLatex *labelScalingDirGamma2040 = new TLatex(11.2,6.0E-5,"x 10^{1}");
//         SetStyleTLatex( labelScalingDirGamma2040, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
//         labelScalingDirGamma2040->Draw();
//
//         TLatex *labelScalingDirGamma4080 = new TLatex(11.2,7.5E-7,"x 10^{0}");
//         SetStyleTLatex( labelScalingDirGamma4080, 0.85*textsizeLabelsDirGamma,4,colorComb4080,42,kFALSE);
//         labelScalingDirGamma4080->Draw();

        if (graphInvYieldDirGammaSysPbPb0020){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSysPbPb0020, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSysPbPb0020->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatPbPb0020){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPb0020);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatPbPb0020, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);
            graphInvYieldDirGammaStatPbPb0020->Draw("p,E1Z,same");
        }

        if (graphInvYieldDirGammaSysPbPb2040){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSysPbPb2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSysPbPb2040->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatPbPb2040){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPb2040);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatPbPb2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphInvYieldDirGammaStatPbPb2040->Draw("p,E1Z,same");
        }
        if (graphInvYieldDirGammaTotArPbPb2040){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArPbPb2040 , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
            graphInvYieldDirGammaTotArPbPb2040->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPb2040);
        }
        if (graphInvYieldDirGammaSysPbPb4080){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSysPbPb4080, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSysPbPb4080->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatPbPb4080){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPb4080);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatPbPb4080, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);
            graphInvYieldDirGammaStatPbPb4080->Draw("p,E1Z,same");
        }
        if (graphInvYieldDirGammaTotArPbPb4080){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArPbPb4080 , 1, 3, colorComb4080, colorComb4080, 1.8, kTRUE);
            graphInvYieldDirGammaTotArPbPb4080->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPb4080);
        }

        if (graphInvYieldDirGammaSyspp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspp2760GeV , markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV, colorCombpp2760GeV, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspp2760GeV->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatpp2760GeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatpp2760GeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV, colorCombpp2760GeV);
            graphInvYieldDirGammaStatpp2760GeV->Draw("p,E1Z,same");
        }

        if (graphInvYieldDirGammaTotArpp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpp2760GeV , 1, 3, colorCombpp2760GeV, colorCombpp2760GeV, 1.8, kTRUE);
            graphInvYieldDirGammaTotArpp2760GeV->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeV);
        }
        if (graphInvYieldDirGammaSyspPb5TeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspPb5TeV, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspPb5TeV->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatpPb5TeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatpPb5TeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpPb5TeV, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb);
            graphInvYieldDirGammaStatpPb5TeV->Draw("p,E1Z,same");
        }
        if (graphInvYieldDirGammaTotArpPb5TeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpPb5TeV , 1, 3, colorCombpPb, colorCombpPb, 1.8, kTRUE);
            graphInvYieldDirGammaTotArpPb5TeV->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpPb5TeV);
        }

        //         TLatex *labelScalingDirGamma4080 = new TLatex(11.2,7.5E-7,"x 10^{0}");
        //         SetStyleTLatex( labelScalingDirGamma4080, 0.85*textsizeLabelsDirGamma,4,colorComb4080,42,kFALSE);
        //         labelScalingDirGamma4080->Draw();


        TLatex *labelALICEDirGamma = new TLatex(0.95,0.94,textALICE);
        SetStyleTLatex( labelALICEDirGamma, 42, 4, 1, 43, kTRUE, 31);

//         labelALICEDirGamma->Draw();

        TLegend* legendDirGammaPP = GetAndSetLegend2(0.21, 0.10+(6*textsizeLabelsDirGamma*0.85), 0.21+0.21, 0.10+(8*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1,
                                                     collisionSystempp2760GeV.Data(), 42, 0.25);
//         legendDirGammaPP->AddEntry((TObject*)0,"ALICE","");
        legendDirGammaPP->AddEntry(graphInvYieldDirGammaSyspp2760GeV,"ALICE","pf");

        legendDirGammaPP->Draw();
        TLegend* legendDirGammaPPb = GetAndSetLegend2(0.21, 0.10+(4*textsizeLabelsDirGamma*0.85), 0.21+0.21, 0.10+(6*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1,
                                                      collisionSystempPb5TeV.Data(), 42, 0.25);
        legendDirGammaPPb->AddEntry(graphInvYieldDirGammaSyspPb5TeV,"ALICE","pf");
        legendDirGammaPPb->Draw();

        TLegend* legendDirGamma = GetAndSetLegend2(0.21, 0.10, 0.21+0.21, 0.10+(4*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, collisionSystemPbPb760GeV.Data(), 42, 0.25);
        legendDirGamma->AddEntry(graphInvYieldDirGammaSysPbPb0020,"  0-20% ALICE","pf");
        legendDirGamma->AddEntry(graphInvYieldDirGammaSysPbPb2040,"20-40% ALICE","pf");
        legendDirGamma->AddEntry(graphInvYieldDirGammaSysPbPb4080,"40-80% ALICE","pf");
        legendDirGamma->Draw();

//         TGraphAsymmErrors* dummyForLegend    = new TGraphAsymmErrors(1);
//         dummyForLegend->SetPoint(0,0.315,1e-6);
//         dummyForLegend->SetPointError(0,0,0,4e-7,0);
//         DrawGammaSetMarkerTGraphAsym(dummyForLegend , 1, 3, colorCombpp2760GeV, colorCombpp2760GeV, 1.8, kTRUE);
//         dummyForLegend->Draw(">,same");
//         PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend,0.01);


    canvasDirGamma->Print(Form("%s/DirGammaSpectra_Unscaled.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectra_Unscaled.pdf",outputDir.Data()));
    
    dummyDirGamma->GetYaxis()->SetRangeUser( 1.2e-10,8.5);
    dummyDirGamma->DrawCopy();

        if (graphInvYieldDirGammaSyspp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspp2760GeV , markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspp2760GeV->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatpp2760GeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatpp2760GeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpp2760GeV, markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV);
            graphInvYieldDirGammaStatpp2760GeV->Draw("p,E1Z,same");
        }

        if (graphInvYieldDirGammaTotArpp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpp2760GeV , 1, 3, colorEpp2760GeV, colorEpp2760GeV, 1.8, kTRUE);
            graphInvYieldDirGammaTotArpp2760GeV->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeV);
        }
        if (graphInvYieldDirGammaSyspp8TeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspp8TeV , markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspp8TeV->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatpp8TeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatpp8TeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpp8TeV, markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV);
            graphInvYieldDirGammaStatpp8TeV->Draw("p,E1Z,same");
        }

        if (graphInvYieldDirGammaTotArpp8TeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpp8TeV , 1, 3, colorEpp8TeV, colorEpp8TeV, 1.8, kTRUE);
            graphInvYieldDirGammaTotArpp8TeV->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp8TeV);
        }
      
        // TLatex *labelEnergyInvYieldPaperAll = new TLatex(0.22, 0.13+0.04*2, collisionSystempp2760GeV.Data());
        // SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        // labelEnergyInvYieldPaperAll->Draw();
        // TLatex *labelALICENormUnPaperAll    = new TLatex(0.22,0.135+0.07*0,"norm. unc. 2.5%");
        // SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixelDirGam*0.85,4, 1, 43, kTRUE, 11);
        // labelALICENormUnPaperAll->Draw();
        // 
        // TLegend* legendDirGammaPPonly = GetAndSetLegend2(0.21, 0.13+(1*textsizeLabelsDirGamma*0.85), 0.21+0.21, 0.13+(2*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, "", 42, 0.25);
        // legendDirGammaPPonly->AddEntry(graphInvYieldDirGammaSyspp2760GeV,"#gamma_{dir} ALICE","pf");
        // legendDirGammaPPonly->Draw();
        // 
        // TLatex *labelEnergyInvYieldPaperAll2 = new TLatex(0.94, 0.83+0.04*2, collisionSystempp8TeV.Data());
        // SetStyleTLatex( labelEnergyInvYieldPaperAll2, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 31);
        // labelEnergyInvYieldPaperAll2->Draw();
        // TLatex *labelALICENormUnPaperAll2    = new TLatex(0.94,0.835+0.07*0,"norm. unc. 2.6%");
        // SetStyleTLatex( labelALICENormUnPaperAll2, textSizeLabelsPixelDirGam*0.85,4, 1, 43, kTRUE, 31);
        // labelALICENormUnPaperAll2->Draw();
        // 
        // TLegend* legendDirGammaPPonly2 = GetAndSetLegend2(0.71, 0.83+(1*textsizeLabelsDirGamma*0.85), 0.75+0.21, 0.83+(2*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, "", 42, 0.25);
        // legendDirGammaPPonly2->AddEntry(graphInvYieldDirGammaSyspp8TeV,"#gamma_{dir} ALICE","pf");
        // legendDirGammaPPonly2->Draw();
        
        TLatex *labelEnergyInvYieldPaperAll = new TLatex(0.23, 0.11+0.04*4, "#gamma_{dir} ALICE");
        SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        labelEnergyInvYieldPaperAll->Draw();
        TLegend* legendDirGammaPPonly2 = GetAndSetLegend2(0.22, 0.11+(0*textsizeLabelsDirGamma*0.85), 0.21+0.25, 0.11+(4*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, "", 42, 0.25);
        legendDirGammaPPonly2->AddEntry(graphInvYieldDirGammaSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        legendDirGammaPPonly2->AddEntry((TObject*)0,"norm. unc. 2.5%","");
        legendDirGammaPPonly2->AddEntry(graphInvYieldDirGammaSyspp8TeV,collisionSystempp8TeV.Data(),"pf");
        legendDirGammaPPonly2->AddEntry((TObject*)0,"norm. unc. 2.6%","");
        legendDirGammaPPonly2->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectraPP.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectraPP.pdf",outputDir.Data()));
    
    TH1D* dummyIncGamma = new TH1D("dummyIncGamma", "dummyIncGamma", 1000, 0., 22.);
    SetStyleHistoTH1ForGraphs( dummyIncGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{inel}} #frac{d^{2}#it{N}_{#gamma_{inc}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.8);
    dummyIncGamma->GetYaxis()->SetRangeUser( 1.2e-9,1.5e1);
    dummyIncGamma->GetXaxis()->SetLabelOffset(-0.015);
    dummyIncGamma->GetXaxis()->SetTickLength(0.025);
    dummyIncGamma->GetYaxis()->SetTickLength(0.025);
    dummyIncGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyIncGamma->DrawCopy();
        if (graphInvYieldIncGammaSyspp2760GeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldIncGammaSyspp2760GeV , markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV, widthLinesBoxes, kTRUE);
            graphInvYieldIncGammaSyspp2760GeV->Draw("E2same");
        }
        if (graphInvYieldIncGammaStatpp2760GeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldIncGammaStatpp2760GeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldIncGammaStatpp2760GeV, markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV, colorEpp2760GeV);
            graphInvYieldIncGammaStatpp2760GeV->Draw("p,E1Z,same");
        }

        if (graphInvYieldIncGammaSyspp8TeV){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldIncGammaSyspp8TeV , markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV, widthLinesBoxes, kTRUE);
            graphInvYieldIncGammaSyspp8TeV->Draw("E2same");
        }
        if (graphInvYieldIncGammaStatpp8TeV){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldIncGammaStatpp8TeV);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldIncGammaStatpp8TeV, markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV, colorEpp8TeV);
            graphInvYieldIncGammaStatpp8TeV->Draw("p,E1Z,same");
        }

        TLatex *labelEnergyIncInvYieldPaperAll = new TLatex(0.23, 0.11+0.04*4, "#gamma_{inc} ALICE");
        SetStyleTLatex( labelEnergyIncInvYieldPaperAll, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        labelEnergyIncInvYieldPaperAll->Draw();
        legendDirGammaPPonly2->Draw();
        
        
    canvasDirGamma->Print(Form("%s/IncGammaSpectraPP.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/IncGammaSpectraPP.pdf",outputDir.Data()));
    
    dummyIncGamma->DrawCopy();
    
        if (graphInvYieldIncGammaSyspp2760GeV){
            graphInvYieldIncGammaSyspp2760GeV->Draw("E2same");
        }
        if (graphInvYieldIncGammaStatpp2760GeV){
            graphInvYieldIncGammaStatpp2760GeV->Draw("p,E1Z,same");
        }

        if (graphInvYieldIncGammaSyspp8TeV){
            graphInvYieldIncGammaSyspp8TeV->Draw("E2same");
        }
        if (graphInvYieldIncGammaStatpp8TeV){
            graphInvYieldIncGammaStatpp8TeV->Draw("p,E1Z,same");
        }
        if(fitInvYieldIncGammapp2760GeV){
          DrawGammaSetMarkerTF1( fitInvYieldIncGammapp2760GeV, 3, 2, colorEpp2760GeV);
          fitInvYieldIncGammapp2760GeV->SetRange(0.45,9.7);
          fitInvYieldIncGammapp2760GeV->Draw("same");
        }
        if(fitInvYieldIncGammapp8TeV){
          DrawGammaSetMarkerTF1( fitInvYieldIncGammapp8TeV, 3, 2, colorEpp8TeV);
          fitInvYieldIncGammapp8TeV->SetRange(0.3,18.);
          fitInvYieldIncGammapp8TeV->Draw("same");
        }
        TF1* fitInvYieldIncGammaDummy = (TF1*)fitInvYieldIncGammapp2760GeV->Clone("fitInvYieldIncGammaDummy");
        DrawGammaSetMarkerTF1( fitInvYieldIncGammaDummy, 3, 2, kGray+2);
        TLatex *labelEnergyIncInvYieldPaperAllwFit = new TLatex(0.23, 0.11+0.04*5, "#gamma_{inc} ALICE");
        SetStyleTLatex( labelEnergyIncInvYieldPaperAllwFit, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        labelEnergyIncInvYieldPaperAllwFit->Draw();
        TLegend* legendDirGammaPPonlywFits = GetAndSetLegend2(0.22, 0.11+(0*textsizeLabelsDirGamma*0.85), 0.21+0.25, 0.11+(5*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, "", 42, 0.25);
        legendDirGammaPPonlywFits->AddEntry(graphInvYieldDirGammaSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        legendDirGammaPPonlywFits->AddEntry((TObject*)0,"norm. unc. 2.5%","");
        legendDirGammaPPonlywFits->AddEntry(graphInvYieldDirGammaSyspp8TeV,collisionSystempp8TeV.Data(),"pf");
        legendDirGammaPPonlywFits->AddEntry((TObject*)0,"norm. unc. 2.6%","");
        legendDirGammaPPonlywFits->AddEntry(fitInvYieldIncGammaDummy,"TCM fit","l");
        legendDirGammaPPonlywFits->Draw();
        
        
    canvasDirGamma->Print(Form("%s/IncGammaSpectraPP_wFit.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/IncGammaSpectraPP_wFit.pdf",outputDir.Data()));

    dummyDirGamma->GetYaxis()->SetRangeUser( 1.2e-10,8.5e4);
    dummyDirGamma->DrawCopy();

        TLatex *labelScalingDirGamma0020 = new TLatex(12.2,1.2E-3,"x 10^{2}");
        SetStyleTLatex( labelScalingDirGamma0020, 0.85*textsizeLabelsDirGamma,4,colorComb0020,42,kFALSE);
        labelScalingDirGamma0020->Draw();
        TLatex *labelScalingDirGamma2040 = new TLatex(12.2,6.0E-6,"x 10^{0}");
        SetStyleTLatex( labelScalingDirGamma2040, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
        labelScalingDirGamma2040->Draw();
        TLatex *labelScalingDirGamma4080 = new TLatex(12.2,7.5E-9,"x 10^{-2}");
        SetStyleTLatex( labelScalingDirGamma4080, 0.85*textsizeLabelsDirGamma,4,colorComb4080,42,kFALSE);
        labelScalingDirGamma4080->Draw();

        if (graphInvYieldDirGammaSysPbPb0020Plot){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSysPbPb0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSysPbPb0020Plot->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatPbPb0020Plot){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPb0020Plot);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatPbPb0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);
            graphInvYieldDirGammaStatPbPb0020Plot->Draw("p,E1Z,same");
        }

        if (graphInvYieldDirGammaSysPbPb2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSysPbPb2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSysPbPb2040Plot->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatPbPb2040Plot){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPb2040Plot);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatPbPb2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphInvYieldDirGammaStatPbPb2040Plot->Draw("p,E1Z,same");
        }
        if (graphInvYieldDirGammaTotArPbPb2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArPbPb2040Plot , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
            graphInvYieldDirGammaTotArPbPb2040Plot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPb2040Plot);
        }
        if (graphInvYieldDirGammaSysPbPb4080Plot){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSysPbPb4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSysPbPb4080Plot->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatPbPb4080Plot){
            ProduceGraphAsymmWithoutXErrors(graphInvYieldDirGammaStatPbPb4080Plot);
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatPbPb4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);
            graphInvYieldDirGammaStatPbPb4080Plot->Draw("p,E1Z,same");
        }
        if (graphInvYieldDirGammaTotArPbPb4080Plot){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArPbPb4080Plot , 1, 3, colorComb4080, colorComb4080, 1.8, kTRUE);
            graphInvYieldDirGammaTotArPbPb4080Plot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPb4080Plot);
        }
        TGraphAsymmErrors* graphInvYieldDirGammaTotArpp2760GeVScaled0020    = ScaleGraph(graphInvYieldDirGammaTotArpp2760GeV,nColl0020*100);
        if (graphInvYieldDirGammaTotArpp2760GeVScaled0020){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpp2760GeVScaled0020 , 1, 3, colorComb0020Ar, colorComb0020Ar, 2.0, kTRUE);
            graphInvYieldDirGammaTotArpp2760GeVScaled0020->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeVScaled0020);
        }
        TGraphAsymmErrors* graphInvYieldDirGammaTotArpp2760GeVScaled2040    = ScaleGraph(graphInvYieldDirGammaTotArpp2760GeV,nColl2040*1);
        if (graphInvYieldDirGammaTotArpp2760GeVScaled2040){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpp2760GeVScaled2040 , 1, 3, colorComb2040Ar, colorComb2040Ar, 2.0, kTRUE);
            graphInvYieldDirGammaTotArpp2760GeVScaled2040->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeVScaled2040);
        }

        TGraphAsymmErrors* graphInvYieldDirGammaTotArpp2760GeVScaled4080    = ScaleGraph(graphInvYieldDirGammaTotArpp2760GeV,nColl4080*1e-2);
        TGraph* graphInvYieldDirGammaTotArpp2760GeVScaled4080Excl           = new TGraph(graphInvYieldDirGammaTotArpp2760GeVScaled4080->GetN());
        for (Int_t i = 0; i < graphInvYieldDirGammaTotArpp2760GeVScaled4080->GetN(); i++){
            graphInvYieldDirGammaTotArpp2760GeVScaled4080Excl->SetPoint(i,graphInvYieldDirGammaTotArpp2760GeVScaled4080->GetX()[i],graphInvYieldDirGammaTotArpp2760GeVScaled4080->GetY()[i]);
        }
        if (graphInvYieldDirGammaTotArpp2760GeVScaled4080){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaTotArpp2760GeVScaled4080 , 1, 3, colorComb4080Ar, colorComb4080Ar, 2.0, kTRUE);
            graphInvYieldDirGammaTotArpp2760GeVScaled4080->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeVScaled4080);
        }


        TGraphAsymmErrors* graphInvYieldDirGammaSyspp2760GeVScaled0020    = ScaleGraph(graphInvYieldDirGammaSyspp2760GeV,nColl0020*100);
        if (graphInvYieldDirGammaSyspp2760GeVScaled0020){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspp2760GeVScaled0020 , markerStyleCombpp2760GeV+4, markerSizeCombpp2760GeV, colorComb0020Ar, colorComb0020Ar, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspp2760GeVScaled0020->Draw("E2,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaSyspp2760GeVScaled0020);
        }
        TGraphAsymmErrors* graphInvYieldDirGammaSyspp2760GeVScaled2040    = ScaleGraph(graphInvYieldDirGammaSyspp2760GeV,nColl2040*1);
        if (graphInvYieldDirGammaSyspp2760GeVScaled2040){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspp2760GeVScaled2040 ,markerStyleCombpp2760GeV+4, markerSizeCombpp2760GeV, colorComb2040Ar, colorComb2040Ar, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspp2760GeVScaled2040->Draw("E2,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaSyspp2760GeVScaled2040);
        }
        TGraphAsymmErrors* graphInvYieldDirGammaSyspp2760GeVScaled4080    = ScaleGraph(graphInvYieldDirGammaSyspp2760GeV,nColl4080*1e-2);
        if (graphInvYieldDirGammaSyspp2760GeVScaled4080){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaSyspp2760GeVScaled4080 , markerStyleCombpp2760GeV+4, markerSizeCombpp2760GeV, colorComb4080Ar, colorComb4080Ar, widthLinesBoxes, kTRUE);
            graphInvYieldDirGammaSyspp2760GeVScaled4080->Draw("E2,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaSyspp2760GeVScaled4080);
        }

        TGraphAsymmErrors* graphInvYieldDirGammaStatpp2760GeVScaled0020    = ScaleGraph(graphInvYieldDirGammaStatpp2760GeV,nColl0020*100);
        if (graphInvYieldDirGammaStatpp2760GeVScaled0020){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpp2760GeVScaled0020 , markerStyleCombpp2760GeV+4, markerSizeCombpp2760GeV, colorComb0020Ar, colorComb0020Ar);
            graphInvYieldDirGammaStatpp2760GeVScaled0020->Draw("p,E1Z,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaStatpp2760GeVScaled0020);
        }
        TGraphAsymmErrors* graphInvYieldDirGammaStatpp2760GeVScaled2040    = ScaleGraph(graphInvYieldDirGammaStatpp2760GeV,nColl2040*1);
        if (graphInvYieldDirGammaStatpp2760GeVScaled2040){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpp2760GeVScaled2040 ,markerStyleCombpp2760GeV+4, markerSizeCombpp2760GeV, colorComb2040Ar, colorComb2040Ar);
            graphInvYieldDirGammaStatpp2760GeVScaled2040->Draw("p,E1Z,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaStatpp2760GeVScaled2040);
        }
        TGraphAsymmErrors* graphInvYieldDirGammaStatpp2760GeVScaled4080    = ScaleGraph(graphInvYieldDirGammaStatpp2760GeV,nColl4080*1e-2);
        if (graphInvYieldDirGammaStatpp2760GeVScaled4080){
            DrawGammaSetMarkerTGraphAsym(graphInvYieldDirGammaStatpp2760GeVScaled4080 , markerStyleCombpp2760GeV+4, markerSizeCombpp2760GeV, colorComb4080Ar, colorComb4080Ar);
            graphInvYieldDirGammaStatpp2760GeVScaled4080->Draw("p,E1Z,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaStatpp2760GeVScaled4080);
        }

//         labelALICEDirGamma->Draw();

        TLegend* legendDirGammaPP2 = GetAndSetLegend2(0.21, 0.10+(4*textsizeLabelsDirGamma*0.85), 0.21+0.21, 0.10+(8*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1,
                                                     collisionSystempp2760GeV.Data(), 42, 0.25);
        legendDirGammaPP2->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled0020,Form("ALICE x %1.1f",nColl0020),"pf");
        legendDirGammaPP2->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled2040,Form("ALICE x %1.1f",nColl2040),"pf");
        legendDirGammaPP2->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled4080,Form("ALICE x %1.1f",nColl4080),"pf");
//         legendDirGammaPP2->AddEntry((TObject*)0,Form("ALICE x %1.1f",nColl0020),"");
//         legendDirGammaPP2->AddEntry((TObject*)0,Form("ALICE x %1.1f",nColl2040),"");
//         legendDirGammaPP2->AddEntry((TObject*)0,Form("ALICE x %1.1f",nColl4080),"");
        legendDirGammaPP2->Draw();
        legendDirGamma->Draw();

//         TGraphAsymmErrors* dummyForLegend2   = new TGraphAsymmErrors(1);
//         dummyForLegend2->SetPoint(0,0.312,2.2e-6);
//         dummyForLegend2->SetPointError(0,0,0,1.2e-6,0);
//         DrawGammaSetMarkerTGraphAsym(dummyForLegend2 , 1, 3, colorComb0020Ar, colorComb0020Ar, 2.0, kTRUE);
//         dummyForLegend2->Draw(">,same");
//         PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend2,0.01);
//
//         TGraphAsymmErrors* dummyForLegend3    = new TGraphAsymmErrors(1);
//         dummyForLegend3->SetPoint(0,0.312,6.3e-7);
//         dummyForLegend3->SetPointError(0,0,0,3.3e-7,0);
//         DrawGammaSetMarkerTGraphAsym(dummyForLegend3 , 1, 3, colorComb2040Ar, colorComb2040Ar, 2.0, kTRUE);
//         dummyForLegend3->Draw(">,same");
//         PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend3,0.01);
//
//         TGraphAsymmErrors* dummyForLegend4    = new TGraphAsymmErrors(1);
//         dummyForLegend4->SetPoint(0,0.312,1.8e-7);
//         dummyForLegend4->SetPointError(0,0,0,0.92e-7,0);
//         DrawGammaSetMarkerTGraphAsym(dummyForLegend4 , 1, 3, colorComb4080Ar, colorComb4080Ar, 2.0, kTRUE);
//         dummyForLegend4->Draw(">,same");
//         PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend4,0.01);

    canvasDirGamma->Print(Form("%s/DirGammaSpectra_WithScaledPP.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectra_WithScaledPP.pdf",outputDir.Data()));

    dummyDirGamma->GetYaxis()->SetRangeUser( 1.2e-10,8.5e5);
    dummyDirGamma->DrawCopy();

        labelScalingDirGamma0020->Draw();
        labelScalingDirGamma2040->Draw();
        labelScalingDirGamma4080->Draw();

        if (graphInvYieldDirGammaSysPbPb0020Plot){
            graphInvYieldDirGammaSysPbPb0020Plot->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatPbPb0020Plot){
            graphInvYieldDirGammaStatPbPb0020Plot->Draw("p,E1Z,same");
        }

        if (graphInvYieldDirGammaSysPbPb2040Plot){
            graphInvYieldDirGammaSysPbPb2040Plot->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatPbPb2040Plot){
            graphInvYieldDirGammaStatPbPb2040Plot->Draw("p,E1Z,same");
        }
        if (graphInvYieldDirGammaTotArPbPb2040Plot){
            graphInvYieldDirGammaTotArPbPb2040Plot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPb2040Plot);
        }
        if (graphInvYieldDirGammaSysPbPb4080Plot){
            graphInvYieldDirGammaSysPbPb4080Plot->Draw("E2same");
        }
        if (graphInvYieldDirGammaStatPbPb4080Plot){
            graphInvYieldDirGammaStatPbPb4080Plot->Draw("p,E1Z,same");
        }
        if (graphInvYieldDirGammaTotArPbPb4080Plot){
            graphInvYieldDirGammaTotArPbPb4080Plot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArPbPb4080Plot);
        }
        if (graphInvYieldDirGammaTotArpp2760GeVScaled0020){
            graphInvYieldDirGammaTotArpp2760GeVScaled0020->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeVScaled0020);
        }
        if (graphInvYieldDirGammaTotArpp2760GeVScaled2040){
            graphInvYieldDirGammaTotArpp2760GeVScaled2040->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeVScaled2040);
        }
        if (graphInvYieldDirGammaTotArpp2760GeVScaled4080){
            graphInvYieldDirGammaTotArpp2760GeVScaled4080->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaTotArpp2760GeVScaled4080);
        }
        if (graphInvYieldDirGammaSyspp2760GeVScaled0020){
            graphInvYieldDirGammaSyspp2760GeVScaled0020->Draw("E2,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaSyspp2760GeVScaled0020);
        }
        if (graphInvYieldDirGammaSyspp2760GeVScaled2040){
            graphInvYieldDirGammaSyspp2760GeVScaled2040->Draw("E2,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaSyspp2760GeVScaled2040);
        }
        if (graphInvYieldDirGammaSyspp2760GeVScaled4080){
            graphInvYieldDirGammaSyspp2760GeVScaled4080->Draw("E2,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaSyspp2760GeVScaled4080);
        }
        if (graphInvYieldDirGammaStatpp2760GeVScaled0020){
            graphInvYieldDirGammaStatpp2760GeVScaled0020->Draw("p,E1Z,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaStatpp2760GeVScaled0020);
        }
        if (graphInvYieldDirGammaStatpp2760GeVScaled2040){
            graphInvYieldDirGammaStatpp2760GeVScaled2040->Draw("p,E1Z,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaStatpp2760GeVScaled2040);
        }
        if (graphInvYieldDirGammaStatpp2760GeVScaled4080){
            graphInvYieldDirGammaStatpp2760GeVScaled4080->Draw("p,E1Z,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphInvYieldDirGammaStatpp2760GeVScaled4080);
        }

        // TLegend* legendDirGammaPP2Split = GetAndSetLegend2(0.21, 0.10+(4*textsizeLabelsDirGamma*0.85), 0.21+0.21, 0.10+(8*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, collisionSystempp2760GeV.Data(), 42, 0.25);
        TLatex *labelEnergyScaledPP = new TLatex(0.22, 0.095+0.04*21, "ALICE");
        SetStyleTLatex( labelEnergyScaledPP, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        labelEnergyScaledPP->Draw();
        TLatex *labelEnergyScaledPP2 = new TLatex(0.22, 0.10+0.04*3, "pp #times <#it{N}_{coll}> , #sqrt{#it{s}} = 2.76 TeV");
        SetStyleTLatex( labelEnergyScaledPP2, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 11);
        labelEnergyScaledPP2->Draw();
        TLegend* legendDirGammaPP2Split = GetAndSetLegend2(0.21, 0.10, 0.21+0.22, 0.10+(3*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, "", 42, 0.25);
        legendDirGammaPP2Split->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled0020,Form("#gamma_{dir} #times %1.1f",nColl0020),"pf");
        legendDirGammaPP2Split->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled2040,Form("#gamma_{dir} #times %1.1f",nColl2040),"pf");
        legendDirGammaPP2Split->AddEntry(graphInvYieldDirGammaSyspp2760GeVScaled4080,Form("#gamma_{dir} #times %1.1f",nColl4080),"pf");
//         legendDirGammaPP2->AddEntry((TObject*)0,Form("ALICE x %1.1f",nColl0020),"");
//         legendDirGammaPP2->AddEntry((TObject*)0,Form("ALICE x %1.1f",nColl2040),"");
//         legendDirGammaPP2->AddEntry((TObject*)0,Form("ALICE x %1.1f",nColl4080),"");
        legendDirGammaPP2Split->Draw();
        
        // TLegend* legendDirGammaSplit = GetAndSetLegend2(0.21, 0.10, 0.21+0.21, 0.10+(4*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, collisionSystemPbPb760GeV.Data(), 42, 0.25);
        TLatex *labelEnergyScaledPBPB2 = new TLatex(0.95, 0.095+0.04*21, Form("%s",collisionSystemPbPb760GeV.Data()));
        SetStyleTLatex( labelEnergyScaledPBPB2, textSizeLabelsPixelDirGam,4, 1, 43, kTRUE, 31);
        labelEnergyScaledPBPB2->Draw();
        TLegend* legendDirGammaSplit = GetAndSetLegend2(0.71, 0.10+(21*textsizeLabelsDirGamma*0.85), 0.71+0.21, 0.10+(24*textsizeLabelsDirGamma*0.85) ,0.85*textsizeLabelsDirGamma, 1, "", 42, 0.25);
        legendDirGammaSplit->AddEntry(graphInvYieldDirGammaSysPbPb0020,"  0-20% #gamma_{dir}","pf");
        legendDirGammaSplit->AddEntry(graphInvYieldDirGammaSysPbPb2040,"20-40% #gamma_{dir}","pf");
        legendDirGammaSplit->AddEntry(graphInvYieldDirGammaSysPbPb4080,"40-80% #gamma_{dir}","pf");
        legendDirGammaSplit->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectra_WithScaledPP_Split.%s",outputDir.Data(),suffix.Data()));
    canvasDirGamma->Print(Form("%s/DirGammaSpectra_WithScaledPP_Split.pdf",outputDir.Data()));


    //*******************************************************************************************************************************************
    //******************************************* DR plot with pp and PbPb measurements *********************************************************
    //*******************************************************************************************************************************************
    TGraphAsymmErrors* graphDRStatPbPb0020Plot = (TGraphAsymmErrors*)graphDRStatPbPb0020->Clone("graphDRStatPbPb0020Plot");
    ProduceGraphAsymmWithoutXErrors(graphDRStatPbPb0020Plot);
    TGraphAsymmErrors* graphDRStatPbPb2040Plot = (TGraphAsymmErrors*)graphDRStatPbPb2040->Clone("graphDRStatPbPb2040Plot");
    ProduceGraphAsymmWithoutXErrors(graphDRStatPbPb2040Plot);
    TGraphAsymmErrors* graphDRStatPbPb4080Plot = (TGraphAsymmErrors*)graphDRStatPbPb4080->Clone("graphDRStatPbPb4080Plot");
    ProduceGraphAsymmWithoutXErrors(graphDRStatPbPb4080Plot);
    TGraphAsymmErrors* graphDRStatpp2760GeVPlot = (TGraphAsymmErrors*)graphDRStatpp2760GeV->Clone("graphDRStatpp2760GeVPlot");
    ProduceGraphAsymmWithoutXErrors(graphDRStatpp2760GeVPlot);

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

    //_______________________________________________________________ 0-20% dummy upper panel ___________________________________________________
    TH2D *dummyDR1 ;
    dummyDR1 = new TH2D("dummyDR1", "dummyDR1", 1000, 0., 22, 1000., doubleRatio[0], doubleRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyDR1, "#it{p}_{T} (GeV/#it{c})", "",
                            0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.95,0.10/(textsizeFacPad1*margin), 510, 505);
    dummyDR1->GetXaxis()->SetLabelOffset(-0.015);
    dummyDR1->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDR1->GetXaxis()->SetTickLength(0.06);
    dummyDR1->GetYaxis()->SetTickLength(0.028);

    //_______________________________________________________________ 20-40% dummy middle panel _________________________________________________
    TH2D *dummyDR2 ;
    dummyDR2 = new TH2D("dummyDR2", "dummyDR2", 1000, 0., 22, 1000., doubleRatio[0], doubleRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyDR2, "#it{p}_{T} (GeV/#it{c})", "#it{R}_{#gamma}", // = (#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})
                            0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.95,0.10/(textsizeFacPad2*margin), 510, 505);
    dummyDR2->GetXaxis()->SetLabelOffset(-0.015);
    dummyDR2->GetYaxis()->CenterTitle(kTRUE);
    dummyDR2->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDR2->GetXaxis()->SetTickLength(0.06);
    dummyDR2->GetYaxis()->SetTickLength(0.028);

    //_______________________________________________________________ 40-80% dummy lower panel __________________________________________________
    TH2D *dummyDR3 ;
    dummyDR3 = new TH2D("dummyDR3", "dummyDR3", 1000, 0., 22, 1000., doubleRatio[0], doubleRatio[1]);
    SetStyleHistoTH2ForGraphs( dummyDR3, "#it{p}_{T} (GeV/#it{c})", "",
                            0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.92,0.10/(textsizeFacPad3*margin), 510, 505);
    dummyDR3->GetXaxis()->SetLabelOffset(-0.015);
    dummyDR3->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDR3->GetXaxis()->SetTickLength(0.055);
    dummyDR3->GetYaxis()->SetTickLength(0.035);

    //_______________________________________________________________ 0-20% panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphDRSysPbPb0020, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);
        graphDRSysPbPb0020->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRSyspp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRSyspp2760GeV->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatPbPb0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);
        graphDRStatPbPb0020Plot->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatpp2760GeVPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");

        TLegend* legendDR0020 = GetAndSetLegend2(0.12, 0.9-(3*textsizeLabelsPad1), 0.12+0.21, 0.9,textsizeLabelsPad1, 1,
                                                 textALICE, 42, 0.3);
        legendDR0020->AddEntry(graphDRSysPbPb0020,Form("%s %s", "0-20%", collisionSystemPbPb760GeV.Data()),"pf");
        legendDR0020->AddEntry(graphDRSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        legendDR0020->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphDRSysPbPb2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphDRSysPbPb2040->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRSyspp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRSyspp2760GeV->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatPbPb2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        graphDRStatPbPb2040Plot->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatpp2760GeVPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");

        TLegend* legendDR2040 = GetAndSetLegend2(0.12, 0.94-(2*textsizeLabelsPad2), 0.12+0.21, 0.94,textsizeLabelsPad2, 1,
                                                 "", 42, 0.3);
        legendDR2040->AddEntry(graphDRSysPbPb2040,Form("%s %s", "20-40%", collisionSystemPbPb760GeV.Data()),"pf");
        legendDR2040->AddEntry(graphDRSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        legendDR2040->Draw();

    //_______________________________________________________________ 40-80% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);

        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphDRSysPbPb4080, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);
        graphDRSysPbPb4080->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRSyspp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRSyspp2760GeV->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatPbPb4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);
        graphDRStatPbPb4080Plot->Draw("p,E1Z,same");
        DrawGammaSetMarkerTGraphAsym(graphDRStatpp2760GeVPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
        graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");

        TLegend* legendDR4080 = GetAndSetLegend2(0.12, 0.94-(2*textsizeLabelsPad3), 0.12+0.21, 0.94,textsizeLabelsPad3, 1,
                                                      "", 42, 0.3);
        legendDR4080->AddEntry(graphDRSysPbPb4080,Form("%s %s","40-80%",collisionSystemPbPb760GeV.Data()),"pf");
        legendDR4080->AddEntry(graphDRSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
        legendDR4080->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_PbPbWithRefPP.%s", outputDir.Data(), suffix.Data()));
    canvasRatioIndDR->SaveAs(Form("%s/DR_PbPbWithRefPP.pdf", outputDir.Data()));
    
    
    
    TGraphAsymmErrors* graphDRStatpp8TeVPlot = (TGraphAsymmErrors*)graphDRStatpp8TeV->Clone("graphDRStatpp8TeVPlot");
    ProduceGraphAsymmWithoutXErrors(graphDRStatpp8TeVPlot);
    // double ratio combined
    TCanvas *canvasDoubleRatio = new TCanvas("canvasDoubleRatio","",0.095,0.09,1000,815);
    DrawGammaCanvasSettings( canvasDoubleRatio, 0.086, 0.01, 0.01, 0.105);
    canvasDoubleRatio->cd();
    canvasDoubleRatio->SetLogx();
    
    widthLinesBoxes                         = 1.4;
    Double_t textSizeSinglePad               = 0.05;
    TH2F * hist2DDRDummySingle       = new TH2F("hist2DDRDummySingle","hist2DDRDummySingle",1000,doubleRatioX[0], doubleRatioX[1],1000,0.72, 1.55);
    SetStyleHistoTH2ForGraphs(hist2DDRDummySingle, "#it{p}_{T} (GeV/#it{c})","#it{R}_{#gamma}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
    hist2DDRDummySingle->GetXaxis()->SetNoExponent();
    hist2DDRDummySingle->GetXaxis()->SetMoreLogLabels(kTRUE);
    hist2DDRDummySingle->DrawCopy();

      TLegend* legendDRSingle = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*3.3,0.43,0.953, textSizeSinglePad, 1, "ALICE", 42, 0.3);
      DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

      DrawGammaSetMarkerTGraphAsym(graphDRSyspp2760GeV, markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV , colorEpp2760GeV,widthLinesBoxes, kTRUE);
      graphDRSyspp2760GeV->Draw("E2same");
      legendDRSingle->AddEntry(graphDRSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");
  
      DrawGammaSetMarkerTGraphAsym(graphDRSyspp8TeV, markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV , colorEpp8TeV,widthLinesBoxes, kTRUE);
      graphDRSyspp8TeV->Draw("E2same");
      legendDRSingle->AddEntry(graphDRSyspp8TeV,collisionSystempp8TeV.Data(),"pf");

      DrawGammaSetMarkerTGraphAsym(graphDRStatpp2760GeVPlot,  markerStyleEpp2760GeV, markerSizeEpp2760GeV, colorEpp2760GeV , colorEpp2760GeV);
      graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");
      DrawGammaSetMarkerTGraphAsym(graphDRStatpp8TeVPlot,  markerStyleEpp8TeV, markerSizeEpp8TeV, colorEpp8TeV , colorEpp8TeV);
      graphDRStatpp8TeVPlot->Draw("pp,E1Z,same");

      legendDRSingle->Draw();
      hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_PP.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_PP.pdf", outputDir.Data()));


    TFile* fileTheory                               = new TFile("ExternalInput/Theory/TheoryCompilationPP.root");
    TGraph* graphTheoryNLODRpp7TeVPromptThermalLiuWerner  = (TGraph*) fileTheory->Get("DirectPhoton/graphRGammaThermalAndPromptDirectPhotonLiuWerner_pp7TeV_ALICECocktail");
    while(graphTheoryNLODRpp7TeVPromptThermalLiuWerner->GetX()[graphTheoryNLODRpp7TeVPromptThermalLiuWerner->GetN()-1] > 20.) graphTheoryNLODRpp7TeVPromptThermalLiuWerner->RemovePoint(graphTheoryNLODRpp7TeVPromptThermalLiuWerner->GetN()-1);

    hist2DDRDummySingle->DrawCopy();

          TLegend* legendDRSingle2 = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*5.3,0.43,0.953, textSizeSinglePad, 1, "ALICE", 42, 0.3);
          DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

          graphDRSyspp2760GeV->Draw("E2same");
          legendDRSingle2->AddEntry(graphDRSyspp2760GeV,collisionSystempp2760GeV.Data(),"pf");

          graphDRSyspp8TeV->Draw("E2same");
          legendDRSingle2->AddEntry(graphDRSyspp8TeV,collisionSystempp8TeV.Data(),"pf");

          graphDRStatpp2760GeVPlot->Draw("p,E1Z,same");
          graphDRStatpp8TeVPlot->Draw("pp,E1Z,same");

          if (graphTheoryNLODRpp7TeVPromptThermalLiuWerner){
              DrawGammaNLOTGraph( graphTheoryNLODRpp7TeVPromptThermalLiuWerner, 2, 5, kPink+2 );
              graphTheoryNLODRpp7TeVPromptThermalLiuWerner->Draw("lc,same");
              legendDRSingle2->AddEntry(graphTheoryNLODRpp7TeVPromptThermalLiuWerner,"NLO pQCD, #scale[0.75]{Prompt + Thermal}","l");
              legendDRSingle2->AddEntry((TObject*)0,"for pp, #sqrt{#it{s}} = 7 TeV","");
          }

          legendDRSingle2->Draw();
          hist2DDRDummySingle->Draw("same,axis");

        canvasDoubleRatio->Print(Form("%s/DR_PP_wThermal.%s", outputDir.Data(), suffix.Data()));
        canvasDoubleRatio->Print(Form("%s/DR_PP_wThermal.pdf", outputDir.Data()));
}
