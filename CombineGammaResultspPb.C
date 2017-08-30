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
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"


void CombineGammaResultspPb(    TString inputFileNamePCM    = "",
                                TString inputFileNamePHOS   = "",
                                TString inputFileNameEMC    = "",
                                TString inputFileNamePCMEMC = "",
                                TString suffix              = "eps",
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
    TString outputDir                                           = Form("%s/%s/CombineGammaMeasurementspPb",suffix.Data(),dateForOutput.Data());
    TString fileNameTheorypPb                                   = "ExternalInputpPb/Theory/TheoryCompilationPPb.root";

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPHOSGammapPb.root", inputFileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMGammapPb.root", inputFileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEMCGammapPb.root", inputFileNamePCMEMC.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCGammapPb.root", inputFileNameEMC.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/TheorypPb.root", fileNameTheorypPb.Data(), outputDir.Data()));

    TString nameFinalResDat                                     = Form("%s/CombinedResultsGamma_FitResults_%s.dat",outputDir.Data(),dateForOutput.Data());
    fstream fileFinalResults;
    fileFinalResults.open(nameFinalResDat.Data(), ios::out);

    //*******************************************************************************************************************************************
    //******************************************************* set ranges for plotting ***********************************************************
    //*******************************************************************************************************************************************
    Double_t doubleRatio[2];
    Double_t indMeasRatio[2];
//     Double_t incRatio[2];
    Double_t doubleRatioX[2];
    Double_t doubleRatioXpp[2];
    doubleRatio[0]      = 0.75;     doubleRatio[1]      = 1.65;
    indMeasRatio[0]     = 0.65;     indMeasRatio[1]     = 1.45;
//  incRatio[0]         = 0.0;      incRatio[1]         = 1.7;
    doubleRatioX[0]     = 0.7;      doubleRatioX[1]     = 16;
    doubleRatioXpp[0]   = 0.23;     doubleRatioXpp[1]   = 30;

    Color_t colorCocktailPi0                        = kRed+2;
    Color_t colorCocktailEta                        = kBlue+1;
    Color_t colorCocktailEtaP                       = kOrange+1;
    Color_t colorCocktailOmega                      = kYellow+2;
    Color_t colorCocktailPhi                        = kViolet;
    Color_t colorCocktailRho0                       = kAzure-2;
    Color_t colorCocktailSigma0                     = kGray+1;

    Color_t colorNLOcalc                            = kBlue-7; // kBlack
    Style_t fillStyleNLO                            = 1001;
    Color_t colorEPS09calc                          = kGray;
    Color_t colorEPS09calc2                         = kGray+1;
    Style_t fillStyleEPS09                          = 3008;
    Color_t colorCT10calc                           = kAzure-4;
    Style_t fillStyleCT10                           = 3001;
    Color_t colorNLOMcGill                          = kOrange+1;
    Style_t styleNLOMcGill                          = 7;
    Color_t colorPHSD                               = kGray+2;
    Style_t stylePHSD                               = 2;
    Color_t colorChatterjee                         = kBlue+1;
    Style_t styleChatterjee                         = 4;
    Color_t colorHees                               = kBlack;
    Style_t styleHees                               = 3;
    Color_t colorHe                                 = kBlue-7;
    Style_t styleHe                                 = 8;
    Style_t styleFit                                = 1;

    Color_t colorPHENIX                             = kGray+2;
    Style_t markerStylePHENIX                       = 24;
    Size_t markerSizePHENIX                         = 2;

    TString  nameMeasGlobal[11]                     = { "PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal",
                                                        "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged",
                                                        "PCMOtherDataset"};
    TString  nameMeasGlobalLabel[11]                = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC",
                                                        "PCM-Dalitz", "PHOS-Dalitz", "EMC-Dalitz", "EMC high pT", "mEMC",
                                                        "PCMOtherDataset"};


    Color_t  colorDet[11];
    Color_t  colorDetMC[11];
    Style_t  markerStyleDet[11];
    Style_t  markerStyleDetMC[11];
    Size_t   markerSizeDet[11];
    Size_t   markerSizeDetMC[11];


    for (Int_t i = 0; i < 11; i++){
        colorDet[i]                                 = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
        colorDetMC[i]                               = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
        markerStyleDet[i]                           = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
        markerStyleDetMC[i]                         = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
        markerSizeDet[i]                            = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
        markerSizeDetMC[i]                          = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
    }

    Color_t colorCombpPb                            = GetColorDefaultColor("pPb_5.023TeV", "", "");
    Color_t colorCombpPbBox                         = GetColorDefaultColor("pPb_5.023TeV", "", "", kTRUE);
    Style_t markerStyleCombpPb                      = GetDefaultMarkerStyle("pPb_5.023TeV", "", "");
    Size_t markerSizeCombpPb                        = GetDefaultMarkerSize("pPb_5.023TeV", "", "");

    Width_t widthLinesBoxes                         = 1.4;
    Width_t widthCommonFit                          = 2.4;

    TString collisionSystempPb                      = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";


    TH1D* histoDRPi0FitStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors*  graphDRPi0FitSysErr[11]     = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoIncGammaStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors* graphIncGammaSysErr[11]      = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from pPb PCM file ******************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMGammapPb                          = new TFile( inputFileNamePCM.Data());
    //________________________________________________ Load PCM pPb 5.023TeV _________________________________________________________________________
    TDirectory* directoryPCMGammapPb                = (TDirectory*)filePCMGammapPb->Get("Gamma_pPb5TeV");
        histoDRPi0FitStatErr[0]                         = (TH1D*) directoryPCMGammapPb->Get("DoubleRatioPi0FitStatError");
        graphDRPi0FitSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("DoubleRatioPi0FitSystError");
        histoIncGammaStatErr[0]                         = (TH1D*) directoryPCMGammapPb->Get("IncRatioStatError");
        graphIncGammaSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncRatioSystError");

    //*******************************************************************************************************************************************
    //*********************************************** Load PCMEMC histograms from pPb PCM file ******************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMEMCGammapPb                       = new TFile( inputFileNamePCMEMC.Data());
    //________________________________________________ Load PCM-EMC pPb 5.023TeV _________________________________________________________________________
    TDirectory* directoryPCMEMCGammapPb             = (TDirectory*)filePCMEMCGammapPb->Get("Gamma_pPb5TeV");
        histoDRPi0FitStatErr[4]                         = (TH1D*) directoryPCMEMCGammapPb->Get("DoubleRatioPi0FitStatError");
        graphDRPi0FitSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("DoubleRatioPi0FitSystError");
        histoIncGammaStatErr[4]                         = (TH1D*) directoryPCMEMCGammapPb->Get("IncRatioStatError");
        graphIncGammaSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("IncRatioSystError");

    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from pPb EMC file ******************************************************
    //*******************************************************************************************************************************************
    TFile* fileEMCGammapPb                          = new TFile( inputFileNameEMC.Data());
    //________________________________________________ Load EMC pPb 5.023TeV _________________________________________________________________________
    TDirectory* directoryEMCGammapPb                = (TDirectory*)fileEMCGammapPb->Get("Gamma_pPb5TeV");
        histoDRPi0FitStatErr[2]                         = (TH1D*) directoryEMCGammapPb->Get("DoubleRatioPi0FitStatError");
        if (histoDRPi0FitStatErr[2]) histoDRPi0FitStatErr[2]->GetXaxis()->SetRangeUser(2.5,14);
        graphDRPi0FitSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("DoubleRatioPi0FitSystError");
        histoIncGammaStatErr[2]                         = (TH1D*) directoryEMCGammapPb->Get("IncRatioStatError");
        graphIncGammaSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("IncRatioSystError");


    //*******************************************************************************************************************************************
    //*********************************************** Load PHOS histograms from PHOS file *******************************************************
    //*******************************************************************************************************************************************
    TFile* filePHOSGamma                            = new TFile( inputFileNamePHOS.Data());
        histoDRPi0FitStatErr[1]                         = (TH1D*) filePHOSGamma->Get("hGamma_PbPb_cen6_Stat");
        TH1D* histoPHOSDRPi0FitSysErrpPb                = (TH1D*) filePHOSGamma->Get("hGamma_PbPb_cen6_SystRatio");
        graphDRPi0FitSysErr[1]                          = new TGraphAsymmErrors(histoPHOSDRPi0FitSysErrpPb);

    //*******************************************************************************************************************************************
    //************************************************ Load theory curves from external input ***************************************************
    //*******************************************************************************************************************************************
    TFile* fileTheory                                       = new TFile( fileNameTheorypPb.Data());
        TGraphAsymmErrors* graphTheoryNLODRpPb                  = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10_ALICECocktail");
        TGraph* graphTheoryNLODRpPbCenter                       = (TGraph*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10_ALICECocktail_Center");
        TGraphAsymmErrors* graphTheoryNLOpPb                    = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphDirectPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10");
        TGraphAsymmErrors* graphTheoryMCGillDRpPb               = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonSpecMcGill5023GeV_ALICECocktail");
        TGraph* graphTheoryMCGillDRpPbCenter                    = (TGraph*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonSpecMcGill5023GeV_ALICECocktail_Center");
        TGraphAsymmErrors* graphTheoryMCGillpPb                 = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphDirectPhotonSpecMcGill5023GeV");



    // double ratio combined
    TCanvas *canvasDoubleRatio = new TCanvas("canvasDoubleRatio","",0.095,0.09,1000,815);
    DrawGammaCanvasSettings( canvasDoubleRatio, 0.086, 0.01, 0.01, 0.09);
    canvasDoubleRatio->cd();
    canvasDoubleRatio->SetLogx();

    Double_t minY                            = 0.85;
    Double_t maxY                            = 1.65;
    Double_t textSizeSinglePad               = 0.05;
    TH2F * hist2DDRDummySingle       = new TH2F("hist2DDRDummySingle","hist2DDRDummySingle",1000,doubleRatioXpp[0], doubleRatioXpp[1],1000,doubleRatio[0], doubleRatio[1]);
    SetStyleHistoTH2ForGraphs(hist2DDRDummySingle, "#it{p}_{T} (GeV/#it{c})","#it{R}_{#gamma}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.75,0.8);
    hist2DDRDummySingle->GetXaxis()->SetLabelOffset(-1e-2);
    hist2DDRDummySingle->DrawCopy();

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        TLegend* legendDRPCMNLOpPbSingle = GetAndSetLegend2(0.12,0.90-1.05*0.85*textSizeSinglePad*5,0.5,0.90, 0.85*textSizeSinglePad, 2, "", 42, 0.25);
        if (graphTheoryNLODRpPb) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryNLODRpPb, 0, 0, kAzure-9, kAzure-9, 0.2, kTRUE, kAzure-9);
            graphTheoryNLODRpPb->Draw("3,same");
            legendDRPCMNLOpPbSingle->AddEntry(graphTheoryNLODRpPb,"NLO pQCD PDF: CT10 FF: GRV ","f");
            legendDRPCMNLOpPbSingle->AddEntry((TObject*)0,"","");
        }
        if (graphTheoryNLODRpPbCenter){
            DrawGammaNLOTGraph( graphTheoryNLODRpPbCenter, 2, 7, kAzure+2);
            graphTheoryNLODRpPbCenter->Draw("lc,same");
        }
        if (graphTheoryMCGillDRpPb) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryMCGillDRpPb, 0, 0, colorNLOMcGill, colorNLOMcGill, 0.2, kTRUE, colorNLOMcGill, kTRUE);
            graphTheoryMCGillDRpPb->Draw("3,same");
            legendDRPCMNLOpPbSingle->AddEntry(graphTheoryMCGillDRpPb,"McGill","f");
            legendDRPCMNLOpPbSingle->AddEntry((TObject*)0,"","");
        }
        if (graphTheoryMCGillDRpPbCenter){
            DrawGammaNLOTGraph( graphTheoryMCGillDRpPbCenter, 2, 7, colorNLOMcGill);
            graphTheoryMCGillDRpPbCenter->Draw("lc,same");
        }

        for (Int_t i = 10; i > -1; i--){
            if (graphDRPi0FitSysErr[i]){
                DrawGammaSetMarkerTGraphAsym(graphDRPi0FitSysErr[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphDRPi0FitSysErr[i]->Draw("E2same");
                legendDRPCMNLOpPbSingle->AddEntry(graphDRPi0FitSysErr[i],nameMeasGlobalLabel[i],"pf");

            }
            if (histoDRPi0FitStatErr[i]){
                DrawGammaSetMarker(histoDRPi0FitStatErr[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoDRPi0FitStatErr[i]->Draw("p,same,e0,X0");
                if (!graphDRPi0FitSysErr[i])legendDRPCMNLOpPbSingle->AddEntry(histoDRPi0FitStatErr[i],nameMeasGlobalLabel[i],"p");
            }
        }

        if (histoDRPi0FitStatErr[4]) histoDRPi0FitStatErr[4]->Draw("p,same,e0,X0");
        legendDRPCMNLOpPbSingle->Draw();


        TLatex *labelALICEDRSingle = new TLatex(0.95,0.91,"ALICE this thesis");
        SetStyleTLatex( labelALICEDRSingle, 0.85*textSizeSinglePad,4, 1, 42, kTRUE, 31);
        labelALICEDRSingle->Draw();

        TLatex *labelDRCentpPbSingle = new TLatex(0.12,0.91,collisionSystempPb.Data());
        SetStyleTLatex( labelDRCentpPbSingle, 0.85*textSizeSinglePad,4, 1, 42, kTRUE, 11);
        labelDRCentpPbSingle->Draw();
    hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_PCMMeasurementTheory_pPb.%s", outputDir.Data(), suffix.Data()));
}