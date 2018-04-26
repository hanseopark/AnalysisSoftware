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



void CombineGammaResultspPb(    TString inputFileNamePCM        = "",
                                TString inputFileNamePHOS       = "",
                                TString inputFileNameEMC        = "",
                                TString inputFileNamePCMEMC     = "",
                                TString suffix                  = "eps",
                                TString fileNameCorrelations    = "",
                                Bool_t enablepValueCalc         = kFALSE,
                                Bool_t isThesis                 = kFALSE
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
    TString outputDirFit                                        = Form("%s/%s/CombineGammaMeasurementspPb/Pi0Fitted",suffix.Data(),dateForOutput.Data());
    TString fileNameTheorypPb                                   = "ExternalInputpPb/Theory/TheoryCompilationPPb.root";
    TString fileNameIsolatedGamma                               = "ExternalInputpPb/EMCAL/isolated_photon_invariant_yield_20180417.root";


    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDirFit);
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
    doubleRatioX[0]     = 0.23;     doubleRatioX[1]     = 52;
    doubleRatioXpp[0]   = 0.23;     doubleRatioXpp[1]   = 52;

    Color_t colorCocktailPi0                        = kRed+2;
    Color_t colorCocktailEta                        = kBlue+1;
    Color_t colorCocktailEtaP                       = kOrange+1;
    Color_t colorCocktailOmega                      = kYellow+2;
    Color_t colorCocktailPhi                        = kViolet;
    Color_t colorCocktailRho0                       = kAzure-2;
    Color_t colorCocktailSigma0                     = kGray+1;


    Color_t colorPHSD                               = kGray+2;
    Style_t stylePHSD                               = 2;
    Color_t  colorNLOWerner                         = kAzure+2;
    Color_t  colorNLOWernerBand                     = kAzure-9;
    Color_t  colorNLONCTEQ                          = kOrange+4;
    Color_t  colorNLONCTEQBand                      = kOrange+5;
    Color_t  colorNLOEPPS                           = kViolet+3;
    Color_t  colorNLOEPPSBand                       = kViolet-8;
    Color_t  colorNLOMcGill                         = kGreen+2;
    Color_t  colorNLOMcGillBand                     = kGreen-6;
    Style_t  styleMarkerNLOWerner                   = 24;
    Style_t  styleMarkerNLONCTEQ                    = 27;
    Style_t  styleMarkerNLOEPPS                     = 30;
    Style_t  styleLineNLOWerner                     = 5;
    Style_t  styleLineNLONCTEQ                      = 8;
    Style_t  styleLineNLOEPPS                       = 6;
    Style_t  styleLineMcGill                        = 7;
    Width_t  widthLineNLO                           = 2.;

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
    Color_t colorComb                               = GetColorDefaultColor("pPb_5.023TeV", "", "");
    Style_t markerStyleComb                         = GetDefaultMarkerStyle("pPb_5.023TeV", "", "");
    Size_t markerSizeComb                           = GetDefaultMarkerSize("pPb_5.023TeV", "", "");

    Color_t colorCombpPb                            = kBlack; //GetColorDefaultColor("pPb_5.023TeV", "", "");
    Color_t colorCombpPbBox                         = kGray+2; //GetColorDefaultColor("pPb_5.023TeV", "", "", kTRUE);
    Style_t markerStyleCombpPb                      = 20; //GetDefaultMarkerStyle("pPb_5.023TeV", "", "");
    Size_t markerSizeCombpPb                        = 1.8; //GetDefaultMarkerSize("pPb_5.023TeV", "", "");

    Width_t widthLinesBoxes                         = 1.4;
    Width_t widthCommonFit                          = 2.4;
    Double_t scalingToNSD                           = 0.964;

    TString collisionSystempPb                      = "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    TString collisionSystempPbNSD                   = "NSD p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";

    TString textALICE                               = "ALICE preliminary";
    TString textALICEPerf                           = "ALICE performance";
    if (isThesis) textALICE                         = "ALICE this thesis";

    cout << "Setting Gamma binning" << endl;
    Double_t xPtLimitsGamma[100]                    = {0};
    Int_t maxNBinsGammaAbs                          = 0;
    Int_t maxNBinsGamma                             = GetBinning( xPtLimitsGamma, maxNBinsGammaAbs, "Gamma", "pPb_5.023TeV", 23 );
    for (Int_t i = 0; i< maxNBinsGamma; i++){
        cout << i << ": "<< xPtLimitsGamma[i] <<" - " << xPtLimitsGamma[i+1]<< ", " <<endl;
    }

    TH1D* histoDRPi0FitStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors*  graphDRPi0FitSysErr[11]     = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoDRNonFitStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors*  graphDRNonFitSysErr[11]     = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoIncGammaRatioStatErr[11]             = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors* graphIncGammaRatioSysErr[11] = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoIncGammaStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors* graphIncGammaSysErr[11]      = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoRawGamma[11]                         = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoConvProb[11]                         = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoEffi[11]                             = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoEffiMCPt[11]                         = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoResolCorr[11]                        = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoPurity[11]                           = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoEffSecCorr[4][11]                    = { {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL} };
    TH1D* histoPileupCorr[11]                       = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from pPb PCM file ******************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMGammapPb                          = new TFile( inputFileNamePCM.Data());
    //________________________________________________ Load PCM pPb 5.023TeV _________________________________________________________________________
    TDirectory* directoryPCMGammapPb                = (TDirectory*)filePCMGammapPb->Get("Gamma_pPb5TeV");
        histoDRPi0FitStatErr[0]                         = (TH1D*) directoryPCMGammapPb->Get("DoubleRatioPi0FitStatError");
        graphDRPi0FitSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("DoubleRatioPi0FitSystError");
        histoDRNonFitStatErr[0]                         = (TH1D*) directoryPCMGammapPb->Get("DoubleRatioStatError");
        graphDRNonFitSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("DoubleRatioSystError");
        histoIncGammaRatioStatErr[0]                    = (TH1D*) directoryPCMGammapPb->Get("IncRatioStatError");
        graphIncGammaRatioSysErr[0]                     = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncRatioSystError");
        histoConvProb[0]                                = (TH1D*) directoryPCMGammapPb->Get("GammaConversionProbability");
        histoEffi[0]                                    = (TH1D*) directoryPCMGammapPb->Get("GammaRecoEfficiency");
        histoEffiMCPt[0]                                = (TH1D*) directoryPCMGammapPb->Get("GammaRecoEfficiencyMCPt");
        histoResolCorr[0]                               = (TH1D*) directoryPCMGammapPb->Get("GammaResolCorr");
        histoPurity[0]                                  = (TH1D*) directoryPCMGammapPb->Get("GammaTruePurity");
        histoIncGammaStatErr[0]                         = (TH1D*) directoryPCMGammapPb->Get("IncGammaStatError");
        if (histoIncGammaStatErr[0])
            histoIncGammaStatErr[0]->Scale(scalingToNSD);
        graphIncGammaSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncGammaSystError");
        if (graphIncGammaSysErr[0])
            graphIncGammaSysErr[0]                      = ScaleGraph(graphIncGammaSysErr[0],scalingToNSD);

        histoPileupCorr[0]                              = (TH1D*) directoryPCMGammapPb->Get("PileUpCorrectionFactor");
        histoRawGamma[0]                                = (TH1D*) directoryPCMGammapPb->Get("GammaRawYields");
        histoEffSecCorr[0][0]                           = (TH1D*) directoryPCMGammapPb->Get("GammaEffectiveSecondaryCorr_K0s");
        histoEffSecCorr[1][0]                           = (TH1D*) directoryPCMGammapPb->Get("GammaEffectiveSecondaryCorr_K0l");
        histoEffSecCorr[2][0]                           = (TH1D*) directoryPCMGammapPb->Get("GammaEffectiveSecondaryCorr_Lambda");
        histoEffSecCorr[3][0]                           = (TH1D*) directoryPCMGammapPb->Get("GammaEffectiveSecondaryCorr_Rest");
    //*******************************************************************************************************************************************
    //*********************************************** Load PCMEMC histograms from pPb PCM file ******************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMEMCGammapPb                       = new TFile( inputFileNamePCMEMC.Data());
    //________________________________________________ Load PCM-EMC pPb 5.023TeV _________________________________________________________________________
    TDirectory* directoryPCMEMCGammapPb             = (TDirectory*)filePCMEMCGammapPb->Get("Gamma_pPb5TeV");
        histoDRPi0FitStatErr[4]                         = (TH1D*) directoryPCMEMCGammapPb->Get("DoubleRatioPi0FitStatError");
        graphDRPi0FitSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("DoubleRatioPi0FitSystError");
        histoDRNonFitStatErr[4]                         = (TH1D*) directoryPCMEMCGammapPb->Get("DoubleRatioStatError");
        graphDRNonFitSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("DoubleRatioSystError");
        histoIncGammaRatioStatErr[4]                    = (TH1D*) directoryPCMEMCGammapPb->Get("IncRatioStatError");
        graphIncGammaRatioSysErr[4]                     = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("IncRatioSystError");
        histoConvProb[4]                                = (TH1D*) directoryPCMEMCGammapPb->Get("GammaConversionProbability");
        histoEffi[4]                                    = (TH1D*) directoryPCMEMCGammapPb->Get("GammaRecoEfficiency");
        histoEffiMCPt[4]                                = (TH1D*) directoryPCMEMCGammapPb->Get("GammaRecoEfficiencyMCPt");
        histoResolCorr[4]                               = (TH1D*) directoryPCMEMCGammapPb->Get("GammaResolCorr");
        histoPurity[4]                                  = (TH1D*) directoryPCMEMCGammapPb->Get("GammaTruePurity");
        histoIncGammaStatErr[4]                         = (TH1D*) directoryPCMEMCGammapPb->Get("IncGammaStatError");
        if (histoIncGammaStatErr[4])
            histoIncGammaStatErr[4]->Scale(scalingToNSD);
        graphIncGammaSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("IncGammaSystError");
        if (graphIncGammaSysErr[4])
            graphIncGammaSysErr[4]                      = ScaleGraph(graphIncGammaSysErr[4],scalingToNSD);
        histoPileupCorr[4]                              = (TH1D*) directoryPCMEMCGammapPb->Get("PileUpCorrectionFactor");
        if (histoPileupCorr[4]) histoPileupCorr[4]->GetXaxis()->SetRangeUser(0.8,14);
        histoRawGamma[4]                                = (TH1D*) directoryPCMEMCGammapPb->Get("GammaRawYields");
        histoEffSecCorr[0][4]                           = (TH1D*) directoryPCMEMCGammapPb->Get("GammaEffectiveSecondaryCorr_K0s");
        histoEffSecCorr[1][4]                           = (TH1D*) directoryPCMEMCGammapPb->Get("GammaEffectiveSecondaryCorr_K0l");
        histoEffSecCorr[2][4]                           = (TH1D*) directoryPCMEMCGammapPb->Get("GammaEffectiveSecondaryCorr_Lambda");
        histoEffSecCorr[3][4]                           = (TH1D*) directoryPCMEMCGammapPb->Get("GammaEffectiveSecondaryCorr_Rest");
        //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from pPb EMC file ******************************************************
    //*******************************************************************************************************************************************
    TFile* fileEMCGammapPb                          = new TFile( inputFileNameEMC.Data());
    //________________________________________________ Load EMC pPb 5.023TeV _________________________________________________________________________
    TDirectory* directoryEMCGammapPb                = (TDirectory*)fileEMCGammapPb->Get("Gamma_pPb5TeV");
        histoDRPi0FitStatErr[2]                         = (TH1D*) directoryEMCGammapPb->Get("DoubleRatioPi0FitStatError");
        graphDRPi0FitSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("DoubleRatioPi0FitSystError");
        histoDRNonFitStatErr[2]                         = (TH1D*) directoryEMCGammapPb->Get("DoubleRatioStatError");
        graphDRNonFitSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("DoubleRatioSystError");
        histoIncGammaRatioStatErr[2]                    = (TH1D*) directoryEMCGammapPb->Get("IncRatioStatError");
        graphIncGammaRatioSysErr[2]                     = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("IncRatioSystError");
        histoEffi[2]                                    = (TH1D*) directoryEMCGammapPb->Get("GammaRecoEfficiency");
        histoEffiMCPt[2]                                = (TH1D*) directoryEMCGammapPb->Get("GammaRecoEfficiencyMCPt");
        histoResolCorr[2]                               = (TH1D*) directoryEMCGammapPb->Get("GammaResolCorr");
        histoPurity[2]                                  = (TH1D*) directoryEMCGammapPb->Get("GammaTruePurity");
        histoIncGammaStatErr[2]                         = (TH1D*) directoryEMCGammapPb->Get("IncGammaStatError");
        if (histoIncGammaStatErr[2])
            histoIncGammaStatErr[2]->Scale(scalingToNSD);
        graphIncGammaSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("IncGammaSystError");
        if (graphIncGammaSysErr[2])
            graphIncGammaSysErr[2]                      = ScaleGraph(graphIncGammaSysErr[2],scalingToNSD);
        histoRawGamma[2]                                = (TH1D*) directoryEMCGammapPb->Get("GammaRawYields");
        histoEffSecCorr[0][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_K0s");
        histoEffSecCorr[1][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_K0l");
        histoEffSecCorr[2][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_Lambda");
        histoEffSecCorr[3][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_Rest");

        Int_t j = 1;
        while (histoEffi[2]->GetBinCenter(j) < 2.05){
          histoEffi[2]->SetBinContent(j,1e-10);
          histoEffiMCPt[2]->SetBinContent(j,1e-10);
          histoResolCorr[2]->SetBinContent(j,1e10);
          histoPurity[2]->SetBinContent(j,1e10);
          histoIncGammaStatErr[2]->SetBinContent(j,0);
          histoDRPi0FitStatErr[2]->SetBinContent(j,0);
          histoDRNonFitStatErr[2]->SetBinContent(j,0);
          histoIncGammaRatioStatErr[2]->SetBinContent(j,0);
          histoEffSecCorr[0][2]->SetBinContent(j,1e10);
          histoEffSecCorr[1][2]->SetBinContent(j,1e10);
          histoEffSecCorr[2][2]->SetBinContent(j,1e10);
          histoEffSecCorr[3][2]->SetBinContent(j,1e10);
          j++;
        }
        while (graphDRPi0FitSysErr[2]->GetX()[0] < 2.05) graphDRPi0FitSysErr[2]->RemovePoint(0);
        while (graphDRNonFitSysErr[2]->GetX()[0] < 2.05) graphDRNonFitSysErr[2]->RemovePoint(0);
        while (graphIncGammaRatioSysErr[2]->GetX()[0] < 2.05) graphIncGammaRatioSysErr[2]->RemovePoint(0);
        while (graphIncGammaSysErr[2]->GetX()[0] < 2.05) graphIncGammaSysErr[2]->RemovePoint(0);
    //*******************************************************************************************************************************************
    //*********************************************** Load PHOS histograms from PHOS file *******************************************************
    //*******************************************************************************************************************************************
    TFile* filePHOSGamma                            = new TFile( inputFileNamePHOS.Data());
    TDirectory* directoryPHOSGammapPb                   = (TDirectory*)filePHOSGamma->Get("PHOS_pPb_5020_Centrality_00-100");
    histoDRPi0FitStatErr[1]                         = (TH1D*) directoryPHOSGammapPb->Get("hPHOS_DoubleRatio_pPb_cen00-100_FitPi0_Stat");
        TH1D* histoPHOSDRPi0FitSysErrpPb                = (TH1D*) directoryPHOSGammapPb->Get("hPHOS_DoubleRatio_pPb_cen00-100_FitPi0_Syst");
        graphDRPi0FitSysErr[1]                          = new TGraphAsymmErrors(histoPHOSDRPi0FitSysErrpPb);
        histoDRNonFitStatErr[1]                         = (TH1D*) directoryPHOSGammapPb->Get("hPHOS_DoubleRatio_pPb_cen00-100_Stat");
        TH1D* histoPHOSDRNonFitSysErrpPb                = (TH1D*) directoryPHOSGammapPb->Get("hPHOS_DoubleRatio_pPb_cen00-100_Syst");
        graphDRNonFitSysErr[1]                          = new TGraphAsymmErrors(histoPHOSDRNonFitSysErrpPb);
        TF1* fitEffiTimesAccPHOSINT7                    = (TF1*) directoryPHOSGammapPb->Get("hPHOS_gammaIncl_pPb_cen00-100_Efficiency");
        TF1* fitEffiTimesAccPHOSPHI7                    = (TF1*) directoryPHOSGammapPb->Get("hPHOS_gammaIncl_pPb_PHI7_cen00-100_Efficiency");
        TF1* fitPurityPHOS                       = (TF1*) directoryPHOSGammapPb->Get("hPHOS_gammaIncl_pPb_cen00-100_Contam");
        histoPurity[1]                                  = (TH1D*) histoPHOSDRPi0FitSysErrpPb->Clone("GammaTruePurity");
        histoEffi[1]                                    = (TH1D*) histoPHOSDRPi0FitSysErrpPb->Clone("GammaRecoEfficiency");
        TH1D* histoEffiPHI7                             = (TH1D*) histoPHOSDRPi0FitSysErrpPb->Clone("GammaRecoEfficiencyPHI7");
        for (Int_t i = 1; i< histoPHOSDRPi0FitSysErrpPb->GetNbinsX()+1; i++){
            if (histoPHOSDRPi0FitSysErrpPb->GetBinContent(i) > 0){
                histoPurity[1]->SetBinContent(i,fitPurityPHOS->Eval(histoPurity[1]->GetBinCenter(i)));
                histoPurity[1]->SetBinError(i,0.005*fitPurityPHOS->Eval(histoPurity[1]->GetBinCenter(i)));
                histoEffi[1]->SetBinContent(i,(fitEffiTimesAccPHOSINT7->Eval(histoEffi[1]->GetBinCenter(i))/(0.125/0.5 *40./360)));
                histoEffi[1]->SetBinError(i,(0.005*fitEffiTimesAccPHOSINT7->Eval(histoEffi[1]->GetBinCenter(i))/(0.125/0.5 *40./360)));
                if (histoEffiPHI7->GetBinCenter(i) > 7){
                    histoEffiPHI7->SetBinContent(i,(fitEffiTimesAccPHOSPHI7->Eval(histoEffiPHI7->GetBinCenter(i))/(0.125/0.5 *40./360)));
                    histoEffiPHI7->SetBinError(i,(0.005*fitEffiTimesAccPHOSPHI7->Eval(histoEffiPHI7->GetBinCenter(i))/(0.125/0.5 *40./360)));
                } else {
                    histoEffiPHI7->SetBinContent(i,0);
                    histoEffiPHI7->SetBinError(i,0);
                }
            } else {
                histoPurity[1]->SetBinContent(i,0);
                histoPurity[1]->SetBinError(i,0);
                histoEffi[1]->SetBinContent(i,0);
                histoEffi[1]->SetBinError(i,0);
                histoEffiPHI7->SetBinContent(i,0);
                histoEffiPHI7->SetBinError(i,0);
            }
        }
        histoIncGammaStatErr[1]                         = (TH1D*) directoryPHOSGammapPb->Get("hPHOS_gammaIncl_pPb_cen00-100_Stat");
        if (histoIncGammaStatErr[1])
            histoIncGammaStatErr[1]->Scale(scalingToNSD);
        TH1D* histoPHOSIncGammaPi0FitSysErrpPb          = (TH1D*) directoryPHOSGammapPb->Get("hPHOS_gammaIncl_pPb_cen00-100_Syst");
        if (histoPHOSIncGammaPi0FitSysErrpPb)
            histoPHOSIncGammaPi0FitSysErrpPb->Scale(scalingToNSD);
        graphIncGammaSysErr[1]                          = new TGraphAsymmErrors(histoPHOSIncGammaPi0FitSysErrpPb);

    //*******************************************************************************************************************************************
    //************************************************ Load theory curves from external input ***************************************************
    //*******************************************************************************************************************************************
    TFile* fileTheory                               = new TFile( fileNameTheorypPb.Data());
        TGraphAsymmErrors* graphTheoryNLODRpPb          = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10_ALICECocktail");
        TGraph* graphTheoryNLODRpPbCenter               = (TGraph*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10_ALICECocktail_Center");
        while (graphTheoryNLODRpPb->GetX()[0] < 3)  graphTheoryNLODRpPb->RemovePoint(0);
        while (graphTheoryNLODRpPbCenter->GetX()[0] < 3)    graphTheoryNLODRpPbCenter->RemovePoint(0);
        TGraphAsymmErrors* graphTheoryNLOpPb            = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphDirectPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10");
        while (graphTheoryNLOpPb->GetX()[0] < 3)  graphTheoryNLOpPb->RemovePoint(0);
        TGraphAsymmErrors* graphTheoryNLOpPbPrompt      = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphPromptPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10");
        while (graphTheoryNLOpPbPrompt->GetX()[0] < 3.5)  graphTheoryNLOpPbPrompt->RemovePoint(0);

        TGraphAsymmErrors* graphTheoryMCGillDRpPb       = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonSpecMcGill5023GeV_ALICECocktail");
        TGraph* graphTheoryMCGillDRpPbCenter            = (TGraph*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonSpecMcGill5023GeV_ALICECocktail_Center");
        TGraphAsymmErrors* graphTheoryMCGillpPb         = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphDirectPhotonSpecMcGill5023GeV");

        TGraphAsymmErrors* graphTheoryPowhegDRnCTEQpPb  = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonPoweheg5023GeV_nCTEQ15_ALICECocktail");
        TGraph* graphTheoryPowhegDRnCTEQpPbCenter       = (TGraph*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonPoweheg5023GeV_nCTEQ15_ALICECocktail_Center");
        TGraphAsymmErrors* graphTheoryPowhegDREPPS16pPb = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonPoweheg5023GeV_EPPS16_ALICECocktail");
        TGraph* graphTheoryPowhegDREPPS16pPbCenter      = (TGraph*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonPoweheg5023GeV_EPPS16_ALICECocktail_Center");

        TGraphAsymmErrors* graphTheoryPowhegnCTEQpPb    = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphPowhegDirectPhotonInvYieldINT7_pPb5020_nCTEQ15");
        TGraphAsymmErrors* graphTheoryPowhegEPPS16pPb   = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphPowhegDirectPhotonInvYieldINT7_pPb5020_EPPS16");

        TGraph* graphTheoryDRAlwinaPeTer                = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonPeTerAlwina");
        TGraph* graphTheoryDRAlwinaPWGGA                = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonPWGGAAlwina");
        TGraph* graphTheoryDRAlwinaVogelsang            = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonVogelsangAlwina");
        TGraph* graphTheoryDRAlwinaJetPHOX              = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonJETPHOXAlwina");

    TFile* fileIsolatedGamma                       = new TFile(fileNameIsolatedGamma);
        TGraphErrors* graphIsolatedGammaStat            = (TGraphErrors*) fileIsolatedGamma->Get("invariant_yield_with_stat");
        TGraphErrors* graphIsolatedGammaSys             = (TGraphErrors*) fileIsolatedGamma->Get("invariant_yield_with_syst");

    //*******************************************************************************************************************************************
    //*********************************************** Combining Rgamma ratios  ******************************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsGamma[11]          = { 2,  8,  2,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
//     Int_t offSetsGammaSys[11]       = { 3,  8,  10,  0,  8,
//                                         0,  0,  0,  0,  0,
//                                         0};
    Int_t offSetsGammaSys[11]       = { 3,  8,  14,  0,  8,
                                        0,  0,  0,  0,  0,
                                        0};

    TGraphAsymmErrors* statErrorRelCollectionDR[11];
    TGraphAsymmErrors* statErrorRelCollectionDRNonFit[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionDR[i]        = NULL;
        statErrorRelCollectionDRNonFit[i]  = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (histoDRPi0FitStatErr[i]){
            statErrorRelCollectionDR[i]    = new TGraphAsymmErrors(histoDRPi0FitStatErr[i]);
            while (statErrorRelCollectionDR[i]->GetY()[0] == 0) statErrorRelCollectionDR[i]->RemovePoint(0);
            while (statErrorRelCollectionDR[i]->GetY()[statErrorRelCollectionDR[i]->GetN()-1] == 0) statErrorRelCollectionDR[i]->RemovePoint(statErrorRelCollectionDR[i]->GetN()-1);
            statErrorRelCollectionDR[i]    = CalculateRelErrUpAsymmGraph( statErrorRelCollectionDR[i], Form("relativeStatErrorDR_%s", nameMeasGlobal[i].Data()));
        }
        if (histoDRNonFitStatErr[i]){
            statErrorRelCollectionDRNonFit[i]   = new TGraphAsymmErrors(histoDRNonFitStatErr[i]);
            while (statErrorRelCollectionDRNonFit[i]->GetY()[0] == 0) statErrorRelCollectionDRNonFit[i]->RemovePoint(0);
            while (statErrorRelCollectionDRNonFit[i]->GetY()[statErrorRelCollectionDRNonFit[i]->GetN()-1] == 0)
                statErrorRelCollectionDRNonFit[i]->RemovePoint(statErrorRelCollectionDRNonFit[i]->GetN()-1);
            statErrorRelCollectionDRNonFit[i]    = CalculateRelErrUpAsymmGraph( statErrorRelCollectionDRNonFit[i], Form("relativeStatErrorDRNonFit_%s", nameMeasGlobal[i].Data()));
        }
    }

    TGraphAsymmErrors* sysErrorRelCollectionDR[11];
    TGraphAsymmErrors* sysErrorRelCollectionDRNonFit[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollectionDR[i]         = NULL;
        sysErrorRelCollectionDRNonFit[i]   = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        cout << i << endl;
        if (graphDRPi0FitSysErr[i]){
            sysErrorRelCollectionDR[i]     = (TGraphAsymmErrors*)graphDRPi0FitSysErr[i]->Clone(Form("relativeSysErrorDR_%s", nameMeasGlobal[i].Data()));
            sysErrorRelCollectionDR[i]->Print();
            while (sysErrorRelCollectionDR[i]->GetY()[0] == 0) sysErrorRelCollectionDR[i]->RemovePoint(0);
            while (sysErrorRelCollectionDR[i]->GetY()[sysErrorRelCollectionDR[i]->GetN()-1] == 0) sysErrorRelCollectionDR[i]->RemovePoint(sysErrorRelCollectionDR[i]->GetN()-1);
            sysErrorRelCollectionDR[i]     = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionDR[i], Form("relativeSysErrorDR_%s", nameMeasGlobal[i].Data()));
            cout << "after" << endl;
            sysErrorRelCollectionDR[i]->Print();
        }
        if (graphDRNonFitSysErr[i]){
            sysErrorRelCollectionDRNonFit[i]     = (TGraphAsymmErrors*)graphDRNonFitSysErr[i]->Clone(Form("relativeSysErrorDRNonFit_%s", nameMeasGlobal[i].Data()));
            sysErrorRelCollectionDRNonFit[i]->Print();
            while (sysErrorRelCollectionDRNonFit[i]->GetY()[0] == 0) sysErrorRelCollectionDRNonFit[i]->RemovePoint(0);
            while (sysErrorRelCollectionDRNonFit[i]->GetY()[sysErrorRelCollectionDRNonFit[i]->GetN()-1] == 0) sysErrorRelCollectionDRNonFit[i]->RemovePoint(sysErrorRelCollectionDRNonFit[i]->GetN()-1);
            sysErrorRelCollectionDRNonFit[i]     = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionDRNonFit[i], Form("relativeSysErrorDRNonFit_%s", nameMeasGlobal[i].Data()));
            cout << "after" << endl;
            sysErrorRelCollectionDRNonFit[i]->Print();
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraph* graphWeightsDR[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsDR[i]                   = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNameDROutputWeighting       = Form("%s/DR_WeightingMethod.dat",outputDirFit.Data());
    TGraphAsymmErrors* graphCombDRStat      = NULL;
    TGraphAsymmErrors* graphCombDRSys       = NULL;
    TGraphAsymmErrors* graphCombDRTot       = CombinePtPointsSpectraFullCorrMat(    histoDRPi0FitStatErr,    graphDRPi0FitSysErr,
                                                                                    xPtLimitsGamma, maxNBinsGamma,
                                                                                    offSetsGamma, offSetsGammaSys,
                                                                                    graphCombDRStat, graphCombDRSys,
                                                                                    fileNameDROutputWeighting, "pPb_5.023TeV", "RGamma", kTRUE,
                                                                                    NULL, fileNameCorrelations );


    if (graphCombDRTot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    while (graphCombDRStat->GetX()[0] < 0.3){
        graphCombDRStat->RemovePoint(0);
    }
    while (graphCombDRTot->GetX()[0] < 0.3){
        graphCombDRTot->RemovePoint(0);
    }
    while (graphCombDRSys->GetX()[0] < 0.3){
        graphCombDRSys->RemovePoint(0);
    }
    graphCombDRTot->Print();

    // Reading weights from output file for plotting
    ifstream fileWeightsDRRead;
    fileWeightsDRRead.open(fileNameDROutputWeighting,ios_base::in);
    cout << "reading" << fileNameDROutputWeighting << endl;
    Double_t xValuesDRRead[50];
    Double_t weightsDRRead[11][50];
    Int_t availableDRMeas[11]      = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetDR               = 4;
    Int_t nPtBinsDRRead            = 0;
    while(!fileWeightsDRRead.eof() && nPtBinsDRRead < 50){
        TString garbage             = "";
        if (nPtBinsDRRead == 0){
            fileWeightsDRRead >> garbage ;//>> availableDRMeas[0] >> availableDRMeas[1] >> availableDRMeas[2] >> availableDRMeas[3];
            for (Int_t i = 0; i < nMeasSetDR; i++){
                fileWeightsDRRead >> availableDRMeas[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableDRMeas[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsDRRead >> xValuesDRRead[nPtBinsDRRead-1];
            for (Int_t i = 0; i < nMeasSetDR; i++){
                fileWeightsDRRead >> weightsDRRead[availableDRMeas[i]][nPtBinsDRRead-1] ;
            }
            cout << "read: "<<  nPtBinsDRRead << "\t"<< xValuesDRRead[nPtBinsDRRead-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetDR; i++){
                cout << weightsDRRead[availableDRMeas[i]][nPtBinsDRRead-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsDRRead++;
    }
    nPtBinsDRRead                  = nPtBinsDRRead-2 ;
    fileWeightsDRRead.close();

    for (Int_t i = 0; i < nMeasSetDR; i++){
        graphWeightsDR[availableDRMeas[i]]                        = new TGraph(nPtBinsDRRead,xValuesDRRead,weightsDRRead[availableDRMeas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsDRRead; n++){
            if (graphWeightsDR[availableDRMeas[i]]->GetY()[bin] == 0) graphWeightsDR[availableDRMeas[i]]->RemovePoint(bin);
            else bin++;
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plotting weights method only EMC *****************************************
    // **********************************************************************************************************************
    Int_t textSizeLabelsPixel                 = 900*0.04;

    TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.01, 0.01, 0.09);
    canvasWeights->SetLogx();

    TH1F * histo2DDRWeights;
    histo2DDRWeights = new TH1F("histo2DDRWeights","histo2DDRWeights",11000,doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(histo2DDRWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DDRWeights->GetYaxis()->SetRangeUser(-0.7,1.3);
    histo2DDRWeights->GetXaxis()->SetMoreLogLabels();
    histo2DDRWeights->GetXaxis()->SetLabelOffset(-0.01);
    histo2DDRWeights->Draw("copy");

    histo2DDRWeights->Draw("copy");
    TLegend* legendWeightsDR   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetDR+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDR; i++){
        DrawGammaSetMarkerTGraph(graphWeightsDR[availableDRMeas[i]], markerStyleDet[availableDRMeas[i]], markerSizeDet[availableDRMeas[i]], colorDet[availableDRMeas[i]] , colorDet[availableDRMeas[i]]);
        graphWeightsDR[availableDRMeas[i]]->Draw("p,same,z");
        legendWeightsDR->AddEntry(graphWeightsDR[availableDRMeas[i]],nameMeasGlobalLabel[availableDRMeas[i]],"p");
    }
    legendWeightsDR->Draw();

    TLatex *labelWeightsEnergy      = new TLatex(0.95,0.20,collisionSystempPbNSD.Data());
    SetStyleTLatex( labelWeightsEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelWeightsEnergy->Draw();
    TLatex *labelWeightsDR         = new TLatex(0.95,0.15,"R_{#gamma}");
    SetStyleTLatex( labelWeightsDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelWeightsDR->Draw();

    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/DR_Weights.%s",outputDirFit.Data(),suffix.Data()));
    canvasWeights->SaveAs(Form("%s/DR_Weights.pdf",outputDirFit.Data()));
    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();

    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0,50.0);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,24.5);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelSysErr->Draw("copy");

    TLegend* legendRelSysErr2       = GetAndSetLegend2(0.62, 0.92-(0.04*(nMeasSetDR+1)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDR; i++){
        cout << "sys\t" << nameMeasGlobalLabel[availableDRMeas[i]] << endl;
        DrawGammaSetMarkerTGraph(sysErrorRelCollectionDR[availableDRMeas[i]], markerStyleDet[availableDRMeas[i]], markerSizeDet[availableDRMeas[i]], colorDet[availableDRMeas[i]],
                                 colorDet[availableDRMeas[i]]);
        sysErrorRelCollectionDR[availableDRMeas[i]]->Draw("p,same,z");
        sysErrorRelCollectionDR[availableDRMeas[i]]->Print();
        legendRelSysErr2->AddEntry(sysErrorRelCollectionDR[availableDRMeas[i]],nameMeasGlobalLabel[availableDRMeas[i]],"p");
    }
    legendRelSysErr2->Draw();

    TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystempPbNSD.Data());
    SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrEnergy->Draw();
    TLatex *labelRelSysErrDR       = new TLatex(0.15,0.85,"R_{#gamma}");
    SetStyleTLatex( labelRelSysErrDR, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DR_RelSysErr.%s",outputDirFit.Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/DR_RelSysErr.pdf",outputDirFit.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();

    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0,50.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetYaxis()->SetRangeUser(0,40.5);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelStatErr->Draw("copy");
    TLegend* legendRelStatErr2       = GetAndSetLegend2(0.14, 0.92-(0.04*(nMeasSetDR+1)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDR; i++){
        DrawGammaSetMarkerTGraph(statErrorRelCollectionDR[availableDRMeas[i]], markerStyleDet[availableDRMeas[i]], markerSizeDet[availableDRMeas[i]], colorDet[availableDRMeas[i]],
                                 colorDet[availableDRMeas[i]]);
        statErrorRelCollectionDR[availableDRMeas[i]]->Draw("p,same,z");
        legendRelStatErr2->AddEntry(statErrorRelCollectionDR[availableDRMeas[i]],nameMeasGlobalLabel[availableDRMeas[i]],"p");
    }
    legendRelStatErr2->Draw();

    TLatex *labelRelStatErrEnergy   = new TLatex(0.95,0.89,collisionSystempPbNSD.Data());
    SetStyleTLatex( labelRelStatErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrEnergy->Draw();
    TLatex *labelRelStatErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
    SetStyleTLatex( labelRelStatErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrDR->Draw();

    canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr.%s",outputDirFit.Data(),suffix.Data()));
    canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr.pdf",outputDirFit.Data()));



    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TGraphAsymmErrors* graphCombDRRelStat       = CalculateRelErrUpAsymmGraph( graphCombDRStat, "relativeStatErrorDR");
    TGraphAsymmErrors* graphCombDRRelSys        = CalculateRelErrUpAsymmGraph( graphCombDRSys, "relativeSysErrorDR");
    TGraphAsymmErrors* graphCombDRRelTot        = CalculateRelErrUpAsymmGraph( graphCombDRTot, "relativeTotalErrorDR");

    canvasRelSysErr->cd();
        TH2F * histo2DRelErr;
        histo2DRelErr                    = new TH2F("histo2DRelErr","histo2DRelErr",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0,50.0);
        SetStyleHistoTH2ForGraphs(histo2DRelErr, "#it{p}_{T} (GeV/#it{c})","Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelErr->GetYaxis()->SetRangeUser(0,40.5);
        histo2DRelErr->GetXaxis()->SetMoreLogLabels();
        histo2DRelErr->GetXaxis()->SetLabelOffset(-0.01);
        histo2DRelErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombDRRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombDRRelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombDRRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombDRRelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombDRRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombDRRelSys->SetLineStyle(7);
        graphCombDRRelSys->Draw("l,x0,same,e1");

        TLegend* legendRelTotErr       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
        legendRelTotErr->AddEntry(graphCombDRRelTot,"tot","p");
        legendRelTotErr->AddEntry(graphCombDRRelStat,"stat","l");
        legendRelTotErr->AddEntry(graphCombDRRelSys,"sys","l");
        legendRelTotErr->Draw();

        TLatex *labelRelTotErrEnergy   = new TLatex(0.95,0.89,collisionSystempPbNSD.Data());
        SetStyleTLatex( labelRelTotErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
        SetStyleTLatex( labelRelTotErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DR_Reldecomp.%s",outputDirFit.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraph* graphWeightsDRNonFit[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsDRNonFit[i]             = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNameDRNonFitOutputWeighting = Form("%s/DR_WeightingMethodNonFit.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombDRNonFitStat= NULL;
    TGraphAsymmErrors* graphCombDRNonFitSys = NULL;
    TGraphAsymmErrors* graphCombDRNonFitTot = CombinePtPointsSpectraFullCorrMat(    histoDRNonFitStatErr,    graphDRNonFitSysErr,
                                                                                    xPtLimitsGamma, maxNBinsGamma,
                                                                                    offSetsGamma, offSetsGammaSys,
                                                                                    graphCombDRNonFitStat, graphCombDRNonFitSys,
                                                                                    fileNameDRNonFitOutputWeighting, "pPb_5.023TeV", "RGamma", kTRUE,
                                                                                    NULL, fileNameCorrelations );


    if (graphCombDRNonFitTot == NULL || graphCombDRNonFitStat== NULL || graphCombDRNonFitSys == NULL ) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        cout << graphCombDRNonFitTot << "\t" << graphCombDRNonFitSys << "\t" << graphCombDRNonFitStat << endl;
        graphCombDRNonFitTot->Print();
        return;
    }
    while (graphCombDRNonFitStat->GetX()[0] < 0.3){
        graphCombDRNonFitStat->RemovePoint(0);
    }
    while (graphCombDRNonFitSys->GetX()[0] < 0.3){
        graphCombDRNonFitSys->RemovePoint(0);
    }
    while (graphCombDRNonFitTot->GetX()[0] < 0.3){
        graphCombDRNonFitTot->RemovePoint(0);
    }
    graphCombDRNonFitTot->Print();

    // Reading weights from output file for plotting
    ifstream fileWeightsDRNonFitRead;
    fileWeightsDRNonFitRead.open(fileNameDRNonFitOutputWeighting,ios_base::in);
    cout << "reading" << fileNameDRNonFitOutputWeighting << endl;
    Double_t xValuesDRNonFitRead[50];
    Double_t weightsDRNonFitRead[11][50];
    Int_t availableDRNonFitMeas[11]      = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetDRNonFit               = 4;
    Int_t nPtBinsDRNonFitRead            = 0;
    while(!fileWeightsDRNonFitRead.eof() && nPtBinsDRNonFitRead < 50){
        TString garbage             = "";
        if (nPtBinsDRNonFitRead == 0){
            fileWeightsDRNonFitRead >> garbage ;//>> availableDRNonFitMeas[0] >> availableDRNonFitMeas[1] >> availableDRNonFitMeas[2] >> availableDRNonFitMeas[3];
            for (Int_t i = 0; i < nMeasSetDRNonFit; i++){
                fileWeightsDRNonFitRead >> availableDRNonFitMeas[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableDRNonFitMeas[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsDRNonFitRead >> xValuesDRNonFitRead[nPtBinsDRNonFitRead-1];
            for (Int_t i = 0; i < nMeasSetDRNonFit; i++){
                fileWeightsDRNonFitRead >> weightsDRNonFitRead[availableDRNonFitMeas[i]][nPtBinsDRNonFitRead-1] ;
            }
            cout << "read: "<<  nPtBinsDRNonFitRead << "\t"<< xValuesDRNonFitRead[nPtBinsDRNonFitRead-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetDRNonFit; i++){
                cout << weightsDRNonFitRead[availableDRNonFitMeas[i]][nPtBinsDRNonFitRead-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsDRNonFitRead++;
    }
    nPtBinsDRNonFitRead                  = nPtBinsDRNonFitRead-2 ;
    fileWeightsDRNonFitRead.close();

    for (Int_t i = 0; i < nMeasSetDRNonFit; i++){
        graphWeightsDRNonFit[availableDRNonFitMeas[i]]                        = new TGraph(nPtBinsDRNonFitRead,xValuesDRNonFitRead,weightsDRNonFitRead[availableDRNonFitMeas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsDRNonFitRead; n++){
            if (graphWeightsDRNonFit[availableDRNonFitMeas[i]]->GetY()[bin] == 0) graphWeightsDRNonFit[availableDRNonFitMeas[i]]->RemovePoint(bin);
            else bin++;
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plotting weights method only EMC *****************************************
    // **********************************************************************************************************************

    canvasWeights->cd();
    histo2DDRWeights->Draw("copy");

    histo2DDRWeights->Draw("copy");
    TLegend* legendWeightsDRNonFit   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetDRNonFit+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDRNonFit; i++){
        DrawGammaSetMarkerTGraph(graphWeightsDRNonFit[availableDRNonFitMeas[i]], markerStyleDet[availableDRNonFitMeas[i]], markerSizeDet[availableDRNonFitMeas[i]], colorDet[availableDRNonFitMeas[i]] , colorDet[availableDRNonFitMeas[i]]);
        graphWeightsDRNonFit[availableDRNonFitMeas[i]]->Draw("p,same,z");
        legendWeightsDRNonFit->AddEntry(graphWeightsDRNonFit[availableDRNonFitMeas[i]],nameMeasGlobalLabel[availableDRNonFitMeas[i]],"p");
    }
    legendWeightsDRNonFit->Draw();

    labelWeightsEnergy->Draw();
    labelWeightsDR->Draw();

    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/DRNonFit_Weights.%s",outputDir.Data(),suffix.Data()));
    canvasWeights->SaveAs(Form("%s/DRNonFit_Weights.pdf",outputDir.Data()));
    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelSysErr->cd();
    histo2DRelSysErr->Draw("copy");

    for (Int_t i = 0; i < nMeasSetDRNonFit; i++){
        cout << "sys\t" << nameMeasGlobalLabel[availableDRNonFitMeas[i]] << endl;
        DrawGammaSetMarkerTGraph(sysErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]], markerStyleDet[availableDRNonFitMeas[i]], markerSizeDet[availableDRNonFitMeas[i]], colorDet[availableDRNonFitMeas[i]],
                                 colorDet[availableDRNonFitMeas[i]]);
        sysErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]]->Draw("p,same,z");
        sysErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]]->Print();
    }
    legendRelSysErr2->Draw();

    labelRelSysErrEnergy->Draw();
    labelRelSysErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DRNonFit_RelSysErr.%s",outputDir.Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/DRNonFit_RelSysErr.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    canvasRelStatErr->cd();
    histo2DRelStatErr->Draw("copy");

    for (Int_t i = 0; i < nMeasSetDRNonFit; i++){
        cout << "stat\t" << nameMeasGlobalLabel[availableDRNonFitMeas[i]] << endl;
        DrawGammaSetMarkerTGraph(statErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]], markerStyleDet[availableDRNonFitMeas[i]], markerSizeDet[availableDRNonFitMeas[i]],
                                 colorDet[availableDRNonFitMeas[i]],colorDet[availableDRNonFitMeas[i]]);
        statErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]]->Draw("p,same,z");
    }
    legendRelStatErr2->Draw();
    labelRelStatErrEnergy->Draw();
    labelRelStatErrDR->Draw();

    canvasRelStatErr->SaveAs(Form("%s/DRNonFit_RelStatErr.%s",outputDir.Data(),suffix.Data()));
    canvasRelStatErr->SaveAs(Form("%s/DRNonFit_RelStatErr.pdf",outputDir.Data()));


    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TGraphAsymmErrors* graphCombDRNonFitRelStat       = CalculateRelErrUpAsymmGraph( graphCombDRNonFitStat, "relativeStatErrorDRNonFit");
    TGraphAsymmErrors* graphCombDRNonFitRelSys        = CalculateRelErrUpAsymmGraph( graphCombDRNonFitSys, "relativeSysErrorDRNonFit");
    TGraphAsymmErrors* graphCombDRNonFitRelTot        = CalculateRelErrUpAsymmGraph( graphCombDRNonFitTot, "relativeTotalErrorDRNonFit");

    canvasRelSysErr->cd();
        histo2DRelErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombDRNonFitRelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombDRNonFitRelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombDRNonFitRelSys->SetLineStyle(7);
        graphCombDRNonFitRelSys->Draw("l,x0,same,e1");

        legendRelTotErr->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DRNonFit_Reldecomp.%s",outputDir.Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/DRNonFit_Reldecomp.pdf",outputDir.Data()));


    //*******************************************************************************************************************************************
    //*********************************************** Combining Rgamma ratios  ******************************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsIncGamma[11]       = { 2,  8,  2,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
//     Int_t offSetsIncGammaSys[11]    = { 3,  8,  10,  0,  8,
//                                         0,  0,  0,  0,  0,
//                                         0};
    Int_t offSetsIncGammaSys[11]    = { 3,  8,  14,  0,  8,
                                        0,  0,  0,  0,  0,
                                        0};

    TGraphAsymmErrors* statErrorGraphCollectionIncGamma[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorGraphCollectionIncGamma[i]   = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (histoIncGammaStatErr[i]){
            statErrorGraphCollectionIncGamma[i]    = new TGraphAsymmErrors(histoIncGammaStatErr[i]);
            while (statErrorGraphCollectionIncGamma[i]->GetY()[0] == 0) statErrorGraphCollectionIncGamma[i]->RemovePoint(0);
            while (statErrorGraphCollectionIncGamma[i]->GetY()[statErrorGraphCollectionIncGamma[i]->GetN()-1] == 0) statErrorGraphCollectionIncGamma[i]->RemovePoint(statErrorGraphCollectionIncGamma[i]->GetN()-1);
            statErrorGraphCollectionIncGamma[i]->SetName(Form("statErrorIncGamma_%s", nameMeasGlobal[i].Data()));
        }
    }


    TGraphAsymmErrors* statErrorRelCollectionIncGamma[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionIncGamma[i]        = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (histoIncGammaStatErr[i]){
            statErrorRelCollectionIncGamma[i]    = new TGraphAsymmErrors(histoIncGammaStatErr[i]);
            while (statErrorRelCollectionIncGamma[i]->GetY()[0] == 0) statErrorRelCollectionIncGamma[i]->RemovePoint(0);
            while (statErrorRelCollectionIncGamma[i]->GetY()[statErrorRelCollectionIncGamma[i]->GetN()-1] == 0) statErrorRelCollectionIncGamma[i]->RemovePoint(statErrorRelCollectionIncGamma[i]->GetN()-1);
            statErrorRelCollectionIncGamma[i]    = CalculateRelErrUpAsymmGraph( statErrorRelCollectionIncGamma[i], Form("relativeStatErrorIncGamma_%s", nameMeasGlobal[i].Data()));
        }
    }


    TGraphAsymmErrors* sysErrorRelCollectionIncGamma[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollectionIncGamma[i]         = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        cout << i << endl;
        if (graphIncGammaSysErr[i]){
            sysErrorRelCollectionIncGamma[i]     = (TGraphAsymmErrors*)graphIncGammaSysErr[i]->Clone(Form("relativeSysErrorIncGamma_%s", nameMeasGlobal[i].Data()));
            sysErrorRelCollectionIncGamma[i]->Print();
            while (sysErrorRelCollectionIncGamma[i]->GetY()[0] == 0) sysErrorRelCollectionIncGamma[i]->RemovePoint(0);
            while (sysErrorRelCollectionIncGamma[i]->GetY()[sysErrorRelCollectionIncGamma[i]->GetN()-1] == 0) sysErrorRelCollectionIncGamma[i]->RemovePoint(sysErrorRelCollectionIncGamma[i]->GetN()-1);
            sysErrorRelCollectionIncGamma[i]     = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionIncGamma[i], Form("relativeSysErrorIncGamma_%s", nameMeasGlobal[i].Data()));
            cout << "after" << endl;
            sysErrorRelCollectionIncGamma[i]->Print();
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraph* graphWeightsIncGamma[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsIncGamma[i]                   = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNameIncGammaOutputWeighting       = Form("%s/IncGamma_WeightingMethod.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombIncGammaStat      = NULL;
    TGraphAsymmErrors* graphCombIncGammaSys       = NULL;
    TGraphAsymmErrors* graphCombIncGammaTot       = CombinePtPointsSpectraFullCorrMat(      histoIncGammaStatErr,    graphIncGammaSysErr,
                                                                                            xPtLimitsGamma, maxNBinsGamma,
                                                                                            offSetsIncGamma, offSetsIncGammaSys,
                                                                                            graphCombIncGammaStat, graphCombIncGammaSys,
                                                                                            fileNameIncGammaOutputWeighting, "pPb_5.023TeV", "GammaInc", kTRUE,
                                                                                            NULL, fileNameCorrelations );


    if (graphCombIncGammaTot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    while (graphCombIncGammaStat->GetX()[0] < 0.3){
        graphCombIncGammaStat->RemovePoint(0);
    }
    while (graphCombIncGammaTot->GetX()[0] < 0.3){
        graphCombIncGammaTot->RemovePoint(0);
    }
    while (graphCombIncGammaSys->GetX()[0] < 0.3){
        graphCombIncGammaSys->RemovePoint(0);
    }
    graphCombIncGammaTot->Print();

    // Reading weights from output file for plotting
    ifstream fileWeightsIncGammaRead;
    fileWeightsIncGammaRead.open(fileNameIncGammaOutputWeighting,ios_base::in);
    cout << "reading" << fileNameIncGammaOutputWeighting << endl;
    Double_t xValuesIncGammaRead[50];
    Double_t weightsIncGammaRead[11][50];
    Int_t availableIncGammaMeas[11]      = { -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1,
                                            -1};
    Int_t nMeasSetIncGamma               = 4;
    Int_t nPtBinsIncGammaRead            = 0;
    while(!fileWeightsIncGammaRead.eof() && nPtBinsIncGammaRead < 50){
        TString garbage             = "";
        if (nPtBinsIncGammaRead == 0){
            fileWeightsIncGammaRead >> garbage ;//>> availableIncGammaMeas[0] >> availableIncGammaMeas[1] >> availableIncGammaMeas[2] >> availableIncGammaMeas[3];
            for (Int_t i = 0; i < nMeasSetIncGamma; i++){
                fileWeightsIncGammaRead >> availableIncGammaMeas[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableIncGammaMeas[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsIncGammaRead >> xValuesIncGammaRead[nPtBinsIncGammaRead-1];
            for (Int_t i = 0; i < nMeasSetIncGamma; i++){
                fileWeightsIncGammaRead >> weightsIncGammaRead[availableIncGammaMeas[i]][nPtBinsIncGammaRead-1] ;
            }
            cout << "read: "<<  nPtBinsIncGammaRead << "\t"<< xValuesIncGammaRead[nPtBinsIncGammaRead-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetIncGamma; i++){
                cout << weightsIncGammaRead[availableIncGammaMeas[i]][nPtBinsIncGammaRead-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsIncGammaRead++;
    }
    nPtBinsIncGammaRead                  = nPtBinsIncGammaRead-2 ;
    fileWeightsIncGammaRead.close();

    for (Int_t i = 0; i < nMeasSetIncGamma; i++){
        graphWeightsIncGamma[availableIncGammaMeas[i]]                        = new TGraph(nPtBinsIncGammaRead,xValuesIncGammaRead,weightsIncGammaRead[availableIncGammaMeas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsIncGammaRead; n++){
            if (graphWeightsIncGamma[availableIncGammaMeas[i]]->GetY()[bin] == 0) graphWeightsIncGamma[availableIncGammaMeas[i]]->RemovePoint(bin);
            else bin++;
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plotting weights method only EMC *****************************************
    // **********************************************************************************************************************
    canvasWeights->cd();
    TH1F * histo2DIncGammaWeights;
    histo2DIncGammaWeights = new TH1F("histo2DIncGammaWeights","histo2DIncGammaWeights",11000,doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(histo2DIncGammaWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DIncGammaWeights->GetYaxis()->SetRangeUser(-3.7,4.3);
    histo2DIncGammaWeights->GetXaxis()->SetMoreLogLabels();
    histo2DIncGammaWeights->GetXaxis()->SetLabelOffset(-0.01);
    histo2DIncGammaWeights->Draw("copy");

    TLegend* legendWeightsIncGamma   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetIncGamma+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetIncGamma; i++){
        DrawGammaSetMarkerTGraph(graphWeightsIncGamma[availableIncGammaMeas[i]], markerStyleDet[availableIncGammaMeas[i]], markerSizeDet[availableIncGammaMeas[i]], colorDet[availableIncGammaMeas[i]] , colorDet[availableIncGammaMeas[i]]);
        graphWeightsIncGamma[availableIncGammaMeas[i]]->Draw("p,same,z");
        legendWeightsIncGamma->AddEntry(graphWeightsIncGamma[availableIncGammaMeas[i]],nameMeasGlobalLabel[availableIncGammaMeas[i]],"p");
    }
    legendWeightsIncGamma->Draw();

    labelWeightsEnergy->Draw();
    TLatex *labelWeightsIncGamma         = new TLatex(0.95,0.15,"#gamma_{inc}");
    SetStyleTLatex( labelWeightsIncGamma, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelWeightsIncGamma->Draw();

    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/IncGamma_Weights.%s",outputDir.Data(),suffix.Data()));
    canvasWeights->SaveAs(Form("%s/IncGamma_Weights.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelSysErr->cd();

    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,17.5);
    histo2DRelSysErr->Draw("copy");

    TLegend* legendRelSysErrIncGamma       = GetAndSetLegend2(0.62, 0.92-(0.04*(nMeasSetIncGamma+1)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetIncGamma; i++){
        cout << "sys\t" << nameMeasGlobalLabel[availableIncGammaMeas[i]] << endl;
        DrawGammaSetMarkerTGraph(sysErrorRelCollectionIncGamma[availableIncGammaMeas[i]], markerStyleDet[availableIncGammaMeas[i]], markerSizeDet[availableIncGammaMeas[i]], colorDet[availableIncGammaMeas[i]],
                                 colorDet[availableIncGammaMeas[i]]);
        sysErrorRelCollectionIncGamma[availableIncGammaMeas[i]]->Draw("p,same,z");
        sysErrorRelCollectionIncGamma[availableIncGammaMeas[i]]->Print();
        legendRelSysErrIncGamma->AddEntry(sysErrorRelCollectionIncGamma[availableIncGammaMeas[i]],nameMeasGlobalLabel[availableIncGammaMeas[i]],"p");
    }
    legendRelSysErrIncGamma->Draw();

    labelRelSysErrEnergy->Draw();
    TLatex *labelRelSysErrIncGamma       = new TLatex(0.15,0.85,"#gamma_{inc}");
    SetStyleTLatex( labelRelSysErrIncGamma, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrIncGamma->Draw();

    canvasRelSysErr->SaveAs(Form("%s/IncGamma_RelSysErr.%s",outputDir.Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/IncGamma_RelSysErr.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    canvasRelStatErr->cd();

    histo2DRelStatErr->GetYaxis()->SetRangeUser(-0.2,30.0);
    histo2DRelStatErr->Draw("copy");
    TLegend* legendRelStatErrIncGamma       = GetAndSetLegend2(0.14, 0.92-(0.04*(nMeasSetIncGamma+1)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetIncGamma; i++){
        DrawGammaSetMarkerTGraph(statErrorRelCollectionIncGamma[availableIncGammaMeas[i]], markerStyleDet[availableIncGammaMeas[i]], markerSizeDet[availableIncGammaMeas[i]], colorDet[availableIncGammaMeas[i]],
                                 colorDet[availableIncGammaMeas[i]]);
        statErrorRelCollectionIncGamma[availableIncGammaMeas[i]]->Draw("p,same,z");
        legendRelStatErrIncGamma->AddEntry(statErrorRelCollectionIncGamma[availableIncGammaMeas[i]],nameMeasGlobalLabel[availableIncGammaMeas[i]],"p");
    }
    legendRelStatErrIncGamma->Draw();

    labelRelStatErrEnergy->Draw();
    TLatex *labelRelStatErrIncGamma      = new TLatex(0.95,0.85,"#gamma_{inc}");
    SetStyleTLatex( labelRelStatErrIncGamma, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrIncGamma->Draw();

    canvasRelStatErr->SaveAs(Form("%s/IncGamma_RelStatErr.%s",outputDir.Data(),suffix.Data()));
    canvasRelStatErr->SaveAs(Form("%s/IncGamma_RelStatErr.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TGraphAsymmErrors* graphCombIncGammaRelStat     = CalculateRelErrUpAsymmGraph( graphCombIncGammaStat, "relativeStatErrorDR");
    TGraphAsymmErrors* graphCombIncGammaRelSys      = CalculateRelErrUpAsymmGraph( graphCombIncGammaSys, "relativeSysErrorDR");
    TGraphAsymmErrors* graphCombIncGammaRelTot      = CalculateRelErrUpAsymmGraph( graphCombIncGammaTot, "relativeTotalErrorDR");

        canvasRelSysErr->cd();
        histo2DRelErr->GetYaxis()->SetRangeUser(-0.2,30.5);
        histo2DRelErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombIncGammaRelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombIncGammaRelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombIncGammaRelSys->SetLineStyle(7);
        graphCombIncGammaRelSys->Draw("l,x0,same,e1");

        legendRelTotErr->Draw();
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrIncGamma      = new TLatex(0.95,0.85,"#gamma_{inc}");
        SetStyleTLatex( labelRelTotErrIncGamma, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrIncGamma->Draw();

    canvasRelSysErr->SaveAs(Form("%s/IncGamma_Reldecomp.%s",outputDir.Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/IncGamma_Reldecomp.pdf",outputDir.Data()));

    //*******************************************************************************************************************************************
    //************************************************* Fitting gamma spectrum ******************************************************************
    //*******************************************************************************************************************************************

    Double_t paramGraphoHag[5]      = {82,-0.22,-0.01,0.57,6.28};
    TF1* fitHagGammaComb            = FitObject("oHag","fitHagGammaComb","Gamma",graphCombIncGammaTot,graphCombIncGammaTot->GetX()[2],
                                                graphCombIncGammaTot->GetX()[graphCombIncGammaTot->GetN()],paramGraphoHag,"QNRME+");
    Double_t paramTCM[5] = {graphCombIncGammaTot->GetY()[1],0.1,graphCombIncGammaTot->GetY()[4],0.6,3};
    TF1* fitTCMGammaComb            = FitObject("tcm", "fitTCMGammaComb","Gamma", graphCombIncGammaTot, graphCombIncGammaTot->GetX()[0],
                                                graphCombIncGammaTot->GetX()[graphCombIncGammaTot->GetN()], paramTCM,"QNRME+");
    TString forOutput               = WriteParameterToFile(fitTCMGammaComb);
    cout << forOutput.Data() << endl;

    TF1* fitTsallisGammaComb        = FitObject("l","fitTsallisGammaComb","Gamma",graphCombIncGammaTot,graphCombIncGammaTot->GetX()[0],
                                                graphCombIncGammaTot->GetX()[graphCombIncGammaTot->GetN()],NULL,"QNRME+");
    forOutput                       = WriteParameterToFile(fitTsallisGammaComb);
    cout << forOutput.Data() << endl;


    // **********************************************************************************************************************
    // ************************************* Calculating bin shifted spectra & fitting **************************************
    // **********************************************************************************************************************

    // Cloning spectra
    TGraphAsymmErrors* graphCombIncGammaTotUnshi         = (TGraphAsymmErrors*)graphCombIncGammaTot->Clone("GammaUnshifted");
    TGraphAsymmErrors* graphCombIncGammaStatUnshi        = (TGraphAsymmErrors*)graphCombIncGammaStat->Clone("GammaUnshiftedStat");
    TGraphAsymmErrors* graphCombIncGammaSysUnshi         = (TGraphAsymmErrors*)graphCombIncGammaSys->Clone("GammaUnshiftedSys");

    TGraphAsymmErrors* graphIndGammaIncStatUnshi[11]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndGammaIncSysUnshi[11]      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndGammaIncStat[11]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndGammaIncSys[11]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndGammaIncStat_yShifted[11] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndGammaIncSys_yShifted[11]  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    for (Int_t i = 0; i< 11; i++){
        if (statErrorGraphCollectionIncGamma[i]){
            graphIndGammaIncStatUnshi[i]                 = (TGraphAsymmErrors*)statErrorGraphCollectionIncGamma[i]->Clone(Form("GammaUnshiftedStat%s",nameMeasGlobalLabel[i].Data()));
            graphIndGammaIncStat[i]                      = (TGraphAsymmErrors*)statErrorGraphCollectionIncGamma[i]->Clone(Form("GammaStat%s",nameMeasGlobalLabel[i].Data()));
            graphIndGammaIncStat_yShifted[i]             = (TGraphAsymmErrors*)statErrorGraphCollectionIncGamma[i]->Clone(Form("GammaYShiftedStat%s",nameMeasGlobalLabel[i].Data()));
        }
        if (graphIncGammaSysErr[i]){
            graphIndGammaIncSysUnshi[i]                  = (TGraphAsymmErrors*)graphIncGammaSysErr[i]->Clone(Form("GammaUnshiftedSys%s",nameMeasGlobalLabel[i].Data()));
            graphIndGammaIncSys[i]                       = (TGraphAsymmErrors*)graphIncGammaSysErr[i]->Clone(Form("GammaSys%s",nameMeasGlobalLabel[i].Data()));
            graphIndGammaIncSys_yShifted[i]              = (TGraphAsymmErrors*)graphIncGammaSysErr[i]->Clone(Form("GammaYShiftedSys%s",nameMeasGlobalLabel[i].Data()));
        }
    }

    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************
    TF1* fitShiftingGamma            = FitObject("tmpt","ShiftingGamma","Gamma");
    fitShiftingGamma->SetParameters(fitTsallisGammaComb->GetParameter(0),fitTsallisGammaComb->GetParameter(1), fitTsallisGammaComb->GetParameter(2));

    TGraphAsymmErrors* graphCombIncGammaTotNoShift = (TGraphAsymmErrors*) graphCombIncGammaTot->Clone("Gamma_NoShift");

//     graphCombIncGammaTot            = ApplyXshift(graphCombIncGammaTot, fitShiftingGamma);
    cout << "comb" << endl;
//     graphCombIncGammaStat->Print();
//     graphCombIncGammaStat           = ApplyXshiftIndividualSpectra( graphCombIncGammaTot,
//                                                                     graphCombIncGammaStat,
//                                                                     fitShiftingGamma,
//                                                                     0, graphCombIncGammaStat->GetN());
//     graphCombIncGammaSys            = ApplyXshiftIndividualSpectra( graphCombIncGammaTot,
//                                                                     graphCombIncGammaSys,
//                                                                     fitShiftingGamma,
//                                                                     0, graphCombIncGammaSys->GetN());
    Int_t offSetGammaShifting[11]   = { 0,  0,  7,  0,  4,
                                        0,  0,  0,  0,  0,
                                        0 };
    Int_t nComBinsGammaShifting[11] = { 24, 0, 17, 0,  20,
                                        0,  0,  0,  0,  0,
                                        0 };

//     for (Int_t i = 0; i< 11; i++){
//         if (graphIndGammaIncStat[i]){
//             cout << "shiting stat err of " << nameMeasGlobalLabel[i].Data();
//             graphIndGammaIncStat[i]  = ApplyXshiftIndividualSpectra(    graphCombIncGammaTot,
//                                                                         graphIndGammaIncStat[i],
//                                                                         fitShiftingGamma,
//                                                                         offSetGammaShifting[i], nComBinsGammaShifting[i]);
//
//         }
//         if (graphIndGammaIncSys[i]){
//             cout << "shiting sys err of " << nameMeasGlobalLabel[i].Data();
//             graphIndGammaIncSys[i]   = ApplyXshiftIndividualSpectra(    graphCombIncGammaTot,
//                                                                         graphIndGammaIncSys[i],
//                                                                         fitShiftingGamma,
//                                                                         offSetGammaShifting[i], nComBinsGammaShifting[i]);
//         }
//     }

    TGraphAsymmErrors* graphRatioGammaIndCombFitStat[11] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioGammaIndCombFitSys[11]  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    TGraphAsymmErrors* graphRatioGammaCombCombFitTot     = (TGraphAsymmErrors*)graphCombIncGammaTot->Clone();
    graphRatioGammaCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitTot, fitTCMGammaComb);
    TGraphAsymmErrors* graphRatioGammaCombCombFitStat    = (TGraphAsymmErrors*)graphCombIncGammaStat->Clone();
    graphRatioGammaCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitStat, fitTCMGammaComb);
    TGraphAsymmErrors* graphRatioGammaCombCombFitSys     = (TGraphAsymmErrors*)graphCombIncGammaSys->Clone();
    graphRatioGammaCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitSys, fitTCMGammaComb);

    for (Int_t i= 0; i< 11; i++){
        if (graphIndGammaIncStat[i]){
            graphRatioGammaIndCombFitStat[i]             = (TGraphAsymmErrors*)graphIndGammaIncStat[i]->Clone(Form("RatioGamma%sToCombFitStat", nameMeasGlobalLabel[i].Data()));
            graphRatioGammaIndCombFitStat[i]             = CalculateGraphErrRatioToFit(graphRatioGammaIndCombFitStat[i], fitTCMGammaComb);
        }
        if (graphIndGammaIncSys[i]){
            graphRatioGammaIndCombFitSys[i]              = (TGraphAsymmErrors*)graphIndGammaIncSys[i]->Clone(Form("RatioGamma%sToCombFitSyst", nameMeasGlobalLabel[i].Data()));
            graphRatioGammaIndCombFitSys[i]              = CalculateGraphErrRatioToFit(graphRatioGammaIndCombFitSys[i], fitTCMGammaComb);
        }
    }

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
    histo2DGammaRatioToCombFit               = new TH2F("histo2DGammaRatioToCombFit","histo2DGammaRatioToCombFit",1000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DGammaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPPb, textsizeLabelsPPb,
                              0.85*textsizeLabelsPPb,textsizeLabelsPPb, 0.9, 0.65, 510, 505);
    histo2DGammaRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DGammaRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
    //  histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(0.75,1.33);
    histo2DGammaRatioToCombFit->Draw("copy");

    ProduceGraphAsymmWithoutXErrors(graphRatioGammaCombCombFitStat);

    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitSys, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE);
    graphRatioGammaCombCombFitSys->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitStat, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb);
    graphRatioGammaCombCombFitStat->Draw("p,same,z");

    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1.,0.1, kGray+2);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1.05, 1.05,0.1, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.95, 0.95,0.1, kGray, 7);

    TLatex *labelRatioToFitEnergy2      = new TLatex(0.95, 0.22, collisionSystempPbNSD.Data());
    SetStyleTLatex( labelRatioToFitEnergy2, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergy2->Draw();
    TLatex *labelRatioToFitALICE2       = new TLatex(0.95, 0.16, textALICE.Data());
    SetStyleTLatex( labelRatioToFitALICE2, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitALICE2->Draw();
    TLatex *labelRatioToFitGamma        = new TLatex(0.15, 0.92, "#gamma_{inc}");
    SetStyleTLatex( labelRatioToFitGamma, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitGamma->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_PPb5023GeV.%s",outputDir.Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_PPb5023GeV.pdf",outputDir.Data()));
    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DGammaRatioToCombFit->Draw("copy");

    for (Int_t i = 10; i > -1 ; i--){
        if (graphRatioGammaIndCombFitSys[i] && i!=4){ //
            DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
            graphRatioGammaIndCombFitSys[i]->Draw("E2same");
        }
        if (graphRatioGammaIndCombFitStat[i] && i!=4){ //
            ProduceGraphAsymmWithoutXErrors(graphRatioGammaIndCombFitStat[i]);
            DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitStat[i], markerStyleDet[i] ,markerSizeDet[i], colorDet[i], colorDet[i]);
            graphRatioGammaIndCombFitStat[i]->Draw("p,same,z");
        }
    }
    graphRatioGammaIndCombFitStat[4]->Draw("p,same,z");

    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1] , 1., 1.,0.5, kGray+2);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1] , 1.05, 1.05,0.5, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1] , 0.95, 0.95,0.5, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1] , 1.1, 1.1,0.5, kGray, 9);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1] , 0.9, 0.9,0.5, kGray, 9);

    TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, collisionSystempPbNSD.Data());
    SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergy->Draw();
    TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICEPerf.Data());
    SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitALICE->Draw();
    TLatex *labelRatioToFitGamma2        = new TLatex(0.95, 0.80, "#gamma_{inc}");
    SetStyleTLatex( labelRatioToFitGamma2, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitGamma2->Draw();
    histo2DGammaRatioToCombFit->Draw("same,axis");


    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    Double_t rowsLegendGammaRatio[4]          = {0.92, 0.86, 0.80, 0.74};
    Double_t rowsLegendGammaRatioAbs[4]       = {1.37, 1.26, 1.22, 1.18 };
    Double_t columnsLegendGammaRatio[6]       = {0.115, 0.232, 0.32, 0.48, 0.7, 0.8};
    Double_t columnsLegendGammaRatioAbs[6]    = {0.215, 0.69, 1.15, 2, 6.6, 9.6};
    Double_t lengthBox                        = 0.2;
    Double_t heightBox                        = 0.03/2;
    //****************** first Column **************************************************
    TLatex *textPCMRatioGamma                 = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[1],nameMeasGlobalLabel[0]);
    SetStyleTLatex( textPCMRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textPCMRatioGamma->Draw();
    TLatex *textEMCALRatioGamma               = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[2],nameMeasGlobalLabel[2]);
    SetStyleTLatex( textEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textEMCALRatioGamma->Draw();
//     TLatex *textPCMEMCALRatioGamma            = new TLatex(columnsLegendGammaRatio[3],rowsLegendGammaRatio[1],nameMeasGlobalLabel[4]);
//     SetStyleTLatex( textPCMEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
//     textPCMEMCALRatioGamma->Draw();
//     TLatex *textPHOSRatioGamma            = new TLatex(columnsLegendGammaRatio[3],rowsLegendGammaRatio[2],nameMeasGlobalLabel[1]);
    TLatex *textPHOSRatioGamma            = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[3],nameMeasGlobalLabel[1]);
    SetStyleTLatex( textPHOSRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textPHOSRatioGamma->Draw();

    //****************** second Column *************************************************
    TLatex *textStatRatioGamma                = new TLatex(columnsLegendGammaRatio[1],rowsLegendGammaRatio[0] ,"stat");
    SetStyleTLatex( textStatRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textStatRatioGamma->Draw();
    TLatex *textSysRatioGamma                 = new TLatex(columnsLegendGammaRatio[2] ,rowsLegendGammaRatio[0],"syst");
    SetStyleTLatex( textSysRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textSysRatioGamma->Draw();
//     TLatex *textStatRatioGamma2               = new TLatex(columnsLegendGammaRatio[4],rowsLegendGammaRatio[0] ,"stat");
//     SetStyleTLatex( textStatRatioGamma2, textSizeLabelsPixel,4, 1, 43);
//     textStatRatioGamma2->Draw();
//     TLatex *textSysRatioGamma2                = new TLatex(columnsLegendGammaRatio[5] ,rowsLegendGammaRatio[0],"syst");
//     SetStyleTLatex( textSysRatioGamma2, textSizeLabelsPixel,4, 1, 43);
//     textSysRatioGamma2->Draw();

    TMarker* markerPCMGammaRatio           = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[0],columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[1],1);
    markerPCMGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[1]);
    TMarker* markerEMCALGammaRatio         = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[2], columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[2],1);
    markerEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[2]);
//     TMarker* markerPCMEMCALGammaRatio      = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[4], columnsLegendGammaRatio[4] ,rowsLegendGammaRatio[1],1);
//     markerPCMEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[4] ,rowsLegendGammaRatioAbs[1]);
//     TMarker* markerPHOSammaRatio      = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[1], columnsLegendGammaRatio[4] ,rowsLegendGammaRatio[2],1);
//     markerPHOSammaRatio->DrawMarker(columnsLegendGammaRatioAbs[4] ,rowsLegendGammaRatioAbs[2]);
    TMarker* markerPHOSammaRatio      = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[1], columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[3],1);
    markerPHOSammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[3]);

    TBox* boxPCMGammaRatio                 = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[0], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[1]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[1]+ heightBox);
    boxPCMGammaRatio->Draw("l");
    TBox* boxEMCALGammaRatio               = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[2], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[2]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[2]+ heightBox);
    boxEMCALGammaRatio->Draw("l");
//     TBox* boxPCMEMCALGammaRatio            = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[4], columnsLegendGammaRatioAbs[5]-0.5*lengthBox , rowsLegendGammaRatioAbs[1]- heightBox, columnsLegendGammaRatioAbs[5]+ 18*lengthBox, rowsLegendGammaRatioAbs[1]+ heightBox);
//     boxPCMEMCALGammaRatio->Draw("l");
//     TBox* boxPHOSGammaRatio            = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[1], columnsLegendGammaRatioAbs[5]-0.5*lengthBox , rowsLegendGammaRatioAbs[2]- heightBox,
    TBox* boxPHOSGammaRatio            = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[1], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[3]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[3]+ heightBox);
    boxPHOSGammaRatio->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit.%s",outputDir.Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit.pdf",outputDir.Data()));

    //*******************************************************************************************************
    //************************** Calculating combined direct photon spectrum ********************************
    //*******************************************************************************************************
    Double_t xArrayCombDR[graphCombDRNonFitStat->GetN()+1];
    xArrayCombDR[0] = graphCombDRNonFitStat->GetX()[0] - graphCombDRNonFitStat->GetEXhigh()[0];
    for (Int_t i = 1; i<graphCombDRNonFitStat->GetN()+1;i++){
        xArrayCombDR[i] = graphCombDRNonFitStat->GetX()[i-1] + graphCombDRNonFitStat->GetEXhigh()[i-1];
    }
    //_______________________ copy inclusive photon spectra _____________________
    TH1D *histoCombDirGammaSpectrumErrSum                   = new TH1D("histoCombDirGammaSpectrumErrSum","",graphCombDRNonFitStat->GetN(),xArrayCombDR);
    TH1D *histoCombDirGammaSpectrumErrSys                   = new TH1D("histoCombDirGammaSpectrumErrSys","",graphCombDRNonFitStat->GetN(),xArrayCombDR);
    TH1D *histoCombDirGammaSpectrumErrStat                  = new TH1D("histoCombDirGammaSpectrumErrStat","",graphCombDRNonFitStat->GetN(),xArrayCombDR);

    //_______________________ get arrays of double ratio errors __________________
    Double_t *SystErrorsCombDR                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *sumErrorsCombDR                               = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *StatErrorsCombDR                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *xErrorsDR                                     = new Double_t[graphCombIncGammaStat->GetN()];
    for (Int_t i = 0; i< graphCombDRStat->GetN(); i++){
        SystErrorsCombDR[i]                                 = graphCombDRNonFitSys->GetEYhigh()[i]/graphCombDRNonFitSys->GetY()[i] *100;
        StatErrorsCombDR[i]                                 = graphCombDRNonFitStat->GetEYhigh()[i]/graphCombDRNonFitStat->GetY()[i] *100;
        sumErrorsCombDR[i]                                  = graphCombDRNonFitTot->GetEYhigh()[i]/graphCombDRNonFitTot->GetY()[i] *100;
        //cout << i << "\t" << graphCombDRSys->GetY()[i] << "\t" << graphCombDRSys->GetEYhigh()[i] << "\t" <<SystErrorsCombDR[i] << endl;
    }
    xErrorsDR                                               = graphCombDRStat->GetX();

    cout << __LINE__ << endl;
    //graphCombDRTot->Print();

    //_______________________ copy inclusive photon spectra _____________________
    TH1D* histoCombErrorsForDRSum                           = new TH1D("histoCombErrorsForDRSum","",graphCombDRNonFitStat->GetN(),xArrayCombDR);
    TH1D* histoCombErrorsForDRStat                          = new TH1D("histoCombErrorsForDRStat","",graphCombDRNonFitStat->GetN(),xArrayCombDR);
    TH1D* histoCombErrorsForDRSys                           = new TH1D("histoCombErrorsForDRSys","",graphCombDRNonFitStat->GetN(),xArrayCombDR);

    for(Int_t i = 1; i<graphCombDRStat->GetN()+1;i++){
        //cout<< i << "\t"<<xErrorsDR[i-1]<<"  "<<histoCombErrorsForDRSum->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum->GetBinWidth(i) <<endl;
        Double_t binErrorSummed                                 = sumErrorsCombDR[i-1];
        Double_t binErrorSyst                                   = SystErrorsCombDR[i-1];
        Double_t binErrorStat                                   = StatErrorsCombDR[i-1];
        Double_t DR                                             = graphCombDRNonFitStat->GetY()[i-1];

        //cout << DR << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
        histoCombErrorsForDRSum->SetBinContent(i,DR);
        histoCombErrorsForDRSys->SetBinContent(i,DR);
        histoCombErrorsForDRStat->SetBinContent(i,DR);
        histoCombErrorsForDRSum->SetBinError(i,(binErrorSummed/100)*DR);
        histoCombErrorsForDRSys->SetBinError(i,(binErrorSyst/100)*DR);
        histoCombErrorsForDRStat->SetBinError(i,(binErrorStat/100)*DR);
    }

    for(Int_t i = 1; i<histoCombErrorsForDRSum->GetNbinsX()+1;i++){
        histoCombDirGammaSpectrumErrSum->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSys->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrStat->SetBinContent(i+1,-1);

        histoCombDirGammaSpectrumErrSum->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSys->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrStat->SetBinError(i+1,0);
    }

    // get the binning of the direct photons from the DR
    TH1D *histoCombDirGammaSpecSysErr                        = new TH1D(*histoCombErrorsForDRSys);
    TH1D *histoCombDirGammaSpecStatErr                       = new TH1D(*histoCombErrorsForDRStat);
    TH1D *histoCombDirGammaSpecSumErr                        = new TH1D(*histoCombErrorsForDRSum);


    cout << __LINE__ << endl;
    //graphCombDRStat->Print();
    for(Int_t i = 1; i<graphCombDRStat->GetN()+1; i++){
        // obtain common quantities
        Double_t Rgamma                 = histoCombErrorsForDRSys->GetBinContent(i);
        Double_t nIncGamma              = graphCombIncGammaStat->GetY()[i-1];

        // calculating Systematics graph
        Double_t errRgamma              = histoCombErrorsForDRSys->GetBinError(i);
        Double_t errNIncGam             = graphCombIncGammaSys->GetEYhigh()[i-1];
        Double_t q1                     = 1 - 1/ Rgamma;

        Double_t q1Error                = errRgamma/(Rgamma*Rgamma);
        Double_t content                = nIncGamma * ( 1 - 1/ Rgamma);
        Double_t error                  = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        Double_t errDR                  = content - error;
        histoCombDirGammaSpecSysErr->SetBinError(i, error);
        histoCombDirGammaSpecSysErr->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSys->SetBinContent(i, errDR);

        // calculating Stat graphs
        errRgamma                       = histoCombErrorsForDRStat->GetBinError(i);
        errNIncGam                      = graphCombIncGammaStat->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecStatErr->SetBinError(i, error);
        histoCombDirGammaSpecStatErr->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrStat->SetBinContent(i, errDR);

        // calculating summed error graphs
        errRgamma                       = histoCombErrorsForDRSum->GetBinError(i);
        errNIncGam                      = graphCombIncGammaTot->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecSumErr->SetBinError(i, error);
        histoCombDirGammaSpecSumErr->SetBinContent(i, content);
        histoCombDirGammaSpectrumErrSum->SetBinContent(i, errDR);
    }

    // purely calculating points based on all Systematic errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSys,histoCombDirGammaSpecStatErr,0,0.5);
    if(graphCombDirGammaSpectrumSystErr)graphCombDirGammaSpectrumSystErr->SetName("graphCombDirGammaSpectrumSystErr");
    if(graphCombDirGammaSpectrumSystErr) cout << "sys has been found" << endl;
    if(graphCombDirGammaSpectrumSystErr)graphCombDirGammaSpectrumSystErr->Print();

    // purely calculating points based on Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumStatErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrStat,histoCombDirGammaSpecStatErr,0,0.5);
    if(graphCombDirGammaSpectrumStatErr)graphCombDirGammaSpectrumStatErr->SetName("graphCombDirGammaSpectrumStatErr");
    if(graphCombDirGammaSpectrumStatErr) cout << "stat has been found" << endl;
    if(graphCombDirGammaSpectrumStatErr)graphCombDirGammaSpectrumStatErr->Print();
    // purely calculating points based on all Systematic + Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum,histoCombDirGammaSpecStatErr,0,0.5);
    if(graphCombDirGammaSpectrumSumErr)graphCombDirGammaSpectrumSumErr->SetName("graphCombDirGammaSpectrumSumErr");
    if(graphCombDirGammaSpectrumSumErr) cout << "tot has been found" << endl;
    if(graphCombDirGammaSpectrumSumErr)graphCombDirGammaSpectrumSumErr->Print();
    // calculate points above confidence level summed errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErrConfi = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum,histoCombDirGammaSpecStatErr,2,0.5);
    if(graphCombDirGammaSpectrumSumErrConfi)graphCombDirGammaSpectrumSumErrConfi->SetName("graphCombDirGammaSpectrumSumErrConfi");
    if(graphCombDirGammaSpectrumSumErrConfi) cout << "confi has been found" << endl;
    if(graphCombDirGammaSpectrumSumErrConfi)graphCombDirGammaSpectrumSumErrConfi->Print();
    // calculate arrows for points with 0, error summed
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErrAr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum,histoCombDirGammaSpecStatErr,5,0.5);
    if(graphCombDirGammaSpectrumSumErrAr)graphCombDirGammaSpectrumSumErrAr->SetName("graphCombDirGammaSpectrumSumErrAr");
    if(graphCombDirGammaSpectrumSumErrAr) cout << "Ar has been found" << endl;
    if(graphCombDirGammaSpectrumSumErrAr)graphCombDirGammaSpectrumSumErrAr->Print();

    //*******************************************************************************************************
    //************************** Calculating combined direct photon spectrum using pi0 fit in DR ************
    //*******************************************************************************************************
    Double_t xArrayCombPi0FitDR[graphCombDRStat->GetN()+1];
    xArrayCombPi0FitDR[0] = graphCombDRStat->GetX()[0] - graphCombDRStat->GetEXhigh()[0];
    for (Int_t i = 1; i<graphCombDRStat->GetN()+1;i++){
        xArrayCombPi0FitDR[i] = graphCombDRStat->GetX()[i-1] + graphCombDRStat->GetEXhigh()[i-1];
    }
    //_______________________ copy inclusive photon spectra _____________________
    TH1D *histoCombPi0FitDirGammaSpectrumErrSum                   = new TH1D("histoCombPi0FitDirGammaSpectrumErrSum","",graphCombDRStat->GetN(),xArrayCombPi0FitDR);
    TH1D *histoCombPi0FitDirGammaSpectrumErrSys                   = new TH1D("histoCombPi0FitDirGammaSpectrumErrSys","",graphCombDRStat->GetN(),xArrayCombPi0FitDR);
    TH1D *histoCombPi0FitDirGammaSpectrumErrStat                  = new TH1D("histoCombPi0FitDirGammaSpectrumErrStat","",graphCombDRStat->GetN(),xArrayCombPi0FitDR);

    //_______________________ get arrays of double ratio errors __________________
    Double_t *SystErrorsCombPi0FitDR                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *sumErrorsCombPi0FitDR                               = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *StatErrorsCombPi0FitDR                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *xErrorsPi0FitDR                                     = new Double_t[graphCombIncGammaStat->GetN()];
    for (Int_t i = 0; i< graphCombDRStat->GetN(); i++){
        SystErrorsCombPi0FitDR[i]                                 = graphCombDRSys->GetEYhigh()[i]/graphCombDRSys->GetY()[i] *100;
        StatErrorsCombPi0FitDR[i]                                 = graphCombDRStat->GetEYhigh()[i]/graphCombDRStat->GetY()[i] *100;
        sumErrorsCombPi0FitDR[i]                                  = graphCombDRTot->GetEYhigh()[i]/graphCombDRTot->GetY()[i] *100;
        //cout << i << "\t" << graphCombPi0FitDRSys->GetY()[i] << "\t" << graphCombPi0FitDRSys->GetEYhigh()[i] << "\t" <<SystErrorsCombPi0FitDR[i] << endl;
    }
    xErrorsPi0FitDR                                               = graphCombDRStat->GetX();

    cout << __LINE__ << endl;
    //graphCombPi0FitDRTot->Print();

    //_______________________ copy inclusive photon spectra _____________________
    TH1D* histoCombErrorsForPi0FitDRSum                           = new TH1D("histoCombErrorsForPi0FitDRSum","",graphCombDRStat->GetN(),xArrayCombPi0FitDR);
    TH1D* histoCombErrorsForPi0FitDRStat                          = new TH1D("histoCombErrorsForPi0FitDRStat","",graphCombDRStat->GetN(),xArrayCombPi0FitDR);
    TH1D* histoCombErrorsForPi0FitDRSys                           = new TH1D("histoCombErrorsForPi0FitDRSys","",graphCombDRStat->GetN(),xArrayCombPi0FitDR);

    for(Int_t i = 1; i<graphCombDRStat->GetN()+1;i++){
        //cout<< i << "\t"<<xErrorsPi0FitDR[i-1]<<"  "<<histoCombErrorsForPi0FitDRSum->GetBinCenter(i)<< "\t"<<histoCombErrorsForPi0FitDRSum->GetBinWidth(i) <<endl;
        Double_t binErrorSummed                                 = sumErrorsCombPi0FitDR[i-1];
        Double_t binErrorSyst                                   = SystErrorsCombPi0FitDR[i-1];
        Double_t binErrorStat                                   = StatErrorsCombPi0FitDR[i-1];
        Double_t DR                                             = graphCombDRStat->GetY()[i-1];

        //cout << DR << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
        histoCombErrorsForPi0FitDRSum->SetBinContent(i,DR);
        histoCombErrorsForPi0FitDRSys->SetBinContent(i,DR);
        histoCombErrorsForPi0FitDRStat->SetBinContent(i,DR);
        histoCombErrorsForPi0FitDRSum->SetBinError(i,(binErrorSummed/100)*DR);
        histoCombErrorsForPi0FitDRSys->SetBinError(i,(binErrorSyst/100)*DR);
        histoCombErrorsForPi0FitDRStat->SetBinError(i,(binErrorStat/100)*DR);
    }

    for(Int_t i = 1; i<histoCombErrorsForPi0FitDRSum->GetNbinsX()+1;i++){
        histoCombPi0FitDirGammaSpectrumErrSum->SetBinContent(i+1,-1);
        histoCombPi0FitDirGammaSpectrumErrSys->SetBinContent(i+1,-1);
        histoCombPi0FitDirGammaSpectrumErrStat->SetBinContent(i+1,-1);

        histoCombPi0FitDirGammaSpectrumErrSum->SetBinError(i+1,0);
        histoCombPi0FitDirGammaSpectrumErrSys->SetBinError(i+1,0);
        histoCombPi0FitDirGammaSpectrumErrStat->SetBinError(i+1,0);
    }

    // get the binning of the direct photons from the DR
    TH1D *histoCombPi0FitDirGammaSpecSysErr                        = new TH1D(*histoCombErrorsForPi0FitDRSys);
    TH1D *histoCombPi0FitDirGammaSpecStatErr                       = new TH1D(*histoCombErrorsForPi0FitDRStat);
    TH1D *histoCombPi0FitDirGammaSpecSumErr                        = new TH1D(*histoCombErrorsForPi0FitDRSum);


    cout << __LINE__ << endl;
    //graphCombDRStat->Print();
    for(Int_t i = 1; i<graphCombDRStat->GetN()+1; i++){
        // obtain common quantities
        Double_t Rgamma                 = histoCombErrorsForPi0FitDRSys->GetBinContent(i);
        Double_t nIncGamma              = graphCombIncGammaStat->GetY()[i-1];

        // calculating Systematics graph
        Double_t errRgamma              = histoCombErrorsForPi0FitDRSys->GetBinError(i);
        Double_t errNIncGam             = graphCombIncGammaSys->GetEYhigh()[i-1];
        Double_t q1                     = 1 - 1/ Rgamma;

        Double_t q1Error                = errRgamma/(Rgamma*Rgamma);
        Double_t content                = nIncGamma * ( 1 - 1/ Rgamma);
        Double_t error                  = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        Double_t errDR                  = content - error;
        histoCombPi0FitDirGammaSpecSysErr->SetBinError(i, error);
        histoCombPi0FitDirGammaSpecSysErr->SetBinContent(i, content);
        histoCombPi0FitDirGammaSpectrumErrSys->SetBinContent(i, errDR);

        // calculating Stat graphs
        errRgamma                       = histoCombErrorsForPi0FitDRStat->GetBinError(i);
        errNIncGam                      = graphCombIncGammaStat->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombPi0FitDirGammaSpecStatErr->SetBinError(i, error);
        histoCombPi0FitDirGammaSpecStatErr->SetBinContent(i, content);
        histoCombPi0FitDirGammaSpectrumErrStat->SetBinContent(i, errDR);

        // calculating summed error graphs
        errRgamma                       = histoCombErrorsForPi0FitDRSum->GetBinError(i);
        errNIncGam                      = graphCombIncGammaTot->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombPi0FitDirGammaSpecSumErr->SetBinError(i, error);
        histoCombPi0FitDirGammaSpecSumErr->SetBinContent(i, content);
        histoCombPi0FitDirGammaSpectrumErrSum->SetBinContent(i, errDR);
    }

    // purely calculating points based on all Systematic errors
    TGraphAsymmErrors *graphCombPi0FitDirGammaSpectrumSystErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombPi0FitDirGammaSpectrumErrSys,histoCombPi0FitDirGammaSpecStatErr,0,0.5);
    if(graphCombPi0FitDirGammaSpectrumSystErr)graphCombPi0FitDirGammaSpectrumSystErr->SetName("graphCombPi0FitDirGammaSpectrumSystErr");
    if(graphCombPi0FitDirGammaSpectrumSystErr) cout << "sys has been found" << endl;
    if(graphCombPi0FitDirGammaSpectrumSystErr)graphCombPi0FitDirGammaSpectrumSystErr->Print();

    // purely calculating points based on Statistical errors
    TGraphAsymmErrors *graphCombPi0FitDirGammaSpectrumStatErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombPi0FitDirGammaSpectrumErrStat,histoCombPi0FitDirGammaSpecStatErr,0,0.5);
    if(graphCombPi0FitDirGammaSpectrumStatErr)graphCombPi0FitDirGammaSpectrumStatErr->SetName("graphCombPi0FitDirGammaSpectrumStatErr");
    if(graphCombPi0FitDirGammaSpectrumStatErr) cout << "stat has been found" << endl;
    if(graphCombPi0FitDirGammaSpectrumStatErr)graphCombPi0FitDirGammaSpectrumStatErr->Print();
    // purely calculating points based on all Systematic + Statistical errors
    TGraphAsymmErrors *graphCombPi0FitDirGammaSpectrumSumErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombPi0FitDirGammaSpectrumErrSum,histoCombPi0FitDirGammaSpecStatErr,0,0.5);
    if(graphCombPi0FitDirGammaSpectrumSumErr)graphCombPi0FitDirGammaSpectrumSumErr->SetName("graphCombPi0FitDirGammaSpectrumSumErr");
    if(graphCombPi0FitDirGammaSpectrumSumErr) cout << "tot has been found" << endl;
    if(graphCombPi0FitDirGammaSpectrumSumErr)graphCombPi0FitDirGammaSpectrumSumErr->Print();
    // calculate points above confidence level summed errors
    TGraphAsymmErrors *graphCombPi0FitDirGammaSpectrumSumErrConfi = CalculateDirectPhotonPointsAndUpperLimits(histoCombPi0FitDirGammaSpectrumErrSum,histoCombPi0FitDirGammaSpecStatErr,2,0.5);
    if(graphCombPi0FitDirGammaSpectrumSumErrConfi)graphCombPi0FitDirGammaSpectrumSumErrConfi->SetName("graphCombPi0FitDirGammaSpectrumSumErrConfi");
    if(graphCombPi0FitDirGammaSpectrumSumErrConfi) cout << "confi has been found" << endl;
    if(graphCombPi0FitDirGammaSpectrumSumErrConfi)graphCombPi0FitDirGammaSpectrumSumErrConfi->Print();
    // calculate arrows for points with 0, error summed
    TGraphAsymmErrors *graphCombPi0FitDirGammaSpectrumSumErrAr = CalculateDirectPhotonPointsAndUpperLimits(histoCombPi0FitDirGammaSpectrumErrSum,histoCombPi0FitDirGammaSpecStatErr,5,0.5);
    if(graphCombPi0FitDirGammaSpectrumSumErrAr)graphCombPi0FitDirGammaSpectrumSumErrAr->SetName("graphCombPi0FitDirGammaSpectrumSumErrAr");
    if(graphCombPi0FitDirGammaSpectrumSumErrAr) cout << "Ar has been found" << endl;
    if(graphCombPi0FitDirGammaSpectrumSumErrAr)graphCombPi0FitDirGammaSpectrumSumErrAr->Print();

    //*******************************************************************************************************************************************
    //********************************************** Plotting individual Rgamma ratios **********************************************************
    //*******************************************************************************************************************************************

    // double ratio combined
    TCanvas *canvasDoubleRatio = new TCanvas("canvasDoubleRatio","",0.095,0.09,1000,815);
    DrawGammaCanvasSettings( canvasDoubleRatio, 0.086, 0.01, 0.01, 0.105);
    canvasDoubleRatio->cd();
    canvasDoubleRatio->SetLogx();

    Double_t minY                            = 0.85;
    Double_t maxY                            = 1.65;
    Double_t textSizeSinglePad               = 0.05;
    TH2F * hist2DDRDummySingle       = new TH2F("hist2DDRDummySingle","hist2DDRDummySingle",1000,doubleRatioXpp[0], doubleRatioXpp[1],1000,doubleRatio[0], doubleRatio[1]);
    SetStyleHistoTH2ForGraphs(hist2DDRDummySingle, "#it{p}_{T} (GeV/#it{c})","#it{R}_{#gamma}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
    hist2DDRDummySingle->GetXaxis()->SetLabelOffset(-0.01);
    hist2DDRDummySingle->GetXaxis()->SetMoreLogLabels(kTRUE);
    hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRSingle = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*4,0.35,0.96, textSizeSinglePad, 1, "", 42, 0.3);
        legendDRSingle->SetTextAlign(11);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 0; i < 10; i++){
            if (graphDRPi0FitSysErr[i]){
                DrawGammaSetMarkerTGraphAsym(graphDRPi0FitSysErr[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphDRPi0FitSysErr[i]->Draw("E2same");
                legendDRSingle->AddEntry(graphDRPi0FitSysErr[i],nameMeasGlobalLabel[i],"pf");

            }
            if (histoDRPi0FitStatErr[i]){
                DrawGammaSetMarker(histoDRPi0FitStatErr[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoDRPi0FitStatErr[i]->Draw("p,same,e0,X0");
                if (!graphDRPi0FitSysErr[i])legendDRSingle->AddEntry(histoDRPi0FitStatErr[i],nameMeasGlobalLabel[i],"p");
            }
        }

        if (histoDRPi0FitStatErr[2]) histoDRPi0FitStatErr[2]->Draw("p,same,e0,X0");
        if (histoDRPi0FitStatErr[0]) histoDRPi0FitStatErr[0]->Draw("p,same,e0,X0");
        if (histoDRPi0FitStatErr[4]) histoDRPi0FitStatErr[4]->Draw("p,same,e0,X0");
        legendDRSingle->Draw();

        TLatex *labelDRSingle = new TLatex(0.95,0.92,collisionSystempPbNSD.Data());
        SetStyleTLatex( labelDRSingle, textSizeSinglePad,4, 1, 42, kTRUE, 31);
        labelDRSingle->Draw();
        TLatex *labelALICEDRSingle = new TLatex(0.95,0.87,textALICE.Data());
        SetStyleTLatex( labelALICEDRSingle, textSizeSinglePad,4, 1, 42, kTRUE, 31);
        labelALICEDRSingle->Draw();


        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pPb5TeV.%s", outputDirFit.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pPb5TeV.pdf", outputDirFit.Data()));

    hist2DDRDummySingle->DrawCopy();
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 0; i < 10; i++){
            if (graphDRNonFitSysErr[i]){
                DrawGammaSetMarkerTGraphAsym(graphDRNonFitSysErr[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphDRNonFitSysErr[i]->Draw("E2same");
            }
            if (histoDRNonFitStatErr[i]){
                DrawGammaSetMarker(histoDRNonFitStatErr[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoDRNonFitStatErr[i]->Draw("p,same,e0,X0");
            }
        }

        if (histoDRNonFitStatErr[2]) histoDRNonFitStatErr[2]->Draw("p,same,e0,X0");
        if (histoDRNonFitStatErr[0]) histoDRNonFitStatErr[0]->Draw("p,same,e0,X0");
        if (histoDRNonFitStatErr[4]) histoDRNonFitStatErr[4]->Draw("p,same,e0,X0");
        legendDRSingle->Draw();

        labelDRSingle->Draw();
        labelALICEDRSingle->Draw();

        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DRNonFit_IndMeasurements_pPb5TeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DRNonFit_IndMeasurements_pPb5TeV.pdf", outputDir.Data()));


    hist2DDRDummySingle->DrawCopy();

        TGraphAsymmErrors* graphCombDRStatPlot    = (TGraphAsymmErrors*)graphCombDRStat->Clone("graphCombDRStatPlot");
        ProduceGraphAsymmWithoutXErrors(graphCombDRStatPlot);

        TLegend* legendDRComb = GetAndSetLegend2(0.12,0.95-textSizeSinglePad*1,0.5,0.95, textSizeSinglePad, 1, "", 42, 0.15);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        DrawGammaSetMarkerTGraphAsym(graphCombDRSys, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRStatPlot, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes);
        legendDRComb->AddEntry(graphCombDRSys,"ALICE prel.","pf");
        graphCombDRSys->Draw("E2same");
        graphCombDRStatPlot->Draw("z,p,same");
//         legendDRComb->Draw();

        labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_Comb_pPb5TeV.%s", outputDirFit.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_Comb_pPb5TeV.pdf", outputDirFit.Data()));

    hist2DDRDummySingle->DrawCopy();
        TGraphAsymmErrors* graphCombDRNonFitStatPlot    = (TGraphAsymmErrors*)graphCombDRNonFitStat->Clone("graphCombDRNonFitStatPlot");
        ProduceGraphAsymmWithoutXErrors(graphCombDRNonFitStatPlot);

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitSys, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitStatPlot, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes);
        graphCombDRNonFitSys->Draw("E2same");
        graphCombDRNonFitStatPlot->Draw("z,p,same");

        labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DRNonFit_Comb_pPb5TeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DRNonFit_Comb_pPb5TeV.pdf", outputDir.Data()));

        hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRTheoryComb     = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*2,0.5,0.96, textSizeSinglePad, 1, "", 42, 0.15);
        TLegend* legendDRTheoryComb2    = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*6,0.5,0.96-textSizeSinglePad*2, textSizeSinglePad, 1, "NLO pQCD: ", 42, 0.15);
        TGraphAsymmErrors* dummyNLO         = new TGraphAsymmErrors(1);
        TGraphAsymmErrors* dummyMcGill      = new TGraphAsymmErrors(1);
        TGraphAsymmErrors* dummyNLOnCTEQ    = new TGraphAsymmErrors(1);
        TGraphAsymmErrors* dummyNLOEPPS     = new TGraphAsymmErrors(1);

        legendDRTheoryComb->AddEntry(graphCombDRSys,"ALICE prel.","pf");

        if (graphTheoryNLODRpPb) {
            DrawGammaSetMarkerTGraphAsym(dummyNLO, 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
            dummyNLO->SetLineStyle(styleLineNLOWerner);
            dummyNLO->SetLineWidth(widthLineNLO);
            dummyNLO->SetLineColor(colorNLOWerner);
            DrawGammaSetMarkerTGraphAsym(graphTheoryNLODRpPb, 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
            graphTheoryNLODRpPb->Draw("3,same");
            legendDRTheoryComb2->AddEntry(dummyNLO,"PDF: CT10, FF: GRV","fl");
        }

        if (graphTheoryMCGillDRpPb) {
            DrawGammaSetMarkerTGraphAsym(dummyMcGill, 0, 0, colorNLOMcGill, colorNLOMcGill, 0.6, kTRUE, colorNLOMcGill, kTRUE);
            dummyMcGill->SetLineStyle(styleLineMcGill);
            dummyMcGill->SetLineWidth(widthLineNLO);
            dummyMcGill->SetLineColor(colorNLOMcGill);
            DrawGammaSetMarkerTGraphAsym(graphTheoryMCGillDRpPb, 0, 0, colorNLOMcGill, colorNLOMcGill, 0.6, kTRUE, colorNLOMcGill, kTRUE);
            graphTheoryMCGillDRpPb->Draw("3,same");
            legendDRTheoryComb->AddEntry(dummyMcGill,"Shen #it{et al.}","fl");
        }
        if (graphTheoryMCGillDRpPbCenter){
            DrawGammaNLOTGraph( graphTheoryMCGillDRpPbCenter, widthLineNLO, styleLineMcGill, colorNLOMcGill);
            graphTheoryMCGillDRpPbCenter->Draw("lc,same");
        }
        if (graphTheoryNLODRpPbCenter){
            DrawGammaNLOTGraph( graphTheoryNLODRpPbCenter, widthLineNLO, styleLineNLOWerner, colorNLOWerner);
            graphTheoryNLODRpPbCenter->Draw("lc,same");
        }
        if (graphTheoryPowhegDRnCTEQpPb) {
            DrawGammaSetMarkerTGraphAsym(dummyNLOnCTEQ, 0, 0, colorNLONCTEQBand, colorNLONCTEQBand, 0.6, kTRUE, colorNLONCTEQBand, kTRUE);
            dummyNLOnCTEQ->SetLineStyle(styleLineNLONCTEQ);
            dummyNLOnCTEQ->SetLineWidth(widthLineNLO);
            dummyNLOnCTEQ->SetLineColor(colorNLONCTEQ);
            DrawGammaSetMarkerTGraphAsym(graphTheoryPowhegDRnCTEQpPb, 0, 0, colorNLONCTEQBand, colorNLONCTEQBand, 0.6, kTRUE, colorNLONCTEQBand, kTRUE);
            graphTheoryPowhegDRnCTEQpPb->Draw("3,same");
            legendDRTheoryComb2->AddEntry(dummyNLOnCTEQ,"nPDF: nCTEQ15, FF: GRV","fl");
        }
        if (graphTheoryPowhegDRnCTEQpPbCenter){
            DrawGammaNLOTGraph( graphTheoryPowhegDRnCTEQpPbCenter, widthLineNLO, styleLineNLONCTEQ, colorNLONCTEQ);
            graphTheoryPowhegDRnCTEQpPbCenter->Draw("lc,same");
        }
        if (graphTheoryPowhegDREPPS16pPb) {
            DrawGammaSetMarkerTGraphAsym(dummyNLOEPPS, 0, 0, colorNLOEPPSBand, colorNLOEPPSBand, 0.6, kTRUE, colorNLOEPPSBand, kTRUE);
            dummyNLOEPPS->SetLineStyle(styleLineNLOEPPS);
            dummyNLOEPPS->SetLineWidth(widthLineNLO);
            dummyNLOEPPS->SetLineColor(colorNLOEPPS);
            DrawGammaSetMarkerTGraphAsym(graphTheoryPowhegDREPPS16pPb, 0, 0, colorNLOEPPSBand, colorNLOEPPSBand, 0.6, kTRUE, colorNLOEPPSBand, kTRUE);
            graphTheoryPowhegDREPPS16pPb->Draw("3,same");
            legendDRTheoryComb2->AddEntry(dummyNLOEPPS,"nPDF: EPPS16, FF: GRV","fl");
        }
        if (graphTheoryPowhegDREPPS16pPbCenter){
            DrawGammaNLOTGraph( graphTheoryPowhegDREPPS16pPbCenter, widthLineNLO, styleLineNLOEPPS, colorNLOEPPS);
            graphTheoryPowhegDREPPS16pPbCenter->Draw("lc,same");
        }

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRSys->Draw("E2same");
        graphCombDRStatPlot->Draw("p,z,same");
        legendDRTheoryComb->Draw();
        legendDRTheoryComb2->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pPb5TeV.%s", outputDirFit.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pPb5TeV.pdf", outputDirFit.Data()));

        hist2DDRDummySingle->DrawCopy();
        if (graphTheoryNLODRpPb) {
            graphTheoryNLODRpPb->Draw("3,same");
        }

        if (graphTheoryMCGillDRpPb) {
            graphTheoryMCGillDRpPb->Draw("3,same");
        }
        if (graphTheoryMCGillDRpPbCenter){
            graphTheoryMCGillDRpPbCenter->Draw("lc,same");
        }
        if (graphTheoryNLODRpPbCenter){
            graphTheoryNLODRpPbCenter->Draw("lc,same");
        }
        if (graphTheoryPowhegDRnCTEQpPb) {
            graphTheoryPowhegDRnCTEQpPb->Draw("3,same");
        }
        if (graphTheoryPowhegDRnCTEQpPbCenter){
            graphTheoryPowhegDRnCTEQpPbCenter->Draw("lc,same");
        }
        if (graphTheoryPowhegDREPPS16pPb) {
            graphTheoryPowhegDREPPS16pPb->Draw("3,same");
        }
        if (graphTheoryPowhegDREPPS16pPbCenter){
            graphTheoryPowhegDREPPS16pPbCenter->Draw("lc,same");
        }

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRNonFitSys->Draw("E2same");
        graphCombDRNonFitStatPlot->Draw("p,z,same");
        legendDRTheoryComb->Draw();
        legendDRTheoryComb2->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DRNonFit_CombAndTheory_pPb5TeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DRNonFit_CombAndTheory_pPb5TeV.pdf", outputDir.Data()));


    hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRTheoryComb5        = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*1,0.5,0.96, textSizeSinglePad, 1, "", 42, 0.15);
        TLegend* legendDRTheoryComb6        = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*5,0.5,0.96-textSizeSinglePad*1, textSizeSinglePad, 1, "", 42, 0.15);

        legendDRTheoryComb5->AddEntry(graphCombDRNonFitSys,"ALICE prel.","pf");

        if (graphTheoryDRAlwinaPeTer){
            DrawGammaNLOTGraph( graphTheoryDRAlwinaPeTer, widthLineNLO, styleLineMcGill, colorNLOMcGill);
            graphTheoryDRAlwinaPeTer->Draw("lc,same");
            legendDRTheoryComb6->AddEntry(graphTheoryDRAlwinaPeTer,"PeTer","l");
        }
        if (graphTheoryDRAlwinaPWGGA){
            DrawGammaNLOTGraph( graphTheoryDRAlwinaPWGGA, widthLineNLO, styleLineNLONCTEQ, colorNLONCTEQ);
            graphTheoryDRAlwinaPWGGA->Draw("lc,same");
            legendDRTheoryComb6->AddEntry(graphTheoryDRAlwinaPWGGA,"PWGGA","l");
        }
        if (graphTheoryDRAlwinaJetPHOX){
            DrawGammaNLOTGraph( graphTheoryDRAlwinaJetPHOX, widthLineNLO, styleLineNLOEPPS, colorNLOEPPS);
            graphTheoryDRAlwinaJetPHOX->Draw("lc,same");
            legendDRTheoryComb6->AddEntry(graphTheoryDRAlwinaJetPHOX,"JetPHOX","l");
        }
        if (graphTheoryDRAlwinaVogelsang){
            DrawGammaNLOTGraph( graphTheoryDRAlwinaVogelsang, widthLineNLO, styleLineNLOWerner, colorNLOWerner);
            graphTheoryDRAlwinaVogelsang->Draw("lc,same");
            legendDRTheoryComb6->AddEntry(graphTheoryDRAlwinaVogelsang,"W. Vogelsang","l");
        }

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRNonFitSys->Draw("E2same");
        graphCombDRNonFitStatPlot->Draw("p,z,same");
        legendDRTheoryComb5->Draw();
        legendDRTheoryComb6->Draw();

        labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DRNonFit_CombAndTheory_pPb5TeV_Alwina.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DRNonFit_CombAndTheory_pPb5TeV_Alwina.pdf", outputDir.Data()));


        hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRTheoryComb3    = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*1,0.5,0.96, textSizeSinglePad, 1, "", 42, 0.15);
        TLegend* legendDRTheoryComb4    = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*5,0.5,0.96-textSizeSinglePad*1, textSizeSinglePad, 1, "NLO pQCD: ", 42, 0.15);

        if (graphTheoryNLODRpPb) {
            graphTheoryNLODRpPb->Draw("3,same");
            legendDRTheoryComb4->AddEntry(dummyNLO,"PDF: CT10, FF: GRV","fl");
        }

        if (graphTheoryMCGillDRpPb) {
            graphTheoryMCGillDRpPb->Draw("3,same");
            legendDRTheoryComb3->AddEntry(dummyMcGill,"Shen #it{et al.}","fl");
        }
        if (graphTheoryMCGillDRpPbCenter){
            graphTheoryMCGillDRpPbCenter->Draw("lc,same");
        }
        if (graphTheoryNLODRpPbCenter){
            graphTheoryNLODRpPbCenter->Draw("lc,same");
        }
        if (graphTheoryPowhegDRnCTEQpPb) {
            graphTheoryPowhegDRnCTEQpPb->Draw("3,same");
            legendDRTheoryComb4->AddEntry(dummyNLOnCTEQ,"nPDF: nCTEQ15, FF: GRV","fl");
        }
        if (graphTheoryPowhegDRnCTEQpPbCenter){
            graphTheoryPowhegDRnCTEQpPbCenter->Draw("lc,same");
        }
        if (graphTheoryPowhegDREPPS16pPb) {
            graphTheoryPowhegDREPPS16pPb->Draw("3,same");
            legendDRTheoryComb4->AddEntry(dummyNLOEPPS,"nPDF: EPPS16, FF: GRV","fl");
        }
        if (graphTheoryPowhegDREPPS16pPbCenter){
            graphTheoryPowhegDREPPS16pPbCenter->Draw("lc,same");
        }

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        legendDRTheoryComb3->Draw();
        legendDRTheoryComb4->Draw();

//         labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_Theory_pPb5TeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_Theory_pPb5TeV.pdf", outputDir.Data()));


    // **********************************************************************************************************************
    // ******************************** Efficiency for gamma individual measurements ****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 55;
    Double_t textSizeLabelsRel          = 55./1200;
    cout << textSizeLabelsRel << endl;

    TCanvas* canvasEff   = new TCanvas("canvasEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEff,  0.09, 0.01, 0.015, 0.095);
//     canvasEff->SetLogy(1);
    canvasEff->SetLogx(1);

    TH1F * histo1DEff            = new TH1F("histo1DEff", "histo1DEff",1000, doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(  histo1DEff, "#it{p}_{T} (GeV/#it{c})","#it{#varepsilon}_{rec}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 0.98);//(#times #epsilon_{pur})
    histo1DEff->GetYaxis()->SetRangeUser(0.18, 0.82 );
    histo1DEff->GetYaxis()->SetLabelOffset(0.001);
    histo1DEff->GetXaxis()->SetLabelOffset(-0.01);
    histo1DEff->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo1DEff->DrawCopy();

        TLegend* legendEffiGamma           = GetAndSetLegend2(0.57, 0.13, 0.95, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel, 2);
        for (Int_t i = 0; i < 11; i++){
            if (histoEffi[i]){
                DrawGammaSetMarker(histoEffi[i],   markerStyleDetMC[i], markerSizeDetMC[i], colorDetMC[i] , colorDetMC[i]);
                histoEffi[i]->Draw("p,same,e");
                legendEffiGamma->AddEntry(histoEffi[i],"   ","p");
//                 if (i == 1){
//                     DrawGammaSetMarker(histoEffiPHI7,   markerStyleDetMC[2], markerSizeDetMC[i], colorDetMC[i] , colorDetMC[i]);
//                     histoEffiPHI7->Draw("p,same,e");
//                 }
            } else if (histoEffiMCPt[i]){
                legendEffiGamma->AddEntry((TObject*)0,"   ","");
            }
            if (histoEffiMCPt[i]){
                DrawGammaSetMarker(histoEffiMCPt[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoEffiMCPt[i]->Draw("p,same,e");
                legendEffiGamma->AddEntry(histoEffiMCPt[i],"    "+nameMeasGlobalLabel[i],"p");
            } else if (histoEffi[i]){
                legendEffiGamma->AddEntry((TObject*)0,"    "+nameMeasGlobalLabel[i],"");
            }
        }
        legendEffiGamma->Draw();

        TLatex *labelPerfEffi           = new TLatex(0.13,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
        labelPerfEffi->Draw();
        TLatex *labelEnergyEffi         = new TLatex(0.13,0.87,collisionSystempPbNSD.Data());
        SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
        labelEnergyEffi->Draw();
        TLatex *labelPerfEffiPTrec      = new TLatex(0.57,0.145+(4*textSizeLabelsRel),"#it{p}_{T}^{rec}");
        SetStyleTLatex( labelPerfEffiPTrec, textSizeLabelsRel,4);
        labelPerfEffiPTrec->Draw();
        TLatex *labelPerfEffiPTtrue      = new TLatex(0.665,0.145+(4*textSizeLabelsRel),"#it{p}_{T}^{true}");
        SetStyleTLatex( labelPerfEffiPTtrue, textSizeLabelsRel,4);
        labelPerfEffiPTtrue->Draw();

    canvasEff->Update();
    canvasEff->Print(Form("%s/Gamma_Effiency.%s",outputDir.Data(),suffix.Data()));
    canvasEff->Print(Form("%s/Gamma_Effiency.pdf",outputDir.Data()));

    // **********************************************************************************************************************
    // ******************************** ResolutionCorr for gamma individual measurements ****************************************
    // **********************************************************************************************************************
    TCanvas* canvasResolCor   = new TCanvas("canvasResolCor", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasResolCor,  0.1, 0.01, 0.015, 0.095);
//     canvasResolCor->SetLogy(1);
    canvasResolCor->SetLogx(1);

    TH2F* histo2DResCor            = new TH2F("histo1DResCor", "histo1DResCor",1000, doubleRatioXpp[0], doubleRatioXpp[1], 1000, -1, 5);
    SetStyleHistoTH2ForGraphs(  histo2DResCor, "#it{p}_{T} (GeV/#it{c})","#it{#varepsilon}_{resol}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
    histo2DResCor->GetYaxis()->SetRangeUser(0.2, 2.6 );
    histo2DResCor->GetYaxis()->SetLabelOffset(0.001);
    histo2DResCor->GetXaxis()->SetLabelOffset(-0.01);
    histo2DResCor->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DResCor->DrawCopy();

        TLegend* legendResolCorGamma           = GetAndSetLegend2(0.15, 0.85-(3*textSizeLabelsRel), 0.33, 0.85,textSizeLabelsPixel);
        for (Int_t i = 0; i < 11; i++){
            if (histoResolCorr[i]){
                DrawGammaSetMarker(histoResolCorr[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoResolCorr[i]->GetXaxis()->SetRangeUser(graphIndGammaIncStat[i]->GetX()[0]-graphIndGammaIncStat[i]->GetEXlow()[0],
                                                            graphIndGammaIncStat[i]->GetX()[graphIndGammaIncStat[i]->GetN()-1]+graphIndGammaIncStat[i]->GetEXhigh()[graphIndGammaIncStat[i]->GetN()-1] );

                histoResolCorr[i]->Draw("p,same,e");
                legendResolCorGamma->AddEntry(histoResolCorr[i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendResolCorGamma->Draw();
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        TLatex *labelPerfResolCor           = new TLatex(0.15,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfResolCor, textSizeLabelsRel,4);
        labelPerfResolCor->Draw();
        TLatex *labelEnergyResolCor         = new TLatex(0.15,0.87,collisionSystempPbNSD.Data());
        SetStyleTLatex( labelEnergyResolCor, textSizeLabelsRel,4);
        labelEnergyResolCor->Draw();

    histo2DResCor->Draw("same,axis");
    canvasResolCor->Update();
    canvasResolCor->Print(Form("%s/Gamma_ResolutionCorrection.%s",outputDir.Data(),suffix.Data()));
    canvasResolCor->Print(Form("%s/Gamma_ResolutionCorrection.pdf",outputDir.Data()));

    // **********************************************************************************************************************
    // ******************************** Purity for gamma individual measurements ****************************************
    // **********************************************************************************************************************
    TCanvas* canvasPurity   = new TCanvas("canvasPurity", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasPurity,  0.1, 0.01, 0.015, 0.095);
    //     canvasPurity->SetLogy(1);
    canvasPurity->SetLogx(1);

    TH1F * histo1DPurity            = new TH1F("histo1DPurity", "histo1DPurity",1000, doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(  histo1DPurity, "#it{p}_{T} (GeV/#it{c})","#it{P}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
    histo1DPurity->GetYaxis()->SetRangeUser(0.8, 1.07 );
    histo1DPurity->GetYaxis()->SetLabelOffset(0.001);
    histo1DPurity->GetXaxis()->SetLabelOffset(-0.01);
    histo1DPurity->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo1DPurity->DrawCopy();

        TLegend* legendPurityGamma           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 0; i < 11; i++){
            if (histoPurity[i]){
                DrawGammaSetMarker(histoPurity[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoPurity[i]->Draw("p,same,e");
                legendPurityGamma->AddEntry(histoPurity[i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendPurityGamma->Draw();

        TLatex *labelPerfPurity           = new TLatex(0.15,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfPurity, textSizeLabelsRel,4);
        labelPerfPurity->Draw();
        TLatex *labelEnergyPurity         = new TLatex(0.15,0.87,collisionSystempPbNSD.Data());
        SetStyleTLatex( labelEnergyPurity, textSizeLabelsRel,4);
        labelEnergyPurity->Draw();

    canvasPurity->Update();
    canvasPurity->Print(Form("%s/Gamma_Purity.%s",outputDir.Data(),suffix.Data()));
    canvasPurity->Print(Form("%s/Gamma_Purity.pdf",outputDir.Data()));

    // **********************************************************************************************************************
    // ******************************** ConvProb for gamma individual measurements ******************************************
    // **********************************************************************************************************************
    TF1* fitConvProbPCM                     = new TF1("fitConvProbPCM","[0]");
    if (histoConvProb[0]) histoConvProb[0]->Fit(fitConvProbPCM,"NRMEX0+","",4,10.);
    else if  (histoConvProb[4]) histoConvProb[4]->Fit(fitConvProbPCM,"NRMEX0+","",4,10.);
    cout << WriteParameterToFile(fitConvProbPCM)<< endl;

    TCanvas* canvasConvProb   = new TCanvas("canvasConvProb", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasConvProb,  0.11, 0.01, 0.015, 0.095);
    //     canvasConvProb->SetLogy(1);
    canvasConvProb->SetLogx(1);

    TH1F * histo1DConvProb            = new TH1F("histo1DConvProb", "histo1DConvProb",1000, doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(  histo1DConvProb, "#it{p}_{T} (GeV/#it{c})","#it{P}_{conv}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.2);//(#times #epsilon_{pur})
    histo1DConvProb->GetYaxis()->SetRangeUser(0.04, 0.12 );
    histo1DConvProb->GetYaxis()->SetLabelOffset(0.003);
    histo1DConvProb->GetXaxis()->SetLabelOffset(-0.01);
    histo1DConvProb->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo1DConvProb->DrawCopy();

        TBox* boxConvProbPCM                = CreateBoxConv( kGray, 4,fitConvProbPCM->GetParameter(0)-fitConvProbPCM->GetParError(0), doubleRatioXpp[1],
                                                            fitConvProbPCM->GetParameter(0)+fitConvProbPCM->GetParError(0));
        boxConvProbPCM->Draw();

        DrawGammaSetMarkerTF1( fitConvProbPCM, 7, 2, kGray+2);
        fitConvProbPCM->SetRange(4,doubleRatioXpp[1]);
        fitConvProbPCM->Draw("same");


        TLegend* legendConvProbGamma           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(2*textSizeLabelsRel),textSizeLabelsPixel);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 0; i < 11; i++){
            if (histoConvProb[i]){
                DrawGammaSetMarker(histoConvProb[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoConvProb[i]->Draw("p,same,e");
                legendConvProbGamma->AddEntry(histoConvProb[i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendConvProbGamma->Draw();

        TLatex *labelPerfConvProb           = new TLatex(0.15,0.92,"ALICE simulation");
        SetStyleTLatex( labelPerfConvProb, textSizeLabelsRel,4);
        labelPerfConvProb->Draw();
        TLatex *labelEnergyConvProb         = new TLatex(0.15,0.87,collisionSystempPbNSD.Data());
        SetStyleTLatex( labelEnergyConvProb, textSizeLabelsRel,4);
        labelEnergyConvProb->Draw();

    histo1DConvProb->Draw("same,axis");

    canvasConvProb->Update();
    canvasConvProb->Print(Form("%s/Gamma_ConvProb.%s",outputDir.Data(),suffix.Data()));
    canvasConvProb->Print(Form("%s/Gamma_ConvProb.pdf",outputDir.Data()));
    // **********************************************************************************************************************
    // ******************************** PileupCorr for gamma individual measurements ****************************************
    // **********************************************************************************************************************

    TCanvas* canvasPileUp   = new TCanvas("canvasPileUp", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasPileUp,  0.115, 0.01, 0.015, 0.095);
    //     canvasPileUp->SetLogy(1);
    canvasPileUp->SetLogx(1);

    TH1F * histo1DPileUp            = new TH1F("histo1DPileUp", "histo1DPileUp",1000, doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(  histo1DPileUp, "#it{p}_{T} (GeV/#it{c})","#it{C}_{pile-up}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.25, 510,505);
    histo1DPileUp->GetYaxis()->SetRangeUser(0.94, 1.02 );
    histo1DPileUp->GetYaxis()->SetLabelOffset(0.005);
    histo1DPileUp->GetXaxis()->SetLabelOffset(-0.01);
    histo1DPileUp->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo1DPileUp->DrawCopy();

        TLegend* legendPileUpGamma           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(2*textSizeLabelsRel),textSizeLabelsPixel);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 0; i < 11; i++){
            if (histoPileupCorr[i]){
                DrawGammaSetMarker(histoPileupCorr[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoPileupCorr[i]->Draw("p,same,e");
                legendPileUpGamma->AddEntry(histoPileupCorr[i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendPileUpGamma->Draw();

        TLatex *labelPerfPileUp           = new TLatex(0.15,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfPileUp, textSizeLabelsRel,4);
        labelPerfPileUp->Draw();
        TLatex *labelEnergyPileUp         = new TLatex(0.15,0.87,collisionSystempPbNSD.Data());
        SetStyleTLatex( labelEnergyPileUp, textSizeLabelsRel,4);
        labelEnergyPileUp->Draw();

    canvasPileUp->Update();
    canvasPileUp->Print(Form("%s/Gamma_PileUp.%s",outputDir.Data(),suffix.Data()));
    canvasPileUp->Print(Form("%s/Gamma_PileUp.pdf",outputDir.Data()));

    // **********************************************************************************************************************
    // ******************************** EffectiveSecCorr for gamma individual measurements ****************************************
    // **********************************************************************************************************************
    TCanvas* canvasEffectiveSecCorr   = new TCanvas("canvasEffectiveSecCorr", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEffectiveSecCorr,  0.11, 0.01, 0.04, 0.095);
    //     canvasEffectiveSecCorr->SetLogy(1);
    canvasEffectiveSecCorr->SetLogx(1);

    TH1F * histo1DEffectiveSecCorr            = new TH1F("histo1DEffectiveSecCorr", "histo1DEffectiveSecCorr",1000, doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(  histo1DEffectiveSecCorr, "#it{p}_{T} (GeV/#it{c})","#it{C}_{sec}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.08);//(#times #epsilon_{pur})
    histo1DEffectiveSecCorr->GetYaxis()->SetLabelOffset(0.003);
    histo1DEffectiveSecCorr->GetXaxis()->SetLabelOffset(-0.01);
    histo1DEffectiveSecCorr->GetXaxis()->SetMoreLogLabels(kTRUE);


    Double_t minYSecCorr[4]             = {0.0, 0.0, 0.0, 0.0};
    Double_t maxYSecCorr[4]             = {0.05, 0.004, 0.8e-3, 0.041};
    TString nameLabelSec[4]             = {"K^{0}_{s}", "K^{0}_{l}", "#Lambda", "Rest"};
    TString nameOutputSec[4]            = {"K0s", "K0l", "Lambda", "Rest"};
    for (Int_t k = 0; k < 4; k++){
        histo1DEffectiveSecCorr->GetYaxis()->SetTitle(Form("C_{sec,%s}", nameLabelSec[k].Data()));
        histo1DEffectiveSecCorr->GetYaxis()->SetRangeUser(minYSecCorr[k], maxYSecCorr[k] );
        histo1DEffectiveSecCorr->DrawCopy();
        TLegend* legendEffectiveSecCorrGamma            = GetAndSetLegend2(0.65, 0.925-(3*textSizeLabelsRel), 0.93, 0.925,textSizeLabelsPixel);
        for (Int_t i = 0; i < 11; i++){
            if (histoEffSecCorr[k][i]){
                DrawGammaSetMarker(histoEffSecCorr[k][i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoEffSecCorr[k][i]->Draw("p,same,e");
                legendEffectiveSecCorrGamma->AddEntry(histoEffSecCorr[k][i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendEffectiveSecCorrGamma->Draw();

        TLatex *labelPerfEffectiveSecCorr           = new TLatex(0.15,0.89,"ALICE performance");
        SetStyleTLatex( labelPerfEffectiveSecCorr, textSizeLabelsRel,4);
        labelPerfEffectiveSecCorr->Draw();
        TLatex *labelEnergyEffectiveSecCorr         = new TLatex(0.15,0.84,collisionSystempPbNSD.Data());
        SetStyleTLatex( labelEnergyEffectiveSecCorr, textSizeLabelsRel,4);
        labelEnergyEffectiveSecCorr->Draw();

        canvasEffectiveSecCorr->Update();
        canvasEffectiveSecCorr->Print(Form("%s/Gamma_EffectiveSecCorr_%s.%s",outputDir.Data(), nameOutputSec[k].Data(), suffix.Data()));
        canvasEffectiveSecCorr->Print(Form("%s/Gamma_EffectiveSecCorr_%s.pdf",outputDir.Data(), nameOutputSec[k].Data()));
    }

    // **********************************************************************************************************************
    // ******************************** InvYield for combined inc gamma measurement *****************************************
    // **********************************************************************************************************************
    TCanvas* canvasInvYieldGamma          = new TCanvas("canvasInvYieldGamma","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasInvYieldGamma, 0.16, 0.02, 0.02, 0.08);
    canvasInvYieldGamma->SetLogx();
    canvasInvYieldGamma->SetLogy();


    TH1F * histo2DYieldGamma              = new TH1F("histo2DYieldGamma","histo2DYieldGamma",11000,doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(histo2DYieldGamma, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                              0.035,0.04, 0.035,0.04, 0.9, 1.7);
    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-11,10e1);
    histo2DYieldGamma->GetXaxis()->SetMoreLogLabels();
    histo2DYieldGamma->GetXaxis()->SetLabelOffset(-0.01);
    histo2DYieldGamma->Draw("copy");

        TGraphAsymmErrors* graphCombIncGammaStatPlot    = (TGraphAsymmErrors*)graphCombIncGammaStat->Clone("graphCombIncGammaStatPlot");
        ProduceGraphAsymmWithoutXErrors(graphCombIncGammaStatPlot);

        TLegend* legendYieldIncGamma           = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(4*textSizeLabelsRel*0.85),textSizeLabelsPixel);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes);
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");
        legendYieldIncGamma->AddEntry(graphCombIncGammaSys,"ALICE preliminary","pf");
        legendYieldIncGamma->Draw();

        DrawGammaSetMarkerTF1( fitHagGammaComb, 7, 2, colorCombpPb);
        legendYieldIncGamma->AddEntry(fitHagGammaComb,"Hagedorn fit","l");
        fitHagGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitHagGammaComb->Draw("same");
        DrawGammaSetMarkerTF1( fitTCMGammaComb, 9, 2, kGray+1);
        legendYieldIncGamma->AddEntry(fitTCMGammaComb,"TCM fit","l");
        fitTCMGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTCMGammaComb->Draw("same");
        DrawGammaSetMarkerTF1( fitTsallisGammaComb, 5, 2, colorCombpPbBox);
        legendYieldIncGamma->AddEntry(fitTsallisGammaComb,"Tsalis fit","l");
        fitTsallisGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTsallisGammaComb->Draw("same");

        TLatex *labelEnergyInvYieldPaperAll = new TLatex(0.20, 0.20+0.04*3, collisionSystempPbNSD.Data());
        TLatex *labelEnergyInvYieldPaperAll4 = new TLatex(0.20, 0.20+0.04*4, collisionSystempPbNSD.Data());
        SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        SetStyleTLatex( labelEnergyInvYieldPaperAll4, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelEnergyInvYieldPaperAll->Draw();
        TLatex *labelALICEInvYieldPaperAll  = new TLatex(0.20,0.20+0.04*3,textALICE.Data());
        SetStyleTLatex( labelALICEInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
//         labelALICEInvYieldPaperAll->Draw();
        TLatex *labelALICENormUnPaperAll    = new TLatex(0.20,0.20+0.05+0.04,"Norm. unc. 3.1%");
        SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICENormUnPaperAll->Draw();

        histo2DYieldGamma->Draw("same,axis");

    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma.pdf",outputDir.Data()));
    // **********************************************************************************************************************
    // ******************************** InvYield for individual inc gamma measurements **************************************
    // **********************************************************************************************************************
    histo2DYieldGamma->Draw("copy");

        TLegend* legendYieldIncGammaInd       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(4*textSizeLabelsRel*0.85),textSizeLabelsPixel);
        for (Int_t i = 0; i < 11; i++){
//             if (i == 4) continue;
            if (graphIndGammaIncSys[i]){
                DrawGammaSetMarkerTGraphAsym(graphIndGammaIncSys[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphIndGammaIncSys[i]->Draw("E2same");
                legendYieldIncGammaInd->AddEntry(graphIndGammaIncSys[i], Form("#gamma_{inc} %s", nameMeasGlobalLabel[i].Data()),"pf");
            }
            if (graphIndGammaIncStat[i]){
                DrawGammaSetMarkerTGraphAsym(graphIndGammaIncStat[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                graphIndGammaIncStat[i]->Draw("Epsame,x0");
                if (!graphIncGammaSysErr[i])legendYieldIncGammaInd->AddEntry(graphIndGammaIncStat[i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendYieldIncGammaInd->Draw();

        labelEnergyInvYieldPaperAll->Draw();
        labelALICEInvYieldPaperAll->Draw();
        labelALICENormUnPaperAll->Draw();

        histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_IndMeas.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_IndMeas.pdf",outputDir.Data()));

    // **********************************************************************************************************************
    // ******************************** InvYield for combined dir gamma measurement **************************************
    // **********************************************************************************************************************

    // prep plot graphs
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErrPlot = NULL;
    if (graphCombDirGammaSpectrumStatErr) graphCombDirGammaSpectrumStatErrPlot       = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr->Clone("graphCombDirGammaSpectrumStatErrPlot");
    if (graphCombDirGammaSpectrumStatErrPlot) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErrPlot);

    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-12,10e1);
    histo2DYieldGamma->Draw("copy");

//     TLegend* legendYieldIncGammaInd       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(3*textSizeLabelsRel*0.85),textSizeLabelsPixel);
//     legendYieldIncGammaInd->Draw();

    if (graphCombDirGammaSpectrumSystErr){
        DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE);
        graphCombDirGammaSpectrumSystErr->Draw("E2same");
    }
    if (graphCombDirGammaSpectrumStatErrPlot){
        DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlot, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb);
        graphCombDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
    }
    if (graphCombDirGammaSpectrumSumErrAr){
        DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrAr , 1, 3, colorCombpPb, colorCombpPb, 1.8, kTRUE);
        graphCombDirGammaSpectrumSumErrAr->Draw(">,same");
        PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr);
    }


    TLatex *labelEnergyDGInvYieldPaperAll = new TLatex(0.94, 0.965-0.04*1, collisionSystempPbNSD.Data());
    SetStyleTLatex( labelEnergyDGInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
    labelEnergyDGInvYieldPaperAll->Draw();
    TLatex *labelALICEDGInvYieldPaperAll  = new TLatex(0.94,0.965-0.04*2,textALICE.Data());
    SetStyleTLatex( labelALICEDGInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
    labelALICEDGInvYieldPaperAll->Draw();
    TLatex *labelALICEDGNormUnPaperAll    = new TLatex(0.94,0.965-(0.04*2+0.04*0.8),"Norm. unc. 3.1%");
    SetStyleTLatex( labelALICEDGNormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 31);
    labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma.pdf",outputDir.Data()));
    histo2DYieldGamma->Draw("copy");

    TLegend* legendYieldDirGamma        = GetAndSetLegend2(0.70, 0.84-(3*textSizeLabelsRel*0.85), 0.93, 0.84,textSizeLabelsPixel, 1, "", 43, 0.3);

        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys, markerStyleCombpPb+4, markerSizeCombpPb+0.2, colorCombpPbBox , colorCombpPbBox,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot, markerStyleCombpPb+4, markerSizeCombpPb+0.2, colorCombpPbBox , colorCombpPbBox, widthLinesBoxes);
        legendYieldDirGamma->AddEntry(graphCombIncGammaSys, "#gamma_{inc}","pf");
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");

        fitTCMGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTCMGammaComb->Draw("same");
        legendYieldDirGamma->AddEntry(fitTCMGammaComb,"TCM fit","l");

        if (graphCombDirGammaSpectrumSystErr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr->Draw("E2same");
            legendYieldDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr, "#gamma_{dir}","pf");
        }
        if (graphCombDirGammaSpectrumStatErrPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlot, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb);
            graphCombDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrAr , 1, 3, colorCombpPb, colorCombpPb, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr);
        }
        legendYieldDirGamma->Draw();

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma.pdf",outputDir.Data()));
    histo2DYieldGamma->Draw("copy");

        TLegend* legendYieldDirGamma2           = GetAndSetLegend2(0.70, 0.84-(2.4*textSizeLabelsRel*0.85), 0.93, 0.84,textSizeLabelsPixel, 1, "", 43, 0.3);
        TLegend* legendYieldDirGamma3           = GetAndSetLegend2(0.20, 0.11+(5.2*textSizeLabelsRel*0.85), 0.5, 0.11+(6.2*textSizeLabelsRel*0.85),textSizeLabelsPixel,1,
                                                                   "", 43, 0.23);
        TLegend* legendYieldDirGammaTheo        = GetAndSetLegend2(0.20, 0.11+(4.2*textSizeLabelsRel*0.85), 0.5, 0.11+(5.2*textSizeLabelsRel*0.85),textSizeLabelsPixel,1,
                                                                  "", 43, 0.23);
        TLegend* legendYieldDirGammaTheo2       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(4*textSizeLabelsRel*0.85),textSizeLabelsPixel,1, "#gamma_{dir} NLO pQCD:", 43, 0.23);
        legendYieldDirGammaTheo2->SetTextAlign(12);
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");
        legendYieldDirGamma2->AddEntry(graphCombIncGammaSys, "#gamma_{inc}","pf");

        fitTCMGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTCMGammaComb->Draw("same");
        legendYieldDirGamma2->AddEntry(fitTCMGammaComb,"TCM fit","l");


        if (graphCombDirGammaSpectrumSystErr){
            graphCombDirGammaSpectrumSystErr->Draw("E2same");
            legendYieldDirGamma3->AddEntry(graphCombDirGammaSpectrumSystErr, "#gamma_{dir}","pf");
        }
        legendYieldDirGamma2->Draw();
        legendYieldDirGamma3->Draw();
        if (graphTheoryNLOpPb) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryNLOpPb, 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
            graphTheoryNLOpPb->Draw("3,same");
            graphTheoryNLOpPb->SetLineWidth(widthLineNLO*1.5);
            legendYieldDirGammaTheo2->AddEntry(graphTheoryNLOpPb,"PDF: CT10, FF: GRV","l");
        }

        if (graphTheoryMCGillpPb) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryMCGillpPb, 0, 0, colorNLOMcGill, colorNLOMcGill, 0.6, kTRUE, colorNLOMcGill, kTRUE);
            graphTheoryMCGillpPb->Draw("3,same");
            graphTheoryMCGillpPb->SetLineWidth(widthLineNLO*1.5);
            legendYieldDirGammaTheo->AddEntry(graphTheoryMCGillpPb,"#gamma_{dir} Shen #it{et al.}","l");
        }
        if (graphTheoryPowhegnCTEQpPb) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryPowhegnCTEQpPb, 0, 0, colorNLONCTEQBand, colorNLONCTEQBand, 0.6, kTRUE, colorNLONCTEQBand, kTRUE);
            graphTheoryPowhegnCTEQpPb->Draw("3,same");
            legendYieldDirGammaTheo2->AddEntry(graphTheoryPowhegnCTEQpPb,"nPDF: nCTEQ15, FF: GRV","f");
        }

        if (graphTheoryPowhegEPPS16pPb) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryPowhegEPPS16pPb, 0, 0, colorNLOEPPSBand, colorNLOEPPSBand, 0.6, kTRUE, colorNLOEPPSBand, kTRUE);
            graphTheoryPowhegEPPS16pPb->Draw("3,same");
            legendYieldDirGammaTheo2->AddEntry(graphTheoryPowhegEPPS16pPb,"nPDF: EPPS16, FF: GRV","f");
        }
        legendYieldDirGammaTheo2->Draw();
        legendYieldDirGammaTheo->Draw();

        if (graphCombDirGammaSpectrumStatErrPlot){
            graphCombDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr){
            graphCombDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr);
        }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory.pdf",outputDir.Data()));

    //***************************************************************************************************
    //********************************* include isolated photon inv yield *******************************
    //***************************************************************************************************
    TH1F * histo2DYieldGamma2              = new TH1F("histo2DYieldGamma","histo2DYieldGamma",11000,doubleRatioXpp[0], 70);
    SetStyleHistoTH1ForGraphs(histo2DYieldGamma2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                              0.035,0.04, 0.035,0.04, 0.9, 1.7);
    histo2DYieldGamma2->GetYaxis()->SetRangeUser(7e-13,10e1);
    histo2DYieldGamma2->GetXaxis()->SetMoreLogLabels();
    histo2DYieldGamma2->GetXaxis()->SetLabelOffset(-0.01);
    histo2DYieldGamma2->Draw("copy");

        TLegend* legendYieldDirGamma4           = GetAndSetLegend2(0.20, 0.11+(5.2*textSizeLabelsRel*0.85), 0.5, 0.11+(7.2*textSizeLabelsRel*0.85),textSizeLabelsPixel,1,
                                                                   "", 43, 0.23);
        legendYieldDirGammaTheo2->SetTextAlign(12);
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");

        fitTCMGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTCMGammaComb->Draw("same");


        if (graphCombDirGammaSpectrumSystErr){
            graphCombDirGammaSpectrumSystErr->Draw("E2same");
            legendYieldDirGamma4->AddEntry(graphCombDirGammaSpectrumSystErr, "#gamma_{dir}","pf");
        }
        // prep plot graphs
        TGraphErrors* graphIsolatedGammaStatPlot = NULL;
        if (graphIsolatedGammaStat) graphIsolatedGammaStatPlot       = (TGraphErrors*)graphIsolatedGammaStat->Clone("graphCombDirGammaSpectrumStatErrPlot");
        if (graphIsolatedGammaStatPlot) ProduceGraphErrWithoutXErrors(graphIsolatedGammaStatPlot);

        if (graphIsolatedGammaSys){
            DrawGammaSetMarkerTGraphErr(graphIsolatedGammaSys, markerStyleCombpPb+4, markerSizeCombpPb+0.3, 807 , 807, widthLinesBoxes, kTRUE);
            graphIsolatedGammaSys->Draw("E2same");
            legendYieldDirGamma4->AddEntry(graphIsolatedGammaSys, "#gamma_{iso}","pf");
        }
        if (graphIsolatedGammaStatPlot){
            DrawGammaSetMarkerTGraphErr(graphIsolatedGammaStatPlot, markerStyleCombpPb+4, markerSizeCombpPb+0.3, 807 , 807);
            graphIsolatedGammaStatPlot->Draw("p,E1Z,same");
        }

        legendYieldDirGamma2->Draw();
        legendYieldDirGamma4->Draw();

        if (graphTheoryNLOpPb) {
            graphTheoryNLOpPb->Draw("3,same");
        }

        if (graphTheoryMCGillpPb) {
            graphTheoryMCGillpPb->Draw("3,same");
        }
        if (graphTheoryPowhegnCTEQpPb) {
            graphTheoryPowhegnCTEQpPb->Draw("3,same");
        }

        if (graphTheoryPowhegEPPS16pPb) {
            graphTheoryPowhegEPPS16pPb->Draw("3,same");
        }
        legendYieldDirGammaTheo2->Draw();
        legendYieldDirGammaTheo->Draw();

        if (graphCombDirGammaSpectrumStatErrPlot){
            graphCombDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr){
            graphCombDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr);
        }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma2->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_IsoGamma_Theory.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_IsoGamma_Theory.pdf",outputDir.Data()));


        // **********************************************************************************************************************
    // ******************************** InvYield for combined dir gamma measurement **************************************
    // **********************************************************************************************************************

    // prep plot graphs
    TGraphAsymmErrors* graphCombPi0FitDirGammaSpectrumStatErrPlot = NULL;
    if (graphCombPi0FitDirGammaSpectrumStatErr) graphCombPi0FitDirGammaSpectrumStatErrPlot       = (TGraphAsymmErrors*)graphCombPi0FitDirGammaSpectrumStatErr->Clone("graphCombPi0FitDirGammaSpectrumStatErrPlot");
    if (graphCombPi0FitDirGammaSpectrumStatErrPlot) ProduceGraphAsymmWithoutXErrors(graphCombPi0FitDirGammaSpectrumStatErrPlot);

    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-12,10e1);
    histo2DYieldGamma->Draw("copy");

//     TLegend* legendYieldIncGammaInd       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(3*textSizeLabelsRel*0.85),textSizeLabelsPixel);
//     legendYieldIncGammaInd->Draw();

    if (graphCombPi0FitDirGammaSpectrumSystErr){
        DrawGammaSetMarkerTGraphAsym(graphCombPi0FitDirGammaSpectrumSystErr, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE);
        graphCombPi0FitDirGammaSpectrumSystErr->Draw("E2same");
    }
    if (graphCombPi0FitDirGammaSpectrumStatErrPlot){
        DrawGammaSetMarkerTGraphAsym(graphCombPi0FitDirGammaSpectrumStatErrPlot, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb);
        graphCombPi0FitDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
    }
    if (graphCombPi0FitDirGammaSpectrumSumErrAr){
        DrawGammaSetMarkerTGraphAsym(graphCombPi0FitDirGammaSpectrumSumErrAr , 1, 3, colorCombpPb, colorCombpPb, 1.8, kTRUE);
        graphCombPi0FitDirGammaSpectrumSumErrAr->Draw(">,same");
        PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombPi0FitDirGammaSpectrumSumErrAr);
    }


    labelEnergyDGInvYieldPaperAll->Draw();
    labelALICEDGInvYieldPaperAll->Draw();
    labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/Pi0Fitted/InvYield_DirGamma.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/Pi0Fitted/InvYield_DirGamma.pdf",outputDir.Data()));
    histo2DYieldGamma->Draw("copy");


        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys, markerStyleCombpPb+4, markerSizeCombpPb+0.2, colorCombpPbBox , colorCombpPbBox,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot, markerStyleCombpPb+4, markerSizeCombpPb+0.2, colorCombpPbBox , colorCombpPbBox, widthLinesBoxes);
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");

        fitTCMGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTCMGammaComb->Draw("same");

        if (graphCombPi0FitDirGammaSpectrumSystErr){
            DrawGammaSetMarkerTGraphAsym(graphCombPi0FitDirGammaSpectrumSystErr, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE);
            graphCombPi0FitDirGammaSpectrumSystErr->Draw("E2same");
        }
        if (graphCombPi0FitDirGammaSpectrumStatErrPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombPi0FitDirGammaSpectrumStatErrPlot, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb);
            graphCombPi0FitDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombPi0FitDirGammaSpectrumSumErrAr){
            DrawGammaSetMarkerTGraphAsym(graphCombPi0FitDirGammaSpectrumSumErrAr , 1, 3, colorCombpPb, colorCombpPb, 1.8, kTRUE);
            graphCombPi0FitDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombPi0FitDirGammaSpectrumSumErrAr);
        }
        legendYieldDirGamma->Draw();

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/Pi0Fitted/InvYield_DirGamma_IncGamma.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/Pi0Fitted/InvYield_DirGamma_IncGamma.pdf",outputDir.Data()));
    histo2DYieldGamma->Draw("copy");

        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");

        fitTCMGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTCMGammaComb->Draw("same");

        if (graphCombPi0FitDirGammaSpectrumSystErr){
            graphCombPi0FitDirGammaSpectrumSystErr->Draw("E2same");
        }
        legendYieldDirGamma2->Draw();
        legendYieldDirGamma3->Draw();
        if (graphTheoryNLOpPb) {
            graphTheoryNLOpPb->Draw("3,same");
        }

        if (graphTheoryMCGillpPb) {
            graphTheoryMCGillpPb->Draw("3,same");
        }
        if (graphTheoryPowhegnCTEQpPb) {
            graphTheoryPowhegnCTEQpPb->Draw("3,same");
        }

        if (graphTheoryPowhegEPPS16pPb) {
            graphTheoryPowhegEPPS16pPb->Draw("3,same");
        }
        legendYieldDirGammaTheo2->Draw();
        legendYieldDirGammaTheo->Draw();

        if (graphCombPi0FitDirGammaSpectrumStatErrPlot){
            graphCombPi0FitDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombPi0FitDirGammaSpectrumSumErrAr){
            graphCombPi0FitDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombPi0FitDirGammaSpectrumSumErrAr);
        }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/Pi0Fitted/InvYield_DirGamma_IncGamma_Theory.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/Pi0Fitted/InvYield_DirGamma_IncGamma_Theory.pdf",outputDir.Data()));

    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-12,0.8e0);
    histo2DYieldGamma->Draw("copy");

        if (graphTheoryNLOpPb) {
            graphTheoryNLOpPb->Draw("3,same");
        }

        if (graphTheoryMCGillpPb) {
            graphTheoryMCGillpPb->Draw("3,same");
        }
        if (graphTheoryPowhegnCTEQpPb) {
            graphTheoryPowhegnCTEQpPb->Draw("3,same");
        }

        if (graphTheoryPowhegEPPS16pPb) {
            graphTheoryPowhegEPPS16pPb->Draw("3,same");
        }
        legendYieldDirGammaTheo2->Draw();
        legendYieldDirGammaTheo->Draw();

        labelEnergyDGInvYieldPaperAll->Draw();
//         labelALICEDGInvYieldPaperAll->Draw();
//         labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_Theory.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_Theory.pdf",outputDir.Data()));

    //***************************************************************************************************
    //********************************* include isolated photon inv yield *******************************
    //***************************************************************************************************
    histo2DYieldGamma2->Draw("copy");

        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");

        fitTCMGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTCMGammaComb->Draw("same");


        if (graphCombPi0FitDirGammaSpectrumSystErr){
            graphCombPi0FitDirGammaSpectrumSystErr->Draw("E2same");
        }
        // prep plot graphs
        if (graphIsolatedGammaSys){
            graphIsolatedGammaSys->Draw("E2same");
        }
        if (graphIsolatedGammaStatPlot){
            graphIsolatedGammaStatPlot->Draw("p,E1Z,same");
        }

        legendYieldDirGamma2->Draw();
        legendYieldDirGamma4->Draw();

        if (graphTheoryNLOpPb) {
            graphTheoryNLOpPb->Draw("3,same");
        }
        if (graphTheoryNLOpPbPrompt) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryNLOpPbPrompt, 0, 0, 807, 807, 0.2, kTRUE, 807);
            graphTheoryNLOpPbPrompt->Draw("3,same");
        }
        if (graphTheoryMCGillpPb) {
            graphTheoryMCGillpPb->Draw("3,same");
        }
        if (graphTheoryPowhegnCTEQpPb) {
            graphTheoryPowhegnCTEQpPb->Draw("3,same");
        }

        if (graphTheoryPowhegEPPS16pPb) {
            graphTheoryPowhegEPPS16pPb->Draw("3,same");
        }
        legendYieldDirGammaTheo2->Draw();
        legendYieldDirGammaTheo->Draw();

        if (graphCombPi0FitDirGammaSpectrumStatErrPlot){
            graphCombPi0FitDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombPi0FitDirGammaSpectrumSumErrAr){
            graphCombPi0FitDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombPi0FitDirGammaSpectrumSumErrAr);
        }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma2->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/Pi0Fitted/InvYield_DirGamma_IncGamma_IsoGamma_Theory.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/Pi0Fitted/InvYield_DirGamma_IncGamma_IsoGamma_Theory.pdf",outputDir.Data()));


    cout << "************ rel stat err iso *****************************************" << endl;
    TGraphErrors* graphRelStatIsolatedGamma    = (TGraphErrors*) graphIsolatedGammaStat->Clone("relStatIso");
    while (graphRelStatIsolatedGamma->GetY()[0] == 0) graphRelStatIsolatedGamma->RemovePoint(0);
    while (graphRelStatIsolatedGamma->GetY()[graphRelStatIsolatedGamma->GetN()-1] == 0) graphRelStatIsolatedGamma->RemovePoint(graphRelStatIsolatedGamma->GetN()-1);
    graphRelStatIsolatedGamma    = CalculateRelErrGraph( graphRelStatIsolatedGamma, "relativeStatErrorIso");
    graphRelStatIsolatedGamma->Print();
    cout << "************ rel sys err iso *****************************************" << endl;
    TGraphErrors* graphRelSysIsolatedGamma    = (TGraphErrors*) graphIsolatedGammaSys->Clone("relSysIso");
    while (graphRelSysIsolatedGamma->GetY()[0] == 0) graphRelSysIsolatedGamma->RemovePoint(0);
    while (graphRelSysIsolatedGamma->GetY()[graphRelSysIsolatedGamma->GetN()-1] == 0) graphRelSysIsolatedGamma->RemovePoint(graphRelSysIsolatedGamma->GetN()-1);
    graphRelSysIsolatedGamma    = CalculateRelErrGraph( graphRelSysIsolatedGamma, "relativeSysErrorIso");
    graphRelSysIsolatedGamma->Print();

    // **********************************************************************************************************************
    // ******************************** GammaRawYields for individual inc gamma measurements ********************************
    // **********************************************************************************************************************
    histo2DYieldGamma->GetYaxis()->SetTitle("#frac{#it{N}_{raw,#gamma}}{#it{N}_{ev}}");
    histo2DYieldGamma->Draw("copy");
    for (Int_t i = 0; i < 11; i++){
        if (histoRawGamma[i]){
            cout << nameMeasGlobalLabel[i].Data() << endl;
            TGraphErrors* dummy         = new TGraphErrors(histoRawGamma[i]);
            dummy->Print();
            DrawGammaSetMarker(histoRawGamma[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
            histoRawGamma[i]->Draw("p,same,e0,X0");
            delete dummy;
        }
    }
    legendYieldIncGammaInd->Draw();

    labelEnergyInvYieldPaperAll->Draw();
    labelALICEInvYieldPaperAll->Draw();
    labelALICENormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/RawYield_IncGamma_IndMeas.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/RawYield_IncGamma_IndMeas.pdf",outputDir.Data()));

    // ****************************************************************************************************************
    // ************************** Store final results including corr factors in 1 file ********************************
    // ****************************************************************************************************************
    TString nameOutputCommonFile    = Form("CombinedGammaResultPPb5TeV_%s.root", dateForOutput.Data());

    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    fCombResults.mkdir("Gamma_pPb5TeV");
    TDirectoryFile* directoryGamma = (TDirectoryFile*)fCombResults.Get("Gamma_pPb5TeV");
    fCombResults.cd("Gamma_pPb5TeV");

    // writing main results
    if (graphCombIncGammaStat) graphCombIncGammaStat->Write("graphInvYieldIncGammaStatErr");
    if (graphCombIncGammaSys) graphCombIncGammaSys->Write("graphInvYieldIncGammaSysErr");
    if (graphCombIncGammaTot) graphCombIncGammaTot->Write("graphInvYieldIncGammaTotErr");
    if (graphCombIncGammaStatUnshi) graphCombIncGammaStatUnshi->Write("graphInvYieldIncGammaStatErr_Unshifted");
    if (graphCombIncGammaSysUnshi) graphCombIncGammaSysUnshi->Write("graphInvYieldIncGammaSysErr_Unshifted");
    if (graphCombIncGammaTotUnshi) graphCombIncGammaTotUnshi->Write("graphInvYieldIncGammaTotErr_Unshifted");

    if (graphCombDRNonFitStat) graphCombDRNonFitStat->Write("graphRGammaCombNonFitStatErr");
    if (graphCombDRNonFitSys) graphCombDRNonFitSys->Write("graphRGammaCombNonFitSysErr");
    if (graphCombDRNonFitTot) graphCombDRNonFitTot->Write("graphRGammaCombNonFitTotErr");
    if (graphCombDirGammaSpectrumSystErr) graphCombDirGammaSpectrumSystErr->Write("graphInvYieldNonFitDirGammaSysErr");
    if (graphCombDirGammaSpectrumStatErr) graphCombDirGammaSpectrumStatErr->Write("graphInvYieldNonFitDirGammaStatErr");
    if (graphCombDirGammaSpectrumSumErrAr) graphCombDirGammaSpectrumSumErrAr->Write("graphInvYieldNonFitDirGammaSumErrAr");
    if (graphCombDRStat) graphCombDRStat->Write("graphRGammaCombStatErr");
    if (graphCombDRSys) graphCombDRSys->Write("graphRGammaCombSysErr");
    if (graphCombDRTot) graphCombDRTot->Write("graphRGammaCombTotErr");
    if (graphCombPi0FitDirGammaSpectrumSystErr) graphCombPi0FitDirGammaSpectrumSystErr->Write("graphInvYieldDirGammaSysErr");
    if (graphCombPi0FitDirGammaSpectrumStatErr) graphCombPi0FitDirGammaSpectrumStatErr->Write("graphInvYieldDirGammaStatErr");
    if (graphCombPi0FitDirGammaSpectrumSumErrAr) graphCombPi0FitDirGammaSpectrumSumErrAr->Write("graphInvYieldDirGammaSumErrAr");

    for (Int_t i = 0; i < 11; i++){
        if (graphIndGammaIncSys[i]) graphIndGammaIncSys[i]->Write(Form("graphInvYieldIncGamma%sSysErr",nameMeasGlobal[i].Data()));
        if (graphIndGammaIncStat[i]) graphIndGammaIncStat[i]->Write(Form("graphInvYieldIncGamma%sStatErr",nameMeasGlobal[i].Data()));
        if (histoIncGammaStatErr[i]) histoIncGammaStatErr[i]->Write(Form("histoInvYieldIncGamma%sStatErr_Unshifted",nameMeasGlobal[i].Data()));
        if (histoDRPi0FitStatErr[i]) histoDRPi0FitStatErr[i]->Write(Form("histoRGamma%sStatErr",nameMeasGlobal[i].Data()));
        if (graphDRPi0FitSysErr[i]) graphDRPi0FitSysErr[i]->Write(Form("graphRGamma%sSysErr",nameMeasGlobal[i].Data()));
        if (histoDRNonFitStatErr[i]) histoDRNonFitStatErr[i]->Write(Form("histoNonFitRGamma%sStatErr",nameMeasGlobal[i].Data()));
        if (graphDRNonFitSysErr[i]) graphDRNonFitSysErr[i]->Write(Form("graphNonfitRGamma%sSysErr",nameMeasGlobal[i].Data()));
    }
    directoryGamma->mkdir("Supporting");
    directoryGamma->cd("Supporting");
    // Writing full correction factors
    for (Int_t i = 0; i < 11; i++){
        if (histoPileupCorr[i]) histoPileupCorr[i]->Write(Form("histoPileupCorr%s",nameMeasGlobal[i].Data()));
        if (histoConvProb[i]) histoConvProb[i]->Write(Form("histoConvProb%s",nameMeasGlobal[i].Data()));
        if (histoEffi[i]) histoEffi[i]->Write(Form("histoEffiEffective%s",nameMeasGlobal[i].Data()));
        if (histoEffiMCPt[i]) histoEffiMCPt[i]->Write(Form("histoEffiWithoutResol%s",nameMeasGlobal[i].Data()));
        if (histoPurity[i]) histoPurity[i]->Write(Form("histoPurity%s",nameMeasGlobal[i].Data()));
        for (Int_t k = 0; k<4; k++){
            if (histoEffSecCorr[k][i]) histoEffSecCorr[k][i]->Write(Form("histoEffectiveSecCorrFrom%s_%s",nameOutputSec[i].Data(), nameMeasGlobal[i].Data()));
        }
    }

    fCombResults.Close();
}