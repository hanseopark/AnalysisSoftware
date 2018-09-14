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



void CombineGammaResultspPb_V2(    TString inputFileNamePCM        = "",
                                TString inputFileNamePHOS       = "",
                                TString inputFileNameEMC        = "",
                                TString inputFileNamePCMEMC     = "",
                                TString suffix                  = "eps",
                                Bool_t enablepValueCalc         = kFALSE,
                                Bool_t isThesis                 = kFALSE
                            ){

    Int_t nMeasSetDR[5]               = {4,4,4,4,4};
    Int_t nMeasSetIncGamma[5]         = {4,4,4,4,4};
    Bool_t combineCent[5]             = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
    TString fileNameCorrelations[5]    = {"/media/nschmidt/local/ANALYSIS_RESULTS/pPb/5.023TeV/dirGammaPHOS/pPb5TeV_CF_20180913.root","/media/nschmidt/local/ANALYSIS_RESULTS/pPb/5.023TeV/dirGammaPHOS/pPb5TeV_CF_20180913.root","/media/nschmidt/local/ANALYSIS_RESULTS/pPb/5.023TeV/dirGammaPHOS/pPb5TeV_CF_20180913.root","/media/nschmidt/local/ANALYSIS_RESULTS/pPb/5.023TeV/dirGammaPHOS/pPb5TeV_CF_20180913.root","/media/nschmidt/local/ANALYSIS_RESULTS/pPb/5.023TeV/dirGammaPHOS/pPb5TeV_CF_20180913.root"};
    //*******************************************************************************************************************************************
    //*********************************************************** Set main style choices ********************************************************
    //*******************************************************************************************************************************************
    StyleSettingsThesis();
    SetPlotStyle();

    //*******************************************************************************************************************************************
    //********************************************* Create output directory and copy input files there ******************************************
    //*******************************************************************************************************************************************
    TString dateForOutput                                       = ReturnDateStringForOutput();
    TString outputDir                                           = Form("%s/%s/CombineGammaResultspPb_V2",suffix.Data(),dateForOutput.Data());
    TString fileNameTheorypPb                                   = "ExternalInputpPb/Theory/TheoryCompilationPPb.root";
    TString fileNameIsolatedGamma                               = "ExternalInputpPb/EMCAL/isolated_photon_invariant_yield_20180417.root";


    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("mkdir -p %s/RelUncAndWeights",outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOSGammapPb.root", inputFileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMGammapPb.root", inputFileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEMCGammapPb.root", inputFileNamePCMEMC.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCGammapPb.root", inputFileNameEMC.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/TheorypPb.root", fileNameTheorypPb.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/CorrelationFactors_0020.root", fileNameCorrelations[0].Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/CorrelationFactors_2040.root", fileNameCorrelations[1].Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/CorrelationFactors_4060.root", fileNameCorrelations[2].Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/CorrelationFactors_60100.root", fileNameCorrelations[3].Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/CorrelationFactors_NSD.root", fileNameCorrelations[4].Data(), outputDir.Data()));

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
    doubleRatioX[0]     = 0.33;     doubleRatioX[1]     = 49;
    doubleRatioXpp[0]   = 0.33;     doubleRatioXpp[1]   = 49;

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
    Color_t  colorMult[5];
    Style_t  markerStyleDet[11];
    Style_t  markerStyleDetMC[11];
    Style_t  markerStyleMult[5];
    Style_t  markerStyleMultOpen[5] = {24,27,24,28,24};
    Size_t   markerSizeDet[11];
    Size_t   markerSizeDetMC[11];
    Size_t   markerSizeMult[5];

    TString multbins[5]              = {"0-20%","20-40%","40-60%","60-100%",""};
    TString multbinsCorrF[5]              = {"0020/","2040/","4060/","60100/",""};
    TString multbinsPHOS[5]              = {"00-20","20-40","40-60","60-100","00-100"};
    for (Int_t i = 0; i < 5; i++){
        colorMult[i]                                 = GetColorDefaultColor("pPb_5.023TeV", "", multbins[i].Data());
        markerStyleMult[i]                           = GetDefaultMarkerStyle("pPb_5.023TeV", "", multbins[i].Data());
        markerSizeMult[i]                            = GetDefaultMarkerSize("pPb_5.023TeV", "", multbins[i].Data());
    }
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

    TString textALICE                               = "ALICE work-in-progress";
    TString textALICEPerf                           = "ALICE performance";
    if (isThesis) textALICE                         = "ALICE this thesis";

    cout << "Setting Gamma binning" << endl;
    Double_t xPtLimitsGamma[5][100]                    = {0};
    Int_t maxNBinsGammaAbs[5]                       = {0};
    Int_t maxNBinsGamma[5];
    for(Int_t ncent = 0; ncent < 5; ncent++){
        if(ncent<4)
            maxNBinsGamma[ncent]                             = GetBinning( xPtLimitsGamma[ncent], maxNBinsGammaAbs[ncent], "Gamma", "pPb_5.023TeV", 24, kFALSE, multbins[ncent].Data() );
        else
            maxNBinsGamma[ncent]                             = GetBinning( xPtLimitsGamma[ncent], maxNBinsGammaAbs[ncent], "Gamma", "pPb_5.023TeV", 23 );
        cout << "combination binning used for " << multbins[ncent].Data() << " is:" << endl;
        for (Int_t i = 0; i< maxNBinsGamma[ncent]; i++){
            cout << i << ": "<< xPtLimitsGamma[ncent][i] <<" - " << xPtLimitsGamma[ncent][i+1]<< ", " <<endl;
        }
    }

    TH1D* histoDRNonFitStatErr[5][11]                  = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TGraphAsymmErrors*  graphDRNonFitSysErr[5][11]     = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TH1D* histoIncGammaRatioStatErr[5][11]             = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TGraphAsymmErrors* graphIncGammaRatioSysErr[5][11] = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TH1D* histoIncGammaStatErr[5][11]                  = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TGraphAsymmErrors* graphIncGammaSysErr[5][11]      = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TH1D* histoRawGamma[5][11]                         = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TH1D* histoConvProb[5][11]                         = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TH1D* histoEffi[5][11]                             = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TH1D* histoEffiMCPt[5][11]                         = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TH1D* histoResolCorr[5][11]                        = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TH1D* histoPurity[5][11]                           = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    // TH1D* histoEffSecCorr[4][5][11]                    = { {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
    //                                                     {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
    //                                                     {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
    //                                                     {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL} };
    TH1D* histoTotalCorrFactor[5][11]                  = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    TH1D* histoPileupCorr[5][11]                       = {  {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                            {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL}};
    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from pPb PCM file ******************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMGammapPb                          = new TFile( inputFileNamePCM.Data());
    //________________________________________________ Load PCM pPb 5.023TeV _________________________________________________________________________
    if (!filePCMGammapPb->IsZombie()){
        for(Int_t i = 0; i<5; i++){
            TDirectory* directoryPCMGammapPb                = (TDirectory*)filePCMGammapPb->Get(Form("Gamma_pPb5TeV%s",multbins[i].Data()));
            if(directoryPCMGammapPb)
                cout << "found PCM " << multbins[i].Data() << "input folder... loading spectra:" << endl;
                histoDRNonFitStatErr[i][0]                         = (TH1D*) directoryPCMGammapPb->Get("DoubleRatioStatError");
                graphDRNonFitSysErr[i][0]                          = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("DoubleRatioSystError");
                histoIncGammaRatioStatErr[i][0]                    = (TH1D*) directoryPCMGammapPb->Get("IncRatioStatError");
                graphIncGammaRatioSysErr[i][0]                     = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncRatioSystError");
                histoConvProb[i][0]                                = (TH1D*) directoryPCMGammapPb->Get("GammaConversionProbability");
                histoEffi[i][0]                                    = (TH1D*) directoryPCMGammapPb->Get("GammaRecoEfficiency");
                histoEffiMCPt[i][0]                                = (TH1D*) directoryPCMGammapPb->Get("GammaRecoEfficiencyMCPt");
                histoResolCorr[i][0]                               = (TH1D*) directoryPCMGammapPb->Get("GammaResolCorr");
                histoPurity[i][0]                                  = (TH1D*) directoryPCMGammapPb->Get("GammaTruePurity");
                histoIncGammaStatErr[i][0]                         = (TH1D*) directoryPCMGammapPb->Get("IncGammaStatError");
                graphIncGammaSysErr[i][0]                          = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncGammaSystError");

                histoPileupCorr[i][0]                              = (TH1D*) directoryPCMGammapPb->Get("PileUpCorrectionFactor");
                histoRawGamma[i][0]                                = (TH1D*) directoryPCMGammapPb->Get("GammaRawYields");
                if(i==4){
                    if (histoIncGammaStatErr[i][0])
                        histoIncGammaStatErr[i][0]->Scale(scalingToNSD);
                    if (graphIncGammaSysErr[i][0])
                        graphIncGammaSysErr[i][0]                      = ScaleGraph(graphIncGammaSysErr[i][0],scalingToNSD);
                }
                // histoEffSecCorr[i][0][i][0]                           = (TH1D*) directoryPCMGammapPb->Get("GammaEffectiveSecondaryCorr_K0s");
                // histoEffSecCorr[1][i][0]                           = (TH1D*) directoryPCMGammapPb->Get("GammaEffectiveSecondaryCorr_K0l");
                // histoEffSecCorr[2][i][0]                           = (TH1D*) directoryPCMGammapPb->Get("GammaEffectiveSecondaryCorr_Lambda");
                // histoEffSecCorr[3][i][0]                           = (TH1D*) directoryPCMGammapPb->Get("GammaEffectiveSecondaryCorr_Rest");
                histoTotalCorrFactor[i][0]                         = (TH1D*) histoEffi[i][0]->Clone("GammaTotalCorrFactorPCM");
                histoTotalCorrFactor[i][0]->Sumw2();
                histoTotalCorrFactor[i][0]->Divide(histoPurity[i][0]);
                histoTotalCorrFactor[i][0]->Multiply(histoConvProb[i][0]);
                cout << "Successfully loaded all inputs" << endl;
        }
    }

    //*******************************************************************************************************************************************
    //*********************************************** Load PCMEMC histograms from pPb PCM file ******************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMEMCGammapPb                       = new TFile( inputFileNamePCMEMC.Data());
    //________________________________________________ Load PCM-EMC pPb 5.023TeV _________________________________________________________________________
    if (!filePCMEMCGammapPb->IsZombie()){
        for(Int_t i = 0; i<5; i++){
            TDirectory* directoryPCMEMCGammapPb                = (TDirectory*)filePCMEMCGammapPb->Get(Form("Gamma_pPb5TeV%s",multbins[i].Data()));
            if(directoryPCMEMCGammapPb)
                cout << "found PCMEMC " << multbins[i].Data() << "input folder... loading spectra:" << endl;
                histoDRNonFitStatErr[i][4]                         = (TH1D*) directoryPCMEMCGammapPb->Get("DoubleRatioStatError");
                graphDRNonFitSysErr[i][4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("DoubleRatioSystError");
                histoIncGammaRatioStatErr[i][4]                    = (TH1D*) directoryPCMEMCGammapPb->Get("IncRatioStatError");
                graphIncGammaRatioSysErr[i][4]                     = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("IncRatioSystError");
                histoConvProb[i][4]                                = (TH1D*) directoryPCMEMCGammapPb->Get("GammaConversionProbability");
                histoEffi[i][4]                                    = (TH1D*) directoryPCMEMCGammapPb->Get("GammaRecoEfficiency");
                histoEffiMCPt[i][4]                                = (TH1D*) directoryPCMEMCGammapPb->Get("GammaRecoEfficiencyMCPt");
                histoResolCorr[i][4]                               = (TH1D*) directoryPCMEMCGammapPb->Get("GammaResolCorr");
                histoPurity[i][4]                                  = (TH1D*) directoryPCMEMCGammapPb->Get("GammaTruePurity");
                histoIncGammaStatErr[i][4]                         = (TH1D*) directoryPCMEMCGammapPb->Get("IncGammaStatError");
                graphIncGammaSysErr[i][4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("IncGammaSystError");

                histoPileupCorr[i][4]                              = (TH1D*) directoryPCMEMCGammapPb->Get("PileUpCorrectionFactor");
                histoRawGamma[i][4]                                = (TH1D*) directoryPCMEMCGammapPb->Get("GammaRawYields");
                if(i==4){
                    if (histoIncGammaStatErr[i][4])
                        histoIncGammaStatErr[i][4]->Scale(scalingToNSD);
                    if (graphIncGammaSysErr[i][4])
                        graphIncGammaSysErr[i][4]                      = ScaleGraph(graphIncGammaSysErr[i][4],scalingToNSD);
                }
                // histoEffSecCorr[i][4][i][4]                           = (TH1D*) directoryPCMEMCGammapPb->Get("GammaEffectiveSecondaryCorr_K0s");
                // histoEffSecCorr[1][i][4]                           = (TH1D*) directoryPCMEMCGammapPb->Get("GammaEffectiveSecondaryCorr_K0l");
                // histoEffSecCorr[2][i][4]                           = (TH1D*) directoryPCMEMCGammapPb->Get("GammaEffectiveSecondaryCorr_Lambda");
                // histoEffSecCorr[3][i][4]                           = (TH1D*) directoryPCMEMCGammapPb->Get("GammaEffectiveSecondaryCorr_Rest");
                histoTotalCorrFactor[i][4]                         = (TH1D*) histoEffi[i][4]->Clone("GammaTotalCorrFactorPCMEMC");
                histoTotalCorrFactor[i][4]->Sumw2();
                histoTotalCorrFactor[i][4]->Divide(histoPurity[i][4]);
                histoTotalCorrFactor[i][4]->Multiply(histoConvProb[i][4]);
                cout << "Successfully loaded all inputs" << endl;
        }
    }

    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from pPb EMC file ******************************************************
    //*******************************************************************************************************************************************
    TFile* fileEMCGammapPb                          = new TFile( inputFileNameEMC.Data());
    //________________________________________________ Load EMC pPb 5.023TeV _________________________________________________________________________
    if (!fileEMCGammapPb->IsZombie()){
        for(Int_t i = 0; i<5; i++){
            TDirectory* directoryEMCGammapPb                = (TDirectory*)fileEMCGammapPb->Get(Form("Gamma_pPb5TeV%s",multbins[i].Data()));
                cout << "found EMC " << multbins[i].Data() << "input folder... loading spectra:" << endl;
                histoDRNonFitStatErr[i][2]                         = (TH1D*) directoryEMCGammapPb->Get("DoubleRatioStatError");
                graphDRNonFitSysErr[i][2]                          = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("DoubleRatioSystError");
                histoIncGammaRatioStatErr[i][2]                    = (TH1D*) directoryEMCGammapPb->Get("IncRatioStatError");
                graphIncGammaRatioSysErr[i][2]                     = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("IncRatioSystError");
                histoEffi[i][2]                                    = (TH1D*) directoryEMCGammapPb->Get("GammaRecoEfficiency");
                histoEffiMCPt[i][2]                                = (TH1D*) directoryEMCGammapPb->Get("GammaRecoEfficiencyMCPt");
                histoResolCorr[i][2]                               = (TH1D*) directoryEMCGammapPb->Get("GammaResolCorr");
                histoPurity[i][2]                                  = (TH1D*) directoryEMCGammapPb->Get("GammaTruePurity");
                histoIncGammaStatErr[i][2]                         = (TH1D*) directoryEMCGammapPb->Get("IncGammaStatError");
                graphIncGammaSysErr[i][2]                          = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("IncGammaSystError");
                histoRawGamma[i][2]                                = (TH1D*) directoryEMCGammapPb->Get("GammaRawYields");
                // histoEffSecCorr[0][i][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_K0s");
                // histoEffSecCorr[1][i][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_K0l");
                // histoEffSecCorr[i][2][i][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_Lambda");
                // histoEffSecCorr[3][i][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_Rest");
                if(i==4){
                    if (histoIncGammaStatErr[i][2])
                        histoIncGammaStatErr[i][2]->Scale(scalingToNSD);
                    if (graphIncGammaSysErr[i][2])
                        graphIncGammaSysErr[i][2]                      = ScaleGraph(graphIncGammaSysErr[i][2],scalingToNSD);

                    while (graphDRNonFitSysErr[i][2]->GetX()[0] < 2.05) graphDRNonFitSysErr[i][2]->RemovePoint(0);
                    while (graphIncGammaRatioSysErr[i][2]->GetX()[0] < 2.05) graphIncGammaRatioSysErr[i][2]->RemovePoint(0);
                    while (graphIncGammaSysErr[i][2]->GetX()[0] < 2.05) graphIncGammaSysErr[i][2]->RemovePoint(0);

                    Int_t j = 1;
                        while (histoEffi[i][2]->GetBinCenter(j) < 2.05){
                        histoEffi[i][2]->SetBinContent(j,1e-10);
                        histoEffiMCPt[i][2]->SetBinContent(j,1e-10);
                        histoResolCorr[i][2]->SetBinContent(j,1e10);
                        histoPurity[i][2]->SetBinContent(j,1e10);
                        histoIncGammaStatErr[i][2]->SetBinContent(j,0);
                        histoDRNonFitStatErr[i][2]->SetBinContent(j,0);
                        histoIncGammaRatioStatErr[i][2]->SetBinContent(j,0);
                        // histoEffSecCorr[0][i][2]->SetBinContent(j,1e10);
                        // histoEffSecCorr[1][i][2]->SetBinContent(j,1e10);
                        // histoEffSecCorr[2][i][2]->SetBinContent(j,1e10);
                        // histoEffSecCorr[3][i][2]->SetBinContent(j,1e10);
                        j++;
                        }
                }

                histoTotalCorrFactor[i][2]                         = (TH1D*) histoEffi[i][2]->Clone("GammaTotalCorrFactorEMC");
                histoTotalCorrFactor[i][2]->Sumw2();
                histoTotalCorrFactor[i][2]->Divide(histoPurity[i][2]);
                histoTotalCorrFactor[i][2]->Scale(100/360.*0.66/0.9);
                histoTotalCorrFactor[i][2]->GetYaxis()->SetRangeUser(0.01,2);
                cout << "Successfully loaded all inputs" << endl;
        }
    }
    //*******************************************************************************************************************************************
    //*********************************************** Load PHOS histograms from PHOS file *******************************************************
    //*******************************************************************************************************************************************
    TFile* filePHOSGamma                            = new TFile( inputFileNamePHOS.Data());
    if (!filePHOSGamma->IsZombie()){
        for(Int_t i = 0; i<5; i++){
            TDirectory* directoryPHOSGammapPb                   = (TDirectory*)filePHOSGamma->Get(Form("PHOS_pPb_2760_Centrality_%s",multbinsPHOS[i].Data()));
                cout << "found PHOS " << multbins[i].Data() << "input folder... loading spectra:" << endl;
                histoDRNonFitStatErr[i][1]                         = (TH1D*) directoryPHOSGammapPb->Get(Form("hPHOS_DoubleRatio_pPb_cen%s_Stat",multbinsPHOS[i].Data()));
                TH1D* histoPHOSDRNonFitSysErrpPb                = (TH1D*) directoryPHOSGammapPb->Get(Form("hPHOS_DoubleRatio_pPb_cen%s_Syst",multbinsPHOS[i].Data()));
                graphDRNonFitSysErr[i][1]                          = new TGraphAsymmErrors(histoPHOSDRNonFitSysErrpPb);
                TF1* fitEffiTimesAccPHOSINT7                    = (TF1*) directoryPHOSGammapPb->Get(Form("hPHOS_gammaIncl_pPb_cen%s_Efficiency",multbinsPHOS[i].Data()));
                TF1* fitEffiTimesAccPHOSPHI7                    = (TF1*) directoryPHOSGammapPb->Get(Form("hPHOS_gammaIncl_pPb_PHI7_cen%s_Efficiency",multbinsPHOS[i].Data()));
                histoPurity[i][1]                                  = (TH1D*) histoPHOSDRNonFitSysErrpPb->Clone("GammaTruePurity");
                histoEffi[i][1]                                    = (TH1D*) histoPHOSDRNonFitSysErrpPb->Clone("GammaRecoEfficiency");

                histoIncGammaStatErr[i][1]                         = (TH1D*) directoryPHOSGammapPb->Get(Form("hPHOS_gammaIncl_pPb_cen%s_Stat",multbinsPHOS[i].Data()));
                TH1D* histoPHOSIncGammaPi0FitSysErrpPb          = (TH1D*) directoryPHOSGammapPb->Get(Form("hPHOS_gammaIncl_pPb_cen%s_Syst",multbinsPHOS[i].Data()));
                graphIncGammaSysErr[i][1]                          = new TGraphAsymmErrors(histoPHOSIncGammaPi0FitSysErrpPb);

                histoTotalCorrFactor[i][1]                         = (TH1D*) histoEffi[i][1]->Clone("GammaTotalCorrFactorEMC");
                histoTotalCorrFactor[i][1]->Sumw2();
                histoTotalCorrFactor[i][1]->Scale(60/360.*0.1/0.9);
                histoTotalCorrFactor[i][1]->Divide(histoPurity[i][1]);
                histoTotalCorrFactor[i][1]->GetYaxis()->SetRangeUser(0.01,2);
                cout << "Successfully loaded all inputs" << endl;
        }
    }
    multbins[4]="NSD";

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



        // Declaration & calculation of combined spectrum
                                              //pcm, phos, emc, pcmphos, pcmemc,   dalitz, phosdalitz,emcdalit
        Int_t offSetsGamma[5][11]              = {  { 0,  6,  0,  0,  0,       0,  0,  0,  0,  0,      0},
                                                    { 0,  6,  0,  0,  0,       0,  0,  0,  0,  0,      0},
                                                    { 0,  6,  0,  0,  0,       0,  0,  0,  0,  0,      0},
                                                    { 0,  6,  0,  0,  0,       0,  0,  0,  0,  0,      0},
                                                    { 2, 14,  2,  0,  0,       0,  0,  0,  0,  0,      0}};
        Int_t offSetsGammaSys[5][11]           = {  { 2,  6,  8,  0,  6,       0,  0,  0,  0,  0,      0},
                                                    { 2,  6,  8,  0,  6,       0,  0,  0,  0,  0,      0},
                                                    { 2,  6,  8,  0,  6,       0,  0,  0,  0,  0,      0},
                                                    { 2,  6,  8,  0,  6,       0,  0,  0,  0,  0,      0},
                                                    { 3, 14, 14,  0,  8,       0,  0,  0,  0,  0,      0}};
        TGraphAsymmErrors* statErrorRelCollectionDRNonFit[5][11];
        TGraphAsymmErrors* sysErrorRelCollectionDRNonFit[5][11];
        TGraph* graphWeightsDR[5][11];

        TString fileNameDROutputWeighting[5];
        TGraphAsymmErrors* graphCombDRStat[5];
        TGraphAsymmErrors* graphCombDRSys[5];
        TGraphAsymmErrors* graphCombDRTot[5];


        Double_t xValuesDRRead[5][50];
        Double_t weightsDRRead[5][11][50];
        Int_t availableDRMeas[5][11]      = {   { -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1},
                                                { -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1},
                                                { -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1},
                                                { -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1},
                                                { -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1}};
        Int_t nPtBinsDRRead[5]            = {0,0,0,0,0};

        Int_t textSizeLabelsPixel                 = 900*0.04;



                                              //pcm, phos, emc, pcmphos, pcmemc,   dalitz, phosdalitz,emcdalit
        Int_t offSetsIncGamma[5][11]           = {  { 0,  6,  0,  0,  0,       0,  0,  0,  0,  0,      0},
                                                    { 0,  6,  0,  0,  0,       0,  0,  0,  0,  0,      0},
                                                    { 0,  6,  0,  0,  0,       0,  0,  0,  0,  0,      0},
                                                    { 0,  6,  0,  0,  0,       0,  0,  0,  0,  0,      0},
                                                    { 2, 14,  2,  0,  0,       0,  0,  0,  0,  0,      0}};
        Int_t offSetsIncGammaSys[5][11]        = {  { 2,  6,  8,  0,  6,       0,  0,  0,  0,  0,      0},
                                                    { 2,  6,  8,  0,  6,       0,  0,  0,  0,  0,      0},
                                                    { 2,  6,  8,  0,  6,       0,  0,  0,  0,  0,      0},
                                                    { 2,  6,  8,  0,  6,       0,  0,  0,  0,  0,      0},
                                                    { 3, 14, 14,  0,  8,       0,  0,  0,  0,  0,      0}};
        TGraphAsymmErrors* statErrorGraphCollectionIncGamma[5][11];
        TGraphAsymmErrors* statErrorRelCollectionIncGamma[5][11];
        TGraphAsymmErrors* sysErrorRelCollectionIncGamma[5][11];
        TGraph* graphWeightsIncGamma[5][11];
        TString fileNameIncGammaOutputWeighting[5];
        TGraphAsymmErrors* graphCombIncGammaStat[5];
        TGraphAsymmErrors* graphCombIncGammaSys[5];
        TGraphAsymmErrors* graphCombIncGammaTot[5];

        Double_t xValuesIncGammaRead[5][50];
        Double_t weightsIncGammaRead[5][11][50];
        Int_t availableIncGammaMeas[5][11]      = {   { -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1},
                                                { -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1},
                                                { -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1},
                                                { -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1},
                                                { -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1}};
        Int_t nPtBinsIncGammaRead[5]            = {0,0,0,0,0};

        TF1* fitHagGammaComb[5];
        TF1* fitTCMGammaComb[5];
        TF1* fitTsallisGammaComb[5];
        TGraphAsymmErrors* graphCombIncGammaTotUnshi[5];
        TGraphAsymmErrors* graphCombIncGammaStatUnshi[5];
        TGraphAsymmErrors* graphCombIncGammaSysUnshi[5];

        TGraphAsymmErrors* graphIndGammaIncStatUnshi[5][11];
        TGraphAsymmErrors* graphIndGammaIncSysUnshi[5][11];
        TGraphAsymmErrors* graphIndGammaIncStat[5][11];
        TGraphAsymmErrors* graphIndGammaIncSys[5][11];
        TGraphAsymmErrors* graphIndGammaIncStat_yShifted[5][11];
        TGraphAsymmErrors* graphIndGammaIncSys_yShifted[5][11];


        TGraphAsymmErrors* graphRatioGammaIndCombFitStat[5][11];
        TGraphAsymmErrors* graphRatioGammaIndCombFitSys[5][11];

        TGraphAsymmErrors* graphRatioGammaCombCombFitTot[5];
        TGraphAsymmErrors* graphRatioGammaCombCombFitStat[5];
        TGraphAsymmErrors* graphRatioGammaCombCombFitSys[5];

        TGraphAsymmErrors *graphCombDirGammaSpectrumSystErr[5];
        TGraphAsymmErrors *graphCombDirGammaSpectrumStatErr[5];
        TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr[5];
        TGraphAsymmErrors *graphCombDirGammaSpectrumSumErrConfi[5];
        TGraphAsymmErrors *graphCombDirGammaSpectrumSumErrAr[5];

        TGraphAsymmErrors* graphCombDRStatPlot[5];
        TGraphAsymmErrors* graphCombIncGammaStatPlot[5];
        TGraphAsymmErrors* graphCombDirGammaSpectrumStatErrPlot[5];
        //*******************************************************************************************************************************************
        //*********************************************** Combining Rgamma ratios  ******************************************************************
        //*******************************************************************************************************************************************
        // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
        // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    for(Int_t ncent = 0; ncent<5; ncent++){
        if(!combineCent[ncent])
            continue;

        for (Int_t i = 0; i < 11; i++){
            statErrorRelCollectionDRNonFit[ncent][i]  = NULL;
            if (histoDRNonFitStatErr[ncent][i]){
                statErrorRelCollectionDRNonFit[ncent][i]   = new TGraphAsymmErrors(histoDRNonFitStatErr[ncent][i]);
                while (statErrorRelCollectionDRNonFit[ncent][i]->GetY()[0] == 0) statErrorRelCollectionDRNonFit[ncent][i]->RemovePoint(0);
                while (statErrorRelCollectionDRNonFit[ncent][i]->GetY()[statErrorRelCollectionDRNonFit[ncent][i]->GetN()-1] == 0)
                    statErrorRelCollectionDRNonFit[ncent][i]->RemovePoint(statErrorRelCollectionDRNonFit[ncent][i]->GetN()-1);
                statErrorRelCollectionDRNonFit[ncent][i]    = CalculateRelErrUpAsymmGraph( statErrorRelCollectionDRNonFit[ncent][i], Form("relativeStatErrorDRNonFit_%s", nameMeasGlobal[i].Data()));
            }
        }


        for (Int_t i = 0; i < 11; i++){
            sysErrorRelCollectionDRNonFit[ncent][i]   = NULL;
            cout << i << endl;
            if (graphDRNonFitSysErr[ncent][i]){
                sysErrorRelCollectionDRNonFit[ncent][i]     = (TGraphAsymmErrors*)graphDRNonFitSysErr[ncent][i]->Clone(Form("relativeSysErrorDRNonFit_%s", nameMeasGlobal[i].Data()));
                sysErrorRelCollectionDRNonFit[ncent][i]->Print();
                while (sysErrorRelCollectionDRNonFit[ncent][i]->GetY()[0] == 0) sysErrorRelCollectionDRNonFit[ncent][i]->RemovePoint(0);
                while (sysErrorRelCollectionDRNonFit[ncent][i]->GetY()[sysErrorRelCollectionDRNonFit[ncent][i]->GetN()-1] == 0) sysErrorRelCollectionDRNonFit[ncent][i]->RemovePoint(sysErrorRelCollectionDRNonFit[ncent][i]->GetN()-1);
                sysErrorRelCollectionDRNonFit[ncent][i]     = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionDRNonFit[ncent][i], Form("relativeSysErrorDRNonFit_%s", nameMeasGlobal[i].Data()));
                cout << "after" << endl;
                sysErrorRelCollectionDRNonFit[ncent][i]->Print();
            }
        }

        // **********************************************************************************************************************
        // ******************************************* Assuming maximal correlation *********************************************
        // **********************************************************************************************************************


        for (Int_t i = 0; i< 11; i++){
            graphWeightsDR[ncent][i]                   = NULL;
        }

        // Declaration & calculation of combined spectrum
        fileNameDROutputWeighting[ncent]        = Form("%s/DR_WeightingMethod_%s.dat",outputDir.Data(),multbinsPHOS[ncent].Data());
        graphCombDRStat[ncent]      = NULL;
        graphCombDRSys[ncent]        = NULL;
        graphCombDRTot[ncent]        = CombinePtPointsSpectraFullCorrMat(    histoDRNonFitStatErr[ncent] ,    graphDRNonFitSysErr[ncent] ,
                                                                                        xPtLimitsGamma[ncent] , maxNBinsGamma[ncent] ,
                                                                                        offSetsGamma[ncent] , offSetsGammaSys[ncent] ,
                                                                                        graphCombDRStat[ncent] , graphCombDRSys[ncent] ,
                                                                                        fileNameDROutputWeighting[ncent] , "pPb_5.023TeV", "RGamma", kTRUE,
                                                                                        NULL, fileNameCorrelations[ncent]  , multbinsCorrF[ncent].Data());
                                                                                        // NULL, ""  );


        if (graphCombDRTot[ncent]  == NULL) {
            cout << "Aborting: something went wrong during the combination of the new spectra for cent " << multbins[ncent].Data() << endl;
            return;
        }
        while (graphCombDRStat[ncent]->GetX()[0] < 0.4){
            graphCombDRStat[ncent]->RemovePoint(0);
        }
        while (graphCombDRTot[ncent]->GetX()[0] < 0.4){
            graphCombDRTot[ncent]->RemovePoint(0);
        }
        while (graphCombDRSys[ncent]->GetX()[0] < 0.4){
            graphCombDRSys[ncent]->RemovePoint(0);
        }
        graphCombDRTot[ncent]->Print();
// return;

    // Reading weights from output file for plotting
    ifstream fileWeightsDRRead;
    fileWeightsDRRead.open(fileNameDROutputWeighting[ncent],ios_base::in);
    cout << "reading" << fileNameDROutputWeighting[ncent] << endl;
    while(!fileWeightsDRRead.eof() && nPtBinsDRRead[ncent] < 50){
        TString garbage             = "";
        if (nPtBinsDRRead[ncent] == 0){
            fileWeightsDRRead >> garbage ;//>> availableDRMeas[0] >> availableDRMeas[1] >> availableDRMeas[2] >> availableDRMeas[3];
            for (Int_t i = 0; i < nMeasSetDR[ncent]; i++){
                fileWeightsDRRead >> availableDRMeas[ncent][i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableDRMeas[ncent][i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsDRRead >> xValuesDRRead[ncent][nPtBinsDRRead[ncent]-1];
            for (Int_t i = 0; i < nMeasSetDR[ncent]; i++){
                fileWeightsDRRead >> weightsDRRead[ncent][availableDRMeas[ncent][i]][nPtBinsDRRead[ncent]-1] ;
            }
            cout << "read: "<<  nPtBinsDRRead[ncent] << "\t"<< xValuesDRRead[ncent][nPtBinsDRRead[ncent]-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetDR[ncent]; i++){
                cout << weightsDRRead[ncent][availableDRMeas[ncent][i]][nPtBinsDRRead[ncent]-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsDRRead[ncent]++;
    }
    nPtBinsDRRead[ncent]                  = nPtBinsDRRead[ncent]-2 ;
    fileWeightsDRRead.close();

    for (Int_t i = 0; i < nMeasSetDR[ncent]; i++){
        graphWeightsDR[ncent][availableDRMeas[ncent][i]]                        = new TGraph(nPtBinsDRRead[ncent],xValuesDRRead[ncent],weightsDRRead[ncent][availableDRMeas[ncent][i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsDRRead[ncent]; n++){
            if (graphWeightsDR[ncent][availableDRMeas[ncent][i]]->GetY()[bin] == 0) graphWeightsDR[ncent][availableDRMeas[ncent][i]]->RemovePoint(bin);
            else bin++;
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plotting weights method only EMC *****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 900*0.04;

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
        TLegend* legendWeightsDR   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetDR[ncent]+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
        for (Int_t i = 0; i < nMeasSetDR[ncent]; i++){
            DrawGammaSetMarkerTGraph(graphWeightsDR[ncent][availableDRMeas[ncent][i]], markerStyleDet[availableDRMeas[ncent][i]], markerSizeDet[availableDRMeas[ncent][i]], colorDet[availableDRMeas[ncent][i]] , colorDet[availableDRMeas[ncent][i]]);
            graphWeightsDR[ncent][availableDRMeas[ncent][i]]->Draw("p,same,z");
            legendWeightsDR->AddEntry(graphWeightsDR[ncent][availableDRMeas[ncent][i]],nameMeasGlobalLabel[availableDRMeas[ncent][i]],"p");
        }
        legendWeightsDR->Draw();

        TLatex *labelWeightsEnergy      = new TLatex(0.95,0.20,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelWeightsEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelWeightsEnergy->Draw();
        TLatex *labelWeightsDR         = new TLatex(0.95,0.15,"R_{#gamma}");
        SetStyleTLatex( labelWeightsDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelWeightsDR->Draw();

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/RelUncAndWeights/DR_Weights_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasWeights->SaveAs(Form("%s/RelUncAndWeights/DR_Weights_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));
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

        TLegend* legendRelSysErr2       = GetAndSetLegend2(0.62, 0.92-(0.04*(nMeasSetDR[ncent]+1)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
        for (Int_t i = 0; i < nMeasSetDR[ncent]; i++){
            cout << "sys\t" << nameMeasGlobalLabel[availableDRMeas[ncent][i]] << endl;
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionDRNonFit[ncent][availableDRMeas[ncent][i]], markerStyleDet[availableDRMeas[ncent][i]], markerSizeDet[availableDRMeas[ncent][i]], colorDet[availableDRMeas[ncent][i]],
                                    colorDet[availableDRMeas[ncent][i]]);
            sysErrorRelCollectionDRNonFit[ncent][availableDRMeas[ncent][i]]->Draw("p,same,z");
            sysErrorRelCollectionDRNonFit[ncent][availableDRMeas[ncent][i]]->Print();
            legendRelSysErr2->AddEntry(sysErrorRelCollectionDRNonFit[ncent][availableDRMeas[ncent][i]],nameMeasGlobalLabel[availableDRMeas[ncent][i]],"p");
        }
        legendRelSysErr2->Draw();

        TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel, 4, 1, 43);
        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrDR       = new TLatex(0.15,0.85,"R_{#gamma}");
        SetStyleTLatex( labelRelSysErrDR, textSizeLabelsPixel, 4, 1, 43);
        labelRelSysErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/RelUncAndWeights/DR_RelSysErr_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/RelUncAndWeights/DR_RelSysErr_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));

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
        TLegend* legendRelStatErr2       = GetAndSetLegend2(0.14, 0.92-(0.04*(nMeasSetDR[ncent]+1)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
        for (Int_t i = 0; i < nMeasSetDR[ncent]; i++){
            DrawGammaSetMarkerTGraph(statErrorRelCollectionDRNonFit[ncent][availableDRMeas[ncent][i]], markerStyleDet[availableDRMeas[ncent][i]], markerSizeDet[availableDRMeas[ncent][i]], colorDet[availableDRMeas[ncent][i]],
                                    colorDet[availableDRMeas[ncent][i]]);
            statErrorRelCollectionDRNonFit[ncent][availableDRMeas[ncent][i]]->Draw("p,same,z");
            legendRelStatErr2->AddEntry(statErrorRelCollectionDRNonFit[ncent][availableDRMeas[ncent][i]],nameMeasGlobalLabel[availableDRMeas[ncent][i]],"p");
        }
        legendRelStatErr2->Draw();

        TLatex *labelRelStatErrEnergy   = new TLatex(0.95,0.89,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelRelStatErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
        SetStyleTLatex( labelRelStatErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelStatErrDR->Draw();

    canvasRelStatErr->SaveAs(Form("%s/RelUncAndWeights/DR_RelStatErr_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasRelStatErr->SaveAs(Form("%s/RelUncAndWeights/DR_RelStatErr_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));


    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TGraphAsymmErrors* graphCombDRRelStat       = CalculateRelErrUpAsymmGraph( graphCombDRStat[ncent], "relativeStatErrorDR");
    TGraphAsymmErrors* graphCombDRRelSys        = CalculateRelErrUpAsymmGraph( graphCombDRSys[ncent], "relativeSysErrorDR");
    TGraphAsymmErrors* graphCombDRRelTot        = CalculateRelErrUpAsymmGraph( graphCombDRTot[ncent], "relativeTotalErrorDR");

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

        TLatex *labelRelTotErrEnergy   = new TLatex(0.95,0.89,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelRelTotErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
        SetStyleTLatex( labelRelTotErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/RelUncAndWeights/DR_Reldecomp_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));





    //*******************************************************************************************************************************************
    //*********************************************** Combining Inc gamma spectra  **************************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};

    for (Int_t i = 0; i< 11; i++){
        statErrorGraphCollectionIncGamma[ncent][i]   = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (histoIncGammaStatErr[ncent][i]){
            statErrorGraphCollectionIncGamma[ncent][i]    = new TGraphAsymmErrors(histoIncGammaStatErr[ncent][i]);
            while (statErrorGraphCollectionIncGamma[ncent][i]->GetY()[0] == 0) statErrorGraphCollectionIncGamma[ncent][i]->RemovePoint(0);
            while (statErrorGraphCollectionIncGamma[ncent][i]->GetY()[statErrorGraphCollectionIncGamma[ncent][i]->GetN()-1] == 0) statErrorGraphCollectionIncGamma[ncent][i]->RemovePoint(statErrorGraphCollectionIncGamma[ncent][i]->GetN()-1);
            statErrorGraphCollectionIncGamma[ncent][i]->SetName(Form("statErrorIncGamma_%s_%s", nameMeasGlobal[i].Data(),multbinsPHOS[ncent].Data()));
        }
    }

    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionIncGamma[ncent][i]        = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (histoIncGammaStatErr[ncent][i]){
            statErrorRelCollectionIncGamma[ncent][i]    = new TGraphAsymmErrors(histoIncGammaStatErr[ncent][i]);
            while (statErrorRelCollectionIncGamma[ncent][i]->GetY()[0] == 0) statErrorRelCollectionIncGamma[ncent][i]->RemovePoint(0);
            while (statErrorRelCollectionIncGamma[ncent][i]->GetY()[statErrorRelCollectionIncGamma[ncent][i]->GetN()-1] == 0) statErrorRelCollectionIncGamma[ncent][i]->RemovePoint(statErrorRelCollectionIncGamma[ncent][i]->GetN()-1);
            statErrorRelCollectionIncGamma[ncent][i]    = CalculateRelErrUpAsymmGraph( statErrorRelCollectionIncGamma[ncent][i], Form("relativeStatErrorIncGamma_%s_%s", nameMeasGlobal[i].Data(),multbinsPHOS[ncent].Data()));
        }
    }


    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollectionIncGamma[ncent][i]         = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        cout << i << endl;
        if (graphIncGammaSysErr[ncent][i]){
            sysErrorRelCollectionIncGamma[ncent][i]     = (TGraphAsymmErrors*)graphIncGammaSysErr[ncent][i]->Clone(Form("relativeSysErrorIncGamma_%s_%s", nameMeasGlobal[i].Data(),multbinsPHOS[ncent].Data()));
            sysErrorRelCollectionIncGamma[ncent][i]->Print();
            while (sysErrorRelCollectionIncGamma[ncent][i]->GetY()[0] == 0) sysErrorRelCollectionIncGamma[ncent][i]->RemovePoint(0);
            while (sysErrorRelCollectionIncGamma[ncent][i]->GetY()[sysErrorRelCollectionIncGamma[ncent][i]->GetN()-1] == 0) sysErrorRelCollectionIncGamma[ncent][i]->RemovePoint(sysErrorRelCollectionIncGamma[ncent][i]->GetN()-1);
            sysErrorRelCollectionIncGamma[ncent][i]     = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionIncGamma[ncent][i], Form("relativeSysErrorIncGamma_%s_%s", nameMeasGlobal[i].Data(),multbinsPHOS[ncent].Data()));
            cout << "after" << endl;
            sysErrorRelCollectionIncGamma[ncent][i]->Print();
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    for (Int_t i = 0; i< 11; i++){
        graphWeightsIncGamma[ncent][i]                   = NULL;
    }

    // Declaration & calculation of combined spectrum
    fileNameIncGammaOutputWeighting[ncent]       = Form("%s/IncGamma_WeightingMethod_%s.dat",outputDir.Data(),multbins[ncent].Data());
    graphCombIncGammaStat[ncent]      = NULL;
    graphCombIncGammaSys[ncent]       = NULL;
    graphCombIncGammaTot[ncent]       = CombinePtPointsSpectraFullCorrMat(      histoIncGammaStatErr[ncent],    graphIncGammaSysErr[ncent],
                                                                                            xPtLimitsGamma[ncent], maxNBinsGamma[ncent],
                                                                                            offSetsIncGamma[ncent], offSetsIncGammaSys[ncent],
                                                                                            graphCombIncGammaStat[ncent], graphCombIncGammaSys[ncent],
                                                                                            fileNameIncGammaOutputWeighting[ncent], "pPb_5.023TeV", "GammaInc", kTRUE,
                                                                                            // NULL, "" );
                                                                                           NULL, fileNameCorrelations[ncent]  , multbinsCorrF[ncent].Data());


    if (graphCombIncGammaTot[ncent] == NULL) {
        cout << "Aborting: something went wrong during the combination of the new inc gamma spectra for cent " << multbins[ncent].Data() << endl;
        return;
    }
    while (graphCombIncGammaStat[ncent]->GetX()[0] < 0.4){
        graphCombIncGammaStat[ncent]->RemovePoint(0);
    }
    while (graphCombIncGammaTot[ncent]->GetX()[0] < 0.4){
        graphCombIncGammaTot[ncent]->RemovePoint(0);
    }
    while (graphCombIncGammaSys[ncent]->GetX()[0] < 0.4){
        graphCombIncGammaSys[ncent]->RemovePoint(0);
    }
    graphCombIncGammaTot[ncent]->Print();


    // Reading weights from output file for plotting
    ifstream fileWeightsIncGammaRead;
    fileWeightsIncGammaRead.open(fileNameIncGammaOutputWeighting[ncent],ios_base::in);
    cout << "reading" << fileNameIncGammaOutputWeighting[ncent] << endl;
    while(!fileWeightsIncGammaRead.eof() && nPtBinsIncGammaRead[ncent] < 50){
        TString garbage             = "";
        if (nPtBinsIncGammaRead[ncent] == 0){
            fileWeightsIncGammaRead >> garbage ;//>> availableIncGammaMeas[0] >> availableIncGammaMeas[1] >> availableIncGammaMeas[2] >> availableIncGammaMeas[3];
            for (Int_t i = 0; i < nMeasSetIncGamma[ncent]; i++){
                fileWeightsIncGammaRead >> availableIncGammaMeas[ncent][i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableIncGammaMeas[ncent][i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsIncGammaRead >> xValuesIncGammaRead[ncent][nPtBinsIncGammaRead[ncent]-1];
            for (Int_t i = 0; i < nMeasSetIncGamma[ncent]; i++){
                fileWeightsIncGammaRead >> weightsIncGammaRead[ncent][availableIncGammaMeas[ncent][i]][nPtBinsIncGammaRead[ncent]-1] ;
            }
            cout << "read: "<<  nPtBinsIncGammaRead[ncent] << "\t"<< xValuesIncGammaRead[ncent][nPtBinsIncGammaRead[ncent]-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetIncGamma[ncent]; i++){
                cout << weightsIncGammaRead[ncent][availableIncGammaMeas[ncent][i]][nPtBinsIncGammaRead[ncent]-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsIncGammaRead[ncent]++;
    }
    nPtBinsIncGammaRead[ncent]                  = nPtBinsIncGammaRead[ncent]-2 ;
    fileWeightsIncGammaRead.close();

    for (Int_t i = 0; i < nMeasSetIncGamma[ncent]; i++){
        graphWeightsIncGamma[ncent][availableIncGammaMeas[ncent][i]]                        = new TGraph(nPtBinsIncGammaRead[ncent],xValuesIncGammaRead[ncent],weightsIncGammaRead[ncent][availableIncGammaMeas[ncent][i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsIncGammaRead[ncent]; n++){
            if (graphWeightsIncGamma[ncent][availableIncGammaMeas[ncent][i]]->GetY()[bin] == 0) graphWeightsIncGamma[ncent][availableIncGammaMeas[ncent][i]]->RemovePoint(bin);
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

    TLegend* legendWeightsIncGamma   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetIncGamma[ncent]+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetIncGamma[ncent]; i++){
        DrawGammaSetMarkerTGraph(graphWeightsIncGamma[ncent][availableIncGammaMeas[ncent][i]], markerStyleDet[availableIncGammaMeas[ncent][i]], markerSizeDet[availableIncGammaMeas[ncent][i]], colorDet[availableIncGammaMeas[ncent][i]] , colorDet[availableIncGammaMeas[ncent][i]]);
        graphWeightsIncGamma[ncent][availableIncGammaMeas[ncent][i]]->Draw("p,same,z");
        legendWeightsIncGamma->AddEntry(graphWeightsIncGamma[ncent][availableIncGammaMeas[ncent][i]],nameMeasGlobalLabel[availableIncGammaMeas[ncent][i]],"p");
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

    canvasWeights->SaveAs(Form("%s/RelUncAndWeights/IncGamma_Weights_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasWeights->SaveAs(Form("%s/RelUncAndWeights/IncGamma_Weights_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelSysErr->cd();

    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,17.5);
    histo2DRelSysErr->Draw("copy");

    TLegend* legendRelSysErrIncGamma       = GetAndSetLegend2(0.62, 0.92-(0.04*(nMeasSetIncGamma[ncent]+1)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetIncGamma[ncent]; i++){
        cout << "sys\t" << nameMeasGlobalLabel[availableIncGammaMeas[ncent][i]] << endl;
        DrawGammaSetMarkerTGraph(sysErrorRelCollectionIncGamma[ncent][availableIncGammaMeas[ncent][i]], markerStyleDet[availableIncGammaMeas[ncent][i]], markerSizeDet[availableIncGammaMeas[ncent][i]], colorDet[availableIncGammaMeas[ncent][i]],
                                 colorDet[availableIncGammaMeas[ncent][i]]);
        sysErrorRelCollectionIncGamma[ncent][availableIncGammaMeas[ncent][i]]->Draw("p,same,z");
        sysErrorRelCollectionIncGamma[ncent][availableIncGammaMeas[ncent][i]]->Print();
        legendRelSysErrIncGamma->AddEntry(sysErrorRelCollectionIncGamma[ncent][availableIncGammaMeas[ncent][i]],nameMeasGlobalLabel[availableIncGammaMeas[ncent][i]],"p");
    }
    legendRelSysErrIncGamma->Draw();

    labelRelSysErrEnergy->Draw();
    TLatex *labelRelSysErrIncGamma       = new TLatex(0.15,0.85,"#gamma_{inc}");
    SetStyleTLatex( labelRelSysErrIncGamma, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrIncGamma->Draw();

    canvasRelSysErr->SaveAs(Form("%s/RelUncAndWeights/IncGamma_RelSysErr_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/RelUncAndWeights/IncGamma_RelSysErr_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    canvasRelStatErr->cd();

    histo2DRelStatErr->GetYaxis()->SetRangeUser(-0.2,30.0);
    histo2DRelStatErr->Draw("copy");
    TLegend* legendRelStatErrIncGamma       = GetAndSetLegend2(0.14, 0.92-(0.04*(nMeasSetIncGamma[ncent]+1)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetIncGamma[ncent]; i++){
        DrawGammaSetMarkerTGraph(statErrorRelCollectionIncGamma[ncent][availableIncGammaMeas[ncent][i]], markerStyleDet[availableIncGammaMeas[ncent][i]], markerSizeDet[availableIncGammaMeas[ncent][i]], colorDet[availableIncGammaMeas[ncent][i]],
                                 colorDet[availableIncGammaMeas[ncent][i]]);
        statErrorRelCollectionIncGamma[ncent][availableIncGammaMeas[ncent][i]]->Draw("p,same,z");
        legendRelStatErrIncGamma->AddEntry(statErrorRelCollectionIncGamma[ncent][availableIncGammaMeas[ncent][i]],nameMeasGlobalLabel[availableIncGammaMeas[ncent][i]],"p");
    }
    legendRelStatErrIncGamma->Draw();

    labelRelStatErrEnergy->Draw();
    TLatex *labelRelStatErrIncGamma      = new TLatex(0.95,0.85,"#gamma_{inc}");
    SetStyleTLatex( labelRelStatErrIncGamma, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrIncGamma->Draw();

    canvasRelStatErr->SaveAs(Form("%s/RelUncAndWeights/IncGamma_RelStatErr_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasRelStatErr->SaveAs(Form("%s/RelUncAndWeights/IncGamma_RelStatErr_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TGraphAsymmErrors* graphCombIncGammaRelStat     = CalculateRelErrUpAsymmGraph( graphCombIncGammaStat[ncent], Form("relativeStatErrorDR_%s",multbinsPHOS[ncent].Data()));
    TGraphAsymmErrors* graphCombIncGammaRelSys      = CalculateRelErrUpAsymmGraph( graphCombIncGammaSys[ncent], Form("relativeSysErrorDR_%s",multbinsPHOS[ncent].Data()));
    TGraphAsymmErrors* graphCombIncGammaRelTot      = CalculateRelErrUpAsymmGraph( graphCombIncGammaTot[ncent], Form("relativeTotalErrorDR_%s",multbinsPHOS[ncent].Data()));

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

    canvasRelSysErr->SaveAs(Form("%s/RelUncAndWeights/IncGamma_Reldecomp_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/RelUncAndWeights/IncGamma_Reldecomp_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));

    //*******************************************************************************************************************************************
    //************************************************* Fitting gamma spectrum ******************************************************************
    //*******************************************************************************************************************************************

    Double_t paramGraphoHag[5]      = {10,-0.22,-0.01,0.57,6.28};
    if(ncent==2){paramGraphoHag[0]=65;paramGraphoHag[1]=-0.17;paramGraphoHag[2]=0.008;paramGraphoHag[3]=1.1;paramGraphoHag[4]=7.1;}
    fitHagGammaComb[ncent]            = FitObject("oHag",Form("fitHagGammaComb%s",multbinsPHOS[ncent].Data()),"Gamma",graphCombIncGammaTot[ncent],ncent==2 ? graphCombIncGammaTot[ncent]->GetX()[3] : graphCombIncGammaTot[ncent]->GetX()[2],
                                                graphCombIncGammaTot[ncent]->GetX()[graphCombIncGammaTot[ncent]->GetN()],paramGraphoHag,"QNRME+");
    TString forOutput               = WriteParameterToFile(fitHagGammaComb[ncent]);
    cout << forOutput.Data() << endl;
    Double_t paramTCM[5] = {graphCombIncGammaTot[ncent]->GetY()[1],0.1,graphCombIncGammaTot[ncent]->GetY()[4],0.6,3};
    fitTCMGammaComb[ncent]            = FitObject("tcm", Form("fitTCMGammaComb%s",multbinsPHOS[ncent].Data()),"Gamma", graphCombIncGammaTot[ncent], graphCombIncGammaTot[ncent]->GetX()[0],
                                                graphCombIncGammaTot[ncent]->GetX()[graphCombIncGammaTot[ncent]->GetN()], paramTCM,"QNRME+");
    forOutput               = WriteParameterToFile(fitTCMGammaComb[ncent]);
    cout << forOutput.Data() << endl;

    fitTsallisGammaComb[ncent]        = FitObject("l",Form("fitTsallisGammaComb%s",multbinsPHOS[ncent].Data()),"Gamma",graphCombIncGammaTot[ncent],graphCombIncGammaTot[ncent]->GetX()[0],
                                                graphCombIncGammaTot[ncent]->GetX()[graphCombIncGammaTot[ncent]->GetN()],NULL,"QNRME+");
    forOutput                       = WriteParameterToFile(fitTsallisGammaComb[ncent]);
    cout << forOutput.Data() << endl;


    // **********************************************************************************************************************
    // ************************************* Calculating bin shifted spectra & fitting **************************************
    // **********************************************************************************************************************

    // Cloning spectra
    graphCombIncGammaTotUnshi[ncent]         = (TGraphAsymmErrors*)graphCombIncGammaTot[ncent]->Clone(Form("GammaUnshifted%s",multbinsPHOS[ncent].Data()));
    graphCombIncGammaStatUnshi[ncent]        = (TGraphAsymmErrors*)graphCombIncGammaStat[ncent]->Clone(Form("GammaUnshiftedStat%s",multbinsPHOS[ncent].Data()));
    graphCombIncGammaSysUnshi[ncent]         = (TGraphAsymmErrors*)graphCombIncGammaSys[ncent]->Clone(Form("GammaUnshiftedSys%s",multbinsPHOS[ncent].Data()));


    for (Int_t i = 0; i< 11; i++){
        graphIndGammaIncStatUnshi[ncent][i]         = NULL;
        graphIndGammaIncStat[ncent][i]              = NULL;
        graphIndGammaIncStat_yShifted[ncent][i]     = NULL;
        if (statErrorGraphCollectionIncGamma[ncent][i]){
            graphIndGammaIncStatUnshi[ncent][i]                 = (TGraphAsymmErrors*)statErrorGraphCollectionIncGamma[ncent][i]->Clone(Form("GammaUnshiftedStat%s_%s",nameMeasGlobalLabel[i].Data(),multbinsPHOS[ncent].Data()));
            graphIndGammaIncStat[ncent][i]                      = (TGraphAsymmErrors*)statErrorGraphCollectionIncGamma[ncent][i]->Clone(Form("GammaStat%s_%s",nameMeasGlobalLabel[i].Data(),multbinsPHOS[ncent].Data()));
            graphIndGammaIncStat_yShifted[ncent][i]             = (TGraphAsymmErrors*)statErrorGraphCollectionIncGamma[ncent][i]->Clone(Form("GammaYShiftedStat%s_%s",nameMeasGlobalLabel[i].Data(),multbinsPHOS[ncent].Data()));
        }
        graphIndGammaIncSysUnshi[ncent][i]          = NULL;
        graphIndGammaIncSys[ncent][i]               = NULL;
        graphIndGammaIncSys_yShifted[ncent][i]      = NULL;
        if (graphIncGammaSysErr[ncent][i]){
            graphIndGammaIncSysUnshi[ncent][i]                  = (TGraphAsymmErrors*)graphIncGammaSysErr[ncent][i]->Clone(Form("GammaUnshiftedSys%s_%s",nameMeasGlobalLabel[i].Data(),multbinsPHOS[ncent].Data()));
            graphIndGammaIncSys[ncent][i]                       = (TGraphAsymmErrors*)graphIncGammaSysErr[ncent][i]->Clone(Form("GammaSys%s_%s",nameMeasGlobalLabel[i].Data(),multbinsPHOS[ncent].Data()));
            graphIndGammaIncSys_yShifted[ncent][i]              = (TGraphAsymmErrors*)graphIncGammaSysErr[ncent][i]->Clone(Form("GammaYShiftedSys%s_%s",nameMeasGlobalLabel[i].Data(),multbinsPHOS[ncent].Data()));
        }
    }

    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************
    TF1* fitShiftingGamma            = FitObject("tmpt","ShiftingGamma","Gamma");
    fitShiftingGamma->SetParameters(fitTsallisGammaComb[ncent]->GetParameter(0),fitTsallisGammaComb[ncent]->GetParameter(1), fitTsallisGammaComb[ncent]->GetParameter(2));

    TGraphAsymmErrors* graphCombIncGammaTotNoShift = (TGraphAsymmErrors*) graphCombIncGammaTot[ncent]->Clone(Form("Gamma_NoShift_%s",multbins[ncent].Data()));

    // graphCombIncGammaTot[ncent]            = ApplyXshift(graphCombIncGammaTot[ncent], fitShiftingGamma);
    // cout << "comb" << endl;
    // graphCombIncGammaStat[ncent]->Print();
    // graphCombIncGammaStat[ncent]           = ApplyXshiftIndividualSpectra( graphCombIncGammaTot[ncent],
    //                                                                 graphCombIncGammaStat[ncent],
    //                                                                 fitShiftingGamma,
    //                                                                 0, graphCombIncGammaStat[ncent]->GetN());
    // graphCombIncGammaSys[ncent]            = ApplyXshiftIndividualSpectra( graphCombIncGammaTot[ncent],
    //                                                                 graphCombIncGammaSys[ncent],
    //                                                                 fitShiftingGamma,
    //                                                                 0, graphCombIncGammaSys[ncent]->GetN());
    Int_t offSetGammaShifting[11]   = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0 };
    Int_t nComBinsGammaShifting[11] = { 0, 0, 0, 0,  0,
                                        0,  0,  0,  0,  0,
                                        0 };

    // for (Int_t i = 0; i< 11; i++){
    //     if (graphIndGammaIncStat[ncent][i]){
    //         cout << "shiting stat err of " << nameMeasGlobalLabel[i].Data();
    //         graphIndGammaIncStat[ncent][i]  = ApplyXshiftIndividualSpectra(    graphCombIncGammaTot[ncent],
    //                                                                     graphIndGammaIncStat[ncent][i],
    //                                                                     fitShiftingGamma,
    //                                                                     offSetGammaShifting[i], nComBinsGammaShifting[i]);

    //     }
    //     if (graphIndGammaIncSys[ncent][i]){
    //         cout << "shiting sys err of " << nameMeasGlobalLabel[i].Data();
    //         graphIndGammaIncSys[ncent][i]   = ApplyXshiftIndividualSpectra(    graphCombIncGammaTot[ncent],
    //                                                                     graphIndGammaIncSys[ncent][i],
    //                                                                     fitShiftingGamma,
    //                                                                     offSetGammaShifting[i], nComBinsGammaShifting[i]);
    //     }
    // }


    graphRatioGammaCombCombFitTot[ncent]     = (TGraphAsymmErrors*)graphCombIncGammaTot[ncent]->Clone(Form("graphRatioGammaCombCombFitTot_%s",multbinsPHOS[ncent].Data()));
    graphRatioGammaCombCombFitTot[ncent]                        = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitTot[ncent], fitHagGammaComb[ncent]);
    graphRatioGammaCombCombFitStat[ncent]    = (TGraphAsymmErrors*)graphCombIncGammaStat[ncent]->Clone(Form("graphRatioGammaCombCombFitStat_%s",multbinsPHOS[ncent].Data()));
    graphRatioGammaCombCombFitStat[ncent]                       = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitStat[ncent], fitHagGammaComb[ncent]);
    graphRatioGammaCombCombFitSys[ncent]     = (TGraphAsymmErrors*)graphCombIncGammaSys[ncent]->Clone(Form("graphRatioGammaCombCombFitSys_%s",multbinsPHOS[ncent].Data()));
    graphRatioGammaCombCombFitSys[ncent]                        = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitSys[ncent], fitHagGammaComb[ncent]);

    for (Int_t i= 0; i< 11; i++){
        graphRatioGammaIndCombFitStat[ncent][i]            = NULL;
        if (graphIndGammaIncStat[ncent][i]){
            graphRatioGammaIndCombFitStat[ncent][i]             = (TGraphAsymmErrors*)graphIndGammaIncStat[ncent][i]->Clone(Form("RatioGamma%sToCombFitStat_%s", nameMeasGlobalLabel[i].Data(),multbinsPHOS[ncent].Data()));
            graphRatioGammaIndCombFitStat[ncent][i]             = CalculateGraphErrRatioToFit(graphRatioGammaIndCombFitStat[ncent][i], fitHagGammaComb[ncent]);
        }
        graphRatioGammaIndCombFitSys[ncent][i]             = NULL;
        if (graphIndGammaIncSys[ncent][i]){
            graphRatioGammaIndCombFitSys[ncent][i]              = (TGraphAsymmErrors*)graphIndGammaIncSys[ncent][i]->Clone(Form("RatioGamma%sToCombFitSyst_%s", nameMeasGlobalLabel[i].Data(),multbinsPHOS[ncent].Data()));
            graphRatioGammaIndCombFitSys[ncent][i]              = CalculateGraphErrRatioToFit(graphRatioGammaIndCombFitSys[ncent][i], fitHagGammaComb[ncent]);
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

    ProduceGraphAsymmWithoutXErrors(graphRatioGammaCombCombFitStat[ncent]);

    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitSys[ncent],  markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent], widthLinesBoxes, kTRUE);
    graphRatioGammaCombCombFitSys[ncent]->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitStat[ncent],  markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent]);
    graphRatioGammaCombCombFitStat[ncent]->Draw("p,same,z");

    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1.,0.1, kGray+2);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1.05, 1.05,0.1, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 0.95, 0.95,0.1, kGray, 7);

    TLatex *labelRatioToFitEnergy2      = new TLatex(0.95, 0.22, Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
    SetStyleTLatex( labelRatioToFitEnergy2, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergy2->Draw();
    TLatex *labelRatioToFitALICE2       = new TLatex(0.95, 0.16, textALICE.Data());
    SetStyleTLatex( labelRatioToFitALICE2, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitALICE2->Draw();
    TLatex *labelRatioToFitGamma        = new TLatex(0.15, 0.92, "#gamma_{inc}");
    SetStyleTLatex( labelRatioToFitGamma, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitGamma->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_PPb5023GeV_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_PPb5023GeV_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));
    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DGammaRatioToCombFit->Draw("copy");

    for (Int_t i = 10; i > -1 ; i--){
        if (graphRatioGammaIndCombFitSys[ncent][i]){ //
            DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitSys[ncent][i], markerStyleDet[i] ,markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
            graphRatioGammaIndCombFitSys[ncent][i]->Draw("E2same");
        }
        if (graphRatioGammaIndCombFitStat[ncent][i]){ //
            ProduceGraphAsymmWithoutXErrors(graphRatioGammaIndCombFitStat[ncent][i]);
            DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitStat[ncent][i], markerStyleDet[i] ,markerSizeDet[i], colorDet[i], colorDet[i]);
            graphRatioGammaIndCombFitStat[ncent][i]->Draw("p,same,z");
        }
    }
    // graphRatioGammaIndCombFitStat[ncent][4]->Draw("p,same,z");

    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1] , 1., 1.,0.5, kGray+2);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1] , 1.05, 1.05,0.5, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1] , 0.95, 0.95,0.5, kGray, 7);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1] , 1.1, 1.1,0.5, kGray, 9);
    DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1] , 0.9, 0.9,0.5, kGray, 9);

    labelRatioToFitEnergy2->Draw();
    labelRatioToFitALICE2->Draw();
    TLatex *labelRatioToFitGamma2        = new TLatex(0.95, 0.92, "#gamma_{inc}");
    SetStyleTLatex( labelRatioToFitGamma2, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitGamma2->Draw();
    histo2DGammaRatioToCombFit->Draw("same,axis");


    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    Double_t rowsLegendGammaRatio[4]          = {0.92, 0.86, 0.80, 0.74};
    Double_t rowsLegendGammaRatioAbs[4]       = {1.37, 1.26, 1.22, 1.18 };
    Double_t columnsLegendGammaRatio[6]       = {0.115, 0.232, 0.32, 0.48, 0.7, 0.8};
    Double_t columnsLegendGammaRatioAbs[6]    = {0.215, 0.89, 1.45, 2, 11.5, 17.5};
    Double_t lengthBox                        = 0.27;
    Double_t heightBox                        = 0.03/2;
    //****************** first Column **************************************************
    TLatex *textPCMRatioGamma                 = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[1],nameMeasGlobalLabel[0]);
    SetStyleTLatex( textPCMRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textPCMRatioGamma->Draw();
    TLatex *textEMCALRatioGamma               = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[2],nameMeasGlobalLabel[2]);
    SetStyleTLatex( textEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textEMCALRatioGamma->Draw();
    TLatex *textPCMEMCALRatioGamma            = new TLatex(columnsLegendGammaRatio[3],rowsLegendGammaRatio[1],nameMeasGlobalLabel[4]);
    SetStyleTLatex( textPCMEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textPCMEMCALRatioGamma->Draw();
    TLatex *textPHOSRatioGamma            = new TLatex(columnsLegendGammaRatio[3],rowsLegendGammaRatio[2],nameMeasGlobalLabel[1]);
    // TLatex *textPHOSRatioGamma            = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[3],nameMeasGlobalLabel[1]);
    SetStyleTLatex( textPHOSRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textPHOSRatioGamma->Draw();

    //****************** second Column *************************************************
    TLatex *textStatRatioGamma                = new TLatex(columnsLegendGammaRatio[1],rowsLegendGammaRatio[0] ,"stat");
    SetStyleTLatex( textStatRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textStatRatioGamma->Draw();
    TLatex *textSysRatioGamma                 = new TLatex(columnsLegendGammaRatio[2] ,rowsLegendGammaRatio[0],"syst");
    SetStyleTLatex( textSysRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textSysRatioGamma->Draw();
    TLatex *textStatRatioGamma2               = new TLatex(columnsLegendGammaRatio[4],rowsLegendGammaRatio[0] ,"stat");
    SetStyleTLatex( textStatRatioGamma2, textSizeLabelsPixel,4, 1, 43);
    textStatRatioGamma2->Draw();
    TLatex *textSysRatioGamma2                = new TLatex(columnsLegendGammaRatio[5] ,rowsLegendGammaRatio[0],"syst");
    SetStyleTLatex( textSysRatioGamma2, textSizeLabelsPixel,4, 1, 43);
    textSysRatioGamma2->Draw();

    TMarker* markerPCMGammaRatio           = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[ncent][0],columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[1],1);
    markerPCMGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[1]);
    TMarker* markerEMCALGammaRatio         = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[ncent][2], columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[2],1);
    markerEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[2]);
    TMarker* markerPCMEMCALGammaRatio      = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[ncent][4], columnsLegendGammaRatio[4] ,rowsLegendGammaRatio[1],1);
    markerPCMEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[4] ,rowsLegendGammaRatioAbs[1]);
    TMarker* markerPHOSammaRatio      = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[ncent][1], columnsLegendGammaRatio[4] ,rowsLegendGammaRatio[2],1);
    markerPHOSammaRatio->DrawMarker(columnsLegendGammaRatioAbs[4] ,rowsLegendGammaRatioAbs[2]);
    // TMarker* markerPHOSammaRatio      = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[ncent][1], columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[3],1);
    // markerPHOSammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[3]);

    TBox* boxPCMGammaRatio                 = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[ncent][0], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[1]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[1]+ heightBox);
    boxPCMGammaRatio->Draw("l");
    // TBox* boxEMCALGammaRatio               = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[ncent][2], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[2]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[2]+ heightBox);
    // boxEMCALGammaRatio->Draw("l");
    TBox* boxPCMEMCALGammaRatio            = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[ncent][4], columnsLegendGammaRatioAbs[5]-0.5*lengthBox , rowsLegendGammaRatioAbs[1]- heightBox, columnsLegendGammaRatioAbs[5]+ 25*lengthBox, rowsLegendGammaRatioAbs[1]+ heightBox);
    boxPCMEMCALGammaRatio->Draw("l");
    TBox* boxPHOSGammaRatio            = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[ncent][1], columnsLegendGammaRatioAbs[5]-0.5*lengthBox , rowsLegendGammaRatioAbs[2]- heightBox, columnsLegendGammaRatioAbs[5]+ 25*lengthBox, rowsLegendGammaRatioAbs[2]+ heightBox);
    // TBox* boxPHOSGammaRatio            = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[ncent][1], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[3]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[3]+ heightBox);
    boxPHOSGammaRatio->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));

    //*******************************************************************************************************
    //************************** Calculating combined direct photon spectrum ********************************
    //*******************************************************************************************************
    Double_t xArrayCombDR[graphCombDRStat[ncent]->GetN()+1];
    xArrayCombDR[0] = graphCombDRStat[ncent]->GetX()[0] - graphCombDRStat[ncent]->GetEXhigh()[0];
    for (Int_t i = 1; i<graphCombDRStat[ncent]->GetN()+1;i++){
        xArrayCombDR[i] = graphCombDRStat[ncent]->GetX()[i-1] + graphCombDRStat[ncent]->GetEXhigh()[i-1];
    }
    //_______________________ copy inclusive photon spectra _____________________
    TH1D *histoCombDirGammaSpectrumErrSum                   = new TH1D("histoCombDirGammaSpectrumErrSum","",graphCombDRStat[ncent]->GetN(),xArrayCombDR);
    TH1D *histoCombDirGammaSpectrumErrSys                   = new TH1D("histoCombDirGammaSpectrumErrSys","",graphCombDRStat[ncent]->GetN(),xArrayCombDR);
    TH1D *histoCombDirGammaSpectrumErrStat                  = new TH1D("histoCombDirGammaSpectrumErrStat","",graphCombDRStat[ncent]->GetN(),xArrayCombDR);

    //_______________________ get arrays of double ratio errors __________________
    Double_t *SystErrorsCombDR                              = new Double_t[graphCombIncGammaStat[ncent]->GetN()];
    Double_t *sumErrorsCombDR                               = new Double_t[graphCombIncGammaStat[ncent]->GetN()];
    Double_t *StatErrorsCombDR                              = new Double_t[graphCombIncGammaStat[ncent]->GetN()];
    Double_t *xErrorsDR                                     = new Double_t[graphCombIncGammaStat[ncent]->GetN()];
    for (Int_t i = 0; i< graphCombDRStat[ncent]->GetN(); i++){
        SystErrorsCombDR[i]                                 = graphCombDRSys[ncent]->GetEYhigh()[i]/graphCombDRSys[ncent]->GetY()[i] *100;
        StatErrorsCombDR[i]                                 = graphCombDRStat[ncent]->GetEYhigh()[i]/graphCombDRStat[ncent]->GetY()[i] *100;
        sumErrorsCombDR[i]                                  = graphCombDRTot[ncent]->GetEYhigh()[i]/graphCombDRTot[ncent]->GetY()[i] *100;
        //cout << i << "\t" << graphCombDRSys[ncent]->GetY()[i] << "\t" << graphCombDRSys[ncent]->GetEYhigh()[i] << "\t" <<SystErrorsCombDR[i] << endl;
    }
    xErrorsDR                                               = graphCombDRStat[ncent]->GetX();

    cout << __LINE__ << endl;
    //graphCombDRTot[ncent]->Print();

    //_______________________ copy inclusive photon spectra _____________________
    TH1D* histoCombErrorsForDRSum                           = new TH1D("histoCombErrorsForDRSum","",graphCombDRStat[ncent]->GetN(),xArrayCombDR);
    TH1D* histoCombErrorsForDRStat                          = new TH1D("histoCombErrorsForDRStat","",graphCombDRStat[ncent]->GetN(),xArrayCombDR);
    TH1D* histoCombErrorsForDRSys                           = new TH1D("histoCombErrorsForDRSys","",graphCombDRStat[ncent]->GetN(),xArrayCombDR);

    for(Int_t i = 1; i<graphCombDRStat[ncent]->GetN()+1;i++){
        //cout<< i << "\t"<<xErrorsDR[i-1]<<"  "<<histoCombErrorsForDRSum->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum->GetBinWidth(i) <<endl;
        Double_t binErrorSummed                                 = sumErrorsCombDR[i-1];
        Double_t binErrorSyst                                   = SystErrorsCombDR[i-1];
        Double_t binErrorStat                                   = StatErrorsCombDR[i-1];
        Double_t DR                                             = graphCombDRStat[ncent]->GetY()[i-1];

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
    //graphCombDRStat[ncent]->Print();
    for(Int_t i = 1; i<graphCombDRStat[ncent]->GetN()+1; i++){
        // obtain common quantities
        Double_t Rgamma                 = histoCombErrorsForDRSys->GetBinContent(i);
        Double_t nIncGamma              = graphCombIncGammaStat[ncent]->GetY()[i-1];

        // calculating Systematics graph
        Double_t errRgamma              = histoCombErrorsForDRSys->GetBinError(i);
        Double_t errNIncGam             = graphCombIncGammaSys[ncent]->GetEYhigh()[i-1];
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
        errNIncGam                      = graphCombIncGammaStat[ncent]->GetEYhigh()[i-1];
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
        errNIncGam                      = graphCombIncGammaTot[ncent]->GetEYhigh()[i-1];
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
    graphCombDirGammaSpectrumSystErr[ncent] = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSys,histoCombDirGammaSpecStatErr,0,0.5);
    if(graphCombDirGammaSpectrumSystErr[ncent])graphCombDirGammaSpectrumSystErr[ncent]->SetName(Form("graphCombDirGammaSpectrumSystErr_%s",multbinsPHOS[ncent].Data()));
    if(graphCombDirGammaSpectrumSystErr[ncent]) cout << "sys has been found" << endl;
    if(graphCombDirGammaSpectrumSystErr[ncent])graphCombDirGammaSpectrumSystErr[ncent]->Print();

    // purely calculating points based on Statistical errors
    graphCombDirGammaSpectrumStatErr[ncent] = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrStat,histoCombDirGammaSpecStatErr,0,0.5);
    if(graphCombDirGammaSpectrumStatErr[ncent])graphCombDirGammaSpectrumStatErr[ncent]->SetName(Form("graphCombDirGammaSpectrumStatErr_%s",multbinsPHOS[ncent].Data()));
    if(graphCombDirGammaSpectrumStatErr[ncent]) cout << "stat has been found" << endl;
    if(graphCombDirGammaSpectrumStatErr[ncent])graphCombDirGammaSpectrumStatErr[ncent]->Print();
    // purely calculating points based on all Systematic + Statistical errors
    graphCombDirGammaSpectrumSumErr[ncent] = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum,histoCombDirGammaSpecStatErr,0,0.5);
    if(graphCombDirGammaSpectrumSumErr[ncent])graphCombDirGammaSpectrumSumErr[ncent]->SetName(Form("graphCombDirGammaSpectrumSumErr_%s",multbinsPHOS[ncent].Data()));
    if(graphCombDirGammaSpectrumSumErr[ncent]) cout << "tot has been found" << endl;
    if(graphCombDirGammaSpectrumSumErr[ncent])graphCombDirGammaSpectrumSumErr[ncent]->Print();
    // calculate points above confidence level summed errors
    graphCombDirGammaSpectrumSumErrConfi[ncent] = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum,histoCombDirGammaSpecStatErr,2,0.5);
    if(graphCombDirGammaSpectrumSumErrConfi[ncent])graphCombDirGammaSpectrumSumErrConfi[ncent]->SetName(Form("graphCombDirGammaSpectrumSumErrConfi_%s",multbinsPHOS[ncent].Data()));
    if(graphCombDirGammaSpectrumSumErrConfi[ncent]) cout << "confi has been found" << endl;
    if(graphCombDirGammaSpectrumSumErrConfi[ncent])graphCombDirGammaSpectrumSumErrConfi[ncent]->Print();
    // calculate arrows for points with 0, error summed
    graphCombDirGammaSpectrumSumErrAr[ncent] = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum,histoCombDirGammaSpecStatErr,5,0.5);
    if(graphCombDirGammaSpectrumSumErrAr[ncent])graphCombDirGammaSpectrumSumErrAr[ncent]->SetName(Form("graphCombDirGammaSpectrumSumErrAr_%s",multbinsPHOS[ncent].Data()));
    if(graphCombDirGammaSpectrumSumErrAr[ncent]) cout << "Ar has been found" << endl;
    if(graphCombDirGammaSpectrumSumErrAr[ncent])graphCombDirGammaSpectrumSumErrAr[ncent]->Print();

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

        TLegend* legendDRSingle = GetAndSetLegend2(0.12,0.958-textSizeSinglePad*5,0.35,0.958, textSizeSinglePad, 1, Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()), 42, 0.3);
        legendDRSingle->SetTextAlign(11);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 0; i < 10; i++){
            if (graphDRNonFitSysErr[ncent][i]){
                DrawGammaSetMarkerTGraphAsym(graphDRNonFitSysErr[ncent][i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphDRNonFitSysErr[ncent][i]->Draw("E2same");
                legendDRSingle->AddEntry(graphDRNonFitSysErr[ncent][i],nameMeasGlobalLabel[i],"pf");

            }
            if (histoDRNonFitStatErr[ncent][i]){
                DrawGammaSetMarker(histoDRNonFitStatErr[ncent][i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoDRNonFitStatErr[ncent][i]->Draw("p,same,e0,X0");
                if (!graphDRNonFitSysErr[ncent][i])legendDRSingle->AddEntry(histoDRNonFitStatErr[ncent][i],nameMeasGlobalLabel[i],"p");
            }
        }

        if (histoDRNonFitStatErr[ncent][2]) histoDRNonFitStatErr[ncent][2]->Draw("p,same,e0,X0");
        if (histoDRNonFitStatErr[ncent][0]) histoDRNonFitStatErr[ncent][0]->Draw("p,same,e0,X0");
        if (histoDRNonFitStatErr[ncent][4]) histoDRNonFitStatErr[ncent][4]->Draw("p,same,e0,X0");
        legendDRSingle->Draw();

        TLatex *labelDRSingle       = new TLatex(0.95,0.19,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelDRSingle, textSizeSinglePad,4, 1, 42, kTRUE, 31);
//         labelDRSingle->Draw();
        TLatex *labelALICEDRSingle  = new TLatex(0.95,0.14,textALICE.Data());
        SetStyleTLatex( labelALICEDRSingle, textSizeSinglePad,4, 1, 42, kTRUE, 31);
        labelALICEDRSingle->Draw();

        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pPb5TeV_%s.%s", outputDir.Data(), multbinsPHOS[ncent].Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pPb5TeV_%s.pdf", outputDir.Data(), multbinsPHOS[ncent].Data()));



    hist2DDRDummySingle->DrawCopy();

        graphCombDRStatPlot[ncent]    = (TGraphAsymmErrors*)graphCombDRStat[ncent]->Clone(Form("graphCombDRStatPlot_%s",multbinsPHOS[ncent].Data()));
        ProduceGraphAsymmWithoutXErrors(graphCombDRStatPlot[ncent]);

        TLegend* legendDRComb = GetAndSetLegend2(0.12,0.95-textSizeSinglePad*1,0.5,0.95, textSizeSinglePad, 1, "", 42, 0.15);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        DrawGammaSetMarkerTGraphAsym(graphCombDRSys[ncent], markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent],widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRStatPlot[ncent], markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent], widthLinesBoxes);
        legendDRComb->AddEntry(graphCombDRSys[ncent],textALICE.Data(),"pf");
        graphCombDRSys[ncent]->Draw("E2same");
        graphCombDRStatPlot[ncent]->Draw("z,p,same");
//         legendDRComb->Draw();

        labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_Comb_pPb5TeV_%s.%s", outputDir.Data(), multbinsPHOS[ncent].Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_Comb_pPb5TeV_%s.pdf", outputDir.Data(), multbinsPHOS[ncent].Data()));


        hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRTheoryComb     = GetAndSetLegend2(0.12,0.958-textSizeSinglePad*3,0.5,0.958, textSizeSinglePad, 1, Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()), 42, 0.15);
        TLegend* legendDRTheoryComb2    = GetAndSetLegend2(0.12,0.958-textSizeSinglePad*7,0.5,0.958-textSizeSinglePad*3, textSizeSinglePad, 1, "NLO pQCD: ", 42, 0.15);
        TGraphAsymmErrors* dummyNLO         = new TGraphAsymmErrors(1);
        TGraphAsymmErrors* dummyMcGill      = new TGraphAsymmErrors(1);
        TGraphAsymmErrors* dummyNLOnCTEQ    = new TGraphAsymmErrors(1);
        TGraphAsymmErrors* dummyNLOEPPS     = new TGraphAsymmErrors(1);

        legendDRTheoryComb->AddEntry(graphCombDRSys[ncent],textALICE.Data(),"pf");

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

        graphCombDRSys[ncent]->Draw("E2same");
        graphCombDRStatPlot[ncent]->Draw("p,z,same");
        legendDRTheoryComb->Draw();
        legendDRTheoryComb2->Draw();

        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pPb5TeV_%s.%s", outputDir.Data(), multbinsPHOS[ncent].Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pPb5TeV_%s.pdf", outputDir.Data(), multbinsPHOS[ncent].Data()));

    // **********************************************************************************************************************
    // ******************************** Efficiency for gamma individual measurements ****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 55;
    Double_t textSizeLabelsRel          = 55./1200;
    cout << textSizeLabelsRel << endl;

    TCanvas* canvasEff   = new TCanvas("canvasEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEff,  0.09, 0.01, 0.015, 0.095);
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
            if (histoEffi[ncent][i]){
                DrawGammaSetMarker(histoEffi[ncent][i],   markerStyleDetMC[i], markerSizeDetMC[i], colorDetMC[i] , colorDetMC[i]);
                histoEffi[ncent][i]->Draw("p,same,e");
                legendEffiGamma->AddEntry(histoEffi[ncent][i],"   ","p");
            } else if (histoEffiMCPt[ncent][i]){
                legendEffiGamma->AddEntry((TObject*)0,"   ","");
            }
            if (histoEffiMCPt[ncent][i]){
                DrawGammaSetMarker(histoEffiMCPt[ncent][i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoEffiMCPt[ncent][i]->Draw("p,same,e");
                legendEffiGamma->AddEntry(histoEffiMCPt[ncent][i],"    "+nameMeasGlobalLabel[i],"p");
            } else if (histoEffi[ncent][i]){
                legendEffiGamma->AddEntry((TObject*)0,"    "+nameMeasGlobalLabel[i],"");
            }
        }
        legendEffiGamma->Draw();

        TLatex *labelPerfEffi           = new TLatex(0.13,0.92,textALICE.Data());
        SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
        labelPerfEffi->Draw();
        TLatex *labelEnergyEffi         = new TLatex(0.13,0.87,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
        labelEnergyEffi->Draw();
        TLatex *labelPerfEffiPTrec      = new TLatex(0.57,0.145+(4*textSizeLabelsRel),"#it{p}_{T}^{rec}");
        SetStyleTLatex( labelPerfEffiPTrec, textSizeLabelsRel,4);
        labelPerfEffiPTrec->Draw();
        TLatex *labelPerfEffiPTtrue      = new TLatex(0.665,0.145+(4*textSizeLabelsRel),"#it{p}_{T}^{true}");
        SetStyleTLatex( labelPerfEffiPTtrue, textSizeLabelsRel,4);
        labelPerfEffiPTtrue->Draw();

    canvasEff->Update();
    canvasEff->Print(Form("%s/Gamma_Effiency_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasEff->Print(Form("%s/Gamma_Effiency_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));

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
            if (histoResolCorr[ncent][i]){
                DrawGammaSetMarker(histoResolCorr[ncent][i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoResolCorr[ncent][i]->GetXaxis()->SetRangeUser(graphIndGammaIncStat[ncent][i]->GetX()[0]-graphIndGammaIncStat[ncent][i]->GetEXlow()[0],
                                                            graphIndGammaIncStat[ncent][i]->GetX()[graphIndGammaIncStat[ncent][i]->GetN()-1]+graphIndGammaIncStat[ncent][i]->GetEXhigh()[graphIndGammaIncStat[ncent][i]->GetN()-1] );

                histoResolCorr[ncent][i]->Draw("p,same,e");
                legendResolCorGamma->AddEntry(histoResolCorr[ncent][i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendResolCorGamma->Draw();
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        TLatex *labelPerfResolCor           = new TLatex(0.15,0.92,textALICE.Data());
        SetStyleTLatex( labelPerfResolCor, textSizeLabelsRel,4);
        labelPerfResolCor->Draw();
        TLatex *labelEnergyResolCor         = new TLatex(0.15,0.87,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelEnergyResolCor, textSizeLabelsRel,4);
        labelEnergyResolCor->Draw();

    histo2DResCor->Draw("same,axis");
    canvasResolCor->Update();
    canvasResolCor->Print(Form("%s/Gamma_ResolutionCorrection_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasResolCor->Print(Form("%s/Gamma_ResolutionCorrection_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));

    // **********************************************************************************************************************
    // ******************************** Purity for gamma individual measurements ****************************************
    // **********************************************************************************************************************
    TCanvas* canvasPurity   = new TCanvas("canvasPurity", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasPurity,  0.1, 0.01, 0.015, 0.095);
    canvasPurity->SetLogx(1);

    TH1F * histo1DPurity            = new TH1F("histo1DPurity", "histo1DPurity",1000, doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(  histo1DPurity, "#it{p}_{T} (GeV/#it{c})","#it{#varepsilon}_{pur}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
    histo1DPurity->GetYaxis()->SetRangeUser(0.8, 1.07 );
    histo1DPurity->GetYaxis()->SetLabelOffset(0.001);
    histo1DPurity->GetXaxis()->SetLabelOffset(-0.01);
    histo1DPurity->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo1DPurity->DrawCopy();

        TLegend* legendPurityGamma           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 0; i < 11; i++){
            if (histoPurity[ncent][i]){
                DrawGammaSetMarker(histoPurity[ncent][i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoPurity[ncent][i]->Draw("p,same,e");
                legendPurityGamma->AddEntry(histoPurity[ncent][i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendPurityGamma->Draw();

        TLatex *labelPerfPurity           = new TLatex(0.15,0.92,textALICE.Data());
        SetStyleTLatex( labelPerfPurity, textSizeLabelsRel,4);
        labelPerfPurity->Draw();
        TLatex *labelEnergyPurity         = new TLatex(0.15,0.87,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelEnergyPurity, textSizeLabelsRel,4);
        labelEnergyPurity->Draw();

    canvasPurity->Update();
    canvasPurity->Print(Form("%s/Gamma_Purity_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasPurity->Print(Form("%s/Gamma_Purity_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));

    // **********************************************************************************************************************
    // ******************************** ConvProb for gamma individual measurements ******************************************
    // **********************************************************************************************************************
    TF1* fitConvProbPCM                     = new TF1("fitConvProbPCM","[0]");
    if (histoConvProb[ncent][0]) histoConvProb[ncent][0]->Fit(fitConvProbPCM,"NRMEX0+","",4,10.);
    else if  (histoConvProb[ncent][4]) histoConvProb[ncent][4]->Fit(fitConvProbPCM,"NRMEX0+","",4,10.);
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
            if (histoConvProb[ncent][i]){
                DrawGammaSetMarker(histoConvProb[ncent][i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoConvProb[ncent][i]->Draw("p,same,e");
                legendConvProbGamma->AddEntry(histoConvProb[ncent][i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendConvProbGamma->Draw();

        TLatex *labelPerfConvProb           = new TLatex(0.15,0.92,textALICE.Data());
        SetStyleTLatex( labelPerfConvProb, textSizeLabelsRel,4);
        labelPerfConvProb->Draw();
        TLatex *labelEnergyConvProb         = new TLatex(0.15,0.87,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelEnergyConvProb, textSizeLabelsRel,4);
        labelEnergyConvProb->Draw();

    histo1DConvProb->Draw("same,axis");

    canvasConvProb->Update();
    canvasConvProb->Print(Form("%s/Gamma_ConvProb_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasConvProb->Print(Form("%s/Gamma_ConvProb_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));
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
            if (histoPileupCorr[ncent][i]){
                DrawGammaSetMarker(histoPileupCorr[ncent][i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoPileupCorr[ncent][i]->Draw("p,same,e");
                legendPileUpGamma->AddEntry(histoPileupCorr[ncent][i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendPileUpGamma->Draw();

        TLatex *labelPerfPileUp           = new TLatex(0.15,0.92,textALICE.Data());
        SetStyleTLatex( labelPerfPileUp, textSizeLabelsRel,4);
        labelPerfPileUp->Draw();
        TLatex *labelEnergyPileUp         = new TLatex(0.15,0.87,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelEnergyPileUp, textSizeLabelsRel,4);
        labelEnergyPileUp->Draw();

    canvasPileUp->Update();
    canvasPileUp->Print(Form("%s/Gamma_PileUp_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasPileUp->Print(Form("%s/Gamma_PileUp_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));


    // **********************************************************************************************************************
    // ******************************** Total Corr for gamma individual measurements ****************************************
    // **********************************************************************************************************************

    TCanvas* canvasTotalCorr   = new TCanvas("canvasTotalCorr", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasTotalCorr,  0.115, 0.01, 0.015, 0.095);
    canvasTotalCorr->SetLogy(1);
    canvasTotalCorr->SetLogx(1);

    TH1F * histo1DTotalCorr            = new TH1F("histo1DTotalCorr", "histo1DTotalCorr",1000, doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(  histo1DTotalCorr, "#it{p}_{T} (GeV/#it{c})","#it{A}#upoint#it{#varepsilon}_{rec}#upoint#it{P}_{conv}/#it{#varepsilon}_{pur}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.25, 510,505);
    histo1DTotalCorr->GetYaxis()->SetRangeUser(0.005, 0.49 );
    histo1DTotalCorr->GetYaxis()->SetLabelOffset(0.005);
    histo1DTotalCorr->GetXaxis()->SetNoExponent();
    histo1DTotalCorr->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo1DTotalCorr->DrawCopy();

    TLegend* legendTotalCorrGamma           = GetAndSetLegend2(0.62, 0.87, 0.96, 0.952, textSizeLabelsPixel, 2, "", 43, 0.32);
    legendTotalCorrGamma->SetTextAlign(13);
    for (Int_t i = 0; i < 11; i++){
        if( i == 4) continue;
        if (i == 2) legendTotalCorrGamma->AddEntry((TObject*)0,"","");
        if (histoTotalCorrFactor[ncent][i]){
            DrawGammaSetMarker(histoTotalCorrFactor[ncent][i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
            histoTotalCorrFactor[ncent][i]->Draw("p,same,e");
            legendTotalCorrGamma->AddEntry(histoTotalCorrFactor[ncent][i],nameMeasGlobalLabel[i],"p");
        }
    }
    legendTotalCorrGamma->Draw();

    TLatex *labelPerfTotalCorr           = new TLatex(0.15,0.92,"ALICE simulation");
    SetStyleTLatex( labelPerfTotalCorr, textSizeLabelsRel,4);
    labelPerfTotalCorr->Draw();
    TLatex *labelEnergyTotalCorr         = new TLatex(0.15,0.87,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
    SetStyleTLatex( labelEnergyTotalCorr, textSizeLabelsRel,4);
    labelEnergyTotalCorr->Draw();

    canvasTotalCorr->Update();
    canvasTotalCorr->Print(Form("%s/Gamma_TotalCorrFactor_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasTotalCorr->Print(Form("%s/Gamma_TotalCorrFactor_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));


    // **********************************************************************************************************************
    // ******************************** EffectiveSecCorr for gamma individual measurements ****************************************
    // **********************************************************************************************************************
    // TCanvas* canvasEffectiveSecCorr   = new TCanvas("canvasEffectiveSecCorr", "", 200, 10, 1200, 1100);  // gives the page size
    // DrawGammaCanvasSettings( canvasEffectiveSecCorr,  0.11, 0.01, 0.04, 0.095);
    // canvasEffectiveSecCorr->SetLogx(1);

    // TH1F * histo1DEffectiveSecCorr            = new TH1F("histo1DEffectiveSecCorr", "histo1DEffectiveSecCorr",1000, doubleRatioXpp[0], doubleRatioXpp[1]);
    // SetStyleHistoTH1ForGraphs(  histo1DEffectiveSecCorr, "#it{p}_{T} (GeV/#it{c})","#it{C}_{sec}",
    //                             0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.08);//(#times #epsilon_{pur})
    // histo1DEffectiveSecCorr->GetYaxis()->SetLabelOffset(0.003);
    // histo1DEffectiveSecCorr->GetXaxis()->SetLabelOffset(-0.01);
    // histo1DEffectiveSecCorr->GetXaxis()->SetMoreLogLabels(kTRUE);


    // Double_t minYSecCorr[4]             = {0.0, 0.0, 0.0, 0.0};
    // Double_t maxYSecCorr[4]             = {0.05, 0.004, 0.8e-3, 0.041};
    // TString nameLabelSec[4]             = {"K^{0}_{s}", "K^{0}_{l}", "#Lambda", "Rest"};
    // TString nameOutputSec[4]            = {"K0s", "K0l", "Lambda", "Rest"};
    // for (Int_t k = 0; k < 4; k++){
    //     histo1DEffectiveSecCorr->GetYaxis()->SetTitle(Form("C_{sec,%s}", nameLabelSec[k].Data()));
    //     histo1DEffectiveSecCorr->GetYaxis()->SetRangeUser(minYSecCorr[k], maxYSecCorr[k] );
    //     histo1DEffectiveSecCorr->DrawCopy();
    //     TLegend* legendEffectiveSecCorrGamma            = GetAndSetLegend2(0.65, 0.925-(3*textSizeLabelsRel), 0.93, 0.925,textSizeLabelsPixel);
    //     for (Int_t i = 0; i < 11; i++){
    //         if (histoEffSecCorr[ncent][k][i]){
    //             DrawGammaSetMarker(histoEffSecCorr[ncent][k][i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
    //             histoEffSecCorr[ncent][k][i]->Draw("p,same,e");
    //             legendEffectiveSecCorrGamma->AddEntry(histoEffSecCorr[ncent][k][i],nameMeasGlobalLabel[i],"p");
    //         }
    //     }
    //     legendEffectiveSecCorrGamma->Draw();

    //     TLatex *labelPerfEffectiveSecCorr           = new TLatex(0.15,0.89,"ALICE performance");
    //     SetStyleTLatex( labelPerfEffectiveSecCorr, textSizeLabelsRel,4);
    //     labelPerfEffectiveSecCorr->Draw();
    //     TLatex *labelEnergyEffectiveSecCorr         = new TLatex(0.15,0.84,Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
    //     SetStyleTLatex( labelEnergyEffectiveSecCorr, textSizeLabelsRel,4);
    //     labelEnergyEffectiveSecCorr->Draw();

    //     canvasEffectiveSecCorr->Update();
    //     canvasEffectiveSecCorr->Print(Form("%s/Gamma_EffectiveSecCorr_%s.%s",outputDir.Data(), nameOutputSec[k].Data(), suffix.Data()));
    //     canvasEffectiveSecCorr->Print(Form("%s/Gamma_EffectiveSecCorr_%s.pdf",outputDir.Data(), nameOutputSec[k].Data()));
    // }

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

        graphCombIncGammaStatPlot[ncent]    = (TGraphAsymmErrors*)graphCombIncGammaStat[ncent]->Clone(Form("graphCombIncGammaStatPlot_%s",multbins[ncent].Data()));
        ProduceGraphAsymmWithoutXErrors(graphCombIncGammaStatPlot[ncent]);

        TLegend* legendYieldIncGamma           = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(4*textSizeLabelsRel*0.85),textSizeLabelsPixel);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys[ncent], markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot[ncent], markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes);
        graphCombIncGammaSys[ncent]->Draw("E2same");
        graphCombIncGammaStatPlot[ncent]->Draw("Epsame");
        legendYieldIncGamma->AddEntry(graphCombIncGammaSys[ncent],textALICE.Data(),"pf");
        legendYieldIncGamma->Draw();

        DrawGammaSetMarkerTF1( fitHagGammaComb[ncent], 7, 2, colorCombpPb);
        legendYieldIncGamma->AddEntry(fitHagGammaComb[ncent],"Hagedorn fit","l");
        fitHagGammaComb[ncent]->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitHagGammaComb[ncent]->Draw("same");
        DrawGammaSetMarkerTF1( fitTCMGammaComb[ncent], 9, 2, kGray+1);
        legendYieldIncGamma->AddEntry(fitTCMGammaComb[ncent],"TCM fit","l");
        fitTCMGammaComb[ncent]->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTCMGammaComb[ncent]->Draw("same");
        DrawGammaSetMarkerTF1( fitTsallisGammaComb[ncent], 5, 2, colorCombpPbBox);
        legendYieldIncGamma->AddEntry(fitTsallisGammaComb[ncent],"Tsalis fit","l");
        fitTsallisGammaComb[ncent]->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTsallisGammaComb[ncent]->Draw("same");

        TLatex *labelEnergyInvYieldPaperAll = new TLatex(0.20, 0.20+0.04*3, Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        TLatex *labelEnergyInvYieldPaperAll4 = new TLatex(0.20, 0.20+0.04*4, Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
        SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        SetStyleTLatex( labelEnergyInvYieldPaperAll4, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelEnergyInvYieldPaperAll->Draw();
        TLatex *labelALICEInvYieldPaperAll  = new TLatex(0.20,0.20+0.04*4,textALICE.Data());
        SetStyleTLatex( labelALICEInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
//         labelALICEInvYieldPaperAll->Draw();
        TLatex *labelALICENormUnPaperAll    = new TLatex(0.20,0.20+0.05+0.04,"Norm. unc. 3.1%");
        SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICENormUnPaperAll->Draw();

        histo2DYieldGamma->Draw("same,axis");

    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(), suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));
    // **********************************************************************************************************************
    // ******************************** InvYield for individual inc gamma measurements **************************************
    // **********************************************************************************************************************
    histo2DYieldGamma->Draw("copy");

        TLegend* legendYieldIncGammaInd       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(4*textSizeLabelsRel*0.85),textSizeLabelsPixel);
        for (Int_t i = 0; i < 11; i++){
//             if (i == 4) continue;
            if (graphIndGammaIncSys[ncent][i]){
                DrawGammaSetMarkerTGraphAsym(graphIndGammaIncSys[ncent][i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphIndGammaIncSys[ncent][i]->Draw("E2same");
                legendYieldIncGammaInd->AddEntry(graphIndGammaIncSys[ncent][i], Form("#gamma_{inc} %s", nameMeasGlobalLabel[i].Data()),"pf");
            }
            if (graphIndGammaIncStat[ncent][i]){
                DrawGammaSetMarkerTGraphAsym(graphIndGammaIncStat[ncent][i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                graphIndGammaIncStat[ncent][i]->Draw("Epsame,x0");
                if (!graphIncGammaSysErr[i])legendYieldIncGammaInd->AddEntry(graphIndGammaIncStat[ncent][i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendYieldIncGammaInd->Draw();

        labelEnergyInvYieldPaperAll->Draw();
        labelALICEInvYieldPaperAll->Draw();
        labelALICENormUnPaperAll->Draw();

        histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_IndMeas_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_IndMeas_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));

    // **********************************************************************************************************************
    // ******************************** InvYield for combined dir gamma measurement **************************************
    // **********************************************************************************************************************

    // prep plot graphs
    if (graphCombDirGammaSpectrumStatErr[ncent]) graphCombDirGammaSpectrumStatErrPlot[ncent]       = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr[ncent]->Clone(Form("graphCombDirGammaSpectrumStatErrPlot_%s",multbinsPHOS[ncent].Data()));
    if (graphCombDirGammaSpectrumStatErrPlot[ncent]) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErrPlot[ncent]);

    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-12,10e1);
    histo2DYieldGamma->Draw("copy");

//     TLegend* legendYieldIncGammaInd       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(3*textSizeLabelsRel*0.85),textSizeLabelsPixel);
//     legendYieldIncGammaInd->Draw();

    if (graphCombDirGammaSpectrumSystErr[ncent]){
        DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr[ncent], markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent], widthLinesBoxes, kTRUE);
        graphCombDirGammaSpectrumSystErr[ncent]->Draw("E2same");
    }
    if (graphCombDirGammaSpectrumStatErrPlot[ncent]){
        DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlot[ncent], markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent]);
        graphCombDirGammaSpectrumStatErrPlot[ncent]->Draw("p,E1Z,same");
    }
    if (graphCombDirGammaSpectrumSumErrAr[ncent]){
        DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrAr[ncent] , 1, 3,  colorMult[ncent] , colorMult[ncent], 1.8, kTRUE);
        graphCombDirGammaSpectrumSumErrAr[ncent]->Draw(">,same");
        PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr[ncent]);
    }


    TLatex *labelEnergyDGInvYieldPaperAll = new TLatex(0.94, 0.965-0.04*1, Form("%s %s",multbins[ncent].Data(),collisionSystempPb.Data()));
    SetStyleTLatex( labelEnergyDGInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
    labelEnergyDGInvYieldPaperAll->Draw();
    TLatex *labelALICEDGInvYieldPaperAll  = new TLatex(0.94,0.965-0.04*2,textALICE.Data());
    SetStyleTLatex( labelALICEDGInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
    labelALICEDGInvYieldPaperAll->Draw();
    TLatex *labelALICEDGNormUnPaperAll    = new TLatex(0.94,0.965-(0.04*2+0.04*0.8),"Norm. unc. 3.1%");
    SetStyleTLatex( labelALICEDGNormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 31);
    labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));
    histo2DYieldGamma->Draw("copy");

    TLegend* legendYieldDirGamma        = GetAndSetLegend2(0.70, 0.84-(3.03*textSizeLabelsRel), 0.93, 0.84,textSizeLabelsPixel, 1, "", 43, 0.3);
        legendYieldDirGamma->SetTextAlign(12);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys[ncent],  markerStyleMultOpen[ncent], markerSizeMult[ncent]+0.2, colorMult[ncent] , colorMult[ncent],widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot[ncent], markerStyleMultOpen[ncent], markerSizeMult[ncent]+0.2, colorMult[ncent] , colorMult[ncent], widthLinesBoxes);
        legendYieldDirGamma->AddEntry(graphCombIncGammaSys[ncent], "#gamma_{inc}","pf");
        graphCombIncGammaSys[ncent]->Draw("E2same");
        graphCombIncGammaStatPlot[ncent]->Draw("Epsame");

        fitTCMGammaComb[ncent]->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTCMGammaComb[ncent]->Draw("same");
        legendYieldDirGamma->AddEntry(fitTCMGammaComb[ncent],"TCM fit","l");

        if (graphCombDirGammaSpectrumSystErr[ncent]){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr[ncent], markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent], widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr[ncent]->Draw("E2same");
            legendYieldDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr[ncent], "#gamma_{dir}","pf");
        }
        if (graphCombDirGammaSpectrumStatErrPlot[ncent]){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlot[ncent], markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent]);
            graphCombDirGammaSpectrumStatErrPlot[ncent]->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr[ncent]){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrAr[ncent] , 1, 3,  colorMult[ncent] , colorMult[ncent], 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErrAr[ncent]->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr[ncent]);
        }
        legendYieldDirGamma->Draw();

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));


     histo2DYieldGamma->Draw("copy");

        TLegend* legendYieldDirGamma2           = GetAndSetLegend2(0.70, 0.84-(2.4*textSizeLabelsRel*0.85), 0.93, 0.84,textSizeLabelsPixel, 1, "", 43, 0.3);
        TLegend* legendYieldDirGamma3           = GetAndSetLegend2(0.20, 0.11+(5.2*textSizeLabelsRel*0.85), 0.5, 0.11+(6.2*textSizeLabelsRel*0.85),textSizeLabelsPixel,1,
                                                                   "", 43, 0.23);
        TLegend* legendYieldDirGammaTheo        = GetAndSetLegend2(0.20, 0.11+(4.2*textSizeLabelsRel*0.85), 0.5, 0.11+(5.2*textSizeLabelsRel*0.85),textSizeLabelsPixel,1,
                                                                  "", 43, 0.23);
        TLegend* legendYieldDirGammaTheo2       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(4*textSizeLabelsRel*0.85),textSizeLabelsPixel,1, "#gamma_{dir} NLO pQCD:", 43, 0.23);
        legendYieldDirGammaTheo2->SetTextAlign(12);
        graphCombIncGammaSys[ncent]->Draw("E2same");
        graphCombIncGammaStatPlot[ncent]->Draw("Epsame");
        legendYieldDirGamma2->AddEntry(graphCombIncGammaSys[ncent], "#gamma_{inc}","pf");

        fitTCMGammaComb[ncent]->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTCMGammaComb[ncent]->Draw("same");
        legendYieldDirGamma2->AddEntry(fitTCMGammaComb[ncent],"TCM fit","l");


        if (graphCombDirGammaSpectrumSystErr[ncent]){
            graphCombDirGammaSpectrumSystErr[ncent]->Draw("E2same");
            legendYieldDirGamma3->AddEntry(graphCombDirGammaSpectrumSystErr[ncent], "#gamma_{dir}","pf");
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

        if (graphCombDirGammaSpectrumStatErrPlot[ncent]){
            graphCombDirGammaSpectrumStatErrPlot[ncent]->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr[ncent]){
            graphCombDirGammaSpectrumSumErrAr[ncent]->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr[ncent]);
        }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));



    // **********************************************************************************************************************
    // ******************************** GammaRawYields for individual inc gamma measurements ********************************
    // **********************************************************************************************************************
    histo2DYieldGamma->GetYaxis()->SetTitle("#frac{#it{N}_{raw,#gamma}}{#it{N}_{ev}}");
    histo2DYieldGamma->Draw("copy");
    for (Int_t i = 0; i < 11; i++){
        if (histoRawGamma[ncent][i]){
            cout << nameMeasGlobalLabel[i].Data() << endl;
            TGraphErrors* dummy         = new TGraphErrors(histoRawGamma[ncent][i]);
            dummy->Print();
            DrawGammaSetMarker(histoRawGamma[ncent][i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
            histoRawGamma[ncent][i]->Draw("p,same,e0,X0");
            delete dummy;
        }
    }
    legendYieldIncGammaInd->Draw();

    labelEnergyInvYieldPaperAll->Draw();
    labelALICEInvYieldPaperAll->Draw();
    labelALICENormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/RawYield_IncGamma_IndMeas_%s.%s",outputDir.Data(),multbinsPHOS[ncent].Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/RawYield_IncGamma_IndMeas_%s.pdf",outputDir.Data(),multbinsPHOS[ncent].Data()));
}/*
    // ****************************************************************************************************************
    // ************************** Store final results including corr factors in 1 file ********************************
    // ****************************************************************************************************************
    TString nameOutputCommonFile    = Form("CombinedGammaResultPPb5TeV_%s.root", dateForOutput.Data());

    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    fCombResults.mkdir("Gamma_pPb5TeV");
    TDirectoryFile* directoryGamma = (TDirectoryFile*)fCombResults.Get("Gamma_pPb5TeV");
    fCombResults.cd("Gamma_pPb5TeV");

    // writing main results
    if (graphCombIncGammaStat[ncent]) graphCombIncGammaStat[ncent]->Write("graphInvYieldIncGammaStatErr");
    if (graphCombIncGammaSys[ncent]) graphCombIncGammaSys[ncent]->Write("graphInvYieldIncGammaSysErr");
    if (graphCombIncGammaTot[ncent]) graphCombIncGammaTot[ncent]->Write("graphInvYieldIncGammaTotErr");
    if (graphCombIncGammaStatUnshi[ncent]) graphCombIncGammaStatUnshi[ncent]->Write("graphInvYieldIncGammaStatErr_Unshifted");
    if (graphCombIncGammaSysUnshi[ncent]) graphCombIncGammaSysUnshi[ncent]->Write("graphInvYieldIncGammaSysErr_Unshifted");
    if (graphCombIncGammaTotUnshi[ncent]) graphCombIncGammaTotUnshi[ncent]->Write("graphInvYieldIncGammaTotErr_Unshifted");

    if (graphCombDRStat[ncent]) graphCombDRStat[ncent]->Write("graphRGammaCombNonFitStatErr");
    if (graphCombDRSys[ncent]) graphCombDRSys[ncent]->Write("graphRGammaCombNonFitSysErr");
    if (graphCombDRTot[ncent]) graphCombDRTot[ncent]->Write("graphRGammaCombNonFitTotErr");
    if (graphCombDirGammaSpectrumSystErr[ncent]) graphCombDirGammaSpectrumSystErr[ncent]->Write("graphInvYieldNonFitDirGammaSysErr");
    if (graphCombDirGammaSpectrumStatErr[ncent]) graphCombDirGammaSpectrumStatErr[ncent]->Write("graphInvYieldNonFitDirGammaStatErr");
    if (graphCombDirGammaSpectrumSumErrAr[ncent]) graphCombDirGammaSpectrumSumErrAr[ncent]->Write("graphInvYieldNonFitDirGammaSumErrAr");
    if (graphCombDRStat[ncent]) graphCombDRStat[ncent]->Write("graphRGammaCombStatErr");
    if (graphCombDRSys[ncent]) graphCombDRSys[ncent]->Write("graphRGammaCombSysErr");
    if (graphCombDRTot[ncent]) graphCombDRTot[ncent]->Write("graphRGammaCombTotErr");
    if (graphCombDirGammaSpectrumSystErr) graphCombDirGammaSpectrumSystErr->Write("graphInvYieldDirGammaSysErr");
    if (graphCombDirGammaSpectrumStatErr) graphCombDirGammaSpectrumStatErr->Write("graphInvYieldDirGammaStatErr");
    if (graphCombDirGammaSpectrumSumErrAr) graphCombDirGammaSpectrumSumErrAr->Write("graphInvYieldDirGammaSumErrAr");

    for (Int_t i = 0; i < 11; i++){
        if (graphIndGammaIncSys[ncent][i]) graphIndGammaIncSys[ncent][i]->Write(Form("graphInvYieldIncGamma%sSysErr",nameMeasGlobal[i].Data()));
        if (graphIndGammaIncStat[ncent][i]) graphIndGammaIncStat[ncent][i]->Write(Form("graphInvYieldIncGamma%sStatErr",nameMeasGlobal[i].Data()));
        if (histoIncGammaStatErr[i]) histoIncGammaStatErr[i]->Write(Form("histoInvYieldIncGamma%sStatErr_Unshifted",nameMeasGlobal[i].Data()));
        if (histoDRStatErr[i]) histoDRStatErr[i]->Write(Form("histoRGamma%sStatErr",nameMeasGlobal[i].Data()));
        if (graphDRSysErr[i]) graphDRSysErr[i]->Write(Form("graphRGamma%sSysErr",nameMeasGlobal[i].Data()));
        if (histoDRNonFitStatErr[i]) histoDRNonFitStatErr[i]->Write(Form("histoNonFitRGamma%sStatErr",nameMeasGlobal[i].Data()));
        if (graphDRNonFitSysErr[i]) graphDRNonFitSysErr[i]->Write(Form("graphNonfitRGamma%sSysErr",nameMeasGlobal[i].Data()));
    }
    directoryGamma->mkdir("Supporting");
    directoryGamma->cd("Supporting");
    // Writing full correction factors
    for (Int_t i = 0; i < 11; i++){
        if (histoPileupCorr[ncent][i]) histoPileupCorr[ncent][i]->Write(Form("histoPileupCorr%s",nameMeasGlobal[i].Data()));
        if (histoConvProb[ncent][i]) histoConvProb[ncent][i]->Write(Form("histoConvProb%s",nameMeasGlobal[i].Data()));
        if (histoEffi[ncent][i]) histoEffi[ncent][i]->Write(Form("histoEffi[ncent]Effective%s",nameMeasGlobal[i].Data()));
        if (histoEffiMCPt[ncent][i]) histoEffiMCPt[ncent][i]->Write(Form("histoEffiWithoutResol%s",nameMeasGlobal[i].Data()));
        if (histoPurity[ncent][i]) histoPurity[ncent][i]->Write(Form("histoPurity%s",nameMeasGlobal[i].Data()));
        for (Int_t k = 0; k<4; k++){
            if (histoEffSecCorr[ncent][k][i]) histoEffSecCorr[ncent][k][i]->Write(Form("histoEffectiveSecCorrFrom%s_%s",nameOutputSec[i].Data(), nameMeasGlobal[i].Data()));
        }
    }

    fCombResults.Close();
    */


    // **********************************************************************************************************************
    // ******************************** InvYield for combined inc gamma measurement *****************************************
    // **********************************************************************************************************************
    TGraphAsymmErrors* graphCombIncGammaSysScaled[5];
    TGraphAsymmErrors* graphCombIncGammaStatPlotScaled[5];
    TF1* fitHagGammaCombScaled[5];
        
    TCanvas* canvasSpectraGamma          = new TCanvas("canvasSpectraGamma","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasSpectraGamma, 0.16, 0.02, 0.02, 0.08);
    canvasSpectraGamma->SetLogx();
    canvasSpectraGamma->SetLogy();
    textSizeLabelsPixel                 = 55;
    Double_t textSizeLabelsRel          = 55./1200;

    TH1F * histo2DSpectraGamma              = new TH1F("histo2DSpectraGamma","histo2DSpectraGamma",11000,doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(histo2DSpectraGamma, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                              0.035,0.04, 0.035,0.04, 0.9, 1.7);
    histo2DSpectraGamma->GetYaxis()->SetRangeUser(2e-10,5e4);
    histo2DSpectraGamma->GetXaxis()->SetMoreLogLabels();
    histo2DSpectraGamma->GetXaxis()->SetLabelOffset(-0.01);
    histo2DSpectraGamma->Draw("copy");
    TLegend* legendYieldIncGamma           = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(4*textSizeLabelsRel*0.85),textSizeLabelsPixel);
    for(Int_t ncent = 0; ncent < 4; ncent++){
        graphCombIncGammaSysScaled[ncent]    = ScaleGraph(graphCombIncGammaSys[ncent],pow(10,3-ncent));
        graphCombIncGammaStatPlotScaled[ncent]    = ScaleGraph(graphCombIncGammaStatPlot[ncent],pow(10,3-ncent));
        ProduceGraphAsymmWithoutXErrors(graphCombIncGammaStatPlotScaled[ncent]);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSysScaled[ncent], markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent],widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlotScaled[ncent], markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent],widthLinesBoxes);
        graphCombIncGammaSysScaled[ncent]->Draw("E2same");
        graphCombIncGammaStatPlotScaled[ncent]->Draw("Epsame");
        legendYieldIncGamma->AddEntry(graphCombIncGammaSysScaled[ncent],Form("#gamma_{inc}, %s (#times 10^{%d})",multbins[ncent].Data(),3-ncent),"pf");
        fitHagGammaCombScaled[ncent] = ScaleTF1(fitHagGammaComb[ncent],pow(10,3-ncent),Form("fitHagGammaCombScaled_%s",multbinsPHOS[ncent].Data()));
        DrawGammaSetMarkerTF1( fitHagGammaCombScaled[ncent], 7, 2, colorCombpPb);
        fitHagGammaCombScaled[ncent]->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        // legendYieldIncGamma->AddEntry(fitHagGammaComb[ncent],"Hagedorn fit","l");
        fitHagGammaCombScaled[ncent]->Draw("same");
        // DrawGammaSetMarkerTF1( fitTCMGammaComb[ncent], 9, 2, kGray+1);
        // legendYieldIncGamma->AddEntry(fitTCMGammaComb[ncent],"TCM fit","l");
        // fitTCMGammaComb[ncent]->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        // fitTCMGammaComb[ncent]->Draw("same");
        // DrawGammaSetMarkerTF1( fitTsallisGammaComb[ncent], 5, 2, colorCombpPbBox);
        // legendYieldIncGamma->AddEntry(fitTsallisGammaComb[ncent],"Tsalis fit","l");
        // fitTsallisGammaComb[ncent]->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        // fitTsallisGammaComb[ncent]->Draw("same");
    }
        legendYieldIncGamma->Draw();
        TLatex *labelEnergyInvYieldPaperAll = new TLatex(0.20, 0.20+0.04*3, Form("%s",collisionSystempPb.Data()));
        TLatex *labelEnergyInvYieldPaperAll4 = new TLatex(0.20, 0.20+0.04*4, Form("%s",collisionSystempPb.Data()));
        SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        SetStyleTLatex( labelEnergyInvYieldPaperAll4, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelEnergyInvYieldPaperAll->Draw();
        TLatex *labelALICEInvYieldPaperAll  = new TLatex(0.20,0.20+0.04*4,textALICE.Data());
        SetStyleTLatex( labelALICEInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelALICEInvYieldPaperAll->Draw();
        TLatex *labelALICENormUnPaperAll    = new TLatex(0.20,0.20+0.05+0.04,"Norm. unc. 3.1%");
        SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICENormUnPaperAll->Draw();

        histo2DSpectraGamma->Draw("same,axis");

    canvasSpectraGamma->SaveAs(Form("%s/InvYield_IncGamma_AllCent.%s",outputDir.Data(), suffix.Data()));
    canvasSpectraGamma->SaveAs(Form("%s/InvYield_IncGamma_AllCent.pdf",outputDir.Data()));


    histo2DSpectraGamma->GetYaxis()->SetRangeUser(4e-12,10e6);
    // histo2DSpectraGamma->GetYaxis()->SetRangeUser(3e-10,10e5);
    histo2DSpectraGamma->Draw("copy");
    legendYieldIncGamma           = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(4*textSizeLabelsRel*0.85),textSizeLabelsPixel);
    TGraphAsymmErrors* graphCombDirGammaSpectrumSystErrScaled[5];
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErrPlotScaled[5];
    TGraphAsymmErrors* graphCombDirGammaSpectrumSumErrArScaled[5];
    for(Int_t ncent = 0; ncent < 4; ncent++){

        if (graphCombDirGammaSpectrumSystErr[ncent]){
            graphCombDirGammaSpectrumSystErrScaled[ncent]    = ScaleGraph(graphCombDirGammaSpectrumSystErr[ncent],pow(10,6-(ncent*2)));
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErrScaled[ncent], markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent], widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErrScaled[ncent]->Draw("E2same");
            legendYieldIncGamma->AddEntry(graphCombDirGammaSpectrumSystErrScaled[ncent],Form("#gamma_{dir}, %s (#times 10^{%d})",multbins[ncent].Data(),6-(ncent*2)),"pf");
        }
        if (graphCombDirGammaSpectrumStatErrPlot[ncent]){
            graphCombDirGammaSpectrumStatErrPlotScaled[ncent]    = ScaleGraph(graphCombDirGammaSpectrumStatErrPlot[ncent],pow(10,6-(ncent*2)));
            ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErrPlotScaled[ncent]);
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlotScaled[ncent], markerStyleMult[ncent], markerSizeMult[ncent], colorMult[ncent] , colorMult[ncent]);
            graphCombDirGammaSpectrumStatErrPlotScaled[ncent]->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr[ncent]){
            graphCombDirGammaSpectrumSumErrArScaled[ncent]    = ScaleGraph(graphCombDirGammaSpectrumSumErrAr[ncent],pow(10,6-(ncent*2)));
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrArScaled[ncent] , 1, 3,  colorMult[ncent] , colorMult[ncent], 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErrArScaled[ncent]->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrArScaled[ncent]);
        }
    }
        legendYieldIncGamma->Draw();
        TLatex *labelEnergyInvYieldPaperAllTop = new TLatex(0.945, 0.76+0.04*3, Form("%s",collisionSystempPb.Data()));
        SetStyleTLatex( labelEnergyInvYieldPaperAllTop, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
        labelEnergyInvYieldPaperAllTop->Draw();
        TLatex *labelALICEInvYieldPaperAllTop  = new TLatex(0.945,0.76+0.04*4,textALICE.Data());
        SetStyleTLatex( labelALICEInvYieldPaperAllTop, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
        labelALICEInvYieldPaperAllTop->Draw();
        TLatex *labelALICENormUnPaperAllTop    = new TLatex(0.945,0.76+0.05+0.04,"Norm. unc. 3.1%");
        SetStyleTLatex( labelALICENormUnPaperAllTop, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 31);
        labelALICENormUnPaperAllTop->Draw();

        histo2DSpectraGamma->Draw("same,axis");

    canvasSpectraGamma->SaveAs(Form("%s/InvYield_DirGamma_AllCent.%s",outputDir.Data(), suffix.Data()));
    canvasSpectraGamma->SaveAs(Form("%s/InvYield_DirGamma_AllCent.pdf",outputDir.Data()));

    //*******************************************************************************************************************************************
    //******************************************* DR plot with pp and PbPb measurements *********************************************************
    //*******************************************************************************************************************************************

    textSizeLabelsPixel = 43;

    Double_t arrayBoundsXIndMeasRatio[2];
    Double_t arrayBoundsYIndMeasRatio[5];
    Double_t relativeMarginsIndMeasRatioX[3];
    Double_t relativeMarginsIndMeasRatioY[3];
    ReturnCorrectValuesForCanvasScaling(1200, 1400, 1, 4, 0.09, 0.01, 0.01, 0.065, arrayBoundsXIndMeasRatio, arrayBoundsYIndMeasRatio, relativeMarginsIndMeasRatioX, relativeMarginsIndMeasRatioY);
    Double_t margin = relativeMarginsIndMeasRatioX[0]*1200;

    TCanvas * canvasRatioIndDR = new TCanvas("canvasRatioIndDR","",10,10,1200,1800);  // gives the page size
    canvasRatioIndDR->cd();
    TPad* padPartRatio[4];
    Double_t textsizeLabelsPad[4];
    Double_t textsizeFacPad[4];
    TH2D *dummyDRCentPlot[4];
    TLegend* legendDRCentPlot;
    TLatex* labelDRCentPlot[4];
    TLatex* labelDRCentPlotALICE;

    for(Int_t ncent = 0; ncent < 4; ncent++){
        cout << arrayBoundsYIndMeasRatio[ncent+1] << endl;
        padPartRatio[ncent] = new TPad(Form("padPartRatio%d",ncent), "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[ncent+1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[ncent],-1, -1, -2);
        DrawGammaPadSettings( padPartRatio[ncent], relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[ncent==0 ? 0 : 1], relativeMarginsIndMeasRatioY[ncent==3 ? 2 : 1]);
        padPartRatio[ncent]->Draw();
    }
    for(Int_t ncent = 0; ncent < 4; ncent++){
        if (padPartRatio[ncent]->XtoPixel(padPartRatio[ncent]->GetX2()) < padPartRatio[ncent]->YtoPixel(padPartRatio[ncent]->GetY1())){
            textsizeLabelsPad[ncent] = (Double_t)textSizeLabelsPixel/padPartRatio[ncent]->XtoPixel(padPartRatio[ncent]->GetX2()) ;
            textsizeFacPad[ncent] = (Double_t)1./padPartRatio[ncent]->XtoPixel(padPartRatio[ncent]->GetX2()) ;
        } else {
            textsizeLabelsPad[ncent] = (Double_t)textSizeLabelsPixel/padPartRatio[ncent]->YtoPixel(padPartRatio[ncent]->GetY1());
            textsizeFacPad[ncent] = (Double_t)1./padPartRatio[ncent]->YtoPixel(padPartRatio[ncent]->GetY1());
        }
        dummyDRCentPlot[ncent] = new TH2D(Form("dummyDRCentPlot[%d]",ncent), Form("dummyDRCentPlot[%d]",ncent), 1000, doubleRatioX[0], doubleRatioX[1], 1000., doubleRatio[0], doubleRatio[1]);
        SetStyleHistoTH2ForGraphs( dummyDRCentPlot[ncent], "#it{p}_{T} (GeV/#it{c})", ncent==2 ? "#it{R}_{#gamma}" : "",
                                0.85*textsizeLabelsPad[ncent], textsizeLabelsPad[ncent], 0.85*textsizeLabelsPad[ncent], textsizeLabelsPad[ncent], 0.95,0.10/(textsizeFacPad[ncent]*margin), 510, 505);
        dummyDRCentPlot[ncent]->GetXaxis()->SetLabelOffset(-0.015);
        dummyDRCentPlot[ncent]->GetXaxis()->SetTickLength(ncent==3 ? 0.05 : 0.055);
        dummyDRCentPlot[ncent]->GetYaxis()->SetTickLength(ncent==3 ? 0.032 : 0.025);

        padPartRatio[ncent]->cd();
        padPartRatio[ncent]->SetLogx(1);
        dummyDRCentPlot[ncent]->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        if(ncent==1)
            legendDRCentPlot = GetAndSetLegend2(0.12, (0.92-3*textsizeLabelsPad[ncent]), 0.12+0.42,  0.94-textsizeLabelsPad[ncent],textsizeLabelsPad[ncent], 2,"", 42, 0.3);

        for(Int_t i = 0; i < 11; i++){
            if(graphDRNonFitSysErr[ncent][i] && histoDRNonFitStatErr[ncent][i]){
                graphDRNonFitSysErr[ncent][i]->Draw("E2same");
                histoDRNonFitStatErr[ncent][i]->Draw("p,same,e0,X0");
                if(ncent==1)
                    legendDRCentPlot->AddEntry(graphDRNonFitSysErr[ncent][i],Form("%s", nameMeasGlobalLabel[i].Data()),"pf");

            }
        }
        labelDRCentPlot[ncent] = new TLatex(0.12, ncent==0 ? 0.82 : (ncent==3 ? 0.88 : 0.86 ), Form("%s %s", multbins[ncent].Data(), collisionSystempPb.Data()));
        SetStyleTLatex( labelDRCentPlot[ncent], textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelDRCentPlot[ncent]->Draw();
        if(ncent==0){
            labelDRCentPlotALICE = new TLatex(0.12, 0.82-textsizeLabelsPad[ncent] , Form("%s", textALICE.Data()));
            SetStyleTLatex( labelDRCentPlotALICE, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
            labelDRCentPlotALICE->Draw();
        }
        if(ncent==1)
            legendDRCentPlot->Draw();
    }


    canvasRatioIndDR->SaveAs(Form("%s/DR_IndMeas_pPb5TeVRun1_allCent.%s", outputDir.Data(), suffix.Data()));
    canvasRatioIndDR->SaveAs(Form("%s/DR_IndMeas_pPb5TeVRun1_allCent.pdf", outputDir.Data()));
    



    for(Int_t ncent = 0; ncent < 4; ncent++){
        padPartRatio[ncent]->cd();
        padPartRatio[ncent]->SetLogx(1);
        dummyDRCentPlot[ncent]->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        if(graphCombDRSys[ncent] && graphCombDRStat[ncent]){
            graphCombDRSys[ncent]->Draw("E2same");
            graphCombDRStatPlot[ncent]->Draw("p,same,z");

        }
        labelDRCentPlot[ncent]->Draw();
        if(ncent==0){
            labelDRCentPlotALICE->Draw();
        }
    }


    canvasRatioIndDR->SaveAs(Form("%s/DR_Comb_pPb5TeVRun1_allCent.%s", outputDir.Data(), suffix.Data()));
    canvasRatioIndDR->SaveAs(Form("%s/DR_Comb_pPb5TeVRun1_allCent.pdf", outputDir.Data()));
    


    for(Int_t ncent = 0; ncent < 4; ncent++){

        padPartRatio[ncent]->cd();
        padPartRatio[ncent]->SetLogx(1);
        // if(ncent==2)
            dummyDRCentPlot[ncent]->GetYaxis()->SetTitle("Data/Fit");
            dummyDRCentPlot[ncent]->GetYaxis()->SetRangeUser( doubleRatio[0], 1.55);
        dummyDRCentPlot[ncent]->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        if(ncent==1)
            legendDRCentPlot = GetAndSetLegend2(0.12, (0.92-3*textsizeLabelsPad[ncent]), 0.12+0.42,  0.94-textsizeLabelsPad[ncent],textsizeLabelsPad[ncent], 2,"", 42, 0.3);

        for (Int_t i = 0; i < 10 ; i++){
            if (graphRatioGammaIndCombFitSys[ncent][i]){ //
                graphRatioGammaIndCombFitSys[ncent][i]->Draw("E2same");
                if(ncent==1)
                    legendDRCentPlot->AddEntry(graphRatioGammaIndCombFitSys[ncent][i],Form("%s", nameMeasGlobalLabel[i].Data()),"pf");
            }
            if (graphRatioGammaIndCombFitStat[ncent][i]){ //
                graphRatioGammaIndCombFitStat[ncent][i]->Draw("p,same,z");
            }
        }

        labelDRCentPlot[ncent]->Draw();
        if(ncent==0){
            labelDRCentPlotALICE->Draw();
            TLatex* labelIncGammaAllCent = new TLatex(0.95, 0.82 , "#gamma_{inc}");
            SetStyleTLatex( labelIncGammaAllCent, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
            labelIncGammaAllCent->Draw();
        }
        if(ncent==1)
            legendDRCentPlot->Draw();
    }


    canvasRatioIndDR->SaveAs(Form("%s/IncGamma_IndMeas_pPb5TeVRun1_allCent.%s", outputDir.Data(), suffix.Data()));
    canvasRatioIndDR->SaveAs(Form("%s/IncGamma_IndMeas_pPb5TeVRun1_allCent.pdf", outputDir.Data()));


    for(Int_t ncent = 0; ncent < 4; ncent++){

        padPartRatio[ncent]->cd();
        padPartRatio[ncent]->SetLogx(1);
        dummyDRCentPlot[ncent]->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        if (graphRatioGammaCombCombFitSys[ncent])
            graphRatioGammaCombCombFitSys[ncent]->Draw("E2same");
        if (graphRatioGammaCombCombFitStat[ncent])
            graphRatioGammaCombCombFitStat[ncent]->Draw("p,same,z");

        labelDRCentPlot[ncent]->Draw();
        if(ncent==0){
            labelDRCentPlotALICE->Draw();
            TLatex* labelIncGammaAllCent = new TLatex(0.95, 0.82 , "#gamma_{inc}");
            SetStyleTLatex( labelIncGammaAllCent, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
            labelIncGammaAllCent->Draw();
        }
    }


    canvasRatioIndDR->SaveAs(Form("%s/IncGamma_Comb_pPb5TeVRun1_allCent.%s", outputDir.Data(), suffix.Data()));
    canvasRatioIndDR->SaveAs(Form("%s/IncGamma_Comb_pPb5TeVRun1_allCent.pdf", outputDir.Data()));
}