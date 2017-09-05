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
    doubleRatioXpp[0]   = 0.23;     doubleRatioXpp[1]   = 24;

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
    TString collisionSystempPbNSD                   = "NSD p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    TString textALICE                               = "ALICE this thesis";
    cout << "Setting Gamma binning" << endl;
    Double_t xPtLimitsGamma[100]                    = {0};
    Int_t maxNBinsGamma                             = GetBinning( xPtLimitsGamma, "Gamma", "pPb_5.023TeV", 20 );
    for (Int_t i = 0; i< maxNBinsGamma; i++){
        cout << i << ": "<< xPtLimitsGamma[i] <<" - " << xPtLimitsGamma[i+1]<< ", " <<endl;
    }

    TH1D* histoDRPi0FitStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors*  graphDRPi0FitSysErr[11]     = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoIncGammaRatioStatErr[11]             = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors* graphIncGammaRatioSysErr[11] = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoIncGammaStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors* graphIncGammaSysErr[11]      = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoRawGamma[11]                         = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoConvProb[11]                         = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoEffi[11]                             = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
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
        histoIncGammaRatioStatErr[0]                    = (TH1D*) directoryPCMGammapPb->Get("IncRatioStatError");
        graphIncGammaRatioSysErr[0]                     = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncRatioSystError");
        histoConvProb[0]                                = (TH1D*) directoryPCMGammapPb->Get("GammaConversionProbability");
        histoEffi[0]                                    = (TH1D*) directoryPCMGammapPb->Get("GammaRecoEfficiency");
        histoPurity[0]                                  = (TH1D*) directoryPCMGammapPb->Get("GammaTruePurity");
        histoIncGammaStatErr[0]                         = (TH1D*) directoryPCMGammapPb->Get("IncGammaStatError");
        graphIncGammaSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncGammaSystError");
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
        histoIncGammaRatioStatErr[4]                    = (TH1D*) directoryPCMEMCGammapPb->Get("IncRatioStatError");
        graphIncGammaRatioSysErr[4]                     = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("IncRatioSystError");
        histoConvProb[4]                                = (TH1D*) directoryPCMEMCGammapPb->Get("GammaConversionProbability");
        histoEffi[4]                                    = (TH1D*) directoryPCMEMCGammapPb->Get("GammaRecoEfficiency");
        histoPurity[4]                                  = (TH1D*) directoryPCMEMCGammapPb->Get("GammaTruePurity");
        histoIncGammaStatErr[4]                         = (TH1D*) directoryPCMEMCGammapPb->Get("IncGammaStatError");
        graphIncGammaSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapPb->Get("IncGammaSystError");
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
//         if (histoDRPi0FitStatErr[2]){
//             histoDRPi0FitStatErr[2]->GetXaxis()->SetRangeUser(2.5,14);
//             for (Int_t i = 1; i < histoDRPi0FitStatErr[2]->GetXaxis()->FindBin(2.4); i++)
//                 histoDRPi0FitStatErr[2]->SetBinContent(i, 0);
//
//         }
        graphDRPi0FitSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("DoubleRatioPi0FitSystError");
//         if (graphDRPi0FitSysErr[2]) {
//             while (graphDRPi0FitSysErr[2]->GetX()[0] < 2.5) graphDRPi0FitSysErr[2]->RemovePoint(0);
//             while (graphDRPi0FitSysErr[2]->GetX()[graphDRPi0FitSysErr[2]->GetN()-1] > 14) graphDRPi0FitSysErr[2]->RemovePoint(graphDRPi0FitSysErr[2]->GetN()-1);
//         }
        histoIncGammaRatioStatErr[2]                    = (TH1D*) directoryEMCGammapPb->Get("IncRatioStatError");
        graphIncGammaRatioSysErr[2]                     = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("IncRatioSystError");
        histoEffi[2]                                    = (TH1D*) directoryEMCGammapPb->Get("GammaRecoEfficiency");
        histoPurity[2]                                  = (TH1D*) directoryEMCGammapPb->Get("GammaTruePurity");
        histoIncGammaStatErr[2]                         = (TH1D*) directoryEMCGammapPb->Get("IncGammaStatError");
//         if (histoIncGammaStatErr[2]){
//             histoIncGammaStatErr[2]->GetXaxis()->SetRangeUser(1.6,14);
//             for (Int_t i = 1; i < histoIncGammaStatErr[2]->GetXaxis()->FindBin(2.4); i++)
//                 histoIncGammaStatErr[2]->SetBinContent(i, 0);
//
//         }
        graphIncGammaSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("IncGammaSystError");
//         if (graphIncGammaSysErr[2]) {
//             while (graphIncGammaSysErr[2]->GetX()[0] < 2.5) graphIncGammaSysErr[2]->RemovePoint(0);
//             while (graphIncGammaSysErr[2]->GetX()[graphIncGammaSysErr[2]->GetN()-1] > 14) graphIncGammaSysErr[2]->RemovePoint(graphIncGammaSysErr[2]->GetN()-1);
//         }
        histoRawGamma[2]                                = (TH1D*) directoryEMCGammapPb->Get("GammaRawYields");
        histoEffSecCorr[0][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_K0s");
        histoEffSecCorr[1][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_K0l");
        histoEffSecCorr[2][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_Lambda");
        histoEffSecCorr[3][2]                           = (TH1D*) directoryEMCGammapPb->Get("GammaEffectiveSecondaryCorr_Rest");

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
    TFile* fileTheory                               = new TFile( fileNameTheorypPb.Data());
        TGraphAsymmErrors* graphTheoryNLODRpPb          = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10_ALICECocktail");
        TGraph* graphTheoryNLODRpPbCenter               = (TGraph*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10_ALICECocktail_Center");
        TGraphAsymmErrors* graphTheoryNLOpPb            = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphDirectPhotonNLOVogelsangInvYieldINT7_pPb5TeV_CT10");
        TGraphAsymmErrors* graphTheoryMCGillDRpPb       = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonSpecMcGill5023GeV_ALICECocktail");
        TGraph* graphTheoryMCGillDRpPbCenter            = (TGraph*) fileTheory->Get("pPb_5.023TeV/graphRGammaDirectPhotonSpecMcGill5023GeV_ALICECocktail_Center");
        TGraphAsymmErrors* graphTheoryMCGillpPb         = (TGraphAsymmErrors*) fileTheory->Get("pPb_5.023TeV/graphDirectPhotonSpecMcGill5023GeV");


    //*******************************************************************************************************************************************
    //*********************************************** Combining Rgamma ratios  ******************************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsGamma[11]          = { 3,  8,  3,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetsGammaSys[11]       = { 4,  8,  11,  0,  8,
                                        0,  0,  0,  0,  0,
                                        0};


    TGraphAsymmErrors* statErrorRelCollectionDR[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionDR[i]        = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (histoDRPi0FitStatErr[i]){
            statErrorRelCollectionDR[i]    = new TGraphAsymmErrors(histoDRPi0FitStatErr[i]);
            while (statErrorRelCollectionDR[i]->GetY()[0] == 0) statErrorRelCollectionDR[i]->RemovePoint(0);
            while (statErrorRelCollectionDR[i]->GetY()[statErrorRelCollectionDR[i]->GetN()-1] == 0) statErrorRelCollectionDR[i]->RemovePoint(statErrorRelCollectionDR[i]->GetN()-1);
            statErrorRelCollectionDR[i]    = CalculateRelErrUpAsymmGraph( statErrorRelCollectionDR[i], Form("relativeStatErrorDR_%s", nameMeasGlobal[i].Data()));
        }
    }

    TGraphAsymmErrors* sysErrorRelCollectionDR[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollectionDR[i]         = NULL;
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
    }

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraph* graphWeightsDR[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsDR[i]                   = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNameDROutputWeighting       = Form("%s/DR_WeightingMethod.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombDRStat      = NULL;
    TGraphAsymmErrors* graphCombDRSys       = NULL;
    TGraphAsymmErrors* graphCombDRTot       = CombinePtPointsSpectraFullCorrMat(    histoDRPi0FitStatErr,    graphDRPi0FitSysErr,
                                                                                    xPtLimitsGamma, maxNBinsGamma,
                                                                                    offSetsGamma, offSetsGammaSys,
                                                                                    graphCombDRStat, graphCombDRSys,
                                                                                    fileNameDROutputWeighting, "pPb_5.023TeV", "DR", kTRUE,
                                                                                    NULL, "" );


    if (graphCombDRTot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    while (graphCombDRStat->GetX()[0] < 0.4){
        graphCombDRStat->RemovePoint(0);
    }
    while (graphCombDRTot->GetX()[0] < 0.4){
        graphCombDRTot->RemovePoint(0);
    }
    while (graphCombDRSys->GetX()[0] < 0.4){
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
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    canvasWeights->SetLogx();

    TH2F * histo2DDRWeights;
    histo2DDRWeights = new TH2F("histo2DDRWeights","histo2DDRWeights",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DDRWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
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

    TLatex *labelWeightsEnergy      = new TLatex(0.95,0.20,collisionSystempPb.Data());
    SetStyleTLatex( labelWeightsEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelWeightsEnergy->Draw();
    TLatex *labelWeightsDR         = new TLatex(0.95,0.15,"R_{#gamma}");
    SetStyleTLatex( labelWeightsDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelWeightsDR->Draw();

    DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/DR_Weights.%s",outputDir.Data(),suffix.Data()));

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

    TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystempPb.Data());
    SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrEnergy->Draw();
    TLatex *labelRelSysErrDR       = new TLatex(0.15,0.85,"R_{#gamma}");
    SetStyleTLatex( labelRelSysErrDR, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DR_RelSysErr.%s",outputDir.Data(),suffix.Data()));


    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();

    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0,50.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetYaxis()->SetRangeUser(0,24.5);
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

    TLatex *labelRelStatErrEnergy   = new TLatex(0.95,0.89,collisionSystempPb.Data());
    SetStyleTLatex( labelRelStatErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrEnergy->Draw();
    TLatex *labelRelStatErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
    SetStyleTLatex( labelRelStatErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrDR->Draw();

    canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr.%s",outputDir.Data(),suffix.Data()));


    //*******************************************************************************************************************************************
    //*********************************************** Combining Rgamma ratios  ******************************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsIncGamma[11]       = { 3,  8,  3,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetsIncGammaSys[11]    = { 4,  8,  11,  0,  8,
                                        0,  0,  0,  0,  0,
                                        0};
//     Int_t offSetGammaShifting[11]   = { 0,  6,  8,  2,  5,
//                                         3,  0,  0,  0,  21,
//                                         0 };
//     Int_t nComBinsGammaShifting[11] = { 30, 31, 30, 0,  31,
//                                         17,  0,  0,  0,  0,
//                                         0 };

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
                                                                                            offSetsGamma, offSetsGammaSys,
                                                                                            graphCombIncGammaStat, graphCombIncGammaSys,
                                                                                            fileNameIncGammaOutputWeighting, "pPb_5.023TeV", "IncGamma", kTRUE,
                                                                                            NULL, "" );


    if (graphCombIncGammaTot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    while (graphCombIncGammaStat->GetX()[0] < 0.4){
        graphCombIncGammaStat->RemovePoint(0);
    }
    while (graphCombIncGammaTot->GetX()[0] < 0.4){
        graphCombIncGammaTot->RemovePoint(0);
    }
    while (graphCombIncGammaSys->GetX()[0] < 0.4){
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
    Int_t nMeasSetIncGamma               = 3;
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

    TH2F * histo2DIncGammaWeights;
    histo2DIncGammaWeights = new TH2F("histo2DIncGammaWeights","histo2DIncGammaWeights",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DIncGammaWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
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

    DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/IncGamma_Weights.%s",outputDir.Data(),suffix.Data()));

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


    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    canvasRelStatErr->cd();

    histo2DRelStatErr->GetYaxis()->SetRangeUser(-0.2,7.5);
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

    //*******************************************************************************************************************************************
    //************************************************* Fitting gamma spectrum ******************************************************************
    //*******************************************************************************************************************************************
    TF1* fitHagGammaComb            = FitObject("h","fitHagGammaComb","Gamma",graphCombIncGammaTot,graphCombIncGammaTot->GetX()[0],
                                                    graphCombIncGammaTot->GetX()[graphCombIncGammaTot->GetN()],NULL,"QNRME+");
    TString forOutput               = WriteParameterToFile(fitHagGammaComb);
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

    graphCombIncGammaTot            = ApplyXshift(graphCombIncGammaTot, fitShiftingGamma);
    cout << "comb" << endl;
    graphCombIncGammaStat->Print();
    graphCombIncGammaStat           = ApplyXshiftIndividualSpectra( graphCombIncGammaTot,
                                                                    graphCombIncGammaStat,
                                                                    fitShiftingGamma,
                                                                    0, graphCombIncGammaStat->GetN());
    graphCombIncGammaSys            = ApplyXshiftIndividualSpectra( graphCombIncGammaTot,
                                                                    graphCombIncGammaSys,
                                                                    fitShiftingGamma,
                                                                    0, graphCombIncGammaSys->GetN());
    Int_t offSetGammaShifting[11]   = { 0,  0,  7,  0,  4,
                                        0,  0,  0,  0,  0,
                                        0 };
    Int_t nComBinsGammaShifting[11] = { 24, 0, 17, 0,  20,
                                        0,  0,  0,  0,  0,
                                        0 };

    for (Int_t i = 0; i< 11; i++){
        if (graphIndGammaIncStat[i]){
            cout << "shiting stat err of " << nameMeasGlobalLabel[i].Data();
            graphIndGammaIncStat[i]  = ApplyXshiftIndividualSpectra(    graphCombIncGammaTot,
                                                                        graphIndGammaIncStat[i],
                                                                        fitShiftingGamma,
                                                                        offSetGammaShifting[i], nComBinsGammaShifting[i]);

        }
        if (graphIndGammaIncSys[i]){
            cout << "shiting sys err of " << nameMeasGlobalLabel[i].Data();
            graphIndGammaIncSys[i]   = ApplyXshiftIndividualSpectra(    graphCombIncGammaTot,
                                                                        graphIndGammaIncSys[i],
                                                                        fitShiftingGamma,
                                                                        offSetGammaShifting[i], nComBinsGammaShifting[i]);
        }
    }

    TGraphAsymmErrors* graphRatioGammaIndCombFitStat[11] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioGammaIndCombFitSys[11]  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    TGraphAsymmErrors* graphRatioGammaCombCombFitTot     = (TGraphAsymmErrors*)graphCombIncGammaTot->Clone();
    graphRatioGammaCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitTot, fitTsallisGammaComb);
    TGraphAsymmErrors* graphRatioGammaCombCombFitStat    = (TGraphAsymmErrors*)graphCombIncGammaStat->Clone();
    graphRatioGammaCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitStat, fitTsallisGammaComb);
    TGraphAsymmErrors* graphRatioGammaCombCombFitSys     = (TGraphAsymmErrors*)graphCombIncGammaSys->Clone();
    graphRatioGammaCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitSys, fitTsallisGammaComb);

    for (Int_t i= 0; i< 11; i++){
        if (graphIndGammaIncStat[i]){
            graphRatioGammaIndCombFitStat[i]             = (TGraphAsymmErrors*)graphIndGammaIncStat[i]->Clone(Form("RatioGamma%sToCombFitStat", nameMeasGlobalLabel[i].Data()));
            graphRatioGammaIndCombFitStat[i]             = CalculateGraphErrRatioToFit(graphRatioGammaIndCombFitStat[i], fitTsallisGammaComb);
        }
        if (graphIndGammaIncSys[i]){
            graphRatioGammaIndCombFitSys[i]              = (TGraphAsymmErrors*)graphIndGammaIncSys[i]->Clone(Form("RatioGamma%sToCombFitSyst", nameMeasGlobalLabel[i].Data()));
            graphRatioGammaIndCombFitSys[i]              = CalculateGraphErrRatioToFit(graphRatioGammaIndCombFitSys[i], fitTsallisGammaComb);
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
    histo2DGammaRatioToCombFit               = new TH2F("histo2DGammaRatioToCombFit","histo2DGammaRatioToCombFit",1000,0.23, 25.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DGammaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPPb, textsizeLabelsPPb,
                              0.85*textsizeLabelsPPb,textsizeLabelsPPb, 0.9, 0.65, 510, 505);
    histo2DGammaRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DGammaRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
    //  histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(0.59,1.42);
    histo2DGammaRatioToCombFit->Draw("copy");

    ProduceGraphAsymmWithoutXErrors(graphRatioGammaCombCombFitStat);

    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitSys, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE);
    graphRatioGammaCombCombFitSys->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitStat, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb);
    graphRatioGammaCombCombFitStat->Draw("p,same,z");

    DrawGammaLines(0.23, 25. , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0.23, 25. , 1.05, 1.05,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.95, 0.95,0.1, kGray, 7);

    TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, collisionSystempPbNSD.Data());
    SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergy->Draw();
    TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICE.Data());
    SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitALICE->Draw();
    TLatex *labelRatioToFitGamma      = new TLatex(0.12, 0.92, "#gamma_{inc}");
    SetStyleTLatex( labelRatioToFitGamma, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
    labelRatioToFitGamma->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_PPb5023GeV.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DGammaRatioToCombFit->Draw("copy");

    for (Int_t i = 10; i > -1 ; i--){
        if (graphRatioGammaIndCombFitSys[i]){
            DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
            graphRatioGammaIndCombFitSys[i]->Draw("E2same");
        }
        if (graphRatioGammaIndCombFitStat[i]){
            ProduceGraphAsymmWithoutXErrors(graphRatioGammaIndCombFitStat[i]);
            DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitStat[i], markerStyleDet[i] ,markerSizeDet[i], colorDet[i], colorDet[i]);
            graphRatioGammaIndCombFitStat[i]->Draw("p,same,z");
        }
    }
    graphRatioGammaIndCombFitStat[4]->Draw("p,same,z");

    DrawGammaLines(0.23, 25. , 1., 1.,0.5, kGray+2);
    DrawGammaLines(0.23, 25. , 1.05, 1.05,0.5, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.95, 0.95,0.5, kGray, 7);
    DrawGammaLines(0.23, 25. , 1.1, 1.1,0.5, kGray, 9);
    DrawGammaLines(0.23, 25. , 0.9, 0.9,0.5, kGray, 9);

    labelRatioToFitEnergy->Draw();
    labelRatioToFitALICE->Draw();
    labelRatioToFitGamma->Draw();
    histo2DGammaRatioToCombFit->Draw("same,axis");


    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    Double_t rowsLegendGammaRatio[4]          = {0.26, 0.21, 0.16, 0.16};
    Double_t rowsLegendGammaRatioAbs[4]       = {0.92, 0.69, 0.64, 0.47 };
    Double_t columnsLegendGammaRatio[6]       = {0.14, 0.26, 0.36, 0.48, 0.7, 0.8};
    Double_t columnsLegendGammaRatioAbs[6]    = {0.215, 0.69, 1.15, 2, 6.6, 9.6};
    Double_t lengthBox                          = 0.2;
    Double_t heightBox                          = 0.03/2;
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

    TMarker* markerPCMGammaRatio           = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[0],columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[1],1);
    markerPCMGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[1]);
    TMarker* markerEMCALGammaRatio         = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[2], columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[2],1);
    markerEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[2]);
    TMarker* markerPCMEMCALGammaRatio      = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[4], columnsLegendGammaRatio[4] ,rowsLegendGammaRatio[1],1);
    markerPCMEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[4] ,rowsLegendGammaRatioAbs[1]);

    TBox* boxPCMGammaRatio                 = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[0], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[1]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[1]+ heightBox);
    boxPCMGammaRatio->Draw("l");
    TBox* boxEMCALGammaRatio               = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[2], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[2]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[2]+ heightBox);
    boxEMCALGammaRatio->Draw("l");
    TBox* boxPCMEMCALGammaRatio            = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[4], columnsLegendGammaRatioAbs[5]-0.5*lengthBox , rowsLegendGammaRatioAbs[1]- heightBox, columnsLegendGammaRatioAbs[5]+ 18*lengthBox, rowsLegendGammaRatioAbs[1]+ heightBox);
    boxPCMEMCALGammaRatio->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit.%s",outputDir.Data(),suffix.Data()));

    //*******************************************************************************************************************************************
    //********************************************** Plotting individual Rgamma ratios **********************************************************
    //*******************************************************************************************************************************************

    // double ratio combined
    TCanvas *canvasDoubleRatio = new TCanvas("canvasDoubleRatio","",0.095,0.09,1000,815);
    DrawGammaCanvasSettings( canvasDoubleRatio, 0.086, 0.01, 0.01, 0.092);
    canvasDoubleRatio->cd();
    canvasDoubleRatio->SetLogx();

    Double_t minY                            = 0.85;
    Double_t maxY                            = 1.65;
    Double_t textSizeSinglePad               = 0.05;
    TH2F * hist2DDRDummySingle       = new TH2F("hist2DDRDummySingle","hist2DDRDummySingle",1000,doubleRatioXpp[0], doubleRatioXpp[1],1000,doubleRatio[0], doubleRatio[1]);
    SetStyleHistoTH2ForGraphs(hist2DDRDummySingle, "#it{p}_{T} (GeV/#it{c})","#it{R}_{#gamma}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.8,0.8);
    hist2DDRDummySingle->GetXaxis()->SetLabelOffset(-0.01);
    hist2DDRDummySingle->GetXaxis()->SetMoreLogLabels(kTRUE);
    hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRpPbSingle = GetAndSetLegend2(0.12,0.90-1.05*0.85*textSizeSinglePad*2,0.5,0.90, 0.85*textSizeSinglePad, 2, "", 42, 0.25);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 10; i > -1; i--){
            if (graphDRPi0FitSysErr[i]){
                DrawGammaSetMarkerTGraphAsym(graphDRPi0FitSysErr[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphDRPi0FitSysErr[i]->Draw("E2same");
                legendDRpPbSingle->AddEntry(graphDRPi0FitSysErr[i],nameMeasGlobalLabel[i],"pf");

            }
            if (histoDRPi0FitStatErr[i]){
                DrawGammaSetMarker(histoDRPi0FitStatErr[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoDRPi0FitStatErr[i]->Draw("p,same,e0,X0");
                if (!graphDRPi0FitSysErr[i])legendDRpPbSingle->AddEntry(histoDRPi0FitStatErr[i],nameMeasGlobalLabel[i],"p");
            }
        }

        if (histoDRPi0FitStatErr[2]) histoDRPi0FitStatErr[2]->Draw("p,same,e0,X0");
        if (histoDRPi0FitStatErr[0]) histoDRPi0FitStatErr[0]->Draw("p,same,e0,X0");
        if (histoDRPi0FitStatErr[4]) histoDRPi0FitStatErr[4]->Draw("p,same,e0,X0");
        legendDRpPbSingle->Draw();


        TLatex *labelALICEDRSingle = new TLatex(0.95,0.91,textALICE.Data());
        SetStyleTLatex( labelALICEDRSingle, 0.85*textSizeSinglePad,4, 1, 42, kTRUE, 31);
        labelALICEDRSingle->Draw();

        TLatex *labelDRCentpPbSingle = new TLatex(0.12,0.91,collisionSystempPbNSD.Data());
        SetStyleTLatex( labelDRCentpPbSingle, 0.85*textSizeSinglePad,4, 1, 42, kTRUE, 11);
        labelDRCentpPbSingle->Draw();
    hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pPb.%s", outputDir.Data(), suffix.Data()));

    hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRPCMNLOpPbComb = GetAndSetLegend2(0.12,0.90-1.05*0.85*textSizeSinglePad*3,0.5,0.90, 0.85*textSizeSinglePad, 2, "", 42, 0.25);
        if (graphTheoryNLODRpPb) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryNLODRpPb, 0, 0, kAzure-9, kAzure-9, 0.2, kTRUE, kAzure-9);
            graphTheoryNLODRpPb->Draw("3,same");
            legendDRPCMNLOpPbComb->AddEntry(graphTheoryNLODRpPb,"NLO pQCD PDF: CT10 FF: GRV ","f");
            legendDRPCMNLOpPbComb->AddEntry((TObject*)0,"","");
        }
        if (graphTheoryNLODRpPbCenter){
            DrawGammaNLOTGraph( graphTheoryNLODRpPbCenter, 2, 7, kAzure+2);
            graphTheoryNLODRpPbCenter->Draw("lc,same");
        }
        if (graphTheoryMCGillDRpPb) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryMCGillDRpPb, 0, 0, colorNLOMcGill, colorNLOMcGill, 0.2, kTRUE, colorNLOMcGill, kTRUE);
            graphTheoryMCGillDRpPb->Draw("3,same");
            legendDRPCMNLOpPbComb->AddEntry(graphTheoryMCGillDRpPb,"McGill","f");
            legendDRPCMNLOpPbComb->AddEntry((TObject*)0,"","");
        }
        if (graphTheoryMCGillDRpPbCenter){
            DrawGammaNLOTGraph( graphTheoryMCGillDRpPbCenter, 2, 7, colorNLOMcGill);
            graphTheoryMCGillDRpPbCenter->Draw("lc,same");
        }

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        DrawGammaSetMarkerTGraphAsym(graphCombDRSys, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRStat, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes);
        legendDRPCMNLOpPbComb->AddEntry(graphCombDRSys,"ALICE","pf");
        graphCombDRSys->Draw("E2same");
        graphCombDRStat->Draw("Epsame");
        legendDRPCMNLOpPbComb->Draw();

        labelALICEDRSingle->Draw();

        labelDRCentpPbSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pPb.%s", outputDir.Data(), suffix.Data()));


    // **********************************************************************************************************************
    // ******************************** Efficiency for gamma individual measurements ****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 55;
    Double_t textSizeLabelsRel          = 55./1200;
    cout << textSizeLabelsRel << endl;

    TCanvas* canvasEff   = new TCanvas("canvasEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEff,  0.1, 0.01, 0.015, 0.095);
//     canvasEff->SetLogy(1);
    canvasEff->SetLogx(1);

    TH1F * histo1DEff            = new TH1F("histo1DEff", "histo1DEff",1000, doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(  histo1DEff, "#it{p}_{T} (GeV/#it{c})","#it{#varepsilon}_{rec}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
    histo1DEff->GetYaxis()->SetRangeUser(0.1, 1.02 );
    histo1DEff->GetYaxis()->SetLabelOffset(0.001);
    histo1DEff->GetXaxis()->SetLabelOffset(-0.01);
    histo1DEff->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo1DEff->DrawCopy();

        TLegend* legendEffiGamma           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
        for (Int_t i = 0; i < 11; i++){
            if (histoEffi[i]){
                DrawGammaSetMarker(histoEffi[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoEffi[i]->Draw("p,same,e");
                legendEffiGamma->AddEntry(histoEffi[i],nameMeasGlobalLabel[i],"p");
            }
        }
        legendEffiGamma->Draw();

        TLatex *labelPerfEffi           = new TLatex(0.15,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
        labelPerfEffi->Draw();
        TLatex *labelEnergyEffi         = new TLatex(0.15,0.87,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
        labelEnergyEffi->Draw();

    canvasEff->Update();
    canvasEff->Print(Form("%s/Gamma_Effiency.%s",outputDir.Data(),suffix.Data()));

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

        TLegend* legendPurityGamma           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
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
        TLatex *labelEnergyPurity         = new TLatex(0.15,0.87,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyPurity, textSizeLabelsRel,4);
        labelEnergyPurity->Draw();

    canvasPurity->Update();
    canvasPurity->Print(Form("%s/Gamma_Purity.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** ConvProb for gamma individual measurements ******************************************
    // **********************************************************************************************************************
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
        TLatex *labelEnergyConvProb         = new TLatex(0.15,0.87,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyConvProb, textSizeLabelsRel,4);
        labelEnergyConvProb->Draw();

    canvasConvProb->Update();
    canvasConvProb->Print(Form("%s/Gamma_ConvProb.%s",outputDir.Data(),suffix.Data()));

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
        TLatex *labelEnergyPileUp         = new TLatex(0.15,0.87,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyPileUp, textSizeLabelsRel,4);
        labelEnergyPileUp->Draw();

    canvasPileUp->Update();
    canvasPileUp->Print(Form("%s/Gamma_PileUp.%s",outputDir.Data(),suffix.Data()));


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
    histo1DEffectiveSecCorr->GetYaxis()->SetLabelOffset(0.001);
    histo1DEffectiveSecCorr->GetXaxis()->SetLabelOffset(-0.01);
    histo1DEffectiveSecCorr->GetXaxis()->SetMoreLogLabels(kTRUE);


    Double_t minYSecCorr[4]             = {0.0, 0.0, 0.0, 0.0};
    Double_t maxYSecCorr[4]             = {0.05, 0.003, 1.0e-3, 0.025};
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
        TLatex *labelEnergyEffectiveSecCorr         = new TLatex(0.15,0.84,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyEffectiveSecCorr, textSizeLabelsRel,4);
        labelEnergyEffectiveSecCorr->Draw();

        canvasEffectiveSecCorr->Update();
        canvasEffectiveSecCorr->Print(Form("%s/Gamma_EffectiveSecCorr_%s.%s",outputDir.Data(), nameOutputSec[k].Data(), suffix.Data()));
    }

    // **********************************************************************************************************************
    // ******************************** InvYield for combined inc gamma measurement *****************************************
    // **********************************************************************************************************************
    TCanvas* canvasInvYieldGamma          = new TCanvas("canvasInvYieldGamma","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasInvYieldGamma, 0.16, 0.02, 0.02, 0.08);
    canvasInvYieldGamma->SetLogx();
    canvasInvYieldGamma->SetLogy();


    TH2F * histo2DYieldGamma              = new TH2F("histo2DYieldGamma","histo2DYieldGamma",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,7e-9,10e1);
    SetStyleHistoTH2ForGraphs(histo2DYieldGamma, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.65);
    histo2DYieldGamma->GetXaxis()->SetMoreLogLabels();
    histo2DYieldGamma->GetXaxis()->SetLabelOffset(-0.01);
    histo2DYieldGamma->Draw("copy");

        TLegend* legendYieldIncGamma           = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(3*textSizeLabelsRel*0.85),textSizeLabelsPixel);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStat, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes);
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStat->Draw("Epsame");
        legendYieldIncGamma->AddEntry(graphCombIncGammaSys,"ALICE","pf");
        legendYieldIncGamma->Draw();

        DrawGammaSetMarkerTF1( fitHagGammaComb, 7, 2, colorCombpPb);
        legendYieldIncGamma->AddEntry(fitHagGammaComb,"Hagedorn fit","l");
        fitHagGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitHagGammaComb->Draw("same");
        DrawGammaSetMarkerTF1( fitTsallisGammaComb, 5, 2, colorCombpPbBox);
        legendYieldIncGamma->AddEntry(fitTsallisGammaComb,"Tsalis fit","l");
        fitTsallisGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTsallisGammaComb->Draw("same");

        TLatex *labelEnergyInvYieldPaperAll = new TLatex(0.20, 0.20+0.04*3, collisionSystempPbNSD.Data());
        SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelEnergyInvYieldPaperAll->Draw();
        TLatex *labelALICEInvYieldPaperAll  = new TLatex(0.20,0.20+0.04*2,textALICE.Data());
        SetStyleTLatex( labelALICEInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelALICEInvYieldPaperAll->Draw();
        TLatex *labelALICENormUnPaperAll    = new TLatex(0.20,0.20+0.05*1,"Norm. unc. 3.1%");
        SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICENormUnPaperAll->Draw();

        histo2DYieldGamma->Draw("same,axis");

    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** InvYield for individual inc gamma measurements **************************************
    // **********************************************************************************************************************
    histo2DYieldGamma->Draw("copy");

        TLegend* legendYieldIncGammaInd       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(3*textSizeLabelsRel*0.85),textSizeLabelsPixel);
        for (Int_t i = 0; i < 11; i++){
            if (graphIndGammaIncSys[i]){
                DrawGammaSetMarkerTGraphAsym(graphIndGammaIncSys[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphIndGammaIncSys[i]->Draw("E2same");
                legendYieldIncGammaInd->AddEntry(graphIndGammaIncSys[i],nameMeasGlobalLabel[i],"pf");
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

}