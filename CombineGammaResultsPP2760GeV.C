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
#include "CommonHeaders/ExtractSignalBinning.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"


void CombineGammaResultsPP2760GeV(  TString inputFileNamePCM        = "",
                                    Bool_t havePHOS                 = kFALSE,
                                    TString inputFileNamePHOS       = "",
                                    Bool_t haveEMC                  = kFALSE,
                                    TString inputFileNameEMC        = "",
                                    Bool_t havePCMEMC               = kFALSE,
                                    TString inputFileNamePCMEMC     = "",
                                    TString suffix                  = "eps",
                                    TString fileNameCorrelations    = "",
                                    Bool_t enablepValueCalc         = kFALSE,
                                    TString inputFileNameCocktail   = "",
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
    TString outputDir                                           = Form("%s/%s/CombineGammaMeasurementspp2760GeV",suffix.Data(),dateForOutput.Data());
    TString fileNameTheorypp2760GeV                             = "ExternalInput/Theory/TheoryCompilationPP.root";

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCMGammapp2760GeV.root", inputFileNamePCM.Data(), outputDir.Data()));
    if (havePHOS) gSystem->Exec(Form("cp %s %s/InputPHOSGammapp2760GeV.root", inputFileNamePHOS.Data(), outputDir.Data()));
    if (havePCMEMC) gSystem->Exec(Form("cp %s %s/InputPCMEMCGammapp2760GeV.root", inputFileNamePCMEMC.Data(), outputDir.Data()));
    if (haveEMC) gSystem->Exec(Form("cp %s %s/InputEMCGammapp2760GeV.root", inputFileNameEMC.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/Theorypp2760GeV.root", fileNameTheorypp2760GeV.Data(), outputDir.Data()));

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
    doubleRatioX[0]     = 0.7;      doubleRatioX[1]     = 14.5;
    doubleRatioXpp[0]   = 0.23;     doubleRatioXpp[1]   = 14.5  ;

    Color_t colorCocktailPi0                        = kRed+2;
    Color_t colorCocktailEta                        = kBlue+1;
    Color_t colorCocktailEtaP                       = kOrange+1;
    Color_t colorCocktailOmega                      = kYellow+2;
    Color_t colorCocktailPhi                        = kViolet;
    Color_t colorCocktailRho0                       = kAzure-2;
    Color_t colorCocktailSigma0                     = kGray+1;
    Color_t  colorJETPHOX                           = kRed+2;
    Color_t  colorJETPHOXBand                       = kRed-9;

    Color_t  colorNLOWerner                         = kAzure+2;
    Color_t  colorNLOWernerBand                     = kAzure-9;
    Color_t  colorNLOMcGill                         = kGreen+2;
    Color_t  colorNLOMcGillBand                     = kGreen-6;
    Style_t  styleMarkerNLOWerner                   = 24;
    Style_t  styleLineNLOWerner                     = 5;
    Style_t  styleMarkerJETPHOX                     = 23;
    Style_t  styleLineJETPHOX                       = 8;
    Style_t  styleLineMcGill                        = 7;
    Width_t  widthLineNLO                           = 2.;

    // Cocktail names and labels
    Int_t nParticles                                = 14;
    TString fParticle[14]                           = { "Pi0", "Eta", "omega", "EtaPrim", "rho0", "rho+", "rho-", "phi", "Delta0", "Delta+", "Sigma0","K0s", "K0l", "Lambda"};
    TString fParticleLatex[14]                      = { "#pi^{0}", "#eta", "#omega", "#eta'", "#rho^{0}", "#rho^{+}", "#rho^{-}", "#phi", "#Delta^{0}", "#Delta^{+}", "#Sigma^{0}", "K^{0}_{s}","K^{0}_{l}", "#Lambda"};
    TString fParticleLatexPartialSums[14]           = { "#pi^{0}", "#eta", "#omega", "#eta'", "#rho^{0}", "#rho^{#pm}", "", "#phi", "#Delta^{0/+}", "", "#Sigma^{0}", "K^{0}_{s}","K^{0}_{l}", "#Lambda"};
    Style_t     cocktailColor[14]                   = {kRed+2,kBlue+1,kYellow+2,kOrange+1,kAzure-2,kGreen+2,kRed-2,kViolet, kBlue-3, kTeal+9,kMagenta+2,kCyan+4,kViolet+4,kAzure-4};
//     Style_t     cocktailColorPartialSums[14]        = {kRed+2,kBlue+1,kYellow+2,kOrange+1,kAzure-2,kRed-2,kViolet,kGreen+2,kOrange+3, kTeal+9,kMagenta+2,kCyan+4,kViolet+4,kAzure-4};
    Color_t     cocktailMarker[14]                  = {20,    21, 25, 24,  20,      21,      24,25, 24,      25,21,        24,     25,       20};



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

    Color_t colorComb                               = GetColorDefaultColor("2.76TeV", "", "");
    Style_t markerStyleComb                         = GetDefaultMarkerStyle("2.76TeV", "", "");
    Size_t markerSizeComb                           = GetDefaultMarkerSize("2.76TeV", "", "");
    Color_t colorCombpp2760GeV                      = kBlack; // GetColorDefaultColor("2.76TeV", "", "");
    Color_t colorCombpp2760GeVBox                   = kGray+2; //GetColorDefaultColor("2.76TeV", "", "", kTRUE);
    Style_t markerStyleCombpp2760GeV                = 20; //GetDefaultMarkerStyle("2.76TeV", "", "");
    Size_t markerSizeCombpp2760GeV                  = 1.8; //GetDefaultMarkerSize("2.76TeV", "", "");

    Width_t widthLinesBoxes                         = 1.4;
    Width_t widthCommonFit                          = 2.4;

    TString collisionSystempp2760GeV                = "pp #sqrt{#it{s}} = 2.76 TeV";
    TString textALICE                               = "ALICE";
    if (isThesis)   textALICE                       = "ALICE this thesis";
    cout << "Setting Gamma binning" << endl;
    Double_t xPtLimitsGamma[100]                    = {0};
    Int_t maxNBinsGamma                             = GetBinning( xPtLimitsGamma, "Gamma", "2.76TeV", 20 );
    for (Int_t i = 0; i< maxNBinsGamma; i++){
        cout << i << ": "<< xPtLimitsGamma[i] <<" - " << xPtLimitsGamma[i+1]<< ", " <<endl;
    }


    //*******************************************************************************************************************************************
    //*********************************************** Create histogram arrays *************************************************
    //*******************************************************************************************************************************************
    TH1D* histoDRPi0FitStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors*  graphDRPi0FitSysErr[11]     = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoDRNonFitStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors*  graphDRNonFitSysErr[11]     = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoIncGammaRatioStatErr[11]             = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors* graphIncGammaRatioSysErr[11] = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoIncGammaStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors* graphIncGammaSysErr[11]      = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoConvProb[11]                         = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoEffi[11]                             = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoEffiMCPt[11]                         = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoPurity[11]                           = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoResolCorr[11]                        = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoEffSecCorr[4][11]                    = { {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL},
                                                        {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL} };
    TH1D* histoPileupCorr[11]                       = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from 2.76TeV PCM file *************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMGammapp2760GeV                    = new TFile( inputFileNamePCM.Data());
    //________________________________________________ Load PCM 2.76TeV _________________________________________________________________________
    TDirectory* directoryPCMGammapp2760GeV          = (TDirectory*)filePCMGammapp2760GeV->Get("Gamma_pp2760GeV");
        histoDRPi0FitStatErr[0]                         = (TH1D*) directoryPCMGammapp2760GeV->Get("DoubleRatioPi0FitStatError");
        graphDRPi0FitSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("DoubleRatioPi0FitSystError");
        histoDRNonFitStatErr[0]                         = (TH1D*) directoryPCMGammapp2760GeV->Get("DoubleRatioStatError");
        graphDRNonFitSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("DoubleRatioSystError");
        histoIncGammaRatioStatErr[0]                    = (TH1D*) directoryPCMGammapp2760GeV->Get("IncRatioStatError");
        graphIncGammaRatioSysErr[0]                     = (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("IncRatioSystError");
        histoConvProb[0]                                = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaConversionProbability");
        histoEffi[0]                                    = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaRecoEfficiency");
        histoEffiMCPt[0]                                = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaRecoEfficiencyMCPt");
        histoPurity[0]                                  = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaTruePurity");
        histoResolCorr[0]                               = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaResolCorr");
        histoIncGammaStatErr[0]                         = (TH1D*) directoryPCMGammapp2760GeV->Get("IncGammaStatError");
        graphIncGammaSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("IncGammaSystError");
        histoPileupCorr[0]                              = (TH1D*) directoryPCMGammapp2760GeV->Get("PileUpCorrectionFactor");
        histoEffSecCorr[0][0]                           = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0s");
        histoEffSecCorr[1][0]                           = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0l");
        histoEffSecCorr[2][0]                           = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Lambda");
        histoEffSecCorr[3][0]                           = (TH1D*) directoryPCMGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Rest");

    //*******************************************************************************************************************************************
    //*********************************************** Load PCMEMC histograms from 2.76TeV PCM file **********************************************
    //*******************************************************************************************************************************************
    if (havePCMEMC){
        TFile* filePCMEMCGammapp2760GeV                 = new TFile( inputFileNamePCMEMC.Data());
        //________________________________________________ Load PCM-EMC 2.76TeV _________________________________________________________________________
        TDirectory* directoryPCMEMCGammapp2760GeV       = (TDirectory*)filePCMEMCGammapp2760GeV->Get("Gamma_pp2760GeV");
            histoDRPi0FitStatErr[4]                         = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("DoubleRatioPi0FitStatError");
            graphDRPi0FitSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapp2760GeV->Get("DoubleRatioPi0FitSystError");
            histoDRNonFitStatErr[4]                         = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("DoubleRatioStatError");
            graphDRNonFitSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapp2760GeV->Get("DoubleRatioSystError");
            histoIncGammaRatioStatErr[4]                    = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("IncRatioStatError");
            graphIncGammaRatioSysErr[4]                     = (TGraphAsymmErrors*) directoryPCMEMCGammapp2760GeV->Get("IncRatioSystError");
            histoConvProb[4]                                = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaConversionProbability");
            histoEffi[4]                                    = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaRecoEfficiency");
            histoEffiMCPt[4]                                = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaRecoEfficiencyMCPt");
            histoPurity[4]                                  = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaTruePurity");
            histoResolCorr[4]                               = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaResolCorr");
            histoIncGammaStatErr[4]                         = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("IncGammaStatError");
            graphIncGammaSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapp2760GeV->Get("IncGammaSystError");
            histoPileupCorr[4]                              = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("PileUpCorrectionFactor");
            if (histoPileupCorr[4]) histoPileupCorr[4]->GetXaxis()->SetRangeUser(0.8,8);
            histoEffSecCorr[0][4]                           = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0s");
            histoEffSecCorr[1][4]                           = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0l");
            histoEffSecCorr[2][4]                           = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Lambda");
            histoEffSecCorr[3][4]                           = (TH1D*) directoryPCMEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Rest");
    }
    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from 2.76TeV EMC file *************************************************
    //*******************************************************************************************************************************************
    if (haveEMC){
        TFile* fileEMCGammapp2760GeV                    = new TFile( inputFileNameEMC.Data());
        //________________________________________________ Load EMC 2.76TeV _________________________________________________________________________
        TDirectory* directoryEMCGammapp2760GeV          = (TDirectory*)fileEMCGammapp2760GeV->Get("Gamma_pp2760GeV");
            histoDRPi0FitStatErr[2]                         = (TH1D*) directoryEMCGammapp2760GeV->Get("DoubleRatioPi0FitStatError");
            graphDRPi0FitSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapp2760GeV->Get("DoubleRatioPi0FitSystError");
            histoDRNonFitStatErr[2]                         = (TH1D*) directoryEMCGammapp2760GeV->Get("DoubleRatioStatError");
            graphDRNonFitSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapp2760GeV->Get("DoubleRatioSystError");
            histoIncGammaRatioStatErr[2]                    = (TH1D*) directoryEMCGammapp2760GeV->Get("IncRatioStatError");
            graphIncGammaRatioSysErr[2]                     = (TGraphAsymmErrors*) directoryEMCGammapp2760GeV->Get("IncRatioSystError");
            histoEffi[2]                                    = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaRecoEfficiency");
            histoEffiMCPt[2]                                = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaRecoEfficiencyMCPt");
            histoPurity[2]                                  = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaTruePurity");
            histoResolCorr[2]                               = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaResolCorr");
            histoIncGammaStatErr[2]                         = (TH1D*) directoryEMCGammapp2760GeV->Get("IncGammaStatError");
            graphIncGammaSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapp2760GeV->Get("IncGammaSystError");
            histoEffSecCorr[0][2]                           = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0s");
            histoEffSecCorr[1][2]                           = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_K0l");
            histoEffSecCorr[2][2]                           = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Lambda");
            histoEffSecCorr[3][2]                           = (TH1D*) directoryEMCGammapp2760GeV->Get("GammaEffectiveSecondaryCorr_Rest");
    }
    //*******************************************************************************************************************************************
    //************************************************ Load theory curves from external input ***************************************************
    //*******************************************************************************************************************************************
    TFile* fileTheory                               = new TFile( fileNameTheorypp2760GeV.Data());
        TGraphAsymmErrors* graphTheoryNLODRpp2760GeV    = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonNLOVogelsangInvYieldINT1_pp2760GeV_CT10_ALICECocktail");
        TGraph* graphTheoryNLODRpp2760GeVCenter         = (TGraph*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonNLOVogelsangInvYieldINT1_pp2760GeV_CT10_ALICECocktail_Center");
        TGraphAsymmErrors* graphTheoryNLOpp2760GeV      = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphDirectPhotonNLOVogelsangInvYieldINT1_2760GeV");
        TGraph* graphTheoryNLODRpp2760GeVPaquettCenter  = (TGraph*) fileTheory->Get("DirectPhoton/graphDRNLOPaquett_2760GeV_ALICECocktail");
        TGraph* graphTheoryNLOpp2760GeVPaquettCenter    = (TGraph*) fileTheory->Get("DirectPhoton/graphNLOPaquett_2760GeV");
        
        TGraphAsymmErrors* graphTheoryJETPHOXDRpp2760GeV    = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonJETPHOXInvYieldINT1_pp2760GeV_ALICECocktail");
        TGraph* graphTheoryJETPHOXDRpp2760GeVCenter         = (TGraph*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonJETPHOXInvYieldINT1_pp2760GeV_ALICECocktail_Center");
        TH1D* histoTheoryJETPHOXpp2760GeV      = (TH1D*) fileTheory->Get("DirectPhoton/graphDirectPhotonJETPHOXInvYield_2760GeV");
        TGraph* graphTheoryJETPHOXpp2760GeV      = new TGraph(histoTheoryJETPHOXpp2760GeV);
        
        while(graphTheoryJETPHOXpp2760GeV->GetX()[0] < 1.5) graphTheoryJETPHOXpp2760GeV->RemovePoint(0);
        // while(graphTheoryJETPHOXDRpp2760GeV->GetX()[0] < 1.5) graphTheoryJETPHOXDRpp2760GeV->RemovePoint(0);
        while(graphTheoryJETPHOXDRpp2760GeVCenter->GetX()[0] < 1.5) graphTheoryJETPHOXDRpp2760GeVCenter->RemovePoint(0);
        
        TGraphAsymmErrors* dummyJETPHOXforLegend    = new TGraphAsymmErrors(1);
        DrawGammaSetMarkerTGraphAsym(dummyJETPHOXforLegend , 2, styleLineJETPHOX, colorJETPHOX, colorJETPHOX, 0.2, kTRUE, colorJETPHOXBand);
        dummyJETPHOXforLegend->SetLineStyle(styleLineJETPHOX);
        dummyJETPHOXforLegend->SetLineWidth(2);
        TGraphAsymmErrors* dummyVogelsangforLegend    = new TGraphAsymmErrors(1);
        DrawGammaSetMarkerTGraphAsym(dummyVogelsangforLegend , 2, styleLineNLOWerner, colorNLOWerner, colorNLOWerner, 0.2, kTRUE, colorNLOWernerBand);
        dummyVogelsangforLegend->SetLineStyle(styleLineNLOWerner);
        dummyVogelsangforLegend->SetLineWidth(2);

    //*******************************************************************************************************************************************
    //*********************************************** Combining Rgamma ratios  ******************************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsGamma[11]          = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetsGammaSys[11]       = { 1,  0,  6,  0,  3,
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
                                                                                    fileNameDROutputWeighting, "2.76TeV", "RGamma", kTRUE,
                                                                                    NULL, fileNameCorrelations );


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
    Int_t nMeasSetDR               = 3;
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
    // ******************************************* Plotting weights *********************************************************
    // **********************************************************************************************************************
    Int_t textSizeLabelsPixel                 = 900*0.04;

    TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.01, 0.01, 0.09);
    canvasWeights->SetLogx();

    TH2F * histo2DDRWeights;
    histo2DDRWeights = new TH2F("histo2DDRWeights","histo2DDRWeights",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DDRWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DDRWeights->GetXaxis()->SetMoreLogLabels();
    histo2DDRWeights->GetXaxis()->SetNoExponent();
    histo2DDRWeights->Draw("copy");

    histo2DDRWeights->Draw("copy");
    TLegend* legendWeightsDR   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetDR+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDR; i++){
        DrawGammaSetMarkerTGraph(graphWeightsDR[availableDRMeas[i]], markerStyleDet[availableDRMeas[i]], markerSizeDet[availableDRMeas[i]], colorDet[availableDRMeas[i]] , colorDet[availableDRMeas[i]]);
        graphWeightsDR[availableDRMeas[i]]->Draw("p,same,z");
        legendWeightsDR->AddEntry(graphWeightsDR[availableDRMeas[i]],nameMeasGlobalLabel[availableDRMeas[i]],"p");
    }
    legendWeightsDR->Draw();

    TLatex *labelWeightsEnergy      = new TLatex(0.95,0.20,collisionSystempp2760GeV.Data());
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
    canvasWeights->SaveAs(Form("%s/DR_Weights.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();

    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0,50.0);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,30.5);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetNoExponent();
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

    TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystempp2760GeV.Data());
    SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrEnergy->Draw();
    TLatex *labelRelSysErrDR       = new TLatex(0.15,0.85,"R_{#gamma}");
    SetStyleTLatex( labelRelSysErrDR, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DR_RelSysErr.%s",outputDir.Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/DR_RelSysErr.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();

    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0,50.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetYaxis()->SetRangeUser(0,30.5);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetNoExponent();
    histo2DRelStatErr->Draw("copy");
    TLegend* legendRelStatErr2       = GetAndSetLegend2(0.14, 0.92-(0.04*(nMeasSetDR+1)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDR; i++){
        DrawGammaSetMarkerTGraph(statErrorRelCollectionDR[availableDRMeas[i]], markerStyleDet[availableDRMeas[i]], markerSizeDet[availableDRMeas[i]], colorDet[availableDRMeas[i]],
                                 colorDet[availableDRMeas[i]]);
        statErrorRelCollectionDR[availableDRMeas[i]]->Draw("p,same,z");
        legendRelStatErr2->AddEntry(statErrorRelCollectionDR[availableDRMeas[i]],nameMeasGlobalLabel[availableDRMeas[i]],"p");
    }
    legendRelStatErr2->Draw();

    TLatex *labelRelStatErrEnergy   = new TLatex(0.95,0.89,collisionSystempp2760GeV.Data());
    SetStyleTLatex( labelRelStatErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrEnergy->Draw();
    TLatex *labelRelStatErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
    SetStyleTLatex( labelRelStatErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrDR->Draw();

    canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr.%s",outputDir.Data(),suffix.Data()));
    canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr.pdf",outputDir.Data()));


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
        histo2DRelErr->GetYaxis()->SetRangeUser(0,30.5);
        histo2DRelErr->GetXaxis()->SetMoreLogLabels();
        histo2DRelErr->GetXaxis()->SetNoExponent();
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

        TLatex *labelRelTotErrEnergy   = new TLatex(0.95,0.89,collisionSystempp2760GeV.Data());
        SetStyleTLatex( labelRelTotErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
        SetStyleTLatex( labelRelTotErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DR_Reldecomp.%s",outputDir.Data(),suffix.Data()));

    //*******************************************************************************************************************************************
    //*********************************************** Combining Rgamma ratios  ******************************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    TGraphAsymmErrors* statErrorRelCollectionDRNonFit[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionDRNonFit[i]        = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (histoDRNonFitStatErr[i]){
            statErrorRelCollectionDRNonFit[i]    = new TGraphAsymmErrors(histoDRNonFitStatErr[i]);
            while (statErrorRelCollectionDRNonFit[i]->GetY()[0] == 0) statErrorRelCollectionDRNonFit[i]->RemovePoint(0);
            while (statErrorRelCollectionDRNonFit[i]->GetY()[statErrorRelCollectionDRNonFit[i]->GetN()-1] == 0) statErrorRelCollectionDRNonFit[i]->RemovePoint(statErrorRelCollectionDRNonFit[i]->GetN()-1);
            statErrorRelCollectionDRNonFit[i]    = CalculateRelErrUpAsymmGraph( statErrorRelCollectionDRNonFit[i], Form("relativeStatErrorDRNonFit_%s", nameMeasGlobal[i].Data()));
        }
    }

    TGraphAsymmErrors* sysErrorRelCollectionDRNonFit[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollectionDRNonFit[i]         = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        cout << i << endl;
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

    TGraph* graphWeightsDRNonFit[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsDRNonFit[i]                   = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNameDRNonFitOutputWeighting       = Form("%s/DRNonFit_WeightingMethod.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombDRNonFitStat      = NULL;
    TGraphAsymmErrors* graphCombDRNonFitSys       = NULL;
    TGraphAsymmErrors* graphCombDRNonFitTot       = CombinePtPointsSpectraFullCorrMat(    histoDRNonFitStatErr,    graphDRNonFitSysErr,
                                                                                    xPtLimitsGamma, maxNBinsGamma,
                                                                                    offSetsGamma, offSetsGammaSys,
                                                                                    graphCombDRNonFitStat, graphCombDRNonFitSys,
                                                                                    fileNameDRNonFitOutputWeighting, "2.76TeV", "RGamma", kTRUE,
                                                                                    NULL, fileNameCorrelations );


    if (graphCombDRNonFitTot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    while (graphCombDRNonFitStat->GetX()[0] < 0.4){
        graphCombDRNonFitStat->RemovePoint(0);
    }
    while (graphCombDRNonFitTot->GetX()[0] < 0.4){
        graphCombDRNonFitTot->RemovePoint(0);
    }
    while (graphCombDRNonFitSys->GetX()[0] < 0.4){
        graphCombDRNonFitSys->RemovePoint(0);
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
    Int_t nMeasSetDRNonFit               = 3;
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
    // ******************************************* Plotting weights *********************************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 900*0.04;
    canvasWeights->cd();

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

    DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/DRNonFit_Weights.%s",outputDir.Data(),suffix.Data()));
    canvasWeights->SaveAs(Form("%s/DRNonFit_Weights.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelSysErr->cd();
    histo2DRelSysErr->Draw("copy");

    TLegend* legenDRNonFitelSysErr2       = GetAndSetLegend2(0.62, 0.92-(0.04*(nMeasSetDRNonFit+1)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDRNonFit; i++){
        cout << "sys\t" << nameMeasGlobalLabel[availableDRNonFitMeas[i]] << endl;
        DrawGammaSetMarkerTGraph(sysErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]], markerStyleDet[availableDRNonFitMeas[i]], markerSizeDet[availableDRNonFitMeas[i]], colorDet[availableDRNonFitMeas[i]],
                                 colorDet[availableDRNonFitMeas[i]]);
        sysErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]]->Draw("p,same,z");
        sysErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]]->Print();
        legenDRNonFitelSysErr2->AddEntry(sysErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]],nameMeasGlobalLabel[availableDRNonFitMeas[i]],"p");
    }
    legenDRNonFitelSysErr2->Draw();

    labelRelSysErrEnergy->Draw();
    labelRelSysErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DRNonFit_RelSysErr.%s",outputDir.Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/DRNonFit_RelSysErr.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    canvasRelStatErr->cd();
    
    histo2DRelStatErr->Draw("copy");
    TLegend* legenDRNonFitelStatErr2       = GetAndSetLegend2(0.14, 0.92-(0.04*(nMeasSetDRNonFit+1)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDRNonFit; i++){
        DrawGammaSetMarkerTGraph(statErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]], markerStyleDet[availableDRNonFitMeas[i]], markerSizeDet[availableDRNonFitMeas[i]], colorDet[availableDRNonFitMeas[i]],
                                 colorDet[availableDRNonFitMeas[i]]);
        statErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]]->Draw("p,same,z");
        legenDRNonFitelStatErr2->AddEntry(statErrorRelCollectionDRNonFit[availableDRNonFitMeas[i]],nameMeasGlobalLabel[availableDRNonFitMeas[i]],"p");
    }
    legenDRNonFitelStatErr2->Draw();
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

        TLegend* legenDRNonFitelTotErr       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
        legenDRNonFitelTotErr->AddEntry(graphCombDRNonFitRelTot,"tot","p");
        legenDRNonFitelTotErr->AddEntry(graphCombDRNonFitRelStat,"stat","l");
        legenDRNonFitelTotErr->AddEntry(graphCombDRNonFitRelSys,"sys","l");
        legenDRNonFitelTotErr->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DR_Reldecomp.%s",outputDir.Data(),suffix.Data()));


    //*******************************************************************************************************************************************
    //*********************************************** Combining IncGamma spectra ****************************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};

    Int_t offSetsIncGamma[11]       = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetsIncGammaSys[11]    = { 1,  0,  6,  0,  3,
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
                                                                                            offSetsGamma, offSetsGammaSys,
                                                                                            graphCombIncGammaStat, graphCombIncGammaSys,
                                                                                            fileNameIncGammaOutputWeighting, "2.76TeV", "GammaInc", kTRUE,
                                                                                            NULL, fileNameCorrelations );


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
    histo2DIncGammaWeights->GetXaxis()->SetNoExponent();
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

    histo2DRelStatErr->GetYaxis()->SetRangeUser(-0.2,17.5);
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
        histo2DRelErr->GetYaxis()->SetRangeUser(0,20.5);
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

//     graphCombIncGammaTot            = ApplyXshift(graphCombIncGammaTot, fitShiftingGamma);
//     cout << "comb" << endl;
//     graphCombIncGammaStat->Print();
//     graphCombIncGammaStat           = ApplyXshiftIndividualSpectra( graphCombIncGammaTot,
//                                                                     graphCombIncGammaStat,
//                                                                     fitShiftingGamma,
//                                                                     0, graphCombIncGammaStat->GetN());
//     graphCombIncGammaSys            = ApplyXshiftIndividualSpectra( graphCombIncGammaTot,
//                                                                     graphCombIncGammaSys,
//                                                                     fitShiftingGamma,
//                                                                     0, graphCombIncGammaSys->GetN());
    Int_t offSetGammaShifting[11]   = { 0,  0,  5,  0,  2,
                                        0,  0,  0,  0,  0,
                                        0 };
    Int_t nComBinsGammaShifting[11] = { 17, 0, 12, 0,  15,
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
    graphRatioGammaCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitTot, fitHagGammaComb);
    TGraphAsymmErrors* graphRatioGammaCombCombFitStat    = (TGraphAsymmErrors*)graphCombIncGammaStat->Clone();
    graphRatioGammaCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitStat, fitHagGammaComb);
    TGraphAsymmErrors* graphRatioGammaCombCombFitSys     = (TGraphAsymmErrors*)graphCombIncGammaSys->Clone();
    graphRatioGammaCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitSys, fitHagGammaComb);

    for (Int_t i= 0; i< 11; i++){
        if (graphIndGammaIncStat[i]){
            graphRatioGammaIndCombFitStat[i]             = (TGraphAsymmErrors*)graphIndGammaIncStat[i]->Clone(Form("RatioGamma%sToCombFitStat", nameMeasGlobalLabel[i].Data()));
            graphRatioGammaIndCombFitStat[i]             = CalculateGraphErrRatioToFit(graphRatioGammaIndCombFitStat[i], fitHagGammaComb);
        }
        if (graphIndGammaIncSys[i]){
            graphRatioGammaIndCombFitSys[i]              = (TGraphAsymmErrors*)graphIndGammaIncSys[i]->Clone(Form("RatioGamma%sToCombFitSyst", nameMeasGlobalLabel[i].Data()));
            graphRatioGammaIndCombFitSys[i]              = CalculateGraphErrRatioToFit(graphRatioGammaIndCombFitSys[i], fitHagGammaComb);
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
    histo2DGammaRatioToCombFit->GetXaxis()->SetNoExponent();
    //  histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(0.59,1.42);
    histo2DGammaRatioToCombFit->Draw("copy");

    ProduceGraphAsymmWithoutXErrors(graphRatioGammaCombCombFitStat);

    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitSys, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV, widthLinesBoxes, kTRUE);
    graphRatioGammaCombCombFitSys->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitStat, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
    graphRatioGammaCombCombFitStat->Draw("p,same,z");

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

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_PPb5023GeV.%s",outputDir.Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_PPb5023GeV.pdf",outputDir.Data()));

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
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit.pdf",outputDir.Data()));

    //*******************************************************************************************************
    //************************** Calculating combined direct photon spectrum ********************************
    //*******************************************************************************************************
    Double_t xArrayCombDR[graphCombDRStat->GetN()+1];
    xArrayCombDR[0] = graphCombDRStat->GetX()[0] - graphCombDRStat->GetEXhigh()[0];
    for (Int_t i = 1; i<graphCombDRStat->GetN()+1;i++){
        xArrayCombDR[i] = graphCombDRStat->GetX()[i-1] + graphCombDRStat->GetEXhigh()[i-1];
    }
    //_______________________ copy inclusive photon spectra _____________________
    TH1D *histoCombDirGammaSpectrumErrSum                   = new TH1D("histoCombDirGammaSpectrumErrSum","",graphCombDRStat->GetN(),xArrayCombDR);
    TH1D *histoCombDirGammaSpectrumErrSys                   = new TH1D("histoCombDirGammaSpectrumErrSys","",graphCombDRStat->GetN(),xArrayCombDR);
    TH1D *histoCombDirGammaSpectrumErrStat                  = new TH1D("histoCombDirGammaSpectrumErrStat","",graphCombDRStat->GetN(),xArrayCombDR);

    //_______________________ get arrays of double ratio errors __________________
    Double_t *SystErrorsCombDR                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *sumErrorsCombDR                               = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *StatErrorsCombDR                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *xErrorsDR                                     = new Double_t[graphCombIncGammaStat->GetN()];
    for (Int_t i = 0; i< graphCombDRStat->GetN(); i++){
        SystErrorsCombDR[i]                                 = graphCombDRSys->GetEYhigh()[i]/graphCombDRSys->GetY()[i] *100;
        StatErrorsCombDR[i]                                 = graphCombDRStat->GetEYhigh()[i]/graphCombDRStat->GetY()[i] *100;
        sumErrorsCombDR[i]                                  = graphCombDRTot->GetEYhigh()[i]/graphCombDRTot->GetY()[i] *100;
        //cout << i << "\t" << graphCombDRSys->GetY()[i] << "\t" << graphCombDRSys->GetEYhigh()[i] << "\t" <<SystErrorsCombDR[i] << endl;
    }
    xErrorsDR                                               = graphCombDRStat->GetX();

    cout << __LINE__ << endl;
    //graphCombDRTot->Print();

    //_______________________ copy inclusive photon spectra _____________________
    TH1D* histoCombErrorsForDRSum                           = new TH1D("histoCombErrorsForDRSum","",graphCombDRStat->GetN(),xArrayCombDR);
    TH1D* histoCombErrorsForDRStat                          = new TH1D("histoCombErrorsForDRStat","",graphCombDRStat->GetN(),xArrayCombDR);
    TH1D* histoCombErrorsForDRSys                           = new TH1D("histoCombErrorsForDRSys","",graphCombDRStat->GetN(),xArrayCombDR);

    for(Int_t i = 1; i<graphCombDRStat->GetN()+1;i++){
        //cout<< i << "\t"<<xErrorsDR[i-1]<<"  "<<histoCombErrorsForDRSum->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum->GetBinWidth(i) <<endl;
        Double_t binErrorSummed                                 = sumErrorsCombDR[i-1];
        Double_t binErrorSyst                                   = SystErrorsCombDR[i-1];
        Double_t binErrorStat                                   = StatErrorsCombDR[i-1];
        Double_t DR                                             = graphCombDRStat->GetY()[i-1];

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
    //************************** Calculating combined direct photon spectrum for NonFitted DR****************
    //*******************************************************************************************************
    Double_t xArrayCombDRNonFit[graphCombDRNonFitStat->GetN()+1];
    xArrayCombDRNonFit[0] = graphCombDRNonFitStat->GetX()[0] - graphCombDRNonFitStat->GetEXhigh()[0];
    for (Int_t i = 1; i<graphCombDRNonFitStat->GetN()+1;i++){
        xArrayCombDRNonFit[i] = graphCombDRNonFitStat->GetX()[i-1] + graphCombDRNonFitStat->GetEXhigh()[i-1];
    }
    //_______________________ copy inclusive photon spectra _____________________
    TH1D *histoCombDirGammaSpectrumNonFitErrSum                   = new TH1D("histoCombDirGammaSpectrumNonFitErrSum","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);
    TH1D *histoCombDirGammaSpectrumNonFitErrSys                   = new TH1D("histoCombDirGammaSpectrumNonFitErrSys","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);
    TH1D *histoCombDirGammaSpectrumNonFitErrStat                  = new TH1D("histoCombDirGammaSpectrumNonFitErrStat","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);

    //_______________________ get arrays of double ratio errors __________________
    Double_t *SystErrorsCombDRNonFit                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *sumErrorsCombDRNonFit                               = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *StatErrorsCombDRNonFit                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *xErrorsDRNonFit                                     = new Double_t[graphCombIncGammaStat->GetN()];
    for (Int_t i = 0; i< graphCombDRNonFitStat->GetN(); i++){
        SystErrorsCombDRNonFit[i]                                 = graphCombDRNonFitSys->GetEYhigh()[i]/graphCombDRNonFitSys->GetY()[i] *100;
        StatErrorsCombDRNonFit[i]                                 = graphCombDRNonFitStat->GetEYhigh()[i]/graphCombDRNonFitStat->GetY()[i] *100;
        sumErrorsCombDRNonFit[i]                                  = graphCombDRNonFitTot->GetEYhigh()[i]/graphCombDRNonFitTot->GetY()[i] *100;
        //cout << i << "\t" << graphCombDRNonFitSys->GetY()[i] << "\t" << graphCombDRNonFitSys->GetEYhigh()[i] << "\t" <<SystErrorsCombDRNonFit[i] << endl;
    }
    xErrorsDRNonFit                                               = graphCombDRNonFitStat->GetX();

    cout << __LINE__ << endl;
    //graphCombDRNonFitTot->Print();

    //_______________________ copy inclusive photon spectra _____________________
    TH1D* histoCombErrorsForDRNonFitSum                           = new TH1D("histoCombErrorsForDRNonFitSum","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);
    TH1D* histoCombErrorsForDRNonFitStat                          = new TH1D("histoCombErrorsForDRNonFitStat","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);
    TH1D* histoCombErrorsForDRNonFitSys                           = new TH1D("histoCombErrorsForDRNonFitSys","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);

    for(Int_t i = 1; i<graphCombDRNonFitStat->GetN()+1;i++){
        //cout<< i << "\t"<<xErrorsDRNonFit[i-1]<<"  "<<histoCombErrorsForDRNonFitSum->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRNonFitSum->GetBinWidth(i) <<endl;
        Double_t binErrorSummed                                 = sumErrorsCombDRNonFit[i-1];
        Double_t binErrorSyst                                   = SystErrorsCombDRNonFit[i-1];
        Double_t binErrorStat                                   = StatErrorsCombDRNonFit[i-1];
        Double_t DRNonFit                                             = graphCombDRNonFitStat->GetY()[i-1];

        //cout << DRNonFit << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
        histoCombErrorsForDRNonFitSum->SetBinContent(i,DRNonFit);
        histoCombErrorsForDRNonFitSys->SetBinContent(i,DRNonFit);
        histoCombErrorsForDRNonFitStat->SetBinContent(i,DRNonFit);
        histoCombErrorsForDRNonFitSum->SetBinError(i,(binErrorSummed/100)*DRNonFit);
        histoCombErrorsForDRNonFitSys->SetBinError(i,(binErrorSyst/100)*DRNonFit);
        histoCombErrorsForDRNonFitStat->SetBinError(i,(binErrorStat/100)*DRNonFit);
    }

    for(Int_t i = 1; i<histoCombErrorsForDRNonFitSum->GetNbinsX()+1;i++){
        histoCombDirGammaSpectrumNonFitErrSum->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumNonFitErrSys->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumNonFitErrStat->SetBinContent(i+1,-1);

        histoCombDirGammaSpectrumNonFitErrSum->SetBinError(i+1,0);
        histoCombDirGammaSpectrumNonFitErrSys->SetBinError(i+1,0);
        histoCombDirGammaSpectrumNonFitErrStat->SetBinError(i+1,0);
    }

    // get the binning of the direct photons from the DRNonFit
    TH1D *histoCombDirGammaSpecNonFitSysErr                        = new TH1D(*histoCombErrorsForDRNonFitSys);
    TH1D *histoCombDirGammaSpecNonFitStatErr                       = new TH1D(*histoCombErrorsForDRNonFitStat);
    TH1D *histoCombDirGammaSpecNonFitSumErr                        = new TH1D(*histoCombErrorsForDRNonFitSum);


    cout << __LINE__ << endl;
    //graphCombDRNonFitStat->Print();
    for(Int_t i = 1; i<graphCombDRNonFitStat->GetN()+1; i++){
        // obtain common quantities
        Double_t Rgamma                 = histoCombErrorsForDRNonFitSys->GetBinContent(i);
        Double_t nIncGamma              = graphCombIncGammaStat->GetY()[i-1];

        // calculating Systematics graph
        Double_t errRgamma              = histoCombErrorsForDRNonFitSys->GetBinError(i);
        Double_t errNIncGam             = graphCombIncGammaSys->GetEYhigh()[i-1];
        Double_t q1                     = 1 - 1/ Rgamma;

        Double_t q1Error                = errRgamma/(Rgamma*Rgamma);
        Double_t content                = nIncGamma * ( 1 - 1/ Rgamma);
        Double_t error                  = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        Double_t errDRNonFit                  = content - error;
        histoCombDirGammaSpecNonFitSysErr->SetBinError(i, error);
        histoCombDirGammaSpecNonFitSysErr->SetBinContent(i, content);
        histoCombDirGammaSpectrumNonFitErrSys->SetBinContent(i, errDRNonFit);

        // calculating Stat graphs
        errRgamma                       = histoCombErrorsForDRNonFitStat->GetBinError(i);
        errNIncGam                      = graphCombIncGammaStat->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDRNonFit                           = content - error;
        histoCombDirGammaSpecNonFitStatErr->SetBinError(i, error);
        histoCombDirGammaSpecNonFitStatErr->SetBinContent(i, content);
        histoCombDirGammaSpectrumNonFitErrStat->SetBinContent(i, errDRNonFit);

        // calculating summed error graphs
        errRgamma                       = histoCombErrorsForDRNonFitSum->GetBinError(i);
        errNIncGam                      = graphCombIncGammaTot->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDRNonFit                           = content - error;
        histoCombDirGammaSpecNonFitSumErr->SetBinError(i, error);
        histoCombDirGammaSpecNonFitSumErr->SetBinContent(i, content);
        histoCombDirGammaSpectrumNonFitErrSum->SetBinContent(i, errDRNonFit);
    }

    // purely calculating points based on all Systematic errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumNonFitSystErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumNonFitErrSys,histoCombDirGammaSpecNonFitStatErr,0,0.5);
    if(graphCombDirGammaSpectrumNonFitSystErr)graphCombDirGammaSpectrumNonFitSystErr->SetName("graphCombDirGammaSpectrumNonFitSystErr");
    if(graphCombDirGammaSpectrumNonFitSystErr) cout << "sys has been found" << endl;
    if(graphCombDirGammaSpectrumNonFitSystErr)graphCombDirGammaSpectrumNonFitSystErr->Print();

    // purely calculating points based on Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumNonFitStatErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumNonFitErrStat,histoCombDirGammaSpecNonFitStatErr,0,0.5);
    if(graphCombDirGammaSpectrumNonFitStatErr)graphCombDirGammaSpectrumNonFitStatErr->SetName("graphCombDirGammaSpectrumNonFitStatErr");
    if(graphCombDirGammaSpectrumNonFitStatErr) cout << "stat has been found" << endl;
    if(graphCombDirGammaSpectrumNonFitStatErr)graphCombDirGammaSpectrumNonFitStatErr->Print();
    // purely calculating points based on all Systematic + Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumNonFitSumErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumNonFitErrSum,histoCombDirGammaSpecNonFitStatErr,0,0.5);
    if(graphCombDirGammaSpectrumNonFitSumErr)graphCombDirGammaSpectrumNonFitSumErr->SetName("graphCombDirGammaSpectrumNonFitSumErr");
    if(graphCombDirGammaSpectrumNonFitSumErr) cout << "tot has been found" << endl;
    if(graphCombDirGammaSpectrumNonFitSumErr)graphCombDirGammaSpectrumNonFitSumErr->Print();
    // calculate points above confidence level summed errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumNonFitSumErrConfi = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumNonFitErrSum,histoCombDirGammaSpecNonFitStatErr,2,0.5);
    if(graphCombDirGammaSpectrumNonFitSumErrConfi)graphCombDirGammaSpectrumNonFitSumErrConfi->SetName("graphCombDirGammaSpectrumNonFitSumErrConfi");
    if(graphCombDirGammaSpectrumNonFitSumErrConfi) cout << "confi has been found" << endl;
    if(graphCombDirGammaSpectrumNonFitSumErrConfi)graphCombDirGammaSpectrumNonFitSumErrConfi->Print();
    // calculate arrows for points with 0, error summed
    TGraphAsymmErrors *graphCombDirGammaSpectrumNonFitSumErrAr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumNonFitErrSum,histoCombDirGammaSpecNonFitStatErr,5,0.5);
    if(graphCombDirGammaSpectrumNonFitSumErrAr)graphCombDirGammaSpectrumNonFitSumErrAr->SetName("graphCombDirGammaSpectrumNonFitSumErrAr");
    if(graphCombDirGammaSpectrumNonFitSumErrAr) cout << "Ar has been found" << endl;
    if(graphCombDirGammaSpectrumNonFitSumErrAr)graphCombDirGammaSpectrumNonFitSumErrAr->Print();



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
    hist2DDRDummySingle->GetXaxis()->SetNoExponent();
    hist2DDRDummySingle->GetXaxis()->SetMoreLogLabels(kTRUE);
    hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRSingle = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*2,0.5,0.953, textSizeSinglePad, 2, "", 42, 0.3);
        legendDRSingle->SetTextAlign(11);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 10; i > -1; i--){
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


        TLatex *labelDRSingle = new TLatex(0.95,0.92,collisionSystempp2760GeV.Data());
        SetStyleTLatex( labelDRSingle, textSizeSinglePad,4, 1, 42, kTRUE, 31);
        labelDRSingle->Draw();
        TLatex *labelALICEDRSingle = new TLatex(0.95,0.87,textALICE.Data());
        SetStyleTLatex( labelALICEDRSingle, textSizeSinglePad,4, 1, 42, kTRUE, 31);
        labelALICEDRSingle->Draw();


        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pp2760GeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pp2760GeV.pdf", outputDir.Data()));
    
    hist2DDRDummySingle->DrawCopy();

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 10; i > -1; i--){
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

    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pp2760GeV_NonFit.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pp2760GeV_NonFit.pdf", outputDir.Data()));

    hist2DDRDummySingle->DrawCopy();
        TGraphAsymmErrors* graphCombDRStatPlot    = (TGraphAsymmErrors*)graphCombDRStat->Clone("graphCombDRStatPlot");
        ProduceGraphAsymmWithoutXErrors(graphCombDRStatPlot);


        TLegend* legendDRComb = GetAndSetLegend2(0.12,0.95-textSizeSinglePad*1,0.5,0.95, textSizeSinglePad, 1, "", 42, 0.15);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        DrawGammaSetMarkerTGraphAsym(graphCombDRSys, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRStatPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV, widthLinesBoxes);
        legendDRComb->AddEntry(graphCombDRSys,"ALICE","pf");
        graphCombDRSys->Draw("E2same");
        graphCombDRStatPlot->Draw("Epsame");
//         legendDRComb->Draw();

        labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_Comb_pp2760GeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_Comb_pp2760GeV.pdf", outputDir.Data()));
    
    hist2DDRDummySingle->DrawCopy();
        TGraphAsymmErrors* graphCombDRNonFitStatPlot    = (TGraphAsymmErrors*)graphCombDRNonFitStat->Clone("graphCombDRNonFitStatPlot");
        ProduceGraphAsymmWithoutXErrors(graphCombDRNonFitStatPlot);


        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitSys, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitStatPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV, widthLinesBoxes);
        graphCombDRNonFitSys->Draw("E2same");
        graphCombDRNonFitStatPlot->Draw("Epsame");

        labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_Comb_pp2760GeV_NonFit.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_Comb_pp2760GeV_NonFit.pdf", outputDir.Data()));

        hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRTheoryComb = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*1,0.5,0.96, textSizeSinglePad, 1, "", 42, 0.15);
        legendDRTheoryComb->SetTextAlign(11);
        TLegend* legendDRTheoryComb2 = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*5,0.5,0.96-textSizeSinglePad*1, textSizeSinglePad, 1, "NLO pQCD:", 42, 0.15);
        legendDRTheoryComb2->SetTextAlign(11);

        DrawGammaSetMarkerTGraphAsym(graphCombDRSys, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRStatPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV, widthLinesBoxes);
        legendDRTheoryComb->AddEntry(graphCombDRSys,"ALICE","pf");

        if (graphTheoryNLODRpp2760GeV) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryNLODRpp2760GeV, 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
            graphTheoryNLODRpp2760GeV->Draw("3,same");
            legendDRTheoryComb2->AddEntry(dummyVogelsangforLegend,"PDF: CT10, FF: GRV","fl");
        }
        if (graphTheoryNLODRpp2760GeVCenter){
            DrawGammaNLOTGraph( graphTheoryNLODRpp2760GeVCenter, 2, styleLineNLOWerner, colorNLOWerner);
            graphTheoryNLODRpp2760GeVCenter->Draw("lc,same");
        }
        if (graphTheoryNLODRpp2760GeVPaquettCenter){
            DrawGammaNLOTGraph( graphTheoryNLODRpp2760GeVPaquettCenter, 2, styleLineMcGill, colorNLOMcGill );
            graphTheoryNLODRpp2760GeVPaquettCenter->RemovePoint(0);
            graphTheoryNLODRpp2760GeVPaquettCenter->Draw("lc,same");
            legendDRTheoryComb2->AddEntry(graphTheoryNLODRpp2760GeVPaquettCenter,"PDF: CTEQ6.1M, FF: BFG2","l");
        }
        if (graphTheoryJETPHOXDRpp2760GeV) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryJETPHOXDRpp2760GeV, 0, 0, colorJETPHOXBand, colorJETPHOXBand, 0.2, kTRUE, colorJETPHOXBand);
            graphTheoryJETPHOXDRpp2760GeV->Draw("3,same");
            legendDRTheoryComb2->AddEntry(dummyJETPHOXforLegend,"JETPHOX","l");
        }
        if (graphTheoryJETPHOXDRpp2760GeVCenter){
            DrawGammaNLOTGraph( graphTheoryJETPHOXDRpp2760GeVCenter, 2, styleLineJETPHOX, colorJETPHOX);
            graphTheoryJETPHOXDRpp2760GeVCenter->Draw("lc,same");
        }
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRSys->Draw("E2same");
        graphCombDRStatPlot->Draw("Epsame");
        legendDRTheoryComb->Draw();
        legendDRTheoryComb2->Draw();

        labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pp2760GeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pp2760GeV.pdf", outputDir.Data()));
    
        hist2DDRDummySingle->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitSys, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitStatPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV, widthLinesBoxes);

        if (graphTheoryNLODRpp2760GeV) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryNLODRpp2760GeV, 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
            graphTheoryNLODRpp2760GeV->Draw("3,same");
        }
        if (graphTheoryNLODRpp2760GeVCenter){
            DrawGammaNLOTGraph( graphTheoryNLODRpp2760GeVCenter, 2, styleLineNLOWerner, colorNLOWerner);
            graphTheoryNLODRpp2760GeVCenter->Draw("lc,same");
        }
        if (graphTheoryNLODRpp2760GeVPaquettCenter){
            DrawGammaNLOTGraph( graphTheoryNLODRpp2760GeVPaquettCenter, 2, styleLineMcGill, colorNLOMcGill );
            graphTheoryNLODRpp2760GeVPaquettCenter->RemovePoint(0);
            graphTheoryNLODRpp2760GeVPaquettCenter->Draw("lc,same");
        }
        if (graphTheoryJETPHOXDRpp2760GeV) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryJETPHOXDRpp2760GeV, 0, 0, colorJETPHOXBand, colorJETPHOXBand, 0.2, kTRUE, colorJETPHOXBand);
            graphTheoryJETPHOXDRpp2760GeV->Draw("3,same");
        }
        if (graphTheoryJETPHOXDRpp2760GeVCenter){
            DrawGammaNLOTGraph( graphTheoryJETPHOXDRpp2760GeVCenter, 2, styleLineJETPHOX, colorJETPHOX);
            graphTheoryJETPHOXDRpp2760GeVCenter->Draw("lc,same");
        }
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRNonFitSys->Draw("E2same");
        graphCombDRNonFitStatPlot->Draw("Epsame");
        legendDRTheoryComb->Draw();
        legendDRTheoryComb2->Draw();

        labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pp2760GeV_NonFit.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pp2760GeV_NonFit.pdf", outputDir.Data()));



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
    histo1DEff->GetYaxis()->SetRangeUser(0.22, 0.82 );
    histo1DEff->GetYaxis()->SetLabelOffset(0.001);
    histo1DEff->GetXaxis()->SetNoExponent();
    histo1DEff->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo1DEff->DrawCopy();

        TLegend* legendEffiGamma           = GetAndSetLegend2(0.57, 0.13, 0.95, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel, 2);
        for (Int_t i = 0; i < 11; i++){
            if (histoEffi[i]){
                DrawGammaSetMarker(histoEffi[i],  markerStyleDetMC[i], markerSizeDetMC[i], colorDetMC[i] , colorDetMC[i]);
                histoEffi[i]->Draw("p,same,e");
                legendEffiGamma->AddEntry(histoEffi[i],"   ","p");
            } else if (histoEffiMCPt[i]){
                legendEffiGamma->AddEntry((TObject*)0,"   ","");
            }
            if (histoEffiMCPt[i]){
                DrawGammaSetMarker(histoEffiMCPt[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
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
        TLatex *labelEnergyEffi         = new TLatex(0.13,0.87,collisionSystempp2760GeV.Data());
        SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
        labelEnergyEffi->Draw();
        TLatex *labelPerfEffiPTrec      = new TLatex(0.57,0.145+(3*textSizeLabelsRel),"#it{p}_{T}^{rec}");
        SetStyleTLatex( labelPerfEffiPTrec, textSizeLabelsRel,4);
        labelPerfEffiPTrec->Draw();
        TLatex *labelPerfEffiPTtrue      = new TLatex(0.665,0.145+(3*textSizeLabelsRel),"#it{p}_{T}^{true}");
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
    histo2DResCor->GetXaxis()->SetNoExponent();
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
        TLatex *labelEnergyResolCor         = new TLatex(0.15,0.87,collisionSystempp2760GeV.Data());
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
    histo1DPurity->GetXaxis()->SetNoExponent();
    histo1DPurity->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo1DPurity->DrawCopy();

        TLegend* legendPurityGamma           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
        DrawGammaLines(0.23, 31., 1., 1., 1.2, kGray+2, 7);
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
        TLatex *labelEnergyPurity         = new TLatex(0.15,0.87,collisionSystempp2760GeV.Data());
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


    TCanvas* canvasConvProb                 = new TCanvas("canvasConvProb", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasConvProb,  0.11, 0.01, 0.015, 0.095);
    //     canvasConvProb->SetLogy(1);
    canvasConvProb->SetLogx(1);

    TH1F * histo1DConvProb                  = new TH1F("histo1DConvProb", "histo1DConvProb",1000, doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(  histo1DConvProb, "#it{p}_{T} (GeV/#it{c})","#it{P}_{conv}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.2);//(#times #epsilon_{pur})
    histo1DConvProb->GetYaxis()->SetRangeUser(0.04, 0.12 );
    histo1DConvProb->GetYaxis()->SetLabelOffset(0.003);
    histo1DConvProb->GetXaxis()->SetNoExponent();
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
        TLatex *labelEnergyConvProb         = new TLatex(0.15,0.87,collisionSystempp2760GeV.Data());
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
    histo1DPileUp->GetYaxis()->SetRangeUser(0.83, 1.03 );
    histo1DPileUp->GetYaxis()->SetLabelOffset(0.005);
    histo1DPileUp->GetXaxis()->SetNoExponent();
    histo1DPileUp->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo1DPileUp->DrawCopy();

        TLegend* legendPileUpGamma           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(2*textSizeLabelsRel),textSizeLabelsPixel);
        DrawGammaLines(0.23, 31., 1., 1., 1.2, kGray+2, 7);
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
        TLatex *labelEnergyPileUp         = new TLatex(0.15,0.87,collisionSystempp2760GeV.Data());
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
    histo1DEffectiveSecCorr->GetXaxis()->SetNoExponent();
    histo1DEffectiveSecCorr->GetXaxis()->SetMoreLogLabels(kTRUE);


    Double_t minYSecCorr[4]             = {0.0, 0.0, 0.0, 0.0};
    Double_t maxYSecCorr[4]             = {0.06, 0.004, 0.8e-3, 0.041};
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
        TLatex *labelEnergyEffectiveSecCorr         = new TLatex(0.15,0.84,collisionSystempp2760GeV.Data());
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
    SetStyleHistoTH1ForGraphs(histo2DYieldGamma, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.7);
    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-9,3.5e0);
    histo2DYieldGamma->GetXaxis()->SetMoreLogLabels();
    histo2DYieldGamma->GetXaxis()->SetNoExponent();
    histo2DYieldGamma->Draw("copy");

        TGraphAsymmErrors* graphCombIncGammaStatPlot    = (TGraphAsymmErrors*)graphCombIncGammaStat->Clone("graphCombIncGammaStatPlot");
        ProduceGraphAsymmWithoutXErrors(graphCombIncGammaStatPlot);

        TLegend* legendYieldIncGamma           = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(3*textSizeLabelsRel*0.85),textSizeLabelsPixel);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV,widthLinesBoxes);
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");
        legendYieldIncGamma->AddEntry(graphCombIncGammaSys,"ALICE","pf");
        legendYieldIncGamma->Draw();

        DrawGammaSetMarkerTF1( fitHagGammaComb, 7, 2, colorCombpp2760GeV);
        legendYieldIncGamma->AddEntry(fitHagGammaComb,"mod. Hagedorn fit","l");
        fitHagGammaComb->SetRange(0.4, 10);
        fitHagGammaComb->Draw("same");
        DrawGammaSetMarkerTF1( fitTsallisGammaComb, 5, 2, colorCombpp2760GeVBox);
        legendYieldIncGamma->AddEntry(fitTsallisGammaComb,"Tsalis fit","l");
        fitTsallisGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTsallisGammaComb->Draw("same");

        TLatex *labelEnergyInvYieldPaperAll = new TLatex(0.20, 0.20+0.04*3, collisionSystempp2760GeV.Data());
        SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelEnergyInvYieldPaperAll->Draw();
        TLatex *labelALICEInvYieldPaperAll  = new TLatex(0.20,0.20+0.04*2,textALICE.Data());
        SetStyleTLatex( labelALICEInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelALICEInvYieldPaperAll->Draw();
        TLatex *labelALICENormUnPaperAll    = new TLatex(0.20,0.20+0.05*1,"Norm. unc. 2.5%");
        SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICENormUnPaperAll->Draw();

        histo2DYieldGamma->Draw("same,axis");

    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma.pdf",outputDir.Data()));


    // **********************************************************************************************************************
    // ******************************** InvYield for individual inc gamma measurements **************************************
    // **********************************************************************************************************************
    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-9,3.5e0);
    histo2DYieldGamma->Draw("copy");
    TLegend* legendYieldIncGammaInd       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(3*textSizeLabelsRel*0.85),textSizeLabelsPixel);
    for (Int_t i = 0; i < 11; i++){
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

    // **********************************************************************************************************************
    // ******************************** InvYield for combined dir gamma measurement **************************************
    // **********************************************************************************************************************

    // prep plot graphs
    TGraphAsymmErrors* graphCombDirGammaSpectrumStatErrPlot = NULL;
    if (graphCombDirGammaSpectrumStatErr) graphCombDirGammaSpectrumStatErrPlot       = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr->Clone("graphCombDirGammaSpectrumStatErrPlot");
    if (graphCombDirGammaSpectrumStatErrPlot) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErrPlot);

    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-11,3.5e0);
    histo2DYieldGamma->Draw("copy");


        if (graphCombDirGammaSpectrumSystErr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErrPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
            graphCombDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrAr , 1, 3, colorCombpp2760GeV, colorCombpp2760GeV, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr);
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

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma.pdf",outputDir.Data()));
    
    // prep plot graphs
    TGraphAsymmErrors* graphCombDirGammaSpectrumNonFitStatErrPlot = NULL;
    if (graphCombDirGammaSpectrumNonFitStatErr) graphCombDirGammaSpectrumNonFitStatErrPlot       = (TGraphAsymmErrors*)graphCombDirGammaSpectrumNonFitStatErr->Clone("graphCombDirGammaSpectrumNonFitStatErrPlot");
    if (graphCombDirGammaSpectrumNonFitStatErrPlot) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumNonFitStatErrPlot);

    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-11,3.5e0);
    histo2DYieldGamma->Draw("copy");


        if (graphCombDirGammaSpectrumNonFitSystErr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumNonFitSystErr, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumNonFitSystErr->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumNonFitStatErrPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumNonFitStatErrPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
            graphCombDirGammaSpectrumNonFitStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumNonFitSumErrAr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumNonFitSumErrAr , 1, 3, colorCombpp2760GeV, colorCombpp2760GeV, 1.8, kTRUE);
            graphCombDirGammaSpectrumNonFitSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumNonFitSumErrAr);
        }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_NonFit.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_NonFit.pdf",outputDir.Data()));
    histo2DYieldGamma->Draw("copy");

        TLegend* legendYieldDirGamma        = GetAndSetLegend2(0.70, 0.84-(2*textSizeLabelsRel*0.85), 0.93, 0.84,textSizeLabelsPixel, 1, "", 43, 0.3);

        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys, markerStyleCombpp2760GeV+4, markerSizeCombpp2760GeV+0.2, colorCombpp2760GeVBox , colorCombpp2760GeVBox,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot, markerStyleCombpp2760GeV+4, markerSizeCombpp2760GeV+0.2, colorCombpp2760GeVBox , colorCombpp2760GeVBox, widthLinesBoxes);
        legendYieldDirGamma->AddEntry(graphCombIncGammaSys, "#gamma_{inc} ALICE","pf");
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");

        if (graphCombDirGammaSpectrumSystErr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErrPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
            graphCombDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrAr , 1, 3, colorCombpp2760GeV, colorCombpp2760GeV, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr);
            legendYieldDirGamma->AddEntry((TObject*)0, "#gamma_{dir} ALICE","");
        }

        TGraphAsymmErrors* dummyForLegend    = new TGraphAsymmErrors(1);
        dummyForLegend->SetPoint(0,4.2,0.02);
        dummyForLegend->SetPointError(0,0,0,0.009,0);
        DrawGammaSetMarkerTGraphAsym(dummyForLegend , 1, 3, colorCombpp2760GeV, colorCombpp2760GeV, 1.8, kTRUE);
        if (graphCombDirGammaSpectrumSumErrAr){
            dummyForLegend->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend);
        }
        legendYieldDirGamma->Draw();

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma.pdf",outputDir.Data()));
    
    histo2DYieldGamma->Draw("copy");


        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");

        if (graphCombDirGammaSpectrumNonFitSystErr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumNonFitSystErr, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumNonFitSystErr->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumNonFitStatErrPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumNonFitStatErrPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
            graphCombDirGammaSpectrumNonFitStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumNonFitSumErrAr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumNonFitSumErrAr , 1, 3, colorCombpp2760GeV, colorCombpp2760GeV, 1.8, kTRUE);
            graphCombDirGammaSpectrumNonFitSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumNonFitSumErrAr);
        }

        
        if (graphCombDirGammaSpectrumNonFitSumErrAr){
            dummyForLegend->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend);
        }
        legendYieldDirGamma->Draw();

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_NonFit.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_NonFit.pdf",outputDir.Data()));
    histo2DYieldGamma->Draw("copy");

    TLegend* legendYieldDirGammaTheo2      = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(4*textSizeLabelsRel*0.85),textSizeLabelsPixel,1, "#gamma_{dir} NLO pQCD:", 43, 0.23);
    legendYieldDirGammaTheo2->SetTextAlign(12);
    graphCombIncGammaSys->Draw("E2same");
    graphCombIncGammaStatPlot->Draw("Epsame");

        if (graphCombDirGammaSpectrumSystErr){
            graphCombDirGammaSpectrumSystErr->Draw("E2same");
        }
        legendYieldDirGamma->Draw();
        if (graphTheoryNLOpp2760GeV) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryNLOpp2760GeV, 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
            graphTheoryNLOpp2760GeV->Draw("3,same");
            graphTheoryNLOpp2760GeV->SetLineWidth(widthLineNLO*1.5);
            legendYieldDirGammaTheo2->AddEntry(graphTheoryNLOpp2760GeV,"PDF: CT10, FF: GRV","l");
        }

        if (graphTheoryNLOpp2760GeVPaquettCenter) {
            DrawGammaNLOTGraph( graphTheoryNLOpp2760GeVPaquettCenter, 2, styleLineMcGill, colorNLOMcGill );
            graphTheoryNLOpp2760GeVPaquettCenter->RemovePoint(0);
            graphTheoryNLOpp2760GeVPaquettCenter->Draw("lc,same");
            legendYieldDirGammaTheo2->AddEntry(graphTheoryNLOpp2760GeVPaquettCenter,"PDF: CTEQ6.1M, FF: BFG2","l");
        }
        
        if (graphTheoryJETPHOXpp2760GeV) {
          DrawGammaNLOTGraph( graphTheoryJETPHOXpp2760GeV, 2, styleLineJETPHOX, colorJETPHOX );
            graphTheoryJETPHOXpp2760GeV->Draw("l,same");
            graphTheoryJETPHOXpp2760GeV->SetLineWidth(widthLineNLO*1.5);
            legendYieldDirGammaTheo2->AddEntry(graphTheoryJETPHOXpp2760GeV,"JETPHOX","l");
        }
        legendYieldDirGammaTheo2->Draw();

        if (graphCombDirGammaSpectrumStatErrPlot){
            graphCombDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr){
            graphCombDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr);
        }
        if (graphCombDirGammaSpectrumSumErrAr){
            dummyForLegend->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend);
        }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory.pdf",outputDir.Data()));
    
    histo2DYieldGamma->Draw("copy");

    graphCombIncGammaSys->Draw("E2same");
    graphCombIncGammaStatPlot->Draw("Epsame");

        if (graphCombDirGammaSpectrumNonFitSystErr){
            graphCombDirGammaSpectrumNonFitSystErr->Draw("E2same");
        }
        legendYieldDirGamma->Draw();
        if (graphTheoryNLOpp2760GeV) {
            graphTheoryNLOpp2760GeV->Draw("3,same");
            graphTheoryNLOpp2760GeV->SetLineWidth(widthLineNLO*1.5);
        }

        if (graphTheoryNLOpp2760GeVPaquettCenter) {
            graphTheoryNLOpp2760GeVPaquettCenter->Draw("lc,same");
        }
        if (graphTheoryJETPHOXpp2760GeV) {
            graphTheoryJETPHOXpp2760GeV->Draw("l,same");
        }
        legendYieldDirGammaTheo2->Draw();

        if (graphCombDirGammaSpectrumNonFitStatErrPlot){
            graphCombDirGammaSpectrumNonFitStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumNonFitSumErrAr){
            graphCombDirGammaSpectrumNonFitSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumNonFitSumErrAr);
        }
        if (graphCombDirGammaSpectrumNonFitSumErrAr){
            dummyForLegend->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend);
        }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory_NonFit.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory_NonFit.pdf",outputDir.Data()));


    //***************************** Plot cocktail gammas to all gammas ratio ****************************************
    Style_t     cocktailColorPartialSums[14]        = {kRed+2,kBlue+1,kYellow+2,kOrange+1,kAzure-2,kRed-2,kViolet,kGreen-3,kOrange+6, kTeal+9,kMagenta-3,kCyan+4,kViolet+4,kAzure-4};

    TDirectory* directoryCocktailpp2760GeV[11]     = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoGammaFractionsCocktail[14][11];
    TGraphAsymmErrors* graphGammaFractionsCocktail[14][11];
    TH1D* histoGammaFractionsCocktailPartialSums[14][11];
    TGraphAsymmErrors* graphGammaFractionsCocktailPartialSums[14][11];
    TGraphAsymmErrors* graphGammaFractionsCocktailBand[14];
    TGraphAsymmErrors* graphGammaFractionsCocktailPartialSumsBand[14];
    for(Int_t i=0;i<14;i++){
        for(Int_t j=0;j<11;j++){
            histoGammaFractionsCocktail[i][j]= NULL;
            graphGammaFractionsCocktail[i][j]= NULL;
            histoGammaFractionsCocktailPartialSums[i][j]= NULL;
            graphGammaFractionsCocktailPartialSums[i][j]= NULL;
        }
        graphGammaFractionsCocktailBand[i]= NULL;
        graphGammaFractionsCocktailPartialSumsBand[i]= NULL;
    }
    TFile* fileCocktailGamma2760GeV                    = new TFile( inputFileNameCocktail.Data());
    directoryCocktailpp2760GeV[0]       = (TDirectory*)fileCocktailGamma2760GeV->Get("2.76TeV_PCM");
    directoryCocktailpp2760GeV[2]       = (TDirectory*)fileCocktailGamma2760GeV->Get("2.76TeV_EMCal");
    directoryCocktailpp2760GeV[4]       = (TDirectory*)fileCocktailGamma2760GeV->Get("2.76TeV_PCMEMCal");
    directoryCocktailpp2760GeV[5]       = (TDirectory*)fileCocktailGamma2760GeV->Get("2.76TeV_Comb");
    for(Int_t k=0;k<11;k++){
        if(directoryCocktailpp2760GeV[k]){
            for(Int_t i=0;i<14;i++){
                histoGammaFractionsCocktail[i][k]                 = (TH1D*) directoryCocktailpp2760GeV[k]->Get(Form("Gamma_From_%s_Pt_OrBin_RatioToAll",fParticle[i].Data()));
                graphGammaFractionsCocktail[i][k] = new TGraphAsymmErrors(histoGammaFractionsCocktail[i][k]);
                if(k==5)
                    graphGammaFractionsCocktailBand[i] = new TGraphAsymmErrors(histoGammaFractionsCocktail[i][k]);
                histoGammaFractionsCocktailPartialSums[i][k] = histoGammaFractionsCocktail[i][k];
            }
        }
    }


    Double_t newBinErrorLow  = 0.;
    Double_t newBinErrorHigh  = 0.;
    Bool_t lowChanged=kFALSE;
    Bool_t highChanged=kFALSE;

    for(Int_t i=0;i<14;i++){
        for(Int_t k = 0; k <= graphGammaFractionsCocktail[i][5]->GetN()-1; k++){
            newBinErrorLow        = graphGammaFractionsCocktail[i][5]->GetY()[k]-graphGammaFractionsCocktail[i][5]->GetEYlow()[k];
            newBinErrorHigh       = graphGammaFractionsCocktail[i][5]->GetY()[k]+graphGammaFractionsCocktail[i][5]->GetEYhigh()[k];
            lowChanged=kFALSE;
            highChanged=kFALSE;
            for(Int_t j=0;j<11;j++){
                if(graphGammaFractionsCocktail[i][j]){
                    if((newBinErrorLow-graphGammaFractionsCocktail[i][j]->GetY()[k]-graphGammaFractionsCocktail[i][j]->GetEYlow()[k])>0){
                        newBinErrorLow = graphGammaFractionsCocktail[i][j]->GetY()[k]-graphGammaFractionsCocktail[i][j]->GetEYlow()[k];
                        lowChanged = kTRUE;
                    }
                    if((newBinErrorHigh-(graphGammaFractionsCocktail[i][j]->GetY()[k]+graphGammaFractionsCocktail[i][j]->GetEYlow()[k]))<0){
                        newBinErrorHigh = graphGammaFractionsCocktail[i][j]->GetY()[k]+graphGammaFractionsCocktail[i][j]->GetEYhigh()[k];
                        highChanged=kTRUE;
                    }
                }
            }
            if(lowChanged)
                graphGammaFractionsCocktailBand[i]->SetPointEYlow(k,     graphGammaFractionsCocktail[i][5]->GetY()[k]-newBinErrorLow);
            if(highChanged)
                graphGammaFractionsCocktailBand[i]->SetPointEYhigh(k,     newBinErrorHigh-graphGammaFractionsCocktail[i][5]->GetY()[k]);
        }
        while(graphGammaFractionsCocktailBand[i]->GetX()[0] < 0.2) graphGammaFractionsCocktailBand[i]->RemovePoint(0);
    }

    TCanvas *canvasGammasRatio2                                 = new TCanvas("canvasGammasRatio2","",1100,1200);
    DrawGammaCanvasSettings(canvasGammasRatio2, 0.12, 0.01, 0.01, 0.075);
    canvasGammasRatio2->SetLogy();
    canvasGammasRatio2->SetLogx();
    TLegend* legendGammasRatio2                                 = GetAndSetLegend2(0.2, 0.96-(40*1.15*(2+1)/1200), 0.95, 0.96, 40, 7);
    legendGammasRatio2->SetHeader("#gamma from");
    TH1D* dummyHist                                                   = new TH1D("dummyHist", "", 1000,  0.25, 19.9);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#gamma_{source} / #gamma_{decay}", 1e-6, 5e1, 1.0, 1.3);
    dummyHist->SetLabelOffset(-0.015, "X");
    dummyHist->SetTitleOffset(0.8, "X");
    dummyHist->GetYaxis()->SetTitleFont(62);
    dummyHist->GetXaxis()->SetTitleFont(62);
    dummyHist->Draw();

    for (Int_t i=0; i<14; i++) {
        if (histoGammaFractionsCocktail[i][0]) {
            DrawGammaSetMarker(             histoGammaFractionsCocktail[i][5], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendGammasRatio2->AddEntry(   histoGammaFractionsCocktail[i][5], fParticleLatex[i].Data(), "l");
            histoGammaFractionsCocktail[i][5]->SetLineWidth(3);
            if(i<4)
                graphGammaFractionsCocktailBand[i]->SetFillColor(cocktailColor[i]);
            else
                graphGammaFractionsCocktailBand[i]->SetFillColorAlpha(cocktailColor[i],0.7);
            graphGammaFractionsCocktailBand[i]->SetFillStyle(4050);
            graphGammaFractionsCocktailBand[i]->Draw("3same");
            if(i<4)
                graphGammaFractionsCocktailBand[i]->Draw("LXsame");

        }
    }
    legendGammasRatio2->Draw("same");
    PutProcessLabelAndEnergyOnPlot(                 0.17, 0.22, 40, "", collisionSystempp2760GeV.Data(), "", 43, 0.03);
    PutALICESimulationLabel(                   0.17, 0.15, 40, 0.03, 1.25, 43);
    dummyHist->Draw("same,axis");

    canvasGammasRatio2->SaveAs(Form("%s/CocktailGammasRatioToAll_8.%s",outputDir.Data(),suffix.Data()));

    cout << __LINE__ << endl;


    for(Int_t j=0;j<11;j++){
        if(histoGammaFractionsCocktail[0][j]){
            for(Int_t i=0;i<14;i++){
                if(i==5||i==8){
                    histoGammaFractionsCocktailPartialSums[i][j]->Sumw2();
                    histoGammaFractionsCocktailPartialSums[i+1][j]->Sumw2();
                    histoGammaFractionsCocktailPartialSums[i][j]->Add(histoGammaFractionsCocktail[i+1][j]);
                }
                if(i==6||i==9)
                    histoGammaFractionsCocktailPartialSums[i][j] = NULL;
                if(histoGammaFractionsCocktailPartialSums[i][j])
                    graphGammaFractionsCocktailPartialSums[i][j] = new TGraphAsymmErrors(histoGammaFractionsCocktailPartialSums[i][j]);
                if(histoGammaFractionsCocktailPartialSums[i][j]&&j==5)
                    graphGammaFractionsCocktailPartialSumsBand[i] = new TGraphAsymmErrors(histoGammaFractionsCocktailPartialSums[i][j]);

            }
        }
    }
    for(Int_t i=0;i<14;i++){
        if(graphGammaFractionsCocktailPartialSums[i][5]){
            for(Int_t k = 0; k <= graphGammaFractionsCocktailPartialSums[i][5]->GetN()-1; k++){
                newBinErrorLow        = graphGammaFractionsCocktailPartialSums[i][5]->GetY()[k]-graphGammaFractionsCocktailPartialSums[i][5]->GetEYlow()[k];
                newBinErrorHigh       = graphGammaFractionsCocktailPartialSums[i][5]->GetY()[k]+graphGammaFractionsCocktailPartialSums[i][5]->GetEYhigh()[k];
                lowChanged=kFALSE;
                highChanged=kFALSE;
                for(Int_t j=0;j<11;j++){
                    if(graphGammaFractionsCocktailPartialSums[i][j]){
                        if((newBinErrorLow-graphGammaFractionsCocktailPartialSums[i][j]->GetY()[k]-graphGammaFractionsCocktailPartialSums[i][j]->GetEYlow()[k])>0){
                            newBinErrorLow = graphGammaFractionsCocktailPartialSums[i][j]->GetY()[k]-graphGammaFractionsCocktailPartialSums[i][j]->GetEYlow()[k];
                            lowChanged = kTRUE;
                        }
                        if((newBinErrorHigh-(graphGammaFractionsCocktailPartialSums[i][j]->GetY()[k]+graphGammaFractionsCocktailPartialSums[i][j]->GetEYlow()[k]))<0){
                            newBinErrorHigh = graphGammaFractionsCocktailPartialSums[i][j]->GetY()[k]+graphGammaFractionsCocktailPartialSums[i][j]->GetEYhigh()[k];
                            highChanged=kTRUE;
                        }
                    }
                }
                if(lowChanged)
                    graphGammaFractionsCocktailPartialSumsBand[i]->SetPointEYlow(k,     graphGammaFractionsCocktailPartialSums[i][5]->GetY()[k]-newBinErrorLow);
                if(highChanged)
                    graphGammaFractionsCocktailPartialSumsBand[i]->SetPointEYhigh(k,     newBinErrorHigh-graphGammaFractionsCocktailPartialSums[i][5]->GetY()[k]);
            }
            while(graphGammaFractionsCocktailPartialSumsBand[i]->GetX()[0] < 0.2) graphGammaFractionsCocktailPartialSumsBand[i]->RemovePoint(0);
        }
    }


    dummyHist->Draw();
    legendGammasRatio2->Clear();
    legendGammasRatio2                                 = GetAndSetLegend2(0.2, 0.96-(40*1.15*(2+1)/1200), 0.95, 0.96, 40, 6);
    legendGammasRatio2->SetHeader("#gamma from");
    for (Int_t i=0; i<14; i++) {
        if (histoGammaFractionsCocktailPartialSums[i][0]) {
            cout << i << " being plotted - particle " << fParticleLatexPartialSums[i].Data() << endl;
            DrawGammaSetMarker(             histoGammaFractionsCocktailPartialSums[i][5], cocktailMarker[i], 1, cocktailColorPartialSums[i],  cocktailColorPartialSums[i]);
            histoGammaFractionsCocktailPartialSums[i][5]->SetLineWidth(3);
            if(i!=8)
                legendGammasRatio2->AddEntry(   histoGammaFractionsCocktailPartialSums[i][5], fParticleLatexPartialSums[i].Data(), "l");
            // if(i<4)
            graphGammaFractionsCocktailPartialSumsBand[i]->SetFillColor(cocktailColorPartialSums[i]);
            // else
            // graphGammaFractionsCocktailPartialSumsBand[i]->SetFillColorAlpha(cocktailColorPartialSums[i],0.7);
            graphGammaFractionsCocktailPartialSumsBand[i]->SetFillStyle(4050);

            if(i==5)  // continue to draw rho above phi
                continue;

            graphGammaFractionsCocktailPartialSumsBand[i]->Draw("3same");
            if(i==7)
                graphGammaFractionsCocktailPartialSumsBand[5]->Draw("3same");

            if(i<4)
                graphGammaFractionsCocktailPartialSumsBand[i]->Draw("LXsame");
        }
    }
    legendGammasRatio2->AddEntry(   histoGammaFractionsCocktailPartialSums[8][5], fParticleLatexPartialSums[8].Data(), "l");
    legendGammasRatio2->Draw("same");
    PutProcessLabelAndEnergyOnPlot(                 0.17, 0.22, 40, "", collisionSystempp2760GeV.Data(), "", 43, 0.03);
    PutALICESimulationLabel(                   0.17, 0.15, 40, 0.03, 1.25, 43);
    dummyHist->Draw("same,axis");

    canvasGammasRatio2->SaveAs(Form("%s/CocktailGammasRatioToAllPartialSums_8.%s",outputDir.Data(),suffix.Data()));

    // ****************************************************************************************************************
    // ************************** Store final results including corr factors in 1 file ********************************
    // ****************************************************************************************************************
    TString nameOutputCommonFile    = Form("CombinedGammaResultPP2760GeV_%s.root", dateForOutput.Data());

    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    fCombResults.mkdir("Gamma2.76TeV");
    TDirectoryFile* directoryGamma = (TDirectoryFile*)fCombResults.Get("Gamma2.76TeV");
    fCombResults.cd("Gamma2.76TeV");

        // writing main results
        if (graphCombDRStat) graphCombDRStat->Write("graphRGammaCombStatErr");
        if (graphCombDRSys) graphCombDRSys->Write("graphRGammaCombSysErr");
        if (graphCombDRTot) graphCombDRTot->Write("graphRGammaCombTotErr");
        if (graphCombIncGammaStat) graphCombIncGammaStat->Write("graphInvYieldIncGammaStatErr");
        if (graphCombIncGammaSys) graphCombIncGammaSys->Write("graphInvYieldIncGammaSysErr");
        if (graphCombIncGammaTot) graphCombIncGammaTot->Write("graphInvYieldIncGammaTotErr");
        if (graphCombIncGammaStatUnshi) graphCombIncGammaStatUnshi->Write("graphInvYieldIncGammaStatErr_Unshifted");
        if (graphCombIncGammaSysUnshi) graphCombIncGammaSysUnshi->Write("graphInvYieldIncGammaSysErr_Unshifted");
        if (graphCombIncGammaTotUnshi) graphCombIncGammaTotUnshi->Write("graphInvYieldIncGammaTotErr_Unshifted");
        if (graphCombDirGammaSpectrumSystErr) graphCombDirGammaSpectrumSystErr->Write("graphInvYieldDirGammaSysErr");
        if (graphCombDirGammaSpectrumStatErr) graphCombDirGammaSpectrumStatErr->Write("graphInvYieldDirGammaStatErr");
        if (graphCombDirGammaSpectrumSumErrAr) graphCombDirGammaSpectrumSumErrAr->Write("graphInvYieldDirGammaSumErrAr");

        for (Int_t i = 0; i < 11; i++){
            if (graphIndGammaIncSys[i]) graphIndGammaIncSys[i]->Write(Form("graphInvYieldIncGamma%sSysErr",nameMeasGlobal[i].Data()));
            if (graphIndGammaIncStat[i]) graphIndGammaIncStat[i]->Write(Form("graphInvYieldIncGamma%sStatErr",nameMeasGlobal[i].Data()));
            if (histoIncGammaStatErr[i]) histoIncGammaStatErr[i]->Write(Form("histoInvYieldIncGamma%sStatErr_Unshifted",nameMeasGlobal[i].Data()));
            if (histoDRPi0FitStatErr[i]) histoDRPi0FitStatErr[i]->Write(Form("histoRGamma%sStatErr",nameMeasGlobal[i].Data()));
            if (graphDRPi0FitSysErr[i]) graphDRPi0FitSysErr[i]->Write(Form("graphRGamma%sSysErr",nameMeasGlobal[i].Data()));
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
