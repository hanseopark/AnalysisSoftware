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


void CombineGammaResultsPP8TeV(     TString inputFileNamePCM        = "CombinationInputGammaPP/InputPCMGammapp8000GeV.root",
                                    Bool_t havePHOS                 = kFALSE,
                                    TString inputFileNamePHOS       = "",
                                    Bool_t haveEMC                  = kTRUE,
                                    TString inputFileNameEMC        = "CombinationInputGammaPP/InputEMCGammapp8000GeV.root",
                                    Bool_t havePCMEMC               = kTRUE,
                                    TString inputFileNamePCMEMC     = "CombinationInputGammaPP/InputPCMEMCGammapp8000GeV.root",
                                    TString inputFileNameCocktail   = "CombinationInputGammaPP/GammaCocktailRatios.root",
                                    TString suffix                  = "pdf",
                                    TString fileNameCorrelations    = "",
                                    Bool_t enablepValueCalc         = kFALSE,
                                    Double_t confidenceLevelnSigma  = 1.28, //1.64 95% C.L. //1.28 //90% C.L.
                                    Bool_t useBinShifted            = kFALSE,
                                    TString inputFileFits           = ""
                            ){

    //*******************************************************************************************************************************************
    //*********************************************************** Set main style choices ********************************************************
    //*******************************************************************************************************************************************
    StyleSettingsThesis();
    SetPlotStyle();

    //*******************************************************************************************************************************************
    //********************************************* Create output directory and copy input files there ******************************************
    //*******************************************************************************************************************************************
    TString strBinShifted                           = "";
    if(useBinShifted){
      strBinShifted                           = "_BinShifted";
      cout << "using bin shifted input spectra" << endl;
    }
    TString dateForOutput                                       = ReturnDateStringForOutput();
    TString outputDir                                           = Form("%s/%s/CombineGammaMeasurementspp8TeV%s",suffix.Data(),dateForOutput.Data(),strBinShifted.Data());
    TString fileNameTheorypp8TeV                                = "ExternalInput/Theory/TheoryCompilationPP.root";

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCMGammapp8TeV.root", inputFileNamePCM.Data(), outputDir.Data()));
    if (havePHOS) gSystem->Exec(Form("cp %s %s/InputPHOSGammapp8TeV.root", inputFileNamePHOS.Data(), outputDir.Data()));
    if (havePCMEMC) gSystem->Exec(Form("cp %s %s/InputPCMEMCGammapp8TeV.root", inputFileNamePCMEMC.Data(), outputDir.Data()));
    if (haveEMC) gSystem->Exec(Form("cp %s %s/InputEMCGammapp8TeV.root", inputFileNameEMC.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/Theorypp8TeV.root", fileNameTheorypp8TeV.Data(), outputDir.Data()));

    TString nameFinalResDat                                     = Form("%s/CombinedResultsGamma_FitResults_%s.dat",outputDir.Data(),dateForOutput.Data());
    fstream fileFinalResults;
    fileFinalResults.open(nameFinalResDat.Data(), ios::out);
    
    TString addNamePi0Output                      = "";
    Bool_t writeFitsToPi0File                     = kFALSE;
    if (inputFileFits.CompareTo("") != 0){
      TFile* filePi0Comb                           = new TFile(inputFileFits);
      if (filePi0Comb->IsZombie()){
        writeFitsToPi0File                        = kFALSE;
      } else {
        writeFitsToPi0File                        = kTRUE;
        addNamePi0Output                          = Form("_Gamma_%s",dateForOutput.Data());
        inputFileFits.Resize(inputFileFits.Length()-5);
        gSystem->Exec(Form("cp %s.root %s%s.root",inputFileFits.Data(),inputFileFits.Data(),addNamePi0Output.Data()));
      }
    }
    //*******************************************************************************************************************************************
    //******************************************************* set ranges for plotting ***********************************************************
    //*******************************************************************************************************************************************
    Double_t doubleRatio[2];
    Double_t indMeasRatio[2];
//     Double_t incRatio[2];
    Double_t doubleRatioX[2];
    Double_t doubleRatioXpp[2];
    doubleRatio[0]      = 0.72 ;     doubleRatio[1]      = 1.55;
    indMeasRatio[0]     = 0.65;     indMeasRatio[1]     = 1.45;
//  incRatio[0]         = 0.0;      incRatio[1]         = 1.7;
    doubleRatioX[0]     = 0.7;      doubleRatioX[1]     = 24.5;
    doubleRatioXpp[0]   = 0.23;     doubleRatioXpp[1]   = 24.5  ;

    Color_t colorCocktailPi0                        = kRed+2;
    Color_t colorCocktailEta                        = kBlue+1;
    Color_t colorCocktailEtaP                       = kOrange+1;
    Color_t colorCocktailOmega                      = kYellow+2;
    Color_t colorCocktailPhi                        = kViolet;
    Color_t colorCocktailRho0                       = kAzure-2;
    Color_t colorCocktailSigma0                     = kGray+1;

    Color_t  colorNLOWerner                         = kAzure+2;
    Color_t  colorNLOWernerBand                     = kAzure-9;
    Color_t  colorJETPHOX                           = 807;
    Color_t  colorJETPHOXBand                       = 807;
    Color_t  colorPOWHEG                            = kRed+2;
    Color_t  colorPOWHEGBand                        = kRed-9;
    Color_t  colorNLOMcGill                         = kGreen+2;
    Color_t  colorNLOMcGillBand                     = kGreen-6;
    Style_t  styleMarkerNLOWerner                   = 24;
    Style_t  styleLineNLOWerner                     = 5;
    Style_t  styleMarkerJETPHOX                     = 23;
    Style_t  styleLineJETPHOX                       = 6;
    Style_t  styleMarkerPOWHEG                      = 25;
    Style_t  styleLinePOWHEG                        = 4;
    Style_t  styleLineMcGill                        = 7;
    Width_t  widthLineNLO                           = 2.;
    Int_t nParticles                                = 14;
    TString fParticle[14]                           = { "Pi0", "Eta", "omega", "EtaPrim", "rho0", "rho+", "rho-", "phi", "Delta0", "Delta+", "Sigma0","K0s", "K0l", "Lambda"};
    TString fParticleLatex[14]                      = { "#pi^{0}", "#eta", "#omega", "#eta'", "#rho^{0}", "#rho^{+}", "#rho^{-}", "#phi", "#Delta^{0}", "#Delta^{+}", "#Sigma^{0}", "K^{0}_{s}","K^{0}_{l}", "#Lambda"};
    TString fParticleLatexPartialSums[14]           = { "#pi^{0}", "#eta", "#omega", "#eta'", "#rho^{0}", "#rho^{#pm}", "", "#phi", "#Delta^{0/+}", "", "#Sigma^{0}", "K^{0}_{s}","K^{0}_{l}", "#Lambda"};
    Style_t     cocktailColor[14]                   = {kRed+2,kBlue+1,kYellow+2,kOrange+1,kAzure-2,kGreen-3,kRed-2,kViolet, kBlue-3, kTeal+9,kMagenta+2,kCyan+4,kViolet+4,kAzure-4};
    Color_t     cocktailMarker[14]                  = {20,    21, 25, 24,  20,      21,      24,25, 24,      25,21,        24,     25,       20};


    TString  nameMeasGlobal[11]                     = { "PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal",
                                                        "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged",
                                                        "PCMOtherDataset"};
    TString  nameMeasGlobalLabel[11]                = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC",
                                                        "PCM-Dalitz", "PHOS-Dalitz", "EMC-Dalitz", "EMC high pT", "mEMC",
                                                        "PCMOtherDataset"};
    TString  nameMeasGlobalLabelGamma[11]                = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM*",
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

    Color_t colorComb                               = GetColorDefaultColor("8TeV", "", "");
    Style_t markerStyleComb                         = GetDefaultMarkerStyle("8TeV", "", "");
    Size_t markerSizeComb                           = GetDefaultMarkerSize("8TeV", "", "");
    Color_t colorCombpp8TeV                         = kBlack; // GetColorDefaultColor("8TeV", "", "");
    Color_t colorCombpp8TeVBox                      = kGray+2; //GetColorDefaultColor("8TeV", "", "", kTRUE);
    Style_t markerStyleCombpp8TeV                   = 20; //GetDefaultMarkerStyle("8TeV", "", "");
    Size_t markerSizeCombpp8TeV                     = 1.8; //GetDefaultMarkerSize("8TeV", "", "");

    Width_t widthLinesBoxes                         = 1.4;
    Width_t widthCommonFit                          = 2.4;

    TString collisionSystempp8TeV                   = "pp, #sqrt{#it{s}} = 8 TeV";
    TString textALICE                               = "ALICE";
      
    cout << "Setting Gamma binning" << endl;
    Double_t xPtLimitsGamma[100]                    = {0};
    Int_t maxNBinsGamma                             = GetBinning( xPtLimitsGamma, "Gamma", "8TeV", 0 );
    for (Int_t i = 0; i< maxNBinsGamma; i++){
        cout << i << ": "<< xPtLimitsGamma[i] <<" - " << xPtLimitsGamma[i+1]<< ", " <<endl;
    }


    TH1D* histoDRPi0FitStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors*  graphDRPi0FitSysErr[11]     = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TH1D* histoDRStatErr[11]                  = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
    TGraphAsymmErrors*  graphDRSysErr[11]     = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
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
    TH1D* histoPileupCorr[11]                       = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};;
    
    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from 8TeV PCM file *************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMGammapp8TeV                    = new TFile( inputFileNamePCM.Data());
    //________________________________________________ Load PCM 8TeV _________________________________________________________________________
    TDirectory* directoryPCMGammapp8TeV          = (TDirectory*)filePCMGammapp8TeV->Get("Gamma_pp8TeV");
        histoDRPi0FitStatErr[0]                         = (TH1D*) directoryPCMGammapp8TeV->Get("DoubleRatioPi0FitStatError");
        graphDRPi0FitSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapp8TeV->Get("DoubleRatioPi0FitSystError");
        histoDRStatErr[0]                               = (TH1D*) directoryPCMGammapp8TeV->Get(Form("DoubleRatioStatError%s",strBinShifted.Data()));
        histoDRStatErr[0]->SetBinContent(2,histoDRStatErr[0]->GetBinContent(2)*1.22);
        histoDRStatErr[0]->SetBinError(2,histoDRStatErr[0]->GetBinError(2)*1.22);
        graphDRSysErr[0]                                = (TGraphAsymmErrors*) directoryPCMGammapp8TeV->Get(Form("DoubleRatioSystError%s",strBinShifted.Data()));
        graphDRSysErr[0]->SetPoint(0,graphDRSysErr[0]->GetX()[0],graphDRSysErr[0]->GetY()[0]*1.22);
        graphDRSysErr[0]->SetPointEYlow(0,graphDRSysErr[0]->GetEYlow()[0]*1.22);
        graphDRSysErr[0]->SetPointEYhigh(0,graphDRSysErr[0]->GetEYhigh()[0]*1.22);        
        histoIncGammaRatioStatErr[0]                    = (TH1D*) directoryPCMGammapp8TeV->Get(Form("IncRatioStatError%s",strBinShifted.Data()));
        graphIncGammaRatioSysErr[0]                     = (TGraphAsymmErrors*) directoryPCMGammapp8TeV->Get(Form("IncRatioSystError%s",strBinShifted.Data()));
        histoConvProb[0]                                = (TH1D*) directoryPCMGammapp8TeV->Get("GammaConversionProbability");
        histoEffi[0]                                    = (TH1D*) directoryPCMGammapp8TeV->Get("GammaRecoEfficiency");
        histoEffiMCPt[0]                                = (TH1D*) directoryPCMGammapp8TeV->Get("GammaRecoEfficiencyMCPt");
        histoPurity[0]                                  = (TH1D*) directoryPCMGammapp8TeV->Get("GammaTruePurity");
        histoResolCorr[0]                               = (TH1D*) directoryPCMGammapp8TeV->Get("GammaResolCorr");
        histoIncGammaStatErr[0]                         = (TH1D*) directoryPCMGammapp8TeV->Get(Form("IncGammaStatError%s",strBinShifted.Data()));
        graphIncGammaSysErr[0]                          = (TGraphAsymmErrors*) directoryPCMGammapp8TeV->Get(Form("IncGammaSystError%s",strBinShifted.Data()));
        histoPileupCorr[0]                              = (TH1D*) directoryPCMGammapp8TeV->Get("PileUpCorrectionFactor");
        histoEffSecCorr[0][0]                           = (TH1D*) directoryPCMGammapp8TeV->Get("GammaEffectiveSecondaryCorr_K0s");
        histoEffSecCorr[1][0]                           = (TH1D*) directoryPCMGammapp8TeV->Get("GammaEffectiveSecondaryCorr_K0l");
        histoEffSecCorr[2][0]                           = (TH1D*) directoryPCMGammapp8TeV->Get("GammaEffectiveSecondaryCorr_Lambda");
        histoEffSecCorr[3][0]                           = (TH1D*) directoryPCMGammapp8TeV->Get("GammaEffectiveSecondaryCorr_Rest");
        
    //*******************************************************************************************************************************************
    //*********************************************** Load PCMEMC histograms from 8TeV PCM file **********************************************
    //*******************************************************************************************************************************************
    if (havePCMEMC){
        TFile* filePCMEMCGammapp8TeV                 = new TFile( inputFileNamePCMEMC.Data());
        //________________________________________________ Load PCM-EMC 8TeV _________________________________________________________________________
        TDirectory* directoryPCMEMCGammapp8TeV       = (TDirectory*)filePCMEMCGammapp8TeV->Get("Gamma_pp8TeV");
            histoDRPi0FitStatErr[4]                         = (TH1D*) directoryPCMEMCGammapp8TeV->Get("DoubleRatioPi0FitStatError");
            graphDRPi0FitSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapp8TeV->Get("DoubleRatioPi0FitSystError");
            histoDRStatErr[4]                         = (TH1D*) directoryPCMEMCGammapp8TeV->Get(Form("DoubleRatioStatError%s",strBinShifted.Data()));
            graphDRSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapp8TeV->Get(Form("DoubleRatioSystError%s",strBinShifted.Data()));
            histoIncGammaRatioStatErr[4]                    = (TH1D*) directoryPCMEMCGammapp8TeV->Get(Form("IncRatioStatError%s",strBinShifted.Data()));
            graphIncGammaRatioSysErr[4]                     = (TGraphAsymmErrors*) directoryPCMEMCGammapp8TeV->Get(Form("IncRatioSystError%s",strBinShifted.Data()));
            histoConvProb[4]                                = (TH1D*) directoryPCMEMCGammapp8TeV->Get("GammaConversionProbability");
            histoEffi[4]                                    = (TH1D*) directoryPCMEMCGammapp8TeV->Get("GammaRecoEfficiency");
            histoEffiMCPt[4]                                = (TH1D*) directoryPCMEMCGammapp8TeV->Get("GammaRecoEfficiencyMCPt");
            histoPurity[4]                                  = (TH1D*) directoryPCMEMCGammapp8TeV->Get("GammaTruePurity");
            histoResolCorr[4]                               = (TH1D*) directoryPCMEMCGammapp8TeV->Get("GammaResolCorr");
            histoIncGammaStatErr[4]                         = (TH1D*) directoryPCMEMCGammapp8TeV->Get(Form("IncGammaStatError%s",strBinShifted.Data()));
            graphIncGammaSysErr[4]                          = (TGraphAsymmErrors*) directoryPCMEMCGammapp8TeV->Get(Form("IncGammaSystError%s",strBinShifted.Data()));
            histoPileupCorr[4]                              = (TH1D*) directoryPCMEMCGammapp8TeV->Get("PileUpCorrectionFactor");
            histoEffSecCorr[0][4]                           = (TH1D*) directoryPCMEMCGammapp8TeV->Get("GammaEffectiveSecondaryCorr_K0s");
            histoEffSecCorr[1][4]                           = (TH1D*) directoryPCMEMCGammapp8TeV->Get("GammaEffectiveSecondaryCorr_K0l");
            histoEffSecCorr[2][4]                           = (TH1D*) directoryPCMEMCGammapp8TeV->Get("GammaEffectiveSecondaryCorr_Lambda");
            histoEffSecCorr[3][4]                           = (TH1D*) directoryPCMEMCGammapp8TeV->Get("GammaEffectiveSecondaryCorr_Rest");
    }
    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from 8TeV EMC file *************************************************
    //*******************************************************************************************************************************************
    if (haveEMC){
        TFile* fileEMCGammapp8TeV                    = new TFile( inputFileNameEMC.Data());
        //________________________________________________ Load EMC 8TeV _________________________________________________________________________
        TDirectory* directoryEMCGammapp8TeV          = (TDirectory*)fileEMCGammapp8TeV->Get("Gamma_pp8TeV");
            histoDRPi0FitStatErr[2]                         = (TH1D*) directoryEMCGammapp8TeV->Get("DoubleRatioPi0FitStatError");
            graphDRPi0FitSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapp8TeV->Get("DoubleRatioPi0FitSystError");
            histoDRStatErr[2]                         = (TH1D*) directoryEMCGammapp8TeV->Get(Form("DoubleRatioStatError%s",strBinShifted.Data()));
            graphDRSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapp8TeV->Get(Form("DoubleRatioSystError%s",strBinShifted.Data()));
            histoIncGammaRatioStatErr[2]                    = (TH1D*) directoryEMCGammapp8TeV->Get(Form("IncRatioStatError%s",strBinShifted.Data()));
            graphIncGammaRatioSysErr[2]                     = (TGraphAsymmErrors*) directoryEMCGammapp8TeV->Get(Form("IncRatioSystError%s",strBinShifted.Data()));
            histoEffi[2]                                    = (TH1D*) directoryEMCGammapp8TeV->Get("GammaRecoEfficiency");
            histoEffiMCPt[2]                                = (TH1D*) directoryEMCGammapp8TeV->Get("GammaRecoEfficiencyMCPt");
            histoPurity[2]                                  = (TH1D*) directoryEMCGammapp8TeV->Get("GammaTruePurity");
            histoResolCorr[2]                               = (TH1D*) directoryEMCGammapp8TeV->Get("GammaResolCorr");
            histoIncGammaStatErr[2]                         = (TH1D*) directoryEMCGammapp8TeV->Get(Form("IncGammaStatError%s",strBinShifted.Data()));
            graphIncGammaSysErr[2]                          = (TGraphAsymmErrors*) directoryEMCGammapp8TeV->Get(Form("IncGammaSystError%s",strBinShifted.Data()));
            histoEffSecCorr[0][2]                           = (TH1D*) directoryEMCGammapp8TeV->Get("GammaEffectiveSecondaryCorr_K0s");
            histoEffSecCorr[1][2]                           = (TH1D*) directoryEMCGammapp8TeV->Get("GammaEffectiveSecondaryCorr_K0l");
            histoEffSecCorr[2][2]                           = (TH1D*) directoryEMCGammapp8TeV->Get("GammaEffectiveSecondaryCorr_Lambda");
            histoEffSecCorr[3][2]                           = (TH1D*) directoryEMCGammapp8TeV->Get("GammaEffectiveSecondaryCorr_Rest");
    }
    //*******************************************************************************************************************************************
    //************************************************ Load theory curves from external input ***************************************************
    //*******************************************************************************************************************************************
    Double_t xmaxNLO = 20;
    TFile* fileTheory                               = new TFile( fileNameTheorypp8TeV.Data());
        TGraphAsymmErrors* graphTheoryNLODRpp8TeV    = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonNLOVogelsangInvYieldINT7_pp8TeV_CT10_ALICECocktail");
        TGraph* graphTheoryNLODRpp8TeVCenter         = (TGraph*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonNLOVogelsangInvYieldINT7_pp8TeV_CT10_ALICECocktail_Center");
        TGraphAsymmErrors* graphTheoryNLOpp8TeV      = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphDirectPhotonNLOVogelsangInvYield_8TeV");
        TGraphAsymmErrors* graphTheoryNLOpp8TeVCenter      = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphDirectPhotonNLOVogelsangInvYield_8TeV");
        TGraphAsymmErrors* graphTheoryJETPHOXDRpp8TeV    = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonJETPHOXInvYieldINT7_pp8TeV_ALICECocktail");
        TGraph* graphTheoryJETPHOXDRpp8TeVCenter         = (TGraph*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonJETPHOXInvYieldINT7_pp8TeV_ALICECocktail_Center");
        TGraphAsymmErrors* graphTheoryJETPHOXpp8TeV      = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphDirectPhotonJETPHOXInvYield_8TeV");
        TGraphAsymmErrors* graphTheoryJETPHOXpp8TeVCenter      = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphDirectPhotonJETPHOXInvYield_8TeV");
        TGraphAsymmErrors* graphTheoryPOWHEGDRpp8TeV    = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonPOWHEGInvYieldINT7_pp8TeV_ALICECocktail");
        TGraph* graphTheoryPOWHEGDRpp8TeVCenter         = (TGraph*) fileTheory->Get("DirectPhoton/graphRGammaDirectPhotonPOWHEGInvYieldINT7_pp8TeV_ALICECocktail_Center");
        TGraphAsymmErrors* graphTheoryPOWHEGpp8TeV      = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphDirectPhotonPOWHEGInvYield_8TeV");
        TGraphAsymmErrors* graphTheoryPOWHEGpp8TeVCenter      = (TGraphAsymmErrors*) fileTheory->Get("DirectPhoton/graphDirectPhotonPOWHEGInvYield_8TeV");
        TGraph* graphTheoryNLODRpp8TeVPaquett  = (TGraph*) fileTheory->Get("DirectPhoton/graphRGammaPaquett_pp8TeV_ALICECocktail");
        TGraph* graphTheoryNLOp8TeVPaquett    = (TGraph*) fileTheory->Get("DirectPhoton/graphPromptPhotonNLOPaquett_8TeV");
        while(graphTheoryNLODRpp8TeVPaquett->GetX()[graphTheoryNLODRpp8TeVPaquett->GetN()-1] > xmaxNLO) graphTheoryNLODRpp8TeVPaquett->RemovePoint(graphTheoryNLODRpp8TeVPaquett->GetN()-1);
        while(graphTheoryNLOp8TeVPaquett->GetX()[graphTheoryNLOp8TeVPaquett->GetN()-1] > xmaxNLO) graphTheoryNLOp8TeVPaquett->RemovePoint(graphTheoryNLOp8TeVPaquett->GetN()-1);

        while(graphTheoryNLODRpp8TeV->GetX()[0] < 1.5) graphTheoryNLODRpp8TeV->RemovePoint(0);
        while(graphTheoryNLODRpp8TeVCenter->GetX()[0] < 1.5) graphTheoryNLODRpp8TeVCenter->RemovePoint(0);
        while(graphTheoryNLOpp8TeV->GetX()[0] < 1.5) graphTheoryNLOpp8TeV->RemovePoint(0);
        while(graphTheoryNLOpp8TeVCenter->GetX()[0] < 1.5) graphTheoryNLOpp8TeVCenter->RemovePoint(0);
        while(graphTheoryNLODRpp8TeV->GetX()[graphTheoryNLODRpp8TeV->GetN()-1] > xmaxNLO) graphTheoryNLODRpp8TeV->RemovePoint(graphTheoryNLODRpp8TeV->GetN()-1);
        while(graphTheoryNLODRpp8TeVCenter->GetX()[graphTheoryNLODRpp8TeVCenter->GetN()-1] > xmaxNLO) graphTheoryNLODRpp8TeVCenter->RemovePoint(graphTheoryNLODRpp8TeVCenter->GetN()-1);
        while(graphTheoryNLOpp8TeV->GetX()[graphTheoryNLOpp8TeV->GetN()-1] > xmaxNLO) graphTheoryNLOpp8TeV->RemovePoint(graphTheoryNLOpp8TeV->GetN()-1);
        while(graphTheoryNLOpp8TeVCenter->GetX()[graphTheoryNLOpp8TeVCenter->GetN()-1] > xmaxNLO) graphTheoryNLOpp8TeVCenter->RemovePoint(graphTheoryNLOpp8TeVCenter->GetN()-1);
        
        while(graphTheoryJETPHOXpp8TeV->GetX()[0] < 1.5) graphTheoryJETPHOXpp8TeV->RemovePoint(0);
        while(graphTheoryJETPHOXpp8TeVCenter->GetX()[0] < 1.5) graphTheoryJETPHOXpp8TeVCenter->RemovePoint(0);
        while(graphTheoryJETPHOXDRpp8TeV->GetX()[0] < 1.5) graphTheoryJETPHOXDRpp8TeV->RemovePoint(0);
        while(graphTheoryJETPHOXDRpp8TeVCenter->GetX()[0] < 1.5) graphTheoryJETPHOXDRpp8TeVCenter->RemovePoint(0);
        while(graphTheoryJETPHOXpp8TeV->GetX()[graphTheoryJETPHOXpp8TeV->GetN()-1] > xmaxNLO) graphTheoryJETPHOXpp8TeV->RemovePoint(graphTheoryJETPHOXpp8TeV->GetN()-1);
        while(graphTheoryJETPHOXpp8TeVCenter->GetX()[graphTheoryJETPHOXpp8TeVCenter->GetN()-1] > xmaxNLO) graphTheoryJETPHOXpp8TeVCenter->RemovePoint(graphTheoryJETPHOXpp8TeVCenter->GetN()-1);
        while(graphTheoryJETPHOXDRpp8TeV->GetX()[graphTheoryJETPHOXDRpp8TeV->GetN()-1] > xmaxNLO) graphTheoryJETPHOXDRpp8TeV->RemovePoint(graphTheoryJETPHOXDRpp8TeV->GetN()-1);
        while(graphTheoryJETPHOXDRpp8TeVCenter->GetX()[graphTheoryJETPHOXDRpp8TeVCenter->GetN()-1] > xmaxNLO) graphTheoryJETPHOXDRpp8TeVCenter->RemovePoint(graphTheoryJETPHOXDRpp8TeVCenter->GetN()-1);
        
        while(graphTheoryPOWHEGpp8TeV->GetX()[0] < 1.5) graphTheoryPOWHEGpp8TeV->RemovePoint(0);
        while(graphTheoryPOWHEGpp8TeVCenter->GetX()[0] < 1.5) graphTheoryPOWHEGpp8TeVCenter->RemovePoint(0);
        while(graphTheoryPOWHEGDRpp8TeV->GetX()[0] < 1.5) graphTheoryPOWHEGDRpp8TeV->RemovePoint(0);
        while(graphTheoryPOWHEGDRpp8TeVCenter->GetX()[0] < 1.5) graphTheoryPOWHEGDRpp8TeVCenter->RemovePoint(0);
        while(graphTheoryPOWHEGpp8TeV->GetX()[graphTheoryPOWHEGpp8TeV->GetN()-1] > xmaxNLO) graphTheoryPOWHEGpp8TeV->RemovePoint(graphTheoryPOWHEGpp8TeV->GetN()-1);
        while(graphTheoryPOWHEGpp8TeVCenter->GetX()[graphTheoryPOWHEGpp8TeVCenter->GetN()-1] > xmaxNLO) graphTheoryPOWHEGpp8TeVCenter->RemovePoint(graphTheoryPOWHEGpp8TeVCenter->GetN()-1);
        while(graphTheoryPOWHEGDRpp8TeV->GetX()[graphTheoryPOWHEGDRpp8TeV->GetN()-1] > xmaxNLO) graphTheoryPOWHEGDRpp8TeV->RemovePoint(graphTheoryPOWHEGDRpp8TeV->GetN()-1);
        while(graphTheoryPOWHEGDRpp8TeVCenter->GetX()[graphTheoryPOWHEGDRpp8TeVCenter->GetN()-1] > xmaxNLO) graphTheoryPOWHEGDRpp8TeVCenter->RemovePoint(graphTheoryPOWHEGDRpp8TeVCenter->GetN()-1);

        
        TGraphAsymmErrors* dummyJETPHOXforLegend    = new TGraphAsymmErrors(1);
        DrawGammaSetMarkerTGraphAsym(dummyJETPHOXforLegend , 2, styleLineJETPHOX, colorJETPHOX, colorJETPHOX, widthLinesBoxes, kTRUE, colorJETPHOX,kTRUE);
        dummyJETPHOXforLegend->SetLineStyle(styleLineJETPHOX);
        dummyJETPHOXforLegend->SetLineWidth(2);
        TGraphAsymmErrors* dummyPOWHEGforLegend    = new TGraphAsymmErrors(1);
        DrawGammaSetMarkerTGraphAsym(dummyPOWHEGforLegend , 2, styleLinePOWHEG, colorPOWHEG, colorPOWHEG, widthLinesBoxes, kTRUE, colorPOWHEG,kTRUE);
        dummyPOWHEGforLegend->SetLineStyle(styleLinePOWHEG);
        dummyPOWHEGforLegend->SetLineWidth(2);
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
    Int_t offSetsGammaSys[11]       = { 1,  0,  6,  0,  4,
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
            cout << "DR rel stat for input  " << i << endl;
            statErrorRelCollectionDR[i]->Print();
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
            cout << "DR rel sys for input  " << i << endl;
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
                                                                                    fileNameDROutputWeighting, "8TeV", "RGamma", kFALSE,
                                                                                    NULL, fileNameCorrelations );
// return;

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

    TLatex *labelWeightsEnergy      = new TLatex(0.95,0.20,collisionSystempp8TeV.Data());
    SetStyleTLatex( labelWeightsEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelWeightsEnergy->Draw();
    TLatex *labelWeightsDR         = new TLatex(0.95,0.15,"R_{#gamma}");
    SetStyleTLatex( labelWeightsDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelWeightsDR->Draw();

    DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/DR_Weights_8.%s",outputDir.Data(),suffix.Data()));
    canvasWeights->SaveAs(Form("%s/DR_Weights_8.pdf",outputDir.Data()));

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

    TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystempp8TeV.Data());
    SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrEnergy->Draw();
    TLatex *labelRelSysErrDR       = new TLatex(0.15,0.85,"R_{#gamma}");
    SetStyleTLatex( labelRelSysErrDR, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DR_RelSysErr_8.%s",outputDir.Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/DR_RelSysErr_8.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();

    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0,50.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetYaxis()->SetRangeUser(0,35.5);
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

    TLatex *labelRelStatErrEnergy   = new TLatex(0.95,0.89,collisionSystempp8TeV.Data());
    SetStyleTLatex( labelRelStatErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrEnergy->Draw();
    TLatex *labelRelStatErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
    SetStyleTLatex( labelRelStatErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrDR->Draw();

    canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr_8.%s",outputDir.Data(),suffix.Data()));
    canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr_8.pdf",outputDir.Data()));


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
        histo2DRelErr->GetYaxis()->SetRangeUser(0,35.5);
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

        TLegend* legendRelTotErr       = GetAndSetLegend2(0.22, 0.92-(0.035*3), 0.45, 0.92, 32);
        legendRelTotErr->AddEntry(graphCombDRRelTot,"tot","p");
        legendRelTotErr->AddEntry(graphCombDRRelStat,"stat","l");
        legendRelTotErr->AddEntry(graphCombDRRelSys,"sys","l");
        legendRelTotErr->Draw();

        TLatex *labelRelTotErrEnergy   = new TLatex(0.95,0.89,collisionSystempp8TeV.Data());
        SetStyleTLatex( labelRelTotErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
        SetStyleTLatex( labelRelTotErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DR_Reldecomp_8.%s",outputDir.Data(),suffix.Data()));

    //*******************************************************************************************************************************************
    //*********************************************** Combining Rgamma ratios unfitted **********************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsGammaNonFit[11]        = { 0,  0,  0,  0,  0,
                                            0,  0,  0,  0,  0,
                                            0};
    Int_t offSetsGammaNonFitSys[11]     = { 1,  0,  6,  0,  4,
                                            0,  0,  0,  0,  0,
                                            0};


    TGraphAsymmErrors* statErrorRelCollectionDRNonFit[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionDRNonFit[i]        = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (histoDRStatErr[i]){
            statErrorRelCollectionDRNonFit[i]    = new TGraphAsymmErrors(histoDRStatErr[i]);
            while (statErrorRelCollectionDRNonFit[i]->GetY()[0] == 0) statErrorRelCollectionDRNonFit[i]->RemovePoint(0);
            while (statErrorRelCollectionDRNonFit[i]->GetY()[statErrorRelCollectionDRNonFit[i]->GetN()-1] == 0) statErrorRelCollectionDRNonFit[i]->RemovePoint(statErrorRelCollectionDRNonFit[i]->GetN()-1);
            statErrorRelCollectionDRNonFit[i]    = CalculateRelErrUpAsymmGraph( statErrorRelCollectionDRNonFit[i], Form("relativeStatErrorDR_%s", nameMeasGlobal[i].Data()));
            cout << "DR nonfit rel stat for input " << i << endl;
            statErrorRelCollectionDRNonFit[i]->Print();
        }
    }

    TGraphAsymmErrors* sysErrorRelCollectionDRNonFit[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollectionDRNonFit[i]         = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        cout << i << endl;
        if (graphDRSysErr[i]){
            sysErrorRelCollectionDRNonFit[i]     = (TGraphAsymmErrors*)graphDRSysErr[i]->Clone(Form("relativeSysErrorDR_%s", nameMeasGlobal[i].Data()));
            sysErrorRelCollectionDRNonFit[i]->Print();
            while (sysErrorRelCollectionDRNonFit[i]->GetY()[0] == 0) sysErrorRelCollectionDRNonFit[i]->RemovePoint(0);
            while (sysErrorRelCollectionDRNonFit[i]->GetY()[sysErrorRelCollectionDRNonFit[i]->GetN()-1] == 0) sysErrorRelCollectionDRNonFit[i]->RemovePoint(sysErrorRelCollectionDRNonFit[i]->GetN()-1);
            sysErrorRelCollectionDRNonFit[i]     = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionDRNonFit[i], Form("relativeSysErrorDR_%s", nameMeasGlobal[i].Data()));
            cout << "DR nonfit rel sys for input " << i << endl;
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
    TString fileNameDROutputWeightingNonFit       = Form("%s/DR_WeightingMethod.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombDRNonFitStat      = NULL;
    TGraphAsymmErrors* graphCombDRNonFitSys       = NULL;
    TGraphAsymmErrors* graphCombDRNonFitTot       = CombinePtPointsSpectraFullCorrMat(    histoDRStatErr,    graphDRSysErr,
                                                                                    xPtLimitsGamma, maxNBinsGamma,
                                                                                    offSetsGammaNonFit, offSetsGammaNonFitSys,
                                                                                    graphCombDRNonFitStat, graphCombDRNonFitSys,
                                                                                    fileNameDROutputWeightingNonFit, "8TeV", "RGamma", kFALSE,
                                                                                    NULL, fileNameCorrelations );
// return;

    if (graphCombDRNonFitTot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    while (graphCombDRNonFitStat->GetX()[0] < 0.3){
        graphCombDRNonFitStat->RemovePoint(0);
    }
    while (graphCombDRNonFitTot->GetX()[0] < 0.3){
        graphCombDRNonFitTot->RemovePoint(0);
    }
    while (graphCombDRNonFitSys->GetX()[0] < 0.3){
        graphCombDRNonFitSys->RemovePoint(0);
    }
    graphCombDRNonFitTot->Print();

    // Reading weights from output file for plotting
    // ifstream fileWeightsDRRead;
    fileWeightsDRRead.open(fileNameDROutputWeightingNonFit,ios_base::in);
    cout << "reading" << fileNameDROutputWeightingNonFit << endl;
    // xValuesDRRead[50];
    // weightsDRRead[11][50];
    nMeasSetDR               = 3;
    nPtBinsDRRead            = 0;
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
        graphWeightsDRNonFit[availableDRMeas[i]]                        = new TGraph(nPtBinsDRRead,xValuesDRRead,weightsDRRead[availableDRMeas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsDRRead; n++){
            if (graphWeightsDRNonFit[availableDRMeas[i]]->GetY()[bin] == 0) graphWeightsDRNonFit[availableDRMeas[i]]->RemovePoint(bin);
            else bin++;
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plotting weights *********************************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 900*0.04;

    canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.01, 0.01, 0.09);
    canvasWeights->SetLogx();

    histo2DDRWeights;
    histo2DDRWeights = new TH2F("histo2DDRWeights","histo2DDRWeights",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DDRWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DDRWeights->GetXaxis()->SetMoreLogLabels();
    histo2DDRWeights->GetXaxis()->SetNoExponent();
    histo2DDRWeights->Draw("copy");

    histo2DDRWeights->Draw("copy");
    legendWeightsDR   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetDR+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDR; i++){
        DrawGammaSetMarkerTGraph(graphWeightsDRNonFit[availableDRMeas[i]], markerStyleDet[availableDRMeas[i]], markerSizeDet[availableDRMeas[i]], colorDet[availableDRMeas[i]] , colorDet[availableDRMeas[i]]);
        graphWeightsDRNonFit[availableDRMeas[i]]->Draw("p,same,z");
        legendWeightsDR->AddEntry(graphWeightsDRNonFit[availableDRMeas[i]],nameMeasGlobalLabel[availableDRMeas[i]],"p");
    }
    legendWeightsDR->Draw();

    labelWeightsEnergy      = new TLatex(0.95,0.20,collisionSystempp8TeV.Data());
    SetStyleTLatex( labelWeightsEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelWeightsEnergy->Draw();
    labelWeightsDR         = new TLatex(0.95,0.15,"R_{#gamma}");
    SetStyleTLatex( labelWeightsDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelWeightsDR->Draw();

    DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/DR_Weights_NonFit_8.%s",outputDir.Data(),suffix.Data()));
    canvasWeights->SaveAs(Form("%s/DR_Weights_NonFit_8.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();

    histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0,50.0);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,30.5);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetNoExponent();
    histo2DRelSysErr->Draw("copy");

    legendRelSysErr2       = GetAndSetLegend2(0.62, 0.92-(0.04*(nMeasSetDR+1)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDR; i++){
        cout << "sys\t" << nameMeasGlobalLabel[availableDRMeas[i]] << endl;
        DrawGammaSetMarkerTGraph(sysErrorRelCollectionDRNonFit[availableDRMeas[i]], markerStyleDet[availableDRMeas[i]], markerSizeDet[availableDRMeas[i]], colorDet[availableDRMeas[i]],
                                 colorDet[availableDRMeas[i]]);
        sysErrorRelCollectionDRNonFit[availableDRMeas[i]]->Draw("p,same,z");
        sysErrorRelCollectionDRNonFit[availableDRMeas[i]]->Print();
        legendRelSysErr2->AddEntry(sysErrorRelCollectionDRNonFit[availableDRMeas[i]],nameMeasGlobalLabel[availableDRMeas[i]],"p");
    }
    legendRelSysErr2->Draw();

    labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystempp8TeV.Data());
    SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrEnergy->Draw();
    labelRelSysErrDR       = new TLatex(0.15,0.85,"R_{#gamma}");
    SetStyleTLatex( labelRelSysErrDR, textSizeLabelsPixel, 4, 1, 43);
    labelRelSysErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DR_RelSysErr_NonFit_8.%s",outputDir.Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/DR_RelSysErr_NonFit_8.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();

    histo2DRelStatErr;
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0,50.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetYaxis()->SetRangeUser(0,30.5);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetNoExponent();
    histo2DRelStatErr->Draw("copy");
    legendRelStatErr2       = GetAndSetLegend2(0.14, 0.92-(0.04*(nMeasSetDR+1)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetDR; i++){
        DrawGammaSetMarkerTGraph(statErrorRelCollectionDRNonFit[availableDRMeas[i]], markerStyleDet[availableDRMeas[i]], markerSizeDet[availableDRMeas[i]], colorDet[availableDRMeas[i]],
                                 colorDet[availableDRMeas[i]]);
        statErrorRelCollectionDRNonFit[availableDRMeas[i]]->Draw("p,same,z");
        legendRelStatErr2->AddEntry(statErrorRelCollectionDRNonFit[availableDRMeas[i]],nameMeasGlobalLabel[availableDRMeas[i]],"p");
    }
    legendRelStatErr2->Draw();

    labelRelStatErrEnergy   = new TLatex(0.95,0.89,collisionSystempp8TeV.Data());
    SetStyleTLatex( labelRelStatErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrEnergy->Draw();
    labelRelStatErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
    SetStyleTLatex( labelRelStatErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRelStatErrDR->Draw();

    canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr_NonFit_8.%s",outputDir.Data(),suffix.Data()));
    canvasRelStatErr->SaveAs(Form("%s/DR_RelStatErr_NonFit_8.pdf",outputDir.Data()));


    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TGraphAsymmErrors* graphCombDRNonFitRelStat       = CalculateRelErrUpAsymmGraph( graphCombDRNonFitStat, "relativeStatErrorDR");
    TGraphAsymmErrors* graphCombDRNonFitRelSys        = CalculateRelErrUpAsymmGraph( graphCombDRNonFitSys, "relativeSysErrorDR");
    TGraphAsymmErrors* graphCombDRNonFitRelTot        = CalculateRelErrUpAsymmGraph( graphCombDRNonFitTot, "relativeTotalErrorDR");

    canvasRelSysErr->cd();
        histo2DRelErr                    = new TH2F("histo2DRelErr","histo2DRelErr",11000,doubleRatioXpp[0], doubleRatioXpp[1],1000,0,50.0);
        SetStyleHistoTH2ForGraphs(histo2DRelErr, "#it{p}_{T} (GeV/#it{c})","Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DRelErr->GetYaxis()->SetRangeUser(0,30.5);
        histo2DRelErr->GetXaxis()->SetMoreLogLabels();
        histo2DRelErr->GetXaxis()->SetNoExponent();
        histo2DRelErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombDRNonFitRelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombDRNonFitRelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombDRNonFitRelSys->SetLineStyle(7);
        graphCombDRNonFitRelSys->Draw("l,x0,same,e1");

        legendRelTotErr       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
        legendRelTotErr->AddEntry(graphCombDRNonFitRelTot,"tot","p");
        legendRelTotErr->AddEntry(graphCombDRNonFitRelStat,"stat","l");
        legendRelTotErr->AddEntry(graphCombDRNonFitRelSys,"sys","l");
        legendRelTotErr->Draw();

        labelRelTotErrEnergy   = new TLatex(0.95,0.89,collisionSystempp8TeV.Data());
        SetStyleTLatex( labelRelTotErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrEnergy->Draw();
        labelRelTotErrDR      = new TLatex(0.95,0.85,"R_{#gamma}");
        SetStyleTLatex( labelRelTotErrDR, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrDR->Draw();

    canvasRelSysErr->SaveAs(Form("%s/DR_Reldecomp_NonFit_8.%s",outputDir.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    cout << "DR - tot" << endl;
    graphCombDRNonFitRelTot->Print();
    cout << "DR - stat" << endl;
    graphCombDRNonFitRelStat->Print();
    cout << "DR - sys" << endl;
    graphCombDRNonFitRelSys->Print();

    //*******************************************************************************************************************************************
    //*********************************************** Combining IncGamma spectra ****************************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};

    Int_t offSetsIncGamma[11]       = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetsIncGammaSys[11]    = { 1,  0,  6,  0,  4,
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
            cout << "gamma rel stat for input " << i << endl;
            statErrorRelCollectionIncGamma[i]->Print();
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
            cout << "gamma rel sys for input " << i << endl;
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
                                                                                            fileNameIncGammaOutputWeighting, "8TeV", "GammaInc", kTRUE,
                                                                                            NULL, fileNameCorrelations );


    if (graphCombIncGammaTot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    
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
    //*******************************************************************************************************************************************
    //*********************************************** Combining PCM IncGamma spectra ************************************************************
    //*******************************************************************************************************************************************
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};

    Int_t offSetsGammaPCM[11]       = { 0,  0,  9,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetsGammaPCMSys[11]    = { 1,  0,  9,  0,  4,
                                        0,  0,  0,  0,  0,
                                        0};

    TGraph* graphWeightsPCMIncGamma[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsPCMIncGamma[i]                   = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNamePCMIncGammaOutputWeighting       = Form("%s/PCMIncGamma_WeightingMethod.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombPCMIncGammaStat      = NULL;
    TGraphAsymmErrors* graphCombPCMIncGammaSys       = NULL;
    TGraphAsymmErrors* graphCombPCMIncGammaTot       = CombinePtPointsSpectraFullCorrMat(      histoIncGammaStatErr,    graphIncGammaSysErr,
                                                                                            xPtLimitsGamma, maxNBinsGamma,
                                                                                            offSetsGammaPCM, offSetsGammaPCMSys,
                                                                                            graphCombPCMIncGammaStat, graphCombPCMIncGammaSys,
                                                                                            fileNamePCMIncGammaOutputWeighting, "8TeV", "GammaInc", kTRUE,
                                                                                            NULL, fileNameCorrelations );
    if (graphCombIncGammaTot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    while (graphCombPCMIncGammaStat->GetX()[0] < 0.3){
        graphCombPCMIncGammaStat->RemovePoint(0);
    }
    while (graphCombPCMIncGammaTot->GetX()[0] < 0.3){
        graphCombPCMIncGammaTot->RemovePoint(0);
    }
    while (graphCombPCMIncGammaSys->GetX()[0] < 0.3){
        graphCombPCMIncGammaSys->RemovePoint(0);
    }  

    // Reading weights from output file for plotting
    ifstream fileWeightsPCMIncGammaRead;
    fileWeightsPCMIncGammaRead.open(fileNamePCMIncGammaOutputWeighting,ios_base::in);
    cout << "reading" << fileNamePCMIncGammaOutputWeighting << endl;
    Double_t xValuesPCMIncGammaRead[50];
    Double_t weightsPCMIncGammaRead[11][50];
    Int_t availablePCMIncGammaMeas[11]      = { -1, -1, -1, -1, -1,
                                                -1, -1, -1, -1, -1,
                                                -1};
    Int_t nMeasSetPCMIncGamma               = 3;
    Int_t nPtBinsPCMIncGammaRead            = 0;
    while(!fileWeightsPCMIncGammaRead.eof() && nPtBinsPCMIncGammaRead < 50){
        TString garbage             = "";
        if (nPtBinsPCMIncGammaRead == 0){
            fileWeightsPCMIncGammaRead >> garbage ;
            for (Int_t i = 0; i < nMeasSetPCMIncGamma; i++){
                fileWeightsPCMIncGammaRead >> availablePCMIncGammaMeas[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availablePCMIncGammaMeas[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsPCMIncGammaRead >> xValuesPCMIncGammaRead[nPtBinsPCMIncGammaRead-1];
            for (Int_t i = 0; i < nMeasSetPCMIncGamma; i++){
                fileWeightsPCMIncGammaRead >> weightsPCMIncGammaRead[availablePCMIncGammaMeas[i]][nPtBinsPCMIncGammaRead-1] ;
            }
            cout << "read: "<<  nPtBinsPCMIncGammaRead << "\t"<< xValuesPCMIncGammaRead[nPtBinsPCMIncGammaRead-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetPCMIncGamma; i++){
                cout << weightsPCMIncGammaRead[availablePCMIncGammaMeas[i]][nPtBinsPCMIncGammaRead-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsPCMIncGammaRead++;
    }
    nPtBinsPCMIncGammaRead                  = nPtBinsPCMIncGammaRead-2 ;
    fileWeightsPCMIncGammaRead.close();

    for (Int_t i = 0; i < nMeasSetPCMIncGamma; i++){
        graphWeightsPCMIncGamma[availablePCMIncGammaMeas[i]]                        = new TGraph(nPtBinsPCMIncGammaRead,xValuesPCMIncGammaRead,weightsPCMIncGammaRead[availablePCMIncGammaMeas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsPCMIncGammaRead; n++){
            if (graphWeightsPCMIncGamma[availablePCMIncGammaMeas[i]]->GetY()[bin] == 0) graphWeightsPCMIncGamma[availablePCMIncGammaMeas[i]]->RemovePoint(bin);
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

    canvasWeights->SaveAs(Form("%s/IncGamma_Weights_8.%s",outputDir.Data(),suffix.Data()));
    canvasWeights->SaveAs(Form("%s/IncGamma_Weights_8.pdf",outputDir.Data()));

    histo2DIncGammaWeights->Draw("copy");

    TLegend* legendWeightsPCMIncGamma   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetPCMIncGamma+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
    for (Int_t i = 0; i < nMeasSetPCMIncGamma; i++){
        DrawGammaSetMarkerTGraph(graphWeightsPCMIncGamma[availablePCMIncGammaMeas[i]], markerStyleDet[availablePCMIncGammaMeas[i]], markerSizeDet[availablePCMIncGammaMeas[i]], colorDet[availablePCMIncGammaMeas[i]] , colorDet[availablePCMIncGammaMeas[i]]);
        graphWeightsPCMIncGamma[availablePCMIncGammaMeas[i]]->Draw("p,same,z");
        legendWeightsPCMIncGamma->AddEntry(graphWeightsPCMIncGamma[availablePCMIncGammaMeas[i]],nameMeasGlobalLabelGamma[availablePCMIncGammaMeas[i]],"p");
    }
    legendWeightsPCMIncGamma->Draw();

    labelWeightsEnergy->Draw();
    labelWeightsIncGamma->Draw();

    DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
    DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/IncGamma_Weights_8_PCM.%s",outputDir.Data(),suffix.Data()));
    canvasWeights->SaveAs(Form("%s/IncGamma_Weights_8_PCM.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelSysErr->cd();

    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,13.5);
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

    canvasRelSysErr->SaveAs(Form("%s/IncGamma_RelSysErr_8.%s",outputDir.Data(),suffix.Data()));
    canvasRelSysErr->SaveAs(Form("%s/IncGamma_RelSysErr_8.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    canvasRelStatErr->cd();

    histo2DRelStatErr->GetYaxis()->SetRangeUser(-0.2,13.5);
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

    canvasRelStatErr->SaveAs(Form("%s/IncGamma_RelStatErr_8.%s",outputDir.Data(),suffix.Data()));
    canvasRelStatErr->SaveAs(Form("%s/IncGamma_RelStatErr_8.pdf",outputDir.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TGraphAsymmErrors* graphCombIncGammaRelStat     = CalculateRelErrUpAsymmGraph( graphCombIncGammaStat, "relativeStatErrorDR");
    TGraphAsymmErrors* graphCombIncGammaRelSys      = CalculateRelErrUpAsymmGraph( graphCombIncGammaSys, "relativeSysErrorDR");
    TGraphAsymmErrors* graphCombIncGammaRelTot      = CalculateRelErrUpAsymmGraph( graphCombIncGammaTot, "relativeTotalErrorDR");

    canvasRelSysErr->cd();
        histo2DRelErr->GetYaxis()->SetRangeUser(0,12.5);
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

    canvasRelSysErr->SaveAs(Form("%s/IncGamma_Reldecomp_8.%s",outputDir.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    cout << "IncGamma - tot" << endl;
    graphCombIncGammaRelTot->Print();
    cout << "IncGamma - stat" << endl;
    graphCombIncGammaRelStat->Print();
    cout << "IncGamma - sys" << endl;
    graphCombIncGammaRelSys->Print();


    TGraphAsymmErrors* graphCombPCMIncGammaRelStat     = CalculateRelErrUpAsymmGraph( graphCombPCMIncGammaStat, "relativeStatErrorIncGammaPCM");
    TGraphAsymmErrors* graphCombPCMIncGammaRelSys      = CalculateRelErrUpAsymmGraph( graphCombPCMIncGammaSys, "relativeSysErrorIncGammaPCM");
    TGraphAsymmErrors* graphCombPCMIncGammaRelTot      = CalculateRelErrUpAsymmGraph( graphCombPCMIncGammaTot, "relativeTotalErrorIncGammaPCM");

    canvasRelSysErr->cd();
        histo2DRelErr->GetYaxis()->SetRangeUser(0,12.5);
        histo2DRelErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPCMIncGammaRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPCMIncGammaRelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombPCMIncGammaRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombPCMIncGammaRelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombPCMIncGammaRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombPCMIncGammaRelSys->SetLineStyle(7);
        graphCombPCMIncGammaRelSys->Draw("l,x0,same,e1");

        legendRelTotErr->Draw();
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrPCMIncGamma      = new TLatex(0.95,0.85,"#gamma_{inc} PCM");
        SetStyleTLatex( labelRelTotErrPCMIncGamma, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRelTotErrPCMIncGamma->Draw();

    canvasRelSysErr->SaveAs(Form("%s/IncGamma_Reldecomp_8_PCM.%s",outputDir.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    cout << "PCMIncGamma - tot" << endl;
    graphCombPCMIncGammaRelTot->Print();
    cout << "PCMIncGamma - stat" << endl;
    graphCombPCMIncGammaRelStat->Print();
    cout << "PCMIncGamma - sys" << endl;
    graphCombPCMIncGammaRelSys->Print();

    //*******************************************************************************************************************************************
    //************************************************* Fitting gamma spectrum ******************************************************************
    //*******************************************************************************************************************************************
    //NOTE: dirty hack to use pcm histo for fitting!
    // TF1* fitHagGammaComb            = FitObject("h","fitHagGammaComb","Gamma",graphCombIncGammaTot,graphCombIncGammaTot->GetX()[0],
                                                    // graphCombIncGammaTot->GetX()[graphCombIncGammaTot->GetN()],NULL,"QNRME+");
    Double_t paramGraphoHag[5]                      = {82,-0.22,-0.01,0.57,6.28};
    TF1* fitHagGammaComb            = FitObject("oHag","fitHagGammaComb","Gamma",graphCombIncGammaTot,graphCombIncGammaTot->GetX()[2],
                                                     graphCombIncGammaTot->GetX()[graphCombIncGammaTot->GetN()],paramGraphoHag,"QNRME+");
    Double_t paramTCM[5] = {graphCombIncGammaTot->GetY()[1],0.1,graphCombIncGammaTot->GetY()[4],0.6,3};
    TF1* fitTCMGammaComb            = FitObject("tcm", "fitTCMGammaComb","Gamma", graphCombIncGammaTot, graphCombIncGammaTot->GetX()[0], graphCombIncGammaTot->GetX()[graphCombIncGammaTot->GetN()], paramTCM,"QNRME+");
    TString forOutput               = WriteParameterToFile(fitTCMGammaComb);
    cout << forOutput.Data() << endl;

    TF1* fitTsallisGammaComb        = FitObject("l","fitTsallisGammaComb","Gamma",graphCombIncGammaTot,graphCombIncGammaTot->GetX()[0],
                                                 graphCombIncGammaTot->GetX()[graphCombIncGammaTot->GetN()],NULL,"QNRME+");
    forOutput                       = WriteParameterToFile(fitTsallisGammaComb);
    cout << forOutput.Data() << endl;


    // Cloning spectra

    TGraphAsymmErrors* graphIndGammaIncStat[11]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndGammaIncSys[11]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    for (Int_t i = 0; i< 11; i++){
        if (statErrorGraphCollectionIncGamma[i]){
            graphIndGammaIncStat[i]                      = (TGraphAsymmErrors*)statErrorGraphCollectionIncGamma[i]->Clone(Form("GammaStat%s",nameMeasGlobalLabel[i].Data()));
        }
        if (graphIncGammaSysErr[i]){
            graphIndGammaIncSys[i]                       = (TGraphAsymmErrors*)graphIncGammaSysErr[i]->Clone(Form("GammaSys%s",nameMeasGlobalLabel[i].Data()));
        }
    }

    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************

    TGraphAsymmErrors* graphRatioGammaIndCombFitStat[11] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioGammaIndCombFitSys[11]  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    TGraphAsymmErrors* graphRatioGammaCombCombFitTot     = (TGraphAsymmErrors*)graphCombIncGammaTot->Clone();
    graphRatioGammaCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitTot, fitTCMGammaComb);
    TGraphAsymmErrors* graphRatioGammaCombCombFitStat    = (TGraphAsymmErrors*)graphCombIncGammaStat->Clone();
    graphRatioGammaCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitStat, fitTCMGammaComb);
    TGraphAsymmErrors* graphRatioGammaCombCombFitSys     = (TGraphAsymmErrors*)graphCombIncGammaSys->Clone();
    graphRatioGammaCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioGammaCombCombFitSys, fitTCMGammaComb);
    
    TGraphAsymmErrors* graphRatioGammaPCMCombCombFitTot     = (TGraphAsymmErrors*)graphCombPCMIncGammaTot->Clone();
    graphRatioGammaPCMCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioGammaPCMCombCombFitTot, fitTCMGammaComb);
    TGraphAsymmErrors* graphRatioGammaPCMCombCombFitStat    = (TGraphAsymmErrors*)graphCombPCMIncGammaStat->Clone();
    graphRatioGammaPCMCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioGammaPCMCombCombFitStat, fitTCMGammaComb);
    TGraphAsymmErrors* graphRatioGammaPCMCombCombFitSys     = (TGraphAsymmErrors*)graphCombPCMIncGammaSys->Clone();
    graphRatioGammaPCMCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioGammaPCMCombCombFitSys, fitTCMGammaComb);

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
    histo2DGammaRatioToCombFit               = new TH2F("histo2DGammaRatioToCombFit","histo2DGammaRatioToCombFit",1000,0.23, 25.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DGammaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPPb, textsizeLabelsPPb,
                              0.85*textsizeLabelsPPb,textsizeLabelsPPb, 0.9, 0.65, 510, 505);
    histo2DGammaRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DGammaRatioToCombFit->GetXaxis()->SetNoExponent();
    //  histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(0.75,1.46);
    histo2DGammaRatioToCombFit->Draw("copy");

    ProduceGraphAsymmWithoutXErrors(graphRatioGammaCombCombFitStat);

    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitSys, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV, widthLinesBoxes, kTRUE);
    graphRatioGammaCombCombFitSys->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphRatioGammaCombCombFitStat, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV);
    graphRatioGammaCombCombFitStat->Draw("p,same,z");

    DrawGammaLines(0.23, 25. , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0.23, 25. , 1.05, 1.05,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.95, 0.95,0.1, kGray, 7);

    TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.91, Form("%s, %s",textALICE.Data(),collisionSystempp8TeV.Data()));
    SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergy->Draw();
    TLatex *labelRatioToFitGamma      = new TLatex(0.12, 0.92, "#gamma_{inc}");
    SetStyleTLatex( labelRatioToFitGamma, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
    labelRatioToFitGamma->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_8TeV.%s",outputDir.Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfCombToCombFit_8TeV.pdf",outputDir.Data()));

    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(0.59,1.47);
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
    labelRatioToFitGamma->Draw();
    histo2DGammaRatioToCombFit->Draw("same,axis");


    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    Double_t rowsLegendGammaRatio[5]          = {0.31, 0.26, 0.21, 0.16, 0.16};
    Double_t rowsLegendGammaRatioAbs[5]       = {0.74, 0.69, 0.64, 0.59 };
    Double_t columnsLegendGammaRatio[7]       = {0.14, 0.28, 0.38, 0.48, 0.7, 0.8};
    Double_t columnsLegendGammaRatioAbs[7]    = {0.215, 0.74, 1.24, 2, 6.6, 9.6};
    Double_t lengthBox                          = 0.2;
    Double_t heightBox                          = 0.03/2;
    //****************** first Column **************************************************
    TLatex *textPCMRatioGamma                 = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[1],nameMeasGlobalLabelGamma[0]);
    SetStyleTLatex( textPCMRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textPCMRatioGamma->Draw();
    TLatex *textEMCALRatioGamma               = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[2],nameMeasGlobalLabelGamma[4]);
    SetStyleTLatex( textEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textEMCALRatioGamma->Draw();
    // TLatex *textPCMEMCALRatioGamma            = new TLatex(columnsLegendGammaRatio[3],rowsLegendGammaRatio[1],nameMeasGlobalLabel[4]);
    // SetStyleTLatex( textPCMEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    // textPCMEMCALRatioGamma->Draw();
    TLatex *textPCMEMCALRatioGamma            = new TLatex(columnsLegendGammaRatio[0],rowsLegendGammaRatio[3],nameMeasGlobalLabelGamma[2]);
    SetStyleTLatex( textPCMEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textPCMEMCALRatioGamma->Draw();

    //****************** second Column *************************************************
    TLatex *textStatRatioGamma                = new TLatex(columnsLegendGammaRatio[1],rowsLegendGammaRatio[0] ,"stat");
    SetStyleTLatex( textStatRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textStatRatioGamma->Draw();
    TLatex *textSysRatioGamma                 = new TLatex(columnsLegendGammaRatio[2] ,rowsLegendGammaRatio[0],"syst");
    SetStyleTLatex( textSysRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textSysRatioGamma->Draw();

    TMarker* markerPCMGammaRatio           = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[0],columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[1],1);
    markerPCMGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[0]);
    TMarker* markerEMCALGammaRatio         = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[4], columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[2],1);
    markerEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[1]);
    TMarker* markerPCMEMCALGammaRatio      = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[2], columnsLegendGammaRatio[1] ,rowsLegendGammaRatio[3],1);
    markerPCMEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbs[1] ,rowsLegendGammaRatioAbs[2]);

    TBox* boxPCMGammaRatio                 = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[0], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[0]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[0]+ heightBox);
    boxPCMGammaRatio->Draw("l");
    TBox* boxEMCALGammaRatio               = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[4], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[1]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[1]+ heightBox);
    boxEMCALGammaRatio->Draw("l");
    TBox* boxPCMEMCALGammaRatio            = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[2], columnsLegendGammaRatioAbs[2]-0.8*lengthBox , rowsLegendGammaRatioAbs[2]- heightBox, columnsLegendGammaRatioAbs[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbs[2]+ heightBox);
    boxPCMEMCALGammaRatio->Draw("l");

    // canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit_0.pdf",outputDir.Data()));

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit_8.%s",outputDir.Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit_8.pdf",outputDir.Data()));
    
    histo2DGammaRatioToCombFit->GetYaxis()->SetRangeUser(0.73,1.27);
    histo2DGammaRatioToCombFit->Draw("copy");
    if (graphRatioGammaPCMCombCombFitSys){
        DrawGammaSetMarkerTGraphAsym(graphRatioGammaPCMCombCombFitSys, markerStyleDet[0] ,markerSizeDet[0], colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        graphRatioGammaPCMCombCombFitSys->Draw("E2same");
    }
    if (graphRatioGammaPCMCombCombFitStat){
        ProduceGraphAsymmWithoutXErrors(graphRatioGammaPCMCombCombFitStat);
        DrawGammaSetMarkerTGraphAsym(graphRatioGammaPCMCombCombFitStat, markerStyleDet[0] ,markerSizeDet[0], colorDet[0], colorDet[0]);
        graphRatioGammaPCMCombCombFitStat->Draw("p,same,z");
    }
    if (graphRatioGammaIndCombFitSys[2]){
        DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitSys[2], markerStyleDet[2] ,markerSizeDet[2], colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        graphRatioGammaIndCombFitSys[2]->Draw("E2same");
    }
    if (graphRatioGammaIndCombFitStat[2]){
        ProduceGraphAsymmWithoutXErrors(graphRatioGammaIndCombFitStat[2]);
        DrawGammaSetMarkerTGraphAsym(graphRatioGammaIndCombFitStat[2], markerStyleDet[2] ,markerSizeDet[2], colorDet[2], colorDet[2]);
        graphRatioGammaIndCombFitStat[2]->Draw("p,same,z");
    }

    DrawGammaLines(0.23, 25. , 1., 1.,0.5, kGray+2);
    DrawGammaLines(0.23, 25. , 1.05, 1.05,0.5, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.95, 0.95,0.5, kGray, 7);
    DrawGammaLines(0.23, 25. , 1.1, 1.1,0.5, kGray, 9);
    DrawGammaLines(0.23, 25. , 0.9, 0.9,0.5, kGray, 9);

    labelRatioToFitEnergy->Draw();
    labelRatioToFitGamma->Draw();
    histo2DGammaRatioToCombFit->Draw("same,axis");


    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    Double_t rowsLegendGammaRatioPCMcomb[5]          = {0.31, 0.27, 0.21, 0.16, 0.16};
    Double_t rowsLegendGammaRatioAbsPCMcomb[5]       = {0.80, 0.77, 0.64, 0.59 };
    Double_t columnsLegendGammaRatioPCMcomb[7]       = {0.14, 0.28, 0.38, 0.48, 0.7, 0.8};
    Double_t columnsLegendGammaRatioAbsPCMcomb[7]    = {0.215, 0.74, 1.24, 2, 6.6, 9.6};
    Double_t lengthBoxPCMcomb                          = 0.2;
    Double_t heightBoxPCMcomb                          = 0.016/2;
    //****************** first Column **************************************************
    textPCMRatioGamma                 = new TLatex(columnsLegendGammaRatioPCMcomb[0],rowsLegendGammaRatioPCMcomb[2],nameMeasGlobalLabelGamma[0]);
    SetStyleTLatex( textPCMRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textPCMRatioGamma->Draw();
    textEMCALRatioGamma               = new TLatex(columnsLegendGammaRatioPCMcomb[0],rowsLegendGammaRatioPCMcomb[3],nameMeasGlobalLabelGamma[2]);
    SetStyleTLatex( textEMCALRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textEMCALRatioGamma->Draw();

    //****************** second Column *************************************************
    textStatRatioGamma                = new TLatex(columnsLegendGammaRatioPCMcomb[1],rowsLegendGammaRatioPCMcomb[1] ,"stat");
    SetStyleTLatex( textStatRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textStatRatioGamma->Draw();
    textSysRatioGamma                 = new TLatex(columnsLegendGammaRatioPCMcomb[2] ,rowsLegendGammaRatioPCMcomb[1],"syst");
    SetStyleTLatex( textSysRatioGamma, textSizeLabelsPixel,4, 1, 43);
    textSysRatioGamma->Draw();

    markerPCMGammaRatio           = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[0],columnsLegendGammaRatioPCMcomb[1] ,rowsLegendGammaRatioPCMcomb[1],1);
    markerPCMGammaRatio->DrawMarker(columnsLegendGammaRatioAbsPCMcomb[1] ,rowsLegendGammaRatioAbsPCMcomb[0]);
    markerEMCALGammaRatio         = CreateMarkerFromGraph(graphRatioGammaIndCombFitSys[2], columnsLegendGammaRatioPCMcomb[1] ,rowsLegendGammaRatioPCMcomb[2],1);
    markerEMCALGammaRatio->DrawMarker(columnsLegendGammaRatioAbsPCMcomb[1] ,rowsLegendGammaRatioAbsPCMcomb[1]);

    boxPCMGammaRatio                 = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[0], columnsLegendGammaRatioAbsPCMcomb[2]-0.8*lengthBox , rowsLegendGammaRatioAbsPCMcomb[0]- heightBox, columnsLegendGammaRatioAbsPCMcomb[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbsPCMcomb[0]+ heightBoxPCMcomb);
    boxPCMGammaRatio->Draw("l");
    boxEMCALGammaRatio               = CreateBoxFromGraph(graphRatioGammaIndCombFitSys[2], columnsLegendGammaRatioAbsPCMcomb[2]-0.8*lengthBox , rowsLegendGammaRatioAbsPCMcomb[1]- heightBox, columnsLegendGammaRatioAbsPCMcomb[2]+ 1.1*lengthBox, rowsLegendGammaRatioAbsPCMcomb[1]+ heightBoxPCMcomb);
    boxEMCALGammaRatio->Draw("l");

    // canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit_0.pdf",outputDir.Data()));

    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit_combPCM_8.%s",outputDir.Data(),suffix.Data()));
    canvasRatioToCombFit->SaveAs(Form("%s/Gamma_RatioOfIndividualMeasToCombFit_combPCM_8.pdf",outputDir.Data()));

    //*******************************************************************************************************
    //************************** Calculating combined direct photon spectrum ********************************
    //*******************************************************************************************************
    Double_t xArrayCombDR[graphCombDRStat->GetN()+1];
    xArrayCombDR[0] = graphCombDRStat->GetX()[0] - graphCombDRStat->GetEXhigh()[0];
    for (Int_t i = 1; i<graphCombDRStat->GetN()+1;i++){
        xArrayCombDR[i] = graphCombDRStat->GetX()[i-1] + graphCombDRStat->GetEXhigh()[i-1];
    }
    Double_t xArrayCombDRNonFit[graphCombDRNonFitStat->GetN()+1];
    xArrayCombDRNonFit[0] = graphCombDRNonFitStat->GetX()[0] - graphCombDRNonFitStat->GetEXhigh()[0];
    for (Int_t i = 1; i<graphCombDRNonFitStat->GetN()+1;i++){
        xArrayCombDRNonFit[i] = graphCombDRNonFitStat->GetX()[i-1] + graphCombDRNonFitStat->GetEXhigh()[i-1];
    }
    //_______________________ copy inclusive photon spectra _____________________
    TH1D *histoCombDirGammaSpectrumErrSum                   = new TH1D("histoCombDirGammaSpectrumErrSum","",graphCombDRStat->GetN(),xArrayCombDR);
    TH1D *histoCombDirGammaSpectrumErrSys                   = new TH1D("histoCombDirGammaSpectrumErrSys","",graphCombDRStat->GetN(),xArrayCombDR);
    TH1D *histoCombDirGammaSpectrumErrStat                  = new TH1D("histoCombDirGammaSpectrumErrStat","",graphCombDRStat->GetN(),xArrayCombDR);
    TH1D *histoCombDirGammaNonFitSpectrumErrSum                   = new TH1D("histoCombDirGammaSpectrumErrSum","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);
    TH1D *histoCombDirGammaNonFitSpectrumErrSys                   = new TH1D("histoCombDirGammaSpectrumErrSys","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);
    TH1D *histoCombDirGammaNonFitSpectrumErrStat                  = new TH1D("histoCombDirGammaSpectrumErrStat","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);

    //_______________________ get arrays of double ratio errors __________________
    Double_t *SystErrorsCombDR                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *sumErrorsCombDR                               = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *StatErrorsCombDR                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *xErrorsDR                                     = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *SystErrorsCombDRNonFit                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *sumErrorsCombDRNonFit                               = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *StatErrorsCombDRNonFit                              = new Double_t[graphCombIncGammaStat->GetN()];
    Double_t *xErrorsDRNonFit                                     = new Double_t[graphCombIncGammaStat->GetN()];
    for (Int_t i = 0; i< graphCombDRStat->GetN(); i++){
        SystErrorsCombDR[i]                                 = graphCombDRSys->GetEYhigh()[i]/graphCombDRSys->GetY()[i] *100;
        StatErrorsCombDR[i]                                 = graphCombDRStat->GetEYhigh()[i]/graphCombDRStat->GetY()[i] *100;
        sumErrorsCombDR[i]                                  = graphCombDRTot->GetEYhigh()[i]/graphCombDRTot->GetY()[i] *100;
    }
    for (Int_t i = 0; i< graphCombDRNonFitStat->GetN(); i++){
        SystErrorsCombDRNonFit[i]                                 = graphCombDRNonFitSys->GetEYhigh()[i]/graphCombDRNonFitSys->GetY()[i] *100;
        StatErrorsCombDRNonFit[i]                                 = graphCombDRNonFitStat->GetEYhigh()[i]/graphCombDRNonFitStat->GetY()[i] *100;
        sumErrorsCombDRNonFit[i]                                  = graphCombDRNonFitTot->GetEYhigh()[i]/graphCombDRNonFitTot->GetY()[i] *100;
    }
    xErrorsDR                                               = graphCombDRStat->GetX();
    xErrorsDRNonFit                                               = graphCombDRStat->GetX();

    cout << __LINE__ << endl;
    //graphCombDRTot->Print();

    //_______________________ copy inclusive photon spectra _____________________
    TH1D* histoCombErrorsForDRSum                           = new TH1D("histoCombErrorsForDRSum","",graphCombDRStat->GetN(),xArrayCombDR);
    TH1D* histoCombErrorsForDRStat                          = new TH1D("histoCombErrorsForDRStat","",graphCombDRStat->GetN(),xArrayCombDR);
    TH1D* histoCombErrorsForDRSys                           = new TH1D("histoCombErrorsForDRSys","",graphCombDRStat->GetN(),xArrayCombDR);
    TH1D* histoCombErrorsForDRNonFitSum                           = new TH1D("histoCombErrorsForDRNonFitSum","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);
    TH1D* histoCombErrorsForDRNonFitStat                          = new TH1D("histoCombErrorsForDRNonFitStat","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);
    TH1D* histoCombErrorsForDRNonFitSys                           = new TH1D("histoCombErrorsForDRNonFitSys","",graphCombDRNonFitStat->GetN(),xArrayCombDRNonFit);

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
    for(Int_t i = 1; i<graphCombDRNonFitStat->GetN()+1;i++){
        //cout<< i << "\t"<<xErrorsDR[i-1]<<"  "<<histoCombErrorsForDRSum->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum->GetBinWidth(i) <<endl;
        Double_t binErrorSummed                                 = sumErrorsCombDRNonFit[i-1];
        Double_t binErrorSyst                                   = SystErrorsCombDRNonFit[i-1];
        Double_t binErrorStat                                   = StatErrorsCombDRNonFit[i-1];
        Double_t DR                                             = graphCombDRNonFitStat->GetY()[i-1];

        //cout << DR << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
        histoCombErrorsForDRNonFitSum->SetBinContent(i,DR);
        histoCombErrorsForDRNonFitSys->SetBinContent(i,DR);
        histoCombErrorsForDRNonFitStat->SetBinContent(i,DR);
        histoCombErrorsForDRNonFitSum->SetBinError(i,(binErrorSummed/100)*DR);
        histoCombErrorsForDRNonFitSys->SetBinError(i,(binErrorSyst/100)*DR);
        histoCombErrorsForDRNonFitStat->SetBinError(i,(binErrorStat/100)*DR);
    }

    for(Int_t i = 1; i<histoCombErrorsForDRSum->GetNbinsX()+1;i++){
        histoCombDirGammaSpectrumErrSum->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrSys->SetBinContent(i+1,-1);
        histoCombDirGammaSpectrumErrStat->SetBinContent(i+1,-1);

        histoCombDirGammaSpectrumErrSum->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrSys->SetBinError(i+1,0);
        histoCombDirGammaSpectrumErrStat->SetBinError(i+1,0);
    }
    for(Int_t i = 1; i<histoCombErrorsForDRNonFitSum->GetNbinsX()+1;i++){
        histoCombDirGammaNonFitSpectrumErrSum->SetBinContent(i+1,-1);
        histoCombDirGammaNonFitSpectrumErrSys->SetBinContent(i+1,-1);
        histoCombDirGammaNonFitSpectrumErrStat->SetBinContent(i+1,-1);

        histoCombDirGammaNonFitSpectrumErrSum->SetBinError(i+1,0);
        histoCombDirGammaNonFitSpectrumErrSys->SetBinError(i+1,0);
        histoCombDirGammaNonFitSpectrumErrStat->SetBinError(i+1,0);
    }

    // get the binning of the direct photons from the DR
    TH1D *histoCombDirGammaSpecSysErr                        = new TH1D(*histoCombErrorsForDRSys);
    TH1D *histoCombDirGammaSpecStatErr                       = new TH1D(*histoCombErrorsForDRStat);
    TH1D *histoCombDirGammaSpecSumErr                        = new TH1D(*histoCombErrorsForDRSum);
    TH1D *histoCombDirGammaSpecNonFitSysErr                        = new TH1D(*histoCombErrorsForDRNonFitSys);
    TH1D *histoCombDirGammaSpecNonFitStatErr                       = new TH1D(*histoCombErrorsForDRNonFitStat);
    TH1D *histoCombDirGammaSpecNonFitSumErr                        = new TH1D(*histoCombErrorsForDRNonFitSum);


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
        Double_t errDR                  = content - error;
        histoCombDirGammaSpecNonFitSysErr->SetBinError(i, error);
        histoCombDirGammaSpecNonFitSysErr->SetBinContent(i, content);
        histoCombDirGammaNonFitSpectrumErrSys->SetBinContent(i, errDR);

        // calculating Stat graphs
        errRgamma                       = histoCombErrorsForDRNonFitStat->GetBinError(i);
        errNIncGam                      = graphCombIncGammaStat->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecNonFitStatErr->SetBinError(i, error);
        histoCombDirGammaSpecNonFitStatErr->SetBinContent(i, content);
        histoCombDirGammaNonFitSpectrumErrStat->SetBinContent(i, errDR);

        // calculating summed error graphs
        errRgamma                       = histoCombErrorsForDRNonFitSum->GetBinError(i);
        errNIncGam                      = graphCombIncGammaTot->GetEYhigh()[i-1];
        q1                              = 1 - 1/ Rgamma;
        q1Error                         = errRgamma/(Rgamma*Rgamma);
        content                         = nIncGamma * ( 1 - 1/ Rgamma);
        error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR                           = content - error;
        histoCombDirGammaSpecNonFitSumErr->SetBinError(i, error);
        histoCombDirGammaSpecNonFitSumErr->SetBinContent(i, content);
        histoCombDirGammaNonFitSpectrumErrSum->SetBinContent(i, errDR);
    }

    // purely calculating points based on all Systematic errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSystErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSys,histoCombDirGammaSpecStatErr,0,0.5,0,confidenceLevelnSigma);
    if(graphCombDirGammaSpectrumSystErr)graphCombDirGammaSpectrumSystErr->SetName("graphCombDirGammaSpectrumSystErr");
    if(graphCombDirGammaSpectrumSystErr) cout << "sys has been found" << endl;
    if(graphCombDirGammaSpectrumSystErr)graphCombDirGammaSpectrumSystErr->Print();

    // purely calculating points based on Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumStatErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrStat,histoCombDirGammaSpecStatErr,0,0.5,0,confidenceLevelnSigma);
    if(graphCombDirGammaSpectrumStatErr)graphCombDirGammaSpectrumStatErr->SetName("graphCombDirGammaSpectrumStatErr");
    if(graphCombDirGammaSpectrumStatErr) cout << "stat has been found" << endl;
    if(graphCombDirGammaSpectrumStatErr)graphCombDirGammaSpectrumStatErr->Print();
    // purely calculating points based on all Systematic + Statistical errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum,histoCombDirGammaSpecStatErr,0,0.5,0,confidenceLevelnSigma);
    if(graphCombDirGammaSpectrumSumErr)graphCombDirGammaSpectrumSumErr->SetName("graphCombDirGammaSpectrumSumErr");
    if(graphCombDirGammaSpectrumSumErr) cout << "tot has been found" << endl;
    if(graphCombDirGammaSpectrumSumErr)graphCombDirGammaSpectrumSumErr->Print();
    // calculate points above confidence level summed errors
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErrConfi = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum,histoCombDirGammaSpecStatErr,2,0.5,0,confidenceLevelnSigma);
    if(graphCombDirGammaSpectrumSumErrConfi)graphCombDirGammaSpectrumSumErrConfi->SetName("graphCombDirGammaSpectrumSumErrConfi");
    if(graphCombDirGammaSpectrumSumErrConfi) cout << "confi has been found" << endl;
    if(graphCombDirGammaSpectrumSumErrConfi)graphCombDirGammaSpectrumSumErrConfi->Print();
    // calculate arrows for points with 0, error summed
    TGraphAsymmErrors *graphCombDirGammaSpectrumSumErrAr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum,histoCombDirGammaSpecStatErr,5,0.5,0,confidenceLevelnSigma);
    if(graphCombDirGammaSpectrumSumErrAr)graphCombDirGammaSpectrumSumErrAr->SetName("graphCombDirGammaSpectrumSumErrAr");
    if(graphCombDirGammaSpectrumSumErrAr) cout << "Ar has been found" << endl;
    if(graphCombDirGammaSpectrumSumErrAr)graphCombDirGammaSpectrumSumErrAr->Print();



    // purely calculating points based on all Systematic errors
    TGraphAsymmErrors *graphCombDirGammaNonFitSpectrumSystErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaNonFitSpectrumErrSys,histoCombDirGammaSpecNonFitStatErr,0,0.5,0,confidenceLevelnSigma);
    if(graphCombDirGammaNonFitSpectrumSystErr)graphCombDirGammaNonFitSpectrumSystErr->SetName("graphCombDirGammaNonFitSpectrumSystErr");
    if(graphCombDirGammaNonFitSpectrumSystErr) cout << "sys has been found" << endl;
    if(graphCombDirGammaNonFitSpectrumSystErr)graphCombDirGammaNonFitSpectrumSystErr->Print();

    // purely calculating points based on Statistical errors
    TGraphAsymmErrors *graphCombDirGammaNonFitSpectrumStatErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaNonFitSpectrumErrStat,histoCombDirGammaSpecNonFitStatErr,0,0.5,0,confidenceLevelnSigma);
    if(graphCombDirGammaNonFitSpectrumStatErr)graphCombDirGammaNonFitSpectrumStatErr->SetName("graphCombDirGammaNonFitSpectrumStatErr");
    if(graphCombDirGammaNonFitSpectrumStatErr) cout << "stat has been found" << endl;
    if(graphCombDirGammaNonFitSpectrumStatErr)graphCombDirGammaNonFitSpectrumStatErr->Print();
    // purely calculating points based on all Systematic + Statistical errors
    TGraphAsymmErrors *graphCombDirGammaNonFitSpectrumSumErr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaNonFitSpectrumErrSum,histoCombDirGammaSpecNonFitStatErr,0,0.5,0,confidenceLevelnSigma);
    if(graphCombDirGammaNonFitSpectrumSumErr)graphCombDirGammaNonFitSpectrumSumErr->SetName("graphCombDirGammaNonFitSpectrumSumErr");
    if(graphCombDirGammaNonFitSpectrumSumErr) cout << "tot has been found" << endl;
    if(graphCombDirGammaNonFitSpectrumSumErr)graphCombDirGammaNonFitSpectrumSumErr->Print();
    // calculate points above confidence level summed errors
    TGraphAsymmErrors *graphCombDirGammaNonFitSpectrumSumErrConfi = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaNonFitSpectrumErrSum,histoCombDirGammaSpecNonFitStatErr,2,0.5,0,confidenceLevelnSigma);
    if(graphCombDirGammaNonFitSpectrumSumErrConfi)graphCombDirGammaNonFitSpectrumSumErrConfi->SetName("graphCombDirGammaNonFitSpectrumSumErrConfi");
    if(graphCombDirGammaNonFitSpectrumSumErrConfi) cout << "confi has been found" << endl;
    if(graphCombDirGammaNonFitSpectrumSumErrConfi)graphCombDirGammaNonFitSpectrumSumErrConfi->Print();
    // calculate arrows for points with 0, error summed
    TGraphAsymmErrors *graphCombDirGammaNonFitSpectrumSumErrAr = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaNonFitSpectrumErrSum,histoCombDirGammaSpecNonFitStatErr,5,0.5,0,confidenceLevelnSigma);
    if(graphCombDirGammaNonFitSpectrumSumErrAr)graphCombDirGammaNonFitSpectrumSumErrAr->SetName("graphCombDirGammaNonFitSpectrumSumErrAr");
    if(graphCombDirGammaNonFitSpectrumSumErrAr) cout << "Ar has been found" << endl;
    if(graphCombDirGammaNonFitSpectrumSumErrAr)graphCombDirGammaNonFitSpectrumSumErrAr->Print();



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

      TLegend* legendDRSingle = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*3,0.31,0.953, textSizeSinglePad, 1, "", 42, 0.3);
        legendDRSingle->SetTextAlign(11);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 0; i < 11; i++){
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


        TLatex *labelDRSingle = new TLatex(0.95,0.92,Form("%s, %s",textALICE.Data(),collisionSystempp8TeV.Data()));
        SetStyleTLatex( labelDRSingle, textSizeSinglePad,4, 1, 42, kTRUE, 31);
        labelDRSingle->Draw();
        // TLatex *labelALICEDRSingle = new TLatex(0.95,0.87,textALICE.Data());
        // SetStyleTLatex( labelALICEDRSingle, textSizeSinglePad,4, 1, 42, kTRUE, 31);
        // labelALICEDRSingle->Draw();


        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pp8TeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_pp8TeV.pdf", outputDir.Data()));

    hist2DDRDummySingle->DrawCopy();

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 10; i > -1; i--){
            if (graphDRSysErr[i]){
                DrawGammaSetMarkerTGraphAsym(graphDRSysErr[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphDRSysErr[i]->Draw("E2same");
            }
        
            if (histoDRStatErr[i]){
                DrawGammaSetMarker(histoDRStatErr[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoDRStatErr[i]->Draw("p,same,e0,X0");
            }
        }

        if (histoDRStatErr[2]) histoDRStatErr[2]->Draw("p,same,e0,X0");
        if (histoDRStatErr[0]) histoDRStatErr[0]->Draw("p,same,e0,X0");
        if (histoDRStatErr[4]) histoDRStatErr[4]->Draw("p,same,e0,X0");
        legendDRSingle->Draw();

        labelDRSingle->Draw();;
        // labelALICEDRSingle->Draw();


        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_NonFit_pp8TeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_NonFit_pp8TeV.pdf", outputDir.Data()));
    
    TF1* pol0DRfitslow[11];
    TF1* pol0DRfitshigh[11];
    for (Int_t i = 0; i < 11; i++){
      pol0DRfitslow[i]            = NULL;
      pol0DRfitshigh[i]           = NULL;
      if (histoDRStatErr[i]){
        pol0DRfitslow[i]          = new TF1(Form("fitPol0low_%d",i),"[0]",1.,3.);
        pol0DRfitshigh[i]         = new TF1(Form("fitPol0high_%d",i),"[0]",1.,6.);
        histoDRStatErr[i]->Fit(pol0DRfitslow[i],"NRME+","",1.,3.);
        histoDRStatErr[i]->Fit(pol0DRfitshigh[i],"NRME+","",1.,6.);
      }
    }

    hist2DDRDummySingle->DrawCopy();
    legendDRSingle->Clear();
    legendDRSingle = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*4,0.31,0.953, textSizeSinglePad, 1, "Pol0 1-3 GeV/#it{c} fit:", 42, 0.3);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 0; i < 11; i++){
            if (graphDRSysErr[i]){
                DrawGammaSetMarkerTGraphAsym(graphDRSysErr[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphDRSysErr[i]->Draw("E2same");
                if (graphDRSysErr[i])legendDRSingle->AddEntry(graphDRSysErr[i],Form("%s: %1.3f #pm %1.3f",nameMeasGlobalLabel[i].Data(),pol0DRfitslow[i]->GetParameter(0),pol0DRfitslow[i]->GetParError(0)),"pf");
    
            }
            if (histoDRStatErr[i]){
                DrawGammaSetMarker(histoDRStatErr[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoDRStatErr[i]->Draw("p,same,e0,X0");
            }
        }
    
        if (histoDRStatErr[2]) histoDRStatErr[2]->Draw("p,same,e0,X0");
        if (histoDRStatErr[0]) histoDRStatErr[0]->Draw("p,same,e0,X0");
        if (histoDRStatErr[4]) histoDRStatErr[4]->Draw("p,same,e0,X0");
        for (Int_t i = 0; i < 11; i++){
            if(pol0DRfitslow[i]){
              DrawGammaSetMarkerTF1( pol0DRfitslow[i], 0, 2, colorDet[i]);
              pol0DRfitslow[i]->Draw("same");
            }
        }
        legendDRSingle->Draw();
        labelDRSingle->Draw();
    
        hist2DDRDummySingle->Draw("same,axis");
    
    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_NonFit_pp8TeV_lowpol0fit.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_NonFit_pp8TeV_lowpol0fit.pdf", outputDir.Data()));
    
    hist2DDRDummySingle->DrawCopy();
    legendDRSingle->Clear();
    legendDRSingle = GetAndSetLegend2(0.12,0.953-textSizeSinglePad*4,0.31,0.953, textSizeSinglePad, 1, "Pol0 1-6 GeV/#it{c} fit:", 42, 0.3);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        for (Int_t i = 0; i < 11; i++){
            if (graphDRSysErr[i]){
                DrawGammaSetMarkerTGraphAsym(graphDRSysErr[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
                graphDRSysErr[i]->Draw("E2same");
                if (graphDRSysErr[i])legendDRSingle->AddEntry(graphDRSysErr[i],Form("%s: %1.3f #pm %1.3f",nameMeasGlobalLabel[i].Data(),pol0DRfitshigh[i]->GetParameter(0),pol0DRfitshigh[i]->GetParError(0)),"pf");
    
            }
            if (histoDRStatErr[i]){
                DrawGammaSetMarker(histoDRStatErr[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
                histoDRStatErr[i]->Draw("p,same,e0,X0");
            }
        }
    
        if (histoDRStatErr[2]) histoDRStatErr[2]->Draw("p,same,e0,X0");
        if (histoDRStatErr[0]) histoDRStatErr[0]->Draw("p,same,e0,X0");
        if (histoDRStatErr[4]) histoDRStatErr[4]->Draw("p,same,e0,X0");
        for (Int_t i = 0; i < 11; i++){
            if(pol0DRfitshigh[i]){
              DrawGammaSetMarkerTF1( pol0DRfitshigh[i], 0, 2, colorDet[i]);
              pol0DRfitshigh[i]->Draw("same");
            }
        }
        legendDRSingle->Draw();
        labelDRSingle->Draw();
    
        hist2DDRDummySingle->Draw("same,axis");
    
    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_NonFit_pp8TeV_highpol0fit.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_IndMeasurements_NonFit_pp8TeV_highpol0fit.pdf", outputDir.Data()));

    hist2DDRDummySingle->DrawCopy();
        TGraphAsymmErrors* graphCombDRStatPlot    = (TGraphAsymmErrors*)graphCombDRStat->Clone("graphCombDRStatPlot");
        ProduceGraphAsymmWithoutXErrors(graphCombDRStatPlot);


        TLegend* legendDRComb = GetAndSetLegend2(0.12,0.95-textSizeSinglePad*1,0.5,0.95, textSizeSinglePad, 1, "", 42, 0.15);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        DrawGammaSetMarkerTGraphAsym(graphCombDRSys, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRStatPlot, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV, widthLinesBoxes);
        legendDRComb->AddEntry(graphCombDRSys,"ALICE","pf");
        graphCombDRSys->Draw("E2same");
        graphCombDRStatPlot->Draw("p,same,z");
//         legendDRComb->Draw();

        // labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_Comb_pp8TeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_Comb_pp8TeV.pdf", outputDir.Data()));

    hist2DDRDummySingle->DrawCopy();
        TGraphAsymmErrors* graphCombDRNonFitStatPlot    = (TGraphAsymmErrors*)graphCombDRNonFitStat->Clone("graphCombDRNonFitStatPlot");
        ProduceGraphAsymmWithoutXErrors(graphCombDRNonFitStatPlot);


        legendDRComb = GetAndSetLegend2(0.12,0.95-textSizeSinglePad*1,0.5,0.95, textSizeSinglePad, 1, "", 42, 0.15);
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitSys, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitStatPlot, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV, widthLinesBoxes);
        legendDRComb->AddEntry(graphCombDRNonFitSys,"ALICE","pf");
        graphCombDRNonFitSys->Draw("E2same");
        graphCombDRNonFitStatPlot->Draw("p,same,z");
//         legendDRComb->Draw();

        // labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_Comb_NonFit_pp8TeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_Comb_NonFit_pp8TeV.pdf", outputDir.Data()));

        hist2DDRDummySingle->DrawCopy();

        TLegend* legendDRTheoryComb         = GetAndSetLegend2(0.12,0.96-textSizeSinglePad*1.2,0.5,0.96, textSizeSinglePad, 1, "", 42, 0.15);
        TLegend* legendDRTheoryComb2        = GetAndSetLegend2(0.12,0.94-textSizeSinglePad*5.8,0.5,0.94-textSizeSinglePad*1, textSizeSinglePad, 1, "", 42, 0.15);
        legendDRTheoryComb2->SetTextAlign(12);

        DrawGammaSetMarkerTGraphAsym(graphCombDRSys, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRStatPlot, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV, widthLinesBoxes);
        legendDRTheoryComb->AddEntry(graphCombDRSys,"ALICE","pf");

        if (graphTheoryNLODRpp8TeV) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryNLODRpp8TeV, 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
            graphTheoryNLODRpp8TeV->Draw("3,same");
            legendDRTheoryComb2->AddEntry(dummyVogelsangforLegend,"NLO pQCD, #scale[0.75]{PDF: CT10, FF: GRV}","fl");
        }
        if (graphTheoryNLODRpp8TeVCenter){
            DrawGammaNLOTGraph( graphTheoryNLODRpp8TeVCenter, 2, styleLineNLOWerner, colorNLOWerner);
            graphTheoryNLODRpp8TeVCenter->Draw("lc,same");
        }
        if (graphTheoryNLODRpp8TeVPaquett){
            DrawGammaNLOTGraph( graphTheoryNLODRpp8TeVPaquett, 2, styleLineMcGill, colorNLOMcGill );
            graphTheoryNLODRpp8TeVPaquett->Draw("lc,same");
            legendDRTheoryComb2->AddEntry(graphTheoryNLODRpp8TeVPaquett,"NLO pQCD, #scale[0.75]{PDF: CTEQ6.1M, FF: BFG2}","l");
        }
        if (graphTheoryJETPHOXDRpp8TeV) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryJETPHOXDRpp8TeV, 0, 0, colorJETPHOX, colorJETPHOX, widthLinesBoxes, kTRUE, colorJETPHOX,kTRUE);
            graphTheoryJETPHOXDRpp8TeV->Draw("3,same");
        }
        if (graphTheoryJETPHOXDRpp8TeVCenter){
            DrawGammaNLOTGraph( graphTheoryJETPHOXDRpp8TeVCenter, 2, styleLineJETPHOX, colorJETPHOX);
            graphTheoryJETPHOXDRpp8TeVCenter->Draw("lc,same");
        }
        if (graphTheoryJETPHOXDRpp8TeV) {
            legendDRTheoryComb2->AddEntry(dummyJETPHOXforLegend,"JETPHOX, #scale[0.75]{PDF: NNPDF2.3QED, FF: BFG2}","fl");
            // legendDRTheoryComb2->AddEntry((TObject*)0,"PDF: NNPDF2.3QED, FF: BFG2","");
        }
        if (graphTheoryPOWHEGDRpp8TeV) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryPOWHEGDRpp8TeV, 0, 0, colorPOWHEG, colorPOWHEG, widthLinesBoxes, kTRUE, colorPOWHEG,kTRUE);
            graphTheoryPOWHEGDRpp8TeV->Draw("3,same");
        }
        if (graphTheoryPOWHEGDRpp8TeVCenter){
            DrawGammaNLOTGraph( graphTheoryPOWHEGDRpp8TeVCenter, 2, styleLinePOWHEG, colorPOWHEG);
            graphTheoryPOWHEGDRpp8TeVCenter->Draw("lc,same");
        }
        if (graphTheoryPOWHEGDRpp8TeV) {
            legendDRTheoryComb2->AddEntry(dummyPOWHEGforLegend,"POWHEG, #scale[0.75]{PDF: NNPDF2.3QED + PYTHIA8 PS}","fl");
            // legendDRTheoryComb2->AddEntry((TObject*)0,"PDF: NNPDF2.3QED, FF: BFG2","");
        }

        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRSys->Draw("E2same");
        graphCombDRStatPlot->Draw("p,same,z");
        legendDRTheoryComb->Draw();
        legendDRTheoryComb2->Draw();

        // labelALICEDRSingle->Draw();
        labelDRSingle = new TLatex(0.95,0.92,Form("%s",collisionSystempp8TeV.Data()));
        SetStyleTLatex( labelDRSingle, textSizeSinglePad,4, 1, 42, kTRUE, 31);
        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pp8TeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_pp8TeV.pdf", outputDir.Data()));


        hist2DDRDummySingle->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitSys, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombDRNonFitStatPlot, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV, widthLinesBoxes);

        if (graphTheoryNLODRpp8TeV) {
            DrawGammaSetMarkerTGraphAsym(graphTheoryNLODRpp8TeV, 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
            graphTheoryNLODRpp8TeV->Draw("3,same");
        }
        if (graphTheoryNLODRpp8TeVCenter){
            DrawGammaNLOTGraph( graphTheoryNLODRpp8TeVCenter, 2, styleLineNLOWerner, colorNLOWerner);
            graphTheoryNLODRpp8TeVCenter->Draw("lc,same");
        }
        if (graphTheoryNLODRpp8TeVPaquett){
            graphTheoryNLODRpp8TeVPaquett->Draw("lc,same");
        }
        if (graphTheoryJETPHOXDRpp8TeV) {
            graphTheoryJETPHOXDRpp8TeV->Draw("3,same");
        }
        if (graphTheoryJETPHOXDRpp8TeVCenter){
            graphTheoryJETPHOXDRpp8TeVCenter->Draw("lc,same");
        }
        if (graphTheoryPOWHEGDRpp8TeV) {
            graphTheoryPOWHEGDRpp8TeV->Draw("3,same");
        }
        if (graphTheoryPOWHEGDRpp8TeVCenter){
            graphTheoryPOWHEGDRpp8TeVCenter->Draw("lc,same");
        }
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        graphCombDRNonFitSys->Draw("E2same");
        graphCombDRNonFitStatPlot->Draw("p,same,z");
        legendDRTheoryComb->Draw();
        legendDRTheoryComb2->Draw();

        // labelALICEDRSingle->Draw();

        labelDRSingle->Draw();
        hist2DDRDummySingle->Draw("same,axis");

    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_NonFit_pp8TeV.%s", outputDir.Data(), suffix.Data()));
    canvasDoubleRatio->Print(Form("%s/DR_CombAndTheory_NonFit_pp8TeV.pdf", outputDir.Data()));



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

        TLegend* legendEffiGamma           = GetAndSetLegend2(0.54, 0.13, 0.92, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel, 2);
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
                legendEffiGamma->AddEntry(histoEffiMCPt[i],"    "+nameMeasGlobalLabelGamma[i],"p");
            } else if (histoEffi[i]){
                legendEffiGamma->AddEntry((TObject*)0,"    "+nameMeasGlobalLabelGamma[i],"");
            }
        }
        legendEffiGamma->Draw();

        TLatex *labelPerfEffi           = new TLatex(0.13,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
        labelPerfEffi->Draw();
        TLatex *labelEnergyEffi         = new TLatex(0.13,0.87,collisionSystempp8TeV.Data());
        SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
        labelEnergyEffi->Draw();
        TLatex *labelPerfEffiPTrec      = new TLatex(0.55,0.145+(3*textSizeLabelsRel),"#it{p}_{T}^{rec}");
        SetStyleTLatex( labelPerfEffiPTrec, textSizeLabelsRel,4);
        labelPerfEffiPTrec->Draw();
        TLatex *labelPerfEffiPTtrue      = new TLatex(0.655,0.145+(3*textSizeLabelsRel),"#it{p}_{T}^{true}");
        SetStyleTLatex( labelPerfEffiPTtrue, textSizeLabelsRel,4);
        labelPerfEffiPTtrue->Draw();

    canvasEff->Update();
    canvasEff->Print(Form("%s/Gamma_Effiency_8.%s",outputDir.Data(),suffix.Data()));
    canvasEff->Print(Form("%s/Gamma_Effiency_8.pdf",outputDir.Data()));

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
                legendResolCorGamma->AddEntry(histoResolCorr[i],nameMeasGlobalLabelGamma[i],"p");
            }
        }
        legendResolCorGamma->Draw();
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

        TLatex *labelPerfResolCor           = new TLatex(0.15,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfResolCor, textSizeLabelsRel,4);
        labelPerfResolCor->Draw();
        TLatex *labelEnergyResolCor         = new TLatex(0.15,0.87,collisionSystempp8TeV.Data());
        SetStyleTLatex( labelEnergyResolCor, textSizeLabelsRel,4);
        labelEnergyResolCor->Draw();

    histo2DResCor->Draw("same,axis");
    canvasResolCor->Update();
    canvasResolCor->Print(Form("%s/Gamma_ResolutionCorrection_8.%s",outputDir.Data(),suffix.Data()));
    canvasResolCor->Print(Form("%s/Gamma_ResolutionCorrection_8.pdf",outputDir.Data()));

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
                legendPurityGamma->AddEntry(histoPurity[i],nameMeasGlobalLabelGamma[i],"p");
            }
        }
        legendPurityGamma->Draw();

        TLatex *labelPerfPurity           = new TLatex(0.15,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfPurity, textSizeLabelsRel,4);
        labelPerfPurity->Draw();
        TLatex *labelEnergyPurity         = new TLatex(0.15,0.87,collisionSystempp8TeV.Data());
        SetStyleTLatex( labelEnergyPurity, textSizeLabelsRel,4);
        labelEnergyPurity->Draw();

    canvasPurity->Update();
    canvasPurity->Print(Form("%s/Gamma_Purity_8.%s",outputDir.Data(),suffix.Data()));
    canvasPurity->Print(Form("%s/Gamma_Purity_8.pdf",outputDir.Data()));

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
                legendConvProbGamma->AddEntry(histoConvProb[i],nameMeasGlobalLabelGamma[i],"p");
            }
        }
        legendConvProbGamma->Draw();

        TLatex *labelPerfConvProb           = new TLatex(0.15,0.92,"ALICE simulation");
        SetStyleTLatex( labelPerfConvProb, textSizeLabelsRel,4);
        labelPerfConvProb->Draw();
        TLatex *labelEnergyConvProb         = new TLatex(0.15,0.87,collisionSystempp8TeV.Data());
        SetStyleTLatex( labelEnergyConvProb, textSizeLabelsRel,4);
        labelEnergyConvProb->Draw();

    histo1DConvProb->Draw("same,axis");

    canvasConvProb->Update();
    canvasConvProb->Print(Form("%s/Gamma_ConvProb_8.%s",outputDir.Data(),suffix.Data()));
    canvasConvProb->Print(Form("%s/Gamma_ConvProb_8.pdf",outputDir.Data()));

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
                legendPileUpGamma->AddEntry(histoPileupCorr[i],nameMeasGlobalLabelGamma[i],"p");
            }
        }
        legendPileUpGamma->Draw();

        TLatex *labelPerfPileUp           = new TLatex(0.15,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfPileUp, textSizeLabelsRel,4);
        labelPerfPileUp->Draw();
        TLatex *labelEnergyPileUp         = new TLatex(0.15,0.87,collisionSystempp8TeV.Data());
        SetStyleTLatex( labelEnergyPileUp, textSizeLabelsRel,4);
        labelEnergyPileUp->Draw();

    canvasPileUp->Update();
    canvasPileUp->Print(Form("%s/Gamma_PileUp_8.%s",outputDir.Data(),suffix.Data()));
    canvasPileUp->Print(Form("%s/Gamma_PileUp_8.pdf",outputDir.Data()));

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
    histo1DEffectiveSecCorr->GetXaxis()->SetNoExponent();
    histo1DEffectiveSecCorr->GetXaxis()->SetMoreLogLabels(kTRUE);


    Double_t minYSecCorr[4]             = {0.0, 0.0, 0.0, 0.0};
    Double_t maxYSecCorr[4]             = {0.066, 0.0046, 1.0e-3, 0.054};
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
                legendEffectiveSecCorrGamma->AddEntry(histoEffSecCorr[k][i],nameMeasGlobalLabelGamma[i],"p");
            }
        }
        legendEffectiveSecCorrGamma->Draw();

        TLatex *labelPerfEffectiveSecCorr           = new TLatex(0.15,0.89,"ALICE performance");
        SetStyleTLatex( labelPerfEffectiveSecCorr, textSizeLabelsRel,4);
        labelPerfEffectiveSecCorr->Draw();
        TLatex *labelEnergyEffectiveSecCorr         = new TLatex(0.15,0.84,collisionSystempp8TeV.Data());
        SetStyleTLatex( labelEnergyEffectiveSecCorr, textSizeLabelsRel,4);
        labelEnergyEffectiveSecCorr->Draw();

        canvasEffectiveSecCorr->Update();
        canvasEffectiveSecCorr->Print(Form("%s/Gamma_EffectiveSecCorr_%s_8.%s",outputDir.Data(), nameOutputSec[k].Data(), suffix.Data()));
        canvasEffectiveSecCorr->Print(Form("%s/Gamma_EffectiveSecCorr_%s_8.pdf",outputDir.Data(), nameOutputSec[k].Data()));
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
    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-9,7e0);
    histo2DYieldGamma->GetXaxis()->SetMoreLogLabels();
    histo2DYieldGamma->GetXaxis()->SetNoExponent();
    histo2DYieldGamma->Draw("copy");

        TGraphAsymmErrors* graphCombIncGammaStatPlot    = (TGraphAsymmErrors*)graphCombIncGammaStat->Clone("graphCombIncGammaStatPlot");
        ProduceGraphAsymmWithoutXErrors(graphCombIncGammaStatPlot);

        TLegend* legendYieldIncGamma           = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(3*textSizeLabelsRel*0.85),textSizeLabelsPixel);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV,widthLinesBoxes);
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");
        legendYieldIncGamma->AddEntry(graphCombIncGammaSys,"ALICE","pf");
        legendYieldIncGamma->Draw();

        
        DrawGammaSetMarkerTF1( fitHagGammaComb, 7, 2, colorCombpp8TeV);
        legendYieldIncGamma->AddEntry(fitHagGammaComb,"mod. Hagedorn fit","l");
        fitHagGammaComb->SetRange(0.28, 20);
        fitHagGammaComb->Draw("same");
        DrawGammaSetMarkerTF1( fitTsallisGammaComb, 5, 2, colorCombpp8TeVBox);
        legendYieldIncGamma->AddEntry(fitTsallisGammaComb,"Tsallis fit","l");
        fitTsallisGammaComb->SetRange(doubleRatioXpp[0], doubleRatioXpp[1]);
        fitTsallisGammaComb->Draw("same");
        DrawGammaSetMarkerTF1( fitTCMGammaComb, 4, 2, kBlue+2);
        legendYieldIncGamma->AddEntry(fitTCMGammaComb,"TCM fit","l");
        fitTCMGammaComb->SetRange(0.28, 20);
        fitTCMGammaComb->Draw("same");

        TLatex *labelEnergyInvYieldPaperAll = new TLatex(0.20, 0.20+0.04*3, collisionSystempp8TeV.Data());
        SetStyleTLatex( labelEnergyInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelEnergyInvYieldPaperAll->Draw();
        TLatex *labelALICEInvYieldPaperAll  = new TLatex(0.20,0.20+0.04*2,textALICE.Data());
        SetStyleTLatex( labelALICEInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelALICEInvYieldPaperAll->Draw();
        TLatex *labelALICENormUnPaperAll    = new TLatex(0.20,0.20+0.05*1,"norm. unc. 2.6%");
        SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICENormUnPaperAll->Draw();

        histo2DYieldGamma->Draw("same,axis");

    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_8.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_8.pdf",outputDir.Data()));

    //********************************************************************************************************
    // Plotting simple comparison of data vs fit to inc gamma spec
    //********************************************************************************************************

    TF1* fitTCMDecomposedGammaL                 = new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",0.,0.,0.));
    fitTCMDecomposedGammaL->SetParameters(fitTCMGammaComb->GetParameter(0),fitTCMGammaComb->GetParameter(1));
    fitTCMDecomposedGammaL->SetRange(0.28, 20.);
    TF1 *fitTCMDecomposedGammaH                 = new TF1("twoCompModel_DecH","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]))");
   //      graphCombEtaInvXSectionTotA->Fit(fitTCMDecomposedH,"QNRMEX0+","",5,20);
    fitTCMDecomposedGammaH->SetParameters(fitTCMGammaComb->GetParameter(2),fitTCMGammaComb->GetParameter(3), fitTCMGammaComb->GetParameter(4));
    fitTCMDecomposedGammaH->SetRange(0.28, 20.);

    histo2DYieldGamma->DrawCopy();

    graphCombIncGammaStatPlot->Draw("pEsame");
    graphCombIncGammaSys->Draw("E2same");

    fitTCMGammaComb->SetLineStyle(0);
    fitTCMGammaComb->SetLineColor(kRed+2);
    fitTCMGammaComb->SetRange(0.28, 20.);
    fitTCMGammaComb->Draw("same");

    fitTCMDecomposedGammaL->SetLineColor(kAzure);
    fitTCMDecomposedGammaL->SetLineStyle(2);
    fitTCMDecomposedGammaL->Draw("same");
    fitTCMDecomposedGammaH->SetLineColor(kGreen+2);
    fitTCMDecomposedGammaH->SetLineStyle(8);
    fitTCMDecomposedGammaH->Draw("same");

    TLatex *labelTCMGamma1= new TLatex(0.43, 0.94, Form("TCM low:"));
    TLatex *labelTCMGamma2= new TLatex(0.43, 0.90, Form("A_{1}: (%.1e #pm %.1e) - T_{e}: (%.3f #pm %.3f)",fitTCMGammaComb->GetParameter(0),fitTCMGammaComb->GetParError(0),fitTCMGammaComb->GetParameter(1),fitTCMGammaComb->GetParError(1)));
    TLatex *labelTCMGamma3= new TLatex(0.43, 0.86, Form("TCM high:"));
    TLatex *labelTCMGamma4= new TLatex(0.43, 0.82, Form("A_{2}: (%.1e #pm %.1e) - T: (%.3f #pm %.3f) - n: (%.3f #pm %.3f)",fitTCMGammaComb->GetParameter(2),fitTCMGammaComb->GetParError(2),abs(fitTCMGammaComb->GetParameter(3)),fitTCMGammaComb->GetParError(3),fitTCMGammaComb->GetParameter(4),fitTCMGammaComb->GetParError(4)));

    TLatex *labelTCMGamma5= new TLatex(0.52, 0.75, Form("Bylinkin-Rostovtsev:"));
    TLatex *labelTCMGamma6= new TLatex(0.52, 0.71, Form("#it{A}_{1} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}_{2}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}"));

    SetStyleTLatex( labelTCMGamma1, 0.03,4);
    labelTCMGamma1->Draw();
    SetStyleTLatex( labelTCMGamma2, 0.02,4);
    labelTCMGamma2->Draw();
    SetStyleTLatex( labelTCMGamma3, 0.03,4);
    labelTCMGamma3->Draw();
    SetStyleTLatex( labelTCMGamma4, 0.02,4);
    labelTCMGamma4->Draw();
    SetStyleTLatex( labelTCMGamma5, 0.03,4);
    labelTCMGamma5->Draw();
    SetStyleTLatex( labelTCMGamma6, 0.03,4);
    labelTCMGamma6->Draw();


    TLegend* legendWithFit   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
    legendWithFit->AddEntry(fitTCMDecomposedGammaL,"TCM low","l");
    legendWithFit->AddEntry(fitTCMDecomposedGammaH,"TCM high","l");
    legendWithFit->AddEntry(fitTCMGammaComb,"Bylinkin-Rostovtsev (TCM)","l");
    legendWithFit->Draw();

    canvasInvYieldGamma->Update();
    canvasInvYieldGamma->Print(Form("%s/InvYield_IncGamma_WithFit_8.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_WithFit_8.pdf",outputDir.Data()));

    // **********************************************************************************************************************
    // ******************************** InvYield for individual inc gamma measurements **************************************
    // **********************************************************************************************************************
    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-9,7e0);
    histo2DYieldGamma->Draw("copy");
    TLegend* legendYieldIncGammaInd       = GetAndSetLegend2(0.20, 0.11, 0.5, 0.11+(3*textSizeLabelsRel*0.85),textSizeLabelsPixel);
    for (Int_t i = 0; i < 11; i++){
        if (graphIndGammaIncSys[i]){
            DrawGammaSetMarkerTGraphAsym(graphIndGammaIncSys[i], markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i],widthLinesBoxes, kTRUE);
            graphIndGammaIncSys[i]->Draw("E2same");
            legendYieldIncGammaInd->AddEntry(graphIndGammaIncSys[i], Form("#gamma_{inc} %s", nameMeasGlobalLabelGamma[i].Data()),"pf");
        }
        if (graphIndGammaIncStat[i]){
            DrawGammaSetMarkerTGraphAsym(graphIndGammaIncStat[i],  markerStyleDet[i], markerSizeDet[i], colorDet[i] , colorDet[i]);
            graphIndGammaIncStat[i]->Draw("Epsame,x0");
            if (!graphIncGammaSysErr[i])legendYieldIncGammaInd->AddEntry(graphIndGammaIncStat[i],nameMeasGlobalLabelGamma[i],"p");
        }
    }
    legendYieldIncGammaInd->Draw();

    labelEnergyInvYieldPaperAll->Draw();
    labelALICEInvYieldPaperAll->Draw();
    labelALICENormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_IncGamma_IndMeas_8.%s",outputDir.Data(),suffix.Data()));

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
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErrPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlot, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV);
            graphCombDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrAr , 1, 3, colorCombpp8TeV, colorCombpp8TeV, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr,0.05,kFALSE,graphCombDRStat);
        }


        TLatex *labelEnergyDGInvYieldPaperAll = new TLatex(0.94, 0.965-0.04*1, Form("%s, %s",textALICE.Data(),collisionSystempp8TeV.Data()));
        SetStyleTLatex( labelEnergyDGInvYieldPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
        labelEnergyDGInvYieldPaperAll->Draw();
        TLatex *labelALICEDGNormUnPaperAll    = new TLatex(0.94,0.965-(0.04*2),"norm. unc. 2.6%");
        SetStyleTLatex( labelALICEDGNormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 31);
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_8.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_8.pdf",outputDir.Data()));

    TGraphAsymmErrors* graphCombDirGammaNonFitSpectrumStatErrPlot = NULL;
    if (graphCombDirGammaNonFitSpectrumStatErr) graphCombDirGammaNonFitSpectrumStatErrPlot       = (TGraphAsymmErrors*)graphCombDirGammaNonFitSpectrumStatErr->Clone("graphCombDirGammaNonFitSpectrumStatErrPlot");
    if (graphCombDirGammaNonFitSpectrumStatErrPlot) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaNonFitSpectrumStatErrPlot);

    histo2DYieldGamma->GetYaxis()->SetRangeUser(7e-11,9.9e0);
    histo2DYieldGamma->Draw("copy");


        if (graphCombDirGammaNonFitSpectrumSystErr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaNonFitSpectrumSystErr, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV, widthLinesBoxes, kTRUE);
            graphCombDirGammaNonFitSpectrumSystErr->Draw("E2same");
        }
        if (graphCombDirGammaNonFitSpectrumStatErrPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaNonFitSpectrumStatErrPlot, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV);
            graphCombDirGammaNonFitSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaNonFitSpectrumSumErrAr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaNonFitSpectrumSumErrAr , 1, 3, colorCombpp8TeV, colorCombpp8TeV, 1.8, kTRUE);
            graphCombDirGammaNonFitSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaNonFitSpectrumSumErrAr,0.05,kFALSE,graphCombDRStat);
        }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_NonFit_8.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_NonFit_8.pdf",outputDir.Data()));
    histo2DYieldGamma->Draw("copy");

        TLegend* legendYieldDirGamma        = GetAndSetLegend2(0.70, 0.84-(2*textSizeLabelsRel*0.85), 0.93, 0.84,textSizeLabelsPixel, 1, "", 43, 0.3);

        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys, markerStyleCombpp8TeV+4, markerSizeCombpp8TeV+0.2, colorCombpp8TeVBox , colorCombpp8TeVBox,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot, markerStyleCombpp8TeV+4, markerSizeCombpp8TeV+0.2, colorCombpp8TeVBox , colorCombpp8TeVBox, widthLinesBoxes);
        legendYieldDirGamma->AddEntry(graphCombIncGammaSys, "#gamma_{inc} ALICE","pf");
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");

        if (graphCombDirGammaSpectrumSystErr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr->Draw("E2same");
            legendYieldDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr, "#gamma_{dir} ALICE","pf");
        }
        if (graphCombDirGammaSpectrumStatErrPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErrPlot, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV);
            graphCombDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErrAr , 1, 3, colorCombpp8TeV, colorCombpp8TeV, 1.8, kTRUE);
            graphCombDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr,0.05,kFALSE,graphCombDRStat);
        }

        // TGraphAsymmErrors* dummyForLegend    = new TGraphAsymmErrors(1);
        // dummyForLegend->SetPoint(0,6.2,0.045);
        // dummyForLegend->SetPointError(0,0,0,0.022,0);
        // DrawGammaSetMarkerTGraphAsym(dummyForLegend , 1, 3, colorCombpp8TeV, colorCombpp8TeV, 1.8, kTRUE);
        // if (graphCombDirGammaSpectrumSumErrAr){
        //     dummyForLegend->Draw(">,same");
        //     PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend);
        // }
        // legendYieldDirGamma->Draw();

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_8.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_8.pdf",outputDir.Data()));
    histo2DYieldGamma->Draw("copy");

        legendYieldDirGamma        = GetAndSetLegend2(0.70, 0.84-(2*textSizeLabelsRel*0.85), 0.93, 0.84,textSizeLabelsPixel, 1, "", 43, 0.3);

        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSys, markerStyleCombpp8TeV+4, markerSizeCombpp8TeV+0.2, colorCombpp8TeVBox , colorCombpp8TeVBox,widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatPlot, markerStyleCombpp8TeV+4, markerSizeCombpp8TeV+0.2, colorCombpp8TeVBox , colorCombpp8TeVBox, widthLinesBoxes);
        legendYieldDirGamma->AddEntry(graphCombIncGammaSys, "#gamma_{inc} ALICE","pf");
        graphCombIncGammaSys->Draw("E2same");
        graphCombIncGammaStatPlot->Draw("Epsame");

        if (graphCombDirGammaNonFitSpectrumSystErr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaNonFitSpectrumSystErr, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV, widthLinesBoxes, kTRUE);
            graphCombDirGammaNonFitSpectrumSystErr->Draw("E2same");
            legendYieldDirGamma->AddEntry(graphCombDirGammaNonFitSpectrumSystErr, "#gamma_{dir} ALICE","pf");
        }
        if (graphCombDirGammaNonFitSpectrumStatErrPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaNonFitSpectrumStatErrPlot, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV);
            graphCombDirGammaNonFitSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaNonFitSpectrumSumErrAr){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaNonFitSpectrumSumErrAr , 1, 3, colorCombpp8TeV, colorCombpp8TeV, 1.8, kTRUE);
            graphCombDirGammaNonFitSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaNonFitSpectrumSumErrAr,0.05,kFALSE,graphCombDRStat);
        }


        // if (graphCombDirGammaNonFitSpectrumSumErrAr){
        //     dummyForLegend->Draw(">,same");
        //     PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend);
        // }
        legendYieldDirGamma->Draw();

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_NonFit_8.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_NonFit_8.pdf",outputDir.Data()));
  
    histo2DYieldGamma->Draw("copy");

    TLegend* legendYieldDirGammaTheo2      = GetAndSetLegend2(0.20, 0.1, 0.5, 0.1+(6*textSizeLabelsRel*0.87),textSizeLabelsPixel,1, "", 43, 0.23);
    legendYieldDirGammaTheo2->SetTextAlign(12);
    graphCombIncGammaSys->Draw("E2same");
    graphCombIncGammaStatPlot->Draw("Epsame");

        if (graphCombDirGammaSpectrumSystErr){
            graphCombDirGammaSpectrumSystErr->Draw("E2same");
        }
        legendYieldDirGamma->Draw();
        legendYieldDirGammaTheo2->AddEntry((TObject*)0,"NLO pQCD","");
        if (graphTheoryNLOpp8TeV) {
          // DrawGammaSetMarkerTGraphAsym(graphTheoryNLOpp8TeV, 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
          // graphTheoryNLOpp8TeV->Draw("3,same");
            DrawGammaNLOTGraph( graphTheoryNLOpp8TeVCenter, 2, styleLineNLOWerner, colorNLOWerner );
            graphTheoryNLOpp8TeVCenter->Draw("lX,same");
            legendYieldDirGammaTheo2->AddEntry(dummyVogelsangforLegend,"#scale[0.75]{PDF: CT10, FF: GRV}","l");
        }

        if (graphTheoryNLOp8TeVPaquett) {
            DrawGammaNLOTGraph( graphTheoryNLOp8TeVPaquett, 2, styleLineMcGill, colorNLOMcGill );
            graphTheoryNLOp8TeVPaquett->Draw("lc,same");
            legendYieldDirGammaTheo2->AddEntry(graphTheoryNLOp8TeVPaquett,"#scale[0.75]{PDF: CTEQ6.1M, FF: BFG2}","l");
        }
        if (graphTheoryJETPHOXpp8TeV) {
          DrawGammaSetMarkerTGraphAsym(graphTheoryJETPHOXpp8TeV, 0, 0, colorJETPHOX, colorJETPHOX, widthLinesBoxes, kTRUE, colorJETPHOX,kTRUE);
            graphTheoryJETPHOXpp8TeV->Draw("3,same");
            graphTheoryJETPHOXpp8TeV->SetLineWidth(widthLineNLO*1.5);
        }
        if (graphTheoryJETPHOXpp8TeVCenter) {
          DrawGammaNLOTGraph( graphTheoryJETPHOXpp8TeVCenter, 2, styleLineJETPHOX, colorJETPHOX );
            graphTheoryJETPHOXpp8TeVCenter->Draw("lX,same");
            graphTheoryJETPHOXpp8TeVCenter->SetLineWidth(widthLineNLO*1.5);
        }
        if (graphTheoryJETPHOXpp8TeV) {
            legendYieldDirGammaTheo2->AddEntry((TObject*)0,"JETPHOX","");
            legendYieldDirGammaTheo2->AddEntry(dummyJETPHOXforLegend,"#scale[0.75]{PDF: NNPDF2.3QED, FF: BFG2}","fl");
        }
        if (graphTheoryPOWHEGpp8TeV) {
          DrawGammaSetMarkerTGraphAsym(graphTheoryPOWHEGpp8TeV, 0, 0, colorPOWHEG, colorPOWHEG, widthLinesBoxes, kTRUE, colorPOWHEG,kTRUE);
            graphTheoryPOWHEGpp8TeV->Draw("3,same");
            graphTheoryPOWHEGpp8TeV->SetLineWidth(widthLineNLO*1.5);
        }
        if (graphTheoryPOWHEGpp8TeVCenter) {
          DrawGammaNLOTGraph( graphTheoryPOWHEGpp8TeVCenter, 2, styleLinePOWHEG, colorPOWHEG );
            graphTheoryPOWHEGpp8TeVCenter->Draw("lX,same");
            graphTheoryPOWHEGpp8TeVCenter->SetLineWidth(widthLineNLO*1.5);
        }
        if (graphTheoryPOWHEGpp8TeV) {
            legendYieldDirGammaTheo2->AddEntry((TObject*)0,"POWHEG","");
            legendYieldDirGammaTheo2->AddEntry(dummyPOWHEGforLegend,"#scale[0.75]{PDF: NNPDF2.3QED + PYTHIA8 PS}","fl");
        }
        legendYieldDirGammaTheo2->Draw();

        if (graphCombDirGammaSpectrumStatErrPlot){
            graphCombDirGammaSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaSpectrumSumErrAr){
            graphCombDirGammaSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErrAr,0.05,kFALSE,graphCombDRStat);
        }
        // if (graphCombDirGammaSpectrumSumErrAr){
        //     dummyForLegend->Draw(">,same");
        //     PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend);
        // }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory_8.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory_8.pdf",outputDir.Data()));
    histo2DYieldGamma->Draw("copy");

    
    graphCombIncGammaSys->Draw("E2same");
    graphCombIncGammaStatPlot->Draw("Epsame");

        if (graphCombDirGammaNonFitSpectrumSystErr){
            graphCombDirGammaNonFitSpectrumSystErr->Draw("E2same");
        }
        legendYieldDirGamma->Draw();

        if (graphTheoryNLOpp8TeV) {
            // graphTheoryNLOpp8TeV->Draw("3,same");
            graphTheoryNLOpp8TeVCenter->Draw("lX,same");
        }
        if (graphTheoryNLOp8TeVPaquett) {
            graphTheoryNLOp8TeVPaquett->Draw("lc,same");
        }
        if (graphTheoryJETPHOXpp8TeV) {
            graphTheoryJETPHOXpp8TeV->Draw("3,same");
            graphTheoryJETPHOXpp8TeVCenter->Draw("lX,same");
        }
        if (graphTheoryPOWHEGpp8TeV) {
            graphTheoryPOWHEGpp8TeV->Draw("3,same");
            graphTheoryPOWHEGpp8TeVCenter->Draw("lX,same");
        }
        legendYieldDirGammaTheo2->Draw();

        if (graphCombDirGammaNonFitSpectrumStatErrPlot){
            graphCombDirGammaNonFitSpectrumStatErrPlot->Draw("p,E1Z,same");
        }
        if (graphCombDirGammaNonFitSpectrumSumErrAr){
            graphCombDirGammaNonFitSpectrumSumErrAr->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaNonFitSpectrumSumErrAr,0.05,kFALSE,graphCombDRStat);
        }
        // if (graphCombDirGammaNonFitSpectrumSumErrAr){
        //     dummyForLegend->Draw(">,same");
        //     PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend);
        // }

        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();

    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory_NonFit_8.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvYield_DirGamma_IncGamma_Theory_NonFit_8.pdf",outputDir.Data()));

    TGraphAsymmErrors* graphXSecCombDirGammaNonFitSpectrumSystErr       = ScaleGraph( graphCombDirGammaNonFitSpectrumSystErr,xSection8TeVV0AND*recalcBarn);
    TGraphAsymmErrors* graphXSecCombDirGammaNonFitSpectrumStatErrPlot   = ScaleGraph( graphCombDirGammaNonFitSpectrumStatErrPlot,xSection8TeVV0AND*recalcBarn);
    TGraphAsymmErrors* graphXSecCombDirGammaNonFitSpectrumSumErrAr      = ScaleGraph( graphCombDirGammaNonFitSpectrumSumErrAr,xSection8TeVV0AND*recalcBarn);
    TGraphAsymmErrors* graphXSecCombIncGammaSys                         = ScaleGraph( graphCombIncGammaSys,xSection8TeVV0AND*recalcBarn);
    TGraphAsymmErrors* graphXSecCombIncGammaStatPlot                    = ScaleGraph( graphCombIncGammaStatPlot,xSection8TeVV0AND*recalcBarn);
    TGraphAsymmErrors* graphXSecTheoryNLOpp8TeV                         = ScaleGraph( graphTheoryNLOpp8TeV,xSection8TeVV0AND*recalcBarn);
    TGraphAsymmErrors* graphXSecTheoryNLOpp8TeVCenter                   = ScaleGraph( graphTheoryNLOpp8TeV,xSection8TeVV0AND*recalcBarn);
    TGraphAsymmErrors* graphXSecTheoryJETPHOXpp8TeV                     = ScaleGraph( graphTheoryJETPHOXpp8TeV,xSection8TeVV0AND*recalcBarn);
    TGraphAsymmErrors* graphXSecTheoryJETPHOXpp8TeVCenter               = ScaleGraph( graphTheoryJETPHOXpp8TeVCenter,xSection8TeVV0AND*recalcBarn);
    TGraphAsymmErrors* graphXSecTheoryPOWHEGpp8TeV                     = ScaleGraph( graphTheoryPOWHEGpp8TeV,xSection8TeVV0AND*recalcBarn);
    TGraphAsymmErrors* graphXSecTheoryPOWHEGpp8TeVCenter               = ScaleGraph( graphTheoryPOWHEGpp8TeVCenter,xSection8TeVV0AND*recalcBarn);
    TGraph* graphXSecTheoryPaquettpp8TeV                                = ScaleGraph( graphTheoryNLOp8TeVPaquett,xSection8TeVV0AND*recalcBarn);
    
    TH1F * histo2DXSecGamma              = new TH1F("histo2DXSecGamma","histo2DXSecGamma",11000,doubleRatioXpp[0], doubleRatioXpp[1]);
    SetStyleHistoTH1ForGraphs(histo2DXSecGamma, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.7);
    histo2DXSecGamma->GetYaxis()->SetRangeUser(2e1,7e11);
    histo2DXSecGamma->GetXaxis()->SetMoreLogLabels();
    histo2DXSecGamma->GetXaxis()->SetNoExponent();
    histo2DXSecGamma->Draw("copy");
    
    DrawGammaSetMarkerTGraphAsym(graphXSecCombIncGammaSys, markerStyleCombpp8TeV+4, markerSizeCombpp8TeV+0.2, colorCombpp8TeVBox , colorCombpp8TeVBox,widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphXSecCombIncGammaStatPlot, markerStyleCombpp8TeV+4, markerSizeCombpp8TeV+0.2, colorCombpp8TeVBox , colorCombpp8TeVBox, widthLinesBoxes);
    graphXSecCombIncGammaSys->Draw("E2same");
    graphXSecCombIncGammaStatPlot->Draw("Epsame");

    TLegend* legendULXSecDirGamma      = GetAndSetLegend2(0.2, 0.12+(6.2*textSizeLabelsRel*0.85), 0.5, 0.12+(7.2*textSizeLabelsRel*0.85),textSizeLabelsPixel,1, "", 43, 0.23);
    legendULXSecDirGamma->SetTextAlign(12);
    
    TLegend* legendTCMfitXSec      = GetAndSetLegend2(0.68, 0.87-(2*textSizeLabelsRel*0.85), 0.93, 0.87,textSizeLabelsPixel, 1, "", 43, 0.3);
    legendTCMfitXSec->AddEntry(graphCombIncGammaSys, "#gamma_{inc} data","pf");
    
    TF1* fitXSecTCMGammaComb = ScaleTF1(fitTCMGammaComb,xSection8TeVV0AND*recalcBarn,"xsecTCMfit");
    DrawGammaSetMarkerTF1( fitXSecTCMGammaComb, 4, 2, kBlue+2);
    legendTCMfitXSec->AddEntry(fitXSecTCMGammaComb,"#gamma_{inc} TCM fit","l");
    fitXSecTCMGammaComb->SetRange(0.28, 17);
    fitXSecTCMGammaComb->Draw("same");
    legendTCMfitXSec->Draw();

    if (graphXSecTheoryNLOpp8TeVCenter) {
        // DrawGammaSetMarkerTGraphAsym(graphXSecTheoryNLOpp8TeV, 0, 0, colorNLOWernerBand, colorNLOWernerBand, 0.2, kTRUE, colorNLOWernerBand);
        // graphXSecTheoryNLOpp8TeV->Draw("3,same");
        // graphXSecTheoryNLOpp8TeV->SetLineWidth(widthLineNLO*1.5);
        DrawGammaNLOTGraph( graphXSecTheoryNLOpp8TeVCenter, 2, styleLineNLOWerner, colorNLOWerner );
        graphXSecTheoryNLOpp8TeVCenter->Draw("lX,same");
    }
    if (graphXSecTheoryPaquettpp8TeV) {
        DrawGammaNLOTGraph( graphXSecTheoryPaquettpp8TeV, 2, styleLineMcGill, colorNLOMcGill );
        graphXSecTheoryPaquettpp8TeV->Draw("lc,same");
    }
    if (graphXSecTheoryJETPHOXpp8TeV) {
      DrawGammaSetMarkerTGraphAsym(graphXSecTheoryJETPHOXpp8TeV, 0, 0, colorJETPHOX, colorJETPHOX, widthLinesBoxes, kTRUE, colorJETPHOX,kTRUE);
        graphXSecTheoryJETPHOXpp8TeV->Draw("3,same");
        graphXSecTheoryJETPHOXpp8TeV->SetLineWidth(widthLineNLO*1.5);
    }
    if (graphXSecTheoryJETPHOXpp8TeVCenter) {
        DrawGammaNLOTGraph( graphXSecTheoryJETPHOXpp8TeVCenter, 2, styleLineJETPHOX, colorJETPHOX );
          graphXSecTheoryJETPHOXpp8TeVCenter->Draw("lX,same");
          graphXSecTheoryJETPHOXpp8TeVCenter->SetLineWidth(widthLineNLO*1.5);
    }
    if (graphXSecTheoryPOWHEGpp8TeV) {
      DrawGammaSetMarkerTGraphAsym(graphXSecTheoryPOWHEGpp8TeV, 0, 0, colorPOWHEG, colorPOWHEG, widthLinesBoxes, kTRUE, colorPOWHEG,kTRUE);
        graphXSecTheoryPOWHEGpp8TeV->Draw("3,same");
        graphXSecTheoryPOWHEGpp8TeV->SetLineWidth(widthLineNLO*1.5);
    }
    if (graphXSecTheoryPOWHEGpp8TeVCenter) {
        DrawGammaNLOTGraph( graphXSecTheoryPOWHEGpp8TeVCenter, 2, styleLinePOWHEG, colorPOWHEG );
          graphXSecTheoryPOWHEGpp8TeVCenter->Draw("lX,same");
          graphXSecTheoryPOWHEGpp8TeVCenter->SetLineWidth(widthLineNLO*1.5);
    }
    if (graphXSecCombDirGammaNonFitSpectrumSystErr){
        DrawGammaSetMarkerTGraphAsym(graphXSecCombDirGammaNonFitSpectrumSystErr, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV, widthLinesBoxes, kTRUE);
        graphXSecCombDirGammaNonFitSpectrumSystErr->Draw("E2same");
        legendULXSecDirGamma->AddEntry(graphXSecCombDirGammaNonFitSpectrumSystErr, "#gamma_{dir} data","pf");
  }
    if (graphXSecCombDirGammaNonFitSpectrumStatErrPlot){
        DrawGammaSetMarkerTGraphAsym(graphXSecCombDirGammaNonFitSpectrumStatErrPlot, markerStyleCombpp8TeV, markerSizeCombpp8TeV, colorCombpp8TeV , colorCombpp8TeV);
        graphXSecCombDirGammaNonFitSpectrumStatErrPlot->Draw("p,same");//Draw("p,E1Z,same");
    }
    if (graphXSecCombDirGammaNonFitSpectrumSumErrAr){
        DrawGammaSetMarkerTGraphAsym(graphXSecCombDirGammaNonFitSpectrumSumErrAr , 1, 3, colorCombpp8TeV, colorCombpp8TeV, 1.8, kTRUE);
        graphXSecCombDirGammaNonFitSpectrumSumErrAr->Draw(">,same");
        PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphXSecCombDirGammaNonFitSpectrumSumErrAr,0.05,kFALSE,graphCombDRStat);
    }

    // dummyForLegend->SetPoint(0,0.32,1e4);
    // dummyForLegend->SetPointError(0,0,0,5e3,0);
    //     if (graphXSecCombDirGammaNonFitSpectrumSumErrAr){
    //         dummyForLegend->Draw(">,same");
    //         PlotErrorBarAtUpperEdgeOfTGraphAsymErr(dummyForLegend,0.015);
    //     }
        legendULXSecDirGamma->Draw();
        legendYieldDirGammaTheo2->Draw();
        labelEnergyDGInvYieldPaperAll->Draw();
        labelALICEDGNormUnPaperAll->Draw();
        
        
    histo2DYieldGamma->Draw("same,axis");
    canvasInvYieldGamma->SaveAs(Form("%s/InvXsection_DirGamma_IncGamma_NonFit_Theory_8.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldGamma->SaveAs(Form("%s/InvXsection_DirGamma_IncGamma_NonFit_Theory_8.pdf",outputDir.Data()));
    //***************************** Plot cocktail gammas to all gammas ratio ****************************************
    Style_t     cocktailColorPartialSums[14]        = {kRed+2,kBlue+1,kYellow+2,kOrange+1,kAzure-2,kRed-2,kViolet,kGreen-3,kOrange+6, kTeal+9,kMagenta-3,kCyan+4,kViolet+4,kAzure-4};

    TDirectory* directoryCocktailpp8TeV[11]     = {NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL};
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
    TFile* fileCocktailGamma8TeV                    = new TFile( inputFileNameCocktail.Data());
    directoryCocktailpp8TeV[0]       = (TDirectory*)fileCocktailGamma8TeV->Get("8TeV_PCM");
    directoryCocktailpp8TeV[2]       = (TDirectory*)fileCocktailGamma8TeV->Get("8TeV_EMCal");
    directoryCocktailpp8TeV[4]       = (TDirectory*)fileCocktailGamma8TeV->Get("8TeV_PCMEMCal");
    directoryCocktailpp8TeV[5]       = (TDirectory*)fileCocktailGamma8TeV->Get("8TeV_Comb");
    for(Int_t k=0;k<11;k++){
        if(directoryCocktailpp8TeV[k]){
            for(Int_t i=0;i<14;i++){
                histoGammaFractionsCocktail[i][k]                 = (TH1D*) directoryCocktailpp8TeV[k]->Get(Form("Gamma_From_%s_Pt_OrBin_RatioToAll",fParticle[i].Data()));
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
            graphGammaFractionsCocktailBand[i]->SetFillColor(cocktailColor[i]);
            graphGammaFractionsCocktailBand[i]->SetFillStyle(4050);
            graphGammaFractionsCocktailBand[i]->Draw("3same");
            if(i<4)
            graphGammaFractionsCocktailBand[i]->Draw("LXsame");

        }
    }
    legendGammasRatio2->Draw("same");
    PutProcessLabelAndEnergyOnPlot(                 0.17, 0.22, 40, "", collisionSystempp8TeV.Data(), "", 43, 0.03);
    PutALICESimulationLabel(                   0.17, 0.15, 40, 0.03, 1.25, 43);
    dummyHist->Draw("same,axis");

    canvasGammasRatio2->SaveAs(Form("%s/CocktailGammasRatioToAll_8.%s",outputDir.Data(),suffix.Data()));
    canvasGammasRatio2->SaveAs(Form("%s/CocktailGammasRatioToAll_8.eps",outputDir.Data()));

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
    PutProcessLabelAndEnergyOnPlot(                 0.17, 0.22, 40, "", collisionSystempp8TeV.Data(), "", 43, 0.03);
    PutALICESimulationLabel(                   0.17, 0.15, 40, 0.03, 1.25, 43);
    dummyHist->Draw("same,axis");

    canvasGammasRatio2->SaveAs(Form("%s/CocktailGammasRatioToAllPartialSums_8.%s",outputDir.Data(),suffix.Data()));
    canvasGammasRatio2->SaveAs(Form("%s/CocktailGammasRatioToAllPartialSums_8.eps",outputDir.Data()));

    // ****************************************************************************************************************
    // ************************** Store final results including corr factors in 1 file ********************************
    // ****************************************************************************************************************
    TString nameOutputCommonFile    = Form("CombinedGammaResultPP8TeV_%s.root", dateForOutput.Data());

    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    fCombResults.mkdir("Gamma8TeV");
    TDirectoryFile* directoryGamma = (TDirectoryFile*)fCombResults.Get("Gamma8TeV");
    fCombResults.cd("Gamma8TeV");

        // writing main results
        if (graphCombDRStat) graphCombDRStat->Write("graphRGammaCombStatErr");
        if (graphCombDRSys) graphCombDRSys->Write("graphRGammaCombSysErr");
        if (graphCombDRTot) graphCombDRTot->Write("graphRGammaCombTotErr");
        if (graphCombDRNonFitStat) graphCombDRNonFitStat->Write("graphRGammaCombNonFitStatErr");
        if (graphCombDRNonFitSys) graphCombDRNonFitSys->Write("graphRGammaCombNonFitSysErr");
        if (graphCombDRNonFitTot) graphCombDRNonFitTot->Write("graphRGammaCombNonFitTotErr");
        if (graphCombIncGammaStat) graphCombIncGammaStat->Write("graphInvYieldIncGammaStatErr");
        if (graphCombIncGammaSys) graphCombIncGammaSys->Write("graphInvYieldIncGammaSysErr");
        if (graphCombIncGammaTot) graphCombIncGammaTot->Write("graphInvYieldIncGammaTotErr");
        if (graphCombDirGammaSpectrumSystErr) graphCombDirGammaSpectrumSystErr->Write("graphInvYieldDirGammaSysErr");
        if (graphCombDirGammaSpectrumStatErr) graphCombDirGammaSpectrumStatErr->Write("graphInvYieldDirGammaStatErr");
        if (graphCombDirGammaSpectrumSumErrAr) graphCombDirGammaSpectrumSumErrAr->Write("graphInvYieldDirGammaSumErrAr");
        if (graphCombDirGammaNonFitSpectrumSystErr) graphCombDirGammaNonFitSpectrumSystErr->Write("graphInvYieldDirGammaNonFitSysErr");
        if (graphCombDirGammaNonFitSpectrumStatErr) graphCombDirGammaNonFitSpectrumStatErr->Write("graphInvYieldDirGammaNonFitStatErr");
        if (graphCombDirGammaNonFitSpectrumSumErrAr) graphCombDirGammaNonFitSpectrumSumErrAr->Write("graphInvYieldDirGammaNonFitSumErrAr");

        for (Int_t i = 0; i < 11; i++){
            if (graphIndGammaIncSys[i]) graphIndGammaIncSys[i]->Write(Form("graphInvYieldIncGamma%sSysErr",nameMeasGlobal[i].Data()));
            if (graphIndGammaIncStat[i]) graphIndGammaIncStat[i]->Write(Form("graphInvYieldIncGamma%sStatErr",nameMeasGlobal[i].Data()));
            if (histoIncGammaStatErr[i]) histoIncGammaStatErr[i]->Write(Form("histoInvYieldIncGamma%sStatErr_Unshifted",nameMeasGlobal[i].Data()));
            if (histoDRPi0FitStatErr[i]) histoDRPi0FitStatErr[i]->Write(Form("histoRGamma%sStatErr",nameMeasGlobal[i].Data()));
            if (graphDRPi0FitSysErr[i]) graphDRPi0FitSysErr[i]->Write(Form("graphRGamma%sSysErr",nameMeasGlobal[i].Data()));
            if (histoDRStatErr[i]) histoDRStatErr[i]->Write(Form("histoRGammaNonFit%sStatErr",nameMeasGlobal[i].Data()));
            if (graphDRSysErr[i]) graphDRSysErr[i]->Write(Form("graphRGammaNonFit%sSysErr",nameMeasGlobal[i].Data()));
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
        directoryGamma->mkdir("Fits");
        directoryGamma->cd("Fits");
            // Writing inclusive gamma fit
                if (fitTCMGammaComb) fitTCMGammaComb->Write("TwoComponentModelFitGamma");
                if (fitHagGammaComb) fitHagGammaComb->Write("ModHagedornFitGamma");
                if (fitTsallisGammaComb) fitTsallisGammaComb->Write("TsallisFitGamma");

    fCombResults.Close();

  if(writeFitsToPi0File){
    TFile fCombResultsPi0(Form("%s%s.root",inputFileFits.Data(),addNamePi0Output.Data()), "UPDATE");
      if (fitTCMGammaComb) fitTCMGammaComb->Write("TwoComponentModelFitGamma");
      if (fitHagGammaComb) fitHagGammaComb->Write("ModHagedornFitGamma");
      if (fitTsallisGammaComb) fitTsallisGammaComb->Write("TsallisFitGamma");
    fCombResultsPi0.Close();
  }
}
