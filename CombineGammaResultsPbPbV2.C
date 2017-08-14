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
                                TString inputFileNamePHOS   = "ExternalInputPbPb/PHOS/PHOS_GammaDir_final_09022015.root",
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
    TString outputDir                                           = Form("%s/%s/CombineGammaMeasurementsPbPb",suffix.Data(),dateForOutput.Data());
    TString fileNameTheoryPbPb                                  = "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";
    TString fileNameExperimentPbPb                              = "ExternalInputPbPb/OtherExperiments/phenix_200.root";
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
    Double_t doubleRatioX[2];
    Double_t doubleRatioXpp[2];
    doubleRatio[0]      = 0.75;     doubleRatio[1]      = 1.95;
    indMeasRatio[0]     = 0.65;     indMeasRatio[1]     = 1.45;
//  incRatio[0]         = 0.0;      incRatio[1]         = 1.7;
    doubleRatioX[0]     = 0.7;      doubleRatioX[1]     = 16;
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
    Color_t colorNLOMcGill                                      = kOrange+1;
    Style_t styleNLOMcGill                                      = 7;
    Color_t colorPHSD                                           = kGray+2;
    Style_t stylePHSD                                           = 2;
    Color_t colorChatterjee                                     = kBlue+1;
    Style_t styleChatterjee                                     = 4;
    Color_t colorHees                                           = kBlack;
    Style_t styleHees                                           = 3;
    Color_t colorHe                                             = kBlue-7;
    Style_t styleHe                                             = 8;
    Style_t styleFit                                            = 1;

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

    TString collisionSystem                                     = "Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString collisionSystemCent0010                             = "0-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString collisionSystemCent2040                             = "20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString collisionSystemCent2050                             = "20-50% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString collisionSystempPb                                  = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    TString collisionSystempp7TeV                               = "pp, #sqrt{#it{s}} = 7 TeV";
    TString collisionSystempp2760GeV                            = "pp, #sqrt{#it{s}} = 2.76 TeV";
    TString collisionSystemRHIC0010                             = "0-10% Au-Au #sqrt{#it{s}_{_{NN}}} = 0.2 TeV";
    TString collisionSystemRHIC2040                             = "20-40% Au-Au #sqrt{#it{s}_{_{NN}}} = 0.2 TeV";


    TString PubCombGammaFileName = "FinalResults/Gamma_CombResults_PbPb_2.76TeV_20150729.root";
    TString PubPCMGammaFileName  = "FinalResults/Gamma_PCMResults_PbPb_2.76TeV_20150729.root";
    TFile *filePubGammaPCM       = new TFile(PubPCMGammaFileName.Data());
    TFile *filePubGammaComb      = new TFile(PubCombGammaFileName.Data());
    if (filePubGammaPCM->IsZombie() || filePubGammaComb->IsZombie()) {
        cout << "published gamma file couldn't be read, aborting....";
        return;
    }

    TDirectoryFile* directoryCombGamma_0020                         = (TDirectoryFile*)filePubGammaComb->Get("Gamma_PbPb_2.76TeV_0-20%");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumStat_0020    = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DirGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumSyst_0020   = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DirGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumStat_0020   = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("IncGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumSyst_0020  = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("IncGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioStat_0020         = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DR_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioSyst_0020         = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DR_comb_SysErr");

    TDirectoryFile* directoryCombGamma_2040                         = (TDirectoryFile*)filePubGammaComb->Get("Gamma_PbPb_2.76TeV_20-40%");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumStat_2040    = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DirGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumSyst_2040   = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DirGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumStat_2040   = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("IncGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumSyst_2040  = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("IncGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioStat_2040         = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DR_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioSyst_2040         = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DR_comb_SysErr");

    TDirectoryFile* directoryPCMGamma_0020                          = (TDirectoryFile*)filePubGammaPCM->Get("Gamma_PbPb_2.76TeV_0-20%");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumStat_0020     = (TGraphAsymmErrors*)directoryPCMGamma_0020->Get("graphDirGammaSpectrumStat");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumSyst_0020     = (TGraphAsymmErrors*)directoryPCMGamma_0020->Get("graphDirGammaSpectrumSyst");
        TH1D *histoPubPCMInclGammaSpectrumStat_0020                 = (TH1D*)directoryPCMGamma_0020->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumStat_0020    = new TGraphAsymmErrors(histoPubPCMInclGammaSpectrumStat_0020);
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumSyst_0020    = (TGraphAsymmErrors*)directoryPCMGamma_0020->Get("IncGammaSystError");
        TH1D *histoPubPCMInclRatioSpectrumStat_0020                 = (TH1D*)directoryPCMGamma_0020->Get("IncRatioStatError");
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumStat_0020    = new TGraphAsymmErrors(histoPubPCMInclRatioSpectrumStat_0020);
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumSyst_0020    = (TGraphAsymmErrors*)directoryPCMGamma_0020->Get("IncRatioSystError");
        TH1D *histoPubPCMDoubleRatioStat_0020                       = (TH1D*)directoryPCMGamma_0020->Get("DoubleRatioStatError");
        TGraphAsymmErrors* graphPubPCMDoubleRatioStat_0020          = new TGraphAsymmErrors(histoPubPCMDoubleRatioStat_0020);
        TGraphAsymmErrors* graphPubPCMDoubleRatioSyst_0020          = (TGraphAsymmErrors*)directoryPCMGamma_0020->Get("DoubleRatioSystError");

    TDirectoryFile* directoryPCMGamma_0010                         = (TDirectoryFile*)filePubGammaPCM->Get("Gamma_PbPb_2.76TeV_0-10%");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumStat_0010    = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("graphDirGammaSpectrumStat");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumSyst_0010   = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("graphDirGammaSpectrumSyst");
        TH1D *histoPubPCMInclGammaSpectrumStat_0010                 = (TH1D*)directoryPCMGamma_0010->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumStat_0010    = new TGraphAsymmErrors(histoPubPCMInclGammaSpectrumStat_0010);
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumSyst_0010    = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("IncGammaSystError");
        TH1D *histoPubPCMInclRatioSpectrumStat_0010                 = (TH1D*)directoryPCMGamma_0010->Get("IncRatioStatError");
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumStat_0010    = new TGraphAsymmErrors(histoPubPCMInclRatioSpectrumStat_0010);
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumSyst_0010    = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("IncRatioSystError");
        TH1D *histoPubPCMDoubleRatioStat_0010                       = (TH1D*)directoryPCMGamma_0010->Get("DoubleRatioStatError");
        TGraphAsymmErrors* graphPubPCMDoubleRatioStat_0010          = new TGraphAsymmErrors(histoPubPCMDoubleRatioStat_0010);
        TGraphAsymmErrors* graphPubPCMDoubleRatioSyst_0010          = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("DoubleRatioSystError");

    TDirectoryFile* directoryPCMGamma_2040                         = (TDirectoryFile*)filePubGammaPCM->Get("Gamma_PbPb_2.76TeV_20-40%");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumStat_2040    = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("graphDirGammaSpectrumStat");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumSyst_2040   = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("graphDirGammaSpectrumSyst");
        TH1D *histoPubPCMInclGammaSpectrumStat_2040                 = (TH1D*)directoryPCMGamma_2040->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumStat_2040    = new TGraphAsymmErrors(histoPubPCMInclGammaSpectrumStat_2040);
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumSyst_2040    = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("IncGammaSystError");
        TH1D *histoPubPCMInclRatioSpectrumStat_2040                 = (TH1D*)directoryPCMGamma_2040->Get("IncRatioStatError");
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumStat_2040    = new TGraphAsymmErrors(histoPubPCMInclRatioSpectrumStat_2040);
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumSyst_2040    = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("IncRatioSystError");
        TH1D *histoPubPCMDoubleRatioStat_2040                       = (TH1D*)directoryPCMGamma_2040->Get("DoubleRatioStatError");
        TGraphAsymmErrors* graphPubPCMDoubleRatioStat_2040          = new TGraphAsymmErrors(histoPubPCMDoubleRatioStat_2040);
        TGraphAsymmErrors* graphPubPCMDoubleRatioSyst_2040          = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("DoubleRatioSystError");



    //*******************************************************************************************************************************************
    //***************************************************** Load theroy predictions *************************************************************
    //*******************************************************************************************************************************************
    TFile* fileTheoryPbPb                           = new TFile( fileNameTheoryPbPb.Data());
    TDirectory* directoryTheoryGamma                = (TDirectory*)fileTheoryPbPb->Get("DirectPhoton");

    //______________________________________________________McGill_________________________________________
        TGraphErrors* graphTheoryMcGill0020                     = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_McGill_0020");
        TGraphErrors* graphTheoryMcGill2040                     = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_McGill_2040");
        TGraphErrors* graphTheoryPromptMcGill0020               = (TGraphErrors*) directoryTheoryGamma->Get("graphPromptPhotonYield_McGill_0020");
        TF1* fitTheoryPromptMcGill0020  = new TF1 ("fitTheoryPromptMcGill0020",
                                                "[0]/64.*1e-9*TMath::Exp(16.20-3.94*TMath::Log(x)-0.269*TMath::Log(x)**2)*1./(0.865779/TMath::Power(4,0.0694875))",
                                                0.1, 20);
        fitTheoryPromptMcGill0020->SetParameter(0,nColl0020);
        TGraphErrors* graphTheoryPromptMcGill2040               = (TGraphErrors*) directoryTheoryGamma->Get("graphPromptPhotonYield_McGill_2040");
        TF1* fitTheoryPromptMcGill2040  = new TF1 ("fitTheoryPromptMcGill2040",
                                                "[0]/64.*1e-9*TMath::Exp(16.20-3.94*TMath::Log(x)-0.269*TMath::Log(x)**2)*1./(0.865779/TMath::Power(4,0.0694875))",
                                                0.1, 20);
        fitTheoryPromptMcGill2040->SetParameter(0,nColl2040);
        TF1* fitTheoryPromptMcGill2050  = new TF1 ("fitTheoryPromptMcGill2050",
                                                "[0]/64.*1e-9*TMath::Exp(16.20-3.94*TMath::Log(x)-0.269*TMath::Log(x)**2)*1./(0.865779/TMath::Power(4,0.0694875))",
                                                0.1, 20);
        fitTheoryPromptMcGill2050->SetParameter(0,nColl2050);

    //______________________________________________________PHSD_________________________________________
        TGraphErrors* graphTheoryPHSD0020                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_PHSD_0020");
        TGraphErrors* graphTheoryPHSD2040                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_PHSD_2040");

    //______________________________________________________Chatterjee_________________________________________
        TGraphErrors* graphTheoryChatterjee0020                 = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_0020_2");
        TGraphErrors* graphTheoryChatterjee2040                 = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_2040_2");
        TGraphErrors* graphTheoryChatterjee0010                 = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_0010");
        TGraphErrors* graphTheoryChatterjee2050                 = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_2050");
        TGraphErrors* graphTheoryChatterjee4060                 = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_4060");

    //______________________________________________________new from Rapp_________________________________________
        TGraphErrors* graphTheoryRapp0010                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonRapp276GeV_0010");
        TGraphErrors* graphTheoryRapp2040                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonRapp276GeV_2040");

    //______________________________________________________VanHees_________________________________________
        TGraphErrors* graphTheoryHees0020                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_VanHees_0020");
        TGraphErrors* graphTheoryHees2040                       = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_VanHees_2040");
        TGraphErrors* graphTheoryHe0020                         = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_He_0020");
        TGraphErrors* graphTheoryHe2040                         = (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_He_2040");


    //*******************************************************************************************************************************************
    //**************************************************** Load results from PHENIX *************************************************************
    //*******************************************************************************************************************************************
    TFile* fileExperimentPbPb                       = new TFile( fileNameExperimentPbPb.Data() );
        TGraphAsymmErrors* graphPHENIXAuAuStat0010              = (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_0_20_StatErr");
        TGraphAsymmErrors* graphPHENIXAuAuSys0010               = (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_0_20_SysErr");
        TGraphAsymmErrors* graphPHENIXAuAuStat2040              = (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_20_40_StatErr");
        TGraphAsymmErrors* graphPHENIXAuAuSys2040               = (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_20_40_SysErr");



    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from PCM file *********************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMGamma                             = new TFile( inputFileNamePCM.Data());
    //________________________________________________ Load PCM 0-10% ___________________________________________________________________________
    TDirectory* directoryPCMGamma0010               = (TDirectory*)filePCMGamma->Get("Gamma_PbPb_2.76TeV_0-10%");
        TH1D* histoPCMDRPi0FitStatErr0010                       = (TH1D*) directoryPCMGamma0010->Get("DoubleRatioPi0FitStatError");
        TGraphAsymmErrors* graphPCMDRPi0FitSysErr0010           = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("DoubleRatioPi0FitSystError");
        TGraphAsymmErrors* graphPCMDRPi0FitSysAErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("DoubleRatioPi0FitSystErrorA");
        TGraphAsymmErrors* graphPCMDRPi0FitSysBErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("DoubleRatioPi0FitSystErrorB");
        TGraphAsymmErrors* graphPCMDRPi0FitSysCErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("DoubleRatioPi0FitSystErrorC");
//         TH1D* histoPCMIncRStatErr0010                        = (TH1D*) directoryPCMGamma0010->Get("IncRatioStatError");
//         TGraphAsymmErrors* graphPCMIncRSysErr0010            = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioSystError");
//         TGraphAsymmErrors* graphPCMIncRSysAErr0010           = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRSysBErr0010           = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRSysCErr0010           = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioSystErrorC");
//         TH1D* histoPCMIncRPi0FitStatErr0010                  = (TH1D*) directoryPCMGamma0010->Get("IncRatioPi0FitStatError");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysErr0010      = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioPi0FitSystError");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysAErr0010     = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioPi0FitSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysBErr0010     = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioPi0FitSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysCErr0010     = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncRatioPi0FitSystErrorC");
        TH1D* histoPCMIncGammaStatErr0010                       = (TH1D*) directoryPCMGamma0010->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPCMIncGammaSysErr0010           = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncGammaSystError");
        TGraphAsymmErrors* graphPCMIncGammaSysAErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncGammaSystErrorA");
        TGraphAsymmErrors* graphPCMIncGammaSysBErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncGammaSystErrorB");
        TGraphAsymmErrors* graphPCMIncGammaSysCErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("IncGammaSystErrorC");
//         TH1D* histoPCMPi0StatErr0010                         = (TH1D*) directoryPCMGamma0010->Get("Pi0StatError");
//         TGraphAsymmErrors* graphPCMPi0SysErr0010             = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("Pi0SystError");
//         TH1D* histoPCMPi0FitStatErr0010                      = (TH1D*) directoryPCMGamma0010->Get("Pi0FitStatError");
//         TGraphAsymmErrors* graphPCMPi0FitSysErr0010          = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("Pi0FitSystError");
        TH1D*  histoCocktailSumGamma                            = (TH1D*) directoryPCMGamma0010->Get("CocktailSumGamma");
        TH1D*  histoCocktailPi0Gamma                            = (TH1D*) directoryPCMGamma0010->Get("CocktailPi0Gamma");
        TH1D*  histoCocktailEtaGamma                            = (TH1D*) directoryPCMGamma0010->Get("CocktailEtaGamma");
        TH1D*  histoCocktailOmegaGamma                          = (TH1D*) directoryPCMGamma0010->Get("CocktailOmegaGamma");
        TH1D*  histoCocktailPhiGamma                            = (TH1D*) directoryPCMGamma0010->Get("CocktailPhiGamma");
        TH1D*  histoCocktailEtaPGamma                           = (TH1D*) directoryPCMGamma0010->Get("CocktailEtapGamma");
        TH1D*  histoCocktailRhoGamma                            = (TH1D*) directoryPCMGamma0010->Get("CocktailRhoGamma");
        TH1D*  histoCocktailSigmaGamma                          = (TH1D*) directoryPCMGamma0010->Get("CocktailSigmaGamma");

        TGraphAsymmErrors* graphPCMDirGammaStatErr0010          = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSysErr0010           = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSumErrAr0010         = NULL;
        graphPCMDirGammaStatErr0010                             = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("graphDirGammaSpectrumStat");
        graphPCMDirGammaSysErr0010                              = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("graphDirGammaSpectrumSyst");
        graphPCMDirGammaSumErrAr0010                            = (TGraphAsymmErrors*) directoryPCMGamma0010->Get("graphDirGammaSpectrumSummedAr");
    //________________________________________________ Load PCM 20-40% __________________________________________________________________________
    TDirectory* directoryPCMGamma2040               = (TDirectory*)filePCMGamma->Get("Gamma_PbPb_2.76TeV_20-40%");
//         TH1D* histoPCMDRStatErr2040                          = (TH1D*) directoryPCMGamma2040->Get("DoubleRatioStatError");
//         TGraphAsymmErrors* graphPCMDRSysErr2040              = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystError");
//         TGraphAsymmErrors* graphPCMDRSysAErr2040             = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystErrorA");
//         TGraphAsymmErrors* graphPCMDRSysBErr2040             = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystErrorB");
//         TGraphAsymmErrors* graphPCMDRSysCErr2040             = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystErrorC");
        TH1D* histoPCMDRPi0FitStatErr2040                       = (TH1D*) directoryPCMGamma2040->Get("DoubleRatioPi0FitStatError");
        TGraphAsymmErrors* graphPCMDRPi0FitSysErr2040           = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystError");
        TGraphAsymmErrors* graphPCMDRPi0FitSysAErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystErrorA");
        TGraphAsymmErrors* graphPCMDRPi0FitSysBErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystErrorB");
        TGraphAsymmErrors* graphPCMDRPi0FitSysCErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystErrorC");
//         TH1D* histoPCMIncRStatErr2040                        = (TH1D*) directoryPCMGamma2040->Get("IncRatioStatError");
//         TGraphAsymmErrors* graphPCMIncRSysErr2040            = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystError");
//         TGraphAsymmErrors* graphPCMIncRSysAErr2040           = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRSysBErr2040           = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRSysCErr2040           = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystErrorC");
//         TH1D* histoPCMIncRPi0FitStatErr2040                  = (TH1D*) directoryPCMGamma2040->Get("IncRatioPi0FitStatError");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysErr2040      = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystError");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysAErr2040     = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysBErr2040     = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysCErr2040     = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystErrorC");
        TH1D* histoPCMIncGammaStatErr2040                       = (TH1D*) directoryPCMGamma2040->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPCMIncGammaSysErr2040           = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystError");
        TGraphAsymmErrors* graphPCMIncGammaSysAErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystErrorA");
        TGraphAsymmErrors* graphPCMIncGammaSysBErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystErrorB");
        TGraphAsymmErrors* graphPCMIncGammaSysCErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystErrorC");
//         TH1D* histoPCMPi0StatErr2040                         = (TH1D*) directoryPCMGamma2040->Get("Pi0StatError");
//         TGraphAsymmErrors* graphPCMPi0SysErr2040             = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("Pi0SystError");
//         TH1D* histoPCMPi0FitStatErr2040                      = (TH1D*) directoryPCMGamma2040->Get("Pi0FitStatError");
//         TGraphAsymmErrors* graphPCMPi0FitSysErr2040          = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("Pi0FitSystError");
//
        TGraphAsymmErrors* graphPCMDirGammaStatErr2040          = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSysErr2040           = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSumErrAr2040         = NULL;
        graphPCMDirGammaStatErr2040                             = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("graphDirGammaSpectrumStat");
        graphPCMDirGammaSysErr2040                              = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("graphDirGammaSpectrumSyst");
        graphPCMDirGammaSumErrAr2040                            = (TGraphAsymmErrors*) directoryPCMGamma2040->Get("graphDirGammaSpectrumSummedAr");
    //________________________________________________ Load PCM 20-50% __________________________________________________________________________
    TDirectory* directoryPCMGamma2050               = (TDirectory*)filePCMGamma->Get("Gamma_PbPb_2.76TeV_20-50%");
//         TH1D* histoPCMDRStatErr2050                          = (TH1D*) directoryPCMGamma2050->Get("DoubleRatioStatError");
//         TGraphAsymmErrors* graphPCMDRSysErr2050              = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioSystError");
//         TGraphAsymmErrors* graphPCMDRSysAErr2050             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioSystErrorA");
//         TGraphAsymmErrors* graphPCMDRSysBErr2050             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioSystErrorB");
//         TGraphAsymmErrors* graphPCMDRSysCErr2050             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioSystErrorC");
        TH1D* histoPCMDRPi0FitStatErr2050                       = (TH1D*) directoryPCMGamma2050->Get("DoubleRatioPi0FitStatError");
        TGraphAsymmErrors* graphPCMDRPi0FitSysErr2050           = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioPi0FitSystError");
        TGraphAsymmErrors* graphPCMDRPi0FitSysAErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioPi0FitSystErrorA");
        TGraphAsymmErrors* graphPCMDRPi0FitSysBErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioPi0FitSystErrorB");
        TGraphAsymmErrors* graphPCMDRPi0FitSysCErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("DoubleRatioPi0FitSystErrorC");
//         TH1D* histoPCMIncRStatErr2050                        = (TH1D*) directoryPCMGamma2050->Get("IncRatioStatError");
//         TGraphAsymmErrors* graphPCMIncRSysErr2050            = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioSystError");
//         TGraphAsymmErrors* graphPCMIncRSysAErr2050           = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRSysBErr2050           = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRSysCErr2050           = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioSystErrorC");
//         TH1D* histoPCMIncRPi0FitStatErr2050                  = (TH1D*) directoryPCMGamma2050->Get("IncRatioPi0FitStatError");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysErr2050      = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioPi0FitSystError");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysAErr2050     = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioPi0FitSystErrorA");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysBErr2050     = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioPi0FitSystErrorB");
//         TGraphAsymmErrors* graphPCMIncRFitPi0SysCErr2050     = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncRatioPi0FitSystErrorC");
        TH1D* histoPCMIncGammaStatErr2050                       = (TH1D*) directoryPCMGamma2050->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPCMIncGammaSysErr2050           = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncGammaSystError");
        TGraphAsymmErrors* graphPCMIncGammaSysAErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncGammaSystErrorA");
        TGraphAsymmErrors* graphPCMIncGammaSysBErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncGammaSystErrorB");
        TGraphAsymmErrors* graphPCMIncGammaSysCErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("IncGammaSystErrorC");
//         TH1D* histoPCMPi0StatErr2050                         = (TH1D*) directoryPCMGamma2050->Get("Pi0StatError");
//         TGraphAsymmErrors* graphPCMPi0SysErr2050             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("Pi0SystError");
//         TH1D* histoPCMPi0FitStatErr2050                      = (TH1D*) directoryPCMGamma2050->Get("Pi0FitStatError");
//         TGraphAsymmErrors* graphPCMPi0FitSysErr2050          = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("Pi0FitSystError");
//
        TGraphAsymmErrors* graphPCMDirGammaStatErr2050          = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSysErr2050           = NULL;
        TGraphAsymmErrors* graphPCMDirGammaSumErrAr2050         = NULL;
        graphPCMDirGammaStatErr2050                             = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("graphDirGammaSpectrumStat");
        graphPCMDirGammaSysErr2050                              = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("graphDirGammaSpectrumSyst");
        graphPCMDirGammaSumErrAr2050                            = (TGraphAsymmErrors*) directoryPCMGamma2050->Get("graphDirGammaSpectrumSummedAr");

    TDirectory* directoryTheory                     = (TDirectory*)filePCMGamma->Get("Theory");
        TGraphAsymmErrors* graphTheoryNLODR0010                 = (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatio_0-10%");
        TGraphAsymmErrors* graphTheoryNLODR2040                 = (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatio_20-40%");
        TGraphAsymmErrors* graphTheoryNLODR2050                 = (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatio_20-50%");
        TGraphAsymmErrors* graphTheoryNLO0010                   = (TGraphAsymmErrors*) directoryTheory->Get("NLO_0-10%");
        TGraphAsymmErrors* graphTheoryNLO2040                   = (TGraphAsymmErrors*) directoryTheory->Get("NLO_20-40%");
        TGraphAsymmErrors* graphTheoryNLO2050                   = (TGraphAsymmErrors*) directoryTheory->Get("NLO_20-50%");
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

        TF1* fitTheoryGammaEPSO90010                            = FitObject("l","fitTheoryGammaEPSO90010","Photon",graphTheoryEPS090010,1.6,12.5,NULL,"QNRMEX0+");
        cout << WriteParameterToFile(fitTheoryGammaEPSO90010)<< endl;
        fitTheoryGammaEPSO90010->SetRange(0.9,14);
        TF1* fitTheoryGammaEPSO92040                            = FitObject("l","fitTheoryGammaEPSO90010","Photon",graphTheoryEPS092040,1.6,12.5,NULL,"QNRMEX0+");
        cout << WriteParameterToFile(fitTheoryGammaEPSO92040)<< endl;
        fitTheoryGammaEPSO92040->SetRange(0.9,14);
        TF1* fitTheoryGammaEPSO92050                            = FitObject("l","fitTheoryGammaEPSO90010","Photon",graphTheoryEPS092050,1.6,12.5,NULL,"QNRMEX0+");
        cout << WriteParameterToFile(fitTheoryGammaEPSO92050)<< endl;
        fitTheoryGammaEPSO92050->SetRange(0.9,14);


    //*******************************************************************************************************************************************
    //*********************************************** Load PHOS histograms from PHOS file *******************************************************
    //*******************************************************************************************************************************************
    TFile* filePHOSGamma                            = new TFile( inputFileNamePHOS.Data());
    //________________________________________________ Load PHOS 0-10% __________________________________________________________________________
    TDirectory* directoryPHOSGamma0010              = (TDirectory*)filePHOSGamma->Get("PHOS_PbPb_2760_Centrality_00-10");
        TH1D* histoPHOSDRPi0FitStatErr0010                      = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_Stat");
        TH1D* histoPHOSDRPi0FitSysErr0010                       = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_Syst");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysErr0010          = new TGraphAsymmErrors(histoPHOSDRPi0FitSysErr0010);
        TH1D* histoPHOSDRPi0FitSysAErr0010                      = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_SystA");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysAErr0010         = new TGraphAsymmErrors(histoPHOSDRPi0FitSysAErr0010);
        TH1D* histoPHOSDRPi0FitSysBcontErr0010                  = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_SystBcont");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBcontErr0010     = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcontErr0010);
//         cout << "B Syst. Cont" << endl;
//         graphPHOSDRPi0FitSysBcontErr0010->Print();
        TH1D* histoPHOSDRPi0FitSysBpi0Err0010                   = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_SystBpi0");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBpi0Err0010      = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBpi0Err0010);
//         cout << "B Sys pi0" << endl;
//         graphPHOSDRPi0FitSysBpi0Err0010->Print();
        TH1D* histoPHOSDRPi0FitSysBcocktErr0010                 = (TH1D*) directoryPHOSGamma0010->Get("hPHOS_DoubleRatio_PbPb_cen00-10_SystBcockt");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysBcocktErr0010    = new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcocktErr0010);
//         cout << "B Sys cockt" << endl;
//         graphPHOSDRPi0FitSysBcocktErr0010->Print();
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

    TGraphAsymmErrors* graphPCMIncGammaStatErr0010              = new TGraphAsymmErrors(histoPCMIncGammaStatErr0010);
        while(graphPCMIncGammaStatErr0010->GetX()[0]<=1.) graphPCMIncGammaStatErr0010->RemovePoint(0);
    TGraphAsymmErrors* graphPCMIncGammaStatSysAErr0010          = AddErrorsOfGraphsQuadratically (graphPCMIncGammaStatErr0010, graphPCMIncGammaSysAErr0010);

    TGraphAsymmErrors* graphPHOSIncGammaStatErr0010             = new TGraphAsymmErrors(histoPHOSIncGammaStatErr0010);
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

    TGraphAsymmErrors* graphPCMIncGammaStatErr2040              = new TGraphAsymmErrors(histoPCMIncGammaStatErr2040);
        while(graphPCMIncGammaStatErr2040->GetX()[0]<=1.) graphPCMIncGammaStatErr2040->RemovePoint(0);
    TGraphAsymmErrors* graphPCMIncGammaStatSysAErr2040          = AddErrorsOfGraphsQuadratically (graphPCMIncGammaStatErr2040, graphPCMIncGammaSysAErr2040);

    TGraphAsymmErrors* graphPHOSIncGammaStatErr2040             = new TGraphAsymmErrors(histoPHOSIncGammaStatErr2040);
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

    TGraphAsymmErrors* graphPCMIncGammaStatErr2050              = new TGraphAsymmErrors(histoPCMIncGammaStatErr2050);
        while(graphPCMIncGammaStatErr2050->GetX()[0]<=1.) graphPCMIncGammaStatErr2050->RemovePoint(0);
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
    TF1* fitFullDirGamma0010Sys                                 = FitObject("qcd","fitFullDirGamma0010Sys","Photon",graphCombDirGammaSpectrumSystErr0010,0.9,14,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitFullDirGamma0010Sys)<< endl;
    TF1* fitFullDirGamma0010Stat                                = FitObject("qcd","fitFullDirGamma0010Stat","Photon",graphCombDirGammaSpectrumStatErr0010,0.9,14,NULL,"QNRMEX0+");
    fileFinalResults << WriteParameterToFile(fitFullDirGamma0010Stat)<< endl;

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
        graphCombThermalGammaSpectrumStatErr0010                    = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumStatErr0010);
        graphCombThermalGammaSpectrumStatErr0010->Print();
        fitPureThermalGamma0010Stat                                 = FitObject("e","fitPureThermalGamma0010Stat","Photon",graphCombThermalGammaSpectrumStatErr0010,0.9,2.1,NULL,"QNRMEX0+");
        fileFinalResults << WriteParameterToFile(fitPureThermalGamma0010Stat)<< endl;
        fitFullThermalGamma0010Stat                                 = FitObject("qcd","fitFullDirGamma0010Stat","Photon",graphCombThermalGammaSpectrumStatErr0010,0.9,14,NULL,"QNRMEX0+");
        fileFinalResults << WriteParameterToFile(fitFullThermalGamma0010Stat)<< endl;
        graphRatioCombFitThermalGammaStatErr0010                    = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumStatErr0010->Clone();
        graphRatioCombFitThermalGammaStatErr0010                    = CalculateGraphErrRatioToFit(graphRatioCombFitThermalGammaStatErr0010, fitFullThermalGamma0010Stat);
    }
    if (graphCombDirGammaSpectrumSystErr0010){
        graphCombThermalGammaSpectrumSysErr0010                     = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSystErr0010);
        graphCombThermalGammaSpectrumSysErr0010->Print();
        fitFullThermalGamma0010Sys                                  = FitObject("qcd","fitFullDirGamma0010Sys","Photon",graphCombThermalGammaSpectrumSysErr0010,0.9,14,NULL,"QNRMEX0+");
        fileFinalResults << WriteParameterToFile(fitFullThermalGamma0010Sys)<< endl;
        graphRatioCombFitThermalGammaSysErr0010                     = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumSysErr0010->Clone();
        graphRatioCombFitThermalGammaSysErr0010                     = CalculateGraphErrRatioToFit(graphRatioCombFitThermalGammaSysErr0010, fitFullThermalGamma0010Sys);

    }
    if (graphCombDirGammaSpectrumSumErr0010Ar){
        graphCombThermalGammaSpectrumSumErr0010Ar                   = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSumErr0010Ar, newBinsComb, 20);
    }

    // Calculate RAA
    cout << __LINE__ << endl;
    cout << endl << "Calculating RAA" << endl;
    TGraphAsymmErrors* graphCombRAADirGammaStat0010                 = NULL;
    TGraphAsymmErrors* graphCombRAADirGammaSys0010                  = NULL;
    TGraphAsymmErrors* graphCombRAADirGammaSum0010Ar                = NULL;
    if (graphCombDirGammaSpectrumStatErr0010)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumStatErr0010, &graphCombRAADirGammaStat0010);
    if (graphCombDirGammaSpectrumSystErr0010)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSystErr0010, &graphCombRAADirGammaSys0010);
    if (graphCombDirGammaSpectrumSumErr0010Ar)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSumErr0010Ar, &graphCombRAADirGammaSum0010Ar, newBinsComb, 20);
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
    if (graphCombDirGammaSpectrumStatErr2050)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2050, graphCombDirGammaSpectrumStatErr2050, &graphCombRAADirGammaStat2050);
    if (graphCombDirGammaSpectrumSystErr2050)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2050, graphCombDirGammaSpectrumSystErr2050, &graphCombRAADirGammaSys2050);
    if (graphCombDirGammaSpectrumSumErr2050Ar)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill2050, graphCombDirGammaSpectrumSumErr2050Ar, &graphCombRAADirGammaSum2050Ar, newBinsComb, 20);
    //if (graphCombRAADirGammaStat2050) graphCombRAADirGammaStat2050->Print();
    //if (graphCombRAADirGammaSys2050) graphCombRAADirGammaSys2050->Print();
    //if (graphCombRAADirGammaSum2050Ar) graphCombRAADirGammaSum2050Ar->Print();


    cout << "Plotting DR individual measurements at " << __LINE__ << endl;
    //*******************************************************************************************************************************************
    //******************************************* DR plot with individual measurements **********************************************************
    //*******************************************************************************************************************************************
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

    //_______________________________________________________________ 20-50% dummy lower panel __________________________________________________
    TH2D *dummyDR3 ;
    dummyDR3 = new TH2D("dummyDR3", "dummyDR3", 1000, 0., 22, 1000., doubleRatio[0], doubleRatio[1]);
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


    cout << "Plotting PCM DR with NLO at " << __LINE__ << endl;
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
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr0010, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR0010, 3, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        graphTheoryNLODR0010->Draw("p3lsame");
        graphPCMDRPi0FitSysErr0010->Draw("E2same");
        histoPCMDRPi0FitStatErr0010->Draw("p,same,e0,X0");

        TLegend* legendDRPCMNLO0010 = new TLegend(0.15,0.82-1.1*0.85*textsizeLabelsPad1*3,0.5,0.82);
        legendDRPCMNLO0010->SetFillStyle(0);
        legendDRPCMNLO0010->SetFillColor(0);
        legendDRPCMNLO0010->SetLineColor(0);
        legendDRPCMNLO0010->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRPCMNLO0010->SetMargin(0.2);
        legendDRPCMNLO0010->SetTextFont(42);
        legendDRPCMNLO0010->AddEntry(graphPCMDRPi0FitSysErr0010,"PCM","pf");
        legendDRPCMNLO0010->AddEntry(graphTheoryNLODR0010,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
        legendDRPCMNLO0010->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
        legendDRPCMNLO0010->Draw();

        labelDRCent0010->Draw();

    //_______________________________________________________________ 20-40% panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR2040, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        graphTheoryNLODR2040->Draw("p3lsame");
        graphPCMDRPi0FitSysErr2040->Draw("E2same");
        histoPCMDRPi0FitStatErr2040->Draw("p,same,e0,X0");

        TLegend* legendDRPCMNLO2040 = new TLegend(0.15,0.85-1.1*0.85*textsizeLabelsPad2*3,0.5,0.85);
        legendDRPCMNLO2040->SetFillStyle(0);
        legendDRPCMNLO2040->SetFillColor(0);
        legendDRPCMNLO2040->SetLineColor(0);
        legendDRPCMNLO2040->SetTextSize(0.85*textsizeLabelsPad2);
        legendDRPCMNLO2040->SetMargin(0.2);
        legendDRPCMNLO2040->SetTextFont(42);
        legendDRPCMNLO2040->AddEntry(graphPCMDRPi0FitSysErr2040,"PCM","pf");
        legendDRPCMNLO2040->AddEntry(graphTheoryNLODR2040,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
        legendDRPCMNLO2040->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
        legendDRPCMNLO2040->Draw();

        labelDRCent2040->Draw();

    //_______________________________________________________________ 20-50% panel _______________________________________________________________
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErr2050, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR2050, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        graphTheoryNLODR2050->Draw("p3lsame");
        graphPCMDRPi0FitSysErr2050->Draw("E2same");
        histoPCMDRPi0FitStatErr2050->Draw("p,same,e0,X0");

        TLegend* legendDRPCMNLO2050 = new TLegend(0.15,0.87-1.1*0.85*textsizeLabelsPad3*3,0.5,0.87);
        legendDRPCMNLO2050->SetFillStyle(0);
        legendDRPCMNLO2050->SetFillColor(0);
        legendDRPCMNLO2050->SetLineColor(0);
        legendDRPCMNLO2050->SetTextSize(0.85*textsizeLabelsPad3);
        legendDRPCMNLO2050->SetMargin(0.2);
        legendDRPCMNLO2050->SetTextFont(42);
        legendDRPCMNLO2050->AddEntry(graphPCMDRPi0FitSysErr2050,"PCM","pf");
        legendDRPCMNLO2050->AddEntry(graphTheoryNLODR2050,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
        legendDRPCMNLO2050->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
        legendDRPCMNLO2050->Draw();

        labelDRCent2050->Draw();

    canvasRatioIndDR->SaveAs(Form("%s/DR_PCMMeasurementTheory.%s", outputDir.Data(), suffix.Data()));


    //*******************************************************************************************************************************************
    //********************************************** DR plot with PCM measurements all cents in one plot  ***************************************
    //*******************************************************************************************************************************************
    TCanvas *canvasDoubleRatio = GetAndSetCanvas("canvasDoubleRatioFinal");
    canvasDoubleRatio->SetLogx();

    TH2D *dummyDR ;
    dummyDR = new TH2D("dummyDR", "dummyDR", 1000, 0., 16, 1000., doubleRatio[0], doubleRatio[1]);
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


        TLegend* legendDRAllCentPCM = GetAndSetLegend(0.15,0.75,4);
        legendDRAllCentPCM->AddEntry(graphPCMDRPi0FitSysErr0010,"  0-10% PCM","pf");
        legendDRAllCentPCM->AddEntry(graphPCMDRPi0FitSysErr2040,"20-40% PCM","pf");
        legendDRAllCentPCM->AddEntry(graphPCMDRPi0FitSysErr2050,"20-50% PCM","pf");

        legendDRAllCentPCM->Draw();

    canvasDoubleRatio->Print(Form("%s/DR_PCMMeasurementAllCentsInOne.%s", outputDir.Data(), suffix.Data()));


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



    cout << "Plotting DR combined with theory at " << __LINE__ << endl;
    //*******************************************************************************************************************************************
    //********************************************** DR plot with combined measurement **********************************************************
    //*******************************************************************************************************************************************
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

        TBox* boxSysCErrPCM0010 = CreateBoxConv(colorPCMBox, 0.76, 1.-(SysCPCMIncGamma0010) , 0.815, 1.+(SysCPCMIncGamma0010));
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


        TBox* boxSysCErrPCM2040 = CreateBoxConv(colorPCMBox, 0.76, 1.-(SysCPCMIncGamma2040) , 0.815, 1.+(SysCPCMIncGamma2040));
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
        TBox* boxSysCErrPCM2050 = CreateBoxConv(colorPCMBox, 0.76, 1.-(SysCPCMIncGamma2050) , 0.815, 1.+(SysCPCMIncGamma2050));
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

    TGraphAsymmErrors* graphCombIncGammaSysErr0010Plot = ScaleGraph(graphPCMIncGammaSysAErr0010,100);
    TGraphAsymmErrors* graphCombIncGammaStatErr0010Plot = ScaleGraph(graphPCMIncGammaStatErr0010,100);

    TGraphAsymmErrors* graphCombIncGammaSysErr2040Plot = ScaleGraph(graphPCMIncGammaSysAErr2040,10);
    TGraphAsymmErrors* graphCombIncGammaStatErr2040Plot = ScaleGraph(graphPCMIncGammaStatErr2040,10);

    TGraphAsymmErrors* graphCombIncGammaSysErr2050Plot = ScaleGraph(graphPCMIncGammaSysAErr2050,1);
    TGraphAsymmErrors* graphCombIncGammaStatErr2050Plot = ScaleGraph(graphPCMIncGammaStatErr2050,1);

    TH1D* histoFitDummyPlotting = (TH1D*) histoFitQCDIncGammaComb2040->Clone("histoFitDummyPlotting");

    TCanvas *canvasIncGamma = new TCanvas("canvasIncGamma","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasIncGamma, 0.16, 0.01, 0.01, 0.07);
    canvasIncGamma->SetLogy();
    canvasIncGamma->SetLogx();

    Int_t textSizeLabelsPixelIncGam = 48;
    Double_t textsizeLabelsIncGamma = 0;
    if (canvasIncGamma->XtoPixel(canvasIncGamma->GetX2()) < canvasIncGamma->YtoPixel(canvasIncGamma->GetY1())){
        textsizeLabelsIncGamma = (Double_t)textSizeLabelsPixelIncGam/canvasIncGamma->XtoPixel(canvasIncGamma->GetX2()) ;
    } else {
        textsizeLabelsIncGamma = (Double_t)textSizeLabelsPixelIncGam/canvasIncGamma->YtoPixel(canvasIncGamma->GetY1());
    }

    TH2D *dummyGamma ;
    dummyGamma = new TH2D("dummyGamma", "dummyGamma", 1000, 0., 22, 1000., 6e-8,8e3);
    SetStyleHistoTH2ForGraphs( dummyGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{inc}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                            0.85*textsizeLabelsIncGamma, textsizeLabelsIncGamma, 0.85*textsizeLabelsIncGamma, textsizeLabelsIncGamma, 0.75, 1.6);
    dummyGamma->GetXaxis()->SetLabelOffset(-0.015);
    dummyGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
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

    TLatex *labelScalingIncGamma0010 = new TLatex(11.25,1.8E-3,"x 10^{2}");
    SetStyleTLatex( labelScalingIncGamma0010, 0.85*textsizeLabelsIncGamma,4,colorComb0010,42,kFALSE);
    labelScalingIncGamma0010->Draw();

    TLatex *labelScalingIncGamma2040 = new TLatex(11.25,1E-4,"x 10^{1}");
    SetStyleTLatex( labelScalingIncGamma2040, 0.85*textsizeLabelsIncGamma,4,colorComb2040,42,kFALSE);
    labelScalingIncGamma2040->Draw();

    TLatex *labelScalingIncGamma2050 = new TLatex(11.25,2E-6,"x 10^{0}");
    SetStyleTLatex( labelScalingIncGamma2050, 0.85*textsizeLabelsIncGamma,4,colorComb2050,42,kFALSE);
    labelScalingIncGamma2050->Draw();

    TLatex *labelIncGammaColl = new TLatex(0.6,0.91,collisionSystem.Data());
    SetStyleTLatex( labelIncGammaColl, 0.85*textsizeLabelsIncGamma,4);
    labelIncGammaColl->Draw();
    SetStyleHisto(histoFitDummyPlotting, widthCommonFit, 5, kGray+1);

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
    legendIncGamma->AddEntry(histoFitDummyPlotting,"Fits to #gamma_{inc}","l");
    legendIncGamma->Draw();

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
    TGraphAsymmErrors* graphTheoryNLO0010Plot                                                = ScaleGraph(graphTheoryNLO0010,100);
    TGraphAsymmErrors* graphTheoryEPS090010Plot                                              = ScaleGraph(graphTheoryEPS090010,100);
    TGraphAsymmErrors* graphTheoryCT100010Plot                                               = ScaleGraph(graphTheoryCT100010,100);
    TGraphErrors* graphTheoryPromptMcGill0020Plot                                            = ScaleGraph(graphTheoryPromptMcGill0020,100);
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
    TGraphAsymmErrors* graphTheoryNLO2040Plot                                                = ScaleGraph(graphTheoryNLO2040,10);
    TGraphAsymmErrors* graphTheoryEPS092040Plot                                              = ScaleGraph(graphTheoryEPS092040,10);
    TGraphAsymmErrors* graphTheoryCT102040Plot                                               = ScaleGraph(graphTheoryCT102040,10);
    TGraphErrors* graphTheoryPromptMcGill2040Plot                                            = ScaleGraph(graphTheoryPromptMcGill2040,10);
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
    TGraphAsymmErrors* graphTheoryNLO2050Plot                                                = ScaleGraph(graphTheoryNLO2050,1);
    TGraphAsymmErrors* graphTheoryEPS092050Plot                                              = ScaleGraph(graphTheoryEPS092050,1);
    TGraphAsymmErrors* graphTheoryCT102050Plot                                               = ScaleGraph(graphTheoryCT102050,1);


    TCanvas *canvasDirGamma = new TCanvas("canvasDirGamma","",10,10,1200,1400);  // gives the page size
    DrawGammaCanvasSettings( canvasDirGamma, 0.165, 0.01, 0.01, 0.07);
    canvasDirGamma->SetLogy();
    canvasDirGamma->SetLogx();

    Int_t textSizeLabelsPixelDirGam = 48;
    Double_t textsizeLabelsDirGamma = 0;
    if (canvasDirGamma->XtoPixel(canvasDirGamma->GetX2()) < canvasDirGamma->YtoPixel(canvasDirGamma->GetY1())){
        textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGamma->XtoPixel(canvasDirGamma->GetX2()) ;
    } else {
        textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGamma->YtoPixel(canvasDirGamma->GetY1());
    }

        TH2D *dummyDirGamma ;
        dummyDirGamma = new TH2D("dummyDirGamma", "dummyDirGamma", 1000, 0., 22, 1000., 1.2e-8,1.5e5);
        SetStyleHistoTH2ForGraphs( dummyDirGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                                0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
        dummyDirGamma->GetXaxis()->SetLabelOffset(-0.015);
        dummyDirGamma->GetXaxis()->SetTickLength(0.025);
        dummyDirGamma->GetYaxis()->SetTickLength(0.025);
        dummyDirGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    //     dummyDirGamma->GetXaxis()->SetRangeUser(0,16);
        dummyDirGamma->DrawCopy();

        TLatex *labelScalingDirGamma0010 = new TLatex(11.2,1.3E-3,"x 10^{2}");
        SetStyleTLatex( labelScalingDirGamma0010, 0.85*textsizeLabelsDirGamma,4,colorComb0010,42,kFALSE);
        labelScalingDirGamma0010->Draw();

        TLatex *labelScalingDirGamma2040 = new TLatex(11.2,6.0E-5,"x 10^{1}");
        SetStyleTLatex( labelScalingDirGamma2040, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
        labelScalingDirGamma2040->Draw();

        TLatex *labelScalingDirGamma2050 = new TLatex(11.2,7.5E-7,"x 10^{0}");
        SetStyleTLatex( labelScalingDirGamma2050, 0.85*textsizeLabelsDirGamma,4,colorComb2050,42,kFALSE);
        labelScalingDirGamma2050->Draw();

        TLatex *labelDirGammaColl = new TLatex(0.25,0.94,Form("%s",collisionSystem.Data()));
        SetStyleTLatex( labelDirGammaColl, 0.85*textsizeLabelsDirGamma,4);
        labelDirGammaColl->Draw();
        SetStyleHisto(histoFitDummyPlotting, widthCommonFit, 5, kGray+1);
    //
        TLegend* legendDirGamma = new TLegend(0.24,0.93-1.*0.85*textsizeLabelsDirGamma*3,0.24+0.21,0.93);
        legendDirGamma->SetFillStyle(0);
        legendDirGamma->SetFillColor(0);
        legendDirGamma->SetLineColor(0);
        legendDirGamma->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGamma->SetMargin(0.25);
        legendDirGamma->SetTextFont(42);
        legendDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot,"  0-10% ALICE","pf");
        legendDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
        legendDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr2050Plot,"20-50% ALICE","pf");
        legendDirGamma->Draw();

        if (graphCombDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphCombDirGammaSpectrumSystErr0010Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphCombDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
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

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum.%s",outputDir.Data(),suffix.Data()));


    cout << "Plotting direct gamma spectra with models at " << __LINE__ << endl;
    SetStyleGammaNLOTGraphWithBand( graphTheoryEPS090010Plot, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
    graphTheoryEPS090010Plot->Draw("p3lsame");
    SetStyleGammaNLOTGraphWithBand( graphTheoryCT100010Plot, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
    graphTheoryCT100010Plot->Draw("p3lsame");

    SetStyleGammaNLOTGraphWithBand( graphTheoryNLO0010Plot, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
    graphTheoryNLO0010Plot->Draw("p3lsame");
    SetStyleGammaNLOTGraphWithBand( graphTheoryPromptMcGill0020Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
    graphTheoryPromptMcGill0020Plot->Draw("lsame");

    if (graphCombDirGammaSpectrumSystErr0010Plot){
        DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
        graphCombDirGammaSpectrumSystErr0010Plot->Draw("E2same");
    }
    if (graphCombDirGammaSpectrumStatErr0010Plot){
        DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        graphCombDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
    }

    SetStyleGammaNLOTGraphWithBand( graphTheoryEPS092040Plot, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
    graphTheoryEPS092040Plot->Draw("p3lsame");

    SetStyleGammaNLOTGraphWithBand( graphTheoryCT102040Plot, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
    graphTheoryCT102040Plot->Draw("p3lsame");

    SetStyleGammaNLOTGraphWithBand( graphTheoryNLO2040Plot, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);
        graphTheoryNLO2040Plot->Draw("p3lsame");
        SetStyleGammaNLOTGraphWithBand( graphTheoryPromptMcGill2040Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
        graphTheoryPromptMcGill2040Plot->Draw("lsame");

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

        SetStyleGammaNLOTGraphWithBand( graphTheoryEPS092050Plot, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
        graphTheoryEPS092050Plot->Draw("p3lsame");

        SetStyleGammaNLOTGraphWithBand( graphTheoryCT102050Plot, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
        graphTheoryCT102050Plot->Draw("p3lsame");

        SetStyleGammaNLOTGraphWithBand( graphTheoryNLO2050Plot, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);
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

        TLegend* legendDirGammaTheory = new TLegend(0.525,0.93-1.*0.85*textsizeLabelsDirGamma*3,0.525+0.21,0.93);
        legendDirGammaTheory->SetFillStyle(0);
        legendDirGammaTheory->SetFillColor(0);
        legendDirGammaTheory->SetLineColor(0);
        legendDirGammaTheory->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaTheory->SetMargin(0.25);
        legendDirGammaTheory->SetTextFont(42);
        legendDirGammaTheory->AddEntry(graphTheoryNLO0010Plot,"PDF: CTEQ6M5, FF: GRV ","l");
        legendDirGammaTheory->AddEntry(graphTheoryPromptMcGill2040Plot,"(n)PDF: CTEQ6.1M/EPS09,","l");
        legendDirGammaTheory->AddEntry((TObject*)0,"FF: BFG2","");
        legendDirGammaTheory->Draw();

        TLegend* legendDirGammaTheory2 = new TLegend(0.525,0.93-1.*0.85*textsizeLabelsDirGamma*7-0.02,0.525+0.21,0.93-1.*0.85*textsizeLabelsDirGamma*3-0.02);
        legendDirGammaTheory2->SetFillStyle(0);
        legendDirGammaTheory2->SetFillColor(0);
        legendDirGammaTheory2->SetLineColor(0);
        legendDirGammaTheory2->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaTheory2->SetMargin(0.25);
        legendDirGammaTheory2->SetTextFont(42);
        legendDirGammaTheory2->AddEntry((TObject*)0,"#it{JETPHOX}","");
        legendDirGammaTheory2->AddEntry(graphTheoryCT100010Plot,"PDF: CT10, FF: BFG2","f");
        legendDirGammaTheory2->AddEntry(graphTheoryEPS090010Plot,"nPDF: EPS09, FF: BFG2","f");
        legendDirGammaTheory2->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
        legendDirGammaTheory2->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum_withNLO.%s",outputDir.Data(),suffix.Data()));

    cout << __LINE__ << endl;
    //************************************************************************************************************************************
    TGraphErrors* graphTheoryMcGill0020Plot    = ScaleGraph(graphTheoryMcGill0020,100);
    TGraph* graphTheoryPHSD0020Plot            = ScaleGraph(graphTheoryPHSD0020,100);
    graphTheoryPHSD0020Plot->RemovePoint(0);
    TGraph* graphTheoryChatterjee0010Plot      = ScaleGraph(graphTheoryChatterjee0010,100);
    TGraph* graphTheoryRapp0010Plot            = ScaleGraph(graphTheoryRapp0010,100);
    TGraph* graphTheoryHees0020Plot            = ScaleGraph(graphTheoryHees0020,100);
    TGraph* graphTheoryHe0020Plot              = ScaleGraph(graphTheoryHe0020,100);

    TGraphErrors* graphTheoryMcGill2040Plot    = ScaleGraph(graphTheoryMcGill2040,10);
    TGraph* graphTheoryPHSD2040Plot            = ScaleGraph(graphTheoryPHSD2040,10);
    graphTheoryPHSD2040Plot->RemovePoint(0);
    TGraph* graphTheoryChatterjee2040Plot      = ScaleGraph(graphTheoryChatterjee2040,10);
    TGraph* graphTheoryRapp2040Plot            = ScaleGraph(graphTheoryRapp2040,10);
    TGraph* graphTheoryHees2040Plot            = ScaleGraph(graphTheoryHees2040,10);
    TGraph* graphTheoryHe2040Plot              = ScaleGraph(graphTheoryHe2040,10);

    TH2D *dummyDirGammaTheory = new TH2D("dummyDirGammaTheory", "dummyDirGamma", 1000, 0., 22, 1000., 4e-8,9e2);
    SetStyleHistoTH2ForGraphs( dummyDirGammaTheory, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})", 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummyDirGammaTheory->GetXaxis()->SetLabelOffset(-0.015);
    dummyDirGammaTheory->GetXaxis()->SetTickLength(0.025);
    dummyDirGammaTheory->GetYaxis()->SetTickLength(0.025);
    dummyDirGammaTheory->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
    dummyDirGammaTheory->DrawCopy();

        SetStyleGammaNLOTGraphWithBand( graphTheoryMcGill0020Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryPHSD0020Plot, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryChatterjee0010Plot, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryHees0020Plot, 3.0, styleHees, colorHees, 3015, colorHees, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryHe0020Plot, 3.0, styleHe, colorHe, 3015, colorHe, 0);
        if (graphCombDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
        }
        SetStyleGammaNLOTGraphWithBand( graphTheoryMcGill2040Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryPHSD2040Plot, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryChatterjee2040Plot, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryHees2040Plot, 3.0, styleHees, colorHees, 3015, colorHees, 0);
        SetStyleGammaNLOTGraphWithBand( graphTheoryHe2040Plot, 3.0, styleHe, colorHe, 3015, colorHe, 0);
        if (graphCombDirGammaSpectrumSystErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
        }
        if (graphCombDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
        }
        if (graphCombDirGammaSpectrumSumErr2040ArPlot){
            DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , styleNLOMcGill, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
        }

        labelScalingDirGamma0010->Draw();
        labelScalingDirGamma2040->Draw();
        labelScalingDirGamma2050->Draw();
        TLatex *labelDirGammaCollRedX = new TLatex(0.585,0.93,Form("%s",collisionSystem.Data()));
        SetStyleTLatex( labelDirGammaCollRedX, 0.85*textsizeLabelsDirGamma,4);
        labelDirGammaCollRedX->Draw();


        TLegend* legendDirGammaRedX = new TLegend(0.58,0.92-1.1*0.85*textsizeLabelsDirGamma*3,0.585+0.21,0.92);
        legendDirGammaRedX->SetFillStyle(0);
        legendDirGammaRedX->SetFillColor(0);
        legendDirGammaRedX->SetLineColor(0);
        legendDirGammaRedX->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaRedX->SetMargin(0.3);
        legendDirGammaRedX->SetTextFont(42);
        legendDirGammaRedX->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot,"  0-10% ALICE","pf");
        legendDirGammaRedX->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
        legendDirGammaRedX->AddEntry(graphCombDirGammaSpectrumSystErr2050Plot,"20-50% ALICE","pf");
        legendDirGammaRedX->Draw();

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

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum_LogX.%s",outputDir.Data(),suffix.Data()));

    TLegend* legendDirGammaTheoryPlusMcGillDiffY = new TLegend(0.20,0.382-0.95*0.85*textsizeLabelsDirGamma*9,0.5,0.382);
    legendDirGammaTheoryPlusMcGillDiffY->SetFillStyle(0);
    legendDirGammaTheoryPlusMcGillDiffY->SetFillColor(0);
    legendDirGammaTheoryPlusMcGillDiffY->SetLineColor(0);
//     legendDirGammaTheoryPlusMcGillDiffY->SetNColumns(2);
    legendDirGammaTheoryPlusMcGillDiffY->SetTextSize(0.85*textsizeLabelsDirGamma);
    legendDirGammaTheoryPlusMcGillDiffY->SetMargin(0.19);
    legendDirGammaTheoryPlusMcGillDiffY->SetTextFont(42);
    legendDirGammaTheoryPlusMcGillDiffY->AddEntry(graphTheoryMcGill0020Plot,"Paquet et al.","l");
    legendDirGammaTheoryPlusMcGillDiffY->AddEntry((TObject*)0,"arXiv:1509.06738","");
    legendDirGammaTheoryPlusMcGillDiffY->AddEntry(graphTheoryPHSD0020Plot,"Linnyk et al.","l");
    legendDirGammaTheoryPlusMcGillDiffY->AddEntry((TObject*)0,"arXiv:1504.05699","");
    legendDirGammaTheoryPlusMcGillDiffY->AddEntry(graphTheoryHe0020Plot,"v. Hees et al.","l");
    legendDirGammaTheoryPlusMcGillDiffY->AddEntry((TObject*)0,"NPA 933(2015) 256","");
    legendDirGammaTheoryPlusMcGillDiffY->AddEntry(graphTheoryChatterjee0010Plot,"Chatterjee et al.","l");
    legendDirGammaTheoryPlusMcGillDiffY->AddEntry((TObject*)0,"PRC 85(2012) 064910","");
    legendDirGammaTheoryPlusMcGillDiffY->AddEntry((TObject*)0,"+ JHEP 1305(2013) 030","");
    legendDirGammaTheoryPlusMcGillDiffY->Draw();

        graphTheoryMcGill0020Plot->Draw("p3lsame");
        graphTheoryPHSD0020Plot->Draw("p3lsame");
        graphTheoryChatterjee0010Plot->Draw("p3lsame");
    //     graphTheoryHees0020Plot->Draw("p3lsame");
        graphTheoryHe0020Plot->Draw("p3lsame");
        if (graphCombDirGammaSpectrumSystErr0010Plot){
            graphCombDirGammaSpectrumSystErr0010Plot->Draw("E2same");
        }
        if (graphCombDirGammaSpectrumStatErr0010Plot){
            graphCombDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }
        graphTheoryMcGill2040Plot->Draw("p3lsame");
        graphTheoryPHSD2040Plot->Draw("p3lsame");
        graphTheoryChatterjee2040Plot->Draw("p3lsame");
    //     graphTheoryHees2040Plot->Draw("p3lsame");
        graphTheoryHe2040Plot->Draw("p3lsame");

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

    //     graphTheoryChatterjee4060Plot->Draw("p3lsame");
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
    TH2D *dummDirGammaRedX ;
    dummDirGammaRedX = new TH2D("dummDirGammaRedX", "dummDirGammaRedX", 100000, 0., 4.25, 1000., 9e-6,9e2);
    SetStyleHistoTH2ForGraphs( dummDirGammaRedX, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                            0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummDirGammaRedX->GetXaxis()->SetLabelOffset(-0.002);
    dummDirGammaRedX->GetXaxis()->SetTickLength(0.025);
    dummDirGammaRedX->GetYaxis()->SetTickLength(0.025);
    dummDirGammaRedX->GetXaxis()->SetRangeUser(0.5,doubleRatioX[1]);
    dummDirGammaRedX->DrawCopy();


        SetStyleHisto(histoFitThermalGamma0010Stat, 3, styleFit, colorComb0010+1 );
        SetStyleHisto(histoFitThermalGamma2040Stat, 3, styleFit, colorComb2040+1 );

        TLatex *labelScalingDirGamma0010_3 = new TLatex(3.65,4E-1,"x 10^{2}");
        SetStyleTLatex( labelScalingDirGamma0010_3, 0.85*textsizeLabelsDirGamma,4,colorComb0010,42,kFALSE);
        labelScalingDirGamma0010_3->Draw();

        TLatex *labelScalingDirGamma2040_3 = new TLatex(3.65,1.7E-2,"x 10^{1}");
        SetStyleTLatex( labelScalingDirGamma2040_3, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
        labelScalingDirGamma2040_3->Draw();

        TLatex *labelScalingDirGamma2050_3 = new TLatex(3.65,4E-4,"x 10^{0}");
        SetStyleTLatex( labelScalingDirGamma2050_3, 0.85*textsizeLabelsDirGamma,4,colorComb2050,42,kFALSE);
        labelScalingDirGamma2050_3->Draw();
        TLatex *labelDirGammaCollRedX2 = new TLatex(0.41,0.93,Form("%s",collisionSystem.Data()));
        SetStyleTLatex( labelDirGammaCollRedX2, 0.85*textsizeLabelsDirGamma,4);
        labelDirGammaCollRedX2->Draw();


        TLegend* legendDirGammaRedXWithFit = new TLegend(0.455,0.915-0.9*0.85*textsizeLabelsDirGamma*3,0.555+0.21,0.915);
        legendDirGammaRedXWithFit->SetFillStyle(0);
        legendDirGammaRedXWithFit->SetFillColor(0);
        legendDirGammaRedXWithFit->SetLineColor(0);
        legendDirGammaRedXWithFit->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaRedXWithFit->SetMargin(0.2);
        legendDirGammaRedXWithFit->SetTextFont(42);
        legendDirGammaRedXWithFit->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot,"  0-10% ALICE","pf");
        legendDirGammaRedXWithFit->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
        legendDirGammaRedXWithFit->AddEntry(graphCombDirGammaSpectrumSystErr2050Plot,"20-50% ALICE","pf");
        legendDirGammaRedXWithFit->Draw();


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


        TLatex *labelDirGammaFitFunc = new TLatex(0.755,0.93,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})");
        SetStyleTLatex( labelDirGammaFitFunc, 0.85*textsizeLabelsDirGamma,4);
        labelDirGammaFitFunc->Draw();
        TLegend* legendDirGammaExpFit = new TLegend(0.795,0.915-0.9*0.85*textsizeLabelsDirGamma*2,0.96,0.915);
        legendDirGammaExpFit->SetFillStyle(0);
        legendDirGammaExpFit->SetFillColor(0);
        legendDirGammaExpFit->SetLineColor(0);
        legendDirGammaExpFit->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaExpFit->SetMargin(0.3);
        legendDirGammaExpFit->SetTextFont(42);
        legendDirGammaExpFit->AddEntry(histoFitThermalGamma0010Stat,"  0-10%","l");
        legendDirGammaExpFit->AddEntry(histoFitThermalGamma2040Stat,"20-40%","l");
        legendDirGammaExpFit->Draw();

        histoFitThermalGamma0010Stat->Draw("same,l");
        histoFitThermalGamma2040Stat->Draw("same,l");

    canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusFit_ReducedX.%s",outputDir.Data(),suffix.Data()));

        TLegend* legendDirGammaTheoryPlusMcGillRedX = new TLegend(0.19,0.295-1.*0.85*textsizeLabelsDirGamma*6,0.94,0.295);
        legendDirGammaTheoryPlusMcGillRedX->SetFillStyle(0);
        legendDirGammaTheoryPlusMcGillRedX->SetFillColor(0);
        legendDirGammaTheoryPlusMcGillRedX->SetLineColor(0);
        legendDirGammaTheoryPlusMcGillRedX->SetNColumns(2);
        legendDirGammaTheoryPlusMcGillRedX->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaTheoryPlusMcGillRedX->SetMargin(0.12);
        legendDirGammaTheoryPlusMcGillRedX->SetTextFont(42);
        legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryMcGill0020Plot,"Paquet et al.","l");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"arXiv:1509.06738","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryPHSD0020Plot,"Linnyk et al.","l");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"arXiv:1504.05699","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryChatterjee0010Plot,"Chatterjee et al.","l");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryHe0020Plot,"v. Hees et al.","l");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"PRC 85(2012) 064910","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"NPA 933(2015) 256","");
        legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"+ JHEP 1305(2013) 030","");
        legendDirGammaTheoryPlusMcGillRedX->Draw();

        dummDirGammaRedX->DrawCopy("same,axis");
        labelScalingDirGamma2050_3->Draw();
        graphTheoryMcGill0020Plot->Draw("p3lsame");
        graphTheoryPHSD0020Plot->Draw("p3lsame");
        graphTheoryChatterjee0010Plot->Draw("p3lsame");
        //     graphTheoryHees0020Plot->Draw("p3lsame");
        graphTheoryHe0020Plot->Draw("p3lsame");
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
        //     graphTheoryHees2040Plot->Draw("p3lsame");
        graphTheoryHe2040Plot->Draw("p3lsame");

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

        //     graphTheoryChatterjee4060Plot->Draw("p3lsame");
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
//     graphPHENIXAuAuStat0010->RemovePoint(graphPHENIXAuAuStat0010->GetN()-1);
//     graphPHENIXAuAuSys0010->RemovePoint(graphPHENIXAuAuSys0010->GetN()-1);
    TF1* fitThermalGamma0010PHENIXStat                = FitObject("e","fitThermalGamma0010PHENIXStat","Photon",graphPHENIXAuAuStat0010,0.6,2,NULL,"QNRMEX0+");
    fitThermalGamma0010PHENIXStat->FixParameter(1,0.239);
    graphPHENIXAuAuStat0010->Fit(fitThermalGamma0010PHENIXStat,"QNRMEX0+","",0.6,2.);
    cout << WriteParameterToFile(fitThermalGamma0010PHENIXStat)<< endl;

    TH1D* histoFitThermalGamma0010PHENIXStat        = (TH1D*)fitThermalGamma0010PHENIXStat->GetHistogram();
    SetStyleHisto(histoFitThermalGamma0010PHENIXStat, 3, styleFit, colorPHENIX );

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

    TH1D *dummDirGammaPHENIX0010 = new TH1D("dummDirGammaPHENIX0010", "dummDirGammaPHENIX0010", 1000, 0., 5.6);
    SetStyleHistoTH1ForGraphs( dummDirGammaPHENIX0010, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})", 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
    dummDirGammaPHENIX0010->GetXaxis()->SetLabelOffset(-0.002);
    dummDirGammaPHENIX0010->GetXaxis()->SetTickLength(0.025);
    dummDirGammaPHENIX0010->GetYaxis()->SetTickLength(0.025);
    dummDirGammaPHENIX0010->GetYaxis()->SetRangeUser(0.9e-5,50);
    dummDirGammaPHENIX0010->GetXaxis()->SetRangeUser(0.,doubleRatioX[1]);
    dummDirGammaPHENIX0010->DrawCopy();
    dummDirGammaPHENIX0010->DrawCopy("same,axis");

        DrawGammaSetMarkerTGraphAsym(graphPHENIXAuAuSys0010, markerStylePHENIX, markerSizePHENIX, colorPHENIX , colorPHENIX, widthLinesBoxes,  kTRUE);
        graphPHENIXAuAuSys0010->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPHENIXAuAuStat0010, markerStylePHENIX, markerSizePHENIX, colorPHENIX , colorPHENIX);
        graphPHENIXAuAuStat0010->Draw("p,E1Z,same");
        histoFitThermalGamma0010PHENIXStat->Draw("same,l");

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

        TLegend* legendDirGammaWithPHENIX0010 = new TLegend(0.455,0.935-1.15*0.85*textsizeLabelsDirGamma*8,0.555+0.21,0.935);
        legendDirGammaWithPHENIX0010->SetFillStyle(0);
        legendDirGammaWithPHENIX0010->SetFillColor(0);
        legendDirGammaWithPHENIX0010->SetLineColor(0);
        legendDirGammaWithPHENIX0010->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaWithPHENIX0010->SetMargin(0.2);
        legendDirGammaWithPHENIX0010->SetTextFont(42);
        legendDirGammaWithPHENIX0010->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot3,"ALICE","pf");
        legendDirGammaWithPHENIX0010->AddEntry((TObject*)0,Form("%s",collisionSystemCent0010.Data()),"");
        legendDirGammaWithPHENIX0010->AddEntry(histoFitThermalGamma0010Stat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
    //     legendDirGammaWithPHENIX0010->AddEntry((TObject*)0,"#it{T}_{eff} = 297 #pm 12^{#it{stat}} #pm 41^{#it{sys}}",""); // numbers for thermal
        legendDirGammaWithPHENIX0010->AddEntry((TObject*)0,"#it{T}_{eff} = 304 #pm 11^{ #it{stat}} #pm 40^{#it{sys}} MeV","");
        legendDirGammaWithPHENIX0010->AddEntry(graphPHENIXAuAuSys0010,"PHENIX","pf");
        legendDirGammaWithPHENIX0010->AddEntry((TObject*)0,Form("%s",collisionSystemRHIC0010.Data()),"");
        legendDirGammaWithPHENIX0010->AddEntry(histoFitThermalGamma0010PHENIXStat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
        legendDirGammaWithPHENIX0010->AddEntry((TObject*)0,"#it{T}_{eff} = 239 #pm 25^{#it{stat}} #pm 7^{ #it{sys}} MeV","");
        legendDirGammaWithPHENIX0010->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusPHENIX0010_ReducedX.%s",outputDir.Data(),suffix.Data()));

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
    dummDirGammaPHENIX2040->GetXaxis()->SetRangeUser(0.,doubleRatioX[1]);
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
        legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,"#it{T}_{eff} = 407 #pm 61^{ #it{stat}} #pm 96^{#it{sys}} MeV","");
        //     legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,"#it{T}_{eff} = 410 #pm 84^{#it{stat}} #pm 140^{#it{sys}}",""); // numbers for thermal
        legendDirGammaWithPHENIX2040->AddEntry(graphPHENIXAuAuSys2040,"PHENIX","pf");
        legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,Form("%s",collisionSystemRHIC2040.Data()),"");
        legendDirGammaWithPHENIX2040->AddEntry(histoFitThermalGamma2040PHENIXStat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
        legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,"#it{T}_{eff} = 260 #pm 33^{#it{ stat}} #pm 8^{#it{sys}} MeV","");
        legendDirGammaWithPHENIX2040->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusPHENIX2040_ReducedX.%s",outputDir.Data(),suffix.Data()));


    cout << "Plotting direct gamma PCM only " << __LINE__ << endl;
    //*******************************************************************************************************************************************
    //*************************************************** Plotting direct Gamma Spectrum PCM only ***********************************************
    //*******************************************************************************************************************************************

    TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr0010Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr0010Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr0010ArPlot = NULL;
    if (graphPCMDirGammaStatErr0010) graphPCMDirGammaSpectrumStatErr0010Plot         = ScaleGraph(graphPCMDirGammaStatErr0010,100);
    if (graphPCMDirGammaSysErr0010) graphPCMDirGammaSpectrumSystErr0010Plot         = ScaleGraph(graphPCMDirGammaSysErr0010,100);
    if (graphPCMDirGammaSumErrAr0010) graphPCMDirGammaSpectrumSumErr0010ArPlot     = ScaleGraph(graphPCMDirGammaSumErrAr0010,100);

    TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr2040Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr2040Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr2040ArPlot = NULL;
    if (graphPCMDirGammaStatErr2040) graphPCMDirGammaSpectrumStatErr2040Plot         = ScaleGraph(graphPCMDirGammaStatErr2040,10);
    if (graphPCMDirGammaSysErr2040) graphPCMDirGammaSpectrumSystErr2040Plot         = ScaleGraph(graphPCMDirGammaSysErr2040,10);
    if (graphPCMDirGammaSumErrAr2040) graphPCMDirGammaSpectrumSumErr2040ArPlot     = ScaleGraph(graphPCMDirGammaSumErrAr2040,10);
    TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr2050Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr2050Plot = NULL;
    TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr2050ArPlot = NULL;
    if (graphPCMDirGammaStatErr2050) graphPCMDirGammaSpectrumStatErr2050Plot         = ScaleGraph(graphPCMDirGammaStatErr2050,1);
    if (graphPCMDirGammaSysErr2050) graphPCMDirGammaSpectrumSystErr2050Plot         = ScaleGraph(graphPCMDirGammaSysErr2050,1);
    if (graphPCMDirGammaSumErrAr2050) graphPCMDirGammaSpectrumSumErr2050ArPlot     = ScaleGraph(graphPCMDirGammaSumErrAr2050,1);

    canvasDirGamma->SetLogx(1);
    canvasDirGamma->cd();
    dummyDirGamma->DrawCopy();

        graphTheoryEPS090010Plot->Draw("p3lsame");
        graphTheoryCT100010Plot->Draw("p3lsame");
        graphTheoryNLO0010Plot->Draw("p3lsame");
        graphTheoryPromptMcGill0020Plot->Draw("p3lsame");

        if (graphPCMDirGammaSpectrumSystErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr0010Plot->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr0010Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr0010Plot, markerStyleComb0010, markerSizeComb0010, colorComb0010 , colorComb0010);
            graphPCMDirGammaSpectrumStatErr0010Plot->Draw("p,E1Z,same");
        }
    //     if (graphPCMDirGammaSpectrumSumErr0010ArPlot){
    //         DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr0010ArPlot , 1, 3, colorComb0010, colorComb0010, 1.8, kTRUE);
    //         graphPCMDirGammaSpectrumSumErr0010ArPlot->Draw(">,same");
    //     }

        graphTheoryEPS092040Plot->Draw("p3lsame");
        graphTheoryCT102040Plot->Draw("p3lsame");
        graphTheoryNLO2040Plot->Draw("p3lsame");
        graphTheoryPromptMcGill2040Plot->Draw("p3lsame");

        if (graphPCMDirGammaSpectrumSystErr2040Plot){
            graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
            graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
            graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
            graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
            graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
            graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
            graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
            graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(5);
            graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(6);
            graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(6);
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
            graphPCMDirGammaSpectrumSystErr2040Plot->Draw("E2same");
        }
        if (graphPCMDirGammaSpectrumStatErr2040Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
            graphPCMDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSumErr2040ArPlot){
            graphPCMDirGammaSpectrumSumErr2040ArPlot->RemovePoint(0);
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr2040ArPlot , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
            graphPCMDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2040ArPlot);
        }
        graphTheoryEPS092050Plot->Draw("p3lsame");
        graphTheoryCT102050Plot->Draw("p3lsame");
        graphTheoryNLO2050Plot->Draw("p3lsame");

    //     if (graphPCMDirGammaSpectrumSystErr2050Plot){
    //         DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050, widthLinesBoxes, kTRUE);
    //         graphPCMDirGammaSpectrumSystErr2050Plot->Draw("E2same");
    //     }
        if (graphPCMDirGammaSpectrumStatErr2050Plot){
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2050Plot, markerStyleComb2050, markerSizeComb2050, colorComb2050 , colorComb2050);
            graphPCMDirGammaSpectrumStatErr2050Plot->Draw("p,E1Z,same");
        }
        if (graphPCMDirGammaSpectrumSumErr2050ArPlot){

            graphPCMDirGammaSpectrumSumErr2050ArPlot->RemovePoint(0);
            graphPCMDirGammaSpectrumSumErr2050ArPlot->Print();
            DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr2050ArPlot , 1, 3, colorComb2050, colorComb2050, 1.8, kTRUE);
            graphPCMDirGammaSpectrumSumErr2050ArPlot->Draw(">,same");
            PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2050ArPlot);
        }


        labelScalingDirGamma0010->Draw();
        labelScalingDirGamma2040->Draw();
        labelScalingDirGamma2050->Draw();
        labelDirGammaColl->Draw();

        legendDirGamma->Draw();
        legendDirGammaTheory->Draw();
        legendDirGammaTheory2->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum_PCMonly.%s",outputDir.Data(),suffix.Data()));


    cout << "Plotting direct gamma PHOS only " << __LINE__ << endl;
    //*******************************************************************************************************************************************
    //*************************************************** Plotting direct Gamma Spectrum PCM only ***********************************************
    //*******************************************************************************************************************************************
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

    canvasDirGamma->SetLogx(1);
    canvasDirGamma->cd();
    dummyDirGamma->DrawCopy();

        graphTheoryEPS090010Plot->Draw("p3lsame");
        graphTheoryCT100010Plot->Draw("p3lsame");
        graphTheoryNLO0010Plot->Draw("p3lsame");
        graphTheoryPromptMcGill0020Plot->Draw("p3lsame");

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
        graphTheoryPromptMcGill2040Plot->Draw("p3lsame");

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

        legendDirGamma->Draw();
        legendDirGammaTheory->Draw();
        legendDirGammaTheory2->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum_PHOSonly.%s",outputDir.Data(),suffix.Data()));


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
        labelScalingDirGamma2050_2->Draw();

        TLatex *labelDirGammaColl_2 = new TLatex(0.26,0.94,Form("%s",collisionSystem.Data()));
        SetStyleTLatex( labelDirGammaColl_2, 0.85*textsizeLabelsDirGamma,4);
        labelDirGammaColl_2->Draw();

        TLegend* legendDirGammaLinX = new TLegend(0.26,0.93-1.1*0.83*textsizeLabelsDirGamma*3,0.26+0.21,0.93);
        legendDirGammaLinX->SetFillStyle(0);
        legendDirGammaLinX->SetFillColor(0);
        legendDirGammaLinX->SetLineColor(0);
        legendDirGammaLinX->SetTextSize(0.85*textsizeLabelsDirGamma);
        legendDirGammaLinX->SetMargin(0.25);
        legendDirGammaLinX->SetTextFont(42);
        legendDirGammaLinX->AddEntry(graphCombDirGammaSpectrumSystErr0010Plot,"  0-10% ALICE","pf");
        legendDirGammaLinX->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
        legendDirGammaLinX->AddEntry(graphCombDirGammaSpectrumSystErr2050Plot,"20-50% ALICE","pf");
        legendDirGammaLinX->Draw();

    canvasDirGammaLinX->Print(Form("%s/DirGammaSpectrum_LinX_withoutFit.%s",outputDir.Data(),suffix.Data()));

        graphTheoryEPS090010Plot->Draw("p3lsame");
        graphTheoryCT100010Plot->Draw("p3lsame");
        graphTheoryNLO0010Plot->Draw("p3lsame");
        graphTheoryPromptMcGill0020Plot->Draw("p3lsame");
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

        TLegend* legendDirGammaTheoryLinX = new TLegend(0.53,0.93-1.*0.83*textsizeLabelsDirGamma*3,0.53+0.21,0.93);
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

        TLegend* legendDirGammaTheoryLinX2 = new TLegend(0.53,0.93-1.*0.83*textsizeLabelsDirGamma*7-0.02,0.53+0.21,0.93-1.*0.85*textsizeLabelsDirGamma*3-0.02);
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


//     cout << "Plotting cocktial contributions " << __LINE__ << endl;
    // ************************************************************************************************************************************************
    // ************************************************************ Plot Cocktail contributions *******************************************************
    // ************************************************************************************************************************************************
//     histoCocktailSumGamma->Rebin(4);
//     histoCocktailPi0Gamma->Rebin(4);
//     histoCocktailEtaGamma->Rebin(4);
//     histoCocktailOmegaGamma->Rebin(4);
//     histoCocktailEtaPGamma->Rebin(4);
//     histoCocktailRhoGamma->Rebin(4);
//     histoCocktailPhiGamma->Rebin(4);
//     histoCocktailSigmaGamma->Rebin(4);
//
//     TH1D* histoCocktailRatioPi0GammaSumGamma = (TH1D*)histoCocktailPi0Gamma->Clone("histoCocktailRatioPi0GammaSumGamma");
//     histoCocktailRatioPi0GammaSumGamma->Divide(histoCocktailRatioPi0GammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// //     histoCocktailRatioPi0GammaSumGamma->Scale(100);
//     TH1D* histoCocktailRatioEtaGammaSumGamma = (TH1D*)histoCocktailEtaGamma->Clone("histoCocktailRatioEtaGammaSumGamma");
//     histoCocktailRatioEtaGammaSumGamma->Divide(histoCocktailRatioEtaGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// //     histoCocktailRatioEtaGammaSumGamma->Scale(100);
//     TH1D* histoCocktailRatioOmegaGammaSumGamma = (TH1D*)histoCocktailOmegaGamma->Clone("histoCocktailRatioOmegaGammaSumGamma");
//     histoCocktailRatioOmegaGammaSumGamma->Divide(histoCocktailRatioOmegaGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// //     histoCocktailRatioOmegaGammaSumGamma->Scale(100);
//     TH1D* histoCocktailRatioPhiGammaSumGamma = (TH1D*)histoCocktailPhiGamma->Clone("histoCocktailRatioPhiGammaSumGamma");
//     histoCocktailRatioPhiGammaSumGamma->Divide(histoCocktailRatioPhiGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// //     histoCocktailRatioPhiGammaSumGamma->Scale(100);
//     TH1D* histoCocktailRatioEtaPGammaSumGamma = (TH1D*)histoCocktailEtaPGamma->Clone("histoCocktailRatioEtaPGammaSumGamma");
//     histoCocktailRatioEtaPGammaSumGamma->Divide(histoCocktailRatioEtaPGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// //     histoCocktailRatioEtaPGammaSumGamma->Scale(100);
//     TH1D* histoCocktailRatioRhoGammaSumGamma = (TH1D*)histoCocktailRhoGamma->Clone("histoCocktailRatioRhoGammaSumGamma");
//     histoCocktailRatioRhoGammaSumGamma->Divide(histoCocktailRatioRhoGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// //     histoCocktailRatioRhoGammaSumGamma->Scale(100);
//     TH1D* histoCocktailRatioSigmaGammaSumGamma = (TH1D*)histoCocktailSigmaGamma->Clone("histoCocktailRatioSigmaGammaSumGamma");
//     histoCocktailRatioSigmaGammaSumGamma->Divide(histoCocktailRatioSigmaGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// //     histoCocktailRatioSigmaGammaSumGamma->Scale(100);
//
//     TCanvas *canvasCocktailRatio = new TCanvas("canvasCocktailRatio","",10,10,1200,1200);  // gives the page size
//     DrawGammaCanvasSettings( canvasCocktailRatio, 0.12, 0.01, 0.01, 0.08);
//     canvasCocktailRatio->SetLogy();
//     canvasCocktailRatio->SetLogx();
//
//     Int_t textSizeLabelsPixelRatioCock = 48;
//     Double_t textsizeLabelsRatioCock = 0;
//     if (canvasCocktailRatio->XtoPixel(canvasCocktailRatio->GetX2()) < canvasCocktailRatio->YtoPixel(canvasCocktailRatio->GetY1())){
//         textsizeLabelsRatioCock = (Double_t)textSizeLabelsPixelRatioCock/canvasCocktailRatio->XtoPixel(canvasCocktailRatio->GetX2()) ;
//     } else {
//         textsizeLabelsRatioCock = (Double_t)textSizeLabelsPixelRatioCock/canvasCocktailRatio->YtoPixel(canvasCocktailRatio->GetY1());
//     }
//
//
//     TH2D *dummyCocktailRatio ;
//     dummyCocktailRatio = new TH2D("dummyCocktailRatio", "dummyCocktailRatio", 1000, 0., 22, 1000., 4e-5,14);
//     SetStyleHistoTH2ForGraphs( dummyCocktailRatio, "#it{p}_{T} (GeV/#it{c})", "#gamma_{source}/#gamma_{decay}",
//                             0.85*textsizeLabelsRatioCock, textsizeLabelsRatioCock, 0.85*textsizeLabelsRatioCock, textsizeLabelsRatioCock, 0.82, 1.3);
//     dummyCocktailRatio->GetXaxis()->SetLabelOffset(-0.011);
// //     dummyCocktailRatio->GetXaxis()->SetMoreLogLabels();
//     dummyCocktailRatio->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
//
//     dummyCocktailRatio->DrawCopy();
//         SetStyleHisto(histoCocktailRatioPi0GammaSumGamma, widthCommonFit*1.5, 1, colorCocktailPi0 );
//         SetStyleHisto(histoCocktailRatioEtaGammaSumGamma, widthCommonFit*1.5, 7, colorCocktailEta );
//         SetStyleHisto(histoCocktailRatioEtaPGammaSumGamma, widthCommonFit*1.5, 2, colorCocktailEtaP );
//         SetStyleHisto(histoCocktailRatioOmegaGammaSumGamma, widthCommonFit*1.5, 4, colorCocktailOmega );
//         SetStyleHisto(histoCocktailRatioPhiGammaSumGamma, widthCommonFit*1.5, 5, colorCocktailPhi );
//         SetStyleHisto(histoCocktailRatioRhoGammaSumGamma, widthCommonFit*1.5, 8, colorCocktailRho0 );
// //         SetStyleHisto(histoCocktailRatioSigmaGammaSumGamma, widthCommonFit*1.5, 3, colorCocktailSigma0 );
//
//         histoCocktailRatioPi0GammaSumGamma->Draw("chistsame");
//         histoCocktailRatioEtaGammaSumGamma->Draw("chistsame");
//         histoCocktailRatioEtaPGammaSumGamma->Draw("chistsame");
//         histoCocktailRatioOmegaGammaSumGamma->Draw("chistsame");
//         histoCocktailRatioPhiGammaSumGamma->Draw("chistsame");
//         histoCocktailRatioRhoGammaSumGamma->Draw("chistsame");
// //         histoCocktailRatioSigmaGammaSumGamma->Draw("chistsame");
//
//     //     TLatex* tpi = new TLatex(0.18,0.92,"#pi^{0} #rightarrow #gamma#gamma (e^{+}e^{-}#gamma)");
//     //     SetStyleTLatex( tpi, 0.85*textsizeLabelsRatioCock,4,colorCocktailPi0,42);
//     //     tpi->Draw();
//     //     TLatex* teta = new TLatex(0.18,0.88,"#eta #rightarrow #gamma#gamma (#pi^{+}#pi^{-}#gamma,e^{+}e^{-}#gamma,#pi^{0}#gamma#gamma)");
//     //     SetStyleTLatex( teta, 00.85*textsizeLabelsRatioCock,4,colorCocktailEta,42);
//     //     teta->Draw();
//     //     TLatex* tomega = new TLatex(0.18,0.84,"#omega #rightarrow #pi^{0}#gamma (#eta#gamma)");
//     //     SetStyleTLatex( tomega, 0.85*textsizeLabelsRatioCock,4,colorCocktailOmega,42);
//     //     tomega->Draw();
//     //     TLatex* tetaprime = new TLatex(0.18,0.80,"#eta' #rightarrow #rho#gamma (#omega#gamma, #gamma#gamma)");
//     //     SetStyleTLatex( tetaprime, 0.85*textsizeLabelsRatioCock,4,colorCocktailEtaP,42);
//     //     tetaprime->Draw();
//     //     TLatex* tphi = new TLatex(0.65,0.92,"#phi #rightarrow #eta#gamma (#pi^{0}#gamma, #omega#gamma)");
//     //     SetStyleTLatex( tphi, 0.85*textsizeLabelsRatioCock,4,colorCocktailPhi,42);
//     //     tphi->Draw();
//     //     TLatex* trho = new TLatex(0.65,0.88,"#rho^{0} #rightarrow #pi^{+}#pi^{-}#gamma (#pi^{0}#gamma, #eta#gamma)");
//     //     SetStyleTLatex( trho, 0.85*textsizeLabelsRatioCock,4,colorCocktailRho0,42);
//     //     trho->Draw();
//     //     TLatex* tSigma = new TLatex(0.65,0.84,"#Sigma^{0} #rightarrow #Lambda#gamma (#Lambda#gamma#gamma)");
//     //     SetStyleTLatex( tSigma, 0.85*textsizeLabelsRatioCock,4,colorCocktailSigma0,42);
//     //     tSigma->Draw();
//
//         TLegend* legendCocktail = new TLegend(0.15,0.96-1.1*0.85*textsizeLabelsRatioCock*4,0.18+0.75,0.96);
//         legendCocktail->SetFillStyle(0);
//         legendCocktail->SetFillColor(0);
//         legendCocktail->SetLineColor(0);
//         legendCocktail->SetTextSize(0.85*textsizeLabelsRatioCock);
//         legendCocktail->SetMargin(0.15);
//         legendCocktail->SetTextFont(42);
//         legendCocktail->SetNColumns(2);
//         legendCocktail->AddEntry(histoCocktailRatioPi0GammaSumGamma,"#pi^{0} #rightarrow #gamma#gamma (e^{+}e^{-}#gamma)","l");
//         legendCocktail->AddEntry(histoCocktailRatioEtaPGammaSumGamma,"#eta' #rightarrow #rho#gamma (#omega#gamma, #gamma#gamma)","l");
//         legendCocktail->AddEntry(histoCocktailRatioEtaGammaSumGamma,"#eta #rightarrow #gamma#gamma (#pi^{+}#pi^{-}#gamma,e^{+}e^{-}#gamma,#pi^{0}#gamma#gamma)","l");
//         legendCocktail->AddEntry(histoCocktailRatioPhiGammaSumGamma,"#phi #rightarrow #eta#gamma (#pi^{0}#gamma, #omega#gamma)","l");
//         legendCocktail->AddEntry(histoCocktailRatioOmegaGammaSumGamma,"#omega #rightarrow #pi^{0}#gamma (#eta#gamma)","l");
//         legendCocktail->AddEntry(histoCocktailRatioRhoGammaSumGamma,"#rho^{0} #rightarrow #pi^{+}#pi^{-}#gamma (#pi^{0}#gamma, #eta#gamma)","l");
// //         legendCocktail->AddEntry(histoCocktailRatioSigmaGammaSumGamma,"#Sigma^{0} #rightarrow #Lambda#gamma (#Lambda#gamma#gamma)","l");
//
//         legendCocktail->Draw();
//
//         TLatex *labelCocktailRatioEnergy = new TLatex(0.16,0.13,collisionSystemCent0010.Data());
//         SetStyleTLatex( labelCocktailRatioEnergy, 0.85*textsizeLabelsRatioCock,4);
//         labelCocktailRatioEnergy->Draw();
//
//         TLatex *labelCocktailSim = new TLatex(0.625,0.13,"Monte Carlo simulation");
//         SetStyleTLatex( labelCocktailSim, 0.85*textsizeLabelsRatioCock,4);
//         labelCocktailSim->Draw();
//
//
//     canvasCocktailRatio->Print(Form("%s/CocktailRatioToSumPCM_0010.%s",outputDir.Data(),suffix.Data()));

    cout << "Plotting direct gamma RAA " << __LINE__ << endl;
    // *****************************************************************************************************
    // ***************************************** Plotting RAA 0-10% ****************************************
    // *****************************************************************************************************
    TGraphAsymmErrors* graphCombRAADirGammaSys0010Plot = NULL;
    if (graphCombRAADirGammaSys0010){
        graphCombRAADirGammaSys0010Plot = (TGraphAsymmErrors*)graphCombRAADirGammaSys0010->Clone("graphCombRAADirGammaSys0010Plot");
    }
    TGraphAsymmErrors* graphCombRAADirGammaStat0010Plot = NULL;
    if (graphCombRAADirGammaStat0010){
        graphCombRAADirGammaStat0010Plot = (TGraphAsymmErrors*)graphCombRAADirGammaStat0010->Clone("graphCombRAADirGammaStat0010Plot");
        ProduceGraphAsymmWithoutXErrors(graphCombRAADirGammaStat0010Plot);
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
    histo2DRAADummy->GetXaxis()->SetMoreLogLabels();
    histo2DRAADummy->DrawCopy("");

        TLatex *labelRAAEnergy0010 = new TLatex(0.48,0.92,collisionSystemCent0010.Data());
        SetStyleTLatex( labelRAAEnergy0010, 0.85*textsizeLabelsRAA,4);
        labelRAAEnergy0010->Draw();

        TBox* boxErrorNorm0010_Single = CreateBoxConv(colorComb0010Box, 0.75, 1.-normErr0010 , 0.8, 1.+normErr0010);
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

    TCanvas* canvasRAA_2050 = new TCanvas("canvasRAA_2050","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRAA_2050,  0.1, 0.01, 0.015, 0.1);
    canvasRAA_2050->SetLogx(1);

    histo2DRAADummy->DrawCopy("");

        TBox* boxErrorNorm2050_Single = CreateBoxConv(colorComb2050Box, 0.75, 1.-normErr2050 , 0.8, 1.+normErr2050);
        boxErrorNorm2050_Single->Draw();

        TLatex *labelRAAEnergy2050 = new TLatex(0.48,0.92,collisionSystemCent2050.Data());
        SetStyleTLatex( labelRAAEnergy2050, 0.85*textsizeLabelsRAA,4);
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
        Int_t maxBinSig0010 = 5;
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