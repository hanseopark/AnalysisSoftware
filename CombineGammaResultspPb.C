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


void CombineGammaResultspPb(    TString inputFileNamePCM    = "", 
                                TString inputFileNamePHOS   = "", 
                                TString inputFileNameEMC    = "", 
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
    TString fileNameTheorypPb                                   = "ExternalInputPbPb/Theory/TheoryCompilationpPb.root";
//     TString fileNameExperimentpPb                               = "ExternalInputPbPb/OtherExperiments/phenix_200.root";
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPHOSGammapPb.root", inputFileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMGammapPb.root", inputFileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCGammapPb.root", inputFileNameEMC.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/TheorypPb.root", fileNameTheorypPb.Data(), outputDir.Data()));
//     gSystem->Exec(Form("cp %s %s/OtherExperimentsPbPb.root", fileNameExperimentPbPb.Data(), outputDir.Data()));
    
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
    doubleRatioXpp[0]   = 0.23;     doubleRatioXpp[1]   = 30;
        
    Color_t colorCocktailPi0                                    = kRed+2;
    Color_t colorCocktailEta                                    = kBlue+1;
    Color_t colorCocktailEtaP                                   = kOrange+1;
    Color_t colorCocktailOmega                                  = kYellow+2;
    Color_t colorCocktailPhi                                    = kViolet;
    Color_t colorCocktailRho0                                   = kAzure-2;
    Color_t colorCocktailSigma0                                 = kGray+1;
    
    Color_t colorNLOcalc                                        = kBlue-7; // kBlack
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
    Color_t colorEMC                                            = GetDefaultColorDiffDetectors("EMCal", kFALSE, kFALSE, kTRUE);
    Color_t colorPCMBox                                         = GetDefaultColorDiffDetectors("PCM", kFALSE, kTRUE, kTRUE);
    Color_t colorPHOSBox                                        = GetDefaultColorDiffDetectors("PHOS", kFALSE, kTRUE, kTRUE);
    Color_t colorEMCBox                                         = GetDefaultColorDiffDetectors("EMCal", kFALSE, kTRUE, kTRUE);
    Style_t markerStylePCM                                      = GetDefaultMarkerStyleDiffDetectors("PCM", kFALSE);
    Style_t markerStylePHOS                                     = GetDefaultMarkerStyleDiffDetectors("PHOS", kFALSE);
    Style_t markerStyleEMC                                      = GetDefaultMarkerStyleDiffDetectors("EMCal", kFALSE);
    Size_t markerSizePCM                                        = GetDefaultMarkerSizeDiffDetectors("PCM", kFALSE);
    Size_t markerSizePHOS                                       = GetDefaultMarkerSizeDiffDetectors("PHOS", kFALSE);
    Size_t markerSizeEMC                                        = GetDefaultMarkerSizeDiffDetectors("EMCal", kFALSE);
    
    Color_t colorCombpPb                                        = GetColorDefaultColor("pPb_5.023TeV", "", "");
    Color_t colorCombpPbBox                                     = GetColorDefaultColor("pPb_5.023TeV", "", "", kTRUE);
    
    Style_t markerStyleCombpPb                                  = GetDefaultMarkerStyle("pPb_5.023TeV", "", "");
    Size_t markerSizeCombpPb                                    = GetDefaultMarkerSize("pPb_5.023TeV", "", "");

    Width_t widthLinesBoxes                                     = 1.4;
    Width_t widthCommonFit                                      = 2.4;
    
    TString collisionSystempPb                                  = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    

    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from pPb PCM file ******************************************************
    //*******************************************************************************************************************************************
    TFile* filePCMGammapPb                                  = new TFile( inputFileNamePCM.Data());
//     //________________________________________________ Load PCM pPb 5.023TeV _________________________________________________________________________
    TDirectory* directoryPCMGammapPb                = (TDirectory*)filePCMGammapPb->Get("Gamma_pPb5TeV");    
        TH1D* histoPCMDRPi0FitStatErrpPb                        = (TH1D*) directoryPCMGammapPb->Get("DoubleRatioPi0FitStatError");
        TGraphAsymmErrors* graphPCMDRPi0FitSysErrpPb            = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("DoubleRatioPi0FitSystError");
        TH1D* histoPCMIncGammaStatErrpPb                        = (TH1D*) directoryPCMGammapPb->Get("IncRatioStatError");
        TGraphAsymmErrors* graphPCMIncGammaSysErrpPb            = (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncRatioSystError");

    // loading theory
    TDirectory* directoryTheorypPb                              = (TDirectory*)filePCMGammapPb->Get("Theory"); 
        TGraphAsymmErrors* graphTheoryNLODRpPb                  = (TGraphAsymmErrors*) directoryTheorypPb->Get("NLODoubleRatio_pPb5TeV");
        TGraphAsymmErrors* graphTheoryNLOpPb                    = (TGraphAsymmErrors*) directoryTheorypPb->Get("NLODirGamma_pPb5TeV");

    //*******************************************************************************************************************************************
    //*********************************************** Load PCM histograms from pPb PCM file ******************************************************
    //*******************************************************************************************************************************************
    TFile* fileEMCGammapPb                                  = new TFile( inputFileNameEMC.Data());
//     //________________________________________________ Load EMC pPb 5.023TeV _________________________________________________________________________
    TDirectory* directoryEMCGammapPb                = (TDirectory*)fileEMCGammapPb->Get("Gamma_pPb5TeV");    
        TH1D* histoEMCDRPi0FitStatErrpPb                        = (TH1D*) directoryEMCGammapPb->Get("DoubleRatioPi0FitStatError");
        TGraphAsymmErrors* graphEMCDRPi0FitSysErrpPb            = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("DoubleRatioPi0FitSystError");
        TH1D* histoEMCIncGammaStatErrpPb                        = (TH1D*) directoryEMCGammapPb->Get("IncRatioStatError");
        TGraphAsymmErrors* graphEMCIncGammaSysErrpPb            = (TGraphAsymmErrors*) directoryEMCGammapPb->Get("IncRatioSystError");

        
//     //*******************************************************************************************************************************************
//     //*********************************************** Load PHOS histograms from PHOS file *******************************************************
//     //*******************************************************************************************************************************************        
    TFile* filePHOSGamma                            = new TFile( inputFileNamePHOS.Data());
//     TDirectorY* directoryPHOSGamma0020              = (TDirectory*)filePHOSGamma->Get("PHOS_PbPb_2760_Centrality_00-20"); 
        TH1D* histoPHOSDRPi0FitStatErrpPb                      = (TH1D*) filePHOSGamma->Get("hGamma_PbPb_cen6_Stat");
        TH1D* histoPHOSDRPi0FitSysErrpPb                       = (TH1D*) filePHOSGamma->Get("hGamma_PbPb_cen6_SystRatio");
        TGraphAsymmErrors* graphPHOSDRPi0FitSysErrpPb          = new TGraphAsymmErrors(histoPHOSDRPi0FitSysErrpPb);    
// 
//     //*******************************************************************************************************************************************
//     //***************************************************** Combine DR of PCM and PHOS **********************************************************
//     //*******************************************************************************************************************************************        
//     Double_t newBinsComb[21]                                    = { 0.9, 1.1, 1.3, 1.5, 1.7, 
//                                                                     1.9, 2.1, 2.3, 2.5, 2.7, 
//                                                                     3.0, 3.3, 3.7, 4.1, 4.6, 
//                                                                     5.4, 6.2, 7.0, 8.0, 11.0, 
//                                                                     14.0};
//     TGraphAsymmErrors *graphCombDRPi0FitSysErr0020;
//     TGraphAsymmErrors *graphCombDRPi0FitSysAErr0020;
//     TGraphAsymmErrors *graphCombDRPi0FitSysBErr0020;
//     TGraphAsymmErrors *graphCombDRPi0FitSysCErr0020;
//     TGraphAsymmErrors *graphCombDRPi0FitStatErr0020;
//     TGraphAsymmErrors *graphCombDRPi0FitSumErr0020;        
//     graphCombDRPi0FitSumErr0020     = CombinePtPointsSpectraAdv(    histoPCMDRPi0FitStatErr0020, graphPCMDRPi0FitSysErr0020, 
//                                                                     graphPCMDRPi0FitSysAErr0020, graphPCMDRPi0FitSysBErr0020, graphPCMDRPi0FitSysCErr0020,
//                                                                     histoPHOSDRPi0FitStatErrpPb, graphPHOSDRPi0FitSysErrpPb ,
//                                                                     graphPHOSDRPi0FitSysAErr0020, graphPHOSDRPi0FitSysBErr0020, graphPHOSDRPi0FitSysCErr0020,
//                                                                     graphCombDRPi0FitStatErr0020, graphCombDRPi0FitSysErr0020, 
//                                                                     graphCombDRPi0FitSysAErr0020, graphCombDRPi0FitSysBErr0020, graphCombDRPi0FitSysCErr0020,
//                                                                     newBinsComb, 21, 0, 0, -1);
//     graphCombDRPi0FitSumErr0020->Print();
//     TGraphAsymmErrors* graphCombDRPi0FitStatSysAErr0020         = AddErrorsOfGraphsQuadratically (graphCombDRPi0FitStatErr0020, graphCombDRPi0FitSysAErr0020);
//     Double_t SysCCombDRPi0Fit0020                               = graphCombDRPi0FitSysCErr0020->GetErrorYlow(4)/graphCombDRPi0FitSysCErr0020->GetY()[4];
// 
//     
//     //*******************************************************************************************************************************************
//     //******************************************** Combine Inclusive gamma of PCM and PHOS ******************************************************
//     //*******************************************************************************************************************************************        
//     //__________________________________________________ 0-20% combine spectra __________________________________________________________________
//     TGraphAsymmErrors *graphCombIncGammaSysErr0020;
//     TGraphAsymmErrors *graphCombIncGammaSysAErr0020;
//     TGraphAsymmErrors *graphCombIncGammaSysBErr0020;
//     TGraphAsymmErrors *graphCombIncGammaSysCErr0020;
//     TGraphAsymmErrors *graphCombIncGammaStatErr0020;
//     TGraphAsymmErrors *graphCombIncGammaSumErr0020;        
//     graphCombIncGammaSumErr0020     = CombinePtPointsSpectraAdv(    histoPCMIncGammaStatErr0020, graphPCMIncGammaSysErr0020,
//                                                                     graphPCMIncGammaSysAErr0020, graphPCMIncGammaSysBErr0020, graphPCMIncGammaSysCErr0020,
//                                                                     histoPHOSIncGammaStatErr0020, graphPHOSIncGammaSysErr0020 ,
//                                                                     graphPHOSIncGammaSysAErr0020, graphPHOSIncGammaSysBErr0020, graphPHOSIncGammaSysCErr0020,
//                                                                     graphCombIncGammaStatErr0020, graphCombIncGammaSysErr0020,
//                                                                     graphCombIncGammaSysAErr0020, graphCombIncGammaSysBErr0020, graphCombIncGammaSysCErr0020,
//                                                                     newBinsComb, 21, 0, 0, -1);
//     
//     //__________________________________________ 0-20% fit combined and build ratio of individual to fit ________________________________________
//     TF1* fitIncGammaCombQCD0020                                 = FitObject("qcd","fitIncGammaCombQCD0020","Photon",graphCombIncGammaSumErr0020,0.9,14,NULL,"QNRMEX0+");
//     cout << WriteParameterToFile(fitIncGammaCombQCD0020)<< endl;   
//     TH1D* histoFitQCDIncGammaComb0020                           = (TH1D*)fitIncGammaCombQCD0020->GetHistogram();
//     
//     TGraphAsymmErrors* graphPCMIncGammaStatErr0020              = new TGraphAsymmErrors(histoPCMIncGammaStatErr0020);
//     graphPCMIncGammaStatErr0020->RemovePoint(0);
//     TGraphAsymmErrors* graphPCMIncGammaStatSysAErr0020          = AddErrorsOfGraphsQuadratically (graphPCMIncGammaStatErr0020, graphPCMIncGammaSysAErr0020);
//     
//     TGraphAsymmErrors* graphPHOSIncGammaStatErr0020             = new TGraphAsymmErrors(histoPHOSIncGammaStatErr0020);
//     graphPHOSIncGammaStatErr0020->RemovePoint(0);
//     graphPHOSIncGammaSysErr0020->RemovePoint(0);
//     graphPHOSIncGammaSysAErr0020->RemovePoint(0);
//     graphPHOSIncGammaSysBErr0020->RemovePoint(0);
//     graphPHOSIncGammaSysCErr0020->RemovePoint(0);
//     TGraphAsymmErrors* graphPHOSIncGammaStatSysAErr0020         = AddErrorsOfGraphsQuadratically (graphPHOSIncGammaStatErr0020, graphPHOSIncGammaSysAErr0020);
//     
//     
//     
//     TGraphAsymmErrors* graphRatioCombPHOSIncGammaStatErr0020    = (TGraphAsymmErrors*) graphPHOSIncGammaStatErr0020->Clone();    
//     TGraphAsymmErrors* graphRatioCombPHOSIncGammaSysErr0020     = (TGraphAsymmErrors*) graphPHOSIncGammaSysErr0020->Clone();    
//     TGraphAsymmErrors* graphRatioCombPCMIncGammaStatErr0020     = (TGraphAsymmErrors*) graphPCMIncGammaStatErr0020->Clone();    
//     TGraphAsymmErrors* graphRatioCombPCMIncGammaSysErr0020      = (TGraphAsymmErrors*) graphPCMIncGammaSysErr0020->Clone();    
//     graphRatioCombPHOSIncGammaStatErr0020                       = CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaStatErr0020, fitIncGammaCombQCD0020); 
//     graphRatioCombPHOSIncGammaSysErr0020                        = CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaSysErr0020, fitIncGammaCombQCD0020); 
//     graphRatioCombPCMIncGammaStatErr0020                        = CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaStatErr0020, fitIncGammaCombQCD0020); 
//     graphRatioCombPCMIncGammaSysErr0020                         = CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaSysErr0020, fitIncGammaCombQCD0020); 
//     
//     TGraphAsymmErrors* graphRatioCombCombFitIncGammaStatErr0020 = (TGraphAsymmErrors*)graphCombIncGammaStatErr0020->Clone();
//     TGraphAsymmErrors* graphRatioCombCombFitIncGammaSysErr0020  = (TGraphAsymmErrors*)graphCombIncGammaSysErr0020->Clone();
//     graphRatioCombCombFitIncGammaStatErr0020                    = CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaStatErr0020, fitIncGammaCombQCD0020); 
//     graphRatioCombCombFitIncGammaSysErr0020                     = CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaSysErr0020, fitIncGammaCombQCD0020);     
// 
//     
//     //*******************************************************************************************************************************************
//     //**************************************************** Calculate direct photon spectrum 0020 ************************************************
//     //*******************************************************************************************************************************************
//     cout << endl;
//     graphCombIncGammaStatErr0020->Print();
//     Double_t xArrayCombined[graphCombDRPi0FitStatErr0020->GetN()+1];
//     xArrayCombined[0] = graphCombDRPi0FitStatErr0020->GetX()[0] - graphCombDRPi0FitStatErr0020->GetEXhigh()[0];
//     cout << "Binning \n" << xArrayCombined[0] << endl;
//     for (Int_t i = 1; i<graphCombDRPi0FitStatErr0020->GetN()+1;i++){
//         xArrayCombined[i] = graphCombDRPi0FitStatErr0020->GetX()[i-1] + graphCombDRPi0FitStatErr0020->GetEXhigh()[i-1];
//         cout << xArrayCombined[i] << endl;
//     }    
//     
//     //_______________________ copy inclusive photon spectra _____________________
//     TH1D *histoCombDirGammaSpectrumErrSum0020                   = new TH1D("histoCombDirGammaSpectrumErrSum0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     TH1D *histoCombDirGammaSpectrumErrSys0020                   = new TH1D("histoCombDirGammaSpectrumErrSys0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     TH1D *histoCombDirGammaSpectrumErrSysA0020                  = new TH1D("histoCombDirGammaSpectrumErrSysA0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     TH1D *histoCombDirGammaSpectrumErrSysB0020                  = new TH1D("histoCombDirGammaSpectrumErrSysB0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     TH1D *histoCombDirGammaSpectrumErrSysC0020                  = new TH1D("histoCombDirGammaSpectrumErrSysC0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     TH1D *histoCombDirGammaSpectrumErrStat0020                  = new TH1D("histoCombDirGammaSpectrumErrStat0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
// 
//     //_______________________ get arrays of double ratio errors __________________
//     Double_t *SystErrorsCombDR0020                              = new Double_t[graphCombIncGammaStatErr0020->GetN()];
//     Double_t *SystAErrorsCombDR0020                             = new Double_t[graphCombIncGammaStatErr0020->GetN()];
//     Double_t *SystBErrorsCombDR0020                             = new Double_t[graphCombIncGammaStatErr0020->GetN()];
//     Double_t *SystCErrorsCombDR0020                             = new Double_t[graphCombIncGammaStatErr0020->GetN()];
//     Double_t *sumErrorsCombDR0020                               = new Double_t[graphCombIncGammaStatErr0020->GetN()];
//     Double_t *StatErrorsCombDR0020                              = new Double_t[graphCombIncGammaStatErr0020->GetN()];
//     Double_t *xErrorsDR0020                                     = new Double_t[graphCombIncGammaStatErr0020->GetN()];
//     graphCombIncGammaSysErr0020->Print();
//     for (Int_t i = 0; i< graphCombDRPi0FitStatErr0020->GetN(); i++){
//         SystErrorsCombDR0020[i]                                 = graphCombDRPi0FitSysErr0020->GetEYhigh()[i]/graphCombDRPi0FitSysErr0020->GetY()[i] *100;
//         SystAErrorsCombDR0020[i]                                = graphCombDRPi0FitSysAErr0020->GetEYhigh()[i]/graphCombDRPi0FitSysAErr0020->GetY()[i] *100;
//         SystBErrorsCombDR0020[i]                                = graphCombDRPi0FitSysBErr0020->GetEYhigh()[i]/graphCombDRPi0FitSysBErr0020->GetY()[i] *100;
//         SystCErrorsCombDR0020[i]                                = graphCombDRPi0FitSysCErr0020->GetEYhigh()[i]/graphCombDRPi0FitSysCErr0020->GetY()[i] *100;
//         StatErrorsCombDR0020[i]                                 = graphCombDRPi0FitStatErr0020->GetEYhigh()[i]/graphCombDRPi0FitStatErr0020->GetY()[i] *100;
//         sumErrorsCombDR0020[i]                                  = graphCombDRPi0FitSumErr0020->GetEYhigh()[i]/graphCombDRPi0FitSumErr0020->GetY()[i] *100;
// //         cout << i << "\t" << graphCombDRPi0FitSysErr0020->GetY()[i] << "\t" << graphCombDRPi0FitSysErr0020->GetEYhigh()[i] << "\t" <<SystErrorsCombDR0020[i] << endl;
//     }
//     xErrorsDR0020                                               = graphCombDRPi0FitStatErr0020->GetX();
// 
//     cout << "here !!! \n\n" << endl;
//     graphCombDRPi0FitSumErr0020->Print();
    
    
//     //_______________________ copy inclusive photon spectra _____________________    
//     TH1D* histoCombErrorsForDRSum0020                           = new TH1D("histoCombErrorsForDRSum0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     TH1D* histoCombErrorsForDRStat0020                          = new TH1D("histoCombErrorsForDRStat0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     TH1D* histoCombErrorsForDRSys0020                           = new TH1D("histoCombErrorsForDRSys0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     TH1D* histoCombErrorsForDRSysA0020                          = new TH1D("histoCombErrorsForDRSysA0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     TH1D* histoCombErrorsForDRSysB0020                          = new TH1D("histoCombErrorsForDRSysB0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     TH1D* histoCombErrorsForDRSysC0020                          = new TH1D("histoCombErrorsForDRSysC0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
//     
//     for(Int_t i = 1; i<graphCombDRPi0FitStatErr0020->GetN()+1;i++){
//         cout<< i << "\t"<<xErrorsDR0020[i-1]<<"  "<<histoCombErrorsForDRSum0020->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum0020->GetBinWidth(i) <<endl;
//         Double_t binErrorSummed                                 = sumErrorsCombDR0020[i-1];
//         Double_t binErrorSyst                                   = SystErrorsCombDR0020[i-1];
//         Double_t binErrorSystA                                  = SystAErrorsCombDR0020[i-1];
//         Double_t binErrorSystB                                  = SystBErrorsCombDR0020[i-1];
//         Double_t binErrorSystC                                  = SystCErrorsCombDR0020[i-1];
//         Double_t binErrorStat                                   = StatErrorsCombDR0020[i-1];
//         Double_t DR                                             = graphCombDRPi0FitStatErr0020->GetY()[i-1];
// 
//         cout << DR << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
//         histoCombErrorsForDRSum0020->SetBinContent(i,DR);
//         histoCombErrorsForDRSys0020->SetBinContent(i,DR);
//         histoCombErrorsForDRSysA0020->SetBinContent(i,DR);
//         histoCombErrorsForDRSysB0020->SetBinContent(i,DR);
//         histoCombErrorsForDRSysC0020->SetBinContent(i,DR);
//         histoCombErrorsForDRStat0020->SetBinContent(i,DR);
//         histoCombErrorsForDRSum0020->SetBinError(i,(binErrorSummed/100)*DR);
//         histoCombErrorsForDRSys0020->SetBinError(i,(binErrorSyst/100)*DR);
//         histoCombErrorsForDRSysA0020->SetBinError(i,(binErrorSystA/100)*DR);
//         histoCombErrorsForDRSysB0020->SetBinError(i,(binErrorSystB/100)*DR);
//         histoCombErrorsForDRSysC0020->SetBinError(i,(binErrorSystC/100)*DR);
//         histoCombErrorsForDRStat0020->SetBinError(i,(binErrorStat/100)*DR);
//     }
// 
//     for(Int_t i = 1; i<histoCombErrorsForDRSum0020->GetNbinsX()+1;i++){
//         histoCombDirGammaSpectrumErrSum0020->SetBinContent(i+1,-1);
//         histoCombDirGammaSpectrumErrSys0020->SetBinContent(i+1,-1);
//         histoCombDirGammaSpectrumErrSysA0020->SetBinContent(i+1,-1);
//         histoCombDirGammaSpectrumErrSysB0020->SetBinContent(i+1,-1);
//         histoCombDirGammaSpectrumErrSysC0020->SetBinContent(i+1,-1);
//         histoCombDirGammaSpectrumErrStat0020->SetBinContent(i+1,-1);
// 
//         histoCombDirGammaSpectrumErrSum0020->SetBinError(i+1,0);
//         histoCombDirGammaSpectrumErrSys0020->SetBinError(i+1,0);
//         histoCombDirGammaSpectrumErrSysA0020->SetBinError(i+1,0);
//         histoCombDirGammaSpectrumErrSysB0020->SetBinError(i+1,0);
//         histoCombDirGammaSpectrumErrSysC0020->SetBinError(i+1,0);
//         histoCombDirGammaSpectrumErrStat0020->SetBinError(i+1,0);
//     }
// 
//     // get the binning of the direct photons from the DR
//     TH1D *hisoCombDirGammaSpecSysErr0020                        = new TH1D(*histoCombErrorsForDRSys0020);
//     TH1D *hisoCombDirGammaSpecSysAErr0020                       = new TH1D(*histoCombErrorsForDRSysA0020);
//     TH1D *hisoCombDirGammaSpecSysBErr0020                       = new TH1D(*histoCombErrorsForDRSysB0020);
//     TH1D *hisoCombDirGammaSpecSysCErr0020                       = new TH1D(*histoCombErrorsForDRSysC0020);
//     TH1D *hisoCombDirGammaSpecStatErr0020                       = new TH1D(*histoCombErrorsForDRStat0020);
//     TH1D *hisoCombDirGammaSpecSumErr0020                        = new TH1D(*histoCombErrorsForDRSum0020);
// 
//     for(Int_t i = 1; i<graphCombDRPi0FitStatErr0020->GetN()+1; i++){
//         // obtain common quantities
//         Double_t Rgamma                 = histoCombErrorsForDRSys0020->GetBinContent(i);
//         Double_t nIncGamma              = graphCombIncGammaStatErr0020->GetY()[i-1];
//         
//         // calculating Systematics graph
//         Double_t errRgamma              = histoCombErrorsForDRSys0020->GetBinError(i);
//         Double_t errNIncGam             = graphCombIncGammaSysErr0020->GetEYhigh()[i-1];
//         Double_t q1                     = 1 - 1/ Rgamma;
//         
//         Double_t q1Error                = errRgamma/(Rgamma*Rgamma);
//         Double_t content                = nIncGamma * ( 1 - 1/ Rgamma);
//         Double_t error                  = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
//         Double_t errDR                  = content - error;
//         hisoCombDirGammaSpecSysErr0020->SetBinError(i, error);
//         hisoCombDirGammaSpecSysErr0020->SetBinContent(i, content);
//         histoCombDirGammaSpectrumErrSys0020->SetBinContent(i, errDR);
// 
//         // calculating Systematics A graph
//         errRgamma                       = histoCombErrorsForDRSysA0020->GetBinError(i);
//         errNIncGam                      = graphCombIncGammaSysAErr0020->GetEYhigh()[i-1];
//         q1                              = 1 - 1/ Rgamma;
//         q1Error                         = errRgamma/(Rgamma*Rgamma);
//         content                         = nIncGamma * ( 1 - 1/ Rgamma);
//         error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
//         errDR                           = content - error;
//         hisoCombDirGammaSpecSysAErr0020->SetBinError(i, error);
//         hisoCombDirGammaSpecSysAErr0020->SetBinContent(i, content);
//         histoCombDirGammaSpectrumErrSysA0020->SetBinContent(i, errDR);
// 
//         // calculating Systematics B graph
//         errRgamma                       = histoCombErrorsForDRSysB0020->GetBinError(i);
//         errNIncGam                      = graphCombIncGammaSysBErr0020->GetEYhigh()[i-1];
//         q1                              = 1 - 1/ Rgamma;
//         q1Error                         = errRgamma/(Rgamma*Rgamma);
//         content                         = nIncGamma * ( 1 - 1/ Rgamma);
//         error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
//         errDR                           = content - error;
//         hisoCombDirGammaSpecSysBErr0020->SetBinError(i, error);
//         hisoCombDirGammaSpecSysBErr0020->SetBinContent(i, content);
//         histoCombDirGammaSpectrumErrSysB0020->SetBinContent(i, errDR);
// 
//         // calculating Systematics C graph
//         errRgamma                       = histoCombErrorsForDRSysC0020->GetBinError(i);
//         errNIncGam                      = graphCombIncGammaSysCErr0020->GetEYhigh()[i-1];
//         q1                              = 1 - 1/ Rgamma;
//         q1Error                         = errRgamma/(Rgamma*Rgamma);
//         content                         = nIncGamma * ( 1 - 1/ Rgamma);
//         error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
//         errDR                           = content - error;
//         hisoCombDirGammaSpecSysCErr0020->SetBinError(i, error);
//         hisoCombDirGammaSpecSysCErr0020->SetBinContent(i, content);
//         histoCombDirGammaSpectrumErrSysC0020->SetBinContent(i, errDR);        
//         
//         // calculating Stat graphs
//         errRgamma                       = histoCombErrorsForDRStat0020->GetBinError(i);
//         errNIncGam                      = graphCombIncGammaStatErr0020->GetEYhigh()[i-1];
//         q1                              = 1 - 1/ Rgamma;
//         q1Error                         = errRgamma/(Rgamma*Rgamma);
//         content                         = nIncGamma * ( 1 - 1/ Rgamma);
//         error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
//         errDR                           = content - error;
//         hisoCombDirGammaSpecStatErr0020->SetBinError(i, error);
//         hisoCombDirGammaSpecStatErr0020->SetBinContent(i, content);
//         histoCombDirGammaSpectrumErrStat0020->SetBinContent(i, errDR);
//         
//         // calculating summed error graphs
//         errRgamma                       = hisoCombDirGammaSpecSumErr0020->GetBinError(i);
//         errNIncGam                      = graphCombIncGammaSumErr0020->GetEYhigh()[i-1];
//         q1                              = 1 - 1/ Rgamma;
//         q1Error                         = errRgamma/(Rgamma*Rgamma);
//         content                         = nIncGamma * ( 1 - 1/ Rgamma);
//         error                           = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
//         errDR                           = content - error;
//         hisoCombDirGammaSpecSumErr0020->SetBinError(i, error);
//         hisoCombDirGammaSpecSumErr0020->SetBinContent(i, content);
//         histoCombDirGammaSpectrumErrSum0020->SetBinContent(i, errDR);
//     }
// 
//     // purely calculating points based on all Systematic errors
//     TGraphAsymmErrors *graphCombDirGammaSpectrumSystErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSys0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
//     if(graphCombDirGammaSpectrumSystErr0020)graphCombDirGammaSpectrumSystErr0020->SetName("graphCombDirGammaSpectrumSystErr0020");
//     if(graphCombDirGammaSpectrumSystErr0020)graphCombDirGammaSpectrumSystErr0020->Print();
//     if(graphCombDirGammaSpectrumSystErr0020)cout << "graph has been found" << endl;
//     // purely calculating points based on all Systematic errors A
//     TGraphAsymmErrors *graphCombDirGammaSpectrumSystAErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysA0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
//     if(graphCombDirGammaSpectrumSystAErr0020)graphCombDirGammaSpectrumSystAErr0020->SetName("graphCombDirGammaSpectrumSystAErr0020");
//     if(graphCombDirGammaSpectrumSystAErr0020)graphCombDirGammaSpectrumSystAErr0020->Print();
//     // purely calculating points based on all Systematic errors B
//     TGraphAsymmErrors *graphCombDirGammaSpectrumSystBErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysB0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
//     if(graphCombDirGammaSpectrumSystBErr0020)graphCombDirGammaSpectrumSystBErr0020->SetName("graphCombDirGammaSpectrumSystBErr0020");
//     if(graphCombDirGammaSpectrumSystBErr0020)graphCombDirGammaSpectrumSystBErr0020->Print();
//     // purely calculating points based on all Systematic errors C
//     TGraphAsymmErrors *graphCombDirGammaSpectrumSystCErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysC0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
//     if(graphCombDirGammaSpectrumSystCErr0020)graphCombDirGammaSpectrumSystCErr0020->SetName("graphCombDirGammaSpectrumSystCErr0020");
//     if(graphCombDirGammaSpectrumSystCErr0020)graphCombDirGammaSpectrumSystCErr0020->Print();
// 
//     // purely calculating points based on Statistical errors
//     TGraphAsymmErrors *graphCombDirGammaSpectrumStatErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrStat0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
//     if(graphCombDirGammaSpectrumStatErr0020)graphCombDirGammaSpectrumStatErr0020->SetName("graphCombDirGammaSpectrumStatErr0020");
//     // purely calculating points based on all Systematic + Statistical errors
//     TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
//     if(graphCombDirGammaSpectrumSumErr0020)graphCombDirGammaSpectrumSumErr0020->SetName("graphCombDirGammaSpectrumSumErr0020");
//     // calculate points above confidence level summed errors
//     TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr0020Confi = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum0020,hisoCombDirGammaSpecStatErr0020,2,0.5);
//     if(graphCombDirGammaSpectrumSumErr0020Confi)graphCombDirGammaSpectrumSumErr0020Confi->SetName("graphCombDirGammaSpectrumSumErr0020Confi");
//     // calculate arrows for points with 0, error summed
//     TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr0020Ar = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum0020,hisoCombDirGammaSpecStatErr0020,5,0.5);
//     if(graphCombDirGammaSpectrumSumErr0020Ar)graphCombDirGammaSpectrumSumErr0020Ar->SetName("graphCombDirGammaSpectrumSumErr0020Ar");
// 
// 
//     TF1* fitThermalGamma0020Sum                                 = FitObject("e","fitThermalGamma0020Sum","Photon",histoCombDirGammaSpectrumErrSum0020,0.9,2.1,NULL,"QNRMEX0+");
//     fileFinalResults << WriteParameterToFile(fitThermalGamma0020Sum)<< endl;       
//     TF1* fitThermalGamma0020Sum23                               = FitObject("e","fitThermalGamma0020Sum23","Photon",histoCombDirGammaSpectrumErrSum0020,0.9,2.3,NULL,"QNRMEX0+");
//     fileFinalResults << WriteParameterToFile(fitThermalGamma0020Sum23)<< endl;       
//     TF1* fitThermalGamma0020Stat                                = FitObject("e","fitThermalGamma0020Stat","Photon",graphCombDirGammaSpectrumStatErr0020,0.9,2.1,NULL,"QNRMEX0+");
//     fileFinalResults << WriteParameterToFile(fitThermalGamma0020Stat)<< endl;   
//     
//     TF1* fitThermalGamma0020Sys                                 = FitObject("e","fitThermalGamma0020Sys","Photon",graphCombDirGammaSpectrumSystErr0020,0.9,2.1,NULL,"QNRMEX0+");
//     fileFinalResults << WriteParameterToFile(fitThermalGamma0020Sys)<< endl;   
//     TF1* fitThermalGamma0020SysA                                = FitObject("e","fitThermalGamma0020SysA","Photon",graphCombDirGammaSpectrumSystAErr0020,0.9,2.1,NULL,"QNRMEX0+");
//     fileFinalResults << WriteParameterToFile(fitThermalGamma0020SysA)<< endl;   
//     TF1* fitThermalGamma0020SysB                                = FitObject("e","fitThermalGamma0020SysB","Photon",graphCombDirGammaSpectrumSystBErr0020,0.9,2.1,NULL,"QNRMEX0+");
//     fileFinalResults << WriteParameterToFile(fitThermalGamma0020SysB)<< endl;   
//     TF1* fitThermalGamma0020SysC                                = FitObject("e","fitThermalGamma0020SysC","Photon",graphCombDirGammaSpectrumSystCErr0020,0.9,2.1,NULL,"QNRMEX0+");
//     fileFinalResults << WriteParameterToFile(fitThermalGamma0020SysC)<< endl;   
//     
//     TGraphAsymmErrors* graphCombDirGammaSpectrumErrSum0020      = new TGraphAsymmErrors(histoCombDirGammaSpectrumErrSum0020);
//     TF1* fitFullDirGamma0020Sys                                 = FitObject("qcd","fitFullDirGamma0020Sys","Photon",graphCombDirGammaSpectrumSystErr0020,0.9,14,NULL,"QNRMEX0+");
//     fileFinalResults << WriteParameterToFile(fitFullDirGamma0020Sys)<< endl;       
//     TF1* fitFullDirGamma0020Stat                                = FitObject("qcd","fitFullDirGamma0020Stat","Photon",graphCombDirGammaSpectrumStatErr0020,0.9,14,NULL,"QNRMEX0+");
//     fileFinalResults << WriteParameterToFile(fitFullDirGamma0020Stat)<< endl;       
//     
//     TGraphAsymmErrors* graphRatioCombFitDirGammaStatErr0020         = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr0020->Clone();
//     TGraphAsymmErrors* graphRatioCombFitDirGammaSysErr0020          = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr0020->Clone();
//     graphRatioCombFitDirGammaStatErr0020                            = CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaStatErr0020, fitFullDirGamma0020Stat); 
//     graphRatioCombFitDirGammaSysErr0020                             = CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaSysErr0020, fitFullDirGamma0020Sys); 
// 
//     // Calculate thermal spectrum
//     TGraphAsymmErrors* graphCombThermalGammaSpectrumStatErr0020     = NULL;
//     TF1* fitPureThermalGamma0020Stat                                = NULL;
//     TGraphAsymmErrors* graphCombThermalGammaSpectrumSysErr0020      = NULL;
//     TGraphAsymmErrors* graphCombThermalGammaSpectrumSumErr0020Ar    = NULL;
//     TF1* fitFullThermalGamma0020Sys                                 = NULL;
//     TF1* fitFullThermalGamma0020Stat                                = NULL;
//     TGraphAsymmErrors* graphRatioCombFitThermalGammaStatErr0020     = NULL;
//     TGraphAsymmErrors* graphRatioCombFitThermalGammaSysErr0020      = NULL;
//     if (graphCombDirGammaSpectrumStatErr0020){
//         graphCombThermalGammaSpectrumStatErr0020                    = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumStatErr0020);
//         graphCombThermalGammaSpectrumStatErr0020->Print();
//         fitPureThermalGamma0020Stat                                 = FitObject("e","fitPureThermalGamma0020Stat","Photon",graphCombThermalGammaSpectrumStatErr0020,0.9,2.1,NULL,"QNRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitPureThermalGamma0020Stat)<< endl;   
//         fitFullThermalGamma0020Stat                                 = FitObject("qcd","fitFullDirGamma0020Stat","Photon",graphCombThermalGammaSpectrumStatErr0020,0.9,14,NULL,"QNRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitFullThermalGamma0020Stat)<< endl;       
//         graphRatioCombFitThermalGammaStatErr0020                    = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumStatErr0020->Clone();
//         graphRatioCombFitThermalGammaStatErr0020                    = CalculateGraphErrRatioToFit(graphRatioCombFitThermalGammaStatErr0020, fitFullThermalGamma0020Stat); 
//     }
//     if (graphCombDirGammaSpectrumSystErr0020){
//         graphCombThermalGammaSpectrumSysErr0020                     = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSystErr0020);
//         graphCombThermalGammaSpectrumSysErr0020->Print();
//         fitFullThermalGamma0020Sys                                  = FitObject("qcd","fitFullDirGamma0020Sys","Photon",graphCombThermalGammaSpectrumSysErr0020,0.9,14,NULL,"QNRMEX0+");
//         fileFinalResults << WriteParameterToFile(fitFullThermalGamma0020Sys)<< endl;       
//         graphRatioCombFitThermalGammaSysErr0020                     = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumSysErr0020->Clone();
//         graphRatioCombFitThermalGammaSysErr0020                     = CalculateGraphErrRatioToFit(graphRatioCombFitThermalGammaSysErr0020, fitFullThermalGamma0020Sys); 
//         
//     }
//     
//     if (graphCombDirGammaSpectrumSumErr0020Ar){
//         graphCombThermalGammaSpectrumSumErr0020Ar                   = SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSumErr0020Ar, newBinsComb, 20);
//     }
//     
//     // Calculate RAA
//     cout << endl << "Calculating RAA" << endl;
//     TGraphAsymmErrors* graphCombRAADirGammaStat0020                 = NULL;
//     TGraphAsymmErrors* graphCombRAADirGammaSys0020                  = NULL;
//     TGraphAsymmErrors* graphCombRAADirGammaSum0020Ar                = NULL;
//     if (graphCombDirGammaSpectrumStatErr0020)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumStatErr0020, &graphCombRAADirGammaStat0020);
//     if (graphCombDirGammaSpectrumSystErr0020)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSystErr0020, &graphCombRAADirGammaSys0020);
//     if (graphCombDirGammaSpectrumSumErr0020Ar)     CalcRaaWithTheoryFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSumErr0020Ar, &graphCombRAADirGammaSum0020Ar, newBinsComb, 20);
//     if (graphCombRAADirGammaStat0020) graphCombRAADirGammaStat0020->Print();
//     if (graphCombRAADirGammaSys0020) graphCombRAADirGammaSys0020->Print();
//     if (graphCombRAADirGammaSum0020Ar) graphCombRAADirGammaSum0020Ar->Print();

    
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
    
    //*******************************************************************************************************************************************
    //********************************************** DR plot with PCM pp & pPb measurement + NLO  ***********************************************
    //*******************************************************************************************************************************************    
    canvasRatioIndDR->cd();
    
    padPartRatioInDR1->Draw();
    padPartRatioInDR2->Draw();
    padPartRatioInDR3->Draw();

    
    //_______________________________________________________________ pPb panel _______________________________________________________________
    padPartRatioInDR1->cd();
    padPartRatioInDR1->SetLogx(1);
        dummyDR1->GetXaxis()->SetRangeUser(doubleRatioXpp[0],doubleRatioXpp[1]);
        dummyDR1->Draw("");
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
        
        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErrpPb, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErrpPb, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODRpPb, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        graphTheoryNLODRpPb->Draw("pE3lsame");
        graphPCMDRPi0FitSysErrpPb->Draw("E2same");
        histoPCMDRPi0FitStatErrpPb->Draw("p,same,e0,X0");
        
        TLegend* legendDRPCMNLOpPb = new TLegend(0.15,0.82-1.1*0.85*textsizeLabelsPad1*3,0.5,0.82);
        legendDRPCMNLOpPb->SetFillStyle(0);
        legendDRPCMNLOpPb->SetFillColor(0);
        legendDRPCMNLOpPb->SetLineColor(0);
        legendDRPCMNLOpPb->SetTextSize(0.85*textsizeLabelsPad1);
        legendDRPCMNLOpPb->SetMargin(0.2);
        legendDRPCMNLOpPb->SetTextFont(42);
        legendDRPCMNLOpPb->AddEntry(graphPCMDRPi0FitSysErrpPb,"PCM","pf");
        legendDRPCMNLOpPb->AddEntry(graphTheoryNLODRpPb,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
        legendDRPCMNLOpPb->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
        legendDRPCMNLOpPb->Draw();

        TLatex *labelDRCentpPb = new TLatex(0.12,0.85,collisionSystempPb.Data());
        SetStyleTLatex( labelDRCentpPb, 0.85*textsizeLabelsPad1,4);
        labelDRCentpPb->Draw();
        
    //_______________________________________________________________ pp 7TeV panel _______________________________________________________________
    padPartRatioInDR2->cd();
    padPartRatioInDR2->SetLogx(1);
        dummyDR2->GetXaxis()->SetRangeUser(doubleRatioXpp[0],doubleRatioXpp[1]);
        dummyDR2->Draw("");
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
// 
//         DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErrpp7TeV, markerStyleCombpp7TeV, markerSizeCombpp7TeV, colorCombpp7TeV , colorCombpp7TeV,widthLinesBoxes, kTRUE);
//         DrawGammaSetMarker(histoPCMDRPi0FitStatErrpp7TeV, markerStyleCombpp7TeV, markerSizeCombpp7TeV, colorCombpp7TeV , colorCombpp7TeV);
//         SetStyleGammaNLOTGraphWithBand( graphTheoryNLODRpp7TeV, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
//         graphTheoryNLODRpp7TeV->Draw("p3lsame");
//         graphPCMDRPi0FitSysErrpp7TeV->Draw("E2same");
//         histoPCMDRPi0FitStatErrpp7TeV->Draw("p,same,e0,X0");
// 
//         TLegend* legendDRPCMNLOpp7TeV = new TLegend(0.15,0.85-1.1*0.85*textsizeLabelsPad2*3,0.5,0.85);
//         legendDRPCMNLOpp7TeV->SetFillStyle(0);
//         legendDRPCMNLOpp7TeV->SetFillColor(0);
//         legendDRPCMNLOpp7TeV->SetLineColor(0);
//         legendDRPCMNLOpp7TeV->SetTextSize(0.85*textsizeLabelsPad2);
//         legendDRPCMNLOpp7TeV->SetMargin(0.2);
//         legendDRPCMNLOpp7TeV->SetTextFont(42);
//         legendDRPCMNLOpp7TeV->AddEntry(graphPCMDRPi0FitSysErrpp7TeV,"PCM","pf");
//         legendDRPCMNLOpp7TeV->AddEntry(graphTheoryNLODRpp7TeV,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
//         legendDRPCMNLOpp7TeV->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
//         legendDRPCMNLOpp7TeV->Draw();        
// 
//         TLatex *labelDRCentpp7TeV = new TLatex(0.12,0.88,collisionSystempp7TeV.Data());
//         SetStyleTLatex( labelDRCentpp7TeV, 0.85*textsizeLabelsPad2,4);        
//         labelDRCentpp7TeV->Draw();
    
    //_______________________________________________________________ pp 2.76TeV panel _______________________________________________________________    
    padPartRatioInDR3->cd();
    padPartRatioInDR3->SetLogx(1);
        dummyDR3->GetXaxis()->SetRangeUser(doubleRatioXpp[0],doubleRatioXpp[1]);
        dummyDR3->Draw("");
        DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);    
//         
//         DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErrpp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV,widthLinesBoxes, kTRUE);
//         DrawGammaSetMarker(histoPCMDRPi0FitStatErrpp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
//         SetStyleGammaNLOTGraphWithBand( graphTheoryNLODRpp2760GeV, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
//         graphTheoryNLODRpp2760GeV->Draw("p3lsame");
//         graphPCMDRPi0FitSysErrpp2760GeV->Draw("E2same");
//         histoPCMDRPi0FitStatErrpp2760GeV->Draw("p,same,e0,X0");
// 
//         TLegend* legendDRPCMNLOpp2760GeV = new TLegend(0.15,0.87-1.1*0.85*textsizeLabelsPad3*3,0.5,0.87);
//         legendDRPCMNLOpp2760GeV->SetFillStyle(0);
//         legendDRPCMNLOpp2760GeV->SetFillColor(0);
//         legendDRPCMNLOpp2760GeV->SetLineColor(0);
//         legendDRPCMNLOpp2760GeV->SetTextSize(0.85*textsizeLabelsPad3);
//         legendDRPCMNLOpp2760GeV->SetMargin(0.2);
//         legendDRPCMNLOpp2760GeV->SetTextFont(42);
//         legendDRPCMNLOpp2760GeV->AddEntry(graphPCMDRPi0FitSysErrpp2760GeV,"PCM","pf");
//         legendDRPCMNLOpp2760GeV->AddEntry(graphTheoryNLODRpp2760GeV,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
//         legendDRPCMNLOpp2760GeV->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
//         legendDRPCMNLOpp2760GeV->Draw();
// 
//         TLatex *labelDRCentpp2760GeV = new TLatex(0.12,0.9,collisionSystempp2760GeV.Data());
//         SetStyleTLatex( labelDRCentpp2760GeV, 0.85*textsizeLabelsPad3,4);
//         labelDRCentpp2760GeV->Draw();
        
    canvasRatioIndDR->SaveAs(Form("%s/DR_PCMMeasurementTheory_multCent_pPb.%s", outputDir.Data(), suffix.Data()));    

//     //*******************************************************************************************************************************************
//     //*************************************************** Plotting inclusive Gamma Spectrum *****************************************************
//     //*******************************************************************************************************************************************
//     
//     TGraphAsymmErrors* graphCombIncGammaSysErr0020Plot = ScaleGraph(graphCombIncGammaSysErr0020,100);
//     TGraphAsymmErrors* graphCombIncGammaStatErr0020Plot = ScaleGraph(graphCombIncGammaStatErr0020,100);
//     histoFitQCDIncGammaComb0020->Scale(100);
//     
//     TGraphAsymmErrors* graphCombIncGammaSysErr2040Plot = ScaleGraph(graphCombIncGammaSysErr2040,10);
//     TGraphAsymmErrors* graphCombIncGammaStatErr2040Plot = ScaleGraph(graphCombIncGammaStatErr2040,10);
//     histoFitQCDIncGammaComb2040->Scale(10);    
//     
//     TGraphAsymmErrors* graphCombIncGammaSysErr4080Plot = ScaleGraph(graphCombIncGammaSysErr4080,1);
//     TGraphAsymmErrors* graphCombIncGammaStatErr4080Plot = ScaleGraph(graphCombIncGammaStatErr4080,1);
//     histoFitQCDIncGammaComb4080->Scale(1);
//     TH1D* histoFitDummyPlotting = (TH1D*) histoFitQCDIncGammaComb4080->Clone("histoFitDummyPlotting");
//     
//     
//     TCanvas *canvasIncGamma = new TCanvas("canvasIncGamma","",10,10,1200,1400);  // gives the page size        
//     DrawGammaCanvasSettings( canvasIncGamma, 0.16, 0.01, 0.01, 0.07);
//     canvasIncGamma->SetLogy();
//     canvasIncGamma->SetLogx();
// 
//     Int_t textSizeLabelsPixelIncGam = 48;
//     Double_t textsizeLabelsIncGamma = 0;
//     if (canvasIncGamma->XtoPixel(canvasIncGamma->GetX2()) < canvasIncGamma->YtoPixel(canvasIncGamma->GetY1())){
//         textsizeLabelsIncGamma = (Double_t)textSizeLabelsPixelIncGam/canvasIncGamma->XtoPixel(canvasIncGamma->GetX2()) ;
//     } else {
//         textsizeLabelsIncGamma = (Double_t)textSizeLabelsPixelIncGam/canvasIncGamma->YtoPixel(canvasIncGamma->GetY1());
//     }
//         
//     TH2D *dummyGamma ;
//     dummyGamma = new TH2D("dummyGamma", "dummyGamma", 1000, 0., 22, 1000., 6e-8,8e3);
//     SetStyleHistoTH2ForGraphs( dummyGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{inc}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
//                             0.85*textsizeLabelsIncGamma, textsizeLabelsIncGamma, 0.85*textsizeLabelsIncGamma, textsizeLabelsIncGamma, 0.75, 1.6);
//     dummyGamma->GetXaxis()->SetLabelOffset(-0.015);
//     dummyGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
// //     dummyGamma->GetXaxis()->SetRangeUser(0,16);
//     dummyGamma->DrawCopy();
// 
//     DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSysErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020,widthLinesBoxes, kTRUE);
//     DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);    
//     graphCombIncGammaSysErr0020Plot->Draw("E2same");
//     graphCombIncGammaStatErr0020Plot->Draw("p,E1Z,same");
// 
//     DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSysErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
//     DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);    
//     graphCombIncGammaSysErr2040Plot->Draw("E2same");
//     graphCombIncGammaStatErr2040Plot->Draw("p,E1Z,same");
//     
//     DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSysErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
//     DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);    
//     graphCombIncGammaSysErr4080Plot->Draw("E2same");
//     graphCombIncGammaStatErr4080Plot->Draw("p,E1Z,same");
// 
//     TLatex *labelScalingIncGamma0020 = new TLatex(11.25,1.8E-3,"x 10^{2}");
//     SetStyleTLatex( labelScalingIncGamma0020, 0.85*textsizeLabelsIncGamma,4,colorComb0020,42,kFALSE);
//     labelScalingIncGamma0020->Draw();
//     
//     TLatex *labelScalingIncGamma2040 = new TLatex(11.25,1E-4,"x 10^{1}");
//     SetStyleTLatex( labelScalingIncGamma2040, 0.85*textsizeLabelsIncGamma,4,colorComb2040,42,kFALSE);
//     labelScalingIncGamma2040->Draw();
//     
//     TLatex *labelScalingIncGamma4080 = new TLatex(11.25,2E-6,"x 10^{0}");
//     SetStyleTLatex( labelScalingIncGamma4080, 0.85*textsizeLabelsIncGamma,4,colorComb4080,42,kFALSE);
//     labelScalingIncGamma4080->Draw();
//     
//     TLatex *labelIncGammaColl = new TLatex(0.6,0.91,collisionSystem.Data());
//     SetStyleTLatex( labelIncGammaColl, 0.85*textsizeLabelsIncGamma,4);
//     labelIncGammaColl->Draw();
//     SetStyleHisto(histoFitDummyPlotting, widthCommonFit, 5, kGray+1);
//     
//     TLegend* legendIncGamma = new TLegend(0.6,0.9-1.*0.85*textsizeLabelsIncGamma*3,0.9,0.9);
//     legendIncGamma->SetFillStyle(0);
//     legendIncGamma->SetFillColor(0);
//     legendIncGamma->SetLineColor(0);
//     legendIncGamma->SetTextSize(0.85*textsizeLabelsIncGamma);
//     legendIncGamma->SetMargin(0.2);
//     legendIncGamma->SetTextFont(42);
//     legendIncGamma->AddEntry(graphCombIncGammaSysErr0020Plot,"  0-20% ALICE","pf");
//     legendIncGamma->AddEntry(graphCombIncGammaSysErr2040Plot,"20-40% ALICE","pf");
//     legendIncGamma->AddEntry(graphCombIncGammaSysErr4080Plot,"40-80% ALICE","pf");
// //     legendIncGamma->AddEntry(histoFitDummyPlotting,"Fits to #gamma_{inc}","l");
//     legendIncGamma->Draw();
// 
//     canvasIncGamma->Print(Form("%s/IncGammaSpectrum.%s",outputDir.Data(),suffix.Data()));
// 
//     
    
    // double ratio combined
    TCanvas *canvasDoubleRatio = new TCanvas("canvasDoubleRatio","",0.095,0.09,1000,815);
    DrawGammaCanvasSettings( canvasDoubleRatio, 0.09, 0.02, 0.02, 0.09);
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
        
        DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErrpPb, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPHOSDRPi0FitStatErrpPb,  markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
        graphPHOSDRPi0FitSysErrpPb->Draw("E2same");
        histoPHOSDRPi0FitStatErrpPb->Draw("p,same,e0,X0");
        
//         DrawGammaSetMarkerTGraphAsym(graphEMCDRPi0FitSysErrpPb, markerStyleEMC, markerSizeEMC, colorEMC , colorEMC,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoEMCDRPi0FitStatErrpPb,  markerStyleEMC, markerSizeEMC, colorEMC , colorEMC);
//         graphEMCDRPi0FitSysErrpPb->Draw("E2same");
        histoEMCDRPi0FitStatErrpPb->Draw("p,same,e0,X0");
        
        DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErrpPb, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
        DrawGammaSetMarker(histoPCMDRPi0FitStatErrpPb,  markerStylePCM, markerSizePCM, colorPCM , colorPCM);
        SetStyleGammaNLOTGraphWithBand( graphTheoryNLODRpPb, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
        graphTheoryNLODRpPb->Draw("pE3lsame");
        graphPCMDRPi0FitSysErrpPb->Draw("E2same");
        histoPCMDRPi0FitStatErrpPb->Draw("p,same,e0,X0");
        
//      
        
        TLegend* legendDRPCMNLOpPbSingle = GetAndSetLegend2(0.15,0.90-1.05*0.85*textSizeSinglePad*5,0.5,0.90, 0.85*textSizeSinglePad, 1, "", 42, 0.2);
        legendDRPCMNLOpPbSingle->AddEntry(graphPCMDRPi0FitSysErrpPb,"PCM","pf");
        legendDRPCMNLOpPbSingle->AddEntry(graphPHOSDRPi0FitSysErrpPb,"PHOS","pf");
        legendDRPCMNLOpPbSingle->AddEntry(histoEMCDRPi0FitStatErrpPb,"EMC","p");
        legendDRPCMNLOpPbSingle->AddEntry(graphTheoryNLODRpPb,"NLO pQCD PDF: CT10 FF: GRV ","l");
        legendDRPCMNLOpPbSingle->Draw();

        TLatex *labelALICEDRSingle = new TLatex(0.95,0.91,"ALICE work in progress");
        SetStyleTLatex( labelALICEDRSingle, 0.85*textSizeSinglePad,4, 1, 42, kTRUE, 31); 
        labelALICEDRSingle->Draw();
       
        TLatex *labelDRCentpPbSingle = new TLatex(0.12,0.91,collisionSystempPb.Data());
        SetStyleTLatex( labelDRCentpPbSingle, 0.85*textSizeSinglePad,4, 1, 42, kTRUE, 11); 
        labelDRCentpPbSingle->Draw();
    hist2DDRDummySingle->Draw("same,axis");
        
    canvasDoubleRatio->Print(Form("%s/DR_PCMMeasurementTheory_pPb.%s", outputDir.Data(), suffix.Data()));    
    
    
}