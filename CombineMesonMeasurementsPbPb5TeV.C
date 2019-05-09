/****************************************************************************************************************************
******         provided by Gamma Conversion Group, PWG4,                                                     *****
******        Friederike Bock, friederike.bock@cern.ch                                                    *****
*****************************************************************************************************************************/

#include <Riostream.h>
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
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ExtractSignalBinning.h"
// #include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;

//____________________________________________________________________________________________________________________________________________
void CombineMesonMeasurementsPbPb5TeV(  TString centralityString        = "20-40%",
                                        TString fileNamePCM             = "/home/mike/1_PbPb_EMC/0_analysis/190227_PbPb_firstcombination/data_PCMResultsFullCorrection_PbPb_EtaOldBinning.root",
                                        TString fileNameEMCAL           = "/home/mike/1_PbPb_EMC/0_analysis/190311_EDC_new/pdf/PbPb_5.02TeV/2019_04_18/data_EMCAL-EMCALResultsFullCorrection_PbPb.root",
                                        TString fileNamePHOS            = "/home/mike/1_PbPb_EMC/0_analysis/190227_PbPb_firstcombination/20190228_PHOS_PbPb_5TeV.root",
                                        TString fileNamePCMEMCAL        = "/home/mike/1_PbPb_EMC/0_analysis/190311_PCMEDC_new/pdf/PbPb_5.02TeV/2019_04_25/data_PCM-EMCALResultsFullCorrection_PbPb.root",
                                        TString fileNamePCMPHOS         = "/home/mike/1_PbPb_EMC/0_analysis/190318_PCMPHOS_new/pdf/PbPb_5.02TeV/2019_04_25/data_PCM-PHOSResultsFullCorrection_PbPb.root",
                                        TString suffix                  = "pdf",
                                        TString fileInputCorrFactors    = "",
                                        Bool_t  enableEta               = kTRUE
                                    ){

    TString date = ReturnDateString();
    cout << "date: " << date << endl;
    Bool_t doOutput                                 = kTRUE;

    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();
    SetPlotStyle();
    
    gStyle->SetEndErrorSize(0);
    
    Bool_t plotMassAndWidth = kTRUE;
    Bool_t plotInvMassBins = kTRUE;
    Bool_t saveResults = kFALSE;

    Bool_t havePCM = kFALSE;
    Bool_t haveEMCAL = kFALSE;
    Bool_t havePHOS = kFALSE;
    Bool_t havePCMEMCAL = kFALSE;
    Bool_t havePCMPHOS = kFALSE;

    if(fileNamePCM.CompareTo("") != 0)         havePCM = kTRUE;
    if(fileNameEMCAL.CompareTo("") != 0)       haveEMCAL = kTRUE;
    if(fileNamePHOS.CompareTo("") != 0)        havePHOS = kTRUE;
    if(fileNamePCMEMCAL.CompareTo("") != 0)    havePCMEMCAL = kTRUE;
    if(fileNamePCMPHOS.CompareTo("") != 0)     havePCMPHOS = kTRUE;

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << "date for output: " << dateForOutput << endl;

    TString CentNameShort = "";
    if(centralityString.CompareTo("0-10%") == 0 ) CentNameShort  = "0_10";
    if(centralityString.CompareTo("10-20%") == 0 ) CentNameShort = "10_20";
    if(centralityString.CompareTo("20-40%") == 0 ) CentNameShort = "20_40";
    if(centralityString.CompareTo("40-60%") == 0 ) CentNameShort = "40_60";
    if(centralityString.CompareTo("60-80%") == 0 ) CentNameShort = "60_80";

    TString collisionSystemPbPb5TeV;
    collisionSystemPbPb5TeV.Form("PbPb, #sqrt{#it{s}_{NN}} = 5.02 TeV, %s", centralityString.Data());

    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements5TeVPbPb/%s",suffix.Data(),dateForOutput.Data(),CentNameShort.Data());
    cout << "output directory: "<< outputDir.Data() << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    if(havePCM)               gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    if(haveEMCAL)             gSystem->Exec(Form("cp %s %s/InputEMCAL.root", fileNameEMCAL.Data(), outputDir.Data()));
    if(havePHOS)              gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    if(havePCMEMCAL)          gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDir.Data()));
    if(havePCMPHOS)           gSystem->Exec(Form("cp %s %s/InputPCMPHOS.root", fileNamePCMPHOS.Data(), outputDir.Data()));
    
    TFile *inputpp        = new TFile("/home/mike/2_pp_EMC/0_analysis/190202_pp_combination/CombineMesonMeasurements5TeV_V2_april15/CombinedResultsPaperPP5TeV_2019_04_15.root");
    TDirectoryFile* inputppDir = (TDirectoryFile*)inputpp->Get("Pi05TeV");
    TDirectoryFile* inputppDirEta = (TDirectoryFile*)inputpp->Get("Eta5TeV");
    TF1* ppFit = (TF1*)inputppDir->Get("TsallisFitPi0");
    TF1* ppFitEta = (TF1*)inputppDirEta->Get("TsallisFitEta");

    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    
    Int_t textSizeLabelsPixel                   = 0;
    Double_t textsizeLabelsPP                   = 0.07;
 
    Width_t  widthLinesBoxes                    = 2;
    Width_t  widthCommonFit                     = 2;
    
    // Definition of colors, styles and markers sizes
    Color_t  colorComb                          = kGray+2;
    Style_t  markerStyleComb                    = 20;
    Size_t   markerSizeComb                     = 2;
    
    Color_t  colorCombLowPt                     = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kFALSE);
    Color_t  colorCombHighPt                    = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kTRUE);
    Style_t  markerStyleCombLowPt               = 20;
    Style_t  markerStyleCombHighPt              = 20;
    Size_t   markerSizeComparison               = 2;
    
    Color_t colorTrigg      [10]                = {kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2};
    Color_t colorTriggShade [10]                = {kGray+1, kGray, kRed-6, kBlue-6, kGreen-8, kCyan-6, kViolet-8, kMagenta-8,  kRed-8, kBlue-8};
    Marker_t markerTrigg    [10]                = {20, 20, 21, 34, 29, 33, 21, 27, 28, 30 };
    Marker_t markerTriggMC  [10]                = {24, 24, 25, 28, 30, 27, 25, 27, 28, 30 };

    Size_t sizeTrigg        [10]                = {1.5, 1.5, 1.5, 2, 2.2, 2., 1.5, 2., 2.5, 1.5 };

    Int_t totalNSets                            = 5;
    TString  nameMeasGlobal[5]                  = { "PCM", "EMCal", "PHOS", "PCM-EMCal", "PCM-PHOS"};
    
    Color_t  colorDet[5];
    Color_t  colorDetMC[5];
    Style_t  markerStyleDet[5];
    Style_t  markerStyleDetMC[5];
    Size_t   markerSizeDet[5];
    Size_t   markerSizeDetMC[5];

    
    for (Int_t i = 0; i < totalNSets; i++){
        colorDet[i]                             = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
        colorDetMC[i]                           = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
        markerStyleDet[i]                       = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
        markerStyleDetMC[i]                     = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
        markerSizeDet[i]                        = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
        markerSizeDetMC[i]                      = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
    }
    
    cout << "loaded plotting styles" << endl;
    TFile* inputFile[5]                            =  {NULL,NULL,NULL,NULL,NULL};
    inputFile[0]                                = new TFile(fileNamePCM.Data());
    inputFile[1]                                = new TFile(fileNameEMCAL.Data());
    inputFile[2]                                = new TFile(fileNamePHOS.Data());
    inputFile[3]                                = new TFile(fileNamePCMEMCAL.Data());
    inputFile[4]                                = new TFile(fileNamePCMPHOS.Data());

    TDirectory* directoryPi0[5]                    =  {NULL,NULL,NULL,NULL,NULL};
    for(Int_t i=0;i<totalNSets;i++){
      if(inputFile[i] && i!=2){
        if(!inputFile[i]->IsZombie()){
          cout << "loading directories for " <<  nameMeasGlobal[i] << endl;
          directoryPi0[i]                           = (TDirectory*)inputFile[i]->Get(Form("Pi0%sPbPb_5.02TeV_V0M",centralityString.Data()));
        }
      }
      if(inputFile[2] && i==2){
        cout << "loading directories for " <<  nameMeasGlobal[i] << endl;
        directoryPi0[i]                           = (TDirectory*)inputFile[i]->Get(Form("Pi05TeV_%s_V0M",CentNameShort.Data()));
      }
    }
    
    TDirectory* directoryEta[5]                    =  {NULL,NULL,NULL,NULL,NULL};
    for(Int_t i=0;i<totalNSets;i++){
      if(inputFile[i] && i!=2){
        if(!inputFile[i]->IsZombie()){
          cout << "loading directories for " <<  nameMeasGlobal[i] << endl;
          directoryEta[i]                           = (TDirectory*)inputFile[i]->Get(Form("Eta%sPbPb_5.02TeV_V0M",centralityString.Data()));
        }
      }
      if(inputFile[2] && i==2){
        cout << "loading directories for " <<  nameMeasGlobal[i] << endl;
        directoryEta[i]                           = (TDirectory*)inputFile[i]->Get(Form("Eta5TeV_%s_V0M",CentNameShort.Data()));
      }
    }

  cout << __LINE__<<endl;
  TH1D* histoPi0Mass[5]                                 = {NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0FWHMMeV[5]                              = {NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0TrueMass[5]                             = {NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0TrueFWHMMeV[5]                          = {NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0Acc[5]                                  = {NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0TrueEffPt[5]                            = {NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0AccTimesEff[5]                          = {NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0InvYield[5]                             = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0InvYieldStat[5]            = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0InvYieldSys[5]             = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphCombPi0Weights[5]             = {NULL,NULL,NULL,NULL,NULL};

  TH1D* histoEtaInvYield[5]                             = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaInvYieldStat[5]            = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaInvYieldSys[5]             = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphCombEtaWeights[5]             = {NULL,NULL,NULL,NULL,NULL};

  TH1D* histoEtaToPi0[5]                             = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaToPi0Stat[5]            = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaToPi0Sys[5]             = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphCombEtaToPi0Weights[5]     = {NULL,NULL,NULL,NULL,NULL};

  Double_t rapidityMeas[5]                       = {1.6, 1.6, 1., 1.6, 1.6};

// *******************************************************************************************************
// ************************** read input files ***********************************************************
// *******************************************************************************************************

  for (Int_t i = 0; i < totalNSets; i++){
    cout << "reading from " << nameMeasGlobal[i].Data() << " measurement..." << endl;
    if(directoryPi0[i]){
      if(i!=2){
      cout << "loading pi0 inputs" << endl;
      // load mass/width/effi plots
      histoPi0Mass[i]                             = (TH1D*)directoryPi0[i]->Get("Pi0_Mass_data_INT7");
      histoPi0FWHMMeV[i]                          = (TH1D*)directoryPi0[i]->Get("Pi0_Width_data_INT7");
      histoPi0TrueMass[i]                         = (TH1D*)directoryPi0[i]->Get("Pi0_Mass_MC_INT7");
      histoPi0TrueFWHMMeV[i]                      = (TH1D*)directoryPi0[i]->Get("Pi0_Width_MC_INT7");

      if(histoPi0Mass[i]) histoPi0Mass[i]->Scale(1000.);
      if(histoPi0FWHMMeV[i]) histoPi0FWHMMeV[i]->Scale(1000.);
      if(histoPi0TrueMass[i]) histoPi0TrueMass[i]->Scale(1000.);
      if(histoPi0TrueFWHMMeV[i]) histoPi0TrueFWHMMeV[i]->Scale(1000.);
      histoPi0Acc[i]                              = (TH1D*)directoryPi0[i]->Get("AcceptancePi0_INT7");
      histoPi0TrueEffPt[i]                        = (TH1D*)directoryPi0[i]->Get("EfficiencyPi0_INT7");
      if(histoPi0TrueEffPt[i] && histoPi0Acc[i]){
        histoPi0AccTimesEff[i]                      = (TH1D*)histoPi0TrueEffPt[i]->Clone(Form("histoPi0AccTimesEff%s",nameMeasGlobal[i].Data()));
        histoPi0AccTimesEff[i]->Multiply(histoPi0Acc[i]);
        histoPi0AccTimesEff[i]->Scale(2*TMath::Pi()*rapidityMeas[i]);
      }

      // load cross section systematics and datapoints
      histoPi0InvYield[i]                      = (TH1D*)directoryPi0[i]->Get("CorrectedYieldPi0");
      graphPi0InvYieldStat[i]                  = new TGraphAsymmErrors(histoPi0InvYield[i]);
      graphPi0InvYieldSys[i]                   = (TGraphAsymmErrors*)directoryPi0[i]->Get("Pi0SystError");
      
      if(enableEta){
        histoEtaInvYield[i]                      = (TH1D*)directoryEta[i]->Get("CorrectedYieldEta");
        graphEtaInvYieldStat[i]                  = new TGraphAsymmErrors(histoEtaInvYield[i]);
        graphEtaInvYieldSys[i]                   = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaSystError");

        histoEtaToPi0[i]                      = (TH1D*)directoryEta[i]->Get("EtaToPi0StatError");
        graphEtaToPi0Stat[i]                  = new TGraphAsymmErrors(histoEtaToPi0[i]);
        graphEtaToPi0Sys[i]                   = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0SystError");
      }

      cout << nameMeasGlobal[i].Data() << " pi0 stat:" << graphPi0InvYieldStat[i] << endl;
      if(doOutput) graphPi0InvYieldStat[i]->Print();
      cout << nameMeasGlobal[i].Data() << " pi0 sys:" << graphPi0InvYieldSys[i] << endl;
      if(doOutput) graphPi0InvYieldSys[i]->Print();
      }
      if(i==2){
        histoPi0AccTimesEff[2]                   = (TH1D*)directoryPi0[2]->Get("EffTimesAcc");
        histoPi0Mass[2]                          = (TH1D*)directoryPi0[2]->Get("Pi0_Mass_data");
        histoPi0FWHMMeV[2]                       = (TH1D*)directoryPi0[2]->Get("Pi0_Width_data");
        histoPi0TrueMass[2]                      = (TH1D*)directoryPi0[2]->Get("Pi0_Mass_MC");
        histoPi0TrueFWHMMeV[2]                   = (TH1D*)directoryPi0[2]->Get("Pi0_Width_MC");
        if(histoPi0Mass[2]) histoPi0Mass[2]->Scale(1000.);
        if(histoPi0FWHMMeV[2]) histoPi0FWHMMeV[2]->Scale(1000.);
        if(histoPi0TrueMass[2]) histoPi0TrueMass[2]->Scale(1000.);
        if(histoPi0TrueFWHMMeV[2]) histoPi0TrueFWHMMeV[2]->Scale(1000.);
        histoPi0AccTimesEff[2]->Scale(2*TMath::Pi()*rapidityMeas[i]);
        
        histoPi0InvYield[2]                      = (TH1D*)directoryPi0[2]->Get("CorrectedYieldPi0");
        graphPi0InvYieldStat[2]                  = new TGraphAsymmErrors(histoPi0InvYield[2]);
        graphPi0InvYieldSys[2]                   = (TGraphAsymmErrors*)directoryPi0[i]->Get("CorrectedYieldPi0Sys");

        if(enableEta){
          histoEtaInvYield[i]                      = (TH1D*)directoryEta[i]->Get("CorrectedYieldEta");
          graphEtaInvYieldStat[i]                  = new TGraphAsymmErrors(histoEtaInvYield[i]);
          graphEtaInvYieldSys[i]                   = (TGraphAsymmErrors*)directoryEta[i]->Get("CorrectedYieldEtaSys");

          histoEtaToPi0[i]                      = (TH1D*)directoryEta[i]->Get("EtaToPi0StatError");
          graphEtaToPi0Stat[i]                  = new TGraphAsymmErrors(histoEtaToPi0[i]);
          graphEtaToPi0Sys[i]                   = (TGraphAsymmErrors*)directoryEta[i]->Get("graphEtaToPi0SystError");
        }
        cout << nameMeasGlobal[i].Data() << " pi0 stat:" << graphPi0InvYieldStat[i] << endl;
        if(doOutput) graphPi0InvYieldStat[i]->Print();
        cout << nameMeasGlobal[i].Data() << " pi0 sys:" << graphPi0InvYieldSys[i] << endl;
        if(doOutput) graphPi0InvYieldSys[i]->Print();
      }
    }
  }
  
  Double_t minPtPi0                               = 0.4;
  Double_t maxPtPi0                               = 35.0; //75
  Double_t delaPtPi0                              = maxPtPi0 - minPtPi0;
  Double_t prodPtPi0                              = maxPtPi0 * minPtPi0;
  Double_t minYieldPi0                            = 10e-10;
  Double_t maxYieldPi0                            = 10e2;
  
  Double_t minPtEta                               = 0.4;
  Double_t maxPtEta                               = 30.0; //75
  Double_t delaPtEta                              = maxPtPi0 - minPtPi0;
  Double_t prodPtEta                              = maxPtPi0 * minPtPi0;
  Double_t minYieldEta                            = 10e-10;
  Double_t maxYieldEta                            = 10e2;

  Double_t pTPi0PbPb5TeV[35]                =  {  0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                                                          2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                                          4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10., 12., 14.,
                                                          16., 20., 25., 30., 35.};

  Double_t pTEtaPbPb5TeV[14]              = { 0.4, 0.8, 1.4, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 25.0, 30.0};

  const Int_t nBinsPi0 = 34;
  const Int_t nBinsEta = 13;
  //                      {"PCM","EMCal","PHOS" ,"PCM-EMCal","PCM-PHOS"};
  //offsets Pi0
  Double_t startpT[5]   = {0.4  ,1.4    ,0.4    ,1.0        ,0.8};
  Double_t endpT[5]     = {12.  ,20.    ,30.    ,20.        ,12.};
  if(centralityString.CompareTo("0-10%") == 0 )  {startpT[1] = 6.0; startpT[3] = 2.2; endpT[2] = 35.;}
  if(centralityString.CompareTo("10-20%") == 0 ) {startpT[1] = 4.0; startpT[3] = 2.2; endpT[2] = 35.;}
  if(centralityString.CompareTo("40-60%") == 0 ) {endpT[2] = 20.;}
  if(centralityString.CompareTo("60-80%") == 0 ) {endpT[2] = 20.;}

  Int_t      offSetMethod[5] = {1,6,1,4,3};
  if(centralityString.CompareTo("0-10%") == 0 )  {offSetMethod[1] = 14; offSetMethod[3] = 10;}
  if(centralityString.CompareTo("10-20%") == 0 ) {offSetMethod[1] = 9; offSetMethod[3] = 10;}
  if(centralityString.CompareTo("40-60%") == 0 ) {}
  if(centralityString.CompareTo("60-80%") == 0 ) {}

  //offsets Eta
  Double_t startpTEta[5]   = {1.4  ,2.0    ,2.0    ,1.4        ,1.4};
  Double_t endpTEta[5]     = {10.  ,20.    ,20.    ,16.        ,8.};
  if(centralityString.CompareTo("0-10%") == 0 )  { startpTEta[0] = 0.8; startpTEta[1] = 4.0; startpTEta[2] = 3.0; startpTEta[3] = 4.0; endpTEta[3] = 10.0;  endpTEta[4] = 6.0;}
  if(centralityString.CompareTo("10-20%") == 0 ) { startpTEta[1] = 3.0; startpTEta[2] = 3.0; startpTEta[3] = 2.0; endpTEta[3] = 10.0;}
  if(centralityString.CompareTo("40-60%") == 0 ) { startpTEta[0] = 0.8; endpTEta[2] = 16.0;  }
  if(centralityString.CompareTo("60-80%") == 0 ) { startpTEta[0] = 0.8; endpTEta[0] = 8.0; endpTEta[2] = 12.0;  endpTEta[4] = 6.0;}
  Int_t    offSetMethodEta[5] = {1,3,3,2,2};
  Int_t    offSetStatPHOS = 3;
  if(centralityString.CompareTo("0-10%") == 0 || centralityString.CompareTo("10-20%") == 0)  { offSetStatPHOS = 4; offSetMethodEta[2] = 4;}
  if(centralityString.CompareTo("0-10%") == 0 )  {offSetMethodEta[1] = 5; offSetMethodEta[3] = 5;}
  if(centralityString.CompareTo("10-20%") == 0 ) {offSetMethodEta[1] = 4; offSetMethodEta[3] = 3;}
  if(centralityString.CompareTo("40-60%") == 0 ) {}
  if(centralityString.CompareTo("60-80%") == 0 ) {}

  // **********************************************************************************************************************
  // ******************************************* Calculate simple combination      ****************************************
  // **********************************************************************************************************************
  cout << endl;
  cout << "*******************************************" << endl;
  cout << "*******************************************" << endl;
  cout << "*****************  Pi0  *******************" << endl;
  cout << "*******************************************" << endl;
  cout << "*******************************************" << endl;
  cout << endl;

  TGraphAsymmErrors* graphCombPi0InvYieldStat  = NULL;
  TGraphAsymmErrors* graphCombPi0InvYieldSys   = NULL;
  TGraphAsymmErrors* graphCombPi0InvYieldStat_Norm  = NULL;
  TGraphAsymmErrors* graphCombPi0InvYieldSys_Norm   = NULL;
  TGraphAsymmErrors* graphCombPi0RAAStat  = NULL;
  TGraphAsymmErrors* graphCombPi0RAASys   = NULL;
  
  Double_t  pTPi0BinCenters[nBinsPi0];
  Double_t  pTPi0BinWidths[nBinsPi0];
  Double_t  y_Comb[nBinsPi0];
  Double_t  y_E_stat_Comb[nBinsPi0];
  Double_t  y_E_syst_Comb[nBinsPi0];
  
  Double_t  y_Comb_Norm[nBinsPi0];
  Double_t  y_E_stat_Comb_Norm[nBinsPi0];
  Double_t  y_E_syst_Comb_Norm[nBinsPi0];
  Double_t  y_Comb_Weights[5][nBinsPi0];
  
  Double_t temp_y[5] = {0,0,0,0,0};
  Double_t temp_y_E_stat[5] = {0,0,0,0,0};
  Double_t temp_y_E_syst[5] = {0,0,0,0,0};
  Double_t temp_y_comb;
  Double_t temp_y_comb_E_tot;
  Int_t    temp_Ny = 0;
  Double_t EStat0, EStat1, EStat2, EStat3, EStat4;
  
  Int_t      n_Meas[5];
  Double_t*  x_Meas[5];
  Double_t*  x_Meas_syst[5];
  Double_t*  y_Meas[5];
  Double_t*  y_E_stat_Meas[5];
  Double_t*  y_E_syst_Meas[5];
  for (Int_t j = 0; j < totalNSets; j++){
      if(directoryPi0[j] && graphPi0InvYieldStat[j]){
        cout << nameMeasGlobal[j].Data() << "\t\t getting values" << endl;
        n_Meas[j] = graphPi0InvYieldStat[j]->GetN();
        x_Meas[j] = graphPi0InvYieldStat[j]->GetX();
        x_Meas_syst[j] = graphPi0InvYieldSys[j]->GetX();
        y_Meas[j] = graphPi0InvYieldStat[j]->GetY();
        y_E_stat_Meas[j] = graphPi0InvYieldStat[j]->GetEYlow();
        y_E_syst_Meas[j] = graphPi0InvYieldSys[j] ->GetEYlow();
      }
  }
  
  cout << "*******************************************" << endl;
  cout << "printing the arrays" << endl;
  cout << "*******************************************" << endl;
  for (Int_t j = 0; j < totalNSets; j++){
    if(directoryPi0[j] && graphPi0InvYieldStat[j]){
      cout << nameMeasGlobal[j].Data() << "\t" << n_Meas[j] << endl;
      cout << "first x value: " << x_Meas[j][0] << endl;
      cout << "first x value with systematics: " << x_Meas_syst[j][0] << endl;
    }
//     if(havePHOS && j==2 ){
//       cout << nameMeasGlobal[j].Data() << "\t" << n_Meas[j] << endl;
//       cout << "first x value: " << x_Meas[j][0] << endl;
//       cout << "first x value with systematics: " << x_Meas_syst[j][0] << endl;
//     }
  }

  cout << "*******************************************" << endl;
  cout << "PHOS syst" << endl;
  cout << "*******************************************" << endl;
  for(Int_t i=0; i<n_Meas[2]; i++){
    cout << y_E_syst_Meas[2][i]/y_Meas[2][i] << endl;
  }
  cout << "*******************************************" << endl;
  cout << "starting to calculate a simple combination of measurements" << endl;
  cout << "*******************************************" << endl;
  for(Int_t i=0; i<nBinsPi0; i++){
    pTPi0BinCenters[i] = ( pTPi0PbPb5TeV[i] + pTPi0PbPb5TeV[i+1] ) / 2. ;
    pTPi0BinWidths[i] = ( pTPi0PbPb5TeV[i+1] - pTPi0PbPb5TeV[i] ) / 2. ;
    cout << endl;
    cout << "pT bin : " << pTPi0PbPb5TeV[i] << "-" << pTPi0PbPb5TeV[i+1] << "\t center :" << pTPi0BinCenters[i] << endl;
    temp_Ny = 0;
    for (Int_t j = 0; j < totalNSets; j++){
      temp_y[j] = 0;
      temp_y_E_stat[j] = 0;
      temp_y_E_syst[j] = 0;
      if(directoryPi0[j] && graphPi0InvYieldStat[j] && j!=2 && (pTPi0BinCenters[i]>startpT[j] && pTPi0BinCenters[i]<endpT[j]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        temp_y[j] = y_Meas[j][i];
        temp_y_E_stat[j] = y_E_stat_Meas[j][i];
        temp_y_E_syst[j] = y_E_syst_Meas[j][i-offSetMethod[j]];
//         if(j==0) temp_y_E_syst[j] = temp_y[j] * 0.1;
        cout << temp_y[j] << endl;
        temp_Ny++;
      }
      if(havePHOS && j==2 && (pTPi0BinCenters[i]>startpT[2] && pTPi0BinCenters[i]<endpT[2]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        temp_y[j] = y_Meas[j][i-1];
        temp_y_E_stat[j] = y_E_stat_Meas[j][i-1];
        temp_y_E_syst[j] = y_E_syst_Meas[j][i-offSetMethod[j]];
        cout << temp_y[j] << endl;
        temp_Ny++;
      }
    }
    cout << "This pt bin has " << temp_Ny << " points to average " << endl;
    temp_y_comb      = 0;
    temp_y_comb_E_tot = 0;
    for (Int_t j = 0; j < totalNSets; j++){
      y_Comb_Weights[j][i] = 0;
      if(temp_y_E_stat[j]>0 && temp_y_E_syst[j]>0){
        temp_y_comb += ( temp_y[j] / (pow(temp_y_E_stat[j],2) + pow(temp_y_E_syst[j],2) ) );
        temp_y_comb_E_tot += ( 1 / (pow(temp_y_E_stat[j],2) + pow(temp_y_E_syst[j],2) ) );
        y_Comb_Weights[j][i] = ( 1 / (pow(temp_y_E_stat[j],2) + pow(temp_y_E_syst[j],2) ) );
      }
    }
    if(temp_y_comb_E_tot>0){
      y_Comb[i]         = temp_y_comb / temp_y_comb_E_tot ;
      y_Comb_Norm[i]    = 1 ;
      for (Int_t j = 0; j < totalNSets; j++){
        if(pTPi0BinCenters[i]>startpT[j] && pTPi0BinCenters[i]<endpT[j]){
          y_Comb_Weights[j][i] /= temp_y_comb_E_tot;
        } else {
          y_Comb_Weights[j][i] = 0.;
        }
      }
    } else {
      y_Comb[i]         = 1e-40;
      y_Comb_Norm[i]    = -10. ;
      for (Int_t j = 0; j < totalNSets; j++){
        y_Comb_Weights[j][i] = 0. ;
      }
    }
    //go fully correlated -> if they are uncorrelated, take RMS
//     temp_Ny = 1.;
    cout << "y_Comb = \t " << y_Comb[i] << endl;
    EStat0 = 0; EStat1 = 0; EStat2 = 0; EStat3 = 0; EStat4 = 0;
    if(temp_y_E_stat[0]>0) EStat0 = pow(1./temp_y_E_stat[0],2);
    if(temp_y_E_stat[1]>0) EStat1 = pow(1./temp_y_E_stat[1],2);
    if(temp_y_E_stat[2]>0) EStat2 = pow(1./temp_y_E_stat[2],2);
    if(temp_y_E_stat[3]>0) EStat3 = pow(1./temp_y_E_stat[3],2);
    if(temp_y_E_stat[4]>0) EStat4 = pow(1./temp_y_E_stat[4],2);
    y_E_stat_Comb[i]  = 1./sqrt(EStat0+EStat1+EStat2+EStat3+EStat4);
    y_E_syst_Comb[i]  = sqrt(((pow(y_Comb_Weights[0][i]*temp_y_E_syst[0],2)+pow(y_Comb_Weights[1][i]*temp_y_E_syst[1],2)+pow(y_Comb_Weights[2][i]*temp_y_E_syst[2],2)+pow(y_Comb_Weights[3][i]*temp_y_E_syst[3],2)+pow(y_Comb_Weights[4][i]*temp_y_E_syst[4],2))/temp_Ny)/(pow(y_Comb_Weights[0][i],2)+pow(y_Comb_Weights[1][i],2)+pow(y_Comb_Weights[2][i],2)+pow(y_Comb_Weights[3][i],2)+pow(y_Comb_Weights[4][i],2)));
    if(y_Comb[i] == 1e-40){
      y_E_stat_Comb[i]  = 0;
      y_E_syst_Comb[i]  = 0;
    }
    y_E_stat_Comb_Norm[i]  = y_E_stat_Comb[i] / y_Comb[i];
    y_E_syst_Comb_Norm[i]  = y_E_syst_Comb[i] / y_Comb[i];
    cout << y_E_stat_Comb[i]  << "\t" << pow(temp_y_E_stat[0],2) << "\t" << pow(temp_y_E_stat[1],2) << "\t" << pow(temp_y_E_stat[2],2) << "\t" << pow(temp_y_E_stat[3],2) << "\t" << pow(temp_y_E_stat[4],2) << endl;
    cout << y_E_syst_Comb[i]  << "\t" << pow(temp_y_E_syst[0],2) << "\t" << pow(temp_y_E_syst[1],2) << "\t" << pow(temp_y_E_syst[2],2) << "\t" << pow(temp_y_E_syst[3],2) << "\t" << pow(temp_y_E_syst[4],2) << endl;
    for (Int_t j = 0; j < totalNSets; j++){
      if( y_Comb_Weights[j][i] == 0 ) y_Comb_Weights[j][i] = -10. ;
    }
  }
  
  graphCombPi0InvYieldStat = new TGraphAsymmErrors(nBinsPi0,pTPi0BinCenters,y_Comb,0,0,y_E_stat_Comb,y_E_stat_Comb);
  graphCombPi0InvYieldSys =  new TGraphAsymmErrors(nBinsPi0,pTPi0BinCenters,y_Comb,pTPi0BinWidths,pTPi0BinWidths,y_E_syst_Comb,y_E_syst_Comb);
  
  graphCombPi0InvYieldStat_Norm = new TGraphAsymmErrors(nBinsPi0,pTPi0BinCenters,y_Comb_Norm,0,0,y_E_stat_Comb_Norm,y_E_stat_Comb_Norm);
  graphCombPi0InvYieldSys_Norm =  new TGraphAsymmErrors(nBinsPi0,pTPi0BinCenters,y_Comb_Norm,pTPi0BinWidths,pTPi0BinWidths,y_E_syst_Comb_Norm,y_E_syst_Comb_Norm);
  
  for (Int_t j = 0; j < totalNSets; j++){
    graphCombPi0Weights[j] = new TGraphAsymmErrors(nBinsPi0,pTPi0BinCenters,y_Comb_Weights[j],0,0,0,0);
  }
  
  cout << "*******************************************" << endl;
  cout << "going to calculate the ratio of indiv measurements to the simple combination..." << endl;
  cout << "*******************************************" << endl;
  Double_t  y_RatioToComb[5][nBinsPi0];
  Double_t  y_RatioToCombTotError[5][nBinsPi0];
  for(Int_t i=0; i<nBinsPi0; i++){
    for (Int_t j = 0; j < totalNSets; j++){
      y_RatioToComb[j][i] = -100;
    }
  }
  for(Int_t i=0; i<nBinsPi0; i++){
    cout << endl;
    cout << "pT bin : " << pTPi0PbPb5TeV[i] << "-" << pTPi0PbPb5TeV[i+1] << "\t center :" << pTPi0BinCenters[i] << endl;
    for (Int_t j = 0; j < totalNSets; j++){
      if(directoryPi0[j] && graphPi0InvYieldStat[j] && j!=2 && (pTPi0BinCenters[i]>startpT[j] && pTPi0BinCenters[i]<endpT[j]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        y_RatioToComb[j][i] = y_Meas[j][i]/y_Comb[i];
        y_RatioToCombTotError[j][i] = sqrt(pow(y_E_stat_Meas[j][i],2)+pow(y_E_syst_Meas[j][i-offSetMethod[j]],2)) / y_Comb[i];
        cout << y_RatioToComb[j][i] << endl;
      }
      if(havePHOS && j==2 && (pTPi0BinCenters[i]>startpT[2] && pTPi0BinCenters[i]<endpT[2]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        y_RatioToComb[j][i] = y_Meas[j][i-1]/y_Comb[i];
        y_RatioToCombTotError[j][i] = sqrt(pow(y_E_stat_Meas[j][i-1],2)+pow(y_E_syst_Meas[j][i-offSetMethod[j]],2)) / y_Comb[i];
        cout << y_RatioToComb[j][i] << endl;
      }
    }
  }
  TGraphAsymmErrors* graphPi0InvYieldStatToComb[5]      = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0InvYieldSysToComb[5]       = {NULL,NULL,NULL,NULL,NULL};
  for (Int_t j = 0; j < totalNSets; j++){
    if(directoryPi0[j] || ( j==2 && havePHOS ) ){
      graphPi0InvYieldStatToComb[j] = new TGraphAsymmErrors(nBinsPi0,pTPi0BinCenters,y_RatioToComb[j],0,0,y_RatioToCombTotError[j],y_RatioToCombTotError[j]);
    }
  }
  
  cout << "*******************************************" << endl;
  cout << "going to calculate the Simple RAA..." << endl;
  cout << "*******************************************" << endl;
  
  Double_t recalcBarn             = 1e12;
  Double_t TAA = 1;
  if(centralityString.CompareTo("0-10%") == 0 )  TAA = 23.07;
  if(centralityString.CompareTo("10-20%") == 0 ) TAA = 14.27;
  if(centralityString.CompareTo("20-40%") == 0 ) TAA = 6.872;
  if(centralityString.CompareTo("40-60%") == 0 ) TAA = 2.046;
  if(centralityString.CompareTo("60-80%") == 0 ) TAA = 0.4173;
  
  Double_t y_RAA[nBinsPi0];
  Double_t y_RAA_error_stat[nBinsPi0];
  Double_t y_RAA_error_syst[nBinsPi0];
  for(Int_t i=0; i<nBinsPi0; i++){
      y_RAA[i] = -100;
      y_RAA_error_stat[i] = 0;
      y_RAA_error_syst[i] = 0;
  }
  for(Int_t i=0; i<nBinsPi0; i++){
    cout << endl;
    cout << "pT bin : " << pTPi0PbPb5TeV[i] << "-" << pTPi0PbPb5TeV[i+1] << "\t center :" << pTPi0BinCenters[i] << endl;
    y_RAA[i] = (y_Comb[i] * recalcBarn * 1e-3 ) / ( TAA * ppFit->Eval(pTPi0BinCenters[i]) );
    y_RAA_error_stat[i] = (y_E_stat_Comb[i] * recalcBarn * 1e-3 ) / ( TAA * ppFit->Eval(pTPi0BinCenters[i]) );
    y_RAA_error_syst[i] = 1.10 * (y_E_syst_Comb[i] * recalcBarn * 1e-3 ) / ( TAA * ppFit->Eval(pTPi0BinCenters[i]) );
    cout << y_RAA[i] << endl;
  }
  
  graphCombPi0RAAStat = new TGraphAsymmErrors(nBinsPi0,pTPi0BinCenters,y_RAA,0,0,y_RAA_error_stat,y_RAA_error_stat);
  graphCombPi0RAASys  = new TGraphAsymmErrors(nBinsPi0,pTPi0BinCenters,y_RAA,pTPi0BinWidths,pTPi0BinWidths,y_RAA_error_syst,y_RAA_error_syst);
  
  // **********************************************************************************************************************
  // ******************************************* Calculate simple combination      ****************************************
  // **********************************************************************************************************************
  cout << endl;
  cout << "*******************************************" << endl;
  cout << "*******************************************" << endl;
  cout << "*****************  Eta  *******************" << endl;
  cout << "*******************************************" << endl;
  cout << "*******************************************" << endl;
  cout << endl;

  TGraphAsymmErrors* graphCombEtaInvYieldStat  = NULL;
  TGraphAsymmErrors* graphCombEtaInvYieldSys   = NULL;
  TGraphAsymmErrors* graphCombEtaInvYieldStat_Norm  = NULL;
  TGraphAsymmErrors* graphCombEtaInvYieldSys_Norm   = NULL;
  TGraphAsymmErrors* graphCombEtaRAAStat  = NULL;
  TGraphAsymmErrors* graphCombEtaRAASys   = NULL;

  Double_t  pTEtaBinCenters[nBinsEta];
  Double_t  pTEtaBinWidths[nBinsEta];
  Double_t  y_CombEta[nBinsEta];
  Double_t  y_E_stat_CombEta[nBinsEta];
  Double_t  y_E_syst_CombEta[nBinsEta];

  Double_t  y_Comb_NormEta[nBinsEta];
  Double_t  y_E_stat_Comb_NormEta[nBinsEta];
  Double_t  y_E_syst_Comb_NormEta[nBinsEta];
  Double_t  y_Comb_WeightsEta[5][nBinsEta];

  for (Int_t j = 0; j < totalNSets; j++){
      if(directoryEta[j] && graphEtaInvYieldStat[j]){
        cout << nameMeasGlobal[j].Data() << "\t\t getting values" << endl;
        n_Meas[j] = graphEtaInvYieldStat[j]->GetN();
        x_Meas[j] = graphEtaInvYieldStat[j]->GetX();
        x_Meas_syst[j] = graphEtaInvYieldSys[j]->GetX();
        y_Meas[j] = graphEtaInvYieldStat[j]->GetY();
        y_E_stat_Meas[j] = graphEtaInvYieldStat[j]->GetEYlow();
        y_E_syst_Meas[j] = graphEtaInvYieldSys[j] ->GetEYlow();
      }
  }

  cout << "*******************************************" << endl;
  cout << "printing the arrays" << endl;
  cout << "*******************************************" << endl;
  for (Int_t j = 0; j < totalNSets; j++){
    if(directoryEta[j] && graphEtaInvYieldStat[j]){
      cout << nameMeasGlobal[j].Data() << "\t" << n_Meas[j] << endl;
      cout << "first x value: " << x_Meas[j][0] << endl;
      cout << "first x value with systematics: " << x_Meas_syst[j][0] << endl;
    }
  }

  cout << "*******************************************" << endl;
  cout << "starting to calculate a simple combination of measurements" << endl;
  cout << "*******************************************" << endl;
  for(Int_t i=0; i<nBinsEta; i++){
    pTEtaBinCenters[i] = ( pTEtaPbPb5TeV[i] + pTEtaPbPb5TeV[i+1] ) / 2. ;
    pTEtaBinWidths[i] = ( pTEtaPbPb5TeV[i+1] - pTEtaPbPb5TeV[i] ) / 2. ;
    cout << endl;
    cout << "pT bin : " << pTEtaPbPb5TeV[i] << "-" << pTEtaPbPb5TeV[i+1] << "\t center :" << pTEtaBinCenters[i] << endl;
    temp_Ny = 0;
    for (Int_t j = 0; j < totalNSets; j++){
      temp_y[j] = 0;
      temp_y_E_stat[j] = 0;
      temp_y_E_syst[j] = 0;
      if(directoryEta[j] && graphEtaInvYieldStat[j] && j!=2 && (pTEtaBinCenters[i]>startpTEta[j] && pTEtaBinCenters[i]<endpTEta[j]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        temp_y[j] = y_Meas[j][i];
        temp_y_E_stat[j] = y_E_stat_Meas[j][i];
        temp_y_E_syst[j] = y_E_syst_Meas[j][i-offSetMethodEta[j]];
        if(j==0) temp_y_E_syst[j] = temp_y[j] * 0.2;
        cout << temp_y[j] << endl;
        temp_Ny++;
      }
      if(havePHOS && j==2 && (pTEtaBinCenters[i]>startpTEta[2] && pTEtaBinCenters[i]<endpTEta[2]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        temp_y[j] = y_Meas[j][i-offSetStatPHOS];
        temp_y_E_stat[j] = y_E_stat_Meas[j][i-offSetStatPHOS];
        temp_y_E_syst[j] = y_E_syst_Meas[j][i-offSetMethodEta[j]];
        cout << temp_y[j] << endl;
        temp_Ny++;
      }
    }
    cout << "This pt bin has " << temp_Ny << " points to average " << endl;
    temp_y_comb      = 0;
    temp_y_comb_E_tot = 0;
    for (Int_t j = 0; j < totalNSets; j++){
      y_Comb_WeightsEta[j][i] = 0;
      if(temp_y_E_stat[j]>0 && temp_y_E_syst[j]>0){
        temp_y_comb += ( temp_y[j] / (pow(temp_y_E_stat[j],2) + pow(temp_y_E_syst[j],2) ) );
        temp_y_comb_E_tot += ( 1 / (pow(temp_y_E_stat[j],2) + pow(temp_y_E_syst[j],2) ) );
        y_Comb_WeightsEta[j][i] = ( 1 / (pow(temp_y_E_stat[j],2) + pow(temp_y_E_syst[j],2) ) );
      }
    }
    if(temp_y_comb_E_tot>0){
      y_CombEta[i]         = temp_y_comb / temp_y_comb_E_tot ;
      y_Comb_NormEta[i]    = 1 ;
      for (Int_t j = 0; j < totalNSets; j++){
        if(pTEtaBinCenters[i]>startpTEta[j] && pTEtaBinCenters[i]<endpTEta[j]){
          y_Comb_WeightsEta[j][i] /= temp_y_comb_E_tot;
        } else {
          y_Comb_WeightsEta[j][i] = 0.;
        }
      }
    } else {
      y_CombEta[i]         = 1e-40;
      y_Comb_NormEta[i]    = -10. ;
      for (Int_t j = 0; j < totalNSets; j++){
        y_Comb_WeightsEta[j][i] = 0. ;
      }
    }
    //go fully correlated -> if they are uncorrelated, take RMS
//     temp_Ny = 1.;
    cout << "y_CombEta = \t " << y_CombEta[i] << endl;
    EStat0 = 0; EStat1 = 0; EStat2 = 0; EStat3 = 0; EStat4 = 0;
    if(temp_y_E_stat[0]>0) EStat0 = pow(1./temp_y_E_stat[0],2);
    if(temp_y_E_stat[1]>0) EStat1 = pow(1./temp_y_E_stat[1],2);
    if(temp_y_E_stat[2]>0) EStat2 = pow(1./temp_y_E_stat[2],2);
    if(temp_y_E_stat[3]>0) EStat3 = pow(1./temp_y_E_stat[3],2);
    if(temp_y_E_stat[4]>0) EStat4 = pow(1./temp_y_E_stat[4],2);
    y_E_stat_CombEta[i]  = 1./sqrt(EStat0+EStat1+EStat2+EStat3+EStat4);
    y_E_syst_CombEta[i]  = sqrt(((pow(y_Comb_WeightsEta[0][i]*temp_y_E_syst[0],2)+pow(y_Comb_WeightsEta[1][i]*temp_y_E_syst[1],2)+pow(y_Comb_WeightsEta[2][i]*temp_y_E_syst[2],2)+pow(y_Comb_WeightsEta[3][i]*temp_y_E_syst[3],2)+pow(y_Comb_WeightsEta[4][i]*temp_y_E_syst[4],2))/temp_Ny)/(pow(y_Comb_WeightsEta[0][i],2)+pow(y_Comb_WeightsEta[1][i],2)+pow(y_Comb_WeightsEta[2][i],2)+pow(y_Comb_WeightsEta[3][i],2)+pow(y_Comb_WeightsEta[4][i],2)));
    if(y_CombEta[i] == 1e-40){
      y_E_stat_CombEta[i]  = 0;
      y_E_syst_CombEta[i]  = 0;
    }
    y_E_stat_Comb_NormEta[i]  = y_E_stat_CombEta[i] / y_CombEta[i];
    y_E_syst_Comb_NormEta[i]  = y_E_syst_CombEta[i] / y_CombEta[i];
    cout << y_E_stat_CombEta[i]  << "\t" << pow(temp_y_E_stat[0],2) << "\t" << pow(temp_y_E_stat[1],2) << "\t" << pow(temp_y_E_stat[2],2) << "\t" << pow(temp_y_E_stat[3],2) << "\t" << pow(temp_y_E_stat[4],2) << endl;
    cout << y_E_syst_CombEta[i]  << "\t" << pow(temp_y_E_syst[0],2) << "\t" << pow(temp_y_E_syst[1],2) << "\t" << pow(temp_y_E_syst[2],2) << "\t" << pow(temp_y_E_syst[3],2) << "\t" << pow(temp_y_E_syst[4],2) << endl;
    for (Int_t j = 0; j < totalNSets; j++){
      if( y_Comb_WeightsEta[j][i] == 0 ) y_Comb_WeightsEta[j][i] = -10. ;
    }
  }

  graphCombEtaInvYieldStat = new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_CombEta,0,0,y_E_stat_CombEta,y_E_stat_CombEta);
  graphCombEtaInvYieldSys =  new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_CombEta,pTEtaBinWidths,pTEtaBinWidths,y_E_syst_CombEta,y_E_syst_CombEta);

  graphCombEtaInvYieldStat_Norm = new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_Comb_NormEta,0,0,y_E_stat_Comb_NormEta,y_E_stat_Comb_NormEta);
  graphCombEtaInvYieldSys_Norm =  new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_Comb_NormEta,pTEtaBinWidths,pTEtaBinWidths,y_E_syst_Comb_NormEta,y_E_syst_Comb_NormEta);

  for (Int_t j = 0; j < totalNSets; j++){
    graphCombEtaWeights[j] = new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_Comb_WeightsEta[j],0,0,0,0);
  }

  cout << "*******************************************" << endl;
  cout << "going to calculate the ratio of indiv measurements to the simple combination..." << endl;
  cout << "*******************************************" << endl;
  Double_t  y_RatioToCombEta[5][nBinsEta];
  Double_t  y_RatioToCombTotErrorEta[5][nBinsEta];
  for(Int_t i=0; i<nBinsEta; i++){
    for (Int_t j = 0; j < totalNSets; j++){
      y_RatioToCombEta[j][i] = -100;
    }
  }
  for(Int_t i=0; i<nBinsEta; i++){
    cout << endl;
    cout << "pT bin : " << pTEtaPbPb5TeV[i] << "-" << pTEtaPbPb5TeV[i+1] << "\t center :" << pTEtaBinCenters[i] << endl;
    for (Int_t j = 0; j < totalNSets; j++){
      if(directoryEta[j] && graphEtaInvYieldStat[j] && j!=2 && (pTEtaBinCenters[i]>startpTEta[j] && pTEtaBinCenters[i]<endpTEta[j]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        y_RatioToCombEta[j][i] = y_Meas[j][i]/y_CombEta[i];
        y_RatioToCombTotErrorEta[j][i] = sqrt(pow(y_E_stat_Meas[j][i],2)+pow(y_E_syst_Meas[j][i-offSetMethodEta[j]],2)) / y_CombEta[i];
        cout << y_RatioToCombEta[j][i] << endl;
      }
      if(havePHOS && j==2 && (pTEtaBinCenters[i]>startpTEta[2] && pTEtaBinCenters[i]<endpTEta[2]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        y_RatioToCombEta[j][i] = y_Meas[j][i-offSetStatPHOS]/y_CombEta[i];
        y_RatioToCombTotErrorEta[j][i] = sqrt(pow(y_E_stat_Meas[j][i-offSetStatPHOS],2)+pow(y_E_syst_Meas[j][i-offSetMethodEta[j]],2)) / y_CombEta[i];
        cout << y_RatioToCombEta[j][i] << endl;
      }
    }
  }
  TGraphAsymmErrors* graphEtaInvYieldStatToComb[5]      = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaInvYieldSysToComb[5]       = {NULL,NULL,NULL,NULL,NULL};
  for (Int_t j = 0; j < totalNSets; j++){
    if(directoryEta[j] || ( j==2 && havePHOS ) ){
      graphEtaInvYieldStatToComb[j] = new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_RatioToCombEta[j],0,0,y_RatioToCombTotErrorEta[j],y_RatioToCombTotErrorEta[j]);
    }
  }

  cout << "*******************************************" << endl;
  cout << "going to calculate the Simple RAA..." << endl;
  cout << "*******************************************" << endl;

  Double_t y_RAAEta[nBinsEta];
  Double_t y_RAA_error_statEta[nBinsEta];
  Double_t y_RAA_error_systEta[nBinsEta];
  for(Int_t i=0; i<nBinsEta; i++){
      y_RAAEta[i] = -100;
      y_RAA_error_statEta[i] = 0;
      y_RAA_error_systEta[i] = 0;
  }
  for(Int_t i=0; i<nBinsEta; i++){
    cout << endl;
    cout << "pT bin : " << pTEtaPbPb5TeV[i] << "-" << pTEtaPbPb5TeV[i+1] << "\t center :" << pTEtaBinCenters[i] << endl;
    y_RAAEta[i] = (y_CombEta[i] * recalcBarn * 1e-3 ) / ( TAA * ppFitEta->Eval(pTEtaBinCenters[i]) );
    y_RAA_error_statEta[i] = (y_E_stat_CombEta[i] * recalcBarn * 1e-3 ) / ( TAA * ppFitEta->Eval(pTEtaBinCenters[i]) );
    y_RAA_error_systEta[i] = 1.10 * (y_E_syst_CombEta[i] * recalcBarn * 1e-3 ) / ( TAA * ppFitEta->Eval(pTEtaBinCenters[i]) );
    cout << y_RAAEta[i] << endl;
  }

  graphCombEtaRAAStat = new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_RAAEta,0,0,y_RAA_error_statEta,y_RAA_error_statEta);
  graphCombEtaRAASys  = new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_RAAEta,pTEtaBinWidths,pTEtaBinWidths,y_RAA_error_systEta,y_RAA_error_systEta);

  // **********************************************************************************************************************
  // ******************************************* Calculate simple combination      ****************************************
  // **********************************************************************************************************************
  cout << endl;
  cout << "*******************************************" << endl;
  cout << "*******************************************" << endl;
  cout << "*****************  EtaToPi0  *******************" << endl;
  cout << "*******************************************" << endl;
  cout << "*******************************************" << endl;
  cout << endl;

  TGraphAsymmErrors* graphCombEtaToPi0Stat  = NULL;
  TGraphAsymmErrors* graphCombEtaToPi0Sys   = NULL;
  TGraphAsymmErrors* graphCombEtaToPi0Stat_Norm  = NULL;
  TGraphAsymmErrors* graphCombEtaToPi0Sys_Norm   = NULL;

  for (Int_t j = 0; j < totalNSets; j++){
      if(directoryEta[j] && graphEtaToPi0Stat[j]){
        cout << nameMeasGlobal[j].Data() << "\t\t getting values" << endl;
        n_Meas[j] = graphEtaToPi0Stat[j]->GetN();
        x_Meas[j] = graphEtaToPi0Stat[j]->GetX();
        x_Meas_syst[j] = graphEtaToPi0Sys[j]->GetX();
        y_Meas[j] = graphEtaToPi0Stat[j]->GetY();
        y_E_stat_Meas[j] = graphEtaToPi0Stat[j]->GetEYlow();
        y_E_syst_Meas[j] = graphEtaToPi0Sys[j] ->GetEYlow();
      }
  }

  cout << "*******************************************" << endl;
  cout << "printing the arrays" << endl;
  cout << "*******************************************" << endl;
  for (Int_t j = 0; j < totalNSets; j++){
    if(directoryEta[j] && graphEtaToPi0Stat[j]){
      cout << nameMeasGlobal[j].Data() << "\t" << n_Meas[j] << endl;
      cout << "first x value: " << x_Meas[j][0] << endl;
      cout << "first x value with systematics: " << x_Meas_syst[j][0] << endl;
    }
  }

  cout << "*******************************************" << endl;
  cout << "starting to calculate a simple combination of measurements" << endl;
  cout << "*******************************************" << endl;
  for(Int_t i=0; i<nBinsEta; i++){
    pTEtaBinCenters[i] = ( pTEtaPbPb5TeV[i] + pTEtaPbPb5TeV[i+1] ) / 2. ;
    pTEtaBinWidths[i] = ( pTEtaPbPb5TeV[i+1] - pTEtaPbPb5TeV[i] ) / 2. ;
    cout << endl;
    cout << "pT bin : " << pTEtaPbPb5TeV[i] << "-" << pTEtaPbPb5TeV[i+1] << "\t center :" << pTEtaBinCenters[i] << endl;
    temp_Ny = 0;
    for (Int_t j = 0; j < totalNSets; j++){
      temp_y[j] = 0;
      temp_y_E_stat[j] = 0;
      temp_y_E_syst[j] = 0;
      if(directoryEta[j] && graphEtaToPi0Stat[j] && j!=2 && (pTEtaBinCenters[i]>startpTEta[j] && pTEtaBinCenters[i]<endpTEta[j]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        temp_y[j] = y_Meas[j][i];
        temp_y_E_stat[j] = y_E_stat_Meas[j][i];
        temp_y_E_syst[j] = y_E_syst_Meas[j][i-offSetMethodEta[j]];
        if(j==0) temp_y_E_syst[j] = temp_y[j] * 0.2;
        cout << temp_y[j] << endl;
        temp_Ny++;
      }
      if(havePHOS && j==2 && (pTEtaBinCenters[i]>startpTEta[2] && pTEtaBinCenters[i]<endpTEta[2]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        temp_y[j] = y_Meas[j][i-offSetStatPHOS];
        temp_y_E_stat[j] = y_E_stat_Meas[j][i-offSetStatPHOS];
        temp_y_E_syst[j] = y_E_syst_Meas[j][i-offSetMethodEta[j]];
        cout << temp_y[j] << endl;
        temp_Ny++;
      }
    }
    cout << "This pt bin has " << temp_Ny << " points to average " << endl;
    temp_y_comb      = 0;
    temp_y_comb_E_tot = 0;
    for (Int_t j = 0; j < totalNSets; j++){
      y_Comb_WeightsEta[j][i] = 0;
      if(temp_y_E_stat[j]>0 && temp_y_E_syst[j]>0){
        temp_y_comb += ( temp_y[j] / (pow(temp_y_E_stat[j],2) + pow(temp_y_E_syst[j],2) ) );
        temp_y_comb_E_tot += ( 1 / (pow(temp_y_E_stat[j],2) + pow(temp_y_E_syst[j],2) ) );
        y_Comb_WeightsEta[j][i] = ( 1 / (pow(temp_y_E_stat[j],2) + pow(temp_y_E_syst[j],2) ) );
      }
    }
    if(temp_y_comb_E_tot>0){
      y_CombEta[i]         = temp_y_comb / temp_y_comb_E_tot ;
      y_Comb_NormEta[i]    = 1 ;
      for (Int_t j = 0; j < totalNSets; j++){
        if(pTEtaBinCenters[i]>startpTEta[j] && pTEtaBinCenters[i]<endpTEta[j]){
          y_Comb_WeightsEta[j][i] /= temp_y_comb_E_tot;
        } else {
          y_Comb_WeightsEta[j][i] = 0.;
        }
      }
    } else {
      y_CombEta[i]         = 1e-40;
      y_Comb_NormEta[i]    = -10. ;
      for (Int_t j = 0; j < totalNSets; j++){
        y_Comb_WeightsEta[j][i] = 0. ;
      }
    }
    //go fully correlated -> if they are uncorrelated, take RMS
//     temp_Ny = 1.;
    cout << "y_CombEta = \t " << y_CombEta[i] << endl;
    EStat0 = 0; EStat1 = 0; EStat2 = 0; EStat3 = 0; EStat4 = 0;
    if(temp_y_E_stat[0]>0) EStat0 = pow(1./temp_y_E_stat[0],2);
    if(temp_y_E_stat[1]>0) EStat1 = pow(1./temp_y_E_stat[1],2);
    if(temp_y_E_stat[2]>0) EStat2 = pow(1./temp_y_E_stat[2],2);
    if(temp_y_E_stat[3]>0) EStat3 = pow(1./temp_y_E_stat[3],2);
    if(temp_y_E_stat[4]>0) EStat4 = pow(1./temp_y_E_stat[4],2);
    y_E_stat_CombEta[i]  = 1./sqrt(EStat0+EStat1+EStat2+EStat3+EStat4);
    y_E_syst_CombEta[i]  = sqrt(((pow(y_Comb_WeightsEta[0][i]*temp_y_E_syst[0],2)+pow(y_Comb_WeightsEta[1][i]*temp_y_E_syst[1],2)+pow(y_Comb_WeightsEta[2][i]*temp_y_E_syst[2],2)+pow(y_Comb_WeightsEta[3][i]*temp_y_E_syst[3],2)+pow(y_Comb_WeightsEta[4][i]*temp_y_E_syst[4],2))/temp_Ny)/(pow(y_Comb_WeightsEta[0][i],2)+pow(y_Comb_WeightsEta[1][i],2)+pow(y_Comb_WeightsEta[2][i],2)+pow(y_Comb_WeightsEta[3][i],2)+pow(y_Comb_WeightsEta[4][i],2)));
    if(y_CombEta[i] == 1e-40){
      y_E_stat_CombEta[i]  = 0;
      y_E_syst_CombEta[i]  = 0;
    }
    y_E_stat_Comb_NormEta[i]  = y_E_stat_CombEta[i] / y_CombEta[i];
    y_E_syst_Comb_NormEta[i]  = y_E_syst_CombEta[i] / y_CombEta[i];
    cout << y_E_stat_CombEta[i]  << "\t" << pow(temp_y_E_stat[0],2) << "\t" << pow(temp_y_E_stat[1],2) << "\t" << pow(temp_y_E_stat[2],2) << "\t" << pow(temp_y_E_stat[3],2) << "\t" << pow(temp_y_E_stat[4],2) << endl;
    cout << y_E_syst_CombEta[i]  << "\t" << pow(temp_y_E_syst[0],2) << "\t" << pow(temp_y_E_syst[1],2) << "\t" << pow(temp_y_E_syst[2],2) << "\t" << pow(temp_y_E_syst[3],2) << "\t" << pow(temp_y_E_syst[4],2) << endl;
    for (Int_t j = 0; j < totalNSets; j++){
      if( y_Comb_WeightsEta[j][i] == 0 ) y_Comb_WeightsEta[j][i] = -10. ;
    }
  }

  graphCombEtaToPi0Stat = new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_CombEta,0,0,y_E_stat_CombEta,y_E_stat_CombEta);
  graphCombEtaToPi0Sys =  new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_CombEta,pTEtaBinWidths,pTEtaBinWidths,y_E_syst_CombEta,y_E_syst_CombEta);

  graphCombEtaToPi0Stat_Norm = new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_Comb_NormEta,0,0,y_E_stat_Comb_NormEta,y_E_stat_Comb_NormEta);
  graphCombEtaToPi0Sys_Norm =  new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_Comb_NormEta,pTEtaBinWidths,pTEtaBinWidths,y_E_syst_Comb_NormEta,y_E_syst_Comb_NormEta);

  for (Int_t j = 0; j < totalNSets; j++){
    graphCombEtaToPi0Weights[j] = new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_Comb_WeightsEta[j],0,0,0,0);
  }

  cout << "*******************************************" << endl;
  cout << "going to calculate the ratio of indiv measurements to the simple combination..." << endl;
  cout << "*******************************************" << endl;
  for(Int_t i=0; i<nBinsEta; i++){
    for (Int_t j = 0; j < totalNSets; j++){
      y_RatioToCombEta[j][i] = -100;
    }
  }
  for(Int_t i=0; i<nBinsEta; i++){
    cout << endl;
    cout << "pT bin : " << pTEtaPbPb5TeV[i] << "-" << pTEtaPbPb5TeV[i+1] << "\t center :" << pTEtaBinCenters[i] << endl;
    for (Int_t j = 0; j < totalNSets; j++){
      if(directoryEta[j] && graphEtaToPi0Stat[j] && j!=2 && (pTEtaBinCenters[i]>startpTEta[j] && pTEtaBinCenters[i]<endpTEta[j]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        y_RatioToCombEta[j][i] = y_Meas[j][i]/y_CombEta[i];
        y_RatioToCombTotErrorEta[j][i] = sqrt(pow(y_E_stat_Meas[j][i],2)+pow(y_E_syst_Meas[j][i-offSetMethodEta[j]],2)) / y_CombEta[i];
        cout << y_RatioToCombEta[j][i] << endl;
      }
      if(havePHOS && j==2 && (pTEtaBinCenters[i]>startpTEta[2] && pTEtaBinCenters[i]<endpTEta[2]) ){
        cout << nameMeasGlobal[j].Data() << "\t\t";
        y_RatioToCombEta[j][i] = y_Meas[j][i-offSetStatPHOS]/y_CombEta[i];
        y_RatioToCombTotErrorEta[j][i] = sqrt(pow(y_E_stat_Meas[j][i-offSetStatPHOS],2)+pow(y_E_syst_Meas[j][i-offSetMethodEta[j]],2)) / y_CombEta[i];
        cout << y_RatioToCombEta[j][i] << endl;
      }
    }
  }
  TGraphAsymmErrors* graphEtaToPi0StatToComb[5]      = {NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaToPi0SysToComb[5]       = {NULL,NULL,NULL,NULL,NULL};
  for (Int_t j = 0; j < totalNSets; j++){
    if(directoryEta[j] || ( j==2 && havePHOS ) ){
      graphEtaToPi0StatToComb[j] = new TGraphAsymmErrors(nBinsEta,pTEtaBinCenters,y_RatioToCombEta[j],0,0,y_RatioToCombTotErrorEta[j],y_RatioToCombTotErrorEta[j]);
    }
  }

  // **********************************************************************************************************************
  // ******************************************* Mass and width for pi0            ****************************************
  // **********************************************************************************************************************

  Double_t arrayBoundariesX1_4[2]; Double_t arrayBoundariesY1_4[3]; Double_t relativeMarginsX[3]; Double_t relativeMarginsY[3];
  textSizeLabelsPixel             = 50;
  ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

  TCanvas* canvasMassWidthPi0     = new TCanvas("canvasMassWidthPi0","",0,0,1350,1250);  // gives the page size
  DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

  TPad* padWidthPi0               = new TPad("padWidthPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
  DrawGammaPadSettings( padWidthPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
  padWidthPi0->Draw();

  TPad* padMassPi0                = new TPad("padMassPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
  DrawGammaPadSettings( padMassPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
  padMassPi0->Draw();

  TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.32, 0.52, 0.52,-1, -1, -2);
  DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
  padMassLegend1->SetFillStyle(0);
  padMassLegend1->Draw();

  padWidthPi0->cd();
  padWidthPi0->SetLogx();

  Double_t margin                 = relativeMarginsX[0]*2.7*1350;
  Double_t textsizeLabelsWidth    = 0;
  Double_t textsizeFacWidth       = 0;
  if (padWidthPi0->XtoPixel(padWidthPi0->GetX2()) < padWidthPi0->YtoPixel(padWidthPi0->GetY1())){
      textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
      textsizeFacWidth            = (Double_t)1./padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
  } else {
      textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->YtoPixel(padWidthPi0->GetY1());
      textsizeFacWidth            = (Double_t)1./padWidthPi0->YtoPixel(padWidthPi0->GetY1());
  }

  TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, minPtPi0,maxPtPi0 ,1000., -30, 60);
  SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
  histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,49.5);//24.5);
  histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
  histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
  histo2DAllPi0FWHM->DrawCopy();

  for (Int_t i = 0; i < totalNSets; i++){
      if(histoPi0FWHMMeV[i] && histoPi0TrueFWHMMeV[i]){
          DrawGammaSetMarker(histoPi0FWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
          histoPi0FWHMMeV[i]->Draw("p,same,e");
          DrawGammaSetMarker(histoPi0TrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
          histoPi0TrueFWHMMeV[i]->Draw("p,same,e");
      }
  }

  TLatex *labelLegendAMass    = new TLatex(0.13,0.06,"a)");
  SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4);
  labelLegendAMass->SetTextFont(43);
  labelLegendAMass->Draw();

  TLatex *labelMassPerf       = new TLatex(0.13,0.87,"ALICE performance");
  SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4);
  labelMassPerf->SetTextFont(43);
  labelMassPerf->Draw();
  TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystemPbPb5TeV.Data());
  SetStyleTLatex( labelMassEnergy, textSizeLabelsPixel,4);
  labelMassEnergy->SetTextFont(43);
  labelMassEnergy->Draw();
  TLatex *labelMassPi0        = new TLatex(0.13,0.69,"#pi^{0} #rightarrow #gamma#gamma");
  SetStyleTLatex( labelMassPi0, textSizeLabelsPixel,4);
  labelMassPi0->SetTextFont(43);
  labelMassPi0->Draw();

  padMassPi0->cd();
  padMassPi0->SetLogx();

  Double_t textsizeLabelsMass         = 0;
  Double_t textsizeFacMass            = 0;
  if (padMassPi0->XtoPixel(padMassPi0->GetX2()) <padMassPi0->YtoPixel(padMassPi0->GetY1()) ){
      textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
      textsizeFacMass                 = (Double_t)1./padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
  } else {
      textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->YtoPixel(padMassPi0->GetY1());
      textsizeFacMass                 = (Double_t)1./padMassPi0->YtoPixel(padMassPi0->GetY1());
  }

  TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, minPtPi0,maxPtPi0, 1000., 100.1, 200.);//125.1, 155.9);
  SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
  histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
  histo2DAllPi0Mass->GetYaxis()->SetRangeUser(121.1, 179.9);//125.1, 155.9);
  histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0Mass->GetXaxis()->SetNoExponent();
  histo2DAllPi0Mass->DrawCopy();

  for (Int_t i = 0; i < totalNSets; i++){
      if(histoPi0Mass[i] && histoPi0TrueMass[i]){
          DrawGammaSetMarker(histoPi0Mass[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
          histoPi0Mass[i]->Draw("p,same,e");
          DrawGammaSetMarker(histoPi0TrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
          histoPi0TrueMass[i]->Draw("p,same,e");
      }
  }

  DrawGammaLines(minPtPi0,maxPtPi0 , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

  TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
  SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4);
  labelLegendBMass->SetTextFont(43);
  labelLegendBMass->Draw();

  //********************************** Defintion of the Legend **************************************************
  Double_t columnsLegendMass2[3]      = {0.,0.57,0.84};
  Double_t  rowsLegendMass2[14]= {0.84,0.66,0.50,0.33,0.01,0.16,0.16,0.16,0.16,0.16,0.16,0.16,0.16};
  //******************* Offsets ***********************
  Double_t offsetMarkerXMass2         = 0.1;
  Double_t offsetMarkerYMass2         = 0.1;
  //****************** Scale factors ******************
  Double_t scaleMarkerMass2           = 1.2;

  padMassLegend1->cd();
  //****************** first Column **************************************************
  TLatex *textMassPCM[10];
  for (Int_t i = 0; i < totalNSets; i++){
      if(histoPi0Mass[i] && histoPi0TrueMass[i] && histoPi0FWHMMeV[i] && histoPi0TrueFWHMMeV[i]){
          textMassPCM[i]                  = new TLatex(columnsLegendMass2[0],rowsLegendMass2[i+1],nameMeasGlobal[i].Data());
          SetStyleTLatex( textMassPCM[i], textSizeLabelsPixel,4);
          textMassPCM[i]->SetTextFont(43);
          textMassPCM[i]->Draw();
      }
  }
  //****************** second Column *************************************************
  TLatex *textMassData                = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"Data");
  SetStyleTLatex( textMassData, textSizeLabelsPixel,4);
  textMassData->SetTextFont(43);
  textMassData->Draw();
  TLatex *textMassMC                  = new TLatex(columnsLegendMass2[2] ,rowsLegendMass2[0],"MC");
  SetStyleTLatex( textMassMC, textSizeLabelsPixel,4);
  textMassMC->SetTextFont(43);
  textMassMC->Draw();

  TMarker* markerPCMPi0Mass[10];
  TMarker* markerPCMPi0MassMC[10];
  for (Int_t i = 0; i < totalNSets; i++){
      if(histoPi0Mass[i] && histoPi0TrueMass[i]){
          markerPCMPi0Mass[i]             = CreateMarkerFromHisto(histoPi0Mass[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
          markerPCMPi0Mass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
          markerPCMPi0MassMC[i]           = CreateMarkerFromHisto(histoPi0TrueMass[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
          markerPCMPi0MassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
      }
  }

  canvasMassWidthPi0->Update();
  canvasMassWidthPi0->SaveAs(Form("%s/Pi0_MassAndWidth.%s",outputDir.Data(),suffix.Data()));
                                                          
  // **********************************************************************************************************************
  // ******************************** Acceptance * Efficiency for pi0 single measurement **************************
  // **********************************************************************************************************************
  textSizeLabelsPixel             = 55;
  Double_t textSizeLabelsRel      = 55./1200;
  cout << textSizeLabelsRel << endl;

  TCanvas* canvasAcceptanceTimesEff       = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.1, 0.01, 0.015, 0.095);
  canvasAcceptanceTimesEff->SetLogy(1);
  canvasAcceptanceTimesEff->SetLogx(1);

  TH2F * histo2DAccEff;
  histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, minPtPi0, maxPtPi0, 1000, 1.01e-4, 10 );
  SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);//(#times #epsilon_{pur})
  histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
  histo2DAccEff->GetXaxis()->SetNoExponent();
  histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAccEff->DrawCopy();

  for (Int_t i = 0; i < totalNSets; i++){
      if(histoPi0AccTimesEff[i]){
          DrawGammaSetMarker(histoPi0AccTimesEff[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
          histoPi0AccTimesEff[i]->Draw("p,same,e");
      }
  }

  TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  for (Int_t i = 0; i < totalNSets; i++){
      if(histoPi0AccTimesEff[i]){
          legendEffiAccPi0->AddEntry(histoPi0AccTimesEff[i],nameMeasGlobal[i].Data(),"p");
      }
  }
  legendEffiAccPi0->Draw();

  drawLatexAdd("ALICE simulation",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystemPbPb5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasAcceptanceTimesEff->Update();
  canvasAcceptanceTimesEff->SaveAs(Form("%s/Pi0_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));
  
  // **********************************************************************************************************************
  // ******************************** Yields for pi0 single measurement                ************************************
  // **********************************************************************************************************************

  TCanvas* canvasYieldPi0  = new TCanvas("canvasYieldPi0","",200,10,1350,1350*1.15);  // gives the page size
  DrawGammaCanvasSettings( canvasYieldPi0, 0.14, 0.02, 0.02, 0.09);
  canvasYieldPi0->SetLogx();
  canvasYieldPi0->SetLogy();

  TH2F * histo2DYieldPi0;
  histo2DYieldPi0          = new TH2F("histo2DYieldPi0","histo2DYieldPi0",11000,minPtPi0,maxPtPi0,10000,minYieldPi0,maxYieldPi0);
  SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
  histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
  histo2DYieldPi0->GetXaxis()->SetNoExponent(kTRUE);
  histo2DYieldPi0->Draw("copy");

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryPi0[i]){
      DrawGammaSetMarkerTGraphAsym(graphPi0InvYieldStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphPi0InvYieldStat[i]->Draw("pEsame");
      DrawGammaSetMarkerTGraphAsym(graphPi0InvYieldSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
      graphPi0InvYieldSys[i]->Draw("E2same");
    }
  }

  TLatex *labelEnergyYieldPi0      = new TLatex(0.44,0.92,collisionSystemPbPb5TeV.Data());
  SetStyleTLatex( labelEnergyYieldPi0, 0.035,4);
  labelEnergyYieldPi0->Draw();
  TLatex *labelDetSysYieldPi0      = new TLatex(0.64,0.88,"#pi^{0} #rightarrow #gamma#gamma");
  SetStyleTLatex( labelDetSysYieldPi0, 0.035,4);
  labelDetSysYieldPi0->Draw();

  TLegend* legendYieldPi0          = new TLegend(0.62,0.62,0.9,0.86);
  legendYieldPi0->SetFillColor(0);
  legendYieldPi0->SetLineColor(0);
  legendYieldPi0->SetTextFont(42);
  legendYieldPi0->SetTextSize(0.035);
  for (Int_t i = 0; i < totalNSets; i++){
      if(graphPi0InvYieldSys[i]){
          legendYieldPi0->AddEntry(graphPi0InvYieldSys[i],nameMeasGlobal[i].Data(),"fp");
      }
  }
  
  legendYieldPi0->Draw();

  canvasYieldPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems.%s",outputDir.Data(),suffix.Data()));

  canvasYieldPi0->cd();
  histo2DYieldPi0->Draw("copy");

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryPi0[i]){
      graphPi0InvYieldStat[i]->Draw("pEsame");
      graphPi0InvYieldSys[i]->Draw("E2same");
    }
  }
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphCombPi0InvYieldSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombPi0InvYieldStat->Draw("p,same");

  labelEnergyYieldPi0->Draw();
  labelDetSysYieldPi0->Draw();

  legendYieldPi0->AddEntry(graphCombPi0InvYieldSys,"comb","fp");
  legendYieldPi0->Draw();

  canvasYieldPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));
  
  // **********************************************************************************************************************
  // ******************************** simple ratios to simple combination                        **************************
  // **********************************************************************************************************************

  TCanvas* canvasSimpleRatio       = new TCanvas("canvasSimpleRatio", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSimpleRatio,  0.1, 0.01, 0.015, 0.095);
  canvasSimpleRatio->SetLogx(1);

  TH2F * histo2DSimpleRatio;
  histo2DSimpleRatio                = new TH2F("histo2DSimpleRatio", "histo2DSimpleRatio",1000, minPtPi0, maxPtPi0, 1000, 0, 10 );
  SetStyleHistoTH2ForGraphs( histo2DSimpleRatio, "#it{p}_{T} (GeV/#it{c})", "Ratio",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);
  histo2DSimpleRatio->GetYaxis()->SetLabelOffset(0.001);
  histo2DSimpleRatio->GetYaxis()->SetRangeUser(0.01, 1.99);
  histo2DSimpleRatio->GetXaxis()->SetNoExponent();
  histo2DSimpleRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DSimpleRatio->DrawCopy();

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryPi0[i]){
      DrawGammaSetMarkerTGraphAsym(graphPi0InvYieldStatToComb[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphPi0InvYieldStatToComb[i]->Draw("pEsame");
    }
  }
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys_Norm, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphCombPi0InvYieldSys_Norm->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat_Norm, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombPi0InvYieldStat_Norm->Draw("p,same");

  TLegend* legendSimpleRatio         = GetAndSetLegend2(0.15, 0.13, 0.43, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  for (Int_t i = 0; i < totalNSets; i++){
      if(graphPi0InvYieldStatToComb[i]){
          legendSimpleRatio->AddEntry(graphPi0InvYieldStatToComb[i],nameMeasGlobal[i].Data(),"p");
      }
  }
  legendSimpleRatio->Draw();

  drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystemPbPb5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasSimpleRatio->Update();
  canvasSimpleRatio->SaveAs(Form("%s/Pi0_SimpleRatio.%s",outputDir.Data(),suffix.Data()));
  
  // **********************************************************************************************************************
  // ******************************** Plotting weights of combination                            **************************
  // **********************************************************************************************************************

  TCanvas* canvasSimpleWeights       = new TCanvas("canvasSimpleWeights", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSimpleWeights,  0.1, 0.01, 0.015, 0.095);
  canvasSimpleWeights->SetLogx(1);

  TH2F * histo2DSimpleWeights;
  histo2DSimpleWeights                = new TH2F("histo2DSimpleWeights", "histo2DSimpleWeights",1000, minPtPi0, maxPtPi0, 1000, -10, 10 );
  SetStyleHistoTH2ForGraphs( histo2DSimpleWeights, "#it{p}_{T} (GeV/#it{c})", "Weight",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);
  histo2DSimpleWeights->GetYaxis()->SetLabelOffset(0.001);
  histo2DSimpleWeights->GetYaxis()->SetRangeUser(-1., 1.99);
  histo2DSimpleWeights->GetXaxis()->SetNoExponent();
  histo2DSimpleWeights->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DSimpleWeights->DrawCopy();

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryPi0[i]){
      DrawGammaSetMarkerTGraphAsym(graphCombPi0Weights[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphCombPi0Weights[i]->Draw("pEsame");
    }
  }

  TLegend* legendSimpleWeights         = GetAndSetLegend2(0.15, 0.13, 0.43, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  for (Int_t i = 0; i < totalNSets; i++){
      if(graphCombPi0Weights[i]){
          legendSimpleWeights->AddEntry(graphCombPi0Weights[i],nameMeasGlobal[i].Data(),"p");
      }
  }
  legendSimpleWeights->Draw();

  drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystemPbPb5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasSimpleWeights->Update();
  canvasSimpleWeights->SaveAs(Form("%s/Pi0_SimpleWeights.%s",outputDir.Data(),suffix.Data()));
  
  // **********************************************************************************************************************
  // ******************************** simple RAA                                                 **************************
  // **********************************************************************************************************************

  TCanvas* canvasSimpleRAA       = new TCanvas("canvasSimpleRAA", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSimpleRAA,  0.1, 0.01, 0.015, 0.095);
  canvasSimpleRAA->SetLogx(1);

  TH2F * histo2DSimpleRAA;
  histo2DSimpleRAA                = new TH2F("histo2DSimpleRAA", "histo2DSimpleRAA",1000, minPtPi0, maxPtPi0, 1000, 0, 10 );
  SetStyleHistoTH2ForGraphs( histo2DSimpleRAA, "#it{p}_{T} (GeV/#it{c})", "R_{AA}",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);
  histo2DSimpleRAA->GetYaxis()->SetLabelOffset(0.001);
  histo2DSimpleRAA->GetYaxis()->SetRangeUser(0.01, 1.19);
  histo2DSimpleRAA->GetXaxis()->SetNoExponent();
  histo2DSimpleRAA->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DSimpleRAA->DrawCopy();

  DrawGammaSetMarkerTGraphAsym(graphCombPi0RAASys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphCombPi0RAASys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0RAAStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombPi0RAAStat->Draw("p,same,z");

  TLegend* legendSimpleRAA         = GetAndSetLegend2(0.15, 0.13, 0.43, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  legendSimpleRAA->AddEntry(graphCombPi0RAASys,"comb","fp");
//   legendSimpleRAA->Draw();

  drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystemPbPb5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasSimpleRAA->Update();
  canvasSimpleRAA->SaveAs(Form("%s/Pi0_SimpleRAA.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************** Yields for Eta single measurement                ************************************
  // **********************************************************************************************************************

  TCanvas* canvasYieldEta  = new TCanvas("canvasYieldEta","",200,10,1350,1350*1.15);  // gives the page size
  DrawGammaCanvasSettings( canvasYieldEta, 0.14, 0.02, 0.02, 0.09);
  canvasYieldEta->SetLogx();
  canvasYieldEta->SetLogy();

  TH2F * histo2DYieldEta;
  histo2DYieldEta          = new TH2F("histo2DYieldEta","histo2DYieldEta",11000,minPtEta,maxPtEta,10000,minYieldEta,maxYieldEta);
  SetStyleHistoTH2ForGraphs(histo2DYieldEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
  histo2DYieldEta->GetXaxis()->SetMoreLogLabels();
  histo2DYieldEta->GetXaxis()->SetNoExponent(kTRUE);
  histo2DYieldEta->Draw("copy");

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryEta[i]){
      DrawGammaSetMarkerTGraphAsym(graphEtaInvYieldStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphEtaInvYieldStat[i]->Draw("pEsame");
      DrawGammaSetMarkerTGraphAsym(graphEtaInvYieldSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
      graphEtaInvYieldSys[i]->Draw("E2same");
    }
  }

  TLatex *labelEnergyYieldEta      = new TLatex(0.44,0.92,collisionSystemPbPb5TeV.Data());
  SetStyleTLatex( labelEnergyYieldEta, 0.035,4);
  labelEnergyYieldEta->Draw();
  TLatex *labelDetSysYieldEta      = new TLatex(0.64,0.88,"#eta #rightarrow #gamma#gamma");
  SetStyleTLatex( labelDetSysYieldEta, 0.035,4);
  labelDetSysYieldEta->Draw();

  TLegend* legendYieldEta          = new TLegend(0.62,0.62,0.9,0.86);
  legendYieldEta->SetFillColor(0);
  legendYieldEta->SetLineColor(0);
  legendYieldEta->SetTextFont(42);
  legendYieldEta->SetTextSize(0.035);
  for (Int_t i = 0; i < totalNSets; i++){
      if(graphEtaInvYieldSys[i]){
          legendYieldEta->AddEntry(graphEtaInvYieldSys[i],nameMeasGlobal[i].Data(),"fp");
      }
  }

  legendYieldEta->Draw();

  canvasYieldEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems.%s",outputDir.Data(),suffix.Data()));

  canvasYieldEta->cd();
  histo2DYieldEta->Draw("copy");

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryEta[i]){
      graphEtaInvYieldStat[i]->Draw("pEsame");
      graphEtaInvYieldSys[i]->Draw("E2same");
    }
  }
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphCombEtaInvYieldSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombEtaInvYieldStat->Draw("p,same");

  labelEnergyYieldEta->Draw();
  labelDetSysYieldEta->Draw();

  legendYieldEta->AddEntry(graphCombEtaInvYieldSys,"comb","fp");
  legendYieldEta->Draw();

  canvasYieldEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************** simple ratios to simple combination                        **************************
  // **********************************************************************************************************************

  TCanvas* canvasSimpleRatioEta       = new TCanvas("canvasSimpleRatioEta", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSimpleRatioEta,  0.1, 0.01, 0.015, 0.095);
  canvasSimpleRatioEta->SetLogx(1);

  TH2F * histo2DSimpleRatioEta;
  histo2DSimpleRatioEta                = new TH2F("histo2DSimpleRatioEta", "histo2DSimpleRatioEta",1000, minPtEta, maxPtEta, 1000, -10, 10 );
  SetStyleHistoTH2ForGraphs( histo2DSimpleRatioEta, "#it{p}_{T} (GeV/#it{c})", "Ratio",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);
  histo2DSimpleRatioEta->GetYaxis()->SetLabelOffset(0.001);
  histo2DSimpleRatioEta->GetYaxis()->SetRangeUser(-0.99, 2.99);
  histo2DSimpleRatioEta->GetXaxis()->SetNoExponent();
  histo2DSimpleRatioEta->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DSimpleRatioEta->DrawCopy();

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryEta[i]){
      DrawGammaSetMarkerTGraphAsym(graphEtaInvYieldStatToComb[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphEtaInvYieldStatToComb[i]->Draw("pEsame");
    }
  }
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys_Norm, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphCombEtaInvYieldSys_Norm->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat_Norm, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombEtaInvYieldStat_Norm->Draw("p,same");

  TLegend* legendSimpleRatioEta         = GetAndSetLegend2(0.15, 0.13, 0.43, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  for (Int_t i = 0; i < totalNSets; i++){
      if(graphEtaInvYieldStatToComb[i]){
          legendSimpleRatioEta->AddEntry(graphEtaInvYieldStatToComb[i],nameMeasGlobal[i].Data(),"p");
      }
  }
  legendSimpleRatioEta->Draw();

  drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystemPbPb5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#eta #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasSimpleRatioEta->Update();
  canvasSimpleRatioEta->SaveAs(Form("%s/Eta_SimpleRatio.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************** Plotting weights of combination                            **************************
  // **********************************************************************************************************************

  TCanvas* canvasSimpleWeightsEta       = new TCanvas("canvasSimpleWeightsEta", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSimpleWeightsEta,  0.1, 0.01, 0.015, 0.095);
  canvasSimpleWeightsEta->SetLogx(1);

  TH2F * histo2DSimpleWeightsEta;
  histo2DSimpleWeightsEta                = new TH2F("histo2DSimpleWeightsEta", "histo2DSimpleWeightsEta",1000, minPtPi0, maxPtPi0, 1000, -10, 10 );
  SetStyleHistoTH2ForGraphs( histo2DSimpleWeightsEta, "#it{p}_{T} (GeV/#it{c})", "Weight",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);
  histo2DSimpleWeightsEta->GetYaxis()->SetLabelOffset(0.001);
  histo2DSimpleWeightsEta->GetYaxis()->SetRangeUser(-1., 1.99);
  histo2DSimpleWeightsEta->GetXaxis()->SetNoExponent();
  histo2DSimpleWeightsEta->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DSimpleWeightsEta->DrawCopy();

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryEta[i]){
      DrawGammaSetMarkerTGraphAsym(graphCombEtaWeights[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphCombEtaWeights[i]->Draw("pEsame");
    }
  }

  TLegend* legendSimpleWeightsEta         = GetAndSetLegend2(0.15, 0.13, 0.43, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  for (Int_t i = 0; i < totalNSets; i++){
      if(graphCombEtaWeights[i]){
          legendSimpleWeightsEta->AddEntry(graphCombEtaWeights[i],nameMeasGlobal[i].Data(),"p");
      }
  }
  legendSimpleWeightsEta->Draw();

  drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystemPbPb5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#eta #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasSimpleWeightsEta->Update();
  canvasSimpleWeightsEta->SaveAs(Form("%s/Eta_SimpleWeights.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************** simple RAA                                                 **************************
  // **********************************************************************************************************************

  TCanvas* canvasSimpleRAAEta       = new TCanvas("canvasSimpleRAAEta", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSimpleRAAEta,  0.1, 0.01, 0.015, 0.095);
  canvasSimpleRAAEta->SetLogx(1);

  TH2F * histo2DSimpleRAAEta;
  histo2DSimpleRAAEta                = new TH2F("histo2DSimpleRAAEta", "histo2DSimpleRAAEta",1000, minPtEta, maxPtEta, 1000, 0, 10 );
  SetStyleHistoTH2ForGraphs( histo2DSimpleRAAEta, "#it{p}_{T} (GeV/#it{c})", "R_{AA}",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);
  histo2DSimpleRAAEta->GetYaxis()->SetLabelOffset(0.001);
  histo2DSimpleRAAEta->GetYaxis()->SetRangeUser(0.01, 1.19);
  histo2DSimpleRAAEta->GetXaxis()->SetNoExponent();
  histo2DSimpleRAAEta->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DSimpleRAAEta->DrawCopy();

  DrawGammaSetMarkerTGraphAsym(graphCombEtaRAASys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphCombEtaRAASys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaRAAStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombEtaRAAStat->Draw("p,same,z");

  TLegend* legendSimpleRAAEta         = GetAndSetLegend2(0.15, 0.13, 0.43, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  legendSimpleRAAEta->AddEntry(graphCombEtaRAASys,"comb","fp");
//   legendSimpleRAAEta->Draw();

  drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystemPbPb5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#eta #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasSimpleRAAEta->Update();
  canvasSimpleRAAEta->SaveAs(Form("%s/Eta_SimpleRAA.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************** Yields for EtaToPi0 single measurement           ************************************
  // **********************************************************************************************************************

  TCanvas* canvasEtaToPi0  = new TCanvas("canvasEtaToPi0","",200,10,1350,1350);  // gives the page size
  DrawGammaCanvasSettings( canvasEtaToPi0, 0.14, 0.02, 0.02, 0.09);
  canvasEtaToPi0->SetLogx();
//   canvasEtaToPi0->SetLogy();

  TH2F * histo2DEtaToPi0;
  histo2DEtaToPi0          = new TH2F("histo2DEtaToPi0","histo2DEtaToPi0",11000,minPtEta,maxPtEta,1000,0.01,1.14);
  SetStyleHistoTH2ForGraphs(histo2DEtaToPi0, "#it{p}_{T} (GeV/#it{c})","#eta / #pi^{0}",0.035,0.04, 0.035,0.04, 0.9,1.45);
  histo2DEtaToPi0->GetXaxis()->SetMoreLogLabels();
  histo2DEtaToPi0->GetXaxis()->SetNoExponent(kTRUE);
  histo2DEtaToPi0->Draw("copy");

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryEta[i]){
      DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Stat[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphEtaToPi0Stat[i]->Draw("pEsame");
      DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Sys[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
      graphEtaToPi0Sys[i]->Draw("E2same");
    }
  }

  TLatex *labelEnergyEtaToPi0      = new TLatex(0.44,0.92,collisionSystemPbPb5TeV.Data());
  SetStyleTLatex( labelEnergyEtaToPi0, 0.035,4);
  labelEnergyEtaToPi0->Draw();
  TLatex *labelDetSysEtaToPi0      = new TLatex(0.64,0.88,"#eta / #pi^{0}");
  SetStyleTLatex( labelDetSysEtaToPi0, 0.035,4);
//   labelDetSysEtaToPi0->Draw();

  TLegend* legendEtaToPi0          = new TLegend(0.22,0.67,0.5,0.91);
  legendEtaToPi0->SetFillColor(0);
  legendEtaToPi0->SetLineColor(0);
  legendEtaToPi0->SetTextFont(42);
  legendEtaToPi0->SetTextSize(0.035);
  for (Int_t i = 0; i < totalNSets; i++){
      if(graphEtaToPi0Sys[i]){
          legendEtaToPi0->AddEntry(graphEtaToPi0Sys[i],nameMeasGlobal[i].Data(),"fp");
      }
  }

  legendEtaToPi0->Draw();

  canvasEtaToPi0->SaveAs(Form("%s/EtaToPi0_CompAllSystems.%s",outputDir.Data(),suffix.Data()));

  canvasEtaToPi0->cd();
  histo2DEtaToPi0->Draw("copy");

  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphCombEtaToPi0Sys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Stat, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombEtaToPi0Stat->Draw("p,same");

  labelEnergyEtaToPi0->Draw();
//   labelDetSysEtaToPi0->Draw();

  canvasEtaToPi0->SaveAs(Form("%s/EtaToPi0_OnlyComb.%s",outputDir.Data(),suffix.Data()));

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryEta[i]){
      graphEtaToPi0Stat[i]->Draw("pEsame");
      graphEtaToPi0Sys[i]->Draw("E2same");
    }
  }
  graphCombEtaToPi0Sys->Draw("E2same");
  graphCombEtaToPi0Stat->Draw("p,same");

  legendEtaToPi0->AddEntry(graphCombEtaToPi0Sys,"comb","fp");
  legendEtaToPi0->Draw();

  canvasEtaToPi0->SaveAs(Form("%s/EtaToPi0_CompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************** simple ratios to simple combination                        **************************
  // **********************************************************************************************************************

  TCanvas* canvasSimpleRatioEtaToPi0       = new TCanvas("canvasSimpleRatioEtaToPi0", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSimpleRatioEtaToPi0,  0.1, 0.01, 0.015, 0.095);
  canvasSimpleRatioEtaToPi0->SetLogx(1);

  TH2F * histo2DSimpleRatioEtaToPi0;
  histo2DSimpleRatioEtaToPi0                = new TH2F("histo2DSimpleRatioEtaToPi0", "histo2DSimpleRatioEtaToPi0",1000, minPtEta, maxPtEta, 1000, -10, 10 );
  SetStyleHistoTH2ForGraphs( histo2DSimpleRatioEtaToPi0, "#it{p}_{T} (GeV/#it{c})", "Ratio",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);
  histo2DSimpleRatioEtaToPi0->GetYaxis()->SetLabelOffset(0.001);
  histo2DSimpleRatioEtaToPi0->GetYaxis()->SetRangeUser(-0.99, 2.99);
  histo2DSimpleRatioEtaToPi0->GetXaxis()->SetNoExponent();
  histo2DSimpleRatioEtaToPi0->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DSimpleRatioEtaToPi0->DrawCopy();

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryEta[i]){
      DrawGammaSetMarkerTGraphAsym(graphEtaToPi0StatToComb[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphEtaToPi0StatToComb[i]->Draw("pEsame");
    }
  }
  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys_Norm, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphCombEtaToPi0Sys_Norm->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Stat_Norm, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombEtaToPi0Stat_Norm->Draw("p,same");

  TLegend* legendSimpleRatioEtaToPi0         = GetAndSetLegend2(0.15, 0.13, 0.43, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  for (Int_t i = 0; i < totalNSets; i++){
      if(graphEtaToPi0StatToComb[i]){
          legendSimpleRatioEtaToPi0->AddEntry(graphEtaToPi0StatToComb[i],nameMeasGlobal[i].Data(),"p");
      }
  }
  legendSimpleRatioEtaToPi0->Draw();

  drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystemPbPb5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#eta / #pi^{0}",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasSimpleRatioEtaToPi0->Update();
  canvasSimpleRatioEtaToPi0->SaveAs(Form("%s/EtaToPi0_SimpleRatio.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************** Plotting weights of combination                            **************************
  // **********************************************************************************************************************

  TCanvas* canvasSimpleWeightsEtaToPi0       = new TCanvas("canvasSimpleWeightsEtaToPi0", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSimpleWeightsEtaToPi0,  0.1, 0.01, 0.015, 0.095);
  canvasSimpleWeightsEtaToPi0->SetLogx(1);

  TH2F * histo2DSimpleWeightsEtaToPi0;
  histo2DSimpleWeightsEtaToPi0                = new TH2F("histo2DSimpleWeightsEtaToPi0", "histo2DSimpleWeightsEtaToPi0",1000, minPtPi0, maxPtPi0, 1000, -10, 10 );
  SetStyleHistoTH2ForGraphs( histo2DSimpleWeightsEtaToPi0, "#it{p}_{T} (GeV/#it{c})", "Weight",
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);
  histo2DSimpleWeightsEtaToPi0->GetYaxis()->SetLabelOffset(0.001);
  histo2DSimpleWeightsEtaToPi0->GetYaxis()->SetRangeUser(-1., 1.99);
  histo2DSimpleWeightsEtaToPi0->GetXaxis()->SetNoExponent();
  histo2DSimpleWeightsEtaToPi0->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DSimpleWeightsEtaToPi0->DrawCopy();

  for (Int_t i = 0; i < totalNSets; i++){
    if(directoryEta[i]){
      DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Weights[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphCombEtaToPi0Weights[i]->Draw("pEsame");
    }
  }

  TLegend* legendSimpleWeightsEtaToPi0         = GetAndSetLegend2(0.15, 0.13, 0.43, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  for (Int_t i = 0; i < totalNSets; i++){
      if(graphCombEtaToPi0Weights[i]){
          legendSimpleWeightsEtaToPi0->AddEntry(graphCombEtaToPi0Weights[i],nameMeasGlobal[i].Data(),"p");
      }
  }
  legendSimpleWeightsEtaToPi0->Draw();

  drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystemPbPb5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#eta / #pi^{0}",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasSimpleWeightsEtaToPi0->Update();
  canvasSimpleWeightsEtaToPi0->SaveAs(Form("%s/EtaToPi0_SimpleWeights.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************** writing the file                                           **************************
  // **********************************************************************************************************************

  TFile *fileOutput = new TFile(Form("%s/%s/PbPb_5.02TeV_Combination.root",suffix.Data(),dateForOutput.Data()), "UPDATE" );
  graphCombPi0InvYieldSys->Write(Form("%s_graphCombPi0InvYieldSys",CentNameShort.Data()));
  graphCombPi0InvYieldStat->Write(Form("%s_graphCombPi0InvYieldStat",CentNameShort.Data()));
  graphCombPi0RAASys->Write(Form("%s_graphCombPi0RAASys",CentNameShort.Data()));
  graphCombPi0RAAStat->Write(Form("%s_graphCombPi0RAAStat",CentNameShort.Data()));
  graphCombEtaInvYieldSys->Write(Form("%s_graphCombEtaInvYieldSys",CentNameShort.Data()));
  graphCombEtaInvYieldStat->Write(Form("%s_graphCombEtaInvYieldStat",CentNameShort.Data()));
  graphCombEtaRAASys->Write(Form("%s_graphCombEtaRAASys",CentNameShort.Data()));
  graphCombEtaRAAStat->Write(Form("%s_graphCombEtaRAAStat",CentNameShort.Data()));
  graphCombEtaToPi0Sys->Write(Form("%s_graphCombEtaToPi0Sys",CentNameShort.Data()));
  graphCombEtaToPi0Stat->Write(Form("%s_graphCombEtaToPi0Stat",CentNameShort.Data()));
  fileOutput->Close();

   
}
