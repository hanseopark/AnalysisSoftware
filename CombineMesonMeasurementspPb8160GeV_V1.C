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
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;

TGraphAsymmErrors* ScaleGraphAsym (TGraphAsymmErrors* graph, Double_t scaleFac){
	TGraphAsymmErrors* dummyGraph = (TGraphAsymmErrors*)graph->Clone(Form("%s_Scaled",graph->GetName()));
	Double_t * xValue         = dummyGraph->GetX();
	Double_t * yValue         = dummyGraph->GetY();
	Double_t* xErrorLow       = dummyGraph->GetEXlow();
	Double_t* xErrorHigh      = dummyGraph->GetEXhigh();
	Double_t* yErrorLow       = dummyGraph->GetEYlow();
	Double_t* yErrorHigh      = dummyGraph->GetEYhigh();
	Int_t nPoints             = dummyGraph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i]               = yValue[i]*scaleFac;
		yErrorLow[i]            = yErrorLow[i]*scaleFac;
		yErrorHigh[i]           = yErrorHigh[i]*scaleFac;
	}
	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
	return returnGraph;
}

TGraphAsymmErrors* ConvertYieldGraph(TGraphAsymmErrors* inputGraph, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt) {
  if (!inputGraph) {
      cout << "Error: Graph is NULL" << endl;
      return NULL;
  }
  if (DivideBy2pi) inputGraph                 = ScaleGraph(inputGraph, 1/(2*TMath::Pi()));
  if (MultiplyBy2pi) inputGraph               = ScaleGraph(inputGraph, 2*TMath::Pi());
  Double_t* xValue                            = inputGraph->GetX();
  Double_t* yValue                            = inputGraph->GetY();
  Double_t* xErrorLow                         = inputGraph->GetEXlow();
  Double_t* xErrorHigh                        = inputGraph->GetEXhigh();
  Double_t* yErrorLow                         = inputGraph->GetEYlow();
  Double_t* yErrorHigh                        = inputGraph->GetEYhigh();
  Int_t nPoints                               = inputGraph->GetN();
  if (DivideByPt || MultiplyByPt) {
    Double_t correctionValue                  = 1;
    for (Int_t i=0; i<nPoints; i++) {
      if (DivideByPt) correctionValue         = 1/xValue[i];
      if (MultiplyByPt) correctionValue       = xValue[i];
      yValue[i]                               = yValue[i]*correctionValue;
      yErrorLow[i]                            = yErrorLow[i]*correctionValue;
      yErrorHigh[i]                           = yErrorHigh[i]*correctionValue;
    }
  }
  inputGraph                                  = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
  return inputGraph;
}

void CombineMesonMeasurementspPb8160GeV_V1(
  TString fileNamePCM     = "",
  TString fileNamePHOS    = "",
  TString fileNameEMCal   = "",
  TString fileNamePCMPHOS = "",
  TString fileNamePCMEMCal= "",
  TString fileNamemEMC    = "",
  Int_t numbersofmeasPi0     = 3,
  Int_t numbersofmeasEta     = 2,
  Double_t minPtPi0Meas   = 0.4,
  Double_t maxPtPi0Meas   = 60,
  Double_t minPtEtaMeas   = 0.5,
  Double_t maxPtEtaMeas   = 60,
  TString bWCorrection    = "",
  TString suffix          = "pdf",
  TString fileNameReference  = "",
  Bool_t isNSD = kTRUE
){

  TString date                                    = ReturnDateString();
  Bool_t doOutput                                 = kTRUE;
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  StyleSettingsThesis();
  SetPlotStyle();

  TString dateForOutput                           = ReturnDateStringForOutput();
  cout << dateForOutput.Data() << endl;
  //___________________________________ Declaration of files _____________________________________________
  TString collisionSystempPb5TeV                     = "p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV";
  TString collisionSystempPb8TeV                     = "p-Pb, #sqrt{#it{s}_{NN}} = 8.16 TeV";
  TString collisionSystempp8TeV                     = "pp, #sqrt{#it{s}} = 8 TeV";
  TString labelALICEPlot                            = "ALICE, work-in-progress";
  TString outputDir                               = Form("%s/%s/CombineMesonMeasurementspPb8160GeV_V1",suffix.Data(),dateForOutput.Data());
  cout << outputDir.Data() << endl;

  TString fileNameTheory                          = "ExternalInput/Theory/TheoryCompilationPP.root";
  TString fileNameTheoryPPb                          = "ExternalInputpPb/Theory/TheoryCompilationPPb.root";
  TString fileNameCHadrons                        = "/home/nschmidt/AnalysisSoftware/CombinationInput5TeV/20180409_ChargedParticle_pp_5.02TeV.root";
  TString fileNameCParticles                      = "/home/nschmidt/AnalysisSoftware/CombinationInput5TeV/Spectra_ppLHC15n_Combined_Histograms.root";
  TString fileNameEtaToPi0                        = "ExternalInput/WorldDataPi0Eta.root";
  // TString fileNamePP8TeV                         = "/media/nschmidt/local/ANALYSIS/pp_8TeV_mEMC/pdf/8TeV/2019_01_29/FinalResultsTriggersPatchedLM1EMC-merged/data_EMC-mergedResultsFullCorrection_PP.root";
  // TString fileNamePP8TeV                         = "/media/nschmidt/local/ANALYSIS/pp_8TeV_mEMC_pTweights/pdf/8TeV/2019_07_19/FinalResultsTriggersPatchedLM1EMC-merged/data_EMC-mergedResultsFullCorrection_PP.root";
  TString fileNamePP8TeV                         = "/media/nschmidt/local/ANALYSIS/pp_8TeV_mEMC/pdf/8TeV/2019_06_05/FinalResultsTriggersPatchedLM1EMC-merged/data_EMC-mergedResultsFullCorrection_PP.root";

  gSystem->Exec("mkdir -p "+outputDir);
  if(fileNamePCM.CompareTo(""))
    gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
  if(fileNamePHOS.CompareTo(""))
    gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
  if(fileNameEMCal.CompareTo(""))
    gSystem->Exec(Form("cp %s %s/InputEMC.root", fileNameEMCal.Data(), outputDir.Data()));
  if(fileNamePCMPHOS.CompareTo(""))
    gSystem->Exec(Form("cp %s %s/InputPCMPHOS.root", fileNamePCMPHOS.Data(), outputDir.Data()));
  if(fileNamePCMEMCal.CompareTo(""))
    gSystem->Exec(Form("cp %s %s/InputPCMEMC.root", fileNamePCMEMCal.Data(), outputDir.Data()));
  if(fileNamemEMC.CompareTo(""))
    gSystem->Exec(Form("cp %s %s/InputmEMC.root", fileNamemEMC.Data(), outputDir.Data()));
  if(fileNameCHadrons.CompareTo(""))
    gSystem->Exec(Form("cp %s %s/ChargedHadrons.root", fileNameCHadrons.Data(), outputDir.Data()));
  if(fileNameCParticles.CompareTo(""))
    gSystem->Exec(Form("cp %s %s/ChardedParticles.root", fileNameCParticles.Data(), outputDir.Data()));
  if(fileNameReference.CompareTo(""))
    gSystem->Exec(Form("cp %s %s/ResultsPP8TeV.root", fileNameReference.Data(), outputDir.Data()));
  cout << "copied input files" << endl;

  fstream fLog;
  fLog.open(Form("%s/CombineMesonpPb8TeV.log",outputDir.Data()), ios::out);
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << dateForOutput.Data() << endl;

  Double_t nCollpPb                           = GetNCollFromName("00100", "pPb_8TeV");
  Double_t nCollErrpPb                        = GetNCollErrFromName("00100", "pPb_8TeV");
  Double_t tpPb                               = GetTAAFromName("", "pPb_8TeV");
  Double_t tpPbErr                            = GetTAAErrFromName("", "pPb_8TeV");

  Double_t mesonMassExpectPi0                     = TDatabasePDG::Instance()->GetParticle(111)->Mass();
  Double_t mesonMassExpectEta                     = TDatabasePDG::Instance()->GetParticle(221)->Mass();

  Width_t widthLinesBoxes                         = 1;
  Width_t widthCommonFit                          = 2;

  // Definition of colors, styles and markers sizes
  Color_t colorComb                               = kBlue+2;
  Style_t markerStyleComb                         = 20;
  Size_t  markerSizeComb                          = 2;

  Color_t  colorCGC                               = kCyan-8;
  Color_t  colorNLO                               = kAzure-4;

  Color_t colorTrigg[10]                          = {kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2};

  Color_t colorCombLowPt                          = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kFALSE);
  Color_t colorCombHighPt                         = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kTRUE);
  Style_t markerStyleCombLowPt                    = 20;
  Style_t markerStyleCombHighPt                   = 20;
  Size_t  markerSizeComparison                    = 2;
  Style_t styleMarkerNLOMuHalf                    = 24;
  Style_t styleMarkerNLOMuOne                     = 27;
  Style_t styleMarkerNLOMuTwo                     = 30;
  Style_t styleLineNLOMuHalf                      = 8;
  Style_t styleLineNLOMuOne                       = 7;
  Style_t styleLineNLOMuTwo                       = 4;
  Style_t styleLineNLOMuTwoBKK                    = 3;
  Style_t styleLineNLOMuTwoDSS                    = 6;
  Size_t  sizeMarkerNLO                           = 1;
  Width_t widthLineNLO                            = 2.;
  vector<TString>* ptSysRemNames                  = new vector<TString>[11];

  Int_t totalNSets                                = 13;
  TString nameMeasGlobal[13]                      = {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged", "PCMOtherDataset", "DCal", "PCM-DCal"};
  TString nameMeasGlobalPlot[13]                  = {"PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "mEMC", "PCMOtherDataset", "DMC", "PCM-DMC"};
  TString nameMeasReference[13]                   = {"PCM", "PHOS", "EMCAL", "PCMPHOS", "PCMEMCAL", "PCMDalitz", "PHOSDalitz", "EMCalDalitz", "spare", "EMCALMerged", "PCMOtherDataset", "DMC", "PCM-DMC"};
  Color_t colorDet[13];
  Color_t colorDetMC[13];
  Style_t markerStyleDet[13];
  Style_t markerStyleDetMC[13];
  Size_t  markerSizeDet[13];
  Size_t  markerSizeDetMC[13];

  for (Int_t i = 0; i < totalNSets; i++){
    colorDet[i]                                 = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
    colorDetMC[i]                               = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
    markerStyleDet[i]                           = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
    markerStyleDetMC[i]                         = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
    markerSizeDet[i]                            = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
    markerSizeDetMC[i]                          = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
  }
  cout << "loaded plotting styles" << endl;
  TFile* inputFile[13]                            =  {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    inputFile[0]                                = new TFile(fileNamePCM.Data());
    inputFile[1]                                = new TFile(fileNamePHOS.Data());
    inputFile[2]                                = new TFile(fileNameEMCal.Data());
    inputFile[3]                                = new TFile(fileNamePCMPHOS.Data());
    inputFile[4]                                = new TFile(fileNamePCMEMCal.Data());
    inputFile[9]                                = new TFile(fileNamemEMC.Data());

  TDirectory* directoryPi0[13]                    =  {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TDirectory* directoryEta[13]                    =  {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  for(Int_t i=0;i<13;i++){
    if(inputFile[i]){
      if(!inputFile[i]->IsZombie()){
        cout << "loading directories for " <<  nameMeasGlobalPlot[i] << endl;
        directoryPi0[i]                           = (TDirectory*)inputFile[i]->Get("Pi00-100%pPb_8TeV_V0A");
        directoryEta[i]                           = (TDirectory*)inputFile[i]->Get("Eta0-100%pPb_8TeV_V0A");
      }
    }
  }
  cout << __LINE__<<endl;
  TGraphAsymmErrors* graphPi0Mass[13]                          = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0FWHMMeV[13]                       = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0TrueMass[13]                      = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0TrueFWHMMeV[13]                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0Acc[13]                           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0TrueEffPt[13]                     = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0AccTimesEff[13]                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0InvYield[13]                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0InvYieldSys[13]   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0InvYieldStat[13]  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

  Bool_t havePi0NuclModFac[13]                   = {kFALSE};
  TH1D* histoPi0NuclModFac[13]                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0NuclModFacSys[13]   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphPi0NuclModFacStat[13]  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

  TGraphAsymmErrors* graphEtaMass[13]                          = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaFWHMMeV[13]                       = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaTrueMass[13]                      = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaTrueFWHMMeV[13]                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaAcc[13]                           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaTrueEffPt[13]                     = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaAccTimesEff[13]                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

  TH1D* histoEtaInvYield[13]                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaInvYieldSys[13]   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaInvYieldStat[13]  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  Bool_t haveEtaNuclModFac[13]                   = {kFALSE};
  TH1D* histoEtaNuclModFac[13]                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaNuclModFacSys[13]   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaNuclModFacStat[13]  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

  TH1D* histoEtaToPi0Stat[13]                     = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaToPi0Stat[13]        = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TGraphAsymmErrors* graphEtaToPi0Sys[13]         = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

  TH1D* histoPi0InvMassSigPlusBG[13]              = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0InvMassSig[13]                    = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0InvMassSigRemBGSub[13]            = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0InvMassBG[13]                     = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0InvMassRemBG[13]                  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TH1D* histoPi0InvMassBGTot[13]                  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TF1* fitPi0InvMassSig[13]                       = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  TF1* fitPi0InvMassBG[13]                        = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

  Bool_t haveAllPi0InvMass[13]                    = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
  TString strTriggerExampleBin[13]                = {"INT7", "INT7", "INT7", "INT7","INT7","INT7", "INT7", "INT7", "INT7","INT7","INT7", "INT7", "INT7"};
  Double_t rapidityMeas[13]                       = {1.6, 1., 1.6, 1.6, 1.6, 1., 1., 1., 1., 1., 1., 1., 1.};

  Double_t minPtPi0                               = 0.25;
  Double_t maxPtPi0                               = 250.0; //75
  Double_t maxPtPi0woMerged                       = 46.0; //75
  Double_t maxPtPi0Theory                         = 200.0; //75
  Double_t delaPtPi0                              = maxPtPi0 - minPtPi0;
  Double_t prodPtPi0                              = maxPtPi0 * minPtPi0;
  Double_t minInvYieldPi0                         = 6e-15;
  Double_t maxInvYieldPi0                         = 3e1;

  Double_t minPtEta                               = 0.3;
  Double_t maxPtEta                               = 60.0;//60
  Double_t maxPtEtaTheory                         = 45.0;//60
  Double_t delaPtEta                              = maxPtEta - minPtEta;
  Double_t prodPtEta                              = maxPtEta * minPtEta;
  Double_t minInvYieldEta                         = 6e-8;
  Double_t maxInvYieldEta                         = 9e-1;

  Double_t scalingToNSD                       = 0.97;
  if (!isNSD)
      scalingToNSD                            = 1;

  // *******************************************************************************************************
  // ************************** read input files ***********************************************************
  // *******************************************************************************************************
  Int_t nNuclModFacAvailPi0 = 0;
  Int_t nNuclModFacAvailEta = 0;
  for (Int_t i = 0; i < totalNSets; i++){
    cout << "reading from " << nameMeasGlobalPlot[i].Data() << " measurement..." << endl;
    if(directoryPi0[i]){
      cout << "loading pi0 inputs" << endl;
      // load mass/width/effi plots
      graphPi0Mass[i]                           = (TGraphAsymmErrors*)directoryPi0[i]->Get("Pi0_Mass_data");
      graphPi0FWHMMeV[i]                        = (TGraphAsymmErrors*)directoryPi0[i]->Get("Pi0_Width_data");
      graphPi0TrueMass[i]                       = (TGraphAsymmErrors*)directoryPi0[i]->Get("Pi0_Mass_MC");
      graphPi0TrueFWHMMeV[i]                    = (TGraphAsymmErrors*)directoryPi0[i]->Get("Pi0_Width_MC");
      if(graphPi0Mass[i]) graphPi0Mass[i]= ScaleGraphAsym(graphPi0Mass[i],1000.);
      if(graphPi0FWHMMeV[i]) graphPi0FWHMMeV[i]= ScaleGraphAsym(graphPi0FWHMMeV[i],1000.);
      if(graphPi0TrueMass[i]) graphPi0TrueMass[i]= ScaleGraphAsym(graphPi0TrueMass[i],1000.);
      if(graphPi0TrueFWHMMeV[i]) graphPi0TrueFWHMMeV[i]= ScaleGraphAsym(graphPi0TrueFWHMMeV[i],1000.);

      graphPi0Acc[i]                              = (TGraphAsymmErrors*)directoryPi0[i]->Get("AcceptancePi0");
      graphPi0TrueEffPt[i]                        = (TGraphAsymmErrors*)directoryPi0[i]->Get("EfficiencyPi0");

      graphPi0AccTimesEff[i]                     = (TGraphAsymmErrors*)directoryPi0[i]->Get("EffTimesAccPi0");
      // graphPi0AccTimesEff[i]= ScaleGraphAsym(graphPi0AccTimesEff[i],2*TMath::Pi()*rapidityMeas[i]);

      // load cross section systematics and datapoints
      histoPi0InvYield[i]                      = (TH1D*)directoryPi0[i]->Get("CorrectedYieldPi0");
      graphPi0InvYieldStat[i]                  = new TGraphAsymmErrors(histoPi0InvYield[i]);
      graphPi0InvYieldSys[i]                   = (TGraphAsymmErrors*)directoryPi0[i]->Get("Pi0SystError");

      if (isNSD) {
        histoPi0InvYield[i]->Scale(scalingToNSD);
        graphPi0InvYieldStat[i]                 = ScaleGraph(graphPi0InvYieldStat[i],scalingToNSD);
        graphPi0InvYieldSys[i]                  = ScaleGraph(graphPi0InvYieldSys[i],scalingToNSD);
      }
      cout << nameMeasGlobalPlot[i].Data() << " pi0 stat:" << graphPi0InvYieldStat[i] << endl;
      if(graphPi0InvYieldStat[i]) graphPi0InvYieldStat[i]->Print();
      else {cout << "stat graph missing. returning..." << endl; return;}
      cout << nameMeasGlobalPlot[i].Data() << " pi0 sys:" << graphPi0InvYieldSys[i] << endl;
      if(graphPi0InvYieldSys[i]) graphPi0InvYieldSys[i]->Print();
      else {cout << "sys graph missing. returning..." << endl; return;}

      // load cross section systematics and datapoints
      histoPi0NuclModFac[i]                      = (TH1D*)directoryPi0[i]->Get("histoNuclearModFactorPi0Stat");
      graphPi0NuclModFacStat[i]                  = (TGraphAsymmErrors*)directoryPi0[i]->Get("graphNuclearModFactorPi0Stat");
      graphPi0NuclModFacSys[i]                   = (TGraphAsymmErrors*)directoryPi0[i]->Get("graphNuclearModFactorPi0Sys");
      if(histoPi0NuclModFac[i] && graphPi0NuclModFacStat[i] && graphPi0NuclModFacSys[i]){
        nNuclModFacAvailPi0+=1;
        havePi0NuclModFac[i] = kTRUE;
        cout << nameMeasGlobalPlot[i].Data() << " Pi0 RpA stat:" << graphPi0NuclModFacStat[i] << endl;
        if(doOutput) graphPi0NuclModFacStat[i]->Print();
        cout << nameMeasGlobalPlot[i].Data() << " Pi0 RpA sys:" << graphPi0NuclModFacSys[i] << endl;
        if(doOutput) graphPi0NuclModFacSys[i]->Print();
      }

      // load invariant mass example bins
      histoPi0InvMassSig[i]                       = (TH1D*)directoryPi0[i]->Get(Form("Pi0_InvMassSig_Example_%s",strTriggerExampleBin[i].Data()));
      histoPi0InvMassSigPlusBG[i]                 = (TH1D*)directoryPi0[i]->Get(Form("Pi0_InvMassSigPlusBG_Example_%s",strTriggerExampleBin[i].Data()));
      histoPi0InvMassBG[i]                        = (TH1D*)directoryPi0[i]->Get(Form("Pi0_InvMassBG_Example_%s",strTriggerExampleBin[i].Data()));
      fitPi0InvMassSig[i]                         = (TF1*)directoryPi0[i]->Get(Form("Pi0_InvMassSigFit_Example_%s",strTriggerExampleBin[i].Data()));

      if (histoPi0InvMassSig[i] && histoPi0InvMassSigPlusBG[i] && histoPi0InvMassBG[i] && fitPi0InvMassSig[i]){
        haveAllPi0InvMass[i]                      = kTRUE;
      }
    }
    if(directoryEta[i]){
      cout << "loading eta inputs" << endl;
      // load mass/width/effi plots

      graphEtaMass[i]                             = (TGraphAsymmErrors*)directoryEta[i]->Get("Eta_Mass_data");
      graphEtaFWHMMeV[i]                          = (TGraphAsymmErrors*)directoryEta[i]->Get("Eta_Width_data");
      graphEtaTrueMass[i]                         = (TGraphAsymmErrors*)directoryEta[i]->Get("Eta_Mass_MC");
      graphEtaTrueFWHMMeV[i]                      = (TGraphAsymmErrors*)directoryEta[i]->Get("Eta_Width_MC");

      if(graphEtaMass[i]) graphEtaMass[i]= ScaleGraphAsym(graphEtaMass[i],1000.);
      if(graphEtaFWHMMeV[i]) graphEtaFWHMMeV[i]= ScaleGraphAsym(graphEtaFWHMMeV[i],1000.);
      if(graphEtaTrueMass[i]) graphEtaTrueMass[i]= ScaleGraphAsym(graphEtaTrueMass[i],1000.);
      if(graphEtaTrueFWHMMeV[i]) graphEtaTrueFWHMMeV[i]= ScaleGraphAsym(graphEtaTrueFWHMMeV[i],1000.);

      graphEtaAcc[i]                              = (TGraphAsymmErrors*)directoryEta[i]->Get("AcceptanceEta");
      graphEtaTrueEffPt[i]                        = (TGraphAsymmErrors*)directoryEta[i]->Get("EfficiencyEta");

      graphEtaAccTimesEff[i]                     = (TGraphAsymmErrors*)directoryEta[i]->Get("EffTimesAccEta");
      // graphEtaAccTimesEff[i]= ScaleGraphAsym(graphEtaAccTimesEff[i],2*TMath::Pi()*rapidityMeas[i]);

      // load cross section systematics and datapoints
      histoEtaInvYield[i]                      = (TH1D*)directoryEta[i]->Get("CorrectedYieldEta");
      graphEtaInvYieldStat[i]                  = new TGraphAsymmErrors(histoEtaInvYield[i]);
      graphEtaInvYieldSys[i]                   = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaSystError");

      if (isNSD) {
        histoEtaInvYield[i]->Scale(scalingToNSD);
        graphEtaInvYieldStat[i]                 = ScaleGraph(graphEtaInvYieldStat[i],scalingToNSD);
        graphEtaInvYieldSys[i]                  = ScaleGraph(graphEtaInvYieldSys[i],scalingToNSD);
      }

      cout << nameMeasGlobalPlot[i].Data() << " eta stat:" << graphEtaInvYieldStat[i] << endl;
      if(doOutput) graphEtaInvYieldStat[i]->Print();
      cout << nameMeasGlobalPlot[i].Data() << " eta sys:" << graphEtaInvYieldSys[i] << endl;
      if(doOutput) graphEtaInvYieldSys[i]->Print();

      // load cross section systematics and datapoints
      histoEtaNuclModFac[i]                      = (TH1D*)directoryEta[i]->Get("histoNuclearModFactorEtaStat");
      graphEtaNuclModFacStat[i]                  = (TGraphAsymmErrors*)directoryEta[i]->Get("graphNuclearModFactorEtaStat");
      graphEtaNuclModFacSys[i]                   = (TGraphAsymmErrors*)directoryEta[i]->Get("graphNuclearModFactorEtaSys");
      if(histoEtaNuclModFac[i] && graphEtaNuclModFacStat[i] && graphEtaNuclModFacSys[i]){
        nNuclModFacAvailEta+=1;
        haveEtaNuclModFac[i] = kTRUE;
        cout << nameMeasGlobalPlot[i].Data() << " eta RpA stat:" << graphEtaNuclModFacStat[i] << endl;
        if(doOutput) graphEtaNuclModFacStat[i]->Print();
        cout << nameMeasGlobalPlot[i].Data() << " eta RpA sys:" << graphEtaNuclModFacSys[i] << endl;
        if(doOutput) graphEtaNuclModFacSys[i]->Print();
      }
      histoEtaToPi0Stat[i]                        = (TH1D*)directoryEta[i]->Get("EtaToPi0StatError");
      graphEtaToPi0Stat[i]                        = (TGraphAsymmErrors*)directoryEta[i]->Get("graphEtaToPi0StatError");
      graphEtaToPi0Sys[i]                         = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0SystError");

      cout << nameMeasGlobalPlot[i].Data() << " eta/pi0 stat:" << histoEtaToPi0Stat[i] << endl;
      cout << nameMeasGlobalPlot[i].Data() << " eta/pi0 stat:" << graphEtaToPi0Stat[i] << endl;
      if(doOutput) if(graphEtaToPi0Stat[i])graphEtaToPi0Stat[i]->Print();
      cout << nameMeasGlobalPlot[i].Data() << " eta/pi0 sys:" << graphEtaToPi0Sys[i] << endl;
      if(doOutput) if(graphEtaToPi0Sys[i]) graphEtaToPi0Sys[i]->Print();
    }


    for (Int_t i = 0; i < totalNSets; i++){
      if (haveAllPi0InvMass[i]){
        histoPi0InvMassBGTot[i]                   = (TH1D*)histoPi0InvMassBG[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameMeasGlobalPlot[i].Data()));
        histoPi0InvMassSigRemBGSub[i]             = (TH1D*)histoPi0InvMassSig[i]->Clone(Form("Pi0_InvMassSigRemBGSub_Example_%s",nameMeasGlobalPlot[i].Data()));
      }
    }
  }

  cout << __LINE__ << endl;
  // *******************************************************************************************************
  // ************************** Loading other energies      ************************************************
  // *******************************************************************************************************
  TFile* filePP8TeV                               = new TFile(fileNamePP8TeV.Data());
  TH1D* histoPP8TeVPi0Stat                        = (TH1D*) filePP8TeV->Get("Pi08TeV/CorrectedYieldPi0");
  TGraphAsymmErrors* graphPP8TeVPi0Stat           = (TGraphAsymmErrors*) filePP8TeV->Get("Pi08TeV/graphCorrectedYieldPi0");
  TGraphAsymmErrors* graphPP8TeVPi0Sys            = (TGraphAsymmErrors*) filePP8TeV->Get("Pi08TeV/Pi0SystError");

  TFile* fileTheoryCompilation                    = new TFile(fileNameTheory.Data());
  TFile* fileTheoryCompilationPPb                    = new TFile(fileNameTheoryPPb.Data());

  // *******************************************************************************************************
  // ************************** Loading reference spectra ***********************************************
  // *******************************************************************************************************
  TGraphAsymmErrors* statErrorCollectionPi0PP[11];
  TGraphAsymmErrors* systErrorCollectionPi0PP[11];
  TGraphAsymmErrors* systErrorUnCorrCollectionPi0PP[11];
  TGraphAsymmErrors* systErrorInterCollectionPi0PP[11];
  TGraphAsymmErrors* statErrorCollectionPi0PPFitRatio_woXErr[11];
  TGraphAsymmErrors* systErrorCollectionPi0PPFitRatio[11];

  TGraphAsymmErrors* statErrorCollectionEtaPP[11];
  TGraphAsymmErrors* systErrorCollectionEtaPP[11];
  TGraphAsymmErrors* systErrorUnCorrCollectionEtaPP[11];
  TGraphAsymmErrors* systErrorInterCollectionEtaPP[11];
  Bool_t haveRefPPPi0[11]                                 = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                              kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                              kFALSE };
  Bool_t haveRefPPEta[11]                                 = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                              kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                              kFALSE };
  for (Int_t i = 0; i< 11; i++){
    statErrorCollectionPi0PP[i]         = NULL;
    systErrorCollectionPi0PP[i]         = NULL;
    systErrorUnCorrCollectionPi0PP[i]   = NULL;
    systErrorInterCollectionPi0PP[i]    = NULL;
    statErrorCollectionEtaPP[i]         = NULL;
    systErrorCollectionEtaPP[i]         = NULL;
    systErrorUnCorrCollectionEtaPP[i]   = NULL;
    systErrorInterCollectionEtaPP[i]    = NULL;
  }

  TFile* fileReference                        = new TFile(fileNameReference.Data());
  TF1* fitPPInvXsectionPi0                        = (TF1*)fileReference->Get("Pi08TeV/TwoComponentModelFitPi0");
  for (Int_t i = 0; i< 11; i++){
    statErrorCollectionPi0PP[i]                 = (TGraphAsymmErrors*)fileReference->Get(Form("Pi08TeV/graphInvCrossSectionPi0%s8TeVStatErr_yShifted",nameMeasReference[i].Data()));
    if(statErrorCollectionPi0PP[i]){
      statErrorCollectionPi0PPFitRatio_woXErr[i] = CalculateGraphErrRatioToFit(statErrorCollectionPi0PP[i],fitPPInvXsectionPi0);
      ProduceGraphAsymmWithoutXErrors(statErrorCollectionPi0PPFitRatio_woXErr[i]);
    }
    systErrorCollectionPi0PP[i]                 = (TGraphAsymmErrors*)fileReference->Get(Form("Pi08TeV/graphInvCrossSectionPi0%s8TeVSysErr_yShifted",nameMeasReference[i].Data()));
    if(systErrorCollectionPi0PP[i])
      systErrorCollectionPi0PPFitRatio[i] = CalculateGraphErrRatioToFit(systErrorCollectionPi0PP[i],fitPPInvXsectionPi0);
    systErrorUnCorrCollectionPi0PP[i]           = (TGraphAsymmErrors*)fileReference->Get(Form("Pi08TeV/graphInvCrossSectionPi0%s8TeVSysErr_yShifted",nameMeasReference[i].Data())); //sUnCorrSysErr
    systErrorInterCollectionPi0PP[i]            = (TGraphAsymmErrors*)fileReference->Get(Form("Pi08TeV/graphInvCrossSectionPi0%s8TeVSysErr_yShifted",nameMeasReference[i].Data())); //sInterpolSysErr
    if (statErrorCollectionPi0PP[i] && systErrorCollectionPi0PP[i] && systErrorUnCorrCollectionPi0PP[i] && systErrorInterCollectionPi0PP[i]){
      haveRefPPPi0[i]                         = kTRUE;
      // statErrorCollectionPi0PP[i] = ScaleGraphAsym(statErrorCollectionPi0PP[i],xSection8TeVINEL/xSection8TeVV0AND);
      // systErrorCollectionPi0PP[i] = ScaleGraphAsym(systErrorCollectionPi0PP[i],xSection8TeVINEL/xSection8TeVV0AND);
      // systErrorUnCorrCollectionPi0PP[i] = ScaleGraphAsym(systErrorUnCorrCollectionPi0PP[i],xSection8TeVINEL/xSection8TeVV0AND);
      // systErrorInterCollectionPi0PP[i] = ScaleGraphAsym(systErrorInterCollectionPi0PP[i],xSection8TeVINEL/xSection8TeVV0AND);
    }

    statErrorCollectionEtaPP[i]                 = (TGraphAsymmErrors*)fileReference->Get(Form("Eta8TeV/graphInvCrossSectionEta%s8TeVStatErr_yShifted",nameMeasReference[i].Data()));
    systErrorCollectionEtaPP[i]                 = (TGraphAsymmErrors*)fileReference->Get(Form("Eta8TeV/graphInvCrossSectionEta%s8TeVSysErr_yShifted",nameMeasReference[i].Data()));
    systErrorUnCorrCollectionEtaPP[i]           = (TGraphAsymmErrors*)fileReference->Get(Form("Eta8TeV/graphInvCrossSectionEta%s8TeVSysErr_yShifted",nameMeasReference[i].Data())); //sUnCorrSysErr
    systErrorInterCollectionEtaPP[i]            = (TGraphAsymmErrors*)fileReference->Get(Form("Eta8TeV/graphInvCrossSectionEta%s8TeVSysErr_yShifted",nameMeasReference[i].Data())); //sInterpolSysErr
    if (statErrorCollectionEtaPP[i] && systErrorCollectionEtaPP[i] && systErrorUnCorrCollectionEtaPP[i] && systErrorInterCollectionEtaPP[i])
      haveRefPPEta[i]                         = kTRUE;
    if (haveRefPPPi0[i])
      cout << "found pi0 pp reference for " << nameMeasReference[i].Data() << endl;
    if (haveRefPPEta[i])
      cout << "found eta pp reference for " << nameMeasReference[i].Data() << endl;
  }
  Double_t xSection8TeV                       = ReturnCorrectXSection("8TeV", 1);; // option is wrong fix when possible
  Double_t xSection8TeVErr                    = xSection8TeVErrDown*1e-3;
  TGraphAsymmErrors* graphPPCombPi0StatFitRatio_woXErr = NULL;
  TGraphAsymmErrors* graphPPCombPi0Stat           = (TGraphAsymmErrors*)fileReference->Get("Pi08TeV/graphInvCrossSectionPi0Comb8TeVAStatErr_yShifted");
  if(graphPPCombPi0Stat){
    graphPPCombPi0StatFitRatio_woXErr = CalculateGraphErrRatioToFit(graphPPCombPi0Stat,fitPPInvXsectionPi0);
    ProduceGraphAsymmWithoutXErrors(graphPPCombPi0StatFitRatio_woXErr);
  }
  TGraphAsymmErrors* graphPPCombPi0UncorrSys      = (TGraphAsymmErrors*)fileReference->Get("Pi08TeV/graphInvCrossSectionPi0Comb8TeVASysErr_yShifted"); // uncorr sys
  TGraphAsymmErrors* graphPPCombPi0FullSysFitRatio = NULL;
  TGraphAsymmErrors* graphPPCombPi0FullSys        = (TGraphAsymmErrors*)fileReference->Get("Pi08TeV/graphInvCrossSectionPi0Comb8TeVASysErr_yShifted"); 
  if(graphPPCombPi0FullSys)
      graphPPCombPi0FullSysFitRatio = CalculateGraphErrRatioToFit(graphPPCombPi0FullSys,fitPPInvXsectionPi0);
  TGraphAsymmErrors* graphPPCombPi0InterSys       = (TGraphAsymmErrors*)fileReference->Get("Pi08TeV/graphInvCrossSectionPi0Comb8TeVASysErr_yShifted");

  TGraphAsymmErrors* graphPPInvYieldCombPi0Stat   = (TGraphAsymmErrors*)graphPPCombPi0Stat->Clone("graphInvYieldStatErrComb");
  TGraphAsymmErrors* graphPPInvYieldCombPi0Sys    = (TGraphAsymmErrors*)graphPPCombPi0FullSys->Clone("graphInvYieldSystErrComb");
  graphPPInvYieldCombPi0Stat                      = ScaleGraph(graphPPCombPi0Stat,1/(xSection8TeV*recalcBarn));
  graphPPInvYieldCombPi0Sys                       = ScaleGraph(graphPPInvYieldCombPi0Sys,1/(xSection8TeV*recalcBarn));
  TGraphAsymmErrors* graphPPInvYieldCombPi0Stat_woXErr   = (TGraphAsymmErrors*)graphPPInvYieldCombPi0Stat->Clone("graphPPInvYieldCombPi0Stat_woXErr");
  ProduceGraphAsymmWithoutXErrors(graphPPInvYieldCombPi0Stat_woXErr);


  TF1* fitPPInvXsectionEta                        = (TF1*)fileReference->Get("Eta8TeV/TwoComponentModelFitEta");
  TGraphAsymmErrors* graphPPCombEtaStat           = (TGraphAsymmErrors*)fileReference->Get("Eta8TeV/graphInvCrossSectionEtaComb8TeVAStatErr_yShifted");
  TGraphAsymmErrors* graphPPCombEtaUncorrSys      = (TGraphAsymmErrors*)fileReference->Get("Eta8TeV/graphInvCrossSectionEtaComb8TeVASysErr_yShifted");
  TGraphAsymmErrors* graphPPCombEtaFullSys        = (TGraphAsymmErrors*)fileReference->Get("Eta8TeV/graphInvCrossSectionEtaComb8TeVASysErr_yShifted");
  TGraphAsymmErrors* graphPPCombEtaInterSys       = (TGraphAsymmErrors*)fileReference->Get("Eta8TeV/graphInvCrossSectionEtaComb8TeVASysErr_yShifted");

  TGraphAsymmErrors* graphPPInvYieldCombEtaStat   = (TGraphAsymmErrors*)graphPPCombEtaStat->Clone("graphInvYieldStatErrComb");
  TGraphAsymmErrors* graphPPInvYieldCombEtaSys    = (TGraphAsymmErrors*)graphPPCombEtaFullSys->Clone("graphInvYieldSystErrComb");
  graphPPInvYieldCombEtaStat                      = ScaleGraph(graphPPCombEtaStat,1/(xSection8TeV*recalcBarn));
  graphPPInvYieldCombEtaSys                       = ScaleGraph(graphPPInvYieldCombEtaSys,1/(xSection8TeV*recalcBarn));
  TGraphAsymmErrors* graphPPInvYieldCombEtaStat_woXErr   = (TGraphAsymmErrors*)graphPPInvYieldCombEtaStat->Clone("graphPPInvYieldCombEtaStat_woXErr");
  ProduceGraphAsymmWithoutXErrors(graphPPInvYieldCombEtaStat_woXErr);

  TFile* fileReferencePub                        = new TFile("ExternalInput/CombNeutralMesons/CombinedResultsPaperPP8TeV_2017_11_16.root");
  TGraphAsymmErrors* graphPPCombPi0StatFitRatioPub_woXErr           = (TGraphAsymmErrors*)fileReferencePub->Get("Pi08TeV/graphInvCrossSectionPi0Comb8TeVAStatErr_yShifted");
    graphPPCombPi0StatFitRatioPub_woXErr = CalculateGraphErrRatioToFit(graphPPCombPi0StatFitRatioPub_woXErr,fitPPInvXsectionPi0);
    ProduceGraphAsymmWithoutXErrors(graphPPCombPi0StatFitRatioPub_woXErr);
  TGraphAsymmErrors* graphPPCombPi0FullSysPub        = (TGraphAsymmErrors*)fileReferencePub->Get("Pi08TeV/graphInvCrossSectionPi0Comb8TeVASysErr_yShifted"); 
  if(graphPPCombPi0FullSysPub)
      graphPPCombPi0FullSysPub = CalculateGraphErrRatioToFit(graphPPCombPi0FullSysPub,fitPPInvXsectionPi0);

  // *******************************************************************************************************
  // ************************** Loading charged spectra     ************************************************
  // *******************************************************************************************************
  TFile* fileCHadrons                             = new TFile(fileNameCHadrons.Data());
  TFile* fileCParticles                           = new TFile(fileNameCParticles.Data());
  TGraphAsymmErrors* graphChargedHadronsStatpPb8TeV= (TGraphAsymmErrors*) fileCHadrons->Get("g1stat");
  graphChargedHadronsStatpPb8TeV                   = ScaleGraph(graphChargedHadronsStatpPb8TeV,67.6*1e-3*recalcBarn);
  TGraphAsymmErrors* graphChargedHadronsSyspPb8TeV = (TGraphAsymmErrors*) fileCHadrons->Get("g1syst");
  graphChargedHadronsSyspPb8TeV                    = ScaleGraph(graphChargedHadronsSyspPb8TeV,67.6*1e-3*recalcBarn);

  cout << __LINE__ << endl;
  TList* listStatChargedPionpPb8TeVprelim            = (TList*)fileCParticles->Get("Summed_Pion");
  TList* listSysChargedPionpPb8TeVprelim             = (TList*)fileCParticles->Get("Summed_Pion_Sys");
  TList* listStatChargedKaonpPb8TeVprelim            = (TList*)fileCParticles->Get("Summed_Kaon");
  TList* listSysChargedKaonpPb8TeVprelim             = (TList*)fileCParticles->Get("Summed_Kaon_Sys");
  TGraphAsymmErrors* graphChargedPionStatpPb8TeV   = new TGraphAsymmErrors( (TH1D*)listStatChargedPionpPb8TeVprelim->FindObject("hSpectraSummedPion_pp_Combined_MB"));
  graphChargedPionStatpPb8TeV                      = ConvertYieldGraph(graphChargedPionStatpPb8TeV, kTRUE, kTRUE, kFALSE, kFALSE);
  graphChargedPionStatpPb8TeV                      = ScaleGraph(graphChargedPionStatpPb8TeV,0.5*67.6*1e-3*recalcBarn);
  TGraphAsymmErrors* graphChargedPionSyspPb8TeV    = new TGraphAsymmErrors( (TH1D*)listSysChargedPionpPb8TeVprelim->FindObject("hSpectraSummedPion_pp_Combined_MB"));
  graphChargedPionSyspPb8TeV                       = ConvertYieldGraph(graphChargedPionSyspPb8TeV, kTRUE, kTRUE, kFALSE, kFALSE);
  graphChargedPionSyspPb8TeV                       = ScaleGraph(graphChargedPionSyspPb8TeV,0.5*67.6*1e-3*recalcBarn);
  TGraphAsymmErrors* graphChargedKaonStatpPb8TeV   = new TGraphAsymmErrors( (TH1D*)listStatChargedKaonpPb8TeVprelim->FindObject("hSpectraSummedKaon_pp_Combined_MB"));
  graphChargedKaonStatpPb8TeV                      = ConvertYieldGraph(graphChargedKaonStatpPb8TeV, kTRUE, kTRUE, kFALSE, kFALSE);
  graphChargedKaonStatpPb8TeV                      = ScaleGraph(graphChargedKaonStatpPb8TeV,0.5*67.6*1e-3*recalcBarn);
  TGraphAsymmErrors* graphChargedKaonSyspPb8TeV    = new TGraphAsymmErrors( (TH1D*)listSysChargedKaonpPb8TeVprelim->FindObject("hSpectraSummedKaon_pp_Combined_MB"));
  graphChargedKaonSyspPb8TeV                       = ConvertYieldGraph(graphChargedKaonSyspPb8TeV, kTRUE, kTRUE, kFALSE, kFALSE);
  graphChargedKaonSyspPb8TeV                       = ScaleGraph(graphChargedKaonSyspPb8TeV,0.5*67.6*1e-3*recalcBarn);


  // *******************************************************************************************************
  // ************************** Loading theory calculations ************************************************
  // *******************************************************************************************************
  TFile* fileEtaToPi                              = new TFile(fileNameEtaToPi0.Data());

    TH1F* histoHIJINGPi0                                = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoPi0SpecHIJING_MCGenpPb8TeV");
    TH1F* histoHIJINGEta                                = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoEtaSpecHIJING_MCGenpPb8TeV");
    TH1F* histoHIJINGEtaToPi0                           = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoEtaToPi0HIJING_MCGenpPb8TeV");
    TH1F* histoHIJINGPi0ToPiCh                          = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoPi0ToPiChHIJING_MCGenpPb8TeV");

      // EPOS JJ From LHC18b9bc
      TH1F* histoEPOSJJInvYield               = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoPi0SpecEPOSJJpPb8TeV_Reb");
      histoEPOSJJInvYield->GetXaxis()->SetRangeUser(10,maxPtPi0Theory-5);
      TH1F* histoEPOSJJInvYieldEta            = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoEtaSpecEPOSJJpPb8TeV_Reb");
      histoEPOSJJInvYieldEta->GetXaxis()->SetRangeUser(10,maxPtEta-1);

      TGraphErrors* graphEPOSJJInvYield       = new TGraphErrors((TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoPi0SpecEPOSJJpPb8TeV_Reb"));
      while(graphEPOSJJInvYield->GetX()[0] < 10) graphEPOSJJInvYield->RemovePoint(0);
      while(graphEPOSJJInvYield->GetX()[graphEPOSJJInvYield->GetN()-1] > maxPtPi0Theory-5) graphEPOSJJInvYield->RemovePoint(graphEPOSJJInvYield->GetN()-1);
      TGraphErrors* graphEPOSJJInvYieldEta    = new TGraphErrors((TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoEtaSpecEPOSJJpPb8TeV_Reb"));
      while(graphEPOSJJInvYieldEta->GetX()[0] < 10) graphEPOSJJInvYieldEta->RemovePoint(0);
      while(graphEPOSJJInvYieldEta->GetX()[graphEPOSJJInvYieldEta->GetN()-1] > maxPtEta-1) graphEPOSJJInvYieldEta->RemovePoint(graphEPOSJJInvYieldEta->GetN()-1);
      TH1F* histoEPOSJJEtaToPi0                  = (TH1F*) histoEPOSJJInvYieldEta->Clone("Pythia8EtaToPi0_EPOS");
      histoEPOSJJEtaToPi0->Divide(histoEPOSJJInvYield);
      histoEPOSJJEtaToPi0->GetXaxis()->SetRangeUser(0.4,maxPtEta-1);
      TGraphErrors* graphEPOSJJEtaToPi0           = new TGraphErrors(histoEPOSJJEtaToPi0);
      while(graphEPOSJJEtaToPi0->GetX()[0] < 10) graphEPOSJJEtaToPi0->RemovePoint(0);
      while(graphEPOSJJEtaToPi0->GetX()[graphEPOSJJEtaToPi0->GetN()-1] > maxPtEta-1) graphEPOSJJEtaToPi0->RemovePoint(graphEPOSJJEtaToPi0->GetN()-1);

      // DPMJET from LHC18f3bc
      TH1F* histoDPMJETInvYield               = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoPi0SpecDPMJETpPb8TeV_Reb");
      histoDPMJETInvYield->GetXaxis()->SetRangeUser(0.3,maxPtPi0Theory-5);
      TH1F* histoDPMJETInvYieldEta            = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoEtaSpecDPMJETpPb8TeV_Reb");
      histoDPMJETInvYieldEta->GetXaxis()->SetRangeUser(0.5,maxPtEta-1);

      TGraphErrors* graphDPMJETInvYield       = new TGraphErrors((TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoPi0SpecDPMJETpPb8TeV_Reb"));
      while(graphDPMJETInvYield->GetX()[0] < 0.3) graphDPMJETInvYield->RemovePoint(0);
      while(graphDPMJETInvYield->GetX()[graphDPMJETInvYield->GetN()-1] > maxPtPi0Theory-5) graphDPMJETInvYield->RemovePoint(graphDPMJETInvYield->GetN()-1);
      TGraphErrors* graphDPMJETInvYieldEta    = new TGraphErrors((TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoEtaSpecDPMJETpPb8TeV_Reb"));
      while(graphDPMJETInvYieldEta->GetX()[0] < 0.5) graphDPMJETInvYieldEta->RemovePoint(0);
      while(graphDPMJETInvYieldEta->GetX()[graphDPMJETInvYieldEta->GetN()-1] > maxPtEta-1) graphDPMJETInvYieldEta->RemovePoint(graphDPMJETInvYieldEta->GetN()-1);
      TH1F* histoDPMJETEtaToPi0                  = (TH1F*) histoDPMJETInvYieldEta->Clone("Pythia8EtaToPi0_DPMJET");
      histoDPMJETEtaToPi0->Divide(histoDPMJETInvYield);
      histoDPMJETEtaToPi0->GetXaxis()->SetRangeUser(0.4,maxPtEta-1);
      TGraphErrors* graphDPMJETEtaToPi0           = new TGraphErrors(histoDPMJETEtaToPi0);
      while(graphDPMJETEtaToPi0->GetX()[0] < 0.4) graphDPMJETEtaToPi0->RemovePoint(0);
      while(graphDPMJETEtaToPi0->GetX()[graphDPMJETEtaToPi0->GetN()-1] > maxPtEta-1) graphDPMJETEtaToPi0->RemovePoint(graphDPMJETEtaToPi0->GetN()-1);


      TH1F* histoPythia8EPPS16InvYield               = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoPi0Pythia8EPPS16pPb8TeV");
      histoPythia8EPPS16InvYield->Scale(208);
      histoPythia8EPPS16InvYield->GetXaxis()->SetRangeUser(0.3,maxPtPi0Theory-5);
      TH1F* histoPythia8EPPS16InvYieldEta            = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoEtaPythia8EPPS16pPb8TeV");
      histoPythia8EPPS16InvYieldEta->Scale(208);
      histoPythia8EPPS16InvYieldEta->GetXaxis()->SetRangeUser(0.5,maxPtEta-1);

      TGraphErrors* graphPythia8EPPS16InvYield       = new TGraphErrors((TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoPi0Pythia8EPPS16pPb8TeV"));
      graphPythia8EPPS16InvYield=ScaleGraph(graphPythia8EPPS16InvYield,208);
      while(graphPythia8EPPS16InvYield->GetX()[0] < 0.3) graphPythia8EPPS16InvYield->RemovePoint(0);
      while(graphPythia8EPPS16InvYield->GetX()[graphPythia8EPPS16InvYield->GetN()-1] > maxPtPi0Theory-5) graphPythia8EPPS16InvYield->RemovePoint(graphPythia8EPPS16InvYield->GetN()-1);
      TGraphErrors* graphPythia8EPPS16InvYieldEta    = new TGraphErrors((TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoEtaPythia8EPPS16pPb8TeV"));
      graphPythia8EPPS16InvYieldEta=ScaleGraph(graphPythia8EPPS16InvYieldEta,208);
      while(graphPythia8EPPS16InvYieldEta->GetX()[0] < 0.5) graphPythia8EPPS16InvYieldEta->RemovePoint(0);
      while(graphPythia8EPPS16InvYieldEta->GetX()[graphPythia8EPPS16InvYieldEta->GetN()-1] > maxPtEta-1) graphPythia8EPPS16InvYieldEta->RemovePoint(graphPythia8EPPS16InvYieldEta->GetN()-1);

      TH1F* histoPythia8nCTEQ15InvYield               = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoPi0Pythia8nCTEQ15pPb8TeV");
      histoPythia8nCTEQ15InvYield->Scale(208);
      histoPythia8nCTEQ15InvYield->GetXaxis()->SetRangeUser(0.3,maxPtPi0Theory-5);
      TH1F* histoPythia8nCTEQ15InvYieldEta            = (TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoEtaPythia8nCTEQ15pPb8TeV");
      histoPythia8nCTEQ15InvYieldEta->Scale(208);
      histoPythia8nCTEQ15InvYieldEta->GetXaxis()->SetRangeUser(0.5,maxPtEta-1);

      TGraphErrors* graphPythia8nCTEQ15InvYield       = new TGraphErrors((TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoPi0Pythia8nCTEQ15pPb8TeV"));
      graphPythia8nCTEQ15InvYield=ScaleGraph(graphPythia8nCTEQ15InvYield,208);
      while(graphPythia8nCTEQ15InvYield->GetX()[0] < 0.3) graphPythia8nCTEQ15InvYield->RemovePoint(0);
      while(graphPythia8nCTEQ15InvYield->GetX()[graphPythia8nCTEQ15InvYield->GetN()-1] > maxPtPi0Theory-5) graphPythia8nCTEQ15InvYield->RemovePoint(graphPythia8nCTEQ15InvYield->GetN()-1);
      TGraphErrors* graphPythia8nCTEQ15InvYieldEta    = new TGraphErrors((TH1F*) fileTheoryCompilationPPb->Get("pPb_8.16TeV/histoEtaPythia8nCTEQ15pPb8TeV"));
      graphPythia8nCTEQ15InvYieldEta=ScaleGraph(graphPythia8nCTEQ15InvYieldEta,208);
      while(graphPythia8nCTEQ15InvYieldEta->GetX()[0] < 0.5) graphPythia8nCTEQ15InvYieldEta->RemovePoint(0);
      while(graphPythia8nCTEQ15InvYieldEta->GetX()[graphPythia8nCTEQ15InvYieldEta->GetN()-1] > maxPtEta-1) graphPythia8nCTEQ15InvYieldEta->RemovePoint(graphPythia8nCTEQ15InvYieldEta->GetN()-1);


      TH1F* histoPythia8Monash2013InvSecpTh               = (TH1F*) fileTheoryCompilation->Get("histoInvSecpThPythia8Monash2013LegoPi08160GeV");
      histoPythia8Monash2013InvSecpTh->Scale(1/(2.1*1e12));
      histoPythia8Monash2013InvSecpTh->GetXaxis()->SetRangeUser(0.3,maxPtPi0Theory-5);
      TH1F* histoPythia8Monash2013InvSecpThEta            = (TH1F*) fileTheoryCompilation->Get("histoInvSecpThPythia8Monash2013LegoEta8160GeV");
      histoPythia8Monash2013InvSecpThEta->Scale(1/(2.1*1e12));
      histoPythia8Monash2013InvSecpThEta->GetXaxis()->SetRangeUser(0.5,maxPtEta-1);

      TGraphErrors* graphPythia8Monash2013InvSecpTh       = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecpThPythia8Monash2013LegoPi08160GeV"));
      while(graphPythia8Monash2013InvSecpTh->GetX()[0] < 0.3) graphPythia8Monash2013InvSecpTh->RemovePoint(0);
      while(graphPythia8Monash2013InvSecpTh->GetX()[graphPythia8Monash2013InvSecpTh->GetN()-1] > maxPtPi0Theory-5) graphPythia8Monash2013InvSecpTh->RemovePoint(graphPythia8Monash2013InvSecpTh->GetN()-1);
      TGraphErrors* graphPythia8Monash2013InvSecpThEta    = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecpThPythia8Monash2013LegoEta8160GeV"));
      while(graphPythia8Monash2013InvSecpThEta->GetX()[0] < 0.5) graphPythia8Monash2013InvSecpThEta->RemovePoint(0);
      while(graphPythia8Monash2013InvSecpThEta->GetX()[graphPythia8Monash2013InvSecpThEta->GetN()-1] > maxPtEta-1) graphPythia8Monash2013InvSecpThEta->RemovePoint(graphPythia8Monash2013InvSecpThEta->GetN()-1);


    // pythia for pp reference (monash 2013)
        TH1F* histoPythia8Monash2013InvXSection                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi08TeV");
        histoPythia8Monash2013InvXSection->Scale(1/(xSection8TeV*recalcBarn));
        histoPythia8Monash2013InvXSection->GetXaxis()->SetRangeUser(0.3,200);
        TH1F* histoPythia8Monash2013InvXSectionEta                    = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta8TeV");
        histoPythia8Monash2013InvXSectionEta->Scale(1/(xSection8TeV*recalcBarn));
        histoPythia8Monash2013InvXSectionEta->GetXaxis()->SetRangeUser(0.4,35);
        TGraphErrors* graphPythia8Monash2013InvXSection               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi08TeV"));
        graphPythia8Monash2013InvXSection = ScaleGraph(graphPythia8Monash2013InvXSection, 1/(xSection8TeV*recalcBarn));
        while(graphPythia8Monash2013InvXSection->GetX()[0] < 0.3) graphPythia8Monash2013InvXSection->RemovePoint(0);
        // while(graphPythia8Monash2013InvXSection->GetX()[graphPythia8Monash2013InvXSection->GetN()-1] > 35.) graphPythia8Monash2013InvXSection->RemovePoint(graphPythia8Monash2013InvXSection->GetN()-1);
        TGraphErrors* graphPythia8Monash2013InvXSectionEta            = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta8TeV"));
        graphPythia8Monash2013InvXSectionEta = ScaleGraph(graphPythia8Monash2013InvXSectionEta, 1/(xSection8TeV*recalcBarn));
        while(graphPythia8Monash2013InvXSectionEta->GetX()[0] < 0.4) graphPythia8Monash2013InvXSectionEta->RemovePoint(0);
        while(graphPythia8Monash2013InvXSectionEta->GetX()[graphPythia8Monash2013InvXSectionEta->GetN()-1] > 35.) graphPythia8Monash2013InvXSectionEta->RemovePoint(graphPythia8Monash2013InvXSectionEta->GetN()-1);
      // *******************************************************************************************************
      // NLO calc
      TGraphAsymmErrors* graphPi0DSS14            = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcDSS14InvSecPi05023GeV");
      while (graphPi0DSS14->GetX()[graphPi0DSS14->GetN()-1] > maxPtPi0 ) graphPi0DSS14->RemovePoint(graphPi0DSS14->GetN()-1);
      while(graphPi0DSS14->GetX()[0] < 2.0) graphPi0DSS14->RemovePoint(0);
      // graphPi0DSS14                       = ScaleGraph(graphPi0DSS14,7.118/(xSection8TeVINEL*recalcBarn));
      graphPi0DSS14                       = ScaleGraph(graphPi0DSS14,7.118/(xSection8TeVV0AND*recalcBarn));

      TGraphAsymmErrors* graphEtaAESSS            = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcAESSSInvSecEta5023GeV");

      TGraphAsymmErrors* graphEtaToPi02760GeV     = (TGraphAsymmErrors*) fileEtaToPi->Get("Alice2760GeV");
      ProduceGraphAsymmWithoutXErrors(graphEtaToPi02760GeV);
      TGraphAsymmErrors* graphEtaToPi08000GeV     = (TGraphAsymmErrors*) fileEtaToPi->Get("Alice8TeV");
      ProduceGraphAsymmWithoutXErrors(graphEtaToPi08000GeV);
  cout << __LINE__ << endl;

      // *******************************************************************************************************

      TGraph* graphNLOCalcPi0MuHalf               = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSS14InvSecPi0MuHalf5023GeV");
      graphNLOCalcPi0MuHalf                       = ScaleGraph(graphNLOCalcPi0MuHalf,7.118/(xSection8TeVINEL*recalcBarn));
      TGraph* graphNLOCalcPi0MuOne                = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSS14InvSecPi0MuOne5023GeV");
      graphNLOCalcPi0MuOne                       = ScaleGraph(graphNLOCalcPi0MuOne,7.118/(xSection8TeVINEL*recalcBarn));
      TGraph* graphNLOCalcPi0MuTwo                = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSS14InvSecPi0MuTwo5023GeV");
      graphNLOCalcPi0MuTwo                       = ScaleGraph(graphNLOCalcPi0MuTwo,7.118/(xSection8TeVINEL*recalcBarn));
      TGraph* graphNLOCalcEtaMuHalf               = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcAESSSInvSecEtaMuHalf5023GeV");
      TGraph* graphNLOCalcEtaMuOne                = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcAESSSInvSecEtaMuOne5023GeV");
      TGraph* graphNLOCalcEtaMuTwo                = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcAESSSInvSecEtaMuTwo5023GeV");

      while (graphNLOCalcPi0MuHalf->GetX()[graphNLOCalcPi0MuHalf->GetN()-1] > maxPtPi0-5 )
          graphNLOCalcPi0MuHalf->RemovePoint(graphNLOCalcPi0MuHalf->GetN()-1);
      while (graphNLOCalcPi0MuOne->GetX()[graphNLOCalcPi0MuOne->GetN()-1] > maxPtPi0-5 )
          graphNLOCalcPi0MuOne->RemovePoint(graphNLOCalcPi0MuOne->GetN()-1);
      while (graphNLOCalcPi0MuTwo->GetX()[graphNLOCalcPi0MuTwo->GetN()-1] > maxPtPi0-5 )
          graphNLOCalcPi0MuTwo->RemovePoint(graphNLOCalcPi0MuTwo->GetN()-1);
      while (graphNLOCalcEtaMuHalf->GetX()[graphNLOCalcEtaMuHalf->GetN()-1] > maxPtEta-1 )
          graphNLOCalcEtaMuHalf->RemovePoint(graphNLOCalcEtaMuHalf->GetN()-1);
      while (graphNLOCalcEtaMuOne->GetX()[graphNLOCalcEtaMuOne->GetN()-1] > maxPtEta-1 )
          graphNLOCalcEtaMuOne->RemovePoint(graphNLOCalcEtaMuOne->GetN()-1);
      while (graphNLOCalcEtaMuTwo->GetX()[graphNLOCalcEtaMuTwo->GetN()-1] > maxPtEta-1 )
          graphNLOCalcEtaMuTwo->RemovePoint(graphNLOCalcEtaMuTwo->GetN()-1);

      TGraphAsymmErrors* graphNLOEtaToPi0         = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcEtaOverPi05023GeV_AESSS_DSS14");
      while (graphNLOEtaToPi0->GetX()[graphNLOEtaToPi0->GetN()-1] > maxPtEta-1 ) graphNLOEtaToPi0->RemovePoint(graphNLOEtaToPi0->GetN()-1);
      while(graphNLOEtaToPi0->GetX()[0] < 2.0) graphNLOEtaToPi0->RemovePoint(0);
  cout << __LINE__ << endl;
  // *******************************************************************************************************
  // ************************** Combination of different measurements **************************************
  // *******************************************************************************************************

  // Int_t nBinsPi0 = 34;//36
  Int_t nBinsPi0 = 58;//36
  Double_t xPtLimits[59]=  {
    0.0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,
    1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6,
    3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,
    8.5, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
    18.0, 20.0, 22.0, 26.0, 30.0, 35.0, 40.0 ,45.0, 50.0,
    55.0, 60.0, 65.0, 70.0, 80.0, 100.0, 125.0, 150.0, 175.0, 200.0
    // 0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4,
    // 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
    // 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
    // 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 24.0,
    // 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 45.0, 50.0,
    // 55.0, 60.0, 65.0, 70.0, 80.0, 100.0, 125.0, 150.0, 175.0, 200.0
  };
  
  // Int_t nBinsPi0 = 60;//36
  // Double_t xPtLimits[61]=  {
  //   0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4,
  //   1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
  //   3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0, 12.0,
  //   13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 24.0,
  //   26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 45.0, 50.0,
  //   55.0, 60.0, 65.0, 70.0, 80.0, 100.0, 125.0, 150.0, 175.0, 200.0
  // };


  Int_t nBinsEta = 25;//15
  Double_t xPtLimitsEta[29] =  {
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5,
    4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0,
    20.0, 25.0, 30.0, 35.0, 40.0,50.,60.
  };
  Double_t xPtLimitsEtaToPi0[29] =  {
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5,
    4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0,
    20.0, 25.0, 30.0, 35.0, 40.0,50.,60.
  };
  // Double_t xPtLimitsEta[29] =  {
  //   0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.4, 2.8, 3.2, 3.6,
  //   4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
  //   18.0, 20.0, 25.0, 30.0, 35.0, 40.0,50.,60.
  // };
  // Double_t xPtLimitsEtaToPi0[29] =  {
  //   0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.4, 2.8, 3.2, 3.6,
  //   4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
  //   18.0, 20.0, 25.0, 30.0, 35.0, 40.0,50.,60.
  // };


  // *******************************************************************************************************
  // ************************** Combination of different measurements **************************************
  // *******************************************************************************************************
  // REMARKS:
  //       - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
  //       - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
  //       - currently only PCM-EMCAL vs others fully implemeted energy independent
  //       - extendable to other energies
  //       - offsets have to be determined manually, see cout's in shell from combination function


  // definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
  // for statistical error and final value from respective method
  cout << __LINE__ << endl;

  TH1D* statErrorCollection[13];
  TH1D* statErrorCollectionEta[13];
  TH1D* statErrorCollectionEtaToPi0[13];
  for (Int_t i = 0; i< totalNSets; i++){
    statErrorCollection[i]                        = NULL;
    statErrorCollectionEta[i]                     = NULL;
    statErrorCollectionEtaToPi0[i]                = NULL;
  }
  for (Int_t i = 0; i< totalNSets; i++){
    if(directoryPi0[i]&&histoPi0InvYield[i]) statErrorCollection[i]    = (TH1D*)histoPi0InvYield[i]->Clone(Form("statErr%sPi0",nameMeasGlobalPlot[i].Data()));
    if(directoryEta[i]&&histoEtaInvYield[i]) statErrorCollectionEta[i] = (TH1D*)histoEtaInvYield[i]->Clone(Form("statErr%sEta",nameMeasGlobalPlot[i].Data()));
    if(directoryEta[i]&&histoEtaToPi0Stat[i]) statErrorCollectionEtaToPi0[i]= (TH1D*)histoEtaToPi0Stat[i]->Clone(Form("statErr%sEtaToPi0",nameMeasGlobalPlot[i].Data()));
  }
  cout << __LINE__ << endl;
  // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
  // for systematic error from respective method
  TGraphAsymmErrors* sysErrorCollection[13];
  TGraphAsymmErrors* sysErrorCollectionEta[13];
  TGraphAsymmErrors* sysErrorCollectionEtaToPi0[13];
  for (Int_t i = 0; i< totalNSets; i++){
      sysErrorCollection[i]                       = NULL;
      sysErrorCollectionEta[i]                    = NULL;
      sysErrorCollectionEtaToPi0[i]               = NULL;
  }
  for (Int_t i = 0; i< totalNSets; i++){
      if(directoryPi0[i]&&graphPi0InvYieldSys[i]) sysErrorCollection[i]   = (TGraphAsymmErrors*)graphPi0InvYieldSys[i]->Clone(Form("sysErr%sPi0",nameMeasGlobalPlot[i].Data()));
      if(directoryEta[i]&&graphEtaInvYieldSys[i]) sysErrorCollectionEta[i]= (TGraphAsymmErrors*)graphEtaInvYieldSys[i]->Clone(Form("sysErr%sEta",nameMeasGlobalPlot[i].Data()));
      if(directoryEta[i]&&graphEtaToPi0Sys[i]) sysErrorCollectionEtaToPi0[i] = (TGraphAsymmErrors*)graphEtaToPi0Sys[i]->Clone(Form("sysErr%sEtaToPi0",nameMeasGlobalPlot[i].Data()));
  }

cout << __LINE__ << endl;

  // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match
  //                                                 PCM,PHOS,EMC,PCMPHOS,PCMEMC,        mEMC
  Int_t offSets[13]                               = { 0,    2,  0,     2,     0, 0,0,0,8,   0,0,0,0};
  Int_t offSetsSys[13]                            = { 1,    2,  7,     3,     5, 0,0,0,6,  38,0,0,0};
  // Int_t offSets[13]                               = { 0,    2,  7,     2,     5, 0,0,0,   8,3,0,0,0};
  // Int_t offSetsSys[13]                            = { 1,    2,  9,     3,     6, 0,0,0,   9,36,0,0,0};

  Int_t offSetPi0Shifting[13]                     = { 0,  0,  0,  0,  0,
                                                      0,  0,  0,  0,  0,
                                                      0,  0,  0 };
  Int_t nComBinsPi0Shifting[13]                   = { 0,  0,  0,  0,  0,
                                                      0,  0,  0,  0,  0,
                                                      0,  0,  0 };
  //                                                 PCM,PHOS,EMC,PCMPHOS,PCMEMC,         mEMC
  Int_t offSetsEta[13]                            = { -1,  1,  0,     2,      0, 0,0,0,   4,0,0,0,0};
  Int_t offSetsSysEta[13]                         = { 0,   1,  4,     3,      1, 0,0,0,   9,0,0,0,0};
  // Int_t offSetsEtaRpA[13]                         = { -4,  1,  0,     2,     -4, 0,0,0,   4,0,0,0,0};
  // Int_t offSetsSysEtaRpA[13]                      = { 0,   1,  5,     3,      2, 0,0,0,   9,0,0,0,0};
  // Int_t offSetsEta[13]                            = { -3,   1,  2,     2,      -3, 0,0,0,   4,0,0,0,0};
  // Int_t offSetsSysEta[13]                         = { 0,    1,  5,     3,      3, 0,0,0,   9,0,0,0,0};
  Int_t offSetEtaShifting[13]                     = { 0,  0,  0,  0,  0,
                                                      0,  0,  0,  0,  0,
                                                      0,  0,  0 };
  Int_t nComBinsEtaShifting[13]                   = { 0,  0,  0,  0,  0,
                                                      0,  0,  0,  0,  0,
                                                      0,  0,  0 };
  //                                                 PCM,PHOS,EMC,PCMPHOS,PCMEMC,         mEMC
  Int_t offSetsEtaToPi0[13]                       = { -1,  1,  0,     2,      0, 0,0,0,   4,0,0,0,0};
  Int_t offSetsSysEtaToPi0[13]                    = { 0,   1,  4,     3,      1, 0,0,0,   9,0,0,0,0};
  // Int_t offSetsEtaToPi0[13]                       = { -3,   1,  2,     2,      -3, 0,0,0,   4,0,0,0,0};
  // Int_t offSetsSysEtaToPi0[13]                    = { 0,    1,  5,     3,      3, 0,0,0,   9,0,0,0,0};
  
  //**********************************************************************************************************************
  //**********************************************************************************************************************
  //**********************************************************************************************************************

  TH1D* statErrorRelCollection[13];
  TH1D* statErrorRelCollectionEta[13];
  TH1D* statErrorRelCollectionEtaToPi0[13];

  for (Int_t i = 0; i < totalNSets; i++){
      if (statErrorCollection[i])
        statErrorRelCollection[i]                 = CalculateRelErrUpTH1D( statErrorCollection[i], Form("relativeStatError_%s", nameMeasGlobalPlot[i].Data()));
      else statErrorRelCollection[i]              = NULL;
      if (statErrorCollectionEta[i])
        statErrorRelCollectionEta[i]              = CalculateRelErrUpTH1D( statErrorCollectionEta[i], Form("relativeStatError_%s", nameMeasGlobalPlot[i].Data()));
      else statErrorRelCollectionEta[i]           = NULL;
      if (statErrorCollectionEtaToPi0[i])
        statErrorRelCollectionEtaToPi0[i]         = CalculateRelErrUpTH1D( statErrorCollectionEtaToPi0[i], Form("relativeStatError_%s", nameMeasGlobalPlot[i].Data()));
      else statErrorRelCollectionEtaToPi0[i]      = NULL;
  }
  cout << __LINE__ << endl;

  TGraphAsymmErrors* sysErrorRelCollection[13];
  TGraphAsymmErrors* sysErrorRelCollectionEta[13];
  TGraphAsymmErrors* sysErrorRelCollectionEtaToPi0[13];

  for (Int_t i = 0; i < totalNSets; i++){
      if (sysErrorCollection[i])
        sysErrorRelCollection[i]                  = CalculateRelErrUpAsymmGraph( sysErrorCollection[i], Form("relativeSysError_%s", nameMeasGlobalPlot[i].Data()));
      else sysErrorRelCollection[i]               = NULL;
      if (sysErrorCollectionEta[i])
        sysErrorRelCollectionEta[i]               = CalculateRelErrUpAsymmGraph( sysErrorCollectionEta[i], Form("relativeSysError_%s", nameMeasGlobalPlot[i].Data()));
      else sysErrorRelCollectionEta[i]            = NULL;
      if (sysErrorCollectionEtaToPi0[i])
        sysErrorRelCollectionEtaToPi0[i]          = CalculateRelErrUpAsymmGraph( sysErrorCollectionEtaToPi0[i], Form("relativeSysError_%s", nameMeasGlobalPlot[i].Data()));
      else sysErrorRelCollectionEtaToPi0[i]       = NULL;
  }
  cout << __LINE__ << endl;

  TGraph* graphWeights[13];
  TGraph* graphWeightsEta[13];
  TGraph* graphWeightsEtaToPi0[13];
  for (Int_t i = 0; i< totalNSets; i++){
      graphWeights[i]                             = NULL;
      graphWeightsEta[i]                          = NULL;
      graphWeightsEtaToPi0[i]                     = NULL;
  }

  // Declaration & calculation of combined spectrum
  TString fileNameOutputWeighting                 = Form("%s/Weighting.dat",outputDir.Data());
  TString fileNameOutputWeightingEta              = Form("%s/WeightingEta.dat",outputDir.Data());
  TString fileNameOutputWeightingEtaToPi0         = Form("%s/WeightingEtaToPi0.dat",outputDir.Data());

  TGraphAsymmErrors* graphCombPi0InvYieldStat  = NULL;
  TGraphAsymmErrors* graphCombPi0InvYieldSys   = NULL;
  TGraphAsymmErrors* graphCombPi0InvYieldTot   = CombinePtPointsSpectraFullCorrMat( statErrorCollection, sysErrorCollection,
                                                                                        xPtLimits, nBinsPi0,
                                                                                        offSets, offSetsSys,
                                                                                        graphCombPi0InvYieldStat, graphCombPi0InvYieldSys,
                                                                                        fileNameOutputWeighting,"pPb_8TeV", "Pi0", kTRUE,
                                                                                        0x0, ""
                                                                                      );

  // return;
  cout << __LINE__ << endl;
  while (graphCombPi0InvYieldStat->GetX()[0] < minPtPi0Meas ) graphCombPi0InvYieldStat->RemovePoint(0);
  while (graphCombPi0InvYieldStat->GetX()[graphCombPi0InvYieldStat->GetN()-1] > maxPtPi0Meas) graphCombPi0InvYieldStat->RemovePoint(graphCombPi0InvYieldStat->GetN()-1);
  while (graphCombPi0InvYieldSys->GetX()[0] < minPtPi0Meas ) graphCombPi0InvYieldSys->RemovePoint(0);
  while (graphCombPi0InvYieldSys->GetX()[graphCombPi0InvYieldSys->GetN()-1] > maxPtPi0Meas) graphCombPi0InvYieldSys->RemovePoint(graphCombPi0InvYieldSys->GetN()-1);
  while (graphCombPi0InvYieldTot->GetX()[0] < minPtPi0Meas ) graphCombPi0InvYieldTot->RemovePoint(0);
  while (graphCombPi0InvYieldTot->GetX()[graphCombPi0InvYieldTot->GetN()-1] > maxPtPi0Meas) graphCombPi0InvYieldTot->RemovePoint(graphCombPi0InvYieldTot->GetN()-1);
  if(doOutput) graphCombPi0InvYieldStat->Print();
  
  TGraphAsymmErrors* graphCombEtaInvYieldStat  = NULL;
  TGraphAsymmErrors* graphCombEtaInvYieldSys   = NULL;
  TGraphAsymmErrors* graphCombEtaInvYieldTot   = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEta, sysErrorCollectionEta,
                                                                                        xPtLimitsEta, nBinsEta,
                                                                                        offSetsEta, offSetsSysEta,
                                                                                        graphCombEtaInvYieldStat, graphCombEtaInvYieldSys,
                                                                                        fileNameOutputWeightingEta,"pPb_8TeV", "Eta", kTRUE,
                                                                                        0x0, ""
                                                                                      );
  // return;
  cout << __LINE__ << endl;
  while (graphCombEtaInvYieldStat->GetX()[0] < minPtEtaMeas ) graphCombEtaInvYieldStat->RemovePoint(0);
  while (graphCombEtaInvYieldStat->GetX()[graphCombEtaInvYieldStat->GetN()-1] > maxPtEtaMeas) graphCombEtaInvYieldStat->RemovePoint(graphCombEtaInvYieldStat->GetN()-1);
  while (graphCombEtaInvYieldSys->GetX()[0] < minPtEtaMeas ) graphCombEtaInvYieldSys->RemovePoint(0);
  while (graphCombEtaInvYieldSys->GetX()[graphCombEtaInvYieldSys->GetN()-1] > maxPtEtaMeas) graphCombEtaInvYieldSys->RemovePoint(graphCombEtaInvYieldSys->GetN()-1);
  while (graphCombEtaInvYieldTot->GetX()[0] < minPtEtaMeas ) graphCombEtaInvYieldTot->RemovePoint(0);
  while (graphCombEtaInvYieldTot->GetX()[graphCombEtaInvYieldTot->GetN()-1] > maxPtEtaMeas) graphCombEtaInvYieldTot->RemovePoint(graphCombEtaInvYieldTot->GetN()-1);
  if(doOutput) graphCombEtaInvYieldStat->Print();
  
  
  TGraphAsymmErrors* graphCombEtaToPi0Stat        = NULL;
  TGraphAsymmErrors* graphCombEtaToPi0Sys         = NULL;
  TGraphAsymmErrors* graphCombEtaToPi0Tot         = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEtaToPi0, sysErrorCollectionEtaToPi0,
                                                                                        xPtLimitsEtaToPi0, nBinsEta,
                                                                                        offSetsEtaToPi0, offSetsSysEtaToPi0,
                                                                                        graphCombEtaToPi0Stat, graphCombEtaToPi0Sys,
                                                                                        fileNameOutputWeightingEtaToPi0,"pPb_8TeV", "EtaToPi0", kTRUE,
                                                                                        0x0, ""
                                                                                      );
  // return;
  cout << __LINE__ << endl;
  while (graphCombEtaToPi0Stat->GetX()[0] < minPtEtaMeas ) graphCombEtaToPi0Stat->RemovePoint(0);
  while (graphCombEtaToPi0Stat->GetX()[graphCombEtaToPi0Stat->GetN()-1] > maxPtEtaMeas) graphCombEtaToPi0Stat->RemovePoint(graphCombEtaToPi0Stat->GetN()-1);
  while (graphCombEtaToPi0Sys->GetX()[0] < minPtEtaMeas ) graphCombEtaToPi0Sys->RemovePoint(0);
  while (graphCombEtaToPi0Sys->GetX()[graphCombEtaToPi0Sys->GetN()-1] > maxPtEtaMeas) graphCombEtaToPi0Sys->RemovePoint(graphCombEtaToPi0Sys->GetN()-1);
  while (graphCombEtaToPi0Tot->GetX()[0] < minPtEtaMeas ) graphCombEtaToPi0Tot->RemovePoint(0);
  while (graphCombEtaToPi0Tot->GetX()[graphCombEtaToPi0Tot->GetN()-1] > maxPtEtaMeas) graphCombEtaToPi0Tot->RemovePoint(graphCombEtaToPi0Tot->GetN()-1);
  if(doOutput) graphCombEtaToPi0Stat->Print();
  
  
  
  //**********************************************************************************************************************
  //**********************************************************************************************************************
  //**********************************************************************************************************************
  // plot weights + unc. for pi0
  //**********************************************************************************************************************
  //**********************************************************************************************************************
  //**********************************************************************************************************************

  // Reading weights from output file for plotting
  ifstream fileWeightsRead;
  fileWeightsRead.open(fileNameOutputWeighting,ios_base::in);
  cout << "reading" << fileNameOutputWeighting << endl;
  Double_t xValuesRead[80];
  Double_t weightsRead[13][80];
  Int_t availableMeas[13]                         = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  Int_t nMeasSetPi0                                  = numbersofmeasPi0;
  Int_t nPtBinsRead                               = 0;
  while(!fileWeightsRead.eof() && nPtBinsRead < 80){
      TString garbage                             = "";
      if (nPtBinsRead == 0){
          fileWeightsRead >> garbage ;
          for (Int_t i = 0; i < nMeasSetPi0; i++){
              fileWeightsRead >> availableMeas[i] ;
          }
          cout << "read following measurements: ";
          for (Int_t i = 0; i < 13; i++){
              cout << availableMeas[i] << "\t" ;
          }
          cout << endl;
      } else {
          fileWeightsRead >> xValuesRead[nPtBinsRead-1];
          for (Int_t i = 0; i < nMeasSetPi0; i++){
              fileWeightsRead >> weightsRead[availableMeas[i]][nPtBinsRead-1] ;
          }
          cout << "read: "<<  nPtBinsRead << "\t"<< xValuesRead[nPtBinsRead-1] << "\t" ;
          for (Int_t i = 0; i < nMeasSetPi0; i++){
              cout << weightsRead[availableMeas[i]][nPtBinsRead-1] << "\t";
          }
          cout << endl;
      }
      nPtBinsRead++;
  }
  nPtBinsRead = nPtBinsRead-2 ;
  fileWeightsRead.close();

  for (Int_t i = 0; i < nMeasSetPi0; i++){
      graphWeights[availableMeas[i]]              = new TGraph(nPtBinsRead,xValuesRead,weightsRead[availableMeas[i]]);
      Int_t bin = 0;
      for (Int_t n = 0; n< nPtBinsRead; n++){
          if (graphWeights[availableMeas[i]]->GetY()[bin] == 0) graphWeights[availableMeas[i]]->RemovePoint(bin);
          else bin++;
      }
  }

  //**********************************************************************************************************************
  //******************************************* Plotting weights for pi0  ************************************************
  //**********************************************************************************************************************

  Int_t textSizeLabelsPixel                   = 900*0.04;

  TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
  canvasWeights->SetLogx();

  TH2F * histo2DWeights = new TH2F("histo2DWeights","histo2DWeights",11000,minPtPi0,maxPtPi0,1000,-0.5,1.1);
  SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DWeights->GetXaxis()->SetMoreLogLabels();
  histo2DWeights->GetXaxis()->SetNoExponent();
  canvasWeights->cd();
  histo2DWeights->Draw("copy");

      TLegend* legendAccWeights                   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0*1.35), 32);
      for (Int_t i = 0; i < nMeasSetPi0; i++){
          DrawGammaSetMarkerTGraph(graphWeights[availableMeas[i]],
                                  markerStyleDet[availableMeas[i]],
                                  markerSizeDet[availableMeas[i]]*0.5,
                                  colorDet[availableMeas[i]] ,
                                  colorDet[availableMeas[i]]);
          graphWeights[availableMeas[i]]->Draw("p,same,e1");
          legendAccWeights->AddEntry(graphWeights[availableMeas[i]],nameMeasGlobalPlot[availableMeas[i]],"p");
      }
      legendAccWeights->Draw();
      TLatex *labelWeightsEnergy                  = new TLatex(0.7,0.20,collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
      labelWeightsEnergy->SetTextFont(43);
      labelWeightsEnergy->Draw();
      TLatex *labelWeightsPi0                     = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
      labelWeightsPi0->SetTextFont(43);
      labelWeightsPi0->Draw();

      DrawGammaLines(minPtPi0, maxPtPi0 , 0.5, 0.5,0.1, kGray, 7);
      DrawGammaLines(minPtPi0, maxPtPi0 , 0.4, 0.4,0.1, kGray, 1);
      DrawGammaLines(minPtPi0, maxPtPi0 , 0.3, 0.3,0.1, kGray, 7);
      DrawGammaLines(minPtPi0, maxPtPi0 , 0.2, 0.2,0.1, kGray, 3);

  canvasWeights->SaveAs(Form("%s/Pi0_WeightsMethods.%s",outputDir.Data(),suffix.Data()));


  //*********************************************************************************************************************
  //************************************ Visualize relative errors ******************************************************
  //*********************************************************************************************************************

  TCanvas* canvasRelSysErr                        = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
  canvasRelSysErr->SetLogx();

  TH2F * histo2DRelSysErr = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,minPtPi0,maxPtPi0,1000,0,80.5);
  SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
  histo2DRelSysErr->GetXaxis()->SetNoExponent();
  histo2DRelSysErr->GetYaxis()->SetRangeUser(0,45.5);
  histo2DRelSysErr->Draw("copy");

      TLegend* legendRelSysErr                    = GetAndSetLegend2(0.62, 0.94-(0.035*nMeasSetPi0*1.35), 0.95, 0.94, 32);
      for (Int_t i = 0; i < nMeasSetPi0; i++){
          DrawGammaSetMarkerTGraph(sysErrorRelCollection[availableMeas[i]], markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]],
                                  colorDet[availableMeas[i]]);
          sysErrorRelCollection[availableMeas[i]]->Draw("p,same,e1");
          legendRelSysErr->AddEntry(sysErrorRelCollection[availableMeas[i]],nameMeasGlobalPlot[availableMeas[i]],"p");
      }
      legendRelSysErr->Draw();

      TLatex *labelRelSysErrEnergy                = new TLatex(0.15,0.89,collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel,4);
      labelRelSysErrEnergy->SetTextFont(43);
      labelRelSysErrEnergy->Draw();
      TLatex *labelRelSysErrPi0                   = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelSysErrPi0, 0.85*textSizeLabelsPixel,4);
      labelRelSysErrPi0->SetTextFont(43);
      labelRelSysErrPi0->Draw();

  canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErr.%s",outputDir.Data(),suffix.Data()));

  //*********************************************************************************************************************
  //************************************ Visualize relative errors ******************************************************
  //*********************************************************************************************************************

  TCanvas* canvasRelStatErr                   = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
  canvasRelStatErr->SetLogx();

  TH2F * histo2DRelStatErr = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,minPtPi0,maxPtPi0,1000,0,80.5);
  SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
  histo2DRelStatErr->GetXaxis()->SetNoExponent();
  histo2DRelStatErr->GetYaxis()->SetRangeUser(0,45.5);
  histo2DRelStatErr->Draw("copy");

      TLegend* legendRelStatErr               = GetAndSetLegend2(0.14, 0.94-(0.035*nMeasSetPi0*1.35), 0.45, 0.94, 32);
      for (Int_t i = 0; i < nMeasSetPi0; i++){
          if (availableMeas[i]== 2){
              DrawGammaSetMarker(statErrorRelCollection[availableMeas[i]], markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]] ,
                          colorDet[availableMeas[i]]);
              TGraphAsymmErrors* graphDummy   = new TGraphAsymmErrors(statErrorRelCollection[availableMeas[i]]);
              DrawGammaSetMarkerTGraphAsym(graphDummy, markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]],
                                  colorDet[availableMeas[i]]);
              graphDummy->Draw("same,p,x0");
              legendRelStatErr->AddEntry(graphDummy,nameMeasGlobalPlot[availableMeas[i]],"p");

          } else {
              DrawGammaSetMarker(statErrorRelCollection[availableMeas[i]], markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]] ,
                          colorDet[availableMeas[i]]);
              statErrorRelCollection[availableMeas[i]]->Draw("p,same,e1");
              legendRelStatErr->AddEntry(statErrorRelCollection[availableMeas[i]],nameMeasGlobalPlot[availableMeas[i]],"p");

          }
      }
      legendRelStatErr->Draw();

      TLatex *labelRelStatErrEnergy           = new TLatex(0.75,0.89,collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelRelStatErrEnergy, 0.85*textSizeLabelsPixel,4);
      labelRelStatErrEnergy->SetTextFont(43);
      labelRelStatErrEnergy->Draw();
      TLatex *labelRelStatErrPi0              = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelStatErrPi0, 0.85*textSizeLabelsPixel,4);
      labelRelStatErrPi0->SetTextFont(43);
      labelRelStatErrPi0->Draw();

  canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErr.%s",outputDir.Data(),suffix.Data()));


  //**********************************************************************************************************************
  //**********************************************************************************************************************

  TGraphAsymmErrors* graphCombPi0InvYieldRelStat     = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldStat, "graphCombPi0InvYieldRelStat");
  TGraphAsymmErrors* graphCombPi0InvYieldRelSys      = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldSys, "graphCombPi0InvYieldRelSys");
  TGraphAsymmErrors* graphCombPi0InvYieldRelTot      = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldTot, "graphCombPi0InvYieldRelTot");

  TCanvas* canvasRelTotErr                    = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
  canvasRelTotErr->SetLogx();

  TH2F * histo2DRelTotErr;
  histo2DRelTotErr                            = new TH2F("histo2DRelTotErr","histo2DRelTotErr",11000,minPtPi0,maxPtPi0,1000,0,57.5);
  SetStyleHistoTH2ForGraphs(histo2DRelTotErr, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DRelTotErr->GetXaxis()->SetMoreLogLabels();
  histo2DRelTotErr->GetXaxis()->SetNoExponent();

  histo2DRelTotErr->GetYaxis()->SetRangeUser(0,57.5);
  histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
  histo2DRelTotErr->Draw("copy");

  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombPi0InvYieldRelTot->Draw("p,same,e1");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
  graphCombPi0InvYieldRelStat->Draw("l,x0,same,e1");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
  graphCombPi0InvYieldRelSys->SetLineStyle(7);
  graphCombPi0InvYieldRelSys->Draw("l,x0,same,e1");

  TLegend* legendRelTotErr2 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
  legendRelTotErr2->AddEntry(graphCombPi0InvYieldRelTot,"tot","p");
  legendRelTotErr2->AddEntry(graphCombPi0InvYieldRelStat,"stat","l");
  legendRelTotErr2->AddEntry(graphCombPi0InvYieldRelSys,"sys","l");
  legendRelTotErr2->Draw();

  TLatex *labelRelTotErrEnergy            = new TLatex(0.75,0.89,collisionSystempPb8TeV.Data());
  SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel,4);
  labelRelTotErrEnergy->SetTextFont(43);
  labelRelTotErrEnergy->Draw();
  TLatex *labelRelTotErrPi0               = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
  SetStyleTLatex( labelRelTotErrPi0, 0.85*textSizeLabelsPixel,4);
  labelRelTotErrPi0->SetTextFont(43);
  labelRelTotErrPi0->Draw();

  canvasRelTotErr->SaveAs(Form("%s/Pi0_Reldecomp.%s",outputDir.Data(),suffix.Data()));

  //**********************************************************************************************************************
  //**********************************************************************************************************************
  //**********************************************************************************************************************
  // plot weights + unc. for eta
  //**********************************************************************************************************************
  //**********************************************************************************************************************
  //**********************************************************************************************************************

  // Reading weights from output file for plotting
  ifstream fileWeightsReadEta;
  fileWeightsReadEta.open(fileNameOutputWeightingEta,ios_base::in);
  cout << "reading" << fileNameOutputWeightingEta << endl;
  Double_t xValuesReadEta[50];
  Double_t weightsReadEta[13][50];
  Int_t availableMeasEta[13]                  = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  Int_t nPtBinsReadEta                        = 0;
  Int_t nMeasSetEta                           = numbersofmeasEta;
  while(!fileWeightsReadEta.eof() && nPtBinsReadEta < 50){
      TString garbage                         = "";
      if (nPtBinsReadEta == 0){
          fileWeightsReadEta >> garbage ;
          for (Int_t i = 0; i < nMeasSetEta; i++){
              fileWeightsReadEta >> availableMeasEta[i] ;
          }
          cout << "read following measurements: ";
          for (Int_t i = 0; i < 13; i++){
              cout << availableMeasEta[i] << "\t" ;
          }
          cout << endl;
      } else {
          fileWeightsReadEta >> xValuesReadEta[nPtBinsReadEta-1];
          for (Int_t i = 0; i < nMeasSetEta; i++){
              fileWeightsReadEta >> weightsReadEta[availableMeasEta[i]][nPtBinsReadEta-1] ;
          }
          cout << "read: "<<  nPtBinsReadEta << "\t"<< xValuesReadEta[nPtBinsReadEta-1] << "\t" ;
          for (Int_t i = 0; i < nMeasSetEta; i++){
              cout << weightsReadEta[availableMeasEta[i]][nPtBinsReadEta-1] << "\t";
          }
          cout << endl;
      }
      nPtBinsReadEta++;
  }
  nPtBinsReadEta = nPtBinsReadEta-2 ;
  fileWeightsReadEta.close();

  for (Int_t i = 0; i < nMeasSetEta; i++){
      graphWeightsEta[availableMeasEta[i]]  = new TGraph(nPtBinsReadEta,xValuesReadEta,weightsReadEta[availableMeasEta[i]]);
      Int_t bin = 0;
      for (Int_t n = 0; n< nPtBinsReadEta; n++){
          if (graphWeightsEta[availableMeasEta[i]]->GetY()[bin] == 0) graphWeightsEta[availableMeasEta[i]]->RemovePoint(bin);
          else bin++;
      }
  }

  //**********************************************************************************************************************
  //******************************************* Plotting weights for Eta  ************************************************
  //**********************************************************************************************************************

  TCanvas* canvasWeightsEta = new TCanvas("canvasWeightsEta","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasWeightsEta, 0.08, 0.02, 0.035, 0.09);
  canvasWeightsEta->SetLogx();

  TH2F * histo2DWeightsEta;
  histo2DWeightsEta                              = new TH2F("histo2DWeightsEta","histo2DWeightsEta",11000,minPtEta,maxPtEta,1000,-0.5,1.1);
  SetStyleHistoTH2ForGraphs(histo2DWeightsEta, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DWeightsEta->GetXaxis()->SetMoreLogLabels();
  histo2DWeightsEta->GetXaxis()->SetNoExponent();
  canvasWeightsEta->cd();
  histo2DWeightsEta->Draw("copy");

      TLegend* legendAccWeightsEta               = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEta*1.35), 32);
      for (Int_t i = 0; i < nMeasSetEta; i++){
          DrawGammaSetMarkerTGraph(graphWeightsEta[availableMeasEta[i]],
                                  markerStyleDet[availableMeasEta[i]],
                                  markerSizeDet[availableMeasEta[i]]*0.5,
                                  colorDet[availableMeasEta[i]] ,
                                  colorDet[availableMeasEta[i]]);
          graphWeightsEta[availableMeasEta[i]]->Draw("p,same,e1");
          legendAccWeightsEta->AddEntry(graphWeightsEta[availableMeasEta[i]],nameMeasGlobalPlot[availableMeasEta[i]],"p");
      }
      legendAccWeightsEta->Draw();
      TLatex *labelWeightsEnergyEta              = new TLatex(0.7,0.20,collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelWeightsEnergyEta, 0.85*textSizeLabelsPixel,4);
      labelWeightsEnergyEta->SetTextFont(43);
      labelWeightsEnergyEta->Draw();
      TLatex *labelWeightsEta                 = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
      SetStyleTLatex( labelWeightsEta, 0.85*textSizeLabelsPixel,4);
      labelWeightsEta->SetTextFont(43);
      labelWeightsEta->Draw();

      DrawGammaLines(minPtEta, maxPtEta , 0.5, 0.5,0.1, kGray, 7);
      DrawGammaLines(minPtEta, maxPtEta , 0.4, 0.4,0.1, kGray, 1);
      DrawGammaLines(minPtEta, maxPtEta , 0.3, 0.3,0.1, kGray, 7);
      DrawGammaLines(minPtEta, maxPtEta , 0.2, 0.2,0.1, kGray, 3);

  canvasWeightsEta->SaveAs(Form("%s/Eta_WeightsMethods.%s",outputDir.Data(),suffix.Data()));


  //*********************************************************************************************************************
  //************************************ Visualize relative errors ******************************************************
  //*********************************************************************************************************************

  TCanvas* canvasRelSysErrEta                    = new TCanvas("canvasRelSysErrEta","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRelSysErrEta, 0.08, 0.02, 0.035, 0.09);
  canvasRelSysErrEta->SetLogx();

  TH2F * histo2DRelSysErrEta = new TH2F("histo2DRelSysErrEta","histo2DRelSysErrEta",11000,minPtEta,maxPtEta,1000,0,80.5);
  SetStyleHistoTH2ForGraphs(histo2DRelSysErrEta, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
  histo2DRelSysErr->GetXaxis()->SetNoExponent();
  histo2DRelSysErr->GetYaxis()->SetRangeUser(0,55.5);
  histo2DRelSysErr->Draw("copy");

      TLegend* legendRelSysErrEta                = GetAndSetLegend2(0.62, 0.94-(0.035*nMeasSetEta*1.35), 0.95, 0.94, 32);
      for (Int_t i = 0; i < nMeasSetEta; i++){
          DrawGammaSetMarkerTGraph(sysErrorRelCollectionEta[availableMeasEta[i]], markerStyleDet[availableMeasEta[i]], markerSizeDet[availableMeasEta[i]]*0.5, colorDet[availableMeasEta[i]],
                                  colorDet[availableMeasEta[i]]);
          sysErrorRelCollectionEta[availableMeasEta[i]]->Draw("p,same,e1");
          legendRelSysErrEta->AddEntry(sysErrorRelCollectionEta[availableMeasEta[i]],nameMeasGlobalPlot[availableMeasEta[i]],"p");
      }
      legendRelSysErrEta->Draw();

      TLatex *labelRelSysErrEnergyEta            = new TLatex(0.15,0.89,collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelRelSysErrEnergyEta, 0.85*textSizeLabelsPixel,4);
      labelRelSysErrEnergyEta->SetTextFont(43);
      labelRelSysErrEnergyEta->Draw();
      TLatex *labelRelSysErrEta               = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelSysErrEta, 0.85*textSizeLabelsPixel,4);
      labelRelSysErrEta->SetTextFont(43);
      labelRelSysErrEta->Draw();

  canvasRelSysErrEta->SaveAs(Form("%s/Eta_RelSysErr.%s",outputDir.Data(),suffix.Data()));

  //*********************************************************************************************************************
  //************************************ Visualize relative errors ******************************************************
  //*********************************************************************************************************************

  TCanvas* canvasRelStatErrEta                   = new TCanvas("canvasRelStatErrEta","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRelStatErrEta, 0.08, 0.02, 0.035, 0.09);
  canvasRelStatErrEta->SetLogx();

  TH2F * histo2DRelStatErrEta = new TH2F("histo2DRelStatErrEta","histo2DRelStatErrEta",11000,minPtEta,maxPtEta,1000,0,80.5);
  SetStyleHistoTH2ForGraphs(histo2DRelStatErrEta, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DRelStatErrEta->GetXaxis()->SetMoreLogLabels();
  histo2DRelStatErrEta->GetXaxis()->SetNoExponent();
  histo2DRelStatErrEta->GetYaxis()->SetRangeUser(0,55.5);
  histo2DRelStatErrEta->Draw("copy");

      TLegend* legendRelStatErrEta               = GetAndSetLegend2(0.14, 0.94-(0.035*nMeasSetEta*1.35), 0.45, 0.94, 32);
      for (Int_t i = 0; i < nMeasSetEta; i++){
          if (availableMeasEta[i]== 2){
              DrawGammaSetMarker(statErrorRelCollectionEta[availableMeasEta[i]], markerStyleDet[availableMeasEta[i]], markerSizeDet[availableMeasEta[i]]*0.5, colorDet[availableMeasEta[i]] ,
                          colorDet[availableMeasEta[i]]);
              TGraphAsymmErrors* graphDummy   = new TGraphAsymmErrors(statErrorRelCollectionEta[availableMeasEta[i]]);
              DrawGammaSetMarkerTGraphAsym(graphDummy, markerStyleDet[availableMeasEta[i]], markerSizeDet[availableMeasEta[i]]*0.5, colorDet[availableMeasEta[i]],
                                  colorDet[availableMeasEta[i]]);
              graphDummy->Draw("same,p,x0");
              legendRelStatErrEta->AddEntry(graphDummy,nameMeasGlobalPlot[availableMeasEta[i]],"p");

          } else {
              DrawGammaSetMarker(statErrorRelCollectionEta[availableMeasEta[i]], markerStyleDet[availableMeasEta[i]], markerSizeDet[availableMeasEta[i]]*0.5, colorDet[availableMeasEta[i]] ,
                          colorDet[availableMeasEta[i]]);
              statErrorRelCollectionEta[availableMeasEta[i]]->Draw("p,same,e1");
              legendRelStatErrEta->AddEntry(statErrorRelCollectionEta[availableMeasEta[i]],nameMeasGlobalPlot[availableMeasEta[i]],"p");

          }
      }
      legendRelStatErrEta->Draw();

      TLatex *labelRelStatErrEnergyEta           = new TLatex(0.75,0.89,collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelRelStatErrEnergyEta, 0.85*textSizeLabelsPixel,4);
      labelRelStatErrEnergyEta->SetTextFont(43);
      labelRelStatErrEnergyEta->Draw();
      TLatex *labelRelStatErrEta              = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelStatErrEta, 0.85*textSizeLabelsPixel,4);
      labelRelStatErrEta->SetTextFont(43);
      labelRelStatErrEta->Draw();

  canvasRelStatErrEta->SaveAs(Form("%s/Eta_RelStatErr.%s",outputDir.Data(),suffix.Data()));


  //************************************************************************************************************************
  //************************************************************************************************************************

  TGraphAsymmErrors* graphCombEtaInvYieldRelStat     = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldStat, "graphCombEtaInvYieldRelStat");
  TGraphAsymmErrors* graphCombEtaInvYieldRelSys      = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldSys, "graphCombEtaInvYieldRelSys");
  TGraphAsymmErrors* graphCombEtaInvYieldRelTot      = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldTot, "graphCombEtaInvYieldRelTot");

  TCanvas* canvasRelTotErrEta                    = new TCanvas("canvasRelTotErrEta","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRelTotErrEta, 0.08, 0.02, 0.035, 0.09);
  canvasRelTotErrEta->SetLogx();

  TH2F * histo2DRelTotErrEta = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,minPtEta,maxPtEta,1000,0,57.5);
  SetStyleHistoTH2ForGraphs(histo2DRelTotErrEta, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DRelTotErrEta->GetXaxis()->SetMoreLogLabels();
  histo2DRelTotErrEta->GetXaxis()->SetNoExponent();

  histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,57.5);
  histo2DRelTotErrEta->GetYaxis()->SetTitle("Err (%)");
  histo2DRelTotErrEta->Draw("copy");

  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombEtaInvYieldRelTot->Draw("p,same,e1");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
  graphCombEtaInvYieldRelStat->Draw("l,x0,same,e1");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
  graphCombEtaInvYieldRelSys->SetLineStyle(7);
  graphCombEtaInvYieldRelSys->Draw("l,x0,same,e1");

  TLegend* legendRelTotErr2Eta = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
  legendRelTotErr2Eta->AddEntry(graphCombEtaInvYieldRelTot,"tot","p");
  legendRelTotErr2Eta->AddEntry(graphCombEtaInvYieldRelStat,"stat","l");
  legendRelTotErr2Eta->AddEntry(graphCombEtaInvYieldRelSys,"sys","l");
  legendRelTotErr2Eta->Draw();

  TLatex *labelRelTotErrEnergyEta            = new TLatex(0.75,0.89,collisionSystempPb8TeV.Data());
  SetStyleTLatex( labelRelTotErrEnergyEta, 0.85*textSizeLabelsPixel,4);
  labelRelTotErrEnergyEta->SetTextFont(43);
  labelRelTotErrEnergyEta->Draw();
  TLatex *labelRelTotErrEta               = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
  SetStyleTLatex( labelRelTotErrEta, 0.85*textSizeLabelsPixel,4);
  labelRelTotErrEta->SetTextFont(43);
  labelRelTotErrEta->Draw();

  canvasRelTotErrEta->SaveAs(Form("%s/Eta_Reldecomp.%s",outputDir.Data(),suffix.Data()));

  //**********************************************************************************************************************
  //**********************************************************************************************************************
  //**********************************************************************************************************************
  // plot weights + unc. for eta/pi0
  //**********************************************************************************************************************
  //**********************************************************************************************************************
  //**********************************************************************************************************************

  // Reading weights from output file for plotting
  ifstream fileWeightsReadEtaToPi0;
  fileWeightsReadEtaToPi0.open(fileNameOutputWeightingEtaToPi0,ios_base::in);
  cout << "reading" << fileNameOutputWeightingEtaToPi0 << endl;
  Double_t xValuesReadEtaToPi0[50];
  Double_t weightsReadEtaToPi0[13][50];
  Int_t availableMeasEtaToPi0[13]                  = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  Int_t nPtBinsReadEtaToPi0                        = 0;
  while(!fileWeightsReadEtaToPi0.eof() && nPtBinsReadEtaToPi0 < 50){
      TString garbage                         = "";
      if (nPtBinsReadEtaToPi0 == 0){
          fileWeightsReadEtaToPi0 >> garbage ;
          for (Int_t i = 0; i < nMeasSetEta; i++){
              fileWeightsReadEtaToPi0 >> availableMeasEtaToPi0[i] ;
          }
          cout << "read following measurements: ";
          for (Int_t i = 0; i < 13; i++){
              cout << availableMeasEtaToPi0[i] << "\t" ;
          }
          cout << endl;
      } else {
          fileWeightsReadEtaToPi0 >> xValuesReadEtaToPi0[nPtBinsReadEtaToPi0-1];
          for (Int_t i = 0; i < nMeasSetEta; i++){
              fileWeightsReadEtaToPi0 >> weightsReadEtaToPi0[availableMeasEtaToPi0[i]][nPtBinsReadEtaToPi0-1] ;
          }
          cout << "read: "<<  nPtBinsReadEtaToPi0 << "\t"<< xValuesReadEtaToPi0[nPtBinsReadEtaToPi0-1] << "\t" ;
          for (Int_t i = 0; i < nMeasSetEta; i++){
              cout << weightsReadEtaToPi0[availableMeasEtaToPi0[i]][nPtBinsReadEtaToPi0-1] << "\t";
          }
          cout << endl;
      }
      nPtBinsReadEtaToPi0++;
  }
  nPtBinsReadEtaToPi0 = nPtBinsReadEtaToPi0-2 ;
  fileWeightsReadEtaToPi0.close();

  for (Int_t i = 0; i < nMeasSetEta; i++){
      graphWeightsEtaToPi0[availableMeasEtaToPi0[i]]  = new TGraph(nPtBinsReadEtaToPi0,xValuesReadEtaToPi0,weightsReadEtaToPi0[availableMeasEtaToPi0[i]]);
      Int_t bin = 0;
      for (Int_t n = 0; n< nPtBinsReadEtaToPi0; n++){
          if (graphWeightsEtaToPi0[availableMeasEtaToPi0[i]]->GetY()[bin] == 0) graphWeightsEtaToPi0[availableMeasEtaToPi0[i]]->RemovePoint(bin);
          else bin++;
      }
  }

  //**********************************************************************************************************************
  //******************************************* Plotting weights for EtaToPi0  ************************************************
  //**********************************************************************************************************************

  TCanvas* canvasWeightsEtaToPi0 = new TCanvas("canvasWeightsEtaToPi0","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasWeightsEtaToPi0, 0.08, 0.02, 0.035, 0.09);
  canvasWeightsEtaToPi0->SetLogx();

  TH2F * histo2DWeightsEtaToPi0;
  histo2DWeightsEtaToPi0                              = new TH2F("histo2DWeightsEtaToPi0","histo2DWeightsEtaToPi0",11000,minPtEta,maxPtEta,1000,-0.5,1.1);
  SetStyleHistoTH2ForGraphs(histo2DWeightsEtaToPi0, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DWeightsEtaToPi0->GetXaxis()->SetMoreLogLabels();
  histo2DWeightsEtaToPi0->GetXaxis()->SetNoExponent();
  canvasWeightsEtaToPi0->cd();
  histo2DWeightsEtaToPi0->Draw("copy");

      TLegend* legendAccWeightsEtaToPi0               = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEta*1.35), 32);
      for (Int_t i = 0; i < nMeasSetEta; i++){
          DrawGammaSetMarkerTGraph(graphWeightsEtaToPi0[availableMeasEtaToPi0[i]],
                                  markerStyleDet[availableMeasEtaToPi0[i]],
                                  markerSizeDet[availableMeasEtaToPi0[i]]*0.5,
                                  colorDet[availableMeasEtaToPi0[i]] ,
                                  colorDet[availableMeasEtaToPi0[i]]);
          graphWeightsEtaToPi0[availableMeasEtaToPi0[i]]->Draw("p,same,e1");
          legendAccWeightsEtaToPi0->AddEntry(graphWeightsEtaToPi0[availableMeasEtaToPi0[i]],nameMeasGlobalPlot[availableMeasEtaToPi0[i]],"p");
      }
      legendAccWeightsEtaToPi0->Draw();
      TLatex *labelWeightsEnergyEtaToPi0              = new TLatex(0.7,0.20,collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelWeightsEnergyEtaToPi0, 0.85*textSizeLabelsPixel,4);
      labelWeightsEnergyEtaToPi0->SetTextFont(43);
      labelWeightsEnergyEtaToPi0->Draw();
      TLatex *labelWeightsEtaToPi0                 = new TLatex(0.7,0.16,"#eta/#pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelWeightsEtaToPi0, 0.85*textSizeLabelsPixel,4);
      labelWeightsEtaToPi0->SetTextFont(43);
      labelWeightsEtaToPi0->Draw();

      DrawGammaLines(minPtEta, maxPtEta , 0.5, 0.5,0.1, kGray, 7);
      DrawGammaLines(minPtEta, maxPtEta , 0.4, 0.4,0.1, kGray, 1);
      DrawGammaLines(minPtEta, maxPtEta , 0.3, 0.3,0.1, kGray, 7);
      DrawGammaLines(minPtEta, maxPtEta , 0.2, 0.2,0.1, kGray, 3);

  canvasWeightsEtaToPi0->SaveAs(Form("%s/EtaToPi0_WeightsMethods.%s",outputDir.Data(),suffix.Data()));


  //********************************************************************************************************************
  //************************************ Visualize relative errors ******************************************************
  //*********************************************************************************************************************

  TCanvas* canvasRelSysErrEtaToPi0                    = new TCanvas("canvasRelSysErrEtaToPi0","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRelSysErrEtaToPi0, 0.08, 0.02, 0.035, 0.09);
  canvasRelSysErrEtaToPi0->SetLogx();

  TH2F * histo2DRelSysErrEtaToPi0 = new TH2F("histo2DRelSysErrEtaToPi0","histo2DRelSysErrEtaToPi0",11000,minPtEta,maxPtEta,1000,0,80.5);
  SetStyleHistoTH2ForGraphs(histo2DRelSysErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
  histo2DRelSysErr->GetXaxis()->SetNoExponent();
  histo2DRelSysErr->GetYaxis()->SetRangeUser(0,45.5);
  histo2DRelSysErr->Draw("copy");

      TLegend* legendRelSysErrEtaToPi0                = GetAndSetLegend2(0.62, 0.94-(0.035*nMeasSetEta*1.35), 0.95, 0.94, 32);
      for (Int_t i = 0; i < nMeasSetEta; i++){
          DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]], markerStyleDet[availableMeasEtaToPi0[i]], markerSizeDet[availableMeasEtaToPi0[i]]*0.5, colorDet[availableMeasEtaToPi0[i]],
                                  colorDet[availableMeasEtaToPi0[i]]);
          sysErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]]->Draw("p,same,e1");
          legendRelSysErrEtaToPi0->AddEntry(sysErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]],nameMeasGlobalPlot[availableMeasEtaToPi0[i]],"p");
      }
      legendRelSysErrEtaToPi0->Draw();

      TLatex *labelRelSysErrEnergyEtaToPi0            = new TLatex(0.15,0.89,collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelRelSysErrEnergyEtaToPi0, 0.85*textSizeLabelsPixel,4);
      labelRelSysErrEnergyEtaToPi0->SetTextFont(43);
      labelRelSysErrEnergyEtaToPi0->Draw();
      TLatex *labelRelSysErrEtaToPi0               = new TLatex(0.15,0.85,"#eta/#pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelSysErrEtaToPi0, 0.85*textSizeLabelsPixel,4);
      labelRelSysErrEtaToPi0->SetTextFont(43);
      labelRelSysErrEtaToPi0->Draw();

  canvasRelSysErrEtaToPi0->SaveAs(Form("%s/EtaToPi0_RelSysErr.%s",outputDir.Data(),suffix.Data()));

  //*********************************************************************************************************************
  //************************************ Visualize relative errors ******************************************************
  //*********************************************************************************************************************

  TCanvas* canvasRelStatErrEtaToPi0                   = new TCanvas("canvasRelStatErrEtaToPi0","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRelStatErrEtaToPi0, 0.08, 0.02, 0.035, 0.09);
  canvasRelStatErrEtaToPi0->SetLogx();

  TH2F * histo2DRelStatErrEtaToPi0 = new TH2F("histo2DRelStatErrEtaToPi0","histo2DRelStatErrEtaToPi0",11000,minPtEta,maxPtEta,1000,0,80.5);
  SetStyleHistoTH2ForGraphs(histo2DRelStatErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DRelStatErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
  histo2DRelStatErrEtaToPi0->GetXaxis()->SetNoExponent();
  histo2DRelStatErrEtaToPi0->GetYaxis()->SetRangeUser(0,45.5);
  histo2DRelStatErrEtaToPi0->Draw("copy");

  TLegend* legendRelStatErrEtaToPi0               = GetAndSetLegend2(0.14, 0.94-(0.035*nMeasSetEta*1.35), 0.45, 0.94, 32);
  for (Int_t i = 0; i < nMeasSetEta; i++){
      if (availableMeasEtaToPi0[i]== 2){
          DrawGammaSetMarker(statErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]], markerStyleDet[availableMeasEtaToPi0[i]], markerSizeDet[availableMeasEtaToPi0[i]]*0.5, colorDet[availableMeasEtaToPi0[i]] ,
                      colorDet[availableMeasEtaToPi0[i]]);
          TGraphAsymmErrors* graphDummy   = new TGraphAsymmErrors(statErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]]);
          DrawGammaSetMarkerTGraphAsym(graphDummy, markerStyleDet[availableMeasEtaToPi0[i]], markerSizeDet[availableMeasEtaToPi0[i]]*0.5, colorDet[availableMeasEtaToPi0[i]],
                              colorDet[availableMeasEtaToPi0[i]]);
          graphDummy->Draw("same,p,x0");
          legendRelStatErrEtaToPi0->AddEntry(graphDummy,nameMeasGlobalPlot[availableMeasEtaToPi0[i]],"p");

      } else {
          DrawGammaSetMarker(statErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]], markerStyleDet[availableMeasEtaToPi0[i]], markerSizeDet[availableMeasEtaToPi0[i]]*0.5, colorDet[availableMeasEtaToPi0[i]] ,
                      colorDet[availableMeasEtaToPi0[i]]);
          statErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]]->Draw("p,same,e1");
          legendRelStatErrEtaToPi0->AddEntry(statErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]],nameMeasGlobalPlot[availableMeasEtaToPi0[i]],"p");

      }
  }
  legendRelStatErrEtaToPi0->Draw();

  TLatex *labelRelStatErrEnergyEtaToPi0           = new TLatex(0.75,0.89,collisionSystempPb8TeV.Data());
  SetStyleTLatex( labelRelStatErrEnergyEtaToPi0, 0.85*textSizeLabelsPixel,4);
  labelRelStatErrEnergyEtaToPi0->SetTextFont(43);
  labelRelStatErrEnergyEtaToPi0->Draw();
  TLatex *labelRelStatErrEtaToPi0              = new TLatex(0.75,0.85,"#eta/#pi^{0} #rightarrow #gamma#gamma");
  SetStyleTLatex( labelRelStatErrEtaToPi0, 0.85*textSizeLabelsPixel,4);
  labelRelStatErrEtaToPi0->SetTextFont(43);
  labelRelStatErrEtaToPi0->Draw();

  canvasRelStatErrEtaToPi0->SaveAs(Form("%s/EtaToPi0_RelStatErr.%s",outputDir.Data(),suffix.Data()));


  //************************************************************************************************************************
  //************************************************************************************************************************

  TGraphAsymmErrors* graphCombEtaToPi0RelStat     = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Stat, "graphCombEtaToPi0RelStat");
  TGraphAsymmErrors* graphCombEtaToPi0RelSys      = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Sys, "graphCombEtaToPi0RelSys");
  TGraphAsymmErrors* graphCombEtaToPi0RelTot      = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Tot, "graphCombEtaToPi0RelTot");

  TCanvas* canvasRelTotErrEtaToPi0                    = new TCanvas("canvasRelTotErrEtaToPi0","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRelTotErrEtaToPi0, 0.08, 0.02, 0.035, 0.09);
  canvasRelTotErrEtaToPi0->SetLogx();

  TH2F * histo2DRelTotErrEtaToPi0 = new TH2F("histo2DRelTotErrEtaToPi0","histo2DRelTotErrEtaToPi0",11000,minPtEta,maxPtEta,1000,0,57.5);
  SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DRelTotErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
  histo2DRelTotErrEtaToPi0->GetXaxis()->SetNoExponent();

  histo2DRelTotErrEtaToPi0->GetYaxis()->SetRangeUser(0,57.5);
  histo2DRelTotErrEtaToPi0->GetYaxis()->SetTitle("Err (%)");
  histo2DRelTotErrEtaToPi0->Draw("copy");

  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombEtaToPi0RelTot->Draw("p,same,e1");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
  graphCombEtaToPi0RelStat->Draw("l,x0,same,e1");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
  graphCombEtaToPi0RelSys->SetLineStyle(7);
  graphCombEtaToPi0RelSys->Draw("l,x0,same,e1");

  TLegend* legendRelTotErr2EtaToPi0 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
  legendRelTotErr2EtaToPi0->AddEntry(graphCombEtaToPi0RelTot,"tot","p");
  legendRelTotErr2EtaToPi0->AddEntry(graphCombEtaToPi0RelStat,"stat","l");
  legendRelTotErr2EtaToPi0->AddEntry(graphCombEtaToPi0RelSys,"sys","l");
  legendRelTotErr2EtaToPi0->Draw();

  TLatex *labelRelTotErrEnergyEtaToPi0            = new TLatex(0.75,0.89,collisionSystempPb8TeV.Data());
  SetStyleTLatex( labelRelTotErrEnergyEtaToPi0, 0.85*textSizeLabelsPixel,4);
  labelRelTotErrEnergyEtaToPi0->SetTextFont(43);
  labelRelTotErrEnergyEtaToPi0->Draw();
  TLatex *labelRelTotErrEtaToPi0               = new TLatex(0.75,0.85,"#eta/#pi^{0} #rightarrow #gamma#gamma");
  SetStyleTLatex( labelRelTotErrEtaToPi0, 0.85*textSizeLabelsPixel,4);
  labelRelTotErrEtaToPi0->SetTextFont(43);
  labelRelTotErrEtaToPi0->Draw();

  canvasRelTotErrEtaToPi0->SaveAs(Form("%s/EtaToPi0_Reldecomp.%s",outputDir.Data(),suffix.Data()));
  
  
  
  
  //**********************************************************************************************************************
  //************************************* Calculating bin shifted spectra & fitting **************************************
  //**********************************************************************************************************************

  // Cloning spectra
  TGraphAsymmErrors* graphCombPi0InvYieldTotUnShifted    = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone("Unshifted");
  TGraphAsymmErrors* graphCombPi0InvYieldStatUnShifted   = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone("UnshiftedStat");
  TGraphAsymmErrors* graphCombPi0InvYieldSysUnShifted    = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone("UnshiftedSys");

  TGraphAsymmErrors* graphPi0InvYieldStatUnShifted[13];
  TGraphAsymmErrors* graphPi0InvYieldSysUnShifted[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  graphPi0InvYieldStatUnShifted[i] = (TGraphAsymmErrors*)graphPi0InvYieldStat[i]->Clone(Form("UnshiftedStat%s",nameMeasGlobalPlot[i].Data()));
  graphPi0InvYieldSysUnShifted[i]  = (TGraphAsymmErrors*)graphPi0InvYieldSys[i] ->Clone(Form("UnshiftedSys%s",nameMeasGlobalPlot[i].Data()));
  }
  }

  // Calculating binshifts
  Double_t paramGraph[3]                      = {1.0e1, 8., 0.13};
  TF1* fitInvXSectionPi0              = FitObject("l","fitInvXSectionPi0","Pi0",graphCombPi0InvYieldTot,0.3,25.,paramGraph,"QNRMEX0+");

  if(bWCorrection.Contains("X")){
  TF1* fitTsallisPi0PtMult        = FitObject("tmpt","TsallisMultWithPtPi0pPb8TeV","Pi0");
  fitTsallisPi0PtMult->SetParameters(paramGraph[0],paramGraph[1], paramGraph[2]) ; // standard parameter optimize if necessary

  graphCombPi0InvYieldTot      = ApplyXshift(graphCombPi0InvYieldTot, fitTsallisPi0PtMult, "Pi0", kTRUE);

  cout << "comb" << endl;
  graphCombPi0InvYieldStat     = ApplyXshiftIndividualSpectra (graphCombPi0InvYieldTot,
                                                              graphCombPi0InvYieldStat,
                                                              fitTsallisPi0PtMult,
                                                              0, graphCombPi0InvYieldStat->GetN());
  graphCombPi0InvYieldSys      = ApplyXshiftIndividualSpectra (graphCombPi0InvYieldTot,
                                                              graphCombPi0InvYieldSys,
                                                              fitTsallisPi0PtMult,
                                                              0, graphCombPi0InvYieldSys->GetN());

  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  cout << nameMeasGlobalPlot[i].Data() << endl;
  graphPi0InvYieldStat[i]      = ApplyXshiftIndividualSpectra (graphCombPi0InvYieldTot,
                                                                  graphPi0InvYieldStat[i],
                                                                  fitTsallisPi0PtMult,
                                                                  offSetPi0Shifting[i], nComBinsPi0Shifting[i]);
  graphPi0InvYieldSys[i]       = ApplyXshiftIndividualSpectra (graphCombPi0InvYieldTot,
                                                                  graphPi0InvYieldSys[i],
                                                                  fitTsallisPi0PtMult,
                                                                  offSetPi0Shifting[i], nComBinsPi0Shifting[i]);
  }
  }

  TF1* fitTsallisPi0PtMultFromShift                 = FitObject("tmpt","TsallisMultWithPtPi0FromShift","Pi0");
  fitTsallisPi0PtMultFromShift->SetRange(0.3,25.);
  fitTsallisPi0PtMultFromShift->SetParameters(fitTsallisPi0PtMult->GetParameter(0),fitTsallisPi0PtMult->GetParameter(1), fitTsallisPi0PtMult->GetParameter(2));

  TF1* fitTsallisPi0PtMultFromShiftScaled = new TF1("TsallisMultWithPtPi0FromShiftScaled","(1/x)*TsallisMultWithPtPi0FromShift",0.3,25.);

  //***************************************************************************************************************
  //************************************Plotting binshift corrections *********************************************
  //***************************************************************************************************************

  TCanvas* canvasShift = new TCanvas("canvasShift","",0,0,1000,900);// gives the page size
  DrawGammaCanvasSettings( canvasShift, 0.10, 0.017, 0.015, 0.1);
  canvasShift->SetLogx(1);

  Size_t textSizeSpectra          = 0.04;
  TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,minPtPi0, maxPtPi0);
  SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (X) / no shift",
                      0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 1.1, 1.2);
  histoBinShift->GetXaxis()->SetMoreLogLabels(1);
  histoBinShift->GetXaxis()->SetNoExponent(kTRUE);
  histoBinShift->GetYaxis()->SetRangeUser(0.95,1.05);
  histoBinShift->DrawCopy();

  TGraphAsymmErrors* graphCombPi0InvYieldTotUnShifted_clone = (TGraphAsymmErrors*) graphCombPi0InvYieldTotUnShifted->Clone("graphCombPi0InvYieldTotUnShifted_clone");

  Int_t numberPoints   = graphCombPi0InvYieldTotUnShifted_clone->GetN();
  Double_t *xPoint     = graphCombPi0InvYieldTotUnShifted_clone->GetX();
  Double_t* xvalueErrUp  = graphCombPi0InvYieldTotUnShifted_clone->GetEXhigh();
  Double_t* xvalueErrLow = graphCombPi0InvYieldTotUnShifted_clone->GetEXlow();
  Double_t *xPointShift= graphCombPi0InvYieldTot->GetX();
  for (Int_t i=0; i<numberPoints; i++) {
  graphCombPi0InvYieldTotUnShifted_clone->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
  graphCombPi0InvYieldTotUnShifted_clone->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
  }
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldTotUnShifted_clone, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombPi0InvYieldTotUnShifted_clone->Draw("p same");

  TLatex *labelRatioToFitBinShift   = new TLatex(0.72, 0.91, collisionSystempPb8TeV.Data());
  SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel,4);
  labelRatioToFitBinShift->SetTextFont(43);
  labelRatioToFitBinShift->Draw();
  TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.852, 0.86, labelALICEPlot.Data());
  SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel,4);
  labelRatioToFitALICEBinShift->SetTextFont(43);
  labelRatioToFitALICEBinShift->Draw();
  TLatex *labelRatioToFitPi0BinShift      = new TLatex(0.826, 0.807, "#pi^{0} #rightarrow #gamma#gamma");
  SetStyleTLatex( labelRatioToFitPi0BinShift, textSizeLabelsPixel,4);
  labelRatioToFitPi0BinShift->SetTextFont(43);
  labelRatioToFitPi0BinShift->Draw();

  canvasShift->Update();
  canvasShift->SaveAs(Form("%s/BinShiftCorrection_Pi0.%s",outputDir.Data(),suffix.Data()));
  canvasShift->SetLogx(0);

  // *************************************************************************************************************
  // Plot control graphs
  // *************************************************************************************************************

  TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
  DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
  canvasDummy2->SetLogy();
  canvasDummy2->SetLogx();
  TH2F * histo2DDummy2;
  histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,minPtPi0, maxPtPi0,1000,minInvYieldPi0,maxInvYieldPi0);
  SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 0.8,1.55);
  histo2DDummy2->DrawCopy();


  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  DrawGammaSetMarkerTGraphAsym(graphPi0InvYieldStat[i], markerStyleDet[i] ,markerSizeDet[i]/2, colorDet[i], colorDet[i]);
  graphPi0InvYieldStat[i]->Draw("pEsame");
  DrawGammaSetMarkerTGraphAsym(graphPi0InvYieldSys[i], markerStyleDet[i] ,markerSizeDet[i]/2, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
  graphPi0InvYieldSys[i]->Draw("pEsame");
  }
  }

  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatUnShifted, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
  graphCombPi0InvYieldStatUnShifted->Draw("pEsame");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
  graphCombPi0InvYieldStat->Draw("pEsame");

  fitInvXSectionPi0->SetLineColor(kBlue+2);
  fitInvXSectionPi0->Draw("same");

  fitTsallisPi0PtMultFromShiftScaled->SetLineColor(kRed+2);
  fitTsallisPi0PtMultFromShiftScaled->Draw("same");

  canvasDummy2->Update();
  canvasDummy2->SaveAs(Form("%s/ComparisonShiftedPi0_pPb8TeV.%s",outputDir.Data(),suffix.Data()));
  delete canvasDummy2;
  }

  TGraphAsymmErrors* graphCombPi0InvYieldStat_WOXErr = (TGraphAsymmErrors*) graphCombPi0InvYieldStat->Clone("graphCombPi0InvYieldStatA_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphCombPi0InvYieldStat_WOXErr);

  TGraphAsymmErrors* graphPi0InvYieldStat_WOXErr[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  graphPi0InvYieldStat_WOXErr[i] = (TGraphAsymmErrors*) graphPi0InvYieldStat[i]->Clone(Form("graphPi0InvYieldStat_%i_WOXErr",i));
  ProduceGraphAsymmWithoutXErrors(graphPi0InvYieldStat_WOXErr[i]);
  }
  }


  //**********************************************************************************************************************
  //************************************* Calculating bin shifted spectra & fitting **************************************
  //**********************************************************************************************************************

  // Cloning spectra
  TGraphAsymmErrors* graphCombEtaInvYieldTotUnShifted    = (TGraphAsymmErrors*)graphCombEtaInvYieldTot->Clone("Unshifted_Eta");
  TGraphAsymmErrors* graphCombEtaInvYieldStatUnShifted   = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone("UnshiftedStat_Eta");
  TGraphAsymmErrors* graphCombEtaInvYieldSysUnShifted    = (TGraphAsymmErrors*)graphCombEtaInvYieldSys->Clone("UnshiftedSys_Eta");

  TGraphAsymmErrors* graphEtaInvYieldStatUnShifted[13];
  TGraphAsymmErrors* graphEtaInvYieldSysUnShifted[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  graphEtaInvYieldStatUnShifted[i] = (TGraphAsymmErrors*)graphEtaInvYieldStat[i]->Clone(Form("UnshiftedStatEta%s",nameMeasGlobalPlot[i].Data()));
  graphEtaInvYieldSysUnShifted[i]  = (TGraphAsymmErrors*)graphEtaInvYieldSys[i] ->Clone(Form("UnshiftedSysEta%s",nameMeasGlobalPlot[i].Data()));
  }
  }

  // Calculating binshifts
  Double_t paramGraphEta[3]           = {1.0e1, 8., 0.13};
  TF1* fitInvXSectionEta              = FitObject("l","fitInvXSectionEta","Eta",graphCombEtaInvYieldTot,0.4,18.,paramGraphEta,"QNRMEX0+");

  if(bWCorrection.Contains("X")){
  TF1* fitTsallisEtaPtMult        = FitObject("tmpt","TsallisMultWithPtEtapPb8TeV","Eta");
  fitTsallisEtaPtMult->SetParameters(paramGraphEta[0],paramGraphEta[1], paramGraphEta[2]) ; // standard parameter optimize if necessary

  graphCombEtaInvYieldTot      = ApplyXshift(graphCombEtaInvYieldTot, fitTsallisEtaPtMult, "Eta", kTRUE);

  cout << "comb" << endl;
  graphCombEtaInvYieldStat     = ApplyXshiftIndividualSpectra (graphCombEtaInvYieldTot,
                                                              graphCombEtaInvYieldStat,
                                                              fitTsallisEtaPtMult,
                                                              0, graphCombEtaInvYieldStat->GetN());
  graphCombEtaInvYieldSys      = ApplyXshiftIndividualSpectra (graphCombEtaInvYieldTot,
                                                              graphCombEtaInvYieldSys,
                                                              fitTsallisEtaPtMult,
                                                              0, graphCombEtaInvYieldSys->GetN());

  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  cout << nameMeasGlobalPlot[i].Data() << endl;
  graphEtaInvYieldStat[i]      = ApplyXshiftIndividualSpectra (graphCombEtaInvYieldTot,
                                                                  graphEtaInvYieldStat[i],
                                                                  fitTsallisEtaPtMult,
                                                                  offSetEtaShifting[i], nComBinsEtaShifting[i]);
  graphEtaInvYieldSys[i]       = ApplyXshiftIndividualSpectra (graphCombEtaInvYieldTot,
                                                                  graphEtaInvYieldSys[i],
                                                                  fitTsallisEtaPtMult,
                                                                  offSetEtaShifting[i], nComBinsEtaShifting[i]);
  }
  }

  TF1* fitTsallisEtaPtMultFromShift                 = FitObject("tmpt","TsallisMultWithPtEtaFromShift","Eta");
  fitTsallisEtaPtMultFromShift->SetRange(0.4,18.);
  fitTsallisEtaPtMultFromShift->SetParameters(fitTsallisEtaPtMult->GetParameter(0),fitTsallisEtaPtMult->GetParameter(1), fitTsallisEtaPtMult->GetParameter(2));

  TF1* fitTsallisEtaPtMultFromShiftScaled = new TF1("TsallisMultWithPtEtaFromShiftScaled","(1/x)*TsallisMultWithPtEtaFromShift",0.4,18.);

  //***************************************************************************************************************
  //************************************Plotting binshift corrections *********************************************
  //***************************************************************************************************************

  TCanvas* canvasShift = new TCanvas("canvasShift","",0,0,1000,900);// gives the page size
  DrawGammaCanvasSettings( canvasShift, 0.10, 0.017, 0.015, 0.1);
  canvasShift->SetLogx(1);

  Size_t textSizeSpectra          = 0.04;
  TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,minPtEta, maxPtEta);
  SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (X) / no shift",
                      0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 1.1, 1.2);
  histoBinShift->GetXaxis()->SetMoreLogLabels(1);
  histoBinShift->GetXaxis()->SetNoExponent(kTRUE);
  histoBinShift->GetYaxis()->SetRangeUser(0.95,1.05);
  histoBinShift->DrawCopy();

  TGraphAsymmErrors* graphCombEtaInvYieldTotUnShifted_clone = (TGraphAsymmErrors*) graphCombEtaInvYieldTotUnShifted->Clone("graphCombEtaInvYieldTotUnShifted_clone");

  Int_t numberPoints   = graphCombEtaInvYieldTotUnShifted_clone->GetN();
  Double_t *xPoint     = graphCombEtaInvYieldTotUnShifted_clone->GetX();
  Double_t* xvalueErrUp  = graphCombEtaInvYieldTotUnShifted_clone->GetEXhigh();
  Double_t* xvalueErrLow = graphCombEtaInvYieldTotUnShifted_clone->GetEXlow();
  Double_t *xPointShift= graphCombEtaInvYieldTot->GetX();
  for (Int_t i=0; i<numberPoints; i++) {
  graphCombEtaInvYieldTotUnShifted_clone->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
  graphCombEtaInvYieldTotUnShifted_clone->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
  }
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldTotUnShifted_clone, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombEtaInvYieldTotUnShifted_clone->Draw("p same");

  TLatex *labelRatioToFitBinShift   = new TLatex(0.72, 0.91, collisionSystempPb8TeV.Data());
  SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel,4);
  labelRatioToFitBinShift->SetTextFont(43);
  labelRatioToFitBinShift->Draw();
  TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.852, 0.86, labelALICEPlot.Data());
  SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel,4);
  labelRatioToFitALICEBinShift->SetTextFont(43);
  labelRatioToFitALICEBinShift->Draw();
  TLatex *labelRatioToFitEtaBinShift      = new TLatex(0.826, 0.807, "#eta #rightarrow #gamma#gamma");
  SetStyleTLatex( labelRatioToFitEtaBinShift, textSizeLabelsPixel,4);
  labelRatioToFitEtaBinShift->SetTextFont(43);
  labelRatioToFitEtaBinShift->Draw();

  canvasShift->Update();
  canvasShift->SaveAs(Form("%s/BinShiftCorrection_Eta.%s",outputDir.Data(),suffix.Data()));
  canvasShift->SetLogx(0);

  // *************************************************************************************************************
  // Plot control graphs
  // *************************************************************************************************************

  TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
  DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
  canvasDummy2->SetLogy();
  canvasDummy2->SetLogx();
  TH2F * histo2DDummy2;
  histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,minPtEta,maxPtEta,1000,minInvYieldEta,maxInvYieldEta);
  SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 0.8,1.55);
  histo2DDummy2->DrawCopy();


  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  DrawGammaSetMarkerTGraphAsym(graphEtaInvYieldStat[i], markerStyleDet[i] ,markerSizeDet[i]/2, colorDet[i], colorDet[i]);
  graphEtaInvYieldStat[i]->Draw("pEsame");
  DrawGammaSetMarkerTGraphAsym(graphEtaInvYieldSys[i], markerStyleDet[i] ,markerSizeDet[i]/2, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
  graphEtaInvYieldSys[i]->Draw("pEsame");
  }
  }

  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatUnShifted, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
  graphCombEtaInvYieldStatUnShifted->Draw("pEsame");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
  graphCombEtaInvYieldStat->Draw("pEsame");

  fitInvXSectionEta->SetLineColor(kBlue+2);
  fitInvXSectionEta->Draw("same");

  fitTsallisEtaPtMultFromShiftScaled->SetLineColor(kRed+2);
  fitTsallisEtaPtMultFromShiftScaled->Draw("same");

  canvasDummy2->Update();
  canvasDummy2->SaveAs(Form("%s/ComparisonShiftedEta_pPb8TeV.%s",outputDir.Data(),suffix.Data()));
  delete canvasDummy2;
  }

  TGraphAsymmErrors* graphCombEtaInvYieldStat_WOXErr = (TGraphAsymmErrors*) graphCombEtaInvYieldStat->Clone("graphCombEtaInvYieldStatA_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphCombEtaInvYieldStat_WOXErr);

  TGraphAsymmErrors* graphEtaInvYieldStat_WOXErr[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  graphEtaInvYieldStat_WOXErr[i] = (TGraphAsymmErrors*) graphEtaInvYieldStat[i]->Clone(Form("graphEtaInvYieldStat_%i_WOXErr",i));
  ProduceGraphAsymmWithoutXErrors(graphEtaInvYieldStat_WOXErr[i]);
  }
  }

  // *************************************************************************************************************
  // Shift graphs in Y direction as well if desired
  // *************************************************************************************************************
  TGraphAsymmErrors* graphCombPi0InvYieldTot_yShifted         = NULL;
  TGraphAsymmErrors* graphCombPi0InvYieldStat_yShifted        = NULL;
  TGraphAsymmErrors* graphCombPi0InvYieldSys_yShifted         = NULL;

  TGraphAsymmErrors* graphPi0InvYieldStat_yShifted[13];
  TGraphAsymmErrors* graphPi0InvYieldSys_yShifted[13];

  if(bWCorrection.Contains("Y") ){
  graphCombPi0InvYieldTot_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvYieldTotUnShifted->Clone("Pi0YShiftedCombTot");
  graphCombPi0InvYieldTot_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldTot_yShifted, fitInvXSectionPi0);
  graphCombPi0InvYieldStat_yShifted       = (TGraphAsymmErrors*)graphCombPi0InvYieldStatUnShifted->Clone("Pi0YShiftedCombStat");
  graphCombPi0InvYieldStat_yShifted       =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldStat_yShifted, fitInvXSectionPi0);
  graphCombPi0InvYieldSys_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvYieldSysUnShifted->Clone("Pi0YShiftedCombSys");
  graphCombPi0InvYieldSys_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldSys_yShifted, fitInvXSectionPi0);

  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  graphPi0InvYieldStat_yShifted[i]       = (TGraphAsymmErrors*)graphPi0InvYieldStatUnShifted[i]->Clone("Pi0YShiftedCombStat");
  graphPi0InvYieldStat_yShifted[i]       =  ApplyYshiftIndividualSpectra( graphPi0InvYieldStat_yShifted[i], fitInvXSectionPi0);
  graphPi0InvYieldSys_yShifted[i]        = (TGraphAsymmErrors*)graphPi0InvYieldSysUnShifted[i]->Clone("Pi0YShiftedCombSys");
  graphPi0InvYieldSys_yShifted[i]        =  ApplyYshiftIndividualSpectra( graphPi0InvYieldSys_yShifted[i], fitInvXSectionPi0);
  }
  }
  }

  // *************************************************************************************************************
  // Shift spectra in Y  direction as well if desired
  // *************************************************************************************************************

  TGraphAsymmErrors* graphCombEtaInvYieldTot_yShifted         = NULL;
  TGraphAsymmErrors* graphCombEtaInvYieldStat_yShifted        = NULL;
  TGraphAsymmErrors* graphCombEtaInvYieldSys_yShifted         = NULL;

  TGraphAsymmErrors* graphEtaInvYieldStat_yShifted[13];
  TGraphAsymmErrors* graphEtaInvYieldSys_yShifted[13];

  if(bWCorrection.Contains("Y") ){
  graphCombEtaInvYieldTot_yShifted        = (TGraphAsymmErrors*)graphCombEtaInvYieldTotUnShifted->Clone("EtaYShiftedCombTot");
  graphCombEtaInvYieldTot_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldTot_yShifted, fitInvXSectionEta);
  graphCombEtaInvYieldStat_yShifted       = (TGraphAsymmErrors*)graphCombEtaInvYieldStatUnShifted->Clone("EtaYShiftedCombStat");
  graphCombEtaInvYieldStat_yShifted       =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldStat_yShifted, fitInvXSectionEta);
  graphCombEtaInvYieldSys_yShifted        = (TGraphAsymmErrors*)graphCombEtaInvYieldSysUnShifted->Clone("EtaYShiftedCombSys");
  graphCombEtaInvYieldSys_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldSys_yShifted, fitInvXSectionEta);

  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  graphEtaInvYieldStat_yShifted[i]       = (TGraphAsymmErrors*)graphEtaInvYieldStatUnShifted[i]->Clone("EtaYShiftedCombStat");
  graphEtaInvYieldStat_yShifted[i]       =  ApplyYshiftIndividualSpectra( graphEtaInvYieldStat_yShifted[i], fitInvXSectionEta);
  graphEtaInvYieldSys_yShifted[i]        = (TGraphAsymmErrors*)graphEtaInvYieldSysUnShifted[i]->Clone("EtaYShiftedCombSys");
  graphEtaInvYieldSys_yShifted[i]        =  ApplyYshiftIndividualSpectra( graphEtaInvYieldSys_yShifted[i], fitInvXSectionEta);
  }
  }
  }



   // *******************************************************************************************************
  // *************************** RpPb calculation for pi0 for individual measurements **********************
  // *******************************************************************************************************
  TGraphAsymmErrors* graphRpPbIndStatPi0[11];
  TH1D* histRpPbIndStatPi0[11];
  TGraphAsymmErrors* graphRpPbIndSystPi0[11];
  TGraphAsymmErrors* graphRpPbIndCombPi0[11];
  TGraph* graphWeightsPi0RpPb[11];
  TGraphAsymmErrors* statErrorRelCollectionPi0RpA[11];
  TGraphAsymmErrors* sysErrorRelCollectionPi0RpA[11];


  for (Int_t i = 0; i < 11; i++){
    graphWeightsPi0RpPb[i]                              = NULL;
    graphRpPbIndStatPi0[i]                              = NULL;
    histRpPbIndStatPi0[i]                               = NULL;
    graphRpPbIndSystPi0[i]                              = NULL;
    graphRpPbIndCombPi0[i]                              = NULL;
    statErrorRelCollectionPi0RpA[i]                     = NULL;
    sysErrorRelCollectionPi0RpA[i]                      = NULL;
    if (haveRefPPPi0[i] &&  directoryPi0[i] && havePi0NuclModFac[i]){
      histRpPbIndStatPi0[i]                       = (TH1D*)histoPi0NuclModFac[i]->Clone(Form("histRpPbInd%sStatPi0",nameMeasGlobalPlot[i].Data()));
      graphRpPbIndStatPi0[i]                       = (TGraphAsymmErrors*)graphPi0NuclModFacStat[i]->Clone(Form("graphRpPbInd%sStatPi0",nameMeasGlobalPlot[i].Data()));
      graphRpPbIndSystPi0[i]                       = (TGraphAsymmErrors*)graphPi0NuclModFacSys[i]->Clone(Form("graphRpPbInd%sSystPi0",nameMeasGlobalPlot[i].Data()));

      if (histRpPbIndStatPi0[i]){
        statErrorRelCollectionPi0RpA[i]     = new TGraphAsymmErrors(histRpPbIndStatPi0[i]);
        while (statErrorRelCollectionPi0RpA[i]->GetY()[0] == 0 && statErrorRelCollectionPi0RpA[i]->GetN() > 0 )
          statErrorRelCollectionPi0RpA[i]->RemovePoint(0);
        while (statErrorRelCollectionPi0RpA[i]->GetY()[statErrorRelCollectionPi0RpA[i]->GetN()-1] == 0 && statErrorRelCollectionPi0RpA[i]->GetN() > 0 )
          statErrorRelCollectionPi0RpA[i]->RemovePoint(statErrorRelCollectionPi0RpA[i]->GetN()-1);
        if (statErrorRelCollectionPi0RpA[i]->GetN() > 0 )
          statErrorRelCollectionPi0RpA[i] = CalculateRelErrUpAsymmGraph( statErrorRelCollectionPi0RpA[i], Form("relativeStatErrorPi0RpPb_%s", nameMeasGlobal[i].Data()));
        else
          statErrorRelCollectionPi0RpA[i] = NULL;
      }
      if (graphRpPbIndSystPi0[i]){
        sysErrorRelCollectionPi0RpA[i]      = (TGraphAsymmErrors*)graphRpPbIndSystPi0[i]->Clone(Form("relativeSysErrorPi0RpPb_%s", nameMeasGlobal[i].Data()));
        while (sysErrorRelCollectionPi0RpA[i]->GetY()[0]== 0 && sysErrorRelCollectionPi0RpA[i]->GetN() > 0 )
          sysErrorRelCollectionPi0RpA[i]->RemovePoint(0);
        while (sysErrorRelCollectionPi0RpA[i]->GetY()[sysErrorRelCollectionPi0RpA[i]->GetN()-1] == 0 && sysErrorRelCollectionPi0RpA[i]->GetN() > 0 )
          sysErrorRelCollectionPi0RpA[i]->RemovePoint(sysErrorRelCollectionPi0RpA[i]->GetN()-1);
        if (sysErrorRelCollectionPi0RpA[i]->GetN() > 0 )
          sysErrorRelCollectionPi0RpA[i]  = CalculateRelErrUpAsymmGraph( graphRpPbIndSystPi0[i], Form("relativeSysErrorPi0RpPb_%s", nameMeasGlobal[i].Data()));
        else
          sysErrorRelCollectionPi0RpA[i]  = NULL;

      }
    }
  }


  // *******************************************************************************************************
  // *************************** RpPb calculation for pi0 **************************************************
  // *******************************************************************************************************
  TGraphAsymmErrors* graphRpPbCombStatPi0         = NULL;
  TGraphAsymmErrors* graphRpPbCombSystPi0         = NULL;
  TGraphAsymmErrors* graphRpPbCombCombPi0         = NULL;

  if(bWCorrection.Contains("Y") ){
    TString fileNamePi0RpPbOutputWeighting       = Form("%s/Pi0RpPb_WeightingMethod.dat",outputDir.Data());
    graphRpPbCombCombPi0                         = CombinePtPointsSpectraFullCorrMat(   histRpPbIndStatPi0,    graphRpPbIndSystPi0,
                                                                                        xPtLimits, nBinsPi0,
                                                                                        offSets, offSetsSys,
                                                                                        graphRpPbCombStatPi0, graphRpPbCombSystPi0,
                                                                                        fileNamePi0RpPbOutputWeighting, "pPb_8TeV", "Pi0RpPb", kTRUE,
                                                                                        // NULL, fileNameCorrFactors
                                                                                        NULL, ""
                                                                                        );
    if (graphRpPbCombCombPi0 == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    while (graphRpPbCombStatPi0->GetX()[0] < 0.3){
        graphRpPbCombStatPi0->RemovePoint(0);
    }
    while (graphRpPbCombSystPi0->GetX()[0] < 0.3){
        graphRpPbCombSystPi0->RemovePoint(0);
    }
    while (graphRpPbCombCombPi0->GetX()[0] < 0.3){
        graphRpPbCombCombPi0->RemovePoint(0);
    }
    // graphRpPbCombCombPi0->Print();
    if(1){
      cout << "comb graph pi0 BEFORE:" << endl;
      graphRpPbCombCombPi0->Print();
      graphRpPbCombCombPi0->SetPoint(33,graphRpPbCombCombPi0->GetX()[33],graphRpPbCombCombPi0->GetY()[33]/1.25);
      graphRpPbCombCombPi0->SetPointError(33,graphRpPbCombCombPi0->GetEXlow()[33],graphRpPbCombCombPi0->GetEXhigh()[33],graphRpPbCombCombPi0->GetEYlow()[33]/1.25,graphRpPbCombCombPi0->GetEYhigh()[33]/1.25);
      graphRpPbCombCombPi0->SetPoint(34,graphRpPbCombCombPi0->GetX()[34],graphRpPbCombCombPi0->GetY()[34]/1.25);
      graphRpPbCombCombPi0->SetPointError(34,graphRpPbCombCombPi0->GetEXlow()[34],graphRpPbCombCombPi0->GetEXhigh()[34],graphRpPbCombCombPi0->GetEYlow()[34]/1.25,graphRpPbCombCombPi0->GetEYhigh()[34]/1.25);
      graphRpPbCombCombPi0->SetPoint(35,graphRpPbCombCombPi0->GetX()[35],graphRpPbCombCombPi0->GetY()[35]/1.25);
      graphRpPbCombCombPi0->SetPointError(35,graphRpPbCombCombPi0->GetEXlow()[35],graphRpPbCombCombPi0->GetEXhigh()[35],graphRpPbCombCombPi0->GetEYlow()[35]/1.25,graphRpPbCombCombPi0->GetEYhigh()[35]/1.25);

      graphRpPbCombSystPi0->SetPoint(33,graphRpPbCombSystPi0->GetX()[33],graphRpPbCombSystPi0->GetY()[33]/1.25);
      graphRpPbCombSystPi0->SetPointError(33,graphRpPbCombSystPi0->GetEXlow()[33],graphRpPbCombSystPi0->GetEXhigh()[33],graphRpPbCombSystPi0->GetEYlow()[33]/1.25,graphRpPbCombSystPi0->GetEYhigh()[33]/1.25);
      graphRpPbCombSystPi0->SetPoint(34,graphRpPbCombSystPi0->GetX()[34],graphRpPbCombSystPi0->GetY()[34]/1.25);
      graphRpPbCombSystPi0->SetPointError(34,graphRpPbCombSystPi0->GetEXlow()[34],graphRpPbCombSystPi0->GetEXhigh()[34],graphRpPbCombSystPi0->GetEYlow()[34]/1.25,graphRpPbCombSystPi0->GetEYhigh()[34]/1.25);
      graphRpPbCombSystPi0->SetPoint(35,graphRpPbCombSystPi0->GetX()[35],graphRpPbCombSystPi0->GetY()[35]/1.25);
      graphRpPbCombSystPi0->SetPointError(35,graphRpPbCombSystPi0->GetEXlow()[35],graphRpPbCombSystPi0->GetEXhigh()[35],graphRpPbCombSystPi0->GetEYlow()[35]/1.25,graphRpPbCombSystPi0->GetEYhigh()[35]/1.25);
      
      graphRpPbCombStatPi0->SetPoint(33,graphRpPbCombStatPi0->GetX()[33],graphRpPbCombStatPi0->GetY()[33]/1.25);
      graphRpPbCombStatPi0->SetPointError(33,graphRpPbCombStatPi0->GetEXlow()[33],graphRpPbCombStatPi0->GetEXhigh()[33],graphRpPbCombStatPi0->GetEYlow()[33]/1.25,graphRpPbCombStatPi0->GetEYhigh()[33]/1.25);
      graphRpPbCombStatPi0->SetPoint(34,graphRpPbCombStatPi0->GetX()[34],graphRpPbCombStatPi0->GetY()[34]/1.25);
      graphRpPbCombStatPi0->SetPointError(34,graphRpPbCombStatPi0->GetEXlow()[34],graphRpPbCombStatPi0->GetEXhigh()[34],graphRpPbCombStatPi0->GetEYlow()[34]/1.25,graphRpPbCombStatPi0->GetEYhigh()[34]/1.25);
      graphRpPbCombStatPi0->SetPoint(35,graphRpPbCombStatPi0->GetX()[35],graphRpPbCombStatPi0->GetY()[35]/1.25);
      graphRpPbCombStatPi0->SetPointError(35,graphRpPbCombStatPi0->GetEXlow()[35],graphRpPbCombStatPi0->GetEXhigh()[35],graphRpPbCombStatPi0->GetEYlow()[35]/1.25,graphRpPbCombStatPi0->GetEYhigh()[35]/1.25);
      cout << "comb graph pi0 AFTER:" << endl;
      graphRpPbCombCombPi0->Print();
      // return;
    }
    // return;
    // Reading weights from output file for plotting
    ifstream fileWeightsPi0RpPbRead;
    fileWeightsPi0RpPbRead.open(fileNamePi0RpPbOutputWeighting,ios_base::in);
    cout << "reading" << fileNamePi0RpPbOutputWeighting << endl;
    Double_t xValuesPi0RpPbRead[60];
    Double_t weightsPi0RpPbRead[11][60];
    Int_t availablePi0RpPbMeas[11]      = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetPi0RpPb               = nNuclModFacAvailPi0;
    Int_t nPtBinsPi0RpPbRead            = 0;
    while(!fileWeightsPi0RpPbRead.eof() && nPtBinsPi0RpPbRead < 60){
        TString garbage             = "";
        if (nPtBinsPi0RpPbRead == 0){
            fileWeightsPi0RpPbRead >> garbage ;//>> availablePi0RpPbMeas[0] >> availablePi0RpPbMeas[1] >> availablePi0RpPbMeas[2] >> availablePi0RpPbMeas[3];
            for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                fileWeightsPi0RpPbRead >> availablePi0RpPbMeas[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availablePi0RpPbMeas[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsPi0RpPbRead >> xValuesPi0RpPbRead[nPtBinsPi0RpPbRead-1];
            for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                fileWeightsPi0RpPbRead >> weightsPi0RpPbRead[availablePi0RpPbMeas[i]][nPtBinsPi0RpPbRead-1] ;
            }
            cout << "read: "<<  nPtBinsPi0RpPbRead << "\t"<< xValuesPi0RpPbRead[nPtBinsPi0RpPbRead-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                cout << weightsPi0RpPbRead[availablePi0RpPbMeas[i]][nPtBinsPi0RpPbRead-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsPi0RpPbRead++;
    }
    nPtBinsPi0RpPbRead                  = nPtBinsPi0RpPbRead-2 ;
    fileWeightsPi0RpPbRead.close();

    for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
        graphWeightsPi0RpPb[availablePi0RpPbMeas[i]]                        = new TGraph(nPtBinsPi0RpPbRead,xValuesPi0RpPbRead,weightsPi0RpPbRead[availablePi0RpPbMeas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsPi0RpPbRead; n++){
            if (graphWeightsPi0RpPb[availablePi0RpPbMeas[i]]->GetY()[bin] == 0) graphWeightsPi0RpPb[availablePi0RpPbMeas[i]]->RemovePoint(bin);
            else bin++;
        }
    }
  cout << __LINE__ << endl;
    // **********************************************************************************************************************
    // ******************************************* Plotting weights method only EMC *****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel           = 900*0.04;
    canvasWeights->cd();

    histo2DWeights->Draw("copy");
        TLegend* legendWeightsPi0RpPb   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetPi0RpPb+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
        for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
            DrawGammaSetMarkerTGraph(graphWeightsPi0RpPb[availablePi0RpPbMeas[i]], markerStyleDet[availablePi0RpPbMeas[i]], markerSizeDet[availablePi0RpPbMeas[i]]*0.5, colorDet[availablePi0RpPbMeas[i]] , colorDet[availablePi0RpPbMeas[i]]);
            graphWeightsPi0RpPb[availablePi0RpPbMeas[i]]->Draw("p,same,z");
            legendWeightsPi0RpPb->AddEntry(graphWeightsPi0RpPb[availablePi0RpPbMeas[i]],nameMeasGlobalPlot[availablePi0RpPbMeas[i]],"p");
        }
        legendWeightsPi0RpPb->Draw();

        labelWeightsEnergy->Draw();
        TLatex *labelWeightsPi0RpPb         = new TLatex(0.95,0.15,"#it{R}_{pA}: #pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsPi0RpPb, 0.85*textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelWeightsPi0RpPb->Draw();

        DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/Pi0RpPb_Weights.%s",outputDir.Data(),suffix.Data()));
  cout << __LINE__ << endl;

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelSysErr->cd();
  cout << __LINE__ << endl;

      histo2DRelSysErr->GetYaxis()->SetRangeUser(0,49.5);
      histo2DRelSysErr->Draw("copy");
      TLegend* legendRelSysErrRpPb       = GetAndSetLegend2(0.60, 0.92-(0.04*(nMeasSetPi0RpPb)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
      for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
          if (!sysErrorRelCollectionPi0RpA[availablePi0RpPbMeas[i]]) continue;
          cout << i << endl;
          DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0RpA[availablePi0RpPbMeas[i]], markerStyleDet[availablePi0RpPbMeas[i]], markerSizeDet[availablePi0RpPbMeas[i]]*0.5,
                                    colorDet[availablePi0RpPbMeas[i]], colorDet[availablePi0RpPbMeas[i]]);
          sysErrorRelCollectionPi0RpA[availablePi0RpPbMeas[i]]->Draw("p,same,z");
          legendRelSysErrRpPb->AddEntry(sysErrorRelCollectionPi0RpA[availablePi0RpPbMeas[i]],nameMeasGlobalPlot[availablePi0RpPbMeas[i]],"p");
      }
      legendRelSysErrRpPb->Draw();

      labelRelSysErrEnergy->Draw();
      TLatex *labelRelSysErrPi0pPb       = new TLatex(0.15,0.85,"#it{R}_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelSysErrPi0pPb, textSizeLabelsPixel, 4, 1, 43);
      labelRelSysErrPi0pPb->Draw();

    canvasRelSysErr->SaveAs(Form("%s/Pi0RpPb_RelSysErr_All.%s",outputDir.Data(),suffix.Data()));
  cout << __LINE__ << endl;

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelStatErr->cd();

      histo2DRelStatErr->GetYaxis()->SetRangeUser(0,49.5);
      histo2DRelStatErr->Draw("copy");
      TLegend* legendRelStatErrRpPb       = GetAndSetLegend2(0.12, 0.92-(0.04*(nMeasSetPi0RpPb)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
      for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
          if (!statErrorRelCollectionPi0RpA[availablePi0RpPbMeas[i]]) continue;
          DrawGammaSetMarkerTGraph(statErrorRelCollectionPi0RpA[availablePi0RpPbMeas[i]], markerStyleDet[availablePi0RpPbMeas[i]], markerSizeDet[availablePi0RpPbMeas[i]]*0.5,
                                    colorDet[availablePi0RpPbMeas[i]], colorDet[availablePi0RpPbMeas[i]]);
          statErrorRelCollectionPi0RpA[availablePi0RpPbMeas[i]]->Draw("p,same,z");
          legendRelStatErrRpPb->AddEntry(statErrorRelCollectionPi0RpA[availablePi0RpPbMeas[i]],nameMeasGlobalPlot[availablePi0RpPbMeas[i]],"p");
      }
      legendRelStatErrRpPb->Draw();

      labelRelStatErrEnergy->Draw();
      TLatex *labelRelStatErrPi0RpPb      = new TLatex(0.95,0.85,"#it{R}_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelStatErrPi0RpPb, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
      labelRelStatErrPi0RpPb->Draw();

    canvasRelStatErr->SaveAs(Form("%s/Pi0RpPb_RelStatErr_All.%s",outputDir.Data(),suffix.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TGraphAsymmErrors* graphCombPi0RpPbRelStat      = CalculateRelErrUpAsymmGraph( graphRpPbCombStatPi0, "relativeStatErrorPi0RpPb_MethodWPP");
    TGraphAsymmErrors* graphCombPi0RpPbRelSys       = CalculateRelErrUpAsymmGraph( graphRpPbCombSystPi0, "relativeSysErrorPi0RpPb_MethodWPP");
    TGraphAsymmErrors* graphCombPi0RpPbRelTot       = CalculateRelErrUpAsymmGraph( graphRpPbCombCombPi0, "relativeTotalErrorPi0RpPb_MethodWPP");
    canvasRelTotErr->cd();

      histo2DRelTotErr->GetYaxis()->SetRangeUser(0,49.5);
      histo2DRelTotErr->Draw("copy");
      DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
      graphCombPi0RpPbRelTot->Draw("p,same,z");
      DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
      graphCombPi0RpPbRelStat->Draw("l,x0,same,e1");
      DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
      graphCombPi0RpPbRelSys->SetLineStyle(7);
      graphCombPi0RpPbRelSys->Draw("l,x0,same,e1");
      legendRelTotErr2->Draw();

      labelRelTotErrEnergy->Draw();
      TLatex *labelRelTotErrPi0RpPb       = new TLatex(0.95,0.85,"#it{R}_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelTotErrPi0RpPb, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
      labelRelTotErrPi0RpPb->Draw();

    canvasRelTotErr->SaveAs(Form("%s/Pi0RpPb_Reldecomp_All.%s",outputDir.Data(),suffix.Data()));
  }




   // *******************************************************************************************************
  // *************************** RpPb calculation for eta for individual measurements **********************
  // *******************************************************************************************************
  TGraphAsymmErrors* graphRpPbIndStatEta[11];
  TH1D* histRpPbIndStatEta[11];
  TGraphAsymmErrors* graphRpPbIndSystEta[11];
  TGraphAsymmErrors* systErrorUnCorrInterCollectionEtaPP[11];
  TGraphAsymmErrors* graphRpPbIndCombEta[11];
  TGraph* graphWeightsEtaRpPb[11];
  TGraphAsymmErrors* statErrorRelCollectionEtaRpA[11];
  TGraphAsymmErrors* sysErrorRelCollectionEtaRpA[11];


  for (Int_t i = 0; i < 11; i++){
    graphWeightsEtaRpPb[i]                              = NULL;
    graphRpPbIndStatEta[i]                              = NULL;
    histRpPbIndStatEta[i]                               = NULL;
    graphRpPbIndSystEta[i]                              = NULL;
    graphRpPbIndCombEta[i]                              = NULL;
    systErrorUnCorrInterCollectionEtaPP[i]              = NULL;
    statErrorRelCollectionEtaRpA[i]                     = NULL;
    sysErrorRelCollectionEtaRpA[i]                      = NULL;
    if (haveRefPPEta[i] &&  directoryEta[i] && haveEtaNuclModFac[i]){
      histRpPbIndStatEta[i]                       = (TH1D*)histoEtaNuclModFac[i]->Clone(Form("histRpPbInd%sStatEta",nameMeasGlobalPlot[i].Data()));
      graphRpPbIndStatEta[i]                       = (TGraphAsymmErrors*)graphEtaNuclModFacStat[i]->Clone(Form("graphRpPbInd%sStatEta",nameMeasGlobalPlot[i].Data()));
      graphRpPbIndSystEta[i]                       = (TGraphAsymmErrors*)graphEtaNuclModFacSys[i]->Clone(Form("graphRpPbInd%sSystEta",nameMeasGlobalPlot[i].Data()));

        cout << "printing cleaned combined meas RpA for method: " << i << endl;
        histRpPbIndStatEta[i]                       = (TH1D*)statErrorCollectionEta[i]->Clone(Form("histRpPb%sStatEta",nameMeasGlobalPlot[i].Data()));
        for (Int_t j = 0; j< histRpPbIndStatEta[i]->GetNbinsX()+1; j++){
          histRpPbIndStatEta[i]->SetBinContent(j,0);
          histRpPbIndStatEta[i]->SetBinError(j,0);
        }
        for (Int_t j = 0; j < graphRpPbIndStatEta[i]->GetN(); j++){
          Int_t bin                               = histRpPbIndStatEta[i]->GetXaxis()->FindBin(graphRpPbIndStatEta[i]->GetX()[j]);
          histRpPbIndStatEta[i]->SetBinContent(bin,graphRpPbIndStatEta[i]->GetY()[j]);
          histRpPbIndStatEta[i]->SetBinError(bin,graphRpPbIndStatEta[i]->GetEYlow()[j]);
        }

      if (histRpPbIndStatEta[i]){
        statErrorRelCollectionEtaRpA[i]     = new TGraphAsymmErrors(histRpPbIndStatEta[i]);
        while (statErrorRelCollectionEtaRpA[i]->GetY()[0] == 0 && statErrorRelCollectionEtaRpA[i]->GetN() > 0 )
          statErrorRelCollectionEtaRpA[i]->RemovePoint(0);
        while (statErrorRelCollectionEtaRpA[i]->GetY()[statErrorRelCollectionEtaRpA[i]->GetN()-1] == 0 && statErrorRelCollectionEtaRpA[i]->GetN() > 0 )
          statErrorRelCollectionEtaRpA[i]->RemovePoint(statErrorRelCollectionEtaRpA[i]->GetN()-1);
        if (statErrorRelCollectionEtaRpA[i]->GetN() > 0 )
          statErrorRelCollectionEtaRpA[i] = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEtaRpA[i], Form("relativeStatErrorEtaRpPb_%s", nameMeasGlobal[i].Data()));
        else
          statErrorRelCollectionEtaRpA[i] = NULL;
      }
      if (graphRpPbIndSystEta[i]){
        sysErrorRelCollectionEtaRpA[i]      = (TGraphAsymmErrors*)graphRpPbIndSystEta[i]->Clone(Form("relativeSysErrorEtaRpPb_%s", nameMeasGlobal[i].Data()));
        while (sysErrorRelCollectionEtaRpA[i]->GetY()[0]== 0 && sysErrorRelCollectionEtaRpA[i]->GetN() > 0 )
          sysErrorRelCollectionEtaRpA[i]->RemovePoint(0);
        while (sysErrorRelCollectionEtaRpA[i]->GetY()[sysErrorRelCollectionEtaRpA[i]->GetN()-1] == 0 && sysErrorRelCollectionEtaRpA[i]->GetN() > 0 )
          sysErrorRelCollectionEtaRpA[i]->RemovePoint(sysErrorRelCollectionEtaRpA[i]->GetN()-1);
        if (sysErrorRelCollectionEtaRpA[i]->GetN() > 0 )
          sysErrorRelCollectionEtaRpA[i]  = CalculateRelErrUpAsymmGraph( graphRpPbIndSystEta[i], Form("relativeSysErrorEtaRpPb_%s", nameMeasGlobal[i].Data()));
        else
          sysErrorRelCollectionEtaRpA[i]  = NULL;

      }
    }
  }
  // return;


  // *******************************************************************************************************
  // *************************** RpPb calculation for Eta **************************************************
  // *******************************************************************************************************

  TGraphAsymmErrors* graphRpPbCombStatEta         = NULL;
  TGraphAsymmErrors* graphRpPbCombSystEta         = NULL;
  TGraphAsymmErrors* graphRpPbCombCombEta         = NULL;

  if(bWCorrection.Contains("Y") ){
    TString fileNameEtaRpPbOutputWeighting       = Form("%s/EtaRpPb_WeightingMethod.dat",outputDir.Data());
    graphRpPbCombCombEta                         = CombinePtPointsSpectraFullCorrMat(   histRpPbIndStatEta,    graphRpPbIndSystEta,
                                                                                        xPtLimitsEta, nBinsEta,
                                                                                        offSetsEta, offSetsSysEta,
                                                                                        graphRpPbCombStatEta, graphRpPbCombSystEta,
                                                                                        fileNameEtaRpPbOutputWeighting, "pPb_8TeV", "EtaRpPb", kTRUE,
                                                                                        // NULL, fileNameCorrFactors
                                                                                        NULL, ""
                                                                                        );
    if (graphRpPbCombCombEta == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    while (graphRpPbCombStatEta->GetX()[0] < 0.3){
        graphRpPbCombStatEta->RemovePoint(0);
    }
    while (graphRpPbCombSystEta->GetX()[0] < 0.3){
        graphRpPbCombSystEta->RemovePoint(0);
    }
    while (graphRpPbCombCombEta->GetX()[0] < 0.3){
        graphRpPbCombCombEta->RemovePoint(0);
    }
    graphRpPbCombCombEta->Print();
    // return;
    // Reading weights from output file for plotting
    ifstream fileWeightsEtaRpPbRead;
    fileWeightsEtaRpPbRead.open(fileNameEtaRpPbOutputWeighting,ios_base::in);
    cout << "reading" << fileNameEtaRpPbOutputWeighting << endl;
    Double_t xValuesEtaRpPbRead[50];
    Double_t weightsEtaRpPbRead[11][50];
    Int_t availableEtaRpPbMeas[11]      = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetEtaRpPb               = nNuclModFacAvailEta;
    Int_t nPtBinsEtaRpPbRead            = 0;
    while(!fileWeightsEtaRpPbRead.eof() && nPtBinsEtaRpPbRead < 50){
        TString garbage             = "";
        if (nPtBinsEtaRpPbRead == 0){
            fileWeightsEtaRpPbRead >> garbage ;//>> availableEtaRpPbMeas[0] >> availableEtaRpPbMeas[1] >> availableEtaRpPbMeas[2] >> availableEtaRpPbMeas[3];
            for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                fileWeightsEtaRpPbRead >> availableEtaRpPbMeas[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableEtaRpPbMeas[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsEtaRpPbRead >> xValuesEtaRpPbRead[nPtBinsEtaRpPbRead-1];
            for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                fileWeightsEtaRpPbRead >> weightsEtaRpPbRead[availableEtaRpPbMeas[i]][nPtBinsEtaRpPbRead-1] ;
            }
            cout << "read: "<<  nPtBinsEtaRpPbRead << "\t"<< xValuesEtaRpPbRead[nPtBinsEtaRpPbRead-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                cout << weightsEtaRpPbRead[availableEtaRpPbMeas[i]][nPtBinsEtaRpPbRead-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsEtaRpPbRead++;
    }
    nPtBinsEtaRpPbRead                  = nPtBinsEtaRpPbRead-2 ;
    fileWeightsEtaRpPbRead.close();

    for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
        graphWeightsEtaRpPb[availableEtaRpPbMeas[i]]                        = new TGraph(nPtBinsEtaRpPbRead,xValuesEtaRpPbRead,weightsEtaRpPbRead[availableEtaRpPbMeas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsEtaRpPbRead; n++){
            if (graphWeightsEtaRpPb[availableEtaRpPbMeas[i]]->GetY()[bin] == 0) graphWeightsEtaRpPb[availableEtaRpPbMeas[i]]->RemovePoint(bin);
            else bin++;
        }
    }
  cout << __LINE__ << endl;
    // **********************************************************************************************************************
    // ******************************************* Plotting weights method only EMC *****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel           = 900*0.04;
    canvasWeights->cd();

    histo2DWeights->Draw("copy");
        TLegend* legendWeightsEtaRpPb   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetEtaRpPb+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
        for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaRpPb[availableEtaRpPbMeas[i]], markerStyleDet[availableEtaRpPbMeas[i]], markerSizeDet[availableEtaRpPbMeas[i]]*0.5, colorDet[availableEtaRpPbMeas[i]] , colorDet[availableEtaRpPbMeas[i]]);
            graphWeightsEtaRpPb[availableEtaRpPbMeas[i]]->Draw("p,same,z");
            legendWeightsEtaRpPb->AddEntry(graphWeightsEtaRpPb[availableEtaRpPbMeas[i]],nameMeasGlobalPlot[availableEtaRpPbMeas[i]],"p");
        }
        legendWeightsEtaRpPb->Draw();

        labelWeightsEnergy->Draw();
        TLatex *labelWeightsEtaRpPb         = new TLatex(0.95,0.15,"#it{R}_{pA}: #pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsEtaRpPb, 0.85*textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelWeightsEtaRpPb->Draw();

        DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/EtaRpPb_Weights.%s",outputDir.Data(),suffix.Data()));
  cout << __LINE__ << endl;

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelSysErr->cd();
  cout << __LINE__ << endl;

      histo2DRelSysErr->GetYaxis()->SetRangeUser(0,49.5);
      histo2DRelSysErr->Draw("copy");
      TLegend* legendRelSysErrRpPb       = GetAndSetLegend2(0.60, 0.92-(0.04*(nMeasSetEtaRpPb)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
      for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
          if (!sysErrorRelCollectionEtaRpA[availableEtaRpPbMeas[i]]) continue;
          cout << i << endl;
          DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaRpA[availableEtaRpPbMeas[i]], markerStyleDet[availableEtaRpPbMeas[i]], markerSizeDet[availableEtaRpPbMeas[i]]*0.5,
                                    colorDet[availableEtaRpPbMeas[i]], colorDet[availableEtaRpPbMeas[i]]);
          sysErrorRelCollectionEtaRpA[availableEtaRpPbMeas[i]]->Draw("p,same,z");
          legendRelSysErrRpPb->AddEntry(sysErrorRelCollectionEtaRpA[availableEtaRpPbMeas[i]],nameMeasGlobalPlot[availableEtaRpPbMeas[i]],"p");
      }
      legendRelSysErrRpPb->Draw();

      labelRelSysErrEnergy->Draw();
      TLatex *labelRelSysErrEtapPb       = new TLatex(0.15,0.85,"#it{R}_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelSysErrEtapPb, textSizeLabelsPixel, 4, 1, 43);
      labelRelSysErrEtapPb->Draw();

    canvasRelSysErr->SaveAs(Form("%s/EtaRpPb_RelSysErr_All.%s",outputDir.Data(),suffix.Data()));
  cout << __LINE__ << endl;

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelStatErr->cd();

      histo2DRelStatErr->GetYaxis()->SetRangeUser(0,49.5);
      histo2DRelStatErr->Draw("copy");
      TLegend* legendRelStatErrRpPb       = GetAndSetLegend2(0.12, 0.92-(0.04*(nMeasSetEtaRpPb)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
      for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
          if (!statErrorRelCollectionEtaRpA[availableEtaRpPbMeas[i]]) continue;
          DrawGammaSetMarkerTGraph(statErrorRelCollectionEtaRpA[availableEtaRpPbMeas[i]], markerStyleDet[availableEtaRpPbMeas[i]], markerSizeDet[availableEtaRpPbMeas[i]]*0.5,
                                    colorDet[availableEtaRpPbMeas[i]], colorDet[availableEtaRpPbMeas[i]]);
          statErrorRelCollectionEtaRpA[availableEtaRpPbMeas[i]]->Draw("p,same,z");
          legendRelStatErrRpPb->AddEntry(statErrorRelCollectionEtaRpA[availableEtaRpPbMeas[i]],nameMeasGlobalPlot[availableEtaRpPbMeas[i]],"p");
      }
      legendRelStatErrRpPb->Draw();

      labelRelStatErrEnergy->Draw();
      TLatex *labelRelStatErrEtaRpPb      = new TLatex(0.95,0.85,"#it{R}_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelStatErrEtaRpPb, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
      labelRelStatErrEtaRpPb->Draw();

    canvasRelStatErr->SaveAs(Form("%s/EtaRpPb_RelStatErr_All.%s",outputDir.Data(),suffix.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TGraphAsymmErrors* graphCombEtaRpPbRelStat      = CalculateRelErrUpAsymmGraph( graphRpPbCombStatEta, "relativeStatErrorEtaRpPb_MethodWPP");
    TGraphAsymmErrors* graphCombEtaRpPbRelSys       = CalculateRelErrUpAsymmGraph( graphRpPbCombSystEta, "relativeSysErrorEtaRpPb_MethodWPP");
    TGraphAsymmErrors* graphCombEtaRpPbRelTot       = CalculateRelErrUpAsymmGraph( graphRpPbCombCombEta, "relativeTotalErrorEtaRpPb_MethodWPP");
    canvasRelTotErr->cd();

      histo2DRelTotErr->GetYaxis()->SetRangeUser(0,49.5);
      histo2DRelTotErr->Draw("copy");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
      graphCombEtaRpPbRelTot->Draw("p,same,z");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
      graphCombEtaRpPbRelStat->Draw("l,x0,same,e1");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
      graphCombEtaRpPbRelSys->SetLineStyle(7);
      graphCombEtaRpPbRelSys->Draw("l,x0,same,e1");
      legendRelTotErr2->Draw();

      labelRelTotErrEnergy->Draw();
      TLatex *labelRelTotErrEtaRpPb       = new TLatex(0.95,0.85,"#it{R}_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelTotErrEtaRpPb, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
      labelRelTotErrEtaRpPb->Draw();

    canvasRelTotErr->SaveAs(Form("%s/EtaRpPb_Reldecomp_All.%s",outputDir.Data(),suffix.Data()));
  }


  // *************************************************************************************************************
  // redo fitting after binshifts
  // *************************************************************************************************************
  // Tsallis function
  graphCombPi0InvYieldTot->Fit(fitInvXSectionPi0,"QNRMEX0+","",0.4,200.);
  fitInvXSectionPi0           = FitObject("l","fitInvXSectionPi0pPb8TeV","Pi0",graphCombPi0InvYieldTot,0.3,200.,paramGraph,"QNRMEX0+");
  cout << WriteParameterToFile(fitInvXSectionPi0)<< endl;

  //Two component model from Bylinkin
  Double_t paramTCMPi0New[5]  = { graphCombPi0InvYieldTot->GetY()[0],0.1,
                          graphCombPi0InvYieldTot->GetY()[4],0.6,3.0};
  TF1* fitTCMInvYieldPi0        = FitObject("tcm","fitTCMInvYieldPi0pPb8TeV","Pi0",graphCombPi0InvYieldTot,0.4,200. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
  cout << WriteParameterToFile(fitTCMInvYieldPi0)<< endl;

  Double_t paramPi0Power[3] = {1E2,0.5,7.5};
  TF1* fitPowInvXSectionPi0   = FitObject("powPure","fitPowInvXSectionPi0pPb8TeV","Pi0",graphCombPi0InvYieldTot,3.5,200. ,NULL,"QNRMEX0+","", kFALSE);
  cout << WriteParameterToFile(fitPowInvXSectionPi0)<< endl;

  TF1* fitPowInvXSectionPi0Stat   = FitObject("powPure","fitPowInvXSectionPi0pPb8TeVStat","Pi0",graphCombPi0InvYieldStat,3.5,200. ,paramPi0Power,"QNRMEX0+","", kFALSE);
  cout << WriteParameterToFile(fitPowInvXSectionPi0Stat)<< endl;

  Double_t paramPi0HageDorn[5] = {1E1,0.3,-0.1,0.5,5.95};
  TF1* fitOHagInvYieldPi0Tot   = FitObject("oHag","fitOHagInvYieldPi0pPb8TeV","Pi0",graphCombPi0InvYieldTot,0.4,200. ,paramPi0HageDorn,"QNRMEX0+","", kFALSE);
  cout << WriteParameterToFile(fitOHagInvYieldPi0Tot)<< endl;

  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Pi0 - Tsallis" << endl;
  fLog << WriteParameterToFile(fitInvXSectionPi0)<< endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Pi0 - TCM" << endl;
  fLog << WriteParameterToFile(fitTCMInvYieldPi0) << endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Pi0 - Hagedorn" << endl;
  fLog << WriteParameterToFile(fitOHagInvYieldPi0Tot) << endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Pi0 - PowerLaw" << endl;
  fLog << WriteParameterToFile(fitPowInvXSectionPi0) << endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Pi0 - PowerLaw - Stat" << endl;
  fLog << WriteParameterToFile(fitPowInvXSectionPi0Stat) << endl;

  // Tsallis function
  graphCombEtaInvYieldTot->Fit(fitInvXSectionEta,"QNRMEX0+","",0.4,18.);
  fitInvXSectionEta        = FitObject("l","fitInvXSectionEtapPb8TeV","Eta",graphCombEtaInvYieldTot,0.4,40.,paramGraphEta,"QNRMEX0+");
  cout << WriteParameterToFile(fitInvXSectionEta)<< endl;

  Double_t paramTCMEta[5]  = {graphCombEtaInvYieldTot->GetY()[0],0.2,graphCombEtaInvYieldTot->GetY()[3],0.75,3.};
  // Two component model by Bylinkin
  TF1* fitTCMInvYieldEta= FitObject("tcm","fitTCMInvYieldEtapPb8TeV","Eta",graphCombEtaInvYieldTot,0.5,40.,paramTCMEta,"QNRMEX0+","", kFALSE);
  fitTCMInvYieldEta     = FitObject("tcm","fitTCMInvYieldEtapPb8TeV","Eta",graphCombEtaInvYieldTot,0.5,40.,paramTCMEta,"QNRMEX0+","", kFALSE);

  TF1* fitTCMDecomposedEtaL                 = new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mesonMassExpectEta,mesonMassExpectEta,mesonMassExpectEta));
  fitTCMDecomposedEtaL->SetParameters(fitTCMInvYieldEta->GetParameter(0),fitTCMInvYieldEta->GetParameter(1));
  fitTCMDecomposedEtaL->SetRange(minPtEta,maxPtEta);
  TF1 *fitTCMDecomposedEtaH                 = new TF1("twoCompModel_DecH","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]))");
  fitTCMDecomposedEtaH->SetParameters(fitTCMInvYieldEta->GetParameter(2),fitTCMInvYieldEta->GetParameter(3), fitTCMInvYieldEta->GetParameter(4));
  fitTCMDecomposedEtaH->SetRange(minPtEta,maxPtEta);
  cout << WriteParameterToFile(fitTCMInvYieldEta)<< endl;

  Double_t paramEtaPower[3] = {1E11,0.5,6.5};
  TF1* fitPowInvXSectionEta   = FitObject("powPure","fitPowInvXSectionEtapPb8TeV","Eta",graphCombEtaInvYieldTot,3.5,40. ,NULL,"QNRMEX0+","", kFALSE);
  cout << WriteParameterToFile(fitPowInvXSectionEta)<< endl;

  Double_t paramEtaHageDorn[5] = {1E11,0.3,-0.1,0.5,5.95};
  TF1* fitOHagInvYieldEtaTot   = FitObject("oHag","fitOHagInvYieldEtapPb8TeV","Eta",graphCombEtaInvYieldTot,0.4,40. ,paramEtaHageDorn,"QNRMEX0+","", kFALSE);
  cout << WriteParameterToFile(fitOHagInvYieldEtaTot)<< endl;

  TF1* fitPowInvXSectionEtaStat   = FitObject("powPure","fitPowInvXSectionEtapPb8TeVStat","Eta",graphCombEtaInvYieldStat,3.5,40. ,NULL,"QNRMEX0+","", kFALSE);
  cout << WriteParameterToFile(fitPowInvXSectionEtaStat)<< endl;

  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Eta - Tsallis" << endl;
  fLog << WriteParameterToFile(fitInvXSectionEta)<< endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Eta - TCM" << endl;
  fLog << WriteParameterToFile(fitTCMInvYieldEta) << endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Eta - Hagedorn" << endl;
  fLog << WriteParameterToFile(fitOHagInvYieldEtaTot) << endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Eta - PowerLaw" << endl;
  fLog << WriteParameterToFile(fitPowInvXSectionEta) << endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Eta - PowerLaw - Stat" << endl;
  fLog << WriteParameterToFile(fitPowInvXSectionEtaStat) << endl;


  //********************************************************************************************************
  // Plotting simple comparison of data vs fit to eta meson spec
  //********************************************************************************************************
  TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
  DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
  canvasDummy2->SetLogy();
  canvasDummy2->SetLogx();
  TH2F* histo2DDummy3;
  histo2DDummy3               = new TH2F("histo2DDummy3","histo2DDummy3",1000,minPtEta,maxPtEta,1000,minInvYieldEta,maxInvYieldEta);
  SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 0.8,1.55);
  histo2DDummy3->DrawCopy();

  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatUnShifted, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
  graphCombEtaInvYieldStatUnShifted->Draw("pEsame");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
  graphCombEtaInvYieldStat->Draw("pEsame");

  //    fitInvXSectionEta->SetLineColor(kBlue+2);
  //    fitInvXSectionEta->Draw("same");
  fitTCMInvYieldEta->SetLineColor(kRed+2);
  fitTCMInvYieldEta->Draw("same");

  fitTCMDecomposedEtaL->SetLineColor(kAzure);
  fitTCMDecomposedEtaL->SetLineStyle(2);
  fitTCMDecomposedEtaL->Draw("same");
  fitTCMDecomposedEtaH->SetLineColor(kGreen+2);
  fitTCMDecomposedEtaH->SetLineStyle(8);
  fitTCMDecomposedEtaH->Draw("same");

  TLatex *labelTCMEta1= new TLatex(0.48, 0.94, Form("TCM low:"));
  TLatex *labelTCMEta2= new TLatex(0.48, 0.90, Form("A_{1}: (%.1e #pm %.1e) - T_{e}: (%.3f #pm %.3f)",fitTCMInvYieldEta->GetParameter(0),fitTCMInvYieldEta->GetParError(0),fitTCMInvYieldEta->GetParameter(1),fitTCMInvYieldEta->GetParError(1)));
  TLatex *labelTCMEta3= new TLatex(0.48, 0.86, Form("TCM high:"));
  TLatex *labelTCMEta4= new TLatex(0.48, 0.82, Form("A_{2}: (%.1e #pm %.1e) - T: (%.3f #pm %.3f) - n: (%.3f #pm %.3f)",fitTCMInvYieldEta->GetParameter(2),fitTCMInvYieldEta->GetParError(2),fitTCMInvYieldEta->GetParameter(3),fitTCMInvYieldEta->GetParError(3),fitTCMInvYieldEta->GetParameter(4),fitTCMInvYieldEta->GetParError(4)));

  TLatex *labelTCMEta5= new TLatex(0.55, 0.75, Form("Bylinkin-Rostovtsev:"));
  TLatex *labelTCMEta6= new TLatex(0.55, 0.71, Form("#it{A}_{1} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}_{2}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}"));

  SetStyleTLatex( labelTCMEta1, 0.03,4);
  labelTCMEta1->Draw();
  SetStyleTLatex( labelTCMEta2, 0.02,4);
  labelTCMEta2->Draw();
  SetStyleTLatex( labelTCMEta3, 0.03,4);
  labelTCMEta3->Draw();
  SetStyleTLatex( labelTCMEta4, 0.02,4);
  labelTCMEta4->Draw();
  SetStyleTLatex( labelTCMEta5, 0.03,4);
  labelTCMEta5->Draw();
  SetStyleTLatex( labelTCMEta6, 0.03,4);
  labelTCMEta6->Draw();

  TLatex *labelRelSysErrEnergyC    = new TLatex(0.18,0.94,collisionSystempPb8TeV.Data());
  SetStyleTLatex( labelRelSysErrEnergyC, 0.85*textSizeLabelsPixel,4);
  labelRelSysErrEnergyC->SetTextFont(43);
  labelRelSysErrEnergyC->Draw();
  TLatex *labelRelSysErrEtaC       = new TLatex(0.18,0.9,"#eta #rightarrow #gamma#gamma");
  SetStyleTLatex( labelRelSysErrEtaC, 0.85*textSizeLabelsPixel,4);
  labelRelSysErrEtaC->SetTextFont(43);
  labelRelSysErrEtaC->Draw();

  TLegend* legendWithFitEta   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
  legendWithFitEta->AddEntry(fitTCMDecomposedEtaL,"TCM low","l");
  legendWithFitEta->AddEntry(fitTCMDecomposedEtaH,"TCM high","l");
  legendWithFitEta->AddEntry(fitTCMInvYieldEta,"Bylinkin-Rostovtsev (TCM)","l");
  legendWithFitEta->Draw();

  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/ComparisonWithFitEta_pPb8TeV.%s",outputDir.Data(),suffix.Data()));
  //********************************************************************************************************
  canvasDummy2->Clear();
  histo2DDummy3->DrawCopy();

  graphCombEtaInvYieldStatUnShifted->Draw("pEsame");
  graphCombEtaInvYieldStat->Draw("pEsame");

  fitInvXSectionEta->SetLineColor(kRed+2);
  fitInvXSectionEta->Draw("same");

  TLatex *labelTCMEta10 = new TLatex(0.35, 0.90, Form("dN/dy: (%.1e #pm %.1e) - n: (%.3f #pm %.3f) - T_{Levy} (GeV/c): (%.3f #pm %.3f)",fitInvXSectionEta->GetParameter(0),fitInvXSectionEta->GetParError(0),fitInvXSectionEta->GetParameter(1),fitInvXSectionEta->GetParError(1),fitInvXSectionEta->GetParameter(2),fitInvXSectionEta->GetParError(2)));
  SetStyleTLatex( labelTCMEta10, 0.02,4);
  labelTCMEta10->Draw();

  labelRelSysErrEnergyC->Draw();
  labelRelSysErrEtaC->Draw();

  TLegend* legendWithFitEta2   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
  legendWithFitEta2->AddEntry(fitInvXSectionEta,"Levy-Tsallis","l");
  legendWithFitEta2->Draw();

  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/ComparisonWithFit_Tsallis_Eta_pPb8TeV.%s",outputDir.Data(),suffix.Data()));

  delete histo2DDummy3;
  canvasDummy2->Clear();

  //********************************************************************************************************
  // Plotting simple comparison of data vs fit to pi0 meson spec
  //********************************************************************************************************

  TF1* fitTCMDecomposedPi0L                 = new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
  fitTCMDecomposedPi0L->SetParameters(fitTCMInvYieldPi0->GetParameter(0),fitTCMInvYieldPi0->GetParameter(1));
  fitTCMDecomposedPi0L->SetRange(minPtPi0,maxPtPi0);
  TF1 *fitTCMDecomposedPi0H                 = new TF1("twoCompModel_DecH","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]))");
  //      graphCombEtaInvYieldTotA->Fit(fitTCMDecomposedH,"QNRMEX0+","",5,20);
  fitTCMDecomposedPi0H->SetParameters(fitTCMInvYieldPi0->GetParameter(2),fitTCMInvYieldPi0->GetParameter(3), fitTCMInvYieldPi0->GetParameter(4));
  fitTCMDecomposedPi0H->SetRange(minPtPi0,maxPtPi0);

  histo2DDummy3               = new TH2F("histo2DDummy2","histo2DDummy2",1000,minPtPi0,maxPtPi0,1000,minInvYieldPi0,maxInvYieldPi0);
  SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 0.8,1.55);
  histo2DDummy3->DrawCopy();

  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatUnShifted, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
  graphCombPi0InvYieldStatUnShifted->Draw("pEsame");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
  graphCombPi0InvYieldStat->Draw("pEsame");

  fitTCMInvYieldPi0->SetLineColor(kRed+2);
  fitTCMInvYieldPi0->SetRange(minPtPi0,maxPtPi0);
  fitTCMInvYieldPi0->Draw("same");

  fitTCMDecomposedPi0L->SetLineColor(kAzure);
  fitTCMDecomposedPi0L->SetLineStyle(2);
  fitTCMDecomposedPi0L->Draw("same");
  fitTCMDecomposedPi0H->SetLineColor(kGreen+2);
  fitTCMDecomposedPi0H->SetLineStyle(8);
  fitTCMDecomposedPi0H->Draw("same");

  TLatex *labelTCMPi01= new TLatex(0.48, 0.94, Form("TCM low:"));
  TLatex *labelTCMPi02= new TLatex(0.48, 0.90, Form("A_{1}: (%.1e #pm %.1e) - T_{e}: (%.3f #pm %.3f)",fitTCMInvYieldPi0->GetParameter(0),fitTCMInvYieldPi0->GetParError(0),fitTCMInvYieldPi0->GetParameter(1),fitTCMInvYieldPi0->GetParError(1)));
  TLatex *labelTCMPi03= new TLatex(0.48, 0.86, Form("TCM high:"));
  TLatex *labelTCMPi04= new TLatex(0.48, 0.82, Form("A_{2}: (%.1e #pm %.1e) - T: (%.3f #pm %.3f) - n: (%.3f #pm %.3f)",fitTCMInvYieldPi0->GetParameter(2),fitTCMInvYieldPi0->GetParError(2),abs(fitTCMInvYieldPi0->GetParameter(3)),fitTCMInvYieldPi0->GetParError(3),fitTCMInvYieldPi0->GetParameter(4),fitTCMInvYieldPi0->GetParError(4)));

  TLatex *labelTCMPi05= new TLatex(0.55, 0.75, Form("Bylinkin-Rostovtsev:"));
  TLatex *labelTCMPi06= new TLatex(0.55, 0.71, Form("#it{A}_{1} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}_{2}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}"));

  SetStyleTLatex( labelTCMPi01, 0.03,4);
  labelTCMPi01->Draw();
  SetStyleTLatex( labelTCMPi02, 0.02,4);
  labelTCMPi02->Draw();
  SetStyleTLatex( labelTCMPi03, 0.03,4);
  labelTCMPi03->Draw();
  SetStyleTLatex( labelTCMPi04, 0.02,4);
  labelTCMPi04->Draw();
  SetStyleTLatex( labelTCMPi05, 0.03,4);
  labelTCMPi05->Draw();
  SetStyleTLatex( labelTCMPi06, 0.03,4);
  labelTCMPi06->Draw();

  labelRelSysErrEnergyC->Draw();
  TLatex *labelRelSysErrPi0C       = new TLatex(0.18,0.9,"#pi^{0} #rightarrow #gamma#gamma");
  SetStyleTLatex( labelRelSysErrPi0C, 0.85*textSizeLabelsPixel,4);
  labelRelSysErrPi0C->SetTextFont(43);
  labelRelSysErrPi0C->Draw();

  TLegend* legendWithFit   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
  legendWithFit->AddEntry(fitTCMDecomposedPi0L,"TCM low","l");
  legendWithFit->AddEntry(fitTCMDecomposedPi0H,"TCM high","l");
  legendWithFit->AddEntry(fitTCMInvYieldPi0,"Bylinkin-Rostovtsev (TCM)","l");
  legendWithFit->Draw();

  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/ComparisonWithFitPi0_pPb8TeV.%s",outputDir.Data(),suffix.Data()));
  //********************************************************************************************************
  canvasDummy2->Clear();
  histo2DDummy3->DrawCopy();

  graphCombPi0InvYieldStatUnShifted->Draw("pEsame");
  graphCombPi0InvYieldStat->Draw("pEsame");

  fitInvXSectionPi0->SetLineColor(kRed+2);
  fitInvXSectionPi0->Draw("same");

  TLatex *labelTCMEta20 = new TLatex(0.35, 0.90, Form("dN/dy: (%.1e #pm %.1e) - n: (%.3f #pm %.3f) - T_{Levy} (GeV/c): (%.3f #pm %.3f)",fitInvXSectionPi0->GetParameter(0),fitInvXSectionPi0->GetParError(0),fitInvXSectionPi0->GetParameter(1),fitInvXSectionPi0->GetParError(1),fitInvXSectionPi0->GetParameter(2),fitInvXSectionPi0->GetParError(2)));
  SetStyleTLatex( labelTCMEta20, 0.02,4);
  labelTCMEta20->Draw();

  labelRelSysErrEnergyC->Draw();
  labelRelSysErrPi0C->Draw();

  TLegend* legendWithFitPi02   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
  legendWithFitPi02->AddEntry(fitInvXSectionPi0,"Levy-Tsallis","l");
  legendWithFitPi02->Draw();

  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/ComparisonWithFit_Tsallis_Pi0_pPb8TeV.%s",outputDir.Data(),suffix.Data()));

  delete canvasDummy2;
  delete histo2DDummy3;


  // *************************************************************************************************************
  // Calculate ratios to combined fit
  // *************************************************************************************************************
  TH1D* histoRatioPi0HIJINGToFit                      = (TH1D*) histoHIJINGPi0->Clone("histoRatioPi0HIJINGToFit");
  histoRatioPi0HIJINGToFit                            = CalculateHistoRatioToFit (histoRatioPi0HIJINGToFit, fitTCMInvYieldPi0);
  // histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(0.3, 20);

  TH1D* histoRatioEPOSJJToFit                     = (TH1D*) histoEPOSJJInvYield->Clone();
  histoRatioEPOSJJToFit                           = CalculateHistoRatioToFit (histoRatioEPOSJJToFit, fitTCMInvYieldPi0);
  TGraphErrors* graphRatioEPOSJJToFit             = (TGraphErrors*) graphEPOSJJInvYield->Clone();
  graphRatioEPOSJJToFit                           = CalculateGraphErrRatioToFit (graphRatioEPOSJJToFit, fitTCMInvYieldPi0);

  TH1D* histoRatioDPMJETToFit                     = (TH1D*) histoDPMJETInvYield->Clone();
  histoRatioDPMJETToFit                           = CalculateHistoRatioToFit (histoRatioDPMJETToFit, fitTCMInvYieldPi0);
  TGraphErrors* graphRatioDPMJETToFit             = (TGraphErrors*) graphDPMJETInvYield->Clone();
  graphRatioDPMJETToFit                           = CalculateGraphErrRatioToFit (graphRatioDPMJETToFit, fitTCMInvYieldPi0);

  TH1D* histoRatioPythia8EPPS16ToFit                     = (TH1D*) histoPythia8EPPS16InvYield->Clone();
  histoRatioPythia8EPPS16ToFit                           = CalculateHistoRatioToFit (histoRatioPythia8EPPS16ToFit, fitTCMInvYieldPi0);
  TGraphErrors* graphRatioPythia8EPPS16ToFit             = (TGraphErrors*) graphPythia8EPPS16InvYield->Clone();
  graphRatioPythia8EPPS16ToFit                           = CalculateGraphErrRatioToFit (graphRatioPythia8EPPS16ToFit, fitTCMInvYieldPi0);

  TH1D* histoRatioPythia8nCTEQ15ToFit                     = (TH1D*) histoPythia8nCTEQ15InvYield->Clone();
  histoRatioPythia8nCTEQ15ToFit                           = CalculateHistoRatioToFit (histoRatioPythia8nCTEQ15ToFit, fitTCMInvYieldPi0);
  TGraphErrors* graphRatioPythia8nCTEQ15ToFit             = (TGraphErrors*) graphPythia8nCTEQ15InvYield->Clone();
  graphRatioPythia8nCTEQ15ToFit                           = CalculateGraphErrRatioToFit (graphRatioPythia8nCTEQ15ToFit, fitTCMInvYieldPi0);

  TH1D* histoRatioPythia8Monash2013ToFit                     = (TH1D*) histoPythia8Monash2013InvSecpTh->Clone();
  histoRatioPythia8Monash2013ToFit                           = CalculateHistoRatioToFit (histoRatioPythia8Monash2013ToFit, fitTCMInvYieldPi0);
  TGraphErrors* graphRatioPythia8Monash2013ToFit             = (TGraphErrors*) graphPythia8Monash2013InvSecpTh->Clone();
  graphRatioPythia8Monash2013ToFit                           = CalculateGraphErrRatioToFit (graphRatioPythia8Monash2013ToFit, fitTCMInvYieldPi0);

  TGraph* graphRatioPi0CombNLOMuHalf               = (TGraph*)graphNLOCalcPi0MuHalf->Clone();cout << __LINE__ << endl;
  TGraph* graphRatioPi0CombNLOMuOne                = (TGraph*)graphNLOCalcPi0MuOne->Clone();cout << __LINE__ << endl;
  TGraph* graphRatioPi0CombNLOMuTwo                = (TGraph*)graphNLOCalcPi0MuTwo->Clone();cout << __LINE__ << endl;
  TGraphAsymmErrors* graphRatioPi0DSS14            = (TGraphAsymmErrors*)graphPi0DSS14->Clone();cout << __LINE__ << endl;

  graphRatioPi0CombNLOMuHalf                       = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuHalf, fitTCMInvYieldPi0); cout << __LINE__ << endl;
  graphRatioPi0CombNLOMuOne                        = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuOne, fitTCMInvYieldPi0); cout << __LINE__ << endl;
  graphRatioPi0CombNLOMuTwo                        = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuTwo, fitTCMInvYieldPi0); cout << __LINE__ << endl;
  graphRatioPi0DSS14                               = CalculateGraphErrRatioToFit (graphRatioPi0DSS14, fitTCMInvYieldPi0); cout << __LINE__ << endl;

  TH1D* histoRatioEtaHIJINGToFit                      = (TH1D*) histoHIJINGEta->Clone("histoRatioEtaHIJINGToFit");
  histoRatioEtaHIJINGToFit                            = CalculateHistoRatioToFit (histoRatioEtaHIJINGToFit, fitTCMInvYieldEta);
  // histoRatioEtaHIJINGToFit->GetXaxis()->SetRangeUser(0.6,20);

  TH1D* histoRatioEPOSJJToFitEta                  = (TH1D*) histoEPOSJJInvYieldEta->Clone();
  histoRatioEPOSJJToFitEta                        = CalculateHistoRatioToFit (histoRatioEPOSJJToFitEta, fitTCMInvYieldEta);
  histoRatioEPOSJJToFitEta->GetXaxis()->SetRangeUser(0.7,45);
  TGraphErrors* graphRatioEPOSJJToFitEta             = (TGraphErrors*) graphEPOSJJInvYieldEta->Clone();
  graphRatioEPOSJJToFitEta                           = CalculateGraphErrRatioToFit (graphRatioEPOSJJToFitEta, fitTCMInvYieldEta);
  while(graphRatioEPOSJJToFitEta->GetX()[0] < 0.7) graphRatioEPOSJJToFitEta->RemovePoint(0);

  TH1D* histoRatioDPMJETToFitEta                  = (TH1D*) histoDPMJETInvYieldEta->Clone();
  histoRatioDPMJETToFitEta                        = CalculateHistoRatioToFit (histoRatioDPMJETToFitEta, fitTCMInvYieldEta);
  histoRatioDPMJETToFitEta->GetXaxis()->SetRangeUser(0.7,45);
  TGraphErrors* graphRatioDPMJETToFitEta             = (TGraphErrors*) graphDPMJETInvYieldEta->Clone();
  graphRatioDPMJETToFitEta                           = CalculateGraphErrRatioToFit (graphRatioDPMJETToFitEta, fitTCMInvYieldEta);
  while(graphRatioDPMJETToFitEta->GetX()[0] < 0.7) graphRatioDPMJETToFitEta->RemovePoint(0);

  TH1D* histoRatioPythia8EPPS16ToFitEta                  = (TH1D*) histoPythia8EPPS16InvYieldEta->Clone();
  histoRatioPythia8EPPS16ToFitEta                        = CalculateHistoRatioToFit (histoRatioPythia8EPPS16ToFitEta, fitTCMInvYieldEta);
  histoRatioPythia8EPPS16ToFitEta->GetXaxis()->SetRangeUser(0.7,45);
  TGraphErrors* graphRatioPythia8EPPS16ToFitEta             = (TGraphErrors*) graphPythia8EPPS16InvYieldEta->Clone();
  graphRatioPythia8EPPS16ToFitEta                           = CalculateGraphErrRatioToFit (graphRatioPythia8EPPS16ToFitEta, fitTCMInvYieldEta);
  while(graphRatioPythia8EPPS16ToFitEta->GetX()[0] < 0.7) graphRatioPythia8EPPS16ToFitEta->RemovePoint(0);

  TH1D* histoRatioPythia8nCTEQ15ToFitEta                  = (TH1D*) histoPythia8nCTEQ15InvYieldEta->Clone();
  histoRatioPythia8nCTEQ15ToFitEta                        = CalculateHistoRatioToFit (histoRatioPythia8nCTEQ15ToFitEta, fitTCMInvYieldEta);
  histoRatioPythia8nCTEQ15ToFitEta->GetXaxis()->SetRangeUser(0.7,45);
  TGraphErrors* graphRatioPythia8nCTEQ15ToFitEta             = (TGraphErrors*) graphPythia8nCTEQ15InvYieldEta->Clone();
  graphRatioPythia8nCTEQ15ToFitEta                           = CalculateGraphErrRatioToFit (graphRatioPythia8nCTEQ15ToFitEta, fitTCMInvYieldEta);
  while(graphRatioPythia8nCTEQ15ToFitEta->GetX()[0] < 0.7) graphRatioPythia8nCTEQ15ToFitEta->RemovePoint(0);

  TH1D* histoRatioPythia8Monash2013ToFitEta                  = (TH1D*) histoPythia8Monash2013InvSecpThEta->Clone();
  histoRatioPythia8Monash2013ToFitEta                        = CalculateHistoRatioToFit (histoRatioPythia8Monash2013ToFitEta, fitTCMInvYieldEta);
  histoRatioPythia8Monash2013ToFitEta->GetXaxis()->SetRangeUser(0.7,45);
  TGraphErrors* graphRatioPythia8Monash2013ToFitEta             = (TGraphErrors*) graphPythia8Monash2013InvSecpThEta->Clone();
  graphRatioPythia8Monash2013ToFitEta                           = CalculateGraphErrRatioToFit (graphRatioPythia8Monash2013ToFitEta, fitTCMInvYieldEta);
  while(graphRatioPythia8Monash2013ToFitEta->GetX()[0] < 0.7) graphRatioPythia8Monash2013ToFitEta->RemovePoint(0);

  TGraph* graphRatioEtaCombNLOMuHalf                  = (TGraph*)graphNLOCalcEtaMuHalf->Clone();
  TGraph* graphRatioEtaCombNLOMuOne                   = (TGraph*)graphNLOCalcEtaMuOne->Clone();
  TGraph* graphRatioEtaCombNLOMuTwo                   = (TGraph*)graphNLOCalcEtaMuTwo->Clone();
  graphRatioEtaCombNLOMuHalf                          = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuHalf, fitTCMInvYieldEta);
  graphRatioEtaCombNLOMuOne                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuOne, fitTCMInvYieldEta);
  graphRatioEtaCombNLOMuTwo                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuTwo, fitTCMInvYieldEta);

  TGraphAsymmErrors* graphRatioEtaAESSS               = (TGraphAsymmErrors*) graphEtaAESSS->Clone();
  graphRatioEtaAESSS                                  = CalculateGraphErrRatioToFit (graphRatioEtaAESSS, fitTCMInvYieldEta);

  // *************************************************************************************************************

  TGraphAsymmErrors* graphRatioCombCombFitTot     = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone();
  graphRatioCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitTot, fitTCMInvYieldPi0);
  TGraphAsymmErrors* graphRatioCombCombFitStat    = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone();
  graphRatioCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioCombCombFitStat, fitTCMInvYieldPi0);
  TGraphAsymmErrors* graphRatioCombCombFitSys     = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone();
  graphRatioCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitSys, fitTCMInvYieldPi0);

  TGraphAsymmErrors* graphRatioCombCombFitTotEta     = (TGraphAsymmErrors*)graphCombEtaInvYieldTot->Clone();
  graphRatioCombCombFitTotEta                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitTotEta, fitTCMInvYieldEta);
  TGraphAsymmErrors* graphRatioCombCombFitStatEta    = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone();
  graphRatioCombCombFitStatEta                       = CalculateGraphErrRatioToFit(graphRatioCombCombFitStatEta, fitTCMInvYieldEta);
  TGraphAsymmErrors* graphRatioCombCombFitSysEta     = (TGraphAsymmErrors*)graphCombEtaInvYieldSys->Clone();
  graphRatioCombCombFitSysEta                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitSysEta, fitTCMInvYieldEta);

  TGraphAsymmErrors* graphRatioCombCombFitTotEtaToPi0     = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone();
  graphRatioCombCombFitTotEtaToPi0                        = CalculateFitRatioToGraphErr(fitTCMInvYieldEta, graphRatioCombCombFitTotEtaToPi0);
  TGraphAsymmErrors* graphRatioCombCombFitStatEtaToPi0    = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone();
  graphRatioCombCombFitStatEtaToPi0                       = CalculateFitRatioToGraphErr(fitTCMInvYieldEta, graphRatioCombCombFitStatEtaToPi0);
  TGraphAsymmErrors* graphRatioCombCombFitSysEtaToPi0     = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone();
  graphRatioCombCombFitSysEtaToPi0                        = CalculateFitRatioToGraphErr(fitTCMInvYieldEta, graphRatioCombCombFitSysEtaToPi0);

  TGraphAsymmErrors* graphRatioEtaCombCombPi0FititTotEtaToPi0     = (TGraphAsymmErrors*)graphCombEtaInvYieldTot->Clone();
  graphRatioEtaCombCombPi0FititTotEtaToPi0                        = CalculateGraphErrRatioToFit(graphRatioEtaCombCombPi0FititTotEtaToPi0, fitTCMInvYieldPi0);
  TGraphAsymmErrors* graphRatioEtaCombCombPi0FitStatEtaToPi0    = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone();
  graphRatioEtaCombCombPi0FitStatEtaToPi0                       = CalculateGraphErrRatioToFit(graphRatioEtaCombCombPi0FitStatEtaToPi0, fitTCMInvYieldPi0);
  TGraphAsymmErrors* graphRatioEtaCombCombPi0FitFitSysEtaToPi0     = (TGraphAsymmErrors*)graphCombEtaInvYieldSys->Clone();
  graphRatioEtaCombCombPi0FitFitSysEtaToPi0                        = CalculateGraphErrRatioToFit(graphRatioEtaCombCombPi0FitFitSysEtaToPi0, fitTCMInvYieldPi0);

  // if(i==0){
  while (graphPi0InvYieldStat[0]->GetX()[graphPi0InvYieldStat[0]->GetN()-1] > 12) graphPi0InvYieldStat[0]->RemovePoint(graphPi0InvYieldStat[0]->GetN()-1);  
  while (graphPi0InvYieldSys[0]->GetX()[graphPi0InvYieldSys[0]->GetN()-1] > 12) graphPi0InvYieldSys[0]->RemovePoint(graphPi0InvYieldSys[0]->GetN()-1);  
  // }
  
  TGraphAsymmErrors* graphRatioCombFitStat[13];
  TGraphAsymmErrors* graphRatioCombFitSys[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  graphRatioCombFitStat[i]                = (TGraphAsymmErrors*)graphPi0InvYieldStat[i]->Clone();
  graphRatioCombFitStat[i]                = CalculateGraphErrRatioToFit(graphRatioCombFitStat[i], fitTCMInvYieldPi0);
  graphRatioCombFitSys[i]                 = (TGraphAsymmErrors*)graphPi0InvYieldSys[i]->Clone();
  graphRatioCombFitSys[i]                 = CalculateGraphErrRatioToFit(graphRatioCombFitSys[i], fitTCMInvYieldPi0);
  }
  }
  TGraphAsymmErrors* graphRatioCombFitStatEta[13];
  TGraphAsymmErrors* graphRatioCombFitSysEta[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  graphRatioCombFitStatEta[i]                = (TGraphAsymmErrors*)graphEtaInvYieldStat[i]->Clone();
  graphRatioCombFitStatEta[i]                = CalculateGraphErrRatioToFit(graphRatioCombFitStatEta[i], fitTCMInvYieldEta);
  graphRatioCombFitSysEta[i]                 = (TGraphAsymmErrors*)graphEtaInvYieldSys[i]->Clone();
  graphRatioCombFitSysEta[i]                 = CalculateGraphErrRatioToFit(graphRatioCombFitSysEta[i], fitTCMInvYieldEta);
  }
  }

  TGraphAsymmErrors* graphRatioCombCombFitTot_WOXErr = (TGraphAsymmErrors*) graphRatioCombCombFitTot->Clone("graphRatioCombCombFitTot_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitTot_WOXErr);
  TGraphAsymmErrors* graphRatioCombCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioCombCombFitStat->Clone("graphRatioCombCombFitStat_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitStat_WOXErr);
  TGraphAsymmErrors* graphRatioCombCombFitSys_WOXErr = (TGraphAsymmErrors*) graphRatioCombCombFitSys->Clone("graphRatioCombCombFitSys_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitSys_WOXErr);

  TGraphAsymmErrors* graphRatioCombCombFitTotEta_WOXErr = (TGraphAsymmErrors*) graphRatioCombCombFitTotEta->Clone("graphRatioCombCombFitTotEta_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitTotEta_WOXErr);
  TGraphAsymmErrors* graphRatioCombCombFitStatEta_WOXErr = (TGraphAsymmErrors*) graphRatioCombCombFitStatEta->Clone("graphRatioCombCombFitStatEta_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitStatEta_WOXErr);
  TGraphAsymmErrors* graphRatioCombCombFitSysEta_WOXErr = (TGraphAsymmErrors*) graphRatioCombCombFitSysEta->Clone("graphRatioCombCombFitSysEta_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitSysEta_WOXErr);

  TGraphAsymmErrors* graphRatioCombFitStat_WOXErr[13];
  TGraphAsymmErrors* graphRatioCombFitSys_WOXErr[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  graphRatioCombFitStat_WOXErr[i] = (TGraphAsymmErrors*) graphRatioCombFitStat[i]->Clone(Form("graphRatioCombFitStat_WOXErr_%i",i));
  ProduceGraphAsymmWithoutXErrors(graphRatioCombFitStat_WOXErr[i]);
  graphRatioCombFitSys_WOXErr[i] = (TGraphAsymmErrors*) graphRatioCombFitSys[i]->Clone(Form("graphRatioCombCombFitStat_WOXErr_%i",i));
  ProduceGraphAsymmWithoutXErrors(graphRatioCombFitSys_WOXErr[i]);
  }
  }
  TGraphAsymmErrors* graphRatioCombFitStatEta_WOXErr[13];
  TGraphAsymmErrors* graphRatioCombFitSysEta_WOXErr[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  graphRatioCombFitStatEta_WOXErr[i] = (TGraphAsymmErrors*) graphRatioCombFitStatEta[i]->Clone(Form("graphRatioCombFitStatEta_WOXErr_%i",i));
  ProduceGraphAsymmWithoutXErrors(graphRatioCombFitStatEta_WOXErr[i]);
  graphRatioCombFitSysEta_WOXErr[i] = (TGraphAsymmErrors*) graphRatioCombFitSysEta[i]->Clone(Form("graphRatioCombCombFitStatEta_WOXErr_%i",i));
  ProduceGraphAsymmWithoutXErrors(graphRatioCombFitSysEta_WOXErr[i]);
  }
  }

  TGraphAsymmErrors* graphRatioPi0CombTsallisFitStat  = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone();
  graphRatioPi0CombTsallisFitStat                     = CalculateGraphErrRatioToFit(graphRatioPi0CombTsallisFitStat, fitInvXSectionPi0);
  TGraphAsymmErrors* graphRatioPi0CombHagedornFitStat = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone();
  graphRatioPi0CombHagedornFitStat                    = CalculateGraphErrRatioToFit(graphRatioPi0CombHagedornFitStat, fitOHagInvYieldPi0Tot);
  TGraphAsymmErrors* graphRatioPi0CombPowerFitStat    = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone();
  graphRatioPi0CombPowerFitStat                       = CalculateGraphErrRatioToFit(graphRatioPi0CombPowerFitStat, fitPowInvXSectionPi0);
  TGraphAsymmErrors* graphRatioPi0CombTsallisFitSys   = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone();
  graphRatioPi0CombTsallisFitSys                      = CalculateGraphErrRatioToFit(graphRatioPi0CombTsallisFitSys, fitInvXSectionPi0);
  TGraphAsymmErrors* graphRatioPi0CombHagedornFitSys  = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone();
  graphRatioPi0CombHagedornFitSys                     = CalculateGraphErrRatioToFit(graphRatioPi0CombHagedornFitSys, fitOHagInvYieldPi0Tot);
  TGraphAsymmErrors* graphRatioPi0CombPowerFitSys     = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone();
  graphRatioPi0CombPowerFitSys                        = CalculateGraphErrRatioToFit(graphRatioPi0CombPowerFitSys, fitPowInvXSectionPi0);

  TGraphAsymmErrors* graphRatioPi0CombTsallisFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombTsallisFitStat->Clone("graphRatioPi0CombTsallisFitStat_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombTsallisFitStat_WOXErr);
  TGraphAsymmErrors* graphRatioPi0CombHagedornFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombHagedornFitStat->Clone("graphRatioPi0CombHagedornFitStat_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombHagedornFitStat_WOXErr);
  TGraphAsymmErrors* graphRatioPi0CombPowerFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombPowerFitStat->Clone("graphRatioPi0CombPowerFitStat_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombPowerFitStat_WOXErr);

  TGraphAsymmErrors* graphRatioEtaCombTsallisFitStat  = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone();
  graphRatioEtaCombTsallisFitStat                     = CalculateGraphErrRatioToFit(graphRatioEtaCombTsallisFitStat, fitInvXSectionEta);
  TGraphAsymmErrors* graphRatioEtaCombHagedornFitStat = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone();
  graphRatioEtaCombHagedornFitStat                    = CalculateGraphErrRatioToFit(graphRatioEtaCombHagedornFitStat, fitOHagInvYieldEtaTot);
  TGraphAsymmErrors* graphRatioEtaCombPowerFitStat    = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone();
  graphRatioEtaCombPowerFitStat                       = CalculateGraphErrRatioToFit(graphRatioEtaCombPowerFitStat, fitPowInvXSectionEta);
  TGraphAsymmErrors* graphRatioEtaCombTsallisFitSys   = (TGraphAsymmErrors*)graphCombEtaInvYieldSys->Clone();
  graphRatioEtaCombTsallisFitSys                      = CalculateGraphErrRatioToFit(graphRatioEtaCombTsallisFitSys, fitInvXSectionEta);
  TGraphAsymmErrors* graphRatioEtaCombHagedornFitSys  = (TGraphAsymmErrors*)graphCombEtaInvYieldSys->Clone();
  graphRatioEtaCombHagedornFitSys                     = CalculateGraphErrRatioToFit(graphRatioEtaCombHagedornFitSys, fitOHagInvYieldEtaTot);
  TGraphAsymmErrors* graphRatioEtaCombPowerFitSys     = (TGraphAsymmErrors*)graphCombEtaInvYieldSys->Clone();
  graphRatioEtaCombPowerFitSys                        = CalculateGraphErrRatioToFit(graphRatioEtaCombPowerFitSys, fitPowInvXSectionEta);

  TGraphAsymmErrors* graphRatioEtaCombTsallisFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaCombTsallisFitStat->Clone("graphRatioEtaCombTsallisFitStat_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombTsallisFitStat_WOXErr);
  TGraphAsymmErrors* graphRatioEtaCombHagedornFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaCombHagedornFitStat->Clone("graphRatioEtaCombHagedornFitStat_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombHagedornFitStat_WOXErr);
  TGraphAsymmErrors* graphRatioEtaCombPowerFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaCombPowerFitStat->Clone("graphRatioEtaCombPowerFitStat_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombPowerFitStat_WOXErr);



  // **********************************************************************************************************************
  // ******************************************* Plot Ratio of Comb to Fit ****************************************
  // **********************************************************************************************************************
  textSizeLabelsPixel                 = 54;
  TCanvas* canvasRatioToCombFit       = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRatioToCombFit, 0.08, 0.01, 0.01, 0.125);
  canvasRatioToCombFit->SetLogx();

  Double_t textsizeLabelsPP       = 0;
  Double_t textsizeFacPP          = 0;
  if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
  textsizeLabelsPP            = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
  textsizeFacPP               = (Double_t)1./canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
  } else {
  textsizeLabelsPP            = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
  textsizeFacPP               = (Double_t)1./canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
  }
  cout << textsizeLabelsPP << endl;

  TH2F * histo2DPi0RatioToCombFit;
  histo2DPi0RatioToCombFit               = new TH2F("histo2DPi0RatioToCombFit","histo2DPi0RatioToCombFit",1000,minPtPi0, maxPtPi0,1000,0.2,3.0);
  SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/TCM fit", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                    0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
  histo2DPi0RatioToCombFit->GetXaxis()->SetMoreLogLabels();
  histo2DPi0RatioToCombFit->GetXaxis()->SetNoExponent(kTRUE);
  //  histo2DPi0RatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
  //  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.6,2.1);
  histo2DPi0RatioToCombFit->Draw("copy");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphRatioCombCombFitSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphRatioCombCombFitStat_WOXErr->Draw("p,same,z");

  DrawGammaLines(minPtPi0, maxPtPi0 , 1., 1.,0.1, kGray+2);
  DrawGammaLines(minPtPi0, maxPtPi0 , 1.1, 1.1,0.1, kGray, 7);
  DrawGammaLines(minPtPi0, maxPtPi0 , 0.9, 0.9,0.1, kGray, 7);

  // TLatex *labelRatioToFitEnergy   = new TLatex(0.72, 0.91, Form("ALICE, %s",collisionSystempPb8TeV.Data()));
  TLatex *labelRatioToFitEnergy   = new TLatex(0.955, 0.905, Form("%s",collisionSystempPb8TeV.Data()));
  SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel,4);
  labelRatioToFitEnergy->SetTextFont(43);
  labelRatioToFitEnergy->SetTextAlign(31);
  labelRatioToFitEnergy->Draw();
  TLatex *labelRatioToFitALICE    = new TLatex(0.955, 0.847, labelALICEPlot.Data());
  SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel,4);
  labelRatioToFitALICE->SetTextFont(43);
  labelRatioToFitALICE->SetTextAlign(31);
  labelRatioToFitALICE->Draw();
  // TLatex *labelRatioToFitPi0      = new TLatex(0.826, 0.807, "#pi^{0} #rightarrow #gamma#gamma");
  TLatex *labelRatioToFitPi0      = new TLatex(0.955, 0.795, "#pi^{0} #rightarrow #gamma#gamma");
  SetStyleTLatex( labelRatioToFitPi0, textSizeLabelsPixel,4);
  labelRatioToFitPi0->SetTextFont(43);
  labelRatioToFitPi0->SetTextAlign(31);
  labelRatioToFitPi0->Draw();

  canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_PP8TeV.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // *******************************************Plot different ratios to fits *********************************************
  // **********************************************************************************************************************

  histo2DPi0RatioToCombFit->SetYTitle("Data/fit");
  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
  histo2DPi0RatioToCombFit->Draw("copy");


  while (graphRatioPi0CombPowerFitSys->GetX()[0] < 3.0) graphRatioPi0CombPowerFitSys->RemovePoint(0);
  while (graphRatioPi0CombPowerFitStat_WOXErr->GetX()[0] < 3.0) graphRatioPi0CombPowerFitStat_WOXErr->RemovePoint(0);
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombPowerFitSys, 33, markerSizeComb*2, kAzure+2 , kAzure+2, widthLinesBoxes, kTRUE);
  graphRatioPi0CombPowerFitSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombPowerFitStat_WOXErr, 33, markerSizeComb*2, kAzure+2 , kAzure+2);
  graphRatioPi0CombPowerFitStat_WOXErr->Draw("p,same,z");

  DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombTsallisFitSys, 25, markerSizeComb, colorTrigg[1], colorTrigg[1], widthLinesBoxes, kTRUE);
  graphRatioPi0CombTsallisFitSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombTsallisFitStat_WOXErr, 25, markerSizeComb, colorTrigg[1] , colorTrigg[1]);
  graphRatioPi0CombTsallisFitStat_WOXErr->Draw("p,same,z");

  graphRatioCombCombFitSys->SetMarkerColor(kGray+2);
  graphRatioCombCombFitSys->Draw("E2same");
  graphRatioCombCombFitStat_WOXErr->SetMarkerColor(kGray+2);
  graphRatioCombCombFitStat_WOXErr->Draw("p,same,z");

  DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombHagedornFitSys, 24, markerSizeComb+0.2, colorTrigg[2], colorTrigg[2], widthLinesBoxes, kTRUE);
  graphRatioPi0CombHagedornFitSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombHagedornFitStat_WOXErr, 24, markerSizeComb+0.2, colorTrigg[2] , colorTrigg[2]);
  graphRatioPi0CombHagedornFitStat_WOXErr->Draw("p,same,z");

  DrawGammaLines(minPtPi0, maxPtPi0, 1., 1.,0.1, kGray+2);
  DrawGammaLines(minPtPi0, maxPtPi0 , 1.1, 1.1,0.1, kGray, 7);
  DrawGammaLines(minPtPi0, maxPtPi0 , 0.9, 0.9,0.1, kGray, 7);

  labelRatioToFitEnergy->Draw();
  labelRatioToFitALICE->Draw();
  labelRatioToFitPi0->Draw();

  TLegend* legendRatioPi0Fits= GetAndSetLegend2(0.12,0.95-4*1.05*textsizeLabelsPP,0.37,0.95, textSizeLabelsPixel, 1, "", 43, 0.2);
  legendRatioPi0Fits->AddEntry(graphRatioCombCombFitSys,"TCM","p");
  legendRatioPi0Fits->AddEntry(graphRatioPi0CombTsallisFitSys,"Levy-Tsallis","p");
  legendRatioPi0Fits->AddEntry(graphRatioPi0CombHagedornFitSys,"mod. Hagedorn","p");
  legendRatioPi0Fits->AddEntry(graphRatioPi0CombPowerFitSys,"pure powerlaw, 3.5-200 GeV/#it{c}","p");
  legendRatioPi0Fits->Draw();

  canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToDifferentFits_pPb8TeV.%s",outputDir.Data(),suffix.Data()));
  histo2DPi0RatioToCombFit->SetYTitle("Data/TCM fit");

  // **********************************************************************************************************************
  // *******************************************Plot Ratio of Individual meas to Fit ******************************************
  // **********************************************************************************************************************

  canvasRatioToCombFit->cd();
  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
  histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(minPtPi0,maxPtPi0);
  histo2DPi0RatioToCombFit->Draw("copy");
  DrawGammaLines(minPtPi0,maxPtPi0 , 1., 1.,1, kGray+2);
  DrawGammaLines(minPtPi0,maxPtPi0 , 1.1, 1.1,1, kGray, 7);
  DrawGammaLines(minPtPi0,maxPtPi0 , 0.9, 0.9,1, kGray, 7);
  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  DrawGammaSetMarkerTGraphAsym(graphRatioCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphRatioCombFitStat_WOXErr[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);

  graphRatioCombFitSys[i]->Draw("E2same");
  graphRatioCombFitStat_WOXErr[i]->Draw("p,same,z");
  }
  }



  labelRatioToFitEnergy->Draw();
  labelRatioToFitALICE->Draw();
  labelRatioToFitPi0->Draw();

  //****************************** Definition of the Legend ******************************************
  //**************** Row def ************************
  Double_t rowsLegendOnlyPi0Ratio[13]      = {0.92,0.88,0.84,0.80,0.76, 0.72,0.68,0.64,0.60,0.56,  0.52,0.48,0.44};
  Double_t rowsLegendOnlyPi0RatioAbs[13]   = {0.91,2.2,2.1,2.0,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1};
  Double_t columnsLegendOnlyPi0Ratio[3]   = {0.14,0.34, 0.40};
  Double_t columnsLegendOnlyPi0RatioAbs[3]= {0.14,pow(prodPtPi0,0.16), pow(prodPtPi0,0.31)};
  Double_t lengthBox                      = pow(prodPtPi0,0.35)-pow(prodPtPi0,0.33);
  Double_t heightBox                      = 0.08/2;
  //****************** first Column **************************************************
  TLatex *textSingleMeasRatioPi0[13];
  Int_t legendIndex = 0;
    for (Int_t i = 0; i < 13; i++){
    if(directoryPi0[i]){
      textSingleMeasRatioPi0[i]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[legendIndex+1],nameMeasGlobalPlot[i].Data());
      SetStyleTLatex( textSingleMeasRatioPi0[i], 0.85*textSizeLabelsPixel,4);
      textSingleMeasRatioPi0[i]->SetTextFont(43);
      textSingleMeasRatioPi0[i]->Draw();
      legendIndex++;
    }
  }

  //****************** second Column *************************************************
  TLatex *textStatOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[1]-0.01,rowsLegendOnlyPi0Ratio[0] ,"stat");
  SetStyleTLatex( textStatOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
  textStatOnlyRatioPi0->SetTextFont(43);
  textStatOnlyRatioPi0->Draw();
  TLatex *textSysOnlyRatioPi0             = new TLatex(columnsLegendOnlyPi0Ratio[2]+0.02 ,rowsLegendOnlyPi0Ratio[0],"syst");
  SetStyleTLatex( textSysOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
  textSysOnlyRatioPi0->SetTextFont(43);
  textSysOnlyRatioPi0->Draw();

  TMarker* markerPi0OnlyRatio[13];
  TBox* boxPi0OnlyRatio[13];
  legendIndex = 0;
  for (Int_t i = 0; i < 13; i++){
    if(directoryPi0[i]){
      markerPi0OnlyRatio[i]               = CreateMarkerFromGraph(graphRatioCombFitSys[i],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[legendIndex+1],1);
      markerPi0OnlyRatio[i]->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[legendIndex+1]);
      boxPi0OnlyRatio[i]                  = CreateBoxFromGraph(graphRatioCombFitSys[i], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[legendIndex+1]- heightBox,
      columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox+0.1, rowsLegendOnlyPi0RatioAbs[legendIndex+1]+ heightBox);
      boxPi0OnlyRatio[i]->Draw("l");
      legendIndex++;
    }
  }

  canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_pPb8TeV.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************************* Ratio of Comb to Fit ****************************************
  // **********************************************************************************************************************
  textSizeLabelsPixel                 = 48;
  canvasRatioToCombFit->cd();
  TH2F * histo2DEtaRatioToCombFit;
  histo2DEtaRatioToCombFit               = new TH2F("histo2DEtaRatioToCombFit","histo2DEtaRatioToCombFit",1000,minPtEta,maxPtEta,1000,0.,7.    );
  SetStyleHistoTH2ForGraphs(histo2DEtaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/TCM fit", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                    0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
  histo2DEtaRatioToCombFit->GetXaxis()->SetMoreLogLabels();
  histo2DEtaRatioToCombFit->GetXaxis()->SetNoExponent(kTRUE);
  //    histo2DEtaRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
  //  histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
  histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.3,1.8);
  histo2DEtaRatioToCombFit->Draw("copy");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysEta, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphRatioCombCombFitSysEta->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatEta_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphRatioCombCombFitStatEta_WOXErr->Draw("p,same,z");

  DrawGammaLines(minPtEta, maxPtEta , 1., 1.,1, kGray+2);
  DrawGammaLines(minPtEta, maxPtEta , 1.1, 1.1,1, kGray, 7);
  DrawGammaLines(minPtEta, maxPtEta , 0.9, 0.9,1, kGray, 7);

  TLatex *labelRatioToFitEnergy2   = new TLatex(0.955, 0.905, Form("%s",collisionSystempPb8TeV.Data()));
  SetStyleTLatex( labelRatioToFitEnergy2, textSizeLabelsPixel,4);
  labelRatioToFitEnergy2->SetTextFont(43);
  labelRatioToFitEnergy2->SetTextAlign(31);
  labelRatioToFitEnergy2->Draw();
  TLatex *labelRatioToFitALICE2    = new TLatex(0.955, 0.847, labelALICEPlot.Data());
  SetStyleTLatex( labelRatioToFitALICE2, textSizeLabelsPixel,4);
  labelRatioToFitALICE2->SetTextFont(43);
  labelRatioToFitALICE2->SetTextAlign(31);
  labelRatioToFitALICE2->Draw();

  // TLatex *labelRatioToFitEta      = new TLatex(0.84,0.82,"#eta #rightarrow #gamma#gamma");
  TLatex *labelRatioToFitEta      = new TLatex(0.955,0.807,"#eta #rightarrow #gamma#gamma");
  SetStyleTLatex( labelRatioToFitEta, textSizeLabelsPixel,4);
  labelRatioToFitEta->SetTextFont(43);
  labelRatioToFitEta->SetTextAlign(31);
  labelRatioToFitEta->Draw();

  canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_pPb8TeV.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // *******************************************Plot different ratios to fits *********************************************
  // **********************************************************************************************************************

  histo2DEtaRatioToCombFit->SetYTitle("Data/fit");
  histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.1,2.45);
  histo2DEtaRatioToCombFit->Draw("copy");

  while (graphRatioEtaCombPowerFitSys->GetX()[0] < 3.0) graphRatioEtaCombPowerFitSys->RemovePoint(0);
  while (graphRatioEtaCombPowerFitStat_WOXErr->GetX()[0] < 3.0) graphRatioEtaCombPowerFitStat_WOXErr->RemovePoint(0);

  DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombPowerFitSys, 33, markerSizeComb*2, kAzure+2 , kAzure+2, widthLinesBoxes, kTRUE);
  graphRatioEtaCombPowerFitSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombPowerFitStat_WOXErr, 33, markerSizeComb*2, kAzure+2 , kAzure+2);
  graphRatioEtaCombPowerFitStat_WOXErr->Draw("p,same,z");

  DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombTsallisFitSys, 25, markerSizeComb, colorTrigg[1], colorTrigg[1], widthLinesBoxes, kTRUE);
  graphRatioEtaCombTsallisFitSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombTsallisFitStat_WOXErr, 25, markerSizeComb, colorTrigg[1] , colorTrigg[1]);
  graphRatioEtaCombTsallisFitStat_WOXErr->Draw("p,same,z");

  graphRatioCombCombFitSysEta->SetMarkerColor(kGray+2);
  graphRatioCombCombFitSysEta->Draw("E2same");
  graphRatioCombCombFitStatEta_WOXErr->SetMarkerColor(kGray+2);
  graphRatioCombCombFitStatEta_WOXErr->Draw("p,same,z");

  DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombHagedornFitSys, 24, markerSizeComb+0.2, colorTrigg[2], colorTrigg[2], widthLinesBoxes, kTRUE);
  graphRatioEtaCombHagedornFitSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombHagedornFitStat_WOXErr, 24, markerSizeComb+0.2, colorTrigg[2] , colorTrigg[2]);
  graphRatioEtaCombHagedornFitStat_WOXErr->Draw("p,same,z");

  DrawGammaLines(minPtEta, maxPtEta , 1., 1.,0.1, kGray+2);
  DrawGammaLines(minPtEta, maxPtEta , 1.1, 1.1,0.1, kGray, 7);
  DrawGammaLines(minPtEta, maxPtEta , 0.9, 0.9,0.1, kGray, 7);

  labelRatioToFitEnergy2->Draw();
  labelRatioToFitALICE2->Draw();
  labelRatioToFitEta->Draw();

  TLegend* legendRatioEtaFits= GetAndSetLegend2(0.12,0.95-4*1.05*textsizeLabelsPP,0.37,0.95, textSizeLabelsPixel, 1, "", 43, 0.2);
  legendRatioEtaFits->AddEntry(graphRatioCombCombFitSysEta,"TCM","p");
  legendRatioEtaFits->AddEntry(graphRatioEtaCombTsallisFitSys,"Levy-Tsallis","p");
  legendRatioEtaFits->AddEntry(graphRatioEtaCombHagedornFitSys,"mod. Hagedorn","p");
  legendRatioEtaFits->AddEntry(graphRatioEtaCombPowerFitSys,"pure powerlaw, 3.5-18 GeV/#it{c}","p");
  legendRatioEtaFits->Draw();

  canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToDifferentFits_pPb8TeV.%s",outputDir.Data(),suffix.Data()));
  histo2DEtaRatioToCombFit->SetYTitle("Data/TCM fit");

  // **********************************************************************************************************************
  // ******************************************* Ratio of Individual meas to Fit ******************************************
  // **********************************************************************************************************************
  histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
  // canvasRatioToCombFit->cd();
  histo2DEtaRatioToCombFit->Draw("copy");

  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  DrawGammaSetMarkerTGraphAsym(graphRatioCombFitSysEta[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphRatioCombFitStatEta_WOXErr[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);

  graphRatioCombFitSysEta[i]->Draw("E2same");
  graphRatioCombFitStatEta_WOXErr[i]->Draw("p,same,z");
  }
  }

  DrawGammaLines(minPtEta,maxPtEta , 1., 1.,1, kGray+2);
  DrawGammaLines(minPtEta,maxPtEta , 1.1, 1.1,1, kGray, 7);
  DrawGammaLines(minPtEta,maxPtEta , 0.9, 0.9,1, kGray, 7);

  labelRatioToFitEnergy->Draw();
  labelRatioToFitALICE->Draw();
  labelRatioToFitEta->Draw();

  //****************************** Definition of the Legend ******************************************
  //**************** Row def ************************
  Double_t rowsLegendOnlyEtaRatio[11]      = {0.92,0.88,0.84,0.80,0.76,0.72,0.68,0.64,0.60};
  Double_t rowsLegendOnlyEtaRatioAbs[11]   = {0.91,2.18,2.07,1.96,1.85,1.74,1.63,1.52,1.41,1.3};
  Double_t columnsLegendOnlyEtaRatio[3]   = {0.15,0.34, 0.40};
  Double_t columnsLegendOnlyEtaRatioAbs[3]= {0.15,pow(prodPtEta,0.10), pow(prodPtEta,0.335)};
  Double_t lengthBoxEta                   = pow(prodPtEta,0.28)-pow(prodPtEta,0.32);
  Int_t rowcounter = 0;
  //****************** first Column **************************************************
  TLatex *textSingleMeasRatioEta[13];
  for (Int_t i = 0; i < 13; i++){
    if(directoryEta[i]&&graphEtaInvYieldSys[i]){
      textSingleMeasRatioEta[i]           = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[rowcounter+1],nameMeasGlobalPlot[i].Data());
      SetStyleTLatex( textSingleMeasRatioEta[i], 0.85*textSizeLabelsPixel,4);
      textSingleMeasRatioEta[i]->SetTextFont(43);
      textSingleMeasRatioEta[i]->Draw();
      rowcounter++;
    }
  }

  //****************** second Column *************************************************
  TLatex *textStatOnlyRatioEta            = new TLatex(columnsLegendOnlyEtaRatio[1]-0.02,rowsLegendOnlyEtaRatio[0] ,"stat");
  SetStyleTLatex( textStatOnlyRatioEta, 0.85*textSizeLabelsPixel,4);
  textStatOnlyRatioEta->SetTextFont(43);
  textStatOnlyRatioEta->Draw();
  TLatex *textSysOnlyRatioEta             = new TLatex(columnsLegendOnlyEtaRatio[2]+0.00 ,rowsLegendOnlyEtaRatio[0],"syst");
  SetStyleTLatex( textSysOnlyRatioEta, 0.85*textSizeLabelsPixel,4);
  textSysOnlyRatioEta->SetTextFont(43);
  textSysOnlyRatioEta->Draw();

  TMarker* markerEtaOnlyRatio[13];
  TBox* boxEtaOnlyRatio[13];
  rowcounter = 0;
  for (Int_t i = 0; i < 13; i++){
    if(directoryEta[i]&&graphEtaInvYieldSys[i]){
      markerEtaOnlyRatio[i]               = CreateMarkerFromGraph(graphRatioCombFitSysEta[i],columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[rowcounter+1],1);
      markerEtaOnlyRatio[i]->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[rowcounter+1]);
      boxEtaOnlyRatio[i]                  = CreateBoxFromGraph(graphRatioCombFitSysEta[i], columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBoxEta , rowsLegendOnlyEtaRatioAbs[rowcounter+1]- heightBox,
      columnsLegendOnlyEtaRatioAbs[2]+ 3*lengthBoxEta+0.1, rowsLegendOnlyEtaRatioAbs[rowcounter+1]+ heightBox);
      boxEtaOnlyRatio[i]->Draw("l");
      rowcounter++;
    }
  }

  canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_pPb8TeV.%s",outputDir.Data(),suffix.Data()));

  // *************************************************************************************************************
  // Calculate ratios to combined TSALLIS fit
  // *************************************************************************************************************

  TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitTot     = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone();
  graphRatioPi0TsallisCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioPi0TsallisCombCombFitTot, fitInvXSectionPi0);
  TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitStat    = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone();
  graphRatioPi0TsallisCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioPi0TsallisCombCombFitStat, fitInvXSectionPi0);
  TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitSys     = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone();
  graphRatioPi0TsallisCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioPi0TsallisCombCombFitSys, fitInvXSectionPi0);

  TGraphAsymmErrors* graphRatioPi0TsallisCombFitSys[13];
  TGraphAsymmErrors* graphRatioPi0TsallisCombFitStat[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  graphRatioPi0TsallisCombFitSys[i]                    = (TGraphAsymmErrors*)graphPi0InvYieldSys[i]->Clone();
  graphRatioPi0TsallisCombFitSys[i]                    = CalculateGraphErrRatioToFit(graphRatioPi0TsallisCombFitSys[i], fitInvXSectionPi0);
  graphRatioPi0TsallisCombFitStat[i]                   = (TGraphAsymmErrors*)graphPi0InvYieldStat[i]->Clone();
  graphRatioPi0TsallisCombFitStat[i]                   = CalculateGraphErrRatioToFit(graphRatioPi0TsallisCombFitStat[i], fitInvXSectionPi0);
  }
  }

  TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitSys_WOXErr = (TGraphAsymmErrors*) graphRatioPi0TsallisCombCombFitSys->Clone("graphRatioPi0TsallisCombCombFitSys_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioPi0TsallisCombCombFitSys_WOXErr);
  TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0TsallisCombCombFitStat->Clone("graphRatioPi0TsallisCombCombFitStat_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioPi0TsallisCombCombFitStat_WOXErr);
  TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitTot_WOXErr = (TGraphAsymmErrors*) graphRatioPi0TsallisCombCombFitTot->Clone("graphRatioPi0TsallisCombCombFitTot_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioPi0TsallisCombCombFitTot_WOXErr);

  TGraphAsymmErrors* graphRatioPi0TsallisCombFitSys_WOXErr[13];
  TGraphAsymmErrors* graphRatioPi0TsallisCombFitStat_WOXErr[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  graphRatioPi0TsallisCombFitSys_WOXErr[i]                    = (TGraphAsymmErrors*)graphRatioPi0TsallisCombFitSys[i]->Clone();
  ProduceGraphAsymmWithoutXErrors(graphRatioPi0TsallisCombFitSys_WOXErr[i]);
  graphRatioPi0TsallisCombFitStat_WOXErr[i]                   = (TGraphAsymmErrors*)graphRatioPi0TsallisCombFitStat[i]->Clone();
  ProduceGraphAsymmWithoutXErrors(graphRatioPi0TsallisCombFitStat_WOXErr[i]);
  }
  }

  TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitTot     = (TGraphAsymmErrors*)graphCombEtaInvYieldTot->Clone();
  graphRatioEtaTsallisCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombCombFitTot, fitInvXSectionEta);
  TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitStat    = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone();
  graphRatioEtaTsallisCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombCombFitStat, fitInvXSectionEta);
  TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitSys     = (TGraphAsymmErrors*)graphCombEtaInvYieldSys->Clone();
  graphRatioEtaTsallisCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombCombFitSys, fitInvXSectionEta);

  TGraphAsymmErrors* graphRatioEtaTsallisCombFitSys[13];
  TGraphAsymmErrors* graphRatioEtaTsallisCombFitStat[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  graphRatioEtaTsallisCombFitSys[i]                    = (TGraphAsymmErrors*)graphEtaInvYieldSys[i]->Clone();
  graphRatioEtaTsallisCombFitSys[i]                    = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombFitSys[i], fitInvXSectionEta);
  graphRatioEtaTsallisCombFitStat[i]                   = (TGraphAsymmErrors*)graphEtaInvYieldStat[i]->Clone();
  graphRatioEtaTsallisCombFitStat[i]                   = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombFitStat[i], fitInvXSectionEta);
  }
  }

  TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitSys_WOXErr = (TGraphAsymmErrors*) graphRatioEtaTsallisCombCombFitSys->Clone("graphRatioEtaTsallisCombCombFitSys_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisCombCombFitSys_WOXErr);
  TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaTsallisCombCombFitStat->Clone("graphRatioEtaTsallisCombCombFitStat_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisCombCombFitStat_WOXErr);
  TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitTot_WOXErr = (TGraphAsymmErrors*) graphRatioEtaTsallisCombCombFitTot->Clone("graphRatioEtaTsallisCombCombFitTot_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisCombCombFitTot_WOXErr);

  TGraphAsymmErrors* graphRatioEtaTsallisCombFitSys_WOXErr[13];
  TGraphAsymmErrors* graphRatioEtaTsallisCombFitStat_WOXErr[13];
  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  graphRatioEtaTsallisCombFitSys_WOXErr[i]                    = (TGraphAsymmErrors*)graphRatioEtaTsallisCombFitSys[i]->Clone();
  ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisCombFitSys_WOXErr[i]);
  graphRatioEtaTsallisCombFitStat_WOXErr[i]                   = (TGraphAsymmErrors*)graphRatioEtaTsallisCombFitStat[i]->Clone();
  ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisCombFitStat_WOXErr[i]);
  }
  }

  // **********************************************************************************************************************
  // ******************************************* Plot Ratio of Comb to Fit ****************************************
  // **********************************************************************************************************************

  canvasRatioToCombFit->SetLogx();
  // histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(0.23,50);
  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.6,2.1);
  histo2DPi0RatioToCombFit->SetYTitle("Data/Tsallis fit");
  histo2DPi0RatioToCombFit->Draw("copy");

  DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisCombCombFitSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphRatioPi0TsallisCombCombFitSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphRatioPi0TsallisCombCombFitStat_WOXErr->Draw("p,same,z");

  DrawGammaLines(minPtPi0, maxPtPi0 , 1., 1.,0.1, kGray+2);
  DrawGammaLines(minPtPi0, maxPtPi0 , 1.1, 1.1,0.1, kGray, 7);
  DrawGammaLines(minPtPi0, maxPtPi0 , 0.9, 0.9,0.1, kGray, 7);

  labelRatioToFitEnergy->Draw();
  // labelRatioToFitALICE->Draw();
  labelRatioToFitPi0->Draw();

  canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombTsallisFit_pPb8TeV.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // *******************************************Plot Ratio of Individual meas to Fit ******************************************
  // **********************************************************************************************************************

  canvasRatioToCombFit->cd();
  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
  histo2DPi0RatioToCombFit->Draw("copy");

  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisCombFitStat_WOXErr[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);

  graphRatioPi0TsallisCombFitSys[i]->Draw("E2same");
  graphRatioPi0TsallisCombFitStat_WOXErr[i]->Draw("p,same,z");
  }
  }

  DrawGammaLines(minPtPi0,maxPtPi0 , 1., 1.,1, kGray+2);
  DrawGammaLines(minPtPi0,maxPtPi0 , 1.1, 1.1,1, kGray, 7);
  DrawGammaLines(minPtPi0,maxPtPi0 , 0.9, 0.9,1, kGray, 7);

  labelRatioToFitEnergy->Draw();
  // labelRatioToFitALICE->Draw();
  labelRatioToFitPi0->Draw();

  //****************************** Definition of the Legend ******************************************
  //**************** Row def ************************

  //****************** first Column **************************************************
  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  textSingleMeasRatioPi0[i]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[i+1],nameMeasGlobalPlot[i].Data());
  SetStyleTLatex( textSingleMeasRatioPi0[i], 0.85*textSizeLabelsPixel,4);
  textSingleMeasRatioPi0[i]->SetTextFont(43);
  textSingleMeasRatioPi0[i]->Draw();
  }
  }

  //****************** second Column *************************************************
  textStatOnlyRatioPi0->Draw();
  textSysOnlyRatioPi0->Draw();

  for (Int_t i = 0; i < 13; i++){
  if(directoryPi0[i]){
  markerPi0OnlyRatio[i]               = CreateMarkerFromGraph(graphRatioPi0TsallisCombFitSys[i],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[i+1],1);
  markerPi0OnlyRatio[i]->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[i+1]);
  boxPi0OnlyRatio[i]                  = CreateBoxFromGraph(graphRatioPi0TsallisCombFitSys[i], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[i+1]- heightBox,
  columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox+0.1, rowsLegendOnlyPi0RatioAbs[i+1]+ heightBox);
  boxPi0OnlyRatio[i]->Draw("l");
  }
  }

  canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToTsallisFit_pPb8TeV.%s",outputDir.Data(),suffix.Data()));
  histo2DPi0RatioToCombFit->SetYTitle("Data/TCM fit");

  // **********************************************************************************************************************
  // ******************************************* Ratio of Comb to Fit ****************************************
  // **********************************************************************************************************************
  canvasRatioToCombFit->cd();
  histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.3,1.8);
  histo2DEtaRatioToCombFit->SetYTitle("Data/Tsallis fit");
  histo2DEtaRatioToCombFit->Draw("copy");

  DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisCombCombFitSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphRatioEtaTsallisCombCombFitSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphRatioEtaTsallisCombCombFitStat_WOXErr->Draw("p,same,z");

  DrawGammaLines(minPtEta,maxPtEta , 1., 1.,1, kGray+2);
  DrawGammaLines(minPtEta,maxPtEta , 1.1, 1.1,1, kGray, 7);
  DrawGammaLines(minPtEta,maxPtEta , 0.9, 0.9,1, kGray, 7);

  labelRatioToFitEnergy2->Draw();
  // labelRatioToFitALICE2->Draw();
  labelRatioToFitEta->Draw();

  canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombTsallisFit_pPb8TeV.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************************* Ratio of Individual meas to Fit ******************************************
  // **********************************************************************************************************************

  canvasRatioToCombFit->cd();
  histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);

  histo2DEtaRatioToCombFit->Draw("copy");

  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisCombFitStat_WOXErr[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);

  graphRatioEtaTsallisCombFitSys[i]->Draw("E2same");
  graphRatioEtaTsallisCombFitStat_WOXErr[i]->Draw("p,same,e");
  }
  }

  DrawGammaLines(minPtEta,maxPtEta , 1., 1.,1, kGray+2);
  DrawGammaLines(minPtEta,maxPtEta , 1.1, 1.1,1, kGray, 7);
  DrawGammaLines(minPtEta,maxPtEta , 0.9, 0.9,1, kGray, 7);

  labelRatioToFitEnergy->Draw();
  // labelRatioToFitALICE2->Draw();
  labelRatioToFitEta->Draw();

  //****************************** Definition of the Legend ******************************************
  //**************** Row def ************************

  //****************** first Column **************************************************
  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  textSingleMeasRatioEta[i]           = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[i+1],nameMeasGlobalPlot[i].Data());
  SetStyleTLatex( textSingleMeasRatioEta[i], 0.85*textSizeLabelsPixel,4);
  textSingleMeasRatioEta[i]->SetTextFont(43);
  textSingleMeasRatioEta[i]->Draw();
  }
  }

  //****************** second Column *************************************************
  textStatOnlyRatioEta->Draw();
  textSysOnlyRatioEta->Draw();

  for (Int_t i = 0; i < 13; i++){
  if(directoryEta[i]){
  markerEtaOnlyRatio[i]               = CreateMarkerFromGraph(graphRatioEtaTsallisCombFitSys[i],columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[i+1],1);
  markerEtaOnlyRatio[i]->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[i+1]);
  boxEtaOnlyRatio[i]                  = CreateBoxFromGraph(graphRatioEtaTsallisCombFitSys[i], columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBoxEta , rowsLegendOnlyEtaRatioAbs[i+1]- heightBox,
  columnsLegendOnlyEtaRatioAbs[2]+ 3*lengthBoxEta+0.1, rowsLegendOnlyEtaRatioAbs[i+1]+ heightBox);
  boxEtaOnlyRatio[i]->Draw("l");
  }
  }

  canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToTsallisFit_pPb8TeV.%s",outputDir.Data(),suffix.Data()));





  // **********************************************************************************************************************
  // ******************************************* Mass and width for pi0 at 8TeV ****************************************
  // **********************************************************************************************************************

  Double_t arrayBoundariesX1_4[2];
  Double_t arrayBoundariesY1_4[3];
  Double_t relativeMarginsX[3];
  Double_t relativeMarginsY[3];
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

  TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, minPtPi0,maxPtPi0woMerged ,1000., -30, 40);
  SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
  histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,25.5);//24.5);
  histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
  histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
  histo2DAllPi0FWHM->DrawCopy();

  for (Int_t i = 0; i < 13; i++){
      if(graphPi0FWHMMeV[i] && graphPi0TrueFWHMMeV[i]){
          DrawGammaSetMarkerTGraphAsym(graphPi0FWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
          graphPi0FWHMMeV[i]->Draw("p,same,e");
          DrawGammaSetMarkerTGraphAsym(graphPi0TrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
          graphPi0TrueFWHMMeV[i]->Draw("p,same,e");
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
  TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystempPb8TeV.Data());
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

  TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, minPtPi0,maxPtPi0woMerged, 1000., 100.1, 160.9);//125.1, 155.9);
  SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
  histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
  histo2DAllPi0Mass->GetYaxis()->SetRangeUser(131.1, 147.9);//125.1, 155.9);
  histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
  histo2DAllPi0Mass->GetXaxis()->SetNoExponent();
  histo2DAllPi0Mass->DrawCopy();

  for (Int_t i = 0; i < 13; i++){
      if(graphPi0Mass[i] && graphPi0TrueMass[i]){
          DrawGammaSetMarkerTGraphAsym(graphPi0Mass[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
          graphPi0Mass[i]->Draw("p,same,e");
          DrawGammaSetMarkerTGraphAsym(graphPi0TrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
          graphPi0TrueMass[i]->Draw("p,same,e");
      }
  }

  DrawGammaLines(minPtPi0,maxPtPi0woMerged , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1, kGray);

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
  Int_t legendRunningIndex = 0;
  for (Int_t i = 0; i < 13; i++){
      if(graphPi0Mass[i] && graphPi0TrueMass[i] && graphPi0FWHMMeV[i] && graphPi0TrueFWHMMeV[i]){
          textMassPCM[i]                  = new TLatex(columnsLegendMass2[0],rowsLegendMass2[legendRunningIndex+1],nameMeasGlobalPlot[i].Data());
          SetStyleTLatex( textMassPCM[i], textSizeLabelsPixel,4);
          textMassPCM[i]->SetTextFont(43);
          textMassPCM[i]->Draw();
          legendRunningIndex++;
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
  legendRunningIndex = 0;
  for (Int_t i = 0; i < 13; i++){
      if(graphPi0Mass[i] && graphPi0TrueMass[i]){
          markerPCMPi0Mass[i]             = CreateMarkerFromGraph(graphPi0Mass[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
          markerPCMPi0Mass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2);
          markerPCMPi0MassMC[i]           = CreateMarkerFromGraph(graphPi0TrueMass[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
          markerPCMPi0MassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[legendRunningIndex+1]+ offsetMarkerYMass2);
          legendRunningIndex++;
      }
  }

  canvasMassWidthPi0->Update();
  canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************************* Mass and width for Eta at 8TeV ****************************************
  // **********************************************************************************************************************


  textSizeLabelsPixel             = 50;
  ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

  TCanvas* canvasMassWidthEta     = new TCanvas("canvasMassWidthEta","",0,0,1350,1250);  // gives the page size
  DrawGammaCanvasSettings( canvasMassWidthEta,  0.13, 0.02, 0.03, 0.06);

  TPad* padWidthEta               = new TPad("padWidthEta", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
  DrawGammaPadSettings( padWidthEta, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
  padWidthEta->Draw();

  TPad* padMassEta                = new TPad("padMassEta", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
  DrawGammaPadSettings( padMassEta, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
  padMassEta->Draw();

  TPad* padMassLegend1Eta            = new TPad("padMassLegend1Eta", "", 0.13, 0.34, 0.52, 0.52,-1, -1, -2);
  DrawGammaPadSettings( padMassLegend1Eta, 0., 0., 0., 0.);
  padMassLegend1Eta->SetFillStyle(0);
  padMassLegend1Eta->Draw();

  padWidthEta->cd();
  padWidthEta->SetLogx();

  margin                 = relativeMarginsX[0]*2.7*1350;
  textsizeLabelsWidth    = 0;
  textsizeFacWidth       = 0;
  if (padWidthEta->XtoPixel(padWidthEta->GetX2()) < padWidthEta->YtoPixel(padWidthEta->GetY1())){
      textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthEta->XtoPixel(padWidthEta->GetX2()) ;
      textsizeFacWidth            = (Double_t)1./padWidthEta->XtoPixel(padWidthEta->GetX2()) ;
  } else {
      textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthEta->YtoPixel(padWidthEta->GetY1());
      textsizeFacWidth            = (Double_t)1./padWidthEta->YtoPixel(padWidthEta->GetY1());
  }

  TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, minPtEta,maxPtEta ,1000., -4, 105.5);
  SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
  histo2DAllEtaFWHM->GetYaxis()->SetRangeUser(-4.,59.5); //69.5
  histo2DAllEtaFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllEtaFWHM->GetYaxis()->SetNdivisions(505);
  histo2DAllEtaFWHM->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllEtaFWHM->GetXaxis()->SetTickLength(0.05);
  histo2DAllEtaFWHM->GetYaxis()->SetTickLength(0.026);
  histo2DAllEtaFWHM->DrawCopy();

  for (Int_t i = 0; i <11; i++){
      if(graphEtaFWHMMeV[i] && graphEtaTrueFWHMMeV[i]){
          DrawGammaSetMarkerTGraphAsym(graphEtaFWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
          graphEtaFWHMMeV[i]->Draw("p,same,e");
          DrawGammaSetMarkerTGraphAsym(graphEtaTrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
          graphEtaTrueFWHMMeV[i]->Draw("p,same,e");
      }
  }

  labelLegendAMass->Draw();

  labelMassPerf->Draw();
  labelMassEnergy->Draw();
  TLatex *labelMassEta        = new TLatex(0.13,0.69,"#eta #rightarrow #gamma#gamma");
  SetStyleTLatex( labelMassEta, textSizeLabelsPixel,4);
  labelMassEta->SetTextFont(43);
  labelMassEta->Draw();

  padMassEta->cd();
  padMassEta->SetLogx();

  Double_t textsizeLabelsMassEta         = 0;
  Double_t textsizeFacMassEta            = 0;
  if (padMassEta->XtoPixel(padMassEta->GetX2()) <padMassEta->YtoPixel(padMassEta->GetY1()) ){
      textsizeLabelsMassEta              = (Double_t)textSizeLabelsPixel/padMassEta->XtoPixel(padMassEta->GetX2()) ;
      textsizeFacMassEta                 = (Double_t)1./padMassEta->XtoPixel(padMassEta->GetX2()) ;
  } else {
      textsizeLabelsMassEta              = (Double_t)textSizeLabelsPixel/padMassEta->YtoPixel(padMassEta->GetY1());
      textsizeFacMassEta                 = (Double_t)1./padMassEta->YtoPixel(padMassEta->GetY1());
  }

  TH2F * histo2DAllEtaMass            = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",20, minPtEta,maxPtEta, 1000., 521.1, 579.9); //485.1, 639.9);
  SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMassEta, textsizeLabelsMassEta, 0.85*textsizeLabelsMassEta,
                            textsizeLabelsMassEta, 0.9, 0.28/(textsizeFacMassEta*margin), 512, 505);
  histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAllEtaMass->GetYaxis()->SetNdivisions(505);
  histo2DAllEtaMass->GetYaxis()->SetNoExponent(kTRUE);
  histo2DAllEtaMass->GetXaxis()->SetTickLength(0.05);
  histo2DAllEtaMass->GetXaxis()->SetNoExponent();
  histo2DAllEtaMass->DrawCopy();

  for (Int_t i = 0; i < 13; i++){
      if(graphEtaMass[i] && graphEtaTrueMass[i]){
          DrawGammaSetMarkerTGraphAsym(graphEtaMass[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
          graphEtaMass[i]->Draw("p,same,e");
          DrawGammaSetMarkerTGraphAsym(graphEtaTrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
          graphEtaTrueMass[i]->Draw("p,same,e");
      }
  }

  DrawGammaLines(minPtEta,maxPtEta , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.1, kGray);

  labelLegendBMass->Draw();

  //********************************** Defintion of the Legend **************************************************
  Double_t columnsLegendMass2Eta[3]      = {0.,0.57,0.84};
  Double_t rowsLegendMass2Eta[5]         = {0.8,0.6,0.4,0.2,0.01};
  //******************* Offsets ***********************
  Double_t offsetMarkerXMass2Eta         = 0.1;
  Double_t offsetMarkerYMass2Eta         = 0.1;
  //****************** Scale factors ******************
  Double_t scaleMarkerMass2Eta           = 1.2;

  padMassLegend1Eta->cd();
  //****************** first Column **************************************************
  TLatex *textMassPCMEta[10];
  Int_t counterEta = 1;
  for (Int_t i = 0; i < 13; i++){
      if(graphEtaMass[i] && graphEtaTrueMass[i] && graphEtaFWHMMeV[i] && graphEtaTrueFWHMMeV[i]){
          textMassPCMEta[i]                  = new TLatex(columnsLegendMass2Eta[0],rowsLegendMass2Eta[counterEta],nameMeasGlobalPlot[i].Data());
          SetStyleTLatex( textMassPCMEta[i], textSizeLabelsPixel,4);
          textMassPCMEta[i]->SetTextFont(43);
          textMassPCMEta[i]->Draw();
          counterEta+=1;
      }
  }
  //****************** second Column *************************************************
  TLatex *textMassDataEta                = new TLatex(columnsLegendMass2Eta[1],rowsLegendMass2Eta[0] ,"Data");
  SetStyleTLatex( textMassDataEta, textSizeLabelsPixel,4);
  textMassDataEta->SetTextFont(43);
  textMassDataEta->Draw();
  TLatex *textMassMCEta                  = new TLatex(columnsLegendMass2Eta[2] ,rowsLegendMass2Eta[0],"MC");
  SetStyleTLatex( textMassMCEta, textSizeLabelsPixel,4);
  textMassMCEta->SetTextFont(43);
  textMassMCEta->Draw();

  TMarker* markerPCMEtaMass[10];
  TMarker* markerPCMEtaMassMC[10];
  counterEta = 1;
  for (Int_t i = 0; i < 13; i++){
      if(graphEtaMass[i] && graphEtaTrueMass[i]){
          markerPCMEtaMass[i]             = CreateMarkerFromGraph(graphEtaMass[i],columnsLegendMass2Eta[1]+ offsetMarkerXMass2Eta ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta ,scaleMarkerMass2Eta);
          markerPCMEtaMass[i]->DrawMarker(columnsLegendMass2Eta[1]+ offsetMarkerXMass2Eta ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta);
          markerPCMEtaMassMC[i]           = CreateMarkerFromGraph(graphEtaTrueMass[i],columnsLegendMass2Eta[2]+ offsetMarkerXMass2Eta ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta ,scaleMarkerMass2Eta);
          markerPCMEtaMassMC[i]->DrawMarker(columnsLegendMass2Eta[2]+ offsetMarkerXMass2Eta-0.04 ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta);
          counterEta+=1;
      }
  }

  canvasMassWidthEta->Update();
  canvasMassWidthEta->Print(Form("%s/Eta_MassAndWidth.%s",outputDir.Data(),suffix.Data()));


  // **********************************************************************************************************************
  // ******************************** Acceptance * Efficiency for pi0 single measurement 8TeV **************************
  // **********************************************************************************************************************
  textSizeLabelsPixel             = 55;
  Double_t textSizeLabelsRel      = 55./1200;
  cout << textSizeLabelsRel << endl;

  TCanvas* canvasAcceptanceTimesEff       = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
  DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.1, 0.01, 0.015, 0.095);
  canvasAcceptanceTimesEff->SetLogy(1);
  canvasAcceptanceTimesEff->SetLogx(1);

  TH2F * histo2DAccEff;
  histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, minPtPi0, maxPtPi0, 1000, 1.01e-4, 14 );
  SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),
                          0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);//(#times #epsilon_{pur})
  histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
  histo2DAccEff->GetXaxis()->SetNoExponent();
  histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
  histo2DAccEff->DrawCopy();

  for (Int_t i = 0; i < 13; i++){
      if(graphPi0AccTimesEff[i]){
          DrawGammaSetMarkerTGraphAsym(graphPi0AccTimesEff[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
          graphPi0AccTimesEff[i]->Draw("p,same,e");
      }
  }

  TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  for (Int_t i = 0; i < 13; i++){
      if(graphPi0AccTimesEff[i]){
          legendEffiAccPi0->AddEntry(graphPi0AccTimesEff[i],nameMeasGlobalPlot[i].Data(),"p");
      }
  }
  legendEffiAccPi0->Draw();

  drawLatexAdd("ALICE simulation",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystempPb8TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasAcceptanceTimesEff->Update();
  canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));
  
  histo2DAccEff->DrawCopy();

  for (Int_t i = 0; i < 13; i++){
      if(graphEtaAccTimesEff[i]){
          DrawGammaSetMarkerTGraphAsym(graphEtaAccTimesEff[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
          graphEtaAccTimesEff[i]->Draw("p,same,e");
      }
  }

  legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
  for (Int_t i = 0; i < 13; i++){
      if(graphEtaAccTimesEff[i]){
          legendEffiAccPi0->AddEntry(graphEtaAccTimesEff[i],nameMeasGlobalPlot[i].Data(),"p");
      }
  }
  legendEffiAccPi0->Draw();

  drawLatexAdd("ALICE simulation",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(collisionSystempPb8TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#eta #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

  canvasAcceptanceTimesEff->Update();
  canvasAcceptanceTimesEff->Print(Form("%s/Eta_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));

  // **********************************************************************************************************************
  // ******************************** Cross section for pi0 single measurement pPb8TeV ************************************
  // **********************************************************************************************************************

  TCanvas* canvasXSectionPi0  = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
  DrawGammaCanvasSettings( canvasXSectionPi0, 0.16, 0.01, 0.01, 0.08);
  canvasXSectionPi0->SetLogx();
  canvasXSectionPi0->SetLogy();

  TH2F * histo2DXSectionPi0;
  histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,minPtPi0,maxPtPi0,1000,minInvYieldPi0,maxInvYieldPi0);
  SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.75);
  histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
  histo2DXSectionPi0->GetXaxis()->SetNoExponent(kTRUE);
  histo2DXSectionPi0->Draw("copy");

  TLegend* legendXSectionPi0          = new TLegend(0.2,0.15,0.48,0.34);
  legendXSectionPi0->SetFillColor(0);
  legendXSectionPi0->SetLineColor(0);
  legendXSectionPi0->SetTextFont(42);
  legendXSectionPi0->SetTextSize(0.035);

  for (Int_t i = 0; i < 13; i++){
    if(directoryPi0[i]){
      // DrawGammaSetMarkerTGraphAsym(graphPi0InvYieldStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      // graphPi0InvYieldStat[i]->Draw("pEsame");
      DrawGammaSetMarkerTGraphAsym(graphPi0InvYieldStat_WOXErr[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphPi0InvYieldStat_WOXErr[i]->Draw("p,same,z");
      DrawGammaSetMarkerTGraphAsym(graphPi0InvYieldSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
      graphPi0InvYieldSys[i]->Draw("E2same");
      legendXSectionPi0->AddEntry(graphPi0InvYieldSys[i],nameMeasGlobalPlot[i].Data(),"fp");
    }
  }

  TLatex *labelEnergyXSectionPi0      = new TLatex(0.94,0.92,collisionSystempPb8TeV.Data());
  SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4,1,42,kTRUE,31);
  labelEnergyXSectionPi0->Draw();
  TLatex *labelDetSysXSectionPi0      = new TLatex(0.94,0.88,"#pi^{0} #rightarrow #gamma#gamma");
  SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4,1,42,kTRUE,31);
  labelDetSysXSectionPi0->Draw();

  legendXSectionPi0->Draw();

  canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems.%s",outputDir.Data(),suffix.Data()));

  canvasXSectionPi0->cd();
  histo2DXSectionPi0->Draw("copy");

  for (Int_t i = 0; i < 13; i++){
    if(directoryPi0[i]){
      graphPi0InvYieldStat_WOXErr[i]->Draw("p,same,z");
      graphPi0InvYieldSys[i]->Draw("E2same");
    }
  }
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphCombPi0InvYieldSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");

  labelEnergyXSectionPi0->Draw();
  labelDetSysXSectionPi0->Draw();

  legendXSectionPi0->AddEntry(graphCombPi0InvYieldSys,"comb","fp");
  legendXSectionPi0->Draw();

  canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));



  TCanvas* canvasXSectionEta      = new TCanvas("canvasXSectionEta","",200,10,1350,1350*1.15);  // gives the page size
  DrawGammaCanvasSettings( canvasXSectionEta, 0.14, 0.02, 0.02, 0.09);
  canvasXSectionEta->SetLogx();
  canvasXSectionEta->SetLogy();

  TH2F * histo2DXSectionEta;
  histo2DXSectionEta              = new TH2F("histo2DXSectionEta","histo2DXSectionEta",11000,minPtEta, maxPtEta,1000,minInvYieldEta,maxInvYieldEta);
  SetStyleHistoTH2ForGraphs(histo2DXSectionEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 1.,1.45);
  histo2DXSectionEta->GetXaxis()->SetMoreLogLabels();
  histo2DXSectionEta->GetXaxis()->SetNoExponent(kTRUE);
  histo2DXSectionEta->Draw("copy");

  TLegend* legendXSectionEta         = new TLegend(0.2,0.15,0.48,0.302);
  legendXSectionEta->SetFillColor(0);
  legendXSectionEta->SetLineColor(0);
  legendXSectionEta->SetTextFont(42);
  legendXSectionEta->SetTextSize(0.035);
  for (Int_t i = 0; i < 13; i++){
    if(directoryEta[i]){
      DrawGammaSetMarkerTGraphAsym(graphEtaInvYieldStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphEtaInvYieldStat[i]->Draw("pEsame");
      DrawGammaSetMarkerTGraphAsym(graphEtaInvYieldStat_WOXErr[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphEtaInvYieldStat_WOXErr[i]->Draw("p,same,z");
      DrawGammaSetMarkerTGraphAsym(graphEtaInvYieldSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
      graphEtaInvYieldSys[i]->Draw("E2same");
      legendXSectionEta->AddEntry(graphEtaInvYieldSys[i],nameMeasGlobalPlot[i].Data(),"fp");
   }
  }

  TLatex *labelEnergyXSectionEta      = new TLatex(0.64,0.92,collisionSystempPb8TeV.Data());
  SetStyleTLatex( labelEnergyXSectionEta, 0.035,4);
  labelEnergyXSectionEta->Draw();
  TLatex *labelDetSysXSectionEta      = new TLatex(0.64,0.88,"#eta #rightarrow #gamma#gamma");
  SetStyleTLatex( labelDetSysXSectionEta, 0.035,4);
  labelDetSysXSectionEta->Draw();



  legendXSectionEta->Draw();

  canvasXSectionEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems.%s",outputDir.Data(),suffix.Data()));
  histo2DXSectionEta->Draw("copy");

  for (Int_t i = 0; i < 13; i++){
    if(directoryEta[i]){
      graphEtaInvYieldStat_WOXErr[i]->Draw("p,same,z");
      graphEtaInvYieldSys[i]->Draw("E2same");
    }
  }

  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
  graphCombEtaInvYieldSys->Draw("E2same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
  graphCombEtaInvYieldStat_WOXErr->Draw("p,same,z");


  labelEnergyXSectionEta->Draw();
  labelDetSysXSectionEta->Draw();

  legendXSectionEta->AddEntry(graphCombEtaInvYieldSys,"comb","fp");
  legendXSectionEta->Draw();

  canvasXSectionEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));



  canvasXSectionPi0->cd();
  histo2DXSectionPi0->Draw("copy");
  for (Int_t i = 0; i < 13; i++){
    if(directoryPi0[i]){
      DrawGammaSetMarkerTGraphAsym(graphPi0InvYieldStat_WOXErr[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphPi0InvYieldStat_WOXErr[i]->Draw("p,same,z");
      graphPi0InvYieldSys[i]->Draw("E2same");
    }
    if(directoryEta[i]){
      DrawGammaSetMarkerTGraphAsym(graphEtaInvYieldStat_WOXErr[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
      graphEtaInvYieldStat_WOXErr[i]->Draw("p,same,z");
      graphEtaInvYieldSys[i]->Draw("E2same");
    }
  }

  labelEnergyXSectionPi0->Draw();
  labelDetSysXSectionPi0->Draw();

  legendXSectionPi0->Draw();

  canvasXSectionPi0->SaveAs(Form("%s/Pi0_Eta_InvYieldCompAllSystems.%s",outputDir.Data(),suffix.Data()));


  // **********************************************************************************************************************
  // ******************************************* Comparison to theory calculations Pi0 ************************************
  // **********************************************************************************************************************
  textSizeLabelsPixel                     = 48;

  TCanvas* canvasRatioPP                  = new TCanvas("canvasRatioPP","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRatioPP,  0.12, 0.01, 0.01, 0.11);
  canvasRatioPP->cd();
  canvasRatioPP->SetLogx();

      textsizeLabelsPP                    = 0;
      textsizeFacPP                       = 0;
      if (canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) <canvasRatioPP->YtoPixel(canvasRatioPP->GetY1()) ){
          textsizeLabelsPP                = (Double_t)textSizeLabelsPixel/canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
          textsizeFacPP                   = (Double_t)1./canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
      } else {
      textsizeLabelsPP                    = (Double_t)textSizeLabelsPixel/canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
      textsizeFacPP                       = (Double_t)1./canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
      }
      cout << textsizeLabelsPP << endl;

      TH2F * ratio2DTheoryPP       = new TH2F("ratio2DTheoryPP","ratio2DTheoryPP",1000,minPtPi0,maxPtPi0,1000,0.6,2.2);
      SetStyleHistoTH2ForGraphs(ratio2DTheoryPP, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                                0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
      ratio2DTheoryPP->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2DTheoryPP->GetYaxis()->SetNdivisions(505);
      ratio2DTheoryPP->GetYaxis()->SetNoExponent(kTRUE);
      ratio2DTheoryPP->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2DTheoryPP->GetXaxis()->SetNoExponent(kTRUE);
      ratio2DTheoryPP->GetXaxis()->SetLabelFont(42);
      ratio2DTheoryPP->GetYaxis()->SetLabelFont(42);
      ratio2DTheoryPP->DrawCopy();

      graphRatioPi0DSS14->SetLineWidth(widthCommonFit);
      graphRatioPi0DSS14->SetLineColor(colorNLO);
      graphRatioPi0DSS14->SetLineStyle(1);
      graphRatioPi0DSS14->SetFillStyle(1001);
      graphRatioPi0DSS14->SetFillColor(colorNLO);
      graphRatioPi0DSS14->Draw("same,e4");

      // DrawGammaSetMarkerTGraphErr(graphRatioEPOSJJToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphRatioEPOSJJToFit->Draw("3,same");
      DrawGammaSetMarker(histoRatioPi0HIJINGToFit, 24, 1.5, kGreen+2 , kGreen+2);
      histoRatioPi0HIJINGToFit->SetLineWidth(widthCommonFit);
      histoRatioPi0HIJINGToFit->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphRatioEPOSJJToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      graphRatioEPOSJJToFit->Draw("3,same");
      DrawGammaSetMarker(histoRatioEPOSJJToFit, 24, 1.5, kRed+2 , kRed+2);
      histoRatioEPOSJJToFit->SetLineWidth(widthCommonFit);
      histoRatioEPOSJJToFit->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphRatioDPMJETToFit, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      graphRatioDPMJETToFit->Draw("3,same");
      DrawGammaSetMarker(histoRatioDPMJETToFit, 24, 1.5, kOrange+2 , kOrange+2);
      histoRatioDPMJETToFit->SetLineWidth(widthCommonFit);
      histoRatioDPMJETToFit->Draw("same,hist,l");

      DrawGammaSetMarker(histoRatioPythia8EPPS16ToFit, 24, 1.5, kPink+2 , kPink+2);
      histoRatioPythia8EPPS16ToFit->SetLineWidth(widthCommonFit);
      histoRatioPythia8EPPS16ToFit->Draw("same,hist,l");

      DrawGammaSetMarker(histoRatioPythia8nCTEQ15ToFit, 24, 1.5, kCyan+2 , kCyan+2);
      histoRatioPythia8nCTEQ15ToFit->SetLineWidth(widthCommonFit);
      histoRatioPythia8nCTEQ15ToFit->Draw("same,hist,l");

      DrawGammaSetMarker(histoRatioPythia8Monash2013ToFit, 24, 1.5, kGreen+2 , kGreen+2);
      histoRatioPythia8Monash2013ToFit->SetLineWidth(widthCommonFit);
      // histoRatioPythia8Monash2013ToFit->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphRatioCombCombFitSys->SetLineWidth(1);
      graphRatioCombCombFitSys->Draw("2,same");
      graphRatioCombCombFitStat_WOXErr->Draw("p,same");

      // TBox* boxErrorSigmaRatio = CreateBoxConv(kGray+2, 0.25, 1.-(0.035 ), 0.3, 1.+(0.035));
      TBox* boxErrorSigmaRatio = CreateBoxConv(kGray+2, 0.3, 1.-(0.024 ), 0.35, 1.+(0.024));
      boxErrorSigmaRatio->SetLineWidth(8);
      boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0,1., 1.,0.5,kGray+2);

      TLegend* legendRatioTheorypp_3Parted = GetAndSetLegend2(0.15,0.76,0.4,0.96, 0.85* textSizeLabelsPixel);
      legendRatioTheorypp_3Parted->AddEntry(graphRatioCombCombFitSys,"Data","pf");
      // legendRatioTheorypp_3Parted->AddEntry(histoRatioEPOSJJToFit,  "EPOS", "l");
      legendRatioTheorypp_3Parted->AddEntry(histoRatioPi0HIJINGToFit,  "HIJING", "l");
      legendRatioTheorypp_3Parted->AddEntry(histoRatioDPMJETToFit,  "DPMJET", "l");
      legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0DSS14,  "NLO, PDF:CT10 - FF:DSS14", "f");
      legendRatioTheorypp_3Parted->Draw();

      TLegend* legendRatioTheoryNormUnc = GetAndSetLegend2(0.29,0.9,0.54,0.95, 0.85* textSizeLabelsPixel);
      legendRatioTheoryNormUnc->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.4%","l");
      legendRatioTheoryNormUnc->Draw();

      TLatex *labelRatioTheoryPPP   = new TLatex(0.218,0.73,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
      SetStyleTLatex( labelRatioTheoryPPP, 0.85*textsizeLabelsPP,4);
      labelRatioTheoryPPP->Draw();

      // TLatex *labelRatioTheoryPP   = new TLatex(0.76,0.925,collisionSystempPb8TeV.Data());
      TLatex *labelRatioTheoryPP   = new TLatex(0.96,0.925,collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelRatioTheoryPP, 0.85*textsizeLabelsPP,4);
      labelRatioTheoryPP->SetTextAlign(31);
      labelRatioTheoryPP->Draw();
      // TLatex *labelRatioTheoryPP1P = new TLatex(0.863,0.875,labelALICEPlot.Data());
      TLatex *labelRatioTheoryPP1P = new TLatex(0.96,0.875,labelALICEPlot.Data());
      SetStyleTLatex( labelRatioTheoryPP1P, 0.85*textsizeLabelsPP,4);
      labelRatioTheoryPP1P->SetTextAlign(31);
      labelRatioTheoryPP1P->Draw();
      // TLatex *labelRatioTheoryPP2P= new TLatex(0.843,0.83,"#pi^{0} #rightarrow #gamma#gamma");
      TLatex *labelRatioTheoryPP2P= new TLatex(0.96,0.83,"#pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRatioTheoryPP2P, 0.85*textsizeLabelsPP,4);
      labelRatioTheoryPP2P->SetTextAlign(31);
      labelRatioTheoryPP2P->Draw();


  canvasRatioPP->Update();
  canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PPB.%s",outputDir.Data(),suffix.Data()));

      ratio2DTheoryPP->DrawCopy();

      // graphRatioEPOSJJToFit->Draw("3,same");
      // histoRatioEPOSJJToFit->Draw("same,hist,l");

      // graphRatioDPMJETToFit->Draw("3,same");
      // histoRatioDPMJETToFit->Draw("same,hist,l");

      histoRatioPi0HIJINGToFit->Draw("same,hist,l");
      histoRatioPythia8EPPS16ToFit->Draw("same,hist,l");
      histoRatioPythia8nCTEQ15ToFit->Draw("same,hist,l");
      // histoRatioPythia8Monash2013ToFit->Draw("same,hist,l");

      graphRatioCombCombFitSys->Draw("2,same");
      graphRatioCombCombFitStat_WOXErr->Draw("p,same");

      boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0,1., 1.,0.5,kGray+2);

      // legendRatioTheorypp_3Parted = GetAndSetLegend2(0.15,0.76,0.4,0.96, 0.85* textSizeLabelsPixel);
      legendRatioTheorypp_3Parted = GetAndSetLegend2(0.15,0.76,0.4,0.96, 0.85* textSizeLabelsPixel);
      legendRatioTheorypp_3Parted->AddEntry(graphRatioCombCombFitSys,"Data","pf");
      legendRatioTheorypp_3Parted->AddEntry(histoRatioPi0HIJINGToFit,  "HIJING", "l");
      // legendRatioTheorypp_3Parted->AddEntry(histoRatioDPMJETToFit,  "DPMJET", "l");
      legendRatioTheorypp_3Parted->AddEntry(histoRatioPythia8EPPS16ToFit,  "PYTHIA8 EPPS16", "l");
      legendRatioTheorypp_3Parted->AddEntry(histoRatioPythia8nCTEQ15ToFit,  "PYTHIA8 nCTEQ15", "l");
      // legendRatioTheorypp_3Parted->AddEntry(histoRatioPythia8Monash2013ToFit,  "PYTHIA8 Monash2013", "l");
      legendRatioTheorypp_3Parted->Draw();

      legendRatioTheoryNormUnc->Draw();

      labelRatioTheoryPP->Draw();
      labelRatioTheoryPP1P->Draw();
      labelRatioTheoryPP2P->Draw();


  canvasRatioPP->Update();
  canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_noNLO.%s",outputDir.Data(),suffix.Data()));


  TH2F * ratio2DTheoryPP2       = new TH2F("ratio2DTheoryPP2","ratio2DTheoryPP2",1000,minPtPi0,maxPtPi0,1000,0.5,3.6);
  SetStyleHistoTH2ForGraphs(ratio2DTheoryPP2, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
  ratio2DTheoryPP2->GetYaxis()->SetMoreLogLabels(kTRUE);
  ratio2DTheoryPP2->GetYaxis()->SetNdivisions(505);
  ratio2DTheoryPP2->GetYaxis()->SetNoExponent(kTRUE);
  ratio2DTheoryPP2->GetXaxis()->SetMoreLogLabels(kTRUE);
  ratio2DTheoryPP2->GetXaxis()->SetNoExponent(kTRUE);
  ratio2DTheoryPP2->GetXaxis()->SetLabelFont(42);
  ratio2DTheoryPP2->GetYaxis()->SetLabelFont(42);
  ratio2DTheoryPP2->DrawCopy();

  graphRatioPi0DSS14->SetLineWidth(widthCommonFit);
  graphRatioPi0DSS14->SetLineColor(colorNLO);
  graphRatioPi0DSS14->SetLineStyle(1);
  graphRatioPi0DSS14->SetFillStyle(1001);
  graphRatioPi0DSS14->SetFillColor(colorNLO);
  graphRatioPi0DSS14->Draw("same,e4");

  DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
  graphRatioPi0CombNLOMuHalf->Draw("same,c");
  DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
  graphRatioPi0CombNLOMuOne->Draw("same,c");
  DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
  graphRatioPi0CombNLOMuTwo->Draw("same,c");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
  graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
  graphRatioCombCombFitSys->SetLineWidth(1);
  graphRatioCombCombFitSys->Draw("2,same");
  graphRatioCombCombFitStat_WOXErr->Draw("p,same");

  boxErrorSigmaRatio->Draw();
  DrawGammaLines(minPtPi0, maxPtPi0,1., 1.,0.5,kGray+2);

  TLegend* legendRatioTheorypp_3Parted2= GetAndSetLegend2(0.15,0.65,0.4,0.96, 0.85* textSizeLabelsPixel);
  legendRatioTheorypp_3Parted2->AddEntry(graphRatioCombCombFitSys,"Data","pf");
  legendRatioTheorypp_3Parted2->AddEntry((TObject*)0,"NLO, PDF:CT10 - FF:DSS14", "");
  legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuHalf, "#mu = 0.5 #it{p}_{T}", "l");
  legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuOne,  "#mu = #it{p}_{T}", "l");
  legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuTwo,  "#mu = 2 #it{p}_{T}", "l");
  legendRatioTheorypp_3Parted2->Draw();

  TLegend* legendRatioTheoryNormUnc2 = GetAndSetLegend2(0.34,0.902,0.59,0.952, 0.85* textSizeLabelsPixel);
  legendRatioTheoryNormUnc2->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.4%","l");
  legendRatioTheoryNormUnc2->Draw();

  TLegend* legendRatioTheorypp_3Parted22= GetAndSetLegend2(0.15,0.14,0.4,0.18, 0.85* textSizeLabelsPixel);
  legendRatioTheorypp_3Parted22->AddEntry(graphRatioPi0DSS14,  "NLO, PDF:CT10 - FF:DSS14, 0.5#it{p}_{T} < #mu < 2#it{p}_{T}", "f");
  legendRatioTheorypp_3Parted22->Draw();

  labelRatioTheoryPP->Draw();
  labelRatioTheoryPP1P->Draw();
  labelRatioTheoryPP2P->Draw();

  canvasRatioPP->Update();
  canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP2.%s",outputDir.Data(),suffix.Data()));

  TH2F * ratio2DTheoryPP3       = new TH2F("ratio2DTheoryPP3","ratio2DTheoryPP3",1000,minPtPi0,maxPtPi0,1000,0.71,2.49);
  SetStyleHistoTH2ForGraphs(ratio2DTheoryPP3, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
  ratio2DTheoryPP3->GetYaxis()->SetMoreLogLabels(kTRUE);
  ratio2DTheoryPP3->GetYaxis()->SetNdivisions(505);
  ratio2DTheoryPP3->GetYaxis()->SetNoExponent(kTRUE);
  ratio2DTheoryPP3->GetXaxis()->SetMoreLogLabels(kTRUE);
  ratio2DTheoryPP3->GetXaxis()->SetNoExponent(kTRUE);
  ratio2DTheoryPP3->GetXaxis()->SetLabelFont(42);
  ratio2DTheoryPP3->GetYaxis()->SetLabelFont(42);
  ratio2DTheoryPP3->DrawCopy();

  DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
  graphRatioPi0CombNLOMuHalf->Draw("same,c");
  DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
  graphRatioPi0CombNLOMuOne->Draw("same,c");
  DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
  graphRatioPi0CombNLOMuTwo->Draw("same,c");

  graphRatioEPOSJJToFit->Draw("3,same");
  histoRatioEPOSJJToFit->Draw("same,hist,l");

  graphRatioDPMJETToFit->Draw("3,same");
  histoRatioDPMJETToFit->Draw("same,hist,l");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
  graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
  graphRatioCombCombFitSys->SetLineWidth(1);
  graphRatioCombCombFitSys->Draw("2,same");
  graphRatioCombCombFitStat_WOXErr->Draw("p,same");

  boxErrorSigmaRatio->Draw();
  DrawGammaLines(minPtPi0, maxPtPi0,1., 1.,0.5,kGray+2);

  TLegend* legendRatioTheorypp_3Parted3= GetAndSetLegend2(0.15,0.62,0.4,0.96, 0.85* textSizeLabelsPixel);
  legendRatioTheorypp_3Parted3->AddEntry(graphRatioCombCombFitSys,"Data","pf");
  legendRatioTheorypp_3Parted3->AddEntry(histoRatioEPOSJJToFit,  "EPOS", "l");
  legendRatioTheorypp_3Parted3->AddEntry(histoRatioDPMJETToFit,  "DPMJET", "l");
  legendRatioTheorypp_3Parted3->AddEntry((TObject*)0,"NLO, PDF:CT10 - FF:DSS14", "");
  legendRatioTheorypp_3Parted3->AddEntry(graphRatioPi0CombNLOMuHalf, "#mu = 0.5 #it{p}_{T}", "l");
  legendRatioTheorypp_3Parted3->AddEntry(graphRatioPi0CombNLOMuOne,  "#mu = #it{p}_{T}", "l");
  legendRatioTheorypp_3Parted3->AddEntry(graphRatioPi0CombNLOMuTwo,  "#mu = 2 #it{p}_{T}", "l");
  legendRatioTheorypp_3Parted3->Draw();

  TLegend* legendRatioTheoryNormUnc3 = GetAndSetLegend2(0.34,0.902,0.59,0.952, 0.85* textSizeLabelsPixel);
  legendRatioTheoryNormUnc3->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.4%","l");
  legendRatioTheoryNormUnc3->Draw();

  labelRatioTheoryPP->Draw();
  labelRatioTheoryPP1P->Draw();
  labelRatioTheoryPP2P->Draw();

  canvasRatioPP->Update();
  canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP3.%s",outputDir.Data(),suffix.Data()));


  // **********************************************************************************************************************
  // ******************************************* Comparison to theory calculations Eta ************************************
  // **********************************************************************************************************************
  canvasRatioPP->cd();
  canvasRatioPP->SetLogx();

      TH2F * ratio2DTheoryPPEta       = new TH2F("ratio2DTheoryPPEta","ratio2DTheoryPPEta",1000,minPtEta,maxPtEta,1000,0.3,3.95);
      SetStyleHistoTH2ForGraphs(ratio2DTheoryPPEta, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                                0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
      ratio2DTheoryPPEta->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2DTheoryPPEta->GetYaxis()->SetNdivisions(505);
      ratio2DTheoryPPEta->GetYaxis()->SetNoExponent(kTRUE);
      ratio2DTheoryPPEta->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2DTheoryPPEta->GetXaxis()->SetNoExponent(kTRUE);
      ratio2DTheoryPPEta->GetXaxis()->SetLabelFont(42);
      ratio2DTheoryPPEta->GetYaxis()->SetLabelFont(42);
      ratio2DTheoryPPEta->DrawCopy();

      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      graphRatioEtaCombNLOMuHalf->Draw("same,c");
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
      graphRatioEtaCombNLOMuOne->Draw("same,c");
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      graphRatioEtaCombNLOMuTwo->Draw("same,c");

      // DrawGammaSetMarkerTGraphErr(graphRatioEPOSJJToFitEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphRatioEPOSJJToFitEta->Draw("3,same");
      DrawGammaSetMarker(histoRatioEtaHIJINGToFit, 24, 1.5, kGreen+2 , kGreen+2);
      histoRatioEtaHIJINGToFit->SetLineWidth(widthCommonFit);
      histoRatioEtaHIJINGToFit->Draw("same,hist,l");

      // DrawGammaSetMarkerTGraphErr(graphRatioEPOSJJToFitEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphRatioEPOSJJToFitEta->Draw("3,same");
      // DrawGammaSetMarker(histoRatioEPOSJJToFitEta, 24, 1.5, kRed+2 , kRed+2);
      // histoRatioEPOSJJToFitEta->SetLineWidth(widthCommonFit);
      // histoRatioEPOSJJToFitEta->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphRatioDPMJETToFitEta, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      graphRatioDPMJETToFitEta->Draw("3,same");
      DrawGammaSetMarker(histoRatioDPMJETToFitEta, 24, 1.5, kOrange+2 , kOrange+2);
      histoRatioDPMJETToFitEta->SetLineWidth(widthCommonFit);
      histoRatioDPMJETToFitEta->Draw("same,hist,l");

      DrawGammaSetMarker(histoRatioPythia8EPPS16ToFitEta, 24, 1.5, kPink+2 , kPink+2);
      histoRatioPythia8EPPS16ToFitEta->SetLineWidth(widthCommonFit);
      histoRatioPythia8EPPS16ToFitEta->Draw("same,hist,l");

      DrawGammaSetMarker(histoRatioPythia8nCTEQ15ToFitEta, 24, 1.5, kCyan+2 , kCyan+2);
      histoRatioPythia8nCTEQ15ToFitEta->SetLineWidth(widthCommonFit);
      histoRatioPythia8nCTEQ15ToFitEta->Draw("same,hist,l");

      DrawGammaSetMarker(histoRatioPythia8Monash2013ToFitEta, 24, 1.5, kGreen+2 , kGreen+2);
      histoRatioPythia8Monash2013ToFitEta->SetLineWidth(widthCommonFit);
      // histoRatioPythia8Monash2013ToFitEta->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatEta_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysEta, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphRatioCombCombFitSysEta->SetLineWidth(1);
      graphRatioCombCombFitSysEta->Draw("2,same");
      graphRatioCombCombFitStatEta_WOXErr->Draw("p,same");

      TBox* boxErrorSigmaRatioEta = CreateBoxConv(kGray+2, 0.4, 1.-(0.035 ), 0.45, 1.+(0.035));
      boxErrorSigmaRatioEta->SetLineWidth(8);
      boxErrorSigmaRatioEta->Draw();
      DrawGammaLines(minPtEta, maxPtEta, 1., 1.,1, kGray+2);

      TLegend* legendRatioTheoryppEta_3Parted= GetAndSetLegend2(0.15,0.84-(0.85*textsizeLabelsPP*5),0.40,0.96, 0.85* textSizeLabelsPixel);
      legendRatioTheoryppEta_3Parted->AddEntry(graphRatioCombCombFitSysEta,"Data","pf");
      legendRatioTheoryppEta_3Parted->AddEntry(histoRatioEtaHIJINGToFit,  "HIJING", "l");
      // legendRatioTheoryppEta_3Parted->AddEntry(histoRatioEPOSJJToFitEta,  "EPOS", "l");
      legendRatioTheoryppEta_3Parted->AddEntry(histoRatioDPMJETToFitEta,  "DPMJET", "l");
      legendRatioTheoryppEta_3Parted->AddEntry((TObject*)0, "NLO, PDF:CTEQ6M5 - FF:AESSS", "");
      legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuHalf, "#mu = 0.5 #it{p}_{T}", "l");
      legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuOne,  "#mu = #it{p}_{T}", "l");
      legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuTwo,  "#mu = 2 #it{p}_{T}", "l");
      legendRatioTheoryppEta_3Parted->Draw();

      TLegend* legendRatioTheoryNormUncEta = GetAndSetLegend2(0.31,0.908,0.56,0.958, 0.85* textSizeLabelsPixel);
      legendRatioTheoryNormUncEta->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.4%","l");
      legendRatioTheoryNormUncEta->Draw();

      TLatex *labelRatioTheoryPP2   = new TLatex(0.78,0.925,collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelRatioTheoryPP2, 0.85*textsizeLabelsPP,4);
      labelRatioTheoryPP2->Draw();
      // TLatex *labelRatioTheoryPP221 = new TLatex(0.883,0.875,labelALICEPlot.Data());
      TLatex *labelRatioTheoryPP221 = new TLatex(0.965,0.875,labelALICEPlot.Data());
      SetStyleTLatex( labelRatioTheoryPP221, 0.85*textsizeLabelsPP,4);
      labelRatioTheoryPP221->SetTextAlign(31);
      labelRatioTheoryPP221->Draw();
      TLatex *labelRatioTheoryPP222= new TLatex(0.873,0.83,"#eta #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRatioTheoryPP222, 0.85*textsizeLabelsPP,4);
      labelRatioTheoryPP222->Draw();

  canvasRatioPP->Update();
  canvasRatioPP->Print(Form("%s/Eta_RatioTheoryToData_PPB.%s",outputDir.Data(),suffix.Data()));
      ratio2DTheoryPPEta->DrawCopy();

      // graphRatioEPOSJJToFitEta->Draw("3,same");
      // histoRatioEPOSJJToFitEta->Draw("same,hist,l");
      histoRatioEtaHIJINGToFit->Draw("same,hist,l");

      graphRatioDPMJETToFitEta->Draw("3,same");
      histoRatioDPMJETToFitEta->Draw("same,hist,l");

      histoRatioPythia8EPPS16ToFitEta->Draw("same,hist,l");

      histoRatioPythia8nCTEQ15ToFitEta->Draw("same,hist,l");
      // histoRatioPythia8Monash2013ToFitEta->Draw("same,hist,l");

      graphRatioCombCombFitSysEta->Draw("2,same");
      graphRatioCombCombFitStatEta_WOXErr->Draw("p,same");

      boxErrorSigmaRatioEta->Draw();
      DrawGammaLines(minPtEta, maxPtEta, 1., 1.,1, kGray+2);

      legendRatioTheoryppEta_3Parted= GetAndSetLegend2(0.15,0.84-(0.85*textsizeLabelsPP*2),0.40,0.96, 0.85* textSizeLabelsPixel);
      legendRatioTheoryppEta_3Parted->AddEntry(graphRatioCombCombFitSysEta,"Data","pf");
      // legendRatioTheoryppEta_3Parted->AddEntry(histoRatioEPOSJJToFitEta,  "EPOS", "l");
      legendRatioTheoryppEta_3Parted->AddEntry(histoRatioEtaHIJINGToFit,  "HIJING", "l");
      legendRatioTheoryppEta_3Parted->AddEntry(histoRatioDPMJETToFitEta,  "DPMJET", "l");
      legendRatioTheoryppEta_3Parted->AddEntry(histoRatioPythia8EPPS16ToFitEta,  "PYTHIA8 EPPS16", "l");
      legendRatioTheoryppEta_3Parted->AddEntry(histoRatioPythia8nCTEQ15ToFitEta,  "PYTHIA8 nCTEQ15", "l");
      // legendRatioTheoryppEta_3Parted->AddEntry(histoRatioPythia8Monash2013ToFitEta,  "PYTHIA8 Monash2013", "l");
      legendRatioTheoryppEta_3Parted->Draw();

      legendRatioTheoryNormUncEta->Draw();

      labelRatioTheoryPP2->Draw();
      labelRatioTheoryPP221->Draw();
      labelRatioTheoryPP222->Draw();

  canvasRatioPP->Update();
  canvasRatioPP->Print(Form("%s/Eta_RatioTheoryToData_PPB_noNLO.%s",outputDir.Data(),suffix.Data()));


  TGraphAsymmErrors* graphRatioChargedHadronsStatpPb8TeV     = (TGraphAsymmErrors*)graphChargedHadronsStatpPb8TeV->Clone();
  graphRatioChargedHadronsStatpPb8TeV                        = CalculateGraphErrRatioToFit(graphRatioChargedHadronsStatpPb8TeV, fitTCMInvYieldPi0);
  TGraphAsymmErrors* graphRatioChargedHadronsStatpPb8TeV_WOXErr = (TGraphAsymmErrors*) graphRatioChargedHadronsStatpPb8TeV->Clone("graphRatioChargedHadronsStatpPb8TeV_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioChargedHadronsStatpPb8TeV_WOXErr);
  TGraphAsymmErrors* graphRatioChargedHadronsSyspPb8TeV     = (TGraphAsymmErrors*)graphChargedHadronsSyspPb8TeV->Clone();
  graphRatioChargedHadronsSyspPb8TeV                        = CalculateGraphErrRatioToFit(graphRatioChargedHadronsSyspPb8TeV, fitTCMInvYieldPi0);
  
  TGraphAsymmErrors* graphRatioChargedPionStatpPb8TeV     = (TGraphAsymmErrors*)graphChargedPionStatpPb8TeV->Clone();
  graphRatioChargedPionStatpPb8TeV                        = CalculateGraphErrRatioToFit(graphRatioChargedPionStatpPb8TeV, fitTCMInvYieldPi0);
  TGraphAsymmErrors* graphRatioChargedPionStatpPb8TeV_WOXErr = (TGraphAsymmErrors*) graphRatioChargedPionStatpPb8TeV->Clone("graphRatioChargedPionStatpPb8TeV_WOXErr");
  ProduceGraphAsymmWithoutXErrors(graphRatioChargedPionStatpPb8TeV_WOXErr);
  TGraphAsymmErrors* graphRatioChargedPionSyspPb8TeV     = (TGraphAsymmErrors*)graphChargedPionSyspPb8TeV->Clone();
  graphRatioChargedPionSyspPb8TeV                        = CalculateGraphErrRatioToFit(graphRatioChargedPionSyspPb8TeV, fitTCMInvYieldPi0);

      TH2F * ratio2DChargedPP       = new TH2F("ratio2DChargedPP","ratio2DChargedPP",1000,minPtPi0,maxPtPi0,1000,0.0,5.6);
      SetStyleHistoTH2ForGraphs(ratio2DChargedPP, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{#pi^{0} TCM fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                                0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
      ratio2DChargedPP->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2DChargedPP->GetYaxis()->SetRangeUser(0.61,5.6);
      ratio2DChargedPP->GetYaxis()->SetNdivisions(505);
      ratio2DChargedPP->GetYaxis()->SetNoExponent(kTRUE);
      ratio2DChargedPP->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2DChargedPP->GetXaxis()->SetNoExponent(kTRUE);
      ratio2DChargedPP->GetXaxis()->SetLabelFont(42);
      ratio2DChargedPP->GetYaxis()->SetLabelFont(42);
      ratio2DChargedPP->DrawCopy();

      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, 29, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, 29, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphRatioCombCombFitSys->SetLineWidth(1);
      graphRatioCombCombFitSys->Draw("2,same");
      graphRatioCombCombFitStat_WOXErr->Draw("p,same");
      
      DrawGammaSetMarkerTGraphAsym(graphRatioChargedHadronsStatpPb8TeV_WOXErr, markerStyleComb, markerSizeComb, kAzure+2, kAzure+2, widthLinesBoxes, kFALSE);
      graphRatioChargedHadronsStatpPb8TeV_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioChargedHadronsSyspPb8TeV, markerStyleComb, markerSizeComb, kAzure+2, kAzure+2, widthLinesBoxes, kTRUE, 0);
      graphRatioChargedHadronsSyspPb8TeV->SetLineWidth(1);
      graphRatioChargedHadronsSyspPb8TeV->Draw("2,same");
      graphRatioChargedHadronsStatpPb8TeV_WOXErr->Draw("p,same");
      
      DrawGammaSetMarkerTGraphAsym(graphRatioChargedPionStatpPb8TeV_WOXErr, 24, markerSizeComb, kMagenta+2, kMagenta+2, widthLinesBoxes, kFALSE);
      graphRatioChargedPionStatpPb8TeV_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioChargedPionSyspPb8TeV, 24, markerSizeComb, kMagenta+2, kMagenta+2, widthLinesBoxes, kTRUE, 0);
      graphRatioChargedPionSyspPb8TeV->SetLineWidth(1);
      graphRatioChargedPionSyspPb8TeV->Draw("2,same");
      graphRatioChargedPionStatpPb8TeV_WOXErr->Draw("p,same");

      boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0,1., 1.,0.5,kGray+2);

      TLegend* legendRatioChargedpp_3Parted = GetAndSetLegend2(0.15,0.76,0.4,0.96, 0.85* textSizeLabelsPixel);
      legendRatioChargedpp_3Parted->AddEntry(graphRatioCombCombFitSys,"Data","pf");
      legendRatioChargedpp_3Parted->AddEntry(graphRatioChargedHadronsSyspPb8TeV,"Charged hadrons","pf");
      legendRatioChargedpp_3Parted->AddEntry(graphRatioChargedPionSyspPb8TeV,"#pi^{#pm}/2","pf");
      legendRatioChargedpp_3Parted->Draw();

      legendRatioTheoryNormUnc->Draw();

      labelRatioTheoryPP->Draw();
      labelRatioTheoryPP1P->Draw();
      labelRatioTheoryPP2P->Draw();


  canvasRatioPP->Update();
  canvasRatioPP->Print(Form("%s/Pi0_RatioChargedToData_PPB.%s",outputDir.Data(),suffix.Data()));
  
      ratio2DChargedPP->GetYaxis()->SetRangeUser(0.65,1.55);
      ratio2DChargedPP->DrawCopy();
      graphRatioCombCombFitSys->Draw("2,same");
      graphRatioCombCombFitStat_WOXErr->Draw("p,same");

      graphRatioChargedPionSyspPb8TeV->Draw("2,same");
      graphRatioChargedPionStatpPb8TeV_WOXErr->Draw("p,same");

      boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0,1., 1.,0.5,kGray+2);

      TLegend* legendRatioChargedpp_3Parted_2 = GetAndSetLegend2(0.15,0.83,0.4,0.96, 0.85* textSizeLabelsPixel);
      legendRatioChargedpp_3Parted_2->AddEntry(graphRatioCombCombFitSys,"Data","pf");
      // legendRatioChargedpp_3Parted->AddEntry(graphRatioChargedHadronsSyspPb8TeV,"Charged hadrons","pf");
      legendRatioChargedpp_3Parted_2->AddEntry(graphRatioChargedPionSyspPb8TeV,"#pi^{#pm}/2","pf");
      legendRatioChargedpp_3Parted_2->Draw();

      legendRatioTheoryNormUnc->Draw();

      labelRatioTheoryPP->Draw();
      labelRatioTheoryPP1P->Draw();
      labelRatioTheoryPP2P->Draw();


  canvasRatioPP->Update();
  canvasRatioPP->Print(Form("%s/Pi0_RatioChargedPionToData_PPB.%s",outputDir.Data(),suffix.Data()));
  
  TGraphErrors* graphRatioBinByBinCPionToPi0AStat = NULL;
  TGraphErrors* graphRatioBinByBinCPionToPi0ASys  = NULL;
  TGraphErrors* graphRatioBinByBinCPionToPi0BStat = NULL;
  TGraphErrors* graphRatioBinByBinCPionToPi0BSys  = NULL;
  TGraphAsymmErrors* graphCPionToPi0inputASys = (TGraphAsymmErrors*) graphChargedPionSyspPb8TeV->Clone();
  TGraphAsymmErrors* graphCPionToPi0inputBSys = (TGraphAsymmErrors*) graphCombPi0InvYieldSys->Clone();

  TGraphErrors* graphRatioBinByBinCPionToPi0 = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                      graphChargedPionStatpPb8TeV->Clone(), graphCPionToPi0inputASys,
                                                                                                      graphCombPi0InvYieldStat->Clone(), graphCPionToPi0inputBSys,
                                                                                                      kTRUE,  kTRUE,
                                                                                                      &graphRatioBinByBinCPionToPi0AStat, &graphRatioBinByBinCPionToPi0ASys,
                                                                                                      &graphRatioBinByBinCPionToPi0BStat, &graphRatioBinByBinCPionToPi0BSys )    ;
  TGraphAsymmErrors* graphCPionToPi0inputASys_ForStat = (TGraphAsymmErrors*) graphChargedPionSyspPb8TeV->Clone();
  TGraphAsymmErrors* graphCPionToPi0inputBSys_ForStat = (TGraphAsymmErrors*) graphCombPi0InvYieldSys->Clone();
    for(Int_t i=0; i<graphCPionToPi0inputASys_ForStat->GetN();i++){
      graphCPionToPi0inputASys_ForStat->GetEYlow()[i] = (1./1000.)*graphCPionToPi0inputASys_ForStat->GetY()[i];
      graphCPionToPi0inputASys_ForStat->GetEYhigh()[i] = (1./1000.)*graphCPionToPi0inputASys_ForStat->GetY()[i];
    }

    for(Int_t i=0; i<graphCPionToPi0inputBSys_ForStat->GetN();i++){
      graphCPionToPi0inputBSys_ForStat->GetEYlow()[i] = (1./1000.)*graphCPionToPi0inputBSys_ForStat->GetY()[i];
      graphCPionToPi0inputBSys_ForStat->GetEYhigh()[i] = (1./1000.)*graphCPionToPi0inputBSys_ForStat->GetY()[i];
    }
  TGraphErrors* graphRatioBinByBinCPionToPi0stat = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                      graphChargedPionStatpPb8TeV->Clone(), graphCPionToPi0inputASys_ForStat,
                                                                                                      graphCombPi0InvYieldStat->Clone(), graphCPionToPi0inputBSys_ForStat,
                                                                                                      kTRUE,  kTRUE,
                                                                                                      &graphRatioBinByBinCPionToPi0AStat, &graphRatioBinByBinCPionToPi0ASys,
                                                                                                      &graphRatioBinByBinCPionToPi0BStat, &graphRatioBinByBinCPionToPi0BSys )    ;
    TGraphErrors* graphRatioBinByBinCPionToPi0stat_WOXErr = (TGraphErrors*) graphRatioBinByBinCPionToPi0stat->Clone("graphRatioBinByBinCPionToPi0stat_WOXErr");
    ProduceGraphErrWithoutXErrors(graphRatioBinByBinCPionToPi0stat_WOXErr);
      ratio2DChargedPP->GetYaxis()->SetRangeUser(0.65,1.55);
      ratio2DChargedPP->GetYaxis()->SetTitle("#pi^{#pm}/#pi^{0}");
      ratio2DChargedPP->DrawCopy();

      DrawGammaSetMarkerTGraphErr(graphRatioBinByBinCPionToPi0, 20, 2, kBlack, kBlack, widthLinesBoxes, kTRUE);
      graphRatioBinByBinCPionToPi0->Draw("E2same");
      for(Int_t i=0; i<graphRatioBinByBinCPionToPi0stat_WOXErr->GetN(); i++) graphRatioBinByBinCPionToPi0stat_WOXErr->GetEX()[i] = 0.;
      DrawGammaSetMarkerTGraphErr(graphRatioBinByBinCPionToPi0stat_WOXErr, 20, 2, kBlack, kBlack);
      graphRatioBinByBinCPionToPi0stat_WOXErr->Draw("p,same,z");

      // boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0,1., 1.,0.5,kGray+2);


      labelRatioTheoryPP->Draw();
      labelRatioTheoryPP1P->Draw();
      labelRatioTheoryPP2P->Draw();


  canvasRatioPP->Update();
  canvasRatioPP->Print(Form("%s/Pi0_RatioChargedPionToData_PP_BinByBin.%s",outputDir.Data(),suffix.Data()));
  
  TGraphErrors* graphRatioBinByBinEtaToCKaonAStat = NULL;
  TGraphErrors* graphRatioBinByBinEtaToCKaonASys  = NULL;
  TGraphErrors* graphRatioBinByBinEtaToCKaonBStat = NULL;
  TGraphErrors* graphRatioBinByBinEtaToCKaonBSys  = NULL;
  TGraphAsymmErrors* graphEtaToCKaoninputASys = (TGraphAsymmErrors*) graphCombEtaInvYieldSys->Clone();
  TGraphAsymmErrors* graphEtaToCKaoninputBSys = (TGraphAsymmErrors*) graphChargedKaonSyspPb8TeV->Clone();

  TGraphErrors* graphRatioBinByBinEtaToCKaon = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                      graphCombEtaInvYieldStat->Clone(), graphEtaToCKaoninputASys,
                                                                                                      graphChargedKaonStatpPb8TeV->Clone(), graphEtaToCKaoninputBSys,
                                                                                                      kTRUE,  kTRUE,
                                                                                                      &graphRatioBinByBinEtaToCKaonAStat, &graphRatioBinByBinEtaToCKaonASys,
                                                                                                      &graphRatioBinByBinEtaToCKaonBStat, &graphRatioBinByBinEtaToCKaonBSys )    ;
  TGraphAsymmErrors* graphEtaToCKaoninputASys_ForStat = (TGraphAsymmErrors*) graphCombEtaInvYieldSys->Clone();
  TGraphAsymmErrors* graphEtaToCKaoninputBSys_ForStat = (TGraphAsymmErrors*) graphChargedKaonSyspPb8TeV->Clone();
    for(Int_t i=0; i<graphEtaToCKaoninputASys_ForStat->GetN();i++){
      graphEtaToCKaoninputASys_ForStat->GetEYlow()[i] = (1./1000.)*graphEtaToCKaoninputASys_ForStat->GetY()[i];
      graphEtaToCKaoninputASys_ForStat->GetEYhigh()[i] = (1./1000.)*graphEtaToCKaoninputASys_ForStat->GetY()[i];
    }

    for(Int_t i=0; i<graphEtaToCKaoninputBSys_ForStat->GetN();i++){
      graphEtaToCKaoninputBSys_ForStat->GetEYlow()[i] = (1./1000.)*graphEtaToCKaoninputBSys_ForStat->GetY()[i];
      graphEtaToCKaoninputBSys_ForStat->GetEYhigh()[i] = (1./1000.)*graphEtaToCKaoninputBSys_ForStat->GetY()[i];
    }
  TGraphErrors* graphRatioBinByBinEtaToCKaonstat = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                      graphCombEtaInvYieldStat->Clone(), graphEtaToCKaoninputASys_ForStat,
                                                                                                      graphChargedKaonStatpPb8TeV->Clone(), graphEtaToCKaoninputBSys_ForStat,
                                                                                                      kTRUE,  kTRUE,
                                                                                                      &graphRatioBinByBinEtaToCKaonAStat, &graphRatioBinByBinEtaToCKaonASys,
                                                                                                      &graphRatioBinByBinEtaToCKaonBStat, &graphRatioBinByBinEtaToCKaonBSys )    ;
    TGraphErrors* graphRatioBinByBinEtaToCKaonstat_WOXErr = (TGraphErrors*) graphRatioBinByBinEtaToCKaonstat->Clone("graphRatioBinByBinEtaToCKaonstat_WOXErr");
    ProduceGraphErrWithoutXErrors(graphRatioBinByBinEtaToCKaonstat_WOXErr);
      ratio2DChargedPP->GetYaxis()->SetRangeUser(0.0,1.01);
      ratio2DChargedPP->GetYaxis()->SetTitle("#eta/K^{#pm}");
      ratio2DChargedPP->DrawCopy();

      DrawGammaSetMarkerTGraphErr(graphRatioBinByBinEtaToCKaon, 20, 2, kBlack, kBlack, widthLinesBoxes, kTRUE);
      graphRatioBinByBinEtaToCKaon->Draw("E2same");
      for(Int_t i=0; i<graphRatioBinByBinEtaToCKaonstat_WOXErr->GetN(); i++) graphRatioBinByBinEtaToCKaonstat_WOXErr->GetEX()[i] = 0.;
      DrawGammaSetMarkerTGraphErr(graphRatioBinByBinEtaToCKaonstat_WOXErr, 20, 2, kBlack, kBlack);
      graphRatioBinByBinEtaToCKaonstat_WOXErr->Draw("p,same,z");

      // boxErrorSigmaRatio->Draw();
      // DrawGammaLines(minPtPi0, maxPtPi0,1., 1.,0.5,kGray+2);


      labelRatioTheoryPP->Draw();
      labelRatioTheoryPP1P->Draw();
      labelRatioTheoryPP2P->Draw();


  canvasRatioPP->Update();
  canvasRatioPP->Print(Form("%s/Eta_RatioEtaToChargedKaon_PP_BinByBin.%s",outputDir.Data(),suffix.Data()));

  //*************************************************************************************************************
  //***************************** Paper plot X-section and ratios ***********************************************
  //*************************************************************************************************************

  Double_t arrayBoundariesX1_XSec[2];
  Double_t arrayBoundariesY1_XSec[6];
  Double_t relativeMarginsXXSec[3];
  Double_t relativeMarginsYXSec[3];
  textSizeLabelsPixel = 48;
  ReturnCorrectValuesForCanvasScaling(1250,2000, 1, 5,0.135, 0.005, 0.003,0.05,arrayBoundariesX1_XSec,arrayBoundariesY1_XSec,relativeMarginsXXSec,relativeMarginsYXSec);

  TCanvas* canvasInvSectionPaper      = new TCanvas("canvasInvSectionPaper","",0,0,1250,2000);  // gives the page size
  DrawGammaCanvasSettings( canvasInvSectionPaper,  0.13, 0.02, 0.03, 0.06);

  TPad* padInvSectionSpec             = new TPad("padInvSectionSpec", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[3], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[0],-1, -1, -2);
  DrawGammaPadSettings( padInvSectionSpec, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[0], relativeMarginsYXSec[1]);
  padInvSectionSpec->Draw();
  Double_t marginXSec                 = relativeMarginsXXSec[0]*1250;
  Double_t textsizeLabelsXSecUp       = 0;
  Double_t textsizeFacXSecUp          = 0;
  if (padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) < padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1())){
      textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
      textsizeFacXSecUp               = (Double_t)1./padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
  } else {
      textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
      textsizeFacXSecUp               = (Double_t)1./padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
  }

  TPad* padInvSectionNLORatio         = new TPad("padInvSectionNLORatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[4], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[3],-1, -1, -2);
  DrawGammaPadSettings( padInvSectionNLORatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
  padInvSectionNLORatio->Draw();
  Double_t textsizeLabelsXSecMiddle   = 0;
  Double_t textsizeFacXSecMiddle      = 0;
  if (padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) < padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1())){
      textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
      textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
  } else {
      textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
      textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
  }

  TPad* padInvSectionPythiaRatio      = new TPad("padInvSectionPythiaRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[5], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[4],-1, -1, -2);
  DrawGammaPadSettings( padInvSectionPythiaRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[2]);
  padInvSectionPythiaRatio->Draw();
  Double_t textsizeLabelsXSecDown     = 0;
  Double_t textsizeFacXSecDown        = 0;
  if (padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) < padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1())){
      textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
      textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
  } else {
      textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
      textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
  }


  padInvSectionSpec->cd();
  padInvSectionSpec->SetLogy(1);
  padInvSectionSpec->SetLogx(1);
      SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                              0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
      histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
      histo2DXSectionPi0->GetXaxis()->SetLabelOffset(+0.01);
      histo2DXSectionPi0->GetYaxis()->SetRangeUser(minInvYieldPi0,6E11);
      histo2DXSectionPi0->Draw();

      graphPi0DSS14->SetLineWidth(widthCommonFit);
      graphPi0DSS14->SetLineColor(colorNLO);
      graphPi0DSS14->SetLineStyle(1);
      graphPi0DSS14->SetFillStyle(1001);
      graphPi0DSS14->SetFillColor(colorNLO);
      graphPi0DSS14->Draw("same,e3");

      // DrawGammaSetMarkerTGraphErr(graphEPOSJJInvYield, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphEPOSJJInvYield->Draw("3,same");
      DrawGammaSetMarker(histoHIJINGPi0, 24, 1.5, kGreen+2 , kGreen+2);
      histoHIJINGPi0->SetLineWidth(widthCommonFit);
      histoHIJINGPi0->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphEPOSJJInvYield, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      graphEPOSJJInvYield->Draw("3,same");
      DrawGammaSetMarker(histoEPOSJJInvYield, 24, 1.5, kRed+2 , kRed+2);
      histoEPOSJJInvYield->SetLineWidth(widthCommonFit);
      histoEPOSJJInvYield->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphDPMJETInvYield, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      graphDPMJETInvYield->Draw("3,same");
      DrawGammaSetMarker(histoDPMJETInvYield, 24, 1.5, kOrange+2 , kOrange+2);
      histoDPMJETInvYield->SetLineWidth(widthCommonFit);
      histoDPMJETInvYield->Draw("same,hist,l");

      DrawGammaSetMarkerTF1( fitTCMInvYieldPi0, 7, 2, kGray+2);
      fitTCMInvYieldPi0->Draw("same");

      DrawGammaSetMarkerTF1( fitInvXSectionPi0, 3, 2, kGray+1);
      fitInvXSectionPi0->Draw("same");

      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
      graphCombPi0InvYieldSys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
      graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");

      // TLatex *labelEnergyXSectionPaper= new TLatex(0.72, 0.91, collisionSystempPb8TeV.Data());
      TLatex *labelEnergyXSectionPaper= new TLatex(0.95, 0.91, collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4);
      labelEnergyXSectionPaper->SetTextAlign(31);
      labelEnergyXSectionPaper->Draw();
      // TLatex *labelALICEXSectionPaper= new TLatex(0.848,0.86,labelALICEPlot.Data());
      TLatex *labelALICEXSectionPaper= new TLatex(0.95,0.86,labelALICEPlot.Data());
      SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4);
      labelALICEXSectionPaper->SetTextAlign(31);
      labelALICEXSectionPaper->Draw();
      // TLatex *labelDetSysXSectionPaper= new TLatex(0.824,0.815,"#pi^{0} #rightarrow #gamma#gamma");
      TLatex *labelDetSysXSectionPaper= new TLatex(0.95,0.815,"#pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4);
      labelDetSysXSectionPaper->SetTextAlign(31);
      labelDetSysXSectionPaper->Draw();

      TLegend* legendXsectionPaper    = GetAndSetLegend2(0.17, 0.08, 0.5, 0.18+0.05*5, textSizeLabelsPixel);
      legendXsectionPaper->SetNColumns(1);
      legendXsectionPaper->SetMargin(0.2);
      legendXsectionPaper->AddEntry(graphCombPi0InvYieldSys,"Data","pf");
      // legendXsectionPaper->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.4%","l");
      legendXsectionPaper->AddEntry(fitTCMInvYieldPi0,"TCM fit","l");
      legendXsectionPaper->AddEntry(fitInvXSectionPi0,"Levy-Tsallis fit","l");
      legendXsectionPaper->AddEntry(histoEPOSJJInvYield,"EPOS","l");
      legendXsectionPaper->AddEntry(histoHIJINGPi0,"HIJING","l");
      legendXsectionPaper->AddEntry(histoDPMJETInvYield,"DPMJET","l");
      legendXsectionPaper->AddEntry(graphPi0DSS14,  "NLO, PDF:CT10 - FF:DSS14", "f");
      legendXsectionPaper->Draw();

      TLatex *labelRatioTheoryPP_Paper   = new TLatex(0.24,0.055,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
      SetStyleTLatex( labelRatioTheoryPP_Paper, 0.8*textsizeLabelsPP,4);
      labelRatioTheoryPP_Paper->Draw();

  padInvSectionNLORatio->cd();
  padInvSectionNLORatio->SetLogx(1);
  Double_t minRatio=0.51;
  Double_t maxRatio=1.65;
      TH2F * ratio2DNLO               = new TH2F("ratio2DNLO","ratio2DNLO",1000,minPtPi0,maxPtPi0,1000,minRatio,maxRatio);
      SetStyleHistoTH2ForGraphs(ratio2DNLO, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{TCM fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
      ratio2DNLO->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2DNLO->GetYaxis()->SetNdivisions(505);
      ratio2DNLO->GetYaxis()->SetNoExponent(kTRUE);
      ratio2DNLO->GetYaxis()->CenterTitle(kTRUE);
      ratio2DNLO->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2DNLO->GetXaxis()->SetNoExponent(kTRUE);
      ratio2DNLO->GetXaxis()->SetLabelFont(42);
      ratio2DNLO->GetYaxis()->SetLabelFont(42);
      ratio2DNLO->GetYaxis()->SetLabelOffset(+0.01);
      ratio2DNLO->GetXaxis()->SetTickLength(0.07);
      ratio2DNLO->DrawCopy();

//        TLegend* legendXsectionPaperPi02     = GetAndSetLegend2(0.17, 0.8, 0.4, 0.83+0.05*1, textSizeLabelsPixel*0.8);
//        legendXsectionPaperPi02->AddEntry(boxErrorSigmaRatio,"norm. unc. 3.5%","f");
//        legendXsectionPaperPi02->Draw();

//        graphRatioPi0CombNLODSS14->RemovePoint(0);
//        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
//        graphRatioPi0CombNLODSS14->Draw("3,same");
//        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
//        graphRatioPi0CombNLOMuHalf->Draw("same,c");
//        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
//        graphRatioPi0CombNLOMuOne->Draw("same,c");
//        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
//        graphRatioPi0CombNLOMuTwo->Draw("same,c");
      graphRatioPi0DSS14->SetLineWidth(widthCommonFit);
      graphRatioPi0DSS14->SetLineColor(colorNLO);
      graphRatioPi0DSS14->SetLineStyle(1);
      graphRatioPi0DSS14->SetFillStyle(1001);
      graphRatioPi0DSS14->SetFillColor(colorNLO);
      graphRatioPi0DSS14->Draw("same,e4");

      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphRatioCombCombFitSys->SetLineWidth(1);
      graphRatioCombCombFitSys->Draw("2,same");
      graphRatioCombCombFitStat_WOXErr->Draw("p,same");

      boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0, 1., 1.,1, kGray+2);

  padInvSectionPythiaRatio->cd();
  padInvSectionPythiaRatio->SetLogx(1);
      TH2F * ratio2DPythia            = new TH2F("ratio2DPythia","ratio2DPythia",1000,minPtPi0,maxPtPi0,1000,minRatio,maxRatio);
      SetStyleHistoTH2ForGraphs(ratio2DPythia, "#it{p}_{T} (GeV/#it{c})","#frac{Pythia, Data}{TCM fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 505);
      ratio2DPythia->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2DPythia->GetYaxis()->SetNdivisions(505);
      ratio2DPythia->GetYaxis()->SetNoExponent(kTRUE);
      ratio2DPythia->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2DPythia->GetXaxis()->SetNoExponent(kTRUE);
      ratio2DPythia->GetXaxis()->SetLabelFont(42);
      ratio2DPythia->GetYaxis()->SetLabelFont(42);
      ratio2DPythia->GetYaxis()->CenterTitle(kTRUE);
      ratio2DPythia->GetYaxis()->SetLabelOffset(+0.01);
      ratio2DPythia->GetXaxis()->SetTickLength(0.06);
      ratio2DPythia->GetYaxis()->SetTickLength(0.04);
      ratio2DPythia->DrawCopy();

      // DrawGammaSetMarkerTGraphErr(graphRatioEPOSJJToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphRatioEPOSJJToFit->Draw("3,same");
      DrawGammaSetMarker(histoRatioPi0HIJINGToFit, 24, 1.5, kGreen+2 , kGreen+2);
      histoRatioPi0HIJINGToFit->SetLineWidth(widthCommonFit);
      histoRatioPi0HIJINGToFit->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphRatioEPOSJJToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      graphRatioEPOSJJToFit->Draw("3,same");
      DrawGammaSetMarker(histoRatioEPOSJJToFit, 24, 1.5, kRed+2 , kRed+2);
      histoRatioEPOSJJToFit->SetLineWidth(widthCommonFit);
      histoRatioEPOSJJToFit->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphRatioDPMJETToFit, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      graphRatioDPMJETToFit->Draw("3,same");
      DrawGammaSetMarker(histoRatioDPMJETToFit, 24, 1.5, kOrange+2 , kOrange+2);
      histoRatioDPMJETToFit->SetLineWidth(widthCommonFit);
      histoRatioDPMJETToFit->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphRatioCombCombFitSys->SetLineWidth(1);
      graphRatioCombCombFitSys->Draw("2,same");
      graphRatioCombCombFitStat_WOXErr->Draw("p,same");

      boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0 , 1., 1.,1, kGray+2);

  canvasInvSectionPaper->Print(Form("%s/Pi0_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));

  padInvSectionSpec->cd();
  padInvSectionSpec->SetLogy(1);
  padInvSectionSpec->SetLogx(1);
      histo2DXSectionPi0->Draw();

      // graphPi0DSS14->Draw("same,e3");

      graphEPOSJJInvYield->Draw("3,same");
      histoEPOSJJInvYield->Draw("same,hist,l");

      graphDPMJETInvYield->Draw("3,same");
      histoDPMJETInvYield->Draw("same,hist,l");

      fitTCMInvYieldPi0->Draw("same");

      fitInvXSectionPi0->Draw("same");

      graphCombPi0InvYieldSys->Draw("E2same");
      graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");

      labelEnergyXSectionPaper->Draw();
      labelALICEXSectionPaper->Draw();
      labelDetSysXSectionPaper->Draw();
      
      DrawGammaNLOTGraph( graphNLOCalcPi0MuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      graphNLOCalcPi0MuHalf->Draw("same,c");
      DrawGammaNLOTGraph( graphNLOCalcPi0MuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
      graphNLOCalcPi0MuOne->Draw("same,c");
      DrawGammaNLOTGraph( graphNLOCalcPi0MuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      graphNLOCalcPi0MuTwo->Draw("same,c");

      TLegend* legendXsectionPaperLines     = GetAndSetLegend2(0.17, 0.03, 0.5, 0.13+0.05*7, textSizeLabelsPixel);
      legendXsectionPaperLines->SetNColumns(1);
      legendXsectionPaperLines->SetMargin(0.2);
      legendXsectionPaperLines->AddEntry(graphCombPi0InvYieldSys,"Data","pf");
      // legendXsectionPaperLines->AddEntry(boxErrorSigmaRatio, "norm. unc. 2.4%", "l");
      legendXsectionPaperLines->AddEntry(fitTCMInvYieldPi0,"TCM fit","l");
      legendXsectionPaperLines->AddEntry(fitInvXSectionPi0,"Levy-Tsallis fit","l");
      legendXsectionPaperLines->AddEntry(histoEPOSJJInvYield,"EPOS","l");
      legendXsectionPaperLines->AddEntry(histoDPMJETInvYield,"DPMJET","l");
      legendXsectionPaperLines->AddEntry((TObject*)0, "", "");
      legendXsectionPaperLines->AddEntry(graphNLOCalcPi0MuHalf,"#mu = 0.5 #it{p}_{T}","l");
      legendXsectionPaperLines->AddEntry(graphNLOCalcPi0MuOne,"#mu = #it{p}_{T}","l");
      legendXsectionPaperLines->AddEntry(graphNLOCalcPi0MuTwo,"#mu = 2 #it{p}_{T}","l");
      legendXsectionPaperLines->Draw();

      TLatex *labelPi0NLO = new TLatex(0.185, 0.21,"NLO, PDF:CT10 - FF:DSS14");
      SetStyleTLatex( labelPi0NLO, 0.75*textsizeLabelsPP,4);
      labelPi0NLO->Draw();


  padInvSectionNLORatio->cd();

      ratio2DNLO->DrawCopy();

//        TLegend* legendXsectionPaperPi02     = GetAndSetLegend2(0.17, 0.8, 0.4, 0.83+0.05*1, textSizeLabelsPixel*0.8);
//        legendXsectionPaperPi02->AddEntry(boxErrorSigmaRatio,"norm. unc. 3.5%","f");
//        legendXsectionPaperPi02->Draw();

//        graphRatioPi0CombNLODSS14->RemovePoint(0);
//        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
//        graphRatioPi0CombNLODSS14->Draw("3,same");
       DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
       graphRatioPi0CombNLOMuHalf->Draw("same,c");
       DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
       graphRatioPi0CombNLOMuOne->Draw("same,c");
       DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
       graphRatioPi0CombNLOMuTwo->Draw("same,c");
      // graphRatioPi0DSS14->SetLineWidth(widthCommonFit);
      // graphRatioPi0DSS14->SetLineColor(colorNLO);
      // graphRatioPi0DSS14->SetLineStyle(1);
      // graphRatioPi0DSS14->SetFillStyle(1001);
      // graphRatioPi0DSS14->SetFillColor(colorNLO);
      // graphRatioPi0DSS14->Draw("same,e4");

      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphRatioCombCombFitSys->SetLineWidth(1);
      graphRatioCombCombFitSys->Draw("2,same");
      graphRatioCombCombFitStat_WOXErr->Draw("p,same");

      boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0, 1., 1.,1, kGray+2);

  padInvSectionPythiaRatio->cd();
      ratio2DPythia->DrawCopy();

      graphRatioEPOSJJToFit->Draw("3,same");
      histoRatioEPOSJJToFit->Draw("same,hist,l");

      graphRatioDPMJETToFit->Draw("3,same");
      histoRatioDPMJETToFit->Draw("same,hist,l");

      graphRatioCombCombFitSys->Draw("2,same");
      graphRatioCombCombFitStat_WOXErr->Draw("p,same");

      boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0 , 1., 1.,1, kGray+2);

  canvasInvSectionPaper->Print(Form("%s/Pi0_InvYieldWithRatios_Lines_Paper.%s",outputDir.Data(),suffix.Data()));

  padInvSectionSpec->cd();
  padInvSectionSpec->SetLogy(1);
  padInvSectionSpec->SetLogx(1);
      histo2DXSectionPi0->Draw();

      graphEPOSJJInvYield->Draw("3,same");
      histoEPOSJJInvYield->Draw("same,hist,l");

      graphDPMJETInvYield->Draw("3,same");
      histoDPMJETInvYield->Draw("same,hist,l");

      fitTCMInvYieldPi0->Draw("same");

      fitInvXSectionPi0->Draw("same");

      graphCombPi0InvYieldSys->Draw("E2same");
      graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");

      labelEnergyXSectionPaper->Draw();
      labelALICEXSectionPaper->Draw();
      labelDetSysXSectionPaper->Draw();
      
      // DrawGammaNLOTGraph( graphNLOCalcPi0MuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      // graphNLOCalcPi0MuHalf->Draw("same,c");
      // DrawGammaNLOTGraph( graphNLOCalcPi0MuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
      // graphNLOCalcPi0MuOne->Draw("same,c");
      // DrawGammaNLOTGraph( graphNLOCalcPi0MuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      // graphNLOCalcPi0MuTwo->Draw("same,c");

      TLegend* legendXsectionPaperIndMeas     = GetAndSetLegend2(0.17, 0.03, 0.5, 0.13+0.05*3, textSizeLabelsPixel);
      legendXsectionPaperIndMeas->SetNColumns(1);
      legendXsectionPaperIndMeas->SetMargin(0.2);
      legendXsectionPaperIndMeas->AddEntry(graphCombPi0InvYieldSys,"Data","pf");
      // legendXsectionPaperIndMeas->AddEntry(boxErrorSigmaRatio, "norm. unc. 2.4%", "l");
      legendXsectionPaperIndMeas->AddEntry(fitTCMInvYieldPi0,"TCM fit","l");
      legendXsectionPaperIndMeas->AddEntry(fitInvXSectionPi0,"Levy-Tsallis fit","l");
      legendXsectionPaperIndMeas->AddEntry(histoEPOSJJInvYield,"EPOS","l");
      legendXsectionPaperIndMeas->AddEntry(histoDPMJETInvYield,"DPMJET","l");
      // legendXsectionPaperIndMeas->AddEntry((TObject*)0, "", "");
      // legendXsectionPaperIndMeas->AddEntry(graphNLOCalcPi0MuHalf,"#mu = 0.5 #it{p}_{T}","l");
      // legendXsectionPaperIndMeas->AddEntry(graphNLOCalcPi0MuOne,"#mu = #it{p}_{T}","l");
      // legendXsectionPaperIndMeas->AddEntry(graphNLOCalcPi0MuTwo,"#mu = 2 #it{p}_{T}","l");
      legendXsectionPaperIndMeas->Draw();


  padInvSectionNLORatio->cd();
        TH2F * ratio2Dindmeasppp               = new TH2F("ratio2Dindmeasppp","ratio2Dindmeasppp",1000,minPtPi0,maxPtPi0,1000,minRatio,maxRatio);
      SetStyleHistoTH2ForGraphs(ratio2Dindmeasppp, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{TCM fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
      ratio2Dindmeasppp->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2Dindmeasppp->GetYaxis()->SetNdivisions(505);
      ratio2Dindmeasppp->GetYaxis()->SetNoExponent(kTRUE);
      ratio2Dindmeasppp->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2Dindmeasppp->GetYaxis()->CenterTitle(kTRUE);
      ratio2Dindmeasppp->GetXaxis()->SetNoExponent(kTRUE);
      ratio2Dindmeasppp->GetXaxis()->SetLabelFont(42);
      ratio2Dindmeasppp->GetYaxis()->SetLabelFont(42);
      ratio2Dindmeasppp->GetYaxis()->SetLabelOffset(+0.01);
      ratio2Dindmeasppp->GetXaxis()->SetTickLength(0.07);
      ratio2Dindmeasppp->DrawCopy();
      ratio2Dindmeasppp->DrawCopy();

      TLegend* legendXsectionPaperIndMeas2     = GetAndSetLegend2(0.17, 0.8, 0.95, 0.95, textSizeLabelsPixel,4);
    for (Int_t i = 0; i < 13; i++){
      if(directoryPi0[i]){
        graphRatioCombFitSys[i]->Draw("E2same");
        graphRatioCombFitStat_WOXErr[i]->Draw("p,same,z");
        legendXsectionPaperIndMeas2->AddEntry(graphRatioCombFitSys[i],nameMeasGlobalPlot[i].Data(),"pf");
      }
    }

      // legendXsectionPaperIndMeas2->SetMargin(0.2);

      legendXsectionPaperIndMeas2->Draw();

      DrawGammaLines(minPtPi0, maxPtPi0, 1., 1.,1, kGray+2);

  padInvSectionPythiaRatio->cd();
      ratio2DPythia->DrawCopy();

      graphRatioEPOSJJToFit->Draw("3,same");
      histoRatioEPOSJJToFit->Draw("same,hist,l");

      graphRatioDPMJETToFit->Draw("3,same");
      histoRatioDPMJETToFit->Draw("same,hist,l");

      graphRatioCombCombFitSys->Draw("2,same");
      graphRatioCombCombFitStat_WOXErr->Draw("p,same");

      DrawGammaLines(minPtPi0, maxPtPi0 , 1., 1.,1, kGray+2);

  canvasInvSectionPaper->Print(Form("%s/Pi0_InvYieldWithRatios_incIndMeas_Paper.%s",outputDir.Data(),suffix.Data()));


  padInvSectionSpec->cd();
  padInvSectionSpec->SetLogy(1);
  padInvSectionSpec->SetLogx(1);
      SetStyleHistoTH2ForGraphs(histo2DXSectionEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                              0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
      histo2DXSectionEta->GetXaxis()->SetMoreLogLabels();
      histo2DXSectionEta->GetXaxis()->SetLabelOffset(+0.01);
      histo2DXSectionEta->GetYaxis()->SetRangeUser(minInvYieldEta,5E10);
      histo2DXSectionEta->Draw();

      DrawGammaNLOTGraph( graphNLOCalcEtaMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      graphNLOCalcEtaMuHalf->Draw("same,c");
      DrawGammaNLOTGraph( graphNLOCalcEtaMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
      graphNLOCalcEtaMuOne->Draw("same,c");
      DrawGammaNLOTGraph( graphNLOCalcEtaMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      graphNLOCalcEtaMuTwo->Draw("same,c");

      // DrawGammaSetMarkerTGraphErr(graphEPOSJJInvYieldEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphEPOSJJInvYieldEta->Draw("3,same");
      DrawGammaSetMarker(histoHIJINGEta, 24, 1.5, kGreen+2 , kGreen+2);
      histoHIJINGEta->SetLineWidth(widthCommonFit);
      histoHIJINGEta->Draw("same,hist,l");

      // DrawGammaSetMarkerTGraphErr(graphEPOSJJInvYieldEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphEPOSJJInvYieldEta->Draw("3,same");
      // DrawGammaSetMarker(histoEPOSJJInvYieldEta, 24, 1.5, kRed+2 , kRed+2);
      // histoEPOSJJInvYieldEta->SetLineWidth(widthCommonFit);
      // histoEPOSJJInvYieldEta->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphDPMJETInvYieldEta, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      graphDPMJETInvYieldEta->Draw("3,same");
      DrawGammaSetMarker(histoDPMJETInvYieldEta, 24, 1.5, kOrange+2 , kOrange+2);
      histoDPMJETInvYieldEta->SetLineWidth(widthCommonFit);
      histoDPMJETInvYieldEta->Draw("same,hist,l");

      DrawGammaSetMarkerTF1( fitTCMInvYieldEta, 7, 2, kGray+2);
      fitTCMInvYieldEta->Draw("same");

      DrawGammaSetMarkerTF1( fitInvXSectionEta, 3, 2, kGray+1);
      fitInvXSectionEta->Draw("same");

      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
      graphCombEtaInvYieldSys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
      graphCombEtaInvYieldStat_WOXErr->Draw("p,same,z");

      labelEnergyXSectionPaper->Draw();
      labelALICEXSectionPaper->Draw();
      TLatex *labelDetSysXSectionPaperEta = new TLatex(0.84,0.815,"#eta #rightarrow #gamma#gamma");
      SetStyleTLatex( labelDetSysXSectionPaperEta, textsizeLabelsXSecUp,4);
      labelDetSysXSectionPaperEta->Draw();

      TLegend* legendXsectionPaperEta     = GetAndSetLegend2(0.17, 0.03, 0.5, 0.13+0.05*8, textSizeLabelsPixel);
      legendXsectionPaperEta->SetNColumns(1);
      legendXsectionPaperEta->SetMargin(0.2);
      legendXsectionPaperEta->AddEntry(graphCombPi0InvYieldSys,"Data","pf");
      legendXsectionPaperEta->AddEntry(boxErrorSigmaRatio, "norm. unc. 2.4%", "l");
      legendXsectionPaperEta->AddEntry(fitTCMInvYieldEta,"TCM fit","l");
      legendXsectionPaperEta->AddEntry(fitInvXSectionEta,"Levy-Tsallis fit","l");
      legendXsectionPaperEta->AddEntry(histoHIJINGEta,"HIJING","l");
      // legendXsectionPaperEta->AddEntry(histoEPOSJJInvYieldEta,"EPOS","l");
      legendXsectionPaperEta->AddEntry(histoDPMJETInvYieldEta,"DPMJET","l");
      legendXsectionPaperEta->AddEntry((TObject*)0, "", "");
      legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuHalf,"#mu = 0.5 #it{p}_{T}","l");
      legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuOne,"#mu = #it{p}_{T}","l");
      legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuTwo,"#mu = 2 #it{p}_{T}","l");
      legendXsectionPaperEta->Draw();

      TLatex *labelEta = new TLatex(0.185, 0.192,"NLO, PDF:CTEQ6M5 - FF:AESSS");
      SetStyleTLatex( labelEta, 0.75*textsizeLabelsPP,4);
      labelEta->Draw();

  padInvSectionNLORatio->cd();
  padInvSectionNLORatio->SetLogx(1);
      TH2F * ratio2DNLOEta                = new TH2F("ratio2DNLOEta","ratio2DNLOEta",1000,minPtEta,maxPtEta,1000,0.25,2.9);
      SetStyleHistoTH2ForGraphs(ratio2DNLOEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{TCM fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
      ratio2DNLOEta->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2DNLOEta->GetYaxis()->SetNdivisions(505);
      ratio2DNLOEta->GetYaxis()->SetNoExponent(kTRUE);
      ratio2DNLOEta->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2DNLOEta->GetXaxis()->SetNoExponent(kTRUE);
      ratio2DNLOEta->GetXaxis()->SetLabelFont(42);
      ratio2DNLOEta->GetYaxis()->SetLabelFont(42);
      ratio2DNLOEta->GetYaxis()->SetLabelOffset(+0.01);
      ratio2DNLOEta->GetXaxis()->SetTickLength(0.07);
      ratio2DNLOEta->DrawCopy();

      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      graphRatioEtaCombNLOMuHalf->Draw("same,c");
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
      graphRatioEtaCombNLOMuOne->Draw("same,c");
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      graphRatioEtaCombNLOMuTwo->Draw("same,c");

      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatEta_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphRatioCombCombFitStatEta_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysEta, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphRatioCombCombFitSysEta->SetLineWidth(1);
      graphRatioCombCombFitSysEta->Draw("2,same");
      graphRatioCombCombFitStatEta_WOXErr->Draw("p,same");

      boxErrorSigmaRatioEta->Draw();
      DrawGammaLines(minPtEta, maxPtEta, 1., 1.,1, kGray+2);

  padInvSectionPythiaRatio->cd();
  padInvSectionPythiaRatio->SetLogx(1);
      TH2F * ratio2DPythiaEta             = new TH2F("ratio2DPythiaEta","ratio2DPythiaEta",1000,minPtEta,maxPtEta,1000,0.4,1.65);
      SetStyleHistoTH2ForGraphs(ratio2DPythiaEta, "#it{p}_{T} (GeV/#it{c})","#frac{Pythia, Data}{TCM fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 505);
      ratio2DPythiaEta->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2DPythiaEta->GetYaxis()->SetNdivisions(505);
      ratio2DPythiaEta->GetYaxis()->SetNoExponent(kTRUE);
      ratio2DPythiaEta->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2DPythiaEta->GetXaxis()->SetNoExponent(kTRUE);
      ratio2DPythiaEta->GetXaxis()->SetLabelFont(42);
      ratio2DPythiaEta->GetYaxis()->SetLabelFont(42);
      ratio2DPythiaEta->GetYaxis()->SetLabelOffset(+0.01);
      ratio2DPythiaEta->GetXaxis()->SetTickLength(0.06);
      ratio2DPythiaEta->GetYaxis()->SetTickLength(0.04);
      ratio2DPythiaEta->DrawCopy();

      // DrawGammaSetMarkerTGraphErr(graphRatioEPOSJJToFitEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphRatioEPOSJJToFitEta->Draw("3,same");
      DrawGammaSetMarker(histoRatioEtaHIJINGToFit, 24, 1.5, kGreen+2 , kGreen+2);
      histoRatioEtaHIJINGToFit->SetLineWidth(widthCommonFit);
      histoRatioEtaHIJINGToFit->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphRatioEPOSJJToFitEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      graphRatioEPOSJJToFitEta->Draw("3,same");
      DrawGammaSetMarker(histoRatioEPOSJJToFitEta, 24, 1.5, kRed+2 , kRed+2);
      histoRatioEPOSJJToFitEta->SetLineWidth(widthCommonFit);
      histoRatioEPOSJJToFitEta->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphRatioDPMJETToFitEta, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      graphRatioDPMJETToFitEta->Draw("3,same");
      DrawGammaSetMarker(histoRatioDPMJETToFitEta, 24, 1.5, kOrange+2 , kOrange+2);
      histoRatioDPMJETToFitEta->SetLineWidth(widthCommonFit);
      histoRatioDPMJETToFitEta->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatEta_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphRatioCombCombFitStatEta_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysEta, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphRatioCombCombFitSysEta->SetLineWidth(1);
      graphRatioCombCombFitSysEta->Draw("2,same");
      graphRatioCombCombFitStatEta_WOXErr->Draw("p,same");

      boxErrorSigmaRatioEta->Draw();
      DrawGammaLines(minPtEta, maxPtEta, 1., 1.,1, kGray+2);

  canvasInvSectionPaper->Print(Form("%s/Eta_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));


  //*************************************************************************************************************
  //***************************** PRL plot X-section and ratios ***********************************************
  //*************************************************************************************************************
  Int_t maxDivisions = 9;
  Double_t arrayBoundariesX1_XSecPRL[2];
  Double_t arrayBoundariesY1_XSecPRL[maxDivisions];
  Double_t relativeMarginsXXSecPRL[3];
  Double_t relativeMarginsYXSecPRL[3];
  textSizeLabelsPixel = 48;
  ReturnCorrectValuesForCanvasScaling(1250,300*maxDivisions, 1, maxDivisions-1,0.135, 0.005, 0.003,0.05,arrayBoundariesX1_XSecPRL,arrayBoundariesY1_XSecPRL,relativeMarginsXXSecPRL,relativeMarginsYXSecPRL);

  TCanvas* canvasInvSectionPRL      = new TCanvas("canvasInvSectionPRL","",0,0,1250,300*maxDivisions);  // gives the page size
  DrawGammaCanvasSettings( canvasInvSectionPRL,  0.13, 0.02, 0.03, 0.06);

  TPad* padPRLInvSectionSpec             = new TPad("padPRLInvSectionSpec", "", arrayBoundariesX1_XSecPRL[0], arrayBoundariesY1_XSecPRL[maxDivisions-5], arrayBoundariesX1_XSecPRL[1], arrayBoundariesY1_XSecPRL[0],-1, -1, -2);
  DrawGammaPadSettings( padPRLInvSectionSpec, relativeMarginsXXSecPRL[0], relativeMarginsXXSecPRL[2], relativeMarginsYXSecPRL[0], relativeMarginsYXSecPRL[1]);
  padPRLInvSectionSpec->Draw();
  marginXSec                 = relativeMarginsXXSecPRL[0]*1250;
  textsizeLabelsXSecUp       = 0;
  textsizeFacXSecUp          = 0;
  if (padPRLInvSectionSpec->XtoPixel(padPRLInvSectionSpec->GetX2()) < padPRLInvSectionSpec->YtoPixel(padPRLInvSectionSpec->GetY1())){
      textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padPRLInvSectionSpec->XtoPixel(padPRLInvSectionSpec->GetX2()) ;
      textsizeFacXSecUp               = (Double_t)1./padPRLInvSectionSpec->XtoPixel(padPRLInvSectionSpec->GetX2()) ;
  } else {
      textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padPRLInvSectionSpec->YtoPixel(padPRLInvSectionSpec->GetY1());
      textsizeFacXSecUp               = (Double_t)1./padPRLInvSectionSpec->YtoPixel(padPRLInvSectionSpec->GetY1());
  }

  TPad* padPRLInvMeaspPb         = new TPad("padPRLInvMeaspPb", "", arrayBoundariesX1_XSecPRL[0], arrayBoundariesY1_XSecPRL[maxDivisions-4], arrayBoundariesX1_XSecPRL[1], arrayBoundariesY1_XSecPRL[maxDivisions-5],-1, -1, -2);
  DrawGammaPadSettings( padPRLInvMeaspPb, relativeMarginsXXSecPRL[0], relativeMarginsXXSecPRL[2], relativeMarginsYXSecPRL[1], relativeMarginsYXSecPRL[1]);
  padPRLInvMeaspPb->Draw();
  textsizeLabelsXSecMiddle   = 0;
  textsizeFacXSecMiddle      = 0;
  if (padPRLInvMeaspPb->XtoPixel(padPRLInvMeaspPb->GetX2()) < padPRLInvMeaspPb->YtoPixel(padPRLInvMeaspPb->GetY1())){
      textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padPRLInvMeaspPb->XtoPixel(padPRLInvMeaspPb->GetX2()) ;
      textsizeFacXSecMiddle           = (Double_t)1./padPRLInvMeaspPb->XtoPixel(padPRLInvMeaspPb->GetX2()) ;
  } else {
      textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padPRLInvMeaspPb->YtoPixel(padPRLInvMeaspPb->GetY1());
      textsizeFacXSecMiddle           = (Double_t)1./padPRLInvMeaspPb->YtoPixel(padPRLInvMeaspPb->GetY1());
  }

  TPad* padPRLInvCombpPb      = new TPad("padPRLInvCombpPb", "", arrayBoundariesX1_XSecPRL[0], arrayBoundariesY1_XSecPRL[maxDivisions-3], arrayBoundariesX1_XSecPRL[1], arrayBoundariesY1_XSecPRL[maxDivisions-4],-1, -1, -2);
  DrawGammaPadSettings( padPRLInvCombpPb, relativeMarginsXXSecPRL[0], relativeMarginsXXSecPRL[2], relativeMarginsYXSecPRL[1], relativeMarginsYXSecPRL[1]);
  padPRLInvCombpPb->Draw();

  TPad* padPRLInvMeasppRef         = new TPad("padPRLInvMeasppRef", "", arrayBoundariesX1_XSecPRL[0], arrayBoundariesY1_XSecPRL[maxDivisions-2], arrayBoundariesX1_XSecPRL[1], arrayBoundariesY1_XSecPRL[maxDivisions-3],-1, -1, -2);
  DrawGammaPadSettings( padPRLInvMeasppRef, relativeMarginsXXSecPRL[0], relativeMarginsXXSecPRL[2], relativeMarginsYXSecPRL[1], relativeMarginsYXSecPRL[1]);
  padPRLInvMeasppRef->Draw();

  TPad* padPRLInvCombppRef      = new TPad("padPRLInvCombppRef", "", arrayBoundariesX1_XSecPRL[0], arrayBoundariesY1_XSecPRL[maxDivisions-1], arrayBoundariesX1_XSecPRL[1], arrayBoundariesY1_XSecPRL[maxDivisions-2],-1, -1, -2);
  DrawGammaPadSettings( padPRLInvCombppRef, relativeMarginsXXSecPRL[0], relativeMarginsXXSecPRL[2], relativeMarginsYXSecPRL[1], relativeMarginsYXSecPRL[2]);
  padPRLInvCombppRef->Draw();
  textsizeLabelsXSecDown     = 0;
  textsizeFacXSecDown        = 0;
  if (padPRLInvCombppRef->XtoPixel(padPRLInvCombppRef->GetX2()) < padPRLInvCombppRef->YtoPixel(padPRLInvCombppRef->GetY1())){
      textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padPRLInvCombppRef->XtoPixel(padPRLInvCombppRef->GetX2()) ;
      textsizeFacXSecDown             = (Double_t)1./padPRLInvCombppRef->XtoPixel(padPRLInvCombppRef->GetX2()) ;
  } else {
      textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padPRLInvCombppRef->YtoPixel(padPRLInvCombppRef->GetY1());
      textsizeFacXSecDown             = (Double_t)1./padPRLInvCombppRef->YtoPixel(padPRLInvCombppRef->GetY1());
  }


  padPRLInvSectionSpec->cd();
  padPRLInvSectionSpec->SetLogy(1);
  padPRLInvSectionSpec->SetLogx(1);
      SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                              0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
      histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
      histo2DXSectionPi0->GetXaxis()->SetLabelOffset(+0.01);
      histo2DXSectionPi0->GetYaxis()->SetRangeUser(minInvYieldPi0,6E11);
      histo2DXSectionPi0->Draw();

      graphPi0DSS14->SetLineWidth(widthCommonFit);
      graphPi0DSS14->SetLineColor(colorNLO);
      graphPi0DSS14->SetLineStyle(1);
      graphPi0DSS14->SetFillStyle(1001);
      graphPi0DSS14->SetFillColor(colorNLO);
      // graphPi0DSS14->Draw("same,e3");

      DrawGammaSetMarkerTGraphErr(graphEPOSJJInvYield, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphEPOSJJInvYield->Draw("3,same");
      DrawGammaSetMarker(histoEPOSJJInvYield, 24, 1.5, kRed+2 , kRed+2);
      histoEPOSJJInvYield->SetLineWidth(widthCommonFit);
      // histoEPOSJJInvYield->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphDPMJETInvYield, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      // graphDPMJETInvYield->Draw("3,same");
      DrawGammaSetMarker(histoDPMJETInvYield, 24, 1.5, kOrange+2 , kOrange+2);
      histoDPMJETInvYield->SetLineWidth(widthCommonFit);
      // histoDPMJETInvYield->Draw("same,hist,l");

      DrawGammaSetMarkerTF1( fitTCMInvYieldPi0, 7, 2, kGray+2);
      fitTCMInvYieldPi0->Draw("same");

      DrawGammaSetMarkerTF1( fitInvXSectionPi0, 3, 2, kGray+1);
      fitInvXSectionPi0->Draw("same");

      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
      graphCombPi0InvYieldSys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
      graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");

      DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0Sys, 24, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
      graphPPInvYieldCombPi0Sys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0Stat_woXErr, 24, markerSizeComb, kBlack, kBlack);
      graphPPInvYieldCombPi0Stat_woXErr->Draw("p,same,z");

      // TLatex *labelEnergyXSectionPaper= new TLatex(0.72, 0.91, collisionSystempPb8TeV.Data());
      labelEnergyXSectionPaper= new TLatex(0.95, 0.91, collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4);
      labelEnergyXSectionPaper->SetTextAlign(31);
      // labelEnergyXSectionPaper->Draw();
      // TLatex *labelALICEXSectionPaper= new TLatex(0.848,0.86,labelALICEPlot.Data());
      labelALICEXSectionPaper= new TLatex(0.95,0.91,labelALICEPlot.Data());
      SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4);
      labelALICEXSectionPaper->SetTextAlign(31);
      labelALICEXSectionPaper->Draw();
      // TLatex *labelDetSysXSectionPaper= new TLatex(0.824,0.815,"#pi^{0} #rightarrow #gamma#gamma");
      labelDetSysXSectionPaper= new TLatex(0.95,0.86,"#pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4);
      labelDetSysXSectionPaper->SetTextAlign(31);
      labelDetSysXSectionPaper->Draw();

      legendXsectionPaper    = GetAndSetLegend2(0.17, 0.08, 0.5, 0.18+0.05*5, textSizeLabelsPixel);
      legendXsectionPaper->SetNColumns(1);
      legendXsectionPaper->SetMargin(0.2);
      legendXsectionPaper->AddEntry(graphCombPi0InvYieldSys,"p-Pb, #sqrt{#it{s}_{NN}} = 8.16 TeV","pf");
      legendXsectionPaper->AddEntry(graphPPInvYieldCombPi0Sys,"pp, #sqrt{#it{s}} = 8 TeV","pf");
      // legendXsectionPaper->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.4%","l");
      legendXsectionPaper->AddEntry(fitTCMInvYieldPi0,"TCM fit","l");
      legendXsectionPaper->AddEntry(fitInvXSectionPi0,"Levy-Tsallis fit","l");
      legendXsectionPaper->AddEntry(histoEPOSJJInvYield,"EPOS","l");
      legendXsectionPaper->AddEntry(histoDPMJETInvYield,"DPMJET","l");
      legendXsectionPaper->AddEntry(graphPi0DSS14,  "NLO, PDF:CT10 - FF:DSS14", "f");
      legendXsectionPaper->Draw();

      labelRatioTheoryPP_Paper   = new TLatex(0.24,0.055,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
      SetStyleTLatex( labelRatioTheoryPP_Paper, 0.8*textsizeLabelsPP,4);
      labelRatioTheoryPP_Paper->Draw();

  padPRLInvMeaspPb->cd();
  padPRLInvMeaspPb->SetLogx(1);
  minRatio=0.51;
  maxRatio=1.65;
      ratio2Dindmeasppp               = new TH2F("ratio2Dindmeasppp","ratio2Dindmeasppp",1000,minPtPi0,maxPtPi0,1000,minRatio,maxRatio);
      SetStyleHistoTH2ForGraphs(ratio2Dindmeasppp, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{TCM fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
      ratio2Dindmeasppp->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2Dindmeasppp->GetYaxis()->SetNdivisions(505);
      ratio2Dindmeasppp->GetYaxis()->SetNoExponent(kTRUE);
      ratio2Dindmeasppp->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2Dindmeasppp->GetYaxis()->CenterTitle(kTRUE);
      ratio2Dindmeasppp->GetXaxis()->SetNoExponent(kTRUE);
      ratio2Dindmeasppp->GetXaxis()->SetLabelFont(42);
      ratio2Dindmeasppp->GetYaxis()->SetLabelFont(42);
      ratio2Dindmeasppp->GetYaxis()->SetLabelOffset(+0.01);
      ratio2Dindmeasppp->GetXaxis()->SetTickLength(0.07);
      ratio2Dindmeasppp->DrawCopy();

      TLegend* legendXsectionPaperIndMeas2PRL     = GetAndSetLegend2(0.17, 0.77, 0.95, 0.95, textSizeLabelsPixel,4);
    for (Int_t i = 0; i < 13; i++){
      if(directoryPi0[i]){
        graphRatioCombFitSys[i]->Draw("E2same");
        graphRatioCombFitStat_WOXErr[i]->Draw("p,same,z");
        legendXsectionPaperIndMeas2PRL->AddEntry(graphRatioCombFitSys[i],nameMeasGlobalPlot[i].Data(),"pf");
      }
    }

      legendXsectionPaperIndMeas2PRL->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0, 1., 1.,1, kGray+2);

  padPRLInvCombpPb->cd();
  padPRLInvCombpPb->SetLogx(1);
      ratio2DPythia            = new TH2F("ratio2DPythia","ratio2DPythia",1000,minPtPi0,maxPtPi0,1000,minRatio,maxRatio);
      SetStyleHistoTH2ForGraphs(ratio2DPythia, "#it{p}_{T} (GeV/#it{c})","#frac{Pythia, Data}{TCM fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 0.9,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
      ratio2DPythia->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2DPythia->GetYaxis()->SetNdivisions(505);
      ratio2DPythia->GetYaxis()->SetNoExponent(kTRUE);
      ratio2DPythia->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2DPythia->GetXaxis()->SetNoExponent(kTRUE);
      ratio2DPythia->GetXaxis()->SetLabelFont(42);
      ratio2DPythia->GetYaxis()->SetLabelFont(42);
      ratio2DPythia->GetYaxis()->CenterTitle(kTRUE);
      ratio2DPythia->GetYaxis()->SetLabelOffset(+0.01);
      ratio2DPythia->GetXaxis()->SetTickLength(0.06);
      ratio2DPythia->GetYaxis()->SetTickLength(0.04);
      ratio2DPythia->DrawCopy();

      //       graphRatioPi0DSS14->SetLineWidth(widthCommonFit);
      // graphRatioPi0DSS14->SetLineColor(colorNLO);
      // graphRatioPi0DSS14->SetLineStyle(1);
      // graphRatioPi0DSS14->SetFillStyle(1001);
      // graphRatioPi0DSS14->SetFillColor(colorNLO);
      // graphRatioPi0DSS14->Draw("same,e3");

      DrawGammaSetMarkerTGraphErr(graphRatioEPOSJJToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphRatioEPOSJJToFit->Draw("3,same");
      DrawGammaSetMarker(histoRatioEPOSJJToFit, 24, 1.5, kRed+2 , kRed+2);
      histoRatioEPOSJJToFit->SetLineWidth(widthCommonFit);
      // histoRatioEPOSJJToFit->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphRatioDPMJETToFit, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      // graphRatioDPMJETToFit->Draw("3,same");
      DrawGammaSetMarker(histoRatioDPMJETToFit, 24, 1.5, kOrange+2 , kOrange+2);
      histoRatioDPMJETToFit->SetLineWidth(widthCommonFit);
      // histoRatioDPMJETToFit->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphRatioCombCombFitSys->SetLineWidth(1);
      graphRatioCombCombFitSys->Draw("2,same");
      graphRatioCombCombFitStat_WOXErr->Draw("p,same");

      boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0 , 1., 1.,1, kGray+2);

  padPRLInvMeasppRef->cd();
  padPRLInvMeasppRef->SetLogx(1);
    
    ratio2Dindmeasppp->DrawCopy();

      // TLegend* legendXsectionPaperIndMeas2PRLpp     = GetAndSetLegend2(0.17, 0.8, 0.95, 0.95, textSizeLabelsPixel,4);
      TLegend* legendXsectionPaperIndMeas2PRLpp     = GetAndSetLegend2(0.23, 0.77, 0.93, 0.95, textSizeLabelsPixel,4);
    // for (Int_t i = 0; i < 11; i++){
    //   if(haveRefPPPi0[i]){
    //     DrawGammaSetMarkerTGraphAsym(systErrorCollectionPi0PPFitRatio[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
    //     DrawGammaSetMarkerTGraphAsym(statErrorCollectionPi0PPFitRatio_woXErr[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);
    //     systErrorCollectionPi0PPFitRatio[i]->Draw("E2same");
    //     statErrorCollectionPi0PPFitRatio_woXErr[i]->Draw("p,same,z");
    //     legendXsectionPaperIndMeas2PRLpp->AddEntry(systErrorCollectionPi0PPFitRatio[i],nameMeasGlobalPlot[i].Data(),"pf");
    //   }
    // }
      DrawGammaSetMarkerTGraphAsym(graphPPCombPi0StatFitRatioPub_woXErr, 24, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphPPCombPi0StatFitRatioPub_woXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphPPCombPi0FullSysPub, 24, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphPPCombPi0FullSysPub->SetLineWidth(1);
      graphPPCombPi0FullSysPub->Draw("2,same");
      graphPPCombPi0StatFitRatioPub_woXErr->Draw("p,same");
      legendXsectionPaperIndMeas2PRLpp->AddEntry(graphPPCombPi0FullSysPub,"EPJC 78 (2018) 263","pf");
      if(haveRefPPPi0[9]){
        DrawGammaSetMarkerTGraphAsym(systErrorCollectionPi0PPFitRatio[9], 30 ,markerSizeDet[9]*0.5, colorDet[9], colorDet[9], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(statErrorCollectionPi0PPFitRatio_woXErr[9], 30 ,markerSizeDet[9]*0.5, colorDet[9], colorDet[9]);
        systErrorCollectionPi0PPFitRatio[9]->Draw("E2same");
        statErrorCollectionPi0PPFitRatio_woXErr[9]->Draw("p,same,z");
        legendXsectionPaperIndMeas2PRLpp->AddEntry(systErrorCollectionPi0PPFitRatio[9],nameMeasGlobalPlot[9].Data(),"pf");
      }

      legendXsectionPaperIndMeas2PRLpp->Draw();

      DrawGammaLines(minPtPi0, maxPtPi0, 1., 1.,1, kGray+2);

  padPRLInvCombppRef->cd();
  padPRLInvCombppRef->SetLogx(1);
      ratio2DPythia            = new TH2F("ratio2DPythia","ratio2DPythia",1000,minPtPi0,maxPtPi0,1000,minRatio,maxRatio);
      SetStyleHistoTH2ForGraphs(ratio2DPythia, "#it{p}_{T} (GeV/#it{c})","#frac{Pythia, Data}{TCM fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 505);
      ratio2DPythia->GetYaxis()->SetMoreLogLabels(kTRUE);
      ratio2DPythia->GetYaxis()->SetNdivisions(505);
      ratio2DPythia->GetYaxis()->SetNoExponent(kTRUE);
      ratio2DPythia->GetXaxis()->SetMoreLogLabels(kTRUE);
      ratio2DPythia->GetXaxis()->SetNoExponent(kTRUE);
      ratio2DPythia->GetXaxis()->SetLabelFont(42);
      ratio2DPythia->GetYaxis()->SetLabelFont(42);
      ratio2DPythia->GetYaxis()->CenterTitle(kTRUE);
      ratio2DPythia->GetYaxis()->SetLabelOffset(+0.01);
      ratio2DPythia->GetXaxis()->SetTickLength(0.06);
      ratio2DPythia->GetYaxis()->SetTickLength(0.04);
      ratio2DPythia->DrawCopy();

      // DrawGammaSetMarkerTGraphErr(graphRatioEPOSJJToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      // graphRatioEPOSJJToFit->Draw("3,same");
      // DrawGammaSetMarker(histoRatioEPOSJJToFit, 24, 1.5, kRed+2 , kRed+2);
      // histoRatioEPOSJJToFit->SetLineWidth(widthCommonFit);
      // histoRatioEPOSJJToFit->Draw("same,hist,l");

      // DrawGammaSetMarkerTGraphErr(graphRatioDPMJETToFit, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      // graphRatioDPMJETToFit->Draw("3,same");
      // DrawGammaSetMarker(histoRatioDPMJETToFit, 24, 1.5, kOrange+2 , kOrange+2);
      // histoRatioDPMJETToFit->SetLineWidth(widthCommonFit);
      // histoRatioDPMJETToFit->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphAsym(graphPPCombPi0StatFitRatio_woXErr, 24, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphPPCombPi0StatFitRatio_woXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphPPCombPi0FullSysFitRatio, 24, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphPPCombPi0FullSysFitRatio->SetLineWidth(1);
      graphPPCombPi0FullSysFitRatio->Draw("2,same");
      graphPPCombPi0StatFitRatio_woXErr->Draw("p,same");

      boxErrorSigmaRatio->Draw();
      DrawGammaLines(minPtPi0, maxPtPi0 , 1., 1.,1, kGray+2);

  canvasInvSectionPRL->Print(Form("%s/Pi0_InvYieldWithRatios_PRL.%s",outputDir.Data(),suffix.Data()));

  // ***************************************************************************************************************
  // ******************************* pi0+eta combined plot  ********************************************************
  // ***************************************************************************************************************

  canvasXSectionPi0->cd();

  TH2F * histo2DXSectionWithEtaAndPi0;
  histo2DXSectionWithEtaAndPi0          = new TH2F("histo2DXSectionWithEtaAndPi0","histo2DXSectionWithEtaAndPi0",11000,0.9*minPtPi0,1.1*maxPtPi0,1000,minInvYieldPi0*0.1,maxInvYieldPi0*3.5);
  SetStyleHistoTH2ForGraphs(histo2DXSectionWithEtaAndPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.75);
  histo2DXSectionWithEtaAndPi0->GetXaxis()->SetMoreLogLabels();
  histo2DXSectionWithEtaAndPi0->GetXaxis()->SetNoExponent();
  histo2DXSectionWithEtaAndPi0->Draw("copy");

      // scale eta graphs
      Double_t scaleFacEtaForCombPlot                              = 1e-2;
      TGraphAsymmErrors* graphCombEtaInvYieldStat_WOXErrCopy   = (TGraphAsymmErrors*) graphCombEtaInvYieldStat_WOXErr->Clone("graphCombEtaInvYieldStatAWOXErrCopy");
      TGraphAsymmErrors* graphCombEtaInvYieldSys_Copy          = (TGraphAsymmErrors*) graphCombEtaInvYieldSys->Clone("graphCombEtaInvYieldSysA_Copy");
      graphCombEtaInvYieldStat_WOXErrCopy                      = ScaleGraph(graphCombEtaInvYieldStat_WOXErrCopy,scaleFacEtaForCombPlot);
      graphCombEtaInvYieldSys_Copy                             = ScaleGraph(graphCombEtaInvYieldSys_Copy,scaleFacEtaForCombPlot);

      TF1* fitTCMInvYieldEtaScaled                               = FitObject("tcm","fitTCMInvYieldEtaScaled","Eta",graphCombEtaInvYieldStat_WOXErrCopy,0.5,35.,paramTCMEta,"QNRMEX0+","", kFALSE);
      fitTCMInvYieldEtaScaled->SetRange(0.5,40);
      TF1* fitInvXSectionEtaScaled                               = FitObject("l","fitInvXSectionEtaScaled","Eta",graphCombEtaInvYieldStat_WOXErrCopy,0.5,18.,paramGraphEta,"QNRMEX0+");

      histoEPOSJJInvYieldEta->Scale(scaleFacEtaForCombPlot);
      TGraphErrors* graphEPOSJJInvYieldEtaScaled              = ScaleGraph(graphEPOSJJInvYieldEta,scaleFacEtaForCombPlot);

      histoDPMJETInvYieldEta->Scale(scaleFacEtaForCombPlot);
      TGraphErrors* graphDPMJETInvYieldEtaScaled              = ScaleGraph(graphDPMJETInvYieldEta,scaleFacEtaForCombPlot);

      TGraphAsymmErrors* graphEtaAESSSCopy                        = (TGraphAsymmErrors*)graphEtaAESSS->Clone("graphEtaAESSSCopy");
      graphEtaAESSSCopy                                           = ScaleGraph(graphEtaAESSSCopy,scaleFacEtaForCombPlot);

      // plotting NLO calcs pi0
      graphPi0DSS14->SetLineWidth(widthCommonFit);
      graphPi0DSS14->SetLineColor(colorNLO);
      graphPi0DSS14->SetLineStyle(1);
      graphPi0DSS14->SetFillStyle(1001);
      graphPi0DSS14->SetFillColor(colorNLO);
      graphPi0DSS14->Draw("same,e3");

      // plotting NLO calcs eta
      DrawGammaSetMarkerTGraphAsym(graphEtaAESSSCopy, 0, 0, colorCGC, colorCGC, widthLinesBoxes, kTRUE, colorCGC);
      graphEtaAESSSCopy->Draw("3,same");

      // plotting EPOS
      DrawGammaSetMarkerTGraphErr(graphEPOSJJInvYield, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      graphEPOSJJInvYield->Draw("3,same");
      DrawGammaSetMarker(histoEPOSJJInvYield, 24, 1.5, kRed+2 , kRed+2);
      histoEPOSJJInvYield->SetLineWidth(widthCommonFit);
      histoEPOSJJInvYield->Draw("same,hist,l");
      // plotting DPMJET
      DrawGammaSetMarkerTGraphErr(graphDPMJETInvYield, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      graphDPMJETInvYield->Draw("3,same");
      DrawGammaSetMarker(histoDPMJETInvYield, 24, 1.5, kOrange+2 , kOrange+2);
      histoDPMJETInvYield->SetLineWidth(widthCommonFit);
      histoDPMJETInvYield->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphPythia8EPPS16InvYield, 0, 0, kPink+2 , kPink+2, widthLinesBoxes, kTRUE, kPink+2);
      graphPythia8EPPS16InvYield->Draw("3,same");
      DrawGammaSetMarker(histoPythia8EPPS16InvYield, 24, 1.5, kPink+2 , kPink+2);
      histoPythia8EPPS16InvYield->SetLineWidth(widthCommonFit);
      histoPythia8EPPS16InvYield->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphPythia8nCTEQ15InvYield, 0, 0, kCyan+2 , kCyan+2, widthLinesBoxes, kTRUE, kCyan+2);
      graphPythia8nCTEQ15InvYield->Draw("3,same");
      DrawGammaSetMarker(histoPythia8nCTEQ15InvYield, 24, 1.5, kCyan+2 , kCyan+2);
      histoPythia8nCTEQ15InvYield->SetLineWidth(widthCommonFit);
      histoPythia8nCTEQ15InvYield->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphPythia8Monash2013InvSecpTh, 0, 0, kGreen+2 , kGreen+2, widthLinesBoxes, kTRUE, kGreen+2);
      // graphPythia8Monash2013InvSecpTh->Draw("3,same");
      DrawGammaSetMarker(histoPythia8Monash2013InvSecpTh, 24, 1.5, kGreen+2 , kGreen+2);
      histoPythia8Monash2013InvSecpTh->SetLineWidth(widthCommonFit);
      // histoPythia8Monash2013InvSecpTh->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphEPOSJJInvYieldEtaScaled, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      graphEPOSJJInvYieldEtaScaled->Draw("3,same");
      DrawGammaSetMarker(histoEPOSJJInvYieldEta, 24, 1.5, kRed+2 , kRed+2);
      histoEPOSJJInvYieldEta->SetLineWidth(widthCommonFit);
      histoEPOSJJInvYieldEta->Draw("same,hist,l");

      DrawGammaSetMarkerTGraphErr(graphDPMJETInvYieldEtaScaled, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
      graphDPMJETInvYieldEtaScaled->Draw("3,same");
      DrawGammaSetMarker(histoDPMJETInvYieldEta, 24, 1.5, kOrange+2 , kOrange+2);
      histoDPMJETInvYieldEta->SetLineWidth(widthCommonFit);
      histoDPMJETInvYieldEta->Draw("same,hist,l");

      // plots fits
      fitTCMInvYieldPi0->Draw("same");
      fitInvXSectionPi0->Draw("same");

      DrawGammaSetMarkerTF1(fitTCMInvYieldEtaScaled, 7, 2, kGray+2);
      fitTCMInvYieldEtaScaled->Draw("same");
      DrawGammaSetMarkerTF1( fitInvXSectionEtaScaled, 3, 2, kGray+1);
      fitInvXSectionEtaScaled->Draw("same");

      // plot data
      graphCombPi0InvYieldSys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys_Copy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
      graphCombEtaInvYieldSys_Copy->Draw("E2same");

      graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat_WOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
      graphCombEtaInvYieldStat_WOXErrCopy->Draw("p,same,z");

      // labels lower left corner
      TLegend* legendXsectionPaperAll    = GetAndSetLegend2(0.17, 0.12, 0.5, 0.11+0.04*4, textSizeLabelsPixel, 1, "", 43, 0.2);
      legendXsectionPaperAll->AddEntry(graphCombPi0InvYieldSys,"#pi^{0}","pf");
      legendXsectionPaperAll->AddEntry(graphCombEtaInvYieldSys_Copy,"#eta (x 10^{-2})","pf");
      legendXsectionPaperAll->AddEntry(fitTCMInvYieldPi0,"TCM fit","l");
      legendXsectionPaperAll->AddEntry(fitInvXSectionPi0,"Levy-Tsallis fit","l");
      legendXsectionPaperAll->Draw();

      TLatex *labelEnergyXSectionPaperAll = new TLatex(0.18, 0.12+0.04*6, collisionSystempPb8TeV.Data());
      SetStyleTLatex( labelEnergyXSectionPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
      labelEnergyXSectionPaperAll->Draw();
      TLatex *labelALICEXSectionPaperAll  = new TLatex(0.18,0.12+0.04*5,labelALICEPlot.Data());
      SetStyleTLatex( labelALICEXSectionPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
      labelALICEXSectionPaperAll->Draw();
      TLatex *labelALICENormUnPaperAll    = new TLatex(0.18,0.12+0.04*4+0.003,"norm. unc. 2.4%");
      SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
      labelALICENormUnPaperAll->Draw();

      // labels upper right corner
      TLegend* legendXsectionPaperPyBoth  = GetAndSetLegend2(0.5, 0.95-0.04*5, 0.54+0.33, 0.95, textSizeLabelsPixel, 1, "", 43, 0.18);
      legendXsectionPaperPyBoth->AddEntry(histoEPOSJJInvYieldEta,"EPOS","l");
      legendXsectionPaperPyBoth->AddEntry(histoDPMJETInvYieldEta,"DPMJET","l");
      legendXsectionPaperPyBoth->AddEntry(graphPi0DSS14,"#pi^{0} pQCD NLO ","f");
      legendXsectionPaperPyBoth->AddEntry((TObject*)0,"#scale[0.75]{PDF: CT10, FF: DSS14}","");
      legendXsectionPaperPyBoth->AddEntry(graphEtaAESSSCopy,"#eta pQCD NLO ","f");
      legendXsectionPaperPyBoth->AddEntry((TObject*)0,"#scale[0.75]{PDF: CTEQ6M5, FF: AESSS}","");
      legendXsectionPaperPyBoth->Draw();

  canvasXSectionPi0->SaveAs(Form("%s/InvYield_Pi0_Eta_Theory.%s",outputDir.Data(),suffix.Data()));

  histo2DXSectionWithEtaAndPi0->Draw("copy");

      // plotting EPOS
      graphEPOSJJInvYield->Draw("3,same");
      histoEPOSJJInvYield->Draw("same,hist,l");
      // plotting DPMJET
      graphDPMJETInvYield->Draw("3,same");
      histoDPMJETInvYield->Draw("same,hist,l");

      graphPythia8EPPS16InvYield->Draw("3,same");
      histoPythia8EPPS16InvYield->Draw("same,hist,l");

      graphPythia8nCTEQ15InvYield->Draw("3,same");
      histoPythia8nCTEQ15InvYield->Draw("same,hist,l");

      // graphPythia8Monash2013InvSecpTh->Draw("3,same");
      // histoPythia8Monash2013InvSecpTh->Draw("same,hist,l");

      graphEPOSJJInvYieldEtaScaled->Draw("3,same");
      histoEPOSJJInvYieldEta->Draw("same,hist,l");

      graphDPMJETInvYieldEtaScaled->Draw("3,same");
      histoDPMJETInvYieldEta->Draw("same,hist,l");

      // plots fits
      fitTCMInvYieldPi0->Draw("same");
      fitInvXSectionPi0->Draw("same");

      fitTCMInvYieldEtaScaled->Draw("same");
      fitInvXSectionEtaScaled->Draw("same");

      // plot data
      graphCombPi0InvYieldSys->Draw("E2same");
      graphCombEtaInvYieldSys_Copy->Draw("E2same");

      graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");
      graphCombEtaInvYieldStat_WOXErrCopy->Draw("p,same,z");

      // labels lower left corner
      legendXsectionPaperAll->Draw();

      labelEnergyXSectionPaperAll->Draw();
      labelALICEXSectionPaperAll->Draw();
      labelALICENormUnPaperAll->Draw();

      // labels upper right corner
      legendXsectionPaperPyBoth  = GetAndSetLegend2(0.5, 0.95-0.04*4, 0.54+0.33, 0.95, textSizeLabelsPixel, 1, "", 43, 0.18);
      legendXsectionPaperPyBoth->AddEntry(histoEPOSJJInvYieldEta,"EPOS","l");
      legendXsectionPaperPyBoth->AddEntry(histoDPMJETInvYieldEta,"DPMJET","l");
      legendXsectionPaperPyBoth->AddEntry(graphPythia8EPPS16InvYield,"PYTHIA8 EPPS16","l");
      legendXsectionPaperPyBoth->AddEntry(graphPythia8nCTEQ15InvYield,"PYTHIA8 nCTEQ15","l");
      // legendXsectionPaperPyBoth->AddEntry(graphPythia8Monash2013InvSecpTh,"PYTHIA8 Monash2013","l");
      legendXsectionPaperPyBoth->Draw();

  canvasXSectionPi0->SaveAs(Form("%s/InvYield_Pi0_Eta_Theory_noNLO.%s",outputDir.Data(),suffix.Data()));

  histo2DXSectionWithEtaAndPi0->Draw("copy");

      // plots fits
      fitTCMInvYieldPi0->Draw("same");
      fitInvXSectionPi0->Draw("same");

      DrawGammaSetMarkerTF1(fitTCMInvYieldEtaScaled, 7, 2, kGray+2);
      fitTCMInvYieldEtaScaled->Draw("same");
      DrawGammaSetMarkerTF1( fitInvXSectionEtaScaled, 3, 2, kGray+1);
      fitInvXSectionEtaScaled->Draw("same");

      // plot data
      graphCombPi0InvYieldSys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys_Copy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
      graphCombEtaInvYieldSys_Copy->Draw("E2same");

      graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat_WOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
      graphCombEtaInvYieldStat_WOXErrCopy->Draw("p,same,z");

      // labels lower left corner
      legendXsectionPaperAll->Draw();

      labelEnergyXSectionPaperAll->Draw();
      labelALICEXSectionPaperAll->Draw();
      labelALICENormUnPaperAll->Draw();

  canvasXSectionPi0->SaveAs(Form("%s/InvYield_Pi0_Eta_Fits.%s",outputDir.Data(),suffix.Data()));

  histo2DXSectionWithEtaAndPi0->Draw("copy");

      // plot data
      graphCombPi0InvYieldSys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys_Copy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
      graphCombEtaInvYieldSys_Copy->Draw("E2same");

      graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat_WOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
      graphCombEtaInvYieldStat_WOXErrCopy->Draw("p,same,z");

      // labels lower left corner
      TLegend* legendXsectionPaperAll2    = GetAndSetLegend2(0.17, 0.20, 0.5, 0.19+0.04*2, textSizeLabelsPixel, 1, "", 43, 0.2);
      legendXsectionPaperAll2->AddEntry(graphCombPi0InvYieldSys,"#pi^{0}","pf");
      legendXsectionPaperAll2->AddEntry(graphCombEtaInvYieldSys_Copy,"#eta (x 10^{-2})","pf");
      legendXsectionPaperAll2->Draw();

      labelEnergyXSectionPaperAll->Draw();
      labelALICEXSectionPaperAll->Draw();
      labelALICENormUnPaperAll->Draw();

  canvasXSectionPi0->SaveAs(Form("%s/InvYield_Pi0_Eta.%s",outputDir.Data(),suffix.Data()));

    histo2DXSectionWithEtaAndPi0->Draw("copy");
      // plot data
      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, kBlue+2, kBlue+2, widthLinesBoxes, kTRUE);
      graphCombPi0InvYieldSys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat_WOXErr, markerStyleComb, markerSizeComb, kBlue+2, kBlue+2);
      graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");
      Color_t colorReference = kGray+1;
      Color_t colorReference2 = kBlue-5;
      DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0Sys, 43, 1.6*markerSizeComb, colorReference2, colorReference2, widthLinesBoxes, kTRUE);
      graphPPInvYieldCombPi0Sys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0Stat_woXErr, 43, 1.6*markerSizeComb, colorReference2, colorReference2);
      graphPPInvYieldCombPi0Stat_woXErr->Draw("p,same,z");


      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys_Copy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
      graphCombEtaInvYieldSys_Copy->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat_WOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
      graphCombEtaInvYieldStat_WOXErrCopy->Draw("p,same,z");


      TGraphAsymmErrors* graphPPInvYieldCombEtaStat_woXErrCopy   = (TGraphAsymmErrors*) graphPPInvYieldCombEtaStat_woXErr->Clone("graphCombEtaInvYieldStatAWOXErrCopy");
      TGraphAsymmErrors* graphPPInvYieldCombEtaSys_Copy          = (TGraphAsymmErrors*) graphPPInvYieldCombEtaSys->Clone("graphCombEtaInvYieldSysA_Copy");
      graphPPInvYieldCombEtaStat_woXErrCopy                      = ScaleGraph(graphPPInvYieldCombEtaStat_woXErrCopy,scaleFacEtaForCombPlot);
      graphPPInvYieldCombEtaSys_Copy                             = ScaleGraph(graphPPInvYieldCombEtaSys_Copy,scaleFacEtaForCombPlot);
      DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombEtaSys_Copy, 42, 1.4*markerSizeComb, colorReference, colorReference, widthLinesBoxes, kTRUE);
      graphPPInvYieldCombEtaSys_Copy->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombEtaStat_woXErrCopy, 42, 1.4*markerSizeComb, colorReference, colorReference);
      graphPPInvYieldCombEtaStat_woXErrCopy->Draw("p,same,z");


      drawLatexAdd(collisionSystempPb8TeV.Data(),0.20, 0.11+0.04*5,0.8*textSizeLabelsRel,kFALSE);
      TLegend* legendXsectionPaperAll22    = GetAndSetLegend2(0.20, 0.11+0.04*3, 0.53, 0.10+0.04*5, textSizeLabelsPixel, 1, "", 43, 0.2);
      legendXsectionPaperAll22->AddEntry(graphCombPi0InvYieldSys,"#pi^{0}","pf");
      legendXsectionPaperAll22->AddEntry(graphCombEtaInvYieldSys_Copy,"#eta (x 10^{-2})","pf");
      legendXsectionPaperAll22->Draw();
      // labels lower left corner
      drawLatexAdd(collisionSystempp8TeV.Data(),0.20, 0.11+0.04*2,0.8*textSizeLabelsRel,kFALSE);
      legendXsectionPaperAll2    = GetAndSetLegend2(0.20, 0.11, 0.53, 0.10+0.04*2, textSizeLabelsPixel, 1, "", 43, 0.2);
      legendXsectionPaperAll2->AddEntry(graphPPInvYieldCombPi0Sys,"#pi^{0}","pf");
      legendXsectionPaperAll2->AddEntry(graphPPInvYieldCombEtaSys_Copy,"#eta (x 10^{-2})","pf");
      legendXsectionPaperAll2->Draw();



  drawLatexAdd("ALICE work-in-progress",0.95,0.92,0.8*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd("#pi^{0}, #eta #rightarrow #gamma#gamma",0.95,0.88,0.8*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      // labelEnergyXSectionPaperAll->Draw();
      // labelALICEXSectionPaperAll->Draw();
      // labelALICENormUnPaperAll->Draw();

  canvasXSectionPi0->SaveAs(Form("%s/InvYield_Pi0_Eta_incReference.%s",outputDir.Data(),suffix.Data()));

    histo2DXSectionWithEtaAndPi0->Draw("copy");


      // graphPythia8EPPS16InvYield->Draw("3,same");
      // histoPythia8EPPS16InvYield->Draw("same,hist,l");

      // graphPythia8nCTEQ15InvYield->Draw("3,same");
      // histoPythia8nCTEQ15InvYield->Draw("same,hist,l");

      // pPb pythia
      Color_t pythiaColor = kRed+2;
      DrawGammaSetMarkerTGraphErr(graphPythia8EPPS16InvYield, 0, 0, pythiaColor , pythiaColor, widthLinesBoxes, kTRUE, pythiaColor);
      graphPythia8EPPS16InvYield->Draw("3,same");
      DrawGammaSetMarker(histoPythia8EPPS16InvYield, 24, 1.5, pythiaColor , pythiaColor);
      histoPythia8EPPS16InvYield->SetLineWidth(widthCommonFit);
      histoPythia8EPPS16InvYield->Draw("same,hist,l");

      TGraphErrors* graphPythia8EPPS16InvYieldEta_Copy = (TGraphErrors*)graphPythia8EPPS16InvYieldEta->Clone("graphPythia8EPPS16InvYieldEta_Copy");
      graphPythia8EPPS16InvYieldEta_Copy = ScaleGraph(graphPythia8EPPS16InvYieldEta_Copy, 1e-2);
      DrawGammaSetMarkerTGraphErr(graphPythia8EPPS16InvYieldEta_Copy, 0, 0, pythiaColor , pythiaColor, widthLinesBoxes, kTRUE, pythiaColor);
      graphPythia8EPPS16InvYieldEta_Copy->Draw("3,same");

      TH1D* histoPythia8EPPS16InvYieldEta_Copy = (TH1D*)histoPythia8EPPS16InvYieldEta->Clone("histoPythia8EPPS16InvYieldEta_Copy");
      histoPythia8EPPS16InvYieldEta_Copy->Scale(1e-2);
      DrawGammaSetMarker(histoPythia8EPPS16InvYieldEta_Copy, 24, 1.5, pythiaColor , pythiaColor);
      histoPythia8EPPS16InvYieldEta_Copy->SetLineWidth(widthCommonFit);
      histoPythia8EPPS16InvYieldEta_Copy->Draw("same,hist,l");

      
      // reference pythia
      pythiaColor = kPink+2;
      DrawGammaSetMarkerTGraphErr(graphPythia8Monash2013InvXSection, 0, 0, pythiaColor , pythiaColor, widthLinesBoxes, kTRUE, pythiaColor);
      graphPythia8Monash2013InvXSection->Draw("3,same");
      DrawGammaSetMarker(histoPythia8Monash2013InvXSection, 24, 1.5, pythiaColor , pythiaColor);
      histoPythia8Monash2013InvXSection->SetLineWidth(widthCommonFit);
      histoPythia8Monash2013InvXSection->Draw("same,hist,l");

      TGraphErrors* graphPythia8Monash2013InvXSectionEta_Copy = (TGraphErrors*)graphPythia8Monash2013InvXSectionEta->Clone("graphPythia8Monash2013InvXSectionEta_Copy");
      graphPythia8Monash2013InvXSectionEta_Copy = ScaleGraph(graphPythia8Monash2013InvXSectionEta_Copy, 1e-2);
      DrawGammaSetMarkerTGraphErr(graphPythia8Monash2013InvXSectionEta_Copy, 0, 0, pythiaColor , pythiaColor, widthLinesBoxes, kTRUE, pythiaColor);
      graphPythia8Monash2013InvXSectionEta_Copy->Draw("3,same");

      TH1D* histoPythia8Monash2013InvXSectionEta_Copy = (TH1D*)histoPythia8Monash2013InvXSectionEta->Clone("histoPythia8Monash2013InvXSectionEta_Copy");
      histoPythia8Monash2013InvXSectionEta_Copy->Scale(1e-2);
      DrawGammaSetMarker(histoPythia8Monash2013InvXSectionEta_Copy, 24, 1.5, pythiaColor , pythiaColor);
      histoPythia8Monash2013InvXSectionEta_Copy->SetLineWidth(widthCommonFit);
      histoPythia8Monash2013InvXSectionEta_Copy->Draw("same,hist,l");


      // plot data
      graphCombPi0InvYieldSys->Draw("E2same");
      graphCombPi0InvYieldStat_WOXErr->Draw("p,same,z");
      graphPPInvYieldCombPi0Sys->Draw("E2same");
      graphPPInvYieldCombPi0Stat_woXErr->Draw("p,same,z");


      graphCombEtaInvYieldSys_Copy->Draw("E2same");
      graphCombEtaInvYieldStat_WOXErrCopy->Draw("p,same,z");
      graphPPInvYieldCombEtaSys_Copy->Draw("E2same");
      graphPPInvYieldCombEtaStat_woXErrCopy->Draw("p,same,z");

      fitTCMInvYieldPi0->Draw("same");
      TF1* fitPPInvYieldPi0                               = FitObject("tcm","fitPPInvYieldPi0","Pi0",graphPPInvYieldCombPi0Stat_woXErr,0.3,200.,paramTCMPi0New,"QNRMEX0+","", kFALSE);
      DrawGammaSetMarkerTF1(fitPPInvYieldPi0, 7, 2, kGray+1);
      fitPPInvYieldPi0->Draw("same");

      fitTCMInvYieldEtaScaled->Draw("same");
      TF1* fitTCMInvYieldEtaScaledPP                               = FitObject("tcm","fitTCMInvYieldEtaScaledPP","Eta",graphPPInvYieldCombEtaStat_woXErrCopy,0.5,35.,paramTCMEta,"QNRMEX0+","", kFALSE);
      DrawGammaSetMarkerTF1(fitTCMInvYieldEtaScaledPP, 7, 2, kGray+1);
      fitTCMInvYieldEtaScaledPP->SetRange(0.5, 40);
      fitTCMInvYieldEtaScaledPP->Draw("same");

      drawLatexAdd(collisionSystempp8TeV.Data(),0.20, 0.11+0.04*9,0.8*textSizeLabelsRel,kFALSE);
      legendXsectionPaperAll22    = GetAndSetLegend2(0.20, 0.11+0.04*5, 0.53, 0.10+0.04*9, textSizeLabelsPixel, 1, "", 43, 0.2);
      legendXsectionPaperAll22->AddEntry(graphPPInvYieldCombPi0Sys,"#pi^{0}","pf");
      legendXsectionPaperAll22->AddEntry(graphPPInvYieldCombEtaSys_Copy,"#eta (x 10^{-2})","pf");
      legendXsectionPaperAll22->AddEntry(fitPPInvYieldPi0,"TCM fit","l");
      legendXsectionPaperAll22->AddEntry(histoPythia8Monash2013InvXSection,"Pythia 8, Monash 2013","l");
      legendXsectionPaperAll22->Draw();
      // labels lower left corner
      drawLatexAdd(collisionSystempPb8TeV.Data(),0.20, 0.11+0.04*4,0.8*textSizeLabelsRel,kFALSE);
      legendXsectionPaperAll2    = GetAndSetLegend2(0.20, 0.11, 0.53, 0.10+0.04*4, textSizeLabelsPixel, 1, "", 43, 0.2);
      legendXsectionPaperAll2->AddEntry(graphCombPi0InvYieldSys,"#pi^{0}","pf");
      legendXsectionPaperAll2->AddEntry(graphCombEtaInvYieldSys_Copy,"#eta (x 10^{-2})","pf");
      legendXsectionPaperAll2->AddEntry(fitTCMInvYieldPi0,"TCM fit","l");
      legendXsectionPaperAll2->AddEntry(histoPythia8EPPS16InvYield,"Pythia 8, EPPS16","l");
      legendXsectionPaperAll2->Draw();



  drawLatexAdd("ALICE work-in-progress",0.95,0.92,0.8*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
  drawLatexAdd("#pi^{0}, #eta #rightarrow #gamma#gamma",0.95,0.88,0.8*textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

  canvasXSectionPi0->SaveAs(Form("%s/InvYield_Pi0_Eta_incReference_Theory.%s",outputDir.Data(),suffix.Data()));

  histoEPOSJJInvYieldEta->Scale(1/scaleFacEtaForCombPlot);
  histoDPMJETInvYieldEta->Scale(1/scaleFacEtaForCombPlot);

// ***************************************************************************************************************
// ******************************** fitting eta/pi0 **************************************************************
// ***************************************************************************************************************

TF1 *fitEtaToPi0 = new TF1("fitEtaToPi0","[0]",3.5,25.);
fitEtaToPi0->SetParameter(0,0.48);

TGraphAsymmErrors* comEtaPi0 = (TGraphAsymmErrors*) graphCombEtaToPi0Stat->Clone();
comEtaPi0->Fit(fitEtaToPi0,"QNRMEX0+","",3.5,25.);
cout << "\n\n\n\n\n++++++++++++++++++++++++++++++++" << endl;
cout << fitEtaToPi0->GetParameter(0) << ", +- " << fitEtaToPi0->GetParError(0) << endl;
cout << "++++++++++++++++++++++++++++++++\n\n\n\n\n" << endl;

fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
fLog << "Eta/Pi0 - Fit pol0: 3.5 < pT < 25.0" << endl;
fLog << fitEtaToPi0->GetParameter(0) << ", +- " << fitEtaToPi0->GetParError(0) << endl;

TGraphAsymmErrors* comEtaPi0Tot = (TGraphAsymmErrors*) graphCombEtaToPi0Tot->Clone();
comEtaPi0Tot->Fit(fitEtaToPi0,"QNRMEX0+","",3.5,25.);
cout << "\n\n\n\n\n++++++++++++++++++++++++++++++++" << endl;
cout << fitEtaToPi0->GetParameter(0) << ", +- " << fitEtaToPi0->GetParError(0) << endl;
cout << "++++++++++++++++++++++++++++++++\n\n\n\n\n" << endl;

fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
fLog << "Eta/Pi0 - Fit pol0: 3.5 < pT < 25.0" << endl;
fLog << fitEtaToPi0->GetParameter(0) << ", +- " << fitEtaToPi0->GetParError(0) << endl;


// ***************************************************************************************************************
// ******************************* eta/pi0 graphs without x-error  ***********************************************
// ***************************************************************************************************************

TGraphAsymmErrors* graphCombEtaToPi0Stat_WOXErr = (TGraphAsymmErrors*) graphCombEtaToPi0Stat->Clone("graphCombEtaToPi0Stat_WOXErr");
ProduceGraphAsymmWithoutXErrors(graphCombEtaToPi0Stat_WOXErr);

TGraphAsymmErrors* graphEtaToPi0Stat_WOXErr[11];
for (Int_t i = 0; i < 13; i++){
 if(directoryEta[i]){
   if(graphEtaToPi0Stat[i]){
     graphEtaToPi0Stat_WOXErr[i] = (TGraphAsymmErrors*) graphEtaToPi0Stat[i]->Clone(Form("graphEtaToPi0Stat_%i_WOXErr",i));
     ProduceGraphAsymmWithoutXErrors(graphEtaToPi0Stat_WOXErr[i]);
   }
 }
}

// ***************************************************************************************************************
// ******************************* Plotting eta/pi0 ratio for single measurements ********************************
// ***************************************************************************************************************
textSizeLabelsPixel                 = 54;
TCanvas* canvasEtatoPi0combo       = new TCanvas("canvasEtatoPi0combo","",200,10,1350,900);  // gives the page size
DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.1, 0.01, 0.01, 0.125);
canvasEtatoPi0combo->SetLogx();

   Double_t textsizeLabelsEtaToPi0 = 0;
   Double_t textsizeFacEtaToPi0    = 0;
   if (canvasEtatoPi0combo->XtoPixel(canvasEtatoPi0combo->GetX2()) <canvasEtatoPi0combo->YtoPixel(canvasEtatoPi0combo->GetY1()) ){
       textsizeLabelsEtaToPi0      = (Double_t)textSizeLabelsPixel/canvasEtatoPi0combo->XtoPixel(canvasEtatoPi0combo->GetX2()) ;
       textsizeFacEtaToPi0         = (Double_t)1./canvasEtatoPi0combo->XtoPixel(canvasEtatoPi0combo->GetX2()) ;
   } else {
       textsizeLabelsEtaToPi0      = (Double_t)textSizeLabelsPixel/canvasEtatoPi0combo->YtoPixel(canvasEtatoPi0combo->GetY1());
       textsizeFacEtaToPi0         = (Double_t)1./canvasEtatoPi0combo->YtoPixel(canvasEtatoPi0combo->GetY1());
   }
   cout << textsizeLabelsEtaToPi0 << endl;

   textsizeLabelsEtaToPi0*=0.9;

TH2F * histo2DEtatoPi0combo;
histo2DEtatoPi0combo               = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",1000,minPtEta,300,1000,-0.05,1.23    );
SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                         0.85*textsizeLabelsEtaToPi0,1.1*textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
histo2DEtatoPi0combo->GetXaxis()->SetNoExponent();
histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(minPtEta,maxPtEta);
histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(-0.05,0.89);
// histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(-0.05,1.15);
histo2DEtatoPi0combo->Draw("copy");

   // plotting systematics graphs
   for (Int_t i = 0; i < 13; i++){
       if(graphEtaToPi0Sys[i]){
           DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Sys[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
           graphEtaToPi0Sys[i]->Draw("E2same");
       }
   }
    DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Stat, markerStyleDet[1], markerSizeDet[1]*0.75, kPink , kPink, widthLinesBoxes, kTRUE);
   // plotting statistics graphs
   for (Int_t i = 0; i < 13; i++){
       if(graphEtaToPi0Stat[i]){
           DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Stat[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i]);
           graphEtaToPi0Stat[i]->Draw("p,same,e");
       }
   }
   DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys, markerStyleDet[1], markerSizeDet[1]*0.75, kPink , kPink);

   TLegend* legendEtaToPi0 = GetAndSetLegend2(0.67, 0.15, 0.9, 0.15+(textsizeLabelsEtaToPi0*4*0.9), textSizeLabelsPixel);
   for (Int_t i = 0; i < 13; i++){
       if(graphEtaToPi0Sys[i]){
           legendEtaToPi0->AddEntry(graphEtaToPi0Sys[i],nameMeasGlobalPlot[i],"pf");
       }
   }
   legendEtaToPi0->Draw();

   drawLatexAdd(collisionSystempPb8TeV.Data(),0.13, 0.92,0.85*textsizeLabelsEtaToPi0,kFALSE);
   drawLatexAdd(labelALICEPlot.Data(),0.13, 0.92-(1*textsizeLabelsEtaToPi0*0.85),0.85*textsizeLabelsEtaToPi0,kFALSE);

canvasEtatoPi0combo->Update();
canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_differentSystems.%s",outputDir.Data(), suffix.Data()));

//*************************************************************************************************************
//*************************************************************************************************************

// histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(minPtEta,32.);
histo2DEtatoPi0combo->Draw("copy");

TLegend* legendXsectionPaperEtaToPi0     = GetAndSetLegend2(0.12, 0.8, 0.45, 0.96, 0.85*textSizeLabelsPixel);
legendXsectionPaperEtaToPi0->SetNColumns(1);
legendXsectionPaperEtaToPi0->SetMargin(0.2);
legendXsectionPaperEtaToPi0->AddEntry(graphCombPi0InvYieldSys,"Data","pf");
legendXsectionPaperEtaToPi0->AddEntry(graphEtaToPi02760GeV,"ALICE pp, #sqrt{#it{s}} = 2.76 TeV","p");
legendXsectionPaperEtaToPi0->AddEntry(graphEtaToPi08000GeV,"ALICE pp, #sqrt{#it{s}} = 8 TeV","p");
legendXsectionPaperEtaToPi0->Draw();

DrawGammaSetMarkerTGraphAsym(graphEtaToPi08000GeV, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[1] , colorDet[1], widthLinesBoxes, kTRUE);
graphEtaToPi08000GeV->Draw("same,p");
DrawGammaSetMarkerTGraphAsym(graphEtaToPi02760GeV, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
graphEtaToPi02760GeV->Draw("same,p");

// plotting labels
// TLatex *labelEnergyEtaToPi02 = new TLatex(0.75, 0.92,collisionSystempPb8TeV.Data());
TLatex *labelEnergyEtaToPi02 = new TLatex(0.96, 0.92,collisionSystempPb8TeV.Data());
SetStyleTLatex( labelEnergyEtaToPi02, 0.85*textsizeLabelsEtaToPi0,4);
labelEnergyEtaToPi02->SetTextAlign(31);
labelEnergyEtaToPi02->Draw();

// TLatex *labelALICEEtaToPi02 = new TLatex(0.852, 0.92-(1*textsizeLabelsEtaToPi0*0.85),labelALICEPlot.Data());
TLatex *labelALICEEtaToPi02 = new TLatex(0.96, 0.92-(1*textsizeLabelsEtaToPi0),labelALICEPlot.Data());
SetStyleTLatex( labelALICEEtaToPi02, 0.85*textsizeLabelsEtaToPi0,4);
labelALICEEtaToPi02->SetTextAlign(31);
labelALICEEtaToPi02->Draw();

// plotting data
graphCombEtaToPi0Stat_WOXErr->Print();
DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Stat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
graphCombEtaToPi0Stat_WOXErr->SetLineWidth(widthLinesBoxes);
DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
graphCombEtaToPi0Sys->SetLineWidth(1);
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->Update();
canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Paper.%s",outputDir.Data(), suffix.Data()));

//*************************************************************************************************************
//*************************************************************************************************************

histo2DEtatoPi0combo->Draw("copy");

// eta/pi0 mt-scaled
TH1F *eta2pi0MtScaled = new TH1F("eta2pi0MtScaled","#eta/#pi^{0} from m_{T} scaling",5000,0.4,45.);
eta2pi0MtScaled->SetLineColor(kBlue+2);
eta2pi0MtScaled->SetLineWidth(2.);

TH1F *eta2pi0MtScaledTCM = new TH1F("eta2pi0MtScaledTCM","#eta/#pi^{0} from m_{T} scaling",5000,0.4,45.);
eta2pi0MtScaledTCM->SetLineColor(kBlue+2);
eta2pi0MtScaledTCM->SetLineWidth(2.);

Double_t eta2Pi0Const = 0.470685;
Double_t mPi0 = 0.134977;
Double_t mEta = 0.547853;
for (Int_t i=1; i<=eta2pi0MtScaled->GetNbinsX(); i++) {
 Double_t ptPi0 = eta2pi0MtScaled->GetBinCenter(i);
 if (ptPi0 < 0.3) continue;
 Double_t mtEta = TMath::Sqrt(mEta*mEta + ptPi0*ptPi0);
 Double_t ptEta = TMath::Sqrt(mtEta*mtEta - mPi0*mPi0);
 Double_t Reta2pi0 = fitInvXSectionPi0->Eval(ptEta) / fitInvXSectionPi0->Eval(ptPi0) * eta2Pi0Const;
 eta2pi0MtScaled->SetBinContent(i,Reta2pi0);

 Double_t Reta2pi0TCM = fitTCMInvYieldPi0->Eval(ptEta) / fitTCMInvYieldPi0->Eval(ptPi0) * eta2Pi0Const;
 eta2pi0MtScaledTCM->SetBinContent(i,Reta2pi0TCM);
}

TGraphAsymmErrors* graphRatioForMt_stat     = (TGraphAsymmErrors*)graphCombEtaToPi0Stat_WOXErr->Clone();
TGraphAsymmErrors* graphRatioForMt_sys      = (TGraphAsymmErrors*)graphCombEtaToPi0Sys->Clone();

Int_t n_stat                = graphRatioForMt_stat->GetN();
Double_t* yValue_stat       = graphRatioForMt_stat->GetY();
Double_t* yErrorLow_stat    = graphRatioForMt_stat->GetEYlow();
Double_t* yErrorHigh_stat   = graphRatioForMt_stat->GetEYhigh();
Double_t* yValue_sys        = graphRatioForMt_sys->GetY();
Double_t* yErrorLow_sys     = graphRatioForMt_sys->GetEYlow();
Double_t* yErrorHigh_sys    = graphRatioForMt_sys->GetEYhigh();
for (Int_t i = 0; i < n_stat; i++){
   Double_t ptPi0 = graphCombEtaToPi0Stat_WOXErr->GetX()[i];
   Double_t mtEta = TMath::Sqrt(mEta*mEta + ptPi0*ptPi0);
   Double_t ptEta = TMath::Sqrt(mtEta*mtEta - mPi0*mPi0);
   Double_t Reta2pi0 = fitInvXSectionPi0->Eval(ptEta) / fitInvXSectionPi0->Eval(ptPi0) * eta2Pi0Const;

   yValue_stat[i]       = yValue_stat[i] / Reta2pi0;
   yErrorLow_stat[i]    = yErrorLow_stat[i] / Reta2pi0;
   yErrorHigh_stat[i]   = yErrorHigh_stat[i] / Reta2pi0;

   yValue_sys[i]        = yValue_sys[i] / Reta2pi0;
   yErrorLow_sys[i]     = yErrorLow_sys[i] / Reta2pi0;
   yErrorHigh_sys[i]    = yErrorHigh_sys[i]/ Reta2pi0;
}

textSizeLabelsPixel = 48;
TLegend* legendXsectionPaperEtaToPi03     = GetAndSetLegend2(0.11, 0.86, 0.96, 0.96, 0.85*textSizeLabelsPixel);
legendXsectionPaperEtaToPi03->SetNColumns(2);
legendXsectionPaperEtaToPi03->SetMargin(0.15);
legendXsectionPaperEtaToPi03->AddEntry(graphCombPi0InvYieldSys,"ALICE pp, #sqrt{#it{s}} = 5 TeV","pf");
legendXsectionPaperEtaToPi03->AddEntry(graphEtaToPi02760GeV,"ALICE pp, #sqrt{#it{s}} = 2.76 TeV","p");
legendXsectionPaperEtaToPi03->AddEntry(graphEtaToPi08000GeV,"ALICE pp, #sqrt{#it{s}} = 8 TeV","p");
legendXsectionPaperEtaToPi03->AddEntry(eta2pi0MtScaled,"ALICE pp m_{T}-scaled, #sqrt{#it{s}} = 5.02 TeV","l");
legendXsectionPaperEtaToPi03->Draw();

DrawGammaSetMarkerTGraphAsym(graphEtaToPi08000GeV, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[1] , colorDet[1], widthLinesBoxes, kTRUE);
graphEtaToPi08000GeV->Draw("same,p");
DrawGammaSetMarkerTGraphAsym(graphEtaToPi02760GeV, 30, 3., kBlue-6, kBlue-6, 2., kFALSE);
graphEtaToPi02760GeV->Draw("same,p");

eta2pi0MtScaled->Draw("][ c same");

// plotting data
DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Stat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
graphCombEtaToPi0Stat_WOXErr->SetLineWidth(widthLinesBoxes);
DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
graphCombEtaToPi0Sys->SetLineWidth(1);
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->Update();
canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_mT.%s",outputDir.Data(), suffix.Data()));

//*************************************************************************************************************
//*************************************************************************************************************

histo2DEtatoPi0combo->Draw("copy");

TLegend* legendXsectionPaperEtaToPi03TCM     = GetAndSetLegend2(0.11, 0.86, 0.96, 0.96, 0.85*textSizeLabelsPixel);
legendXsectionPaperEtaToPi03TCM->SetNColumns(2);
legendXsectionPaperEtaToPi03TCM->SetMargin(0.15);
legendXsectionPaperEtaToPi03TCM->AddEntry(graphCombPi0InvYieldSys,"ALICE pp, #sqrt{#it{s}} = 5 TeV","pf");
legendXsectionPaperEtaToPi03TCM->AddEntry(graphEtaToPi02760GeV,"ALICE pp, #sqrt{#it{s}} = 2.76 TeV","p");
legendXsectionPaperEtaToPi03TCM->AddEntry(graphEtaToPi08000GeV,"ALICE pp, #sqrt{#it{s}} = 8 TeV","p");
legendXsectionPaperEtaToPi03TCM->AddEntry(eta2pi0MtScaledTCM,"ALICE pp m_{T}-scaled, #sqrt{#it{s}} = 5 TeV","l");
legendXsectionPaperEtaToPi03TCM->Draw();

graphEtaToPi08000GeV->Draw("same,p");
graphEtaToPi02760GeV->Draw("same,p");

eta2pi0MtScaledTCM->Draw("][ c same");

// plotting data
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->Update();
canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_mT_TCM.%s",outputDir.Data(), suffix.Data()));

//*************************************************************************************************************
//*************************************************************************************************************

TH2F * histo2DEtatoPi0ratio;
histo2DEtatoPi0ratio               = new TH2F("histo2DEtatoPi0ratio","histo2DEtatoPi0ratio",1000,minPtEta,maxPtEta,1000,-0.05,1.99    );
SetStyleHistoTH2ForGraphs(histo2DEtatoPi0ratio, "#it{p}_{T} (GeV/#it{c})","ratio", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                         0.85*textsizeLabelsEtaToPi0,textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
histo2DEtatoPi0ratio->GetXaxis()->SetMoreLogLabels();
histo2DEtatoPi0ratio->GetXaxis()->SetNoExponent(kTRUE);
histo2DEtatoPi0ratio->Draw("copy");

DrawGammaLines(minPtEta,maxPtEta , 1., 1.,1, kGray+2);

// plotting data
DrawGammaSetMarkerTGraphAsym(graphRatioForMt_stat, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
graphRatioForMt_stat->SetLineWidth(widthLinesBoxes);
DrawGammaSetMarkerTGraphAsym(graphRatioForMt_sys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
graphRatioForMt_sys->SetLineWidth(1);
graphRatioForMt_sys->Draw("2,same");
graphRatioForMt_stat->Draw("p,same");

TLegend* legendXsectionPaperEtaToPiRatio     = GetAndSetLegend2(0.13, 0.88, 0.96, 0.96, 0.85*textSizeLabelsPixel);
legendXsectionPaperEtaToPiRatio->SetNColumns(2);
legendXsectionPaperEtaToPiRatio->SetMargin(0.15);
legendXsectionPaperEtaToPiRatio->AddEntry(graphRatioForMt_sys,"ALICE pp, #sqrt{#it{s}} = 5 TeV - (#eta/#pi^{0})_{data}/(#eta/#pi^{0})_{m_{T}}","pf");
legendXsectionPaperEtaToPiRatio->Draw();


histo2DEtatoPi0ratio->Draw("axis,same");

canvasEtatoPi0combo->Update();
canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_mT_ratio.%s",outputDir.Data(), suffix.Data()));

//*************************************************************************************************************
//*************************************************************************************************************

histo2DEtatoPi0combo->Draw("copy");

textSizeLabelsPixel = 48;
TLegend* legendXsectionPaperEtaToPi03n     = GetAndSetLegend2(0.11, 0.81, 0.5, 0.96, 0.85*textSizeLabelsPixel);
legendXsectionPaperEtaToPi03n->SetNColumns(1);
legendXsectionPaperEtaToPi03n->SetMargin(0.15);
legendXsectionPaperEtaToPi03n->AddEntry(graphCombPi0InvYieldSys,"ALICE pp, #sqrt{#it{s}} = 5 TeV","pf");
legendXsectionPaperEtaToPi03n->AddEntry(graphEtaToPi02760GeV,"ALICE pp, #sqrt{#it{s}} = 2.76 TeV","p");
legendXsectionPaperEtaToPi03n->AddEntry(graphEtaToPi08000GeV,"ALICE pp, #sqrt{#it{s}} = 8 TeV","p");
legendXsectionPaperEtaToPi03n->Draw();

graphEtaToPi08000GeV->Draw("same,p");
graphEtaToPi02760GeV->Draw("same,p");

// plotting data
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->Update();
canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Comparison_no_mT.%s",outputDir.Data(), suffix.Data()));

//*************************************************************************************************************
//*************************************************************************************************************

histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(minPtEta,maxPtEta);
histo2DEtatoPi0combo->Draw("copy");

// DrawGammaSetMarkerTGraphErr(graphEPOSJJEtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
// graphEPOSJJEtaToPi0->Draw("3,same");
DrawGammaSetMarker(histoHIJINGEtaToPi0, 24, 1.5, kGreen+2 , kGreen+2);
histoHIJINGEtaToPi0->SetLineWidth(widthCommonFit);
histoHIJINGEtaToPi0->Draw("same,hist,l");

// DrawGammaSetMarkerTGraphErr(graphEPOSJJEtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
// graphEPOSJJEtaToPi0->Draw("3,same");
// DrawGammaSetMarker(histoEPOSJJEtaToPi0, 24, 1.5, kRed+2 , kRed+2);
// histoEPOSJJEtaToPi0->SetLineWidth(widthCommonFit);
// histoEPOSJJEtaToPi0->Draw("same,hist,l");

DrawGammaSetMarkerTGraphErr(graphDPMJETEtaToPi0, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
graphDPMJETEtaToPi0->Draw("3,same");
DrawGammaSetMarker(histoDPMJETEtaToPi0, 24, 1.5, kOrange+2 , kOrange+2);
histoDPMJETEtaToPi0->SetLineWidth(widthCommonFit);
histoDPMJETEtaToPi0->Draw("same,hist,l");

textSizeLabelsPixel = 48;
TLegend* legendXsectionPaperEtaToPi02     = GetAndSetLegend2(0.12, 0.69, 0.45, 0.69+0.045*6, 0.85*textSizeLabelsPixel);
legendXsectionPaperEtaToPi02->SetNColumns(1);
legendXsectionPaperEtaToPi02->SetMargin(0.2);
legendXsectionPaperEtaToPi02->AddEntry(graphCombPi0InvYieldSys,"Data","pf");
// legendXsectionPaperEtaToPi02->AddEntry(histoEPOSJJEtaToPi0,"EPOS","l");
legendXsectionPaperEtaToPi02->AddEntry(histoHIJINGEtaToPi0,"HIJING","l");
legendXsectionPaperEtaToPi02->AddEntry(histoDPMJETEtaToPi0,"DPMJET","l");
legendXsectionPaperEtaToPi02->AddEntry(graphNLOEtaToPi0,"NLO, PDF:CT10","f");
legendXsectionPaperEtaToPi02->AddEntry((TObject*)0,"#pi^{0} FF:DSS14, #eta FF:AESSS","");
legendXsectionPaperEtaToPi02->AddEntry((TObject*)0,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}","");
legendXsectionPaperEtaToPi02->Draw();

// plotting NLO
graphNLOEtaToPi0->SetLineWidth(widthCommonFit);
graphNLOEtaToPi0->SetLineColor(colorNLO);
graphNLOEtaToPi0->SetLineStyle(1);
graphNLOEtaToPi0->SetFillStyle(1001);
graphNLOEtaToPi0->SetFillColor(colorNLO);
graphNLOEtaToPi0->Draw("same,e4");

// plotting data
DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Stat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
graphCombEtaToPi0Stat_WOXErr->SetLineWidth(widthLinesBoxes);
DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
graphCombEtaToPi0Sys->SetLineWidth(1);
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

// plotting labels
labelEnergyEtaToPi02->Draw();
labelALICEEtaToPi02->Draw();

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->Update();
canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper.%s",outputDir.Data(), suffix.Data()));

//*************************************************************************************************************
//*************************************************************************************************************

histo2DEtatoPi0combo->Draw("copy");

// DrawGammaSetMarkerTGraphErr(graphEPOSJJEtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
// graphEPOSJJEtaToPi0->Draw("3,same");
DrawGammaSetMarker(histoHIJINGEtaToPi0, 24, 1.5, kGreen+2 , kGreen+2);
histoHIJINGEtaToPi0->SetLineWidth(widthCommonFit);
histoHIJINGEtaToPi0->Draw("same,hist,l");

// DrawGammaSetMarkerTGraphErr(graphEPOSJJEtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
// graphEPOSJJEtaToPi0->Draw("3,same");
// DrawGammaSetMarker(histoEPOSJJEtaToPi0, 24, 1.5, kRed+2 , kRed+2);
// histoEPOSJJEtaToPi0->SetLineWidth(widthCommonFit);
// histoEPOSJJEtaToPi0->Draw("same,hist,l");

DrawGammaSetMarkerTGraphErr(graphDPMJETEtaToPi0, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
graphDPMJETEtaToPi0->Draw("3,same");
DrawGammaSetMarker(histoDPMJETEtaToPi0, 24, 1.5, kOrange+2 , kOrange+2);
histoDPMJETEtaToPi0->SetLineWidth(widthCommonFit);
histoDPMJETEtaToPi0->Draw("same,hist,l");

TLegend* legendXsectionPaperEtaToPi05     = GetAndSetLegend2(0.12, 0.645, 0.45, 0.645+0.045*7, 0.85*textSizeLabelsPixel);
legendXsectionPaperEtaToPi05->SetNColumns(1);
legendXsectionPaperEtaToPi05->SetMargin(0.2);
legendXsectionPaperEtaToPi05->AddEntry(graphCombPi0InvYieldSys,"Data","pf");
legendXsectionPaperEtaToPi05->AddEntry(eta2pi0MtScaled,"ALICE pp, #sqrt{#it{s}} = 5 TeV from m_{T} scaling","l");
legendXsectionPaperEtaToPi05->AddEntry(histoHIJINGEtaToPi0,"HIJING","l");
// legendXsectionPaperEtaToPi05->AddEntry(histoEPOSJJEtaToPi0,"EPOS","l");
legendXsectionPaperEtaToPi05->AddEntry(histoDPMJETEtaToPi0,"DPMJET","l");
legendXsectionPaperEtaToPi05->AddEntry(graphNLOEtaToPi0,"NLO, PDF:CT10","f");
legendXsectionPaperEtaToPi05->AddEntry((TObject*)0,"#pi^{0} FF:DSS14, #eta FF:AESSS","");
legendXsectionPaperEtaToPi05->AddEntry((TObject*)0,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}","");
legendXsectionPaperEtaToPi05->Draw();

// plotting NLO
graphNLOEtaToPi0->Draw("same,e4");

// plotting data
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

// plotting labels
labelEnergyEtaToPi02->Draw();
labelALICEEtaToPi02->Draw();

eta2pi0MtScaled->Draw("][ c same");

// plotting data
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->Update();
canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper_mT.%s",outputDir.Data(), suffix.Data()));

//*************************************************************************************************************
//*************************************************************************************************************

histo2DEtatoPi0combo->Draw("copy");

// graphEPOSJJEtaToPi0->Draw("3,same");
// histoEPOSJJEtaToPi0->Draw("same,hist,l");
histoHIJINGEtaToPi0->Draw("same,hist,l");

graphDPMJETEtaToPi0->Draw("3,same");
histoDPMJETEtaToPi0->Draw("same,hist,l");

TLegend* legendXsectionPaperEtaToPi06     = GetAndSetLegend2(0.12, 0.645+0.045*3.5, 0.45, 0.645+0.045*7, 0.85*textSizeLabelsPixel);
legendXsectionPaperEtaToPi06->SetNColumns(1);
legendXsectionPaperEtaToPi06->SetMargin(0.2);
legendXsectionPaperEtaToPi06->AddEntry(graphCombPi0InvYieldSys,"Data","pf");
legendXsectionPaperEtaToPi06->AddEntry(eta2pi0MtScaled,"ALICE pp, #sqrt{#it{s}} = 5 TeV from m_{T} scaling","l");
// legendXsectionPaperEtaToPi06->AddEntry(histoEPOSJJEtaToPi0,"EPOS","l");
legendXsectionPaperEtaToPi06->AddEntry(histoHIJINGEtaToPi0,"HIJING","l");
legendXsectionPaperEtaToPi06->AddEntry(histoDPMJETEtaToPi0,"DPMJET","l");
legendXsectionPaperEtaToPi06->Draw();

// plotting data
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

// plotting labels
labelEnergyEtaToPi02->Draw();
labelALICEEtaToPi02->Draw();

eta2pi0MtScaled->Draw("][ c same");

// plotting data
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->Update();
canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Pythia_Paper_mT.%s",outputDir.Data(), suffix.Data()));

//*************************************************************************************************************
//*************************************************************************************************************

histo2DEtatoPi0combo->Draw("copy");

// plotting data
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

// plotting labels
labelEnergyEtaToPi02->Draw();
labelALICEEtaToPi02->Draw();

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Combined.%s",outputDir.Data(), suffix.Data()));
//*************************************************************************************************************
histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(minPtEta,250);

TGraphAsymmErrors* graphRatioCombCombFitStatEtaToPi0_WOXErr = (TGraphAsymmErrors*) graphRatioCombCombFitStatEtaToPi0->Clone("graphRatioCombCombFitStatEtaToPi0_WOXErr");
ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitStatEtaToPi0_WOXErr);
DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatEtaToPi0_WOXErr, 24, markerSizeComb, kOrange+2, kOrange+2, widthLinesBoxes, kFALSE);
graphRatioCombCombFitStatEtaToPi0_WOXErr->SetLineWidth(widthLinesBoxes);
DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysEtaToPi0, 24, markerSizeComb, kOrange+2, kOrange+2, widthLinesBoxes, kTRUE, 0);


histo2DEtatoPi0combo->Draw("copy");

// plotting data
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

graphRatioCombCombFitSysEtaToPi0->Draw("2,same");
graphRatioCombCombFitStatEtaToPi0_WOXErr->Draw("p,same");

// plotting labels
labelEnergyEtaToPi02->Draw();
labelALICEEtaToPi02->Draw();

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Combined_withFitRatio.%s",outputDir.Data(), suffix.Data()));

TGraphAsymmErrors* graphRatioEtaCombCombPi0FitStatEtaToPi0_WOXErr = (TGraphAsymmErrors*) graphRatioEtaCombCombPi0FitStatEtaToPi0->Clone("graphRatioEtaCombCombPi0FitStatEtaToPi0_WOXErr");
ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombPi0FitStatEtaToPi0_WOXErr);
DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombPi0FitStatEtaToPi0_WOXErr, 25, markerSizeComb, kMagenta+2, kMagenta+2, widthLinesBoxes, kFALSE);
graphRatioEtaCombCombPi0FitStatEtaToPi0_WOXErr->SetLineWidth(widthLinesBoxes);
DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombPi0FitFitSysEtaToPi0, 25, markerSizeComb, kMagenta+2, kMagenta+2, widthLinesBoxes, kTRUE, 0);


histo2DEtatoPi0combo->Draw("copy");

    TBox* boxHighpTetaToPi0Ratio = CreateBoxConv(kGray, 4.0, 0.459-TMath::Sqrt(TMath::Power(0.006,2)+TMath::Power(0.011,2)), 250, 0.459+TMath::Sqrt(TMath::Power(0.006,2)+TMath::Power(0.011,2)));
    boxHighpTetaToPi0Ratio->SetLineWidth(8);
    boxHighpTetaToPi0Ratio->Draw();
    DrawGammaLines(4.0,  250, 0.459,0.459,2,kBlack,7);

// plotting data
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

graphRatioCombCombFitSysEtaToPi0->Draw("2,same");
graphRatioCombCombFitStatEtaToPi0_WOXErr->Draw("p,same");

graphRatioEtaCombCombPi0FitFitSysEtaToPi0->Draw("2,same");
graphRatioEtaCombCombPi0FitStatEtaToPi0_WOXErr->Draw("p,same");

TF1* fitRatioEtaToPi0 = DivideTF1(fitTCMInvYieldEta, fitTCMInvYieldPi0,"fitRatioEtaToPi0");
DrawGammaSetMarkerTF1(fitRatioEtaToPi0, 7, 2, kBlue+2);
fitRatioEtaToPi0->SetRange(0.3,250);
cout << fitRatioEtaToPi0->GetExpFormula() << endl;
fitRatioEtaToPi0->Draw("same");
// plotting labels
labelEnergyEtaToPi02->Draw();
labelALICEEtaToPi02->Draw();

textSizeLabelsPixel = 48;
TLegend* legendXsectionPaperEtaToPi03nxx     = GetAndSetLegend2(0.13, 0.96-0.05*4, 0.53, 0.96, 0.85*textSizeLabelsPixel);
legendXsectionPaperEtaToPi03nxx->SetNColumns(1);
legendXsectionPaperEtaToPi03nxx->SetMargin(0.15);
legendXsectionPaperEtaToPi03nxx->AddEntry(graphCombEtaToPi0Sys,"combined spectrum","pf");
legendXsectionPaperEtaToPi03nxx->AddEntry(graphRatioCombCombFitSysEtaToPi0,"#eta fit / #pi^{0} data","pf");
legendXsectionPaperEtaToPi03nxx->AddEntry(graphRatioEtaCombCombPi0FitFitSysEtaToPi0,"#eta data / #pi^{0} fit","pf");
legendXsectionPaperEtaToPi03nxx->AddEntry(fitRatioEtaToPi0,"#eta fit / #pi^{0} fit","l");
legendXsectionPaperEtaToPi03nxx->Draw();

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Combined_withFitRatio_v2.%s",outputDir.Data(), suffix.Data()));
histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(minPtEta,maxPtEta);

//*************************************************************************************************************
//*************************************************************************************************************

canvasEtatoPi0combo->SetRightMargin(0.02);
histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(0.);
histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(0.,13.0);
histo2DEtatoPi0combo->Draw("copy");
legendXsectionPaperEtaToPi02->Draw();

//plotting MC
// DrawGammaSetMarkerTGraphErr(graphEPOSJJEtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
// graphEPOSJJEtaToPi0->Draw("3,same");
DrawGammaSetMarker(histoHIJINGEtaToPi0, 24, 1.5, kRed+2 , kRed+2);
histoHIJINGEtaToPi0->SetLineWidth(widthCommonFit);
histoHIJINGEtaToPi0->Draw("same,hist,l");
// DrawGammaSetMarkerTGraphErr(graphEPOSJJEtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
// graphEPOSJJEtaToPi0->Draw("3,same");
// DrawGammaSetMarker(histoEPOSJJEtaToPi0, 24, 1.5, kRed+2 , kRed+2);
// histoEPOSJJEtaToPi0->SetLineWidth(widthCommonFit);
// histoEPOSJJEtaToPi0->Draw("same,hist,l");
DrawGammaSetMarkerTGraphErr(graphDPMJETEtaToPi0, 0, 0, kOrange+2 , kOrange+2, widthLinesBoxes, kTRUE, kOrange+2);
graphDPMJETEtaToPi0->Draw("3,same");
DrawGammaSetMarker(histoDPMJETEtaToPi0, 24, 1.5, kOrange+2 , kOrange+2);
histoDPMJETEtaToPi0->SetLineWidth(widthCommonFit);
histoDPMJETEtaToPi0->Draw("same,hist,l");

// plotting NLO
graphNLOEtaToPi0->Draw("same,e4");

// plotting data
graphCombEtaToPi0Sys->Draw("2,same");
graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

// plotting labels
labelEnergyEtaToPi02->Draw();
labelALICEEtaToPi02->Draw();

histo2DEtatoPi0combo->Draw("axis,same");

canvasEtatoPi0combo->SetLogx(kFALSE);
canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper_LIN.%s",outputDir.Data(), suffix.Data()));




 //*****************************************************************************************************************
  // Plotting RpA s for different mesons
  //*****************************************************************************************************************
  textSizeLabelsPixel                         = 54;
  TCanvas* canvasRpPb = new TCanvas("canvasRpPb","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasRpPb,   0.08, 0.01, 0.01, 0.12);
  canvasRpPb->SetLogx();

  if (bWCorrection.Contains("Y")){
      // TH2F * histo2DRpPb  = new TH2F("histo2DRpPb","histo2DRpPb",1000,minPtPi0,maxPtPi0,1000,0.3,1.65);
      TH2F * histo2DRpPb  = new TH2F("histo2DRpPb","histo2DRpPb",1000,-10,300,1000,0.3,1.65);
      SetStyleHistoTH2ForGraphs(histo2DRpPb, "#it{p}_{T} (GeV/#it{c})","#it{R}_{pA}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                    0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
      histo2DRpPb->GetXaxis()->SetRangeUser(minPtPi0,maxPtPi0);
      histo2DRpPb->GetXaxis()->SetMoreLogLabels();
      histo2DRpPb->GetXaxis()->SetNoExponent(kTRUE);
      // histo2DRpPb->GetXaxis()->SetLabelOffset(-0.01);
      histo2DRpPb->DrawCopy();

        TBox* boxErrorNormRpAPi0          = CreateBoxConv(kGray+1, minPtPi0, 1.-(0.0472), minPtPi0*1.1, 1.+(0.0472));
        boxErrorNormRpAPi0->Draw();

        DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystPi0, markerStyleComb, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphRpPbCombSystPi0->Draw("E2same");
        // graphRpPbCombSystEta->RemovePoint(graphRpPbCombSystEta->GetN());
        // DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystEta, markerStyleComb+4, markerSizeComb, kGray+1 , kGray+1, widthLinesBoxes, kTRUE);
        // graphRpPbCombSystEta->Draw("E2same");

        DrawGammaLines(minPtPi0, maxPtPi0 , 1, 1 ,1, kGray, 7);


        TGraphAsymmErrors* graphRpPbCombStatPi0WOXErr   = (TGraphAsymmErrors*)graphRpPbCombStatPi0->Clone("graphRpPbCombStatPi0WOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatPi0WOXErr);
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatPi0WOXErr, markerStyleComb, markerSizeComb, kBlack , kBlack);
        graphRpPbCombStatPi0WOXErr->Draw("p,same,z");
        // TGraphAsymmErrors* graphRpPbCombStatEtaWOXErr   = (TGraphAsymmErrors*)graphRpPbCombStatEta->Clone("graphRpPbCombStatEtaWOXErr");
        // ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatEtaWOXErr);
        // graphRpPbCombStatEtaWOXErr->RemovePoint(graphRpPbCombStatEtaWOXErr->GetN());
        // DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatEtaWOXErr, markerStyleComb+4, markerSizeComb, kGray+1 , kGray+1);
        // graphRpPbCombStatEtaWOXErr->Draw("p,same,z");


        TLatex *labelEnergyRpPb     = new TLatex(0.12, 0.95-0.04*1, collisionSystempPb8TeV.Data());
        SetStyleTLatex( labelEnergyRpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelEnergyRpPb->Draw();
        TLatex *labelALICERpPb  = new TLatex(0.12,0.947-0.04*2.1,labelALICEPlot.Data());
        SetStyleTLatex( labelALICERpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICERpPb->Draw();

        TLegend* legendRpPbComb     = GetAndSetLegend2(0.12, 0.95-0.04*1.25*4, 0.32 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
        legendRpPbComb->AddEntry(graphRpPbCombSystPi0,"#pi^{0}","pf");
        // legendRpPbComb->AddEntry(graphRpPbCombSystEta,"#eta","pf");
        legendRpPbComb->Draw();

      histo2DRpPb->Draw("same,axis");
      canvasRpPb->Update();
      canvasRpPb->Print(Form("%s/Pi0_RpPb.%s",outputDir.Data(),suffix.Data()));



      //*****************************************************************************************************************
      // Plotting RpA pi0 diff rec techniques
      //*****************************************************************************************************************
      histo2DRpPb->GetYaxis()->SetRangeUser(0,1.8);
      histo2DRpPb->DrawCopy();

          boxErrorNormRpAPi0->Draw();
          for (Int_t i = 10; i> -1; i--){
              if (graphRpPbIndSystPi0[i]){
                // if( i == 0 ){
                //   // while(graphRpPbIndStatPi0WOXErr->GetX()[0] < 10) graphRpPbIndStatPi0WOXErr->RemovePoint(0);
                //   while(graphRpPbIndSystPi0[i]->GetX()[graphRpPbIndSystPi0[i]->GetN()-1] > 0.6) graphRpPbIndSystPi0[i]->RemovePoint(graphRpPbIndSystPi0[i]->GetN()-1);
                // }
                  DrawGammaSetMarkerTGraphAsym(graphRpPbIndSystPi0[i], markerStyleDet[i], markerSizeDet[i]*0.5, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
                  graphRpPbIndSystPi0[i]->Draw("E2same");
              }
          }
          TGraphAsymmErrors* graphRpPbIndStatPi0WOXErr[11]    = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

          for (Int_t i = 10; i> -1; i--){
              if (graphRpPbIndStatPi0[i]){
                // if( i == 0 ){
                //     while(graphRpPbIndStatPi0[i]->GetX()[graphRpPbIndStatPi0[i]->GetN()-1] > 0.6) graphRpPbIndStatPi0[i]->RemovePoint(graphRpPbIndStatPi0[i]->GetN()-1);
                // }
                  graphRpPbIndStatPi0WOXErr[i]                = (TGraphAsymmErrors*)graphRpPbIndStatPi0[i]->Clone(Form("graphRpPb%sStatPi0WOXErr", nameMeasGlobalPlot[i].Data()));
                  ProduceGraphAsymmWithoutXErrors(graphRpPbIndStatPi0WOXErr[i]);
                  DrawGammaSetMarkerTGraphAsym(graphRpPbIndStatPi0WOXErr[i], markerStyleDet[i], markerSizeDet[i]*0.5, colorDet[i] , colorDet[i]);
                  graphRpPbIndStatPi0WOXErr[i]->Draw("p,same,z");
              }
          }

          DrawGammaLines(minPtPi0, maxPtPi0 , 1, 1 ,1, kGray, 7);

          labelEnergyRpPb->Draw();
          labelALICERpPb->Draw();
          TLatex *labelPi0RpPb  = new TLatex(0.12,0.94-0.04*3,"#pi^{0} #rightarrow #gamma#gamma");
          SetStyleTLatex( labelPi0RpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
          labelPi0RpPb->Draw();

          // TLegend* legendRpPbInd     = GetAndSetLegend2(0.55, 0.15, 0.95 , 0.15+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
          TLegend* legendRpPbInd     = GetAndSetLegend2(0.12, 0.95-0.04*9, 0.37 , 0.95-0.04*4, textSizeLabelsPixel*0.85,1, "", 43, 0.2);
          for (Int_t i = 0; i< 11; i++){
              if (graphRpPbIndSystPi0[i])legendRpPbInd->AddEntry(graphRpPbIndSystPi0[i],nameMeasGlobalPlot[i].Data(),"pf");
          }
          legendRpPbInd->Draw();

      histo2DRpPb->Draw("same,axis");
      canvasRpPb->Update();
      canvasRpPb->Print(Form("%s/Pi0_RpPb_IndividualMeasurements.%s",outputDir.Data(),suffix.Data()));

  }


  if (bWCorrection.Contains("Y")){
      TH2F * histo2DRpPb  = new TH2F("histo2DRpPb","histo2DRpPb",1000,minPtEta,maxPtEta,1000,0.0,1.85);
      SetStyleHistoTH2ForGraphs(histo2DRpPb, "#it{p}_{T} (GeV/#it{c})","#it{R}_{pA}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                    0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
      histo2DRpPb->GetXaxis()->SetMoreLogLabels();
      histo2DRpPb->GetXaxis()->SetNoExponent(kTRUE);
      // histo2DRpPb->GetXaxis()->SetLabelOffset(-0.01);
      histo2DRpPb->DrawCopy();

        TBox* boxErrorNormRpAEta          = CreateBoxConv(kGray+1, minPtEta, 1.-(0.0472), minPtEta*1.1, 1.+(0.0472));
        boxErrorNormRpAEta->Draw();

        DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystEta, markerStyleComb, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphRpPbCombSystEta->Draw("E2same");
        // graphRpPbCombSystEta->RemovePoint(graphRpPbCombSystEta->GetN());
        // DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystEta, markerStyleComb+4, markerSizeComb, kGray+1 , kGray+1, widthLinesBoxes, kTRUE);
        // graphRpPbCombSystEta->Draw("E2same");

        DrawGammaLines(minPtEta, maxPtEta , 1, 1 ,1, kGray, 7);


        TGraphAsymmErrors* graphRpPbCombStatEtaWOXErr   = (TGraphAsymmErrors*)graphRpPbCombStatEta->Clone("graphRpPbCombStatEtaWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatEtaWOXErr);
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatEtaWOXErr, markerStyleComb, markerSizeComb, kBlack , kBlack);
        graphRpPbCombStatEtaWOXErr->Draw("p,same,z");
        // TGraphAsymmErrors* graphRpPbCombStatEtaWOXErr   = (TGraphAsymmErrors*)graphRpPbCombStatEta->Clone("graphRpPbCombStatEtaWOXErr");
        // ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatEtaWOXErr);
        // graphRpPbCombStatEtaWOXErr->RemovePoint(graphRpPbCombStatEtaWOXErr->GetN());
        // DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatEtaWOXErr, markerStyleComb+4, markerSizeComb, kGray+1 , kGray+1);
        // graphRpPbCombStatEtaWOXErr->Draw("p,same,z");


        TLatex *labelEnergyRpPb     = new TLatex(0.12, 0.95-0.04*1, collisionSystempPb8TeV.Data());
        SetStyleTLatex( labelEnergyRpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelEnergyRpPb->Draw();
        TLatex *labelALICERpPb  = new TLatex(0.12,0.947-0.04*2.1,labelALICEPlot.Data());
        SetStyleTLatex( labelALICERpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICERpPb->Draw();

        TLegend* legendRpPbComb     = GetAndSetLegend2(0.12, 0.95-0.04*1.25*4, 0.32 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
        // legendRpPbComb->AddEntry(graphRpPbCombSystEta,"#pi^{0}","pf");
        legendRpPbComb->AddEntry(graphRpPbCombSystEta,"#eta","pf");
        legendRpPbComb->Draw();

      histo2DRpPb->Draw("same,axis");
      canvasRpPb->Update();
      canvasRpPb->Print(Form("%s/Eta_RpPb.%s",outputDir.Data(),suffix.Data()));



      //*****************************************************************************************************************
      // Plotting RpA Eta diff rec techniques
      //*****************************************************************************************************************
      // histo2DRpPb->GetYaxis()->SetRangeUser(0,1.8);
      histo2DRpPb->DrawCopy();

          boxErrorNormRpAEta->Draw();
          for (Int_t i = 10; i> -1; i--){
              if (graphRpPbIndSystEta[i]){
                // if( i == 0 ){
                //   // while(graphRpPbIndStatEtaWOXErr->GetX()[0] < 10) graphRpPbIndStatEtaWOXErr->RemovePoint(0);
                //   while(graphRpPbIndSystEta[i]->GetX()[graphRpPbIndSystEta[i]->GetN()-1] > 0.6) graphRpPbIndSystEta[i]->RemovePoint(graphRpPbIndSystEta[i]->GetN()-1);
                // }
                  DrawGammaSetMarkerTGraphAsym(graphRpPbIndSystEta[i], markerStyleDet[i], markerSizeDet[i]*0.5, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
                  graphRpPbIndSystEta[i]->Draw("E2same");
              }
          }
          TGraphAsymmErrors* graphRpPbIndStatEtaWOXErr[11]    = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

          for (Int_t i = 10; i> -1; i--){
              if (graphRpPbIndStatEta[i]){
                // if( i == 0 ){
                //     while(graphRpPbIndStatEta[i]->GetX()[graphRpPbIndStatEta[i]->GetN()-1] > 0.6) graphRpPbIndStatEta[i]->RemovePoint(graphRpPbIndStatEta[i]->GetN()-1);
                // }
                  graphRpPbIndStatEtaWOXErr[i]                = (TGraphAsymmErrors*)graphRpPbIndStatEta[i]->Clone(Form("graphRpPb%sStatEtaWOXErr", nameMeasGlobalPlot[i].Data()));
                  ProduceGraphAsymmWithoutXErrors(graphRpPbIndStatEtaWOXErr[i]);
                  DrawGammaSetMarkerTGraphAsym(graphRpPbIndStatEtaWOXErr[i], markerStyleDet[i], markerSizeDet[i]*0.5, colorDet[i] , colorDet[i]);
                  graphRpPbIndStatEtaWOXErr[i]->Draw("p,same,z");
              }
          }

          DrawGammaLines(minPtEta, maxPtEta , 1, 1 ,1, kGray, 7);

          labelEnergyRpPb->Draw();
          labelALICERpPb->Draw();
          TLatex *labelEtaRpPb  = new TLatex(0.12,0.94-0.04*3,"#eta #rightarrow #gamma#gamma");
          SetStyleTLatex( labelEtaRpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
          labelEtaRpPb->Draw();

          // TLegend* legendRpPbInd     = GetAndSetLegend2(0.55, 0.15, 0.95 , 0.15+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
          TLegend* legendRpPbInd     = GetAndSetLegend2(0.12, 0.95-0.04*9, 0.37 , 0.95-0.04*4, textSizeLabelsPixel*0.85,1, "", 43, 0.2);
          for (Int_t i = 0; i< 11; i++){
              if (graphRpPbIndSystEta[i])legendRpPbInd->AddEntry(graphRpPbIndSystEta[i],nameMeasGlobalPlot[i].Data(),"pf");
          }
          legendRpPbInd->Draw();

      histo2DRpPb->Draw("same,axis");
      canvasRpPb->Update();
      canvasRpPb->Print(Form("%s/Eta_RpPb_IndividualMeasurements.%s",outputDir.Data(),suffix.Data()));

  }

  // Pi0 and Eta RpA plot
  if (bWCorrection.Contains("Y")){
      TH2F * histo2DRpPb  = new TH2F("histo2DRpPb","histo2DRpPb",1000,minPtPi0*0.9,maxPtPi0*1.1,1000,0.0,1.85);
      SetStyleHistoTH2ForGraphs(histo2DRpPb, "#it{p}_{T} (GeV/#it{c})","#it{R}_{pA}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                    0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
      histo2DRpPb->GetXaxis()->SetMoreLogLabels();
      histo2DRpPb->GetXaxis()->SetNoExponent(kTRUE);
      // histo2DRpPb->GetXaxis()->SetLabelOffset(-0.01);
      histo2DRpPb->DrawCopy();

        TBox* boxErrorNormRpAEta          = CreateBoxConv(kGray+1, 0.9*minPtPi0, 1.-(0.0472), 0.9*minPtPi0*1.1, 1.+(0.0472));
        boxErrorNormRpAEta->Draw();

        Color_t colorRpAPi0 = kBlue+2;
        TGraphAsymmErrors* graphRpPbCombStatPi0WOXErr   = (TGraphAsymmErrors*)graphRpPbCombStatPi0->Clone("graphRpPbCombStatPi0WOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatPi0WOXErr);
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatPi0WOXErr, markerStyleComb, markerSizeComb, colorRpAPi0 , colorRpAPi0);
        graphRpPbCombStatPi0WOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystPi0, markerStyleComb, markerSizeComb, colorRpAPi0 , colorRpAPi0, widthLinesBoxes, kTRUE);
        graphRpPbCombSystPi0->Draw("E2same");

        Color_t colorRpAEta = kOrange+2;
        TGraphAsymmErrors* graphRpPbCombStatEtaWOXErr   = (TGraphAsymmErrors*)graphRpPbCombStatEta->Clone("graphRpPbCombStatEtaWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatEtaWOXErr);
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatEtaWOXErr, 43, 1.5*markerSizeComb, colorRpAEta , colorRpAEta);
        graphRpPbCombStatEtaWOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystEta, 43, 1.5*markerSizeComb, colorRpAEta , colorRpAEta, widthLinesBoxes, kTRUE);
        graphRpPbCombSystEta->Draw("E2same");

        DrawGammaLines(minPtPi0*0.9,maxPtPi0*1.1 , 1, 1 ,1, kGray+1, 7);





        TLatex *labelEnergyRpPb     = new TLatex(0.12, 0.92, labelALICEPlot.Data());
        SetStyleTLatex( labelEnergyRpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelEnergyRpPb->Draw();
        TLatex *labelALICERpPb  = new TLatex(0.12,0.87,collisionSystempPb8TeV.Data());
        SetStyleTLatex( labelALICERpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICERpPb->Draw();

        TLegend* legendRpPbComb     = GetAndSetLegend2(0.12, 0.95-0.04*1.25*4, 0.32 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
        legendRpPbComb->AddEntry(graphRpPbCombSystPi0,"#pi^{0}","pf");
        legendRpPbComb->AddEntry(graphRpPbCombSystEta,"#eta","pf");
        legendRpPbComb->Draw();

      histo2DRpPb->Draw("same,axis");
      canvasRpPb->Update();
      canvasRpPb->Print(Form("%s/Pi0_Eta_RpPb.%s",outputDir.Data(),suffix.Data()));

  }

  // TFile* inputRpA5TeVFredi = new TFile("/home/nschmidt/Downloads/CombinedResultsPaperPPb5023GeV_2017_10_02.root");
  // TGraphAsymmErrors* graphPi0RpPbComb5TeVFrediStat = (TGraphAsymmErrors*) inputRpA5TeVFredi->Get("Pi0RpPb/CombinedPi0RpPbStatErr");
  // TGraphAsymmErrors* graphPi0RpPbComb5TeVFrediSyst = (TGraphAsymmErrors*) inputRpA5TeVFredi->Get("Pi0RpPb/CombinedPi0RpPbSystErr");
  // TGraphAsymmErrors* graphEtaRpPbComb5TeVFrediStat = (TGraphAsymmErrors*) inputRpA5TeVFredi->Get("EtaRpPb/CombinedEtaRpPbStatErr");
  // TGraphAsymmErrors* graphEtaRpPbComb5TeVFrediSyst = (TGraphAsymmErrors*) inputRpA5TeVFredi->Get("EtaRpPb/CombinedEtaRpPbSystErr");

  TFile* inputRpA5TeV = new TFile("ExternalInputpPb/CombNeutralMesons/CombinedResultsPaper_pPb_5023GeV_2018_05_25_Published.root");
  TGraphAsymmErrors* graphPi0RpPbComb5TeVStat = (TGraphAsymmErrors*) inputRpA5TeV->Get("Pi0RpPb/CombinedPi0RpPbStatErr");
  TGraphAsymmErrors* graphPi0RpPbComb5TeVSyst = (TGraphAsymmErrors*) inputRpA5TeV->Get("Pi0RpPb/CombinedPi0RpPbSystErr");
  TGraphAsymmErrors* graphEtaRpPbComb5TeVStat = (TGraphAsymmErrors*) inputRpA5TeV->Get("EtaRpPb/CombinedEtaRpPbStatErr");
  TGraphAsymmErrors* graphEtaRpPbComb5TeVSyst = (TGraphAsymmErrors*) inputRpA5TeV->Get("EtaRpPb/CombinedEtaRpPbSystErr");

  TFile* inputRpA5TeVCMS = new TFile("/home/nschmidt/Downloads/HEPData-ins1496050-v2_chargedPartRpA_CMS.root");
  TGraphAsymmErrors* graphPi0RpPbComb5TeVCMSTot = (TGraphAsymmErrors*) inputRpA5TeVCMS->Get("Table 16/Graph1D_y1");
  // TGraphAsymmErrors* graphPi0RpPbComb5TeVCMSSyst = (TGraphAsymmErrors*) inputRpA5TeV->Get("Pi0RpPb/CombinedPi0RpPbSystErr");

  TFile* inputRpA5TeVATLAS = new TFile("/home/nschmidt/Downloads/hepdata_ATLAS_chargedpartRpA/HEPData-ins1463284-v1-root.root");
  TGraphAsymmErrors* graphPi0RpPbComb5TeVATLASTot = (TGraphAsymmErrors*) inputRpA5TeVATLAS->Get("Table 5/Graph1D_y1");

  while(graphPi0RpPbComb5TeVATLASTot->GetX()[0] < 0.25) graphPi0RpPbComb5TeVATLASTot->RemovePoint(0);

  // Comparison to ppb5TeV
  if (bWCorrection.Contains("Y")){
      TH2F * histo2DRpPb  = new TH2F("histo2DRpPb","histo2DRpPb",1000,minPtPi0*0.9,maxPtPi0*1.1,1000,0.0,1.85);
      // TH2F * histo2DRpPb  = new TH2F("histo2DRpPb","histo2DRpPb",1000,-30,300,1000,0.0,2.85);
      SetStyleHistoTH2ForGraphs(histo2DRpPb, "#it{p}_{T} (GeV/#it{c})","#it{R}_{pA}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                    0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
      // histo2DRpPb->GetYaxis()->SetRangeUser(0.0,1.85);
      // histo2DRpPb->GetXaxis()->SetRangeUser(minPtPi0*0.9,maxPtPi0*1.1);
      histo2DRpPb->GetXaxis()->SetMoreLogLabels();
      histo2DRpPb->GetXaxis()->SetNoExponent(kTRUE);
      // histo2DRpPb->GetXaxis()->SetLabelOffset(-0.01);
      histo2DRpPb->DrawCopy();

        TBox* boxErrorNormRpAEta          = CreateBoxConv(kGray+1, 0.9*minPtPi0, 1.-(0.0472), 0.9*minPtPi0*1.1, 1.+(0.0472));
        boxErrorNormRpAEta->Draw();

        Color_t colorRpAPi0 = kBlue+2;
        TGraphAsymmErrors* graphRpPbCombStatPi0WOXErr   = (TGraphAsymmErrors*)graphRpPbCombStatPi0->Clone("graphRpPbCombStatPi0WOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatPi0WOXErr);
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatPi0WOXErr, markerStyleComb, markerSizeComb, colorRpAPi0 , colorRpAPi0);
        graphRpPbCombStatPi0WOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystPi0, markerStyleComb, markerSizeComb, colorRpAPi0 , colorRpAPi0, widthLinesBoxes, kTRUE);
        graphRpPbCombSystPi0->Draw("E2same");

        Color_t colorRpA5TeV = kGray+1;
        TGraphAsymmErrors* graphPi0RpPbComb5TeVStatWOXErr   = (TGraphAsymmErrors*)graphPi0RpPbComb5TeVStat->Clone("graphPi0RpPbComb5TeVStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphPi0RpPbComb5TeVStatWOXErr);
        DrawGammaSetMarkerTGraphAsym(graphPi0RpPbComb5TeVStatWOXErr, 30, 1.4*markerSizeComb, colorRpA5TeV , colorRpA5TeV);
        graphPi0RpPbComb5TeVStatWOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPi0RpPbComb5TeVSyst, 30, 1.4*markerSizeComb, colorRpA5TeV , colorRpA5TeV, widthLinesBoxes, kTRUE);
        graphPi0RpPbComb5TeVSyst->Draw("E2same");



        DrawGammaLines(minPtPi0*0.9,maxPtPi0*1.1 , 1, 1 ,1, kGray+1, 7);
        DrawGammaLines(minPtPi0*0.9,maxPtPi0*1.1 , 0.9, 0.9 ,1, kGray+1, 8);

        // TLatex *labelEnergyRpPb     = new TLatex(0.12, 0.92, Form("p-Pb, %s", labelALICEPlot.Data()));
        TLatex *labelEnergyRpPb     = new TLatex(0.12, 0.92, Form("ALICE, p-Pb, work-in-progress"));
        SetStyleTLatex( labelEnergyRpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelEnergyRpPb->Draw();
        TLatex *labelALICERpPb  = new TLatex(0.12,0.87,"#pi^{0} #rightarrow #gamma #gamma");
        SetStyleTLatex( labelALICERpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICERpPb->Draw();

        TLegend* legendRpPbComb     = GetAndSetLegend2(0.12, 0.93-0.04*1.25*4, 0.32 , 0.93-0.04*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
        legendRpPbComb->AddEntry(graphPi0RpPbComb5TeVSyst,Form("%s (EPJC 78 (2018) 624)",collisionSystempPb5TeV.Data()),"pf");
        legendRpPbComb->AddEntry(graphRpPbCombSystPi0,collisionSystempPb8TeV.Data(),"pf");
        legendRpPbComb->Draw();

      histo2DRpPb->Draw("same,axis");
      canvasRpPb->Update();
      canvasRpPb->Print(Form("%s/Pi0_RpPb_Comparison5TeV.%s",outputDir.Data(),suffix.Data()));



      histo2DRpPb->DrawCopy();
        boxErrorNormRpAEta->Draw();

        colorRpAPi0 = kBlue+2;
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatPi0WOXErr, markerStyleComb, markerSizeComb, colorRpAPi0 , colorRpAPi0);
        graphRpPbCombStatPi0WOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystPi0, markerStyleComb, markerSizeComb, colorRpAPi0 , colorRpAPi0, widthLinesBoxes, kTRUE);
        graphRpPbCombSystPi0->Draw("E2same");

        colorRpA5TeV = kGray+1;
        DrawGammaSetMarkerTGraphAsym(graphPi0RpPbComb5TeVStatWOXErr, 30, 1.4*markerSizeComb, colorRpA5TeV , colorRpA5TeV);
        graphPi0RpPbComb5TeVStatWOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPi0RpPbComb5TeVSyst, 30, 1.4*markerSizeComb, colorRpA5TeV , colorRpA5TeV, widthLinesBoxes, kTRUE);
        graphPi0RpPbComb5TeVSyst->Draw("E2same");

        Color_t colorRpA5TeVCMS = kRed+1;
        TGraphAsymmErrors* graphPi0RpPbComb5TeVCMSTotWOXErr   = (TGraphAsymmErrors*)graphPi0RpPbComb5TeVCMSTot->Clone("graphPi0RpPbComb5TeVCMSTotWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphPi0RpPbComb5TeVCMSTotWOXErr);
        DrawGammaSetMarkerTGraphAsym(graphPi0RpPbComb5TeVCMSTotWOXErr, 42, 1.4*markerSizeComb, colorRpA5TeVCMS , colorRpA5TeVCMS);
        graphPi0RpPbComb5TeVCMSTotWOXErr->Draw("p,same,z");



        DrawGammaLines(minPtPi0*0.9,maxPtPi0*1.1 , 1, 1 ,1, kGray+1, 7);
        labelEnergyRpPb->Draw();
        labelALICERpPb->Draw();

        legendRpPbComb     = GetAndSetLegend2(0.12, 0.93-0.04*1.25*4.2, 0.32 , 0.93-0.04*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
        legendRpPbComb->AddEntry(graphPi0RpPbComb5TeVSyst,Form("#sqrt{#it{s}_{NN}} = 5.02 TeV #scale[0.8]{(EPJC 78 (2018) 624)}"),"pf");
        // legendRpPbComb->AddEntry((TObject*)0,Form("#scale[0.8]{EPJC 78 (2018) 624}"),"");
        legendRpPbComb->AddEntry(graphRpPbCombSystPi0,Form("#sqrt{#it{s}_{NN}} = 8.16 TeV"),"pf");
        legendRpPbComb->Draw();

        TLegend * legendRpPbComb2     = GetAndSetLegend2(0.25, 0.31-0.04*1.25*3, 0.45 , 0.31-0.04*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
        legendRpPbComb2->AddEntry(graphPi0RpPbComb5TeVCMSTotWOXErr,Form("h^{#pm}, CMS, #sqrt{#it{s}_{NN}} = 5.02 TeV #scale[0.8]{(JHEP 1704 (2017) 039)}"),"pe");
        legendRpPbComb2->Draw();

        // TLatex *labelcollsys     = new TLatex(0.90, 0.92, "p-Pb");
        // SetStyleTLatex( labelcollsys, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 31);
        // labelcollsys->Draw();

      histo2DRpPb->Draw("same,axis");
      canvasRpPb->Update();
      canvasRpPb->Print(Form("%s/Pi0_RpPb_Comparison5TeV_CMS.%s",outputDir.Data(),suffix.Data()));


      histo2DRpPb->DrawCopy();
        boxErrorNormRpAEta->Draw();

        colorRpAPi0 = kBlue+2;
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatPi0WOXErr, markerStyleComb, markerSizeComb, colorRpAPi0 , colorRpAPi0);
        graphRpPbCombStatPi0WOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystPi0, markerStyleComb, markerSizeComb, colorRpAPi0 , colorRpAPi0, widthLinesBoxes, kTRUE);
        graphRpPbCombSystPi0->Draw("E2same");

        colorRpA5TeV = kGray+1;
        DrawGammaSetMarkerTGraphAsym(graphPi0RpPbComb5TeVStatWOXErr, 30, 1.4*markerSizeComb, colorRpA5TeV , colorRpA5TeV);
        graphPi0RpPbComb5TeVStatWOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPi0RpPbComb5TeVSyst, 30, 1.4*markerSizeComb, colorRpA5TeV , colorRpA5TeV, widthLinesBoxes, kTRUE);
        graphPi0RpPbComb5TeVSyst->Draw("E2same");

        colorRpA5TeVCMS = kRed+1;
        DrawGammaSetMarkerTGraphAsym(graphPi0RpPbComb5TeVCMSTotWOXErr, 42, 1.4*markerSizeComb, colorRpA5TeVCMS , colorRpA5TeVCMS);
        graphPi0RpPbComb5TeVCMSTotWOXErr->Draw("p,same,z");

        Color_t colorRpA5TeVATLAS = kOrange+1;
        TGraphAsymmErrors* graphPi0RpPbComb5TeVATLASTotWOXErr   = (TGraphAsymmErrors*)graphPi0RpPbComb5TeVATLASTot->Clone("graphPi0RpPbComb5TeVATLASTotWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphPi0RpPbComb5TeVATLASTotWOXErr);
        DrawGammaSetMarkerTGraphAsym(graphPi0RpPbComb5TeVATLASTotWOXErr, 28, 1.1*markerSizeComb, colorRpA5TeVATLAS , colorRpA5TeVATLAS);
        graphPi0RpPbComb5TeVATLASTotWOXErr->Draw("p,same,z");



        DrawGammaLines(minPtPi0*0.9,maxPtPi0*1.1 , 1, 1 ,1, kGray+1, 7);
        labelEnergyRpPb->Draw();
        labelALICERpPb->Draw();

        legendRpPbComb     = GetAndSetLegend2(0.12, 0.93-0.04*1.25*4.2, 0.32 , 0.93-0.04*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
        legendRpPbComb->AddEntry(graphPi0RpPbComb5TeVSyst,Form("#sqrt{#it{s}_{NN}} = 5.02 TeV #scale[0.8]{(EPJC 78 (2018) 624)}"),"pf");
        legendRpPbComb->AddEntry(graphRpPbCombSystPi0,Form("#sqrt{#it{s}_{NN}} = 8.16 TeV"),"pf");
        legendRpPbComb->Draw();

        legendRpPbComb2     = GetAndSetLegend2(0.25, 0.31-0.04*1.35*3, 0.45 , 0.31-0.04*2+0.04*1.25, textSizeLabelsPixel*0.85,1,"",43,0.3);
        legendRpPbComb2->AddEntry(graphPi0RpPbComb5TeVCMSTotWOXErr,Form("h^{#pm}, CMS, #sqrt{#it{s}_{NN}} = 5.02 TeV #scale[0.8]{(JHEP 1704 (2017) 039)}"),"pe");
        legendRpPbComb2->AddEntry(graphPi0RpPbComb5TeVATLASTotWOXErr,Form("h^{#pm}, ATLAS, #sqrt{#it{s}_{NN}} = 5.02 TeV #scale[0.8]{(PLB 763 (2016) 313)}"),"pe");
        legendRpPbComb2->Draw();

      histo2DRpPb->Draw("same,axis");
      canvasRpPb->Update();
      canvasRpPb->Print(Form("%s/Pi0_RpPb_Comparison5TeV_CMS_ATLAS.%s",outputDir.Data(),suffix.Data()));

      canvasRpPb->SetLogx(kFALSE);
      TH2F * histo2DRpPb_LinX  = new TH2F("histo2DRpPb_LinX","histo2DRpPb_LinX",1000,-30,300,1000,0.0,2.85);
      SetStyleHistoTH2ForGraphs(histo2DRpPb_LinX, "#it{p}_{T} (GeV/#it{c})","#it{R}_{pA}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                    0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
      histo2DRpPb_LinX->GetXaxis()->SetMoreLogLabels();
      histo2DRpPb_LinX->GetXaxis()->SetNoExponent(kTRUE);
      histo2DRpPb_LinX->GetYaxis()->SetRangeUser(0.1,1.99);
      histo2DRpPb_LinX->GetXaxis()->SetRangeUser(-8,209);
      histo2DRpPb_LinX->DrawCopy();
        boxErrorNormRpAEta->Draw();

        graphRpPbCombStatPi0WOXErr->Draw("p,same,z");
        graphRpPbCombSystPi0->Draw("E2same");

        graphPi0RpPbComb5TeVStatWOXErr->Draw("p,same,z");
        graphPi0RpPbComb5TeVSyst->Draw("E2same");

        graphPi0RpPbComb5TeVCMSTotWOXErr->Draw("p,same,z");

        graphPi0RpPbComb5TeVATLASTotWOXErr->Draw("p,same,z");



        DrawGammaLines(minPtPi0*0.9,maxPtPi0*1.1 , 1, 1 ,1, kGray+1, 7);
        labelEnergyRpPb->Draw();
        labelALICERpPb->Draw();

        legendRpPbComb->Draw();

        legendRpPbComb2->Draw();

      histo2DRpPb_LinX->Draw("same,axis");
      canvasRpPb->Update();
      canvasRpPb->Print(Form("%s/Pi0_RpPb_Comparison5TeV_CMS_ATLAS_linX.%s",outputDir.Data(),suffix.Data()));

  }








// **********************************************************************************************************************
// ************************* Saving of final results ********************************************************************
// **********************************************************************************************************************
cout << "#######################################################################"<< endl;
cout << "             creating output root file" << endl;
cout << "#######################################################################"<< endl;
  TString nameOutputCommonFile    = Form("%s/CombinedResultsPaperpPb8TeV_%s.root",outputDir.Data(), dateForOutput.Data());
  TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

  fCombResults.mkdir("Pi0pPb8TeV");
  TDirectoryFile* fileDirectoryPi0 = (TDirectoryFile*)fCombResults.Get("Pi0pPb8TeV");
  fCombResults.cd("Pi0pPb8TeV");
      if(graphCombPi0InvYieldTot)graphCombPi0InvYieldTot->Write("graphInvCrossSectionPi0CombpPb8TeV");
      if(graphCombPi0InvYieldStat)graphCombPi0InvYieldStat->Write("graphInvCrossSectionPi0CombpPb8TeVStatErr");
      if(graphCombPi0InvYieldSys)graphCombPi0InvYieldSys->Write("graphInvCrossSectionPi0CombpPb8TeVSysErr");
      if(graphCombPi0InvYieldTot)graphCombPi0InvYieldTot->Write("graphInvCrossSectionPi0CombpPb8TeVATotErr");
      if(graphCombPi0InvYieldStat)graphCombPi0InvYieldStat->Write("graphInvCrossSectionPi0CombpPb8TeVAStatErr");
      if(graphCombPi0InvYieldSys)graphCombPi0InvYieldSys->Write("graphInvCrossSectionPi0CombpPb8TeVASysErr");

      for (Int_t i = 0; i < 13; i++){
        if(directoryPi0[i]){
          if(graphPi0InvYieldStat[i])graphPi0InvYieldStat[i]          ->Write(Form("graphInvCrossSectionPi0%spPb8TeVStatErr",nameMeasGlobalPlot[i].Data()));
          if(graphPi0InvYieldSys[i])graphPi0InvYieldSys[i]           ->Write(Form("graphInvCrossSectionPi0%spPb8TeVSysErr",nameMeasGlobalPlot[i].Data()));
        }
      }

       // fits for pi0
      if(fitInvXSectionPi0)fitInvXSectionPi0->Write("TsallisFitPi0");
      if(fitTCMInvYieldPi0)fitTCMInvYieldPi0->Write("TwoComponentModelFitPi0");

      if (bWCorrection.Contains("Y")){
          if(graphCombPi0InvYieldTot_yShifted)graphCombPi0InvYieldTot_yShifted->Write("graphInvCrossSectionPi0CombpPb8TeV_yShifted");
          if(graphCombPi0InvYieldStat_yShifted)graphCombPi0InvYieldStat_yShifted->Write("graphInvCrossSectionPi0CombpPb8TeVStatErr_yShifted");
          if(graphCombPi0InvYieldSys_yShifted)graphCombPi0InvYieldSys_yShifted->Write("graphInvCrossSectionPi0CombpPb8TeVSysErr_yShifted");

          for (Int_t i = 0; i < 13; i++){
            if(directoryPi0[i]){
              if(graphPi0InvYieldStat_yShifted[i])graphPi0InvYieldStat_yShifted[i]          ->Write(Form("graphInvCrossSectionPi0%spPb8TeVStatErr_yShifted",nameMeasGlobalPlot[i].Data()));
              if(graphPi0InvYieldSys_yShifted[i])graphPi0InvYieldSys_yShifted[i]           ->Write(Form("graphInvCrossSectionPi0%spPb8TeVSysErr_yShifted",nameMeasGlobalPlot[i].Data()));
            }
          }

      }

      fileDirectoryPi0->mkdir("Supporting");
      fileDirectoryPi0->cd("Supporting");
      for (Int_t i = 0; i < 13; i++){
        if(directoryPi0[i]){
          // Writing full correction factors
          if(graphPi0AccTimesEff[i])graphPi0AccTimesEff[i]          ->Write(Form("Pi0CorrectionFactor%s",nameMeasGlobalPlot[i].Data()));
          if(graphPi0Mass[i])graphPi0Mass[i]                 ->Write(Form("Pi0MassData%s",nameMeasGlobalPlot[i].Data()));
          if(graphPi0TrueMass[i])graphPi0TrueMass[i]             ->Write(Form("Pi0MassMC%s",nameMeasGlobalPlot[i].Data()));
          if(graphPi0FWHMMeV[i])graphPi0FWHMMeV[i]              ->Write(Form("Pi0WidthData%s",nameMeasGlobalPlot[i].Data()));
          if(graphPi0TrueFWHMMeV[i])graphPi0TrueFWHMMeV[i]          ->Write(Form("Pi0WidthMC%s",nameMeasGlobalPlot[i].Data()));
        }
      }

  fCombResults.mkdir("EtapPb8TeV");
  TDirectoryFile* fileDirectoryEta = (TDirectoryFile*)fCombResults.Get("EtapPb8TeV");
  fCombResults.cd("EtapPb8TeV");
      if(graphCombEtaInvYieldTot)graphCombEtaInvYieldTot->Write("graphInvCrossSectionEtaCombpPb8TeV");
      if(graphCombEtaInvYieldStat)graphCombEtaInvYieldStat->Write("graphInvCrossSectionEtaCombpPb8TeVStatErr");
      if(graphCombEtaInvYieldSys)graphCombEtaInvYieldSys->Write("graphInvCrossSectionEtaCombpPb8TeVSysErr");
      if(graphCombEtaInvYieldTot)graphCombEtaInvYieldTot->Write("graphInvCrossSectionEtaCombpPb8TeVATotErr");
      if(graphCombEtaInvYieldStat)graphCombEtaInvYieldStat->Write("graphInvCrossSectionEtaCombpPb8TeVAStatErr");
      if(graphCombEtaInvYieldSys)graphCombEtaInvYieldSys->Write("graphInvCrossSectionEtaCombpPb8TeVASysErr");

      for (Int_t i = 0; i < 13; i++){
        if(directoryEta[i]){
          if(graphEtaInvYieldStat[i])graphEtaInvYieldStat[i]          ->Write(Form("graphInvCrossSectionEta%spPb8TeVStatErr",nameMeasGlobalPlot[i].Data()));
          if(graphEtaInvYieldSys[i])graphEtaInvYieldSys[i]           ->Write(Form("graphInvCrossSectionEta%spPb8TeVSysErr",nameMeasGlobalPlot[i].Data()));
        }
      }

       // fits for Eta
      if(fitInvXSectionEta)fitInvXSectionEta->Write("TsallisFitEta");
      if(fitTCMInvYieldEta)fitTCMInvYieldEta->Write("TwoComponentModelFitEta");

      if (bWCorrection.Contains("Y")){
          if(graphCombEtaInvYieldTot_yShifted)graphCombEtaInvYieldTot_yShifted->Write("graphInvCrossSectionEtaCombpPb8TeV_yShifted");
          if(graphCombEtaInvYieldStat_yShifted)graphCombEtaInvYieldStat_yShifted->Write("graphInvCrossSectionEtaCombpPb8TeVStatErr_yShifted");
          if(graphCombEtaInvYieldSys_yShifted)graphCombEtaInvYieldSys_yShifted->Write("graphInvCrossSectionEtaCombpPb8TeVSysErr_yShifted");

          for (Int_t i = 0; i < 13; i++){
            if(directoryEta[i]){
              if(graphEtaInvYieldStat_yShifted[i])graphEtaInvYieldStat_yShifted[i]          ->Write(Form("graphInvCrossSectionEta%spPb8TeVStatErr_yShifted",nameMeasGlobalPlot[i].Data()));
              if(graphEtaInvYieldSys_yShifted[i])graphEtaInvYieldSys_yShifted[i]           ->Write(Form("graphInvCrossSectionEta%spPb8TeVSysErr_yShifted",nameMeasGlobalPlot[i].Data()));
            }
          }

      }


      if(graphCombEtaToPi0Tot){
        while(graphCombEtaToPi0Tot->GetX()[graphCombEtaToPi0Tot->GetN()-1] >20) graphCombEtaToPi0Tot->RemovePoint(graphCombEtaToPi0Tot->GetN()-1);
        graphCombEtaToPi0Tot->Write("graphRatioEtaToPi0CombpPb8TeVTotErr");
      }
      if(graphCombEtaToPi0Stat){
        while(graphCombEtaToPi0Stat->GetX()[graphCombEtaToPi0Stat->GetN()-1] >20) graphCombEtaToPi0Stat->RemovePoint(graphCombEtaToPi0Stat->GetN()-1);
        graphCombEtaToPi0Stat->Write("graphRatioEtaToPi0CombpPb8TeVStatErr");
      }
      if(graphCombEtaToPi0Sys){
        while(graphCombEtaToPi0Sys->GetX()[graphCombEtaToPi0Sys->GetN()-1] >20) graphCombEtaToPi0Sys->RemovePoint(graphCombEtaToPi0Sys->GetN()-1);
        graphCombEtaToPi0Sys->Write("graphRatioEtaToPi0CombpPb8TeVSysErr");
      }

      for (Int_t i = 0; i < 13; i++){
        if(directoryEta[i]){
          if(graphEtaToPi0Stat[i])graphEtaToPi0Stat[i]          ->Write(Form("graphEtaToPi0%spPb8TeVStatErr",nameMeasGlobalPlot[i].Data()));
          if(graphEtaToPi0Sys[i])graphEtaToPi0Sys[i]           ->Write(Form("graphEtaToPi0%spPb8TeVSysErr",nameMeasGlobalPlot[i].Data()));
        }
      }

      if(graphRatioForMt_stat)graphRatioForMt_stat->Write("graphRatioMtScalingToEtaToPi0StatErr");
      if(graphRatioForMt_sys)graphRatioForMt_sys->Write("graphRatioMtScalingToEtaToPi0SysErr");

      fileDirectoryEta->mkdir("Supporting");
      fileDirectoryEta->cd("Supporting");
      for (Int_t i = 0; i < 13; i++){
        if(directoryEta[i]){
          // Writing full correction factors
          if(graphEtaAccTimesEff[i])graphEtaAccTimesEff[i]          ->Write(Form("EtaCorrectionFactor%s",nameMeasGlobalPlot[i].Data()));
          if(graphEtaMass[i])graphEtaMass[i]                 ->Write(Form("EtaMassData%s",nameMeasGlobalPlot[i].Data()));
          if(graphEtaTrueMass[i])graphEtaTrueMass[i]             ->Write(Form("EtaMassMC%s",nameMeasGlobalPlot[i].Data()));
          if(graphEtaFWHMMeV[i])graphEtaFWHMMeV[i]              ->Write(Form("EtaWidthData%s",nameMeasGlobalPlot[i].Data()));
          if(graphEtaTrueFWHMMeV[i])graphEtaTrueFWHMMeV[i]          ->Write(Form("EtaWidthMC%s",nameMeasGlobalPlot[i].Data()));
        }
      }
  fCombResults.Close();

// **********************************************************************************************************************
// ************************* Saving only fits to final results **********************************************************
// **********************************************************************************************************************

  TString nameOutputCommonFileFitsOnly    = Form("FitsPaperpPb8TeV_%s.root", dateForOutput.Data());
  TFile fFitsResults(nameOutputCommonFileFitsOnly.Data(), "RECREATE");

       // fits for pi0
      if(fitInvXSectionPi0)fitInvXSectionPi0->Write("TsallisFitPi0");
      if(fitTCMInvYieldPi0)fitTCMInvYieldPi0->Write("TwoComponentModelFitPi0");

      // fits for eta
      if(fitInvXSectionEta)fitInvXSectionEta->Write("TsallisFitEta");
      if(fitTCMInvYieldEta)fitTCMInvYieldEta->Write("TwoComponentModelFitEta");

  fFitsResults.Close();



  
}