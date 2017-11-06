/****************************************************************************************************************************
******        provided by Gamma Conversion Group                                                          *****
******        Nicolas Schmidt, n.schmidt@cern.ch                                                          *****
******        Daniel Mühlheim, d.muehlheim@cern.ch                                                         *****
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
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //    TString name;
};

void CombineMesonMeasurements7TeV_V3(   TString fileNamePCM          = "CombinationInput7TeV/data_PCMResultsFullCorrection_PP_7TeV_20170718.root",
                                        TString fileNamePHOS         = "/home/nschmidt/AnalysisResults/pp/7TeV/PHOS/pp7TeV_pass4_ppareek_PHOSResultsFullCorrection_10092017.root",
                                        //TString fileNameEMCAL        = "/home/nschmidt/AnalysisSoftware/CombinationInput7TeV/data_EMCAL-EMCALResultsFullCorrection_PP_20170714.root",
                                        TString fileNameEMCAL        = "/home/nschmidt/AnalysisResults/pp/7TeV/EMCal/Evi/mesonSpecrta7TeV_2011EMCAL_14Oct2017.root",
                                        TString fileNamePCMEMCAL     = "/home/nschmidt/AnalysisSoftware/CombinationInput7TeV/data_PCM-EMCALResultsFullCorrection_PP_20170714.root",
                                        TString fileNamePCMPHOS      = "CombinationInput7TeV/data_PCM-PHOSResultsFullCorrection_PP_NoBinShifting_v2.root",
                                        TString fileNameEMCAL2       = "/home/nschmidt/AnalysisSoftware/CombinationInput7TeV/data_EMCAL-EMCALResultsFullCorrection_PP_20170714.root",
                                        TString fileInputCorrFactors = "/home/nschmidt/AnalysisResults/pp/7TeV/Comb/correlationInput/ComputeCorrelationFactors_pp7TeV/pp7TeV.root",
                                        TString suffix               = "pdf",
                                        TString isMC                 = "",
                                        TString thesisPlots          = "",
                                        TString bWCorrection         = "",
                                        Int_t numbersofmeas          = 5,
                                        Bool_t useDanielmeas         = kFALSE
                                    ){

    TString date = ReturnDateString();

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;

    //___________________________________ Declaration of files _____________________________________________
    TString collisionSystem7TeV                 = "pp, #sqrt{#it{s}} = 7 TeV";

    TString fileNameTheory                      = "ExternalInput/Theory/TheoryCompilationPP.root";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements7TeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
    cout << outputDir.Data() << endl;

    cout << "------------------------------------------" << endl;
    cout << "input files: " << endl;
    cout << fileNamePCM.Data() << endl;
    cout << fileNamePHOS.Data() << endl;
    cout << fileNameEMCAL.Data() << endl;
    cout << fileNamePCMEMCAL.Data() << endl;
fileNamePCMPHOS="";
    cout << fileNamePCMPHOS.Data() << endl;
fileNameEMCAL2="";
    cout << fileNameEMCAL2.Data() << endl;
    cout << "------------------------------------------" << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMPHOS.root", fileNamePCMPHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEMC.root", fileNamePCMEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCAL.root", fileNameEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/Theory.root", fileNameTheory.Data(), outputDir.Data()));

    Bool_t thesis                               = kFALSE;
    if(thesisPlots.CompareTo("thesis") == 0){
        thesis                                  = kTRUE;
    }

    TString prefix2                             = "";
    if (isMC.CompareTo("kTRUE")==0){
        prefix2                                 = "MC";
    } else {
        prefix2                                 = "Data";
    }

    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();

    Double_t xSection7TeV                       = 62.22*1e-3;
    Double_t recalcBarn                         = 1e12; //NLO in pbarn!!!!

    Width_t widthLinesBoxes                     = 1.4;
    Width_t widthCommonFit                      = 2;

    // Definition of colors, styles and markers sizes
    Color_t colorComb                           = kBlue+2;
    Style_t markerStyleComb                     = 20;
    Size_t  markerSizeComb                      = 2;

    Color_t colorCombLowPt                      = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kFALSE);
    Color_t colorCombHighPt                     = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kTRUE);
    Style_t markerStyleCombLowPt                = 20;
    Style_t markerStyleCombHighPt               = 20;
    Size_t  markerSizeComparison                = 2;

    TString nameMeasGlobal[11]                  = {"PCM",
                                                   "PHOS",
                                                   "EMCal",
                                                   "PCM-PHOS",
                                                   "PCM-EMCal",
                                                   "PCM-Dalitz",
                                                   "PHOS-Dalitz",
                                                   "EMCal-Dalitz",
                                                   "EMCal",
                                                   "EMCAL merged",
                                                   "PCMOtherDataset"};

    //location of alternative EMCal measurement
    Int_t iEviEMCal = 8;
    if(!useDanielmeas){
      nameMeasGlobal[2] = "EMCal high pT";
      nameMeasGlobal[8] = "EMCal";
      iEviEMCal = 2;
    }else{
      nameMeasGlobal[2] = "EMCal";
      nameMeasGlobal[8] = "EMCal high pT";
    }

    Color_t colorDet[11];
    Color_t colorDetMC[11];
    Style_t markerStyleDet[11];
    Style_t markerStyleDetMC[11];
    Size_t  markerSizeDet[11];
    Size_t  markerSizeDetMC[11];

    Style_t styleMarkerNLOMuHalf                = 24;
    Style_t styleMarkerNLOMuOne                 = 27;
    Style_t styleMarkerNLOMuTwo                 = 30;
    Style_t styleLineNLOMuHalf                  = 8;
    Style_t styleLineNLOMuOne                   = 7;
    Style_t styleLineNLOMuTwo                   = 4;
    Style_t styleLineNLOMuTwoBKK                = 3;
    Style_t styleLineNLOMuTwoDSS                = 6;
    Size_t  sizeMarkerNLO                       = 1;
    Width_t widthLineNLO                        = 2.;

    for (Int_t i = 0; i < 11; i++){
        colorDet[i]                             = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
        colorDetMC[i]                           = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
        markerStyleDet[i]                       = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
        markerStyleDetMC[i]                     = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
        markerSizeDet[i]                        = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
        markerSizeDetMC[i]                      = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
    }

    TFile* inputFile[11]                     = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TDirectory* directoryPi0[11]             = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TDirectory* directoryEta[11]             = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    inputFile[0]                            = new TFile(fileNamePCM.Data());
    inputFile[1]                            = new TFile(fileNamePHOS.Data());
    inputFile[2]                            = new TFile(fileNameEMCAL.Data());
    inputFile[3]                            = new TFile(fileNamePCMPHOS.Data());
    inputFile[4]                            = new TFile(fileNamePCMEMCAL.Data());
    inputFile[8]                            = new TFile(fileNameEMCAL2.Data());

    for (Int_t i = 0; i < 11; i++){
      if(inputFile[i]){
        if(inputFile[i]->IsZombie()){
          cout << endl;
          cout << "*********************************************" << endl;
          cout << "inputFile #'" << i << "' could not be opened!" << endl;
          cout << "*********************************************" << endl;
          cout << endl;
          directoryPi0[i] = NULL;
          directoryEta[i] = NULL;
        }else{
          directoryPi0[i]                     = (TDirectory*)inputFile[i]->Get("Pi07TeV");
          directoryEta[i]                     = (TDirectory*)inputFile[i]->Get("Eta7TeV");
        }
      }else{
        directoryPi0[i] = NULL;
        directoryEta[i] = NULL;
      }
    }
        
    TH1D* histoPi0Mass[11]                                  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0FWHMMeV[11]                               = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0TrueMass[11]                              = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0TrueFWHMMeV[11]                           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0Acc[11]                                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0TrueEffPt[11]                             = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0AccTimesEff[11]                           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    TH1D* histoPi0InvCrossSection[11]                       = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphPi0InvCrossSectionSys[11]       = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphPi0InvCrossSectionStat[11]      = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    TH1D* histoEtaMass[11]                                  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaFWHMMeV[11]                               = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaTrueMass[11]                              = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaTrueFWHMMeV[11]                           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaAcc[11]                                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaTrueEffPt[11]                             = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaAccTimesEff[11]                           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    TH1D* histoEtaInvCrossSection[11]                       = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaInvCrossSectionSys[11]       = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaInvCrossSectionStat[11]      = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    TH1D* histoEtaToPi0Stat[11]                             = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaToPi0Stat[11]                = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaToPi0Sys[11]                 = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    TH1D* histoPi0InvMassSigPlusBG[11]                      = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvMassSig[11]                            = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvMassSigRemBGSub[11]                    = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvMassBG[11]                             = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvMassRemBG[11]                          = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvMassBGTot[11]                          = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TF1* fitPi0InvMassSig[11]                               = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TF1* fitPi0InvMassBG[11]                                = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    Bool_t haveAllPi0InvMass[11]                            = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
    TString strInvMassBin[11]                               = {"04", "22", "3To3_2", "04",""};
//     TString strInvMassBin[11]                              = {"04", "22", "12", "04",""}; //DanielEMCAL
    Double_t rapidityMeas[11]                               = {1.6, 1., 1.6, 1.6, 1.6, 1., 1., 1., 1., 1.};

    if(!useDanielmeas){
      rapidityMeas[2] = 1.;
      rapidityMeas[8] = 1.6;
      nameMeasGlobal[2] = "EMCal2011";
      nameMeasGlobal[8] = "EMCal2010";
    }else{
      nameMeasGlobal[2] = "EMCal2010";
      nameMeasGlobal[8] = "EMCal2011";
    }

    if(numbersofmeas < 5) nameMeasGlobal[2] = "EMCal";

    Double_t minPtPi0 = 0.23;
    Double_t maxPtPi0 = 32.0;
    Double_t minPtEta = 0.33;
    Double_t maxPtEta = 18.0;

    Bool_t doOutput = kTRUE;

    // *******************************************************************************************************
    // ************************** read input files ***********************************************************
    // *******************************************************************************************************

    for (Int_t i = 0; i < 11; i++){
      cout << "reading from " << nameMeasGlobal[i].Data() << " measurement..." << endl;
      if(directoryPi0[i]){
        // load mass/width/effi plots
        histoPi0Mass[i]                     = (TH1D*)directoryPi0[i]->Get("Pi0_Mass_data_INT1");
          if(histoPi0Mass[i]) histoPi0Mass[i]->Scale(1000.);
        histoPi0FWHMMeV[i]                  = (TH1D*)directoryPi0[i]->Get("Pi0_Width_data_INT1");
          if(histoPi0FWHMMeV[i]) histoPi0FWHMMeV[i]->Scale(1000.);
        histoPi0TrueMass[i]                 = (TH1D*)directoryPi0[i]->Get("Pi0_Mass_MC_INT1");
          if(histoPi0TrueMass[i]) histoPi0TrueMass[i]->Scale(1000.);
        histoPi0TrueFWHMMeV[i]              = (TH1D*)directoryPi0[i]->Get("Pi0_Width_MC_INT1");
          if(histoPi0TrueFWHMMeV[i]) histoPi0TrueFWHMMeV[i]->Scale(1000.);
        histoPi0Acc[i]                      = (TH1D*)directoryPi0[i]->Get("AcceptancePi0_INT1");
        histoPi0TrueEffPt[i]                = (TH1D*)directoryPi0[i]->Get("EfficiencyPi0_INT1");

        if(i==1 || i==3){ // PHOS & PCM-PHOS
          histoPi0Mass[i]                     = (TH1D*)directoryPi0[i]->Get("MassPi0");
          histoPi0FWHMMeV[i]                  = (TH1D*)directoryPi0[i]->Get("FWHMPi0MeV");
          histoPi0TrueMass[i]                 = (TH1D*)directoryPi0[i]->Get("TrueMassPi0");
          histoPi0TrueFWHMMeV[i]              = (TH1D*)directoryPi0[i]->Get("TrueFWHMPi0MeV");
          histoPi0Acc[i]                      = (TH1D*)directoryPi0[i]->Get("AcceptancePi0");
          histoPi0TrueEffPt[i]                = (TH1D*)directoryPi0[i]->Get("EfficiencyPi0");
        }

        if(i!=iEviEMCal){
          histoPi0AccTimesEff[i]            = (TH1D*)histoPi0TrueEffPt[i]->Clone(Form("histoPi0AccTimesEff%s",nameMeasGlobal[i].Data()));
          histoPi0AccTimesEff[i]            ->Multiply(histoPi0Acc[i]);
          histoPi0AccTimesEff[i]            ->Scale(2*TMath::Pi()*rapidityMeas[i]);
        }else{
          histoPi0AccTimesEff[i]            = (TH1D*)directoryPi0[i]->Get("EfficiencyPi0");
        }

        // load cross section systematics and datapoints
        histoPi0InvCrossSection[i]          = (TH1D*)directoryPi0[i]->Get("InvCrossSectionPi0");
        graphPi0InvCrossSectionStat[i]      = new TGraphAsymmErrors(histoPi0InvCrossSection[i]);
          while (graphPi0InvCrossSectionStat[i]->GetY()[0] <= 1E-50 ) graphPi0InvCrossSectionStat[i]->RemovePoint(0);
          while (graphPi0InvCrossSectionStat[i]->GetY()[graphPi0InvCrossSectionStat[i]->GetN()-1] <= 1E-50 ) graphPi0InvCrossSectionStat[i]->RemovePoint(graphPi0InvCrossSectionStat[i]->GetN()-1);
        graphPi0InvCrossSectionSys[i]       = (TGraphAsymmErrors*)directoryPi0[i]->Get("InvCrossSectionPi0Sys");
          while (graphPi0InvCrossSectionSys[i]->GetY()[0] <= 1E-50 ) graphPi0InvCrossSectionSys[i]->RemovePoint(0);
          while (graphPi0InvCrossSectionSys[i]->GetY()[graphPi0InvCrossSectionSys[i]->GetN()-1] <= 1E-50 ) graphPi0InvCrossSectionSys[i]->RemovePoint(graphPi0InvCrossSectionSys[i]->GetN()-1);

        cout << nameMeasGlobal[i].Data() << " pi0 stat:" << graphPi0InvCrossSectionStat[i] << endl;
        if(doOutput) graphPi0InvCrossSectionStat[i]->Print();
        cout << nameMeasGlobal[i].Data() << " pi0 sys:" << graphPi0InvCrossSectionSys[i] << endl;
        if(doOutput) graphPi0InvCrossSectionSys[i]->Print();

        // load invariant mass example bins
        histoPi0InvMassSig[i]               = (TH1D*)directoryPi0[i]->Get(Form("InvMassSig_PtBin%s",strInvMassBin[i].Data()));
        histoPi0InvMassSigPlusBG[i]         = (TH1D*)directoryPi0[i]->Get(Form("InvMassSigPlusBG_PtBin%s",strInvMassBin[i].Data()));
        histoPi0InvMassBG[i]                = (TH1D*)directoryPi0[i]->Get(Form("InvMassBG_PtBin%s",strInvMassBin[i].Data()));
        fitPi0InvMassSig[i]                 = (TF1*)directoryPi0[i]->Get(Form("FitInvMassSig_PtBin%s",strInvMassBin[i].Data()));

        if (histoPi0InvMassSig[i] && histoPi0InvMassSigPlusBG[i] && histoPi0InvMassBG[i] && fitPi0InvMassSig[i]){
          haveAllPi0InvMass[i]              = kTRUE;
        }
      }
      if(directoryEta[i]){
        // load mass/width/effi plots
        histoEtaMass[i]                     = (TH1D*)directoryEta[i]->Get("Eta_Mass_data_INT1");
          if(histoEtaMass[i]) histoEtaMass[i]->Scale(1000.);
        histoEtaFWHMMeV[i]                  = (TH1D*)directoryEta[i]->Get("Eta_Width_data_INT1");
          if(histoEtaFWHMMeV[i]) histoEtaFWHMMeV[i]->Scale(1000.);
        histoEtaTrueMass[i]                 = (TH1D*)directoryEta[i]->Get("Eta_Mass_MC_INT1");
          if(histoEtaTrueMass[i]) histoEtaTrueMass[i]->Scale(1000.);
        histoEtaTrueFWHMMeV[i]              = (TH1D*)directoryEta[i]->Get("Eta_Width_MC_INT1");
          if(histoEtaTrueFWHMMeV[i]) histoEtaTrueFWHMMeV[i]->Scale(1000.);
        histoEtaAcc[i]                      = (TH1D*)directoryEta[i]->Get("AcceptanceEta_INT1");
        histoEtaTrueEffPt[i]                = (TH1D*)directoryEta[i]->Get("EfficiencyEta_INT1");

        if(i==1 || i==3){ // PHOS & PCM-PHOS
          histoEtaMass[i]                     = (TH1D*)directoryEta[i]->Get("MassEta");
          histoEtaFWHMMeV[i]                  = (TH1D*)directoryEta[i]->Get("FWHMEtaMeV");
          histoEtaTrueMass[i]                 = (TH1D*)directoryEta[i]->Get("TrueMassEta");
          histoEtaTrueFWHMMeV[i]              = (TH1D*)directoryEta[i]->Get("TrueFWHMEtaMeV");
          histoEtaAcc[i]                      = (TH1D*)directoryEta[i]->Get("AcceptanceEta");
          histoEtaTrueEffPt[i]                = (TH1D*)directoryEta[i]->Get("EfficiencyEta");
        }

        if(i!=iEviEMCal){
          histoEtaAccTimesEff[i]            = (TH1D*)histoEtaTrueEffPt[i]->Clone(Form("histoEtaAccTimesEff%s",nameMeasGlobal[i].Data()));
          histoEtaAccTimesEff[i]            ->Multiply(histoEtaAcc[i]);
          histoEtaAccTimesEff[i]            ->Scale(2*TMath::Pi()*rapidityMeas[i]);
        }else{
          histoEtaAccTimesEff[i]            = (TH1D*)directoryEta[i]->Get("EfficiencyEta");
        }

        // load cross section systematics and datapoints
        histoEtaInvCrossSection[i]          = (TH1D*)directoryEta[i]->Get("InvCrossSectionEta");
        graphEtaInvCrossSectionStat[i]      = new TGraphAsymmErrors(histoEtaInvCrossSection[i]);
        graphEtaInvCrossSectionSys[i]       = (TGraphAsymmErrors*)directoryEta[i]->Get("InvCrossSectionEtaSys");
          for(Int_t iB=0;iB<histoEtaInvCrossSection[i]->GetNbinsX();iB++){if(histoEtaInvCrossSection[i]->GetBinContent(iB) <= 1E-50) histoEtaInvCrossSection[i]->SetBinContent(iB,0.);}
          while (graphEtaInvCrossSectionStat[i]->GetY()[0] <= 1E-50 ) graphEtaInvCrossSectionStat[i]->RemovePoint(0);
          while (graphEtaInvCrossSectionStat[i]->GetY()[graphEtaInvCrossSectionStat[i]->GetN()-1] <= 1E-50 ) graphEtaInvCrossSectionStat[i]->RemovePoint(graphEtaInvCrossSectionStat[i]->GetN()-1);
          while (graphEtaInvCrossSectionSys[i]->GetY()[0] <= 1E-50 ) graphEtaInvCrossSectionSys[i]->RemovePoint(0);
          while (graphEtaInvCrossSectionSys[i]->GetY()[graphEtaInvCrossSectionSys[i]->GetN()-1] <= 1E-50 ) graphEtaInvCrossSectionSys[i]->RemovePoint(graphEtaInvCrossSectionSys[i]->GetN()-1);

        cout << nameMeasGlobal[i].Data() << " eta stat:" << graphEtaInvCrossSectionStat[i] << endl;
        if(doOutput) graphEtaInvCrossSectionStat[i]->Print();
        cout << nameMeasGlobal[i].Data() << " eta sys:" << graphEtaInvCrossSectionSys[i] << endl;
        if(doOutput) graphEtaInvCrossSectionSys[i]->Print();


        histoEtaToPi0Stat[i]               = (TH1D*)directoryEta[i]->Get("EtaToPi0YShiftedStatError");
        graphEtaToPi0Stat[i]               = (TGraphAsymmErrors*)directoryEta[i]->Get("graphEtaToPi0YShiftedStatError");
        graphEtaToPi0Sys[i]                = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0YShiftedSystError");
        if(i==1){
            histoEtaToPi0Stat[i]           = (TH1D*)directoryEta[i]->Get("RatioEtaPi0");
            graphEtaToPi0Stat[i]           = new TGraphAsymmErrors(histoEtaToPi0Stat[i]);
            graphEtaToPi0Sys[i]            = (TGraphAsymmErrors*)directoryEta[i]->Get("RatioEtaPi0Sys");
        }
        if(i==3){
            histoEtaToPi0Stat[i]           = (TH1D*)directoryEta[i]->Get("EtatoPi0RatioConversion");
            graphEtaToPi0Stat[i]           = new TGraphAsymmErrors(histoEtaToPi0Stat[i]);
            graphEtaToPi0Sys[i]            = (TGraphAsymmErrors*)directoryEta[i]->Get("EtatoPi0RatioConversionSys");
        }
        if(i==iEviEMCal){
            histoEtaToPi0Stat[i]           = (TH1D*)directoryEta[i]->Get("EtaToPi0RatioEMCal");
            graphEtaToPi0Stat[i]           = new TGraphAsymmErrors(histoEtaToPi0Stat[i]);
            graphEtaToPi0Sys[i]            = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0RatioEMCalSys");
        }

        for(Int_t iB=0;iB<histoEtaToPi0Stat[i]->GetNbinsX();iB++){if(histoEtaToPi0Stat[i]->GetBinContent(iB) <= 1E-50) histoEtaToPi0Stat[i]->SetBinContent(iB,0.);}
        while (graphEtaToPi0Stat[i]->GetY()[0] <= 1E-50 ) graphEtaToPi0Stat[i]->RemovePoint(0);
        while (graphEtaToPi0Stat[i]->GetY()[graphEtaToPi0Stat[i]->GetN()-1] <= 1E-50 ) graphEtaToPi0Stat[i]->RemovePoint(graphEtaToPi0Stat[i]->GetN()-1);
        while (graphEtaToPi0Sys[i]->GetY()[0] <= 1E-50 ) graphEtaToPi0Sys[i]->RemovePoint(0);
        while (graphEtaToPi0Sys[i]->GetY()[graphEtaToPi0Sys[i]->GetN()-1] <= 1E-50 ) graphEtaToPi0Sys[i]->RemovePoint(graphEtaToPi0Sys[i]->GetN()-1);

        cout << nameMeasGlobal[i].Data() << " eta/pi0 stat:" << histoEtaToPi0Stat[i] << endl;
        cout << nameMeasGlobal[i].Data() << " eta/pi0 stat:" << graphEtaToPi0Stat[i] << endl;
        if(doOutput) graphEtaToPi0Stat[i]->Print();
        cout << nameMeasGlobal[i].Data() << " eta/pi0 sys:" << graphEtaToPi0Sys[i] << endl;
        if(doOutput) graphEtaToPi0Sys[i]->Print();
      }


      for (Int_t i = 0; i < 11; i++){
          if (haveAllPi0InvMass[i]){
              histoPi0InvMassBGTot[i]      = (TH1D*)histoPi0InvMassBG[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameMeasGlobal[i].Data()));
              histoPi0InvMassSigRemBGSub[i]= (TH1D*)histoPi0InvMassSig[i]->Clone(Form("Pi0_InvMassSigRemBGSub_Example_%s",nameMeasGlobal[i].Data()));
          }
      }
    }

    // calculate necessary example bin histograms for PCM
//     for (Int_t i = 0; i < numbersofmeas; i++){
//         if (i!=1&&i!=2&&i!=4){
//             cout << "found all example invariant mass components for " << nameMeasGlobal[i].Data() << endl;
//             histoPi0InvMassSig[i]->Fit(fitPi0InvMassSig[i],"QRME0");
//             for (Int_t l=0; l < 6; l++){
//                 cout << fitPi0InvMassSig[i]->GetParameter(l) << "\t +- " << fitPi0InvMassSig[i]->GetParError(l) << endl;
//             }
//             fitPi0InvMassBG[i]                  = new TF1("Linearpp","[0]+[1]*x",0.02,0.25);
//             fitPi0InvMassBG[i]->SetParameter(0, fitPi0InvMassSig[i]->GetParameter(4));
//             fitPi0InvMassBG[i]->SetParameter(1, fitPi0InvMassSig[i]->GetParameter(5));
//             TVirtualFitter * fitter             = TVirtualFitter::GetFitter();
//             Int_t nFreePar                      = fitPi0InvMassSig[i]->GetNumberFreeParameters();
//             double * covMatrix                  = fitter->GetCovarianceMatrix();
//
//             histoPi0InvMassRemBG[i]             = (TH1D*)histoPi0InvMassBG[i]->Clone(Form("Pi0_InvMassRemBG_Example_%s",nameMeasGlobal[i].Data()));
//             for (Int_t j = 1; j < histoPi0InvMassRemBG[i]->GetNbinsX()+1; j++){
//                 histoPi0InvMassRemBG[i]->SetBinContent(j,0);
//                 histoPi0InvMassRemBG[i]->SetBinError(j,0);
//             }
//             for (Int_t j = histoPi0InvMassSig[i]->GetXaxis()->FindBin(0.01); j < histoPi0InvMassSig[i]->GetXaxis()->FindBin(0.30)+1; j++){
//                 Double_t startBinEdge           = histoPi0InvMassSig[i]->GetXaxis()->GetBinLowEdge(j);
//                 Double_t endBinEdge             = histoPi0InvMassSig[i]->GetXaxis()->GetBinUpEdge(j);
//                 Double_t intLinearBack          = fitPi0InvMassBG[i]->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
//                 Double_t errorLinearBck         = pow(( pow( (endBinEdge-startBinEdge)*fitPi0InvMassSig[i]->GetParError(4),2) +
//                                                   pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPi0InvMassSig[i]->GetParError(5),2)
//                                                   +2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*
//                                                   (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
//                 histoPi0InvMassRemBG[i]->SetBinContent(j,intLinearBack);
//                 histoPi0InvMassRemBG[i]->SetBinError(j,errorLinearBck);
//             }
//             histoPi0InvMassBGTot[i]             = (TH1D*)histoPi0InvMassBG[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameMeasGlobal[i].Data()));
//             histoPi0InvMassBGTot[i]->Sumw2();
//             histoPi0InvMassBGTot[i]->Add(histoPi0InvMassRemBG[i]);
//             histoPi0InvMassSigRemBGSub[i]       = (TH1D*)histoPi0InvMassSig[i]->Clone(Form("Pi0_InvMassSigRemBGSub_Example_%s",nameMeasGlobal[i].Data()));
//             histoPi0InvMassSigRemBGSub[i]->Sumw2();
//             histoPi0InvMassSigRemBGSub[i]->Add(histoPi0InvMassRemBG[i],-1);
//             fitPi0InvMassSig[i]->SetParameter(4, 0);
//             fitPi0InvMassSig[i]->SetParameter(5, 0);
//         }
//     }


    Int_t nBinsPi0 = 48;
    Double_t xPtLimits[49]                      =  { 0.0, 0.3, 0.4, 0.5, 0.6,
                                                     0.7, 0.8, 0.9, 1.0, 1.1,
                                                     1.2, 1.3, 1.4, 1.5, 1.6,
                                                     1.7, 1.8, 1.9, 2.0, 2.1,
                                                     2.2, 2.3, 2.4, 2.6, 2.8,
                                                     3.0, 3.2, 3.4, 3.6, 3.8,
                                                     4.0, 4.3, 4.6, 5.0, 5.5,
                                                     6.0, 6.5, 7.0, 8.0, 9.0,
                                                     10.0, 11.0, 12.0, 13.0, 14.0,
                                                     16.0, 18.0, 20.0, 25.0
                                                    };

    Int_t nBinsEta = 17;
    Double_t xPtLimitsEta[49]                   =  { 0.0, 0.4, 0.6, 0.8, 1.0,
                                                     1.4, 1.8, 2.2, 2.6, 3.0,
                                                     3.5, 4.0, 5.0, 6.0, 8.0,
                                                     10.0, 12.0, 14.0, 15.0, 16.0 ,
                                                     18.0, 20.0, 25.0, 35.0
                                                   };
    Double_t xPtLimitsEtaToPi0[49]              =  { 0.0, 0.4, 0.6, 0.8, 1.0,
                                                     1.4, 1.8, 2.2, 2.6, 3.0,
                                                     3.5, 4.0, 5.0, 6.0, 8.0,
                                                     10.0, 12.0, 14.0, 15.0, 16.0,
                                                     18.0, 20.0, 25.0, 35.0
                                                   };
    if(!useDanielmeas){
      // with Evi
      nBinsEta = 24;
      Double_t xPtLimitsEta_Evi[49]              =  { 0.0, 0.4, 0.6, 0.8, 1.0,
                                                      1.4, 1.8, 2.2, 2.6, 3.0,
                                                      3.5, 4.0, 5.0, 6.0, 8.0,
                                                      10.0, 12.0, 13.0, 14.0, 15.0,
                                                      16.0 ,18.0, 20.0, 25.0, 35.0
                                                    };
      for(Int_t i=0; i<=nBinsEta; i++) xPtLimitsEta[i] = xPtLimitsEta_Evi[i];

      Double_t xPtLimitsEtaToPi0_Evi[49]         =  { 0.0, 0.4, 0.6, 0.8, 1.0,
                                                      1.4, 1.8, 2.2, 2.6, 3.0,
                                                      3.5, 4.0, 5.0, 6.0, 8.0,
                                                      10.0, 12.0, 13.0, 14.0, 15.0,
                                                      16.0, 18.0, 20.0, 25.0
                                                    };
      for(Int_t i=0; i<nBinsEta; i++) xPtLimitsEtaToPi0[i] = xPtLimitsEtaToPi0_Evi[i];
    }

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

    TH1D* statErrorCollection[11];
    TH1D* statErrorCollectionEta[11];
    TH1D* statErrorCollectionEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorCollection[i]                  = NULL;
        statErrorCollectionEta[i]               = NULL;
        statErrorCollectionEtaToPi0[i]          = NULL;
    }
    for (Int_t i = 0; i< 11; i++){
        if(directoryPi0[i]) statErrorCollection[i]                  = (TH1D*)histoPi0InvCrossSection[i]->Clone(Form("statErr%sPi0",nameMeasGlobal[i].Data()));
        if(directoryEta[i]) statErrorCollectionEta[i]               = (TH1D*)histoEtaInvCrossSection[i]->Clone(Form("statErr%sEta",nameMeasGlobal[i].Data()));
        if(directoryEta[i]) statErrorCollectionEtaToPi0[i]          = (TH1D*)histoEtaToPi0Stat[i]->Clone(Form("statErr%sEtaToPi0",nameMeasGlobal[i].Data()));
    }
    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollection[11];
    TGraphAsymmErrors* sysErrorCollectionEta[11];
    TGraphAsymmErrors* sysErrorCollectionEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollection[i]                   = NULL;
        sysErrorCollectionEta[i]                = NULL;
        sysErrorCollectionEtaToPi0[i]           = NULL;
    }
    for (Int_t i = 0; i< 11; i++){
        if(directoryPi0[i]) sysErrorCollection[i]                   = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i]->Clone(Form("sysErr%sPi0",nameMeasGlobal[i].Data()));
        if(directoryEta[i]) sysErrorCollectionEta[i]                = (TGraphAsymmErrors*)graphEtaInvCrossSectionSys[i]->Clone(Form("sysErr%sEta",nameMeasGlobal[i].Data()));
        if(directoryEta[i]) sysErrorCollectionEtaToPi0[i]           = (TGraphAsymmErrors*)graphEtaToPi0Sys[i]->Clone(Form("sysErr%sEtaToPi0",nameMeasGlobal[i].Data()));
    }



    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match
    //                                            PCM,PHOS,EMC,PCMPHOS,PCMEMC,        EMC
    Int_t offSets[11]                           =  {0,    8,  1,     2,     0, 0,0,0,   6,0,0};
    Int_t offSetsSys[11]                        =  {1,    8, 11,     3,     6, 0,0,0,   6,0,0};
    if(!useDanielmeas){
      offSets[2]    = 6;
      offSetsSys[2] = 6;
      offSets[8]    = 1;
      offSetsSys[8] = 11;
    }

    //                                            PCM,PHOS,EMC,PCMPHOS,PCMEMC,         EMC
    Int_t offSetsEta[11]                        =  {0,    4,  1,     2,      1, 0,0,0,   4,0,0};
    Int_t offSetsSysEta[11]                     =  {1,    4,  7,     3,      5, 0,0,0,   9,0,0};
    if(!useDanielmeas){
      offSetsEta[2]    = 4;
      offSetsSysEta[2] = 9;
      offSetsEta[8]    = 1;
      offSetsSysEta[8] = 7;
    }

    //                                            PCM,PHOS,EMC,PCMPHOS,PCMEMC,         EMC
    Int_t offSetsEtaToPi0[11]                   =  {0,    4,  1,     2,      1, 0,0,0,   4,0,0};
    Int_t offSetsSysEtaToPi0[11]                =  {1,    4,  7,     3,      5, 0,0,0,   9,0,0};
    if(!useDanielmeas){
      offSetsEtaToPi0[2]    = 4;
      offSetsSysEtaToPi0[2] = 9;
      offSetsEtaToPi0[8]    = 1;
      offSetsSysEtaToPi0[8] = 7;
    }
    
    //**********************************************************************************************************************
    //**********************************************************************************************************************
    //**********************************************************************************************************************

    TH1D* statErrorRelCollection[11];
    TH1D* statErrorRelCollectionEta[11];
    TH1D* statErrorRelCollectionEtaToPi0[11];

    for (Int_t i = 0; i < 11; i++){
        if (statErrorCollection[i]) statErrorRelCollection[i] = CalculateRelErrUpTH1D( statErrorCollection[i], Form("relativeStatError_%s", nameMeasGlobal[i].Data()));
        else statErrorRelCollection[i] = NULL;
        if (statErrorCollectionEta[i]) statErrorRelCollectionEta[i] = CalculateRelErrUpTH1D( statErrorCollectionEta[i], Form("relativeStatError_%s", nameMeasGlobal[i].Data()));
        else statErrorRelCollectionEta[i] = NULL;
        if (statErrorCollectionEtaToPi0[i]) statErrorRelCollectionEtaToPi0[i] = CalculateRelErrUpTH1D( statErrorCollectionEtaToPi0[i], Form("relativeStatError_%s", nameMeasGlobal[i].Data()));
        else statErrorRelCollectionEtaToPi0[i] = NULL;
    }

    TGraphAsymmErrors* sysErrorRelCollection[11];
    TGraphAsymmErrors* sysErrorRelCollectionEta[11];
    TGraphAsymmErrors* sysErrorRelCollectionEtaToPi0[11];

    for (Int_t i = 0; i < 11; i++){
        if (sysErrorCollection[i]) sysErrorRelCollection[i] = CalculateRelErrUpAsymmGraph( sysErrorCollection[i], Form("relativeSysError_%s", nameMeasGlobal[i].Data()));
        else sysErrorRelCollection[i] = NULL;
        if (sysErrorCollectionEta[i]) sysErrorRelCollectionEta[i] = CalculateRelErrUpAsymmGraph( sysErrorCollectionEta[i], Form("relativeSysError_%s", nameMeasGlobal[i].Data()));
        else sysErrorRelCollectionEta[i] = NULL;
        if (sysErrorCollectionEtaToPi0[i]) sysErrorRelCollectionEtaToPi0[i] = CalculateRelErrUpAsymmGraph( sysErrorCollectionEtaToPi0[i], Form("relativeSysError_%s", nameMeasGlobal[i].Data()));
        else sysErrorRelCollectionEtaToPi0[i] = NULL;
    }

    TGraph* graphWeights[11];
    TGraph* graphWeightsEta[11];
    TGraph* graphWeightsEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeights[i] = NULL;
        graphWeightsEta[i] = NULL;
        graphWeightsEtaToPi0[i] = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNameOutputWeighting                       = Form("%s/Weighting.dat",outputDir.Data());
    TString fileNameOutputWeightingEta                    = Form("%s/WeightingEta.dat",outputDir.Data());
    TString fileNameOutputWeightingEtaToPi0               = Form("%s/WeightingEtaToPi0.dat",outputDir.Data());

    TGraphAsymmErrors* graphCombPi0InvCrossSectionStatPCMEMCPHOS= NULL;
    TGraphAsymmErrors* graphCombPi0InvCrossSectionSysPCMEMCPHOS = NULL;
    TGraphAsymmErrors* graphCombPi0InvCrossSectionTot = CombinePtPointsSpectraFullCorrMat( statErrorCollection, sysErrorCollection,
                                                                                           xPtLimits, nBinsPi0,
                                                                                           offSets, offSetsSys,
                                                                                           graphCombPi0InvCrossSectionStatPCMEMCPHOS, graphCombPi0InvCrossSectionSysPCMEMCPHOS,
                                                                                           fileNameOutputWeighting,"7TeV", "Pi0", kTRUE,
                                                                                           0x0, fileInputCorrFactors
                                                                                          );

    //return;
    if(doOutput) graphCombPi0InvCrossSectionStatPCMEMCPHOS->Print();
    
    TGraphAsymmErrors* graphCombEtaInvCrossSectionStatPCMEMCPHOS= NULL;
    TGraphAsymmErrors* graphCombEtaInvCrossSectionSysPCMEMCPHOS = NULL;
    TGraphAsymmErrors* graphCombEtaInvCrossSectionTot = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEta, sysErrorCollectionEta,
                                                                                           xPtLimitsEta, nBinsEta,
                                                                                           offSetsEta, offSetsSysEta,
                                                                                           graphCombEtaInvCrossSectionStatPCMEMCPHOS, graphCombEtaInvCrossSectionSysPCMEMCPHOS,
                                                                                           fileNameOutputWeightingEta,"7TeV", "Eta", kTRUE,
                                                                                           0x0, fileInputCorrFactors
                                                                                         );
    //return;
    if(doOutput) graphCombEtaInvCrossSectionStatPCMEMCPHOS->Print();
    
    
    TGraphAsymmErrors* graphCombEtaToPi0StatPCMEMCPHOS= NULL;
    TGraphAsymmErrors* graphCombEtaToPi0SysPCMEMCPHOS = NULL;
    TGraphAsymmErrors* graphCombEtaToPi0Tot = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEtaToPi0, sysErrorCollectionEtaToPi0,
                                                                                 xPtLimitsEtaToPi0, nBinsEta-1,
                                                                                 offSetsEtaToPi0, offSetsSysEtaToPi0,
                                                                                 graphCombEtaToPi0StatPCMEMCPHOS, graphCombEtaToPi0SysPCMEMCPHOS,
                                                                                 fileNameOutputWeightingEtaToPi0,"7TeV", "EtaToPi0", kTRUE,
                                                                                 0x0, fileInputCorrFactors
                                                                               );
    //return;
    if(doOutput) graphCombEtaToPi0StatPCMEMCPHOS->Print();
    
    
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
    Double_t xValuesRead[50];
    Double_t weightsRead[11][50];
    Int_t availableMeas[11]                     = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    Int_t nMeasSet                              = numbersofmeas;
    Int_t nPtBinsRead                           = 0;
    while(!fileWeightsRead.eof() && nPtBinsRead < 50){
        TString garbage                         = "";
        if (nPtBinsRead == 0){
            fileWeightsRead >> garbage ;
            for (Int_t i = 0; i < nMeasSet; i++){
                fileWeightsRead >> availableMeas[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableMeas[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsRead >> xValuesRead[nPtBinsRead-1];
            for (Int_t i = 0; i < nMeasSet; i++){
                fileWeightsRead >> weightsRead[availableMeas[i]][nPtBinsRead-1] ;
            }
            cout << "read: "<<  nPtBinsRead << "\t"<< xValuesRead[nPtBinsRead-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSet; i++){
                cout << weightsRead[availableMeas[i]][nPtBinsRead-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsRead++;
    }
    nPtBinsRead = nPtBinsRead-2 ;
    fileWeightsRead.close();

    for (Int_t i = 0; i < nMeasSet; i++){
        graphWeights[availableMeas[i]]  = new TGraph(nPtBinsRead,xValuesRead,weightsRead[availableMeas[i]]);
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

        TLegend* legendAccWeights               = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSet*1.35), 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            DrawGammaSetMarkerTGraph(graphWeights[availableMeas[i]],
                                    markerStyleDet[availableMeas[i]],
                                    markerSizeDet[availableMeas[i]]*0.5,
                                    colorDet[availableMeas[i]] ,
                                    colorDet[availableMeas[i]]);
            graphWeights[availableMeas[i]]->Draw("p,same,e1");
            legendAccWeights->AddEntry(graphWeights[availableMeas[i]],nameMeasGlobal[availableMeas[i]],"p");
        }
        legendAccWeights->Draw();
        TLatex *labelWeightsEnergy              = new TLatex(0.7,0.20,collisionSystem7TeV.Data());
        SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
        labelWeightsEnergy->SetTextFont(43);
        labelWeightsEnergy->Draw();
        TLatex *labelWeightsPi0                 = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
        labelWeightsPi0->SetTextFont(43);
        labelWeightsPi0->Draw();

        DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/Pi0_WeightsMethods.%s",outputDir.Data(),suffix.Data()));


    //*********************************************************************************************************************
    //************************************ Visualize relative errors ******************************************************
    //*********************************************************************************************************************

    TCanvas* canvasRelSysErr                    = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();

    TH2F * histo2DRelSysErr = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,minPtPi0,maxPtPi0,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetNoExponent();
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,45.5);
    histo2DRelSysErr->Draw("copy");

        TLegend* legendRelSysErr                = GetAndSetLegend2(0.62, 0.94-(0.035*nMeasSet*1.35), 0.95, 0.94, 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollection[availableMeas[i]], markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]],
                                    colorDet[availableMeas[i]]);
            sysErrorRelCollection[availableMeas[i]]->Draw("p,same,e1");
            legendRelSysErr->AddEntry(sysErrorRelCollection[availableMeas[i]],nameMeasGlobal[availableMeas[i]],"p");
        }
        legendRelSysErr->Draw();

        TLatex *labelRelSysErrEnergy            = new TLatex(0.15,0.89,collisionSystem7TeV.Data());
        SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrEnergy->SetTextFont(43);
        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrPi0               = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
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

        TLegend* legendRelStatErr               = GetAndSetLegend2(0.14, 0.94-(0.035*nMeasSet*1.35), 0.45, 0.94, 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            if (availableMeas[i]== 2){
                DrawGammaSetMarker(statErrorRelCollection[availableMeas[i]], markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]] ,
                            colorDet[availableMeas[i]]);
                TGraphAsymmErrors* graphDummy   = new TGraphAsymmErrors(statErrorRelCollection[availableMeas[i]]);
                DrawGammaSetMarkerTGraphAsym(graphDummy, markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]],
                                    colorDet[availableMeas[i]]);
                graphDummy->Draw("same,p,x0");
                legendRelStatErr->AddEntry(graphDummy,nameMeasGlobal[availableMeas[i]],"p");

            } else {
                DrawGammaSetMarker(statErrorRelCollection[availableMeas[i]], markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]] ,
                            colorDet[availableMeas[i]]);
                statErrorRelCollection[availableMeas[i]]->Draw("p,same,e1");
                legendRelStatErr->AddEntry(statErrorRelCollection[availableMeas[i]],nameMeasGlobal[availableMeas[i]],"p");

            }
        }
        legendRelStatErr->Draw();

        TLatex *labelRelStatErrEnergy           = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
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

    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelStat7TeV     = CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionStatPCMEMCPHOS, "graphCombPi0InvCrossSectionRelStat7TeV");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelSys7TeV      = CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionSysPCMEMCPHOS, "graphCombPi0InvCrossSectionRelSys7TeV");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelTot7TeV      = CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionTot, "graphCombPi0InvCrossSectionRelTot7TeV");

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

    DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelTot7TeV, markerStyleComb, markerSizeComb, colorComb , colorComb);
    graphCombPi0InvCrossSectionRelTot7TeV->Draw("p,same,e1");
    DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelStat7TeV, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
    graphCombPi0InvCrossSectionRelStat7TeV->Draw("l,x0,same,e1");
    DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelSys7TeV, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
    graphCombPi0InvCrossSectionRelSys7TeV->SetLineStyle(7);
    graphCombPi0InvCrossSectionRelSys7TeV->Draw("l,x0,same,e1");

    TLegend* legendRelTotErr2 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
    legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelTot7TeV,"tot","p");
    legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelStat7TeV,"stat","l");
    legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelSys7TeV,"sys","l");
    legendRelTotErr2->Draw();

    TLatex *labelRelTotErrEnergy            = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
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
    Double_t weightsReadEta[11][50];
    Int_t availableMeasEta[11]                  = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    Int_t nPtBinsReadEta                        = 0;
    while(!fileWeightsReadEta.eof() && nPtBinsReadEta < 50){
        TString garbage                         = "";
        if (nPtBinsReadEta == 0){
            fileWeightsReadEta >> garbage ;
            for (Int_t i = 0; i < nMeasSet; i++){
                fileWeightsReadEta >> availableMeasEta[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableMeasEta[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsReadEta >> xValuesReadEta[nPtBinsReadEta-1];
            for (Int_t i = 0; i < nMeasSet; i++){
                fileWeightsReadEta >> weightsReadEta[availableMeasEta[i]][nPtBinsReadEta-1] ;
            }
            cout << "read: "<<  nPtBinsReadEta << "\t"<< xValuesReadEta[nPtBinsReadEta-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSet; i++){
                cout << weightsReadEta[availableMeasEta[i]][nPtBinsReadEta-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsReadEta++;
    }
    nPtBinsReadEta = nPtBinsReadEta-2 ;
    fileWeightsReadEta.close();

    for (Int_t i = 0; i < nMeasSet; i++){
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

        TLegend* legendAccWeightsEta               = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSet*1.35), 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEta[availableMeasEta[i]],
                                    markerStyleDet[availableMeasEta[i]],
                                    markerSizeDet[availableMeasEta[i]]*0.5,
                                    colorDet[availableMeasEta[i]] ,
                                    colorDet[availableMeasEta[i]]);
            graphWeightsEta[availableMeasEta[i]]->Draw("p,same,e1");
            legendAccWeightsEta->AddEntry(graphWeightsEta[availableMeasEta[i]],nameMeasGlobal[availableMeasEta[i]],"p");
        }
        legendAccWeightsEta->Draw();
        TLatex *labelWeightsEnergyEta              = new TLatex(0.7,0.20,collisionSystem7TeV.Data());
        SetStyleTLatex( labelWeightsEnergyEta, 0.85*textSizeLabelsPixel,4);
        labelWeightsEnergyEta->SetTextFont(43);
        labelWeightsEnergyEta->Draw();
        TLatex *labelWeightsEta                 = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsEta, 0.85*textSizeLabelsPixel,4);
        labelWeightsEta->SetTextFont(43);
        labelWeightsEta->Draw();

        DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);

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
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,45.5);
    histo2DRelSysErr->Draw("copy");

        TLegend* legendRelSysErrEta                = GetAndSetLegend2(0.62, 0.94-(0.035*nMeasSet*1.35), 0.95, 0.94, 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionEta[availableMeasEta[i]], markerStyleDet[availableMeasEta[i]], markerSizeDet[availableMeasEta[i]]*0.5, colorDet[availableMeasEta[i]],
                                    colorDet[availableMeasEta[i]]);
            sysErrorRelCollectionEta[availableMeasEta[i]]->Draw("p,same,e1");
            legendRelSysErrEta->AddEntry(sysErrorRelCollectionEta[availableMeasEta[i]],nameMeasGlobal[availableMeasEta[i]],"p");
        }
        legendRelSysErrEta->Draw();

        TLatex *labelRelSysErrEnergyEta            = new TLatex(0.15,0.89,collisionSystem7TeV.Data());
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
    histo2DRelStatErrEta->GetYaxis()->SetRangeUser(0,45.5);
    histo2DRelStatErrEta->Draw("copy");

        TLegend* legendRelStatErrEta               = GetAndSetLegend2(0.14, 0.94-(0.035*nMeasSet*1.35), 0.45, 0.94, 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            if (availableMeasEta[i]== 2){
                DrawGammaSetMarker(statErrorRelCollectionEta[availableMeasEta[i]], markerStyleDet[availableMeasEta[i]], markerSizeDet[availableMeasEta[i]]*0.5, colorDet[availableMeasEta[i]] ,
                            colorDet[availableMeasEta[i]]);
                TGraphAsymmErrors* graphDummy   = new TGraphAsymmErrors(statErrorRelCollectionEta[availableMeasEta[i]]);
                DrawGammaSetMarkerTGraphAsym(graphDummy, markerStyleDet[availableMeasEta[i]], markerSizeDet[availableMeasEta[i]]*0.5, colorDet[availableMeasEta[i]],
                                    colorDet[availableMeasEta[i]]);
                graphDummy->Draw("same,p,x0");
                legendRelStatErrEta->AddEntry(graphDummy,nameMeasGlobal[availableMeasEta[i]],"p");

            } else {
                DrawGammaSetMarker(statErrorRelCollectionEta[availableMeasEta[i]], markerStyleDet[availableMeasEta[i]], markerSizeDet[availableMeasEta[i]]*0.5, colorDet[availableMeasEta[i]] ,
                            colorDet[availableMeasEta[i]]);
                statErrorRelCollectionEta[availableMeasEta[i]]->Draw("p,same,e1");
                legendRelStatErrEta->AddEntry(statErrorRelCollectionEta[availableMeasEta[i]],nameMeasGlobal[availableMeasEta[i]],"p");

            }
        }
        legendRelStatErrEta->Draw();

        TLatex *labelRelStatErrEnergyEta           = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
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

    TGraphAsymmErrors* graphCombEtaInvCrossSectionRelStat7TeV     = CalculateRelErrUpAsymmGraph( graphCombEtaInvCrossSectionStatPCMEMCPHOS, "graphCombEtaInvCrossSectionRelStat7TeV");
    TGraphAsymmErrors* graphCombEtaInvCrossSectionRelSys7TeV      = CalculateRelErrUpAsymmGraph( graphCombEtaInvCrossSectionSysPCMEMCPHOS, "graphCombEtaInvCrossSectionRelSys7TeV");
    TGraphAsymmErrors* graphCombEtaInvCrossSectionRelTot7TeV      = CalculateRelErrUpAsymmGraph( graphCombEtaInvCrossSectionTot, "graphCombEtaInvCrossSectionRelTot7TeV");

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

    DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionRelTot7TeV, markerStyleComb, markerSizeComb, colorComb , colorComb);
    graphCombEtaInvCrossSectionRelTot7TeV->Draw("p,same,e1");
    DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionRelStat7TeV, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
    graphCombEtaInvCrossSectionRelStat7TeV->Draw("l,x0,same,e1");
    DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionRelSys7TeV, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
    graphCombEtaInvCrossSectionRelSys7TeV->SetLineStyle(7);
    graphCombEtaInvCrossSectionRelSys7TeV->Draw("l,x0,same,e1");

    TLegend* legendRelTotErr2Eta = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
    legendRelTotErr2Eta->AddEntry(graphCombEtaInvCrossSectionRelTot7TeV,"tot","p");
    legendRelTotErr2Eta->AddEntry(graphCombEtaInvCrossSectionRelStat7TeV,"stat","l");
    legendRelTotErr2Eta->AddEntry(graphCombEtaInvCrossSectionRelSys7TeV,"sys","l");
    legendRelTotErr2Eta->Draw();

    TLatex *labelRelTotErrEnergyEta            = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
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
    Double_t weightsReadEtaToPi0[11][50];
    Int_t availableMeasEtaToPi0[11]                  = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    Int_t nPtBinsReadEtaToPi0                        = 0;
    while(!fileWeightsReadEtaToPi0.eof() && nPtBinsReadEtaToPi0 < 50){
        TString garbage                         = "";
        if (nPtBinsReadEtaToPi0 == 0){
            fileWeightsReadEtaToPi0 >> garbage ;
            for (Int_t i = 0; i < nMeasSet; i++){
                fileWeightsReadEtaToPi0 >> availableMeasEtaToPi0[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableMeasEtaToPi0[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsReadEtaToPi0 >> xValuesReadEtaToPi0[nPtBinsReadEtaToPi0-1];
            for (Int_t i = 0; i < nMeasSet; i++){
                fileWeightsReadEtaToPi0 >> weightsReadEtaToPi0[availableMeasEtaToPi0[i]][nPtBinsReadEtaToPi0-1] ;
            }
            cout << "read: "<<  nPtBinsReadEtaToPi0 << "\t"<< xValuesReadEtaToPi0[nPtBinsReadEtaToPi0-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSet; i++){
                cout << weightsReadEtaToPi0[availableMeasEtaToPi0[i]][nPtBinsReadEtaToPi0-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsReadEtaToPi0++;
    }
    nPtBinsReadEtaToPi0 = nPtBinsReadEtaToPi0-2 ;
    fileWeightsReadEtaToPi0.close();

    for (Int_t i = 0; i < nMeasSet; i++){
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

        TLegend* legendAccWeightsEtaToPi0               = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSet*1.35), 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaToPi0[availableMeasEtaToPi0[i]],
                                    markerStyleDet[availableMeasEtaToPi0[i]],
                                    markerSizeDet[availableMeasEtaToPi0[i]]*0.5,
                                    colorDet[availableMeasEtaToPi0[i]] ,
                                    colorDet[availableMeasEtaToPi0[i]]);
            graphWeightsEtaToPi0[availableMeasEtaToPi0[i]]->Draw("p,same,e1");
            legendAccWeightsEtaToPi0->AddEntry(graphWeightsEtaToPi0[availableMeasEtaToPi0[i]],nameMeasGlobal[availableMeasEtaToPi0[i]],"p");
        }
        legendAccWeightsEtaToPi0->Draw();
        TLatex *labelWeightsEnergyEtaToPi0              = new TLatex(0.7,0.20,collisionSystem7TeV.Data());
        SetStyleTLatex( labelWeightsEnergyEtaToPi0, 0.85*textSizeLabelsPixel,4);
        labelWeightsEnergyEtaToPi0->SetTextFont(43);
        labelWeightsEnergyEtaToPi0->Draw();
        TLatex *labelWeightsEtaToPi0                 = new TLatex(0.7,0.16,"#EtaToPi0 #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsEtaToPi0, 0.85*textSizeLabelsPixel,4);
        labelWeightsEtaToPi0->SetTextFont(43);
        labelWeightsEtaToPi0->Draw();

        DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);

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

        TLegend* legendRelSysErrEtaToPi0                = GetAndSetLegend2(0.62, 0.94-(0.035*nMeasSet*1.35), 0.95, 0.94, 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]], markerStyleDet[availableMeasEtaToPi0[i]], markerSizeDet[availableMeasEtaToPi0[i]]*0.5, colorDet[availableMeasEtaToPi0[i]],
                                    colorDet[availableMeasEtaToPi0[i]]);
            sysErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]]->Draw("p,same,e1");
            legendRelSysErrEtaToPi0->AddEntry(sysErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]],nameMeasGlobal[availableMeasEtaToPi0[i]],"p");
        }
        legendRelSysErrEtaToPi0->Draw();

        TLatex *labelRelSysErrEnergyEtaToPi0            = new TLatex(0.15,0.89,collisionSystem7TeV.Data());
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

        TLegend* legendRelStatErrEtaToPi0               = GetAndSetLegend2(0.14, 0.94-(0.035*nMeasSet*1.35), 0.45, 0.94, 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            if (availableMeasEtaToPi0[i]== 2){
                DrawGammaSetMarker(statErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]], markerStyleDet[availableMeasEtaToPi0[i]], markerSizeDet[availableMeasEtaToPi0[i]]*0.5, colorDet[availableMeasEtaToPi0[i]] ,
                            colorDet[availableMeasEtaToPi0[i]]);
                TGraphAsymmErrors* graphDummy   = new TGraphAsymmErrors(statErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]]);
                DrawGammaSetMarkerTGraphAsym(graphDummy, markerStyleDet[availableMeasEtaToPi0[i]], markerSizeDet[availableMeasEtaToPi0[i]]*0.5, colorDet[availableMeasEtaToPi0[i]],
                                    colorDet[availableMeasEtaToPi0[i]]);
                graphDummy->Draw("same,p,x0");
                legendRelStatErrEtaToPi0->AddEntry(graphDummy,nameMeasGlobal[availableMeasEtaToPi0[i]],"p");

            } else {
                DrawGammaSetMarker(statErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]], markerStyleDet[availableMeasEtaToPi0[i]], markerSizeDet[availableMeasEtaToPi0[i]]*0.5, colorDet[availableMeasEtaToPi0[i]] ,
                            colorDet[availableMeasEtaToPi0[i]]);
                statErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]]->Draw("p,same,e1");
                legendRelStatErrEtaToPi0->AddEntry(statErrorRelCollectionEtaToPi0[availableMeasEtaToPi0[i]],nameMeasGlobal[availableMeasEtaToPi0[i]],"p");

            }
        }
        legendRelStatErrEtaToPi0->Draw();

        TLatex *labelRelStatErrEnergyEtaToPi0           = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
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

    TGraphAsymmErrors* graphCombEtaToPi0InvCrossSectionRelStat7TeV     = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0StatPCMEMCPHOS, "graphCombEtaToPi0InvCrossSectionRelStat7TeV");
    TGraphAsymmErrors* graphCombEtaToPi0InvCrossSectionRelSys7TeV      = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0SysPCMEMCPHOS, "graphCombEtaToPi0InvCrossSectionRelSys7TeV");
    TGraphAsymmErrors* graphCombEtaToPi0InvCrossSectionRelTot7TeV      = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Tot, "graphCombEtaToPi0InvCrossSectionRelTot7TeV");

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

    DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0InvCrossSectionRelTot7TeV, markerStyleComb, markerSizeComb, colorComb , colorComb);
    graphCombEtaToPi0InvCrossSectionRelTot7TeV->Draw("p,same,e1");
    DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0InvCrossSectionRelStat7TeV, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
    graphCombEtaToPi0InvCrossSectionRelStat7TeV->Draw("l,x0,same,e1");
    DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0InvCrossSectionRelSys7TeV, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
    graphCombEtaToPi0InvCrossSectionRelSys7TeV->SetLineStyle(7);
    graphCombEtaToPi0InvCrossSectionRelSys7TeV->Draw("l,x0,same,e1");

    TLegend* legendRelTotErr2EtaToPi0 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
    legendRelTotErr2EtaToPi0->AddEntry(graphCombEtaToPi0InvCrossSectionRelTot7TeV,"tot","p");
    legendRelTotErr2EtaToPi0->AddEntry(graphCombEtaToPi0InvCrossSectionRelStat7TeV,"stat","l");
    legendRelTotErr2EtaToPi0->AddEntry(graphCombEtaToPi0InvCrossSectionRelSys7TeV,"sys","l");
    legendRelTotErr2EtaToPi0->Draw();

    TLatex *labelRelTotErrEnergyEtaToPi0            = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
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
    TGraphAsymmErrors* graphCombPi0InvCrossSectionTotUnShifted              = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTot->Clone("Unshifted");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionStatPCMEMCPHOSUnShifted   = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStatPCMEMCPHOS->Clone("UnshiftedStat");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionSysPCMEMCPHOSUnShifted    = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSysPCMEMCPHOS->Clone("UnshiftedSys");

    TGraphAsymmErrors* graphPi0InvCrossSectionStatUnShifted[11];
    TGraphAsymmErrors* graphPi0InvCrossSectionSysUnShifted[11];
    for (Int_t i = 0; i < 11; i++){
      if(directoryPi0[i]){
        graphPi0InvCrossSectionStatUnShifted[i] = (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[i]->Clone(Form("UnshiftedStat%s",nameMeasGlobal[i].Data()));
        graphPi0InvCrossSectionSysUnShifted[i]  = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i] ->Clone(Form("UnshiftedSys%s",nameMeasGlobal[i].Data()));
      }
    }

    // Calculating binshifts
    Double_t paramGraph[3]                      = {1.0e12, 8., 0.13};
    Double_t paramGraphEta[3]                      = {1.0e11, 8., 0.13};
    TF1* fitInvCrossSectionPi07TeV              = FitObject("l","fitInvCrossSectionPi07TeV","Pi0",histoPi0InvCrossSection[2],0.3,25.,paramGraph,"QNRMEX0+");
    TF1* fitInvCrossSectionEta7TeV              = FitObject("l","fitInvCrossSectionEta7TeV","Eta",histoEtaInvCrossSection[2],0.3,25.,paramGraphEta,"QNRMEX0+");
    TF1* fitInvCrossSectionPi07TeVGraph         = (TF1*)histoPi0InvCrossSection[2]->Clone("fitInvCrossSectionPi07TeVGraph");

    if(bWCorrection.CompareTo("X")==0 ){
        TF1* fitTsallisPi07TeVPtMult            = FitObject("tmpt","TsallisMultWithPtPi07TeV","Pi0");
        fitTsallisPi07TeVPtMult->SetParameters(paramGraph[0],paramGraph[1], paramGraph[2]) ; // standard parameter optimize if necessary
        graphCombPi0InvCrossSectionTot          = ApplyXshift(graphCombPi0InvCrossSectionTot, fitTsallisPi07TeVPtMult);
        for (Int_t i = 0; i < 11; i++){
          if(directoryPi0[i]){
            graphPi0InvCrossSectionStat[i]      = ApplyXshiftIndividualSpectra (graphCombPi0InvCrossSectionTot,
                                                                                        graphPi0InvCrossSectionStat[i],
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        0, graphPi0InvCrossSectionStat[i]->GetN());
            graphPi0InvCrossSectionSys[i]       = ApplyXshiftIndividualSpectra (graphCombPi0InvCrossSectionTot,
                                                                                        graphPi0InvCrossSectionSys[i],
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        0, graphPi0InvCrossSectionSys[i]->GetN());
          }
        }

        TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2,  0.13, 0.01, 0.015, 0.08);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F * histo2DDummy2;
        histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,minPtPi0,maxPtPi0,1000,1e-1,10e11);
        SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
        histo2DDummy2->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatPCMEMCPHOS, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
        graphCombPi0InvCrossSectionStatPCMEMCPHOS->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTot, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvCrossSectionTot->Draw("pEsame");

        fitInvCrossSectionPi07TeV->SetLineColor(kBlue+2);
        fitInvCrossSectionPi07TeV->Draw("same");

        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_7TeV.%s",outputDir.Data(),suffix.Data()));
    }

    graphCombPi0InvCrossSectionTot->Fit(fitInvCrossSectionPi07TeV,"QNRMEX0+","",0.4,50.);

//     fitInvCrossSectionPi07TeV = FitObject("l","fitInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot,0.4,50.,paramGraph,"QNRMEX0+");
    fitInvCrossSectionPi07TeV = FitObject("l","fitInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot,0.4,50. ,paramGraph,"QNRMEX0+");
    fitInvCrossSectionEta7TeV = FitObject("l","fitInvCrossSectionEta7TeV","Eta",graphCombEtaInvCrossSectionTot,0.4,30. ,paramGraphEta,"QNRMEX0+");

    cout << "pi0 fit:" << endl << WriteParameterToFile(fitInvCrossSectionPi07TeV)<< endl;
    cout << "eta fit:" << endl << WriteParameterToFile(fitInvCrossSectionEta7TeV)<< endl;
//     cout << "1: " << graphCombPi0InvCrossSectionTot->GetY()[2] << endl;0.1, csGraphs[2]->GetY()[4],0.6,3.0};
    Double_t paramTCM[5] = {graphCombPi0InvCrossSectionTot->GetY()[1],0.1,graphCombPi0InvCrossSectionTot->GetY()[4],0.6,3};
    Double_t paramTCMEta[5] = {graphCombEtaInvCrossSectionTot->GetY()[1],0.1,graphCombEtaInvCrossSectionTot->GetY()[4],0.6,3};
//     fitInvCrossSectionPi07TeV = FitObject("tcm","fitTCMInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot,0.4,50.,paramTCM,"QNRMEX0+");
//     fitInvCrossSectionPi07TeV = FitObject("l","fitPtLevy","Pi0",graphCombPi0InvCrossSectionTot,0.3,25,NULL,"QNRME+");
    TF1* fitTCMInvCrossSectionPi07TeV = FitObject("tcm","fitTCMInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot,0.3,25. ,paramTCM,"QNRMEX0+","", kFALSE);
    TF1* fitTCMInvCrossSectionEta7TeV = FitObject("tcm","graphCombEtaInvCrossSectionTotFit","Pi0",graphCombEtaInvCrossSectionTot,0.4,25.,paramTCMEta,"QNRMEX0+","", kFALSE);
//     fitTCMInvCrossSectionPi07TeV = FitObject("l","fitPtLevy","Pi0",graphCombPi0InvCrossSectionTot,0.4,25,NULL,"QNRME+");

//     cout << WriteParameterToFile(fitTCMInvCrossSectionPi07TeV)<< endl;

    TString forOutput= WriteParameterToFile(fitInvCrossSectionPi07TeV);
//     cout<< forOutput.Data()<< endl;


    TGraphAsymmErrors* graphRatioCombCombFitTot7TeV     = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTot->Clone();
    graphRatioCombCombFitTot7TeV                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitTot7TeV, fitTCMInvCrossSectionPi07TeV);
    TGraphAsymmErrors* graphRatioCombCombFitStat7TeV    = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStatPCMEMCPHOS->Clone();
    graphRatioCombCombFitStat7TeV                       = CalculateGraphErrRatioToFit(graphRatioCombCombFitStat7TeV, fitTCMInvCrossSectionPi07TeV);
    TGraphAsymmErrors* graphRatioCombCombFitSys7TeV     = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSysPCMEMCPHOS->Clone();
    graphRatioCombCombFitSys7TeV                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitSys7TeV, fitTCMInvCrossSectionPi07TeV);

    TGraphAsymmErrors* graphRatioCombCombFitTot7TeVEta     = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionTot->Clone();
    graphRatioCombCombFitTot7TeVEta                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitTot7TeVEta, fitInvCrossSectionEta7TeV);
    TGraphAsymmErrors* graphRatioCombCombFitStat7TeVEta    = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionStatPCMEMCPHOS->Clone();
    graphRatioCombCombFitStat7TeVEta                       = CalculateGraphErrRatioToFit(graphRatioCombCombFitStat7TeVEta, fitInvCrossSectionEta7TeV);
    TGraphAsymmErrors* graphRatioCombCombFitSys7TeVEta     = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionSysPCMEMCPHOS->Clone();
    graphRatioCombCombFitSys7TeVEta                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitSys7TeVEta, fitInvCrossSectionEta7TeV);

    TGraphAsymmErrors* graphRatioCombFitStat[11];
    TGraphAsymmErrors* graphRatioCombFitSys[11];
    for (Int_t i = 0; i < 11; i++){
      if(directoryPi0[i]){
        graphRatioCombFitStat[i]                = (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[i]->Clone();
        graphRatioCombFitStat[i]                = CalculateGraphErrRatioToFit(graphRatioCombFitStat[i], fitTCMInvCrossSectionPi07TeV);
        graphRatioCombFitSys[i]                 = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i]->Clone();
        graphRatioCombFitSys[i]                 = CalculateGraphErrRatioToFit(graphRatioCombFitSys[i], fitTCMInvCrossSectionPi07TeV);
      }
    }
    TGraphAsymmErrors* graphRatioCombFitStatEta[11];
    TGraphAsymmErrors* graphRatioCombFitSysEta[11];
    for (Int_t i = 0; i < 11; i++){
       if(directoryEta[i]){
            graphRatioCombFitStatEta[i]                = (TGraphAsymmErrors*)graphEtaInvCrossSectionStat[i]->Clone();
            graphRatioCombFitStatEta[i]                = CalculateGraphErrRatioToFit(graphRatioCombFitStatEta[i], fitInvCrossSectionEta7TeV);
            graphRatioCombFitSysEta[i]                 = (TGraphAsymmErrors*)graphEtaInvCrossSectionSys[i]->Clone();
            graphRatioCombFitSysEta[i]                 = CalculateGraphErrRatioToFit(graphRatioCombFitSysEta[i], fitInvCrossSectionEta7TeV);
       }
    }


    //  **********************************************************************************************************************
    //  ******************************************* Plot with Fit ****************************************
    //  **********************************************************************************************************************
    textSizeLabelsPixel                         = 48;
    TCanvas* canvasPlotwithFit                  = new TCanvas("canvasPlotwithFit","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasPlotwithFit, 0.12, 0.01, 0.01, 0.11);
    canvasPlotwithFit->SetLogx();
    canvasPlotwithFit->SetLogy();

        Double_t textsizeLabelsPP2              = 0;
        Double_t textsizeFacPP2                 = 0;
        if (canvasPlotwithFit->XtoPixel(canvasPlotwithFit->GetX2()) <canvasPlotwithFit->YtoPixel(canvasPlotwithFit->GetY1()) ){
            textsizeLabelsPP2 = (Double_t)textSizeLabelsPixel/canvasPlotwithFit->XtoPixel(canvasPlotwithFit->GetX2()) ;
            textsizeFacPP2 = (Double_t)1./canvasPlotwithFit->XtoPixel(canvasPlotwithFit->GetX2()) ;
        } else {
            textsizeLabelsPP2 = (Double_t)textSizeLabelsPixel/canvasPlotwithFit->YtoPixel(canvasPlotwithFit->GetY1());
            textsizeFacPP2 = (Double_t)1./canvasPlotwithFit->YtoPixel(canvasPlotwithFit->GetY1());
        }
        cout << textsizeLabelsPP2 << endl;

    TH2F * histo2DCombWithFit;
    histo2DCombWithFit                          = new TH2F("histo2DCombWithFit","histo2DCombWithFit",1000,minPtPi0,maxPtPi0,1000,0.01,9e13  );
    SetStyleHistoTH2ForGraphs(histo2DCombWithFit, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.85*textsizeLabelsPP2, textsizeLabelsPP2,
                            0.85*textsizeLabelsPP2,textsizeLabelsPP2, 0.9, 0.95, 510, 505);
    histo2DCombWithFit->GetXaxis()->SetMoreLogLabels();
    histo2DCombWithFit->GetXaxis()->SetNoExponent();
    histo2DCombWithFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTot, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombPi0InvCrossSectionTot->Draw("E2same");
//         fitInvCrossSectionPi07TeV->Draw("p,same,e1");
        fitTCMInvCrossSectionPi07TeV->Draw("p,same,e1");

    canvasPlotwithFit->SaveAs(Form("%s/CombWithFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));
    
    // ETA PLOT
    histo2DCombWithFit->Draw();

        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionTot, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombEtaInvCrossSectionTot->Draw("E2same");
//         fitInvCrossSectionEta7TeV->Draw("p,same,e1");
        fitTCMInvCrossSectionEta7TeV->Draw("p,same,e1");

    canvasPlotwithFit->SaveAs(Form("%s/Eta_CombWithFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));

    //    **********************************************************************************************************************
    //    ******************************************* Ratio of Comb to Fit ****************************************
    //    **********************************************************************************************************************
    textSizeLabelsPixel                         = 48;
    TCanvas* canvasRatioToCombFit               = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToCombFit, 0.12, 0.01, 0.01, 0.11);
    canvasRatioToCombFit->SetLogx();

        Double_t textsizeLabelsPP = 0;
        Double_t textsizeFacPP= 0;
        if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
            textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
            textsizeFacPP = (Double_t)1./canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
        } else {
            textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
            textsizeFacPP = (Double_t)1./canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
        }
        cout << textsizeLabelsPP << endl;

    TH2F * histo2DRatioToCombFit;
    histo2DRatioToCombFit                       = new TH2F("histo2DRatioToCombFit","histo2DRatioToCombFit",1000,minPtPi0,maxPtPi0,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DRatioToCombFit->GetXaxis()->SetNoExponent();
    histo2DRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DRatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys7TeV, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioCombCombFitSys7TeV->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat7TeV, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioCombCombFitStat7TeV->Draw("p,same,e1");

        DrawGammaLines(minPtPi0,maxPtPi0 , 1., 1.,0.1, kGray+2);
        DrawGammaLines(minPtPi0,maxPtPi0 , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(minPtPi0,maxPtPi0 , 0.9, 0.9,0.1, kGray, 7);

        TLatex *labelRatioToFitEnergy           = new TLatex(0.73,0.92,collisionSystem7TeV.Data());
        SetStyleTLatex( labelRatioToFitEnergy, 0.85*textSizeLabelsPixel,4);
        labelRatioToFitEnergy->SetTextFont(43);
        labelRatioToFitEnergy->Draw();
        TLatex *labelRatioToFitPi0              = new TLatex(0.73,0.87,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitPi0, 0.85*textSizeLabelsPixel,4);
        labelRatioToFitPi0->SetTextFont(43);
        labelRatioToFitPi0->Draw();


    canvasRatioToCombFit->SaveAs(Form("%s/RatioOfCombToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));
    
    // ETA PLOT
    histo2DRatioToCombFit->GetXaxis()->SetRangeUser(minPtEta,maxPtEta);
    histo2DRatioToCombFit->Draw();

        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys7TeVEta, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioCombCombFitSys7TeVEta->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat7TeVEta, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioCombCombFitStat7TeVEta->Draw("p,same,e1");

        DrawGammaLines(minPtEta,maxPtEta , 1., 1.,0.1, kGray+2);
        DrawGammaLines(minPtEta,maxPtEta , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(minPtEta,maxPtEta , 0.9, 0.9,0.1, kGray, 7);


        labelRatioToFitEnergy->Draw();
        TLatex *labelRatioToFitEta              = new TLatex(0.73,0.87,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitEta, 0.85*textSizeLabelsPixel,4);
        labelRatioToFitEta->SetTextFont(43);
        labelRatioToFitEta->Draw();


    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));

    //    **********************************************************************************************************************
    //    ******************************************* Ratio of Individual meas to Fit ******************************************
    //    **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DRatioToCombFit->GetXaxis()->SetRangeUser(minPtPi0,maxPtPi0);
    histo2DRatioToCombFit->Draw("copy");

    for (Int_t i = 0; i < 11; i++){
      if(directoryPi0[i]){
        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);

        graphRatioCombFitSys[i]->Draw("E2same");
        graphRatioCombFitStat[i]->Draw("p,same,e");
      }
    }

    DrawGammaLines(minPtPi0,maxPtPi0 , 1., 1.,0.5, kGray+2);
    DrawGammaLines(minPtPi0,maxPtPi0 , 1.1, 1.1,0.5, kGray, 7);
    DrawGammaLines(minPtPi0,maxPtPi0 , 0.9, 0.9,0.5, kGray, 7);

    labelRatioToFitEnergy->Draw();
    labelRatioToFitPi0->Draw();

    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    //         Double_t rowsLegendOnlyPi0Ratio[6]      = {0.92,0.88,0.84,0.80,0.76,0.72};
    Double_t rowsLegendOnlyPi0Ratio[6]      = {0.92,0.88,0.84,0.80,0.79,0.76};
    //         Double_t rowsLegendOnlyPi0RatioAbs[6]   = {0.91,2.2,2.1,2.0,1.9,1.8};
    Double_t rowsLegendOnlyPi0RatioAbs[6]   = {0.91,2.2,2.1,2.0,1.95,1.9};
    Double_t columnsLegendOnlyPi0Ratio[3]   = {0.15,0.32, 0.38};
    Double_t columnsLegendOnlyPi0RatioAbs[3]= {0.15,1.04, 1.37};
    Double_t lengthBox                      = 0.2/2;
    Double_t heightBox                      = 0.08/2;
    //****************** first Column **************************************************
    TLatex *textSingleMeasRatioPi0[11];
    for (Int_t i = 0; i < 11; i++){
      if(directoryPi0[i]){
        textSingleMeasRatioPi0[i]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[i+1],nameMeasGlobal[i].Data());
        SetStyleTLatex( textSingleMeasRatioPi0[i], 0.85*textSizeLabelsPixel,4);
        textSingleMeasRatioPi0[i]->SetTextFont(43);
        textSingleMeasRatioPi0[i]->Draw();
      }
    }

    //****************** second Column *************************************************
    TLatex *textStatOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[1]+0.03,rowsLegendOnlyPi0Ratio[0] ,"stat");
    SetStyleTLatex( textStatOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
    textStatOnlyRatioPi0->SetTextFont(43);
    textStatOnlyRatioPi0->Draw();
    TLatex *textSysOnlyRatioPi0             = new TLatex(columnsLegendOnlyPi0Ratio[2]+0.04 ,rowsLegendOnlyPi0Ratio[0],"syst");
    SetStyleTLatex( textSysOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
    textSysOnlyRatioPi0->SetTextFont(43);
    textSysOnlyRatioPi0->Draw();

    TMarker* markerPi0OnlyRatio[11];
    TBox* boxPi0OnlyRatio[11];

    for (Int_t i = 0; i < 11; i++){
      if(directoryPi0[i]){
        markerPi0OnlyRatio[i]               = CreateMarkerFromGraph(graphRatioCombFitSys[i],columnsLegendOnlyPi0Ratio[1]-0.1 ,rowsLegendOnlyPi0Ratio[i+1],1);
        markerPi0OnlyRatio[i]->DrawMarker(columnsLegendOnlyPi0RatioAbs[1]-0.1 ,rowsLegendOnlyPi0RatioAbs[i+1]);
        boxPi0OnlyRatio[i]                  = CreateBoxFromGraph(graphRatioCombFitSys[i], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox-0.1 , rowsLegendOnlyPi0RatioAbs[i+1]- heightBox,
        columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox-0.1, rowsLegendOnlyPi0RatioAbs[i+1]+ heightBox);
        boxPi0OnlyRatio[i]->Draw("l");
      }
    }

    canvasRatioToCombFit->SaveAs(Form("%s/RatioOfIndividualMeasToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));

    //    **********************************************************************************************************************
    //    ******************************************* Ratio of Individual meas to Fit ******************************************
    //    **********************************************************************************************************************

//    if(graphRatioCombFitSys[2] && graphRatioCombFitSys[8] && graphRatioCombFitStat[2] && graphRatioCombFitStat[8]){
//      canvasRatioToCombFit->cd();
//      TH2F * histo2DRatioToCombFitEviDaniel;
//      histo2DRatioToCombFitEviDaniel                       = new TH2F("histo2DRatioToCombFitEviDaniel","histo2DRatioToCombFitEviDaniel",1000,0.51,29.9,1000,0.2,4.    );
//      SetStyleHistoTH2ForGraphs(histo2DRatioToCombFitEviDaniel, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
//                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
//      histo2DRatioToCombFitEviDaniel->GetXaxis()->SetMoreLogLabels();
//      histo2DRatioToCombFitEviDaniel->GetXaxis()->SetNoExponent();
//      histo2DRatioToCombFitEviDaniel->GetYaxis()->SetRangeUser(0.61,1.65);
//      histo2DRatioToCombFitEviDaniel->Draw("copy");

//      graphRatioCombFitSys[2]->Draw("E2same");
//      graphRatioCombFitSys[8]->Draw("E2same");

//      graphRatioCombFitStat[2]->Draw("p,same,e");
//      graphRatioCombFitStat[8]->Draw("p,same,e");

//      DrawGammaLines(0.51, 29.9 , 1., 1.,0.5, kGray+2);
//      DrawGammaLines(0.51, 29.9 , 1.1, 1.1,0.5, kGray, 7);
//      DrawGammaLines(0.51, 29.9 , 0.9, 0.9,0.5, kGray, 7);

//      labelRatioToFitEnergy->Draw();
//      labelRatioToFitPi0->Draw();

//      textSingleMeasRatioPi0[2]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[0+1],"EMCal LHC10");
//      textSingleMeasRatioPi0[8]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1+1],"EMCal LHC11");
//      SetStyleTLatex( textSingleMeasRatioPi0[2], 0.85*textSizeLabelsPixel,4);
//      textSingleMeasRatioPi0[2]->SetTextFont(43);
//      SetStyleTLatex( textSingleMeasRatioPi0[8], 0.85*textSizeLabelsPixel,4);
//      textSingleMeasRatioPi0[8]->SetTextFont(43);
//      textSingleMeasRatioPi0[2]->Draw();
//      textSingleMeasRatioPi0[8]->Draw();

//      //****************** second Column *************************************************;
//      textStatOnlyRatioPi0->Draw();
//      textSysOnlyRatioPi0->Draw();
//      Double_t rowsLegendOnlyPi0RatioAbsEMC[6]   = {0.91,1.53,1.48,2.0,1.95,1.9};
//      Double_t columnsLegendOnlyPi0RatioAbsEMC[3]= {0.15,1.8, 2.2};
//      Double_t lengthBoxEMC                      = 0.3/2;
//      Double_t heightBoxEMC                      = 0.04/2;
//      markerPi0OnlyRatio[7]               = CreateMarkerFromGraph(graphRatioCombFitSys[2],columnsLegendOnlyPi0Ratio[1]-0.1 ,rowsLegendOnlyPi0Ratio[0+1],1);
//      markerPi0OnlyRatio[7]->DrawMarker(columnsLegendOnlyPi0RatioAbsEMC[1]-0.1 ,rowsLegendOnlyPi0RatioAbsEMC[0+1]);
//      markerPi0OnlyRatio[8]               = CreateMarkerFromGraph(graphRatioCombFitSys[8],columnsLegendOnlyPi0Ratio[1]-0.1 ,rowsLegendOnlyPi0Ratio[1+1],1);
//      markerPi0OnlyRatio[8]->DrawMarker(columnsLegendOnlyPi0RatioAbsEMC[1]-0.1 ,rowsLegendOnlyPi0RatioAbsEMC[1+1]-0.01);

//      boxPi0OnlyRatio[7]                  = CreateBoxFromGraph(graphRatioCombFitSys[2], columnsLegendOnlyPi0RatioAbsEMC[2]-0.5*lengthBoxEMC-0.1 , rowsLegendOnlyPi0RatioAbsEMC[0+1]- heightBoxEMC,
//          columnsLegendOnlyPi0RatioAbsEMC[2]+ 3*lengthBoxEMC-0.1, rowsLegendOnlyPi0RatioAbsEMC[0+1]+ heightBoxEMC);
//      boxPi0OnlyRatio[7]->Draw("l");
//      boxPi0OnlyRatio[8]                  = CreateBoxFromGraph(graphRatioCombFitSys[8], columnsLegendOnlyPi0RatioAbsEMC[2]-0.5*lengthBoxEMC-0.1 , rowsLegendOnlyPi0RatioAbsEMC[1+1]- heightBoxEMC,
//          columnsLegendOnlyPi0RatioAbsEMC[2]+ 3*lengthBoxEMC-0.1, rowsLegendOnlyPi0RatioAbsEMC[1+1]+ heightBoxEMC);
//      boxPi0OnlyRatio[8]->Draw("l");

//      canvasRatioToCombFit->SaveAs(Form("%s/RatioOfEMCMeasToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));
//    }

    //    **********************************************************************************************************************
    //    ******************* ETA ETA  ************* Ratio of Individual meas to Fit ***************** ETA  ETA ****************
    //    **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DRatioToCombFit->GetXaxis()->SetRangeUser(minPtEta,maxPtEta);
    histo2DRatioToCombFit->Draw("copy");

    for (Int_t i = 0; i < 11; i++){
      if(directoryEta[i]){
        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitSysEta[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitStatEta[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);

        graphRatioCombFitSysEta[i]->Draw("E2same");
        graphRatioCombFitStatEta[i]->Draw("p,same,e");
      }
    }

    DrawGammaLines(minPtEta,maxPtEta , 1., 1.,0.5, kGray+2);
    DrawGammaLines(minPtEta,maxPtEta , 1.1, 1.1,0.5, kGray, 7);
    DrawGammaLines(minPtEta,maxPtEta , 0.9, 0.9,0.5, kGray, 7);

    labelRatioToFitEnergy->Draw();
    labelRatioToFitEta->Draw();

    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************

    //****************** first Column **************************************************
    for (Int_t i = 0; i < 11; i++){
      if(directoryEta[i]){
        textSingleMeasRatioPi0[i]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[i+1],nameMeasGlobal[i].Data());
        SetStyleTLatex( textSingleMeasRatioPi0[i], 0.85*textSizeLabelsPixel,4);
        textSingleMeasRatioPi0[i]->SetTextFont(43);
        textSingleMeasRatioPi0[i]->Draw();
      }
    }

    //****************** second Column *************************************************
    textStatOnlyRatioPi0->Draw();
    textSysOnlyRatioPi0->Draw();

    for (Int_t i = 0; i < 11; i++){
      if(directoryEta[i]){
        markerPi0OnlyRatio[i]               = CreateMarkerFromGraph(graphRatioCombFitSysEta[i],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[i+1],1);
        markerPi0OnlyRatio[i]->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[i+1]);
        boxPi0OnlyRatio[i]                  = CreateBoxFromGraph(graphRatioCombFitSysEta[i], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[i+1]- heightBox,
        columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox+0.1, rowsLegendOnlyPi0RatioAbs[i+1]+ heightBox);
        boxPi0OnlyRatio[i]->Draw("l");
      }
    }

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));


//    if(graphRatioCombFitSysEta[2] && graphRatioCombFitSysEta[8] && graphRatioCombFitStatEta[2] && graphRatioCombFitStatEta[8]){
//      canvasRatioToCombFit->cd();
//      TH2F * histo2DRatioToCombFitEviDaniel;
//      histo2DRatioToCombFitEviDaniel                       = new TH2F("histo2DRatioToCombFitEviDaniel","histo2DRatioToCombFitEviDaniel",1000,0.51,29.9,1000,0.2,4.    );
//      SetStyleHistoTH2ForGraphs(histo2DRatioToCombFitEviDaniel, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
//                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
//      histo2DRatioToCombFitEviDaniel->GetXaxis()->SetMoreLogLabels();
//      histo2DRatioToCombFitEviDaniel->GetXaxis()->SetNoExponent();
//      histo2DRatioToCombFitEviDaniel->GetYaxis()->SetRangeUser(0.61,1.65);
//      histo2DRatioToCombFitEviDaniel->Draw("copy");

//      graphRatioCombFitSysEta[2]->Draw("E2same");
//      graphRatioCombFitSysEta[8]->Draw("E2same");

//      graphRatioCombFitStatEta[2]->Draw("p,same,e");
//      graphRatioCombFitStatEta[8]->Draw("p,same,e");

//      DrawGammaLines(0.51, 29.9 , 1., 1.,0.5, kGray+2);
//      DrawGammaLines(0.51, 29.9 , 1.1, 1.1,0.5, kGray, 7);
//      DrawGammaLines(0.51, 29.9 , 0.9, 0.9,0.5, kGray, 7);

//      labelRatioToFitEnergy->Draw();
//      labelRatioToFitEta->Draw();

//      textSingleMeasRatioPi0[2]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[0+1],"EMCal LHC10");
//      textSingleMeasRatioPi0[8]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1+1],"EMCal LHC11");
//      SetStyleTLatex( textSingleMeasRatioPi0[2], 0.85*textSizeLabelsPixel,4);
//      textSingleMeasRatioPi0[2]->SetTextFont(43);
//      SetStyleTLatex( textSingleMeasRatioPi0[8], 0.85*textSizeLabelsPixel,4);
//      textSingleMeasRatioPi0[8]->SetTextFont(43);
//      textSingleMeasRatioPi0[2]->Draw();
//      textSingleMeasRatioPi0[8]->Draw();

//      //****************** second Column *************************************************;
//      textStatOnlyRatioPi0->Draw();
//      textSysOnlyRatioPi0->Draw();
//      Double_t rowsLegendOnlyPi0RatioAbsEMC[6]   = {0.91,1.53,1.48,2.0,1.95,1.9};
//      Double_t columnsLegendOnlyPi0RatioAbsEMC[3]= {0.15,1.8, 2.2};
//      Double_t lengthBoxEMC                      = 0.3/2;
//      Double_t heightBoxEMC                      = 0.04/2;
//      markerPi0OnlyRatio[7]               = CreateMarkerFromGraph(graphRatioCombFitSysEta[2],columnsLegendOnlyPi0Ratio[1]-0.1 ,rowsLegendOnlyPi0Ratio[0+1],1);
//      markerPi0OnlyRatio[7]->DrawMarker(columnsLegendOnlyPi0RatioAbsEMC[1]-0.1 ,rowsLegendOnlyPi0RatioAbsEMC[0+1]);
//      markerPi0OnlyRatio[8]               = CreateMarkerFromGraph(graphRatioCombFitSysEta[8],columnsLegendOnlyPi0Ratio[1]-0.1 ,rowsLegendOnlyPi0Ratio[1+1],1);
//      markerPi0OnlyRatio[8]->DrawMarker(columnsLegendOnlyPi0RatioAbsEMC[1]-0.1 ,rowsLegendOnlyPi0RatioAbsEMC[1+1]-0.01);

//      boxPi0OnlyRatio[7]                  = CreateBoxFromGraph(graphRatioCombFitSysEta[2], columnsLegendOnlyPi0RatioAbsEMC[2]-0.5*lengthBoxEMC-0.1 , rowsLegendOnlyPi0RatioAbsEMC[0+1]- heightBoxEMC,
//          columnsLegendOnlyPi0RatioAbsEMC[2]+ 3*lengthBoxEMC-0.1, rowsLegendOnlyPi0RatioAbsEMC[0+1]+ heightBoxEMC);
//      boxPi0OnlyRatio[7]->Draw("l");
//      boxPi0OnlyRatio[8]                  = CreateBoxFromGraph(graphRatioCombFitSysEta[8], columnsLegendOnlyPi0RatioAbsEMC[2]-0.5*lengthBoxEMC-0.1 , rowsLegendOnlyPi0RatioAbsEMC[1+1]- heightBoxEMC,
//          columnsLegendOnlyPi0RatioAbsEMC[2]+ 3*lengthBoxEMC-0.1, rowsLegendOnlyPi0RatioAbsEMC[1+1]+ heightBoxEMC);
//      boxPi0OnlyRatio[8]->Draw("l");

//      canvasRatioToCombFit->SaveAs(Form("%s/RatioOfEMCMeasToCombFit_Eta_PP7TeV.%s",outputDir.Data(),suffix.Data()));
//    }

    TString nameOutputCommonFile                = Form("CombinedResultsPaperPP7TeV_%s.root", dateForOutput.Data());
    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    // PI0 MESON
    fCombResults.mkdir("Pi07TeV");
    fCombResults.cd("Pi07TeV");

    graphCombPi0InvCrossSectionStatPCMEMCPHOS   ->Write("graphInvCrossSectionPi0CombStat");
    graphCombPi0InvCrossSectionSysPCMEMCPHOS    ->Write("graphInvCrossSectionPi0CombSys");
    graphCombPi0InvCrossSectionStatPCMEMCPHOS   ->Write("graphInvCrossSectionPi0Comb7TeVAStatErr");
    graphCombPi0InvCrossSectionSysPCMEMCPHOS    ->Write("graphInvCrossSectionPi0Comb7TeVASysErr");
    
    TGraphAsymmErrors* graphCombPi0InvYieldStat = ScaleGraph( graphCombPi0InvCrossSectionStatPCMEMCPHOS,1/xSection7TeV/recalcBarn);
    TGraphAsymmErrors* graphCombPi0InvYieldSys = ScaleGraph( graphCombPi0InvCrossSectionSysPCMEMCPHOS,1/xSection7TeV/recalcBarn);
    graphCombPi0InvYieldStat   ->Write("graphCombPi0InvYieldStat");
    graphCombPi0InvYieldSys    ->Write("graphCombPi0InvYieldSys");
    
    for (Int_t i = 0; i < 11; i++){
      if(directoryPi0[i]){
        graphPi0InvCrossSectionStat[i]          ->Write(Form("graphInvCrossSectionPi0%sStat",nameMeasGlobal[i].Data()));
        graphPi0InvCrossSectionSys[i]           ->Write(Form("graphInvCrossSectionPi0%sSys",nameMeasGlobal[i].Data()));
        graphPi0InvCrossSectionStat[i]          ->Write(Form("graphInvCrossSectionPi0%s7TeVStatErr",nameMeasGlobal[i].Data()));
        graphPi0InvCrossSectionSys[i]           ->Write(Form("graphInvCrossSectionPi0%s7TeVSysErr",nameMeasGlobal[i].Data()));
      }
    }
    fitInvCrossSectionPi07TeV                   ->Write("TsallisFitPi0");
    fitTCMInvCrossSectionPi07TeV                ->Write("TwoComponentModelFitPi0");

    // ETA MESON
    fCombResults.mkdir("Eta7TeV");
    fCombResults.cd("Eta7TeV");

    graphCombEtaInvCrossSectionStatPCMEMCPHOS   ->Write("graphInvCrossSectionEtaCombStat");
    graphCombEtaInvCrossSectionSysPCMEMCPHOS    ->Write("graphInvCrossSectionEtaCombSys");
    graphCombEtaInvCrossSectionStatPCMEMCPHOS   ->Write("graphInvCrossSectionEtaComb7TeVAStatErr");
    graphCombEtaInvCrossSectionSysPCMEMCPHOS    ->Write("graphInvCrossSectionEtaComb7TeVASysErr");
    TGraphAsymmErrors* graphCombEtaInvYieldStat = ScaleGraph( graphCombEtaInvCrossSectionStatPCMEMCPHOS,1/xSection7TeV/recalcBarn);
    TGraphAsymmErrors* graphCombEtaInvYieldSys = ScaleGraph( graphCombEtaInvCrossSectionSysPCMEMCPHOS,1/xSection7TeV/recalcBarn);

    graphCombEtaInvYieldStat   ->Write("graphCombEtaInvYieldStat");
    graphCombEtaInvYieldSys    ->Write("graphCombEtaInvYieldSys");
    for (Int_t i = 0; i < 11; i++){
      if(directoryEta[i]){
       graphEtaInvCrossSectionStat[i]      ->Write(Form("graphInvCrossSectionEta%sStat",nameMeasGlobal[i].Data()));
       graphEtaInvCrossSectionSys[i]       ->Write(Form("graphInvCrossSectionEta%sSys",nameMeasGlobal[i].Data()));
       graphEtaInvCrossSectionStat[i]      ->Write(Form("graphInvCrossSectionEta%s7TeVStatErr",nameMeasGlobal[i].Data()));
       graphEtaInvCrossSectionSys[i]       ->Write(Form("graphInvCrossSectionEta%s7TeVSysErr",nameMeasGlobal[i].Data()));
      }
    }
    
    graphCombEtaToPi0StatPCMEMCPHOS             ->Write("graphEtaToPi0CombStat");
    graphCombEtaToPi0SysPCMEMCPHOS              ->Write("graphEtaToPi0CombSys");
    graphCombEtaToPi0StatPCMEMCPHOS             ->Write("graphRatioEtaToPi0Comb7TeVStatErr");
    graphCombEtaToPi0SysPCMEMCPHOS              ->Write("graphRatioEtaToPi0Comb7TeVSysErr");
    for (Int_t i = 0; i < 11; i++){
      if(directoryEta[i]){
        graphEtaToPi0Stat[i]                ->Write(Form("graphEtaToPi0%sStat",nameMeasGlobal[i].Data()));
        graphEtaToPi0Sys[i]                 ->Write(Form("graphEtaToPi0%sSys",nameMeasGlobal[i].Data()));
        graphEtaToPi0Stat[i]                ->Write(Form("graphRatioEtaToPi0%s7TeVStatErr",nameMeasGlobal[i].Data()));
        graphEtaToPi0Sys[i]                 ->Write(Form("graphRatioEtaToPi0%s7TeVSysErr",nameMeasGlobal[i].Data()));
      }
    }
    fitInvCrossSectionEta7TeV                   ->Write("TsallisFitEta");
    fitTCMInvCrossSectionPi07TeV                ->Write("TwoComponentModelFitEta");
    
    fCombResults.Close();
    

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

        TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, minPtPi0,maxPtPi0 ,1000., -30, 40);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
        histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,24.5);
        histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
        histo2DAllPi0FWHM->DrawCopy();

        for (Int_t i = 0; i < 11; i++){
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
        TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystem7TeV.Data());
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

        TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, minPtPi0,maxPtPi0, 1000., 125.1, 155.9);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                                  textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
        histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
        histo2DAllPi0Mass->DrawCopy();

        for (Int_t i = 0; i < 11; i++){
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
//         Double_t rowsLegendMass2[5] = {0.8,0.6,0.4,0.2,0.01};
//         Double_t rowsLegendMass2[6] = {0.84,0.66,0.50,0.33,0.16,0.01};
          Double_t  rowsLegendMass2[7]= {0.84,0.66,0.50,0.33,0.01,0.16};
        //******************* Offsets ***********************
        Double_t offsetMarkerXMass2         = 0.1;
        Double_t offsetMarkerYMass2         = 0.1;
        //****************** Scale factors ******************
        Double_t scaleMarkerMass2           = 1.2;

        padMassLegend1->cd();
        //****************** first Column **************************************************
        TLatex *textMassPCM[10];
        for (Int_t i = 0; i < 11; i++){
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
        for (Int_t i = 0; i < 11; i++){
            if(histoPi0Mass[i] && histoPi0TrueMass[i]){
                markerPCMPi0Mass[i]             = CreateMarkerFromHisto(histoPi0Mass[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
                markerPCMPi0Mass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
                markerPCMPi0MassMC[i]           = CreateMarkerFromHisto(histoPi0TrueMass[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
                markerPCMPi0MassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
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

        TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, minPtEta,maxPtEta ,1000., -1, 55.5);
        SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
//         histo2DAllEtaFWHM->GetYaxis()->SetRangeUser(-1.,45.5);
        histo2DAllEtaFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllEtaFWHM->GetYaxis()->SetNdivisions(505);
        histo2DAllEtaFWHM->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaFWHM->GetXaxis()->SetTickLength(0.05);
        histo2DAllEtaFWHM->GetYaxis()->SetTickLength(0.026);
        histo2DAllEtaFWHM->DrawCopy();

        for (Int_t i = 0; i <11; i++){
            if(histoEtaFWHMMeV[i] && histoEtaTrueFWHMMeV[i]){
                DrawGammaSetMarker(histoEtaFWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
                histoEtaFWHMMeV[i]->Draw("p,same,e");
                DrawGammaSetMarker(histoEtaTrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
                histoEtaTrueFWHMMeV[i]->Draw("p,same,e");
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

        TH2F * histo2DAllEtaMass            = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",20, minPtEta,maxPtEta, 1000., 515.1, 589.9);
        SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMassEta, textsizeLabelsMassEta, 0.85*textsizeLabelsMassEta,
                                  textsizeLabelsMassEta, 0.9, 0.28/(textsizeFacMassEta*margin), 512, 505);
        histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllEtaMass->GetYaxis()->SetNdivisions(505);
        histo2DAllEtaMass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaMass->GetXaxis()->SetTickLength(0.05);
        histo2DAllEtaMass->GetXaxis()->SetLabelOffset(-0.015);
        histo2DAllEtaMass->DrawCopy();

        for (Int_t i = 0; i < 11; i++){
            if(histoEtaMass[i] && histoEtaTrueMass[i]){
                DrawGammaSetMarker(histoEtaMass[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
                histoEtaMass[i]->Draw("p,same,e");
                DrawGammaSetMarker(histoEtaTrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
                histoEtaTrueMass[i]->Draw("p,same,e");
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
        for (Int_t i = 0; i < 11; i++){
            if(histoEtaMass[i] && histoEtaTrueMass[i] && histoEtaFWHMMeV[i] && histoEtaTrueFWHMMeV[i]){
                textMassPCMEta[i]                  = new TLatex(columnsLegendMass2Eta[0],rowsLegendMass2Eta[counterEta],nameMeasGlobal[i].Data());
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
        for (Int_t i = 0; i < 11; i++){
            if(histoEtaMass[i] && histoEtaTrueMass[i]){
                markerPCMEtaMass[i]             = CreateMarkerFromHisto(histoEtaMass[i],columnsLegendMass2Eta[1]+ offsetMarkerXMass2Eta ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta ,scaleMarkerMass2Eta);
                markerPCMEtaMass[i]->DrawMarker(columnsLegendMass2Eta[1]+ offsetMarkerXMass2Eta ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta);
                markerPCMEtaMassMC[i]           = CreateMarkerFromHisto(histoEtaTrueMass[i],columnsLegendMass2Eta[2]+ offsetMarkerXMass2Eta ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta ,scaleMarkerMass2Eta);
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
        histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, minPtPi0, maxPtPi0, 1000, 3e-6, 2 );
        SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);//(#times #epsilon_{pur})
        histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
        histo2DAccEff->GetXaxis()->SetNoExponent();
        histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAccEff->DrawCopy();

        for (Int_t i = 0; i < 11; i++){
            if(histoPi0AccTimesEff[i]){
                DrawGammaSetMarker(histoPi0AccTimesEff[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
                histoPi0AccTimesEff[i]->Draw("p,same,e");
            }
        }

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        for (Int_t i = 0; i < 11; i++){
            if(histoPi0AccTimesEff[i]){
                legendEffiAccPi0->AddEntry(histoPi0AccTimesEff[i],nameMeasGlobal[i].Data(),"p");
            }
        }
        legendEffiAccPi0->Draw();

        drawLatexAdd("ALICE performance",0.15,0.92,textSizeLabelsRel,kFALSE);
        drawLatexAdd(collisionSystem7TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
        drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));



    // ***************************************************************************************************************
    // ******************************* Plotting eta/pi0 ratio for single measurements ********************************
    // ***************************************************************************************************************
    cout << "PLOTTING: Eta/Pi0 ratio" << endl;
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

    TH2F * histo2DEtatoPi0combo;
    histo2DEtatoPi0combo               = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",1000,0.33,25.,1000,0.,1.05    );
    SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                              0.85*textsizeLabelsEtaToPi0,1.1*textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
    histo2DEtatoPi0combo->GetXaxis()->SetNoExponent();
    histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.0,1.05);
    histo2DEtatoPi0combo->Draw();
        // plotting systematics graphs
        for (Int_t i = 0; i < 11; i++){
            if(graphEtaToPi0Sys[i]){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Sys[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
                graphEtaToPi0Sys[i]->Draw("E2same");
            }
        }
         DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatPCMEMCPHOS, markerStyleDet[1], markerSizeDet[1]*0.75, kPink , kPink, widthLinesBoxes, kTRUE);
//                 graphCombEtaToPi0StatPCMEMCPHOS->Draw("E2same");
        // plotting statistics graphs
        for (Int_t i = 0; i < 11; i++){
            if(graphEtaToPi0Stat[i]){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Stat[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i]);
                graphEtaToPi0Stat[i]->Draw("p,same,e");
            }
        }
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0SysPCMEMCPHOS, markerStyleDet[1], markerSizeDet[1]*0.75, kPink , kPink);
//                 graphCombEtaToPi0SysPCMEMCPHOS->Draw("p,same,e");

        TLegend* legendEtaToPi0 = GetAndSetLegend2(0.67, 0.15, 0.9, 0.15+(textsizeLabelsEtaToPi0*4*0.9), textSizeLabelsPixel);
        for (Int_t i = 0; i < 11; i++){
            if(graphEtaToPi0Sys[i]){
                legendEtaToPi0->AddEntry(graphEtaToPi0Sys[i],nameMeasGlobal[i],"pf");
            }
        }
        legendEtaToPi0->Draw();

        drawLatexAdd(collisionSystem7TeV.Data(),0.13, 0.92,0.85*textsizeLabelsEtaToPi0,kFALSE);
        drawLatexAdd("ALICE",0.13, 0.92-(1*textsizeLabelsEtaToPi0*0.85),0.85*textsizeLabelsEtaToPi0,kFALSE);

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_differentSystems.%s",outputDir.Data(), suffix.Data()));



    // **********************************************************************************************************************
    // **************************Plot example invariant mass bins ***********************************************************
    // **********************************************************************************************************************

    textSizeLabelsPixel                         = 100*3/5;
    TCanvas* canvasInvMassSamplePlot            = new TCanvas("canvasInvMassSamplePlot","",0,0,1500,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasInvMassSamplePlot,  0.09, 0.01, 0.035, 0.08);

    Style_t markerStyleInvMassSGBG              = 0;
    Size_t markerSizeInvMassSGBG                = 0;
    Color_t markerColorInvMassSGBG              = kBlack;
    Style_t markerStyleInvMassMBG               = 24;
    Size_t markerSizeInvMassMBG                 = 1.5;
    Color_t markerColorInvMassMBG               = kGray+2;
    Style_t markerStyleInvMassBG                = 20;
    Size_t markerSizeInvMassBG                  = 2;
    Color_t markerColorInvMassBG                = kBlack;
    Style_t markerStyleInvMassSG                = 20;
    Size_t markerSizeInvMassSG                  = 3;
    Color_t markerColorInvMassSG                = kRed+2;
    Color_t fitColorInvMassSG                   = kAzure+2;

    Double_t marginInvMass                      = 0.1*1500;
    Double_t textsizeLabelsInvMass              = 0;
    Double_t textsizeFacInvMass                 = 0;

    TString strLowerEdgeExamplePi0[10]          = {"0.6","6.0","1.7","0.8","1.6"};
    TString strUpperEdgeExamplePi0[10]          = {"0.8","7.0","1.8","0.9","1.8"};
    if (canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) < canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1())){
        textsizeLabelsInvMass                   = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
        textsizeFacInvMass                      = (Double_t)1./canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
    } else {
        textsizeLabelsInvMass                   = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
        textsizeFacInvMass                      = (Double_t)1./canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
    }
    cout << textsizeLabelsInvMass << endl;

    TH2F * histo2DPi0InvMassDummy;
    histo2DPi0InvMassDummy                      = new TH2F("histo2DPi0InvMassDummy","histo2DPi0InvMassDummy",11000,0.05,0.235,21000,-1000,20000);
    SetStyleHistoTH2ForGraphs(histo2DPi0InvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                            0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));

    TH2F * histo2DEtaInvMassDummy;
    histo2DEtaInvMassDummy                      = new TH2F("histo2DEtaInvMassDummy","histo2DEtaInvMassDummy",11000,0.35,0.695,21000,-1000,20000);
    SetStyleHistoTH2ForGraphs(histo2DEtaInvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                            0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));

//     for (Int_t i =0 ; i < numbersofmeas; i++){
//         if (haveAllPi0InvMass[i]&&0){
//             canvasInvMassSamplePlot->cd();
//             histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.255);
//             histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSub[i]->GetMinimum(),1.15*histoPi0InvMassSigPlusBG[i]->GetMaximum());
//             if(i==2)
//                 histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSub[i]->GetMinimum(),1.35*histoPi0InvMassSigPlusBG[i]->GetMaximum());
//             histo2DPi0InvMassDummy->DrawCopy();
// 
//             TLatex *labelInvMassPtRange = new TLatex(0.945,0.9,Form("#pi^{0}: %s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}",strLowerEdgeExamplePi0[i].Data(),strUpperEdgeExamplePi0[i].Data()));
// 
//             DrawGammaSetMarker(histoPi0InvMassSigPlusBG[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
//             histoPi0InvMassSigPlusBG[i]->SetLineWidth(1);
//             histoPi0InvMassSigPlusBG[i]->Draw("hist,e,same");
//             DrawGammaSetMarker(histoPi0InvMassBGTot[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
//             histoPi0InvMassBGTot[i]->Draw("same");
// 
//             DrawGammaSetMarker(histoPi0InvMassSigRemBGSub[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
//             histoPi0InvMassSigRemBGSub[i]->Draw("same");
//             fitPi0InvMassSig[i]->SetNpx(1000);
//             fitPi0InvMassSig[i]->SetRange(0,0.255);
//             fitPi0InvMassSig[i]->SetLineColor(fitColorInvMassSG);
//             fitPi0InvMassSig[i]->SetLineWidth(1);
//             fitPi0InvMassSig[i]->Draw("same");
// 
//             TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9,collisionSystem7TeV.Data());
//             SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
//             labelInvMassEnergy->SetTextFont(43);
//             labelInvMassEnergy->Draw();
// 
//             TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.8*textsizeLabelsPP,"MinBias");
//             SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
//             labelInvMassTrigger->SetTextFont(43);
//             labelInvMassTrigger->Draw();
// 
//             TLatex *labelInvMassReco  = new TLatex(0.135,0.9-2*0.8*textsizeLabelsPP,Form("%s",nameMeasGlobal[i].Data()));
//             SetStyleTLatex( labelInvMassReco, 0.85*textSizeLabelsPixel,4);
//             labelInvMassReco->SetTextFont(43);
//             labelInvMassReco->Draw();
// 
//             SetStyleTLatex( labelInvMassPtRange, 0.85*textSizeLabelsPixel,4);
//             labelInvMassPtRange->SetTextAlign(31);
//             labelInvMassPtRange->SetTextFont(43);
//             labelInvMassPtRange->Draw();
// 
//             TLegend* legendInvMass  = GetAndSetLegend2(0.67, 0.88-5*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
//             legendInvMass->SetMargin(0.25);
//             legendInvMass->AddEntry(histoPi0InvMassSigPlusBG[i],"Raw real events","l");
//             if(i!=2){
//             legendInvMass->AddEntry(histoPi0InvMassBGTot[i],"Mixed event +","p");
//             legendInvMass->AddEntry((TObject*)0,"corr. BG","");
//             } else{
//             legendInvMass->AddEntry(histoPi0InvMassBGTot[i],"Mixed event","p");
//             }
//             legendInvMass->AddEntry(histoPi0InvMassSigRemBGSub[i],"BG subtracted","p");
//             legendInvMass->AddEntry(fitPi0InvMassSig[i], "Fit","l");
//             legendInvMass->Draw();
//             canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBin_%s.%s",outputDir.Data(), nameMeasGlobal[i].Data(), suffix.Data()));
//         } else {
//             cout << "missing partial input for invariant mass bin for  for trigger: " << nameMeasGlobal[i].Data() << endl;
//         }
//     }
    
   
}

















