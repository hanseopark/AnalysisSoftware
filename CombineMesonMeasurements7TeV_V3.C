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

void CombineMesonMeasurements7TeV_V3(   TString fileNamePCM     = "CombinationInput7TeV/data_PCMResultsFullCorrection_PP_7TeV_20170718.root",
                                        // TString fileNameEMCAL   = "/home/nschmidt/AnalysisSoftware/CombinationInput7TeV/data_EMCAL-EMCALResultsFullCorrection_PP_20170714.root",
                                        TString fileNameEMCAL   = "/home/nschmidt/AnalysisResults/pp/7TeV/EMCal/Evi/mesonSpecrta7TeV_2011EMCAL_14Oct2017.root",
                                        TString fileNameEMCAL2   = "/home/nschmidt/AnalysisSoftware/CombinationInput7TeV/data_EMCAL-EMCALResultsFullCorrection_PP_20170714.root",
                                        TString fileNamePHOS    = "/home/nschmidt/AnalysisResults/pp/7TeV/PHOS/pp7TeV_pass4_ppareek_PHOSResultsFullCorrection_10092017.root",
                                        TString fileNamePCMPHOS = "CombinationInput7TeV/data_PCM-PHOSResultsFullCorrection_PP_NoBinShifting_v2.root",
//                                         TString fileNamePCMEMCAL = "CombinationInput7TeV/data_PCM-EMCALResultsFullCorrection_PP.2.root",
                                        TString fileNamePCMEMCAL = "/home/nschmidt/AnalysisSoftware/CombinationInput7TeV/data_PCM-EMCALResultsFullCorrection_PP_20170714.root",
                                        TString fileInputCorrFactors = "/home/nschmidt/AnalysisResults/pp/7TeV/Comb/correlationInput/ComputeCorrelationFactors_pp7TeV/pp7TeV.root",
                                        TString suffix          = "pdf",
                                        TString isMC            = "",
                                        TString thesisPlots     = "",
                                        TString bWCorrection    ="",
                                        Int_t numbersofmeas     = 5,
                                        Bool_t useDanielmeas    = kFALSE
                                    ){

    TString date                                = ReturnDateString();

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
    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
    cout << outputDir.Data() << endl;
    cout << fileNamePCM.Data() << endl;

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
    Double_t xSection7TeVV0AND                  = 54.31*1e-3;
    Double_t xSection7TeVErrUp                  = 2.18;
    Double_t xSection7TeVErrDown                = 2.18;
    Double_t xSection7TeVppINEL                 = 73.2*1e9;
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

    TString nameMeasGlobal[11]                  = {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCMOtherDataset"};
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

    TFile* inputFile[10];
        inputFile[0]                            = new TFile(fileNamePCM.Data());
        inputFile[1]                            = new TFile(fileNamePHOS.Data());
        inputFile[2]                            = new TFile(fileNameEMCAL.Data());
        inputFile[3]                            = new TFile(fileNamePCMPHOS.Data());
        inputFile[4]                            = new TFile(fileNamePCMEMCAL.Data());
        inputFile[5]                            = new TFile(fileNameEMCAL2.Data());

    TDirectory* directoryPi0[10];
    TDirectory* directoryEta[10];
        for (Int_t i = 0; i < numbersofmeas+1; i++){
            directoryPi0[i]                     = (TDirectory*)inputFile[i]->Get("Pi07TeV");
            directoryEta[i]                     = (TDirectory*)inputFile[i]->Get("Eta7TeV");
        }
        
    TH1D* histoNumberOfEvents[10];
    TH1D* histoPi0Mass[10];
    TH1D* histoPi0FWHMMeV[10];
    TH1D* histoPi0TrueMass[10];
    TH1D* histoPi0TrueFWHMMeV[10];
    TH1D* histoEtaMass[10];
    TH1D* histoEtaFWHMMeV[10];
    TH1D* histoEtaTrueMass[10];
    TH1D* histoEtaTrueFWHMMeV[10];
    TH1D* histoPi0Acc[10];
    TH1D* histoPi0TrueEffPt[10];
    TH1D* histoPi0AccTimesEff[10];
    TH1D* histoPi0InvCrossSection[10];
    TGraphAsymmErrors* graphPi0InvCrossSectionSys[10];
    TGraphAsymmErrors* graphPi0InvCrossSectionStat[10];
    TH1D* histoEtaInvCrossSection[10];
    TGraphAsymmErrors* graphEtaInvCrossSectionSys[10];
    TGraphAsymmErrors* graphEtaInvCrossSectionStat[10];
    TH1D* histoPi0InvMassSigPlusBG[10];
    TH1D* histoPi0InvMassSig[10];
    TH1D* histoPi0InvMassSigRemBGSub[10];
    TH1D* histoPi0InvMassBG[10];
    TH1D* histoPi0InvMassRemBG[10];
    TH1D* histoPi0InvMassBGTot[10];
    TF1* fitPi0InvMassSig[10];
    TF1* fitPi0InvMassBG[10];
    Bool_t haveAllPi0InvMass[10]                    = {kFALSE, kFALSE, kFALSE,kFALSE,kFALSE};
    TString strInvMassBin[10]                       = {"04", "22", "3To3_2", "04",""};
//     TString strInvMassBin[10]                       = {"04", "22", "12", "04",""}; //DanielEMCAL

    TH1D* histoEtaToPi0Stat[10];
    TGraphAsymmErrors* graphEtaToPi0Stat[10]        = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaToPi0Sys[10]         = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    Double_t rapidityMeas[10]                       = {1.6, 1,1, 1.6,1.6,1,1};

    for (Int_t i = 0; i < numbersofmeas; i++){
        // LOAD PCM AND PHOS INPUT
        if(i!=0 && i!=4 && !(useDanielmeas&&i==2)){
            histoNumberOfEvents[i]              = (TH1D*)inputFile[i]->Get("histoNumberOfEvents7TeV");
            // load invariant mass peak positions and widths
            histoPi0Mass[i]                     = (TH1D*)directoryPi0[i]->Get("MassPi0");
            histoPi0FWHMMeV[i]                  = (TH1D*)directoryPi0[i]->Get("FWHMPi0MeV");
            histoPi0TrueMass[i]                 = (TH1D*)directoryPi0[i]->Get("TrueMassPi0");
            histoPi0TrueFWHMMeV[i]              = (TH1D*)directoryPi0[i]->Get("TrueFWHMPi0MeV");
            histoPi0InvCrossSection[i]          = (TH1D*)directoryPi0[i]->Get("InvCrossSectionPi0");
            if(i==3)
              histoPi0InvCrossSection[5]          = (TH1D*)directoryPi0[5]->Get("InvCrossSectionPi0");

            // load acceptance and efficiency and calculate acc*eff*y*2pi
            histoPi0Acc[i]                      = (TH1D*)directoryPi0[i]->Get("AcceptancePi0");
            histoPi0TrueEffPt[i]                = (TH1D*)directoryPi0[i]->Get("EfficiencyPi0");
            if(i!=2){
                histoPi0AccTimesEff[i]          = (TH1D*)histoPi0TrueEffPt[i]->Clone(Form("histoPi0AccTimesEff%s",nameMeasGlobal[i].Data()));
                histoPi0AccTimesEff[i]          ->Multiply(histoPi0Acc[i]);
                histoPi0AccTimesEff[i]          ->Scale(2*TMath::Pi()*rapidityMeas[i]);
            }else{
                histoPi0AccTimesEff[i]          = (TH1D*)directoryPi0[i]->Get("EfficiencyPi0");
            }
            // load invariant mass example bins
            histoPi0InvMassSig[i]               = (TH1D*)directoryPi0[i]->Get(Form("InvMassSig_PtBin%s",strInvMassBin[i].Data()));
            histoPi0InvMassSigPlusBG[i]         = (TH1D*)directoryPi0[i]->Get(Form("InvMassSigPlusBG_PtBin%s",strInvMassBin[i].Data()));
            histoPi0InvMassBG[i]                = (TH1D*)directoryPi0[i]->Get(Form("InvMassBG_PtBin%s",strInvMassBin[i].Data()));
            fitPi0InvMassSig[i]                 = (TF1*)directoryPi0[i]->Get(Form("FitInvMassSig_PtBin%s",strInvMassBin[i].Data()));
            if (histoPi0InvMassSig[i] && histoPi0InvMassSigPlusBG[i] && histoPi0InvMassBG[i] && fitPi0InvMassSig[i]){
                haveAllPi0InvMass[i]            = kTRUE;
            }
            // load cross section systematics and datapoints
            graphPi0InvCrossSectionSys[i]       = (TGraphAsymmErrors*)directoryPi0[i]->Get("InvCrossSectionPi0Sys");
            cout << nameMeasGlobal[i].Data() << " sys:" << endl;
//                 graphPi0InvCrossSectionSys[i]      ->Print();

            graphPi0InvCrossSectionStat[i]      = new TGraphAsymmErrors(histoPi0InvCrossSection[i]);
            if(i==3){
              graphPi0InvCrossSectionSys[5]       = (TGraphAsymmErrors*)directoryPi0[5]->Get("InvCrossSectionPi0Sys");
              graphPi0InvCrossSectionStat[5]      = new TGraphAsymmErrors(histoPi0InvCrossSection[5]);
                              graphPi0InvCrossSectionStat[i]      ->Print();
                              graphPi0InvCrossSectionSys[i]      ->Print();
            }
            if((i==0 || i==3) && graphPi0InvCrossSectionStat[i])
                graphPi0InvCrossSectionStat[i]  ->RemovePoint(0);
            cout << nameMeasGlobal[i].Data() << " stat:" << endl;
//                 graphPi0InvCrossSectionStat[i]      ->Print();

            // load eta invariant mass peak positions and widths
            if(directoryEta[i]){
                if(i==2||i==1){
                    histoEtaInvCrossSection[i]          = (TH1D*)directoryEta[i]->Get("InvCrossSectionEta");
                    // if(i==2 && !useDanielmeas)
                    //   histoEtaInvCrossSection[i]          = (TH1D*)directoryEta[i]->Get("InvariantCrossSectionEta");
                    graphEtaInvCrossSectionStat[i]      = new TGraphAsymmErrors(histoEtaInvCrossSection[i]);
                    cout << nameMeasGlobal[i].Data() << " stat:" << endl;
//                     graphEtaInvCrossSectionStat[i]      ->Print();
                    graphEtaInvCrossSectionSys[i]       = (TGraphAsymmErrors*)directoryEta[i]->Get("InvCrossSectionEtaSys");
                    cout << nameMeasGlobal[i].Data() << " sys:" << endl;
//                     graphEtaInvCrossSectionSys[i]      ->Print();
                    cout << __LINE__ << endl;
                }
                if(i==2){
                  cout << __LINE__ << endl;
                    histoEtaInvCrossSection[5]          = (TH1D*)directoryEta[5]->Get("InvCrossSectionEta");
                    graphEtaInvCrossSectionStat[5]      = new TGraphAsymmErrors(histoEtaInvCrossSection[5]);
                    graphEtaInvCrossSectionSys[5]       = (TGraphAsymmErrors*)directoryEta[5]->Get("InvCrossSectionEtaSys");
                    while(graphEtaInvCrossSectionSys[i]->GetX()[0] < 3.) graphEtaInvCrossSectionSys[i]->RemovePoint(0);
                    while(graphEtaInvCrossSectionStat[i]->GetX()[0] < 3.) graphEtaInvCrossSectionStat[i]->RemovePoint(0);
                    cout << __LINE__ << endl;
                }
            
                histoEtaMass[i]                 = (TH1D*)directoryEta[i]->Get("MassEta");
                histoEtaFWHMMeV[i]              = (TH1D*)directoryEta[i]->Get("FWHMEtaMeV");
                histoEtaTrueMass[i]             = (TH1D*)directoryEta[i]->Get("TrueMassEta");
                histoEtaTrueFWHMMeV[i]          = (TH1D*)directoryEta[i]->Get("TrueFWHMEtaMeV");
                if(i==2){
                    histoEtaToPi0Stat[i]            = (TH1D*)directoryEta[i]->Get("EtaToPi0RatioEMCal");
                    graphEtaToPi0Stat[i]            = new TGraphAsymmErrors(histoEtaToPi0Stat[i]);
                    cout << __LINE__ << endl;
//                     graphEtaToPi0Sys[i]            ->RemovePoint(0);
                    graphEtaToPi0Sys[i]             = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0RatioEMCalSys");
                    while(graphEtaToPi0Sys[i]->GetX()[graphEtaToPi0Sys[i]->GetN()-1] > 30.) graphEtaToPi0Sys[i]->RemovePoint(graphEtaToPi0Sys[i]->GetN()-1);
                    
                    
                    
                    graphEtaToPi0Stat[i]            ->Print();
                    graphEtaToPi0Sys[i]             ->Print();
                    cout << __LINE__ << endl;
                }
                if(i==1){
                    histoEtaToPi0Stat[i]            = (TH1D*)directoryEta[i]->Get("RatioEtaPi0");
                    graphEtaToPi0Stat[i]            = new TGraphAsymmErrors(histoEtaToPi0Stat[i]);
//                     graphEtaToPi0Sys[i]            ->RemovePoint(0);
                    graphEtaToPi0Sys[i]             = (TGraphAsymmErrors*)directoryEta[i]->Get("RatioEtaPi0Sys");
                    graphEtaToPi0Stat[i]            ->Print();
                    graphEtaToPi0Sys[i]             ->Print();
                }
                    
            }else{
                histoEtaMass[i]                 = NULL;
                histoEtaFWHMMeV[i]              = NULL;
                histoEtaTrueMass[i]             = NULL;
                histoEtaTrueFWHMMeV[i]          = NULL;
            }
        } else {
            // LOAD EMCAL AND PCM-EMCAL INPUT
            cout << "using daniel output" << endl;
            histoNumberOfEvents[i]              = (TH1D*)inputFile[i]->Get("histoNumberOfEvents7TeV");
            // load invariant mass peak positions and widths
            histoPi0Mass[i]                     = (TH1D*)directoryPi0[i]->Get("Pi0_Mass_data_INT1");
            histoPi0FWHMMeV[i]                  = (TH1D*)directoryPi0[i]->Get("Pi0_Width_data_INT1");
            histoPi0TrueMass[i]                 = (TH1D*)directoryPi0[i]->Get("Pi0_Mass_MC_INT1");
            histoPi0TrueFWHMMeV[i]              = (TH1D*)directoryPi0[i]->Get("Pi0_Width_MC_INT1");
            histoPi0Mass[i]                     ->Scale(1000);
            histoPi0FWHMMeV[i]                  ->Scale(1000);
            histoPi0TrueMass[i]                 ->Scale(1000);
            histoPi0TrueFWHMMeV[i]              ->Scale(1000);
            histoPi0InvCrossSection[i]          = (TH1D*)directoryPi0[i]->Get("InvCrossSectionPi0");

            // load acceptance and efficiency and calculate acc*eff*y*2pi
            histoPi0Acc[i]                      = (TH1D*)directoryPi0[i]->Get("AcceptancePi0_INT1");
            histoPi0TrueEffPt[i]                = (TH1D*)directoryPi0[i]->Get("EfficiencyPi0_INT1");
            histoPi0AccTimesEff[i]              = (TH1D*)histoPi0TrueEffPt[i]->Clone(Form("histoPi0AccTimesEff%s",nameMeasGlobal[i].Data()));
            histoPi0AccTimesEff[i]->Multiply(histoPi0Acc[i]);
            histoPi0AccTimesEff[i]->Scale(2*TMath::Pi()*rapidityMeas[i]);

            // load invariant mass example bins
            histoPi0InvMassSig[i]               = (TH1D*)directoryPi0[i]->Get("Pi0_InvMassSig_Example_INT1");
            histoPi0InvMassSigPlusBG[i]         = (TH1D*)directoryPi0[i]->Get("Pi0_InvMassSigPlusBG_Example_INT1");
            histoPi0InvMassBG[i]                = (TH1D*)directoryPi0[i]->Get("Pi0_InvMassBG_Example_INT1");
            fitPi0InvMassSig[i]                 = (TF1*)directoryPi0[i]->Get("Pi0_InvMassSigFit_Example_INT1");
            if (histoPi0InvMassSig[i] && histoPi0InvMassSigPlusBG[i] && histoPi0InvMassBG[i] && fitPi0InvMassSig[i]){
                haveAllPi0InvMass[i]            = kTRUE;
            }
            // load cross section systematics and datapoints
            graphPi0InvCrossSectionSys[i]       = (TGraphAsymmErrors*)directoryPi0[i]->Get("InvCrossSectionPi0Sys");
            cout << nameMeasGlobal[i].Data() << " sys:" << endl;
//                 graphPi0InvCrossSectionSys[i]      ->Print();

            graphPi0InvCrossSectionStat[i]      = new TGraphAsymmErrors(histoPi0InvCrossSection[i]);
            if((i==0 || i==3) && graphPi0InvCrossSectionStat[i])
                graphPi0InvCrossSectionStat[i]  ->RemovePoint(0);
            cout << nameMeasGlobal[i].Data() << " stat:" << endl;
//                 graphPi0InvCrossSectionStat[i]      ->Print();

            // load eta invariant mass peak positions and widths
            if(directoryEta[i]){
                histoEtaInvCrossSection[i]          = (TH1D*)directoryEta[i]->Get("InvCrossSectionEta");
                graphEtaInvCrossSectionStat[i]      = new TGraphAsymmErrors(histoEtaInvCrossSection[i]);
                if(graphEtaInvCrossSectionStat[i])
                    graphEtaInvCrossSectionStat[i]  ->RemovePoint(0);
                if((i==4) && graphEtaInvCrossSectionStat[i]){
                    graphEtaInvCrossSectionStat[i]  ->RemovePoint(0);
                    graphEtaInvCrossSectionStat[i]  ->RemovePoint(0);
                }
                cout << nameMeasGlobal[i].Data() << " stat:" << endl;
//                 graphEtaInvCrossSectionStat[i]      ->Print();
                graphEtaInvCrossSectionSys[i]       = (TGraphAsymmErrors*)directoryEta[i]->Get("InvCrossSectionEtaSys");
                cout << nameMeasGlobal[i].Data() << " sys:" << endl;
//                 graphEtaInvCrossSectionSys[i]      ->Print();
                histoEtaMass[i]                 = (TH1D*)directoryEta[i]->Get("Eta_Mass_data_INT1");
                histoEtaFWHMMeV[i]              = (TH1D*)directoryEta[i]->Get("Eta_Width_data_INT1");
                histoEtaTrueMass[i]             = (TH1D*)directoryEta[i]->Get("Eta_Mass_MC_INT1");
                histoEtaTrueFWHMMeV[i]          = (TH1D*)directoryEta[i]->Get("Eta_Width_MC_INT1");
                histoEtaMass[i]                 ->Scale(1000);
                histoEtaFWHMMeV[i]              ->Scale(1000);
                histoEtaTrueMass[i]             ->Scale(1000);
                histoEtaTrueFWHMMeV[i]          ->Scale(1000);
                histoEtaToPi0Stat[i]             = (TH1D*)directoryEta[i]->Get("EtaToPi0StatError");
                graphEtaToPi0Stat[i]             = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0StatError_INT1");
                graphEtaToPi0Sys[i]             = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0SystError_INT1");
                if(i==0||i==4||i==2){
                    histoEtaToPi0Stat[i]             = (TH1D*)directoryEta[i]->Get("EtaToPi0YShiftedStatError");
                    graphEtaToPi0Stat[i]             = (TGraphAsymmErrors*)directoryEta[i]->Get("graphEtaToPi0YShiftedStatError");
                    graphEtaToPi0Sys[i]             = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0YShiftedSystError_INT1");
//                     graphEtaToPi0Stat[i]            ->Print();
//                     graphEtaToPi0Sys[i]             ->Print();
                }
            }else{
                histoEtaMass[i]                 = NULL;
                histoEtaFWHMMeV[i]              = NULL;
                histoEtaTrueMass[i]             = NULL;
                histoEtaTrueFWHMMeV[i]          = NULL;
            }
        }
    }
    for (Int_t i = 0; i < numbersofmeas; i++){
        if (haveAllPi0InvMass[i]){
            histoPi0InvMassBGTot[i]                 = (TH1D*)histoPi0InvMassBG[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameMeasGlobal[i].Data()));
            histoPi0InvMassSigRemBGSub[i]           = (TH1D*)histoPi0InvMassSig[i]->Clone(Form("Pi0_InvMassSigRemBGSub_Example_%s",nameMeasGlobal[i].Data()));
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
        statErrorCollectionEta[i]                  = NULL;
        statErrorCollectionEtaToPi0[i]                  = NULL;
    }
    for (Int_t i = 0; i< numbersofmeas; i++){
        if(i!=3)
        statErrorCollection[i]                  = (TH1D*)histoPi0InvCrossSection[i]->Clone(Form("statErr%sPi0",nameMeasGlobal[i].Data()));
        if(i!=3)
        statErrorCollectionEta[i]                  = (TH1D*)histoEtaInvCrossSection[i]->Clone(Form("statErr%sEta",nameMeasGlobal[i].Data()));
        if(i!=3)
        statErrorCollectionEtaToPi0[i]                  = (TH1D*)histoEtaToPi0Stat[i]->Clone(Form("statErr%sEtaToPi0",nameMeasGlobal[i].Data()));
    }
    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollection[11];
    TGraphAsymmErrors* sysErrorCollectionEta[11];
    TGraphAsymmErrors* sysErrorCollectionEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollection[i]                   = NULL;
        sysErrorCollectionEta[i]                   = NULL;
        sysErrorCollectionEtaToPi0[i]                   = NULL;
    }
    for (Int_t i = 0; i< numbersofmeas; i++){
        if(i!=3 )
        sysErrorCollection[i]                   = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i]->Clone(Form("sysErr%sPi0",nameMeasGlobal[i].Data()));
        if(i!=3 )
        sysErrorCollectionEta[i]                   = (TGraphAsymmErrors*)graphEtaInvCrossSectionSys[i]->Clone(Form("sysErr%sEta",nameMeasGlobal[i].Data()));
        if(i!=3 )
        sysErrorCollectionEtaToPi0[i]                   = (TGraphAsymmErrors*)graphEtaToPi0Sys[i]->Clone(Form("sysErr%sEtaToPi0",nameMeasGlobal[i].Data()));
    }


    // Definition of final pt binning (has to be set manually)
    Double_t xPtLimits[51]                      =  { 0.0, 0.3, 0.4, 0.5, 0.6,
                                                     0.7, 0.8, 0.9, 1.0, 1.1,
                                                     1.2, 1.3, 1.4, 1.5, 1.6,
                                                     1.7, 1.8, 1.9, 2.0, 2.1,
                                                     2.2, 2.3, 2.4, 2.6, 2.8,
                                                     3.0, 3.2, 3.4, 3.6, 3.8,
                                                     4.0, 4.3, 4.6, 5.0, 5.5,
                                                     6.0, 6.5, 7.0, 8.0, 9.0,
                                                    10.0,11.0,12.0,13.0,14.0,16.0,18,
                                                    20.0,25.0
                                                    };
  // with daniel
    // Double_t xPtLimitsEta[51]                      =  { 0.0, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0,
    //                                                 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,14.,15., 16.0 ,18,
    //                                                 20.0,25.0,35.
    //                                                 };
  // with Evi
    Double_t xPtLimitsEta[51]                      =  { 0.0, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0,
                                                    3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,13.,14,15., 16.0 ,18,
                                                    20.0,25.0,35.
                                                    };
  // with daniel
    // Double_t xPtLimitsEtaToPi0[51]                      =  { 0.0, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0,
    //                                                 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,14.,15., 16.0 ,18,
    //                                                 20.0,25.0,35.
    //                                                 };
  // with Evi
    Double_t xPtLimitsEtaToPi0[50]                      =  { 0.0, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0,
                                                    3.5, 4.0, 5.0, 6.0, 8.0,10,12.0,13.0,14.0,15.,16.0,18,
                                                    20.0,25.0
                                                    };


    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match
//     Int_t offSets[11]                           =  {0,   8,     1,      2,        0,           0,   0,     0,      0,        0, 0};
//     Int_t offSetsSys[11]                        =  {1,   8,    11,      3,        6,           0,   0,     0,      2,        0, 0};
    Int_t offSets[11]                           =  {0,   8,     1,      2,        0,           0,   0,     0,      0,        0, 0};
    Int_t offSetsSys[11]                        =  {1,   8,    11,      3,        6,           0,   0,     0,      2,        0, 0};
    if(!useDanielmeas){
      offSets[2]    = 6;//4
      offSetsSys[2] = 6;//4
    }
                                                // pcm, phos, emcal, pcmphos, pcmemcal
                                                //  0    1      2       3         4
    Int_t offSetsEta[11]                           =  {0,   4,     1,      2,        1,           0,   0,     0,      0,        0, 0};
    Int_t offSetsSysEta[11]                        =  {1,   4,    7,      3,        4,           0,   0,     0,      2,        0, 0};
    if(!useDanielmeas){
      offSetsEta[2]    = 4;
      offSetsSysEta[2] = 9;
    }
                                                    // pcm, phos, emcal, pcmphos, pcmemcal
                                                    //  0    1      2       3         4
    Int_t offSetsEtaToPi0[11]                      =  {0,   4,     1,      2,        1,           0,   0,     0,      0,        0, 0};
    Int_t offSetsSysEtaToPi0[11]                   =  {1,   4,    7,      3,        4,           0,   0,     0,      2,        0, 0};
    if(!useDanielmeas){
      offSetsEtaToPi0[2]    = 4;
      offSetsSysEtaToPi0[2] = 4;
    }
    
    Int_t nBinsPi0 = 48;
    if(!useDanielmeas)
      nBinsPi0 = 48;
    Int_t nBinsEta = 17;
    if(!useDanielmeas)
      nBinsEta = 24;

    //    **********************************************************************************************************************
    //    ******************************************* Calculation of spectrum including EMCal only *****************************
    //    **********************************************************************************************************************

    TH1D* statErrorRelCollection[11];
    TH1D* statErrorRelCollectionEta[11];
    TH1D* statErrorRelCollectionEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollection[i] = NULL;
        statErrorRelCollectionEta[i] = NULL;
        statErrorRelCollectionEtaToPi0[i] = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (statErrorCollection[i]) statErrorRelCollection[i] = CalculateRelErrUpTH1D( statErrorCollection[i], Form("relativeStatError_%s", nameMeasGlobal[i].Data()));
        if (statErrorCollectionEta[i]) statErrorRelCollection[i] = CalculateRelErrUpTH1D( statErrorCollectionEta[i], Form("relativeStatError_%s", nameMeasGlobal[i].Data()));
        if (statErrorCollectionEtaToPi0[i]) statErrorRelCollection[i] = CalculateRelErrUpTH1D( statErrorCollectionEtaToPi0[i], Form("relativeStatError_%s", nameMeasGlobal[i].Data()));
    }

    TGraphAsymmErrors* sysErrorRelCollection[11];
    TGraphAsymmErrors* sysErrorRelCollectionEta[11];
    TGraphAsymmErrors* sysErrorRelCollectionEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollection[i] = NULL;
        sysErrorRelCollectionEta[i] = NULL;
        sysErrorRelCollectionEtaToPi0[i] = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (sysErrorCollection[i]) sysErrorRelCollection[i] = CalculateRelErrUpAsymmGraph( sysErrorCollection[i], Form("relativeSysError_%s", nameMeasGlobal[i].Data()));
        if (sysErrorCollectionEta[i]) sysErrorRelCollectionEta[i] = CalculateRelErrUpAsymmGraph( sysErrorCollectionEta[i], Form("relativeSysError_%s", nameMeasGlobal[i].Data()));
        if (sysErrorCollectionEtaToPi0[i]) sysErrorRelCollectionEta[i] = CalculateRelErrUpAsymmGraph( sysErrorCollectionEtaToPi0[i], Form("relativeSysError_%s", nameMeasGlobal[i].Data()));
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
    TString fileNameOutputWeightingEta                       = Form("%s/WeightingEta.dat",outputDir.Data());
    TString fileNameOutputWeightingEtaToPi0                       = Form("%s/WeightingEtaToPi0.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombPi0InvCrossSectionStatPCMEMCPHOS= NULL;
    TGraphAsymmErrors* graphCombPi0InvCrossSectionSysPCMEMCPHOS = NULL;
//     TGraphAsymmErrors* graphCombPi0InvCrossSectionTot = NULL;
    TGraphAsymmErrors* graphCombPi0InvCrossSectionTot = CombinePtPointsSpectraFullCorrMat( statErrorCollection, sysErrorCollection,
                                                                                                   xPtLimits, nBinsPi0,
                                                                                                   offSets, offSetsSys,
                                                                                                   graphCombPi0InvCrossSectionStatPCMEMCPHOS, graphCombPi0InvCrossSectionSysPCMEMCPHOS,
//                                                                                                    fileNameOutputWeighting,1
                                                                                                   fileNameOutputWeighting,"7TeV", "Pi0", kTRUE,
                                                                                                0x0, fileInputCorrFactors
                                                                                                );

    // return;

    graphCombPi0InvCrossSectionStatPCMEMCPHOS->Print();
    
    TGraphAsymmErrors* graphCombEtaInvCrossSectionStatPCMEMCPHOS= NULL;
    TGraphAsymmErrors* graphCombEtaInvCrossSectionSysPCMEMCPHOS = NULL;
    TGraphAsymmErrors* graphCombEtaInvCrossSectionTot = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEta, sysErrorCollectionEta,
                                                                                                   xPtLimitsEta, nBinsEta,
                                                                                                   offSetsEta, offSetsSysEta,
                                                                                                   graphCombEtaInvCrossSectionStatPCMEMCPHOS, graphCombEtaInvCrossSectionSysPCMEMCPHOS,
//                                                                                                    fileNameOutputWeightingEta,1
                                                                                                   fileNameOutputWeightingEta,"7TeV", "Eta", kTRUE,
                                                                                                   0x0, fileInputCorrFactors
                                                                                                );
    // return;

    graphCombEtaInvCrossSectionStatPCMEMCPHOS->Print();
    
    
    TGraphAsymmErrors* graphCombEtaToPi0StatPCMEMCPHOS= NULL;
    TGraphAsymmErrors* graphCombEtaToPi0SysPCMEMCPHOS = NULL;
    TGraphAsymmErrors* graphCombEtaToPi0Tot = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEtaToPi0, sysErrorCollectionEtaToPi0,
                                                                                                   xPtLimitsEtaToPi0, nBinsEta-1,
                                                                                                   offSetsEtaToPi0, offSetsSysEtaToPi0,
                                                                                                   graphCombEtaToPi0StatPCMEMCPHOS, graphCombEtaToPi0SysPCMEMCPHOS,
//                                                                                                    fileNameOutputWeightingEtaToPi0,1
                                                                                                   fileNameOutputWeightingEtaToPi0,"7TeV", "EtaToPi0", kTRUE,
                                                                                                   0x0, fileInputCorrFactors
                                                                                                );
    // return;

    graphCombEtaToPi0StatPCMEMCPHOS->Print();
    
    
    

    // Reading weights from output file for plotting
    ifstream fileWeightsRead;
    fileWeightsRead.open(fileNameOutputWeighting,ios_base::in);
    cout << "reading" << fileNameOutputWeighting << endl;
    Double_t xValuesRead[50];
    Double_t weightsRead[11][50];
    Int_t availableMeas[11]                     = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    Int_t nMeasSet                              = 4;
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

    //    **********************************************************************************************************************
    //    ******************************************* Plotting weights method only EMC *****************************************
    //    **********************************************************************************************************************
    Int_t textSizeLabelsPixel                   = 900*0.04;

    TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    canvasWeights->SetLogx();

    TH2F * histo2DWeights;
    histo2DWeights                              = new TH2F("histo2DWeights","histo2DWeights",11000,0.23,70.,1000,-0.5,1.1);
    SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DWeights->GetXaxis()->SetMoreLogLabels();
    histo2DWeights->GetXaxis()->SetLabelOffset(-0.01);
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

    canvasWeights->SaveAs(Form("%s/WeightsMethods.%s",outputDir.Data(),suffix.Data()));



    //    **********************************************************************************************************************
    //    ******************************************* Ratio of Comb to Fit ****************************************
    //    **********************************************************************************************************************

    TCanvas* canvasRatioToOldCombined           = new TCanvas("canvasRatioToOldCombined","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToOldCombined, 0.1, 0.02, 0.035, 0.09);
    canvasRatioToOldCombined->SetLogx();

    TH2F * histo2DRatioToOldCombined;
    histo2DRatioToOldCombined                   = new TH2F("histo2DRatioToOldCombined","histo2DRatioToOldCombined",11000,0.23,10.,1000,0.75,1.25);
    SetStyleHistoTH2ForGraphs(histo2DRatioToOldCombined, "#it{p}_{T} (GeV/#it{c})","#frac{Comb}{Comb Old}",0.035,0.04, 0.035,0.04, 1.,1.,510,505);
    histo2DRatioToOldCombined->GetXaxis()->SetMoreLogLabels();
    histo2DRatioToOldCombined->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRatioToOldCombined->Draw("copy");

        DrawGammaLines(0.23, 10. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 10. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 10. , 0.9, 0.9,0.1, kGray, 7);
        DrawGammaLines(0.23, 10. , 1.05, 1.05,0.1, kGray, 3);
        DrawGammaLines(0.23, 10. , 0.95, 0.95,0.1, kGray, 3);

        TLatex *labelRatioToOldEnergy           = new TLatex(0.15,0.89,collisionSystem7TeV.Data());
        SetStyleTLatex( labelRatioToOldEnergy, 0.85*textSizeLabelsPixel,4);
        labelRatioToOldEnergy->SetTextFont(43);
        labelRatioToOldEnergy->Draw();
        TLatex *labelRatioToOldPi0              = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToOldPi0, 0.85*textSizeLabelsPixel,4);
        labelRatioToOldPi0->SetTextFont(43);
        labelRatioToOldPi0->Draw();

    canvasRatioToOldCombined->SaveAs(Form("%s/RatioOfCombToCombOld_PP7TeV.%s",outputDir.Data(),suffix.Data()));


    //     *********************************************************************************************************************
    //     ************************************ Visualize relative errors ******************************************************
    //     *********************************************************************************************************************

    TCanvas* canvasRelSysErr                    = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();

    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                            = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23,70.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
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

    canvasRelSysErr->SaveAs(Form("%s/RelSysErr.%s",outputDir.Data(),suffix.Data()));

    //     *********************************************************************************************************************
    //     ************************************ Visualize relative errors ******************************************************
    //     *********************************************************************************************************************

    TCanvas* canvasRelStatErr                   = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();

    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                           = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23,70.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
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

    canvasRelStatErr->SaveAs(Form("%s/RelStatErr.%s",outputDir.Data(),suffix.Data()));


    //************************************************************************************************************************
    //************************************** Comparison sys and stat for new and old combined ********************************
    //************************************************************************************************************************
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelStat7TeV     = CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionStatPCMEMCPHOS, "relativeStatError_");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelSys7TeV      = CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionSysPCMEMCPHOS, "relativeSysError_");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelTot7TeV      = CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionTot, "relativeTotalError_");

    TCanvas* canvasRelTotErr                    = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();

    TH2F * histo2DRelTotErr;
    histo2DRelTotErr                            = new TH2F("histo2DRelTotErr","histo2DRelTotErr",11000,0.23,70.,1000,0,57.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErr, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelTotErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelTot7TeV, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
        graphCombPi0InvCrossSectionRelTot7TeV->Draw("p,same,e1");

        TLegend* legendRelTotErr                = GetAndSetLegend2(0.14, 0.94-(0.035*4*1.35), 0.45, 0.94, 32);
        legendRelTotErr->AddEntry(graphCombPi0InvCrossSectionRelTot7TeV,"PCM, PHOS, EMCAL","p");
        legendRelTotErr->Draw();

        TLatex *labelRelTotErrEnergy            = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
        SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrEnergy->SetTextFont(43);
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrPi0               = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelTotErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrPi0->SetTextFont(43);
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/RelTotErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelSysErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelSys7TeV, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
        graphCombPi0InvCrossSectionRelSys7TeV->Draw("p,same,e1");

        legendRelTotErr->Draw();
        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/RelSysErr_diffComb.%s",outputDir.Data(),suffix.Data()));

    histo2DRelStatErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelStat7TeV, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
        graphCombPi0InvCrossSectionRelStat7TeV->Draw("p,same,e1");

        legendRelTotErr->Draw();
        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/RelStatErr_diffComb.%s",outputDir.Data(),suffix.Data()));

    histo2DRelTotErr->GetYaxis()->SetRangeUser(0,57.5);
    histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErr->Draw("copy");

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

        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/Reldecomp.%s",outputDir.Data(),suffix.Data()));



    //    **********************************************************************************************************************
    //    ************************************* Calculating bin shifted spectra & fitting **************************************
    //    **********************************************************************************************************************

    // Cloning spectra
    TGraphAsymmErrors* graphCombPi0InvCrossSectionTotUnShifted              = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTot->Clone("Unshifted");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionStatPCMEMCPHOSUnShifted   = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStatPCMEMCPHOS->Clone("UnshiftedStat");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionSysPCMEMCPHOSUnShifted    = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSysPCMEMCPHOS->Clone("UnshiftedSys");

    TGraphAsymmErrors* graphPi0InvCrossSectionStatUnShifted[10];
    TGraphAsymmErrors* graphPi0InvCrossSectionSysUnShifted[10];
    for (Int_t i = 0; i < 4; i++){
        graphPi0InvCrossSectionStatUnShifted[i] = (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[i]->Clone(Form("UnshiftedStat%s",nameMeasGlobal[i].Data()));
        graphPi0InvCrossSectionSysUnShifted[i]  = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i] ->Clone(Form("UnshiftedSys%s",nameMeasGlobal[i].Data()));
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
        for (Int_t i = 0; i < 4; i++){
            graphPi0InvCrossSectionStat[i]      = ApplyXshiftIndividualSpectra (graphCombPi0InvCrossSectionTot,
                                                                                        graphPi0InvCrossSectionStat[i],
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        0, graphPi0InvCrossSectionStat[i]->GetN());
            graphPi0InvCrossSectionSys[i]       = ApplyXshiftIndividualSpectra (graphCombPi0InvCrossSectionTot,
                                                                                        graphPi0InvCrossSectionSys[i],
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        0, graphPi0InvCrossSectionSys[i]->GetN());
        }

        TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2,  0.13, 0.01, 0.015, 0.08);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F * histo2DDummy2;
        histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,70.,1000,1e-1,10e11);
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
    Double_t paramTCM[5] = {graphCombPi0InvCrossSectionStatPCMEMCPHOS->GetY()[1],0.1,graphCombPi0InvCrossSectionStatPCMEMCPHOS->GetY()[4],0.6,3};
    Double_t paramTCMEta[5] = {graphCombEtaInvCrossSectionStatPCMEMCPHOS->GetY()[1],0.1,graphCombEtaInvCrossSectionStatPCMEMCPHOS->GetY()[4],0.6,3};
//     fitInvCrossSectionPi07TeV = FitObject("tcm","fitTCMInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot,0.4,50.,paramTCM,"QNRMEX0+");
//     fitInvCrossSectionPi07TeV = FitObject("l","fitPtLevy","Pi0",graphCombPi0InvCrossSectionTot,0.3,25,NULL,"QNRME+");
    TF1* fitTCMInvCrossSectionPi07TeV = FitObject("tcm","fitTCMInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionStatPCMEMCPHOS,0.3,25. ,paramTCM,"QNRMEX0+","", kFALSE);
    TF1* fitTCMInvCrossSectionEta7TeV = FitObject("tcm","graphCombEtaInvCrossSectionTotFit","Pi0",graphCombEtaInvCrossSectionStatPCMEMCPHOS,0.4,25.,paramTCMEta,"QNRMEX0+","", kFALSE);
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

    TGraphAsymmErrors* graphRatioCombFitStat[10];
    TGraphAsymmErrors* graphRatioCombFitSys[10];
    for (Int_t i = 0; i < numbersofmeas+1; i++){
        graphRatioCombFitStat[i]                = (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[i]->Clone();
        graphRatioCombFitStat[i]                = CalculateGraphErrRatioToFit(graphRatioCombFitStat[i], fitTCMInvCrossSectionPi07TeV);
        graphRatioCombFitSys[i]                 = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i]->Clone();
        graphRatioCombFitSys[i]                 = CalculateGraphErrRatioToFit(graphRatioCombFitSys[i], fitTCMInvCrossSectionPi07TeV);
    }
    TGraphAsymmErrors* graphRatioCombFitStatEta[10];
    TGraphAsymmErrors* graphRatioCombFitSysEta[10];
    for (Int_t i = 0; i < numbersofmeas+1; i++){
        if(i!=3){
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
    histo2DCombWithFit                          = new TH2F("histo2DCombWithFit","histo2DCombWithFit",1000,0.23,70.,1000,0.01,9e13  );
    SetStyleHistoTH2ForGraphs(histo2DCombWithFit, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.85*textsizeLabelsPP2, textsizeLabelsPP2,
                            0.85*textsizeLabelsPP2,textsizeLabelsPP2, 0.9, 0.95, 510, 505);
    histo2DCombWithFit->GetXaxis()->SetMoreLogLabels();
    histo2DCombWithFit->GetXaxis()->SetLabelOffset(-0.01);
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
    histo2DRatioToCombFit                       = new TH2F("histo2DRatioToCombFit","histo2DRatioToCombFit",1000,0.23,70.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DRatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys7TeV, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioCombCombFitSys7TeV->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat7TeV, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioCombCombFitStat7TeV->Draw("p,same,e1");

        DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);

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
    histo2DRatioToCombFit->GetXaxis()->SetRangeUser(0.23,29.9);
    histo2DRatioToCombFit->Draw();

        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys7TeVEta, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioCombCombFitSys7TeVEta->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat7TeVEta, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioCombCombFitStat7TeVEta->Draw("p,same,e1");

        DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);


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
    histo2DRatioToCombFit->Draw("copy");

    for (Int_t i = 0; i < numbersofmeas; i++){
        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);
    }
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3)
            graphRatioCombFitSys[i]->Draw("E2same");
    }
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3)
        graphRatioCombFitStat[i]->Draw("p,same,e");
    }

        DrawGammaLines(0.23, 29.9 , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 29.9 , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 29.9 , 0.9, 0.9,0.5, kGray, 7);

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
        TLatex *textSingleMeasRatioPi0[10];
        for (Int_t i = 0; i < numbersofmeas; i++){
            textSingleMeasRatioPi0[i]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[i+1],nameMeasGlobal[i].Data());
            SetStyleTLatex( textSingleMeasRatioPi0[i], 0.85*textSizeLabelsPixel,4);
            textSingleMeasRatioPi0[i]->SetTextFont(43);
            if(i!=3)
            textSingleMeasRatioPi0[i]->Draw();
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

        TMarker* markerPi0OnlyRatio[10];
        for (Int_t i = 0; i < numbersofmeas; i++){
            markerPi0OnlyRatio[i]               = CreateMarkerFromGraph(graphRatioCombFitSys[i],columnsLegendOnlyPi0Ratio[1]-0.1 ,rowsLegendOnlyPi0Ratio[i+1],1);
            if(i!=3)
            markerPi0OnlyRatio[i]->DrawMarker(columnsLegendOnlyPi0RatioAbs[1]-0.1 ,rowsLegendOnlyPi0RatioAbs[i+1]);
        }

        TBox* boxPi0OnlyRatio[10];
        for (Int_t i = 0; i < numbersofmeas; i++){
            boxPi0OnlyRatio[i]                  = CreateBoxFromGraph(graphRatioCombFitSys[i], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox-0.1 , rowsLegendOnlyPi0RatioAbs[i+1]- heightBox,
                                                        columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox-0.1, rowsLegendOnlyPi0RatioAbs[i+1]+ heightBox);
            if(i!=3)
            boxPi0OnlyRatio[i]->Draw("l");
        }

    canvasRatioToCombFit->SaveAs(Form("%s/RatioOfIndividualMeasToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));
    //    **********************************************************************************************************************
    //    ******************************************* Ratio of Individual meas to Fit ******************************************
    //    **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    TH2F * histo2DRatioToCombFitEviDaniel;
    histo2DRatioToCombFitEviDaniel                       = new TH2F("histo2DRatioToCombFitEviDaniel","histo2DRatioToCombFitEviDaniel",1000,0.51,29.9,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DRatioToCombFitEviDaniel, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DRatioToCombFitEviDaniel->GetXaxis()->SetMoreLogLabels();
    histo2DRatioToCombFitEviDaniel->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRatioToCombFitEviDaniel->GetYaxis()->SetRangeUser(0.61,1.65);
    histo2DRatioToCombFitEviDaniel->Draw("copy");
    Color_t colorEMCDaniel = kGreen-8;
    DrawGammaSetMarkerTGraphAsym(graphRatioCombFitSys[5], markerStyleDet[5] ,markerSizeDet[5]*0.5, colorEMCDaniel, colorEMCDaniel, widthLinesBoxes, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphRatioCombFitStat[5], markerStyleDet[5] ,markerSizeDet[5]*0.5, colorEMCDaniel, colorEMCDaniel);
            graphRatioCombFitSys[2]->Draw("E2same");
            graphRatioCombFitSys[5]->Draw("E2same");

        graphRatioCombFitStat[2]->Draw("p,same,e");
        graphRatioCombFitStat[5]->Draw("p,same,e");

        DrawGammaLines(0.51, 29.9 , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.51, 29.9 , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.51, 29.9 , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitPi0->Draw();

      
            textSingleMeasRatioPi0[2]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[0+1],"EMCal LHC11");
              textSingleMeasRatioPi0[5]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1+1],"EMCal LHC10");
            SetStyleTLatex( textSingleMeasRatioPi0[2], 0.85*textSizeLabelsPixel,4);
            textSingleMeasRatioPi0[2]->SetTextFont(43);
            SetStyleTLatex( textSingleMeasRatioPi0[5], 0.85*textSizeLabelsPixel,4);
            textSingleMeasRatioPi0[5]->SetTextFont(43);
            textSingleMeasRatioPi0[2]->Draw();
            textSingleMeasRatioPi0[5]->Draw();

        //****************** second Column *************************************************;
        textStatOnlyRatioPi0->Draw();
        textSysOnlyRatioPi0->Draw();
        Double_t rowsLegendOnlyPi0RatioAbsEMC[6]   = {0.91,1.53,1.48,2.0,1.95,1.9};
        Double_t columnsLegendOnlyPi0RatioAbsEMC[3]= {0.15,1.8, 2.2};
        Double_t lengthBoxEMC                      = 0.3/2;
        Double_t heightBoxEMC                      = 0.04/2;
            markerPi0OnlyRatio[7]               = CreateMarkerFromGraph(graphRatioCombFitSys[2],columnsLegendOnlyPi0Ratio[1]-0.1 ,rowsLegendOnlyPi0Ratio[0+1],1);
            markerPi0OnlyRatio[7]->DrawMarker(columnsLegendOnlyPi0RatioAbsEMC[1]-0.1 ,rowsLegendOnlyPi0RatioAbsEMC[0+1]);
            markerPi0OnlyRatio[8]               = CreateMarkerFromGraph(graphRatioCombFitSys[5],columnsLegendOnlyPi0Ratio[1]-0.1 ,rowsLegendOnlyPi0Ratio[1+1],1);
            markerPi0OnlyRatio[8]->DrawMarker(columnsLegendOnlyPi0RatioAbsEMC[1]-0.1 ,rowsLegendOnlyPi0RatioAbsEMC[1+1]-0.01);

            boxPi0OnlyRatio[7]                  = CreateBoxFromGraph(graphRatioCombFitSys[2], columnsLegendOnlyPi0RatioAbsEMC[2]-0.5*lengthBoxEMC-0.1 , rowsLegendOnlyPi0RatioAbsEMC[0+1]- heightBoxEMC,
              columnsLegendOnlyPi0RatioAbsEMC[2]+ 3*lengthBoxEMC-0.1, rowsLegendOnlyPi0RatioAbsEMC[0+1]+ heightBoxEMC);
            boxPi0OnlyRatio[7]->Draw("l");
            boxPi0OnlyRatio[8]                  = CreateBoxFromGraph(graphRatioCombFitSys[5], columnsLegendOnlyPi0RatioAbsEMC[2]-0.5*lengthBoxEMC-0.1 , rowsLegendOnlyPi0RatioAbsEMC[1+1]- heightBoxEMC,
              columnsLegendOnlyPi0RatioAbsEMC[2]+ 3*lengthBoxEMC-0.1, rowsLegendOnlyPi0RatioAbsEMC[1+1]+ heightBoxEMC);
            boxPi0OnlyRatio[8]->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/RatioOfEMCMeasToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));
    
    //    **********************************************************************************************************************
    //    ******************* ETA ETA  ************* Ratio of Individual meas to Fit ***************** ETA  ETA ****************
    //    **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DRatioToCombFit->GetXaxis()->SetRangeUser(0.23,39.9);
    histo2DRatioToCombFit->Draw("copy");

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3){
        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitSysEta[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioCombFitStatEta[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);
        }
    }
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3)
            graphRatioCombFitSysEta[i]->Draw("E2same");
    }
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3)
        graphRatioCombFitStatEta[i]->Draw("p,same,e");
    }

        DrawGammaLines(0.23, 39.9 , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 39.9 , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 39.9 , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitEta->Draw();

        //****************************** Definition of the Legend ******************************************
        //**************** Row def ************************

        //****************** first Column **************************************************
        for (Int_t i = 0; i < numbersofmeas; i++){
            textSingleMeasRatioPi0[i]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[i+1],nameMeasGlobal[i].Data());
            SetStyleTLatex( textSingleMeasRatioPi0[i], 0.85*textSizeLabelsPixel,4);
            textSingleMeasRatioPi0[i]->SetTextFont(43);
            if(i!=3)
            textSingleMeasRatioPi0[i]->Draw();
        }

        //****************** second Column *************************************************
        textStatOnlyRatioPi0->Draw();
        textSysOnlyRatioPi0->Draw();

        for (Int_t i = 0; i < numbersofmeas; i++){
            markerPi0OnlyRatio[i]               = CreateMarkerFromGraph(graphRatioCombFitSys[i],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[i+1],1);
            if(i!=3)
            markerPi0OnlyRatio[i]->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[i+1]);
        }

        for (Int_t i = 0; i < numbersofmeas; i++){
            boxPi0OnlyRatio[i]                  = CreateBoxFromGraph(graphRatioCombFitSys[i], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[i+1]- heightBox,
                                                        columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox+0.1, rowsLegendOnlyPi0RatioAbs[i+1]+ heightBox);
            if(i!=3)
            boxPi0OnlyRatio[i]->Draw("l");
        }

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));


    TString nameMeasGlobal22[11]                = {"PCM", "PHOS", "EMCAL", "PCMPHOS", "PCMEMCAL", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCMOtherDataset"};
    
    TString nameOutputCommonFile                = Form("%s/CombinedResultsPaperPP7TeV_%s.root", outputDir.Data(), dateForOutput.Data());    
    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    // PI0 MESON
    fCombResults.mkdir("Pi07TeV");
    TDirectoryFile* directoryPi02               = (TDirectoryFile*)fCombResults.Get("Pi07TeV");
    fCombResults.cd("Pi07TeV");

    graphCombPi0InvCrossSectionStatPCMEMCPHOS   ->RemovePoint(0);
    graphCombPi0InvCrossSectionSysPCMEMCPHOS    ->RemovePoint(0);
    graphCombPi0InvCrossSectionStatPCMEMCPHOS   ->Write("graphInvCrossSectionPi0CombStat");
    graphCombPi0InvCrossSectionSysPCMEMCPHOS    ->Write("graphInvCrossSectionPi0CombSys");
    graphCombPi0InvCrossSectionStatPCMEMCPHOS   ->Write("graphInvCrossSectionPi0Comb7TeVAStatErr");
    graphCombPi0InvCrossSectionSysPCMEMCPHOS    ->Write("graphInvCrossSectionPi0Comb7TeVASysErr");
    

    TGraphAsymmErrors* graphCombPi0InvYieldStat = ScaleGraph( graphCombPi0InvCrossSectionStatPCMEMCPHOS,1/xSection7TeV/recalcBarn);
    TGraphAsymmErrors* graphCombPi0InvYieldSys = ScaleGraph( graphCombPi0InvCrossSectionSysPCMEMCPHOS,1/xSection7TeV/recalcBarn);;
    graphCombPi0InvYieldStat   ->Write("graphCombPi0InvYieldStat");
    graphCombPi0InvYieldSys    ->Write("graphCombPi0InvYieldSys");
    
    while(graphPi0InvCrossSectionStat[2]->GetX()[0] < 1.2) graphPi0InvCrossSectionStat[2]->RemovePoint(0);
    while(graphPi0InvCrossSectionStat[4]->GetX()[0] < 0.8) graphPi0InvCrossSectionStat[4]->RemovePoint(0);
    for (Int_t i = 0; i < numbersofmeas; i++){
        while(graphPi0InvCrossSectionStat[i]->GetY()[0] == 0) graphPi0InvCrossSectionStat[i]->RemovePoint(0);
        while(graphPi0InvCrossSectionSys[i]->GetY()[0] == 0) graphPi0InvCrossSectionStat[i]->RemovePoint(0);
        graphPi0InvCrossSectionStat[i]          ->Write(Form("graphInvCrossSectionPi0%sStat",nameMeasGlobal[i].Data()));
        graphPi0InvCrossSectionSys[i]           ->Write(Form("graphInvCrossSectionPi0%sSys",nameMeasGlobal[i].Data()));
        graphPi0InvCrossSectionStat[i]          ->Write(Form("graphInvCrossSectionPi0%s7TeVStatErr",nameMeasGlobal22[i].Data()));
        graphPi0InvCrossSectionSys[i]           ->Write(Form("graphInvCrossSectionPi0%s7TeVSysErr",nameMeasGlobal22[i].Data()));
    }
    fitInvCrossSectionPi07TeV                   ->Write("TsallisFitPi0");
    fitTCMInvCrossSectionPi07TeV                ->Write("TwoComponentModelFitPi0");

    // ETA MESON
    fCombResults.mkdir("Eta7TeV");
    TDirectoryFile* directoryEta2               = (TDirectoryFile*)fCombResults.Get("Eta7TeV");
    fCombResults.cd("Eta7TeV");

    graphCombEtaInvCrossSectionStatPCMEMCPHOS   ->RemovePoint(0);
    graphCombEtaInvCrossSectionSysPCMEMCPHOS    ->RemovePoint(0);
    graphCombEtaInvCrossSectionStatPCMEMCPHOS   ->Write("graphInvCrossSectionEtaCombStat");
    graphCombEtaInvCrossSectionSysPCMEMCPHOS    ->Write("graphInvCrossSectionEtaCombSys");
    graphCombEtaInvCrossSectionStatPCMEMCPHOS   ->Write("graphInvCrossSectionEtaComb7TeVAStatErr");
    graphCombEtaInvCrossSectionSysPCMEMCPHOS    ->Write("graphInvCrossSectionEtaComb7TeVASysErr");
    TGraphAsymmErrors* graphCombEtaInvYieldStat = ScaleGraph( graphCombEtaInvCrossSectionStatPCMEMCPHOS,1/xSection7TeV/recalcBarn);
    TGraphAsymmErrors* graphCombEtaInvYieldSys = ScaleGraph( graphCombEtaInvCrossSectionSysPCMEMCPHOS,1/xSection7TeV/recalcBarn);

    graphCombEtaInvYieldStat   ->Write("graphCombEtaInvYieldStat");
    graphCombEtaInvYieldSys    ->Write("graphCombEtaInvYieldSys");
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3){
            while(graphEtaInvCrossSectionStat[i]->GetY()[0] < 1e-50) graphEtaInvCrossSectionStat[i]->RemovePoint(0);
            while(graphEtaInvCrossSectionSys[i]->GetY()[0] < 1e-50) graphEtaInvCrossSectionSys[i]->RemovePoint(0);
            graphEtaInvCrossSectionStat[i]      ->Write(Form("graphInvCrossSectionEta%sStat",nameMeasGlobal[i].Data()));
            graphEtaInvCrossSectionSys[i]       ->Write(Form("graphInvCrossSectionEta%sSys",nameMeasGlobal[i].Data()));
            graphEtaInvCrossSectionStat[i]      ->Write(Form("graphInvCrossSectionEta%s7TeVStatErr",nameMeasGlobal22[i].Data()));
            graphEtaInvCrossSectionSys[i]       ->Write(Form("graphInvCrossSectionEta%s7TeVSysErr",nameMeasGlobal22[i].Data()));
        }
    }
    
    graphCombEtaToPi0StatPCMEMCPHOS   ->RemovePoint(0);
    graphCombEtaToPi0SysPCMEMCPHOS    ->RemovePoint(0);
    graphCombEtaToPi0StatPCMEMCPHOS             ->Write("graphEtaToPi0CombStat");
    graphCombEtaToPi0SysPCMEMCPHOS              ->Write("graphEtaToPi0CombSys");
    graphCombEtaToPi0StatPCMEMCPHOS             ->Write("graphRatioEtaToPi0Comb7TeVStatErr");
    graphCombEtaToPi0SysPCMEMCPHOS              ->Write("graphRatioEtaToPi0Comb7TeVSysErr");
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3){
            while(graphEtaToPi0Stat[i]->GetY()[0] < 1e-50)graphEtaToPi0Stat[i]->RemovePoint(0);
            while(graphEtaToPi0Sys[i]->GetY()[0] < 1e-50) graphEtaToPi0Sys[i]->RemovePoint(0);
            graphEtaToPi0Stat[i]                ->Write(Form("graphEtaToPi0%sStat",nameMeasGlobal[i].Data()));
            graphEtaToPi0Sys[i]                 ->Write(Form("graphEtaToPi0%sSys",nameMeasGlobal[i].Data()));
            graphEtaToPi0Stat[i]                ->Write(Form("graphRatioEtaToPi0%s7TeVStatErr",nameMeasGlobal22[i].Data()));
            graphEtaToPi0Sys[i]                 ->Write(Form("graphRatioEtaToPi0%s7TeVSysErr",nameMeasGlobal22[i].Data()));
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

        TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, 0.23, 50. ,1000., -30, 40);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
        histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,24.5);
        histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
        histo2DAllPi0FWHM->DrawCopy();

        for (Int_t i = 0; i < numbersofmeas; i++){
            if(histoPi0FWHMMeV[i] && histoPi0TrueFWHMMeV[i]&&i!=3){
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

        TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, 0.23, 50., 1000., 125.1, 155.9);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                                  textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
        histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
        histo2DAllPi0Mass->DrawCopy();

        for (Int_t i = 0; i < numbersofmeas; i++){
            if(histoPi0Mass[i] && histoPi0TrueMass[i]&&i!=3){
                DrawGammaSetMarker(histoPi0Mass[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
                histoPi0Mass[i]->Draw("p,same,e");
                DrawGammaSetMarker(histoPi0TrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
                histoPi0TrueMass[i]->Draw("p,same,e");
            }
        }

        DrawGammaLines(0.23, 50. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

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
        for (Int_t i = 0; i < numbersofmeas; i++){
            if(histoPi0Mass[i] && histoPi0TrueMass[i] && histoPi0FWHMMeV[i] && histoPi0TrueFWHMMeV[i]&&i!=3){
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
        for (Int_t i = 0; i < numbersofmeas; i++){
            if(histoPi0Mass[i] && histoPi0TrueMass[i]&&i!=3){
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

        TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, 0.23, 50. ,1000., -1, 55.5);
        SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
//         histo2DAllEtaFWHM->GetYaxis()->SetRangeUser(-1.,45.5);
        histo2DAllEtaFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllEtaFWHM->GetYaxis()->SetNdivisions(505);
        histo2DAllEtaFWHM->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaFWHM->GetXaxis()->SetTickLength(0.05);
        histo2DAllEtaFWHM->GetYaxis()->SetTickLength(0.026);
        histo2DAllEtaFWHM->DrawCopy();

        for (Int_t i = 0; i <numbersofmeas; i++){
            if(histoEtaFWHMMeV[i] && histoEtaTrueFWHMMeV[i]&&i!=3){
                DrawGammaSetMarker(histoEtaFWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
                histoEtaFWHMMeV[i]->Draw("p,same,e");
                DrawGammaSetMarker(histoEtaTrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
                histoEtaTrueFWHMMeV[i]->Draw("p,same,e");
            }
        }

        drawLatexAdd("b)",0.13,0.06,textSizeLabelsPixel,kTRUE);
        drawLatexAdd("ALICE performance",0.13,0.87,textSizeLabelsPixel,kTRUE);
        drawLatexAdd(collisionSystem7TeV.Data(),0.13,0.78,textSizeLabelsPixel,kTRUE);
        drawLatexAdd("#eta #rightarrow #gamma#gamma",0.13,0.69,textSizeLabelsPixel,kTRUE);

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

        TH2F * histo2DAllEtaMass            = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",20, 0.23, 50., 1000., 515.1, 589.9);
        SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMassEta, textsizeLabelsMassEta, 0.85*textsizeLabelsMassEta,
                                  textsizeLabelsMassEta, 0.9, 0.28/(textsizeFacMassEta*margin), 512, 505);
        histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllEtaMass->GetYaxis()->SetNdivisions(505);
        histo2DAllEtaMass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaMass->GetXaxis()->SetTickLength(0.05);
        histo2DAllEtaMass->GetXaxis()->SetLabelOffset(-0.015);
        histo2DAllEtaMass->DrawCopy();

        for (Int_t i = 0; i < numbersofmeas; i++){
            if(histoEtaMass[i] && histoEtaTrueMass[i]&&i!=3){
                DrawGammaSetMarker(histoEtaMass[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
                histoEtaMass[i]->Draw("p,same,e");
                DrawGammaSetMarker(histoEtaTrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
                histoEtaTrueMass[i]->Draw("p,same,e");
            }
        }

        DrawGammaLines(0.23, 50. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.1, kGray);

        drawLatexAdd("b)",0.13,0.22,textSizeLabelsPixel,kTRUE);


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
        for (Int_t i = 0; i < numbersofmeas; i++){
            if(histoEtaMass[i] && histoEtaTrueMass[i] && histoEtaFWHMMeV[i] && histoEtaTrueFWHMMeV[i]&&i!=3){
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
        for (Int_t i = 0; i < numbersofmeas; i++){
            if(histoEtaMass[i] && histoEtaTrueMass[i]&&i!=3){
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
        histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.23,  39, 1000, 8e-6, 2 );
        SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);//(#times #epsilon_{pur})
        histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
        histo2DAccEff->GetXaxis()->SetLabelOffset(-0.01);
        histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAccEff->DrawCopy();

        for (Int_t i = 0; i < numbersofmeas; i++){
            if(histoPi0AccTimesEff[i]&&i!=3){
                DrawGammaSetMarker(histoPi0AccTimesEff[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
                histoPi0AccTimesEff[i]->Draw("p,same,e");
            }
        }

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        for (Int_t i = 0; i < numbersofmeas; i++){
            if(histoPi0AccTimesEff[i]&&i!=3){
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
    histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
    histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.0,1.05);
    histo2DEtatoPi0combo->Draw();
        // plotting systematics graphs
        for (Int_t i = 0; i < numbersofmeas; i++){
            if(graphEtaToPi0Sys[i]&&i!=3){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Sys[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
                graphEtaToPi0Sys[i]->Draw("E2same");
            }
        }
         DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatPCMEMCPHOS, markerStyleDet[1], markerSizeDet[1]*0.75, kPink , kPink, widthLinesBoxes, kTRUE);
//                 graphCombEtaToPi0StatPCMEMCPHOS->Draw("E2same");
        // plotting statistics graphs
        for (Int_t i = 0; i < numbersofmeas; i++){
            if(graphEtaToPi0Stat[i]&&i!=3){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Stat[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i]);
                graphEtaToPi0Stat[i]->Draw("p,same,e");
            }
        }
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0SysPCMEMCPHOS, markerStyleDet[1], markerSizeDet[1]*0.75, kPink , kPink);
//                 graphCombEtaToPi0SysPCMEMCPHOS->Draw("p,same,e");

        TLegend* legendEtaToPi0 = GetAndSetLegend2(0.67, 0.15, 0.9, 0.15+(textsizeLabelsEtaToPi0*4*0.9), textSizeLabelsPixel);
        for (Int_t i = 0; i < numbersofmeas; i++){
            if(graphEtaToPi0Sys[i]&&i!=3){
                legendEtaToPi0->AddEntry(graphEtaToPi0Sys[i],nameMeasGlobal[i],"pf");
            }
        }
        legendEtaToPi0->Draw();

        drawLatexAdd(collisionSystem7TeV.Data(),0.13, 0.92,0.85*textsizeLabelsEtaToPi0,kFALSE);
        drawLatexAdd("ALICE",0.13, 0.92-(1*textsizeLabelsEtaToPi0*0.85),0.85*textsizeLabelsEtaToPi0,kFALSE);
        drawLatexAdd("#eta/#pi^{0}",0.13, 0.92-(2*textsizeLabelsEtaToPi0*0.9),textsizeLabelsEtaToPi0,kFALSE);

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

















