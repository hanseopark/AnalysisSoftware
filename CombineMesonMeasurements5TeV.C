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


void CombineMesonMeasurements5TeV(      TString fileNamePCM     = "",
                                        TString fileNamePHOS    = "",
                                        TString fileNameEMCal   = "/home/nschmidt/AnalysisSoftware/pdf/5TeV/2017_08_31/FinalResultsTriggersPatched_EMC_27/data_EMCAL-EMCALResultsFullCorrection_PP.root",
                                        TString fileNamePCMPHOS = "",
                                        TString fileNamePCMEMCal= "/home/nschmidt/AnalysisSoftware/pdf/5TeV/2017_08_31/FinalResultsTriggersPatched_PCMEMC_27/data_PCM-EMCALResultsFullCorrection_PP.root",
                                        TString fileNameInterpolation = "/home/nschmidt/AnalysisSoftware/CombinationInput5TeV/Interpolation.root",
                                        TString suffix          = "pdf",
                                        Int_t numbersofmeas     = 5
                                    ){

    TString date                                = ReturnDateString();

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    //___________________________________ Declaration of files _____________________________________________
    TString collisionSystem5TeV               = "pp, #sqrt{#it{s}} = 5 TeV";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements5TeV",suffix.Data(),dateForOutput.Data());
    cout << outputDir.Data() << endl;
    cout << fileNamePCM.Data() << endl;

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
    if(fileNameInterpolation.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputInterpolation.root", fileNameInterpolation.Data(), outputDir.Data()));

    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();

    Width_t widthLinesBoxes                     = 1.4;

    // Definition of colors, styles and markers sizes
    Color_t colorComb                           = kBlue+2;
    Style_t markerStyleComb                     = 20;
    Size_t  markerSizeComb                      = 2;

    TString nameMeasGlobal[11]                  = {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCMOtherDataset"};
    Color_t colorDet[11];
    Color_t colorDetMC[11];
    Style_t markerStyleDet[11];
    Style_t markerStyleDetMC[11];
    Size_t  markerSizeDet[11];
    Size_t  markerSizeDetMC[11];

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
        inputFile[2]                            = new TFile(fileNameEMCal.Data());
        inputFile[3]                            = new TFile(fileNamePCMPHOS.Data());
        inputFile[4]                            = new TFile(fileNamePCMEMCal.Data());
    TFile* inputFileInterpolation               = new TFile(fileNameInterpolation.Data());

    TDirectory* directoryPi0[10];
    TDirectory* directoryEta[10];
    for(Int_t i=0;i<numbersofmeas;i++){
      if(!inputFile[i]->IsZombie()){
        cout << "loading directories for " <<  nameMeasGlobal[i] << endl;
        directoryPi0[i]                     = (TDirectory*)inputFile[i]->Get("Pi05TeV");
        directoryEta[i]                     = (TDirectory*)inputFile[i]->Get("Eta5TeV");
      }
    }
    cout << __LINE__<<endl;
    TH1D* histoNumberOfEvents[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0Mass[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvCrossSectionSys[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0FWHMMeV[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0TrueMass[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0TrueFWHMMeV[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaMass[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaFWHMMeV[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaTrueMass[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaTrueFWHMMeV[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0Acc[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0TrueEffPt[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0AccTimesEff[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaAcc[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaTrueEffPt[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaAccTimesEff[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvCrossSection[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaInvCrossSection[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0RawYields[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaRawYields[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphPi0RawYields[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaRawYields[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphPi0InvCrossSectionSys[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphPi0InvCrossSectionStat[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaInvCrossSectionStat[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaInvCrossSectionSys[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaToPi0PCMStat;
    TGraphAsymmErrors* graphEtaToPi0PCMSys;
    TH1D* histoPi0InvMassSigPlusBG[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvMassSig[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvMassSigRemBGSub[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvMassBG[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvMassRemBG[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoPi0InvMassBGTot[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaInvMassSigPlusBG[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaInvMassSig[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaInvMassSigRemBGSub[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaInvMassBG[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaInvMassRemBG[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1D* histoEtaInvMassBGTot[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TF1* fitPi0InvMassSig[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TF1* fitPi0InvMassBG[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TF1* fitEtaInvMassSig[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TF1* fitEtaInvMassBG[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    Bool_t haveAllPi0InvMass[10]                = {kFALSE, kFALSE, kFALSE,kFALSE,kFALSE};
    Bool_t haveAllEtaInvMass[10]                = {kFALSE, kFALSE, kFALSE,kFALSE,kFALSE};
    TH1D* histoEtaToPi0Stat[10]= {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaToPi0Stat[10]    = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaToPi0Sys[10]     = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    Double_t rapidityMeas[10]                   = {1.6, 1,1, 1.6,1,1,1};
    Double_t availableMeas[10]                  = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
    Int_t nMeasSet                              = 0;

    for (Int_t i = 0; i < 5; i++){
      if(inputFile[i]->IsZombie())
        continue;
      else
        cout << "loading  histos for " << nameMeasGlobal[i] << endl;
        
    //______________________________ Neutral pion cross section and raw yield
      histoPi0RawYields[i]          = (TH1D*)directoryPi0[i]->Get("RAWYieldPerEventsPi0_INT7");
      graphPi0RawYields[i]      = new TGraphAsymmErrors(histoPi0RawYields[i]);
      graphPi0InvCrossSectionSys[i]       = (TGraphAsymmErrors*)directoryPi0[i]->Get("InvCrossSectionPi0Sys");
      histoPi0InvCrossSection[i]          = (TH1D*)directoryPi0[i]->Get("InvCrossSectionPi0");
      graphPi0InvCrossSectionStat[i]      = new TGraphAsymmErrors(histoPi0InvCrossSection[i]);
      if(!graphPi0InvCrossSectionSys[i] || !histoPi0InvCrossSection[i]){
        cout << "missing invariant cross section histograms of pi0... returning!" << endl;
        return;
      } else {
        availableMeas[i]                  = kTRUE;
        nMeasSet++;
        cout << "loaded " << nameMeasGlobal[i] << " pi0 invariant cross section" << endl;
      }
      
      //______________________________ Neutral pion invariant mass peak pos and FWHM
      histoPi0Mass[i]                     = (TH1D*)directoryPi0[i]->Get("Pi0_Mass_data_INT7");
      histoPi0FWHMMeV[i]                  = (TH1D*)directoryPi0[i]->Get("Pi0_Width_data_INT7");
      histoPi0TrueMass[i]                 = (TH1D*)directoryPi0[i]->Get("Pi0_Mass_MC_INT7");
      histoPi0TrueFWHMMeV[i]              = (TH1D*)directoryPi0[i]->Get("Pi0_Width_MC_INT7");
      if(!histoPi0Mass[i] || !histoPi0FWHMMeV[i] || !histoPi0FWHMMeV[i] || !histoPi0FWHMMeV[i]){
        cout << "missing mass or width histograms... returning!" << endl;
        return;
      } else {
        // scaling peak pos and FWHM to go from GeV/c2 to MeV/c2
        histoPi0Mass[i]                   ->Scale(1000);
        histoPi0FWHMMeV[i]                ->Scale(1000);
        histoPi0TrueMass[i]               ->Scale(1000);
        histoPi0TrueFWHMMeV[i]            ->Scale(1000);
        cout << "loaded " << nameMeasGlobal[i] << " mass and width" << endl;
      }
      
      //______________________________ Neutral pion acceptance and efficiency and calculate acc*eff*y*2pi
      histoPi0Acc[i]                      = (TH1D*)directoryPi0[i]->Get("AcceptancePi0_INT7");
      histoPi0TrueEffPt[i]                = (TH1D*)directoryPi0[i]->Get("EfficiencyPi0_INT7");
      if(!histoPi0Acc[i] || !histoPi0TrueEffPt[i]){
        cout << "missing acceptance or efficiency histograms... returning!" << endl;
        return;
      } else {
        // calculating acceptance times efficiency
        histoPi0AccTimesEff[i]            = (TH1D*)histoPi0TrueEffPt[i]->Clone(Form("histoPi0AccTimesEff%s",nameMeasGlobal[i].Data()));
        histoPi0AccTimesEff[i]->Multiply(histoPi0Acc[i]);
        histoPi0AccTimesEff[i]->Scale(2*TMath::Pi()*rapidityMeas[i]);
        cout << "loaded " << nameMeasGlobal[i] << "acceptance and efficiency" << endl;
      }
      
      //______________________________ Neutral pion invariant mass example bin
      histoPi0InvMassSig[i]               = (TH1D*)directoryPi0[i]->Get("Pi0_InvMassSig_Example_INT7");
      histoPi0InvMassSigPlusBG[i]         = (TH1D*)directoryPi0[i]->Get("Pi0_InvMassSigPlusBG_Example_INT7");
      histoPi0InvMassBG[i]                = (TH1D*)directoryPi0[i]->Get("Pi0_InvMassBG_Example_INT7");
      fitPi0InvMassSig[i]                 = (TF1*)directoryPi0[i]->Get("Pi0_InvMassSigFit_Example_INT7");
      if (histoPi0InvMassSig[i] && histoPi0InvMassSigPlusBG[i] && histoPi0InvMassBG[i] && fitPi0InvMassSig[i]){
        haveAllPi0InvMass[i]              = kTRUE;
        histoPi0InvMassBGTot[i]           = (TH1D*)histoPi0InvMassBG[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameMeasGlobal[i].Data()));
        histoPi0InvMassSigRemBGSub[i]     = (TH1D*)histoPi0InvMassSig[i]->Clone(Form("Pi0_InvMassSigRemBGSub_Example_%s",nameMeasGlobal[i].Data()));
        cout << "loaded " << nameMeasGlobal[i] << " example mass bin" << endl;
      } else {
        cout << "missing pi0 invariant mass example bin histograms... they will not be plotted!" << endl;
      }

      //______________________________ Eta meson invariant cross section
        histoEtaRawYields[i]          = (TH1D*)directoryEta[i]->Get("RAWYieldPerEventsEta_INT7");
        graphEtaRawYields[i]      = new TGraphAsymmErrors(histoEtaRawYields[i]);
        histoEtaInvCrossSection[i]        = (TH1D*)directoryEta[i]->Get("InvCrossSectionEta");
        graphEtaInvCrossSectionStat[i]    = (TGraphAsymmErrors*)directoryEta[i]->Get("graphInvCrossSectionEta");
        graphEtaInvCrossSectionSys[i]     = (TGraphAsymmErrors*)directoryEta[i]->Get("InvCrossSectionEtaSys");
        graphEtaToPi0PCMStat              = (TGraphAsymmErrors*)directoryEta[i]->Get("graphEtaToPi0StatError");
        graphEtaToPi0PCMSys               = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0SystError");
        if(!graphEtaInvCrossSectionSys[i] || !histoEtaInvCrossSection[i]){
          cout << "missing invariant cross section histograms of eta... returning!" << endl;
          return;
        } else {
          cout << "loaded " << nameMeasGlobal[i] << " eta invariant cross section" << endl;
        }
        
        //______________________________ Eta meson invariant mass peak pos and FWHM
        histoEtaMass[i]                     = (TH1D*)directoryEta[i]->Get("Eta_Mass_data_INT7");
        histoEtaFWHMMeV[i]                  = (TH1D*)directoryEta[i]->Get("Eta_Width_data_INT7");
        histoEtaTrueMass[i]                 = (TH1D*)directoryEta[i]->Get("Eta_Mass_MC_INT7");
        histoEtaTrueFWHMMeV[i]              = (TH1D*)directoryEta[i]->Get("Eta_Width_MC_INT7");
        if(!histoEtaMass[i] || !histoEtaFWHMMeV[i] || !histoEtaFWHMMeV[i] || !histoEtaFWHMMeV[i]){
          cout << "missing mass or width histograms... returning!" << endl;
          return;
        } else {
          // scaling peak pos and FWHM to go from GeV/c2 to MeV/c2
          histoEtaMass[i]                   ->Scale(1000);
          histoEtaFWHMMeV[i]                ->Scale(1000);
          histoEtaTrueMass[i]               ->Scale(1000);
          histoEtaTrueFWHMMeV[i]            ->Scale(1000);
          cout << "loaded " << nameMeasGlobal[i] << " mass and width" << endl;
        }
        
        //______________________________ Neutral pion acceptance and efficiency and calculate acc*eff*y*2pi
        histoEtaAcc[i]                      = (TH1D*)directoryEta[i]->Get("AcceptanceEta_INT7");
        histoEtaTrueEffPt[i]                = (TH1D*)directoryEta[i]->Get("EfficiencyEta_INT7");
        if(!histoEtaAcc[i] || !histoEtaTrueEffPt[i]){
          cout << "missing acceptance or efficiency histograms... returning!" << endl;
          return;
        } else {
          // calculating acceptance times efficiency
          histoEtaAccTimesEff[i]            = (TH1D*)histoEtaTrueEffPt[i]->Clone(Form("histoEtaAccTimesEff%s",nameMeasGlobal[i].Data()));
          histoEtaAccTimesEff[i]->Multiply(histoEtaAcc[i]);
          histoEtaAccTimesEff[i]->Scale(2*TMath::Pi()*rapidityMeas[i]);
          cout << "loaded " << nameMeasGlobal[i] << "acceptance and efficiency" << endl;
        }
        
        //______________________________ Eta meson invariant mass example bin
        histoEtaInvMassSig[i]               = (TH1D*)directoryEta[i]->Get("Eta_InvMassSig_Example_INT7");
        histoEtaInvMassSigPlusBG[i]         = (TH1D*)directoryEta[i]->Get("Eta_InvMassSigPlusBG_Example_INT7");
        histoEtaInvMassBG[i]                = (TH1D*)directoryEta[i]->Get("Eta_InvMassBG_Example_INT7");
        fitEtaInvMassSig[i]                 = (TF1*)directoryEta[i]->Get("Eta_InvMassSigFit_Example_INT7");
        if (histoEtaInvMassSig[i] && histoEtaInvMassSigPlusBG[i] && histoEtaInvMassBG[i] && fitEtaInvMassSig[i]){
          haveAllEtaInvMass[i]            = kTRUE;
          histoEtaInvMassBGTot[i]         = (TH1D*)histoEtaInvMassBG[i]->Clone(Form("Eta_InvMassTotBG_Example_%s",nameMeasGlobal[i].Data()));
          histoEtaInvMassSigRemBGSub[i]   = (TH1D*)histoEtaInvMassSig[i]->Clone(Form("Eta_InvMassSigRemBGSub_Example_%s",nameMeasGlobal[i].Data()));
          cout << "loaded " << nameMeasGlobal[i] << " example mass bin" << endl;
        } else {
          cout << "missing eta invariant mass example bin histograms... they will not be plotted!" << endl;
        }
      
    }
    // loading interpolation inputs
    TGraphAsymmErrors* graphPi0InvCrossSectionSysInterpolation[10];
    TGraphAsymmErrors* graphPi0InvCrossSectionStatInterpolation[10];
    TGraphAsymmErrors* graphEtaInvCrossSectionStatInterpolation[10];
    TGraphAsymmErrors* graphEtaInvCrossSectionSysInterpolation[10];
    graphPi0InvCrossSectionStatInterpolation[4]      = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionStatErrPCM-EMC_Pi0_5.023TeV");
    graphPi0InvCrossSectionSysInterpolation[4]       = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionSystErrPCM-EMC_Pi0_5.023TeV");
    graphPi0InvCrossSectionStatInterpolation[2]      = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionStatErrEMC_Pi0_5.023TeV");
    graphPi0InvCrossSectionSysInterpolation[2]       = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionSystErrEMC_Pi0_5.023TeV");
    graphPi0InvCrossSectionStatInterpolation[9]      = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionStatErrComb_Pi0_5.023TeV");
    graphPi0InvCrossSectionSysInterpolation[9]       = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionSystErrComb_Pi0_5.023TeV");

    graphEtaInvCrossSectionStatInterpolation[4]      = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionStatErrPCM-EMC_Eta_5.023TeV");
    graphEtaInvCrossSectionSysInterpolation[4]       = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionSystErrPCM-EMC_Eta_5.023TeV");
    graphEtaInvCrossSectionStatInterpolation[2]      = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionStatErrEMC_Eta_5.023TeV");
    graphEtaInvCrossSectionSysInterpolation[2]       = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionSystErrEMC_Eta_5.023TeV");
    graphEtaInvCrossSectionStatInterpolation[9]      = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionStatErrComb_Eta_5.023TeV");
    graphEtaInvCrossSectionSysInterpolation[9]       = (TGraphAsymmErrors*)inputFileInterpolation->Get("graphInvXSectionSystErrComb_Eta_5.023TeV");
    
    cout << "Finished loading inputs" << endl;
      
    // calculation of relative statistical and systematic uncertainties
    // FOR pi0:
    TH1D* statErrorCollection[11];
    TGraphAsymmErrors* sysErrorCollection[11];
    TH1D* statErrorRelCollection[11];
    TGraphAsymmErrors* sysErrorRelCollection[11];
    for (Int_t i = 0; i< 11; i++){
      // initialize all histograms and graphs as NULL
        statErrorCollection[i]                  = NULL;
        sysErrorCollection[i]                   = NULL;
        statErrorRelCollection[i]               = NULL;
        sysErrorRelCollection[i]                = NULL;
        // add available measurements to the collection
        if(i<numbersofmeas && availableMeas[i]){
            statErrorCollection[i]              = (TH1D*)histoPi0InvCrossSection[i]->Clone(Form("statErr%sPi0",nameMeasGlobal[i].Data()));
            sysErrorCollection[i]               = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i]->Clone(Form("sysErr%sPi0",nameMeasGlobal[i].Data()));
            cout << "calculating pi0 relative uncertainties for " << nameMeasGlobal[i] << endl;
        }
        // calculate relative errors for the available measurements
        if (statErrorCollection[i]) statErrorRelCollection[i] = CalculateRelErrUpTH1D( statErrorCollection[i], Form("relativeStatErrorPi0_%s", nameMeasGlobal[i].Data()));
        if (sysErrorCollection[i]) sysErrorRelCollection[i]   = CalculateRelErrUpAsymmGraph( sysErrorCollection[i], Form("relativeSysErrorPi0_%s", nameMeasGlobal[i].Data()));  
    }
    // FOR eta:
    TH1D* statErrorCollectionEta[11];
    TGraphAsymmErrors* sysErrorCollectionEta[11];
    TH1D* statErrorRelCollectionEta[11];
    TGraphAsymmErrors* sysErrorRelCollectionEta[11];
    for (Int_t i = 0; i< 11; i++){
      // initialize all histograms and graphs as NULL
        statErrorCollectionEta[i]                  = NULL;
        sysErrorCollectionEta[i]                   = NULL;
        statErrorRelCollectionEta[i]               = NULL;
        sysErrorRelCollectionEta[i]                = NULL;
        // add available measurements to the collection
        if(i<numbersofmeas && availableMeas[i]){
            statErrorCollectionEta[i]              = (TH1D*)histoEtaInvCrossSection[i]->Clone(Form("statErr%sEta",nameMeasGlobal[i].Data()));
            sysErrorCollectionEta[i]               = (TGraphAsymmErrors*)graphEtaInvCrossSectionSys[i]->Clone(Form("sysErr%sEta",nameMeasGlobal[i].Data()));
            cout << "calculating eta relative uncertainties for " << nameMeasGlobal[i] << endl;
        }
        // calculate relative errors for the available measurements
        if (statErrorCollectionEta[i]) statErrorRelCollectionEta[i] = CalculateRelErrUpTH1D( statErrorCollectionEta[i], Form("relativeStatErrorEta_%s", nameMeasGlobal[i].Data()));
        if (sysErrorCollectionEta[i]) sysErrorRelCollectionEta[i]   = CalculateRelErrUpAsymmGraph( sysErrorCollectionEta[i], Form("relativeSysErrorEta_%s", nameMeasGlobal[i].Data()));  
    }


    TString nameMeasGlobal22[11]                = {"PCM", "PHOS", "EMCAL", "PCMPHOS", "PCMEMCAL", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCMOtherDataset"};
    
    TString nameOutputCommonFile                = Form("%s/CombinedResultsPaperPP5TeV_%s.root", outputDir.Data(), dateForOutput.Data());    
    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    // PI0 MESON
    fCombResults.mkdir("Pi05TeV");
    TDirectoryFile* directoryPi02               = (TDirectoryFile*)fCombResults.Get("Pi05TeV");
    fCombResults.cd("Pi05TeV");

    // graphCombPi0InvCrossSectionStat   ->Write("graphInvCrossSectionPi0Comb5TeVAStatErr");
    // graphCombPi0InvCrossSectionSysS    ->Write("graphInvCrossSectionPi0Comb5TeVASysErr");
    cout << __LINE__ << endl;
    
    for (Int_t i = 0; i < numbersofmeas; i++){
      if(graphPi0InvCrossSectionStat[i]&&graphPi0InvCrossSectionSys[i]){
        cout << i << endl;
        while(graphPi0InvCrossSectionStat[i]->GetY()[0] == 0) graphPi0InvCrossSectionStat[i]->RemovePoint(0);
        while(graphPi0InvCrossSectionSys[i]->GetY()[0] == 0) graphPi0InvCrossSectionStat[i]->RemovePoint(0);
        graphPi0InvCrossSectionStat[i]          ->Write(Form("graphInvCrossSectionPi0%s5TeVStatErr",nameMeasGlobal22[i].Data()));
        graphPi0InvCrossSectionSys[i]           ->Write(Form("graphInvCrossSectionPi0%s5TeVSysErr",nameMeasGlobal22[i].Data()));
        cout << i << endl;
      }
    }
cout << __LINE__ << endl;
    // ETA MESON
    fCombResults.mkdir("Eta5TeV");
    TDirectoryFile* directoryEta2               = (TDirectoryFile*)fCombResults.Get("Eta5TeV");
    fCombResults.cd("Eta5TeV");

    // graphCombEtaInvCrossSectionStat   ->Write("graphInvCrossSectionEtaComb5TeVAStatErr");
    // graphCombEtaInvCrossSectionSys    ->Write("graphInvCrossSectionEtaComb5TeVASysErr");
    cout << __LINE__ << endl;
    
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(graphEtaInvCrossSectionStat[i]&&graphEtaInvCrossSectionSys[i]){
          cout << i << endl;
            while(graphEtaInvCrossSectionStat[i]->GetY()[0] < 1e-50) graphEtaInvCrossSectionStat[i]->RemovePoint(0);
            while(graphEtaInvCrossSectionSys[i]->GetY()[0] < 1e-50) graphEtaInvCrossSectionSys[i]->RemovePoint(0);
            graphEtaInvCrossSectionStat[i]      ->Write(Form("graphInvCrossSectionEta%s5TeVStatErr",nameMeasGlobal22[i].Data()));
            graphEtaInvCrossSectionSys[i]       ->Write(Form("graphInvCrossSectionEta%s5TeVSysErr",nameMeasGlobal22[i].Data()));
            cout << i << endl;
        }
      }
      cout << __LINE__ << endl;
  

    // graphCombEtaToPi0Stat             ->Write("graphRatioEtaToPi0Comb5TeVStatErr");
    // graphCombEtaToPi0Sys              ->Write("graphRatioEtaToPi0Comb5TeVSysErr");
    // for (Int_t i = 0; i < numbersofmeas; i++){
    //   if(graphEtaToPi0Stat[i]&&graphEtaToPi0Sys[i]){
    //     cout << i << endl;
            // while(graphEtaToPi0Stat[i]->GetY()[0] < 1e-50)graphEtaToPi0Stat[i]->RemovePoint(0);
            // while(graphEtaToPi0Sys[i]->GetY()[0] < 1e-50) graphEtaToPi0Sys[i]->RemovePoint(0);
    //         graphEtaToPi0Stat[i]                ->Write(Form("graphRatioEtaToPi0%s5TeVStatErr",nameMeasGlobal22[i].Data()));
    //         graphEtaToPi0Sys[i]                 ->Write(Form("graphRatioEtaToPi0%s5TeVSysErr",nameMeasGlobal22[i].Data()));
    //         cout << i << endl;
    //   }
    // }
    
    fCombResults.Close();
    cout << __LINE__ << endl;
    


    Double_t minX                               = 0.61;
    Double_t maxX                               = 29.9;

    // **********************************************************************************************************************
    // ******************************************* Pi0 invariant cross section            ****************************************
    // **********************************************************************************************************************
    Double_t textSizeLabelsPixel             = 50;
    Double_t textSizeLabelsRel      = 50./1200;
    TCanvas* canvasCrossSectionPi0       = new TCanvas("canvasCrossSectionPi0", "", 0,0,1250,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasCrossSectionPi0,  0.14, 0.01, 0.01, 0.08);
    canvasCrossSectionPi0->SetLogy(1);
    canvasCrossSectionPi0->SetLogx(1);
        TH2F * histoDummyCrossSection;
            histoDummyCrossSection                = new TH2F("histoDummyCrossSection", "histoDummyCrossSection",1000, minX,  maxX, 1000, 3, 4e12 );
        SetStyleHistoTH2ForGraphs( histoDummyCrossSection, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2}#it{c}^{3})",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.4);//(#times #epsilon_{pur})
        
        histoDummyCrossSection->GetYaxis()->SetLabelOffset(0.001);
        histoDummyCrossSection->GetXaxis()->SetLabelOffset(-0.01);
        histoDummyCrossSection->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummyCrossSection->DrawCopy();
        Int_t exampleInteger = 0;
        Int_t totalMeasAvail = 0;
        for(Int_t i=0;i<numbersofmeas;i++)
          if(availableMeas[i]){
            exampleInteger = i;
            totalMeasAvail++;
          }
        TH1D * histoBlack               = (TH1D*)histoPi0InvCrossSection[exampleInteger]->Clone("histoBlack") ;
        histoBlack->SetLineColor(kBlack);
        histoBlack->SetMarkerStyle(21) ;
        histoBlack->SetMarkerColor(kBlack) ;
        histoBlack->SetMarkerSize(1.) ;

        TGraphAsymmErrors * graphGrey   = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[exampleInteger]->Clone("graphGrey") ;
        graphGrey   ->SetFillColor(kGray);
        graphGrey   ->SetLineColor(kGray);
        graphGrey   ->SetFillStyle(1001);
        
        TGraphAsymmErrors * graphPi0InvCrossSectionSysClone[10];
        TGraphAsymmErrors * graphPi0InvCrossSectionStatClone[10];
        
        TLegend* legendCrossSectionPi0           = GetAndSetLegend2(0.18, 0.13, 0.43, 0.13+(3.5*textSizeLabelsRel),textSizeLabelsPixel);
        
        Int_t scalingFactor = pow(10,totalMeasAvail-1);
        Int_t legendScalingFactor = totalMeasAvail-1;
        // Int_t legendScalingFactor = 0;
        for(Int_t i=numbersofmeas;i>-1;i--){
          if(availableMeas[i]){
            graphPi0InvCrossSectionSysClone[i]= (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i]->Clone(Form("csSys%d",i)) ;
            graphPi0InvCrossSectionStatClone[i]= (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[i]->Clone(Form("csStat%d",i)) ;
            for (int j=0;j<graphPi0InvCrossSectionSysClone[i]->GetN();j++){
                graphPi0InvCrossSectionSysClone[i]->GetY()[j] *= scalingFactor;
                graphPi0InvCrossSectionSysClone[i]->GetEYhigh()[j] *= scalingFactor;
                graphPi0InvCrossSectionSysClone[i]->GetEYlow()[j] *= scalingFactor;
            }
            // DrawGammaSetMarkerTGraphAsym(graphPi0InvCrossSectionSysInterpolation[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDetMC[i] , colorDetMC[i], widthLinesBoxes, kTRUE);
            // graphPi0InvCrossSectionSysInterpolation[i]     ->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphPi0InvCrossSectionSysClone[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
            graphPi0InvCrossSectionSysClone[i]     ->Draw("E2same");
            legendCrossSectionPi0->AddEntry(graphPi0InvCrossSectionSysClone[i],Form("%s x10^{%d}",nameMeasGlobal[i].Data(),legendScalingFactor),"pf");

            // statistics
            for (int j=0;j<graphPi0InvCrossSectionStatClone[i]->GetN();j++){
                graphPi0InvCrossSectionStatClone[i]->GetY()[j] *= scalingFactor;
                graphPi0InvCrossSectionStatClone[i]->GetEYhigh()[j] *= scalingFactor;
                graphPi0InvCrossSectionStatClone[i]->GetEYlow()[j] *= scalingFactor;
            }
            // DrawGammaSetMarkerTGraph(graphPi0InvCrossSectionStatInterpolation[i],  markerStyleDet[i], markerSizeDet[i]*0.55, colorDetMC[i] , colorDetMC[i]);
            // graphPi0InvCrossSectionStatInterpolation[i]->Draw("p,same,z");
            DrawGammaSetMarkerTGraph(graphPi0InvCrossSectionStatClone[i],  markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            graphPi0InvCrossSectionStatClone[i]->Draw("p,same,z");
            scalingFactor/=10;
            legendScalingFactor--;
          }
        }
        legendCrossSectionPi0->Draw();


        // TLegend* legendPi0Err2 = GetAndSetLegend2(0.72, 0.72, 0.98, 0.72+(2*textSizeLabelsRel),0.85*textSizeLabelsPixel);
        // legendPi0Err2->AddEntry(histoBlack, "stat. Err.","ple");
        // legendPi0Err2->AddEntry(graphGrey,  "syst. Err.","f");
        // legendPi0Err2->Draw();

        drawLatexAdd("ALICE",0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem5TeV,0.93,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        histoDummyCrossSection->Draw("sameaxis");
    canvasCrossSectionPi0->Update();
    canvasCrossSectionPi0->Print(Form("%s/Pi0_InvariantCrossSectionMeas.%s",outputDir.Data(),suffix.Data()));
    
    

    // **********************************************************************************************************************
    // ******************************************* Pi0 invariant cross section            ****************************************
    // **********************************************************************************************************************
        TH2F * histoDummyCrossSectionEta;
            histoDummyCrossSectionEta                = new TH2F("histoDummyCrossSectionEta", "histoDummyCrossSectionEta",1000, minX,  maxX, 1000, 3, 4e10 );
        SetStyleHistoTH2ForGraphs( histoDummyCrossSectionEta, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2}#it{c}^{3})",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.4);//(#times #epsilon_{pur})
        
        histoDummyCrossSectionEta->GetYaxis()->SetLabelOffset(0.001);
        histoDummyCrossSectionEta->GetXaxis()->SetLabelOffset(-0.01);
        histoDummyCrossSectionEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummyCrossSectionEta->DrawCopy();
        exampleInteger = 0;
        totalMeasAvail = 0;
        for(Int_t i=0;i<numbersofmeas;i++)
          if(availableMeas[i]){
            exampleInteger = i;
            totalMeasAvail++;
          }
        
        TGraphAsymmErrors * graphEtaInvCrossSectionSysClone[10];
        TGraphAsymmErrors * graphEtaInvCrossSectionStatClone[10];
        
        TLegend* legendCrossSectionEta           = GetAndSetLegend2(0.18, 0.13, 0.43, 0.13+(3.5*textSizeLabelsRel),textSizeLabelsPixel);
        
        scalingFactor = pow(10,totalMeasAvail-1);
        legendScalingFactor = totalMeasAvail-1;
        // legendScalingFactor = 0;
        for(Int_t i=numbersofmeas;i>-1;i--){
          if(availableMeas[i]){
            graphEtaInvCrossSectionSysClone[i]= (TGraphAsymmErrors*)graphEtaInvCrossSectionSys[i]->Clone(Form("csSys%d",i)) ;
            graphEtaInvCrossSectionStatClone[i]= (TGraphAsymmErrors*)graphEtaInvCrossSectionStat[i]->Clone(Form("csStat%d",i)) ;
            for (int j=0;j<graphEtaInvCrossSectionSysClone[i]->GetN();j++){
                graphEtaInvCrossSectionSysClone[i]->GetY()[j] *= scalingFactor;
                graphEtaInvCrossSectionSysClone[i]->GetEYhigh()[j] *= scalingFactor;
                graphEtaInvCrossSectionSysClone[i]->GetEYlow()[j] *= scalingFactor;
            }
            // DrawGammaSetMarkerTGraphAsym(graphEtaInvCrossSectionSysInterpolation[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDetMC[i] , colorDetMC[i], widthLinesBoxes, kTRUE);
            // graphEtaInvCrossSectionSysInterpolation[i]     ->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphEtaInvCrossSectionSysClone[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
            graphEtaInvCrossSectionSysClone[i]     ->Draw("E2same");
            legendCrossSectionEta->AddEntry(graphEtaInvCrossSectionSysClone[i],Form("%s x10^{%d}",nameMeasGlobal[i].Data(),legendScalingFactor),"pf");
              
            // statistics
            for (int j=0;j<graphEtaInvCrossSectionStatClone[i]->GetN();j++){
                graphEtaInvCrossSectionStatClone[i]->GetY()[j] *= scalingFactor;
                graphEtaInvCrossSectionStatClone[i]->GetEYhigh()[j] *= scalingFactor;
                graphEtaInvCrossSectionStatClone[i]->GetEYlow()[j] *= scalingFactor;
            }
            DrawGammaSetMarkerTGraph(graphEtaInvCrossSectionStatClone[i],  markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            graphEtaInvCrossSectionStatClone[i]->Draw("p,same,z");
            scalingFactor/=10;
            // legendScalingFactor--;
          }
        }
        // graphEtaInvCrossSectionSysInterpolation[9]     ->Draw("p,same,z");
        legendCrossSectionEta->Draw();


        // TLegend* legendEtaErr2 = GetAndSetLegend2(0.72, 0.72, 0.98, 0.72+(2*textSizeLabelsRel),0.85*textSizeLabelsPixel);
        // legendEtaErr2->AddEntry(histoBlack, "stat. Err.","ple");
        // legendEtaErr2->AddEntry(graphGrey,  "syst. Err.","f");
        // legendEtaErr2->Draw();

        drawLatexAdd("ALICE",0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem5TeV,0.93,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("#eta #rightarrow #gamma#gamma",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        histoDummyCrossSectionEta->Draw("sameaxis");
    canvasCrossSectionPi0->Update();
    canvasCrossSectionPi0->Print(Form("%s/Eta_InvariantCrossSectionMeas.%s",outputDir.Data(),suffix.Data()));
    
    // **********************************************************************************************************************
    // ******************************************* RAW YIELDS                        ****************************************
    // **********************************************************************************************************************
    
    TH2F * histoDummyPi0RawYields;
        histoDummyPi0RawYields                = new TH2F("histoDummyPi0RawYields", "histoDummyPi0RawYields",1000, minX,  maxX, 1000, 2e-9, 8e-2 );
    SetStyleHistoTH2ForGraphs( histoDummyPi0RawYields, "#it{p}_{T} (GeV/#it{c})", "RAW Yield per Event",
                            0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.4);//(#times #epsilon_{pur})
    
    histoDummyPi0RawYields->GetYaxis()->SetLabelOffset(0.001);
    histoDummyPi0RawYields->GetXaxis()->SetLabelOffset(-0.01);
    histoDummyPi0RawYields->GetXaxis()->SetMoreLogLabels(kTRUE);
    histoDummyPi0RawYields->DrawCopy();
    scalingFactor = pow(10,totalMeasAvail-1);
    legendScalingFactor = totalMeasAvail-1;
    for(Int_t i=numbersofmeas;i>-1;i--){
      if(availableMeas[i]){
        // statistics
        for (int j=0;j<graphPi0RawYields[i]->GetN();j++){
            graphPi0RawYields[i]->GetY()[j] *= scalingFactor;
            graphPi0RawYields[i]->GetEYhigh()[j] *= scalingFactor;
            graphPi0RawYields[i]->GetEYlow()[j] *= scalingFactor;
        }
        DrawGammaSetMarkerTGraph(graphPi0RawYields[i],  markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
        graphPi0RawYields[i]->Draw("p,same,z");
        scalingFactor/=10;
        legendScalingFactor--;
      }
    }
    legendCrossSectionPi0->Draw();
    drawLatexAdd("ALICE",0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(collisionSystem5TeV,0.93,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    histoDummyPi0RawYields->Draw("sameaxis");
    canvasCrossSectionPi0->Update();
    canvasCrossSectionPi0->Print(Form("%s/Pi0_RawYields.%s",outputDir.Data(),suffix.Data()));
    
    TH2F * histoDummyEtaRawYields;
        histoDummyEtaRawYields                = new TH2F("histoDummyEtaRawYields", "histoDummyEtaRawYields",1000, minX,  maxX, 1000, 2e-8, 2e-3 );
    SetStyleHistoTH2ForGraphs( histoDummyEtaRawYields, "#it{p}_{T} (GeV/#it{c})", "RAW Yield per Event",
                            0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.4);//(#times #epsilon_{pur})
    
    histoDummyEtaRawYields->GetYaxis()->SetLabelOffset(0.001);
    histoDummyEtaRawYields->GetXaxis()->SetLabelOffset(-0.01);
    histoDummyEtaRawYields->GetXaxis()->SetMoreLogLabels(kTRUE);
    histoDummyEtaRawYields->DrawCopy();
    scalingFactor = pow(10,totalMeasAvail-1);
    legendScalingFactor = totalMeasAvail-1;
    for(Int_t i=numbersofmeas;i>-1;i--){
      if(availableMeas[i]){
        // statistics
        for (int j=0;j<graphEtaRawYields[i]->GetN();j++){
            graphEtaRawYields[i]->GetY()[j] *= scalingFactor;
            graphEtaRawYields[i]->GetEYhigh()[j] *= scalingFactor;
            graphEtaRawYields[i]->GetEYlow()[j] *= scalingFactor;
        }
        DrawGammaSetMarkerTGraph(graphEtaRawYields[i],  markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
        graphEtaRawYields[i]->Draw("p,same,z");
        scalingFactor/=10;
        legendScalingFactor--;
      }
    }
    legendCrossSectionPi0->Draw();
    drawLatexAdd("ALICE",0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd(collisionSystem5TeV,0.93,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    histoDummyEtaRawYields->Draw("sameaxis");
    canvasCrossSectionPi0->Update();
    canvasCrossSectionPi0->Print(Form("%s/Eta_RawYields.%s",outputDir.Data(),suffix.Data()));
    
    // **********************************************************************************************************************
    // ******************************************* Compare to Interpolation           ****************************************
    // **********************************************************************************************************************
    TCanvas* canvasRatioToCombFit               = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToCombFit, 0.12, 0.01, 0.01, 0.12);
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
    histo2DRatioToCombFit                       = new TH2F("histo2DRatioToCombFit","histo2DRatioToCombFit",1000,minX,maxX,1000,0.01,5.    );
    SetStyleHistoTH2ForGraphs(histo2DRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","#frac{measurement}{interpolation}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRatioToCombFit->GetYaxis()->SetRangeUser(0.71,1.49);
    histo2DRatioToCombFit->Draw("copy");
    TGraphAsymmErrors* measDummyGraph              = (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[2]->Clone("dummypi0EMC");
    while(measDummyGraph->GetX()[0] < 1.5) measDummyGraph->RemovePoint(0);
    TGraphAsymmErrors* interpolationRatioPCMEMCStatErrBand = CalculateAsymGraphRatioToGraph(measDummyGraph,measDummyGraph);
    DrawGammaSetMarkerTGraphAsym(interpolationRatioPCMEMCStatErrBand, markerStyleDet[2], markerSizeDet[2]*0.55, colorDetMC[2] , colorDetMC[2], widthLinesBoxes, kTRUE);
    interpolationRatioPCMEMCStatErrBand     ->Draw("E2same");
    
    TGraphAsymmErrors* interpolationRatioPCMEMC = CalculateAsymGraphRatioToGraph(measDummyGraph,graphPi0InvCrossSectionStatInterpolation[2]);
    DrawGammaSetMarkerTGraph(interpolationRatioPCMEMC,  markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
    interpolationRatioPCMEMC->Draw("p,same,e1");

    DrawGammaLines(minX,maxX , 1., 1.,0.1, kGray+2);
    DrawGammaLines(minX,maxX , 1.1, 1.1,0.1, kGray, 7);
    DrawGammaLines(minX,maxX , 0.9, 0.9,0.1, kGray, 7);


    drawLatexAdd(collisionSystem5TeV.Data(),0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE);
    drawLatexAdd("EMC",0.15,0.82,textSizeLabelsRel,kFALSE,kFALSE);

    canvasRatioToCombFit->SaveAs(Form("%s/InterpolationRatioPi0EMC.%s",outputDir.Data(),suffix.Data()));
    
    histo2DRatioToCombFit->GetYaxis()->SetRangeUser(0.67,1.44);
      histo2DRatioToCombFit->Draw("copy");
      measDummyGraph              = (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[4]->Clone("dummypi0EMC");
      while(measDummyGraph->GetX()[0] < 0.9)measDummyGraph->RemovePoint(0);
      TGraphAsymmErrors* interpolationRatioEMCStatErrBand = CalculateAsymGraphRatioToGraph(measDummyGraph,measDummyGraph);
      DrawGammaSetMarkerTGraphAsym(interpolationRatioEMCStatErrBand, markerStyleDet[4], markerSizeDet[4]*0.55, colorDetMC[4] , colorDetMC[4], widthLinesBoxes, kTRUE);
      interpolationRatioEMCStatErrBand     ->Draw("E2same");
      TGraphAsymmErrors* interpolationRatioEMC = CalculateAsymGraphRatioToGraph(measDummyGraph,graphPi0InvCrossSectionStatInterpolation[4]);
      DrawGammaSetMarkerTGraph(interpolationRatioEMC,  markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
      interpolationRatioEMC->Draw("p,same,e1");
      
      DrawGammaLines(minX,maxX , 1., 1.,0.1, kGray+2);
      DrawGammaLines(minX,maxX , 1.1, 1.1,0.1, kGray, 7);
      DrawGammaLines(minX,maxX , 0.9, 0.9,0.1, kGray, 7);
      
      drawLatexAdd(collisionSystem5TeV.Data(),0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.93,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("PCM-EMC",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      
      canvasRatioToCombFit->SaveAs(Form("%s/InterpolationRatioPi0PCMEMC.%s",outputDir.Data(),suffix.Data()));
      
      
      TH2F * histo2DRatioToCombFitEta;
      histo2DRatioToCombFitEta                       = new TH2F("histo2DRatioToCombFitEta","histo2DRatioToCombFitEta",1000,minX,maxX,1000,0.31,1.99    );
      SetStyleHistoTH2ForGraphs(histo2DRatioToCombFitEta, "#it{p}_{T} (GeV/#it{c})","#frac{measurement}{interpolation}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
      histo2DRatioToCombFitEta->GetXaxis()->SetMoreLogLabels();
      histo2DRatioToCombFitEta->GetXaxis()->SetLabelOffset(-0.01);
  
      histo2DRatioToCombFitEta->Draw("copy");
      measDummyGraph              = (TGraphAsymmErrors*)graphEtaInvCrossSectionStat[4]->Clone("dummypi0EMC");
      // while(measDummyGraph->GetX()[0] < 1.9)measDummyGraph->RemovePoint(0);
      while(measDummyGraph->GetX()[measDummyGraph->GetN()-1] > 11.) measDummyGraph->RemovePoint(measDummyGraph->GetN()-1);
      TGraphAsymmErrors* interpolationRatioEMCStatErrBandEta = CalculateAsymGraphRatioToGraph(measDummyGraph,measDummyGraph);
      DrawGammaSetMarkerTGraphAsym(interpolationRatioEMCStatErrBandEta, markerStyleDet[4], markerSizeDet[4]*0.55, colorDetMC[4] , colorDetMC[4], widthLinesBoxes, kTRUE);
      interpolationRatioEMCStatErrBandEta     ->Draw("E2same");
      TGraphAsymmErrors* interpolationRatioPCMEMCEta = CalculateAsymGraphRatioToGraph(measDummyGraph,graphEtaInvCrossSectionStatInterpolation[4]);
      DrawGammaSetMarkerTGraph(interpolationRatioPCMEMCEta,  markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
interpolationRatioPCMEMCEta->Draw("p,same,e1");
      
      DrawGammaLines(minX,maxX , 1., 1.,0.1, kGray+2);
      DrawGammaLines(minX,maxX , 1.1, 1.1,0.1, kGray, 7);
      DrawGammaLines(minX,maxX , 0.9, 0.9,0.1, kGray, 7);
      
      drawLatexAdd(collisionSystem5TeV.Data(),0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("#eta #rightarrow #gamma#gamma",0.93,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("PCM-EMC",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      canvasRatioToCombFit->SaveAs(Form("%s/InterpolationRatioEtaPCMEMC.%s",outputDir.Data(),suffix.Data()));
      
      
      
      histo2DRatioToCombFitEta->Draw("copy");
      measDummyGraph              = (TGraphAsymmErrors*)graphEtaInvCrossSectionStatInterpolation[2]->Clone("dummypi0EMC");
      while(measDummyGraph->GetX()[0] < 2.0)measDummyGraph->RemovePoint(0);
      TGraphAsymmErrors* interpolationRatioPCMEMCStatErrBandEta = CalculateAsymGraphRatioToGraph(graphEtaInvCrossSectionStat[2],graphEtaInvCrossSectionStat[2]);
      DrawGammaSetMarkerTGraphAsym(interpolationRatioPCMEMCStatErrBandEta, markerStyleDet[2], markerSizeDet[2]*0.55, colorDetMC[2] , colorDetMC[2], widthLinesBoxes, kTRUE);
      interpolationRatioPCMEMCStatErrBandEta     ->Draw("E2same");
      
      TGraphAsymmErrors* interpolationRatioEMCEta = CalculateAsymGraphRatioToGraph(graphEtaInvCrossSectionStat[2],measDummyGraph);
      DrawGammaSetMarkerTGraph(interpolationRatioEMCEta,  markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
      interpolationRatioEMCEta->Draw("p,same,e1");
      
      DrawGammaLines(minX,maxX , 1., 1.,0.1, kGray+2);
      DrawGammaLines(minX,maxX , 1.1, 1.1,0.1, kGray, 7);
      DrawGammaLines(minX,maxX , 0.9, 0.9,0.1, kGray, 7);
      
      drawLatexAdd(collisionSystem5TeV.Data(),0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("#eta #rightarrow #gamma#gamma",0.93,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      drawLatexAdd("EMC",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
      canvasRatioToCombFit->SaveAs(Form("%s/InterpolationRatioEtaEMC.%s",outputDir.Data(),suffix.Data()));
      
      
      
      // RATIO TO COMBINED INTERPOLATION Pi0
      histo2DRatioToCombFit->GetYaxis()->SetTitle("#frac{indiv. measurement}{comb. interpolation}");
      histo2DRatioToCombFit->Draw("copy");
      TGraphAsymmErrors* interpolationRatioEMCStatErrBandComb = CalculateAsymGraphRatioToGraph(graphPi0InvCrossSectionStat[2],graphPi0InvCrossSectionStat[2]);
      DrawGammaSetMarkerTGraphAsym(interpolationRatioEMCStatErrBandComb, markerStyleDet[2], markerSizeDet[2]*0.55, colorDetMC[2] , colorDetMC[2], widthLinesBoxes, kTRUE);
      interpolationRatioEMCStatErrBandComb     ->Draw("E2same");
      TGraphAsymmErrors* interpolationRatioPCMEMCStatErrBandComb = CalculateAsymGraphRatioToGraph(graphPi0InvCrossSectionStat[4],graphPi0InvCrossSectionStat[4]);
      DrawGammaSetMarkerTGraphAsym(interpolationRatioPCMEMCStatErrBandComb, markerStyleDet[4], markerSizeDet[4]*0.55, colorDetMC[4] , colorDetMC[4], widthLinesBoxes, kTRUE);
      interpolationRatioPCMEMCStatErrBandComb     ->Draw("E2same");
      
      TGraphAsymmErrors* measDummyGraphEMC              = (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[2]->Clone("dummypi0combint");
      while(measDummyGraphEMC->GetX()[0] < 0.35)measDummyGraphEMC->RemovePoint(0);
      TGraphAsymmErrors* interpolationRatioEMCToComb = CalculateAsymGraphRatioToGraph(measDummyGraphEMC,graphPi0InvCrossSectionStatInterpolation[9]);
      DrawGammaSetMarkerTGraph(interpolationRatioEMCToComb,  markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
      interpolationRatioEMCToComb->Draw("p,same,e1");
      
      TGraphAsymmErrors* measDummyGraphPCMEMC              = (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[4]->Clone("dummypi0combint");
      while(measDummyGraphPCMEMC->GetX()[0] < 0.35)measDummyGraphPCMEMC->RemovePoint(0);
      TGraphAsymmErrors* interpolationRatioPCMEMCToComb = CalculateAsymGraphRatioToGraph(measDummyGraphPCMEMC,graphPi0InvCrossSectionStatInterpolation[9]);
      DrawGammaSetMarkerTGraph(interpolationRatioPCMEMCToComb,  markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
      interpolationRatioPCMEMCToComb->Draw("p,same,e1");
      
      DrawGammaLines(minX,maxX , 1., 1.,0.1, kGray+2);
      DrawGammaLines(minX,maxX , 1.1, 1.1,0.1, kGray, 7);
      DrawGammaLines(minX,maxX , 0.9, 0.9,0.1, kGray, 7);
      
      TLegend* legendCrossSectionPi0Compare           = GetAndSetLegend2(0.35, 0.86, 0.63, 0.86+(2*textSizeLabelsRel),0.8*textSizeLabelsPixel);
      legendCrossSectionPi0Compare->AddEntry(interpolationRatioEMCToComb,"EMC","p");
      legendCrossSectionPi0Compare->AddEntry(interpolationRatioPCMEMCToComb,"PCM-EMC","p");
      legendCrossSectionPi0Compare->Draw();
      
      drawLatexAdd(collisionSystem5TeV.Data(),0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE);
      drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE);
      // drawLatexAdd("EMC",0.15,0.82,textSizeLabelsRel,kFALSE,kFALSE);
      
      canvasRatioToCombFit->SaveAs(Form("%s/InterpolationRatioPi0ToCombined.%s",outputDir.Data(),suffix.Data()));
      
      // RATIO TO COMBINED INTERPOLATION Eta
      histo2DRatioToCombFitEta->Draw("copy");
      histo2DRatioToCombFitEta->GetYaxis()->SetTitle("#frac{indiv. measurement}{comb. interpolation}");
      TGraphAsymmErrors* interpolationRatioEMCStatErrBandCombEta = CalculateAsymGraphRatioToGraph(graphEtaInvCrossSectionStat[2],graphEtaInvCrossSectionStat[2]);
      DrawGammaSetMarkerTGraphAsym(interpolationRatioEMCStatErrBandCombEta, markerStyleDet[2], markerSizeDet[2]*0.55, colorDetMC[2] , colorDetMC[2], widthLinesBoxes, kTRUE);
      interpolationRatioEMCStatErrBandCombEta     ->Draw("E2same");
      while(graphEtaInvCrossSectionStat[4]->GetX()[graphEtaInvCrossSectionStat[4]->GetN()-1] > 11.) graphEtaInvCrossSectionStat[4]->RemovePoint(graphEtaInvCrossSectionStat[4]->GetN()-1);
      TGraphAsymmErrors* interpolationRatioPCMEMCStatErrBandCombEta = CalculateAsymGraphRatioToGraph(graphEtaInvCrossSectionStat[4],graphEtaInvCrossSectionStat[4]);
      DrawGammaSetMarkerTGraphAsym(interpolationRatioPCMEMCStatErrBandCombEta, markerStyleDet[4], markerSizeDet[4]*0.55, colorDetMC[4] , colorDetMC[4], widthLinesBoxes, kTRUE);
      interpolationRatioPCMEMCStatErrBandCombEta     ->Draw("E2same");
      
      TGraphAsymmErrors* measDummyGraphEMCEta              = (TGraphAsymmErrors*)graphEtaInvCrossSectionStatInterpolation[9]->Clone("dummyEtacombint");
      while(measDummyGraphEMCEta->GetX()[0] < 2)measDummyGraphEMCEta->RemovePoint(0);
      TGraphAsymmErrors* interpolationRatioEMCToCombEta = CalculateAsymGraphRatioToGraph(graphEtaInvCrossSectionStat[2],measDummyGraphEMCEta);
      DrawGammaSetMarkerTGraph(interpolationRatioEMCToCombEta,  markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
      interpolationRatioEMCToCombEta->Draw("p,same,e1");
      
      TGraphAsymmErrors* measDummyGraphPCMEMCEta              = (TGraphAsymmErrors*)graphEtaInvCrossSectionStatInterpolation[9]->Clone("dummyEtacombint");
      while(measDummyGraphPCMEMCEta->GetX()[0] < 1.25)measDummyGraphPCMEMCEta->RemovePoint(0);
      
      TGraphAsymmErrors* interpolationRatioPCMEMCToCombEta = CalculateAsymGraphRatioToGraph(graphEtaInvCrossSectionStat[4],measDummyGraphPCMEMCEta);
      DrawGammaSetMarkerTGraph(interpolationRatioPCMEMCToCombEta,  markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
      interpolationRatioPCMEMCToCombEta->Draw("p,same,e1");
      
      DrawGammaLines(minX,maxX , 1., 1.,0.1, kGray+2);
      DrawGammaLines(minX,maxX , 1.1, 1.1,0.1, kGray, 7);
      DrawGammaLines(minX,maxX , 0.9, 0.9,0.1, kGray, 7);
      
      legendCrossSectionPi0Compare->Draw();
      drawLatexAdd(collisionSystem5TeV.Data(),0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE);
      drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE);
      // drawLatexAdd("EMC",0.15,0.82,textSizeLabelsRel,kFALSE,kFALSE);
      
      canvasRatioToCombFit->SaveAs(Form("%s/InterpolationRatioEtaToCombined.%s",outputDir.Data(),suffix.Data()));
      
      
    // **********************************************************************************************************************
    // ******************************************* NONLINEARITY COMPARISONS ***************************************
    // **********************************************************************************************************************
  
    Double_t minXNL = 0.45;
    Double_t maxXNL = 29;
    cout << "plotting Nonlinearity parameterizations" << endl;
    canvasRatioToCombFit->cd();
    TH2F * histo2DNonlinearityComp;
    histo2DNonlinearityComp                       = new TH2F("histo2DNonlinearityComp","histo2DNonlinearityComp",1000,minXNL,maxXNL,1000,0.955,1.055    );
    SetStyleHistoTH2ForGraphs(histo2DNonlinearityComp, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (data)} #GT / M_{#pi^{0} (PDG)}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DNonlinearityComp->GetXaxis()->SetMoreLogLabels();
    histo2DNonlinearityComp->GetXaxis()->SetLabelOffset(-0.01);

    histo2DNonlinearityComp->Draw("copy");
    TF1* fitDataMass2760GeV       = new TF1("fitDataMass2760GeV", "[0] + [1]*pow(x,[2])" ,minXNL,maxXNL);
    fitDataMass2760GeV->SetParameter(0, 1.1673716264 );
    fitDataMass2760GeV->SetParameter(1, -0.1853095466 );
    fitDataMass2760GeV->SetParameter(2, -0.0848801702 );
    // fitDataMass2760GeV->SetParameter(0, 1.0443938253 );
    // fitDataMass2760GeV->SetParameter(1, -0.0691830812 );
    // fitDataMass2760GeV->SetParameter(2, -0.1247555443 );
    DrawGammaSetMarkerTF1( fitDataMass2760GeV, 7, 2, kMagenta+2);
    fitDataMass2760GeV->Draw("same");
    
    TF1* fitDataMass5023GeV       = new TF1("fitDataMass5023GeV", "[0] + [1]*pow(x,[2])" ,minXNL,maxXNL);
    fitDataMass5023GeV->SetParameter(0, 1.1237705555 );
    fitDataMass5023GeV->SetParameter(1, -0.1328806845 );
    fitDataMass5023GeV->SetParameter(2, -0.1064755375 );
    // fitDataMass5023GeV->SetParameter(0, 1.0825118468 );
    // fitDataMass5023GeV->SetParameter(1, -0.1079736424 );
    // fitDataMass5023GeV->SetParameter(2, -0.0800008954 );
    // fitDataMass5023GeV->SetParameter(0, 1.1818406492 ); // after fix
    // fitDataMass5023GeV->SetParameter(1, -0.1999998957 );// after fix
    // fitDataMass5023GeV->SetParameter(2, -0.1434322871 );// after fix
    // fitDataMass5023GeV->SetParameter(0, 1.1497456392 );// after fix
    // fitDataMass5023GeV->SetParameter(1, -0.1999999732 );// after fix
    // fitDataMass5023GeV->SetParameter(2, -0.0839303140 );// after fix
    DrawGammaSetMarkerTF1( fitDataMass5023GeV, 8, 2, kOrange+2);
    fitDataMass5023GeV->Draw("same");
    
    TF1* fitDataMass7000GeV       = new TF1("fitDataMass7000GeV", "[0] + [1]*pow(x,[2])" ,minXNL,maxXNL);
    fitDataMass7000GeV->SetParameter(0, 1.1850179319 );
    fitDataMass7000GeV->SetParameter(1, -0.1999999950 );
    fitDataMass7000GeV->SetParameter(2, -0.0863054172 );
    // fitDataMass7000GeV->SetParameter(0, 1.1082846035 );
    // fitDataMass7000GeV->SetParameter(1, -0.1369968318 );
    // fitDataMass7000GeV->SetParameter(2, -0.0800000002 );
    DrawGammaSetMarkerTF1( fitDataMass7000GeV, 3, 2, kGreen+2);
    fitDataMass7000GeV->Draw("same");
    
    TF1* fitDataMass8000GeV       = new TF1("fitDataMass8000GeV", "[0] + [1]*pow(x,[2])" ,minXNL,maxXNL);
    fitDataMass8000GeV->SetParameter(0, 1.1814766150 );
    fitDataMass8000GeV->SetParameter(1, -0.1980098061 );
    fitDataMass8000GeV->SetParameter(2, -0.0854569214 );
    // fitDataMass8000GeV->SetParameter(0, 1.0654169768 );
    // fitDataMass8000GeV->SetParameter(1, -0.0935785719 );
    // fitDataMass8000GeV->SetParameter(2, -0.1137883054 );
    DrawGammaSetMarkerTF1( fitDataMass8000GeV, 2, 2, kBlue+2);
    fitDataMass8000GeV->Draw("same");

    TF1* fitDataMass5000GeVpPb       = new TF1("fitDataMass5000GeVpPb", "[0]-TMath::Exp(-[1]*x+[2])" ,minXNL,maxXNL);
    fitDataMass5000GeVpPb->SetParameter(0, 1.0255088817 );
    fitDataMass5000GeVpPb->SetParameter(1, 0.3070452373 );
    fitDataMass5000GeVpPb->SetParameter(2, -2.9149185308 );
    // fitDataMass5000GeVpPb->SetParameter(0, 0.9910691195 );
    // fitDataMass5000GeVpPb->SetParameter(1, 0.4901455923 );
    // fitDataMass5000GeVpPb->SetParameter(2, -3.6647921806 );
    DrawGammaSetMarkerTF1( fitDataMass5000GeVpPb, 4, 2, kCyan+2);
    fitDataMass5000GeVpPb->Draw("same");
    
    
  
    TLegend* legendNLPCMEMC           = GetAndSetLegend2(0.15, 0.76-textSizeLabelsRel, 0.53, 0.76+(4.5*textSizeLabelsRel),0.8*textSizeLabelsPixel);
    legendNLPCMEMC->AddEntry(fitDataMass8000GeV,"pp, #sqrt{#it{s}} = 8 TeV","l");
    legendNLPCMEMC->AddEntry(fitDataMass7000GeV,"pp, #sqrt{#it{s}} = 7 TeV","l");
    legendNLPCMEMC->AddEntry(fitDataMass5023GeV,"pp, #sqrt{#it{s}} = 5.02TeV","l");
    legendNLPCMEMC->AddEntry(fitDataMass2760GeV,"pp, #sqrt{#it{s}} = 2.76 TeV","l");
    legendNLPCMEMC->AddEntry(fitDataMass5000GeVpPb,"p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV","l");
    legendNLPCMEMC->Draw();
    
    
    DrawGammaLines(minXNL,maxXNL , 1., 1.,0.1, kGray+2, 7);
    // DrawGammaLines(minXNL,maxXNL , 1.1, 1.1,0.1, kGray, 7);
    // DrawGammaLines(minXNL,maxXNL , 0.9, 0.9,0.1, kGray, 7);
    
    drawLatexAdd("ALICE",0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.93,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("PCM-EMC",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    canvasRatioToCombFit->SaveAs(Form("%s/NonlinearityComparisonsPCMEMC.%s",outputDir.Data(),suffix.Data()));
    
    
    histo2DNonlinearityComp                       = new TH2F("histo2DNonlinearityComp","histo2DNonlinearityComp",1000,minXNL,maxXNL,1000,0.905,1.105    );
    SetStyleHistoTH2ForGraphs(histo2DNonlinearityComp, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (data)} #GT / M_{#pi^{0} (PDG)}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DNonlinearityComp->GetXaxis()->SetMoreLogLabels();
    histo2DNonlinearityComp->GetXaxis()->SetLabelOffset(-0.01);
    histo2DNonlinearityComp->Draw("copy");
    
    fitDataMass2760GeV->SetParameter(0, 1.0383412435 );
    fitDataMass2760GeV->SetParameter(1, -0.0851830429 );
    fitDataMass2760GeV->SetParameter(2, -0.4999999996 );
    // fitDataMass2760GeV->SetParameter(0, 0.9980625418 );
    // fitDataMass2760GeV->SetParameter(1, -0.0564782662 );
    // fitDataMass2760GeV->SetParameter(2, -0.5 );
    DrawGammaSetMarkerTF1( fitDataMass2760GeV, 7, 2, kMagenta+2);
    fitDataMass2760GeV->Draw("same");
    
    fitDataMass5023GeV->SetParameter(0, 1.0434825769 );
    fitDataMass5023GeV->SetParameter(1, -0.0788025699 );
    fitDataMass5023GeV->SetParameter(2, -0.5000000000 );
    // fitDataMass5023GeV->SetParameter(0, 1.0194962801 );
    // fitDataMass5023GeV->SetParameter(1, -0.0830205281 );
    // fitDataMass5023GeV->SetParameter(2, -0.4628036411 );
    DrawGammaSetMarkerTF1( fitDataMass5023GeV, 8, 2, kOrange+2);
    fitDataMass5023GeV->Draw("same");
    
    fitDataMass7000GeV->SetParameter(0, 1.1224162203 );
    fitDataMass7000GeV->SetParameter(1, -0.1586806096 );
    fitDataMass7000GeV->SetParameter(2, -0.2458351112 );
    // fitDataMass7000GeV->SetParameter(0, 1.0074002842 );
    // fitDataMass7000GeV->SetParameter(1, -0.0682543971 );
    // fitDataMass7000GeV->SetParameter(2, -0.4509341085 );
    DrawGammaSetMarkerTF1( fitDataMass7000GeV, 3, 2, kGreen+2);
    fitDataMass7000GeV->Draw("same");
    
    fitDataMass8000GeV->SetParameter(0, 1.1603460704 );
    fitDataMass8000GeV->SetParameter(1, -0.1999999989 );
    fitDataMass8000GeV->SetParameter(2, -0.2194447313 );
    // fitDataMass8000GeV->SetParameter(0, 1.1389201636 );
    // fitDataMass8000GeV->SetParameter(1, -0.1999994717 );
    // fitDataMass8000GeV->SetParameter(2, -0.1622237979 );
    DrawGammaSetMarkerTF1( fitDataMass8000GeV, 2, 2, kBlue+2);
    fitDataMass8000GeV->Draw("same");
    
    fitDataMass5000GeVpPb->SetParameter(0, 1.0165873637 );
    fitDataMass5000GeVpPb->SetParameter(1, 0.6999387334 );
    fitDataMass5000GeVpPb->SetParameter(2, -2.1324782465 );
    // fitDataMass5000GeVpPb->SetParameter(0, 0.9795532189 );
    // fitDataMass5000GeVpPb->SetParameter(1, 0.8578583955 );
    // fitDataMass5000GeVpPb->SetParameter(2, -2.3447892540 );
    DrawGammaSetMarkerTF1( fitDataMass5000GeVpPb, 4, 2, kCyan+2);
    fitDataMass5000GeVpPb->Draw("same");
    
    
    
    legendNLPCMEMC->Draw();
    
    
    DrawGammaLines(minXNL,maxXNL , 1., 1.,0.1, kGray+2, 7);
    // DrawGammaLines(minXNL,maxXNL , 1.1, 1.1,0.1, kGray, 7);
    // DrawGammaLines(minXNL,maxXNL , 0.9, 0.9,0.1, kGray, 7);
    
    drawLatexAdd("ALICE",0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.93,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("EMCal",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
    canvasRatioToCombFit->SaveAs(Form("%s/NonlinearityComparisonsEMC.%s",outputDir.Data(),suffix.Data()));
    // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0            ****************************************
    // **********************************************************************************************************************
    cout << "plotting pi0 mass" << endl;
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

    TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, minX, maxX ,1000., -30, 40);
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
        if(histoPi0FWHMMeV[i] && histoPi0TrueFWHMMeV[i] && availableMeas[i]){
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

    TLatex *labelMassPerf       = new TLatex(0.13,0.87,"ALICE");
    SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4);
    labelMassPerf->SetTextFont(43);
    labelMassPerf->Draw();
    TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystem5TeV.Data());
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

    TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, minX, maxX, 1000., 120.1, 160.9);
    SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
    histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
    histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
    histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
    histo2DAllPi0Mass->DrawCopy();

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0Mass[i] && histoPi0TrueMass[i] && availableMeas[i]){
            DrawGammaSetMarker(histoPi0Mass[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            histoPi0Mass[i]->Draw("p,same,e");
            DrawGammaSetMarker(histoPi0TrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
            histoPi0TrueMass[i]->Draw("p,same,e");
        }
    }

    DrawGammaLines(minX, maxX , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

    TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
    SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4);
    labelLegendBMass->SetTextFont(43);
    labelLegendBMass->Draw();

    //********************************** Defintion of the Legend **************************************************
    Double_t columnsLegendMass2[3]      = {0.,0.57,0.84};
    //         Double_t rowsLegendMass2[5] = {0.8,0.6,0.4,0.2,0.01};
    //         Double_t rowsLegendMass2[6] = {0.84,0.66,0.50,0.33,0.16,0.01};
    //           Double_t  rowsLegendMass2[7]= {0.84,0.66,0.50,0.33,0.01,0.16};
    Double_t  rowsLegendMass2[9]= {0.84,0.66,0.51,0.50,0.331,0.33,0.01,0.16}; //setting for use without PHOS and PCM-PHOS peak positions
    //******************* Offsets ***********************
    Double_t offsetMarkerXMass2         = 0.1;
    Double_t offsetMarkerYMass2         = 0.1;
    //****************** Scale factors ******************
    Double_t scaleMarkerMass2           = 1.2;

    padMassLegend1->cd();
    //****************** first Column **************************************************
    TLatex *textMassPCM[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0Mass[i] && histoPi0TrueMass[i] && histoPi0FWHMMeV[i] && histoPi0TrueFWHMMeV[i] && availableMeas[i]){
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
        if(histoPi0Mass[i] && histoPi0TrueMass[i]&& availableMeas[i]){
            markerPCMPi0Mass[i]             = CreateMarkerFromHisto(histoPi0Mass[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
            markerPCMPi0Mass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
            markerPCMPi0MassMC[i]           = CreateMarkerFromHisto(histoPi0TrueMass[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
            markerPCMPi0MassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
        }
    }

    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ******************************************* Mass and width for Eta            ****************************************
    // **********************************************************************************************************************
    cout << "plotting eta mass" << endl;
    textSizeLabelsPixel             = 50;
    canvasMassWidthPi0->cd();
    padWidthPi0->Draw();
    padMassPi0->Draw();
    padMassLegend1->Draw();

    padWidthPi0->cd();
    padWidthPi0->SetLogx();

    TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, minX, maxX ,1000., 1, 89);
    SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
    histo2DAllEtaFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllEtaFWHM->GetYaxis()->SetNdivisions(505);
    histo2DAllEtaFWHM->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllEtaFWHM->GetXaxis()->SetTickLength(0.05);
    histo2DAllEtaFWHM->GetYaxis()->SetTickLength(0.026);
    histo2DAllEtaFWHM->DrawCopy();

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoEtaFWHMMeV[i] && histoEtaTrueFWHMMeV[i] && availableMeas[i]){
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

    padMassPi0->cd();
    padMassPi0->SetLogx();

    TH2F * histo2DAllEtaMass            = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",20, minX, maxX, 1000.,  500.1, 584.9);
    SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
    histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllEtaMass->GetYaxis()->SetNdivisions(505);
    histo2DAllEtaMass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllEtaMass->GetXaxis()->SetTickLength(0.05);
    histo2DAllEtaMass->GetXaxis()->SetLabelOffset(-0.015);
    histo2DAllEtaMass->DrawCopy();

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoEtaMass[i] && histoEtaTrueMass[i] && availableMeas[i]){
            DrawGammaSetMarker(histoEtaMass[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            histoEtaMass[i]->Draw("p,same,e");
            DrawGammaSetMarker(histoEtaTrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
            histoEtaTrueMass[i]->Draw("p,same,e");
        }
    }

    DrawGammaLines(minX, maxX , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.1, kGray);

    labelLegendBMass->Draw();


    padMassLegend1->cd();
    //****************** first Column **************************************************
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoEtaMass[i] && histoEtaTrueMass[i] && histoEtaFWHMMeV[i] && histoEtaTrueFWHMMeV[i] && availableMeas[i]){
            textMassPCM[i]->Draw();
        }
    }
    //****************** second Column *************************************************
    textMassData->Draw();
    textMassMC->Draw();

    TMarker* markerPCMEtaMass[10];
    TMarker* markerPCMEtaMassMC[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoEtaMass[i] && histoEtaTrueMass[i]&& availableMeas[i]){
            markerPCMEtaMass[i]             = CreateMarkerFromHisto(histoEtaMass[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
            markerPCMEtaMass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
            markerPCMEtaMassMC[i]           = CreateMarkerFromHisto(histoEtaTrueMass[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
            markerPCMEtaMassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
        }
    }

    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Eta_MassAndWidth.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for pi0 single measurement  *********************************
    // **********************************************************************************************************************
    textSizeLabelsPixel             = 55;
    textSizeLabelsRel      = 55./1200;
    cout << textSizeLabelsRel << endl;

    TCanvas* canvasAcceptanceTimesEff       = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.1, 0.01, 0.015, 0.095);
    canvasAcceptanceTimesEff->SetLogy(1);
    canvasAcceptanceTimesEff->SetLogx(1);

    TH2F * histo2DAccEff;
    histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, minX,  maxX, 1000, 8e-5, 3 );
    SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),
                            0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);//(#times #epsilon_{pur})
                            histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
                            histo2DAccEff->GetXaxis()->SetLabelOffset(-0.01);
                            histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
                            histo2DAccEff->DrawCopy();

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0AccTimesEff[i] && availableMeas[i]){
            DrawGammaSetMarker(histoPi0AccTimesEff[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            histoPi0AccTimesEff[i]->Draw("p,same,e");
        }
    }

    TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.59, 0.13, 0.93, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0AccTimesEff[i] && availableMeas[i]){
            legendEffiAccPi0->AddEntry(histoPi0AccTimesEff[i],nameMeasGlobal[i].Data(),"p");
        }
    }
    legendEffiAccPi0->Draw();

    drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(collisionSystem5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));
    
    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for eta single measurement  *********************************
    // **********************************************************************************************************************
    TH2F * histo2DAccEffEta;
    histo2DAccEffEta                = new TH2F("histo2DAccEffEta", "histo2DAccEffEta",1000, minX,  maxX, 1000, 8e-5, 3 );
    SetStyleHistoTH2ForGraphs( histo2DAccEffEta, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),
                            0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);//(#times #epsilon_{pur})
                            histo2DAccEffEta->GetYaxis()->SetLabelOffset(0.001);
                            histo2DAccEffEta->GetXaxis()->SetLabelOffset(-0.01);
                            histo2DAccEffEta->GetXaxis()->SetMoreLogLabels(kTRUE);
                            histo2DAccEffEta->DrawCopy();

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoEtaAccTimesEff[i] && availableMeas[i]){
            DrawGammaSetMarker(histoEtaAccTimesEff[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            histoEtaAccTimesEff[i]->Draw("p,same,e");
        }
    }

    TLegend* legendEffiAccEta           = GetAndSetLegend2(0.59, 0.13, 0.93, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoEtaAccTimesEff[i] && availableMeas[i]){
            legendEffiAccEta->AddEntry(histoEtaAccTimesEff[i],nameMeasGlobal[i].Data(),"p");
        }
    }
    TLegend* legendEffiAccPi0Eta           = GetAndSetLegend2(0.59, 0.13, 0.93, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0AccTimesEff[i] && availableMeas[i]){
            legendEffiAccPi0Eta->AddEntry(histoEtaAccTimesEff[i],nameMeasGlobal[i].Data(),"p");
        }
    }
    legendEffiAccEta->Draw();

    drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(collisionSystem5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Eta_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));


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

    for (Int_t i =0 ; i < numbersofmeas; i++){
        if (haveAllPi0InvMass[i] && availableMeas[i]){
            canvasInvMassSamplePlot->cd();
            histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.255);
            histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSub[i]->GetMinimum(),1.15*histoPi0InvMassSigPlusBG[i]->GetMaximum());
            if(i==2)
                histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSub[i]->GetMinimum(),1.35*histoPi0InvMassSigPlusBG[i]->GetMaximum());
            histo2DPi0InvMassDummy->DrawCopy();

            // TLatex *labelInvMassPtRange = new TLatex(0.945,0.9,Form("#pi^{0}: %s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}",strLowerEdgeExamplePi0[i].Data(),strUpperEdgeExamplePi0[i].Data()));
            TLatex *labelInvMassPtRange = new TLatex(0.945,0.9,Form("#pi^{0}: %s",histoPi0InvMassSigRemBGSub[i]->GetTitle()));

            DrawGammaSetMarker(histoPi0InvMassSigPlusBG[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
            histoPi0InvMassSigPlusBG[i]->SetLineWidth(1);
            histoPi0InvMassSigPlusBG[i]->Draw("hist,e,same");
            DrawGammaSetMarker(histoPi0InvMassBGTot[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
            histoPi0InvMassBGTot[i]->Draw("same");

            DrawGammaSetMarker(histoPi0InvMassSigRemBGSub[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
            histoPi0InvMassSigRemBGSub[i]->Draw("same");
            fitPi0InvMassSig[i]->SetNpx(1000);
            fitPi0InvMassSig[i]->SetRange(0,0.255);
            fitPi0InvMassSig[i]->SetLineColor(fitColorInvMassSG);
            fitPi0InvMassSig[i]->SetLineWidth(1);
            fitPi0InvMassSig[i]->Draw("same");

            TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9,collisionSystem5TeV.Data());
            SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
            labelInvMassEnergy->SetTextFont(43);
            labelInvMassEnergy->Draw();

            TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.8*textsizeLabelsInvMass,"MinBias");
            SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
            labelInvMassTrigger->SetTextFont(43);
            labelInvMassTrigger->Draw();

            TLatex *labelInvMassReco  = new TLatex(0.135,0.9-2*0.8*textsizeLabelsInvMass,Form("%s",nameMeasGlobal[i].Data()));
            SetStyleTLatex( labelInvMassReco, 0.85*textSizeLabelsPixel,4);
            labelInvMassReco->SetTextFont(43);
            labelInvMassReco->Draw();

            SetStyleTLatex( labelInvMassPtRange, 0.85*textSizeLabelsPixel,4);
            labelInvMassPtRange->SetTextAlign(31);
            labelInvMassPtRange->SetTextFont(43);
            labelInvMassPtRange->Draw();

            TLegend* legendInvMass  = GetAndSetLegend2(0.67, 0.88-5*0.75*textsizeLabelsInvMass, 0.9, 0.88, 0.85*textSizeLabelsPixel);
            legendInvMass->SetMargin(0.25);
            legendInvMass->AddEntry(histoPi0InvMassSigPlusBG[i],"Raw real events","l");
            if(i!=2){
                legendInvMass->AddEntry(histoPi0InvMassBGTot[i],"Mixed event +","p");
                legendInvMass->AddEntry((TObject*)0,"corr. BG","");
            } else{
                legendInvMass->AddEntry(histoPi0InvMassBGTot[i],"Mixed event","p");
            }
            legendInvMass->AddEntry(histoPi0InvMassSigRemBGSub[i],"BG subtracted","p");
            legendInvMass->AddEntry(fitPi0InvMassSig[i], "Fit","l");
            legendInvMass->Draw();
            canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBin_%s.%s",outputDir.Data(), nameMeasGlobal[i].Data(), suffix.Data()));
        } else {
            cout << "missing partial input for invariant mass bin for pi0 for meas: " << nameMeasGlobal[i].Data() << endl;
        }
        if (haveAllEtaInvMass[i] && availableMeas[i]){
            canvasInvMassSamplePlot->cd();
            histo2DEtaInvMassDummy->GetXaxis()->SetRangeUser(0.02,0.255);
            histo2DEtaInvMassDummy->GetYaxis()->SetRangeUser(histoEtaInvMassSigRemBGSub[i]->GetMinimum(),1.15*histoEtaInvMassSigPlusBG[i]->GetMaximum());
            if(i==2)
                histo2DEtaInvMassDummy->GetYaxis()->SetRangeUser(histoEtaInvMassSigRemBGSub[i]->GetMinimum(),1.35*histoEtaInvMassSigPlusBG[i]->GetMaximum());
            histo2DEtaInvMassDummy->DrawCopy();

            // TLatex *labelInvMassPtRange = new TLatex(0.945,0.9,Form("#pi^{0}: %s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}",strLowerEdgeExampleEta[i].Data(),strUpperEdgeExampleEta[i].Data()));
            TLatex *labelInvMassPtRange = new TLatex(0.945,0.9,Form("#eta: %s",histoEtaInvMassSigRemBGSub[i]->GetTitle()));

            DrawGammaSetMarker(histoEtaInvMassSigPlusBG[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
            histoEtaInvMassSigPlusBG[i]->SetLineWidth(1);
            histoEtaInvMassSigPlusBG[i]->Draw("hist,e,same");
            DrawGammaSetMarker(histoEtaInvMassBGTot[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
            histoEtaInvMassBGTot[i]->Draw("same");

            DrawGammaSetMarker(histoEtaInvMassSigRemBGSub[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
            histoEtaInvMassSigRemBGSub[i]->Draw("same");
            fitEtaInvMassSig[i]->SetNpx(1000);
            fitEtaInvMassSig[i]->SetRange(0.35,0.695);
            fitEtaInvMassSig[i]->SetLineColor(fitColorInvMassSG);
            fitEtaInvMassSig[i]->SetLineWidth(1);
            fitEtaInvMassSig[i]->Draw("same");

            TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9,collisionSystem5TeV.Data());
            SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
            labelInvMassEnergy->SetTextFont(43);
            labelInvMassEnergy->Draw();

            TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.8*textsizeLabelsInvMass,"MinBias");
            SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
            labelInvMassTrigger->SetTextFont(43);
            labelInvMassTrigger->Draw();

            TLatex *labelInvMassReco  = new TLatex(0.135,0.9-2*0.8*textsizeLabelsInvMass,Form("%s",nameMeasGlobal[i].Data()));
            SetStyleTLatex( labelInvMassReco, 0.85*textSizeLabelsPixel,4);
            labelInvMassReco->SetTextFont(43);
            labelInvMassReco->Draw();

            SetStyleTLatex( labelInvMassPtRange, 0.85*textSizeLabelsPixel,4);
            labelInvMassPtRange->SetTextAlign(31);
            labelInvMassPtRange->SetTextFont(43);
            labelInvMassPtRange->Draw();

            TLegend* legendInvMass  = GetAndSetLegend2(0.67, 0.88-5*0.75*textsizeLabelsInvMass, 0.9, 0.88, 0.85*textSizeLabelsPixel);
            legendInvMass->SetMargin(0.25);
            legendInvMass->AddEntry(histoEtaInvMassSigPlusBG[i],"Raw real events","l");
            if(i!=2){
                legendInvMass->AddEntry(histoEtaInvMassBGTot[i],"Mixed event +","p");
                legendInvMass->AddEntry((TObject*)0,"corr. BG","");
            } else{
                legendInvMass->AddEntry(histoEtaInvMassBGTot[i],"Mixed event","p");
            }
            legendInvMass->AddEntry(histoEtaInvMassSigRemBGSub[i],"BG subtracted","p");
            legendInvMass->AddEntry(fitEtaInvMassSig[i], "Fit","l");
            legendInvMass->Draw();
            canvasInvMassSamplePlot->SaveAs(Form("%s/Eta_InvMassBin_%s.%s",outputDir.Data(), nameMeasGlobal[i].Data(), suffix.Data()));
        } else {
            cout << "missing partial input for invariant mass bin for eta for meas: " << nameMeasGlobal[i].Data() << endl;
        }
    }
    


}