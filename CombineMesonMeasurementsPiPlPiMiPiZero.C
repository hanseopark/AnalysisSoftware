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

/*void drawLatexAdd(TString latextext, Double_t textcolumn, Double_t textrow, Double_t textSizePixel,Bool_t setFont = kFALSE, Bool_t setFont2 = kFALSE, Bool_t alignRight = kFALSE, Color_t textcolor = kBlack){
    TLatex *latexDummy                  = new TLatex(textcolumn ,textrow,latextext);
    SetStyleTLatex( latexDummy, textSizePixel,4);
    if(setFont)
        latexDummy->SetTextFont(62);
    if(setFont2)
        latexDummy->SetTextFont(43);
    if(alignRight)
        latexDummy->SetTextAlign(31);
    latexDummy->SetTextColor(textcolor);
    latexDummy->Draw();
}*/


void CombineMesonMeasurementsPiPlPiMiPiZero(      TString fileNamePCM     = "",
                                        TString fileNamePHOS    = "",
                                        TString fileNameEMCal   = "/home/nschmidt/AnalysisSoftware/pdf/5TeV/2017_08_19/FinalResultsTriggersPatched_EMC/data_EMCAL-EMCALResultsFullCorrection_PP.root",
                                        TString fileNamePCMPHOS = "",
                                        TString fileNamePCMEMCal= "/home/nschmidt/AnalysisSoftware/pdf/5TeV/2017_08_20/FinalResultsTriggersPatched/data_PCM-EMCALResultsFullCorrection_PP.root",
                                        TString fileNamePi0PCM    = "",
                                        TString fileNamePi0PHOS    = "",
                                        TString fileNamePi0EMCal   = "/home/nschmidt/AnalysisSoftware/pdf/5TeV/2017_08_19/FinalResultsTriggersPatched_EMC/data_EMCAL-EMCALResultsFullCorrection_PP.root",
                                        TString fileNamePi0PCMPHOS = "",
                                        TString fileNamePi0PCMEMCal= "/home/nschmidt/AnalysisSoftware/pdf/5TeV/2017_08_20/FinalResultsTriggersPatched/data_PCM-EMCALResultsFullCorrection_PP.root",
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
    TString collisionSystem5TeV               = "pp, #sqrt{#it{s}} = 7 TeV";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements7TeV",suffix.Data(),dateForOutput.Data());
    cout << outputDir.Data() << endl;
    cout << fileNamePCM.Data() << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    if(fileNamePCM.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    if(fileNamePHOS.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    if(fileNameEMCal.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputEMCal.root", fileNameEMCal.Data(), outputDir.Data()));
    if(fileNamePCMPHOS.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputPCMPHOS.root", fileNamePCMPHOS.Data(), outputDir.Data()));
    if(fileNamePCMEMCal.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputPCMEMCal.root", fileNamePCMEMCal.Data(), outputDir.Data()));

    if(fileNamePi0PCM.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputPi0PCM.root", fileNamePi0PCM.Data(), outputDir.Data()));
    if(fileNamePi0PHOS.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputPi0PHOS.root", fileNamePi0PHOS.Data(), outputDir.Data()));
    if(fileNamePi0EMCal.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputPi0EMCal.root", fileNamePi0EMCal.Data(), outputDir.Data()));
    if(fileNamePi0PCMPHOS.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputPi0PCMPHOS.root", fileNamePi0PCMPHOS.Data(), outputDir.Data()));
    if(fileNamePi0PCMEMCal.CompareTo(""))
      gSystem->Exec(Form("cp %s %s/InputPi0PCMEMCal.root", fileNamePi0PCMEMCal.Data(), outputDir.Data()));


    Double_t mesonMassExpectOmega                = TDatabasePDG::Instance()->GetParticle(223)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();

    Width_t widthLinesBoxes                     = 1.4;

    // Definition of colors, styles and markers sizes
    Color_t colorComb                           = kBlue+2;
    Style_t markerStyleComb                     = 20;
    Size_t  markerSizeComb                      = 2;

    TString nameMeasGlobal[11]                  = {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCMOtherDataset"};
    TString nameMeasGlobalPi0[11]                  = {"PCM (#pi^{0})", "PHOS (#pi^{0})", "EMCal (#pi^{0})", "PCM-PHOS (#pi^{0})", "PCM-EMCal (#pi^{0})", "PCM-Dalitz (#pi^{0})",
                                                      "PHOS-Dalitz (#pi^{0})", "EMCal-Dalitz (#pi^{0})", "EMCal high pT (#pi^{0})", "EMCal merged (#pi^{0})", "PCMOtherDataset (#pi^{0})"};
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
    TFile* inputFilePi0[10];
        inputFile[0]                            = new TFile(fileNamePCM.Data());
        inputFile[1]                            = new TFile(fileNamePHOS.Data());
        inputFile[2]                            = new TFile(fileNameEMCal.Data());
        inputFile[3]                            = new TFile(fileNamePCMPHOS.Data());
        inputFile[4]                            = new TFile(fileNamePCMEMCal.Data());

        // Pi0 input files for comparison
        inputFilePi0[0]                            = new TFile(fileNamePi0PCM.Data());
        inputFilePi0[1]                            = new TFile(fileNamePi0PHOS.Data());
        inputFilePi0[2]                            = new TFile(fileNamePi0EMCal.Data());
        inputFilePi0[3]                            = new TFile(fileNamePi0PCMPHOS.Data());
        inputFilePi0[4]                            = new TFile(fileNamePi0PCMEMCal.Data());

    TDirectory* directoryOmega[10];
    TDirectory* directoryEta[10];
    TDirectory* directoryPi0[10];
    for(Int_t i=0;i<numbersofmeas;i++){
      if(!inputFile[i]->IsZombie()){
        cout << "loading directories for " <<  nameMeasGlobal[i] << endl;
        directoryOmega[i]                     = (TDirectory*)inputFile[i]->Get("Omega7TeV");
        directoryEta[i]                     = (TDirectory*)inputFile[i]->Get("Eta7TeV");
        directoryPi0[i]                     = (TDirectory*)inputFilePi0[i]->Get("Pi07TeV");
      }
    }
    cout << __LINE__<<endl;
    TH1D* histoOmegaMass[10];
    TH1D* histoOmegaFWHMMeV[10];
    TH1D* histoOmegaTrueMass[10];
    TH1D* histoOmegaTrueFWHMMeV[10];
    TH1D* histoEtaMass[10];
    TH1D* histoEtaFWHMMeV[10];
    TH1D* histoEtaTrueMass[10];
    TH1D* histoEtaTrueFWHMMeV[10];
    TH1D* histoOmegaAcc[10];
    TH1D* histoOmegaTrueEffPt[10];
    TH1D* histoOmegaAccTimesEff[10];
    TH1D* histoEtaAcc[10];
    TH1D* histoEtaTrueEffPt[10];
    TH1D* histoEtaAccTimesEff[10];
    TH1D* histoOmegaInvCrossSection[10];
    TH1D* histoEtaInvCrossSection[10];
    TH1D* histoPi0Acc[10];
    TH1D* histoPi0TrueEffPt[10];
    TH1D* histoPi0AccTimesEff[10];
    TGraphAsymmErrors* graphOmegaInvCrossSectionSys[10];
    TGraphAsymmErrors* graphOmegaInvCrossSectionStat[10];
    TGraphAsymmErrors* graphEtaInvCrossSectionStat[10];
    TGraphAsymmErrors* graphEtaInvCrossSectionSys[10];
    TGraphAsymmErrors* graphEtaToOmegaPCMStat;
    TGraphAsymmErrors* graphEtaToOmegaPCMSys;

    Double_t rapidityMeas[10]                   = {1.6, 1,1, 1.6,1,1,1};
    Double_t availableMeas[10]                  = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
    Double_t availableMeasPi0[10]                  = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
    Int_t nMeasSet                              = 0;

    for (Int_t i = 0; i < 5; i++){
      if(inputFile[i]->IsZombie())
        continue;
      else
        cout << "loading  histos for " << nameMeasGlobal[i] << endl;
        
    //______________________________ Neutral pion cross section
      graphOmegaInvCrossSectionSys[i]       = (TGraphAsymmErrors*)directoryOmega[i]->Get("InvCrossSectionOmegaSys");
      histoOmegaInvCrossSection[i]          = (TH1D*)directoryOmega[i]->Get("InvCrossSectionOmega");
      graphOmegaInvCrossSectionStat[i]      = new TGraphAsymmErrors(histoOmegaInvCrossSection[i]);
      if(!graphOmegaInvCrossSectionSys[i] || !histoOmegaInvCrossSection[i]){
        cout << "missing invariant cross section histograms of Omega... returning!" << endl;
        return;
      } else {
        availableMeas[i]                  = kTRUE;
        nMeasSet++;
        cout << "loaded " << nameMeasGlobal[i] << " Omega invariant cross section" << endl;
      }
      
      //______________________________ Neutral pion invariant mass peak pos and FWHM
      histoOmegaMass[i]                     = (TH1D*)directoryOmega[i]->Get("Omega_Mass_data_INT1");
      histoOmegaFWHMMeV[i]                  = (TH1D*)directoryOmega[i]->Get("Omega_Width_data_INT1");
      histoOmegaTrueMass[i]                 = (TH1D*)directoryOmega[i]->Get("Omega_Mass_MC_INT1");
      histoOmegaTrueFWHMMeV[i]              = (TH1D*)directoryOmega[i]->Get("Omega_Width_MC_INT1");
      if(!histoOmegaMass[i] || !histoOmegaFWHMMeV[i] || !histoOmegaFWHMMeV[i] || !histoOmegaFWHMMeV[i]){
        cout << "missing mass or width histograms... returning!" << endl;
        return;
      } else {
        // scaling peak pos and FWHM to go from GeV/c2 to MeV/c2
        histoOmegaMass[i]                   ->Scale(1000);
        histoOmegaFWHMMeV[i]                ->Scale(1000);
        histoOmegaTrueMass[i]               ->Scale(1000);
        histoOmegaTrueFWHMMeV[i]            ->Scale(1000);
        cout << "loaded " << nameMeasGlobal[i] << " mass and width" << endl;
      }
      
      //______________________________ Neutral pion acceptance and efficiency and calculate acc*eff*y*2pi
      histoOmegaAcc[i]                      = (TH1D*)directoryOmega[i]->Get("AcceptanceOmega_INT1");
      histoOmegaTrueEffPt[i]                = (TH1D*)directoryOmega[i]->Get("EfficiencyOmega_INT1");
      // calculating for pi0 for later comparison

      if(!histoOmegaAcc[i] || !histoOmegaTrueEffPt[i]){
        cout << "missing acceptance or efficiency histograms... returning!" << endl;
        return;
      } else {
        // calculating acceptance times efficiency
        histoOmegaAccTimesEff[i]            = (TH1D*)histoOmegaTrueEffPt[i]->Clone(Form("histoOmegaAccTimesEff%s",nameMeasGlobal[i].Data()));
        histoOmegaAccTimesEff[i]->Multiply(histoOmegaAcc[i]);
        histoOmegaAccTimesEff[i]->Scale(2*TMath::Pi()*rapidityMeas[i]);
        cout << "loaded " << nameMeasGlobal[i] << "acceptance and efficiency" << endl;
      }

      if(!inputFilePi0[i]->IsZombie()){
          histoPi0Acc[i]                      = (TH1D*)directoryPi0[i]->Get("AcceptancePi0_INT1");
          histoPi0TrueEffPt[i]                = (TH1D*)directoryPi0[i]->Get("EfficiencyPi0_INT1");

          if(!histoPi0Acc[i] || !histoPi0TrueEffPt[i]){
              cout << "missing acceptance or efficiency histograms of pi0 ... returning!" << endl;
              return;
          } else {
              // calculating acceptance times efficiency
              availableMeasPi0[i]=kTRUE;
              histoPi0AccTimesEff[i]            = (TH1D*)histoPi0TrueEffPt[i]->Clone(Form("histoPi0AccTimesEff%s",nameMeasGlobalPi0[i].Data()));
              histoPi0AccTimesEff[i]->Multiply(histoPi0Acc[i]);
              histoPi0AccTimesEff[i]->Scale(2*TMath::Pi()*rapidityMeas[i]);
              cout << "loaded " << nameMeasGlobalPi0[i] << "acceptance and efficiency" << endl;
          }
      }
      
      //______________________________ Eta meson invariant cross section
        histoEtaInvCrossSection[i]        = (TH1D*)directoryEta[i]->Get("InvCrossSectionEta");
        graphEtaInvCrossSectionStat[i]    = (TGraphAsymmErrors*)directoryEta[i]->Get("graphInvCrossSectionEta");
        graphEtaInvCrossSectionSys[i]     = (TGraphAsymmErrors*)directoryEta[i]->Get("InvCrossSectionEtaSys");
        graphEtaToOmegaPCMStat              = (TGraphAsymmErrors*)directoryEta[i]->Get("graphEtaToOmegaStatError");
        graphEtaToOmegaPCMSys               = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToOmegaSystError");
        if(!graphEtaInvCrossSectionSys[i] || !histoEtaInvCrossSection[i]){
          cout << "missing invariant cross section histograms of eta... returning!" << endl;
          return;
        } else {
          cout << "loaded " << nameMeasGlobal[i] << " eta invariant cross section" << endl;
        }
        
        //______________________________ Eta meson invariant mass peak pos and FWHM
        histoEtaMass[i]                     = (TH1D*)directoryEta[i]->Get("Eta_Mass_data_INT1");
        histoEtaFWHMMeV[i]                  = (TH1D*)directoryEta[i]->Get("Eta_Width_data_INT1");
        histoEtaTrueMass[i]                 = (TH1D*)directoryEta[i]->Get("Eta_Mass_MC_INT1");
        histoEtaTrueFWHMMeV[i]              = (TH1D*)directoryEta[i]->Get("Eta_Width_MC_INT1");
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
        histoEtaAcc[i]                      = (TH1D*)directoryEta[i]->Get("AcceptanceEta_INT1");
        histoEtaTrueEffPt[i]                = (TH1D*)directoryEta[i]->Get("EfficiencyEta_INT1");
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
    }

    // calculation of relative statistical and systematic uncertainties
    // FOR Omega:
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
            statErrorCollection[i]              = (TH1D*)histoOmegaInvCrossSection[i]->Clone(Form("statErr%sOmega",nameMeasGlobal[i].Data()));
            sysErrorCollection[i]               = (TGraphAsymmErrors*)graphOmegaInvCrossSectionSys[i]->Clone(Form("sysErr%sOmega",nameMeasGlobal[i].Data()));
            cout << "calculating Omega relative uncertainties for " << nameMeasGlobal[i] << endl;
        }
        // calculate relative errors for the available measurements
        if (statErrorCollection[i]) statErrorRelCollection[i] = CalculateRelErrUpTH1D( statErrorCollection[i], Form("relativeStatErrorOmega_%s", nameMeasGlobal[i].Data()));
        if (sysErrorCollection[i]) sysErrorRelCollection[i]   = CalculateRelErrUpAsymmGraph( sysErrorCollection[i], Form("relativeSysErrorOmega_%s", nameMeasGlobal[i].Data()));
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


    Double_t minX                               = 1;
    Double_t maxX                               = 39.;

    // **********************************************************************************************************************
    // ******************************************* Omega invariant cross section            ****************************************
    // **********************************************************************************************************************
    Double_t textSizeLabelsPixel             = 50;
    Double_t textSizeLabelsRel      = 50./1200;
    TCanvas* canvasCrossSectionOmega       = new TCanvas("canvasCrossSectionOmega", "", 0,0,1250,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasCrossSectionOmega,  0.14, 0.01, 0.01, 0.08);
    canvasCrossSectionOmega->SetLogy(1);
    canvasCrossSectionOmega->SetLogx(1);
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
        TH1D * histoBlack               = (TH1D*)histoOmegaInvCrossSection[exampleInteger]->Clone("histoBlack") ;
        histoBlack->SetLineColor(kBlack);
        histoBlack->SetMarkerStyle(21) ;
        histoBlack->SetMarkerColor(kBlack) ;
        histoBlack->SetMarkerSize(1.) ;

        TGraphAsymmErrors * graphGrey   = (TGraphAsymmErrors*)graphOmegaInvCrossSectionSys[exampleInteger]->Clone("graphGrey") ;
        graphGrey   ->SetFillColor(kGray);
        graphGrey   ->SetLineColor(kGray);
        graphGrey   ->SetFillStyle(1001);
        
        TGraphAsymmErrors * graphOmegaInvCrossSectionSysClone[10];
        TGraphAsymmErrors * graphOmegaInvCrossSectionStatClone[10];
        
        TLegend* legendCrossSectionOmega           = GetAndSetLegend2(0.18, 0.13, 0.43, 0.13+(3.5*textSizeLabelsRel),textSizeLabelsPixel);
        
        Int_t scalingFactor = pow(10,totalMeasAvail-1);
        // Int_t legendScalingFactor = totalMeasAvail-1;
        Int_t legendScalingFactor = 0;
        for(Int_t i=numbersofmeas;i>-1;i--){
          if(availableMeas[i]){
            graphOmegaInvCrossSectionSysClone[i]= (TGraphAsymmErrors*)graphOmegaInvCrossSectionSys[i]->Clone(Form("csSys%d",i)) ;
            graphOmegaInvCrossSectionStatClone[i]= (TGraphAsymmErrors*)graphOmegaInvCrossSectionStat[i]->Clone(Form("csStat%d",i)) ;
            // for (int j=0;j<graphOmegaInvCrossSectionSysClone[i]->GetN();j++){
            //     graphOmegaInvCrossSectionSysClone[i]->GetY()[j] *= scalingFactor;
            //     graphOmegaInvCrossSectionSysClone[i]->GetEYhigh()[j] *= scalingFactor;
            //     graphOmegaInvCrossSectionSysClone[i]->GetEYlow()[j] *= scalingFactor;
            // }
            // DrawGammaSetMarkerTGraphAsym(graphOmegaInvCrossSectionSysInterpolation[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDetMC[i] , colorDetMC[i], widthLinesBoxes, kTRUE);
            // graphOmegaInvCrossSectionSysInterpolation[i]     ->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphOmegaInvCrossSectionSysClone[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
            graphOmegaInvCrossSectionSysClone[i]     ->Draw("E2same");
            legendCrossSectionOmega->AddEntry(graphOmegaInvCrossSectionSysClone[i],Form("%s x10^{%d}",nameMeasGlobal[i].Data(),legendScalingFactor),"pf");

            // statistics
            // for (int j=0;j<graphOmegaInvCrossSectionStatClone[i]->GetN();j++){
            //     graphOmegaInvCrossSectionStatClone[i]->GetY()[j] *= scalingFactor;
            //     graphOmegaInvCrossSectionStatClone[i]->GetEYhigh()[j] *= scalingFactor;
            //     graphOmegaInvCrossSectionStatClone[i]->GetEYlow()[j] *= scalingFactor;
            // }
            // DrawGammaSetMarkerTGraph(graphOmegaInvCrossSectionStatInterpolation[i],  markerStyleDet[i], markerSizeDet[i]*0.55, colorDetMC[i] , colorDetMC[i]);
            // graphOmegaInvCrossSectionStatInterpolation[i]->Draw("p,same,z");
            DrawGammaSetMarkerTGraph(graphOmegaInvCrossSectionStatClone[i],  markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            graphOmegaInvCrossSectionStatClone[i]->Draw("p,same,z");
            scalingFactor/=10;
            // legendScalingFactor--;
          }
        }
        legendCrossSectionOmega->Draw();


        // TLegend* legendOmegaErr2 = GetAndSetLegend2(0.72, 0.72, 0.98, 0.72+(2*textSizeLabelsRel),0.85*textSizeLabelsPixel);
        // legendOmegaErr2->AddEntry(histoBlack, "stat. Err.","ple");
        // legendOmegaErr2->AddEntry(graphGrey,  "syst. Err.","f");
        // legendOmegaErr2->Draw();

        drawLatexAdd("ALICE",0.93,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd(collisionSystem5TeV,0.93,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        histoDummyCrossSection->Draw("sameaxis");
    canvasCrossSectionOmega->Update();
    canvasCrossSectionOmega->Print(Form("%s/Omega_InvariantCrossSectionMeas.%s",outputDir.Data(),suffix.Data()));
    
    // **********************************************************************************************************************
    // ******************************************* Omega invariant cross section            ****************************************
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
        // Int_t legendScalingFactor = totalMeasAvail-1;
        legendScalingFactor = 0;
        for(Int_t i=numbersofmeas;i>-1;i--){
          if(availableMeas[i]){
            graphEtaInvCrossSectionSysClone[i]= (TGraphAsymmErrors*)graphEtaInvCrossSectionSys[i]->Clone(Form("csSys%d",i)) ;
            graphEtaInvCrossSectionStatClone[i]= (TGraphAsymmErrors*)graphEtaInvCrossSectionStat[i]->Clone(Form("csStat%d",i)) ;
            // for (int j=0;j<graphEtaInvCrossSectionSysClone[i]->GetN();j++){
            //     graphEtaInvCrossSectionSysClone[i]->GetY()[j] *= scalingFactor;
            //     graphEtaInvCrossSectionSysClone[i]->GetEYhigh()[j] *= scalingFactor;
            //     graphEtaInvCrossSectionSysClone[i]->GetEYlow()[j] *= scalingFactor;
            // }
            DrawGammaSetMarkerTGraphAsym(graphEtaInvCrossSectionSysClone[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
            graphEtaInvCrossSectionSysClone[i]     ->Draw("E2same");
            legendCrossSectionEta->AddEntry(graphEtaInvCrossSectionSysClone[i],Form("%s x10^{%d}",nameMeasGlobal[i].Data(),legendScalingFactor),"pf");
              
            // statistics
            // for (int j=0;j<graphEtaInvCrossSectionStatClone[i]->GetN();j++){
            //     graphEtaInvCrossSectionStatClone[i]->GetY()[j] *= scalingFactor;
            //     graphEtaInvCrossSectionStatClone[i]->GetEYhigh()[j] *= scalingFactor;
            //     graphEtaInvCrossSectionStatClone[i]->GetEYlow()[j] *= scalingFactor;
            // }
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
        drawLatexAdd("#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}",0.93,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        histoDummyCrossSectionEta->Draw("sameaxis");
    canvasCrossSectionOmega->Update();
    canvasCrossSectionOmega->Print(Form("%s/Eta_InvariantCrossSectionMeas.%s",outputDir.Data(),suffix.Data()));
     
    // **********************************************************************************************************************
    // ******************************************* Mass and width for Omega            ****************************************
    // **********************************************************************************************************************
    cout << "plotting Omega mass" << endl;
    Double_t arrayBoundariesX1_4[2];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    textSizeLabelsPixel             = 35;
    ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasMassWidthOmega     = new TCanvas("canvasMassWidthOmega","",0,0,1350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthOmega,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthOmega               = new TPad("padWidthOmega", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthOmega, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthOmega->Draw();

    TPad* padMassOmega                = new TPad("padMassOmega", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassOmega, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassOmega->Draw();

    TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.32, 0.52, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
    padMassLegend1->SetFillStyle(0);
    padMassLegend1->Draw();

    padWidthOmega->cd();
    padWidthOmega->SetLogx();
    Double_t margin                 = relativeMarginsX[0]*2.7*1350;
    Double_t textsizeLabelsWidth    = 0;
    Double_t textsizeFacWidth       = 0;
    if (padWidthOmega->XtoPixel(padWidthOmega->GetX2()) < padWidthOmega->YtoPixel(padWidthOmega->GetY1())){
        textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthOmega->XtoPixel(padWidthOmega->GetX2()) ;
        textsizeFacWidth            = (Double_t)1./padWidthOmega->XtoPixel(padWidthOmega->GetX2()) ;
    } else {
        textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthOmega->YtoPixel(padWidthOmega->GetY1());
        textsizeFacWidth            = (Double_t)1./padWidthOmega->YtoPixel(padWidthOmega->GetY1());
    }

    TH2F * histo2DAllOmegaFWHM    = new TH2F("histo2DAllOmegaFWHM","histo2DAllOmegaFWHM", 20, minX, maxX ,1000., -30, 80);
    SetStyleHistoTH2ForGraphs(histo2DAllOmegaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
    histo2DAllOmegaFWHM->GetYaxis()->SetRangeUser(-1.,79.5);
    histo2DAllOmegaFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllOmegaFWHM->GetYaxis()->SetNdivisions(505);
    histo2DAllOmegaFWHM->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllOmegaFWHM->GetXaxis()->SetTickLength(0.05);
    histo2DAllOmegaFWHM->GetYaxis()->SetTickLength(0.026);
    histo2DAllOmegaFWHM->DrawCopy();

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoOmegaFWHMMeV[i] && histoOmegaTrueFWHMMeV[i] && availableMeas[i]){
            DrawGammaSetMarker(histoOmegaFWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            histoOmegaFWHMMeV[i]->Draw("p,same,e");
            DrawGammaSetMarker(histoOmegaTrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
            histoOmegaTrueFWHMMeV[i]->Draw("p,same,e");
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
    TLatex *labelMassOmega        = new TLatex(0.13,0.69,"#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}");
    SetStyleTLatex( labelMassOmega, textSizeLabelsPixel,4);
    labelMassOmega->SetTextFont(43);
    labelMassOmega->Draw();

    padMassOmega->cd();
    padMassOmega->SetLogx();

    Double_t textsizeLabelsMass         = 0;
    Double_t textsizeFacMass            = 0;
    if (padMassOmega->XtoPixel(padMassOmega->GetX2()) <padMassOmega->YtoPixel(padMassOmega->GetY1()) ){
        textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassOmega->XtoPixel(padMassOmega->GetX2()) ;
        textsizeFacMass                 = (Double_t)1./padMassOmega->XtoPixel(padMassOmega->GetX2()) ;
    } else {
        textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassOmega->YtoPixel(padMassOmega->GetY1());
        textsizeFacMass                 = (Double_t)1./padMassOmega->YtoPixel(padMassOmega->GetY1());
    }

    TH2F * histo2DAllOmegaMass            = new TH2F("histo2DAllOmegaMass","histo2DAllOmegaMass",20, minX, maxX, 1000., 730, 860);
    SetStyleHistoTH2ForGraphs(histo2DAllOmegaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
    histo2DAllOmegaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllOmegaMass->GetYaxis()->SetNdivisions(505);
    histo2DAllOmegaMass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllOmegaMass->GetXaxis()->SetTickLength(0.05);
    histo2DAllOmegaMass->GetXaxis()->SetLabelOffset(-0.015);
    histo2DAllOmegaMass->DrawCopy();

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoOmegaMass[i] && histoOmegaTrueMass[i] && availableMeas[i]){
            DrawGammaSetMarker(histoOmegaMass[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            histoOmegaMass[i]->Draw("p,same,e");
            DrawGammaSetMarker(histoOmegaTrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
            histoOmegaTrueMass[i]->Draw("p,same,e");
        }
    }

    DrawGammaLines(minX, maxX , mesonMassExpectOmega*1000., mesonMassExpectOmega*1000.,0.1, kGray);

    TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
    SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4);
    labelLegendBMass->SetTextFont(43);
    labelLegendBMass->Draw();

    //********************************** Defintion of the Legend **************************************************
    Double_t columnsLegendMass2[3]      = {0.,0.57,0.84};
    //   Double_t rowsLegendMass2[5] = {0.8,0.6,0.4,0.2,0.01};
    Double_t rowsLegendMass2[6] = {0.84,0.66,0.50,0.33,0.16,0.01};
    //           Double_t  rowsLegendMass2[7]= {0.84,0.66,0.50,0.33,0.01,0.16};
    //Double_t  rowsLegendMass2[9]= {0.84,0.66,0.51,0.50,0.331,0.33,0.01,0.16}; //setting for use without PHOS and PCM-PHOS peak positions
    //******************* Offsets ***********************
    Double_t offsetMarkerXMass2         = 0.1;
    Double_t offsetMarkerYMass2         = 0.1;
    //****************** Scale factors ******************
    Double_t scaleMarkerMass2           = 1.2;

    padMassLegend1->cd();
    //****************** first Column **************************************************
    TLatex *textMassPCM[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoOmegaMass[i] && histoOmegaTrueMass[i] && histoOmegaFWHMMeV[i] && histoOmegaTrueFWHMMeV[i] && availableMeas[i]){
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

    TMarker* markerPCMOmegaMass[10];
    TMarker* markerPCMOmegaMassMC[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoOmegaMass[i] && histoOmegaTrueMass[i]&& availableMeas[i]){
            markerPCMOmegaMass[i]             = CreateMarkerFromHisto(histoOmegaMass[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
            markerPCMOmegaMass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
            markerPCMOmegaMassMC[i]           = CreateMarkerFromHisto(histoOmegaTrueMass[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
            markerPCMOmegaMassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
        }
    }

    canvasMassWidthOmega->Update();
    canvasMassWidthOmega->Print(Form("%s/Omega_MassAndWidth.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ******************************************* Mass and width for Eta            ****************************************
    // **********************************************************************************************************************
    cout << "plotting eta mass" << endl;
    textSizeLabelsPixel             = 50;
    canvasMassWidthOmega->cd();
    padWidthOmega->Draw();
    padMassOmega->Draw();
    padMassLegend1->Draw();

    padWidthOmega->cd();
    padWidthOmega->SetLogx();

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
    TLatex *labelMassEta        = new TLatex(0.13,0.69,"#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}");
    SetStyleTLatex( labelMassEta, textSizeLabelsPixel,4);
    labelMassEta->SetTextFont(43);
    labelMassEta->Draw();

    padMassOmega->cd();
    padMassOmega->SetLogx();

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

    canvasMassWidthOmega->Update();
    canvasMassWidthOmega->Print(Form("%s/Eta_MassAndWidth.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for Omega single measurement  *********************************
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
        if(histoOmegaAccTimesEff[i] && availableMeas[i]){
            DrawGammaSetMarker(histoOmegaAccTimesEff[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            // Draw Pi0 measurment for comparison. For now use MC markers
            if(availableMeasPi0[i]) DrawGammaSetMarker(histoPi0AccTimesEff[i], markerStyleDetMC[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            histoOmegaAccTimesEff[i]->Draw("p,same,e");
            if(availableMeasPi0[i]) histoPi0AccTimesEff[i]->Draw("p,same,e");
        }
    }

    TLegend* legendEffiAccOmega           = GetAndSetLegend2(0.59, 0.13, 0.96, 0.13+(numbersofmeas*textSizeLabelsRel),textSizeLabelsPixel);
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoOmegaAccTimesEff[i] && availableMeas[i]){
            legendEffiAccOmega->AddEntry(histoOmegaAccTimesEff[i],nameMeasGlobal[i].Data(),"p");
            if(availableMeasPi0[i]) legendEffiAccOmega->AddEntry(histoPi0AccTimesEff[i],nameMeasGlobalPi0[i].Data(),"p");
        }
    }
    legendEffiAccOmega->Draw();

    drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(collisionSystem5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
    drawLatexAdd("#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}",0.15,0.82,textSizeLabelsRel,kFALSE);

    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Omega_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));
    
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
    TLegend* legendEffiAccOmegaEta           = GetAndSetLegend2(0.59, 0.13, 0.93, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoOmegaAccTimesEff[i] && availableMeas[i]){
            legendEffiAccOmegaEta->AddEntry(histoEtaAccTimesEff[i],nameMeasGlobal[i].Data(),"p");
        }
    }
    legendEffiAccEta->Draw();

    drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(collisionSystem5TeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
    drawLatexAdd("#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}",0.15,0.82,textSizeLabelsRel,kFALSE);

    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Eta_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));

}
