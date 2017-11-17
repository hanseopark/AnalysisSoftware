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

  void CombineMesonMeasurements900GeV(    TString fileNamePCM     = "/home/nschmidt/AnalysisResults/pp/900GeV/PCM/data_PCMResultsFullCorrection_PP_20170817.root",
                                          TString fileNamePHOS    = "/home/nschmidt/AnalysisResults/pp/900GeV/PHOS/PHOS_pp_pi0_900GeV_noBWcorr_K0Scorr_20111206.root",
                                          TString fileNameEMCAL   = "/home/nschmidt/AnalysisResults/pp/900GeV/EMCal/data_EMCAL-EMCALResultsFullCorrection_PP_20170420.root",
                                          TString fileNamePCMEMCAL= "/home/nschmidt/AnalysisResults/pp/900GeV/PCM-EMC/data_PCM-EMCALResultsFullCorrection_PP_20170420.root",
                                          TString fileInputCorrFactors = "/home/nschmidt/AnalysisResults/pp/7TeV/Comb/correlationInput/ComputeCorrelationFactors_pp900GeV/pp7TeV.root",
                                          TString suffix          = "pdf",
                                          TString isMC            = "",
                                          TString thesisPlots     = "",
                                          TString bWCorrection    ="",
                                          Int_t numbersofmeas     = 5
                                      ){

      TString date = ReturnDateString();

      gROOT->Reset();
      gROOT->SetStyle("Plain");

      StyleSettingsThesis();
      SetPlotStyle();

      gStyle->SetEndErrorSize(0);

      TString dateForOutput                       = ReturnDateStringForOutput();
      cout << dateForOutput.Data() << endl;

      //___________________________________ Declaration of files _____________________________________________
      TString collisionSystem900GeV               = "pp, #sqrt{#it{s}} = 900 GeV";

      if(thesisPlots.Contains("Y")) collisionSystem900GeV = "ALICE this thesis";

      TString fileNameTheory                      = "ExternalInput/Theory/TheoryCompilationPP.root";
      TString fileNameEtaToPi0                    = "ExternalInput/WorldDataPi0Eta.root";
      TString outputDir                           = Form("%s/%s/CombineMesonMeasurements900GeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
      cout << outputDir.Data() << endl;

      cout << "------------------------------------------" << endl;
      cout << "input files: " << endl;
      cout << fileNamePCM.Data() << endl;
      cout << fileNamePHOS.Data() << endl;
      cout << fileNameEMCAL.Data() << endl;
      cout << fileNamePCMEMCAL.Data() << endl;
      cout << "------------------------------------------" << endl;

      gSystem->Exec("mkdir -p "+outputDir);
      gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
      gSystem->Exec(Form("cp %s %s/InputPCMEMC.root", fileNamePCMEMCAL.Data(), outputDir.Data()));
      gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
      gSystem->Exec(Form("cp %s %s/InputEMCAL.root", fileNameEMCAL.Data(), outputDir.Data()));
      gSystem->Exec(Form("cp %s %s/Theory.root", fileNameTheory.Data(), outputDir.Data()));

      fstream fLog;
      fLog.open(Form("%s/CombineMeson900GeV.log",outputDir.Data()), ios::out);
      fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      fLog << dateForOutput.Data() << endl;

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

      Double_t xSection900GeVINEL                 = 52.5*1e-3;
      Double_t xSection900GeV                     = 47.78*1e-3;
      Double_t recalcBarn                         = 1e12; //NLO in pbarn!!!!

      Width_t widthLinesBoxes                     = 1.4;
      Width_t widthCommonFit                      = 2;

      // Definition of colors, styles and markers sizes
      Color_t colorComb                           = kBlue+2;
      Style_t markerStyleComb                     = 20;
      Size_t  markerSizeComb                      = 2;

      Color_t  colorCGC                           = kCyan-8;
      Color_t  colorNLO                           = kAzure-4;

      Color_t colorTrigg      [10]                = {kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2};

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
      TString nameMeasGlobalWriteToFile[11]       = {"PCM",
                                                     "PHOS",
                                                     "EMCAL",
                                                     "PCMPHOS",
                                                     "PCMEMCAL",
                                                     "PCMDalitz",
                                                     "PHOSDalitz",
                                                     "EMCALDalitz",
                                                     "EMCAL",
                                                     "EMCALmerged",
                                                     "PCMOtherDataset"};
      TString  nameSecPi0SourceRead[4]            = {"K0S", "K0L", "Lambda", "Rest"};
      TString  nameSecPi0SourceLabel[4]           = {"K^{0}_{s}", "K^{0}_{l}", "#Lambda", "had. int."};
      Double_t maxSecCorr[4]                      = { 0.06, 0.007, 0.00005, 0.04};

      Bool_t haveEffSecCorr[4][11]                    = { {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE},
                                                          {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE},
                                                          {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE},
                                                          {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE} };
      TGraphAsymmErrors* graphPi0EffSecCorrFromX[4][11]   = { { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                              { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                              { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                              { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} };

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
      inputFile[4]                            = new TFile(fileNamePCMEMCAL.Data());

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
            directoryPi0[i]                     = (TDirectory*)inputFile[i]->Get("Pi0900GeV");
            directoryEta[i]                     = (TDirectory*)inputFile[i]->Get("Eta900GeV");
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

      TH1D* histoPi0InvXSection[11]                           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
      TGraphAsymmErrors* graphPi0InvXSectionSys[11]           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
      TGraphAsymmErrors* graphPi0InvXSectionStat[11]          = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

      TH1D* histoEtaMass[11]                                  = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
      TH1D* histoEtaFWHMMeV[11]                               = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
      TH1D* histoEtaTrueMass[11]                              = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
      TH1D* histoEtaTrueFWHMMeV[11]                           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
      TH1D* histoEtaAcc[11]                                   = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
      TH1D* histoEtaTrueEffPt[11]                             = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
      TH1D* histoEtaAccTimesEff[11]                           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

      TH1D* histoEtaInvXSection[11]                           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
      TGraphAsymmErrors* graphEtaInvXSectionSys[11]           = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
      TGraphAsymmErrors* graphEtaInvXSectionStat[11]          = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

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

      Double_t minPtPi0 = 0.33;
      Double_t maxPtPi0 = 12.0;
      Double_t minXSectionPi0 = 2e2;
      Double_t maxXSectionPi0 = 9e12;

      Double_t minPtEta = 0.7;
      Double_t minPtEtaToPi0 = 0.33;
      Double_t maxPtEta = 5.1;
      Double_t minXSectionEta = 1e5;
      Double_t maxXSectionEta = 9e9;

      Bool_t doOutput = kTRUE;

      // *******************************************************************************************************
      // ************************** define binnings ************************************************************
      // *******************************************************************************************************

       Int_t nBinsPi0 = 15;
       Double_t xPtLimits[51]                      =  {0.0, 0.4, 0.6, 0.8, 1.0,
                                                       1.2, 1.4, 1.6, 2.0, 2.5,
                                                       3.0, 3.5, 4.0, 5.0, 7.0,
                                                       10.
                                                      };
       Int_t nBinsEta = 4;
       Double_t xPtLimitsEta[5]                    =  {0.0, 0.9, 1.8, 3.0, 4.0
                                                      };

      // *******************************************************************************************************
      // ************************** read input files ***********************************************************
      // *******************************************************************************************************

      for (Int_t i = 0; i < 11; i++){
        cout << "reading from " << nameMeasGlobal[i].Data() << " measurement..." << endl;

        if(i==1){ // PHOS
          directoryPi0[i] = (TDirectory*)inputFile[0]->Get("Pi0900GeV");
          histoPi0InvXSection[i]          = (TH1D*)inputFile[i]->Get("hPi0900GeVStat");
          histoPi0InvXSection[i]->Sumw2();
          histoPi0InvXSection[i]->Scale(xSection900GeV*recalcBarn);
          graphPi0InvXSectionStat[i]      = new TGraphAsymmErrors(histoPi0InvXSection[i]);
          TH1D* histoPi0InvXSectionSys    = (TH1D*)inputFile[i]->Get("hPi0900GeVSys");
          histoPi0InvXSectionSys->Sumw2();
          histoPi0InvXSectionSys->Scale(xSection900GeV*recalcBarn);
          graphPi0InvXSectionSys[i]       = new TGraphAsymmErrors(histoPi0InvXSectionSys);

          cout << nameMeasGlobal[i].Data() << " pi0 stat:" << graphPi0InvXSectionStat[i] << endl;
          if(doOutput) graphPi0InvXSectionStat[i]->Print();
          cout << nameMeasGlobal[i].Data() << " pi0 sys:" << graphPi0InvXSectionSys[i] << endl;
          if(doOutput) graphPi0InvXSectionSys[i]->Print();

          continue;
        }

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

          histoPi0AccTimesEff[i]            = (TH1D*)histoPi0TrueEffPt[i]->Clone(Form("histoPi0AccTimesEff%s",nameMeasGlobal[i].Data()));
          histoPi0AccTimesEff[i]            ->Multiply(histoPi0Acc[i]);
          histoPi0AccTimesEff[i]            ->Scale(2*TMath::Pi()*rapidityMeas[i]);

          // load cross section systematics and datapoints
          histoPi0InvXSection[i]          = (TH1D*)directoryPi0[i]->Get("InvCrossSectionPi0");
          graphPi0InvXSectionStat[i]      = new TGraphAsymmErrors(histoPi0InvXSection[i]);
            while (graphPi0InvXSectionStat[i]->GetY()[0] <= 1E-50 ) graphPi0InvXSectionStat[i]->RemovePoint(0);
            while (graphPi0InvXSectionStat[i]->GetY()[graphPi0InvXSectionStat[i]->GetN()-1] <= 1E-50 ) graphPi0InvXSectionStat[i]->RemovePoint(graphPi0InvXSectionStat[i]->GetN()-1);
          graphPi0InvXSectionSys[i]       = (TGraphAsymmErrors*)directoryPi0[i]->Get("InvCrossSectionPi0Sys");
            while (graphPi0InvXSectionSys[i]->GetY()[0] <= 1E-50 ) graphPi0InvXSectionSys[i]->RemovePoint(0);
            while (graphPi0InvXSectionSys[i]->GetY()[graphPi0InvXSectionSys[i]->GetN()-1] <= 1E-50 ) graphPi0InvXSectionSys[i]->RemovePoint(graphPi0InvXSectionSys[i]->GetN()-1);

          cout << nameMeasGlobal[i].Data() << " pi0 stat:" << graphPi0InvXSectionStat[i] << endl;
          if(doOutput) graphPi0InvXSectionStat[i]->Print();
          cout << nameMeasGlobal[i].Data() << " pi0 sys:" << graphPi0InvXSectionSys[i] << endl;
          if(doOutput) graphPi0InvXSectionSys[i]->Print();

          // load invariant mass example bins
          histoPi0InvMassSig[i]               = (TH1D*)directoryPi0[i]->Get(Form("InvMassSig_PtBin%s",strInvMassBin[i].Data()));
          histoPi0InvMassSigPlusBG[i]         = (TH1D*)directoryPi0[i]->Get(Form("InvMassSigPlusBG_PtBin%s",strInvMassBin[i].Data()));
          histoPi0InvMassBG[i]                = (TH1D*)directoryPi0[i]->Get(Form("InvMassBG_PtBin%s",strInvMassBin[i].Data()));
          fitPi0InvMassSig[i]                 = (TF1*)directoryPi0[i]->Get(Form("FitInvMassSig_PtBin%s",strInvMassBin[i].Data()));

          if (histoPi0InvMassSig[i] && histoPi0InvMassSigPlusBG[i] && histoPi0InvMassBG[i] && fitPi0InvMassSig[i]){
            haveAllPi0InvMass[i]              = kTRUE;
          }

          for (Int_t k = 0; k < 4; k++){
              graphPi0EffSecCorrFromX[k][i]           = (TGraphAsymmErrors*)directoryPi0[i]->Get(Form("EffectiveSecondaryPi0CorrFrom%s",nameSecPi0SourceRead[k].Data()));
              if (graphPi0EffSecCorrFromX[k][i]){
                  cout << nameSecPi0SourceRead[k].Data() << endl;
                  if(doOutput) graphPi0EffSecCorrFromX[k][i]->Print();
                  Int_t nAboveZero                    = 0;
                  for (Int_t j = 0; j< graphPi0EffSecCorrFromX[k][i]->GetN(); j++){
                      if(graphPi0EffSecCorrFromX[k][i]->GetY()[j] > 0) nAboveZero++;
                  }
                  if (nAboveZero>0){
                      haveEffSecCorr[k][i]            = kTRUE;
                  } else {
                      graphPi0EffSecCorrFromX[k][i]   = NULL;
                  }
              }
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

          histoEtaAccTimesEff[i]            = (TH1D*)histoEtaTrueEffPt[i]->Clone(Form("histoEtaAccTimesEff%s",nameMeasGlobal[i].Data()));
          histoEtaAccTimesEff[i]            ->Multiply(histoEtaAcc[i]);
          histoEtaAccTimesEff[i]            ->Scale(2*TMath::Pi()*rapidityMeas[i]);

          // load cross section systematics and datapoints
          histoEtaInvXSection[i]          = new TH1D("InvCrossSectionEta","InvCrossSectionEta",4,xPtLimitsEta);
          graphEtaInvXSectionStat[i]      = (TGraphAsymmErrors*)directoryEta[i]->Get("graphInvCrossSectionEta");
          for(Int_t iB=0; iB < graphEtaInvXSectionStat[i]->GetN(); iB++) {
            histoEtaInvXSection[i]->SetBinContent(histoEtaInvXSection[i]->FindBin(graphEtaInvXSectionStat[i]->GetX()[iB]),graphEtaInvXSectionStat[i]->GetY()[iB]);
            histoEtaInvXSection[i]->SetBinError(histoEtaInvXSection[i]->FindBin(graphEtaInvXSectionStat[i]->GetX()[iB]),graphEtaInvXSectionStat[i]->GetErrorYhigh(iB));
          }

          graphEtaInvXSectionSys[i]       = (TGraphAsymmErrors*)directoryEta[i]->Get("InvCrossSectionEtaSys");
            for(Int_t iB=0;iB<histoEtaInvXSection[i]->GetNbinsX();iB++){if(histoEtaInvXSection[i]->GetBinContent(iB) <= 1E-50) histoEtaInvXSection[i]->SetBinContent(iB,0.);}
            while (graphEtaInvXSectionStat[i]->GetY()[0] <= 1E-50 ) graphEtaInvXSectionStat[i]->RemovePoint(0);
            while (graphEtaInvXSectionStat[i]->GetY()[graphEtaInvXSectionStat[i]->GetN()-1] <= 1E-50 ) graphEtaInvXSectionStat[i]->RemovePoint(graphEtaInvXSectionStat[i]->GetN()-1);
            while (graphEtaInvXSectionSys[i]->GetY()[0] <= 1E-50 ) graphEtaInvXSectionSys[i]->RemovePoint(0);
            while (graphEtaInvXSectionSys[i]->GetY()[graphEtaInvXSectionSys[i]->GetN()-1] <= 1E-50 ) graphEtaInvXSectionSys[i]->RemovePoint(graphEtaInvXSectionSys[i]->GetN()-1);

          cout << nameMeasGlobal[i].Data() << " eta stat:" << graphEtaInvXSectionStat[i] << endl;
          if(doOutput) graphEtaInvXSectionStat[i]->Print();
          cout << nameMeasGlobal[i].Data() << " eta sys:" << graphEtaInvXSectionSys[i] << endl;
          if(doOutput) graphEtaInvXSectionSys[i]->Print();


          histoEtaToPi0Stat[i]               = new TH1D("EtaToPi0StatErr","EtaToPi0StatErr",4,xPtLimitsEta);
          graphEtaToPi0Stat[i]               = (TGraphAsymmErrors*)directoryEta[i]->Get("graphEtaToPi0StatError");
          for(Int_t iB=0; iB < graphEtaToPi0Stat[i]->GetN(); iB++) {
            histoEtaToPi0Stat[i]->SetBinContent(histoEtaToPi0Stat[i]->FindBin(graphEtaToPi0Stat[i]->GetX()[iB]),graphEtaToPi0Stat[i]->GetY()[iB]);
            histoEtaToPi0Stat[i]->SetBinError(histoEtaToPi0Stat[i]->FindBin(graphEtaToPi0Stat[i]->GetX()[iB]),graphEtaToPi0Stat[i]->GetErrorYhigh(iB));
          }
          graphEtaToPi0Sys[i]                = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0SystError");

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

      // *******************************************************************************************************
      // ************************** Loading theory calculations ************************************************
      // *******************************************************************************************************
      TFile* fileEtaToPi                              = new TFile(fileNameEtaToPi0.Data());
      TFile* fileTheoryCompilation                    = new TFile(fileNameTheory.Data());

          // Pythia8 Monash2013:
          TH1F* histoPythia8InvXSection                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi0900GeV");
          histoPythia8InvXSection->GetXaxis()->SetRangeUser(minPtPi0,maxPtPi0);
          TH1F* histoPythia8InvXSectionEta                    = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta900GeV");
          histoPythia8InvXSectionEta->GetXaxis()->SetRangeUser(minPtEta,maxPtEta);
          TGraphErrors* graphPythia8InvXSection               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi0900GeV"));
          while(graphPythia8InvXSection->GetX()[0] < minPtPi0) graphPythia8InvXSection->RemovePoint(0);
          while(graphPythia8InvXSection->GetX()[graphPythia8InvXSection->GetN()-1] > 12.) graphPythia8InvXSection->RemovePoint(graphPythia8InvXSection->GetN()-1);
          TGraphErrors* graphPythia8InvXSectionEta            = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta900GeV"));
          while(graphPythia8InvXSectionEta->GetX()[0] < minPtEta) graphPythia8InvXSectionEta->RemovePoint(0);
          while(graphPythia8InvXSectionEta->GetX()[graphPythia8InvXSectionEta->GetN()-1] > 5.) graphPythia8InvXSectionEta->RemovePoint(graphPythia8InvXSectionEta->GetN()-1);
          TH1F* histoPythia8EtaToPi0                          = (TH1F*) histoPythia8InvXSectionEta->Clone("Pythia8EtaToPi0");
          histoPythia8EtaToPi0->Divide(histoPythia8InvXSection);
          histoPythia8EtaToPi0->GetXaxis()->SetRangeUser(minPtEta,maxPtEta);
          TGraphErrors* graphPythia8EtaToPi0                  = new TGraphErrors(histoPythia8EtaToPi0);
          while(graphPythia8EtaToPi0->GetX()[0] < 0.4) graphPythia8EtaToPi0->RemovePoint(0);
          while(graphPythia8EtaToPi0->GetX()[graphPythia8EtaToPi0->GetN()-1] > 5.) graphPythia8EtaToPi0->RemovePoint(graphPythia8EtaToPi0->GetN()-1);

          // *******************************************************************************************************
          // NLO calc
          TGraphAsymmErrors* graphPi0DSS07                    = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcDSS07InvSecPi0900GeV");
          while (graphPi0DSS07->GetX()[graphPi0DSS07->GetN()-1] > 12. ) graphPi0DSS07->RemovePoint(graphPi0DSS07->GetN()-1);

          TGraphAsymmErrors* graphEtaAESSS                    = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcAESSSInvSecEta900GeV");
          while (graphEtaAESSS->GetX()[graphEtaAESSS->GetN()-1] > 5. ) graphEtaAESSS->RemovePoint(graphEtaAESSS->GetN()-1);

          TGraphAsymmErrors* graphEtaToPi02760GeV             = (TGraphAsymmErrors*) fileEtaToPi->Get("Alice2760GeV");
          ProduceGraphAsymmWithoutXErrors(graphEtaToPi02760GeV);
          TGraphAsymmErrors* graphEtaToPi08000GeV             = (TGraphAsymmErrors*) fileEtaToPi->Get("Alice8TeV");
          ProduceGraphAsymmWithoutXErrors(graphEtaToPi08000GeV);

          // *******************************************************************************************************

          TGraph* graphNLOCalcPi0MuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf900GeV");
          TGraph* graphNLOCalcPi0MuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne900GeV");
          TGraph* graphNLOCalcPi0MuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo900GeV");
          TGraph* graphNLOCalcEtaMuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf900GeV");
          TGraph* graphNLOCalcEtaMuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne900GeV");
          TGraph* graphNLOCalcEtaMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo900GeV");

          while (graphNLOCalcEtaMuHalf->GetX()[graphNLOCalcEtaMuHalf->GetN()-1] > 5. )
              graphNLOCalcEtaMuHalf->RemovePoint(graphNLOCalcEtaMuHalf->GetN()-1);
          while (graphNLOCalcEtaMuOne->GetX()[graphNLOCalcEtaMuOne->GetN()-1] > 5. )
              graphNLOCalcEtaMuOne->RemovePoint(graphNLOCalcEtaMuOne->GetN()-1);
          while (graphNLOCalcEtaMuTwo->GetX()[graphNLOCalcEtaMuTwo->GetN()-1] > 5. )
              graphNLOCalcEtaMuTwo->RemovePoint(graphNLOCalcEtaMuTwo->GetN()-1);

          TGraphAsymmErrors* graphNLOEtaToPi0                = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0900GeV_AESSS_DSS07");
          while (graphNLOEtaToPi0->GetX()[graphNLOEtaToPi0->GetN()-1] > 6. ) graphNLOEtaToPi0->RemovePoint(graphNLOEtaToPi0->GetN()-1);

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
          if(directoryPi0[i]) statErrorCollection[i]                  = (TH1D*)histoPi0InvXSection[i]->Clone(Form("statErr%sPi0",nameMeasGlobal[i].Data()));
          if(directoryEta[i]) statErrorCollectionEta[i]               = (TH1D*)histoEtaInvXSection[i]->Clone(Form("statErr%sEta",nameMeasGlobal[i].Data()));
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
          if(directoryPi0[i]) sysErrorCollection[i]                   = (TGraphAsymmErrors*)graphPi0InvXSectionSys[i]->Clone(Form("sysErr%sPi0",nameMeasGlobal[i].Data()));
          if(directoryEta[i]) sysErrorCollectionEta[i]                = (TGraphAsymmErrors*)graphEtaInvXSectionSys[i]->Clone(Form("sysErr%sEta",nameMeasGlobal[i].Data()));
          if(directoryEta[i]) sysErrorCollectionEtaToPi0[i]           = (TGraphAsymmErrors*)graphEtaToPi0Sys[i]->Clone(Form("sysErr%sEtaToPi0",nameMeasGlobal[i].Data()));
      }



      // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match
      Int_t offSets[11]                           =  {0,   2,    5,      2,        3,           0,   0,     0,      0,        0, 0};
      Int_t offSetsSys[11]                        =  {1,   2,    6,      3,        4,           0,   0,     0,      2,        0, 0};
      // pcm, phos, emcal, pcmphos, pcmemcal
      //  0    1      2       3         4

      Int_t offSetPi0Shifting[11]     = { 0,  0,  0,  0,  0,
                                          0,  0,  0,  0,  0,
                                          0 };
      Int_t nComBinsPi0Shifting[11]   = { 0,  0,  0,  0,  0,
                                          0,  0,  0,  0,  0,
                                          0 };

      //                                            PCM,PHOS,EMC,PCMPHOS,PCMEMC,         EMC
      Int_t offSetsEta[11]                        =  {0,    0,  0,     0,      0, 0,0,0,   0,0,0};
      Int_t offSetsSysEta[11]                     =  {1,    0,  0,     0,      0, 0,0,0,   0,0,0};

      Int_t offSetEtaShifting[11]     = { 0,  0,  0,  0,  0,
                                          0,  0,  0,  0,  0,
                                          0 };
      Int_t nComBinsEtaShifting[11]   = { 0,  0,  0,  0,  0,
                                          0,  0,  0,  0,  0,
                                          0 };

      //                                            PCM,PHOS,EMC,PCMPHOS,PCMEMC,         EMC
      Int_t offSetsEtaToPi0[11]                   =  {0,    0,  0,     0,      0, 0,0,0,   0,0,0};
      Int_t offSetsSysEtaToPi0[11]                =  {1,    0,  0,     0,      0, 0,0,0,   0,0,0};

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

      TGraphAsymmErrors* graphCombPi0InvXSectionStat= NULL;
      TGraphAsymmErrors* graphCombPi0InvXSectionSys = NULL;
      TGraphAsymmErrors* graphCombPi0InvXSectionTot = CombinePtPointsSpectraFullCorrMat( statErrorCollection, sysErrorCollection,
                                                                                         xPtLimits, nBinsPi0,
                                                                                         offSets, offSetsSys,
                                                                                         graphCombPi0InvXSectionStat, graphCombPi0InvXSectionSys,
                                                                                         fileNameOutputWeighting,"900GeV", "Pi0", kTRUE,
                                                                                         0x0, fileInputCorrFactors
                                                                                        );

      //return;
      if(doOutput) graphCombPi0InvXSectionStat->Print();

      TGraphAsymmErrors* graphCombEtaInvXSectionStat= NULL;
      TGraphAsymmErrors* graphCombEtaInvXSectionSys = NULL;
      TGraphAsymmErrors* graphCombEtaInvXSectionTot = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEta, sysErrorCollectionEta,
                                                                                         xPtLimitsEta, nBinsEta-1,
                                                                                         offSetsEta, offSetsSysEta,
                                                                                         graphCombEtaInvXSectionStat, graphCombEtaInvXSectionSys,
                                                                                         fileNameOutputWeightingEta,"900GeV", "Eta", kTRUE,
                                                                                         0x0, fileInputCorrFactors
                                                                                       );
      //return;
      if(doOutput) graphCombEtaInvXSectionStat->Print();


      TGraphAsymmErrors* graphCombEtaToPi0Stat= NULL;
      TGraphAsymmErrors* graphCombEtaToPi0Sys = NULL;
      TGraphAsymmErrors* graphCombEtaToPi0Tot = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEtaToPi0, sysErrorCollectionEtaToPi0,
                                                                                   xPtLimitsEta, nBinsEta-1,
                                                                                   offSetsEtaToPi0, offSetsSysEtaToPi0,
                                                                                   graphCombEtaToPi0Stat, graphCombEtaToPi0Sys,
                                                                                   fileNameOutputWeightingEtaToPi0,"900GeV", "EtaToPi0", kTRUE,
                                                                                   0x0, fileInputCorrFactors
                                                                                 );
      //return;
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
          TLatex *labelWeightsEnergy              = new TLatex(0.7,0.20,collisionSystem900GeV.Data());
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

          TLatex *labelRelSysErrEnergy            = new TLatex(0.15,0.89,collisionSystem900GeV.Data());
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

      TH2F * histo2DRelStatErr = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,minPtPi0,maxPtPi0,1000,0,87.5);
      SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
      histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
      histo2DRelStatErr->GetXaxis()->SetNoExponent();
      histo2DRelStatErr->GetYaxis()->SetRangeUser(0,87.5);
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

          TLatex *labelRelStatErrEnergy           = new TLatex(0.75,0.89,collisionSystem900GeV.Data());
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

      TGraphAsymmErrors* graphCombPi0InvXSectionRelStat     = CalculateRelErrUpAsymmGraph( graphCombPi0InvXSectionStat, "graphCombPi0InvXSectionRelStat");
      TGraphAsymmErrors* graphCombPi0InvXSectionRelSys      = CalculateRelErrUpAsymmGraph( graphCombPi0InvXSectionSys, "graphCombPi0InvXSectionRelSys");
      TGraphAsymmErrors* graphCombPi0InvXSectionRelTot      = CalculateRelErrUpAsymmGraph( graphCombPi0InvXSectionTot, "graphCombPi0InvXSectionRelTot");

      TCanvas* canvasRelTotErr                    = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
      DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
      canvasRelTotErr->SetLogx();

      TH2F * histo2DRelTotErr;
      histo2DRelTotErr                            = new TH2F("histo2DRelTotErr","histo2DRelTotErr",11000,minPtPi0,maxPtPi0,1000,0,87.5);
      SetStyleHistoTH2ForGraphs(histo2DRelTotErr, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
      histo2DRelTotErr->GetXaxis()->SetMoreLogLabels();
      histo2DRelTotErr->GetXaxis()->SetNoExponent();

      histo2DRelTotErr->GetYaxis()->SetRangeUser(0,87.5);
      histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
      histo2DRelTotErr->Draw("copy");

      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
      graphCombPi0InvXSectionRelTot->Draw("p,same,e1");
      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
      graphCombPi0InvXSectionRelStat->Draw("l,x0,same,e1");
      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
      graphCombPi0InvXSectionRelSys->SetLineStyle(7);
      graphCombPi0InvXSectionRelSys->Draw("l,x0,same,e1");

      TLegend* legendRelTotErr2 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
      legendRelTotErr2->AddEntry(graphCombPi0InvXSectionRelTot,"tot","p");
      legendRelTotErr2->AddEntry(graphCombPi0InvXSectionRelStat,"stat","l");
      legendRelTotErr2->AddEntry(graphCombPi0InvXSectionRelSys,"sys","l");
      legendRelTotErr2->Draw();

      TLatex *labelRelTotErrEnergy            = new TLatex(0.75,0.89,collisionSystem900GeV.Data());
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

      nMeasSet = 1;

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
          TLatex *labelWeightsEnergyEta              = new TLatex(0.7,0.20,collisionSystem900GeV.Data());
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

          TLatex *labelRelSysErrEnergyEta            = new TLatex(0.15,0.89,collisionSystem900GeV.Data());
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
      histo2DRelStatErrEta->GetYaxis()->SetRangeUser(0,65.5);
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

          TLatex *labelRelStatErrEnergyEta           = new TLatex(0.75,0.89,collisionSystem900GeV.Data());
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

      TGraphAsymmErrors* graphCombEtaInvXSectionRelStat     = CalculateRelErrUpAsymmGraph( graphCombEtaInvXSectionStat, "graphCombEtaInvXSectionRelStat");
      TGraphAsymmErrors* graphCombEtaInvXSectionRelSys      = CalculateRelErrUpAsymmGraph( graphCombEtaInvXSectionSys, "graphCombEtaInvXSectionRelSys");
      TGraphAsymmErrors* graphCombEtaInvXSectionRelTot      = CalculateRelErrUpAsymmGraph( graphCombEtaInvXSectionTot, "graphCombEtaInvXSectionRelTot");

      TCanvas* canvasRelTotErrEta                    = new TCanvas("canvasRelTotErrEta","",200,10,1350,900);  // gives the page size
      DrawGammaCanvasSettings( canvasRelTotErrEta, 0.08, 0.02, 0.035, 0.09);
      canvasRelTotErrEta->SetLogx();

      TH2F * histo2DRelTotErrEta = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,minPtEta,maxPtEta,1000,0,77.5);
      SetStyleHistoTH2ForGraphs(histo2DRelTotErrEta, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
      histo2DRelTotErrEta->GetXaxis()->SetMoreLogLabels();
      histo2DRelTotErrEta->GetXaxis()->SetNoExponent();

      histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,77.5);
      histo2DRelTotErrEta->GetYaxis()->SetTitle("Err (%)");
      histo2DRelTotErrEta->Draw("copy");

      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
      graphCombEtaInvXSectionRelTot->Draw("p,same,e1");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
      graphCombEtaInvXSectionRelStat->Draw("l,x0,same,e1");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
      graphCombEtaInvXSectionRelSys->SetLineStyle(7);
      graphCombEtaInvXSectionRelSys->Draw("l,x0,same,e1");

      TLegend* legendRelTotErr2Eta = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
      legendRelTotErr2Eta->AddEntry(graphCombEtaInvXSectionRelTot,"tot","p");
      legendRelTotErr2Eta->AddEntry(graphCombEtaInvXSectionRelStat,"stat","l");
      legendRelTotErr2Eta->AddEntry(graphCombEtaInvXSectionRelSys,"sys","l");
      legendRelTotErr2Eta->Draw();

      TLatex *labelRelTotErrEnergyEta            = new TLatex(0.75,0.89,collisionSystem900GeV.Data());
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
          TLatex *labelWeightsEnergyEtaToPi0              = new TLatex(0.7,0.20,collisionSystem900GeV.Data());
          SetStyleTLatex( labelWeightsEnergyEtaToPi0, 0.85*textSizeLabelsPixel,4);
          labelWeightsEnergyEtaToPi0->SetTextFont(43);
          labelWeightsEnergyEtaToPi0->Draw();
          TLatex *labelWeightsEtaToPi0                 = new TLatex(0.7,0.16,"#eta/#pi^{0} #rightarrow #gamma#gamma");
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

          TLatex *labelRelSysErrEnergyEtaToPi0            = new TLatex(0.15,0.89,collisionSystem900GeV.Data());
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
      histo2DRelStatErrEtaToPi0->GetYaxis()->SetRangeUser(0,75.5);
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

          TLatex *labelRelStatErrEnergyEtaToPi0           = new TLatex(0.75,0.89,collisionSystem900GeV.Data());
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

      TH2F * histo2DRelTotErrEtaToPi0 = new TH2F("histo2DRelTotErrEtaToPi0","histo2DRelTotErrEtaToPi0",11000,minPtEta,maxPtEta,1000,0,87.5);
      SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
      histo2DRelTotErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
      histo2DRelTotErrEtaToPi0->GetXaxis()->SetNoExponent();

      histo2DRelTotErrEtaToPi0->GetYaxis()->SetRangeUser(0,77.5);
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

      TLatex *labelRelTotErrEnergyEtaToPi0            = new TLatex(0.75,0.89,collisionSystem900GeV.Data());
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
      TGraphAsymmErrors* graphCombPi0InvXSectionTotUnShifted    = (TGraphAsymmErrors*)graphCombPi0InvXSectionTot->Clone("Unshifted");
      TGraphAsymmErrors* graphCombPi0InvXSectionStatUnShifted   = (TGraphAsymmErrors*)graphCombPi0InvXSectionStat->Clone("UnshiftedStat");
      TGraphAsymmErrors* graphCombPi0InvXSectionSysUnShifted    = (TGraphAsymmErrors*)graphCombPi0InvXSectionSys->Clone("UnshiftedSys");

      TGraphAsymmErrors* graphPi0InvXSectionStatUnShifted[11];
      TGraphAsymmErrors* graphPi0InvXSectionSysUnShifted[11];
      for (Int_t i = 0; i < 11; i++){
        if(directoryPi0[i]){
          graphPi0InvXSectionStatUnShifted[i] = (TGraphAsymmErrors*)graphPi0InvXSectionStat[i]->Clone(Form("UnshiftedStat%s",nameMeasGlobal[i].Data()));
          graphPi0InvXSectionSysUnShifted[i]  = (TGraphAsymmErrors*)graphPi0InvXSectionSys[i] ->Clone(Form("UnshiftedSys%s",nameMeasGlobal[i].Data()));
        }
      }

      // Calculating binshifts
      Double_t paramGraph[3]                      = {1.0e11, 8., 0.13};
      TF1* fitInvXSectionPi0              = FitObject("l","fitInvXSectionPi0","Pi0",graphCombPi0InvXSectionTot,0.4,10.,paramGraph,"QNRMEX0+");

      if(bWCorrection.Contains("X")){
          TF1* fitTsallisPi0PtMult        = FitObject("tmpt","TsallisMultWithPtPi0900GeV","Pi0");
          fitTsallisPi0PtMult->SetParameters(paramGraph[0],paramGraph[1], paramGraph[2]) ; // standard parameter optimize if necessary

          graphCombPi0InvXSectionTot      = ApplyXshift(graphCombPi0InvXSectionTot, fitTsallisPi0PtMult, "Pi0", kTRUE);

          cout << "comb" << endl;
          graphCombPi0InvXSectionStat     = ApplyXshiftIndividualSpectra (graphCombPi0InvXSectionTot,
                                                                          graphCombPi0InvXSectionStat,
                                                                          fitTsallisPi0PtMult,
                                                                          0, graphCombPi0InvXSectionStat->GetN());
          graphCombPi0InvXSectionSys      = ApplyXshiftIndividualSpectra (graphCombPi0InvXSectionTot,
                                                                          graphCombPi0InvXSectionSys,
                                                                          fitTsallisPi0PtMult,
                                                                          0, graphCombPi0InvXSectionSys->GetN());

          for (Int_t i = 0; i < 11; i++){
            if(directoryPi0[i]){
              cout << nameMeasGlobal[i].Data() << endl;
              graphPi0InvXSectionStat[i]      = ApplyXshiftIndividualSpectra (graphCombPi0InvXSectionTot,
                                                                              graphPi0InvXSectionStat[i],
                                                                              fitTsallisPi0PtMult,
                                                                              offSetPi0Shifting[i], nComBinsPi0Shifting[i]);
              graphPi0InvXSectionSys[i]       = ApplyXshiftIndividualSpectra (graphCombPi0InvXSectionTot,
                                                                              graphPi0InvXSectionSys[i],
                                                                              fitTsallisPi0PtMult,
                                                                              offSetPi0Shifting[i], nComBinsPi0Shifting[i]);
            }
          }

          TF1* fitTsallisPi0PtMultFromShift                 = FitObject("tmpt","TsallisMultWithPtPi0FromShift","Pi0");
          fitTsallisPi0PtMultFromShift->SetRange(0.4,10.);
          fitTsallisPi0PtMultFromShift->SetParameters(fitTsallisPi0PtMult->GetParameter(0),fitTsallisPi0PtMult->GetParameter(1), fitTsallisPi0PtMult->GetParameter(2));

          TF1* fitTsallisPi0PtMultFromShiftScaled = new TF1("TsallisMultWithPtPi0FromShiftScaled","(1/x)*TsallisMultWithPtPi0FromShift",0.4,10.);

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

          TGraphAsymmErrors* graphCombPi0InvXSectionTotUnShifted_clone = (TGraphAsymmErrors*) graphCombPi0InvXSectionTotUnShifted->Clone("graphCombPi0InvXSectionTotUnShifted_clone");

          Int_t numberPoints   = graphCombPi0InvXSectionTotUnShifted_clone->GetN();
          Double_t *xPoint     = graphCombPi0InvXSectionTotUnShifted_clone->GetX();
          Double_t* xvalueErrUp  = graphCombPi0InvXSectionTotUnShifted_clone->GetEXhigh();
          Double_t* xvalueErrLow = graphCombPi0InvXSectionTotUnShifted_clone->GetEXlow();
          Double_t *xPointShift= graphCombPi0InvXSectionTot->GetX();
          for (Int_t i=0; i<numberPoints; i++) {
            graphCombPi0InvXSectionTotUnShifted_clone->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
            graphCombPi0InvXSectionTotUnShifted_clone->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
          }
          DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionTotUnShifted_clone, markerStyleComb, markerSizeComb, colorComb , colorComb);
          graphCombPi0InvXSectionTotUnShifted_clone->Draw("p same");

          TLatex *labelRatioToFitBinShift   = new TLatex(0.685, 0.91, collisionSystem900GeV.Data());
          SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel,4);
          labelRatioToFitBinShift->SetTextFont(43);
          labelRatioToFitBinShift->Draw();
          TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.852, 0.86, "ALICE");
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
          histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,minPtPi0, maxPtPi0,1000,minXSectionPi0,maxXSectionPi0);
          SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);
          histo2DDummy2->DrawCopy();


          for (Int_t i = 0; i < 11; i++){
            if(directoryPi0[i]){
              DrawGammaSetMarkerTGraphAsym(graphPi0InvXSectionStat[i], markerStyleDet[i] ,markerSizeDet[i]/2, colorDet[i], colorDet[i]);
              graphPi0InvXSectionStat[i]->Draw("pEsame");
              DrawGammaSetMarkerTGraphAsym(graphPi0InvXSectionSys[i], markerStyleDet[i] ,markerSizeDet[i]/2, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
              graphPi0InvXSectionSys[i]->Draw("pEsame");
            }
          }

          DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatUnShifted, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
          graphCombPi0InvXSectionStatUnShifted->Draw("pEsame");
          DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
          graphCombPi0InvXSectionStat->Draw("pEsame");

          fitInvXSectionPi0->SetLineColor(kBlue+2);
          fitInvXSectionPi0->Draw("same");

          fitTsallisPi0PtMultFromShiftScaled->SetLineColor(kRed+2);
          fitTsallisPi0PtMultFromShiftScaled->Draw("same");

          canvasDummy2->Update();
          canvasDummy2->SaveAs(Form("%s/ComparisonShiftedPi0_900GeV.%s",outputDir.Data(),suffix.Data()));
          delete canvasDummy2;
      }

      TGraphAsymmErrors* graphCombPi0InvXSectionStat_WOXErr = (TGraphAsymmErrors*) graphCombPi0InvXSectionStat->Clone("graphCombPi0InvXSectionStatA_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphCombPi0InvXSectionStat_WOXErr);

      TGraphAsymmErrors* graphPi0InvXSectionStat_WOXErr[11];
      for (Int_t i = 0; i < 11; i++){
        if(directoryPi0[i]){
          graphPi0InvXSectionStat_WOXErr[i] = (TGraphAsymmErrors*) graphPi0InvXSectionStat[i]->Clone(Form("graphPi0InvXSectionStat_%i_WOXErr",i));
          ProduceGraphAsymmWithoutXErrors(graphPi0InvXSectionStat_WOXErr[i]);
        }
      }


      //**********************************************************************************************************************
      //************************************* Calculating bin shifted spectra & fitting **************************************
      //**********************************************************************************************************************

      // Cloning spectra
      TGraphAsymmErrors* graphCombEtaInvXSectionTotUnShifted    = (TGraphAsymmErrors*)graphCombEtaInvXSectionTot->Clone("Unshifted_Eta");
      TGraphAsymmErrors* graphCombEtaInvXSectionStatUnShifted   = (TGraphAsymmErrors*)graphCombEtaInvXSectionStat->Clone("UnshiftedStat_Eta");
      TGraphAsymmErrors* graphCombEtaInvXSectionSysUnShifted    = (TGraphAsymmErrors*)graphCombEtaInvXSectionSys->Clone("UnshiftedSys_Eta");

      TGraphAsymmErrors* graphEtaInvXSectionStatUnShifted[11];
      TGraphAsymmErrors* graphEtaInvXSectionSysUnShifted[11];
      for (Int_t i = 0; i < 11; i++){
        if(directoryEta[i]){
          graphEtaInvXSectionStatUnShifted[i] = (TGraphAsymmErrors*)graphEtaInvXSectionStat[i]->Clone(Form("UnshiftedStatEta%s",nameMeasGlobal[i].Data()));
          graphEtaInvXSectionSysUnShifted[i]  = (TGraphAsymmErrors*)graphEtaInvXSectionSys[i] ->Clone(Form("UnshiftedSysEta%s",nameMeasGlobal[i].Data()));
        }
      }

      Double_t paramGraphEta[3]           = {1.0e9, 7.9, 0.13};
      TF1* fitInvXSectionEta              = FitObject("l","fitInvXSectionEta","Eta",NULL,0.9,3.,paramGraphEta,"QNRMEX0+");
      fitInvXSectionEta->FixParameter(1,7.9);
      fitInvXSectionEta              = FitObject("l","fitInvXSectionEta","Eta",graphCombEtaInvXSectionTot,0.9,3.,paramGraphEta,"QNRMEX0+");

      TGraphAsymmErrors* graphCombEtaInvXSectionStat_WOXErr = (TGraphAsymmErrors*) graphCombEtaInvXSectionStat->Clone("graphCombEtaInvXSectionStatA_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphCombEtaInvXSectionStat_WOXErr);

      TGraphAsymmErrors* graphEtaInvXSectionStat_WOXErr[11];
      for (Int_t i = 0; i < 11; i++){
        if(directoryEta[i]){
          graphEtaInvXSectionStat_WOXErr[i] = (TGraphAsymmErrors*) graphEtaInvXSectionStat[i]->Clone(Form("graphEtaInvXSectionStat_%i_WOXErr",i));
          ProduceGraphAsymmWithoutXErrors(graphEtaInvXSectionStat_WOXErr[i]);
        }
      }

      // *************************************************************************************************************
      // Shift graphs in Y direction as well if desired
      // *************************************************************************************************************
      TGraphAsymmErrors* graphCombPi0InvXSectionTot_yShifted         = NULL;
      TGraphAsymmErrors* graphCombPi0InvXSectionStat_yShifted        = NULL;
      TGraphAsymmErrors* graphCombPi0InvXSectionSys_yShifted         = NULL;

      TGraphAsymmErrors* graphPi0InvXSectionStat_yShifted[11];
      TGraphAsymmErrors* graphPi0InvXSectionSys_yShifted[11];

      if(bWCorrection.Contains("Y") ){
          graphCombPi0InvXSectionTot_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotUnShifted->Clone("Pi0YShiftedCombTot");
          graphCombPi0InvXSectionTot_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvXSectionTot_yShifted, fitInvXSectionPi0);
          graphCombPi0InvXSectionStat_yShifted       = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatUnShifted->Clone("Pi0YShiftedCombStat");
          graphCombPi0InvXSectionStat_yShifted       =  ApplyYshiftIndividualSpectra( graphCombPi0InvXSectionStat_yShifted, fitInvXSectionPi0);
          graphCombPi0InvXSectionSys_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysUnShifted->Clone("Pi0YShiftedCombSys");
          graphCombPi0InvXSectionSys_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvXSectionSys_yShifted, fitInvXSectionPi0);

          for (Int_t i = 0; i < 11; i++){
            if(directoryPi0[i]){
              graphPi0InvXSectionStat_yShifted[i]       = (TGraphAsymmErrors*)graphPi0InvXSectionStatUnShifted[i]->Clone("Pi0YShiftedCombStat");
              graphPi0InvXSectionStat_yShifted[i]       =  ApplyYshiftIndividualSpectra( graphPi0InvXSectionStat_yShifted[i], fitInvXSectionPi0);
              graphPi0InvXSectionSys_yShifted[i]        = (TGraphAsymmErrors*)graphPi0InvXSectionSysUnShifted[i]->Clone("Pi0YShiftedCombSys");
              graphPi0InvXSectionSys_yShifted[i]        =  ApplyYshiftIndividualSpectra( graphPi0InvXSectionSys_yShifted[i], fitInvXSectionPi0);
            }
          }
      }

      // *************************************************************************************************************
      // Shift spectra in Y  direction as well if desired
      // *************************************************************************************************************

      TGraphAsymmErrors* graphCombEtaInvXSectionTot_yShifted         = NULL;
      TGraphAsymmErrors* graphCombEtaInvXSectionStat_yShifted        = NULL;
      TGraphAsymmErrors* graphCombEtaInvXSectionSys_yShifted         = NULL;

      TGraphAsymmErrors* graphEtaInvXSectionStat_yShifted[11];
      TGraphAsymmErrors* graphEtaInvXSectionSys_yShifted[11];

      if(bWCorrection.Contains("Y") ){
          graphCombEtaInvXSectionTot_yShifted        = (TGraphAsymmErrors*)graphCombEtaInvXSectionTotUnShifted->Clone("EtaYShiftedCombTot");
          graphCombEtaInvXSectionTot_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaInvXSectionTot_yShifted, fitInvXSectionEta);
          graphCombEtaInvXSectionStat_yShifted       = (TGraphAsymmErrors*)graphCombEtaInvXSectionStatUnShifted->Clone("EtaYShiftedCombStat");
          graphCombEtaInvXSectionStat_yShifted       =  ApplyYshiftIndividualSpectra( graphCombEtaInvXSectionStat_yShifted, fitInvXSectionEta);
          graphCombEtaInvXSectionSys_yShifted        = (TGraphAsymmErrors*)graphCombEtaInvXSectionSysUnShifted->Clone("EtaYShiftedCombSys");
          graphCombEtaInvXSectionSys_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaInvXSectionSys_yShifted, fitInvXSectionEta);

          for (Int_t i = 0; i < 11; i++){
            if(directoryEta[i]){
              graphEtaInvXSectionStat_yShifted[i]       = (TGraphAsymmErrors*)graphEtaInvXSectionStatUnShifted[i]->Clone("EtaYShiftedCombStat");
              graphEtaInvXSectionStat_yShifted[i]       =  ApplyYshiftIndividualSpectra( graphEtaInvXSectionStat_yShifted[i], fitInvXSectionEta);
              graphEtaInvXSectionSys_yShifted[i]        = (TGraphAsymmErrors*)graphEtaInvXSectionSysUnShifted[i]->Clone("EtaYShiftedCombSys");
              graphEtaInvXSectionSys_yShifted[i]        =  ApplyYshiftIndividualSpectra( graphEtaInvXSectionSys_yShifted[i], fitInvXSectionEta);
            }
          }
      }

      // *************************************************************************************************************
      // redo fitting after binshifts
      // *************************************************************************************************************
      // Tsallis function
      fitInvXSectionPi0           = FitObject("l","fitInvXSectionPi0900GeV","Pi0",graphCombPi0InvXSectionTot,0.3,25.,paramGraph,"QNRMEX0+");
      cout << WriteParameterToFile(fitInvXSectionPi0)<< endl;

      //Two component model from Bylinkin
      Double_t paramTCMPi0New[5]  = { 5E10,0.2,
                                      1E9,0.8,3.8};
      TF1* fitTCMInvXSectionPi0        = FitObject("tcm","fitTCMInvXSectionPi0900GeV","Pi0",graphCombPi0InvXSectionTot,0.4,10. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
      fitTCMInvXSectionPi0        = FitObject("tcm","fitTCMInvXSectionPi0900GeV","Pi0",graphCombPi0InvXSectionTot,0.4,10. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
      cout << WriteParameterToFile(fitTCMInvXSectionPi0)<< endl;

      Double_t paramPi0Power[3] = {1E11,0.5,6.5};
      TF1* fitPowInvXSectionPi0   = FitObject("powPure","fitPowInvXSectionPi0900GeV","Pi0",graphCombPi0InvXSectionTot,3.5,10. ,paramPi0Power,"QNRMEX0+","", kFALSE);
      cout << WriteParameterToFile(fitPowInvXSectionPi0)<< endl;

      TF1* fitPowInvXSectionPi0Stat   = FitObject("powPure","fitPowInvXSectionPi0900GeVStat","Pi0",graphCombPi0InvXSectionStat,3.5,10. ,paramPi0Power,"QNRMEX0+","", kFALSE);
      cout << WriteParameterToFile(fitPowInvXSectionPi0Stat)<< endl;

      Double_t paramPi0HageDorn[5] = {1E11,0.3,-0.1,0.5,5.95};
      TF1* fitOHagInvYieldPi0Tot   = FitObject("oHag","fitOHagInvYieldPi0900GeV","Pi0",graphCombPi0InvXSectionTot,0.4,10. ,paramPi0HageDorn,"QNRMEX0+","", kFALSE);
      cout << WriteParameterToFile(fitOHagInvYieldPi0Tot)<< endl;

      fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      fLog << "Pi0 - Tsallis" << endl;
      fLog << WriteParameterToFile(fitInvXSectionPi0)<< endl;
      fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      fLog << "Pi0 - TCM" << endl;
      fLog << WriteParameterToFile(fitTCMInvXSectionPi0) << endl;
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
      cout << WriteParameterToFile(fitInvXSectionEta)<< endl;

      fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      fLog << "Eta - Tsallis" << endl;
      fLog << WriteParameterToFile(fitInvXSectionEta)<< endl;


      //********************************************************************************************************
      // Plotting simple comparison of data vs fit to eta meson spec
      //********************************************************************************************************
      TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
      DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
      canvasDummy2->SetLogy();
      canvasDummy2->SetLogx();
      TH2F* histo2DDummy3;
      histo2DDummy3               = new TH2F("histo2DDummy3","histo2DDummy3",1000,minPtEta,maxPtEta,1000,minXSectionEta,maxXSectionEta);
      SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);
      histo2DDummy3->DrawCopy();

      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatUnShifted, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);

      TLatex *labelRelSysErrEnergyC    = new TLatex(0.18,0.94,collisionSystem900GeV.Data());
      SetStyleTLatex( labelRelSysErrEnergyC, 0.85*textSizeLabelsPixel,4);
      labelRelSysErrEnergyC->SetTextFont(43);
      labelRelSysErrEnergyC->Draw();
      TLatex *labelRelSysErrEtaC       = new TLatex(0.18,0.9,"#eta #rightarrow #gamma#gamma");
      SetStyleTLatex( labelRelSysErrEtaC, 0.85*textSizeLabelsPixel,4);
      labelRelSysErrEtaC->SetTextFont(43);
      labelRelSysErrEtaC->Draw();

  //********************************************************************************************************
      canvasDummy2->Clear();
      histo2DDummy3->DrawCopy();

      graphCombEtaInvXSectionStatUnShifted->Draw("pEsame");
      graphCombEtaInvXSectionStat->Draw("pEsame");

      fitInvXSectionEta->SetLineColor(kRed+2);
      fitInvXSectionEta->Draw("same");

      TLatex *labelTCMEta10 = new TLatex(0.35, 0.90, Form("dN/dy: (%.1e #pm %.1e) - n: (%.3f #pm %.3f) - T_{Levy} (GeV/c): (%.3f #pm %.3f)",fitInvXSectionEta->GetParameter(0),fitInvXSectionEta->GetParError(0),fitInvXSectionEta->GetParameter(1),fitInvXSectionEta->GetParError(1),fitInvXSectionEta->GetParameter(2),fitInvXSectionEta->GetParError(2)));
      SetStyleTLatex( labelTCMEta10, 0.02,4);
      labelTCMEta10->Draw();

      labelRelSysErrEnergyC->Draw();
      labelRelSysErrEtaC->Draw();

      TLegend* legendWithFitEta2   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
      legendWithFitEta2->AddEntry(fitInvXSectionEta,"Tsallis","l");
      legendWithFitEta2->Draw();

      canvasDummy2->Update();
      canvasDummy2->Print(Form("%s/ComparisonWithFit_Tsallis_Eta_900GeV.%s",outputDir.Data(),suffix.Data()));

      delete histo2DDummy3;
      canvasDummy2->Clear();

      //********************************************************************************************************
      // Plotting simple comparison of data vs fit to pi0 meson spec
      //********************************************************************************************************

      TF1* fitTCMDecomposedPi0L                 = new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
      fitTCMDecomposedPi0L->SetParameters(fitTCMInvXSectionPi0->GetParameter(0),fitTCMInvXSectionPi0->GetParameter(1));
      fitTCMDecomposedPi0L->SetRange(minPtPi0,maxPtPi0);
      TF1 *fitTCMDecomposedPi0H                 = new TF1("twoCompModel_DecH","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]))");
     //      graphCombEtaInvXSectionTotA->Fit(fitTCMDecomposedH,"QNRMEX0+","",5,20);
      fitTCMDecomposedPi0H->SetParameters(fitTCMInvXSectionPi0->GetParameter(2),fitTCMInvXSectionPi0->GetParameter(3), fitTCMInvXSectionPi0->GetParameter(4));
      fitTCMDecomposedPi0H->SetRange(minPtPi0,maxPtPi0);

      histo2DDummy3               = new TH2F("histo2DDummy2","histo2DDummy2",1000,minPtPi0,maxPtPi0,1000,minXSectionPi0,maxXSectionPi0);
      SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);
      histo2DDummy3->DrawCopy();

      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatUnShifted, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
      graphCombPi0InvXSectionStatUnShifted->Draw("pEsame");
      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
      graphCombPi0InvXSectionStat->Draw("pEsame");

      fitTCMInvXSectionPi0->SetLineColor(kRed+2);
      fitTCMInvXSectionPi0->SetRange(minPtPi0,maxPtPi0);
      fitTCMInvXSectionPi0->Draw("same");

      fitTCMDecomposedPi0L->SetLineColor(kAzure);
      fitTCMDecomposedPi0L->SetLineStyle(2);
      fitTCMDecomposedPi0L->Draw("same");
      fitTCMDecomposedPi0H->SetLineColor(kGreen+2);
      fitTCMDecomposedPi0H->SetLineStyle(8);
      fitTCMDecomposedPi0H->Draw("same");

      TLatex *labelTCMPi01= new TLatex(0.48, 0.94, Form("TCM low:"));
      TLatex *labelTCMPi02= new TLatex(0.48, 0.90, Form("A_{1}: (%.1e #pm %.1e) - T_{e}: (%.3f #pm %.3f)",fitTCMInvXSectionPi0->GetParameter(0),fitTCMInvXSectionPi0->GetParError(0),fitTCMInvXSectionPi0->GetParameter(1),fitTCMInvXSectionPi0->GetParError(1)));
      TLatex *labelTCMPi03= new TLatex(0.48, 0.86, Form("TCM high:"));
      TLatex *labelTCMPi04= new TLatex(0.48, 0.82, Form("A_{2}: (%.1e #pm %.1e) - T: (%.3f #pm %.3f) - n: (%.3f #pm %.3f)",fitTCMInvXSectionPi0->GetParameter(2),fitTCMInvXSectionPi0->GetParError(2),abs(fitTCMInvXSectionPi0->GetParameter(3)),fitTCMInvXSectionPi0->GetParError(3),fitTCMInvXSectionPi0->GetParameter(4),fitTCMInvXSectionPi0->GetParError(4)));

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
      legendWithFit->AddEntry(fitTCMInvXSectionPi0,"Bylinkin-Rostovtsev (TCM)","l");
      legendWithFit->Draw();

      canvasDummy2->Update();
      canvasDummy2->Print(Form("%s/ComparisonWithFitPi0_900GeV.%s",outputDir.Data(),suffix.Data()));
      //********************************************************************************************************
      canvasDummy2->Clear();
      histo2DDummy3->DrawCopy();

      graphCombPi0InvXSectionStatUnShifted->Draw("pEsame");
      graphCombPi0InvXSectionStat->Draw("pEsame");

      fitInvXSectionPi0->SetLineColor(kRed+2);
      fitInvXSectionPi0->Draw("same");

      TLatex *labelTCMEta20 = new TLatex(0.35, 0.90, Form("dN/dy: (%.1e #pm %.1e) - n: (%.3f #pm %.3f) - T_{Levy} (GeV/c): (%.3f #pm %.3f)",fitInvXSectionPi0->GetParameter(0),fitInvXSectionPi0->GetParError(0),fitInvXSectionPi0->GetParameter(1),fitInvXSectionPi0->GetParError(1),fitInvXSectionPi0->GetParameter(2),fitInvXSectionPi0->GetParError(2)));
      SetStyleTLatex( labelTCMEta20, 0.02,4);
      labelTCMEta20->Draw();

      labelRelSysErrEnergyC->Draw();
      labelRelSysErrPi0C->Draw();

      TLegend* legendWithFitPi02   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
      legendWithFitPi02->AddEntry(fitInvXSectionPi0,"Tsallis","l");
      legendWithFitPi02->Draw();

      canvasDummy2->Update();
      canvasDummy2->Print(Form("%s/ComparisonWithFit_Tsallis_Pi0_900GeV.%s",outputDir.Data(),suffix.Data()));

      delete canvasDummy2;
      delete histo2DDummy3;


      // *************************************************************************************************************
      // Calculate ratios to combined fit
      // *************************************************************************************************************
      TH1D* histoRatioPythia8ToFit                     = (TH1D*) histoPythia8InvXSection->Clone();
      histoRatioPythia8ToFit                           = CalculateHistoRatioToFit (histoRatioPythia8ToFit, fitTCMInvXSectionPi0);
      TGraphErrors* graphRatioPythia8ToFit             = (TGraphErrors*) graphPythia8InvXSection->Clone();
      graphRatioPythia8ToFit                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFit, fitTCMInvXSectionPi0);

      TGraph* graphRatioPi0CombNLOMuHalf               = (TGraph*)graphNLOCalcPi0MuHalf->Clone();cout << __LINE__ << endl;
      TGraph* graphRatioPi0CombNLOMuOne                = (TGraph*)graphNLOCalcPi0MuOne->Clone();cout << __LINE__ << endl;
      TGraph* graphRatioPi0CombNLOMuTwo                = (TGraph*)graphNLOCalcPi0MuTwo->Clone();cout << __LINE__ << endl;
      TGraphAsymmErrors* graphRatioPi0DSS07            = (TGraphAsymmErrors*)graphPi0DSS07->Clone();cout << __LINE__ << endl;

      graphRatioPi0CombNLOMuHalf                       = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuHalf, fitTCMInvXSectionPi0); cout << __LINE__ << endl;
      graphRatioPi0CombNLOMuOne                        = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuOne, fitTCMInvXSectionPi0); cout << __LINE__ << endl;
      graphRatioPi0CombNLOMuTwo                        = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuTwo, fitTCMInvXSectionPi0); cout << __LINE__ << endl;
      graphRatioPi0DSS07                               = CalculateGraphErrRatioToFit (graphRatioPi0DSS07, fitTCMInvXSectionPi0); cout << __LINE__ << endl;

      TH1D* histoRatioPythia8ToFitEta                  = (TH1D*) histoPythia8InvXSectionEta->Clone();
      histoRatioPythia8ToFitEta                        = CalculateHistoRatioToFit (histoRatioPythia8ToFitEta, fitInvXSectionEta);
      histoRatioPythia8ToFitEta->GetXaxis()->SetRangeUser(minPtEta,maxPtEta);

      TGraphErrors* graphRatioPythia8ToFitEta             = (TGraphErrors*) graphPythia8InvXSectionEta->Clone();
      graphRatioPythia8ToFitEta                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFitEta, fitInvXSectionEta);
      while(graphRatioPythia8ToFitEta->GetX()[0] < minPtEta) graphRatioPythia8ToFitEta->RemovePoint(0);

      TGraph* graphRatioEtaCombNLOMuHalf                  = (TGraph*)graphNLOCalcEtaMuHalf->Clone();
      TGraph* graphRatioEtaCombNLOMuOne                   = (TGraph*)graphNLOCalcEtaMuOne->Clone();
      TGraph* graphRatioEtaCombNLOMuTwo                   = (TGraph*)graphNLOCalcEtaMuTwo->Clone();
      graphRatioEtaCombNLOMuHalf                          = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuHalf, fitInvXSectionEta);
      graphRatioEtaCombNLOMuOne                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuOne, fitInvXSectionEta);
      graphRatioEtaCombNLOMuTwo                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuTwo, fitInvXSectionEta);

      TGraphAsymmErrors* graphRatioEtaAESSS               = (TGraphAsymmErrors*) graphEtaAESSS->Clone();
      graphRatioEtaAESSS                                  = CalculateGraphErrRatioToFit (graphRatioEtaAESSS, fitInvXSectionEta);

      // *************************************************************************************************************

      TGraphAsymmErrors* graphRatioCombCombFitTot     = (TGraphAsymmErrors*)graphCombPi0InvXSectionTot->Clone();
      graphRatioCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitTot, fitTCMInvXSectionPi0);
      TGraphAsymmErrors* graphRatioCombCombFitStat    = (TGraphAsymmErrors*)graphCombPi0InvXSectionStat->Clone();
      graphRatioCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioCombCombFitStat, fitTCMInvXSectionPi0);
      TGraphAsymmErrors* graphRatioCombCombFitSys     = (TGraphAsymmErrors*)graphCombPi0InvXSectionSys->Clone();
      graphRatioCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitSys, fitTCMInvXSectionPi0);

      TGraphAsymmErrors* graphRatioCombCombFitTotEta     = (TGraphAsymmErrors*)graphCombEtaInvXSectionTot->Clone();
      graphRatioCombCombFitTotEta                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitTotEta, fitInvXSectionEta);
      TGraphAsymmErrors* graphRatioCombCombFitStatEta    = (TGraphAsymmErrors*)graphCombEtaInvXSectionStat->Clone();
      graphRatioCombCombFitStatEta                       = CalculateGraphErrRatioToFit(graphRatioCombCombFitStatEta, fitInvXSectionEta);
      TGraphAsymmErrors* graphRatioCombCombFitSysEta     = (TGraphAsymmErrors*)graphCombEtaInvXSectionSys->Clone();
      graphRatioCombCombFitSysEta                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitSysEta, fitInvXSectionEta);

      TGraphAsymmErrors* graphRatioCombFitStat[11];
      TGraphAsymmErrors* graphRatioCombFitSys[11];
      for (Int_t i = 0; i < 11; i++){
        if(directoryPi0[i]){
          graphRatioCombFitStat[i]                = (TGraphAsymmErrors*)graphPi0InvXSectionStat[i]->Clone();
          graphRatioCombFitStat[i]                = CalculateGraphErrRatioToFit(graphRatioCombFitStat[i], fitTCMInvXSectionPi0);
          graphRatioCombFitSys[i]                 = (TGraphAsymmErrors*)graphPi0InvXSectionSys[i]->Clone();
          graphRatioCombFitSys[i]                 = CalculateGraphErrRatioToFit(graphRatioCombFitSys[i], fitTCMInvXSectionPi0);
        }
      }
      TGraphAsymmErrors* graphRatioCombFitStatEta[11];
      TGraphAsymmErrors* graphRatioCombFitSysEta[11];
      for (Int_t i = 0; i < 11; i++){
         if(directoryEta[i]){
              graphRatioCombFitStatEta[i]                = (TGraphAsymmErrors*)graphEtaInvXSectionStat[i]->Clone();
              graphRatioCombFitStatEta[i]                = CalculateGraphErrRatioToFit(graphRatioCombFitStatEta[i], fitInvXSectionEta);
              graphRatioCombFitSysEta[i]                 = (TGraphAsymmErrors*)graphEtaInvXSectionSys[i]->Clone();
              graphRatioCombFitSysEta[i]                 = CalculateGraphErrRatioToFit(graphRatioCombFitSysEta[i], fitInvXSectionEta);
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

      TGraphAsymmErrors* graphRatioCombFitStat_WOXErr[11];
      TGraphAsymmErrors* graphRatioCombFitSys_WOXErr[11];
      for (Int_t i = 0; i < 11; i++){
        if(directoryPi0[i]){
          graphRatioCombFitStat_WOXErr[i] = (TGraphAsymmErrors*) graphRatioCombFitStat[i]->Clone(Form("graphRatioCombFitStat_WOXErr_%i",i));
          ProduceGraphAsymmWithoutXErrors(graphRatioCombFitStat_WOXErr[i]);
          graphRatioCombFitSys_WOXErr[i] = (TGraphAsymmErrors*) graphRatioCombFitSys[i]->Clone(Form("graphRatioCombCombFitStat_WOXErr_%i",i));
          ProduceGraphAsymmWithoutXErrors(graphRatioCombFitSys_WOXErr[i]);
        }
      }
      TGraphAsymmErrors* graphRatioCombFitStatEta_WOXErr[11];
      TGraphAsymmErrors* graphRatioCombFitSysEta_WOXErr[11];
      for (Int_t i = 0; i < 11; i++){
         if(directoryEta[i]){
           graphRatioCombFitStatEta_WOXErr[i] = (TGraphAsymmErrors*) graphRatioCombFitStatEta[i]->Clone(Form("graphRatioCombFitStatEta_WOXErr_%i",i));
           ProduceGraphAsymmWithoutXErrors(graphRatioCombFitStatEta_WOXErr[i]);
           graphRatioCombFitSysEta_WOXErr[i] = (TGraphAsymmErrors*) graphRatioCombFitSysEta[i]->Clone(Form("graphRatioCombCombFitStatEta_WOXErr_%i",i));
           ProduceGraphAsymmWithoutXErrors(graphRatioCombFitSysEta_WOXErr[i]);
         }
      }

      TGraphAsymmErrors* graphRatioPi0CombTsallisFitStat  = (TGraphAsymmErrors*)graphCombPi0InvXSectionStat->Clone();
      graphRatioPi0CombTsallisFitStat                     = CalculateGraphErrRatioToFit(graphRatioPi0CombTsallisFitStat, fitInvXSectionPi0);
      TGraphAsymmErrors* graphRatioPi0CombHagedornFitStat = (TGraphAsymmErrors*)graphCombPi0InvXSectionStat->Clone();
      graphRatioPi0CombHagedornFitStat                    = CalculateGraphErrRatioToFit(graphRatioPi0CombHagedornFitStat, fitOHagInvYieldPi0Tot);
      TGraphAsymmErrors* graphRatioPi0CombPowerFitStat    = (TGraphAsymmErrors*)graphCombPi0InvXSectionStat->Clone();
      graphRatioPi0CombPowerFitStat                       = CalculateGraphErrRatioToFit(graphRatioPi0CombPowerFitStat, fitPowInvXSectionPi0);
      TGraphAsymmErrors* graphRatioPi0CombTsallisFitSys   = (TGraphAsymmErrors*)graphCombPi0InvXSectionSys->Clone();
      graphRatioPi0CombTsallisFitSys                      = CalculateGraphErrRatioToFit(graphRatioPi0CombTsallisFitSys, fitInvXSectionPi0);
      TGraphAsymmErrors* graphRatioPi0CombHagedornFitSys  = (TGraphAsymmErrors*)graphCombPi0InvXSectionSys->Clone();
      graphRatioPi0CombHagedornFitSys                     = CalculateGraphErrRatioToFit(graphRatioPi0CombHagedornFitSys, fitOHagInvYieldPi0Tot);
      TGraphAsymmErrors* graphRatioPi0CombPowerFitSys     = (TGraphAsymmErrors*)graphCombPi0InvXSectionSys->Clone();
      graphRatioPi0CombPowerFitSys                        = CalculateGraphErrRatioToFit(graphRatioPi0CombPowerFitSys, fitPowInvXSectionPi0);

      TGraphAsymmErrors* graphRatioPi0CombTsallisFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombTsallisFitStat->Clone("graphRatioPi0CombTsallisFitStat_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombTsallisFitStat_WOXErr);
      TGraphAsymmErrors* graphRatioPi0CombHagedornFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombHagedornFitStat->Clone("graphRatioPi0CombHagedornFitStat_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombHagedornFitStat_WOXErr);
      TGraphAsymmErrors* graphRatioPi0CombPowerFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombPowerFitStat->Clone("graphRatioPi0CombPowerFitStat_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombPowerFitStat_WOXErr);

      TGraphAsymmErrors* graphRatioEtaCombTsallisFitStat  = (TGraphAsymmErrors*)graphCombEtaInvXSectionStat->Clone();
      graphRatioEtaCombTsallisFitStat                     = CalculateGraphErrRatioToFit(graphRatioEtaCombTsallisFitStat, fitInvXSectionEta);
      TGraphAsymmErrors* graphRatioEtaCombTsallisFitSys   = (TGraphAsymmErrors*)graphCombEtaInvXSectionSys->Clone();
      graphRatioEtaCombTsallisFitSys                      = CalculateGraphErrRatioToFit(graphRatioEtaCombTsallisFitSys, fitInvXSectionEta);

      TGraphAsymmErrors* graphRatioEtaCombTsallisFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaCombTsallisFitStat->Clone("graphRatioEtaCombTsallisFitStat_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombTsallisFitStat_WOXErr);

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
      histo2DPi0RatioToCombFit               = new TH2F("histo2DPi0RatioToCombFit","histo2DPi0RatioToCombFit",1000,minPtPi0, maxPtPi0,1000,0.1,3.5);
      SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/TCM fit", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                                0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
      histo2DPi0RatioToCombFit->GetXaxis()->SetMoreLogLabels();
      histo2DPi0RatioToCombFit->GetXaxis()->SetNoExponent(kTRUE);
  //  histo2DPi0RatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
  //  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
      histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.2,2.3);
      histo2DPi0RatioToCombFit->Draw("copy");

          DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
          graphRatioCombCombFitSys->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
          graphRatioCombCombFitStat_WOXErr->Draw("p,same,z");

          DrawGammaLines(minPtPi0,maxPtPi0 , 1., 1.,0.1, kGray+2);
          DrawGammaLines(minPtPi0,maxPtPi0 , 1.1, 1.1,0.1, kGray, 7);
          DrawGammaLines(minPtPi0,maxPtPi0 , 0.9, 0.9,0.1, kGray, 7);

          TLatex *labelRatioToFitEnergy   = new TLatex(0.67, 0.91, collisionSystem900GeV.Data());
          SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel,4);
          labelRatioToFitEnergy->SetTextFont(43);
          labelRatioToFitEnergy->Draw();
          TLatex *labelRatioToFitALICE    = new TLatex(0.852, 0.86, "ALICE");
          SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel,4);
          labelRatioToFitALICE->SetTextFont(43);
          labelRatioToFitALICE->Draw();
          TLatex *labelRatioToFitPi0      = new TLatex(0.826, 0.807, "#pi^{0} #rightarrow #gamma#gamma");
          SetStyleTLatex( labelRatioToFitPi0, textSizeLabelsPixel,4);
          labelRatioToFitPi0->SetTextFont(43);
          labelRatioToFitPi0->Draw();

      canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_PP8TeV.%s",outputDir.Data(),suffix.Data()));

      // **********************************************************************************************************************
      // *******************************************Plot different ratios to fits *********************************************
      // **********************************************************************************************************************

      histo2DPi0RatioToCombFit->SetYTitle("Data/fit");
      histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,3.4);
      histo2DPi0RatioToCombFit->Draw("copy");


      while (graphRatioPi0CombPowerFitSys->GetX()[0] < 1.6) graphRatioPi0CombPowerFitSys->RemovePoint(0);
      while (graphRatioPi0CombPowerFitStat_WOXErr->GetX()[0] < 1.6) graphRatioPi0CombPowerFitStat_WOXErr->RemovePoint(0);
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
          legendRatioPi0Fits->AddEntry(graphRatioPi0CombTsallisFitSys,"Tsallis","p");
          legendRatioPi0Fits->AddEntry(graphRatioPi0CombHagedornFitSys,"mod. Hagedorn","p");
          legendRatioPi0Fits->AddEntry(graphRatioPi0CombPowerFitSys,"pure powerlaw, 3.5-25 GeV/#it{c}","p");
          legendRatioPi0Fits->Draw();

      canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToDifferentFits_PP900GeV.%s",outputDir.Data(),suffix.Data()));
      histo2DPi0RatioToCombFit->SetYTitle("Data/TCM fit");

      // **********************************************************************************************************************
      // *******************************************Plot Ratio of Individual meas to Fit ******************************************
      // **********************************************************************************************************************

      canvasRatioToCombFit->cd();
      histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,3.05);
      histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(minPtPi0,maxPtPi0);
      histo2DPi0RatioToCombFit->Draw("copy");

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
      labelRatioToFitALICE->Draw();
      labelRatioToFitPi0->Draw();

      //****************************** Definition of the Legend ******************************************
      //**************** Row def ************************
      Double_t rowsLegendOnlyPi0Ratio[7]      = {0.92,0.88,0.84,0.80,0.76,0.76,0.68};
      Double_t rowsLegendOnlyPi0RatioAbs[7]   = {0.91,2.75,2.61,2.47,2.33,2.33,2.05};
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
      TLatex *textStatOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[1]+0.02,rowsLegendOnlyPi0Ratio[0]+0.01 ,"stat");
      SetStyleTLatex( textStatOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
      textStatOnlyRatioPi0->SetTextFont(43);
      textStatOnlyRatioPi0->Draw();
      TLatex *textSysOnlyRatioPi0             = new TLatex(columnsLegendOnlyPi0Ratio[2]+0.06 ,rowsLegendOnlyPi0Ratio[0]+0.01,"syst");
      SetStyleTLatex( textSysOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
      textSysOnlyRatioPi0->SetTextFont(43);
      textSysOnlyRatioPi0->Draw();

      TMarker* markerPi0OnlyRatio[11];
      TBox* boxPi0OnlyRatio[11];

      for (Int_t i = 0; i < 11; i++){
        if(directoryPi0[i]){
          markerPi0OnlyRatio[i]               = CreateMarkerFromGraph(graphRatioCombFitSys[i],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[i+1],1);
          markerPi0OnlyRatio[i]->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[i+1]);
          boxPi0OnlyRatio[i]                  = CreateBoxFromGraph(graphRatioCombFitSys[i], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[i+1]- heightBox,
          columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox+0.1, rowsLegendOnlyPi0RatioAbs[i+1]+ heightBox);
          boxPi0OnlyRatio[i]->Draw("l");
        }
      }

      canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_PP900GeV.%s",outputDir.Data(),suffix.Data()));

      // **********************************************************************************************************************
      // ******************************************* Ratio of Comb to Fit ****************************************
      // **********************************************************************************************************************
      textSizeLabelsPixel                 = 48;
      canvasRatioToCombFit->cd();
      TH2F * histo2DEtaRatioToCombFit;
      histo2DEtaRatioToCombFit               = new TH2F("histo2DEtaRatioToCombFit","histo2DEtaRatioToCombFit",1000,minPtEta,maxPtEta,1000,0.2,7.    );
      SetStyleHistoTH2ForGraphs(histo2DEtaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Tsallis fit", 0.85*textsizeLabelsPP, textsizeLabelsPP,
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

          DrawGammaLines(0.33, 50. , 1., 1.,0.5, kGray+2);
          DrawGammaLines(0.33, 50. , 1.1, 1.1,0.5, kGray, 7);
          DrawGammaLines(0.33, 50. , 0.9, 0.9,0.5, kGray, 7);

          TLatex *labelRatioToFitEnergy2   = new TLatex(0.695, 0.91, collisionSystem900GeV.Data());
          SetStyleTLatex( labelRatioToFitEnergy2, textSizeLabelsPixel,4);
          labelRatioToFitEnergy2->SetTextFont(43);
          labelRatioToFitEnergy2->Draw();
          TLatex *labelRatioToFitALICE2    = new TLatex(0.85, 0.86, "ALICE");
          SetStyleTLatex( labelRatioToFitALICE2, textSizeLabelsPixel,4);
          labelRatioToFitALICE2->SetTextFont(43);
          labelRatioToFitALICE2->Draw();

          TLatex *labelRatioToFitEta      = new TLatex(0.84,0.82,"#eta #rightarrow #gamma#gamma");
          SetStyleTLatex( labelRatioToFitEta, textSizeLabelsPixel,4);
          labelRatioToFitEta->SetTextFont(43);
          labelRatioToFitEta->Draw();

      canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_PP900GeV.%s",outputDir.Data(),suffix.Data()));

      // **********************************************************************************************************************
      // ******************************************* Ratio of Individual meas to Fit ******************************************
      // **********************************************************************************************************************

      canvasRatioToCombFit->cd();
      histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.1,2.45);
      histo2DEtaRatioToCombFit->GetXaxis()->SetRangeUser(minPtEta,maxPtEta);
      histo2DEtaRatioToCombFit->Draw("copy");

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
      labelRatioToFitALICE->Draw();
      labelRatioToFitEta->Draw();

      //****************************** Definition of the Legend ******************************************
      //**************** Row def ************************
      Double_t rowsLegendOnlyEtaRatio[6]      = {0.92,0.88,0.84,0.80,0.79,0.76};
      Double_t rowsLegendOnlyEtaRatioAbs[6]   = {0.91,2.2,2.1,2.0,1.95,1.9};
      Double_t columnsLegendOnlyEtaRatio[3]   = {0.15,0.32, 0.38};
      Double_t columnsLegendOnlyEtaRatioAbs[3]= {0.15,1.04, 1.37};
      //****************** first Column **************************************************
      TLatex *textSingleMeasRatioEta[11];
      for (Int_t i = 0; i < 11; i++){
        if(directoryEta[i]){
          textSingleMeasRatioEta[i]           = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[i+1],nameMeasGlobal[i].Data());
          SetStyleTLatex( textSingleMeasRatioEta[i], 0.85*textSizeLabelsPixel,4);
          textSingleMeasRatioEta[i]->SetTextFont(43);
          textSingleMeasRatioEta[i]->Draw();
        }
      }

      //****************** second Column *************************************************
      TLatex *textStatOnlyRatioEta            = new TLatex(columnsLegendOnlyEtaRatio[1]-0.08,rowsLegendOnlyEtaRatio[0] ,"stat");
      SetStyleTLatex( textStatOnlyRatioEta, 0.85*textSizeLabelsPixel,4);
      textStatOnlyRatioEta->SetTextFont(43);
      textStatOnlyRatioEta->Draw();
      TLatex *textSysOnlyRatioEta             = new TLatex(columnsLegendOnlyEtaRatio[2]+0.03 ,rowsLegendOnlyEtaRatio[0],"syst");
      SetStyleTLatex( textSysOnlyRatioEta, 0.85*textSizeLabelsPixel,4);
      textSysOnlyRatioEta->SetTextFont(43);
      textSysOnlyRatioEta->Draw();

      TMarker* markerEtaOnlyRatio[11];
      TBox* boxEtaOnlyRatio[11];
      for (Int_t i = 0; i < 11; i++){
        if(directoryEta[i]){
          markerEtaOnlyRatio[i]               = CreateMarkerFromGraph(graphRatioCombFitSysEta[i],columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[i+1],1);
          markerEtaOnlyRatio[i]->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[i+1]);
          boxEtaOnlyRatio[i]                  = CreateBoxFromGraph(graphRatioCombFitSysEta[i], columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[i+1]- heightBox,
          columnsLegendOnlyEtaRatioAbs[2]+ 3*lengthBox+0.1, rowsLegendOnlyEtaRatioAbs[i+1]+ heightBox);
          boxEtaOnlyRatio[i]->Draw("l");
        }
      }


      canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_PP900GeV.%s",outputDir.Data(),suffix.Data()));

      // *************************************************************************************************************
      // Calculate ratios to combined TSALLIS fit
      // *************************************************************************************************************

      TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitTot     = (TGraphAsymmErrors*)graphCombPi0InvXSectionTot->Clone();
      graphRatioPi0TsallisCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioPi0TsallisCombCombFitTot, fitInvXSectionPi0);
      TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitStat    = (TGraphAsymmErrors*)graphCombPi0InvXSectionStat->Clone();
      graphRatioPi0TsallisCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioPi0TsallisCombCombFitStat, fitInvXSectionPi0);
      TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitSys     = (TGraphAsymmErrors*)graphCombPi0InvXSectionSys->Clone();
      graphRatioPi0TsallisCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioPi0TsallisCombCombFitSys, fitInvXSectionPi0);

      TGraphAsymmErrors* graphRatioPi0TsallisCombFitSys[11];
      TGraphAsymmErrors* graphRatioPi0TsallisCombFitStat[11];
      for (Int_t i = 0; i < 11; i++){
         if(directoryPi0[i]){
           graphRatioPi0TsallisCombFitSys[i]                    = (TGraphAsymmErrors*)graphPi0InvXSectionSys[i]->Clone();
           graphRatioPi0TsallisCombFitSys[i]                    = CalculateGraphErrRatioToFit(graphRatioPi0TsallisCombFitSys[i], fitInvXSectionPi0);
           graphRatioPi0TsallisCombFitStat[i]                   = (TGraphAsymmErrors*)graphPi0InvXSectionStat[i]->Clone();
           graphRatioPi0TsallisCombFitStat[i]                   = CalculateGraphErrRatioToFit(graphRatioPi0TsallisCombFitStat[i], fitInvXSectionPi0);
         }
      }

      TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitSys_WOXErr = (TGraphAsymmErrors*) graphRatioPi0TsallisCombCombFitSys->Clone("graphRatioPi0TsallisCombCombFitSys_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphRatioPi0TsallisCombCombFitSys_WOXErr);
      TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0TsallisCombCombFitStat->Clone("graphRatioPi0TsallisCombCombFitStat_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphRatioPi0TsallisCombCombFitStat_WOXErr);
      TGraphAsymmErrors* graphRatioPi0TsallisCombCombFitTot_WOXErr = (TGraphAsymmErrors*) graphRatioPi0TsallisCombCombFitTot->Clone("graphRatioPi0TsallisCombCombFitTot_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphRatioPi0TsallisCombCombFitTot_WOXErr);

      TGraphAsymmErrors* graphRatioPi0TsallisCombFitSys_WOXErr[11];
      TGraphAsymmErrors* graphRatioPi0TsallisCombFitStat_WOXErr[11];
      for (Int_t i = 0; i < 11; i++){
         if(directoryPi0[i]){
           graphRatioPi0TsallisCombFitSys_WOXErr[i]                    = (TGraphAsymmErrors*)graphRatioPi0TsallisCombFitSys[i]->Clone();
           ProduceGraphAsymmWithoutXErrors(graphRatioPi0TsallisCombFitSys_WOXErr[i]);
           graphRatioPi0TsallisCombFitStat_WOXErr[i]                   = (TGraphAsymmErrors*)graphRatioPi0TsallisCombFitStat[i]->Clone();
           ProduceGraphAsymmWithoutXErrors(graphRatioPi0TsallisCombFitStat_WOXErr[i]);
         }
      }

      TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitTot     = (TGraphAsymmErrors*)graphCombEtaInvXSectionTot->Clone();
      graphRatioEtaTsallisCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombCombFitTot, fitInvXSectionEta);
      TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitStat    = (TGraphAsymmErrors*)graphCombEtaInvXSectionStat->Clone();
      graphRatioEtaTsallisCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombCombFitStat, fitInvXSectionEta);
      TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitSys     = (TGraphAsymmErrors*)graphCombEtaInvXSectionSys->Clone();
      graphRatioEtaTsallisCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombCombFitSys, fitInvXSectionEta);

      TGraphAsymmErrors* graphRatioEtaTsallisCombFitSys[11];
      TGraphAsymmErrors* graphRatioEtaTsallisCombFitStat[11];
      for (Int_t i = 0; i < 11; i++){
         if(directoryEta[i]){
           graphRatioEtaTsallisCombFitSys[i]                    = (TGraphAsymmErrors*)graphEtaInvXSectionSys[i]->Clone();
           graphRatioEtaTsallisCombFitSys[i]                    = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombFitSys[i], fitInvXSectionEta);
           graphRatioEtaTsallisCombFitStat[i]                   = (TGraphAsymmErrors*)graphEtaInvXSectionStat[i]->Clone();
           graphRatioEtaTsallisCombFitStat[i]                   = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombFitStat[i], fitInvXSectionEta);
         }
      }

      TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitSys_WOXErr = (TGraphAsymmErrors*) graphRatioEtaTsallisCombCombFitSys->Clone("graphRatioEtaTsallisCombCombFitSys_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisCombCombFitSys_WOXErr);
      TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaTsallisCombCombFitStat->Clone("graphRatioEtaTsallisCombCombFitStat_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisCombCombFitStat_WOXErr);
      TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitTot_WOXErr = (TGraphAsymmErrors*) graphRatioEtaTsallisCombCombFitTot->Clone("graphRatioEtaTsallisCombCombFitTot_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisCombCombFitTot_WOXErr);

      TGraphAsymmErrors* graphRatioEtaTsallisCombFitSys_WOXErr[11];
      TGraphAsymmErrors* graphRatioEtaTsallisCombFitStat_WOXErr[11];
      for (Int_t i = 0; i < 11; i++){
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
      histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(minPtPi0,maxPtPi0);
      histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.3);
      histo2DPi0RatioToCombFit->SetYTitle("Data/Tsallis fit");
      histo2DPi0RatioToCombFit->Draw("copy");

          DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisCombCombFitSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
          graphRatioPi0TsallisCombCombFitSys->Draw("E2same");
          DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
          graphRatioPi0TsallisCombCombFitStat_WOXErr->Draw("p,same,z");

          DrawGammaLines(minPtPi0,maxPtPi0 , 1., 1.,0.1, kGray+2);
          DrawGammaLines(minPtPi0,maxPtPi0 , 1.1, 1.1,0.1, kGray, 7);
          DrawGammaLines(minPtPi0,maxPtPi0 , 0.9, 0.9,0.1, kGray, 7);

          labelRatioToFitEnergy->Draw();
          labelRatioToFitALICE->Draw();
          labelRatioToFitPi0->Draw();

      canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombTsallisFit_PP900GeV.%s",outputDir.Data(),suffix.Data()));

      // **********************************************************************************************************************
      // *******************************************Plot Ratio of Individual meas to Fit ******************************************
      // **********************************************************************************************************************

      canvasRatioToCombFit->cd();
      histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,3.05);
      histo2DPi0RatioToCombFit->Draw("copy");

      for (Int_t i = 0; i < 11; i++){
        if(directoryPi0[i]){
          DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
          DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisCombFitStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);

          graphRatioPi0TsallisCombFitSys[i]->Draw("E2same");
          graphRatioPi0TsallisCombFitStat[i]->Draw("p,same,e");
        }
      }

      DrawGammaLines(minPtPi0,maxPtPi0 , 1., 1.,0.5, kGray+2);
      DrawGammaLines(minPtPi0,maxPtPi0 , 1.1, 1.1,0.5, kGray, 7);
      DrawGammaLines(minPtPi0,maxPtPi0 , 0.9, 0.9,0.5, kGray, 7);

      labelRatioToFitEnergy->Draw();
      labelRatioToFitALICE->Draw();
      labelRatioToFitPi0->Draw();

      //****************************** Definition of the Legend ******************************************
      //**************** Row def ************************

      //****************** first Column **************************************************
      for (Int_t i = 0; i < 11; i++){
        if(directoryPi0[i]){
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
        if(directoryPi0[i]){
          markerPi0OnlyRatio[i]               = CreateMarkerFromGraph(graphRatioPi0TsallisCombFitSys[i],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[i+1],1);
          markerPi0OnlyRatio[i]->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[i+1]);
          boxPi0OnlyRatio[i]                  = CreateBoxFromGraph(graphRatioPi0TsallisCombFitSys[i], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[i+1]- heightBox,
          columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox+0.1, rowsLegendOnlyPi0RatioAbs[i+1]+ heightBox);
          boxPi0OnlyRatio[i]->Draw("l");
        }
      }

      canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToTsallisFit_PP900GeV.%s",outputDir.Data(),suffix.Data()));
      histo2DPi0RatioToCombFit->SetYTitle("Data/TCM fit");

      // **********************************************************************************************************************
      // ******************************************* Ratio of Individual meas to Fit ******************************************
      // **********************************************************************************************************************

      canvasRatioToCombFit->cd();
      histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
      histo2DEtaRatioToCombFit->GetXaxis()->SetRangeUser(minPtEta,maxPtEta);
      histo2DEtaRatioToCombFit->Draw("copy");

      for (Int_t i = 0; i < 11; i++){
        if(directoryEta[i]){
          DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
          DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisCombFitStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);

          graphRatioEtaTsallisCombFitSys[i]->Draw("E2same");
          graphRatioEtaTsallisCombFitStat[i]->Draw("p,same,e");
        }
      }

      DrawGammaLines(minPtEta,maxPtEta , 1., 1.,0.5, kGray+2);
      DrawGammaLines(minPtEta,maxPtEta , 1.1, 1.1,0.5, kGray, 7);
      DrawGammaLines(minPtEta,maxPtEta , 0.9, 0.9,0.5, kGray, 7);

      labelRatioToFitEnergy->Draw();
      labelRatioToFitALICE2->Draw();
      labelRatioToFitEta->Draw();

      //****************************** Definition of the Legend ******************************************
      //**************** Row def ************************

      //****************** first Column **************************************************
      for (Int_t i = 0; i < 11; i++){
        if(directoryEta[i]){
          textSingleMeasRatioEta[i]           = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[i+1],nameMeasGlobal[i].Data());
          SetStyleTLatex( textSingleMeasRatioEta[i], 0.85*textSizeLabelsPixel,4);
          textSingleMeasRatioEta[i]->SetTextFont(43);
          textSingleMeasRatioEta[i]->Draw();
        }
      }

      //****************** second Column *************************************************
      textStatOnlyRatioEta->Draw();
      textSysOnlyRatioEta->Draw();

      for (Int_t i = 0; i < 11; i++){
        if(directoryEta[i]){
          markerEtaOnlyRatio[i]               = CreateMarkerFromGraph(graphRatioEtaTsallisCombFitSys[i],columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[i+1],1);
          markerEtaOnlyRatio[i]->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[i+1]);
          boxEtaOnlyRatio[i]                  = CreateBoxFromGraph(graphRatioEtaTsallisCombFitSys[i], columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[i+1]- heightBox,
          columnsLegendOnlyEtaRatioAbs[2]+ 3*lengthBox+0.1, rowsLegendOnlyEtaRatioAbs[i+1]+ heightBox);
          boxEtaOnlyRatio[i]->Draw("l");
        }
      }

      canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToTsallisFit_PP900GeV.%s",outputDir.Data(),suffix.Data()));


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
          TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystem900GeV.Data());
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

          TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, minPtPi0,maxPtPi0, 1000., 120.1, 152.9);
          SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                                    textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
          histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
          histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
          histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
          histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
          histo2DAllPi0Mass->GetXaxis()->SetNoExponent();
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
          Double_t  rowsLegendMass2[7]= {0.84,0.66,0.66,0.50,0.01,0.33};
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

          TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, minPtEta,maxPtEta ,1000., -4, 55.5);
          SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                    0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
          histo2DAllEtaFWHM->GetYaxis()->SetRangeUser(-4.,45.5);
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

          TH2F * histo2DAllEtaMass            = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",20, minPtEta,maxPtEta, 1000., 505.1, 599.9);
          SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMassEta, textsizeLabelsMassEta, 0.85*textsizeLabelsMassEta,
                                    textsizeLabelsMassEta, 0.9, 0.28/(textsizeFacMassEta*margin), 512, 505);
          histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
          histo2DAllEtaMass->GetYaxis()->SetNdivisions(505);
          histo2DAllEtaMass->GetYaxis()->SetNoExponent(kTRUE);
          histo2DAllEtaMass->GetXaxis()->SetTickLength(0.05);
          histo2DAllEtaMass->GetXaxis()->SetNoExponent();
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
          histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, minPtPi0, maxPtPi0, 1000, 3e-5, 2 );
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
          drawLatexAdd(collisionSystem900GeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
          drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

      canvasAcceptanceTimesEff->Update();
      canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));

      // **********************************************************************************************************************
      // ******************************** Cross section for pi0 single measurement 900GeV ************************************
      // **********************************************************************************************************************

      TCanvas* canvasXSectionPi0  = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
      DrawGammaCanvasSettings( canvasXSectionPi0, 0.14, 0.02, 0.02, 0.09);
      canvasXSectionPi0->SetLogx();
      canvasXSectionPi0->SetLogy();

      TH2F * histo2DXSectionPi0;
      histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,minPtPi0,maxPtPi0,1000,minXSectionPi0,maxXSectionPi0);
      SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
      histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
      histo2DXSectionPi0->GetXaxis()->SetNoExponent(kTRUE);
      histo2DXSectionPi0->Draw("copy");

      for (Int_t i = 0; i < 11; i++){
        if(directoryPi0[i]){
          DrawGammaSetMarkerTGraphAsym(graphPi0InvXSectionStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
          graphPi0InvXSectionStat[i]->Draw("pEsame");
          DrawGammaSetMarkerTGraphAsym(graphPi0InvXSectionSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
          graphPi0InvXSectionSys[i]->Draw("E2same");
        }
      }

      TLatex *labelEnergyXSectionPi0      = new TLatex(0.64,0.92,collisionSystem900GeV.Data());
      SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4);
      labelEnergyXSectionPi0->Draw();
      TLatex *labelDetSysXSectionPi0      = new TLatex(0.64,0.88,"#pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4);
      labelDetSysXSectionPi0->Draw();

      TLegend* legendXSectionPi0          = new TLegend(0.62,0.62,0.9,0.86);
      legendXSectionPi0->SetFillColor(0);
      legendXSectionPi0->SetLineColor(0);
      legendXSectionPi0->SetTextFont(42);
      legendXSectionPi0->SetTextSize(0.035);
      legendXSectionPi0->AddEntry(graphPi0InvXSectionSys[0],"PCM","fp");
      legendXSectionPi0->AddEntry(graphPi0InvXSectionSys[1],"PHOS","fp");
      legendXSectionPi0->AddEntry(graphPi0InvXSectionSys[2],"EMC","fp");
      legendXSectionPi0->AddEntry(graphPi0InvXSectionSys[4],"PCM-EMC","fp");
      legendXSectionPi0->Draw();

      canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvXSectionCompAllSystems.%s",outputDir.Data(),suffix.Data()));

      canvasXSectionPi0->cd();
      histo2DXSectionPi0->Draw("copy");

      for (Int_t i = 0; i < 11; i++){
        if(directoryPi0[i]){
          graphPi0InvXSectionStat[i]->Draw("pEsame");
          graphPi0InvXSectionSys[i]->Draw("E2same");
        }
      }
      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
      graphCombPi0InvXSectionSys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStat_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
      graphCombPi0InvXSectionStat_WOXErr->Draw("p,same,z");

      labelEnergyXSectionPi0->Draw();
      labelDetSysXSectionPi0->Draw();

      legendXSectionPi0->AddEntry(graphCombPi0InvXSectionSys,"comb","fp");
      legendXSectionPi0->Draw();

      canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvXSectionCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));



      TCanvas* canvasXSectionEta      = new TCanvas("canvasXSectionEta","",200,10,1350,1350*1.15);  // gives the page size
      DrawGammaCanvasSettings( canvasXSectionEta, 0.14, 0.02, 0.02, 0.09);
      canvasXSectionEta->SetLogx();
      canvasXSectionEta->SetLogy();

      TH2F * histo2DXSectionEta;
      histo2DXSectionEta              = new TH2F("histo2DXSectionEta","histo2DXSectionEta",11000,minPtEta, maxPtEta,1000,minXSectionEta,maxXSectionEta);
      SetStyleHistoTH2ForGraphs(histo2DXSectionEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 1.,1.45);
      histo2DXSectionEta->GetXaxis()->SetMoreLogLabels();
      histo2DXSectionEta->GetXaxis()->SetNoExponent(kTRUE);
      histo2DXSectionEta->Draw("copy");


      for (Int_t i = 0; i < 11; i++){
        if(directoryEta[i]){
          DrawGammaSetMarkerTGraphAsym(graphEtaInvXSectionStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i]);
          graphEtaInvXSectionStat[i]->Draw("pEsame");
          DrawGammaSetMarkerTGraphAsym(graphEtaInvXSectionSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.75, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
          graphEtaInvXSectionSys[i]->Draw("E2same");
        }
      }

      TLatex *labelEnergyXSectionEta      = new TLatex(0.64,0.92,collisionSystem900GeV.Data());
      SetStyleTLatex( labelEnergyXSectionEta, 0.035,4);
      labelEnergyXSectionEta->Draw();
      TLatex *labelDetSysXSectionEta      = new TLatex(0.64,0.88,"#eta #rightarrow #gamma#gamma");
      SetStyleTLatex( labelDetSysXSectionEta, 0.035,4);
      labelDetSysXSectionEta->Draw();

      TLegend* legendXSectionEta          = new TLegend(0.62,0.8,0.9,0.86);
      legendXSectionEta->SetFillColor(0);
      legendXSectionEta->SetLineColor(0);
      legendXSectionEta->SetTextFont(42);
      legendXSectionEta->SetTextSize(0.035);
      legendXSectionEta->AddEntry(graphEtaInvXSectionSys[0],"PCM","fp");
      legendXSectionEta->Draw();

      canvasXSectionEta->SaveAs(Form("%s/Eta_InvXSectionCompAllSystems.%s",outputDir.Data(),suffix.Data()));
      histo2DXSectionEta->Draw("copy");

      for (Int_t i = 0; i < 11; i++){
        if(directoryEta[i]){
          graphEtaInvXSectionStat[i]->Draw("pEsame");
          graphEtaInvXSectionSys[i]->Draw("E2same");
        }
      }

      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
      graphCombEtaInvXSectionSys->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStat_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
      graphCombEtaInvXSectionStat_WOXErr->Draw("p,same,z");


      labelEnergyXSectionEta->Draw();
      labelDetSysXSectionEta->Draw();

      legendXSectionEta->AddEntry(graphCombEtaInvXSectionSys,"comb","fp");
      legendXSectionEta->Draw();

     canvasXSectionEta->SaveAs(Form("%s/Eta_InvXSectionCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));


     // **********************************************************************************************************************
     // ******************************** effective secondary correction drawing for different methods ************************
     // **********************************************************************************************************************
     TCanvas* canvasEffectiveSecCorr       = new TCanvas("canvasEffectiveSecCorr", "", 200, 10, 1200, 1100);  // gives the page size
     DrawGammaCanvasSettings( canvasEffectiveSecCorr,  0.1, 0.01, 0.04, 0.095);
     canvasEffectiveSecCorr->SetLogx(1);

         TH1F * histo1DEffSecCorr;
         histo1DEffSecCorr                = new TH1F("histo1DEffSecCorr", "histo1DEffSecCorr",1000, minPtPi0, maxPtPi0);
         SetStyleHistoTH1ForGraphs( histo1DEffSecCorr, "#it{p}_{T} (GeV/#it{c})","R_{K^{0}_{s}}",
                                 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
         histo1DEffSecCorr->GetYaxis()->SetRangeUser(0, 10 );
         histo1DEffSecCorr->GetXaxis()->SetLabelOffset(-0.01);
         histo1DEffSecCorr->GetXaxis()->SetMoreLogLabels(kTRUE);

         for (Int_t k = 0; k < 4; k++){
             Bool_t plotCorr     = kFALSE;
             Int_t nCorrAvail    = 0;
             for (Int_t i = 0; i < 11; i++){
                 if (haveEffSecCorr[k][i]){
                     nCorrAvail++;
                     plotCorr    = kTRUE;
                 }
             }
             TLegend* legendEffSecCorrPi0           = GetAndSetLegend2(0.65, 0.925-(nCorrAvail*textSizeLabelsRel), 0.93, 0.925,textSizeLabelsPixel);
             if (plotCorr){
                 histo1DEffSecCorr->GetYaxis()->SetTitle(Form("R_{%s}",nameSecPi0SourceLabel[k].Data()));
                 histo1DEffSecCorr->GetYaxis()->SetRangeUser(0, maxSecCorr[k]);
                 histo1DEffSecCorr->DrawCopy();
                 for (Int_t i = 0; i < 11; i++){
                     if (haveEffSecCorr[k][i]){
                         DrawGammaSetMarkerTGraphAsym(graphPi0EffSecCorrFromX[k][i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
                         graphPi0EffSecCorrFromX[k][i]->Draw("p,same,z");
                         legendEffSecCorrPi0->AddEntry(graphPi0EffSecCorrFromX[k][i],nameMeasGlobal[i],"p");
                     }
                 }
                 legendEffSecCorrPi0->Draw();
                 Double_t xAlign = 0.15;
                 if(k==2) xAlign += 0.05;
                 TLatex *labelPerfSecCorr               = new TLatex(xAlign,0.90,"ALICE performance");
                 SetStyleTLatex( labelPerfSecCorr, textSizeLabelsRel,4);
                 labelPerfSecCorr->Draw();
                 TLatex *labelEnergySecCorr             = new TLatex(xAlign,0.85,collisionSystem900GeV.Data());
                 SetStyleTLatex( labelEnergySecCorr, textSizeLabelsRel,4);
                 labelEnergySecCorr->Draw();
                 TLatex *labelDetSysSecCorrPi0          = new TLatex(xAlign,0.80,"#pi^{0} #rightarrow #gamma#gamma");
                 SetStyleTLatex( labelDetSysSecCorrPi0, textSizeLabelsRel,4);
                 labelDetSysSecCorrPi0->Draw();

                 canvasEffectiveSecCorr->Update();
                 canvasEffectiveSecCorr->Print(Form("%s/Pi0_EffectiveSecCorr_%s.%s",outputDir.Data(), nameSecPi0SourceRead[k].Data() , suffix.Data()));
             }
         }

     delete canvasEffectiveSecCorr;


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

             TH2F * ratio2DTheoryPP       = new TH2F("ratio2DTheoryPP","ratio2DTheoryPP",1000,minPtPi0,maxPtPi0,1000,0.4,5.1);
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

             graphRatioPi0DSS07->SetLineWidth(widthCommonFit);
             graphRatioPi0DSS07->SetLineColor(colorNLO);
             graphRatioPi0DSS07->SetLineStyle(1);
             graphRatioPi0DSS07->SetFillStyle(1001);
             graphRatioPi0DSS07->SetFillColor(colorNLO);
             graphRatioPi0DSS07->Draw("same,e4");

             DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
             graphRatioPythia8ToFit->Draw("3,same");
             DrawGammaSetMarker(histoRatioPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);
             histoRatioPythia8ToFit->SetLineWidth(widthCommonFit);
             histoRatioPythia8ToFit->Draw("same,hist,l");

             DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
             graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
             DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
             graphRatioCombCombFitSys->SetLineWidth(0);
             graphRatioCombCombFitSys->Draw("2,same");
             graphRatioCombCombFitStat_WOXErr->Draw("p,same");

             TBox* boxErrorSigmaRatio = CreateBoxConv(kGray+2, 0.35, 1.-(0.035 ), 0.4, 1.+(0.035));
             boxErrorSigmaRatio->SetLineWidth(8);
             boxErrorSigmaRatio->Draw();
             DrawGammaLines(minPtPi0, maxPtPi0,1., 1.,0.5,kGray+2);

             TLegend* legendRatioTheorypp_3Parted = GetAndSetLegend2(0.2,0.82,0.45,0.96, 0.85* textSizeLabelsPixel);
             legendRatioTheorypp_3Parted->AddEntry(graphRatioCombCombFitSys,"Data","pf");
             legendRatioTheorypp_3Parted->AddEntry(histoRatioPythia8ToFit,  "PYTHIA 8.2, Monash 2013", "l");
             legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0DSS07,  "NLO, PDF:MSTW08 - FF:DSS07", "f");
             legendRatioTheorypp_3Parted->Draw();

             TLegend* legendRatioTheoryNormUnc = GetAndSetLegend2(0.34,0.91,0.59,0.96, 0.85* textSizeLabelsPixel);
             legendRatioTheoryNormUnc->AddEntry(boxErrorSigmaRatio,"norm. unc. 3.5%","l");
             legendRatioTheoryNormUnc->Draw();

             TLatex *labelRatioTheoryPPP   = new TLatex(0.268,0.79,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
             SetStyleTLatex( labelRatioTheoryPPP, 0.85*textsizeLabelsPP,4);
             labelRatioTheoryPPP->Draw();

             TLatex *labelRatioTheoryPP   = new TLatex(0.72,0.925,collisionSystem900GeV.Data());
             SetStyleTLatex( labelRatioTheoryPP, 0.85*textsizeLabelsPP,4);
             labelRatioTheoryPP->Draw();
             TLatex *labelRatioTheoryPP1P = new TLatex(0.863,0.875,"ALICE");
             SetStyleTLatex( labelRatioTheoryPP1P, 0.85*textsizeLabelsPP,4);
             labelRatioTheoryPP1P->Draw();
             TLatex *labelRatioTheoryPP2P= new TLatex(0.843,0.83,"#pi^{0} #rightarrow #gamma#gamma");
             SetStyleTLatex( labelRatioTheoryPP2P, 0.85*textsizeLabelsPP,4);
             labelRatioTheoryPP2P->Draw();


         canvasRatioPP->Update();
         canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP.%s",outputDir.Data(),suffix.Data()));


         TH2F * ratio2DTheoryPP2       = new TH2F("ratio2DTheoryPP2","ratio2DTheoryPP2",1000,minPtPi0,maxPtPi0,1000,0.4,5.1);
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

         DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, kGray+1);
         graphRatioPi0CombNLOMuHalf->Draw("same,c");
         DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, kGray+1);
         graphRatioPi0CombNLOMuOne->Draw("same,c");
         DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, kGray+2);
         graphRatioPi0CombNLOMuTwo->Draw("same,c");

         DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
         graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
         DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
         graphRatioCombCombFitSys->SetLineWidth(0);
         graphRatioCombCombFitSys->Draw("2,same");
         graphRatioCombCombFitStat_WOXErr->Draw("p,same");

         boxErrorSigmaRatio->Draw();
         DrawGammaLines(minPtPi0, maxPtPi0,1., 1.,0.5,kGray+2);

         TLegend* legendRatioTheorypp_3Parted2= GetAndSetLegend2(0.15,0.65,0.4,0.96, 0.85* textSizeLabelsPixel);
         legendRatioTheorypp_3Parted2->AddEntry(graphRatioCombCombFitSys,"Data","pf");
         legendRatioTheorypp_3Parted2->AddEntry((TObject*)0,"NLO, PDF:CTEQ6M5 - FF:DSS07", "");
         legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuHalf, "#mu = 0.5 #it{p}_{T}", "l");
         legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuOne,  "#mu = #it{p}_{T}", "l");
         legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuTwo,  "#mu = 2 #it{p}_{T}", "l");
         legendRatioTheorypp_3Parted2->Draw();

         TLegend* legendRatioTheoryNormUnc2 = GetAndSetLegend2(0.34,0.902,0.59,0.952, 0.85* textSizeLabelsPixel);
         legendRatioTheoryNormUnc2->AddEntry(boxErrorSigmaRatio,"norm. unc. 3.5%","l");
         legendRatioTheoryNormUnc2->Draw();

         labelRatioTheoryPP->Draw();
         labelRatioTheoryPP1P->Draw();
         labelRatioTheoryPP2P->Draw();

         canvasRatioPP->Update();
         canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP2.%s",outputDir.Data(),suffix.Data()));


         // **********************************************************************************************************************
         // ******************************************* Comparison to theory calculations Eta ************************************
         // **********************************************************************************************************************
         canvasRatioPP->cd();
         canvasRatioPP->SetLogx();

             TH2F * ratio2DTheoryPPEta       = new TH2F("ratio2DTheoryPPEta","ratio2DTheoryPPEta",1000,minPtEta,maxPtEta,1000,0.2,3.9);
             SetStyleHistoTH2ForGraphs(ratio2DTheoryPPEta, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{Tsallis fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
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

             DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFitEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
             graphRatioPythia8ToFitEta->Draw("3,same");
             DrawGammaSetMarker(histoRatioPythia8ToFitEta, 24, 1.5, kRed+2 , kRed+2);
             histoRatioPythia8ToFitEta->SetLineWidth(widthCommonFit);
             histoRatioPythia8ToFitEta->Draw("same,hist,l");

             DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatEta_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
             graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
             DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysEta, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
             graphRatioCombCombFitSysEta->SetLineWidth(0);
             graphRatioCombCombFitSysEta->Draw("2,same");
             graphRatioCombCombFitStatEta_WOXErr->Draw("p,same");

             TBox* boxErrorSigmaRatioEta = CreateBoxConv(kGray+2, 0.35, 1.-(0.035 ), 0.4, 1.+(0.035));
             boxErrorSigmaRatioEta->SetLineWidth(8);
             boxErrorSigmaRatioEta->Draw();
             DrawGammaLines(minPtEta, maxPtEta, 1., 1.,0.5, kGray+2);

             TLegend* legendRatioTheoryppEta_3Parted= GetAndSetLegend2(0.15,0.84-(0.85*textsizeLabelsPP*5),0.40,0.96, 0.85* textSizeLabelsPixel);
             legendRatioTheoryppEta_3Parted->AddEntry(graphRatioCombCombFitSysEta,"Data","pf");
             legendRatioTheoryppEta_3Parted->AddEntry(histoRatioPythia8ToFitEta,  "PYTHIA 8.2, Monash 2013", "l");
             legendRatioTheoryppEta_3Parted->AddEntry((TObject*)0, "NLO, PDF:CTEQ6M5 - FF:AESSS", "");
             legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuHalf, "#mu = 0.5 #it{p}_{T}", "l");
             legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuOne,  "#mu = #it{p}_{T}", "l");
             legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuTwo,  "#mu = 2 #it{p}_{T}", "l");
             legendRatioTheoryppEta_3Parted->Draw();

             TLegend* legendRatioTheoryNormUncEta = GetAndSetLegend2(0.31,0.9,0.56,0.96, 0.85* textSizeLabelsPixel);
             legendRatioTheoryNormUncEta->AddEntry(boxErrorSigmaRatio,"norm. unc. 3.5%","l");
             legendRatioTheoryNormUncEta->Draw();

             TLatex *labelRatioTheoryPP2   = new TLatex(0.75,0.925,collisionSystem900GeV.Data());
             SetStyleTLatex( labelRatioTheoryPP2, 0.85*textsizeLabelsPP,4);
             labelRatioTheoryPP2->Draw();
             TLatex *labelRatioTheoryPP221 = new TLatex(0.883,0.875,"ALICE");
             SetStyleTLatex( labelRatioTheoryPP221, 0.85*textsizeLabelsPP,4);
             labelRatioTheoryPP221->Draw();
             TLatex *labelRatioTheoryPP222= new TLatex(0.873,0.83,"#eta #rightarrow #gamma#gamma");
             SetStyleTLatex( labelRatioTheoryPP222, 0.85*textsizeLabelsPP,4);
             labelRatioTheoryPP222->Draw();

         canvasRatioPP->Update();
         canvasRatioPP->Print(Form("%s/Eta_RatioTheoryToData_PP.%s",outputDir.Data(),suffix.Data()));

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
             SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",
                                     0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
             histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
             histo2DXSectionPi0->GetYaxis()->SetRangeUser(minXSectionPi0,maxXSectionPi0/20.);
             histo2DXSectionPi0->GetXaxis()->SetLabelOffset(+0.01);
             histo2DXSectionPi0->Draw();

             graphPi0DSS07->SetLineWidth(widthCommonFit);
             graphPi0DSS07->SetLineColor(colorNLO);
             graphPi0DSS07->SetLineStyle(1);
             graphPi0DSS07->SetFillStyle(1001);
             graphPi0DSS07->SetFillColor(colorNLO);
             graphPi0DSS07->Draw("same,e3");

             DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
             graphPythia8InvXSection->Draw("3,same");
             DrawGammaSetMarker(histoPythia8InvXSection, 24, 1.5, kRed+2 , kRed+2);
             histoPythia8InvXSection->SetLineWidth(widthCommonFit);
             histoPythia8InvXSection->Draw("same,hist,l");

             DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0, 7, 2, kGray+2);
             fitTCMInvXSectionPi0->Draw("same");

             DrawGammaSetMarkerTF1( fitInvXSectionPi0, 3, 2, kGray+1);
             fitInvXSectionPi0->Draw("same");

             DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
             graphCombPi0InvXSectionSys->Draw("E2same");
             DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
             graphCombPi0InvXSectionStat_WOXErr->Draw("p,same,z");

             TLatex *labelEnergyXSectionPaper= new TLatex(0.678, 0.91, collisionSystem900GeV.Data());
             SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4);
             labelEnergyXSectionPaper->Draw();
             TLatex *labelALICEXSectionPaper= new TLatex(0.848,0.86,"ALICE");
             SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4);
             labelALICEXSectionPaper->Draw();
             TLatex *labelDetSysXSectionPaper= new TLatex(0.824,0.815,"#pi^{0} #rightarrow #gamma#gamma");
             SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4);
             labelDetSysXSectionPaper->Draw();

             TLegend* legendXsectionPaper    = GetAndSetLegend2(0.17, 0.08, 0.5, 0.18+0.05*6, textSizeLabelsPixel);
             legendXsectionPaper->SetNColumns(1);
             legendXsectionPaper->SetMargin(0.2);
             legendXsectionPaper->AddEntry(graphCombPi0InvXSectionSys,"Data","pf");
             legendXsectionPaper->AddEntry(boxErrorSigmaRatio,"norm. unc. 3.5%","l");
             legendXsectionPaper->AddEntry(fitTCMInvXSectionPi0,"TCM fit","l");
             legendXsectionPaper->AddEntry(fitInvXSectionPi0,"Tsallis fit","l");
             legendXsectionPaper->AddEntry(histoPythia8InvXSection,"PYTHIA 8.2, Monash 2013","l");
             legendXsectionPaper->AddEntry(graphPi0DSS07,  "NLO, PDF:MSTW08 - FF:DSS07", "f");
             legendXsectionPaper->Draw();

             TLatex *labelRatioTheoryPP_Paper   = new TLatex(0.24,0.055,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
             SetStyleTLatex( labelRatioTheoryPP_Paper, 0.8*textsizeLabelsPP,4);
             labelRatioTheoryPP_Paper->Draw();

         padInvSectionNLORatio->cd();
         padInvSectionNLORatio->SetLogx(1);
             TH2F * ratio2DNLO               = new TH2F("ratio2DNLO","ratio2DNLO",1000,minPtPi0,maxPtPi0,1000,0.4,4.95);
             SetStyleHistoTH2ForGraphs(ratio2DNLO, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{TCM fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                       0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
             ratio2DNLO->GetYaxis()->SetMoreLogLabels(kTRUE);
             ratio2DNLO->GetYaxis()->SetNdivisions(505);
             ratio2DNLO->GetYaxis()->SetNoExponent(kTRUE);
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

     //        graphRatioPi0CombNLODSS07->RemovePoint(0);
     //        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS07, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
     //        graphRatioPi0CombNLODSS07->Draw("3,same");
     //        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
     //        graphRatioPi0CombNLOMuHalf->Draw("same,c");
     //        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
     //        graphRatioPi0CombNLOMuOne->Draw("same,c");
     //        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
     //        graphRatioPi0CombNLOMuTwo->Draw("same,c");
             graphRatioPi0DSS07->SetLineWidth(widthCommonFit);
             graphRatioPi0DSS07->SetLineColor(colorNLO);
             graphRatioPi0DSS07->SetLineStyle(1);
             graphRatioPi0DSS07->SetFillStyle(1001);
             graphRatioPi0DSS07->SetFillColor(colorNLO);
             graphRatioPi0DSS07->Draw("same,e4");

             DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
             graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
             DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
             graphRatioCombCombFitSys->SetLineWidth(0);
             graphRatioCombCombFitSys->Draw("2,same");
             graphRatioCombCombFitStat_WOXErr->Draw("p,same");

             boxErrorSigmaRatio->Draw();
             DrawGammaLines(minPtPi0, maxPtPi0, 1., 1.,0.5, kGray+2);

         padInvSectionPythiaRatio->cd();
         padInvSectionPythiaRatio->SetLogx(1);
             TH2F * ratio2DPythia            = new TH2F("ratio2DPythia","ratio2DPythia",1000,minPtPi0,maxPtPi0,1000,0.4,2.45);
             SetStyleHistoTH2ForGraphs(ratio2DPythia, "#it{p}_{T} (GeV/#it{c})","#frac{Pythia, Data}{TCM fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                       0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 505);
             ratio2DPythia->GetYaxis()->SetMoreLogLabels(kTRUE);
             ratio2DPythia->GetYaxis()->SetNdivisions(505);
             ratio2DPythia->GetYaxis()->SetNoExponent(kTRUE);
             ratio2DPythia->GetXaxis()->SetMoreLogLabels(kTRUE);
             ratio2DPythia->GetXaxis()->SetNoExponent(kTRUE);
             ratio2DPythia->GetXaxis()->SetLabelFont(42);
             ratio2DPythia->GetYaxis()->SetLabelFont(42);
             ratio2DPythia->GetYaxis()->SetLabelOffset(+0.01);
             ratio2DPythia->GetXaxis()->SetTickLength(0.06);
             ratio2DPythia->GetYaxis()->SetTickLength(0.04);
             ratio2DPythia->DrawCopy();

             DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
             graphRatioPythia8ToFit->Draw("3,same");
             DrawGammaSetMarker(histoRatioPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);
             histoRatioPythia8ToFit->SetLineWidth(widthCommonFit);
             histoRatioPythia8ToFit->Draw("same,hist,l");

             DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
             graphRatioCombCombFitStat_WOXErr->SetLineWidth(widthLinesBoxes);
             DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
             graphRatioCombCombFitSys->SetLineWidth(0);
             graphRatioCombCombFitSys->Draw("2,same");
             graphRatioCombCombFitStat_WOXErr->Draw("p,same");

             boxErrorSigmaRatio->Draw();
             DrawGammaLines(minPtPi0, maxPtPi0 , 1., 1.,0.5, kGray+2);

         canvasInvSectionPaper->Print(Form("%s/Pi0_InvXSectionWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));

         padInvSectionSpec->cd();
         padInvSectionSpec->SetLogy(1);
         padInvSectionSpec->SetLogx(1);
             SetStyleHistoTH2ForGraphs(histo2DXSectionEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",
                                     0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
             histo2DXSectionEta->GetXaxis()->SetMoreLogLabels();
             histo2DXSectionEta->GetXaxis()->SetLabelOffset(+0.01);
             histo2DXSectionEta->Draw();

             DrawGammaNLOTGraph( graphNLOCalcEtaMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
             graphNLOCalcEtaMuHalf->Draw("same,c");
             DrawGammaNLOTGraph( graphNLOCalcEtaMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
             graphNLOCalcEtaMuOne->Draw("same,c");
             DrawGammaNLOTGraph( graphNLOCalcEtaMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
             graphNLOCalcEtaMuTwo->Draw("same,c");

             DrawGammaSetMarkerTGraphErr(graphPythia8InvXSectionEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
             graphPythia8InvXSectionEta->Draw("3,same");
             DrawGammaSetMarker(histoPythia8InvXSectionEta, 24, 1.5, kRed+2 , kRed+2);
             histoPythia8InvXSectionEta->SetLineWidth(widthCommonFit);
             histoPythia8InvXSectionEta->Draw("same,hist,l");

             DrawGammaSetMarkerTF1( fitInvXSectionEta, 7, 2, kGray+2);
             fitInvXSectionEta->Draw("same");

             DrawGammaSetMarkerTF1( fitInvXSectionEta, 3, 2, kGray+1);
             fitInvXSectionEta->Draw("same");

             DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
             graphCombEtaInvXSectionSys->Draw("E2same");
             DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
             graphCombEtaInvXSectionStat_WOXErr->Draw("p,same,z");

             labelEnergyXSectionPaper->Draw();
             labelALICEXSectionPaper->Draw();
             TLatex *labelDetSysXSectionPaperEta = new TLatex(0.84,0.815,"#eta #rightarrow #gamma#gamma");
             SetStyleTLatex( labelDetSysXSectionPaperEta, textsizeLabelsXSecUp,4);
             labelDetSysXSectionPaperEta->Draw();

             TLegend* legendXsectionPaperEta     = GetAndSetLegend2(0.17, 0.03, 0.5, 0.13+0.05*6, textSizeLabelsPixel);
             legendXsectionPaperEta->SetNColumns(1);
             legendXsectionPaperEta->SetMargin(0.2);
             legendXsectionPaperEta->AddEntry(graphCombPi0InvXSectionSys,"Data","pf");
             legendXsectionPaperEta->AddEntry(boxErrorSigmaRatio, "norm. unc. 3.5%", "l");
             legendXsectionPaperEta->AddEntry(fitInvXSectionEta,"Tsallis fit","l");
             legendXsectionPaperEta->AddEntry(histoPythia8InvXSectionEta,"PYTHIA 8.2, Monash 2013","l");
             legendXsectionPaperEta->AddEntry((TObject*)0, "", "");
             legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuHalf,"#mu = 0.5 #it{p}_{T}","l");
             legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuOne,"#mu = #it{p}_{T}","l");
             legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuTwo,"#mu = 2 #it{p}_{T}","l");
             legendXsectionPaperEta->Draw();

             TLatex *labelEta = new TLatex(0.17, 0.19,"NLO, PDF:CTEQ6M5 - FF:AESSS");
             SetStyleTLatex( labelEta, 0.75*textsizeLabelsPP,4);
             labelEta->Draw();

         padInvSectionNLORatio->cd();
         padInvSectionNLORatio->SetLogx(1);
             TH2F * ratio2DNLOEta                = new TH2F("ratio2DNLOEta","ratio2DNLOEta",1000,minPtEta,maxPtEta,1000,0.35,3.65);
             SetStyleHistoTH2ForGraphs(ratio2DNLOEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{Tsallis fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
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
             graphRatioCombCombFitSysEta->SetLineWidth(0);
             graphRatioCombCombFitSysEta->Draw("2,same");
             graphRatioCombCombFitStatEta_WOXErr->Draw("p,same");

             boxErrorSigmaRatioEta->Draw();
             DrawGammaLines(minPtEta, maxPtEta, 1., 1.,0.5, kGray+2);

         padInvSectionPythiaRatio->cd();
         padInvSectionPythiaRatio->SetLogx(1);
             TH2F * ratio2DPythiaEta             = new TH2F("ratio2DPythiaEta","ratio2DPythiaEta",1000,minPtEta,maxPtEta,1000,0.4,1.95);
             SetStyleHistoTH2ForGraphs(ratio2DPythiaEta, "#it{p}_{T} (GeV/#it{c})","#frac{Pythia, Data}{Tsallis fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
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

             DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFitEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
             graphRatioPythia8ToFitEta->Draw("3,same");
             DrawGammaSetMarker(histoRatioPythia8ToFitEta, 24, 1.5, kRed+2 , kRed+2);
             histoRatioPythia8ToFitEta->SetLineWidth(widthCommonFit);
             histoRatioPythia8ToFitEta->Draw("same,hist,l");

             DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStatEta_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
             graphRatioCombCombFitStatEta_WOXErr->SetLineWidth(widthLinesBoxes);
             DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSysEta, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
             graphRatioCombCombFitSysEta->SetLineWidth(0);
             graphRatioCombCombFitSysEta->Draw("2,same");
             graphRatioCombCombFitStatEta_WOXErr->Draw("p,same");

             boxErrorSigmaRatioEta->Draw();
             DrawGammaLines(minPtEta, maxPtEta, 1., 1.,0.5, kGray+2);

         canvasInvSectionPaper->Print(Form("%s/Eta_InvXSectionWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));


         // ***************************************************************************************************************
         // ******************************* pi0+eta combined plot  ********************************************************
         // ***************************************************************************************************************

         canvasXSectionPi0->cd();

         TH2F * histo2DXSectionWithEtaAndPi0;
         histo2DXSectionWithEtaAndPi0          = new TH2F("histo2DXSectionWithEtaAndPi0","histo2DXSectionWithEtaAndPi0",11000,minPtPi0,maxPtPi0,1000,minXSectionPi0,maxXSectionPi0/10.);
         SetStyleHistoTH2ForGraphs(histo2DXSectionWithEtaAndPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
         histo2DXSectionWithEtaAndPi0->GetXaxis()->SetMoreLogLabels();
         histo2DXSectionWithEtaAndPi0->GetXaxis()->SetNoExponent();
         histo2DXSectionWithEtaAndPi0->Draw("copy");

             // scale eta graphs
             Double_t scaleFacEtaForCombPlot                              = 1e-2;
             TGraphAsymmErrors* graphCombEtaInvXSectionStat_WOXErrCopy   = (TGraphAsymmErrors*) graphCombEtaInvXSectionStat_WOXErr->Clone("graphCombEtaInvXSectionStatAWOXErrCopy");
             TGraphAsymmErrors* graphCombEtaInvXSectionSys_Copy          = (TGraphAsymmErrors*) graphCombEtaInvXSectionSys->Clone("graphCombEtaInvXSectionSysA_Copy");
             graphCombEtaInvXSectionStat_WOXErrCopy                      = ScaleGraph(graphCombEtaInvXSectionStat_WOXErrCopy,scaleFacEtaForCombPlot);
             graphCombEtaInvXSectionSys_Copy                             = ScaleGraph(graphCombEtaInvXSectionSys_Copy,scaleFacEtaForCombPlot);

             TH1D* histfitInvXSectionEta                               = (TH1D*)fitInvXSectionEta->GetHistogram();
             histfitInvXSectionEta->Scale(scaleFacEtaForCombPlot);

             TF1* tf1FitInvXSectionEta = new TF1("tf1FitInvXSectionEta","(1/100) * fitInvXSectionEta", 0.9, 3.);

             histoPythia8InvXSectionEta->Scale(scaleFacEtaForCombPlot);
             TGraphErrors* graphPythia8InvXSectionEtaScaled              = ScaleGraph(graphPythia8InvXSectionEta,scaleFacEtaForCombPlot);

             TGraphAsymmErrors* graphEtaAESSSCopy                        = (TGraphAsymmErrors*)graphEtaAESSS->Clone("graphEtaAESSSCopy");
             graphEtaAESSSCopy                                           = ScaleGraph(graphEtaAESSSCopy,scaleFacEtaForCombPlot);

             // plotting NLO calcs pi0
             graphPi0DSS07->SetLineWidth(widthCommonFit);
             graphPi0DSS07->SetLineColor(colorNLO);
             graphPi0DSS07->SetLineStyle(1);
             graphPi0DSS07->SetFillStyle(1001);
             graphPi0DSS07->SetFillColor(colorNLO);
             graphPi0DSS07->Draw("same,e3");

             // plotting NLO calcs eta
             DrawGammaSetMarkerTGraphAsym(graphEtaAESSSCopy, 0, 0, colorCGC, colorCGC, widthLinesBoxes, kTRUE, colorCGC);
             graphEtaAESSSCopy->Draw("3,same");

             // plotting Pythia 8.2 Monash
             DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
             graphPythia8InvXSection->Draw("3,same");
             DrawGammaSetMarker(histoPythia8InvXSection, 24, 1.5, kRed+2 , kRed+2);
             histoPythia8InvXSection->SetLineWidth(widthCommonFit);
             histoPythia8InvXSection->Draw("same,hist,l");

             DrawGammaSetMarkerTGraphErr(graphPythia8InvXSectionEtaScaled, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
             graphPythia8InvXSectionEtaScaled->Draw("3,same");
             DrawGammaSetMarker(histoPythia8InvXSectionEta, 24, 1.5, kRed+2 , kRed+2);
             histoPythia8InvXSectionEta->SetLineWidth(widthCommonFit);
             histoPythia8InvXSectionEta->Draw("same,hist,l");

             // plots fits
             fitTCMInvXSectionPi0->Draw("same");
             fitInvXSectionPi0->Draw("same");

             SetStyleHisto(histfitInvXSectionEta, 2, 7, kGray+2);
             histfitInvXSectionEta->Draw("same,c");
             DrawGammaSetMarkerTF1( tf1FitInvXSectionEta, 3, 2, kGray+1);
             tf1FitInvXSectionEta->Draw("same");

             // plot data
             graphCombPi0InvXSectionSys->Draw("E2same");
             DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSys_Copy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
             graphCombEtaInvXSectionSys_Copy->Draw("E2same");

             graphCombPi0InvXSectionStat_WOXErr->Draw("p,same,z");
             DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStat_WOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
             graphCombEtaInvXSectionStat_WOXErrCopy->Draw("p,same,z");

             // labels lower left corner
             TLegend* legendXsectionPaperAll    = GetAndSetLegend2(0.17, 0.12, 0.5, 0.11+0.04*4, textSizeLabelsPixel, 1, "", 43, 0.2);
             legendXsectionPaperAll->AddEntry(graphCombPi0InvXSectionSys,"#pi^{0}","pf");
             legendXsectionPaperAll->AddEntry(graphCombEtaInvXSectionSys_Copy,"#eta (x 10^{-2})","pf");
             legendXsectionPaperAll->AddEntry(fitTCMInvXSectionPi0,"TCM fit","l");
             legendXsectionPaperAll->AddEntry(fitInvXSectionPi0,"Tsallis fit","l");
             legendXsectionPaperAll->Draw();

             TLatex *labelEnergyXSectionPaperAll = new TLatex(0.18, 0.12+0.04*6, collisionSystem900GeV.Data());
             SetStyleTLatex( labelEnergyXSectionPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
             labelEnergyXSectionPaperAll->Draw();
             TLatex *labelALICEXSectionPaperAll  = new TLatex(0.18,0.12+0.04*5,"ALICE");
             SetStyleTLatex( labelALICEXSectionPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
             labelALICEXSectionPaperAll->Draw();
             TLatex *labelALICENormUnPaperAll    = new TLatex(0.18,0.12+0.04*4+0.003,"norm. unc. 3.5%");
             SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
             labelALICENormUnPaperAll->Draw();

             // labels upper right corner
             TLegend* legendXsectionPaperPyBoth  = GetAndSetLegend2(0.5, 0.95-0.04*5, 0.54+0.33, 0.95, textSizeLabelsPixel, 1, "", 43, 0.18);
             legendXsectionPaperPyBoth->AddEntry(histoPythia8InvXSectionEta,"PYTHIA 8.2, Monash 2013","l");
             legendXsectionPaperPyBoth->AddEntry(graphPi0DSS07,"#pi^{0} pQCD NLO ","f");
             legendXsectionPaperPyBoth->AddEntry((TObject*)0,"#scale[0.75]{PDF: MSTW, FF: DSS07}","");
             legendXsectionPaperPyBoth->AddEntry(graphEtaAESSSCopy,"#eta pQCD NLO ","f");
             legendXsectionPaperPyBoth->AddEntry((TObject*)0,"#scale[0.75]{PDF: CTEQ6M5, FF: AESSS}","");
             legendXsectionPaperPyBoth->Draw();

         canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Pi0_Eta_Theory.%s",outputDir.Data(),suffix.Data()));

         histo2DXSectionWithEtaAndPi0->Draw("copy");

             // plots fits
             fitTCMInvXSectionPi0->Draw("same");
             fitInvXSectionPi0->Draw("same");

             SetStyleHisto(histfitInvXSectionEta, 2, 7, kGray+2);
             histfitInvXSectionEta->Draw("same,c");
             DrawGammaSetMarkerTF1( tf1FitInvXSectionEta, 3, 2, kGray+1);
             tf1FitInvXSectionEta->Draw("same");

             // plot data
             graphCombPi0InvXSectionSys->Draw("E2same");
             DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSys_Copy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
             graphCombEtaInvXSectionSys_Copy->Draw("E2same");

             graphCombPi0InvXSectionStat_WOXErr->Draw("p,same,z");
             DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStat_WOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
             graphCombEtaInvXSectionStat_WOXErrCopy->Draw("p,same,z");

             // labels lower left corner
             legendXsectionPaperAll->Draw();

             labelEnergyXSectionPaperAll->Draw();
             labelALICEXSectionPaperAll->Draw();
             labelALICENormUnPaperAll->Draw();

         canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Pi0_Eta_Fits.%s",outputDir.Data(),suffix.Data()));

         histo2DXSectionWithEtaAndPi0->Draw("copy");

             // plot data
             graphCombPi0InvXSectionSys->Draw("E2same");
             DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSys_Copy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
             graphCombEtaInvXSectionSys_Copy->Draw("E2same");

             graphCombPi0InvXSectionStat_WOXErr->Draw("p,same,z");
             DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStat_WOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
             graphCombEtaInvXSectionStat_WOXErrCopy->Draw("p,same,z");

             // labels lower left corner
             TLegend* legendXsectionPaperAll2    = GetAndSetLegend2(0.17, 0.20, 0.5, 0.19+0.04*2, textSizeLabelsPixel, 1, "", 43, 0.2);
             legendXsectionPaperAll2->AddEntry(graphCombPi0InvXSectionSys,"#pi^{0}","pf");
             legendXsectionPaperAll2->AddEntry(graphCombEtaInvXSectionSys_Copy,"#eta (x 10^{-2})","pf");
             legendXsectionPaperAll2->Draw();

             labelEnergyXSectionPaperAll->Draw();
             labelALICEXSectionPaperAll->Draw();
             labelALICENormUnPaperAll->Draw();

         canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Pi0_Eta.%s",outputDir.Data(),suffix.Data()));

         histoPythia8InvXSectionEta->Scale(1/scaleFacEtaForCombPlot);

      // ***************************************************************************************************************
      // ******************************* eta/pi0 graphs without x-error  ***********************************************
      // ***************************************************************************************************************

      TGraphAsymmErrors* graphCombEtaToPi0Stat_WOXErr = (TGraphAsymmErrors*) graphCombEtaToPi0Stat->Clone("graphCombEtaToPi0Stat_WOXErr");
      ProduceGraphAsymmWithoutXErrors(graphCombEtaToPi0Stat_WOXErr);

      TGraphAsymmErrors* graphEtaToPi0Stat_WOXErr[11];
      for (Int_t i = 0; i < 11; i++){
        if(directoryEta[i]){
          graphEtaToPi0Stat_WOXErr[i] = (TGraphAsymmErrors*) graphEtaToPi0Stat[i]->Clone(Form("graphEtaToPi0Stat_%i_WOXErr",i));
          ProduceGraphAsymmWithoutXErrors(graphEtaToPi0Stat_WOXErr[i]);
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
      histo2DEtatoPi0combo               = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",1000,minPtEtaToPi0,35.,1000,-0.05,1.05    );
      SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                                0.85*textsizeLabelsEtaToPi0,1.1*textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
      histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
      histo2DEtatoPi0combo->GetXaxis()->SetNoExponent();
      histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(minPtEtaToPi0,maxPtEta);
      histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(-0.05,1.05);
      histo2DEtatoPi0combo->Draw("copy");

          // plotting systematics graphs
          for (Int_t i = 0; i < 11; i++){
              if(graphEtaToPi0Sys[i]){
                  DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Sys[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
                  graphEtaToPi0Sys[i]->Draw("E2same");
              }
          }
           DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Stat, markerStyleDet[1], markerSizeDet[1]*0.75, kPink , kPink, widthLinesBoxes, kTRUE);
          // plotting statistics graphs
          for (Int_t i = 0; i < 11; i++){
              if(graphEtaToPi0Stat[i]){
                  DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Stat[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i]);
                  graphEtaToPi0Stat[i]->Draw("p,same,e");
              }
          }
          DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys, markerStyleDet[1], markerSizeDet[1]*0.75, kPink , kPink);

          TLegend* legendEtaToPi0 = GetAndSetLegend2(0.67, 0.9, 0.9, 0.9+(textsizeLabelsEtaToPi0*1*0.9), textSizeLabelsPixel);
          for (Int_t i = 0; i < 11; i++){
              if(graphEtaToPi0Sys[i]){
                  legendEtaToPi0->AddEntry(graphEtaToPi0Sys[i],nameMeasGlobal[i],"pf");
              }
          }
          legendEtaToPi0->Draw();

          drawLatexAdd(collisionSystem900GeV.Data(),0.13, 0.92,0.85*textsizeLabelsEtaToPi0,kFALSE);
          drawLatexAdd("ALICE",0.13, 0.92-(1*textsizeLabelsEtaToPi0*0.85),0.85*textsizeLabelsEtaToPi0,kFALSE);

      canvasEtatoPi0combo->Update();
      canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_differentSystems.%s",outputDir.Data(), suffix.Data()));

      //*************************************************************************************************************
      //*************************************************************************************************************

      histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(minPtEtaToPi0,32.);
      histo2DEtatoPi0combo->Draw("copy");

      TLegend* legendXsectionPaperEtaToPi0     = GetAndSetLegend2(0.12, 0.8, 0.45, 0.96, 0.85*textSizeLabelsPixel);
      legendXsectionPaperEtaToPi0->SetNColumns(1);
      legendXsectionPaperEtaToPi0->SetMargin(0.2);
      legendXsectionPaperEtaToPi0->AddEntry(graphCombPi0InvXSectionSys,"Data","pf");
      legendXsectionPaperEtaToPi0->AddEntry(graphEtaToPi02760GeV,"ALICE pp, #sqrt{#it{s}} = 2.76 TeV","p");
      legendXsectionPaperEtaToPi0->AddEntry(graphEtaToPi08000GeV,"ALICE pp, #sqrt{#it{s}} = 8 TeV","p");
      legendXsectionPaperEtaToPi0->Draw();

      DrawGammaSetMarkerTGraphAsym(graphEtaToPi08000GeV, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[1] , colorDet[1], widthLinesBoxes, kTRUE);
      graphEtaToPi08000GeV->Draw("same,p");
      DrawGammaSetMarkerTGraphAsym(graphEtaToPi02760GeV, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
      graphEtaToPi02760GeV->Draw("same,p");

      // plotting labels
      TLatex *labelEnergyEtaToPi02 = new TLatex(0.71, 0.92,collisionSystem900GeV.Data());
      SetStyleTLatex( labelEnergyEtaToPi02, 0.85*textsizeLabelsEtaToPi0,4);
      labelEnergyEtaToPi02->Draw();

      TLatex *labelALICEEtaToPi02 = new TLatex(0.852, 0.92-(1*textsizeLabelsEtaToPi0*0.85),"ALICE");
      SetStyleTLatex( labelALICEEtaToPi02, 0.85*textsizeLabelsEtaToPi0,4);
      labelALICEEtaToPi02->Draw();

      // plotting data
      graphCombEtaToPi0Stat_WOXErr->Print();
      DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Stat_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
      graphCombEtaToPi0Stat_WOXErr->SetLineWidth(widthLinesBoxes);
      DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
      graphCombEtaToPi0Sys->SetLineWidth(0);
      graphCombEtaToPi0Sys->Draw("2,same");
      graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

      histo2DEtatoPi0combo->Draw("axis,same");

      canvasEtatoPi0combo->Update();
      canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Paper.%s",outputDir.Data(), suffix.Data()));

      //*************************************************************************************************************
      //*************************************************************************************************************

      histo2DEtatoPi0combo->Draw("copy");

      textSizeLabelsPixel = 48;
      TLegend* legendXsectionPaperEtaToPi03n     = GetAndSetLegend2(0.11, 0.81, 0.5, 0.96, 0.85*textSizeLabelsPixel);
      legendXsectionPaperEtaToPi03n->SetNColumns(1);
      legendXsectionPaperEtaToPi03n->SetMargin(0.15);
      legendXsectionPaperEtaToPi03n->AddEntry(graphCombPi0InvXSectionSys,"ALICE pp, #sqrt{#it{s}} = 900 GeV","pf");
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

      histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(minPtEtaToPi0,maxPtEta);
      histo2DEtatoPi0combo->Draw("copy");

      DrawGammaSetMarkerTGraphErr(graphPythia8EtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      graphPythia8EtaToPi0->Draw("3,same");
      DrawGammaSetMarker(histoPythia8EtaToPi0, 24, 1.5, kRed+2 , kRed+2);
      histoPythia8EtaToPi0->SetLineWidth(widthCommonFit);
      histoPythia8EtaToPi0->Draw("same,hist,l");

      textSizeLabelsPixel = 48;
      TLegend* legendXsectionPaperEtaToPi02     = GetAndSetLegend2(0.12, 0.69, 0.45, 0.69+0.045*6, 0.85*textSizeLabelsPixel);
      legendXsectionPaperEtaToPi02->SetNColumns(1);
      legendXsectionPaperEtaToPi02->SetMargin(0.2);
      legendXsectionPaperEtaToPi02->AddEntry(graphCombPi0InvXSectionSys,"Data","pf");
      legendXsectionPaperEtaToPi02->AddEntry(histoPythia8EtaToPi0,"PYTHIA 8.2, Monash 2013","l");
      legendXsectionPaperEtaToPi02->AddEntry(graphNLOEtaToPi0,"NLO, PDF:CTEQ6M5","f");
      legendXsectionPaperEtaToPi02->AddEntry((TObject*)0,"#pi^{0} FF:DSS07, #eta FF:AESSS","");
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
      graphCombEtaToPi0Sys->SetLineWidth(0);
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

      // plotting data
      graphCombEtaToPi0Sys->Draw("2,same");
      graphCombEtaToPi0Stat_WOXErr->Draw("p,same");

      // plotting labels
      labelEnergyEtaToPi02->Draw();
      labelALICEEtaToPi02->Draw();

      histo2DEtatoPi0combo->Draw("axis,same");

      canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Combined.%s",outputDir.Data(), suffix.Data()));

      //*************************************************************************************************************
      //*************************************************************************************************************

      canvasEtatoPi0combo->SetRightMargin(0.02);
      histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(0.);
      histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(0.,5.0);
      histo2DEtatoPi0combo->Draw("copy");
      legendXsectionPaperEtaToPi02->Draw();

      //plotting MC
      DrawGammaSetMarkerTGraphErr(graphPythia8EtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
      graphPythia8EtaToPi0->Draw("3,same");
      DrawGammaSetMarker(histoPythia8EtaToPi0, 24, 1.5, kRed+2 , kRed+2);
      histoPythia8EtaToPi0->SetLineWidth(widthCommonFit);
      histoPythia8EtaToPi0->Draw("same,hist,l");

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

      canvasEtatoPi0combo->SetLogx(kFALSE);
      canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Combined_LIN.%s",outputDir.Data(), suffix.Data()));

      // **********************************************************************************************************************
      // ************************* Saving of final results ********************************************************************
      // **********************************************************************************************************************

         TString nameOutputCommonFile    = Form("CombinedResultsPaperPP900GeV_%s.root", dateForOutput.Data());
         TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

         fCombResults.mkdir("Pi0900GeV");
         TDirectoryFile* fileDirectoryPi0 = (TDirectoryFile*)fCombResults.Get("Pi0900GeV");
         fCombResults.cd("Pi0900GeV");
             graphCombPi0InvXSectionTot->Write("graphInvCrossSectionPi0Comb900GeVA");
             graphCombPi0InvXSectionStat->Write("graphInvCrossSectionPi0Comb900GeVAStatErr");
             graphCombPi0InvXSectionSys->Write("graphInvCrossSectionPi0Comb900GeVASysErr");

             for (Int_t i = 0; i < 11; i++){
               if(directoryPi0[i]){
                 graphPi0InvXSectionStat[i]          ->Write(Form("graphInvCrossSectionPi0%s900GeVStatErr",nameMeasGlobalWriteToFile[i].Data()));
                 graphPi0InvXSectionSys[i]           ->Write(Form("graphInvCrossSectionPi0%s900GeVSysErr",nameMeasGlobalWriteToFile[i].Data()));
               }
             }

              // fits for pi0
             fitInvXSectionPi0->Write("TsallisFitPi0");
             fitTCMInvXSectionPi0->Write("TwoComponentModelFitPi0");

             if (bWCorrection.Contains("Y")){
                 if(graphCombPi0InvXSectionTot_yShifted)graphCombPi0InvXSectionTot_yShifted->Write("graphInvCrossSectionPi0Comb900GeV_yShifted");
                 if(graphCombPi0InvXSectionStat_yShifted)graphCombPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0Comb900GeVStatErr_yShifted");
                 if(graphCombPi0InvXSectionSys_yShifted)graphCombPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0Comb900GeVSysErr_yShifted");

                 for (Int_t i = 0; i < 11; i++){
                   if(directoryPi0[i]){
                     graphPi0InvXSectionStat_yShifted[i]          ->Write(Form("graphInvCrossSectionPi0%s900GeVStatErr_yShifted",nameMeasGlobalWriteToFile[i].Data()));
                     graphPi0InvXSectionSys_yShifted[i]           ->Write(Form("graphInvCrossSectionPi0%s900GeVSysErr_yShifted",nameMeasGlobalWriteToFile[i].Data()));
                   }
                 }

             }

             fileDirectoryPi0->mkdir("Supporting");
             fileDirectoryPi0->cd("Supporting");
             for (Int_t i = 0; i < 11; i++){
               if(directoryPi0[i] && i!=1){
                 // Writing full correction factors
                 histoPi0AccTimesEff[i]          ->Write(Form("Pi0CorrectionFactor%s",nameMeasGlobalWriteToFile[i].Data()));
                 histoPi0Mass[i]                 ->Write(Form("Pi0MassData%s",nameMeasGlobalWriteToFile[i].Data()));
                 histoPi0TrueMass[i]             ->Write(Form("Pi0MassMC%s",nameMeasGlobalWriteToFile[i].Data()));
                 histoPi0FWHMMeV[i]              ->Write(Form("Pi0WidthData%s",nameMeasGlobalWriteToFile[i].Data()));
                 histoPi0TrueFWHMMeV[i]          ->Write(Form("Pi0WidthMC%s",nameMeasGlobalWriteToFile[i].Data()));
               }
             }

         fCombResults.mkdir("Eta900GeV");
         TDirectoryFile* fileDirectoryEta = (TDirectoryFile*)fCombResults.Get("Eta900GeV");
         fCombResults.cd("Eta900GeV");
             graphCombEtaInvXSectionTot->Write("graphInvCrossSectionEtaComb900GeVA");
             graphCombEtaInvXSectionStat->Write("graphInvCrossSectionEtaComb900GeVAStatErr");
             graphCombEtaInvXSectionSys->Write("graphInvCrossSectionEtaComb900GeVASysErr");

             for (Int_t i = 0; i < 11; i++){
               if(directoryEta[i]){
                 graphEtaInvXSectionStat[i]          ->Write(Form("graphInvCrossSectionEta%s900GeVStatErr",nameMeasGlobalWriteToFile[i].Data()));
                 graphEtaInvXSectionSys[i]           ->Write(Form("graphInvCrossSectionEta%s900GeVSysErr",nameMeasGlobalWriteToFile[i].Data()));
               }
             }

              // fits for Eta
             fitInvXSectionEta->Write("TsallisFitEta");
             fitInvXSectionEta->Write("TwoComponentModelFitEta");

             if (bWCorrection.Contains("Y")){
                 if(graphCombEtaInvXSectionTot_yShifted)graphCombEtaInvXSectionTot_yShifted->Write("graphInvCrossSectionEtaComb900GeV_yShifted");
                 if(graphCombEtaInvXSectionStat_yShifted)graphCombEtaInvXSectionStat_yShifted->Write("graphInvCrossSectionEtaComb900GeVStatErr_yShifted");
                 if(graphCombEtaInvXSectionSys_yShifted)graphCombEtaInvXSectionSys_yShifted->Write("graphInvCrossSectionEtaComb900GeVSysErr_yShifted");

                 for (Int_t i = 0; i < 11; i++){
                   if(directoryEta[i]){
                     graphEtaInvXSectionStat_yShifted[i]          ->Write(Form("graphInvCrossSectionEta%s900GeVStatErr_yShifted",nameMeasGlobalWriteToFile[i].Data()));
                     graphEtaInvXSectionSys_yShifted[i]           ->Write(Form("graphInvCrossSectionEta%s900GeVSysErr_yShifted",nameMeasGlobalWriteToFile[i].Data()));
                   }
                 }

             }


             graphCombEtaToPi0Tot->Write("graphRatioEtaToPi0Comb900GeVTotErr");
             graphCombEtaToPi0Stat->Write("graphRatioEtaToPi0Comb900GeVStatErr");
             graphCombEtaToPi0Sys->Write("graphRatioEtaToPi0Comb900GeVSysErr");

             for (Int_t i = 0; i < 11; i++){
               if(directoryEta[i]){
                 graphEtaToPi0Stat[i]          ->Write(Form("graphRatioEtaToPi0%s900GeVStatErr",nameMeasGlobalWriteToFile[i].Data()));
                 graphEtaToPi0Sys[i]           ->Write(Form("graphRatioEtaToPi0%s900GeVSysErr",nameMeasGlobalWriteToFile[i].Data()));
               }
             }

             fileDirectoryEta->mkdir("Supporting");
             fileDirectoryEta->cd("Supporting");
             for (Int_t i = 0; i < 11; i++){
               if(directoryEta[i]){
                 // Writing full correction factors
                 histoEtaAccTimesEff[i]          ->Write(Form("EtaCorrectionFactor%s",nameMeasGlobalWriteToFile[i].Data()));
                 histoEtaMass[i]                 ->Write(Form("EtaMassData%s",nameMeasGlobalWriteToFile[i].Data()));
                 histoEtaTrueMass[i]             ->Write(Form("EtaMassMC%s",nameMeasGlobalWriteToFile[i].Data()));
                 histoEtaFWHMMeV[i]              ->Write(Form("EtaWidthData%s",nameMeasGlobalWriteToFile[i].Data()));
                 histoEtaTrueFWHMMeV[i]          ->Write(Form("EtaWidthMC%s",nameMeasGlobalWriteToFile[i].Data()));
               }
             }
         fCombResults.Close();

      // **********************************************************************************************************************
      // ************************* Saving only fits to final results **********************************************************
      // **********************************************************************************************************************

         TString nameOutputCommonFileFitsOnly    = Form("FitsPaperPP900GeV_%s.root", dateForOutput.Data());
         TFile fFitsResults(nameOutputCommonFileFitsOnly.Data(), "RECREATE");

              // fits for pi0
             fitInvXSectionPi0->Write("TsallisFitPi0");
             fitTCMInvXSectionPi0->Write("TwoComponentModelFitPi0");

             // fits for eta
             fitInvXSectionEta->Write("TsallisFitEta");

         fFitsResults.Close();

      // **********************************************************************************************************************
      // **************************Plot example invariant mass bins ***********************************************************
      // **********************************************************************************************************************
}
