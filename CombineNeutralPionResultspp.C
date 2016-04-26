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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"

extern TRandom*   gRandom;
extern TBenchmark*   gBenchmark;
extern TSystem*   gSystem;
extern TMinuit*   gMinuit;

void ScaleMCYield(TH1D* histoCorrectedToBeScaled, Double_t deltaRapid, Double_t scaling, Double_t nEvtMC, TString nameMeson, TString optionDalitz ){
   histoCorrectedToBeScaled->Sumw2();
   histoCorrectedToBeScaled->Scale(1./deltaRapid);
   histoCorrectedToBeScaled->Scale(scaling);
   histoCorrectedToBeScaled->Scale(1./nEvtMC);
   for (Int_t i = 1; i < histoCorrectedToBeScaled->GetNbinsX()+1 ; i++){
      Double_t newBinContent = histoCorrectedToBeScaled->GetBinContent(i)/histoCorrectedToBeScaled->GetBinCenter(i);
      Double_t newBinError = histoCorrectedToBeScaled->GetBinError(i)/histoCorrectedToBeScaled->GetBinCenter(i);
      histoCorrectedToBeScaled->SetBinContent(i,newBinContent);
      histoCorrectedToBeScaled->SetBinError(i,newBinError);
   }
   if (nameMeson.CompareTo("Pi0") == 0 ||nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
      if (optionDalitz.CompareTo("kFALSE")==0){
         histoCorrectedToBeScaled->Scale(1./0.98798);
      } else {
         histoCorrectedToBeScaled->Scale(1./0.01198);
      }
   }else{
      if (optionDalitz.CompareTo("kFALSE")==0){
         histoCorrectedToBeScaled->Scale(1./0.3931);
      } else {
         histoCorrectedToBeScaled->Scale(1./6.8e-5);
      }
      
   }
}


void CombineNeutralPionResultspp(TString suffix = "pdf", 
				 TString nameFilepp1 = "Pythia", 
				 TString nameFilepp2 = "PythiaAddSign",
				 TString nameFilepp3 = "Phojet", 
				 Bool_t runDrawReweighted = kTRUE){


   gROOT->Reset();   
   gROOT->SetStyle("Plain");
      
   TString dateForOutput = ReturnDateStringForOutput();
   TString outputDir = Form("%s/%s/CombineNeutralPionResultspp",suffix.Data(),dateForOutput.Data());
   gSystem->Exec("mkdir -p "+outputDir);
   gSystem->Exec(Form("cp -r %s %s/InputFilePCMPPPythia.root ",nameFilepp1.Data(),outputDir.Data())); 
   gSystem->Exec(Form("cp -r %s %s/InputFilePCMPPPythiaPlus.root ",nameFilepp2.Data(),outputDir.Data())); 
   gSystem->Exec(Form("cp -r %s %s/InputFilePCMPPPhojet.root ",nameFilepp3.Data(),outputDir.Data())); 


   gSystem->Exec("mkdir -p "+outputDir);
   
   StyleSettingsThesis();  
   SetPlotStyle();

   Color_t  colorCombPP             = kBlack;
   Color_t  colorCombMCPP           = kRed+3;
   Style_t  markerStylePP     = 33 ;   
   Style_t  markerStylePPMC   = 24 ;   
   Size_t   markerSizePP      = 2.5;
   
   Color_t  colorPi0900GeV          = kRed +2;
   Color_t  colorPi02760GeV         = kMagenta+2;
   Color_t  colorPi07TeV            = kBlue+2;
   Color_t  colorPi08TeV            = kGreen+2;
   Color_t  colorPi0900GeVBox = colorPi0900GeV-10;
   Color_t  colorPi02760GeVBox = colorPi02760GeV-10;
   Color_t  colorPi07TeVBox = colorPi07TeV-10;

   Color_t  colorMCPythiaPP900GeV   = colorPi0900GeV-4;
   Color_t  colorMCPythiaPP2760GeV = colorPi02760GeV+2;
   Color_t  colorMCPythiaPP7TeV  = colorPi07TeV+3;
   Color_t  colorMCPythiaPP8TeV  = colorPi08TeV+2;
   Color_t  colorMCPhojetPP900GeV   = colorPi0900GeV+2;
   Color_t  colorMCPhojetPP2760GeV = colorPi02760GeV-4;
   Color_t  colorMCPhojetPP7TeV  = colorPi07TeV-3;
   Color_t  colorMCPhojetPP8TeV  = colorPi08TeV-2;

   Style_t  markerStyleSpectrum7TeVMC  = 24 ;
   Style_t  markerStyleSpectrum900GeVMC = 25 ;
   Style_t  markerStyleSpectrum8TeVMC = 22 ;
   Style_t  markerStyleSpectrum2760GeVMC = 30 ;
   Style_t  markerStyleSpectrum7TeV    = 20 ;
   Style_t  markerStyleSpectrum900GeV = 21 ;
   Style_t  markerStyleSpectrum8TeV = 24 ;
   Style_t  markerStyleSpectrum2760GeV = 29 ;
   
   Double_t xSection7TeVppINEL = 73.2*1e9;
   Double_t xSection2760GeVppINEL = 62.8*1e9;
   Double_t xSection900GeVppINEL = 52.5*1e9;

   Style_t  markerStyleMCPP8TeV  = 24 ;   
   Style_t  markerStyleMCPP7TeV  = 24 ;
   Style_t  markerStyleMCPP900GeV   = 25 ;
   Style_t  markerStyleMCPP2760GeV  = 30 ;

   Size_t   markerSizePi0PP8TeV  = 1.8;   
   Size_t   markerSizePi0PP7TeV  = 1.8;
   Size_t   markerSizePi0PP900GeV = 1.8;
   Size_t   markerSizePi0PP2760GeV  = 2.2;

   TString collisionSystemPP2760GeV = "pp #sqrt{#it{s}} = 2.76 TeV";      
   TString collisionSystemPP7TeV = "pp #sqrt{#it{s}} = 7 TeV";      
   TString collisionSystemPP8TeV = "pp #sqrt{#it{s}} = 8 TeV";         
   TString collisionSystemPP900GeV = "pp #sqrt{#it{s}} = 0.9 TeV";     
  
 
   Size_t markerSizeComparison = 0.5;
   Double_t maxPtMesonEffFit = 12.;
   Double_t minPtMesonEffFit = 1.2;
   Int_t offsetCorrectionHighPt= 1;
   TF1* fitTrueEffi = new TF1("EffiFitDummy","1 - [0]*exp([1]*x)+[1]");
   fitTrueEffi->SetRange(minPtMesonEffFit,maxPtMesonEffFit);


   //******************** Retrieving plots from PP file ************************
   //**************************   PYTHIA only   *************************************
   
   TFile*   filePythia =              new TFile(nameFilepp1.Data());

   //*********   8 TeV   *****************************************************************************************************
   TDirectory* directoryPythiaPi08TeV =  (TDirectory*)filePythia->Get("Pi08TeV"); 
   
   TH1D* histoPi0CorrectedSpecPythia8TeV =   (TH1D*)directoryPythiaPi08TeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPythia8TeV =        (TGraphAsymmErrors*)directoryPythiaPi08TeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPythia8TeV =        (TH1D*)directoryPythiaPi08TeV->Get("MassPi0");
   TH1D* histoPi0MassMCPythia8TeV =          (TH1D*)directoryPythiaPi08TeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPythia8TeV =       (TH1D*)directoryPythiaPi08TeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPythia8TeV =         (TH1D*)directoryPythiaPi08TeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPythia8TeV =        (TH1D*)directoryPythiaPi08TeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPythia8TeV =          (TH1D*)directoryPythiaPi08TeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPythia8TeV =         (TH1D*)directoryPythiaPi08TeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPythia8TeVWOWeights = (TH1D*)directoryPythiaPi08TeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPythia8TeV =        (TH1D*)directoryPythiaPi08TeV->Get("Pi0_Weights");
   TH1D* histoPi0RawYieldPythia8TeV =          (TH1D*)directoryPythiaPi08TeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPythia8TeV->Scale(1000.);
   histoPi0MassMCPythia8TeV->Scale(1000.);

   
   TDirectory* directoryPythiaEta8TeV =  (TDirectory*)filePythia->Get("Eta_pp_8TeV"); 
   
   TH1D* histoEtaCorrectedSpecPythia8TeV =   (TH1D*)directoryPythiaEta8TeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPythia8TeV =        (TGraphAsymmErrors*)directoryPythiaEta8TeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPythia8TeV =        (TH1D*)directoryPythiaEta8TeV->Get("MassEta");
   TH1D* histoEtaMassMCPythia8TeV =          (TH1D*)directoryPythiaEta8TeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPythia8TeV =       (TH1D*)directoryPythiaEta8TeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPythia8TeV =         (TH1D*)directoryPythiaEta8TeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPythia8TeV =        (TH1D*)directoryPythiaEta8TeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPythia8TeV =          (TH1D*)directoryPythiaEta8TeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPythia8TeV =         (TH1D*)directoryPythiaEta8TeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPythia8TeVWOWeights = (TH1D*)directoryPythiaEta8TeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPythia8TeV =        (TH1D*)directoryPythiaEta8TeV->Get("Eta_Weights");
   TH1D* histoEtaRawYieldPythia8TeV =          (TH1D*)directoryPythiaEta8TeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPythia8TeV->Scale(1000.);
   histoEtaMassMCPythia8TeV->Scale(1000.);
   
   
   //*********   7 TeV   *****************************************************************************************************
   TDirectory* directoryPythiaPi07TeV =  (TDirectory*)filePythia->Get("Pi0_pp_7TeV"); 
   
   TH1D* histoPi0CorrectedSpecPythia7TeV =   (TH1D*)directoryPythiaPi07TeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPythia7TeV =        (TGraphAsymmErrors*)directoryPythiaPi07TeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPythia7TeV =        (TH1D*)directoryPythiaPi07TeV->Get("MassPi0");
   TH1D* histoPi0MassMCPythia7TeV =          (TH1D*)directoryPythiaPi07TeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPythia7TeV =       (TH1D*)directoryPythiaPi07TeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPythia7TeV =         (TH1D*)directoryPythiaPi07TeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPythia7TeV =        (TH1D*)directoryPythiaPi07TeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPythia7TeV =          (TH1D*)directoryPythiaPi07TeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPythia7TeV =         (TH1D*)directoryPythiaPi07TeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPythia7TeVWOWeights = (TH1D*)directoryPythiaPi07TeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPythia7TeV =        (TH1D*)directoryPythiaPi07TeV->Get("Pi0_Weights");
   TH1D* histoPi0RawYieldPythia7TeV =          (TH1D*)directoryPythiaPi07TeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPythia7TeV->Scale(1000.);
   histoPi0MassMCPythia7TeV->Scale(1000.);

   
   TDirectory* directoryPythiaEta7TeV =  (TDirectory*)filePythia->Get("Eta_pp_7TeV"); 
   
   TH1D* histoEtaCorrectedSpecPythia7TeV =   (TH1D*)directoryPythiaEta7TeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPythia7TeV =        (TGraphAsymmErrors*)directoryPythiaEta7TeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPythia7TeV =        (TH1D*)directoryPythiaEta7TeV->Get("MassEta");
   TH1D* histoEtaMassMCPythia7TeV =          (TH1D*)directoryPythiaEta7TeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPythia7TeV =       (TH1D*)directoryPythiaEta7TeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPythia7TeV =         (TH1D*)directoryPythiaEta7TeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPythia7TeV =        (TH1D*)directoryPythiaEta7TeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPythia7TeV =          (TH1D*)directoryPythiaEta7TeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPythia7TeV =         (TH1D*)directoryPythiaEta7TeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPythia7TeVWOWeights = (TH1D*)directoryPythiaEta7TeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPythia7TeV =        (TH1D*)directoryPythiaEta7TeV->Get("Eta_Weights");
   TH1D* histoEtaRawYieldPythia7TeV =          (TH1D*)directoryPythiaEta7TeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPythia7TeV->Scale(1000.);
   histoEtaMassMCPythia7TeV->Scale(1000.);
   
   
   //*********   2.76 TeV   *****************************************************************************************************
   TDirectory* directoryPythiaPi02760GeV =  (TDirectory*)filePythia->Get("Pi0_pp_2760GeV"); 
   
   TH1D* histoPi0CorrectedSpecPythia2760GeV =   (TH1D*)directoryPythiaPi02760GeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPythia2760GeV =        (TGraphAsymmErrors*)directoryPythiaPi02760GeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPythia2760GeV =        (TH1D*)directoryPythiaPi02760GeV->Get("MassPi0");
   TH1D* histoPi0MassMCPythia2760GeV =          (TH1D*)directoryPythiaPi02760GeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPythia2760GeV =       (TH1D*)directoryPythiaPi02760GeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPythia2760GeV =         (TH1D*)directoryPythiaPi02760GeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPythia2760GeV =        (TH1D*)directoryPythiaPi02760GeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPythia2760GeV =          (TH1D*)directoryPythiaPi02760GeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPythia2760GeV =         (TH1D*)directoryPythiaPi02760GeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPythia2760GeVWOWeights = (TH1D*)directoryPythiaPi02760GeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPythia2760GeV =        (TH1D*)directoryPythiaPi02760GeV->Get("Pi0_Weights");
   TH1D* histoPi0RawYieldPythia2760GeV =          (TH1D*)directoryPythiaPi02760GeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPythia2760GeV->Scale(1000.);
   histoPi0MassMCPythia2760GeV->Scale(1000.);

   
   TDirectory* directoryPythiaEta2760GeV =  (TDirectory*)filePythia->Get("Eta_pp_2760GeV"); 
   
   TH1D* histoEtaCorrectedSpecPythia2760GeV =   (TH1D*)directoryPythiaEta2760GeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPythia2760GeV =        (TGraphAsymmErrors*)directoryPythiaEta2760GeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPythia2760GeV =        (TH1D*)directoryPythiaEta2760GeV->Get("MassEta");
   TH1D* histoEtaMassMCPythia2760GeV =          (TH1D*)directoryPythiaEta2760GeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPythia2760GeV =       (TH1D*)directoryPythiaEta2760GeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPythia2760GeV =         (TH1D*)directoryPythiaEta2760GeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPythia2760GeV =        (TH1D*)directoryPythiaEta2760GeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPythia2760GeV =          (TH1D*)directoryPythiaEta2760GeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPythia2760GeV =         (TH1D*)directoryPythiaEta2760GeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPythia2760GeVWOWeights = (TH1D*)directoryPythiaEta2760GeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPythia2760GeV =        (TH1D*)directoryPythiaEta2760GeV->Get("Eta_Weights");
   TH1D* histoEtaRawYieldPythia2760GeV =          (TH1D*)directoryPythiaEta2760GeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPythia2760GeV->Scale(1000.);
   histoEtaMassMCPythia2760GeV->Scale(1000.);
   
   
   //*********  0.9 TeV   *****************************************************************************************************
   TDirectory* directoryPythiaPi0900GeV =  (TDirectory*)filePythia->Get("Pi0_pp_900GeV"); 
   
   TH1D* histoPi0CorrectedSpecPythia900GeV =   (TH1D*)directoryPythiaPi0900GeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPythia900GeV =        (TGraphAsymmErrors*)directoryPythiaPi0900GeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPythia900GeV =        (TH1D*)directoryPythiaPi0900GeV->Get("MassPi0");
   TH1D* histoPi0MassMCPythia900GeV =          (TH1D*)directoryPythiaPi0900GeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPythia900GeV =       (TH1D*)directoryPythiaPi0900GeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPythia900GeV =         (TH1D*)directoryPythiaPi0900GeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPythia900GeV =        (TH1D*)directoryPythiaPi0900GeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPythia900GeV =          (TH1D*)directoryPythiaPi0900GeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPythia900GeV =         (TH1D*)directoryPythiaPi0900GeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPythia900GeVWOWeights = (TH1D*)directoryPythiaPi0900GeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPythia900GeV =        (TH1D*)directoryPythiaPi0900GeV->Get("Pi0_Weights");
   TH1D* histoPi0RawYieldPythia900GeV =          (TH1D*)directoryPythiaPi0900GeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPythia900GeV->Scale(1000.);
   histoPi0MassMCPythia900GeV->Scale(1000.);

   
   TDirectory* directoryPythiaEta900GeV =  (TDirectory*)filePythia->Get("Eta_pp_900GeV"); 
   
   TH1D* histoEtaCorrectedSpecPythia900GeV =   (TH1D*)directoryPythiaEta900GeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPythia900GeV =        (TGraphAsymmErrors*)directoryPythiaEta900GeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPythia900GeV =        (TH1D*)directoryPythiaEta900GeV->Get("MassEta");
   TH1D* histoEtaMassMCPythia900GeV =          (TH1D*)directoryPythiaEta900GeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPythia900GeV =       (TH1D*)directoryPythiaEta900GeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPythia900GeV =         (TH1D*)directoryPythiaEta900GeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPythia900GeV =        (TH1D*)directoryPythiaEta900GeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPythia900GeV =          (TH1D*)directoryPythiaEta900GeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPythia900GeV =         (TH1D*)directoryPythiaEta900GeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPythia900GeVWOWeights = (TH1D*)directoryPythiaEta900GeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPythia900GeV =        (TH1D*)directoryPythiaEta900GeV->Get("Eta_Weights");
   TH1D* histoEtaRawYieldPythia900GeV =          (TH1D*)directoryPythiaEta900GeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPythia900GeV->Scale(1000.);
   histoEtaMassMCPythia900GeV->Scale(1000.);
   
   
   
   
  //*********************************************************************************************
  //**************************   PYTHIA plus Added Signal   *************************************
      TFile*   filePythiaPlus =              new TFile(nameFilepp2.Data());

   //*********   8 TeV   *****************************************************************************************************
   TDirectory* directoryPythiaPlusPi08TeV =  (TDirectory*)filePythiaPlus->Get("Pi0_pp_8TeV"); 
   
   TH1D* histoPi0CorrectedSpecPythiaPlus8TeV =   (TH1D*)directoryPythiaPlusPi08TeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPythiaPlus8TeV =        (TGraphAsymmErrors*)directoryPythiaPlusPi08TeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPythiaPlus8TeV =        (TH1D*)directoryPythiaPlusPi08TeV->Get("MassPi0");
   TH1D* histoPi0MassMCPythiaPlus8TeV =          (TH1D*)directoryPythiaPlusPi08TeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPythiaPlus8TeV =       (TH1D*)directoryPythiaPlusPi08TeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPythiaPlus8TeV =         (TH1D*)directoryPythiaPlusPi08TeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPythiaPlus8TeV =        (TH1D*)directoryPythiaPlusPi08TeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPythiaPlus8TeV =          (TH1D*)directoryPythiaPlusPi08TeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPythiaPlus8TeV =         (TH1D*)directoryPythiaPlusPi08TeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPythiaPlus8TeVWOWeights = (TH1D*)directoryPythiaPlusPi08TeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPythiaPlus8TeV =        (TH1D*)directoryPythiaPlusPi08TeV->Get("Pi0_Weights");
   TH1D* histoPi0RawYieldPythiaPlus8TeV =          (TH1D*)directoryPythiaPlusPi08TeV->Get("RAWYieldPerEventsPi0");
   TH1D* histoMCPi0YieldPtPythiaPlus8TeVAddedSig =         (TH1D*)directoryPythiaPlusPi08TeV->Get("Pi0_Input_Reweighted_AddedSig");
   TH1D* histoMCPi0YieldPtPythiaPlus8TeVAddedSigWOWeights = (TH1D*)directoryPythiaPlusPi08TeV->Get("Pi0_Input_AddedSig");


   histoPi0MassDataPythiaPlus8TeV->Scale(1000.);
   histoPi0MassMCPythiaPlus8TeV->Scale(1000.);

   
   TDirectory* directoryPythiaPlusEta8TeV =  (TDirectory*)filePythiaPlus->Get("Eta_pp_8TeV"); 
   
   TH1D* histoEtaCorrectedSpecPythiaPlus8TeV =   (TH1D*)directoryPythiaPlusEta8TeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPythiaPlus8TeV =        (TGraphAsymmErrors*)directoryPythiaPlusEta8TeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPythiaPlus8TeV =        (TH1D*)directoryPythiaPlusEta8TeV->Get("MassEta");
   TH1D* histoEtaMassMCPythiaPlus8TeV =          (TH1D*)directoryPythiaPlusEta8TeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPythiaPlus8TeV =       (TH1D*)directoryPythiaPlusEta8TeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPythiaPlus8TeV =         (TH1D*)directoryPythiaPlusEta8TeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPythiaPlus8TeV =        (TH1D*)directoryPythiaPlusEta8TeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPythiaPlus8TeV =          (TH1D*)directoryPythiaPlusEta8TeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPythiaPlus8TeV =         (TH1D*)directoryPythiaPlusEta8TeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPythiaPlus8TeVWOWeights = (TH1D*)directoryPythiaPlusEta8TeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPythiaPlus8TeV =        (TH1D*)directoryPythiaPlusEta8TeV->Get("Eta_Weights");
   TH1D* histoMCEtaYieldPtPythiaPlus8TeVAddedSig =         (TH1D*)directoryPythiaPlusEta8TeV->Get("Eta_Input_Reweighted_AddedSig");
   TH1D* histoMCEtaYieldPtPythiaPlus8TeVAddedSigWOWeights = (TH1D*)directoryPythiaPlusEta8TeV->Get("Eta_Input_AddedSig");

   //TH1D* histoEtaRawYieldPythiaPlus8TeV =          (TH1D*)directoryPythiaPlusEta8TeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPythiaPlus8TeV->Scale(1000.);
   histoEtaMassMCPythiaPlus8TeV->Scale(1000.);
   
   
   //*********   7 TeV   *****************************************************************************************************
   TDirectory* directoryPythiaPlusPi07TeV =  (TDirectory*)filePythiaPlus->Get("Pi0_pp_7TeV"); 
   
   TH1D* histoPi0CorrectedSpecPythiaPlus7TeV =   (TH1D*)directoryPythiaPlusPi07TeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPythiaPlus7TeV =        (TGraphAsymmErrors*)directoryPythiaPlusPi07TeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPythiaPlus7TeV =        (TH1D*)directoryPythiaPlusPi07TeV->Get("MassPi0");
   TH1D* histoPi0MassMCPythiaPlus7TeV =          (TH1D*)directoryPythiaPlusPi07TeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPythiaPlus7TeV =       (TH1D*)directoryPythiaPlusPi07TeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPythiaPlus7TeV =         (TH1D*)directoryPythiaPlusPi07TeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPythiaPlus7TeV =        (TH1D*)directoryPythiaPlusPi07TeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPythiaPlus7TeV =          (TH1D*)directoryPythiaPlusPi07TeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPythiaPlus7TeV =         (TH1D*)directoryPythiaPlusPi07TeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPythiaPlus7TeVWOWeights = (TH1D*)directoryPythiaPlusPi07TeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPythiaPlus7TeV =        (TH1D*)directoryPythiaPlusPi07TeV->Get("Pi0_Weights");
   TH1D* histoMCPi0YieldPtPythiaPlus7TeVAddedSig =         (TH1D*)directoryPythiaPlusPi07TeV->Get("Pi0_Input_Reweighted_AddedSig");
   TH1D* histoMCPi0YieldPtPythiaPlus7TeVAddedSigWOWeights = (TH1D*)directoryPythiaPlusPi07TeV->Get("Pi0_Input_AddedSig");

   TH1D* histoPi0RawYieldPythiaPlus7TeV =          (TH1D*)directoryPythiaPlusPi07TeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPythiaPlus7TeV->Scale(1000.);
   histoPi0MassMCPythiaPlus7TeV->Scale(1000.);

   
   TDirectory* directoryPythiaPlusEta7TeV =  (TDirectory*)filePythiaPlus->Get("Eta_pp_7TeV"); 
   
   TH1D* histoEtaCorrectedSpecPythiaPlus7TeV =   (TH1D*)directoryPythiaPlusEta7TeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPythiaPlus7TeV =        (TGraphAsymmErrors*)directoryPythiaPlusEta7TeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPythiaPlus7TeV =        (TH1D*)directoryPythiaPlusEta7TeV->Get("MassEta");
   TH1D* histoEtaMassMCPythiaPlus7TeV =          (TH1D*)directoryPythiaPlusEta7TeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPythiaPlus7TeV =       (TH1D*)directoryPythiaPlusEta7TeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPythiaPlus7TeV =         (TH1D*)directoryPythiaPlusEta7TeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPythiaPlus7TeV =        (TH1D*)directoryPythiaPlusEta7TeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPythiaPlus7TeV =          (TH1D*)directoryPythiaPlusEta7TeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPythiaPlus7TeV =         (TH1D*)directoryPythiaPlusEta7TeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPythiaPlus7TeVWOWeights = (TH1D*)directoryPythiaPlusEta7TeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPythiaPlus7TeV =        (TH1D*)directoryPythiaPlusEta7TeV->Get("Eta_Weights");
   TH1D* histoMCEtaYieldPtPythiaPlus7TeVAddedSig =         (TH1D*)directoryPythiaPlusEta7TeV->Get("Eta_Input_Reweighted_AddedSig");
   TH1D* histoMCEtaYieldPtPythiaPlus7TeVAddedSigWOWeights = (TH1D*)directoryPythiaPlusEta7TeV->Get("Eta_Input_AddedSig");

   TH1D* histoEtaRawYieldPythiaPlus7TeV =          (TH1D*)directoryPythiaPlusEta7TeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPythiaPlus7TeV->Scale(1000.);
   histoEtaMassMCPythiaPlus7TeV->Scale(1000.);
   
   
   //*********   2.76 TeV   *****************************************************************************************************
   TDirectory* directoryPythiaPlusPi02760GeV =  (TDirectory*)filePythiaPlus->Get("Pi0_pp_2760GeV"); 
   
   TH1D* histoPi0CorrectedSpecPythiaPlus2760GeV =   (TH1D*)directoryPythiaPlusPi02760GeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPythiaPlus2760GeV =        (TGraphAsymmErrors*)directoryPythiaPlusPi02760GeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPythiaPlus2760GeV =        (TH1D*)directoryPythiaPlusPi02760GeV->Get("MassPi0");
   TH1D* histoPi0MassMCPythiaPlus2760GeV =          (TH1D*)directoryPythiaPlusPi02760GeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPythiaPlus2760GeV =       (TH1D*)directoryPythiaPlusPi02760GeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPythiaPlus2760GeV =         (TH1D*)directoryPythiaPlusPi02760GeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPythiaPlus2760GeV =        (TH1D*)directoryPythiaPlusPi02760GeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPythiaPlus2760GeV =          (TH1D*)directoryPythiaPlusPi02760GeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPythiaPlus2760GeV =         (TH1D*)directoryPythiaPlusPi02760GeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPythiaPlus2760GeVWOWeights = (TH1D*)directoryPythiaPlusPi02760GeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPythiaPlus2760GeV =        (TH1D*)directoryPythiaPlusPi02760GeV->Get("Pi0_Weights");
   TH1D* histoMCPi0YieldPtPythiaPlus2760GeVAddedSig =         (TH1D*)directoryPythiaPlusPi02760GeV->Get("Pi0_Input_Reweighted_AddedSig");
   TH1D* histoMCPi0YieldPtPythiaPlus2760GeVAddedSigWOWeights = (TH1D*)directoryPythiaPlusPi02760GeV->Get("Pi0_Input_AddedSig");

   TH1D* histoPi0RawYieldPythiaPlus2760GeV =          (TH1D*)directoryPythiaPlusPi02760GeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPythiaPlus2760GeV->Scale(1000.);
   histoPi0MassMCPythiaPlus2760GeV->Scale(1000.);

   
   TDirectory* directoryPythiaPlusEta2760GeV =  (TDirectory*)filePythiaPlus->Get("Eta_pp_2760GeV"); 
   
   TH1D* histoEtaCorrectedSpecPythiaPlus2760GeV =   (TH1D*)directoryPythiaPlusEta2760GeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPythiaPlus2760GeV =        (TGraphAsymmErrors*)directoryPythiaPlusEta2760GeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPythiaPlus2760GeV =        (TH1D*)directoryPythiaPlusEta2760GeV->Get("MassEta");
   TH1D* histoEtaMassMCPythiaPlus2760GeV =          (TH1D*)directoryPythiaPlusEta2760GeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPythiaPlus2760GeV =       (TH1D*)directoryPythiaPlusEta2760GeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPythiaPlus2760GeV =         (TH1D*)directoryPythiaPlusEta2760GeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPythiaPlus2760GeV =        (TH1D*)directoryPythiaPlusEta2760GeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPythiaPlus2760GeV =          (TH1D*)directoryPythiaPlusEta2760GeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPythiaPlus2760GeV =         (TH1D*)directoryPythiaPlusEta2760GeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPythiaPlus2760GeVWOWeights = (TH1D*)directoryPythiaPlusEta2760GeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPythiaPlus2760GeV =        (TH1D*)directoryPythiaPlusEta2760GeV->Get("Eta_Weights");
   TH1D* histoMCEtaYieldPtPythiaPlus2760GeVAddedSig =         (TH1D*)directoryPythiaPlusEta2760GeV->Get("Eta_Input_Reweighted_AddedSig");
   TH1D* histoMCEtaYieldPtPythiaPlus2760GeVAddedSigWOWeights = (TH1D*)directoryPythiaPlusEta2760GeV->Get("Eta_Input_AddedSig");

   TH1D* histoEtaRawYieldPythiaPlus2760GeV =          (TH1D*)directoryPythiaPlusEta2760GeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPythiaPlus2760GeV->Scale(1000.);
   histoEtaMassMCPythiaPlus2760GeV->Scale(1000.);
   
   
   //*********  0.9 TeV   *****************************************************************************************************
   TDirectory* directoryPythiaPlusPi0900GeV =  (TDirectory*)filePythiaPlus->Get("Pi0_pp_900GeV"); 
   
   TH1D* histoPi0CorrectedSpecPythiaPlus900GeV =   (TH1D*)directoryPythiaPlusPi0900GeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPythiaPlus900GeV =        (TGraphAsymmErrors*)directoryPythiaPlusPi0900GeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPythiaPlus900GeV =        (TH1D*)directoryPythiaPlusPi0900GeV->Get("MassPi0");
   TH1D* histoPi0MassMCPythiaPlus900GeV =          (TH1D*)directoryPythiaPlusPi0900GeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPythiaPlus900GeV =       (TH1D*)directoryPythiaPlusPi0900GeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPythiaPlus900GeV =         (TH1D*)directoryPythiaPlusPi0900GeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPythiaPlus900GeV =        (TH1D*)directoryPythiaPlusPi0900GeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPythiaPlus900GeV =          (TH1D*)directoryPythiaPlusPi0900GeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPythiaPlus900GeV =         (TH1D*)directoryPythiaPlusPi0900GeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPythiaPlus900GeVWOWeights = (TH1D*)directoryPythiaPlusPi0900GeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPythiaPlus900GeV =        (TH1D*)directoryPythiaPlusPi0900GeV->Get("Pi0_Weights");
   TH1D* histoMCPi0YieldPtPythiaPlus900GeVAddedSig =         (TH1D*)directoryPythiaPlusPi0900GeV->Get("Pi0_Input_Reweighted_AddedSig");
   TH1D* histoMCPi0YieldPtPythiaPlus900GeVAddedSigWOWeights = (TH1D*)directoryPythiaPlusPi0900GeV->Get("Pi0_Input_AddedSig");

   TH1D* histoPi0RawYieldPythiaPlus900GeV =          (TH1D*)directoryPythiaPlusPi0900GeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPythiaPlus900GeV->Scale(1000.);
   histoPi0MassMCPythiaPlus900GeV->Scale(1000.);

   
   TDirectory* directoryPythiaPlusEta900GeV =  (TDirectory*)filePythiaPlus->Get("Eta_pp_900GeV"); 
   
   TH1D* histoEtaCorrectedSpecPythiaPlus900GeV =   (TH1D*)directoryPythiaPlusEta900GeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPythiaPlus900GeV =        (TGraphAsymmErrors*)directoryPythiaPlusEta900GeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPythiaPlus900GeV =        (TH1D*)directoryPythiaPlusEta900GeV->Get("MassEta");
   TH1D* histoEtaMassMCPythiaPlus900GeV =          (TH1D*)directoryPythiaPlusEta900GeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPythiaPlus900GeV =       (TH1D*)directoryPythiaPlusEta900GeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPythiaPlus900GeV =         (TH1D*)directoryPythiaPlusEta900GeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPythiaPlus900GeV =        (TH1D*)directoryPythiaPlusEta900GeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPythiaPlus900GeV =          (TH1D*)directoryPythiaPlusEta900GeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPythiaPlus900GeV =         (TH1D*)directoryPythiaPlusEta900GeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPythiaPlus900GeVWOWeights = (TH1D*)directoryPythiaPlusEta900GeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPythiaPlus900GeV =        (TH1D*)directoryPythiaPlusEta900GeV->Get("Eta_Weights");
   TH1D* histoMCEtaYieldPtPythiaPlus900GeVAddedSig =         (TH1D*)directoryPythiaPlusEta900GeV->Get("Eta_Input_Reweighted_AddedSig");
   TH1D* histoMCEtaYieldPtPythiaPlus900GeVAddedSigWOWeights = (TH1D*)directoryPythiaPlusEta900GeV->Get("Eta_Input_AddedSig");

   TH1D* histoEtaRawYieldPythiaPlus900GeV =          (TH1D*)directoryPythiaPlusEta900GeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPythiaPlus900GeV->Scale(1000.);
   histoEtaMassMCPythiaPlus900GeV->Scale(1000.);
   
   
   
   //***************************************************************************
   //**************************   PHOJET   *************************************
      TFile*   filePhojet =              new TFile(nameFilepp2.Data());

   //*********   8 TeV   *****************************************************************************************************
   TDirectory* directoryPhojetPi08TeV =  (TDirectory*)filePhojet->Get("Pi0_pp_8TeV"); 
   
   TH1D* histoPi0CorrectedSpecPhojet8TeV =   (TH1D*)directoryPhojetPi08TeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPhojet8TeV =        (TGraphAsymmErrors*)directoryPhojetPi08TeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPhojet8TeV =        (TH1D*)directoryPhojetPi08TeV->Get("MassPi0");
   TH1D* histoPi0MassMCPhojet8TeV =          (TH1D*)directoryPhojetPi08TeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPhojet8TeV =       (TH1D*)directoryPhojetPi08TeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPhojet8TeV =         (TH1D*)directoryPhojetPi08TeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPhojet8TeV =        (TH1D*)directoryPhojetPi08TeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPhojet8TeV =          (TH1D*)directoryPhojetPi08TeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPhojet8TeV =         (TH1D*)directoryPhojetPi08TeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPhojet8TeVWOWeights = (TH1D*)directoryPhojetPi08TeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPhojet8TeV =        (TH1D*)directoryPhojetPi08TeV->Get("Pi0_Weights");
   TH1D* histoPi0RawYieldPhojet8TeV =          (TH1D*)directoryPhojetPi08TeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPhojet8TeV->Scale(1000.);
   histoPi0MassMCPhojet8TeV->Scale(1000.);

   
   TDirectory* directoryPhojetEta8TeV =  (TDirectory*)filePhojet->Get("Eta_pp_8TeV"); 
   
   TH1D* histoEtaCorrectedSpecPhojet8TeV =   (TH1D*)directoryPhojetEta8TeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPhojet8TeV =        (TGraphAsymmErrors*)directoryPhojetEta8TeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPhojet8TeV =        (TH1D*)directoryPhojetEta8TeV->Get("MassEta");
   TH1D* histoEtaMassMCPhojet8TeV =          (TH1D*)directoryPhojetEta8TeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPhojet8TeV =       (TH1D*)directoryPhojetEta8TeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPhojet8TeV =         (TH1D*)directoryPhojetEta8TeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPhojet8TeV =        (TH1D*)directoryPhojetEta8TeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPhojet8TeV =          (TH1D*)directoryPhojetEta8TeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPhojet8TeV =         (TH1D*)directoryPhojetEta8TeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPhojet8TeVWOWeights = (TH1D*)directoryPhojetEta8TeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPhojet8TeV =        (TH1D*)directoryPhojetEta8TeV->Get("Eta_Weights");
   TH1D* histoEtaRawYieldPhojet8TeV =          (TH1D*)directoryPhojetEta8TeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPhojet8TeV->Scale(1000.);
   histoEtaMassMCPhojet8TeV->Scale(1000.);
   
   
   //*********   7 TeV   *****************************************************************************************************
   TDirectory* directoryPhojetPi07TeV =  (TDirectory*)filePhojet->Get("Pi0_pp_7TeV"); 
   
   TH1D* histoPi0CorrectedSpecPhojet7TeV =   (TH1D*)directoryPhojetPi07TeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPhojet7TeV =        (TGraphAsymmErrors*)directoryPhojetPi07TeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPhojet7TeV =        (TH1D*)directoryPhojetPi07TeV->Get("MassPi0");
   TH1D* histoPi0MassMCPhojet7TeV =          (TH1D*)directoryPhojetPi07TeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPhojet7TeV =       (TH1D*)directoryPhojetPi07TeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPhojet7TeV =         (TH1D*)directoryPhojetPi07TeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPhojet7TeV =        (TH1D*)directoryPhojetPi07TeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPhojet7TeV =          (TH1D*)directoryPhojetPi07TeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPhojet7TeV =         (TH1D*)directoryPhojetPi07TeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPhojet7TeVWOWeights = (TH1D*)directoryPhojetPi07TeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPhojet7TeV =        (TH1D*)directoryPhojetPi07TeV->Get("Pi0_Weights");
   TH1D* histoPi0RawYieldPhojet7TeV =          (TH1D*)directoryPhojetPi07TeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPhojet7TeV->Scale(1000.);
   histoPi0MassMCPhojet7TeV->Scale(1000.);

   
   TDirectory* directoryPhojetEta7TeV =  (TDirectory*)filePhojet->Get("Eta_pp_7TeV"); 
   
   TH1D* histoEtaCorrectedSpecPhojet7TeV =   (TH1D*)directoryPhojetEta7TeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPhojet7TeV =        (TGraphAsymmErrors*)directoryPhojetEta7TeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPhojet7TeV =        (TH1D*)directoryPhojetEta7TeV->Get("MassEta");
   TH1D* histoEtaMassMCPhojet7TeV =          (TH1D*)directoryPhojetEta7TeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPhojet7TeV =       (TH1D*)directoryPhojetEta7TeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPhojet7TeV =         (TH1D*)directoryPhojetEta7TeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPhojet7TeV =        (TH1D*)directoryPhojetEta7TeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPhojet7TeV =          (TH1D*)directoryPhojetEta7TeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPhojet7TeV =         (TH1D*)directoryPhojetEta7TeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPhojet7TeVWOWeights = (TH1D*)directoryPhojetEta7TeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPhojet7TeV =        (TH1D*)directoryPhojetEta7TeV->Get("Eta_Weights");
   TH1D* histoEtaRawYieldPhojet7TeV =          (TH1D*)directoryPhojetEta7TeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPhojet7TeV->Scale(1000.);
   histoEtaMassMCPhojet7TeV->Scale(1000.);
   
   
   //*********   2.76 TeV   *****************************************************************************************************
   TDirectory* directoryPhojetPi02760GeV =  (TDirectory*)filePhojet->Get("Pi0_pp_2760GeV"); 
   
   TH1D* histoPi0CorrectedSpecPhojet2760GeV =   (TH1D*)directoryPhojetPi02760GeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPhojet2760GeV =        (TGraphAsymmErrors*)directoryPhojetPi02760GeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPhojet2760GeV =        (TH1D*)directoryPhojetPi02760GeV->Get("MassPi0");
   TH1D* histoPi0MassMCPhojet2760GeV =          (TH1D*)directoryPhojetPi02760GeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPhojet2760GeV =       (TH1D*)directoryPhojetPi02760GeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPhojet2760GeV =         (TH1D*)directoryPhojetPi02760GeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPhojet2760GeV =        (TH1D*)directoryPhojetPi02760GeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPhojet2760GeV =          (TH1D*)directoryPhojetPi02760GeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPhojet2760GeV =         (TH1D*)directoryPhojetPi02760GeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPhojet2760GeVWOWeights = (TH1D*)directoryPhojetPi02760GeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPhojet2760GeV =        (TH1D*)directoryPhojetPi02760GeV->Get("Pi0_Weights");
   TH1D* histoPi0RawYieldPhojet2760GeV =          (TH1D*)directoryPhojetPi02760GeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPhojet2760GeV->Scale(1000.);
   histoPi0MassMCPhojet2760GeV->Scale(1000.);

   
   TDirectory* directoryPhojetEta2760GeV =  (TDirectory*)filePhojet->Get("Eta_pp_2760GeV"); 
   
   TH1D* histoEtaCorrectedSpecPhojet2760GeV =   (TH1D*)directoryPhojetEta2760GeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPhojet2760GeV =        (TGraphAsymmErrors*)directoryPhojetEta2760GeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPhojet2760GeV =        (TH1D*)directoryPhojetEta2760GeV->Get("MassEta");
   TH1D* histoEtaMassMCPhojet2760GeV =          (TH1D*)directoryPhojetEta2760GeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPhojet2760GeV =       (TH1D*)directoryPhojetEta2760GeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPhojet2760GeV =         (TH1D*)directoryPhojetEta2760GeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPhojet2760GeV =        (TH1D*)directoryPhojetEta2760GeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPhojet2760GeV =          (TH1D*)directoryPhojetEta2760GeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPhojet2760GeV =         (TH1D*)directoryPhojetEta2760GeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPhojet2760GeVWOWeights = (TH1D*)directoryPhojetEta2760GeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPhojet2760GeV =        (TH1D*)directoryPhojetEta2760GeV->Get("Eta_Weights");
   TH1D* histoEtaRawYieldPhojet2760GeV =          (TH1D*)directoryPhojetEta2760GeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPhojet2760GeV->Scale(1000.);
   histoEtaMassMCPhojet2760GeV->Scale(1000.);
   
   
   //*********  0.9 TeV   *****************************************************************************************************
   TDirectory* directoryPhojetPi0900GeV =  (TDirectory*)filePhojet->Get("Pi0_pp_900GeV"); 
   
   TH1D* histoPi0CorrectedSpecPhojet900GeV =   (TH1D*)directoryPhojetPi0900GeV->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPi0CorrectedSpecSysPhojet900GeV =        (TGraphAsymmErrors*)directoryPhojetPi0900GeV->Get("Pi0SystError"); 
   TH1D* histoPi0MassDataPhojet900GeV =        (TH1D*)directoryPhojetPi0900GeV->Get("MassPi0");
   TH1D* histoPi0MassMCPhojet900GeV =          (TH1D*)directoryPhojetPi0900GeV->Get("TrueMassPi0");
   TH1D* histoPi0WidthDataPhojet900GeV =       (TH1D*)directoryPhojetPi0900GeV->Get("FWHMPi0MeV");
   TH1D* histoPi0WidthMCPhojet900GeV =         (TH1D*)directoryPhojetPi0900GeV->Get("TrueFWHMPi0MeV");
   TH1D* histoPi0TrueEffiPtPhojet900GeV =        (TH1D*)directoryPhojetPi0900GeV->Get("EfficiencyPi0");
   TH1D* histoPi0AcceptPtPhojet900GeV =          (TH1D*)directoryPhojetPi0900GeV->Get("AcceptancePi0");
   TH1D* histoMCPi0YieldPtPhojet900GeV =         (TH1D*)directoryPhojetPi0900GeV->Get("Pi0_Input_Reweighted");
   TH1D* histoMCPi0YieldPtPhojet900GeVWOWeights = (TH1D*)directoryPhojetPi0900GeV->Get("Pi0_Input");
   TH1D* histoPi0WeightsPhojet900GeV =        (TH1D*)directoryPhojetPi0900GeV->Get("Pi0_Weights");
   TH1D* histoPi0RawYieldPhojet900GeV =          (TH1D*)directoryPhojetPi0900GeV->Get("RAWYieldPerEventsPi0");

   histoPi0MassDataPhojet900GeV->Scale(1000.);
   histoPi0MassMCPhojet900GeV->Scale(1000.);

   
   TDirectory* directoryPhojetEta900GeV =  (TDirectory*)filePhojet->Get("Eta_pp_900GeV"); 
   
   TH1D* histoEtaCorrectedSpecPhojet900GeV =   (TH1D*)directoryPhojetEta900GeV->Get("CorrectedYieldEta");   
   TGraphAsymmErrors* graphEtaCorrectedSpecSysPhojet900GeV =        (TGraphAsymmErrors*)directoryPhojetEta900GeV->Get("EtaSystError"); 
   TH1D* histoEtaMassDataPhojet900GeV =        (TH1D*)directoryPhojetEta900GeV->Get("MassEta");
   TH1D* histoEtaMassMCPhojet900GeV =          (TH1D*)directoryPhojetEta900GeV->Get("TrueMassEta");
   TH1D* histoEtaWidthDataPhojet900GeV =       (TH1D*)directoryPhojetEta900GeV->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMCPhojet900GeV =         (TH1D*)directoryPhojetEta900GeV->Get("TrueFWHMEtaMeV");
   TH1D* histoEtaTrueEffiPtPhojet900GeV =        (TH1D*)directoryPhojetEta900GeV->Get("EfficiencyEta");
   TH1D* histoEtaAcceptPtPhojet900GeV =          (TH1D*)directoryPhojetEta900GeV->Get("EfficiencyEta");
   TH1D* histoMCEtaYieldPtPhojet900GeV =         (TH1D*)directoryPhojetEta900GeV->Get("Eta_Input_Reweighted");
   TH1D* histoMCEtaYieldPtPhojet900GeVWOWeights = (TH1D*)directoryPhojetEta900GeV->Get("Eta_Input");
   TH1D* histoEtaWeightsPhojet900GeV =        (TH1D*)directoryPhojetEta900GeV->Get("Eta_Weights");
   TH1D* histoEtaRawYieldPhojet900GeV =          (TH1D*)directoryPhojetEta900GeV->Get("RAWYieldPerEventsEta");

   histoEtaMassDataPhojet900GeV->Scale(1000.);
   histoEtaMassMCPhojet900GeV->Scale(1000.);
   
   
   
   

   
   TFile* fileNeutralPionDataPP = new TFile(nameFilePPOld.Data());
   
   //*************    7 TeV    *****************
   TDirectory* directoryPi07TeV =      (TDirectory*)fileNeutralPionDataPP->Get("Pi07TeV"); 
   TH1D* histoAccPi07TeV =             (TH1D*)directoryPi07TeV->Get("AcceptancePi0");
   TH1D* histoTrueEffPtPi07TeV =       (TH1D*)directoryPi07TeV->Get("EfficiencyPi0");
   TH1D* histoRawYieldPi07TeV =        (TH1D*)directoryPi07TeV->Get("RAWYieldPerEventsPi0");
   TDirectory* directoryEta7TeV =      (TDirectory*)fileNeutralPionDataPP->Get("Eta7TeV");
   TH1D* histoCorrectedYieldEta7TeV = (TH1D*)directoryEta7TeV->Get("CorrectedYieldEta");
   TGraphAsymmErrors* graphCorrectedYieldSysEta7TeV =   (TGraphAsymmErrors*)directoryEta7TeV->Get("EtaSystError");
   TH1D* histoAccEta7TeV =             (TH1D*)directoryEta7TeV->Get("AcceptanceEta");
   TH1D* histoTrueEffPtEta7TeV =       (TH1D*)directoryEta7TeV->Get("EfficiencyEta");
   TH1D* histoRawYieldEta7TeV =        (TH1D*)directoryEta7TeV->Get("RAWYieldPerEventsEta");
   TH1D* histoEtaMassData7TeV =        (TH1D*)directoryEta7TeV->Get("MassEta");
   TH1D* histoEtaWidthData7TeV =       (TH1D*)directoryEta7TeV->Get("FWHMEtaMeV");
   TH1D* histoEtaMassMC7TeV =          (TH1D*)directoryEta7TeV->Get("TrueMassEta");
   TH1D* histoEtaWidthMC7TeV =         (TH1D*)directoryEta7TeV->Get("TrueFWHMEtaMeV");
   TH1D* histoRatioEtaPi07TeV=        (TH1D*)directoryEta7TeV->Get("EtatoPi0RatioConversion");   
   TGraphAsymmErrors* graphRatioEtaPi0SystErr7TeV=             (TGraphAsymmErrors*)directoryEta7TeV->Get("EtatoPi0RatioConversionSys");

   histoEtaMassData7TeV->Scale(1000.);
   histoEtaMassMC7TeV->Scale(1000.);
   histoEtaMassData7TeV->SetBinContent(histoEtaMassData7TeV->GetNbinsX(),0);
   histoEtaMassMC7TeV->SetBinContent(histoEtaMassMC7TeV->GetNbinsX(),0);
   histoEtaWidthData7TeV->SetBinContent(histoEtaWidthData7TeV->GetNbinsX(),10000.);
   histoEtaWidthMC7TeV->SetBinContent(histoEtaWidthMC7TeV->GetNbinsX(),10000.);

   //*************    0.9 TeV    *****************
   TDirectory* directoryPi0900GeV =    (TDirectory*)fileNeutralPionDataPP->Get("Pi0900GeV"); 
   TH1D* histoAccPi0900GeV =           (TH1D*)directoryPi0900GeV->Get("AcceptancePi0");
   TH1D* histoTrueEffPtPi0900GeV =  (TH1D*)directoryPi0900GeV->Get("EfficiencyPi0");
   TH1D* histoRawYieldPi0900GeV =      (TH1D*)directoryPi0900GeV->Get("RAWYieldPerEventsPi0");
   TDirectory* directoryEta900GeV =    (TDirectory*)fileNeutralPionDataPP->Get("Eta900GeV"); 
   TH1D* histoCorrectedYieldEta900GeV =     (TH1D*)directoryEta900GeV->Get("CorrectedYieldEta");
   TGraphAsymmErrors* graphCorrectedYieldSysEta900GeV =  (TGraphAsymmErrors*)directoryEta900GeV->Get("EtaSystError");
   TH1D* histoAccEta900GeV =           (TH1D*)directoryEta900GeV->Get("AcceptanceEta");
   TH1D* histoTrueEffPtEta900GeV =  (TH1D*)directoryEta900GeV->Get("EfficiencyEta");
   TH1D* histoRawYieldEta900GeV =      (TH1D*)directoryEta900GeV->Get("RAWYieldPerEventsEta");
   TH1D* histoRatioEtaPi0900GeV=      (TH1D*)directoryEta900GeV->Get("EtatoPi0RatioConversion");   
 
   //*************    2.76 TeV    *****************
   TGraphAsymmErrors* graphRatioEtaPi0SystErr900GeV=           (TGraphAsymmErrors*)directoryEta900GeV->Get("EtatoPi0RatioConversionSys");
   TDirectory* directoryPi02760GeV =   (TDirectory*)fileNeutralPionDataPP->Get("Pi02.76TeV"); 
   TH1D* histoAccPi02760GeV =          (TH1D*)directoryPi02760GeV->Get("AcceptancePi0");
   TH1D* histoTrueEffPtPi02760GeV =    (TH1D*)directoryPi02760GeV->Get("EfficiencyPi0");
   TH1D* histoRawYieldPi02760GeV =  (TH1D*)directoryPi02760GeV->Get("RAWYieldPerEventsPi0");
   TH1D* histoCorrectedYieldPi02760GeV =  (TH1D*)directoryPi02760GeV->Get("CorrectedYieldPi0");
   TGraphAsymmErrors* graphCorrectedYieldSysPi02760GeV = (TGraphAsymmErrors*)directoryPi02760GeV->Get("Pi0SystError"); 
   TDirectory* directoryEta2760GeV =   (TDirectory*)fileNeutralPionDataPP->Get("Eta2.76TeV"); 
   TH1D* histoAccEta2760GeV =          (TH1D*)directoryEta2760GeV->Get("AcceptanceEta");
   TH1D* histoTrueEffPtEta2760GeV =    (TH1D*)directoryEta2760GeV->Get("EfficiencyEta");
   TH1D* histoCorrectedYieldEta2760GeV =  (TH1D*)directoryEta2760GeV->Get("CorrectedYieldEta");
   TGraphAsymmErrors* graphCorrectedYieldSysEta2760GeV = (TGraphAsymmErrors*)directoryEta2760GeV->Get("EtaSystError");
   TH1D* histoRawYieldEta2760GeV =  (TH1D*)directoryEta2760GeV->Get("RAWYieldPerEventsEta");
   TH1D* histoEtaMassData2760GeV =     (TH1D*)directoryEta2760GeV->Get("MassEta");
   TH1D* histoEtaWidthData2760GeV =    (TH1D*)directoryEta2760GeV->Get("FWHMEtaMeV");
   TH1D* histoEtaMassMC2760GeV =       (TH1D*)directoryEta2760GeV->Get("TrueMassEta");
   TH1D* histoEtaWidthMC2760GeV =      (TH1D*)directoryEta2760GeV->Get("TrueFWHMEtaMeV");
   TH1D* histoRatioEtaPi02760GeV=     (TH1D*)directoryEta2760GeV->Get("EtatoPi0RatioConversion");   
   TGraphAsymmErrors* graphRatioEtaPi0SystErr2760GeV=          (TGraphAsymmErrors*)directoryEta2760GeV->Get("EtatoPi0RatioConversionSys");

   histoEtaMassData2760GeV->Scale(1000.);
   histoEtaMassMC2760GeV->Scale(1000.);
   histoEtaMassData2760GeV->SetBinContent(histoEtaMassData2760GeV->GetNbinsX(),0);
   histoEtaMassMC2760GeV->SetBinContent(histoEtaMassMC2760GeV->GetNbinsX(),0);
   histoEtaWidthData2760GeV->SetBinContent(histoEtaWidthData2760GeV->GetNbinsX(),10000.);
   histoEtaWidthMC2760GeV->SetBinContent(histoEtaWidthMC2760GeV->GetNbinsX(),10000.);

   TH1D* histoPCMMassDataPP2760GeV =      (TH1D*)directoryPi02760GeV->Get("MassPi0");
   TH1D* histoPCMWidthDataPP2760GeV =     (TH1D*)directoryPi02760GeV->Get("FWHMPi0MeV");
   TH1D* histoPCMMassMCPP2760GeV =        (TH1D*)directoryPi02760GeV->Get("TrueMassPi0");
   TH1D* histoPCMWidthMCPP2760GeV =          (TH1D*)directoryPi02760GeV->Get("TrueFWHMPi0MeV");
   histoPCMMassDataPP2760GeV->Scale(1000.);
   histoPCMMassMCPP2760GeV->Scale(1000.);
   histoPCMMassDataPP2760GeV->SetBinContent(histoPCMMassDataPP2760GeV->GetNbinsX(),0);
   histoPCMMassMCPP2760GeV->SetBinContent(histoPCMMassMCPP2760GeV->GetNbinsX(),0);
   histoPCMWidthDataPP2760GeV->SetBinContent(histoPCMWidthDataPP2760GeV->GetNbinsX(),10000.);
   histoPCMWidthMCPP2760GeV->SetBinContent(histoPCMWidthMCPP2760GeV->GetNbinsX(),10000.);


   cout << "here" << endl;
   
   
   //********************* processing old pp files *************************************
/*
   TFile*   fileCocktail =                new TFile("ExternalInputPbPb/cocktail_allCentpluspp.root");
   TDirectory* directoryCocktailpp2760GeV =           (TDirectory*)fileCocktail->Get("cocktail_pp_2760GeV_qcd"); 
   TH1D* histoEtaFromCocktailpp2760GeV = (TH1D*)directoryCocktailpp2760GeV->Get("ptEta"); 

   
   TFile* fileNeutralPionCombDataPP = new TFile("CombinedResultsPaperX_18_Feb_2014.root");
   TGraphAsymmErrors* graphInvYieldPi0Comb7TeV= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeV");
   graphInvYieldPi0Comb7TeV = ScaleGraph(graphInvYieldPi0Comb7TeV,1./xSection7TeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0Comb7TeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeVStatErr");
   graphInvYieldPi0Comb7TeVStatErr = ScaleGraph(graphInvYieldPi0Comb7TeVStatErr,1./xSection7TeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0Comb7TeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeVSysErr");
   graphInvYieldPi0Comb7TeVSysErr = ScaleGraph(graphInvYieldPi0Comb7TeVSysErr,1./xSection7TeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PCM7TeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMStat7TeV");
   graphInvYieldPi0PCM7TeVStatErr = ScaleGraph(graphInvYieldPi0PCM7TeVStatErr,1./xSection7TeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PCM7TeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMSys7TeV");
   graphInvYieldPi0PCM7TeVSysErr = ScaleGraph(graphInvYieldPi0PCM7TeVSysErr,1./xSection7TeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PHOS7TeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSStat7TeV");
   graphInvYieldPi0PHOS7TeVStatErr = ScaleGraph(graphInvYieldPi0PHOS7TeVStatErr,1./xSection7TeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PHOS7TeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSSys7TeV");
   graphInvYieldPi0PHOS7TeVSysErr = ScaleGraph(graphInvYieldPi0PHOS7TeVSysErr,1./xSection7TeVppINEL);
   cout << "here 2.76TeV" << endl;
   TGraphAsymmErrors* graphInvYieldPi0Comb2760GeV= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb2760GeV");
   graphInvYieldPi0Comb2760GeV = ScaleGraph(graphInvYieldPi0Comb2760GeV,1./xSection2760GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0Comb2760GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb2760GeVStatErr");
   graphInvYieldPi0Comb2760GeVStatErr = ScaleGraph(graphInvYieldPi0Comb2760GeVStatErr,1./xSection2760GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0Comb2760GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb2760GeVSysErr");
   graphInvYieldPi0Comb2760GeVSysErr = ScaleGraph(graphInvYieldPi0Comb2760GeVSysErr,1./xSection2760GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PCM2760GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
   graphInvYieldPi0PCM2760GeVStatErr = ScaleGraph(graphInvYieldPi0PCM2760GeVStatErr,1./xSection2760GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PCM2760GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCM2760GeVSysErr");
   graphInvYieldPi0PCM2760GeVSysErr = ScaleGraph(graphInvYieldPi0PCM2760GeVSysErr,1./xSection2760GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PHOS2760GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOS2760GeVStatErr");
   graphInvYieldPi0PHOS2760GeVStatErr = ScaleGraph(graphInvYieldPi0PHOS2760GeVStatErr,1./xSection2760GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PHOS2760GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOS2760GeVSysErr");
   graphInvYieldPi0PHOS2760GeVSysErr = ScaleGraph(graphInvYieldPi0PHOS2760GeVSysErr,1./xSection2760GeVppINEL);
   cout << "here 900 GeV" << endl;
   TGraphAsymmErrors* graphInvYieldPi0Comb900GeV= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeV");
   graphInvYieldPi0Comb900GeV = ScaleGraph(graphInvYieldPi0Comb900GeV,1./xSection900GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0Comb900GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeVStatErr");
   graphInvYieldPi0Comb900GeVStatErr = ScaleGraph(graphInvYieldPi0Comb900GeVStatErr,1./xSection900GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0Comb900GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeVSysErr");
   graphInvYieldPi0Comb900GeVSysErr = ScaleGraph(graphInvYieldPi0Comb900GeVSysErr,1./xSection900GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PCM900GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMStat900GeV");
   graphInvYieldPi0PCM900GeVStatErr = ScaleGraph(graphInvYieldPi0PCM900GeVStatErr,1./xSection900GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PCM900GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMSys900GeV");
   graphInvYieldPi0PCM900GeVSysErr = ScaleGraph(graphInvYieldPi0PCM900GeVSysErr,1./xSection900GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PHOS900GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSStat900GeV");
   graphInvYieldPi0PHOS900GeVStatErr = ScaleGraph(graphInvYieldPi0PHOS900GeVStatErr,1./xSection900GeVppINEL);
   TGraphAsymmErrors* graphInvYieldPi0PHOS900GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSSys900GeV");
   graphInvYieldPi0PHOS900GeVSysErr = ScaleGraph(graphInvYieldPi0PHOS900GeVSysErr,1./xSection900GeVppINEL);



   cout << "efficiencies Added Sig 2.76TeV GeV" << endl;
   TFile* fileNeutralPion2760GeVDataPPEffiAddedSig = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistosAddSig_000001100209366300380000000_01631031009000.root");
   TH1D* histoEffiAddedSigPP = (TH1D*)fileNeutralPion2760GeVDataPPEffiAddedSig->Get("TrueMesonEffiPt");
   
   TFile* fileNeutralPion7TeVDataPPEffi2760GeVCut = new TFile("ppAdditionalInput/0000011002093663003800000_01631031009/7TeV/Pi0_MC_GammaConvV1CorrectionHistosD_0000011002093663003800000_01631031009.root");
   TH1D* histoEffi7TeV2760GeVCut = (TH1D*)fileNeutralPion7TeVDataPPEffi2760GeVCut->Get("TrueMesonEffiPt");
   
   TFile* fileNeutralPion2760GeVDataPPEffiWithSDD = new TFile("ppAdditionalInput/0002011002093663003800000_01631031009/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_0002011002093663003800000_01631031009.root");
   TH1D* histoEffi2760GeVWithSDD = (TH1D*)fileNeutralPion2760GeVDataPPEffiWithSDD->Get("TrueMesonEffiPt");
   
   cout << "efficiencies MinBias 2.76TeV GeV" << endl;
   TFile* fileNeutralPion2760GeVDataPPEffiMinBias = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistosMinBias_000001100209366300380000000_01631031009000.root");
   TH1D* histoEffiMinBiasPP = (TH1D*)fileNeutralPion2760GeVDataPPEffiMinBias->Get("TrueMesonEffiPt");

   
   cout << "MC spectra 2.76TeV GeV" << endl;
   TFile* fileNeutralPion2760GeVDataPP = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Pi0_data_GammaConvV1Correction_000001100209366300380000000_01631031009000.root");
   TH1D* histoMCYieldPi02760GeVFinal = (TH1D*)fileNeutralPion2760GeVDataPP->Get("MCYield_Meson_oldBin");
   TFile* fileEta2760GeVDataPP = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Eta_data_GammaConvV1Correction_000001100209366300380000000_01631031009000.root");
   TH1D* histoMCYieldEta2760GeVFinal = (TH1D*)fileEta2760GeVDataPP->Get("MCYield_Meson_oldBin");
   
   cout << "MC spectra added signal 2.76TeV GeV" << endl;
   TFile* fileNeutralPion2760GeVDataPPAddSig = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistosAddSig_000001100209366300380000000_01631031009000.root");
   TH1D* histoMCYieldPi02760GeVAddSig = (TH1D*)fileNeutralPion2760GeVDataPPAddSig->Get("MC_Meson_genPt_oldBin");
   TH1F *histoEventQualityMC =         (TH1F*)fileNeutralPion2760GeVDataPPAddSig->Get("NEvents");
   Float_t nEvtMC = GetNEvents(histoEventQualityMC);
   TString rapidityRange = "";
   Double_t deltaRapid =  ReturnRapidityStringAndDouble("01631031009", rapidityRange);
   Double_t scaling = 1./(2.*TMath::Pi());
   ScaleMCYield(histoMCYieldPi02760GeVAddSig,  deltaRapid,  scaling,  nEvtMC,  "Pi0" ,"kFALSE");
   cout << "MC spectra withough weighting added signal 2.76TeV GeV" << endl;
   TFile* fileNeutralPion2760GeVDataPPAddSigWOWeighting = new TFile("ppAdditionalInput/0000012002093663003800000_01631031009/2.76TeV/Pi0_data_GammaConvV1Correction_0000012002093663003800000_01631031009.root");
   TH1D* histoMCYieldPi02760GeVAddSigWOWeighting = (TH1D*)fileNeutralPion2760GeVDataPPAddSigWOWeighting->Get("MCYield_Meson_oldBin");
   
   cout << "MC spectra without weighting 2.76TeV GeV" << endl;
   TFile* fileNeutralPion2760GeVDataPPWOWeigthing = new TFile("ppAdditionalInput/0000011002093663003800000_01631031009/2.76TeV/Pi0_data_GammaConvV1Correction_0000011002093663003800000_01631031009.root");
   TH1D* histoMCYieldPi02760GeVWOWeighting = (TH1D*)fileNeutralPion2760GeVDataPPWOWeigthing->Get("MCYield_Meson_oldBin");
   cout << "MC spectra 7TeV GeV" << endl;
   TFile* fileNeutralPion7TeVDataPP = new TFile("900366208010033211360000000900/7TeV/Pi0_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
   TH1D* histoMCYieldPi07TeV = (TH1D*)fileNeutralPion7TeVDataPP->Get("MCYield_Meson_oldBin");
   TFile* fileEta7TeVDataPP = new TFile("900366208010033211360000000900/7TeV/Eta_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
   TH1D* histoMCYieldEta7TeV = (TH1D*)fileEta7TeVDataPP->Get("MCYield_Meson_oldBin");
   cout << "MC spectra 0.9TeV GeV" << endl;
   TFile* fileNeutralPion900GeVDataPP = new TFile("900366208010033211360000000900/900GeV/Pi0_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
   TH1D* histoMCYieldPi0900GeV = (TH1D*)fileNeutralPion900GeVDataPP->Get("MCYield_Meson_oldBin");
   TFile* fileEta900GeVDataPP = new TFile("900366208010033211360000000900/900GeV/Eta_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
   TH1D* histoMCYieldEta900GeV = (TH1D*)fileEta900GeVDataPP->Get("MCYield_Meson_oldBin");
*/

   TFile* filePi0Pythia = new TFile("0000011_002092570028250400000_01521065000000/8TeV/Pi0_data_GammaConvV1Correction_0000011_002092570028250400000_01521065000000.root");
   TH1D* histoMCYieldPi0Pythia = (TH1D*)filePi0Pythia->Get("MCYield_Meson_oldBin");
   TFile* filePi0PythiaPlus = new TFile("0000012_002092570028250400000_01521065000000/8TeV/Pi0_MC_GammaConvV1Correction_0000012_002092570028250400000_01521065000000.root");
   TH1D* histoMCYieldPi0PythiaPlus = (TH1D*)filePi0PythiaPlus->Get("MCYield_Meson_oldBin");

   TFile* fileEtaPythia = new TFile("0000011_002092570028250400000_01521065000000/8TeV/Eta_data_GammaConvV1Correction_0000011_002092570028250400000_01521065000000.root");
   TH1D* histoMCYieldEtaPythia = (TH1D*)fileEtaPythia->Get("MCYield_Meson_oldBin");
   TFile* fileEtaPythiaPlus = new TFile("0000012_002092570028250400000_01521065000000/8TeV/Eta_MC_GammaConvV1Correction_0000012_002092570028250400000_01521065000000.root");
   TH1D* histoMCYieldEtaPythiaPlus = (TH1D*)fileEtaPythiaPlus->Get("MCYield_Meson_oldBin");


   
   Width_t  widthLinesBoxes            = 1.4;
   Width_t  widthCommonFit             = 2.;
   Width_t  widthStatErrBars           = 1.5;
   Width_t  widthCommonErrors          = 1.1;
   Width_t  widthCommonSpectrumBoxes         = 0.99;
   if (suffix.CompareTo("eps")==0){
      widthLinesBoxes            = 1.4;
      widthCommonFit             = 2.;
      widthStatErrBars           = 1.5;
      widthCommonErrors          = 1.1;
      widthCommonSpectrumBoxes         = 0.99;
   } else {
      widthLinesBoxes            = 2.3;
      widthCommonFit             = 2.6;
      widthStatErrBars           = 2.6;
      widthCommonErrors          = 2.;
      widthCommonSpectrumBoxes         = 2.3;
   }


   cout << "Pi0 efficiency in different system - Pythia" << endl;   
   
   TCanvas* canvasEfficiencyPi0Pythia = new TCanvas("canvasEfficiencyPi0Pythia","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyPi0Pythia, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffPi0Pythia;
   histo2DEffPi0Pythia = new TH2F("histo2DEffPi0Pythia","histo2DEffPi0Pythia",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffPi0Pythia, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffPi0Pythia->Draw("copy");

//   DrawGammaSetMarker(histoTrueEffPtPi07TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
//   histoTrueEffPtPi07TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueEffPtPi02760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
//   histoTrueEffPtPi02760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueEffPtPi0900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
//   histoTrueEffPtPi0900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0TrueEffiPtPythia8TeV, markerStyleSpectrum8TeVMC, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
   histoPi0TrueEffiPtPythia8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0TrueEffiPtPythia7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
   histoPi0TrueEffiPtPythia7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0TrueEffiPtPythia2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
   histoPi0TrueEffiPtPythia2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoPi0TrueEffiPtPythia900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0TrueEffiPtPythia900GeV->DrawCopy("e1,same");    

   TLegend* legendEffiPi0Pythia = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiPi0Pythia->SetFillColor(0);
   legendEffiPi0Pythia->SetLineColor(0);
   legendEffiPi0Pythia->SetTextSize(0.027);
   legendEffiPi0Pythia->SetNColumns(2);
   legendEffiPi0Pythia->AddEntry(histoPi0TrueEffiPtPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendEffiPi0Pythia->AddEntry(histoPi0TrueEffiPtPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendEffiPi0Pythia->AddEntry(histoPi0TrueEffiPtPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendEffiPi0Pythia->AddEntry(histoPi0TrueEffiPtPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendEffiPi0Pythia->AddEntry(histoTrueEffPtPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendEffiPi0Pythia->AddEntry(histoTrueEffPtPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendEffiPi0Pythia->AddEntry(histoTrueEffPtPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendEffiPi0Pythia->Draw();
   
   canvasEfficiencyPi0Pythia->SaveAs(Form("%s/EfficiencyPi0Pythia.%s",outputDir.Data(),suffix.Data()));
   histo2DEffPi0Pythia->Draw("copy");
   
   
   
   cout << "Pi0 efficiency in different system - Pythia with added signals" << endl;   
   
   TCanvas* canvasEfficiencyPi0PythiaPlus = new TCanvas("canvasEfficiencyPi0PythiaPlus","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyPi0PythiaPlus, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffPi0PythiaPlus;
   histo2DEffPi0PythiaPlus = new TH2F("histo2DEffPi0PythiaPlus","histo2DEffPi0PythiaPlus",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffPi0PythiaPlus, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffPi0PythiaPlus->Draw("copy");

//   DrawGammaSetMarker(histoTrueEffPtPi07TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
//   histoTrueEffPtPi07TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueEffPtPi02760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
//   histoTrueEffPtPi02760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueEffPtPi0900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
//   histoTrueEffPtPi0900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0TrueEffiPtPythiaPlus8TeV, markerStyleSpectrum8TeVMC, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
   histoPi0TrueEffiPtPythiaPlus8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0TrueEffiPtPythiaPlus7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
   histoPi0TrueEffiPtPythiaPlus7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0TrueEffiPtPythiaPlus2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
   histoPi0TrueEffiPtPythiaPlus2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoPi0TrueEffiPtPythiaPlus900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0TrueEffiPtPythiaPlus900GeV->DrawCopy("e1,same");    

   TLegend* legendEffiPi0PythiaPlus = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiPi0PythiaPlus->SetFillColor(0);
   legendEffiPi0PythiaPlus->SetLineColor(0);
   legendEffiPi0PythiaPlus->SetTextSize(0.027);
   legendEffiPi0PythiaPlus->SetNColumns(2);
   legendEffiPi0PythiaPlus->AddEntry(histoPi0TrueEffiPtPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendEffiPi0PythiaPlus->AddEntry(histoPi0TrueEffiPtPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendEffiPi0PythiaPlus->AddEntry(histoPi0TrueEffiPtPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendEffiPi0PythiaPlus->AddEntry(histoPi0TrueEffiPtPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendEffiPi0PythiaPlus->AddEntry(histoTrueEffPtPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendEffiPi0PythiaPlus->AddEntry(histoTrueEffPtPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendEffiPi0PythiaPlus->AddEntry(histoTrueEffPtPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendEffiPi0PythiaPlus->Draw();
   
   canvasEfficiencyPi0PythiaPlus->SaveAs(Form("%s/EfficiencyPi0PythiaPlus.%s",outputDir.Data(),suffix.Data()));
   histo2DEffPi0PythiaPlus->Draw("copy");
   
   
   
   cout << "Pi0 efficiency in different system - Phojet" << endl;   
   
   TCanvas* canvasEfficiencyPi0Phojet = new TCanvas("canvasEfficiencyPi0Phojet","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyPi0Phojet, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffPi0Phojet;
   histo2DEffPi0Phojet = new TH2F("histo2DEffPi0Phojet","histo2DEffPi0Phojet",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffPi0Phojet, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffPi0Phojet->Draw("copy");

//   DrawGammaSetMarker(histoTrueEffPtPi07TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
//   histoTrueEffPtPi07TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueEffPtPi02760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
//   histoTrueEffPtPi02760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueEffPtPi0900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
//   histoTrueEffPtPi0900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0TrueEffiPtPhojet8TeV, markerStyleSpectrum8TeVMC, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
   histoPi0TrueEffiPtPhojet8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0TrueEffiPtPhojet7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
   histoPi0TrueEffiPtPhojet7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0TrueEffiPtPhojet2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
   histoPi0TrueEffiPtPhojet2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoPi0TrueEffiPtPhojet900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0TrueEffiPtPhojet900GeV->DrawCopy("e1,same");    

   TLegend* legendEffiPi0Phojet = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiPi0Phojet->SetFillColor(0);
   legendEffiPi0Phojet->SetLineColor(0);
   legendEffiPi0Phojet->SetTextSize(0.027);
   legendEffiPi0Phojet->SetNColumns(2);
   legendEffiPi0Phojet->AddEntry(histoPi0TrueEffiPtPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendEffiPi0Phojet->AddEntry(histoPi0TrueEffiPtPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendEffiPi0Phojet->AddEntry(histoPi0TrueEffiPtPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendEffiPi0Phojet->AddEntry(histoPi0TrueEffiPtPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendEffiPi0Phojet->AddEntry(histoTrueEffPtPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendEffiPi0Phojet->AddEntry(histoTrueEffPtPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendEffiPi0Phojet->AddEntry(histoTrueEffPtPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendEffiPi0Phojet->Draw();
   
   canvasEfficiencyPi0Phojet->SaveAs(Form("%s/EfficiencyPi0Phojet.%s",outputDir.Data(),suffix.Data()));
   histo2DEffPi0Phojet->Draw("copy");

   
   //******************************************************************************************************************************************************
   
   cout << "Pi0 efficiency 8TeV" << endl;   
   
   TCanvas* canvasEfficiencyPi08TeV = new TCanvas("canvasEfficiencyPi08TeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyPi08TeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffPi08TeV;
   histo2DEffPi08TeV = new TH2F("histo2DEffPi08TeV","histo2DEffPi08TeV",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffPi08TeV, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffPi08TeV->Draw("copy");

   DrawGammaSetMarker(histoPi0TrueEffiPtPythia8TeV, markerStyleSpectrum8TeVMC, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
   histoPi0TrueEffiPtPythia8TeV->DrawCopy("e1,same");   
//    DrawGammaSetMarker(histoPi0TrueEffiPtPythiaPlus8TeV, markerStyleSpectrum8TeVMC+1, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
//    histoPi0TrueEffiPtPythiaPlus8TeV->DrawCopy("e1,same");
//    DrawGammaSetMarker(histoPi0TrueEffiPtPhojet8TeV, markerStyleSpectrum8TeVMC+2, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
//    histoPi0TrueEffiPtPhojet8TeV->DrawCopy("e1,same");
   
   TLegend* legendEffiPi08TeV = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiPi08TeV->SetFillColor(0);
   legendEffiPi08TeV->SetLineColor(0);
   legendEffiPi08TeV->SetTextSize(0.027);
   legendEffiPi08TeV->SetNColumns(2);
   legendEffiPi08TeV->AddEntry(histoPi0TrueEffiPtPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia","p");
//    legendEffiPi08TeV->AddEntry(histoPi0TrueEffiPtPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia+AddSign","p");
//    legendEffiPi08TeV->AddEntry(histoPi0TrueEffiPtPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV - Phojet","p");
   legendEffiPi08TeV->Draw();
   
   canvasEfficiencyPi08TeV->SaveAs(Form("%s/EfficiencyPi08TeV.%s",outputDir.Data(),suffix.Data()));
   histo2DEffPi08TeV->Draw("copy");
   
 
   cout << "Pi0 efficiency 7TeV" << endl;   
   
   TCanvas* canvasEfficiencyPi07TeV = new TCanvas("canvasEfficiencyPi07TeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyPi07TeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffPi07TeV;
   histo2DEffPi07TeV = new TH2F("histo2DEffPi07TeV","histo2DEffPi07TeV",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffPi07TeV, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffPi07TeV->Draw("copy");

   DrawGammaSetMarker(histoPi0TrueEffiPtPythia7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV); 
   histoPi0TrueEffiPtPythia7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0TrueEffiPtPythiaPlus7TeV, markerStyleSpectrum7TeVMC+1, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV); 
   histoPi0TrueEffiPtPythiaPlus7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0TrueEffiPtPhojet7TeV, markerStyleSpectrum7TeVMC+2, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV); 
   histoPi0TrueEffiPtPhojet7TeV->DrawCopy("e1,same");
   
   TLegend* legendEffiPi07TeV = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiPi07TeV->SetFillColor(0);
   legendEffiPi07TeV->SetLineColor(0);
   legendEffiPi07TeV->SetTextSize(0.027);
   legendEffiPi07TeV->SetNColumns(2);
   legendEffiPi07TeV->AddEntry(histoPi0TrueEffiPtPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia","p");
   legendEffiPi07TeV->AddEntry(histoPi0TrueEffiPtPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia+AddSign","p");
   legendEffiPi07TeV->AddEntry(histoPi0TrueEffiPtPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV - Phojet","p");
   legendEffiPi07TeV->Draw();
   
   canvasEfficiencyPi07TeV->SaveAs(Form("%s/EfficiencyPi07TeV.%s",outputDir.Data(),suffix.Data()));
   histo2DEffPi07TeV->Draw("copy");
 
   
   cout << "Pi0 efficiency 2760GeV" << endl;   
   
   TCanvas* canvasEfficiencyPi02760GeV = new TCanvas("canvasEfficiencyPi02760GeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyPi02760GeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffPi02760GeV;
   histo2DEffPi02760GeV = new TH2F("histo2DEffPi02760GeV","histo2DEffPi02760GeV",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffPi02760GeV, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffPi02760GeV->Draw("copy");

   DrawGammaSetMarker(histoPi0TrueEffiPtPythia2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV); 
   histoPi0TrueEffiPtPythia2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0TrueEffiPtPythiaPlus2760GeV, markerStyleSpectrum2760GeVMC+1, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV); 
   histoPi0TrueEffiPtPythiaPlus2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0TrueEffiPtPhojet2760GeV, markerStyleSpectrum2760GeVMC+2, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV); 
   histoPi0TrueEffiPtPhojet2760GeV->DrawCopy("e1,same");
   
   TLegend* legendEffiPi02760GeV = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiPi02760GeV->SetFillColor(0);
   legendEffiPi02760GeV->SetLineColor(0);
   legendEffiPi02760GeV->SetTextSize(0.027);
   legendEffiPi02760GeV->SetNColumns(2);
   legendEffiPi02760GeV->AddEntry(histoPi0TrueEffiPtPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia","p");
   legendEffiPi02760GeV->AddEntry(histoPi0TrueEffiPtPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia+AddSign","p");
   legendEffiPi02760GeV->AddEntry(histoPi0TrueEffiPtPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Phojet","p");
   legendEffiPi02760GeV->Draw();
   
   canvasEfficiencyPi02760GeV->SaveAs(Form("%s/EfficiencyPi02760GeV.%s",outputDir.Data(),suffix.Data()));
   histo2DEffPi02760GeV->Draw("copy");
 
   
   cout << "Pi0 efficiency 900GeV" << endl;   
   
   TCanvas* canvasEfficiencyPi0900GeV = new TCanvas("canvasEfficiencyPi0900GeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyPi0900GeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffPi0900GeV;
   histo2DEffPi0900GeV = new TH2F("histo2DEffPi0900GeV","histo2DEffPi0900GeV",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffPi0900GeV->Draw("copy");

   DrawGammaSetMarker(histoPi0TrueEffiPtPythia900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0TrueEffiPtPythia900GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0TrueEffiPtPythiaPlus900GeV, markerStyleSpectrum900GeVMC+1, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0TrueEffiPtPythiaPlus900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0TrueEffiPtPhojet900GeV, markerStyleSpectrum900GeVMC+2, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0TrueEffiPtPhojet900GeV->DrawCopy("e1,same");
   
   TLegend* legendEffiPi0900GeV = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiPi0900GeV->SetFillColor(0);
   legendEffiPi0900GeV->SetLineColor(0);
   legendEffiPi0900GeV->SetTextSize(0.027);
   legendEffiPi0900GeV->SetNColumns(2);
   legendEffiPi0900GeV->AddEntry(histoPi0TrueEffiPtPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia","p");
   legendEffiPi0900GeV->AddEntry(histoPi0TrueEffiPtPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia+AddSign","p");
   legendEffiPi0900GeV->AddEntry(histoPi0TrueEffiPtPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Phojet","p");
   legendEffiPi0900GeV->Draw();
   
   canvasEfficiencyPi0900GeV->SaveAs(Form("%s/EfficiencyPi0900GeV.%s",outputDir.Data(),suffix.Data()));
   histo2DEffPi0900GeV->Draw("copy");
 
   
   
//**********************************************************************************************************************
   //!!!! histo names for old pp file may not be right !!!!!
   
   cout << "Pi0 acceptance in different system - Pythia" << endl;   
   
   TCanvas* canvasAcceptancePi0Pythia = new TCanvas("canvasAcceptancePi0Pythia","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptancePi0Pythia, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccPi0Pythia;
   histo2DAccPi0Pythia = new TH2F("histo2DAccPi0Pythia","histo2DAccPi0Pythia",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccPi0Pythia, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccPi0Pythia->Draw("copy");

//   DrawGammaSetMarker(histoTrueAccPtPi07TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
//   histoTrueAccPtPi07TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueAccPtPi02760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
//   histoTrueAccPtPi02760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueAccPtPi0900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
//   histoTrueAccPtPi0900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0AcceptPtPythia8TeV, markerStyleSpectrum8TeVMC, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
   histoPi0AcceptPtPythia8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0AcceptPtPythia7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
   histoPi0AcceptPtPythia7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0AcceptPtPythia2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
   histoPi0AcceptPtPythia2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoPi0AcceptPtPythia900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0AcceptPtPythia900GeV->DrawCopy("e1,same");    

   TLegend* legendAccPi0Pythia = new TLegend(0.34,0.13,0.93,0.43);
   legendAccPi0Pythia->SetFillColor(0);
   legendAccPi0Pythia->SetLineColor(0);
   legendAccPi0Pythia->SetTextSize(0.027);
   legendAccPi0Pythia->SetNColumns(2);
   legendAccPi0Pythia->AddEntry(histoPi0AcceptPtPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendAccPi0Pythia->AddEntry(histoPi0AcceptPtPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendAccPi0Pythia->AddEntry(histoPi0AcceptPtPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendAccPi0Pythia->AddEntry(histoPi0AcceptPtPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendAccPi0Pythia->AddEntry(histoTrueAccPtPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendAccPi0Pythia->AddEntry(histoTrueAccPtPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendAccPi0Pythia->AddEntry(histoTrueAccPtPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendAccPi0Pythia->Draw();
   
   canvasAcceptancePi0Pythia->SaveAs(Form("%s/AcceptancePi0Pythia.%s",outputDir.Data(),suffix.Data()));
   histo2DAccPi0Pythia->Draw("copy");
   
   
   
   cout << "Pi0 acceptance in different system - Pythia with added signals" << endl;   
   
   TCanvas* canvasAcceptancePi0PythiaPlus = new TCanvas("canvasAcceptancePi0PythiaPlus","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptancePi0PythiaPlus, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccPi0PythiaPlus;
   histo2DAccPi0PythiaPlus = new TH2F("histo2DAccPi0PythiaPlus","histo2DAccPi0PythiaPlus",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccPi0PythiaPlus, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccPi0PythiaPlus->Draw("copy");

//   DrawGammaSetMarker(histoTrueAccPtPi07TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
//   histoTrueAccPtPi07TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueAccPtPi02760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
//   histoTrueAccPtPi02760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueAccPtPi0900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
//   histoTrueAccPtPi0900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0AcceptPtPythiaPlus8TeV, markerStyleSpectrum8TeVMC, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
   histoPi0AcceptPtPythiaPlus8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0AcceptPtPythiaPlus7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
   histoPi0AcceptPtPythiaPlus7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0AcceptPtPythiaPlus2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
   histoPi0AcceptPtPythiaPlus2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoPi0AcceptPtPythiaPlus900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0AcceptPtPythiaPlus900GeV->DrawCopy("e1,same");    

   TLegend* legendAccPi0PythiaPlus = new TLegend(0.34,0.13,0.93,0.43);
   legendAccPi0PythiaPlus->SetFillColor(0);
   legendAccPi0PythiaPlus->SetLineColor(0);
   legendAccPi0PythiaPlus->SetTextSize(0.027);
   legendAccPi0PythiaPlus->SetNColumns(2);
   legendAccPi0PythiaPlus->AddEntry(histoPi0AcceptPtPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendAccPi0PythiaPlus->AddEntry(histoPi0AcceptPtPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendAccPi0PythiaPlus->AddEntry(histoPi0AcceptPtPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendAccPi0PythiaPlus->AddEntry(histoPi0AcceptPtPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendAccPi0PythiaPlus->AddEntry(histoTrueAccPtPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendAccPi0PythiaPlus->AddEntry(histoTrueAccPtPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendAccPi0PythiaPlus->AddEntry(histoTrueAccPtPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendAccPi0PythiaPlus->Draw();
   
   canvasAcceptancePi0PythiaPlus->SaveAs(Form("%s/AcceptancePi0PythiaPlus.%s",outputDir.Data(),suffix.Data()));
   histo2DAccPi0PythiaPlus->Draw("copy");
   
   
   
   cout << "Pi0 acceptance in different system - Phojet" << endl;   
   
   TCanvas* canvasAcceptancePi0Phojet = new TCanvas("canvasAcceptancePi0Phojet","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptancePi0Phojet, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccPi0Phojet;
   histo2DAccPi0Phojet = new TH2F("histo2DAccPi0Phojet","histo2DAccPi0Phojet",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccPi0Phojet, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccPi0Phojet->Draw("copy");

//   DrawGammaSetMarker(histoTrueAccPtPi07TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
//   histoTrueAccPtPi07TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueAccPtPi02760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
//   histoTrueAccPtPi02760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueAccPtPi0900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
//   histoTrueAccPtPi0900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0AcceptPtPhojet8TeV, markerStyleSpectrum8TeVMC, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
   histoPi0AcceptPtPhojet8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0AcceptPtPhojet7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);  
   histoPi0AcceptPtPhojet7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0AcceptPtPhojet2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
   histoPi0AcceptPtPhojet2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoPi0AcceptPtPhojet900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0AcceptPtPhojet900GeV->DrawCopy("e1,same");    

   TLegend* legendAccPi0Phojet = new TLegend(0.34,0.13,0.93,0.43);
   legendAccPi0Phojet->SetFillColor(0);
   legendAccPi0Phojet->SetLineColor(0);
   legendAccPi0Phojet->SetTextSize(0.027);
   legendAccPi0Phojet->SetNColumns(2);
   legendAccPi0Phojet->AddEntry(histoPi0AcceptPtPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendAccPi0Phojet->AddEntry(histoPi0AcceptPtPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendAccPi0Phojet->AddEntry(histoPi0AcceptPtPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendAccPi0Phojet->AddEntry(histoPi0AcceptPtPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendAccPi0Phojet->AddEntry(histoTrueAccPtPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendAccPi0Phojet->AddEntry(histoTrueAccPtPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendAccPi0Phojet->AddEntry(histoTrueAccPtPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendAccPi0Phojet->Draw();
   
   canvasAcceptancePi0Phojet->SaveAs(Form("%s/AcceptancePi0Phojet.%s",outputDir.Data(),suffix.Data()));
   histo2DAccPi0Phojet->Draw("copy");

   
   //******************************************************************************************************************************************************
   
   cout << "Pi0 acceptance 8TeV" << endl;   
   
   TCanvas* canvasAcceptancePi08TeV = new TCanvas("canvasAcceptancePi08TeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptancePi08TeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccPi08TeV;
   histo2DAccPi08TeV = new TH2F("histo2DAccPi08TeV","histo2DAccPi08TeV",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccPi08TeV, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccPi08TeV->Draw("copy");

   DrawGammaSetMarker(histoPi0AcceptPtPythia8TeV, markerStyleSpectrum8TeVMC, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
   histoPi0AcceptPtPythia8TeV->DrawCopy("e1,same");   
//    DrawGammaSetMarker(histoPi0AcceptPtPythiaPlus8TeV, markerStyleSpectrum8TeVMC+1, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
//    histoPi0AcceptPtPythiaPlus8TeV->DrawCopy("e1,same");
//    DrawGammaSetMarker(histoPi0AcceptPtPhojet8TeV, markerStyleSpectrum8TeVMC+2, markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
//    histoPi0AcceptPtPhojet8TeV->DrawCopy("e1,same");
   
   TLegend* legendAccPi08TeV = new TLegend(0.34,0.13,0.93,0.43);
   legendAccPi08TeV->SetFillColor(0);
   legendAccPi08TeV->SetLineColor(0);
   legendAccPi08TeV->SetTextSize(0.027);
   legendAccPi08TeV->SetNColumns(2);
   legendAccPi08TeV->AddEntry(histoPi0AcceptPtPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia","p");
//    legendAccPi08TeV->AddEntry(histoPi0AcceptPtPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia+AddSign","p");
//    legendAccPi08TeV->AddEntry(histoPi0AcceptPtPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV - Phojet","p");
   legendAccPi08TeV->Draw();
   
   canvasAcceptancePi08TeV->SaveAs(Form("%s/AcceptancePi08TeV.%s",outputDir.Data(),suffix.Data()));
   histo2DAccPi08TeV->Draw("copy");
   
   
   cout << "Pi0 acceptance 7TeV" << endl;   
   
   TCanvas* canvasAcceptancePi07TeV = new TCanvas("canvasAcceptancePi07TeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptancePi07TeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccPi07TeV;
   histo2DAccPi07TeV = new TH2F("histo2DAccPi07TeV","histo2DAccPi07TeV",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccPi07TeV, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccPi07TeV->Draw("copy");

   DrawGammaSetMarker(histoPi0AcceptPtPythia7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV); 
   histoPi0AcceptPtPythia7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0AcceptPtPythiaPlus7TeV, markerStyleSpectrum7TeVMC+1, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV); 
   histoPi0AcceptPtPythiaPlus7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0AcceptPtPhojet7TeV, markerStyleSpectrum7TeVMC+2, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV); 
   histoPi0AcceptPtPhojet7TeV->DrawCopy("e1,same");
   
   TLegend* legendAccPi07TeV = new TLegend(0.34,0.13,0.93,0.43);
   legendAccPi07TeV->SetFillColor(0);
   legendAccPi07TeV->SetLineColor(0);
   legendAccPi07TeV->SetTextSize(0.027);
   legendAccPi07TeV->SetNColumns(2);
   legendAccPi07TeV->AddEntry(histoPi0AcceptPtPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia","p");
   legendAccPi07TeV->AddEntry(histoPi0AcceptPtPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia+AddSign","p");
   legendAccPi07TeV->AddEntry(histoPi0AcceptPtPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV - Phojet","p");
   legendAccPi07TeV->Draw();
   
   canvasAcceptancePi07TeV->SaveAs(Form("%s/AcceptancePi07TeV.%s",outputDir.Data(),suffix.Data()));
   histo2DAccPi07TeV->Draw("copy");
 
   
   cout << "Pi0 acceptance 2760GeV" << endl;   
   
   TCanvas* canvasAcceptancePi02760GeV = new TCanvas("canvasAcceptancePi02760GeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptancePi02760GeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccPi02760GeV;
   histo2DAccPi02760GeV = new TH2F("histo2DAccPi02760GeV","histo2DAccPi02760GeV",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccPi02760GeV, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccPi02760GeV->Draw("copy");

   DrawGammaSetMarker(histoPi0AcceptPtPythia2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV); 
   histoPi0AcceptPtPythia2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0AcceptPtPythiaPlus2760GeV, markerStyleSpectrum2760GeVMC+1, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV); 
   histoPi0AcceptPtPythiaPlus2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0AcceptPtPhojet2760GeV, markerStyleSpectrum2760GeVMC+2, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV); 
   histoPi0AcceptPtPhojet2760GeV->DrawCopy("e1,same");
   
   TLegend* legendAccPi02760GeV = new TLegend(0.34,0.13,0.93,0.43);
   legendAccPi02760GeV->SetFillColor(0);
   legendAccPi02760GeV->SetLineColor(0);
   legendAccPi02760GeV->SetTextSize(0.027);
   legendAccPi02760GeV->SetNColumns(2);
   legendAccPi02760GeV->AddEntry(histoPi0AcceptPtPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia","p");
   legendAccPi02760GeV->AddEntry(histoPi0AcceptPtPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia+AddSign","p");
   legendAccPi02760GeV->AddEntry(histoPi0AcceptPtPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Phojet","p");
   legendAccPi02760GeV->Draw();
   
   canvasAcceptancePi02760GeV->SaveAs(Form("%s/AcceptancePi02760GeV.%s",outputDir.Data(),suffix.Data()));
   histo2DAccPi02760GeV->Draw("copy");
 
   
   cout << "Pi0 acceptance 900GeV" << endl;   
   
   TCanvas* canvasAcceptancePi0900GeV = new TCanvas("canvasAcceptancePi0900GeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptancePi0900GeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccPi0900GeV;
   histo2DAccPi0900GeV = new TH2F("histo2DAccPi0900GeV","histo2DAccPi0900GeV",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccPi0900GeV, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccPi0900GeV->Draw("copy");

   DrawGammaSetMarker(histoPi0AcceptPtPythia900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0AcceptPtPythia900GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0AcceptPtPythiaPlus900GeV, markerStyleSpectrum900GeVMC+1, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0AcceptPtPythiaPlus900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0AcceptPtPhojet900GeV, markerStyleSpectrum900GeVMC+2, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
   histoPi0AcceptPtPhojet900GeV->DrawCopy("e1,same");
   
   TLegend* legendAccPi0900GeV = new TLegend(0.34,0.13,0.93,0.43);
   legendAccPi0900GeV->SetFillColor(0);
   legendAccPi0900GeV->SetLineColor(0);
   legendAccPi0900GeV->SetTextSize(0.027);
   legendAccPi0900GeV->SetNColumns(2);
   legendAccPi0900GeV->AddEntry(histoPi0AcceptPtPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia","p");
   legendAccPi0900GeV->AddEntry(histoPi0AcceptPtPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia+AddSign","p");
   legendAccPi0900GeV->AddEntry(histoPi0AcceptPtPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 Te/*V - Phojet","p");
   legendAccPi0900GeV->Draw();
   
   canvasAcceptancePi0900GeV->SaveAs(Form("%s/AcceptancePi0900GeV.%s",outputDir.Data(),suffix.Data()));
   histo2DAccPi0900GeV->Draw("copy");
 
   
//************************************************************************************************************************************
   
   
   
   cout << "Pi0 raw yield in different system - Pythia" << endl;   

   TCanvas* canvasRawYieldPi0Pythia = new TCanvas("canvasRawYieldPi0Pythia","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldPi0Pythia, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawPi0Pythia;
   canvasRawYieldPi0Pythia->SetLogy();
   histo2DRawPi0Pythia = new TH2F("histo2DRawPi0Pythia","histo2DRawPi0Pythia",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawPi0Pythia, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawPi0Pythia->Draw("copy");
   
//   DrawGammaSetMarker(histoRawYieldPi07TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
//   histoRawYieldPi07TeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldPi02760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
//   histoRawYieldPi02760GeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldPi0900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
//   histoRawYieldPi0900GeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0RawYieldPythia8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
   histoPi0RawYieldPythia8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0RawYieldPythia7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
   histoPi0RawYieldPythia7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0RawYieldPythia2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
   histoPi0RawYieldPythia2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0RawYieldPythia900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
   histoPi0RawYieldPythia900GeV->DrawCopy("e1,same");     

   TLatex *labelRawPi0Pythia = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawPi0Pythia, 0.038,4);
   labelRawPi0Pythia->Draw();

   TLegend* legendRawYieldPi0Pythia = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldPi0Pythia->SetFillColor(0);
   legendRawYieldPi0Pythia->SetLineColor(0);
   legendRawYieldPi0Pythia->SetTextSize(0.04);
//   legendRawYieldPi0Pythia->AddEntry(histoRawYieldPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
//   legendRawYieldPi0Pythia->AddEntry(histoRawYieldPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//   legendRawYieldPi0Pythia->AddEntry(histoRawYieldPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 
   legendRawYieldPi0Pythia->AddEntry(histoPi0RawYieldPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendRawYieldPi0Pythia->AddEntry(histoPi0RawYieldPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendRawYieldPi0Pythia->AddEntry(histoPi0RawYieldPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendRawYieldPi0Pythia->AddEntry(histoPi0RawYieldPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendRawYieldPi0Pythia->Draw();
   
   canvasRawYieldPi0Pythia->SaveAs(Form("%s/RawYieldCompPi0Pythia.%s",outputDir.Data(),suffix.Data()));


   cout << "Pi0 raw yield in different system - PythiaPlus" << endl;   

   TCanvas* canvasRawYieldPi0PythiaPlus = new TCanvas("canvasRawYieldPi0PythiaPlus","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldPi0PythiaPlus, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawPi0PythiaPlus;
   canvasRawYieldPi0PythiaPlus->SetLogy();
   histo2DRawPi0PythiaPlus = new TH2F("histo2DRawPi0PythiaPlus","histo2DRawPi0PythiaPlus",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawPi0PythiaPlus, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawPi0PythiaPlus->Draw("copy");
   
//   DrawGammaSetMarker(histoRawYieldPi07TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
//   histoRawYieldPi07TeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldPi02760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
//   histoRawYieldPi02760GeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldPi0900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
//   histoRawYieldPi0900GeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0RawYieldPythiaPlus8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
   histoPi0RawYieldPythiaPlus8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0RawYieldPythiaPlus7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
   histoPi0RawYieldPythiaPlus7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0RawYieldPythiaPlus2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
   histoPi0RawYieldPythiaPlus2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0RawYieldPythiaPlus900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
   histoPi0RawYieldPythiaPlus900GeV->DrawCopy("e1,same");     

   TLatex *labelRawPi0PythiaPlus = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawPi0PythiaPlus, 0.038,4);
   labelRawPi0PythiaPlus->Draw();

   TLegend* legendRawYieldPi0PythiaPlus = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldPi0PythiaPlus->SetFillColor(0);
   legendRawYieldPi0PythiaPlus->SetLineColor(0);
   legendRawYieldPi0PythiaPlus->SetTextSize(0.04);
//   legendRawYieldPi0PythiaPlus->AddEntry(histoRawYieldPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
//   legendRawYieldPi0PythiaPlus->AddEntry(histoRawYieldPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//   legendRawYieldPi0PythiaPlus->AddEntry(histoRawYieldPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 
   legendRawYieldPi0PythiaPlus->AddEntry(histoPi0RawYieldPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendRawYieldPi0PythiaPlus->AddEntry(histoPi0RawYieldPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendRawYieldPi0PythiaPlus->AddEntry(histoPi0RawYieldPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendRawYieldPi0PythiaPlus->AddEntry(histoPi0RawYieldPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendRawYieldPi0PythiaPlus->Draw();
   
   canvasRawYieldPi0PythiaPlus->SaveAs(Form("%s/RawYieldCompPi0PythiaPlus.%s",outputDir.Data(),suffix.Data()));


   
   cout << "Pi0 raw yield in different system - Phojet" << endl;   

   TCanvas* canvasRawYieldPi0Phojet = new TCanvas("canvasRawYieldPi0Phojet","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldPi0Phojet, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawPi0Phojet;
   canvasRawYieldPi0Phojet->SetLogy();
   histo2DRawPi0Phojet = new TH2F("histo2DRawPi0Phojet","histo2DRawPi0Phojet",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawPi0Phojet, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawPi0Phojet->Draw("copy");
   
//   DrawGammaSetMarker(histoRawYieldPi07TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
//   histoRawYieldPi07TeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldPi02760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
//   histoRawYieldPi02760GeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldPi0900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
//   histoRawYieldPi0900GeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0RawYieldPhojet8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
   histoPi0RawYieldPhojet8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0RawYieldPhojet7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
   histoPi0RawYieldPhojet7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0RawYieldPhojet2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
   histoPi0RawYieldPhojet2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0RawYieldPhojet900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
   histoPi0RawYieldPhojet900GeV->DrawCopy("e1,same");     

   TLatex *labelRawPi0Phojet = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawPi0Phojet, 0.038,4);
   labelRawPi0Phojet->Draw();

   TLegend* legendRawYieldPi0Phojet = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldPi0Phojet->SetFillColor(0);
   legendRawYieldPi0Phojet->SetLineColor(0);
   legendRawYieldPi0Phojet->SetTextSize(0.04);
//   legendRawYieldPi0Phojet->AddEntry(histoRawYieldPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
//   legendRawYieldPi0Phojet->AddEntry(histoRawYieldPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//   legendRawYieldPi0Phojet->AddEntry(histoRawYieldPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 
   legendRawYieldPi0Phojet->AddEntry(histoPi0RawYieldPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendRawYieldPi0Phojet->AddEntry(histoPi0RawYieldPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendRawYieldPi0Phojet->AddEntry(histoPi0RawYieldPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendRawYieldPi0Phojet->AddEntry(histoPi0RawYieldPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendRawYieldPi0Phojet->Draw();
   
   canvasRawYieldPi0Phojet->SaveAs(Form("%s/RawYieldCompPi0Phojet.%s",outputDir.Data(),suffix.Data()));

   
//*************************************************************************************************************************
   

   cout << "Pi0 raw yield 8TeV " << endl;   

   TCanvas* canvasRawYieldPi08TeV = new TCanvas("canvasRawYieldPi08TeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldPi08TeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawPi08TeV;
   canvasRawYieldPi08TeV->SetLogy();
   histo2DRawPi08TeV = new TH2F("histo2DRawPi08TeV","histo2DRawPi08TeV",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawPi08TeV, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawPi08TeV->Draw("copy");
   
   DrawGammaSetMarker(histoPi0RawYieldPythia8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
   histoPi0RawYieldPythia8TeV->DrawCopy("e1,same");
//    DrawGammaSetMarker(histoPi0RawYieldPythiaPlus8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
//    histoPi0RawYieldPythiaPlus8TeV->DrawCopy("e1,same");
//    DrawGammaSetMarker(histoPi0RawYieldPhojet8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
//    histoPi0RawYieldPhojet8TeV->DrawCopy("e1,same");
   
   TLatex *labelRawPi08TeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawPi08TeV, 0.038,4);
   labelRawPi08TeV->Draw();

   TLegend* legendRawYieldPi08TeV = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldPi08TeV->SetFillColor(0);
   legendRawYieldPi08TeV->SetLineColor(0);
   legendRawYieldPi08TeV->SetTextSize(0.04);
   legendRawYieldPi08TeV->AddEntry(histoPi0RawYieldPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia","p");      
//    legendRawYieldPi08TeV->AddEntry(histoPi0RawYieldPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia+AddSign","p");      
//    legendRawYieldPi08TeV->AddEntry(histoPi0RawYieldPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV - Phojet","p");      
   legendRawYieldPi08TeV->Draw();
   
   canvasRawYieldPi08TeV->SaveAs(Form("%s/RawYieldCompPi08TeV.%s",outputDir.Data(),suffix.Data()));

  
   cout << "Pi0 raw yield 7TeV " << endl;   

   TCanvas* canvasRawYieldPi07TeV = new TCanvas("canvasRawYieldPi07TeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldPi07TeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawPi07TeV;
   canvasRawYieldPi07TeV->SetLogy();
   histo2DRawPi07TeV = new TH2F("histo2DRawPi07TeV","histo2DRawPi07TeV",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawPi07TeV, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawPi07TeV->Draw("copy");
   
   DrawGammaSetMarker(histoPi0RawYieldPythia7TeV , markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);  
   histoPi0RawYieldPythia7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0RawYieldPythiaPlus7TeV , markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);  
   histoPi0RawYieldPythiaPlus7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0RawYieldPhojet7TeV , markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);  
   histoPi0RawYieldPhojet7TeV->DrawCopy("e1,same");
   
   TLatex *labelRawPi07TeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawPi07TeV, 0.038,4);
   labelRawPi07TeV->Draw();

   TLegend* legendRawYieldPi07TeV = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldPi07TeV->SetFillColor(0);
   legendRawYieldPi07TeV->SetLineColor(0);
   legendRawYieldPi07TeV->SetTextSize(0.04);
   legendRawYieldPi07TeV->AddEntry(histoPi0RawYieldPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia","p");      
   legendRawYieldPi07TeV->AddEntry(histoPi0RawYieldPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia+AddSign","p");      
   legendRawYieldPi07TeV->AddEntry(histoPi0RawYieldPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV - Phojet","p");      
   legendRawYieldPi07TeV->Draw();
   
   canvasRawYieldPi07TeV->SaveAs(Form("%s/RawYieldCompPi07TeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Pi0 raw yield 2760GeV " << endl;   

   TCanvas* canvasRawYieldPi02760GeV = new TCanvas("canvasRawYieldPi02760GeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldPi02760GeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawPi02760GeV;
   canvasRawYieldPi02760GeV->SetLogy();
   histo2DRawPi02760GeV = new TH2F("histo2DRawPi02760GeV","histo2DRawPi02760GeV",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawPi02760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawPi02760GeV->Draw("copy");
   
   DrawGammaSetMarker(histoPi0RawYieldPythia2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);  
   histoPi0RawYieldPythia2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0RawYieldPythiaPlus2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);  
   histoPi0RawYieldPythiaPlus2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0RawYieldPhojet2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);  
   histoPi0RawYieldPhojet2760GeV->DrawCopy("e1,same");
   
   TLatex *labelRawPi02760GeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawPi02760GeV, 0.038,4);
   labelRawPi02760GeV->Draw();

   TLegend* legendRawYieldPi02760GeV = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldPi02760GeV->SetFillColor(0);
   legendRawYieldPi02760GeV->SetLineColor(0);
   legendRawYieldPi02760GeV->SetTextSize(0.04);
   legendRawYieldPi02760GeV->AddEntry(histoPi0RawYieldPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia","p");      
   legendRawYieldPi02760GeV->AddEntry(histoPi0RawYieldPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia+AddSign","p");      
   legendRawYieldPi02760GeV->AddEntry(histoPi0RawYieldPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Phojet","p");      
   legendRawYieldPi02760GeV->Draw();
   
   canvasRawYieldPi02760GeV->SaveAs(Form("%s/RawYieldCompPi02760GeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Pi0 raw yield 900GeV " << endl;   

   TCanvas* canvasRawYieldPi0900GeV = new TCanvas("canvasRawYieldPi0900GeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldPi0900GeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawPi0900GeV;
   canvasRawYieldPi0900GeV->SetLogy();
   histo2DRawPi0900GeV = new TH2F("histo2DRawPi0900GeV","histo2DRawPi0900GeV",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawPi0900GeV->Draw("copy");
   
   DrawGammaSetMarker(histoPi0RawYieldPythia900GeV , markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);  
   histoPi0RawYieldPythia900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0RawYieldPythiaPlus900GeV , markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);  
   histoPi0RawYieldPythiaPlus900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0RawYieldPhojet900GeV , markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);  
   histoPi0RawYieldPhojet900GeV->DrawCopy("e1,same");
   
   TLatex *labelRawPi0900GeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawPi0900GeV, 0.038,4);
   labelRawPi0900GeV->Draw();

   TLegend* legendRawYieldPi0900GeV = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldPi0900GeV->SetFillColor(0);
   legendRawYieldPi0900GeV->SetLineColor(0);
   legendRawYieldPi0900GeV->SetTextSize(0.04);
   legendRawYieldPi0900GeV->AddEntry(histoPi0RawYieldPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia","p");      
   legendRawYieldPi0900GeV->AddEntry(histoPi0RawYieldPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia+AddSign","p");      
   legendRawYieldPi0900GeV->AddEntry(histoPi0RawYieldPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Phojet","p");      
   legendRawYieldPi0900GeV->Draw();
   
   canvasRawYieldPi0900GeV->SaveAs(Form("%s/RawYieldCompPi0900GeV.%s",outputDir.Data(),suffix.Data()));

   
//////////////////////////////// spectra all together //////////////////////////////////////////////////////////////////////////

//************************************************************************************************************************************
   
   
   
   cout << "Pi0 corrected  yield in different system - Pythia" << endl;   

   TCanvas* canvasCorrectedSpecPi0Pythia = new TCanvas("canvasCorrectedSpecPi0Pythia","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecPi0Pythia, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrPi0Pythia;
   canvasCorrectedSpecPi0Pythia->SetLogy();
   histo2DCorrPi0Pythia = new TH2F("histo2DCorrPi0Pythia","histo2DCorrPi0Pythia",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrPi0Pythia, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrPi0Pythia->Draw("copy");
   
   DrawGammaSetMarker(histoPi0CorrectedSpecPythia8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
   histoPi0CorrectedSpecPythia8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0CorrectedSpecPythia7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
   histoPi0CorrectedSpecPythia7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0CorrectedSpecPythia2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
   histoPi0CorrectedSpecPythia2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0CorrectedSpecPythia900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
   histoPi0CorrectedSpecPythia900GeV->DrawCopy("e1,same");     

   TLatex *labelCorrPi0Pythia = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrPi0Pythia, 0.038,4);
   labelCorrPi0Pythia->Draw();

   TLegend* legendCorrectedSpecPi0Pythia = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecPi0Pythia->SetFillColor(0);
   legendCorrectedSpecPi0Pythia->SetLineColor(0);
   legendCorrectedSpecPi0Pythia->SetTextSize(0.04); 
   legendCorrectedSpecPi0Pythia->AddEntry(histoPi0CorrectedSpecPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendCorrectedSpecPi0Pythia->AddEntry(histoPi0CorrectedSpecPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendCorrectedSpecPi0Pythia->AddEntry(histoPi0CorrectedSpecPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendCorrectedSpecPi0Pythia->AddEntry(histoPi0CorrectedSpecPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendCorrectedSpecPi0Pythia->Draw();
   
   canvasCorrectedSpecPi0Pythia->SaveAs(Form("%s/CorrectedSpecCompPi0Pythia.%s",outputDir.Data(),suffix.Data()));


   cout << "Pi0 corrected  yield in different system - PythiaPlus" << endl;   

   TCanvas* canvasCorrectedSpecPi0PythiaPlus = new TCanvas("canvasCorrectedSpecPi0PythiaPlus","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecPi0PythiaPlus, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrPi0PythiaPlus;
   canvasCorrectedSpecPi0PythiaPlus->SetLogy();
   histo2DCorrPi0PythiaPlus = new TH2F("histo2DCorrPi0PythiaPlus","histo2DCorrPi0PythiaPlus",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrPi0PythiaPlus, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrPi0PythiaPlus->Draw("copy");
   
//   DrawGammaSetMarker(histoCorrectedSpecPi07TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
//   histoCorrectedSpecPi07TeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoCorrectedSpecPi02760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
//   histoCorrectedSpecPi02760GeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoCorrectedSpecPi0900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
//   histoCorrectedSpecPi0900GeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0CorrectedSpecPythiaPlus8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
   histoPi0CorrectedSpecPythiaPlus8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0CorrectedSpecPythiaPlus7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
   histoPi0CorrectedSpecPythiaPlus7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0CorrectedSpecPythiaPlus2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
   histoPi0CorrectedSpecPythiaPlus2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0CorrectedSpecPythiaPlus900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
   histoPi0CorrectedSpecPythiaPlus900GeV->DrawCopy("e1,same");     

   TLatex *labelCorrPi0PythiaPlus = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrPi0PythiaPlus, 0.038,4);
   labelCorrPi0PythiaPlus->Draw();

   TLegend* legendCorrectedSpecPi0PythiaPlus = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecPi0PythiaPlus->SetFillColor(0);
   legendCorrectedSpecPi0PythiaPlus->SetLineColor(0);
   legendCorrectedSpecPi0PythiaPlus->SetTextSize(0.04);
//   legendCorrectedSpecPi0PythiaPlus->AddEntry(histoCorrectedSpecPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
//   legendCorrectedSpecPi0PythiaPlus->AddEntry(histoCorrectedSpecPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//   legendCorrectedSpecPi0PythiaPlus->AddEntry(histoCorrectedSpecPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 
   legendCorrectedSpecPi0PythiaPlus->AddEntry(histoPi0CorrectedSpecPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendCorrectedSpecPi0PythiaPlus->AddEntry(histoPi0CorrectedSpecPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendCorrectedSpecPi0PythiaPlus->AddEntry(histoPi0CorrectedSpecPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendCorrectedSpecPi0PythiaPlus->AddEntry(histoPi0CorrectedSpecPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendCorrectedSpecPi0PythiaPlus->Draw();
   
   canvasCorrectedSpecPi0PythiaPlus->SaveAs(Form("%s/CorrectedSpecCompPi0PythiaPlus.%s",outputDir.Data(),suffix.Data()));


   
   cout << "Pi0 corrected yield in different system - Phojet" << endl;   

   TCanvas* canvasCorrectedSpecPi0Phojet = new TCanvas("canvasCorrectedSpecPi0Phojet","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecPi0Phojet, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrPi0Phojet;
   canvasCorrectedSpecPi0Phojet->SetLogy();
   histo2DCorrPi0Phojet = new TH2F("histo2DCorrPi0Phojet","histo2DCorrPi0Phojet",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrPi0Phojet, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrPi0Phojet->Draw("copy");
   
//   DrawGammaSetMarker(histoCorrectedSpecPi07TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
//   histoCorrectedSpecPi07TeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoCorrectedSpecPi02760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
//   histoCorrectedSpecPi02760GeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoCorrectedSpecPi0900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
//   histoCorrectedSpecPi0900GeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoPi0CorrectedSpecPhojet8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
   histoPi0CorrectedSpecPhojet8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0CorrectedSpecPhojet7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV); 
   histoPi0CorrectedSpecPhojet7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0CorrectedSpecPhojet2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV); 
   histoPi0CorrectedSpecPhojet2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoPi0CorrectedSpecPhojet900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);   
   histoPi0CorrectedSpecPhojet900GeV->DrawCopy("e1,same");     

   TLatex *labelCorrPi0Phojet = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrPi0Phojet, 0.038,4);
   labelCorrPi0Phojet->Draw();

   TLegend* legendCorrectedSpecPi0Phojet = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecPi0Phojet->SetFillColor(0);
   legendCorrectedSpecPi0Phojet->SetLineColor(0);
   legendCorrectedSpecPi0Phojet->SetTextSize(0.04);
//   legendCorrectedSpecPi0Phojet->AddEntry(histoCorrectedSpecPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
//   legendCorrectedSpecPi0Phojet->AddEntry(histoCorrectedSpecPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//   legendCorrectedSpecPi0Phojet->AddEntry(histoCorrectedSpecPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 
   legendCorrectedSpecPi0Phojet->AddEntry(histoPi0CorrectedSpecPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendCorrectedSpecPi0Phojet->AddEntry(histoPi0CorrectedSpecPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendCorrectedSpecPi0Phojet->AddEntry(histoPi0CorrectedSpecPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendCorrectedSpecPi0Phojet->AddEntry(histoPi0CorrectedSpecPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendCorrectedSpecPi0Phojet->Draw();
   
   canvasCorrectedSpecPi0Phojet->SaveAs(Form("%s/CorrectedSpecCompPi0Phojet.%s",outputDir.Data(),suffix.Data()));

   
//*************************************************************************************************************************
   

   cout << "Pi0 corrected yield 8TeV " << endl;   

   TCanvas* canvasCorrectedSpecPi08TeV = new TCanvas("canvasCorrectedSpecPi08TeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecPi08TeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrPi08TeV;
   canvasCorrectedSpecPi08TeV->SetLogy();
   histo2DCorrPi08TeV = new TH2F("histo2DCorrPi08TeV","histo2DCorrPi08TeV",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrPi08TeV, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrPi08TeV->Draw("copy");
   
   DrawGammaSetMarker(histoPi0CorrectedSpecPythia8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
   histoPi0CorrectedSpecPythia8TeV->DrawCopy("e1,same");
//    DrawGammaSetMarker(histoPi0CorrectedSpecPythiaPlus8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
//    histoPi0CorrectedSpecPythiaPlus8TeV->DrawCopy("e1,same");
//    DrawGammaSetMarker(histoPi0CorrectedSpecPhojet8TeV , markerStyleSpectrum8TeV, markerSizePi0PP8TeV*1.5, colorPi08TeV, colorPi08TeV);  
//    histoPi0CorrectedSpecPhojet8TeV->DrawCopy("e1,same");
   
   TLatex *labelCorrPi08TeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrPi08TeV, 0.038,4);
   labelCorrPi08TeV->Draw();

   TLegend* legendCorrectedSpecPi08TeV = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecPi08TeV->SetFillColor(0);
   legendCorrectedSpecPi08TeV->SetLineColor(0);
   legendCorrectedSpecPi08TeV->SetTextSize(0.04);
   legendCorrectedSpecPi08TeV->AddEntry(histoPi0CorrectedSpecPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia","p");      
//    legendCorrectedSpecPi08TeV->AddEntry(histoPi0CorrectedSpecPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia+AddSign","p");      
//    legendCorrectedSpecPi08TeV->AddEntry(histoPi0CorrectedSpecPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV - Phojet","p");      
   legendCorrectedSpecPi08TeV->Draw();
   
   canvasCorrectedSpecPi08TeV->SaveAs(Form("%s/CorrectedSpecCompPi08TeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Pi0 corrected yield 7TeV " << endl;   

   TCanvas* canvasCorrectedSpecPi07TeV = new TCanvas("canvasCorrectedSpecPi07TeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecPi07TeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrPi07TeV;
   canvasCorrectedSpecPi07TeV->SetLogy();
   histo2DCorrPi07TeV = new TH2F("histo2DCorrPi07TeV","histo2DCorrPi07TeV",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrPi07TeV, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrPi07TeV->Draw("copy");
   
   DrawGammaSetMarker(histoPi0CorrectedSpecPythia7TeV , markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);  
   histoPi0CorrectedSpecPythia7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0CorrectedSpecPythiaPlus7TeV , markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);  
   histoPi0CorrectedSpecPythiaPlus7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0CorrectedSpecPhojet7TeV , markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);  
   histoPi0CorrectedSpecPhojet7TeV->DrawCopy("e1,same");
   
   TLatex *labelCorrPi07TeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrPi07TeV, 0.038,4);
   labelCorrPi07TeV->Draw();

   TLegend* legendCorrectedSpecPi07TeV = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecPi07TeV->SetFillColor(0);
   legendCorrectedSpecPi07TeV->SetLineColor(0);
   legendCorrectedSpecPi07TeV->SetTextSize(0.04);
   legendCorrectedSpecPi07TeV->AddEntry(histoPi0CorrectedSpecPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia","p");      
   legendCorrectedSpecPi07TeV->AddEntry(histoPi0CorrectedSpecPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia+AddSign","p");      
   legendCorrectedSpecPi07TeV->AddEntry(histoPi0CorrectedSpecPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV - Phojet","p");      
   legendCorrectedSpecPi07TeV->Draw();
   
   canvasCorrectedSpecPi07TeV->SaveAs(Form("%s/CorrectedSpecCompPi07TeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Pi0 corrected yield 2760GeV " << endl;   

   TCanvas* canvasCorrectedSpecPi02760GeV = new TCanvas("canvasCorrectedSpecPi02760GeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecPi02760GeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrPi02760GeV;
   canvasCorrectedSpecPi02760GeV->SetLogy();
   histo2DCorrPi02760GeV = new TH2F("histo2DCorrPi02760GeV","histo2DCorrPi02760GeV",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrPi02760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrPi02760GeV->Draw("copy");
   
   DrawGammaSetMarker(histoPi0CorrectedSpecPythia2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);  
   histoPi0CorrectedSpecPythia2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0CorrectedSpecPythiaPlus2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);  
   histoPi0CorrectedSpecPythiaPlus2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0CorrectedSpecPhojet2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);  
   histoPi0CorrectedSpecPhojet2760GeV->DrawCopy("e1,same");
   
   TLatex *labelCorrPi02760GeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrPi02760GeV, 0.038,4);
   labelCorrPi02760GeV->Draw();

   TLegend* legendCorrectedSpecPi02760GeV = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecPi02760GeV->SetFillColor(0);
   legendCorrectedSpecPi02760GeV->SetLineColor(0);
   legendCorrectedSpecPi02760GeV->SetTextSize(0.04);
   legendCorrectedSpecPi02760GeV->AddEntry(histoPi0CorrectedSpecPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia","p");      
   legendCorrectedSpecPi02760GeV->AddEntry(histoPi0CorrectedSpecPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia+AddSign","p");      
   legendCorrectedSpecPi02760GeV->AddEntry(histoPi0CorrectedSpecPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Phojet","p");      
   legendCorrectedSpecPi02760GeV->Draw();
   
   canvasCorrectedSpecPi02760GeV->SaveAs(Form("%s/CorrectedSpecCompPi02760GeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Pi0 corrected yield 900GeV " << endl;   

   TCanvas* canvasCorrectedSpecPi0900GeV = new TCanvas("canvasCorrectedSpecPi0900GeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecPi0900GeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrPi0900GeV;
   canvasCorrectedSpecPi0900GeV->SetLogy();
   histo2DCorrPi0900GeV = new TH2F("histo2DCorrPi0900GeV","histo2DCorrPi0900GeV",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrPi0900GeV->Draw("copy");
   
   DrawGammaSetMarker(histoPi0CorrectedSpecPythia900GeV , markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);  
   histoPi0CorrectedSpecPythia900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0CorrectedSpecPythiaPlus900GeV , markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);  
   histoPi0CorrectedSpecPythiaPlus900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoPi0CorrectedSpecPhojet900GeV , markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);  
   histoPi0CorrectedSpecPhojet900GeV->DrawCopy("e1,same");
   
   TLatex *labelCorrPi0900GeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrPi0900GeV, 0.038,4);
   labelCorrPi0900GeV->Draw();

   TLegend* legendCorrectedSpecPi0900GeV = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecPi0900GeV->SetFillColor(0);
   legendCorrectedSpecPi0900GeV->SetLineColor(0);
   legendCorrectedSpecPi0900GeV->SetTextSize(0.04);
   legendCorrectedSpecPi0900GeV->AddEntry(histoPi0CorrectedSpecPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia","p");      
   legendCorrectedSpecPi0900GeV->AddEntry(histoPi0CorrectedSpecPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia+AddSign","p");      
   legendCorrectedSpecPi0900GeV->AddEntry(histoPi0CorrectedSpecPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Phojet","p");      
   legendCorrectedSpecPi0900GeV->Draw();
   
   canvasCorrectedSpecPi0900GeV->SaveAs(Form("%s/CorrectedSpecCompPi0900GeV.%s",outputDir.Data(),suffix.Data()));

   
   


 	TCanvas* canvasFraction8TeV = new TCanvas("canvasFraction8TeV","",1550,1200);  // gives the page size
 	canvasFraction8TeV->SetTickx();
 	canvasFraction8TeV->SetTicky();
 	canvasFraction8TeV->SetGridx(0);
 	canvasFraction8TeV->SetGridy(0);
 	canvasFraction8TeV->SetLogy(0);
 	canvasFraction8TeV->SetLeftMargin(0.13);
 	canvasFraction8TeV->SetRightMargin(0.02);
 	canvasFraction8TeV->SetTopMargin(0.02);
 	canvasFraction8TeV->SetFillColor(0);
   
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation PP ***********
   //**********************************************************************************


   canvasCorrectedSpecPi08TeV->cd();
   canvasCorrectedSpecPi08TeV->SetLogx(1);
   histo2DCorrPi08TeV->GetXaxis()->SetRangeUser(0.01,30);
   histo2DCorrPi08TeV->DrawCopy(); 
   

   TString  nameFinalResDatPi08TeV = Form("%s/FitResultsMCPi08TeV.dat",outputDir.Data());
   TString forOutputPi08TeV;
   fstream fileFinalResultsPi08TeV;
   fileFinalResultsPi08TeV.open(nameFinalResDatPi08TeV.Data(), ios::out);

   TF1* fitYieldDataQCDPi08TeV = NULL;
   if (directoryPythiaPi08TeV){
      canvasCorrectedSpecPi08TeV->cd();
      canvasCorrectedSpecPi08TeV->SetLogy(1);
      canvasCorrectedSpecPi08TeV->SetLogx(1);
      histo2DCorrPi08TeV->GetXaxis()->SetRangeUser(0.01,30);
      histo2DCorrPi08TeV->DrawCopy(); 
   
      histoMCYieldPi0Pythia->SetMarkerStyle(markerStyleMCPP8TeV);
      histoMCYieldPi0Pythia->Draw("hist,pe1,same");
   
      histoPi0CorrectedSpecPythia8TeV->Draw("hist,pe1,same");

      fitYieldDataQCDPi08TeV = FitObject("qcd","fitYieldDataQCDPi08TeV","Pi0",histoPi0CorrectedSpecPythia8TeV,0.4,14,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDPi08TeV, 1, 1.5, colorPi08TeV);
      fitYieldDataQCDPi08TeV->Draw("same");

      forOutputPi08TeV = WriteParameterToFile(fitYieldDataQCDPi08TeV);

      fileFinalResultsPi08TeV << forOutputPi08TeV.Data() << endl;  

      canvasCorrectedSpecPi08TeV->Update();
      canvasCorrectedSpecPi08TeV->Print(Form("%s/Pi0_MCInputSpectraFittedPP.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDPi08TeV = CalculateHistoRatioToFit(histoPi0CorrectedSpecPythia8TeV, fitYieldDataQCDPi08TeV);
      TH1D* histoRatioMCtoDataFitQCDPi08TeV = CalculateHistoRatioToFit(histoMCYieldPi0Pythia, fitYieldDataQCDPi08TeV);
      TH1D* histoRatioMCUnweightedtoDataFitQCDPi08TeV = NULL;

      if (histoMCPi0YieldPtPythia8TeVWOWeights) histoRatioMCUnweightedtoDataFitQCDPi08TeV = CalculateHistoRatioToFit(histoMCPi0YieldPtPythia8TeVWOWeights, fitYieldDataQCDPi08TeV);
      canvasFraction8TeV->cd();
      if (histoRatioMCUnweightedtoDataFitQCDPi08TeV) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPi08TeV, markerStyleMCPP8TeV,markerSizePi0PP8TeV, colorPi08TeV, colorPi08TeV); 
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDPi08TeV,  markerStyleMCPP8TeV+1,markerSizePi0PP8TeV, colorPi08TeV+2, colorPi08TeV+2); 
      DrawGammaSetMarker(histoRatioDatatoFitQCDPi08TeV, markerStyleMCPP8TeV,markerSizePi0PP8TeV, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPi08TeV,
                  "", "#it{p}_{T} (GeV/#it{c})", "Pi0 Spectrum/ fit to Spectrum",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 13.5);
      histoRatioDatatoFitQCDPi08TeV->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDPi08TeV->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDPi08TeV) histoRatioMCUnweightedtoDataFitQCDPi08TeV->Draw("same,e,p");  
      TLegend* legendFit = new TLegend(0.16,0.81,0.4,0.9);
      legendFit->SetFillColor(0);
      legendFit->SetLineColor(0);
      legendFit->SetTextSize(0.025);
      legendFit->SetMargin(0.2);
      legendFit->AddEntry(histoRatioDatatoFitQCDPi08TeV,"Data/QCD fit to Data","p");
      if (runDrawReweighted) legendFit->AddEntry(histoRatioMCtoDataFitQCDPi08TeV,"MC weighted/QCD fit to Data","p");
      if (histoRatioMCUnweightedtoDataFitQCDPi08TeV) legendFit->AddEntry(histoRatioMCUnweightedtoDataFitQCDPi08TeV,"MC/QCD fit to Data","p");
      legendFit->Draw();
      TLatex *labelRatioMCData = new TLatex(0.2,0.92,collisionSystemPP8TeV.Data());
      SetStyleTLatex(labelRatioMCData, 0.04,4);
      labelRatioMCData->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFraction8TeV->Update();
      canvasFraction8TeV->SaveAs(Form("%s/Pi0_Ratio_pp_MCToDataFit_8TeV.%s",outputDir.Data(),suffix.Data()));
   }
   

 	TCanvas* canvasFraction7TeV = new TCanvas("canvasFraction7TeV","",1550,1200);  // gives the page size
 	canvasFraction7TeV->SetTickx();
 	canvasFraction7TeV->SetTicky();
 	canvasFraction7TeV->SetGridx(0);
 	canvasFraction7TeV->SetGridy(0);
 	canvasFraction7TeV->SetLogy(0);
 	canvasFraction7TeV->SetLeftMargin(0.13);
 	canvasFraction7TeV->SetRightMargin(0.02);
 	canvasFraction7TeV->SetTopMargin(0.02);
 	canvasFraction7TeV->SetFillColor(0);
   
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation PP ***********
   //**********************************************************************************


   canvasCorrectedSpecPi07TeV->cd();
   canvasCorrectedSpecPi07TeV->SetLogx(1);
   histo2DCorrPi07TeV->GetXaxis()->SetRangeUser(0.01,30);
   histo2DCorrPi07TeV->DrawCopy(); 
   

   TString  nameFinalResDatPi07TeV = Form("%s/FitResultsMCPi07TeV.dat",outputDir.Data());
   TString forOutputPi07TeV;
   fstream fileFinalResultsPi07TeV;
   fileFinalResultsPi07TeV.open(nameFinalResDatPi07TeV.Data(), ios::out);

   TF1* fitYieldDataQCDPi07TeV = NULL;
   if (directoryPythiaPi07TeV){
      canvasCorrectedSpecPi07TeV->cd();
      canvasCorrectedSpecPi07TeV->SetLogy(1);
      canvasCorrectedSpecPi07TeV->SetLogx(1);
      histo2DCorrPi07TeV->GetXaxis()->SetRangeUser(0.01,30);
      histo2DCorrPi07TeV->DrawCopy(); 
   
      histoMCYieldPi0Pythia->SetMarkerStyle(markerStyleMCPP7TeV);
      histoMCYieldPi0Pythia->Draw("hist,pe1,same");
   
      histoPi0CorrectedSpecPythia7TeV->Draw("hist,pe1,same");

      fitYieldDataQCDPi07TeV = FitObject("qcd","fitYieldDataQCDPi07TeV","Pi0",histoPi0CorrectedSpecPythia7TeV,0.4,14,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDPi07TeV, 1, 1.5, colorPi07TeV);
      fitYieldDataQCDPi07TeV->Draw("same");

      forOutputPi07TeV = WriteParameterToFile(fitYieldDataQCDPi07TeV);

      fileFinalResultsPi07TeV << forOutputPi07TeV.Data() << endl;  

      canvasCorrectedSpecPi07TeV->Update();
      canvasCorrectedSpecPi07TeV->Print(Form("%s/Pi0_MCInputSpectraFittedPP.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDPi07TeV = CalculateHistoRatioToFit(histoPi0CorrectedSpecPythia7TeV, fitYieldDataQCDPi0Pi07TeV);
      TH1D* histoRatioMCtoDataFitQCDPi07TeV = CalculateHistoRatioToFit(histoMCYieldPi0Pythia, fitYieldDataQCDPi0Pi07TeV);
      TH1D* histoRatioMCUnweightedtoDataFitQCDPi07TeV = NULL;

      if (histoMCPi0YieldPtPythia7TeVWOWeights) histoRatioMCUnweightedtoDataFitQCDPi07TeV = CalculateHistoRatioToFit(histoMCPi0YieldPtPythia7TeVWOWeights, fitYieldDataQCDPi0Pi07TeV);
      canvasFraction7TeV->cd();
      if (histoRatioMCUnweightedtoDataFitQCDPi07TeV) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPi07TeV, markerStyleMCPP7TeV,markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV); 
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDPi07TeV,  markerStyleMCPP7TeV+1,markerSizePi0PP7TeV, colorPi07TeV+2, colorPi07TeV+2); 
      DrawGammaSetMarker(histoRatioDatatoFitQCDPi07TeV, markerStyleMCPP7TeV,markerSizePi0PP7TeV, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPi07TeV,
                  "", "#it{p}_{T} (GeV/#it{c})", "Pi0 Spectrum/ fit to Spectrum",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 13.5);
      histoRatioDatatoFitQCDPi07TeV->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDPi07TeV->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDPi07TeV) histoRatioMCUnweightedtoDataFitQCDPi07TeV->Draw("same,e,p");  
      TLegend* legendFit = new TLegend(0.16,0.81,0.4,0.9);
      legendFit->SetFillColor(0);
      legendFit->SetLineColor(0);
      legendFit->SetTextSize(0.025);
      legendFit->SetMargin(0.2);
      legendFit->AddEntry(histoRatioDatatoFitQCDPi07TeV,"Data/QCD fit to Data","p");
      if (runDrawReweighted) legendFit->AddEntry(histoRatioMCtoDataFitQCDPi07TeV,"MC weighted/QCD fit to Data","p");
      if (histoRatioMCUnweightedtoDataFitQCDPi07TeV) legendFit->AddEntry(histoRatioMCUnweightedtoDataFitQCDPi07TeV,"MC/QCD fit to Data","p");
      legendFit->Draw();
      TLatex *labelRatioMCData = new TLatex(0.2,0.92,collisionSystemPP7TeV.Data());
      SetStyleTLatex(labelRatioMCData, 0.04,4);
      labelRatioMCData->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFraction7TeV->Update();
      canvasFraction7TeV->SaveAs(Form("%s/Pi0_Ratio_pp_MCToDataFit_7TeV.%s",outputDir.Data(),suffix.Data()));
   }
   
   
   
 	TCanvas* canvasFraction2760GeV = new TCanvas("canvasFraction2760GeV","",1550,1200);  // gives the page size
 	canvasFraction2760GeV->SetTickx();
 	canvasFraction2760GeV->SetTicky();
 	canvasFraction2760GeV->SetGridx(0);
 	canvasFraction2760GeV->SetGridy(0);
 	canvasFraction2760GeV->SetLogy(0);
 	canvasFraction2760GeV->SetLeftMargin(0.13);
 	canvasFraction2760GeV->SetRightMargin(0.02);
 	canvasFraction2760GeV->SetTopMargin(0.02);
 	canvasFraction2760GeV->SetFillColor(0);
   
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation PP ***********
   //**********************************************************************************


   canvasCorrectedSpecPi02760GeV->cd();
   canvasCorrectedSpecPi02760GeV->SetLogx(1);
   histo2DCorrPi02760GeV->GetXaxis()->SetRangeUser(0.01,30);
   histo2DCorrPi02760GeV->DrawCopy(); 
   

   TString  nameFinalResDatPi02760GeV = Form("%s/FitResultsMCPi02760GeV.dat",outputDir.Data());
   TString forOutputPi02760GeV;
   fstream fileFinalResultsPi02760GeV;
   fileFinalResultsPi02760GeV.open(nameFinalResDatPi02760GeV.Data(), ios::out);

   TF1* fitYieldDataQCDPi02760GeV = NULL;
   if (directoryPythiaPi02760GeV){
      canvasCorrectedSpecPi02760GeV->cd();
      canvasCorrectedSpecPi02760GeV->SetLogy(1);
      canvasCorrectedSpecPi02760GeV->SetLogx(1);
      histo2DCorrPi02760GeV->GetXaxis()->SetRangeUser(0.01,30);
      histo2DCorrPi02760GeV->DrawCopy(); 
   
      histoMCYieldPi0Pythia->SetMarkerStyle(markerStyleMCPP2760GeV);
      histoMCYieldPi0Pythia->Draw("hist,pe1,same");
   
      histoPi0CorrectedSpecPythia2760GeV->Draw("hist,pe1,same");

      fitYieldDataQCDPi02760GeV = FitObject("qcd","fitYieldDataQCDPi02760GeV","Pi0",histoPi0CorrectedSpecPythia2760GeV,0.4,14,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDPi02760GeV, 1, 1.5, colorPi02760GeV);
      fitYieldDataQCDPi02760GeV->Draw("same");

      forOutputPi02760GeV = WriteParameterToFile(fitYieldDataQCDPi02760GeV);

      fileFinalResultsPi02760GeV << forOutputPi02760GeV.Data() << endl;  

      canvasCorrectedSpecPi02760GeV->Update();
      canvasCorrectedSpecPi02760GeV->Print(Form("%s/Pi0_MCInputSpectraFittedPP.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDPi02760GeV = CalculateHistoRatioToFit(histoPi0CorrectedSpecPythia2760GeV, fitYieldDataQCDPi0Pi02760GeV);
      TH1D* histoRatioMCtoDataFitQCDPi02760GeV = CalculateHistoRatioToFit(histoMCYieldPi0Pythia, fitYieldDataQCDPi0Pi02760GeV);
      TH1D* histoRatioMCUnweightedtoDataFitQCDPi02760GeV = NULL;

      if (histoMCPi0YieldPtPythia2760GeVWOWeights) histoRatioMCUnweightedtoDataFitQCDPi02760GeV = CalculateHistoRatioToFit(histoMCPi0YieldPtPythia2760GeVWOWeights, fitYieldDataQCDPi0Pi02760GeV);
      canvasFraction2760GeV->cd();
      if (histoRatioMCUnweightedtoDataFitQCDPi02760GeV) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPi02760GeV, markerStyleMCPP2760GeV,markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV); 
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDPi02760GeV,  markerStyleMCPP2760GeV+1,markerSizePi0PP2760GeV, colorPi02760GeV+2, colorPi02760GeV+2); 
      DrawGammaSetMarker(histoRatioDatatoFitQCDPi02760GeV, markerStyleMCPP2760GeV,markerSizePi0PP2760GeV, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPi02760GeV,
                  "", "#it{p}_{T} (GeV/#it{c})", "Pi0 Spectrum/ fit to Spectrum",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 13.5);
      histoRatioDatatoFitQCDPi02760GeV->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDPi02760GeV->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDPi02760GeV) histoRatioMCUnweightedtoDataFitQCDPi02760GeV->Draw("same,e,p");  
      TLegend* legendFit = new TLegend(0.16,0.81,0.4,0.9);
      legendFit->SetFillColor(0);
      legendFit->SetLineColor(0);
      legendFit->SetTextSize(0.025);
      legendFit->SetMargin(0.2);
      legendFit->AddEntry(histoRatioDatatoFitQCDPi02760GeV,"Data/QCD fit to Data","p");
      if (runDrawReweighted) legendFit->AddEntry(histoRatioMCtoDataFitQCDPi02760GeV,"MC weighted/QCD fit to Data","p");
      if (histoRatioMCUnweightedtoDataFitQCDPi02760GeV) legendFit->AddEntry(histoRatioMCUnweightedtoDataFitQCDPi02760GeV,"MC/QCD fit to Data","p");
      legendFit->Draw();
      TLatex *labelRatioMCData = new TLatex(0.2,0.92,collisionSystemPP2760GeV.Data());
      SetStyleTLatex(labelRatioMCData, 0.04,4);
      labelRatioMCData->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFraction2760GeV->Update();
      canvasFraction2760GeV->SaveAs(Form("%s/Pi0_Ratio_pp_MCToDataFit_2760GeV.%s",outputDir.Data(),suffix.Data()));
   }
   
   
 	TCanvas* canvasFraction900GeV = new TCanvas("canvasFraction900GeV","",1550,1200);  // gives the page size
 	canvasFraction900GeV->SetTickx();
 	canvasFraction900GeV->SetTicky();
 	canvasFraction900GeV->SetGridx(0);
 	canvasFraction900GeV->SetGridy(0);
 	canvasFraction900GeV->SetLogy(0);
 	canvasFraction900GeV->SetLeftMargin(0.13);
 	canvasFraction900GeV->SetRightMargin(0.02);
 	canvasFraction900GeV->SetTopMargin(0.02);
 	canvasFraction900GeV->SetFillColor(0);
   
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation PP ***********
   //**********************************************************************************


   canvasCorrectedSpecPi0900GeV->cd();
   canvasCorrectedSpecPi0900GeV->SetLogx(1);
   histo2DCorrPi0900GeV->GetXaxis()->SetRangeUser(0.01,30);
   histo2DCorrPi0900GeV->DrawCopy(); 
   

   TString  nameFinalResDatPi0900GeV = Form("%s/FitResultsMCPi0900GeV.dat",outputDir.Data());
   TString forOutputPi0900GeV;
   fstream fileFinalResultsPi0900GeV;
   fileFinalResultsPi0900GeV.open(nameFinalResDatPi0900GeV.Data(), ios::out);

   TF1* fitYieldDataQCDPi0900GeV = NULL;
   if (directoryPythiaPi0900GeV){
      canvasCorrectedSpecPi0900GeV->cd();
      canvasCorrectedSpecPi0900GeV->SetLogy(1);
      canvasCorrectedSpecPi0900GeV->SetLogx(1);
      histo2DCorrPi0900GeV->GetXaxis()->SetRangeUser(0.01,30);
      histo2DCorrPi0900GeV->DrawCopy(); 
   
      histoMCYieldPi0Pythia->SetMarkerStyle(markerStyleMCPP900GeV);
      histoMCYieldPi0Pythia->Draw("hist,pe1,same");
   
      histoPi0CorrectedSpecPythia900GeV->Draw("hist,pe1,same");

      fitYieldDataQCDPi0900GeV = FitObject("qcd","fitYieldDataQCDPi0900GeV","Pi0",histoPi0CorrectedSpecPythia900GeV,0.4,14,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDPi0900GeV, 1, 1.5, colorPi0900GeV);
      fitYieldDataQCDPi0900GeV->Draw("same");

      forOutputPi0900GeV = WriteParameterToFile(fitYieldDataQCDPi0900GeV);

      fileFinalResultsPi0900GeV << forOutputPi0900GeV.Data() << endl;  

      canvasCorrectedSpecPi0900GeV->Update();
      canvasCorrectedSpecPi0900GeV->Print(Form("%s/Pi0_MCInputSpectraFittedPP.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDPi0900GeV = CalculateHistoRatioToFit(histoPi0CorrectedSpecPythia900GeV, fitYieldDataQCDPi0Pi0900GeV);
      TH1D* histoRatioMCtoDataFitQCDPi0900GeV = CalculateHistoRatioToFit(histoMCYieldPi0Pythia, fitYieldDataQCDPi0Pi0900GeV);
      TH1D* histoRatioMCUnweightedtoDataFitQCDPi0900GeV = NULL;

      if (histoMCPi0YieldPtPythia900GeVWOWeights) histoRatioMCUnweightedtoDataFitQCDPi0900GeV = CalculateHistoRatioToFit(histoMCPi0YieldPtPythia900GeVWOWeights, fitYieldDataQCDPi0Pi0900GeV);
      canvasFraction900GeV->cd();
      if (histoRatioMCUnweightedtoDataFitQCDPi0900GeV) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPi0900GeV, markerStyleMCPP900GeV,markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV); 
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDPi0900GeV,  markerStyleMCPP900GeV+1,markerSizePi0PP900GeV, colorPi0900GeV+2, colorPi0900GeV+2); 
      DrawGammaSetMarker(histoRatioDatatoFitQCDPi0900GeV, markerStyleMCPP900GeV,markerSizePi0PP900GeV, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPi0900GeV,
                  "", "#it{p}_{T} (GeV/#it{c})", "Pi0 Spectrum/ fit to Spectrum",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 13.5);
      histoRatioDatatoFitQCDPi0900GeV->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDPi0900GeV->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDPi0900GeV) histoRatioMCUnweightedtoDataFitQCDPi0900GeV->Draw("same,e,p");  
      TLegend* legendFit = new TLegend(0.16,0.81,0.4,0.9);
      legendFit->SetFillColor(0);
      legendFit->SetLineColor(0);
      legendFit->SetTextSize(0.025);
      legendFit->SetMargin(0.2);
      legendFit->AddEntry(histoRatioDatatoFitQCDPi0900GeV,"Data/QCD fit to Data","p");
      if (runDrawReweighted) legendFit->AddEntry(histoRatioMCtoDataFitQCDPi0900GeV,"MC weighted/QCD fit to Data","p");
      if (histoRatioMCUnweightedtoDataFitQCDPi0900GeV) legendFit->AddEntry(histoRatioMCUnweightedtoDataFitQCDPi0900GeV,"MC/QCD fit to Data","p");
      legendFit->Draw();
      TLatex *labelRatioMCData = new TLatex(0.2,0.92,collisionSystemPP900GeV.Data());
      SetStyleTLatex(labelRatioMCData, 0.04,4);
      labelRatioMCData->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFraction900GeV->Update();
      canvasFraction900GeV->SaveAs(Form("%s/Pi0_Ratio_pp_MCToDataFit_900GeV.%s",outputDir.Data(),suffix.Data()));
   }
   
   
   
   cout << "Eta efficiency in different system - Pythia" << endl;   
   
   TCanvas* canvasEfficiencyEtaPythia = new TCanvas("canvasEfficiencyEtaPythia","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyEtaPythia, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffEtaPythia;
   histo2DEffEtaPythia = new TH2F("histo2DEffEtaPythia","histo2DEffEtaPythia",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffEtaPythia, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffEtaPythia->Draw("copy");

//   DrawGammaSetMarker(histoTrueEffPtEta7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
//   histoTrueEffPtEta7TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueEffPtEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
//   histoTrueEffPtEta2760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueEffPtEta900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
//   histoTrueEffPtEta900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaTrueEffiPtPythia8TeV, markerStyleSpectrum8TeVMC, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaTrueEffiPtPythia8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaTrueEffiPtPythia7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
   histoEtaTrueEffiPtPythia7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaTrueEffiPtPythia2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
   histoEtaTrueEffiPtPythia2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoEtaTrueEffiPtPythia900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaTrueEffiPtPythia900GeV->DrawCopy("e1,same");    

   TLegend* legendEffiEtaPythia = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiEtaPythia->SetFillColor(0);
   legendEffiEtaPythia->SetLineColor(0);
   legendEffiEtaPythia->SetTextSize(0.027);
   legendEffiEtaPythia->SetNColumns(2);
   legendEffiEtaPythia->AddEntry(histoEtaTrueEffiPtPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendEffiEtaPythia->AddEntry(histoEtaTrueEffiPtPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendEffiEtaPythia->AddEntry(histoEtaTrueEffiPtPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendEffiEtaPythia->AddEntry(histoEtaTrueEffiPtPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendEffiEtaPythia->AddEntry(histoTrueEffPtEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendEffiEtaPythia->AddEntry(histoTrueEffPtEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendEffiEtaPythia->AddEntry(histoTrueEffPtEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendEffiEtaPythia->Draw();
   
   canvasEfficiencyEtaPythia->SaveAs(Form("%s/EfficiencyEtaPythia.%s",outputDir.Data(),suffix.Data()));
   histo2DEffEtaPythia->Draw("copy");
   
   
   
   cout << "Eta efficiency in different system - Pythia with added signals" << endl;   
   
   TCanvas* canvasEfficiencyEtaPythiaPlus = new TCanvas("canvasEfficiencyEtaPythiaPlus","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyEtaPythiaPlus, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffEtaPythiaPlus;
   histo2DEffEtaPythiaPlus = new TH2F("histo2DEffEtaPythiaPlus","histo2DEffEtaPythiaPlus",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffEtaPythiaPlus, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffEtaPythiaPlus->Draw("copy");

//   DrawGammaSetMarker(histoTrueEffPtEta7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
//   histoTrueEffPtEta7TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueEffPtEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
//   histoTrueEffPtEta2760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueEffPtEta900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
//   histoTrueEffPtEta900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaTrueEffiPtPythiaPlus8TeV, markerStyleSpectrum8TeVMC, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaTrueEffiPtPythiaPlus8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaTrueEffiPtPythiaPlus7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
   histoEtaTrueEffiPtPythiaPlus7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaTrueEffiPtPythiaPlus2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
   histoEtaTrueEffiPtPythiaPlus2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoEtaTrueEffiPtPythiaPlus900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaTrueEffiPtPythiaPlus900GeV->DrawCopy("e1,same");    

   TLegend* legendEffiEtaPythiaPlus = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiEtaPythiaPlus->SetFillColor(0);
   legendEffiEtaPythiaPlus->SetLineColor(0);
   legendEffiEtaPythiaPlus->SetTextSize(0.027);
   legendEffiEtaPythiaPlus->SetNColumns(2);
   legendEffiEtaPythiaPlus->AddEntry(histoEtaTrueEffiPtPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendEffiEtaPythiaPlus->AddEntry(histoEtaTrueEffiPtPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendEffiEtaPythiaPlus->AddEntry(histoEtaTrueEffiPtPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendEffiEtaPythiaPlus->AddEntry(histoEtaTrueEffiPtPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendEffiEtaPythiaPlus->AddEntry(histoTrueEffPtEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendEffiEtaPythiaPlus->AddEntry(histoTrueEffPtEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendEffiEtaPythiaPlus->AddEntry(histoTrueEffPtEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendEffiEtaPythiaPlus->Draw();
   
   canvasEfficiencyEtaPythiaPlus->SaveAs(Form("%s/EfficiencyEtaPythiaPlus.%s",outputDir.Data(),suffix.Data()));
   histo2DEffEtaPythiaPlus->Draw("copy");
   
   
   
   cout << "Eta efficiency in different system - Phojet" << endl;   
   
   TCanvas* canvasEfficiencyEtaPhojet = new TCanvas("canvasEfficiencyEtaPhojet","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyEtaPhojet, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffEtaPhojet;
   histo2DEffEtaPhojet = new TH2F("histo2DEffEtaPhojet","histo2DEffEtaPhojet",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffEtaPhojet, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffEtaPhojet->Draw("copy");

//   DrawGammaSetMarker(histoTrueEffPtEta7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
//   histoTrueEffPtEta7TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueEffPtEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
//   histoTrueEffPtEta2760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueEffPtEta900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
//   histoTrueEffPtEta900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaTrueEffiPtPhojet8TeV, markerStyleSpectrum8TeVMC, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaTrueEffiPtPhojet8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaTrueEffiPtPhojet7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
   histoEtaTrueEffiPtPhojet7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaTrueEffiPtPhojet2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
   histoEtaTrueEffiPtPhojet2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoEtaTrueEffiPtPhojet900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaTrueEffiPtPhojet900GeV->DrawCopy("e1,same");    

   TLegend* legendEffiEtaPhojet = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiEtaPhojet->SetFillColor(0);
   legendEffiEtaPhojet->SetLineColor(0);
   legendEffiEtaPhojet->SetTextSize(0.027);
   legendEffiEtaPhojet->SetNColumns(2);
   legendEffiEtaPhojet->AddEntry(histoEtaTrueEffiPtPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendEffiEtaPhojet->AddEntry(histoEtaTrueEffiPtPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendEffiEtaPhojet->AddEntry(histoEtaTrueEffiPtPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendEffiEtaPhojet->AddEntry(histoEtaTrueEffiPtPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendEffiEtaPhojet->AddEntry(histoTrueEffPtEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendEffiEtaPhojet->AddEntry(histoTrueEffPtEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendEffiEtaPhojet->AddEntry(histoTrueEffPtEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendEffiEtaPhojet->Draw();
   
   canvasEfficiencyEtaPhojet->SaveAs(Form("%s/EfficiencyEtaPhojet.%s",outputDir.Data(),suffix.Data()));
   histo2DEffEtaPhojet->Draw("copy");

   
   //******************************************************************************************************************************************************
   
   cout << "Eta efficiency 8TeV" << endl;   
   
   TCanvas* canvasEfficiencyEta8TeV = new TCanvas("canvasEfficiencyEta8TeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyEta8TeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffEta8TeV;
   histo2DEffEta8TeV = new TH2F("histo2DEffEta8TeV","histo2DEffEta8TeV",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffEta8TeV, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffEta8TeV->Draw("copy");

   DrawGammaSetMarker(histoEtaTrueEffiPtPythia8TeV, markerStyleSpectrum8TeVMC, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaTrueEffiPtPythia8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaTrueEffiPtPythiaPlus8TeV, markerStyleSpectrum8TeVMC+1, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaTrueEffiPtPythiaPlus8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaTrueEffiPtPhojet8TeV, markerStyleSpectrum8TeVMC+2, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaTrueEffiPtPhojet8TeV->DrawCopy("e1,same");
   
   TLegend* legendEffiEta8TeV = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiEta8TeV->SetFillColor(0);
   legendEffiEta8TeV->SetLineColor(0);
   legendEffiEta8TeV->SetTextSize(0.027);
   legendEffiEta8TeV->SetNColumns(2);
   legendEffiEta8TeV->AddEntry(histoEtaTrueEffiPtPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia","p");
   legendEffiEta8TeV->AddEntry(histoEtaTrueEffiPtPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia+AddSign","p");
   legendEffiEta8TeV->AddEntry(histoEtaTrueEffiPtPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV - Phojet","p");
   legendEffiEta8TeV->Draw();
   
   canvasEfficiencyEta8TeV->SaveAs(Form("%s/EfficiencyEta8TeV.%s",outputDir.Data(),suffix.Data()));
   histo2DEffEta8TeV->Draw("copy");
   
   
   cout << "Eta efficiency 7TeV" << endl;   
   
   TCanvas* canvasEfficiencyEta7TeV = new TCanvas("canvasEfficiencyEta7TeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyEta7TeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffEta7TeV;
   histo2DEffEta7TeV = new TH2F("histo2DEffEta7TeV","histo2DEffEta7TeV",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffEta7TeV, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffEta7TeV->Draw("copy");

   DrawGammaSetMarker(histoEtaTrueEffiPtPythia7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV); 
   histoEtaTrueEffiPtPythia7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaTrueEffiPtPythiaPlus7TeV, markerStyleSpectrum7TeVMC+1, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV); 
   histoEtaTrueEffiPtPythiaPlus7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaTrueEffiPtPhojet7TeV, markerStyleSpectrum7TeVMC+2, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV); 
   histoEtaTrueEffiPtPhojet7TeV->DrawCopy("e1,same");
   
   TLegend* legendEffiEta7TeV = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiEta7TeV->SetFillColor(0);
   legendEffiEta7TeV->SetLineColor(0);
   legendEffiEta7TeV->SetTextSize(0.027);
   legendEffiEta7TeV->SetNColumns(2);
   legendEffiEta7TeV->AddEntry(histoEtaTrueEffiPtPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia","p");
   legendEffiEta7TeV->AddEntry(histoEtaTrueEffiPtPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia+AddSign","p");
   legendEffiEta7TeV->AddEntry(histoEtaTrueEffiPtPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV - Phojet","p");
   legendEffiEta7TeV->Draw();
   
   canvasEfficiencyEta7TeV->SaveAs(Form("%s/EfficiencyEta7TeV.%s",outputDir.Data(),suffix.Data()));
   histo2DEffEta7TeV->Draw("copy");
 
   
   cout << "Eta efficiency 2760GeV" << endl;   
   
   TCanvas* canvasEfficiencyEta2760GeV = new TCanvas("canvasEfficiencyEta2760GeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyEta2760GeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffEta2760GeV;
   histo2DEffEta2760GeV = new TH2F("histo2DEffEta2760GeV","histo2DEffEta2760GeV",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffEta2760GeV, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffEta2760GeV->Draw("copy");

   DrawGammaSetMarker(histoEtaTrueEffiPtPythia2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV); 
   histoEtaTrueEffiPtPythia2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaTrueEffiPtPythiaPlus2760GeV, markerStyleSpectrum2760GeVMC+1, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV); 
   histoEtaTrueEffiPtPythiaPlus2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaTrueEffiPtPhojet2760GeV, markerStyleSpectrum2760GeVMC+2, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV); 
   histoEtaTrueEffiPtPhojet2760GeV->DrawCopy("e1,same");
   
   TLegend* legendEffiEta2760GeV = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiEta2760GeV->SetFillColor(0);
   legendEffiEta2760GeV->SetLineColor(0);
   legendEffiEta2760GeV->SetTextSize(0.027);
   legendEffiEta2760GeV->SetNColumns(2);
   legendEffiEta2760GeV->AddEntry(histoEtaTrueEffiPtPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia","p");
   legendEffiEta2760GeV->AddEntry(histoEtaTrueEffiPtPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia+AddSign","p");
   legendEffiEta2760GeV->AddEntry(histoEtaTrueEffiPtPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Phojet","p");
   legendEffiEta2760GeV->Draw();
   
   canvasEfficiencyEta2760GeV->SaveAs(Form("%s/EfficiencyEta2760GeV.%s",outputDir.Data(),suffix.Data()));
   histo2DEffEta2760GeV->Draw("copy");
 
   
   cout << "Eta efficiency 900GeV" << endl;   
   
   TCanvas* canvasEfficiencyEta900GeV = new TCanvas("canvasEfficiencyEta900GeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasEfficiencyEta900GeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DEffEta900GeV;
   histo2DEffEta900GeV = new TH2F("histo2DEffEta900GeV","histo2DEffEta900GeV",1000,0,8,2000,0.e-3,3.7e-3 );
   SetStyleHistoTH2ForGraphs(histo2DEffEta900GeV, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DEffEta900GeV->Draw("copy");

   DrawGammaSetMarker(histoEtaTrueEffiPtPythia900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaTrueEffiPtPythia900GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaTrueEffiPtPythiaPlus900GeV, markerStyleSpectrum900GeVMC+1, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaTrueEffiPtPythiaPlus900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaTrueEffiPtPhojet900GeV, markerStyleSpectrum900GeVMC+2, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaTrueEffiPtPhojet900GeV->DrawCopy("e1,same");
   
   TLegend* legendEffiEta900GeV = new TLegend(0.12,0.78,0.83,0.93);
   legendEffiEta900GeV->SetFillColor(0);
   legendEffiEta900GeV->SetLineColor(0);
   legendEffiEta900GeV->SetTextSize(0.027);
   legendEffiEta900GeV->SetNColumns(2);
   legendEffiEta900GeV->AddEntry(histoEtaTrueEffiPtPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia","p");
   legendEffiEta900GeV->AddEntry(histoEtaTrueEffiPtPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia+AddSign","p");
   legendEffiEta900GeV->AddEntry(histoEtaTrueEffiPtPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Phojet","p");
   legendEffiEta900GeV->Draw();
   
   canvasEfficiencyEta900GeV->SaveAs(Form("%s/EfficiencyEta900GeV.%s",outputDir.Data(),suffix.Data()));
   histo2DEffEta900GeV->Draw("copy");
 
   
   
//**********************************************************************************************************************
   //!!!! histo names for old pp file may not be right !!!!!
   
   cout << "Eta acceptance in different system - Pythia" << endl;   
   
   TCanvas* canvasAcceptanceEtaPythia = new TCanvas("canvasAcceptanceEtaPythia","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptanceEtaPythia, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccEtaPythia;
   histo2DAccEtaPythia = new TH2F("histo2DAccEtaPythia","histo2DAccEtaPythia",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccEtaPythia, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccEtaPythia->Draw("copy");

//   DrawGammaSetMarker(histoTrueAccPtEta7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
//   histoTrueAccPtEta7TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueAccPtEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
//   histoTrueAccPtEta2760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueAccPtEta900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
//   histoTrueAccPtEta900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaAcceptPtPythia8TeV, markerStyleSpectrum8TeVMC, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaAcceptPtPythia8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaAcceptPtPythia7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
   histoEtaAcceptPtPythia7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaAcceptPtPythia2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
   histoEtaAcceptPtPythia2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoEtaAcceptPtPythia900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaAcceptPtPythia900GeV->DrawCopy("e1,same");    

   TLegend* legendAccEtaPythia = new TLegend(0.34,0.13,0.93,0.43);
   legendAccEtaPythia->SetFillColor(0);
   legendAccEtaPythia->SetLineColor(0);
   legendAccEtaPythia->SetTextSize(0.027);
   legendAccEtaPythia->SetNColumns(2);
   legendAccEtaPythia->AddEntry(histoEtaAcceptPtPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendAccEtaPythia->AddEntry(histoEtaAcceptPtPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendAccEtaPythia->AddEntry(histoEtaAcceptPtPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendAccEtaPythia->AddEntry(histoEtaAcceptPtPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendAccEtaPythia->AddEntry(histoTrueAccPtEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendAccEtaPythia->AddEntry(histoTrueAccPtEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendAccEtaPythia->AddEntry(histoTrueAccPtEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendAccEtaPythia->Draw();
   
   canvasAcceptanceEtaPythia->SaveAs(Form("%s/AcceptanceEtaPythia.%s",outputDir.Data(),suffix.Data()));
   histo2DAccEtaPythia->Draw("copy");
   
   
   
   cout << "Eta acceptance in different system - Pythia with added signals" << endl;   
   
   TCanvas* canvasAcceptanceEtaPythiaPlus = new TCanvas("canvasAcceptanceEtaPythiaPlus","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptanceEtaPythiaPlus, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccEtaPythiaPlus;
   histo2DAccEtaPythiaPlus = new TH2F("histo2DAccEtaPythiaPlus","histo2DAccEtaPythiaPlus",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccEtaPythiaPlus, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccEtaPythiaPlus->Draw("copy");

//   DrawGammaSetMarker(histoTrueAccPtEta7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
//   histoTrueAccPtEta7TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueAccPtEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
//   histoTrueAccPtEta2760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueAccPtEta900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
//   histoTrueAccPtEta900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaAcceptPtPythiaPlus8TeV, markerStyleSpectrum8TeVMC, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaAcceptPtPythiaPlus8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaAcceptPtPythiaPlus7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
   histoEtaAcceptPtPythiaPlus7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaAcceptPtPythiaPlus2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
   histoEtaAcceptPtPythiaPlus2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoEtaAcceptPtPythiaPlus900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaAcceptPtPythiaPlus900GeV->DrawCopy("e1,same");    

   TLegend* legendAccEtaPythiaPlus = new TLegend(0.34,0.13,0.93,0.43);
   legendAccEtaPythiaPlus->SetFillColor(0);
   legendAccEtaPythiaPlus->SetLineColor(0);
   legendAccEtaPythiaPlus->SetTextSize(0.027);
   legendAccEtaPythiaPlus->SetNColumns(2);
   legendAccEtaPythiaPlus->AddEntry(histoEtaAcceptPtPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendAccEtaPythiaPlus->AddEntry(histoEtaAcceptPtPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendAccEtaPythiaPlus->AddEntry(histoEtaAcceptPtPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendAccEtaPythiaPlus->AddEntry(histoEtaAcceptPtPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendAccEtaPythiaPlus->AddEntry(histoTrueAccPtEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendAccEtaPythiaPlus->AddEntry(histoTrueAccPtEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendAccEtaPythiaPlus->AddEntry(histoTrueAccPtEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendAccEtaPythiaPlus->Draw();
   
   canvasAcceptanceEtaPythiaPlus->SaveAs(Form("%s/AcceptanceEtaPythiaPlus.%s",outputDir.Data(),suffix.Data()));
   histo2DAccEtaPythiaPlus->Draw("copy");
   
   
   
   cout << "Eta acceptance in different system - Phojet" << endl;   
   
   TCanvas* canvasAcceptanceEtaPhojet = new TCanvas("canvasAcceptanceEtaPhojet","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptanceEtaPhojet, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccEtaPhojet;
   histo2DAccEtaPhojet = new TH2F("histo2DAccEtaPhojet","histo2DAccEtaPhojet",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccEtaPhojet, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccEtaPhojet->Draw("copy");

//   DrawGammaSetMarker(histoTrueAccPtEta7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
//   histoTrueAccPtEta7TeV->DrawCopy("e1,same");  
//   DrawGammaSetMarker(histoTrueAccPtEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
//   histoTrueAccPtEta2760GeV->DrawCopy("pe1,same");    
//   DrawGammaSetMarker(histoTrueAccPtEta900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
//   histoTrueAccPtEta900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaAcceptPtPhojet8TeV, markerStyleSpectrum8TeVMC, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaAcceptPtPhojet8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaAcceptPtPhojet7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV);  
   histoEtaAcceptPtPhojet7TeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaAcceptPtPhojet2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV);  
   histoEtaAcceptPtPhojet2760GeV->DrawCopy("pe1,same");    
   DrawGammaSetMarker(histoEtaAcceptPtPhojet900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaAcceptPtPhojet900GeV->DrawCopy("e1,same");    

   TLegend* legendAccEtaPhojet = new TLegend(0.34,0.13,0.93,0.43);
   legendAccEtaPhojet->SetFillColor(0);
   legendAccEtaPhojet->SetLineColor(0);
   legendAccEtaPhojet->SetTextSize(0.027);
   legendAccEtaPhojet->SetNColumns(2);
   legendAccEtaPhojet->AddEntry(histoEtaAcceptPtPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");
   legendAccEtaPhojet->AddEntry(histoEtaAcceptPtPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   legendAccEtaPhojet->AddEntry(histoEtaAcceptPtPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   legendAccEtaPhojet->AddEntry(histoEtaAcceptPtPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   //legendAccEtaPhojet->AddEntry(histoTrueAccPtEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
   //legendAccEtaPhojet->AddEntry(histoTrueAccPtEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
   //legendAccEtaPhojet->AddEntry(histoTrueAccPtEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
   legendAccEtaPhojet->Draw();
   
   canvasAcceptanceEtaPhojet->SaveAs(Form("%s/AcceptanceEtaPhojet.%s",outputDir.Data(),suffix.Data()));
   histo2DAccEtaPhojet->Draw("copy");

   
   //******************************************************************************************************************************************************
   
   cout << "Eta acceptance 8TeV" << endl;   
   
   TCanvas* canvasAcceptanceEta8TeV = new TCanvas("canvasAcceptanceEta8TeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptanceEta8TeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccEta8TeV;
   histo2DAccEta8TeV = new TH2F("histo2DAccEta8TeV","histo2DAccEta8TeV",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccEta8TeV, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccEta8TeV->Draw("copy");

   DrawGammaSetMarker(histoEtaAcceptPtPythia8TeV, markerStyleSpectrum8TeVMC, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaAcceptPtPythia8TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaAcceptPtPythiaPlus8TeV, markerStyleSpectrum8TeVMC+1, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaAcceptPtPythiaPlus8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaAcceptPtPhojet8TeV, markerStyleSpectrum8TeVMC+2, markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
   histoEtaAcceptPtPhojet8TeV->DrawCopy("e1,same");
   
   TLegend* legendAccEta8TeV = new TLegend(0.34,0.13,0.93,0.43);
   legendAccEta8TeV->SetFillColor(0);
   legendAccEta8TeV->SetLineColor(0);
   legendAccEta8TeV->SetTextSize(0.027);
   legendAccEta8TeV->SetNColumns(2);
   legendAccEta8TeV->AddEntry(histoEtaAcceptPtPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia","p");
   legendAccEta8TeV->AddEntry(histoEtaAcceptPtPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia+AddSign","p");
   legendAccEta8TeV->AddEntry(histoEtaAcceptPtPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV - Phojet","p");
   legendAccEta8TeV->Draw();
   
   canvasAcceptanceEta8TeV->SaveAs(Form("%s/AcceptanceEta8TeV.%s",outputDir.Data(),suffix.Data()));
   histo2DAccEta8TeV->Draw("copy");
   
   
   cout << "Eta acceptance 7TeV" << endl;   
   
   TCanvas* canvasAcceptanceEta7TeV = new TCanvas("canvasAcceptanceEta7TeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptanceEta7TeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccEta7TeV;
   histo2DAccEta7TeV = new TH2F("histo2DAccEta7TeV","histo2DAccEta7TeV",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccEta7TeV, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccEta7TeV->Draw("copy");

   DrawGammaSetMarker(histoEtaAcceptPtPythia7TeV, markerStyleSpectrum7TeVMC, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV); 
   histoEtaAcceptPtPythia7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaAcceptPtPythiaPlus7TeV, markerStyleSpectrum7TeVMC+1, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV); 
   histoEtaAcceptPtPythiaPlus7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaAcceptPtPhojet7TeV, markerStyleSpectrum7TeVMC+2, markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV); 
   histoEtaAcceptPtPhojet7TeV->DrawCopy("e1,same");
   
   TLegend* legendAccEta7TeV = new TLegend(0.34,0.13,0.93,0.43);
   legendAccEta7TeV->SetFillColor(0);
   legendAccEta7TeV->SetLineColor(0);
   legendAccEta7TeV->SetTextSize(0.027);
   legendAccEta7TeV->SetNColumns(2);
   legendAccEta7TeV->AddEntry(histoEtaAcceptPtPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia","p");
   legendAccEta7TeV->AddEntry(histoEtaAcceptPtPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia+AddSign","p");
   legendAccEta7TeV->AddEntry(histoEtaAcceptPtPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV - Phojet","p");
   legendAccEta7TeV->Draw();
   
   canvasAcceptanceEta7TeV->SaveAs(Form("%s/AcceptanceEta7TeV.%s",outputDir.Data(),suffix.Data()));
   histo2DAccEta7TeV->Draw("copy");
 
   
   cout << "Eta acceptance 2760GeV" << endl;   
   
   TCanvas* canvasAcceptanceEta2760GeV = new TCanvas("canvasAcceptanceEta2760GeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptanceEta2760GeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccEta2760GeV;
   histo2DAccEta2760GeV = new TH2F("histo2DAccEta2760GeV","histo2DAccEta2760GeV",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccEta2760GeV, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccEta2760GeV->Draw("copy");

   DrawGammaSetMarker(histoEtaAcceptPtPythia2760GeV, markerStyleSpectrum2760GeVMC, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV); 
   histoEtaAcceptPtPythia2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaAcceptPtPythiaPlus2760GeV, markerStyleSpectrum2760GeVMC+1, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV); 
   histoEtaAcceptPtPythiaPlus2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaAcceptPtPhojet2760GeV, markerStyleSpectrum2760GeVMC+2, markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV); 
   histoEtaAcceptPtPhojet2760GeV->DrawCopy("e1,same");
   
   TLegend* legendAccEta2760GeV = new TLegend(0.34,0.13,0.93,0.43);
   legendAccEta2760GeV->SetFillColor(0);
   legendAccEta2760GeV->SetLineColor(0);
   legendAccEta2760GeV->SetTextSize(0.027);
   legendAccEta2760GeV->SetNColumns(2);
   legendAccEta2760GeV->AddEntry(histoEtaAcceptPtPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia","p");
   legendAccEta2760GeV->AddEntry(histoEtaAcceptPtPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia+AddSign","p");
   legendAccEta2760GeV->AddEntry(histoEtaAcceptPtPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Phojet","p");
   legendAccEta2760GeV->Draw();
   
   canvasAcceptanceEta2760GeV->SaveAs(Form("%s/AcceptanceEta2760GeV.%s",outputDir.Data(),suffix.Data()));
   histo2DAccEta2760GeV->Draw("copy");
 
   
   cout << "Eta acceptance 900GeV" << endl;   
   
   TCanvas* canvasAcceptanceEta900GeV = new TCanvas("canvasAcceptanceEta900GeV","",200,10,1350,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasAcceptanceEta900GeV, 0.1, 0.02, 0.035, 0.09);
   TH2F * histo2DAccEta900GeV;
   histo2DAccEta900GeV = new TH2F("histo2DAccEta900GeV","histo2DAccEta900GeV",1000,0,8.,2000,0.4,1.02 );
   SetStyleHistoTH2ForGraphs(histo2DAccEta900GeV, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
   histo2DAccEta900GeV->Draw("copy");

   DrawGammaSetMarker(histoEtaAcceptPtPythia900GeV, markerStyleSpectrum900GeVMC, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaAcceptPtPythia900GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaAcceptPtPythiaPlus900GeV, markerStyleSpectrum900GeVMC+1, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaAcceptPtPythiaPlus900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaAcceptPtPhojet900GeV, markerStyleSpectrum900GeVMC+2, markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
   histoEtaAcceptPtPhojet900GeV->DrawCopy("e1,same");
   
   TLegend* legendAccEta900GeV = new TLegend(0.34,0.13,0.93,0.43);
   legendAccEta900GeV->SetFillColor(0);
   legendAccEta900GeV->SetLineColor(0);
   legendAccEta900GeV->SetTextSize(0.027);
   legendAccEta900GeV->SetNColumns(2);
   legendAccEta900GeV->AddEntry(histoEtaAcceptPtPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia","p");
   legendAccEta900GeV->AddEntry(histoEtaAcceptPtPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia+AddSign","p");
   legendAccEta900GeV->AddEntry(histoEtaAcceptPtPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Phojet","p");
   legendAccEta900GeV->Draw();
   
   canvasAcceptanceEta900GeV->SaveAs(Form("%s/AcceptanceEta900GeV.%s",outputDir.Data(),suffix.Data()));
   histo2DAccEta900GeV->Draw("copy");
 
   
//************************************************************************************************************************************
   
   
   
   cout << "Eta raw yield in different system - Pythia" << endl;   

   TCanvas* canvasRawYieldEtaPythia = new TCanvas("canvasRawYieldEtaPythia","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldEtaPythia, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawEtaPythia;
   canvasRawYieldEtaPythia->SetLogy();
   histo2DRawEtaPythia = new TH2F("histo2DRawEtaPythia","histo2DRawEtaPythia",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawEtaPythia, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawEtaPythia->Draw("copy");
   
//   DrawGammaSetMarker(histoRawYieldEta7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
//   histoRawYieldEta7TeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldEta2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
//   histoRawYieldEta2760GeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldEta900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
//   histoRawYieldEta900GeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaRawYieldPythia8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaRawYieldPythia8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPythia7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
   histoEtaRawYieldPythia7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaRawYieldPythia2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
   histoEtaRawYieldPythia2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaRawYieldPythia900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
   histoEtaRawYieldPythia900GeV->DrawCopy("e1,same");     

   TLatex *labelRawEtaPythia = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawEtaPythia, 0.038,4);
   labelRawEtaPythia->Draw();

   TLegend* legendRawYieldEtaPythia = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldEtaPythia->SetFillColor(0);
   legendRawYieldEtaPythia->SetLineColor(0);
   legendRawYieldEtaPythia->SetTextSize(0.04);
//   legendRawYieldEtaPythia->AddEntry(histoRawYieldEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
//   legendRawYieldEtaPythia->AddEntry(histoRawYieldEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//   legendRawYieldEtaPythia->AddEntry(histoRawYieldEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 
   legendRawYieldEtaPythia->AddEntry(histoEtaRawYieldPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendRawYieldEtaPythia->AddEntry(histoEtaRawYieldPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendRawYieldEtaPythia->AddEntry(histoEtaRawYieldPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendRawYieldEtaPythia->AddEntry(histoEtaRawYieldPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendRawYieldEtaPythia->Draw();
   
   canvasRawYieldEtaPythia->SaveAs(Form("%s/RawYieldCompEtaPythia.%s",outputDir.Data(),suffix.Data()));


   cout << "Eta raw yield in different system - PythiaPlus" << endl;   

   TCanvas* canvasRawYieldEtaPythiaPlus = new TCanvas("canvasRawYieldEtaPythiaPlus","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldEtaPythiaPlus, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawEtaPythiaPlus;
   canvasRawYieldEtaPythiaPlus->SetLogy();
   histo2DRawEtaPythiaPlus = new TH2F("histo2DRawEtaPythiaPlus","histo2DRawEtaPythiaPlus",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawEtaPythiaPlus, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawEtaPythiaPlus->Draw("copy");
   
//   DrawGammaSetMarker(histoRawYieldEta7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
//   histoRawYieldEta7TeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldEta2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
//   histoRawYieldEta2760GeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldEta900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
//   histoRawYieldEta900GeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaRawYieldPythiaPlus8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaRawYieldPythiaPlus8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPythiaPlus7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
   histoEtaRawYieldPythiaPlus7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaRawYieldPythiaPlus2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
   histoEtaRawYieldPythiaPlus2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaRawYieldPythiaPlus900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
   histoEtaRawYieldPythiaPlus900GeV->DrawCopy("e1,same");     

   TLatex *labelRawEtaPythiaPlus = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawEtaPythiaPlus, 0.038,4);
   labelRawEtaPythiaPlus->Draw();

   TLegend* legendRawYieldEtaPythiaPlus = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldEtaPythiaPlus->SetFillColor(0);
   legendRawYieldEtaPythiaPlus->SetLineColor(0);
   legendRawYieldEtaPythiaPlus->SetTextSize(0.04);
//   legendRawYieldEtaPythiaPlus->AddEntry(histoRawYieldEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
//   legendRawYieldEtaPythiaPlus->AddEntry(histoRawYieldEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//   legendRawYieldEtaPythiaPlus->AddEntry(histoRawYieldEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 
   legendRawYieldEtaPythiaPlus->AddEntry(histoEtaRawYieldPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendRawYieldEtaPythiaPlus->AddEntry(histoEtaRawYieldPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendRawYieldEtaPythiaPlus->AddEntry(histoEtaRawYieldPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendRawYieldEtaPythiaPlus->AddEntry(histoEtaRawYieldPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendRawYieldEtaPythiaPlus->Draw();
   
   canvasRawYieldEtaPythiaPlus->SaveAs(Form("%s/RawYieldCompEtaPythiaPlus.%s",outputDir.Data(),suffix.Data()));


   
   cout << "Eta raw yield in different system - Phojet" << endl;   

   TCanvas* canvasRawYieldEtaPhojet = new TCanvas("canvasRawYieldEtaPhojet","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldEtaPhojet, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawEtaPhojet;
   canvasRawYieldEtaPhojet->SetLogy();
   histo2DRawEtaPhojet = new TH2F("histo2DRawEtaPhojet","histo2DRawEtaPhojet",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawEtaPhojet, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawEtaPhojet->Draw("copy");
   
//   DrawGammaSetMarker(histoRawYieldEta7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
//   histoRawYieldEta7TeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldEta2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
//   histoRawYieldEta2760GeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoRawYieldEta900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
//   histoRawYieldEta900GeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaRawYieldPhojet8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaRawYieldPhojet8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPhojet7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
   histoEtaRawYieldPhojet7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaRawYieldPhojet2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
   histoEtaRawYieldPhojet2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaRawYieldPhojet900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
   histoEtaRawYieldPhojet900GeV->DrawCopy("e1,same");     

   TLatex *labelRawEtaPhojet = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawEtaPhojet, 0.038,4);
   labelRawEtaPhojet->Draw();

   TLegend* legendRawYieldEtaPhojet = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldEtaPhojet->SetFillColor(0);
   legendRawYieldEtaPhojet->SetLineColor(0);
   legendRawYieldEtaPhojet->SetTextSize(0.04);
//   legendRawYieldEtaPhojet->AddEntry(histoRawYieldEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
//   legendRawYieldEtaPhojet->AddEntry(histoRawYieldEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//   legendRawYieldEtaPhojet->AddEntry(histoRawYieldEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 
   legendRawYieldEtaPhojet->AddEntry(histoEtaRawYieldPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendRawYieldEtaPhojet->AddEntry(histoEtaRawYieldPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendRawYieldEtaPhojet->AddEntry(histoEtaRawYieldPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendRawYieldEtaPhojet->AddEntry(histoEtaRawYieldPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendRawYieldEtaPhojet->Draw();
   
   canvasRawYieldEtaPhojet->SaveAs(Form("%s/RawYieldCompEtaPhojet.%s",outputDir.Data(),suffix.Data()));

   
//*************************************************************************************************************************
   

   cout << "Eta raw yield 8TeV " << endl;   

   TCanvas* canvasRawYieldEta8TeV = new TCanvas("canvasRawYieldEta8TeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldEta8TeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawEta8TeV;
   canvasRawYieldEta8TeV->SetLogy();
   histo2DRawEta8TeV = new TH2F("histo2DRawEta8TeV","histo2DRawEta8TeV",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawEta8TeV, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawEta8TeV->Draw("copy");
   
   DrawGammaSetMarker(histoEtaRawYieldPythia8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaRawYieldPythia8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPythiaPlus8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaRawYieldPythiaPlus8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPhojet8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaRawYieldPhojet8TeV->DrawCopy("e1,same");
   
   TLatex *labelRawEta8TeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawEta8TeV, 0.038,4);
   labelRawEta8TeV->Draw();

   TLegend* legendRawYieldEta8TeV = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldEta8TeV->SetFillColor(0);
   legendRawYieldEta8TeV->SetLineColor(0);
   legendRawYieldEta8TeV->SetTextSize(0.04);
   legendRawYieldEta8TeV->AddEntry(histoEtaRawYieldPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia","p");      
   legendRawYieldEta8TeV->AddEntry(histoEtaRawYieldPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia+AddSign","p");      
   legendRawYieldEta8TeV->AddEntry(histoEtaRawYieldPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV - Phojet","p");      
   legendRawYieldEta8TeV->Draw();
   
   canvasRawYieldEta8TeV->SaveAs(Form("%s/RawYieldCompEta8TeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Eta raw yield 7TeV " << endl;   

   TCanvas* canvasRawYieldEta7TeV = new TCanvas("canvasRawYieldEta7TeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldEta7TeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawEta7TeV;
   canvasRawYieldEta7TeV->SetLogy();
   histo2DRawEta7TeV = new TH2F("histo2DRawEta7TeV","histo2DRawEta7TeV",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawEta7TeV, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawEta7TeV->Draw("copy");
   
   DrawGammaSetMarker(histoEtaRawYieldPythia7TeV , markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV);  
   histoEtaRawYieldPythia7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPythiaPlus7TeV , markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV);  
   histoEtaRawYieldPythiaPlus7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPhojet7TeV , markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV);  
   histoEtaRawYieldPhojet7TeV->DrawCopy("e1,same");
   
   TLatex *labelRawEta7TeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawEta7TeV, 0.038,4);
   labelRawEta7TeV->Draw();

   TLegend* legendRawYieldEta7TeV = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldEta7TeV->SetFillColor(0);
   legendRawYieldEta7TeV->SetLineColor(0);
   legendRawYieldEta7TeV->SetTextSize(0.04);
   legendRawYieldEta7TeV->AddEntry(histoEtaRawYieldPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia","p");      
   legendRawYieldEta7TeV->AddEntry(histoEtaRawYieldPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia+AddSign","p");      
   legendRawYieldEta7TeV->AddEntry(histoEtaRawYieldPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV - Phojet","p");      
   legendRawYieldEta7TeV->Draw();
   
   canvasRawYieldEta7TeV->SaveAs(Form("%s/RawYieldCompEta7TeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Eta raw yield 2760GeV " << endl;   

   TCanvas* canvasRawYieldEta2760GeV = new TCanvas("canvasRawYieldEta2760GeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldEta2760GeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawEta2760GeV;
   canvasRawYieldEta2760GeV->SetLogy();
   histo2DRawEta2760GeV = new TH2F("histo2DRawEta2760GeV","histo2DRawEta2760GeV",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawEta2760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawEta2760GeV->Draw("copy");
   
   DrawGammaSetMarker(histoEtaRawYieldPythia2760GeV , markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV);  
   histoEtaRawYieldPythia2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPythiaPlus2760GeV , markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV);  
   histoEtaRawYieldPythiaPlus2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPhojet2760GeV , markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV);  
   histoEtaRawYieldPhojet2760GeV->DrawCopy("e1,same");
   
   TLatex *labelRawEta2760GeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawEta2760GeV, 0.038,4);
   labelRawEta2760GeV->Draw();

   TLegend* legendRawYieldEta2760GeV = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldEta2760GeV->SetFillColor(0);
   legendRawYieldEta2760GeV->SetLineColor(0);
   legendRawYieldEta2760GeV->SetTextSize(0.04);
   legendRawYieldEta2760GeV->AddEntry(histoEtaRawYieldPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia","p");      
   legendRawYieldEta2760GeV->AddEntry(histoEtaRawYieldPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia+AddSign","p");      
   legendRawYieldEta2760GeV->AddEntry(histoEtaRawYieldPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Phojet","p");      
   legendRawYieldEta2760GeV->Draw();
   
   canvasRawYieldEta2760GeV->SaveAs(Form("%s/RawYieldCompEta2760GeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Eta raw yield 900GeV " << endl;   

   TCanvas* canvasRawYieldEta900GeV = new TCanvas("canvasRawYieldEta900GeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasRawYieldEta900GeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DRawEta900GeV;
   canvasRawYieldEta900GeV->SetLogy();
   histo2DRawEta900GeV = new TH2F("histo2DRawEta900GeV","histo2DRawEta900GeV",1000,0.,20,2000,1.e-8,1e-3 );
   SetStyleHistoTH2ForGraphs(histo2DRawEta900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
   histo2DRawEta900GeV->Draw("copy");
   
   DrawGammaSetMarker(histoEtaRawYieldPythia900GeV , markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);  
   histoEtaRawYieldPythia900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPythiaPlus900GeV , markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);  
   histoEtaRawYieldPythiaPlus900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaRawYieldPhojet900GeV , markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);  
   histoEtaRawYieldPhojet900GeV->DrawCopy("e1,same");
   
   TLatex *labelRawEta900GeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelRawEta900GeV, 0.038,4);
   labelRawEta900GeV->Draw();

   TLegend* legendRawYieldEta900GeV = new TLegend(0.65,0.73,0.93,0.93);
   legendRawYieldEta900GeV->SetFillColor(0);
   legendRawYieldEta900GeV->SetLineColor(0);
   legendRawYieldEta900GeV->SetTextSize(0.04);
   legendRawYieldEta900GeV->AddEntry(histoEtaRawYieldPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia","p");      
   legendRawYieldEta900GeV->AddEntry(histoEtaRawYieldPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia+AddSign","p");      
   legendRawYieldEta900GeV->AddEntry(histoEtaRawYieldPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Phojet","p");      
   legendRawYieldEta900GeV->Draw();
   
   canvasRawYieldEta900GeV->SaveAs(Form("%s/RawYieldCompEta900GeV.%s",outputDir.Data(),suffix.Data()));

   
//////////////////////////////// spectra all together //////////////////////////////////////////////////////////////////////////

//************************************************************************************************************************************
   
   
   
   cout << "Eta corrected  yield in different system - Pythia" << endl;   

   TCanvas* canvasCorrectedSpecEtaPythia = new TCanvas("canvasCorrectedSpecEtaPythia","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecEtaPythia, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrEtaPythia;
   canvasCorrectedSpecEtaPythia->SetLogy();
   histo2DCorrEtaPythia = new TH2F("histo2DCorrEtaPythia","histo2DCorrEtaPythia",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrEtaPythia, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrEtaPythia->Draw("copy");
   
   DrawGammaSetMarker(histoEtaCorrectedSpecPythia8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaCorrectedSpecPythia8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPythia7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
   histoEtaCorrectedSpecPythia7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaCorrectedSpecPythia2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
   histoEtaCorrectedSpecPythia2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaCorrectedSpecPythia900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
   histoEtaCorrectedSpecPythia900GeV->DrawCopy("e1,same");     

   TLatex *labelCorrEtaPythia = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrEtaPythia, 0.038,4);
   labelCorrEtaPythia->Draw();

   TLegend* legendCorrectedSpecEtaPythia = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecEtaPythia->SetFillColor(0);
   legendCorrectedSpecEtaPythia->SetLineColor(0);
   legendCorrectedSpecEtaPythia->SetTextSize(0.04); 
   legendCorrectedSpecEtaPythia->AddEntry(histoEtaCorrectedSpecPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendCorrectedSpecEtaPythia->AddEntry(histoEtaCorrectedSpecPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendCorrectedSpecEtaPythia->AddEntry(histoEtaCorrectedSpecPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendCorrectedSpecEtaPythia->AddEntry(histoEtaCorrectedSpecPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendCorrectedSpecEtaPythia->Draw();
   
   canvasCorrectedSpecEtaPythia->SaveAs(Form("%s/CorrectedSpecCompEtaPythia.%s",outputDir.Data(),suffix.Data()));


   cout << "Eta corrected  yield in different system - PythiaPlus" << endl;   

   TCanvas* canvasCorrectedSpecEtaPythiaPlus = new TCanvas("canvasCorrectedSpecEtaPythiaPlus","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecEtaPythiaPlus, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrEtaPythiaPlus;
   canvasCorrectedSpecEtaPythiaPlus->SetLogy();
   histo2DCorrEtaPythiaPlus = new TH2F("histo2DCorrEtaPythiaPlus","histo2DCorrEtaPythiaPlus",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrEtaPythiaPlus, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrEtaPythiaPlus->Draw("copy");
   
//   DrawGammaSetMarker(histoCorrectedSpecEta7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
//   histoCorrectedSpecEta7TeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoCorrectedSpecEta2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
//   histoCorrectedSpecEta2760GeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoCorrectedSpecEta900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
//   histoCorrectedSpecEta900GeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaCorrectedSpecPythiaPlus8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaCorrectedSpecPythiaPlus8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPythiaPlus7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
   histoEtaCorrectedSpecPythiaPlus7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaCorrectedSpecPythiaPlus2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
   histoEtaCorrectedSpecPythiaPlus2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaCorrectedSpecPythiaPlus900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
   histoEtaCorrectedSpecPythiaPlus900GeV->DrawCopy("e1,same");     

   TLatex *labelCorrEtaPythiaPlus = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrEtaPythiaPlus, 0.038,4);
   labelCorrEtaPythiaPlus->Draw();

   TLegend* legendCorrectedSpecEtaPythiaPlus = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecEtaPythiaPlus->SetFillColor(0);
   legendCorrectedSpecEtaPythiaPlus->SetLineColor(0);
   legendCorrectedSpecEtaPythiaPlus->SetTextSize(0.04);
//   legendCorrectedSpecEtaPythiaPlus->AddEntry(histoCorrectedSpecEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
//   legendCorrectedSpecEtaPythiaPlus->AddEntry(histoCorrectedSpecEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//   legendCorrectedSpecEtaPythiaPlus->AddEntry(histoCorrectedSpecEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 
   legendCorrectedSpecEtaPythiaPlus->AddEntry(histoEtaCorrectedSpecPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendCorrectedSpecEtaPythiaPlus->AddEntry(histoEtaCorrectedSpecPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendCorrectedSpecEtaPythiaPlus->AddEntry(histoEtaCorrectedSpecPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendCorrectedSpecEtaPythiaPlus->AddEntry(histoEtaCorrectedSpecPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendCorrectedSpecEtaPythiaPlus->Draw();
   
   canvasCorrectedSpecEtaPythiaPlus->SaveAs(Form("%s/CorrectedSpecCompEtaPythiaPlus.%s",outputDir.Data(),suffix.Data()));


   
   cout << "Eta corrected yield in different system - Phojet" << endl;   

   TCanvas* canvasCorrectedSpecEtaPhojet = new TCanvas("canvasCorrectedSpecEtaPhojet","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecEtaPhojet, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrEtaPhojet;
   canvasCorrectedSpecEtaPhojet->SetLogy();
   histo2DCorrEtaPhojet = new TH2F("histo2DCorrEtaPhojet","histo2DCorrEtaPhojet",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrEtaPhojet, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrEtaPhojet->Draw("copy");
   
//   DrawGammaSetMarker(histoCorrectedSpecEta7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
//   histoCorrectedSpecEta7TeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoCorrectedSpecEta2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
//   histoCorrectedSpecEta2760GeV->DrawCopy("e1,same");   
//   DrawGammaSetMarker(histoCorrectedSpecEta900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
//   histoCorrectedSpecEta900GeV->DrawCopy("e1,same");  
   DrawGammaSetMarker(histoEtaCorrectedSpecPhojet8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaCorrectedSpecPhojet8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPhojet7TeV, markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV); 
   histoEtaCorrectedSpecPhojet7TeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaCorrectedSpecPhojet2760GeV, markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV); 
   histoEtaCorrectedSpecPhojet2760GeV->DrawCopy("e1,same");   
   DrawGammaSetMarker(histoEtaCorrectedSpecPhojet900GeV, markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);   
   histoEtaCorrectedSpecPhojet900GeV->DrawCopy("e1,same");     

   TLatex *labelCorrEtaPhojet = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrEtaPhojet, 0.038,4);
   labelCorrEtaPhojet->Draw();

   TLegend* legendCorrectedSpecEtaPhojet = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecEtaPhojet->SetFillColor(0);
   legendCorrectedSpecEtaPhojet->SetLineColor(0);
   legendCorrectedSpecEtaPhojet->SetTextSize(0.04);
//   legendCorrectedSpecEtaPhojet->AddEntry(histoCorrectedSpecEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
//   legendCorrectedSpecEtaPhojet->AddEntry(histoCorrectedSpecEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//   legendCorrectedSpecEtaPhojet->AddEntry(histoCorrectedSpecEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 
   legendCorrectedSpecEtaPhojet->AddEntry(histoEtaCorrectedSpecPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV","p");      
   legendCorrectedSpecEtaPhojet->AddEntry(histoEtaCorrectedSpecPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");      
   legendCorrectedSpecEtaPhojet->AddEntry(histoEtaCorrectedSpecPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");      
   legendCorrectedSpecEtaPhojet->AddEntry(histoEtaCorrectedSpecPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");      
   legendCorrectedSpecEtaPhojet->Draw();
   
   canvasCorrectedSpecEtaPhojet->SaveAs(Form("%s/CorrectedSpecCompEtaPhojet.%s",outputDir.Data(),suffix.Data()));

   
//*************************************************************************************************************************
   

   cout << "Eta corrected yield 8TeV " << endl;   

   TCanvas* canvasCorrectedSpecEta8TeV = new TCanvas("canvasCorrectedSpecEta8TeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecEta8TeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrEta8TeV;
   canvasCorrectedSpecEta8TeV->SetLogy();
   histo2DCorrEta8TeV = new TH2F("histo2DCorrEta8TeV","histo2DCorrEta8TeV",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrEta8TeV, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrEta8TeV->Draw("copy");
   
   DrawGammaSetMarker(histoEtaCorrectedSpecPythia8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaCorrectedSpecPythia8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPythiaPlus8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaCorrectedSpecPythiaPlus8TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPhojet8TeV , markerStyleSpectrum8TeV, markerSizeEtaPP8TeV*1.5, colorEta8TeV, colorEta8TeV);  
   histoEtaCorrectedSpecPhojet8TeV->DrawCopy("e1,same");
   
   TLatex *labelCorrEta8TeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrEta8TeV, 0.038,4);
   labelCorrEta8TeV->Draw();

   TLegend* legendCorrectedSpecEta8TeV = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecEta8TeV->SetFillColor(0);
   legendCorrectedSpecEta8TeV->SetLineColor(0);
   legendCorrectedSpecEta8TeV->SetTextSize(0.04);
   legendCorrectedSpecEta8TeV->AddEntry(histoEtaCorrectedSpecPythia8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia","p");      
   legendCorrectedSpecEta8TeV->AddEntry(histoEtaCorrectedSpecPythiaPlus8TeV,"pp #sqrt{#it{s}} = 8 TeV - Pythia+AddSign","p");      
   legendCorrectedSpecEta8TeV->AddEntry(histoEtaCorrectedSpecPhojet8TeV,"pp #sqrt{#it{s}} = 8 TeV - Phojet","p");      
   legendCorrectedSpecEta8TeV->Draw();
   
   canvasCorrectedSpecEta8TeV->SaveAs(Form("%s/CorrectedSpecCompEta8TeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Eta corrected yield 7TeV " << endl;   

   TCanvas* canvasCorrectedSpecEta7TeV = new TCanvas("canvasCorrectedSpecEta7TeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecEta7TeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrEta7TeV;
   canvasCorrectedSpecEta7TeV->SetLogy();
   histo2DCorrEta7TeV = new TH2F("histo2DCorrEta7TeV","histo2DCorrEta7TeV",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrEta7TeV, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrEta7TeV->Draw("copy");
   
   DrawGammaSetMarker(histoEtaCorrectedSpecPythia7TeV , markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV);  
   histoEtaCorrectedSpecPythia7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPythiaPlus7TeV , markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV);  
   histoEtaCorrectedSpecPythiaPlus7TeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPhojet7TeV , markerStyleSpectrum7TeV, markerSizeEtaPP7TeV*1.5, colorEta7TeV, colorEta7TeV);  
   histoEtaCorrectedSpecPhojet7TeV->DrawCopy("e1,same");
   
   TLatex *labelCorrEta7TeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrEta7TeV, 0.038,4);
   labelCorrEta7TeV->Draw();

   TLegend* legendCorrectedSpecEta7TeV = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecEta7TeV->SetFillColor(0);
   legendCorrectedSpecEta7TeV->SetLineColor(0);
   legendCorrectedSpecEta7TeV->SetTextSize(0.04);
   legendCorrectedSpecEta7TeV->AddEntry(histoEtaCorrectedSpecPythia7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia","p");      
   legendCorrectedSpecEta7TeV->AddEntry(histoEtaCorrectedSpecPythiaPlus7TeV,"pp #sqrt{#it{s}} = 7 TeV - Pythia+AddSign","p");      
   legendCorrectedSpecEta7TeV->AddEntry(histoEtaCorrectedSpecPhojet7TeV,"pp #sqrt{#it{s}} = 7 TeV - Phojet","p");      
   legendCorrectedSpecEta7TeV->Draw();
   
   canvasCorrectedSpecEta7TeV->SaveAs(Form("%s/CorrectedSpecCompEta7TeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Eta corrected yield 2760GeV " << endl;   

   TCanvas* canvasCorrectedSpecEta2760GeV = new TCanvas("canvasCorrectedSpecEta2760GeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecEta2760GeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrEta2760GeV;
   canvasCorrectedSpecEta2760GeV->SetLogy();
   histo2DCorrEta2760GeV = new TH2F("histo2DCorrEta2760GeV","histo2DCorrEta2760GeV",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrEta2760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrEta2760GeV->Draw("copy");
   
   DrawGammaSetMarker(histoEtaCorrectedSpecPythia2760GeV , markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV);  
   histoEtaCorrectedSpecPythia2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPythiaPlus2760GeV , markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV);  
   histoEtaCorrectedSpecPythiaPlus2760GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPhojet2760GeV , markerStyleSpectrum2760GeV, markerSizeEtaPP2760GeV*1.5, colorEta2760GeV, colorEta2760GeV);  
   histoEtaCorrectedSpecPhojet2760GeV->DrawCopy("e1,same");
   
   TLatex *labelCorrEta2760GeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrEta2760GeV, 0.038,4);
   labelCorrEta2760GeV->Draw();

   TLegend* legendCorrectedSpecEta2760GeV = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecEta2760GeV->SetFillColor(0);
   legendCorrectedSpecEta2760GeV->SetLineColor(0);
   legendCorrectedSpecEta2760GeV->SetTextSize(0.04);
   legendCorrectedSpecEta2760GeV->AddEntry(histoEtaCorrectedSpecPythia2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia","p");      
   legendCorrectedSpecEta2760GeV->AddEntry(histoEtaCorrectedSpecPythiaPlus2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Pythia+AddSign","p");      
   legendCorrectedSpecEta2760GeV->AddEntry(histoEtaCorrectedSpecPhojet2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV - Phojet","p");      
   legendCorrectedSpecEta2760GeV->Draw();
   
   canvasCorrectedSpecEta2760GeV->SaveAs(Form("%s/CorrectedSpecCompEta2760GeV.%s",outputDir.Data(),suffix.Data()));

   
   cout << "Eta corrected yield 900GeV " << endl;   

   TCanvas* canvasCorrectedSpecEta900GeV = new TCanvas("canvasCorrectedSpecEta900GeV","",200,10,1350*1.4,1350);  // gives the page size
   DrawGammaCanvasSettings(canvasCorrectedSpecEta900GeV, 0.12, 0.02, 0.035, 0.09);
   TH2F * histo2DCorrEta900GeV;
   canvasCorrectedSpecEta900GeV->SetLogy();
   histo2DCorrEta900GeV = new TH2F("histo2DCorrEta900GeV","histo2DCorrEta900GeV",1000,0.23,20.,1000,1e-8,2e2 );
   SetStyleHistoTH2ForGraphs(histo2DCorrEta900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
   histo2DCorrEta900GeV->Draw("copy");
   
   DrawGammaSetMarker(histoEtaCorrectedSpecPythia900GeV , markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);  
   histoEtaCorrectedSpecPythia900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPythiaPlus900GeV , markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);  
   histoEtaCorrectedSpecPythiaPlus900GeV->DrawCopy("e1,same");
   DrawGammaSetMarker(histoEtaCorrectedSpecPhojet900GeV , markerStyleSpectrum900GeV, markerSizeEtaPP900GeV*1.5, colorEta900GeV, colorEta900GeV);  
   histoEtaCorrectedSpecPhojet900GeV->DrawCopy("e1,same");
   
   TLatex *labelCorrEta900GeV = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelCorrEta900GeV, 0.038,4);
   labelCorrEta900GeV->Draw();

   TLegend* legendCorrectedSpecEta900GeV = new TLegend(0.65,0.73,0.93,0.93);
   legendCorrectedSpecEta900GeV->SetFillColor(0);
   legendCorrectedSpecEta900GeV->SetLineColor(0);
   legendCorrectedSpecEta900GeV->SetTextSize(0.04);
   legendCorrectedSpecEta900GeV->AddEntry(histoEtaCorrectedSpecPythia900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia","p");      
   legendCorrectedSpecEta900GeV->AddEntry(histoEtaCorrectedSpecPythiaPlus900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Pythia+AddSign","p");      
   legendCorrectedSpecEta900GeV->AddEntry(histoEtaCorrectedSpecPhojet900GeV,"pp #sqrt{#it{s}} = 0.9 TeV - Phojet","p");      
   legendCorrectedSpecEta900GeV->Draw();
   
   canvasCorrectedSpecEta900GeV->SaveAs(Form("%s/CorrectedSpecCompEta900GeV.%s",outputDir.Data(),suffix.Data()));

   
   


 	TCanvas* canvasFractionEta8TeV = new TCanvas("canvasFractionEta8TeV","",1550,1200);  // gives the page size
 	canvasFractionEta8TeV->SetTickx();
 	canvasFractionEta8TeV->SetTicky();
 	canvasFractionEta8TeV->SetGridx(0);
 	canvasFractionEta8TeV->SetGridy(0);
 	canvasFractionEta8TeV->SetLogy(0);
 	canvasFractionEta8TeV->SetLeftMargin(0.13);
 	canvasFractionEta8TeV->SetRightMargin(0.02);
 	canvasFractionEta8TeV->SetTopMargin(0.02);
 	canvasFractionEta8TeV->SetFillColor(0);
   
   //**********************************************************************************
   //**************************** Eta reweighting evalulation PP ***********
   //**********************************************************************************


   canvasCorrectedSpecEta8TeV->cd();
   canvasCorrectedSpecEta8TeV->SetLogx(1);
   histo2DCorrEta8TeV->GetXaxis()->SetRangeUser(0.01,30);
   histo2DCorrEta8TeV->DrawCopy(); 
   

   TString  nameFinalResDatEta8TeV = Form("%s/FitResultsMCEta8TeV.dat",outputDir.Data());
   TString forOutputEta8TeV;
   fstream fileFinalResultsEta8TeV;
   fileFinalResultsEta8TeV.open(nameFinalResDatEta8TeV.Data(), ios::out);

   TF1* fitYieldDataQCDEta8TeV = NULL;
   if (directoryPythiaEta8TeV){
      canvasCorrectedSpecEta8TeV->cd();
      canvasCorrectedSpecEta8TeV->SetLogy(1);
      canvasCorrectedSpecEta8TeV->SetLogx(1);
      histo2DCorrEta8TeV->GetXaxis()->SetRangeUser(0.01,30);
      histo2DCorrEta8TeV->DrawCopy(); 
   
      histoMCYieldEtaPythia->SetMarkerStyle(markerStyleMCPP8TeV);
      histoMCYieldEtaPythia->Draw("hist,pe1,same");
   
      histoEtaCorrectedSpecPythia8TeV->Draw("hist,pe1,same");

      fitYieldDataQCDEta8TeV = FitObject("qcd","fitYieldDataQCDEta8TeV","Eta",histoEtaCorrectedSpecPythia8TeV,0.4,14,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDEta8TeV, 1, 1.5, colorEta8TeV);
      fitYieldDataQCDEta8TeV->Draw("same");

      forOutputEta8TeV = WriteParameterToFile(fitYieldDataQCDEta8TeV);

      fileFinalResultsEta8TeV << forOutputEta8TeV.Data() << endl;  

      canvasCorrectedSpecEta8TeV->Update();
      canvasCorrectedSpecEta8TeV->Print(Form("%s/Eta_MCInputSpectraFittedPP.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDEta8TeV = CalculateHistoRatioToFit(histoEtaCorrectedSpecPythia8TeV, fitYieldDataQCDEtaEta8TeV);
      TH1D* histoRatioMCtoDataFitQCDEta8TeV = CalculateHistoRatioToFit(histoMCYieldEtaPythia, fitYieldDataQCDEtaEta8TeV);
      TH1D* histoRatioMCUnweightedtoDataFitQCDEta8TeV = NULL;

      if (histoMCEtaYieldPtPythiaWOWeights) histoRatioMCUnweightedtoDataFitQCDEta8TeV = CalculateHistoRatioToFit(histoMCEtaYieldPtPythiaWOWeights, fitYieldDataQCDEtaEta8TeV);
      canvasFractionEta8TeV->cd();
      if (histoRatioMCUnweightedtoDataFitQCDEta8TeV) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDEta8TeV, markerStyleMCPP8TeV,markerSizeEtaPP8TeV, colorEta8TeV, colorEta8TeV); 
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDEta8TeV,  markerStyleMCPP8TeV+1,markerSizeEtaPP8TeV, colorEta8TeV+2, colorEta8TeV+2); 
      DrawGammaSetMarker(histoRatioDatatoFitQCDEta8TeV, markerStyleMCPP8TeV,markerSizeEtaPP8TeV, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDEta8TeV,
                  "", "#it{p}_{T} (GeV/#it{c})", "Eta Spectrum/ fit to Spectrum",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 13.5);
      histoRatioDatatoFitQCDEta8TeV->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDEta8TeV->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDEta8TeV) histoRatioMCUnweightedtoDataFitQCDEta8TeV->Draw("same,e,p");  
      TLegend* legendFit = new TLegend(0.16,0.81,0.4,0.9);
      legendFit->SetFillColor(0);
      legendFit->SetLineColor(0);
      legendFit->SetTextSize(0.025);
      legendFit->SetMargin(0.2);
      legendFit->AddEntry(histoRatioDatatoFitQCDEta8TeV,"Data/QCD fit to Data","p");
      if (runDrawReweighted) legendFit->AddEntry(histoRatioMCtoDataFitQCDEta8TeV,"MC weighted/QCD fit to Data","p");
      if (histoRatioMCUnweightedtoDataFitQCDEta8TeV) legendFit->AddEntry(histoRatioMCUnweightedtoDataFitQCDEta8TeV,"MC/QCD fit to Data","p");
      legendFit->Draw();
      TLatex *labelRatioMCData = new TLatex(0.2,0.92,collisionSystemPP8TeV.Data());
      SetStyleTLatex(labelRatioMCData, 0.04,4);
      labelRatioMCData->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFractionEta8TeV->Update();
      canvasFractionEta8TeV->SaveAs(Form("%s/Eta_Ratio_pp_MCToDataFit_8TeV.%s",outputDir.Data(),suffix.Data()));
   }
   

 	TCanvas* canvasFractionEta7TeV = new TCanvas("canvasFractionEta7TeV","",1550,1200);  // gives the page size
 	canvasFractionEta7TeV->SetTickx();
 	canvasFractionEta7TeV->SetTicky();
 	canvasFractionEta7TeV->SetGridx(0);
 	canvasFractionEta7TeV->SetGridy(0);
 	canvasFractionEta7TeV->SetLogy(0);
 	canvasFractionEta7TeV->SetLeftMargin(0.13);
 	canvasFractionEta7TeV->SetRightMargin(0.02);
 	canvasFractionEta7TeV->SetTopMargin(0.02);
 	canvasFractionEta7TeV->SetFillColor(0);
   
   //**********************************************************************************
   //**************************** Eta reweighting evalulation PP ***********
   //**********************************************************************************


   canvasCorrectedSpecEta7TeV->cd();
   canvasCorrectedSpecEta7TeV->SetLogx(1);
   histo2DCorrEta7TeV->GetXaxis()->SetRangeUser(0.01,30);
   histo2DCorrEta7TeV->DrawCopy(); 
   

   TString  nameFinalResDatEta7TeV = Form("%s/FitResultsMCEta7TeV.dat",outputDir.Data());
   TString forOutputEta7TeV;
   fstream fileFinalResultsEta7TeV;
   fileFinalResultsEta7TeV.open(nameFinalResDatEta7TeV.Data(), ios::out);

   TF1* fitYieldDataQCDEta7TeV = NULL;
   if (directoryPythiaEta7TeV){
      canvasCorrectedSpecEta7TeV->cd();
      canvasCorrectedSpecEta7TeV->SetLogy(1);
      canvasCorrectedSpecEta7TeV->SetLogx(1);
      histo2DCorrEta7TeV->GetXaxis()->SetRangeUser(0.01,30);
      histo2DCorrEta7TeV->DrawCopy(); 
   
      histoMCYieldEtaPythia->SetMarkerStyle(markerStyleMCPP7TeV);
      histoMCYieldEtaPythia->Draw("hist,pe1,same");
   
      histoEtaCorrectedSpecPythia7TeV->Draw("hist,pe1,same");

      fitYieldDataQCDEta7TeV = FitObject("qcd","fitYieldDataQCDEta7TeV","Eta",histoEtaCorrectedSpecPythia7TeV,0.4,14,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDEta7TeV, 1, 1.5, colorEta7TeV);
      fitYieldDataQCDEta7TeV->Draw("same");

      forOutputEta7TeV = WriteParameterToFile(fitYieldDataQCDEta7TeV);

      fileFinalResultsEta7TeV << forOutputEta7TeV.Data() << endl;  

      canvasCorrectedSpecEta7TeV->Update();
      canvasCorrectedSpecEta7TeV->Print(Form("%s/Eta_MCInputSpectraFittedPP.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDEta7TeV = CalculateHistoRatioToFit(histoEtaCorrectedSpecPythia7TeV, fitYieldDataQCDEtaEta7TeV);
      TH1D* histoRatioMCtoDataFitQCDEta7TeV = CalculateHistoRatioToFit(histoMCYieldEtaPythia, fitYieldDataQCDEtaEta7TeV);
      TH1D* histoRatioMCUnweightedtoDataFitQCDEta7TeV = NULL;

      if (histoMCEtaYieldPtPythiaWOWeights) histoRatioMCUnweightedtoDataFitQCDEta7TeV = CalculateHistoRatioToFit(histoMCEtaYieldPtPythiaWOWeights, fitYieldDataQCDEtaEta7TeV);
      canvasFractionEta7TeV->cd();
      if (histoRatioMCUnweightedtoDataFitQCDEta7TeV) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDEta7TeV, markerStyleMCPP7TeV,markerSizeEtaPP7TeV, colorEta7TeV, colorEta7TeV); 
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDEta7TeV,  markerStyleMCPP7TeV+1,markerSizeEtaPP7TeV, colorEta7TeV+2, colorEta7TeV+2); 
      DrawGammaSetMarker(histoRatioDatatoFitQCDEta7TeV, markerStyleMCPP7TeV,markerSizeEtaPP7TeV, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDEta7TeV,
                  "", "#it{p}_{T} (GeV/#it{c})", "Eta Spectrum/ fit to Spectrum",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 13.5);
      histoRatioDatatoFitQCDEta7TeV->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDEta7TeV->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDEta7TeV) histoRatioMCUnweightedtoDataFitQCDEta7TeV->Draw("same,e,p");  
      TLegend* legendFit = new TLegend(0.16,0.81,0.4,0.9);
      legendFit->SetFillColor(0);
      legendFit->SetLineColor(0);
      legendFit->SetTextSize(0.025);
      legendFit->SetMargin(0.2);
      legendFit->AddEntry(histoRatioDatatoFitQCDEta7TeV,"Data/QCD fit to Data","p");
      if (runDrawReweighted) legendFit->AddEntry(histoRatioMCtoDataFitQCDEta7TeV,"MC weighted/QCD fit to Data","p");
      if (histoRatioMCUnweightedtoDataFitQCDEta7TeV) legendFit->AddEntry(histoRatioMCUnweightedtoDataFitQCDEta7TeV,"MC/QCD fit to Data","p");
      legendFit->Draw();
      TLatex *labelRatioMCData = new TLatex(0.2,0.92,collisionSystemPP7TeV.Data());
      SetStyleTLatex(labelRatioMCData, 0.04,4);
      labelRatioMCData->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFractionEta7TeV->Update();
      canvasFractionEta7TeV->SaveAs(Form("%s/Eta_Ratio_pp_MCToDataFit_7TeV.%s",outputDir.Data(),suffix.Data()));
   }
   
   
   
 	TCanvas* canvasFractionEta2760GeV = new TCanvas("canvasFractionEta2760GeV","",1550,1200);  // gives the page size
 	canvasFractionEta2760GeV->SetTickx();
 	canvasFractionEta2760GeV->SetTicky();
 	canvasFractionEta2760GeV->SetGridx(0);
 	canvasFractionEta2760GeV->SetGridy(0);
 	canvasFractionEta2760GeV->SetLogy(0);
 	canvasFractionEta2760GeV->SetLeftMargin(0.13);
 	canvasFractionEta2760GeV->SetRightMargin(0.02);
 	canvasFractionEta2760GeV->SetTopMargin(0.02);
 	canvasFractionEta2760GeV->SetFillColor(0);
   
   //**********************************************************************************
   //**************************** Eta reweighting evalulation PP ***********
   //**********************************************************************************


   canvasCorrectedSpecEta2760GeV->cd();
   canvasCorrectedSpecEta2760GeV->SetLogx(1);
   histo2DCorrEta2760GeV->GetXaxis()->SetRangeUser(0.01,30);
   histo2DCorrEta2760GeV->DrawCopy(); 
   

   TString  nameFinalResDatEta2760GeV = Form("%s/FitResultsMCEta2760GeV.dat",outputDir.Data());
   TString forOutputEta2760GeV;
   fstream fileFinalResultsEta2760GeV;
   fileFinalResultsEta2760GeV.open(nameFinalResDatEta2760GeV.Data(), ios::out);

   TF1* fitYieldDataQCDEta2760GeV = NULL;
   if (directoryPythiaEta2760GeV){
      canvasCorrectedSpecEta2760GeV->cd();
      canvasCorrectedSpecEta2760GeV->SetLogy(1);
      canvasCorrectedSpecEta2760GeV->SetLogx(1);
      histo2DCorrEta2760GeV->GetXaxis()->SetRangeUser(0.01,30);
      histo2DCorrEta2760GeV->DrawCopy(); 
   
      histoMCYieldEtaPythia->SetMarkerStyle(markerStyleMCPP2760GeV);
      histoMCYieldEtaPythia->Draw("hist,pe1,same");
   
      histoEtaCorrectedSpecPythia2760GeV->Draw("hist,pe1,same");

      fitYieldDataQCDEta2760GeV = FitObject("qcd","fitYieldDataQCDEta2760GeV","Eta",histoEtaCorrectedSpecPythia2760GeV,0.4,14,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDEta2760GeV, 1, 1.5, colorEta2760GeV);
      fitYieldDataQCDEta2760GeV->Draw("same");

      forOutputEta2760GeV = WriteParameterToFile(fitYieldDataQCDEta2760GeV);

      fileFinalResultsEta2760GeV << forOutputEta2760GeV.Data() << endl;  

      canvasCorrectedSpecEta2760GeV->Update();
      canvasCorrectedSpecEta2760GeV->Print(Form("%s/Eta_MCInputSpectraFittedPP.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDEta2760GeV = CalculateHistoRatioToFit(histoEtaCorrectedSpecPythia2760GeV, fitYieldDataQCDEtaEta2760GeV);
      TH1D* histoRatioMCtoDataFitQCDEta2760GeV = CalculateHistoRatioToFit(histoMCYieldEtaPythia, fitYieldDataQCDEtaEta2760GeV);
      TH1D* histoRatioMCUnweightedtoDataFitQCDEta2760GeV = NULL;

      if (histoMCEtaYieldPtPythiaWOWeights) histoRatioMCUnweightedtoDataFitQCDEta2760GeV = CalculateHistoRatioToFit(histoMCEtaYieldPtPythiaWOWeights, fitYieldDataQCDEtaEta2760GeV);
      canvasFractionEta2760GeV->cd();
      if (histoRatioMCUnweightedtoDataFitQCDEta2760GeV) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDEta2760GeV, markerStyleMCPP2760GeV,markerSizeEtaPP2760GeV, colorEta2760GeV, colorEta2760GeV); 
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDEta2760GeV,  markerStyleMCPP2760GeV+1,markerSizeEtaPP2760GeV, colorEta2760GeV+2, colorEta2760GeV+2); 
      DrawGammaSetMarker(histoRatioDatatoFitQCDEta2760GeV, markerStyleMCPP2760GeV,markerSizeEtaPP2760GeV, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDEta2760GeV,
                  "", "#it{p}_{T} (GeV/#it{c})", "Eta Spectrum/ fit to Spectrum",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 13.5);
      histoRatioDatatoFitQCDEta2760GeV->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDEta2760GeV->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDEta2760GeV) histoRatioMCUnweightedtoDataFitQCDEta2760GeV->Draw("same,e,p");  
      TLegend* legendFit = new TLegend(0.16,0.81,0.4,0.9);
      legendFit->SetFillColor(0);
      legendFit->SetLineColor(0);
      legendFit->SetTextSize(0.025);
      legendFit->SetMargin(0.2);
      legendFit->AddEntry(histoRatioDatatoFitQCDEta2760GeV,"Data/QCD fit to Data","p");
      if (runDrawReweighted) legendFit->AddEntry(histoRatioMCtoDataFitQCDEta2760GeV,"MC weighted/QCD fit to Data","p");
      if (histoRatioMCUnweightedtoDataFitQCDEta2760GeV) legendFit->AddEntry(histoRatioMCUnweightedtoDataFitQCDEta2760GeV,"MC/QCD fit to Data","p");
      legendFit->Draw();
      TLatex *labelRatioMCData = new TLatex(0.2,0.92,collisionSystemPP2760GeV.Data());
      SetStyleTLatex(labelRatioMCData, 0.04,4);
      labelRatioMCData->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFractionEta2760GeV->Update();
      canvasFractionEta2760GeV->SaveAs(Form("%s/Eta_Ratio_pp_MCToDataFit_2760GeV.%s",outputDir.Data(),suffix.Data()));
   }
   
   
 	TCanvas* canvasFractionEta900GeV = new TCanvas("canvasFractionEta900GeV","",1550,1200);  // gives the page size
 	canvasFractionEta900GeV->SetTickx();
 	canvasFractionEta900GeV->SetTicky();
 	canvasFractionEta900GeV->SetGridx(0);
 	canvasFractionEta900GeV->SetGridy(0);
 	canvasFractionEta900GeV->SetLogy(0);
 	canvasFractionEta900GeV->SetLeftMargin(0.13);
 	canvasFractionEta900GeV->SetRightMargin(0.02);
 	canvasFractionEta900GeV->SetTopMargin(0.02);
 	canvasFractionEta900GeV->SetFillColor(0);
   
   //**********************************************************************************
   //**************************** Eta reweighting evalulation PP ***********
   //**********************************************************************************


   canvasCorrectedSpecEta900GeV->cd();
   canvasCorrectedSpecEta900GeV->SetLogx(1);
   histo2DCorrEta900GeV->GetXaxis()->SetRangeUser(0.01,30);
   histo2DCorrEta900GeV->DrawCopy(); 
   

   TString  nameFinalResDatEta900GeV = Form("%s/FitResultsMCEta900GeV.dat",outputDir.Data());
   TString forOutputEta900GeV;
   fstream fileFinalResultsEta900GeV;
   fileFinalResultsEta900GeV.open(nameFinalResDatEta900GeV.Data(), ios::out);

   TF1* fitYieldDataQCDEta900GeV = NULL;
   if (directoryPythiaEta900GeV){
      canvasCorrectedSpecEta900GeV->cd();
      canvasCorrectedSpecEta900GeV->SetLogy(1);
      canvasCorrectedSpecEta900GeV->SetLogx(1);
      histo2DCorrEta900GeV->GetXaxis()->SetRangeUser(0.01,30);
      histo2DCorrEta900GeV->DrawCopy(); 
   
      histoMCYieldEtaPythia->SetMarkerStyle(markerStyleMCPP900GeV);
      histoMCYieldEtaPythia->Draw("hist,pe1,same");
   
      histoEtaCorrectedSpecPythia900GeV->Draw("hist,pe1,same");

      fitYieldDataQCDEta900GeV = FitObject("qcd","fitYieldDataQCDEta900GeV","Eta",histoEtaCorrectedSpecPythia900GeV,0.4,14,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDEta900GeV, 1, 1.5, colorEta900GeV);
      fitYieldDataQCDEta900GeV->Draw("same");

      forOutputEta900GeV = WriteParameterToFile(fitYieldDataQCDEta900GeV);

      fileFinalResultsEta900GeV << forOutputEta900GeV.Data() << endl;  

      canvasCorrectedSpecEta900GeV->Update();
      canvasCorrectedSpecEta900GeV->Print(Form("%s/Eta_MCInputSpectraFittedPP.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDEta900GeV = CalculateHistoRatioToFit(histoEtaCorrectedSpecPythia900GeV, fitYieldDataQCDEtaEta900GeV);
      TH1D* histoRatioMCtoDataFitQCDEta900GeV = CalculateHistoRatioToFit(histoMCYieldEtaPythia, fitYieldDataQCDEtaEta900GeV);
      TH1D* histoRatioMCUnweightedtoDataFitQCDEta900GeV = NULL;

      if (histoMCEtaYieldPtPythiaWOWeights) histoRatioMCUnweightedtoDataFitQCDEta900GeV = CalculateHistoRatioToFit(histoMCEtaYieldPtPythiaWOWeights, fitYieldDataQCDEtaEta900GeV);
      canvasFractionEta900GeV->cd();
      if (histoRatioMCUnweightedtoDataFitQCDEta900GeV) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDEta900GeV, markerStyleMCPP900GeV,markerSizeEtaPP900GeV, colorEta900GeV, colorEta900GeV); 
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDEta900GeV,  markerStyleMCPP900GeV+1,markerSizeEtaPP900GeV, colorEta900GeV+2, colorEta900GeV+2); 
      DrawGammaSetMarker(histoRatioDatatoFitQCDEta900GeV, markerStyleMCPP900GeV,markerSizeEtaPP900GeV, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDEta900GeV,
                  "", "#it{p}_{T} (GeV/#it{c})", "Eta Spectrum/ fit to Spectrum",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 13.5);
      histoRatioDatatoFitQCDEta900GeV->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDEta900GeV->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDEta900GeV) histoRatioMCUnweightedtoDataFitQCDEta900GeV->Draw("same,e,p");  
      TLegend* legendFit = new TLegend(0.16,0.81,0.4,0.9);
      legendFit->SetFillColor(0);
      legendFit->SetLineColor(0);
      legendFit->SetTextSize(0.025);
      legendFit->SetMargin(0.2);
      legendFit->AddEntry(histoRatioDatatoFitQCDEta900GeV,"Data/QCD fit to Data","p");
      if (runDrawReweighted) legendFit->AddEntry(histoRatioMCtoDataFitQCDEta900GeV,"MC weighted/QCD fit to Data","p");
      if (histoRatioMCUnweightedtoDataFitQCDEta900GeV) legendFit->AddEntry(histoRatioMCUnweightedtoDataFitQCDEta900GeV,"MC/QCD fit to Data","p");
      legendFit->Draw();
      TLatex *labelRatioMCData = new TLatex(0.2,0.92,collisionSystemPP900GeV.Data());
      SetStyleTLatex(labelRatioMCData, 0.04,4);
      labelRatioMCData->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFractionEta900GeV->Update();
      canvasFractionEta900GeV->SaveAs(Form("%s/Eta_Ratio_pp_MCToDataFit_900GeV.%s",outputDir.Data(),suffix.Data()));
   }
   

   
   TFile fMCSpectraInput("MCSpectraInputpp.root","UPDATE");

      if (fitYieldDataQCDPi08TeV){
         fitYieldDataQCDPi08TeV->SetRange(0,30);
         fitYieldDataQCDPi08TeV->Write("Pi0_Fit_Data_PP_8TeV",TObject::kOverwrite);
      }
      if (fitYieldDataQCDPi07TeV){
         fitYieldDataQCDPi07TeV->SetRange(0,30);
         fitYieldDataQCDPi07TeV->Write("Pi0_Fit_Data_PP_7TeV",TObject::kOverwrite);
      }
      if (fitYieldDataQCDPi02760GeV){
         fitYieldDataQCDPi02760GeV->SetRange(0,30);
         fitYieldDataQCDPi02760GeV->Write("Pi0_Fit_Data_PP_2760GeV",TObject::kOverwrite);
      }
      if (fitYieldDataQCDPi0900GeV){
         fitYieldDataQCDPi0900GeV->SetRange(0,30);
         fitYieldDataQCDPi0900GeV->Write("Pi0_Fit_Data_PP_900GeV",TObject::kOverwrite);
      }
      if (fitYieldDataQCDEta8TeV){
         fitYieldDataQCDEta8TeV->SetRange(0,30);
         fitYieldDataQCDEta8TeV->Write("Eta_Fit_Data_PP_8TeV",TObject::kOverwrite);
      }
      if (fitYieldDataQCDEta7TeV){
         fitYieldDataQCDEta7TeV->SetRange(0,30);
         fitYieldDataQCDEta7TeV->Write("Eta_Fit_Data_PP_7TeV",TObject::kOverwrite);
      }
      if (fitYieldDataQCDEta2760GeV){
         fitYieldDataQCDEta2760GeV->SetRange(0,30);
         fitYieldDataQCDEta2760GeV->Write("Eta_Fit_Data_PP_2760GeV",TObject::kOverwrite);
      }
      if (fitYieldDataQCDEta900GeV){
         fitYieldDataQCDEta900GeV->SetRange(0,30);
         fitYieldDataQCDEta900GeV->Write("Eta_Fit_Data_PP_900GeV",TObject::kOverwrite);
      }

      
      histoMCPi0YieldPtPythia8TeVWOWeights->Write("Pi0_PP_8TeV",TObject::kOverwrite);
      histoMCPi0YieldPtPythiaPlus8AddedSigWOWeights->Write("Pi0_addSig_PP_8TeV",TObject::kOverwrite);
      histoMCPi0YieldPtPythia7TeVWOWeights->Write("Pi0_PP_7TeV",TObject::kOverwrite);
      histoMCPi0YieldPtPythiaPlus7AddedSigWOWeights->Write("Pi0_addSig_PP_7TeV",TObject::kOverwrite);
      histoMCPi0YieldPtPythia2760GeVWOWeights->Write("Pi0_PP_2760GeV",TObject::kOverwrite);
      histoMCPi0YieldPtPythiaPlus2760GeVAddedSigWOWeights->Write("Pi0_addSig_PP_2760GeV",TObject::kOverwrite);
      histoMCPi0YieldPtPythia900GeVWOWeights->Write("Pi0_PP_900GeV",TObject::kOverwrite);
      histoMCPi0YieldPtPythiaPlus900GeVAddedSigWOWeights->Write("Pi0_addSig_PP_900GeV",TObject::kOverwrite);

      histoMCEtaYieldPtPythia8TeVWOWeights->Write("Eta_PP_8TeV",TObject::kOverwrite);
      histoMCEtaYieldPtPythiaPlus8AddedSigWOWeights->Write("Eta_addSig_PP_8TeV",TObject::kOverwrite);
      histoMCEtaYieldPtPythia7TeVWOWeights->Write("Eta_PP_7TeV",TObject::kOverwrite);
      histoMCEtaYieldPtPythiaPlus7AddedSigWOWeights->Write("Eta_addSig_PP_7TeV",TObject::kOverwrite);
      histoMCEtaYieldPtPythia2760GeVWOWeights->Write("Eta_PP_2760GeV",TObject::kOverwrite);
      histoMCEtaYieldPtPythiaPlus2760GeVAddedSigWOWeights->Write("Eta_addSig_PP_2760GeV",TObject::kOverwrite);
      histoMCEtaYieldPtPythia900GeVWOWeights->Write("Eta_PP_900GeV",TObject::kOverwrite);
      histoMCEtaYieldPtPythiaPlus900GeVAddedSigWOWeights->Write("Eta_addSig_PP_900GeV",TObject::kOverwrite);

      
      
   fMCSpectraInput.Close();
}

