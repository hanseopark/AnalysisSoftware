/****************************************************************************************************************************
******      Friederike Bock, friederike.bock@cern.ch                                      *****
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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

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


void CombineNeutralPionResultspPb(TString suffix = "pdf", TString nameFilepPb = "data_PCMResults_pPb.root", Bool_t runDrawReweighted = kTRUE){


	gROOT->Reset();   
	gROOT->SetStyle("Plain");

	TString dateForOutput = ReturnDateStringForOutput();
	TString outputDir = Form("%s/%s/CombineNeutralPionResultspPb",suffix.Data(),dateForOutput.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	gSystem->Exec(Form("cp %s %s/InputFilePCMPionpPb.root ",nameFilepPb.Data(),outputDir.Data() ));

	StyleSettingsThesis();  
	SetPlotStyle();

	Color_t  colorCombpPb         	   = kBlack;
	Color_t  colorCombpPb0020          = kRed+1;
	Color_t  colorCombpPb2040          = kGreen+2;
	Color_t  colorCombpPb4060          = kCyan+2;
	Color_t  colorCombpPb6080          = kBlue+1;
	Color_t  colorCombpPb60100         = kViolet+3;

	Color_t  colorCombMCpPb 	      = kGray+1;
	Color_t  colorCombMCpPb0020       = kRed-6;
	Color_t  colorCombMCpPb2040       = kGreen-6;
	Color_t  colorCombMCpPb4060       = kCyan-6;
	Color_t  colorCombMCpPb6080       = kBlue-6;
	Color_t  colorCombMCpPb60100      = kViolet+6;

	Color_t  colorCombMCReweightedpPb 	      = kGray+2;
	Color_t  colorCombMCReweightedAddSigpPb   = kRed+2;
	Color_t  colorCombMCReweightedpPb0020       = kRed-6;
	Color_t  colorCombMCReweightedpPb2040       = kGreen-6;
	Color_t  colorCombMCReweightedpPb4060       = kCyan-6;
	Color_t  colorCombMCReweightedpPb6080       = kBlue-6;
	Color_t  colorCombMCReweightedpPb60100      = kViolet+6;

	Style_t  markerStylepPb  = 20 ;
	Style_t  markerStylepPb0020   = 20 ;
	Style_t  markerStylepPb2040   = 21 ;
	Style_t  markerStylepPb4060   = 29 ;
	Style_t  markerStylepPb6080   = 33 ;
	Style_t  markerStylepPb60100 = 34 ;

	Style_t  lineStyleMCpPb       = 1;
	Style_t  lineStyleMCpPb0020   = 1;
	Style_t  lineStyleMCpPb2040   = 1 ;
	Style_t  lineStyleMCpPb4060   = 1 ;
	Style_t  lineStyleMCpPb6080   = 1 ;
	Style_t  lineStyleMCpPb60100  = 1 ;

	Style_t  lineStyleMCAddSigpPb       = 2;
	Style_t  lineStyleMCAddSigpPb0020   = 3;
	Style_t  lineStyleMCAddSigpPb2040   = 4 ;
	Style_t  lineStyleMCAddSigpPb4060   = 5 ;
	Style_t  lineStyleMCAddSigpPb6080   = 6 ;
	Style_t  lineStyleMCAddSigpPb60100  = 7 ;

	Size_t   markerSizepPb  = 2.;
	Size_t   markerSizepPb0020  = 2.;
	Size_t   markerSizepPb2040  = 2.;
	Size_t   markerSizepPb4060  = 2.5;
	Size_t   markerSizepPb6080  = 2.5;
	Size_t   markerSizepPb60100  = 2.;

	TString collisionSystempPb = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.023 TeV";     
	TString collisionSystempPb0020 = "0-20% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.023 TeV";		
	TString collisionSystempPb2040 = "20-40% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.023 TeV";      
	TString collisionSystempPb4060 = "40-60% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.023 TeV";    
	TString collisionSystempPb60100 = "60-100% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.023 TeV";     
	TString collisionSystempPb6080 = "60-80% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.023 TeV";     

	TString nameHistoPCM = "CorrectedYieldPi0";
	TString nameHistoPCMEta = "CorrectedYieldEta";
	TString nameGraphPCM = "Pi0SystError";
	TString nameGraphPCMEta = "EtaSystError";

	cout << "PCM" << endl;
	cout << "0-100%" << endl;
	TFile* filePCMpPb = 					new TFile(nameFilepPb);
	TDirectory*	directoryPCMPi0pPb = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
	TH1D* histoPCMNumberOfEventspPb= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV0-100%");
	TH1D* histoPCMYieldPi0pPb = 				(TH1D*)directoryPCMPi0pPb->Get(nameHistoPCM.Data());
	TH1D* histoPCMRAWYieldPi0pPb =            (TH1D*)directoryPCMPi0pPb->Get("RAWYieldPerEventsPi0");
	TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb= 	(TGraphAsymmErrors*)directoryPCMPi0pPb->Get(nameGraphPCM.Data());	
	TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb= 	(TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystErrorA");	
	TH1D* histoPCMMassPi0DatapPb = 			(TH1D*)directoryPCMPi0pPb->Get("MassPi0");
	TH1D* histoPCMMassPi0MCpPb = 			(TH1D*)directoryPCMPi0pPb->Get("TrueMassPi0");
	TH1D* histoPCMWidthPi0DatapPb = 		(TH1D*)directoryPCMPi0pPb->Get("FWHMPi0MeV");
	TH1D* histoPCMWidthPi0MCpPb = 			(TH1D*)directoryPCMPi0pPb->Get("TrueFWHMPi0MeV");
	TH1D* histoPi0InputMCWOWeightspPb = 		(TH1D*)directoryPCMPi0pPb->Get("Pi0MCInputWOWeights");
	TH1D* histoPi0InputMCWOWeightsAddedSigpPb = (TH1D*)directoryPCMPi0pPb->Get("Pi0MCInputWOWeightsAddedSig");
	TH1D* histoPi0InputMCpPb = 					(TH1D*)directoryPCMPi0pPb->Get("Pi0MCInput");
	TH1D* histoPi0InputMCAddedSigpPb = 			(TH1D*)directoryPCMPi0pPb->Get("Pi0MCInputAddedSig");
	TH1D* histoPi0EfficiencypPb = (TH1D*)directoryPCMPi0pPb->Get("EfficiencyPi0");
	TH1D* histoPi0AcceptancepPb = (TH1D*)directoryPCMPi0pPb->Get("AcceptancePi0");
	
	histoPCMMassPi0DatapPb->Scale(1000.);
	histoPCMMassPi0MCpPb->Scale(1000.);
	Double_t nPCMEventpPb = 	histoPCMNumberOfEventspPb->GetBinContent(1);
   
	/*cout << "0-20%" << endl;
	TDirectory*	directoryPCMPi0pPb0020 = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_0-20%"); 
	TH1D* histoPCMNumberOfEvents0020= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV0-20%");
	TH1D* histoPCMYieldPi0pPb0020 = 				(TH1D*)directoryPCMPi0pPb0020->Get(nameHistoPCM.Data());
	TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb0020= 	(TGraphAsymmErrors*)directoryPCMPi0pPb0020->Get(nameGraphPCM.Data());	
	TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb0020= 	(TGraphAsymmErrors*)directoryPCMPi0pPb0020->Get("Pi0SystErrorA");	
	Double_t nPCMEventpPb0020 = 							histoPCMNumberOfEvents0020->GetBinContent(1);
	TH1D* histoPCMMassPi0DatapPb0020 =        (TH1D*)directoryPCMPi0pPb0020->Get("MassPi0");
	TH1D* histoPCMMassPi0MCpPb0020 =          (TH1D*)directoryPCMPi0pPb0020->Get("TrueMassPi0");
	TH1D* histoPCMWidthPi0DatapPb0020 =       (TH1D*)directoryPCMPi0pPb0020->Get("FWHMPi0MeV");
	TH1D* histoPCMWidthPi0MCpPb0020 =         (TH1D*)directoryPCMPi0pPb0020->Get("TrueFWHMPi0MeV");
	TH1D* histoPi0InputMCWOWeightspPb0020 = 		(TH1D*)directoryPCMPi0pPb0020->Get("Pi0MCInputWOWeights");
	TH1D* histoPi0InputMCWOWeightsAddedSigpPb0020 = (TH1D*)directoryPCMPi0pPb0020->Get("Pi0MCInputWOWeightsAddedSig");
	TH1D* histoPi0InputMCpPb0020 = 					(TH1D*)directoryPCMPi0pPb0020->Get("Pi0MCInput");
	TH1D* histoPi0InputMCAddedSigpPb0020 = 			(TH1D*)directoryPCMPi0pPb0020->Get("Pi0MCInputAddedSig");
	TH1D* histoPi0EfficiencypPb0020 = (TH1D*)directoryPCMPi0pPb0020->Get("EfficiencyPi0");
	TH1D* histoPi0AcceptancepPb0020 = (TH1D*)directoryPCMPi0pPb0020->Get("AcceptancePi0");

	histoPCMMassPi0DatapPb0020->Scale(1000.);
	histoPCMMassPi0MCpPb0020->Scale(1000.);
   
	cout << "20-40%" << endl;
	TDirectory*	directoryPCMPi0pPb2040 = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_20-40%"); 
	TH1D* histoPCMNumberOfEvents2040= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV20-40%");
	TH1D* histoPCMYieldPi0pPb2040 = 				(TH1D*)directoryPCMPi0pPb2040->Get(nameHistoPCM.Data());
	TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb2040= 	(TGraphAsymmErrors*)directoryPCMPi0pPb2040->Get(nameGraphPCM.Data());	
	TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb2040= 	(TGraphAsymmErrors*)directoryPCMPi0pPb2040->Get("Pi0SystErrorA");		
	Double_t nPCMEventpPb2040 = 							histoPCMNumberOfEvents2040->GetBinContent(1);
	TH1D* histoPCMMassPi0DatapPb2040 =        (TH1D*)directoryPCMPi0pPb2040->Get("MassPi0");
	TH1D* histoPCMMassPi0MCpPb2040 =          (TH1D*)directoryPCMPi0pPb2040->Get("TrueMassPi0");
	TH1D* histoPCMWidthPi0DatapPb2040 =       (TH1D*)directoryPCMPi0pPb2040->Get("FWHMPi0MeV");
	TH1D* histoPCMWidthPi0MCpPb2040 =         (TH1D*)directoryPCMPi0pPb2040->Get("TrueFWHMPi0MeV");
	TH1D* histoPi0InputMCWOWeightspPb2040 = 		(TH1D*)directoryPCMPi0pPb2040->Get("Pi0MCInputWOWeights");
	TH1D* histoPi0InputMCWOWeightsAddedSigpPb2040 = (TH1D*)directoryPCMPi0pPb2040->Get("Pi0MCInputWOWeightsAddedSig");
	TH1D* histoPi0InputMCpPb2040 = 					(TH1D*)directoryPCMPi0pPb2040->Get("Pi0MCInput");
	TH1D* histoPi0InputMCAddedSigpPb2040 = 			(TH1D*)directoryPCMPi0pPb2040->Get("Pi0MCInputAddedSig");
	TH1D* histoPi0EfficiencypPb2040 = (TH1D*)directoryPCMPi0pPb2040->Get("EfficiencyPi0");
	TH1D* histoPi0AcceptancepPb2040 = (TH1D*)directoryPCMPi0pPb2040->Get("AcceptancePi0");

	
	histoPCMMassPi0DatapPb2040->Scale(1000.);
	histoPCMMassPi0MCpPb2040->Scale(1000.);
		
	cout << "40-60%" << endl;
	TDirectory*	directoryPCMPi0pPb4060 = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_40-60%"); 
	TH1D* histoPCMNumberOfEvents4060= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV40-60%");
	TH1D* histoPCMYieldPi0pPb4060 = 				(TH1D*)directoryPCMPi0pPb4060->Get(nameHistoPCM.Data());
	TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb4060= 	(TGraphAsymmErrors*)directoryPCMPi0pPb4060->Get(nameGraphPCM.Data());	
	TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb4060= 	(TGraphAsymmErrors*)directoryPCMPi0pPb4060->Get("Pi0SystErrorA");	
	Double_t nPCMEventpPb4060 = 							histoPCMNumberOfEvents4060->GetBinContent(1);
	TH1D* histoPCMMassPi0DatapPb4060 =        (TH1D*)directoryPCMPi0pPb4060->Get("MassPi0");
	TH1D* histoPCMMassPi0MCpPb4060 =          (TH1D*)directoryPCMPi0pPb4060->Get("TrueMassPi0");
	TH1D* histoPCMWidthPi0DatapPb4060 =       (TH1D*)directoryPCMPi0pPb4060->Get("FWHMPi0MeV");
	TH1D* histoPCMWidthPi0MCpPb4060 =         (TH1D*)directoryPCMPi0pPb4060->Get("TrueFWHMPi0MeV");
	TH1D* histoPi0InputMCWOWeightspPb4060 = 		(TH1D*)directoryPCMPi0pPb4060->Get("Pi0MCInputWOWeights");
	TH1D* histoPi0InputMCWOWeightsAddedSigpPb4060 = (TH1D*)directoryPCMPi0pPb4060->Get("Pi0MCInputWOWeightsAddedSig");
	TH1D* histoPi0InputMCpPb4060 = 					(TH1D*)directoryPCMPi0pPb4060->Get("Pi0MCInput");
	TH1D* histoPi0InputMCAddedSigpPb4060 = 			(TH1D*)directoryPCMPi0pPb4060->Get("Pi0MCInputAddedSig");
	TH1D* histoPi0EfficiencypPb4060 = (TH1D*)directoryPCMPi0pPb4060->Get("EfficiencyPi0");
	TH1D* histoPi0AcceptancepPb4060 = (TH1D*)directoryPCMPi0pPb4060->Get("AcceptancePi0");

	histoPCMMassPi0DatapPb4060->Scale(1000.);
	histoPCMMassPi0MCpPb4060->Scale(1000.);

	cout << "60-80%" << endl;
	TDirectory*	directoryPCMPi0pPb6080 = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_60-80%"); 
	TH1D* histoPCMNumberOfEvents6080= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV60-80%");
	TH1D* histoPCMYieldPi0pPb6080 = 				(TH1D*)directoryPCMPi0pPb6080->Get(nameHistoPCM.Data());
	TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb6080= 	(TGraphAsymmErrors*)directoryPCMPi0pPb6080->Get(nameGraphPCM.Data());	
	TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb6080= 	(TGraphAsymmErrors*)directoryPCMPi0pPb6080->Get("Pi0SystErrorA");	
	TH1D* histoPCMMassPi0DatapPb6080 = 			(TH1D*)directoryPCMPi0pPb6080->Get("MassPi0");
	TH1D* histoPCMMassPi0MCpPb6080 = 			(TH1D*)directoryPCMPi0pPb6080->Get("TrueMassPi0");
	TH1D* histoPCMWidthPi0DatapPb6080 = 		(TH1D*)directoryPCMPi0pPb6080->Get("FWHMPi0MeV");
	TH1D* histoPCMWidthPi0MCpPb6080 = 			(TH1D*)directoryPCMPi0pPb6080->Get("TrueFWHMPi0MeV");
	TH1D* histoPi0InputMCWOWeightspPb6080 = 		(TH1D*)directoryPCMPi0pPb6080->Get("Pi0MCInputWOWeights");
	TH1D* histoPi0InputMCWOWeightsAddedSigpPb6080 = (TH1D*)directoryPCMPi0pPb6080->Get("Pi0MCInputWOWeightsAddedSig");
	TH1D* histoPi0InputMCpPb6080 = 					(TH1D*)directoryPCMPi0pPb6080->Get("Pi0MCInput");
	TH1D* histoPi0InputMCAddedSigpPb6080 = 			(TH1D*)directoryPCMPi0pPb6080->Get("Pi0MCInputAddedSig");
	TH1D* histoPi0EfficiencypPb6080 = (TH1D*)directoryPCMPi0pPb6080->Get("EfficiencyPi0");
	TH1D* histoPi0AcceptancepPb6080 = (TH1D*)directoryPCMPi0pPb6080->Get("AcceptancePi0");

	histoPCMMassPi0DatapPb6080->Scale(1000.);
	histoPCMMassPi0MCpPb6080->Scale(1000.);
	Double_t nPCMEventpPb6080 = 							histoPCMNumberOfEvents6080->GetBinContent(1);
	
   
	cout << "60-100%" << endl;
	TDirectory*	directoryPCMPi0pPb60100 =            (TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_60-0%"); 
	TH1D* histoPCMNumberOfEvents60100=         (TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV60-0%");
	TH1D* histoPCMYieldPi0pPb60100 =           (TH1D*)directoryPCMPi0pPb60100->Get(nameHistoPCM.Data());
	TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb60100=   (TGraphAsymmErrors*)directoryPCMPi0pPb60100->Get(nameGraphPCM.Data());   
	TGraphAsymmErrors* graphPCMYieldPi0SysErrRAApPb60100=   (TGraphAsymmErrors*)directoryPCMPi0pPb60100->Get("Pi0SystErrorA"); 
	TH1D* histoPCMMassPi0DatapPb60100 =        (TH1D*)directoryPCMPi0pPb60100->Get("MassPi0");
	TH1D* histoPCMMassPi0MCpPb60100 =          (TH1D*)directoryPCMPi0pPb60100->Get("TrueMassPi0");
	TH1D* histoPCMWidthPi0DatapPb60100 =       (TH1D*)directoryPCMPi0pPb60100->Get("FWHMPi0MeV");
	TH1D* histoPCMWidthPi0MCpPb60100 =         (TH1D*)directoryPCMPi0pPb60100->Get("TrueFWHMPi0MeV");
	TH1D* histoPi0InputMCWOWeightspPb60100 = 		(TH1D*)directoryPCMPi0pPb60100->Get("Pi0MCInputWOWeights");
	TH1D* histoPi0InputMCWOWeightsAddedSigpPb60100 = (TH1D*)directoryPCMPi0pPb60100->Get("Pi0MCInputWOWeightsAddedSig");
	TH1D* histoPi0InputMCpPb60100 = 					(TH1D*)directoryPCMPi0pPb60100->Get("Pi0MCInput");
	TH1D* histoPi0InputMCAddedSigpPb60100 = 			(TH1D*)directoryPCMPi0pPb60100->Get("Pi0MCInputAddedSig");
	TH1D* histoPi0EfficiencypPb60100 = (TH1D*)directoryPCMPi0pPb60100->Get("EfficiencyPi0");
	TH1D* histoPi0AcceptancepPb60100 = (TH1D*)directoryPCMPi0pPb60100->Get("AcceptancePi0");

	histoPCMMassPi0DatapPb60100->Scale(1000.);
	histoPCMMassPi0MCpPb60100->Scale(1000.);
	Double_t nPCMEventpPb60100 =                     histoPCMNumberOfEvents60100->GetBinContent(1);	
	*///****************************************************************************************************
	//************************** Read Eta data for PCM *******************************************************
	//****************************************************************************************************
	cout << "PCM Eta" << endl;
	cout << "0-100%" << endl;
	TDirectory*	directoryPCMEtapPb = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_0-100%"); 
	TH1D* histoPCMYieldEtapPb = 				(TH1D*)directoryPCMEtapPb->Get(nameHistoPCMEta.Data());
	TH1D* histoPCMRAWYieldEtapPb =            (TH1D*)directoryPCMEtapPb->Get("RAWYieldPerEventsEta");
	TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb= 	(TGraphAsymmErrors*)directoryPCMEtapPb->Get(nameGraphPCMEta.Data());	
	TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb= 	(TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtaSystErrorA");	
	TH1D* histoPCMEtaPi0RatiopPb =            (TH1D*)directoryPCMEtapPb->Get("EtatoPi0Ratio");
	TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb=    (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtatoPi0RatioSys"); 
	TH1D* histoPCMMassEtaDatapPb = 			(TH1D*)directoryPCMEtapPb->Get("MassEta");
	TH1D* histoPCMMassEtaMCpPb = 			(TH1D*)directoryPCMEtapPb->Get("TrueMassEta");
	TH1D* histoPCMWidthEtaDatapPb = 		(TH1D*)directoryPCMEtapPb->Get("FWHMEtaMeV");
	TH1D* histoPCMWidthEtaMCpPb = 			(TH1D*)directoryPCMEtapPb->Get("TrueFWHMEtaMeV");
	TH1D* histoEtaInputMCWOWeightspPb = 			(TH1D*)directoryPCMEtapPb->Get("EtaMCInputWOWeights");
	TH1D* histoEtaInputMCWOWeightsAddedSigpPb = 			(TH1D*)directoryPCMEtapPb->Get("EtaMCInputWOWeightsAddedSig");
	TH1D* histoEtaInputMCpPb = 			(TH1D*)directoryPCMEtapPb->Get("EtaMCInput");
	TH1D* histoEtaInputMCAddedSigpPb = 			(TH1D*)directoryPCMEtapPb->Get("EtaMCInputAddedSig");
	TH1D* histoEtaEfficiencypPb = (TH1D*)directoryPCMEtapPb->Get("EfficiencyEta");
	TH1D* histoEtaAcceptancepPb = (TH1D*)directoryPCMEtapPb->Get("AcceptanceEta");

	histoPCMMassEtaDatapPb->Scale(1000.);
	histoPCMMassEtaMCpPb->Scale(1000.);

// 	cout << "0-20%" << endl;
// 	TDirectory*	directoryPCMEtapPb0020 = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_0-20%"); 
// 	TH1D* histoPCMYieldEtapPb0020 = 				(TH1D*)directoryPCMEtapPb0020->Get(nameHistoPCMEta.Data());
// 	TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb0020= 	(TGraphAsymmErrors*)directoryPCMEtapPb0020->Get(nameGraphPCMEta.Data());	
// 	TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb0020= 	(TGraphAsymmErrors*)directoryPCMEtapPb0020->Get("EtaSystErrorA");
// 	TH1D* histoPCMEtaPi0RatiopPb0020 =            (TH1D*)directoryPCMEtapPb0020->Get("EtatoPi0Ratio");
// 	TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb0020=    (TGraphAsymmErrors*)directoryPCMEtapPb0020->Get("EtatoPi0RatioSys"); 
// 	TH1D* histoPCMMassEtaDatapPb0020 =        (TH1D*)directoryPCMEtapPb0020->Get("MassEta");
// 	TH1D* histoPCMMassEtaMCpPb0020 =          (TH1D*)directoryPCMEtapPb0020->Get("TrueMassEta");
// 	TH1D* histoPCMWidthEtaDatapPb0020 =       (TH1D*)directoryPCMEtapPb0020->Get("FWHMEtaMeV");
// 	TH1D* histoPCMWidthEtaMCpPb0020 =         (TH1D*)directoryPCMEtapPb0020->Get("TrueFWHMEtaMeV");
// 	TH1D* histoEtaInputMCWOWeightspPb0020 = 		(TH1D*)directoryPCMEtapPb0020->Get("EtaMCInputWOWeights");
// 	TH1D* histoEtaInputMCWOWeightsAddedSigpPb0020 = (TH1D*)directoryPCMEtapPb0020->Get("EtaMCInputWOWeightsAddedSig");
// 	TH1D* histoEtaInputMCpPb0020 = 					(TH1D*)directoryPCMEtapPb0020->Get("EtaMCInput");
// 	TH1D* histoEtaInputMCAddedSigpPb0020 = 			(TH1D*)directoryPCMEtapPb0020->Get("EtaMCInputAddedSig");
// 	TH1D* histoEtaEfficiencypPb0020 = (TH1D*)directoryPCMEtapPb0020->Get("EfficiencyEta");
// 	TH1D* histoEtaAcceptancepPb0020 = (TH1D*)directoryPCMEtapPb0020->Get("AcceptanceEta");
// 
// 	histoPCMMassEtaDatapPb0020->Scale(1000.);
// 	histoPCMMassEtaMCpPb0020->Scale(1000.);
// 
// 	cout << "20-40%" << endl;
// 	TDirectory*	directoryPCMEtapPb2040 = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_20-40%"); 
// 	TH1D* histoPCMYieldEtapPb2040 = 				(TH1D*)directoryPCMEtapPb2040->Get(nameHistoPCMEta.Data());
// 	TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb2040= 	(TGraphAsymmErrors*)directoryPCMEtapPb2040->Get(nameGraphPCMEta.Data());	
// 	TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb2040= 	(TGraphAsymmErrors*)directoryPCMEtapPb2040->Get("EtaSystErrorA");		
// 	TH1D* histoPCMEtaPi0RatiopPb2040 =            (TH1D*)directoryPCMEtapPb2040->Get("EtatoPi0Ratio");
// 	TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb2040=    (TGraphAsymmErrors*)directoryPCMEtapPb2040->Get("EtatoPi0RatioSys"); 
// 	TH1D* histoPCMMassEtaDatapPb2040 =        (TH1D*)directoryPCMEtapPb2040->Get("MassEta");
// 	TH1D* histoPCMMassEtaMCpPb2040 =          (TH1D*)directoryPCMEtapPb2040->Get("TrueMassEta");
// 	TH1D* histoPCMWidthEtaDatapPb2040 =       (TH1D*)directoryPCMEtapPb2040->Get("FWHMEtaMeV");
// 	TH1D* histoPCMWidthEtaMCpPb2040 =         (TH1D*)directoryPCMEtapPb2040->Get("TrueFWHMEtaMeV");
// 	TH1D* histoEtaInputMCWOWeightspPb2040 = 		(TH1D*)directoryPCMEtapPb2040->Get("EtaMCInputWOWeights");
// 	TH1D* histoEtaInputMCWOWeightsAddedSigpPb2040 = (TH1D*)directoryPCMEtapPb2040->Get("EtaMCInputWOWeightsAddedSig");
// 	TH1D* histoEtaInputMCpPb2040 = 					(TH1D*)directoryPCMEtapPb2040->Get("EtaMCInput");
// 	TH1D* histoEtaInputMCAddedSigpPb2040 = 			(TH1D*)directoryPCMEtapPb2040->Get("EtaMCInputAddedSig");
// 	TH1D* histoEtaEfficiencypPb2040 = (TH1D*)directoryPCMEtapPb2040->Get("EfficiencyEta");
// 	TH1D* histoEtaAcceptancepPb2040 = (TH1D*)directoryPCMEtapPb2040->Get("AcceptanceEta");
// 	histoPCMMassEtaDatapPb2040->Scale(1000.);
// 	histoPCMMassEtaMCpPb2040->Scale(1000.);
// 		
// 	cout << "40-60%" << endl;
// 	TDirectory*	directoryPCMEtapPb4060 = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_40-60%"); 
// 	TH1D* histoPCMYieldEtapPb4060 = 				(TH1D*)directoryPCMEtapPb4060->Get(nameHistoPCMEta.Data());
// 	TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb4060= 	(TGraphAsymmErrors*)directoryPCMEtapPb4060->Get(nameGraphPCMEta.Data());	
// 	TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb4060= 	(TGraphAsymmErrors*)directoryPCMEtapPb4060->Get("EtaSystErrorA");	
// 	TH1D* histoPCMEtaPi0RatiopPb4060 =            (TH1D*)directoryPCMEtapPb4060->Get("EtatoPi0Ratio");
// 	TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb4060=    (TGraphAsymmErrors*)directoryPCMEtapPb4060->Get("EtatoPi0RatioSys"); 
// 	TH1D* histoPCMMassEtaDatapPb4060 =        (TH1D*)directoryPCMEtapPb4060->Get("MassEta");
// 	TH1D* histoPCMMassEtaMCpPb4060 =          (TH1D*)directoryPCMEtapPb4060->Get("TrueMassEta");
// 	TH1D* histoPCMWidthEtaDatapPb4060 =       (TH1D*)directoryPCMEtapPb4060->Get("FWHMEtaMeV");
// 	TH1D* histoPCMWidthEtaMCpPb4060 =         (TH1D*)directoryPCMEtapPb4060->Get("TrueFWHMEtaMeV");
// 	TH1D* histoEtaInputMCWOWeightspPb4060 = 		(TH1D*)directoryPCMEtapPb4060->Get("EtaMCInputWOWeights");
// 	TH1D* histoEtaInputMCWOWeightsAddedSigpPb4060 = (TH1D*)directoryPCMEtapPb4060->Get("EtaMCInputWOWeightsAddedSig");
// 	TH1D* histoEtaInputMCpPb4060 = 					(TH1D*)directoryPCMEtapPb4060->Get("EtaMCInput");
// 	TH1D* histoEtaInputMCAddedSigpPb4060 = 			(TH1D*)directoryPCMEtapPb4060->Get("EtaMCInputAddedSig");
// 	TH1D* histoEtaEfficiencypPb4060 = (TH1D*)directoryPCMEtapPb4060->Get("EfficiencyEta");
// 	TH1D* histoEtaAcceptancepPb4060 = (TH1D*)directoryPCMEtapPb4060->Get("AcceptanceEta");
// 
// 	histoPCMMassEtaDatapPb4060->Scale(1000.);
// 	histoPCMMassEtaMCpPb4060->Scale(1000.);
// 
// 	cout << "60-80%" << endl;
// 	TDirectory*	directoryPCMEtapPb6080 = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_60-80%"); 
// 	TH1D* histoPCMYieldEtapPb6080 = 				(TH1D*)directoryPCMEtapPb6080->Get(nameHistoPCMEta.Data());
// 	TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb6080= 	(TGraphAsymmErrors*)directoryPCMEtapPb6080->Get(nameGraphPCMEta.Data());	
// 	TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb6080= 	(TGraphAsymmErrors*)directoryPCMEtapPb6080->Get("EtaSystErrorA");	
// 	TH1D* histoPCMEtaPi0RatiopPb6080 =            (TH1D*)directoryPCMEtapPb6080->Get("EtatoPi0Ratio");
// 	TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb6080=    (TGraphAsymmErrors*)directoryPCMEtapPb6080->Get("EtatoPi0RatioSys"); 
// 	TH1D* histoPCMMassEtaDatapPb6080 = 			(TH1D*)directoryPCMEtapPb6080->Get("MassEta");
// 	TH1D* histoPCMMassEtaMCpPb6080 = 			(TH1D*)directoryPCMEtapPb6080->Get("TrueMassEta");
// 	TH1D* histoPCMWidthEtaDatapPb6080 = 		(TH1D*)directoryPCMEtapPb6080->Get("FWHMEtaMeV");
// 	TH1D* histoPCMWidthEtaMCpPb6080 = 			(TH1D*)directoryPCMEtapPb6080->Get("TrueFWHMEtaMeV");
// 	TH1D* histoEtaInputMCWOWeightspPb6080 = 		(TH1D*)directoryPCMEtapPb6080->Get("EtaMCInputWOWeights");
// 	TH1D* histoEtaInputMCWOWeightsAddedSigpPb6080 = (TH1D*)directoryPCMEtapPb6080->Get("EtaMCInputWOWeightsAddedSig");
// 	TH1D* histoEtaInputMCpPb6080 = 					(TH1D*)directoryPCMEtapPb6080->Get("EtaMCInput");
// 	TH1D* histoEtaInputMCAddedSigpPb6080 = 			(TH1D*)directoryPCMEtapPb6080->Get("EtaMCInputAddedSig");
// 	TH1D* histoEtaEfficiencypPb6080 = (TH1D*)directoryPCMEtapPb6080->Get("EfficiencyEta");
// 	TH1D* histoEtaAcceptancepPb6080 = (TH1D*)directoryPCMEtapPb6080->Get("AcceptanceEta");
// 	histoPCMMassEtaDatapPb6080->Scale(1000.);
// 	histoPCMMassEtaMCpPb6080->Scale(1000.);
// 
// 	cout << "60-100%" << endl;
// 	TDirectory*	directoryPCMEtapPb60100 = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_60-0%"); 
// 	TH1D* histoPCMYieldEtapPb60100 = 				(TH1D*)directoryPCMEtapPb60100->Get(nameHistoPCMEta.Data());
// 	TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb60100= 	(TGraphAsymmErrors*)directoryPCMEtapPb60100->Get(nameGraphPCMEta.Data());	
// 	TGraphAsymmErrors* graphPCMYieldEtaSysErrRAApPb60100= 	(TGraphAsymmErrors*)directoryPCMEtapPb60100->Get("EtaSystErrorA");	
// 	TH1D* histoPCMEtaPi0RatiopPb60100 =            (TH1D*)directoryPCMEtapPb60100->Get("EtatoPi0Ratio");
// 	TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb60100=    (TGraphAsymmErrors*)directoryPCMEtapPb60100->Get("EtatoPi0RatioSys"); 
// 	TH1D* histoPCMMassEtaDatapPb60100 = 			(TH1D*)directoryPCMEtapPb60100->Get("MassEta");
// 	TH1D* histoPCMMassEtaMCpPb60100 = 			(TH1D*)directoryPCMEtapPb60100->Get("TrueMassEta");
// 	TH1D* histoPCMWidthEtaDatapPb60100 = 		(TH1D*)directoryPCMEtapPb60100->Get("FWHMEtaMeV");
// 	TH1D* histoPCMWidthEtaMCpPb60100 = 			(TH1D*)directoryPCMEtapPb60100->Get("TrueFWHMEtaMeV");
// 	TH1D* histoEtaInputMCWOWeightspPb60100 = 		(TH1D*)directoryPCMEtapPb60100->Get("EtaMCInputWOWeights");
// 	TH1D* histoEtaInputMCWOWeightsAddedSigpPb60100 = (TH1D*)directoryPCMEtapPb60100->Get("EtaMCInputWOWeightsAddedSig");
// 	TH1D* histoEtaInputMCpPb60100 = 					(TH1D*)directoryPCMEtapPb60100->Get("EtaMCInput");
// 	TH1D* histoEtaInputMCAddedSigpPb60100 = 			(TH1D*)directoryPCMEtapPb60100->Get("EtaMCInputAddedSig");
// 	TH1D* histoEtaEfficiencypPb60100 = (TH1D*)directoryPCMEtapPb60100->Get("EfficiencyEta");
// 	TH1D* histoEtaAcceptancepPb60100 = (TH1D*)directoryPCMEtapPb60100->Get("AcceptanceEta");
// 	histoPCMMassEtaDatapPb60100->Scale(1000.);
// 	histoPCMMassEtaMCpPb60100->Scale(1000.);

	TFile* fileGammaPCMpPb = new TFile("8000011002092170008260400000_01621035009000/pPb_5.023TeV/Gamma_Pi0_data_GammaConvV1Correction_8000011002092170008260400000_01621035009000.root");
	TH1D* histoGammaPuritypPb = (TH1D*)fileGammaPCMpPb->Get("MCGammaTruePurity");
	TH1D* histoGammaRecoEffpPb = (TH1D*)fileGammaPCMpPb->Get("MCGammaPrimaryRecoEff");
	
	TFile* fileGammaPCMpp7TeV = new TFile("900366208010033211360000000900/7TeV/Pi0_MC_AnalysisResultsCorrectionHistos_900366208010033211360000000900.root");
	TH1D* histoGammaConvProbpp7TeV = (TH1D*)fileGammaPCMpp7TeV->Get("fMCGammaConvProb");
	
	TFile* fileGammaPCMpPbOnlyCorr = new TFile("8000011002092170008260400000_01621035009000/pPb_5.023TeV/Pi0_MC_GammaConvV1CorrectionHistos_8000011002092170008260400000_01621035009000.root");
	TH1D* histoGammaConvProbpPb = (TH1D*)fileGammaPCMpPbOnlyCorr->Get("MCGammaConvProb");
	
// 	TFile* filePi0PCMpPb_DPMJET = 					new TFile("DPMJET/8000011002092170008260400000_01621035009000/pPb_5.023TeV/Pi0_data_GammaConvV1Correction_8000011002092170008260400000_01621035009000.root");
// 	TH1D* histoPi0InputMCWOWeightspPb_DPMJET = 	(TH1D*)filePi0PCMpPb_DPMJET->Get("MCYield_Meson_oldBinWOWeights");
// 	TH1D* histoPi0InputMCpP_DPMJET = 			(TH1D*)filePi0PCMpPb_DPMJET->Get("MCYield_Meson_oldBin");
// 	TH1D* histoPi0EfficiencypPb_DPMJET = 		(TH1D*)filePi0PCMpPb_DPMJET->Get("TrueMesonEffiPt");
// 	TH1D* histoPi0AcceptancepPb_DPMJET = 		(TH1D*)filePi0PCMpPb_DPMJET->Get("fMCMesonAccepPt");
// 
// 	TFile* fileEtaPCMpPb_DPMJET = 					new TFile("DPMJET/8000011002092170008260400000_01621035009000/pPb_5.023TeV/Eta_data_GammaConvV1Correction_8000011002092170008260400000_01621035009000.root");
// 	TH1D* histoEtaInputMCWOWeightspPb_DPMJET = 	(TH1D*)fileEtaPCMpPb_DPMJET->Get("MCYield_Meson_oldBinWOWeights");
// 	TH1D* histoEtaInputMCpP_DPMJET = 			(TH1D*)fileEtaPCMpPb_DPMJET->Get("MCYield_Meson_oldBin");
// 	TH1D* histoEtaEfficiencypPb_DPMJET = 		(TH1D*)fileEtaPCMpPb_DPMJET->Get("TrueMesonEffiPt");
// 	TH1D* histoEtaAcceptancepPb_DPMJET = 		(TH1D*)fileEtaPCMpPb_DPMJET->Get("fMCMesonAccepPt");
	
// 	TFile* fileGammaPCMpPb_DPMJET = new TFile("DPMJET/8000011002092170008260400000_01621035009000/pPb_5.023TeV/Gamma_Pi0_data_GammaConvV1Correction_8000011002092170008260400000_01621035009000.root");
// 	TH1D* histoGammaPuritypPb_DPMJET = (TH1D*)fileGammaPCMpPb_DPMJET->Get("MCGammaTruePurity");
// 	TH1D* histoGammaRecoEffpPb_DPMJET = (TH1D*)fileGammaPCMpPb_DPMJET->Get("MCGammaPrimaryRecoEff");
// 	TFile* fileGammaPCMpPbOnlyCorr_DPMJET = new TFile("DPMJET/8000011002092170008260400000_01621035009000/pPb_5.023TeV/Pi0_MC_GammaConvV1CorrectionHistos_8000011002092170008260400000_01621035009000.root");
// 	TH1D* histoGammaConvProbpPb_DPMJET = (TH1D*)fileGammaPCMpPbOnlyCorr_DPMJET->Get("MCGammaConvProb");

	
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

		
	//	**********************************************************************************************************************
	//	****************************************Pi0 Spectra compared to MC****************************************************
	//	**********************************************************************************************************************

	
// 	TCanvas* canvasPi0SpectraAllTogether = new TCanvas("canvasPi0SpectraAllTogether","",200,10,1200,1100);  // gives the page size
// 	DrawGammaCanvasSettings( canvasPi0SpectraAllTogether,  0.13, 0.01, 0.015, 0.08);
// 	
// 	canvasPi0SpectraAllTogether->SetLogy();
// 	canvasPi0SpectraAllTogether->SetLogx();
// 	TH2F * histo2DPi0SpectraAll;
// 	histo2DPi0SpectraAll = new TH2F("histo2DPi0SpectraAll","histo2DPi0SpectraAll",1000,0.23,20.,1000,1e-8,2e2 );
// 	SetStyleHistoTH2ForGraphs(histo2DPi0SpectraAll, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{d#it{N}^{AA}}{d#it{p}_{T} d#it{y}}}{ #frac{1}{N_{evt}^{pp}}#frac{d#it{N}^{pp}}{d#it{p}_{T} d#it{y}}}
// 	histo2DPi0SpectraAll->GetXaxis()->SetLabelOffset(-0.01);
// 	histo2DPi0SpectraAll->GetYaxis()->SetLabelOffset(0.01);
// 	histo2DPi0SpectraAll->DrawCopy(); 
// 	
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrpPb, markerStylepPb,markerSizepPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE);
// 	graphPCMYieldPi0SysErrpPb->Draw("E2same");
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrpPb0020, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020, widthLinesBoxes, kTRUE);
// 	graphPCMYieldPi0SysErrpPb0020->Draw("E2same");
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrpPb2040, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040, widthLinesBoxes, kTRUE);
// 	graphPCMYieldPi0SysErrpPb2040->Draw("E2same");
// 	
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrpPb4060, markerStylepPb4060,markerSizepPb4060, colorCombpPb4060 , colorCombpPb4060, widthLinesBoxes, kTRUE);
// 	graphPCMYieldPi0SysErrpPb4060->Draw("E2same");
// 	
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrpPb6080, markerStylepPb6080,markerSizepPb6080, colorCombpPb6080 , colorCombpPb6080, widthLinesBoxes, kTRUE);
// 	graphPCMYieldPi0SysErrpPb6080->Draw("E2same");
// 	
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrpPb60100, markerStylepPb60100,markerSizepPb60100, colorCombpPb60100 , colorCombpPb60100, widthLinesBoxes, kTRUE);
// 	graphPCMYieldPi0SysErrpPb60100->Draw("E2same");
// 	
// 	DrawGammaSetMarker(histoPCMYieldPi0pPb, markerStylepPb,markerSizepPb, colorCombpPb , colorCombpPb);
// 	histoPCMYieldPi0pPb->Draw("p,same,e1");
// 	DrawGammaSetMarker(histoPCMYieldPi0pPb0020, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020);
// 	histoPCMYieldPi0pPb0020->Draw("p,same,e1");
// 	
// 	DrawGammaSetMarker(histoPCMYieldPi0pPb2040, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	histoPCMYieldPi0pPb2040->Draw("p,same,e1");
// 	
// 	DrawGammaSetMarker(histoPCMYieldPi0pPb4060, markerStylepPb4060,markerSizepPb4060, colorCombpPb4060 , colorCombpPb4060);
// 	histoPCMYieldPi0pPb4060->Draw("p,same,e1");
// 	
// 	DrawGammaSetMarker(histoPCMYieldPi0pPb6080, markerStylepPb6080,markerSizepPb6080, colorCombpPb6080 , colorCombpPb6080);
// 	histoPCMYieldPi0pPb6080->Draw("p,same,e1");
// 	
// 	DrawGammaSetMarker(histoPCMYieldPi0pPb60100, markerStylepPb60100,markerSizepPb60100, colorCombpPb60100 , colorCombpPb60100);
// 	histoPCMYieldPi0pPb60100->Draw("same");
// 
// 	TLatex *labelSpectraPi0Label = new TLatex(0.55,0.92,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
// 	SetStyleTLatex( labelSpectraPi0Label, 0.035  ,4);
// 	labelSpectraPi0Label->Draw();
// 	
// 	TLegend* legendSpectraPi0 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraPi0->SetFillColor(0);
// 	legendSpectraPi0->SetLineColor(0);
// 	legendSpectraPi0->SetTextSize(0.025);
// 	legendSpectraPi0->SetNColumns(2);
// 	legendSpectraPi0->SetMargin(0.2);
// 	legendSpectraPi0->AddEntry(graphPCMYieldPi0SysErrpPb,"PCM","pf");
// 	legendSpectraPi0->AddEntry((TObject*)0, collisionSystempPb.Data(),"");
// 	legendSpectraPi0->AddEntry(graphPCMYieldPi0SysErrpPb0020,"PCM","pf");
// 	legendSpectraPi0->AddEntry((TObject*)0, collisionSystempPb0020.Data(),"");
// 	legendSpectraPi0->AddEntry(graphPCMYieldPi0SysErrpPb2040,"PCM","pf");
// 	legendSpectraPi0->AddEntry((TObject*)0, collisionSystempPb2040.Data(),"");
// 	legendSpectraPi0->AddEntry(graphPCMYieldPi0SysErrpPb4060,"PCM","pf");
// 	legendSpectraPi0->AddEntry((TObject*)0, collisionSystempPb4060.Data(),"");
// 	legendSpectraPi0->AddEntry(graphPCMYieldPi0SysErrpPb6080,"PCM","pf");
// 	legendSpectraPi0->AddEntry((TObject*)0, collisionSystempPb6080.Data(),"");
// 	legendSpectraPi0->AddEntry(graphPCMYieldPi0SysErrpPb60100,"PCM","pf");
// 	legendSpectraPi0->AddEntry((TObject*)0, collisionSystempPb60100.Data(),"");
// 	legendSpectraPi0->Draw();
// 	
// 	canvasPi0SpectraAllTogether->Update();
// 	canvasPi0SpectraAllTogether->Print(Form("%s/Pi0_Spectra_All.%s",outputDir.Data(),suffix.Data()));
// 
// 	canvasPi0SpectraAllTogether->cd();
// 	histo2DPi0SpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldPi0SysErrpPb->Draw("E2same");
// 	histoPCMYieldPi0pPb->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoPi0InputMCWOWeightspPb, 1., lineStyleMCpPb, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightspPb->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCWOWeightsAddedSigpPb, 1., lineStyleMCAddSigpPb, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightsAddedSigpPb->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCpPb, 1., lineStyleMCpPb, colorCombMCReweightedpPb);  
// 	histoPi0InputMCpPb->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCAddedSigpPb, 1., lineStyleMCAddSigpPb, colorCombMCReweightedAddSigpPb);  
// 	histoPi0InputMCAddedSigpPb->Draw("same,hist,c");    
// 
// 	labelSpectraPi0Label->Draw();
// 	TLatex *labelSpectraMinBias = new TLatex(0.55,0.88,collisionSystempPb.Data());
// 	SetStyleTLatex( labelSpectraMinBias, 0.035  ,4);
// 	labelSpectraMinBias->Draw(); 
// 	
// 	TLegend* legendSpectraPi0MC = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraPi0MC->SetFillColor(0);
// 	legendSpectraPi0MC->SetLineColor(0);
// 	legendSpectraPi0MC->SetTextSize(0.025);
// 	legendSpectraPi0MC->SetNColumns(1);
// 	legendSpectraPi0MC->SetMargin(0.2);
// 	legendSpectraPi0MC->AddEntry(graphPCMYieldPi0SysErrpPb,"measured","pf");
// 	legendSpectraPi0MC->AddEntry(histoPi0InputMCWOWeightspPb,"HIJING","l");
// 	legendSpectraPi0MC->AddEntry(histoPi0InputMCWOWeightsAddedSigpPb,"added Signals","l");
// 	legendSpectraPi0MC->AddEntry(histoPi0InputMCpPb,"HIJING reweighted","l");
// 	legendSpectraPi0MC->AddEntry(histoPi0InputMCAddedSigpPb,"added Signals reweighted","l");
// 	
// 	legendSpectraPi0MC->Draw();
// 	
// 	canvasPi0SpectraAllTogether->Update();
// 	canvasPi0SpectraAllTogether->Print(Form("%s/Pi0_Spectra_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 	canvasPi0SpectraAllTogether->cd();
// 	histo2DPi0SpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldPi0SysErrpPb->Draw("E2same");
// 	histoPCMYieldPi0pPb->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoPi0InputMCWOWeightspPb, 1., lineStyleMCpPb, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightspPb->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCWOWeightsAddedSigpPb, 1., lineStyleMCAddSigpPb, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightsAddedSigpPb->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCWOWeightspPb_DPMJET, 1., lineStyleMCpPb, kRed-7);  
// 	histoPi0InputMCWOWeightspPb_DPMJET->Draw("same,hist,c");    
// 	
// 	labelSpectraPi0Label->Draw();
// 	labelSpectraMinBias->Draw(); 
// 	
// 	TLegend* legendSpectraPi0MCAll = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraPi0MCAll->SetFillColor(0);
// 	legendSpectraPi0MCAll->SetLineColor(0);
// 	legendSpectraPi0MCAll->SetTextSize(0.025);
// 	legendSpectraPi0MCAll->SetNColumns(1);
// 	legendSpectraPi0MCAll->SetMargin(0.2);
// 	legendSpectraPi0MCAll->AddEntry(graphPCMYieldPi0SysErrpPb,"measured","pf");
// 	legendSpectraPi0MCAll->AddEntry(histoPi0InputMCWOWeightspPb,"HIJING","l");
// 	legendSpectraPi0MCAll->AddEntry(histoPi0InputMCWOWeightsAddedSigpPb,"added Signals","l");
// 	legendSpectraPi0MCAll->AddEntry(histoPi0InputMCWOWeightspPb_DPMJET,"DPMJET","l");
// 	legendSpectraPi0MCAll->Draw();
// 	
// 	canvasPi0SpectraAllTogether->Update();
// 	canvasPi0SpectraAllTogether->Print(Form("%s/Pi0_Spectra_WithMC_HIJINGAndDPMJET.%s",outputDir.Data(),suffix.Data()));
// 
// 	canvasPi0SpectraAllTogether->cd();
// 	histo2DPi0SpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldPi0SysErrpPb0020->Draw("E2same");
// 	histoPCMYieldPi0pPb0020->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoPi0InputMCWOWeightspPb0020, 1., lineStyleMCpPb0020, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightspPb0020->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCWOWeightsAddedSigpPb0020, 1., lineStyleMCAddSigpPb0020, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightsAddedSigpPb0020->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCpPb0020, 1., lineStyleMCpPb0020, colorCombMCReweightedpPb);  
// 	histoPi0InputMCpPb0020->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCAddedSigpPb0020, 1., lineStyleMCAddSigpPb0020, colorCombMCReweightedAddSigpPb);  
// 	histoPi0InputMCAddedSigpPb0020->Draw("same,hist,c");    
// 
// 	labelSpectraPi0Label->Draw();
// 	TLatex *labelSpectra0020 = new TLatex(0.55,0.88,collisionSystempPb0020.Data());
// 	SetStyleTLatex( labelSpectra0020, 0.035  ,4);
// 	labelSpectra0020->Draw(); 
// 	
// 	TLegend* legendSpectraPi0MC0020 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraPi0MC0020->SetFillColor(0);
// 	legendSpectraPi0MC0020->SetLineColor(0);
// 	legendSpectraPi0MC0020->SetTextSize(0.025);
// 	legendSpectraPi0MC0020->SetNColumns(1);
// 	legendSpectraPi0MC0020->SetMargin(0.2);
// 	legendSpectraPi0MC0020->AddEntry(graphPCMYieldPi0SysErrpPb0020,"measured","pf");
// 	legendSpectraPi0MC0020->AddEntry(histoPi0InputMCWOWeightspPb0020,"HIJING","l");
// 	legendSpectraPi0MC0020->AddEntry(histoPi0InputMCWOWeightsAddedSigpPb0020,"added Signals","l");
// 	legendSpectraPi0MC0020->AddEntry(histoPi0InputMCpPb0020,"HIJING reweighted","l");
// 	legendSpectraPi0MC0020->AddEntry(histoPi0InputMCAddedSigpPb0020,"added Signals reweighted","l");
// 	
// 	legendSpectraPi0MC0020->Draw();
// 	
// 	canvasPi0SpectraAllTogether->Update();
// 	canvasPi0SpectraAllTogether->Print(Form("%s/Pi0_Spectra_pPb0020_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 	
// 	canvasPi0SpectraAllTogether->cd();
// 	histo2DPi0SpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldPi0SysErrpPb2040->Draw("E2same");
// 	histoPCMYieldPi0pPb2040->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoPi0InputMCWOWeightspPb2040, 1., lineStyleMCpPb2040, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightspPb2040->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCWOWeightsAddedSigpPb2040, 1., lineStyleMCAddSigpPb2040, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightsAddedSigpPb2040->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCpPb2040, 1., lineStyleMCpPb2040, colorCombMCReweightedpPb);  
// 	histoPi0InputMCpPb2040->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCAddedSigpPb2040, 1., lineStyleMCAddSigpPb2040, colorCombMCReweightedAddSigpPb);  
// 	histoPi0InputMCAddedSigpPb2040->Draw("same,hist,c");    
// 
// 	labelSpectraPi0Label->Draw();
// 	TLatex *labelSpectra2040 = new TLatex(0.55,0.88,collisionSystempPb2040.Data());
// 	SetStyleTLatex( labelSpectra2040, 0.035  ,4);
// 	labelSpectra2040->Draw(); 
// 	
// 	TLegend* legendSpectraPi0MC2040 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraPi0MC2040->SetFillColor(0);
// 	legendSpectraPi0MC2040->SetLineColor(0);
// 	legendSpectraPi0MC2040->SetTextSize(0.025);
// 	legendSpectraPi0MC2040->SetNColumns(1);
// 	legendSpectraPi0MC2040->SetMargin(0.2);
// 	legendSpectraPi0MC2040->AddEntry(graphPCMYieldPi0SysErrpPb2040,"measured","pf");
// 	legendSpectraPi0MC2040->AddEntry(histoPi0InputMCWOWeightspPb2040,"HIJING","l");
// 	legendSpectraPi0MC2040->AddEntry(histoPi0InputMCWOWeightsAddedSigpPb2040,"added Signals","l");
// 	legendSpectraPi0MC2040->AddEntry(histoPi0InputMCpPb2040,"HIJING reweighted","l");
// 	legendSpectraPi0MC2040->AddEntry(histoPi0InputMCAddedSigpPb2040,"added Signals reweighted","l");
// 	
// 	legendSpectraPi0MC2040->Draw();
// 	
// 	canvasPi0SpectraAllTogether->Update();
// 	canvasPi0SpectraAllTogether->Print(Form("%s/Pi0_Spectra_pPb2040_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 		canvasPi0SpectraAllTogether->cd();
// 	histo2DPi0SpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldPi0SysErrpPb4060->Draw("E2same");
// 	histoPCMYieldPi0pPb4060->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoPi0InputMCWOWeightspPb4060, 1., lineStyleMCpPb4060, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightspPb4060->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCWOWeightsAddedSigpPb4060, 1., lineStyleMCAddSigpPb4060, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightsAddedSigpPb4060->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCpPb4060, 1., lineStyleMCpPb4060, colorCombMCReweightedpPb);  
// 	histoPi0InputMCpPb4060->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCAddedSigpPb4060, 1., lineStyleMCAddSigpPb4060, colorCombMCReweightedAddSigpPb);  
// 	histoPi0InputMCAddedSigpPb4060->Draw("same,hist,c");    
// 
// 	labelSpectraPi0Label->Draw();
// 	TLatex *labelSpectra4060 = new TLatex(0.55,0.88,collisionSystempPb4060.Data());
// 	SetStyleTLatex( labelSpectra4060, 0.035  ,4);
// 	labelSpectra4060->Draw(); 
// 	
// 	TLegend* legendSpectraPi0MC4060 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraPi0MC4060->SetFillColor(0);
// 	legendSpectraPi0MC4060->SetLineColor(0);
// 	legendSpectraPi0MC4060->SetTextSize(0.025);
// 	legendSpectraPi0MC4060->SetNColumns(1);
// 	legendSpectraPi0MC4060->SetMargin(0.2);
// 	legendSpectraPi0MC4060->AddEntry(graphPCMYieldPi0SysErrpPb4060,"measured","pf");
// 	legendSpectraPi0MC4060->AddEntry(histoPi0InputMCWOWeightspPb4060,"HIJING","l");
// 	legendSpectraPi0MC4060->AddEntry(histoPi0InputMCWOWeightsAddedSigpPb4060,"added Signals","l");
// 	legendSpectraPi0MC4060->AddEntry(histoPi0InputMCpPb4060,"HIJING reweighted","l");
// 	legendSpectraPi0MC4060->AddEntry(histoPi0InputMCAddedSigpPb4060,"added Signals reweighted","l");
// 	
// 	legendSpectraPi0MC4060->Draw();
// 	
// 	canvasPi0SpectraAllTogether->Update();
// 	canvasPi0SpectraAllTogether->Print(Form("%s/Pi0_Spectra_pPb4060_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 		canvasPi0SpectraAllTogether->cd();
// 	histo2DPi0SpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldPi0SysErrpPb6080->Draw("E2same");
// 	histoPCMYieldPi0pPb6080->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoPi0InputMCWOWeightspPb6080, 1., lineStyleMCpPb6080, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightspPb6080->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCWOWeightsAddedSigpPb6080, 1., lineStyleMCAddSigpPb6080, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightsAddedSigpPb6080->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCpPb6080, 1., lineStyleMCpPb6080, colorCombMCReweightedpPb);  
// 	histoPi0InputMCpPb6080->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCAddedSigpPb6080, 1., lineStyleMCAddSigpPb6080, colorCombMCReweightedAddSigpPb);  
// 	histoPi0InputMCAddedSigpPb6080->Draw("same,hist,c");    
// 
// 	labelSpectraPi0Label->Draw();
// 	TLatex *labelSpectra6080 = new TLatex(0.55,0.88,collisionSystempPb6080.Data());
// 	SetStyleTLatex( labelSpectra6080, 0.035  ,4);
// 	labelSpectra6080->Draw(); 
// 	
// 	TLegend* legendSpectraPi0MC6080 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraPi0MC6080->SetFillColor(0);
// 	legendSpectraPi0MC6080->SetLineColor(0);
// 	legendSpectraPi0MC6080->SetTextSize(0.025);
// 	legendSpectraPi0MC6080->SetNColumns(1);
// 	legendSpectraPi0MC6080->SetMargin(0.2);
// 	legendSpectraPi0MC6080->AddEntry(graphPCMYieldPi0SysErrpPb6080,"measured","pf");
// 	legendSpectraPi0MC6080->AddEntry(histoPi0InputMCWOWeightspPb6080,"HIJING","l");
// 	legendSpectraPi0MC6080->AddEntry(histoPi0InputMCWOWeightsAddedSigpPb6080,"added Signals","l");
// 	legendSpectraPi0MC6080->AddEntry(histoPi0InputMCpPb6080,"HIJING reweighted","l");
// 	legendSpectraPi0MC6080->AddEntry(histoPi0InputMCAddedSigpPb6080,"added Signals reweighted","l");
// 	
// 	legendSpectraPi0MC6080->Draw();
// 	
// 	canvasPi0SpectraAllTogether->Update();
// 	canvasPi0SpectraAllTogether->Print(Form("%s/Pi0_Spectra_pPb6080_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 		canvasPi0SpectraAllTogether->cd();
// 	histo2DPi0SpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldPi0SysErrpPb60100->Draw("E2same");
// 	histoPCMYieldPi0pPb60100->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoPi0InputMCWOWeightspPb60100, 1., lineStyleMCpPb60100, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightspPb60100->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCWOWeightsAddedSigpPb60100, 1., lineStyleMCAddSigpPb60100, colorCombMCpPb);  
// 	histoPi0InputMCWOWeightsAddedSigpPb60100->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCpPb60100, 1., lineStyleMCpPb60100, colorCombMCReweightedpPb);  
// 	histoPi0InputMCpPb60100->Draw("same,hist,c");    
// 	SetStyleHisto(histoPi0InputMCAddedSigpPb60100, 1., lineStyleMCAddSigpPb60100, colorCombMCReweightedAddSigpPb);  
// 	histoPi0InputMCAddedSigpPb60100->Draw("same,hist,c");    
// 
// 	labelSpectraPi0Label->Draw();
// 	TLatex *labelSpectra60100 = new TLatex(0.55,0.88,collisionSystempPb60100.Data());
// 	SetStyleTLatex( labelSpectra60100, 0.035  ,4);
// 	labelSpectra60100->Draw(); 
// 	
// 	TLegend* legendSpectraPi0MC60100 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraPi0MC60100->SetFillColor(0);
// 	legendSpectraPi0MC60100->SetLineColor(0);
// 	legendSpectraPi0MC60100->SetTextSize(0.025);
// 	legendSpectraPi0MC60100->SetNColumns(1);
// 	legendSpectraPi0MC60100->SetMargin(0.2);
// 	legendSpectraPi0MC60100->AddEntry(graphPCMYieldPi0SysErrpPb60100,"measured","pf");
// 	legendSpectraPi0MC60100->AddEntry(histoPi0InputMCWOWeightspPb60100,"HIJING","l");
// 	legendSpectraPi0MC60100->AddEntry(histoPi0InputMCWOWeightsAddedSigpPb60100,"added Signals","l");
// 	legendSpectraPi0MC60100->AddEntry(histoPi0InputMCpPb60100,"HIJING reweighted","l");
// 	legendSpectraPi0MC60100->AddEntry(histoPi0InputMCAddedSigpPb60100,"added Signals reweighted","l");
// 	
// 	legendSpectraPi0MC60100->Draw();
// 	
// 	canvasPi0SpectraAllTogether->Update();
// 	canvasPi0SpectraAllTogether->Print(Form("%s/Pi0_Spectra_pPb60100_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 	
// 	//	**********************************************************************************************************************
// 	//	****************************************Eta Spectra compared to MC****************************************************
// 	//	**********************************************************************************************************************
// 	
// 	TCanvas* canvasEtaSpectraAllTogether = new TCanvas("canvasEtaSpectraAllTogether","",200,10,1200,1100);  // gives the page size
// 	DrawGammaCanvasSettings( canvasEtaSpectraAllTogether,  0.13, 0.01, 0.015, 0.08);
// 	
// 	canvasEtaSpectraAllTogether->SetLogy();
// 	canvasEtaSpectraAllTogether->SetLogx();
// 	TH2F * histo2DEtaSpectraAll;
// 	histo2DEtaSpectraAll = new TH2F("histo2DEtaSpectraAll","histo2DEtaSpectraAll",1000,0.4,15.,1000,1e-8,2e2 );
// 	SetStyleHistoTH2ForGraphs(histo2DEtaSpectraAll, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{d#it{N}^{AA}}{d#it{p}_{T} d#it{y}}}{ #frac{1}{N_{evt}^{pp}}#frac{d#it{N}^{pp}}{d#it{p}_{T} d#it{y}}}
// 	histo2DEtaSpectraAll->GetXaxis()->SetLabelOffset(-0.01);
// 	histo2DEtaSpectraAll->GetYaxis()->SetLabelOffset(0.01);
// 	histo2DEtaSpectraAll->DrawCopy(); 
// 
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldEtaSysErrpPb, markerStylepPb,markerSizepPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE);
// 	graphPCMYieldEtaSysErrpPb->Draw("E2same");
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldEtaSysErrpPb0020, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020, widthLinesBoxes, kTRUE);
// 	graphPCMYieldEtaSysErrpPb0020->Draw("E2same");
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldEtaSysErrpPb2040, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040, widthLinesBoxes, kTRUE);
// 	graphPCMYieldEtaSysErrpPb2040->Draw("E2same");
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldEtaSysErrpPb4060, markerStylepPb4060,markerSizepPb4060, colorCombpPb4060 , colorCombpPb4060, widthLinesBoxes, kTRUE);
// 	graphPCMYieldEtaSysErrpPb4060->Draw("E2same");
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldEtaSysErrpPb6080, markerStylepPb6080,markerSizepPb6080, colorCombpPb6080 , colorCombpPb6080, widthLinesBoxes, kTRUE);
// 	graphPCMYieldEtaSysErrpPb6080->Draw("E2same"); 
// 	DrawGammaSetMarkerTGraphAsym(graphPCMYieldEtaSysErrpPb60100, markerStylepPb60100,markerSizepPb60100, colorCombpPb60100 , colorCombpPb60100, widthLinesBoxes, kTRUE);
// 	graphPCMYieldEtaSysErrpPb60100->Draw("E2same");
// 	
// 	DrawGammaSetMarker(histoPCMYieldEtapPb, markerStylepPb,markerSizepPb, colorCombpPb , colorCombpPb);
// 	histoPCMYieldEtapPb->Draw("p,same,e1");
// 	DrawGammaSetMarker(histoPCMYieldEtapPb0020, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020);
// 	histoPCMYieldEtapPb0020->Draw("p,same,e1");
// 	DrawGammaSetMarker(histoPCMYieldEtapPb2040, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	histoPCMYieldEtapPb2040->Draw("p,same,e1");
// 	DrawGammaSetMarker(histoPCMYieldEtapPb4060, markerStylepPb4060,markerSizepPb4060, colorCombpPb4060 , colorCombpPb4060);
// 	histoPCMYieldEtapPb4060->Draw("p,same,e1");
// 	DrawGammaSetMarker(histoPCMYieldEtapPb6080, markerStylepPb6080,markerSizepPb6080, colorCombpPb6080 , colorCombpPb6080);
// 	histoPCMYieldEtapPb6080->Draw("same");
// 	DrawGammaSetMarker(histoPCMYieldEtapPb60100, markerStylepPb60100,markerSizepPb60100, colorCombpPb60100 , colorCombpPb60100);
// 	histoPCMYieldEtapPb60100->Draw("same");
// 	
// 
// 	TLatex *labelSpectraEtaLabel = new TLatex(0.55,0.92,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
// 	SetStyleTLatex( labelSpectraEtaLabel, 0.03  ,4);
// 	labelSpectraEtaLabel->Draw();
// 	
// 	TLegend* legendSpectra = new TLegend(0.16,0.09,0.63,0.3);
// 	legendSpectra->SetFillColor(0);
// 	legendSpectra->SetLineColor(0);
// 	legendSpectra->SetTextSize(0.025);
// 	legendSpectra->SetNColumns(2);
// 	legendSpectra->SetMargin(0.2);
// 	legendSpectra->AddEntry(graphPCMYieldEtaSysErrpPb,"PCM","pf");
// 	legendSpectra->AddEntry((TObject*)0, collisionSystempPb.Data(),"");
// 	legendSpectra->AddEntry(graphPCMYieldEtaSysErrpPb0020,"PCM","pf");
// 	legendSpectra->AddEntry((TObject*)0, collisionSystempPb0020.Data(),"");
// 	legendSpectra->AddEntry(graphPCMYieldEtaSysErrpPb2040,"PCM","pf");
// 	legendSpectra->AddEntry((TObject*)0, collisionSystempPb2040.Data(),"");
// 	legendSpectra->AddEntry(graphPCMYieldEtaSysErrpPb4060,"PCM","pf");
// 	legendSpectra->AddEntry((TObject*)0, collisionSystempPb4060.Data(),"");
// 	legendSpectra->AddEntry(graphPCMYieldEtaSysErrpPb6080,"PCM","pf");
// 	legendSpectra->AddEntry((TObject*)0, collisionSystempPb6080.Data(),"");
// 	legendSpectra->AddEntry(graphPCMYieldEtaSysErrpPb60100,"PCM","pf");
// 	legendSpectra->AddEntry((TObject*)0, collisionSystempPb60100.Data(),"");
// 	legendSpectra->Draw();
// 	
// 	canvasEtaSpectraAllTogether->Update();
// 	canvasEtaSpectraAllTogether->Print(Form("%s/Eta_Spectra_All.%s",outputDir.Data(),suffix.Data()));
// 
// 	
// 	canvasEtaSpectraAllTogether->cd();
// 	histo2DEtaSpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldEtaSysErrpPb->Draw("E2same");
// 	histoPCMYieldEtapPb->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoEtaInputMCWOWeightspPb, 1., lineStyleMCpPb, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightspPb->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCWOWeightsAddedSigpPb, 1., lineStyleMCAddSigpPb, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightsAddedSigpPb->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCpPb, 1., lineStyleMCpPb, colorCombMCReweightedpPb);  
// 	histoEtaInputMCpPb->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCAddedSigpPb, 1., lineStyleMCAddSigpPb, colorCombMCReweightedAddSigpPb);  
// 	histoEtaInputMCAddedSigpPb->Draw("same,hist,c");    
// 
// 	labelSpectraEtaLabel->Draw();
// 	labelSpectraMinBias->Draw(); 
// 	
// 	TLegend* legendSpectraEtaMC = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraEtaMC->SetFillColor(0);
// 	legendSpectraEtaMC->SetLineColor(0);
// 	legendSpectraEtaMC->SetTextSize(0.025);
// 	legendSpectraEtaMC->SetNColumns(1);
// 	legendSpectraEtaMC->SetMargin(0.2);
// 	legendSpectraEtaMC->AddEntry(graphPCMYieldEtaSysErrpPb,"measured","pf");
// 	legendSpectraEtaMC->AddEntry(histoEtaInputMCWOWeightspPb,"HIJING","l");
// 	legendSpectraEtaMC->AddEntry(histoEtaInputMCWOWeightsAddedSigpPb,"added Signals","l");
// 	legendSpectraEtaMC->AddEntry(histoEtaInputMCpPb,"HIJING reweighted","l");
// 	legendSpectraEtaMC->AddEntry(histoEtaInputMCAddedSigpPb,"added Signals reweighted","l");
// 	
// 	legendSpectraEtaMC->Draw();
// 	
// 	canvasEtaSpectraAllTogether->Update();
// 	canvasEtaSpectraAllTogether->Print(Form("%s/Eta_Spectra_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 	canvasEtaSpectraAllTogether->cd();
// 	histo2DEtaSpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldEtaSysErrpPb->Draw("E2same");
// 	histoPCMYieldEtapPb->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoEtaInputMCWOWeightspPb, 1., lineStyleMCpPb, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightspPb->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCWOWeightsAddedSigpPb, 1., lineStyleMCAddSigpPb, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightsAddedSigpPb->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCWOWeightspPb_DPMJET, 1., lineStyleMCpPb, kRed-7);  
// 	histoEtaInputMCWOWeightspPb_DPMJET->Draw("same,hist,c");    
// 	
// 	labelSpectraEtaLabel->Draw();
// 	labelSpectraMinBias->Draw(); 
// 	
// 	TLegend* legendSpectraEtaMCAll = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraEtaMCAll->SetFillColor(0);
// 	legendSpectraEtaMCAll->SetLineColor(0);
// 	legendSpectraEtaMCAll->SetTextSize(0.025);
// 	legendSpectraEtaMCAll->SetNColumns(1);
// 	legendSpectraEtaMCAll->SetMargin(0.2);
// 	legendSpectraEtaMCAll->AddEntry(graphPCMYieldEtaSysErrpPb,"measured","pf");
// 	legendSpectraEtaMCAll->AddEntry(histoEtaInputMCWOWeightspPb,"HIJING","l");
// 	legendSpectraEtaMCAll->AddEntry(histoEtaInputMCWOWeightsAddedSigpPb,"added Signals","l");
// 	legendSpectraEtaMCAll->AddEntry(histoEtaInputMCWOWeightspPb_DPMJET,"DPMJET","l");
// 	legendSpectraEtaMCAll->Draw();
// 	
// 	canvasEtaSpectraAllTogether->Update();
// 	canvasEtaSpectraAllTogether->Print(Form("%s/Eta_Spectra_WithMC_HIJINGAndDPMJET.%s",outputDir.Data(),suffix.Data()));
// 
// 	
// 	canvasEtaSpectraAllTogether->cd();
// 	histo2DEtaSpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldEtaSysErrpPb0020->Draw("E2same");
// 	histoPCMYieldEtapPb0020->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoEtaInputMCWOWeightspPb0020, 1., lineStyleMCpPb0020, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightspPb0020->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCWOWeightsAddedSigpPb0020, 1., lineStyleMCAddSigpPb0020, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightsAddedSigpPb0020->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCpPb0020, 1., lineStyleMCpPb0020, colorCombMCReweightedpPb);  
// 	histoEtaInputMCpPb0020->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCAddedSigpPb0020, 1., lineStyleMCAddSigpPb0020, colorCombMCReweightedAddSigpPb);  
// 	histoEtaInputMCAddedSigpPb0020->Draw("same,hist,c");    
// 
// 	labelSpectraEtaLabel->Draw();
// 	labelSpectra0020->Draw(); 
// 	
// 	TLegend* legendSpectraEtaMC0020 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraEtaMC0020->SetFillColor(0);
// 	legendSpectraEtaMC0020->SetLineColor(0);
// 	legendSpectraEtaMC0020->SetTextSize(0.025);
// 	legendSpectraEtaMC0020->SetNColumns(1);
// 	legendSpectraEtaMC0020->SetMargin(0.2);
// 	legendSpectraEtaMC0020->AddEntry(graphPCMYieldEtaSysErrpPb0020,"measured","pf");
// 	legendSpectraEtaMC0020->AddEntry(histoEtaInputMCWOWeightspPb0020,"HIJING","l");
// 	legendSpectraEtaMC0020->AddEntry(histoEtaInputMCWOWeightsAddedSigpPb0020,"added Signals","l");
// 	legendSpectraEtaMC0020->AddEntry(histoEtaInputMCpPb0020,"HIJING reweighted","l");
// 	legendSpectraEtaMC0020->AddEntry(histoEtaInputMCAddedSigpPb0020,"added Signals reweighted","l");
// 	
// 	legendSpectraEtaMC0020->Draw();
// 	
// 	canvasEtaSpectraAllTogether->Update();
// 	canvasEtaSpectraAllTogether->Print(Form("%s/Eta_Spectra_pPb0020_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 	canvasEtaSpectraAllTogether->cd();
// 	histo2DEtaSpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldEtaSysErrpPb2040->Draw("E2same");
// 	histoPCMYieldEtapPb2040->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoEtaInputMCWOWeightspPb2040, 1., lineStyleMCpPb2040, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightspPb2040->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCWOWeightsAddedSigpPb2040, 1., lineStyleMCAddSigpPb2040, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightsAddedSigpPb2040->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCpPb2040, 1., lineStyleMCpPb2040, colorCombMCReweightedpPb);  
// 	histoEtaInputMCpPb2040->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCAddedSigpPb2040, 1., lineStyleMCAddSigpPb2040, colorCombMCReweightedAddSigpPb);  
// 	histoEtaInputMCAddedSigpPb2040->Draw("same,hist,c");    
// 
// 	labelSpectraEtaLabel->Draw();
// 	labelSpectra2040->Draw(); 
// 	
// 	TLegend* legendSpectraEtaMC2040 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraEtaMC2040->SetFillColor(0);
// 	legendSpectraEtaMC2040->SetLineColor(0);
// 	legendSpectraEtaMC2040->SetTextSize(0.025);
// 	legendSpectraEtaMC2040->SetNColumns(1);
// 	legendSpectraEtaMC2040->SetMargin(0.2);
// 	legendSpectraEtaMC2040->AddEntry(graphPCMYieldEtaSysErrpPb2040,"measured","pf");
// 	legendSpectraEtaMC2040->AddEntry(histoEtaInputMCWOWeightspPb2040,"HIJING","l");
// 	legendSpectraEtaMC2040->AddEntry(histoEtaInputMCWOWeightsAddedSigpPb2040,"added Signals","l");
// 	legendSpectraEtaMC2040->AddEntry(histoEtaInputMCpPb2040,"HIJING reweighted","l");
// 	legendSpectraEtaMC2040->AddEntry(histoEtaInputMCAddedSigpPb2040,"added Signals reweighted","l");
// 	
// 	legendSpectraEtaMC2040->Draw();
// 	
// 	canvasEtaSpectraAllTogether->Update();
// 	canvasEtaSpectraAllTogether->Print(Form("%s/Eta_Spectra_pPb2040_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 		canvasEtaSpectraAllTogether->cd();
// 	histo2DEtaSpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldEtaSysErrpPb4060->Draw("E2same");
// 	histoPCMYieldEtapPb4060->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoEtaInputMCWOWeightspPb4060, 1., lineStyleMCpPb4060, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightspPb4060->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCWOWeightsAddedSigpPb4060, 1., lineStyleMCAddSigpPb4060, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightsAddedSigpPb4060->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCpPb4060, 1., lineStyleMCpPb4060, colorCombMCReweightedpPb);  
// 	histoEtaInputMCpPb4060->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCAddedSigpPb4060, 1., lineStyleMCAddSigpPb4060, colorCombMCReweightedAddSigpPb);  
// 	histoEtaInputMCAddedSigpPb4060->Draw("same,hist,c");    
// 
// 	labelSpectraEtaLabel->Draw();
// 	labelSpectra4060->Draw(); 
// 	
// 	TLegend* legendSpectraEtaMC4060 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraEtaMC4060->SetFillColor(0);
// 	legendSpectraEtaMC4060->SetLineColor(0);
// 	legendSpectraEtaMC4060->SetTextSize(0.025);
// 	legendSpectraEtaMC4060->SetNColumns(1);
// 	legendSpectraEtaMC4060->SetMargin(0.2);
// 	legendSpectraEtaMC4060->AddEntry(graphPCMYieldEtaSysErrpPb4060,"measured","pf");
// 	legendSpectraEtaMC4060->AddEntry(histoEtaInputMCWOWeightspPb4060,"HIJING","l");
// 	legendSpectraEtaMC4060->AddEntry(histoEtaInputMCWOWeightsAddedSigpPb4060,"added Signals","l");
// 	legendSpectraEtaMC4060->AddEntry(histoEtaInputMCpPb4060,"HIJING reweighted","l");
// 	legendSpectraEtaMC4060->AddEntry(histoEtaInputMCAddedSigpPb4060,"added Signals reweighted","l");
// 	
// 	legendSpectraEtaMC4060->Draw();
// 	
// 	canvasEtaSpectraAllTogether->Update();
// 	canvasEtaSpectraAllTogether->Print(Form("%s/Eta_Spectra_pPb4060_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 		canvasEtaSpectraAllTogether->cd();
// 	histo2DEtaSpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldEtaSysErrpPb6080->Draw("E2same");
// 	histoPCMYieldEtapPb6080->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoEtaInputMCWOWeightspPb6080, 1., lineStyleMCpPb6080, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightspPb6080->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCWOWeightsAddedSigpPb6080, 1., lineStyleMCAddSigpPb6080, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightsAddedSigpPb6080->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCpPb6080, 1., lineStyleMCpPb6080, colorCombMCReweightedpPb);  
// 	histoEtaInputMCpPb6080->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCAddedSigpPb6080, 1., lineStyleMCAddSigpPb6080, colorCombMCReweightedAddSigpPb);  
// 	histoEtaInputMCAddedSigpPb6080->Draw("same,hist,c");    
// 
// 	labelSpectraEtaLabel->Draw();
// 	labelSpectra6080->Draw(); 
// 	
// 	TLegend* legendSpectraEtaMC6080 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraEtaMC6080->SetFillColor(0);
// 	legendSpectraEtaMC6080->SetLineColor(0);
// 	legendSpectraEtaMC6080->SetTextSize(0.025);
// 	legendSpectraEtaMC6080->SetNColumns(1);
// 	legendSpectraEtaMC6080->SetMargin(0.2);
// 	legendSpectraEtaMC6080->AddEntry(graphPCMYieldEtaSysErrpPb6080,"measured","pf");
// 	legendSpectraEtaMC6080->AddEntry(histoEtaInputMCWOWeightspPb6080,"HIJING","l");
// 	legendSpectraEtaMC6080->AddEntry(histoEtaInputMCWOWeightsAddedSigpPb6080,"added Signals","l");
// 	legendSpectraEtaMC6080->AddEntry(histoEtaInputMCpPb6080,"HIJING reweighted","l");
// 	legendSpectraEtaMC6080->AddEntry(histoEtaInputMCAddedSigpPb6080,"added Signals reweighted","l");
// 	
// 	legendSpectraEtaMC6080->Draw();
// 	
// 	canvasEtaSpectraAllTogether->Update();
// 	canvasEtaSpectraAllTogether->Print(Form("%s/Eta_Spectra_pPb6080_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 		canvasEtaSpectraAllTogether->cd();
// 	histo2DEtaSpectraAll->DrawCopy(); 
// 	
// 	graphPCMYieldEtaSysErrpPb60100->Draw("E2same");
// 	histoPCMYieldEtapPb60100->Draw("p,same,e1");
// 
// 	SetStyleHisto(histoEtaInputMCWOWeightspPb60100, 1., lineStyleMCpPb60100, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightspPb60100->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCWOWeightsAddedSigpPb60100, 1., lineStyleMCAddSigpPb60100, colorCombMCpPb);  
// 	histoEtaInputMCWOWeightsAddedSigpPb60100->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCpPb60100, 1., lineStyleMCpPb60100, colorCombMCReweightedpPb);  
// 	histoEtaInputMCpPb60100->Draw("same,hist,c");    
// 	SetStyleHisto(histoEtaInputMCAddedSigpPb60100, 1., lineStyleMCAddSigpPb60100, colorCombMCReweightedAddSigpPb);  
// 	histoEtaInputMCAddedSigpPb60100->Draw("same,hist,c");    
// 
// 	labelSpectraEtaLabel->Draw();
// 	labelSpectra60100->Draw(); 
// 	
// 	TLegend* legendSpectraEtaMC60100 = new TLegend(0.16,0.09,0.73,0.3);
// 	legendSpectraEtaMC60100->SetFillColor(0);
// 	legendSpectraEtaMC60100->SetLineColor(0);
// 	legendSpectraEtaMC60100->SetTextSize(0.025);
// 	legendSpectraEtaMC60100->SetNColumns(1);
// 	legendSpectraEtaMC60100->SetMargin(0.2);
// 	legendSpectraEtaMC60100->AddEntry(graphPCMYieldEtaSysErrpPb60100,"measured","pf");
// 	legendSpectraEtaMC60100->AddEntry(histoEtaInputMCWOWeightspPb60100,"HIJING","l");
// 	legendSpectraEtaMC60100->AddEntry(histoEtaInputMCWOWeightsAddedSigpPb60100,"added Signals","l");
// 	legendSpectraEtaMC60100->AddEntry(histoEtaInputMCpPb60100,"HIJING reweighted","l");
// 	legendSpectraEtaMC60100->AddEntry(histoEtaInputMCAddedSigpPb60100,"added Signals reweighted","l");
// 	
// 	legendSpectraEtaMC60100->Draw();
// 	
// 	canvasEtaSpectraAllTogether->Update();
// 	canvasEtaSpectraAllTogether->Print(Form("%s/Eta_Spectra_pPb60100_WithMC_HIJING.%s",outputDir.Data(),suffix.Data()));
// 
// 	
// 	
// 	      
// 	TCanvas* canvasFraction2 = new TCanvas("canvasFraction2","",1550,1200);  // gives the page size
// 	canvasFraction2->SetTickx();
// 	canvasFraction2->SetTicky();
// 	canvasFraction2->SetGridx(0);
// 	canvasFraction2->SetGridy(0);
// 	canvasFraction2->SetLogy(0);
// 	canvasFraction2->SetLeftMargin(0.13);
// 	canvasFraction2->SetRightMargin(0.02);
// 	canvasFraction2->SetTopMargin(0.02);
// 	canvasFraction2->SetFillColor(0);
//  
// 	
// 	TF1* fitYieldDataPi0pPb = NULL;
// 	fitYieldDataPi0pPb = FitObject("l","fitYieldDataPi0pPb","Pi0",histoPCMYieldPi0pPb,0.3,14,NULL,"QNRME+");
// 	TH1D* histoPi0RatioDatatoFitpPb = CalculateHistoRatioToFit (histoPCMYieldPi0pPb, fitYieldDataPi0pPb);
// 	TH1D* histoPi0RatioMCtoDataFitpPb = CalculateHistoRatioToFit (histoPi0InputMCpPb, fitYieldDataPi0pPb);
// 	TH1D* histoPi0RatioMCUnweightedtoDataFitpPb = NULL;
// 	if (histoPi0InputMCWOWeightspPb) histoPi0RatioMCUnweightedtoDataFitpPb = CalculateHistoRatioToFit (histoPi0InputMCWOWeightspPb, fitYieldDataPi0pPb);
// 	canvasFraction2->cd();
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb) DrawGammaSetMarker(histoPi0RatioMCUnweightedtoDataFitpPb, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoPi0RatioMCtoDataFitpPb, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoPi0RatioDatatoFitpPb, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoPi0RatioDatatoFitpPb,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 15);
// 	histoPi0RatioDatatoFitpPb->Draw("same,e,p");  
// 	if (runDrawReweighted) histoPi0RatioMCtoDataFitpPb->Draw("same,hist,c");  
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb) histoPi0RatioMCUnweightedtoDataFitpPb->Draw("same,hist,c");  
// 	
// 	TLegend* legendPi0RatiopPb = new TLegend(0.16,0.81,0.4,0.9);
// 	legendPi0RatiopPb->SetFillColor(0);
// 	legendPi0RatiopPb->SetLineColor(0);
// 	legendPi0RatiopPb->SetTextSize(0.025);
// //    legendPi0RatiopPb->SetNColumns(3);
// 	legendPi0RatiopPb->SetMargin(0.2);
// 	legendPi0RatiopPb->AddEntry(histoPi0RatioDatatoFitpPb,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendPi0RatiopPb->AddEntry(histoPi0RatioMCtoDataFitpPb,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb) legendPi0RatiopPb->AddEntry(histoPi0RatioMCUnweightedtoDataFitpPb,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendPi0RatiopPb->Draw();
// 	TLatex *labelRatioMCDatapPb = new TLatex(0.2,0.92,collisionSystempPb.Data());
// 	SetStyleTLatex( labelRatioMCDatapPb, 0.04,4);
// 	labelRatioMCDatapPb->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Pi0_RatioMCToDataFit_pPb.%s",outputDir.Data(),suffix.Data()));
// 
// 	TF1* fitYieldDataPi0pPb0020 = NULL;
// 	fitYieldDataPi0pPb0020 = FitObject("l","fitYieldDataPi0pPb0020","Pi0",histoPCMYieldPi0pPb0020,0.3,14,NULL,"QNRME+");
// 	TH1D* histoPi0RatioDatatoFitpPb0020 = CalculateHistoRatioToFit (histoPCMYieldPi0pPb0020, fitYieldDataPi0pPb0020);
// 	TH1D* histoPi0RatioMCtoDataFitpPb0020 = CalculateHistoRatioToFit (histoPi0InputMCpPb0020, fitYieldDataPi0pPb0020);
// 	TH1D* histoPi0RatioMCUnweightedtoDataFitpPb0020 = NULL;
// 	if (histoPi0InputMCWOWeightspPb0020) histoPi0RatioMCUnweightedtoDataFitpPb0020 = CalculateHistoRatioToFit (histoPi0InputMCWOWeightspPb0020, fitYieldDataPi0pPb0020);
// 	canvasFraction2->cd();
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb0020) DrawGammaSetMarker(histoPi0RatioMCUnweightedtoDataFitpPb0020, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoPi0RatioMCtoDataFitpPb0020, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoPi0RatioDatatoFitpPb0020, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoPi0RatioDatatoFitpPb0020,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 15);
// 	histoPi0RatioDatatoFitpPb0020->Draw("same,e,p");  
// 	if (runDrawReweighted) histoPi0RatioMCtoDataFitpPb0020->Draw("same,hist,c");  
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb0020) histoPi0RatioMCUnweightedtoDataFitpPb0020->Draw("same,hist,c");  
// 	
// 	TLegend* legendPi0RatiopPb0020 = new TLegend(0.16,0.81,0.4,0.9);
// 	legendPi0RatiopPb0020->SetFillColor(0);
// 	legendPi0RatiopPb0020->SetLineColor(0);
// 	legendPi0RatiopPb0020->SetTextSize(0.025);
// //    legendPi0RatiopPb0020->SetNColumns(3);
// 	legendPi0RatiopPb0020->SetMargin(0.2);
// 	legendPi0RatiopPb0020->AddEntry(histoPi0RatioDatatoFitpPb0020,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendPi0RatiopPb0020->AddEntry(histoPi0RatioMCtoDataFitpPb0020,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb0020) legendPi0RatiopPb0020->AddEntry(histoPi0RatioMCUnweightedtoDataFitpPb0020,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendPi0RatiopPb0020->Draw();
// 	TLatex *labelRatioMCDatapPb0020 = new TLatex(0.2,0.92,collisionSystempPb0020.Data());
// 	SetStyleTLatex( labelRatioMCDatapPb0020, 0.04,4);
// 	labelRatioMCDatapPb0020->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Pi0_RatioMCToDataFit_pPb0020.%s",outputDir.Data(),suffix.Data()));
// 
// 	TF1* fitYieldDataPi0pPb2040 = NULL;
// 	fitYieldDataPi0pPb2040 = FitObject("l","fitYieldDataPi0pPb2040","Pi0",histoPCMYieldPi0pPb2040,0.3,14,NULL,"QNRME+");
// 	TH1D* histoPi0RatioDatatoFitpPb2040 = CalculateHistoRatioToFit (histoPCMYieldPi0pPb2040, fitYieldDataPi0pPb2040);
// 	TH1D* histoPi0RatioMCtoDataFitpPb2040 = CalculateHistoRatioToFit (histoPi0InputMCpPb2040, fitYieldDataPi0pPb2040);
// 	TH1D* histoPi0RatioMCUnweightedtoDataFitpPb2040 = NULL;
// 	if (histoPi0InputMCWOWeightspPb2040) histoPi0RatioMCUnweightedtoDataFitpPb2040 = CalculateHistoRatioToFit (histoPi0InputMCWOWeightspPb2040, fitYieldDataPi0pPb2040);
// 	canvasFraction2->cd();
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb2040) DrawGammaSetMarker(histoPi0RatioMCUnweightedtoDataFitpPb2040, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoPi0RatioMCtoDataFitpPb2040, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoPi0RatioDatatoFitpPb2040, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoPi0RatioDatatoFitpPb2040,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 15);
// 	histoPi0RatioDatatoFitpPb2040->Draw("same,e,p");  
// 	if (runDrawReweighted) histoPi0RatioMCtoDataFitpPb2040->Draw("same,hist,c");  
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb2040) histoPi0RatioMCUnweightedtoDataFitpPb2040->Draw("same,hist,c");  
// 	
// 	TLegend* legendPi0RatiopPb2040 = new TLegend(0.16,0.81,0.4,0.9);
// 	legendPi0RatiopPb2040->SetFillColor(0);
// 	legendPi0RatiopPb2040->SetLineColor(0);
// 	legendPi0RatiopPb2040->SetTextSize(0.025);
// //    legendPi0RatiopPb2040->SetNColumns(3);
// 	legendPi0RatiopPb2040->SetMargin(0.2);
// 	legendPi0RatiopPb2040->AddEntry(histoPi0RatioDatatoFitpPb2040,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendPi0RatiopPb2040->AddEntry(histoPi0RatioMCtoDataFitpPb2040,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb2040) legendPi0RatiopPb2040->AddEntry(histoPi0RatioMCUnweightedtoDataFitpPb2040,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendPi0RatiopPb2040->Draw();
// 	TLatex *labelRatioMCDatapPb2040 = new TLatex(0.2,0.92,collisionSystempPb2040.Data());
// 	SetStyleTLatex( labelRatioMCDatapPb2040, 0.04,4);
// 	labelRatioMCDatapPb2040->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Pi0_RatioMCToDataFit_pPb2040.%s",outputDir.Data(),suffix.Data()));
// 
// 	TF1* fitYieldDataPi0pPb4060 = NULL;
// 	fitYieldDataPi0pPb4060 = FitObject("l","fitYieldDataPi0pPb4060","Pi0",histoPCMYieldPi0pPb4060,0.3,14,NULL,"QNRME+");
// 	TH1D* histoPi0RatioDatatoFitpPb4060 = CalculateHistoRatioToFit (histoPCMYieldPi0pPb4060, fitYieldDataPi0pPb4060);
// 	TH1D* histoPi0RatioMCtoDataFitpPb4060 = CalculateHistoRatioToFit (histoPi0InputMCpPb4060, fitYieldDataPi0pPb4060);
// 	TH1D* histoPi0RatioMCUnweightedtoDataFitpPb4060 = NULL;
// 	if (histoPi0InputMCWOWeightspPb4060) histoPi0RatioMCUnweightedtoDataFitpPb4060 = CalculateHistoRatioToFit (histoPi0InputMCWOWeightspPb4060, fitYieldDataPi0pPb4060);
// 	canvasFraction2->cd();
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb4060) DrawGammaSetMarker(histoPi0RatioMCUnweightedtoDataFitpPb4060, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoPi0RatioMCtoDataFitpPb4060, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoPi0RatioDatatoFitpPb4060, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoPi0RatioDatatoFitpPb4060,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 15);
// 	histoPi0RatioDatatoFitpPb4060->Draw("same,e,p");  
// 	if (runDrawReweighted) histoPi0RatioMCtoDataFitpPb4060->Draw("same,hist,c");  
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb4060) histoPi0RatioMCUnweightedtoDataFitpPb4060->Draw("same,hist,c");  
// 	
// 	TLegend* legendPi0RatiopPb4060 = new TLegend(0.16,0.81,0.4,0.9);
// 	legendPi0RatiopPb4060->SetFillColor(0);
// 	legendPi0RatiopPb4060->SetLineColor(0);
// 	legendPi0RatiopPb4060->SetTextSize(0.025);
// //    legendPi0RatiopPb4060->SetNColumns(3);
// 	legendPi0RatiopPb4060->SetMargin(0.2);
// 	legendPi0RatiopPb4060->AddEntry(histoPi0RatioDatatoFitpPb4060,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendPi0RatiopPb4060->AddEntry(histoPi0RatioMCtoDataFitpPb4060,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb4060) legendPi0RatiopPb4060->AddEntry(histoPi0RatioMCUnweightedtoDataFitpPb4060,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendPi0RatiopPb4060->Draw();
// 	TLatex *labelRatioMCDatapPb4060 = new TLatex(0.2,0.92,collisionSystempPb4060.Data());
// 	SetStyleTLatex( labelRatioMCDatapPb4060, 0.04,4);
// 	labelRatioMCDatapPb4060->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Pi0_RatioMCToDataFit_pPb4060.%s",outputDir.Data(),suffix.Data()));
// 
// 	TF1* fitYieldDataPi0pPb6080 = NULL;
// 	fitYieldDataPi0pPb6080 = FitObject("l","fitYieldDataPi0pPb6080","Pi0",histoPCMYieldPi0pPb6080,0.3,14,NULL,"QNRME+");
// 	TH1D* histoPi0RatioDatatoFitpPb6080 = CalculateHistoRatioToFit (histoPCMYieldPi0pPb6080, fitYieldDataPi0pPb6080);
// 	TH1D* histoPi0RatioMCtoDataFitpPb6080 = CalculateHistoRatioToFit (histoPi0InputMCpPb6080, fitYieldDataPi0pPb6080);
// 	TH1D* histoPi0RatioMCUnweightedtoDataFitpPb6080 = NULL;
// 	if (histoPi0InputMCWOWeightspPb6080) histoPi0RatioMCUnweightedtoDataFitpPb6080 = CalculateHistoRatioToFit (histoPi0InputMCWOWeightspPb6080, fitYieldDataPi0pPb6080);
// 	canvasFraction2->cd();
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb6080) DrawGammaSetMarker(histoPi0RatioMCUnweightedtoDataFitpPb6080, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoPi0RatioMCtoDataFitpPb6080, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoPi0RatioDatatoFitpPb6080, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoPi0RatioDatatoFitpPb6080,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 15);
// 	histoPi0RatioDatatoFitpPb6080->Draw("same,e,p");  
// 	if (runDrawReweighted) histoPi0RatioMCtoDataFitpPb6080->Draw("same,hist,c");  
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb6080) histoPi0RatioMCUnweightedtoDataFitpPb6080->Draw("same,hist,c");  
// 	
// 	TLegend* legendPi0RatiopPb6080 = new TLegend(0.16,0.81,0.4,0.9);
// 	legendPi0RatiopPb6080->SetFillColor(0);
// 	legendPi0RatiopPb6080->SetLineColor(0);
// 	legendPi0RatiopPb6080->SetTextSize(0.025);
// //    legendPi0RatiopPb6080->SetNColumns(3);
// 	legendPi0RatiopPb6080->SetMargin(0.2);
// 	legendPi0RatiopPb6080->AddEntry(histoPi0RatioDatatoFitpPb6080,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendPi0RatiopPb6080->AddEntry(histoPi0RatioMCtoDataFitpPb6080,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb6080) legendPi0RatiopPb6080->AddEntry(histoPi0RatioMCUnweightedtoDataFitpPb6080,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendPi0RatiopPb6080->Draw();
// 	TLatex *labelRatioMCDatapPb6080 = new TLatex(0.2,0.92,collisionSystempPb6080.Data());
// 	SetStyleTLatex( labelRatioMCDatapPb6080, 0.04,4);
// 	labelRatioMCDatapPb6080->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Pi0_RatioMCToDataFit_pPb6080.%s",outputDir.Data(),suffix.Data()));
// 
// 	TF1* fitYieldDataPi0pPb60100 = NULL;
// 	fitYieldDataPi0pPb60100 = FitObject("l","fitYieldDataPi0pPb60100","Pi0",histoPCMYieldPi0pPb60100,0.3,14,NULL,"QNRME+");
// 	TH1D* histoPi0RatioDatatoFitpPb60100 = CalculateHistoRatioToFit (histoPCMYieldPi0pPb60100, fitYieldDataPi0pPb60100);
// 	TH1D* histoPi0RatioMCtoDataFitpPb60100 = CalculateHistoRatioToFit (histoPi0InputMCpPb60100, fitYieldDataPi0pPb60100);
// 	TH1D* histoPi0RatioMCUnweightedtoDataFitpPb60100 = NULL;
// 	if (histoPi0InputMCWOWeightspPb60100) histoPi0RatioMCUnweightedtoDataFitpPb60100 = CalculateHistoRatioToFit (histoPi0InputMCWOWeightspPb60100, fitYieldDataPi0pPb60100);
// 	canvasFraction2->cd();
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb60100) DrawGammaSetMarker(histoPi0RatioMCUnweightedtoDataFitpPb60100, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoPi0RatioMCtoDataFitpPb60100, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoPi0RatioDatatoFitpPb60100, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoPi0RatioDatatoFitpPb60100,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 15);
// 	histoPi0RatioDatatoFitpPb60100->Draw("same,e,p");  
// 	if (runDrawReweighted) histoPi0RatioMCtoDataFitpPb60100->Draw("same,hist,c");  
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb60100) histoPi0RatioMCUnweightedtoDataFitpPb60100->Draw("same,hist,c");  
// 	
// 	TLegend* legendPi0RatiopPb60100 = new TLegend(0.16,0.81,0.4,0.9);
// 	legendPi0RatiopPb60100->SetFillColor(0);
// 	legendPi0RatiopPb60100->SetLineColor(0);
// 	legendPi0RatiopPb60100->SetTextSize(0.025);
// //    legendPi0RatiopPb60100->SetNColumns(3);
// 	legendPi0RatiopPb60100->SetMargin(0.2);
// 	legendPi0RatiopPb60100->AddEntry(histoPi0RatioDatatoFitpPb60100,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendPi0RatiopPb60100->AddEntry(histoPi0RatioMCtoDataFitpPb60100,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoPi0RatioMCUnweightedtoDataFitpPb60100) legendPi0RatiopPb60100->AddEntry(histoPi0RatioMCUnweightedtoDataFitpPb60100,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendPi0RatiopPb60100->Draw();
// 	TLatex *labelRatioMCDatapPb60100 = new TLatex(0.2,0.92,collisionSystempPb60100.Data());
// 	SetStyleTLatex( labelRatioMCDatapPb60100, 0.04,4);
// 	labelRatioMCDatapPb60100->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Pi0_RatioMCToDataFit_pPb60100.%s",outputDir.Data(),suffix.Data()));
// 	
// 	//***********************************************************************************************************************************
// 	//******************************************Eta weighting comparison*****************************************************************
// 	//***********************************************************************************************************************************
// 	TF1* fitYieldDataEtapPb = NULL;
// 	fitYieldDataEtapPb = FitObject("l","fitYieldDataEtapPb","Eta",histoPCMYieldEtapPb,0.7,10,NULL,"QNRME+");
// 	TH1D* histoEtaRatioDatatoFitpPb = CalculateHistoRatioToFit (histoPCMYieldEtapPb, fitYieldDataEtapPb);
// 	TH1D* histoEtaRatioMCtoDataFitpPb = CalculateHistoRatioToFit (histoEtaInputMCpPb, fitYieldDataEtapPb);
// 	TH1D* histoEtaRatioMCUnweightedtoDataFitpPb = NULL;
// 	if (histoEtaInputMCWOWeightspPb) histoEtaRatioMCUnweightedtoDataFitpPb = CalculateHistoRatioToFit (histoEtaInputMCWOWeightspPb, fitYieldDataEtapPb);
// 	canvasFraction2->cd();
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitpPb, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoEtaRatioMCtoDataFitpPb, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoEtaRatioDatatoFitpPb, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitpPb,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 10);
// 	histoEtaRatioDatatoFitpPb->Draw("same,e,p");  
// 	if (runDrawReweighted) histoEtaRatioMCtoDataFitpPb->Draw("same,hist,c");  
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb) histoEtaRatioMCUnweightedtoDataFitpPb->Draw("same,hist,c");  
// 	
// 	TLegend* legendEtaRatiopPb = new TLegend(0.16,0.81,0.4,0.9);
// 	legendEtaRatiopPb->SetFillColor(0);
// 	legendEtaRatiopPb->SetLineColor(0);
// 	legendEtaRatiopPb->SetTextSize(0.025);
// //    legendEtaRatiopPb->SetNColumns(3);
// 	legendEtaRatiopPb->SetMargin(0.2);
// 	legendEtaRatiopPb->AddEntry(histoEtaRatioDatatoFitpPb,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendEtaRatiopPb->AddEntry(histoEtaRatioMCtoDataFitpPb,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb) legendEtaRatiopPb->AddEntry(histoEtaRatioMCUnweightedtoDataFitpPb,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendEtaRatiopPb->Draw();
// 	labelRatioMCDatapPb->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Eta_RatioMCToDataFit_pPb.%s",outputDir.Data(),suffix.Data()));
// 
// 	TF1* fitYieldDataEtapPb0020 = NULL;
// 	fitYieldDataEtapPb0020 = FitObject("l","fitYieldDataEtapPb0020","Eta",histoPCMYieldEtapPb0020,0.7,10,NULL,"QNRME+");
// 	TH1D* histoEtaRatioDatatoFitpPb0020 = CalculateHistoRatioToFit (histoPCMYieldEtapPb0020, fitYieldDataEtapPb0020);
// 	TH1D* histoEtaRatioMCtoDataFitpPb0020 = CalculateHistoRatioToFit (histoEtaInputMCpPb0020, fitYieldDataEtapPb0020);
// 	TH1D* histoEtaRatioMCUnweightedtoDataFitpPb0020 = NULL;
// 	if (histoEtaInputMCWOWeightspPb0020) histoEtaRatioMCUnweightedtoDataFitpPb0020 = CalculateHistoRatioToFit (histoEtaInputMCWOWeightspPb0020, fitYieldDataEtapPb0020);
// 	canvasFraction2->cd();
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb0020) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitpPb0020, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoEtaRatioMCtoDataFitpPb0020, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoEtaRatioDatatoFitpPb0020, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitpPb0020,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 10);
// 	histoEtaRatioDatatoFitpPb0020->Draw("same,e,p");  
// 	if (runDrawReweighted) histoEtaRatioMCtoDataFitpPb0020->Draw("same,hist,c");  
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb0020) histoEtaRatioMCUnweightedtoDataFitpPb0020->Draw("same,hist,c");  
// 	
// 	TLegend* legendEtaRatiopPb0020 = new TLegend(0.16,0.81,0.4,0.9);
// 	legendEtaRatiopPb0020->SetFillColor(0);
// 	legendEtaRatiopPb0020->SetLineColor(0);
// 	legendEtaRatiopPb0020->SetTextSize(0.025);
// //    legendEtaRatiopPb0020->SetNColumns(3);
// 	legendEtaRatiopPb0020->SetMargin(0.2);
// 	legendEtaRatiopPb0020->AddEntry(histoEtaRatioDatatoFitpPb0020,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendEtaRatiopPb0020->AddEntry(histoEtaRatioMCtoDataFitpPb0020,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb0020) legendEtaRatiopPb0020->AddEntry(histoEtaRatioMCUnweightedtoDataFitpPb0020,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendEtaRatiopPb0020->Draw();
// 	labelRatioMCDatapPb0020->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Eta_RatioMCToDataFit_pPb0020.%s",outputDir.Data(),suffix.Data()));
// 
// 	TF1* fitYieldDataEtapPb2040 = NULL;
// 	fitYieldDataEtapPb2040 = FitObject("l","fitYieldDataEtapPb2040","Eta",histoPCMYieldEtapPb2040,0.7,10,NULL,"QNRME+");
// 	TH1D* histoEtaRatioDatatoFitpPb2040 = CalculateHistoRatioToFit (histoPCMYieldEtapPb2040, fitYieldDataEtapPb2040);
// 	TH1D* histoEtaRatioMCtoDataFitpPb2040 = CalculateHistoRatioToFit (histoEtaInputMCpPb2040, fitYieldDataEtapPb2040);
// 	TH1D* histoEtaRatioMCUnweightedtoDataFitpPb2040 = NULL;
// 	if (histoEtaInputMCWOWeightspPb2040) histoEtaRatioMCUnweightedtoDataFitpPb2040 = CalculateHistoRatioToFit (histoEtaInputMCWOWeightspPb2040, fitYieldDataEtapPb2040);
// 	canvasFraction2->cd();
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb2040) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitpPb2040, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoEtaRatioMCtoDataFitpPb2040, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoEtaRatioDatatoFitpPb2040, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitpPb2040,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 10);
// 	histoEtaRatioDatatoFitpPb2040->Draw("same,e,p");  
// 	if (runDrawReweighted) histoEtaRatioMCtoDataFitpPb2040->Draw("same,hist,c");  
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb2040) histoEtaRatioMCUnweightedtoDataFitpPb2040->Draw("same,hist,c");  
// 	
// 	TLegend* legendEtaRatiopPb2040 = new TLegend(0.16,0.81,0.4,0.9);
// 	legendEtaRatiopPb2040->SetFillColor(0);
// 	legendEtaRatiopPb2040->SetLineColor(0);
// 	legendEtaRatiopPb2040->SetTextSize(0.025);
// //    legendEtaRatiopPb2040->SetNColumns(3);
// 	legendEtaRatiopPb2040->SetMargin(0.2);
// 	legendEtaRatiopPb2040->AddEntry(histoEtaRatioDatatoFitpPb2040,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendEtaRatiopPb2040->AddEntry(histoEtaRatioMCtoDataFitpPb2040,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb2040) legendEtaRatiopPb2040->AddEntry(histoEtaRatioMCUnweightedtoDataFitpPb2040,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendEtaRatiopPb2040->Draw();
// 	labelRatioMCDatapPb2040->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Eta_RatioMCToDataFit_pPb2040.%s",outputDir.Data(),suffix.Data()));
// 
// 	TF1* fitYieldDataEtapPb4060 = NULL;
// 	fitYieldDataEtapPb4060 = FitObject("l","fitYieldDataEtapPb4060","Eta",histoPCMYieldEtapPb4060,0.7,10,NULL,"QNRME+");
// 	TH1D* histoEtaRatioDatatoFitpPb4060 = CalculateHistoRatioToFit (histoPCMYieldEtapPb4060, fitYieldDataEtapPb4060);
// 	TH1D* histoEtaRatioMCtoDataFitpPb4060 = CalculateHistoRatioToFit (histoEtaInputMCpPb4060, fitYieldDataEtapPb4060);
// 	TH1D* histoEtaRatioMCUnweightedtoDataFitpPb4060 = NULL;
// 	if (histoEtaInputMCWOWeightspPb4060) histoEtaRatioMCUnweightedtoDataFitpPb4060 = CalculateHistoRatioToFit (histoEtaInputMCWOWeightspPb4060, fitYieldDataEtapPb4060);
// 	canvasFraction2->cd();
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb4060) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitpPb4060, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoEtaRatioMCtoDataFitpPb4060, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoEtaRatioDatatoFitpPb4060, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitpPb4060,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 10);
// 	histoEtaRatioDatatoFitpPb4060->Draw("same,e,p");  
// 	if (runDrawReweighted) histoEtaRatioMCtoDataFitpPb4060->Draw("same,hist,c");  
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb4060) histoEtaRatioMCUnweightedtoDataFitpPb4060->Draw("same,hist,c");  
// 	
// 	TLegend* legendEtaRatiopPb4060 = new TLegend(0.16,0.81,0.4,0.9);
// 	legendEtaRatiopPb4060->SetFillColor(0);
// 	legendEtaRatiopPb4060->SetLineColor(0);
// 	legendEtaRatiopPb4060->SetTextSize(0.025);
// //    legendEtaRatiopPb4060->SetNColumns(3);
// 	legendEtaRatiopPb4060->SetMargin(0.2);
// 	legendEtaRatiopPb4060->AddEntry(histoEtaRatioDatatoFitpPb4060,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendEtaRatiopPb4060->AddEntry(histoEtaRatioMCtoDataFitpPb4060,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb4060) legendEtaRatiopPb4060->AddEntry(histoEtaRatioMCUnweightedtoDataFitpPb4060,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendEtaRatiopPb4060->Draw();
// 	labelRatioMCDatapPb4060->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Eta_RatioMCToDataFit_pPb4060.%s",outputDir.Data(),suffix.Data()));
// 
// 	TF1* fitYieldDataEtapPb6080 = NULL;
// 	fitYieldDataEtapPb6080 = FitObject("l","fitYieldDataEtapPb6080","Eta",histoPCMYieldEtapPb6080,0.7,10,NULL,"QNRME+");
// 	TH1D* histoEtaRatioDatatoFitpPb6080 = CalculateHistoRatioToFit (histoPCMYieldEtapPb6080, fitYieldDataEtapPb6080);
// 	TH1D* histoEtaRatioMCtoDataFitpPb6080 = CalculateHistoRatioToFit (histoEtaInputMCpPb6080, fitYieldDataEtapPb6080);
// 	TH1D* histoEtaRatioMCUnweightedtoDataFitpPb6080 = NULL;
// 	if (histoEtaInputMCWOWeightspPb6080) histoEtaRatioMCUnweightedtoDataFitpPb6080 = CalculateHistoRatioToFit (histoEtaInputMCWOWeightspPb6080, fitYieldDataEtapPb6080);
// 	canvasFraction2->cd();
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb6080) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitpPb6080, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoEtaRatioMCtoDataFitpPb6080, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoEtaRatioDatatoFitpPb6080, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitpPb6080,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 10);
// 	histoEtaRatioDatatoFitpPb6080->Draw("same,e,p");  
// 	if (runDrawReweighted) histoEtaRatioMCtoDataFitpPb6080->Draw("same,hist,c");  
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb6080) histoEtaRatioMCUnweightedtoDataFitpPb6080->Draw("same,hist,c");  
// 	
// 	TLegend* legendEtaRatiopPb6080 = new TLegend(0.16,0.81,0.4,0.9);
// 	legendEtaRatiopPb6080->SetFillColor(0);
// 	legendEtaRatiopPb6080->SetLineColor(0);
// 	legendEtaRatiopPb6080->SetTextSize(0.025);
// //    legendEtaRatiopPb6080->SetNColumns(3);
// 	legendEtaRatiopPb6080->SetMargin(0.2);
// 	legendEtaRatiopPb6080->AddEntry(histoEtaRatioDatatoFitpPb6080,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendEtaRatiopPb6080->AddEntry(histoEtaRatioMCtoDataFitpPb6080,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb6080) legendEtaRatiopPb6080->AddEntry(histoEtaRatioMCUnweightedtoDataFitpPb6080,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendEtaRatiopPb6080->Draw();
// 	labelRatioMCDatapPb6080->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Eta_RatioMCToDataFit_pPb6080.%s",outputDir.Data(),suffix.Data()));
// 
// 	TF1* fitYieldDataEtapPb60100 = NULL;
// 	fitYieldDataEtapPb60100 = FitObject("l","fitYieldDataEtapPb60100","Eta",histoPCMYieldEtapPb60100,0.7,10,NULL,"QNRME+");
// 	TH1D* histoEtaRatioDatatoFitpPb60100 = CalculateHistoRatioToFit (histoPCMYieldEtapPb60100, fitYieldDataEtapPb60100);
// 	TH1D* histoEtaRatioMCtoDataFitpPb60100 = CalculateHistoRatioToFit (histoEtaInputMCpPb60100, fitYieldDataEtapPb60100);
// 	TH1D* histoEtaRatioMCUnweightedtoDataFitpPb60100 = NULL;
// 	if (histoEtaInputMCWOWeightspPb60100) histoEtaRatioMCUnweightedtoDataFitpPb60100 = CalculateHistoRatioToFit (histoEtaInputMCWOWeightspPb60100, fitYieldDataEtapPb60100);
// 	canvasFraction2->cd();
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb60100) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitpPb60100, markerStylepPb0020,markerSizepPb0020, colorCombpPb0020 , colorCombpPb0020 );
// 	DrawGammaSetMarker(histoEtaRatioMCtoDataFitpPb60100, markerStylepPb2040,markerSizepPb2040, colorCombpPb2040 , colorCombpPb2040);
// 	DrawGammaSetMarker(histoEtaRatioDatatoFitpPb60100, markerStylepPb,markerSizepPb, kBlack , kBlack);
// 	DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitpPb60100,
// 				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
// 				kFALSE, 1.5, 0, kTRUE,
// 				kTRUE, 0., 3.,
// 				kTRUE, 0., 10);
// 	histoEtaRatioDatatoFitpPb60100->Draw("same,e,p");  
// 	if (runDrawReweighted) histoEtaRatioMCtoDataFitpPb60100->Draw("same,hist,c");  
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb60100) histoEtaRatioMCUnweightedtoDataFitpPb60100->Draw("same,hist,c");  
// 	
// 	TLegend* legendEtaRatiopPb60100 = new TLegend(0.16,0.81,0.4,0.9);
// 	legendEtaRatiopPb60100->SetFillColor(0);
// 	legendEtaRatiopPb60100->SetLineColor(0);
// 	legendEtaRatiopPb60100->SetTextSize(0.025);
// //    legendEtaRatiopPb60100->SetNColumns(3);
// 	legendEtaRatiopPb60100->SetMargin(0.2);
// 	legendEtaRatiopPb60100->AddEntry(histoEtaRatioDatatoFitpPb60100,"Data/fit to Data (0.3 <pT<14)","p");
// 	if (runDrawReweighted) legendEtaRatiopPb60100->AddEntry(histoEtaRatioMCtoDataFitpPb60100,"MC weighted/fit to Data (0.3 <pT<14)","p");
// 	if (histoEtaRatioMCUnweightedtoDataFitpPb60100) legendEtaRatiopPb60100->AddEntry(histoEtaRatioMCUnweightedtoDataFitpPb60100,"MC/fit to Data (0.4 <pT<14)","p");
// 	legendEtaRatiopPb60100->Draw();
// 	labelRatioMCDatapPb60100->Draw();
// 
// 	DrawGammaLines(0., 30.,1., 1.,0.1);
// 	canvasFraction2->Update();
// 	canvasFraction2->SaveAs(Form("%s/Eta_RatioMCToDataFit_pPb60100.%s",outputDir.Data(),suffix.Data()));
// 
// 	
// 	TCanvas* canvasEffAllEta = new TCanvas("canvasEffAllEta","",200,10,1350,1350);  // gives the page size
// 	DrawGammaCanvasSettings( canvasEffAllEta, 0.1, 0.02, 0.035, 0.09);
// 		TH2F * histo2DEffEta;
// 		histo2DEffEta = new TH2F("histo2DEffEta","histo2DEffEta",1000,0,10,2000,0.e-3,3.7e-3 );
// 	//       histo2DEffEta->GetXaxis()->SetRangeUser(0.,12.);
// 		SetStyleHistoTH2ForGraphs(histo2DEffEta, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #eta}",0.03,0.04, 0.03,0.04, 1.,1.);
// 		histo2DEffEta->Draw("copy");
// 
// 		DrawGammaSetMarker(histoEtaEfficiencypPb60100, markerStylepPb60100, markerSizepPb60100, colorCombpPb60100, colorCombpPb60100); 
// 		histoEtaEfficiencypPb60100->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoEtaEfficiencypPb6080, markerStylepPb6080, markerSizepPb6080, colorCombpPb6080, colorCombpPb6080); 
// 		histoEtaEfficiencypPb6080->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoEtaEfficiencypPb2040, markerStylepPb2040, markerSizepPb2040, colorCombpPb2040, colorCombpPb2040); 
// 		histoEtaEfficiencypPb2040->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoEtaEfficiencypPb4060, markerStylepPb4060, markerSizepPb4060, colorCombpPb4060, colorCombpPb4060); 
// 		histoEtaEfficiencypPb4060->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoEtaEfficiencypPb0020, markerStylepPb0020, markerSizepPb0020, colorCombpPb0020, colorCombpPb0020); 
// 		histoEtaEfficiencypPb0020->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoEtaEfficiencypPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
// 		histoEtaEfficiencypPb->DrawCopy("e1,same");    
// 
// 		
// 		TLegend* legendEffiEta = new TLegend(0.12,0.75,0.43,0.93);
// 		legendEffiEta->SetFillColor(0);
// 		legendEffiEta->SetLineColor(0);
// 		legendEffiEta->SetTextSize(0.027);
// 		legendEffiEta->SetNColumns(1);
// 		legendEffiEta->AddEntry(histoEtaEfficiencypPb,collisionSystempPb.Data(),"p");
// 		legendEffiEta->AddEntry(histoEtaEfficiencypPb0020,collisionSystempPb0020.Data(),"p");
// 		legendEffiEta->AddEntry(histoEtaEfficiencypPb2040,collisionSystempPb2040.Data(),"p");
// 		legendEffiEta->AddEntry(histoEtaEfficiencypPb4060,collisionSystempPb4060.Data(),"p");
// 		legendEffiEta->AddEntry(histoEtaEfficiencypPb6080,collisionSystempPb6080.Data(),"p");
// 		legendEffiEta->AddEntry(histoEtaEfficiencypPb60100,collisionSystempPb60100.Data(),"p");
// 
// 		legendEffiEta->Draw();
// 	
// 	canvasEffAllEta->SaveAs(Form("%s/EfficiencyCompEta_pPb.%s",outputDir.Data(),suffix.Data()));
// 
// 	TCanvas* canvasEffAllPi0 = new TCanvas("canvasEffAllPi0","",200,10,1350,1350);  // gives the page size
// 	DrawGammaCanvasSettings( canvasEffAllPi0, 0.1, 0.02, 0.035, 0.09);
// 		TH2F * histo2DEffPi0;
// 		histo2DEffPi0 = new TH2F("histo2DEffPi0","histo2DEffPi0",1000,0,15,2000,0.e-3,3.7e-3 );
// 	//       histo2DEffPi0->GetXaxis()->SetRangeUser(0.,12.);
// 		SetStyleHistoTH2ForGraphs(histo2DEffPi0, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
// 		histo2DEffPi0->Draw("copy");
// 
// 		DrawGammaSetMarker(histoPi0EfficiencypPb60100, markerStylepPb60100, markerSizepPb60100, colorCombpPb60100, colorCombpPb60100); 
// 		histoPi0EfficiencypPb60100->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoPi0EfficiencypPb6080, markerStylepPb6080, markerSizepPb6080, colorCombpPb6080, colorCombpPb6080); 
// 		histoPi0EfficiencypPb6080->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoPi0EfficiencypPb2040, markerStylepPb2040, markerSizepPb2040, colorCombpPb2040, colorCombpPb2040); 
// 		histoPi0EfficiencypPb2040->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoPi0EfficiencypPb4060, markerStylepPb4060, markerSizepPb4060, colorCombpPb4060, colorCombpPb4060); 
// 		histoPi0EfficiencypPb4060->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoPi0EfficiencypPb0020, markerStylepPb0020, markerSizepPb0020, colorCombpPb0020, colorCombpPb0020); 
// 		histoPi0EfficiencypPb0020->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoPi0EfficiencypPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
// 		histoPi0EfficiencypPb->DrawCopy("e1,same");    
// 
// 		
// 		TLegend* legendEffiPi0 = new TLegend(0.12,0.75,0.43,0.93);
// 		legendEffiPi0->SetFillColor(0);
// 		legendEffiPi0->SetLineColor(0);
// 		legendEffiPi0->SetTextSize(0.027);
// 		legendEffiPi0->SetNColumns(1);
// 		legendEffiPi0->AddEntry(histoPi0EfficiencypPb,collisionSystempPb.Data(),"p");
// 		legendEffiPi0->AddEntry(histoPi0EfficiencypPb0020,collisionSystempPb0020.Data(),"p");
// 		legendEffiPi0->AddEntry(histoPi0EfficiencypPb2040,collisionSystempPb2040.Data(),"p");
// 		legendEffiPi0->AddEntry(histoPi0EfficiencypPb4060,collisionSystempPb4060.Data(),"p");
// 		legendEffiPi0->AddEntry(histoPi0EfficiencypPb6080,collisionSystempPb6080.Data(),"p");
// 		legendEffiPi0->AddEntry(histoPi0EfficiencypPb60100,collisionSystempPb60100.Data(),"p");
// 
// 		legendEffiPi0->Draw();
// 	
// 	canvasEffAllPi0->SaveAs(Form("%s/EfficiencyCompPi0_pPb.%s",outputDir.Data(),suffix.Data()));
// 
// 	
// 	canvasEffAllPi0->cd();
// 		histo2DEffPi0->GetYaxis()->SetTitle("#epsilon_{reco}");
// 		histo2DEffPi0->GetYaxis()->SetRangeUser(0,2.5e-3);
// 		histo2DEffPi0->Draw("copy");
// 
// 		DrawGammaSetMarker(histoPi0EfficiencypPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
// 		histoPi0EfficiencypPb->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoEtaEfficiencypPb, markerStylepPb+4, markerSizepPb, kGray+1, kGray+1); 
// 		histoEtaEfficiencypPb->DrawCopy("e1,same");    
// 
// 		
// 		TLegend* legendEffiPi0Eta = new TLegend(0.12,0.85,0.43,0.93);
// 		legendEffiPi0Eta->SetFillColor(0);
// 		legendEffiPi0Eta->SetLineColor(0);
// 		legendEffiPi0Eta->SetTextSize(0.027);
// 		legendEffiPi0Eta->SetNColumns(1);
// 		legendEffiPi0Eta->AddEntry(histoPi0EfficiencypPb,"#pi^{0}","p");
// 		legendEffiPi0Eta->AddEntry(histoEtaEfficiencypPb,"#eta","p");
// 		legendEffiPi0Eta->Draw();
// 	
// 	canvasEffAllPi0->SaveAs(Form("%s/EfficiencyPi0AndEta_pPbMB.%s",outputDir.Data(),suffix.Data()));
// 	
// 	canvasEffAllPi0->cd();
// 		histo2DEffPi0->GetYaxis()->SetTitle("#epsilon_{reco, #pi^{0}}");
// 		histo2DEffPi0->GetYaxis()->SetRangeUser(0,2.5e-3);
// 		histo2DEffPi0->Draw("copy");
// 
// 		DrawGammaSetMarker(histoPi0EfficiencypPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
// 		histoPi0EfficiencypPb->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoPi0EfficiencypPb_DPMJET, markerStylepPb+4, markerSizepPb, kRed-7, kRed-7); 
// 		histoPi0EfficiencypPb_DPMJET->DrawCopy("e1,same");    
// 
// 		
// 		TLegend* legendEffiPi0DiffMC = new TLegend(0.12,0.85,0.43,0.93);
// 		legendEffiPi0DiffMC->SetFillColor(0);
// 		legendEffiPi0DiffMC->SetLineColor(0);
// 		legendEffiPi0DiffMC->SetTextSize(0.027);
// 		legendEffiPi0DiffMC->SetNColumns(1);
// 		legendEffiPi0DiffMC->AddEntry(histoPi0EfficiencypPb,"HIJING","p");
// 		legendEffiPi0DiffMC->AddEntry(histoPi0EfficiencypPb_DPMJET,"DPMJET","p");
// 		legendEffiPi0DiffMC->Draw();
// 	
// 	canvasEffAllPi0->SaveAs(Form("%s/EfficiencyPi0_DiffMC_pPbMB.%s",outputDir.Data(),suffix.Data()));
// 	
// 	canvasEffAllPi0->cd();
// 		histo2DEffPi0->GetYaxis()->SetTitle("#epsilon_{reco, #eta}");
// 		histo2DEffPi0->GetYaxis()->SetRangeUser(0,2.5e-3);
// 		histo2DEffPi0->Draw("copy");
// 
// 		DrawGammaSetMarker(histoEtaEfficiencypPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
// 		histoEtaEfficiencypPb->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoEtaEfficiencypPb_DPMJET, markerStylepPb+4, markerSizepPb, kRed-7, kRed-7); 
// 		histoEtaEfficiencypPb_DPMJET->DrawCopy("e1,same");    
// 
// 		
// 		TLegend* legendEffiEtaDiffMC = new TLegend(0.12,0.85,0.43,0.93);
// 		legendEffiEtaDiffMC->SetFillColor(0);
// 		legendEffiEtaDiffMC->SetLineColor(0);
// 		legendEffiEtaDiffMC->SetTextSize(0.027);
// 		legendEffiEtaDiffMC->SetNColumns(1);
// 		legendEffiEtaDiffMC->AddEntry(histoEtaEfficiencypPb,"HIJING","p");
// 		legendEffiEtaDiffMC->AddEntry(histoEtaEfficiencypPb_DPMJET,"DPMJET","p");
// 		legendEffiEtaDiffMC->Draw();
// 	
// 	canvasEffAllPi0->SaveAs(Form("%s/EfficiencyEta_DiffMC_pPbMB.%s",outputDir.Data(),suffix.Data()));
// 	
// 	TCanvas* canvasAccepAllCent = new TCanvas("canvasAccepAllCent","",200,10,1350,1350);  // gives the page size
// 	DrawGammaCanvasSettings( canvasAccepAllCent, 0.1, 0.02, 0.035, 0.09);
// 	TH2F * histo2DAccepAllEta;
// 	histo2DAccepAllEta = new TH2F("histo2DAccepAllEta","histo2DAccepAllEta",1000,0,15.,2000,0.4,1.02 );
// 	histo2DAccepAllEta->GetXaxis()->SetRangeUser(0.,15.);
// 	SetStyleHistoTH2ForGraphs(histo2DAccepAllEta, "#it{p}_{T} (GeV/#it{c})","A",0.03,0.04, 0.03,0.04, 1.,1.2);
// 	histo2DAccepAllEta->Draw("copy");
// 		
// 		DrawGammaSetMarker(histoPi0AcceptancepPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
// 		histoPi0AcceptancepPb->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoEtaAcceptancepPb, markerStylepPb+4, markerSizepPb, kGray+1, kGray+1); 
// 		histoEtaAcceptancepPb->DrawCopy("e1,same");    
// 
// 		TLegend* legendAcceptAllEtaPi0 = new TLegend(0.34,0.13,0.63,0.23);
// 		legendAcceptAllEtaPi0->SetFillColor(0);
// 		legendAcceptAllEtaPi0->SetLineColor(0);
// 		legendAcceptAllEtaPi0->SetTextSize(0.03);
// 		legendAcceptAllEtaPi0->SetNColumns(2);
// 		legendAcceptAllEtaPi0->AddEntry(histoPi0AcceptancepPb,"#pi^{0}","p");
// 		legendAcceptAllEtaPi0->AddEntry((TObject*)0, " |y| < 0.8","");
// 		legendAcceptAllEtaPi0->AddEntry(histoEtaAcceptancepPb,"#eta","p");
// 		legendAcceptAllEtaPi0->AddEntry((TObject*)0, " |y| < 0.8","");
// 		legendAcceptAllEtaPi0->Draw();
// 	
// 	canvasAccepAllCent->SaveAs(Form("%s/AcceptanceEtaAndPi0_pPbMB.%s",outputDir.Data(),suffix.Data()));
// 	
// 	canvasAccepAllCent->cd();
// 	histo2DAccepAllEta->GetYaxis()->SetTitle("A_{#pi^{0}}");
// 	histo2DAccepAllEta->Draw("copy");
// 		
// 		DrawGammaSetMarker(histoPi0AcceptancepPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
// 		histoPi0AcceptancepPb->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoPi0AcceptancepPb_DPMJET, markerStylepPb+4, markerSizepPb, kRed-7, kRed-7); 
// 		histoPi0AcceptancepPb_DPMJET->DrawCopy("e1,same");    
// 
// 		TLegend* legendAcceptAllPi0DiffMC = new TLegend(0.34,0.13,0.63,0.23);
// 		legendAcceptAllPi0DiffMC->SetFillColor(0);
// 		legendAcceptAllPi0DiffMC->SetLineColor(0);
// 		legendAcceptAllPi0DiffMC->SetTextSize(0.03);
// 		legendAcceptAllPi0DiffMC->SetNColumns(2);
// 		legendAcceptAllPi0DiffMC->AddEntry(histoPi0AcceptancepPb,"HIJING","p");
// 		legendAcceptAllPi0DiffMC->AddEntry((TObject*)0, " |y| < 0.8","");
// 		legendAcceptAllPi0DiffMC->AddEntry(histoPi0AcceptancepPb_DPMJET,"DPMJET","p");
// 		legendAcceptAllPi0DiffMC->AddEntry((TObject*)0, " |y| < 0.8","");
// 		legendAcceptAllPi0DiffMC->Draw();
// 	
// 	canvasAccepAllCent->SaveAs(Form("%s/AcceptancePi0_DiffMC_pPbMB.%s",outputDir.Data(),suffix.Data()));
// 
// 	canvasAccepAllCent->cd();
// 	histo2DAccepAllEta->GetYaxis()->SetTitle("A_{#eta}");
// 	histo2DAccepAllEta->Draw("copy");
// 		
// 		DrawGammaSetMarker(histoEtaAcceptancepPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
// 		histoEtaAcceptancepPb->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoEtaAcceptancepPb_DPMJET, markerStylepPb+4, markerSizepPb, kRed-7, kRed-7); 
// 		histoEtaAcceptancepPb_DPMJET->DrawCopy("e1,same");    
// 
// 		TLegend* legendAcceptAllEtaDiffMC = new TLegend(0.34,0.13,0.63,0.23);
// 		legendAcceptAllEtaDiffMC->SetFillColor(0);
// 		legendAcceptAllEtaDiffMC->SetLineColor(0);
// 		legendAcceptAllEtaDiffMC->SetTextSize(0.03);
// 		legendAcceptAllEtaDiffMC->SetNColumns(2);
// 		legendAcceptAllEtaDiffMC->AddEntry(histoEtaAcceptancepPb,"HIJING","p");
// 		legendAcceptAllEtaDiffMC->AddEntry((TObject*)0, " |y| < 0.8","");
// 		legendAcceptAllEtaDiffMC->AddEntry(histoEtaAcceptancepPb_DPMJET,"DPMJET","p");
// 		legendAcceptAllEtaDiffMC->AddEntry((TObject*)0, " |y| < 0.8","");
// 		legendAcceptAllEtaDiffMC->Draw();
// 	
// 	canvasAccepAllCent->SaveAs(Form("%s/AcceptanceEta_DiffMC_pPbMB.%s",outputDir.Data(),suffix.Data()));
// 
// 	delete canvasAccepAllCent;

	TCanvas* canvasGammaRecEff = new TCanvas("canvasGammaRecEff","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasGammaRecEff, 0.1, 0.02, 0.035, 0.09);
	TH2F * histo2DRecEffGamma;
	histo2DRecEffGamma = new TH2F("histo2DRecEffGamma","histo2DRecEffGamma",1000,0,15.,2000,0.2,0.8 );
	histo2DRecEffGamma->GetXaxis()->SetRangeUser(0.,15.);
	SetStyleHistoTH2ForGraphs(histo2DRecEffGamma, "#it{p}_{T} (GeV/#it{c})","#epsilon_{rec, #gamma}",0.03,0.04, 0.03,0.04, 1.,1.2);
	histo2DRecEffGamma->Draw("copy");
		
		DrawGammaSetMarker(histoGammaRecoEffpPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
		histoGammaRecoEffpPb->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoGammaRecoEffpPb_DPMJET, markerStylepPb+4, markerSizepPb, kRed-7, kRed-7); 
// 		histoGammaRecoEffpPb_DPMJET->DrawCopy("e1,same");    

		TLegend* legendRecEffGammaDiffMC = new TLegend(0.34,0.13,0.63,0.18);
		legendRecEffGammaDiffMC->SetFillColor(0);
		legendRecEffGammaDiffMC->SetLineColor(0);
		legendRecEffGammaDiffMC->SetTextSize(0.03);
		legendRecEffGammaDiffMC->SetNColumns(2);
		legendRecEffGammaDiffMC->AddEntry(histoGammaRecoEffpPb,"HIJING","p");
		legendRecEffGammaDiffMC->AddEntry((TObject*)0, " |#eta| < 0.9","");
// 		legendRecEffGammaDiffMC->AddEntry(histoGammaRecoEffpPb_DPMJET,"DPMJET","p");
// 		legendRecEffGammaDiffMC->AddEntry((TObject*)0, " |#eta| < 0.9","");
		legendRecEffGammaDiffMC->Draw();
	
	canvasGammaRecEff->SaveAs(Form("%s/PhotonRecEff_DiffMC_pPbMB.%s",outputDir.Data(),suffix.Data()));

	TCanvas* canvasGammaPurity = new TCanvas("canvasGammaPurity","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasGammaPurity, 0.1, 0.02, 0.035, 0.09);
	TH2F * histo2DPurityGamma;
	histo2DPurityGamma = new TH2F("histo2DPurityGamma","histo2DPurityGamma",1000,0,15.,2000,0.88,1.01 );
	histo2DPurityGamma->GetXaxis()->SetRangeUser(0.,15.);
	SetStyleHistoTH2ForGraphs(histo2DPurityGamma, "#it{p}_{T} (GeV/#it{c})","#epsilon_{pur, #gamma}",0.03,0.04, 0.03,0.04, 1.,1.2);
	histo2DPurityGamma->Draw("copy");
		
		DrawGammaSetMarker(histoGammaPuritypPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
		histoGammaPuritypPb->DrawCopy("e1,same");    

		TLegend* legendPurityGammaDiffMC = new TLegend(0.34,0.13,0.63,0.18);
		legendPurityGammaDiffMC->SetFillColor(0);
		legendPurityGammaDiffMC->SetLineColor(0);
		legendPurityGammaDiffMC->SetTextSize(0.03);
		legendPurityGammaDiffMC->SetNColumns(2);
		legendPurityGammaDiffMC->AddEntry(histoGammaPuritypPb,"HIJING","p");
		legendPurityGammaDiffMC->AddEntry((TObject*)0, " |#eta| < 0.9","");
		legendPurityGammaDiffMC->Draw();
	
	canvasGammaPurity->SaveAs(Form("%s/PhotonPurity_DiffMC_pPbMB.%s",outputDir.Data(),suffix.Data()));
	
	TCanvas* canvasGammaConvProb = new TCanvas("canvasGammaConvProb","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasGammaConvProb, 0.1, 0.02, 0.035, 0.09);
	TH2F * histo2DConvProbGamma;
	histo2DConvProbGamma = new TH2F("histo2DConvProbGamma","histo2DConvProbGamma",1000,0,15.,2000,0.,0.1 );
	histo2DConvProbGamma->GetXaxis()->SetRangeUser(0.,15.);
	SetStyleHistoTH2ForGraphs(histo2DConvProbGamma, "#it{p}_{T} (GeV/#it{c})","P_{conv, #gamma}",0.03,0.04, 0.03,0.04, 1.,1.2);
	histo2DConvProbGamma->Draw("copy");
		
		DrawGammaSetMarker(histoGammaConvProbpPb, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb); 
		histoGammaConvProbpPb->DrawCopy("e1,same");    
// 		DrawGammaSetMarker(histoGammaConvProbpPb_DPMJET, markerStylepPb+4, markerSizepPb, kRed-7, kRed-7); 
// 		histoGammaConvProbpPb_DPMJET->DrawCopy("e1,same");    
		DrawGammaSetMarker(histoGammaConvProbpp7TeV, markerStylepPb+4, markerSizepPb, kRed-7, kRed-7); 
		histoGammaConvProbpp7TeV->DrawCopy("e1,same");    

		
		TLegend* legendConvProbGammaDiffMC = new TLegend(0.34,0.13,0.93,0.18);
		legendConvProbGammaDiffMC->SetFillColor(0);
		legendConvProbGammaDiffMC->SetLineColor(0);
		legendConvProbGammaDiffMC->SetTextSize(0.03);
		legendConvProbGammaDiffMC->SetNColumns(2);
		legendConvProbGammaDiffMC->AddEntry(histoGammaConvProbpPb,"HIJING","p");
		legendConvProbGammaDiffMC->AddEntry((TObject*)0, " |#eta| < 0.9","");
		legendConvProbGammaDiffMC->AddEntry(histoGammaConvProbpp7TeV,"Pythia 6 + PHOJET, pp 7 TeV","p");
		legendConvProbGammaDiffMC->AddEntry((TObject*)0, " |#eta| < 0.9","");
		legendConvProbGammaDiffMC->Draw();
	
	canvasGammaConvProb->SaveAs(Form("%s/PhotonConvProb_DiffMC_pPbMB.%s",outputDir.Data(),suffix.Data()));
	
// 	TFile fMCSpectraInput("MCSpectraInputpPb.root","UPDATE");
// 		if (fitYieldDataPi0pPb){
// 			fitYieldDataPi0pPb->SetRange(0,30);
// 			fitYieldDataPi0pPb->Write("Pi0_Fit_Data_pPb_5023GeV_MBV0A",TObject::kOverwrite);
// 		}
// 		if (fitYieldDataPi0pPb0020){
// 			fitYieldDataPi0pPb0020->SetRange(0,30);
// 			fitYieldDataPi0pPb0020->Write("Pi0_Fit_Data_pPb_5023GeV_0020V0A",TObject::kOverwrite);
// 		}   
// 		if (fitYieldDataPi0pPb2040){
// 			fitYieldDataPi0pPb2040->SetRange(0,30);
// 			fitYieldDataPi0pPb2040->Write("Pi0_Fit_Data_pPb_5023GeV_2040V0A",TObject::kOverwrite);
// 		}
// 		if (fitYieldDataPi0pPb4060){
// 			fitYieldDataPi0pPb4060->SetRange(0,30);
// 			fitYieldDataPi0pPb4060->Write("Pi0_Fit_Data_pPb_5023GeV_4060V0A",TObject::kOverwrite);
// 		}
// 		if (fitYieldDataPi0pPb6080){
// 			fitYieldDataPi0pPb6080->SetRange(0,30);
// 			fitYieldDataPi0pPb6080->Write("Pi0_Fit_Data_pPb_5023GeV_6080V0A",TObject::kOverwrite);
// 		}
// 		if (fitYieldDataPi0pPb60100){
// 			fitYieldDataPi0pPb60100->SetRange(0,30);
// 			fitYieldDataPi0pPb60100->Write("Pi0_Fit_Data_pPb_5023GeV_60100V0A",TObject::kOverwrite);
// 		}
// 		if (fitYieldDataEtapPb){
// 			fitYieldDataEtapPb->SetRange(0,30);
// 			fitYieldDataEtapPb->Write("Eta_Fit_Data_pPb_5023GeV_MBV0A",TObject::kOverwrite);
// 		}
// 		if (fitYieldDataEtapPb0020){
// 			fitYieldDataEtapPb0020->SetRange(0,30);
// 			fitYieldDataEtapPb0020->Write("Eta_Fit_Data_pPb_5023GeV_0020V0A",TObject::kOverwrite);
// 		}
// 		if (fitYieldDataEtapPb2040){
// 			fitYieldDataEtapPb2040->SetRange(0,30);
// 			fitYieldDataEtapPb2040->Write("Eta_Fit_Data_pPb_5023GeV_2040V0A",TObject::kOverwrite);
// 		}
// 		if (fitYieldDataEtapPb4060){
// 			fitYieldDataEtapPb4060->SetRange(0,30);
// 			fitYieldDataEtapPb4060->Write("Eta_Fit_Data_pPb_5023GeV_4060V0A",TObject::kOverwrite);
// 		}
// 		if (fitYieldDataEtapPb6080){
// 			fitYieldDataEtapPb6080->SetRange(0,30);
// 			fitYieldDataEtapPb6080->Write("Eta_Fit_Data_pPb_5023GeV_6080V0A",TObject::kOverwrite);
// 		}   
// 		if (fitYieldDataEtapPb60100){
// 			fitYieldDataEtapPb60100->SetRange(0,30);
// 			fitYieldDataEtapPb60100->Write("Eta_Fit_Data_pPb_5023GeV_60100V0A",TObject::kOverwrite);
// 		}   
// 		histoPi0InputMCWOWeightsAddedSigpPb->Write("Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightspPb->Write("Pi0_Hijing_LHC13e7_pPb_5023GeV_MBV0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightsAddedSigpPb0020->Write("Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_0020V0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightspPb0020->Write("Pi0_Hijing_LHC13e7_pPb_5023GeV_0020V0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightsAddedSigpPb2040->Write("Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_2040V0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightspPb2040->Write("Pi0_Hijing_LHC13e7_pPb_5023GeV_2040V0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightsAddedSigpPb4060->Write("Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_4060V0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightspPb4060->Write("Pi0_Hijing_LHC13e7_pPb_5023GeV_4060V0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightsAddedSigpPb6080->Write("Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_6080V0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightspPb6080->Write("Pi0_Hijing_LHC13e7_pPb_5023GeV_6080V0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightsAddedSigpPb60100->Write("Pi0_Hijing_LHC13e7_addSig_pPb_5023GeV_60100V0A",TObject::kOverwrite);
// 		histoPi0InputMCWOWeightspPb60100->Write("Pi0_Hijing_LHC13e7_pPb_5023GeV_60100V0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightsAddedSigpPb->Write("Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_MBV0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightspPb->Write("Eta_Hijing_LHC13e7_pPb_5023GeV_MBV0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightsAddedSigpPb0020->Write("Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_0020V0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightspPb0020->Write("Eta_Hijing_LHC13e7_pPb_5023GeV_0020V0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightsAddedSigpPb2040->Write("Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_2040V0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightspPb2040->Write("Eta_Hijing_LHC13e7_pPb_5023GeV_2040V0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightsAddedSigpPb4060->Write("Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_4060V0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightspPb4060->Write("Eta_Hijing_LHC13e7_pPb_5023GeV_4060V0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightsAddedSigpPb6080->Write("Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_6080V0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightspPb6080->Write("Eta_Hijing_LHC13e7_pPb_5023GeV_6080V0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightsAddedSigpPb60100->Write("Eta_Hijing_LHC13e7_addSig_pPb_5023GeV_60100V0A",TObject::kOverwrite);
// 		histoEtaInputMCWOWeightspPb60100->Write("Eta_Hijing_LHC13e7_pPb_5023GeV_60100V0A",TObject::kOverwrite);
// 	fMCSpectraInput.Close();
		
}

