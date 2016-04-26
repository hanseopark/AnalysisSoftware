/****************************************************************************************************************************
****** 		provided by Gamma Conversion Group, PWGGA, 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
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

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

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


void CombineNeutralPionResults(TString outputDir = "pdf/EfficiencyComparison", TString suffix = "pdf", TString nameFilePP = "data_PCMResultsFullCorrection_PP_NoBinShifting.root", TString nameFilePbPb = "data_GammaConversionResultsFullCorrection_PbPb_2.76TeV_.root", Bool_t runDrawReweighted = kTRUE){

	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	gSystem->Exec("mkdir -p "+outputDir);
	
	StyleSettingsThesis();	
	SetPlotStyle();

	Color_t	colorCombPP					= kBlack;
	Color_t	colorCombPbPb0005				= kRed+1;
   Color_t  colorCombPbPb0010          = kRed+1;
	Color_t	colorCombPbPb0510				= 807;
	Color_t	colorCombPbPb1020				= 800;
	Color_t	colorCombPbPb2040				= kGreen+2;
	Color_t	colorCombPbPb4060				= kCyan+2;
	Color_t	colorCombPbPb6080				= kBlue+1;

	Color_t	colorCombMCPbPb0005				= kRed+3;
   Color_t  colorCombMCPbPb0010           = kRed+3;
	Color_t	colorCombMCPbPb0510				= 807+2;
	Color_t	colorCombMCPbPb1020				= 800+2;
	Color_t	colorCombMCPbPb2040				= kGreen+4;
	Color_t	colorCombMCPbPb4060				= kCyan+4;
	Color_t	colorCombMCPbPb6080				= kBlue+3;

	Style_t 	markerStylePP	 	= 33 ;
	Style_t 	markerStylePbPb0005 	= 20 ;
	Style_t 	markerStylePbPb0010 	= 20 ;
   Style_t  markerStylePbPb0510  = 21 ;
	Style_t 	markerStylePbPb1020 	= 29 ;
	Style_t 	markerStylePbPb2040 	= 33 ;
	Style_t 	markerStylePbPb4060 	= 20 ;
	Style_t 	markerStylePbPb6080 	= 21 ;

	Style_t 	markerStylePbPb0005MC 	= 24 ;
	Style_t 	markerStylePbPb0010MC 	= 24 ;
   Style_t  markerStylePbPb0510MC   = 25 ;
	Style_t 	markerStylePbPb1020MC 	= 30 ;
	Style_t 	markerStylePbPb2040MC 	= 27 ;
	Style_t 	markerStylePbPb4060MC 	= 24 ;
	Style_t 	markerStylePbPb6080MC 	= 25 ;
   
	Size_t 	markerSizePP	 	= 2.5;
	Size_t 	markerSizePbPb0005 	= 2.;
   Size_t   markerSizePbPb0010   = 2.;
	Size_t 	markerSizePbPb0510 	= 2.;
	Size_t 	markerSizePbPb1020 	= 2.5;
	Size_t 	markerSizePbPb2040 	= 2.5;
	Size_t 	markerSizePbPb4060 	= 2.;
	Size_t 	markerSizePbPb6080 	= 2.;
	
	Color_t 	colorPi0900GeV 			= kRed +2;
	Color_t 	colorPi02760GeV 			= kMagenta+2;
	Color_t	colorPi07TeV				= kBlue+2;
   Color_t  colorPi0900GeVBox = colorPi0900GeV-10;
   Color_t  colorPi02760GeVBox = colorPi02760GeV-10;
   Color_t  colorPi07TeVBox = colorPi07TeV-10;

	Color_t 	colorMCPythiaPP900GeV 	= colorPi0900GeV-4;
	Color_t 	colorMCPythiaPP2760GeV = colorPi02760GeV+2;
	Color_t 	colorMCPythiaPP7TeV 	= colorPi07TeV+3;
	Color_t 	colorMCPhojetPP900GeV 	= colorPi0900GeV+2;
	Color_t 	colorMCPhojetPP2760GeV = colorPi02760GeV-4;
	Color_t 	colorMCPhojetPP7TeV 	= colorPi07TeV-3;

	Style_t 	markerStyleSpectrum7TeVMC 	= 24 ;
	Style_t 	markerStyleSpectrum900GeVMC = 25 ;
	Style_t 	markerStyleSpectrum2760GeVMC = 30 ;
	Style_t 	markerStyleSpectrum7TeV 	= 20 ;
	Style_t 	markerStyleSpectrum900GeV = 21 ;
	Style_t 	markerStyleSpectrum2760GeV = 29 ;
	
	Double_t xSection7TeVppINEL = 73.2*1e9;
	Double_t xSection2760GeVppINEL = 62.8*1e9;
	Double_t xSection900GeVppINEL = 52.5*1e9;
	
	Style_t 	markerStyleMCPP7TeV 	= 24 ;
	Style_t 	markerStyleMCPP900GeV 	= 25 ;
	Style_t 	markerStyleMCPP2760GeV 	= 30 ;
	
	Size_t 	markerSizePi0PP7TeV 	= 1.8;
	Size_t 	markerSizePi0PP900GeV = 1.8;
	Size_t 	markerSizePi0PP2760GeV 	= 2.2;

   Color_t  colorCombpPb          = kBlack;
   Color_t  colorCombpPb0020          = kRed+1;
   Color_t  colorCombpPb2040          = kGreen+2;
   Color_t  colorCombpPb4060          = kCyan+2;
   Color_t  colorCombpPb6080          = kBlue+1;
   Color_t  colorCombpPb80100         = kViolet+3;
   Color_t  colorCombpPb60100         = kBlue+1;

   Color_t  colorCombpPbBox       = kGray+1;
   Color_t  colorCombpPb0020Box       = kRed-6;
   Color_t  colorCombpPb2040Box       = kGreen-6;
   Color_t  colorCombpPb4060Box       = kCyan-6;
   Color_t  colorCombpPb6080Box       = kBlue-6;
   Color_t  colorCombpPb80100Box      = kViolet+6;
   Color_t  colorCombpPb60100Box      = kBlue-6;

   Style_t  markerStylepPb  = 20 ;
   Style_t  markerStylepPb0020   = 20 ;
   Style_t  markerStylepPb2040   = 21 ;
   Style_t  markerStylepPb4060   = 29 ;
   Style_t  markerStylepPb6080   = 33 ;
   Style_t  markerStylepPb80100 = 34 ;
   Style_t  markerStylepPb60100 = 33 ;

   Size_t   markerSizepPb  = 2.;
   Size_t   markerSizepPb0020  = 2.;
   Size_t   markerSizepPb2040  = 2.;
   Size_t   markerSizepPb4060  = 2.5;
   Size_t   markerSizepPb6080  = 2.5;
   Size_t   markerSizepPb80100  = 2.;
   Size_t   markerSizepPb60100  = 2.5;

   
	TString collisionSystemPP2760GeV = "pp #sqrt{#it{s}} = 2.76 TeV";		
	TString collisionSystemPP7TeV = "pp #sqrt{#it{s}} = 7 TeV";		
	TString collisionSystemPP900GeV = "pp #sqrt{#it{s}} = 0.9 TeV";		
	TString collisionSystemPbPb0005 = "0-5% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemPbPb0510 = "5-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemPbPb1020 = "10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemPbPb0010 = "0-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemPbPb2040 = "20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemPbPb4060 = "40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemPbPb6080 = "60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
   TString collisionSystemPbPb0020 = "0-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";     
   TString collisionSystemPbPb0080 = "0-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";     
   TString collisionSystemPbPb0040 = "0-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";     
   TString collisionSystemPbPb4080 = "40-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";     
   TString collisionSystempPb0020 = "0-20% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   TString collisionSystempPb2040 = "20-40% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";      
   TString collisionSystempPb4060 = "40-60% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";    
   TString collisionSystempPb60100 = "60-100% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   // TString collisionSystempPb6080 = "60-80% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   // TString collisionSystempPb80100 = "80-100% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   TString collisionSystempPb = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   

	Size_t markerSizeComparison = 0.5;
	Double_t maxPtMesonEffFit = 12.;
	Double_t minPtMesonEffFit = 1.2;
	Int_t offsetCorrectionHighPt= 1;
	TF1* fitTrueEffi = new TF1("EffiFitDummy","1 - [0]*exp([1]*x)+[1]");
	fitTrueEffi->SetRange(minPtMesonEffFit,maxPtMesonEffFit);


	TFile* fileNeutralPionDataPP = new TFile(nameFilePP.Data());
	
	TDirectory*	directoryPi07TeV = 		(TDirectory*)fileNeutralPionDataPP->Get("Pi07TeV"); 
	TH1D* histoAccPi07TeV = 				(TH1D*)directoryPi07TeV->Get("AcceptancePi0");
	TH1D* histoTrueEffPtPi07TeV = 		(TH1D*)directoryPi07TeV->Get("EfficiencyPi0");
	TH1D* histoRawYieldPi07TeV = 			(TH1D*)directoryPi07TeV->Get("RAWYieldPerEventsPi0");
	TDirectory*	directoryEta7TeV = 		(TDirectory*)fileNeutralPionDataPP->Get("Eta7TeV");
   TH1D* histoCorrectedYieldEta7TeV = (TH1D*)directoryEta7TeV->Get("CorrectedYieldEta");
   TGraphAsymmErrors* graphCorrectedYieldSysEta7TeV =   (TGraphAsymmErrors*)directoryEta7TeV->Get("EtaSystError");
	TH1D* histoAccEta7TeV = 				(TH1D*)directoryEta7TeV->Get("AcceptanceEta");
	TH1D* histoTrueEffPtEta7TeV = 		(TH1D*)directoryEta7TeV->Get("EfficiencyEta");
	TH1D* histoRawYieldEta7TeV = 			(TH1D*)directoryEta7TeV->Get("RAWYieldPerEventsEta");
	TH1D*	histoEtaMassData7TeV = 			(TH1D*)directoryEta7TeV->Get("MassEta");
	TH1D*	histoEtaWidthData7TeV = 		(TH1D*)directoryEta7TeV->Get("FWHMEtaMeV");
	TH1D*	histoEtaMassMC7TeV = 			(TH1D*)directoryEta7TeV->Get("TrueMassEta");
	TH1D*	histoEtaWidthMC7TeV = 			(TH1D*)directoryEta7TeV->Get("TrueFWHMEtaMeV");
   TH1D* histoRatioEtaPi07TeV=        (TH1D*)directoryEta7TeV->Get("EtatoPi0RatioConversion");   
   TGraphAsymmErrors* graphRatioEtaPi0SystErr7TeV=             (TGraphAsymmErrors*)directoryEta7TeV->Get("EtatoPi0RatioConversionSys");

	histoEtaMassData7TeV->Scale(1000.);
	histoEtaMassMC7TeV->Scale(1000.);
	histoEtaMassData7TeV->SetBinContent(histoEtaMassData7TeV->GetNbinsX(),0);
	histoEtaMassMC7TeV->SetBinContent(histoEtaMassMC7TeV->GetNbinsX(),0);
	histoEtaWidthData7TeV->SetBinContent(histoEtaWidthData7TeV->GetNbinsX(),10000.);
	histoEtaWidthMC7TeV->SetBinContent(histoEtaWidthMC7TeV->GetNbinsX(),10000.);

	TDirectory*	directoryPi0900GeV = 	(TDirectory*)fileNeutralPionDataPP->Get("Pi0900GeV"); 
	TH1D* histoAccPi0900GeV = 				(TH1D*)directoryPi0900GeV->Get("AcceptancePi0");
	TH1D* histoTrueEffPtPi0900GeV = 	(TH1D*)directoryPi0900GeV->Get("EfficiencyPi0");
	TH1D* histoRawYieldPi0900GeV = 		(TH1D*)directoryPi0900GeV->Get("RAWYieldPerEventsPi0");
	TDirectory*	directoryEta900GeV = 	(TDirectory*)fileNeutralPionDataPP->Get("Eta900GeV"); 
   TH1D* histoCorrectedYieldEta900GeV =     (TH1D*)directoryEta900GeV->Get("CorrectedYieldEta");
   TGraphAsymmErrors* graphCorrectedYieldSysEta900GeV =  (TGraphAsymmErrors*)directoryEta900GeV->Get("EtaSystError");
	TH1D* histoAccEta900GeV = 				(TH1D*)directoryEta900GeV->Get("AcceptanceEta");
	TH1D* histoTrueEffPtEta900GeV = 	(TH1D*)directoryEta900GeV->Get("EfficiencyEta");
	TH1D* histoRawYieldEta900GeV = 		(TH1D*)directoryEta900GeV->Get("RAWYieldPerEventsEta");
   TH1D* histoRatioEtaPi0900GeV=      (TH1D*)directoryEta900GeV->Get("EtatoPi0RatioConversion");   
   TGraphAsymmErrors* graphRatioEtaPi0SystErr900GeV=           (TGraphAsymmErrors*)directoryEta900GeV->Get("EtatoPi0RatioConversionSys");

	TDirectory*	directoryPi02760GeV = 	(TDirectory*)fileNeutralPionDataPP->Get("Pi02.76TeV"); 
	TH1D* histoAccPi02760GeV = 			(TH1D*)directoryPi02760GeV->Get("AcceptancePi0");
	TH1D* histoTrueEffPtPi02760GeV = 	(TH1D*)directoryPi02760GeV->Get("EfficiencyPi0");
	TH1D* histoRawYieldPi02760GeV = 	(TH1D*)directoryPi02760GeV->Get("RAWYieldPerEventsPi0");
	TH1D* histoCorrectedYieldPi02760GeV = 	(TH1D*)directoryPi02760GeV->Get("CorrectedYieldPi0");
   TGraphAsymmErrors* graphCorrectedYieldSysPi02760GeV = (TGraphAsymmErrors*)directoryPi02760GeV->Get("Pi0SystError"); 
	TDirectory*	directoryEta2760GeV = 	(TDirectory*)fileNeutralPionDataPP->Get("Eta2.76TeV"); 
	TH1D* histoAccEta2760GeV = 			(TH1D*)directoryEta2760GeV->Get("AcceptanceEta");
	TH1D* histoTrueEffPtEta2760GeV = 	(TH1D*)directoryEta2760GeV->Get("EfficiencyEta");
	TH1D* histoCorrectedYieldEta2760GeV = 	(TH1D*)directoryEta2760GeV->Get("CorrectedYieldEta");
	TGraphAsymmErrors* graphCorrectedYieldSysEta2760GeV = (TGraphAsymmErrors*)directoryEta2760GeV->Get("EtaSystError");
	TH1D* histoRawYieldEta2760GeV = 	(TH1D*)directoryEta2760GeV->Get("RAWYieldPerEventsEta");
	TH1D*	histoEtaMassData2760GeV = 		(TH1D*)directoryEta2760GeV->Get("MassEta");
	TH1D*	histoEtaWidthData2760GeV = 	(TH1D*)directoryEta2760GeV->Get("FWHMEtaMeV");
	TH1D*	histoEtaMassMC2760GeV = 		(TH1D*)directoryEta2760GeV->Get("TrueMassEta");
	TH1D*	histoEtaWidthMC2760GeV = 		(TH1D*)directoryEta2760GeV->Get("TrueFWHMEtaMeV");
   TH1D* histoRatioEtaPi02760GeV=     (TH1D*)directoryEta2760GeV->Get("EtatoPi0RatioConversion");   
   TGraphAsymmErrors* graphRatioEtaPi0SystErr2760GeV=          (TGraphAsymmErrors*)directoryEta2760GeV->Get("EtatoPi0RatioConversionSys");

	histoEtaMassData2760GeV->Scale(1000.);
	histoEtaMassMC2760GeV->Scale(1000.);
	histoEtaMassData2760GeV->SetBinContent(histoEtaMassData2760GeV->GetNbinsX(),0);
	histoEtaMassMC2760GeV->SetBinContent(histoEtaMassMC2760GeV->GetNbinsX(),0);
	histoEtaWidthData2760GeV->SetBinContent(histoEtaWidthData2760GeV->GetNbinsX(),10000.);
	histoEtaWidthMC2760GeV->SetBinContent(histoEtaWidthMC2760GeV->GetNbinsX(),10000.);

	
	TH1D*	histoPCMMassDataPP2760GeV = 		(TH1D*)directoryPi02760GeV->Get("MassPi0");
	TH1D*	histoPCMWidthDataPP2760GeV = 		(TH1D*)directoryPi02760GeV->Get("FWHMPi0MeV");
	TH1D*	histoPCMMassMCPP2760GeV = 			(TH1D*)directoryPi02760GeV->Get("TrueMassPi0");
	TH1D*	histoPCMWidthMCPP2760GeV = 			(TH1D*)directoryPi02760GeV->Get("TrueFWHMPi0MeV");
	histoPCMMassDataPP2760GeV->Scale(1000.);
	histoPCMMassMCPP2760GeV->Scale(1000.);
	histoPCMMassDataPP2760GeV->SetBinContent(histoPCMMassDataPP2760GeV->GetNbinsX(),0);
	histoPCMMassMCPP2760GeV->SetBinContent(histoPCMMassMCPP2760GeV->GetNbinsX(),0);
	histoPCMWidthDataPP2760GeV->SetBinContent(histoPCMWidthDataPP2760GeV->GetNbinsX(),10000.);
	histoPCMWidthMCPP2760GeV->SetBinContent(histoPCMWidthMCPP2760GeV->GetNbinsX(),10000.);

	TFile*	filePCMPbPb = 					new TFile(nameFilePbPb.Data());
   cout << "0-5%" << endl;
	TDirectory* directoryPi0PbPb0005 = 	(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_0-5%"); 
   TH1D* histoPCMPi0CorrectedSpec0005 =   (TH1D*)directoryPi0PbPb0005->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys0005 =        (TGraphAsymmErrors*)directoryPi0PbPb0005->Get("Pi0SystError"); 
	TH1D* histoPCMMassData0005 = 			(TH1D*)directoryPi0PbPb0005->Get("MassPi0");
	TH1D* histoPCMMassMC0005 = 			(TH1D*)directoryPi0PbPb0005->Get("TrueMassPi0");
	TH1D* histoPCMWidthData0005 = 		(TH1D*)directoryPi0PbPb0005->Get("FWHMPi0MeV");
	TH1D* histoPCMWidthMC0005 = 			(TH1D*)directoryPi0PbPb0005->Get("TrueFWHMPi0MeV");
   TH1D* histoTrueEffiPt0005 =        (TH1D*)directoryPi0PbPb0005->Get("Pi0_Efficiency");
   TH1D* histoAcceptPt0005 =          (TH1D*)directoryPi0PbPb0005->Get("Pi0_Acceptance");
   TH1D* histoMCYieldPt0005 =         (TH1D*)directoryPi0PbPb0005->Get("Pi0_HIJING_Input_Reweighted");
   TH1D* histoMCYieldPtWOWeights0005 = (TH1D*)directoryPi0PbPb0005->Get("Pi0_HIJING_Input");
   TH1D* histoPi0Weights0005 =        (TH1D*)directoryPi0PbPb0005->Get("Pi0_HIJING_Weights");
   TH1D* histoRawYield0005 =          (TH1D*)directoryPi0PbPb0005->Get("Pi0_RawYieldPerEvent");
   
	histoPCMMassData0005->Scale(1000.);
	histoPCMMassMC0005->Scale(1000.);
	histoPCMMassData0005->SetBinContent(1,0);
	histoPCMMassData0005->SetBinContent(2,0);
	histoPCMMassData0005->SetBinContent(histoPCMMassData0005->GetNbinsX(),0);
	histoPCMMassData0005->SetBinContent(histoPCMMassData0005->GetNbinsX()-1,0);
	histoPCMMassMC0005->SetBinContent(1,0);
	histoPCMMassMC0005->SetBinContent(2,0);
	histoPCMMassMC0005->SetBinContent(histoPCMMassMC0005->GetNbinsX(),0);
	histoPCMMassMC0005->SetBinContent(histoPCMMassMC0005->GetNbinsX()-1,0);
	histoPCMWidthData0005->SetBinContent(1,10000);
	histoPCMWidthData0005->SetBinContent(2,10000);
	histoPCMWidthData0005->SetBinContent(histoPCMWidthData0005->GetNbinsX(),10000.);
	histoPCMWidthData0005->SetBinContent(histoPCMWidthData0005->GetNbinsX()-1,10000.);
	histoPCMWidthMC0005->SetBinContent(1,10000);
	histoPCMWidthMC0005->SetBinContent(2,10000);
	histoPCMWidthMC0005->SetBinContent(histoPCMWidthMC0005->GetNbinsX(),10000.);
	histoPCMWidthMC0005->SetBinContent(histoPCMWidthMC0005->GetNbinsX()-1,10000.);

   cout << "0-10%" << endl;
   TDirectory* directoryPi0PbPb0010 =           (TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_0-10%"); 
   TH1D* histoPCMPi0CorrectedSpec0010 =   (TH1D*)directoryPi0PbPb0010->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys0010 =        (TGraphAsymmErrors*)directoryPi0PbPb0010->Get("Pi0SystError"); 
   TH1D* histoTrueEffiPt0010 =        (TH1D*)directoryPi0PbPb0010->Get("Pi0_Efficiency");
   TH1D* histoAcceptPt0010 =          (TH1D*)directoryPi0PbPb0010->Get("Pi0_Acceptance");
   TH1D* histoMCYieldPt0010 =         (TH1D*)directoryPi0PbPb0010->Get("Pi0_HIJING_Input_Reweighted");
   TH1D* histoMCYieldPtWOWeights0010 = (TH1D*)directoryPi0PbPb0010->Get("Pi0_HIJING_Input");
   TH1D* histoPi0Weights0010 =        (TH1D*)directoryPi0PbPb0010->Get("Pi0_HIJING_Weights");
   TH1D* histoRawYield0010 =          (TH1D*)directoryPi0PbPb0010->Get("Pi0_RawYieldPerEvent");
   
   cout << "0-20%" << endl;
   TDirectory* directoryPi0PbPb0020 =           (TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_0-20%"); 
   TH1D* histoPCMPi0CorrectedSpec0020 =   NULL;
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys0020 = NULL;
   TH1D* histoMCYieldPt0020 =         NULL;
   TH1D* histoMCYieldPtWOWeights0020 = NULL; 
   if (directoryPi0PbPb0020){
      histoPCMPi0CorrectedSpec0020 =   (TH1D*)directoryPi0PbPb0020->Get("CorrectedYieldPi0");   
      graphPCMPi0CorrectedSpecSys0020 =        (TGraphAsymmErrors*)directoryPi0PbPb0020->Get("Pi0SystError"); 
      histoMCYieldPt0020 =         (TH1D*)directoryPi0PbPb0020->Get("Pi0_HIJING_Input_Reweighted");
      histoMCYieldPtWOWeights0020 = (TH1D*)directoryPi0PbPb0020->Get("Pi0_HIJING_Input");
   //    TH1D* histoTrueEffiPt0020 =        (TH1D*)directoryPi0PbPb0020->Get("Pi0_Efficiency");
   //    TH1D* histoAcceptPt0020 =          (TH1D*)directoryPi0PbPb0020->Get("Pi0_Acceptance");      
   //    TH1D* histoPi0Weights0020 =        (TH1D*)directoryPi0PbPb0020->Get("Pi0_HIJING_Weights");
   //    TH1D* histoRawYield0020 =          (TH1D*)directoryPi0PbPb0020->Get("Pi0_RawYieldPerEvent");
   }
   
   cout << "0-40%" << endl;
   TDirectory* directoryPi0PbPb0040 =           (TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_0-40%"); 
   TH1D* histoPCMPi0CorrectedSpec0040 =   NULL;
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys0040 = NULL;
   TH1D* histoMCYieldPt0040 =         NULL;
   TH1D* histoMCYieldPtWOWeights0040 = NULL; 
   if (directoryPi0PbPb0040){
      histoPCMPi0CorrectedSpec0040 =   (TH1D*)directoryPi0PbPb0040->Get("CorrectedYieldPi0");   
      graphPCMPi0CorrectedSpecSys0040 =        (TGraphAsymmErrors*)directoryPi0PbPb0040->Get("Pi0SystError"); 
      histoMCYieldPt0040 =         (TH1D*)directoryPi0PbPb0040->Get("Pi0_HIJING_Input_Reweighted");
      histoMCYieldPtWOWeights0040 = (TH1D*)directoryPi0PbPb0040->Get("Pi0_HIJING_Input");
   //    TH1D* histoTrueEffiPt0040 =        (TH1D*)directoryPi0PbPb0040->Get("Pi0_Efficiency");
   //    TH1D* histoAcceptPt0040 =          (TH1D*)directoryPi0PbPb0040->Get("Pi0_Acceptance");      
   //    TH1D* histoPi0Weights0040 =        (TH1D*)directoryPi0PbPb0040->Get("Pi0_HIJING_Weights");
   //    TH1D* histoRawYield0040 =          (TH1D*)directoryPi0PbPb0040->Get("Pi0_RawYieldPerEvent");
   }
   
   cout << "0-80%" << endl;
   TDirectory* directoryPi0PbPb0080 =           (TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_0-80%"); 
   TH1D* histoPCMPi0CorrectedSpec0080 =   NULL;
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys0080 = NULL;
   TH1D* histoMCYieldPt0080 =         NULL;
   TH1D* histoMCYieldPtWOWeights0080 = NULL; 
   if (directoryPi0PbPb0080){
      histoPCMPi0CorrectedSpec0080 =   (TH1D*)directoryPi0PbPb0080->Get("CorrectedYieldPi0");   
      graphPCMPi0CorrectedSpecSys0080 =        (TGraphAsymmErrors*)directoryPi0PbPb0080->Get("Pi0SystError"); 
      histoMCYieldPt0080 =         (TH1D*)directoryPi0PbPb0080->Get("Pi0_HIJING_Input_Reweighted");
      histoMCYieldPtWOWeights0080 = (TH1D*)directoryPi0PbPb0080->Get("Pi0_HIJING_Input");
   //    TH1D* histoTrueEffiPt0080 =        (TH1D*)directoryPi0PbPb0080->Get("Pi0_Efficiency");
   //    TH1D* histoAcceptPt0080 =          (TH1D*)directoryPi0PbPb0080->Get("Pi0_Acceptance");      
   //    TH1D* histoPi0Weights0080 =        (TH1D*)directoryPi0PbPb0080->Get("Pi0_HIJING_Weights");
   //    TH1D* histoRawYield0080 =          (TH1D*)directoryPi0PbPb0080->Get("Pi0_RawYieldPerEvent");
   }
   
   cout << "40-80%" << endl;
   TDirectory* directoryPi0PbPb4080 =           (TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_40-80%"); 
   TH1D* histoPCMPi0CorrectedSpec4080 =   NULL;
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys4080 = NULL;
   TH1D* histoMCYieldPt4080 =         NULL;
   TH1D* histoMCYieldPtWOWeights4080 = NULL; 
   if (directoryPi0PbPb4080){
      histoPCMPi0CorrectedSpec4080 =   (TH1D*)directoryPi0PbPb4080->Get("CorrectedYieldPi0");   
      graphPCMPi0CorrectedSpecSys4080 =        (TGraphAsymmErrors*)directoryPi0PbPb4080->Get("Pi0SystError"); 
      histoMCYieldPt4080 =         (TH1D*)directoryPi0PbPb4080->Get("Pi0_HIJING_Input_Reweighted");
      histoMCYieldPtWOWeights4080 = (TH1D*)directoryPi0PbPb4080->Get("Pi0_HIJING_Input");
   //    TH1D* histoTrueEffiPt4080 =        (TH1D*)directoryPi0PbPb4080->Get("Pi0_Efficiency");
   //    TH1D* histoAcceptPt4080 =          (TH1D*)directoryPi0PbPb4080->Get("Pi0_Acceptance");      
   //    TH1D* histoPi0Weights4080 =        (TH1D*)directoryPi0PbPb4080->Get("Pi0_HIJING_Weights");
   //    TH1D* histoRawYield4080 =          (TH1D*)directoryPi0PbPb4080->Get("Pi0_RawYieldPerEvent");
   }
   
	cout << "5-10%" << endl;
	TDirectory* directoryPi0PbPb0510 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_5-10%"); 
   TH1D* histoPCMPi0CorrectedSpec0510 =   (TH1D*)directoryPi0PbPb0510->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys0510 =        (TGraphAsymmErrors*)directoryPi0PbPb0510->Get("Pi0SystError"); 
	TH1D* histoPCMMassData0510 = 			(TH1D*)directoryPi0PbPb0510->Get("MassPi0");
	TH1D* histoPCMMassMC0510 = 			(TH1D*)directoryPi0PbPb0510->Get("TrueMassPi0");
	TH1D* histoPCMWidthData0510 = 		(TH1D*)directoryPi0PbPb0510->Get("FWHMPi0MeV");
	TH1D* histoPCMWidthMC0510 = 			(TH1D*)directoryPi0PbPb0510->Get("TrueFWHMPi0MeV");
   TH1D* histoTrueEffiPt0510 =        (TH1D*)directoryPi0PbPb0510->Get("Pi0_Efficiency");
   TH1D* histoAcceptPt0510 =          (TH1D*)directoryPi0PbPb0510->Get("Pi0_Acceptance");
   TH1D* histoMCYieldPt0510 =         (TH1D*)directoryPi0PbPb0510->Get("Pi0_HIJING_Input_Reweighted");
   TH1D* histoMCYieldPtWOWeights0510 = (TH1D*)directoryPi0PbPb0510->Get("Pi0_HIJING_Input");
   TH1D* histoPi0Weights0510 =        (TH1D*)directoryPi0PbPb0510->Get("Pi0_HIJING_Weights");
   TH1D* histoRawYield0510 =          (TH1D*)directoryPi0PbPb0510->Get("Pi0_RawYieldPerEvent");
   
	histoPCMMassData0510->Scale(1000.);
	histoPCMMassMC0510->Scale(1000.);
	histoPCMMassData0510->SetBinContent(histoPCMMassData0510->GetNbinsX(),0);
	histoPCMMassData0510->SetBinContent(histoPCMMassData0510->GetNbinsX()-1,0);
	histoPCMMassMC0510->SetBinContent(histoPCMMassMC0510->GetNbinsX(),0);
	histoPCMMassMC0510->SetBinContent(histoPCMMassMC0510->GetNbinsX()-1,0);
	histoPCMWidthData0510->SetBinContent(histoPCMWidthData0510->GetNbinsX(),10000.);
	histoPCMWidthData0510->SetBinContent(histoPCMWidthData0510->GetNbinsX()-1,10000.);
	histoPCMWidthMC0510->SetBinContent(histoPCMWidthMC0510->GetNbinsX(),10000.);
	histoPCMWidthMC0510->SetBinContent(histoPCMWidthMC0510->GetNbinsX()-1,10000.);

	cout << "10-20%" << endl;
	TDirectory* directoryPi0PbPb1020 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_10-20%"); 
   TH1D* histoPCMPi0CorrectedSpec1020 =   (TH1D*)directoryPi0PbPb1020->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys1020 =        (TGraphAsymmErrors*)directoryPi0PbPb1020->Get("Pi0SystError"); 
	TH1D* histoPCMMassData1020 = 			(TH1D*)directoryPi0PbPb1020->Get("MassPi0");
	TH1D* histoPCMMassMC1020 = 			(TH1D*)directoryPi0PbPb1020->Get("TrueMassPi0");
	TH1D* histoPCMWidthData1020 = 		(TH1D*)directoryPi0PbPb1020->Get("FWHMPi0MeV");
	TH1D* histoPCMWidthMC1020 = 			(TH1D*)directoryPi0PbPb1020->Get("TrueFWHMPi0MeV");
   TH1D* histoTrueEffiPt1020 =        (TH1D*)directoryPi0PbPb1020->Get("Pi0_Efficiency");
   TH1D* histoAcceptPt1020 =          (TH1D*)directoryPi0PbPb1020->Get("Pi0_Acceptance");
   TH1D* histoMCYieldPt1020 =         (TH1D*)directoryPi0PbPb1020->Get("Pi0_HIJING_Input_Reweighted");
   TH1D* histoMCYieldPtWOWeights1020 = (TH1D*)directoryPi0PbPb1020->Get("Pi0_HIJING_Input");
   TH1D* histoPi0Weights1020 =        (TH1D*)directoryPi0PbPb1020->Get("Pi0_HIJING_Weights");
   TH1D* histoRawYield1020 =          (TH1D*)directoryPi0PbPb1020->Get("Pi0_RawYieldPerEvent");
	histoPCMMassData1020->Scale(1000.);
	histoPCMMassMC1020->Scale(1000.);
	histoPCMMassData1020->SetBinContent(histoPCMMassData1020->GetNbinsX(),0);
	histoPCMMassData1020->SetBinContent(histoPCMMassData1020->GetNbinsX()-1,0);
	histoPCMMassMC1020->SetBinContent(histoPCMMassMC1020->GetNbinsX(),0);
	histoPCMMassMC1020->SetBinContent(histoPCMMassMC1020->GetNbinsX()-1,0);
	histoPCMWidthData1020->SetBinContent(histoPCMWidthData1020->GetNbinsX(),10000.);
	histoPCMWidthData1020->SetBinContent(histoPCMWidthData1020->GetNbinsX()-1,10000.);
	histoPCMWidthMC1020->SetBinContent(histoPCMWidthMC1020->GetNbinsX(),10000.);
	histoPCMWidthMC1020->SetBinContent(histoPCMWidthMC1020->GetNbinsX()-1,10000.);
	
	cout << "20-40%" << endl;
	TDirectory* directoryPi0PbPb2040 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_20-40%"); 
	TH1D* histoPCMPi0CorrectedSpec2040 =   (TH1D*)directoryPi0PbPb2040->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys2040 =        (TGraphAsymmErrors*)directoryPi0PbPb2040->Get("Pi0SystError"); 
   TH1D* histoTrueEffiPt2040 =        (TH1D*)directoryPi0PbPb2040->Get("Pi0_Efficiency");
   TH1D* histoAcceptPt2040 =          (TH1D*)directoryPi0PbPb2040->Get("Pi0_Acceptance");
   TH1D* histoMCYieldPt2040 =         (TH1D*)directoryPi0PbPb2040->Get("Pi0_HIJING_Input_Reweighted");
   TH1D* histoMCYieldPtWOWeights2040 = (TH1D*)directoryPi0PbPb2040->Get("Pi0_HIJING_Input");
   TH1D* histoPi0Weights2040 =        (TH1D*)directoryPi0PbPb2040->Get("Pi0_HIJING_Weights");
   TH1D* histoRawYield2040 =          (TH1D*)directoryPi0PbPb2040->Get("Pi0_RawYieldPerEvent");
   
   
	cout << "40-60%" << endl;
	TDirectory* directoryPi0PbPb4060 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_40-60%"); 
   TH1D* histoPCMPi0CorrectedSpec4060 =   (TH1D*)directoryPi0PbPb4060->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys4060 =        (TGraphAsymmErrors*)directoryPi0PbPb4060->Get("Pi0SystError"); 
   TH1D* histoTrueEffiPt4060 =        (TH1D*)directoryPi0PbPb4060->Get("Pi0_Efficiency");
   TH1D* histoAcceptPt4060 =          (TH1D*)directoryPi0PbPb4060->Get("Pi0_Acceptance");
   TH1D* histoMCYieldPt4060 =         (TH1D*)directoryPi0PbPb4060->Get("Pi0_HIJING_Input_Reweighted");
   TH1D* histoMCYieldPtWOWeights4060 = (TH1D*)directoryPi0PbPb4060->Get("Pi0_HIJING_Input");
   TH1D* histoPi0Weights4060 =        (TH1D*)directoryPi0PbPb4060->Get("Pi0_HIJING_Weights");
   TH1D* histoRawYield4060 =          (TH1D*)directoryPi0PbPb4060->Get("Pi0_RawYieldPerEvent");
   
	
	cout << "60-80%" << endl;
	TDirectory* directoryPi0PbPb6080 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_60-80%");
   TH1D* histoPCMPi0CorrectedSpec6080 =   (TH1D*)directoryPi0PbPb6080->Get("CorrectedYieldPi0");   
   TGraphAsymmErrors* graphPCMPi0CorrectedSpecSys6080 =        (TGraphAsymmErrors*)directoryPi0PbPb6080->Get("Pi0SystError"); 
	TH1D* histoPCMMassData6080 = 			(TH1D*)directoryPi0PbPb6080->Get("MassPi0");
	TH1D* histoPCMMassMC6080 = 			(TH1D*)directoryPi0PbPb6080->Get("TrueMassPi0");
	TH1D* histoPCMWidthData6080 = 		(TH1D*)directoryPi0PbPb6080->Get("FWHMPi0MeV");
	TH1D* histoPCMWidthMC6080 = 			(TH1D*)directoryPi0PbPb6080->Get("TrueFWHMPi0MeV");
   TH1D* histoTrueEffiPt6080 =        (TH1D*)directoryPi0PbPb6080->Get("Pi0_Efficiency");
   TH1D* histoAcceptPt6080 =          (TH1D*)directoryPi0PbPb6080->Get("Pi0_Acceptance");
   TH1D* histoMCYieldPt6080 =         (TH1D*)directoryPi0PbPb6080->Get("Pi0_HIJING_Input_Reweighted");
   TH1D* histoMCYieldPtWOWeights6080 = (TH1D*)directoryPi0PbPb6080->Get("Pi0_HIJING_Input");
   TH1D* histoPi0Weights6080 =        (TH1D*)directoryPi0PbPb6080->Get("Pi0_HIJING_Weights");
   TH1D* histoRawYield6080 =          (TH1D*)directoryPi0PbPb6080->Get("Pi0_RawYieldPerEvent");
	histoPCMMassData6080->Scale(1000.);
	histoPCMMassMC6080->Scale(1000.);
	histoPCMMassData6080->SetBinContent(histoPCMMassData6080->GetNbinsX(),0);
	histoPCMMassData6080->SetBinContent(histoPCMMassData6080->GetNbinsX()-1,0);
	histoPCMMassMC6080->SetBinContent(histoPCMMassMC6080->GetNbinsX(),0);
	histoPCMMassMC6080->SetBinContent(histoPCMMassMC6080->GetNbinsX()-1,0);
	histoPCMWidthData6080->SetBinContent(histoPCMWidthData6080->GetNbinsX(),10000.);
	histoPCMWidthData6080->SetBinContent(histoPCMWidthData6080->GetNbinsX()-1,10000.);
	histoPCMWidthMC6080->SetBinContent(histoPCMWidthMC6080->GetNbinsX(),10000.);
	histoPCMWidthMC6080->SetBinContent(histoPCMWidthMC6080->GetNbinsX()-1,10000.);

	TDirectory* directoryEtaPbPb0010 = 				(TDirectory*)filePCMPbPb->Get("Eta_PbPb_2.76TeV_0-10%"); 
	TH1D* histoPCMEtaCorrectedSpec0010 = 			(TH1D*)directoryEtaPbPb0010->Get("CorrectedYieldEta");	
	TGraphAsymmErrors* graphPCMEtaCorrectedSpecSys0010 = 			(TGraphAsymmErrors*)directoryEtaPbPb0010->Get("EtaSystError");	
	TH1D* histoEtaTrueEffiPt0010 =        (TH1D*)directoryEtaPbPb0010->Get("Eta_Efficiency");
   TH1D* histoEtaAcceptPt0010 =          (TH1D*)directoryEtaPbPb0010->Get("Eta_Acceptance");
   TH1D* histoEtaMCYieldPt0010 =         (TH1D*)directoryEtaPbPb0010->Get("Eta_HIJING_Input_Reweighted");
   TH1D* histoEtaRawYield0010 =          (TH1D*)directoryEtaPbPb0010->Get("Eta_RawYieldPerEvent");
   
	TDirectory* directoryEtaPbPb1020 = 				(TDirectory*)filePCMPbPb->Get("Eta_PbPb_2.76TeV_10-20%"); 
	TH1D* histoPCMEtaCorrectedSpec1020 = 			(TH1D*)directoryEtaPbPb1020->Get("CorrectedYieldEta");	
	TGraphAsymmErrors* graphPCMEtaCorrectedSpecSys1020 = 			(TGraphAsymmErrors*)directoryEtaPbPb1020->Get("EtaSystError");	
	TH1D* histoEtaTrueEffiPt1020 =        (TH1D*)directoryEtaPbPb1020->Get("Eta_Efficiency");
   TH1D* histoEtaAcceptPt1020 =          (TH1D*)directoryEtaPbPb1020->Get("Eta_Acceptance");
   TH1D* histoEtaMCYieldPt1020 =         (TH1D*)directoryEtaPbPb1020->Get("Eta_HIJING_Input_Reweighted");
   TH1D* histoEtaRawYield1020 =          (TH1D*)directoryEtaPbPb1020->Get("Eta_RawYieldPerEvent");
   
	TDirectory* directoryEtaPbPb2040 = 				(TDirectory*)filePCMPbPb->Get("Eta_PbPb_2.76TeV_20-40%"); 
	TH1D* histoPCMEtaCorrectedSpec2040 = 			(TH1D*)directoryEtaPbPb2040->Get("CorrectedYieldEta");	
	TGraphAsymmErrors* graphPCMEtaCorrectedSpecSys2040 = 			(TGraphAsymmErrors*)directoryEtaPbPb2040->Get("EtaSystError");	
	TH1D* histoEtaTrueEffiPt2040 =        (TH1D*)directoryEtaPbPb2040->Get("Eta_Efficiency");
   TH1D* histoEtaAcceptPt2040 =          (TH1D*)directoryEtaPbPb2040->Get("Eta_Acceptance");
   TH1D* histoEtaMCYieldPt2040 =         (TH1D*)directoryEtaPbPb2040->Get("Eta_HIJING_Input_Reweighted");
   TH1D* histoEtaRawYield2040 =          (TH1D*)directoryEtaPbPb2040->Get("Eta_RawYieldPerEvent");
   
	TDirectory* directoryEtaPbPb4060 = 				(TDirectory*)filePCMPbPb->Get("Eta_PbPb_2.76TeV_40-60%"); 
	TH1D* histoPCMEtaCorrectedSpec4060 = 			(TH1D*)directoryEtaPbPb4060->Get("CorrectedYieldEta");	
	TGraphAsymmErrors* graphPCMEtaCorrectedSpecSys4060 = 			(TGraphAsymmErrors*)directoryEtaPbPb4060->Get("EtaSystError");	
	TH1D* histoEtaTrueEffiPt4060 =        (TH1D*)directoryEtaPbPb4060->Get("Eta_Efficiency");
   TH1D* histoEtaAcceptPt4060 =          (TH1D*)directoryEtaPbPb4060->Get("Eta_Acceptance");
   TH1D* histoEtaMCYieldPt4060 =         (TH1D*)directoryEtaPbPb4060->Get("Eta_HIJING_Input_Reweighted");
   TH1D* histoEtaRawYield4060 =          (TH1D*)directoryEtaPbPb4060->Get("Eta_RawYieldPerEvent");
   
	TDirectory* directoryEtaPbPb6080 = 				(TDirectory*)filePCMPbPb->Get("Eta_PbPb_2.76TeV_60-80%"); 
	TH1D* histoPCMEtaCorrectedSpec6080 = 			(TH1D*)directoryEtaPbPb6080->Get("CorrectedYieldEta");	
	TGraphAsymmErrors* graphPCMEtaCorrectedSpecSys6080 = 			(TGraphAsymmErrors*)directoryEtaPbPb6080->Get("EtaSystError");	
   TH1D* histoEtaTrueEffiPt6080 =        (TH1D*)directoryEtaPbPb6080->Get("Eta_Efficiency");
   TH1D* histoEtaAcceptPt6080 =          (TH1D*)directoryEtaPbPb6080->Get("Eta_Acceptance");
   TH1D* histoEtaMCYieldPt6080 =         (TH1D*)directoryEtaPbPb6080->Get("Eta_HIJING_Input_Reweighted");
   TH1D* histoEtaRawYield6080 =          (TH1D*)directoryEtaPbPb6080->Get("Eta_RawYieldPerEvent");
   TH1D* histoEtaMassData6080 =          (TH1D*)directoryEtaPbPb6080->Get("MassEta");
   TH1D* histoEtaMassMC6080 =            (TH1D*)directoryEtaPbPb6080->Get("TrueMassEta");
   TH1D* histoEtaWidthData6080 =         (TH1D*)directoryEtaPbPb6080->Get("FWHMEtaMeV");
   TH1D* histoEtaWidthMC6080 =           (TH1D*)directoryEtaPbPb6080->Get("TrueFWHMEtaMeV");
   histoEtaMassData6080->Scale(1000.);
   histoEtaMassMC6080->Scale(1000.);
   
   TFile* filePCMpPb =               new TFile("data_GammaConversionResultsFullCorrection_pPb_MidRap_CentBinning_Lego.root");
   TDirectory* directoryPCMPi0pPb =             (TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
   TH1D* histoPCMYieldPi0pPb =            (TH1D*)directoryPCMPi0pPb->Get("CorrectedYieldPi0");
//    histoPCMRAWYieldPi0pPb =            (TH1D*)directoryPCMPi0pPb->Get("RawYieldPi0");
   TGraphAsymmErrors*graphPCMYieldPi0SysErrpPb=    (TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystError"); 
   TH1D* histoPCMMassPi0DatapPb =         (TH1D*)directoryPCMPi0pPb->Get("MassPi0");
   TH1D* histoPCMMassPi0MCpPb =        (TH1D*)directoryPCMPi0pPb->Get("TrueMassPi0");
   TH1D* histoPCMWidthPi0DatapPb =     (TH1D*)directoryPCMPi0pPb->Get("FWHMPi0MeV");
   TH1D* histoPCMWidthPi0MCpPb =          (TH1D*)directoryPCMPi0pPb->Get("TrueFWHMPi0MeV");
   histoPCMMassPi0DatapPb->Scale(1000.);
   histoPCMMassPi0MCpPb->Scale(1000.);
   TDirectory* directoryPCMEtapPb =             (TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_0-100%"); 
   TH1D* histoPCMYieldEtapPb =            (TH1D*)directoryPCMEtapPb->Get("CorrectedYieldEta");
//    histoPCMRAWYieldEtapPb =            (TH1D*)directoryPCMEtapPb->Get("RawYieldEta");
   TGraphAsymmErrors* graphPCMYieldEtaSysErrpPb=    (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtaSystError"); 
   TH1D* histoPCMEtaPi0RatiopPb =            (TH1D*)directoryPCMEtapPb->Get("EtatoPi0Ratio");
   TGraphAsymmErrors* graphPCMEtaPi0RatioSysErrpPb=    (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtatoPi0RatioSys"); 
   TH1D* histoPCMMassEtaDatapPb =         (TH1D*)directoryPCMEtapPb->Get("MassEta");
   TH1D* histoPCMMassEtaMCpPb =        (TH1D*)directoryPCMEtapPb->Get("TrueMassEta");
   TH1D* histoPCMWidthEtaDatapPb =     (TH1D*)directoryPCMEtapPb->Get("FWHMEtaMeV");
   TH1D* histoPCMWidthEtaMCpPb =          (TH1D*)directoryPCMEtapPb->Get("TrueFWHMEtaMeV");
   histoPCMMassEtaDatapPb->Scale(1000.);
   histoPCMMassEtaMCpPb->Scale(1000.);
   
	TFile*	fileCocktail = 					new TFile("CocktailInput/cocktail_allCentpluspp.root");
	TDirectory* directoryCocktailpp2760GeV = 				(TDirectory*)fileCocktail->Get("cocktail_pp_2760GeV_qcd"); 
	TH1D* histoEtaFromCocktailpp2760GeV = (TH1D*)directoryCocktailpp2760GeV->Get("ptEta");	
	TDirectory* directoryCocktail0010 = 				(TDirectory*)fileCocktail->Get("cocktail_PbPb_0010_qcd"); 
	TH1D* histoEtaFromCocktail0010 = (TH1D*)directoryCocktail0010->Get("ptEta");	
	TDirectory* directoryCocktail0510 = 				(TDirectory*)fileCocktail->Get("cocktail_PbPb_0510_qcd"); 
	TH1D* histoEtaFromCocktail0510 = (TH1D*)directoryCocktail0510->Get("ptEta");	
	TDirectory* directoryCocktail1020 = 				(TDirectory*)fileCocktail->Get("cocktail_PbPb_1020_qcd"); 
	TH1D* histoEtaFromCocktail1020 = (TH1D*)directoryCocktail1020->Get("ptEta");	
	TDirectory* directoryCocktail2040 = 				(TDirectory*)fileCocktail->Get("cocktail_PbPb_2040_qcd"); 
	TH1D* histoEtaFromCocktail2040 = (TH1D*)directoryCocktail2040->Get("ptEta");	
	TDirectory* directoryCocktail4060 = 				(TDirectory*)fileCocktail->Get("cocktail_PbPb_4060_qcd"); 
	TH1D* histoEtaFromCocktail4060 = (TH1D*)directoryCocktail4060->Get("ptEta");	
	TDirectory* directoryCocktail6080 = 				(TDirectory*)fileCocktail->Get("cocktail_PbPb_6080_qcd"); 
	TH1D* histoEtaFromCocktail6080 = (TH1D*)directoryCocktail6080->Get("ptEta");	
	cout << "here 7TeV" << endl;
   
	TFile* fileNeutralPionCombDataPP = new TFile("CombinedResultsPaperX.root");
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
	TFile* fileNeutralPion2760GeVDataPPEffiAddedSig = new TFile("0000011002093663003800000_01631031009/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistosAddSig_0000011002093663003800000_01631031009.root");
	TH1D* histoEffiAddedSigPP = (TH1D*)fileNeutralPion2760GeVDataPPEffiAddedSig->Get("TrueMesonEffiPt");
	
   TFile* fileNeutralPion7TeVDataPPEffi2760GeVCut = new TFile("/home/fredi/Photon/Results/ppTests/0000011002093663003800000_01631031009/7TeV/Pi0_MC_GammaConvV1CorrectionHistosD_0000011002093663003800000_01631031009.root");
   TH1D* histoEffi7TeV2760GeVCut = (TH1D*)fileNeutralPion7TeVDataPPEffi2760GeVCut->Get("TrueMesonEffiPt");
   
   TFile* fileNeutralPion2760GeVDataPPEffiWithSDD = new TFile("/home/fredi/Photon/Results/ppTests/0002011002093663003800000_01631031009/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_0002011002093663003800000_01631031009.root");
   TH1D* histoEffi2760GeVWithSDD = (TH1D*)fileNeutralPion2760GeVDataPPEffiWithSDD->Get("TrueMesonEffiPt");
   
   cout << "efficiencies MinBias 2.76TeV GeV" << endl;
	TFile* fileNeutralPion2760GeVDataPPEffiMinBias = new TFile("0000011002093663003800000_01631031009/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistosMinBias_0000011002093663003800000_01631031009.root");
	TH1D* histoEffiMinBiasPP = (TH1D*)fileNeutralPion2760GeVDataPPEffiMinBias->Get("TrueMesonEffiPt");

	cout << "old spectra 2.76TeV GeV" << endl;
	TFile* fileNeutralPionDataPPOld = new TFile("0000011002093663003800000_01631031009_old/data_GammaConversionResultsFullCorrectionNoBinShifting.root");
	TDirectory*	directoryPi02760GeVOld = 			(TDirectory*)fileNeutralPionDataPPOld->Get("Pi02.76TeV"); 
	TH1D* histoRawYieldPi02760GeVOld = 				(TH1D*)directoryPi02760GeVOld->Get("RAWYieldPerEventsPi0");
	TH1D* histoCorrectedYieldPi02760GeVOld = 				(TH1D*)directoryPi02760GeVOld->Get("CorrectedYieldPi0");
	
   cout << "MC spectra 2.76TeV GeV" << endl;
	TFile* fileNeutralPion2760GeVDataPP = new TFile("0000011002093663003800000_01631031009/2.76TeV/Pi0_data_GammaConvV1Correction_0000011002093663003800000_01631031009.root");
	TH1D*	histoMCYieldPi02760GeVFinal = (TH1D*)fileNeutralPion2760GeVDataPP->Get("MCYield_Meson_oldBin");
   TFile* fileEta2760GeVDataPP = new TFile("0000011002093663003800000_01631031009/2.76TeV/Eta_data_GammaConvV1Correction_0000011002093663003800000_01631031009.root");
   TH1D* histoMCYieldEta2760GeVFinal = (TH1D*)fileEta2760GeVDataPP->Get("MCYield_Meson_oldBin");
   
   cout << "MC spectra added signal 2.76TeV GeV" << endl;
   TFile* fileNeutralPion2760GeVDataPPAddSig = new TFile("0000011002093663003800000_01631031009/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistosAddSig_0000011002093663003800000_01631031009.root");
   TH1D* histoMCYieldPi02760GeVAddSig = (TH1D*)fileNeutralPion2760GeVDataPPAddSig->Get("MC_Meson_genPt_oldBin");
   TH1F *histoEventQualityMC =         (TH1F*)fileNeutralPion2760GeVDataPPAddSig->Get("NEvents");
   Float_t nEvtMC = GetNEvents(histoEventQualityMC);
   TString rapidityRange = "";
   Double_t deltaRapid =  ReturnRapidityStringAndDouble("01631031009", rapidityRange);
   Double_t scaling = 1./(2.*TMath::Pi());
   ScaleMCYield(histoMCYieldPi02760GeVAddSig,  deltaRapid,  scaling,  nEvtMC,  "Pi0" ,"kFALSE");
   cout << "MC spectra withough weighting added signal 2.76TeV GeV" << endl;
   TFile* fileNeutralPion2760GeVDataPPAddSigWOWeighting = new TFile("/home/fredi/Photon/Results/ppTestsWithoutWeighting/0000012002093663003800000_01631031009/2.76TeV/Pi0_data_GammaConvV1Correction_0000012002093663003800000_01631031009.root");
   TH1D* histoMCYieldPi02760GeVAddSigWOWeighting = (TH1D*)fileNeutralPion2760GeVDataPPAddSigWOWeighting->Get("MCYield_Meson_oldBin");
   
   cout << "MC spectra without weighting 2.76TeV GeV" << endl;
   TFile* fileNeutralPion2760GeVDataPPWOWeigthing = new TFile("/home/fredi/Photon/Results/ppTestsWithoutWeighting/0000011002093663003800000_01631031009/2.76TeV/Pi0_data_GammaConvV1Correction_0000011002093663003800000_01631031009.root");
   TH1D* histoMCYieldPi02760GeVWOWeighting = (TH1D*)fileNeutralPion2760GeVDataPPWOWeigthing->Get("MCYield_Meson_oldBin");
   cout << "MC spectra 7TeV GeV" << endl;
	TFile* fileNeutralPion7TeVDataPP = new TFile("900366208010033211360000000900/7TeV/Pi0_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
	TH1D*	histoMCYieldPi07TeV = (TH1D*)fileNeutralPion7TeVDataPP->Get("MCYield_Meson_oldBin");
   TFile* fileEta7TeVDataPP = new TFile("900366208010033211360000000900/7TeV/Eta_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
   TH1D* histoMCYieldEta7TeV = (TH1D*)fileEta7TeVDataPP->Get("MCYield_Meson_oldBin");
   cout << "MC spectra 0.9TeV GeV" << endl;
	TFile* fileNeutralPion900GeVDataPP = new TFile("900366208010033211360000000900/900GeV/Pi0_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
	TH1D*	histoMCYieldPi0900GeV = (TH1D*)fileNeutralPion900GeVDataPP->Get("MCYield_Meson_oldBin");
   TFile* fileEta900GeVDataPP = new TFile("900366208010033211360000000900/900GeV/Eta_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
   TH1D* histoMCYieldEta900GeV = (TH1D*)fileEta900GeVDataPP->Get("MCYield_Meson_oldBin");


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

   
   
	TCanvas* canvasEffAllEta = new TCanvas("canvasEffAllEta","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasEffAllEta, 0.1, 0.02, 0.035, 0.09);
		TH2F * histo2DEffEta;
		histo2DEffEta = new TH2F("histo2DEffEta","histo2DEffEta",1000,0,8,2000,0.e-3,3.7e-3	);
// 		histo2DEffEta->GetXaxis()->SetRangeUser(0.,12.);
		SetStyleHistoTH2ForGraphs(histo2DEffEta, "p_{T} (GeV/c)","#epsilon_{reco, #eta}",0.03,0.04, 0.03,0.04, 1.,1.);
		histo2DEffEta->Draw("copy");

		DrawGammaSetMarker(histoTrueEffPtEta7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);	
		histoTrueEffPtEta7TeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoTrueEffPtEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);	
		histoTrueEffPtEta2760GeV->DrawCopy("pe1,same"); 	
		DrawGammaSetMarker(histoTrueEffPtEta900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV);	
		histoTrueEffPtEta900GeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoEtaTrueEffiPt6080, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);	
		histoEtaTrueEffiPt6080->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoEtaTrueEffiPt0010, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);	
		histoEtaTrueEffiPt0010->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoEtaTrueEffiPt1020, markerStylePbPb1020, markerSizePbPb1020, colorCombPbPb1020, colorCombPbPb1020);	
		histoEtaTrueEffiPt1020->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoEtaTrueEffiPt2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);	
		histoEtaTrueEffiPt2040->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoEtaTrueEffiPt4060, markerStylePbPb4060, markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060);	
		histoEtaTrueEffiPt4060->DrawCopy("e1,same"); 	

		
		TLegend* legendEffiEta = new TLegend(0.12,0.78,0.83,0.93);
		legendEffiEta->SetFillColor(0);
		legendEffiEta->SetLineColor(0);
		legendEffiEta->SetTextSize(0.027);
		legendEffiEta->SetNColumns(2);
		legendEffiEta->AddEntry(histoEtaTrueEffiPt0010,"0-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendEffiEta->AddEntry(histoTrueEffPtEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
		legendEffiEta->AddEntry(histoEtaTrueEffiPt1020,"10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendEffiEta->AddEntry(histoTrueEffPtEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
		legendEffiEta->AddEntry(histoEtaTrueEffiPt2040,"20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendEffiEta->AddEntry(histoTrueEffPtEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
		legendEffiEta->AddEntry(histoEtaTrueEffiPt4060,"40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendEffiEta->AddEntry((TObject*)0, "","");
		legendEffiEta->AddEntry(histoEtaTrueEffiPt6080,"60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");

		legendEffiEta->Draw();
	
	canvasEffAllEta->SaveAs(Form("%s/EffCompEta.%s",outputDir.Data(),suffix.Data()));

		TCanvas* canvasRawYieldEta = new TCanvas("canvasRawYieldEta","",200,10,1350*1.4,1350);  // gives the page size
		DrawGammaCanvasSettings( canvasRawYieldEta, 0.12, 0.02, 0.035, 0.09);
		TH2F * histo2DRawEta;
// 		canvasRawYieldEta->SetLogx();
		canvasRawYieldEta->SetLogy();
		histo2DRawEta = new TH2F("histo2DRawEta","histo2DRawEta",1000,0.,9,2000,1.e-7,1e-4	);
		SetStyleHistoTH2ForGraphs(histo2DRawEta, "p_{T} (GeV/c)","#frac{dN_{raw}^{#eta}}{N_{evt} dp_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
		histo2DRawEta->Draw("copy");
	
		DrawGammaSetMarker(histoRawYieldEta7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);	
		histoRawYieldEta7TeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoRawYieldEta2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);	
		histoRawYieldEta2760GeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoRawYieldEta900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);	
		histoRawYieldEta900GeV->DrawCopy("e1,same"); 	

		TLatex *labelRawEta = new TLatex(0.34,0.88,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
		SetStyleTLatex( labelRawEta, 0.038,4);
		labelRawEta->Draw();

		
		TLegend* legendRawYieldEta = new TLegend(0.65,0.73,0.93,0.93);
		legendRawYieldEta->SetFillColor(0);
		legendRawYieldEta->SetLineColor(0);
		legendRawYieldEta->SetTextSize(0.04);
		legendRawYieldEta->AddEntry(histoRawYieldEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
		legendRawYieldEta->AddEntry(histoRawYieldEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
		legendRawYieldEta->AddEntry(histoRawYieldEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 		
		legendRawYieldEta->Draw();
	
	canvasRawYieldEta->SaveAs(Form("%s/RawYieldCompEta.%s",outputDir.Data(),suffix.Data()));

		TCanvas* canvasRawYieldPi0 = new TCanvas("canvasRawYieldPi0","",200,10,1350*1.4,1350);  // gives the page size
		DrawGammaCanvasSettings( canvasRawYieldPi0, 0.12, 0.02, 0.035, 0.09);
		TH2F * histo2DRawPi0;
// 		canvasRawYieldPi0->SetLogx();
		canvasRawYieldPi0->SetLogy();
		histo2DRawPi0 = new TH2F("histo2DRawPi0","histo2DRawPi0",1000,0.,20,2000,1.e-8,1e-3	);
		SetStyleHistoTH2ForGraphs(histo2DRawPi0, "p_{T} (GeV/c)","#frac{dN_{raw}^{#pi^{0}}}{N_{evt} dp_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
		histo2DRawPi0->Draw("copy");
	
		DrawGammaSetMarker(histoRawYieldPi07TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);	
		histoRawYieldPi07TeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoRawYieldPi02760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);	
		histoRawYieldPi02760GeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoRawYieldPi0900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);	
		histoRawYieldPi0900GeV->DrawCopy("e1,same"); 	

		TLatex *labelRawPi0 = new TLatex(0.34,0.88,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
		SetStyleTLatex( labelRawPi0, 0.038,4);
		labelRawPi0->Draw();

		
		TLegend* legendRawYieldPi0 = new TLegend(0.65,0.73,0.93,0.93);
		legendRawYieldPi0->SetFillColor(0);
		legendRawYieldPi0->SetLineColor(0);
		legendRawYieldPi0->SetTextSize(0.04);
		legendRawYieldPi0->AddEntry(histoRawYieldPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
		legendRawYieldPi0->AddEntry(histoRawYieldPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
		legendRawYieldPi0->AddEntry(histoRawYieldPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p"); 		
		legendRawYieldPi0->Draw();
	
	canvasRawYieldPi0->SaveAs(Form("%s/RawYieldCompPi0.%s",outputDir.Data(),suffix.Data()));

	TCanvas* canvasRawYieldPi02769GeV = new TCanvas("canvasRawYieldPi02769GeV","",200,10,1350*1.4,1350);  // gives the page size
		DrawGammaCanvasSettings( canvasRawYieldPi02769GeV, 0.12, 0.02, 0.035, 0.09);
// 		canvasRawYieldPi02769GeV->SetLogx();
		canvasRawYieldPi02769GeV->SetLogy();
		histo2DRawPi0->GetXaxis()->SetRangeUser(0.,10.1);
		histo2DRawPi0->Draw("copy");
	
		DrawGammaSetMarker(histoRawYieldPi02760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);	
		histoRawYieldPi02760GeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoRawYieldPi02760GeVOld, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);	
		histoRawYieldPi02760GeVOld->DrawCopy("e1,same"); 	

		labelRawPi0->Draw();

		
		TLegend* legendRawYieldPi02760GeV = new TLegend(0.65,0.73,0.93,0.93);
		legendRawYieldPi02760GeV->SetFillColor(0);
		legendRawYieldPi02760GeV->SetLineColor(0);
		legendRawYieldPi02760GeV->SetTextSize(0.04);
		legendRawYieldPi02760GeV->AddEntry(histoRawYieldPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
		legendRawYieldPi02760GeV->AddEntry(histoRawYieldPi0900GeV,"pp #sqrt{#it{s}} = 2.76 TeV, old","p"); 		
		legendRawYieldPi02760GeV->Draw();
	
	canvasRawYieldPi02769GeV->SaveAs(Form("%s/RawYieldCompPi02760GeV.%s",outputDir.Data(),suffix.Data()));
	
	TCanvas* canvasFraction2 = new TCanvas("canvasFraction2","",1550,1200);  // gives the page size
	canvasFraction2->SetTickx();
	canvasFraction2->SetTicky();
	canvasFraction2->SetGridx(0);
	canvasFraction2->SetGridy(0);
	canvasFraction2->SetLogy(0);
	canvasFraction2->SetLeftMargin(0.13);
	canvasFraction2->SetRightMargin(0.02);
	canvasFraction2->SetTopMargin(0.02);
	canvasFraction2->SetFillColor(0);

	TH1D* histoRatioRawSpectra2760NewOld = (TH1D*) histoRawYieldPi02760GeV->Clone("histoRatioRawSpectra2760NewOld");
	histoRatioRawSpectra2760NewOld->Divide(histoRawYieldPi02760GeVOld, histoRawYieldPi02760GeV, 1.,1.,"B");

	TH1D* histoRatioCorrSpectra2760NewOld = (TH1D*) histoCorrectedYieldPi02760GeV->Clone("histoRatioCorrSpectra2760NewOld");
	histoRatioCorrSpectra2760NewOld->Divide(histoCorrectedYieldPi02760GeVOld, histoCorrectedYieldPi02760GeV, 1.,1.,"B");

	DrawGammaSetMarker(histoRatioRawSpectra2760NewOld, 20, 1., kBlue+2, kBlue+2);
	DrawAutoGammaMesonHistos( histoRatioRawSpectra2760NewOld,
					"", "p_{T} (GeV/c)", "#frac{N^{raw}_{#pi^{0},old}}{N^{raw}_{#pi^{0},new}}, pp @ 2.76TeV",
					kFALSE, 5., 10e-10, kTRUE,
					kTRUE, 0.5, 1.5,
					kTRUE, 0., 7.9);
	canvasFraction2->Update();
	DrawGammaLines(0., 8.,1., 1.,0.1);

	canvasFraction2->SaveAs(Form("%s/RatioRawYieldNewOld.%s",outputDir.Data(),suffix.Data()));

	canvasFraction2->cd();
	DrawGammaSetMarker(histoRatioCorrSpectra2760NewOld, 20, 1., kBlue+2, kBlue+2);
	DrawAutoGammaMesonHistos( histoRatioCorrSpectra2760NewOld,
					"", "p_{T} (GeV/c)", "#frac{N^{corr}_{#pi^{0},old}}{N^{corr}_{#pi^{0},new}}, pp @ 2.76TeV",
					kFALSE, 5., 10e-10, kTRUE,
					kTRUE, 0.5, 1.5,
					kTRUE, 0., 7.9);
	canvasFraction2->Update();
	DrawGammaLines(0., 8.,1., 1.,0.1);

	canvasFraction2->SaveAs(Form("%s/RatioCorrYieldNewOld.%s",outputDir.Data(),suffix.Data()));
	
	TCanvas* canvasMCYieldPi0 = new TCanvas("canvasMCYieldPi0","",200,10,1350*1.2,1350);  // gives the page size
		DrawGammaCanvasSettings( canvasMCYieldPi0, 0.12, 0.02, 0.035, 0.09);
		TH2F * histo2DMCYieldPi0;
// 		canvasMCYieldPi0->SetLogx();
		canvasMCYieldPi0->SetLogy();
		histo2DMCYieldPi0 = new TH2F("histo2DMCYieldPi0","histo2DMCYieldPi0",1000,0.,20,2000,1.e-9,1e1	);
		SetStyleHistoTH2ForGraphs(histo2DMCYieldPi0, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}",0.035,0.04, 0.035,0.04, 1.,1.3);
		histo2DMCYieldPi0->Draw("copy");
	
      DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCM7TeVSysErr, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV,widthLinesBoxes, kTRUE);  
      graphInvYieldPi0PCM7TeVSysErr->Draw("e2,same");  
      graphInvYieldPi0PCM2760GeVSysErr= ScaleGraph(graphInvYieldPi0PCM2760GeVSysErr,1e-1);
      DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCM2760GeVSysErr, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV,widthLinesBoxes, kTRUE);
      graphInvYieldPi0PCM2760GeVSysErr->Draw("e2,same");  
      graphInvYieldPi0PCM900GeVSysErr= ScaleGraph(graphInvYieldPi0PCM900GeVSysErr,1e-2);
      DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCM900GeVSysErr, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV,widthLinesBoxes, kTRUE);
      graphInvYieldPi0PCM900GeVSysErr->Draw("e2,same");   

		DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCM7TeVStatErr, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);	
		graphInvYieldPi0PCM7TeVStatErr->Draw("pe1,same"); 	
		graphInvYieldPi0PCM2760GeVStatErr= ScaleGraph(graphInvYieldPi0PCM2760GeVStatErr,1e-1);
		DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCM2760GeVStatErr, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);	
		graphInvYieldPi0PCM2760GeVStatErr->Draw("pe1,same"); 	
		graphInvYieldPi0PCM900GeVStatErr= ScaleGraph(graphInvYieldPi0PCM900GeVStatErr,1e-2);
		DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCM900GeVStatErr, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV);	
		graphInvYieldPi0PCM900GeVStatErr->Draw("pe1,same"); 	

		SetStyleHisto(histoMCYieldPi07TeV, 1., 1, colorPi07TeV+2);	
		histoMCYieldPi07TeV->Draw("same,hist,c"); 	
		histoMCYieldPi02760GeVFinal->Scale(1e-1);
		DrawGammaSetMarker(histoMCYieldPi02760GeVFinal, markerStyleSpectrum2760GeV+1, markerSizePi0PP2760GeV*1.5, colorPi02760GeV+2, colorPi02760GeV+2);	
		histoMCYieldPi02760GeVFinal->Draw("same,hist,c"); 	
		histoMCYieldPi0900GeV->Scale(1e-2);
		DrawGammaSetMarker(histoMCYieldPi0900GeV, markerStyleSpectrum900GeV+4, markerSizePi0PP900GeV*1.5, colorPi0900GeV+2, colorPi0900GeV+2);	
		histoMCYieldPi0900GeV->DrawCopy("same,hist,c"); 	

		TLatex *labelMCPi0PP = new TLatex(0.5,0.73,"#pi^{0} #rightarrow #gamma #gamma (#rightarrow e^{+}e^{-} e^{+}e^{-})");
		SetStyleTLatex( labelMCPi0PP, 0.038,4);
		labelMCPi0PP->Draw();

		TLegend* legendMCYieldPi0 = new TLegend(0.35,0.78,0.93,0.93);
		legendMCYieldPi0->SetFillColor(0);
		legendMCYieldPi0->SetLineColor(0);
		legendMCYieldPi0->SetTextSize(0.035);
		legendMCYieldPi0->SetNColumns(3);
		
		legendMCYieldPi0->AddEntry(graphInvYieldPi0PCM7TeVSysErr,"Data","p");
		legendMCYieldPi0->AddEntry(histoMCYieldPi07TeV,"MC input","l");
		legendMCYieldPi0->AddEntry((TObject*)0,"pp #sqrt{#it{s}} = 7 TeV","");
		
		legendMCYieldPi0->AddEntry(graphInvYieldPi0PCM2760GeVSysErr,"Data","p");
		legendMCYieldPi0->AddEntry(histoMCYieldPi02760GeVFinal,"MC input","l");
		legendMCYieldPi0->AddEntry((TObject*)0,"pp #sqrt{#it{s}} = 2.76 TeV #times 10^{-1}","");
		
		legendMCYieldPi0->AddEntry(graphInvYieldPi0PCM900GeVSysErr,"Data","p"); 		
		legendMCYieldPi0->AddEntry(histoMCYieldPi0900GeV,"MC input","l"); 		
		legendMCYieldPi0->AddEntry((TObject*)0,"pp #sqrt{#it{s}} = 0.9 TeV #times 10^{-2}","");
		legendMCYieldPi0->Draw();
	
	canvasMCYieldPi0->SaveAs(Form("%s/MCYieldCompPi0.%s",outputDir.Data(),suffix.Data()));

      canvasMCYieldPi0->cd();
      canvasMCYieldPi0->SetLogx();
      histo2DMCYieldPi0->GetXaxis()->SetRangeUser(0.23,20);
      histo2DMCYieldPi0->Draw("copy");
   
//       graphInvYieldPi0Comb2760GeV= ScaleGraph(graphInvYieldPi0Comb2760GeV,1e1);
      DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0Comb2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);  
      graphInvYieldPi0Comb2760GeV->Draw("pe1,same");  
      
      DrawGammaSetMarker(histoMCYieldPi02760GeVAddSig, markerStylePP,markerSizePP, colorPi02760GeV-5, colorPi02760GeV-5);//,
      histoMCYieldPi02760GeVAddSig->SetLineStyle(2);
      histoMCYieldPi02760GeVAddSig->Draw("hist,c,same");
      
      histoMCYieldPi02760GeVFinal->Scale(1e1);
      histoMCYieldPi02760GeVFinal->Draw("same,hist,c");    
      DrawGammaSetMarker(histoMCYieldPi02760GeVAddSigWOWeighting, markerStylePP,markerSizePP, kGray+2, kGray+2);//, colorCombPbPb1020-5);
      histoMCYieldPi02760GeVAddSigWOWeighting->SetLineStyle(2);
      histoMCYieldPi02760GeVAddSigWOWeighting->Draw("hist,c,same");
      DrawGammaSetMarker(histoMCYieldPi02760GeVWOWeighting, markerStylePP,markerSizePP, kGray+2, kGray+2);//, colorCombPbPb1020-5);
      histoMCYieldPi02760GeVWOWeighting->Draw("hist,c,same");
      
      TLatex *labelMCPi0PPExpl = new TLatex(0.4,0.9,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}, pp #sqrt{#it{s}} = 2.76 TeV");
      SetStyleTLatex( labelMCPi0PPExpl, 0.038,4);
      labelMCPi0PPExpl->Draw();

      TLegend* legendMCYieldPi0Expl = new TLegend(0.15,0.13,0.73,0.25);
      legendMCYieldPi0Expl->SetFillColor(0);
      legendMCYieldPi0Expl->SetLineColor(0);
      legendMCYieldPi0Expl->SetTextSize(0.035);
      legendMCYieldPi0Expl->SetNColumns(2);
      legendMCYieldPi0Expl->AddEntry(graphInvYieldPi0Comb2760GeV,"Data","p");
      legendMCYieldPi0Expl->AddEntry(histoMCYieldPi02760GeVWOWeighting,"MC min Bias w/o weighting","l");
      legendMCYieldPi0Expl->AddEntry((TObject*)0,"","");
      legendMCYieldPi0Expl->AddEntry(histoMCYieldPi02760GeVFinal,"MC min Bias w weighting","l");
      legendMCYieldPi0Expl->AddEntry((TObject*)0,"","");
      legendMCYieldPi0Expl->AddEntry(histoMCYieldPi02760GeVAddSigWOWeighting,"MC add Sig w/o weighting","l");
      legendMCYieldPi0Expl->AddEntry((TObject*)0,"","");
      legendMCYieldPi0Expl->AddEntry(histoMCYieldPi02760GeVAddSig,"MC add Sig w weighting","l");
      legendMCYieldPi0Expl->Draw();
      
   canvasMCYieldPi0->SaveAs(Form("%s/MCYield2760GeVPi0Expl.%s",outputDir.Data(),suffix.Data()));
	
   TCanvas* canvasMCYieldEta = new TCanvas("canvasMCYieldEta","",200,10,1350*1.2,1350);  // gives the page size
      DrawGammaCanvasSettings( canvasMCYieldEta, 0.12, 0.02, 0.035, 0.09);
      TH2F * histo2DMCYieldEta;
//       canvasMCYieldEta->SetLogx();
      canvasMCYieldEta->SetLogy();
      histo2DMCYieldEta = new TH2F("histo2DMCYieldEta","histo2DMCYieldEta",1000,0.,10,2000,1.e-9,1e1  );
      SetStyleHistoTH2ForGraphs(histo2DMCYieldEta, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}",0.035,0.04, 0.035,0.04, 1.,1.3);
      histo2DMCYieldEta->Draw("copy");
   
      DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldSysEta7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV,widthLinesBoxes, kTRUE);  
      graphCorrectedYieldSysEta7TeV->Draw("e2,same");  
      graphCorrectedYieldSysEta2760GeV= ScaleGraph(graphCorrectedYieldSysEta2760GeV,1e-1);
      DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldSysEta2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV,widthLinesBoxes, kTRUE);
      graphCorrectedYieldSysEta2760GeV->Draw("e2,same");  
      graphCorrectedYieldSysEta900GeV= ScaleGraph(graphCorrectedYieldSysEta900GeV,1e-2);
      DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldSysEta900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV,widthLinesBoxes, kTRUE);
      graphCorrectedYieldSysEta900GeV->Draw("e2,same");   

      DrawGammaSetMarker(histoCorrectedYieldEta7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);  
      histoCorrectedYieldEta7TeV->Draw("pe1,same");  
      histoCorrectedYieldEta2760GeV->Scale(1e-1);
      DrawGammaSetMarker(histoCorrectedYieldEta2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);  
      histoCorrectedYieldEta2760GeV->Draw("pe1,same");  
      histoCorrectedYieldEta900GeV->Scale(1e-2);
      DrawGammaSetMarker(histoCorrectedYieldEta900GeV, markerStyleSpectrum900GeV, markerSizePi0PP900GeV*1.5, colorPi0900GeV, colorPi0900GeV); 
      histoCorrectedYieldEta900GeV->Draw("pe1,same");   

      SetStyleHisto(histoMCYieldEta7TeV, 1., 1, colorPi07TeV+2);  
      histoMCYieldEta7TeV->Draw("same,hist,c");    
      histoMCYieldEta2760GeVFinal->Scale(1e-1);
      DrawGammaSetMarker(histoMCYieldEta2760GeVFinal, markerStyleSpectrum2760GeV+1, markerSizePi0PP2760GeV*1.5, colorPi02760GeV+2, colorPi02760GeV+2);   
      histoMCYieldEta2760GeVFinal->Draw("same,hist,c");  
      histoMCYieldEta900GeV->Scale(1e-2);
      DrawGammaSetMarker(histoMCYieldEta900GeV, markerStyleSpectrum900GeV+4, markerSizePi0PP900GeV*1.5, colorPi0900GeV+2, colorPi0900GeV+2); 
      histoMCYieldEta900GeV->DrawCopy("same,hist,c");    

      TLatex *labelMCEtaPP = new TLatex(0.5,0.73,"#eta #rightarrow #gamma #gamma (#rightarrow e^{+}e^{-} e^{+}e^{-})");
      SetStyleTLatex( labelMCEtaPP, 0.038,4);
      labelMCEtaPP->Draw();

      TLegend* legendMCYieldEta = new TLegend(0.35,0.78,0.93,0.93);
      legendMCYieldEta->SetFillColor(0);
      legendMCYieldEta->SetLineColor(0);
      legendMCYieldEta->SetTextSize(0.035);
      legendMCYieldEta->SetNColumns(3);
      
      legendMCYieldEta->AddEntry(graphCorrectedYieldSysEta7TeV,"Data","p");
      legendMCYieldEta->AddEntry(histoMCYieldEta7TeV,"MC input","l");
      legendMCYieldEta->AddEntry((TObject*)0,"pp #sqrt{#it{s}} = 7 TeV","");
      
      legendMCYieldEta->AddEntry(graphCorrectedYieldSysEta2760GeV,"Data","p");
      legendMCYieldEta->AddEntry(histoMCYieldEta2760GeVFinal,"MC input","l");
      legendMCYieldEta->AddEntry((TObject*)0,"pp #sqrt{#it{s}} = 2.76 TeV #times 10^{-1}","");
      
      legendMCYieldEta->AddEntry(graphCorrectedYieldSysEta900GeV,"Data","p");       
      legendMCYieldEta->AddEntry(histoMCYieldEta900GeV,"MC input","l");       
      legendMCYieldEta->AddEntry((TObject*)0,"pp #sqrt{#it{s}} = 0.9 TeV #times 10^{-2}","");
      legendMCYieldEta->Draw();
   
   canvasMCYieldEta->SaveAs(Form("%s/MCYieldCompEta.%s",outputDir.Data(),suffix.Data()));

   
   
	TCanvas* canvasRawYieldPi0PbPb = new TCanvas("canvasRawYieldPi0PbPb","",200,10,1350*1.4,1350);  // gives the page size
		DrawGammaCanvasSettings( canvasRawYieldPi0PbPb, 0.12, 0.02, 0.035, 0.09);
		TH2F * histo2DRawPi0PbPb;
// 		canvasRawYieldPi0PbPb->SetLogx();
		canvasRawYieldPi0PbPb->SetLogy();
		histo2DRawPi0PbPb = new TH2F("histo2DRawPi0PbPb","histo2DRawPi0PbPb",1000,0.,14,2000,1.e-7,1e-1	);
		SetStyleHistoTH2ForGraphs(histo2DRawPi0PbPb, "p_{T} (GeV/c)","#frac{N_{raw}^{#pi^{0}}}{N_{evt} dp_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
		histo2DRawPi0PbPb->Draw("copy");
	
		DrawGammaSetMarker(histoRawYield6080, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);	
		histoRawYield6080->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoRawYield0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);	
		histoRawYield0005->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoRawYield0510, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0510, colorCombPbPb0510);	
		histoRawYield0510->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoRawYield1020, markerStylePbPb1020, markerSizePbPb1020, colorCombPbPb1020, colorCombPbPb1020);	
		histoRawYield1020->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoRawYield2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);	
		histoRawYield2040->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoRawYield4060, markerStylePbPb4060, markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060);	
		histoRawYield4060->DrawCopy("e1,same"); 	

		TLatex *labelRawPi0PbPb = new TLatex(0.3,0.9,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
		SetStyleTLatex( labelRawPi0PbPb, 0.038,4);
		labelRawPi0PbPb->Draw();

		
		TLegend* legendRawYieldPi0PbPb = new TLegend(0.55,0.63,0.93,0.93);
		legendRawYieldPi0PbPb->SetFillColor(0);
		legendRawYieldPi0PbPb->SetLineColor(0);
		legendRawYieldPi0PbPb->SetTextSize(0.035);
		legendRawYieldPi0PbPb->AddEntry(histoRawYield0005,"0-5% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldPi0PbPb->AddEntry(histoRawYield0510,"5-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldPi0PbPb->AddEntry(histoRawYield1020,"10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldPi0PbPb->AddEntry(histoRawYield2040,"20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldPi0PbPb->AddEntry(histoRawYield4060,"40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldPi0PbPb->AddEntry(histoRawYield6080,"60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldPi0PbPb->Draw();
	
	canvasRawYieldPi0PbPb->SaveAs(Form("%s/RawYieldCompPi0PbPb.%s",outputDir.Data(),suffix.Data()));

	TCanvas* canvasRawYieldEtaPbPb = new TCanvas("canvasRawYieldEtaPbPb","",200,10,1350*1.4,1350);  // gives the page size
		DrawGammaCanvasSettings( canvasRawYieldEtaPbPb, 0.12, 0.02, 0.035, 0.09);
		TH2F * histo2DRawEtaPbPb;
// 		canvasRawYieldEtaPbPb->SetLogx();
		canvasRawYieldEtaPbPb->SetLogy();
		histo2DRawEtaPbPb = new TH2F("histo2DRawEtaPbPb","histo2DRawEtaPbPb",1000,0.,10,2000,1.e-7,1e-1	);
		SetStyleHistoTH2ForGraphs(histo2DRawEtaPbPb, "p_{T} (GeV/c)","#frac{N_{raw}^{#eta}}{N_{evt} dp_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
		histo2DRawEtaPbPb->Draw("copy");
	
		DrawGammaSetMarker(histoEtaRawYield6080, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);	
		histoEtaRawYield6080->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoEtaRawYield0010, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);	
		histoEtaRawYield0010->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoEtaRawYield1020, markerStylePbPb1020, markerSizePbPb1020, colorCombPbPb1020, colorCombPbPb1020);	
		histoEtaRawYield1020->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoEtaRawYield2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);	
		histoEtaRawYield2040->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoEtaRawYield4060, markerStylePbPb4060, markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060);	
		histoEtaRawYield4060->DrawCopy("e1,same"); 	

		TLatex *labelRawEtaPbPb = new TLatex(0.3,0.9,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
		SetStyleTLatex( labelRawEtaPbPb, 0.038,4);
		labelRawEtaPbPb->Draw();

		TLegend* legendRawYieldEtaPbPb = new TLegend(0.55,0.63,0.93,0.93);
		legendRawYieldEtaPbPb->SetFillColor(0);
		legendRawYieldEtaPbPb->SetLineColor(0);
		legendRawYieldEtaPbPb->SetTextSize(0.035);
		legendRawYieldEtaPbPb->AddEntry(histoEtaRawYield0010,"0-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldEtaPbPb->AddEntry(histoEtaRawYield1020,"10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldEtaPbPb->AddEntry(histoEtaRawYield2040,"20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldEtaPbPb->AddEntry(histoEtaRawYield4060,"40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldEtaPbPb->AddEntry(histoEtaRawYield6080,"60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendRawYieldEtaPbPb->Draw();
	
	canvasRawYieldEtaPbPb->SaveAs(Form("%s/RawYieldCompEtaPbPb.%s",outputDir.Data(),suffix.Data()));

	
	TCanvas* canvasEffAllCent = new TCanvas("canvasEffAllCent","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasEffAllCent, 0.1, 0.02, 0.035, 0.09);
		TH2F * histo2DEffAllCent;
		histo2DEffAllCent = new TH2F("histo2DEffAllCent","histo2DEffAllCent",1000,0,16,2000,0.e-3,3.2e-3	);
// 		histo2DEffAllCent->GetXaxis()->SetRangeUser(0.,12.);
		SetStyleHistoTH2ForGraphs(histo2DEffAllCent, "p_{T} (GeV/c)","#epsilon_{reco, #pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
		histo2DEffAllCent->Draw("copy");

		DrawGammaSetMarker(histoTrueEffPtPi07TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);	
		histoTrueEffPtPi07TeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoTrueEffPtPi02760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);	
		histoTrueEffPtPi02760GeV->DrawCopy("pe1,same"); 	
		DrawGammaSetMarker(histoTrueEffPtPi0900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV);	
		histoTrueEffPtPi0900GeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoTrueEffiPt6080, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);	
		histoTrueEffiPt6080->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoTrueEffiPt0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);	
		histoTrueEffiPt0005->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoTrueEffiPt0510, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0510, colorCombPbPb0510);	
		histoTrueEffiPt0510->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoTrueEffiPt1020, markerStylePbPb1020, markerSizePbPb1020, colorCombPbPb1020, colorCombPbPb1020);	
		histoTrueEffiPt1020->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoTrueEffiPt2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);	
		histoTrueEffiPt2040->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoTrueEffiPt4060, markerStylePbPb4060, markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060);	
		histoTrueEffiPt4060->DrawCopy("e1,same"); 	
		cout << "bla"<<endl;
		TLegend* legendEffiCompAllCent = new TLegend(0.27,0.13,0.94,0.32);
		legendEffiCompAllCent->SetFillColor(0);
		legendEffiCompAllCent->SetLineColor(0);
		legendEffiCompAllCent->SetTextSize(0.027);
		legendEffiCompAllCent->SetNColumns(2);
		legendEffiCompAllCent->AddEntry(histoTrueEffiPt0005,"0-5% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendEffiCompAllCent->AddEntry(histoTrueEffPtPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
		legendEffiCompAllCent->AddEntry(histoTrueEffiPt0510,"5-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendEffiCompAllCent->AddEntry(histoTrueEffPtPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
		legendEffiCompAllCent->AddEntry(histoTrueEffiPt1020,"10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendEffiCompAllCent->AddEntry(histoTrueEffPtPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
		legendEffiCompAllCent->AddEntry(histoTrueEffiPt2040,"20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendEffiCompAllCent->AddEntry((TObject*)0, "","");
		legendEffiCompAllCent->AddEntry(histoTrueEffiPt4060,"40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendEffiCompAllCent->AddEntry((TObject*)0, "","");
		legendEffiCompAllCent->AddEntry(histoTrueEffiPt6080,"60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendEffiCompAllCent->Draw();
		
	canvasEffAllCent->SaveAs(Form("%s/EffCompSimpleAllCent.%s",outputDir.Data(),suffix.Data()));

		canvasEffAllCent->cd();
		histo2DEffAllCent->Draw("copy");

      DrawGammaSetMarker(histoEffi2760GeVWithSDD, 24  , markerSizePi0PP7TeV, kGreen+3, kGreen+3);   
      histoEffi2760GeVWithSDD->DrawCopy("pe1,same");
      
      DrawGammaSetMarker(histoEffi7TeV2760GeVCut, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);   
      histoEffi7TeV2760GeVCut->DrawCopy("pe1,same");
      
		DrawGammaSetMarker(histoEffiAddedSigPP, 24, markerSizePi0PP2760GeV, kRed+2, kRed+2);	
		histoEffiAddedSigPP->DrawCopy("pe1,same");
		DrawGammaSetMarker(histoEffiMinBiasPP, 26, markerSizePi0PP2760GeV, colorCombPbPb4060, colorCombPbPb4060);  
		histoEffiMinBiasPP->DrawCopy("pe1,same");
		DrawGammaSetMarker(histoTrueEffPtPi02760GeV, markerStyleSpectrum2760GeVMC-1, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);	
		histoTrueEffPtPi02760GeV->DrawCopy("pe1,same");
		
		TLegend* legendEffComperison = new TLegend(0.24,0.13,0.94,0.32);
		legendEffComperison->SetFillColor(0);
		legendEffComperison->SetLineColor(0);
		legendEffComperison->SetTextSize(0.027);
		legendEffComperison->AddEntry(histoTrueEffPtPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV, final","p");
		legendEffComperison->AddEntry(histoEffiAddedSigPP,"pp #sqrt{#it{s}} = 2.76 TeV, added signal","p");
		legendEffComperison->AddEntry(histoEffiMinBiasPP,"pp #sqrt{#it{s}} = 2.76 TeV, minBias","p");
      legendEffComperison->AddEntry(histoEffi7TeV2760GeVCut,"pp #sqrt{#it{s}} = 7 TeV, 2.76TeV cut","p");
      legendEffComperison->AddEntry(histoEffi2760GeVWithSDD,"pp #sqrt{#it{s}} = 2.76 TeV, with SDD","p");
		legendEffComperison->Draw();
	
	canvasEffAllCent->SaveAs(Form("%s/EffCompSimplePP.%s",outputDir.Data(),suffix.Data()));
	delete canvasEffAllCent;

	
	TCanvas* canvasAccepAllCent = new TCanvas("canvasAccepAllCent","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasAccepAllCent, 0.1, 0.02, 0.035, 0.09);
		TH2F * histo2DAccepAllEta;
		histo2DAccepAllEta = new TH2F("histo2DAccepAllEta","histo2DAccepAllEta",1000,0,8.,2000,0.4,1.02	);
		histo2DAccepAllEta->GetXaxis()->SetRangeUser(0.,12.);
		
		SetStyleHistoTH2ForGraphs(histo2DAccepAllEta, "p_{T} (GeV/c)","A_{#eta}",0.03,0.04, 0.03,0.04, 1.,1.2);
		histo2DAccepAllEta->Draw("copy");

		DrawGammaSetMarker(histoAccEta7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);	
		histoAccEta7TeV->DrawCopy("e1,same"); 	
		
		DrawGammaSetMarker(histoAccEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);	
		histoAccEta2760GeV->DrawCopy("e1,same"); 	
		cout << "bla"<<endl;
		DrawGammaSetMarker(histoAccEta900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV);	
		histoAccEta900GeV->DrawCopy("e1,same"); 	
		cout << "bla"<<endl;
		
		DrawGammaSetMarker(histoEtaAcceptPt6080, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);	
		histoEtaAcceptPt6080->DrawCopy("e1,same"); 	
		cout << "bla"<<endl;
		DrawGammaSetMarker(histoEtaAcceptPt0010, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);	
		histoEtaAcceptPt0010->DrawCopy("e1,same"); 	
		cout << "bla"<<endl;
		DrawGammaSetMarker(histoEtaAcceptPt1020, markerStylePbPb1020, markerSizePbPb1020, colorCombPbPb1020, colorCombPbPb1020);	
		histoEtaAcceptPt1020->DrawCopy("e1,same"); 	
		cout << "bla"<<endl;
		DrawGammaSetMarker(histoEtaAcceptPt2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);	
		histoEtaAcceptPt2040->DrawCopy("e1,same"); 	
		cout << "bla"<<endl;
		DrawGammaSetMarker(histoEtaAcceptPt4060, markerStylePbPb4060, markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060);	
		histoEtaAcceptPt4060->DrawCopy("e1,same"); 	
		cout << "bla"<<endl;

		TLegend* legendAcceptAllEta = new TLegend(0.34,0.13,0.93,0.43);
		legendAcceptAllEta->SetFillColor(0);
		legendAcceptAllEta->SetLineColor(0);
		legendAcceptAllEta->SetTextSize(0.03);
		legendAcceptAllEta->SetNColumns(2);
		legendAcceptAllEta->AddEntry(histoEtaAcceptPt0010,"0-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllEta->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllEta->AddEntry(histoEtaAcceptPt1020,"10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllEta->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllEta->AddEntry(histoEtaAcceptPt2040,"20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllEta->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllEta->AddEntry(histoEtaAcceptPt4060,"40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllEta->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllEta->AddEntry(histoEtaAcceptPt6080,"60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllEta->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllEta->AddEntry(histoAccEta7TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
		legendAcceptAllEta->AddEntry((TObject*)0, " |y| < 0.8","");
		legendAcceptAllEta->AddEntry(histoAccEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
		legendAcceptAllEta->AddEntry((TObject*)0, " |y| < 0.8","");
		legendAcceptAllEta->AddEntry(histoAccEta900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
		legendAcceptAllEta->AddEntry((TObject*)0, " |y| < 0.8","");
		legendAcceptAllEta->Draw();
	
	canvasAccepAllCent->SaveAs(Form("%s/AcceptanceCompEtaAll.%s",outputDir.Data(),suffix.Data()));
	delete canvasAccepAllCent;

	TCanvas* canvasAccepEta = new TCanvas("canvasAccepEta","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasAccepEta, 0.1, 0.02, 0.035, 0.09);
		TH2F * histo2DAccepAllCent;
		histo2DAccepAllCent = new TH2F("histo2DAccepAllCent","histo2DAccepAllCent",1000,0,16.,2000,0.6,1.02	);
// 		histo2DAccepAllCent->GetXaxis()->SetRangeUser(0.,12.);
		
		SetStyleHistoTH2ForGraphs(histo2DAccepAllCent, "p_{T} (GeV/c)","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.2);
		histo2DAccepAllCent->Draw("copy");

		DrawGammaSetMarker(histoAccPi07TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);	
		histoAccPi07TeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoAccPi02760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);	
		histoAccPi02760GeV->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoAccPi0900GeV, markerStyleSpectrum900GeVMC, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV);	
		histoAccPi0900GeV->DrawCopy("e1,same"); 	
		
		DrawGammaSetMarker(histoAcceptPt6080, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);	
		histoAcceptPt6080->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoAcceptPt0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);	
		histoAcceptPt0005->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoAcceptPt0510, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0510, colorCombPbPb0510);	
		histoAcceptPt0510->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoAcceptPt1020, markerStylePbPb1020, markerSizePbPb1020, colorCombPbPb1020, colorCombPbPb1020);	
		histoAcceptPt1020->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoAcceptPt2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);	
		histoAcceptPt2040->DrawCopy("e1,same"); 	
		DrawGammaSetMarker(histoAcceptPt4060, markerStylePbPb4060, markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060);	
		histoAcceptPt4060->DrawCopy("e1,same"); 	

		TLegend* legendAcceptAllCent = new TLegend(0.34,0.13,0.93,0.43);
		legendAcceptAllCent->SetFillColor(0);
		legendAcceptAllCent->SetLineColor(0);
		legendAcceptAllCent->SetTextSize(0.03);
		legendAcceptAllCent->SetNColumns(2);
		legendAcceptAllCent->AddEntry(histoAcceptPt0005,"0-5% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllCent->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllCent->AddEntry(histoAcceptPt0510,"5-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllCent->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllCent->AddEntry(histoAcceptPt1020,"10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllCent->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllCent->AddEntry(histoAcceptPt2040,"20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllCent->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllCent->AddEntry(histoAcceptPt4060,"40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllCent->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllCent->AddEntry(histoAcceptPt6080,"60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
		legendAcceptAllCent->AddEntry((TObject*)0, " |y| < 0.7","");
		legendAcceptAllCent->AddEntry(histoAccPi07TeV,"pp #sqrt{#it{s}} = 7 TeV","p");
		legendAcceptAllCent->AddEntry((TObject*)0, " |y| < 0.8","");
		legendAcceptAllCent->AddEntry(histoAccPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
		legendAcceptAllCent->AddEntry((TObject*)0, " |y| < 0.8","");
		legendAcceptAllCent->AddEntry(histoAccPi0900GeV,"pp #sqrt{#it{s}} = 0.9 TeV","p");
		legendAcceptAllCent->AddEntry((TObject*)0, " |y| < 0.8","");
		legendAcceptAllCent->Draw();
	
	canvasAccepEta->SaveAs(Form("%s/AcceptanceCompAllCent.%s",outputDir.Data(),suffix.Data()));
	delete canvasAccepEta;


	Double_t mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	Double_t mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();
	
	
	//********************************** Defintion of the Legend **************************************************	
	Double_t columnsLegend[4] 	= {0.,0.18,0.47,0.75};
	Double_t rowsLegend[6] 		= {0.88,0.75,0.57,0.4,0.22,0.05}; //with EMCAL {0.88,0.75,0.57,0.4,0.22,0.05};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumn	= 0.13;
	Size_t textSizeTopRow		= 0.13; 
	Size_t textSizeSecondRow 	= 0.11;
	//******************* Offsets ***********************
	Double_t offsetSystColumn 	= 0.15;
	Double_t offsetMarkerX		= 0.1;
	Double_t offsetMarkerY		= 0.05;
	Double_t offsetBoxSizeY		= 0.05;
	Double_t offsetFit			= 0.04;
	//****************** Scale factors ******************
	Double_t scaleWidthLine 		= 0.8;
		
	TCanvas * canvas6PartMassWidth = new TCanvas("canvas6PartMassWidth","",10,10,2400,1300);  // gives the page size		
	canvas6PartMassWidth->cd();
	DrawGammaCanvasSettings( canvas6PartMassWidth, 0.13, 0.0, 0.02, 0.09);
	
	TPad* pad6PartMassWidth1 = new TPad("pad6PartMassWidth1", "", 0., 0.54, 0.37, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth1, 0.16, 0.0, 0.02, 0.);
	pad6PartMassWidth1->Draw();
	TPad* pad6PartMassWidth2 = new TPad("pad6PartMassWidth2", "", 0., 0., 0.37, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth2, 0.16, 0.0, 0., 0.14);
	pad6PartMassWidth2->Draw();
	
	TPad* pad6PartMassWidth3 = new TPad("pad6PartMassWidth3", "", 0.37, 0.54, 0.685, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth3, 0.0, 0.0, 0.02, 0.);
	pad6PartMassWidth3->Draw();
	TPad* pad6PartMassWidth4 = new TPad("pad6PartMassWidth4", "", 0.37, 0., 0.685, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth4, 0.0, 0.0, 0., 0.14);
	pad6PartMassWidth4->Draw();

	TPad* pad6PartMassWidth5 = new TPad("pad6PartMassWidth5", "", 0.685, 0.54, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth5, 0.0, 0.02, 0.02, 0.);
	pad6PartMassWidth5->Draw();
	TPad* pad6PartMassWidth6 = new TPad("pad6PartMassWidth6", "", 0.685, 0., 1., 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth6, 0.0, 0.02, 0., 0.14);
	pad6PartMassWidth6->Draw();

	TPad* padMassLegend1 = new TPad("padMassLegend1", "", 0.07, 0.09, 0.20, 0.185,-1, -1, -2);
	DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
	padMassLegend1->Draw();
	
	TPad* padFWHMLegend1 = new TPad("padFWHMLegend1", "", 0.07, 0.85, 0.28, 0.94,-1, -1, -2);
	DrawGammaPadSettings( padFWHMLegend1, 0., 0., 0., 0.);
	padFWHMLegend1->Draw();

	
	TH2D *histo2DPi0FWHM;
	histo2DPi0FWHM = new TH2D("histo2DPi0FWHM", "histo2DPi0FWHM", 20,0.35,20. ,1000.,-30,40);
	SetStyleHistoTH2ForGraphs(histo2DPi0FWHM, "p_{T} (GeV/c)","FWHM/2.36 (MeV/c^{2})", 0.035,0.05, 0.065,0.08, 1,1, 515, 504); 
	histo2DPi0FWHM->GetYaxis()->SetRangeUser(-0.5,10);
	histo2DPi0FWHM->GetYaxis()->SetLabelOffset(0.01);

	TH2D *histo2DPi0Mass;
	histo2DPi0Mass = new TH2D("histo2DPi0Mass", "histo2DPi0Mass", 20,0.35,20. ,1000.,125.,150);
	SetStyleHistoTH2ForGraphs(histo2DPi0Mass, "p_{T} (GeV/c)","m_{#pi^{0}} (MeV/c^{2})", 0.062,0.076, 0.062,0.078, 0.8,1, 515, 510); 
	histo2DPi0Mass->GetYaxis()->SetRangeUser(130.,140.5);
	histo2DPi0Mass->GetXaxis()->SetLabelOffset(-0.02);
	histo2DPi0Mass->GetYaxis()->SetLabelOffset(0.01);
	pad6PartMassWidth1->cd();
 	pad6PartMassWidth1->SetLogx();
	histo2DPi0FWHM->DrawCopy();
			
	DrawGammaSetMarker(histoPCMWidthDataPP2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);
	histoPCMWidthDataPP2760GeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMWidthMCPP2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);
	histoPCMWidthMCPP2760GeV->DrawCopy("same,e1,p"); 

	TLatex *labelMassPi0PP = new TLatex(0.2,0.9,collisionSystemPP2760GeV.Data());
	SetStyleTLatex( labelMassPi0PP, 0.062,4);
	labelMassPi0PP->Draw();
	TLatex *labelLegendAMass = new TLatex(0.92,0.88,"a)");
	SetStyleTLatex( labelLegendAMass, 0.08,4);
	labelLegendAMass->Draw();
	
	pad6PartMassWidth1->Update();
	pad6PartMassWidth2->cd();
 	pad6PartMassWidth2->SetLogx();
 	histo2DPi0Mass->DrawCopy();
	
	DrawGammaSetMarker(histoPCMMassDataPP2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);					 
	histoPCMMassDataPP2760GeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMMassMCPP2760GeV , markerStyleSpectrum2760GeVMC , markerSizePi0PP2760GeV, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);					 
	histoPCMMassMCPP2760GeV->DrawCopy("same,e1,p"); 

	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,1.);
	TLatex *labelLegendDMass = new TLatex(0.92,0.9,"d)");
	SetStyleTLatex( labelLegendDMass, 0.075,4);
	labelLegendDMass->Draw();
	
	pad6PartMassWidth2->Update();

	pad6PartMassWidth5->cd();
	pad6PartMassWidth5->SetLogx();
	histo2DPi0FWHM->DrawCopy();
			
	DrawGammaSetMarker(histoPCMWidthData0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
	histoPCMWidthData0005->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMWidthMC0005, markerStylePbPb0005MC, markerSizePbPb0005, colorCombMCPbPb0005, colorCombMCPbPb0005);
	histoPCMWidthMC0005->DrawCopy("same,e1,p"); 
	
	TLatex *labelMassPi0PbPb0005 = new TLatex(0.05,0.9,collisionSystemPbPb0005.Data());
	SetStyleTLatex( labelMassPi0PbPb0005, 0.062,4);
	labelMassPi0PbPb0005->Draw();
	TLatex *labelLegendBMass = new TLatex(0.89,0.88,"c)");
	SetStyleTLatex( labelLegendBMass, 0.08,4);
	labelLegendBMass->Draw();
	
	pad6PartMassWidth5->Update();
	pad6PartMassWidth6->cd();
 	pad6PartMassWidth6->SetLogx();
 	histo2DPi0Mass->DrawCopy();
	
	DrawGammaSetMarker(histoPCMMassData0005 , markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);					 
	histoPCMMassData0005->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMMassMC0005 , markerStylePbPb0005MC , markerSizePbPb0005, colorCombMCPbPb0005, colorCombMCPbPb0005);					 
	histoPCMMassMC0005->DrawCopy("same,e1,p"); 
	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,1);
	TLatex *labelLegendFMass = new TLatex(0.89,0.9,"f)");
	SetStyleTLatex( labelLegendFMass, 0.075,4);
	labelLegendFMass->Draw();
	
	pad6PartMassWidth6->Update();

	pad6PartMassWidth3->cd();
 	pad6PartMassWidth3->SetLogx();
	histo2DPi0FWHM->DrawCopy();
		
	DrawGammaSetMarker(histoPCMWidthData6080, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
	histoPCMWidthData6080->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMWidthMC6080, markerStylePbPb6080MC, markerSizePbPb6080, colorCombMCPbPb6080, colorCombMCPbPb6080);
	histoPCMWidthMC6080->DrawCopy("same,e1,p"); 
	
	TLatex *labelMassPi0PbPb6080 = new TLatex(0.05,0.9,collisionSystemPbPb6080.Data());
	SetStyleTLatex( labelMassPi0PbPb6080, 0.062,4);
	labelMassPi0PbPb6080->Draw();
	TLatex *labelLegendCMass = new TLatex(0.91,0.88,"b)");
	SetStyleTLatex( labelLegendCMass, 0.08,4);
	labelLegendCMass->Draw();

	pad6PartMassWidth3->Update();
	pad6PartMassWidth4->cd();
 	pad6PartMassWidth4->SetLogx();
 	histo2DPi0Mass->DrawCopy();
		
	DrawGammaSetMarker(histoPCMMassData6080 , markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);					 
	histoPCMMassData6080->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMMassMC6080 , markerStylePbPb0005MC , markerSizePbPb6080, colorCombMCPbPb6080, colorCombMCPbPb6080);					 
	histoPCMMassMC6080->DrawCopy("same,e1,p"); 
	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,1);
	TLatex *labelLegendEMass = new TLatex(0.91,0.9,"e)");
	SetStyleTLatex( labelLegendEMass, 0.075,4);
	labelLegendEMass->Draw();

	pad6PartMassWidth4->Update();

	canvas6PartMassWidth->Update();	
	canvas6PartMassWidth->SaveAs(Form("%s/MassWidth_Parted.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartMassWidth1;	
	delete pad6PartMassWidth2;	
	delete pad6PartMassWidth3;	
	delete pad6PartMassWidth4;	
	delete canvas6PartMassWidth;	

	TCanvas * canvas6PartEtaMassWidth = new TCanvas("canvas6PartEtaMassWidth","",10,10,2400,1300);  // gives the page size		
	canvas6PartEtaMassWidth->cd();
	DrawGammaCanvasSettings( canvas6PartEtaMassWidth, 0.13, 0.0, 0.02, 0.09);
	
	TPad* pad6PartEtaMassWidth1 = new TPad("pad6PartEtaMassWidth1", "", 0., 0.54, 0.37, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth1, 0.16, 0.0, 0.02, 0.);
	pad6PartEtaMassWidth1->Draw();
	TPad* pad6PartEtaMassWidth2 = new TPad("pad6PartEtaMassWidth2", "", 0., 0., 0.37, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth2, 0.16, 0.0, 0., 0.14);
	pad6PartEtaMassWidth2->Draw();
	
	TPad* pad6PartEtaMassWidth3 = new TPad("pad6PartEtaMassWidth3", "", 0.37, 0.54, 0.685, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth3, 0.0, 0.0, 0.02, 0.);
	pad6PartEtaMassWidth3->Draw();
	TPad* pad6PartEtaMassWidth4 = new TPad("pad6PartEtaMassWidth4", "", 0.37, 0., 0.685, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth4, 0.0, 0.0, 0., 0.14);
	pad6PartEtaMassWidth4->Draw();

	TPad* pad6PartEtaMassWidth5 = new TPad("pad6PartEtaMassWidth5", "", 0.685, 0.54, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth5, 0.0, 0.03, 0.02, 0.);
	pad6PartEtaMassWidth5->Draw();
	TPad* pad6PartEtaMassWidth6 = new TPad("pad6PartEtaMassWidth6", "", 0.685, 0., 1., 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth6, 0.0, 0.03, 0., 0.14);
	pad6PartEtaMassWidth6->Draw();
	
	TH2D *histo2DEtaFWHM;
	histo2DEtaFWHM = new TH2D("histo2DEtaFWHM", "histo2DEtaFWHM", 20,0.35,12. ,1000.,-30,40);
	SetStyleHistoTH2ForGraphs(histo2DEtaFWHM, "p_{T} (GeV/c)","FWHM/2.36 (MeV/c^{2})", 0.035,0.05, 0.065,0.08, 1,1, 515, 504); 
	histo2DEtaFWHM->GetYaxis()->SetRangeUser(-0.5,20);
	histo2DEtaFWHM->GetYaxis()->SetLabelOffset(0.01);

	TH2D *histo2DEtaMass;
	histo2DEtaMass = new TH2D("histo2DEtaMass", "histo2DEtaMass", 20,0.35,12. ,1000.,540.,559.5);
	SetStyleHistoTH2ForGraphs(histo2DEtaMass, "p_{T} (GeV/c)","m_{#eta} (MeV/c^{2})", 0.062,0.076, 0.062,0.078, 0.8,1, 515, 504); 
// 	histo2DEtaMass->GetYaxis()->SetRangeUser(130.,140.5);
	histo2DEtaMass->GetXaxis()->SetLabelOffset(-0.02);
	histo2DEtaMass->GetYaxis()->SetLabelOffset(0.01);
	pad6PartEtaMassWidth1->cd();
 	pad6PartEtaMassWidth1->SetLogx();
	histo2DEtaFWHM->DrawCopy();
			
	DrawGammaSetMarker(histoEtaWidthData7TeV, markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);
	histoEtaWidthData7TeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoEtaWidthMC7TeV, markerStyleSpectrum7TeVMC, markerSizePi0PP7TeV*1.5, colorMCPhojetPP7TeV, colorMCPhojetPP7TeV);
	histoEtaWidthMC7TeV->DrawCopy("same,e1,p"); 

	TLatex *labelMassEtaPP = new TLatex(0.2,0.9,collisionSystemPP7TeV.Data());
	SetStyleTLatex( labelMassEtaPP, 0.062,4);
	labelMassEtaPP->Draw();
	labelLegendAMass->Draw();
	
	pad6PartEtaMassWidth1->Update();
	pad6PartEtaMassWidth2->cd();
 	pad6PartEtaMassWidth2->SetLogx();
 	histo2DEtaMass->DrawCopy();
	
	DrawGammaSetMarker(histoEtaMassData7TeV , markerStyleSpectrum7TeV, markerSizePi0PP7TeV*1.5, colorPi07TeV, colorPi07TeV);
	histoEtaMassData7TeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoEtaMassMC7TeV , markerStyleSpectrum7TeVMC , markerSizePi0PP7TeV*1.5, colorMCPhojetPP7TeV, colorMCPhojetPP7TeV);					 
	histoEtaMassMC7TeV->DrawCopy("same,e1,p"); 

	DrawGammaLines(0.35, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,1.);
	labelLegendDMass->Draw();
	
	pad6PartEtaMassWidth2->Update();

	pad6PartEtaMassWidth5->cd();
	pad6PartEtaMassWidth5->SetLogx();
	histo2DEtaFWHM->DrawCopy();
			
	DrawGammaSetMarker(histoEtaWidthData6080, markerStylePbPb6080, markerSizePbPb6080*1.5, colorCombPbPb6080, colorCombPbPb6080);
	histoEtaWidthData6080->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoEtaWidthMC6080, markerStylePbPb6080MC, markerSizePbPb6080*1.5, colorCombMCPbPb6080, colorCombMCPbPb6080);
	histoEtaWidthMC6080->DrawCopy("same,e1,p"); 
	
	TLatex *labelMassEtaPbPb6080 = new TLatex(0.05,0.9,collisionSystemPbPb6080.Data());
	SetStyleTLatex( labelMassEtaPbPb6080, 0.062,4);
	labelMassEtaPbPb6080->Draw();
	labelLegendBMass->Draw();
	
	pad6PartEtaMassWidth5->Update();
	pad6PartEtaMassWidth6->cd();
 	pad6PartEtaMassWidth6->SetLogx();
 	histo2DEtaMass->DrawCopy();
	
	DrawGammaSetMarker(histoEtaMassData6080 , markerStylePbPb6080, markerSizePbPb6080*1.5, colorCombPbPb6080, colorCombPbPb6080);					 
	histoEtaMassData6080->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoEtaMassMC6080 , markerStylePbPb6080MC , markerSizePbPb6080*1.5, colorCombMCPbPb6080, colorCombMCPbPb6080);				 
	histoEtaMassMC6080->DrawCopy("same,e1,p"); 
	DrawGammaLines(0.35, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,1);
	labelLegendFMass->Draw();
	
	pad6PartEtaMassWidth6->Update();

	pad6PartEtaMassWidth3->cd();
 	pad6PartEtaMassWidth3->SetLogx();
	histo2DEtaFWHM->DrawCopy();
		
	DrawGammaSetMarker(histoEtaWidthData2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);
	histoEtaWidthData2760GeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoEtaWidthMC2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV*1.5, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);
	histoEtaWidthMC2760GeV->DrawCopy("same,e1,p"); 
	
		TLatex *labelMassEtaPP2760GeV = new TLatex(0.05,0.9,collisionSystemPP2760GeV.Data());
	SetStyleTLatex( labelMassEtaPP2760GeV, 0.062,4);
	labelMassEtaPP2760GeV->Draw();
	labelLegendCMass->Draw();

	pad6PartEtaMassWidth3->Update();
	pad6PartEtaMassWidth4->cd();
 	pad6PartEtaMassWidth4->SetLogx();
 	histo2DEtaMass->DrawCopy();
		
	DrawGammaSetMarker(histoEtaMassData2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV*1.5, colorPi02760GeV, colorPi02760GeV);
	histoEtaMassData2760GeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoEtaMassMC2760GeV , markerStyleSpectrum2760GeVMC , markerSizePi0PP2760GeV*1.5, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);					 
	histoEtaMassMC2760GeV->DrawCopy("same,e1,p"); 

	
	DrawGammaLines(0.35, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,1);
	labelLegendEMass->Draw();

	pad6PartEtaMassWidth4->Update();

	canvas6PartEtaMassWidth->Update();	
	canvas6PartEtaMassWidth->SaveAs(Form("%s/MassWidth_Eta.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartEtaMassWidth1;	
	delete pad6PartEtaMassWidth2;	
	delete pad6PartEtaMassWidth3;	
	delete pad6PartEtaMassWidth4;	
	delete canvas6PartEtaMassWidth;	

	TFile* fileCorrectionsSecPi0 = new TFile("ExternalInput/PCM/SecondaryFractionHistogramms7TeV.root");
	TH1D*	histoDefaultTrueSecFracMeson = 				(TH1D*)fileCorrectionsSecPi0->Get("TrueSecFrac");
	TH1D*	histoDefaultTrueSecFracFromK0SMeson = 		(TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracFromK0S");
	TF1* fitDefaultSecFrac = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
	fitDefaultSecFrac->SetRange(0.3,16.);
	TFitResultPtr resultSecFrac = histoDefaultTrueSecFracMeson->Fit(fitDefaultSecFrac,"SINRME+","",0.3,16.);
	TF1* fitDefaultSecFracFromK0 = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
	fitDefaultSecFracFromK0->SetRange(0.3,16.);
	TFitResultPtr resultSecFracFromK0 = histoDefaultTrueSecFracFromK0SMeson->Fit(fitDefaultSecFracFromK0,"SINRME+","",0.3,16.);

	
		TCanvas* canvasSecFrac = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasSecFrac, 0.09, 0.02, 0.04, 0.09);
//       canvasSecFrac->SetLogy(1);	
						
		DrawAutoGammaMesonHistos( histoDefaultTrueSecFracMeson, 
				"", "p_{T} (GeV/c)", "r_{X} = #frac{X->#pi^{0}}{#pi^{0}}", 
				kTRUE, 1.5, 0, kFALSE,
				kFALSE, 0., 0.7, 
				kFALSE, 0., 10.);
		histoDefaultTrueSecFracMeson->GetYaxis()->SetTitleOffset(0.9);		
		DrawGammaSetMarker(histoDefaultTrueSecFracMeson, 20, 1.2, kBlack, kBlack);										 
		fitDefaultSecFrac->SetLineWidth(1.);
		fitDefaultSecFrac->SetLineColor(kBlack);	
		DrawGammaSetMarker(histoDefaultTrueSecFracFromK0SMeson, 25, 1.2, kBlue, kBlue);										 
		fitDefaultSecFracFromK0->SetLineWidth(1.);
		fitDefaultSecFracFromK0->SetLineColor(kBlue);
	
		TLatex *labelFrac = new TLatex(0.15,0.85,collisionSystemPP7TeV.Data());
		SetStyleTLatex( labelFrac, 0.05,4);
// 		labelFrac->Draw();
		
		
		histoDefaultTrueSecFracMeson->DrawCopy("e1"); 	
		fitDefaultSecFrac->Draw("same");
		histoDefaultTrueSecFracFromK0SMeson->DrawCopy("e1,same"); 	
		fitDefaultSecFracFromK0->Draw("same");
		labelFrac->Draw();
		
		TLegend* legendSecFrac = new TLegend(0.7,0.78,0.96,0.94);
		legendSecFrac->SetTextSize(0.03);			
		legendSecFrac->SetFillColor(0);
		legendSecFrac->SetLineColor(0);
		legendSecFrac->AddEntry(histoDefaultTrueSecFracMeson,"X= All Particles");
		legendSecFrac->AddEntry(fitDefaultSecFrac,"fit on total fraction");
		legendSecFrac->AddEntry(histoDefaultTrueSecFracFromK0SMeson,"X = K_{s}^{0}");
		legendSecFrac->AddEntry(fitDefaultSecFracFromK0,"fit on fraction from K_{s}^{0}");
		legendSecFrac->Draw();

/*     if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
*/     canvasSecFrac->Update();
		
		canvasSecFrac->SaveAs(Form("%s/FracSecondaries.%s",outputDir.Data(),suffix.Data()));
		delete canvasSecFrac;

	cout << "here" << endl;
	TCanvas* canvasEtaSpectraAllTogether = new TCanvas("canvasEtaSpectraAllTogether","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasEtaSpectraAllTogether,  0.13, 0.01, 0.015, 0.08);
 	
	canvasEtaSpectraAllTogether->SetLogy();
  	canvasEtaSpectraAllTogether->SetLogx();
	TH2F * histo2DEtaSpectraAll;
	histo2DEtaSpectraAll = new TH2F("histo2DEtaSpectraAll","histo2DEtaSpectraAll",1000,0.3,11.,1000,1e-8,2e2	);
	SetStyleHistoTH2ForGraphs(histo2DEtaSpectraAll, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}
 	histo2DEtaSpectraAll->GetXaxis()->SetLabelOffset(-0.01);
	histo2DEtaSpectraAll->GetYaxis()->SetLabelOffset(0.01);
	histo2DEtaSpectraAll->DrawCopy(); 
	graphCorrectedYieldSysEta2760GeV = ScaleGraph(graphCorrectedYieldSysEta2760GeV,1e1);
	DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldSysEta2760GeV, markerStylePP,markerSizePP, kBlack , kBlack, widthLinesBoxes, kTRUE);//, colorCombPbPb1020-5);
	graphCorrectedYieldSysEta2760GeV->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSys0010, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005, widthLinesBoxes, kTRUE);//, colorCombPbPb0510-5);
	graphPCMEtaCorrectedSpecSys0010->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSys1020, markerStylePbPb1020,markerSizePbPb1020, colorCombPbPb1020 , colorCombPbPb1020, widthLinesBoxes, kTRUE);//, colorCombPbPb1020-5);
	graphPCMEtaCorrectedSpecSys1020->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSys2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, widthLinesBoxes, kTRUE);//, colorCombPbPb2040-5);
	graphPCMEtaCorrectedSpecSys2040->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSys4060, markerStylePbPb4060,markerSizePbPb4060, colorCombPbPb4060 , colorCombPbPb4060, widthLinesBoxes, kTRUE);//, colorCombPbPb4060-5);
	graphPCMEtaCorrectedSpecSys4060->Draw("E2same");
	
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSys6080, markerStylePbPb6080,markerSizePbPb6080, colorCombPbPb6080 , colorCombPbPb6080, widthLinesBoxes, kTRUE);//, colorCombPbPb6080-5);
	graphPCMEtaCorrectedSpecSys6080->Draw("E2same");
	cout << "here" << endl;
   histoCorrectedYieldEta2760GeV->Scale(1e1);
	DrawGammaSetMarker(histoCorrectedYieldEta2760GeV, markerStylePP,markerSizePP, kBlack , kBlack);
	histoCorrectedYieldEta2760GeV->Draw("p,same,e1");
	histoEtaFromCocktailpp2760GeV->SetLineWidth(2.);
	histoEtaFromCocktailpp2760GeV->SetLineColor(kBlack);
	histoEtaFromCocktailpp2760GeV->Draw("hist, l,same");

	DrawGammaSetMarker(histoPCMEtaCorrectedSpec0010, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
	histoPCMEtaCorrectedSpec0010->Draw("p,same,e1");
	histoEtaFromCocktail0010->SetLineWidth(2.);
	histoEtaFromCocktail0010->SetLineColor(colorCombPbPb0005);
	DrawGammaSetMarker(histoEtaFromCocktail0010, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
// 	histoEtaFromCocktail0010->SetFillColor(colorCombPbPb0005);
  	histoEtaFromCocktail0010->Draw("hist,l,same");
// 	histoEtaFromCocktail0510->SetLineWidth(2.);
// 	histoEtaFromCocktail0510->SetLineColor(colorCombPbPb0510);
// 	DrawGammaSetMarker(histoEtaFromCocktail0510, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0510 , colorCombPbPb0510);
// 	histoEtaFromCocktail0510->Draw("hist, l,same");
	
	DrawGammaSetMarker(histoPCMEtaCorrectedSpec1020, markerStylePbPb1020,markerSizePbPb1020, colorCombPbPb1020 , colorCombPbPb1020);
	histoPCMEtaCorrectedSpec1020->Draw("p,same,e1");
	histoEtaFromCocktail1020	->SetLineWidth(2.);
	histoEtaFromCocktail1020->SetLineColor(colorCombPbPb1020);
// 	DrawGammaSetMarker(histoEtaFromCocktail1020, markerStylePbPb1020,markerSizePbPb1020, colorCombPbPb1020 , colorCombPbPb1020);
	histoEtaFromCocktail1020->Draw("hist,l,same");
	
	DrawGammaSetMarker(histoPCMEtaCorrectedSpec2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
	histoPCMEtaCorrectedSpec2040->Draw("p,same,e1");
	histoEtaFromCocktail2040->SetLineWidth(2.);
	histoEtaFromCocktail2040->SetLineColor(colorCombPbPb2040 );
// 	DrawGammaSetMarker(histoEtaFromCocktail2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
	histoEtaFromCocktail2040->Draw("hist,l,same");
	
	DrawGammaSetMarker(histoPCMEtaCorrectedSpec4060, markerStylePbPb4060,markerSizePbPb4060, colorCombPbPb4060 , colorCombPbPb4060);
	histoPCMEtaCorrectedSpec4060->Draw("p,same,e1");
	histoEtaFromCocktail4060->SetLineWidth(2.);
	histoEtaFromCocktail4060->SetLineColor(colorCombPbPb4060 );
// 	DrawGammaSetMarker(histoEtaFromCocktail4060, markerStylePbPb4060,markerSizePbPb4060, colorCombPbPb4060 , colorCombPbPb4060);
	histoEtaFromCocktail4060->Draw("hist,l,same");

	DrawGammaSetMarker(histoPCMEtaCorrectedSpec6080, markerStylePbPb6080,markerSizePbPb6080, colorCombPbPb6080 , colorCombPbPb6080);
	histoPCMEtaCorrectedSpec6080->Draw("same");
	histoEtaFromCocktail6080->SetLineWidth(2.);
	histoEtaFromCocktail6080->SetLineColor(colorCombPbPb6080 );
// 	DrawGammaSetMarker(histoEtaFromCocktail6080, markerStylePbPb6080,markerSizePbPb6080, colorCombPbPb6080 , colorCombPbPb6080);
	histoEtaFromCocktail6080->Draw("hist,l,same");
	cout << "here" << endl;


	TLatex *labelSpectraEtaLabelPbPb = new TLatex(0.6,0.89,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelSpectraEtaLabelPbPb, 0.05	,4);
	labelSpectraEtaLabelPbPb->Draw();
	
	TLegend* legendSpectra = new TLegend(0.16,0.09,0.83,0.3);
	legendSpectra->SetFillColor(0);
	legendSpectra->SetLineColor(0);
	legendSpectra->SetTextSize(0.025);
	legendSpectra->SetNColumns(3);
	legendSpectra->SetMargin(0.2);
	legendSpectra->AddEntry(histoEtaFromCocktail0010,"m_{T} scaled","l");
	legendSpectra->AddEntry(graphPCMEtaCorrectedSpecSys0010,"measured","pf");
	legendSpectra->AddEntry((TObject*)0, collisionSystemPbPb0010.Data(),"");
	legendSpectra->AddEntry(histoEtaFromCocktail1020,"m_{T} scaled","l");
	legendSpectra->AddEntry(graphPCMEtaCorrectedSpecSys1020,"measured","pf");
	legendSpectra->AddEntry((TObject*)0, collisionSystemPbPb1020.Data(),"");
	legendSpectra->AddEntry(histoEtaFromCocktail2040,"m_{T} scaled","l");
	legendSpectra->AddEntry(graphPCMEtaCorrectedSpecSys2040,"measured","pf");
	legendSpectra->AddEntry((TObject*)0, collisionSystemPbPb2040.Data(),"");
	legendSpectra->AddEntry(histoEtaFromCocktail4060,"m_{T} scaled","l");
	legendSpectra->AddEntry(graphPCMEtaCorrectedSpecSys4060,"measured","pf");
	legendSpectra->AddEntry((TObject*)0, collisionSystemPbPb4060.Data(),"");
	legendSpectra->AddEntry(histoEtaFromCocktail6080,"m_{T} scaled","l");
	legendSpectra->AddEntry(graphPCMEtaCorrectedSpecSys6080,"measured","pf");
	legendSpectra->AddEntry((TObject*)0, collisionSystemPbPb6080.Data(),"");
	legendSpectra->AddEntry(histoEtaFromCocktailpp2760GeV,"m_{T} scaled","l");
	legendSpectra->AddEntry(graphCorrectedYieldSysEta2760GeV,"measured","pf");
	legendSpectra->AddEntry((TObject*)0, collisionSystemPP2760GeV.Data(),"");

// 	legendSpectra->AddEntry(graphInvSectionCombSysPi02760GeVPlot,collisionSystemPP.Data(),"pf");
	
// 	legendSpectra->AddEntry(fitInvCrossSectionPi0Comb2760GeV,"Tsallis Fit","l");
// 	legendSpectra->AddEntry(fitInvCrossSectionPi0Comb2760GeVPow,"Powerlaw Fit","l");
	legendSpectra->Draw();
	
	canvasEtaSpectraAllTogether->Update();
	canvasEtaSpectraAllTogether->Print(Form("%s/Eta_Spectra_All.%s",outputDir.Data(),suffix.Data()));

   
   TCanvas* canvasPi0SpectraAllTogether = new TCanvas("canvasPi0SpectraAllTogether","",200,10,1200,1100);  // gives the page size
   DrawGammaCanvasSettings( canvasPi0SpectraAllTogether,  0.13, 0.01, 0.015, 0.08);
   
   canvasPi0SpectraAllTogether->SetLogy();
   canvasPi0SpectraAllTogether->SetLogx();
   TH2F * histo2DPi0SpectraAll;
   histo2DPi0SpectraAll = new TH2F("histo2DPi0SpectraAll","histo2DPi0SpectraAll",1000,0.23,15.,1000,1e-7,2e3 );
   SetStyleHistoTH2ForGraphs(histo2DPi0SpectraAll, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}
   histo2DPi0SpectraAll->GetXaxis()->SetLabelOffset(-0.01);
   histo2DPi0SpectraAll->GetYaxis()->SetLabelOffset(0.01);
   histo2DPi0SpectraAll->DrawCopy(); 
   
   graphCorrectedYieldSysPi02760GeV->RemovePoint(17);
   DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldSysPi02760GeV, markerStylePP,markerSizePP, kBlack , kBlack, widthLinesBoxes, kTRUE);//, colorCombPbPb1020-5);
   graphCorrectedYieldSysPi02760GeV->Draw("E2same");
   graphPCMPi0CorrectedSpecSys0005->RemovePoint(0);
   graphPCMPi0CorrectedSpecSys0005->RemovePoint(16);
   graphPCMPi0CorrectedSpecSys0005->RemovePoint(16);
   
   graphPCMPi0CorrectedSpecSys0005 = ScaleGraph(graphPCMPi0CorrectedSpecSys0005,4.);
   DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSys0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005, widthLinesBoxes, kTRUE);//, colorCombPbPb0510-5);
   graphPCMPi0CorrectedSpecSys0005->Draw("E2same");
   graphPCMPi0CorrectedSpecSys0510->RemovePoint(0);
   graphPCMPi0CorrectedSpecSys0510->RemovePoint(16);
   graphPCMPi0CorrectedSpecSys0510->RemovePoint(16);
   graphPCMPi0CorrectedSpecSys0510 = ScaleGraph(graphPCMPi0CorrectedSpecSys0510,2.);
   DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSys0510, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0510 , colorCombPbPb0510, widthLinesBoxes, kTRUE);//, colorCombPbPb0510-5);
   graphPCMPi0CorrectedSpecSys0510->Draw("E2same");
   graphPCMPi0CorrectedSpecSys1020->RemovePoint(17);
   graphPCMPi0CorrectedSpecSys1020->RemovePoint(17);
   DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSys1020, markerStylePbPb1020,markerSizePbPb1020, colorCombPbPb1020 , colorCombPbPb1020, widthLinesBoxes, kTRUE);//, colorCombPbPb1020-5);
   graphPCMPi0CorrectedSpecSys1020->Draw("E2same");
   graphPCMPi0CorrectedSpecSys2040->RemovePoint(17);
   graphPCMPi0CorrectedSpecSys2040->RemovePoint(17);
   
   DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSys2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, widthLinesBoxes, kTRUE);//, colorCombPbPb2040-5);
   graphPCMPi0CorrectedSpecSys2040->Draw("E2same");
   
   graphPCMPi0CorrectedSpecSys4060->RemovePoint(17);
   graphPCMPi0CorrectedSpecSys4060->RemovePoint(17);
   DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSys4060, markerStylePbPb4060,markerSizePbPb4060, colorCombPbPb4060 , colorCombPbPb4060, widthLinesBoxes, kTRUE);//, colorCombPbPb4060-5);
   graphPCMPi0CorrectedSpecSys4060->Draw("E2same");
   
   graphPCMPi0CorrectedSpecSys6080->RemovePoint(17);
   graphPCMPi0CorrectedSpecSys6080->RemovePoint(17);
   DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSys6080, markerStylePbPb6080,markerSizePbPb6080, colorCombPbPb6080 , colorCombPbPb6080, widthLinesBoxes, kTRUE);//, colorCombPbPb6080-5);
   graphPCMPi0CorrectedSpecSys6080->Draw("E2same");
  
   histoCorrectedYieldPi02760GeV->SetBinContent(19,0);
   DrawGammaSetMarker(histoCorrectedYieldPi02760GeV, markerStylePP,markerSizePP, kBlack , kBlack);
   histoCorrectedYieldPi02760GeV->Draw("p,same,e1");
   histoPCMPi0CorrectedSpec0005->SetBinContent(2,0);
   histoPCMPi0CorrectedSpec0005->SetBinContent(19,0);
   histoPCMPi0CorrectedSpec0005->SetBinContent(20,0);
   histoPCMPi0CorrectedSpec0005->Scale(4.);
   DrawGammaSetMarker(histoPCMPi0CorrectedSpec0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
   histoPCMPi0CorrectedSpec0005->Draw("p,same,e1");
   histoPCMPi0CorrectedSpec0510->SetBinContent(2,0);
   histoPCMPi0CorrectedSpec0510->SetBinContent(19,0);
   histoPCMPi0CorrectedSpec0510->SetBinContent(20,0);
   histoPCMPi0CorrectedSpec0510->Scale(2.);
   DrawGammaSetMarker(histoPCMPi0CorrectedSpec0510, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0510 , colorCombPbPb0510);
   histoPCMPi0CorrectedSpec0510->Draw("p,same,e1");
   
   histoPCMPi0CorrectedSpec1020->SetBinContent(19,0);
   histoPCMPi0CorrectedSpec1020->SetBinContent(20,0);
   DrawGammaSetMarker(histoPCMPi0CorrectedSpec1020, markerStylePbPb1020,markerSizePbPb1020, colorCombPbPb1020 , colorCombPbPb1020);
   histoPCMPi0CorrectedSpec1020->Draw("p,same,e1");
   
   histoPCMPi0CorrectedSpec2040->SetBinContent(19,0);
   histoPCMPi0CorrectedSpec2040->SetBinContent(20,0);
   DrawGammaSetMarker(histoPCMPi0CorrectedSpec2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
   histoPCMPi0CorrectedSpec2040->Draw("p,same,e1");
   
   histoPCMPi0CorrectedSpec4060->SetBinContent(19,0);
   histoPCMPi0CorrectedSpec4060->SetBinContent(20,0);
   DrawGammaSetMarker(histoPCMPi0CorrectedSpec4060, markerStylePbPb4060,markerSizePbPb4060, colorCombPbPb4060 , colorCombPbPb4060);
   histoPCMPi0CorrectedSpec4060->Draw("p,same,e1");
   
   histoPCMPi0CorrectedSpec6080->SetBinContent(19,0);
   histoPCMPi0CorrectedSpec6080->SetBinContent(20,0);
   DrawGammaSetMarker(histoPCMPi0CorrectedSpec6080, markerStylePbPb6080,markerSizePbPb6080, colorCombPbPb6080 , colorCombPbPb6080);
   histoPCMPi0CorrectedSpec6080->Draw("same");

   TLatex *labelSpectraPi0LabelPbPb = new TLatex(0.6,0.89,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
   SetStyleTLatex( labelSpectraPi0LabelPbPb, 0.05  ,4);
   labelSpectraPi0LabelPbPb->Draw();
   
   TLegend* legendSpectraPi0 = new TLegend(0.16,0.09,0.73,0.3);
   legendSpectraPi0->SetFillColor(0);
   legendSpectraPi0->SetLineColor(0);
   legendSpectraPi0->SetTextSize(0.025);
   legendSpectraPi0->SetNColumns(2);
   legendSpectraPi0->SetMargin(0.2);
   legendSpectraPi0->AddEntry(graphPCMPi0CorrectedSpecSys0005,"PCM #times 4","pf");
   legendSpectraPi0->AddEntry((TObject*)0, collisionSystemPbPb0005.Data(),"");
   legendSpectraPi0->AddEntry(graphPCMPi0CorrectedSpecSys0510,"PCM #times 2","pf");
   legendSpectraPi0->AddEntry((TObject*)0, collisionSystemPbPb0510.Data(),"");
   legendSpectraPi0->AddEntry(graphPCMPi0CorrectedSpecSys1020,"PCM","pf");
   legendSpectraPi0->AddEntry((TObject*)0, collisionSystemPbPb1020.Data(),"");
   legendSpectraPi0->AddEntry(graphPCMPi0CorrectedSpecSys2040,"PCM","pf");
   legendSpectraPi0->AddEntry((TObject*)0, collisionSystemPbPb2040.Data(),"");
   legendSpectraPi0->AddEntry(graphPCMPi0CorrectedSpecSys4060,"PCM","pf");
   legendSpectraPi0->AddEntry((TObject*)0, collisionSystemPbPb4060.Data(),"");
   legendSpectraPi0->AddEntry(graphPCMPi0CorrectedSpecSys6080,"PCM","pf");
   legendSpectraPi0->AddEntry((TObject*)0, collisionSystemPbPb6080.Data(),"");
   legendSpectraPi0->AddEntry(graphCorrectedYieldSysPi02760GeV,"PCM","pf");
   legendSpectraPi0->AddEntry((TObject*)0, collisionSystemPP2760GeV.Data(),"");
   legendSpectraPi0->Draw();
   
   canvasPi0SpectraAllTogether->Update();
   canvasPi0SpectraAllTogether->Print(Form("%s/Pi0_Spectra_All.%s",outputDir.Data(),suffix.Data()));

   
   TCanvas* canvasSpectraMCAllTogether = new TCanvas("canvasSpectraMCAllTogether","",200,10,1200,1100);  // gives the page size
   DrawGammaCanvasSettings( canvasSpectraMCAllTogether,  0.13, 0.01, 0.015, 0.08);
   
   canvasSpectraMCAllTogether->SetLogy();
   canvasSpectraMCAllTogether->SetLogx();
   TH2F * histo2DSpectraMCAll;
   histo2DSpectraMCAll = new TH2F("histo2DSpectraMCAll","histo2DSpectraMCAll",1000,0.,30.,1000,1e-8,2e4 );
   SetStyleHistoTH2ForGraphs(histo2DSpectraMCAll, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}
   histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.3,30);
   histo2DSpectraMCAll->GetXaxis()->SetLabelOffset(-0.01);
   histo2DSpectraMCAll->GetYaxis()->SetLabelOffset(0.01);
   histo2DSpectraMCAll->DrawCopy(); 
   
   DrawGammaSetMarker(histoMCYieldPi02760GeVFinal, markerStylePP,markerSizePP, kBlack , kBlack);//, colorCombPbPb1020-5);
   histoMCYieldPi02760GeVFinal->Draw("hist,l,same");
   DrawGammaSetMarker(histoMCYieldPi02760GeVAddSigWOWeighting, markerStylePP,markerSizePP, kGray, kGray);//, colorCombPbPb1020-5);
   histoMCYieldPi02760GeVAddSigWOWeighting->SetLineStyle(2);
   histoMCYieldPi02760GeVAddSigWOWeighting->Draw("hist,l,same");
   TH1D* histoMCYieldPi02760GeVAddSigWOWeighting0005 = (TH1D*)histoMCYieldPi02760GeVAddSigWOWeighting->Clone("histoMCYieldPi02760GeVAddSigWOWeighting0005");
   histoMCYieldPi02760GeVAddSigWOWeighting0005->Scale(6.);
   DrawGammaSetMarker(histoMCYieldPi02760GeVAddSigWOWeighting0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005-4 , colorCombPbPb0005-4);//, colorCombPbPb0510-5);
   histoMCYieldPi02760GeVAddSigWOWeighting0005->SetLineStyle(2);
   histoMCYieldPi02760GeVAddSigWOWeighting0005->Draw("hist,l,same");
   
   TH1D* histoMCYieldPi02760GeVAddSigWOWeighting6080 = (TH1D*)histoMCYieldPi02760GeVAddSigWOWeighting->Clone("histoMCYieldPi02760GeVAddSigWOWeighting6080");
   histoMCYieldPi02760GeVAddSigWOWeighting6080->Scale(3.);
   DrawGammaSetMarker(histoMCYieldPi02760GeVAddSigWOWeighting6080, markerStylePbPb6080,markerSizePbPb6080, colorCombPbPb6080-2 , colorCombPbPb6080-2);//, colorCombPbPb6080-5);
   histoMCYieldPi02760GeVAddSigWOWeighting6080->SetLineStyle(2);
   histoMCYieldPi02760GeVAddSigWOWeighting6080->Draw("hist,l,same");
   
   DrawGammaSetMarker(histoMCYieldPt0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);//, colorCombPbPb0510-5);
   histoMCYieldPt0005->Draw("hist,l,same");
   DrawGammaSetMarker(histoMCYieldPt0510, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0510 , colorCombPbPb0510);//, colorCombPbPb0510-5);
   histoMCYieldPt0510->Draw("hist,l,same");
   
   DrawGammaSetMarker(histoMCYieldPt1020, markerStylePbPb1020,markerSizePbPb1020, colorCombPbPb1020 , colorCombPbPb1020);//, colorCombPbPb1020-5);
   histoMCYieldPt1020->Draw("hist,l,same");
   DrawGammaSetMarker(histoMCYieldPt2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);//, colorCombPbPb2040-5);
   histoMCYieldPt2040->Draw("hist,l,same");
   DrawGammaSetMarker(histoMCYieldPt4060, markerStylePbPb4060,markerSizePbPb4060, colorCombPbPb4060 , colorCombPbPb4060);//, colorCombPbPb4060-5);
   histoMCYieldPt4060->Draw("hist,l,same");
   DrawGammaSetMarker(histoMCYieldPt6080, markerStylePbPb6080,markerSizePbPb6080, colorCombPbPb6080 , colorCombPbPb6080);//, colorCombPbPb6080-5);
   histoMCYieldPt6080->Draw("hist,l,same");
   
   TLatex *labelSpectraPi0InputPbPb = new TLatex(0.5,0.89,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-} MC input");
   SetStyleTLatex( labelSpectraPi0InputPbPb, 0.04  ,4);
   labelSpectraPi0InputPbPb->Draw();
   
   TLegend* legendSpectraMC = new TLegend(0.16,0.09,0.5,0.35);
   legendSpectraMC->SetFillColor(0);
   legendSpectraMC->SetLineColor(0);
   legendSpectraMC->SetTextSize(0.025);
   legendSpectraMC->SetNColumns(1);
   legendSpectraMC->SetMargin(0.2);
   
   legendSpectraMC->AddEntry(histoMCYieldPt0005,Form("%s",collisionSystemPbPb0005.Data()),"l");
   legendSpectraMC->AddEntry(histoMCYieldPi02760GeVAddSigWOWeighting0005,Form("int. add. sig. %s",collisionSystemPbPb0005.Data()),"l");
   legendSpectraMC->AddEntry(histoMCYieldPt0510,Form("%s",collisionSystemPbPb0510.Data()),"l");
   legendSpectraMC->AddEntry(histoMCYieldPt1020,Form("%s",collisionSystemPbPb1020.Data()),"l");
   legendSpectraMC->AddEntry(histoMCYieldPt2040,Form("%s",collisionSystemPbPb2040.Data()),"l");
   legendSpectraMC->AddEntry(histoMCYieldPt4060,Form("%s",collisionSystemPbPb4060.Data()),"l");
   legendSpectraMC->AddEntry(histoMCYieldPt6080,Form("%s",collisionSystemPbPb6080.Data()),"l");
   legendSpectraMC->AddEntry(histoMCYieldPi02760GeVAddSigWOWeighting6080,Form("int. add. sig. %s",collisionSystemPbPb6080.Data()),"l");
   legendSpectraMC->AddEntry(histoMCYieldPi02760GeVFinal,Form("%s",collisionSystemPP2760GeV.Data()),"l");
   legendSpectraMC->AddEntry(histoMCYieldPi02760GeVAddSigWOWeighting,Form("add. sig. %s",collisionSystemPP2760GeV.Data()),"l");
   legendSpectraMC->Draw();
   
   canvasSpectraMCAllTogether->Update();
   canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectra_All.%s",outputDir.Data(),suffix.Data()));

   
   
   canvasPi0SpectraAllTogether->cd();
   canvasPi0SpectraAllTogether->SetLogy();
   canvasPi0SpectraAllTogether->SetLogx();
   TH2F * histo2DPi0SpectraAllPlusMC;
   histo2DPi0SpectraAllPlusMC = new TH2F("histo2DPi0SpectraAllPlusMC","histo2DPi0SpectraAllPlusMC",1000,0.23,15.,1000,1e-9,2e4 );
   SetStyleHistoTH2ForGraphs(histo2DPi0SpectraAllPlusMC, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}
   histo2DPi0SpectraAllPlusMC->GetXaxis()->SetLabelOffset(-0.01);
   histo2DPi0SpectraAllPlusMC->GetYaxis()->SetLabelOffset(0.01);
   histo2DPi0SpectraAllPlusMC->DrawCopy(); 
   
   
   graphCorrectedYieldSysPi02760GeV->Draw("E2same");
   graphPCMPi0CorrectedSpecSys0005->Draw("E2same");
   graphPCMPi0CorrectedSpecSys0510->Draw("E2same");
   graphPCMPi0CorrectedSpecSys1020->Draw("E2same");
   graphPCMPi0CorrectedSpecSys2040->Draw("E2same");
   graphPCMPi0CorrectedSpecSys4060->Draw("E2same");
   graphPCMPi0CorrectedSpecSys6080->Draw("E2same");
   histoCorrectedYieldPi02760GeV->Draw("p,same,e1");
   
   histoPCMPi0CorrectedSpec0005->Draw("p,same,e1");
   histoPCMPi0CorrectedSpec0510->Draw("p,same,e1");
   histoPCMPi0CorrectedSpec1020->Draw("p,same,e1");
   histoPCMPi0CorrectedSpec2040->Draw("p,same,e1");
   histoPCMPi0CorrectedSpec4060->Draw("p,same,e1");
   histoPCMPi0CorrectedSpec6080->Draw("same");

   histoMCYieldPi02760GeVFinal->Draw("hist,l,same");
   histoMCYieldPt0005->Scale(4.);
   DrawGammaSetMarker(histoMCYieldPt0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);//, colorCombPbPb0510-5);
   histoMCYieldPt0005->Draw("hist,l,same");
   histoMCYieldPt0510->Scale(2.);
   DrawGammaSetMarker(histoMCYieldPt0510, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0510 , colorCombPbPb0510);//, colorCombPbPb0510-5);
   histoMCYieldPt0510->Draw("hist,l,same");
   histoMCYieldPt1020->Draw("hist,l,same");
   histoMCYieldPt2040->Draw("hist,l,same");
   histoMCYieldPt4060->Draw("hist,l,same");
   histoMCYieldPt6080->Draw("hist,l,same");
   
   
   labelSpectraPi0LabelPbPb->Draw();
   
   TLegend* legendSpectraPi0PlusMC = new TLegend(0.16,0.09,0.83,0.3);
   legendSpectraPi0PlusMC->SetFillColor(0);
   legendSpectraPi0PlusMC->SetLineColor(0);
   legendSpectraPi0PlusMC->SetTextSize(0.025);
   legendSpectraPi0PlusMC->SetNColumns(3);
   legendSpectraPi0PlusMC->SetMargin(0.2);
   legendSpectraPi0PlusMC->AddEntry(graphPCMPi0CorrectedSpecSys0005,"PCM","pf");
   legendSpectraPi0PlusMC->AddEntry(histoMCYieldPt0005,"HIJING #times 4","l");
   legendSpectraPi0PlusMC->AddEntry((TObject*)0, collisionSystemPbPb0005.Data(),"");
   legendSpectraPi0PlusMC->AddEntry(graphPCMPi0CorrectedSpecSys0510,"PCM","pf");
   legendSpectraPi0PlusMC->AddEntry(histoMCYieldPt0510,"HIJING #times 2","l");
   legendSpectraPi0PlusMC->AddEntry((TObject*)0, collisionSystemPbPb0510.Data(),"");
   legendSpectraPi0PlusMC->AddEntry(graphPCMPi0CorrectedSpecSys1020,"PCM","pf");
   legendSpectraPi0PlusMC->AddEntry(histoMCYieldPt1020,"HIJING","l");
   legendSpectraPi0PlusMC->AddEntry((TObject*)0, collisionSystemPbPb1020.Data(),"");
   legendSpectraPi0PlusMC->AddEntry(graphPCMPi0CorrectedSpecSys2040,"PCM","pf");
   legendSpectraPi0PlusMC->AddEntry(histoMCYieldPt2040,"HIJING","l");
   legendSpectraPi0PlusMC->AddEntry((TObject*)0, collisionSystemPbPb2040.Data(),"");
   legendSpectraPi0PlusMC->AddEntry(graphPCMPi0CorrectedSpecSys4060,"PCM","pf");
   legendSpectraPi0PlusMC->AddEntry(histoMCYieldPt4060,"HIJING","l");
   legendSpectraPi0PlusMC->AddEntry((TObject*)0, collisionSystemPbPb4060.Data(),"");
   legendSpectraPi0PlusMC->AddEntry(graphPCMPi0CorrectedSpecSys6080,"PCM","pf");
   legendSpectraPi0PlusMC->AddEntry(histoMCYieldPt6080,"HIJING","l");
   legendSpectraPi0PlusMC->AddEntry((TObject*)0, collisionSystemPbPb6080.Data(),"");
   legendSpectraPi0PlusMC->AddEntry(graphCorrectedYieldSysPi02760GeV,"PCM","pf");
   legendSpectraPi0PlusMC->AddEntry(histoMCYieldPi02760GeVFinal,"PHO+PYT8","l");
   legendSpectraPi0PlusMC->AddEntry((TObject*)0, collisionSystemPP2760GeV.Data(),"");
   legendSpectraPi0PlusMC->Draw();
   
   canvasPi0SpectraAllTogether->Update();
   canvasPi0SpectraAllTogether->Print(Form("%s/Pi0_SpectraPlusMC_All.%s",outputDir.Data(),suffix.Data()));

   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 0-5 % *******************
   //**********************************************************************************
   
   canvasSpectraMCAllTogether->cd();
   canvasSpectraMCAllTogether->SetLogx(1);
   histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
   histo2DSpectraMCAll->DrawCopy(); 
   
   histoMCYieldPt0005->Scale(1./4.);
   histoMCYieldPt0005->SetMarkerStyle(markerStylePbPb0005MC);
   histoMCYieldPt0005->Draw("hist,pe1,same");
 
   TString  nameFinalResDat = Form("%s/FitResultsMC.dat",outputDir.Data());
   TString forOutput;
   fstream  fileFinalResults;
   fileFinalResults.open(nameFinalResDat.Data(), ios::out);
   histoPCMPi0CorrectedSpec0005->Scale(1/4.);
   histoPCMPi0CorrectedSpec0005->Draw("hist,pe1,same");
   TF1* fitYieldDataQCDPi0PbPb0005 = FitObject("qcd","fitYieldDataQCDPi0PbPb0005","Pi0",histoPCMPi0CorrectedSpec0005,0.6,8,NULL,"QNRME+");
   DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb0005, 1, 1.5, colorCombPbPb0005);
//    fitYieldDataQCDPi0PbPb0005->SetRange(0.1,20);
   fitYieldDataQCDPi0PbPb0005->Draw("same");
   forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb0005);
   fileFinalResults<< forOutput.Data()<< endl;  
   canvasSpectraMCAllTogether->Update();
   canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted0005.%s",outputDir.Data(),suffix.Data()));
   
   TH1D* histoRatioDatatoFitQCDPbPb0005 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec0005, fitYieldDataQCDPi0PbPb0005);
   TH1D* histoRatioMCtoDataFitQCDPbPb0005 = CalculateHistoRatioToFit (histoMCYieldPt0005, fitYieldDataQCDPi0PbPb0005);
   TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb0005 = NULL;
   if (histoMCYieldPtWOWeights0005) histoRatioMCUnweightedtoDataFitQCDPbPb0005 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights0005, fitYieldDataQCDPi0PbPb0005);
   canvasFraction2->cd();
   if (histoRatioMCUnweightedtoDataFitQCDPbPb0005) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb0005, markerStylePbPb0510,markerSizePbPb0005, colorCombPbPb0510 , colorCombPbPb0510 );
   DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
   DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb0005, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
   DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb0005,
               "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 0-5%",
               kFALSE, 1.5, 0, kTRUE,
               kTRUE, -0.5, 8.,
               kTRUE, 0., 7.9);
   histoRatioDatatoFitQCDPbPb0005->Draw("same,e,p");  
   if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb0005->Draw("same,e,p");  
   if (histoRatioMCUnweightedtoDataFitQCDPbPb0005) histoRatioMCUnweightedtoDataFitQCDPbPb0005->Draw("same,e,p");  
   TLegend* legendFit0005 = new TLegend(0.16,0.81,0.4,0.9);
   legendFit0005->SetFillColor(0);
   legendFit0005->SetLineColor(0);
   legendFit0005->SetTextSize(0.025);
//    legendFit0005->SetNColumns(3);
   legendFit0005->SetMargin(0.2);
   legendFit0005->AddEntry(histoRatioDatatoFitQCDPbPb0005,"Data/QCD fit to Data (0.4 <pT<8)","p");
   if (runDrawReweighted) legendFit0005->AddEntry(histoRatioMCtoDataFitQCDPbPb0005,"MC weighted/QCD fit to Data (0.4 <pT<8)","p");
   if (histoRatioMCUnweightedtoDataFitQCDPbPb0005) legendFit0005->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb0005,"MC/QCD fit to Data (0.4 <pT<8)","p");
   legendFit0005->Draw();

   TLatex *labelRatioMCData0005 = new TLatex(0.05,0.9,collisionSystemPbPb0005.Data());
   SetStyleTLatex( labelRatioMCData0005, 0.062,4);
   labelRatioMCData0005->Draw();

   
   DrawGammaLines(0., 30.,1., 1.,0.1);
   canvasFraction2->Update();
   canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit0005.%s",outputDir.Data(),suffix.Data()));

   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 0-10 % ******************
   //**********************************************************************************


   canvasSpectraMCAllTogether->cd();
   canvasSpectraMCAllTogether->SetLogx(1);
   histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
   histo2DSpectraMCAll->DrawCopy(); 
  
   DrawGammaSetMarker(histoMCYieldPt0010, markerStylePbPb0010MC,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
   histoMCYieldPt0010->Draw("hist,pe1,same");
 
   DrawGammaSetMarker(histoPCMPi0CorrectedSpec0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
   histoPCMPi0CorrectedSpec0010->Draw("hist,pe1,same");
   TF1* fitYieldDataQCDPi0PbPb0010 = FitObject("qcd","fitYieldDataQCDPi0PbPb0010","Pi0",histoPCMPi0CorrectedSpec0010,0.6,8,NULL,"QNRME+");
   DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb0010, 1, 1.5, colorCombPbPb0010);
//    fitYieldDataQCDPi0PbPb0010->SetRange(0.1,20);
   fitYieldDataQCDPi0PbPb0010->Draw("same");
   forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb0010);
   fileFinalResults<< forOutput.Data()<< endl;  
   canvasSpectraMCAllTogether->Update();
   canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted0010.%s",outputDir.Data(),suffix.Data()));
   
   TH1D* histoRatioDatatoFitQCDPbPb0010 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec0010, fitYieldDataQCDPi0PbPb0010);
   TH1D* histoRatioMCtoDataFitQCDPbPb0010 = CalculateHistoRatioToFit (histoMCYieldPt0010, fitYieldDataQCDPi0PbPb0010);
   TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb0010 = NULL;
   if (histoMCYieldPtWOWeights0010) histoRatioMCUnweightedtoDataFitQCDPbPb0010 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights0010, fitYieldDataQCDPi0PbPb0010);
   canvasFraction2->cd();
   if (histoRatioMCUnweightedtoDataFitQCDPbPb0010) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb0010, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
   DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb0010, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
   DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb0010, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
   DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb0010,
               "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 0-10%",
               kFALSE, 1.5, 0, kTRUE,
               kTRUE, -0.5, 8.,
               kTRUE, 0., 7.9);
   histoRatioDatatoFitQCDPbPb0010->Draw("same,e,p");  
   if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb0010->Draw("same,e,p");  
   if (histoRatioMCUnweightedtoDataFitQCDPbPb0010) histoRatioMCUnweightedtoDataFitQCDPbPb0010->Draw("same,e,p");  
   TLegend* legendFit0010 = new TLegend(0.16,0.81,0.4,0.9);
   legendFit0010->SetFillColor(0);
   legendFit0010->SetLineColor(0);
   legendFit0010->SetTextSize(0.025);
//    legendFit0010->SetNColumns(3);
   legendFit0010->SetMargin(0.2);
   legendFit0010->AddEntry(histoRatioDatatoFitQCDPbPb0010,"Data/QCD fit to Data (0.6 <pT<8)","p");
   if (runDrawReweighted) legendFit0010->AddEntry(histoRatioMCtoDataFitQCDPbPb0010,"MC weighted/QCD fit to Data (0.6 <pT<8)","p");
   if (histoRatioMCUnweightedtoDataFitQCDPbPb0010) legendFit0010->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb0010,"MC/QCD fit to Data (0.6 <pT<8)","p");
   legendFit0010->Draw();
   
   DrawGammaLines(0., 30.,1., 1.,0.1);
   canvasFraction2->Update();
   canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit0010.%s",outputDir.Data(),suffix.Data()));

   
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 5-10 % ******************
   //**********************************************************************************

   canvasSpectraMCAllTogether->cd();
   canvasSpectraMCAllTogether->SetLogx(1);
   histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
   histo2DSpectraMCAll->DrawCopy(); 
  
   histoMCYieldPt0510->Scale(1./2.);
   histoMCYieldPt0510->SetMarkerStyle(markerStylePbPb0510MC);
   histoMCYieldPt0510->Draw("hist,pe1,same");
   histoPCMPi0CorrectedSpec0510->Scale(1/2.);
   histoPCMPi0CorrectedSpec0510->Draw("hist,pe1,same");
   TF1* fitYieldDataQCDPi0PbPb0510 = FitObject("qcd","fitYieldDataQCDPi0PbPb0510","Pi0",histoPCMPi0CorrectedSpec0510,0.6,8,NULL,"QNRME+");
   DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb0510, 1, 1.5, colorCombPbPb0510);
//    fitYieldDataQCDPi0PbPb0510->SetRange(0.1,20);
   fitYieldDataQCDPi0PbPb0510->Draw("same");
   forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb0510);
   fileFinalResults<< forOutput.Data()<< endl;  
   canvasSpectraMCAllTogether->Update();
   canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted0510.%s",outputDir.Data(),suffix.Data()));
   
   TH1D* histoRatioDatatoFitQCDPbPb0510 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec0510, fitYieldDataQCDPi0PbPb0510);
   TH1D* histoRatioMCtoDataFitQCDPbPb0510 = CalculateHistoRatioToFit (histoMCYieldPt0510, fitYieldDataQCDPi0PbPb0510);
   TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb0510 = NULL;
   if (histoMCYieldPtWOWeights0510) histoRatioMCUnweightedtoDataFitQCDPbPb0510 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights0510, fitYieldDataQCDPi0PbPb0510);
   canvasFraction2->cd();
   if (histoRatioMCUnweightedtoDataFitQCDPbPb0510) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
   DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb0510, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
   DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb0510, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
   DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb0510,
               "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 5-10%",
               kFALSE, 1.5, 0, kTRUE,
               kTRUE, -0.5, 8.,
               kTRUE, 0., 7.9);
   histoRatioDatatoFitQCDPbPb0510->Draw("same,e,p");  
   if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb0510->Draw("same,e,p");  
   if (histoRatioMCUnweightedtoDataFitQCDPbPb0510) histoRatioMCUnweightedtoDataFitQCDPbPb0510->Draw("same,e,p");  
   TLegend* legendFit0510 = new TLegend(0.16,0.81,0.4,0.9);
   legendFit0510->SetFillColor(0);
   legendFit0510->SetLineColor(0);
   legendFit0510->SetTextSize(0.025);
//    legendFit0510->SetNColumns(3);
   legendFit0510->SetMargin(0.2);
   legendFit0510->AddEntry(histoRatioDatatoFitQCDPbPb0510,"Data/QCD fit to Data (0.6 <pT<8)","p");
   if (runDrawReweighted) legendFit0510->AddEntry(histoRatioMCtoDataFitQCDPbPb0510,"MC weighted/QCD fit to Data (0.6 <pT<8)","p");
   if (histoRatioMCUnweightedtoDataFitQCDPbPb0510) legendFit0510->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb0510,"MC/QCD fit to Data (0.6 <pT<8)","p");
   legendFit0510->Draw();
   
   DrawGammaLines(0., 30.,1., 1.,0.1);
   canvasFraction2->Update();
   canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit0510.%s",outputDir.Data(),suffix.Data()));

   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 10-20 % *****************
   //**********************************************************************************
   
   canvasSpectraMCAllTogether->cd();
   canvasSpectraMCAllTogether->SetLogx(1);
   histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
   histo2DSpectraMCAll->DrawCopy(); 
  
   histoMCYieldPt1020->SetMarkerStyle(markerStylePbPb1020MC);
   histoMCYieldPt1020->Draw("hist,pe1,same");
 
   histoPCMPi0CorrectedSpec1020->Draw("hist,pe1,same");
   TF1* fitYieldDataQCDPi0PbPb1020 = FitObject("qcd","fitYieldDataQCDPi0PbPb1020","Pi0",histoPCMPi0CorrectedSpec1020,0.4,8,NULL,"QNRME+");
   DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb1020, 1, 1.5, colorCombPbPb1020);
//    fitYieldDataQCDPi0PbPb1020->SetRange(0.1,20);
   fitYieldDataQCDPi0PbPb1020->Draw("same");
   forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb1020);
   fileFinalResults<< forOutput.Data()<< endl;  
   canvasSpectraMCAllTogether->Update();
   canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted1020.%s",outputDir.Data(),suffix.Data()));
   
   TH1D* histoRatioDatatoFitQCDPbPb1020 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec1020, fitYieldDataQCDPi0PbPb1020);
   TH1D* histoRatioMCtoDataFitQCDPbPb1020 = CalculateHistoRatioToFit (histoMCYieldPt1020, fitYieldDataQCDPi0PbPb1020);
   TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb1020 = NULL;
   if (histoMCYieldPtWOWeights1020) histoRatioMCUnweightedtoDataFitQCDPbPb1020 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights1020, fitYieldDataQCDPi0PbPb1020);
   canvasFraction2->cd();
   if (histoRatioMCUnweightedtoDataFitQCDPbPb1020) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb1020, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
   DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb1020, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
   DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb1020, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
   DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb1020,
               "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 10-20%",
               kFALSE, 1.5, 0, kTRUE,
               kTRUE, -0.5, 8.,
               kTRUE, 0., 7.9);
   histoRatioDatatoFitQCDPbPb1020->Draw("same,e,p");  
   if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb1020->Draw("same,e,p");  
   if (histoRatioMCUnweightedtoDataFitQCDPbPb1020) histoRatioMCUnweightedtoDataFitQCDPbPb1020->Draw("same,e,p");  
   TLegend* legendFit1020 = new TLegend(0.16,0.81,0.4,0.9);
   legendFit1020->SetFillColor(0);
   legendFit1020->SetLineColor(0);
   legendFit1020->SetTextSize(0.025);
//    legendFit1020->SetNColumns(3);
   legendFit1020->SetMargin(0.2);
   legendFit1020->AddEntry(histoRatioDatatoFitQCDPbPb1020,"Data/QCD fit to Data (0.4 <pT<8)","p");
   if (runDrawReweighted) legendFit1020->AddEntry(histoRatioMCtoDataFitQCDPbPb1020,"MC weighted/QCD fit to Data (0.4 <pT<8)","p");
   if (histoRatioMCUnweightedtoDataFitQCDPbPb1020) legendFit1020->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb1020,"MC/QCD fit to Data (0.4 <pT<8)","p");
   legendFit1020->Draw();
   
   DrawGammaLines(0., 30.,1., 1.,0.1);
   canvasFraction2->Update();
   canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit1020.%s",outputDir.Data(),suffix.Data()));
 
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 20-40 % *****************
   //**********************************************************************************
   
   canvasSpectraMCAllTogether->cd();
   canvasSpectraMCAllTogether->SetLogx(1);
   histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
   histo2DSpectraMCAll->DrawCopy(); 
  
   histoMCYieldPt2040->SetMarkerStyle(markerStylePbPb2040MC);
   histoMCYieldPt2040->Draw("hist,pe1,same");
 
   histoPCMPi0CorrectedSpec2040->Draw("hist,pe1,same");
   TF1* fitYieldDataQCDPi0PbPb2040 = FitObject("qcd","fitYieldDataQCDPi0PbPb2040","Pi0",histoPCMPi0CorrectedSpec2040,0.4,8,NULL,"QNRME+");
   DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb2040, 1, 1.5, colorCombPbPb2040);
//    fitYieldDataQCDPi0PbPb2040->SetRange(0.1,20);
   fitYieldDataQCDPi0PbPb2040->Draw("same");
   forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb2040);
   fileFinalResults<< forOutput.Data()<< endl;  
   canvasSpectraMCAllTogether->Update();
   canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted2040.%s",outputDir.Data(),suffix.Data()));
   
   TH1D* histoRatioDatatoFitQCDPbPb2040 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec2040, fitYieldDataQCDPi0PbPb2040);
   TH1D* histoRatioMCtoDataFitQCDPbPb2040 = CalculateHistoRatioToFit (histoMCYieldPt2040, fitYieldDataQCDPi0PbPb2040);
   TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb2040 = NULL;
   if (histoMCYieldPtWOWeights2040) histoRatioMCUnweightedtoDataFitQCDPbPb2040 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights2040, fitYieldDataQCDPi0PbPb2040);
   canvasFraction2->cd();
   if (histoRatioMCUnweightedtoDataFitQCDPbPb2040) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb2040, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
   DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb2040, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
   DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb2040, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
   DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb2040,
               "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 20-40%",
               kFALSE, 1.5, 0, kTRUE,
               kTRUE, -0.5, 8.,
               kTRUE, 0., 7.9);
   histoRatioDatatoFitQCDPbPb2040->Draw("same,e,p");  
   if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb2040->Draw("same,e,p");  
   if (histoRatioMCUnweightedtoDataFitQCDPbPb2040) histoRatioMCUnweightedtoDataFitQCDPbPb2040->Draw("same,e,p");  
   TLegend* legendFit2040 = new TLegend(0.16,0.81,0.4,0.9);
   legendFit2040->SetFillColor(0);
   legendFit2040->SetLineColor(0);
   legendFit2040->SetTextSize(0.025);
//    legendFit2040->SetNColumns(3);
   legendFit2040->SetMargin(0.2);
   legendFit2040->AddEntry(histoRatioDatatoFitQCDPbPb2040,"Data/QCD fit to Data (0.4 <pT<8)","p");
   if (runDrawReweighted) legendFit2040->AddEntry(histoRatioMCtoDataFitQCDPbPb2040,"MC weighted/QCD fit to Data (0.4 <pT<8)","p");
   if (histoRatioMCUnweightedtoDataFitQCDPbPb2040) legendFit2040->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb2040,"MC/QCD fit to Data (0.4 <pT<8)","p");
   legendFit2040->Draw();
   
   DrawGammaLines(0., 30.,1., 1.,0.1);
   canvasFraction2->Update();
   canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit2040.%s",outputDir.Data(),suffix.Data()));
 
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 40-60 % *****************
   //**********************************************************************************
   
   canvasSpectraMCAllTogether->cd();
   canvasSpectraMCAllTogether->SetLogx(1);
   histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
   histo2DSpectraMCAll->DrawCopy(); 
  
   histoMCYieldPt4060->SetMarkerStyle(markerStylePbPb4060MC);
   histoMCYieldPt4060->Draw("hist,pe1,same");
 
   histoPCMPi0CorrectedSpec4060->Draw("hist,pe1,same");
   TF1* fitYieldDataQCDPi0PbPb4060 = FitObject("qcd","fitYieldDataQCDPi0PbPb4060","Pi0",histoPCMPi0CorrectedSpec4060,0.4,8,NULL,"QNRME+");
   DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb4060, 1, 1.5, colorCombPbPb4060);
//    fitYieldDataQCDPi0PbPb4060->SetRange(0.1,20);
   fitYieldDataQCDPi0PbPb4060->Draw("same");
   forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb4060);
   fileFinalResults<< forOutput.Data()<< endl;  
   canvasSpectraMCAllTogether->Update();
   canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted4060.%s",outputDir.Data(),suffix.Data()));
   
   TH1D* histoRatioDatatoFitQCDPbPb4060 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec4060, fitYieldDataQCDPi0PbPb4060);
   TH1D* histoRatioMCtoDataFitQCDPbPb4060 = CalculateHistoRatioToFit (histoMCYieldPt4060, fitYieldDataQCDPi0PbPb4060);
   TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb4060 = NULL;
   if (histoMCYieldPtWOWeights4060) histoRatioMCUnweightedtoDataFitQCDPbPb4060 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights4060, fitYieldDataQCDPi0PbPb4060);
   canvasFraction2->cd();
   if (histoRatioMCUnweightedtoDataFitQCDPbPb4060) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb4060, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
   DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb4060, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
   DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb4060, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
   DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb4060,
               "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 40-60%",
               kFALSE, 1.5, 0, kTRUE,
               kTRUE, -0.5, 8.,
               kTRUE, 0., 7.9);
   histoRatioDatatoFitQCDPbPb4060->Draw("same,e,p");  
   if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb4060->Draw("same,e,p");  
   if (histoRatioMCUnweightedtoDataFitQCDPbPb4060) histoRatioMCUnweightedtoDataFitQCDPbPb4060->Draw("same,e,p");  
   TLegend* legendFit4060 = new TLegend(0.16,0.81,0.4,0.9);
   legendFit4060->SetFillColor(0);
   legendFit4060->SetLineColor(0);
   legendFit4060->SetTextSize(0.025);
//    legendFit4060->SetNColumns(3);
   legendFit4060->SetMargin(0.2);
   legendFit4060->AddEntry(histoRatioDatatoFitQCDPbPb4060,"Data/QCD fit to Data (0.4 <pT<8)","p");
   if (runDrawReweighted) legendFit4060->AddEntry(histoRatioMCtoDataFitQCDPbPb4060,"MC weighted/QCD fit to Data (0.4 <pT<8)","p");
   if (histoRatioMCUnweightedtoDataFitQCDPbPb4060) legendFit4060->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb4060,"MC/QCD fit to Data (0.4 <pT<8)","p");
   legendFit4060->Draw();
   
   DrawGammaLines(0., 30.,1., 1.,0.1);
   canvasFraction2->Update();
   canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit4060.%s",outputDir.Data(),suffix.Data()));
 
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 60-80 % *****************
   //**********************************************************************************
   
   canvasSpectraMCAllTogether->cd();
   canvasSpectraMCAllTogether->SetLogx(1);
   histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
   histo2DSpectraMCAll->DrawCopy(); 
  
   histoMCYieldPt6080->SetMarkerStyle(markerStylePbPb6080MC);
   histoMCYieldPt6080->Draw("hist,pe1,same");
 
   histoPCMPi0CorrectedSpec6080->Draw("hist,pe1,same");
   TF1* fitYieldDataQCDPi0PbPb6080 = FitObject("qcd","fitYieldDataQCDPi0PbPb6080","Pi0",histoPCMPi0CorrectedSpec6080,0.4,8,NULL,"QNRME+");
   DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb6080, 1, 1.5, colorCombPbPb6080);
//    fitYieldDataQCDPi0PbPb6080->SetRange(0.1,20);
   fitYieldDataQCDPi0PbPb6080->Draw("same");
   forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb6080);
   fileFinalResults<< forOutput.Data()<< endl;  
   canvasSpectraMCAllTogether->Update();
   canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted6080.%s",outputDir.Data(),suffix.Data()));
   
   TH1D* histoRatioDatatoFitQCDPbPb6080 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec6080, fitYieldDataQCDPi0PbPb6080);
   TH1D* histoRatioMCtoDataFitQCDPbPb6080 = CalculateHistoRatioToFit (histoMCYieldPt6080, fitYieldDataQCDPi0PbPb6080);
   TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb6080 = NULL;
   if (histoMCYieldPtWOWeights6080) histoRatioMCUnweightedtoDataFitQCDPbPb6080 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights6080, fitYieldDataQCDPi0PbPb6080);
   canvasFraction2->cd();
   if (histoRatioMCUnweightedtoDataFitQCDPbPb6080) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb6080, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
   DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb6080, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
   DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb6080, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
   DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb6080,
               "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 60-80%",
               kFALSE, 1.5, 0, kTRUE,
               kTRUE, -0.5, 8.,
               kTRUE, 0., 7.9);
   histoRatioDatatoFitQCDPbPb6080->Draw("same,e,p");  
   if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb6080->Draw("same,e,p");  
   if (histoRatioMCUnweightedtoDataFitQCDPbPb6080) histoRatioMCUnweightedtoDataFitQCDPbPb6080->Draw("same,e,p");  
   TLegend* legendFit6080 = new TLegend(0.16,0.81,0.4,0.9);
   legendFit6080->SetFillColor(0);
   legendFit6080->SetLineColor(0);
   legendFit6080->SetTextSize(0.025);
//    legendFit6080->SetNColumns(3);
   legendFit6080->SetMargin(0.2);
   legendFit6080->AddEntry(histoRatioDatatoFitQCDPbPb6080,"Data/QCD fit to Data (0.4 <pT<8)","p");
   if (runDrawReweighted) legendFit6080->AddEntry(histoRatioMCtoDataFitQCDPbPb6080,"MC weighted/QCD fit to Data (0.4 <pT<8)","p");
   if (histoRatioMCUnweightedtoDataFitQCDPbPb6080) legendFit6080->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb6080,"MC/QCD fit to Data (0.4 <pT<8)","p");
   legendFit6080->Draw();
   
   DrawGammaLines(0., 30.,1., 1.,0.1);
   canvasFraction2->Update();
   canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit6080.%s",outputDir.Data(),suffix.Data()));

   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 0-20 % *****************
   //**********************************************************************************
   if (directoryPi0PbPb0020){
      canvasSpectraMCAllTogether->cd();
      canvasSpectraMCAllTogether->SetLogx(1);
      histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
      histo2DSpectraMCAll->DrawCopy(); 
   
      histoMCYieldPt0020->SetMarkerStyle(markerStylePbPb6080MC);
      histoMCYieldPt0020->Draw("hist,pe1,same");
   
      histoPCMPi0CorrectedSpec0020->Draw("hist,pe1,same");
      TF1* fitYieldDataQCDPi0PbPb0020 = FitObject("qcd","fitYieldDataQCDPi0PbPb0020","Pi0",histoPCMPi0CorrectedSpec0020,0.4,8,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb0020, 1, 1.5, colorCombPbPb6080);
   //    fitYieldDataQCDPi0PbPb0020->SetRange(0.1,20);
      fitYieldDataQCDPi0PbPb0020->Draw("same");
      forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb0020);
      fileFinalResults<< forOutput.Data()<< endl;  
      canvasSpectraMCAllTogether->Update();
      canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted0020.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDPbPb0020 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec0020, fitYieldDataQCDPi0PbPb0020);
      TH1D* histoRatioMCtoDataFitQCDPbPb0020 = CalculateHistoRatioToFit (histoMCYieldPt0020, fitYieldDataQCDPi0PbPb0020);
      TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb0020 = NULL;
      if (histoMCYieldPtWOWeights0020) histoRatioMCUnweightedtoDataFitQCDPbPb0020 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights0020, fitYieldDataQCDPi0PbPb0020);
      canvasFraction2->cd();
      if (histoRatioMCUnweightedtoDataFitQCDPbPb0020) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb0020, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb0020, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
      DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb0020, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb0020,
                  "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 0-20%",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 7.9);
      histoRatioDatatoFitQCDPbPb0020->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb0020->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDPbPb0020) histoRatioMCUnweightedtoDataFitQCDPbPb0020->Draw("same,e,p");  
      TLegend* legendFit0020 = new TLegend(0.16,0.81,0.4,0.9);
      legendFit0020->SetFillColor(0);
      legendFit0020->SetLineColor(0);
      legendFit0020->SetTextSize(0.025);
   //    legendFit0020->SetNColumns(3);
      legendFit0020->SetMargin(0.2);
      legendFit0020->AddEntry(histoRatioDatatoFitQCDPbPb0020,"Data/QCD fit to Data (0.4 <pT<8)","p");
      if (runDrawReweighted) legendFit0020->AddEntry(histoRatioMCtoDataFitQCDPbPb0020,"MC weighted/QCD fit to Data (0.4 <pT<8)","p");
      if (histoRatioMCUnweightedtoDataFitQCDPbPb0020) legendFit0020->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb0020,"MC/QCD fit to Data (0.4 <pT<8)","p");
      legendFit0020->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFraction2->Update();
      canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit0020.%s",outputDir.Data(),suffix.Data()));
   }
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 0-40 % *****************
   //**********************************************************************************
   
   if (directoryPi0PbPb0040){
      canvasSpectraMCAllTogether->cd();
      canvasSpectraMCAllTogether->SetLogx(1);
      histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
      histo2DSpectraMCAll->DrawCopy(); 
   
      histoMCYieldPt0040->SetMarkerStyle(markerStylePbPb6080MC);
      histoMCYieldPt0040->Draw("hist,pe1,same");
   
      histoPCMPi0CorrectedSpec0040->Draw("hist,pe1,same");
      TF1* fitYieldDataQCDPi0PbPb0040 = FitObject("qcd","fitYieldDataQCDPi0PbPb0040","Pi0",histoPCMPi0CorrectedSpec0040,0.4,8,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb0040, 1, 1.5, colorCombPbPb6080);
   //    fitYieldDataQCDPi0PbPb0040->SetRange(0.1,20);
      fitYieldDataQCDPi0PbPb0040->Draw("same");
      forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb0040);
      fileFinalResults<< forOutput.Data()<< endl;  
      canvasSpectraMCAllTogether->Update();
      canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted0040.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDPbPb0040 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec0040, fitYieldDataQCDPi0PbPb0040);
      TH1D* histoRatioMCtoDataFitQCDPbPb0040 = CalculateHistoRatioToFit (histoMCYieldPt0040, fitYieldDataQCDPi0PbPb0040);
      TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb0040 = NULL;
      if (histoMCYieldPtWOWeights0040) histoRatioMCUnweightedtoDataFitQCDPbPb0040 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights0040, fitYieldDataQCDPi0PbPb0040);
      canvasFraction2->cd();
      if (histoRatioMCUnweightedtoDataFitQCDPbPb0040) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb0040, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb0040, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
      DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb0040, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb0040,
                  "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 0-40%",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 7.9);
      histoRatioDatatoFitQCDPbPb0040->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb0040->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDPbPb0040) histoRatioMCUnweightedtoDataFitQCDPbPb0040->Draw("same,e,p");  
      TLegend* legendFit0040 = new TLegend(0.16,0.81,0.4,0.9);
      legendFit0040->SetFillColor(0);
      legendFit0040->SetLineColor(0);
      legendFit0040->SetTextSize(0.025);
   //    legendFit0040->SetNColumns(3);
      legendFit0040->SetMargin(0.2);
      legendFit0040->AddEntry(histoRatioDatatoFitQCDPbPb0040,"Data/QCD fit to Data (0.4 <pT<8)","p");
      if (runDrawReweighted) legendFit0040->AddEntry(histoRatioMCtoDataFitQCDPbPb0040,"MC weighted/QCD fit to Data (0.4 <pT<8)","p");
      if (histoRatioMCUnweightedtoDataFitQCDPbPb0040) legendFit0040->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb0040,"MC/QCD fit to Data (0.4 <pT<8)","p");
      legendFit0040->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFraction2->Update();
      canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit0040.%s",outputDir.Data(),suffix.Data()));
   }
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 0-80 % *****************
   //**********************************************************************************
   
   if (directoryPi0PbPb0080){
      canvasSpectraMCAllTogether->cd();
      canvasSpectraMCAllTogether->SetLogx(1);
      histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
      histo2DSpectraMCAll->DrawCopy(); 
   
      histoMCYieldPt0080->SetMarkerStyle(markerStylePbPb6080MC);
      histoMCYieldPt0080->Draw("hist,pe1,same");
   
      histoPCMPi0CorrectedSpec0080->Draw("hist,pe1,same");
      TF1* fitYieldDataQCDPi0PbPb0080 = FitObject("qcd","fitYieldDataQCDPi0PbPb0080","Pi0",histoPCMPi0CorrectedSpec0080,0.4,8,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb0080, 1, 1.5, colorCombPbPb6080);
   //    fitYieldDataQCDPi0PbPb0080->SetRange(0.1,20);
      fitYieldDataQCDPi0PbPb0080->Draw("same");
      forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb0080);
      fileFinalResults<< forOutput.Data()<< endl;  
      canvasSpectraMCAllTogether->Update();
      canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted0080.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDPbPb0080 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec0080, fitYieldDataQCDPi0PbPb0080);
      TH1D* histoRatioMCtoDataFitQCDPbPb0080 = CalculateHistoRatioToFit (histoMCYieldPt0080, fitYieldDataQCDPi0PbPb0080);
      TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb0080 = NULL;
      if (histoMCYieldPtWOWeights0080) histoRatioMCUnweightedtoDataFitQCDPbPb0080 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights0080, fitYieldDataQCDPi0PbPb0080);
      canvasFraction2->cd();
      if (histoRatioMCUnweightedtoDataFitQCDPbPb0080) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb0080, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb0080, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
      DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb0080, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb0080,
                  "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 0-80%",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 7.9);
      histoRatioDatatoFitQCDPbPb0080->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb0080->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDPbPb0080) histoRatioMCUnweightedtoDataFitQCDPbPb0080->Draw("same,e,p");  
      TLegend* legendFit0080 = new TLegend(0.16,0.81,0.4,0.9);
      legendFit0080->SetFillColor(0);
      legendFit0080->SetLineColor(0);
      legendFit0080->SetTextSize(0.025);
   //    legendFit0080->SetNColumns(3);
      legendFit0080->SetMargin(0.2);
      legendFit0080->AddEntry(histoRatioDatatoFitQCDPbPb0080,"Data/QCD fit to Data (0.4 <pT<8)","p");
      if (runDrawReweighted) legendFit0080->AddEntry(histoRatioMCtoDataFitQCDPbPb0080,"MC weighted/QCD fit to Data (0.4 <pT<8)","p");
      if (histoRatioMCUnweightedtoDataFitQCDPbPb0080) legendFit0080->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb0080,"MC/QCD fit to Data (0.4 <pT<8)","p");
      legendFit0080->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFraction2->Update();
      canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit0080.%s",outputDir.Data(),suffix.Data()));
   }
   //**********************************************************************************
   //**************************** Pi0 reweighting evalulation 40-80 % *****************
   //**********************************************************************************
   if (directoryPi0PbPb4080){
      canvasSpectraMCAllTogether->cd();
      canvasSpectraMCAllTogether->SetLogx(1);
      histo2DSpectraMCAll->GetXaxis()->SetRangeUser(0.01,30);
      histo2DSpectraMCAll->DrawCopy(); 
   
      histoMCYieldPt4080->SetMarkerStyle(markerStylePbPb6080MC);
      histoMCYieldPt4080->Draw("hist,pe1,same");
   
      histoPCMPi0CorrectedSpec4080->Draw("hist,pe1,same");
      TF1* fitYieldDataQCDPi0PbPb4080 = FitObject("qcd","fitYieldDataQCDPi0PbPb4080","Pi0",histoPCMPi0CorrectedSpec4080,0.4,8,NULL,"QNRME+");
      DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb4080, 1, 1.5, colorCombPbPb6080);
   //    fitYieldDataQCDPi0PbPb4080->SetRange(0.1,20);
      fitYieldDataQCDPi0PbPb4080->Draw("same");
      forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPb4080);
      fileFinalResults<< forOutput.Data()<< endl;  
      canvasSpectraMCAllTogether->Update();
      canvasSpectraMCAllTogether->Print(Form("%s/Pi0_MCInputSpectraFitted4080.%s",outputDir.Data(),suffix.Data()));
      
      TH1D* histoRatioDatatoFitQCDPbPb4080 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpec4080, fitYieldDataQCDPi0PbPb4080);
      TH1D* histoRatioMCtoDataFitQCDPbPb4080 = CalculateHistoRatioToFit (histoMCYieldPt4080, fitYieldDataQCDPi0PbPb4080);
      TH1D* histoRatioMCUnweightedtoDataFitQCDPbPb4080 = NULL;
      if (histoMCYieldPtWOWeights4080) histoRatioMCUnweightedtoDataFitQCDPbPb4080 = CalculateHistoRatioToFit (histoMCYieldPtWOWeights4080, fitYieldDataQCDPi0PbPb4080);
      canvasFraction2->cd();
      if (histoRatioMCUnweightedtoDataFitQCDPbPb4080) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPb4080, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
      DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPb4080, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
      DrawGammaSetMarker(histoRatioDatatoFitQCDPbPb4080, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
      DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPb4080,
                  "", "p_{T} (GeV/c)", "Spectrum/ fit to Spectrum, PbPb 40-80%",
                  kFALSE, 1.5, 0, kTRUE,
                  kTRUE, -0.5, 8.,
                  kTRUE, 0., 7.9);
      histoRatioDatatoFitQCDPbPb4080->Draw("same,e,p");  
      if (runDrawReweighted) histoRatioMCtoDataFitQCDPbPb4080->Draw("same,e,p");  
      if (histoRatioMCUnweightedtoDataFitQCDPbPb4080) histoRatioMCUnweightedtoDataFitQCDPbPb4080->Draw("same,e,p");  
      TLegend* legendFit4080 = new TLegend(0.16,0.81,0.4,0.9);
      legendFit4080->SetFillColor(0);
      legendFit4080->SetLineColor(0);
      legendFit4080->SetTextSize(0.025);
   //    legendFit4080->SetNColumns(3);
      legendFit4080->SetMargin(0.2);
      legendFit4080->AddEntry(histoRatioDatatoFitQCDPbPb4080,"Data/QCD fit to Data (0.4 <pT<8)","p");
      if (runDrawReweighted) legendFit4080->AddEntry(histoRatioMCtoDataFitQCDPbPb4080,"MC weighted/QCD fit to Data (0.4 <pT<8)","p");
      if (histoRatioMCUnweightedtoDataFitQCDPbPb4080) legendFit4080->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPb4080,"MC/QCD fit to Data (0.4 <pT<8)","p");
      legendFit4080->Draw();
      
      DrawGammaLines(0., 30.,1., 1.,0.1);
      canvasFraction2->Update();
      canvasFraction2->SaveAs(Form("%s/RatioMCToMCFit4080.%s",outputDir.Data(),suffix.Data()));
   }
   
   //*************************************************************************************************************
   //********************************** Eta/Pi0 ratio for different collisions Systems ***************************
   //*************************************************************************************************************
   
   TGraphAsymmErrors* graphRatioEtaPi0StatErr7TeV = new TGraphAsymmErrors(histoRatioEtaPi07TeV);
   graphRatioEtaPi0StatErr7TeV->RemovePoint(0);
 
   TGraphAsymmErrors* graphRatioEtaPi0StatErr2760GeV = new TGraphAsymmErrors(histoRatioEtaPi02760GeV);
   graphRatioEtaPi0StatErr2760GeV->RemovePoint(0);
 
   TGraphAsymmErrors* graphRatioEtaPi0StatErr900GeV = new TGraphAsymmErrors(histoRatioEtaPi0900GeV);
   graphRatioEtaPi0StatErr900GeV->RemovePoint(0);
   
   TGraphAsymmErrors* graphRatioEtaPi0StatErrpPb = new TGraphAsymmErrors(histoPCMEtaPi0RatiopPb);
   graphRatioEtaPi0StatErrpPb->RemovePoint(0);
   
   TCanvas* canvasRatioEtaPi0DiffEnergiesSep = new TCanvas("canvasRatioEtaPi0DiffEnergiesSep","",200,10,1200,1100);// gives the page size
   DrawGammaCanvasSettings( canvasRatioEtaPi0DiffEnergiesSep,0.08, 0.01, 0.015, 0.08);
   
   
   TH2D *histo2DRatioEtaPi0DiffEnergiesSep = new TH2D("histo2DRatioEtaPi0DiffEnergiesSep", "histo2DRatioEtaPi0DiffEnergiesSep", 1000,0.,8.,1000,0,1.1   );
   SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0DiffEnergiesSep, "p_{T} (GeV/c)","#eta/#pi^{0}", 0.03,0.04, 0.03,0.04, 0.83,0.95);
   histo2DRatioEtaPi0DiffEnergiesSep->Draw();
   
   TGraphAsymmErrors* graphRatioEtaPi0SysErr900GeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0SystErr900GeV->Clone();
   ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0SysErr900GeVNewXError, 0.5);
   DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0SysErr900GeVNewXError, markerStyleSpectrum900GeV, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeVBox);
   graphRatioEtaPi0SysErr900GeVNewXError->SetFillColor(colorPi0900GeVBox);
   graphRatioEtaPi0SysErr900GeVNewXError->SetFillStyle(1001);
   graphRatioEtaPi0SysErr900GeVNewXError->Draw("2same");

   TGraphAsymmErrors* graphRatioEtaPi0SysErr2760GeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0SystErr2760GeV->Clone();
   ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0SysErr2760GeVNewXError, 0.4);
   DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0SysErr2760GeVNewXError, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeVBox);
   graphRatioEtaPi0SysErr2760GeVNewXError->SetFillColor(colorPi02760GeVBox);
   graphRatioEtaPi0SysErr2760GeVNewXError->SetFillStyle(1001);
   graphRatioEtaPi0SysErr2760GeVNewXError->Draw("2same");

   TGraphAsymmErrors* graphRatioEtaPi0SysErr7TeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0SystErr7TeV->Clone();
   ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0SysErr7TeVNewXError, 0.3);
   DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0SysErr7TeVNewXError, markerStyleSpectrum7TeV, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeVBox);
   graphRatioEtaPi0SysErr7TeVNewXError->SetFillColor(colorPi07TeVBox);
   graphRatioEtaPi0SysErr7TeVNewXError->SetFillStyle(1001);
   graphRatioEtaPi0SysErr7TeVNewXError->Draw("2same");

   TGraphAsymmErrors* graphRatioEtaPi0SysErrpPbNewXError = (TGraphAsymmErrors*)graphPCMEtaPi0RatioSysErrpPb->Clone();
   ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0SysErrpPbNewXError, 0.2);
   DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0SysErrpPbNewXError, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPbBox);
   graphRatioEtaPi0SysErrpPbNewXError->SetFillColor(colorCombpPbBox);
   graphRatioEtaPi0SysErrpPbNewXError->SetFillStyle(1001);
   graphRatioEtaPi0SysErrpPbNewXError->Draw("2same");

   TGraphAsymmErrors* graphRatioEtaPi0StatErr900GeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0StatErr900GeV->Clone();
   ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0StatErr900GeVNewXError, 0.);
   DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0StatErr900GeVNewXError, markerStyleSpectrum900GeV, markerSizePi0PP900GeV, colorPi0900GeV, colorPi0900GeV);
   graphRatioEtaPi0StatErr900GeVNewXError->Draw("p,same");

   TGraphAsymmErrors* graphRatioEtaPi0StatErr2760GeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0StatErr2760GeV->Clone();
   ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0StatErr2760GeVNewXError, 0.);
   DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0StatErr2760GeVNewXError, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);
   graphRatioEtaPi0StatErr2760GeVNewXError->Draw("p,same");

   TGraphAsymmErrors* graphRatioEtaPi0StatErr7TeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0StatErr7TeV->Clone();
   ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0StatErr7TeVNewXError, 0.);
   DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0StatErr7TeVNewXError, markerStyleSpectrum7TeV, markerSizePi0PP7TeV, colorPi07TeV, colorPi07TeV);
   graphRatioEtaPi0StatErr7TeVNewXError->Draw("p,same");

   TGraphAsymmErrors* graphRatioEtaPi0StatErrpPbNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0StatErrpPb->Clone();
   ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0StatErrpPbNewXError, 0.);
   DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0StatErrpPbNewXError, markerStylepPb, markerSizepPb, colorCombpPb, colorCombpPb);
   graphRatioEtaPi0StatErrpPbNewXError->Draw("p,same");

   TLegend* legendRatioEtaToPi0DiffEnergiesSep = new TLegend(0.13,0.75,0.45,0.95);
   legendRatioEtaToPi0DiffEnergiesSep->SetTextSize(0.03);         
   legendRatioEtaToPi0DiffEnergiesSep->SetFillColor(0);
   legendRatioEtaToPi0DiffEnergiesSep->SetBorderSize(0);
   legendRatioEtaToPi0DiffEnergiesSep->AddEntry(graphRatioEtaPi0SysErr7TeVNewXError,Form("PCM, %s", collisionSystemPP7TeV.Data()),"fp");
   legendRatioEtaToPi0DiffEnergiesSep->AddEntry(graphRatioEtaPi0SysErr2760GeVNewXError,Form("PCM, %s", collisionSystemPP2760GeV.Data()),"fp");
   legendRatioEtaToPi0DiffEnergiesSep->AddEntry(graphRatioEtaPi0SysErr900GeVNewXError,Form("PCM, %s", collisionSystemPP900GeV.Data()),"fp");
      legendRatioEtaToPi0DiffEnergiesSep->AddEntry(graphRatioEtaPi0SysErrpPbNewXError,Form("PCM, %s", collisionSystempPb.Data()),"fp");
   legendRatioEtaToPi0DiffEnergiesSep->Draw();
   
   canvasRatioEtaPi0DiffEnergiesSep->Update();
   canvasRatioEtaPi0DiffEnergiesSep->SaveAs(Form("%s/EtaPi0RatioDiffEnergiesSepPCM.%s",outputDir.Data(), suffix.Data()));

   
   TFile* fileRatioK0Correction = new TFile("ExternalInputPbPb/NeutralKaons/ratio_MCLHC11a10a_bis_K0s.root");
   TH1D* histoRatioMCtoDataK0s0005 = (TH1D*)fileRatioK0Correction->Get("ratio_cent0");
   TH1D* histoRatioDatatoMCK0s0005 = InvertHisto(histoRatioMCtoDataK0s0005);
   TH1D* histoRatioMCtoDataK0s0510 = (TH1D*)fileRatioK0Correction->Get("ratio_cent5");
   TH1D* histoRatioDatatoMCK0s0510 = InvertHisto(histoRatioMCtoDataK0s0510);
   TH1D* histoRatioDatatoMCK0s0010 = (TH1D*)histoRatioDatatoMCK0s0005->Clone("histoRatioDatatoMCK0s0010");
   for (Int_t i = 1; i < histoRatioDatatoMCK0s0010->GetNbinsX()+1; i++){
      histoRatioDatatoMCK0s0010->SetBinContent(i,(histoRatioDatatoMCK0s0005->GetBinContent(i) + histoRatioDatatoMCK0s0510->GetBinContent(i))/2);
      histoRatioDatatoMCK0s0010->SetBinError(i,(histoRatioDatatoMCK0s0005->GetBinError(i) + histoRatioDatatoMCK0s0510->GetBinError(i))/2);
   }
   TH1D* histoRatioMCtoDataK0s1020 = (TH1D*)fileRatioK0Correction->Get("ratio_cent10");
   TH1D* histoRatioDatatoMCK0s1020 = InvertHisto(histoRatioMCtoDataK0s1020);
   TH1D* histoRatioDatatoMCK0s0020 = (TH1D*)histoRatioDatatoMCK0s0010->Clone("histoRatioDatatoMCK0s0020");
   for (Int_t i = 1; i < histoRatioDatatoMCK0s0020->GetNbinsX()+1; i++){
      histoRatioDatatoMCK0s0020->SetBinContent(i,(histoRatioDatatoMCK0s0010->GetBinContent(i) + histoRatioDatatoMCK0s1020->GetBinContent(i))/2);
      histoRatioDatatoMCK0s0020->SetBinError(i,(histoRatioDatatoMCK0s0010->GetBinError(i) + histoRatioDatatoMCK0s1020->GetBinError(i))/2);
   }
   
   TH1D* histoRatioMCtoDataK0s2040 = (TH1D*)fileRatioK0Correction->Get("ratio_cent20");
   TH1D* histoRatioDatatoMCK0s2040 = InvertHisto(histoRatioMCtoDataK0s2040);
   TH1D* histoRatioDatatoMCK0s0040 = (TH1D*)histoRatioDatatoMCK0s0020->Clone("histoRatioDatatoMCK0s0040");
   for (Int_t i = 1; i < histoRatioDatatoMCK0s0040->GetNbinsX()+1; i++){
      histoRatioDatatoMCK0s0040->SetBinContent(i,(histoRatioDatatoMCK0s0020->GetBinContent(i) + histoRatioDatatoMCK0s2040->GetBinContent(i))/2);
      histoRatioDatatoMCK0s0040->SetBinError(i,(histoRatioDatatoMCK0s0020->GetBinError(i) + histoRatioDatatoMCK0s2040->GetBinError(i))/2);
   }
   TH1D* histoRatioMCtoDataK0s4060 = (TH1D*)fileRatioK0Correction->Get("ratio_cent40");
   TH1D* histoRatioDatatoMCK0s4060 = InvertHisto(histoRatioMCtoDataK0s4060);
   TH1D* histoRatioMCtoDataK0s6080 = (TH1D*)fileRatioK0Correction->Get("ratio_cent60");
   TH1D* histoRatioDatatoMCK0s6080 = InvertHisto(histoRatioMCtoDataK0s6080);
   TH1D* histoRatioDatatoMCK0s4080 = (TH1D*)histoRatioDatatoMCK0s4060->Clone("histoRatioDatatoMCK0s4080");
   for (Int_t i = 1; i < histoRatioDatatoMCK0s4080->GetNbinsX()+1; i++){
      histoRatioDatatoMCK0s4080->SetBinContent(i,(histoRatioDatatoMCK0s4060->GetBinContent(i) + histoRatioDatatoMCK0s4060->GetBinContent(i))/2);
      histoRatioDatatoMCK0s4080->SetBinError(i,(histoRatioDatatoMCK0s4060->GetBinError(i) + histoRatioDatatoMCK0s4060->GetBinError(i))/2);
   }
   TH1D* histoRatioDatatoMCK0s0080 = (TH1D*)histoRatioDatatoMCK0s4060->Clone("histoRatioDatatoMCK0s0080");
   for (Int_t i = 1; i < histoRatioDatatoMCK0s0080->GetNbinsX()+1; i++){
      histoRatioDatatoMCK0s0080->SetBinContent(i,(histoRatioDatatoMCK0s0040->GetBinContent(i) + histoRatioDatatoMCK0s4080->GetBinContent(i))/2);
      histoRatioDatatoMCK0s0080->SetBinError(i,(histoRatioDatatoMCK0s0040->GetBinError(i) + histoRatioDatatoMCK0s4080->GetBinError(i))/2);
   }
   
   
   TFile fMCSpectraInput("MCSpectraInput.root","RECREATE");
      histoMCYieldPtWOWeights0005->Write("Pi0_Hijing_PbPb_2760GeV_0005");
      histoMCYieldPtWOWeights0510->Write("Pi0_Hijing_PbPb_2760GeV_0510");
      histoMCYieldPtWOWeights0010->Write("Pi0_Hijing_PbPb_2760GeV_0010");
      histoMCYieldPtWOWeights1020->Write("Pi0_Hijing_PbPb_2760GeV_1020");
      histoMCYieldPtWOWeights2040->Write("Pi0_Hijing_PbPb_2760GeV_2040");
      histoMCYieldPtWOWeights4060->Write("Pi0_Hijing_PbPb_2760GeV_4060");
      histoMCYieldPtWOWeights6080->Write("Pi0_Hijing_PbPb_2760GeV_6080");
      if (directoryPi0PbPb0020) histoMCYieldPtWOWeights0020->Write("Pi0_Hijing_PbPb_2760GeV_0020");
      if (directoryPi0PbPb0040) histoMCYieldPtWOWeights0040->Write("Pi0_Hijing_PbPb_2760GeV_0040");
      if (directoryPi0PbPb0080) histoMCYieldPtWOWeights0080->Write("Pi0_Hijing_PbPb_2760GeV_0080");
      if (directoryPi0PbPb4080) histoMCYieldPtWOWeights4080->Write("Pi0_Hijing_PbPb_2760GeV_4080");
      histoRatioDatatoMCK0s0005->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0005");
      histoRatioDatatoMCK0s0010->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0010");
      histoRatioDatatoMCK0s0510->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0510");
      histoRatioDatatoMCK0s1020->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_1020");
      histoRatioDatatoMCK0s2040->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_2040");
      histoRatioDatatoMCK0s4060->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_4060");
      histoRatioDatatoMCK0s6080->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_6080");
      histoRatioDatatoMCK0s0020->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0020");
      histoRatioDatatoMCK0s0040->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0040");
      histoRatioDatatoMCK0s0080->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_0080");
      histoRatioDatatoMCK0s4080->Write("K0s_RatioDataToMC_Hijing_PbPb_2760GeV_4080");
   fMCSpectraInput.Close();
}

//    cout << "gamma corrections 2.76TeV GeV" << endl;
//    TFile* fileGammaDataPP = new TFile("0000011002093663003800000_01631031009/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_0000011002093663003800000_01631031009.root");
//    TH1D* histoTrueGammaPurityPP = (TH1D*)fileGammaDataPP->Get("fMCGammaTruePurity");
//    TH1D* histoTrueGammaEffiPP = (TH1D*)fileGammaDataPP->Get("fMCGammaRecoEff");
//    TH1D* histoConversionProbPP = (TH1D*)fileGammaDataPP->Get("fMCGammaConvProb");
//    
//    cout << "HI 0005" << endl;
//    TFile* fileGammaData0005 = new TFile("3010001042092970023220000_01522045000/PbPb_2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_3010001042092970023220000_01522045000.root");
//    TH1D* histoTrueGammaPurity0005 = (TH1D*)fileGammaData0005->Get("fMCGammaTruePurity");
//    TH1D* histoTrueGammaEffi0005 = (TH1D*)fileGammaData0005->Get("fMCGammaRecoEff");
//    TH1D* histoConversionProb0005 = (TH1D*)fileGammaData0005->Get("fMCGammaConvProb");
// 
//    cout << "HI 0010" << endl;
//    TFile* fileGammaData0010 = new TFile("1010001042092970023220000_01522045000/PbPb_2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_1010001042092970023220000_01522045000.root");
//    TH1D* histoTrueGammaPurity0010 = (TH1D*)fileGammaData0010->Get("fMCGammaTruePurity");
//    TH1D* histoTrueGammaEffi0010 = (TH1D*)fileGammaData0010->Get("fMCGammaRecoEff");
//    TH1D* histoConversionProb0010 = (TH1D*)fileGammaData0010->Get("fMCGammaConvProb");
// 
//    cout << "HI 0510" << endl;
//    TFile* fileGammaData0510 = new TFile("3120001042092970023220000_01522045000/PbPb_2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_3120001042092970023220000_01522045000.root");
//    TH1D* histoTrueGammaPurity0510 = (TH1D*)fileGammaData0510->Get("fMCGammaTruePurity");
//    TH1D* histoTrueGammaEffi0510 = (TH1D*)fileGammaData0510->Get("fMCGammaRecoEff");
//    TH1D* histoConversionProb0510 = (TH1D*)fileGammaData0510->Get("fMCGammaConvProb");
// 
//    cout << "HI 1020" << endl;
//    TFile* fileGammaData1020 = new TFile("1120001042092970023220000_01522045000/PbPb_2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_1120001042092970023220000_01522045000.root");
//    TH1D* histoTrueGammaPurity1020 = (TH1D*)fileGammaData1020->Get("fMCGammaTruePurity");
//    TH1D* histoTrueGammaEffi1020 = (TH1D*)fileGammaData1020->Get("fMCGammaRecoEff");
//    TH1D* histoConversionProb1020 = (TH1D*)fileGammaData1020->Get("fMCGammaConvProb");
// 
//    cout << "HI 2040" << endl;
//    TFile* fileGammaData2040 = new TFile("1240001042092970023220000_01522045000/PbPb_2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_1240001042092970023220000_01522045000.root");
//    TH1D* histoTrueGammaPurity2040 = (TH1D*)fileGammaData2040->Get("fMCGammaTruePurity");
//    TH1D* histoTrueGammaEffi2040 = (TH1D*)fileGammaData2040->Get("fMCGammaRecoEff");
//    TH1D* histoConversionProb2040 = (TH1D*)fileGammaData2040->Get("fMCGammaConvProb");
// 
//    cout << "HI 4060" << endl;
//    TFile* fileGammaData4060 = new TFile("1460001042092970023220000_01522065000/PbPb_2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_1460001042092970023220000_01522065000.root");
//    TH1D* histoTrueGammaPurity4060 = (TH1D*)fileGammaData4060->Get("fMCGammaTruePurity");
//    TH1D* histoTrueGammaEffi4060 = (TH1D*)fileGammaData4060->Get("fMCGammaRecoEff");
//    TH1D* histoConversionProb4060 = (TH1D*)fileGammaData4060->Get("fMCGammaConvProb");
//    
//    cout << "HI 6080" << endl;
//    TFile* fileGammaData6080 = new TFile("1680001042092970023220000_01522065000/PbPb_2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_1680001042092970023220000_01522065000.root");
//    TH1D* histoTrueGammaPurity6080 = (TH1D*)fileGammaData6080->Get("fMCGammaTruePurity");
//    TH1D* histoTrueGammaEffi6080 = (TH1D*)fileGammaData6080->Get("fMCGammaRecoEff");
//    TH1D* histoConversionProb6080 = (TH1D*)fileGammaData6080->Get("fMCGammaConvProb");
// 
// 
   
//       TCanvas* canvasConversionProbAllCent = new TCanvas("canvasConversionProbAllCent","",200,10,1350,1350);  // gives the page size
//    DrawGammaCanvasSettings( canvasConversionProbAllCent, 0.1, 0.02, 0.035, 0.09);
//       TH2F * histo2DConversionProbAllCent;
//       histo2DConversionProbAllCent = new TH2F("histo2DConversionProbAllCent","histo2DConversionProbAllCent",1000,0,10.,2000,0.,0.12 );
// //       histo2DConversionProbAllCent->GetXaxis()->SetRangeUser(0.,12.);
//       SetStyleHistoTH2ForGraphs(histo2DConversionProbAllCent, "p_{T} (GeV/c)","P",0.03,0.04, 0.03,0.04, 1.,1.);
//       histo2DConversionProbAllCent->Draw("copy");
// 
//       DrawGammaSetMarker(histoConversionProbPP, markerStylePP, markerSizePP, colorCombPP, colorCombPP);   
//       histoConversionProbPP->DrawCopy("e1,same");  
//       DrawGammaSetMarker(histoConversionProb6080, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);  
//       histoConversionProb6080->DrawCopy("e1,same");   
//       DrawGammaSetMarker(histoConversionProb0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);  
//       histoConversionProb0005->DrawCopy("e1,same");   
//       DrawGammaSetMarker(histoConversionProb0510, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0510, colorCombPbPb0510);  
//       histoConversionProb0510->DrawCopy("e1,same");   
//       DrawGammaSetMarker(histoConversionProb1020, markerStylePbPb1020, markerSizePbPb1020, colorCombPbPb1020, colorCombPbPb1020);  
//       histoConversionProb1020->DrawCopy("e1,same");   
//       DrawGammaSetMarker(histoConversionProb2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);  
//       histoConversionProb2040->DrawCopy("e1,same");   
//       DrawGammaSetMarker(histoConversionProb4060, markerStylePbPb4060, markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060);  
//       histoConversionProb4060->DrawCopy("e1,same");   
//       histoConversionProbPP->DrawCopy("e1,same");  
//       
//       TLegend* legendConversionProbAllCent = new TLegend(0.5,0.12,0.94,0.35);
//       legendConversionProbAllCent->SetFillColor(0);
//       legendConversionProbAllCent->SetLineColor(0);
//       legendConversionProbAllCent->SetTextSize(0.027);
//       legendConversionProbAllCent->AddEntry(histoConversionProb0005,"0-5% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendConversionProbAllCent->AddEntry(histoConversionProb0510,"5-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendConversionProbAllCent->AddEntry(histoConversionProb1020,"10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendConversionProbAllCent->AddEntry(histoConversionProb2040,"20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendConversionProbAllCent->AddEntry(histoConversionProb4060,"40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendConversionProbAllCent->AddEntry(histoConversionProb6080,"60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendConversionProbAllCent->AddEntry(histoConversionProbPP,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//       legendConversionProbAllCent->Draw();
//    
//    canvasConversionProbAllCent->SaveAs(Form("%s/ConversionProbability.%s",outputDir.Data(),suffix.Data()));
//    delete canvasConversionProbAllCent;
// 
//    
//    
//    TCanvas* canvasGammaEffiAllCent = new TCanvas("canvasGammaEffiAllCent","",200,10,1350,1350);  // gives the page size
//    DrawGammaCanvasSettings( canvasGammaEffiAllCent, 0.1, 0.02, 0.035, 0.09);
//       TH2F * histo2DGammaEffiAllCent;
//       histo2DGammaEffiAllCent = new TH2F("histo2DGammaEffiAllCent","histo2DGammaEffiAllCent",1000,0,10.,2000,0.,0.72 );
// //       histo2DGammaEffiAllCent->GetXaxis()->SetRangeUser(0.,12.);
//       SetStyleHistoTH2ForGraphs(histo2DGammaEffiAllCent, "p_{T} (GeV/c)","#epsilon_{reco}",0.03,0.04, 0.03,0.04, 1.,1.);
//       histo2DGammaEffiAllCent->Draw("copy");
// 
//       DrawGammaSetMarker(histoTrueGammaEffiPP, markerStylePP, markerSizePP, colorCombPP, colorCombPP); 
//       histoTrueGammaEffiPP->DrawCopy("e1,same");   
//       DrawGammaSetMarker(histoTrueGammaEffi6080, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);   
//       histoTrueGammaEffi6080->DrawCopy("e1,same");    
//       DrawGammaSetMarker(histoTrueGammaEffi0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);   
//       histoTrueGammaEffi0005->DrawCopy("e1,same");    
//       DrawGammaSetMarker(histoTrueGammaEffi0510, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0510, colorCombPbPb0510);   
//       histoTrueGammaEffi0510->DrawCopy("e1,same");    
//       DrawGammaSetMarker(histoTrueGammaEffi1020, markerStylePbPb1020, markerSizePbPb1020, colorCombPbPb1020, colorCombPbPb1020);   
//       histoTrueGammaEffi1020->DrawCopy("e1,same");    
//       DrawGammaSetMarker(histoTrueGammaEffi2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);   
//       histoTrueGammaEffi2040->DrawCopy("e1,same");    
//       DrawGammaSetMarker(histoTrueGammaEffi4060, markerStylePbPb4060, markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060);   
//       histoTrueGammaEffi4060->DrawCopy("e1,same");    
//       histoTrueGammaEffiPP->DrawCopy("e1,same");   
//       
//       TLegend* legendGammaEffiAllCent = new TLegend(0.5,0.12,0.94,0.35);
//       legendGammaEffiAllCent->SetFillColor(0);
//       legendGammaEffiAllCent->SetLineColor(0);
//       legendGammaEffiAllCent->SetTextSize(0.027);
//       legendGammaEffiAllCent->AddEntry(histoTrueGammaEffi0005,"0-5% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendGammaEffiAllCent->AddEntry(histoTrueGammaEffi0510,"5-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendGammaEffiAllCent->AddEntry(histoTrueGammaEffi1020,"10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendGammaEffiAllCent->AddEntry(histoTrueGammaEffi2040,"20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendGammaEffiAllCent->AddEntry(histoTrueGammaEffi4060,"40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendGammaEffiAllCent->AddEntry(histoTrueGammaEffi6080,"60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendGammaEffiAllCent->AddEntry(histoTrueGammaEffiPP,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//       legendGammaEffiAllCent->Draw();
//    
//    canvasGammaEffiAllCent->SaveAs(Form("%s/EffiAllCent.%s",outputDir.Data(),suffix.Data()));
//    delete canvasGammaEffiAllCent;
// 
//    
//    TCanvas* canvasPurityAllCent = new TCanvas("canvasPurityAllCent","",200,10,1350,1350);  // gives the page size
//    DrawGammaCanvasSettings( canvasPurityAllCent, 0.1, 0.02, 0.035, 0.09);
//       TH2F * histo2DPurityAllCent;
//       histo2DPurityAllCent = new TH2F("histo2DPurityAllCent","histo2DPurityAllCent",1000,0,10.,2000,0.4,1.05   );
// //       histo2DPurityAllCent->GetXaxis()->SetRangeUser(0.,12.);
//       SetStyleHistoTH2ForGraphs(histo2DPurityAllCent, "p_{T} (GeV/c)","#epsilon_{pur}",0.03,0.04, 0.03,0.04, 1.,1.);
//       histo2DPurityAllCent->Draw("copy");
// 
//       DrawGammaSetMarker(histoTrueGammaPurityPP, markerStylePP, markerSizePP, colorCombPP, colorCombPP);  
//       histoTrueGammaPurityPP->DrawCopy("e1,same");    
//       DrawGammaSetMarker(histoTrueGammaPurity6080, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080); 
//       histoTrueGammaPurity6080->DrawCopy("e1,same");  
//       DrawGammaSetMarker(histoTrueGammaPurity0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005); 
//       histoTrueGammaPurity0005->DrawCopy("e1,same");  
//       DrawGammaSetMarker(histoTrueGammaPurity0510, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0510, colorCombPbPb0510); 
//       histoTrueGammaPurity0510->DrawCopy("e1,same");  
//       DrawGammaSetMarker(histoTrueGammaPurity1020, markerStylePbPb1020, markerSizePbPb1020, colorCombPbPb1020, colorCombPbPb1020); 
//       histoTrueGammaPurity1020->DrawCopy("e1,same");  
//       DrawGammaSetMarker(histoTrueGammaPurity2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040); 
//       histoTrueGammaPurity2040->DrawCopy("e1,same");  
//       DrawGammaSetMarker(histoTrueGammaPurity4060, markerStylePbPb4060, markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060); 
//       histoTrueGammaPurity4060->DrawCopy("e1,same");  
//       histoTrueGammaPurityPP->DrawCopy("e1,same");    
//       
//       TLegend* legendPurityAllCent = new TLegend(0.5,0.12,0.94,0.35);
//       legendPurityAllCent->SetFillColor(0);
//       legendPurityAllCent->SetLineColor(0);
//       legendPurityAllCent->SetTextSize(0.027);
//       legendPurityAllCent->AddEntry(histoTrueGammaPurity0005,"0-5% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendPurityAllCent->AddEntry(histoTrueGammaPurity0510,"5-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendPurityAllCent->AddEntry(histoTrueGammaPurity1020,"10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendPurityAllCent->AddEntry(histoTrueGammaPurity2040,"20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendPurityAllCent->AddEntry(histoTrueGammaPurity4060,"40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendPurityAllCent->AddEntry(histoTrueGammaPurity6080,"60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","p");
//       legendPurityAllCent->AddEntry(histoTrueGammaPurityPP,"pp #sqrt{#it{s}} = 2.76 TeV","p");
//       legendPurityAllCent->Draw();
//    
//    canvasPurityAllCent->SaveAs(Form("%s/PurityAllCent.%s",outputDir.Data(),suffix.Data()));
//    delete canvasPurityAllCent;
   