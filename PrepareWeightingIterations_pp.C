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

extern TRandom*  	 gRandom;
extern TBenchmark*   gBenchmark;
extern TSystem*  	 gSystem;
extern TMinuit*  	 gMinuit;

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


void PrepareWeightingIterations_pp(TString suffix = "pdf", TString nameFilepp2760GeV = "data_PCMResults_pp.root", TString nameFilepp8TeV = "data_PCMResults_pp.root", Bool_t runDrawReweighted = kTRUE, TString stringIterationNumber = ""){

	gROOT->Reset();   
	gROOT->SetStyle("Plain");

	TString dateForOutput 											= ReturnDateStringForOutput();
	TString outputDir 												= Form("%s/%s/PrepareWeightingIterations_pp",suffix.Data(),dateForOutput.Data());
	
	gSystem->Exec("mkdir -p "+outputDir);
	gSystem->Exec(Form("cp %s %s/InputFilePCMPionpp2760GeV.root ",nameFilepp2760GeV.Data(),outputDir.Data() ));
	gSystem->Exec(Form("cp %s %s/InputFilePCMPionpp8TeV.root ",nameFilepp8TeV.Data(),outputDir.Data() ));

	StyleSettingsThesis();  
	SetPlotStyle();
	
// 	Color_t 	color900GeV 										= kRed +2;
	Color_t 	color2760GeV 										= kMagenta+2;
	Color_t 	color8TeV 											= kGreen+2;
// 	Color_t		color7TeV											= kBlue+2;
// 	Color_t 	color900GeVBox 										= color900GeV-10;
// 	Color_t 	color2760GeVBox 									= color2760GeV-10;
// 	Color_t 	color7TeVBox 										= color7TeV-10;

// 	Color_t 	colorMCPythiaPP900GeV 								= color900GeV-4;
	Color_t 	colorMCPythiaPP2760GeV 								= color2760GeV+2;
	Color_t 	colorMCPythia2PP2760GeV 							= color2760GeV+1;
	Color_t 	colorMCPythiaAddSigPP2760GeV 						= color2760GeV-6;
// 	Color_t 	colorMCPythiaPP7TeV 								= color7TeV+3;
// 	Color_t 	colorMCPhojetPP900GeV 								= color900GeV+2;
	Color_t 	colorMCPhojetPP2760GeV 								= color2760GeV-4;
// 	Color_t 	colorMCPhojetPP7TeV 								= color7TeV-3;

	Color_t 	colorMCPythiaPP8TeV 								= color8TeV+2;
	Color_t 	colorMCPythia2PP8TeV 								= color8TeV+1;
	Color_t 	colorMCPythiaAddSigPP8TeV 							= color8TeV-6;
	Color_t 	colorMCPhojetPP8TeV 								= color8TeV-4;
	
	
// 	Style_t 	markerStyleSpectrum7TeV 							= 20 ;
// 	Style_t 	markerStyleSpectrum900GeV 							= 21 ;
	Style_t 	markerStyleSpectrum2760GeV 							= 29 ;
	Style_t 	markerStyleSpectrum8TeV 							= 33 ;

// 	Size_t 		markerSizePP7TeV 									= 1.8;
// 	Size_t 		markerSizePP900GeV 									= 1.8;
	Size_t 		markerSizePP2760GeV 								= 2.2;
	Size_t 		markerSizePP8TeV 									= 2.2;

	TString 	collisionSystemPP2760GeV 							= "pp #sqrt{#it{s}} = 2.76 TeV";      
	TString 	collisionSystemPP8TeV 								= "pp #sqrt{#it{s}} = 8 TeV";      
// 	TString 	collisionSystemPP7TeV 								= "pp #sqrt{#it{s}} = 7 TeV";      
// 	TString 	collisionSystemPP900GeV 							= "pp #sqrt{#it{s}} = 0.9 TeV";     

	TString 	nameHistoPCM 										= "CorrectedYieldPi0";
	TString 	nameHistoPCMEta 									= "CorrectedYieldEta";

	Style_t  	lineStyleMCA										= 1;
	Style_t  	lineStyleMCB										= 2;
	Style_t  	lineStyleMCC										= 3;
	Style_t  	lineStyleMCAddSig									= 4 ;
	
	Style_t		markerStyleMCA										= 20;
	Style_t		markerStyleMCB										= 25;
	Style_t		markerStyleMCC										= 30;
	Style_t		markerStyleMCAddSig									= 29;
	
	// read reconstructed data
	TFile* 		filePCMpp2760GeV 									= new TFile(nameFilepp2760GeV);
	TDirectory*	directoryPCMPi02760GeV 								= (TDirectory*)filePCMpp2760GeV->Get("Pi02.76TeV"); 
	TH1D* 		histoPCMYieldPi02760GeV 							= (TH1D*)directoryPCMPi02760GeV->Get(nameHistoPCM.Data());
	TH1D*		histoPi0InputFullReweighted2760GeV					= (TH1D*)directoryPCMPi02760GeV->Get("Pi0_Input_Reweighted");
	TH1D*		histoPi0InputFull2760GeV							= (TH1D*)directoryPCMPi02760GeV->Get("Pi0_Input");
	TH1D*		histoPi0InputFullAddSigReweighted2760GeV			= (TH1D*)directoryPCMPi02760GeV->Get("Pi0_Input_Reweighted_AddedSig");
	TH1D*		histoPi0InputFullAddSig2760GeV						= (TH1D*)directoryPCMPi02760GeV->Get("Pi0_Input_AddedSig");
	TDirectory*	directoryPCMEta2760GeV 								= (TDirectory*)filePCMpp2760GeV->Get("Eta2.76TeV"); 
	TH1D* 		histoPCMYieldEta2760GeV 							= (TH1D*)directoryPCMEta2760GeV->Get(nameHistoPCMEta.Data());
	TH1D*		histoEtaInputFullReweighted2760GeV					= (TH1D*)directoryPCMEta2760GeV->Get("Eta_Input_Reweighted");
	TH1D*		histoEtaInputFull2760GeV							= (TH1D*)directoryPCMEta2760GeV->Get("Eta_Input");
	TH1D*		histoEtaInputFullAddSigReweighted2760GeV			= (TH1D*)directoryPCMEta2760GeV->Get("Eta_Input_Reweighted_AddedSig");
	TH1D*		histoEtaInputFullAddSig2760GeV						= (TH1D*)directoryPCMEta2760GeV->Get("Eta_Input_AddedSig");
	// Pythia 8
	TFile* 		filePi0PCM2760GeV_LHC12f1a_PYT8 					= new TFile("0000011_002000092970028250400000_01525065000000_LHC12f1a-woSDD/2.76TeV/Pi0_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D* 		histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV 		= (TH1D*)filePi0PCM2760GeV_LHC12f1a_PYT8->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoPi0Efficiency_LHC12f1a_PYT8_2760GeV 			= (TH1D*)filePi0PCM2760GeV_LHC12f1a_PYT8->Get("TrueMesonEffiPt");
	TFile* 		fileEtaPCM2760GeV_LHC12f1a_PYT8 					= new TFile("0000011_002000092970028250400000_01525065000000_LHC12f1a-woSDD/2.76TeV/Eta_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D* 		histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV 		= (TH1D*)fileEtaPCM2760GeV_LHC12f1a_PYT8->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoEtaEfficiency_LHC12f1a_PYT8_2760GeV 			= (TH1D*)fileEtaPCM2760GeV_LHC12f1a_PYT8->Get("TrueMesonEffiPt");

	// Phojet
	TFile* 		filePi0PCM2760GeV_LHC12f1b_PHO 						= new TFile("0000011_002000092970028250400000_01525065000000_LHC12f1b-woSDD/2.76TeV/Pi0_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D* 		histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV 		= (TH1D*)filePi0PCM2760GeV_LHC12f1b_PHO->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoPi0Efficiency_LHC12f1b_PHO_2760GeV 			= (TH1D*)filePi0PCM2760GeV_LHC12f1b_PHO->Get("TrueMesonEffiPt");
	TFile*		fileEtaPCM2760GeV_LHC12f1b_PHO 						= new TFile("0000011_002000092970028250400000_01525065000000_LHC12f1b-woSDD/2.76TeV/Eta_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D*		histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV 		= (TH1D*)fileEtaPCM2760GeV_LHC12f1b_PHO->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoEtaEfficiency_LHC12f1b_PHO_2760GeV 			= (TH1D*)fileEtaPCM2760GeV_LHC12f1b_PHO->Get("TrueMesonEffiPt");

	// Pythia 8, LHC12i3
	TFile*		filePi0PCM2760GeV_LHC12i3_PYT8 						= new TFile("0000011_002000092970028250400000_01525065000000_LHC12i3-woSDD/2.76TeV/Pi0_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D*		histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV 		= (TH1D*)filePi0PCM2760GeV_LHC12i3_PYT8->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoPi0Efficiency_LHC12i3_PYT8_2760GeV 			= (TH1D*)filePi0PCM2760GeV_LHC12i3_PYT8->Get("TrueMesonEffiPt");
	TFile*		fileEtaPCM2760GeV_LHC12i3_PYT8 						= new TFile("0000011_002000092970028250400000_01525065000000_LHC12i3-woSDD/2.76TeV/Eta_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D*		histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV 		= (TH1D*)fileEtaPCM2760GeV_LHC12i3_PYT8->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoEtaEfficiency_LHC12i3_PYT8_2760GeV 			= (TH1D*)fileEtaPCM2760GeV_LHC12i3_PYT8->Get("TrueMesonEffiPt");

	// Pythia 8, LHC12i3 added signals
	TFile*		filePi0PCM2760GeV_LHC12i3_PYT8Addsig 				= new TFile("0000012_002000092970028250400000_01525065000000_LHC12i3-woSDD/2.76TeV/Pi0_MC_GammaConvV1Correction_0000012_002000092970028250400000_01525065000000.root");
	TH1D*		histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV = (TH1D*)filePi0PCM2760GeV_LHC12i3_PYT8Addsig->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV 		= (TH1D*)filePi0PCM2760GeV_LHC12i3_PYT8Addsig->Get("TrueMesonEffiPt");
	TFile*		fileEtaPCM2760GeV_LHC12i3_PYT8Addsig 				= new TFile("0000012_002000092970028250400000_01525065000000_LHC12i3-woSDD/2.76TeV/Eta_MC_GammaConvV1Correction_0000012_002000092970028250400000_01525065000000.root");
	TH1D* 		histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV = (TH1D*)fileEtaPCM2760GeV_LHC12i3_PYT8Addsig->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV 		= (TH1D*)fileEtaPCM2760GeV_LHC12i3_PYT8Addsig->Get("TrueMesonEffiPt");

	// Pythia 8
	TFile* 		filePi0PCM2760GeV_LHC12f1a_PYT8_wSDD 				= new TFile("0000011_002000092970028250400000_01525065000000_LHC12f1a-wSDD/2.76TeV/Pi0_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D* 		histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD = (TH1D*)filePi0PCM2760GeV_LHC12f1a_PYT8_wSDD->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDD 		= (TH1D*)filePi0PCM2760GeV_LHC12f1a_PYT8_wSDD->Get("TrueMesonEffiPt");
	TFile* 		fileEtaPCM2760GeV_LHC12f1a_PYT8_wSDD 				= new TFile("0000011_002000092970028250400000_01525065000000_LHC12f1a-wSDD/2.76TeV/Eta_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D* 		histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD	= (TH1D*)fileEtaPCM2760GeV_LHC12f1a_PYT8_wSDD->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDD		= (TH1D*)fileEtaPCM2760GeV_LHC12f1a_PYT8_wSDD->Get("TrueMesonEffiPt");

	// Phojet
	TFile* 		filePi0PCM2760GeV_LHC12f1b_PHO_wSDD 				= new TFile("0000011_002000092970028250400000_01525065000000_LHC12f1b-wSDD/2.76TeV/Pi0_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D* 		histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD 	= (TH1D*)filePi0PCM2760GeV_LHC12f1b_PHO_wSDD->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDD 		= (TH1D*)filePi0PCM2760GeV_LHC12f1b_PHO_wSDD->Get("TrueMesonEffiPt");
	TFile*		fileEtaPCM2760GeV_LHC12f1b_PHO_wSDD 				= new TFile("0000011_002000092970028250400000_01525065000000_LHC12f1b-wSDD/2.76TeV/Eta_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D*		histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD 	= (TH1D*)fileEtaPCM2760GeV_LHC12f1b_PHO_wSDD->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDD 		= (TH1D*)fileEtaPCM2760GeV_LHC12f1b_PHO_wSDD->Get("TrueMesonEffiPt");

	// Pythia 8, LHC12i3
	TFile*		filePi0PCM2760GeV_LHC12i3_PYT8_wSDD 						= new TFile("0000011_002000092970028250400000_01525065000000_LHC12i3-wSDD/2.76TeV/Pi0_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D*		histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD	= (TH1D*)filePi0PCM2760GeV_LHC12i3_PYT8_wSDD->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDD		= (TH1D*)filePi0PCM2760GeV_LHC12i3_PYT8_wSDD->Get("TrueMesonEffiPt");
	TFile*		fileEtaPCM2760GeV_LHC12i3_PYT8_wSDD 						= new TFile("0000011_002000092970028250400000_01525065000000_LHC12i3-wSDD/2.76TeV/Eta_data_GammaConvV1Correction_0000011_002000092970028250400000_01525065000000.root");
	TH1D*		histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD 	= (TH1D*)fileEtaPCM2760GeV_LHC12i3_PYT8_wSDD->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDD 		= (TH1D*)fileEtaPCM2760GeV_LHC12i3_PYT8_wSDD->Get("TrueMesonEffiPt");

	// Pythia 8, LHC12i3 added signals
	TFile*		filePi0PCM2760GeV_LHC12i3_PYT8Addsig_wSDD 				= new TFile("0000012_002000092970028250400000_01525065000000_LHC12i3-wSDD/2.76TeV/Pi0_MC_GammaConvV1Correction_0000012_002000092970028250400000_01525065000000.root");
	TH1D*		histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD = (TH1D*)filePi0PCM2760GeV_LHC12i3_PYT8Addsig_wSDD->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDD 	= (TH1D*)filePi0PCM2760GeV_LHC12i3_PYT8Addsig_wSDD->Get("TrueMesonEffiPt");
	TFile*		fileEtaPCM2760GeV_LHC12i3_PYT8Addsig_wSDD			= new TFile("0000012_002000092970028250400000_01525065000000_LHC12i3-wSDD/2.76TeV/Eta_MC_GammaConvV1Correction_0000012_002000092970028250400000_01525065000000.root");
	TH1D* 		histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD = (TH1D*)fileEtaPCM2760GeV_LHC12i3_PYT8Addsig_wSDD->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDD 	= (TH1D*)fileEtaPCM2760GeV_LHC12i3_PYT8Addsig_wSDD->Get("TrueMesonEffiPt");
	
// 	Width_t		widthLinesBoxes										= 1.4;
// 	Width_t		widthCommonFit									 	= 2.;
// 	Width_t		widthStatErrBars									= 1.5;
// 	Width_t		widthCommonErrors 									= 1.1;
// 	Width_t		widthCommonSpectrumBoxes 							= 0.99;
// 	if (suffix.CompareTo("eps")!=0){
// 		widthLinesBoxes 											= 2.3;
// 		widthCommonFit 												= 2.6;
// 		widthStatErrBars 											= 2.6;
// 		widthCommonErrors 											= 2.;
// 		widthCommonSpectrumBoxes 									= 2.3;
// 	}

	// read reconstructed data
	TFile* 		filePCMpp8TeV 										= new TFile(nameFilepp8TeV);
	TDirectory*	directoryPCMPi08TeV 								= (TDirectory*)filePCMpp8TeV->Get("Pi08TeV"); 
	TH1D* 		histoPCMYieldPi08TeV 								= (TH1D*)directoryPCMPi08TeV->Get(nameHistoPCM.Data());
	TH1D*		histoPi0InputFullReweighted8TeV						= (TH1D*)directoryPCMPi08TeV->Get("Pi0_Input_Reweighted");
	TH1D*		histoPi0InputFull8TeV								= (TH1D*)directoryPCMPi08TeV->Get("Pi0_Input");
	TH1D*		histoPi0InputFullAddSigReweighted8TeV				= (TH1D*)directoryPCMPi08TeV->Get("Pi0_Input_Reweighted_AddedSig");
	TH1D*		histoPi0InputFullAddSig8TeV							= (TH1D*)directoryPCMPi08TeV->Get("Pi0_Input_AddedSig");

	TDirectory*	directoryPCMEta8TeV 								= (TDirectory*)filePCMpp8TeV->Get("Eta8TeV"); 
	TH1D* 		histoPCMYieldEta8TeV 								= (TH1D*)directoryPCMEta8TeV->Get(nameHistoPCMEta.Data());
	TH1D*		histoEtaInputFullReweighted8TeV						= (TH1D*)directoryPCMEta8TeV->Get("Eta_Input_Reweighted");
	TH1D*		histoEtaInputFull8TeV								= (TH1D*)directoryPCMEta8TeV->Get("Eta_Input");
	TH1D*		histoEtaInputFullAddSigReweighted8TeV				= (TH1D*)directoryPCMEta8TeV->Get("Eta_Input_Reweighted_AddedSig");
	TH1D*		histoEtaInputFullAddSig8TeV							= (TH1D*)directoryPCMEta8TeV->Get("Eta_Input_AddedSig");
	
	// Pythia 8
	TFile* 		filePi0PCM8TeV_LHC14e2a_PYT8 						= new TFile("0000011_002000092273028250400000_01521035000000_LHC14e2a/8TeV/Pi0_data_GammaConvV1Correction_0000011_002000092273028250400000_01521035000000.root");
	TH1D* 		histoPi0InputMCWOWeights_LHC14e2a_PYT8_8TeV 		= (TH1D*)filePi0PCM8TeV_LHC14e2a_PYT8->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoPi0Efficiency_LHC14e2a_PYT8_8TeV 				= (TH1D*)filePi0PCM8TeV_LHC14e2a_PYT8->Get("TrueMesonEffiPt");
	TFile* 		fileEtaPCM8TeV_LHC14e2a_PYT8 						= new TFile("0000011_002000092273028250400000_01521035000000_LHC14e2a/8TeV/Eta_data_GammaConvV1Correction_0000011_002000092273028250400000_01521035000000.root");
	TH1D* 		histoEtaInputMCWOWeights_LHC14e2a_PYT8_8TeV 		= (TH1D*)fileEtaPCM8TeV_LHC14e2a_PYT8->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoEtaEfficiency_LHC14e2a_PYT8_8TeV 				= (TH1D*)fileEtaPCM8TeV_LHC14e2a_PYT8->Get("TrueMesonEffiPt");

	// Phojet
	TFile* 		filePi0PCM8TeV_LHC14e2c_PHO 						= new TFile("0000011_002000092273028250400000_01521035000000_LHC14e2c/8TeV/Pi0_data_GammaConvV1Correction_0000011_002000092273028250400000_01521035000000.root");
	TH1D* 		histoPi0InputMCWOWeights_LHC14e2c_PHO_8TeV 			= (TH1D*)filePi0PCM8TeV_LHC14e2c_PHO->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoPi0Efficiency_LHC14e2c_PHO_8TeV 				= (TH1D*)filePi0PCM8TeV_LHC14e2c_PHO->Get("TrueMesonEffiPt");
	TFile*		fileEtaPCM8TeV_LHC14e2c_PHO 						= new TFile("0000011_002000092273028250400000_01521035000000_LHC14e2c/8TeV/Eta_data_GammaConvV1Correction_0000011_002000092273028250400000_01521035000000.root");
	TH1D*		histoEtaInputMCWOWeights_LHC14e2c_PHO_8TeV 			= (TH1D*)fileEtaPCM8TeV_LHC14e2c_PHO->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoEtaEfficiency_LHC14e2c_PHO_8TeV 				= (TH1D*)fileEtaPCM8TeV_LHC14e2c_PHO->Get("TrueMesonEffiPt");

	// Pythia 8, LHC12i3
	TFile*		filePi0PCM8TeV_LHC14e2b_PYT8 						= new TFile("0000011_002000092273028250400000_01521035000000_LHC14e2b/8TeV/Pi0_data_GammaConvV1Correction_0000011_002000092273028250400000_01521035000000.root");
	TH1D*		histoPi0InputMCWOWeights_LHC14e2b_PYT8_8TeV 		= (TH1D*)filePi0PCM8TeV_LHC14e2b_PYT8->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoPi0Efficiency_LHC14e2b_PYT8_8TeV 				= (TH1D*)filePi0PCM8TeV_LHC14e2b_PYT8->Get("TrueMesonEffiPt");
	TFile*		fileEtaPCM8TeV_LHC14e2b_PYT8 						= new TFile("0000011_002000092273028250400000_01521035000000_LHC14e2b/8TeV/Eta_data_GammaConvV1Correction_0000011_002000092273028250400000_01521035000000.root");
	TH1D*		histoEtaInputMCWOWeights_LHC14e2b_PYT8_8TeV 		= (TH1D*)fileEtaPCM8TeV_LHC14e2b_PYT8->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoEtaEfficiency_LHC14e2b_PYT8_8TeV 				= (TH1D*)fileEtaPCM8TeV_LHC14e2b_PYT8->Get("TrueMesonEffiPt");

	// Pythia 8, LHC12i3 added signals
	TFile*		filePi0PCM8TeV_LHC14e2b_PYT8Addsig 					= new TFile("0000012_002000092273028250400000_01521035000000_LHC14e2b/8TeV/Pi0_MC_GammaConvV1Correction_0000012_002000092273028250400000_01521035000000.root");
	TH1D*		histoPi0InputMCWOWeights_LHC14e2b_PYT8addSig_8TeV 	= (TH1D*)filePi0PCM8TeV_LHC14e2b_PYT8Addsig->Get("MCYield_Meson_oldBinWOWeights");
	TH1D*		histoPi0Efficiency_LHC14e2b_PYT8addSig_8TeV 		= (TH1D*)filePi0PCM8TeV_LHC14e2b_PYT8Addsig->Get("TrueMesonEffiPt");
	TFile*		fileEtaPCM8TeV_LHC14e2b_PYT8Addsig 					= new TFile("0000012_002000092273028250400000_01521035000000_LHC14e2b/8TeV/Eta_MC_GammaConvV1Correction_0000012_002000092273028250400000_01521035000000.root");
	TH1D* 		histoEtaInputMCWOWeights_LHC14e2b_PYT8addSig_8TeV 	= (TH1D*)fileEtaPCM8TeV_LHC14e2b_PYT8Addsig->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* 		histoEtaEfficiency_LHC14e2b_PYT8addSig_8TeV 		= (TH1D*)fileEtaPCM8TeV_LHC14e2b_PYT8Addsig->Get("TrueMesonEffiPt");

	
	//	**********************************************************************************************************************
	//	****************************************Pi0 Spectra 2.76TeV compared to MC********************************************
	//	**********************************************************************************************************************
	TCanvas* canvasPi0Spectra2760GeV 								= new TCanvas("canvasPi0Spectra2760GeV", "", 200, 10, 1200, 1100);  // gives the page size
	DrawGammaCanvasSettings( canvasPi0Spectra2760GeV,  0.13, 0.01, 0.015, 0.08);
	canvasPi0Spectra2760GeV->SetLogy();
	canvasPi0Spectra2760GeV->SetLogx();
	
	TH2F * histo2DPi0Spectra2760GeV;
	histo2DPi0Spectra2760GeV 											= new TH2F("histo2DPi0Spectra2760GeV", "histo2DPi0Spectra2760GeV",1000, 0.23, 20., 1000, 1e-8, 1e1 );
	SetStyleHistoTH2ForGraphs( histo2DPi0Spectra2760GeV, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
							   0.03, 0.04, 0.03, 0.04, 0.83, 1.4);
	histo2DPi0Spectra2760GeV->GetXaxis()->SetLabelOffset(-0.01);
	histo2DPi0Spectra2760GeV->GetYaxis()->SetLabelOffset(0.01);
	histo2DPi0Spectra2760GeV->DrawCopy(); 
	
	DrawGammaSetMarker(histoPCMYieldPi02760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
	histoPCMYieldPi02760GeV->Draw("p,same,e1");
	SetStyleHisto(histoPi0InputFullReweighted2760GeV, 1., lineStyleMCA, kBlue+1);  
	histoPi0InputFullReweighted2760GeV->Draw("same,hist,c");    
	SetStyleHisto(histoPi0InputFullAddSigReweighted2760GeV, 1., lineStyleMCAddSig, kBlue+2);  
	histoPi0InputFullAddSigReweighted2760GeV->Draw("same,hist,c");    

	SetStyleHisto(histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV, 1., lineStyleMCA, colorMCPythiaPP2760GeV);  
	histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV->Draw("same,hist,c");    
	SetStyleHisto(histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV, 1., lineStyleMCB, colorMCPhojetPP2760GeV);  
	histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV->Draw("same,hist,c");    
	SetStyleHisto(histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV, 1., lineStyleMCC, colorMCPythia2PP2760GeV);  
	histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV->Draw("same,hist,c");    
	SetStyleHisto(histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV, 1., lineStyleMCAddSig, colorMCPythiaAddSigPP2760GeV);  
	histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV->Draw("same,hist,c");    

	TLatex *labelSpectraPi0Label 									= new TLatex(0.55,0.92,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelSpectraPi0Label, 0.035  ,4);
	labelSpectraPi0Label->Draw();
	
	TLegend* legendSpectraPi02760GeV 										= new TLegend(0.16,0.09,0.73,0.3);
	legendSpectraPi02760GeV->SetFillColor(0);
	legendSpectraPi02760GeV->SetLineColor(0);
	legendSpectraPi02760GeV->SetTextSize(0.025);
// 	legendSpectraPi02760GeV->SetNColumns(2);
	legendSpectraPi02760GeV->SetMargin(0.2);
	legendSpectraPi02760GeV->AddEntry(histoPCMYieldPi02760GeV,collisionSystemPP2760GeV.Data(),"pf");
	legendSpectraPi02760GeV->AddEntry(histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV,"Phojet, LHC12f1b","l");
	legendSpectraPi02760GeV->AddEntry(histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV,"Pythia 8, LHC12f1a","l");
	legendSpectraPi02760GeV->AddEntry(histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV,"Pythia 8, LHC12i3","l");
	legendSpectraPi02760GeV->AddEntry(histoPi0InputFullReweighted2760GeV,"full MC reweighted","l");
	legendSpectraPi02760GeV->AddEntry(histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV,"Pythia 8, LHC12i3, added Signals","l");
	legendSpectraPi02760GeV->AddEntry(histoPi0InputFullAddSigReweighted2760GeV,"added Signals, reweighted","l");
	legendSpectraPi02760GeV->Draw();
	
	canvasPi0Spectra2760GeV->Update();
	canvasPi0Spectra2760GeV->Print(Form("%s/Pi0_Spectra_2760GeV.%s",outputDir.Data(),suffix.Data()));

	//	**********************************************************************************************************************
	//	****************************************Fit Pi0 Spectra 2.76TeV and plot *********************************************
	//	**********************************************************************************************************************
	
	// fit spectrum with Tsallis function
	TF1* fitInvYieldDataPi0Comb2760GeV 								= FitObject("l","fitInvYieldDataPi0Comb2760GeV","Pi0");
	histoPCMYieldPi02760GeV->Fit(fitInvYieldDataPi0Comb2760GeV,"QNRMEI+","",0.4,10.);
	fitInvYieldDataPi0Comb2760GeV->SetRange(0.1,20.);
	// print fit result to shell
	cout << WriteParameterToFile(fitInvYieldDataPi0Comb2760GeV)<< endl;	

	// plot result with input
	canvasPi0Spectra2760GeV->cd();
	histo2DPi0Spectra2760GeV->DrawCopy(); 

	DrawGammaSetMarker(histoPCMYieldPi02760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
	histoPCMYieldPi02760GeV->Draw("p,same,e1");

	
	fitInvYieldDataPi0Comb2760GeV->SetLineColor(color2760GeV);
	fitInvYieldDataPi0Comb2760GeV->Draw("same");
	
	SetStyleHisto(histoPi0InputFullReweighted2760GeV, 1., lineStyleMCAddSig, kBlue+1);  
	histoPi0InputFullReweighted2760GeV->Draw("same,hist,c");    
	SetStyleHisto(histoPi0InputFullAddSigReweighted2760GeV, 1., lineStyleMCAddSig, kBlue+2);  
	histoPi0InputFullAddSigReweighted2760GeV->Draw("same,hist,c");    

	labelSpectraPi0Label->Draw();
	
	canvasPi0Spectra2760GeV->Print(Form("%s/Pi0_Spectra_WithFit_2760GeV.%s",outputDir.Data(),suffix.Data()));
	
	// **********************************************************************************************************************
    // ******************************* Ratio of data to fit and MC input to fit for pi0 2.76 TeV ****************************
	// **********************************************************************************************************************
	// definition of canvas for ratios
	TCanvas* canvasRatioToFit = new TCanvas("canvasRatioToFit","",1550,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioToFit,  0.08, 0.015, 0.015, 0.08);
	canvasRatioToFit->SetGridx(0);
	canvasRatioToFit->SetGridy(0);

	// Calculation of ratio histograms
	TH1D* histoRatioPi0DatatoFit2760GeV 		= CalculateHistoRatioToFit (histoPCMYieldPi02760GeV, fitInvYieldDataPi0Comb2760GeV,kTRUE);
	TH1D* histoRatioPi0MCtoDataFit2760GeV 		= CalculateHistoRatioToFit (histoPi0InputFullReweighted2760GeV, fitInvYieldDataPi0Comb2760GeV,kTRUE);
	TH1D* histoRatioPi0MCAddSigtoDataFit2760GeV = CalculateHistoRatioToFit (histoPi0InputFullAddSigReweighted2760GeV, fitInvYieldDataPi0Comb2760GeV,kTRUE);
	TH1D* histoRatioPi0MCUnweightedtoDataFit2760GeV = NULL;
	if (histoPi0InputFull2760GeV) histoRatioPi0MCUnweightedtoDataFit2760GeV = CalculateHistoRatioToFit (histoPi0InputFull2760GeV, fitInvYieldDataPi0Comb2760GeV,kTRUE);
	if (histoRatioPi0MCUnweightedtoDataFit2760GeV) SetStyleHisto(histoRatioPi0MCUnweightedtoDataFit2760GeV, 2, lineStyleMCB, 807 );
	canvasRatioToFit->cd();
	SetStyleHisto(histoRatioPi0MCtoDataFit2760GeV, 2, lineStyleMCA, kRed+2 );
	SetStyleHisto(histoRatioPi0MCAddSigtoDataFit2760GeV, 2, lineStyleMCAddSig, kBlue+2 );
	DrawGammaSetMarker(histoRatioPi0DatatoFit2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, kBlack , kBlack);
	DrawAutoGammaMesonHistos( histoRatioPi0DatatoFit2760GeV,
				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
				kFALSE, 1.5, 0, kTRUE,
				kTRUE, 0, 2.1,
				kTRUE, 0.,histoRatioPi0DatatoFit2760GeV->GetXaxis()->GetBinUpEdge(histoRatioPi0DatatoFit2760GeV->GetNbinsX()));
	histoRatioPi0DatatoFit2760GeV->GetYaxis()->SetTitleOffset(0.9);
	histoRatioPi0DatatoFit2760GeV->Draw("e,p");  
	if (runDrawReweighted) histoRatioPi0MCtoDataFit2760GeV->Draw("same,hist,l");  
	if (runDrawReweighted) histoRatioPi0MCAddSigtoDataFit2760GeV->Draw("same,hist,l");  
	if (histoRatioPi0MCUnweightedtoDataFit2760GeV) histoRatioPi0MCUnweightedtoDataFit2760GeV->Draw("same,hist,l");  

	TLatex *labelSpectraPi0LabelRatio 									= new TLatex(0.65,0.9,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelSpectraPi0LabelRatio, 0.035  ,4);
	labelSpectraPi0LabelRatio->Draw();
	TLatex *labelEnergy2760GeVRatio 									= new TLatex(0.65,0.86,collisionSystemPP2760GeV.Data());
	SetStyleTLatex( labelEnergy2760GeVRatio, 0.035  ,4);
	labelEnergy2760GeVRatio->Draw();
	
	TLegend* legendRatioPi02760GeV = new TLegend(0.11,0.12,0.4,0.30);
	legendRatioPi02760GeV->SetFillColor(0);
	legendRatioPi02760GeV->SetLineColor(0);
	legendRatioPi02760GeV->SetTextSize(0.035);
	legendRatioPi02760GeV->SetMargin(0.2);
	legendRatioPi02760GeV->AddEntry(histoRatioPi0DatatoFit2760GeV,"Data/Tsallis fit to Data (0.4 <#it{p}_{T}<10)","p");
	if (runDrawReweighted) legendRatioPi02760GeV->AddEntry(histoRatioPi0MCtoDataFit2760GeV,Form("MC weighted %s/Tsallis  fit to Data (0.4 <#it{p}_{T}<10)",stringIterationNumber.Data()),"l");
	if (runDrawReweighted) legendRatioPi02760GeV->AddEntry(histoRatioPi0MCAddSigtoDataFit2760GeV,Form("MC add Sig weighted %s/Tsallis  fit to Data (0.4 <#it{p}_{T}<10)", 
														   stringIterationNumber.Data()), "l");
	if (histoRatioPi0MCUnweightedtoDataFit2760GeV) legendRatioPi02760GeV->AddEntry(histoRatioPi0MCUnweightedtoDataFit2760GeV,"MC/Tsallis  fit to Data (0.4 <#it{p}_{T}<10)","l");
	legendRatioPi02760GeV->Draw();
	
	DrawGammaLines(0., histoRatioPi0DatatoFit2760GeV->GetXaxis()->GetBinUpEdge(histoRatioPi0DatatoFit2760GeV->GetNbinsX()) ,1., 1.,0.1);
	canvasRatioToFit->Update();
	canvasRatioToFit->SaveAs(Form("%s/Pi0_RatioToDataFit_2760GeV.%s",outputDir.Data(),suffix.Data()));
	
	//	**********************************************************************************************************************
	//	******************************Compare Efficiencies for pi0 at 2.76TeV for different MC *******************************
	//	**********************************************************************************************************************
	TCanvas* canvasPi0Efficiencies2760GeVDiffMC 					= new TCanvas("canvasPi0Efficiencies2760GeVDiffMC", "", 200, 10, 1200, 1100);  // gives the page size
	DrawGammaCanvasSettings( canvasPi0Efficiencies2760GeVDiffMC,  0.09, 0.01, 0.015, 0.08);
	canvasPi0Efficiencies2760GeVDiffMC->SetLogy();
// 	canvasPi0Efficiencies2760GeVDiffMC->SetLogx();
	
	TH2F * histo2DPi0Effi2760GeV;
	histo2DPi0Effi2760GeV 													= new TH2F("histo2DPi0Effi2760GeV", "histo2DPi0Effi2760GeV",1000, 0., 10.5, 1000, 1e-5, 1e-2 );
	SetStyleHistoTH2ForGraphs( histo2DPi0Effi2760GeV, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{#pi^{0}}", 
							   0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
	histo2DPi0Effi2760GeV->GetYaxis()->SetLabelOffset(0.01);
	histo2DPi0Effi2760GeV->DrawCopy(); 

	DrawGammaSetMarker(histoPi0Efficiency_LHC12f1a_PYT8_2760GeV, markerStyleMCA, markerSizePP2760GeV, colorMCPythiaPP2760GeV , colorMCPythiaPP2760GeV);
	histoPi0Efficiency_LHC12f1a_PYT8_2760GeV->Draw("p,same,e1");

	DrawGammaSetMarker(histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV, markerStyleMCAddSig, markerSizePP2760GeV, colorMCPythiaAddSigPP2760GeV , colorMCPythiaAddSigPP2760GeV);
	histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV->Draw("p,same,e1");

	DrawGammaSetMarker(histoPi0Efficiency_LHC12i3_PYT8_2760GeV, markerStyleMCC, markerSizePP2760GeV, colorMCPythia2PP2760GeV , colorMCPythia2PP2760GeV);
	histoPi0Efficiency_LHC12i3_PYT8_2760GeV->Draw("p,same,e1");
	
	DrawGammaSetMarker(histoPi0Efficiency_LHC12f1b_PHO_2760GeV, markerStyleMCB, markerSizePP2760GeV, colorMCPhojetPP2760GeV , colorMCPhojetPP2760GeV);
	histoPi0Efficiency_LHC12f1b_PHO_2760GeV->Draw("p,same,e1");

	TLegend* legendEfficiencyPi02760GeV 										= new TLegend(0.56,0.09,0.93,0.3);
	legendEfficiencyPi02760GeV->SetFillColor(0);
	legendEfficiencyPi02760GeV->SetLineColor(0);
	legendEfficiencyPi02760GeV->SetTextSize(0.025);
	legendEfficiencyPi02760GeV->SetMargin(0.2);
	legendEfficiencyPi02760GeV->AddEntry(histoPi0Efficiency_LHC12f1b_PHO_2760GeV,"Phojet, LHC12f1b","p");
	legendEfficiencyPi02760GeV->AddEntry(histoPi0Efficiency_LHC12f1a_PYT8_2760GeV,"Pythia 8, LHC12f1a","p");
	legendEfficiencyPi02760GeV->AddEntry(histoPi0Efficiency_LHC12i3_PYT8_2760GeV,"Pythia 8, LHC12i3","p");
	legendEfficiencyPi02760GeV->AddEntry(histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV,"Pythia 8, LHC12i3, added Signals","p");
	legendEfficiencyPi02760GeV->Draw();
	
	canvasPi0Efficiencies2760GeVDiffMC->Update();
	canvasPi0Efficiencies2760GeVDiffMC->Print(Form("%s/Pi0_Efficiency_WOSDD_2760GeV.%s",outputDir.Data(),suffix.Data()));

	canvasPi0Efficiencies2760GeVDiffMC->cd();
	histo2DPi0Effi2760GeV->DrawCopy(); 

	DrawGammaSetMarker(histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDD, markerStyleMCA, markerSizePP2760GeV, colorMCPythiaPP2760GeV , colorMCPythiaPP2760GeV);
	histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDD->Draw("p,same,e1");

	DrawGammaSetMarker(histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDD, markerStyleMCAddSig, markerSizePP2760GeV, colorMCPythiaAddSigPP2760GeV , colorMCPythiaAddSigPP2760GeV);
	histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDD->Draw("p,same,e1");

	DrawGammaSetMarker(histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDD, markerStyleMCC, markerSizePP2760GeV, colorMCPythia2PP2760GeV , colorMCPythia2PP2760GeV);
	histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDD->Draw("p,same,e1");
	
	DrawGammaSetMarker(histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDD, markerStyleMCB, markerSizePP2760GeV, colorMCPhojetPP2760GeV , colorMCPhojetPP2760GeV);
	histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDD->Draw("p,same,e1");

	legendEfficiencyPi02760GeV->Draw();
	
	canvasPi0Efficiencies2760GeVDiffMC->Update();
	canvasPi0Efficiencies2760GeVDiffMC->Print(Form("%s/Pi0_Efficiency_WSDD_2760GeV.%s",outputDir.Data(),suffix.Data()));
	
	canvasPi0Efficiencies2760GeVDiffMC->SetLogy(0);
	canvasPi0Efficiencies2760GeVDiffMC->SetTopMargin(0.035);
	
	histo2DPi0Effi2760GeV->GetYaxis()->SetRangeUser(1e-5,3E-3);
	histo2DPi0Effi2760GeV->DrawCopy(); 
	histoPi0Efficiency_LHC12f1a_PYT8_2760GeV->Draw("p,same,e1");
	histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV->Draw("p,same,e1");
	histoPi0Efficiency_LHC12i3_PYT8_2760GeV->Draw("p,same,e1");
	histoPi0Efficiency_LHC12f1b_PHO_2760GeV->Draw("p,same,e1");
	legendEfficiencyPi02760GeV->Draw();
	
	canvasPi0Efficiencies2760GeVDiffMC->Update();
	canvasPi0Efficiencies2760GeVDiffMC->Print(Form("%s/Pi0_Efficiency_WOSDD_2760GeV_LinY.%s",outputDir.Data(),suffix.Data()));
	
	canvasPi0Efficiencies2760GeVDiffMC->cd();
	
	histo2DPi0Effi2760GeV->DrawCopy(); 
	histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDD->Draw("p,same,e1");
	histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDD->Draw("p,same,e1");
	histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDD->Draw("p,same,e1");
	histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDD->Draw("p,same,e1");
	legendEfficiencyPi02760GeV->Draw();
	
	canvasPi0Efficiencies2760GeVDiffMC->Update();
	canvasPi0Efficiencies2760GeVDiffMC->Print(Form("%s/Pi0_Efficiency_WSDD_2760GeV_LinY.%s",outputDir.Data(),suffix.Data()));
	
	//	**********************************************************************************************************************
	//	****************************************Eta Spectra 2.76TeV compared to MC********************************************
	//	**********************************************************************************************************************	
	TCanvas* canvasEtaSpectra2760GeV 								= new TCanvas("canvasEtaSpectra2760GeV", "", 200, 10, 1200, 1100);  // gives the page size
	DrawGammaCanvasSettings( canvasEtaSpectra2760GeV,  0.13, 0.01, 0.015, 0.08);
	canvasEtaSpectra2760GeV->SetLogy();
	canvasEtaSpectra2760GeV->SetLogx();
	
	TH2F * histo2DEtaSpectraAll;
	histo2DEtaSpectraAll 											= new TH2F("histo2DEtaSpectraAll", "histo2DEtaSpectraAll",1000, 0.23, 20., 1000, 1e-8, 1e1 );
	SetStyleHistoTH2ForGraphs( histo2DEtaSpectraAll, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
							   0.03, 0.04, 0.03, 0.04, 0.83, 1.4);
	histo2DEtaSpectraAll->GetXaxis()->SetLabelOffset(-0.01);
	histo2DEtaSpectraAll->GetYaxis()->SetLabelOffset(0.01);
	histo2DEtaSpectraAll->DrawCopy(); 
	
	DrawGammaSetMarker(histoPCMYieldEta2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
	histoPCMYieldEta2760GeV->Draw("p,same,e1");
	SetStyleHisto(histoEtaInputFullReweighted2760GeV, 1., lineStyleMCA, kBlue+1);  
	histoEtaInputFullReweighted2760GeV->Draw("same,hist,c");    
	SetStyleHisto(histoEtaInputFullAddSigReweighted2760GeV, 1., lineStyleMCAddSig, kBlue+2);  
	histoEtaInputFullAddSigReweighted2760GeV->Draw("same,hist,c");    

	SetStyleHisto(histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV, 1., lineStyleMCA, colorMCPythiaPP2760GeV);  
	histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV->Draw("same,hist,c");    
	SetStyleHisto(histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV, 1., lineStyleMCB, colorMCPhojetPP2760GeV);  
	histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV->Draw("same,hist,c");    
	SetStyleHisto(histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV, 1., lineStyleMCC, colorMCPythia2PP2760GeV);  
	histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV->Draw("same,hist,c");    
	SetStyleHisto(histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV, 1., lineStyleMCAddSig, colorMCPythiaAddSigPP2760GeV);  
	histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV->Draw("same,hist,c");    
	
	TLatex *labelSpectraEtaLabel 									= new TLatex(0.55,0.92,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelSpectraEtaLabel, 0.035  ,4);
	labelSpectraEtaLabel->Draw();
	
	TLegend* legendSpectraEta2760GeV 										= new TLegend(0.16,0.09,0.73,0.3);
	legendSpectraEta2760GeV->SetFillColor(0);
	legendSpectraEta2760GeV->SetLineColor(0);
	legendSpectraEta2760GeV->SetTextSize(0.025);
	legendSpectraEta2760GeV->SetMargin(0.2);
	legendSpectraEta2760GeV->AddEntry(histoPCMYieldEta2760GeV,collisionSystemPP2760GeV.Data(),"pf");
	legendSpectraEta2760GeV->AddEntry(histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV,"Phojet, LHC12f1b","l");
	legendSpectraEta2760GeV->AddEntry(histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV,"Pythia 8, LHC12f1a","l");
	legendSpectraEta2760GeV->AddEntry(histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV,"Pythia 8, LHC12i3","l");
	legendSpectraEta2760GeV->AddEntry(histoEtaInputFullReweighted2760GeV,"full MC reweighted","l");
	legendSpectraEta2760GeV->AddEntry(histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV,"Pythia 8, LHC12i3, added Signals","l");
	legendSpectraEta2760GeV->AddEntry(histoEtaInputFullAddSigReweighted2760GeV,"added Signals, reweighted","l");
	legendSpectraEta2760GeV->Draw();
	
	canvasEtaSpectra2760GeV->Update();
	canvasEtaSpectra2760GeV->Print(Form("%s/Eta_Spectra_All.%s",outputDir.Data(),suffix.Data()));

	//	**********************************************************************************************************************
	//	****************************************Fit Pi0 Spectra 2.76TeV and plot *********************************************
	//	**********************************************************************************************************************
	
	// fit spectrum with Tsallis function - "n" fixed by result from pi0 at same energy	
	TF1* fitInvYieldDataEtaComb2760GeV 								= FitObject("l","fitInvYieldDataEtaComb2760GeV","Eta");
// 	fitInvYieldDataEtaComb2760GeV->FixParameter(1,fitInvYieldDataPi0Comb2760GeV->GetParameter(1));
// 	fitInvYieldDataEtaComb2760GeV->SetParameter(1,fitInvYieldDataPi0Comb2760GeV->GetParameter(1));
// 	fitInvYieldDataEtaComb2760GeV->SetParLimits(1,fitInvYieldDataPi0Comb2760GeV->GetParameter(1)*0.9,fitInvYieldDataPi0Comb2760GeV->GetParameter(1)*1.1);
	histoPCMYieldEta2760GeV->Fit(fitInvYieldDataEtaComb2760GeV,"QNRMEI+","",0.6,6.);
	fitInvYieldDataEtaComb2760GeV->SetRange(0.1,20.);
	// print result to shell
	cout << WriteParameterToFile(fitInvYieldDataEtaComb2760GeV)<< endl;	
	
	// plot result with input
	canvasEtaSpectra2760GeV->cd();
	histo2DEtaSpectraAll->DrawCopy(); 

	DrawGammaSetMarker(histoPCMYieldEta2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
	histoPCMYieldEta2760GeV->Draw("p,same,e1");

	fitInvYieldDataEtaComb2760GeV->SetLineColor(color2760GeV);
	fitInvYieldDataEtaComb2760GeV->Draw("same");
	SetStyleHisto(histoEtaInputFullReweighted2760GeV, 1., lineStyleMCAddSig, kBlue+1);  
	histoEtaInputFullReweighted2760GeV->Draw("same,hist,c");    
	SetStyleHisto(histoEtaInputFullAddSigReweighted2760GeV, 1., lineStyleMCAddSig, kBlue+2);  
	histoEtaInputFullAddSigReweighted2760GeV->Draw("same,hist,c");    

	
	labelSpectraEtaLabel->Draw();
	
	canvasEtaSpectra2760GeV->Print(Form("%s/Eta_Spectra_WithFit_2760GeV.%s",outputDir.Data(),suffix.Data()));

	// **********************************************************************************************************************
    // ******************************* Ratio of data to fit and MC input to fit for Eta 2.76 TeV ****************************
	// **********************************************************************************************************************
	
	// Calculation of ratio histograms
	TH1D* histoRatioEtaDatatoFit2760GeV 		= CalculateHistoRatioToFit (histoPCMYieldEta2760GeV, fitInvYieldDataEtaComb2760GeV,kTRUE);
	TH1D* histoRatioEtaMCtoDataFit2760GeV 		= CalculateHistoRatioToFit (histoEtaInputFullReweighted2760GeV, fitInvYieldDataEtaComb2760GeV,kTRUE);
	TH1D* histoRatioEtaMCAddSigtoDataFit2760GeV = CalculateHistoRatioToFit (histoEtaInputFullAddSigReweighted2760GeV, fitInvYieldDataEtaComb2760GeV,kTRUE);
	TH1D* histoRatioEtaMCUnweightedtoDataFit2760GeV = NULL;
	if (histoEtaInputFull2760GeV) histoRatioEtaMCUnweightedtoDataFit2760GeV = CalculateHistoRatioToFit (histoEtaInputFull2760GeV, fitInvYieldDataEtaComb2760GeV,kTRUE);
	if (histoRatioEtaMCUnweightedtoDataFit2760GeV) SetStyleHisto(histoRatioEtaMCUnweightedtoDataFit2760GeV, 2, lineStyleMCB, 807 );
	// plotting ratios
	canvasRatioToFit->cd();
	SetStyleHisto(histoRatioEtaMCtoDataFit2760GeV, 2, lineStyleMCA, kRed+2 );
	SetStyleHisto(histoRatioEtaMCAddSigtoDataFit2760GeV, 2, lineStyleMCAddSig, kBlue+2 );
	DrawGammaSetMarker(histoRatioEtaDatatoFit2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, kBlack , kBlack);
	DrawAutoGammaMesonHistos( histoRatioEtaDatatoFit2760GeV,
				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
				kFALSE, 1.5, 0, kTRUE,
				kTRUE, 0, 2.1,
				kTRUE, 0.,histoRatioEtaDatatoFit2760GeV->GetXaxis()->GetBinUpEdge(histoRatioEtaDatatoFit2760GeV->GetNbinsX()));
	histoRatioEtaDatatoFit2760GeV->GetYaxis()->SetTitleOffset(0.9);
	histoRatioEtaDatatoFit2760GeV->Draw("e,p");  
	if (runDrawReweighted) histoRatioEtaMCtoDataFit2760GeV->Draw("same,hist,l");  
	if (runDrawReweighted) histoRatioEtaMCAddSigtoDataFit2760GeV->Draw("same,hist,l");  
	if (histoRatioEtaMCUnweightedtoDataFit2760GeV) histoRatioEtaMCUnweightedtoDataFit2760GeV->Draw("same,hist,l");  

	// labeling
	TLatex *labelSpectraEtaLabelRatio 									= new TLatex(0.65,0.9,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelSpectraEtaLabelRatio, 0.035  ,4);
	labelSpectraEtaLabelRatio->Draw();
	labelEnergy2760GeVRatio->Draw();
	
	TLegend* legendRatioEta2760GeV = new TLegend(0.11,0.12,0.4,0.30);
	legendRatioEta2760GeV->SetFillStyle(0);
	legendRatioEta2760GeV->SetFillColor(0);
	legendRatioEta2760GeV->SetLineColor(0);
	legendRatioEta2760GeV->SetTextSize(0.035);
	legendRatioEta2760GeV->SetMargin(0.2);
	legendRatioEta2760GeV->AddEntry(histoRatioEtaDatatoFit2760GeV,"Data/Tsallis fit to Data (0.6 <#it{p}_{T}<6)","p");
	if (runDrawReweighted) legendRatioEta2760GeV->AddEntry(histoRatioEtaMCtoDataFit2760GeV,Form("MC weighted %s/Tsallis  fit to Data (0.6 <#it{p}_{T}<6)",stringIterationNumber.Data()),"l");
	if (runDrawReweighted) legendRatioEta2760GeV->AddEntry(histoRatioEtaMCAddSigtoDataFit2760GeV,Form("MC add Sig weighted %s/Tsallis  fit to Data (0.6 <#it{p}_{T}<6)",
														   stringIterationNumber.Data()), "l");
	if (histoRatioEtaMCUnweightedtoDataFit2760GeV) legendRatioEta2760GeV->AddEntry(histoRatioEtaMCUnweightedtoDataFit2760GeV,"MC/Tsallis  fit to Data (0.6 <#it{p}_{T}<6)","l");
	legendRatioEta2760GeV->Draw();
	
	DrawGammaLines(0., histoRatioEtaDatatoFit2760GeV->GetXaxis()->GetBinUpEdge(histoRatioEtaDatatoFit2760GeV->GetNbinsX()) ,1., 1.,0.1);
	canvasRatioToFit->Update();
	canvasRatioToFit->SaveAs(Form("%s/Eta_RatioToDataFit_2760GeV.%s",outputDir.Data(),suffix.Data()));
	
	
	//	**********************************************************************************************************************
	//	******************************Compare Efficiencies for Eta at 2.76TeV for different MC *******************************
	//	**********************************************************************************************************************
	TCanvas* canvasEtaEfficiencies2760GeVDiffMC 					= new TCanvas("canvasEtaEfficiencies2760GeVDiffMC", "", 200, 10, 1200, 1100);  // gives the page size
	DrawGammaCanvasSettings( canvasEtaEfficiencies2760GeVDiffMC,  0.09, 0.01, 0.015, 0.08);
	canvasEtaEfficiencies2760GeVDiffMC->SetLogy();
// 	canvasEtaEfficiencies2760GeVDiffMC->SetLogx();
	
	TH2F * histo2DEtaEffi2760GeV;
	histo2DEtaEffi2760GeV 													= new TH2F("histo2DEtaEffi2760GeV", "histo2DEtaEffi2760GeV",1000, 0., 6.5, 1000, 1e-5, 1e-2 );
	SetStyleHistoTH2ForGraphs( histo2DEtaEffi2760GeV, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{#eta}", 
							   0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
	histo2DEtaEffi2760GeV->GetYaxis()->SetLabelOffset(0.01);
	histo2DEtaEffi2760GeV->DrawCopy(); 

	DrawGammaSetMarker(histoEtaEfficiency_LHC12f1a_PYT8_2760GeV, markerStyleMCA, markerSizePP2760GeV, colorMCPythiaPP2760GeV , colorMCPythiaPP2760GeV);
	histoEtaEfficiency_LHC12f1a_PYT8_2760GeV->Draw("p,same,e1");

	DrawGammaSetMarker(histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV, markerStyleMCAddSig, markerSizePP2760GeV, colorMCPythiaAddSigPP2760GeV , colorMCPythiaAddSigPP2760GeV);
	histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV->Draw("p,same,e1");

	DrawGammaSetMarker(histoEtaEfficiency_LHC12i3_PYT8_2760GeV, markerStyleMCC, markerSizePP2760GeV, colorMCPythia2PP2760GeV , colorMCPythia2PP2760GeV);
	histoEtaEfficiency_LHC12i3_PYT8_2760GeV->Draw("p,same,e1");
	
	DrawGammaSetMarker(histoEtaEfficiency_LHC12f1b_PHO_2760GeV, markerStyleMCB, markerSizePP2760GeV, colorMCPhojetPP2760GeV , colorMCPhojetPP2760GeV);
	histoEtaEfficiency_LHC12f1b_PHO_2760GeV->Draw("p,same,e1");

	TLegend* legendEfficiencyEta2760GeV 										= new TLegend(0.56,0.09,0.93,0.3);
	legendEfficiencyEta2760GeV->SetFillColor(0);
	legendEfficiencyEta2760GeV->SetLineColor(0);
	legendEfficiencyEta2760GeV->SetTextSize(0.025);
	legendEfficiencyEta2760GeV->SetMargin(0.2);
	legendEfficiencyEta2760GeV->AddEntry(histoEtaEfficiency_LHC12f1b_PHO_2760GeV,"Phojet, LHC12f1b","p");
	legendEfficiencyEta2760GeV->AddEntry(histoEtaEfficiency_LHC12f1a_PYT8_2760GeV,"Pythia 8, LHC12f1a","p");
	legendEfficiencyEta2760GeV->AddEntry(histoEtaEfficiency_LHC12i3_PYT8_2760GeV,"Pythia 8, LHC12i3","p");
	legendEfficiencyEta2760GeV->AddEntry(histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV,"Pythia 8, LHC12i3, added Signals","p");
	legendEfficiencyEta2760GeV->Draw();
	
	canvasEtaEfficiencies2760GeVDiffMC->Update();
	canvasEtaEfficiencies2760GeVDiffMC->Print(Form("%s/Eta_Efficiency_WOSDD_2760GeV.%s",outputDir.Data(),suffix.Data()));

	canvasEtaEfficiencies2760GeVDiffMC->cd();
	histo2DEtaEffi2760GeV->DrawCopy(); 

	DrawGammaSetMarker(histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDD, markerStyleMCA, markerSizePP2760GeV, colorMCPythiaPP2760GeV , colorMCPythiaPP2760GeV);
	histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDD->Draw("p,same,e1");

	DrawGammaSetMarker(histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDD, markerStyleMCAddSig, markerSizePP2760GeV, colorMCPythiaAddSigPP2760GeV , colorMCPythiaAddSigPP2760GeV);
	histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDD->Draw("p,same,e1");

	DrawGammaSetMarker(histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDD, markerStyleMCC, markerSizePP2760GeV, colorMCPythia2PP2760GeV , colorMCPythia2PP2760GeV);
	histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDD->Draw("p,same,e1");
	
	DrawGammaSetMarker(histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDD, markerStyleMCB, markerSizePP2760GeV, colorMCPhojetPP2760GeV , colorMCPhojetPP2760GeV);
	histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDD->Draw("p,same,e1");

	legendEfficiencyEta2760GeV->Draw();
	
	canvasEtaEfficiencies2760GeVDiffMC->Update();
	canvasEtaEfficiencies2760GeVDiffMC->Print(Form("%s/Eta_Efficiency_WSDD_2760GeV.%s",outputDir.Data(),suffix.Data()));
	
	canvasEtaEfficiencies2760GeVDiffMC->SetLogy(0);
	canvasEtaEfficiencies2760GeVDiffMC->SetTopMargin(0.035);
	
	histo2DEtaEffi2760GeV->GetYaxis()->SetRangeUser(1e-5,3E-3);
	histo2DEtaEffi2760GeV->DrawCopy(); 
	histoEtaEfficiency_LHC12f1a_PYT8_2760GeV->Draw("p,same,e1");
	histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV->Draw("p,same,e1");
	histoEtaEfficiency_LHC12i3_PYT8_2760GeV->Draw("p,same,e1");
	histoEtaEfficiency_LHC12f1b_PHO_2760GeV->Draw("p,same,e1");
	legendEfficiencyEta2760GeV->Draw();
	
	canvasEtaEfficiencies2760GeVDiffMC->Update();
	canvasEtaEfficiencies2760GeVDiffMC->Print(Form("%s/Eta_Efficiency_WOSDD_2760GeV_LinY.%s",outputDir.Data(),suffix.Data()));

	canvasEtaEfficiencies2760GeVDiffMC->cd();
	
	histo2DEtaEffi2760GeV->DrawCopy(); 
	histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDD->Draw("p,same,e1");
	histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDD->Draw("p,same,e1");
	histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDD->Draw("p,same,e1");
	histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDD->Draw("p,same,e1");
	legendEfficiencyEta2760GeV->Draw();
	
	canvasEtaEfficiencies2760GeVDiffMC->Update();
	canvasEtaEfficiencies2760GeVDiffMC->Print(Form("%s/Eta_Efficiency_WSDD_2760GeV_LinY.%s",outputDir.Data(),suffix.Data()));
	
	//	**********************************************************************************************************************
	//	****************************************Pi0 Spectra 8TeV compared to MC********************************************
	//	**********************************************************************************************************************
	TCanvas* canvasPi0Spectra8TeV 								= new TCanvas("canvasPi0Spectra8TeV", "", 200, 10, 1200, 1100);  // gives the page size
	DrawGammaCanvasSettings( canvasPi0Spectra8TeV,  0.13, 0.01, 0.015, 0.08);
	canvasPi0Spectra8TeV->SetLogy();
	canvasPi0Spectra8TeV->SetLogx();
	
	TH2F * histo2DPi0Spectra8TeV;
	histo2DPi0Spectra8TeV 											= new TH2F("histo2DPi0Spectra8TeV", "histo2DPi0Spectra8TeV",1000, 0.23, 20., 1000, 1e-8, 1e1 );
	SetStyleHistoTH2ForGraphs( histo2DPi0Spectra8TeV, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
							   0.03, 0.04, 0.03, 0.04, 0.83, 1.4);
	histo2DPi0Spectra8TeV->GetXaxis()->SetLabelOffset(-0.01);
	histo2DPi0Spectra8TeV->GetYaxis()->SetLabelOffset(0.01);
	histo2DPi0Spectra8TeV->DrawCopy(); 
	
	DrawGammaSetMarker(histoPCMYieldPi08TeV, markerStyleSpectrum8TeV, markerSizePP8TeV, color8TeV , color8TeV);
	histoPCMYieldPi08TeV->Draw("p,same,e1");
	SetStyleHisto(histoPi0InputMCWOWeights_LHC14e2a_PYT8_8TeV, 1., lineStyleMCA, colorMCPythiaPP8TeV);  
	histoPi0InputMCWOWeights_LHC14e2a_PYT8_8TeV->Draw("same,hist,c");    
	SetStyleHisto(histoPi0InputMCWOWeights_LHC14e2c_PHO_8TeV, 1., lineStyleMCB, colorMCPhojetPP8TeV);  
	histoPi0InputMCWOWeights_LHC14e2c_PHO_8TeV->Draw("same,hist,c");    
	SetStyleHisto(histoPi0InputMCWOWeights_LHC14e2b_PYT8_8TeV, 1., lineStyleMCC, colorMCPythia2PP8TeV);  
	histoPi0InputMCWOWeights_LHC14e2b_PYT8_8TeV->Draw("same,hist,c");    
	SetStyleHisto(histoPi0InputMCWOWeights_LHC14e2b_PYT8addSig_8TeV, 1., lineStyleMCAddSig, colorMCPythiaAddSigPP8TeV);  
	histoPi0InputMCWOWeights_LHC14e2b_PYT8addSig_8TeV->Draw("same,hist,c");    
	SetStyleHisto(histoPi0InputFullReweighted8TeV, 1., lineStyleMCAddSig, kBlue);  
	histoPi0InputFullReweighted8TeV->Draw("same,hist,c");
	labelSpectraPi0Label->Draw();
	
	TLegend* legendSpectraPi08TeV 										= new TLegend(0.16,0.09,0.73,0.32);
	legendSpectraPi08TeV->SetFillColor(0);
	legendSpectraPi08TeV->SetLineColor(0);
	legendSpectraPi08TeV->SetTextSize(0.025);
// 	legendSpectraPi08TeV->SetNColumns(2);
	legendSpectraPi08TeV->SetMargin(0.2);
	legendSpectraPi08TeV->AddEntry(histoPCMYieldPi08TeV,collisionSystemPP8TeV.Data(),"pf");
	legendSpectraPi08TeV->AddEntry(histoPi0InputFullReweighted8TeV,"full MC reweighted","l");
	legendSpectraPi08TeV->AddEntry(histoPi0InputMCWOWeights_LHC14e2c_PHO_8TeV,"Phojet, LHC14e2c","l");
	legendSpectraPi08TeV->AddEntry(histoPi0InputMCWOWeights_LHC14e2a_PYT8_8TeV,"Pythia 8, LHC14e2a","l");
	legendSpectraPi08TeV->AddEntry(histoPi0InputMCWOWeights_LHC14e2b_PYT8_8TeV,"Pythia 8, LHC14e2b","l");
	legendSpectraPi08TeV->AddEntry(histoPi0InputMCWOWeights_LHC14e2b_PYT8addSig_8TeV,"Pythia 8, LHC14e2b, added Signals","l");
	legendSpectraPi08TeV->Draw();
	
	canvasPi0Spectra8TeV->Update();
	canvasPi0Spectra8TeV->Print(Form("%s/Pi0_Spectra_8TeV.%s",outputDir.Data(),suffix.Data()));

	//	**********************************************************************************************************************
	//	****************************************Fit Pi0 Spectra 8TeV and plot *********************************************
	//	**********************************************************************************************************************
	
	// fit spectrum with Tsallis function
	TF1* fitInvYieldDataPi0Comb8TeV 								= FitObject("l","fitInvYieldDataPi0Comb8TeV","Pi0");
	histoPCMYieldPi08TeV->Fit(fitInvYieldDataPi0Comb8TeV,"QNRMEI+","",0.3,12.);
	fitInvYieldDataPi0Comb8TeV->SetRange(0.1,20.);
	// print fit result to shell
	cout << WriteParameterToFile(fitInvYieldDataPi0Comb8TeV)<< endl;	

	// plot result with input
	canvasPi0Spectra8TeV->cd();
	histo2DPi0Spectra8TeV->DrawCopy(); 

	DrawGammaSetMarker(histoPCMYieldPi08TeV, markerStyleSpectrum8TeV, markerSizePP8TeV, color8TeV , color8TeV);
	histoPCMYieldPi08TeV->Draw("p,same,e1");

	fitInvYieldDataPi0Comb8TeV->SetLineColor(color8TeV);
	fitInvYieldDataPi0Comb8TeV->Draw("same");
	labelSpectraPi0Label->Draw();
	histoPi0InputFullReweighted8TeV->Draw("same,hist,c");
	
	canvasPi0Spectra8TeV->Print(Form("%s/Pi0_Spectra_WithFit_8TeV.%s",outputDir.Data(),suffix.Data()));

	// **********************************************************************************************************************
    // ******************************* Ratio of data to fit and MC input to fit for pi0 8 TeV ****************************
	// **********************************************************************************************************************

	// Calculation of ratio histograms
	TH1D* histoRatioPi0DatatoFit8TeV 									= CalculateHistoRatioToFit (histoPCMYieldPi08TeV, fitInvYieldDataPi0Comb8TeV,kTRUE);
	TH1D* histoRatioPi0MCtoDataFit8TeV 									= CalculateHistoRatioToFit (histoPi0InputFullReweighted8TeV, fitInvYieldDataPi0Comb8TeV,kTRUE);
	TH1D* histoRatioPi0MCAddSigtoDataFit8TeV 							= CalculateHistoRatioToFit (histoPi0InputFullAddSigReweighted8TeV, fitInvYieldDataPi0Comb8TeV,kTRUE);
	TH1D* histoRatioPi0MCUnweightedtoDataFit8TeV 						= NULL;
	if (histoPi0InputFull8TeV) histoRatioPi0MCUnweightedtoDataFit8TeV 	= CalculateHistoRatioToFit (histoPi0InputFull8TeV, fitInvYieldDataPi0Comb8TeV,kTRUE);
	
	canvasRatioToFit->cd();
	if (histoRatioPi0MCUnweightedtoDataFit8TeV) SetStyleHisto(histoRatioPi0MCUnweightedtoDataFit8TeV, 2, lineStyleMCB, 807 );
	SetStyleHisto(histoRatioPi0MCtoDataFit8TeV, 2, lineStyleMCA, kRed+2 );
	SetStyleHisto(histoRatioPi0MCAddSigtoDataFit8TeV, 2, lineStyleMCAddSig, kBlue+2 );
	DrawGammaSetMarker(histoRatioPi0DatatoFit8TeV, markerStyleSpectrum8TeV, markerSizePP8TeV, kBlack , kBlack);
	DrawAutoGammaMesonHistos( histoRatioPi0DatatoFit8TeV,
				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
				kFALSE, 1.5, 0, kTRUE,
				kTRUE, 0, 2.1,
				kTRUE, 0.,histoRatioPi0DatatoFit8TeV->GetXaxis()->GetBinUpEdge(histoRatioPi0DatatoFit8TeV->GetNbinsX()));
	histoRatioPi0DatatoFit8TeV->GetYaxis()->SetTitleOffset(0.9);
	histoRatioPi0DatatoFit8TeV->Draw("e,p");  
	if (runDrawReweighted) histoRatioPi0MCtoDataFit8TeV->Draw("same,hist,l");  
	if (runDrawReweighted) histoRatioPi0MCAddSigtoDataFit8TeV->Draw("same,hist,l");  
	if (histoRatioPi0MCUnweightedtoDataFit8TeV) histoRatioPi0MCUnweightedtoDataFit8TeV->Draw("same,hist,l");  

	labelSpectraPi0LabelRatio->Draw();
	TLatex *labelEnergy8TeVRatio 									= new TLatex(0.65,0.86,collisionSystemPP8TeV.Data());
	SetStyleTLatex( labelEnergy8TeVRatio, 0.035  ,4);
	labelEnergy8TeVRatio->Draw();
	
	TLegend* legendRatioPi08TeV = new TLegend(0.11,0.12,0.4,0.30);
	legendRatioPi08TeV->SetFillColor(0);
	legendRatioPi08TeV->SetLineColor(0);
	legendRatioPi08TeV->SetTextSize(0.035);
	legendRatioPi08TeV->SetMargin(0.2);
	legendRatioPi08TeV->AddEntry(histoRatioPi0DatatoFit8TeV,"Data/Tsallis fit to Data (0.3 <#it{p}_{T}<12)","p");
	if (runDrawReweighted) legendRatioPi08TeV->AddEntry(histoRatioPi0MCtoDataFit8TeV,Form("MC weighted %s/Tsallis  fit to Data (0.3 <#it{p}_{T}<12)",stringIterationNumber.Data()),"l");
	if (runDrawReweighted) legendRatioPi08TeV->AddEntry(histoRatioPi0MCAddSigtoDataFit8TeV,Form("MC add Sig weighted %s/Tsallis  fit to Data (0.3 <#it{p}_{T}<12)", 		
														stringIterationNumber.Data()),"l");
	if (histoRatioPi0MCUnweightedtoDataFit8TeV) legendRatioPi08TeV->AddEntry(histoRatioPi0MCUnweightedtoDataFit8TeV,"MC/Tsallis  fit to Data (0.3 <#it{p}_{T}<12)","l");
	legendRatioPi08TeV->Draw();
	
	DrawGammaLines(0., histoRatioPi0DatatoFit8TeV->GetXaxis()->GetBinUpEdge(histoRatioPi0DatatoFit8TeV->GetNbinsX()) ,1., 1.,0.1);
	canvasRatioToFit->Update();
	canvasRatioToFit->SaveAs(Form("%s/Pi0_RatioToDataFit_8TeV.%s",outputDir.Data(),suffix.Data()));
	
	//	**********************************************************************************************************************
	//	******************************Compare Efficiencies for Pi0 at 8TeV for different MC *******************************
	//	**********************************************************************************************************************
	TCanvas* canvasPi0Efficiencies8TeVDiffMC 					= new TCanvas("canvasPi0Efficiencies8TeVDiffMC", "", 200, 10, 1200, 1100);  // gives the page size
	DrawGammaCanvasSettings( canvasPi0Efficiencies8TeVDiffMC,  0.09, 0.01, 0.015, 0.08);
	canvasPi0Efficiencies8TeVDiffMC->SetLogy();
// 	canvasPi0Efficiencies8TeVDiffMC->SetLogx();
	
	TH2F * histo2DPi0Effi8TeV;
	histo2DPi0Effi8TeV 													= new TH2F("histo2DPi0Effi8TeV", "histo2DPi0Effi8TeV",1000, 0., 10.5, 1000, 1e-5, 1e-2 );
	SetStyleHistoTH2ForGraphs( histo2DPi0Effi8TeV, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{#pi^{0}}", 
							   0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
	histo2DPi0Effi8TeV->GetYaxis()->SetLabelOffset(0.01);
	histo2DPi0Effi8TeV->DrawCopy(); 

	DrawGammaSetMarker(histoPi0Efficiency_LHC14e2a_PYT8_8TeV, markerStyleMCA, markerSizePP8TeV, colorMCPythiaPP8TeV , colorMCPythiaPP8TeV);
	histoPi0Efficiency_LHC14e2a_PYT8_8TeV->Draw("p,same,e1");

	DrawGammaSetMarker(histoPi0Efficiency_LHC14e2b_PYT8addSig_8TeV, markerStyleMCAddSig, markerSizePP8TeV, colorMCPythiaAddSigPP8TeV , colorMCPythiaAddSigPP8TeV);
	histoPi0Efficiency_LHC14e2b_PYT8addSig_8TeV->Draw("p,same,e1");

	DrawGammaSetMarker(histoPi0Efficiency_LHC14e2b_PYT8_8TeV, markerStyleMCC, markerSizePP8TeV, colorMCPythia2PP8TeV , colorMCPythia2PP8TeV);
	histoPi0Efficiency_LHC14e2b_PYT8_8TeV->Draw("p,same,e1");
	
	DrawGammaSetMarker(histoPi0Efficiency_LHC14e2c_PHO_8TeV, markerStyleMCB, markerSizePP8TeV, colorMCPhojetPP8TeV , colorMCPhojetPP8TeV);
	histoPi0Efficiency_LHC14e2c_PHO_8TeV->Draw("p,same,e1");

	TLegend* legendEfficiencyPi08TeV 										= new TLegend(0.56,0.09,0.93,0.3);
	legendEfficiencyPi08TeV->SetFillColor(0);
	legendEfficiencyPi08TeV->SetLineColor(0);
	legendEfficiencyPi08TeV->SetTextSize(0.025);
	legendEfficiencyPi08TeV->SetMargin(0.2);
	legendEfficiencyPi08TeV->AddEntry(histoPi0Efficiency_LHC14e2c_PHO_8TeV,"Phojet, LHC14e2c","p");
	legendEfficiencyPi08TeV->AddEntry(histoPi0Efficiency_LHC14e2a_PYT8_8TeV,"Pythia 8, LHC14e2a","p");
	legendEfficiencyPi08TeV->AddEntry(histoPi0Efficiency_LHC14e2b_PYT8_8TeV,"Pythia 8, LHC14e2b","p");
	legendEfficiencyPi08TeV->AddEntry(histoPi0Efficiency_LHC14e2b_PYT8addSig_8TeV,"Pythia 8, LHC14e2b, added Signals","p");
	legendEfficiencyPi08TeV->Draw();
	
	canvasPi0Efficiencies8TeVDiffMC->Update();
	canvasPi0Efficiencies8TeVDiffMC->Print(Form("%s/Pi0_Efficiency_8TeV.%s",outputDir.Data(),suffix.Data()));

	canvasPi0Efficiencies8TeVDiffMC->SetLogy(0);
	canvasPi0Efficiencies8TeVDiffMC->SetTopMargin(0.035);
	
	histo2DPi0Effi8TeV->GetYaxis()->SetRangeUser(1e-5,3E-3);
	histo2DPi0Effi8TeV->DrawCopy(); 
	histoPi0Efficiency_LHC14e2a_PYT8_8TeV->Draw("p,same,e1");
	histoPi0Efficiency_LHC14e2b_PYT8addSig_8TeV->Draw("p,same,e1");
	histoPi0Efficiency_LHC14e2b_PYT8_8TeV->Draw("p,same,e1");
	histoPi0Efficiency_LHC14e2c_PHO_8TeV->Draw("p,same,e1");
	legendEfficiencyPi08TeV->Draw();
	
	canvasPi0Efficiencies8TeVDiffMC->Update();
	canvasPi0Efficiencies8TeVDiffMC->Print(Form("%s/Pi0_Efficiency_8TeV_LinY.%s",outputDir.Data(),suffix.Data()));
	
	//	**********************************************************************************************************************
	//	****************************************Eta Spectra 8TeV compared to MC********************************************
	//	**********************************************************************************************************************
	TCanvas* canvasEtaSpectra8TeV 								= new TCanvas("canvasEtaSpectra8TeV", "", 200, 10, 1200, 1100);  // gives the page size
	DrawGammaCanvasSettings( canvasEtaSpectra8TeV,  0.13, 0.01, 0.015, 0.08);
	canvasEtaSpectra8TeV->SetLogy();
	canvasEtaSpectra8TeV->SetLogx();
	
	TH2F * histo2DEtaSpectra8TeV;
	histo2DEtaSpectra8TeV 											= new TH2F("histo2DEtaSpectra8TeV", "histo2DEtaSpectra8TeV",1000, 0.23, 20., 1000, 1e-8, 1e1 );
	SetStyleHistoTH2ForGraphs( histo2DEtaSpectra8TeV, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
							   0.03, 0.04, 0.03, 0.04, 0.83, 1.4);
	histo2DEtaSpectra8TeV->GetXaxis()->SetLabelOffset(-0.01);
	histo2DEtaSpectra8TeV->GetYaxis()->SetLabelOffset(0.01);
	histo2DEtaSpectra8TeV->DrawCopy(); 
	
	DrawGammaSetMarker(histoPCMYieldEta8TeV, markerStyleSpectrum8TeV, markerSizePP8TeV, color8TeV , color8TeV);
	histoPCMYieldEta8TeV->Draw("p,same,e1");
	SetStyleHisto(histoEtaInputMCWOWeights_LHC14e2a_PYT8_8TeV, 1., lineStyleMCA, colorMCPythiaPP8TeV);  
	histoEtaInputMCWOWeights_LHC14e2a_PYT8_8TeV->Draw("same,hist,c");    
	SetStyleHisto(histoEtaInputMCWOWeights_LHC14e2c_PHO_8TeV, 1., lineStyleMCB, colorMCPhojetPP8TeV);  
	histoEtaInputMCWOWeights_LHC14e2c_PHO_8TeV->Draw("same,hist,c");    
	SetStyleHisto(histoEtaInputMCWOWeights_LHC14e2b_PYT8_8TeV, 1., lineStyleMCC, colorMCPythia2PP8TeV);  
	histoEtaInputMCWOWeights_LHC14e2b_PYT8_8TeV->Draw("same,hist,c");    
	SetStyleHisto(histoEtaInputMCWOWeights_LHC14e2b_PYT8addSig_8TeV, 1., lineStyleMCAddSig, colorMCPythiaAddSigPP8TeV);  
	histoEtaInputMCWOWeights_LHC14e2b_PYT8addSig_8TeV->Draw("same,hist,c");    
	SetStyleHisto(histoEtaInputFullReweighted8TeV, 1., lineStyleMCAddSig, kBlue);  
	histoEtaInputFullReweighted8TeV->Draw("same,hist,c");

	labelSpectraEtaLabel->Draw();
	
	TLegend* legendSpectraEta8TeV 										= new TLegend(0.16,0.09,0.73,0.32);
	legendSpectraEta8TeV->SetFillColor(0);
	legendSpectraEta8TeV->SetLineColor(0);
	legendSpectraEta8TeV->SetTextSize(0.025);
// 	legendSpectraEta8TeV->SetNColumns(2);
	legendSpectraEta8TeV->SetMargin(0.2);
	legendSpectraEta8TeV->AddEntry(histoPCMYieldEta8TeV,collisionSystemPP8TeV.Data(),"pf");
	legendSpectraEta8TeV->AddEntry(histoEtaInputFullReweighted8TeV,"full MC reweighted","l");
	legendSpectraEta8TeV->AddEntry(histoEtaInputMCWOWeights_LHC14e2c_PHO_8TeV,"Phojet, LHC14e2c","l");
	legendSpectraEta8TeV->AddEntry(histoEtaInputMCWOWeights_LHC14e2a_PYT8_8TeV,"Pythia 8, LHC14e2a","l");
	legendSpectraEta8TeV->AddEntry(histoEtaInputMCWOWeights_LHC14e2b_PYT8_8TeV,"Pythia 8, LHC14e2b","l");
	legendSpectraEta8TeV->AddEntry(histoEtaInputMCWOWeights_LHC14e2b_PYT8addSig_8TeV,"Pythia 8, LHC14e2b, added Signals","l");
	legendSpectraEta8TeV->Draw();
	
	canvasEtaSpectra8TeV->Update();
	canvasEtaSpectra8TeV->Print(Form("%s/Eta_Spectra_8TeV.%s",outputDir.Data(),suffix.Data()));

	//	**********************************************************************************************************************
	//	****************************************Fit Eta Spectra 8TeV and plot *********************************************
	//	**********************************************************************************************************************
	
	// fit spectrum with Tsallis function
	TF1* fitInvYieldDataEtaComb8TeV 								= FitObject("l","fitInvYieldDataEtaComb8TeV","Eta");
	histoPCMYieldEta8TeV->Fit(fitInvYieldDataEtaComb8TeV,"QNRMEI+","",0.5,8.);
	fitInvYieldDataEtaComb8TeV->SetRange(0.1,20.);
	// print fit result to shell
	cout << WriteParameterToFile(fitInvYieldDataEtaComb8TeV)<< endl;	

	// plot result with input
	canvasEtaSpectra8TeV->cd();
	histo2DEtaSpectra8TeV->DrawCopy(); 

	DrawGammaSetMarker(histoPCMYieldEta8TeV, markerStyleSpectrum8TeV, markerSizePP8TeV, color8TeV , color8TeV);
	histoPCMYieldEta8TeV->Draw("p,same,e1");

	fitInvYieldDataEtaComb8TeV->SetLineColor(color8TeV);
	fitInvYieldDataEtaComb8TeV->Draw("same");
	labelSpectraEtaLabel->Draw();
	histoEtaInputFullReweighted8TeV->Draw("same,hist,c");
	
	canvasEtaSpectra8TeV->Print(Form("%s/Eta_Spectra_WithFit_8TeV.%s",outputDir.Data(),suffix.Data()));

	// **********************************************************************************************************************
    // ******************************* Ratio of data to fit and MC input to fit for Eta 8 TeV ****************************
	// **********************************************************************************************************************
	
	// Calculation of ratio histograms
	TH1D* histoRatioEtaDatatoFit8TeV 		= CalculateHistoRatioToFit (histoPCMYieldEta8TeV, fitInvYieldDataEtaComb8TeV,kTRUE);
	TH1D* histoRatioEtaMCtoDataFit8TeV 		= CalculateHistoRatioToFit (histoEtaInputFullReweighted8TeV, fitInvYieldDataEtaComb8TeV,kTRUE);
	TH1D* histoRatioEtaMCAddSigtoDataFit8TeV = CalculateHistoRatioToFit (histoEtaInputFullAddSigReweighted8TeV, fitInvYieldDataEtaComb8TeV,kTRUE);
	TH1D* histoRatioEtaMCUnweightedtoDataFit8TeV = NULL;
	if (histoEtaInputFull8TeV) histoRatioEtaMCUnweightedtoDataFit8TeV = CalculateHistoRatioToFit (histoEtaInputFull8TeV, fitInvYieldDataEtaComb8TeV,kTRUE);
	if (histoRatioEtaMCUnweightedtoDataFit8TeV) SetStyleHisto(histoRatioEtaMCUnweightedtoDataFit8TeV, 2, lineStyleMCB, 807 );
	// plotting ratios
	canvasRatioToFit->cd();
	SetStyleHisto(histoRatioEtaMCtoDataFit8TeV, 2, lineStyleMCA, kRed+2 );
	SetStyleHisto(histoRatioEtaMCAddSigtoDataFit8TeV, 2, lineStyleMCAddSig, kBlue+2 );
	DrawGammaSetMarker(histoRatioEtaDatatoFit8TeV, markerStyleSpectrum8TeV, markerSizePP8TeV, kBlack , kBlack);
	DrawAutoGammaMesonHistos( histoRatioEtaDatatoFit8TeV,
				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
				kFALSE, 1.5, 0, kTRUE,
				kTRUE, 0, 2.1,
				kTRUE, 0.,histoRatioEtaDatatoFit8TeV->GetXaxis()->GetBinUpEdge(histoRatioEtaDatatoFit8TeV->GetNbinsX()));
	histoRatioEtaDatatoFit8TeV->GetYaxis()->SetTitleOffset(0.9);
	histoRatioEtaDatatoFit8TeV->Draw("e,p");  
	if (runDrawReweighted) histoRatioEtaMCtoDataFit8TeV->Draw("same,hist,l");  
	if (runDrawReweighted) histoRatioEtaMCAddSigtoDataFit8TeV->Draw("same,hist,l");  
	if (histoRatioEtaMCUnweightedtoDataFit8TeV) histoRatioEtaMCUnweightedtoDataFit8TeV->Draw("same,hist,l");  

	// labeling
	labelSpectraEtaLabelRatio->Draw();
	labelEnergy8TeVRatio->Draw();
	
	TLegend* legendRatioEta8TeV = new TLegend(0.11,0.12,0.4,0.30);
	legendRatioEta8TeV->SetFillStyle(0);
	legendRatioEta8TeV->SetFillColor(0);
	legendRatioEta8TeV->SetLineColor(0);
	legendRatioEta8TeV->SetTextSize(0.035);
	legendRatioEta8TeV->SetMargin(0.2);
	legendRatioEta8TeV->AddEntry(histoRatioEtaDatatoFit8TeV,"Data/Tsallis fit to Data (0.5 <#it{p}_{T}<8)","p");
	if (runDrawReweighted) legendRatioEta8TeV->AddEntry(histoRatioEtaMCtoDataFit8TeV,Form("MC weighted %s/Tsallis  fit to Data (0.5 <#it{p}_{T}<8)",stringIterationNumber.Data()),"l");
	if (runDrawReweighted) legendRatioEta8TeV->AddEntry(histoRatioEtaMCAddSigtoDataFit8TeV,Form("MC add Sig weighted %s/Tsallis  fit to Data (0.5 <#it{p}_{T}<8)",
														   stringIterationNumber.Data()), "l");
	if (histoRatioEtaMCUnweightedtoDataFit8TeV) legendRatioEta8TeV->AddEntry(histoRatioEtaMCUnweightedtoDataFit8TeV,"MC/Tsallis  fit to Data (0.5 <#it{p}_{T}<8)","l");
	legendRatioEta8TeV->Draw();
	
	DrawGammaLines(0., histoRatioEtaDatatoFit8TeV->GetXaxis()->GetBinUpEdge(histoRatioEtaDatatoFit8TeV->GetNbinsX()) ,1., 1.,0.1);
	canvasRatioToFit->Update();
	canvasRatioToFit->SaveAs(Form("%s/Eta_RatioToDataFit_8TeV.%s",outputDir.Data(),suffix.Data()));

	
	//	**********************************************************************************************************************
	//	******************************Compare Efficiencies for Eta at 8TeV for different MC *******************************
	//	**********************************************************************************************************************
	TCanvas* canvasEtaEfficiencies8TeVDiffMC 					= new TCanvas("canvasEtaEfficiencies8TeVDiffMC", "", 200, 10, 1200, 1100);  // gives the page size
	DrawGammaCanvasSettings( canvasEtaEfficiencies8TeVDiffMC,  0.09, 0.01, 0.015, 0.08);
	canvasEtaEfficiencies8TeVDiffMC->SetLogy();
// 	canvasEtaEfficiencies8TeVDiffMC->SetLogx();
	
	TH2F * histo2DEtaEffi8TeV;
	histo2DEtaEffi8TeV 													= new TH2F("histo2DEtaEffi8TeV", "histo2DEtaEffi8TeV",1000, 0., 10.5, 1000, 1e-5, 1e-2 );
	SetStyleHistoTH2ForGraphs( histo2DEtaEffi8TeV, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{#eta}", 
							   0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
	histo2DEtaEffi8TeV->GetYaxis()->SetLabelOffset(0.01);
	histo2DEtaEffi8TeV->DrawCopy(); 

	DrawGammaSetMarker(histoEtaEfficiency_LHC14e2a_PYT8_8TeV, markerStyleMCA, markerSizePP8TeV, colorMCPythiaPP8TeV , colorMCPythiaPP8TeV);
	histoEtaEfficiency_LHC14e2a_PYT8_8TeV->Draw("p,same,e1");

	DrawGammaSetMarker(histoEtaEfficiency_LHC14e2b_PYT8addSig_8TeV, markerStyleMCAddSig, markerSizePP8TeV, colorMCPythiaAddSigPP8TeV , colorMCPythiaAddSigPP8TeV);
	histoEtaEfficiency_LHC14e2b_PYT8addSig_8TeV->Draw("p,same,e1");

	DrawGammaSetMarker(histoEtaEfficiency_LHC14e2b_PYT8_8TeV, markerStyleMCC, markerSizePP8TeV, colorMCPythia2PP8TeV , colorMCPythia2PP8TeV);
	histoEtaEfficiency_LHC14e2b_PYT8_8TeV->Draw("p,same,e1");
	
	DrawGammaSetMarker(histoEtaEfficiency_LHC14e2c_PHO_8TeV, markerStyleMCB, markerSizePP8TeV, colorMCPhojetPP8TeV , colorMCPhojetPP8TeV);
	histoEtaEfficiency_LHC14e2c_PHO_8TeV->Draw("p,same,e1");

	TLegend* legendEfficiencyEta8TeV 										= new TLegend(0.56,0.09,0.93,0.3);
	legendEfficiencyEta8TeV->SetFillColor(0);
	legendEfficiencyEta8TeV->SetLineColor(0);
	legendEfficiencyEta8TeV->SetTextSize(0.025);
	legendEfficiencyEta8TeV->SetMargin(0.2);
	legendEfficiencyEta8TeV->AddEntry(histoEtaEfficiency_LHC14e2c_PHO_8TeV,"Phojet, LHC14e2c","p");
	legendEfficiencyEta8TeV->AddEntry(histoEtaEfficiency_LHC14e2a_PYT8_8TeV,"Pythia 8, LHC14e2a","p");
	legendEfficiencyEta8TeV->AddEntry(histoEtaEfficiency_LHC14e2b_PYT8_8TeV,"Pythia 8, LHC14e2b","p");
	legendEfficiencyEta8TeV->AddEntry(histoEtaEfficiency_LHC14e2b_PYT8addSig_8TeV,"Pythia 8, LHC14e2b, added Signals","p");
	legendEfficiencyEta8TeV->Draw();
	
	canvasEtaEfficiencies8TeVDiffMC->Update();
	canvasEtaEfficiencies8TeVDiffMC->Print(Form("%s/Eta_Efficiency_8TeV.%s",outputDir.Data(),suffix.Data()));

	canvasEtaEfficiencies8TeVDiffMC->SetLogy(0);
	canvasEtaEfficiencies8TeVDiffMC->SetTopMargin(0.035);
	
	histo2DEtaEffi8TeV->GetYaxis()->SetRangeUser(1e-5,3E-3);
	histo2DEtaEffi8TeV->DrawCopy(); 
	histoEtaEfficiency_LHC14e2a_PYT8_8TeV->Draw("p,same,e1");
	histoEtaEfficiency_LHC14e2b_PYT8addSig_8TeV->Draw("p,same,e1");
	histoEtaEfficiency_LHC14e2b_PYT8_8TeV->Draw("p,same,e1");
	histoEtaEfficiency_LHC14e2c_PHO_8TeV->Draw("p,same,e1");
	legendEfficiencyEta8TeV->Draw();
	
	canvasEtaEfficiencies8TeVDiffMC->Update();
	canvasEtaEfficiencies8TeVDiffMC->Print(Form("%s/Eta_Efficiency_8TeV_LinY.%s",outputDir.Data(),suffix.Data()));
	

	
	//	**********************************************************************************************************************
	//	****************************************Write fits & input MC spectra to file ****************************************
	//	**********************************************************************************************************************	
	TFile fMCSpectraInput("MCSpectraInputpp.root","UPDATE");
		if (fitInvYieldDataPi0Comb2760GeV){
			fitInvYieldDataPi0Comb2760GeV->SetRange(0,30);
			fitInvYieldDataPi0Comb2760GeV->Write("Pi0_Fit_Data_2760GeV",TObject::kOverwrite);
		}
		if (fitInvYieldDataEtaComb2760GeV){
			fitInvYieldDataEtaComb2760GeV->SetRange(0,30);
			fitInvYieldDataEtaComb2760GeV->Write("Eta_Fit_Data_2760GeV",TObject::kOverwrite);
		}
		if (histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV){
			histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV->SetTitle("Pi0_Phojet_LHC12f1b_WOSDD_2760GeV");
			histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV->Write("Pi0_Phojet_LHC12f1b_WOSDD_2760GeV",TObject::kOverwrite);
		}
		if (histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV){
			histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV->SetTitle("Pi0_Pythia8_LHC12f1a_WOSDD_2760GeV");
			histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV->Write("Pi0_Pythia8_LHC12f1a_WOSDD_2760GeV",TObject::kOverwrite);
		}	
		if (histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV){
			histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV->SetTitle("Pi0_Pythia8_LHC12i3_WOSDD_2760GeV");
			histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV->Write("Pi0_Pythia8_LHC12i3_WOSDD_2760GeV",TObject::kOverwrite);
		}
		if (histoPi0InputFullAddSig2760GeV){
			histoPi0InputFullAddSig2760GeV->SetTitle("Pi0_Pythia8_LHC12i3_WOSDD_addSig_2760GeV");
			histoPi0InputFullAddSig2760GeV->Write("Pi0_Pythia8_LHC12i3_WOSDD_addSig_2760GeV",TObject::kOverwrite);
		}
		if (histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV){
			histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV->SetTitle("Eta_Phojet_LHC12f1b_WOSDD_2760GeV");
			histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV->Write("Eta_Phojet_LHC12f1b_WOSDD_2760GeV",TObject::kOverwrite);
		}
		if (histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV){
			histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV->SetTitle("Eta_Pythia8_LHC12f1a_WOSDD_2760GeV");
			histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV->Write("Eta_Pythia8_LHC12f1a_WOSDD_2760GeV",TObject::kOverwrite);
		}	
		if (histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV){
			histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV->SetTitle("Eta_Pythia8_LHC12i3_WOSDD_2760GeV");
			histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV->Write("Eta_Pythia8_LHC12i3_WOSDD_2760GeV",TObject::kOverwrite);
		}
		if (histoEtaInputFullAddSig2760GeV){
			histoEtaInputFullAddSig2760GeV->SetTitle("Eta_Pythia8_LHC12i3_WOSDD_addSig_2760GeV");
			histoEtaInputFullAddSig2760GeV->Write("Eta_Pythia8_LHC12i3_WOSDD_addSig_2760GeV",TObject::kOverwrite);
		}
		if (histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD){
			histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD->SetTitle("Pi0_Phojet_LHC12f1b_WSDD_2760GeV");
			histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD->Write("Pi0_Phojet_LHC12f1b_WSDD_2760GeV",TObject::kOverwrite);
		}
		if (histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD){
			histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD->SetTitle("Pi0_Pythia8_LHC12f1a_WSDD_2760GeV");
			histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD->Write("Pi0_Pythia8_LHC12f1a_WSDD_2760GeV",TObject::kOverwrite);
		}	
		if (histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD){
			histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD->SetTitle("Pi0_Pythia8_LHC12i3_WSDD_2760GeV");
			histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD->Write("Pi0_Pythia8_LHC12i3_WSDD_2760GeV",TObject::kOverwrite);
		}
		if (histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD){
			histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD->SetTitle("Pi0_Pythia8_LHC12i3_WSDD_addSig_2760GeV");
			histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD->Write("Pi0_Pythia8_LHC12i3_WSDD_addSig_2760GeV",TObject::kOverwrite);
		}
		if (histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD){
			histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD->SetTitle("Eta_Phojet_LHC12f1b_WSDD_2760GeV");
			histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD->Write("Eta_Phojet_LHC12f1b_WSDD_2760GeV",TObject::kOverwrite);
		}
		if (histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD){
			histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD->SetTitle("Eta_Pythia8_LHC12f1a_WSDD_2760GeV");
			histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD->Write("Eta_Pythia8_LHC12f1a_WSDD_2760GeV",TObject::kOverwrite);
		}	
		if (histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD){
			histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD->SetTitle("Eta_Pythia8_LHC12i3_WSDD_2760GeV");
			histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD->Write("Eta_Pythia8_LHC12i3_WSDD_2760GeV",TObject::kOverwrite);
		}
		if (histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD){
			histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD->SetTitle("Eta_Pythia8_LHC12i3_WSDD_addSig_2760GeV");
			histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD->Write("Eta_Pythia8_LHC12i3_WSDD_addSig_2760GeV",TObject::kOverwrite);
		}

		if (fitInvYieldDataPi0Comb8TeV){
			fitInvYieldDataPi0Comb8TeV->SetRange(0,30);
			fitInvYieldDataPi0Comb8TeV->Write("Pi0_Fit_Data_8TeV",TObject::kOverwrite);
		}
		if (fitInvYieldDataEtaComb8TeV){
			fitInvYieldDataEtaComb8TeV->SetRange(0,30);
			fitInvYieldDataEtaComb8TeV->Write("Eta_Fit_Data_8TeV",TObject::kOverwrite);
		}
		if (histoPi0InputMCWOWeights_LHC14e2c_PHO_8TeV){
			histoPi0InputMCWOWeights_LHC14e2c_PHO_8TeV->SetTitle("Pi0_Phojet_LHC14e2c_8TeV");
			histoPi0InputMCWOWeights_LHC14e2c_PHO_8TeV->Write("Pi0_Phojet_LHC14e2c_8TeV",TObject::kOverwrite);
		}
		if (histoPi0InputMCWOWeights_LHC14e2a_PYT8_8TeV){
			histoPi0InputMCWOWeights_LHC14e2a_PYT8_8TeV->SetTitle("Pi0_Pythia8_LHC14e2a_8TeV");
			histoPi0InputMCWOWeights_LHC14e2a_PYT8_8TeV->Write("Pi0_Pythia8_LHC14e2a_8TeV",TObject::kOverwrite);
		}	
		if (histoPi0InputMCWOWeights_LHC14e2b_PYT8_8TeV){
			histoPi0InputMCWOWeights_LHC14e2b_PYT8_8TeV->SetTitle("Pi0_Pythia8_LHC14e2b_8TeV");
			histoPi0InputMCWOWeights_LHC14e2b_PYT8_8TeV->Write("Pi0_Pythia8_LHC14e2b_8TeV",TObject::kOverwrite);
		}
		if (histoPi0InputFullAddSig8TeV){
			histoPi0InputFullAddSig8TeV->SetTitle("Pi0_Pythia8_LHC14e2b_addSig_8TeV");
			histoPi0InputFullAddSig8TeV->Write("Pi0_Pythia8_LHC14e2b_addSig_8TeV",TObject::kOverwrite);
		}
		if (histoEtaInputMCWOWeights_LHC14e2c_PHO_8TeV){
			histoEtaInputMCWOWeights_LHC14e2c_PHO_8TeV->SetTitle("Eta_Phojet_LHC14e2c_8TeV");
			histoEtaInputMCWOWeights_LHC14e2c_PHO_8TeV->Write("Eta_Phojet_LHC14e2c_8TeV",TObject::kOverwrite);
		}
		if (histoEtaInputMCWOWeights_LHC14e2a_PYT8_8TeV){
			histoEtaInputMCWOWeights_LHC14e2a_PYT8_8TeV->SetTitle("Eta_Pythia8_LHC14e2a_8TeV");
			histoEtaInputMCWOWeights_LHC14e2a_PYT8_8TeV->Write("Eta_Pythia8_LHC14e2a_8TeV",TObject::kOverwrite);
		}	
		if (histoEtaInputMCWOWeights_LHC14e2b_PYT8_8TeV){
			histoEtaInputMCWOWeights_LHC14e2b_PYT8_8TeV->SetTitle("Eta_Pythia8_LHC14e2b_8TeV");
			histoEtaInputMCWOWeights_LHC14e2b_PYT8_8TeV->Write("Eta_Pythia8_LHC14e2b_8TeV",TObject::kOverwrite);
		}
		if (histoEtaInputFullAddSig8TeV){
			histoEtaInputFullAddSig8TeV->SetTitle("Eta_Pythia8_LHC14e2b_addSig_8TeV");
			histoEtaInputFullAddSig8TeV->Write("Eta_Pythia8_LHC14e2b_addSig_8TeV",TObject::kOverwrite);
		}

	fMCSpectraInput.Close();
		
}

