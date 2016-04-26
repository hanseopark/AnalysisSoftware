#include <Riostream.h>
#include <fstream>
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
#include "TGaxis.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
//#include "AliHEPDataParser.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"


void CombineGammaResultsPbPb( TString inputFileNamePCM = "", TString inputFileNamePHOS = "", TString suffix = "eps", Bool_t enablepValueCalc = kFALSE){

	//*******************************************************************************************************************************************
	//*********************************************************** Set main style choices ********************************************************
	//*******************************************************************************************************************************************
	StyleSettingsThesis();
	SetPlotStyle();

	//*******************************************************************************************************************************************
	//*****************************************************Definition pp & pPb input ************************************************************
	//*******************************************************************************************************************************************
	TString inputFileNamePCMpp									= "ExternalInput/PCM/Gamma_PCMResults_pp.root";
	TString inputFileNamePCMpPb									= "ExternalInputpPb/PCM/Gamma_PCMResults_pPb_5.023TeV.root";
	
	//*******************************************************************************************************************************************
	//********************************************* Create output directory and copy input files there ******************************************
	//*******************************************************************************************************************************************
	TString dateForOutput 										= ReturnDateStringForOutput();
	TString outputDir 											= Form("%s/%s/CombineGammaMeasurementsPbPb",suffix.Data(),dateForOutput.Data());
	TString fileNameTheoryPbPb									= "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";
	TString fileNameExperimentPbPb								= "ExternalInputPbPb/OtherExperiments/phenix_200.root";
	gSystem->Exec("mkdir -p "+outputDir);
	gSystem->Exec(Form("cp %s %s/InputPCMGammaPbPb.root", inputFileNamePCM.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPHOSGammaPbPb.root", inputFileNamePHOS.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPCMGammapPb.root", inputFileNamePCMpPb.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPCMGammapp.root", inputFileNamePCMpp.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/TheoryPbPb.root", fileNameTheoryPbPb.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/OtherExperimentsPbPb.root", fileNameExperimentPbPb.Data(), outputDir.Data()));
	
	TString nameFinalResDat 									= Form("%s/CombinedResultsGamma_FitResults_%s.dat",outputDir.Data(),dateForOutput.Data());
	fstream fileFinalResults;
	fileFinalResults.open(nameFinalResDat.Data(), ios::out);	

	//*******************************************************************************************************************************************
	//******************************************************* set ranges for plotting ***********************************************************
	//*******************************************************************************************************************************************
	Double_t doubleRatio[2];
	Double_t indMeasRatio[2];
// 	Double_t incRatio[2];
	Double_t doubleRatioX[2];
	Double_t doubleRatioXpp[2];
	doubleRatio[0] 		= 0.75; 	doubleRatio[1] 		= 1.95;
	indMeasRatio[0] 	= 0.65;	 	indMeasRatio[1]		= 1.45;
// 	incRatio[0] 		= 0.0; 		incRatio[1] 		= 1.7;
	doubleRatioX[0] 	= 0.7; 		doubleRatioX[1] 	= 16;
	doubleRatioXpp[0] 	= 0.23; 		doubleRatioXpp[1] 	= 23;
	
	Double_t nCollErr0020 										= GetNCollErrFromName("0020");
	Double_t nCollErr2040 										= GetNCollErrFromName("2040");
	Double_t nCollErr4080 										= GetNCollErrFromName("4080");

	Double_t nColl0020 											= GetNCollFromName("0020");
	Double_t nColl2040 											= GetNCollFromName("2040");
	Double_t nColl4080 											= GetNCollFromName("4080");
	
	Double_t normErr0020 										= nCollErr0020/nColl0020;
	Double_t normErr2040 										= nCollErr2040/nColl2040;
	Double_t normErr4080 										= nCollErr4080/nColl4080;
	
	Color_t colorCocktailPi0 									= kRed+2;
	Color_t colorCocktailEta 									= kBlue+1;
	Color_t colorCocktailEtaP 									= kOrange+1;
	Color_t colorCocktailOmega 									= kYellow+2;
	Color_t colorCocktailPhi 									= kViolet;
	Color_t colorCocktailRho0 									= kAzure-2;
	Color_t colorCocktailSigma0 								= kGray+1;
	
	Color_t colorNLOcalc 										= kBlack;//kBlue-7;
	Style_t fillStyleNLO										= 1001;
	Color_t colorEPS09calc 										= kGray;
	Color_t colorEPS09calc2 									= kGray+1;
	Style_t fillStyleEPS09										= 3008;
	Color_t colorCT10calc 										= kAzure-4;
	Style_t fillStyleCT10										= 3001;
	Color_t colorNLOMcGill 										= kOrange+1;
	Style_t styleNLOMcGill 										= 7;
	Color_t colorPHSD 											= kGray+2;
	Style_t stylePHSD	 										= 2;
	Color_t colorChatterjee 									= kBlue+1;
	Style_t styleChatterjee	 									= 4;
	Color_t colorHees 											= kBlack;
	Style_t styleHees	 										= 3;
	Color_t colorHe 											= kBlue-7;
	Style_t styleHe	 											= 8;
	Style_t styleFit											= 1;
	
	Color_t colorPHENIX											= kGray+2;
	Style_t markerStylePHENIX									= 24;
	Size_t markerSizePHENIX										= 2;
	
	Color_t colorPCM 											= GetDefaultColorDiffDetectors("PCM", kFALSE, kFALSE, kTRUE);
	Color_t colorPHOS	 										= GetDefaultColorDiffDetectors("PHOS", kFALSE, kFALSE, kTRUE);
	Color_t colorPCMBox 										= GetDefaultColorDiffDetectors("PCM", kFALSE, kTRUE, kTRUE);
	Color_t colorPHOSBox	 									= GetDefaultColorDiffDetectors("PHOS", kFALSE, kTRUE, kTRUE);
	Style_t markerStylePCM 										= GetDefaultMarkerStyleDiffDetectors("PCM", kFALSE);
	Style_t markerStylePHOS 									= GetDefaultMarkerStyleDiffDetectors("PHOS", kFALSE);
	Size_t markerSizePCM 										= GetDefaultMarkerSizeDiffDetectors("PCM", kFALSE);
	Size_t markerSizePHOS 										= GetDefaultMarkerSizeDiffDetectors("PHOS", kFALSE);
	
	Color_t colorComb0020										= GetColorDefaultColor("PbPb_2.76TeV", "", "0-20%");
	Color_t colorComb2040										= GetColorDefaultColor("PbPb_2.76TeV", "", "20-40%");
	Color_t colorComb4080										= GetColorDefaultColor("PbPb_2.76TeV", "", "40-80%");
	Color_t colorComb0020Box									= GetColorDefaultColor("PbPb_2.76TeV", "", "0-20%", kTRUE);
	Color_t colorComb2040Box									= GetColorDefaultColor("PbPb_2.76TeV", "", "20-40%", kTRUE);
	Color_t colorComb4080Box									= GetColorDefaultColor("PbPb_2.76TeV", "", "40-80%", kTRUE);
	Color_t colorCombpPb										= GetColorDefaultColor("pPb_5.023TeV", "", "");
	Color_t colorCombpPbBox										= GetColorDefaultColor("pPb_5.023TeV", "", "", kTRUE);
	Color_t colorCombpp7TeV										= GetColorDefaultColor("7TeV", "", "");
	Color_t colorCombpp7TeVBox									= GetColorDefaultColor("7TeV", "", "", kTRUE);
	Color_t colorCombpp2760GeV									= GetColorDefaultColor("2.76TeV", "", "");
	Color_t colorCombpp2760GeVBox								= GetColorDefaultColor("2.76TeV", "", "", kTRUE);
	
	Style_t markerStyleComb0020									= GetDefaultMarkerStyle("PbPb_2.76TeV", "", "0-20%");
	Style_t markerStyleComb2040									= GetDefaultMarkerStyle("PbPb_2.76TeV", "", "20-40%");
	Style_t markerStyleComb4080									= GetDefaultMarkerStyle("PbPb_2.76TeV", "", "40-80%");
	Style_t markerStyleCombpPb									= GetDefaultMarkerStyle("pPb_5.023TeV", "", "");
	Style_t markerStyleCombpp7TeV								= GetDefaultMarkerStyle("7TeV", "", "");
	Style_t markerStyleCombpp2760GeV							= GetDefaultMarkerStyle("2.76TeV", "", "");
	Size_t markerSizeComb0020									= GetDefaultMarkerSize("PbPb_2.76TeV", "", "0-20%");
	Size_t markerSizeComb2040									= GetDefaultMarkerSize("PbPb_2.76TeV", "", "20-40%");
	Size_t markerSizeComb4080									= GetDefaultMarkerSize("PbPb_2.76TeV", "", "40-80%");
	Size_t markerSizeCombpPb									= GetDefaultMarkerSize("pPb_5.023TeV", "", "");
	Size_t markerSizeCombpp7TeV									= GetDefaultMarkerSize("7TeV", "", "");
	Size_t markerSizeCombpp2760GeV								= GetDefaultMarkerSize("2.76TeV", "", "");

	Width_t widthLinesBoxes										= 1.4;
	Width_t widthCommonFit										= 2.4;
	
	TString collisionSystem			 							= "Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
	TString collisionSystemCent0020 							= "0-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
	TString collisionSystemCent2040 							= "20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
	TString collisionSystemCent4080 							= "40-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
	TString collisionSystempPb 									= "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
	TString collisionSystempp7TeV 								= "pp, #sqrt{#it{s}} = 7 TeV";
	TString collisionSystempp2760GeV 							= "pp, #sqrt{#it{s}} = 2.76 TeV";
	TString collisionSystemRHIC0020			 					= "0-20% Au-Au #sqrt{#it{s}_{_{NN}}} = 0.2 TeV";
	TString collisionSystemRHIC2040			 					= "20-40% Au-Au #sqrt{#it{s}_{_{NN}}} = 0.2 TeV";
	
	TFile*				fileTheoryPbPb 						= new TFile( fileNameTheoryPbPb.Data());
	TDirectory* 		directoryTheoryGamma 				= (TDirectory*)fileTheoryPbPb->Get("DirectPhoton"); 
		TGraphErrors* 	graphTheoryMcGill0020				= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_McGill_0020");
		TGraphErrors* 	graphTheoryMcGill2040				= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_McGill_2040");
		TGraphErrors* 	graphTheoryPromptMcGill0020			= (TGraphErrors*) directoryTheoryGamma->Get("graphPromptPhotonYield_McGill_0020");
		TF1* fitTheoryPromptMcGill0020	= new TF1 ("fitTheoryPromptMcGill0020",
												   "[0]/64.*1e-9*TMath::Exp(16.20-3.94*TMath::Log(x)-0.269*TMath::Log(x)**2)*1./(0.865779/TMath::Power(4,0.0694875))",
												   0.1, 20);  
		fitTheoryPromptMcGill0020->SetParameter(0,nColl0020);		
		TGraphErrors* 	graphTheoryPromptMcGill2040			= (TGraphErrors*) directoryTheoryGamma->Get("graphPromptPhotonYield_McGill_2040");
		TF1* fitTheoryPromptMcGill2040	= new TF1 ("fitTheoryPromptMcGill2040",
												   "[0]/64.*1e-9*TMath::Exp(16.20-3.94*TMath::Log(x)-0.269*TMath::Log(x)**2)*1./(0.865779/TMath::Power(4,0.0694875))",
												   0.1, 20);  
		fitTheoryPromptMcGill2040->SetParameter(0,nColl2040);
		TF1* fitTheoryPromptMcGill4080	= new TF1 ("fitTheoryPromptMcGill4080",
												   "[0]/64.*1e-9*TMath::Exp(16.20-3.94*TMath::Log(x)-0.269*TMath::Log(x)**2)*1./(0.865779/TMath::Power(4,0.0694875))",
												   0.1, 20);  
		fitTheoryPromptMcGill4080->SetParameter(0,nColl4080);

		TGraphErrors* 	graphTheoryPHSD0020					= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_PHSD_0020");
		TGraphErrors* 	graphTheoryPHSD2040					= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_PHSD_2040");
		TGraphErrors* 	graphTheoryPHSD4080					= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_PHSD_4080");
		TGraphErrors* 	graphTheoryChatterjee0020			= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_0020_2");
		TGraphErrors* 	graphTheoryChatterjee2040			= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_2040_2");
		TGraphErrors* 	graphTheoryChatterjee4060			= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_Chatterjee_4060");
		TGraphErrors* 	graphTheoryHees0020					= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_VanHees_0020");
		TGraphErrors* 	graphTheoryHees2040					= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_VanHees_2040");
		TGraphErrors* 	graphTheoryHees4080					= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_VanHees_4080");
		TGraphErrors* 	graphTheoryHe0020					= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_He_0020");
		TGraphErrors* 	graphTheoryHe2040					= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_He_2040");
		TGraphErrors* 	graphTheoryHe4080					= (TGraphErrors*) directoryTheoryGamma->Get("graphDirectPhotonYield_He_4080");
	
	TFile*				fileExperimentPbPb					= new TFile( fileNameExperimentPbPb.Data() );
		TGraphAsymmErrors* graphPHENIXAuAuStat0020			= (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_0_20_StatErr");
		TGraphAsymmErrors* graphPHENIXAuAuSys0020			= (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_0_20_SysErr");
		TGraphAsymmErrors* graphPHENIXAuAuStat2040			= (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_20_40_StatErr");
		TGraphAsymmErrors* graphPHENIXAuAuSys2040			= (TGraphAsymmErrors*)fileExperimentPbPb->Get("phenix_auau_20_40_SysErr");
		
	//*******************************************************************************************************************************************
	//*********************************************** Load PCM histograms from PCM file *********************************************************
	//*******************************************************************************************************************************************
	TFile*				filePCMGamma 						= new TFile( inputFileNamePCM.Data());
	//________________________________________________ Load PCM 0-20% ___________________________________________________________________________
	TDirectory* 		directoryPCMGamma0020 					= (TDirectory*)filePCMGamma->Get("Gamma_PbPb_2.76TeV_0-20%"); 
		TH1D * 				histoPCMDRPi0FitStatErr0020 		= (TH1D*) directoryPCMGamma0020->Get("DoubleRatioPi0FitStatError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysErr0020 			= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("DoubleRatioPi0FitSystError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysAErr0020			= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("DoubleRatioPi0FitSystErrorA");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysBErr0020			= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("DoubleRatioPi0FitSystErrorB");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysCErr0020			= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("DoubleRatioPi0FitSystErrorC");
// 		TH1D * 				histoPCMIncRStatErr0020				= (TH1D*) directoryPCMGamma0020->Get("IncRatioStatError");
// 		TGraphAsymmErrors* 	graphPCMIncRSysErr0020 				= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncRatioSystError");
// 		TGraphAsymmErrors* 	graphPCMIncRSysAErr0020				= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncRatioSystErrorA");
// 		TGraphAsymmErrors* 	graphPCMIncRSysBErr0020				= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncRatioSystErrorB");
// 		TGraphAsymmErrors* 	graphPCMIncRSysCErr0020				= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncRatioSystErrorC");
// 		TH1D * 				histoPCMIncRPi0FitStatErr0020		= (TH1D*) directoryPCMGamma0020->Get("IncRatioPi0FitStatError");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysErr0020		= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncRatioPi0FitSystError");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysAErr0020		= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncRatioPi0FitSystErrorA");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysBErr0020		= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncRatioPi0FitSystErrorB");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysCErr0020		= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncRatioPi0FitSystErrorC");
		TH1D * 				histoPCMIncGammaStatErr0020			= (TH1D*) directoryPCMGamma0020->Get("IncGammaStatError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysErr0020			= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncGammaSystError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysAErr0020			= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncGammaSystErrorA");
		TGraphAsymmErrors* 	graphPCMIncGammaSysBErr0020			= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncGammaSystErrorB");
		TGraphAsymmErrors* 	graphPCMIncGammaSysCErr0020			= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("IncGammaSystErrorC");
// 		TH1D * 				histoPCMPi0StatErr0020				= (TH1D*) directoryPCMGamma0020->Get("Pi0StatError");
// 		TGraphAsymmErrors* 	graphPCMPi0SysErr0020				= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("Pi0SystError");
// 		TH1D * 				histoPCMPi0FitStatErr0020			= (TH1D*) directoryPCMGamma0020->Get("Pi0FitStatError");
// 		TGraphAsymmErrors* 	graphPCMPi0FitSysErr0020			= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("Pi0FitSystError");
		TH1D* 				histoCocktailSumGamma				= (TH1D*) directoryPCMGamma0020->Get("CocktailSumGamma");
		TH1D* 				histoCocktailPi0Gamma				= (TH1D*) directoryPCMGamma0020->Get("CocktailPi0Gamma");
		TH1D* 				histoCocktailEtaGamma				= (TH1D*) directoryPCMGamma0020->Get("CocktailEtaGamma");
		TH1D* 				histoCocktailOmegaGamma				= (TH1D*) directoryPCMGamma0020->Get("CocktailOmegaGamma");
		TH1D* 				histoCocktailPhiGamma				= (TH1D*) directoryPCMGamma0020->Get("CocktailPhiGamma");
		TH1D* 				histoCocktailEtaPGamma				= (TH1D*) directoryPCMGamma0020->Get("CocktailEtapGamma");
		TH1D* 				histoCocktailRhoGamma				= (TH1D*) directoryPCMGamma0020->Get("CocktailRhoGamma");
		TH1D* 				histoCocktailSigmaGamma				= (TH1D*) directoryPCMGamma0020->Get("CocktailSigmaGamma");
		
		TGraphAsymmErrors* 	graphPCMDirGammaStatErr0020 		= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSysErr0020			= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSumErrAr0020		= NULL;
		graphPCMDirGammaStatErr0020								= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("graphDirGammaSpectrumStat");
		graphPCMDirGammaSysErr0020 								= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("graphDirGammaSpectrumSyst");
		graphPCMDirGammaSumErrAr0020							= (TGraphAsymmErrors*) directoryPCMGamma0020->Get("graphDirGammaSpectrumSummedAr");
	//________________________________________________ Load PCM 20-40% __________________________________________________________________________	
	TDirectory* directoryPCMGamma2040 		= (TDirectory*)filePCMGamma->Get("Gamma_PbPb_2.76TeV_20-40%"); 
// 		TH1D * 				histoPCMDRStatErr2040				= (TH1D*) directoryPCMGamma2040->Get("DoubleRatioStatError");
// 		TGraphAsymmErrors* 	graphPCMDRSysErr2040 				= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystError");
// 		TGraphAsymmErrors* 	graphPCMDRSysAErr2040 				= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystErrorA");
// 		TGraphAsymmErrors* 	graphPCMDRSysBErr2040 				= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystErrorB");
// 		TGraphAsymmErrors* 	graphPCMDRSysCErr2040 				= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioSystErrorC");
		TH1D * 				histoPCMDRPi0FitStatErr2040 		= (TH1D*) directoryPCMGamma2040->Get("DoubleRatioPi0FitStatError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysErr2040 			= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysAErr2040			= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystErrorA");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysBErr2040			= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystErrorB");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysCErr2040			= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("DoubleRatioPi0FitSystErrorC");
// 		TH1D * 				histoPCMIncRStatErr2040				= (TH1D*) directoryPCMGamma2040->Get("IncRatioStatError");
// 		TGraphAsymmErrors* 	graphPCMIncRSysErr2040 				= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystError");
// 		TGraphAsymmErrors* 	graphPCMIncRSysAErr2040				= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystErrorA");
// 		TGraphAsymmErrors* 	graphPCMIncRSysBErr2040				= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystErrorB");
// 		TGraphAsymmErrors* 	graphPCMIncRSysCErr2040				= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioSystErrorC");
// 		TH1D * 				histoPCMIncRPi0FitStatErr2040		= (TH1D*) directoryPCMGamma2040->Get("IncRatioPi0FitStatError");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysErr2040		= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystError");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysAErr2040		= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystErrorA");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysBErr2040		= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystErrorB");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysCErr2040		= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncRatioPi0FitSystErrorC");
		TH1D * 				histoPCMIncGammaStatErr2040			= (TH1D*) directoryPCMGamma2040->Get("IncGammaStatError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysErr2040			= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysAErr2040			= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystErrorA");
		TGraphAsymmErrors* 	graphPCMIncGammaSysBErr2040			= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystErrorB");
		TGraphAsymmErrors* 	graphPCMIncGammaSysCErr2040			= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("IncGammaSystErrorC");
// 		TH1D * 				histoPCMPi0StatErr2040				= (TH1D*) directoryPCMGamma2040->Get("Pi0StatError");
// 		TGraphAsymmErrors* 	graphPCMPi0SysErr2040				= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("Pi0SystError");
// 		TH1D * 				histoPCMPi0FitStatErr2040			= (TH1D*) directoryPCMGamma2040->Get("Pi0FitStatError");
// 		TGraphAsymmErrors* 	graphPCMPi0FitSysErr2040			= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("Pi0FitSystError");
// 
		TGraphAsymmErrors* 	graphPCMDirGammaStatErr2040 		= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSysErr2040			= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSumErrAr2040		= NULL;
		graphPCMDirGammaStatErr2040								= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("graphDirGammaSpectrumStat");
		graphPCMDirGammaSysErr2040 								= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("graphDirGammaSpectrumSyst");
		graphPCMDirGammaSumErrAr2040							= (TGraphAsymmErrors*) directoryPCMGamma2040->Get("graphDirGammaSpectrumSummedAr");
	//________________________________________________ Load PCM 40-80% __________________________________________________________________________
	TDirectory* directoryPCMGamma4080 		= (TDirectory*)filePCMGamma->Get("Gamma_PbPb_2.76TeV_40-80%"); 
// 		TH1D * 				histoPCMDRStatErr4080				= (TH1D*) directoryPCMGamma4080->Get("DoubleRatioStatError");
// 		TGraphAsymmErrors* 	graphPCMDRSysErr4080 				= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("DoubleRatioSystError");
// 		TGraphAsymmErrors* 	graphPCMDRSysAErr4080 				= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("DoubleRatioSystErrorA");
// 		TGraphAsymmErrors* 	graphPCMDRSysBErr4080 				= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("DoubleRatioSystErrorB");
// 		TGraphAsymmErrors* 	graphPCMDRSysCErr4080 				= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("DoubleRatioSystErrorC");
		TH1D * 				histoPCMDRPi0FitStatErr4080 		= (TH1D*) directoryPCMGamma4080->Get("DoubleRatioPi0FitStatError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysErr4080 			= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("DoubleRatioPi0FitSystError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysAErr4080			= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("DoubleRatioPi0FitSystErrorA");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysBErr4080			= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("DoubleRatioPi0FitSystErrorB");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysCErr4080			= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("DoubleRatioPi0FitSystErrorC");
// 		TH1D * 				histoPCMIncRStatErr4080				= (TH1D*) directoryPCMGamma4080->Get("IncRatioStatError");
// 		TGraphAsymmErrors* 	graphPCMIncRSysErr4080 				= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncRatioSystError");
// 		TGraphAsymmErrors* 	graphPCMIncRSysAErr4080				= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncRatioSystErrorA");
// 		TGraphAsymmErrors* 	graphPCMIncRSysBErr4080				= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncRatioSystErrorB");
// 		TGraphAsymmErrors* 	graphPCMIncRSysCErr4080				= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncRatioSystErrorC");
// 		TH1D * 				histoPCMIncRPi0FitStatErr4080		= (TH1D*) directoryPCMGamma4080->Get("IncRatioPi0FitStatError");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysErr4080		= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncRatioPi0FitSystError");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysAErr4080		= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncRatioPi0FitSystErrorA");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysBErr4080		= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncRatioPi0FitSystErrorB");
// 		TGraphAsymmErrors* 	graphPCMIncRFitPi0SysCErr4080		= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncRatioPi0FitSystErrorC");
		TH1D * 				histoPCMIncGammaStatErr4080			= (TH1D*) directoryPCMGamma4080->Get("IncGammaStatError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysErr4080			= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncGammaSystError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysAErr4080			= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncGammaSystErrorA");
		TGraphAsymmErrors* 	graphPCMIncGammaSysBErr4080			= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncGammaSystErrorB");
		TGraphAsymmErrors* 	graphPCMIncGammaSysCErr4080			= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("IncGammaSystErrorC");
// 		TH1D * 				histoPCMPi0StatErr4080				= (TH1D*) directoryPCMGamma4080->Get("Pi0StatError");
// 		TGraphAsymmErrors* 	graphPCMPi0SysErr4080				= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("Pi0SystError");
// 		TH1D * 				histoPCMPi0FitStatErr4080			= (TH1D*) directoryPCMGamma4080->Get("Pi0FitStatError");
// 		TGraphAsymmErrors* 	graphPCMPi0FitSysErr4080			= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("Pi0FitSystError");
// 
		TGraphAsymmErrors* 	graphPCMDirGammaStatErr4080 		= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSysErr4080			= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSumErrAr4080		= NULL;
		graphPCMDirGammaStatErr4080								= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("graphDirGammaSpectrumStat");
		graphPCMDirGammaSysErr4080 								= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("graphDirGammaSpectrumSyst");
		graphPCMDirGammaSumErrAr4080							= (TGraphAsymmErrors*) directoryPCMGamma4080->Get("graphDirGammaSpectrumSummedAr");

	TDirectory* directoryTheory 		= (TDirectory*)filePCMGamma->Get("Theory"); 
		TGraphAsymmErrors* graphTheoryNLODR0020 				= (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatio_0-20%");
		TGraphAsymmErrors* graphTheoryNLODR2040 				= (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatio_20-40%");
		TGraphAsymmErrors* graphTheoryNLODR4080 				= (TGraphAsymmErrors*) directoryTheory->Get("NLODoubleRatio_40-80%");
		TGraphAsymmErrors* graphTheoryNLO0020 					= (TGraphAsymmErrors*) directoryTheory->Get("NLO_0-20%");
		TGraphAsymmErrors* graphTheoryNLO2040 					= (TGraphAsymmErrors*) directoryTheory->Get("NLO_20-40%");
		TGraphAsymmErrors* graphTheoryNLO4080 					= (TGraphAsymmErrors*) directoryTheory->Get("NLO_40-80%");
		TGraphAsymmErrors* graphTheoryEPS090020 				= (TGraphAsymmErrors*) directoryTheory->Get("EPS09_0-20%");
		TGraphAsymmErrors* graphTheoryRelErrEPS090020			= CalculateRelErrAsymmGraphAround1(graphTheoryEPS090020, "graphTheoryRelErrEPS090020");
		TGraphAsymmErrors* graphTheoryEPS092040 				= (TGraphAsymmErrors*) directoryTheory->Get("EPS09_20-40%");
		TGraphAsymmErrors* graphTheoryRelErrEPS092040			= CalculateRelErrAsymmGraphAround1(graphTheoryEPS092040, "graphTheoryRelErrEPS092040");
		TGraphAsymmErrors* graphTheoryEPS094080 				= (TGraphAsymmErrors*) directoryTheory->Get("EPS09_40-80%");
		TGraphAsymmErrors* graphTheoryRelErrEPS094080			= CalculateRelErrAsymmGraphAround1(graphTheoryEPS094080, "graphTheoryRelErrEPS094080");
		TGraphAsymmErrors* graphTheoryEPS09DR0020 				= (TGraphAsymmErrors*) directoryTheory->Get("EPS09DoubleRatio_0-20%");
		TGraphAsymmErrors* graphTheoryEPS09DR2040 				= (TGraphAsymmErrors*) directoryTheory->Get("EPS09DoubleRatio_20-40%");
		TGraphAsymmErrors* graphTheoryEPS09DR4080 				= (TGraphAsymmErrors*) directoryTheory->Get("EPS09DoubleRatio_40-80%");
		TGraphAsymmErrors* graphTheoryCT100020 					= (TGraphAsymmErrors*) directoryTheory->Get("CT10BF_0-20%");
		TGraphAsymmErrors* graphTheoryCT102040 					= (TGraphAsymmErrors*) directoryTheory->Get("CT10BF_20-40%");
		TGraphAsymmErrors* graphTheoryCT104080 					= (TGraphAsymmErrors*) directoryTheory->Get("CT10BF_40-80%");
		TGraphAsymmErrors* graphTheoryCT10DR0020 				= (TGraphAsymmErrors*) directoryTheory->Get("CT10BFDoubleRatio_0-20%");
		TGraphAsymmErrors* graphTheoryCT10DR2040 				= (TGraphAsymmErrors*) directoryTheory->Get("CT10BFDoubleRatio_20-40%");
		TGraphAsymmErrors* graphTheoryCT10DR4080 				= (TGraphAsymmErrors*) directoryTheory->Get("CT10BFDoubleRatio_40-80%");
		
		
		TF1* fitTheoryGammaEPSO90020	 		= FitObject("l","fitTheoryGammaEPSO90020","Photon",graphTheoryEPS090020,1.6,12.5,NULL,"QNRMEX0+");
		cout << WriteParameterToFile(fitTheoryGammaEPSO90020)<< endl;   
		fitTheoryGammaEPSO90020->SetRange(0.9,14);
		TF1* fitTheoryGammaEPSO92040	 		= FitObject("l","fitTheoryGammaEPSO90020","Photon",graphTheoryEPS092040,1.6,12.5,NULL,"QNRMEX0+");
		cout << WriteParameterToFile(fitTheoryGammaEPSO92040)<< endl;   
		fitTheoryGammaEPSO92040->SetRange(0.9,14);
		TF1* fitTheoryGammaEPSO94080	 		= FitObject("l","fitTheoryGammaEPSO90020","Photon",graphTheoryEPS094080,1.6,12.5,NULL,"QNRMEX0+");
		cout << WriteParameterToFile(fitTheoryGammaEPSO94080)<< endl;   
		fitTheoryGammaEPSO94080->SetRange(0.9,14);


	//*******************************************************************************************************************************************
	//*********************************************** Load PCM histograms from pp PCM file ******************************************************
	//*******************************************************************************************************************************************
	TFile*				filePCMGammapp 						= new TFile( inputFileNamePCMpp.Data());
	//________________________________________________ Load PCM pp 7TeV _________________________________________________________________________
	TDirectory* 		directoryPCMGammapp7TeV 				= (TDirectory*)filePCMGammapp->Get("Gamma_7TeV_pp"); 
		TH1D * 				histoPCMDRPi0FitStatErrpp7TeV 		= (TH1D*) directoryPCMGammapp7TeV->Get("DoubleRatioPi0FitStatError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysErrpp7TeV 		= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("DoubleRatioPi0FitSystError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysAErrpp7TeV		= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("DoubleRatioPi0FitSystErrorA");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysBErrpp7TeV		= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("DoubleRatioPi0FitSystErrorB");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysCErrpp7TeV		= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("DoubleRatioPi0FitSystErrorC");
		TH1D * 				histoPCMIncGammaStatErrpp7TeV		= (TH1D*) directoryPCMGammapp7TeV->Get("IncGammaStatError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysErrpp7TeV		= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("IncGammaSystError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysAErrpp7TeV		= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("IncGammaSystErrorA");
		TGraphAsymmErrors* 	graphPCMIncGammaSysBErrpp7TeV		= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("IncGammaSystErrorB");
		TGraphAsymmErrors* 	graphPCMIncGammaSysCErrpp7TeV		= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("IncGammaSystErrorC");
		TGraphAsymmErrors* 	graphPCMDirGammaStatErrpp7TeV 		= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSysErrpp7TeV		= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSumErrArpp7TeV		= NULL;
		graphPCMDirGammaStatErrpp7TeV							= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("graphDirGammaSpectrumStat");
		graphPCMDirGammaSysErrpp7TeV 							= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("graphDirGammaSpectrumSyst");
		graphPCMDirGammaSumErrArpp7TeV							= (TGraphAsymmErrors*) directoryPCMGammapp7TeV->Get("graphDirGammaSpectrumSummedAr");
	//________________________________________________ Load PCM pp 2.76TeV _________________________________________________________________________	
	TDirectory* 		directoryPCMGammapp2760GeV 				= (TDirectory*)filePCMGammapp->Get("Gamma_2.76TeV_pp"); 
		TH1D * 				histoPCMDRPi0FitStatErrpp2760GeV 	= (TH1D*) directoryPCMGammapp2760GeV->Get("DoubleRatioPi0FitStatError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysErrpp2760GeV 	= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("DoubleRatioPi0FitSystError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysAErrpp2760GeV	= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("DoubleRatioPi0FitSystErrorA");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysBErrpp2760GeV	= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("DoubleRatioPi0FitSystErrorB");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysCErrpp2760GeV	= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("DoubleRatioPi0FitSystErrorC");
		TH1D * 				histoPCMIncGammaStatErrpp2760GeV	= (TH1D*) directoryPCMGammapp2760GeV->Get("IncGammaStatError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysErrpp2760GeV		= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("IncGammaSystError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysAErrpp2760GeV	= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("IncGammaSystErrorA");
		TGraphAsymmErrors* 	graphPCMIncGammaSysBErrpp2760GeV	= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("IncGammaSystErrorB");
		TGraphAsymmErrors* 	graphPCMIncGammaSysCErrpp2760GeV	= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("IncGammaSystErrorC");
		TGraphAsymmErrors* 	graphPCMDirGammaStatErrpp2760GeV 	= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSysErrpp2760GeV		= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSumErrArpp2760GeV	= NULL;
		graphPCMDirGammaStatErrpp2760GeV						= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("graphDirGammaSpectrumStat");
		graphPCMDirGammaSysErrpp2760GeV 						= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("graphDirGammaSpectrumSyst");
		graphPCMDirGammaSumErrArpp2760GeV						= (TGraphAsymmErrors*) directoryPCMGammapp2760GeV->Get("graphDirGammaSpectrumSummedAr");
	TDirectory* directoryTheorypp 		= (TDirectory*)filePCMGammapp->Get("Theory"); 
		TGraphAsymmErrors* graphTheoryNLODRpp2760GeV 			= (TGraphAsymmErrors*) directoryTheorypp->Get("NLODoubleRatio_2.76TeV_pp");
		TGraphAsymmErrors* graphTheoryNLODRpp7TeV 				= (TGraphAsymmErrors*) directoryTheorypp->Get("NLODoubleRatio_7TeV_pp");
		TGraphAsymmErrors* graphTheoryNLOpp2760GeV 				= (TGraphAsymmErrors*) directoryTheorypp->Get("NLO_2.76TeV_pp");
		TGraphAsymmErrors* graphTheoryNLOpp7TeV 				= (TGraphAsymmErrors*) directoryTheorypp->Get("NLO_7TeV_pp");
		

	//*******************************************************************************************************************************************
	//*********************************************** Load PCM histograms from pPb PCM file ******************************************************
	//*******************************************************************************************************************************************
	TFile*				filePCMGammapPb 						= new TFile( inputFileNamePCMpPb.Data());
	//________________________________________________ Load PCM pPb 5.023TeV _________________________________________________________________________
	TDirectory* 		directoryPCMGammapPb 					= (TDirectory*)filePCMGammapPb->Get("Gamma_pPb_5.023TeV_0-100%"); 
		TH1D * 				histoPCMDRPi0FitStatErrpPb 			= (TH1D*) directoryPCMGammapPb->Get("DoubleRatioPi0FitStatError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysErrpPb 			= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("DoubleRatioPi0FitSystError");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysAErrpPb			= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("DoubleRatioPi0FitSystErrorA");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysBErrpPb			= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("DoubleRatioPi0FitSystErrorB");
		TGraphAsymmErrors* 	graphPCMDRPi0FitSysCErrpPb			= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("DoubleRatioPi0FitSystErrorC");
		TH1D * 				histoPCMIncGammaStatErrpPb				= (TH1D*) directoryPCMGammapPb->Get("IncGammaStatError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysErrpPb			= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncGammaSystError");
		TGraphAsymmErrors* 	graphPCMIncGammaSysAErrpPb			= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncGammaSystErrorA");
		TGraphAsymmErrors* 	graphPCMIncGammaSysBErrpPb			= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncGammaSystErrorB");
		TGraphAsymmErrors* 	graphPCMIncGammaSysCErrpPb			= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("IncGammaSystErrorC");
		TGraphAsymmErrors* 	graphPCMDirGammaStatErrpPb 			= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSysErrpPb			= NULL;
		TGraphAsymmErrors* 	graphPCMDirGammaSumErrArpPb			= NULL;
		graphPCMDirGammaStatErrpPb								= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("graphDirGammaSpectrumStat");
		graphPCMDirGammaSysErrpPb 								= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("graphDirGammaSpectrumSyst");
		graphPCMDirGammaSumErrArpPb								= (TGraphAsymmErrors*) directoryPCMGammapPb->Get("graphDirGammaSpectrumSummedAr");
	TDirectory* directoryTheorypPb 								= (TDirectory*)filePCMGammapPb->Get("Theory"); 
		TGraphAsymmErrors* graphTheoryNLODRpPb 					= (TGraphAsymmErrors*) directoryTheorypPb->Get("NLODoubleRatio_0-100%");
		TGraphAsymmErrors* graphTheoryNLOpPb 					= (TGraphAsymmErrors*) directoryTheorypPb->Get("NLO_0-100%");

		
	//*******************************************************************************************************************************************
	//*********************************************** Load PHOS histograms from PHOS file *******************************************************
	//*******************************************************************************************************************************************		
	TFile*				filePHOSGamma 						= new TFile( inputFileNamePHOS.Data());
	//________________________________________________ Load PHOS 0-20% __________________________________________________________________________
	TDirectory* 		directoryPHOSGamma0020 					= (TDirectory*)filePHOSGamma->Get("PHOS_PbPb_2760_Centrality_00-20"); 
		TH1D * 				histoPHOSDRPi0FitStatErr0020 		= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_DoubleRatio_PbPb_cen00-20_Stat");
		TH1D * 				histoPHOSDRPi0FitSysErr0020 		= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_DoubleRatio_PbPb_cen00-20_Syst");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysErr0020 		= new TGraphAsymmErrors(histoPHOSDRPi0FitSysErr0020);	
		TH1D * 				histoPHOSDRPi0FitSysAErr0020		= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_DoubleRatio_PbPb_cen00-20_SystA");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysAErr0020 		= new TGraphAsymmErrors(histoPHOSDRPi0FitSysAErr0020);	
		cout << "here" << endl;
		TH1D * 				histoPHOSDRPi0FitSysBcontErr0020	= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_DoubleRatio_PbPb_cen00-20_SystBcont");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBcontErr0020 	= new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcontErr0020);	
		cout << "B Syst. Cont" << endl;
		graphPHOSDRPi0FitSysBcontErr0020->Print();
		TH1D * 				histoPHOSDRPi0FitSysBpi0Err0020		= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_DoubleRatio_PbPb_cen00-20_SystBpi0");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBpi0Err0020 	= new TGraphAsymmErrors(histoPHOSDRPi0FitSysBpi0Err0020);	
		cout << "B Sys pi0" << endl;
		graphPHOSDRPi0FitSysBpi0Err0020->Print();
		TH1D * 				histoPHOSDRPi0FitSysBcocktErr0020	= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_DoubleRatio_PbPb_cen00-20_SystBcockt");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBcocktErr0020 	= new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcocktErr0020);	
		cout << "B Sys cockt" << endl;
		graphPHOSDRPi0FitSysBcocktErr0020->Print();
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBErr0020 		= Add3ErrorsOfGraphsQuadratically (graphPHOSDRPi0FitSysBcontErr0020, graphPHOSDRPi0FitSysBpi0Err0020, graphPHOSDRPi0FitSysBcocktErr0020);
		graphPHOSDRPi0FitSysBErr0020->SetName("hPHOS_DoubleRatio_PbPb_cen00-20_SystBtot");
		
		cout << "here" << endl;
		TH1D * 				histoPHOSDRPi0FitSysCErr0020		= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_DoubleRatio_PbPb_cen00-20_SystC");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysCErr0020 		= new TGraphAsymmErrors(histoPHOSDRPi0FitSysCErr0020);	
		cout << "here" << endl;
// 		TH1D * 				histoPHOSPi0StatErr0020				= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_pi0_PbPb_cen00-20_Stat");
// 		TH1D *			 	histoPHOSPi0SysErr0020				= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_pi0_PbPb_cen00-20_Syst");
// 		TH1D *			 	histoPHOSPi0SysAErr0020				= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_pi0_PbPb_cen00-20_SystA");
// 		TH1D *			 	histoPHOSPi0SysBErr0020				= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_pi0_PbPb_cen00-20_SystB");
// 		TH1D *			 	histoPHOSPi0SysCErr0020				= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_pi0_PbPb_cen00-20_SystC");
		
		TH1D * 				histoPHOSIncGammaStatErr0020		= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_gammaIncl_PbPb_cen00-20_Stat");
		TH1D *			 	histoPHOSIncGammaSysErr0020			= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_gammaIncl_PbPb_cen00-20_Syst");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysErr0020 		= new TGraphAsymmErrors(histoPHOSIncGammaSysErr0020);	
		TH1D *			 	histoPHOSIncGammaSysAErr0020		= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_gammaIncl_PbPb_cen00-20_SystA");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysAErr0020 		= new TGraphAsymmErrors(histoPHOSIncGammaSysAErr0020);	
		TH1D *			 	histoPHOSIncGammaSysBErr0020		= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_gammaIncl_PbPb_cen00-20_SystBtot");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysBErr0020 		= new TGraphAsymmErrors(histoPHOSIncGammaSysBErr0020);	
		TH1D *			 	histoPHOSIncGammaSysCErr0020		= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_gammaIncl_PbPb_cen00-20_SystC");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysCErr0020 		= new TGraphAsymmErrors(histoPHOSIncGammaSysCErr0020);	
		
		TH1D* 				histoPHOSDirGammaStatErr0020 		= NULL;
		TH1D* 				histoPHOSDirGammaSysErr0020			= NULL;
// 		TH1D* 				histoPHOSDirGammaSysBErrnl0020		= NULL;
// 		TH1D* 				histoPHOSDirGammaSysBErrglobalE0020	= NULL;
// 		TH1D* 				histoPHOSDirGammaSysCErr0020		= NULL;
		histoPHOSDirGammaStatErr0020							= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_gammaDir_PbPb_cen00-20_Stat");
		histoPHOSDirGammaSysErr0020 							= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_gammaDir_PbPb_cen00-20_Syst");
// 		histoPHOSDirGammaSysBErrnl0020							= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_gammaDir_PbPb_cen00-20_SystBnl");
// 		histoPHOSDirGammaSysBErrglobalE0020						= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_gammaDir_PbPb_cen00-20_SystBglobalE");
// 		histoPHOSDirGammaSysCErr0020							= (TH1D*) directoryPHOSGamma0020->Get("hPHOS_gammaDir_PbPb_cen00-20_SystC");
	//________________________________________________ Load PHOS 20-40% _________________________________________________________________________
	TDirectory* 		directoryPHOSGamma2040 					= (TDirectory*)filePHOSGamma->Get("PHOS_PbPb_2760_Centrality_20-40"); 
		TH1D * 				histoPHOSDRPi0FitStatErr2040 		= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_Stat");
		TH1D * 				histoPHOSDRPi0FitSysErr2040 		= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_Syst");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysErr2040 		= new TGraphAsymmErrors(histoPHOSDRPi0FitSysErr2040);	
		TH1D * 				histoPHOSDRPi0FitSysAErr2040		= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystA");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysAErr2040 		= new TGraphAsymmErrors(histoPHOSDRPi0FitSysAErr2040);	
		TH1D * 				histoPHOSDRPi0FitSysBcontErr2040	= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystBcont");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBcontErr2040 	= new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcontErr2040);	
		TH1D * 				histoPHOSDRPi0FitSysBpi0Err2040		= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystBpi0");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBpi0Err2040 	= new TGraphAsymmErrors(histoPHOSDRPi0FitSysBpi0Err2040);	
		TH1D * 				histoPHOSDRPi0FitSysBcocktErr2040	= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystBcockt");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBcocktErr2040 	= new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcocktErr2040);	
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBErr2040 		= Add3ErrorsOfGraphsQuadratically (graphPHOSDRPi0FitSysBcontErr2040, graphPHOSDRPi0FitSysBpi0Err2040, graphPHOSDRPi0FitSysBcocktErr2040);
		graphPHOSDRPi0FitSysBErr2040->SetName("hPHOS_DoubleRatio_PbPb_cen20-40_SystBtot");
		TH1D * 				histoPHOSDRPi0FitSysCErr2040		= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_DoubleRatio_PbPb_cen20-40_SystC");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysCErr2040 		= new TGraphAsymmErrors(histoPHOSDRPi0FitSysCErr2040);	
	
// 		TH1D * 				histoPHOSPi0StatErr2040				= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_pi0_PbPb_cen20-40_Stat");
// 		TH1D *			 	histoPHOSPi0SysErr2040				= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_pi0_PbPb_cen20-40_Syst");
// 		TH1D *			 	histoPHOSPi0SysAErr2040				= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_pi0_PbPb_cen20-40_SystA");
// 		TH1D *			 	histoPHOSPi0SysBErr2040				= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_pi0_PbPb_cen20-40_SystB");
// 		TH1D *			 	histoPHOSPi0SysCErr2040				= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_pi0_PbPb_cen20-40_SystC");
		
		TH1D * 				histoPHOSIncGammaStatErr2040		= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaIncl_PbPb_cen20-40_Stat");
		TH1D *			 	histoPHOSIncGammaSysErr2040			= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaIncl_PbPb_cen20-40_Syst");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysErr2040 		= new TGraphAsymmErrors(histoPHOSIncGammaSysErr2040);	
		TH1D *			 	histoPHOSIncGammaSysAErr2040		= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaIncl_PbPb_cen20-40_SystA");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysAErr2040 		= new TGraphAsymmErrors(histoPHOSIncGammaSysAErr2040);	
		TH1D *			 	histoPHOSIncGammaSysBErr2040		= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaIncl_PbPb_cen20-40_SystBtot");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysBErr2040 		= new TGraphAsymmErrors(histoPHOSIncGammaSysBErr2040);	
		TH1D *			 	histoPHOSIncGammaSysCErr2040		= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaIncl_PbPb_cen20-40_SystC");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysCErr2040 		= new TGraphAsymmErrors(histoPHOSIncGammaSysCErr2040);	
 	
		TH1D* 				histoPHOSDirGammaStatErr2040 		= NULL;
		TH1D* 				histoPHOSDirGammaSysErr2040			= NULL;
// 		TH1D* 				histoPHOSDirGammaSysBErrnl2040		= NULL;
// 		TH1D* 				histoPHOSDirGammaSysBErrglobalE2040	= NULL;
// 		TH1D* 				histoPHOSDirGammaSysCErr2040		= NULL;
		histoPHOSDirGammaStatErr2040							= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaDir_PbPb_cen20-40_Stat");
		histoPHOSDirGammaSysErr2040 							= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaDir_PbPb_cen20-40_Syst");
// 		histoPHOSDirGammaSysBErrnl2040							= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaDir_PbPb_cen20-40_SystBnl");
// 		histoPHOSDirGammaSysBErrglobalE2040						= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaDir_PbPb_cen20-40_SystBglobalE");
// 		histoPHOSDirGammaSysCErr2040							= (TH1D*) directoryPHOSGamma2040->Get("hPHOS_gammaDir_PbPb_cen20-40_SystC");
	//________________________________________________ Load PHOS 40-80% _________________________________________________________________________
	TDirectory* 		directoryPHOSGamma4080 					= (TDirectory*)filePHOSGamma->Get("PHOS_PbPb_2760_Centrality_40-80"); 
		TH1D * 				histoPHOSDRPi0FitStatErr4080 		= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_DoubleRatio_PbPb_cen40-80_Stat");
		TH1D * 				histoPHOSDRPi0FitSysErr4080 		= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_DoubleRatio_PbPb_cen40-80_Syst");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysErr4080 		= new TGraphAsymmErrors(histoPHOSDRPi0FitSysErr4080);	
		TH1D * 				histoPHOSDRPi0FitSysAErr4080		= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_DoubleRatio_PbPb_cen40-80_SystA");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysAErr4080 		= new TGraphAsymmErrors(histoPHOSDRPi0FitSysAErr4080);	
		TH1D * 				histoPHOSDRPi0FitSysBcontErr4080	= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_DoubleRatio_PbPb_cen40-80_SystBcont");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBcontErr4080 	= new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcontErr4080);	
		TH1D * 				histoPHOSDRPi0FitSysBpi0Err4080		= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_DoubleRatio_PbPb_cen40-80_SystBpi0");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBpi0Err4080 	= new TGraphAsymmErrors(histoPHOSDRPi0FitSysBpi0Err4080);	
		TH1D * 				histoPHOSDRPi0FitSysBcocktErr4080	= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_DoubleRatio_PbPb_cen40-80_SystBcockt");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBcocktErr4080 	= new TGraphAsymmErrors(histoPHOSDRPi0FitSysBcocktErr4080);	
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysBErr4080 		= Add3ErrorsOfGraphsQuadratically (graphPHOSDRPi0FitSysBcontErr4080, graphPHOSDRPi0FitSysBpi0Err4080, graphPHOSDRPi0FitSysBcocktErr4080);
		graphPHOSDRPi0FitSysBErr4080->SetName("hPHOS_DoubleRatio_PbPb_cen40-80_SystBtot");
		TH1D * 				histoPHOSDRPi0FitSysCErr4080		= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_DoubleRatio_PbPb_cen40-80_SystC");
		TGraphAsymmErrors* 	graphPHOSDRPi0FitSysCErr4080 		= new TGraphAsymmErrors(histoPHOSDRPi0FitSysCErr4080);	

// 		TH1D * 				histoPHOSPi0StatErr4080				= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_pi0_PbPb_cen40-80_Stat");
// 		TH1D *			 	histoPHOSPi0SysErr4080				= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_pi0_PbPb_cen40-80_Syst");
// 		TH1D *			 	histoPHOSPi0SysAErr4080				= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_pi0_PbPb_cen40-80_SystA");
// 		TH1D *			 	histoPHOSPi0SysBErr4080				= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_pi0_PbPb_cen40-80_SystB");
// 		TH1D *			 	histoPHOSPi0SysCErr4080				= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_pi0_PbPb_cen40-80_SystC");
		
		TH1D * 				histoPHOSIncGammaStatErr4080		= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_gammaIncl_PbPb_cen40-80_Stat");
		TH1D *			 	histoPHOSIncGammaSysErr4080			= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_gammaIncl_PbPb_cen40-80_Syst");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysErr4080 		= new TGraphAsymmErrors(histoPHOSIncGammaSysErr4080);	
		TH1D *			 	histoPHOSIncGammaSysAErr4080		= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_gammaIncl_PbPb_cen40-80_SystA");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysAErr4080 		= new TGraphAsymmErrors(histoPHOSIncGammaSysAErr4080);	
		TH1D *			 	histoPHOSIncGammaSysBErr4080		= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_gammaIncl_PbPb_cen40-80_SystBtot");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysBErr4080 		= new TGraphAsymmErrors(histoPHOSIncGammaSysBErr4080);	
		TH1D *			 	histoPHOSIncGammaSysCErr4080		= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_gammaIncl_PbPb_cen40-80_SystC");
		TGraphAsymmErrors* 	graphPHOSIncGammaSysCErr4080 		= new TGraphAsymmErrors(histoPHOSIncGammaSysCErr4080);	
// 
		TH1D* 				histoPHOSDirGammaStatErr4080 		= NULL;
		TH1D* 				histoPHOSDirGammaSysErr4080			= NULL;
// 		TH1D* 				histoPHOSDirGammaSysBErrnl4080		= NULL;
// 		TH1D* 				histoPHOSDirGammaSysBErrglobalE4080	= NULL;
// 		TH1D* 				histoPHOSDirGammaSysCErr4080		= NULL;
		histoPHOSDirGammaStatErr4080							= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_gammaDir_PbPb_cen40-80_Stat");
		histoPHOSDirGammaSysErr4080 							= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_gammaDir_PbPb_cen40-80_Syst");
// 		histoPHOSDirGammaSysBErrnl4080							= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_gammaDir_PbPb_cen40-80_SystBnl");
// 		histoPHOSDirGammaSysBErrglobalE4080						= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_gammaDir_PbPb_cen40-80_SystBglobalE");
// 		histoPHOSDirGammaSysCErr4080							= (TH1D*) directoryPHOSGamma4080->Get("hPHOS_gammaDir_PbPb_cen40-80_SystC");
// 	

	//*******************************************************************************************************************************************
	//***************************************************** Combine DR of PCM and PHOS **********************************************************
	//*******************************************************************************************************************************************		
	Double_t newBinsComb[21] 									= {	0.9, 1.1, 1.3, 1.5, 1.7, 
																	1.9, 2.1, 2.3, 2.5, 2.7, 
																	3.0, 3.3, 3.7, 4.1, 4.6, 
																	5.4, 6.2, 7.0, 8.0, 11.0, 
																	14.0};
	TGraphAsymmErrors *graphCombDRPi0FitSysErr0020;
	TGraphAsymmErrors *graphCombDRPi0FitSysAErr0020;
	TGraphAsymmErrors *graphCombDRPi0FitSysBErr0020;
	TGraphAsymmErrors *graphCombDRPi0FitSysCErr0020;
	TGraphAsymmErrors *graphCombDRPi0FitStatErr0020;
	TGraphAsymmErrors *graphCombDRPi0FitSumErr0020;		
	graphCombDRPi0FitSumErr0020 	= CombinePtPointsSpectraAdv(	histoPCMDRPi0FitStatErr0020, graphPCMDRPi0FitSysErr0020, 
																graphPCMDRPi0FitSysAErr0020, graphPCMDRPi0FitSysBErr0020, graphPCMDRPi0FitSysCErr0020,
																histoPHOSDRPi0FitStatErr0020, graphPHOSDRPi0FitSysErr0020 ,
																graphPHOSDRPi0FitSysAErr0020, graphPHOSDRPi0FitSysBErr0020, graphPHOSDRPi0FitSysCErr0020,
																graphCombDRPi0FitStatErr0020, graphCombDRPi0FitSysErr0020, 
																graphCombDRPi0FitSysAErr0020, graphCombDRPi0FitSysBErr0020, graphCombDRPi0FitSysCErr0020,
																newBinsComb, 21, 0, 0, -1);
	graphCombDRPi0FitSumErr0020->Print();
	TGraphAsymmErrors* graphCombDRPi0FitStatSysAErr0020 		= AddErrorsOfGraphsQuadratically (graphCombDRPi0FitStatErr0020, graphCombDRPi0FitSysAErr0020);
	Double_t SysCCombDRPi0Fit0020								= graphCombDRPi0FitSysCErr0020->GetErrorYlow(4)/graphCombDRPi0FitSysCErr0020->GetY()[4];



	TGraphAsymmErrors *graphCombDRPi0FitSysErr2040;
	TGraphAsymmErrors *graphCombDRPi0FitSysAErr2040;
	TGraphAsymmErrors *graphCombDRPi0FitSysBErr2040;
	TGraphAsymmErrors *graphCombDRPi0FitSysCErr2040;
	TGraphAsymmErrors *graphCombDRPi0FitStatErr2040;
	TGraphAsymmErrors *graphCombDRPi0FitSumErr2040;		
	graphCombDRPi0FitSumErr2040 	= CombinePtPointsSpectraAdv(	histoPCMDRPi0FitStatErr2040, graphPCMDRPi0FitSysErr2040, 
																graphPCMDRPi0FitSysAErr2040, graphPCMDRPi0FitSysBErr2040, graphPCMDRPi0FitSysCErr2040,
																histoPHOSDRPi0FitStatErr2040, graphPHOSDRPi0FitSysErr2040 ,
																graphPHOSDRPi0FitSysAErr2040, graphPHOSDRPi0FitSysBErr2040, graphPHOSDRPi0FitSysCErr2040,
																graphCombDRPi0FitStatErr2040, graphCombDRPi0FitSysErr2040, 
																graphCombDRPi0FitSysAErr2040, graphCombDRPi0FitSysBErr2040, graphCombDRPi0FitSysCErr2040,
																newBinsComb, 21, 0, 0, -1);
	graphCombDRPi0FitSumErr2040->Print();
	TGraphAsymmErrors* graphCombDRPi0FitStatSysAErr2040 		= AddErrorsOfGraphsQuadratically (graphCombDRPi0FitStatErr2040, graphCombDRPi0FitSysAErr2040);
	Double_t SysCCombDRPi0Fit2040								= graphCombDRPi0FitSysCErr2040->GetErrorYlow(4)/graphCombDRPi0FitSysCErr2040->GetY()[4];
	
	TGraphAsymmErrors *graphCombDRPi0FitSysErr4080;
	TGraphAsymmErrors *graphCombDRPi0FitSysAErr4080;
	TGraphAsymmErrors *graphCombDRPi0FitSysBErr4080;
	TGraphAsymmErrors *graphCombDRPi0FitSysCErr4080;
	TGraphAsymmErrors *graphCombDRPi0FitStatErr4080;
	TGraphAsymmErrors *graphCombDRPi0FitSumErr4080;		
	graphCombDRPi0FitSumErr4080 	= CombinePtPointsSpectraAdv(	histoPCMDRPi0FitStatErr4080, graphPCMDRPi0FitSysErr4080, 
																graphPCMDRPi0FitSysAErr4080, graphPCMDRPi0FitSysBErr4080, graphPCMDRPi0FitSysCErr4080,
																histoPHOSDRPi0FitStatErr4080, graphPHOSDRPi0FitSysErr4080 ,
																graphPHOSDRPi0FitSysAErr4080, graphPHOSDRPi0FitSysBErr4080, graphPHOSDRPi0FitSysCErr4080,
																graphCombDRPi0FitStatErr4080, graphCombDRPi0FitSysErr4080, 
																graphCombDRPi0FitSysAErr4080, graphCombDRPi0FitSysBErr4080, graphCombDRPi0FitSysCErr4080,
																newBinsComb, 21, 0, 0, -1);
	graphCombDRPi0FitSumErr4080->Print();
	TGraphAsymmErrors* graphCombDRPi0FitStatSysAErr4080 		= AddErrorsOfGraphsQuadratically (graphCombDRPi0FitStatErr4080, graphCombDRPi0FitSysAErr4080);
	Double_t SysCCombDRPi0Fit4080								= graphCombDRPi0FitSysCErr4080->GetErrorYlow(4)/graphCombDRPi0FitSysCErr4080->GetY()[4];

	
	//*******************************************************************************************************************************************
	//******************************************** Combine Inclusive gamma of PCM and PHOS ******************************************************
	//*******************************************************************************************************************************************		
	//__________________________________________________ 0-20% combine spectra __________________________________________________________________
	TGraphAsymmErrors *graphCombIncGammaSysErr0020;
	TGraphAsymmErrors *graphCombIncGammaSysAErr0020;
	TGraphAsymmErrors *graphCombIncGammaSysBErr0020;
	TGraphAsymmErrors *graphCombIncGammaSysCErr0020;
	TGraphAsymmErrors *graphCombIncGammaStatErr0020;
	TGraphAsymmErrors *graphCombIncGammaSumErr0020;		
	graphCombIncGammaSumErr0020 	= CombinePtPointsSpectraAdv(	histoPCMIncGammaStatErr0020, graphPCMIncGammaSysErr0020,
																graphPCMIncGammaSysAErr0020, graphPCMIncGammaSysBErr0020, graphPCMIncGammaSysCErr0020,
																histoPHOSIncGammaStatErr0020, graphPHOSIncGammaSysErr0020 ,
																graphPHOSIncGammaSysAErr0020, graphPHOSIncGammaSysBErr0020, graphPHOSIncGammaSysCErr0020,
																graphCombIncGammaStatErr0020, graphCombIncGammaSysErr0020,
																graphCombIncGammaSysAErr0020, graphCombIncGammaSysBErr0020, graphCombIncGammaSysCErr0020,
																newBinsComb, 21, 0, 0, -1);
	
	//__________________________________________ 0-20% fit combined and build ratio of individual to fit ________________________________________
	TF1* fitIncGammaCombQCD0020	 								= FitObject("qcd","fitIncGammaCombQCD0020","Photon",graphCombIncGammaSumErr0020,0.9,14,NULL,"QNRMEX0+");
	cout << WriteParameterToFile(fitIncGammaCombQCD0020)<< endl;   
	TH1D* histoFitQCDIncGammaComb0020 							= (TH1D*)fitIncGammaCombQCD0020->GetHistogram();
	
	TGraphAsymmErrors* graphPCMIncGammaStatErr0020 				= new TGraphAsymmErrors(histoPCMIncGammaStatErr0020);
	graphPCMIncGammaStatErr0020->RemovePoint(0);
	TGraphAsymmErrors* graphPCMIncGammaStatSysAErr0020 			= AddErrorsOfGraphsQuadratically (graphPCMIncGammaStatErr0020, graphPCMIncGammaSysAErr0020);
	
	TGraphAsymmErrors* graphPHOSIncGammaStatErr0020 			= new TGraphAsymmErrors(histoPHOSIncGammaStatErr0020);
	graphPHOSIncGammaStatErr0020->RemovePoint(0);
	graphPHOSIncGammaSysErr0020->RemovePoint(0);
	graphPHOSIncGammaSysAErr0020->RemovePoint(0);
	graphPHOSIncGammaSysBErr0020->RemovePoint(0);
	graphPHOSIncGammaSysCErr0020->RemovePoint(0);
	TGraphAsymmErrors* graphPHOSIncGammaStatSysAErr0020 		= AddErrorsOfGraphsQuadratically (graphPHOSIncGammaStatErr0020, graphPHOSIncGammaSysAErr0020);
	
	
	
	TGraphAsymmErrors* graphRatioCombPHOSIncGammaStatErr0020 	= (TGraphAsymmErrors*) graphPHOSIncGammaStatErr0020->Clone();	
	TGraphAsymmErrors* graphRatioCombPHOSIncGammaSysErr0020 	= (TGraphAsymmErrors*) graphPHOSIncGammaSysErr0020->Clone();	
	TGraphAsymmErrors* graphRatioCombPCMIncGammaStatErr0020 	= (TGraphAsymmErrors*) graphPCMIncGammaStatErr0020->Clone();	
	TGraphAsymmErrors* graphRatioCombPCMIncGammaSysErr0020 		= (TGraphAsymmErrors*) graphPCMIncGammaSysErr0020->Clone();	
	graphRatioCombPHOSIncGammaStatErr0020 						= CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaStatErr0020, fitIncGammaCombQCD0020); 
	graphRatioCombPHOSIncGammaSysErr0020 						= CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaSysErr0020, fitIncGammaCombQCD0020); 
	graphRatioCombPCMIncGammaStatErr0020 						= CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaStatErr0020, fitIncGammaCombQCD0020); 
	graphRatioCombPCMIncGammaSysErr0020 						= CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaSysErr0020, fitIncGammaCombQCD0020); 
	
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatErr0020 	= (TGraphAsymmErrors*) graphPHOSIncGammaStatErr0020->Clone();		
	graphRatioCombPPHOSIncGammaStatErr0020 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatErr0020, graphCombIncGammaSumErr0020);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysErr0020 	= (TGraphAsymmErrors*) graphPHOSIncGammaSysErr0020->Clone();		
	graphRatioCombPPHOSIncGammaSysErr0020 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysErr0020, graphCombIncGammaSumErr0020);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysBErr0020 	= (TGraphAsymmErrors*) graphPHOSIncGammaSysBErr0020->Clone();		
	graphRatioCombPPHOSIncGammaSysBErr0020 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysBErr0020, graphCombIncGammaSumErr0020);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysCErr0020 	= (TGraphAsymmErrors*) graphPHOSIncGammaSysCErr0020->Clone();		
	graphRatioCombPPHOSIncGammaSysCErr0020 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysCErr0020, graphCombIncGammaSumErr0020);

	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr0020 	= (TGraphAsymmErrors*) graphPHOSIncGammaStatSysAErr0020->Clone();		
	graphRatioCombPPHOSIncGammaStatSysAErr0020 					= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatSysAErr0020, graphCombIncGammaSumErr0020);

	
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatErr0020 	= (TGraphAsymmErrors*) graphPCMIncGammaStatErr0020->Clone();	
	graphRatioCombPPCMIncGammaStatErr0020 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatErr0020, graphCombIncGammaSumErr0020);
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysErr0020 	= (TGraphAsymmErrors*) graphPCMIncGammaSysErr0020->Clone();	
	graphRatioCombPPCMIncGammaSysErr0020 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysErr0020, graphCombIncGammaSumErr0020);	
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysBErr0020 	= (TGraphAsymmErrors*) graphPCMIncGammaSysBErr0020->Clone();		
	graphRatioCombPPCMIncGammaSysBErr0020 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysBErr0020, graphCombIncGammaSumErr0020);
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysCErr0020 	= (TGraphAsymmErrors*) graphPCMIncGammaSysCErr0020->Clone();		
	graphRatioCombPPCMIncGammaSysCErr0020 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysCErr0020, graphCombIncGammaSumErr0020);
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr0020 	= (TGraphAsymmErrors*) graphPCMIncGammaStatSysAErr0020->Clone();		
	graphRatioCombPPCMIncGammaStatSysAErr0020 					= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatSysAErr0020, graphCombIncGammaSumErr0020);

	Double_t SysCPHOSIncGamma0020								= graphPHOSIncGammaSysCErr0020->GetErrorYlow(4)/graphPHOSIncGammaSysCErr0020->GetY()[4];
	Double_t SysCPCMIncGamma0020								= graphPCMIncGammaSysCErr0020->GetErrorYlow(4)/graphPCMIncGammaSysCErr0020->GetY()[4];

	TGraphAsymmErrors* graphRatioCombCombFitIncGammaStatErr0020	= (TGraphAsymmErrors*)graphCombIncGammaStatErr0020->Clone();
	TGraphAsymmErrors* graphRatioCombCombFitIncGammaSysErr0020 	= (TGraphAsymmErrors*)graphCombIncGammaSysErr0020->Clone();
	graphRatioCombCombFitIncGammaStatErr0020 					= CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaStatErr0020, fitIncGammaCombQCD0020); 
	graphRatioCombCombFitIncGammaSysErr0020 					= CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaSysErr0020, fitIncGammaCombQCD0020); 	

	//__________________________________________________ 20-40% combine spectra _________________________________________________________________
	TGraphAsymmErrors *graphCombIncGammaSysErr2040;
	TGraphAsymmErrors *graphCombIncGammaSysAErr2040;
	TGraphAsymmErrors *graphCombIncGammaSysBErr2040;
	TGraphAsymmErrors *graphCombIncGammaSysCErr2040;
	TGraphAsymmErrors *graphCombIncGammaStatErr2040;
	TGraphAsymmErrors *graphCombIncGammaSumErr2040;		
	graphCombIncGammaSumErr2040 	= CombinePtPointsSpectraAdv(	histoPCMIncGammaStatErr2040, graphPCMIncGammaSysErr2040,
																graphPCMIncGammaSysAErr2040, graphPCMIncGammaSysBErr2040, graphPCMIncGammaSysCErr2040,
																histoPHOSIncGammaStatErr2040, graphPHOSIncGammaSysErr2040 ,
																graphPHOSIncGammaSysAErr2040, graphPHOSIncGammaSysBErr2040, graphPHOSIncGammaSysCErr2040,
																graphCombIncGammaStatErr2040, graphCombIncGammaSysErr2040,
																graphCombIncGammaSysAErr2040, graphCombIncGammaSysBErr2040, graphCombIncGammaSysCErr2040,
																newBinsComb, 21, 0, 0, -1);

	//__________________________________________ 20-40% fit combined and build ratio of individual to fit _______________________________________
	TF1* fitIncGammaCombQCD2040	 								= FitObject("qcd","fitIncGammaCombQCD2040","Photon",graphCombIncGammaSumErr2040,0.9,14,NULL,"QNRMEX0+");
	cout << WriteParameterToFile(fitIncGammaCombQCD2040)<< endl;   
	TH1D* histoFitQCDIncGammaComb2040 							= (TH1D*)fitIncGammaCombQCD2040->GetHistogram();
	
	TGraphAsymmErrors* graphPCMIncGammaStatErr2040 				= new TGraphAsymmErrors(histoPCMIncGammaStatErr2040);
	graphPCMIncGammaStatErr2040->RemovePoint(0);
	TGraphAsymmErrors* graphPCMIncGammaStatSysAErr2040 			= AddErrorsOfGraphsQuadratically (graphPCMIncGammaStatErr2040, graphPCMIncGammaSysAErr2040);
	
	TGraphAsymmErrors* graphPHOSIncGammaStatErr2040 			= new TGraphAsymmErrors(histoPHOSIncGammaStatErr2040);
	graphPHOSIncGammaStatErr2040->RemovePoint(0);
	graphPHOSIncGammaSysErr2040->RemovePoint(0);
	graphPHOSIncGammaSysAErr2040->RemovePoint(0);
	graphPHOSIncGammaSysBErr2040->RemovePoint(0);
	graphPHOSIncGammaSysCErr2040->RemovePoint(0);
	TGraphAsymmErrors* graphPHOSIncGammaStatSysAErr2040 		= AddErrorsOfGraphsQuadratically (graphPHOSIncGammaStatErr2040, graphPHOSIncGammaSysAErr2040);
	
	
	TGraphAsymmErrors* graphRatioCombPHOSIncGammaStatErr2040 	= (TGraphAsymmErrors*) graphPHOSIncGammaStatErr2040->Clone();	
	TGraphAsymmErrors* graphRatioCombPHOSIncGammaSysErr2040 	= (TGraphAsymmErrors*) graphPHOSIncGammaSysErr2040->Clone();	
	TGraphAsymmErrors* graphRatioCombPCMIncGammaStatErr2040 	= (TGraphAsymmErrors*) graphPCMIncGammaStatErr2040->Clone();	
	TGraphAsymmErrors* graphRatioCombPCMIncGammaSysErr2040 		= (TGraphAsymmErrors*) graphPCMIncGammaSysErr2040->Clone();	
	graphRatioCombPHOSIncGammaStatErr2040 						= CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaStatErr2040, fitIncGammaCombQCD2040); 
	graphRatioCombPHOSIncGammaSysErr2040 						= CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaSysErr2040, fitIncGammaCombQCD2040); 
	graphRatioCombPCMIncGammaStatErr2040 						= CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaStatErr2040, fitIncGammaCombQCD2040); 
	graphRatioCombPCMIncGammaSysErr2040 						= CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaSysErr2040, fitIncGammaCombQCD2040); 
	TGraphAsymmErrors* graphRatioCombCombFitIncGammaStatErr2040	= (TGraphAsymmErrors*)graphCombIncGammaStatErr2040->Clone();
	TGraphAsymmErrors* graphRatioCombCombFitIncGammaSysErr2040 	= (TGraphAsymmErrors*)graphCombIncGammaSysErr2040->Clone();
	graphRatioCombCombFitIncGammaStatErr2040 					= CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaStatErr2040, fitIncGammaCombQCD2040); 
	graphRatioCombCombFitIncGammaSysErr2040 					= CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaSysErr2040, fitIncGammaCombQCD2040); 

	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatErr2040 	= (TGraphAsymmErrors*) graphPHOSIncGammaStatErr2040->Clone();		
	graphRatioCombPPHOSIncGammaStatErr2040 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatErr2040, graphCombIncGammaSumErr2040);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysErr2040 	= (TGraphAsymmErrors*) graphPHOSIncGammaSysErr2040->Clone();		
	graphRatioCombPPHOSIncGammaSysErr2040 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysErr2040, graphCombIncGammaSumErr2040);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysBErr2040 	= (TGraphAsymmErrors*) graphPHOSIncGammaSysBErr2040->Clone();		
	graphRatioCombPPHOSIncGammaSysBErr2040 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysBErr2040, graphCombIncGammaSumErr2040);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr2040 	= (TGraphAsymmErrors*) graphPHOSIncGammaStatSysAErr2040->Clone();		
	graphRatioCombPPHOSIncGammaStatSysAErr2040 					= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatSysAErr2040, graphCombIncGammaSumErr2040);

	TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatErr2040 	= (TGraphAsymmErrors*) graphPCMIncGammaStatErr2040->Clone();	
	graphRatioCombPPCMIncGammaStatErr2040 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatErr2040, graphCombIncGammaSumErr2040);
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysErr2040 	= (TGraphAsymmErrors*) graphPCMIncGammaSysErr2040->Clone();	
	graphRatioCombPPCMIncGammaSysErr2040 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysErr2040, graphCombIncGammaSumErr2040);	
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysBErr2040 	= (TGraphAsymmErrors*) graphPCMIncGammaSysBErr2040->Clone();	
	graphRatioCombPPCMIncGammaSysBErr2040 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysBErr2040, graphCombIncGammaSumErr2040);	
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr2040 	= (TGraphAsymmErrors*) graphPCMIncGammaStatSysAErr2040->Clone();		
	graphRatioCombPPCMIncGammaStatSysAErr2040 					= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatSysAErr2040, graphCombIncGammaSumErr2040);

	Double_t SysCPHOSIncGamma2040								= graphPHOSIncGammaSysCErr2040->GetErrorYlow(4)/graphPHOSIncGammaSysCErr2040->GetY()[4];
	Double_t SysCPCMIncGamma2040								= graphPCMIncGammaSysCErr2040->GetErrorYlow(4)/graphPCMIncGammaSysCErr2040->GetY()[4];
	
	
	//__________________________________________________ 40-80% combine spectra _________________________________________________________________
	TGraphAsymmErrors *graphCombIncGammaSysErr4080;
	TGraphAsymmErrors *graphCombIncGammaSysAErr4080;
	TGraphAsymmErrors *graphCombIncGammaSysBErr4080;
	TGraphAsymmErrors *graphCombIncGammaSysCErr4080;
	TGraphAsymmErrors *graphCombIncGammaStatErr4080;
	TGraphAsymmErrors *graphCombIncGammaSumErr4080;		
	graphCombIncGammaSumErr4080 	= CombinePtPointsSpectraAdv(	histoPCMIncGammaStatErr4080, graphPCMIncGammaSysErr4080,
																graphPCMIncGammaSysAErr4080, graphPCMIncGammaSysBErr4080, graphPCMIncGammaSysCErr4080,
																histoPHOSIncGammaStatErr4080, graphPHOSIncGammaSysErr4080 ,
																graphPHOSIncGammaSysAErr4080, graphPHOSIncGammaSysBErr4080, graphPHOSIncGammaSysCErr4080,
																graphCombIncGammaStatErr4080, graphCombIncGammaSysErr4080,
																graphCombIncGammaSysAErr4080, graphCombIncGammaSysBErr4080, graphCombIncGammaSysCErr4080,
																newBinsComb, 21, 0, 0, -1);

	//__________________________________________ 40-80% fit combined and build ratio of individual to fit _______________________________________
	TF1* fitIncGammaCombQCD4080	 								= FitObject("qcd","fitIncGammaCombQCD4080","Photon",graphCombIncGammaSumErr4080,0.9,14,NULL,"QNRMEX0+");
	cout << WriteParameterToFile(fitIncGammaCombQCD4080)<< endl;   
	TH1D* histoFitQCDIncGammaComb4080 							= (TH1D*)fitIncGammaCombQCD4080->GetHistogram();
	
	TGraphAsymmErrors* graphPCMIncGammaStatErr4080 				= new TGraphAsymmErrors(histoPCMIncGammaStatErr4080);
	graphPCMIncGammaStatErr4080->RemovePoint(0);
	TGraphAsymmErrors* graphPCMIncGammaStatSysAErr4080 			= AddErrorsOfGraphsQuadratically (graphPCMIncGammaStatErr4080, graphPCMIncGammaSysAErr4080);
	TGraphAsymmErrors* graphPHOSIncGammaStatErr4080 			= new TGraphAsymmErrors(histoPHOSIncGammaStatErr4080);
	graphPHOSIncGammaStatErr4080->RemovePoint(0);
	graphPHOSIncGammaSysErr4080->RemovePoint(0);
	graphPHOSIncGammaSysAErr4080->RemovePoint(0);
	graphPHOSIncGammaSysBErr4080->RemovePoint(0);
	graphPHOSIncGammaSysCErr4080->RemovePoint(0);
	TGraphAsymmErrors* graphPHOSIncGammaStatSysAErr4080 		= AddErrorsOfGraphsQuadratically (graphPHOSIncGammaStatErr4080, graphPHOSIncGammaSysAErr4080);
	
	TGraphAsymmErrors* graphRatioCombPHOSIncGammaStatErr4080 	= (TGraphAsymmErrors*) graphPHOSIncGammaStatErr4080->Clone();	
	TGraphAsymmErrors* graphRatioCombPHOSIncGammaSysErr4080 	= (TGraphAsymmErrors*) graphPHOSIncGammaSysErr4080->Clone();	
	TGraphAsymmErrors* graphRatioCombPCMIncGammaStatErr4080 	= (TGraphAsymmErrors*) graphPCMIncGammaStatErr4080->Clone();	
	TGraphAsymmErrors* graphRatioCombPCMIncGammaSysErr4080 		= (TGraphAsymmErrors*) graphPCMIncGammaSysErr4080->Clone();	
	graphRatioCombPHOSIncGammaStatErr4080 						= CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaStatErr4080, fitIncGammaCombQCD4080); 
	graphRatioCombPHOSIncGammaSysErr4080 						= CalculateGraphErrRatioToFit (graphRatioCombPHOSIncGammaSysErr4080, fitIncGammaCombQCD4080); 
	graphRatioCombPCMIncGammaStatErr4080 						= CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaStatErr4080, fitIncGammaCombQCD4080); 
	graphRatioCombPCMIncGammaSysErr4080 						= CalculateGraphErrRatioToFit (graphRatioCombPCMIncGammaSysErr4080, fitIncGammaCombQCD4080); 
	TGraphAsymmErrors* graphRatioCombCombFitIncGammaStatErr4080	= (TGraphAsymmErrors*)graphCombIncGammaStatErr4080->Clone();
	TGraphAsymmErrors* graphRatioCombCombFitIncGammaSysErr4080 	= (TGraphAsymmErrors*)graphCombIncGammaSysErr4080->Clone();
	graphRatioCombCombFitIncGammaStatErr4080 					= CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaStatErr4080, fitIncGammaCombQCD4080); 
	graphRatioCombCombFitIncGammaSysErr4080 					= CalculateGraphErrRatioToFit(graphRatioCombCombFitIncGammaSysErr4080, fitIncGammaCombQCD4080); 

	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatErr4080 	= (TGraphAsymmErrors*) graphPHOSIncGammaStatErr4080->Clone();		
	graphRatioCombPPHOSIncGammaStatErr4080 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatErr4080, graphCombIncGammaSumErr4080);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysErr4080 	= (TGraphAsymmErrors*) graphPHOSIncGammaSysErr4080->Clone();		
	graphRatioCombPPHOSIncGammaSysErr4080 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysErr4080, graphCombIncGammaSumErr4080);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaSysBErr4080 	= (TGraphAsymmErrors*) graphPHOSIncGammaSysBErr4080->Clone();		
	graphRatioCombPPHOSIncGammaSysBErr4080 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaSysBErr4080, graphCombIncGammaSumErr4080);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr4080 	= (TGraphAsymmErrors*) graphPHOSIncGammaStatSysAErr4080->Clone();		
	graphRatioCombPPHOSIncGammaStatSysAErr4080 					= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPHOSIncGammaStatSysAErr4080, graphCombIncGammaSumErr4080);

	TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatErr4080 	= (TGraphAsymmErrors*) graphPCMIncGammaStatErr4080->Clone();	
	graphRatioCombPPCMIncGammaStatErr4080 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatErr4080, graphCombIncGammaSumErr4080);
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysErr4080 	= (TGraphAsymmErrors*) graphPCMIncGammaSysErr4080->Clone();	
	graphRatioCombPPCMIncGammaSysErr4080 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysErr4080, graphCombIncGammaSumErr4080);	
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaSysBErr4080 	= (TGraphAsymmErrors*) graphPCMIncGammaSysBErr4080->Clone();	
	graphRatioCombPPCMIncGammaSysBErr4080 						= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaSysBErr4080, graphCombIncGammaSumErr4080);	
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr4080 	= (TGraphAsymmErrors*) graphPCMIncGammaStatSysAErr4080->Clone();		
	graphRatioCombPPCMIncGammaStatSysAErr4080 					= CalculateGraphErrRatioToOtherTGraphErr(graphRatioCombPPCMIncGammaStatSysAErr4080, graphCombIncGammaSumErr4080);

	Double_t SysCPHOSIncGamma4080								= graphPHOSIncGammaSysCErr4080->GetErrorYlow(4)/graphPHOSIncGammaSysCErr4080->GetY()[4];
	Double_t SysCPCMIncGamma4080								= graphPCMIncGammaSysCErr4080->GetErrorYlow(4)/graphPCMIncGammaSysCErr4080->GetY()[4];
	
	//*******************************************************************************************************************************************
	//**************************************************** Calculate direct photon spectrum 0020 ************************************************
	//*******************************************************************************************************************************************
	cout << endl;
	graphCombIncGammaStatErr0020->Print();
	Double_t xArrayCombined[graphCombDRPi0FitStatErr0020->GetN()+1];
	xArrayCombined[0] = graphCombDRPi0FitStatErr0020->GetX()[0] - graphCombDRPi0FitStatErr0020->GetEXhigh()[0];
	cout << "Binning \n" << xArrayCombined[0] << endl;
	for (Int_t i = 1; i<graphCombDRPi0FitStatErr0020->GetN()+1;i++){
		xArrayCombined[i] = graphCombDRPi0FitStatErr0020->GetX()[i-1] + graphCombDRPi0FitStatErr0020->GetEXhigh()[i-1];
		cout << xArrayCombined[i] << endl;
	}	
	
	//_______________________ copy inclusive photon spectra _____________________
	TH1D *histoCombDirGammaSpectrumErrSum0020 					= new TH1D("histoCombDirGammaSpectrumErrSum0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSys0020 					= new TH1D("histoCombDirGammaSpectrumErrSys0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSysA0020 					= new TH1D("histoCombDirGammaSpectrumErrSysA0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSysB0020 					= new TH1D("histoCombDirGammaSpectrumErrSysB0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSysC0020 					= new TH1D("histoCombDirGammaSpectrumErrSysC0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrStat0020 					= new TH1D("histoCombDirGammaSpectrumErrStat0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);

	//_______________________ get arrays of double ratio errors __________________
	Double_t *SystErrorsCombDR0020 								= new Double_t[graphCombIncGammaStatErr0020->GetN()];
	Double_t *SystAErrorsCombDR0020 							= new Double_t[graphCombIncGammaStatErr0020->GetN()];
	Double_t *SystBErrorsCombDR0020 							= new Double_t[graphCombIncGammaStatErr0020->GetN()];
	Double_t *SystCErrorsCombDR0020 							= new Double_t[graphCombIncGammaStatErr0020->GetN()];
	Double_t *sumErrorsCombDR0020 								= new Double_t[graphCombIncGammaStatErr0020->GetN()];
	Double_t *StatErrorsCombDR0020 								= new Double_t[graphCombIncGammaStatErr0020->GetN()];
	Double_t *xErrorsDR0020 									= new Double_t[graphCombIncGammaStatErr0020->GetN()];
	graphCombIncGammaSysErr0020->Print();
	for (Int_t i = 0; i< graphCombDRPi0FitStatErr0020->GetN(); i++){
		SystErrorsCombDR0020[i] 								= graphCombDRPi0FitSysErr0020->GetEYhigh()[i]/graphCombDRPi0FitSysErr0020->GetY()[i] *100;
		SystAErrorsCombDR0020[i] 								= graphCombDRPi0FitSysAErr0020->GetEYhigh()[i]/graphCombDRPi0FitSysAErr0020->GetY()[i] *100;
		SystBErrorsCombDR0020[i] 								= graphCombDRPi0FitSysBErr0020->GetEYhigh()[i]/graphCombDRPi0FitSysBErr0020->GetY()[i] *100;
		SystCErrorsCombDR0020[i] 								= graphCombDRPi0FitSysCErr0020->GetEYhigh()[i]/graphCombDRPi0FitSysCErr0020->GetY()[i] *100;
		StatErrorsCombDR0020[i] 								= graphCombDRPi0FitStatErr0020->GetEYhigh()[i]/graphCombDRPi0FitStatErr0020->GetY()[i] *100;
		sumErrorsCombDR0020[i] 									= graphCombDRPi0FitSumErr0020->GetEYhigh()[i]/graphCombDRPi0FitSumErr0020->GetY()[i] *100;
// 		cout << i << "\t" << graphCombDRPi0FitSysErr0020->GetY()[i] << "\t" << graphCombDRPi0FitSysErr0020->GetEYhigh()[i] << "\t" <<SystErrorsCombDR0020[i] << endl;
	}
	xErrorsDR0020 												= graphCombDRPi0FitStatErr0020->GetX();

	cout << "here !!! \n\n" << endl;
	graphCombDRPi0FitSumErr0020->Print();
	
	
	//_______________________ copy inclusive photon spectra _____________________	
	TH1D* histoCombErrorsForDRSum0020 							= new TH1D("histoCombErrorsForDRSum0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRStat0020 							= new TH1D("histoCombErrorsForDRStat0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSys0020 							= new TH1D("histoCombErrorsForDRSys0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSysA0020 							= new TH1D("histoCombErrorsForDRSysA0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSysB0020 							= new TH1D("histoCombErrorsForDRSysB0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSysC0020 							= new TH1D("histoCombErrorsForDRSysC0020","",graphCombDRPi0FitStatErr0020->GetN(),xArrayCombined);
	
	for(Int_t i = 1; i<graphCombDRPi0FitStatErr0020->GetN()+1;i++){
		cout<< i << "\t"<<xErrorsDR0020[i-1]<<"  "<<histoCombErrorsForDRSum0020->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum0020->GetBinWidth(i) <<endl;
		Double_t binErrorSummed 								= sumErrorsCombDR0020[i-1];
		Double_t binErrorSyst 									= SystErrorsCombDR0020[i-1];
		Double_t binErrorSystA 									= SystAErrorsCombDR0020[i-1];
		Double_t binErrorSystB 									= SystBErrorsCombDR0020[i-1];
		Double_t binErrorSystC 									= SystCErrorsCombDR0020[i-1];
		Double_t binErrorStat 									= StatErrorsCombDR0020[i-1];
		Double_t DR 											= graphCombDRPi0FitStatErr0020->GetY()[i-1];

		cout << DR << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
		histoCombErrorsForDRSum0020->SetBinContent(i,DR);
		histoCombErrorsForDRSys0020->SetBinContent(i,DR);
		histoCombErrorsForDRSysA0020->SetBinContent(i,DR);
		histoCombErrorsForDRSysB0020->SetBinContent(i,DR);
		histoCombErrorsForDRSysC0020->SetBinContent(i,DR);
		histoCombErrorsForDRStat0020->SetBinContent(i,DR);
		histoCombErrorsForDRSum0020->SetBinError(i,(binErrorSummed/100)*DR);
		histoCombErrorsForDRSys0020->SetBinError(i,(binErrorSyst/100)*DR);
		histoCombErrorsForDRSysA0020->SetBinError(i,(binErrorSystA/100)*DR);
		histoCombErrorsForDRSysB0020->SetBinError(i,(binErrorSystB/100)*DR);
		histoCombErrorsForDRSysC0020->SetBinError(i,(binErrorSystC/100)*DR);
		histoCombErrorsForDRStat0020->SetBinError(i,(binErrorStat/100)*DR);
	}
   
	for(Int_t i = 1; i<histoCombErrorsForDRSum0020->GetNbinsX()+1;i++){
		histoCombDirGammaSpectrumErrSum0020->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSys0020->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSysA0020->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSysB0020->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSysC0020->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrStat0020->SetBinContent(i+1,-1);

		histoCombDirGammaSpectrumErrSum0020->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSys0020->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSysA0020->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSysB0020->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSysC0020->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrStat0020->SetBinError(i+1,0);
	}

	// get the binning of the direct photons from the DR
	TH1D *hisoCombDirGammaSpecSysErr0020 				= new TH1D(*histoCombErrorsForDRSys0020);
	TH1D *hisoCombDirGammaSpecSysAErr0020 				= new TH1D(*histoCombErrorsForDRSysA0020);
	TH1D *hisoCombDirGammaSpecSysBErr0020 				= new TH1D(*histoCombErrorsForDRSysB0020);
	TH1D *hisoCombDirGammaSpecSysCErr0020 				= new TH1D(*histoCombErrorsForDRSysC0020);
	TH1D *hisoCombDirGammaSpecStatErr0020 				= new TH1D(*histoCombErrorsForDRStat0020);
	TH1D *hisoCombDirGammaSpecSumErr0020 				= new TH1D(*histoCombErrorsForDRSum0020);

	for(Int_t i = 1; i<graphCombDRPi0FitStatErr0020->GetN()+1; i++){
		// obtain common quantities
		Double_t Rgamma 			= histoCombErrorsForDRSys0020->GetBinContent(i);
		Double_t nIncGamma			= graphCombIncGammaStatErr0020->GetY()[i-1];
		
		// calculating Systematics graph
		Double_t errRgamma			= histoCombErrorsForDRSys0020->GetBinError(i);
		Double_t errNIncGam 		= graphCombIncGammaSysErr0020->GetEYhigh()[i-1];
		Double_t q1 				= 1 - 1/ Rgamma;
		
		Double_t q1Error 			= errRgamma/(Rgamma*Rgamma);
		Double_t content 			= nIncGamma * ( 1 - 1/ Rgamma);
		Double_t error 				= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		Double_t errDR				= content - error;
		hisoCombDirGammaSpecSysErr0020->SetBinError(i, error);
		hisoCombDirGammaSpecSysErr0020->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSys0020->SetBinContent(i, errDR);

		// calculating Systematics A graph
		errRgamma					= histoCombErrorsForDRSysA0020->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSysAErr0020->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSysAErr0020->SetBinError(i, error);
		hisoCombDirGammaSpecSysAErr0020->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSysA0020->SetBinContent(i, errDR);

		// calculating Systematics B graph
		errRgamma					= histoCombErrorsForDRSysB0020->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSysBErr0020->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSysBErr0020->SetBinError(i, error);
		hisoCombDirGammaSpecSysBErr0020->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSysB0020->SetBinContent(i, errDR);

		// calculating Systematics C graph
		errRgamma					= histoCombErrorsForDRSysC0020->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSysCErr0020->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSysCErr0020->SetBinError(i, error);
		hisoCombDirGammaSpecSysCErr0020->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSysC0020->SetBinContent(i, errDR);		
		
		// calculating Stat graphs
		errRgamma					= histoCombErrorsForDRStat0020->GetBinError(i);
		errNIncGam 					= graphCombIncGammaStatErr0020->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecStatErr0020->SetBinError(i, error);
		hisoCombDirGammaSpecStatErr0020->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrStat0020->SetBinContent(i, errDR);
		
		// calculating summed error graphs
		errRgamma					= hisoCombDirGammaSpecSumErr0020->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSumErr0020->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSumErr0020->SetBinError(i, error);
		hisoCombDirGammaSpecSumErr0020->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSum0020->SetBinContent(i, errDR);
	}

	// purely calculating points based on all Systematic errors
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSys0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
	if(graphCombDirGammaSpectrumSystErr0020)graphCombDirGammaSpectrumSystErr0020->SetName("graphCombDirGammaSpectrumSystErr0020");
	if(graphCombDirGammaSpectrumSystErr0020)graphCombDirGammaSpectrumSystErr0020->Print();
	if(graphCombDirGammaSpectrumSystErr0020)cout << "graph has been found" << endl;
	// purely calculating points based on all Systematic errors A
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystAErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysA0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
	if(graphCombDirGammaSpectrumSystAErr0020)graphCombDirGammaSpectrumSystAErr0020->SetName("graphCombDirGammaSpectrumSystAErr0020");
	if(graphCombDirGammaSpectrumSystAErr0020)graphCombDirGammaSpectrumSystAErr0020->Print();
	// purely calculating points based on all Systematic errors B
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystBErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysB0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
	if(graphCombDirGammaSpectrumSystBErr0020)graphCombDirGammaSpectrumSystBErr0020->SetName("graphCombDirGammaSpectrumSystBErr0020");
	if(graphCombDirGammaSpectrumSystBErr0020)graphCombDirGammaSpectrumSystBErr0020->Print();
	// purely calculating points based on all Systematic errors C
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystCErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysC0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
	if(graphCombDirGammaSpectrumSystCErr0020)graphCombDirGammaSpectrumSystCErr0020->SetName("graphCombDirGammaSpectrumSystCErr0020");
	if(graphCombDirGammaSpectrumSystCErr0020)graphCombDirGammaSpectrumSystCErr0020->Print();

	// purely calculating points based on Statistical errors
	TGraphAsymmErrors *graphCombDirGammaSpectrumStatErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrStat0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
	if(graphCombDirGammaSpectrumStatErr0020)graphCombDirGammaSpectrumStatErr0020->SetName("graphCombDirGammaSpectrumStatErr0020");
	// purely calculating points based on all Systematic + Statistical errors
	TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr0020 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum0020,hisoCombDirGammaSpecStatErr0020,0,0.5);
	if(graphCombDirGammaSpectrumSumErr0020)graphCombDirGammaSpectrumSumErr0020->SetName("graphCombDirGammaSpectrumSumErr0020");
 	// calculate points above confidence level summed errors
	TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr0020Confi = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum0020,hisoCombDirGammaSpecStatErr0020,2,0.5);
	if(graphCombDirGammaSpectrumSumErr0020Confi)graphCombDirGammaSpectrumSumErr0020Confi->SetName("graphCombDirGammaSpectrumSumErr0020Confi");
 	// calculate arrows for points with 0, error summed
	TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr0020Ar = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum0020,hisoCombDirGammaSpecStatErr0020,5,0.5);
	if(graphCombDirGammaSpectrumSumErr0020Ar)graphCombDirGammaSpectrumSumErr0020Ar->SetName("graphCombDirGammaSpectrumSumErr0020Ar");


	TF1* fitThermalGamma0020Sum						 			= FitObject("e","fitThermalGamma0020Sum","Photon",histoCombDirGammaSpectrumErrSum0020,0.9,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma0020Sum)<< endl;   	
	TF1* fitThermalGamma0020Sum23						 		= FitObject("e","fitThermalGamma0020Sum23","Photon",histoCombDirGammaSpectrumErrSum0020,0.9,2.3,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma0020Sum23)<< endl;   	
	TF1* fitThermalGamma0020Stat						 		= FitObject("e","fitThermalGamma0020Stat","Photon",graphCombDirGammaSpectrumStatErr0020,0.9,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma0020Stat)<< endl;   
	
	TF1* fitThermalGamma0020Sys						 			= FitObject("e","fitThermalGamma0020Sys","Photon",graphCombDirGammaSpectrumSystErr0020,0.9,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma0020Sys)<< endl;   
	TF1* fitThermalGamma0020SysA						 		= FitObject("e","fitThermalGamma0020SysA","Photon",graphCombDirGammaSpectrumSystAErr0020,0.9,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma0020SysA)<< endl;   
	TF1* fitThermalGamma0020SysB						 		= FitObject("e","fitThermalGamma0020SysB","Photon",graphCombDirGammaSpectrumSystBErr0020,0.9,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma0020SysB)<< endl;   
	TF1* fitThermalGamma0020SysC						 		= FitObject("e","fitThermalGamma0020SysC","Photon",graphCombDirGammaSpectrumSystCErr0020,0.9,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma0020SysC)<< endl;   
	
	TGraphAsymmErrors* graphCombDirGammaSpectrumErrSum0020 			= new TGraphAsymmErrors(histoCombDirGammaSpectrumErrSum0020);
	TF1* fitFullDirGamma0020Sys									= FitObject("qcd","fitFullDirGamma0020Sys","Photon",graphCombDirGammaSpectrumSystErr0020,0.9,14,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitFullDirGamma0020Sys)<< endl;   	
	TF1* fitFullDirGamma0020Stat								= FitObject("qcd","fitFullDirGamma0020Stat","Photon",graphCombDirGammaSpectrumStatErr0020,0.9,14,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitFullDirGamma0020Stat)<< endl;   	
	
	TGraphAsymmErrors* graphRatioCombFitDirGammaStatErr0020	= (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr0020->Clone();
	TGraphAsymmErrors* graphRatioCombFitDirGammaSysErr0020 	= (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr0020->Clone();
	graphRatioCombFitDirGammaStatErr0020 					= CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaStatErr0020, fitFullDirGamma0020Stat); 
	graphRatioCombFitDirGammaSysErr0020 					= CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaSysErr0020, fitFullDirGamma0020Sys); 

	// Calculate thermal spectrum
	TGraphAsymmErrors* graphCombThermalGammaSpectrumStatErr0020 	= NULL;
	TF1* fitPureThermalGamma0020Stat								= NULL;
	TGraphAsymmErrors* graphCombThermalGammaSpectrumSysErr0020 		= NULL;
	TGraphAsymmErrors* graphCombThermalGammaSpectrumSumErr0020Ar	= NULL;
	TF1* fitFullThermalGamma0020Sys									= NULL;
	TF1* fitFullThermalGamma0020Stat								= NULL;
	TGraphAsymmErrors* graphRatioCombFitThermalGammaStatErr0020		= NULL;
	TGraphAsymmErrors* graphRatioCombFitThermalGammaSysErr0020 		= NULL;
	if (graphCombDirGammaSpectrumStatErr0020){
		graphCombThermalGammaSpectrumStatErr0020 				= SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumStatErr0020);
		graphCombThermalGammaSpectrumStatErr0020->Print();
		fitPureThermalGamma0020Stat						 		= FitObject("e","fitPureThermalGamma0020Stat","Photon",graphCombThermalGammaSpectrumStatErr0020,0.9,2.1,NULL,"QNRMEX0+");
		fileFinalResults << WriteParameterToFile(fitPureThermalGamma0020Stat)<< endl;   
		fitFullThermalGamma0020Stat								= FitObject("qcd","fitFullDirGamma0020Stat","Photon",graphCombThermalGammaSpectrumStatErr0020,0.9,14,NULL,"QNRMEX0+");
		fileFinalResults << WriteParameterToFile(fitFullThermalGamma0020Stat)<< endl;   	
		graphRatioCombFitThermalGammaStatErr0020				= (TGraphAsymmErrors*)graphCombThermalGammaSpectrumStatErr0020->Clone();
		graphRatioCombFitThermalGammaStatErr0020 				= CalculateGraphErrRatioToFit(graphRatioCombFitThermalGammaStatErr0020, fitFullThermalGamma0020Stat); 
	}
	if (graphCombDirGammaSpectrumSystErr0020){
		graphCombThermalGammaSpectrumSysErr0020 				= SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSystErr0020);
		graphCombThermalGammaSpectrumSysErr0020->Print();
		fitFullThermalGamma0020Sys								= FitObject("qcd","fitFullDirGamma0020Sys","Photon",graphCombThermalGammaSpectrumSysErr0020,0.9,14,NULL,"QNRMEX0+");
		fileFinalResults << WriteParameterToFile(fitFullThermalGamma0020Sys)<< endl;   	
		graphRatioCombFitThermalGammaSysErr0020 				= (TGraphAsymmErrors*)graphCombThermalGammaSpectrumSysErr0020->Clone();
		graphRatioCombFitThermalGammaSysErr0020 				= CalculateGraphErrRatioToFit(graphRatioCombFitThermalGammaSysErr0020, fitFullThermalGamma0020Sys); 
		
	}
	
	if (graphCombDirGammaSpectrumSumErr0020Ar){
		graphCombThermalGammaSpectrumSumErr0020Ar	= SubtractPromptPhotonsViaFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSumErr0020Ar, newBinsComb, 20);
	}
	
	// Calculate RAA
	cout << endl << "Calculating RAA" << endl;
	TGraphAsymmErrors* graphCombRAADirGammaStat0020 	= NULL;
	TGraphAsymmErrors* graphCombRAADirGammaSys0020		= NULL;
	TGraphAsymmErrors* graphCombRAADirGammaSum0020Ar	= NULL;
	if (graphCombDirGammaSpectrumStatErr0020) 	CalcRaaWithTheoryFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumStatErr0020, &graphCombRAADirGammaStat0020);
	if (graphCombDirGammaSpectrumSystErr0020) 	CalcRaaWithTheoryFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSystErr0020, &graphCombRAADirGammaSys0020);
	if (graphCombDirGammaSpectrumSumErr0020Ar) 	CalcRaaWithTheoryFit( fitTheoryPromptMcGill0020, graphCombDirGammaSpectrumSumErr0020Ar, &graphCombRAADirGammaSum0020Ar, newBinsComb, 20);
	if (graphCombRAADirGammaStat0020) graphCombRAADirGammaStat0020->Print();
	if (graphCombRAADirGammaSys0020) graphCombRAADirGammaSys0020->Print();
	if (graphCombRAADirGammaSum0020Ar) graphCombRAADirGammaSum0020Ar->Print();

	//*******************************************************************************************************************************************
	//*********************************************** Calculate direct photon spectrum 2040 *****************************************************
	//*******************************************************************************************************************************************	
	//_______________________ copy inclusive photon spectra _____________________
	TH1D *histoCombDirGammaSpectrumErrSum2040 					= new TH1D("histoCombDirGammaSpectrumErrSum2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSys2040 					= new TH1D("histoCombDirGammaSpectrumErrSys2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSysA2040 					= new TH1D("histoCombDirGammaSpectrumErrSysA2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSysB2040 					= new TH1D("histoCombDirGammaSpectrumErrSysB2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSysC2040 					= new TH1D("histoCombDirGammaSpectrumErrSysC2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrStat2040 					= new TH1D("histoCombDirGammaSpectrumErrStat2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);

	//_______________________ get arrays of double ratio errors __________________
	Double_t *SystErrorsCombDR2040 								= new Double_t[graphCombIncGammaStatErr2040->GetN()];
	Double_t *SystAErrorsCombDR2040 								= new Double_t[graphCombIncGammaStatErr2040->GetN()];
	Double_t *SystBErrorsCombDR2040 								= new Double_t[graphCombIncGammaStatErr2040->GetN()];
	Double_t *SystCErrorsCombDR2040 								= new Double_t[graphCombIncGammaStatErr2040->GetN()];
	Double_t *sumErrorsCombDR2040 								= new Double_t[graphCombIncGammaStatErr2040->GetN()];
	Double_t *StatErrorsCombDR2040 								= new Double_t[graphCombIncGammaStatErr2040->GetN()];
	Double_t *xErrorsDR2040 									= new Double_t[graphCombIncGammaStatErr2040->GetN()];
	graphCombIncGammaSysErr2040->Print();
	for (Int_t i = 0; i< graphCombDRPi0FitStatErr2040->GetN(); i++){
		SystErrorsCombDR2040[i] 								= graphCombDRPi0FitSysErr2040->GetEYhigh()[i]/graphCombDRPi0FitSysErr2040->GetY()[i] *100;
		SystAErrorsCombDR2040[i] 								= graphCombDRPi0FitSysAErr2040->GetEYhigh()[i]/graphCombDRPi0FitSysAErr2040->GetY()[i] *100;
		SystBErrorsCombDR2040[i] 								= graphCombDRPi0FitSysBErr2040->GetEYhigh()[i]/graphCombDRPi0FitSysBErr2040->GetY()[i] *100;
		SystCErrorsCombDR2040[i] 								= graphCombDRPi0FitSysCErr2040->GetEYhigh()[i]/graphCombDRPi0FitSysCErr2040->GetY()[i] *100;
		StatErrorsCombDR2040[i] 								= graphCombDRPi0FitStatErr2040->GetEYhigh()[i]/graphCombDRPi0FitStatErr2040->GetY()[i] *100;
		sumErrorsCombDR2040[i] 									= graphCombDRPi0FitSumErr2040->GetEYhigh()[i]/graphCombDRPi0FitSumErr2040->GetY()[i] *100;
// 		cout << i << "\t" << graphCombDRPi0FitSysErr2040->GetY()[i] << "\t" << graphCombDRPi0FitSysErr2040->GetEYhigh()[i] << "\t" <<SystErrorsCombDR2040[i] << endl;
	}
	xErrorsDR2040 												= graphCombDRPi0FitStatErr2040->GetX();

	cout << "here !!! \n\n" << endl;
	graphCombDRPi0FitSumErr2040->Print();
	
	
	//_______________________ copy inclusive photon spectra _____________________	
	TH1D* histoCombErrorsForDRSum2040 							= new TH1D("histoCombErrorsForDRSum2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRStat2040 							= new TH1D("histoCombErrorsForDRStat2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSys2040 							= new TH1D("histoCombErrorsForDRSys2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSysA2040 							= new TH1D("histoCombErrorsForDRSysA2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSysB2040 							= new TH1D("histoCombErrorsForDRSysB2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSysC2040 							= new TH1D("histoCombErrorsForDRSysC2040","",graphCombDRPi0FitStatErr2040->GetN(),xArrayCombined);
	
	for(Int_t i = 1; i<graphCombDRPi0FitStatErr2040->GetN()+1;i++){
		cout<< i << "\t"<<xErrorsDR2040[i-1]<<"  "<<histoCombErrorsForDRSum2040->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum2040->GetBinWidth(i) <<endl;
		Double_t binErrorSummed 								= sumErrorsCombDR2040[i-1];
		Double_t binErrorSyst 									= SystErrorsCombDR2040[i-1];
		Double_t binErrorSystA 									= SystAErrorsCombDR2040[i-1];
		Double_t binErrorSystB 									= SystBErrorsCombDR2040[i-1];
		Double_t binErrorSystC 									= SystCErrorsCombDR2040[i-1];
		Double_t binErrorStat 									= StatErrorsCombDR2040[i-1];
		Double_t DR 											= graphCombDRPi0FitStatErr2040->GetY()[i-1];

		cout << DR << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
		histoCombErrorsForDRSum2040->SetBinContent(i,DR);
		histoCombErrorsForDRSys2040->SetBinContent(i,DR);
		histoCombErrorsForDRSysA2040->SetBinContent(i,DR);
		histoCombErrorsForDRSysB2040->SetBinContent(i,DR);
		histoCombErrorsForDRSysC2040->SetBinContent(i,DR);
		histoCombErrorsForDRStat2040->SetBinContent(i,DR);
		histoCombErrorsForDRSum2040->SetBinError(i,(binErrorSummed/100)*DR);
		histoCombErrorsForDRSys2040->SetBinError(i,(binErrorSyst/100)*DR);
		histoCombErrorsForDRSysA2040->SetBinError(i,(binErrorSystA/100)*DR);
		histoCombErrorsForDRSysB2040->SetBinError(i,(binErrorSystB/100)*DR);
		histoCombErrorsForDRSysC2040->SetBinError(i,(binErrorSystC/100)*DR);
		histoCombErrorsForDRStat2040->SetBinError(i,(binErrorStat/100)*DR);
	}
   
	for(Int_t i = 1; i<histoCombErrorsForDRSum2040->GetNbinsX()+1;i++){
		histoCombDirGammaSpectrumErrSum2040->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSys2040->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSysA2040->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSysB2040->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSysC2040->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrStat2040->SetBinContent(i+1,-1);

		histoCombDirGammaSpectrumErrSum2040->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSys2040->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSysA2040->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSysB2040->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSysC2040->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrStat2040->SetBinError(i+1,0);
	}

	// get the binning of the direct photons from the DR
	TH1D *hisoCombDirGammaSpecSysErr2040 				= new TH1D(*histoCombErrorsForDRSys2040);
	TH1D *hisoCombDirGammaSpecSysAErr2040 				= new TH1D(*histoCombErrorsForDRSysA2040);
	TH1D *hisoCombDirGammaSpecSysBErr2040 				= new TH1D(*histoCombErrorsForDRSysB2040);
	TH1D *hisoCombDirGammaSpecSysCErr2040 				= new TH1D(*histoCombErrorsForDRSysC2040);
	TH1D *hisoCombDirGammaSpecStatErr2040 				= new TH1D(*histoCombErrorsForDRStat2040);
	TH1D *hisoCombDirGammaSpecSumErr2040 				= new TH1D(*histoCombErrorsForDRSum2040);

	for(Int_t i = 1; i<graphCombDRPi0FitStatErr2040->GetN()+1; i++){
		// obtain common quantities
		Double_t Rgamma 			= histoCombErrorsForDRSys2040->GetBinContent(i);
		Double_t nIncGamma			= graphCombIncGammaStatErr2040->GetY()[i-1];
		
		// calculating Systematics graph
		Double_t errRgamma			= histoCombErrorsForDRSys2040->GetBinError(i);
		Double_t errNIncGam 		= graphCombIncGammaSysErr2040->GetEYhigh()[i-1];
		Double_t q1 				= 1 - 1/ Rgamma;
		
		Double_t q1Error 			= errRgamma/(Rgamma*Rgamma);
		Double_t content 			= nIncGamma * ( 1 - 1/ Rgamma);
		Double_t error 				= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		Double_t errDR				= content - error;
		hisoCombDirGammaSpecSysErr2040->SetBinError(i, error);
		hisoCombDirGammaSpecSysErr2040->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSys2040->SetBinContent(i, errDR);

		// calculating Systematics A graph
		errRgamma					= histoCombErrorsForDRSysA2040->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSysAErr2040->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSysAErr2040->SetBinError(i, error);
		hisoCombDirGammaSpecSysAErr2040->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSysA2040->SetBinContent(i, errDR);

		// calculating Systematics B graph
		errRgamma					= histoCombErrorsForDRSysB2040->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSysBErr2040->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSysBErr2040->SetBinError(i, error);
		hisoCombDirGammaSpecSysBErr2040->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSysB2040->SetBinContent(i, errDR);

		// calculating Systematics C graph
		errRgamma					= histoCombErrorsForDRSysC2040->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSysCErr2040->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSysCErr2040->SetBinError(i, error);
		hisoCombDirGammaSpecSysCErr2040->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSysC2040->SetBinContent(i, errDR);		
		
		// calculating Stat graphs
		errRgamma					= histoCombErrorsForDRStat2040->GetBinError(i);
		errNIncGam 					= graphCombIncGammaStatErr2040->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecStatErr2040->SetBinError(i, error);
		hisoCombDirGammaSpecStatErr2040->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrStat2040->SetBinContent(i, errDR);
		
		// calculating summed error graphs
		errRgamma					= hisoCombDirGammaSpecSumErr2040->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSumErr2040->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSumErr2040->SetBinError(i, error);
		hisoCombDirGammaSpecSumErr2040->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSum2040->SetBinContent(i, errDR);
	}
	// purely calculating points based on all Systematic errors
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystErr2040 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSys2040,hisoCombDirGammaSpecStatErr2040,0,0.5);
	if(graphCombDirGammaSpectrumSystErr2040)graphCombDirGammaSpectrumSystErr2040->SetName("graphCombDirGammaSpectrumSystErr2040");
	if(graphCombDirGammaSpectrumSystErr2040)graphCombDirGammaSpectrumSystErr2040->Print();
	if(graphCombDirGammaSpectrumSystErr2040)cout << "graph has been found" << endl;
	// purely calculating points based on all Systematic errors A
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystAErr2040 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysA2040,hisoCombDirGammaSpecStatErr2040,0,0.5);
	if(graphCombDirGammaSpectrumSystAErr2040)graphCombDirGammaSpectrumSystAErr2040->SetName("graphCombDirGammaSpectrumSystAErr2040");
	if(graphCombDirGammaSpectrumSystAErr2040)graphCombDirGammaSpectrumSystAErr2040->Print();
	// purely calculating points based on all Systematic errors B
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystBErr2040 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysB2040,hisoCombDirGammaSpecStatErr2040,0,0.5);
	if(graphCombDirGammaSpectrumSystBErr2040)graphCombDirGammaSpectrumSystBErr2040->SetName("graphCombDirGammaSpectrumSystBErr2040");
	if(graphCombDirGammaSpectrumSystBErr2040)graphCombDirGammaSpectrumSystBErr2040->Print();
	// purely calculating points based on all Systematic errors C
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystCErr2040 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysC2040,hisoCombDirGammaSpecStatErr2040,0,0.5);
	if(graphCombDirGammaSpectrumSystCErr2040)graphCombDirGammaSpectrumSystCErr2040->SetName("graphCombDirGammaSpectrumSystCErr2040");
	if(graphCombDirGammaSpectrumSystCErr2040)graphCombDirGammaSpectrumSystCErr2040->Print();

	// purely calculating points based on Statistical errors
	TGraphAsymmErrors *graphCombDirGammaSpectrumStatErr2040 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrStat2040,hisoCombDirGammaSpecStatErr2040,0,0.5);
	if(graphCombDirGammaSpectrumStatErr2040)graphCombDirGammaSpectrumStatErr2040->SetName("graphCombDirGammaSpectrumStatErr2040");
	// purely calculating points based on all Systematic + Statistical errors
	TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr2040 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum2040,hisoCombDirGammaSpecStatErr2040,0,0.5);
	if(graphCombDirGammaSpectrumSumErr2040)graphCombDirGammaSpectrumSumErr2040->SetName("graphCombDirGammaSpectrumSumErr2040");
	// calculate arrows for points with 0, error summed
	TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr2040Ar = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum2040,hisoCombDirGammaSpecStatErr2040,5,0.5);
	if(graphCombDirGammaSpectrumSumErr2040Ar)graphCombDirGammaSpectrumSumErr2040Ar->SetName("graphCombDirGammaSpectrumSumErr2040Ar");
	// calculate points below confidence level summed errors with arrows	
	TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr2040ArConfi = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum2040,hisoCombDirGammaSpecStatErr2040,7,0.5);
	if(graphCombDirGammaSpectrumSumErr2040ArConfi)graphCombDirGammaSpectrumSumErr2040ArConfi->SetName("graphCombDirGammaSpectrumSumErr2040Ar");

	TF1* fitThermalGamma2040Sum						 			= FitObject("e","fitThermalGamma2040Sum","Photon",histoCombDirGammaSpectrumErrSum2040,1.1,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma2040Sum)<< endl;   	
	TF1* fitThermalGamma2040Sum23						 		= FitObject("e","fitThermalGamma2040Sum23","Photon",histoCombDirGammaSpectrumErrSum2040,1.1,2.3,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma2040Sum23)<< endl;   	
	TF1* fitThermalGamma2040Stat						 		= FitObject("e","fitThermalGamma2040Stat","Photon",graphCombDirGammaSpectrumStatErr2040,1.1,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma2040Stat)<< endl;   
	TF1* fitThermalGamma2040Sys						 			= FitObject("e","fitThermalGamma2040Sys","Photon",graphCombDirGammaSpectrumSystErr2040,1.1,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma2040Sys)<< endl;   
	TF1* fitThermalGamma2040SysA						 		= FitObject("e","fitThermalGamma2040SysA","Photon",graphCombDirGammaSpectrumSystAErr2040,1.1,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma2040SysA)<< endl;   
	TF1* fitThermalGamma2040SysB						 		= FitObject("e","fitThermalGamma2040SysB","Photon",graphCombDirGammaSpectrumSystBErr2040,1.1,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma2040SysB)<< endl;   
	TF1* fitThermalGamma2040SysC						 		= FitObject("e","fitThermalGamma2040SysC","Photon",graphCombDirGammaSpectrumSystCErr2040,1.1,2.1,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitThermalGamma2040SysC)<< endl;   

	TF1* fitFullDirGamma2040Sys									= FitObject("qcd","fitFullDirGamma2040Sys","Photon",graphCombDirGammaSpectrumSystErr2040,0.9,14,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitFullDirGamma2040Sys)<< endl;   	
	TF1* fitFullDirGamma2040Stat								= FitObject("qcd","fitFullDirGamma2040Stat","Photon",graphCombDirGammaSpectrumStatErr2040,0.9,14,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitFullDirGamma2040Stat)<< endl;   	
	
	TGraphAsymmErrors* graphRatioCombFitDirGammaStatErr2040	= (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr2040->Clone();
	TGraphAsymmErrors* graphRatioCombFitDirGammaSysErr2040 	= (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr2040->Clone();
	graphRatioCombFitDirGammaStatErr2040 					= CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaStatErr2040, fitFullDirGamma2040Stat); 
	graphRatioCombFitDirGammaSysErr2040 					= CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaSysErr2040, fitFullDirGamma2040Sys); 
	
	TGraphAsymmErrors* graphCombThermalGammaSpectrumStatErr2040 	= SubtractPromptPhotonsViaFit( 	fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumStatErr2040);
	TGraphAsymmErrors* graphCombThermalGammaSpectrumSysErr2040 		= SubtractPromptPhotonsViaFit( 	fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumSystErr2040);
	TGraphAsymmErrors* graphCombThermalGammaSpectrumSumErr2040Ar	= SubtractPromptPhotonsViaFit( 	fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumSumErr2040Ar, newBinsComb, 20);
	
	cout << endl << "Calculating RAA" << endl;
	TGraphAsymmErrors* graphCombRAADirGammaStat2040 	= NULL;
	TGraphAsymmErrors* graphCombRAADirGammaSys2040		= NULL;
	TGraphAsymmErrors* graphCombRAADirGammaSum2040Ar	= NULL;
	if (graphCombDirGammaSpectrumStatErr2040) 	CalcRaaWithTheoryFit( fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumStatErr2040, &graphCombRAADirGammaStat2040);
	if (graphCombDirGammaSpectrumSystErr2040) 	CalcRaaWithTheoryFit( fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumSystErr2040, &graphCombRAADirGammaSys2040);
	if (graphCombDirGammaSpectrumSumErr2040Ar) 	CalcRaaWithTheoryFit( fitTheoryPromptMcGill2040, graphCombDirGammaSpectrumSumErr2040Ar, &graphCombRAADirGammaSum2040Ar, newBinsComb, 20);
	if (graphCombRAADirGammaStat2040) graphCombRAADirGammaStat2040->Print();
	if (graphCombRAADirGammaSys2040) graphCombRAADirGammaSys2040->Print();
	if (graphCombRAADirGammaSum2040Ar) graphCombRAADirGammaSum2040Ar->Print();
	
	//*******************************************************************************************************************************************
	//*********************************************** Calculate direct photon spectrum 4080 *****************************************************
	//*******************************************************************************************************************************************	
	//_______________________ copy inclusive photon spectra _____________________
	//_______________________ copy inclusive photon spectra _____________________
	TH1D *histoCombDirGammaSpectrumErrSum4080 					= new TH1D("histoCombDirGammaSpectrumErrSum4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSys4080 					= new TH1D("histoCombDirGammaSpectrumErrSys4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSysA4080 					= new TH1D("histoCombDirGammaSpectrumErrSysA4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSysB4080 					= new TH1D("histoCombDirGammaSpectrumErrSysB4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrSysC4080 					= new TH1D("histoCombDirGammaSpectrumErrSysC4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	TH1D *histoCombDirGammaSpectrumErrStat4080 					= new TH1D("histoCombDirGammaSpectrumErrStat4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);

	//_______________________ get arrays of double ratio errors __________________
	Double_t *SystErrorsCombDR4080 								= new Double_t[graphCombIncGammaStatErr4080->GetN()];
	Double_t *SystAErrorsCombDR4080 								= new Double_t[graphCombIncGammaStatErr4080->GetN()];
	Double_t *SystBErrorsCombDR4080 								= new Double_t[graphCombIncGammaStatErr4080->GetN()];
	Double_t *SystCErrorsCombDR4080 								= new Double_t[graphCombIncGammaStatErr4080->GetN()];
	Double_t *sumErrorsCombDR4080 								= new Double_t[graphCombIncGammaStatErr4080->GetN()];
	Double_t *StatErrorsCombDR4080 								= new Double_t[graphCombIncGammaStatErr4080->GetN()];
	Double_t *xErrorsDR4080 									= new Double_t[graphCombIncGammaStatErr4080->GetN()];
	graphCombIncGammaSysErr4080->Print();
	for (Int_t i = 0; i< graphCombDRPi0FitStatErr4080->GetN(); i++){
		SystErrorsCombDR4080[i] 								= graphCombDRPi0FitSysErr4080->GetEYhigh()[i]/graphCombDRPi0FitSysErr4080->GetY()[i] *100;
		SystAErrorsCombDR4080[i] 								= graphCombDRPi0FitSysAErr4080->GetEYhigh()[i]/graphCombDRPi0FitSysAErr4080->GetY()[i] *100;
		SystBErrorsCombDR4080[i] 								= graphCombDRPi0FitSysBErr4080->GetEYhigh()[i]/graphCombDRPi0FitSysBErr4080->GetY()[i] *100;
		SystCErrorsCombDR4080[i] 								= graphCombDRPi0FitSysCErr4080->GetEYhigh()[i]/graphCombDRPi0FitSysCErr4080->GetY()[i] *100;
		StatErrorsCombDR4080[i] 								= graphCombDRPi0FitStatErr4080->GetEYhigh()[i]/graphCombDRPi0FitStatErr4080->GetY()[i] *100;
		sumErrorsCombDR4080[i] 									= graphCombDRPi0FitSumErr4080->GetEYhigh()[i]/graphCombDRPi0FitSumErr4080->GetY()[i] *100;
// 		cout << i << "\t" << graphCombDRPi0FitSysErr4080->GetY()[i] << "\t" << graphCombDRPi0FitSysErr4080->GetEYhigh()[i] << "\t" <<SystErrorsCombDR4080[i] << endl;
	}
	xErrorsDR4080 												= graphCombDRPi0FitStatErr4080->GetX();

	cout << "here !!! \n\n" << endl;
	graphCombDRPi0FitSumErr4080->Print();
	
	
	//_______________________ copy inclusive photon spectra _____________________	
	TH1D* histoCombErrorsForDRSum4080 							= new TH1D("histoCombErrorsForDRSum4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRStat4080 							= new TH1D("histoCombErrorsForDRStat4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSys4080 							= new TH1D("histoCombErrorsForDRSys4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSysA4080 							= new TH1D("histoCombErrorsForDRSysA4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSysB4080 							= new TH1D("histoCombErrorsForDRSysB4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	TH1D* histoCombErrorsForDRSysC4080 							= new TH1D("histoCombErrorsForDRSysC4080","",graphCombDRPi0FitStatErr4080->GetN(),xArrayCombined);
	
	for(Int_t i = 1; i<graphCombDRPi0FitStatErr4080->GetN()+1;i++){
		cout<< i << "\t"<<xErrorsDR4080[i-1]<<"  "<<histoCombErrorsForDRSum4080->GetBinCenter(i)<< "\t"<<histoCombErrorsForDRSum4080->GetBinWidth(i) <<endl;
		Double_t binErrorSummed 								= sumErrorsCombDR4080[i-1];
		Double_t binErrorSyst 									= SystErrorsCombDR4080[i-1];
		Double_t binErrorSystA 									= SystAErrorsCombDR4080[i-1];
		Double_t binErrorSystB 									= SystBErrorsCombDR4080[i-1];
		Double_t binErrorSystC 									= SystCErrorsCombDR4080[i-1];
		Double_t binErrorStat 									= StatErrorsCombDR4080[i-1];
		Double_t DR 											= graphCombDRPi0FitStatErr4080->GetY()[i-1];

		cout << DR << "\t" << binErrorStat << "\t" << binErrorSyst << "\t" << binErrorSummed << endl;
		histoCombErrorsForDRSum4080->SetBinContent(i,DR);
		histoCombErrorsForDRSys4080->SetBinContent(i,DR);
		histoCombErrorsForDRSysA4080->SetBinContent(i,DR);
		histoCombErrorsForDRSysB4080->SetBinContent(i,DR);
		histoCombErrorsForDRSysC4080->SetBinContent(i,DR);
		histoCombErrorsForDRStat4080->SetBinContent(i,DR);
		histoCombErrorsForDRSum4080->SetBinError(i,(binErrorSummed/100)*DR);
		histoCombErrorsForDRSys4080->SetBinError(i,(binErrorSyst/100)*DR);
		histoCombErrorsForDRSysA4080->SetBinError(i,(binErrorSystA/100)*DR);
		histoCombErrorsForDRSysB4080->SetBinError(i,(binErrorSystB/100)*DR);
		histoCombErrorsForDRSysC4080->SetBinError(i,(binErrorSystC/100)*DR);
		histoCombErrorsForDRStat4080->SetBinError(i,(binErrorStat/100)*DR);
	}
   
	for(Int_t i = 1; i<histoCombErrorsForDRSum4080->GetNbinsX()+1;i++){
		histoCombDirGammaSpectrumErrSum4080->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSys4080->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSysA4080->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSysB4080->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrSysC4080->SetBinContent(i+1,-1);
		histoCombDirGammaSpectrumErrStat4080->SetBinContent(i+1,-1);

		histoCombDirGammaSpectrumErrSum4080->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSys4080->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSysA4080->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSysB4080->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrSysC4080->SetBinError(i+1,0);
		histoCombDirGammaSpectrumErrStat4080->SetBinError(i+1,0);
	}

	// get the binning of the direct photons from the DR
	TH1D *hisoCombDirGammaSpecSysErr4080 				= new TH1D(*histoCombErrorsForDRSys4080);
	TH1D *hisoCombDirGammaSpecSysAErr4080 				= new TH1D(*histoCombErrorsForDRSysA4080);
	TH1D *hisoCombDirGammaSpecSysBErr4080 				= new TH1D(*histoCombErrorsForDRSysB4080);
	TH1D *hisoCombDirGammaSpecSysCErr4080 				= new TH1D(*histoCombErrorsForDRSysC4080);
	TH1D *hisoCombDirGammaSpecStatErr4080 				= new TH1D(*histoCombErrorsForDRStat4080);
	TH1D *hisoCombDirGammaSpecSumErr4080 				= new TH1D(*histoCombErrorsForDRSum4080);

	for(Int_t i = 1; i<graphCombDRPi0FitStatErr4080->GetN()+1; i++){
		// obtain common quantities
		Double_t Rgamma 			= histoCombErrorsForDRSys4080->GetBinContent(i);
		Double_t nIncGamma			= graphCombIncGammaStatErr4080->GetY()[i-1];
		
		// calculating Systematics graph
		Double_t errRgamma			= histoCombErrorsForDRSys4080->GetBinError(i);
		Double_t errNIncGam 		= graphCombIncGammaSysErr4080->GetEYhigh()[i-1];
		Double_t q1 				= 1 - 1/ Rgamma;
		
		Double_t q1Error 			= errRgamma/(Rgamma*Rgamma);
		Double_t content 			= nIncGamma * ( 1 - 1/ Rgamma);
		Double_t error 				= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		Double_t errDR				= content - error;
		hisoCombDirGammaSpecSysErr4080->SetBinError(i, error);
		hisoCombDirGammaSpecSysErr4080->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSys4080->SetBinContent(i, errDR);

		// calculating Systematics A graph
		errRgamma					= histoCombErrorsForDRSysA4080->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSysAErr4080->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSysAErr4080->SetBinError(i, error);
		hisoCombDirGammaSpecSysAErr4080->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSysA4080->SetBinContent(i, errDR);

		// calculating Systematics B graph
		errRgamma					= histoCombErrorsForDRSysB4080->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSysBErr4080->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSysBErr4080->SetBinError(i, error);
		hisoCombDirGammaSpecSysBErr4080->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSysB4080->SetBinContent(i, errDR);

		// calculating Systematics C graph
		errRgamma					= histoCombErrorsForDRSysC4080->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSysCErr4080->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSysCErr4080->SetBinError(i, error);
		hisoCombDirGammaSpecSysCErr4080->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSysC4080->SetBinContent(i, errDR);		
		
		// calculating Stat graphs
		errRgamma					= histoCombErrorsForDRStat4080->GetBinError(i);
		errNIncGam 					= graphCombIncGammaStatErr4080->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecStatErr4080->SetBinError(i, error);
		hisoCombDirGammaSpecStatErr4080->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrStat4080->SetBinContent(i, errDR);
		
		// calculating summed error graphs
		errRgamma					= hisoCombDirGammaSpecSumErr4080->GetBinError(i);
		errNIncGam 					= graphCombIncGammaSumErr4080->GetEYhigh()[i-1];
		q1 							= 1 - 1/ Rgamma;
		q1Error 					= errRgamma/(Rgamma*Rgamma);
		content 					= nIncGamma * ( 1 - 1/ Rgamma);
		error 						= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR						= content - error;
		hisoCombDirGammaSpecSumErr4080->SetBinError(i, error);
		hisoCombDirGammaSpecSumErr4080->SetBinContent(i, content);
		histoCombDirGammaSpectrumErrSum4080->SetBinContent(i, errDR);
	}
	// purely calculating points based on all Systematic errors
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystErr4080 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSys4080,hisoCombDirGammaSpecStatErr4080,0,0.5);
	if(graphCombDirGammaSpectrumSystErr4080)graphCombDirGammaSpectrumSystErr4080->SetName("graphCombDirGammaSpectrumSystErr4080");
	if(graphCombDirGammaSpectrumSystErr4080)graphCombDirGammaSpectrumSystErr4080->Print();
	if(graphCombDirGammaSpectrumSystErr4080)cout << "graph has been found" << endl;
	// purely calculating points based on all Systematic errors A
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystAErr4080 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysA4080,hisoCombDirGammaSpecStatErr4080,0,0.5);
	if(graphCombDirGammaSpectrumSystAErr4080)graphCombDirGammaSpectrumSystAErr4080->SetName("graphCombDirGammaSpectrumSystAErr4080");
	if(graphCombDirGammaSpectrumSystAErr4080)graphCombDirGammaSpectrumSystAErr4080->Print();
	// purely calculating points based on all Systematic errors B
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystBErr4080 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysB4080,hisoCombDirGammaSpecStatErr4080,0,0.5);
	if(graphCombDirGammaSpectrumSystBErr4080)graphCombDirGammaSpectrumSystBErr4080->SetName("graphCombDirGammaSpectrumSystBErr4080");
	if(graphCombDirGammaSpectrumSystBErr4080)graphCombDirGammaSpectrumSystBErr4080->Print();
	// purely calculating points based on all Systematic errors C
	TGraphAsymmErrors *graphCombDirGammaSpectrumSystCErr4080 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSysC4080,hisoCombDirGammaSpecStatErr4080,0,0.5);
	if(graphCombDirGammaSpectrumSystCErr4080)graphCombDirGammaSpectrumSystCErr4080->SetName("graphCombDirGammaSpectrumSystCErr4080");
	if(graphCombDirGammaSpectrumSystCErr4080)graphCombDirGammaSpectrumSystCErr4080->Print();

	// purely calculating points based on Statistical errors
	TGraphAsymmErrors *graphCombDirGammaSpectrumStatErr4080 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrStat4080,hisoCombDirGammaSpecStatErr4080,0,0.5);
	if(graphCombDirGammaSpectrumStatErr4080)graphCombDirGammaSpectrumStatErr4080->SetName("graphCombDirGammaSpectrumStatErr4080");
	// purely calculating points based on all Systematic + Statistical errors
	TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr4080 = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum4080,hisoCombDirGammaSpecStatErr4080,0,0.5);
	if(graphCombDirGammaSpectrumSumErr4080)graphCombDirGammaSpectrumSumErr4080->SetName("graphCombDirGammaSpectrumSumErr4080");
 	// calculate arrows for points with 0, error summed
	TGraphAsymmErrors *graphCombDirGammaSpectrumSumErr4080Ar = CalculateDirectPhotonPointsAndUpperLimits(histoCombDirGammaSpectrumErrSum4080,hisoCombDirGammaSpecStatErr4080,5,0.5);
	if(graphCombDirGammaSpectrumSumErr4080Ar)graphCombDirGammaSpectrumSumErr4080Ar->SetName("graphCombDirGammaSpectrumSumErr4080Ar");

	TF1* fitFullDirGamma4080Sys									= FitObject("qcd","fitFullDirGamma4080Sys","Photon",graphCombDirGammaSpectrumSystErr4080,0.9,14,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitFullDirGamma4080Sys)<< endl;   	
	TF1* fitFullDirGamma4080Stat								= FitObject("qcd","fitFullDirGamma4080Stat","Photon",graphCombDirGammaSpectrumStatErr4080,0.9,14,NULL,"QNRMEX0+");
	fileFinalResults << WriteParameterToFile(fitFullDirGamma4080Stat)<< endl;   	
	
	TGraphAsymmErrors* graphRatioCombFitDirGammaStatErr4080	= (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr4080->Clone();
	TGraphAsymmErrors* graphRatioCombFitDirGammaSysErr4080 	= (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr4080->Clone();
	graphRatioCombFitDirGammaStatErr4080 					= CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaStatErr4080, fitFullDirGamma4080Stat); 
	graphRatioCombFitDirGammaSysErr4080 					= CalculateGraphErrRatioToFit(graphRatioCombFitDirGammaSysErr4080, fitFullDirGamma4080Sys); 

	
	cout << endl << "Calculating RAA" << endl;
	TGraphAsymmErrors* graphCombRAADirGammaStat4080 	= NULL;
	TGraphAsymmErrors* graphCombRAADirGammaSys4080		= NULL;
	TGraphAsymmErrors* graphCombRAADirGammaSum4080Ar	= NULL;
	if (graphCombDirGammaSpectrumStatErr4080) 	CalcRaaWithTheoryFit( fitTheoryPromptMcGill4080, graphCombDirGammaSpectrumStatErr4080, &graphCombRAADirGammaStat4080);
	if (graphCombDirGammaSpectrumSystErr4080) 	CalcRaaWithTheoryFit( fitTheoryPromptMcGill4080, graphCombDirGammaSpectrumSystErr4080, &graphCombRAADirGammaSys4080);
	if (graphCombDirGammaSpectrumSumErr4080Ar) 	CalcRaaWithTheoryFit( fitTheoryPromptMcGill4080, graphCombDirGammaSpectrumSumErr4080Ar, &graphCombRAADirGammaSum4080Ar, newBinsComb, 20);
	if (graphCombRAADirGammaStat4080) graphCombRAADirGammaStat4080->Print();
	if (graphCombRAADirGammaSys4080) graphCombRAADirGammaSys4080->Print();
	if (graphCombRAADirGammaSum4080Ar) graphCombRAADirGammaSum4080Ar->Print();


	
	//*******************************************************************************************************************************************
	//******************************************* DR plot with individual measurements **********************************************************
	//*******************************************************************************************************************************************
	Double_t arrayBoundsXIndMeasRatio[2];
	Double_t arrayBoundsYIndMeasRatio[4];
	Double_t relativeMarginsIndMeasRatioX[3];
	Double_t relativeMarginsIndMeasRatioY[3];
	ReturnCorrectValuesForCanvasScaling(1200, 1400, 1, 3, 0.09, 0.01, 0.01, 0.065, arrayBoundsXIndMeasRatio, arrayBoundsYIndMeasRatio, relativeMarginsIndMeasRatioX, relativeMarginsIndMeasRatioY);
		
	TCanvas * canvasRatioIndDR = new TCanvas("canvasRatioIndDR","",10,10,1200,1400);  // gives the page size		
	canvasRatioIndDR->cd();
	
	TPad* padPartRatioInDR1 = new TPad("padPartRatioInDR1", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( padPartRatioInDR1, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	padPartRatioInDR1->Draw();
	TPad* padPartRatioInDR2 = new TPad("padPartRatioInDR2", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( padPartRatioInDR2, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[1]);
	padPartRatioInDR2->Draw();
	TPad* padPartRatioInDR3 = new TPad("padPartRatioInDR3", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[3], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[2],-1, -1, -2);
	DrawGammaPadSettings( padPartRatioInDR3, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	padPartRatioInDR3->Draw();

	//_______________________________________________________________ define text sizes _________________________________________________________
	Int_t textSizeLabelsPixel = 43;
	Double_t margin = relativeMarginsIndMeasRatioX[0]*1200;
	Double_t textsizeLabelsPad1 = 0;
	Double_t textsizeFacPad1 = 0;
	if (padPartRatioInDR1->XtoPixel(padPartRatioInDR1->GetX2()) < padPartRatioInDR1->YtoPixel(padPartRatioInDR1->GetY1())){
		textsizeLabelsPad1 = (Double_t)textSizeLabelsPixel/padPartRatioInDR1->XtoPixel(padPartRatioInDR1->GetX2()) ;
		textsizeFacPad1 = (Double_t)1./padPartRatioInDR1->XtoPixel(padPartRatioInDR1->GetX2()) ;
	} else {
		textsizeLabelsPad1 = (Double_t)textSizeLabelsPixel/padPartRatioInDR1->YtoPixel(padPartRatioInDR1->GetY1());
		textsizeFacPad1 = (Double_t)1./padPartRatioInDR1->YtoPixel(padPartRatioInDR1->GetY1());
	}
	Double_t textsizeLabelsPad2 = 0;
	Double_t textsizeFacPad2 = 0;
	if (padPartRatioInDR2->XtoPixel(padPartRatioInDR2->GetX2()) <padPartRatioInDR2->YtoPixel(padPartRatioInDR2->GetY1()) ){
		textsizeLabelsPad2 = (Double_t)textSizeLabelsPixel/padPartRatioInDR2->XtoPixel(padPartRatioInDR2->GetX2()) ;
		textsizeFacPad2 = (Double_t)1./padPartRatioInDR2->XtoPixel(padPartRatioInDR2->GetX2()) ;
	} else {
		textsizeLabelsPad2 = (Double_t)textSizeLabelsPixel/padPartRatioInDR2->YtoPixel(padPartRatioInDR2->GetY1());
		textsizeFacPad2 = (Double_t)1./padPartRatioInDR2->YtoPixel(padPartRatioInDR2->GetY1());
	}
	Double_t textsizeLabelsPad3 = 0;
	Double_t textsizeFacPad3= 0;
	if (padPartRatioInDR3->XtoPixel(padPartRatioInDR3->GetX2()) <padPartRatioInDR3->YtoPixel(padPartRatioInDR3->GetY1()) ){
		textsizeLabelsPad3 = (Double_t)textSizeLabelsPixel/padPartRatioInDR3->XtoPixel(padPartRatioInDR3->GetX2()) ;
		textsizeFacPad3 = (Double_t)1./padPartRatioInDR3->XtoPixel(padPartRatioInDR3->GetX2()) ;
	} else {
		textsizeLabelsPad3 = (Double_t)textSizeLabelsPixel/padPartRatioInDR3->YtoPixel(padPartRatioInDR3->GetY1());
		textsizeFacPad3 = (Double_t)1./padPartRatioInDR3->YtoPixel(padPartRatioInDR3->GetY1());
	}

	//_______________________________________________________________ 0-20% dummy upper panel ___________________________________________________
	TH2D *dummyDR1 ;
	dummyDR1 = new TH2D("dummyDR1", "dummyDR1", 1000, 0., 22, 1000., doubleRatio[0], doubleRatio[1]);
	SetStyleHistoTH2ForGraphs( dummyDR1, "#it{p}_{T} (GeV/#it{c})", "", 
							   0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.95,0.10/(textsizeFacPad1*margin), 510, 505);
	dummyDR1->GetXaxis()->SetLabelOffset(-0.015);
	dummyDR1->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyDR1->GetXaxis()->SetTickLength(0.06);
	dummyDR1->GetYaxis()->SetTickLength(0.028);

	//_______________________________________________________________ 20-40% dummy middle panel _________________________________________________
	TH2D *dummyDR2 ;
	dummyDR2 = new TH2D("dummyDR2", "dummyDR2", 1000, 0., 22, 1000., doubleRatio[0], doubleRatio[1]);
	SetStyleHistoTH2ForGraphs( dummyDR2, "#it{p}_{T} (GeV/#it{c})", "#it{R}_{#gamma}", // = (#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}}) 
							   0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.95,0.10/(textsizeFacPad2*margin), 510, 505);
	dummyDR2->GetXaxis()->SetLabelOffset(-0.015);
	dummyDR2->GetYaxis()->CenterTitle(kTRUE);
	dummyDR2->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyDR2->GetXaxis()->SetTickLength(0.06);
	dummyDR2->GetYaxis()->SetTickLength(0.028);

	//_______________________________________________________________ 40-80% dummy lower panel __________________________________________________
	TH2D *dummyDR3 ;
	dummyDR3 = new TH2D("dummyDR3", "dummyDR3", 1000, 0., 22, 1000., doubleRatio[0], doubleRatio[1]);
	SetStyleHistoTH2ForGraphs( dummyDR3, "#it{p}_{T} (GeV/#it{c})", "", 
							   0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.92,0.10/(textsizeFacPad3*margin), 510, 505);
	dummyDR3->GetXaxis()->SetLabelOffset(-0.015);
	dummyDR3->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyDR3->GetXaxis()->SetTickLength(0.055);
	dummyDR3->GetYaxis()->SetTickLength(0.035);
	
	//_______________________________________________________________ 0-20% panel _______________________________________________________________
	padPartRatioInDR1->cd();
	padPartRatioInDR1->SetLogx(1);
		dummyDR1->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr0020, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		graphPCMDRPi0FitSysErr0020->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr0020, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		graphPHOSDRPi0FitSysErr0020->Draw("E2same");

		
		DrawGammaSetMarker(histoPCMDRPi0FitStatErr0020, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		histoPCMDRPi0FitStatErr0020->Draw("p,same,e0,X0");		
		DrawGammaSetMarker(histoPHOSDRPi0FitStatErr0020, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		histoPHOSDRPi0FitStatErr0020->Draw("p,same,e0,X0");
		
		TLatex *labelDRCent0020 = new TLatex(0.12,0.85,collisionSystemCent0020.Data());
		SetStyleTLatex( labelDRCent0020, 0.85*textsizeLabelsPad1,4);
		labelDRCent0020->Draw();

		TLatex *labelALICEInd = new TLatex(0.82,0.85,"ALICE");
		SetStyleTLatex( labelALICEInd, 0.85*textsizeLabelsPad1,4);
		labelALICEInd->Draw();

		TLegend* legendDRIndMeas = new TLegend(0.15,0.6,0.4,0.78);
		legendDRIndMeas->SetFillStyle(0);
		legendDRIndMeas->SetFillColor(0);
		legendDRIndMeas->SetLineColor(0);
		legendDRIndMeas->SetTextSize(0.85*textsizeLabelsPad1);
		legendDRIndMeas->SetMargin(0.2);
		legendDRIndMeas->SetTextFont(42);
		legendDRIndMeas->AddEntry(graphPCMDRPi0FitSysErr0020, "PCM","pf");
		legendDRIndMeas->AddEntry(graphPHOSDRPi0FitSysErr0020,"PHOS","pf");
		legendDRIndMeas->Draw();
	
	//_______________________________________________________________ 20-40% panel _______________________________________________________________
	padPartRatioInDR2->cd();
	padPartRatioInDR2->SetLogx(1);
		dummyDR2->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
	
		DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2040, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		graphPCMDRPi0FitSysErr2040->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		graphPHOSDRPi0FitSysErr2040->Draw("E2same");

		DrawGammaSetMarker(histoPCMDRPi0FitStatErr2040, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		histoPCMDRPi0FitStatErr2040->Draw("p,same,e0,X0");
		DrawGammaSetMarker(histoPHOSDRPi0FitStatErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		histoPHOSDRPi0FitStatErr2040->Draw("p,same,e0,X0");

		TLatex *labelDRCent2040 = new TLatex(0.12,0.88,collisionSystemCent2040.Data());
		SetStyleTLatex( labelDRCent2040, 0.85*textsizeLabelsPad2,4);
		labelDRCent2040->Draw();
	
	//_______________________________________________________________ 40-80% panel _______________________________________________________________	
	padPartRatioInDR3->cd();
	padPartRatioInDR3->SetLogx(1);
		dummyDR3->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);	
		
		DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr4080, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		graphPCMDRPi0FitSysErr4080->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr4080, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		graphPHOSDRPi0FitSysErr4080->Draw("E2same");
		
		DrawGammaSetMarker(histoPCMDRPi0FitStatErr4080, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		histoPCMDRPi0FitStatErr4080->Draw("p,same,e0,X0");
		DrawGammaSetMarker(histoPHOSDRPi0FitStatErr4080, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		histoPHOSDRPi0FitStatErr4080->Draw("p,same,e0,X0");

		TLatex *labelDRCent4080 = new TLatex(0.12,0.9,collisionSystemCent4080.Data());
		SetStyleTLatex( labelDRCent4080, 0.85*textsizeLabelsPad3,4);
		labelDRCent4080->Draw();
		
	canvasRatioIndDR->SaveAs(Form("%s/DR_individualMeasurements.%s", outputDir.Data(), suffix.Data()));	

	//*******************************************************************************************************************************************
	//********************************************** DR plot with PCM measurement + NLO  ********************************************************
	//*******************************************************************************************************************************************	
	canvasRatioIndDR->cd();
	
	padPartRatioInDR1->Draw();
	padPartRatioInDR2->Draw();
	padPartRatioInDR3->Draw();

	
	//_______________________________________________________________ 0-20% panel _______________________________________________________________
	padPartRatioInDR1->cd();
	padPartRatioInDR1->SetLogx(1);
		dummyDR1->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr0020, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020,widthLinesBoxes, kTRUE);
		DrawGammaSetMarker(histoPCMDRPi0FitStatErr0020, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);
		SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR0020, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
		graphTheoryNLODR0020->Draw("pE3lsame");
		graphPCMDRPi0FitSysErr0020->Draw("E2same");
		histoPCMDRPi0FitStatErr0020->Draw("p,same,e0,X0");

		
		
		TLegend* legendDRPCMNLO0020 = new TLegend(0.15,0.82-1.1*0.85*textsizeLabelsPad1*3,0.5,0.82);
		legendDRPCMNLO0020->SetFillStyle(0);
		legendDRPCMNLO0020->SetFillColor(0);
		legendDRPCMNLO0020->SetLineColor(0);
		legendDRPCMNLO0020->SetTextSize(0.85*textsizeLabelsPad1);
		legendDRPCMNLO0020->SetMargin(0.2);
		legendDRPCMNLO0020->SetTextFont(42);
		legendDRPCMNLO0020->AddEntry(graphPCMDRPi0FitSysErr0020,"PCM","pf");
		legendDRPCMNLO0020->AddEntry(graphTheoryNLODR0020,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
		legendDRPCMNLO0020->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
		legendDRPCMNLO0020->Draw();

		labelDRCent0020->Draw();
		
	//_______________________________________________________________ 20-40% panel _______________________________________________________________
	padPartRatioInDR2->cd();
	padPartRatioInDR2->SetLogx(1);
		dummyDR2->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

		DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
		DrawGammaSetMarker(histoPCMDRPi0FitStatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
		SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR2040, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
		graphTheoryNLODR2040->Draw("p3lsame");
		graphPCMDRPi0FitSysErr2040->Draw("E2same");
		histoPCMDRPi0FitStatErr2040->Draw("p,same,e0,X0");

		TLegend* legendDRPCMNLO2040 = new TLegend(0.15,0.85-1.1*0.85*textsizeLabelsPad2*3,0.5,0.85);
		legendDRPCMNLO2040->SetFillStyle(0);
		legendDRPCMNLO2040->SetFillColor(0);
		legendDRPCMNLO2040->SetLineColor(0);
		legendDRPCMNLO2040->SetTextSize(0.85*textsizeLabelsPad2);
		legendDRPCMNLO2040->SetMargin(0.2);
		legendDRPCMNLO2040->SetTextFont(42);
		legendDRPCMNLO2040->AddEntry(graphPCMDRPi0FitSysErr2040,"PCM","pf");
		legendDRPCMNLO2040->AddEntry(graphTheoryNLODR2040,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
		legendDRPCMNLO2040->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
		legendDRPCMNLO2040->Draw();		
		
		labelDRCent2040->Draw();
	
	//_______________________________________________________________ 40-80% panel _______________________________________________________________	
	padPartRatioInDR3->cd();
	padPartRatioInDR3->SetLogx(1);
		dummyDR3->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);	
		
		DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErr4080, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080,widthLinesBoxes, kTRUE);
		DrawGammaSetMarker(histoPCMDRPi0FitStatErr4080, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);
		SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR4080, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
		graphTheoryNLODR4080->Draw("p3lsame");
		graphPCMDRPi0FitSysErr4080->Draw("E2same");
		histoPCMDRPi0FitStatErr4080->Draw("p,same,e0,X0");

		TLegend* legendDRPCMNLO4080 = new TLegend(0.15,0.87-1.1*0.85*textsizeLabelsPad3*3,0.5,0.87);
		legendDRPCMNLO4080->SetFillStyle(0);
		legendDRPCMNLO4080->SetFillColor(0);
		legendDRPCMNLO4080->SetLineColor(0);
		legendDRPCMNLO4080->SetTextSize(0.85*textsizeLabelsPad3);
		legendDRPCMNLO4080->SetMargin(0.2);
		legendDRPCMNLO4080->SetTextFont(42);
		legendDRPCMNLO4080->AddEntry(graphPCMDRPi0FitSysErr4080,"PCM","pf");
		legendDRPCMNLO4080->AddEntry(graphTheoryNLODR4080,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
		legendDRPCMNLO4080->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
		legendDRPCMNLO4080->Draw();

		labelDRCent4080->Draw();
		
	canvasRatioIndDR->SaveAs(Form("%s/DR_PCMMeasurementTheory.%s", outputDir.Data(), suffix.Data()));	

	
	//*******************************************************************************************************************************************
	//********************************************** DR plot with PCM measurements all cents in one plot  ***************************************
	//*******************************************************************************************************************************************	
	TCanvas *canvasDoubleRatio = GetAndSetCanvas("canvasDoubleRatioFinal");
	canvasDoubleRatio->SetLogx();

	TH2D *dummyDR ;
	dummyDR = new TH2D("dummyDR", "dummyDR", 1000, 0., 16, 1000., doubleRatio[0], doubleRatio[1]);
	SetStyleHistoTH2ForGraphs( dummyDR, "#it{p}_{T} (GeV/#it{c})", "#it{R}_{#gamma} = (#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})", 
							   0.045, 0.05, 0.045, 0.05, 0.85, 0.85);
	dummyDR->GetXaxis()->SetLabelOffset(-0.015);
	dummyDR->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyDR->DrawCopy();

		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);	

		graphPCMDRPi0FitSysErr0020->Draw("E2same");
		histoPCMDRPi0FitStatErr0020->Draw("p,same,e0,X0");
		graphPCMDRPi0FitSysErr2040->Draw("E2same");
		histoPCMDRPi0FitStatErr2040->Draw("p,same,e0,X0");
		graphPCMDRPi0FitSysErr4080->Draw("E2same");
		histoPCMDRPi0FitStatErr4080->Draw("p,same,e0,X0");

	
		TLegend* legendDRAllCentPCM = GetAndSetLegend(0.15,0.75,4);
		legendDRAllCentPCM->AddEntry(graphPCMDRPi0FitSysErr0020,"  0-20% PCM","pf");
		legendDRAllCentPCM->AddEntry(graphPCMDRPi0FitSysErr2040,"20-40% PCM","pf");
		legendDRAllCentPCM->AddEntry(graphPCMDRPi0FitSysErr4080,"40-80% PCM","pf");
		
		legendDRAllCentPCM->Draw();

	canvasDoubleRatio->Print(Form("%s/DR_PCMMeasurementAllCentsInOne.%s", outputDir.Data(), suffix.Data()));	


	//*******************************************************************************************************************************************
	//********************************************** DR plot with PHOS measurement + NLO  ********************************************************
	//*******************************************************************************************************************************************	
	canvasRatioIndDR->cd();
	
	padPartRatioInDR1->Draw();
	padPartRatioInDR2->Draw();
	padPartRatioInDR3->Draw();

	
	//_______________________________________________________________ 0-20% panel _______________________________________________________________
	padPartRatioInDR1->cd();
	padPartRatioInDR1->SetLogx(1);
		dummyDR1->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr0020, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020,widthLinesBoxes, kTRUE);
		DrawGammaSetMarker(histoPHOSDRPi0FitStatErr0020, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);
		SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR0020, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
		graphTheoryNLODR0020->Draw("pE3lsame");
		graphPHOSDRPi0FitSysErr0020->Draw("E2same");
		histoPHOSDRPi0FitStatErr0020->Draw("p,same,e0,X0");
		
		TLegend* legendDRPHOSNLO0020 = new TLegend(0.15,0.82-1.1*0.85*textsizeLabelsPad1*3,0.5,0.82);
		legendDRPHOSNLO0020->SetFillStyle(0);
		legendDRPHOSNLO0020->SetFillColor(0);
		legendDRPHOSNLO0020->SetLineColor(0);
		legendDRPHOSNLO0020->SetTextSize(0.85*textsizeLabelsPad1);
		legendDRPHOSNLO0020->SetMargin(0.2);
		legendDRPHOSNLO0020->SetTextFont(42);
		legendDRPHOSNLO0020->AddEntry(graphPHOSDRPi0FitSysErr0020,"PHOS","pf");
		legendDRPHOSNLO0020->AddEntry(graphTheoryNLODR0020,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
		legendDRPHOSNLO0020->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
		legendDRPHOSNLO0020->Draw();

		labelDRCent0020->Draw();
		
	//_______________________________________________________________ 20-40% panel _______________________________________________________________
	padPartRatioInDR2->cd();
	padPartRatioInDR2->SetLogx(1);
		dummyDR2->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

		DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
		DrawGammaSetMarker(histoPHOSDRPi0FitStatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
		SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR2040, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
		graphTheoryNLODR2040->Draw("p3lsame");
		graphPHOSDRPi0FitSysErr2040->Draw("E2same");
		histoPHOSDRPi0FitStatErr2040->Draw("p,same,e0,X0");

		TLegend* legendDRPHOSNLO2040 = new TLegend(0.15,0.85-1.1*0.85*textsizeLabelsPad2*3,0.5,0.85);
		legendDRPHOSNLO2040->SetFillStyle(0);
		legendDRPHOSNLO2040->SetFillColor(0);
		legendDRPHOSNLO2040->SetLineColor(0);
		legendDRPHOSNLO2040->SetTextSize(0.85*textsizeLabelsPad2);
		legendDRPHOSNLO2040->SetMargin(0.2);
		legendDRPHOSNLO2040->SetTextFont(42);
		legendDRPHOSNLO2040->AddEntry(graphPHOSDRPi0FitSysErr2040,"PHOS","pf");
		legendDRPHOSNLO2040->AddEntry(graphTheoryNLODR2040,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
		legendDRPHOSNLO2040->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
		legendDRPHOSNLO2040->Draw();		
		
		labelDRCent2040->Draw();
	
	//_______________________________________________________________ 40-80% panel _______________________________________________________________	
	padPartRatioInDR3->cd();
	padPartRatioInDR3->SetLogx(1);
		dummyDR3->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);	
		
		DrawGammaSetMarkerTGraphAsym(graphPHOSDRPi0FitSysErr4080, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080,widthLinesBoxes, kTRUE);
		DrawGammaSetMarker(histoPHOSDRPi0FitStatErr4080, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);
		SetStyleGammaNLOTGraphWithBand( graphTheoryNLODR4080, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
		graphTheoryNLODR4080->Draw("p3lsame");
		graphPHOSDRPi0FitSysErr4080->Draw("E2same");
		histoPHOSDRPi0FitStatErr4080->Draw("p,same,e0,X0");

		TLegend* legendDRPHOSNLO4080 = new TLegend(0.15,0.87-1.1*0.85*textsizeLabelsPad3*3,0.5,0.87);
		legendDRPHOSNLO4080->SetFillStyle(0);
		legendDRPHOSNLO4080->SetFillColor(0);
		legendDRPHOSNLO4080->SetLineColor(0);
		legendDRPHOSNLO4080->SetTextSize(0.85*textsizeLabelsPad3);
		legendDRPHOSNLO4080->SetMargin(0.2);
		legendDRPHOSNLO4080->SetTextFont(42);
		legendDRPHOSNLO4080->AddEntry(graphPHOSDRPi0FitSysErr4080,"PHOS","pf");
		legendDRPHOSNLO4080->AddEntry(graphTheoryNLODR4080,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
		legendDRPHOSNLO4080->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
		legendDRPHOSNLO4080->Draw();

		labelDRCent4080->Draw();
		
	canvasRatioIndDR->SaveAs(Form("%s/DR_PHOSMeasurementTheory.%s", outputDir.Data(), suffix.Data()));	
	

	//*******************************************************************************************************************************************
	//********************************************** DR plot with PHOS measurements all cents in one plot  ***************************************
	//*******************************************************************************************************************************************	
	canvasDoubleRatio->cd();

	dummyDR->DrawCopy();

		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);	

		graphPHOSDRPi0FitSysErr0020->Draw("E2same");
		histoPHOSDRPi0FitStatErr0020->Draw("p,same,e0,X0");
		graphPHOSDRPi0FitSysErr2040->Draw("E2same");
		histoPHOSDRPi0FitStatErr2040->Draw("p,same,e0,X0");
		graphPHOSDRPi0FitSysErr4080->Draw("E2same");
		histoPHOSDRPi0FitStatErr4080->Draw("p,same,e0,X0");

	
		TLegend* legendDRAllCentPHOS = GetAndSetLegend(0.15,0.75,4);
		legendDRAllCentPHOS->AddEntry(graphPHOSDRPi0FitSysErr0020,"  0-20% PHOS","pf");
		legendDRAllCentPHOS->AddEntry(graphPHOSDRPi0FitSysErr2040,"20-40% PHOS","pf");
		legendDRAllCentPHOS->AddEntry(graphPHOSDRPi0FitSysErr4080,"40-80% PHOS","pf");
		
		legendDRAllCentPHOS->Draw();

	canvasDoubleRatio->Print(Form("%s/DR_PHOSMeasurementAllCentsInOne.%s", outputDir.Data(), suffix.Data()));	


	
	//*******************************************************************************************************************************************
	//********************************************** DR plot with combined measurement **********************************************************
	//*******************************************************************************************************************************************	
	canvasRatioIndDR->cd();
	
	padPartRatioInDR1->Draw();
	padPartRatioInDR2->Draw();
	padPartRatioInDR3->Draw();

	TGraphAsymmErrors* graphCombDRPi0FitStatErr0020Plot = (TGraphAsymmErrors*)graphCombDRPi0FitStatErr0020->Clone("graphCombDRPi0FitStatErr0020Plot");
	ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatErr0020Plot);
	TGraphAsymmErrors* graphCombDRPi0FitStatErr2040Plot = (TGraphAsymmErrors*)graphCombDRPi0FitStatErr2040->Clone("graphCombDRPi0FitStatErr2040Plot");
	ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatErr2040Plot);
	TGraphAsymmErrors* graphCombDRPi0FitStatErr4080Plot = (TGraphAsymmErrors*)graphCombDRPi0FitStatErr4080->Clone("graphCombDRPi0FitStatErr4080Plot");
	ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatErr4080Plot);
	//_______________________________________________________________ 0-20% panel _______________________________________________________________
	padPartRatioInDR1->cd();
	padPartRatioInDR1->SetLogx(1);
		dummyDR1->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysErr0020, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);
		SetStyleGammaNLOTGraphWithBand( graphTheoryEPS09DR0020, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
		graphTheoryEPS09DR0020->Draw("p3lsame");
		SetStyleGammaNLOTGraphWithBand( graphTheoryCT10DR0020, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
		graphTheoryCT10DR0020->Draw("p3lsame");

		graphTheoryNLODR0020->Draw("p3lsame");
		
		graphCombDRPi0FitSysErr0020->Draw("E2same");
		graphCombDRPi0FitStatErr0020Plot->Draw("p,E1Z,same");		

		TLegend* legendDRCombNLO0020 = new TLegend(0.12,0.82-1.02*0.85*textsizeLabelsPad1*5,0.42,0.82);
		legendDRCombNLO0020->SetFillStyle(0);
		legendDRCombNLO0020->SetFillColor(0);
		legendDRCombNLO0020->SetLineColor(0);
		legendDRCombNLO0020->SetTextSize(0.85*textsizeLabelsPad1);
		legendDRCombNLO0020->SetMargin(0.23);
		legendDRCombNLO0020->SetTextFont(42);
		legendDRCombNLO0020->AddEntry(graphCombDRPi0FitSysErr0020,"ALICE","pf");
		legendDRCombNLO0020->AddEntry(graphTheoryNLODR0020,"NLO pQCD PDF: CTEQ6M5 FF: GRV ","l");
		legendDRCombNLO0020->AddEntry(graphTheoryCT10DR0020,"JETPHOX PDF: CT10, FF: BFG2","f");
		legendDRCombNLO0020->AddEntry(graphTheoryEPS09DR0020,"JETPHOX nPDF: EPS09, FF: BFG2","f");
		legendDRCombNLO0020->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
		legendDRCombNLO0020->Draw();
		
		labelDRCent0020->Draw();
		
	//_______________________________________________________________ 20-40% panel _______________________________________________________________
	padPartRatioInDR2->cd();
	padPartRatioInDR2->SetLogx(1);
		dummyDR2->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
		
		SetStyleGammaNLOTGraphWithBand( graphTheoryEPS09DR2040, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
		graphTheoryEPS09DR2040->Draw("p3lsame");
		SetStyleGammaNLOTGraphWithBand( graphTheoryCT10DR2040, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
		graphTheoryCT10DR2040->Draw("p3lsame");

		graphTheoryNLODR2040->Draw("p3lsame");
		
		graphCombDRPi0FitSysErr2040->Draw("E2same");
		graphCombDRPi0FitStatErr2040Plot->Draw("p,E1Z,same");

		TLegend* legendDRCombNLO2040 = new TLegend(0.12,0.85-1.02*0.85*textsizeLabelsPad2*5,0.42,0.85);
		legendDRCombNLO2040->SetFillStyle(0);
		legendDRCombNLO2040->SetFillColor(0);
		legendDRCombNLO2040->SetLineColor(0);
		legendDRCombNLO2040->SetTextSize(0.85*textsizeLabelsPad2);
		legendDRCombNLO2040->SetMargin(0.23);
		legendDRCombNLO2040->SetTextFont(42);
		legendDRCombNLO2040->AddEntry(graphCombDRPi0FitSysErr2040,"ALICE","pf");
		legendDRCombNLO2040->AddEntry(graphTheoryNLODR0020,"NLO pQCD PDF: CTEQ6M5 FF: GRV ","l");
		legendDRCombNLO2040->AddEntry(graphTheoryCT10DR0020,"JETPHOX PDF: CT10, FF: BFG2","f");
		legendDRCombNLO2040->AddEntry(graphTheoryEPS09DR0020,"JETPHOX nPDF: EPS09, FF: BFG2","f");
		legendDRCombNLO2040->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
		legendDRCombNLO2040->Draw();		
		
		labelDRCent2040->Draw();
	
	//_______________________________________________________________ 40-80% panel _______________________________________________________________	
	padPartRatioInDR3->cd();
	padPartRatioInDR3->SetLogx(1);
		dummyDR3->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);	
		
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysErr4080, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);
		SetStyleGammaNLOTGraphWithBand( graphTheoryEPS09DR4080, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
		graphTheoryEPS09DR4080->Draw("p3lsame");
		SetStyleGammaNLOTGraphWithBand( graphTheoryCT10DR4080, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
		graphTheoryCT10DR4080->Draw("p3lsame");
		graphTheoryNLODR4080->Draw("p3lsame");
		
		graphCombDRPi0FitSysErr4080->Draw("E2same");
		graphCombDRPi0FitStatErr4080Plot->Draw("p,E1Z,same");

		TLegend* legendDRCombNLO4080 = new TLegend(0.12,0.87-1.02*0.85*textsizeLabelsPad3*5,0.42,0.87);
		legendDRCombNLO4080->SetFillStyle(0);
		legendDRCombNLO4080->SetFillColor(0);
		legendDRCombNLO4080->SetLineColor(0);
		legendDRCombNLO4080->SetTextSize(0.85*textsizeLabelsPad3);
		legendDRCombNLO4080->SetMargin(0.23);
		legendDRCombNLO4080->SetTextFont(42);
		legendDRCombNLO4080->AddEntry(graphCombDRPi0FitSysErr4080,"ALICE","pf");
		legendDRCombNLO4080->AddEntry(graphTheoryNLODR0020,"NLO pQCD PDF: CTEQ6M5 FF: GRV ","l");
		legendDRCombNLO4080->AddEntry(graphTheoryCT10DR0020,"JETPHOX PDF: CT10, FF: BFG2","f");
		legendDRCombNLO4080->AddEntry(graphTheoryEPS09DR0020,"JETPHOX nPDF: EPS09, FF: BFG2","f");
		legendDRCombNLO4080->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
		legendDRCombNLO4080->Draw();
				
		labelDRCent4080->Draw();
		
	canvasRatioIndDR->SaveAs(Form("%s/DR_combMeasurement_incNLO.%s", outputDir.Data(), suffix.Data()));	

	
	//*******************************************************************************************************************************************
	//********************************************** DR plot with combined measurement inc NLO **************************************************
	//*******************************************************************************************************************************************	
	
	canvasRatioIndDR->cd();
	
	padPartRatioInDR1->Draw();
	padPartRatioInDR2->Draw();
	padPartRatioInDR3->Draw();

	
	//_______________________________________________________________ 0-20% panel _______________________________________________________________
	padPartRatioInDR1->cd();
	padPartRatioInDR1->SetLogx(1);
		dummyDR1->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		graphCombDRPi0FitSysErr0020->Draw("E2same");
		graphCombDRPi0FitStatErr0020Plot->Draw("p,E1Z,same");		
		
		labelDRCent0020->Draw();

		TLatex *labelALICECent0020 = new TLatex(0.87,0.85,"ALICE");
		SetStyleTLatex( labelALICECent0020, 0.85*textsizeLabelsPad1,4);
		labelALICECent0020->Draw();

		
	//_______________________________________________________________ 20-40% panel _______________________________________________________________
	padPartRatioInDR2->cd();
	padPartRatioInDR2->SetLogx(1);
		dummyDR2->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

		graphCombDRPi0FitSysErr2040->Draw("E2same");
		graphCombDRPi0FitStatErr2040Plot->Draw("p,E1Z,same");
				
		labelDRCent2040->Draw();
	
		TLatex *labelALICECent2040 = new TLatex(0.87,0.88,"ALICE");
		SetStyleTLatex( labelALICECent2040, 0.85*textsizeLabelsPad2,4);
		labelALICECent2040->Draw();

		
	//_______________________________________________________________ 40-80% panel _______________________________________________________________	
	padPartRatioInDR3->cd();
	padPartRatioInDR3->SetLogx(1);
		dummyDR3->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);	
		
		graphCombDRPi0FitSysErr4080->Draw("E2same");
		graphCombDRPi0FitStatErr4080Plot->Draw("p,E1Z,same");

		labelDRCent4080->Draw();
		
		TLatex *labelALICECent4080 = new TLatex(0.87,0.9,"ALICE");
		SetStyleTLatex( labelALICECent4080, 0.85*textsizeLabelsPad3,4);
		labelALICECent4080->Draw();

		
	canvasRatioIndDR->SaveAs(Form("%s/DR_combMeasurement.%s", outputDir.Data(), suffix.Data()));	

	//*******************************************************************************************************************************************
	//********************************************** DR plot with combined measurement diff Syst Err ********************************************
	//*******************************************************************************************************************************************	
	
	ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatSysAErr0020);
	ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatSysAErr2040);
	ProduceGraphAsymmWithoutXErrors(graphCombDRPi0FitStatSysAErr4080);
	
	canvasRatioIndDR->cd();
	
	padPartRatioInDR1->Draw();
	padPartRatioInDR2->Draw();
	padPartRatioInDR3->Draw();

	
	//_______________________________________________________________ 0-20% panel _______________________________________________________________
	padPartRatioInDR1->cd();
	padPartRatioInDR1->SetLogx(1);
		dummyDR1->Draw("");
		
		TBox* boxSysCErrComb0020 = CreateBoxConv(colorComb0020Box, 0.74, 1.-(SysCCombDRPi0Fit0020) , 0.83, 1.+(SysCCombDRPi0Fit0020));
		boxSysCErrComb0020->Draw();
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
	
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysCErr0020, markerStyleComb0020, markerSizeComb0020, colorComb0020Box , colorComb0020Box, widthLinesBoxes,  kTRUE, colorComb0020Box);		
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysBErr0020, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatSysAErr0020, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);		
		graphCombDRPi0FitSysBErr0020->Draw("E2same");
		graphCombDRPi0FitStatSysAErr0020->Draw("p,E1Z,same");		
	
		TLegend* legendDRCombDiffErrRepPad1 = new TLegend(0.15,0.82-1.1*0.85*textsizeLabelsPad1*3,0.35,0.82);
		legendDRCombDiffErrRepPad1->SetFillStyle(0);
		legendDRCombDiffErrRepPad1->SetFillColor(0);
		legendDRCombDiffErrRepPad1->SetLineColor(0);
		legendDRCombDiffErrRepPad1->SetTextSize(0.85*textsizeLabelsPad1);
		legendDRCombDiffErrRepPad1->SetMargin(0.25);
		legendDRCombDiffErrRepPad1->SetTextFont(42);
		legendDRCombDiffErrRepPad1->AddEntry(graphCombDRPi0FitStatSysAErr0020,"Stat. #oplus Syst. A","pe");
		legendDRCombDiffErrRepPad1->AddEntry(graphCombDRPi0FitSysBErr0020,"Syst. B","f");
		legendDRCombDiffErrRepPad1->AddEntry(graphCombDRPi0FitSysCErr0020,"Syst. C","f");
		legendDRCombDiffErrRepPad1->Draw();
				
		labelDRCent0020->Draw();
		labelALICECent0020->Draw();
	//_______________________________________________________________ 20-40% panel _______________________________________________________________
	padPartRatioInDR2->cd();
	padPartRatioInDR2->SetLogx(1);
		dummyDR2->Draw("");
		
		TBox* boxSysCErrComb2040 = CreateBoxConv(colorComb2040Box, 0.74, 1.-(SysCCombDRPi0Fit2040) , 0.83, 1.+(SysCCombDRPi0Fit2040));
		boxSysCErrComb2040->Draw();
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysCErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040Box , colorComb2040Box, widthLinesBoxes,  kTRUE, colorComb2040Box);		
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysBErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatSysAErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);		
		graphCombDRPi0FitSysBErr2040->Draw("E2same");
		graphCombDRPi0FitStatSysAErr2040->Draw("p,E1Z,same");		
	
		TLegend* legendDRCombDiffErrRepPad2 = new TLegend(0.15,0.86-1.1*0.85*textsizeLabelsPad2*3,0.35,0.86);
		legendDRCombDiffErrRepPad2->SetFillStyle(0);
		legendDRCombDiffErrRepPad2->SetFillColor(0);
		legendDRCombDiffErrRepPad2->SetLineColor(0);
		legendDRCombDiffErrRepPad2->SetTextSize(0.85*textsizeLabelsPad2);
		legendDRCombDiffErrRepPad2->SetMargin(0.25);
		legendDRCombDiffErrRepPad2->SetTextFont(42);
		legendDRCombDiffErrRepPad2->AddEntry(graphCombDRPi0FitStatSysAErr2040,"Stat. #oplus Syst. A","pe");
		legendDRCombDiffErrRepPad2->AddEntry(graphCombDRPi0FitSysBErr2040,"Syst. B","f");
		legendDRCombDiffErrRepPad2->AddEntry(graphCombDRPi0FitSysCErr2040,"Syst. C","f");
		legendDRCombDiffErrRepPad2->Draw();

		labelDRCent2040->Draw();
		labelALICECent2040->Draw();
	//_______________________________________________________________ 40-80% panel _______________________________________________________________	
	padPartRatioInDR3->cd();
	padPartRatioInDR3->SetLogx(1);
		dummyDR3->Draw("");
		TBox* boxSysCErrComb4080 = CreateBoxConv(colorComb4080Box, 0.74, 1.-(SysCCombDRPi0Fit4080) , 0.83, 1.+(SysCCombDRPi0Fit4080));
		boxSysCErrComb4080->Draw();
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysCErr4080, markerStyleComb4080, markerSizeComb4080, colorComb4080Box , colorComb4080Box, widthLinesBoxes,  kTRUE, colorComb4080Box);		
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitSysBErr4080, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphCombDRPi0FitStatSysAErr4080, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);		
		graphCombDRPi0FitSysBErr4080->Draw("E2same");
		graphCombDRPi0FitStatSysAErr4080->Draw("p,E1Z,same");		

		TLegend* legendDRCombDiffErrRepPad3 = new TLegend(0.15,0.88-1.1*0.85*textsizeLabelsPad3*3,0.35,0.88);
		legendDRCombDiffErrRepPad3->SetFillStyle(0);
		legendDRCombDiffErrRepPad3->SetFillColor(0);
		legendDRCombDiffErrRepPad3->SetLineColor(0);
		legendDRCombDiffErrRepPad3->SetTextSize(0.85*textsizeLabelsPad3);
		legendDRCombDiffErrRepPad3->SetMargin(0.25);
		legendDRCombDiffErrRepPad3->SetTextFont(42);
		legendDRCombDiffErrRepPad3->AddEntry(graphCombDRPi0FitStatSysAErr4080,"Stat. #oplus Syst. A","pe");
		legendDRCombDiffErrRepPad3->AddEntry(graphCombDRPi0FitSysBErr4080,"Syst. B","f");
		legendDRCombDiffErrRepPad3->AddEntry(graphCombDRPi0FitSysCErr4080,"Syst. C","f");
		legendDRCombDiffErrRepPad3->Draw();

		labelDRCent4080->Draw();
		labelALICECent4080->Draw();
		
	canvasRatioIndDR->SaveAs(Form("%s/DR_combMeasurement_diffErrRep.%s", outputDir.Data(), suffix.Data()));	
	
	
	//*******************************************************************************************************************************************
	//************************************ plot ratio of individual measurements to combined fit for inc gamma **********************************
	//*******************************************************************************************************************************************
	canvasRatioIndDR->cd();
	
	padPartRatioInDR1->Draw();
	padPartRatioInDR2->Draw();
	padPartRatioInDR3->Draw();

	//_______________________________________________________________ 0-20% dummy upper panel ___________________________________________________
	TH2D *dummyIncGammaIndMeas1 ;
	dummyIncGammaIndMeas1 = new TH2D("dummyIncGammaIndMeas1", "dummyIncGammaIndMeas1", 1000, 0., 22, 1000., indMeasRatio[0], indMeasRatio[1]);
	SetStyleHistoTH2ForGraphs( dummyIncGammaIndMeas1, "#it{p}_{T} (GeV/#it{c})", "", 
							   0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.85*textsizeLabelsPad1, textsizeLabelsPad1, 0.95,0.105/(textsizeFacPad1*margin), 510, 505);
	dummyIncGammaIndMeas1->GetXaxis()->SetLabelOffset(-0.015);
	dummyIncGammaIndMeas1->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyIncGammaIndMeas1->GetXaxis()->SetTickLength(0.06);
	dummyIncGammaIndMeas1->GetYaxis()->SetTickLength(0.028);

	//_______________________________________________________________ 20-40% dummy middle panel _________________________________________________
	TH2D *dummyIncGammaIndMeas2 ;
	dummyIncGammaIndMeas2 = new TH2D("dummyIncGammaIndMeas2", "dummyIncGammaIndMeas2", 1000, 0., 22, 1000., indMeasRatio[0], indMeasRatio[1]);
	SetStyleHistoTH2ForGraphs( dummyIncGammaIndMeas2, "#it{p}_{T} (GeV/#it{c})", "#gamma_{inc} (Data / Combined)", 
							   0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.85*textsizeLabelsPad2, textsizeLabelsPad2, 0.95,0.105/(textsizeFacPad2*margin), 510, 505);
	dummyIncGammaIndMeas2->GetXaxis()->SetLabelOffset(-0.015);
	dummyIncGammaIndMeas2->GetYaxis()->CenterTitle(kTRUE);
	dummyIncGammaIndMeas2->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyIncGammaIndMeas2->GetXaxis()->SetTickLength(0.06);
	dummyIncGammaIndMeas2->GetYaxis()->SetTickLength(0.028);

	//_______________________________________________________________ 40-80% dummy lower panel __________________________________________________
	TH2D *dummyIncGammaIndMeas3 ;
	dummyIncGammaIndMeas3 = new TH2D("dummyIncGammaIndMeas3", "dummyIncGammaIndMeas3", 1000, 0., 22, 1000., indMeasRatio[0], indMeasRatio[1]);
	SetStyleHistoTH2ForGraphs( dummyIncGammaIndMeas3, "#it{p}_{T} (GeV/#it{c})", "", 
							   0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.85*textsizeLabelsPad3, textsizeLabelsPad3, 0.95,0.105/(textsizeFacPad3*margin), 510, 505);
	dummyIncGammaIndMeas3->GetXaxis()->SetLabelOffset(-0.015);
	dummyIncGammaIndMeas3->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyIncGammaIndMeas3->GetXaxis()->SetTickLength(0.055);
	dummyIncGammaIndMeas3->GetYaxis()->SetTickLength(0.035);
	
	
	//_______________________________________________________________ 0-20% panel _______________________________________________________________
	padPartRatioInDR1->cd();
	padPartRatioInDR1->SetLogx(1);
		dummyIncGammaIndMeas1->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysErr0020, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
		ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatErr0020);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatErr0020, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		graphRatioCombPPHOSIncGammaSysErr0020->Draw("E2same");
		graphRatioCombPPHOSIncGammaStatErr0020->Draw("p,E1Z,same");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysErr0020, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
		ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatErr0020);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatErr0020, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		graphRatioCombPPCMIncGammaSysErr0020->Draw("E2same");
		graphRatioCombPPCMIncGammaStatErr0020->Draw("p,E1Z,same");

		TLatex *labelALICEIncGamma = new TLatex(0.12,0.76,"ALICE");
		SetStyleTLatex( labelALICEIncGamma, 0.85*textsizeLabelsPad1,4);
		labelALICEIncGamma->Draw();
		
		TLegend* legendIncGammaRIndMeas = new TLegend(0.8,0.92-1.1*0.85*textsizeLabelsPad1*2,0.95,0.92);
		legendIncGammaRIndMeas->SetFillStyle(0);
		legendIncGammaRIndMeas->SetFillColor(0);
		legendIncGammaRIndMeas->SetLineColor(0);
		legendIncGammaRIndMeas->SetTextSize(0.85*textsizeLabelsPad1);
		legendIncGammaRIndMeas->SetMargin(0.4);
		legendIncGammaRIndMeas->SetTextFont(42);
		legendIncGammaRIndMeas->AddEntry(graphRatioCombPPCMIncGammaSysErr0020,"PCM","pf");
		legendIncGammaRIndMeas->AddEntry(graphRatioCombPPHOSIncGammaSysErr0020,"PHOS","pf");
		legendIncGammaRIndMeas->Draw();
		
		labelDRCent0020->Draw();
		
	//_______________________________________________________________ 20-40% panel _______________________________________________________________
	padPartRatioInDR2->cd();
	padPartRatioInDR2->SetLogx(1);
		dummyIncGammaIndMeas2->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);

		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
		ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatErr2040);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		graphRatioCombPPHOSIncGammaSysErr2040->Draw("E2same");
		graphRatioCombPPHOSIncGammaStatErr2040->Draw("p,E1Z,same");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysErr2040, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
		ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatErr2040);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatErr2040, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		graphRatioCombPPCMIncGammaSysErr2040->Draw("E2same");
		graphRatioCombPPCMIncGammaStatErr2040->Draw("p,E1Z,same");
				
		labelDRCent2040->Draw();
	
	//_______________________________________________________________ 40-80% panel _______________________________________________________________	
	padPartRatioInDR3->cd();
	padPartRatioInDR3->SetLogx(1);
		dummyIncGammaIndMeas3->Draw("");
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);	
		
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysErr4080, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
		ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatErr4080);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatErr4080, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		graphRatioCombPPHOSIncGammaSysErr4080->Draw("E2same");
		graphRatioCombPPHOSIncGammaStatErr4080->Draw("p,E1Z,same");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysErr4080, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
		ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatErr4080);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatErr4080, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		graphRatioCombPPCMIncGammaSysErr4080->Draw("E2same");
		graphRatioCombPPCMIncGammaStatErr4080->Draw("p,E1Z,same");

		labelDRCent4080->Draw();
		
	canvasRatioIndDR->SaveAs(Form("%s/Ratio_indMeasurement_incGamma.%s", outputDir.Data(), suffix.Data()));	

	//*******************************************************************************************************************************************
	//************************************ plot ratio of individual measurements to combined fit for inc gamma **********************************
	//*******************************************************************************************************************************************
	canvasRatioIndDR->cd();
	
	padPartRatioInDR1->Draw();
	padPartRatioInDR2->Draw();
	padPartRatioInDR3->Draw();

	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr0020Plot = (TGraphAsymmErrors*)graphRatioCombPPHOSIncGammaStatSysAErr0020->Clone("graphRatioCombPPHOSIncGammaStatSysAErr0020Plot");
	ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatSysAErr0020Plot);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr2040Plot = (TGraphAsymmErrors*)graphRatioCombPPHOSIncGammaStatSysAErr2040->Clone("graphRatioCombPPHOSIncGammaStatSysAErr2040Plot");
	ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatSysAErr2040Plot);
	TGraphAsymmErrors* graphRatioCombPPHOSIncGammaStatSysAErr4080Plot = (TGraphAsymmErrors*)graphRatioCombPPHOSIncGammaStatSysAErr4080->Clone("graphRatioCombPPHOSIncGammaStatSysAErr4080Plot");
	ProduceGraphAsymmWithoutXErrors(graphRatioCombPPHOSIncGammaStatSysAErr4080Plot);

	TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr0020Plot = (TGraphAsymmErrors*)graphRatioCombPPCMIncGammaStatSysAErr0020->Clone("graphRatioCombPPCMIncGammaStatSysAErr0020Plot");
	ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatSysAErr0020Plot);
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr2040Plot = (TGraphAsymmErrors*)graphRatioCombPPCMIncGammaStatSysAErr2040->Clone("graphRatioCombPPCMIncGammaStatSysAErr2040Plot");
	ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatSysAErr2040Plot);
	TGraphAsymmErrors* graphRatioCombPPCMIncGammaStatSysAErr4080Plot = (TGraphAsymmErrors*)graphRatioCombPPCMIncGammaStatSysAErr4080->Clone("graphRatioCombPPCMIncGammaStatSysAErr4080Plot");
	ProduceGraphAsymmWithoutXErrors(graphRatioCombPPCMIncGammaStatSysAErr4080Plot);
	
	//_______________________________________________________________ 0-20% panel _______________________________________________________________
	padPartRatioInDR1->cd();
	padPartRatioInDR1->SetLogx(1);
		dummyIncGammaIndMeas1->Draw("");
				
		TBox* boxSysCErrPCM0020 = CreateBoxConv(colorPCMBox, 0.76, 1.-(SysCPCMIncGamma0020) , 0.815, 1.+(SysCPCMIncGamma0020));
		boxSysCErrPCM0020->Draw();
		TBox* boxSysCErrPHOS0020 = CreateBoxConv(colorPHOSBox, 0.815, 1.-(SysCPHOSIncGamma0020) , 0.88, 1.+(SysCPHOSIncGamma0020));
		boxSysCErrPHOS0020->Draw();

		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysBErr0020, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatSysAErr0020Plot, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		graphRatioCombPPHOSIncGammaSysBErr0020->Draw("E2same");
		graphRatioCombPPHOSIncGammaStatSysAErr0020Plot->Draw("p,E1Z,same");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysCErr0020, markerStylePCM, markerSizePCM, colorPCMBox , colorPCMBox,widthLinesBoxes, kTRUE, colorPCMBox);		
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysBErr0020, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatSysAErr0020Plot, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		graphRatioCombPPCMIncGammaSysBErr0020->Draw("E2same");
		graphRatioCombPPCMIncGammaStatSysAErr0020Plot->Draw("p,E1Z,same");

		legendIncGammaRIndMeas->Draw();
		dummyIncGammaIndMeas1->Draw("axis,same");
		labelDRCent0020->Draw();
		labelALICEIncGamma->Draw();
		
	//_______________________________________________________________ 20-40% panel _______________________________________________________________
	padPartRatioInDR2->cd();
	padPartRatioInDR2->SetLogx(1);
		dummyIncGammaIndMeas2->Draw("");
		
		TLegend* legendIncGammaRIndMeas2 = new TLegend(0.65,0.95-1.1*0.85*textsizeLabelsPad2*3,0.9,0.95);
		legendIncGammaRIndMeas2->SetFillStyle(0);
		legendIncGammaRIndMeas2->SetFillColor(0);
		legendIncGammaRIndMeas2->SetLineColor(0);
		legendIncGammaRIndMeas2->SetTextSize(0.85*textsizeLabelsPad2);
		legendIncGammaRIndMeas2->SetMargin(0.45);
		legendIncGammaRIndMeas2->SetTextFont(42);
		legendIncGammaRIndMeas2->SetNColumns(2);
		legendIncGammaRIndMeas2->SetHeader("For both methods:");
		legendIncGammaRIndMeas2->AddEntry(graphRatioCombPPCMIncGammaStatSysAErr0020Plot,"Stat. #oplus Syst. A","pel");
		legendIncGammaRIndMeas2->AddEntry((TObject*)0,"","");
		legendIncGammaRIndMeas2->AddEntry(graphRatioCombPPCMIncGammaSysBErr0020,"Syst. B","f");
		legendIncGammaRIndMeas2->AddEntry(graphRatioCombPPCMIncGammaSysCErr0020,"Syst. C","f");
		legendIncGammaRIndMeas2->Draw();

		
		TBox* boxSysCErrPCM2040 = CreateBoxConv(colorPCMBox, 0.76, 1.-(SysCPCMIncGamma2040) , 0.815, 1.+(SysCPCMIncGamma2040));
		boxSysCErrPCM2040->Draw();
		TBox* boxSysCErrPHOS2040 = CreateBoxConv(colorPHOSBox, 0.815, 1.-(SysCPHOSIncGamma2040) , 0.88, 1.+(SysCPHOSIncGamma2040));
		boxSysCErrPHOS2040->Draw();
		
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
				
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysBErr2040, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatSysAErr2040Plot, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		graphRatioCombPPHOSIncGammaSysBErr2040->Draw("E2same");
		graphRatioCombPPHOSIncGammaStatSysAErr2040Plot->Draw("p,E1Z,same");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysBErr2040, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatSysAErr2040Plot, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		graphRatioCombPPCMIncGammaSysBErr2040->Draw("E2same");
		graphRatioCombPPCMIncGammaStatSysAErr2040Plot->Draw("p,E1Z,same");
		
		dummyIncGammaIndMeas2->Draw("axis,same");
		labelDRCent2040->Draw();
	
	//_______________________________________________________________ 40-80% panel _______________________________________________________________	
	padPartRatioInDR3->cd();
	padPartRatioInDR3->SetLogx(1);
		dummyIncGammaIndMeas3->Draw("");
		TBox* boxSysCErrPCM4080 = CreateBoxConv(colorPCMBox, 0.76, 1.-(SysCPCMIncGamma4080) , 0.815, 1.+(SysCPCMIncGamma4080));
		boxSysCErrPCM4080->Draw();
		TBox* boxSysCErrPHOS4080 = CreateBoxConv(colorPHOSBox, 0.815, 1.-(SysCPHOSIncGamma4080) , 0.88, 1.+(SysCPHOSIncGamma4080));
		boxSysCErrPHOS4080->Draw();
		
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaSysBErr4080, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPHOSIncGammaStatSysAErr4080Plot, markerStylePHOS, markerSizePHOS, colorPHOS , colorPHOS);
		graphRatioCombPPHOSIncGammaSysBErr4080->Draw("E2same");
		graphRatioCombPPHOSIncGammaStatSysAErr4080Plot->Draw("p,E1Z,same");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaSysBErr4080, markerStylePCM, markerSizePCM, colorPCM , colorPCM,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombPPCMIncGammaStatSysAErr4080Plot, markerStylePCM, markerSizePCM, colorPCM , colorPCM);
		graphRatioCombPPCMIncGammaSysBErr4080->Draw("E2same");
		graphRatioCombPPCMIncGammaStatSysAErr4080Plot->Draw("p,E1Z,same");

		dummyIncGammaIndMeas3->Draw("axis,same");

		labelDRCent4080->Draw();
		
	canvasRatioIndDR->SaveAs(Form("%s/Ratio_indMeasurementDiffErrorRep_incGamma.%s", outputDir.Data(), suffix.Data()));	


	//*******************************************************************************************************************************************
	//************************************ plot ratio of dirGamma measurement to combined fit for  **********************************************
	//*******************************************************************************************************************************************
	canvasRatioIndDR->cd();
	
	padPartRatioInDR1->Draw();
	padPartRatioInDR2->Draw();
	padPartRatioInDR3->Draw();
	
	//_______________________________________________________________ 0-20% panel _______________________________________________________________
	padPartRatioInDR1->cd();
	padPartRatioInDR1->SetLogx(1);
		dummyIncGammaIndMeas1->Draw("");
				
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		if (graphRatioCombFitDirGammaSysErr0020){
			DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaSysErr0020, markerStyleComb0020, markerSizeComb2040, colorComb0020 , colorComb0020,widthLinesBoxes, kTRUE);
			graphRatioCombFitDirGammaSysErr0020->Draw("E2same");
		}	
		if (graphRatioCombFitDirGammaStatErr0020){
			DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaStatErr0020, markerStyleComb0020, markerSizeComb2040, colorComb0020 , colorComb0020);
			graphRatioCombFitDirGammaStatErr0020->Draw("p,E1Z,same");
		}
// 		legendIncGammaRIndMeas->Draw();
		dummyIncGammaIndMeas1->Draw("axis,same");
		labelDRCent0020->Draw();
		labelALICEIncGamma->Draw();
		
	//_______________________________________________________________ 20-40% panel _______________________________________________________________
	padPartRatioInDR2->cd();
	padPartRatioInDR2->SetLogx(1);
		dummyIncGammaIndMeas2->Draw("");
		
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
					
		DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaSysErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaStatErr2040, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);
		graphRatioCombFitDirGammaSysErr2040->Draw("E2same");
		graphRatioCombFitDirGammaStatErr2040->Draw("p,E1Z,same");
		
		dummyIncGammaIndMeas2->Draw("axis,same");
		labelDRCent2040->Draw();
	
	//_______________________________________________________________ 40-80% panel _______________________________________________________________	
	padPartRatioInDR3->cd();
	padPartRatioInDR3->SetLogx(1);
		dummyIncGammaIndMeas3->Draw("");
		
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaSysErr4080, markerStyleComb4080, markerSizeComb2040, colorComb4080 , colorComb4080,widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombFitDirGammaStatErr4080, markerStyleComb4080, markerSizeComb2040, colorComb4080 , colorComb4080);
		graphRatioCombFitDirGammaSysErr4080->Draw("E2same");
		graphRatioCombFitDirGammaStatErr4080->Draw("p,E1Z,same");

		dummyIncGammaIndMeas3->Draw("axis,same");

		labelDRCent4080->Draw();
		
	canvasRatioIndDR->SaveAs(Form("%s/Ratio_DirGammaCombMeasurementToFit.%s", outputDir.Data(), suffix.Data()));	

	
	
	//*******************************************************************************************************************************************
	//********************************************** DR plot with PCM pp & pPb measurement + NLO  ***********************************************
	//*******************************************************************************************************************************************	
	canvasRatioIndDR->cd();
	
	padPartRatioInDR1->Draw();
	padPartRatioInDR2->Draw();
	padPartRatioInDR3->Draw();

	
	//_______________________________________________________________ pPb panel _______________________________________________________________
	padPartRatioInDR1->cd();
	padPartRatioInDR1->SetLogx(1);
		dummyDR1->GetXaxis()->SetRangeUser(doubleRatioXpp[0],doubleRatioXpp[1]);
		dummyDR1->Draw("");
		DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);
		
		DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErrpPb, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb,widthLinesBoxes, kTRUE);
		DrawGammaSetMarker(histoPCMDRPi0FitStatErrpPb, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb);
		SetStyleGammaNLOTGraphWithBand( graphTheoryNLODRpPb, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
		graphTheoryNLODRpPb->Draw("pE3lsame");
		graphPCMDRPi0FitSysErrpPb->Draw("E2same");
		histoPCMDRPi0FitStatErrpPb->Draw("p,same,e0,X0");

		
		
		TLegend* legendDRPCMNLOpPb = new TLegend(0.15,0.82-1.1*0.85*textsizeLabelsPad1*3,0.5,0.82);
		legendDRPCMNLOpPb->SetFillStyle(0);
		legendDRPCMNLOpPb->SetFillColor(0);
		legendDRPCMNLOpPb->SetLineColor(0);
		legendDRPCMNLOpPb->SetTextSize(0.85*textsizeLabelsPad1);
		legendDRPCMNLOpPb->SetMargin(0.2);
		legendDRPCMNLOpPb->SetTextFont(42);
		legendDRPCMNLOpPb->AddEntry(graphPCMDRPi0FitSysErrpPb,"PCM","pf");
		legendDRPCMNLOpPb->AddEntry(graphTheoryNLODRpPb,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
		legendDRPCMNLOpPb->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
		legendDRPCMNLOpPb->Draw();

		TLatex *labelDRCentpPb = new TLatex(0.12,0.85,collisionSystempPb.Data());
		SetStyleTLatex( labelDRCentpPb, 0.85*textsizeLabelsPad1,4);
		labelDRCentpPb->Draw();
		
	//_______________________________________________________________ pp 7TeV panel _______________________________________________________________
	padPartRatioInDR2->cd();
	padPartRatioInDR2->SetLogx(1);
		dummyDR2->GetXaxis()->SetRangeUser(doubleRatioXpp[0],doubleRatioXpp[1]);
		dummyDR2->Draw("");
		DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);

		DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErrpp7TeV, markerStyleCombpp7TeV, markerSizeCombpp7TeV, colorCombpp7TeV , colorCombpp7TeV,widthLinesBoxes, kTRUE);
		DrawGammaSetMarker(histoPCMDRPi0FitStatErrpp7TeV, markerStyleCombpp7TeV, markerSizeCombpp7TeV, colorCombpp7TeV , colorCombpp7TeV);
		SetStyleGammaNLOTGraphWithBand( graphTheoryNLODRpp7TeV, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
		graphTheoryNLODRpp7TeV->Draw("p3lsame");
		graphPCMDRPi0FitSysErrpp7TeV->Draw("E2same");
		histoPCMDRPi0FitStatErrpp7TeV->Draw("p,same,e0,X0");

		TLegend* legendDRPCMNLOpp7TeV = new TLegend(0.15,0.85-1.1*0.85*textsizeLabelsPad2*3,0.5,0.85);
		legendDRPCMNLOpp7TeV->SetFillStyle(0);
		legendDRPCMNLOpp7TeV->SetFillColor(0);
		legendDRPCMNLOpp7TeV->SetLineColor(0);
		legendDRPCMNLOpp7TeV->SetTextSize(0.85*textsizeLabelsPad2);
		legendDRPCMNLOpp7TeV->SetMargin(0.2);
		legendDRPCMNLOpp7TeV->SetTextFont(42);
		legendDRPCMNLOpp7TeV->AddEntry(graphPCMDRPi0FitSysErrpp7TeV,"PCM","pf");
		legendDRPCMNLOpp7TeV->AddEntry(graphTheoryNLODRpp7TeV,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
		legendDRPCMNLOpp7TeV->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
		legendDRPCMNLOpp7TeV->Draw();		

		TLatex *labelDRCentpp7TeV = new TLatex(0.12,0.88,collisionSystempp7TeV.Data());
		SetStyleTLatex( labelDRCentpp7TeV, 0.85*textsizeLabelsPad2,4);		
		labelDRCentpp7TeV->Draw();
	
	//_______________________________________________________________ pp 2.76TeV panel _______________________________________________________________	
	padPartRatioInDR3->cd();
	padPartRatioInDR3->SetLogx(1);
		dummyDR3->GetXaxis()->SetRangeUser(doubleRatioXpp[0],doubleRatioXpp[1]);
		dummyDR3->Draw("");
		DrawGammaLines(doubleRatioXpp[0], doubleRatioXpp[1], 1., 1., 1.2, kGray+2, 7);	
		
		DrawGammaSetMarkerTGraphAsym(graphPCMDRPi0FitSysErrpp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV,widthLinesBoxes, kTRUE);
		DrawGammaSetMarker(histoPCMDRPi0FitStatErrpp2760GeV, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);
		SetStyleGammaNLOTGraphWithBand( graphTheoryNLODRpp2760GeV, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
		graphTheoryNLODRpp2760GeV->Draw("p3lsame");
		graphPCMDRPi0FitSysErrpp2760GeV->Draw("E2same");
		histoPCMDRPi0FitStatErrpp2760GeV->Draw("p,same,e0,X0");

		TLegend* legendDRPCMNLOpp2760GeV = new TLegend(0.15,0.87-1.1*0.85*textsizeLabelsPad3*3,0.5,0.87);
		legendDRPCMNLOpp2760GeV->SetFillStyle(0);
		legendDRPCMNLOpp2760GeV->SetFillColor(0);
		legendDRPCMNLOpp2760GeV->SetLineColor(0);
		legendDRPCMNLOpp2760GeV->SetTextSize(0.85*textsizeLabelsPad3);
		legendDRPCMNLOpp2760GeV->SetMargin(0.2);
		legendDRPCMNLOpp2760GeV->SetTextFont(42);
		legendDRPCMNLOpp2760GeV->AddEntry(graphPCMDRPi0FitSysErrpp2760GeV,"PCM","pf");
		legendDRPCMNLOpp2760GeV->AddEntry(graphTheoryNLODRpp2760GeV,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
		legendDRPCMNLOpp2760GeV->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
		legendDRPCMNLOpp2760GeV->Draw();

		TLatex *labelDRCentpp2760GeV = new TLatex(0.12,0.9,collisionSystempp2760GeV.Data());
		SetStyleTLatex( labelDRCentpp2760GeV, 0.85*textsizeLabelsPad3,4);
		labelDRCentpp2760GeV->Draw();
		
	canvasRatioIndDR->SaveAs(Form("%s/DR_PCMMeasurementTheory_ppAndpPb.%s", outputDir.Data(), suffix.Data()));	

	//*******************************************************************************************************************************************
	//*************************************************** Plotting inclusive Gamma Spectrum *****************************************************
	//*******************************************************************************************************************************************
	
	TGraphAsymmErrors* graphCombIncGammaSysErr0020Plot = ScaleGraph(graphCombIncGammaSysErr0020,100);
	TGraphAsymmErrors* graphCombIncGammaStatErr0020Plot = ScaleGraph(graphCombIncGammaStatErr0020,100);
	histoFitQCDIncGammaComb0020->Scale(100);
	
	TGraphAsymmErrors* graphCombIncGammaSysErr2040Plot = ScaleGraph(graphCombIncGammaSysErr2040,10);
	TGraphAsymmErrors* graphCombIncGammaStatErr2040Plot = ScaleGraph(graphCombIncGammaStatErr2040,10);
	histoFitQCDIncGammaComb2040->Scale(10);	
	
	TGraphAsymmErrors* graphCombIncGammaSysErr4080Plot = ScaleGraph(graphCombIncGammaSysErr4080,1);
	TGraphAsymmErrors* graphCombIncGammaStatErr4080Plot = ScaleGraph(graphCombIncGammaStatErr4080,1);
	histoFitQCDIncGammaComb4080->Scale(1);
	TH1D* histoFitDummyPlotting = (TH1D*) histoFitQCDIncGammaComb4080->Clone("histoFitDummyPlotting");
	
	
	TCanvas *canvasIncGamma = new TCanvas("canvasIncGamma","",10,10,1200,1400);  // gives the page size		
	DrawGammaCanvasSettings( canvasIncGamma, 0.16, 0.01, 0.01, 0.07);
	canvasIncGamma->SetLogy();
	canvasIncGamma->SetLogx();

	Int_t textSizeLabelsPixelIncGam = 48;
	Double_t textsizeLabelsIncGamma = 0;
	if (canvasIncGamma->XtoPixel(canvasIncGamma->GetX2()) < canvasIncGamma->YtoPixel(canvasIncGamma->GetY1())){
		textsizeLabelsIncGamma = (Double_t)textSizeLabelsPixelIncGam/canvasIncGamma->XtoPixel(canvasIncGamma->GetX2()) ;
	} else {
		textsizeLabelsIncGamma = (Double_t)textSizeLabelsPixelIncGam/canvasIncGamma->YtoPixel(canvasIncGamma->GetY1());
	}
		
	TH2D *dummyGamma ;
	dummyGamma = new TH2D("dummyGamma", "dummyGamma", 1000, 0., 22, 1000., 6e-8,8e3);
	SetStyleHistoTH2ForGraphs( dummyGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{inc}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.85*textsizeLabelsIncGamma, textsizeLabelsIncGamma, 0.85*textsizeLabelsIncGamma, textsizeLabelsIncGamma, 0.75, 1.6);
	dummyGamma->GetXaxis()->SetLabelOffset(-0.015);
	dummyGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
// 	dummyGamma->GetXaxis()->SetRangeUser(0,16);
	dummyGamma->DrawCopy();

	DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSysErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020,widthLinesBoxes, kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);	
	graphCombIncGammaSysErr0020Plot->Draw("E2same");
	graphCombIncGammaStatErr0020Plot->Draw("p,E1Z,same");

	DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSysErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
	graphCombIncGammaSysErr2040Plot->Draw("E2same");
	graphCombIncGammaStatErr2040Plot->Draw("p,E1Z,same");
	
	DrawGammaSetMarkerTGraphAsym(graphCombIncGammaSysErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphCombIncGammaStatErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);	
	graphCombIncGammaSysErr4080Plot->Draw("E2same");
	graphCombIncGammaStatErr4080Plot->Draw("p,E1Z,same");

	TLatex *labelScalingIncGamma0020 = new TLatex(11.25,1.8E-3,"x 10^{2}");
	SetStyleTLatex( labelScalingIncGamma0020, 0.85*textsizeLabelsIncGamma,4,colorComb0020,42,kFALSE);
	labelScalingIncGamma0020->Draw();
	
	TLatex *labelScalingIncGamma2040 = new TLatex(11.25,1E-4,"x 10^{1}");
	SetStyleTLatex( labelScalingIncGamma2040, 0.85*textsizeLabelsIncGamma,4,colorComb2040,42,kFALSE);
	labelScalingIncGamma2040->Draw();
	
	TLatex *labelScalingIncGamma4080 = new TLatex(11.25,2E-6,"x 10^{0}");
	SetStyleTLatex( labelScalingIncGamma4080, 0.85*textsizeLabelsIncGamma,4,colorComb4080,42,kFALSE);
	labelScalingIncGamma4080->Draw();
	
	TLatex *labelIncGammaColl = new TLatex(0.6,0.91,collisionSystem.Data());
	SetStyleTLatex( labelIncGammaColl, 0.85*textsizeLabelsIncGamma,4);
	labelIncGammaColl->Draw();
	SetStyleHisto(histoFitDummyPlotting, widthCommonFit, 5, kGray+1);
	
	TLegend* legendIncGamma = new TLegend(0.6,0.9-1.*0.85*textsizeLabelsIncGamma*3,0.9,0.9);
	legendIncGamma->SetFillStyle(0);
	legendIncGamma->SetFillColor(0);
	legendIncGamma->SetLineColor(0);
	legendIncGamma->SetTextSize(0.85*textsizeLabelsIncGamma);
	legendIncGamma->SetMargin(0.2);
	legendIncGamma->SetTextFont(42);
	legendIncGamma->AddEntry(graphCombIncGammaSysErr0020Plot,"  0-20% ALICE","pf");
	legendIncGamma->AddEntry(graphCombIncGammaSysErr2040Plot,"20-40% ALICE","pf");
	legendIncGamma->AddEntry(graphCombIncGammaSysErr4080Plot,"40-80% ALICE","pf");
// 	legendIncGamma->AddEntry(histoFitDummyPlotting,"Fits to #gamma_{inc}","l");
	legendIncGamma->Draw();

	canvasIncGamma->Print(Form("%s/IncGammaSpectrum.%s",outputDir.Data(),suffix.Data()));

	//*******************************************************************************************************************************************
	//*************************************************** Plotting direct Gamma Spectrum ********************************************************
	//*******************************************************************************************************************************************
		
	TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr0020Plot;
	TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr0020Plot;
	TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr0020Plot;
	TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr0020ArPlot;
	if (graphCombDirGammaSpectrumSumErr0020) graphCombDirGammaSpectrumSumErr0020Plot 		= ScaleGraph(graphCombDirGammaSpectrumSumErr0020,100);
	if (graphCombDirGammaSpectrumStatErr0020) graphCombDirGammaSpectrumStatErr0020Plot 		= ScaleGraph(graphCombDirGammaSpectrumStatErr0020,100);
	if (graphCombDirGammaSpectrumStatErr0020Plot) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErr0020Plot);
	if (graphCombDirGammaSpectrumSystErr0020) graphCombDirGammaSpectrumSystErr0020Plot 		= ScaleGraph(graphCombDirGammaSpectrumSystErr0020,100);
	if (graphCombDirGammaSpectrumSumErr0020Ar) graphCombDirGammaSpectrumSumErr0020ArPlot 	= ScaleGraph(graphCombDirGammaSpectrumSumErr0020Ar,100);
	TGraphAsymmErrors* graphTheoryNLO0020Plot 												= ScaleGraph(graphTheoryNLO0020,100);
	TGraphAsymmErrors* graphTheoryEPS090020Plot 											= ScaleGraph(graphTheoryEPS090020,100);
	TGraphAsymmErrors* graphTheoryCT100020Plot 												= ScaleGraph(graphTheoryCT100020,100);
	TGraphErrors* graphTheoryPromptMcGill0020Plot 											= ScaleGraph(graphTheoryPromptMcGill0020,100);
	TH1D* histoFitThermalGamma0020Stat														= (TH1D*)fitThermalGamma0020Stat->GetHistogram();
	histoFitThermalGamma0020Stat->Scale(100);
	TH1D* histoFitPureThermalGamma0020Stat													= NULL;
	if (graphCombThermalGammaSpectrumStatErr0020){
		histoFitPureThermalGamma0020Stat													= (TH1D*)fitPureThermalGamma0020Stat->GetHistogram();
		histoFitPureThermalGamma0020Stat->Scale(100);
	}	
	
	TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2040Plot;
	TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr2040Plot;
	TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr2040Plot;
	TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2040ArPlot;
	if (graphCombDirGammaSpectrumSumErr2040) graphCombDirGammaSpectrumSumErr2040Plot 		= ScaleGraph(graphCombDirGammaSpectrumSumErr2040,10);
	if (graphCombDirGammaSpectrumStatErr2040) graphCombDirGammaSpectrumStatErr2040Plot 		= ScaleGraph(graphCombDirGammaSpectrumStatErr2040,10);
	if (graphCombDirGammaSpectrumStatErr2040Plot) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErr2040Plot);
	if (graphCombDirGammaSpectrumSystErr2040) graphCombDirGammaSpectrumSystErr2040Plot 		= ScaleGraph(graphCombDirGammaSpectrumSystErr2040,10);
	if (graphCombDirGammaSpectrumSumErr2040Ar) graphCombDirGammaSpectrumSumErr2040ArPlot	= ScaleGraph(graphCombDirGammaSpectrumSumErr2040Ar,10);
	TGraphAsymmErrors* graphTheoryNLO2040Plot 												= ScaleGraph(graphTheoryNLO2040,10);
	TGraphAsymmErrors* graphTheoryEPS092040Plot 											= ScaleGraph(graphTheoryEPS092040,10);
	TGraphAsymmErrors* graphTheoryCT102040Plot 												= ScaleGraph(graphTheoryCT102040,10);
	TGraphErrors* graphTheoryPromptMcGill2040Plot 											= ScaleGraph(graphTheoryPromptMcGill2040,10);
	TH1D* histoFitThermalGamma2040Stat														= (TH1D*)fitThermalGamma2040Stat->GetHistogram();
	histoFitThermalGamma2040Stat->Scale(10);
	
	TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr4080Plot;
	TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr4080Plot;
	TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr4080Plot;
	TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr4080ArPlot;
	if (graphCombDirGammaSpectrumSumErr4080) graphCombDirGammaSpectrumSumErr4080Plot 		= ScaleGraph(graphCombDirGammaSpectrumSumErr4080,1);
	if (graphCombDirGammaSpectrumStatErr4080) graphCombDirGammaSpectrumStatErr4080Plot 		= ScaleGraph(graphCombDirGammaSpectrumStatErr4080,1);
	if (graphCombDirGammaSpectrumStatErr4080Plot) ProduceGraphAsymmWithoutXErrors(graphCombDirGammaSpectrumStatErr4080Plot);
	if (graphCombDirGammaSpectrumSystErr4080) graphCombDirGammaSpectrumSystErr4080Plot 		= ScaleGraph(graphCombDirGammaSpectrumSystErr4080,1);
	if (graphCombDirGammaSpectrumSumErr4080Ar) graphCombDirGammaSpectrumSumErr4080ArPlot 	= ScaleGraph(graphCombDirGammaSpectrumSumErr4080Ar,1);
	TGraphAsymmErrors* graphTheoryNLO4080Plot 												= ScaleGraph(graphTheoryNLO4080,1);
	TGraphAsymmErrors* graphTheoryEPS094080Plot 											= ScaleGraph(graphTheoryEPS094080,1);
	TGraphAsymmErrors* graphTheoryCT104080Plot 												= ScaleGraph(graphTheoryCT104080,1);

	
	TCanvas *canvasDirGamma = new TCanvas("canvasDirGamma","",10,10,1200,1400);  // gives the page size		
	DrawGammaCanvasSettings( canvasDirGamma, 0.165, 0.01, 0.01, 0.07);
	canvasDirGamma->SetLogy();
	canvasDirGamma->SetLogx();

	Int_t textSizeLabelsPixelDirGam = 48;
	Double_t textsizeLabelsDirGamma = 0;
	if (canvasDirGamma->XtoPixel(canvasDirGamma->GetX2()) < canvasDirGamma->YtoPixel(canvasDirGamma->GetY1())){
		textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGamma->XtoPixel(canvasDirGamma->GetX2()) ;
	} else {
		textsizeLabelsDirGamma = (Double_t)textSizeLabelsPixelDirGam/canvasDirGamma->YtoPixel(canvasDirGamma->GetY1());
	}
		
	TH2D *dummyDirGamma ;
	dummyDirGamma = new TH2D("dummyDirGamma", "dummyDirGamma", 1000, 0., 22, 1000., 1.2e-8,1.5e5);
	SetStyleHistoTH2ForGraphs( dummyDirGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
	dummyDirGamma->GetXaxis()->SetLabelOffset(-0.015);
	dummyDirGamma->GetXaxis()->SetTickLength(0.025);
	dummyDirGamma->GetYaxis()->SetTickLength(0.025);
	dummyDirGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
// 	dummyDirGamma->GetXaxis()->SetRangeUser(0,16);
	dummyDirGamma->DrawCopy();

	TLatex *labelScalingDirGamma0020 = new TLatex(11.2,1.3E-3,"x 10^{2}");
	SetStyleTLatex( labelScalingDirGamma0020, 0.85*textsizeLabelsDirGamma,4,colorComb0020,42,kFALSE);
	labelScalingDirGamma0020->Draw();
	
	TLatex *labelScalingDirGamma2040 = new TLatex(11.2,6.0E-5,"x 10^{1}");
	SetStyleTLatex( labelScalingDirGamma2040, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
	labelScalingDirGamma2040->Draw();
	
	TLatex *labelScalingDirGamma4080 = new TLatex(11.2,7.5E-7,"x 10^{0}");
	SetStyleTLatex( labelScalingDirGamma4080, 0.85*textsizeLabelsDirGamma,4,colorComb4080,42,kFALSE);
	labelScalingDirGamma4080->Draw();
	
	TLatex *labelDirGammaColl = new TLatex(0.25,0.94,Form("%s",collisionSystem.Data()));
	SetStyleTLatex( labelDirGammaColl, 0.85*textsizeLabelsDirGamma,4);
	labelDirGammaColl->Draw();
	SetStyleHisto(histoFitDummyPlotting, widthCommonFit, 5, kGray+1);
// 	
	TLegend* legendDirGamma = new TLegend(0.24,0.93-1.*0.85*textsizeLabelsDirGamma*3,0.24+0.21,0.93);
	legendDirGamma->SetFillStyle(0);
	legendDirGamma->SetFillColor(0);
	legendDirGamma->SetLineColor(0);
	legendDirGamma->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGamma->SetMargin(0.25);
	legendDirGamma->SetTextFont(42);
	legendDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr0020Plot,"  0-20% ALICE","pf");
	legendDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
	legendDirGamma->AddEntry(graphCombDirGammaSpectrumSystErr4080Plot,"40-80% ALICE","pf");
	legendDirGamma->Draw();

	if (graphCombDirGammaSpectrumSystErr0020Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr0020Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr0020Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);	
		graphCombDirGammaSpectrumStatErr0020Plot->Draw("p,E1Z,same");
	}

	if (graphCombDirGammaSpectrumSystErr2040Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
	}	
	if (graphCombDirGammaSpectrumStatErr2040Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
		graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr2040ArPlot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
		graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
	}	
	if (graphCombDirGammaSpectrumSystErr4080Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr4080Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr4080Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);	
		graphCombDirGammaSpectrumStatErr4080Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr4080ArPlot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr4080ArPlot , 1, 3, colorComb4080, colorComb4080, 1.8, kTRUE);
		graphCombDirGammaSpectrumSumErr4080ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr4080ArPlot);
	}	
	
	canvasDirGamma->Print(Form("%s/DirGammaSpectrum.%s",outputDir.Data(),suffix.Data()));
	
	SetStyleGammaNLOTGraphWithBand( graphTheoryEPS090020Plot, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
	graphTheoryEPS090020Plot->Draw("p3lsame");

	SetStyleGammaNLOTGraphWithBand( graphTheoryCT100020Plot, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
	graphTheoryCT100020Plot->Draw("p3lsame");
	cout << "JETPHOX CT10 0-20%" << endl;

	SetStyleGammaNLOTGraphWithBand( graphTheoryNLO0020Plot, 3.0, 1, colorNLOcalc, fillStyleNLO, colorNLOcalc, 0);
	graphTheoryNLO0020Plot->Draw("p3lsame");
	SetStyleGammaNLOTGraphWithBand( graphTheoryPromptMcGill0020Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
	graphTheoryPromptMcGill0020Plot->Draw("lsame");
	
	if (graphCombDirGammaSpectrumSystErr0020Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr0020Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr0020Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);	
		graphCombDirGammaSpectrumStatErr0020Plot->Draw("p,E1Z,same");
	}

	SetStyleGammaNLOTGraphWithBand( graphTheoryEPS092040Plot, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
	graphTheoryEPS092040Plot->Draw("p3lsame");

	SetStyleGammaNLOTGraphWithBand( graphTheoryCT102040Plot, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
	graphTheoryCT102040Plot->Draw("p3lsame");

	SetStyleGammaNLOTGraphWithBand( graphTheoryNLO2040Plot, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);
	graphTheoryNLO2040Plot->Draw("p3lsame");
	SetStyleGammaNLOTGraphWithBand( graphTheoryPromptMcGill2040Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
	graphTheoryPromptMcGill2040Plot->Draw("lsame");

	if (graphCombDirGammaSpectrumSystErr2040Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
	}	
	if (graphCombDirGammaSpectrumStatErr2040Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
		graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr2040ArPlot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
		graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
	}	

	SetStyleGammaNLOTGraphWithBand( graphTheoryEPS094080Plot, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
	graphTheoryEPS094080Plot->Draw("p3lsame");

	SetStyleGammaNLOTGraphWithBand( graphTheoryCT104080Plot, 1.0, 1, colorCT10calc, fillStyleCT10, colorCT10calc, 0);
	graphTheoryCT104080Plot->Draw("p3lsame");

	SetStyleGammaNLOTGraphWithBand( graphTheoryNLO4080Plot, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);
	graphTheoryNLO4080Plot->Draw("p3lsame");

	if (graphCombDirGammaSpectrumSystErr4080Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr4080Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr4080Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);	
		graphCombDirGammaSpectrumStatErr4080Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr4080ArPlot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr4080ArPlot , 1, 3, colorComb4080, colorComb4080, 1.8, kTRUE);
		graphCombDirGammaSpectrumSumErr4080ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr4080ArPlot);
	}	

	TLegend* legendDirGammaTheory = new TLegend(0.525,0.93-1.*0.85*textsizeLabelsDirGamma*3,0.525+0.21,0.93);
	legendDirGammaTheory->SetFillStyle(0);
	legendDirGammaTheory->SetFillColor(0);
	legendDirGammaTheory->SetLineColor(0);
	legendDirGammaTheory->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaTheory->SetMargin(0.25);
	legendDirGammaTheory->SetTextFont(42);
	legendDirGammaTheory->AddEntry(graphTheoryNLO0020Plot,"PDF: CTEQ6M5, FF: GRV ","l");
	legendDirGammaTheory->AddEntry(graphTheoryPromptMcGill2040Plot,"(n)PDF: CTEQ6.1M/EPS09,","l");
	legendDirGammaTheory->AddEntry((TObject*)0,"FF: BFG2","");
	legendDirGammaTheory->Draw();

	TLegend* legendDirGammaTheory2 = new TLegend(0.525,0.93-1.*0.85*textsizeLabelsDirGamma*7-0.02,0.525+0.21,0.93-1.*0.85*textsizeLabelsDirGamma*3-0.02);
	legendDirGammaTheory2->SetFillStyle(0);
	legendDirGammaTheory2->SetFillColor(0);
	legendDirGammaTheory2->SetLineColor(0);
	legendDirGammaTheory2->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaTheory2->SetMargin(0.25);
	legendDirGammaTheory2->SetTextFont(42);
	legendDirGammaTheory2->AddEntry((TObject*)0,"#it{JETPHOX}","");
	legendDirGammaTheory2->AddEntry(graphTheoryCT100020Plot,"PDF: CT10, FF: BFG2","f");
	legendDirGammaTheory2->AddEntry(graphTheoryEPS090020Plot,"nPDF: EPS09, FF: BFG2","f");
	legendDirGammaTheory2->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
	legendDirGammaTheory2->Draw();
	
	canvasDirGamma->Print(Form("%s/DirGammaSpectrum_withNLO.%s",outputDir.Data(),suffix.Data()));

	//************************************************************************************************************************************
	
	TGraphErrors* graphTheoryMcGill0020Plot	= ScaleGraph(graphTheoryMcGill0020,100);
	TGraphErrors* graphTheoryMcGill2040Plot	= ScaleGraph(graphTheoryMcGill2040,10);

	TGraph* graphTheoryPHSD0020Plot	= ScaleGraph(graphTheoryPHSD0020,100);
	graphTheoryPHSD0020Plot->RemovePoint(0);
	TGraph* graphTheoryPHSD2040Plot	= ScaleGraph(graphTheoryPHSD2040,10);
	graphTheoryPHSD2040Plot->RemovePoint(0);
	TGraph* graphTheoryPHSD4080Plot	= ScaleGraph(graphTheoryPHSD4080,1);
	graphTheoryPHSD4080Plot->RemovePoint(0);

	TGraph* graphTheoryChatterjee0020Plot	= ScaleGraph(graphTheoryChatterjee0020,100);
	TGraph* graphTheoryChatterjee2040Plot	= ScaleGraph(graphTheoryChatterjee2040,10);
	TGraph* graphTheoryChatterjee4060Plot	= ScaleGraph(graphTheoryChatterjee4060,1);

	TGraph* graphTheoryHees0020Plot	= ScaleGraph(graphTheoryHees0020,100);
	TGraph* graphTheoryHees2040Plot	= ScaleGraph(graphTheoryHees2040,10);
	TGraph* graphTheoryHees4080Plot	= ScaleGraph(graphTheoryHees4080,1);

	TGraph* graphTheoryHe0020Plot	= ScaleGraph(graphTheoryHe0020,100);
	TGraph* graphTheoryHe2040Plot	= ScaleGraph(graphTheoryHe2040,10);
	TGraph* graphTheoryHe4080Plot	= ScaleGraph(graphTheoryHe4080,1);
	
	TH2D *dummyDirGammaTheory ;
	dummyDirGammaTheory = new TH2D("dummyDirGammaTheory", "dummyDirGamma", 1000, 0., 22, 1000., 4e-8,9e2);
	SetStyleHistoTH2ForGraphs( dummyDirGammaTheory, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
	dummyDirGammaTheory->GetXaxis()->SetLabelOffset(-0.015);
	dummyDirGammaTheory->GetXaxis()->SetTickLength(0.025);
	dummyDirGammaTheory->GetYaxis()->SetTickLength(0.025);
	dummyDirGammaTheory->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyDirGammaTheory->DrawCopy();

	SetStyleGammaNLOTGraphWithBand( graphTheoryMcGill0020Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
	SetStyleGammaNLOTGraphWithBand( graphTheoryPHSD0020Plot, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
	SetStyleGammaNLOTGraphWithBand( graphTheoryChatterjee0020Plot, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
	SetStyleGammaNLOTGraphWithBand( graphTheoryHees0020Plot, 3.0, styleHees, colorHees, 3015, colorHees, 0);
	SetStyleGammaNLOTGraphWithBand( graphTheoryHe0020Plot, 3.0, styleHe, colorHe, 3015, colorHe, 0);
	if (graphCombDirGammaSpectrumSystErr0020Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
	}
	if (graphCombDirGammaSpectrumStatErr0020Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);	
	}
	SetStyleGammaNLOTGraphWithBand( graphTheoryMcGill2040Plot, 3.0, styleNLOMcGill, colorNLOMcGill, 3001, colorNLOMcGill, 0);
	SetStyleGammaNLOTGraphWithBand( graphTheoryPHSD2040Plot, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
	SetStyleGammaNLOTGraphWithBand( graphTheoryChatterjee2040Plot, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
	SetStyleGammaNLOTGraphWithBand( graphTheoryHees2040Plot, 3.0, styleHees, colorHees, 3015, colorHees, 0);
	SetStyleGammaNLOTGraphWithBand( graphTheoryHe2040Plot, 3.0, styleHe, colorHe, 3015, colorHe, 0);
	if (graphCombDirGammaSpectrumSystErr2040Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
	}	
	if (graphCombDirGammaSpectrumStatErr2040Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
	}
	if (graphCombDirGammaSpectrumSumErr2040ArPlot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , styleNLOMcGill, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
	}	
	SetStyleGammaNLOTGraphWithBand( graphTheoryPHSD4080Plot, 3.0, stylePHSD, colorPHSD, 3015, colorPHSD, 0);
// 	SetStyleGammaNLOTGraphWithBand( graphTheoryChatterjee4060Plot, 3.0, styleChatterjee, colorChatterjee, 3015, colorChatterjee, 0);
	SetStyleGammaNLOTGraphWithBand( graphTheoryHees4080Plot, 3.0, styleHees, colorHees, 3015, colorHees, 0);
	SetStyleGammaNLOTGraphWithBand( graphTheoryHe4080Plot, 3.0, styleHe, colorHe, 3015, colorHe, 0);
	if (graphCombDirGammaSpectrumSystErr4080Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
	}
	if (graphCombDirGammaSpectrumStatErr4080Plot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);	
	}
	if (graphCombDirGammaSpectrumSumErr4080ArPlot){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr4080ArPlot , 1, 3, colorComb4080, colorComb4080, 1.8, kTRUE);
	}	

	labelScalingDirGamma0020->Draw();
	labelScalingDirGamma2040->Draw();
	labelScalingDirGamma4080->Draw();
	TLatex *labelDirGammaCollRedX = new TLatex(0.585,0.93,Form("%s",collisionSystem.Data()));
	SetStyleTLatex( labelDirGammaCollRedX, 0.85*textsizeLabelsDirGamma,4);
	labelDirGammaCollRedX->Draw();

	
	TLegend* legendDirGammaRedX = new TLegend(0.58,0.92-1.1*0.85*textsizeLabelsDirGamma*3,0.585+0.21,0.92);
	legendDirGammaRedX->SetFillStyle(0);
	legendDirGammaRedX->SetFillColor(0);
	legendDirGammaRedX->SetLineColor(0);
	legendDirGammaRedX->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaRedX->SetMargin(0.3);
	legendDirGammaRedX->SetTextFont(42);
	legendDirGammaRedX->AddEntry(graphCombDirGammaSpectrumSystErr0020Plot,"  0-20% ALICE","pf");
	legendDirGammaRedX->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
	legendDirGammaRedX->AddEntry(graphCombDirGammaSpectrumSystErr4080Plot,"40-80% ALICE","pf");
	legendDirGammaRedX->Draw();
	
	if (graphCombDirGammaSpectrumSystErr0020Plot){
		graphCombDirGammaSpectrumSystErr0020Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr0020Plot){
		graphCombDirGammaSpectrumStatErr0020Plot->Draw("p,E1Z,same");
	}

	if (graphCombDirGammaSpectrumSystErr2040Plot){
		graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
	}	
	if (graphCombDirGammaSpectrumStatErr2040Plot){
		graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr2040ArPlot){
		graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
	}	
	if (graphCombDirGammaSpectrumSystErr4080Plot){
		graphCombDirGammaSpectrumSystErr4080Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr4080Plot){
		graphCombDirGammaSpectrumStatErr4080Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr4080ArPlot){
		graphCombDirGammaSpectrumSumErr4080ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr4080ArPlot);
	}	

	canvasDirGamma->Print(Form("%s/DirGammaSpectrum_LogX.%s",outputDir.Data(),suffix.Data()));	

	TLegend* legendDirGammaTheoryPlusMcGillDiffY = new TLegend(0.20,0.382-0.95*0.85*textsizeLabelsDirGamma*9,0.5,0.382);
	legendDirGammaTheoryPlusMcGillDiffY->SetFillStyle(0);
	legendDirGammaTheoryPlusMcGillDiffY->SetFillColor(0);
	legendDirGammaTheoryPlusMcGillDiffY->SetLineColor(0);
// 	legendDirGammaTheoryPlusMcGillDiffY->SetNColumns(2);
	legendDirGammaTheoryPlusMcGillDiffY->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaTheoryPlusMcGillDiffY->SetMargin(0.19);
	legendDirGammaTheoryPlusMcGillDiffY->SetTextFont(42);
	legendDirGammaTheoryPlusMcGillDiffY->AddEntry(graphTheoryMcGill0020Plot,"Paquet et al.","l");
	legendDirGammaTheoryPlusMcGillDiffY->AddEntry((TObject*)0,"arXiv:1509.06738","");
	legendDirGammaTheoryPlusMcGillDiffY->AddEntry(graphTheoryPHSD0020Plot,"Linnyk et al.","l");
	legendDirGammaTheoryPlusMcGillDiffY->AddEntry((TObject*)0,"arXiv:1504.05699","");
	legendDirGammaTheoryPlusMcGillDiffY->AddEntry(graphTheoryHe0020Plot,"v. Hees et al.","l");	
	legendDirGammaTheoryPlusMcGillDiffY->AddEntry((TObject*)0,"NPA 933(2015) 256","");
	legendDirGammaTheoryPlusMcGillDiffY->AddEntry(graphTheoryChatterjee0020Plot,"Chatterjee et al.","l");
	legendDirGammaTheoryPlusMcGillDiffY->AddEntry((TObject*)0,"PRC 85(2012) 064910","");
	legendDirGammaTheoryPlusMcGillDiffY->AddEntry((TObject*)0,"+ JHEP 1305(2013) 030","");
	legendDirGammaTheoryPlusMcGillDiffY->Draw();

	graphTheoryMcGill0020Plot->Draw("p3lsame");
	graphTheoryPHSD0020Plot->Draw("p3lsame");
	graphTheoryChatterjee0020Plot->Draw("p3lsame");
// 	graphTheoryHees0020Plot->Draw("p3lsame");
	graphTheoryHe0020Plot->Draw("p3lsame");	
	if (graphCombDirGammaSpectrumSystErr0020Plot){
		graphCombDirGammaSpectrumSystErr0020Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr0020Plot){
		graphCombDirGammaSpectrumStatErr0020Plot->Draw("p,E1Z,same");
	}
	graphTheoryMcGill2040Plot->Draw("p3lsame");
	graphTheoryPHSD2040Plot->Draw("p3lsame");
	graphTheoryChatterjee2040Plot->Draw("p3lsame");
// 	graphTheoryHees2040Plot->Draw("p3lsame");
	graphTheoryHe2040Plot->Draw("p3lsame");

	if (graphCombDirGammaSpectrumSystErr2040Plot){
		graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
	}	
	if (graphCombDirGammaSpectrumStatErr2040Plot){
		graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr2040ArPlot){
		graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
	}	

	graphTheoryPHSD4080Plot->Draw("p3lsame");
// 	graphTheoryChatterjee4060Plot->Draw("p3lsame");
// 	graphTheoryHees4080Plot->Draw("p3lsame");
	graphTheoryHe4080Plot->Draw("p3lsame");
	if (graphCombDirGammaSpectrumSystErr4080Plot){
		graphCombDirGammaSpectrumSystErr4080Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr4080Plot){
		graphCombDirGammaSpectrumStatErr4080Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr4080ArPlot){
		graphCombDirGammaSpectrumSumErr4080ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr4080ArPlot);
	}	
	
	canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusMcGill_LogX.%s",outputDir.Data(),suffix.Data()));

	//*****************************************************************************************************************
	//******************************* Plotting spectrum with reduced x range ******************************************
	//*****************************************************************************************************************
	// Cut plotting graphs
	TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr0020Plot2 = NULL;
	if (graphCombDirGammaSpectrumSystErr0020Plot) {
		graphCombDirGammaSpectrumSystErr0020Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr0020Plot->Clone("graphCombDirGammaSpectrumSystErr0020Plot2");
		while (graphCombDirGammaSpectrumSystErr0020Plot2->GetX()[graphCombDirGammaSpectrumSystErr0020Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumSystErr0020Plot2->GetN()>0)
			graphCombDirGammaSpectrumSystErr0020Plot2->RemovePoint(graphCombDirGammaSpectrumSystErr0020Plot2->GetN()-1);
		if (graphCombDirGammaSpectrumSystErr0020Plot2->GetN()==0) graphCombDirGammaSpectrumSystErr0020Plot2 = NULL;
	}
	TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr0020Plot2 = NULL;
	if (graphCombDirGammaSpectrumStatErr0020Plot) {
		graphCombDirGammaSpectrumStatErr0020Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr0020Plot->Clone("graphCombDirGammaSpectrumStatErr0020Plot2");
		while (graphCombDirGammaSpectrumStatErr0020Plot2->GetX()[graphCombDirGammaSpectrumStatErr0020Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumStatErr0020Plot2->GetN()>0)
			graphCombDirGammaSpectrumStatErr0020Plot2->RemovePoint(graphCombDirGammaSpectrumStatErr0020Plot2->GetN()-1);
		if (graphCombDirGammaSpectrumStatErr0020Plot2->GetN()==0) graphCombDirGammaSpectrumStatErr0020Plot2 = NULL;
	}	
	TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr2040Plot2 = NULL;
	if (graphCombDirGammaSpectrumSystErr2040Plot) {
		graphCombDirGammaSpectrumSystErr2040Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr2040Plot->Clone("graphCombDirGammaSpectrumSystErr2040Plot2");
		while (graphCombDirGammaSpectrumSystErr2040Plot2->GetX()[graphCombDirGammaSpectrumSystErr2040Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumSystErr2040Plot2->GetN()>0)
			graphCombDirGammaSpectrumSystErr2040Plot2->RemovePoint(graphCombDirGammaSpectrumSystErr2040Plot2->GetN()-1);
		if (graphCombDirGammaSpectrumSystErr2040Plot2->GetN()==0) graphCombDirGammaSpectrumSystErr2040Plot2 = NULL;
	}
	TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr2040Plot2 = NULL;
	if (graphCombDirGammaSpectrumStatErr2040Plot) {
		graphCombDirGammaSpectrumStatErr2040Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr2040Plot->Clone("graphCombDirGammaSpectrumStatErr2040Plot2");
		while (graphCombDirGammaSpectrumStatErr2040Plot2->GetX()[graphCombDirGammaSpectrumStatErr2040Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumStatErr2040Plot2->GetN()>0)
			graphCombDirGammaSpectrumStatErr2040Plot2->RemovePoint(graphCombDirGammaSpectrumStatErr2040Plot2->GetN()-1);
		if (graphCombDirGammaSpectrumStatErr2040Plot2->GetN()==0) graphCombDirGammaSpectrumStatErr2040Plot2 = NULL;
	}	

	TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2040ArPlot2 = NULL;
	if (graphCombDirGammaSpectrumSumErr2040ArPlot) {
		graphCombDirGammaSpectrumSumErr2040ArPlot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSumErr2040ArPlot->Clone("graphCombDirGammaSpectrumSumErr2040ArPlot2");
		while (graphCombDirGammaSpectrumSumErr2040ArPlot2->GetX()[graphCombDirGammaSpectrumSumErr2040ArPlot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumSumErr2040ArPlot2->GetN()>0)
			graphCombDirGammaSpectrumSumErr2040ArPlot2->RemovePoint(graphCombDirGammaSpectrumSumErr2040ArPlot2->GetN()-1);
		if (graphCombDirGammaSpectrumSumErr2040ArPlot2->GetN()==0) graphCombDirGammaSpectrumSumErr2040ArPlot2 = NULL;
	}		
	
	TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr4080Plot2 = NULL;
	if (graphCombDirGammaSpectrumSystErr4080Plot) {
		graphCombDirGammaSpectrumSystErr4080Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr4080Plot->Clone("graphCombDirGammaSpectrumSystErr4080Plot2");
		while (graphCombDirGammaSpectrumSystErr4080Plot2->GetX()[graphCombDirGammaSpectrumSystErr4080Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumSystErr4080Plot2->GetN()>0)
			graphCombDirGammaSpectrumSystErr4080Plot2->RemovePoint(graphCombDirGammaSpectrumSystErr4080Plot2->GetN()-1);
		if (graphCombDirGammaSpectrumSystErr4080Plot2->GetN()==0) graphCombDirGammaSpectrumSystErr4080Plot2 = NULL;
	}
	TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr4080Plot2 = NULL;
	if (graphCombDirGammaSpectrumStatErr4080Plot) {
		graphCombDirGammaSpectrumStatErr4080Plot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr4080Plot->Clone("graphCombDirGammaSpectrumStatErr4080Plot2");
		while (graphCombDirGammaSpectrumStatErr4080Plot2->GetX()[graphCombDirGammaSpectrumStatErr4080Plot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumStatErr4080Plot2->GetN()>0)
			graphCombDirGammaSpectrumStatErr4080Plot2->RemovePoint(graphCombDirGammaSpectrumStatErr4080Plot2->GetN()-1);
		if (graphCombDirGammaSpectrumStatErr4080Plot2->GetN()==0) graphCombDirGammaSpectrumStatErr4080Plot2 = NULL;
	}	
	
	TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr4080ArPlot2 = NULL;
	if (graphCombDirGammaSpectrumSumErr4080ArPlot) {
		graphCombDirGammaSpectrumSumErr4080ArPlot2 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSumErr4080ArPlot->Clone("graphCombDirGammaSpectrumSumErr4080ArPlot2");
		while (graphCombDirGammaSpectrumSumErr4080ArPlot2->GetX()[graphCombDirGammaSpectrumSumErr4080ArPlot2->GetN()-1] > 4.1 && graphCombDirGammaSpectrumSumErr4080ArPlot2->GetN()>0)
			graphCombDirGammaSpectrumSumErr4080ArPlot2->RemovePoint(graphCombDirGammaSpectrumSumErr4080ArPlot2->GetN()-1);
		if (graphCombDirGammaSpectrumSumErr4080ArPlot2->GetN()==0) graphCombDirGammaSpectrumSumErr4080ArPlot2 = NULL;
	}	

	
	canvasDirGamma->SetLogx(0);
	TH2D *dummDirGammaRedX ;
	dummDirGammaRedX = new TH2D("dummDirGammaRedX", "dummDirGammaRedX", 100000, 0., 4.25, 1000., 9e-6,9e2);
	SetStyleHistoTH2ForGraphs( dummDirGammaRedX, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
	dummDirGammaRedX->GetXaxis()->SetLabelOffset(-0.002);
	dummDirGammaRedX->GetXaxis()->SetTickLength(0.025);
	dummDirGammaRedX->GetYaxis()->SetTickLength(0.025);
	dummDirGammaRedX->GetXaxis()->SetRangeUser(0.5,doubleRatioX[1]);
	dummDirGammaRedX->DrawCopy();

	
	SetStyleHisto(histoFitThermalGamma0020Stat, 3, styleFit, colorComb0020+1 );
	SetStyleHisto(histoFitThermalGamma2040Stat, 3, styleFit, colorComb2040+1 );

	TLatex *labelScalingDirGamma0020_3 = new TLatex(3.65,4E-1,"x 10^{2}");
	SetStyleTLatex( labelScalingDirGamma0020_3, 0.85*textsizeLabelsDirGamma,4,colorComb0020,42,kFALSE);
	labelScalingDirGamma0020_3->Draw();
	
	TLatex *labelScalingDirGamma2040_3 = new TLatex(3.65,1.7E-2,"x 10^{1}");
	SetStyleTLatex( labelScalingDirGamma2040_3, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
	labelScalingDirGamma2040_3->Draw();
	
	TLatex *labelScalingDirGamma4080_3 = new TLatex(3.65,4E-4,"x 10^{0}");
	SetStyleTLatex( labelScalingDirGamma4080_3, 0.85*textsizeLabelsDirGamma,4,colorComb4080,42,kFALSE);
	labelScalingDirGamma4080_3->Draw();
	TLatex *labelDirGammaCollRedX2 = new TLatex(0.41,0.93,Form("%s",collisionSystem.Data()));
	SetStyleTLatex( labelDirGammaCollRedX2, 0.85*textsizeLabelsDirGamma,4);
	labelDirGammaCollRedX2->Draw();

	
	TLegend* legendDirGammaRedXWithFit = new TLegend(0.455,0.915-0.9*0.85*textsizeLabelsDirGamma*3,0.555+0.21,0.915);
	legendDirGammaRedXWithFit->SetFillStyle(0);
	legendDirGammaRedXWithFit->SetFillColor(0);
	legendDirGammaRedXWithFit->SetLineColor(0);
	legendDirGammaRedXWithFit->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaRedXWithFit->SetMargin(0.2);
	legendDirGammaRedXWithFit->SetTextFont(42);
	legendDirGammaRedXWithFit->AddEntry(graphCombDirGammaSpectrumSystErr0020Plot,"  0-20% ALICE","pf");
	legendDirGammaRedXWithFit->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
	legendDirGammaRedXWithFit->AddEntry(graphCombDirGammaSpectrumSystErr4080Plot,"40-80% ALICE","pf");
	legendDirGammaRedXWithFit->Draw();

	
	if (graphCombDirGammaSpectrumSystErr0020Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr0020Plot2->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr0020Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);	
		graphCombDirGammaSpectrumStatErr0020Plot2->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSystErr2040Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot2, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr2040Plot2->Draw("E2same");
	}	
	if (graphCombDirGammaSpectrumStatErr2040Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
		graphCombDirGammaSpectrumStatErr2040Plot2->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr2040ArPlot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , styleNLOMcGill, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
		graphCombDirGammaSpectrumSumErr2040ArPlot2->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot2);
	}	
	if (graphCombDirGammaSpectrumSystErr4080Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr4080Plot2, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr4080Plot2->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr4080Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr4080Plot2, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);	
		graphCombDirGammaSpectrumStatErr4080Plot2->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr4080ArPlot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr4080ArPlot2 , 1, 3, colorComb4080, colorComb4080, 1.8, kTRUE);
		graphCombDirGammaSpectrumSumErr4080ArPlot2->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr4080ArPlot2);
	}	
	
	canvasDirGamma->Print(Form("%s/DirGammaSpectrum_ReducedX.%s",outputDir.Data(),suffix.Data()));

	
	TLatex *labelDirGammaFitFunc = new TLatex(0.755,0.93,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})");
	SetStyleTLatex( labelDirGammaFitFunc, 0.85*textsizeLabelsDirGamma,4);
	labelDirGammaFitFunc->Draw();
	TLegend* legendDirGammaExpFit = new TLegend(0.795,0.915-0.9*0.85*textsizeLabelsDirGamma*2,0.96,0.915);
	legendDirGammaExpFit->SetFillStyle(0);
	legendDirGammaExpFit->SetFillColor(0);
	legendDirGammaExpFit->SetLineColor(0);
	legendDirGammaExpFit->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaExpFit->SetMargin(0.3);
	legendDirGammaExpFit->SetTextFont(42);
	legendDirGammaExpFit->AddEntry(histoFitThermalGamma0020Stat,"  0-20%","l");
	legendDirGammaExpFit->AddEntry(histoFitThermalGamma2040Stat,"20-40%","l");
	legendDirGammaExpFit->Draw();
	
	histoFitThermalGamma0020Stat->Draw("same,l");
	histoFitThermalGamma2040Stat->Draw("same,l");

	
	canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusFit_ReducedX.%s",outputDir.Data(),suffix.Data()));
	
	TLegend* legendDirGammaTheoryPlusMcGillRedX = new TLegend(0.19,0.295-1.*0.85*textsizeLabelsDirGamma*6,0.94,0.295);
	legendDirGammaTheoryPlusMcGillRedX->SetFillStyle(0);
	legendDirGammaTheoryPlusMcGillRedX->SetFillColor(0);
	legendDirGammaTheoryPlusMcGillRedX->SetLineColor(0);
	legendDirGammaTheoryPlusMcGillRedX->SetNColumns(2);
	legendDirGammaTheoryPlusMcGillRedX->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaTheoryPlusMcGillRedX->SetMargin(0.12);
	legendDirGammaTheoryPlusMcGillRedX->SetTextFont(42);
	legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryMcGill0020Plot,"Paquet et al.","l");
	legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"","");
	legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"arXiv:1509.06738","");
	legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"","");
	legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryPHSD0020Plot,"Linnyk et al.","l");
	legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"","");
	legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"arXiv:1504.05699","");	
	legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryChatterjee0020Plot,"Chatterjee et al.","l");
	legendDirGammaTheoryPlusMcGillRedX->AddEntry(graphTheoryHe0020Plot,"v. Hees et al.","l");
	legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"PRC 85(2012) 064910","");
	legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"NPA 933(2015) 256","");
	legendDirGammaTheoryPlusMcGillRedX->AddEntry((TObject*)0,"+ JHEP 1305(2013) 030","");
	legendDirGammaTheoryPlusMcGillRedX->Draw();
	
	dummDirGammaRedX->DrawCopy("same,axis");
	labelScalingDirGamma4080_3->Draw();
	graphTheoryMcGill0020Plot->Draw("p3lsame");
	graphTheoryPHSD0020Plot->Draw("p3lsame");
	graphTheoryChatterjee0020Plot->Draw("p3lsame");
// 	graphTheoryHees0020Plot->Draw("p3lsame");
	graphTheoryHe0020Plot->Draw("p3lsame");	
	if (graphCombDirGammaSpectrumSystErr0020Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr0020Plot2->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr0020Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);	
		graphCombDirGammaSpectrumStatErr0020Plot2->Draw("p,E1Z,same");
	}
	histoFitThermalGamma0020Stat->Draw("same,l");
	
	graphTheoryMcGill2040Plot->Draw("p3lsame");
	graphTheoryPHSD2040Plot->Draw("p3lsame");
	graphTheoryChatterjee2040Plot->Draw("p3lsame");
// 	graphTheoryHees2040Plot->Draw("p3lsame");
	graphTheoryHe2040Plot->Draw("p3lsame");

	if (graphCombDirGammaSpectrumSystErr2040Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot2, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr2040Plot2->Draw("E2same");
	}	
	if (graphCombDirGammaSpectrumStatErr2040Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
		graphCombDirGammaSpectrumStatErr2040Plot2->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr2040ArPlot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot , styleNLOMcGill, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
		graphCombDirGammaSpectrumSumErr2040ArPlot2->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot2);
	}	
	histoFitThermalGamma2040Stat->Draw("same,l");
	
	graphTheoryPHSD4080Plot->Draw("p3lsame");
// 	graphTheoryChatterjee4060Plot->Draw("p3lsame");
// 	graphTheoryHees4080Plot->Draw("p3lsame");
	graphTheoryHe4080Plot->Draw("p3lsame");
	if (graphCombDirGammaSpectrumSystErr4080Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr4080Plot2, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr4080Plot2->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr4080Plot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr4080Plot2, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);	
		graphCombDirGammaSpectrumStatErr4080Plot2->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr4080ArPlot2){
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr4080ArPlot2 , 1, 3, colorComb4080, colorComb4080, 1.8, kTRUE);
		graphCombDirGammaSpectrumSumErr4080ArPlot2->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr4080ArPlot2);
	}	
	
	canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusMcGill_ReducedX.%s",outputDir.Data(),suffix.Data()));

	//*******************************************************************************************************
	//********************************PHENIX and ALICE together for 0-20% ***********************************
	//*******************************************************************************************************
// 	graphPHENIXAuAuStat0020->RemovePoint(graphPHENIXAuAuStat0020->GetN()-1);
// 	graphPHENIXAuAuSys0020->RemovePoint(graphPHENIXAuAuSys0020->GetN()-1);
	TF1* fitThermalGamma0020PHENIXStat				= FitObject("e","fitThermalGamma0020PHENIXStat","Photon",graphPHENIXAuAuStat0020,0.6,2,NULL,"QNRMEX0+");
	fitThermalGamma0020PHENIXStat->FixParameter(1,0.239);
	graphPHENIXAuAuStat0020->Fit(fitThermalGamma0020PHENIXStat,"QNRMEX0+","",0.6,2.);
	cout << WriteParameterToFile(fitThermalGamma0020PHENIXStat)<< endl;   

	TH1D* histoFitThermalGamma0020PHENIXStat		= (TH1D*)fitThermalGamma0020PHENIXStat->GetHistogram();
	SetStyleHisto(histoFitThermalGamma0020PHENIXStat, 3, styleFit, colorPHENIX );
	
	TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr0020Plot3 = NULL;
	if (graphCombDirGammaSpectrumSystErr0020Plot) {
		graphCombDirGammaSpectrumSystErr0020Plot3 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr0020Plot->Clone("graphCombDirGammaSpectrumSystErr0020Plot3");
		while (graphCombDirGammaSpectrumSystErr0020Plot3->GetX()[graphCombDirGammaSpectrumSystErr0020Plot3->GetN()-1] > 5.1 && graphCombDirGammaSpectrumSystErr0020Plot3->GetN()>0)
			graphCombDirGammaSpectrumSystErr0020Plot3->RemovePoint(graphCombDirGammaSpectrumSystErr0020Plot3->GetN()-1);
		if (graphCombDirGammaSpectrumSystErr0020Plot3->GetN()==0) graphCombDirGammaSpectrumSystErr0020Plot3 = NULL;
	}
	TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr0020Plot3 = NULL;
	if (graphCombDirGammaSpectrumStatErr0020Plot) {
		graphCombDirGammaSpectrumStatErr0020Plot3 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr0020Plot->Clone("graphCombDirGammaSpectrumStatErr0020Plot3");
		while (graphCombDirGammaSpectrumStatErr0020Plot3->GetX()[graphCombDirGammaSpectrumStatErr0020Plot3->GetN()-1] > 5.1 && graphCombDirGammaSpectrumStatErr0020Plot3->GetN()>0)
			graphCombDirGammaSpectrumStatErr0020Plot3->RemovePoint(graphCombDirGammaSpectrumStatErr0020Plot3->GetN()-1);
		if (graphCombDirGammaSpectrumStatErr0020Plot3->GetN()==0) graphCombDirGammaSpectrumStatErr0020Plot3 = NULL;
	}	
	
// 	cout << "here" << endl;
// 	TGraphAsymmErrors* graphCombThermalGammaSpectrumSystErr0020Plot2 = NULL;
// 	if (graphCombThermalGammaSpectrumSysErr0020) {
// 		graphCombThermalGammaSpectrumSystErr0020Plot2 = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumSysErr0020->Clone("graphCombThermalGammaSpectrumSystErr0020Plot2");
// 		while (graphCombThermalGammaSpectrumSystErr0020Plot2->GetX()[graphCombThermalGammaSpectrumSystErr0020Plot2->GetN()-1] > 4.1 && graphCombThermalGammaSpectrumSystErr0020Plot2->GetN()>0)
// 			graphCombThermalGammaSpectrumSystErr0020Plot2->RemovePoint(graphCombThermalGammaSpectrumSystErr0020Plot2->GetN()-1);
// 		if (graphCombThermalGammaSpectrumSystErr0020Plot2->GetN()==0) graphCombThermalGammaSpectrumSystErr0020Plot2 = NULL;
// 	}
// 	cout << "here" << endl;
// 	TGraphAsymmErrors* graphCombThermalGammaSpectrumStatErr0020Plot2 = NULL;
// 	if (graphCombThermalGammaSpectrumStatErr0020) {
// 		graphCombThermalGammaSpectrumStatErr0020Plot2 = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumStatErr0020->Clone("graphCombThermalGammaSpectrumStatErr0020Plot2");
// 		ProduceGraphAsymmWithoutXErrors(graphCombThermalGammaSpectrumStatErr0020Plot2);
// 		while (graphCombThermalGammaSpectrumStatErr0020Plot2->GetX()[graphCombThermalGammaSpectrumStatErr0020Plot2->GetN()-1] > 4.1 && graphCombThermalGammaSpectrumStatErr0020Plot2->GetN()>0)
// 			graphCombThermalGammaSpectrumStatErr0020Plot2->RemovePoint(graphCombThermalGammaSpectrumStatErr0020Plot2->GetN()-1);
// 		if (graphCombThermalGammaSpectrumStatErr0020Plot2->GetN()==0) graphCombThermalGammaSpectrumStatErr0020Plot2 = NULL;
// 	}	
// 	cout << "here" << endl;
	
	TH1D *dummDirGammaPHENIX0020 ;
	dummDirGammaPHENIX0020 = new TH1D("dummDirGammaPHENIX0020", "dummDirGammaPHENIX0020", 1000, 0., 5.6);
	SetStyleHistoTH1ForGraphs( dummDirGammaPHENIX0020, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
	dummDirGammaPHENIX0020->GetXaxis()->SetLabelOffset(-0.002);
	dummDirGammaPHENIX0020->GetXaxis()->SetTickLength(0.025);
	dummDirGammaPHENIX0020->GetYaxis()->SetTickLength(0.025);
	dummDirGammaPHENIX0020->GetYaxis()->SetRangeUser(0.9e-5,50);
	dummDirGammaPHENIX0020->GetXaxis()->SetRangeUser(0.,doubleRatioX[1]);
	dummDirGammaPHENIX0020->DrawCopy();
	
	dummDirGammaPHENIX0020->DrawCopy("same,axis");

	DrawGammaSetMarkerTGraphAsym(graphPHENIXAuAuSys0020, markerStylePHENIX, markerSizePHENIX, colorPHENIX , colorPHENIX, widthLinesBoxes,  kTRUE);		
	graphPHENIXAuAuSys0020->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPHENIXAuAuStat0020, markerStylePHENIX, markerSizePHENIX, colorPHENIX , colorPHENIX);		
	graphPHENIXAuAuStat0020->Draw("p,E1Z,same");
	histoFitThermalGamma0020PHENIXStat->Draw("same,l");

	if (histoFitThermalGamma0020Stat){
		SetStyleHisto(histoFitThermalGamma0020Stat, 3, styleFit, colorComb0020+1 );
		histoFitThermalGamma0020Stat->Scale(1e-2);
		histoFitThermalGamma0020Stat->Draw("same,l");
	}
	
	
	if (graphCombDirGammaSpectrumSystErr0020Plot3){
		graphCombDirGammaSpectrumSystErr0020Plot3= ScaleGraph(graphCombDirGammaSpectrumSystErr0020Plot3,1e-2);
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr0020Plot3, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr0020Plot3->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr0020Plot3){
		graphCombDirGammaSpectrumStatErr0020Plot3= ScaleGraph(graphCombDirGammaSpectrumStatErr0020Plot3,1e-2);
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr0020Plot3, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);	
		graphCombDirGammaSpectrumStatErr0020Plot3->Draw("p,E1Z,same");
	}
	
	
	TLegend* legendDirGammaWithPHENIX0020 = new TLegend(0.455,0.935-1.15*0.85*textsizeLabelsDirGamma*8,0.555+0.21,0.935);
	legendDirGammaWithPHENIX0020->SetFillStyle(0);
	legendDirGammaWithPHENIX0020->SetFillColor(0);
	legendDirGammaWithPHENIX0020->SetLineColor(0);
	legendDirGammaWithPHENIX0020->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaWithPHENIX0020->SetMargin(0.2);
	legendDirGammaWithPHENIX0020->SetTextFont(42);
	legendDirGammaWithPHENIX0020->AddEntry(graphCombDirGammaSpectrumSystErr0020Plot3,"ALICE","pf");
	legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,Form("%s",collisionSystemCent0020.Data()),"");
	legendDirGammaWithPHENIX0020->AddEntry(histoFitThermalGamma0020Stat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
// 	legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,"#it{T}_{eff} = 297 #pm 12^{#it{stat}} #pm 41^{#it{sys}}",""); // numbers for thermal
	legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,"#it{T}_{eff} = 304 #pm 11^{ #it{stat}} #pm 40^{#it{sys}} MeV","");
	legendDirGammaWithPHENIX0020->AddEntry(graphPHENIXAuAuSys0020,"PHENIX","pf");
	legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,Form("%s",collisionSystemRHIC0020.Data()),"");
	legendDirGammaWithPHENIX0020->AddEntry(histoFitThermalGamma0020PHENIXStat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
	legendDirGammaWithPHENIX0020->AddEntry((TObject*)0,"#it{T}_{eff} = 239 #pm 25^{#it{stat}} #pm 7^{ #it{sys}} MeV","");

	legendDirGammaWithPHENIX0020->Draw();
	
	canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusPHENIX0020_ReducedX.%s",outputDir.Data(),suffix.Data()));

	//*******************************************************************************************************
	//********************************PHENIX and ALICE together for 0-20% ***********************************
	//*******************************************************************************************************
// 	graphPHENIXAuAuStat2040->RemovePoint(4);
// 	graphPHENIXAuAuStat2040->RemovePoint(graphPHENIXAuAuStat2040->GetN()-1);
// 	graphPHENIXAuAuSys2040->RemovePoint(graphPHENIXAuAuSys2040->GetN()-1);
	TF1* fitThermalGamma2040PHENIXStat				= FitObject("e","fitThermalGamma2040PHENIXStat","Photon",graphPHENIXAuAuStat2040,0.6,2,NULL,"QNRMEX0+");
	fitThermalGamma2040PHENIXStat->FixParameter(1,0.260);
	fitThermalGamma2040PHENIXStat->SetParLimits(0,0,50);
	graphPHENIXAuAuStat2040->Fit(fitThermalGamma2040PHENIXStat,"QNRMEX0+","",0.6,2.);
	cout << WriteParameterToFile(fitThermalGamma2040PHENIXStat)<< endl;   

	TH1D* histoFitThermalGamma2040PHENIXStat		= (TH1D*)fitThermalGamma2040PHENIXStat->GetHistogram();
	SetStyleHisto(histoFitThermalGamma2040PHENIXStat, 3, styleFit, colorPHENIX );
	
	TGraphAsymmErrors* graphCombDirGammaSpectrumSystErr2040Plot3 = NULL;
	if (graphCombDirGammaSpectrumSystErr2040Plot) {
		graphCombDirGammaSpectrumSystErr2040Plot3 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSystErr2040Plot->Clone("graphCombDirGammaSpectrumSystErr2040Plot3");
		while (graphCombDirGammaSpectrumSystErr2040Plot3->GetX()[graphCombDirGammaSpectrumSystErr2040Plot3->GetN()-1] > 5.1 && graphCombDirGammaSpectrumSystErr2040Plot3->GetN()>0)
			graphCombDirGammaSpectrumSystErr2040Plot3->RemovePoint(graphCombDirGammaSpectrumSystErr2040Plot3->GetN()-1);
		if (graphCombDirGammaSpectrumSystErr2040Plot3->GetN()==0) graphCombDirGammaSpectrumSystErr2040Plot3 = NULL;
	}
	TGraphAsymmErrors* graphCombDirGammaSpectrumStatErr2040Plot3 = NULL;
	if (graphCombDirGammaSpectrumStatErr2040Plot) {
		graphCombDirGammaSpectrumStatErr2040Plot3 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumStatErr2040Plot->Clone("graphCombDirGammaSpectrumStatErr2040Plot3");
		while (graphCombDirGammaSpectrumStatErr2040Plot3->GetX()[graphCombDirGammaSpectrumStatErr2040Plot3->GetN()-1] > 5.1 && graphCombDirGammaSpectrumStatErr2040Plot3->GetN()>0)
			graphCombDirGammaSpectrumStatErr2040Plot3->RemovePoint(graphCombDirGammaSpectrumStatErr2040Plot3->GetN()-1);
		if (graphCombDirGammaSpectrumStatErr2040Plot3->GetN()==0) graphCombDirGammaSpectrumStatErr2040Plot3 = NULL;
	}	

	TGraphAsymmErrors* graphCombDirGammaSpectrumSumErr2040ArPlot3 = NULL;
	if (graphCombDirGammaSpectrumSumErr2040ArPlot) {
		graphCombDirGammaSpectrumSumErr2040ArPlot3 = (TGraphAsymmErrors*)graphCombDirGammaSpectrumSumErr2040ArPlot->Clone("graphCombDirGammaSpectrumSumErr2040ArPlot3");
		while (graphCombDirGammaSpectrumSumErr2040ArPlot3->GetX()[graphCombDirGammaSpectrumSumErr2040ArPlot3->GetN()-1] > 5.1 && graphCombDirGammaSpectrumSumErr2040ArPlot3->GetN()>0)
			graphCombDirGammaSpectrumSumErr2040ArPlot3->RemovePoint(graphCombDirGammaSpectrumSumErr2040ArPlot3->GetN()-1);
		if (graphCombDirGammaSpectrumSumErr2040ArPlot3->GetN()==0) graphCombDirGammaSpectrumSumErr2040ArPlot3 = NULL;
	}		
	
// 	cout << "here" << endl;
// 	TGraphAsymmErrors* graphCombThermalGammaSpectrumSystErr2040Plot2 = NULL;
// 	if (graphCombThermalGammaSpectrumSysErr2040) {
// 		graphCombThermalGammaSpectrumSystErr2040Plot2 = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumSysErr2040->Clone("graphCombThermalGammaSpectrumSystErr2040Plot2");
// 		while (graphCombThermalGammaSpectrumSystErr2040Plot2->GetX()[graphCombThermalGammaSpectrumSystErr2040Plot2->GetN()-1] > 4.1 && graphCombThermalGammaSpectrumSystErr2040Plot2->GetN()>0)
// 			graphCombThermalGammaSpectrumSystErr2040Plot2->RemovePoint(graphCombThermalGammaSpectrumSystErr2040Plot2->GetN()-1);
// 		if (graphCombThermalGammaSpectrumSystErr2040Plot2->GetN()==0) graphCombThermalGammaSpectrumSystErr2040Plot2 = NULL;
// 	}
// 	cout << "here" << endl;
// 	TGraphAsymmErrors* graphCombThermalGammaSpectrumStatErr2040Plot2 = NULL;
// 	if (graphCombThermalGammaSpectrumStatErr2040) {
// 		graphCombThermalGammaSpectrumStatErr2040Plot2 = (TGraphAsymmErrors*)graphCombThermalGammaSpectrumStatErr2040->Clone("graphCombThermalGammaSpectrumStatErr2040Plot2");
// 		ProduceGraphAsymmWithoutXErrors(graphCombThermalGammaSpectrumStatErr2040Plot2);
// 		while (graphCombThermalGammaSpectrumStatErr2040Plot2->GetX()[graphCombThermalGammaSpectrumStatErr2040Plot2->GetN()-1] > 4.1 && graphCombThermalGammaSpectrumStatErr2040Plot2->GetN()>0)
// 			graphCombThermalGammaSpectrumStatErr2040Plot2->RemovePoint(graphCombThermalGammaSpectrumStatErr2040Plot2->GetN()-1);
// 		if (graphCombThermalGammaSpectrumStatErr2040Plot2->GetN()==0) graphCombThermalGammaSpectrumStatErr2040Plot2 = NULL;
// 	}	
// 	cout << "here" << endl;
	canvasDirGamma->cd();
	TH2D *dummDirGammaPHENIX2040 ;
	dummDirGammaPHENIX2040 = new TH2D("dummDirGammaPHENIX2040", "dummDirGammaPHENIX2040", 100, 0., 5.6,100000,2e-5,20);
	SetStyleHistoTH2ForGraphs( dummDirGammaPHENIX2040, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
	dummDirGammaPHENIX2040->GetXaxis()->SetLabelOffset(-0.002);
	dummDirGammaPHENIX2040->GetXaxis()->SetTickLength(0.025);
	dummDirGammaPHENIX2040->GetYaxis()->SetTickLength(0.025);
	dummDirGammaPHENIX2040->GetYaxis()->SetRangeUser(1.1e-5,20);
	dummDirGammaPHENIX2040->GetXaxis()->SetRangeUser(0.,doubleRatioX[1]);
	dummDirGammaPHENIX2040->DrawCopy();

	DrawGammaSetMarkerTGraphAsym(graphPHENIXAuAuSys2040, markerStylePHENIX, markerSizePHENIX, colorPHENIX , colorPHENIX, widthLinesBoxes,  kTRUE);		
	graphPHENIXAuAuSys2040->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPHENIXAuAuStat2040, markerStylePHENIX, markerSizePHENIX, colorPHENIX , colorPHENIX);		
	graphPHENIXAuAuStat2040->Draw("p,E1Z,same");
	histoFitThermalGamma2040PHENIXStat->Draw("same,l");

	if (histoFitThermalGamma2040Stat){
		SetStyleHisto(histoFitThermalGamma2040Stat, 3, styleFit, colorComb2040+1 );
		histoFitThermalGamma2040Stat->Scale(1e-1);
		histoFitThermalGamma2040Stat->Draw("same,l");
	}
		
	if (graphCombDirGammaSpectrumSystErr2040Plot3){
		graphCombDirGammaSpectrumSystErr2040Plot3= ScaleGraph(graphCombDirGammaSpectrumSystErr2040Plot3,1e-1);
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSystErr2040Plot3, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
		graphCombDirGammaSpectrumSystErr2040Plot3->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr2040Plot3){
		graphCombDirGammaSpectrumStatErr2040Plot3= ScaleGraph(graphCombDirGammaSpectrumStatErr2040Plot3,1e-1);
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumStatErr2040Plot3, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
		graphCombDirGammaSpectrumStatErr2040Plot3->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr2040ArPlot3){
		graphCombDirGammaSpectrumSumErr2040ArPlot3= ScaleGraph(graphCombDirGammaSpectrumSumErr2040ArPlot3,1e-1);
		DrawGammaSetMarkerTGraphAsym(graphCombDirGammaSpectrumSumErr2040ArPlot3, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
		graphCombDirGammaSpectrumSumErr2040ArPlot3->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot3);
	}

	TLegend* legendDirGammaWithPHENIX2040 = new TLegend(0.455,0.935-1.15*0.85*textsizeLabelsDirGamma*8,0.555+0.21,0.935);
	legendDirGammaWithPHENIX2040->SetFillStyle(0);
	legendDirGammaWithPHENIX2040->SetFillColor(0);
	legendDirGammaWithPHENIX2040->SetLineColor(0);
	legendDirGammaWithPHENIX2040->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaWithPHENIX2040->SetMargin(0.2);
	legendDirGammaWithPHENIX2040->SetTextFont(42);
	legendDirGammaWithPHENIX2040->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot3,"ALICE","pf");
	legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,Form("%s",collisionSystemCent2040.Data()),"");
	legendDirGammaWithPHENIX2040->AddEntry(histoFitThermalGamma2040Stat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
	legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,"#it{T}_{eff} = 407 #pm 61^{ #it{stat}} #pm 96^{#it{sys}} MeV","");
	// 	legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,"#it{T}_{eff} = 410 #pm 84^{#it{stat}} #pm 140^{#it{sys}}",""); // numbers for thermal
	legendDirGammaWithPHENIX2040->AddEntry(graphPHENIXAuAuSys2040,"PHENIX","pf");
	legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,Form("%s",collisionSystemRHIC2040.Data()),"");
	legendDirGammaWithPHENIX2040->AddEntry(histoFitThermalGamma2040PHENIXStat,"#it{A} exp(-#it{p}_{T}/#it{T}_{eff})","l");
	legendDirGammaWithPHENIX2040->AddEntry((TObject*)0,"#it{T}_{eff} = 260 #pm 33^{#it{ stat}} #pm 8^{#it{sys}} MeV","");
	legendDirGammaWithPHENIX2040->Draw();

	canvasDirGamma->Print(Form("%s/DirGammaSpectrumPlusPHENIX2040_ReducedX.%s",outputDir.Data(),suffix.Data()));
	
	//*******************************************************************************************************************************************
	//*************************************************** Plotting direct Gamma Spectrum PCM only ***********************************************
	//*******************************************************************************************************************************************
	
	TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr0020Plot;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr0020Plot;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr0020ArPlot;
	if (graphPCMDirGammaStatErr0020) graphPCMDirGammaSpectrumStatErr0020Plot 		= ScaleGraph(graphPCMDirGammaStatErr0020,100);
	if (graphPCMDirGammaSysErr0020) graphPCMDirGammaSpectrumSystErr0020Plot 		= ScaleGraph(graphPCMDirGammaSysErr0020,100);
	if (graphPCMDirGammaSumErrAr0020) graphPCMDirGammaSpectrumSumErr0020ArPlot 	= ScaleGraph(graphPCMDirGammaSumErrAr0020,100);

	TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr2040Plot;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr2040Plot;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr2040ArPlot;
	if (graphPCMDirGammaStatErr2040) graphPCMDirGammaSpectrumStatErr2040Plot 		= ScaleGraph(graphPCMDirGammaStatErr2040,10);
	if (graphPCMDirGammaSysErr2040) graphPCMDirGammaSpectrumSystErr2040Plot 		= ScaleGraph(graphPCMDirGammaSysErr2040,10);
	if (graphPCMDirGammaSumErrAr2040) graphPCMDirGammaSpectrumSumErr2040ArPlot 	= ScaleGraph(graphPCMDirGammaSumErrAr2040,10);
	TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErr4080Plot;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErr4080Plot;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErr4080ArPlot;
	if (graphPCMDirGammaStatErr4080) graphPCMDirGammaSpectrumStatErr4080Plot 		= ScaleGraph(graphPCMDirGammaStatErr4080,1);
	if (graphPCMDirGammaSysErr4080) graphPCMDirGammaSpectrumSystErr4080Plot 		= ScaleGraph(graphPCMDirGammaSysErr4080,1);
	if (graphPCMDirGammaSumErrAr4080) graphPCMDirGammaSpectrumSumErr4080ArPlot 	= ScaleGraph(graphPCMDirGammaSumErrAr4080,1);
	
	canvasDirGamma->SetLogx(1);
	canvasDirGamma->cd();
	dummyDirGamma->DrawCopy();

	graphTheoryEPS090020Plot->Draw("p3lsame");
	graphTheoryCT100020Plot->Draw("p3lsame");
	graphTheoryNLO0020Plot->Draw("p3lsame");
	graphTheoryPromptMcGill0020Plot->Draw("p3lsame");
	
	if (graphPCMDirGammaSpectrumSystErr0020Plot){
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
		graphPCMDirGammaSpectrumSystErr0020Plot->Draw("E2same");
	}
	if (graphPCMDirGammaSpectrumStatErr0020Plot){
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);	
		graphPCMDirGammaSpectrumStatErr0020Plot->Draw("p,E1Z,same");
	}
// 	if (graphPCMDirGammaSpectrumSumErr0020ArPlot){
// 		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr0020ArPlot , 1, 3, colorComb0020, colorComb0020, 1.8, kTRUE);
// 		graphPCMDirGammaSpectrumSumErr0020ArPlot->Draw(">,same");
// 	}	

	graphTheoryEPS092040Plot->Draw("p3lsame");
	graphTheoryCT102040Plot->Draw("p3lsame");
	graphTheoryNLO2040Plot->Draw("p3lsame");
	graphTheoryPromptMcGill2040Plot->Draw("p3lsame");

	if (graphPCMDirGammaSpectrumSystErr2040Plot){
		graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
		graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
		graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
		graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
		graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
		graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
		graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(0);
		graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(5);
		graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(6);
		graphPCMDirGammaSpectrumSystErr2040Plot->RemovePoint(6);
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
		graphPCMDirGammaSpectrumSystErr2040Plot->Draw("E2same");
	}	
	if (graphPCMDirGammaSpectrumStatErr2040Plot){
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
		graphPCMDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
	}
	if (graphPCMDirGammaSpectrumSumErr2040ArPlot){
		graphPCMDirGammaSpectrumSumErr2040ArPlot->RemovePoint(0);
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr2040ArPlot , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
		graphPCMDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr2040ArPlot);
	}	
	graphTheoryEPS094080Plot->Draw("p3lsame");
	graphTheoryCT104080Plot->Draw("p3lsame");
	graphTheoryNLO4080Plot->Draw("p3lsame");

// 	if (graphPCMDirGammaSpectrumSystErr4080Plot){
// 		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
// 		graphPCMDirGammaSpectrumSystErr4080Plot->Draw("E2same");
// 	}
	if (graphPCMDirGammaSpectrumStatErr4080Plot){
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);	
		graphPCMDirGammaSpectrumStatErr4080Plot->Draw("p,E1Z,same");
	}
	if (graphPCMDirGammaSpectrumSumErr4080ArPlot){
		
		graphPCMDirGammaSpectrumSumErr4080ArPlot->RemovePoint(0);
		graphPCMDirGammaSpectrumSumErr4080ArPlot->Print();
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErr4080ArPlot , 1, 3, colorComb4080, colorComb4080, 1.8, kTRUE);
		graphPCMDirGammaSpectrumSumErr4080ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErr4080ArPlot);
	}	
	
	
	labelScalingDirGamma0020->Draw();	
	labelScalingDirGamma2040->Draw();
	labelScalingDirGamma4080->Draw();
	labelDirGammaColl->Draw();

	legendDirGamma->Draw();
	legendDirGammaTheory->Draw();
	legendDirGammaTheory2->Draw();
	
	canvasDirGamma->Print(Form("%s/DirGammaSpectrum_PCMonly.%s",outputDir.Data(),suffix.Data()));

	
	//*******************************************************************************************************************************************
	//*************************************************** Plotting direct Gamma Spectrum PCM only ***********************************************
	//*******************************************************************************************************************************************
	
	TGraphAsymmErrors* graphPHOSDirGammaSpectrumStatErr0020Plot;
	TGraphAsymmErrors* graphPHOSDirGammaSpectrumSystErr0020Plot;
	if (histoPHOSDirGammaStatErr0020){
		graphPHOSDirGammaSpectrumStatErr0020Plot 		= new TGraphAsymmErrors(histoPHOSDirGammaStatErr0020);	
		graphPHOSDirGammaSpectrumStatErr0020Plot 		= ScaleGraph(graphPHOSDirGammaSpectrumStatErr0020Plot,100);
	}
	if (histoPHOSDirGammaSysErr0020){
		graphPHOSDirGammaSpectrumSystErr0020Plot 		= new TGraphAsymmErrors(histoPHOSDirGammaSysErr0020);	
		graphPHOSDirGammaSpectrumSystErr0020Plot 		= ScaleGraph(graphPHOSDirGammaSpectrumSystErr0020Plot,100);
	}
	TGraphAsymmErrors* graphPHOSDirGammaSpectrumStatErr2040Plot;
	TGraphAsymmErrors* graphPHOSDirGammaSpectrumSystErr2040Plot;
	if (histoPHOSDirGammaStatErr2040){
		graphPHOSDirGammaSpectrumStatErr2040Plot 		= new TGraphAsymmErrors(histoPHOSDirGammaStatErr2040);	
		graphPHOSDirGammaSpectrumStatErr2040Plot 		= ScaleGraph(graphPHOSDirGammaSpectrumStatErr2040Plot,10);
	}
	if (histoPHOSDirGammaSysErr2040){
		graphPHOSDirGammaSpectrumSystErr2040Plot 		= new TGraphAsymmErrors(histoPHOSDirGammaSysErr2040);	
		graphPHOSDirGammaSpectrumSystErr2040Plot 		= ScaleGraph(graphPHOSDirGammaSpectrumSystErr2040Plot,10);
	}

	TGraphAsymmErrors* graphPHOSDirGammaSpectrumStatErr4080Plot;
	TGraphAsymmErrors* graphPHOSDirGammaSpectrumSystErr4080Plot;
	if (histoPHOSDirGammaStatErr4080){
		graphPHOSDirGammaSpectrumStatErr4080Plot 		= new TGraphAsymmErrors(histoPHOSDirGammaStatErr4080);	
		graphPHOSDirGammaSpectrumStatErr4080Plot 		= ScaleGraph(graphPHOSDirGammaSpectrumStatErr4080Plot,1);
	}	
	if (histoPHOSDirGammaSysErr4080){
		graphPHOSDirGammaSpectrumSystErr4080Plot 		= new TGraphAsymmErrors(histoPHOSDirGammaSysErr4080);	
		graphPHOSDirGammaSpectrumSystErr4080Plot 		= ScaleGraph(graphPHOSDirGammaSpectrumSystErr4080Plot,1);
	}	
	
	canvasDirGamma->SetLogx(1);
	canvasDirGamma->cd();
	dummyDirGamma->DrawCopy();

	graphTheoryEPS090020Plot->Draw("p3lsame");
	graphTheoryCT100020Plot->Draw("p3lsame");
	graphTheoryNLO0020Plot->Draw("p3lsame");
	graphTheoryPromptMcGill0020Plot->Draw("p3lsame");
	
	if (graphPHOSDirGammaSpectrumSystErr0020Plot){
		DrawGammaSetMarkerTGraphAsym(graphPHOSDirGammaSpectrumSystErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
		graphPHOSDirGammaSpectrumSystErr0020Plot->Draw("E2same");
	}
	if (graphPHOSDirGammaSpectrumStatErr0020Plot){
		DrawGammaSetMarkerTGraphAsym(graphPHOSDirGammaSpectrumStatErr0020Plot, markerStyleComb0020, markerSizeComb0020, colorComb0020 , colorComb0020);	
		graphPHOSDirGammaSpectrumStatErr0020Plot->Draw("p,E1Z,same");
	}

	graphTheoryEPS092040Plot->Draw("p3lsame");
	graphTheoryCT102040Plot->Draw("p3lsame");
	graphTheoryNLO2040Plot->Draw("p3lsame");
	graphTheoryPromptMcGill2040Plot->Draw("p3lsame");

	if (graphPHOSDirGammaSpectrumSystErr2040Plot){
		DrawGammaSetMarkerTGraphAsym(graphPHOSDirGammaSpectrumSystErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
		graphPHOSDirGammaSpectrumSystErr2040Plot->Draw("E2same");
	}	
	if (graphPHOSDirGammaSpectrumStatErr2040Plot){
		DrawGammaSetMarkerTGraphAsym(graphPHOSDirGammaSpectrumStatErr2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
		graphPHOSDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
	}
	graphTheoryEPS094080Plot->Draw("p3lsame");
	graphTheoryCT104080Plot->Draw("p3lsame");
	graphTheoryNLO4080Plot->Draw("p3lsame");

	if (graphPHOSDirGammaSpectrumSystErr4080Plot){
		DrawGammaSetMarkerTGraphAsym(graphPHOSDirGammaSpectrumSystErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
		graphPHOSDirGammaSpectrumSystErr4080Plot->Draw("E2same");
	}
	if (graphPHOSDirGammaSpectrumStatErr4080Plot){
		DrawGammaSetMarkerTGraphAsym(graphPHOSDirGammaSpectrumStatErr4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);	
		graphPHOSDirGammaSpectrumStatErr4080Plot->Draw("p,E1Z,same");
	}
	
	labelScalingDirGamma0020->Draw();	
	labelScalingDirGamma2040->Draw();
	labelScalingDirGamma4080->Draw();
	labelDirGammaColl->Draw();

	legendDirGamma->Draw();
	legendDirGammaTheory->Draw();
	legendDirGammaTheory2->Draw();
	
	canvasDirGamma->Print(Form("%s/DirGammaSpectrum_PHOSonly.%s",outputDir.Data(),suffix.Data()));
	
	//*******************************************************************************************************************************************
	//*************************************************** Plotting direct Gamma Spectrum PCM only ***********************************************
	//*******************************************************************************************************************************************
	
	Double_t xSection2760GeVppINEL 	= 62.8;
	Double_t xSection5023GeVppINEL 	= 70;
	Double_t xSection7TeVppINEL    	= 73.2;
	Double_t nCollpPb				= 6.9; 
	TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErrpPbPlot 					= NULL;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErrpPbPlot					= NULL;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErrpPbArPlot					= NULL;
	if (graphPCMDirGammaStatErrpPb) graphPCMDirGammaSpectrumStatErrpPbPlot 		= ScaleGraph(graphPCMDirGammaStatErrpPb,10000);
	if (graphPCMDirGammaSysErrpPb) graphPCMDirGammaSpectrumSystErrpPbPlot 		= ScaleGraph(graphPCMDirGammaSysErrpPb,10000);
	if (graphPCMDirGammaSumErrArpPb) graphPCMDirGammaSpectrumSumErrpPbArPlot 	= ScaleGraph(graphPCMDirGammaSumErrArpPb,10000);
	TGraphAsymmErrors* graphTheoryNLOpPbPlot 									= ScaleGraph(graphTheoryNLOpPb,10000);
	
	TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErrpp7TeVPlot					= NULL;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErrpp7TeVPlot					= NULL;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErrpp7TeVArPlot					= NULL;
	if (graphPCMDirGammaStatErrpp7TeV) graphPCMDirGammaSpectrumStatErrpp7TeVPlot 	= ScaleGraph(graphPCMDirGammaStatErrpp7TeV,100);
	if (graphPCMDirGammaSysErrpp7TeV) graphPCMDirGammaSpectrumSystErrpp7TeVPlot 	= ScaleGraph(graphPCMDirGammaSysErrpp7TeV,100);
	if (graphPCMDirGammaSumErrArpp7TeV) graphPCMDirGammaSpectrumSumErrpp7TeVArPlot 	= ScaleGraph(graphPCMDirGammaSumErrArpp7TeV,100);
	TGraphAsymmErrors* graphTheoryNLOpp7TeVPlot 									= ScaleGraph(graphTheoryNLOpp7TeV,100);
	
	TGraphAsymmErrors* graphPCMDirGammaSpectrumStatErrpp2760GeVPlot;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSystErrpp2760GeVPlot;
	TGraphAsymmErrors* graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot;
	if (graphPCMDirGammaStatErrpp2760GeV) graphPCMDirGammaSpectrumStatErrpp2760GeVPlot 		= ScaleGraph(graphPCMDirGammaStatErrpp2760GeV,1);
	if (graphPCMDirGammaSysErrpp2760GeV) graphPCMDirGammaSpectrumSystErrpp2760GeVPlot 		= ScaleGraph(graphPCMDirGammaSysErrpp2760GeV,1);
	if (graphPCMDirGammaSumErrArpp2760GeV) graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot 	= ScaleGraph(graphPCMDirGammaSumErrArpp2760GeV,1);
	TGraphAsymmErrors* graphTheoryNLOpp2760GeVPlot 											= ScaleGraph(graphTheoryNLOpp2760GeV,1);
		
	TF1* fitNLOMcGill2760GeV = new TF1 ("JFNLOpp2760GeV","TMath::Exp([0] + [1]*TMath::Log(x) + [2]*TMath::Log(x)*TMath::Log(x))*[3]*1./([4]/TMath::Power([6],[5]))",0.5,10);  //
	fitNLOMcGill2760GeV->SetParameter(0,16.586869134650662);
	fitNLOMcGill2760GeV->SetParameter(1,-4.159632883579499);
	fitNLOMcGill2760GeV->SetParameter(2,-0.2320496354057888);
	fitNLOMcGill2760GeV->SetParameter(3,1e-9);
	fitNLOMcGill2760GeV->SetParameter(4,0.865779);
	fitNLOMcGill2760GeV->SetParameter(5,0.0694875);
	fitNLOMcGill2760GeV->SetParameter(6,4);
	TH1D* histoTheoryNLOMcGill2760GeV = (TH1D*)fitNLOMcGill2760GeV->GetHistogram();
	histoTheoryNLOMcGill2760GeV->Scale(1./xSection2760GeVppINEL);

	TF1* fitNLOMcGill5023GeV = new TF1 ("JFNLOpp5023GeV","TMath::Exp([0] + [1]*TMath::Log(x) + [2]*TMath::Log(x)*TMath::Log(x))*[3]*1./([4]/TMath::Power([6],[5]))",0.5,10);  //
	fitNLOMcGill5023GeV->SetParameter(0,16.870532033687542);
	fitNLOMcGill5023GeV->SetParameter(1,-4.05339420384128);
	fitNLOMcGill5023GeV->SetParameter(2,-0.20936877337479987);
	fitNLOMcGill5023GeV->SetParameter(3,1e-9);
	fitNLOMcGill5023GeV->SetParameter(4,0.905908);
	fitNLOMcGill5023GeV->SetParameter(5,0.0577656);
	fitNLOMcGill5023GeV->SetParameter(6,8);
	TH1D* histoTheoryNLOMcGill5023GeV = (TH1D*)fitNLOMcGill5023GeV->GetHistogram();
	histoTheoryNLOMcGill5023GeV->Scale(1./xSection5023GeVppINEL*10000*nCollpPb);
	
	TF1* fitNLOMcGill7TeV = new TF1 ("JFNLOpp7TeV","TMath::Exp([0] + [1]*TMath::Log(x) + [2]*TMath::Log(x)*TMath::Log(x))*[3]*1./([4]/TMath::Power([6],[5]))",0.5,10);  //
	fitNLOMcGill7TeV->SetParameter(0,17.0572573576821);
	fitNLOMcGill7TeV->SetParameter(1,-3.971862631159262);
	fitNLOMcGill7TeV->SetParameter(2,-0.2213356968221152);
	fitNLOMcGill7TeV->SetParameter(3,1e-9);
	fitNLOMcGill7TeV->SetParameter(4,0.890255);
	fitNLOMcGill7TeV->SetParameter(5,0.0541605);
	fitNLOMcGill7TeV->SetParameter(6,4);
	TH1D* histoTheoryNLOMcGill7TeV = (TH1D*)fitNLOMcGill7TeV->GetHistogram();
	histoTheoryNLOMcGill7TeV->Scale(1./xSection7TeVppINEL*100);
	
	
	canvasDirGamma->cd();
	TH2D* dummyDirGammapp;
	dummyDirGammapp = new TH2D("dummyDirGammapp", "dummyDirGammapp", 1000, 0., 25, 1000., 6e-11,4e6);
	SetStyleHistoTH2ForGraphs( dummyDirGammapp, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
	dummyDirGammapp->GetXaxis()->SetLabelOffset(-0.015);
	dummyDirGammapp->GetXaxis()->SetRangeUser(doubleRatioXpp[0],doubleRatioXpp[1]);
// 	dummyDirGamma->GetXaxis()->SetRangeUser(0,16);
	dummyDirGammapp->DrawCopy();

	SetStyleGammaNLOTGraphWithBand( graphTheoryNLOpPbPlot, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);
	graphTheoryNLOpPbPlot->Draw("p3lsame");
	
	if (graphPCMDirGammaSpectrumSystErrpPbPlot){
		graphPCMDirGammaSpectrumSystErrpPbPlot->RemovePoint(1);
		graphPCMDirGammaSpectrumSystErrpPbPlot->RemovePoint(1);
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErrpPbPlot, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb, widthLinesBoxes, kTRUE);
		graphPCMDirGammaSpectrumSystErrpPbPlot->Draw("E2same");
		
	}
	if (graphPCMDirGammaSpectrumStatErrpPbPlot){
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErrpPbPlot, markerStyleCombpPb, markerSizeCombpPb, colorCombpPb , colorCombpPb);	
		graphPCMDirGammaSpectrumStatErrpPbPlot->Draw("p,E1Z,same");
	}
	if (graphPCMDirGammaSpectrumSumErrpPbArPlot){
		graphPCMDirGammaSpectrumSumErrpPbArPlot->RemovePoint(0);
		graphPCMDirGammaSpectrumSumErrpPbArPlot->RemovePoint(graphPCMDirGammaSpectrumSumErrpPbArPlot->GetN()-3);
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErrpPbArPlot , 1, 3, colorCombpPb, colorCombpPb, 1.8, kTRUE);
		graphPCMDirGammaSpectrumSumErrpPbArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErrpPbArPlot);		
	}	
	SetStyleHisto(histoTheoryNLOMcGill5023GeV, widthCommonFit, styleNLOMcGill, colorNLOMcGill);
	histoTheoryNLOMcGill5023GeV->Draw("same,c");

	
	SetStyleGammaNLOTGraphWithBand( graphTheoryNLOpp7TeVPlot, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);
	graphTheoryNLOpp7TeVPlot->Draw("p3lsame");

	if (graphPCMDirGammaSpectrumSystErrpp7TeVPlot){
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErrpp7TeVPlot, markerStyleCombpp7TeV, markerSizeCombpp7TeV, colorCombpp7TeV , colorCombpp7TeV, widthLinesBoxes, kTRUE);
		graphPCMDirGammaSpectrumSystErrpp7TeVPlot->Draw("E2same");
	}	
	if (graphPCMDirGammaSpectrumStatErrpp7TeVPlot){
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErrpp7TeVPlot, markerStyleCombpp7TeV, markerSizeCombpp7TeV, colorCombpp7TeV , colorCombpp7TeV);	
		graphPCMDirGammaSpectrumStatErrpp7TeVPlot->Draw("p,E1Z,same");
	}
	if (graphPCMDirGammaSpectrumSumErrpp7TeVArPlot){
		graphPCMDirGammaSpectrumSumErrpp7TeVArPlot->RemovePoint(0);
		graphPCMDirGammaSpectrumSumErrpp7TeVArPlot->RemovePoint(graphPCMDirGammaSpectrumSumErrpp7TeVArPlot->GetN()-3);
		graphPCMDirGammaSpectrumSumErrpp7TeVArPlot->RemovePoint(graphPCMDirGammaSpectrumSumErrpp7TeVArPlot->GetN()-3);
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErrpp7TeVArPlot , 1, 3, colorCombpp7TeV, colorCombpp7TeV, 1.8, kTRUE);
		graphPCMDirGammaSpectrumSumErrpp7TeVArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErrpp7TeVArPlot);
	}	
	SetStyleHisto(histoTheoryNLOMcGill7TeV, widthCommonFit, styleNLOMcGill, colorNLOMcGill);
	histoTheoryNLOMcGill7TeV->Draw("same,c");

	SetStyleGammaNLOTGraphWithBand( graphTheoryNLOpp2760GeVPlot, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);
	graphTheoryNLOpp2760GeVPlot->Draw("p3lsame");
	
	
	if (graphPCMDirGammaSpectrumSystErrpp2760GeVPlot){
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSystErrpp2760GeVPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV, widthLinesBoxes, kTRUE);
		graphPCMDirGammaSpectrumSystErrpp2760GeVPlot->Draw("E2same");
	}
	if (graphPCMDirGammaSpectrumStatErrpp2760GeVPlot){
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumStatErrpp2760GeVPlot, markerStyleCombpp2760GeV, markerSizeCombpp2760GeV, colorCombpp2760GeV , colorCombpp2760GeV);	
		graphPCMDirGammaSpectrumStatErrpp2760GeVPlot->Draw("p,E1Z,same");
	}
	if (graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot){
		
		graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot->RemovePoint(0);
		graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot->RemovePoint(1);
		graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot->RemovePoint(graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot->GetN()-1);
		graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot->Print();
		DrawGammaSetMarkerTGraphAsym(graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot , 1, 3, colorCombpp2760GeV, colorCombpp2760GeV, 1.8, kTRUE);
		graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphPCMDirGammaSpectrumSumErrpp2760GeVArPlot);
	}	
	SetStyleHisto(histoTheoryNLOMcGill2760GeV, widthCommonFit, styleNLOMcGill, colorNLOMcGill);
	histoTheoryNLOMcGill2760GeV->Draw("same,c");
	
	TLatex *labelScalingDirGammapPb = new TLatex(13.5,1.1E-4,"x 10^{5}");
	SetStyleTLatex( labelScalingDirGammapPb, 0.85*textsizeLabelsDirGamma,4,colorCombpPb,42,kFALSE);
	labelScalingDirGammapPb->Draw();
	
	TLatex *labelScalingDirGammapp7TeV = new TLatex(13.5,9.5E-7,"x 10^{3}");
	SetStyleTLatex( labelScalingDirGammapp7TeV, 0.85*textsizeLabelsDirGamma,4,colorCombpp7TeV,42,kFALSE);
	labelScalingDirGammapp7TeV->Draw();
	
	TLatex *labelScalingDirGammapp2760GeV = new TLatex(13.5,5E-9,"x 10^{0}");
	SetStyleTLatex( labelScalingDirGammapp2760GeV, 0.85*textsizeLabelsDirGamma,4,colorCombpp2760GeV,42,kFALSE);
	labelScalingDirGammapp2760GeV->Draw();
// 	labelDirGammaColl->Draw();

	TLegend* legendDirGammapp = new TLegend(0.59,0.95-1.1*0.85*textsizeLabelsDirGamma*3,0.59+0.21,0.95);
	legendDirGammapp->SetFillStyle(0);
	legendDirGammapp->SetFillColor(0);
	legendDirGammapp->SetLineColor(0);
	legendDirGammapp->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammapp->SetMargin(0.25);
	legendDirGammapp->SetTextFont(42);
	legendDirGammapp->AddEntry(graphPCMDirGammaSpectrumSystErrpPbPlot,collisionSystempPb.Data(),"pf");
	legendDirGammapp->AddEntry(graphPCMDirGammaSpectrumSystErrpp7TeVPlot,collisionSystempp7TeV.Data(),"pf");
	legendDirGammapp->AddEntry(graphPCMDirGammaSpectrumSystErrpp2760GeVPlot,collisionSystempp2760GeV.Data(),"pf");
	legendDirGammapp->Draw();

	TLegend* legendDirGammaTheorypp = new TLegend(0.21,0.19-1.1*0.85*textsizeLabelsDirGamma*2,0.21+0.21,0.19);
	legendDirGammaTheorypp->SetFillStyle(0);
	legendDirGammaTheorypp->SetFillColor(0);
	legendDirGammaTheorypp->SetLineColor(0);
	legendDirGammaTheorypp->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaTheorypp->SetMargin(0.25);
	legendDirGammaTheorypp->SetTextFont(42);
	legendDirGammaTheorypp->AddEntry(graphTheoryNLOpPb,"PDF: CTEQ6M5 FF: GRV ","l");
	legendDirGammaTheorypp->AddEntry(histoTheoryNLOMcGill2760GeV,"nPDF: EPS09, FF: BFG2","l");
	legendDirGammaTheorypp->Draw();
	
	canvasDirGamma->Print(Form("%s/DirGammaSpectrum_PCMonly_ppAndpPb.%s",outputDir.Data(),suffix.Data()));
	
	
	//*******************************************************************************************************************************************
	//********************************************** Plotting direct Gamma Spectrum LinX ********************************************************
	//*******************************************************************************************************************************************
		
	TCanvas *canvasDirGammaLinX = new TCanvas("canvasDirGammaLinX","",10,10,1200,1400);  // gives the page size		
	DrawGammaCanvasSettings( canvasDirGammaLinX, 0.165, 0.01, 0.01, 0.07);
	canvasDirGammaLinX->SetLogy();
		
	TH2D *dummyDirGammaLinX ;
	dummyDirGammaLinX = new TH2D("dummyDirGammaLinX", "dummyDirGammaLinX", 1000, 0., 22, 1000., 6e-9,9e3);
	SetStyleHistoTH2ForGraphs( dummyDirGammaLinX, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}_{#gamma_{dir}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.85*textsizeLabelsDirGamma, textsizeLabelsDirGamma, 0.75, 1.65);
	dummyDirGammaLinX->GetXaxis()->SetLabelOffset(-0.002);
	dummyDirGammaLinX->GetXaxis()->SetRangeUser(0,14.5);
	dummyDirGammaLinX->DrawCopy();

	if (graphCombDirGammaSpectrumSystErr0020Plot){
		graphCombDirGammaSpectrumSystErr0020Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr0020Plot){
		graphCombDirGammaSpectrumStatErr0020Plot->Draw("p,E1Z,same");
	}

	if (graphCombDirGammaSpectrumSystErr2040Plot){
		graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
	}	
	if (graphCombDirGammaSpectrumStatErr2040Plot){
		graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr2040ArPlot){
		graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
	}	

	if (graphCombDirGammaSpectrumSystErr4080Plot){
		graphCombDirGammaSpectrumSystErr4080Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr4080Plot){
		graphCombDirGammaSpectrumStatErr4080Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr4080ArPlot){
		graphCombDirGammaSpectrumSumErr4080ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr4080ArPlot);
	}	

	TLatex *labelScalingDirGamma0020_2 = new TLatex(12.5,6.E-4,"x 10^{2}");
	SetStyleTLatex( labelScalingDirGamma0020_2, 0.85*textsizeLabelsDirGamma,4,colorComb0020,42,kFALSE);
	labelScalingDirGamma0020_2->Draw();
	
	TLatex *labelScalingDirGamma2040_2 = new TLatex(12.5,3.5E-5,"x 10^{1}");
	SetStyleTLatex( labelScalingDirGamma2040_2, 0.85*textsizeLabelsDirGamma,4,colorComb2040,42,kFALSE);
	labelScalingDirGamma2040_2->Draw();
	
	TLatex *labelScalingDirGamma4080_2 = new TLatex(12.5,6.5E-7,"x 10^{0}");
	SetStyleTLatex( labelScalingDirGamma4080_2, 0.85*textsizeLabelsDirGamma,4,colorComb4080,42,kFALSE);
	labelScalingDirGamma4080_2->Draw();
	
	TLatex *labelDirGammaColl_2 = new TLatex(0.26,0.94,Form("%s",collisionSystem.Data()));
	SetStyleTLatex( labelDirGammaColl_2, 0.85*textsizeLabelsDirGamma,4);
	labelDirGammaColl_2->Draw();

	TLegend* legendDirGammaLinX = new TLegend(0.26,0.93-1.1*0.83*textsizeLabelsDirGamma*3,0.26+0.21,0.93);
	legendDirGammaLinX->SetFillStyle(0);
	legendDirGammaLinX->SetFillColor(0);
	legendDirGammaLinX->SetLineColor(0);
	legendDirGammaLinX->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaLinX->SetMargin(0.25);
	legendDirGammaLinX->SetTextFont(42);
	legendDirGammaLinX->AddEntry(graphCombDirGammaSpectrumSystErr0020Plot,"  0-20% ALICE","pf");
	legendDirGammaLinX->AddEntry(graphCombDirGammaSpectrumSystErr2040Plot,"20-40% ALICE","pf");
	legendDirGammaLinX->AddEntry(graphCombDirGammaSpectrumSystErr4080Plot,"40-80% ALICE","pf");
	legendDirGammaLinX->Draw();
	
	canvasDirGammaLinX->Print(Form("%s/DirGammaSpectrum_LinX_withoutFit.%s",outputDir.Data(),suffix.Data()));
	
	graphTheoryEPS090020Plot->Draw("p3lsame");
	graphTheoryCT100020Plot->Draw("p3lsame");
	graphTheoryNLO0020Plot->Draw("p3lsame");
	graphTheoryPromptMcGill0020Plot->Draw("p3lsame");
	if (graphCombDirGammaSpectrumSystErr0020Plot){
		graphCombDirGammaSpectrumSystErr0020Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr0020Plot){
		graphCombDirGammaSpectrumStatErr0020Plot->Draw("p,E1Z,same");
	}
	
	graphTheoryEPS092040Plot->Draw("p3lsame");
	graphTheoryCT102040Plot->Draw("p3lsame");
	graphTheoryNLO2040Plot->Draw("p3lsame");
	graphTheoryPromptMcGill2040Plot->Draw("p3lsame");
	if (graphCombDirGammaSpectrumSystErr2040Plot){
		graphCombDirGammaSpectrumSystErr2040Plot->Draw("E2same");
	}	
	if (graphCombDirGammaSpectrumStatErr2040Plot){
		graphCombDirGammaSpectrumStatErr2040Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr2040ArPlot){
		graphCombDirGammaSpectrumSumErr2040ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr2040ArPlot);
	}	
	
	graphTheoryEPS094080Plot->Draw("p3lsame");
	graphTheoryCT104080Plot->Draw("p3lsame");
	graphTheoryNLO4080Plot->Draw("p3lsame");

	if (graphCombDirGammaSpectrumSystErr4080Plot){
		graphCombDirGammaSpectrumSystErr4080Plot->Draw("E2same");
	}
	if (graphCombDirGammaSpectrumStatErr4080Plot){
		graphCombDirGammaSpectrumStatErr4080Plot->Draw("p,E1Z,same");
	}
	if (graphCombDirGammaSpectrumSumErr4080ArPlot){
		graphCombDirGammaSpectrumSumErr4080ArPlot->Draw(">,same");
		PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombDirGammaSpectrumSumErr4080ArPlot);
	}	

	TLegend* legendDirGammaTheoryLinX = new TLegend(0.53,0.93-1.*0.83*textsizeLabelsDirGamma*3,0.53+0.21,0.93);
	legendDirGammaTheoryLinX->SetFillStyle(0);
	legendDirGammaTheoryLinX->SetFillColor(0);
	legendDirGammaTheoryLinX->SetLineColor(0);
	legendDirGammaTheoryLinX->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaTheoryLinX->SetMargin(0.23);
	legendDirGammaTheoryLinX->SetTextFont(42);
	legendDirGammaTheoryLinX->AddEntry(graphTheoryNLO0020Plot,"PDF: CTEQ6M5, FF: GRV ","l");
	legendDirGammaTheoryLinX->AddEntry(graphTheoryPromptMcGill2040Plot,"(n)PDF: CTEQ6.1M/EPS09,","l");
	legendDirGammaTheoryLinX->AddEntry((TObject*)0,"FF: BFG2","");

	legendDirGammaTheoryLinX->Draw();

	TLegend* legendDirGammaTheoryLinX2 = new TLegend(0.53,0.93-1.*0.83*textsizeLabelsDirGamma*7-0.02,0.53+0.21,0.93-1.*0.85*textsizeLabelsDirGamma*3-0.02);
	legendDirGammaTheoryLinX2->SetFillStyle(0);
	legendDirGammaTheoryLinX2->SetFillColor(0);
	legendDirGammaTheoryLinX2->SetLineColor(0);
	legendDirGammaTheoryLinX2->SetTextSize(0.85*textsizeLabelsDirGamma);
	legendDirGammaTheoryLinX2->SetMargin(0.23);
	legendDirGammaTheoryLinX2->SetTextFont(42);
	legendDirGammaTheoryLinX2->AddEntry((TObject*)0,"#it{JETPHOX}","");
	legendDirGammaTheoryLinX2->AddEntry(graphTheoryCT100020Plot,"PDF: CT10, FF: BFG2","f");
	legendDirGammaTheoryLinX2->AddEntry(graphTheoryEPS090020Plot,"nPDF: EPS09, FF: BFG2","f");
	legendDirGammaTheoryLinX2->AddEntry((TObject*)0,"(all scaled by #it{N}_{coll})","");
	legendDirGammaTheoryLinX2->Draw();
		
	canvasDirGammaLinX->Print(Form("%s/DirGammaSpectrum_LinX_withNLO_withoutFit.%s",outputDir.Data(),suffix.Data()));

	
	// ************************************************************************************************************************************************
	// ************************************************************ Plot Cocktail contributions *******************************************************
	// ************************************************************************************************************************************************
	histoCocktailSumGamma->Rebin(4);
	histoCocktailPi0Gamma->Rebin(4);
	histoCocktailEtaGamma->Rebin(4);
	histoCocktailOmegaGamma->Rebin(4);
	histoCocktailEtaPGamma->Rebin(4);
	histoCocktailRhoGamma->Rebin(4);
	histoCocktailPhiGamma->Rebin(4);
	histoCocktailSigmaGamma->Rebin(4);
	
	TH1D* histoCocktailRatioPi0GammaSumGamma = (TH1D*)histoCocktailPi0Gamma->Clone("histoCocktailRatioPi0GammaSumGamma");
	histoCocktailRatioPi0GammaSumGamma->Divide(histoCocktailRatioPi0GammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// 	histoCocktailRatioPi0GammaSumGamma->Scale(100);
	TH1D* histoCocktailRatioEtaGammaSumGamma = (TH1D*)histoCocktailEtaGamma->Clone("histoCocktailRatioEtaGammaSumGamma");
	histoCocktailRatioEtaGammaSumGamma->Divide(histoCocktailRatioEtaGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// 	histoCocktailRatioEtaGammaSumGamma->Scale(100);
	TH1D* histoCocktailRatioOmegaGammaSumGamma = (TH1D*)histoCocktailOmegaGamma->Clone("histoCocktailRatioOmegaGammaSumGamma");
	histoCocktailRatioOmegaGammaSumGamma->Divide(histoCocktailRatioOmegaGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// 	histoCocktailRatioOmegaGammaSumGamma->Scale(100);
	TH1D* histoCocktailRatioPhiGammaSumGamma = (TH1D*)histoCocktailPhiGamma->Clone("histoCocktailRatioPhiGammaSumGamma");
	histoCocktailRatioPhiGammaSumGamma->Divide(histoCocktailRatioPhiGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// 	histoCocktailRatioPhiGammaSumGamma->Scale(100);
	TH1D* histoCocktailRatioEtaPGammaSumGamma = (TH1D*)histoCocktailEtaPGamma->Clone("histoCocktailRatioEtaPGammaSumGamma");
	histoCocktailRatioEtaPGammaSumGamma->Divide(histoCocktailRatioEtaPGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// 	histoCocktailRatioEtaPGammaSumGamma->Scale(100);
	TH1D* histoCocktailRatioRhoGammaSumGamma = (TH1D*)histoCocktailRhoGamma->Clone("histoCocktailRatioRhoGammaSumGamma");
	histoCocktailRatioRhoGammaSumGamma->Divide(histoCocktailRatioRhoGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// 	histoCocktailRatioRhoGammaSumGamma->Scale(100);
	TH1D* histoCocktailRatioSigmaGammaSumGamma = (TH1D*)histoCocktailSigmaGamma->Clone("histoCocktailRatioSigmaGammaSumGamma");
	histoCocktailRatioSigmaGammaSumGamma->Divide(histoCocktailRatioSigmaGammaSumGamma,histoCocktailSumGamma, 1, 1,"B" );
// 	histoCocktailRatioSigmaGammaSumGamma->Scale(100);
	
	TCanvas *canvasCocktailRatio = new TCanvas("canvasCocktailRatio","",10,10,1200,1200);  // gives the page size		
	DrawGammaCanvasSettings( canvasCocktailRatio, 0.12, 0.01, 0.01, 0.08);
	canvasCocktailRatio->SetLogy();
	canvasCocktailRatio->SetLogx();
	
	Int_t textSizeLabelsPixelRatioCock = 48;
	Double_t textsizeLabelsRatioCock = 0;
	if (canvasCocktailRatio->XtoPixel(canvasCocktailRatio->GetX2()) < canvasCocktailRatio->YtoPixel(canvasCocktailRatio->GetY1())){
		textsizeLabelsRatioCock = (Double_t)textSizeLabelsPixelRatioCock/canvasCocktailRatio->XtoPixel(canvasCocktailRatio->GetX2()) ;
	} else {
		textsizeLabelsRatioCock = (Double_t)textSizeLabelsPixelRatioCock/canvasCocktailRatio->YtoPixel(canvasCocktailRatio->GetY1());
	}

	
	TH2D *dummyCocktailRatio ;
	dummyCocktailRatio = new TH2D("dummyCocktailRatio", "dummyCocktailRatio", 1000, 0., 22, 1000., 4e-5,14);
	SetStyleHistoTH2ForGraphs( dummyCocktailRatio, "#it{p}_{T} (GeV/#it{c})", "#gamma_{source}/#gamma_{decay}",
							   0.85*textsizeLabelsRatioCock, textsizeLabelsRatioCock, 0.85*textsizeLabelsRatioCock, textsizeLabelsRatioCock, 0.82, 1.3);
	dummyCocktailRatio->GetXaxis()->SetLabelOffset(-0.011);
// 	dummyCocktailRatio->GetXaxis()->SetMoreLogLabels();
	dummyCocktailRatio->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	
	dummyCocktailRatio->DrawCopy();
		SetStyleHisto(histoCocktailRatioPi0GammaSumGamma, widthCommonFit*1.5, 1, colorCocktailPi0 );
		SetStyleHisto(histoCocktailRatioEtaGammaSumGamma, widthCommonFit*1.5, 7, colorCocktailEta );
		SetStyleHisto(histoCocktailRatioEtaPGammaSumGamma, widthCommonFit*1.5, 2, colorCocktailEtaP );
		SetStyleHisto(histoCocktailRatioOmegaGammaSumGamma, widthCommonFit*1.5, 4, colorCocktailOmega );
		SetStyleHisto(histoCocktailRatioPhiGammaSumGamma, widthCommonFit*1.5, 5, colorCocktailPhi );
		SetStyleHisto(histoCocktailRatioRhoGammaSumGamma, widthCommonFit*1.5, 8, colorCocktailRho0 );
// 		SetStyleHisto(histoCocktailRatioSigmaGammaSumGamma, widthCommonFit*1.5, 3, colorCocktailSigma0 );
		
		histoCocktailRatioPi0GammaSumGamma->Draw("chistsame");
		histoCocktailRatioEtaGammaSumGamma->Draw("chistsame");
		histoCocktailRatioEtaPGammaSumGamma->Draw("chistsame");
		histoCocktailRatioOmegaGammaSumGamma->Draw("chistsame");
		histoCocktailRatioPhiGammaSumGamma->Draw("chistsame");
		histoCocktailRatioRhoGammaSumGamma->Draw("chistsame");
// 		histoCocktailRatioSigmaGammaSumGamma->Draw("chistsame");

	// 	TLatex* tpi = new TLatex(0.18,0.92,"#pi^{0} #rightarrow #gamma#gamma (e^{+}e^{-}#gamma)");
	// 	SetStyleTLatex( tpi, 0.85*textsizeLabelsRatioCock,4,colorCocktailPi0,42);
	// 	tpi->Draw();   
	// 	TLatex* teta = new TLatex(0.18,0.88,"#eta #rightarrow #gamma#gamma (#pi^{+}#pi^{-}#gamma,e^{+}e^{-}#gamma,#pi^{0}#gamma#gamma)");
	// 	SetStyleTLatex( teta, 00.85*textsizeLabelsRatioCock,4,colorCocktailEta,42);
	// 	teta->Draw();
	// 	TLatex* tomega = new TLatex(0.18,0.84,"#omega #rightarrow #pi^{0}#gamma (#eta#gamma)");
	// 	SetStyleTLatex( tomega, 0.85*textsizeLabelsRatioCock,4,colorCocktailOmega,42);
	// 	tomega->Draw();
	// 	TLatex* tetaprime = new TLatex(0.18,0.80,"#eta' #rightarrow #rho#gamma (#omega#gamma, #gamma#gamma)");
	// 	SetStyleTLatex( tetaprime, 0.85*textsizeLabelsRatioCock,4,colorCocktailEtaP,42);
	// 	tetaprime->Draw();
	// 	TLatex* tphi = new TLatex(0.65,0.92,"#phi #rightarrow #eta#gamma (#pi^{0}#gamma, #omega#gamma)");
	// 	SetStyleTLatex( tphi, 0.85*textsizeLabelsRatioCock,4,colorCocktailPhi,42);
	// 	tphi->Draw();
	// 	TLatex* trho = new TLatex(0.65,0.88,"#rho^{0} #rightarrow #pi^{+}#pi^{-}#gamma (#pi^{0}#gamma, #eta#gamma)");
	// 	SetStyleTLatex( trho, 0.85*textsizeLabelsRatioCock,4,colorCocktailRho0,42);
	// 	trho->Draw();
	// 	TLatex* tSigma = new TLatex(0.65,0.84,"#Sigma^{0} #rightarrow #Lambda#gamma (#Lambda#gamma#gamma)");
	// 	SetStyleTLatex( tSigma, 0.85*textsizeLabelsRatioCock,4,colorCocktailSigma0,42);
	// 	tSigma->Draw();

		TLegend* legendCocktail = new TLegend(0.15,0.96-1.1*0.85*textsizeLabelsRatioCock*4,0.18+0.75,0.96);
		legendCocktail->SetFillStyle(0);
		legendCocktail->SetFillColor(0);
		legendCocktail->SetLineColor(0);
		legendCocktail->SetTextSize(0.85*textsizeLabelsRatioCock);
		legendCocktail->SetMargin(0.15);
		legendCocktail->SetTextFont(42);
		legendCocktail->SetNColumns(2);
		legendCocktail->AddEntry(histoCocktailRatioPi0GammaSumGamma,"#pi^{0} #rightarrow #gamma#gamma (e^{+}e^{-}#gamma)","l");
		legendCocktail->AddEntry(histoCocktailRatioEtaPGammaSumGamma,"#eta' #rightarrow #rho#gamma (#omega#gamma, #gamma#gamma)","l");
		legendCocktail->AddEntry(histoCocktailRatioEtaGammaSumGamma,"#eta #rightarrow #gamma#gamma (#pi^{+}#pi^{-}#gamma,e^{+}e^{-}#gamma,#pi^{0}#gamma#gamma)","l");
		legendCocktail->AddEntry(histoCocktailRatioPhiGammaSumGamma,"#phi #rightarrow #eta#gamma (#pi^{0}#gamma, #omega#gamma)","l");
		legendCocktail->AddEntry(histoCocktailRatioOmegaGammaSumGamma,"#omega #rightarrow #pi^{0}#gamma (#eta#gamma)","l");
		legendCocktail->AddEntry(histoCocktailRatioRhoGammaSumGamma,"#rho^{0} #rightarrow #pi^{+}#pi^{-}#gamma (#pi^{0}#gamma, #eta#gamma)","l");
// 		legendCocktail->AddEntry(histoCocktailRatioSigmaGammaSumGamma,"#Sigma^{0} #rightarrow #Lambda#gamma (#Lambda#gamma#gamma)","l");
		
		legendCocktail->Draw();
		
		TLatex *labelCocktailRatioEnergy = new TLatex(0.16,0.13,collisionSystemCent0020.Data());
		SetStyleTLatex( labelCocktailRatioEnergy, 0.85*textsizeLabelsRatioCock,4);		
		labelCocktailRatioEnergy->Draw();

		TLatex *labelCocktailSim = new TLatex(0.625,0.13,"Monte Carlo simulation");
		SetStyleTLatex( labelCocktailSim, 0.85*textsizeLabelsRatioCock,4);		
		labelCocktailSim->Draw();
		
	
	canvasCocktailRatio->Print(Form("%s/CocktailRatioToSumPCM_0020.%s",outputDir.Data(),suffix.Data()));

	// *****************************************************************************************************
	// ***************************************** Plotting RAA 0-20% ****************************************
	// *****************************************************************************************************
	cout << "here" << endl;
	TGraphAsymmErrors* graphCombRAADirGammaSys0020Plot = NULL;
	if (graphCombRAADirGammaSys0020){
		graphCombRAADirGammaSys0020Plot = (TGraphAsymmErrors*)graphCombRAADirGammaSys0020->Clone("graphCombRAADirGammaSys0020Plot");
	}
	cout << "here" << endl;
	TGraphAsymmErrors* graphCombRAADirGammaStat0020Plot = NULL;
	if (graphCombRAADirGammaStat0020){
		graphCombRAADirGammaStat0020Plot = (TGraphAsymmErrors*)graphCombRAADirGammaStat0020->Clone("graphCombRAADirGammaStat0020Plot");
		ProduceGraphAsymmWithoutXErrors(graphCombRAADirGammaStat0020Plot);
	}
	
	cout << "here" << endl;
	TCanvas* canvasRAA_0020 = new TCanvas("canvasRAA_0020","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasRAA_0020,  0.1, 0.01, 0.015, 0.1);
	canvasRAA_0020->SetLogx(1);

	Int_t textSizeLabelsPixelRAA = 48;
	Double_t textsizeLabelsRAA = 0;
	if (canvasRAA_0020->XtoPixel(canvasRAA_0020->GetX2()) < canvasRAA_0020->YtoPixel(canvasRAA_0020->GetY1())){
		textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAA_0020->XtoPixel(canvasRAA_0020->GetX2()) ;
	} else {
		textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAA_0020->YtoPixel(canvasRAA_0020->GetY1());
	}

	TH2F * histo2DRAADummy;
	histo2DRAADummy = new TH2F("histo2DRAADummy","histo2DRAADummy",1000,0.,20.,1000,0.0,11.5);
	SetStyleHistoTH2ForGraphs(histo2DRAADummy, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA,0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.95,1., 510, 510); //#frac{#frac{1
	// 	histo2DRAAAll3->GetYaxis()->SetRangeUser(0.05,8.);
	histo2DRAADummy->GetYaxis()->SetLabelOffset(0.005);
	histo2DRAADummy->GetXaxis()->SetLabelOffset(-0.005);
	histo2DRAADummy->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	histo2DRAADummy->GetXaxis()->SetMoreLogLabels();
	histo2DRAADummy->DrawCopy("");

		TLatex *labelRAAEnergy0020 = new TLatex(0.48,0.92,collisionSystemCent0020.Data());
		SetStyleTLatex( labelRAAEnergy0020, 0.85*textsizeLabelsRAA,4);		
		labelRAAEnergy0020->Draw();
	
		TBox* boxErrorNorm0020_Single = CreateBoxConv(colorComb0020Box, 0.75, 1.-normErr0020 , 0.8, 1.+normErr0020);
		boxErrorNorm0020_Single->Draw();
	
		SetStyleGammaNLOTGraphWithBand( graphTheoryRelErrEPS090020, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
		graphTheoryRelErrEPS090020->Draw("p3lsame");
	
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1] , 1, 1 ,3,kGray+1,7);
		
		if (graphCombRAADirGammaSys0020Plot){
			DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaSys0020Plot, markerStyleComb0020,markerSizeComb0020, colorComb0020 , colorComb0020, widthLinesBoxes, kTRUE);
			graphCombRAADirGammaSys0020Plot->Draw("E2same");
		}	
		if (graphCombRAADirGammaStat0020Plot){
			DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaStat0020Plot, markerStyleComb0020,markerSizeComb0020, colorComb0020 , colorComb0020);
			graphCombRAADirGammaStat0020Plot->Draw("p,same,e1Z");
		}

		TLegend* legendRAA0020 = new TLegend(0.48,0.91-1.25*0.85*textsizeLabelsRAA*6,0.48+0.21,0.91);
		legendRAA0020->SetFillStyle(0);
		legendRAA0020->SetFillColor(0);
		legendRAA0020->SetLineColor(0);
		legendRAA0020->SetTextSize(0.85*textsizeLabelsRAA);
		legendRAA0020->SetMargin(0.3);
		legendRAA0020->SetTextFont(42);
		if (graphCombRAADirGammaSys0020Plot)legendRAA0020->AddEntry(graphCombRAADirGammaSys0020Plot,"ALICE #it{#gamma}_{dir}","fp");
		legendRAA0020->AddEntry((TObject*)0,"pp reference: ","");
		legendRAA0020->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
		legendRAA0020->AddEntry((TObject*)0,"FF: BFG2","");
		legendRAA0020->AddEntry(graphTheoryRelErrEPS090020,"Rel. error #it{JETPHOX}","fp");
		legendRAA0020->AddEntry((TObject*)0,"nPDF: EPS09, FF: BFG2","");
		legendRAA0020->Draw();

		histo2DRAADummy->Draw("axis,same");

	canvasRAA_0020->Update();
	canvasRAA_0020->Print(Form("%s/DirGammaRAA_0020.%s",outputDir.Data(),suffix.Data()));

	// *****************************************************************************************************
	// ***************************************** Plotting RAA 20-40% ****************************************
	// *****************************************************************************************************
	TGraphAsymmErrors* graphCombRAADirGammaSys2040Plot = NULL;
	if (graphCombRAADirGammaSys2040){
		graphCombRAADirGammaSys2040Plot = (TGraphAsymmErrors*)graphCombRAADirGammaSys2040->Clone("graphCombRAADirGammaSys2040Plot");
	}
	cout << "here" << endl;
	TGraphAsymmErrors* graphCombRAADirGammaStat2040Plot = NULL;
	if (graphCombRAADirGammaStat2040){
		graphCombRAADirGammaStat2040Plot = (TGraphAsymmErrors*)graphCombRAADirGammaStat2040->Clone("graphCombRAADirGammaStat2040Plot");
		ProduceGraphAsymmWithoutXErrors(graphCombRAADirGammaStat2040Plot);
	}

	TCanvas* canvasRAA_2040 = new TCanvas("canvasRAA_2040","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasRAA_2040,  0.1, 0.01, 0.015, 0.1);
	canvasRAA_2040->SetLogx(1);

	histo2DRAADummy->DrawCopy("");

		TBox* boxErrorNorm2040_Single = CreateBoxConv(colorComb2040Box, 0.75, 1.-normErr2040 , 0.8, 1.+normErr2040);
		boxErrorNorm2040_Single->Draw();

		TLatex *labelRAAEnergy2040 = new TLatex(0.48,0.92,collisionSystemCent2040.Data());
		SetStyleTLatex( labelRAAEnergy2040, 0.85*textsizeLabelsRAA,4);		
		labelRAAEnergy2040->Draw();
		
		SetStyleGammaNLOTGraphWithBand( graphTheoryRelErrEPS092040, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
		graphTheoryRelErrEPS092040->Draw("p3lsame");
	
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1] , 1, 1 ,3,kGray+1,7);

		if (graphCombRAADirGammaSys2040Plot){
			DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaSys2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);
			graphCombRAADirGammaSys2040Plot->Draw("E2same");
		}	
		if (graphCombRAADirGammaStat2040Plot){
			DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaStat2040Plot, markerStyleComb2040, markerSizeComb2040, colorComb2040 , colorComb2040);	
			graphCombRAADirGammaStat2040Plot->Draw("p,E1Z,same");
		}
		if (graphCombRAADirGammaSum2040Ar){
// 			graphCombRAADirGammaSum2040Ar->RemovePoint(0);
			DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaSum2040Ar , 1, 3, colorComb2040, colorComb2040, 1.8, kTRUE);
			for (Int_t i = 0; i < graphCombRAADirGammaSum2040Ar->GetN(); i++){
				graphCombRAADirGammaSum2040Ar->SetPointEYhigh(i,0);
			}	
			graphCombRAADirGammaSum2040Ar->Draw(">,same");
			PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombRAADirGammaSum2040Ar);
			graphCombRAADirGammaSum2040Ar->Print();
		}	
		histo2DRAADummy->Draw("axis,same");

		TLegend* legendRAA2040 = new TLegend(0.48,0.91-1.25*0.85*textsizeLabelsRAA*6,0.48+0.21,0.91);
		legendRAA2040->SetFillStyle(0);
		legendRAA2040->SetFillColor(0);
		legendRAA2040->SetLineColor(0);
		legendRAA2040->SetTextSize(0.85*textsizeLabelsRAA);
		legendRAA2040->SetMargin(0.3);
		legendRAA2040->SetTextFont(42);
		if (graphCombRAADirGammaSys2040Plot)legendRAA2040->AddEntry(graphCombRAADirGammaSys2040Plot,"ALICE #it{#gamma}_{dir}","fp");
		legendRAA2040->AddEntry((TObject*)0,"pp reference: ","");
		legendRAA2040->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
		legendRAA2040->AddEntry((TObject*)0,"FF: BFG2","");
		legendRAA2040->AddEntry(graphTheoryRelErrEPS092040,"Rel. error #it{JETPHOX}","fp");
		legendRAA2040->AddEntry((TObject*)0,"nPDF: EPS09, FF: BFG2","");
		legendRAA2040->Draw();
		
		
	canvasRAA_2040->Update();
	canvasRAA_2040->Print(Form("%s/DirGammaRAA_2040.%s",outputDir.Data(),suffix.Data()));
	
	// *****************************************************************************************************
	// ***************************************** Plotting RAA 40-80% ****************************************
	// *****************************************************************************************************
	TGraphAsymmErrors* graphCombRAADirGammaSys4080Plot = NULL;
	if (graphCombRAADirGammaSys4080){
		graphCombRAADirGammaSys4080Plot = (TGraphAsymmErrors*)graphCombRAADirGammaSys4080->Clone("graphCombRAADirGammaSys4080Plot");
	}
	cout << "here" << endl;
	TGraphAsymmErrors* graphCombRAADirGammaStat4080Plot = NULL;
	if (graphCombRAADirGammaStat4080){
		graphCombRAADirGammaStat4080Plot = (TGraphAsymmErrors*)graphCombRAADirGammaStat4080->Clone("graphCombRAADirGammaStat4080Plot");
		ProduceGraphAsymmWithoutXErrors(graphCombRAADirGammaStat4080Plot);
	}

	TCanvas* canvasRAA_4080 = new TCanvas("canvasRAA_4080","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasRAA_4080,  0.1, 0.01, 0.015, 0.1);
	canvasRAA_4080->SetLogx(1);

	histo2DRAADummy->DrawCopy("");

		TBox* boxErrorNorm4080_Single = CreateBoxConv(colorComb4080Box, 0.75, 1.-normErr4080 , 0.8, 1.+normErr4080);
		boxErrorNorm4080_Single->Draw();

		TLatex *labelRAAEnergy4080 = new TLatex(0.48,0.92,collisionSystemCent4080.Data());
		SetStyleTLatex( labelRAAEnergy4080, 0.85*textsizeLabelsRAA,4);		
		labelRAAEnergy4080->Draw();
		
		SetStyleGammaNLOTGraphWithBand( graphTheoryRelErrEPS094080, 1.0, 1, colorEPS09calc, fillStyleEPS09, colorEPS09calc, 0);
		graphTheoryRelErrEPS094080->Draw("p3lsame");
	
		DrawGammaLines(doubleRatioX[0], doubleRatioX[1] , 1, 1 ,3,kGray+1,7);

		if (graphCombRAADirGammaSys4080Plot){
			DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaSys4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080, widthLinesBoxes, kTRUE);
			graphCombRAADirGammaSys4080Plot->Draw("E2same");
		}	
		if (graphCombRAADirGammaStat4080Plot){
			DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaStat4080Plot, markerStyleComb4080, markerSizeComb4080, colorComb4080 , colorComb4080);	
			graphCombRAADirGammaStat4080Plot->Draw("p,E1Z,same");
		}
		if (graphCombRAADirGammaSum4080Ar){
// 			graphCombRAADirGammaSum4080Ar->RemovePoint(0);
			DrawGammaSetMarkerTGraphAsym(graphCombRAADirGammaSum4080Ar , 1, 3, colorComb4080, colorComb4080, 1.8, kTRUE);
			for (Int_t i = 0; i < graphCombRAADirGammaSum4080Ar->GetN(); i++){
				graphCombRAADirGammaSum4080Ar->SetPointEYhigh(i,0);
			}	
			graphCombRAADirGammaSum4080Ar->Draw(">,same");
			PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphCombRAADirGammaSum4080Ar);
			graphCombRAADirGammaSum4080Ar->Print();
		}	
		histo2DRAADummy->Draw("axis,same");

		TLegend* legendRAA4080 = new TLegend(0.48,0.91-1.25*0.85*textsizeLabelsRAA*6,0.48+0.21,0.91);
		legendRAA4080->SetFillStyle(0);
		legendRAA4080->SetFillColor(0);
		legendRAA4080->SetLineColor(0);
		legendRAA4080->SetTextSize(0.85*textsizeLabelsRAA);
		legendRAA4080->SetMargin(0.3);
		legendRAA4080->SetTextFont(42);
		if (graphCombRAADirGammaSys4080Plot)legendRAA4080->AddEntry(graphCombRAADirGammaSys4080Plot,"ALICE #it{#gamma}_{dir}","fp");
		legendRAA4080->AddEntry((TObject*)0,"pp reference: ","");
		legendRAA4080->AddEntry((TObject*)0,"(n)PDF: CTEQ6.1M/EPS09,","");
		legendRAA4080->AddEntry((TObject*)0,"FF: BFG2","");
		legendRAA4080->AddEntry(graphTheoryRelErrEPS094080,"Rel. error #it{JETPHOX}","fp");
		legendRAA4080->AddEntry((TObject*)0,"nPDF: EPS09, FF: BFG2","");
		legendRAA4080->Draw();
		
		
	canvasRAA_4080->Update();
	canvasRAA_4080->Print(Form("%s/DirGammaRAA_4080.%s",outputDir.Data(),suffix.Data()));

	//*******************************************************************************************************************
	//******************************************** Significance calculations ********************************************
	//*******************************************************************************************************************
	// originally written by Klaus Reygers 
	
	if (enablepValueCalc){
		
		//*********************************************************************************
		// Calculating Significance of combined double ratio for 0-20%
		//*********************************************************************************
		Int_t minBinSig0020 = 0;
		Int_t maxBinSig0020 = 5;
		const Int_t nPointSig0020 = maxBinSig0020- minBinSig0020 +1;
		
		// filling of graphs
		TGraph* graphAbsTypeAPlusStatErr0020 					= new TGraph(nPointSig0020);
		TGraph* graphRelTypeAPlusStatErr0020 					= new TGraph(nPointSig0020);
		TGraph* graphRelTypeBErr0020 							= new TGraph(nPointSig0020);
		TGraphErrors* graphDoubleRatioTypeAPlusStatErr0020 		= new TGraphErrors(nPointSig0020);
		for (Int_t i = 0; i < nPointSig0020; i++){
			Double_t pt 			= graphCombDRPi0FitSysBErr0020->GetX()[minBinSig0020+i];
			Double_t R				= graphCombDRPi0FitSysBErr0020->GetY()[minBinSig0020+i];
			Double_t statErr		= graphCombDRPi0FitStatErr0020->GetEYlow()[minBinSig0020+i];
			Double_t sysAErr		= graphCombDRPi0FitSysAErr0020->GetEYlow()[minBinSig0020+i];
			Double_t relErrSysA 	= sysAErr/ R;
			Double_t statSysAErr	= TMath::Sqrt(statErr*statErr+sysAErr*sysAErr);
			Double_t relErrStatSysA = statSysAErr/ R;
			Double_t relErrSysB		= graphCombDRPi0FitSysBErr0020->GetEYlow()[minBinSig0020+i]/ graphCombDRPi0FitSysBErr0020->GetY()[minBinSig0020+i];
			
			graphAbsTypeAPlusStatErr0020->SetPoint(i, pt, statSysAErr);
			graphRelTypeAPlusStatErr0020->SetPoint(i, pt, relErrStatSysA);
			graphRelTypeBErr0020->SetPoint(i, pt, relErrSysB);
			graphDoubleRatioTypeAPlusStatErr0020->SetPoint(i, pt, R);
			graphDoubleRatioTypeAPlusStatErr0020->SetPointError(i,0,statSysAErr);
		}	
		Double_t relErrSysC0020 = graphCombDRPi0FitSysCErr0020->GetEYlow()[minBinSig0020]/ graphCombDRPi0FitSysCErr0020->GetY()[minBinSig0020];
		
		// Generating pseudo data
		TH1F* histoChi2NullHypo0020 = new TH1F("histoChi2NullHypo0020", "histoChi2NullHypo0020", 3000, 0., 1500);
		FillChi2HistForNullHypo(100000000, histoChi2NullHypo0020, graphRelTypeAPlusStatErr0020, graphRelTypeBErr0020, relErrSysC0020, kFALSE);
		
		// calculating significance
		Double_t chi2Data0020 	= Chi2ForNullHypo(graphDoubleRatioTypeAPlusStatErr0020);
		Int_t binChi2Data0020 	= histoChi2NullHypo0020->FindBin(chi2Data0020);
		Int_t binFirst0020 		= 1;
		Int_t binLast0020		= histoChi2NullHypo0020->GetNbinsX();
		Double_t intTot0020		= histoChi2NullHypo0020->Integral(binFirst0020, binLast0020);
		Double_t intData0020	= histoChi2NullHypo0020->Integral(binChi2Data0020, binLast0020);
		Double_t pValue0020		= intData0020/ intTot0020;
		Double_t nSigma0020		= PValueToNSigma(pValue0020);
		
		cout << "0-20%: chi2data =" << chi2Data0020 << "\t, pVal = " << pValue0020 << "\t, nSigma = " << nSigma0020 << endl;
		
		// preparing histo for drawing
		TH1F* histoChi2NullData0020 = (TH1F*)histoChi2NullHypo0020->Clone("histoChi2NullData0020");
		for (Int_t i = 1; i < binChi2Data0020; i++ ){
			histoChi2NullData0020->SetBinContent(i,-1);
		}
		for (Int_t i = 1; i < binLast0020; i++ ){
			histoChi2NullData0020->SetBinError(i,0);
		}
		
		//*********************************************************************************
		// Calculating Significance of combined double ratio for 20-40%
		//*********************************************************************************
		Int_t minBinSig2040 = 0;
		Int_t maxBinSig2040 = 5;
		const Int_t nPointSig2040 = maxBinSig2040- minBinSig2040 +1;

		// filling of graphs		
		TGraph* graphAbsTypeAPlusStatErr2040 					= new TGraph(nPointSig2040);
		TGraph* graphRelTypeAPlusStatErr2040 					= new TGraph(nPointSig2040);
		TGraph* graphRelTypeBErr2040 							= new TGraph(nPointSig2040);
		TGraphErrors* graphDoubleRatioTypeAPlusStatErr2040 		= new TGraphErrors(nPointSig2040);

		for (Int_t i = 0; i < nPointSig2040; i++){
			Double_t pt 			= graphCombDRPi0FitSysBErr2040->GetX()[minBinSig2040+i];
			Double_t R				= graphCombDRPi0FitSysBErr2040->GetY()[minBinSig2040+i];
			Double_t statErr		= graphCombDRPi0FitStatErr2040->GetEYlow()[minBinSig2040+i];
			Double_t sysAErr		= graphCombDRPi0FitSysAErr2040->GetEYlow()[minBinSig2040+i];
			Double_t relErrSysA 	= sysAErr/ R;
			Double_t statSysAErr	= TMath::Sqrt(statErr*statErr+sysAErr*sysAErr);
			Double_t relErrStatSysA = statSysAErr/ R;
			Double_t relErrSysB		= graphCombDRPi0FitSysBErr2040->GetEYlow()[minBinSig2040+i]/ graphCombDRPi0FitSysBErr2040->GetY()[minBinSig2040+i];
			
			graphAbsTypeAPlusStatErr2040->SetPoint(i, pt, statSysAErr);
			graphRelTypeAPlusStatErr2040->SetPoint(i, pt, relErrStatSysA);
			graphRelTypeBErr2040->SetPoint(i, pt, relErrSysB);
			graphDoubleRatioTypeAPlusStatErr2040->SetPoint(i, pt, R);
			graphDoubleRatioTypeAPlusStatErr2040->SetPointError(i,0,statSysAErr);
		}	
		Double_t relErrSysC2040 = graphCombDRPi0FitSysCErr2040->GetEYlow()[minBinSig2040]/ graphCombDRPi0FitSysCErr2040->GetY()[minBinSig2040];
		
		// Generating pseudo data
		TH1F* histoChi2NullHypo2040 = new TH1F("histoChi2NullHypo2040", "histoChi2NullHypo2040", 3000, 0., 1500);
		FillChi2HistForNullHypo(100000000, histoChi2NullHypo2040, graphRelTypeAPlusStatErr2040, graphRelTypeBErr2040, relErrSysC2040, kFALSE);
		
		// calculating significance
		Double_t chi2Data2040 	= Chi2ForNullHypo(graphDoubleRatioTypeAPlusStatErr2040);
		Int_t binChi2Data2040 	= histoChi2NullHypo2040->FindBin(chi2Data2040);
		Int_t binFirst2040 		= 1;
		Int_t binLast2040		= histoChi2NullHypo2040->GetNbinsX();
		Double_t intTot2040		= histoChi2NullHypo2040->Integral(binFirst2040, binLast2040);
		Double_t intData2040	= histoChi2NullHypo2040->Integral(binChi2Data2040, binLast2040);
		Double_t pValue2040		= intData2040/ intTot2040;
		Double_t nSigma2040		= PValueToNSigma(pValue2040);
		
		cout << "20-40%:  chi2data =" << chi2Data2040 << "\t, pVal = " << pValue2040 << "\t, nSigma = " << nSigma2040 << endl;
		
		// preparing histo for drawing
		TH1F* histoChi2NullData2040 = (TH1F*)histoChi2NullHypo2040->Clone("histoChi2NullData2040");
		for (Int_t i = 1; i < binChi2Data2040; i++ ){
			histoChi2NullData2040->SetBinContent(i,-1);
		}
		for (Int_t i = 1; i < binLast2040; i++ ){
			histoChi2NullData2040->SetBinError(i,0);
		}
		
		//*********************************************************************************
		// Calculating Significance of combined double ratio for 40-80%
		//*********************************************************************************

		Int_t minBinSig4080 = 0;
		Int_t maxBinSig4080 = 5;
		const Int_t nPointSig4080 = maxBinSig4080- minBinSig4080 +1;
		
		// filling of graphs		
		TGraph* graphAbsTypeAPlusStatErr4080 					= new TGraph(nPointSig4080);
		TGraph* graphRelTypeAPlusStatErr4080 					= new TGraph(nPointSig4080);
		TGraph* graphRelTypeBErr4080 							= new TGraph(nPointSig4080);
		TGraphErrors* graphDoubleRatioTypeAPlusStatErr4080 		= new TGraphErrors(nPointSig4080);

		for (Int_t i = 0; i < nPointSig4080; i++){
			Double_t pt 			= graphCombDRPi0FitSysBErr4080->GetX()[minBinSig4080+i];
			Double_t R				= graphCombDRPi0FitSysBErr4080->GetY()[minBinSig4080+i];
			Double_t statErr		= graphCombDRPi0FitStatErr4080->GetEYlow()[minBinSig4080+i];
			Double_t sysAErr		= graphCombDRPi0FitSysAErr4080->GetEYlow()[minBinSig4080+i];
			Double_t relErrSysA 	= sysAErr/ R;
			Double_t statSysAErr	= TMath::Sqrt(statErr*statErr+sysAErr*sysAErr);
			Double_t relErrStatSysA = statSysAErr/ R;
			Double_t relErrSysB		= graphCombDRPi0FitSysBErr4080->GetEYlow()[minBinSig4080+i]/ graphCombDRPi0FitSysBErr4080->GetY()[minBinSig4080+i];
			
			graphAbsTypeAPlusStatErr4080->SetPoint(i, pt, statSysAErr);
			graphRelTypeAPlusStatErr4080->SetPoint(i, pt, relErrStatSysA);
			graphRelTypeBErr4080->SetPoint(i, pt, relErrSysB);
			graphDoubleRatioTypeAPlusStatErr4080->SetPoint(i, pt, R);
			graphDoubleRatioTypeAPlusStatErr4080->SetPointError(i,0,statSysAErr);
		}	
		Double_t relErrSysC4080 = graphCombDRPi0FitSysCErr4080->GetEYlow()[minBinSig4080]/ graphCombDRPi0FitSysCErr4080->GetY()[minBinSig4080];
		
		// Generating pseudo data
		TH1F* histoChi2NullHypo4080 = new TH1F("histoChi2NullHypo4080", "histoChi2NullHypo4080", 3000, 0., 1500);
		FillChi2HistForNullHypo(100000000, histoChi2NullHypo4080, graphRelTypeAPlusStatErr4080, graphRelTypeBErr4080, relErrSysC4080, kFALSE);
		
		// calculating significance
		Double_t chi2Data4080 	= Chi2ForNullHypo(graphDoubleRatioTypeAPlusStatErr4080);
		Int_t binChi2Data4080 	= histoChi2NullHypo4080->FindBin(chi2Data4080);
		Int_t binFirst4080 		= 1;
		Int_t binLast4080		= histoChi2NullHypo4080->GetNbinsX();
		Double_t intTot4080		= histoChi2NullHypo4080->Integral(binFirst4080, binLast4080);
		Double_t intData4080	= histoChi2NullHypo4080->Integral(binChi2Data4080, binLast4080);
		Double_t pValue4080		= intData4080/ intTot4080;
		Double_t nSigma4080		= PValueToNSigma(pValue4080);
		cout << "40-80%:  chi2data =" << chi2Data4080 << "\t, pVal = " << pValue4080 << "\t, nSigma = " << nSigma4080 << endl;
		
		// preparing histo for drawing
		TH1F* histoChi2NullData4080 = (TH1F*)histoChi2NullHypo4080->Clone("histoChi2NullData4080");
		for (Int_t i = 1; i < binChi2Data4080; i++ ){
			histoChi2NullData4080->SetBinContent(i,-1);
		}
		for (Int_t i = 1; i < binLast4080; i++ ){
			histoChi2NullData4080->SetBinError(i,0);
		}
	
		//*********************************************************************************
		// Plotting Significance of combined double ratio for 0-20%
		//*********************************************************************************

		TCanvas* canvasSignificance = new TCanvas("canvasSignificance","",200,10,1400,1100);  // gives the page size
		DrawGammaCanvasSettings( canvasSignificance,  0.09, 0.01, 0.015, 0.1);
		canvasSignificance->SetLogy(1);
		
		Int_t textSizeLabelsPixelSignificance = 48;
		Double_t textsizeLabelsSignificance = 0;
		if (canvasSignificance->XtoPixel(canvasSignificance->GetX2()) < canvasSignificance->YtoPixel(canvasSignificance->GetY1())){
			textsizeLabelsSignificance = (Double_t)textSizeLabelsPixelSignificance/canvasSignificance->XtoPixel(canvasSignificance->GetX2()) ;
		} else {
			textsizeLabelsSignificance = (Double_t)textSizeLabelsPixelSignificance/canvasSignificance->YtoPixel(canvasSignificance->GetY1());
		}

		SetStyleHistoTH1ForGraphs(histoChi2NullHypo0020, "Test statistics #it{t}","Counts", 
								0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
								0.85*textsizeLabelsSignificance, textsizeLabelsSignificance, 
								0.95, 1., 510, 510); 
		histoChi2NullHypo0020->GetYaxis()->SetLabelOffset(0.005);
		
		Int_t firstAbove0020 = histoChi2NullHypo0020->FindFirstBinAbove();
		Int_t lastAbove0020 = histoChi2NullHypo0020->FindLastBinAbove();
		
		histoChi2NullHypo0020->GetXaxis()->SetRange(1,lastAbove0020+1);
		DrawGammaSetMarker(histoChi2NullHypo0020, 1, 0.1, kGray+2 , kGray+2);
		histoChi2NullHypo0020->DrawCopy("");

		DrawGammaSetMarker(histoChi2NullData0020, 1, 0.1, colorComb0020Box , colorComb0020Box);
		histoChi2NullData0020->SetFillColor(colorComb0020Box);
		histoChi2NullData0020->SetFillStyle(3356);
		histoChi2NullData0020->Draw("same,lf");
		histoChi2NullHypo0020->DrawCopy("same");
		histoChi2NullHypo0020->DrawCopy("same,axis");
		
			TLatex *labelSignificanceEnergy0020 = new TLatex(0.58,0.92,collisionSystemCent0020.Data());
			SetStyleTLatex( labelSignificanceEnergy0020, 0.85*textsizeLabelsSignificance,4);		
			labelSignificanceEnergy0020->Draw();

			TLatex *labelSignificanceALICE = new TLatex(0.58,0.87,"ALICE pseudo data");
			SetStyleTLatex( labelSignificanceALICE, 0.85*textsizeLabelsSignificance,4);		
			labelSignificanceALICE->Draw();
			
			DrawGammaLines(chi2Data0020, chi2Data0020, 0, histoChi2NullHypo0020->GetMaximum()*0.1, 3, colorComb0020, 7);

			TLatex *labelTData0020 = new TLatex(chi2Data0020+10,histoChi2NullHypo0020->GetMaximum()*0.1,"#it{t}_{data}");
			SetStyleTLatex( labelTData0020, 0.85*textsizeLabelsSignificance,4,colorComb0020,42,kFALSE);
			labelTData0020->Draw();

			TLatex *labelSignificancePValue0020 = new TLatex(0.48,0.57,Form("#it{p}-value = %1.4f",pValue0020));
			SetStyleTLatex( labelSignificancePValue0020, 0.85*textsizeLabelsSignificance,4);		
			labelSignificancePValue0020->Draw();
			TLatex *labelSignificanceNSigma0020 = new TLatex(0.48,0.53,Form("%1.2f #sigma",nSigma0020));
			SetStyleTLatex( labelSignificanceNSigma0020, 0.85*textsizeLabelsSignificance,4);		
			labelSignificanceNSigma0020->Draw();
			
			
		canvasSignificance->Update();
		canvasSignificance->Print(Form("%s/DirGammaSignificanceTest_0020.%s",outputDir.Data(),suffix.Data()));

		//*********************************************************************************
		// Plotting Significance of combined double ratio for 20-40%
		//*********************************************************************************

		SetStyleHistoTH1ForGraphs(histoChi2NullHypo2040, "Test statistics #it{t}","Counts", 
								0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
								0.85*textsizeLabelsSignificance, textsizeLabelsSignificance, 
								0.95, 1., 510, 510); 
		histoChi2NullHypo2040->GetYaxis()->SetLabelOffset(0.005);
		
		Int_t firstAbove2040 = histoChi2NullHypo2040->FindFirstBinAbove();
		Int_t lastAbove2040 = histoChi2NullHypo2040->FindLastBinAbove();
		
		histoChi2NullHypo2040->GetXaxis()->SetRange(1,lastAbove2040+1);
		DrawGammaSetMarker(histoChi2NullHypo2040, 1, 0.1, kGray+2 , kGray+2);
		histoChi2NullHypo2040->DrawCopy("");

		DrawGammaSetMarker(histoChi2NullData2040, 1, 0.1, colorComb2040Box , colorComb2040Box);
		histoChi2NullData2040->SetFillColor(colorComb2040Box);
		histoChi2NullData2040->SetFillStyle(3356);
		histoChi2NullData2040->Draw("same,lf");
		histoChi2NullHypo2040->DrawCopy("same");
		histoChi2NullHypo2040->DrawCopy("same,axis");
		
			TLatex *labelSignificanceEnergy2040 = new TLatex(0.58,0.92,collisionSystemCent2040.Data());
			SetStyleTLatex( labelSignificanceEnergy2040, 0.85*textsizeLabelsSignificance,4);		
			labelSignificanceEnergy2040->Draw();

			labelSignificanceALICE->Draw();
			
			DrawGammaLines(chi2Data2040, chi2Data2040, 0, histoChi2NullHypo2040->GetMaximum()*0.5, 3, colorComb2040, 7);

			TLatex *labelTData2040 = new TLatex(chi2Data2040+10,histoChi2NullHypo2040->GetMaximum()*0.5,"#it{t}_{data}");
			SetStyleTLatex( labelTData2040, 0.85*textsizeLabelsSignificance,4,colorComb2040,42,kFALSE);
			labelTData2040->Draw();

			TLatex *labelSignificancePValue2040 = new TLatex(0.48,0.57,Form("#it{p}-value = %1.4f",pValue2040));
			SetStyleTLatex( labelSignificancePValue2040, 0.85*textsizeLabelsSignificance,4);		
			labelSignificancePValue2040->Draw();
			TLatex *labelSignificanceNSigma2040 = new TLatex(0.48,0.53,Form("%1.2f #sigma",nSigma2040));
			SetStyleTLatex( labelSignificanceNSigma2040, 0.85*textsizeLabelsSignificance,4);		
			labelSignificanceNSigma2040->Draw();
			
			
		canvasSignificance->Update();
		canvasSignificance->Print(Form("%s/DirGammaSignificanceTest_2040.%s",outputDir.Data(),suffix.Data()));

		//*********************************************************************************
		// Plotting Significance of combined double ratio for 40-80%
		//*********************************************************************************

		SetStyleHistoTH1ForGraphs(histoChi2NullHypo4080, "Test statistics #it{t}","Counts", 
								0.85*textsizeLabelsSignificance, textsizeLabelsSignificance,
								0.85*textsizeLabelsSignificance, textsizeLabelsSignificance, 
								0.95, 1., 510, 510); 
		histoChi2NullHypo4080->GetYaxis()->SetLabelOffset(0.005);
		
		Int_t firstAbove4080 = histoChi2NullHypo4080->FindFirstBinAbove();
		Int_t lastAbove4080 = histoChi2NullHypo4080->FindLastBinAbove();
		
		histoChi2NullHypo4080->GetXaxis()->SetRange(1,lastAbove4080+1);
		histoChi2NullHypo4080->GetYaxis()->SetRangeUser(1, 10*histoChi2NullHypo4080->GetMaximum());
		DrawGammaSetMarker(histoChi2NullHypo4080, 1, 0.1, kGray+2 , kGray+2);
		histoChi2NullHypo4080->DrawCopy("");

		DrawGammaSetMarker(histoChi2NullData4080, 1, 0.1, colorComb4080Box , colorComb4080Box);
		histoChi2NullData4080->SetFillColor(colorComb4080Box);
		histoChi2NullData4080->SetFillStyle(3356);
		histoChi2NullData4080->Draw("same,lf");
		histoChi2NullHypo4080->DrawCopy("same");
		histoChi2NullHypo4080->DrawCopy("same,axis");
		
			TLatex *labelSignificanceEnergy4080 = new TLatex(0.58,0.92,collisionSystemCent4080.Data());
			SetStyleTLatex( labelSignificanceEnergy4080, 0.85*textsizeLabelsSignificance,4);		
			labelSignificanceEnergy4080->Draw();

			labelSignificanceALICE->Draw();
			
			DrawGammaLines(chi2Data4080, chi2Data4080, 0, histoChi2NullHypo4080->GetMaximum()*0.2, 3, colorComb4080, 7);

			TLatex *labelTData4080 = new TLatex(chi2Data4080+10,histoChi2NullHypo4080->GetMaximum()*0.2,"#it{t}_{data}");
			SetStyleTLatex( labelTData4080, 0.85*textsizeLabelsSignificance,4,colorComb4080,42,kFALSE);
			labelTData4080->Draw();

			TLatex *labelSignificancePValue4080 = new TLatex(0.48,0.57,Form("#it{p}-value = %1.4f",pValue4080));
			SetStyleTLatex( labelSignificancePValue4080, 0.85*textsizeLabelsSignificance,4);		
			labelSignificancePValue4080->Draw();
			TLatex *labelSignificanceNSigma4080 = new TLatex(0.48,0.53,Form("%1.2f #sigma",nSigma4080));
			SetStyleTLatex( labelSignificanceNSigma4080, 0.85*textsizeLabelsSignificance,4);		
			labelSignificanceNSigma4080->Draw();
			
			
		canvasSignificance->Update();
		canvasSignificance->Print(Form("%s/DirGammaSignificanceTest_4080.%s",outputDir.Data(),suffix.Data()));
		
	}
	
	// ************************************************************************************************************************************************
	// ************************************************************ Write data to output file *********************************************************
	// ************************************************************************************************************************************************
	const char* fileNameOutputComp = "Gamma_CombResults_PbPb_2.76TeV.root";
	TFile* fileGammaSpectrum = new TFile(fileNameOutputComp,"UPDATE");		
		fileGammaSpectrum->mkdir("Gamma_PbPb_2.76TeV_0-20%");
		fileGammaSpectrum->cd("Gamma_PbPb_2.76TeV_0-20%");
			if (graphCombDRPi0FitStatErr0020) graphCombDRPi0FitStatErr0020->Write("DR_comb_StatErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysErr0020) graphCombDRPi0FitSysErr0020->Write("DR_comb_SysErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSumErr0020) graphCombDRPi0FitSumErr0020->Write("DR_comb_totErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysAErr0020) graphCombDRPi0FitSysAErr0020->Write("DR_comb_SysAErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysBErr0020) graphCombDRPi0FitSysBErr0020->Write("DR_comb_SysBErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysCErr0020) graphCombDRPi0FitSysCErr0020->Write("DR_comb_SysCErr",TObject::kOverwrite);

			if (graphCombIncGammaStatErr0020) graphCombIncGammaStatErr0020->Write("IncGammaSpec_comb_StatErr",TObject::kOverwrite);
			if (graphCombIncGammaSysErr0020) graphCombIncGammaSysErr0020->Write("IncGammaSpec_comb_SysErr",TObject::kOverwrite);
			if (graphCombIncGammaSumErr0020) graphCombIncGammaSumErr0020->Write("IncGammaSpec_comb_totErr",TObject::kOverwrite);
			if (graphCombIncGammaSysAErr0020) graphCombIncGammaSysAErr0020->Write("IncGammaSpec_comb_SysAErr",TObject::kOverwrite);
			if (graphCombIncGammaSysBErr0020) graphCombIncGammaSysBErr0020->Write("IncGammaSpec_comb_SysBErr",TObject::kOverwrite);
			if (graphCombIncGammaSysCErr0020) graphCombIncGammaSysCErr0020->Write("IncGammaSpec_comb_SysCErr",TObject::kOverwrite);
			
			if (graphCombDirGammaSpectrumStatErr0020) graphCombDirGammaSpectrumStatErr0020->Write("DirGammaSpec_comb_StatErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystErr0020) graphCombDirGammaSpectrumSystErr0020->Write("DirGammaSpec_comb_SysErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystAErr0020) graphCombDirGammaSpectrumSystAErr0020->Write("DirGammaSpec_comb_SysAErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystBErr0020) graphCombDirGammaSpectrumSystBErr0020->Write("DirGammaSpec_comb_SysBErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystCErr0020) graphCombDirGammaSpectrumSystCErr0020->Write("DirGammaSpec_comb_SysCErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSumErr0020) graphCombDirGammaSpectrumSumErr0020->Write("DirGammaSpec_comb_totErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSumErr0020Ar) graphCombDirGammaSpectrumSumErr0020Ar->Write("DirGammaSpec_comb_upperLimits",TObject::kOverwrite);

			if (graphCombThermalGammaSpectrumStatErr0020) graphCombThermalGammaSpectrumStatErr0020->Write("ThermalGammaSpec_comb_StatErr",TObject::kOverwrite);
			if (graphCombThermalGammaSpectrumSysErr0020) graphCombThermalGammaSpectrumSysErr0020->Write("ThermalGammaSpec_comb_SysErr",TObject::kOverwrite);
			if (fitPureThermalGamma0020Stat) fitPureThermalGamma0020Stat->Write("PureThermal_ExpFit_Stat",TObject::kOverwrite);
			
			if (histoPCMIncGammaStatErr0020) histoPCMIncGammaStatErr0020->Write("hIncGammaSpec_PCM_StatErr",TObject::kOverwrite);
			if (graphPCMIncGammaSysAErr0020) graphPCMIncGammaSysAErr0020->Write("IncGammaSpec_PCM_SysAErr",TObject::kOverwrite);
			if (graphPCMIncGammaSysBErr0020) graphPCMIncGammaSysBErr0020->Write("IncGammaSpec_PCM_SysBErr",TObject::kOverwrite);
			if (graphPCMIncGammaSysCErr0020) graphPCMIncGammaSysCErr0020->Write("IncGammaSpec_PCM_SysCErr",TObject::kOverwrite);

			if (histoPHOSIncGammaStatErr0020) histoPHOSIncGammaStatErr0020->Write("hIncGammaSpec_PHOS_StatErr",TObject::kOverwrite);
			if (graphPHOSIncGammaSysAErr0020) graphPHOSIncGammaSysAErr0020->Write("IncGammaSpec_PHOS_SysAErr",TObject::kOverwrite);
			if (graphPHOSIncGammaSysBErr0020) graphPHOSIncGammaSysBErr0020->Write("IncGammaSpec_PHOS_SysBErr",TObject::kOverwrite);
			if (graphPHOSIncGammaSysCErr0020) graphPHOSIncGammaSysCErr0020->Write("IncGammaSpec_PHOS_SysCErr",TObject::kOverwrite);

			if (histoPCMDRPi0FitStatErr0020) histoPCMDRPi0FitStatErr0020->Write("hDR_PCM_StatErr",TObject::kOverwrite);
			if (graphPCMDRPi0FitSysAErr0020) graphPCMDRPi0FitSysAErr0020->Write("DR_PCM_SysAErr",TObject::kOverwrite);
			if (graphPCMDRPi0FitSysBErr0020) graphPCMDRPi0FitSysBErr0020->Write("DR_PCM_SysBErr",TObject::kOverwrite);
			if (graphPCMDRPi0FitSysCErr0020) graphPCMDRPi0FitSysCErr0020->Write("DR_PCM_SysCErr",TObject::kOverwrite);

			if (histoPHOSDRPi0FitStatErr0020) histoPHOSDRPi0FitStatErr0020->Write("hDR_PHOS_StatErr",TObject::kOverwrite);
			if (graphPHOSDRPi0FitSysAErr0020) graphPHOSDRPi0FitSysAErr0020->Write("DR_PHOS_SysAErr",TObject::kOverwrite);
			if (graphPHOSDRPi0FitSysBErr0020) graphPHOSDRPi0FitSysBErr0020->Write("DR_PHOS_SysBErr",TObject::kOverwrite);
			if (graphPHOSDRPi0FitSysCErr0020) graphPHOSDRPi0FitSysCErr0020->Write("DR_PHOS_SysCErr",TObject::kOverwrite);

            if (graphCombRAADirGammaStat0020) graphCombRAADirGammaStat0020->Write("RAA_comb_StatErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSys0020) graphCombRAADirGammaSys0020->Write("RAA_comb_SysErr",TObject::kOverwrite);
			if (graphCombRAADirGammaSum0020Ar) graphCombRAADirGammaSum0020Ar->Write("RAA_comb_upperLimits",TObject::kOverwrite);
            
		fileGammaSpectrum->mkdir("Gamma_PbPb_2.76TeV_20-40%");
		fileGammaSpectrum->cd("Gamma_PbPb_2.76TeV_20-40%");
			if (graphCombDRPi0FitStatErr2040) graphCombDRPi0FitStatErr2040->Write("DR_comb_StatErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysErr2040) graphCombDRPi0FitSysErr2040->Write("DR_comb_SysErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSumErr2040) graphCombDRPi0FitSumErr2040->Write("DR_comb_totErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysAErr2040) graphCombDRPi0FitSysAErr2040->Write("DR_comb_SysAErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysBErr2040) graphCombDRPi0FitSysBErr2040->Write("DR_comb_SysBErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysCErr2040) graphCombDRPi0FitSysCErr2040->Write("DR_comb_SysCErr",TObject::kOverwrite);

			if (graphCombIncGammaStatErr2040) graphCombIncGammaStatErr2040->Write("IncGammaSpec_comb_StatErr",TObject::kOverwrite);
			if (graphCombIncGammaSysErr2040) graphCombIncGammaSysErr2040->Write("IncGammaSpec_comb_SysErr",TObject::kOverwrite);
			if (graphCombIncGammaSumErr2040) graphCombIncGammaSumErr2040->Write("IncGammaSpec_comb_totErr",TObject::kOverwrite);
			if (graphCombIncGammaSysAErr2040) graphCombIncGammaSysAErr2040->Write("IncGammaSpec_comb_SysAErr",TObject::kOverwrite);
			if (graphCombIncGammaSysBErr2040) graphCombIncGammaSysBErr2040->Write("IncGammaSpec_comb_SysBErr",TObject::kOverwrite);
			if (graphCombIncGammaSysCErr2040) graphCombIncGammaSysCErr2040->Write("IncGammaSpec_comb_SysCErr",TObject::kOverwrite);
			
			if (graphCombDirGammaSpectrumStatErr2040) graphCombDirGammaSpectrumStatErr2040->Write("DirGammaSpec_comb_StatErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystErr2040) graphCombDirGammaSpectrumSystErr2040->Write("DirGammaSpec_comb_SysErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystAErr2040) graphCombDirGammaSpectrumSystAErr2040->Write("DirGammaSpec_comb_SysAErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystBErr2040) graphCombDirGammaSpectrumSystBErr2040->Write("DirGammaSpec_comb_SysBErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystCErr2040) graphCombDirGammaSpectrumSystCErr2040->Write("DirGammaSpec_comb_SysCErr",TObject::kOverwrite);

			if (graphCombDirGammaSpectrumSumErr2040) graphCombDirGammaSpectrumSumErr2040->Write("DirGammaSpec_comb_totErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSumErr2040Ar) graphCombDirGammaSpectrumSumErr2040Ar->Write("DirGammaSpec_comb_upperLimits",TObject::kOverwrite);

			if (histoPCMIncGammaStatErr2040) histoPCMIncGammaStatErr2040->Write("hIncGammaSpec_PCM_StatErr",TObject::kOverwrite);
			if (graphPCMIncGammaSysAErr2040) graphPCMIncGammaSysAErr2040->Write("IncGammaSpec_PCM_SysAErr",TObject::kOverwrite);
			if (graphPCMIncGammaSysBErr2040) graphPCMIncGammaSysBErr2040->Write("IncGammaSpec_PCM_SysBErr",TObject::kOverwrite);
			if (graphPCMIncGammaSysCErr2040) graphPCMIncGammaSysCErr2040->Write("IncGammaSpec_PCM_SysCErr",TObject::kOverwrite);

			if (histoPHOSIncGammaStatErr2040) histoPHOSIncGammaStatErr2040->Write("hIncGammaSpec_PHOS_StatErr",TObject::kOverwrite);
			if (graphPHOSIncGammaSysAErr2040) graphPHOSIncGammaSysAErr2040->Write("IncGammaSpec_PHOS_SysAErr",TObject::kOverwrite);
			if (graphPHOSIncGammaSysBErr2040) graphPHOSIncGammaSysBErr2040->Write("IncGammaSpec_PHOS_SysBErr",TObject::kOverwrite);
			if (graphPHOSIncGammaSysCErr2040) graphPHOSIncGammaSysCErr2040->Write("IncGammaSpec_PHOS_SysCErr",TObject::kOverwrite);

			if (histoPCMDRPi0FitStatErr2040) histoPCMDRPi0FitStatErr2040->Write("hDR_PCM_StatErr",TObject::kOverwrite);
			if (graphPCMDRPi0FitSysAErr2040) graphPCMDRPi0FitSysAErr2040->Write("DR_PCM_SysAErr",TObject::kOverwrite);
			if (graphPCMDRPi0FitSysBErr2040) graphPCMDRPi0FitSysBErr2040->Write("DR_PCM_SysBErr",TObject::kOverwrite);
			if (graphPCMDRPi0FitSysCErr2040) graphPCMDRPi0FitSysCErr2040->Write("DR_PCM_SysCErr",TObject::kOverwrite);

			if (histoPHOSDRPi0FitStatErr2040) histoPHOSDRPi0FitStatErr2040->Write("hDR_PHOS_StatErr",TObject::kOverwrite);
			if (graphPHOSDRPi0FitSysAErr2040) graphPHOSDRPi0FitSysAErr2040->Write("DR_PHOS_SysAErr",TObject::kOverwrite);
			if (graphPHOSDRPi0FitSysBErr2040) graphPHOSDRPi0FitSysBErr2040->Write("DR_PHOS_SysBErr",TObject::kOverwrite);
			if (graphPHOSDRPi0FitSysCErr2040) graphPHOSDRPi0FitSysCErr2040->Write("DR_PHOS_SysCErr",TObject::kOverwrite);

            if (graphCombRAADirGammaStat2040) graphCombRAADirGammaStat2040->Write("RAA_comb_StatErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSys2040) graphCombRAADirGammaSys2040->Write("RAA_comb_SysErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSum2040Ar) graphCombRAADirGammaSum2040Ar->Write("RAA_comb_upperLimits",TObject::kOverwrite);
			
		fileGammaSpectrum->mkdir("Gamma_PbPb_2.76TeV_40-80%");
		fileGammaSpectrum->cd("Gamma_PbPb_2.76TeV_40-80%");
			if (graphCombDRPi0FitStatErr4080) graphCombDRPi0FitStatErr4080->Write("DR_comb_StatErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysErr4080) graphCombDRPi0FitSysErr4080->Write("DR_comb_SysErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSumErr4080) graphCombDRPi0FitSumErr4080->Write("DR_comb_totErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysAErr4080) graphCombDRPi0FitSysAErr4080->Write("DR_comb_SysAErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysBErr4080) graphCombDRPi0FitSysBErr4080->Write("DR_comb_SysBErr",TObject::kOverwrite);
			if (graphCombDRPi0FitSysCErr4080) graphCombDRPi0FitSysCErr4080->Write("DR_comb_SysCErr",TObject::kOverwrite);

			if (graphCombIncGammaStatErr4080) graphCombIncGammaStatErr4080->Write("IncGammaSpec_comb_StatErr",TObject::kOverwrite);
			if (graphCombIncGammaSysErr4080) graphCombIncGammaSysErr4080->Write("IncGammaSpec_comb_SysErr",TObject::kOverwrite);
			if (graphCombIncGammaSumErr4080) graphCombIncGammaSumErr4080->Write("IncGammaSpec_comb_totErr",TObject::kOverwrite);
			if (graphCombIncGammaSysAErr4080) graphCombIncGammaSysAErr4080->Write("IncGammaSpec_comb_SysAErr",TObject::kOverwrite);
			if (graphCombIncGammaSysBErr4080) graphCombIncGammaSysBErr4080->Write("IncGammaSpec_comb_SysBErr",TObject::kOverwrite);
			if (graphCombIncGammaSysCErr4080) graphCombIncGammaSysCErr4080->Write("IncGammaSpec_comb_SysCErr",TObject::kOverwrite);
			
			if (graphCombDirGammaSpectrumStatErr4080) graphCombDirGammaSpectrumStatErr4080->Write("DirGammaSpec_comb_StatErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystErr4080) graphCombDirGammaSpectrumSystErr4080->Write("DirGammaSpec_comb_SysErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystAErr4080) graphCombDirGammaSpectrumSystAErr4080->Write("DirGammaSpec_comb_SysAErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystBErr4080) graphCombDirGammaSpectrumSystBErr4080->Write("DirGammaSpec_comb_SysBErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSystCErr4080) graphCombDirGammaSpectrumSystCErr4080->Write("DirGammaSpec_comb_SysCErr",TObject::kOverwrite);

			if (graphCombDirGammaSpectrumSumErr4080) graphCombDirGammaSpectrumSumErr4080->Write("DirGammaSpec_comb_totErr",TObject::kOverwrite);
			if (graphCombDirGammaSpectrumSumErr4080Ar) graphCombDirGammaSpectrumSumErr4080Ar->Write("DirGammaSpec_comb_upperLimits",TObject::kOverwrite);
			
			if (histoPCMIncGammaStatErr4080) histoPCMIncGammaStatErr4080->Write("hIncGammaSpec_PCM_StatErr",TObject::kOverwrite);
			if (graphPCMIncGammaSysAErr4080) graphPCMIncGammaSysAErr4080->Write("IncGammaSpec_PCM_SysAErr",TObject::kOverwrite);
			if (graphPCMIncGammaSysBErr4080) graphPCMIncGammaSysBErr4080->Write("IncGammaSpec_PCM_SysBErr",TObject::kOverwrite);
			if (graphPCMIncGammaSysCErr4080) graphPCMIncGammaSysCErr4080->Write("IncGammaSpec_PCM_SysCErr",TObject::kOverwrite);

			if (histoPHOSIncGammaStatErr4080) histoPHOSIncGammaStatErr4080->Write("hIncGammaSpec_PHOS_StatErr",TObject::kOverwrite);
			if (graphPHOSIncGammaSysAErr4080) graphPHOSIncGammaSysAErr4080->Write("IncGammaSpec_PHOS_SysAErr",TObject::kOverwrite);
			if (graphPHOSIncGammaSysBErr4080) graphPHOSIncGammaSysBErr4080->Write("IncGammaSpec_PHOS_SysBErr",TObject::kOverwrite);
			if (graphPHOSIncGammaSysCErr4080) graphPHOSIncGammaSysCErr4080->Write("IncGammaSpec_PHOS_SysCErr",TObject::kOverwrite);

			if (histoPCMDRPi0FitStatErr4080) histoPCMDRPi0FitStatErr4080->Write("hDR_PCM_StatErr",TObject::kOverwrite);
			if (graphPCMDRPi0FitSysAErr4080) graphPCMDRPi0FitSysAErr4080->Write("DR_PCM_SysAErr",TObject::kOverwrite);
			if (graphPCMDRPi0FitSysBErr4080) graphPCMDRPi0FitSysBErr4080->Write("DR_PCM_SysBErr",TObject::kOverwrite);
			if (graphPCMDRPi0FitSysCErr4080) graphPCMDRPi0FitSysCErr4080->Write("DR_PCM_SysCErr",TObject::kOverwrite);

			if (histoPHOSDRPi0FitStatErr4080) histoPHOSDRPi0FitStatErr4080->Write("hDR_PHOS_StatErr",TObject::kOverwrite);
			if (graphPHOSDRPi0FitSysAErr4080) graphPHOSDRPi0FitSysAErr4080->Write("DR_PHOS_SysAErr",TObject::kOverwrite);
			if (graphPHOSDRPi0FitSysBErr4080) graphPHOSDRPi0FitSysBErr4080->Write("DR_PHOS_SysBErr",TObject::kOverwrite);
			if (graphPHOSDRPi0FitSysCErr4080) graphPHOSDRPi0FitSysCErr4080->Write("DR_PHOS_SysCErr",TObject::kOverwrite);

            if (graphCombRAADirGammaStat4080) graphCombRAADirGammaStat4080->Write("RAA_comb_StatErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSys4080) graphCombRAADirGammaSys4080->Write("RAA_comb_SysErr",TObject::kOverwrite);
            if (graphCombRAADirGammaSum4080Ar) graphCombRAADirGammaSum4080Ar->Write("RAA_comb_upperLimits",TObject::kOverwrite);
            
	fileGammaSpectrum->Write();
	fileGammaSpectrum->Close();

}