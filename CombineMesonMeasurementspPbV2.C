/****************************************************************************************************************************
****** 		provided by Gamma Conversion Group, PWG4, 													*****
******		Ana Marin, marin@physi.uni-heidelberg.de													*****
******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
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
#include <TPaveText.h>
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
#include "CombineMesonMeasurementspPbV2.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

void CombineMesonMeasurementspPbV2(const char *fileNameConversionpPb = "", TString suffix = "eps", TString bWCorrection="X"){	
	
	date = ReturnDateString();
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();

	//___________________________________ Declaration of files _____________________________________________
	collisionSystempPb 		= "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
	TString dateForOutput 	= ReturnDateStringForOutput();
	TString outputDir 		= Form("%s/%s/CombineMesonMeasurementspPb%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	nameFinalResDat 		= Form("%s/CombinedResultspPb%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
	fileNameCaloEMCALpPb 	= "ExternalInputpPb/EMCAL/EMCALResults_11Sept.root";
	fileNamePCMEMCALpPb 	= "ExternalInputpPb/Hybrid/data_PCM-EMCALResults_pPb.root";
	fileNamePCMPHOSpPb 	= "ExternalInputpPb/Hybrid/data_PCM-PHOSResults_pPb.root";
	fileNameEMCALEMCALpPb 	= "ExternalInputpPb/Hybrid/data_EMCAL-EMCALResults_pPb.root";
	fileNamePHOSPHOSpPb 	= "ExternalInputpPb/Hybrid/data_PHOS-PHOSResults_pPb.root";
	fileNameDalitzpPb	= "ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_2014-09-12.root";
	fileNameCaloPHOSpPb 	= "ExternalInputpPb/PHOS/data_PHOSResultsFullCorrection_pPb-23102014.root";   
	
	TString datePHOS = "2013-05-15";
	TString datePCM = "2013-06-18";

	TString nameHistoPCM = "CorrectedYieldPi0";
	TString nameHistoPCMEta = "CorrectedYieldEta";
	TString nameGraphPCM = "Pi0SystError";
	TString nameGraphPCMEta = "EtaSystError";
	TString nameHistoDalitz = "Pi0SystError";

	minPtForFits = 0.4;
	
	//declaration for printing logo 	
	prefix2 = "data";
	pictDrawingOptions[1] = kFALSE;
	
	mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	cout << mesonMassExpectPi0 << endl;
	mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();
	cout << mesonMassExpectEta << endl;
	
	//****************************************************************************************************
	//************************** Read data for PCM *******************************************************
	//****************************************************************************************************
	cout << "PCM" << endl;
	cout << "0-100%" << endl;
	filePCMpPb 								= new TFile(fileNameConversionpPb);
	directoryPCMPi0pPb 						= (TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
	histoPCMNumberOfEventspPb 				= (TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV0-100%");
	histoPCMYieldPi0pPb 					= (TH1D*)directoryPCMPi0pPb->Get(nameHistoPCM.Data());
	histoPCMRAWYieldPi0pPb 					= (TH1D*)directoryPCMPi0pPb->Get("RawYieldPi0");
	graphPCMYieldPi0SysErrpPb				= (TGraphAsymmErrors*)directoryPCMPi0pPb->Get(nameGraphPCM.Data());	
	graphPCMYieldPi0SysErrRAApPb 			= (TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystErrorA");	
	histoPCMMassPi0DatapPb 					= (TH1D*)directoryPCMPi0pPb->Get("MassPi0");
	histoPCMMassPi0MCpPb 					= (TH1D*)directoryPCMPi0pPb->Get("TrueMassPi0");
	histoPCMWidthPi0DatapPb 				= (TH1D*)directoryPCMPi0pPb->Get("FWHMPi0MeV");
	histoPCMWidthPi0MCpPb 					= (TH1D*)directoryPCMPi0pPb->Get("TrueFWHMPi0MeV");
	nPCMEventpPb 							= histoPCMNumberOfEventspPb->GetBinContent(1);
	histoRatioPCMMassPi0DiffDataMC 			= (TH1D*)histoPCMMassPi0DatapPb->Clone("histoRatioPCMMassPi0DiffDataMC");
	histoRatioPCMMassPi0DiffDataMC->Add(histoPCMMassPi0MCpPb,-1);
	histoRatioPCMMassPi0DiffDataMC->Scale(1/mesonMassExpectPi0*100);
	histoRatioPCMWidthPi0DiffDataMC 		= (TH1D*)histoPCMWidthPi0DatapPb->Clone("histoRatioPCMWidthPi0DiffDataMC");
	histoRatioPCMWidthPi0DiffDataMC->Add(histoPCMWidthPi0MCpPb,-1);
	histoRatioPCMWidthPi0DiffDataMC->Divide(histoRatioPCMWidthPi0DiffDataMC,histoPCMWidthPi0DatapPb,1.,1.,"");
	histoRatioPCMWidthPi0DiffDataMC->Scale(100.);
	
	histoPCMMassPi0DatapPb->Scale(1000.);
	histoPCMMassPi0MCpPb->Scale(1000.);
	
	// remove last points
	for (Int_t i = 0; i < 1; i++){
		histoPCMMassPi0DatapPb->SetBinContent(histoPCMMassPi0DatapPb->GetNbinsX()-i,0);
		histoPCMMassPi0MCpPb->SetBinContent(histoPCMMassPi0MCpPb->GetNbinsX()-i,0);
		histoPCMWidthPi0DatapPb->SetBinContent(histoPCMWidthPi0DatapPb->GetNbinsX()-i,10000.);
		histoPCMWidthPi0MCpPb->SetBinContent(histoPCMWidthPi0MCpPb->GetNbinsX()-i,10000.);
		histoRatioPCMWidthPi0DiffDataMC->SetBinContent(histoRatioPCMWidthPi0DiffDataMC->GetNbinsX()-i,10000.);
		histoRatioPCMMassPi0DiffDataMC->SetBinContent(histoRatioPCMMassPi0DiffDataMC->GetNbinsX()-i,10000.);
	}
   
	//****************************************************************************************************
	//************************** Read Eta data for PCM *******************************************************
	//****************************************************************************************************
	cout << "PCM Eta" << endl;
	cout << "0-100%" << endl;	
	directoryPCMEtapPb 						= (TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_0-100%"); 
	histoPCMYieldEtapPb 					= (TH1D*)directoryPCMEtapPb->Get(nameHistoPCMEta.Data());
	histoPCMRAWYieldEtapPb 					= (TH1D*)directoryPCMEtapPb->Get("RawYieldEta");
	graphPCMYieldEtaSysErrpPb				= (TGraphAsymmErrors*)directoryPCMEtapPb->Get(nameGraphPCMEta.Data());	
	graphPCMYieldEtaSysErrRAApPb 			= (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtaSystErrorA");	
	histoPCMEtaPi0RatiopPb 					= (TH1D*)directoryPCMEtapPb->Get("EtatoPi0Ratio");
	graphPCMEtaPi0RatioSysErrpPb 			= (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtatoPi0RatioSys"); 
	histoPCMMassEtaDatapPb 					= (TH1D*)directoryPCMEtapPb->Get("MassEta");
	histoPCMMassEtaMCpPb 					= (TH1D*)directoryPCMEtapPb->Get("TrueMassEta");
	histoPCMWidthEtaDatapPb 				= (TH1D*)directoryPCMEtapPb->Get("FWHMEtaMeV");
	histoPCMWidthEtaMCpPb 					= (TH1D*)directoryPCMEtapPb->Get("TrueFWHMEtaMeV");
	histoRatioPCMMassEtaDiffDataMC 	= (TH1D*)histoPCMMassEtaDatapPb->Clone("histoRatioPCMMassEtaDiffDataMC");
	histoRatioPCMMassEtaDiffDataMC->Add(histoPCMMassEtaMCpPb,-1);
	histoRatioPCMMassEtaDiffDataMC->Scale(1/(mesonMassExpectEta)*100);	
	histoPCMMassEtaDatapPb->Scale(1000.);
	histoPCMMassEtaMCpPb->Scale(1000.);
   
	//****************************************************************************************************
	//************************** Read data for EMCAL *****************************************************
	//****************************************************************************************************
	cout << "EMCAL" << endl;
	cout << "0-100%" << endl;
	fileEMCALpPb 							= new TFile(fileNameCaloEMCALpPb);
	directoryPi0EMCALpPb 					= (TDirectory*)fileEMCALpPb->Get("Pi05.02TeV_pPb"); 
	histoEMCALYieldPi0pPb 					= (TH1D*)directoryPi0EMCALpPb->Get("CorrectedYieldPi0");
	histoEMCALYieldPi0pPbSys				= (TH1D*)directoryPi0EMCALpPb->Get("Pi0SystError");
	histoEMCALRawYieldPi0pPb 				= (TH1D*)directoryPi0EMCALpPb->Get("RawYieldPerEventPi0");
	histoEMCALMassPi0DatapPb 				= (TH1D*)directoryPi0EMCALpPb->Get("MassPi0");
	histoEMCALMassPi0MCpPb 					= (TH1D*)directoryPi0EMCALpPb->Get("TrueMassPi0");
	histoEMCALWidthPi0DatapPb 				= (TH1D*)directoryPi0EMCALpPb->Get("FWHMPi0");
	histoEMCALWidthPi0MCpPb 				= (TH1D*)directoryPi0EMCALpPb->Get("TrueFWHMPi0MeV");
	histoRatioEMCALMassPi0DiffDataMC 		= (TH1D*)histoEMCALMassPi0DatapPb->Clone("histoRatioEMCALMassPi0DiffDataMC");
	histoRatioEMCALMassPi0DiffDataMC->Add(histoEMCALMassPi0MCpPb,-1);
	histoRatioEMCALMassPi0DiffDataMC->Scale(1/(mesonMassExpectPi0*1000)*100);
	histoRatioEMCALWidthPi0DiffDataMC 		= (TH1D*)histoEMCALWidthPi0DatapPb->Clone("histoRatioEMCALWidthPi0DiffDataMC");
	histoRatioEMCALWidthPi0DiffDataMC->Add(histoEMCALWidthPi0MCpPb,-1);
	histoRatioEMCALWidthPi0DiffDataMC->Divide(histoRatioEMCALWidthPi0DiffDataMC,histoEMCALWidthPi0DatapPb,1.,1.,"");
	histoRatioEMCALWidthPi0DiffDataMC->Scale(100.);
	histoEMCALWidthPi0DatapPb->Scale(1/2.355);
	histoEMCALWidthPi0MCpPb->Scale(1/2.355);
	
	//****************************************************************************************************
	//************************** Read data for PCM-EMCAL *******************************************************
	//****************************************************************************************************
	cout << "PCM-EMCAL" << endl;
	cout << "0-100%" << endl;
	filePCMEMCALpPb 						= new TFile(fileNamePCMEMCALpPb);
	directoryPCMEMCALPi0pPb 				= (TDirectory*)filePCMEMCALpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
	histoPCMEMCALYieldPi0pPb 				= (TH1D*)directoryPCMEMCALPi0pPb->Get("CorrectedYieldNormEffPi0");
	histoPCMEMCALRAWYieldPi0pPb 			= (TH1D*)directoryPCMEMCALPi0pPb->Get("RawYieldPi0PerEventsPi0");
	graphPCMEMCALYieldPi0SysErrpPb 			= (TGraphAsymmErrors*)directoryPCMEMCALPi0pPb->Get(nameGraphPCM.Data());	
	graphPCMEMCALYieldPi0SysErrRAApPb		= (TGraphAsymmErrors*)directoryPCMEMCALPi0pPb->Get("Pi0SystErrorA");	
	histoPCMEMCALMassPi0DatapPb 			= (TH1D*)directoryPCMEMCALPi0pPb->Get("MassPi0");
	histoPCMEMCALMassPi0MCpPb 				= (TH1D*)directoryPCMEMCALPi0pPb->Get("TrueMassPi0");
	histoPCMEMCALWidthPi0DatapPb 			= (TH1D*)directoryPCMEMCALPi0pPb->Get("FWHMPi0MeV");
	histoPCMEMCALWidthPi0MCpPb 				= (TH1D*)directoryPCMEMCALPi0pPb->Get("TrueFWHMPi0MeV");
	histoRatioPCMEMCALMassPi0DiffDataMC 	= (TH1D*)histoPCMEMCALMassPi0DatapPb->Clone("histoRatioPCMEMCALMassPi0DiffDataMC");
	histoRatioPCMEMCALMassPi0DiffDataMC->Add(histoPCMEMCALMassPi0MCpPb,-1);
	histoRatioPCMEMCALMassPi0DiffDataMC->Scale(1/(mesonMassExpectPi0)*100);	
	histoRatioPCMEMCALWidthPi0DiffDataMC 		= (TH1D*)histoPCMEMCALWidthPi0DatapPb->Clone("histoRatioPCMEMCALWidthPi0DiffDataMC");
	histoRatioPCMEMCALWidthPi0DiffDataMC->Add(histoPCMEMCALWidthPi0MCpPb,-1);
	histoRatioPCMEMCALWidthPi0DiffDataMC->Divide(histoRatioPCMEMCALWidthPi0DiffDataMC,histoPCMEMCALWidthPi0DatapPb,1.,1.,"");
	histoRatioPCMEMCALWidthPi0DiffDataMC->Scale(100.);
	histoPCMEMCALMassPi0DatapPb->Scale(1000.);
	histoPCMEMCALMassPi0MCpPb->Scale(1000.);

	// remove last points
	for (Int_t i = 0; i < 3; i++){
		histoPCMEMCALMassPi0DatapPb->SetBinContent(histoPCMEMCALMassPi0DatapPb->GetNbinsX()-i,0);
		histoPCMEMCALMassPi0MCpPb->SetBinContent(histoPCMEMCALMassPi0MCpPb->GetNbinsX()-i,0);
		histoPCMEMCALWidthPi0DatapPb->SetBinContent(histoPCMEMCALWidthPi0DatapPb->GetNbinsX()-i,10000.);
		histoPCMEMCALWidthPi0MCpPb->SetBinContent(histoPCMEMCALWidthPi0MCpPb->GetNbinsX()-i,10000.);
		histoRatioPCMEMCALMassPi0DiffDataMC->SetBinContent(histoRatioPCMEMCALMassPi0DiffDataMC->GetNbinsX()-i,-100);
		histoRatioPCMEMCALWidthPi0DiffDataMC->SetBinContent(histoRatioPCMEMCALWidthPi0DiffDataMC->GetNbinsX()-i,-100);
	}

	directoryPCMEMCALEtapPb 				= (TDirectory*)filePCMEMCALpPb->Get("Eta_pPb_5.023TeV_0-100%"); 
	histoPCMEMCALYieldEtapPb 				= (TH1D*)directoryPCMEMCALEtapPb->Get("CorrectedYieldNormEffEta");
	histoPCMEMCALRAWYieldEtapPb 			= (TH1D*)directoryPCMEMCALEtapPb->Get("RawYieldEta");
	graphPCMEMCALYieldEtaSysErrpPb			= (TGraphAsymmErrors*)directoryPCMEMCALEtapPb->Get("EtaSystError");	
	graphPCMEMCALYieldEtaSysErrRAApPb 		= (TGraphAsymmErrors*)directoryPCMEMCALEtapPb->Get("EtaSystErrorA");	
	histoPCMEMCALEtaPi0RatiopPb 			= (TH1D*)directoryPCMEMCALEtapPb->Get("EtatoPi0RatioNormEff");
	graphPCMEMCALEtaPi0RatioSysErrpPb 		= (TGraphAsymmErrors*)directoryPCMEMCALEtapPb->Get("EtatoPi0RatioSysNormEff"); 
	histoPCMEMCALMassEtaDatapPb 			= (TH1D*)directoryPCMEMCALEtapPb->Get("MassEta");
	histoPCMEMCALMassEtaMCpPb 				= (TH1D*)directoryPCMEMCALEtapPb->Get("TrueMassEta");
	histoPCMEMCALWidthEtaDatapPb 			= (TH1D*)directoryPCMEMCALEtapPb->Get("FWHMEtaMeV");
	histoPCMEMCALWidthEtaMCpPb 				= (TH1D*)directoryPCMEMCALEtapPb->Get("TrueFWHMEtaMeV");
	histoRatioPCMEMCALMassEtaDiffDataMC 	= (TH1D*)histoPCMEMCALMassEtaDatapPb->Clone("histoRatioPCMEMCALMassEtaDiffDataMC");
	histoRatioPCMEMCALMassEtaDiffDataMC->Add(histoPCMEMCALMassEtaMCpPb,-1);
	histoRatioPCMEMCALMassEtaDiffDataMC->Scale(1/(mesonMassExpectEta)*100);	

	histoPCMEMCALMassEtaDatapPb->Scale(1000.);
	histoPCMEMCALMassEtaMCpPb->Scale(1000.);
	
	//****************************************************************************************************
	//************************** Read data for PCM-PHOS *******************************************************
	//****************************************************************************************************
	cout << "PCM-PHOS" << endl;
	cout << "0-100%" << endl;
	filePCMPHOSpPb 							= new TFile(fileNamePCMPHOSpPb);
	directoryPCMPHOSPi0pPb 					= (TDirectory*)filePCMPHOSpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
	histoPCMPHOSYieldPi0pPb 				= (TH1D*)directoryPCMPHOSPi0pPb->Get("CorrectedYieldNormEffPi0");
	histoPCMPHOSRAWYieldPi0pPb 				= (TH1D*)directoryPCMPHOSPi0pPb->Get("RawYieldPi0PerEventsPi0");
	graphPCMPHOSYieldPi0SysErrpPb 			= (TGraphAsymmErrors*)directoryPCMPHOSPi0pPb->Get(nameGraphPCM.Data());	
	graphPCMPHOSYieldPi0SysErrRAApPb		= (TGraphAsymmErrors*)directoryPCMPHOSPi0pPb->Get("Pi0SystErrorA");	
	histoPCMPHOSMassPi0DatapPb 				= (TH1D*)directoryPCMPHOSPi0pPb->Get("MassPi0");
	histoPCMPHOSMassPi0MCpPb 				= (TH1D*)directoryPCMPHOSPi0pPb->Get("TrueMassPi0");
	histoPCMPHOSWidthPi0DatapPb 			= (TH1D*)directoryPCMPHOSPi0pPb->Get("FWHMPi0MeV");
	histoPCMPHOSWidthPi0MCpPb 				= (TH1D*)directoryPCMPHOSPi0pPb->Get("TrueFWHMPi0MeV");
	histoRatioPCMPHOSMassPi0DiffDataMC 		= (TH1D*)histoPCMPHOSMassPi0DatapPb->Clone("histoRatioPCMPHOSMassPi0DiffDataMC");
	histoRatioPCMPHOSMassPi0DiffDataMC->Add(histoPCMPHOSMassPi0MCpPb,-1);
	histoRatioPCMPHOSMassPi0DiffDataMC->Scale(1/(mesonMassExpectPi0)*100);	
	histoRatioPCMPHOSWidthPi0DiffDataMC 		= (TH1D*)histoPCMPHOSWidthPi0DatapPb->Clone("histoRatioPCMPHOSWidthPi0DiffDataMC");
	histoRatioPCMPHOSWidthPi0DiffDataMC->Add(histoPCMPHOSWidthPi0MCpPb,-1);
	histoRatioPCMPHOSWidthPi0DiffDataMC->Divide(histoRatioPCMPHOSWidthPi0DiffDataMC,histoPCMPHOSWidthPi0DatapPb,1.,1.,"");
	histoRatioPCMPHOSWidthPi0DiffDataMC->Scale(100.);
	histoPCMPHOSMassPi0DatapPb->Scale(1000.);
	histoPCMPHOSMassPi0MCpPb->Scale(1000.);

	// remove last points
	for (Int_t i = 0; i < 3; i++){
		histoPCMPHOSMassPi0DatapPb->SetBinContent(histoPCMPHOSMassPi0DatapPb->GetNbinsX()-i,0);
		histoPCMPHOSMassPi0MCpPb->SetBinContent(histoPCMPHOSMassPi0MCpPb->GetNbinsX()-i,0);
		histoPCMPHOSWidthPi0DatapPb->SetBinContent(histoPCMPHOSWidthPi0DatapPb->GetNbinsX()-i,10000.);
		histoPCMPHOSWidthPi0MCpPb->SetBinContent(histoPCMPHOSWidthPi0MCpPb->GetNbinsX()-i,10000.);
		histoRatioPCMPHOSMassPi0DiffDataMC->SetBinContent(histoRatioPCMPHOSMassPi0DiffDataMC->GetNbinsX()-i,-100);
		histoRatioPCMPHOSWidthPi0DiffDataMC->SetBinContent(histoRatioPCMPHOSMassPi0DiffDataMC->GetNbinsX()-i,-100);
	}

	directoryPCMPHOSEtapPb 					= (TDirectory*)filePCMPHOSpPb->Get("Eta_pPb_5.023TeV_0-100%"); 
	histoPCMPHOSYieldEtapPb 				= (TH1D*)directoryPCMPHOSEtapPb->Get("CorrectedYieldNormEffEta");
	histoPCMPHOSRAWYieldEtapPb 				= (TH1D*)directoryPCMPHOSEtapPb->Get("RawYieldEta");
	graphPCMPHOSYieldEtaSysErrpPb			= (TGraphAsymmErrors*)directoryPCMPHOSEtapPb->Get("EtaSystError");	
	graphPCMPHOSYieldEtaSysErrRAApPb 		= (TGraphAsymmErrors*)directoryPCMPHOSEtapPb->Get("EtaSystErrorA");	
	histoPCMPHOSEtaPi0RatiopPb 				= (TH1D*)directoryPCMPHOSEtapPb->Get("EtatoPi0RatioNormEff");
	graphPCMPHOSEtaPi0RatioSysErrpPb 		= (TGraphAsymmErrors*)directoryPCMPHOSEtapPb->Get("EtatoPi0RatioSysNormEff"); 
	histoPCMPHOSMassEtaDatapPb 				= (TH1D*)directoryPCMPHOSEtapPb->Get("MassEta");
	histoPCMPHOSMassEtaMCpPb 				= (TH1D*)directoryPCMPHOSEtapPb->Get("TrueMassEta");
	histoPCMPHOSWidthEtaDatapPb 			= (TH1D*)directoryPCMPHOSEtapPb->Get("FWHMEtaMeV");
	histoPCMPHOSWidthEtaMCpPb 				= (TH1D*)directoryPCMPHOSEtapPb->Get("TrueFWHMEtaMeV");
	histoRatioPCMPHOSMassEtaDiffDataMC 		= (TH1D*)histoPCMPHOSMassEtaDatapPb->Clone("histoRatioPCMPHOSMassEtaDiffDataMC");
	histoRatioPCMPHOSMassEtaDiffDataMC->Add(histoPCMPHOSMassEtaMCpPb,-1);
	histoRatioPCMPHOSMassEtaDiffDataMC->Scale(1/(mesonMassExpectEta)*100);	
	
	histoPCMPHOSMassEtaDatapPb->Scale(1000.);
	histoPCMPHOSMassEtaMCpPb->Scale(1000.);
	
	//****************************************************************************************************
	//************************** Read data for EMCAL-EMCAL ***********************************************
	//****************************************************************************************************
	cout << "EMCAL-EMCAL" << endl;
	cout << "0-100%" << endl;
	fileEMCALEMCALpPb 						= new TFile(fileNameEMCALEMCALpPb);
	directoryEMCALEMCALPi0pPb 				= (TDirectory*)fileEMCALEMCALpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
	histoEMCALEMCALYieldPi0pPb 				= (TH1D*)directoryEMCALEMCALPi0pPb->Get("CorrectedYieldNormEffPi0");
	histoEMCALEMCALRAWYieldPi0pPb 			= (TH1D*)directoryEMCALEMCALPi0pPb->Get("RawYieldPi0PerEventsPi0");
	graphEMCALEMCALYieldPi0SysErrpPb		= (TGraphAsymmErrors*)directoryEMCALEMCALPi0pPb->Get(nameGraphPCM.Data());	
	graphEMCALEMCALYieldPi0SysErrRAApPb		= (TGraphAsymmErrors*)directoryEMCALEMCALPi0pPb->Get("Pi0SystErrorA");	
	histoEMCALEMCALMassPi0DatapPb 			= (TH1D*)directoryEMCALEMCALPi0pPb->Get("MassPi0");
	histoEMCALEMCALMassPi0MCpPb 			= (TH1D*)directoryEMCALEMCALPi0pPb->Get("TrueMassPi0");
	histoEMCALEMCALWidthPi0DatapPb 			= (TH1D*)directoryEMCALEMCALPi0pPb->Get("FWHMPi0MeV");
	histoEMCALEMCALWidthPi0MCpPb 			= (TH1D*)directoryEMCALEMCALPi0pPb->Get("TrueFWHMPi0MeV");
	histoRatioEMCALEMCALMassPi0DiffDataMC 	= (TH1D*)histoEMCALEMCALMassPi0DatapPb->Clone("histoRatioEMCALEMCALMassPi0DiffDataMC");
	histoRatioEMCALEMCALMassPi0DiffDataMC->Add(histoEMCALEMCALMassPi0MCpPb,-1);
	histoRatioEMCALEMCALMassPi0DiffDataMC->Scale(1/(mesonMassExpectPi0)*100);
	histoRatioEMCALEMCALWidthPi0DiffDataMC 		= (TH1D*)histoEMCALEMCALWidthPi0DatapPb->Clone("histoRatioEMCALEMCALWidthPi0DiffDataMC");
	histoRatioEMCALEMCALWidthPi0DiffDataMC->Add(histoEMCALEMCALWidthPi0MCpPb,-1);
	histoRatioEMCALEMCALWidthPi0DiffDataMC->Divide(histoRatioEMCALEMCALWidthPi0DiffDataMC,histoEMCALEMCALWidthPi0DatapPb,1.,1.,"");
	histoRatioEMCALEMCALWidthPi0DiffDataMC->Scale(100.);
	histoEMCALEMCALMassPi0DatapPb->Scale(1000.);
	histoEMCALEMCALMassPi0MCpPb->Scale(1000.);
	// remove first points
	for (Int_t i = 0; i < 8; i++){
		histoEMCALEMCALMassPi0DatapPb->SetBinContent(i,0);
		histoEMCALEMCALMassPi0MCpPb->SetBinContent(i,0);
		histoEMCALEMCALWidthPi0DatapPb->SetBinContent(i,10000.);
		histoEMCALEMCALWidthPi0MCpPb->SetBinContent(i,10000.);
		histoRatioEMCALEMCALMassPi0DiffDataMC->SetBinContent(i,-100.);
		histoRatioEMCALEMCALWidthPi0DiffDataMC->SetBinContent(i,-100.);
	}

	directoryEMCALEMCALEtapPb 				= (TDirectory*)fileEMCALEMCALpPb->Get("Eta_pPb_5.023TeV_0-100%"); 
	histoEMCALEMCALYieldEtapPb 				= (TH1D*)directoryEMCALEMCALEtapPb->Get("CorrectedYieldNormEffEta");
	histoEMCALEMCALRAWYieldEtapPb 			= (TH1D*)directoryEMCALEMCALEtapPb->Get("RawYieldEta");
	graphEMCALEMCALYieldEtaSysErrpPb		= (TGraphAsymmErrors*)directoryEMCALEMCALEtapPb->Get("EtaSystError");	
	graphEMCALEMCALYieldEtaSysErrRAApPb 	= (TGraphAsymmErrors*)directoryEMCALEMCALEtapPb->Get("EtaSystErrorA");	
	histoEMCALEMCALEtaPi0RatiopPb 			= (TH1D*)directoryEMCALEMCALEtapPb->Get("EtatoPi0RatioNormEff");
	graphEMCALEMCALEtaPi0RatioSysErrpPb 	= (TGraphAsymmErrors*)directoryEMCALEMCALEtapPb->Get("EtatoPi0RatioSysNormEff"); 
	histoEMCALEMCALMassEtaDatapPb 			= (TH1D*)directoryEMCALEMCALEtapPb->Get("MassEta");
	histoEMCALEMCALMassEtaMCpPb 			= (TH1D*)directoryEMCALEMCALEtapPb->Get("TrueMassEta");
	histoEMCALEMCALWidthEtaDatapPb 			= (TH1D*)directoryEMCALEMCALEtapPb->Get("FWHMEtaMeV");
	histoEMCALEMCALWidthEtaMCpPb 			= (TH1D*)directoryEMCALEMCALEtapPb->Get("TrueFWHMEtaMeV");
	histoRatioEMCALEMCALMassEtaDiffDataMC 	= (TH1D*)histoEMCALEMCALMassEtaDatapPb->Clone("histoRatioEMCALEMCALMassEtaDiffDataMC");
	histoRatioEMCALEMCALMassEtaDiffDataMC->Add(histoEMCALEMCALMassEtaMCpPb,-1);
	histoRatioEMCALEMCALMassEtaDiffDataMC->Scale(1/(mesonMassExpectEta)*100);	
	histoEMCALEMCALMassEtaDatapPb->Scale(1000.);
	histoEMCALEMCALMassEtaMCpPb->Scale(1000.);

	//****************************************************************************************************
	//************************** Read data for PHOS-PHOS *************************************************
	//****************************************************************************************************
	cout << "PHOS-PHOS" << endl;
	cout << "0-100%" << endl;
	filePHOSPHOSpPb 						= new TFile(fileNamePHOSPHOSpPb);
	directoryPHOSPHOSPi0pPb 				= (TDirectory*)filePHOSPHOSpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
	histoPHOSPHOSYieldPi0pPb 				= (TH1D*)directoryPHOSPHOSPi0pPb->Get("CorrectedYieldNormEffPi0");
	histoPHOSPHOSRAWYieldPi0pPb 			= (TH1D*)directoryPHOSPHOSPi0pPb->Get("RawYieldPi0PerEventsPi0");
	graphPHOSPHOSYieldPi0SysErrpPb			= (TGraphAsymmErrors*)directoryPHOSPHOSPi0pPb->Get(nameGraphPCM.Data());	
	graphPHOSPHOSYieldPi0SysErrRAApPb		= (TGraphAsymmErrors*)directoryPHOSPHOSPi0pPb->Get("Pi0SystErrorA");	
	histoPHOSPHOSMassPi0DatapPb 			= (TH1D*)directoryPHOSPHOSPi0pPb->Get("MassPi0");
	histoPHOSPHOSMassPi0MCpPb 				= (TH1D*)directoryPHOSPHOSPi0pPb->Get("TrueMassPi0");
	histoPHOSPHOSWidthPi0DatapPb 			= (TH1D*)directoryPHOSPHOSPi0pPb->Get("FWHMPi0MeV");
	histoPHOSPHOSWidthPi0MCpPb 				= (TH1D*)directoryPHOSPHOSPi0pPb->Get("TrueFWHMPi0MeV");
	histoRatioPHOSPHOSMassPi0DiffDataMC 	= (TH1D*)histoPHOSPHOSMassPi0DatapPb->Clone("histoRatioPHOSPHOSMassPi0DiffDataMC");
	histoRatioPHOSPHOSMassPi0DiffDataMC->Add(histoPHOSPHOSMassPi0MCpPb,-1);
	histoRatioPHOSPHOSMassPi0DiffDataMC->Scale(1/(mesonMassExpectPi0)*100);
	histoRatioPHOSPHOSWidthPi0DiffDataMC 		= (TH1D*)histoPHOSPHOSWidthPi0DatapPb->Clone("histoRatioPHOSPHOSWidthPi0DiffDataMC");
	histoRatioPHOSPHOSWidthPi0DiffDataMC->Add(histoPHOSPHOSWidthPi0MCpPb,-1);
	histoRatioPHOSPHOSWidthPi0DiffDataMC->Divide(histoRatioPHOSPHOSWidthPi0DiffDataMC,histoPHOSPHOSWidthPi0DatapPb,1.,1.,"");
	histoRatioPHOSPHOSWidthPi0DiffDataMC->Scale(100.);
	histoPHOSPHOSMassPi0DatapPb->Scale(1000.);
	histoPHOSPHOSMassPi0MCpPb->Scale(1000.);
	// remove first points
// 	for (Int_t i = 0; i < 8; i++){
// 		histoPHOSPHOSMassPi0DatapPb->SetBinContent(i,0);
// 		histoPHOSPHOSMassPi0MCpPb->SetBinContent(i,0);
// 		histoPHOSPHOSWidthPi0DatapPb->SetBinContent(i,10000.);
// 		histoPHOSPHOSWidthPi0MCpPb->SetBinContent(i,10000.);
// 		histoRatioPHOSPHOSMassPi0DiffDataMC->SetBinContent(i,-100.);
// 	}

	
	//****************************************************************************************************
	//************************** Read data for Dalitz ****************************************************
	//****************************************************************************************************
	fileDalitzpPb 							= new TFile(fileNameDalitzpPb);
	directoryDalitzPi0pPb 					= (TDirectory*)fileDalitzpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
	TH1D* histoDalitzYieldPi0pPb 			= (TH1D*)directoryDalitzPi0pPb->Get(nameHistoPCM.Data());
	graphDalitzYieldPi0pPb 					= (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get(nameHistoDalitz.Data());
	graphDalitzYieldPi0SysErrpPb			= (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get(nameGraphPCM.Data()); 
	graphDalitzYieldPi0pPb->RemovePoint(0); 
	TGraphAsymmErrors*	graphDalitzYieldPi0SysErrRAApPb= 	(TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0SystErrorA");
	graphDalitzYieldPi0SysErrpPb->RemovePoint(0);
	histoDalitzMassPi0DatapPb 				= (TH1D*)directoryDalitzPi0pPb->Get("MassPi0");
	histoDalitzMassPi0MCpPb 				= (TH1D*)directoryDalitzPi0pPb->Get("TrueMassPi0");
	histoDalitzWidthPi0DatapPb 				= (TH1D*)directoryDalitzPi0pPb->Get("FWHMPi0MeV");
	histoDalitzWidthPi0MCpPb 				= (TH1D*)directoryDalitzPi0pPb->Get("TrueFWHMPi0MeV");
	histoRatioDalitzMassPi0DiffDataMC 		= (TH1D*)histoDalitzMassPi0DatapPb->Clone("histoRatioDalitzMassPi0DiffDataMC");
	histoRatioDalitzMassPi0DiffDataMC->Add(histoDalitzMassPi0MCpPb,-1);
	histoRatioDalitzMassPi0DiffDataMC->Scale(1/(mesonMassExpectPi0)*100);
	histoRatioDalitzWidthPi0DiffDataMC 		= (TH1D*)histoDalitzWidthPi0DatapPb->Clone("histoRatioDalitzWidthPi0DiffDataMC");
	histoRatioDalitzWidthPi0DiffDataMC->Add(histoDalitzWidthPi0MCpPb,-1);
	histoRatioDalitzWidthPi0DiffDataMC->Divide(histoRatioDalitzWidthPi0DiffDataMC,histoDalitzWidthPi0DatapPb,1.,1.,"");
	histoRatioDalitzWidthPi0DiffDataMC->Scale(100.);
	histoDalitzMassPi0DatapPb->Scale(1000.);
	histoDalitzMassPi0MCpPb->Scale(1000.);
	

	//****************************************************************************************************
	//************************** Read data for PHOS ******************************************************
	//****************************************************************************************************
	filePHOSpPb 							= new TFile(fileNameCaloPHOSpPb);
	directoryPHOSPi0pPb 					= (TDirectory*)filePHOSpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
	histoPHOSYieldPi0pPb 					= (TH1D*)directoryPHOSPi0pPb->Get(nameHistoPCM.Data());      		
	histoPHOSYieldPi0pPbSys 				= (TH1D*)directoryPHOSPi0pPb->Get("Pi0SystError");      
	histoPHOSMassPi0DatapPb 				= (TH1D*)directoryPHOSPi0pPb->Get("MassPi0");
	histoPHOSMassPi0MCpPb 					= (TH1D*)directoryPHOSPi0pPb->Get("TrueMassPi0");
	histoPHOSWidthPi0DatapPb 				= (TH1D*)directoryPHOSPi0pPb->Get("FWHMPi0");
	histoPHOSWidthPi0MCpPb 					= (TH1D*)directoryPHOSPi0pPb->Get("TrueFWHMPi0");
	histoRatioPHOSMassPi0DiffDataMC 		= (TH1D*)histoPHOSMassPi0DatapPb->Clone("histoRatioPHOSMassPi0DiffDataMC");
	histoRatioPHOSMassPi0DiffDataMC->Add(histoPHOSMassPi0MCpPb,-1);
	histoRatioPHOSMassPi0DiffDataMC->Scale(1/(mesonMassExpectPi0)*100);
	histoRatioPHOSWidthPi0DiffDataMC 		= (TH1D*)histoPHOSWidthPi0DatapPb->Clone("histoRatioPHOSWidthPi0DiffDataMC");
	histoRatioPHOSWidthPi0DiffDataMC->Add(histoPHOSWidthPi0MCpPb,-1);
	histoRatioPHOSWidthPi0DiffDataMC->Divide(histoRatioPHOSWidthPi0DiffDataMC,histoPHOSWidthPi0DatapPb,1.,1.,"");
	histoRatioPHOSWidthPi0DiffDataMC->Scale(100.);
	histoPHOSMassPi0DatapPb->Scale(1000.);
	histoPHOSMassPi0MCpPb->Scale(1000.);
	histoPHOSWidthPi0DatapPb->Scale(1000.);
	histoPHOSWidthPi0MCpPb->Scale(1000.);
	
	Int_t PHOS_start=histoPHOSYieldPi0pPbSys->FindBin(2.3);
	Int_t PHOS_stopp=histoPHOSYieldPi0pPbSys->FindBin(18.);
	Int_t PHOS_Bin=histoPHOSYieldPi0pPbSys->GetNbinsX();
	cout <<"PHOS range"  <<  PHOS_start << "   "<< PHOS_stopp ;
	for (int i=1; i<PHOS_start;i++ ){
		histoPHOSYieldPi0pPbSys->SetBinError(i,0.);
		histoPHOSYieldPi0pPbSys->SetBinContent(i,0.);
		histoPHOSYieldPi0pPb->SetBinError(i,0.);
		histoPHOSYieldPi0pPb->SetBinContent(i,0.);
	}
	for (int i=PHOS_stopp; i<=PHOS_Bin;i++ ){
		histoPHOSYieldPi0pPbSys->SetBinError(i,0.);
		histoPHOSYieldPi0pPbSys->SetBinContent(i,0.);
		histoPHOSYieldPi0pPb->SetBinError(i,0.);
		histoPHOSYieldPi0pPb->SetBinContent(i,0.);
	}  	

	//****************************************************************************************************
	//*************************** setting up width for plots *********************************************
	//****************************************************************************************************
	if (suffix.CompareTo("eps")==0){
		widthLinesBoxes				= 1.4;
		widthCommonFit					= 2.;
		widthStatErrBars				= 1.5;
		widthCommonErrors				= 1.1;
		widthCommonSpectrumBoxes			= 0.99;
	} else {
		widthLinesBoxes				= 2.3;
		widthCommonFit					= 2.6;
		widthStatErrBars				= 2.6;
		widthCommonErrors				= 2.;
		widthCommonSpectrumBoxes			= 2.3;
	}
	Size_t textSizeMassWidth = 0.039;
   //****************************************************************************************************
   //************************** Combine & Compare PCM & Dalitz minBias spectrum ***************************
   //****************************************************************************************************
   TGraphAsymmErrors* graphPCMYieldPi0SysErrpPbCopy1 = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrpPb->Clone("graphPCMYieldPi0SysErrpPbCopy1");
   TGraphAsymmErrors* graphPCMYieldPi0SysErrpPbCopy2 = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrpPb->Clone("graphPCMYieldPi0SysErrpPbCopy2");
   TGraphAsymmErrors* graphPCMYieldPi0SysErrpPbCopy3 = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrpPb->Clone("graphPCMYieldPi0SysErrpPbCopy3");
   TGraphAsymmErrors* graphPCMYieldPi0SysErrpPbCopy4 = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrpPb->Clone("graphPCMYieldPi0SysErrpPbCopy4");

   TGraphAsymmErrors* graphPCMYieldPi0SysErrApPbCopy1 = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrRAApPb->Clone("graphPCMYieldPi0SysErrApPbCopy1");
  TGraphAsymmErrors* graphDalitzYieldPi0SysErrApPbCopy1 = (TGraphAsymmErrors*) graphDalitzYieldPi0SysErrRAApPb->Clone("graphDalitzYieldPi0SysErrApPbCopy1");
 

 TGraphErrors* bla1 = NULL;
 TGraphErrors* bla2 = NULL;
 TGraphErrors* bla3 = NULL;
 TGraphErrors* bla4 = NULL;
   cout << "PCM Spectrum  - Dalitz" << endl;
   TGraphErrors* graphRatioPCMDalitz = CalculateRatioBetweenSpectraWithDifferentBinning(histoDalitzYieldPi0pPb, graphDalitzYieldPi0SysErrApPbCopy1,histoPCMYieldPi0pPb,graphPCMYieldPi0SysErrApPbCopy1,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
   graphRatioPCMDalitz->SetPointError(0,0.05,0.34515);

   graphRatioPCMDalitz->Print();
   
   cout << "PCM Spectrum  - PHOS " << endl;
   TGraphErrors* graphRatioPCMPHOS = CalculateRatioBetweenSpectraWithDifferentBinning(histoPHOSYieldPi0pPb,histoPHOSYieldPi0pPbSys,histoPCMYieldPi0pPb,graphPCMYieldPi0SysErrpPbCopy1,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
   graphRatioPCMPHOS->Print();
   cout << "PCM Spectrum  - PCM-PHOS " << endl;
   TGraphErrors* graphRatioPCM_PCMPHOS = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMPHOSYieldPi0pPb,graphPCMPHOSYieldPi0SysErrpPb,histoPCMYieldPi0pPb,graphPCMYieldPi0SysErrpPbCopy2,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
   graphRatioPCMPHOS->Print();                                                                                    
   cout << "PCM Spectrum  - EMCal " << endl;
   TGraphErrors* graphRatioPCMEMCal = CalculateRatioBetweenSpectraWithDifferentBinning(histoEMCALYieldPi0pPb,histoEMCALYieldPi0pPbSys,histoPCMYieldPi0pPb,graphPCMYieldPi0SysErrpPbCopy3,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
   graphRatioPCMPHOS->Print();
  cout << "PCM Spectrum  - PCM-EMCal " << endl;
  TGraphErrors* graphRatioPCM_PCMEMCal = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMEMCALYieldPi0pPb,graphPCMEMCALYieldPi0SysErrpPb,histoPCMYieldPi0pPb,graphPCMYieldPi0SysErrpPbCopy4,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4);
   graphRatioPCMPHOS->Print();

   //**************************************************************************************************************************************************
   //************************************* Plotting ratio of corrected spectra MinBias PCM & Dalitz Pi0 & Dalitz Pi0*************************************************
   //**************************************************************************************************************************************************
   
   TCanvas* canvasSpectraDiffDetMinBiasRatio = new TCanvas("canvasSpectraDiffDetMinBiasRatio","",200,10,1200,900);  // gives the page size
   DrawGammaCanvasSettings( canvasSpectraDiffDetMinBiasRatio,  0.1, 0.01, 0.015, 0.13);

   canvasSpectraDiffDetMinBiasRatio->SetLogx();
   TH2F * ratio2DYieldPi0pPbDiffDet;
   ratio2DYieldPi0pPbDiffDet = new TH2F("ratio2DYieldPi0pPbDiffDet","ratio2DYieldPi0pPbDiffDet",1000,0.23,30.,1000,0.25,2.65  );
   SetStyleHistoTH2ForGraphs(ratio2DYieldPi0pPbDiffDet, "#it{p}_{T} (GeV/#it{c})","Detector A/Detector B", 0.05,0.06, 0.05,0.06, 1.,0.8, 512, 505);
   ratio2DYieldPi0pPbDiffDet->DrawCopy(); 

   DrawGammaSetMarkerTGraphErr(graphRatioPCMDalitz, markerStyleDalitz, markerSizeInvYield, colorDalitz, colorDalitz);
   graphRatioPCMDalitz->Draw("p,same,e1");
   
   DrawGammaSetMarkerTGraphErr(graphRatioPCMPHOS , markerStylePHOS, markerSizeInvYield, colorPHOS, colorPHOS);
   graphRatioPCMPHOS->Draw("p,same,e1");
 
   DrawGammaSetMarkerTGraphErr(graphRatioPCMEMCal, markerStyleEMCAL, markerSizeInvYield, colorEMCAL, colorEMCAL);
   graphRatioPCMEMCal->Draw("p,same,e1");
  
   DrawGammaSetMarkerTGraphErr(graphRatioPCM_PCMPHOS, markerStylePCMPHOS, markerSizeInvYield, colorPCMPHOS, colorPCMPHOS);
   graphRatioPCM_PCMPHOS->Draw("p,same,e1");

   DrawGammaSetMarkerTGraphErr(graphRatioPCM_PCMEMCal, markerStylePCMEMCAL, markerSizeInvYield, colorPCMEMCAL, colorPCMEMCAL);
   graphRatioPCM_PCMEMCal->Draw("p,same,e1"); 


   TLatex *labelRatioPi0pPb1 = new TLatex(0.13,0.88,Form("#pi^{0},%s",collisionSystempPb.Data()));
   SetStyleTLatex( labelRatioPi0pPb1, 0.06,4);
   labelRatioPi0pPb1->Draw();
   
   DrawGammaLines(0., 30.,1., 1.,0.1);
      
   TLegend* legendSpectraDiffDetMinBiasRatio = new TLegend(0.65,0.7,0.89,0.95);
   legendSpectraDiffDetMinBiasRatio->SetFillColor(0);
   legendSpectraDiffDetMinBiasRatio->SetLineColor(0);
   legendSpectraDiffDetMinBiasRatio->SetTextSize(0.045);
   legendSpectraDiffDetMinBiasRatio->AddEntry(graphRatioPCMDalitz,Form("Dalitz/PCM"),"p");
   legendSpectraDiffDetMinBiasRatio->AddEntry(graphRatioPCMPHOS,Form("PHOS/PCM"),"p");
   legendSpectraDiffDetMinBiasRatio->AddEntry(graphRatioPCM_PCMPHOS,Form("PCM-PHOS/PCM"),"p");
   legendSpectraDiffDetMinBiasRatio->AddEntry(graphRatioPCMEMCal,Form("EMCal/PCM"),"p");
  legendSpectraDiffDetMinBiasRatio->AddEntry(graphRatioPCM_PCMEMCal,Form("PCM-EMCal/PCM"),"p");
   legendSpectraDiffDetMinBiasRatio->Draw();

   // TPaveText* pt= new TPaveText(0.75,0.82,0.89,0.9);
   // pt->ConvertNDCtoPad();
   // pt->AddText(Form("PCM spectra dated %s",datePCM.Data()));
   // pt->AddText(Form("Dalitz spectra dated %s",dateDalitz.Data()));
   // pt->AddText(Form("PHOS spectra dated %s",datePHOS.Data()));
   // pt->AddText(Form("EMCal spectra dated %s",dateEMCal.Data()));
   // pt->Draw();

   canvasSpectraDiffDetMinBiasRatio->Update();
   canvasSpectraDiffDetMinBiasRatio->Print(Form("%s/Pi0_Ratio_Spectra_MinBias_DiffDetectors.%s",outputDir.Data(),suffix.Data()));

	
	//****************************************************************************************************
	//******************************* Plotting FWHM for diff detectors ***********************************
	//****************************************************************************************************
	
	TCanvas* canvasFWHM = new TCanvas("canvasFWHM","",200,10,1300,1000);  // gives the page size
	DrawGammaCanvasSettings( canvasFWHM, 0.09, 0.02, 0.01, 0.09);   
	canvasFWHM->cd();
		
	canvasFWHM->SetLogx();
	
	TH2D *histo2DAllPi0FWHM;
	histo2DAllPi0FWHM = new TH2D("histo2DAllPi0FWHM", "histo2DAllPi0FWHM", 20,0.25,30. ,1000.,-30,60);
	SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})","peak width (MeV/#it{c}^{2})", 0.85*textSizeMassWidth,textSizeMassWidth, 0.85*textSizeMassWidth,textSizeMassWidth, 0.9,1.05, 510, 510);
	histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,25.5);
	histo2DAllPi0FWHM->GetXaxis()->SetLabelOffset(-0.01);
	histo2DAllPi0FWHM->GetYaxis()->SetLabelOffset(0.01);
	histo2DAllPi0FWHM->SetLineWidth(2);
	histo2DAllPi0FWHM->GetYaxis()->SetDecimals();
	histo2DAllPi0FWHM->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DAllPi0FWHM->DrawCopy();

	DrawGammaSetMarker(histoPHOSWidthPi0DatapPb, markerStylePHOS, markerSizeMass, colorPHOS, colorPHOS);
	histoPHOSWidthPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoPHOSWidthPi0MCpPb, markerStylePHOSMC, markerSizeMass, colorPHOSMC , colorPHOSMC);
	histoPHOSWidthPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	DrawGammaSetMarker(histoDalitzWidthPi0DatapPb, markerStyleDalitz, markerSizeMass, colorDalitz, colorDalitz);
	histoDalitzWidthPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoDalitzWidthPi0MCpPb, markerStyleDalitzMC, markerSizeMass, colorDalitzMC , colorDalitzMC);
	histoDalitzWidthPi0MCpPb->DrawCopy("same,e1,p,x0"); 
	
	DrawGammaSetMarker( histoEMCALWidthPi0DatapPb, markerStyleEMCAL, markerSizeMass*1.5, colorEMCAL, colorEMCAL);
	histoEMCALWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoEMCALWidthPi0MCpPb, markerStyleEMCALMC, markerSizeMass*1.5, colorEMCALMC, colorEMCALMC);
	histoEMCALWidthPi0MCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoPCMEMCALWidthPi0DatapPb, markerStylePCMEMCAL, markerSizeMass, colorPCMEMCAL, colorPCMEMCAL);
	histoPCMEMCALWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoPCMEMCALWidthPi0MCpPb, markerStylePCMEMCALMC, markerSizeMass, colorPCMEMCALMC, colorPCMEMCALMC);
	histoPCMEMCALWidthPi0MCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoPCMPHOSWidthPi0DatapPb, markerStylePCMPHOS, markerSizeMass, colorPCMPHOS, colorPCMPHOS);
	histoPCMPHOSWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoPCMPHOSWidthPi0MCpPb, markerStylePCMPHOSMC, markerSizeMass, colorPCMPHOSMC, colorPCMPHOSMC);
	histoPCMPHOSWidthPi0MCpPb->Draw("same,e1,p,x0"); 

// 	DrawGammaSetMarker( histoEMCALEMCALWidthPi0DatapPb, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
// 	histoEMCALEMCALWidthPi0DatapPb->Draw("same,e1,p,x0"); 
// 	DrawGammaSetMarker( histoEMCALEMCALWidthPi0MCpPb, markerStyleEMCALMC, markerSizeMass, colorEMCALMC+1, colorEMCALMC+1);
// 	histoEMCALEMCALWidthPi0MCpPb->Draw("same,e1,p,x0"); 
	
	DrawGammaSetMarker(histoPCMWidthPi0DatapPb, markerStylePCM, markerSizeMass, colorPCM, colorPCM);
	histoPCMWidthPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoPCMWidthPi0MCpPb, markerStylePCMMC, markerSizeMass, colorPCMMC, colorPCMMC);
	histoPCMWidthPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	//********************************** Defintion of the Legend **************************************************   
	Double_t columnsLegendFWHM2[4]    = {0.13,0.43,0.515,0.39};
	Double_t columnsLegendFWHM2Abs[4]    = {4,1.8,2.7,12};
	//    Double_t rowsLegendFWHM2[3]       = {0.66,0.33,0.0};
	Double_t rowsLegendFWHM2[7]       = {0.93,0.89,0.85,0.81, 0.77, 0.73, 0.69};
	Double_t rowsLegendFWHM2Abs[7]       = {0.2,23,21.8,20.6, 19.4, 18.2, 17.0};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnFWHM2 = textSizeMassWidth*0.85;
	Size_t textSizeTopRowFWHM2  = textSizeMassWidth*0.85; 
	//****************** Scale factors ******************
	Double_t scaleMarkerFWHM2      = 1.0;
	
	//    padFWHMLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textFWHM2PCM = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[1],"PCM (FWHM/2.35)");
	SetStyleTLatex( textFWHM2PCM, textSizeLeftColumnFWHM2,4);
	textFWHM2PCM->SetTextFont(42);
	textFWHM2PCM->Draw();
	TLatex *textFWHM2PCMEMCAL = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[2],"PCM-EMCal (FWHM/2.35)");
	SetStyleTLatex( textFWHM2PCMEMCAL, textSizeLeftColumnFWHM2,4);
	textFWHM2PCMEMCAL->SetTextFont(42);
	textFWHM2PCMEMCAL->Draw();
	TLatex *textFWHM2EMCAL = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[3],"EMCal (#sigma)");
	SetStyleTLatex( textFWHM2EMCAL, textSizeLeftColumnFWHM2,4);
	textFWHM2EMCAL->SetTextFont(42);
	textFWHM2EMCAL->Draw();

	TLatex *textFWHM2Dalitz = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[4],"Dalitz (FWHM/2.35)");
	SetStyleTLatex( textFWHM2Dalitz, textSizeLeftColumnFWHM2,4);
	textFWHM2Dalitz->SetTextFont(42);
	textFWHM2Dalitz->Draw();
	
// 	TLatex *textFWHM2EMCALEMCAL = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[4],"EMCal-EMCal (FWHM/2.35)");
// 	SetStyleTLatex( textFWHM2EMCALEMCAL, textSizeLeftColumnFWHM2,4);
// 	textFWHM2EMCALEMCAL->SetTextFont(42);
// 	textFWHM2EMCALEMCAL->Draw();
	TLatex *textFWHM2PHOS = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[5],"PHOS (#sigma)");
	SetStyleTLatex( textFWHM2PHOS, textSizeLeftColumnFWHM2,4);
	textFWHM2PHOS->SetTextFont(42);
	textFWHM2PHOS->Draw();
	TLatex *textFWHM2PCMPHOS = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[6],"PCM-PHOS (FWHM/2.35)");
	SetStyleTLatex( textFWHM2PCMPHOS, textSizeLeftColumnFWHM2,4);
	textFWHM2PCMPHOS->SetTextFont(42);
	textFWHM2PCMPHOS->Draw();

	
	//****************** second Column *************************************************
	TLatex *textFWHM2Data2 = new TLatex(columnsLegendFWHM2[1],rowsLegendFWHM2[0] ,"data");
	SetStyleTLatex( textFWHM2Data2, textSizeTopRowFWHM2 ,4);
	textFWHM2Data2->SetTextFont(42);
	textFWHM2Data2->Draw();
	TLatex *textFWHM2MC2 = new TLatex(columnsLegendFWHM2[2]
	,rowsLegendFWHM2[0],"MC");
	SetStyleTLatex( textFWHM2MC2, textSizeTopRowFWHM2,4);
	textFWHM2MC2->SetTextFont(42);
	textFWHM2MC2->Draw();
	
	TMarker* markerPCMPi0FWHM2 = CreateMarkerFromHisto(histoPCMWidthPi0DatapPb,columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[1] ,scaleMarkerFWHM2);
	markerPCMPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[1]);
	TMarker* markerPCMEMCALPi0FWHM2 = CreateMarkerFromHisto(histoPCMEMCALWidthPi0DatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[2],scaleMarkerFWHM2);
	markerPCMEMCALPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[2]);
	TMarker* markerEMCALPi0FWHM2 = CreateMarkerFromHisto(histoEMCALWidthPi0DatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[3],scaleMarkerFWHM2);
	markerEMCALPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[3]);
// 	TMarker* markerEMCALEMCALPi0FWHM2 = CreateMarkerFromHisto(histoEMCALEMCALWidthPi0DatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[4],scaleMarkerFWHM2);
// 	markerEMCALEMCALPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[4]);
	TMarker* markerDalitzPi0FWHM2 = CreateMarkerFromHisto(histoDalitzWidthPi0DatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[4],scaleMarkerFWHM2);
	markerDalitzPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[4]);
	TMarker* markerPHOSPi0FWHM2 = CreateMarkerFromHisto(histoPHOSWidthPi0DatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[5],scaleMarkerFWHM2);
	markerPHOSPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[5]);
	TMarker* markerPCMPHOSPi0FWHM2 = CreateMarkerFromHisto(histoPCMPHOSWidthPi0DatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[6],scaleMarkerFWHM2);
	markerPCMPHOSPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[6]);

	
	TMarker* markerPCMPi0FWHM2MC = CreateMarkerFromHisto(histoPCMWidthPi0MCpPb,columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[1],scaleMarkerFWHM2);
	markerPCMPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[1]);
	TMarker* markerPCMEMCALPi0FWHM2MC = CreateMarkerFromHisto(histoPCMEMCALWidthPi0MCpPb,columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[2] ,scaleMarkerFWHM2);
	markerPCMEMCALPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[2]);
	TMarker* markerEMCALPi0FWHM2MC = CreateMarkerFromHisto(histoEMCALWidthPi0MCpPb,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[3] ,scaleMarkerFWHM2);
	markerEMCALPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[3]);
// 	TMarker* markerEMCALEMCALPi0FWHM2MC = CreateMarkerFromHisto(histoEMCALEMCALWidthPi0MCpPb,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[4] ,scaleMarkerFWHM2);
// 	markerEMCALEMCALPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[4]);
	TMarker* markerDalitzPi0FWHM2MC = CreateMarkerFromHisto(histoDalitzWidthPi0MCpPb,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[4] ,scaleMarkerFWHM2);
	markerDalitzPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[4]);
	TMarker* markerPHOSPi0FWHM2MC = CreateMarkerFromHisto(histoPHOSWidthPi0MCpPb,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[5] ,scaleMarkerFWHM2);
	markerPHOSPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[5]);
	TMarker* markerPCMPHOSPi0FWHM2MC = CreateMarkerFromHisto(histoPCMPHOSWidthPi0MCpPb,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[6] ,scaleMarkerFWHM2);
	markerPCMPHOSPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[6]);

	
	TLatex *labelIndMeaspPbMassWidth = new TLatex(0.70,0.93,Form("%s",collisionSystempPb.Data()));
	SetStyleTLatex( labelIndMeaspPbMassWidth, 0.85*textSizeMassWidth,4);	
	labelIndMeaspPbMassWidth->Draw();
	TLatex *labelIndMeasProcMassWidth = new TLatex(0.83,0.89,Form("#pi^{0} #rightarrow #gamma#gamma"));
	SetStyleTLatex( labelIndMeasProcMassWidth, 0.85*textSizeMassWidth,4);	
	labelIndMeasProcMassWidth->Draw();

	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Pi0_WidthDifferentDetectors_pPb.%s",outputDir.Data(), suffix.Data()));

	//****************************************************************************************************
	//************************* Plotting FWHM for diff detectors only new EMCAL **************************
	//****************************************************************************************************

	canvasFWHM->cd();		
	canvasFWHM->SetLogx();
	
	histo2DAllPi0FWHM->DrawCopy();
	
	histoEMCALWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	histoEMCALWidthPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMEMCALWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPCMEMCALWidthPi0MCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoEMCALEMCALWidthPi0DatapPb, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
	histoEMCALEMCALWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoEMCALEMCALWidthPi0MCpPb, markerStyleEMCALMC, markerSizeMass, colorEMCALMC+1, colorEMCALMC+1);
	histoEMCALEMCALWidthPi0MCpPb->Draw("same,e1,p,x0"); 
	
	histoPCMWidthPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	histoPCMWidthPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textFWHM2PCM->Draw();
	textFWHM2PCMEMCAL->Draw();
	textFWHM2EMCAL->Draw();
	
	TLatex *textFWHM2EMCALEMCAL = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[4],"EMCal-EMCal (FWHM/2.35)");
	SetStyleTLatex( textFWHM2EMCALEMCAL, textSizeLeftColumnFWHM2,4);
	textFWHM2EMCALEMCAL->SetTextFont(42);
	textFWHM2EMCALEMCAL->Draw();
	
	//****************** second Column *************************************************
	textFWHM2Data2->Draw();
	textFWHM2MC2->Draw();
	
	markerPCMPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[1]);
	markerPCMEMCALPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[2]);
	markerEMCALPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[3]);
	TMarker* markerEMCALEMCALPi0FWHM2 = CreateMarkerFromHisto(histoEMCALEMCALWidthPi0DatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[4],scaleMarkerFWHM2);
	markerEMCALEMCALPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[4]);
	
	markerPCMPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[1]);
	markerPCMEMCALPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[2]);
	markerEMCALPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[3]);
	TMarker* markerEMCALEMCALPi0FWHM2MC = CreateMarkerFromHisto(histoEMCALEMCALWidthPi0MCpPb,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[4] ,scaleMarkerFWHM2);
	markerEMCALEMCALPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[4]);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Pi0_WidthDifferentDetectorsOnlyNewPlusRefEMCAL_pPb.%s",outputDir.Data(), suffix.Data()));

	canvasFWHM->cd();		
	canvasFWHM->SetLogx();
	
	histo2DAllPi0FWHM->DrawCopy();
	
	histoEMCALWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	histoEMCALWidthPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMEMCALWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPCMEMCALWidthPi0MCpPb->Draw("same,e1,p,x0"); 
	
	histoPCMWidthPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	histoPCMWidthPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textFWHM2PCM->Draw();
	textFWHM2PCMEMCAL->Draw();
	textFWHM2EMCAL->Draw();
		
	//****************** second Column *************************************************
	textFWHM2Data2->Draw();
	textFWHM2MC2->Draw();
	
	markerPCMPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[1]);
	markerPCMEMCALPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[2]);
	markerEMCALPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[3]);
	
	markerPCMPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[1]);
	markerPCMEMCALPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[2]);
	markerEMCALPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[3]);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Pi0_WidthDifferentDetectorsOnlyNewPlusRefEMCAL_stripped_pPb.%s",outputDir.Data(), suffix.Data()));
	
	//****************************************************************************************************
	//************************* Plotting FWHM for diff detectors only new PHOS ***************************
	//****************************************************************************************************

	canvasFWHM->cd();		
	canvasFWHM->SetLogx();
	
	histo2DAllPi0FWHM->DrawCopy();
	
	histoPHOSWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPHOSWidthPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMPHOSWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPCMPHOSWidthPi0MCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoPHOSPHOSWidthPi0DatapPb, markerStylePHOS, markerSizeMass*0.7, colorPHOS+1, colorPHOS+1);
	histoPHOSPHOSWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoPHOSPHOSWidthPi0MCpPb, markerStylePHOSMC, markerSizeMass*0.7, colorPHOSMC+1, colorPHOSMC+1);
	histoPHOSPHOSWidthPi0MCpPb->Draw("same,e1,p,x0"); 
	
	histoPCMWidthPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	histoPCMWidthPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textFWHM2PCM->Draw();

	TLatex *textFWHM2PCMPHOS2 = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[2],"PCM-PHOS (FWHM/2.35)");
	SetStyleTLatex( textFWHM2PCMPHOS2, textSizeLeftColumnFWHM2,4);
	textFWHM2PCMPHOS2->SetTextFont(42);
	textFWHM2PCMPHOS2->Draw();

	TLatex *textFWHM2PHOS2 = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[3],"PHOS (#sigma)");
	SetStyleTLatex( textFWHM2PHOS2, textSizeLeftColumnFWHM2,4);
	textFWHM2PHOS2->SetTextFont(42);
	textFWHM2PHOS2->Draw();
	
	TLatex *textFWHM2PHOSPHOS = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[4],"PHOS-PHOS (FWHM/2.35)");
	SetStyleTLatex( textFWHM2PHOSPHOS, textSizeLeftColumnFWHM2,4);
	textFWHM2PHOSPHOS->SetTextFont(42);
	textFWHM2PHOSPHOS->Draw();
	
	//****************** second Column *************************************************
	textFWHM2Data2->Draw();
	textFWHM2MC2->Draw();
	
	markerPCMPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[1]);
	markerPCMPHOSPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[2]);
	markerPHOSPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[3]);
	TMarker* markerPHOSPHOSPi0FWHM2 = CreateMarkerFromHisto(histoPHOSPHOSWidthPi0DatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[4],scaleMarkerFWHM2);
	markerPHOSPHOSPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[4]);
	
	markerPCMPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[1]);
	markerPCMPHOSPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[2]);
	markerPHOSPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[3]);
	TMarker* markerPHOSPHOSPi0FWHM2MC = CreateMarkerFromHisto(histoPHOSPHOSWidthPi0MCpPb,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[4] ,scaleMarkerFWHM2);
	markerPHOSPHOSPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[4]);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Pi0_WidthDifferentDetectorsOnlyNewPlusRefPHOS_pPb.%s",outputDir.Data(), suffix.Data()));

	canvasFWHM->cd();		
	canvasFWHM->SetLogx();
	
	histo2DAllPi0FWHM->DrawCopy();
	
	histoPHOSWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPHOSWidthPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMPHOSWidthPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPCMPHOSWidthPi0MCpPb->Draw("same,e1,p,x0"); 
	
	histoPCMWidthPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	histoPCMWidthPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textFWHM2PCM->Draw();
	textFWHM2PCMPHOS2->Draw();
	textFWHM2PHOS2->Draw();
	
	//****************** second Column *************************************************
	textFWHM2Data2->Draw();
	textFWHM2MC2->Draw();
	
	markerPCMPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[1]);
	markerPCMPHOSPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[2]);
	markerPHOSPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[3]);
	
	markerPCMPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[1]);
	markerPCMPHOSPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[2]);
	markerPHOSPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[3]);
	
	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Pi0_WidthDifferentDetectorsOnlyNewPlusRefPHOS_stripped_pPb.%s",outputDir.Data(), suffix.Data()));
	
	//****************************************************************************************************
	//******************************* Plotting Eta FWHM for diff detectors ***********************************
	//****************************************************************************************************
	
	canvasFWHM->cd();		
	canvasFWHM->SetLogx();
	
	TH2D *histo2DAllEtaFWHM;
	histo2DAllEtaFWHM = new TH2D("histo2DAllEtaFWHM", "histo2DAllEtaFWHM", 20,0.25,30. ,1000.,-30,100);
	SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})","peak width (MeV/#it{c}^{2})", 0.85*textSizeMassWidth,textSizeMassWidth, 0.85*textSizeMassWidth,textSizeMassWidth, 0.9,1.05, 510, 510);
	histo2DAllEtaFWHM->GetYaxis()->SetRangeUser(-1.,70.5);
	histo2DAllEtaFWHM->GetXaxis()->SetLabelOffset(-0.01);
	histo2DAllEtaFWHM->GetYaxis()->SetLabelOffset(0.01);
	histo2DAllEtaFWHM->SetLineWidth(2);
	histo2DAllEtaFWHM->GetYaxis()->SetDecimals();
	histo2DAllEtaFWHM->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DAllEtaFWHM->DrawCopy();
	
	DrawGammaSetMarker( histoPCMEMCALWidthEtaDatapPb, markerStylePCMEMCAL, markerSizeMass, colorPCMEMCAL, colorPCMEMCAL);
	histoPCMEMCALWidthEtaDatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoPCMEMCALWidthEtaMCpPb, markerStylePCMEMCALMC, markerSizeMass, colorPCMEMCALMC, colorPCMEMCALMC);
	histoPCMEMCALWidthEtaMCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoPCMPHOSWidthEtaDatapPb, markerStylePCMPHOS, markerSizeMass, colorPCMPHOS, colorPCMPHOS);
	histoPCMPHOSWidthEtaDatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoPCMPHOSWidthEtaMCpPb, markerStylePCMPHOSMC, markerSizeMass, colorPCMPHOSMC, colorPCMPHOSMC);
	histoPCMPHOSWidthEtaMCpPb->Draw("same,e1,p,x0"); 
	
	DrawGammaSetMarker( histoEMCALEMCALWidthEtaDatapPb, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
	histoEMCALEMCALWidthEtaDatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoEMCALEMCALWidthEtaMCpPb, markerStyleEMCALMC, markerSizeMass, colorEMCALMC+1, colorEMCALMC+1);
	histoEMCALEMCALWidthEtaMCpPb->Draw("same,e1,p,x0"); 
	
	DrawGammaSetMarker(histoPCMWidthEtaDatapPb, markerStylePCM, markerSizeMass, colorPCM, colorPCM);
	histoPCMWidthEtaDatapPb->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoPCMWidthEtaMCpPb, markerStylePCMMC, markerSizeMass, colorPCMMC, colorPCMMC);
	histoPCMWidthEtaMCpPb->DrawCopy("same,e1,p,x0"); 

	Double_t rowsLegendFWHM2AbsEta[6]       = {0.2,63.5, 60.5, 57, 54, 18.2};

	// 	//****************** first Column **************************************************
	textFWHM2PCM->Draw();
	textFWHM2PCMEMCAL->Draw();
	TLatex *textFWHM2EMCALEMCALEta = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[3],"EMCal-EMCal (FWHM/2.35)");
	SetStyleTLatex( textFWHM2EMCALEMCALEta, textSizeLeftColumnFWHM2,4);
	textFWHM2EMCALEMCALEta->SetTextFont(42);
	textFWHM2EMCALEMCALEta->Draw();
	TLatex *textFWHM2PCMPHOSEta = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[4],"PCM-PHOS (FWHM/2.35)");
	SetStyleTLatex( textFWHM2PCMPHOSEta, textSizeLeftColumnFWHM2,4);
	textFWHM2PCMPHOSEta->SetTextFont(42);
	textFWHM2PCMPHOSEta->Draw();

// 	//****************** second Column *************************************************
	textFWHM2Data2->Draw();
	textFWHM2MC2->Draw();
	TMarker* markerPCMEtaFWHM2 = CreateMarkerFromHisto(histoPCMWidthEtaDatapPb,columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2AbsEta[1] ,scaleMarkerFWHM2);
	markerPCMEtaFWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2AbsEta[1]);
	TMarker* markerPCMEMCALEtaFWHM2 = CreateMarkerFromHisto(histoPCMEMCALWidthEtaDatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2AbsEta[2],scaleMarkerFWHM2);
	markerPCMEMCALEtaFWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2AbsEta[2]);
	TMarker* markerEMCALEMCALEtaFWHM2 = CreateMarkerFromHisto(histoEMCALEMCALWidthEtaDatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2AbsEta[3],scaleMarkerFWHM2);
	markerEMCALEMCALEtaFWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2AbsEta[3]);
	TMarker* markerPCMPHOSEtaFWHM2 = CreateMarkerFromHisto(histoPCMPHOSWidthEtaDatapPb,columnsLegendFWHM2Abs[1],rowsLegendFWHM2AbsEta[4],scaleMarkerFWHM2);
	markerPCMPHOSEtaFWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2AbsEta[4]);
	
	TMarker* markerPCMEtaFWHM2MC = CreateMarkerFromHisto(histoPCMWidthEtaMCpPb,columnsLegendFWHM2Abs[2],rowsLegendFWHM2AbsEta[1],scaleMarkerFWHM2);
	markerPCMEtaFWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2AbsEta[1]);
	TMarker* markerPCMEMCALEtaFWHM2MC = CreateMarkerFromHisto(histoPCMEMCALWidthEtaMCpPb,columnsLegendFWHM2Abs[2],rowsLegendFWHM2AbsEta[2] ,scaleMarkerFWHM2);
	markerPCMEMCALEtaFWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2AbsEta[2]);
	TMarker* markerEMCALEMCALEtaFWHM2MC = CreateMarkerFromHisto(histoEMCALEMCALWidthEtaMCpPb,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2AbsEta[3] ,scaleMarkerFWHM2);
	markerEMCALEMCALEtaFWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2AbsEta[3]);
	TMarker* markerPCMPHOSEtaFWHM2MC = CreateMarkerFromHisto(histoPCMPHOSWidthEtaMCpPb,columnsLegendFWHM2Abs[2],rowsLegendFWHM2AbsEta[4] ,scaleMarkerFWHM2);
	markerPCMPHOSEtaFWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2AbsEta[4]);
	
	labelIndMeaspPbMassWidth->Draw();
	TLatex *labelIndMeasProcEtaMassWidth = new TLatex(0.83,0.89,Form("#eta #rightarrow #gamma#gamma"));
	SetStyleTLatex( labelIndMeasProcEtaMassWidth, 0.85*textSizeMassWidth,4);	
	labelIndMeasProcEtaMassWidth->Draw();
	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Eta_WidthDifferentDetectors_pPb.%s",outputDir.Data(), suffix.Data()));
	
	
	//****************************************************************************************************
	//******************************* Plotting Mass for diff detectors ***********************************
	//****************************************************************************************************
	
	TCanvas* canvasMass = new TCanvas("canvasMass","",200,10,1300,1000);  // gives the page size
	DrawGammaCanvasSettings( canvasMass, 0.09, 0.02, 0.01, 0.09);   
	canvasMass->cd();
		
	canvasMass->SetLogx();
	
	TH2D *histo2DAllPi0Mass;
	histo2DAllPi0Mass = new TH2D("histo2DAllPi0Mass", "histo2DAllPi0Mass", 20,0.25,30. ,1000.,125.,170);

	SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})","peak position (MeV/#it{c}^{2})", 0.85*textSizeMassWidth,textSizeMassWidth, 0.85*textSizeMassWidth,textSizeMassWidth, 0.9,1.05, 515, 510);
	histo2DAllPi0Mass->GetYaxis()->SetRangeUser(129.,149);
	histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.01);
	histo2DAllPi0Mass->GetYaxis()->SetLabelOffset(0.01);
	histo2DAllPi0Mass->SetLineWidth(2);
	histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DAllPi0Mass->GetYaxis()->SetDecimals();
	histo2DAllPi0Mass->DrawCopy();

	DrawGammaSetMarker(histoPHOSMassPi0DatapPb, markerStylePHOS, markerSizeMass, colorPHOS, colorPHOS);
	histoPHOSMassPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoPHOSMassPi0MCpPb, markerStylePHOSMC, markerSizeMass, colorPHOSMC , colorPHOSMC);
	histoPHOSMassPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	DrawGammaSetMarker(histoDalitzMassPi0DatapPb, markerStyleDalitz, markerSizeMass, colorDalitz, colorDalitz);
	histoDalitzMassPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoDalitzMassPi0MCpPb, markerStyleDalitzMC, markerSizeMass, colorDalitzMC , colorDalitzMC);
	histoDalitzMassPi0MCpPb->DrawCopy("same,e1,p,x0"); 
	
	DrawGammaSetMarker( histoEMCALMassPi0DatapPb, markerStyleEMCAL, markerSizeMass*1.5, colorEMCAL, colorEMCAL);
	histoEMCALMassPi0DatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoEMCALMassPi0MCpPb, markerStyleEMCALMC, markerSizeMass*1.5, colorEMCALMC, colorEMCALMC);
	histoEMCALMassPi0MCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoPCMEMCALMassPi0DatapPb, markerStylePCMEMCAL, markerSizeMass, colorPCMEMCAL, colorPCMEMCAL);
	histoPCMEMCALMassPi0DatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoPCMEMCALMassPi0MCpPb, markerStylePCMEMCALMC, markerSizeMass, colorPCMEMCALMC, colorPCMEMCALMC);
	histoPCMEMCALMassPi0MCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoPCMPHOSMassPi0DatapPb, markerStylePCMPHOS, markerSizeMass, colorPCMPHOS, colorPCMPHOS);
	histoPCMPHOSMassPi0DatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoPCMPHOSMassPi0MCpPb, markerStylePCMPHOSMC, markerSizeMass, colorPCMPHOSMC, colorPCMPHOSMC);
	histoPCMPHOSMassPi0MCpPb->Draw("same,e1,p,x0"); 

	
// 	DrawGammaSetMarker( histoEMCALEMCALMassPi0DatapPb, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
// 	histoEMCALEMCALMassPi0DatapPb->Draw("same,e1,p,x0"); 
// 	DrawGammaSetMarker( histoEMCALEMCALMassPi0MCpPb, markerStyleEMCALMC, markerSizeMass, colorEMCALMC+1, colorEMCALMC+1);
// 	histoEMCALEMCALMassPi0MCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker(histoPCMMassPi0DatapPb, markerStylePCM, markerSizeMass, colorPCM, colorPCM);
	histoPCMMassPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoPCMMassPi0MCpPb, markerStylePCMMC, markerSizeMass, colorPCMMC, colorPCMMC);
	histoPCMMassPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	//********************************** Defintion of the Legend **************************************************   
	Double_t columnsLegendMass2[4]    = {0.25,0.415,0.505,0.39};
	Double_t columnsLegendMass2Abs[4]    = {4,1.6,2.55,12};
	//    Double_t rowsLegendMass2[3]       = {0.66,0.33,0.0};
	Double_t rowsLegendMass2[7]       = {0.93,0.89,0.85,0.81, 0.77,0.73, 0.69};
	Double_t rowsLegendMass2Abs[7]       = {0.2,147.0,146.2,145.3, 144.4, 143.5, 142.6};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnMass2 = textSizeMassWidth*0.85;
	Size_t textSizeTopRowMass2  = textSizeMassWidth*0.85; 
	//****************** Scale factors ******************
	Double_t scaleMarkerMass2      = 1.1;
	
	//    padMassLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textMass2PCM = new TLatex(columnsLegendMass2[0],rowsLegendMass2[1],"PCM");
	SetStyleTLatex( textMass2PCM, textSizeLeftColumnMass2,4);
	textMass2PCM->SetTextFont(42);
	textMass2PCM->Draw();
	TLatex *textMass2PCMEMCAL = new TLatex(columnsLegendMass2[0],rowsLegendMass2[2],"PCM-EMCal");
	SetStyleTLatex( textMass2PCMEMCAL, textSizeLeftColumnMass2,4);
	textMass2PCMEMCAL->SetTextFont(42);
	textMass2PCMEMCAL->Draw();
	TLatex *textMass2EMCAL = new TLatex(columnsLegendMass2[0],rowsLegendMass2[3],"EMCal");
	SetStyleTLatex( textMass2EMCAL, textSizeLeftColumnMass2,4);
	textMass2EMCAL->SetTextFont(42);
	textMass2EMCAL->Draw();
// 	TLatex *textMass2EMCALEMCAL = new TLatex(columnsLegendMass2[0],rowsLegendMass2[4],"EMCal-EMCal");
// 	SetStyleTLatex( textMass2EMCALEMCAL, textSizeLeftColumnMass2,4);
// 	textMass2EMCALEMCAL->SetTextFont(42);
// 	textMass2EMCALEMCAL->Draw();
	TLatex *textMass2Dalitz = new TLatex(columnsLegendMass2[0],rowsLegendMass2[4],"Dalitz");
	SetStyleTLatex( textMass2Dalitz, textSizeLeftColumnMass2,4);
	textMass2Dalitz->SetTextFont(42);
	textMass2Dalitz->Draw();
	TLatex *textMass2PHOS = new TLatex(columnsLegendMass2[0],rowsLegendMass2[5],"PHOS");
	SetStyleTLatex( textMass2PHOS, textSizeLeftColumnMass2,4);
	textMass2PHOS->SetTextFont(42);
	textMass2PHOS->Draw();
	TLatex *textMass2PCMPHOS = new TLatex(columnsLegendMass2[0],rowsLegendMass2[6],"PCM-PHOS");
	SetStyleTLatex( textMass2PCMPHOS, textSizeLeftColumnMass2,4);
	textMass2PCMPHOS->SetTextFont(42);
	textMass2PCMPHOS->Draw();

	
	//****************** second Column *************************************************
	TLatex *textMass2Data2 = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"data");
	SetStyleTLatex( textMass2Data2, textSizeTopRowMass2 ,4);
	textMass2Data2->SetTextFont(42);
	textMass2Data2->Draw();
	TLatex *textMass2MC2 = new TLatex(columnsLegendMass2[2]
	,rowsLegendMass2[0],"MC");
	SetStyleTLatex( textMass2MC2, textSizeTopRowMass2,4);
	textMass2MC2->SetTextFont(42);
	textMass2MC2->Draw();
	
	TMarker* markerPCMPi0Mass2 = CreateMarkerFromHisto(histoPCMMassPi0DatapPb,columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[1] ,scaleMarkerMass2);
	markerPCMPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[1]);
	TMarker* markerPCMEMCALPi0Mass2 = CreateMarkerFromHisto(histoPCMEMCALMassPi0DatapPb,columnsLegendMass2Abs[1],rowsLegendMass2Abs[2],scaleMarkerMass2);
	markerPCMEMCALPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[2]);
	TMarker* markerEMCALPi0Mass2 = CreateMarkerFromHisto(histoEMCALMassPi0DatapPb,columnsLegendMass2Abs[1],rowsLegendMass2Abs[3],scaleMarkerMass2);
	markerEMCALPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[3]);
// 	TMarker* markerEMCALEMCALPi0Mass2 = CreateMarkerFromHisto(histoEMCALEMCALMassPi0DatapPb,columnsLegendMass2Abs[1],rowsLegendMass2Abs[4],scaleMarkerMass2);
// 	markerEMCALEMCALPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[4]);
	TMarker* markerDalitzPi0Mass2 = CreateMarkerFromHisto(histoDalitzMassPi0DatapPb,columnsLegendMass2Abs[1],rowsLegendMass2Abs[4],scaleMarkerMass2);
	markerDalitzPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[4]);
	TMarker* markerPHOSPi0Mass2 = CreateMarkerFromHisto(histoPHOSMassPi0DatapPb,columnsLegendMass2Abs[1],rowsLegendMass2Abs[5],scaleMarkerMass2);
	markerPHOSPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[5]);
	TMarker* markerPCMPHOSPi0Mass2 = CreateMarkerFromHisto(histoPCMPHOSMassPi0DatapPb,columnsLegendMass2Abs[1],rowsLegendMass2Abs[5],scaleMarkerMass2);
	markerPCMPHOSPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[6]);
	
	TMarker* markerPCMPi0Mass2MC = CreateMarkerFromHisto(histoPCMMassPi0MCpPb,columnsLegendMass2Abs[2],rowsLegendMass2Abs[1],scaleMarkerMass2);
	markerPCMPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[1]);
	TMarker* markerPCMEMCALPi0Mass2MC = CreateMarkerFromHisto(histoPCMEMCALMassPi0MCpPb,columnsLegendMass2Abs[2],rowsLegendMass2Abs[2] ,scaleMarkerMass2);
	markerPCMEMCALPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[2]);
	TMarker* markerEMCALPi0Mass2MC = CreateMarkerFromHisto(histoEMCALMassPi0MCpPb,columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[3] ,scaleMarkerMass2);
	markerEMCALPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[3]);
// 	TMarker* markerEMCALEMCALPi0Mass2MC = CreateMarkerFromHisto(histoEMCALEMCALMassPi0MCpPb,columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[4] ,scaleMarkerMass2);
// 	markerEMCALEMCALPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[4]);
	TMarker* markerDalitzPi0Mass2MC = CreateMarkerFromHisto(histoDalitzMassPi0MCpPb,columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[4] ,scaleMarkerMass2);
	markerDalitzPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[4]);
	TMarker* markerPHOSPi0Mass2MC = CreateMarkerFromHisto(histoPHOSMassPi0MCpPb,columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[5] ,scaleMarkerMass2);
	markerPHOSPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[5]);
	TMarker* markerPCMPHOSPi0Mass2MC = CreateMarkerFromHisto(histoPCMPHOSMassPi0MCpPb,columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[5] ,scaleMarkerMass2);
	markerPCMPHOSPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[6]);

	DrawGammaLines(0.25, 30. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Pi0_MassDifferentDetectors_pPb.%s",outputDir.Data(), suffix.Data()));

	//****************************************************************************************************
	//*********************** Plotting Mass for diff detectors only new EMCAL ****************************
	//****************************************************************************************************

	canvasMass->cd();
	canvasMass->SetLogx();
	
	histo2DAllPi0Mass->DrawCopy();

	histoEMCALMassPi0DatapPb->Draw("same,e1,p,x0"); 
	histoEMCALMassPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMEMCALMassPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPCMEMCALMassPi0MCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoEMCALEMCALMassPi0DatapPb, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
	histoEMCALEMCALMassPi0DatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoEMCALEMCALMassPi0MCpPb, markerStyleEMCALMC, markerSizeMass, colorEMCALMC+1, colorEMCALMC+1);
	histoEMCALEMCALMassPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMMassPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	histoPCMMassPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textMass2PCM->Draw();
	textMass2PCMEMCAL->Draw();
	textMass2EMCAL->Draw();
	TLatex *textMass2EMCALEMCAL = new TLatex(columnsLegendMass2[0],rowsLegendMass2[4],"EMCal-EMCal");
	SetStyleTLatex( textMass2EMCALEMCAL, textSizeLeftColumnMass2,4);
	textMass2EMCALEMCAL->SetTextFont(42);
	textMass2EMCALEMCAL->Draw();
	
	//****************** second Column *************************************************
	textMass2Data2->Draw();
	textMass2MC2->Draw();
	
	markerPCMPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[1]);
	markerPCMEMCALPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[2]);
	markerEMCALPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[3]);
	TMarker* markerEMCALEMCALPi0Mass2 = CreateMarkerFromHisto(histoEMCALEMCALMassPi0DatapPb,columnsLegendMass2Abs[1],rowsLegendMass2Abs[4],scaleMarkerMass2);
	markerEMCALEMCALPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[4]);
	
	markerPCMPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[1]);
	markerPCMEMCALPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[2]);
	markerEMCALPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[3]);
	TMarker* markerEMCALEMCALPi0Mass2MC = CreateMarkerFromHisto(histoEMCALEMCALMassPi0MCpPb,columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[4] ,scaleMarkerMass2);
	markerEMCALEMCALPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[4]);

	DrawGammaLines(0.25, 30. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Pi0_MassDifferentDetectorsOnlyNewPlusRefEMCAL_pPb.%s",outputDir.Data(), suffix.Data()));

	canvasMass->cd();
	canvasMass->SetLogx();
	
	histo2DAllPi0Mass->DrawCopy();

	histoEMCALMassPi0DatapPb->Draw("same,e1,p,x0"); 
	histoEMCALMassPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMEMCALMassPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPCMEMCALMassPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMMassPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	histoPCMMassPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textMass2PCM->Draw();
	textMass2PCMEMCAL->Draw();
	textMass2EMCAL->Draw();
	
	//****************** second Column *************************************************
	textMass2Data2->Draw();
	textMass2MC2->Draw();
	
	markerPCMPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[1]);
	markerPCMEMCALPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[2]);
	markerEMCALPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[3]);
	
	markerPCMPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[1]);
	markerPCMEMCALPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[2]);
	markerEMCALPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[3]);

	DrawGammaLines(0.25, 30. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Pi0_MassDifferentDetectorsOnlyNewPlusRefEMCAL_stripped_pPb.%s",outputDir.Data(), suffix.Data()));
	
	//****************************************************************************************************
	//*********************** Plotting Mass for diff detectors only new PHOS *****************************
	//****************************************************************************************************

	canvasMass->cd();
	canvasMass->SetLogx();
	
	histo2DAllPi0Mass->DrawCopy();

	histoPHOSMassPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPHOSMassPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMPHOSMassPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPCMPHOSMassPi0MCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoPHOSPHOSMassPi0DatapPb, markerStylePHOS, markerSizeMass*0.7, colorPHOS+1, colorPHOS+1);
	histoPHOSPHOSMassPi0DatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoPHOSPHOSMassPi0MCpPb, markerStylePHOSMC, markerSizeMass*0.7, colorPHOSMC+1, colorPHOSMC+1);
	histoPHOSPHOSMassPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMMassPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	histoPCMMassPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textMass2PCM->Draw();

	TLatex *textMass2PCMPHOS2 = new TLatex(columnsLegendMass2[0],rowsLegendMass2[2],"PCM-PHOS");
	SetStyleTLatex( textMass2PCMPHOS2, textSizeLeftColumnMass2,4);
	textMass2PCMPHOS2->SetTextFont(42);
	textMass2PCMPHOS2->Draw();
	TLatex *textMass2PHOS2 = new TLatex(columnsLegendMass2[0],rowsLegendMass2[3],"PHOS");
	SetStyleTLatex( textMass2PHOS2, textSizeLeftColumnMass2,4);
	textMass2PHOS2->SetTextFont(42);
	textMass2PHOS2->Draw();
	TLatex *textMass2PHOSPHOS = new TLatex(columnsLegendMass2[0],rowsLegendMass2[4],"PHOS-PHOS");
	SetStyleTLatex( textMass2PHOSPHOS, textSizeLeftColumnMass2,4);
	textMass2PHOSPHOS->SetTextFont(42);
	textMass2PHOSPHOS->Draw();
	
	//****************** second Column *************************************************
	textMass2Data2->Draw();
	textMass2MC2->Draw();
	
	markerPCMPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[1]);
	markerPCMPHOSPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[2]);
	markerPHOSPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[3]);
	TMarker* markerPHOSPHOSPi0Mass2 = CreateMarkerFromHisto(histoPHOSPHOSMassPi0DatapPb,columnsLegendMass2Abs[1],rowsLegendMass2Abs[4],scaleMarkerMass2);
	markerPHOSPHOSPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[4]);
	
	markerPCMPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[1]);
	markerPCMPHOSPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[2]);
	markerPHOSPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[3]);
	TMarker* markerPHOSPHOSPi0Mass2MC = CreateMarkerFromHisto(histoPHOSPHOSMassPi0MCpPb,columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[4] ,scaleMarkerMass2);
	markerPHOSPHOSPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[4]);

	DrawGammaLines(0.25, 30. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Pi0_MassDifferentDetectorsOnlyNewPlusRefPHOS_pPb.%s",outputDir.Data(), suffix.Data()));

	canvasMass->cd();
	canvasMass->SetLogx();
	
	histo2DAllPi0Mass->DrawCopy();

	histoPHOSMassPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPHOSMassPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMPHOSMassPi0DatapPb->Draw("same,e1,p,x0"); 
	histoPCMPHOSMassPi0MCpPb->Draw("same,e1,p,x0"); 

	histoPCMMassPi0DatapPb->DrawCopy("same,e1,p,x0"); 
	histoPCMMassPi0MCpPb->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textMass2PCM->Draw();
	textMass2PCMPHOS2->Draw();
	textMass2PHOS2->Draw();
	
	//****************** second Column *************************************************
	textMass2Data2->Draw();
	textMass2MC2->Draw();
	
	markerPCMPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[1]);
	markerPCMPHOSPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[2]);
	markerPHOSPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2Abs[3]);
	
	markerPCMPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2] ,rowsLegendMass2Abs[1]);
	markerPCMPHOSPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[2]);
	markerPHOSPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2Abs[3]);

	DrawGammaLines(0.25, 30. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Pi0_MassDifferentDetectorsOnlyNewPlusRefPHOS_stripped_pPb.%s",outputDir.Data(), suffix.Data()));

	//****************************************************************************************************
	//******************************* Plotting Eta Mass for diff detectors *******************************
	//****************************************************************************************************
	
	canvasMass->cd();		
	canvasMass->SetLogx();
	
	TH2D *histo2DAllEtaMass;
	histo2DAllEtaMass = new TH2D("histo2DAllEtaMass", "histo2DAllEtaMass", 20,0.25,30. ,1000.,500.,600.);

	SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})","peak position (MeV/#it{c}^{2})", 0.85*textSizeMassWidth,textSizeMassWidth, 0.85*textSizeMassWidth,textSizeMassWidth, 0.9,1.05, 515, 510);
	histo2DAllEtaMass->GetYaxis()->SetRangeUser(520.,572.);
	histo2DAllEtaMass->GetXaxis()->SetLabelOffset(-0.01);
	histo2DAllEtaMass->GetYaxis()->SetLabelOffset(0.01);
	histo2DAllEtaMass->SetLineWidth(2);
	histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DAllEtaMass->GetYaxis()->SetDecimals();
	histo2DAllEtaMass->DrawCopy();
	
	DrawGammaSetMarker( histoPCMEMCALMassEtaDatapPb, markerStylePCMEMCAL, markerSizeMass, colorPCMEMCAL, colorPCMEMCAL);
	histoPCMEMCALMassEtaDatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoPCMEMCALMassEtaMCpPb, markerStylePCMEMCALMC, markerSizeMass, colorPCMEMCALMC, colorPCMEMCALMC);
	histoPCMEMCALMassEtaMCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoPCMPHOSMassEtaDatapPb, markerStylePCMPHOS, markerSizeMass, colorPCMPHOS, colorPCMPHOS);
	histoPCMPHOSMassEtaDatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoPCMPHOSMassEtaMCpPb, markerStylePCMPHOSMC, markerSizeMass, colorPCMPHOSMC, colorPCMPHOSMC);
	histoPCMPHOSMassEtaMCpPb->Draw("same,e1,p,x0"); 

	
	DrawGammaSetMarker( histoEMCALEMCALMassEtaDatapPb, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
	histoEMCALEMCALMassEtaDatapPb->Draw("same,e1,p,x0"); 
	DrawGammaSetMarker( histoEMCALEMCALMassEtaMCpPb, markerStyleEMCALMC, markerSizeMass, colorEMCALMC+1, colorEMCALMC+1);
	histoEMCALEMCALMassEtaMCpPb->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker(histoPCMMassEtaDatapPb, markerStylePCM, markerSizeMass, colorPCM, colorPCM);
	histoPCMMassEtaDatapPb->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoPCMMassEtaMCpPb, markerStylePCMMC, markerSizeMass, colorPCMMC, colorPCMMC);
	histoPCMMassEtaMCpPb->DrawCopy("same,e1,p,x0"); 

	//********************************** Defintion of the Legend **************************************************   
	Double_t rowsLegendMass2AbsEta[7]       = {0.2,566.95,564.75,562.25, 560.0, 143.5, 142.6};
	
	//****************** first Column **************************************************
	textMass2PCM->Draw();
	TLatex *textMass2PCMEMCALEta = new TLatex(columnsLegendMass2[0],rowsLegendMass2[2],"PCM-EMCal");
	SetStyleTLatex( textMass2PCMEMCALEta, textSizeLeftColumnMass2,4);
	textMass2PCMEMCALEta->SetTextFont(42);
	textMass2PCMEMCALEta->Draw();
	TLatex *textMass2EMCALEMCALEta = new TLatex(columnsLegendMass2[0],rowsLegendMass2[3],"EMCal-EMCal");
	SetStyleTLatex( textMass2EMCALEMCALEta, textSizeLeftColumnMass2,4);
	textMass2EMCALEMCALEta->SetTextFont(42);
	textMass2EMCALEMCALEta->Draw();
	TLatex *textMass2PCMPHOSEta = new TLatex(columnsLegendMass2[0],rowsLegendMass2[4],"PCM-PHOS");
	SetStyleTLatex( textMass2PCMPHOSEta, textSizeLeftColumnMass2,4);
	textMass2PCMPHOSEta->SetTextFont(42);
	textMass2PCMPHOSEta->Draw();

	
	//****************** second Column *************************************************
	textMass2Data2->Draw();
	textMass2MC2->Draw();
	
	markerPCMPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2AbsEta[1]);
	markerPCMEMCALPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2AbsEta[2]);
	markerEMCALEMCALPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2AbsEta[3]);
	markerPCMPHOSPi0Mass2->DrawMarker(columnsLegendMass2Abs[1] ,rowsLegendMass2AbsEta[4]);
	
	markerPCMPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2] ,rowsLegendMass2AbsEta[1]);
	markerPCMEMCALPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2AbsEta[2]);
	markerEMCALEMCALPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2AbsEta[3]);
	markerPCMPHOSPi0Mass2MC->DrawMarker(columnsLegendMass2Abs[2],rowsLegendMass2AbsEta[4]);

	DrawGammaLines(0.25, 30. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcEtaMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Eta_MassDifferentDetectors_pPb.%s",outputDir.Data(), suffix.Data()));
	
	//****************************************************************************************************
	//********************** Plotting Relative Mass accuracy for diff detectors **************************
	//****************************************************************************************************
	
	DrawGammaCanvasSettings( canvasMass, 0.1, 0.02, 0.01, 0.09);   
	canvasMass->cd();
		
	canvasMass->SetLogx();
	
	TH2D *histo2DAllPi0MassRelAcc;
	histo2DAllPi0MassRelAcc = new TH2D("histo2DAllPi0MassRelAcc", "histo2DAllPi0MassRelAcc", 20,0.25,30. ,1000.,-3.5,4.5);
	SetStyleHistoTH2ForGraphs(histo2DAllPi0MassRelAcc, "#it{p}_{T} (GeV/#it{c})","#frac{#it{m}_{#pi^{0}, data}- #it{m}_{#pi^{0}, MC}}{#it{m}_{#pi^{0}, PDG}} (%)", 0.85*textSizeMassWidth,textSizeMassWidth, 0.85*textSizeMassWidth,textSizeMassWidth, 0.9,1, 515, 510);
	histo2DAllPi0MassRelAcc->GetXaxis()->SetLabelOffset(-0.01);
	histo2DAllPi0MassRelAcc->GetYaxis()->SetLabelOffset(0.01);
	histo2DAllPi0MassRelAcc->SetLineWidth(2);
	histo2DAllPi0MassRelAcc->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DAllPi0MassRelAcc->GetYaxis()->SetDecimals();
	histo2DAllPi0MassRelAcc->DrawCopy();

	DrawGammaSetMarker(histoRatioPHOSMassPi0DiffDataMC, markerStylePHOS, markerSizeMass, colorPHOS, colorPHOS);
	histoRatioPHOSMassPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	DrawGammaSetMarker(histoRatioDalitzMassPi0DiffDataMC, markerStyleDalitz, markerSizeMass, colorDalitz, colorDalitz);
	histoRatioDalitzMassPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 
	
	DrawGammaSetMarker( histoRatioEMCALMassPi0DiffDataMC, markerStyleEMCAL, markerSizeMass*1.5, colorEMCAL, colorEMCAL);
	histoRatioEMCALMassPi0DiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoRatioPCMEMCALMassPi0DiffDataMC, markerStylePCMEMCAL, markerSizeMass, colorPCMEMCAL, colorPCMEMCAL);
	histoRatioPCMEMCALMassPi0DiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoRatioPCMPHOSMassPi0DiffDataMC, markerStylePCMPHOS, markerSizeMass, colorPCMPHOS, colorPCMPHOS);
	histoRatioPCMPHOSMassPi0DiffDataMC->Draw("same,e1,p,x0"); 

// 	DrawGammaSetMarker( histoRatioEMCALEMCALMassPi0DiffDataMC, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
// 	histoRatioEMCALEMCALMassPi0DiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker(histoRatioPCMMassPi0DiffDataMC, markerStylePCM, markerSizeMass, colorPCM, colorPCM);
	histoRatioPCMMassPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	//********************************** Defintion of the Legend **************************************************   
	Double_t columnsLegendMass3[4]    = {0.2,0.405,0.435,0.39};
	Double_t columnsLegendMass3Abs[4]    = {4,1.45,4.8,12};
	//    Double_t rowsLegendMass3[3]       = {0.66,0.33,0.0};
	Double_t rowsLegendMass3[7]       = {0.93,0.89,0.85,0.81, 0.77,0.73, 0.69};
	Double_t rowsLegendMass3Abs[7]       = {0.2, 3.7, 3.35, 3.0, 2.65, 2.3, 1.95 };
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnMass3 = textSizeMassWidth*0.85;
	Size_t textSizeTopRowMass3  = textSizeMassWidth*0.85; 
	//****************** Scale factors ******************
	Double_t scaleMarkerMass3      = 1.1;
// 	
// 	//    padMassLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textMass3PCM = new TLatex(columnsLegendMass3[0],rowsLegendMass3[1],"PCM");
	SetStyleTLatex( textMass3PCM, textSizeLeftColumnMass3,4);
	textMass3PCM->SetTextFont(42);
	textMass3PCM->Draw();
	TLatex *textMass3PCMEMCAL = new TLatex(columnsLegendMass3[0],rowsLegendMass3[2],"PCM-EMCal");
	SetStyleTLatex( textMass3PCMEMCAL, textSizeLeftColumnMass3,4);
	textMass3PCMEMCAL->SetTextFont(42);
	textMass3PCMEMCAL->Draw();
	TLatex *textMass3EMCAL = new TLatex(columnsLegendMass3[0],rowsLegendMass3[3],"EMCal");
	SetStyleTLatex( textMass3EMCAL, textSizeLeftColumnMass3,4);
	textMass3EMCAL->SetTextFont(42);
	textMass3EMCAL->Draw();
// 	TLatex *textMass3EMCALEMCAL = new TLatex(columnsLegendMass3[0],rowsLegendMass3[4],"EMCal-EMCal");
// 	SetStyleTLatex( textMass3EMCALEMCAL, textSizeLeftColumnMass3,4);
// 	textMass3EMCALEMCAL->SetTextFont(42);
// 	textMass3EMCALEMCAL->Draw();
	TLatex *textMass3Dalitz = new TLatex(columnsLegendMass3[0],rowsLegendMass3[4],"Dalitz");
	SetStyleTLatex( textMass3Dalitz, textSizeLeftColumnMass3,4);
	textMass3Dalitz->SetTextFont(42);
	textMass3Dalitz->Draw();
	TLatex *textMass3PHOS = new TLatex(columnsLegendMass3[0],rowsLegendMass3[5],"PHOS");
	SetStyleTLatex( textMass3PHOS, textSizeLeftColumnMass3,4);
	textMass3PHOS->SetTextFont(42);
	textMass3PHOS->Draw();
	TLatex *textMass3PCMPHOS = new TLatex(columnsLegendMass3[0],rowsLegendMass3[6],"PCM-PHOS");
	SetStyleTLatex( textMass3PCMPHOS, textSizeLeftColumnMass3,4);
	textMass3PCMPHOS->SetTextFont(42);
	textMass3PCMPHOS->Draw();

	
	//****************** second Column *************************************************
	TLatex *textMass3Data2 = new TLatex(columnsLegendMass3[1],rowsLegendMass3[0] ,"data");
	SetStyleTLatex( textMass3Data2, textSizeTopRowMass3 ,4);
	textMass3Data2->SetTextFont(42);
	textMass3Data2->Draw();
	TMarker* markerPCMPi0Mass3 = CreateMarkerFromHisto(histoPCMMassPi0DatapPb,columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[1] ,scaleMarkerMass3);
	markerPCMPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[1]);
	TMarker* markerPCMEMCALPi0Mass3 = CreateMarkerFromHisto(histoPCMEMCALMassPi0DatapPb,columnsLegendMass3Abs[1],rowsLegendMass3Abs[2],scaleMarkerMass3);
	markerPCMEMCALPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[2]);
	TMarker* markerEMCALPi0Mass3 = CreateMarkerFromHisto(histoEMCALMassPi0DatapPb,columnsLegendMass3Abs[1],rowsLegendMass3Abs[3],scaleMarkerMass3);
	markerEMCALPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[3]);
// 	TMarker* markerEMCALEMCALPi0Mass3 = CreateMarkerFromHisto(histoEMCALEMCALMassPi0DatapPb,columnsLegendMass3Abs[1],rowsLegendMass3Abs[4],scaleMarkerMass3);
// 	markerEMCALEMCALPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[4]);
	TMarker* markerDalitzPi0Mass3 = CreateMarkerFromHisto(histoDalitzMassPi0DatapPb,columnsLegendMass3Abs[1],rowsLegendMass3Abs[4],scaleMarkerMass3);
	markerDalitzPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[4]);
	TMarker* markerPHOSPi0Mass3 = CreateMarkerFromHisto(histoPHOSMassPi0DatapPb,columnsLegendMass3Abs[1],rowsLegendMass3Abs[5],scaleMarkerMass3);
	markerPHOSPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[5]);
	TMarker* markerPCMPHOSPi0Mass3 = CreateMarkerFromHisto(histoPCMPHOSMassPi0DatapPb,columnsLegendMass3Abs[1],rowsLegendMass3Abs[6],scaleMarkerMass3);
	markerPCMPHOSPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[6]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Pi0_MassRelAccDifferentDetectors_pPb.%s",outputDir.Data(), suffix.Data()));


	//****************************************************************************************************
	//***************** Plotting Relative Mass accuracy for diff detectors only new EMCAL ****************
	//****************************************************************************************************
		
	DrawGammaCanvasSettings( canvasMass, 0.1, 0.02, 0.01, 0.09);   
	canvasMass->cd();
	canvasMass->SetLogx();

	histo2DAllPi0MassRelAcc->DrawCopy();
	
	histoRatioEMCALMassPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMEMCALMassPi0DiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoRatioEMCALEMCALMassPi0DiffDataMC, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
	histoRatioEMCALEMCALMassPi0DiffDataMC->Draw("same,e1,p,x0"); 

	histoRatioPCMMassPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textMass3PCM->Draw();
	textMass3PCMEMCAL->Draw();
	textMass3EMCAL->Draw();
	TLatex *textMass3EMCALEMCAL = new TLatex(columnsLegendMass3[0],rowsLegendMass3[4],"EMCal-EMCal");
	SetStyleTLatex( textMass3EMCALEMCAL, textSizeLeftColumnMass3,4);
	textMass3EMCALEMCAL->SetTextFont(42);
	textMass3EMCALEMCAL->Draw();
	
	//****************** second Column *************************************************
	textMass3Data2->Draw();
	markerPCMPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[1]);
	markerPCMEMCALPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[2]);
	markerEMCALPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[3]);
	TMarker* markerEMCALEMCALPi0Mass3 = CreateMarkerFromHisto(histoEMCALEMCALMassPi0DatapPb,columnsLegendMass3Abs[1],rowsLegendMass3Abs[4],scaleMarkerMass3);
	markerEMCALEMCALPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[4]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Pi0_MassRelAccDifferentDetectorsOnlyNewRefEMCAL_pPb.%s",outputDir.Data(), suffix.Data()));

		
	canvasMass->cd();
	canvasMass->SetLogx();

	histo2DAllPi0MassRelAcc->DrawCopy();
	
	histoRatioEMCALMassPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMEMCALMassPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMMassPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textMass3PCM->Draw();
	textMass3PCMEMCAL->Draw();
	textMass3EMCAL->Draw();
	
	//****************** second Column *************************************************
	textMass3Data2->Draw();
	markerPCMPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[1]);
	markerPCMEMCALPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[2]);
	markerEMCALPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[3]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Pi0_MassRelAccDifferentDetectorsOnlyNewRefEMCAL_stripped_pPb.%s",outputDir.Data(), suffix.Data()));
	
	//****************************************************************************************************
	//***************** Plotting Relative Mass accuracy for diff detectors only new PHOS ****************
	//****************************************************************************************************

	DrawGammaCanvasSettings( canvasMass, 0.1, 0.02, 0.01, 0.09);   
	canvasMass->cd();
	canvasMass->SetLogx();
		
	histo2DAllPi0MassRelAcc->DrawCopy();
	
	histoRatioPHOSMassPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMPHOSMassPi0DiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoRatioPHOSPHOSMassPi0DiffDataMC, markerStylePHOS, markerSizeMass*0.7, colorPHOS+1, colorPHOS+1);
	histoRatioPHOSPHOSMassPi0DiffDataMC->Draw("same,e1,p,x0"); 

	histoRatioPCMMassPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textMass3PCM->Draw();
	TLatex *textMass3PCMPHOS2 = new TLatex(columnsLegendMass3[0],rowsLegendMass3[2],"PCM-PHOS");
	SetStyleTLatex( textMass3PCMPHOS2, textSizeLeftColumnMass3,4);
	textMass3PCMPHOS2->SetTextFont(42);
	textMass3PCMPHOS2->Draw();
	TLatex *textMass3PHOS2 = new TLatex(columnsLegendMass3[0],rowsLegendMass3[3],"PHOS");
	SetStyleTLatex( textMass3PHOS2, textSizeLeftColumnMass3,4);
	textMass3PHOS2->SetTextFont(42);
	textMass3PHOS2->Draw();
	TLatex *textMass3PHOSPHOS = new TLatex(columnsLegendMass3[0],rowsLegendMass3[4],"PHOS-PHOS");
	SetStyleTLatex( textMass3PHOSPHOS, textSizeLeftColumnMass3,4);
	textMass3PHOSPHOS->SetTextFont(42);
	textMass3PHOSPHOS->Draw();
	
	//****************** second Column *************************************************
	textMass3Data2->Draw();
	markerPCMPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[1]);
	markerPCMPHOSPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[2]);
	markerPHOSPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[3]);
	TMarker* markerPHOSPHOSPi0Mass3 = CreateMarkerFromHisto(histoPHOSPHOSMassPi0DatapPb,columnsLegendMass3Abs[1],rowsLegendMass3Abs[4],scaleMarkerMass3);
	markerPHOSPHOSPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[4]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Pi0_MassRelAccDifferentDetectorsOnlyNewRefPHOS_pPb.%s",outputDir.Data(), suffix.Data()));

	canvasMass->cd();
	canvasMass->SetLogx();
		
	histo2DAllPi0MassRelAcc->DrawCopy();
	
	histoRatioPHOSMassPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMPHOSMassPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMMassPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textMass3PCM->Draw();
	textMass3PCMPHOS2->Draw();
	textMass3PHOS2->Draw();
	
	//****************** second Column *************************************************
	textMass3Data2->Draw();
	markerPCMPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[1]);
	markerPCMPHOSPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[2]);
	markerPHOSPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3Abs[3]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Pi0_MassRelAccDifferentDetectorsOnlyNewRefPHOS_stripped_pPb.%s",outputDir.Data(), suffix.Data()));
	
	//****************************************************************************************************
	//********************** Plotting Relative Eta Mass accuracy for diff detectors **************************
	//****************************************************************************************************
	
	DrawGammaCanvasSettings( canvasMass, 0.1, 0.02, 0.01, 0.09);   
	canvasMass->cd();
		
	canvasMass->SetLogx();
	
	TH2D *histo2DAllEtaMassRelAcc;
	histo2DAllEtaMassRelAcc = new TH2D("histo2DAllEtaMassRelAcc", "histo2DAllEtaMassRelAcc", 20,0.25,30. ,1000.,-3.5,4.5);
	SetStyleHistoTH2ForGraphs(histo2DAllEtaMassRelAcc, "#it{p}_{T} (GeV/#it{c})","#frac{#it{m}_{#eta, data}- #it{m}_{#eta, MC}}{#it{m}_{#eta, PDG}} (%)", 0.85*textSizeMassWidth,textSizeMassWidth, 0.85*textSizeMassWidth,textSizeMassWidth, 0.9,1, 515, 510);
	histo2DAllEtaMassRelAcc->GetXaxis()->SetLabelOffset(-0.01);
	histo2DAllEtaMassRelAcc->GetYaxis()->SetLabelOffset(0.01);
	histo2DAllEtaMassRelAcc->SetLineWidth(2);
	histo2DAllEtaMassRelAcc->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DAllEtaMassRelAcc->GetYaxis()->SetDecimals();
	histo2DAllEtaMassRelAcc->DrawCopy();

	DrawGammaSetMarker( histoRatioPCMEMCALMassEtaDiffDataMC, markerStylePCMEMCAL, markerSizeMass, colorPCMEMCAL, colorPCMEMCAL);
	histoRatioPCMEMCALMassEtaDiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoRatioPCMPHOSMassEtaDiffDataMC, markerStylePCMPHOS, markerSizeMass, colorPCMPHOS, colorPCMPHOS);
	histoRatioPCMPHOSMassEtaDiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoRatioEMCALEMCALMassEtaDiffDataMC, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
	histoRatioEMCALEMCALMassEtaDiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker(histoRatioPCMMassEtaDiffDataMC, markerStylePCM, markerSizeMass, colorPCM, colorPCM);
	histoRatioPCMMassEtaDiffDataMC->DrawCopy("same,e1,p,x0"); 

	Double_t rowsLegendMass3AbsEta[7]       = {0.2, 3.7, 3.35, 3.0, 2.65, 2.3, 1.95 };

	//****************** first Column **************************************************
	textMass3PCM->Draw();
	textMass3PCMEMCAL->Draw();
	TLatex *textMass3EMCALEMCAL2 = new TLatex(columnsLegendMass3[0],rowsLegendMass3[3],"EMCal-EMCal");
	SetStyleTLatex( textMass3EMCALEMCAL2, textSizeLeftColumnMass3,4);
	textMass3EMCALEMCAL2->SetTextFont(42);
	textMass3EMCALEMCAL2->Draw();
	TLatex *textMass3PCMPHOS3 = new TLatex(columnsLegendMass3[0],rowsLegendMass3[4],"PCM-PHOS");
	SetStyleTLatex( textMass3PCMPHOS3, textSizeLeftColumnMass3,4);
	textMass3PCMPHOS3->SetTextFont(42);
	textMass3PCMPHOS3->Draw();

	//****************** second Column *************************************************
	textMass3Data2->Draw();
	markerPCMPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3AbsEta[1]);
	markerPCMEMCALPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3AbsEta[2]);
	markerEMCALEMCALPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3AbsEta[3]);
	markerPCMPHOSPi0Mass3->DrawMarker(columnsLegendMass3Abs[1] ,rowsLegendMass3AbsEta[4]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcEtaMassWidth->Draw();
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/Eta_MassRelAccDifferentDetectors_pPb.%s",outputDir.Data(), suffix.Data()));





	//****************************************************************************************************
	//********************** Plotting Relative Width accuracy for diff detectors **************************
	//****************************************************************************************************
	
	DrawGammaCanvasSettings( canvasFWHM, 0.1, 0.02, 0.01, 0.09);   
	canvasFWHM->cd();
		
	canvasFWHM->SetLogx();
	
	TH2D *histo2DAllPi0WidthRelAcc;
	histo2DAllPi0WidthRelAcc = new TH2D("histo2DAllPi0WidthRelAcc", "histo2DAllPi0WidthRelAcc", 20,0.25,30. ,1000.,-130,130);
	SetStyleHistoTH2ForGraphs(histo2DAllPi0WidthRelAcc, "#it{p}_{T} (GeV/#it{c})","#frac{#it{#sigma}_{#pi^{0}, data}- #it{#sigma}_{#pi^{0}, MC}}{#it{#sigma}_{#pi^{0}, data}} (%)", 0.85*textSizeMassWidth,textSizeMassWidth, 0.85*textSizeMassWidth,textSizeMassWidth, 0.9,1, 515, 510);
	histo2DAllPi0WidthRelAcc->GetXaxis()->SetLabelOffset(-0.01);
	histo2DAllPi0WidthRelAcc->GetYaxis()->SetLabelOffset(0.01);
	histo2DAllPi0WidthRelAcc->SetLineWidth(2);
	histo2DAllPi0WidthRelAcc->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DAllPi0WidthRelAcc->GetYaxis()->SetDecimals();
	histo2DAllPi0WidthRelAcc->DrawCopy();

	DrawGammaSetMarker(histoRatioPHOSWidthPi0DiffDataMC, markerStylePHOS, markerSizeMass, colorPHOS, colorPHOS);
	histoRatioPHOSWidthPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	DrawGammaSetMarker(histoRatioDalitzWidthPi0DiffDataMC, markerStyleDalitz, markerSizeMass, colorDalitz, colorDalitz);
	histoRatioDalitzWidthPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 
	
	DrawGammaSetMarker( histoRatioEMCALWidthPi0DiffDataMC, markerStyleEMCAL, markerSizeMass*1.5, colorEMCAL, colorEMCAL);
	histoRatioEMCALWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoRatioPCMEMCALWidthPi0DiffDataMC, markerStylePCMEMCAL, markerSizeMass, colorPCMEMCAL, colorPCMEMCAL);
	histoRatioPCMEMCALWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoRatioPCMPHOSWidthPi0DiffDataMC, markerStylePCMPHOS, markerSizeMass, colorPCMPHOS, colorPCMPHOS);
	histoRatioPCMPHOSWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 

// 	DrawGammaSetMarker( histoRatioEMCALEMCALWidthPi0DiffDataMC, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
// 	histoRatioEMCALEMCALWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker(histoRatioPCMWidthPi0DiffDataMC, markerStylePCM, markerSizeMass, colorPCM, colorPCM);
	histoRatioPCMWidthPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	//********************************** Defintion of the Legend **************************************************   
	Double_t columnsLegendWidth3[4]    = {0.2,0.405,0.435,0.39};
	Double_t columnsLegendWidth3Abs[4]    = {4,1.45,4.8,12};
	//    Double_t rowsLegendWidth3[3]       = {0.66,0.33,0.0};
	Double_t rowsLegendWidth3[7]       = {0.93,0.89,0.85,0.81, 0.77,0.73, 0.69};
	Double_t rowsLegendWidth3Abs[7]       = {0.2, 104, 93, 82, 70, 59, 47 };
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnWidth3 = textSizeMassWidth*0.85;
	Size_t textSizeTopRowWidth3  = textSizeMassWidth*0.85; 
	//****************** Scale factors ******************
	Double_t scaleMarkerWidth3      = 1.1;
// 	
// 	//    padWidthLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textWidth3PCM = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[1],"PCM");
	SetStyleTLatex( textWidth3PCM, textSizeLeftColumnWidth3,4);
	textWidth3PCM->SetTextFont(42);
	textWidth3PCM->Draw();
	TLatex *textWidth3PCMEMCAL = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[2],"PCM-EMCal");
	SetStyleTLatex( textWidth3PCMEMCAL, textSizeLeftColumnWidth3,4);
	textWidth3PCMEMCAL->SetTextFont(42);
	textWidth3PCMEMCAL->Draw();
	TLatex *textWidth3EMCAL = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[3],"EMCal");
	SetStyleTLatex( textWidth3EMCAL, textSizeLeftColumnWidth3,4);
	textWidth3EMCAL->SetTextFont(42);
	textWidth3EMCAL->Draw();
// 	TLatex *textWidth3EMCALEMCAL = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[4],"EMCal-EMCal");
// 	SetStyleTLatex( textWidth3EMCALEMCAL, textSizeLeftColumnWidth3,4);
// 	textWidth3EMCALEMCAL->SetTextFont(42);
// 	textWidth3EMCALEMCAL->Draw();
	TLatex *textWidth3Dalitz = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[4],"Dalitz");
	SetStyleTLatex( textWidth3Dalitz, textSizeLeftColumnWidth3,4);
	textWidth3Dalitz->SetTextFont(42);
	textWidth3Dalitz->Draw();
	TLatex *textWidth3PHOS = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[5],"PHOS");
	SetStyleTLatex( textWidth3PHOS, textSizeLeftColumnWidth3,4);
	textWidth3PHOS->SetTextFont(42);
	textWidth3PHOS->Draw();
	TLatex *textWidth3PCMPHOS = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[6],"PCM-PHOS");
	SetStyleTLatex( textWidth3PCMPHOS, textSizeLeftColumnWidth3,4);
	textWidth3PCMPHOS->SetTextFont(42);
	textWidth3PCMPHOS->Draw();

	
	//****************** second Column *************************************************
	TLatex *textWidth3Data2 = new TLatex(columnsLegendWidth3[1],rowsLegendWidth3[0] ,"data");
	SetStyleTLatex( textWidth3Data2, textSizeTopRowWidth3 ,4);
	textWidth3Data2->SetTextFont(42);
	textWidth3Data2->Draw();
	TMarker* markerPCMPi0Width3 = CreateMarkerFromHisto(histoPCMWidthPi0DatapPb,columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[1] ,scaleMarkerWidth3);
	markerPCMPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[1]);
	TMarker* markerPCMEMCALPi0Width3 = CreateMarkerFromHisto(histoPCMEMCALWidthPi0DatapPb,columnsLegendWidth3Abs[1],rowsLegendWidth3Abs[2],scaleMarkerWidth3);
	markerPCMEMCALPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[2]);
	TMarker* markerEMCALPi0Width3 = CreateMarkerFromHisto(histoEMCALWidthPi0DatapPb,columnsLegendWidth3Abs[1],rowsLegendWidth3Abs[3],scaleMarkerWidth3);
	markerEMCALPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[3]);
// 	TMarker* markerEMCALEMCALPi0Width3 = CreateMarkerFromHisto(histoEMCALEMCALWidthPi0DatapPb,columnsLegendWidth3Abs[1],rowsLegendWidth3Abs[4],scaleMarkerWidth3);
// 	markerEMCALEMCALPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[4]);
	TMarker* markerDalitzPi0Width3 = CreateMarkerFromHisto(histoDalitzWidthPi0DatapPb,columnsLegendWidth3Abs[1],rowsLegendWidth3Abs[4],scaleMarkerWidth3);
	markerDalitzPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[4]);
	TMarker* markerPHOSPi0Width3 = CreateMarkerFromHisto(histoPHOSWidthPi0DatapPb,columnsLegendWidth3Abs[1],rowsLegendWidth3Abs[5],scaleMarkerWidth3);
	markerPHOSPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[5]);
	TMarker* markerPCMPHOSPi0Width3 = CreateMarkerFromHisto(histoPCMPHOSWidthPi0DatapPb,columnsLegendWidth3Abs[1],rowsLegendWidth3Abs[6],scaleMarkerWidth3);
	markerPCMPHOSPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[6]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Pi0_WidthRelAccDifferentDetectors_pPb.%s",outputDir.Data(), suffix.Data()));


	//****************************************************************************************************
	//***************** Plotting Relative Width accuracy for diff detectors only new EMCAL ****************
	//****************************************************************************************************
		
	DrawGammaCanvasSettings( canvasFWHM, 0.1, 0.02, 0.01, 0.09);   
	canvasFWHM->cd();
	canvasFWHM->SetLogx();

	histo2DAllPi0WidthRelAcc->DrawCopy();
	
	histoRatioEMCALWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMEMCALWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoRatioEMCALEMCALWidthPi0DiffDataMC, markerStyleEMCAL, markerSizeMass, colorEMCAL+1, colorEMCAL+1);
	histoRatioEMCALEMCALWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 

	histoRatioPCMWidthPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textWidth3PCM->Draw();
	textWidth3PCMEMCAL->Draw();
	textWidth3EMCAL->Draw();
	TLatex *textWidth3EMCALEMCAL = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[4],"EMCal-EMCal");
	SetStyleTLatex( textWidth3EMCALEMCAL, textSizeLeftColumnWidth3,4);
	textWidth3EMCALEMCAL->SetTextFont(42);
	textWidth3EMCALEMCAL->Draw();
	
	//****************** second Column *************************************************
	textWidth3Data2->Draw();
	markerPCMPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[1]);
	markerPCMEMCALPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[2]);
	markerEMCALPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[3]);
	TMarker* markerEMCALEMCALPi0Width3 = CreateMarkerFromHisto(histoEMCALEMCALWidthPi0DatapPb,columnsLegendWidth3Abs[1],rowsLegendWidth3Abs[4],scaleMarkerWidth3);
	markerEMCALEMCALPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[4]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Pi0_WidthRelAccDifferentDetectorsOnlyNewRefEMCAL_pPb.%s",outputDir.Data(), suffix.Data()));

		
	canvasFWHM->cd();
	canvasFWHM->SetLogx();

	histo2DAllPi0WidthRelAcc->DrawCopy();
	
	histoRatioEMCALWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMEMCALWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMWidthPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textWidth3PCM->Draw();
	textWidth3PCMEMCAL->Draw();
	textWidth3EMCAL->Draw();
	
	//****************** second Column *************************************************
	textWidth3Data2->Draw();
	markerPCMPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[1]);
	markerPCMEMCALPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[2]);
	markerEMCALPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[3]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Pi0_WidthRelAccDifferentDetectorsOnlyNewRefEMCAL_stripped_pPb.%s",outputDir.Data(), suffix.Data()));
	
	//****************************************************************************************************
	//***************** Plotting Relative Width accuracy for diff detectors only new PHOS ****************
	//****************************************************************************************************

	DrawGammaCanvasSettings( canvasFWHM, 0.1, 0.02, 0.01, 0.09);   
	canvasFWHM->cd();
	canvasFWHM->SetLogx();
		
	histo2DAllPi0WidthRelAcc->DrawCopy();
	
	histoRatioPHOSWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMPHOSWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 

	DrawGammaSetMarker( histoRatioPHOSPHOSWidthPi0DiffDataMC, markerStylePHOS, markerSizeMass*0.7, colorPHOS+1, colorPHOS+1);
	histoRatioPHOSPHOSWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 

	histoRatioPCMWidthPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textWidth3PCM->Draw();
	TLatex *textWidth3PCMPHOS2 = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[2],"PCM-PHOS");
	SetStyleTLatex( textWidth3PCMPHOS2, textSizeLeftColumnWidth3,4);
	textWidth3PCMPHOS2->SetTextFont(42);
	textWidth3PCMPHOS2->Draw();
	TLatex *textWidth3PHOS2 = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[3],"PHOS");
	SetStyleTLatex( textWidth3PHOS2, textSizeLeftColumnWidth3,4);
	textWidth3PHOS2->SetTextFont(42);
	textWidth3PHOS2->Draw();
	TLatex *textWidth3PHOSPHOS = new TLatex(columnsLegendWidth3[0],rowsLegendWidth3[4],"PHOS-PHOS");
	SetStyleTLatex( textWidth3PHOSPHOS, textSizeLeftColumnWidth3,4);
	textWidth3PHOSPHOS->SetTextFont(42);
	textWidth3PHOSPHOS->Draw();
	
	//****************** second Column *************************************************
	textWidth3Data2->Draw();
	markerPCMPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[1]);
	markerPCMPHOSPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[2]);
	markerPHOSPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[3]);
	TMarker* markerPHOSPHOSPi0Width3 = CreateMarkerFromHisto(histoPHOSPHOSWidthPi0DatapPb,columnsLegendWidth3Abs[1],rowsLegendWidth3Abs[4],scaleMarkerWidth3);
	markerPHOSPHOSPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[4]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Pi0_WidthRelAccDifferentDetectorsOnlyNewRefPHOS_pPb.%s",outputDir.Data(), suffix.Data()));

	canvasFWHM->cd();
	canvasFWHM->SetLogx();
		
	histo2DAllPi0WidthRelAcc->DrawCopy();
	
	histoRatioPHOSWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMPHOSWidthPi0DiffDataMC->Draw("same,e1,p,x0"); 
	histoRatioPCMWidthPi0DiffDataMC->DrawCopy("same,e1,p,x0"); 

	//****************** first Column **************************************************
	textWidth3PCM->Draw();
	textWidth3PCMPHOS2->Draw();
	textWidth3PHOS2->Draw();
	
	//****************** second Column *************************************************
	textWidth3Data2->Draw();
	markerPCMPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[1]);
	markerPCMPHOSPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[2]);
	markerPHOSPi0Width3->DrawMarker(columnsLegendWidth3Abs[1] ,rowsLegendWidth3Abs[3]);
	
	DrawGammaLines(0.25, 30. , 0., 0.,0.1,colorPCM);

	labelIndMeaspPbMassWidth->Draw();
	labelIndMeasProcMassWidth->Draw();
	
	canvasFWHM->Update();
	canvasFWHM->SaveAs(Form("%s/Pi0_WidthRelAccDifferentDetectorsOnlyNewRefPHOS_stripped_pPb.%s",outputDir.Data(), suffix.Data()));

	
	//****************************************************************************************************
	//*************************** Plotting Pi0 corrected spectra for diff detectors **********************
	//****************************************************************************************************

	TCanvas* canvasSpectraDiffDetMinBias = new TCanvas("canvasSpectraDiffDetMinBias","",200,10,800,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasSpectraDiffDetMinBias,  0.16, 0.01, 0.015, 0.07);
	
	canvasSpectraDiffDetMinBias->SetLogy();
	canvasSpectraDiffDetMinBias->SetLogx();
	TH2F * histo2DPCMSpectraPi0AllCent;
	histo2DPCMSpectraPi0AllCent = new TH2F("histo2DPCMSpectraPi0AllCent","histo2DPCMSpectraPi0AllCent",1000,0.23,25.,1000,1e-9,7e1   );
	SetStyleHistoTH2ForGraphs(histo2DPCMSpectraPi0AllCent, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.85*textSizeMassWidth,textSizeMassWidth, 0.85*textSizeMassWidth,textSizeMassWidth, 0.83,1.7);
	histo2DPCMSpectraPi0AllCent->GetXaxis()->SetLabelOffset(-0.01);
	histo2DPCMSpectraPi0AllCent->GetYaxis()->SetLabelOffset(0.005);
	histo2DPCMSpectraPi0AllCent->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DPCMSpectraPi0AllCent->DrawCopy(); 
	
	DrawGammaSetMarker(histoPCMYieldPi0pPb, markerStylePCM,markerSizeInvYield, colorPCM , colorPCM);
	histoPCMYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoDalitzYieldPi0pPb, markerStyleDalitz,markerSizeInvYield, colorDalitz , colorDalitz);
	histoDalitzYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMEMCALYieldPi0pPb, markerStylePCMEMCAL,markerSizeInvYield, colorPCMEMCAL , colorPCMEMCAL);
	histoPCMEMCALYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoEMCALYieldPi0pPb, markerStyleEMCAL,markerSizeInvYield*1.5, colorEMCAL , colorEMCAL);
	histoEMCALYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoEMCALEMCALYieldPi0pPb, markerStyleEMCAL,markerSizeInvYield, colorEMCAL+1 , colorEMCAL+1);
	//	histoEMCALEMCALYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPHOSYieldPi0pPb, markerStylePHOS,markerSizeInvYield, colorPHOS , colorPHOS);
		histoPHOSYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPHOSPHOSYieldPi0pPb, markerStylePHOS,markerSizeInvYield*0.7, colorPHOS+1 , colorPHOS+1);
	//	histoPHOSPHOSYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMPHOSYieldPi0pPb, markerStylePCMPHOS,markerSizeInvYield, colorPCMPHOS , colorPCMPHOS);
	histoPCMPHOSYieldPi0pPb->Draw("p,same,e1");

	DrawGammaSetMarker(histoPCMYieldEtapPb, markerStylePCM,markerSizeInvYield, colorPCM , colorPCM);
	histoPCMYieldEtapPb->Scale(0.1);
	histoPCMYieldEtapPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMEMCALYieldEtapPb, markerStylePCMEMCAL,markerSizeInvYield, colorPCMEMCAL , colorPCMEMCAL);
	histoPCMEMCALYieldEtapPb->Scale(0.1);
	histoPCMEMCALYieldEtapPb->Draw("p,same,e1");
	//	DrawGammaSetMarker(histoPCMPHOSYieldEtapPb, markerStylePCMPHOS,markerSizeInvYield, colorPCMPHOS , colorPCMPHOS);
	//	histoPCMPHOSYieldEtapPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoEMCALEMCALYieldEtapPb, markerStyleEMCAL,markerSizeInvYield, colorEMCAL+1 , colorEMCAL+1);
	histoEMCALEMCALYieldEtapPb->Scale(0.1);
	histoEMCALEMCALYieldEtapPb->Draw("p,same,e1");

	
	
	TLegend* legendSpectraDiffDetMinBias = new TLegend(0.18,0.13,0.5,0.39);
	legendSpectraDiffDetMinBias->SetFillColor(0);
	legendSpectraDiffDetMinBias->SetLineColor(0);
	legendSpectraDiffDetMinBias->SetTextFont(42);
	legendSpectraDiffDetMinBias->SetTextSize(0.85*textSizeMassWidth);
	legendSpectraDiffDetMinBias->AddEntry(histoPCMYieldPi0pPb,Form("PCM"),"pf");
	legendSpectraDiffDetMinBias->AddEntry(histoPCMEMCALYieldPi0pPb,Form("PCM-EMCal"),"pf");
	legendSpectraDiffDetMinBias->AddEntry(histoPCMPHOSYieldPi0pPb,Form("PCM-PHOS"),"pf");
	legendSpectraDiffDetMinBias->AddEntry(histoDalitzYieldPi0pPb,Form("Dalitz"),"pf");
	legendSpectraDiffDetMinBias->AddEntry(histoEMCALYieldPi0pPb,Form("EMCal"),"pf");
	// legendSpectraDiffDetMinBias->AddEntry(histoEMCALEMCALYieldPi0pPb,Form("EMCal-EMCal"),"pf");
	 legendSpectraDiffDetMinBias->AddEntry(histoPHOSYieldPi0pPb,Form("PHOS"),"pf");
	// legendSpectraDiffDetMinBias->AddEntry(histoPHOSPHOSYieldPi0pPb,Form("PHOS-PHOS"),"pf");
	legendSpectraDiffDetMinBias->Draw();

	TLatex *labelIndMeaspPbSpec = new TLatex(0.65,0.9,Form("%s",collisionSystempPb.Data()));
	SetStyleTLatex( labelIndMeaspPbSpec, 0.85*textSizeMassWidth,4);	
	labelIndMeaspPbSpec->Draw();
	TLatex *labelIndMeasProcSpec = new TLatex(0.6,0.7,Form("#pi^{0} #rightarrow #gamma#gamma"));
	SetStyleTLatex( labelIndMeasProcSpec, 0.85*textSizeMassWidth,4);	
	labelIndMeasProcSpec->Draw();
	TLatex *labelIndMeasProcSpecEta2 = new TLatex(0.38,0.5,Form("#eta #rightarrow #gamma#gamma x 0.1"));
	SetStyleTLatex( labelIndMeasProcSpecEta2, 0.85*textSizeMassWidth,4);	
	labelIndMeasProcSpecEta2->Draw();

	canvasSpectraDiffDetMinBias->Update();
	canvasSpectraDiffDetMinBias->Print(Form("%s/Pi0_Spectra_MinBias_DiffDetectors.%s",outputDir.Data(),suffix.Data()));

	histo2DPCMSpectraPi0AllCent->DrawCopy(); 
	
	DrawGammaSetMarker(histoPCMYieldPi0pPb, markerStylePCM,markerSizeInvYield, colorPCM , colorPCM);
	histoPCMYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoDalitzYieldPi0pPb, markerStyleDalitz,markerSizeInvYield, colorDalitz , colorDalitz);
	histoDalitzYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMPHOSYieldPi0pPb, markerStylePCMPHOS,markerSizeInvYield, colorPCMPHOS , colorPCMPHOS);
	histoPCMPHOSYieldPi0pPb->Draw("p,same,e1");

	DrawGammaSetMarker(histoEMCALYieldPi0pPb, markerStyleEMCAL,markerSizeInvYield*1.5, colorEMCAL , colorEMCAL);
	histoEMCALYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPHOSYieldPi0pPb, markerStylePHOS,markerSizeInvYield, colorPHOS , colorPHOS);
	histoPHOSYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMEMCALYieldPi0pPb, markerStylePCMEMCAL,markerSizeInvYield, colorPCMEMCAL , colorPCMEMCAL);
	histoPCMEMCALYieldPi0pPb->Draw("p,same,e1");

		
	TLegend* legendSpectraDiffDetMinBiasStrip = new TLegend(0.18,0.13,0.5,0.39);
	legendSpectraDiffDetMinBiasStrip->SetFillColor(0);
	legendSpectraDiffDetMinBiasStrip->SetLineColor(0);
	legendSpectraDiffDetMinBiasStrip->SetTextFont(42);
	legendSpectraDiffDetMinBiasStrip->SetTextSize(0.85*textSizeMassWidth);
	legendSpectraDiffDetMinBiasStrip->AddEntry(histoPCMYieldPi0pPb,Form("PCM"),"pf");
	legendSpectraDiffDetMinBiasStrip->AddEntry(histoPCMEMCALYieldPi0pPb,Form("PCM-EMCal"),"pf");
	legendSpectraDiffDetMinBiasStrip->AddEntry(histoPCMPHOSYieldPi0pPb,Form("PCM-PHOS"),"pf");
	legendSpectraDiffDetMinBiasStrip->AddEntry(histoDalitzYieldPi0pPb,Form("Dalitz"),"pf");
	legendSpectraDiffDetMinBiasStrip->AddEntry(histoEMCALYieldPi0pPb,Form("EMCal"),"pf");
	legendSpectraDiffDetMinBiasStrip->AddEntry(histoPHOSYieldPi0pPb,Form("PHOS"),"pf");
	legendSpectraDiffDetMinBiasStrip->Draw();

	labelIndMeaspPbSpec->Draw();
	labelIndMeasProcSpec->Draw();

	canvasSpectraDiffDetMinBias->Update();
	canvasSpectraDiffDetMinBias->Print(Form("%s/Pi0_Spectra_MinBias_DiffDetectors_stripped.%s",outputDir.Data(),suffix.Data()));
	
	//****************************************************************************************************
	//*************************** Plotting Eta corrected spectra for diff detectors **********************
	//****************************************************************************************************

	canvasSpectraDiffDetMinBias->cd();
	canvasSpectraDiffDetMinBias->SetLogy();
	canvasSpectraDiffDetMinBias->SetLogx();
	TH2F * histo2DPCMSpectraEtaAllCent;
	histo2DPCMSpectraEtaAllCent = new TH2F("histo2DPCMSpectraEtaAllCent","histo2DPCMSpectraEtaAllCent",1000,0.23,25.,1000,1e-9,7e1   );
	SetStyleHistoTH2ForGraphs(histo2DPCMSpectraEtaAllCent, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}", 0.85*textSizeMassWidth,textSizeMassWidth, 0.85*textSizeMassWidth,textSizeMassWidth, 0.83,1.4);
	histo2DPCMSpectraEtaAllCent->GetXaxis()->SetLabelOffset(-0.01);
	histo2DPCMSpectraEtaAllCent->GetYaxis()->SetLabelOffset(0.005);
	histo2DPCMSpectraEtaAllCent->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DPCMSpectraEtaAllCent->DrawCopy(); 
	
	DrawGammaSetMarker(histoPCMYieldEtapPb, markerStylePCM,markerSizeInvYield, colorPCM , colorPCM);
	histoPCMYieldEtapPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMEMCALYieldEtapPb, markerStylePCMEMCAL,markerSizeInvYield, colorPCMEMCAL , colorPCMEMCAL);
	histoPCMEMCALYieldEtapPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMPHOSYieldEtapPb, markerStylePCMPHOS,markerSizeInvYield, colorPCMPHOS , colorPCMPHOS);
	histoPCMPHOSYieldEtapPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoEMCALEMCALYieldEtapPb, markerStyleEMCAL,markerSizeInvYield, colorEMCAL+1 , colorEMCAL+1);
	histoEMCALEMCALYieldEtapPb->Draw("p,same,e1");
	
	TLegend* legendSpectraDiffDetMinBiasEta = new TLegend(0.18,0.13,0.5,0.29);
	legendSpectraDiffDetMinBiasEta->SetFillColor(0);
	legendSpectraDiffDetMinBiasEta->SetLineColor(0);
	legendSpectraDiffDetMinBiasEta->SetTextFont(42);
	legendSpectraDiffDetMinBiasEta->SetTextSize(0.85*textSizeMassWidth);
	legendSpectraDiffDetMinBiasEta->AddEntry(histoPCMYieldEtapPb,Form("PCM"),"pf");
	legendSpectraDiffDetMinBiasEta->AddEntry(histoPCMEMCALYieldEtapPb,Form("PCM-EMCal"),"pf");
	legendSpectraDiffDetMinBiasEta->AddEntry(histoPCMPHOSYieldEtapPb,Form("PCM-PHOS"),"pf");
	legendSpectraDiffDetMinBiasEta->AddEntry(histoEMCALEMCALYieldEtapPb,Form("EMCal-EMCal"),"pf");
	legendSpectraDiffDetMinBiasEta->Draw();

	labelIndMeaspPbSpec->Draw();
	TLatex *labelIndMeasProcSpecEta = new TLatex(0.81,0.85,Form("#eta #rightarrow #gamma#gamma"));
	SetStyleTLatex( labelIndMeasProcSpecEta, 0.85*textSizeMassWidth,4);	
	labelIndMeasProcSpecEta->Draw();

	canvasSpectraDiffDetMinBias->Update();
	canvasSpectraDiffDetMinBias->Print(Form("%s/Eta_Spectra_MinBias_DiffDetectors.%s",outputDir.Data(),suffix.Data()));


	//****************************************************************************************************
	//*************************** Comparison to charged pion data for diff detectors **********************
	//****************************************************************************************************

	
	fileChargedPionspPb 			= new TFile("ExternalInputpPb/ChargedPionSpectrapPb_4_Apr_2014.root");
	histoChargedPionSyspPb			= (TH1D*)fileChargedPionspPb->Get("histoChargedPionSpecFullPtSyspPb"); 
	histoChargedPionStatpPb			= (TH1D*)fileChargedPionspPb->Get("histoChargedPionSpecFullPtStatpPb");

	TGraphAsymmErrors* graphPCMYieldPi0SysErrpPbCopy 			= (TGraphAsymmErrors*) graphPCMYieldPi0SysErrpPb->Clone("graphPCMYieldPi0SysErrpPbCopy");
	TGraphAsymmErrors* graphPCMEMCALYieldPi0SysErrpPbCopy 		= (TGraphAsymmErrors*) graphPCMEMCALYieldPi0SysErrpPb->Clone("graphPCMEMCALYieldPi0SysErrpPbCopy");
// 	for (Int_t i = 0; i<2; i++) graphPCMEMCALYieldPi0SysErrpPbCopy->RemovePoint(0);
// 	for (Int_t i = 0; i<2; i++) graphPCMEMCALYieldPi0SysErrpPbCopy->RemovePoint(0);
	
	TGraphAsymmErrors* graphEMCALEMCALYieldPi0SysErrpPbCopy 	= (TGraphAsymmErrors*) graphEMCALEMCALYieldPi0SysErrpPb->Clone("graphEMCALEMCALYieldPi0SysErrpPbCopy");
	for (Int_t i = 0; i<4; i++) graphEMCALEMCALYieldPi0SysErrpPbCopy->RemovePoint(0);
	TGraphAsymmErrors* graphPCMPHOSYieldPi0SysErrpPbCopy 		= (TGraphAsymmErrors*) graphPCMPHOSYieldPi0SysErrpPb->Clone("graphPCMPHOSYieldPi0SysErrpPbCopy");
	TGraphAsymmErrors* graphPHOSPHOSYieldPi0SysErrpPbCopy 		= (TGraphAsymmErrors*) graphPHOSPHOSYieldPi0SysErrpPb->Clone("graphPHOSPHOSYieldPi0SysErrpPbCopy");
	for (Int_t i = 0; i<7; i++) graphPHOSPHOSYieldPi0SysErrpPbCopy->RemovePoint(0);

	TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPbCopy 		= (TGraphAsymmErrors*) graphDalitzYieldPi0SysErrpPb->Clone("graphDalitzYieldPi0SysErrpPbCopy");   
	TGraphAsymmErrors* graphDalitzYieldPi0StatErrpPbCopy 		= (TGraphAsymmErrors*) graphDalitzYieldPi0pPb->Clone("graphDalitzYieldPi0StatErrpPbCopy");   
	cout << "*************************************************************************"<< endl;  
	cout << "******************************  pPb MinBias *****************************"<< endl;
	cout << "*************************************************************************"<< endl;
	// TGraphErrors* bla1 				= NULL;
	// TGraphErrors* bla2 				= NULL;
	// TGraphErrors* bla3 				= NULL;
	// TGraphErrors* bla4 				= NULL;

	cout << "*************************************************************************"<< endl;  
	cout << "******************************  pPb MinBias PCM *************************"<< endl;
	cout << "*************************************************************************"<< endl;
	TGraphErrors* graphRatioChargedPionsPCMpPb 	= CalculateRatioBetweenSpectraWithDifferentBinning(			histoPCMYieldPi0pPb, graphPCMYieldPi0SysErrpPbCopy, 
																											histoChargedPionStatpPb, histoChargedPionSyspPb,  
																											kTRUE,  kTRUE, 
																											&bla1, &bla2, &bla3, &bla4)  ;
	graphRatioChargedPionsPCMpPb->Print();
	
	cout << "*************************************************************************"<< endl;  
	cout << "******************************  pPb MinBias Dalitz **********************"<< endl;
	cout << "*************************************************************************"<< endl;
	TGraphErrors* graphRatioChargedPionsDalitzpPb = CalculateRatioBetweenSpectraWithDifferentBinning(		graphDalitzYieldPi0StatErrpPbCopy, graphDalitzYieldPi0SysErrpPbCopy, 
																											histoChargedPionStatpPb, histoChargedPionSyspPb,  
																											kTRUE,  kTRUE,
																											&bla1, &bla2, &bla3, &bla4)  ;
	graphRatioChargedPionsDalitzpPb->Print();

	cout << "*************************************************************************"<< endl;  
	cout << "******************************  pPb MinBias PHOS ************************"<< endl;
	cout << "*************************************************************************"<< endl;
	TGraphErrors* graphRatioChargedPionsPHOSpPb = CalculateRatioBetweenSpectraWithDifferentBinning(			histoPHOSYieldPi0pPb, histoPHOSYieldPi0pPbSys, 
																											histoChargedPionStatpPb, histoChargedPionSyspPb,  
																											kTRUE,  kTRUE,
																											&bla1, &bla2, &bla3, &bla4)  ;
	graphRatioChargedPionsPHOSpPb->Print();

	cout << "*************************************************************************"<< endl;  
	cout << "******************************  pPb MinBias EMCAL ***********************"<< endl;
	cout << "*************************************************************************"<< endl;	
	cout << "EMCal" << endl;
	TGraphErrors* graphRatioChargedPionsEMCALpPb = CalculateRatioBetweenSpectraWithDifferentBinning(		histoEMCALYieldPi0pPb, histoEMCALYieldPi0pPbSys, 
																											histoChargedPionStatpPb, histoChargedPionSyspPb,  
																											kTRUE,  kTRUE, 
																											&bla1, &bla2, &bla3, &bla4)  ;
	graphRatioChargedPionsEMCALpPb->Print();
	
	cout << "*************************************************************************"<< endl;  
	cout << "******************************  pPb MinBias PCM - EMCAL *****************"<< endl;
	cout << "*************************************************************************"<< endl;
	TGraphErrors* graphRatioChargedPionsPCMEMCALpPb = CalculateRatioBetweenSpectraWithDifferentBinning(		histoPCMEMCALYieldPi0pPb, graphPCMEMCALYieldPi0SysErrpPbCopy, 
																											histoChargedPionStatpPb, histoChargedPionSyspPb,  
																											kTRUE,  kTRUE, 
																											&bla1, &bla2, &bla3, &bla4)  ;
	for (Int_t i = 0; i<3; i++) graphRatioChargedPionsPCMEMCALpPb->RemovePoint(0);
	graphRatioChargedPionsPCMEMCALpPb->Print();

	cout << "*************************************************************************"<< endl;  
	cout << "******************************  pPb MinBias ECMAL - EMCAL ***************"<< endl;
	cout << "*************************************************************************"<< endl;	
	TGraphErrors* graphRatioChargedPionsEMCALEMCALpPb = CalculateRatioBetweenSpectraWithDifferentBinning(	histoEMCALEMCALYieldPi0pPb, graphEMCALEMCALYieldPi0SysErrpPbCopy, 
																											histoChargedPionStatpPb, histoChargedPionSyspPb,  
																											kTRUE,  kTRUE, 
																											&bla1, &bla2, &bla3, &bla4)  ;
	for (Int_t i = 0; i<2; i++) graphRatioChargedPionsEMCALEMCALpPb->RemovePoint(0);
	graphRatioChargedPionsEMCALEMCALpPb->Print();

	cout << "*************************************************************************"<< endl;  
	cout << "******************************  pPb MinBias PCM - PHOS ******************"<< endl;
	cout << "*************************************************************************"<< endl;
	TGraphErrors* graphRatioChargedPionsPCMPHOSpPb = CalculateRatioBetweenSpectraWithDifferentBinning(		histoPCMPHOSYieldPi0pPb, graphPCMPHOSYieldPi0SysErrpPbCopy, 
																											histoChargedPionStatpPb, histoChargedPionSyspPb,  
																											kTRUE,  kTRUE, 
																											&bla1, &bla2, &bla3, &bla4)  ;
// 	for (Int_t i = 0; i<3; i++) graphRatioChargedPionsPCMPHOSpPb->RemovePoint(0);
	graphRatioChargedPionsPCMPHOSpPb->Print();
	

	cout << "*************************************************************************"<< endl;  
	cout << "******************************  pPb MinBias PCHOS - PHOS ****************"<< endl;
	cout << "*************************************************************************"<< endl;	
	TGraphErrors* graphRatioChargedPionsPHOSPHOSpPb = CalculateRatioBetweenSpectraWithDifferentBinning(		histoPHOSPHOSYieldPi0pPb, graphPHOSPHOSYieldPi0SysErrpPbCopy, 
																											histoChargedPionStatpPb, histoChargedPionSyspPb,  
																											kTRUE,  kTRUE, 
																											&bla1, &bla2, &bla3, &bla4)  ;
// 	for (Int_t i = 0; i<2; i++) graphRatioChargedPionsPHOSPHOSpPb->RemovePoint(0);
	graphRatioChargedPionsPHOSPHOSpPb->Print();
	
	TCanvas* canvasCompYieldpPbInd = new TCanvas("canvasCompYieldpPbInd","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvasCompYieldpPbInd,  0.12, 0.02, 0.02, 0.12);
	
	canvasCompYieldpPbInd->SetLogx();
	TH2F * histo2DCompCombinedRatio2;
	histo2DCompCombinedRatio2 = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.23,30.,1000,0.2,4.   );
	SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
	histo2DCompCombinedRatio2->GetXaxis()->SetLabelOffset(-0.02);
	histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.23,30.);
	histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.1,2.1);
	histo2DCompCombinedRatio2->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPHOSpPb,markerStylePHOS,markerSizeComparison,  colorPHOS, colorPHOS);
		graphRatioChargedPionsPHOSpPb->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsEMCALpPb,markerStyleEMCAL,markerSizeComparison*1.2,  colorEMCAL, colorEMCAL);
		graphRatioChargedPionsEMCALpPb->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPCMpPb,markerStylePCM,markerSizeComparison, colorPCM, colorPCM);
		graphRatioChargedPionsPCMpPb->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsDalitzpPb,markerStyleDalitz,markerSizeComparison*1.2, colorDalitz, colorDalitz);
		graphRatioChargedPionsDalitzpPb->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPCMEMCALpPb,markerStylePCMEMCAL,markerSizeComparison,  colorPCMEMCAL, colorPCMEMCAL);
		graphRatioChargedPionsPCMEMCALpPb->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPCMPHOSpPb,markerStylePCMPHOS,markerSizeComparison,  colorPCMPHOS, colorPCMPHOS);
		graphRatioChargedPionsPCMPHOSpPb->Draw("E1psame");

		TLatex *labelRatioPi0pPb = new TLatex(0.16,0.9,collisionSystempPb.Data());
		SetStyleTLatex( labelRatioPi0pPb, 0.06,4);
		labelRatioPi0pPb->Draw();

		TLegend* legendPi0CompIndChargedPionspPb = new TLegend(0.15,0.16,0.95,0.25);
		legendPi0CompIndChargedPionspPb->SetFillColor(0);
		legendPi0CompIndChargedPionspPb->SetFillStyle(0);
		legendPi0CompIndChargedPionspPb->SetLineColor(0);
		legendPi0CompIndChargedPionspPb->SetNColumns(3);
		legendPi0CompIndChargedPionspPb->SetTextSize(0.038);
		legendPi0CompIndChargedPionspPb->SetMargin(0.14);
		legendPi0CompIndChargedPionspPb->AddEntry(graphRatioChargedPionsPCMpPb,Form("#pi^{0}/#pi^{#pm} PCM"),"p");  // -0.8 < y_{lab} < 0.8
		legendPi0CompIndChargedPionspPb->AddEntry(graphRatioChargedPionsDalitzpPb,Form("#pi^{0}/#pi^{#pm} Dalitz"),"p");
		legendPi0CompIndChargedPionspPb->AddEntry(graphRatioChargedPionsPHOSpPb,Form("#pi^{0}/#pi^{#pm} PHOS"),"p");
		legendPi0CompIndChargedPionspPb->AddEntry(graphRatioChargedPionsEMCALpPb,Form("#pi^{0}/#pi^{#pm} EMCal"),"p");  // -0.7 < y_{lab} < 0.7
		legendPi0CompIndChargedPionspPb->AddEntry(graphRatioChargedPionsPCMEMCALpPb,Form("#pi^{0}/#pi^{#pm} PCM-EMCal"),"p");  // -0.7 < y_{lab} < 0.7
		legendPi0CompIndChargedPionspPb->AddEntry(graphRatioChargedPionsPCMPHOSpPb,Form("#pi^{0}/#pi^{#pm} PCM-PHOS"),"p");  // -0.7 < y_{lab} < 0.7

		legendPi0CompIndChargedPionspPb->Draw();
		DrawGammaLines(0., 30. , 1, 1 ,1,kGray);
	
  
	canvasCompYieldpPbInd->Update();
	canvasCompYieldpPbInd->Print(Form("%s/ComparisonChargedToNeutralPions_pPb.%s",outputDir.Data(),suffix.Data()));

	canvasCompYieldpPbInd->cd();
	histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.55,1.55);
	histo2DCompCombinedRatio2->DrawCopy();

// 		graphRatioChargedPionsPHOSpPb->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsEMCALpPb,markerStyleEMCALMC,markerSizeComparison,  colorEMCAL, colorEMCAL);
		graphRatioChargedPionsEMCALpPb->Draw("E1psame");
		graphRatioChargedPionsPCMpPb->DrawClone("E1psame");
// 		graphRatioChargedPionsDalitzpPb->Draw("E1psame");
		
		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsEMCALEMCALpPb,markerStyleEMCAL,markerSizeComparison,  colorEMCAL+1, colorEMCAL+1);
		graphRatioChargedPionsEMCALEMCALpPb->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPCMEMCALpPb,markerStylePCMEMCAL,markerSizeComparison,  colorPCMEMCAL, colorPCMEMCAL);
		graphRatioChargedPionsPCMEMCALpPb->Draw("E1psame");

		labelRatioPi0pPb->Draw();

		TLegend* legendPi0CompIndChargedPionspPb2 = new TLegend(0.13,0.15,0.95,0.25);
		legendPi0CompIndChargedPionspPb2->SetFillColor(0);
		legendPi0CompIndChargedPionspPb2->SetFillStyle(0);
		legendPi0CompIndChargedPionspPb2->SetLineColor(0);
		legendPi0CompIndChargedPionspPb2->SetNColumns(2);
		legendPi0CompIndChargedPionspPb2->SetTextSize(0.038);
		legendPi0CompIndChargedPionspPb2->SetMargin(0.14);
		legendPi0CompIndChargedPionspPb2->AddEntry(graphRatioChargedPionsPCMpPb,Form("#pi^{0}/#pi^{#pm} PCM - Annika/Fredi"),"p");  // -0.8 < y_{lab} < 0.8
		legendPi0CompIndChargedPionspPb2->AddEntry(graphRatioChargedPionsPCMEMCALpPb,Form("#pi^{0}/#pi^{#pm} PCM-EMCal - Fredi"),"p"); // -0.8 < y_{lab} < 0.8 
// 		legendPi0CompIndChargedPionspPb2->AddEntry(graphRatioChargedPionsDalitzpPb,Form("#pi^{0}/#pi^{#pm} Dalitz - Pedro"),"p");
// 		legendPi0CompIndChargedPionspPb2->AddEntry(graphRatioChargedPionsPHOSpPb,Form("#pi^{0}/#pi^{#pm} PHOS - Boris"),"p");
		legendPi0CompIndChargedPionspPb2->AddEntry(graphRatioChargedPionsEMCALpPb,Form("#pi^{0}/#pi^{#pm} EMCal - Jason"),"p");  // -0.7 < y_{lab} < 0.7
		legendPi0CompIndChargedPionspPb2->AddEntry(graphRatioChargedPionsEMCALEMCALpPb,Form("#pi^{0}/#pi^{#pm} EMCal-EMCal - Fredi"),"p");  // -0.7 < y_{lab} < 0.7
		
		legendPi0CompIndChargedPionspPb2->Draw();
		DrawGammaLines(0., 30. , 1, 1 ,1,kGray);
	
  
	canvasCompYieldpPbInd->Update();
	canvasCompYieldpPbInd->Print(Form("%s/ComparisonChargedToNeutralPionsOnlyNewRefEMCAL_pPb.%s",outputDir.Data(),suffix.Data()));

	histo2DCompCombinedRatio2->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPHOSpPb,markerStylePHOSMC,markerSizeComparison,  colorPHOS, colorPHOS);
		graphRatioChargedPionsPHOSpPb->Draw("E1psame");
// 		graphRatioChargedPionsEMCALpPb->Draw("E1psame");
		graphRatioChargedPionsPCMpPb->DrawClone("E1psame");
// 		graphRatioChargedPionsDalitzpPb->Draw("E1psame");
		
		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPHOSPHOSpPb,markerStylePHOS,markerSizeComparison,  colorPHOS+1, colorPHOS+1);
		graphRatioChargedPionsPHOSPHOSpPb->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPCMPHOSpPb,markerStylePCMPHOS,markerSizeComparison,  colorPCMPHOS, colorPCMPHOS);
		graphRatioChargedPionsPCMPHOSpPb->Draw("E1psame");

		labelRatioPi0pPb->Draw();

		TLegend* legendPi0CompIndChargedPionspPb3 = new TLegend(0.13,0.15,0.95,0.25);
		legendPi0CompIndChargedPionspPb3->SetFillColor(0);
		legendPi0CompIndChargedPionspPb3->SetFillStyle(0);
		legendPi0CompIndChargedPionspPb3->SetLineColor(0);
		legendPi0CompIndChargedPionspPb3->SetNColumns(2);
		legendPi0CompIndChargedPionspPb3->SetTextSize(0.038);
		legendPi0CompIndChargedPionspPb3->SetMargin(0.14);
		legendPi0CompIndChargedPionspPb3->AddEntry(graphRatioChargedPionsPCMpPb,Form("#pi^{0}/#pi^{#pm} PCM - Annika/Fredi"),"p");  // -0.8 < y_{lab} < 0.8
		legendPi0CompIndChargedPionspPb3->AddEntry(graphRatioChargedPionsPCMPHOSpPb,Form("#pi^{0}/#pi^{#pm} PCM-PHOS - Fredi"),"p"); // -0.8 < y_{lab} < 0.8 
// 		legendPi0CompIndChargedPionspPb3->AddEntry(graphRatioChargedPionsDalitzpPb,Form("#pi^{0}/#pi^{#pm} Dalitz - Pedro"),"p");
		legendPi0CompIndChargedPionspPb3->AddEntry(graphRatioChargedPionsPHOSpPb,Form("#pi^{0}/#pi^{#pm} PHOS - Boris"),"p");
// 		legendPi0CompIndChargedPionspPb3->AddEntry(graphRatioChargedPionsPHOSpPb,Form("#pi^{0}/#pi^{#pm} EMCal - Jason"),"p");  // -0.7 < y_{lab} < 0.7
		legendPi0CompIndChargedPionspPb3->AddEntry(graphRatioChargedPionsPHOSPHOSpPb,Form("#pi^{0}/#pi^{#pm} PHOS-PHOS - Fredi"),"p");  // -0.7 < y_{lab} < 0.7
		
		legendPi0CompIndChargedPionspPb3->Draw();
		DrawGammaLines(0., 30. , 1, 1 ,1,kGray);
	
  
	canvasCompYieldpPbInd->Update();
	canvasCompYieldpPbInd->Print(Form("%s/ComparisonChargedToNeutralPionsOnlyNewRefPHOS_pPb.%s",outputDir.Data(),suffix.Data()));

	canvasCompYieldpPbInd->cd();
	histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.55,1.55);
	histo2DCompCombinedRatio2->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsEMCALpPb,markerStyleEMCAL,markerSizeComparison,  colorEMCAL, colorEMCAL);
		graphRatioChargedPionsEMCALpPb->Draw("E1psame");
		graphRatioChargedPionsPCMpPb->DrawClone("E1psame");
		
		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPCMEMCALpPb,markerStylePCMEMCAL,markerSizeComparison,  colorPCMEMCAL, colorPCMEMCAL);
		graphRatioChargedPionsPCMEMCALpPb->Draw("E1psame");

		labelRatioPi0pPb->Draw();

		TLegend* legendPi0CompIndChargedPionspPb4 = new TLegend(0.13,0.15,0.95,0.25);
		legendPi0CompIndChargedPionspPb4->SetFillColor(0);
		legendPi0CompIndChargedPionspPb4->SetFillStyle(0);
		legendPi0CompIndChargedPionspPb4->SetLineColor(0);
		legendPi0CompIndChargedPionspPb4->SetNColumns(2);
		legendPi0CompIndChargedPionspPb4->SetTextSize(0.038);
		legendPi0CompIndChargedPionspPb4->SetMargin(0.14);
		legendPi0CompIndChargedPionspPb4->AddEntry(graphRatioChargedPionsPCMpPb,Form("#pi^{0}/#pi^{#pm} PCM"),"p");  // -0.8 < y_{lab} < 0.8
		legendPi0CompIndChargedPionspPb4->AddEntry(graphRatioChargedPionsPCMEMCALpPb,Form("#pi^{0}/#pi^{#pm} PCM-EMCal"),"p"); // -0.8 < y_{lab} < 0.8 
// 		legendPi0CompIndChargedPionspPb4->AddEntry(graphRatioChargedPionsDalitzpPb,Form("#pi^{0}/#pi^{#pm} Dalitz - Pedro"),"p");
// 		legendPi0CompIndChargedPionspPb4->AddEntry(graphRatioChargedPionsPHOSpPb,Form("#pi^{0}/#pi^{#pm} PHOS - Boris"),"p");
		legendPi0CompIndChargedPionspPb4->AddEntry(graphRatioChargedPionsEMCALpPb,Form("#pi^{0}/#pi^{#pm} EMCal"),"p");  // -0.7 < y_{lab} < 0.7
		
		
		legendPi0CompIndChargedPionspPb4->Draw();
		DrawGammaLines(0., 30. , 1, 1 ,1,kGray);
	
  
	canvasCompYieldpPbInd->Update();
	canvasCompYieldpPbInd->Print(Form("%s/ComparisonChargedToNeutralPionsOnlyNewRefEMCAL_stripped_pPb.%s",outputDir.Data(),suffix.Data()));

	histo2DCompCombinedRatio2->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPHOSpPb,markerStylePHOS,markerSizeComparison,  colorPHOS, colorPHOS);
		graphRatioChargedPionsPHOSpPb->Draw("E1psame");
		graphRatioChargedPionsPCMpPb->DrawClone("E1psame");
		
		DrawGammaSetMarkerTGraphErr(graphRatioChargedPionsPCMPHOSpPb,markerStylePCMPHOS,markerSizeComparison,  colorPCMPHOS, colorPCMPHOS);
		graphRatioChargedPionsPCMPHOSpPb->Draw("E1psame");

		labelRatioPi0pPb->Draw();

		TLegend* legendPi0CompIndChargedPionspPb5 = new TLegend(0.13,0.15,0.95,0.25);
		legendPi0CompIndChargedPionspPb5->SetFillColor(0);
		legendPi0CompIndChargedPionspPb5->SetFillStyle(0);
		legendPi0CompIndChargedPionspPb5->SetLineColor(0);
		legendPi0CompIndChargedPionspPb5->SetNColumns(2);
		legendPi0CompIndChargedPionspPb5->SetTextSize(0.038);
		legendPi0CompIndChargedPionspPb5->SetMargin(0.14);
		legendPi0CompIndChargedPionspPb5->AddEntry(graphRatioChargedPionsPCMpPb,Form("#pi^{0}/#pi^{#pm} PCM"),"p");  // -0.8 < y_{lab} < 0.8
		legendPi0CompIndChargedPionspPb5->AddEntry(graphRatioChargedPionsPCMPHOSpPb,Form("#pi^{0}/#pi^{#pm} PCM-PHOS"),"p"); // -0.8 < y_{lab} < 0.8 
		legendPi0CompIndChargedPionspPb5->AddEntry(graphRatioChargedPionsPHOSpPb,Form("#pi^{0}/#pi^{#pm} PHOS"),"p");
		
		legendPi0CompIndChargedPionspPb5->Draw();
		DrawGammaLines(0., 30. , 1, 1 ,1,kGray);
	
  
	canvasCompYieldpPbInd->Update();
	canvasCompYieldpPbInd->Print(Form("%s/ComparisonChargedToNeutralPionsOnlyNewRefPHOS_stripped_pPb.%s",outputDir.Data(),suffix.Data()));

	
	//*********************** Comparison of Measurements in ALICE ***************************************************
	TCanvas* canvasRatioEtaToPi0ALICE = new TCanvas("canvasRatioEtaToPi0ALICE","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaToPi0ALICE, 0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioEtaToPi0ALICE = new TH2D("histo2DRatioEtaToPi0ALICE", "histo2DRatioEtaToPi0ALICE", 20,0.,15.,1000.,-0.4,1.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaToPi0ALICE, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
	histo2DRatioEtaToPi0ALICE->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DRatioEtaToPi0ALICE->Draw();
	
	DrawGammaSetMarker(histoPCMEtaPi0RatiopPb, markerStylePCM, markerSizeComparison*1.5, colorPCM, colorPCM);
	DrawGammaSetMarker(histoEMCALEMCALEtaPi0RatiopPb, markerStyleEMCAL, markerSizeComparison*1.8, colorEMCAL, colorEMCAL);
	DrawGammaSetMarker(histoPCMEMCALEtaPi0RatiopPb,markerStylePCMEMCAL, markerSizeComparison*1.5, colorPCMEMCAL, colorPCMEMCAL);
// 	DrawGammaSetMarker(histoPCMPHOSEtaPi0RatiopPb,markerStylePCMPHOS, markerSizeComparison, colorPCMPHOS, colorPCMPHOS);
	histoPCMEtaPi0RatiopPb->Draw("same,e1p");
	histoPCMEMCALEtaPi0RatiopPb->Draw("same,e1p");
// 	histoPCMPHOSEtaPi0RatiopPb->Draw("same,e1p");
	histoEMCALEMCALEtaPi0RatiopPb->Draw("same,e1p");

	labelRatioPi0pPb->Draw();

	
	TLegend* legendEtaToPi0Ratio = new TLegend(0.3,0.16,0.977,0.34);
	legendEtaToPi0Ratio->SetTextSize(0.04);			
	legendEtaToPi0Ratio->SetFillColor(0);
	legendEtaToPi0Ratio->SetFillStyle(0);
	legendEtaToPi0Ratio->SetBorderSize(0);
	legendEtaToPi0Ratio->SetMargin(0.1);
	legendEtaToPi0Ratio->AddEntry(histoPCMEtaPi0RatiopPb,Form("#pi^{0}/#pi^{#pm} PCM"),"p");  // -0.8 < y_{lab} < 0.8
	legendEtaToPi0Ratio->AddEntry(histoEMCALEMCALEtaPi0RatiopPb,Form("#pi^{0}/#pi^{#pm} EMCal-EMCal"),"p"); // -0.8 < y_{lab} < 0.8 
	legendEtaToPi0Ratio->AddEntry(histoPCMEMCALEtaPi0RatiopPb,Form("#pi^{0}/#pi^{#pm} PCM-EMCAL"),"p");

	legendEtaToPi0Ratio->Draw();
	
	canvasRatioEtaToPi0ALICE->Update();
	canvasRatioEtaToPi0ALICE->SaveAs(Form("%s/EtaToPi0Ratio_DiffSystems_pPb.%s",outputDir.Data(), suffix.Data()));
	
}
	
