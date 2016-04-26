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
#include "CombineMesonMeasurementspPb.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

void CombineMesonMeasurementspPb(const char *fileNameConversionpPb = "", TString suffix = "eps", TString bWCorrection="X"){	
	
	date = ReturnDateString();
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();

	//___________________________________ Declaration of files _____________________________________________
	collisionSystem0020 = "0-20% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";		
	collisionSystem2040 = "20-40% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";      
	collisionSystem4060 = "40-60% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";    
	collisionSystem60100 = "60-100% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
	collisionSystem6080 = "60-80% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
	//	collisionSystem80100 = "80-100% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
	collisionSystempPb = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   
	TString outputDir = Form("%s/CombineMesonMeasurementspPb%s",suffix.Data(),bWCorrection.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	nameFinalResDat = Form("%s/CombinedResultspPb%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
	fileNameCaloEMCALpPb = "ExternalInputpPb/EMCal_LHC13b_AODpass2_300v2.root";
	//   fileNameCaloPHOSpPb = "ExternalInputpPb/PHOS_pPb_08042013.root";
	fileNameCaloPHOSpPb = "ExternalInputpPb/data_PHOSResultsFullCorrection_pPb_20130826.root";   
	TString datePHOS = "2013-05-15";
	TString datePCM = "2013-06-18";


	TString nameHistoPCM = "CorrectedYieldPi0";
	TString nameHistoPCMEta = "CorrectedYieldEta";
	TString nameGraphPCM = "Pi0SystError";
	TString nameGraphPCMEta = "EtaSystError";
	
	minPtForFits = 0.4;
	
	//declaration for printing logo 	
	prefix2 = "data";
	pictDrawingOptions[1] = kFALSE;
	
	mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	cout << mesonMassExpectPi0 << endl;
	
	//****************************************************************************************************
	//************************** Read data for PCM *******************************************************
	//****************************************************************************************************
	cout << "PCM" << endl;
	cout << "0-100%" << endl;
	filePCMpPb = 					new TFile(fileNameConversionpPb);
	directoryPCMPi0pPb = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
	histoPCMNumberOfEventspPb= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV0-100%");
	histoPCMYieldPi0pPb = 				(TH1D*)directoryPCMPi0pPb->Get(nameHistoPCM.Data());
	histoPCMRAWYieldPi0pPb =            (TH1D*)directoryPCMPi0pPb->Get("RawYieldPi0");
	graphPCMYieldPi0SysErrpPb= 	(TGraphAsymmErrors*)directoryPCMPi0pPb->Get(nameGraphPCM.Data());	
	graphPCMYieldPi0SysErrRAApPb= 	(TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystErrorA");	
	histoPCMMassPi0DatapPb = 			(TH1D*)directoryPCMPi0pPb->Get("MassPi0");
	histoPCMMassPi0MCpPb = 			(TH1D*)directoryPCMPi0pPb->Get("TrueMassPi0");
	histoPCMWidthPi0DatapPb = 		(TH1D*)directoryPCMPi0pPb->Get("FWHMPi0MeV");
	histoPCMWidthPi0MCpPb = 			(TH1D*)directoryPCMPi0pPb->Get("TrueFWHMPi0MeV");
	histoPCMMassPi0DatapPb->Scale(1000.);
	histoPCMMassPi0MCpPb->Scale(1000.);
	nPCMEventpPb = 							histoPCMNumberOfEventspPb->GetBinContent(1);
   
	cout << "0-20%" << endl;
	directoryPCMPi0pPb0020 = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_0-20%"); 
	histoPCMNumberOfEvents0020= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV0-20%");
	histoPCMYieldPi0pPb0020 = 				(TH1D*)directoryPCMPi0pPb0020->Get(nameHistoPCM.Data());
	graphPCMYieldPi0SysErrpPb0020= 	(TGraphAsymmErrors*)directoryPCMPi0pPb0020->Get(nameGraphPCM.Data());	
	graphPCMYieldPi0SysErrRAApPb0020= 	(TGraphAsymmErrors*)directoryPCMPi0pPb0020->Get("Pi0SystErrorA");	
	nPCMEventpPb0020 = 							histoPCMNumberOfEvents0020->GetBinContent(1);
	histoPCMMassPi0DatapPb0020 =        (TH1D*)directoryPCMPi0pPb0020->Get("MassPi0");
	histoPCMMassPi0MCpPb0020 =          (TH1D*)directoryPCMPi0pPb0020->Get("TrueMassPi0");
	histoPCMWidthPi0DatapPb0020 =       (TH1D*)directoryPCMPi0pPb0020->Get("FWHMPi0MeV");
	histoPCMWidthPi0MCpPb0020 =         (TH1D*)directoryPCMPi0pPb0020->Get("TrueFWHMPi0MeV");
	histoPCMMassPi0DatapPb0020->Scale(1000.);
	histoPCMMassPi0MCpPb0020->Scale(1000.);
   
	cout << "20-40%" << endl;
	directoryPCMPi0pPb2040 = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_20-40%"); 
	histoPCMNumberOfEvents2040= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV20-40%");
	histoPCMYieldPi0pPb2040 = 				(TH1D*)directoryPCMPi0pPb2040->Get(nameHistoPCM.Data());
	graphPCMYieldPi0SysErrpPb2040= 	(TGraphAsymmErrors*)directoryPCMPi0pPb2040->Get(nameGraphPCM.Data());	
	graphPCMYieldPi0SysErrRAApPb2040= 	(TGraphAsymmErrors*)directoryPCMPi0pPb2040->Get("Pi0SystErrorA");		
	nPCMEventpPb2040 = 							histoPCMNumberOfEvents2040->GetBinContent(1);
	histoPCMMassPi0DatapPb2040 =        (TH1D*)directoryPCMPi0pPb2040->Get("MassPi0");
	histoPCMMassPi0MCpPb2040 =          (TH1D*)directoryPCMPi0pPb2040->Get("TrueMassPi0");
	histoPCMWidthPi0DatapPb2040 =       (TH1D*)directoryPCMPi0pPb2040->Get("FWHMPi0MeV");
	histoPCMWidthPi0MCpPb2040 =         (TH1D*)directoryPCMPi0pPb2040->Get("TrueFWHMPi0MeV");
	histoPCMMassPi0DatapPb2040->Scale(1000.);
	histoPCMMassPi0MCpPb2040->Scale(1000.);
		
	cout << "40-60%" << endl;
	directoryPCMPi0pPb4060 = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_40-60%"); 
	histoPCMNumberOfEvents4060= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV40-60%");
	histoPCMYieldPi0pPb4060 = 				(TH1D*)directoryPCMPi0pPb4060->Get(nameHistoPCM.Data());
	graphPCMYieldPi0SysErrpPb4060= 	(TGraphAsymmErrors*)directoryPCMPi0pPb4060->Get(nameGraphPCM.Data());	
	graphPCMYieldPi0SysErrRAApPb4060= 	(TGraphAsymmErrors*)directoryPCMPi0pPb4060->Get("Pi0SystErrorA");	
	nPCMEventpPb4060 = 							histoPCMNumberOfEvents4060->GetBinContent(1);
	histoPCMMassPi0DatapPb4060 =        (TH1D*)directoryPCMPi0pPb4060->Get("MassPi0");
	histoPCMMassPi0MCpPb4060 =          (TH1D*)directoryPCMPi0pPb4060->Get("TrueMassPi0");
	histoPCMWidthPi0DatapPb4060 =       (TH1D*)directoryPCMPi0pPb4060->Get("FWHMPi0MeV");
	histoPCMWidthPi0MCpPb4060 =         (TH1D*)directoryPCMPi0pPb4060->Get("TrueFWHMPi0MeV");
	histoPCMMassPi0DatapPb4060->Scale(1000.);
	histoPCMMassPi0MCpPb4060->Scale(1000.);

 //   cout << "60-80%" << endl;
// 	directoryPCMPi0pPb6080 = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_60-80%"); 
// 	histoPCMNumberOfEvents6080= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV60-80%");
// 	histoPCMYieldPi0pPb6080 = 				(TH1D*)directoryPCMPi0pPb6080->Get(nameHistoPCM.Data());
// 	graphPCMYieldPi0SysErrpPb6080= 	(TGraphAsymmErrors*)directoryPCMPi0pPb6080->Get(nameGraphPCM.Data());	
// 	graphPCMYieldPi0SysErrRAApPb6080= 	(TGraphAsymmErrors*)directoryPCMPi0pPb6080->Get("Pi0SystErrorA");	
// 	histoPCMMassPi0DatapPb6080 = 			(TH1D*)directoryPCMPi0pPb6080->Get("MassPi0");
// 	histoPCMMassPi0MCpPb6080 = 			(TH1D*)directoryPCMPi0pPb6080->Get("TrueMassPi0");
// 	histoPCMWidthPi0DatapPb6080 = 		(TH1D*)directoryPCMPi0pPb6080->Get("FWHMPi0MeV");
// 	histoPCMWidthPi0MCpPb6080 = 			(TH1D*)directoryPCMPi0pPb6080->Get("TrueFWHMPi0MeV");
// 	histoPCMMassPi0DatapPb6080->Scale(1000.);
// 	histoPCMMassPi0MCpPb6080->Scale(1000.);
// 	nPCMEventpPb6080 = 							histoPCMNumberOfEvents6080->GetBinContent(1);
	
//    cout << "80-100%" << endl;
//    directoryPCMPi0pPb80100 =            (TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_80-100%"); 
//    histoPCMNumberOfEvents80100=         (TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV80-100%");
//    histoPCMYieldPi0pPb80100 =           (TH1D*)directoryPCMPi0pPb80100->Get(nameHistoPCM.Data());
//    graphPCMYieldPi0SysErrpPb80100=   (TGraphAsymmErrors*)directoryPCMPi0pPb80100->Get(nameGraphPCM.Data());   
//    graphPCMYieldPi0SysErrRAApPb80100=   (TGraphAsymmErrors*)directoryPCMPi0pPb80100->Get("Pi0SystErrorA"); 
//    histoPCMMassPi0DatapPb80100 =        (TH1D*)directoryPCMPi0pPb80100->Get("MassPi0");
//    histoPCMMassPi0MCpPb80100 =          (TH1D*)directoryPCMPi0pPb80100->Get("TrueMassPi0");
//    histoPCMWidthPi0DatapPb80100 =       (TH1D*)directoryPCMPi0pPb80100->Get("FWHMPi0MeV");
//    histoPCMWidthPi0MCpPb80100 =         (TH1D*)directoryPCMPi0pPb80100->Get("TrueFWHMPi0MeV");
//    histoPCMMassPi0DatapPb80100->Scale(1000.);
//    histoPCMMassPi0MCpPb80100->Scale(1000.);
//    nPCMEventpPb80100 =                     histoPCMNumberOfEvents80100->GetBinContent(1);
   
	cout << "60-100%" << endl;
	directoryPCMPi0pPb60100 =            (TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_60-100%"); 
	histoPCMNumberOfEvents60100=         (TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV60-100%");
	histoPCMYieldPi0pPb60100 =           (TH1D*)directoryPCMPi0pPb60100->Get(nameHistoPCM.Data());
	graphPCMYieldPi0SysErrpPb60100=   (TGraphAsymmErrors*)directoryPCMPi0pPb60100->Get(nameGraphPCM.Data());   
	graphPCMYieldPi0SysErrRAApPb60100=   (TGraphAsymmErrors*)directoryPCMPi0pPb60100->Get("Pi0SystErrorA"); 
	histoPCMMassPi0DatapPb60100 =        (TH1D*)directoryPCMPi0pPb60100->Get("MassPi0");
	histoPCMMassPi0MCpPb60100 =          (TH1D*)directoryPCMPi0pPb60100->Get("TrueMassPi0");
	histoPCMWidthPi0DatapPb60100 =       (TH1D*)directoryPCMPi0pPb60100->Get("FWHMPi0MeV");
	histoPCMWidthPi0MCpPb60100 =         (TH1D*)directoryPCMPi0pPb60100->Get("TrueFWHMPi0MeV");
	histoPCMMassPi0DatapPb60100->Scale(1000.);
	histoPCMMassPi0MCpPb60100->Scale(1000.);
	nPCMEventpPb60100 =                     histoPCMNumberOfEvents60100->GetBinContent(1);	
	//****************************************************************************************************
	//************************** Read Eta data for PCM *******************************************************
	//****************************************************************************************************
	cout << "PCM Eta" << endl;
	cout << "0-100%" << endl;
	directoryPCMEtapPb = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_0-100%"); 
	histoPCMYieldEtapPb = 				(TH1D*)directoryPCMEtapPb->Get(nameHistoPCMEta.Data());
	histoPCMRAWYieldEtapPb =            (TH1D*)directoryPCMEtapPb->Get("RawYieldEta");
	graphPCMYieldEtaSysErrpPb= 	(TGraphAsymmErrors*)directoryPCMEtapPb->Get(nameGraphPCMEta.Data());	
	graphPCMYieldEtaSysErrRAApPb= 	(TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtaSystErrorA");	
	histoPCMEtaPi0RatiopPb =            (TH1D*)directoryPCMEtapPb->Get("EtatoPi0Ratio");
	graphPCMEtaPi0RatioSysErrpPb=    (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtatoPi0RatioSys"); 
	histoPCMMassEtaDatapPb = 			(TH1D*)directoryPCMEtapPb->Get("MassEta");
	histoPCMMassEtaMCpPb = 			(TH1D*)directoryPCMEtapPb->Get("TrueMassEta");
	histoPCMWidthEtaDatapPb = 		(TH1D*)directoryPCMEtapPb->Get("FWHMEtaMeV");
	histoPCMWidthEtaMCpPb = 			(TH1D*)directoryPCMEtapPb->Get("TrueFWHMEtaMeV");
	histoPCMMassEtaDatapPb->Scale(1000.);
	histoPCMMassEtaMCpPb->Scale(1000.);
   
	cout << "0-20%" << endl;
	directoryPCMEtapPb0020 = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_0-20%"); 
	histoPCMYieldEtapPb0020 = 				(TH1D*)directoryPCMEtapPb0020->Get(nameHistoPCMEta.Data());
	graphPCMYieldEtaSysErrpPb0020= 	(TGraphAsymmErrors*)directoryPCMEtapPb0020->Get(nameGraphPCMEta.Data());	
	graphPCMYieldEtaSysErrRAApPb0020= 	(TGraphAsymmErrors*)directoryPCMEtapPb0020->Get("EtaSystErrorA");
	histoPCMEtaPi0RatiopPb0020 =            (TH1D*)directoryPCMEtapPb0020->Get("EtatoPi0Ratio");
	graphPCMEtaPi0RatioSysErrpPb0020=    (TGraphAsymmErrors*)directoryPCMEtapPb0020->Get("EtatoPi0RatioSys"); 
	histoPCMMassEtaDatapPb0020 =        (TH1D*)directoryPCMEtapPb0020->Get("MassEta");
	histoPCMMassEtaMCpPb0020 =          (TH1D*)directoryPCMEtapPb0020->Get("TrueMassEta");
	histoPCMWidthEtaDatapPb0020 =       (TH1D*)directoryPCMEtapPb0020->Get("FWHMEtaMeV");
	histoPCMWidthEtaMCpPb0020 =         (TH1D*)directoryPCMEtapPb0020->Get("TrueFWHMEtaMeV");
	histoPCMMassEtaDatapPb0020->Scale(1000.);
	histoPCMMassEtaMCpPb0020->Scale(1000.);
   
	cout << "20-40%" << endl;
	directoryPCMEtapPb2040 = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_20-40%"); 
	histoPCMYieldEtapPb2040 = 				(TH1D*)directoryPCMEtapPb2040->Get(nameHistoPCMEta.Data());
	graphPCMYieldEtaSysErrpPb2040= 	(TGraphAsymmErrors*)directoryPCMEtapPb2040->Get(nameGraphPCMEta.Data());	
	graphPCMYieldEtaSysErrRAApPb2040= 	(TGraphAsymmErrors*)directoryPCMEtapPb2040->Get("EtaSystErrorA");		
	histoPCMEtaPi0RatiopPb2040 =            (TH1D*)directoryPCMEtapPb2040->Get("EtatoPi0Ratio");
	graphPCMEtaPi0RatioSysErrpPb2040=    (TGraphAsymmErrors*)directoryPCMEtapPb2040->Get("EtatoPi0RatioSys"); 
	histoPCMMassEtaDatapPb2040 =        (TH1D*)directoryPCMEtapPb2040->Get("MassEta");
	histoPCMMassEtaMCpPb2040 =          (TH1D*)directoryPCMEtapPb2040->Get("TrueMassEta");
	histoPCMWidthEtaDatapPb2040 =       (TH1D*)directoryPCMEtapPb2040->Get("FWHMEtaMeV");
	histoPCMWidthEtaMCpPb2040 =         (TH1D*)directoryPCMEtapPb2040->Get("TrueFWHMEtaMeV");
	histoPCMMassEtaDatapPb2040->Scale(1000.);
	histoPCMMassEtaMCpPb2040->Scale(1000.);
		
	cout << "40-60%" << endl;
	directoryPCMEtapPb4060 = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_40-60%"); 
	histoPCMYieldEtapPb4060 = 				(TH1D*)directoryPCMEtapPb4060->Get(nameHistoPCMEta.Data());
	graphPCMYieldEtaSysErrpPb4060= 	(TGraphAsymmErrors*)directoryPCMEtapPb4060->Get(nameGraphPCMEta.Data());	
	graphPCMYieldEtaSysErrRAApPb4060= 	(TGraphAsymmErrors*)directoryPCMEtapPb4060->Get("EtaSystErrorA");	
	histoPCMEtaPi0RatiopPb4060 =            (TH1D*)directoryPCMEtapPb4060->Get("EtatoPi0Ratio");
	graphPCMEtaPi0RatioSysErrpPb4060=    (TGraphAsymmErrors*)directoryPCMEtapPb4060->Get("EtatoPi0RatioSys"); 
	histoPCMMassEtaDatapPb4060 =        (TH1D*)directoryPCMEtapPb4060->Get("MassEta");
	histoPCMMassEtaMCpPb4060 =          (TH1D*)directoryPCMEtapPb4060->Get("TrueMassEta");
	histoPCMWidthEtaDatapPb4060 =       (TH1D*)directoryPCMEtapPb4060->Get("FWHMEtaMeV");
	histoPCMWidthEtaMCpPb4060 =         (TH1D*)directoryPCMEtapPb4060->Get("TrueFWHMEtaMeV");
	histoPCMMassEtaDatapPb4060->Scale(1000.);
	histoPCMMassEtaMCpPb4060->Scale(1000.);

//    cout << "60-80%" << endl;
// 	directoryPCMEtapPb6080 = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_60-80%"); 
// 	histoPCMYieldEtapPb6080 = 				(TH1D*)directoryPCMEtapPb6080->Get(nameHistoPCMEta.Data());
// 	graphPCMYieldEtaSysErrpPb6080= 	(TGraphAsymmErrors*)directoryPCMEtapPb6080->Get(nameGraphPCM.Data());	
// 	graphPCMYieldEtaSysErrRAApPb6080= 	(TGraphAsymmErrors*)directoryPCMEtapPb6080->Get("EtaSystErrorA");	
//    histoPCMEtaPi0RatiopPb6080 =            (TH1D*)directoryPCMEtapPb6080->Get("EtatoPi0Ratio");
//    graphPCMEtaPi0RatioSysErrpPb6080=    (TGraphAsymmErrors*)directoryPCMEtapPb6080->Get("EtatoPi0RatioSys"); 
// 	histoPCMMassEtaDatapPb6080 = 			(TH1D*)directoryPCMEtapPb6080->Get("MassEta");
// 	histoPCMMassEtaMCpPb6080 = 			(TH1D*)directoryPCMEtapPb6080->Get("TrueMassEta");
// 	histoPCMWidthEtaDatapPb6080 = 		(TH1D*)directoryPCMEtapPb6080->Get("FWHMEtaMeV");
// 	histoPCMWidthEtaMCpPb6080 = 			(TH1D*)directoryPCMEtapPb6080->Get("TrueFWHMEtaMeV");
// 	histoPCMMassEtaDatapPb6080->Scale(1000.);
// 	histoPCMMassEtaMCpPb6080->Scale(1000.);

	
//    cout << "80-100%" << endl;
//    directoryPCMEtapPb80100 =            (TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_80-100%"); 
//    histoPCMYieldEtapPb80100 =           (TH1D*)directoryPCMEtapPb80100->Get(nameHistoPCM.Data());
//    graphPCMYieldEtaSysErrpPb80100=   (TGraphAsymmErrors*)directoryPCMEtapPb80100->Get(nameGraphPCM.Data());   
//    graphPCMYieldEtaSysErrRAApPb80100=   (TGraphAsymmErrors*)directoryPCMEtapPb80100->Get("EtaSystErrorA"); 
//    histoPCMEtaPi0RatiopPb80100 =            (TH1D*)directoryPCMEtapPb80100->Get("EtatoPi0Ratio");
//    graphPCMEtaPi0RatioSysErrpPb80100=    (TGraphAsymmErrors*)directoryPCMEtapPb80100->Get("EtatoPi0RatioSys"); 
//    histoPCMMassEtaDatapPb80100 =        (TH1D*)directoryPCMEtapPb80100->Get("MassEta");
//    histoPCMMassEtaMCpPb80100 =          (TH1D*)directoryPCMEtapPb80100->Get("TrueMassEta");
//    histoPCMWidthEtaDatapPb80100 =       (TH1D*)directoryPCMEtapPb80100->Get("FWHMEtaMeV");
//    histoPCMWidthEtaMCpPb80100 =         (TH1D*)directoryPCMEtapPb80100->Get("TrueFWHMEtaMeV");
//    histoPCMMassEtaDatapPb80100->Scale(1000.);
//    histoPCMMassEtaMCpPb80100->Scale(1000.);

   
	cout << "60-100%" << endl;
	directoryPCMEtapPb60100 =            (TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_60-100%"); 
	histoPCMYieldEtapPb60100 =           (TH1D*)directoryPCMEtapPb60100->Get(nameHistoPCMEta.Data());
	graphPCMYieldEtaSysErrpPb60100=   (TGraphAsymmErrors*)directoryPCMEtapPb60100->Get(nameGraphPCMEta.Data());   
	graphPCMYieldEtaSysErrRAApPb60100=   (TGraphAsymmErrors*)directoryPCMEtapPb60100->Get("EtaSystErrorA"); 
	histoPCMEtaPi0RatiopPb60100 =            (TH1D*)directoryPCMEtapPb60100->Get("EtatoPi0Ratio");
	graphPCMEtaPi0RatioSysErrpPb60100=    (TGraphAsymmErrors*)directoryPCMEtapPb60100->Get("EtatoPi0RatioSys"); 
	histoPCMMassEtaDatapPb60100 =        (TH1D*)directoryPCMEtapPb60100->Get("MassEta");
	histoPCMMassEtaMCpPb60100 =          (TH1D*)directoryPCMEtapPb60100->Get("TrueMassEta");
	histoPCMWidthEtaDatapPb60100 =       (TH1D*)directoryPCMEtapPb60100->Get("FWHMEtaMeV");
	histoPCMWidthEtaMCpPb60100 =         (TH1D*)directoryPCMEtapPb60100->Get("TrueFWHMEtaMeV");
	histoPCMMassEtaDatapPb60100->Scale(1000.);
	histoPCMMassEtaMCpPb60100->Scale(1000.);
		
	//****************************************************************************************************
	//************************** Read data for EMCAL *****************************************************
	//****************************************************************************************************
	cout << "EMCAL" << endl;
	cout << "0-100%" << endl;
	fileEMCALpPb = 		new TFile(fileNameCaloEMCALpPb);
	directoryPi0EMCALpPb =       (TDirectory*)fileEMCALpPb->Get("Pi0HI00-100%"); 
	histoEMCALNumberOfEventspPb= (TH1D*)fileEMCALpPb->Get("histoNumberOfEventsHI00-100%");
	histoEMCALRawYieldPi0pPb =      (TH1D*)directoryPi0EMCALpPb->Get("YieldPi0");
	histoEMCALMassPi0DatapPb =      (TH1D*)directoryPi0EMCALpPb->Get("MassPi0");
	histoEMCALWidthPi0DatapPb =     (TH1D*)directoryPi0EMCALpPb->Get("SigmaPi0");
	histoEMCALMassPi0DatapPb->Scale(1000.);
	histoEMCALMassPi0DatapPb->SetBinContent(histoEMCALMassPi0DatapPb->GetNbinsX(),0);
	histoEMCALMassPi0DatapPb->SetBinContent(histoEMCALMassPi0DatapPb->GetNbinsX()-1,0);
	histoEMCALWidthPi0DatapPb->SetBinContent(histoEMCALWidthPi0DatapPb->GetNbinsX(),10000.);
	histoEMCALWidthPi0DatapPb->SetBinContent(histoEMCALWidthPi0DatapPb->GetNbinsX()-1,10000.);
	nEMCALEventpPb =  histoEMCALNumberOfEventspPb->GetBinContent(1);
	histoEMCALRawYieldPi0pPb->Scale(1./nEMCALEventpPb);

	cout << "0-20%" << endl;
	directoryEMCALPi0pPb0020 =         (TDirectory*)fileEMCALpPb->Get("Pi0HI00-20%"); 
	histoEMCALNumberOfEvents0020=      (TH1D*)fileEMCALpPb->Get("histoNumberOfEventsHI00-20%");
	nEMCALEventpPb0020 =               histoEMCALNumberOfEvents0020->GetBinContent(1);
	histoEMCALRawYieldPi0pPb0020 =        (TH1D*)directoryEMCALPi0pPb0020->Get("YieldPi0");
	histoEMCALRawYieldPi0pPb0020->Scale(1./nEMCALEventpPb0020);
	histoEMCALMassPi0DatapPb0020 =        (TH1D*)directoryEMCALPi0pPb0020->Get("MassPi0");
	histoEMCALWidthPi0DatapPb0020 =       (TH1D*)directoryEMCALPi0pPb0020->Get("SigmaPi0");
	histoEMCALMassPi0DatapPb0020->Scale(1000.);
   
	cout << "20-40%" << endl;
	directoryEMCALPi0pPb2040 =         (TDirectory*)fileEMCALpPb->Get("Pi0HI20-40%"); 
	histoEMCALNumberOfEvents2040=      (TH1D*)fileEMCALpPb->Get("histoNumberOfEventsHI20-40%");
	nEMCALEventpPb2040 =               histoEMCALNumberOfEvents2040->GetBinContent(1);
	histoEMCALRawYieldPi0pPb2040 =        (TH1D*)directoryEMCALPi0pPb2040->Get("YieldPi0");
	histoEMCALRawYieldPi0pPb2040->Scale(1./nEMCALEventpPb2040);
	histoEMCALMassPi0DatapPb2040 =        (TH1D*)directoryEMCALPi0pPb2040->Get("MassPi0");
	histoEMCALWidthPi0DatapPb2040 =       (TH1D*)directoryEMCALPi0pPb2040->Get("SigmaPi0");
	histoEMCALMassPi0DatapPb2040->Scale(1000.);
   
	cout << "40-60%" << endl;
	directoryEMCALPi0pPb4060 =            (TDirectory*)fileEMCALpPb->Get("Pi0HI40-60%"); 
	histoEMCALNumberOfEvents4060=         (TH1D*)fileEMCALpPb->Get("histoNumberOfEventsHI40-60%");
	nEMCALEventpPb4060 =                     histoEMCALNumberOfEvents4060->GetBinContent(1);
	histoEMCALRawYieldPi0pPb4060 =        (TH1D*)directoryEMCALPi0pPb4060->Get("YieldPi0");
	histoEMCALRawYieldPi0pPb4060->Scale(1./nEMCALEventpPb4060);
	histoEMCALMassPi0DatapPb4060 =        (TH1D*)directoryEMCALPi0pPb4060->Get("MassPi0");
	histoEMCALWidthPi0DatapPb4060 =       (TH1D*)directoryEMCALPi0pPb4060->Get("SigmaPi0");
	histoEMCALMassPi0DatapPb4060->Scale(1000.);
      
	cout << "60-80%" << endl;
	directoryEMCALPi0pPb6080 =            (TDirectory*)fileEMCALpPb->Get("Pi0HI60-80%"); 
	histoEMCALNumberOfEvents6080=         (TH1D*)fileEMCALpPb->Get("histoNumberOfEventsHI60-80%");
	histoEMCALMassPi0DatapPb6080 =        (TH1D*)directoryEMCALPi0pPb6080->Get("MassPi0");
	histoEMCALWidthPi0DatapPb6080 =       (TH1D*)directoryEMCALPi0pPb6080->Get("FWHMPi0MeV");
	histoEMCALMassPi0DatapPb6080->Scale(1000.);
	nEMCALEventpPb6080 =                     histoEMCALNumberOfEvents6080->GetBinContent(1);
	histoEMCALRawYieldPi0pPb6080 =        (TH1D*)directoryEMCALPi0pPb6080->Get("YieldPi0");
	histoEMCALRawYieldPi0pPb6080->Scale(1./nEMCALEventpPb6080);

	cout << "80-100%" << endl;
	directoryEMCALPi0pPb80100 =            (TDirectory*)fileEMCALpPb->Get("Pi0HI80-100%"); 
	histoEMCALNumberOfEvents80100=         (TH1D*)fileEMCALpPb->Get("histoNumberOfEventsHI80-100%");
	histoEMCALMassPi0DatapPb80100 =        (TH1D*)directoryEMCALPi0pPb80100->Get("MassPi0");
	histoEMCALWidthPi0DatapPb80100 =       (TH1D*)directoryEMCALPi0pPb80100->Get("FWHMPi0MeV");
	histoEMCALMassPi0DatapPb80100->Scale(1000.);
	nEMCALEventpPb80100 =                     histoEMCALNumberOfEvents80100->GetBinContent(1);
	histoEMCALRawYieldPi0pPb80100 =        (TH1D*)directoryEMCALPi0pPb80100->Get("YieldPi0");
	histoEMCALRawYieldPi0pPb80100->Scale(1./nEMCALEventpPb80100);
  
	//****************************************************************************************************
	//************************** Read data for PHOS ******************************************************
	//****************************************************************************************************
   
	cout << "PHOS" << endl;
	cout << "0-100%" << endl;

	filePHOSpPb = 					new TFile(fileNameCaloPHOSpPb);
	directoryPHOSPi0pPb = 				(TDirectory*)filePHOSpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
	histoPHOSYieldPi0pPb = 				(TH1D*)directoryPHOSPi0pPb->Get(nameHistoPCM.Data());
 
	//   filePHOSpPb =       new TFile(fileNameCaloPHOSpPb);
	//histoPHOSYieldPi0pPb = (TH1D*)filePHOSpPb->Get("yeild1_GS_Disp_cen0");
	
	//****************************************************************************************************
	//************************** Combine & Compare PCM & PHOS minBias spectrum ***************************
	//****************************************************************************************************
	TH1D* ratioYieldPi0pPb_PHOS_PCM = new TH1D("ratioYieldPi0pPb_PHOS_PCM", "ratioYieldPi0pPb_PHOS_PCM", 23, ptBinningPHOSpPb);
	for (Int_t i = 1; i < ratioYieldPi0pPb_PHOS_PCM->GetNbinsX()+1; i++){
		Int_t PCMbin = i+5;
		Int_t PHOSbin = i;
		Double_t ratio = histoPHOSYieldPi0pPb->GetBinContent(PHOSbin)/histoPCMYieldPi0pPb->GetBinContent(PCMbin);
		Double_t ratioErr =  TMath::Sqrt( pow(histoPHOSYieldPi0pPb->GetBinError(PHOSbin)/histoPCMYieldPi0pPb->GetBinContent(PCMbin),2)  + pow( histoPCMYieldPi0pPb->GetBinError(PCMbin)*histoPHOSYieldPi0pPb->GetBinContent(PHOSbin)/pow(histoPCMYieldPi0pPb->GetBinContent(PCMbin),2),2) );
		ratioYieldPi0pPb_PHOS_PCM->SetBinContent(i,ratio);
		ratioYieldPi0pPb_PHOS_PCM->SetBinError(i,ratioErr);
	}
      
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
		
   //***********************************************************************************************************************************************   
	//*********************************** Plotting raw EMCAL yields Pi0 *****************************************************************************
   //***********************************************************************************************************************************************   
	TCanvas* canvasRawSpectraEMCALAllCent = new TCanvas("canvasRawSpectraEMCALAllCent","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasRawSpectraEMCALAllCent,  0.13, 0.01, 0.015, 0.08);
 	
	canvasRawSpectraEMCALAllCent->SetLogy();
  	canvasRawSpectraEMCALAllCent->SetLogx();
	TH2F * histo2DEMCALRawSpectraAllCent;
	histo2DEMCALRawSpectraAllCent = new TH2F("histo2DEMCALRawSpectraAllCent","histo2DEMCALRawSpectraAllCent",1000,0.23,30.,1000,1e-8,7e0	);
	SetStyleHistoTH2ForGraphs(histo2DEMCALRawSpectraAllCent, "p_{T} (GeV/c)","Raw Yield/Evt", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}
 	histo2DEMCALRawSpectraAllCent->GetXaxis()->SetLabelOffset(-0.01);
	histo2DEMCALRawSpectraAllCent->GetYaxis()->SetLabelOffset(0.01);
	histo2DEMCALRawSpectraAllCent->DrawCopy(); 
	
	DrawGammaSetMarker(histoEMCALRawYieldPi0pPb0020, markerStyleCommmonSpectrum0020,markerSizeCommonSpectrum0020, colorComb0020 , colorComb0020);
	histoEMCALRawYieldPi0pPb0020->Draw("p,same,e1");
	DrawGammaSetMarker(histoEMCALRawYieldPi0pPb2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
	histoEMCALRawYieldPi0pPb2040->Draw("p,same,e1");
	DrawGammaSetMarker(histoEMCALRawYieldPi0pPb4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
	histoEMCALRawYieldPi0pPb4060->Draw("p,same,e1");
	DrawGammaSetMarker(histoEMCALRawYieldPi0pPb6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
	histoEMCALRawYieldPi0pPb6080->Draw("p,same,e1");
	DrawGammaSetMarker(histoEMCALRawYieldPi0pPb80100, markerStyleCommmonSpectrum80100,markerSizeCommonSpectrum80100, colorComb80100 , colorComb80100);
	histoEMCALRawYieldPi0pPb80100->Draw("p,same,e1");
   
	TLatex *labelSpectraPi0LabelEMCALpPb = new TLatex(0.75,0.92,"#pi^{0} #rightarrow #gamma #gamma (EMCAL)");
	SetStyleTLatex( labelSpectraPi0LabelEMCALpPb, 0.03,4);
	labelSpectraPi0LabelEMCALpPb->Draw();
	
	TLegend* legendSpectraEMCAL = new TLegend(0.16,0.11,0.5,0.3);
	legendSpectraEMCAL->SetFillColor(0);
	legendSpectraEMCAL->SetLineColor(0);
	legendSpectraEMCAL->SetTextSize(0.025);
	legendSpectraEMCAL->AddEntry(histoEMCALRawYieldPi0pPb0020,collisionSystem0020.Data(),"pf");
	legendSpectraEMCAL->AddEntry(histoEMCALRawYieldPi0pPb2040,collisionSystem2040.Data(),"pf");
	legendSpectraEMCAL->AddEntry(histoEMCALRawYieldPi0pPb4060,collisionSystem4060.Data(),"pf");
	legendSpectraEMCAL->AddEntry(histoEMCALRawYieldPi0pPb6080,collisionSystem6080.Data(),"pf");
	legendSpectraEMCAL->AddEntry(histoEMCALRawYieldPi0pPb80100,collisionSystem80100.Data(),"pf");
	legendSpectraEMCAL->Draw();
	
	canvasRawSpectraEMCALAllCent->Update();
	canvasRawSpectraEMCALAllCent->Print(Form("%s/Pi0_RawSpectra_EMCAL.%s",outputDir.Data(),suffix.Data()));

	//**************************************************************************************************************************************************
	//************************************* Plotting corrected spectra all centralities PCM Pi0 ********************************************************
	//**************************************************************************************************************************************************
	TCanvas* canvasSpectraPCMPi0AllCent = new TCanvas("canvasSpectraPCMPi0AllCent","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasSpectraPCMPi0AllCent,  0.13, 0.01, 0.015, 0.08);

	canvasSpectraPCMPi0AllCent->SetLogy();
	canvasSpectraPCMPi0AllCent->SetLogx();
	TH2F * histo2DPCMSpectraPi0AllCent;
	histo2DPCMSpectraPi0AllCent = new TH2F("histo2DPCMSpectraPi0AllCent","histo2DPCMSpectraPi0AllCent",1000,0.23,30.,1000,1e-8,7e1   );
	SetStyleHistoTH2ForGraphs(histo2DPCMSpectraPi0AllCent, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}
	histo2DPCMSpectraPi0AllCent->GetXaxis()->SetLabelOffset(-0.01);
	histo2DPCMSpectraPi0AllCent->GetYaxis()->SetLabelOffset(0.01);
	histo2DPCMSpectraPi0AllCent->DrawCopy(); 

	DrawGammaSetMarker(histoPCMYieldPi0pPb0020, markerStyleCommmonSpectrum0020,markerSizeCommonSpectrum0020, colorComb0020 , colorComb0020);
	histoPCMYieldPi0pPb0020->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMYieldPi0pPb2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
	histoPCMYieldPi0pPb2040->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMYieldPi0pPb4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
	histoPCMYieldPi0pPb4060->Draw("p,same,e1");
	//    DrawGammaSetMarker(histoPCMYieldPi0pPb6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
	//    histoPCMYieldPi0pPb6080->Draw("p,same,e1");
	//    DrawGammaSetMarker(histoPCMYieldPi0pPb80100, markerStyleCommmonSpectrum80100,markerSizeCommonSpectrum80100, colorComb80100 , colorComb80100);
	//    histoPCMYieldPi0pPb80100->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMYieldPi0pPb60100, markerStyleCommmonSpectrum60100,markerSizeCommonSpectrum60100, colorComb60100 , colorComb60100);
	histoPCMYieldPi0pPb60100->Draw("p,same,e1"); 
	DrawGammaSetMarker(histoPCMYieldPi0pPb, 34,markerSizeCommonSpectrum60100, colorPCM , colorPCM);
	histoPCMYieldPi0pPb->Draw("p,same,e1");  

	TLatex *labelSpectraPi0LabelPCMpPb = new TLatex(0.65,0.92,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-} (PCM)");
	SetStyleTLatex( labelSpectraPi0LabelPCMpPb, 0.03,4);
	labelSpectraPi0LabelPCMpPb->Draw();
	
	TLegend* legendSpectraPi0PCM = new TLegend(0.16,0.11,0.5,0.3);
	legendSpectraPi0PCM->SetFillColor(0);
	legendSpectraPi0PCM->SetLineColor(0);
	legendSpectraPi0PCM->SetTextSize(0.025);
	legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb0020,collisionSystem0020.Data(),"pf");
	legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb2040,collisionSystem2040.Data(),"pf");
	legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb4060,collisionSystem4060.Data(),"pf");  
	// legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb6080,collisionSystem6080.Data(),"pf");
	legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb60100,collisionSystem60100.Data(),"pf");
	legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb,collisionSystempPb.Data(),"pf");
	//   legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb6080,collisionSystem6080.Data(),"pf");
	//  legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb80100,collisionSystem80100.Data(),"pf");
	legendSpectraPi0PCM->Draw();

	canvasSpectraPCMPi0AllCent->Update();
	canvasSpectraPCMPi0AllCent->Print(Form("%s/Pi0_Spectra_PCM.%s",outputDir.Data(),suffix.Data()));
   
    //**************************************************************************************************************************************************
	//************************************* Plotting corrected spectra all centralities PCM Eta ********************************************************
	//**************************************************************************************************************************************************
	TCanvas* canvasSpectraPCMEtaAllCent = new TCanvas("canvasSpectraPCMEtaAllCent","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasSpectraPCMEtaAllCent,  0.13, 0.01, 0.015, 0.08);
	cout <<"hier   !!!!!!!!!!!"<<endl;   
	canvasSpectraPCMEtaAllCent->SetLogy();
	canvasSpectraPCMEtaAllCent->SetLogx();
	TH2F * histo2DPCMSpectraEtaAllCent;
	histo2DPCMSpectraEtaAllCent = new TH2F("histo2DPCMSpectraEtaAllCent","histo2DPCMSpectraEtaAllCent",1000,0.8,10.,1000,1e-7,4e-1   );
	SetStyleHistoTH2ForGraphs(histo2DPCMSpectraEtaAllCent, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}
	histo2DPCMSpectraEtaAllCent->GetXaxis()->SetLabelOffset(-0.01);
	histo2DPCMSpectraEtaAllCent->GetYaxis()->SetLabelOffset(0.01);
	histo2DPCMSpectraEtaAllCent->DrawCopy(); 
	cout <<"hier   !!!!!!!!!!!"<<endl;	      
	DrawGammaSetMarker(histoPCMYieldEtapPb0020, markerStyleCommmonSpectrum0020,markerSizeCommonSpectrum0020, colorComb0020 , colorComb0020);
	histoPCMYieldEtapPb0020->Draw("p,same,e1");    cout <<"hier   !!!!!!!!!!!"<<endl;
	DrawGammaSetMarker(histoPCMYieldEtapPb2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
	histoPCMYieldEtapPb2040->Draw("p,same,e1");    cout <<"hier   !!!!!!!!!!!"<<endl;
	DrawGammaSetMarker(histoPCMYieldEtapPb4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
	histoPCMYieldEtapPb4060->Draw("p,same,e1");    cout <<"hier   !!!!!!!!!!!"<<endl;
// 	DrawGammaSetMarker(histoPCMYieldEtapPb6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
// 	histoPCMYieldEtapPb6080->Draw("p,same,e1");    cout <<"hier   !!!!!!!!!!!"<<endl;
// 	DrawGammaSetMarker(histoPCMYieldEtapPb80100, markerStyleCommmonSpectrum80100,markerSizeCommonSpectrum80100, colorComb80100 , colorComb80100);
// 	histoPCMYieldEtapPb80100->Draw("p,same,e1");
// 	cout <<"hier   !!!!!!!!!!!"<<endl;	      
    DrawGammaSetMarker(histoPCMYieldEtapPb60100, markerStyleCommmonSpectrum60100,markerSizeCommonSpectrum60100, colorComb60100 , colorComb60100);
    histoPCMYieldEtapPb60100->Draw("p,same,e1");  
	DrawGammaSetMarker(histoPCMYieldEtapPb, 34,markerSizeCommonSpectrum60100, colorPCM , colorPCM);
    histoPCMYieldEtapPb->Draw("p,same,e1");  
   
	TLatex *labelSpectraEtaLabelPCMpPb = new TLatex(0.65,0.92,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-} (PCM)");
	SetStyleTLatex( labelSpectraEtaLabelPCMpPb, 0.03,4);
	labelSpectraEtaLabelPCMpPb->Draw();
	
	TLegend* legendSpectraEtaPCM = new TLegend(0.16,0.11,0.5,0.3);
	legendSpectraEtaPCM->SetFillColor(0);
	legendSpectraEtaPCM->SetLineColor(0);
	legendSpectraEtaPCM->SetTextSize(0.025);
	legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb0020,collisionSystem0020.Data(),"pf");
	legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb2040,collisionSystem2040.Data(),"pf");
	legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb4060,collisionSystem4060.Data(),"pf");  
	//  legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb6080,collisionSystem6080.Data(),"pf");
	legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb60100,collisionSystem60100.Data(),"pf");
	legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb,collisionSystempPb.Data(),"pf");
	//   legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb6080,collisionSystem6080.Data(),"pf");
	//  legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb80100,collisionSystem80100.Data(),"pf");
	legendSpectraEtaPCM->Draw();
	cout <<"hier   !!!!!!!!!!!"<<endl;	      
	canvasSpectraPCMEtaAllCent->Update();
	canvasSpectraPCMEtaAllCent->Print(Form("%s/Eta_Spectra_PCM.%s",outputDir.Data(),suffix.Data()));

		//************************************************************************************************************************************************  	   **
	//************************************* Plotting eta/pi0 ratio all centralities PCM Eta ********************************************************
	//**************************************************************************************************************************************************
	TCanvas* canvasEtaToPi0PCMAllCent = new TCanvas("canvasEtaToPi0PCMAllCent","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasEtaToPi0PCMAllCent,  0.08, 0.01, 0.015, 0.08);
	
	//    canvasEtaToPi0PCMAllCent->SetLogy();
	//    canvasEtaToPi0PCMAllCent->SetLogx();
	TH2F * histo2DEtaToPi0PCMAllCent;
	histo2DEtaToPi0PCMAllCent = new TH2F("histo2DEtaToPi0PCMAllCent","histo2DEtaToPi0PCMAllCent",1000,0.,8.,1000,0,1.1   );
	SetStyleHistoTH2ForGraphs(histo2DEtaToPi0PCMAllCent, "p_{T} (GeV/c)","#eta/#pi^{0}", 0.03,0.04, 0.03,0.04, 0.83,0.95);
	//    histo2DEtaToPi0PCMAllCent->GetXaxis()->SetLabelOffset(-0.01);
	//    histo2DEtaToPi0PCMAllCent->GetYaxis()->SetLabelOffset(0.01);
	histo2DEtaToPi0PCMAllCent->DrawCopy(); 
	cout <<"hier   !!!!!!!!!!!"<<endl;  
	DrawGammaSetMarker(histoPCMEtaPi0RatiopPb0020, markerStyleCommmonSpectrum0020,markerSizeCommonSpectrum0020, colorComb0020 , colorComb0020);
	histoPCMEtaPi0RatiopPb0020->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMEtaPi0RatiopPb2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
	histoPCMEtaPi0RatiopPb2040->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMEtaPi0RatiopPb4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
	histoPCMEtaPi0RatiopPb4060->Draw("p,same,e1");
//    DrawGammaSetMarker(histoPCMEtaPi0RatiopPb6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
//    histoPCMEtaPi0RatiopPb6080->Draw("p,same,e1");
	DrawGammaSetMarker(histoPCMEtaPi0RatiopPb60100, markerStyleCommmonSpectrum60100,markerSizeCommonSpectrum60100, colorComb60100 , colorComb60100);
	histoPCMEtaPi0RatiopPb60100->Draw("p,same,e1");
   
//    DrawGammaSetMarker(histoPCMEtaPi0RatiopPb60100, markerStyleCommmonSpectrum60100,markerSizeCommonSpectrum60100, colorComb60100 , colorComb60100);
//    histoPCMEtaPi0RatiopPb60100->Draw("p,same,e1"); 
   
	DrawGammaSetMarker(histoPCMEtaPi0RatiopPb, 34,markerSizeCommonSpectrum60100, colorPCM , colorPCM);
	histoPCMEtaPi0RatiopPb->Draw("p,same,e1");  
	cout <<"hier   !!!!!!!!!!!"<<endl;   
   
	TLegend* legendEtaToPi0PCM = new TLegend(0.16,0.76,0.5,0.95);
	legendEtaToPi0PCM->SetFillColor(0);
	legendEtaToPi0PCM->SetLineColor(0);
	legendEtaToPi0PCM->SetTextSize(0.025);
	legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb0020,collisionSystem0020.Data(),"pf");
	legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb2040,collisionSystem2040.Data(),"pf");
	legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb4060,collisionSystem4060.Data(),"pf"); 
	//  legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb6080,collisionSystem6080.Data(),"pf");
	legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb60100,collisionSystem60100.Data(),"pf");
	legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb,collisionSystempPb.Data(),"pf");
	//   legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb6080,collisionSystem6080.Data(),"pf");
	//  legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb80100,collisionSystem80100.Data(),"pf");
	legendEtaToPi0PCM->Draw();
	
	canvasEtaToPi0PCMAllCent->Update();
	canvasEtaToPi0PCMAllCent->Print(Form("%s/EtaToPi0_Ratio_PCM.%s",outputDir.Data(),suffix.Data()));
	
	cout <<"hier   !!!!!!!!!!!"<<endl;  
	//**************************************************************************************************************************************************
	//************************************* Plotting corrected spectra MinBias PCM & PHOS Pi0 *********************************************************    cout <<"hier   !!!!!!!!!!!"<<endl; //**************************************************************************************************************************************************
	TCanvas* canvasSpectraDiffDetMinBias = new TCanvas("canvasSpectraDiffDetMinBias","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasSpectraDiffDetMinBias,  0.13, 0.01, 0.015, 0.08);
	
	canvasSpectraDiffDetMinBias->SetLogy();
	canvasSpectraDiffDetMinBias->SetLogx();
	histo2DPCMSpectraPi0AllCent->DrawCopy(); 
	
	DrawGammaSetMarker(histoPCMYieldPi0pPb, markerStylePCM,markerSizeInvYield, colorPCM , colorPCM);
	histoPCMYieldPi0pPb->Draw("p,same,e1");
	DrawGammaSetMarker(histoPHOSYieldPi0pPb, markerStylePHOS,markerSizeInvYield, colorPHOS , colorPHOS);
	histoPHOSYieldPi0pPb->Draw("p,same,e1");
	
	TLegend* legendSpectraDiffDetMinBias = new TLegend(0.16,0.11,0.5,0.2);
	legendSpectraDiffDetMinBias->SetFillColor(0);
	legendSpectraDiffDetMinBias->SetLineColor(0);
	legendSpectraDiffDetMinBias->SetTextSize(0.025);
	legendSpectraDiffDetMinBias->AddEntry(histoPCMYieldPi0pPb,Form("#pi^{0} PCM, %s - %s",collisionSystempPb.Data(),datePCM.Data()),"pf");
	legendSpectraDiffDetMinBias->AddEntry(histoPHOSYieldPi0pPb,Form("#pi^{0} PHOS, %s - %s",collisionSystempPb.Data(),datePHOS.Data()),"pf");
	legendSpectraDiffDetMinBias->Draw();
	
	canvasSpectraDiffDetMinBias->Update();
	canvasSpectraDiffDetMinBias->Print(Form("%s/Pi0_Spectra_MinBias_DiffDetectors.%s",outputDir.Data(),suffix.Data()));
  
	//**************************************************************************************************************************************************
	//************************************* Plotting ratio of corrected spectra MinBias PCM & PHOS Pi0 *************************************************
	//**************************************************************************************************************************************************
	
	TCanvas* canvasSpectraDiffDetMinBiasRatio = new TCanvas("canvasSpectraDiffDetMinBiasRatio","",200,10,1200,700);  // gives the page size
	DrawGammaCanvasSettings( canvasSpectraDiffDetMinBiasRatio,  0.08, 0.01, 0.015, 0.13);

	canvasSpectraDiffDetMinBiasRatio->SetLogx();
	TH2F * ratio2DYieldPi0pPbDiffDet;
	ratio2DYieldPi0pPbDiffDet = new TH2F("ratio2DYieldPi0pPbDiffDet","ratio2DYieldPi0pPbDiffDet",1000,0.23,30.,1000,0.25,2.15  );
	SetStyleHistoTH2ForGraphs(ratio2DYieldPi0pPbDiffDet, "p_{T} (GeV/c)","Detector A/Detector B", 0.05,0.064, 0.05,0.064, 0.8,0.6, 512, 505);
	ratio2DYieldPi0pPbDiffDet->DrawCopy(); 
	
	DrawGammaSetMarker(ratioYieldPi0pPb_PHOS_PCM, markerStylePCM,markerSizeInvYield, colorPCM , colorPCM);
	ratioYieldPi0pPb_PHOS_PCM->Draw("p,same,e1");
	
	TLatex *labelRatioPi0pPb = new TLatex(0.13,0.88,Form("#pi^{0},%s",collisionSystempPb.Data()));
	SetStyleTLatex( labelRatioPi0pPb, 0.07,4);
	labelRatioPi0pPb->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
		
	TLegend* legendSpectraDiffDetMinBiasRatio = new TLegend(0.75,0.90,0.89,0.95);
	legendSpectraDiffDetMinBiasRatio->SetFillColor(0);
	legendSpectraDiffDetMinBiasRatio->SetLineColor(0);
	legendSpectraDiffDetMinBiasRatio->SetTextSize(0.045);
	legendSpectraDiffDetMinBiasRatio->AddEntry(ratioYieldPi0pPb_PHOS_PCM,Form("PHOS/PCM"),"pf");
	legendSpectraDiffDetMinBiasRatio->Draw();

	TPaveText* pt= new TPaveText(0.75,0.82,0.89,0.9);
	pt->ConvertNDCtoPad();
	pt->AddText(Form("PCM spectra dated %s",datePCM.Data()));
	pt->AddText(Form("PHOS spectra dated %s",datePCM.Data()));
	pt->Draw();

	canvasSpectraDiffDetMinBiasRatio->Update();
	canvasSpectraDiffDetMinBiasRatio->Print(Form("%s/Pi0_Ratio_Spectra_MinBias_DiffDetectors.%s",outputDir.Data(),suffix.Data()));
   

}
	
