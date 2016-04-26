
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
#include "MultiplicityPlots.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

void MultiplicityPlots(){	

	
  date = ReturnDateString();
  gROOT->Reset();	
  gROOT->SetStyle("Plain");
	
  StyleSettingsThesis();	
  SetPlotStyle();

  TString fileNameConversionpPb = "data_PCMResults_pPb_DPG2014_EtaBinning.root";
  TString suffix = "eps";
  //___________________________________ Declaration of files _____________________________________________
  collisionSystem0020 = "0-20% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";		
  collisionSystem2040 = "20-40% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";      
  collisionSystem4060 = "40-60% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";    
  collisionSystem60100 = "60-100% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
  collisionSystem6080 = "60-80% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     

  collisionSystempPb = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   
  TString outputDir = Form("%s/MultiplicityPlots_DPG",suffix.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  nameFinalResDat = Form("%s/MultResultspPb_FitResults.dat",outputDir.Data());
 
  TString datePCM = "2014-03-04";


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

  cout << "60-80%" << endl;
  directoryPCMPi0pPb6080 = 				(TDirectory*)filePCMpPb->Get("Pi0_pPb_5.023TeV_60-80%"); 
  histoPCMNumberOfEvents6080= 			(TH1D*)filePCMpPb->Get("histoNumberOfEventspPb_5.023TeV60-80%");
  histoPCMYieldPi0pPb6080 = 				(TH1D*)directoryPCMPi0pPb6080->Get(nameHistoPCM.Data());
  graphPCMYieldPi0SysErrpPb6080= 	(TGraphAsymmErrors*)directoryPCMPi0pPb6080->Get(nameGraphPCM.Data());	
  graphPCMYieldPi0SysErrRAApPb6080= 	(TGraphAsymmErrors*)directoryPCMPi0pPb6080->Get("Pi0SystErrorA");	
  histoPCMMassPi0DatapPb6080 = 			(TH1D*)directoryPCMPi0pPb6080->Get("MassPi0");
  histoPCMMassPi0MCpPb6080 = 			(TH1D*)directoryPCMPi0pPb6080->Get("TrueMassPi0");
  histoPCMWidthPi0DatapPb6080 = 		(TH1D*)directoryPCMPi0pPb6080->Get("FWHMPi0MeV");
  histoPCMWidthPi0MCpPb6080 = 			(TH1D*)directoryPCMPi0pPb6080->Get("TrueFWHMPi0MeV");
  histoPCMMassPi0DatapPb6080->Scale(1000.);
  histoPCMMassPi0MCpPb6080->Scale(1000.);
  nPCMEventpPb6080 = 							histoPCMNumberOfEvents6080->GetBinContent(1);

   
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

  cout << "60-80%" << endl;
  directoryPCMEtapPb6080 = 				(TDirectory*)filePCMpPb->Get("Eta_pPb_5.023TeV_60-80%"); 
  histoPCMYieldEtapPb6080 = 				(TH1D*)directoryPCMEtapPb6080->Get(nameHistoPCMEta.Data());
  graphPCMYieldEtaSysErrpPb6080= 	(TGraphAsymmErrors*)directoryPCMEtapPb6080->Get(nameGraphPCM.Data());	
  graphPCMYieldEtaSysErrRAApPb6080= 	(TGraphAsymmErrors*)directoryPCMEtapPb6080->Get("EtaSystErrorA");	
  histoPCMEtaPi0RatiopPb6080 =            (TH1D*)directoryPCMEtapPb6080->Get("EtatoPi0Ratio");
  graphPCMEtaPi0RatioSysErrpPb6080=    (TGraphAsymmErrors*)directoryPCMEtapPb6080->Get("EtatoPi0RatioSys"); 
  histoPCMMassEtaDatapPb6080 = 			(TH1D*)directoryPCMEtapPb6080->Get("MassEta");
  histoPCMMassEtaMCpPb6080 = 			(TH1D*)directoryPCMEtapPb6080->Get("TrueMassEta");
  histoPCMWidthEtaDatapPb6080 = 		(TH1D*)directoryPCMEtapPb6080->Get("FWHMEtaMeV");
  histoPCMWidthEtaMCpPb6080 = 			(TH1D*)directoryPCMEtapPb6080->Get("TrueFWHMEtaMeV");
  histoPCMMassEtaDatapPb6080->Scale(1000.);
  histoPCMMassEtaMCpPb6080->Scale(1000.);


   
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

  //**************************************************************************************************************************************************
  //************************************* Plotting corrected minmum bias spectrum  PCM Pi0 ********************************************************
  //**************************************************************************************************************************************************	


  //---------read CGC calculation-----------------

	const char *fileNameEPS09sPi0CGC = "ExternalInputpPb/CGC_Pi0spectrum.dat";

	ifstream 		inCGC;	

	Int_t nlinesPi0CGC = 0;	
	inCGC.open(fileNameEPS09sPi0CGC,ios_base::in);
	
	Double_t xESPsPi0CGC[100], yESPsPi0CGC[100];
	
		
	while(!inCGC.eof()){
			nlinesPi0CGC++;
			inCGC >> xESPsPi0CGC[nlinesPi0CGC]  >> yESPsPi0CGC[nlinesPi0CGC];
			cout << nlinesPi0CGC << "         "  << xESPsPi0CGC[nlinesPi0CGC] << "         "  <<yESPsPi0CGC[nlinesPi0CGC]<<endl;
	
	}
	inCGC.close();
	
	TGraph* graphPi0CGC = new TGraph(nlinesPi0CGC,xESPsPi0CGC,yESPsPi0CGC);	

	//---------------------------------------
  
TCanvas* canvasSpectraPCMPi0 = new TCanvas("canvasSpectraPCM","",200,10,1200,1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSpectraPCMPi0,  0.13, 0.01, 0.015, 0.08);

  canvasSpectraPCMPi0->SetLogy();
  canvasSpectraPCMPi0->SetLogx();
  TH2F * histo2DPCMSpectraPi0;
  histo2DPCMSpectraPi0 = new TH2F("histo2DPCMSpectraPi0","histo2DPCMSpectraPi0",1000,0.23,20.,1000,1e-8,10000   );

   SetStyleHistoTH2ForGraphs(histo2DPCMSpectraPi0, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}	
  histo2DPCMSpectraPi0->GetXaxis()->SetLabelOffset(-0.01);
  histo2DPCMSpectraPi0->GetYaxis()->SetLabelOffset(0.01);
  histo2DPCMSpectraPi0->DrawCopy(); 

  DrawGammaSetMarkerTGraph(  graphPCMYieldPi0SysErrpPb, 20, 0.2, kBlack, kBlack);
  graphPCMYieldPi0SysErrpPb->SetFillColor(kGray+1);
  graphPCMYieldPi0SysErrpPb->Draw("same,2,p"); 
  DrawGammaSetMarker(histoPCMYieldPi0pPb, 34,markerSizeCommonSpectrum60100, colorPCM , colorPCM);
   histoPCMYieldPi0pPb->Draw("p,same,e1");
   histoPCMYieldPi0pPb->SetMarkerColor(kRed+2);
  histoPCMYieldPi0pPb->SetLineColor(kRed+2);
   graphPi0CGC->SetMarkerStyle(28);
   graphPi0CGC->SetMarkerSize(2);
   graphPi0CGC->SetMarkerColor(kGreen+2);
   graphPi0CGC->SetLineColor(kGreen+2);
  graphPi0CGC->SetLineWidth(2);
  graphPi0CGC->SetLineStyle(2);


 graphPi0CGC->Draw("p,same");

  TLatex *labelSpectraPi0LabelPCMpPbMB = new TLatex(0.65,0.85,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-} (PCM)");
  SetStyleTLatex( labelSpectraPi0LabelPCMpPbMB, 0.03,4);
  labelSpectraPi0LabelPCMpPbMB->Draw();
   //   TLatex *LabelpPb = new TLatex(0.65,0.9,"ALICE work in progress");
//     SetStyleTLatex( LabelpPb, 0.03,4);
//     LabelpPb->Draw();
 
    TLatex *labelSpectraPi0LabelpPbMB = new TLatex(0.65,0.8,"|y_{#pi^{0},lab}| < 0.8");
    SetStyleTLatex( labelSpectraPi0LabelpPbMB, 0.03,4);
   labelSpectraPi0LabelpPbMB->Draw();  
  
  TLegend* legendSpectraPi0PCMMB = new TLegend(0.16,0.25,0.5,0.3);
  legendSpectraPi0PCMMB->SetFillColor(0);
  legendSpectraPi0PCMMB->SetLineColor(0);
  legendSpectraPi0PCMMB->SetTextSize(0.025); 
  legendSpectraPi0PCMMB->AddEntry(histoPCMYieldPi0pPb,"minimum bias p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV","pf");
  legendSpectraPi0PCMMB->Draw();

  canvasSpectraPCMPi0->Update();
  canvasSpectraPCMPi0->Print(Form("%s/Pi0_MinimumBias_Spectrum_PCM.%s",outputDir.Data(),suffix.Data()));
   




  //**************************************************************************************************************************************************
  //************************************* Plotting corrected spectra all centralities PCM Pi0 ********************************************************
  //**************************************************************************************************************************************************
  TCanvas* canvasSpectraPCMPi0AllCent = new TCanvas("canvasSpectraPCMPi0AllCent","",200,10,1200,1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSpectraPCMPi0AllCent,  0.13, 0.01, 0.015, 0.08);

  canvasSpectraPCMPi0AllCent->SetLogy();
  canvasSpectraPCMPi0AllCent->SetLogx();
  TH2F * histo2DPCMSpectraPi0AllCent;
  histo2DPCMSpectraPi0AllCent = new TH2F("histo2DPCMSpectraPi0AllCent","histo2DPCMSpectraPi0AllCent",1000,0.23,20.,1000,1e-8,10000   );
  //  SetStyleHistoTH2ForGraphs(histo2DPCMSpectraPi0AllCent, "p_{T} (GeV/c)","corrected #pi^{0} yield (a. u.)", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 
   SetStyleHistoTH2ForGraphs(histo2DPCMSpectraPi0AllCent, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}
  histo2DPCMSpectraPi0AllCent->GetXaxis()->SetLabelOffset(-0.01);
  histo2DPCMSpectraPi0AllCent->GetYaxis()->SetLabelOffset(0.01);
  histo2DPCMSpectraPi0AllCent->DrawCopy(); 

   DrawGammaSetMarker(histoPCMYieldPi0pPb0020, markerStyleCommmonSpectrum0020,markerSizeCommonSpectrum0020, colorComb0020 , colorComb0020);
   histoPCMYieldPi0pPb0020->Scale(16);
   histoPCMYieldPi0pPb0020->Draw("p,same,e1");
   DrawGammaSetMarker(histoPCMYieldPi0pPb2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
   histoPCMYieldPi0pPb2040->Scale(8);
   histoPCMYieldPi0pPb2040->Draw("p,same,e1");
   DrawGammaSetMarker(histoPCMYieldPi0pPb4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
   histoPCMYieldPi0pPb4060->Scale(4);
   histoPCMYieldPi0pPb4060->Draw("p,same,e1");
   DrawGammaSetMarker(histoPCMYieldPi0pPb6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
    histoPCMYieldPi0pPb6080->Scale(2);
    histoPCMYieldPi0pPb6080->Draw("p,same,e1");
//    DrawGammaSetMarker(histoPCMYieldPi0pPb80100, markerStyleCommmonSpectrum80100,markerSizeCommonSpectrum80100, colorComb80100 , colorComb80100);
//    histoPCMYieldPi0pPb80100->Draw("p,same,e1");
   DrawGammaSetMarker(histoPCMYieldPi0pPb60100, markerStyleCommmonSpectrum60100,markerSizeCommonSpectrum60100, colorComb60100 , colorComb60100);
   histoPCMYieldPi0pPb60100->Scale(1);
   histoPCMYieldPi0pPb60100->Draw("p,same,e1"); 

TGraphAsymmErrors*  graphPCMYieldPi0SysErrpPbScale=ScaleGraph(graphPCMYieldPi0SysErrpPb,128);
	DrawGammaSetMarkerTGraph(  graphPCMYieldPi0SysErrpPbScale, 20, 0.2, kBlack, kBlack);
	  graphPCMYieldPi0SysErrpPbScale->SetFillColor(kGray+1);
	  graphPCMYieldPi0SysErrpPbScale->Draw("same,2,p"); 

  DrawGammaSetMarker(histoPCMYieldPi0pPb, 34,markerSizeCommonSpectrum60100, colorPCM , colorPCM);
   histoPCMYieldPi0pPb->Scale(128); 
   histoPCMYieldPi0pPb->Draw("p,same,e1"); 
//   DrawGammaSetMarker(histoPCMYieldPi0pPb0020, markerStyleCommmonSpectrum0020,markerSizeCommonSpectrum0020, colorComb0020 , colorComb0020);
//   histoPCMYieldPi0pPb0020->Draw("p,same,e1");
//   DrawGammaSetMarker(histoPCMYieldPi0pPb2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
//   histoPCMYieldPi0pPb2040->Draw("p,same,e1");
//   DrawGammaSetMarker(histoPCMYieldPi0pPb4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
//   histoPCMYieldPi0pPb4060->Draw("p,same,e1");
//   DrawGammaSetMarker(histoPCMYieldPi0pPb6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
//   histoPCMYieldPi0pPb6080->Draw("p,same,e1");

//   DrawGammaSetMarker(histoPCMYieldPi0pPb60100, markerStyleCommmonSpectrum60100,markerSizeCommonSpectrum60100, colorComb60100 , colorComb60100);
//   histoPCMYieldPi0pPb60100->Draw("p,same,e1"); 
//   DrawGammaSetMarker(histoPCMYieldPi0pPb, 34,markerSizeCommonSpectrum60100, colorPCM , colorPCM);
//   histoPCMYieldPi0pPb->Draw("p,same,e1");  

  TLatex *labelSpectraPi0LabelPCMpPb = new TLatex(0.65,0.85,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-} (PCM)");
  SetStyleTLatex( labelSpectraPi0LabelPCMpPb, 0.03,4);
  labelSpectraPi0LabelPCMpPb->Draw();
     TLatex *LabelpPb = new TLatex(0.65,0.9,"ALICE work in progress");
    SetStyleTLatex( LabelpPb, 0.03,4);
    LabelpPb->Draw();
 
    TLatex *labelSpectraPi0LabelpPb = new TLatex(0.65,0.8,"|y_{#pi^{0},lab}| < 0.8");
    SetStyleTLatex( labelSpectraPi0LabelpPb, 0.03,4);
   labelSpectraPi0LabelpPb->Draw();  
  
  TLegend* legendSpectraPi0PCM = new TLegend(0.16,0.11,0.5,0.3);
  legendSpectraPi0PCM->SetFillColor(0);
  legendSpectraPi0PCM->SetLineColor(0);
  legendSpectraPi0PCM->SetTextSize(0.025); 
  legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb,"minimum bias p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV(x128)","pf");
  legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb0020,"0-20% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV (x16)","pf");
  legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb2040,"20-40% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV (x8)","pf");
  legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb4060,"40-60% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV (x4)","pf");  
  legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb6080,"60-80% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV (x2)","pf");
  legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb60100,"60-100% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV (x1)","pf");
 //  legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb0020,collisionSystem0020.Data(),"pf");
//   legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb2040,collisionSystem2040.Data(),"pf");
//   legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb4060,collisionSystem4060.Data(),"pf");  
//   legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb6080,collisionSystem6080.Data(),"pf");
//   legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb60100,collisionSystem60100.Data(),"pf");
//   legendSpectraPi0PCM->AddEntry(histoPCMYieldPi0pPb,collisionSystempPb.Data(),"pf");

  legendSpectraPi0PCM->Draw();

  canvasSpectraPCMPi0AllCent->Update();
  canvasSpectraPCMPi0AllCent->Print(Form("%s/Pi0_Spectra_PCM.%s",outputDir.Data(),suffix.Data()));
   
  //**************************************************************************************************************************************************
  //************************************* Plotting corrected spectra all centralities PCM Eta ********************************************************
  //**************************************************************************************************************************************************
  TCanvas* canvasSpectraPCMEtaAllCent = new TCanvas("canvasSpectraPCMEtaAllCent","",200,10,1200,1100);  // gives the page size
  DrawGammaCanvasSettings( canvasSpectraPCMEtaAllCent,  0.13, 0.01, 0.015, 0.08);

  canvasSpectraPCMEtaAllCent->SetLogy();
  canvasSpectraPCMEtaAllCent->SetLogx();
  TH2F * histo2DPCMSpectraEtaAllCent;
  histo2DPCMSpectraEtaAllCent = new TH2F("histo2DPCMSpectraEtaAllCent","histo2DPCMSpectraEtaAllCent",1000,0.4,12.,1000,1e-7,100   );
  SetStyleHistoTH2ForGraphs(histo2DPCMSpectraEtaAllCent, "p_{T} (GeV/c)","corrected #eta yield (a. u.)", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}
   SetStyleHistoTH2ForGraphs(histo2DPCMSpectraEtaAllCent, "p_{T} (GeV/c)","#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV/c)^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{dN^{AA}}{dp_{T} dy}}{ #frac{1}{N_{evt}^{pp}}#frac{dN^{pp}}{dp_{T} dy}}
  histo2DPCMSpectraEtaAllCent->GetXaxis()->SetLabelOffset(-0.01);
  histo2DPCMSpectraEtaAllCent->GetYaxis()->SetLabelOffset(0.01);
  histo2DPCMSpectraEtaAllCent->DrawCopy(); 
  DrawGammaSetMarker(histoPCMYieldEtapPb0020, markerStyleCommmonSpectrum0020,markerSizeCommonSpectrum0020, colorComb0020 , colorComb0020); 
  histoPCMYieldEtapPb0020->Scale(16);
  histoPCMYieldEtapPb0020->Draw("p,same,e1");  
  DrawGammaSetMarker(histoPCMYieldEtapPb2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040); 
  histoPCMYieldEtapPb2040->Scale(8);
  histoPCMYieldEtapPb2040->Draw("p,same,e1");  
  DrawGammaSetMarker(histoPCMYieldEtapPb4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
  histoPCMYieldEtapPb4060->Scale(4);
  histoPCMYieldEtapPb4060->Draw("p,same,e1");  
  DrawGammaSetMarker(histoPCMYieldEtapPb6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080); 
  histoPCMYieldEtapPb6080->Scale(2);
  histoPCMYieldEtapPb6080->Draw("p,same,e1"); 

  DrawGammaSetMarker(histoPCMYieldEtapPb60100, markerStyleCommmonSpectrum60100,markerSizeCommonSpectrum60100, colorComb60100 , colorComb60100); 

  histoPCMYieldEtapPb60100->Draw("p,same,e1");  

 TGraphAsymmErrors*  graphPCMYieldEtaSysErrpPbScale=ScaleGraph(graphPCMYieldEtaSysErrpPb,128);
	DrawGammaSetMarkerTGraph(  graphPCMYieldEtaSysErrpPbScale, 20, 0.2, kBlack, kBlack);
	  graphPCMYieldEtaSysErrpPbScale->SetFillColor(kGray+1);
	  graphPCMYieldEtaSysErrpPbScale->Draw("same,2,p");   
  DrawGammaSetMarker(histoPCMYieldEtapPb, 34,markerSizeCommonSpectrum60100, colorPCM , colorPCM);
  histoPCMYieldEtapPb->Scale(128);
  histoPCMYieldEtapPb->Draw("p,same,e1"); 
  TLatex *labelSpectraEtaLabelPCMpPb = new TLatex(0.65,0.85,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-} (PCM)");
  SetStyleTLatex( labelSpectraEtaLabelPCMpPb, 0.03,4);
  labelSpectraEtaLabelPCMpPb->Draw();
     TLatex *LabelpPbEta = new TLatex(0.65,0.9,"ALICE work in progress");
    SetStyleTLatex( LabelpPbEta, 0.03,4);
    LabelpPbEta->Draw();
 
    TLatex *labelSpectraPi0LabelpPbEta = new TLatex(0.65,0.8,"|y_{#eta,lab}| < 0.8");
    SetStyleTLatex( labelSpectraPi0LabelpPbEta, 0.03,4);
   labelSpectraPi0LabelpPbEta->Draw();  
  TLegend* legendSpectraEtaPCM = new TLegend(0.16,0.11,0.5,0.3);
  legendSpectraEtaPCM->SetFillColor(0);
  legendSpectraEtaPCM->SetLineColor(0);
  legendSpectraEtaPCM->SetTextSize(0.025); 
  legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb,"minimum bias p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV(x128)","pf");
  legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb0020,"0-20% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV (x16)","pf");
  legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb2040,"20-40% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV (x8)","pf");
  legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb4060,"40-60% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV (x4)","pf");  
  legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb6080,"60-80% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV (x2)","pf");
  legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb60100,"60-100% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV (x1)","pf");
//   legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb0020,collisionSystem0020.Data(),"pf");
//   legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb2040,collisionSystem2040.Data(),"pf");
//   legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb4060,collisionSystem4060.Data(),"pf");  
//   legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb6080,collisionSystem6080.Data(),"pf");
//   legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb60100,collisionSystem60100.Data(),"pf");
//   legendSpectraEtaPCM->AddEntry(histoPCMYieldEtapPb,collisionSystempPb.Data(),"pf");

  legendSpectraEtaPCM->Draw();

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

  DrawGammaSetMarker(histoPCMEtaPi0RatiopPb0020, markerStyleCommmonSpectrum0020,markerSizeCommonSpectrum0020, colorComb0020 , colorComb0020);
  histoPCMEtaPi0RatiopPb0020->Draw("p,same,e1");
  DrawGammaSetMarker(histoPCMEtaPi0RatiopPb2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
  histoPCMEtaPi0RatiopPb2040->Draw("p,same,e1");
  DrawGammaSetMarker(histoPCMEtaPi0RatiopPb4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
  histoPCMEtaPi0RatiopPb4060->Draw("p,same,e1");
  DrawGammaSetMarker(histoPCMEtaPi0RatiopPb6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
  histoPCMEtaPi0RatiopPb6080->Draw("p,same,e1");
  DrawGammaSetMarker(histoPCMEtaPi0RatiopPb60100, markerStyleCommmonSpectrum60100,markerSizeCommonSpectrum60100, colorComb60100 , colorComb60100);
  histoPCMEtaPi0RatiopPb60100->Draw("p,same,e1");
   
   
  DrawGammaSetMarker(histoPCMEtaPi0RatiopPb, 34,markerSizeCommonSpectrum60100, colorPCM , colorPCM);
  histoPCMEtaPi0RatiopPb->Draw("p,same,e1"); 
 
      TLatex *LabelpPbEtaPi0 = new TLatex(0.65,0.9,"ALICE work in progress");
    SetStyleTLatex( LabelpPbEtaPi0, 0.03,4);
    LabelpPbEtaPi0->Draw();
 
    TLatex *labelSpectraPi0LabelpPbEtaPi0 = new TLatex(0.65,0.8,"|y_{meson,lab}| < 0.8");
    SetStyleTLatex( labelSpectraPi0LabelpPbEtaPi0, 0.03,4);
   labelSpectraPi0LabelpPbEtaPi0->Draw();  
   
  TLegend* legendEtaToPi0PCM = new TLegend(0.16,0.76,0.5,0.95);
  legendEtaToPi0PCM->SetFillColor(0);
  legendEtaToPi0PCM->SetLineColor(0);
  legendEtaToPi0PCM->SetTextSize(0.025);
  legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb0020,collisionSystem0020.Data(),"pf");
  legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb2040,collisionSystem2040.Data(),"pf");
  legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb4060,collisionSystem4060.Data(),"pf"); 
  legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb6080,collisionSystem6080.Data(),"pf");
  legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb60100,collisionSystem60100.Data(),"pf");
  legendEtaToPi0PCM->AddEntry(histoPCMEtaPi0RatiopPb,collisionSystempPb.Data(),"pf");

  legendEtaToPi0PCM->Draw();
   
  canvasEtaToPi0PCMAllCent->Update();
  canvasEtaToPi0PCMAllCent->Print(Form("%s/EtaToPi0_Ratio_PCM.%s",outputDir.Data(),suffix.Data()));

  //************************************************************************************************************************************************  	   **
  //************************************* Plotting eta/pi0 ratioMin Bias ********************************************************
  //**************************************************************************************************************************************************


  TCanvas* canvasEtaToPi0PCMMB = new TCanvas("canvasEtaToPi0PCMMB","",200,10,1200,1100);  // gives the page size
  DrawGammaCanvasSettings( canvasEtaToPi0PCMMB,  0.08, 0.01, 0.015, 0.08);
   
  //    canvasEtaToPi0PCMAllCent->SetLogy();
  //    canvasEtaToPi0PCMAllCent->SetLogx();
  TH2F * histo2DEtaToPi0PCMMB;
  histo2DEtaToPi0PCMMB = new TH2F("histo2DEtaToPi0PCMMB","histo2DEtaToPi0PCMMB",1000,0.,8.,1000,0,1.1   );
  SetStyleHistoTH2ForGraphs(histo2DEtaToPi0PCMMB, "p_{T} (GeV/c)","#eta/#pi^{0}", 0.03,0.04, 0.03,0.04, 0.83,0.95);
  //    histo2DEtaToPi0PCMMB->GetXaxis()->SetLabelOffset(-0.01);
  //    histo2DEtaToPi0PCMMB->GetYaxis()->SetLabelOffset(0.01);
  histo2DEtaToPi0PCMMB->DrawCopy(); 

    
  DrawGammaSetMarker(histoPCMEtaPi0RatiopPb, 34,markerSizeCommonSpectrum60100, colorPCM , colorPCM);
  histoPCMEtaPi0RatiopPb->Draw("p,same,e1");  
  
      TLatex *LabelpPbEtaPi0MB = new TLatex(0.65,0.9,"ALICE work in progress");
    SetStyleTLatex( LabelpPbEtaPi0MB, 0.03,4);
    LabelpPbEtaPi0MB->Draw();
 
    TLatex *labelSpectraPi0LabelpPbEtaPi0MB = new TLatex(0.65,0.8,"|y_{meson,lab}| < 0.8");
    SetStyleTLatex( labelSpectraPi0LabelpPbEtaPi0MB, 0.03,4);
   labelSpectraPi0LabelpPbEtaPi0MB->Draw(); 
   
  TLegend* legendEtaToPi0PCMMB = new TLegend(0.16,0.76,0.5,0.95);
  legendEtaToPi0PCMMB->SetFillColor(0);
  legendEtaToPi0PCMMB->SetLineColor(0);
  legendEtaToPi0PCMMB->SetTextSize(0.025);
  legendEtaToPi0PCMMB->AddEntry(histoPCMEtaPi0RatiopPb,"minimum bias p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV","pf");

  legendEtaToPi0PCMMB->Draw();
   
  canvasEtaToPi0PCMMB->Update();
  canvasEtaToPi0PCMMB->Print(Form("%s/EtaToPi0_Ratio_PCM_MB.%s",outputDir.Data(),suffix.Data()));
}
