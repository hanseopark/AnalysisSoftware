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
#include "TaskV1/ProduceFinalResults.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"
#include "CombineMesonMeasurements.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;	

struct SysErrorConversion {
	Double_t value;
	Double_t error;
	//	TString name;
};

void CombineMesonMeasurements(TString fileNameConversions = "", const char *fileNameConversionsPrelim = "", TString suffix = "eps", TString isMC= "", TString thesisPlots = "", TString bWCorrection="X"){	

	date = ReturnDateString();
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TString dateForOutput 			= ReturnDateStringForOutput();
	cout << dateForOutput.Data() << endl;
	//___________________________________ Declaration of files _____________________________________________
	collisionSystem7TeV 			= "pp, #sqrt{#it{s}} = 7 TeV";		
	collisionSystem2760GeV 			= "pp, #sqrt{#it{s}} = 2.76 TeV";		
	collisionSystem900GeV 			= "pp, #sqrt{#it{s}} = 0.9 TeV";
	collisionSystemCombined 		= "pp, #sqrt{#it{s}} = 0.9 & 7 TeV";
	collisionSystemCombinedReallyAll = "pp #sqrt{#it{s}} = 0.9, 2.76, 7 TeV";
		
	fileNameCaloPhos7TeV 			= "ExternalInput/PHOS/7TeV/PHOS_pi0_7TeV_20111030_BWcorr.root";
	fileNameCaloPhosOmega7TeV 		= "ExternalInput/PHOS/7TeV/PHOS_pp_omega_7TeV_07082012.root";
	fileNameCaloEmcal7TeV 			= "ExternalInput/EMCAL/7TeV/pi0EMCAL_12052011.root";
	fileNameCaloEmcalEtaToPi07TeV 	= "ExternalInput/EMCAL/7TeV/eta2pi0EMCAL_11052011.root";
	fileNameCaloPhos7TeVEta 		= "ExternalInput/PHOS/7TeV/PHOS_eta_7TeV_06052011.root";
	fileNameChargedSpectra7TeV 		= "ExternalInput/IdentifiedCharged/1PtdNdPtSpectra_INEL_LHC10d_7TeV_Final_301110.txt";
	fileNameChargedSpectra2760GeV 	= "ExternalInput/IdentifiedCharged/charged_dNdPt_pp_276.root";
	fileNameChargedExpectation7TeV	= "ExternalInput/IdentifiedCharged/ChargedPythia22092010.dat";
	fileNamePHOSMassData7TeV 		= "ExternalInput/PHOS/7TeV/PHOS_Pi0MassPosData-CommonBinning-20111109.root";
	fileNamePHOSMassMC7TeV 			= "ExternalInput/PHOS/7TeV/PHOS_Pi0MassPosMC-CommonBinning-20111108.root";
	fileNameEMCALMassData7TeVPi0 	= "ExternalInput/EMCAL/7TeV/EMCALResults_pp_v2Clusterizer.root";
	fileNameEMCALMassData7TeVEta 	= "ExternalInput/EMCAL/7TeV/EMCAL_DataEtaMassWidth.root";
	fileNameEMCALMassMC7TeVPi0 		= "ExternalInput/EMCAL/7TeV/EMCAL_SimPi0MassWidth.root";
	fileNameEMCALMassMC7TeVEta 		= "ExternalInput/EMCAL/7TeV/EMCAL_SimEtaMassWidth.root";
	fileNamePHOSEtaToPi0 			= "ExternalInput/PHOS/7TeV/PHOS_etaPi_ratio_7TeV_20111217.root";
	fileNameCaloPhos900GeV 			= "ExternalInput/PHOS/0.9GeV/PHOS_pp_pi0_900GeV_20110413_K0Scorr.root";
	if(bWCorrection.CompareTo("X")==0){
		fileNameCaloPhos900GeV 		= "ExternalInput/PHOS/0.9GeV/PHOS_pp_pi0_900GeV_noBWcorr_K0Scorr_20111206.root";
		fileNameCaloPhos7TeV 		= "ExternalInput/PHOS/7TeV/PHOS_pp_pi0_7TeV_20120515.root";
		fileNameCaloPhos7TeVEta 	= "ExternalInput/PHOS/7TeV/PHOS_pp_eta_7TeV_20120515.root";
	}
	fileNameCaloPhos2760GeV			= "ExternalInput/PHOS/2.76TeV/LHC11a_PHOS_pi0_pp2760_noBWCorr_FDcorr_20140218.root";
   cout << "USING PHOS FILE 2.76TeV: " << fileNameCaloPhos2760GeV<< endl;
   
	TString outputDir 				= Form("%s/%s/CombineMesonMeasurements%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
	nameFinalResDat 				= Form("%s/CombinedResults%s_FitResults.dat",dateForOutput.Data(),bWCorrection.Data());
	cout << outputDir.Data() << endl;
   cout << fileNameConversions.Data() << endl;
//    TString bla = Form("cp %s %s/InputGammaConv.root", fileNameConversions.Data(), outputDir.Data());
//    cout << bla.Data() << endl;
	gSystem->Exec("mkdir -p "+outputDir);
 	gSystem->Exec(Form("cp %s %s/InputGammaConv.root", fileNameConversions.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPHOS900GeV.root", fileNameCaloPhos900GeV, outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPHOS7TeV.root", fileNameCaloPhos7TeV, outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPHOS7TeVEta.root", fileNameCaloPhos7TeVEta, outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPHOS2760GeV.root", fileNameCaloPhos2760GeV, outputDir.Data()));
	if(thesisPlots.CompareTo("thesis") == 0){// means we want to plot values for the pi0
      thesis = kTRUE;
   }
		
	//declaration for printing logo 
	
	if (isMC.CompareTo("kTRUE")==0){ 
		prefix2 = "MC";
		pictDrawingOptions[1] = kTRUE;
	} else {	
		prefix2 = "data";
		pictDrawingOptions[1] = kFALSE;
	}
	
	mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();
	
	//************************** Read data for conversions **************************************************
	// File definitions
	fileConversions = 		new TFile(fileNameConversions.Data());
	fileConversionsPrelim = 		new TFile(fileNameConversionsPrelim);
	directoryPi07TeV = 			(TDirectory*)fileConversions->Get("Pi07TeV"); 
	directoryEta7TeV = 			(TDirectory*)fileConversions->Get("Eta7TeV");
	histoNumberOfEvents= 					(TH1D*)fileConversions->Get("histoNumberOfEvents7TeV");
	histoMassMesonPi0 = 					(TH1D*)directoryPi07TeV->Get("MassPi0");
	histoFWHMMesonPi0MeV = 					(TH1D*)directoryPi07TeV->Get("FWHMPi0MeV");
	histoTrueMassMesonPi0 = 					(TH1D*)directoryPi07TeV->Get("TrueMassPi0");
	histoTrueFWHMMesonPi0MeV = 				(TH1D*)directoryPi07TeV->Get("TrueFWHMPi0MeV");
	histoAccPi0 = 							(TH1D*)directoryPi07TeV->Get("AcceptancePi0");
	histoTrueEffPtPi0 = 					(TH1D*)directoryPi07TeV->Get("EfficiencyPi0");
	histoPi0CorrYieldBinShifted = 			(TH1D*)directoryPi07TeV->Get("CorrectedYieldPi0BinShifted");
	histoPi0CorrYieldMtBinShifted = 			(TH1D*)directoryPi07TeV->Get("CorrectedYieldPi0MtBinShifted");
	histoPi0CorrYieldXtBinShifted = 			(TH1D*)directoryPi07TeV->Get("CorrectedYieldPi0XtBinShifted");
	histoInvCrossSectionPi0= 				(TH1D*)directoryPi07TeV->Get("InvCrossSectionPi0");
	graphInvCrossSectionSysPi0= 				(TGraphAsymmErrors*)directoryPi07TeV->Get("InvCrossSectionPi0Sys");
	graphInvCrossSectionSysAPi0= 				(TGraphAsymmErrors*)directoryPi07TeV->Get("InvCrossSectionPi0SysA");
	graphCorrectedYieldPi0SysErr= 			(TGraphAsymmErrors*)directoryPi07TeV->Get("Pi0SystError");
	graphCorrectedYieldPi0SysErrBinShifted= 			(TGraphAsymmErrors*)directoryPi07TeV->Get("Pi0SystErrorBinShifted");
	histoEtaToPi0Phojet7TeV= 				(TH1D*)directoryPi07TeV->Get("EtaToPi0MCTruthPhojet");
	histoEtaToPi0Pythia7TeV= 				(TH1D*)directoryPi07TeV->Get("EtaToPi0MCTruthPythia");
	histoPi0ToChargedPhojet7TeV= 				(TH1D*)directoryPi07TeV->Get("Pi0ToChargedPhojet");
	histoPi0ToChargedPythia7TeV= 				(TH1D*)directoryPi07TeV->Get("Pi0ToChargedPythia");
	cout << "bis hier her geschafft" << endl;
	histoCorrectedYieldEta= 					(TH1D*)directoryEta7TeV->Get("CorrectedYieldEta");
	histoUnCorrectedYieldEta= 				(TH1D*)directoryEta7TeV->Get("RAWYieldPerEventsEta");
	histoMassMesonEta=						(TH1D*)directoryEta7TeV->Get("MassEta");
	histoFWHMMesonEtaMeV = 					(TH1D*)directoryEta7TeV->Get("FWHMEtaMeV");
	histoTrueMassMesonEta = 					(TH1D*)directoryEta7TeV->Get("TrueMassEta");
	histoTrueFWHMMesonEtaMeV = 				(TH1D*)directoryEta7TeV->Get("TrueFWHMEtaMeV");
	histoEtaCorrYieldBinShifted= 				(TH1D*)directoryEta7TeV->Get("CorrectedYieldEtaBinShifted");
	histoEtaCorrYieldMtBinShifted= 				(TH1D*)directoryEta7TeV->Get("CorrectedYieldEtaMtBinShifted");
	histoEtaCorrYieldXtBinShifted= 				(TH1D*)directoryEta7TeV->Get("CorrectedYieldEtaXtBinShifted");
	histoInvCrossSectionEta= 				(TH1D*)directoryEta7TeV->Get("InvCrossSectionEta");
	graphInvCrossSectionSysAEta= 				(TGraphAsymmErrors*)directoryEta7TeV->Get("InvCrossSectionEtaSysA");
	graphInvCrossSectionSysEta= 				(TGraphAsymmErrors*)directoryEta7TeV->Get("InvCrossSectionEtaSys");
	graphCorrectedYieldEtaSysErr= 			(TGraphAsymmErrors*)directoryEta7TeV->Get("EtaSystError");
	graphCorrectedYieldEtaSysErrBinShifted= 			(TGraphAsymmErrors*)directoryEta7TeV->Get("EtaSystErrorBinShifted");
	graphCorrectedYieldEtaSysErrMtBinShifted= 			(TGraphAsymmErrors*)directoryEta7TeV->Get("EtaSystErrorBinShiftedMt");
	histoRatioEtaPi0= 						(TH1D*)directoryEta7TeV->Get("EtatoPi0RatioConversionBinShifted");	
	graphSystErrRatio= 						(TGraphAsymmErrors*)directoryEta7TeV->Get("EtatoPi0RatioConversionBinShiftedSys");
	nEvt7TeV = 					histoNumberOfEvents->GetBinContent(1);
	Bool_t kV0ANDConv7TeV = kFALSE;
	if (histoNumberOfEvents->GetNbinsX() > 1 && histoNumberOfEvents->GetBinContent(2) > 0)  kV0ANDConv7TeV = kTRUE;
	cout << "7 TeV data read, kVOAND: " <<kV0ANDConv7TeV  <<endl;
	
	directoryPi02760GeV = 			(TDirectory*)fileConversions->Get("Pi02.76TeV"); 
	if (conference) {
// 		directoryPi02760GeV = 			(TDirectory*)fileConversionsPrelim->Get("Pi02.76TeV"); 
		directoryEta2760GeV = 			(TDirectory*)fileConversionsPrelim->Get("Eta2.76TeV"); 
		cout << "reading preliminary file for 2.76TeV" << endl;
	} else {
		directoryEta2760GeV = 			(TDirectory*)fileConversions->Get("Eta2.76TeV"); 
	}

	histoNumberOfEvents2760GeV= 				(TH1D*)fileConversions->Get("histoNumberOfEvents");
	histoMassMesonPi02760GeV = 				(TH1D*)directoryPi02760GeV->Get("MassPi0");
	histoFWHMMesonPi0MeV2760GeV = 				(TH1D*)directoryPi02760GeV->Get("FWHMPi0MeV");
	histoTrueMassMesonPi02760GeV = 			(TH1D*)directoryPi02760GeV->Get("TrueMassPi0");
	histoTrueFWHMMesonPi0MeV2760GeV = 			(TH1D*)directoryPi02760GeV->Get("TrueFWHMPi0MeV");
	histoAccPi02760GeV = 					(TH1D*)directoryPi02760GeV->Get("AcceptancePi0");
	histoTrueEffPtPi02760GeV = 				(TH1D*)directoryPi02760GeV->Get("EfficiencyPi0");
	histoPi0CorrYieldBinShifted2760GeV = 		(TH1D*)directoryPi02760GeV->Get("CorrectedYieldPi0BinShifted");
	histoPi0CorrYieldMtBinShifted2760GeV = 		(TH1D*)directoryPi02760GeV->Get("CorrectedYieldPi0MtBinShifted");
	histoInvCrossSectionPi02760GeV= 			(TH1D*)directoryPi02760GeV->Get("InvCrossSectionPi0");
	graphInvCrossSectionSysAPi02760GeV= 		(TGraphAsymmErrors*)directoryPi02760GeV->Get("InvCrossSectionPi0SysA");
	graphInvCrossSectionSysPi02760GeV= 			(TGraphAsymmErrors*)directoryPi02760GeV->Get("InvCrossSectionPi0Sys");
	graphCorrectedYieldPi0SysErr2760GeV= 		(TGraphAsymmErrors*)directoryPi02760GeV->Get("Pi0SystError");
	graphCorrectedYieldPi0SysErrBinShifted2760GeV= 		(TGraphAsymmErrors*)directoryPi02760GeV->Get("Pi0SystErrorBinShifted");
	histoEtaToPi0Phojet2760GeV= 				(TH1D*)directoryPi02760GeV->Get("EtaToPi0MCTruthPhojet");
	histoEtaToPi0Pythia2760GeV= 				(TH1D*)directoryPi02760GeV->Get("EtaToPi0MCTruthPythia");
	histoPi0ToChargedPhojet2760GeV= 			(TH1D*)directoryPi02760GeV->Get("Pi0ToChargedPhojet");
	histoPi0ToChargedPythia2760GeV= 			(TH1D*)directoryPi02760GeV->Get("Pi0ToChargedPythia");	
	histoCorrectedYieldEta2760GeV= 			(TH1D*)directoryEta2760GeV->Get("CorrectedYieldEta");
	histoUnCorrectedYieldEta2760GeV= 			(TH1D*)directoryEta2760GeV->Get("RAWYieldPerEventsEta");
	histoMassMesonEta2760GeV=				(TH1D*)directoryEta2760GeV->Get("MassEta");
	histoFWHMMesonEtaMeV2760GeV = 			(TH1D*)directoryEta2760GeV->Get("FWHMEtaMeV");
	histoTrueMassMesonEta2760GeV = 			(TH1D*)directoryEta2760GeV->Get("TrueMassEta");
	histoTrueFWHMMesonEtaMeV2760GeV = 			(TH1D*)directoryEta2760GeV->Get("TrueFWHMEtaMeV");
	histoEtaCorrYieldBinShifted2760GeV= 		(TH1D*)directoryEta2760GeV->Get("CorrectedYieldEtaBinShifted");
	histoEtaCorrYieldMtBinShifted2760GeV= 		(TH1D*)directoryEta2760GeV->Get("CorrectedYieldEtaMtBinShifted");
	histoInvCrossSectionEta2760GeV= 			(TH1D*)directoryEta2760GeV->Get("InvCrossSectionEta");
	graphInvCrossSectionSysAEta2760GeV= 		(TGraphAsymmErrors*)directoryEta2760GeV->Get("InvCrossSectionEtaSysA");
	graphInvCrossSectionSysEta2760GeV= 		(TGraphAsymmErrors*)directoryEta2760GeV->Get("InvCrossSectionEtaSys");
	graphCorrectedYieldEtaSysErr2760GeV= 		(TGraphAsymmErrors*)directoryEta2760GeV->Get("EtaSystError");
	graphCorrectedYieldEtaSysErrBinShifted2760GeV= 		(TGraphAsymmErrors*)directoryEta2760GeV->Get("EtaSystErrorBinShifted");
	graphCorrectedYieldEtaSysErrMtBinShifted2760GeV= 		(TGraphAsymmErrors*)directoryEta2760GeV->Get("EtaSystErrorBinShiftedMt");
	histoRatioEtaPi02760GeV= 				(TH1D*)directoryEta2760GeV->Get("EtatoPi0RatioConversionBinShifted");	
	graphSystErrRatio2760GeV= 				(TGraphAsymmErrors*)directoryEta2760GeV->Get("EtatoPi0RatioConversionBinShiftedSys");
	
	
	nEvt2760GeV = 					histoNumberOfEvents2760GeV->GetBinContent(1);
	Bool_t kV0ANDConv2760GeV = kFALSE;
	if (histoNumberOfEvents2760GeV->GetNbinsX() > 1 && histoNumberOfEvents2760GeV->GetBinContent(2) > 0)  kV0ANDConv2760GeV = kTRUE;
	cout << "2.76 TeV eingelesen, kVOAND: " <<kV0ANDConv2760GeV  <<endl;
	
	directoryPi0900GeV = 			(TDirectory*)fileConversions->Get("Pi0900GeV"); 
	if (conference) {
		directoryEta900GeV = 			(TDirectory*)fileConversionsPrelim->Get("Eta900GeV"); 
		cout << "reading preliminary file for eta,  900GeV" << endl;
	} else {
		directoryEta900GeV = 			(TDirectory*)fileConversions->Get("Eta900GeV"); 
	}	
	
	histoNumberOfEvents900GeV= 				(TH1D*)fileConversions->Get("histoNumberOfEvents");
	histoMassMesonPi0900GeV = 				(TH1D*)directoryPi0900GeV->Get("MassPi0");
	histoFWHMMesonPi0MeV900GeV = 				(TH1D*)directoryPi0900GeV->Get("FWHMPi0MeV");
	histoTrueMassMesonPi0900GeV = 			(TH1D*)directoryPi0900GeV->Get("TrueMassPi0");
	histoTrueFWHMMesonPi0MeV900GeV = 			(TH1D*)directoryPi0900GeV->Get("TrueFWHMPi0MeV");
	histoAccPi0900GeV = 					(TH1D*)directoryPi0900GeV->Get("AcceptancePi0");
	histoTrueEffPtPi0900GeV = 				(TH1D*)directoryPi0900GeV->Get("EfficiencyPi0");
	histoPi0CorrYieldBinShifted900GeV = 		(TH1D*)directoryPi0900GeV->Get("CorrectedYieldPi0BinShifted");
	histoInvCrossSectionPi0900GeV= 			(TH1D*)directoryPi0900GeV->Get("InvCrossSectionPi0");
	graphInvCrossSectionSysAPi0900GeV= 		(TGraphAsymmErrors*)directoryPi0900GeV->Get("InvCrossSectionPi0SysA");
	graphInvCrossSectionSysPi0900GeV= 			(TGraphAsymmErrors*)directoryPi0900GeV->Get("InvCrossSectionPi0Sys");
	graphCorrectedYieldPi0SysErr900GeV= 		(TGraphAsymmErrors*)directoryPi0900GeV->Get("Pi0SystError");
	graphCorrectedYieldPi0SysErrBinShifted900GeV= 		(TGraphAsymmErrors*)directoryPi0900GeV->Get("Pi0SystErrorBinShifted");
	histoEtaToPi0Phojet900GeV= 				(TH1D*)directoryPi0900GeV->Get("EtaToPi0MCTruthPhojet");
	histoEtaToPi0Pythia900GeV= 				(TH1D*)directoryPi0900GeV->Get("EtaToPi0MCTruthPythia");
	histoPi0ToChargedPhojet900GeV= 			(TH1D*)directoryPi0900GeV->Get("Pi0ToChargedPhojet");
	histoPi0ToChargedPythia900GeV= 			(TH1D*)directoryPi0900GeV->Get("Pi0ToChargedPythia");

	histoCorrectedYieldEta900GeV= 				(TH1D*)directoryEta900GeV->Get("CorrectedYieldEta");
	histoUnCorrectedYieldEta900GeV= 				(TH1D*)directoryEta900GeV->Get("RAWYieldPerEventsEta");
	histoMassMesonEta900GeV=						(TH1D*)directoryEta900GeV->Get("MassEta");
	histoFWHMMesonEtaMeV900GeV = 					(TH1D*)directoryEta900GeV->Get("FWHMEtaMeV");
	histoTrueMassMesonEta900GeV = 				(TH1D*)directoryEta900GeV->Get("TrueMassEta");
	histoTrueFWHMMesonEtaMeV900GeV = 				(TH1D*)directoryEta900GeV->Get("TrueFWHMEtaMeV");
	histoEtaCorrYieldBinShifted900GeV= 			(TH1D*)directoryEta900GeV->Get("CorrectedYieldEtaBinShifted");
	histoInvCrossSectionEta900GeV= 				(TH1D*)directoryEta900GeV->Get("InvCrossSectionEta");
	graphInvCrossSectionSysAEta900GeV= 			(TGraphAsymmErrors*)directoryEta900GeV->Get("InvCrossSectionEtaSysA");
	graphInvCrossSectionSysEta900GeV= 				(TGraphAsymmErrors*)directoryEta900GeV->Get("InvCrossSectionEtaSys");
	graphCorrectedYieldEtaSysErr900GeV= 			(TGraphAsymmErrors*)directoryEta900GeV->Get("EtaSystError");
	graphCorrectedYieldEtaSysErrBinShifted900GeV= 			(TGraphAsymmErrors*)directoryEta900GeV->Get("EtaSystErrorBinShifted");
	histoRatioEtaPi0900GeV= 						(TH1D*)directoryEta900GeV->Get("EtatoPi0RatioConversion");	
	graphSystErrRatio900GeV= 					(TGraphAsymmErrors*)directoryEta900GeV->Get("EtatoPi0RatioConversionSys");
	
	nEvt900GeV = 					histoNumberOfEvents900GeV->GetBinContent(1);
	Bool_t kV0ANDConv900GeV = kFALSE;
	if (histoNumberOfEvents900GeV->GetNbinsX() > 1 && histoNumberOfEvents900GeV->GetBinContent(2) > 0)  kV0ANDConv900GeV = kTRUE;
	cout << "900 GeV eingelesen, kVOAND: " <<kV0ANDConv2760GeV  <<endl;

	
	histoMassMesonPi0MinusExp = CalculateMassMinusExpectedMass(histoMassMesonPi0,mesonMassExpectPi0+0.005);
	histoTrueMassMesonPi0MinusExp = CalculateMassMinusExpectedMass(histoTrueMassMesonPi0,mesonMassExpectPi0+0.005);
	histoMassMesonPi0MinusExp900GeV = CalculateMassMinusExpectedMass(histoMassMesonPi0900GeV,mesonMassExpectPi0+0.005);
	histoTrueMassMesonPi0MinusExp900GeV = CalculateMassMinusExpectedMass(histoTrueMassMesonPi0900GeV,mesonMassExpectPi0+0.005);
	histoMassMesonPi0MinusExp2760GeV = CalculateMassMinusExpectedMass(histoMassMesonPi02760GeV,mesonMassExpectPi0+0.005);
	histoTrueMassMesonPi0MinusExp2760GeV = CalculateMassMinusExpectedMass(histoTrueMassMesonPi02760GeV,mesonMassExpectPi0+0.005);
	histoMassMesonPi0->Scale(1000.);
	histoTrueMassMesonPi0->Scale(1000.);
	
	histoMassMesonEtaMinusExp= CalculateMassMinusExpectedMass(histoMassMesonEta,mesonMassExpectEta+0.005);
	histoTrueMassMesonEtaMinusExp = CalculateMassMinusExpectedMass(histoTrueMassMesonEta,mesonMassExpectEta+0.005);
	histoMassMesonEtaMinusExp900GeV= CalculateMassMinusExpectedMass(histoMassMesonEta900GeV,mesonMassExpectEta+0.005);
	histoTrueMassMesonEtaMinusExp900GeV = CalculateMassMinusExpectedMass(histoTrueMassMesonEta900GeV,mesonMassExpectEta+0.005);
	histoMassMesonEtaMinusExp2760GeV= CalculateMassMinusExpectedMass(histoMassMesonEta2760GeV,mesonMassExpectEta+0.005);
	histoTrueMassMesonEtaMinusExp2760GeV = CalculateMassMinusExpectedMass(histoTrueMassMesonEta2760GeV,mesonMassExpectEta+0.005);
	
	
	histoTrueEffPtPi0->Multiply(histoTrueEffPtPi0,histoAccPi0);
	histoTrueEffPtPi0900GeV->Multiply(histoTrueEffPtPi0900GeV,histoAccPi0900GeV);
	histoTrueEffPtPi02760GeV->Multiply(histoTrueEffPtPi02760GeV,histoAccPi02760GeV);

// 	if (patched){
// 		graphCorrectedYieldPi0SysErr->RemovePoint(0);
// 	}
	relSystErrorPi07TeVDown = ExtractRelErrDownAsymmGraph(graphCorrectedYieldPi0SysErr);
	relSystErrorPi07TeVUp = ExtractRelErrUpAsymmGraph(graphCorrectedYieldPi0SysErr);
	nPointsPi0 = graphCorrectedYieldPi0SysErr->GetN();
	cout << "systematic errors Pi0 Conv 7 TeV" << endl;
	for (Int_t i = 0; i < nPointsPi0; i++){
		cout << relSystErrorPi07TeVDown[i] << "\t" << relSystErrorPi07TeVUp[i] << endl;
	}
	
	relSystErrorPi02760GeVDown = ExtractRelErrDownAsymmGraph(graphCorrectedYieldPi0SysErr2760GeV);
	relSystErrorPi02760GeVUp = ExtractRelErrUpAsymmGraph(graphCorrectedYieldPi0SysErr2760GeV);
	nPointsPi02760GeV = graphCorrectedYieldPi0SysErr2760GeV->GetN();
	cout << "systematic errors Pi0 Conv 2760 GeV" << endl;
	for (Int_t i = 0; i < nPointsPi02760GeV; i++){
		cout << relSystErrorPi02760GeVDown[i] << "\t" << relSystErrorPi02760GeVUp[i] << endl;
	}

	relSystErrorPi0900GeVDown = ExtractRelErrDownAsymmGraph(graphCorrectedYieldPi0SysErr900GeV);
	relSystErrorPi0900GeVUp = ExtractRelErrUpAsymmGraph(graphCorrectedYieldPi0SysErr900GeV);
	nPointsPi0900GeV = graphCorrectedYieldPi0SysErr900GeV->GetN();
	cout << "systematic errors Pi0 Conv 900 GeV" << endl;
	for (Int_t i = 0; i < nPointsPi0900GeV; i++){
		cout << relSystErrorPi0900GeVDown[i] << "\t" << relSystErrorPi0900GeVUp[i] << endl;
	}

	relSystErrorEta7TeVDown = ExtractRelErrDownAsymmGraph(graphCorrectedYieldEtaSysErr);
	relSystErrorEta7TeVUp = ExtractRelErrUpAsymmGraph(graphCorrectedYieldEtaSysErr);
	nPointsEta = graphCorrectedYieldEtaSysErr->GetN();
	cout << "systematic errors Eta Conv 7TeV" << endl;
	for (Int_t i = 0; i < nPointsEta; i++){
		cout << relSystErrorEta7TeVDown[i] << "\t" << relSystErrorEta7TeVUp[i] << endl;
	}

	relSystErrorEta2760GeVDown = ExtractRelErrDownAsymmGraph(graphCorrectedYieldEtaSysErr2760GeV);
	relSystErrorEta2760GeVUp = ExtractRelErrUpAsymmGraph(graphCorrectedYieldEtaSysErr2760GeV);
	nPointsEta2760GeV = graphCorrectedYieldEtaSysErr2760GeV->GetN();
	cout << "systematic errors Eta Conv 2.76TeV" << endl;
	for (Int_t i = 0; i < nPointsEta2760GeV; i++){
		cout << relSystErrorEta2760GeVDown[i] << "\t" << relSystErrorEta2760GeVUp[i] << endl;
	}

// 	if (patched){
// 		graphSystErrRatio->RemovePoint(0);
// 		graphSystErrRatio->RemovePoint(0);
// 	}
	TGraphAsymmErrors* dummyGraphSystErr = (TGraphAsymmErrors*)graphSystErrRatio->Clone("dummyGraphSystErr");
	relSystErrorEtaPi07TeVDown = ExtractRelErrDownAsymmGraph(dummyGraphSystErr);
	relSystErrorEtaPi07TeVUp = ExtractRelErrUpAsymmGraph(dummyGraphSystErr);
	nPointsEtaPi07TeV = dummyGraphSystErr->GetN();
	cout << "systematic errors Eta Conv 7TeV" << endl;
	for (Int_t i = 0; i < nPointsEtaPi07TeV; i++){
		cout << relSystErrorEtaPi07TeVDown[i] << "\t" << relSystErrorEtaPi07TeVUp[i] << endl;
	}
	cout << "Ratio Eta/pi0 including syst" << endl;
	relSystErrorEta900GeVDown = ExtractRelErrDownAsymmGraph(graphCorrectedYieldEtaSysErr900GeV);
	relSystErrorEta900GeVUp = ExtractRelErrUpAsymmGraph(graphCorrectedYieldEtaSysErr900GeV);
	nPointsEta900GeV = graphCorrectedYieldEtaSysErr900GeV->GetN();
	cout << "systematic errors Eta Conv 900GeV" << endl;
	for (Int_t i = 0; i < nPointsEta900GeV; i++){
		cout << relSystErrorEta900GeVDown[i] << "\t" << relSystErrorEta900GeVUp[i] << endl;
	}
	
	dummyGraphSystErr = (TGraphAsymmErrors*)graphSystErrRatio2760GeV->Clone("dummyGraphSystErr");
	relSystErrorEtaPi02760GeVDown = ExtractRelErrDownAsymmGraph(dummyGraphSystErr);
	relSystErrorEtaPi02760GeVUp = ExtractRelErrUpAsymmGraph(dummyGraphSystErr);
	nPointsEta2760GeV = dummyGraphSystErr->GetN();
	cout << "systematic errors Eta Conv 2.76TeV" << endl;
	for (Int_t i = 0; i < nPointsEta2760GeV; i++){
		cout << relSystErrorEtaPi02760GeVDown[i] << "\t" << relSystErrorEtaPi02760GeVUp[i] << endl;
	}
	graphRatioEtaPi0ComplErr2760GeV = CalculateSysErrFromRelSysHistoComplete( histoRatioEtaPi02760GeV , "EtaPi0ComplError2760GeV",relSystErrorEtaPi02760GeVDown , relSystErrorEtaPi02760GeVDown, 2, nPointsEta2760GeV);
	graphRatioEtaPi0StatErr2760GeV = new TGraphAsymmErrors(histoRatioEtaPi02760GeV);
	graphRatioEtaPi0StatErr2760GeV->RemovePoint(0);
	cout << "stat Eta/Pi0 2.76 TeV" << endl;
	graphRatioEtaPi0StatErr2760GeV->Print();
	graphRatioEtaPi0SysErr2760GeV = (TGraphAsymmErrors*)graphSystErrRatio2760GeV->Clone("graphRatioEtaPi0SysErr900GeV");

	dummyGraphSystErr = (TGraphAsymmErrors*)graphSystErrRatio900GeV->Clone("dummyGraphSystErr");
	relSystErrorEtaPi0900GeVDown = ExtractRelErrDownAsymmGraph(dummyGraphSystErr);
	relSystErrorEtaPi0900GeVUp = ExtractRelErrUpAsymmGraph(dummyGraphSystErr);
	
	nPointsEta900GeV = dummyGraphSystErr->GetN();
	cout << "systematic errors Eta Conv 900GeV" << endl;
	for (Int_t i = 0; i < nPointsEta900GeV; i++){
		cout << relSystErrorEtaPi0900GeVDown[i] << "\t" << relSystErrorEtaPi0900GeVUp[i] << endl;
	}
	graphRatioEtaPi0ComplErr900GeV = CalculateSysErrFromRelSysHistoComplete( histoRatioEtaPi0900GeV , "EtaPi0ComplError900GeV",relSystErrorEtaPi0900GeVDown , relSystErrorEtaPi0900GeVDown, 2, nPointsEta900GeV);
	graphRatioEtaPi0StatErr900GeV = new TGraphAsymmErrors(histoRatioEtaPi0900GeV);
	graphRatioEtaPi0StatErr900GeV->RemovePoint(0);
	graphRatioEtaPi0SysErr900GeV = (TGraphAsymmErrors*)graphSystErrRatio900GeV->Clone("graphRatioEtaPi0SysErr900GeV");


   TFile* fileTheoryCompilation = new TFile("ExternalInput/TheoryCompilationPP.root");
   TH1F* histoPythia8InvXSection = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Spec2760GeV");
   TH1F* histoPythia8InvXSection_VarBinning = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Spec2760GeVVarBinning");
   graphNLOCalcMuHalf900GeV=    (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf900GeV");
   graphNLOCalcMuOne900GeV=     (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne900GeV");
   graphNLOCalcMuTwo900GeV=     (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo900GeV");
   graphNLOCalcMuHalf2760GeV=    (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf2760GeV");
   graphNLOCalcMuOne2760GeV=     (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne2760GeV");
   graphNLOCalcMuTwo2760GeV=     (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo2760GeV");
   graphNLOCalcMuHalf7TeV=    (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf7000GeV");
   graphNLOCalcMuOne7TeV=     (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne7000GeV");
   graphNLOCalcMuTwo7TeV=     (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo7000GeV");
   graphNLOCalcEtaMuHalf900GeV= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf900GeV");
   graphNLOCalcEtaMuOne900GeV=  (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne900GeV");
   graphNLOCalcEtaMuTwo900GeV=  (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo900GeV"); 
   graphNLOCalcEtaMuHalf2760GeV=   (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf2760GeV");
   graphNLOCalcEtaMuOne2760GeV=    (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne2760GeV");
   graphNLOCalcEtaMuTwo2760GeV=    (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo2760GeV");
   graphNLOCalcEtaMuHalf7TeV=          (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf7000GeV");
   graphNLOCalcEtaMuOne7TeV=           (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne7000GeV");
   graphNLOCalcEtaMuTwo7TeV=           (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo7000GeV");
   graphEtaToPi0NLOMuHalf7TeV =  (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuHalf7000GeV");
   graphEtaToPi0NLOMuOne7TeV =   (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuOne7000GeV");
   graphEtaToPi0NLOMuTwo7TeV =   (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuTwo7000GeV");
   graphEtaToPi0NLOMuHalf2760GeV = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuHalf2760GeV");
   graphEtaToPi0NLOMuOne2760GeV =  (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuOne2760GeV");
   graphEtaToPi0NLOMuTwo2760GeV =  (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuTwo2760GeV");
   graphNLOBKKCalcMuTwo900GeV=     (TGraph*)fileTheoryCompilation->Get("graphNLOCalcBKKInvSecPi0MuTwo900GeV");
   graphNLOBKKCalcMuTwo7TeV =  (TGraph*)fileTheoryCompilation->Get("graphNLOCalcBKKInvSecPi0MuTwo7000GeV");
   graphNLODSSCalcMuTwo900GeV=      (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSSInvSecPi0MuTwo900GeV");
   graphNLODSSCalcMuTwo2760GeV=     (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSSInvSecPi0MuTwo2760GeV");
   graphNLODSSCalcMuTwo7TeV=            (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSSInvSecPi0MuTwo7000GeV");
   
	//************************** Read data for PHOS **************************************************
	filePhos7TeV = 		new TFile(fileNameCaloPhos7TeV);
	histoPi0Phos7TeV = 		(TH1D*)filePhos7TeV->Get("hPi07TeVStat");
	histoPi0PhosSys7TeV = 	(TH1D*)filePhos7TeV->Get("hPi07TeVSys");
	filePhos7TeVEta = 		new TFile(fileNameCaloPhos7TeVEta);
	histoEtaPhos7TeV = 		(TH1D*)filePhos7TeVEta->Get("hEta7TeVStat");
	histoEtaPhosSys7TeV = 	(TH1D*)filePhos7TeVEta->Get("hEta7TeVSys");
	filePhos900GeV = 		new TFile(fileNameCaloPhos900GeV);
	histoPi0Phos900GeV = 	(TH1D*)filePhos900GeV->Get("hPi0900GeVStat");
	histoPi0PhosSys900GeV = 	(TH1D*)filePhos900GeV->Get("hPi0900GeVSys");
   
	filePhos2760GeV = 		new TFile(fileNameCaloPhos2760GeV);
	directoryPHOSPi02760GeV =   (TDirectory*)filePhos2760GeV->Get("pp2760"); 
	histoPi0Phos2760GeV = 		(TH1D*)directoryPHOSPi02760GeV->Get("hPi02760GeVStat");
	histoPi0PhosSys2760GeV = 	(TH1D*)directoryPHOSPi02760GeV->Get("hPi02760GeVSys");
	histoPi0PhosSysRAA2760GeV = 	(TH1D*)directoryPHOSPi02760GeV->Get("hPi02760GeVSysTypeB");
	
	for (Int_t i = 0; i <  histoPi0Phos2760GeV->GetNbinsX()+1 ; i++){
		cout<< histoPi0Phos2760GeV->GetBinCenter(i)<< "\t" << histoPi0Phos2760GeV->GetBinError(i)/histoPi0Phos2760GeV->GetBinContent(i)*100 << "\t" << histoPi0PhosSys2760GeV->GetBinError(i)/histoPi0PhosSys2760GeV->GetBinContent(i)*100 <<"\t"<<  TMath::Sqrt(TMath::Power(histoPi0PhosSys2760GeV->GetBinError(i)/histoPi0PhosSys2760GeV->GetBinContent(i)*100,2) + 4*4) << "\t" << histoPi0PhosSysRAA2760GeV->GetBinError(i)/histoPi0PhosSysRAA2760GeV->GetBinContent(i)*100 << endl;
	}   
   
   
	filePhosEtaToPi0 = 		new TFile(fileNamePHOSEtaToPi0);
	histoEtaToPi0PHOS =		(TH1D*)filePhosEtaToPi0->Get("hEtaPiRatio7TeVStat");
	histoEtaToPi0PHOSSys =	(TH1D*)filePhosEtaToPi0->Get("hEtaPiRatio7TeVSys");
	filePhos7TeVMassData = 	new TFile(fileNamePHOSMassData7TeV);
	filePhos7TeVMassMC = 	new TFile(fileNamePHOSMassMC7TeV);
	histoMassMesonPi0PHOS = 	(TH1D*)filePhos7TeVMassData->Get("Mix_mr1");		
	histoFWHMMesonPi0PHOS= 	(TH1D*)filePhos7TeVMassData->Get("Mix_sr1");		
	histoTrueMassMesonPi0PHOS= (TH1D*)filePhos7TeVMassMC->Get("mass1_GS");		
	histoTrueFWHMMesonPi0PHOS= (TH1D*)filePhos7TeVMassMC->Get("width1_GS");		
	histoMassMesonPi0PHOSMinusExp = CalculateMassMinusExpectedMass(histoMassMesonPi0PHOS,mesonMassExpectPi0+0.010);
	histoTrueMassMesonPi0PHOSMinusExp = CalculateMassMinusExpectedMass(histoTrueMassMesonPi0PHOS,mesonMassExpectPi0+0.010);
	histoFWHMMesonPi0PHOS->Scale(1000.);
	histoTrueFWHMMesonPi0PHOS->Scale(1000.);
	histoMassMesonPi0PHOS->Scale(1000.);
	histoTrueMassMesonPi0PHOS->Scale(1000.);
	cout << "loaded phos" << endl;
	nEvtPHOS7TeV =3.181e+8;
	nEvtPHOS900GeV = 7.3e+6;

	histoForPHOSGraph1 = (TH1D*)histoEtaToPi0PHOSSys->Clone();
	graphSysErrEtaToPi0PHOS = new TGraphAsymmErrors(histoForPHOSGraph1);	
	
	//************************** Read data for EMCAL **************************************************
// 	fileEMCAL7TeV = 		new TFile(fileNameCaloEmcal7TeV);
// 	histoPi0EMCAL7TeV = 	(TH1D*)fileEMCAL7TeV->Get("pi0pro");
// 	histoPi0EMCAL7TeV->SetBinContent(1,NULL);
// 	histoPi0EMCAL7TeV->SetBinError(1,NULL);
	
	fileEMCAL7TeVMassData = 	new TFile(fileNameEMCALMassData7TeVPi0);
	TDirectory* directoryEMCALPi07TeV = (TDirectory*)fileEMCAL7TeVMassData->Get("Pi07TeV");
	histoMassMesonPi0EMCAL = 	(TH1D*)directoryEMCALPi07TeV->Get("MassPi0");		
	histoFWHMMesonPi0EMCAL= 	(TH1D*)directoryEMCALPi07TeV->Get("FWHMPi0MeV");		
	histoTrueMassMesonPi0EMCAL= (TH1D*)directoryEMCALPi07TeV->Get("TrueMassPi0");		
	histoTrueFWHMMesonPi0EMCAL= (TH1D*)directoryEMCALPi07TeV->Get("TrueFWHMPi0MeV");		
	histoMassMesonPi0EMCALMinusExp = CalculateMassMinusExpectedMass(histoMassMesonPi0EMCAL,mesonMassExpectPi0+0.015);
	histoTrueMassMesonPi0EMCALMinusExp = CalculateMassMinusExpectedMass(histoTrueMassMesonPi0EMCAL,mesonMassExpectPi0+0.015);
   for (Int_t  i = 1; i < histoMassMesonPi0EMCAL->GetXaxis()->FindBin(0.8); i++){
      histoMassMesonPi0EMCAL->SetBinContent(i,0.);
      histoFWHMMesonPi0EMCAL->SetBinContent(i,10000.);
      histoTrueFWHMMesonPi0EMCAL->SetBinContent(i,10000.);
      histoTrueMassMesonPi0EMCAL->SetBinContent(i,0.);
   }
   
   for (Int_t  i = histoMassMesonPi0EMCAL->GetXaxis()->FindBin(18.0); i < histoMassMesonPi0EMCAL->GetNbinsX()+1; i++){
      histoMassMesonPi0EMCAL->SetBinContent(i,0.);
      histoFWHMMesonPi0EMCAL->SetBinContent(i,10000.);
      histoTrueFWHMMesonPi0EMCAL->SetBinContent(i,10000.);
      histoTrueMassMesonPi0EMCAL->SetBinContent(i,0.);
   }

//  	histoMassMesonPi0EMCAL->Scale(1000.);
//  	histoTrueMassMesonPi0EMCAL->Scale(1000.);
   
	fileEMCAL7TeVMassDataEta = 	new TFile(fileNameEMCALMassData7TeVEta);
	fileEMCAL7TeVMassMCEta = 	new TFile(fileNameEMCALMassMC7TeVEta);
	histoMassMesonEtaEMCAL = 	(TH1D*)fileEMCAL7TeVMassDataEta->Get("hEtaMass");		
	histoFWHMMesonEtaEMCAL= 		(TH1D*)fileEMCAL7TeVMassDataEta->Get("hEtaWidth");		
	histoTrueMassMesonEtaEMCAL= 	(TH1D*)fileEMCAL7TeVMassMCEta->Get("hEtaMass");		
	histoTrueFWHMMesonEtaEMCAL= 	(TH1D*)fileEMCAL7TeVMassMCEta->Get("hEtaWidth");		
	histoMassMesonEtaEMCALMinusExp = CalculateMassMinusExpectedMass(histoMassMesonEtaEMCAL,mesonMassExpectEta+0.015);
	histoTrueMassMesonEtaEMCALMinusExp = CalculateMassMinusExpectedMass(histoTrueMassMesonEtaEMCAL,mesonMassExpectEta+0.015);
	histoFWHMMesonEtaEMCAL->Scale(1000.);
	histoTrueFWHMMesonEtaEMCAL->Scale(1000.);
	fileFinalResults.open(nameFinalResDat, ios::out);	
	cout << "loaded emcal" << endl;

	filePhosOmega7TeV = 		new TFile(fileNameCaloPhosOmega7TeV);
	graphOmegaPhos7TeV = 		(TGraphAsymmErrors*)filePhosOmega7TeV->Get("graphOmegaStat");
	graphOmegaPhosSys7TeV = 	(TGraphAsymmErrors*)filePhosOmega7TeV->Get("graphOmegaSyst");
	graphOmegaPhosComb7TeV = 	(TGraphAsymmErrors*)filePhosOmega7TeV->Get("graphOmegaComb");
	
	maxPtPi0 = histoPi0Phos7TeV->GetXaxis()->GetBinUpEdge(histoPi0Phos7TeV->GetNbinsX());
	maxPtEta2760GeV = histoInvCrossSectionEta2760GeV->GetXaxis()->GetBinUpEdge(histoInvCrossSectionEta2760GeV->GetNbinsX());
	cout << maxPtPi0 << endl;
   maxPtPi0900GeV= histoPi0Phos900GeV->GetXaxis()->GetBinUpEdge(histoPi0Phos900GeV->GetNbinsX());
   maxPtEta = histoEtaPhos7TeV->GetXaxis()->GetBinUpEdge(histoEtaPhos7TeV->GetNbinsX());
   maxPtEtaPHOS = histoEtaPhos7TeV->GetXaxis()->GetBinUpEdge(histoEtaPhos7TeV->GetNbinsX());
   maxPtPi02760GeV = histoPi0Phos2760GeV->GetXaxis()->GetBinUpEdge(histoPi0Phos2760GeV->GetNbinsX()-1);
	//----------------------- Pi0 7 TeV NLO mu = pt/2 ------------------------------
	graphNLOMuHalfPi07TeV = (TGraph*)graphNLOCalcMuHalf7TeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuHalfPi07TeV, styleMarkerNLOMuHalf, sizeMarkerNLO, colorNLOPi07TeVMuHalf, colorNLOPi07TeVMuHalf );
	//------------------------- Pi0 7 TeV NLO mu = pt ----------------------------
	graphNLOMuOnePi07TeV = (TGraph*)graphNLOCalcMuOne7TeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuOnePi07TeV, styleMarkerNLOMuOne, sizeMarkerNLO, colorNLOPi07TeVMuOne, colorNLOPi07TeVMuOne);
	//------------------------- Pi0 7 TeV NLO mu = 2pt -----------------------------
	graphNLOMuTwoPi07TeV = (TGraph*)graphNLOCalcMuTwo7TeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuTwoPi07TeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOPi07TeVMuTwo, colorNLOPi07TeVMuTwo);
	DrawGammaSetMarkerTGraph(graphNLOBKKCalcMuTwo7TeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOBKKPi07TeVMuTwo, colorNLOBKKPi07TeVMuTwo);
	DrawGammaSetMarkerTGraph(graphNLODSSCalcMuTwo7TeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLODSSPi07TeVMuTwo, colorNLODSSPi07TeVMuTwo);
	
	//----------------------- Pi0 900 TeV NLO mu = pt/2 ------------------------------
	graphNLOMuHalfPi0900GeV = (TGraph*)graphNLOCalcMuHalf900GeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuHalfPi0900GeV, styleMarkerNLOMuHalf, sizeMarkerNLO, colorNLOPi0900GeVMuHalf, colorNLOPi0900GeVMuHalf );
	//------------------------- Pi0 900 TeV NLO mu = pt ----------------------------
	graphNLOMuOnePi0900GeV = (TGraph*)graphNLOCalcMuOne900GeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuOnePi0900GeV, styleMarkerNLOMuOne, sizeMarkerNLO, colorNLOPi0900GeVMuOne, colorNLOPi0900GeVMuOne);
	//------------------------- Pi0 900 TeV NLO mu = 2pt -----------------------------
	graphNLOMuTwoPi0900GeV = (TGraph*)graphNLOCalcMuTwo900GeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuTwoPi0900GeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOPi0900GeVMuTwo, colorNLOPi0900GeVMuTwo);

	Double_t* dummyArrayXValuesNLOBKK = graphNLOBKKCalcMuTwo900GeV->GetX();
	Int_t dummyArrayNBKK = graphNLOBKKCalcMuTwo900GeV->GetN();
	Int_t numberOfBinsBKK = 0;
	Int_t nBKK = 0;
	while (nBKK!=dummyArrayNBKK ){
		if (dummyArrayXValuesNLOBKK[nBKK] <= 10.5){
			nBKK++;
			numberOfBinsBKK++;
			cout << nBKK << endl;
		} else {
			nBKK = dummyArrayNBKK;
		} 
	}
	
	DrawGammaSetMarkerTGraph(graphNLOBKKCalcMuTwo900GeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOBKKPi0900GeVMuTwo, colorNLOBKKPi0900GeVMuTwo);
	for (Int_t i = numberOfBinsBKK; i < dummyArrayNBKK; i++){
		graphNLOBKKCalcMuTwo900GeV->RemovePoint(numberOfBinsBKK);
	}

	DrawGammaSetMarkerTGraph(graphNLODSSCalcMuTwo900GeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLODSSPi07TeVMuTwo, colorNLODSSPi07TeVMuTwo);

	//----------------------- Pi0 2760 GeV NLO mu = pt/2 ------------------------------
	graphNLOMuHalfPi02760GeV = (TGraph*)graphNLOCalcMuHalf2760GeV->Clone();
	cout << "hier" << endl;
	DrawGammaSetMarkerTGraph(graphNLOMuHalfPi02760GeV, styleMarkerNLOMuHalf, sizeMarkerNLO, colorNLOPi02760GeVMuHalf, colorNLOPi02760GeVMuHalf );
	//------------------------- Pi0 2760 GeV NLO mu = pt ----------------------------
	graphNLOMuOnePi02760GeV = (TGraph*)graphNLOCalcMuOne2760GeV->Clone();
	cout << "hier" << endl;
	DrawGammaSetMarkerTGraph(graphNLOMuOnePi02760GeV, styleMarkerNLOMuOne, sizeMarkerNLO, colorNLOPi02760GeVMuOne, colorNLOPi02760GeVMuOne);
	//------------------------- Pi0 2760 GeV NLO mu = 2pt -----------------------------
	graphNLOMuTwoPi02760GeV = (TGraph*)graphNLOCalcMuTwo2760GeV->Clone();
	cout << "hier" << endl;
	DrawGammaSetMarkerTGraph(graphNLOMuTwoPi02760GeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOPi02760GeVMuTwo, colorNLOPi02760GeVMuTwo);
/*
	DrawGammaSetMarkerTGraph(graphNLODSSCalcMuTwo2760GeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLODSSPi07TeVMuTwo, colorNLODSSPi07TeVMuTwo);
	graphNLODSSCalcMuTwo2760GeV->RemovePoint(0);*/



	//----------------------- Eta 7TeV NLO mu = pt/2 ------------------------------
	graphNLOMuHalfEta7TeV = (TGraph*)graphNLOCalcEtaMuHalf7TeV->Clone();
	Double_t* dummyArrayXValuesNLO = graphNLOMuHalfEta7TeV->GetX();
	Int_t dummyArrayN = graphNLOMuHalfEta7TeV->GetN();
	Int_t numberOfBins = 0;
	Int_t n = 0;
	while (n!=dummyArrayN ){
		if (dummyArrayXValuesNLO[n] <= 15.5){
			n++;
			numberOfBins++;
			cout << n << endl;
		} else {
			n = dummyArrayN;
		} 
	}
	cout << "Number of bins Eta 7TeV "<<numberOfBins << endl;
	DrawGammaSetMarkerTGraph(graphNLOMuHalfEta7TeV, styleMarkerNLOMuHalf, sizeMarkerNLO, colorNLOEta7TeVMuHalf, colorNLOEta7TeVMuHalf );
	for (Int_t i = numberOfBins; i < dummyArrayN; i++){
		graphNLOMuHalfEta7TeV->RemovePoint(numberOfBins);
	}

	//------------------------- Eta 7TeV NLO mu = pt ----------------------------
	graphNLOMuOneEta7TeV = (TGraph*)graphNLOCalcEtaMuOne7TeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuOneEta7TeV,  styleMarkerNLOMuOne, sizeMarkerNLO, colorNLOEta7TeVMuOne, colorNLOEta7TeVMuOne);
	for (Int_t i = numberOfBins; i < dummyArrayN; i++){
		graphNLOMuOneEta7TeV->RemovePoint(numberOfBins);
	}
	
	//------------------------- Eta 7TeV NLO mu = 2pt -----------------------------
	graphNLOMuTwoEta7TeV = (TGraph*)graphNLOCalcEtaMuTwo7TeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuTwoEta7TeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOEta7TeVMuTwo, colorNLOEta7TeVMuTwo);
	for (Int_t i = numberOfBins; i < dummyArrayN; i++){
		graphNLOMuTwoEta7TeV->RemovePoint(numberOfBins);
	}
	
	//----------------------- Eta 2760GeV NLO mu = pt/2 ------------------------------
	graphNLOMuHalfEta2760GeV = (TGraph*)graphNLOCalcEtaMuHalf2760GeV->Clone();
	dummyArrayXValuesNLO = graphNLOMuHalfEta2760GeV->GetX();
	dummyArrayN = graphNLOMuHalfEta2760GeV->GetN();
	numberOfBins = 0;
	n = 0;
	while (n!=dummyArrayN ){
		if (dummyArrayXValuesNLO[n] <= 7.5){
			n++;
			numberOfBins++;
			cout << n << endl;
		} else {
			n = dummyArrayN;
		} 
	}
	cout << "Number of bins Eta 2.76TeV "<<numberOfBins << endl;
	graphNLOMuHalfEta2760GeV = ScaleGraph(graphNLOMuHalfEta2760GeV,xSection2760GeV*recalcBarn);
	DrawGammaSetMarkerTGraph(graphNLOMuHalfEta2760GeV, styleMarkerNLOMuHalf, sizeMarkerNLO, colorNLOEta2760GeVMuHalf, colorNLOEta2760GeVMuHalf );
	for (Int_t i = numberOfBins; i < dummyArrayN; i++){
		graphNLOMuHalfEta2760GeV->RemovePoint(numberOfBins);
	}
	graphNLOMuHalfEta2760GeV->RemovePoint(0);
	graphNLOMuHalfEta2760GeV->RemovePoint(0);
	//------------------------- Eta 2760GeV NLO mu = pt ----------------------------
	graphNLOMuOneEta2760GeV = (TGraph*)graphNLOCalcEtaMuOne2760GeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuOneEta2760GeV,  styleMarkerNLOMuOne, sizeMarkerNLO, colorNLOEta2760GeVMuOne, colorNLOEta2760GeVMuOne);
	for (Int_t i = numberOfBins; i < dummyArrayN; i++){
		graphNLOMuOneEta2760GeV->RemovePoint(numberOfBins);
	}
	//------------------------- Eta 2760GeV NLO mu = 2pt -----------------------------
	graphNLOMuTwoEta2760GeV = (TGraph*)graphNLOCalcEtaMuTwo2760GeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuTwoEta2760GeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOEta2760GeVMuTwo, colorNLOEta2760GeVMuTwo);
	for (Int_t i = numberOfBins; i < dummyArrayN; i++){
		graphNLOMuTwoEta2760GeV->RemovePoint(numberOfBins);
	}

	//----------------------- Eta 900GeV NLO mu = pt/2 ------------------------------
	graphNLOMuHalfEta900GeV = (TGraph*)graphNLOCalcEtaMuHalf900GeV->Clone();
	dummyArrayXValuesNLO = graphNLOMuHalfEta900GeV->GetX();
	dummyArrayN = graphNLOMuHalfEta900GeV->GetN();
	numberOfBins = 0;
	n = 0;
	while (n!=dummyArrayN ){
		if (dummyArrayXValuesNLO[n] <= 4.5){
			n++;
			numberOfBins++;
			cout << n << endl;
		} else {
			n = dummyArrayN;
		} 
	}
	cout << "Number of bins Eta 0.9TeV "<<numberOfBins << endl;
   
	DrawGammaSetMarkerTGraph(graphNLOMuHalfEta900GeV, styleMarkerNLOMuHalf, sizeMarkerNLO, colorNLOEta900GeVMuHalf, colorNLOEta900GeVMuHalf );
	for (Int_t i = numberOfBins; i < dummyArrayN; i++){
		graphNLOMuHalfEta900GeV->RemovePoint(numberOfBins);
	}
	//------------------------- Eta 900GeV NLO mu = pt ----------------------------
	graphNLOMuOneEta900GeV = (TGraph*)graphNLOCalcEtaMuOne900GeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuOneEta900GeV,  styleMarkerNLOMuOne, sizeMarkerNLO, colorNLOEta900GeVMuOne, colorNLOEta900GeVMuOne);
	for (Int_t i = numberOfBins; i < dummyArrayN; i++){
		graphNLOMuOneEta900GeV->RemovePoint(numberOfBins);
	}
	//------------------------- Eta 900GeV NLO mu = 2pt -----------------------------
	graphNLOMuTwoEta900GeV = (TGraph*)graphNLOCalcEtaMuTwo900GeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuTwoEta900GeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOEta900GeVMuTwo, colorNLOEta900GeVMuTwo);
	for (Int_t i = numberOfBins; i < dummyArrayN; i++){
		graphNLOMuTwoEta900GeV->RemovePoint(numberOfBins);
	}
	
	
	
	
	//*********************** Comparison of Measuremnts in ALICE ***************************************************
	TCanvas* canvasRatioEtaPi0ALICE = new TCanvas("canvasRatioEtaPi0ALICE","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaPi0ALICE, 0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioEtaPi0ALICE = new TH2D("histo2DRatioEtaPi0ALICE", "histo2DRatioEtaPi0ALICE", 20,0.,maxPtEta,1000.,-0.4,1.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0ALICE, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.8,0.65, 510, 510);
	histo2DRatioEtaPi0ALICE->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DRatioEtaPi0ALICE->Draw();
	
	TGraphAsymmErrors *graphSysErrEtaToPi0PHOSNewXError = (TGraphAsymmErrors*)graphSysErrEtaToPi0PHOS->Clone();
	ProduceGraphAsymmFixedXErrors(graphSysErrEtaToPi0PHOSNewXError, 0.4);
	graphSysErrEtaToPi0PHOSNewXError->SetFillColor(colorPHOSSyst);
     graphSysErrEtaToPi0PHOSNewXError->Draw("same,2,p");

	TGraphAsymmErrors* graphSystErrRatioNewXError = (TGraphAsymmErrors*)graphSystErrRatio->Clone();
	ProduceGraphAsymmFixedXErrors(graphSystErrRatioNewXError, 0.3);
	graphSystErrRatioNewXError->SetFillColor(colorConvSyst);
	graphSystErrRatioNewXError->Draw("same,2,p");

//   DrawGammaSetMarkerTGraphErr(graphEtaToPi0EMCALSys7TeV, markerStyleEMCAL, markerSizeInvYield, colorEMCAL, colorEMCAL);
//   graphEtaToPi0EMCALSys7TeV->SetFillColor(colorEMCALSyst);
//   graphEtaToPi0EMCALSys7TeV->SetFillStyle(fillStyleEMCAL);
//   graphEtaToPi0EMCALSys7TeV->Draw("p,2same");
	
//   DrawGammaSetMarkerTGraphErr(graphEtaToPi0EMCAL7TeV, markerStyleEMCAL,markerSizeInvYield,colorEMCAL,colorEMCAL);
//   graphEtaToPi0EMCAL7TeV->Draw("p,same");
		
	DrawGammaSetMarker(histoEtaToPi0Pythia7TeV, markerStyleMCEtaToPi0, markerSizeInvYield, colorPythiaEtaToPi0, colorPythiaEtaToPi0);
	//histoEtaToPi0Pythia7TeV->Draw("same");
	DrawGammaSetMarker(histoEtaToPi0Phojet7TeV, markerStyleMCEtaToPi0, markerSizeInvYield, colorPhojetEtaToPi0,colorPhojetEtaToPi0);
	//histoEtaToPi0Phojet7TeV->Draw("same");
	
	graphEtaToPi0NLOMuHalf7TeV->SetLineColor(colorNLOPi07TeVMuHalf);
	graphEtaToPi0NLOMuHalf7TeV->SetLineStyle(styleLineNLOMuHalf);
	graphEtaToPi0NLOMuHalf7TeV->SetLineWidth(widthLineNLO);
	graphEtaToPi0NLOMuOne7TeV->SetLineColor(colorNLOPi07TeVMuOne);
	graphEtaToPi0NLOMuOne7TeV->SetLineStyle(styleLineNLOMuOne);
	graphEtaToPi0NLOMuOne7TeV->SetLineWidth(widthLineNLO);	
	graphEtaToPi0NLOMuTwo7TeV->SetLineColor(colorNLOPi07TeVMuTwo);
	graphEtaToPi0NLOMuTwo7TeV->SetLineStyle(styleLineNLOMuTwo);
	graphEtaToPi0NLOMuTwo7TeV->SetLineWidth(widthLineNLO);

	
	DrawGammaSetMarker(histoRatioEtaPi0, markerStyleConv, markerSizeInvYield, colorConv, colorConv);
	histoRatioEtaPi0->DrawCopy("e1,x0,same"); 
		
	TGraphErrors *graphhistoEtaToPi0PHOS =  new TGraphErrors(histoEtaToPi0PHOS);              
	TGraphErrors *graphhistoEtaToPi0PHOSDisplacedX = (TGraphErrors*)graphhistoEtaToPi0PHOS->Clone();
        ProduceGraphErrDisplacedX(graphhistoEtaToPi0PHOSDisplacedX,0.02);
	DrawGammaSetMarkerTGraphErr(graphhistoEtaToPi0PHOSDisplacedX, markerStylePHOS, markerSizeInvYield, colorPHOS, colorPHOS); 
	graphhistoEtaToPi0PHOSDisplacedX->Draw("same,p");  

	TLegend* legendRatioALICE = new TLegend(0.8,0.16,0.977,0.28);
	legendRatioALICE->SetTextSize(0.05);			
	legendRatioALICE->SetFillColor(0);
	legendRatioALICE->SetBorderSize(0);
	legendRatioALICE->AddEntry(histoRatioEtaPi0,Form("PCM"),"p");
	legendRatioALICE->AddEntry(graphhistoEtaToPi0PHOSDisplacedX,Form("PHOS"),"pe");
	legendRatioALICE->Draw();
	
	canvasRatioEtaPi0ALICE->Update();
	canvasRatioEtaPi0ALICE->SaveAs(Form("%s/%s_Pi0EtaRatioALICE_7TeV_Paper.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	if(!thesis)DrawAliceLogoCombinedPreliminary(pictDrawingCoordinatesPi0Eta[0], pictDrawingCoordinatesPi0Eta[1], pictDrawingCoordinatesPi0Eta[2], pictDrawingCoordinatesPi0Eta[3], pictDrawingCoordinatesPi0Eta[4], pictDrawingCoordinatesPi0Eta[5], pictDrawingCoordinatesPi0Eta[6], pictDrawingCoordinatesPi0Eta[7], pictDrawingCoordinates[8],collisionSystem7TeV, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900);
	
	canvasRatioEtaPi0ALICE->Update();
	canvasRatioEtaPi0ALICE->SaveAs(Form("%s/%s_Pi0EtaRatioALICE_7TeV.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	delete canvasRatioEtaPi0ALICE;
	delete legendRatioALICE;
	
	graphCombinedEtaToPi0 = CombinePtPointsSpectra(	histoRatioEtaPi0, 		graphSystErrRatio,
												histoEtaToPi0PHOS, 	graphSysErrEtaToPi0PHOS,
												graphStatErrCombinedEtaToPi0, graphSysErrCombinedEtaToPi0,
												xPtLimitsEta7TeV, 14, 1, 0,2);

	TCanvas* canvasRatioEtaPi0NLO = new TCanvas("canvasRatioEtaPi0NLO","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaPi0NLO, 0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioEtaPi0Theory = new TH2D("histo2DRatioEtaPi0Theory", "histo2DRatioEtaPi0Theory", 20,0.,maxPtEta,1000.,-0.4,1.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0Theory, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.8,0.65, 510, 510);
	histo2DRatioEtaPi0Theory->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DRatioEtaPi0Theory->Draw();
	
	TGraphAsymmErrors* graphCombinedEtaToPi0NewXError = (TGraphAsymmErrors*)graphCombinedEtaToPi0->Clone();
	ProduceGraphAsymmFixedXErrors(graphCombinedEtaToPi0NewXError, 0.3);
	DrawGammaSetMarkerTGraphAsym(graphCombinedEtaToPi0NewXError, markerStyleConv,markerSizeInvYield, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeVBox);
     graphCombinedEtaToPi0NewXError->SetFillColor(colorCommonSpectrumPi07TeVBox);
     graphCombinedEtaToPi0NewXError->SetFillStyle(1001);
	graphCombinedEtaToPi0NewXError->Draw("p,E2same");
	
	DrawGammaSetMarker(histoEtaToPi0Pythia900GeV, markerStyleMCEtaToPi0, markerSizeInvYield-0.2, colorPythiaEtaToPi0900GeV, colorPythiaEtaToPi0900GeV);
	DrawGammaSetMarker(histoEtaToPi0Phojet900GeV, markerStyleMCEtaToPi0, markerSizeInvYield-0.2, colorPhojetEtaToPi0900GeV, colorPhojetEtaToPi0900GeV);
	graphEtaToPi0NLOMuHalf7TeV->Draw("same,c");
	graphEtaToPi0NLOMuOne7TeV->Draw("same,c");
	graphEtaToPi0NLOMuTwo7TeV->Draw("same,c");
	
	histo2DRatioEtaPi0Theory->Draw("AXIS,same");
	TLegend* legendRatioEtaToPi0Theory = new TLegend(0.2,0.16,0.65,0.33);
	legendRatioEtaToPi0Theory->SetTextSize(0.04);			
	legendRatioEtaToPi0Theory->SetFillColor(0);
	legendRatioEtaToPi0Theory->SetBorderSize(0);
	legendRatioEtaToPi0Theory->AddEntry(graphCombinedEtaToPi0NewXError,Form("stat + syst, PCM, PHOS, #sqrt{#it{s}} = 7 TeV "),"fp");
	legendRatioEtaToPi0Theory->AddEntry(graphEtaToPi0NLOMuHalf7TeV,"NLO #mu = 0.5 #it{p}_{T}","l");
	legendRatioEtaToPi0Theory->AddEntry(graphEtaToPi0NLOMuOne7TeV,"NLO #mu = #it{p}_{T}","l");
	legendRatioEtaToPi0Theory->AddEntry(graphEtaToPi0NLOMuTwo7TeV,"NLO #mu = 2 #it{p}_{T}","l");
	legendRatioEtaToPi0Theory->Draw();
	
	canvasRatioEtaPi0NLO->Update();
	canvasRatioEtaPi0NLO->SaveAs(Form("%s/%s_Pi0EtaRatioTheory_Paper.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	
	canvasRatioEtaPi0NLO->SetLogx(0);
	DrawAliceLogoSimple(pictDrawingCoordinatesPi0Eta[4], pictDrawingCoordinatesPi0Eta[5], pictDrawingCoordinatesPi0Eta[6],1350,900);
	canvasRatioEtaPi0NLO->Update();
	canvasRatioEtaPi0NLO->SaveAs(Form("%s/%s_Pi0EtaRatioTheory.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	
	TCanvas* canvasRatioEtaPi0NLOSep = new TCanvas("canvasRatioEtaPi0NLOSep","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaPi0NLOSep, 0.09, 0.01, 0.015, 0.115);
	
	histo2DRatioEtaPi0Theory->Draw();
	
	TGraphAsymmErrors* graphStatCombinedEtaToPi0NewXError = (TGraphAsymmErrors*)graphStatErrCombinedEtaToPi0->Clone();
	ProduceGraphAsymmFixedXErrors(graphStatCombinedEtaToPi0NewXError, 0.);
	DrawGammaSetMarkerTGraphAsym(graphStatCombinedEtaToPi0NewXError, markerStyleConv,markerSizeInvYield, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV);
	TGraphAsymmErrors* graphSysCombinedEtaToPi0NewXError = (TGraphAsymmErrors*)graphSysErrCombinedEtaToPi0->Clone();
	ProduceGraphAsymmFixedXErrors(graphSysCombinedEtaToPi0NewXError, 0.3);
	DrawGammaSetMarkerTGraphAsym(graphSysCombinedEtaToPi0NewXError, markerStyleConv,markerSizeInvYield, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeVBox);
	graphSysCombinedEtaToPi0NewXError->SetFillColor(colorCommonSpectrumPi07TeVBox);
	graphSysCombinedEtaToPi0NewXError->SetFillStyle(1001);
	graphSysCombinedEtaToPi0NewXError->Draw("2same");
	graphStatCombinedEtaToPi0NewXError->Draw("p,same");
	
	graphEtaToPi0NLOMuHalf7TeV->Draw("same,c");
	graphEtaToPi0NLOMuOne7TeV->Draw("same,c");
	graphEtaToPi0NLOMuTwo7TeV->Draw("same,c");
	histo2DRatioEtaPi0Theory->Draw("AXIS,same");
	
	TLegend* legendRatioEtaToPi0TheorySep = new TLegend(0.2,0.16,0.65,0.33);
	legendRatioEtaToPi0TheorySep->SetTextSize(0.04);			
	legendRatioEtaToPi0TheorySep->SetFillColor(0);
	legendRatioEtaToPi0TheorySep->SetBorderSize(0);
	legendRatioEtaToPi0TheorySep->AddEntry(graphCombinedEtaToPi0NewXError,Form("ALICE, #sqrt{#it{s}} = 7 TeV"),"fp");
	legendRatioEtaToPi0TheorySep->AddEntry(graphEtaToPi0NLOMuHalf7TeV,"NLO #mu = 0.5 #it{p}_{T}","l");
	legendRatioEtaToPi0TheorySep->AddEntry(graphEtaToPi0NLOMuOne7TeV,"NLO #mu = #it{p}_{T}","l");
	legendRatioEtaToPi0TheorySep->AddEntry(graphEtaToPi0NLOMuTwo7TeV,"NLO #mu = 2 #it{p}_{T}","l");
	legendRatioEtaToPi0TheorySep->Draw();
	
	canvasRatioEtaPi0NLOSep->Update();
	canvasRatioEtaPi0NLOSep->SaveAs(Form("%s/%s_Pi0EtaRatioTheorySep_Paper.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	canvasRatioEtaPi0NLOSep->SetLogx(0);
	DrawAliceLogoSimple(pictDrawingCoordinatesPi0Eta[4], pictDrawingCoordinatesPi0Eta[5], pictDrawingCoordinatesPi0Eta[6],1350,900);
	canvasRatioEtaPi0NLOSep->Update();
	canvasRatioEtaPi0NLOSep->SaveAs(Form("%s/%s_Pi0EtaRatioTheorySep.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	
	TCanvas* canvasRatioEtaPi0DiffEnergies = new TCanvas("canvasRatioEtaPi0DiffEnergies","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaPi0DiffEnergies,0.09, 0.01, 0.015, 0.115);
	
	
	
	TH2D *histo2DRatioEtaPi0DiffEnergies = new TH2D("histo2DRatioEtaPi0DiffEnergies", "histo2DRatioEtaPi0DiffEnergies", 20,0.,maxPtEta,1000.,-0.4,1.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0DiffEnergies, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.05,0.064, 0.05,0.064, 0.8,0.65, 510, 510);
	histo2DRatioEtaPi0DiffEnergies->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DRatioEtaPi0DiffEnergies->Draw();
	
	TGraphAsymmErrors* graphRatioEtaPi0ComplErr900GeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0ComplErr900GeV->Clone();
	ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0ComplErr900GeVNewXError, 0.5);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0ComplErr900GeVNewXError, markerStyleConv, markerSizeInvYield, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeVBox);
   graphRatioEtaPi0ComplErr900GeVNewXError->SetFillColor(colorCommonSpectrumPi0900GeVBox);
   graphRatioEtaPi0ComplErr900GeVNewXError->SetFillStyle(1001);
	graphRatioEtaPi0ComplErr900GeVNewXError->Draw("p,E2same");


	TGraphAsymmErrors* graphRatioEtaPi0ComplErr2760GeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0ComplErr2760GeV->Clone();
	ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0ComplErr2760GeVNewXError, 0.4);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0ComplErr2760GeVNewXError, markerStyleConv, markerSizeInvYield, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeVBox);
   graphRatioEtaPi0ComplErr2760GeVNewXError->SetFillColor(colorCommonSpectrumPi02760GeVBox);
   graphRatioEtaPi0ComplErr2760GeVNewXError->SetFillStyle(1001);
	graphRatioEtaPi0ComplErr2760GeVNewXError->Draw("p,E2same");

	graphCombinedEtaToPi0NewXError->Draw("p,E2same");

	TGraphAsymmErrors* graphRatioEtaPi0ComplErr900GeVNoError = (TGraphAsymmErrors*)graphRatioEtaPi0ComplErr900GeV->Clone();
	ProduceGraphAsymmWithoutXYErrors(graphRatioEtaPi0ComplErr900GeVNoError);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0ComplErr900GeVNoError, markerStyleConv, markerSizeInvYield, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV);
	
	TGraphAsymmErrors* graphRatioEtaPi0ComplErr2760GeVNoError = (TGraphAsymmErrors*)graphRatioEtaPi0ComplErr2760GeV->Clone();
	ProduceGraphAsymmWithoutXYErrors(graphRatioEtaPi0ComplErr2760GeVNoError);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0ComplErr2760GeVNoError, markerStyleConv, markerSizeInvYield, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV);

	TGraphAsymmErrors* graphCombinedEtaToPi0NoError = (TGraphAsymmErrors*)graphCombinedEtaToPi0->Clone();
	ProduceGraphAsymmWithoutXYErrors(graphCombinedEtaToPi0NoError);
	DrawGammaSetMarkerTGraphAsym(graphCombinedEtaToPi0NoError, markerStyleConv,markerSizeInvYield, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV);


	graphCombinedEtaToPi0NoError->Draw("p,same");
	graphRatioEtaPi0ComplErr2760GeVNoError->Draw("p,same");
	graphRatioEtaPi0ComplErr900GeVNoError->Draw("p,same");

	TLegend* legendRatioEtaToPi0DiffEnergies = new TLegend(0.27,0.15,0.7,0.29);
	legendRatioEtaToPi0DiffEnergies->SetTextSize(0.04);			
	legendRatioEtaToPi0DiffEnergies->SetFillColor(0);
	legendRatioEtaToPi0DiffEnergies->SetBorderSize(0);
	legendRatioEtaToPi0DiffEnergies->AddEntry(graphCombinedEtaToPi0NewXError,Form("stat + syst, ALICE, #sqrt{#it{s}} = 7 TeV (PLB 717 (2012) 162-172)"),"fp");
	legendRatioEtaToPi0DiffEnergies->AddEntry(graphRatioEtaPi0ComplErr2760GeVNewXError,Form("stat + syst, ALICE (PCM), #sqrt{#it{s}} = 2.76 TeV"),"fp");
	legendRatioEtaToPi0DiffEnergies->AddEntry(graphRatioEtaPi0ComplErr900GeVNewXError,Form("stat + syst, ALICE (PCM), #sqrt{#it{s}} = 0.9 TeV"),"fp");
	legendRatioEtaToPi0DiffEnergies->Draw();
	
	canvasRatioEtaPi0DiffEnergies->Update();
	canvasRatioEtaPi0DiffEnergies->SaveAs(Form("%s/%s_Pi0EtaRatioDiffEnergiesALICE_Paper.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));

	TCanvas* canvasRatioEtaPi0DiffEnergiesSep = new TCanvas("canvasRatioEtaPi0DiffEnergiesSep","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaPi0DiffEnergiesSep,0.09, 0.01, 0.015, 0.115);
	
	
	
	TH2D *histo2DRatioEtaPi0DiffEnergiesSep = new TH2D("histo2DRatioEtaPi0DiffEnergiesSep", "histo2DRatioEtaPi0DiffEnergiesSep", 20,0.,maxPtEta,1000.,-0.4,1.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0DiffEnergiesSep, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.05,0.064, 0.05,0.064, 0.8,0.65, 510, 510);
	histo2DRatioEtaPi0DiffEnergiesSep->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DRatioEtaPi0DiffEnergiesSep->Draw();
	
	TGraphAsymmErrors* graphRatioEtaPi0SysErr900GeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0SysErr900GeV->Clone();
	ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0SysErr900GeVNewXError, 0.5);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0SysErr900GeVNewXError, markerStyleConv, markerSizeInvYield, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeVBox);
   graphRatioEtaPi0SysErr900GeVNewXError->SetFillColor(colorCommonSpectrumPi0900GeVBox);
   graphRatioEtaPi0SysErr900GeVNewXError->SetFillStyle(1001);
	graphRatioEtaPi0SysErr900GeVNewXError->Draw("2same");


	TGraphAsymmErrors* graphRatioEtaPi0SysErr2760GeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0SysErr2760GeV->Clone();
	ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0SysErr2760GeVNewXError, 0.4);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0SysErr2760GeVNewXError, markerStyleConv, markerSizeInvYield, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeVBox);
   graphRatioEtaPi0SysErr2760GeVNewXError->SetFillColor(colorCommonSpectrumPi02760GeVBox);
   graphRatioEtaPi0SysErr2760GeVNewXError->SetFillStyle(1001);
	graphRatioEtaPi0SysErr2760GeVNewXError->Draw("2same");

	graphSysCombinedEtaToPi0NewXError->Draw("2same");

	TGraphAsymmErrors* graphRatioEtaPi0StatErr900GeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0StatErr900GeV->Clone();
	ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0StatErr900GeVNewXError, 0.);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0StatErr900GeVNewXError, markerStyleConv, markerSizeInvYield, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV);
   graphRatioEtaPi0StatErr900GeVNewXError->Draw("p,same");

	TGraphAsymmErrors* graphRatioEtaPi0StatErr2760GeVNewXError = (TGraphAsymmErrors*)graphRatioEtaPi0StatErr2760GeV->Clone();
	ProduceGraphAsymmFixedXErrors(graphRatioEtaPi0StatErr2760GeVNewXError, 0.);
	DrawGammaSetMarkerTGraphAsym(graphRatioEtaPi0StatErr2760GeVNewXError, markerStyleConv, markerSizeInvYield, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV);
   graphRatioEtaPi0StatErr2760GeVNewXError->Draw("p,same");

	graphStatCombinedEtaToPi0NewXError->Draw("p,same");
	
	TLegend* legendRatioEtaToPi0DiffEnergiesSep = new TLegend(0.27,0.15,0.7,0.29);
	legendRatioEtaToPi0DiffEnergiesSep->SetTextSize(0.04);			
	legendRatioEtaToPi0DiffEnergiesSep->SetFillColor(0);
	legendRatioEtaToPi0DiffEnergiesSep->SetBorderSize(0);
	legendRatioEtaToPi0DiffEnergiesSep->AddEntry(graphCombinedEtaToPi0NewXError,Form("ALICE, #sqrt{#it{s}} = 7 TeV (PLB 717 (2012) 162-172)"),"fp");
	legendRatioEtaToPi0DiffEnergiesSep->AddEntry(graphRatioEtaPi0ComplErr2760GeVNewXError,Form("ALICE (PCM), #sqrt{#it{s}} = 2.76 TeV"),"fp");
	legendRatioEtaToPi0DiffEnergiesSep->AddEntry(graphRatioEtaPi0ComplErr900GeVNewXError,Form("ALICE (PCM), #sqrt{#it{s}} = 0.9 TeV"),"fp");
	legendRatioEtaToPi0DiffEnergiesSep->Draw();
	
	canvasRatioEtaPi0DiffEnergiesSep->Update();
	canvasRatioEtaPi0DiffEnergiesSep->SaveAs(Form("%s/%s_Pi0EtaRatioDiffEnergiesSepALICE_Paper.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	if(!thesis)DrawAliceLogoCombinedPreliminary(pictDrawingCoordinatesPi0Eta[0], pictDrawingCoordinatesPi0Eta[1], pictDrawingCoordinatesPi0Eta[2], pictDrawingCoordinatesPi0Eta[3], pictDrawingCoordinatesPi0Eta[4], pictDrawingCoordinatesPi0Eta[5], pictDrawingCoordinatesPi0Eta[6], pictDrawingCoordinatesPi0Eta[7], pictDrawingCoordinates[8], collisionSystemCombinedReallyAll, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900);
	
	canvasRatioEtaPi0DiffEnergiesSep->Update();
	canvasRatioEtaPi0DiffEnergiesSep->SaveAs(Form("%s/%s_Pi0EtaRatioDiffEnergiesSepALICE.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	delete canvasRatioEtaPi0DiffEnergiesSep;
	delete legendRatioEtaToPi0DiffEnergiesSep;

	TCanvas* canvasRatioEtaPi0NLO2760GeV = new TCanvas("canvasRatioEtaPi0NLO2760GeV","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaPi0NLO2760GeV,0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioEtaPi0Theory2760GeV = new TH2D("histo2DRatioEtaPi0Theory2760GeV", "histo2DRatioEtaPi0Theory2760GeV", 20,0.,maxPtEta2760GeV,1000.,-0.4,1.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0Theory2760GeV, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.8,0.65, 510, 510);
	histo2DRatioEtaPi0Theory2760GeV->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DRatioEtaPi0Theory2760GeV->Draw();
	
	graphRatioEtaPi0ComplErr2760GeVNewXError->Draw("p,E2same");
	DrawGammaNLOTGraph(graphEtaToPi0NLOMuHalf2760GeV, widthLineNLO, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	DrawGammaNLOTGraph(graphEtaToPi0NLOMuOne2760GeV, widthLineNLO, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	DrawGammaNLOTGraph(graphEtaToPi0NLOMuTwo2760GeV, widthLineNLO, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);
	
	graphEtaToPi0NLOMuHalf2760GeV->Draw("same,c");
	graphEtaToPi0NLOMuOne2760GeV->Draw("same,c");
	graphEtaToPi0NLOMuTwo2760GeV->Draw("same,c");
	
	histo2DRatioEtaPi0Theory2760GeV->Draw("AXIS,same");
	
	TLegend* legendRatioEtaToPi0Theory2760GeV = new TLegend(0.5,0.16,0.977,0.34);
	legendRatioEtaToPi0Theory2760GeV->SetTextSize(0.036);			
	legendRatioEtaToPi0Theory2760GeV->SetFillColor(0);
	legendRatioEtaToPi0Theory2760GeV->SetBorderSize(0);
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphRatioEtaPi0ComplErr2760GeVNewXError,Form("stat + syst, PCM, #sqrt{#it{s}} = 2.76 TeV "),"fp");
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphEtaToPi0NLOMuHalf2760GeV,"NLO #mu = 0.5 #it{p}_{T}","l");
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphEtaToPi0NLOMuOne2760GeV,"NLO #mu = #it{p}_{T}","l");
	legendRatioEtaToPi0Theory2760GeV->AddEntry(graphEtaToPi0NLOMuTwo2760GeV,"NLO #mu = 2 #it{p}_{T}","l");
	legendRatioEtaToPi0Theory2760GeV->Draw();
	
	canvasRatioEtaPi0NLO2760GeV->Update();
	canvasRatioEtaPi0NLO2760GeV->SaveAs(Form("%s/%s_Pi0EtaRatioTheory2760GeV_Paper.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	if(!thesis)DrawAliceLogoCombinedPreliminary(pictDrawingCoordinatesPi0Eta[0], pictDrawingCoordinatesPi0Eta[1], pictDrawingCoordinatesPi0Eta[2], pictDrawingCoordinatesPi0Eta[3], pictDrawingCoordinatesPi0Eta[4], pictDrawingCoordinatesPi0Eta[5], pictDrawingCoordinatesPi0Eta[6], pictDrawingCoordinatesPi0Eta[7], pictDrawingCoordinates[8],collisionSystem2760GeV, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900);
	
	canvasRatioEtaPi0NLO2760GeV->Update();
	canvasRatioEtaPi0NLO2760GeV->SaveAs(Form("%s/%s_Pi0EtaRatioTheory2760GeV.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	TCanvas* canvasRatioEtaPi0NLO2760GeVSep = new TCanvas("canvasRatioEtaPi0NLO2760GeVSep","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaPi0NLO2760GeVSep,0.09, 0.01, 0.015, 0.115);
	
	TH2D *histo2DRatioEtaPi0TheorySep2760GeV = new TH2D("histo2DRatioEtaPi0TheorySep2760GeV", "histo2DRatioEtaPi0TheorySep2760GeV", 20,0.,maxPtEta2760GeV,1000.,-0.4,1.2);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0TheorySep2760GeV, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.8,0.65, 510, 510);
	histo2DRatioEtaPi0TheorySep2760GeV->GetYaxis()->SetRangeUser(0.,1.02);
	histo2DRatioEtaPi0TheorySep2760GeV->Draw();
	
	graphRatioEtaPi0SysErr2760GeVNewXError->Draw("2same");
	graphRatioEtaPi0StatErr2760GeVNewXError->Draw("psame");
	graphEtaToPi0NLOMuHalf2760GeV->Draw("same,c");
	graphEtaToPi0NLOMuOne2760GeV->Draw("same,c");
	graphEtaToPi0NLOMuTwo2760GeV->Draw("same,c");
	
	histo2DRatioEtaPi0TheorySep2760GeV->Draw("AXIS,same");
	
	TLegend* legendRatioEtaToPi0TheorySep2760GeV = new TLegend(0.5,0.16,0.977,0.34);
	legendRatioEtaToPi0TheorySep2760GeV->SetTextSize(0.036);			
	legendRatioEtaToPi0TheorySep2760GeV->SetFillColor(0);
	legendRatioEtaToPi0TheorySep2760GeV->SetBorderSize(0);
	legendRatioEtaToPi0TheorySep2760GeV->AddEntry(graphRatioEtaPi0ComplErr2760GeVNewXError,Form("ALICE (PCM), #sqrt{#it{s}} = 2.76 TeV "),"fp");
	legendRatioEtaToPi0TheorySep2760GeV->AddEntry(graphEtaToPi0NLOMuHalf2760GeV,"NLO #mu = 0.5 #it{p}_{T}","l");
	legendRatioEtaToPi0TheorySep2760GeV->AddEntry(graphEtaToPi0NLOMuOne2760GeV,"NLO #mu = #it{p}_{T}","l");
	legendRatioEtaToPi0TheorySep2760GeV->AddEntry(graphEtaToPi0NLOMuTwo2760GeV,"NLO #mu = 2 #it{p}_{T}","l");
	legendRatioEtaToPi0TheorySep2760GeV->Draw();
	
	canvasRatioEtaPi0NLO2760GeVSep->Update();
	canvasRatioEtaPi0NLO2760GeVSep->SaveAs(Form("%s/%s_Pi0EtaRatioTheorySep2760GeV_Paper.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	if(!thesis)DrawAliceLogoCombinedPreliminary(pictDrawingCoordinatesPi0Eta[0], pictDrawingCoordinatesPi0Eta[1], pictDrawingCoordinatesPi0Eta[2], pictDrawingCoordinatesPi0Eta[3], pictDrawingCoordinatesPi0Eta[4], pictDrawingCoordinatesPi0Eta[5], pictDrawingCoordinatesPi0Eta[6], pictDrawingCoordinatesPi0Eta[7], pictDrawingCoordinates[8],collisionSystem2760GeV, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900);
	
	canvasRatioEtaPi0NLO2760GeVSep->Update();
	canvasRatioEtaPi0NLO2760GeVSep->SaveAs(Form("%s/%s_Pi0EtaRatioTheorySep2760GeV.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	
	//************************** Ratio + World data ******************************************** 
	cout << "vorm WorldDataFile einlesen" << endl;
	fileWorldDataPi0Eta = new TFile("ExternalInput/WorldDataPi0Eta.root");
	graphDonaldson100GeV = (TGraphErrors*)fileWorldDataPi0Eta->Get("donaldson100GeV");
	DrawGammaSetMarkerTGraphErr(graphDonaldson100GeV, 24, 0.8, kGray+3, kGray+3);
	graphDonaldson200GeV = (TGraphErrors*)fileWorldDataPi0Eta->Get("donaldson200GeV");
	DrawGammaSetMarkerTGraphErr(graphDonaldson200GeV, 25, 0.8, kGray+3, kGray+3);
	graphAntille87pp = (TGraphErrors*)fileWorldDataPi0Eta->Get("Antille24.3GeVpp");
	DrawGammaSetMarkerTGraphErr(graphAntille87pp, 26, 0.8, kGray+3, kGray+3);
	graphAguilar400GeV = (TGraphErrors*)fileWorldDataPi0Eta->Get("Aguilar400GeV");
	DrawGammaSetMarkerTGraphErr(graphAguilar400GeV, 27, 0.8, kGray+3, kGray+3);
	graphKourkou79pp = (TGraphErrors*)fileWorldDataPi0Eta->Get("Kourkou30.6GeVpp");
	DrawGammaSetMarkerTGraphErr(graphKourkou79pp, 28, 0.8, kGreen+2, kGreen+2);
	graphApana530GeV = (TGraphErrors*)fileWorldDataPi0Eta->Get("Apana530GeV");
	DrawGammaSetMarkerTGraphErr(graphApana530GeV, 29, 0.8, kGray+3, kGray+3);
	graphKourkou79pp52 = (TGraphErrors*)fileWorldDataPi0Eta->Get("Kourkou52.7GeVpp");
	DrawGammaSetMarkerTGraphErr(graphKourkou79pp52, 30, 0.8, kGreen+2, kGreen+2);
	graphAkesson53GeVpp = (TGraphErrors*)fileWorldDataPi0Eta->Get("Akesson53GeVpp");
	DrawGammaSetMarkerTGraphErr(graphAkesson53GeVpp, 3, 0.8, kGreen+2, kGreen+2);
	graphKourkou79pp62 = (TGraphErrors*)fileWorldDataPi0Eta->Get("Kourkou79pp62");
	DrawGammaSetMarkerTGraphErr(graphKourkou79pp62, 22, 0.6, kGreen+2, kGreen+2);
	graphPhenix200GeV = (TGraphErrors*)fileWorldDataPi0Eta->Get("Phenix200GeV");
	DrawGammaSetMarkerTGraphErr(graphPhenix200GeV, 5, 0.8, kGreen+2, kGreen+2);
	cout << "WorldDataFile fertig" << endl;
	
	TCanvas* canvasRatioEtaPi03 = new TCanvas("canvasRatioEtaPi03","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEtaPi03,0.09, 0.01, 0.015, 0.115);
	canvasRatioEtaPi03->SetLogy(1);
	
	TH2D *histo2DRatioEtaPi0 = new TH2D("histo2DRatioEtaPi0", "histo2DRatioEtaPi0", 20,0.,maxPtEtaPHOS,1000.,-0.4,6.0);
	SetStyleHistoTH2ForGraphs(histo2DRatioEtaPi0, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.8,0.65, 510, 510);
	histo2DRatioEtaPi0->GetYaxis()->SetRangeUser(0.03,2.8);
	histo2DRatioEtaPi0->Draw();

	graphRatioEtaPi0ComplErr900GeVNewXError->Draw("p,E2same");  
	graphRatioEtaPi0ComplErr2760GeVNewXError->Draw("p,E2same");  
	graphCombinedEtaToPi0NewXError->Draw("p,E2same");  

	graphDonaldson100GeV->Draw("p,same");
	graphDonaldson200GeV->Draw("p,same");
	graphAntille87pp->Draw("p,same");
	graphAguilar400GeV->Draw("p,same");
	graphKourkou79pp->Draw("p,same");
	graphApana530GeV->Draw("p,same");
	graphKourkou79pp52->Draw("p,same");
	graphAkesson53GeVpp->Draw("p,same");
	graphKourkou79pp62->Draw("p,same");
	graphPhenix200GeV->Draw("p,same");
	
	graphCombinedEtaToPi0NoError->Draw("p,same"); 
	graphRatioEtaPi0ComplErr2760GeVNoError->Draw("p,same");  
	graphRatioEtaPi0ComplErr900GeVNoError->Draw("p,same");  
		
	TLegend* legendRatio3 = new TLegend(0.22,0.14,0.96,0.28);
	legendRatio3->SetTextSize(0.027);			
	legendRatio3->SetFillColor(0);
	legendRatio3->SetBorderSize(0);
	legendRatio3->SetNColumns(2);
	legendRatio3->AddEntry(graphCombinedEtaToPi0NewXError,Form("p+p ALICE (#sqrt{#it{s}} = 7 TeV, PLB 717 (2012) 162-172)"),"fp");
	legendRatio3->AddEntry(graphRatioEtaPi0ComplErr2760GeVNewXError,Form("p+p ALICE (#sqrt{#it{s}} = 2.76 TeV)"),"fp");
	legendRatio3->AddEntry(graphRatioEtaPi0ComplErr900GeVNewXError,Form("p+p ALICE (#sqrt{#it{s}} = 0.9 TeV)"),"fp");
	legendRatio3->AddEntry(graphDonaldson100GeV,graphDonaldson100GeV->GetTitle(),"p");
	legendRatio3->AddEntry(graphDonaldson200GeV,graphDonaldson200GeV->GetTitle(),"p");
	legendRatio3->AddEntry(graphAntille87pp,graphAntille87pp->GetTitle(),"p");
	legendRatio3->AddEntry(graphAguilar400GeV,graphAguilar400GeV->GetTitle(),"p");
	legendRatio3->AddEntry(graphKourkou79pp,graphKourkou79pp->GetTitle(),"p");
	legendRatio3->AddEntry(graphApana530GeV,graphApana530GeV->GetTitle(),"p");
	legendRatio3->AddEntry(graphKourkou79pp52,graphKourkou79pp52->GetTitle(),"p");
	legendRatio3->AddEntry(graphAkesson53GeVpp,graphAkesson53GeVpp->GetTitle(),"p");
	legendRatio3->AddEntry(graphKourkou79pp62,graphKourkou79pp62->GetTitle(),"p");
	legendRatio3->AddEntry(graphPhenix200GeV,graphPhenix200GeV->GetTitle(),"p");	
	legendRatio3->Draw();
	
	canvasRatioEtaPi03->Update();
	canvasRatioEtaPi03->SaveAs(Form("%s/%s_Pi0EtaRatioWorld_Paper.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	if(!thesis)DrawAliceLogoCombinedPreliminary(pictDrawingCoordinatesPi0Eta[0], pictDrawingCoordinatesPi0Eta[1], pictDrawingCoordinatesPi0Eta[2], pictDrawingCoordinatesPi0Eta[3], pictDrawingCoordinatesPi0Eta[4], pictDrawingCoordinatesPi0Eta[5], pictDrawingCoordinatesPi0Eta[6], pictDrawingCoordinatesPi0Eta[7], pictDrawingCoordinates[8],collisionSystemCombinedReallyAll, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900);
	
	canvasRatioEtaPi03->Update();
	canvasRatioEtaPi03->SaveAs(Form("%s/%s_Pi0EtaRatioWorld.%s",outputDir.Data(), prefix2.Data(), suffix.Data()));
	
	delete canvasRatioEtaPi03;
	delete legendRatio3;
	cout << "Ratio" << endl;

	
	//**************************************************************************************************
	//********************************* Efficiency Pi0 *************************************************
	//**************************************************************************************************
	TCanvas* canvasEffiSimple = new TCanvas("canvasEffiSimple","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasEffiSimple, 0.08, 0.02, 0.02, 0.09);
	canvasEffiSimple->SetLogy(1);	
	
	DrawAutoGammaMesonHistos( histoTrueEffPtPi0, 
								"", "#it{p}_{T} (GeV/#it{c})", "Efficiency x Acceptance", 
								kTRUE, 2.5, 3e-7, kFALSE,
								kFALSE, 0., 0.7, 
								kFALSE, 0., 10.);
	histoTrueEffPtPi0->GetYaxis()->SetTitleOffset(1.);
	histoTrueEffPtPi0->GetXaxis()->SetTitleOffset(1.);
	DrawGammaSetMarker(histoTrueEffPtPi0, markerStyleConv , markerSizeInvYield , kBlack, kBlack);										 
	histoTrueEffPtPi0->DrawCopy("e1"); 	
	
	DrawGammaSetMarker(histoTrueEffPtPi0900GeV, markerStyleConv, markerSizeInvYield, kBlue, kBlue);										 
	histoTrueEffPtPi0900GeV->DrawCopy("e1,same"); 	
	
	TLegend* legendEffi = new TLegend(0.65,0.12,0.95,0.22);
	legendEffi->SetTextSize(0.03);		
	legendEffi->SetFillColor(0);
	legendEffi->AddEntry(histoTrueEffPtPi0,Form("#pi^{0}, #sqrt{#it{s}} = 7 TeV"));
	legendEffi->AddEntry(histoTrueEffPtPi0900GeV,Form("#pi^{0}, #sqrt{#it{s}} = 0.9 TeV"));
	legendEffi->Draw();
	
	if(!thesis)DrawAliceLogoPi0Performance(pictDrawingCoordinatesAccEff[0], pictDrawingCoordinatesAccEff[1], pictDrawingCoordinatesAccEff[2], pictDrawingCoordinatesAccEff[3], pictDrawingCoordinatesAccEff[4], pictDrawingCoordinatesAccEff[5], pictDrawingCoordinatesAccEff[6], pictDrawingCoordinatesAccEff[7], pictDrawingCoordinates[8],collisionSystemCombined, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date);
	canvasEffiSimple->Update();
	
	canvasEffiSimple->SaveAs(Form("%s/EffXAccConvDiffEnergies.%s",outputDir.Data(),suffix.Data()));
	
	delete canvasEffiSimple;

	
	
	
	
	//************************************** Mass and FWHM combined *****************************************************************
	TCanvas* canvasMassPlusFWHM = new TCanvas("canvasMassPlusFWHM","",200,10,2700,1800);  // gives the page size
	DrawGammaCanvasSettings( canvasMassPlusFWHM, 0.07, 0.02, 0.02, 0.10);	
	canvasMassPlusFWHM->cd();
	canvasMassPlusFWHM->SetLogx();
	
	TH2D *histo2DAllPi0MassAndWidthCombined;
	histo2DAllPi0MassAndWidthCombined = new TH2D("histo2DAllPi0MassAndWidthCombined", "histo2DAllPi0MassAndWidthCombined", 20,0.2,maxPtPi0 ,1000.,-30,40);
	histo2DAllPi0MassAndWidthCombined->GetYaxis()->SetRangeUser(-27,20);
	histo2DAllPi0MassAndWidthCombined->SetYTitle("m- m_{PDG}- m_{o}, #sigma (MeV/#it{c}^{2})");
	histo2DAllPi0MassAndWidthCombined->SetXTitle("#it{p}_{T} (GeV/#it{c})");
	histo2DAllPi0MassAndWidthCombined->GetXaxis()->SetNdivisions(515,kTRUE);
	histo2DAllPi0MassAndWidthCombined->GetYaxis()->SetNdivisions(515,kTRUE);
	histo2DAllPi0MassAndWidthCombined->GetYaxis()->SetLabelSize(0.035);
	histo2DAllPi0MassAndWidthCombined->GetYaxis()->SetTitleSize(0.05);	
	histo2DAllPi0MassAndWidthCombined->GetYaxis()->SetDecimals();
	histo2DAllPi0MassAndWidthCombined->GetYaxis()->SetTitleOffset(0.63);
	histo2DAllPi0MassAndWidthCombined->GetXaxis()->SetTitleSize(0.05);	
	histo2DAllPi0MassAndWidthCombined->GetXaxis()->SetLabelSize(0.035);
	histo2DAllPi0MassAndWidthCombined->SetTitle("");
	histo2DAllPi0MassAndWidthCombined->GetXaxis()->SetTitleOffset(0.9);
	histo2DAllPi0MassAndWidthCombined->DrawCopy();
	
	DrawGammaSetMarker(histoFWHMMesonPi0PHOS, markerStylePHOS, markerSizeMass, colorPHOS, colorPHOS);
	histoFWHMMesonPi0PHOS->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoTrueFWHMMesonPi0PHOS, markerStylePHOSMC, markerSizeMass, colorPHOSMC , colorPHOSMC);
	histoTrueFWHMMesonPi0PHOS->DrawCopy("same,e1,p"); 
	
	DrawGammaSetMarker(histoFWHMMesonPi0EMCAL, markerStylePHOS, markerSizeMass, colorEMCAL, colorEMCAL);
	histoFWHMMesonPi0EMCAL->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoTrueFWHMMesonPi0EMCAL, markerStylePHOSMC, markerSizeMass, colorEMCALMC, colorEMCALMC);
	histoTrueFWHMMesonPi0EMCAL->DrawCopy("same,e1,p"); 
	
	DrawGammaSetMarker(histoFWHMMesonPi0MeV, markerStylePHOS, markerSizeMass, colorConv, colorConv);
	histoFWHMMesonPi0MeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoTrueFWHMMesonPi0MeV, markerStylePHOSMC, markerSizeMass, colorConvMC, colorConvMC);
	histoTrueFWHMMesonPi0MeV->DrawCopy("same,e1,p"); 
	
	DrawGammaSetMarker(histoMassMesonPi0PHOSMinusExp, markerStyleConv, markerSizeMass, colorPHOS,colorPHOS);
	histoMassMesonPi0PHOSMinusExp->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoTrueMassMesonPi0PHOSMinusExp, markerStyleConvMC, markerSizeMass, colorPHOSMC, colorPHOSMC);
	histoTrueMassMesonPi0PHOSMinusExp->DrawCopy("same,e1,p"); 
	
	DrawGammaSetMarker(histoMassMesonPi0EMCALMinusExp, markerStyleConv, markerSizeMass, colorEMCAL,colorEMCAL);
	histoMassMesonPi0EMCALMinusExp->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoTrueMassMesonPi0EMCALMinusExp, markerStyleConvMC, markerSizeMass, colorEMCALMC, colorEMCALMC);
	histoTrueMassMesonPi0EMCALMinusExp->DrawCopy("same,e1,p"); 
	
	DrawGammaSetMarker(histoMassMesonPi0MinusExp, markerStyleConv, markerSizeMass, colorConv, colorConv);					 
	histoMassMesonPi0MinusExp->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoTrueMassMesonPi0MinusExp, markerStyleConvMC , markerSizeMass, colorConvMC, colorConvMC);					 
	histoTrueMassMesonPi0MinusExp->DrawCopy("same,e1,p"); 
	
	
	DrawGammaLines(0., maxPtPi0 ,-5., -5.,0.1,colorConv);
	DrawGammaLines(0., maxPtPi0 ,-10., -10.,0.1,colorPHOS );
	DrawGammaLines(0., maxPtPi0,-15.,-15.,0.1, colorEMCAL);
	
	//********************************** Defintion of the Legend **************************************************	
	Double_t columnsLegendMass[7] 	= {0.,0.18,0.29,0.4,0.59,0.70,0.83};
	Double_t rowsLegendMass[5] 		= {0.80,0.60,0.40,0.20,0.0};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnMass	= 0.19;
	Size_t textSizeTopRowMass	= 0.19; 
	Size_t textSizeSecondRowMass 	= 0.17;
	//******************* Offsets ***********************
	Double_t offsetMarkerXMass	= 0.02;
	Double_t offsetMarkerYMass	= 0.05;
	//****************** Scale factors ******************
	Double_t scaleMarkerMass		= 1.5;
	
	TPad* padMassLegend = new TPad("padMassLegend", "", 0.15, 0.12, 0.95, 0.27,-1, -1, -2);
	DrawGammaPadSettings( padMassLegend, 0., 0., 0., 0.);
	padMassLegend->Draw();
	padMassLegend->cd();
	
	//****************** first Column **************************************************
	TLatex *textMassCTS = new TLatex(columnsLegendMass[0],rowsLegendMass[2],"PCM");
	SetStyleTLatex( textMassCTS, textSizeLeftColumnMass,4);
	textMassCTS->Draw();
	TLatex *textMassPHOS = new TLatex(columnsLegendMass[0],rowsLegendMass[3],"PHOS");
	SetStyleTLatex( textMassPHOS, textSizeLeftColumnMass,4);
	textMassPHOS->Draw();
	TLatex *textMassEMCAL = new TLatex(columnsLegendMass[0],rowsLegendMass[4],"EMCAL");
	SetStyleTLatex( textMassEMCAL, textSizeLeftColumnMass,4);
	textMassEMCAL->Draw();
	
	//****************** second Column *************************************************
	TLatex *textMass = new TLatex(columnsLegendMass[1],rowsLegendMass[0],"Mass: m- m_{PDG}- m_{o}");
	SetStyleTLatex( textMass, textSizeTopRowMass,4);
	textMass->Draw();
	TLatex *textMassData = new TLatex(columnsLegendMass[1],rowsLegendMass[1] ,"Data");
	SetStyleTLatex( textMassData, textSizeSecondRowMass,4);
	textMassData->Draw();
	TLatex *textMassMC = new TLatex(columnsLegendMass[2] ,rowsLegendMass[1],"MC");
	SetStyleTLatex( textMassMC, textSizeSecondRowMass,4);
	textMassMC->Draw();
	TLatex *textMassOffset = new TLatex(columnsLegendMass[3] ,rowsLegendMass[1],"m_{0} (MeV/#it{c}^{2})");
	SetStyleTLatex( textMassOffset, textSizeSecondRowMass,4);
	textMassOffset->Draw();
	
	TMarker* markerCTSPi07TeVMass = CreateMarkerFromHisto(histoMassMesonPi0MinusExp,columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass ,scaleMarkerMass);
	markerCTSPi07TeVMass->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
	TMarker* markerPHOSPi07TeVMass = CreateMarkerFromHisto(histoMassMesonPi0PHOSMinusExp,columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass ,scaleMarkerMass);
	markerPHOSPi07TeVMass->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	TMarker* markerEMCALPi07TeVMass = CreateMarkerFromHisto(histoMassMesonPi0EMCALMinusExp,columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass ,scaleMarkerMass);
	markerEMCALPi07TeVMass->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);
	
	TMarker* markerCTSPi07TeVMassMC = CreateMarkerFromHisto(histoTrueMassMesonPi0MinusExp,columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass ,scaleMarkerMass);
	markerCTSPi07TeVMassMC->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
	TMarker* markerPHOSPi07TeVMassMC = CreateMarkerFromHisto(histoTrueMassMesonPi0PHOSMinusExp,columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass ,scaleMarkerMass);
	markerPHOSPi07TeVMassMC->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	TMarker* markerEMCALPi07TeVMassMC = CreateMarkerFromHisto(histoTrueMassMesonPi0EMCALMinusExp,columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass ,scaleMarkerMass);
	markerEMCALPi07TeVMassMC->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);

	TLatex *textOffsetConv = new TLatex(columnsLegendMass[3]+offsetMarkerXMass,rowsLegendMass[2] ,"5");
	SetStyleTLatex( textOffsetConv, textSizeSecondRowMass,4);
	textOffsetConv->Draw();
	TLatex *textOffsetPHOS = new TLatex(columnsLegendMass[3]+offsetMarkerXMass ,rowsLegendMass[3],"10");
	SetStyleTLatex( textOffsetPHOS, textSizeSecondRowMass,4);
	textOffsetPHOS->Draw();
	TLatex *textOffsetEMCAL = new TLatex(columnsLegendMass[3]+offsetMarkerXMass ,rowsLegendMass[4],"15");
	SetStyleTLatex( textOffsetEMCAL, textSizeSecondRowMass,4);
	textOffsetEMCAL->Draw();
	
	//***************** third Column ***************************************************
	TLatex *textWidth = new TLatex(columnsLegendMass[4],rowsLegendMass[0],"Width");
	SetStyleTLatex( textWidth, textSizeTopRowMass,4);
	textWidth->Draw();
	TLatex *textWidthData = new TLatex(columnsLegendMass[4],rowsLegendMass[1] ,"Data");
	SetStyleTLatex( textWidthData, textSizeSecondRowMass,4);
	textWidthData->Draw();
	TLatex *textWidthMC = new TLatex(columnsLegendMass[5] ,rowsLegendMass[1],"MC");
	SetStyleTLatex( textWidthMC, textSizeSecondRowMass,4);
	textWidthMC->Draw();

	TMarker* markerCTSPi07TeVWidth = CreateMarkerFromHisto(histoFWHMMesonPi0MeV,columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass ,scaleMarkerMass);
	markerCTSPi07TeVWidth->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
	TMarker* markerPHOSPi07TeVWidth = CreateMarkerFromHisto(histoFWHMMesonPi0PHOS,columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass ,scaleMarkerMass);
	markerPHOSPi07TeVWidth->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	TMarker* markerEMCALPi07TeVWidth = CreateMarkerFromHisto(histoFWHMMesonPi0EMCAL,columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass ,scaleMarkerMass);
	markerEMCALPi07TeVWidth->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);
	
	TMarker* markerCTSPi07TeVWidthMC = CreateMarkerFromHisto(histoTrueFWHMMesonPi0MeV,columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass ,scaleMarkerMass);
	markerCTSPi07TeVWidthMC->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
	TMarker* markerPHOSPi07TeVWidthMC = CreateMarkerFromHisto(histoTrueFWHMMesonPi0PHOS,columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass ,scaleMarkerMass);
	markerPHOSPi07TeVWidthMC->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	TMarker* markerEMCALPi07TeVWidthMC = CreateMarkerFromHisto(histoTrueFWHMMesonPi0EMCAL,columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass ,scaleMarkerMass);
	markerEMCALPi07TeVWidthMC->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);
	
	TLatex *textWidthConv = new TLatex(columnsLegendMass[6],rowsLegendMass[2] ,"FWHM/2.36");
	SetStyleTLatex( textWidthConv, textSizeSecondRowMass,4);
	textWidthConv->Draw();
	TLatex *textWidthPHOS = new TLatex(columnsLegendMass[6] ,rowsLegendMass[3],"#sigma");
	SetStyleTLatex( textWidthPHOS, textSizeSecondRowMass,4);
	textWidthPHOS->Draw();
	TLatex *textWidthEMCAL = new TLatex(columnsLegendMass[6] ,rowsLegendMass[4],"#sigma");
	SetStyleTLatex( textWidthEMCAL, textSizeSecondRowMass,4);
	textWidthEMCAL->Draw();
	
	canvasMassPlusFWHM->Update();
	canvasMassPlusFWHM->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurements_7TeVLogX_Paper.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));
	
	canvasMassPlusFWHM->SetLogx(0);
	canvasMassPlusFWHM->Update();
	canvasMassPlusFWHM->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurements_7TeVLinX_Paper.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));
	

	TCanvas* canvasMassPlusFWHMLogo = new TCanvas("canvasMassPlusFWHMLogo","",200,10,2700,1800);  // gives the page size
	DrawGammaCanvasSettings( canvasMassPlusFWHMLogo, 0.07, 0.02, 0.02, 0.10);	
	canvasMassPlusFWHMLogo->cd();
	canvasMassPlusFWHMLogo->SetLogx();
	TH2D *histo2DAllPi0MassAndWidthCombinedLogo;
	histo2DAllPi0MassAndWidthCombinedLogo = new TH2D("histo2DAllPi0MassAndWidthCombinedLogo", "histo2DAllPi0MassAndWidthCombinedLogo", 20,0.2,maxPtPi0 ,1000.,-30,40);
	histo2DAllPi0MassAndWidthCombinedLogo->GetYaxis()->SetRangeUser(-27.5,25);
	histo2DAllPi0MassAndWidthCombinedLogo->SetYTitle("(MeV/#it{c}^{2})");
	histo2DAllPi0MassAndWidthCombinedLogo->SetXTitle("#it{p}_{T} (GeV/#it{c})");
   histo2DAllPi0MassAndWidthCombinedLogo->GetXaxis()->SetNdivisions(515,kTRUE);
	histo2DAllPi0MassAndWidthCombinedLogo->GetYaxis()->SetNdivisions(515,kTRUE);
	histo2DAllPi0MassAndWidthCombinedLogo->GetYaxis()->SetLabelSize(0.035);
	histo2DAllPi0MassAndWidthCombinedLogo->GetYaxis()->SetTitleSize(0.05);	
	histo2DAllPi0MassAndWidthCombinedLogo->GetYaxis()->SetDecimals();
	histo2DAllPi0MassAndWidthCombinedLogo->GetYaxis()->SetTitleOffset(0.63);
	histo2DAllPi0MassAndWidthCombinedLogo->GetXaxis()->SetTitleSize(0.05);	
	histo2DAllPi0MassAndWidthCombinedLogo->GetXaxis()->SetLabelSize(0.035);
	histo2DAllPi0MassAndWidthCombinedLogo->SetTitle("");
	histo2DAllPi0MassAndWidthCombinedLogo->GetXaxis()->SetTitleOffset(0.9);
	histo2DAllPi0MassAndWidthCombinedLogo->DrawCopy();
	
	histoFWHMMesonPi0PHOS->DrawCopy("same,e1,p"); 
	histoTrueFWHMMesonPi0PHOS->DrawCopy("same,e1,p"); 
	histoFWHMMesonPi0EMCAL->DrawCopy("same,e1,p"); 
	histoTrueFWHMMesonPi0EMCAL->DrawCopy("same,e1,p"); 	
	histoFWHMMesonPi0MeV->DrawCopy("same,e1,p"); 
	histoTrueFWHMMesonPi0MeV->DrawCopy("same,e1,p"); 
	histoMassMesonPi0PHOSMinusExp->DrawCopy("same,e1,p"); 
	histoTrueMassMesonPi0PHOSMinusExp->DrawCopy("same,e1,p"); 
	histoMassMesonPi0EMCALMinusExp->DrawCopy("same,e1,p"); 
	histoTrueMassMesonPi0EMCALMinusExp->DrawCopy("same,e1,p"); 
	histoMassMesonPi0MinusExp->DrawCopy("same,e1,p"); 
	histoTrueMassMesonPi0MinusExp->DrawCopy("same,e1,p"); 
	DrawGammaLines(0., maxPtPi0 ,-5., -5.,0.1,colorConv);
	DrawGammaLines(0., maxPtPi0,-10., -10.,0.1,colorPHOS );
	DrawGammaLines(0., maxPtPi0,-15.,-15.,0.1, colorEMCAL);
	
	padMassLegend->Draw();

	canvasMassPlusFWHMLogo->Update();
	canvasMassPlusFWHMLogo->cd();
	
	if(!thesis)DrawAliceLogoPi0WithPHOSPerformance(pictDrawingCoordinatesMassFWHM2[0], pictDrawingCoordinatesMassFWHM2[1], pictDrawingCoordinatesMassFWHM2[2], pictDrawingCoordinatesMassFWHM2[3], pictDrawingCoordinatesMassFWHM2[4], pictDrawingCoordinatesMassFWHM2[5], pictDrawingCoordinatesMassFWHM2[6], pictDrawingCoordinatesMassFWHM2[7], pictDrawingCoordinates[8], collisionSystem7TeV, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],2700,1800,date);
	canvasMassPlusFWHMLogo->Update();
	canvasMassPlusFWHMLogo->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurements_7TeVLogX.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));
	
	
	canvasMassPlusFWHMLogo->SetLogx(0);
	canvasMassPlusFWHMLogo->Update();
	canvasMassPlusFWHMLogo->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurements_7TeVLinX.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));

	//************************************** Mass and FWHM combined devided pads*****************************************************************
	TCanvas* canvasMassPlusFWHMDifPads = new TCanvas("canvasMassPlusFWHMDifPads","",200,10,2700,2500);  // gives the page size
	DrawGammaCanvasSettings( canvasMassPlusFWHMDifPads, 0.00, 0.00, 0.00, 0.0);	
	canvasMassPlusFWHMDifPads->cd();
	
	TPad* padFWHM7TeV = new TPad("padFWHM7TeV", "", 0., 0.56, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padFWHM7TeV, 0.11, 0.005, 0.02, 0.);
	padFWHM7TeV->SetLogx();
	padFWHM7TeV->Draw();
	
	TPad* padMass7TeV = new TPad("padMass7TeV", "", 0., 0., 1., 0.56,-1, -1, -2);
	DrawGammaPadSettings( padMass7TeV, 0.11, 0.005, 0., 0.15);
	padMass7TeV->Draw();
	
	TPad* padMassLegend1 = new TPad("padMassLegend1", "", 0.15, 0.41, 0.42, 0.53,-1, -1, -2);
	DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
	padMassLegend1->Draw();
		
	padFWHM7TeV->cd();
	padFWHM7TeV->SetLogx();
	
	TH2D *histo2DAllPi0FWHM7TeV;
	histo2DAllPi0FWHM7TeV = new TH2D("histo2DAllPi0FWHM7TeV", "histo2DAllPi0FWHM7TeV", 20,0.25,35. ,1000.,-30,40);
	histo2DAllPi0FWHM7TeV->GetYaxis()->SetRangeUser(-1.,25.5);
	histo2DAllPi0FWHM7TeV->SetYTitle("peak width (MeV/#it{c}^{2})");
	histo2DAllPi0FWHM7TeV->SetXTitle("#it{p}_{T} (GeV/#it{c})");
   histo2DAllPi0FWHM7TeV->SetLineWidth(2);
   histo2DAllPi0FWHM7TeV->SetTitleFont(42,"Y");
   histo2DAllPi0FWHM7TeV->SetLabelFont(42,"Y");
	histo2DAllPi0FWHM7TeV->GetXaxis()->SetNdivisions(515,kTRUE);
	histo2DAllPi0FWHM7TeV->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DAllPi0FWHM7TeV->GetYaxis()->SetLabelSize(0.074);
   histo2DAllPi0FWHM7TeV->GetYaxis()->SetTitleSize(0.085);  
   histo2DAllPi0FWHM7TeV->GetYaxis()->SetTickLength(0.015);  
   histo2DAllPi0FWHM7TeV->GetYaxis()->SetDecimals();
   histo2DAllPi0FWHM7TeV->GetYaxis()->SetLabelOffset(0.01);
   histo2DAllPi0FWHM7TeV->GetYaxis()->SetTitleOffset(0.67);
   histo2DAllPi0FWHM7TeV->GetXaxis()->SetTickLength(0.04);  
   histo2DAllPi0FWHM7TeV->SetTitle("");
   histo2DAllPi0FWHM7TeV->DrawCopy();

   
	DrawGammaSetMarker(histoFWHMMesonPi0PHOS, markerStylePHOS, markerSizeMass*1.5, colorPHOSMass, colorPHOSMass);
	histoFWHMMesonPi0PHOS->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoTrueFWHMMesonPi0PHOS, markerStylePHOSMC, markerSizeMass*1.5, colorPHOSMCMass , colorPHOSMCMass);
	histoTrueFWHMMesonPi0PHOS->DrawCopy("same,e1,p,x0"); 
	
   DrawGammaSetMarker(histoFWHMMesonPi0EMCAL, markerStyleEMCAL, markerSizeMass*1.5*1.5, colorEMCALMass, colorEMCALMass);
   histoFWHMMesonPi0EMCAL->DrawCopy("same,e1,p,x0"); 
   DrawGammaSetMarker(histoTrueFWHMMesonPi0EMCAL, markerStyleEMCALMC, markerSizeMass*1.5*1.5, colorEMCALMCMass , colorEMCALMCMass);
   histoTrueFWHMMesonPi0EMCAL->DrawCopy("same,e1,p,x0"); 
   
   
	DrawGammaSetMarker(histoFWHMMesonPi0MeV, markerStyleConv, markerSizeMass*1.5, colorConv, colorConv);
	histoFWHMMesonPi0MeV->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoTrueFWHMMesonPi0MeV, markerStyleConvMC, markerSizeMass*1.5, colorConvMC, colorConvMC);
   
	histoTrueFWHMMesonPi0MeV->DrawCopy("same,e1,p,x0"); 
	
	TLatex *labelLegendAMass = new TLatex(0.15,0.04,"a)");
	SetStyleTLatex( labelLegendAMass, 0.07,4);
   labelLegendAMass->SetTextFont(42);
	labelLegendAMass->Draw();
	
   TLatex *labelMassEnergy = new TLatex(0.8,0.9,"pp #sqrt{#it{s}} = 7 TeV");
   SetStyleTLatex( labelMassEnergy, 0.07,4);
   labelMassEnergy->SetTextFont(42);
   labelMassEnergy->Draw();

   //********************************** Defintion of the Legend **************************************************	
	Double_t columnsLegendFWHM[4] 	= {0.15,0.4,0.48,0.39};
   Double_t columnsLegendFWHMAbs[4]    = {4,1.5,2.3,12};
// 	Double_t rowsLegendFWHM[3] 		= {0.66,0.33,0.0};
   Double_t rowsLegendFWHM[4]       = {0.9,0.82,0.76,0.7};
   Double_t rowsLegendFWHMAbs[4]       = {0.2,22,20.,18.};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnFWHM	= 0.07;
	Size_t textSizeTopRowFWHM	= 0.07; 
	//******************* Offsets ***********************
	Double_t offsetMarkerXFWHM	= 0.07;
	Double_t offsetMarkerYFWHM	= 0.07;
	//****************** Scale factors ******************
	Double_t scaleMarkerFWHM		= 1.2;
	
// 	padFWHMLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textFWHMCTS = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM[1],"PCM (FWHM/2.36)");
	SetStyleTLatex( textFWHMCTS, textSizeLeftColumnFWHM,4);
   textFWHMCTS->SetTextFont(42);
	textFWHMCTS->Draw();
	TLatex *textFWHMPHOS = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM[2],"PHOS (#sigma)");
	SetStyleTLatex( textFWHMPHOS, textSizeLeftColumnFWHM,4);
   textFWHMPHOS->SetTextFont(42);
	textFWHMPHOS->Draw();
   TLatex *textFWHMEMCAL = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM[3],"EMCal (#sigma)");
   SetStyleTLatex( textFWHMEMCAL, textSizeLeftColumnFWHM,4);
   textFWHMEMCAL->SetTextFont(42);
   textFWHMEMCAL->Draw();
	
	//****************** second Column *************************************************
	TLatex *textFWHMData2 = new TLatex(columnsLegendFWHM[1],rowsLegendFWHM[0] ,"data");
	SetStyleTLatex( textFWHMData2, textSizeTopRowFWHM ,4);
   textFWHMData2->SetTextFont(42);
	textFWHMData2->Draw();
	TLatex *textFWHMMC2 = new TLatex(columnsLegendFWHM[2] ,rowsLegendFWHM[0],"MC");
	SetStyleTLatex( textFWHMMC2, textSizeTopRowFWHM,4);
   textFWHMMC2->SetTextFont(42);
	textFWHMMC2->Draw();
	
	TMarker* markerCTSPi07TeVFWHM = CreateMarkerFromHisto(histoFWHMMesonPi0MeV,columnsLegendFWHMAbs[1] ,rowsLegendFWHMAbs[1] ,scaleMarkerFWHM);
	markerCTSPi07TeVFWHM->DrawMarker(columnsLegendFWHMAbs[1] ,rowsLegendFWHMAbs[1]);
	TMarker* markerPHOSPi07TeVFWHM = CreateMarkerFromHisto(histoFWHMMesonPi0PHOS,columnsLegendFWHMAbs[1],rowsLegendFWHMAbs[2],scaleMarkerFWHM);
	markerPHOSPi07TeVFWHM->DrawMarker(columnsLegendFWHMAbs[1] ,rowsLegendFWHMAbs[2]);
	TMarker* markerEMCALPi07TeVFWHM = CreateMarkerFromHisto(histoFWHMMesonPi0EMCAL,columnsLegendFWHMAbs[1],rowsLegendFWHMAbs[3],scaleMarkerFWHM);
   markerEMCALPi07TeVFWHM->DrawMarker(columnsLegendFWHMAbs[1] ,rowsLegendFWHMAbs[3]);
   
	TMarker* markerCTSPi07TeVFWHMMC = CreateMarkerFromHisto(histoTrueFWHMMesonPi0MeV,columnsLegendFWHMAbs[2],rowsLegendFWHMAbs[1],scaleMarkerFWHM);
	markerCTSPi07TeVFWHMMC->DrawMarker(columnsLegendFWHMAbs[2] ,rowsLegendFWHMAbs[1]);
	TMarker* markerPHOSPi07TeVFWHMMC = CreateMarkerFromHisto(histoTrueFWHMMesonPi0PHOS,columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
	markerPHOSPi07TeVFWHMMC->DrawMarker(columnsLegendFWHMAbs[2],rowsLegendFWHMAbs[2]);
	TMarker* markerEMCALPi07TeVFWHMMC = CreateMarkerFromHisto(histoTrueFWHMMesonPi0EMCAL,columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[3]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
   markerEMCALPi07TeVFWHMMC->DrawMarker(columnsLegendFWHMAbs[2],rowsLegendFWHMAbs[3]);
   
	padMass7TeV->cd();
	padMass7TeV->SetLogx();
	
	TH2D *histo2DAllPi0Mass7TeV;
	histo2DAllPi0Mass7TeV = new TH2D("histo2DAllPi0Mass7TeV", "histo2DAllPi0Mass7TeV", 20,0.25,35. ,1000.,125.,170);
	histo2DAllPi0Mass7TeV->GetYaxis()->SetRangeUser(125.,156);
	histo2DAllPi0Mass7TeV->SetYTitle("peak position (MeV/#it{c}^{2})");
	histo2DAllPi0Mass7TeV->SetXTitle("#it{p}_{T} (GeV/#it{c})");
   histo2DAllPi0Mass7TeV->SetTitleFont(42,"X");
   histo2DAllPi0Mass7TeV->SetTitleFont(42,"Y");
   histo2DAllPi0Mass7TeV->SetLabelFont(42,"X");
   histo2DAllPi0Mass7TeV->SetLabelFont(42,"Y");
   histo2DAllPi0Mass7TeV->SetLineWidth(2);
	histo2DAllPi0Mass7TeV->GetXaxis()->SetNdivisions(515,kTRUE);
   histo2DAllPi0Mass7TeV->GetYaxis()->SetNdivisions(510,kTRUE);
   histo2DAllPi0Mass7TeV->GetYaxis()->SetLabelSize(0.063);
   histo2DAllPi0Mass7TeV->GetYaxis()->SetTitleSize(0.068);   
   histo2DAllPi0Mass7TeV->GetYaxis()->SetDecimals();
   histo2DAllPi0Mass7TeV->GetYaxis()->SetLabelOffset(0.01);
   histo2DAllPi0Mass7TeV->GetYaxis()->SetTitleOffset(0.82);
   histo2DAllPi0Mass7TeV->GetXaxis()->SetTitleSize(0.068);   
   histo2DAllPi0Mass7TeV->GetXaxis()->SetLabelSize(0.063);
   histo2DAllPi0Mass7TeV->GetYaxis()->SetTickLength(0.02);  
   histo2DAllPi0Mass7TeV->GetXaxis()->SetTickLength(0.04);  
   histo2DAllPi0Mass7TeV->SetTitle("");
   histo2DAllPi0Mass7TeV->GetXaxis()->SetLabelOffset(-0.015);
   histo2DAllPi0Mass7TeV->GetXaxis()->SetTitleOffset(0.85);
   histo2DAllPi0Mass7TeV->DrawCopy();

	DrawGammaSetMarker(histoMassMesonPi0PHOS , markerStylePHOS, markerSizeMass*1.5, colorPHOSMass,colorPHOSMass);
	histoMassMesonPi0PHOS->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoTrueMassMesonPi0PHOS , markerStylePHOSMC, markerSizeMass*1.5, colorPHOSMCMass, colorPHOSMCMass);
	histoTrueMassMesonPi0PHOS->DrawCopy("same,e1,p,x0"); 

   DrawGammaSetMarker(histoMassMesonPi0EMCAL, markerStyleEMCAL, markerSizeMass*1.5*1.5, colorEMCALMass, colorEMCALMass);
   histoMassMesonPi0EMCAL->DrawCopy("same,e1,p,x0"); 
   DrawGammaSetMarker(histoTrueMassMesonPi0EMCAL, markerStyleEMCALMC, markerSizeMass*1.5*1.5, colorEMCALMCMass , colorEMCALMCMass);
   histoTrueMassMesonPi0EMCAL->DrawCopy("same,e1,p,x0"); 

   
	DrawGammaSetMarker(histoMassMesonPi0 , markerStyleConv, markerSizeMass*1.5, colorConv, colorConv);					 
	histoMassMesonPi0->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoTrueMassMesonPi0 , markerStyleConvMC , markerSizeMass*1.5, colorConvMC, colorConvMC);					 
	histoTrueMassMesonPi0->DrawCopy("same,e1,p,x0"); 
	
	DrawGammaLines(0.25, 35. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1,colorConv);

	TLatex *labelLegendBMass = new TLatex(0.15,0.18,"b)");
	SetStyleTLatex( labelLegendBMass, 0.06,4);
   labelLegendBMass->SetTextFont(42);
	labelLegendBMass->Draw();
	
	//********************************** Defintion of the Legend **************************************************	
	Double_t columnsLegendMass2[3] 	= {0.,0.45,0.75};
	Double_t rowsLegendMass2[4] 		= {0.75,0.5,0.25,0.01};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnMass2	= 0.26;
	Size_t textSizeSecondRowMass2 	= 0.25;
	//******************* Offsets ***********************
	Double_t offsetMarkerXMass2	= 0.1;
	Double_t offsetMarkerYMass2	= 0.1;
	//****************** Scale factors ******************
	Double_t scaleMarkerMass2		= 1.2;
	
	padMassLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textMassCTS2 = new TLatex(columnsLegendMass2[0],rowsLegendMass2[1],"PCM");
	SetStyleTLatex( textMassCTS2, textSizeLeftColumnMass2,4);
   textMassCTS2->SetTextFont(42);
	textMassCTS2->Draw();
	TLatex *textMassPHOS2 = new TLatex(columnsLegendMass2[0],rowsLegendMass2[2],"PHOS");
	SetStyleTLatex( textMassPHOS2, textSizeLeftColumnMass2,4);
   textMassPHOS2->SetTextFont(42);
	textMassPHOS2->Draw();
	TLatex *textMassEMCAL2 = new TLatex(columnsLegendMass2[0],rowsLegendMass2[3],"EMCal");
   SetStyleTLatex( textMassEMCAL2, textSizeLeftColumnMass2,4);
   textMassEMCAL2->SetTextFont(42);
   textMassEMCAL2->Draw();
   
	//****************** second Column *************************************************
	TLatex *textMassData2 = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"data");
	SetStyleTLatex( textMassData2, textSizeSecondRowMass2,4);
   textMassData2->SetTextFont(42);
	textMassData2->Draw();
	TLatex *textMassMC2 = new TLatex(columnsLegendMass2[2] ,rowsLegendMass2[0],"MC");
	SetStyleTLatex( textMassMC2, textSizeSecondRowMass2,4);
   textMassMC2->SetTextFont(42);
	textMassMC2->Draw();
	
	TMarker* markerCTSPi07TeVMass2 = CreateMarkerFromHisto(histoMassMesonPi0,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerCTSPi07TeVMass2->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
	TMarker* markerPHOSPi07TeVMass2 = CreateMarkerFromHisto(histoMassMesonPi0PHOS,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerPHOSPi07TeVMass2->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
	TMarker* markerEMCALPi07TeVMass2 = CreateMarkerFromHisto(histoMassMesonPi0EMCAL,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2 ,scaleMarkerMass2);
   markerEMCALPi07TeVMass2->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
   
	TMarker* markerCTSPi07TeVMassMC2 = CreateMarkerFromHisto(histoTrueMassMesonPi0,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerCTSPi07TeVMassMC2->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
	TMarker* markerPHOSPi07TeVMassMC2 = CreateMarkerFromHisto(histoTrueMassMesonPi0PHOS,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerPHOSPi07TeVMassMC2->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
	TMarker* markerEMCALPi07TeVMassMC2 = CreateMarkerFromHisto(histoTrueMassMesonPi0EMCAL,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2 ,scaleMarkerMass2);
   markerEMCALPi07TeVMassMC2->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
   
	canvasMassPlusFWHMDifPads->Update();
	canvasMassPlusFWHMDifPads->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurementsDP_7TeVLogX_Paper.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));
	canvasMassPlusFWHMDifPads->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurementsDP_7TeVLogX_Paper.C",outputDir.Data(), prefix2.Data()));
	canvasMassPlusFWHMDifPads->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSOnlyPerformance(pictDrawingCoordinatesMassFWHM3[0], pictDrawingCoordinatesMassFWHM3[1], pictDrawingCoordinatesMassFWHM3[2], pictDrawingCoordinatesMassFWHM3[3], pictDrawingCoordinatesMassFWHM3[4], pictDrawingCoordinatesMassFWHM3[5], pictDrawingCoordinatesMassFWHM3[6], pictDrawingCoordinatesMassFWHM3[7], pictDrawingCoordinates[8], collisionSystem7TeV, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],2700,2300,date);
	canvasMassPlusFWHMDifPads->Update();
	canvasMassPlusFWHMDifPads->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurementsDP_7TeVLogX.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));
	

	//************************************** Mass and FWHM combined *****************************************************************
	TCanvas* canvasMassPlusFWHM900GeV = new TCanvas("canvasMassPlusFWHM900GeV","",200,10,2700,1800);  // gives the page size
	DrawGammaCanvasSettings( canvasMassPlusFWHM900GeV, 0.07, 0.02, 0.04, 0.10);	
	canvasMassPlusFWHM900GeV->cd();

	TH2D *histo2DAllPi0900GeVMassAndWidthCombined;
	histo2DAllPi0900GeVMassAndWidthCombined = new TH2D("histo2DAllPi0900GeVMassAndWidthCombined", "histo2DAllPi0900GeVMassAndWidthCombined", 20,0.,maxPtPi0900GeV ,1000.,-17,17);
	histo2DAllPi0900GeVMassAndWidthCombined->SetYTitle("(m- m_{PDG}- m_{o}, #sigma (MeV/#it{c}^{2})");
	histo2DAllPi0900GeVMassAndWidthCombined->SetXTitle("#it{p}_{T} (GeV/#it{c})");
	histo2DAllPi0900GeVMassAndWidthCombined->GetXaxis()->SetNdivisions(515,kTRUE);
	histo2DAllPi0900GeVMassAndWidthCombined->GetYaxis()->SetNdivisions(505,kTRUE);
	histo2DAllPi0900GeVMassAndWidthCombined->GetYaxis()->SetLabelSize(0.035);
	histo2DAllPi0900GeVMassAndWidthCombined->GetYaxis()->SetTitleSize(0.05);	
	histo2DAllPi0900GeVMassAndWidthCombined->GetYaxis()->SetDecimals();
	histo2DAllPi0900GeVMassAndWidthCombined->GetYaxis()->SetTitleOffset(0.63);
	histo2DAllPi0900GeVMassAndWidthCombined->GetXaxis()->SetTitleSize(0.05);	
	histo2DAllPi0900GeVMassAndWidthCombined->GetXaxis()->SetLabelSize(0.035);
	histo2DAllPi0900GeVMassAndWidthCombined->SetTitle("");
	histo2DAllPi0900GeVMassAndWidthCombined->GetXaxis()->SetTitleOffset(0.9);
	histo2DAllPi0900GeVMassAndWidthCombined->DrawCopy();

	/*	DrawGammaSetMarker(histoMassMesonPi0PHOSMinusExp, markerStylePHOS, markerSizeMass, colorPHOS,colorPHOS);
	histoMassMesonPi0PHOSMinusExp->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoTrueMassMesonPi0PHOSMinusExp, colorPHOS, markerSizeMass, colorPHOSMC, colorPHOSMC);
	histoTrueMassMesonPi0PHOSMinusExp->GetXaxis()->SetRangeUser(0.4,25.);
	histoTrueMassMesonPi0PHOSMinusExp->DrawCopy("same,e1,p"); */

	/*	DrawGammaSetMarker(histoFWHMMesonPi0PHOS, markerStylePHOS, markerSizeMass, colorPHOS, colorPHOS);
	histoFWHMMesonPi0PHOS->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoTrueFWHMMesonPi0PHOS, markerStylePHOSMC, markerSizeMass, colorPHOSMC , colorPHOSMC);
	histoTrueFWHMMesonPi0PHOS->GetXaxis()->SetRangeUser(0.6,25.);
	histoTrueFWHMMesonPi0PHOS->DrawCopy("same,e1,p"); */

	DrawGammaSetMarker(histoTrueMassMesonPi0MinusExp900GeV , markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);					 
	histoTrueMassMesonPi0MinusExp900GeV->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoMassMesonPi0MinusExp900GeV, markerStyleConv, markerSizeMass, colorConv, colorConv);					 
	histoMassMesonPi0MinusExp900GeV->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoFWHMMesonPi0MeV900GeV, markerStylePHOS , markerSizeMass, colorConv, colorConv);
	histoFWHMMesonPi0MeV900GeV->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoTrueFWHMMesonPi0MeV900GeV, markerStylePHOSMC, markerSizeMass, colorConvMC, colorConvMC);
	histoTrueFWHMMesonPi0MeV900GeV->DrawCopy("same,e1,p"); 

	DrawGammaLines(0., maxPtPi0900GeV ,-5., -5.,0.1,colorConv);
/*	DrawGammaLines(0., maxPtPi0900GeV ,-10., -10.,0.1,colorPHOS );
	DrawGammaLines(0., maxPtPi0900GeV,-15.,-15.,0.1, colorEMCAL);*/
	
	TPad* padMassLegendPi0900GeV = new TPad("padMassLegendPi0900GeV", "", 0.15, 0.12, 0.95, 0.27,-1, -1, -2);
	DrawGammaPadSettings( padMassLegendPi0900GeV, 0., 0., 0., 0.);
	padMassLegendPi0900GeV->Draw();
	padMassLegendPi0900GeV->cd();

	//****************** first Column **************************************************
	textMassCTS->Draw();
/*	textMassPHOS->Draw();
	textMassEMCAL->Draw();*/
	
	//****************** second Column *************************************************
	textMass->Draw();
	textMassData->Draw();
	textMassMC->Draw();
	textMassOffset->Draw();
	
	markerCTSPi07TeVMass->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
/*	markerPHOSPi07TeVMass->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	markerEMCALPi07TeVMass->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);*/
	
	markerCTSPi07TeVMassMC->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
/*	markerPHOSPi07TeVMassMC->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	markerEMCALPi07TeVMassMC->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);*/
	
	textOffsetConv->Draw();
/*	textOffsetPHOS->Draw();
	textOffsetEMCAL->Draw();*/
	
	//***************** third Column ***************************************************
	textWidth->Draw();
	textWidthData->Draw();
	textWidthMC->Draw();
	
	markerCTSPi07TeVWidth->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
/*	markerPHOSPi07TeVWidth->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	markerEMCALPi07TeVWidth->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);*/
	
	markerCTSPi07TeVWidthMC->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
/*	markerPHOSPi07TeVWidthMC->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	markerEMCALPi07TeVWidthMC->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);*/
	
	textWidthConv->Draw();
/*	textWidthPHOS->Draw();
	textWidthEMCAL->Draw();*/
	
	canvasMassPlusFWHM900GeV->Update();
	canvasMassPlusFWHM900GeV->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurements_900GeV_Paper.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));

	canvasMassPlusFWHM900GeV->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSPerformance(pictDrawingCoordinatesMassFWHM[0], pictDrawingCoordinatesMassFWHM[1], pictDrawingCoordinatesMassFWHM[2], pictDrawingCoordinatesMassFWHM[3], pictDrawingCoordinatesMassFWHM[4], pictDrawingCoordinatesMassFWHM[5], pictDrawingCoordinatesMassFWHM[6], pictDrawingCoordinatesMassFWHM[7], pictDrawingCoordinates[8], collisionSystem900GeV, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],2700,1800,date);
	
	canvasMassPlusFWHM900GeV->Update();
	canvasMassPlusFWHM900GeV->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurements_900GeV.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));




	//*******************************************************************************************************
	//************************************ Eta @ 7TeV *******************************************************
	//*******************************************************************************************************


	//************************************** Mass and FWHM combined *****************************************************************
	TCanvas* canvasMassPlusFWHMEta = new TCanvas("canvasMassPlusFWHMEta","",200,10,2700,1800);  // gives the page size
	DrawGammaCanvasSettings( canvasMassPlusFWHMEta, 0.07, 0.02, 0.04, 0.10);	
	canvasMassPlusFWHMEta->cd();
	canvasMassPlusFWHMEta->SetLogx(1);
	
	TH2D *histo2DAllEtaMassAndWidthCombined;
	histo2DAllEtaMassAndWidthCombined = new TH2D("histo2DAllEtaMassAndWidthCombined", "histo2DAllEtaMassAndWidthCombined", 20,0.27,15.,1000.,-60,60);
	histo2DAllEtaMassAndWidthCombined->SetYTitle("m- m_{PDG}- m_{o}, #sigma (MeV/#it{c}^{2})");
	histo2DAllEtaMassAndWidthCombined->SetXTitle("#it{p}_{T} (GeV/#it{c})");
	histo2DAllEtaMassAndWidthCombined->GetXaxis()->SetNdivisions(515,kTRUE);
	histo2DAllEtaMassAndWidthCombined->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DAllEtaMassAndWidthCombined->GetYaxis()->SetLabelSize(0.035);
	histo2DAllEtaMassAndWidthCombined->GetYaxis()->SetTitleSize(0.05);	
	histo2DAllEtaMassAndWidthCombined->GetYaxis()->SetDecimals();
	histo2DAllEtaMassAndWidthCombined->GetYaxis()->SetTitleOffset(0.63);
	histo2DAllEtaMassAndWidthCombined->GetXaxis()->SetTitleSize(0.05);	
	histo2DAllEtaMassAndWidthCombined->GetXaxis()->SetLabelSize(0.035);
	histo2DAllEtaMassAndWidthCombined->SetTitle("");
	histo2DAllEtaMassAndWidthCombined->GetXaxis()->SetTitleOffset(0.9);
	histo2DAllEtaMassAndWidthCombined->DrawCopy();

	DrawGammaSetMarker(histoMassMesonEtaEMCALMinusExp, markerStyleConv, markerSizeMass, colorEMCAL,colorEMCAL);
	histoMassMesonEtaEMCALMinusExp->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoTrueMassMesonEtaEMCALMinusExp, markerStyleConvMC, markerSizeMass, colorEMCALMC, colorEMCALMC);
	histoTrueMassMesonEtaEMCALMinusExp->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoFWHMMesonEtaEMCAL, markerStylePHOS, markerSizeMass, colorEMCAL, colorEMCAL);
	histoFWHMMesonEtaEMCAL->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoTrueFWHMMesonEtaEMCAL, markerStylePHOSMC, markerSizeMass, colorEMCALMC, colorEMCALMC);
	histoTrueFWHMMesonEtaEMCAL->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoTrueMassMesonEtaMinusExp, markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);					 
	histoTrueMassMesonEtaMinusExp->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoMassMesonEtaMinusExp, markerStyleConv, markerSizeMass, colorConv, colorConv);					 
	histoMassMesonEtaMinusExp->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoFWHMMesonEtaMeV, markerStylePHOS, markerSizeMass, colorConv, colorConv);
	histoFWHMMesonEtaMeV->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoTrueFWHMMesonEtaMeV, markerStylePHOSMC, markerSizeMass, colorConvMC, colorConvMC);
	histoTrueFWHMMesonEtaMeV->DrawCopy("same,e1,p"); 

	DrawGammaLines(0., maxPtEta,-5., -5.,0.1,colorConv);
	//DrawGammaLines(0., maxPtEta ,-10., -10.,0.1,colorPHOS );
	DrawGammaLines(0., maxPtEta,-15.,-15.,0.1, colorEMCAL);

	TPad* padMassLegendEta7TeV = new TPad("padMassLegendEta7TeV", "", 0.15, 0.12, 0.95, 0.27,-1, -1, -2);
	DrawGammaPadSettings( padMassLegendEta7TeV, 0., 0., 0., 0.);
	padMassLegendEta7TeV->Draw();
	padMassLegendEta7TeV->cd();
	
	//****************** first Column **************************************************
	textMassCTS->Draw();
	//textMassPHOS->Draw();
	textMassEMCAL->Draw();
	
	//****************** second Column *************************************************
	textMass->Draw();
	textMassData->Draw();
	textMassMC->Draw();
	textMassOffset->Draw();
	
	markerCTSPi07TeVMass->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
	//markerPHOSPi07TeVMass->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	markerEMCALPi07TeVMass->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);
	
	markerCTSPi07TeVMassMC->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
	//markerPHOSPi07TeVMassMC->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	markerEMCALPi07TeVMassMC->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);
	
	textOffsetConv->Draw();
	//textOffsetPHOS->Draw();
	textOffsetEMCAL->Draw();
	
	//***************** third Column ***************************************************
	textWidth->Draw();
	textWidthData->Draw();
	textWidthMC->Draw();
	
	markerCTSPi07TeVWidth->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
		//markerPHOSPi07TeVWidth->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	markerEMCALPi07TeVWidth->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);
	
	markerCTSPi07TeVWidthMC->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
		//markerPHOSPi07TeVWidthMC->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
	markerEMCALPi07TeVWidthMC->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass ,rowsLegendMass[4]+ offsetMarkerYMass);
	
	textWidthConv->Draw();
	//textWidthPHOS->Draw();
	textWidthEMCAL->Draw();

	canvasMassPlusFWHMEta->Update();
	canvasMassPlusFWHMEta->SaveAs(Form("%s/%s_Eta_CombinedMassAndWidthALLMeasurementsLogX_7TeV_Paper.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));
	
	canvasMassPlusFWHMEta->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSPerformance(pictDrawingCoordinatesMassFWHM[0], pictDrawingCoordinatesMassFWHM[1], pictDrawingCoordinatesMassFWHM[2], pictDrawingCoordinatesMassFWHM[3], pictDrawingCoordinatesMassFWHM[4], pictDrawingCoordinatesMassFWHM[5], pictDrawingCoordinatesMassFWHM[6], pictDrawingCoordinatesMassFWHM[7], pictDrawingCoordinates[8], collisionSystem7TeV, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,2700,1800,date);
	
	canvasMassPlusFWHMEta->Update();
	canvasMassPlusFWHMEta->SaveAs(Form("%s/%s_Eta_CombinedMassAndWidthALLMeasurementsLogX_7TeV.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));


	//**************************************************************************************************
	//*********************** Inv Cross sections PHOS, CTS *********************************************
	//**************************************************************************************************

	cout << "hey hier bin ich" << endl;
	histoInvCrossSectionPi0PHOS= (TH1D*)histoPi0Phos7TeV->Clone();
	histoInvCrossSectionPi0PHOS->Scale(xSection7TeV*recalcBarn);
	histoInvCrossSectionPi0SysPHOS = (TH1D*)histoPi0PhosSys7TeV->Clone();
	histoInvCrossSectionPi0SysPHOS->Scale(xSection7TeV*recalcBarn);
	graphSysErrInvCrossSectionPi0PHOS = new TGraphAsymmErrors(histoInvCrossSectionPi0SysPHOS);	
// 	histoInvCrossSectionPi0EMCAL= (TH1D*)histoPi0EMCAL7TeV->Clone();
// 	histoInvCrossSectionPi0EMCAL->Scale(xSection7TeV*recalcBarn);
	
	cout << "hey hier bin ich" << endl;
	
	histoInvCrossSectionEtaPHOS= (TH1D*)histoEtaPhos7TeV->Clone();
	histoInvCrossSectionEtaPHOS->Scale(xSection7TeV*recalcBarn);
	histoInvCrossSectionEtaSysPHOS = (TH1D*)histoEtaPhosSys7TeV->Clone();
	histoInvCrossSectionEtaSysPHOS->Scale(xSection7TeV*recalcBarn);
	graphSysErrInvCrossSectionEtaPHOS = new TGraphAsymmErrors(histoInvCrossSectionEtaSysPHOS);	
	//graphInvCrossSectionEtaEMCAL= (TGraphErrors*)graphEtaEMCAL7TeV->Clone();
	//graphInvCrossSectionEtaEMCAL= ScaleGraph(graphInvCrossSectionEtaEMCAL,xSection7TeV*recalcBarn);

	cout << "hey hier bin ich" << endl;
	
	histoInvCrossSectionPi0PHOS2760GeV= (TH1D*)histoPi0Phos2760GeV->Clone();
	histoInvCrossSectionPi0PHOS2760GeV->Scale(xSection2760GeV*recalcBarn);
	cout << "hier" << endl;
	histoInvCrossSectionPi0SysPHOS2760GeV = (TH1D*)histoPi0PhosSys2760GeV->Clone();
	histoInvCrossSectionPi0SysPHOS2760GeV->Scale(xSection2760GeV*recalcBarn);
	cout << "hier" << endl;
	if (!conference) {
		histoInvCrossSectionPi0SysRAAPHOS2760GeV = (TH1D*)histoPi0PhosSysRAA2760GeV->Clone();
		histoInvCrossSectionPi0SysRAAPHOS2760GeV->Scale(xSection2760GeV*recalcBarn);
		cout << "hier" << endl;
		graphSysErrRAAInvCrossSectionPi0PHOS2760GeV = new TGraphAsymmErrors(histoInvCrossSectionPi0SysRAAPHOS2760GeV);	
	}
	graphSysErrInvCrossSectionPi0PHOS2760GeV = new TGraphAsymmErrors(histoInvCrossSectionPi0SysPHOS2760GeV);	
	
	
	cout << "hey hier bin ich" << endl;
	
	histoInvCrossSectionPi0PHOS900GeV= (TH1D*)histoPi0Phos900GeV->Clone();
	histoInvCrossSectionPi0PHOS900GeV->Scale(xSection900GeV*recalcBarn);
	histoInvCrossSectionPi0SysPHOS900GeV = (TH1D*)histoPi0PhosSys900GeV->Clone();
	histoInvCrossSectionPi0SysPHOS900GeV->Scale(xSection900GeV*recalcBarn);
	graphSysErrInvCrossSectionPi0PHOS900GeV = new TGraphAsymmErrors(histoInvCrossSectionPi0SysPHOS900GeV);	

	cout << "hey hier bin ich" << endl;
	//******************* CombinePoints 7TeV Pi0 **********************
	//mode 1 with material systematic errors, otherwise material errors out
	Double_t paramGraph[3] = {1.0e12,7.,0.13};
	graphInvCrossSectionPi0Comb7TeV = CombinePtPointsSpectra(	histoInvCrossSectionPi0, 		graphInvCrossSectionSysPi0,
												histoInvCrossSectionPi0PHOS, 	graphSysErrInvCrossSectionPi0PHOS,
												graphInvCrossSectionPi0Comb7TeVStatErr, graphInvCrossSectionPi0Comb7TeVSysErr,
												xPtLimits7TeVNewPHOS, 34, 1, 0,4);
	graphInvCrossSectionPi0Comb7TeVStatErr->Print();
	graphInvCrossSectionPi0Comb7TeVSysErr->Print();
	fitInvCrossSectionPi0 = FitObject("l","fitInvCrossSectionPi0","Pi0",histoInvCrossSectionPi0,0.3,maxPtPi0 ,paramGraph,"QNRMEX0+");
	graphInvCrossSectionPi0Comb7TeV->Fit(fitInvCrossSectionPi0,"QNRMEX0+","",0.3,maxPtPi0);
	paramGraph[0] = fitInvCrossSectionPi0->GetParameter(0);
	paramGraph[1] = fitInvCrossSectionPi0->GetParameter(1);
	paramGraph[2] = fitInvCrossSectionPi0->GetParameter(2);
	
	if(bWCorrection.CompareTo("X")==0){
		fitTsallisPi07TeVPtMult = FitObject("tmpt","TsallisMultWithPtPi07TeV","Pi0");
		fitTsallisPi07TeVPtMult->SetParameters(paramGraph[0],paramGraph[1], paramGraph[2]) ; // standard parameter optimize if necessary
		graphInvCrossSectionPi0Comb7TeV   = ApplyXshift(graphInvCrossSectionPi0Comb7TeV ,fitTsallisPi07TeVPtMult);
	}
	graphInvCrossSectionPi0Comb7TeV->SetMarkerStyle(markerStyleConv);
	graphInvCrossSectionPi0Comb7TeV->SetMarkerColor(1);
	graphInvCrossSectionPi0Comb7TeV->SetLineColor(1);
	graphInvCrossSectionPi0Comb7TeV->SetLineWidth(5);
	graphInvCrossSectionPi0Comb7TeV->SetFillColor(0);
	graphInvCrossSectionPi0Comb7TeV->SetFillStyle(0);
	graphInvCrossSectionPi0Comb7TeVUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb7TeV->Clone("graphInvCrossSectionPi0Comb7TeVUnscaled");
	fitInvCrossSectionPi0 = FitObject("l","fitInvCrossSectionPi0","Pi0",graphInvCrossSectionPi0Comb7TeV,0.3,maxPtPi0 ,paramGraph,"QNRMEX0+");
	cout << WriteParameterToFile(fitInvCrossSectionPi0)<< endl;	
	cout << "Pi0 7TeV ____________________________________________" << endl;
	
	graphConversionXSectionPi07TeV = new TGraphAsymmErrors(histoInvCrossSectionPi0);
	graphConversionXSectionPi07TeV->RemovePoint(0);
	graphPHOSXSectionPi07TeV = new TGraphAsymmErrors(histoInvCrossSectionPi0PHOS);
	graphEMCALXSectionPi07TeV = new TGraphAsymmErrors(histoInvCrossSectionPi0EMCAL);
	
	if(bWCorrection.CompareTo("X")==0){
		graphConversionXSectionPi07TeV= ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb7TeV, graphConversionXSectionPi07TeV, fitTsallisPi07TeVPtMult, 0, 28);
		graphInvCrossSectionSysPi0= ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb7TeV,graphInvCrossSectionSysPi0, fitTsallisPi07TeVPtMult, 0, 28);
		graphPHOSXSectionPi07TeV = ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb7TeV, graphPHOSXSectionPi07TeV, fitTsallisPi07TeVPtMult,4,5+28);
		graphSysErrInvCrossSectionPi0PHOS = ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb7TeV, graphSysErrInvCrossSectionPi0PHOS, fitTsallisPi07TeVPtMult,4,5+28);
		graphInvCrossSectionPi0Comb7TeVStatErr = ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb7TeV, graphInvCrossSectionPi0Comb7TeVStatErr, fitTsallisPi07TeVPtMult, 0, graphInvCrossSectionPi0Comb7TeVStatErr->GetN());
		graphInvCrossSectionPi0Comb7TeVSysErr = ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb7TeV, graphInvCrossSectionPi0Comb7TeVSysErr, fitTsallisPi07TeVPtMult, 0, graphInvCrossSectionPi0Comb7TeVSysErr->GetN());
	}
	graphInvCrossSectionPi0Comb7TeVStatErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb7TeVStatErr->Clone("graphInvCrossSectionPi0Comb7TeVStatErrUnscaled");
	graphInvCrossSectionPi0Comb7TeVSysErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb7TeVSysErr->Clone("graphInvCrossSectionPi0Comb7TeVSysErrUnscaled");
	
	graphConversionXSectionPi07TeVUnscaled = (TGraphAsymmErrors*) graphConversionXSectionPi07TeV->Clone();	
	graphInvCrossSectionSysPi0Unscaled =  (TGraphAsymmErrors*) graphInvCrossSectionSysPi0->Clone();	
	graphPHOSXSectionPi07TeVUnscaled =  (TGraphAsymmErrors*) graphPHOSXSectionPi07TeV->Clone();	
	graphSysErrInvCrossSectionPi0PHOSUnscaled = (TGraphAsymmErrors*) graphSysErrInvCrossSectionPi0PHOS->Clone();	
	
	forOutput= WriteParameterToFile(fitInvCrossSectionPi0);
	fileFinalResults<< forOutput.Data()<< endl;
	
	graphRatioCombPHOSPi0 = (TGraphAsymmErrors*) graphPHOSXSectionPi07TeV->Clone();	
	graphRatioCombPHOSPi0Sys = (TGraphAsymmErrors*) graphSysErrInvCrossSectionPi0PHOS->Clone();	
	graphRatioCombConvPi0 = (TGraphAsymmErrors*) graphConversionXSectionPi07TeV->Clone();	
	graphRatioCombConvPi0Sys = (TGraphAsymmErrors*) graphInvCrossSectionSysPi0->Clone();	
// 	histoRatioCombEMCALPi0 = (TH1D*) histoInvCrossSectionPi0EMCAL->Clone();	
	
	graphRatioCombNLOPi07TeVMuHalf= (TGraph*)graphNLOMuHalfPi07TeV->Clone();
	graphRatioCombNLOPi07TeVMuOne= (TGraph*)graphNLOMuOnePi07TeV->Clone();
	graphRatioCombNLOPi07TeVMuTwo= (TGraph*)graphNLOMuTwoPi07TeV->Clone();
	graphRatioCombNLOBKKPi07TeVMuTwo = (TGraph*)graphNLOBKKCalcMuTwo7TeV->Clone();
	graphRatioCombNLODSSPi07TeVMuTwo = (TGraph*)graphNLODSSCalcMuTwo7TeV->Clone();
	
	graphRatioCombPHOSPi0 = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0, fitInvCrossSectionPi0); 
	graphRatioCombPHOSPi0Sys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0Sys, fitInvCrossSectionPi0); 
	
	graphRatioCombConvPi0 = CalculateGraphErrRatioToFit (graphRatioCombConvPi0, fitInvCrossSectionPi0); 
	graphRatioCombConvPi0Sys= CalculateGraphErrRatioToFit (graphRatioCombConvPi0Sys, fitInvCrossSectionPi0); 
	
// 	histoRatioCombEMCALPi0 = CalculateHistoRatioToFit (histoRatioCombEMCALPi0,fitInvCrossSectionPi0); 
	graphRatioCombNLOPi07TeVMuHalf = CalculateGraphRatioToFit (graphRatioCombNLOPi07TeVMuHalf, fitInvCrossSectionPi0); 
	graphRatioCombNLOPi07TeVMuOne = CalculateGraphRatioToFit (graphRatioCombNLOPi07TeVMuOne , fitInvCrossSectionPi0); 
	graphRatioCombNLOPi07TeVMuTwo = CalculateGraphRatioToFit (graphRatioCombNLOPi07TeVMuTwo, fitInvCrossSectionPi0); 
	graphRatioCombNLOBKKPi07TeVMuTwo = CalculateGraphRatioToFit (graphRatioCombNLOBKKPi07TeVMuTwo, fitInvCrossSectionPi0); 
	graphRatioCombNLODSSPi07TeVMuTwo = CalculateGraphRatioToFit (graphRatioCombNLODSSPi07TeVMuTwo, fitInvCrossSectionPi0); 
	cout << "\n\n\n\n\n ratios\n\nVogelsang\n" << endl;
	graphRatioCombNLOPi07TeVMuTwo->Print();
	cout << "\n\n BKK \n" << endl;
	graphRatioCombNLOBKKPi07TeVMuTwo->Print();
	cout << "\n\n DSS \n" << endl;
	graphRatioCombNLODSSPi07TeVMuTwo->Print();
	cout << "\n\n ende \n\n\n" << endl;
	
	graphRatioCombCombFit = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb7TeV->Clone();
	graphRatioCombCombFit = CalculateGraphErrRatioToFit(graphRatioCombCombFit, fitInvCrossSectionPi0); 
	graphRatioCombCombFitStat = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb7TeVStatErr->Clone();
	graphRatioCombCombFitStat = CalculateGraphErrRatioToFit(graphRatioCombCombFitStat, fitInvCrossSectionPi0); 
	graphRatioCombCombFitSys = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb7TeVSysErr->Clone();
	graphRatioCombCombFitSys = CalculateGraphErrRatioToFit(graphRatioCombCombFitSys, fitInvCrossSectionPi0); 
	
	Double_t paramGraphEta[3] = {1.0e10,7.,0.13};
	graphInvCrossSectionEtaComb7TeV =  CombinePtPointsSpectra(	histoInvCrossSectionEta, 		graphInvCrossSectionSysEta,
												histoInvCrossSectionEtaPHOS, 	graphSysErrInvCrossSectionEtaPHOS,
												graphInvCrossSectionEtaComb7TeVStatErr, graphInvCrossSectionEtaComb7TeVSysErr,
												xPtLimitsEta7TeV, 14, 1, 0,2);
	
	if(bWCorrection.CompareTo("X")==0){
		fitTsallisEta7TeVPtMult = FitObject("tmpt","TsallisMultWithPtEta7TeV","Eta");
		fitTsallisEta7TeVPtMult->SetParameters(paramGraphEta[0],paramGraphEta[1], paramGraphEta[2]) ; // standard parameter optimize if necessary
		graphInvCrossSectionEtaComb7TeV   = ApplyXshift(graphInvCrossSectionEtaComb7TeV ,fitTsallisEta7TeVPtMult);
	}
	graphInvCrossSectionEtaComb7TeV->SetMarkerStyle(markerStyleConv);
	graphInvCrossSectionEtaComb7TeV->SetMarkerColor(1);
	graphInvCrossSectionEtaComb7TeV->SetLineColor(1);
	graphInvCrossSectionEtaComb7TeV->SetLineWidth(5);
	graphInvCrossSectionEtaComb7TeV->SetFillColor(0);
	graphInvCrossSectionEtaComb7TeV->SetFillStyle(0);
	graphInvCrossSectionEtaComb7TeVUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb7TeV->Clone("graphInvCrossSectionEtaComb7TeVUnscaled");
	
	fitInvCrossSectionEta7TeV = FitObject("l","fitInvCrossSectionEta7TeV","Eta",histoInvCrossSectionEta,0.4,maxPtEta,paramGraphEta,"QNRMEX0+");
	graphInvCrossSectionEtaComb7TeV->Fit(fitInvCrossSectionEta7TeV,"QNRMEX0+","",0.4,maxPtEta);
	cout << WriteParameterToFile(fitInvCrossSectionEta7TeV)<< endl;
	cout << "Eta 7TeV ____________________________________________" << endl;	
	graphConversionXSectionEta7TeV = new TGraphAsymmErrors(histoInvCrossSectionEta);
	graphConversionXSectionEta7TeV->RemovePoint(0);
	graphPHOSXSectionEta7TeV = new TGraphAsymmErrors(histoInvCrossSectionEtaPHOS);
	
	if(bWCorrection.CompareTo("X")==0){
		graphConversionXSectionEta7TeV= ApplyXshiftIndividualSpectra(graphInvCrossSectionEtaComb7TeV, graphConversionXSectionEta7TeV, fitTsallisEta7TeVPtMult, 0,11);
		graphInvCrossSectionSysEta= ApplyXshiftIndividualSpectra(graphInvCrossSectionEtaComb7TeV, graphInvCrossSectionSysEta, fitTsallisEta7TeVPtMult, 0,11);
		graphPHOSXSectionEta7TeV= ApplyXshiftIndividualSpectra(graphInvCrossSectionEtaComb7TeV, graphPHOSXSectionEta7TeV, fitTsallisEta7TeVPtMult,2,2+11);
		graphSysErrInvCrossSectionEtaPHOS= ApplyXshiftIndividualSpectra(graphInvCrossSectionEtaComb7TeV, graphSysErrInvCrossSectionEtaPHOS, fitTsallisEta7TeVPtMult,2,2+11);
		graphInvCrossSectionEtaComb7TeVStatErr = ApplyXshiftIndividualSpectra(graphInvCrossSectionEtaComb7TeV, graphInvCrossSectionEtaComb7TeVStatErr, fitTsallisEta7TeVPtMult,0,graphInvCrossSectionEtaComb7TeVStatErr->GetN());
		graphInvCrossSectionEtaComb7TeVSysErr = ApplyXshiftIndividualSpectra(graphInvCrossSectionEtaComb7TeV, graphInvCrossSectionEtaComb7TeVSysErr, fitTsallisEta7TeVPtMult,0,graphInvCrossSectionEtaComb7TeVSysErr->GetN());
	}
	graphInvCrossSectionEtaComb7TeVStatErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb7TeVStatErr->Clone("graphInvCrossSectionEtaComb7TeVStatErrUnscaled");
	graphInvCrossSectionEtaComb7TeVSysErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb7TeVSysErr->Clone("graphInvCrossSectionEtaComb7TeVSysErrUnscaled");
	
	graphConversionXSectionEta7TeVUnscaled = (TGraphAsymmErrors*) graphConversionXSectionEta7TeV->Clone();	
	graphInvCrossSectionSysEtaUnscaled =  (TGraphAsymmErrors*) graphInvCrossSectionSysEta->Clone();	
	graphPHOSXSectionEta7TeVUnscaled =  (TGraphAsymmErrors*) graphPHOSXSectionEta7TeV->Clone();	
	graphSysErrInvCrossSectionEtaPHOSUnscaled = (TGraphAsymmErrors*) graphSysErrInvCrossSectionEtaPHOS->Clone();	
	
	forOutput= WriteParameterToFile(fitInvCrossSectionEta7TeV);
	fileFinalResults<< forOutput.Data()<< endl;
	
	graphRatioCombNLOEta7TeVMuHalf= (TGraph*)graphNLOMuHalfEta7TeV->Clone();
	graphRatioCombNLOEta7TeVMuOne= (TGraph*)graphNLOMuOneEta7TeV->Clone();
	graphRatioCombNLOEta7TeVMuTwo= (TGraph*)graphNLOMuTwoEta7TeV->Clone();	
	graphRatioCombPHOSEta7TeV = (TGraphAsymmErrors*) graphPHOSXSectionEta7TeV->Clone();	
	graphRatioCombPHOSEta7TeVSys = (TGraphAsymmErrors*) graphSysErrInvCrossSectionEtaPHOS->Clone();	
	graphRatioCombConvEta7TeV = (TGraphAsymmErrors*) graphConversionXSectionEta7TeV->Clone();	
	graphRatioCombConvEta7TeVSys = (TGraphAsymmErrors*) graphInvCrossSectionSysEta->Clone();	
//	graphRatioCombEMCALEta7TeV = (TGraphErrors*) graphInvCrossSectionEtaEMCAL->Clone();
	
	graphRatioCombPHOSEta7TeV = CalculateGraphErrRatioToFit (graphRatioCombPHOSEta7TeV, fitInvCrossSectionEta7TeV); 
	graphRatioCombPHOSEta7TeVSys = CalculateGraphErrRatioToFit (graphRatioCombPHOSEta7TeVSys, fitInvCrossSectionEta7TeV); 
	graphRatioCombConvEta7TeV = CalculateGraphErrRatioToFit (graphRatioCombConvEta7TeV, fitInvCrossSectionEta7TeV); 
	graphRatioCombConvEta7TeVSys = CalculateGraphErrRatioToFit (graphRatioCombConvEta7TeVSys, fitInvCrossSectionEta7TeV); 
//	graphRatioCombEMCALEta7TeV = CalculateGraphErrRatioToFit (graphRatioCombEMCALEta7TeV, fitInvCrossSectionEta7TeV); 
	
	graphRatioCombNLOEta7TeVMuHalf = CalculateGraphRatioToFit (graphRatioCombNLOEta7TeVMuHalf, fitInvCrossSectionEta7TeV); 
	graphRatioCombNLOEta7TeVMuOne = CalculateGraphRatioToFit (graphRatioCombNLOEta7TeVMuOne, fitInvCrossSectionEta7TeV); 
	graphRatioCombNLOEta7TeVMuTwo = CalculateGraphRatioToFit (graphRatioCombNLOEta7TeVMuTwo, fitInvCrossSectionEta7TeV); 
	
	graphRatioCombCombFitEta7TeV = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb7TeV->Clone();
	graphRatioCombCombFitEta7TeV = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta7TeV, fitInvCrossSectionEta7TeV); 
	graphRatioCombCombFitEta7TeVStat = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb7TeVStatErr->Clone();
	graphRatioCombCombFitEta7TeVStat = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta7TeVStat, fitInvCrossSectionEta7TeV); 
	graphRatioCombCombFitEta7TeVSys = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb7TeVSysErr->Clone();
	graphRatioCombCombFitEta7TeVSys = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta7TeVSys, fitInvCrossSectionEta7TeV); 
	
	histoFitInvCrossSectionEta7TeV = (TH1D*)fitInvCrossSectionEta7TeV->GetHistogram();
	histoFitInvCrossSectionEta7TeV->Scale(1e-3);
	
	//**************** CombinePoints 2760 GeV**************************
	cout << "Input Graph Conv" << endl;
	if (conference){
		graphInvCrossSectionPi0Comb2760GeV = CombinePtPointsSpectra(	histoInvCrossSectionPi02760GeV, 		graphInvCrossSectionSysPi02760GeV,
												histoInvCrossSectionPi0PHOS2760GeV, 	graphSysErrInvCrossSectionPi0PHOS2760GeV,
												graphInvCrossSectionPi0Comb2760GeVStatErr, graphInvCrossSectionPi0Comb2760GeVSysErr,
												xPtLimits2760GeV, 20, 1, 0,1);
	} else {
		histoInvCrossSectionPi0PHOS2760GeV->SetBinContent(1,0.);
		histoInvCrossSectionPi0PHOS2760GeV->SetBinError(1,0.);
		graphInvCrossSectionPi0Comb2760GeV = CombinePtPointsSpectra(	histoInvCrossSectionPi02760GeV, 		graphInvCrossSectionSysPi02760GeV,
												histoInvCrossSectionPi0PHOS2760GeV, 	graphSysErrInvCrossSectionPi0PHOS2760GeV,
												graphInvCrossSectionPi0Comb2760GeVStatErr, graphInvCrossSectionPi0Comb2760GeVSysErr,
												xPtLimits2760GeVNew, 20, 1, 0,1,kTRUE);
	}
 	TGraphAsymmErrors* graphInvCrossSectionPi0Comb2760GeVUnShifted = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeV->Clone("Unshifted"); 
	TGraphAsymmErrors* graphInvCrossSectionPi0Comb2760GeVStatErrUnShifted = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVStatErr->Clone("UnshiftedStat"); 
	TGraphAsymmErrors* graphInvCrossSectionPi0Comb2760GeVSysErrUnShifted = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVSysErr->Clone("Unshifted"); 
	graphConversionXSectionPi02760GeV = new TGraphAsymmErrors(histoInvCrossSectionPi02760GeV);
	graphConversionXSectionPi02760GeV->RemovePoint(0);
	graphConversionXSectionPi02760GeV->RemovePoint(graphConversionXSectionPi02760GeV->GetN()-1);
	TGraphAsymmErrors* graphConversionXSectionPi02760GeVUnShifted = (TGraphAsymmErrors*)graphConversionXSectionPi02760GeV->Clone("UnshiftedConv"); 
	graphInvCrossSectionSysPi02760GeV->RemovePoint(graphInvCrossSectionSysPi02760GeV->GetN()-1);
	TGraphAsymmErrors* graphConversionXSectionSysPi02760GeVUnShifted = (TGraphAsymmErrors*)graphInvCrossSectionSysPi02760GeV->Clone("UnshiftedSysConv"); 
	graphPHOSXSectionPi02760GeV = new TGraphAsymmErrors(histoInvCrossSectionPi0PHOS2760GeV);
	graphPHOSXSectionPi02760GeV->RemovePoint(0);
	TGraphAsymmErrors* graphPHOSXSectionPi02760GeVUnshifted = (TGraphAsymmErrors*)graphPHOSXSectionPi02760GeV->Clone("UnshiftedPhos"); 
	graphSysErrInvCrossSectionPi0PHOS2760GeV->RemovePoint(0);
	TGraphAsymmErrors* graphPHOSXSectionSysPi02760GeVUnshifted = (TGraphAsymmErrors*)graphSysErrInvCrossSectionPi0PHOS2760GeV->Clone("UnshiftedSysPhos"); 
	
	fitInvCrossSectionPi02760GeV = FitObject("l","fitInvCrossSectionPi02760GeV","Pi0",histoInvCrossSectionPi02760GeV,0.4,12.,paramGraph,"QNRMEX0+");
	TF1* fitInvCrossSectionPi02760GeVGraph = (TF1*)fitInvCrossSectionPi02760GeV->Clone("fitInvCrossSectionPi02760GeVGraph"); 
	
	if(bWCorrection.CompareTo("X")==0 && !conference){
		TF1* fitTsallisPi02760GeVPtMult = FitObject("tmpt","TsallisMultWithPtPi02760GeV","Pi0");
		fitTsallisPi02760GeVPtMult->SetParameters(paramGraph[0],paramGraph[1], paramGraph[2]) ; // standard parameter optimize if necessary
// 		fitTsallisPi02760GeVPtMult->SetRange();
		graphInvCrossSectionPi0Comb2760GeV   = ApplyXshift(graphInvCrossSectionPi0Comb2760GeV ,fitTsallisPi02760GeVPtMult);
		
		graphInvCrossSectionPi0Comb2760GeVStatErr = ApplyXshiftIndividualSpectra (graphInvCrossSectionPi0Comb2760GeV, graphInvCrossSectionPi0Comb2760GeVStatErr, fitTsallisPi02760GeVPtMult, 0, graphInvCrossSectionPi0Comb2760GeVStatErr->GetN());
		graphInvCrossSectionPi0Comb2760GeVSysErr = ApplyXshiftIndividualSpectra (graphInvCrossSectionPi0Comb2760GeV, graphInvCrossSectionPi0Comb2760GeVSysErr, fitTsallisPi02760GeVPtMult, 0, graphInvCrossSectionPi0Comb2760GeVSysErr->GetN());
		graphConversionXSectionPi02760GeV= ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb2760GeV, graphConversionXSectionPi02760GeV, fitTsallisPi07TeVPtMult, 0, 17);
		graphInvCrossSectionSysPi02760GeV= ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb2760GeV,graphInvCrossSectionSysPi02760GeV, fitTsallisPi07TeVPtMult, 0, 17);
		graphPHOSXSectionPi02760GeV = ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb2760GeV, graphPHOSXSectionPi02760GeV, fitTsallisPi07TeVPtMult,2,graphInvCrossSectionPi0Comb2760GeV->GetN());
		graphSysErrInvCrossSectionPi0PHOS2760GeV = ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb2760GeV, graphSysErrInvCrossSectionPi0PHOS2760GeV, fitTsallisPi07TeVPtMult,2,graphInvCrossSectionPi0Comb2760GeV->GetN());
		
		TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
		DrawGammaCanvasSettings( canvasDummy2,  0.1, 0.01, 0.015, 0.08);
		canvasDummy2->SetLogy();
		canvasDummy2->SetLogx();
		TH2F * histo2DDummy2;
		histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,30.,1000,10e1,10e11);
		SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
		histo2DDummy2->DrawCopy(); 
	
		DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb2760GeVStatErrUnShifted, 20,markerSizeCommonSpectrumPi0900GeV, kRed, kRed, widthLinesBoxes, kTRUE);
		graphInvCrossSectionPi0Comb2760GeVStatErrUnShifted->Draw("pEsame");
		DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb2760GeV, 24,markerSizeCommonSpectrumPi07TeV, kBlack, kBlack, widthLinesBoxes, kTRUE);
		graphInvCrossSectionPi0Comb2760GeV->Draw("pEsame");
		DrawGammaSetMarkerTGraphAsym(graphConversionXSectionPi02760GeV, 25,markerSizeCommonSpectrumPi07TeV, kBlue, kBlue, widthLinesBoxes, kTRUE);
		graphConversionXSectionPi02760GeV->Draw("pEsame");
		DrawGammaSetMarkerTGraphAsym(graphSysErrInvCrossSectionPi0PHOS2760GeV, 26,markerSizeCommonSpectrumPi07TeV, kGreen+2, kGreen+2, widthLinesBoxes, kTRUE);
		graphSysErrInvCrossSectionPi0PHOS2760GeV->Draw("pEsame");

		graphInvCrossSectionPi0Comb2760GeV->Fit(fitInvCrossSectionPi02760GeVGraph,"QNRMVE+","",0.4,10.);
		
		fitInvCrossSectionPi02760GeVGraph->SetLineColor(kBlue+2);
		fitInvCrossSectionPi02760GeVGraph->Draw("same");
		
// 		fitInvCrossSectionPi02760GeV->SetLineColor(kRed+2);
// 		fitInvCrossSectionPi02760GeV->Draw("same");
		
		canvasDummy2->Update();
		canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_2760GeV.%s",outputDir.Data(),suffix.Data()));

// 		return;
	}
	
	graphInvCrossSectionPi0Comb2760GeV->SetMarkerStyle(markerStyleConv);
	graphInvCrossSectionPi0Comb2760GeV->SetMarkerColor(1);
	graphInvCrossSectionPi0Comb2760GeV->SetLineColor(1);
	graphInvCrossSectionPi0Comb2760GeV->SetLineWidth(5);
	graphInvCrossSectionPi0Comb2760GeV->SetFillColor(0);
	graphInvCrossSectionPi0Comb2760GeV->SetFillStyle(0);
	graphInvCrossSectionPi0Comb2760GeVUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeV->Clone("graphInvCrossSectionPi0Comb2760GeVUnscaled");
	graphInvCrossSectionPi0Comb2760GeVStatErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVStatErr->Clone("graphInvCrossSectionPi0Comb2760GeVStatErrUnscaled");
	graphInvCrossSectionPi0Comb2760GeVSysErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVSysErr->Clone("graphInvCrossSectionPi0Comb2760GeVSysErrUnscaled");
	
	graphInvCrossSectionPi0Comb2760GeV->Fit(fitInvCrossSectionPi02760GeV,"QNRMEX0+","",0.4,10.);
	
// 	fitInvCrossSectionPi02760GeV = FitObject("l","fitInvCrossSectionPi02760GeV","Pi0",graphInvCrossSectionPi0Comb2760GeV,0.4,12.,paramGraph,"QNRMEX0+");
// 	fitInvCrossSectionPi02760GeV = FitObject("l","fitInvCrossSectionPi02760GeV","Pi0",graphInvCrossSectionPi0Comb2760GeV,0.4,12. ,paramGraph,"QNRMEX0+");

	cout << WriteParameterToFile(fitInvCrossSectionPi02760GeV)<< endl;
	forOutput= WriteParameterToFile(fitInvCrossSectionPi02760GeV);
	fileFinalResults<< forOutput.Data()<< endl;
	
	graphRatioCombNLOPi02760GeVMuHalf= (TGraph*)graphNLOMuHalfPi02760GeV->Clone();
	graphRatioCombNLOPi02760GeVMuOne= (TGraph*)graphNLOMuOnePi02760GeV->Clone();
	graphRatioCombNLOPi02760GeVMuTwo= (TGraph*)graphNLOMuTwoPi02760GeV->Clone();
	
	histoInvCrossSectionPi0PHOS2760GeV->SetBinContent(histoInvCrossSectionPi0PHOS2760GeV->GetNbinsX(),0);
	histoRatioCombPHOSPi02760GeV = (TH1D*) histoInvCrossSectionPi0PHOS2760GeV->Clone();	
	histoInvCrossSectionPi0SysPHOS2760GeV->SetBinContent(1,0);
	histoInvCrossSectionPi0SysPHOS2760GeV->SetBinContent(histoInvCrossSectionPi0SysPHOS2760GeV->GetNbinsX(),0);
	histoRatioCombPHOSPi0Sys2760GeV = (TH1D*) histoInvCrossSectionPi0SysPHOS2760GeV->Clone();	
	
	histoRatioCombPHOSPi02760GeV = CalculateHistoRatioToFit (histoRatioCombPHOSPi02760GeV, fitInvCrossSectionPi02760GeV); 
	histoRatioCombPHOSPi0Sys2760GeV = CalculateHistoRatioToFit (histoRatioCombPHOSPi0Sys2760GeV, fitInvCrossSectionPi02760GeV); 
	histoInvCrossSectionPi02760GeV->SetBinContent(histoInvCrossSectionPi02760GeV->GetNbinsX(),0);
	histoInvCrossSectionPi02760GeV->SetBinError(histoInvCrossSectionPi02760GeV->GetNbinsX(),0);
	histoRatioCombConvPi02760GeV = (TH1D*) histoInvCrossSectionPi02760GeV->Clone();		
	histoRatioCombConvPi02760GeV = CalculateHistoRatioToFit (histoRatioCombConvPi02760GeV, fitInvCrossSectionPi02760GeV); 
   TH1D* histoRatioPythia8ToFit2760GeV = (TH1D*) histoPythia8InvXSection->Clone();     
   histoRatioPythia8ToFit2760GeV = CalculateHistoRatioToFit (histoRatioPythia8ToFit2760GeV, fitInvCrossSectionPi02760GeV); 
   TH1D* histoRatioPythia8VarBinningToFit2760GeV = (TH1D*) histoPythia8InvXSection_VarBinning->Clone();     
   histoRatioPythia8VarBinningToFit2760GeV = CalculateHistoRatioToFit (histoRatioPythia8VarBinningToFit2760GeV, fitInvCrossSectionPi02760GeV); 
   
	graphRatioCombNLOPi02760GeVMuHalf = CalculateGraphRatioToFit (graphRatioCombNLOPi02760GeVMuHalf, fitInvCrossSectionPi02760GeV); 
	graphRatioCombNLOPi02760GeVMuOne = CalculateGraphRatioToFit (graphRatioCombNLOPi02760GeVMuOne, fitInvCrossSectionPi02760GeV); 
	graphRatioCombNLOPi02760GeVMuTwo = CalculateGraphRatioToFit (graphRatioCombNLOPi02760GeVMuTwo, fitInvCrossSectionPi02760GeV); 
	
	graphRatioCombCombFit2760GeV = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeV->Clone();
	graphRatioCombCombFit2760GeV = CalculateGraphErrRatioToFit(graphRatioCombCombFit2760GeV, fitInvCrossSectionPi02760GeV); 
	graphRatioCombCombFit2760GeVStat = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVStatErr->Clone();
	graphRatioCombCombFit2760GeVStat = CalculateGraphErrRatioToFit(graphRatioCombCombFit2760GeVStat, fitInvCrossSectionPi02760GeV); 
	graphRatioCombCombFit2760GeVSys = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVSysErr->Clone();
	graphRatioCombCombFit2760GeVSys = CalculateGraphErrRatioToFit(graphRatioCombCombFit2760GeVSys, fitInvCrossSectionPi02760GeV); 
	
	graphSysErrRatioCombPHOSPi02760GeV = new TGraphErrors(histoRatioCombPHOSPi0Sys2760GeV);	
	graphSysErrRatioCombConvPi02760GeV = CalculateSysErrFromRelSysHisto( histoRatioCombConvPi02760GeV, "Pi0SystError_RatioFitComb2760GeV",relSystErrorPi02760GeVDown, relSystErrorPi02760GeVUp, 2, nPointsPi02760GeV);
	histoFitInvCrossSectionPi02760GeV = (TH1D*)fitInvCrossSectionPi02760GeV->GetHistogram();

	paramGraphEta[1] =fitInvCrossSectionPi02760GeV->GetParameter(1);
	graphInvCrossSectionEtaComb2760GeVStatErr = new TGraphAsymmErrors(histoInvCrossSectionEta2760GeV);	
	graphInvCrossSectionEtaComb2760GeVStatErr->RemovePoint(0);
	graphInvCrossSectionEtaComb2760GeV = CalculateSysErrFromRelSysHistoComplete( histoInvCrossSectionEta2760GeV, "graphInvCrossSectionEtaComb2760GeV",relSystErrorEta2760GeVDown , relSystErrorEta2760GeVUp, 2, nPointsEta2760GeV);
	graphInvCrossSectionEtaComb2760GeVSysErr =  (TGraphAsymmErrors*)graphInvCrossSectionSysEta2760GeV->Clone("graphInvCrossSectionEtaComb2760GeVSysErr");
	if(bWCorrection.CompareTo("X")==0 && !conference){
		TF1* fitTsallisEta2760GeVPtMult = FitObject("tmpt","TsallisMultWithPtEta2760GeV","Eta");
		fitTsallisEta2760GeVPtMult->SetParameters(paramGraphEta[0],paramGraphEta[1], paramGraphEta[2]) ; // standard parameter optimize if necessary
		graphInvCrossSectionEtaComb2760GeV   = ApplyXshift(graphInvCrossSectionEtaComb2760GeV ,fitTsallisEta2760GeVPtMult);
		graphInvCrossSectionEtaComb2760GeVStatErr = ApplyXshiftIndividualSpectra (graphInvCrossSectionEtaComb2760GeV, graphInvCrossSectionEtaComb2760GeVStatErr, fitTsallisEta2760GeVPtMult, 0, graphInvCrossSectionEtaComb2760GeVStatErr->GetN());
		graphInvCrossSectionEtaComb2760GeVSysErr = ApplyXshiftIndividualSpectra (graphInvCrossSectionEtaComb2760GeV, graphInvCrossSectionEtaComb2760GeVSysErr, fitTsallisEta2760GeVPtMult, 0, graphInvCrossSectionEtaComb2760GeVSysErr->GetN());
	}
	graphInvCrossSectionEtaComb2760GeV->SetMarkerStyle(markerStyleConv);
	graphInvCrossSectionEtaComb2760GeV->SetMarkerColor(1);
	graphInvCrossSectionEtaComb2760GeV->SetLineColor(1);
	graphInvCrossSectionEtaComb2760GeV->SetLineWidth(5);
	graphInvCrossSectionEtaComb2760GeV->SetFillColor(0);
	graphInvCrossSectionEtaComb2760GeV->SetFillStyle(0);
	graphInvCrossSectionEtaComb2760GeVUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeV->Clone("graphInvCrossSectionEtaComb2760GeVUnscaled");
	graphInvCrossSectionEtaComb2760GeVStatErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVStatErr->Clone("graphInvCrossSectionEtaComb2760GeVStatErrUnscaled");
	graphInvCrossSectionEtaComb2760GeVSysErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVSysErr->Clone("graphInvCrossSectionEtaComb2760GeVSysErrUnscaled");
	
	fitInvCrossSectionEta2760GeV = FitObject("l","fitInvCrossSectionEta2760GeV","Eta",histoInvCrossSectionEta2760GeV,0.6,6.,paramGraphEta,"NQRME+","fixn");
	graphInvCrossSectionEtaComb2760GeV->Fit(fitInvCrossSectionEta2760GeV,"QNRMEX0+","",0.6,6.);
	
	cout << WriteParameterToFile(fitInvCrossSectionEta2760GeV)<< endl;
	forOutput= WriteParameterToFile(fitInvCrossSectionEta2760GeV);
	fileFinalResults<< forOutput.Data()<< endl;
	
	
	graphRatioCombNLOEta2760GeVMuHalf= (TGraph*)graphNLOMuHalfEta2760GeV->Clone();
	graphRatioCombNLOEta2760GeVMuOne= (TGraph*)graphNLOMuOneEta2760GeV->Clone();
	graphRatioCombNLOEta2760GeVMuTwo= (TGraph*)graphNLOMuTwoEta2760GeV->Clone();
	
	histoRatioCombConvEta2760GeV = (TH1D*) histoInvCrossSectionEta2760GeV->Clone();		
	histoRatioCombConvEta2760GeV = CalculateHistoRatioToFit (histoRatioCombConvEta2760GeV, fitInvCrossSectionEta2760GeV); 
	graphRatioCombNLOEta2760GeVMuHalf = CalculateGraphRatioToFit (graphRatioCombNLOEta2760GeVMuHalf, fitInvCrossSectionEta2760GeV); 
	graphRatioCombNLOEta2760GeVMuOne = CalculateGraphRatioToFit (graphRatioCombNLOEta2760GeVMuOne, fitInvCrossSectionEta2760GeV); 
	graphRatioCombNLOEta2760GeVMuTwo = CalculateGraphRatioToFit (graphRatioCombNLOEta2760GeVMuTwo, fitInvCrossSectionEta2760GeV); 
	
	graphRatioCombCombFitEta2760GeV = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeV->Clone();
	graphRatioCombCombFitEta2760GeV = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta2760GeV, fitInvCrossSectionEta2760GeV); 
	graphRatioCombCombFitEta2760GeVStat = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVStatErr->Clone();
	graphRatioCombCombFitEta2760GeVStat = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta2760GeVStat, fitInvCrossSectionEta2760GeV); 
	graphRatioCombCombFitEta2760GeVSys = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb2760GeVSysErr->Clone();
	graphRatioCombCombFitEta2760GeVSys = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta2760GeVSys, fitInvCrossSectionEta2760GeV); 
	
	histoFitInvCrossSectionEta2760GeV = (TH1D*)fitInvCrossSectionEta2760GeV->GetHistogram();
	
	//**************** CombinePoints 900 GeV**************************

	
	graphInvCrossSectionPi0Comb900GeV = CombinePtPointsSpectra(histoInvCrossSectionPi0900GeV, 		graphInvCrossSectionSysPi0900GeV,
												  histoInvCrossSectionPi0PHOS900GeV, 	graphSysErrInvCrossSectionPi0PHOS900GeV,
												  graphInvCrossSectionPi0Comb900GeVStatErr, graphInvCrossSectionPi0Comb900GeVSysErr,
												  xPtLimits900GeV, 14, 1, 0,1);
	if(bWCorrection.CompareTo("X")==0){
		fitTsallisPi0900GeVPtMult = FitObject("tmpt","TsallisMultWithPtPi0900GeV","Pi0");
		fitTsallisPi0900GeVPtMult->SetParameters(paramGraph[0],paramGraph[1], paramGraph[2]) ; // standard parameter optimize if necessary
		graphInvCrossSectionPi0Comb900GeV   = ApplyXshift(graphInvCrossSectionPi0Comb900GeV ,fitTsallisPi0900GeVPtMult);
	}
	
	graphInvCrossSectionPi0Comb900GeV->SetMarkerStyle(markerStyleConv);
	graphInvCrossSectionPi0Comb900GeV->SetMarkerColor(1);
	graphInvCrossSectionPi0Comb900GeV->SetLineColor(1);
	graphInvCrossSectionPi0Comb900GeV->SetLineWidth(5);
	graphInvCrossSectionPi0Comb900GeV->SetFillColor(0);
	graphInvCrossSectionPi0Comb900GeV->SetFillStyle(0);
	graphInvCrossSectionPi0Comb900GeVUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb900GeV->Clone("graphInvCrossSectionPi0Comb900GeVUnscaled");
	fitInvCrossSectionPi0900GeV = FitObject("l","fitInvCrossSectionPi0900GeV","Pi0",graphInvCrossSectionPi0Comb900GeV,0.4,7.,paramGraph,"QNRMEX0+");
	paramGraphEta[1] =fitInvCrossSectionPi0900GeV->GetParameter(1);
	
	cout << WriteParameterToFile(fitInvCrossSectionPi0900GeV)<< endl;
	cout << "Pi0 900GeV ____________________________________________" << endl;	
	graphConversionXSectionPi0900GeV = new TGraphAsymmErrors(histoInvCrossSectionPi0900GeV);
	graphConversionXSectionPi0900GeV->RemovePoint(0);
	graphPHOSXSectionPi0900GeV = new TGraphAsymmErrors(histoInvCrossSectionPi0PHOS900GeV);
	
	if(bWCorrection.CompareTo("X")==0){
		graphConversionXSectionPi0900GeV = ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb900GeV, graphConversionXSectionPi0900GeV, fitTsallisPi0900GeVPtMult,0,8);
		graphInvCrossSectionSysPi0900GeV = ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb900GeV, graphInvCrossSectionSysPi0900GeV, fitTsallisPi0900GeVPtMult,0,8);
		graphPHOSXSectionPi0900GeV= ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb900GeV, graphPHOSXSectionPi0900GeV, fitTsallisPi0900GeVPtMult,1,1+12);
		graphSysErrInvCrossSectionPi0PHOS900GeV= ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb900GeV, graphSysErrInvCrossSectionPi0PHOS900GeV, fitTsallisPi0900GeVPtMult,1,1+12);
		graphInvCrossSectionPi0Comb900GeVStatErr = ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb900GeV, graphInvCrossSectionPi0Comb900GeVStatErr, fitTsallisPi0900GeVPtMult,0,graphInvCrossSectionPi0Comb900GeV->GetN());
		graphInvCrossSectionPi0Comb900GeVSysErr = ApplyXshiftIndividualSpectra(graphInvCrossSectionPi0Comb900GeV, graphInvCrossSectionPi0Comb900GeVSysErr, fitTsallisPi0900GeVPtMult,0,graphInvCrossSectionPi0Comb900GeVSysErr->GetN());
	}
	graphInvCrossSectionPi0Comb900GeVStatErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb900GeVStatErr->Clone("graphInvCrossSectionPi0Comb900GeVStatErrUnscaled");
	graphInvCrossSectionPi0Comb900GeVSysErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb900GeVSysErr->Clone("graphInvCrossSectionPi0Comb900GeVSysErrUnscaled");
	
	graphConversionXSectionPi0900GeVUnscaled = (TGraphAsymmErrors*) graphConversionXSectionPi0900GeV->Clone();
	graphInvCrossSectionSysPi0900GeVUnscaled = (TGraphAsymmErrors*) graphInvCrossSectionSysPi0900GeV->Clone();
	graphPHOSXSectionPi0900GeVUnscaled= (TGraphAsymmErrors*) graphPHOSXSectionPi0900GeV->Clone();
	graphSysErrInvCrossSectionPi0PHOS900GeVUnscaled = (TGraphAsymmErrors*) graphSysErrInvCrossSectionPi0PHOS900GeV->Clone();
	
	forOutput= WriteParameterToFile(fitInvCrossSectionPi0900GeV);
	fileFinalResults<< forOutput.Data()<< endl;
	
	graphRatioCombNLOPi0900GeVMuHalf= (TGraph*)graphNLOMuHalfPi0900GeV->Clone();
	graphRatioCombNLOPi0900GeVMuOne= (TGraph*)graphNLOMuOnePi0900GeV->Clone();
	graphRatioCombNLOPi0900GeVMuTwo= (TGraph*)graphNLOMuTwoPi0900GeV->Clone();
	graphRatioCombNLOBKKPi0900GeVMuTwo= (TGraph*)graphNLOBKKCalcMuTwo900GeV->Clone();
	
	graphRatioCombPHOSPi0900GeV = (TGraphAsymmErrors*) graphPHOSXSectionPi0900GeV->Clone();	
	graphRatioCombPHOSPi0900GeVSys = (TGraphAsymmErrors*) graphSysErrInvCrossSectionPi0PHOS900GeV->Clone();	
	graphRatioCombConvPi0900GeV = (TGraphAsymmErrors*) graphConversionXSectionPi0900GeV->Clone();		
	graphRatioCombConvPi0900GeVSys = (TGraphAsymmErrors*) graphInvCrossSectionSysPi0900GeV->Clone();		
	
	graphRatioCombPHOSPi0900GeV = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0900GeV, fitInvCrossSectionPi0900GeV); 
	graphRatioCombPHOSPi0900GeVSys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0900GeVSys, fitInvCrossSectionPi0900GeV); 
	graphRatioCombConvPi0900GeV = CalculateGraphErrRatioToFit (graphRatioCombConvPi0900GeV, fitInvCrossSectionPi0900GeV); 
	graphRatioCombConvPi0900GeVSys = CalculateGraphErrRatioToFit (graphRatioCombConvPi0900GeVSys, fitInvCrossSectionPi0900GeV); 
	graphRatioCombNLOPi0900GeVMuHalf = CalculateGraphRatioToFit (graphRatioCombNLOPi0900GeVMuHalf, fitInvCrossSectionPi0900GeV); 
	graphRatioCombNLOPi0900GeVMuOne = CalculateGraphRatioToFit (graphRatioCombNLOPi0900GeVMuOne, fitInvCrossSectionPi0900GeV); 
	graphRatioCombNLOPi0900GeVMuTwo = CalculateGraphRatioToFit (graphRatioCombNLOPi0900GeVMuTwo, fitInvCrossSectionPi0900GeV); 
	graphRatioCombNLOBKKPi0900GeVMuTwo = CalculateGraphRatioToFit (graphRatioCombNLOBKKPi0900GeVMuTwo, fitInvCrossSectionPi0900GeV); 
	
	graphRatioCombCombFit900GeV = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb900GeV->Clone();
	graphRatioCombCombFit900GeV = CalculateGraphErrRatioToFit(graphRatioCombCombFit900GeV, fitInvCrossSectionPi0900GeV); 
	graphRatioCombCombFit900GeVStat = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb900GeVStatErr->Clone();
	graphRatioCombCombFit900GeVStat = CalculateGraphErrRatioToFit(graphRatioCombCombFit900GeVStat, fitInvCrossSectionPi0900GeV); 
	graphRatioCombCombFit900GeVSys = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb900GeVSysErr->Clone();
	graphRatioCombCombFit900GeVSys = CalculateGraphErrRatioToFit(graphRatioCombCombFit900GeVSys, fitInvCrossSectionPi0900GeV); 
		
	histoFitInvCrossSectionPi0900GeV = (TH1D*)fitInvCrossSectionPi0900GeV->GetHistogram();
	histoFitInvCrossSectionPi0900GeV->Scale(1e-1);

	graphInvCrossSectionEtaComb900GeVStatErr = new TGraphAsymmErrors(histoInvCrossSectionEta900GeV);	
	graphInvCrossSectionEtaComb900GeVStatErr->RemovePoint(0);
	graphInvCrossSectionEtaComb900GeVSysErr =  (TGraphAsymmErrors*)graphInvCrossSectionSysEta900GeV->Clone("graphInvCrossSectionEtaComb900GeVSysErr");
	graphInvCrossSectionEtaComb900GeV = CalculateSysErrFromRelSysHistoComplete( histoInvCrossSectionEta900GeV, "graphInvCrossSectionEtaComb900GeV",relSystErrorEta900GeVDown, relSystErrorEta900GeVDown, 2, nPointsEta900GeV);
	graphInvCrossSectionEtaComb900GeV->SetMarkerStyle(markerStyleConv);
	graphInvCrossSectionEtaComb900GeV->SetMarkerColor(1);
	graphInvCrossSectionEtaComb900GeV->SetLineColor(1);
	graphInvCrossSectionEtaComb900GeV->SetLineWidth(5);
	graphInvCrossSectionEtaComb900GeV->SetFillColor(0);
	graphInvCrossSectionEtaComb900GeV->SetFillStyle(0);
	graphInvCrossSectionEtaComb900GeVUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb900GeV->Clone("graphInvCrossSectionEtaComb900GeVUnscaled");
	graphInvCrossSectionEtaComb900GeVStatErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb900GeVStatErr->Clone("graphInvCrossSectionEtaComb900GeVStatErrUnscaled");
	graphInvCrossSectionEtaComb900GeVSysErrUnscaled = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb900GeVSysErr->Clone("graphInvCrossSectionEtaComb900GeVSysErrUnscaled");
	fitInvCrossSectionEta900GeV = FitObject("l","fitInvCrossSectionEta900GeV","Eta",graphInvCrossSectionEtaComb900GeV,0.9,3.,paramGraphEta,"NQRME+","fixn");
	cout << WriteParameterToFile(fitInvCrossSectionEta900GeV)<< endl;
	forOutput= WriteParameterToFile(fitInvCrossSectionEta900GeV);
	fileFinalResults<< forOutput.Data()<< endl;
	
	graphRatioCombNLOEta900GeVMuHalf= (TGraph*)graphNLOMuHalfEta900GeV->Clone();
	graphRatioCombNLOEta900GeVMuOne= (TGraph*)graphNLOMuOneEta900GeV->Clone();
	graphRatioCombNLOEta900GeVMuTwo= (TGraph*)graphNLOMuTwoEta900GeV->Clone();
	histoRatioCombConvEta900GeV = (TH1D*) histoInvCrossSectionEta900GeV->Clone();		
	histoRatioCombConvEta900GeV = CalculateHistoRatioToFit (histoRatioCombConvEta900GeV, fitInvCrossSectionEta900GeV); 
	graphRatioCombNLOEta900GeVMuHalf = CalculateGraphRatioToFit (graphRatioCombNLOEta900GeVMuHalf, fitInvCrossSectionEta900GeV); 
	graphRatioCombNLOEta900GeVMuOne = CalculateGraphRatioToFit (graphRatioCombNLOEta900GeVMuOne, fitInvCrossSectionEta900GeV); 
	graphRatioCombNLOEta900GeVMuTwo = CalculateGraphRatioToFit (graphRatioCombNLOEta900GeVMuTwo, fitInvCrossSectionEta900GeV); 
	graphRatioCombCombFitEta900GeV = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb900GeV->Clone();
	graphRatioCombCombFitEta900GeV = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta900GeV, fitInvCrossSectionEta900GeV); 
	graphRatioCombCombFitEta900GeVStat = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb900GeVStatErr->Clone();
	graphRatioCombCombFitEta900GeVStat = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta900GeVStat, fitInvCrossSectionEta900GeV); 
	graphRatioCombCombFitEta900GeVSys = (TGraphAsymmErrors*)graphInvCrossSectionEtaComb900GeVSysErr->Clone();
	graphRatioCombCombFitEta900GeVSys = CalculateGraphErrRatioToFit(graphRatioCombCombFitEta900GeVSys, fitInvCrossSectionEta900GeV); 
	
	histoFitInvCrossSectionEta900GeV = (TH1D*)fitInvCrossSectionEta900GeV->GetHistogram();

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
	
	//****************************************** Inv Cross Section **********************************************************************
	TCanvas* canvasInvXSection = new TCanvas("canvasInvXSection","",200,10,1200,2000);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSection,  0.15, 0.02, 0.03, 0.06);

	TPad* padComparisonXSection = new TPad("padComparisonXSection", "", 0., 0.42, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSection, 0.15, 0.02, 0.03, 0.);
	padComparisonXSection->Draw();
	
	TPad* padXSectionRatioPi07TeV = new TPad("padXSectionRatioPi07TeV", "", 0., 0.30, 1., 0.42,-1, -1, -2);
	DrawGammaPadSettings( padXSectionRatioPi07TeV, 0.15, 0.02, 0., 0.);
	padXSectionRatioPi07TeV->Draw();
	
	TPad* padXSectionRatioPi0900GeV = new TPad("padXSectionRatioPi0900GeV", "", 0., 0.18, 1., 0.30,-1, -1, -2);
	DrawGammaPadSettings( padXSectionRatioPi0900GeV, 0.15, 0.02, 0., 0.);
	padXSectionRatioPi0900GeV->Draw();
	
	TPad* padXSectionRatioEta7TeV = new TPad("padXSectionRatioEta7TeV", "", 0., 0., 1., 0.18,-1, -1, -2);
	DrawGammaPadSettings( padXSectionRatioEta7TeV,  0.15, 0.02, 0., 0.28);
	padXSectionRatioEta7TeV->Draw();
			
	padComparisonXSection->cd();
	padComparisonXSection->SetLogy();		
	padComparisonXSection->SetLogx();		

	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSection;
	histo2DInvXSection = new TH2F("histo2DInvXSection","histo2DInvXSection",1000,0.23,30.,1000,2e-4,10e12);
	SetStyleHistoTH2ForGraphs(histo2DInvXSection, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSection->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphSysErrInvCrossSectionPi0PHOS, markerStylePHOS,markerSizeSpectrum, colorPHOSPi07TeV , colorPHOSPi07TeV, widthLinesBoxes, kTRUE);
//	graphSysErrInvCrossSectionPi0PHOS->Draw("same,p");
	graphSysErrInvCrossSectionPi0PHOS->Draw("E2same");
	
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionSysPi0,markerStyleConv,markerSizeSpectrum, colorConvPi07TeV , colorConvPi07TeV, widthLinesBoxes, kTRUE);
//	graphInvCrossSectionSysPi0->Draw("p,same");
	graphInvCrossSectionSysPi0->Draw("E2same");
	
	graphSysErrInvCrossSectionPi0PHOS900GeV = ScaleGraph(graphSysErrInvCrossSectionPi0PHOS900GeV,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphSysErrInvCrossSectionPi0PHOS900GeV,markerStylePHOS,markerSizeSpectrum, colorPHOSPi0900GeV , colorPHOSPi0900GeV, widthLinesBoxes, kTRUE);
//	graphSysErrInvCrossSectionPi0PHOS900GeV->Draw("same,p");
	graphSysErrInvCrossSectionPi0PHOS900GeV->Draw("E2same");
	
	graphInvCrossSectionSysPi0900GeV = ScaleGraph(graphInvCrossSectionSysPi0900GeV,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionSysPi0900GeV,markerStyleConv,markerSizeSpectrum, colorConvPi0900GeV, colorConvPi0900GeV, widthLinesBoxes, kTRUE);
	//graphInvCrossSectionSysPi0900GeV->Draw("same,p");
	graphInvCrossSectionSysPi0900GeV->Draw("E2same");
	
	graphSysErrInvCrossSectionEtaPHOS = ScaleGraph(graphSysErrInvCrossSectionEtaPHOS,1e-3);
	DrawGammaSetMarkerTGraphAsym(graphSysErrInvCrossSectionEtaPHOS,markerStylePHOS,markerSizeSpectrum, colorPHOSEta7TeV , colorPHOSEta7TeV, widthLinesBoxes, kTRUE);
	//graphSysErrInvCrossSectionEtaPHOS->Draw("same,p");
	graphSysErrInvCrossSectionEtaPHOS->Draw("E2same");
	
	graphInvCrossSectionSysEta = ScaleGraph(graphInvCrossSectionSysEta,1e-3);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionSysEta,markerStyleConv,markerSizeSpectrum, colorConvEta7TeV, colorConvEta7TeV, widthLinesBoxes, kTRUE);
	//graphInvCrossSectionSysEta->Draw("p,same");
	graphInvCrossSectionSysEta->Draw("E2same");
		
	//-------------------------- Drawing Spectrum
	DrawGammaSetMarkerTGraphAsym(graphPHOSXSectionPi07TeV,markerStylePHOS,markerSizeSpectrum, colorPHOSPi07TeV, colorPHOSPi07TeV);
	graphPHOSXSectionPi07TeV->Draw("p,same,e1");

	DrawGammaSetMarkerTGraphAsym(graphConversionXSectionPi07TeV ,markerStyleConv,markerSizeSpectrum, colorConvPi07TeV, colorConvPi07TeV);
	graphConversionXSectionPi07TeV->Draw("p,same,e1");
	
// 	DrawGammaSetMarker(histoInvCrossSectionPi0EMCAL, markerStyleEMCAL, markerSizeSpectrum, colorEMCALPi07TeV,colorEMCALPi07TeV);
//	histoInvCrossSectionPi0EMCAL->Draw("same,e1");
	
	graphConversionXSectionPi0900GeV = ScaleGraph(graphConversionXSectionPi0900GeV,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphConversionXSectionPi0900GeV ,markerStyleConv,markerSizeSpectrum, colorConvPi0900GeV, colorConvPi0900GeV);
	graphConversionXSectionPi0900GeV->Draw("p,same,e1");
	
	graphPHOSXSectionPi0900GeV= ScaleGraph(graphPHOSXSectionPi0900GeV,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphPHOSXSectionPi0900GeV,markerStylePHOS,markerSizeSpectrum, colorPHOSPi0900GeV , colorPHOSPi0900GeV);
	graphPHOSXSectionPi0900GeV->Draw("p,same,e1");

	graphPHOSXSectionEta7TeV = ScaleGraph(graphPHOSXSectionEta7TeV,1e-3);
	DrawGammaSetMarkerTGraphAsym(graphPHOSXSectionEta7TeV,markerStylePHOS,markerSizeSpectrum, colorPHOSEta7TeV , colorPHOSEta7TeV);
	graphPHOSXSectionEta7TeV->Draw("p,same,e1");
	
	graphConversionXSectionEta7TeV = ScaleGraph(graphConversionXSectionEta7TeV,1e-3);
	DrawGammaSetMarkerTGraphAsym(graphConversionXSectionEta7TeV , markerStyleConv,markerSizeSpectrum, colorConvEta7TeV, colorConvEta7TeV);
	graphConversionXSectionEta7TeV->Draw("p,same,e1");
	
/*	graphInvCrossSectionEtaEMCAL = ScaleGraph(graphInvCrossSectionEtaEMCAL,1e-3);
	DrawGammaSetMarkerTGraphErr(graphInvCrossSectionEtaEMCAL, markerStyleEMCAL, markerSizeSpectrum, colorEMCALEta7TeV,colorEMCALEta7TeV,widthStatErrBars);
	graphInvCrossSectionEtaEMCAL->Draw("same,ep");*/
	
	DrawGammaSetMarkerTF1( fitInvCrossSectionPi0, styleFitCommonSpectrum, widthCommonFit, colorCommonSpectrumPi07TeV);
	fitInvCrossSectionPi0->Draw("same");
	
	SetStyleHisto(histoFitInvCrossSectionPi0900GeV, widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi0900GeV);
	histoFitInvCrossSectionPi0900GeV->Draw("same,c");
	
	SetStyleHisto(histoFitInvCrossSectionEta7TeV, widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumEta7TeV);
	histoFitInvCrossSectionEta7TeV->Draw("same,c");
	
	TLatex *labelScalingPi07TeV = new TLatex(0.27,3E11,"x 1");
	SetStyleTLatex( labelScalingPi07TeV, 0.025,4,fitInvCrossSectionPi0->GetLineColor(),62,kFALSE);
	labelScalingPi07TeV->Draw();
	
	TLatex *labelScalingEta7TeV = new TLatex(0.27,1E7,"x 10^{-3}");
	SetStyleTLatex( labelScalingEta7TeV, 0.025,4,histoFitInvCrossSectionEta7TeV->GetLineColor(),62,kFALSE);
	labelScalingEta7TeV->Draw();
	
	TLatex *labelScalingPi0900GeV = new TLatex(0.27,5E9,"x 10^{-1}");
	SetStyleTLatex( labelScalingPi0900GeV, 0.025,4,histoFitInvCrossSectionPi0900GeV->GetLineColor(),62,kFALSE);
	labelScalingPi0900GeV->Draw();
	DrawNormalizationErrorText(normalizationInvX3EnPi0Eta[0],normalizationInvX3EnPi0Eta[1],normalizationInvX3EnPi0Eta[2],normalizationInvX3EnPi0Eta[3],normalizationInvX3EnPi0Eta[4],"No2.76"); 

	//********************************** Defintion of the Legend **************************************************	
		Double_t columnsLegend[4] 	= {0.,0.18,0.47,0.75};
		Double_t rowsLegend[6] 		= {0.88,0.75,0.57,0.4,0.22,0.05}; //with EMCAL {0.88,0.75,0.57,0.4,0.22,0.05};
		//******************* Text sizes *******************
		Size_t textSizeLeftColumn	= 0.13;
		Size_t textSizeTopRow		= 0.13; 
		Size_t textSizeSecondRow 	= 0.11;
		//******************* Offsets ***********************
		Double_t offsetSystColumn 	= 0.15;
		Double_t offsetMarkerX		= 0.04;
		Double_t offsetMarkerY		= 0.05;
		Double_t offsetBoxSizeY		= 0.05;
		Double_t offsetFit			= 0.04;
		//****************** Scale factors ******************
		Double_t scaleWidthLine 		= 0.8;
		
		TPad* padXSectionLegend = new TPad("padXSectionLegend", "", 0.17, 0.005, 0.95, 0.15,-1, -1, -2);  
		//TPad* padXSectionLegend = new TPad("padXSectionLegend", "", 0.17, 0.02, 0.95, 0.18,-1, -1, -2);
		DrawGammaPadSettings( padXSectionLegend, 0., 0., 0., 0.);
		padXSectionLegend->Draw();
		padXSectionLegend->cd();

		//****************** first Column **************************************************
		TLatex *textCTS = new TLatex(columnsLegend[0],rowsLegend[2],"PCM");
		SetStyleTLatex( textCTS, textSizeLeftColumn,4);
		textCTS->Draw();
		TLatex *textPHOS = new TLatex(columnsLegend[0],rowsLegend[3],"PHOS");
		SetStyleTLatex( textPHOS, textSizeLeftColumn,4);
		textPHOS->Draw();
		TLatex *textEMCAL = new TLatex(columnsLegend[0],rowsLegend[4],"EMCAL");
		SetStyleTLatex( textEMCAL, textSizeLeftColumn,4);
// 		textEMCAL->Draw();
		TLatex *textFit = new TLatex(columnsLegend[0] ,rowsLegend[4] ,"Fit");
		SetStyleTLatex( textFit, textSizeLeftColumn,4);
		textFit->Draw();

		//****************** second Column *************************************************
		TLatex *textPi07TeV = new TLatex(columnsLegend[1],rowsLegend[0],"#pi^{0}, #sqrt{#it{s}} = 7 TeV");
		SetStyleTLatex( textPi07TeV, textSizeTopRow,4);
		textPi07TeV->Draw();
		TLatex *textPi07TeVstat = new TLatex(columnsLegend[1],rowsLegend[1] ,"stat");
		SetStyleTLatex( textPi07TeVstat, textSizeSecondRow,4);
		textPi07TeVstat->Draw();
		TLatex *textPi07TeVsys = new TLatex(columnsLegend[1]+ offsetSystColumn ,rowsLegend[1],"syst");
		SetStyleTLatex( textPi07TeVsys, textSizeSecondRow,4);
		textPi07TeVsys->Draw();
		TMarker* markerCTSPi07TeV = CreateMarkerFromGraph(graphConversionXSectionPi07TeV,columnsLegend[1]+ offsetMarkerX ,rowsLegend[2]+ offsetMarkerY ,scaleWidthLine);
		markerCTSPi07TeV->DrawMarker(columnsLegend[1]+ offsetMarkerX ,rowsLegend[2]+ offsetMarkerY);
		TMarker* markerPHOSPi07TeV = CreateMarkerFromGraph(graphPHOSXSectionPi07TeV, columnsLegend[1]+ offsetMarkerX ,rowsLegend[3]+ offsetMarkerY ,scaleWidthLine);
		markerPHOSPi07TeV->DrawMarker(columnsLegend[1]+ offsetMarkerX ,rowsLegend[3]+ offsetMarkerY );
// 		TMarker* markerEMCALPi07TeV = CreateMarkerFromHisto(histoInvCrossSectionPi0EMCAL, columnsLegend[1]+ offsetMarkerX ,rowsLegend[4]+ offsetMarkerY ,scaleWidthLine);
// 		markerEMCALPi07TeV->DrawMarker(columnsLegend[1]+ offsetMarkerX ,rowsLegend[4]+ offsetMarkerY );
		TLine * lineFit7TeV = CreateLineFromFit(fitInvCrossSectionPi0, columnsLegend[1]+ offsetMarkerX- offsetFit, rowsLegend[4]+offsetMarkerY, columnsLegend[1]+ offsetMarkerX+ offsetFit, rowsLegend[4]+offsetMarkerY, scaleWidthLine);
		lineFit7TeV->Draw("same");	
		TBox* boxCTSPi07TeV = CreateBoxFromGraph(graphPHOSXSectionPi07TeV,columnsLegend[1]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend[2]+ offsetMarkerY- offsetBoxSizeY, columnsLegend[1]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend[2]+ offsetMarkerY+ offsetBoxSizeY);
		boxCTSPi07TeV->Draw("l");
		TBox* boxPHOSPi07TeV = CreateBoxFromGraph(graphPHOSXSectionPi07TeV,columnsLegend[1]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend[3]+ offsetMarkerY- offsetBoxSizeY, columnsLegend[1]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend[3]+ offsetMarkerY+ offsetBoxSizeY);
		boxPHOSPi07TeV->Draw("l");
/*		TBox* boxEMCALPi07TeV = CreateBoxFromHisto(histoInvCrossSectionPi0EMCAL,columnsLegend[1]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend[4]+ offsetMarkerY- offsetBoxSizeY, columnsLegend[1]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend[4]+ offsetMarkerY+ offsetBoxSizeY);
		boxEMCALPi07TeV->Draw("l");*/
			
		//***************** third Column ***************************************************
		TLatex *textPi0900GeV = new TLatex(columnsLegend[2],rowsLegend[0],"#pi^{0}, #sqrt{#it{s}} = 0.9 TeV");
		SetStyleTLatex( textPi0900GeV, textSizeTopRow,4);
		textPi0900GeV->Draw();
		TLatex *textPi0900GeVstat = new TLatex(columnsLegend[2],rowsLegend[1],"stat");
		SetStyleTLatex( textPi0900GeVstat, textSizeSecondRow,4);
		textPi0900GeVstat->Draw();
		TLatex *textPi0900GeVsys = new TLatex(columnsLegend[2]+ offsetSystColumn,rowsLegend[1],"syst");
		SetStyleTLatex( textPi0900GeVsys, textSizeSecondRow,4);
		textPi0900GeVsys->Draw();
		TMarker* markerCTSPi0900GeV = CreateMarkerFromGraph(graphConversionXSectionPi0900GeV,columnsLegend[2]+ offsetMarkerX ,rowsLegend[2]+ offsetMarkerY ,scaleWidthLine);
		markerCTSPi0900GeV->DrawMarker(columnsLegend[2]+ offsetMarkerX ,rowsLegend[2]+ offsetMarkerY);
		TMarker* markerPHOSPi0900GeV = CreateMarkerFromGraph(graphPHOSXSectionPi0900GeV, columnsLegend[2]+ offsetMarkerX ,rowsLegend[3]+ offsetMarkerY ,scaleWidthLine);
		markerPHOSPi0900GeV->DrawMarker(columnsLegend[2]+ offsetMarkerX ,rowsLegend[3]+ offsetMarkerY);
/*		TMarker* markerEMCALPi0900GeV = CreateMarkerFromGraph(graphInvCrossSectionPi0EMCAL, columnsLegend[2]+ offsetMarkerX ,rowsLegend[4]+ offsetMarkerY ,scaleWidthLine);
		markerEMCALPi0900GeV->DrawMarker(columnsLegend[2]+ offsetMarkerX ,rowsLegend[4]+ offsetMarkerY);*/
		TLine * lineFit900GeV = CreateLineFromHisto(histoFitInvCrossSectionPi0900GeV, columnsLegend[2]+ offsetMarkerX- offsetFit, rowsLegend[4]+offsetMarkerY, columnsLegend[2]+ offsetMarkerX+ offsetFit, rowsLegend[4]+offsetMarkerY, scaleWidthLine);
		lineFit900GeV->Draw("same");		
		TBox* boxCTSPi0900GeV = CreateBoxFromGraph(graphConversionXSectionPi0900GeV,columnsLegend[2]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend[2]+ offsetMarkerY- offsetBoxSizeY, columnsLegend[2]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend[2]+ offsetMarkerY+ offsetBoxSizeY);
		boxCTSPi0900GeV->Draw("l");		
		TBox* boxPHOSPi0900GeV = CreateBoxFromGraph(graphPHOSXSectionPi0900GeV,columnsLegend[2]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend[3]+ offsetMarkerY- offsetBoxSizeY, columnsLegend[2]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend[3]+ offsetMarkerY+ offsetBoxSizeY);
		boxPHOSPi0900GeV->Draw("l");
// 		TBox* boxEMCALPi0900GeV =  CreateBoxFromGraph(graphInvCrossSectionPi0EMCAL,columnsLegend[2]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend[4]+ offsetMarkerY- offsetBoxSizeY, columnsLegend[2]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend[4]+ offsetMarkerY+ offsetBoxSizeY);
// 		boxEMCALPi0900GeV->Draw("l");

		//***************** forth Column ***************************************************
		TLatex *textEta7TeV = new TLatex(columnsLegend[3],rowsLegend[0],"#eta, #sqrt{#it{s}} = 7 TeV");
		SetStyleTLatex( textEta7TeV, textSizeTopRow,4);
		textEta7TeV->Draw();
		TLatex *textEta7TeVstat = new TLatex(columnsLegend[3],rowsLegend[1],"stat");
		SetStyleTLatex( textEta7TeVstat, textSizeSecondRow,4);
		textEta7TeVstat->Draw();
		TLatex *textEta7TeVsys = new TLatex(columnsLegend[3]+ offsetSystColumn ,rowsLegend[1],"syst");
		SetStyleTLatex( textEta7TeVsys, textSizeSecondRow,4);
		textEta7TeVsys->Draw();
		TMarker* markerCTSEta7TeV = CreateMarkerFromGraph(graphConversionXSectionEta7TeV,columnsLegend[3]+ offsetMarkerX ,rowsLegend[2]+ offsetMarkerY ,scaleWidthLine);
		markerCTSEta7TeV->DrawMarker(columnsLegend[3]+ offsetMarkerX ,rowsLegend[2]+ offsetMarkerY);
		TMarker* markerPHOSEta7TeV = CreateMarkerFromGraph(graphPHOSXSectionEta7TeV, columnsLegend[3]+ offsetMarkerX ,rowsLegend[3]+ offsetMarkerY ,scaleWidthLine);
		markerPHOSEta7TeV->DrawMarker(columnsLegend[3]+ offsetMarkerX ,rowsLegend[3]+ offsetMarkerY );
/*		TMarker* markerEMCALEta7TeV = CreateMarkerFromGraph(graphInvCrossSectionEtaEMCAL, columnsLegend[3]+ offsetMarkerX ,rowsLegend[4]+ offsetMarkerY ,scaleWidthLine);
		markerEMCALEta7TeV->DrawMarker(columnsLegend[3]+ offsetMarkerX ,rowsLegend[4]+ offsetMarkerY);*/
		TLine * lineFitEta7TeV = CreateLineFromHisto(histoFitInvCrossSectionEta7TeV, columnsLegend[3]+ offsetMarkerX- offsetFit, rowsLegend[4]+offsetMarkerY, columnsLegend[3]+ offsetMarkerX+ offsetFit, rowsLegend[4]+offsetMarkerY, scaleWidthLine);
		lineFitEta7TeV->Draw("same");
		TBox* boxCTSEta7TeV = CreateBoxFromGraph(graphConversionXSectionEta7TeV,columnsLegend[3]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend[2]+ offsetMarkerY- offsetBoxSizeY, columnsLegend[3]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend[2]+ offsetMarkerY+ offsetBoxSizeY);
		boxCTSEta7TeV->Draw("l");
		TBox* boxPHOSEta7TeV = CreateBoxFromGraph(graphPHOSXSectionEta7TeV,columnsLegend[3]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend[3]+ offsetMarkerY- offsetBoxSizeY, columnsLegend[3]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend[3]+ offsetMarkerY+ offsetBoxSizeY);
		boxPHOSEta7TeV->Draw("l");
		/*	TBox* boxEMCALEta7TeV = CreateBoxFromGraph(graphInvCrossSectionPi0EMCAL,columnsLegend[3]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend[4]+ offsetMarkerY- offsetBoxSizeY, columnsLegend[3]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend[4]+ offsetMarkerY+ offsetBoxSizeY);
		boxEMCALEta7TeV->Draw("l");*/

	//************************************* End Legend Definition ***************************************************************

	padXSectionRatioPi07TeV->cd();
	padXSectionRatioPi07TeV->SetLogx();
	TH2F * ratio2DInvXSectionPi0;
	ratio2DInvXSectionPi0 = new TH2F("ratio2DInvXSectionPi0","ratio2DInvXSectionPi0",1000,0.23,30.,1000,0.4,1.72);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionPi0->DrawCopy(); 

	TLatex *labelRatioNLOPi07TeV = new TLatex(0.18,0.75,"#pi^{0}, #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( labelRatioNLOPi07TeV, 0.17,4);
	
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0Sys, markerStylePHOS,markerSizeSpectrum, colorPHOSPi07TeV , colorPHOSPi07TeV, widthLinesBoxes, kTRUE);
//	graphRatioCombPHOSPi0Sys->Draw("same,p");
	graphRatioCombPHOSPi0Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0Sys, markerStyleConv,markerSizeSpectrum, colorConvPi07TeV , colorConvPi07TeV, widthLinesBoxes, kTRUE);
//	graphRatioCombConvPi0Sys->Draw("p,same");
	graphRatioCombConvPi0Sys->Draw("E2same");
// 	DrawGammaSetMarker(histoRatioCombEMCALPi0, markerStyleEMCAL, markerSizeSpectrum, colorEMCALPi07TeV,colorEMCALPi07TeV);
// 	histoRatioCombEMCALPi0->Draw("same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0, markerStylePHOS,markerSizeSpectrum, colorPHOSPi07TeV , colorPHOSPi07TeV);
	graphRatioCombPHOSPi0->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0 , markerStyleConv,markerSizeSpectrum, colorConvPi07TeV , colorConvPi07TeV);
	graphRatioCombConvPi0->Draw("p,same,e1");
	TBox* boxErrorSigmaPi07TeVRatio = CreateBoxConv(colorCommonSpectrumPi07TeVBox, 0.25, 1.-(xSection7TeVErrDown /(xSection7TeV*1e3) ), 0.29, 1.+(xSection7TeVErrUp /(xSection7TeV*1e3)));
	boxErrorSigmaPi07TeVRatio->Draw();
	
	labelRatioNLOPi07TeV->Draw();
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionRatioPi0900GeV->cd();
	padXSectionRatioPi0900GeV->SetLogx();
	TH2F * ratio2DInvXSectionPi0900GeV = new TH2F("ratio2DInvXSectionPi0900GeV","ratio2DInvXSectionPi0900GeV",1000,0.23,30.,1000,0.4,1.72);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionPi0900GeV->DrawCopy(); 

	TLatex *labelRatioNLOPi0900GeV = new TLatex(0.18,0.75,"#pi^{0}, #sqrt{#it{s}} = 0.9 TeV");
	SetStyleTLatex( labelRatioNLOPi0900GeV, 0.17,4);
	
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0900GeVSys,  markerStylePHOS,markerSizeSpectrum, colorPHOSPi0900GeV , colorPHOSPi0900GeV, widthLinesBoxes, kTRUE);
//	graphRatioCombPHOSPi0900GeVSys->Draw("same,p");
	graphRatioCombPHOSPi0900GeVSys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0900GeVSys, markerStyleConv,markerSizeSpectrum, colorConvPi0900GeV, colorConvPi0900GeV,widthLinesBoxes, kTRUE);
//	graphRatioCombConvPi0900GeVSys->Draw("same,p");
	graphRatioCombConvPi0900GeVSys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0900GeV ,  markerStyleConv,markerSizeSpectrum, colorConvPi0900GeV, colorConvPi0900GeV);
	graphRatioCombConvPi0900GeV->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0900GeV, markerStylePHOS,markerSizeSpectrum, colorPHOSPi0900GeV , colorPHOSPi0900GeV);
	graphRatioCombPHOSPi0900GeV->Draw("p,same,e1");
	TBox* boxErrorSigmaPi0900GeVRatio = CreateBoxConv(colorCommonSpectrumPi0900GeVBox, 0.25, 1.-(xSection900GeVErrDown/(xSection900GeV*1e3) ), 0.29, 1.+(xSection900GeVErrUp/(xSection900GeV*1e3)));
	boxErrorSigmaPi0900GeVRatio->Draw();
	labelRatioNLOPi0900GeV->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionRatioEta7TeV->cd();
	padXSectionRatioEta7TeV->SetLogx();
	
	TH2F * ratio2DInvXSectionEta= new TH2F("ratio2DInvXSectionEta","ratio2DInvXSectionEta",1000,0.23,30.,1000,0.4,1.72);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEta, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Fit}", 0.12,0.13, 0.125,0.128, 0.9,0.45, 510, 505);	
	ratio2DInvXSectionEta->DrawCopy(); 
	
	TLatex *labelRatioNLOEta7TeV = new TLatex(0.18,0.85,"#eta, #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( labelRatioNLOEta7TeV, 0.115,4);
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSEta7TeVSys, markerStylePHOS,markerSizeSpectrum, colorPHOSEta7TeV, colorPHOSEta7TeV, widthLinesBoxes, kTRUE);
//	graphRatioCombPHOSEta7TeVSys->Draw("same,p");
	graphRatioCombPHOSEta7TeVSys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvEta7TeVSys, markerStyleConv,markerSizeSpectrum, colorConvEta7TeV, colorConvEta7TeV,widthLinesBoxes, kTRUE);
//	graphRatioCombConvEta7TeVSys->Draw("p,same");
	graphRatioCombConvEta7TeVSys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSEta7TeV,  markerStylePHOS,markerSizeSpectrum, colorPHOSEta7TeV, colorPHOSEta7TeV);
	graphRatioCombPHOSEta7TeV->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvEta7TeV , markerStyleConv,markerSizeSpectrum, colorConvEta7TeV, colorConvEta7TeV);
	graphRatioCombConvEta7TeV->Draw("p,same,e1");
/*	DrawGammaSetMarkerTGraphErr(graphRatioCombEMCALEta7TeV, markerStyleEMCAL, markerSizeSpectrum, colorEMCALEta7TeV,colorEMCALEta7TeV,widthStatErrBars);
	graphRatioCombEMCALEta7TeV->Draw("same,ep");*/
	TBox* boxErrorSigmaEta7TeVRatio = CreateBoxConv(colorCommonSpectrumEta7TeVBox,0.25, 1.-(xSection7TeVErrDown/(xSection7TeV*1e3) ), 0.29, 1.+(xSection7TeVErrUp/(xSection7TeV*1e3) ));
	boxErrorSigmaEta7TeVRatio->Draw();
	labelRatioNLOEta7TeV->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	canvasInvXSection->Update();

	canvasInvXSection->Print(Form("%s/InvXSection_Paper.%s",outputDir.Data(),suffix.Data()));

	padComparisonXSection->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSPi0EtaPrelim(pictDrawingCoordinatesCombineMeasCross[0], pictDrawingCoordinatesCombineMeasCross[1], pictDrawingCoordinatesCombineMeasCross[2], pictDrawingCoordinatesCombineMeasCross[3], pictDrawingCoordinatesCombineMeasCross[4], pictDrawingCoordinatesCombineMeasCross[5], pictDrawingCoordinatesCombineMeasCross[6], pictDrawingCoordinatesCombineMeasCross[7], pictDrawingCoordinates[8],collisionSystemCombined, pictDrawingOptions[1], pictDrawingOptions[2], kTRUE,1200,1220);
	canvasInvXSection->Update();
	
	canvasInvXSection->Print(Form("%s/InvXSection.%s",outputDir.Data(),suffix.Data()));
	
	//**************************************************************************************************************************************
	//************************************ Inv Crosssection only ratio Pi0 ***************************************************************
	//**************************************************************************************************************************************
	
	TCanvas* canvasInvXSectionOnlyRatioPi0 = new TCanvas("canvasInvXSectionOnlyRatioPi0","",200,10,1200,700);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionOnlyRatioPi0,  0.08, 0.01, 0.015, 0.13);

	canvasInvXSectionOnlyRatioPi0->SetLogx();
	TH2F * ratio2DInvXSectionOnlyPi0;
	ratio2DInvXSectionOnlyPi0 = new TH2F("ratio2DInvXSectionOnlyPi0","ratio2DInvXSectionOnlyPi0",1000,0.23,30.,1000,0.25,2.15	);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionOnlyPi0, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.064, 0.05,0.064, 0.8,0.6, 512, 505);
	ratio2DInvXSectionOnlyPi0->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0Sys, markerStylePHOS,markerSizeSpectrum, colorPHOSMass , colorPHOSMass, widthLinesBoxes, kTRUE);
	graphRatioCombPHOSPi0Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0Sys, markerStyleConv,markerSizeSpectrum, colorConv , colorConv, widthLinesBoxes, kTRUE);
	graphRatioCombConvPi0Sys->Draw("E2same");
// 	DrawGammaSetMarker(histoRatioCombEMCALPi0, markerStyleEMCAL, markerSizeSpectrum, colorEMCAL,colorEMCAL);
// 	histoRatioCombEMCALPi0->Draw("same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0, markerStylePHOS,markerSizeSpectrum*2., colorPHOSMass , colorPHOSMass);
	graphRatioCombPHOSPi0->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0 , markerStyleConv,markerSizeSpectrum*2., colorConv , colorConv);
	graphRatioCombConvPi0->Draw("p,same,e1");
	
	TLatex *labelRatioPi07TeV = new TLatex(0.13,0.88,"#pi^{0}, #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( labelRatioPi07TeV, 0.07,4);
	labelRatioPi07TeV->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	TPad* padXSectionLegendOnlyPi0Ratio = new TPad("padXSectionLegendOnlyPi0Ratio", "",   0.5, 0.8, 0.8, 0.95,-1, -1, -2);
	DrawGammaPadSettings( padXSectionLegendOnlyPi0Ratio, 0.0, 0.00, 0.0, 0.);
	padXSectionLegendOnlyPi0Ratio->Draw();
	
	padXSectionLegendOnlyPi0Ratio->cd();
	//****************************** Definition of the Legend ******************************************
	//**************** Row def ************************
	Double_t rowsLegendOnlyPi0Ratio[4] 		= {0.78,0.45,0.05,0.00};
	Double_t columnsLegendOnlyPi0Ratio[2] 	= {0.,0.3};
	//****************** first Column **************************************************
	TLatex *textCTSOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],"PCM");
	SetStyleTLatex( textCTSOnlyRatioPi0, textSizeLeftColumn*2.5,4);
	textCTSOnlyRatioPi0->Draw();
	TLatex *textPHOSOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[2],"PHOS");
	SetStyleTLatex( textPHOSOnlyRatioPi0, textSizeLeftColumn*2.5,4);
	textPHOSOnlyRatioPi0->Draw();
	TLatex *textEMCALOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],"EMCAL");
	SetStyleTLatex( textEMCALOnlyRatioPi0, textSizeLeftColumn*2.5,4);
// 	textEMCALOnlyRatioPi0->Draw();
	
	//****************** second Column *************************************************
	TLatex *textPi07TeVstatOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
	SetStyleTLatex( textPi07TeVstatOnlyRatioPi0, textSizeLeftColumn*2.5,4);
	textPi07TeVstatOnlyRatioPi0->Draw();
	TLatex *textPi07TeVsysOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[1]+ offsetSystColumn*1.5 ,rowsLegendOnlyPi0Ratio[0],"syst");
	SetStyleTLatex( textPi07TeVsysOnlyRatioPi0, textSizeLeftColumn*2.5,4);
	textPi07TeVsysOnlyRatioPi0->Draw();
	TMarker* markerCTSPi07TeVOnlyRatioPi0 = CreateMarkerFromGraph(graphRatioCombConvPi0,columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[1]+ offsetMarkerY ,scaleWidthLine);
	markerCTSPi07TeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[1]+ offsetMarkerY);
	TMarker* markerPHOSPi07TeVOnlyRatioPi0 = CreateMarkerFromGraph(graphRatioCombPHOSPi0, columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[2]+ offsetMarkerY ,scaleWidthLine);
	markerPHOSPi07TeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[2]+ offsetMarkerY );
// 	TMarker* markerEMCALPi07TeVOnlyRatioPi0 = CreateMarkerFromHisto(histoRatioCombEMCALPi0, columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[3]+ offsetMarkerY ,scaleWidthLine);
// 	markerEMCALPi07TeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[3]+ offsetMarkerY );
	TBox* boxCTSPi07TeVOnlyRatioPi0 = CreateBoxFromGraph(graphRatioCombConvPi0,columnsLegendOnlyPi0Ratio[1]+offsetSystColumn*1.5+offsetMarkerX-offsetFit, rowsLegendOnlyPi0Ratio[1]+ offsetMarkerY- offsetBoxSizeY, columnsLegendOnlyPi0Ratio[1]+ offsetSystColumn*1.5+ offsetMarkerX*2+ offsetFit*2, rowsLegendOnlyPi0Ratio[1]+ offsetMarkerY+ offsetBoxSizeY);
	boxCTSPi07TeVOnlyRatioPi0->Draw("l");
	TBox* boxPHOSPi07TeVOnlyRatioPi0 = CreateBoxFromGraph(graphRatioCombPHOSPi0,columnsLegendOnlyPi0Ratio[1]+offsetSystColumn*1.5+offsetMarkerX-offsetFit, rowsLegendOnlyPi0Ratio[2]+ offsetMarkerY- offsetBoxSizeY, columnsLegendOnlyPi0Ratio[1]+ offsetSystColumn*1.5+ offsetMarkerX*2+ offsetFit*2, rowsLegendOnlyPi0Ratio[2]+ offsetMarkerY+ offsetBoxSizeY);
	boxPHOSPi07TeVOnlyRatioPi0->Draw("l");
	/*	TBox* boxEMCALPi07TeVOnlyRatio = CreateBoxFromHisto(histoRatioCombEMCALPi0,columnsLegendOnlyPi0Ratio[1]+offsetSystColumn*1.5+offsetMarkerX-offsetFit, rowsLegendOnlyPi0Ratio[3]+ offsetMarkerY- offsetBoxSizeY, columnsLegendOnlyPi0Ratio[1]+ offsetSystColumn*1.5+ offsetMarkerX+ offsetFit, rowsLegendOnlyPi0Ratio[3]+ offsetMarkerY+ offsetBoxSizeY);
	boxEMCALPi07TeVOnlyRatio->Draw("l");*/
	
	
	canvasInvXSectionOnlyRatioPi0->Update();
	canvasInvXSectionOnlyRatioPi0->Print(Form("%s/InvXSection_OnlyRatioPi07TeV_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

	//********************************* Ratio 7 TeV Eta Only********************************************************************
	
		TCanvas* canvasInvXSectionOnlyRatioEta = new TCanvas("canvasInvXSectionOnlyRatioEta","",200,10,1200,700);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionOnlyRatioEta,  0.08, 0.01, 0.015, 0.13);

	canvasInvXSectionOnlyRatioEta->SetLogx();
	TH2F * ratio2DInvXSectionOnlyEta;
	ratio2DInvXSectionOnlyEta = new TH2F("ratio2DInvXSectionOnlyEta","ratio2DInvXSectionOnlyEta",1000,0.23,30.,1000,0.5,1.85);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionOnlyEta, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.064, 0.05,0.064, 0.8,0.6, 512, 505);
	ratio2DInvXSectionOnlyEta->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSEta7TeVSys, markerStylePHOS,markerSizeSpectrum*2, colorPHOSMass , colorPHOSMass, widthLinesBoxes, kTRUE);
//	graphSysErrRatioCombPHOSEta->Draw("same,p");
	graphRatioCombPHOSEta7TeVSys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvEta7TeVSys, markerStyleConv,markerSizeSpectrum*2, colorConv , colorConv, widthLinesBoxes, kTRUE);
//	graphSysErrRatioCombConvEta->Draw("p,same");
	graphRatioCombConvEta7TeVSys->Draw("E2same");
// 	histoRatioCombEMCALEta->Draw("same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSEta7TeV, markerStylePHOS,markerSizeSpectrum*2, colorPHOSMass , colorPHOSMass);
	graphRatioCombPHOSEta7TeV->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvEta7TeV , markerStyleConv,markerSizeSpectrum*2, colorConv , colorConv);
	graphRatioCombConvEta7TeV->Draw("p,same,e1");
	
	TLatex *labelRatioEta7TeV = new TLatex(0.13,0.88,"#eta, #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( labelRatioEta7TeV, 0.07,4);
	labelRatioEta7TeV->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	TPad* padXSectionLegendOnlyEtaRatio = new TPad("padXSectionLegendOnlyEtaRatio", "",  0.5, 0.8, 0.8, 0.95,-1, -1, -2);
	DrawGammaPadSettings( padXSectionLegendOnlyEtaRatio, 0.15, 0.02, 0.03, 0.);
	padXSectionLegendOnlyEtaRatio->Draw();
	
	padXSectionLegendOnlyEtaRatio->cd();
	//****************************** Definition of the Legend ******************************************
	//**************** Row def ************************
	textCTSOnlyRatioPi0->Draw();
	textPHOSOnlyRatioPi0->Draw();
	
	//****************** second Column *************************************************
	textPi07TeVstatOnlyRatioPi0->Draw();
	textPi07TeVsysOnlyRatioPi0->Draw();
	markerCTSPi07TeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[1]+ offsetMarkerY);
	markerPHOSPi07TeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[2]+ offsetMarkerY );
	boxCTSPi07TeVOnlyRatioPi0->Draw("l");
	boxPHOSPi07TeVOnlyRatioPi0->Draw("l");
	
	
	canvasInvXSectionOnlyRatioEta->Update();
	canvasInvXSectionOnlyRatioEta->Print(Form("%s/InvXSection_OnlyRatioEta7TeV_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

	//**************************************************************************************************************************************
	//************************************ Inv Crosssection only ratio Pi0 ***************************************************************
	//**************************************************************************************************************************************
	
	TCanvas* canvasInvXSectionOnlyRatioPi02760GeV = new TCanvas("canvasInvXSectionOnlyRatioPi02760GeV","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionOnlyRatioPi02760GeV,  0.12, 0.02, 0.02, 0.12);
	
	canvasInvXSectionOnlyRatioPi02760GeV->SetLogx();
	TH2F * ratio2DInvXSectionOnlyPi02760GeV;
	ratio2DInvXSectionOnlyPi02760GeV = new TH2F("ratio2DInvXSectionOnlyPi02760GeV","ratio2DInvXSectionOnlyPi02760GeV",1000,0.33,20.,1000,0.05,2.45);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionOnlyPi02760GeV, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.062, 0.05,0.062, 0.82,0.8, 512, 505);
	ratio2DInvXSectionOnlyPi02760GeV->GetXaxis()->SetLabelOffset(-0.015);
	ratio2DInvXSectionOnlyPi02760GeV->GetYaxis()->SetLabelOffset(0.01);
	ratio2DInvXSectionOnlyPi02760GeV->GetYaxis()->SetLabelFont(42);
	ratio2DInvXSectionOnlyPi02760GeV->GetXaxis()->SetLabelFont(42);
	ratio2DInvXSectionOnlyPi02760GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphSysErrRatioCombPHOSPi02760GeV , markerStylePHOS,markerSizeSpectrum, colorPHOSMass , colorPHOSMass, widthLinesBoxes, kTRUE);
//	graphSysErrRatioCombPHOSPi02760GeV->Draw("same,p");
	graphSysErrRatioCombPHOSPi02760GeV->Draw("E2,same");
	DrawGammaSetMarkerTGraphAsym(graphSysErrRatioCombConvPi02760GeV , markerStyleConv,markerSizeSpectrum, colorConv , colorConv, widthLinesBoxes, kTRUE);
//	graphSysErrRatioCombConvPi02760GeV->Draw("p,same");
	graphSysErrRatioCombConvPi02760GeV->Draw("E2,same");
	DrawGammaSetMarker(histoRatioCombPHOSPi02760GeV, markerStylePHOS,markerSizeSpectrum*1.5, colorPHOSMass , colorPHOSMass);
	histoRatioCombPHOSPi02760GeV->Draw("same,e1,x0");
	DrawGammaSetMarker(histoRatioCombConvPi02760GeV , markerStyleConv,markerSizeSpectrum*1.5, colorConv , colorConv);
	histoRatioCombConvPi02760GeV->DrawCopy("same,e1,x0");
	
	TLatex *labelRatioPi02760GeV = new TLatex(0.16,0.88,"pp #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( labelRatioPi02760GeV, 0.064,4);
	labelRatioPi02760GeV->Draw();
	
	DrawGammaLines(0., 20.,1., 1.,0.1,kBlack,2);
	
	TPad* padXSectionLegendOnlyPi0Ratio2760GeV = new TPad("padXSectionLegendOnlyPi0Ratio2760GeV", "",  0.65, 0.17, 0.96, 0.35,-1, -1, -2);
	DrawGammaPadSettings( padXSectionLegendOnlyPi0Ratio2760GeV, 0.15, 0.02, 0.03, 0.);
	padXSectionLegendOnlyPi0Ratio2760GeV->Draw();
	
	padXSectionLegendOnlyPi0Ratio2760GeV->cd();
	Double_t rowsLegendOnlyPi0Ratio2760GeV[4] 		= {0.78,0.45,0.05,0.00};
	Double_t columnsLegendOnlyPi0Ratio2760GeV[3] 	= {0.,0.4,0.7};
	Size_t textSizeLeftColumn2760GeV = 0.3;
	

	//****************************** Definition of the Legend ******************************************
	//****************** first Column **************************************************
	TLatex *textCTSOnlyRatioPi02760GeV = new TLatex(columnsLegendOnlyPi0Ratio2760GeV[0],rowsLegendOnlyPi0Ratio2760GeV[1],"PCM");
	SetStyleTLatex( textCTSOnlyRatioPi02760GeV, textSizeLeftColumn2760GeV,4);
	textCTSOnlyRatioPi02760GeV->Draw();
	TLatex *textPHOSOnlyRatioPi02760GeV = new TLatex(columnsLegendOnlyPi0Ratio2760GeV[0],rowsLegendOnlyPi0Ratio2760GeV[2],"PHOS");
	SetStyleTLatex( textPHOSOnlyRatioPi02760GeV, textSizeLeftColumn2760GeV,4);
	textPHOSOnlyRatioPi02760GeV->Draw();
	
	//****************** second Column *************************************************
	TLatex *textPi02760GeVstatOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio2760GeV[1],rowsLegendOnlyPi0Ratio2760GeV[0] ,"stat");
	SetStyleTLatex( textPi02760GeVstatOnlyRatioPi0, textSizeLeftColumn2760GeV,4);
	textPi02760GeVstatOnlyRatioPi0->Draw();
	TLatex *textPi02760GeVsysOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio2760GeV[2] ,rowsLegendOnlyPi0Ratio2760GeV[0],"syst");
	SetStyleTLatex( textPi02760GeVsysOnlyRatioPi0, textSizeLeftColumn2760GeV,4);
	textPi02760GeVsysOnlyRatioPi0->Draw();
	TMarker* markerCTSPi02760GeVOnlyRatioPi0 = CreateMarkerFromGraph(graphRatioCombConvPi0,columnsLegendOnlyPi0Ratio2760GeV[1]+ 0.1 ,rowsLegendOnlyPi0Ratio2760GeV[1]+ 0.1 ,scaleWidthLine);
	markerCTSPi02760GeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio2760GeV[1]+ 0.1 ,rowsLegendOnlyPi0Ratio2760GeV[1] + 0.1);
	TMarker* markerPHOSPi02760GeVOnlyRatioPi0 = CreateMarkerFromGraph(graphRatioCombPHOSPi0, columnsLegendOnlyPi0Ratio2760GeV[1]+ +0.1 ,rowsLegendOnlyPi0Ratio2760GeV[2]+ 0.1 ,scaleWidthLine);
	markerPHOSPi02760GeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio2760GeV[1]+ 0.1 ,rowsLegendOnlyPi0Ratio2760GeV[2]+ 0.1 );
	TBox* boxCTSPi02760GeVOnlyRatioPi0 = CreateBoxFromGraph(graphRatioCombConvPi0,columnsLegendOnlyPi0Ratio2760GeV[2]+0.1-offsetFit, rowsLegendOnlyPi0Ratio2760GeV[1]+ 0.1- offsetBoxSizeY, columnsLegendOnlyPi0Ratio2760GeV[2]+ 0.1+ offsetFit, rowsLegendOnlyPi0Ratio2760GeV[1]+ 0.1+ offsetBoxSizeY);
	boxCTSPi02760GeVOnlyRatioPi0->Draw("l");
	TBox* boxPHOSPi02760GeVOnlyRatioPi0 = CreateBoxFromGraph(graphRatioCombPHOSPi0,columnsLegendOnlyPi0Ratio2760GeV[2]+0.1-offsetFit, rowsLegendOnlyPi0Ratio2760GeV[2]+ 0.1- offsetBoxSizeY, columnsLegendOnlyPi0Ratio2760GeV[2]+ 0.1+ offsetFit, rowsLegendOnlyPi0Ratio2760GeV[2]+ 0.1+ offsetBoxSizeY);
	boxPHOSPi02760GeVOnlyRatioPi0->Draw("l");
	
	
	canvasInvXSectionOnlyRatioPi02760GeV->Update();
	canvasInvXSectionOnlyRatioPi02760GeV->Print(Form("%s/InvXSection_OnlyRatioPi02760GeV_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
	
	//**************************************************************************************************************************************
	//************************************ Inv Crosssection only ratio Pi0 900***************************************************************
	//**************************************************************************************************************************************
	
	TCanvas* canvasInvXSectionOnlyRatioPi0900GeV = new TCanvas("canvasInvXSectionOnlyRatioPi0900GeV","",200,10,1200,700);  // gives the page size
	DrawGammaCanvasSettings(canvasInvXSectionOnlyRatioPi0900GeV,  0.08, 0.01, 0.015, 0.13);
	
	canvasInvXSectionOnlyRatioPi0900GeV->SetLogx();
	TH2F * ratio2DInvXSectionOnlyPi0900GeV;
	ratio2DInvXSectionOnlyPi0900GeV = new TH2F("ratio2DInvXSectionOnlyPi0900GeV","ratio2DInvXSectionOnlyPi0900GeV",1000,0.23,30.,1000,0.5,1.85);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionOnlyPi0900GeV, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.064, 0.05,0.064, 0.8,0.6, 512, 505);
	ratio2DInvXSectionOnlyPi0900GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0900GeVSys , markerStylePHOS,markerSizeSpectrum*2, colorPHOSMass , colorPHOSMass, widthLinesBoxes, kTRUE);
//	graphRatioCombPHOSPi0900GeVSys->Draw("same,p");
	graphRatioCombPHOSPi0900GeVSys->Draw("E2,same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0900GeVSys , markerStyleConv,markerSizeSpectrum*2, colorConv , colorConv, widthLinesBoxes, kTRUE);
//	graphRatioCombConvPi0900GeVSys->Draw("p,same");
	graphRatioCombConvPi0900GeVSys->Draw("E2,same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0900GeV, markerStylePHOS,markerSizeSpectrum*2, colorPHOSMass , colorPHOSMass);
	graphRatioCombPHOSPi0900GeV->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0900GeV , markerStyleConv,markerSizeSpectrum*2, colorConv , colorConv);
	graphRatioCombConvPi0900GeV->Draw("p,same,e1");
	
	TLatex *labelRatioPi0900GeV = new TLatex(0.13,0.88,"#pi^{0}, #sqrt{#it{s}} = 0.9 TeV");
	SetStyleTLatex( labelRatioPi0900GeV, 0.07,4);
	labelRatioPi0900GeV->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	TPad* padXSectionLegendOnlyPi0Ratio900GeV = new TPad("padXSectionLegendOnlyPi0Ratio900GeV", "",  0.65, 0.8, 0.9, 0.95,-1, -1, -2);
	DrawGammaPadSettings( padXSectionLegendOnlyPi0Ratio900GeV, 0.15, 0.02, 0.03, 0.);
	padXSectionLegendOnlyPi0Ratio900GeV->Draw();
	
	padXSectionLegendOnlyPi0Ratio900GeV->cd();
	//****************************** Definition of the Legend ******************************************
	//****************** first Column **************************************************
	TLatex *textCTSOnlyRatioPi0900GeV = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],"PCM");
	SetStyleTLatex( textCTSOnlyRatioPi0900GeV, textSizeLeftColumn*2.,4);
	textCTSOnlyRatioPi0900GeV->Draw();
	TLatex *textPHOSOnlyRatioPi0900GeV = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[2],"PHOS");
	SetStyleTLatex( textPHOSOnlyRatioPi0900GeV, textSizeLeftColumn*2.,4);
	textPHOSOnlyRatioPi0900GeV->Draw();
	
	//****************** second Column *************************************************
	TLatex *textPi0900GeVstatOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
	SetStyleTLatex( textPi0900GeVstatOnlyRatioPi0, textSizeLeftColumn*2.,4);
	textPi0900GeVstatOnlyRatioPi0->Draw();
	TLatex *textPi0900GeVsysOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[1]+ offsetSystColumn*1.5 ,rowsLegendOnlyPi0Ratio[0],"syst");
	SetStyleTLatex( textPi0900GeVsysOnlyRatioPi0, textSizeLeftColumn*2.,4);
	textPi0900GeVsysOnlyRatioPi0->Draw();
	TMarker* markerCTSPi0900GeVOnlyRatioPi0 = CreateMarkerFromGraph(graphRatioCombConvPi0,columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[1]+ offsetMarkerY ,scaleWidthLine);
	markerCTSPi0900GeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[1]+ offsetMarkerY);
	TMarker* markerPHOSPi0900GeVOnlyRatioPi0 = CreateMarkerFromGraph(graphRatioCombPHOSPi0, columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[2]+ offsetMarkerY ,scaleWidthLine);
	markerPHOSPi0900GeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[2]+ offsetMarkerY );
	TBox* boxCTSPi0900GeVOnlyRatioPi0 = CreateBoxFromGraph(graphRatioCombConvPi0,columnsLegendOnlyPi0Ratio[1]+offsetSystColumn*1.5+offsetMarkerX-offsetFit, rowsLegendOnlyPi0Ratio[1]+ offsetMarkerY- offsetBoxSizeY, columnsLegendOnlyPi0Ratio[1]+ offsetSystColumn*1.5+ offsetMarkerX+ offsetFit, rowsLegendOnlyPi0Ratio[1]+ offsetMarkerY+ offsetBoxSizeY);
	boxCTSPi0900GeVOnlyRatioPi0->Draw("l");
	TBox* boxPHOSPi0900GeVOnlyRatioPi0 = CreateBoxFromGraph(graphRatioCombPHOSPi0,columnsLegendOnlyPi0Ratio[1]+offsetSystColumn*1.5+offsetMarkerX-offsetFit, rowsLegendOnlyPi0Ratio[2]+ offsetMarkerY- offsetBoxSizeY, columnsLegendOnlyPi0Ratio[1]+ offsetSystColumn*1.5+ offsetMarkerX+ offsetFit, rowsLegendOnlyPi0Ratio[2]+ offsetMarkerY+ offsetBoxSizeY);
	boxPHOSPi0900GeVOnlyRatioPi0->Draw("l");
	
	
	canvasInvXSectionOnlyRatioPi0900GeV->Update();
	canvasInvXSectionOnlyRatioPi0900GeV->Print(Form("%s/InvXSection_OnlyRatioPi0900GeV_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
	
	
	TCanvas * canvas4PartRatioIndMeas = new TCanvas("canvas4PartRatioIndMeas","",10,10,1400,1000);  // gives the page size		
	canvas4PartRatioIndMeas->cd();
// 	DrawGammaCanvasSettings( canvas4PartRatioIndMeas, 0.13, 0.0, 0.02, 0.09);
	
	TPad* pad4PartRatioIndMeas1 = new TPad("pad4PartRatioIndMeas1", "", 0., 0.52, 0.52, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad4PartRatioIndMeas1, 0.12, 0.0, 0.02, 0.);
	pad4PartRatioIndMeas1->Draw();
	TPad* pad4PartRatioIndMeas2 = new TPad("pad4PartRatioIndMeas2", "", 0., 0., 0.52, 0.52,-1, -1, -2);
	DrawGammaPadSettings( pad4PartRatioIndMeas2, 0.12, 0.0, 0., 0.12);
	pad4PartRatioIndMeas2->Draw();
	
	TPad* pad4PartRatioIndMeas3 = new TPad("pad4PartRatioIndMeas3", "", 0.52, 0.52, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( pad4PartRatioIndMeas3, 0.0, 0.02, 0.02, 0.);
	pad4PartRatioIndMeas3->Draw();
	TPad* pad4PartRatioIndMeas4 = new TPad("pad4PartRatioIndMeas4", "", 0.52, 0., 1., 0.52,-1, -1, -2);
	DrawGammaPadSettings( pad4PartRatioIndMeas4, 0.0, 0.02, 0., 0.12);
	pad4PartRatioIndMeas4->Draw();
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionOnlyPi0, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.062, 0.05,0.062, 0.8,0.82, 512, 505);
	ratio2DInvXSectionOnlyPi0->GetXaxis()->SetLabelOffset(-0.015);
	ratio2DInvXSectionOnlyPi0->GetYaxis()->SetLabelOffset(0.01);

	pad4PartRatioIndMeas1->cd();
	pad4PartRatioIndMeas1->SetLogx();
	ratio2DInvXSectionOnlyPi0->DrawCopy(); 
		
	graphRatioCombPHOSPi0900GeVSys->Draw("E2,same");
	graphRatioCombConvPi0900GeVSys->Draw("E2,same");
	graphRatioCombPHOSPi0900GeV->Draw("p,same,e1");
	graphRatioCombConvPi0900GeV->Draw("p,same,e1");
	
	TLatex *labelLegendA = new TLatex(0.16,0.88,"a)");
	SetStyleTLatex( labelLegendA, 0.062,4);
	labelLegendA->Draw();
	
	TLatex *labelRatioIndMeas900GeV = new TLatex(0.15,0.05,"#pi^{0}, pp #sqrt{#it{s}} = 0.9 TeV");
	SetStyleTLatex( labelRatioIndMeas900GeV, 0.062,4);
	labelRatioIndMeas900GeV->Draw();

	DrawGammaLines(0., 30.,1., 1.,1.,kBlack,2);
	TPad* padXSectionLegengOnly4Parted = new TPad("padXSectionLegengOnly4Parted", "",  0.5, 0.73, 0.96, 0.95,-1, -1, -2);
	DrawGammaPadSettings( padXSectionLegengOnly4Parted, 0.0, 0.00, 0.0, 0.);
	padXSectionLegengOnly4Parted->Draw();
	
	padXSectionLegengOnly4Parted->cd();
	//****************************** Definition of the Legend ******************************************
	//****************** first Column **************************************************
	textCTSOnlyRatioPi0900GeV->Draw();
	textPHOSOnlyRatioPi0900GeV->Draw();
	
	//****************** second Column *************************************************
	textPi0900GeVstatOnlyRatioPi0->Draw();
	textPi0900GeVsysOnlyRatioPi0->Draw();
	markerCTSPi0900GeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[1]+ offsetMarkerY);
	markerPHOSPi0900GeVOnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0Ratio[1]+ offsetMarkerX ,rowsLegendOnlyPi0Ratio[2]+ offsetMarkerY );
	boxCTSPi0900GeVOnlyRatioPi0->Draw("l");
	boxPHOSPi0900GeVOnlyRatioPi0->Draw("l");
	
	pad4PartRatioIndMeas1->Update();
	pad4PartRatioIndMeas2->cd();
	pad4PartRatioIndMeas2->SetLogx();
 	ratio2DInvXSectionOnlyPi0->DrawCopy(); 
	graphRatioCombPHOSPi0Sys->Draw("E2same");
	graphRatioCombConvPi0Sys->Draw("E2same");
	graphRatioCombPHOSPi0->Draw("p,same,e1");
	graphRatioCombConvPi0->Draw("p,same,e1");
	DrawGammaLines(0., 30.,1., 1.,1.,kBlack,2);
	TLatex *labelRatioIndMeasPbPb7TeV = new TLatex(0.15,0.17,"#pi^{0}, pp #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( labelRatioIndMeasPbPb7TeV, 0.06,4);
	labelRatioIndMeasPbPb7TeV->Draw();
	TLatex *labelLegendC = new TLatex(0.16,0.92,"c)");
	SetStyleTLatex( labelLegendC, 0.062,4);
	labelLegendC->Draw();

	pad4PartRatioIndMeas2->Update();

	pad4PartRatioIndMeas3->cd();
	pad4PartRatioIndMeas3->SetLogx();
	
	ratio2DInvXSectionOnlyPi0->DrawCopy(); 
	graphSysErrRatioCombPHOSPi02760GeV->Draw("E2,same");
	graphSysErrRatioCombConvPi02760GeV->Draw("E2,same");
	histoRatioCombPHOSPi02760GeV->SetMarkerSize(2.);
	histoRatioCombConvPi02760GeV->SetMarkerSize(2.);
	histoRatioCombPHOSPi02760GeV->Draw("same,e1");
	histoRatioCombConvPi02760GeV->DrawCopy("same,e1");
	
	DrawGammaLines(0., 30.,1., 1.,1.,kBlack,2);
	
	TLatex *labelRatioIndMeas2760GeV = new TLatex(0.03,0.05,"#pi^{0}, pp #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( labelRatioIndMeas2760GeV, 0.062,4);
	labelRatioIndMeas2760GeV->Draw();
	TLatex *labelLegendB = new TLatex(0.04,0.88,"b)");
	SetStyleTLatex( labelLegendB, 0.062,4);
	labelLegendB->Draw();
	
	pad4PartRatioIndMeas3->Update();
	pad4PartRatioIndMeas4->cd();
	pad4PartRatioIndMeas4->SetLogx();
	ratio2DInvXSectionOnlyPi0->DrawCopy(); 
	graphRatioCombPHOSEta7TeVSys->Draw("E2same");
	graphRatioCombConvEta7TeVSys->Draw("E2same");
	graphRatioCombPHOSEta7TeV->Draw("p,same,e1");
	graphRatioCombConvEta7TeV->Draw("p,same,e1");
	DrawGammaLines(0., 30.,1., 1.,1.,kBlack,2);
	
	TLatex *labelLegendD = new TLatex(0.04,0.92,"d)");
	SetStyleTLatex( labelLegendD, 0.062,4);
	labelLegendD->Draw();
	
	TLatex *labelRatioIndMeasEta7TeV = new TLatex(0.03,0.17,"#eta, pp #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( labelRatioIndMeasEta7TeV, 0.06,4);
	labelRatioIndMeasEta7TeV->Draw();
	
	pad4PartRatioIndMeas4->Update();

	canvas4PartRatioIndMeas->Update();	
	canvas4PartRatioIndMeas->SaveAs(Form("%s/InvYieldOnlyRatio_All4Parted_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
	delete pad4PartRatioIndMeas1;	
	delete pad4PartRatioIndMeas2;	
	delete pad4PartRatioIndMeas3;	
	delete pad4PartRatioIndMeas4;	
	delete canvas4PartRatioIndMeas;	

	
		
	//****************************** Definition of the Legend ******************************************

		//*************** Size factors ********************
		Double_t scaleMarkerNLO = 1.5;
		Double_t scaleLineWidthNLO = 1.;

	//**************************************************************************************************************
	//***************************** Larger Variant of Paper Plot with NLO ******************************************
	//**************************************************************************************************************
	
	TCanvas* canvasInvXSectionNLOVar2 = new TCanvas("canvasInvXSectionNLOVar2","",200,10,1300,3000);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOVar2,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionNLOVar2 = new TPad("padComparisonXSectionNLOVar2", "", 0., 0.49, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionNLOVar2, 0.15, 0.02, 0.03, 0.);
	padComparisonXSectionNLOVar2->Draw();
	
	TPad* padXSectionNLOVar2RatioPi07TeV = new TPad("padXSectionNLOVar2RatioPi07TeV", "", 0., 0.34, 1., 0.49,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOVar2RatioPi07TeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOVar2RatioPi07TeV->Draw();
	
	TPad* padXSectionNLOVar2RatioPi0900GeV = new TPad("padXSectionNLOVar2RatioPi0900GeV", "", 0., 0.19, 1., 0.34,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOVar2RatioPi0900GeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOVar2RatioPi0900GeV->Draw();
	
	TPad* padXSectionNLOVar2RatioEta7TeV = new TPad("padXSectionNLOVar2RatioEta7TeV", "", 0., 0., 1., 0.19,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOVar2RatioEta7TeV, 0.15, 0.02, 0., 0.2);
	padXSectionNLOVar2RatioEta7TeV->Draw();
	
	
	padComparisonXSectionNLOVar2->cd();
	padComparisonXSectionNLOVar2->SetLogy();		
	padComparisonXSectionNLOVar2->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionNLOVar2 = new TH2F("histo2DInvXSectionNLOVar2","histo2DInvXSectionNLOVar2",1000,0.23,30.,1000,2e-4,10e12);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLOVar2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionNLOVar2->GetYaxis()->SetRangeUser(2e-4,10e11);
	histo2DInvXSectionNLOVar2->DrawCopy(); 

	graphInvCrossSectionPi0Comb900GeV = ScaleGraph(graphInvCrossSectionPi0Comb900GeV,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb900GeV, markerStyleCommmonSpectrumPi0900GeV,markerSizeCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionPi0Comb900GeV->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionPi0Comb900GeV->Draw("p,E2same");	
	
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb7TeV, markerStyleCommmonSpectrumPi07TeV,markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionPi0Comb7TeV->SetLineWidth(widthCommonErrors);	
	graphInvCrossSectionPi0Comb7TeV->Draw("p,E2same");
	
	graphInvCrossSectionEtaComb7TeV = ScaleGraph(graphInvCrossSectionEtaComb7TeV,1e-3);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb7TeV,  markerStyleCommmonSpectrumEta7TeV,markerSizeCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionEtaComb7TeV->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb7TeV->Draw("p,E2same");

	DrawGammaSetMarkerTF1( fitInvCrossSectionPi0, styleFitCommonSpectrum, 0.5*widthCommonFit, colorCommonSpectrumPi07TeV);
	fitInvCrossSectionPi0->Draw("same");
	
	SetStyleHisto(histoFitInvCrossSectionPi0900GeV, 0.5*widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi0900GeV);
	histoFitInvCrossSectionPi0900GeV->Draw("same,c");
	
	SetStyleHisto(histoFitInvCrossSectionEta7TeV, 0.5*widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumEta7TeV);
	histoFitInvCrossSectionEta7TeV->Draw("same,c");

	DrawGammaNLOTGraph( graphNLOMuHalfPi07TeV, 0.5*widthCommonFit, styleLineNLOMuHalf, colorNLOPi07TeVMuHalf);
	graphNLOMuHalfPi07TeV->Draw("same,c");
	DrawGammaNLOTGraph( graphNLOMuOnePi07TeV, 0.5*widthCommonFit, styleLineNLOMuOne, colorNLOPi07TeVMuOne);
	graphNLOMuOnePi07TeV->Draw("same,c");
	DrawGammaNLOTGraph( graphNLOMuTwoPi07TeV,0.5*widthCommonFit, styleLineNLOMuTwo, colorNLOPi07TeVMuTwo);
	graphNLOMuTwoPi07TeV->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOBKKPi07TeVMuTwo, 0.5*widthCommonFit, styleLineNLOMuTwoBKK, colorNLOPi07TeVMuTwo);
	DrawGammaNLOTGraph( graphRatioCombNLODSSPi07TeVMuTwo, 0.5*widthCommonFit, styleLineNLOMuTwoDSS, colorNLOPi07TeVMuTwo);
	
	graphNLOMuHalfPi0900GeV= ScaleGraph(graphNLOMuHalfPi0900GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuHalfPi0900GeV, 0.5*widthCommonFit, styleLineNLOMuHalf, colorNLOPi0900GeVMuHalf);
	graphNLOMuHalfPi0900GeV->Draw("same,c");
	graphNLOMuOnePi0900GeV= ScaleGraph(graphNLOMuOnePi0900GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuOnePi0900GeV, 0.5*widthCommonFit, styleLineNLOMuOne, colorNLOPi0900GeVMuOne);
	graphNLOMuOnePi0900GeV->Draw("same,c");
	graphNLOMuTwoPi0900GeV = ScaleGraph(graphNLOMuTwoPi0900GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuTwoPi0900GeV, 0.5*widthCommonFit, styleLineNLOMuTwo, colorNLOPi0900GeVMuTwo);
	graphNLOMuTwoPi0900GeV->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOBKKPi0900GeVMuTwo, 0.5*widthCommonFit, styleLineNLOMuTwoBKK, colorNLOPi0900GeVMuTwo);
	
	graphNLOMuHalfEta7TeV= ScaleGraph(graphNLOMuHalfEta7TeV,1e-3);
	DrawGammaNLOTGraph( graphNLOMuHalfEta7TeV, 0.5*widthCommonFit, styleLineNLOMuHalf, colorNLOEta7TeVMuHalf);
	graphNLOMuHalfEta7TeV->Draw("same,c");
	graphNLOMuOneEta7TeV= ScaleGraph(graphNLOMuOneEta7TeV,1e-3);
	DrawGammaNLOTGraph( graphNLOMuOneEta7TeV, 0.5*widthCommonFit, styleLineNLOMuOne, colorNLOEta7TeVMuOne);
	graphNLOMuOneEta7TeV->Draw("same,c");
	graphNLOMuTwoEta7TeV = ScaleGraph(graphNLOMuTwoEta7TeV,1e-3);
	DrawGammaNLOTGraph( graphNLOMuTwoEta7TeV,0.5*widthCommonFit, styleLineNLOMuTwo, colorNLOEta7TeVMuTwo);
	graphNLOMuTwoEta7TeV->Draw("same,c");

	TLatex *labelLegendAVar2 = new TLatex(15,2E11,"a)");
	SetStyleTLatex( labelLegendAVar2, 0.045,4,1,kFALSE);
	labelLegendAVar2->Draw();
	
	labelScalingPi07TeV->Draw();
	labelScalingPi0900GeV->Draw();
	labelScalingEta7TeV->Draw();
	DrawNormalizationErrorText(normalizationInvX3En[0],normalizationInvX3En[1],normalizationInvX3En[2],normalizationInvX3En[3],normalizationInvX3En[4],"No2.76"); 
	
	//****************************** Definition of the Legend ******************************************
		//**************** Row def ************************
		Double_t rowsNLOLegendVar2[8] = {0.89,0.79,0.68,0.54,0.43,0.31,0.17,0.05};
		Double_t columnsNLOLegendVar2[4] = {0.,0.27,0.51,0.78};
		//*************** Label sizes *********************
		Double_t textSizeLeftLabelsVar2 = 0.11;
		Double_t textSizeTopLablesVar2 = 0.11;
		Double_t textSizeTopLowerLablesVar2 = 0.105;
		//*************** Size factors ********************
		Double_t scaleMarkerNLOVar2 = 1.;
		Double_t scaleLineWidthNLOVar2 = 1.2;
		//*************** Offsets *************************
		Double_t offsetNLOLegendMarkerXVar2 = 0.09;
		Double_t offsetNLOLegendMarkerYVar2 = 0.03;
		Double_t offsetNLOLegendBoxVar2 = 0.05;
		Double_t offsetNLOLegendLineVar2 = 0.06;

		TPad* padXSectionNLOVar2Legend = new TPad("padXSectionNLOVar2Legend", "", 0.17, 0.02, 0.95, 0.205,-1, -1, -2); 
		//TPad* padXSectionNLOVar2Legend = new TPad("padXSectionNLOVar2Legend", "", 0.17, 0.02, 0.95, 0.23,-1, -1, -2);
		DrawGammaPadSettings( padXSectionNLOVar2Legend, 0., 0., 0., 0.);
		padXSectionNLOVar2Legend->Draw();
		padXSectionNLOVar2Legend->cd();
		
		//*************** first Column **********************************************************
		TLatex *textSpectrumVar = new TLatex(columnsNLOLegendVar2[0],rowsNLOLegendVar2[2],"combined Spec.");
		SetStyleTLatex( textSpectrumVar, textSizeLeftLabelsVar2,4);
		textSpectrumVar->Draw();
		TLatex *textFitCombVar2 = new TLatex(columnsNLOLegendVar2[0],rowsNLOLegendVar2[3],"fit combined");
		SetStyleTLatex( textFitCombVar2, textSizeLeftLabelsVar2,4);
		textFitCombVar2->Draw();
		TLatex *textNLOMuHalfVar2 = new TLatex(columnsNLOLegendVar2[0],rowsNLOLegendVar2[4],"NLO #mu= 0.5 #it{p}_{T}");
		SetStyleTLatex( textNLOMuHalfVar2, textSizeLeftLabelsVar2,4);
		textNLOMuHalfVar2->Draw();
		TLatex *textNLOMuOneVar2 = new TLatex(columnsNLOLegendVar2[0],rowsNLOLegendVar2[5],"NLO #mu= #it{p}_{T}");
		SetStyleTLatex( textNLOMuOneVar2, textSizeLeftLabelsVar2,4);
		textNLOMuOneVar2->Draw();
		TLatex *textNLOMuTwoVar2 = new TLatex(columnsNLOLegendVar2[0],rowsNLOLegendVar2[6],"NLO #mu= 2 #it{p}_{T}");
		SetStyleTLatex( textNLOMuTwoVar2, textSizeLeftLabelsVar2,4);
		textNLOMuTwoVar2->Draw();
		TLatex *textNLOMuTwoBKKVar2 = new TLatex(columnsNLOLegendVar2[0],rowsNLOLegendVar2[7],"NLO #mu= 2 #it{p}_{T} (BKK)");
		SetStyleTLatex( textNLOMuTwoBKKVar2, textSizeLeftLabelsVar2,4);
		textNLOMuTwoBKKVar2->Draw();
		
		
		//*************** second Column **********************************************************
		TLatex *textPi07TeVNLOVar2 = new TLatex(columnsNLOLegendVar2[1],rowsNLOLegendVar2[0],"#pi^{0}, #sqrt{#it{s}} = 7 TeV");
		SetStyleTLatex( textPi07TeVNLOVar2, textSizeTopLablesVar2,4);
		textPi07TeVNLOVar2->Draw();
		TLatex *textPi07TeVNLOsysVar2 = new TLatex(columnsNLOLegendVar2[1],rowsNLOLegendVar2[1],"syst. + stat.");
		SetStyleTLatex( textPi07TeVNLOsysVar2, textSizeTopLowerLablesVar2,4);
		textPi07TeVNLOsysVar2->Draw();
		TBox* boxCombinedPi07TeVVar2 = CreateBoxFromGraph(graphInvCrossSectionPi0Comb7TeV,columnsNLOLegendVar2[1]+offsetNLOLegendMarkerXVar2-offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2- offsetNLOLegendBoxVar2, columnsNLOLegendVar2[1]+offsetNLOLegendMarkerXVar2+offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2+offsetNLOLegendBoxVar2);
		boxCombinedPi07TeVVar2->Draw("l");
		TMarker* markerCombinedPi07TeVVar2 = CreateMarkerFromGraph(graphInvCrossSectionPi0Comb7TeV,columnsNLOLegendVar2[1]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2 ,scaleMarkerNLOVar2);
		markerCombinedPi07TeVVar2->DrawMarker(columnsNLOLegendVar2[1]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		TLine * lineFit7TeVNLOVar2 = CreateLineFromFit(fitInvCrossSectionPi0, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[3]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[3]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2);
		lineFit7TeVNLOVar2->Draw("same");
		TLine * lineNLOPi07TeVMuHalfVar2 = CreateLineFromGraph(graphNLOMuHalfPi07TeV, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[4] + offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2, rowsNLOLegendVar2[4]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOPi07TeVMuHalfVar2->Draw("same");
		TLine * lineNLOPi07TeVMuOneVar2 = CreateLineFromGraph(graphNLOMuOnePi07TeV,columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[5]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[5]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOPi07TeVMuOneVar2->Draw("same");
		TLine * lineNLOPi07TeVMuTwoVar2 = CreateLineFromGraph(graphNLOMuTwoPi07TeV, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[6]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[6]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOPi07TeVMuTwoVar2->Draw("same");
		TLine * lineNLOPi07TeVMuTwoBKKVar2 = CreateLineFromGraph(graphRatioCombNLOBKKPi07TeVMuTwo, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[7]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[7]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOPi07TeVMuTwoBKKVar2->Draw("same");
// 		TLine * lineNLOPi07TeVMuTwoDSSVar2 = CreateLineFromGraph(graphRatioCombNLODSSPi07TeVMuTwo, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[7]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[1]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[7]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
// 		lineNLOPi07TeVMuTwoDSSVar2->Draw("same");		
		//*************** third Column **********************************************************
		TLatex *textPi0900GeVNLOVar2 = new TLatex(columnsNLOLegendVar2[2],rowsNLOLegendVar2[0],"#pi^{0}, #sqrt{#it{s}} = 0.9 TeV");
		SetStyleTLatex( textPi0900GeVNLOVar2, textSizeTopLablesVar2,4);
		textPi0900GeVNLOVar2->Draw();
		TLatex *textPi0900GeVNLOsysVar2 = new TLatex(columnsNLOLegendVar2[2],rowsNLOLegendVar2[1],"syst. + stat.");
		SetStyleTLatex( textPi0900GeVNLOsysVar2, textSizeTopLowerLablesVar2,4);
		textPi0900GeVNLOsysVar2->Draw();
		TBox* boxCombinedPi0900GeVVar2 = CreateBoxFromGraph(graphInvCrossSectionPi0Comb900GeV,columnsNLOLegendVar2[2]+offsetNLOLegendMarkerXVar2-offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2- offsetNLOLegendBoxVar2, columnsNLOLegendVar2[2]+offsetNLOLegendMarkerXVar2+offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2+offsetNLOLegendBoxVar2);
		boxCombinedPi0900GeVVar2->Draw("l");
		TMarker* markerCombinedPi0900GeVVar2 = CreateMarkerFromGraph(graphInvCrossSectionPi0Comb900GeV,columnsNLOLegendVar2[2]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2 ,scaleMarkerNLOVar2);
		markerCombinedPi0900GeVVar2->DrawMarker(columnsNLOLegendVar2[2]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		TLine * lineFitPi0900GeVNLOVar2 =CreateLineFromHisto(histoFitInvCrossSectionPi0900GeV, columnsNLOLegendVar2[2]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[3]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[2]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[3]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2);
		lineFitPi0900GeVNLOVar2->Draw("same");
		TLine * lineNLOPi0900GeVMuHalfVar2 = CreateLineFromGraph(graphNLOMuHalfPi0900GeV,  columnsNLOLegendVar2[2]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[4] + offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[2]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2, rowsNLOLegendVar2[4]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOPi0900GeVMuHalfVar2->Draw("same");
		TLine * lineNLOPi0900GeVMuOneVar2 = CreateLineFromGraph(graphNLOMuOnePi0900GeV, columnsNLOLegendVar2[2]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[5]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[2]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[5]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOPi0900GeVMuOneVar2->Draw("same");
		TLine * lineNLOPi0900GeVMuTwoVar2 = CreateLineFromGraph(graphNLOMuTwoPi0900GeV, columnsNLOLegendVar2[2]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[6]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[2]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[6]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOPi0900GeVMuTwoVar2->Draw("same");
		TLine * lineNLOPi0900GeVMuTwoBKKVar2 = CreateLineFromGraph(graphRatioCombNLOBKKPi0900GeVMuTwo, columnsNLOLegendVar2[2]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[7]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[2]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[7]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOPi0900GeVMuTwoBKKVar2->Draw("same");
		
		
		//**************** forth Column **********************************************************
		TLatex *textEta7TeVNLOVar2 = new TLatex(columnsNLOLegendVar2[3],rowsNLOLegendVar2[0],"#eta, #sqrt{#it{s}} = 7 TeV");
		SetStyleTLatex( textEta7TeVNLOVar2, textSizeTopLablesVar2,4);
		textEta7TeVNLOVar2->Draw();
		TLatex *textEta7TeVNLOsysVar2 = new TLatex(columnsNLOLegendVar2[3],rowsNLOLegendVar2[1],"syst. + stat.");
		SetStyleTLatex( textEta7TeVNLOsysVar2, textSizeTopLowerLablesVar2,4);
		textEta7TeVNLOsysVar2->Draw();
		TBox* boxCombinedEta7TeVVar2 = CreateBoxFromGraph(graphInvCrossSectionEtaComb7TeV,columnsNLOLegendVar2[3]+offsetNLOLegendMarkerXVar2-offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2- offsetNLOLegendBoxVar2, columnsNLOLegendVar2[3]+offsetNLOLegendMarkerXVar2+offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2+offsetNLOLegendBoxVar2);
		boxCombinedEta7TeVVar2->Draw("l");
		TMarker* markerCombinedEta7TeVVar2 = CreateMarkerFromGraph(graphInvCrossSectionEtaComb7TeV,columnsNLOLegendVar2[3]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2 ,scaleMarkerNLOVar2);
		markerCombinedEta7TeVVar2->DrawMarker(columnsNLOLegendVar2[3]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		TLine * lineFitEta7TeVNLOVar2 = CreateLineFromHisto(histoFitInvCrossSectionEta7TeV,columnsNLOLegendVar2[3]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[3]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[3]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[3]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2);
		lineFitEta7TeVNLOVar2->Draw("same");
		TLine * lineNLOEta7TeVMuHalfVar2 =  CreateLineFromGraph(graphNLOMuHalfEta7TeV, columnsNLOLegendVar2[3]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[4] + offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[3]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2, rowsNLOLegendVar2[4]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOEta7TeVMuHalfVar2->Draw("same");
		TLine * lineNLOEta7TeVMuOneVar2 =  CreateLineFromGraph(graphNLOMuOneEta7TeV,columnsNLOLegendVar2[3]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[5]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[3]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[5]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOEta7TeVMuOneVar2->Draw("same");
		TLine * lineNLOEta7TeVMuTwoVar2 =  CreateLineFromGraph(graphNLOMuTwoEta7TeV, columnsNLOLegendVar2[3]+ offsetNLOLegendMarkerXVar2- offsetNLOLegendLineVar2, rowsNLOLegendVar2[6]+offsetNLOLegendMarkerYVar2, columnsNLOLegendVar2[3]+ offsetNLOLegendMarkerXVar2+ offsetNLOLegendLineVar2,rowsNLOLegendVar2[6]+offsetNLOLegendMarkerYVar2, scaleLineWidthNLOVar2); 
		lineNLOEta7TeVMuTwoVar2->Draw("same");
	//************************************************* End Legend ***************************************************
	
	padXSectionNLOVar2RatioPi07TeV->cd();
	padXSectionNLOVar2RatioPi07TeV->SetLogx();
// 	padXSectionNLOVar2RatioPi07TeV->SetLogy();
	
	TH2F * ratio2DInvXSectionNLOVar2Pi0 = new TH2F("ratio2DInvXSectionNLOVar2Pi0","ratio2DInvXSectionNLOVar2Pi0",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOVar2Pi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.11,0.115, 1,0.53, 512, 505);
	ratio2DInvXSectionNLOVar2Pi0->GetYaxis()->SetMoreLogLabels(kTRUE);
	ratio2DInvXSectionNLOVar2Pi0->GetYaxis()->SetNdivisions(515);
	ratio2DInvXSectionNLOVar2Pi0->GetYaxis()->SetNoExponent(kTRUE);
	ratio2DInvXSectionNLOVar2Pi0->GetXaxis()->SetTickLength(0.07);
	ratio2DInvXSectionNLOVar2Pi0->DrawCopy(); 

	DrawGammaNLOTGraph( graphRatioCombNLOPi07TeVMuHalf, 0.5*widthCommonFit, styleLineNLOMuHalf, colorNLOPi07TeVMuHalf);
	graphRatioCombNLOPi07TeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOPi07TeVMuOne, 0.5*widthCommonFit, styleLineNLOMuOne, colorNLOPi07TeVMuOne);
	graphRatioCombNLOPi07TeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOPi07TeVMuTwo, 0.5*widthCommonFit, styleLineNLOMuTwo, colorNLOPi07TeVMuTwo);
	graphRatioCombNLOPi07TeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi07TeVMuTwo->Draw("same,c");
// 	graphRatioCombNLODSSPi07TeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi07TeVRatio->Draw();
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit, markerStyleCommmonSpectrumPi07TeV,markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFit->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFit->Draw("p,E2same");
	
	TLatex *labelLegendBVar2 = new TLatex(0.18,0.75,"b)");
	SetStyleTLatex( labelLegendBVar2, 0.12,4);
	labelLegendBVar2->Draw();
		
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOVar2RatioPi0900GeV->cd();
	padXSectionNLOVar2RatioPi0900GeV->SetLogx();
// 	padXSectionNLOVar2RatioPi0900GeV->SetLogy();
	TH2F * ratio2DInvXSectionNLOVar2Pi0900GeV;
	ratio2DInvXSectionNLOVar2Pi0900GeV = new TH2F("ratio2DInvXSectionNLOVar2Pi0900GeV","ratio2DInvXSectionNLOVar2Pi0900GeV",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOVar2Pi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.11,0.115, 1,0.53, 512, 505);
	ratio2DInvXSectionNLOVar2Pi0900GeV->GetYaxis()->SetMoreLogLabels(kTRUE);
	ratio2DInvXSectionNLOVar2Pi0900GeV->GetYaxis()->SetNdivisions(515);
	ratio2DInvXSectionNLOVar2Pi0900GeV->GetYaxis()->SetNoExponent(kTRUE);

	ratio2DInvXSectionNLOVar2Pi0900GeV->GetXaxis()->SetTickLength(0.07);
	ratio2DInvXSectionNLOVar2Pi0900GeV->DrawCopy(); 

	DrawGammaNLOTGraph( graphRatioCombNLOPi0900GeVMuHalf, 0.5*widthCommonFit, styleLineNLOMuHalf, colorNLOPi0900GeVMuHalf);
	graphRatioCombNLOPi0900GeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOPi0900GeVMuOne, 0.5*widthCommonFit, styleLineNLOMuOne, colorNLOPi0900GeVMuOne);
	graphRatioCombNLOPi0900GeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOPi0900GeVMuTwo, 0.5*widthCommonFit, styleLineNLOMuTwo, colorNLOPi0900GeVMuTwo);
	graphRatioCombNLOPi0900GeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi0900GeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi0900GeVRatio->Draw();
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit900GeV, markerStyleCommmonSpectrumPi0900GeV,markerSizeCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFit900GeV->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFit900GeV->Draw("p,E2same");
	boxErrorSigmaPi0900GeVRatio->Draw();	
	
	TLatex *labelLegendCVar2 = new TLatex(0.18,0.75,"c)");
	SetStyleTLatex( labelLegendCVar2, 0.12,4);
	labelLegendCVar2->Draw();
		
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOVar2RatioEta7TeV->cd();
	padXSectionNLOVar2RatioEta7TeV->SetLogx();
// 	padXSectionNLOVar2RatioEta7TeV->SetLogy();
	TH2F * ratio2DInvXSectionNLOVar2Eta=  new TH2F("ratio2DInvXSectionNLOVar2Eta","ratio2DInvXSectionNLOVar2Eta",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOVar2Eta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.09,0.095, 0.09,0.095, 0.95,0.63, 510, 505);
	ratio2DInvXSectionNLOVar2Eta->GetYaxis()->SetMoreLogLabels(kTRUE);
 	ratio2DInvXSectionNLOVar2Eta->GetYaxis()->SetNdivisions(512);
	ratio2DInvXSectionNLOVar2Eta->GetYaxis()->SetNoExponent(kTRUE);
	ratio2DInvXSectionNLOVar2Eta->GetXaxis()->SetMoreLogLabels(kTRUE);
	ratio2DInvXSectionNLOVar2Eta->GetXaxis()->SetNoExponent(kTRUE);
// 	ratio2DInvXSectionNLOVar2Eta->GetXaxis()->SetLabelOffset(-0.03);
	ratio2DInvXSectionNLOVar2Eta->GetXaxis()->SetTickLength(0.07);
	ratio2DInvXSectionNLOVar2Eta->DrawCopy(); 

	boxErrorSigmaEta7TeVRatio->Draw();
	
	DrawGammaNLOTGraph( graphRatioCombNLOEta7TeVMuHalf, 0.5*widthCommonFit, styleLineNLOMuHalf, colorNLOEta7TeVMuHalf);
	graphRatioCombNLOEta7TeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta7TeVMuOne, 0.5*widthCommonFit, styleLineNLOMuOne, colorNLOEta7TeVMuOne);
	graphRatioCombNLOEta7TeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta7TeVMuTwo, 0.5*widthCommonFit, styleLineNLOMuTwo, colorNLOEta7TeVMuTwo);
	graphRatioCombNLOEta7TeVMuTwo->Draw("same,c");
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta7TeV,markerStyleCommmonSpectrumEta7TeV,markerSizeCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFitEta7TeV->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFitEta7TeV->Draw("p,E2same");

	TLatex *labelLegendDVar2 = new TLatex(0.18,0.85,"d)");
	SetStyleTLatex( labelLegendDVar2, 0.10,4);
	labelLegendDVar2->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	canvasInvXSectionNLOVar2->Update();
	canvasInvXSectionNLOVar2->Print(Form("%s/InvXSectionNLOVar2_Paper.%s",outputDir.Data(),suffix.Data()));


	//**************************************************************************************************************
	//***************************** Larger Variant of Paper Plot with NLO + separated sys & stat errror*************
	//**************************************************************************************************************
	
	TCanvas* canvasInvXSectionNLOVar2Sep = new TCanvas("canvasInvXSectionNLOVar2Sep","",200,10,1300,3000);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOVar2Sep,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionNLOVar2Sep = new TPad("padComparisonXSectionNLOVar2Sep", "", 0., 0.49, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionNLOVar2Sep, 0.15, 0.02, 0.03, 0.);
	padComparisonXSectionNLOVar2Sep->Draw();
	
	TPad* padXSectionNLOVar2SepRatioPi07TeV = new TPad("padXSectionNLOVar2SepRatioPi07TeV", "", 0., 0.34, 1., 0.49,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOVar2SepRatioPi07TeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOVar2SepRatioPi07TeV->Draw();
	
	TPad* padXSectionNLOVar2SepRatioPi0900GeV = new TPad("padXSectionNLOVar2SepRatioPi0900GeV", "", 0., 0.19, 1., 0.34,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOVar2SepRatioPi0900GeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOVar2SepRatioPi0900GeV->Draw();
	
	TPad* padXSectionNLOVar2SepRatioEta7TeV = new TPad("padXSectionNLOVar2SepRatioEta7TeV", "", 0., 0., 1., 0.19,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOVar2SepRatioEta7TeV, 0.15, 0.02, 0., 0.2);
	padXSectionNLOVar2SepRatioEta7TeV->Draw();
	
	
	padComparisonXSectionNLOVar2Sep->cd();
	padComparisonXSectionNLOVar2Sep->SetLogy();		
	padComparisonXSectionNLOVar2Sep->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionNLOVar2Sep = new TH2F("histo2DInvXSectionNLOVar2Sep","histo2DInvXSectionNLOVar2Sep",1000,0.23,30.,1000,2e-4,10e12);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLOVar2Sep, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionNLOVar2Sep->GetYaxis()->SetRangeUser(2e-4,10e11);
	histo2DInvXSectionNLOVar2Sep->DrawCopy(); 

	graphInvCrossSectionPi0Comb900GeVStatErr = ScaleGraph(graphInvCrossSectionPi0Comb900GeVStatErr,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb900GeVStatErr, markerStyleCommmonSpectrumPi0900GeV,markerSizeCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionPi0Comb900GeVStatErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionPi0Comb900GeVSysErr = ScaleGraph(graphInvCrossSectionPi0Comb900GeVSysErr,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb900GeVSysErr, markerStyleCommmonSpectrumPi0900GeV,markerSizeCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi0900GeVBox);
	graphInvCrossSectionPi0Comb900GeVSysErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionPi0Comb900GeVSysErr->Draw("2,same");	
	graphInvCrossSectionPi0Comb900GeVStatErr->Draw("p,same");	
	
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb7TeVStatErr, markerStyleCommmonSpectrumPi07TeV,markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionPi0Comb7TeVStatErr->SetLineWidth(widthCommonErrors);	
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb7TeVSysErr, markerStyleCommmonSpectrumPi07TeV,markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi07TeVBox);
	graphInvCrossSectionPi0Comb7TeVSysErr->SetLineWidth(widthCommonErrors);	
	graphInvCrossSectionPi0Comb7TeVSysErr->Draw("2,same");
	graphInvCrossSectionPi0Comb7TeVStatErr->Draw("p,same");
	
	graphInvCrossSectionEtaComb7TeVStatErr = ScaleGraph(graphInvCrossSectionEtaComb7TeVStatErr,1e-3);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb7TeVStatErr,  markerStyleCommmonSpectrumEta7TeV,markerSizeCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionEtaComb7TeVStatErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb7TeVSysErr = ScaleGraph(graphInvCrossSectionEtaComb7TeVSysErr,1e-3);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb7TeVSysErr,  markerStyleCommmonSpectrumEta7TeV,markerSizeCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumEta7TeVBox);
	graphInvCrossSectionEtaComb7TeVSysErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb7TeVSysErr->Draw("2,same");
	graphInvCrossSectionEtaComb7TeVStatErr->Draw("p,same");

	fitInvCrossSectionPi0->Draw("same");
	histoFitInvCrossSectionPi0900GeV->Draw("same,c");
	histoFitInvCrossSectionEta7TeV->Draw("same,c");
	graphNLOMuHalfPi07TeV->Draw("same,c");
	graphNLOMuOnePi07TeV->Draw("same,c");
	graphNLOMuTwoPi07TeV->Draw("same,c");
	
	graphNLOMuHalfPi0900GeV->Draw("same,c");
	graphNLOMuOnePi0900GeV->Draw("same,c");
	graphNLOMuTwoPi0900GeV->Draw("same,c");
	
	graphNLOMuHalfEta7TeV->Draw("same,c");
	graphNLOMuOneEta7TeV->Draw("same,c");
	graphNLOMuTwoEta7TeV->Draw("same,c");

	labelLegendAVar2->Draw();
	labelScalingPi07TeV->Draw();
	labelScalingPi0900GeV->Draw();
	labelScalingEta7TeV->Draw();
	DrawNormalizationErrorText(normalizationInvX3En[0],normalizationInvX3En[1],normalizationInvX3En[2],normalizationInvX3En[3],normalizationInvX3En[4],"No2.76"); 
	
		TPad* padXSectionNLOVar2LegendSep = new TPad("padXSectionNLOVar2LegendSep", "", 0.17, 0.02, 0.95, 0.205,-1, -1, -2); 
		DrawGammaPadSettings( padXSectionNLOVar2LegendSep, 0., 0., 0., 0.);
		padXSectionNLOVar2LegendSep->Draw();
		padXSectionNLOVar2LegendSep->cd();
		
		//*************** first Column **********************************************************
		textSpectrumVar->Draw();
		textFitCombVar2->Draw();
		textNLOMuHalfVar2->Draw();
		textNLOMuOneVar2->Draw();
		textNLOMuTwoVar2->Draw();
		textNLOMuTwoBKKVar2->Draw();
		
		
		//*************** second Column **********************************************************
		textPi07TeVNLOVar2->Draw();
		TLatex *textPi07TeVNLOsysVar2Sep = new TLatex(columnsNLOLegendVar2[1],rowsNLOLegendVar2[1],"syst., stat.");
		SetStyleTLatex( textPi07TeVNLOsysVar2Sep, textSizeTopLowerLablesVar2,4);
		textPi07TeVNLOsysVar2Sep->Draw();
		TBox* boxCombinedPi07TeVVar2Sep = CreateBoxFromGraphWithFill(graphInvCrossSectionPi0Comb7TeVSysErr,columnsNLOLegendVar2[1]+offsetNLOLegendMarkerXVar2-offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2- offsetNLOLegendBoxVar2, columnsNLOLegendVar2[1]+offsetNLOLegendMarkerXVar2+offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2+offsetNLOLegendBoxVar2);
		boxCombinedPi07TeVVar2Sep->Draw("f");
		markerCombinedPi07TeVVar2->DrawMarker(columnsNLOLegendVar2[1]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		lineFit7TeVNLOVar2->Draw("same");
		lineNLOPi07TeVMuHalfVar2->Draw("same");
		lineNLOPi07TeVMuOneVar2->Draw("same");
		lineNLOPi07TeVMuTwoVar2->Draw("same");
		lineNLOPi07TeVMuTwoBKKVar2->Draw("same");
// 		lineNLOPi07TeVMuTwoDSSVar2->Draw("same");		
		//*************** third Column **********************************************************
		textPi0900GeVNLOVar2->Draw();
		TLatex *textPi0900GeVNLOsysVar2Sep = new TLatex(columnsNLOLegendVar2[2],rowsNLOLegendVar2[1],"syst., stat.");
		SetStyleTLatex( textPi0900GeVNLOsysVar2Sep, textSizeTopLowerLablesVar2,4);
		textPi0900GeVNLOsysVar2Sep->Draw();
		TBox* boxCombinedPi0900GeVVar2Sep = CreateBoxFromGraphWithFill(graphInvCrossSectionPi0Comb900GeVSysErr,columnsNLOLegendVar2[2]+offsetNLOLegendMarkerXVar2-offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2- offsetNLOLegendBoxVar2, columnsNLOLegendVar2[2]+offsetNLOLegendMarkerXVar2+offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2+offsetNLOLegendBoxVar2);
		boxCombinedPi0900GeVVar2Sep->Draw("f");
		markerCombinedPi0900GeVVar2->DrawMarker(columnsNLOLegendVar2[2]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		lineFitPi0900GeVNLOVar2->Draw("same");
		lineNLOPi0900GeVMuHalfVar2->Draw("same");
		lineNLOPi0900GeVMuOneVar2->Draw("same");
		lineNLOPi0900GeVMuTwoVar2->Draw("same");
		lineNLOPi0900GeVMuTwoBKKVar2->Draw("same");
		
		//**************** forth Column **********************************************************
		textEta7TeVNLOVar2->Draw();
		TLatex *textEta7TeVNLOsysVar2Sep = new TLatex(columnsNLOLegendVar2[3],rowsNLOLegendVar2[1],"syst., stat.");
		SetStyleTLatex( textEta7TeVNLOsysVar2Sep, textSizeTopLowerLablesVar2,4);
		textEta7TeVNLOsysVar2Sep->Draw();
		TBox* boxCombinedEta7TeVVar2Sep = CreateBoxFromGraphWithFill(graphInvCrossSectionEtaComb7TeVSysErr,columnsNLOLegendVar2[3]+offsetNLOLegendMarkerXVar2-offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2- offsetNLOLegendBoxVar2, columnsNLOLegendVar2[3]+offsetNLOLegendMarkerXVar2+offsetNLOLegendLineVar2, rowsNLOLegendVar2[2]+ offsetNLOLegendMarkerYVar2+offsetNLOLegendBoxVar2);
		boxCombinedEta7TeVVar2Sep->Draw("f");
		markerCombinedEta7TeVVar2->DrawMarker(columnsNLOLegendVar2[3]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		lineFitEta7TeVNLOVar2->Draw("same");
		lineNLOEta7TeVMuHalfVar2->Draw("same");
		lineNLOEta7TeVMuOneVar2->Draw("same");
		lineNLOEta7TeVMuTwoVar2->Draw("same");
	//************************************************* End Legend ***************************************************
	
	padXSectionNLOVar2SepRatioPi07TeV->cd();
	padXSectionNLOVar2SepRatioPi07TeV->SetLogx();
	
	TH2F * ratio2DInvXSectionNLOVar2SepPi0 = new TH2F("ratio2DInvXSectionNLOVar2SepPi0","ratio2DInvXSectionNLOVar2SepPi0",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOVar2SepPi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.11,0.115, 1,0.53, 512, 505);
	ratio2DInvXSectionNLOVar2SepPi0->GetYaxis()->SetMoreLogLabels(kTRUE);
	ratio2DInvXSectionNLOVar2SepPi0->GetYaxis()->SetNdivisions(515);
	ratio2DInvXSectionNLOVar2SepPi0->GetYaxis()->SetNoExponent(kTRUE);
	ratio2DInvXSectionNLOVar2SepPi0->GetXaxis()->SetTickLength(0.07);
	ratio2DInvXSectionNLOVar2SepPi0->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat, markerStyleCommmonSpectrumPi07TeV,markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kFALSE);
	graphRatioCombCombFitStat->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleCommmonSpectrumPi07TeV,markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE,colorCommonSpectrumPi07TeVBox);
	graphRatioCombCombFitSys->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFitSys->Draw("2,same");
	graphRatioCombCombFitStat->Draw("p,same");

	graphRatioCombNLOPi07TeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi07TeVMuOne->Draw("same,c");
	graphRatioCombNLOPi07TeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi07TeVMuTwo->Draw("same,c");
// 	graphRatioCombNLODSSPi07TeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi07TeVRatio->Draw();
	labelLegendBVar2->Draw();
		
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOVar2SepRatioPi0900GeV->cd();
	padXSectionNLOVar2SepRatioPi0900GeV->SetLogx();
	TH2F * ratio2DInvXSectionNLOVar2SepPi0900GeV;
	ratio2DInvXSectionNLOVar2SepPi0900GeV = new TH2F("ratio2DInvXSectionNLOVar2SepPi0900GeV","ratio2DInvXSectionNLOVar2SepPi0900GeV",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOVar2SepPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.11,0.115, 1,0.53, 512, 505);
	ratio2DInvXSectionNLOVar2SepPi0900GeV->GetYaxis()->SetMoreLogLabels(kTRUE);
	ratio2DInvXSectionNLOVar2SepPi0900GeV->GetYaxis()->SetNdivisions(515);
	ratio2DInvXSectionNLOVar2SepPi0900GeV->GetYaxis()->SetNoExponent(kTRUE);
	ratio2DInvXSectionNLOVar2SepPi0900GeV->GetXaxis()->SetTickLength(0.07);
	ratio2DInvXSectionNLOVar2SepPi0900GeV->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit900GeVStat, markerStyleCommmonSpectrumPi0900GeV,markerSizeCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kFALSE);
	graphRatioCombCombFit900GeVStat->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit900GeVSys, markerStyleCommmonSpectrumPi0900GeV,markerSizeCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi0900GeVBox);
	graphRatioCombCombFit900GeVSys->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFit900GeVSys->Draw("2,same");
	graphRatioCombCombFit900GeVStat->Draw("p,same");

	graphRatioCombNLOPi0900GeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuOne->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi0900GeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi0900GeVRatio->Draw();
	
	boxErrorSigmaPi0900GeVRatio->Draw();	
	labelLegendCVar2->Draw();
		
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOVar2SepRatioEta7TeV->cd();
	padXSectionNLOVar2SepRatioEta7TeV->SetLogx();
// 	padXSectionNLOVar2SepRatioEta7TeV->SetLogy();
	TH2F * ratio2DInvXSectionNLOVar2SepEta=  new TH2F("ratio2DInvXSectionNLOVar2SepEta","ratio2DInvXSectionNLOVar2SepEta",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOVar2SepEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.09,0.095, 0.09,0.095, 0.95,0.63, 510, 505);
	ratio2DInvXSectionNLOVar2SepEta->GetYaxis()->SetMoreLogLabels(kTRUE);
 	ratio2DInvXSectionNLOVar2SepEta->GetYaxis()->SetNdivisions(512);
	ratio2DInvXSectionNLOVar2SepEta->GetYaxis()->SetNoExponent(kTRUE);
	ratio2DInvXSectionNLOVar2SepEta->GetXaxis()->SetMoreLogLabels(kTRUE);
	ratio2DInvXSectionNLOVar2SepEta->GetXaxis()->SetNoExponent(kTRUE);
// 	ratio2DInvXSectionNLOVar2SepEta->GetXaxis()->SetLabelOffset(-0.03);
	ratio2DInvXSectionNLOVar2SepEta->GetXaxis()->SetTickLength(0.07);
	ratio2DInvXSectionNLOVar2SepEta->DrawCopy(); 

	boxErrorSigmaEta7TeVRatio->Draw();
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta7TeVStat,markerStyleCommmonSpectrumEta7TeV,markerSizeCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, widthCommonSpectrumBoxes, kFALSE);
	graphRatioCombCombFitEta7TeVStat->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta7TeVSys,markerStyleCommmonSpectrumEta7TeV,markerSizeCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, colorCommonSpectrumEta7TeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumEta7TeVBox);
	graphRatioCombCombFitEta7TeVSys->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFitEta7TeVSys->Draw("2,same");
	graphRatioCombCombFitEta7TeVStat->Draw("p,same");
	
	graphRatioCombNLOEta7TeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta7TeVMuOne->Draw("same,c");
	graphRatioCombNLOEta7TeVMuTwo->Draw("same,c");
	
	
	labelLegendDVar2->Draw();
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	canvasInvXSectionNLOVar2Sep->Update();
	canvasInvXSectionNLOVar2Sep->Print(Form("%s/InvXSectionNLOVar2Sep_Paper.%s",outputDir.Data(),suffix.Data()));


		//*****************************************************************************************************************
	//************************ Pi0 Meson Combined + NLO just spectrum******************************************************
	//*****************************************************************************************************************

	TCanvas* canvasInvXSectionNLOOnlySpectraMixed = new TCanvas("canvasInvXSectionNLOOnlySpectraMixed","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlySpectraMixed,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionNLOOnlySpectraMixed = new TPad("padComparisonXSectionNLOOnlySpectraMixed", "", 0., 0., 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionNLOOnlySpectraMixed, 0.15, 0.02, 0.02, 0.09);
	padComparisonXSectionNLOOnlySpectraMixed->Draw();
	
	TPad* padXSectionNLOLegendOnlySpectraMixed = new TPad("padXSectionNLOLegendOnlySpectraMixed", "", 0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraMixed, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlySpectraMixed->Draw();
	
	padComparisonXSectionNLOOnlySpectraMixed->cd();
	padComparisonXSectionNLOOnlySpectraMixed->SetLogy();		
	padComparisonXSectionNLOOnlySpectraMixed->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionNLOOnlySpectraMixed = new TH2F("histo2DInvXSectionNLOOnlySpectraMixed","histo2DInvXSectionNLOOnlySpectraMixed",1000,0.23,30.,1000,2e-1,1e12);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLOOnlySpectraMixed, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionNLOOnlySpectraMixed->DrawCopy(); 
	
	graphInvCrossSectionPi0Comb900GeVSysErr->Draw("2,same");	
	graphInvCrossSectionPi0Comb900GeVStatErr->Draw("p,same");	
	
	graphInvCrossSectionPi0Comb7TeVSysErr->Draw("2,same");
	graphInvCrossSectionPi0Comb7TeVStatErr->Draw("p,same");
	
	graphInvCrossSectionEtaComb7TeVSysErr->Draw("2,same");
	graphInvCrossSectionEtaComb7TeVStatErr->Draw("p,same");

	fitInvCrossSectionPi0->SetLineWidth(widthCommonFit*1.5);
	fitInvCrossSectionPi0->Draw("same");
	histoFitInvCrossSectionPi0900GeV->SetLineWidth(widthCommonFit*1.5);
	histoFitInvCrossSectionPi0900GeV->Draw("same,c");
	histoFitInvCrossSectionEta7TeV->SetLineWidth(widthCommonFit*1.5);
	histoFitInvCrossSectionEta7TeV->Draw("same,c");
	graphNLOMuHalfPi07TeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOMuHalfPi07TeV->Draw("same,c");
	graphNLOMuOnePi07TeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOMuOnePi07TeV->Draw("same,c");
	graphNLOMuTwoPi07TeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOMuTwoPi07TeV->Draw("same,c");
	DrawGammaNLOTGraph( graphNLOBKKCalcMuTwo7TeV, 0.5*widthCommonFit, styleLineNLOMuTwoBKK, colorNLOPi07TeVMuTwo);
	graphNLOBKKCalcMuTwo7TeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOBKKCalcMuTwo7TeV->Draw("same,c");
	graphNLOMuHalfPi0900GeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOMuHalfPi0900GeV->Draw("same,c");
	graphNLOMuOnePi0900GeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOMuOnePi0900GeV->Draw("same,c");
	graphNLOMuTwoPi0900GeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOMuTwoPi0900GeV->Draw("same,c");
	graphNLOBKKCalcMuTwo900GeV= ScaleGraph(graphNLOBKKCalcMuTwo900GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOBKKCalcMuTwo900GeV, 0.5*widthCommonFit, styleLineNLOMuTwoBKK, colorNLOPi0900GeVMuTwo);
	graphNLOBKKCalcMuTwo900GeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOBKKCalcMuTwo900GeV->Draw("same,c");
	
	graphNLOMuHalfEta7TeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOMuHalfEta7TeV->Draw("same,c");
	graphNLOMuOneEta7TeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOMuOneEta7TeV->Draw("same,c");
	graphNLOMuTwoEta7TeV->SetLineWidth(widthCommonFit*1.5);
	graphNLOMuTwoEta7TeV->Draw("same,c");

	labelLegendAVar2->Draw();
	labelScalingPi07TeV->Draw();
	labelScalingPi0900GeV->Draw();
	labelScalingEta7TeV->Draw();
	DrawNormalizationErrorText(normalizationInvX3En[0],normalizationInvX3En[1],normalizationInvX3En[2],normalizationInvX3En[3]*1.5,normalizationInvX3En[4]*1.5,"No2.76"); 
	
	
	//************************************************* Begin Legend ***************************************************
	
	padXSectionNLOLegendOnlySpectraMixed->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraMixed, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlySpectraMixed->SetBorderMode(-1);
	padXSectionNLOLegendOnlySpectraMixed->SetBorderSize(3);
	padXSectionNLOLegendOnlySpectraMixed->Draw();
	padXSectionNLOLegendOnlySpectraMixed->cd();

		textSpectrumVar->Draw();
		textFitCombVar2->Draw();
		textNLOMuHalfVar2->Draw();
		textNLOMuOneVar2->Draw();
		textNLOMuTwoVar2->Draw();
		textNLOMuTwoBKKVar2->Draw();
		
		
		//*************** second Column **********************************************************
		textPi07TeVNLOVar2->Draw();
		textPi07TeVNLOsysVar2Sep->Draw();
		boxCombinedPi07TeVVar2Sep->Draw("f");
		markerCombinedPi07TeVVar2->DrawMarker(columnsNLOLegendVar2[1]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		lineFit7TeVNLOVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineFit7TeVNLOVar2->Draw("same");
		lineNLOPi07TeVMuHalfVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOPi07TeVMuHalfVar2->Draw("same");
		lineNLOPi07TeVMuOneVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOPi07TeVMuOneVar2->Draw("same");
		lineNLOPi07TeVMuTwoVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOPi07TeVMuTwoVar2->Draw("same");
		lineNLOPi07TeVMuTwoBKKVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOPi07TeVMuTwoBKKVar2->Draw("same");
// 		lineNLOPi07TeVMuTwoDSSVar2->Draw("same");		
		//*************** third Column **********************************************************
		textPi0900GeVNLOVar2->Draw();
		textPi0900GeVNLOsysVar2Sep->Draw();
		boxCombinedPi0900GeVVar2Sep->Draw("f");
		markerCombinedPi0900GeVVar2->DrawMarker(columnsNLOLegendVar2[2]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		lineFitPi0900GeVNLOVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineFitPi0900GeVNLOVar2->Draw("same");
		lineNLOPi0900GeVMuHalfVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOPi0900GeVMuHalfVar2->Draw("same");
		lineNLOPi0900GeVMuOneVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOPi0900GeVMuOneVar2->Draw("same");
		lineNLOPi0900GeVMuTwoVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOPi0900GeVMuTwoVar2->Draw("same");
		lineNLOPi0900GeVMuTwoBKKVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOPi0900GeVMuTwoBKKVar2->Draw("same");
		
		//**************** forth Column **********************************************************
		textEta7TeVNLOVar2->Draw();
		textEta7TeVNLOsysVar2Sep->Draw();
		boxCombinedEta7TeVVar2Sep->Draw("f");
		markerCombinedEta7TeVVar2->DrawMarker(columnsNLOLegendVar2[3]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		lineFitEta7TeVNLOVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineFitEta7TeVNLOVar2->Draw("same");
		lineNLOEta7TeVMuHalfVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOEta7TeVMuHalfVar2->Draw("same");
		lineNLOEta7TeVMuOneVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOEta7TeVMuOneVar2->Draw("same");
		lineNLOEta7TeVMuTwoVar2->SetLineWidth(widthCommonFit*1.5*scaleLineWidthNLOVar2);
		lineNLOEta7TeVMuTwoVar2->Draw("same");
	//lineNLOPi0900GeVMuTwoBKKALLEnergies->Draw("same");

	canvasInvXSectionNLOOnlySpectraMixed->Update();
	canvasInvXSectionNLOOnlySpectraMixed->Print(Form("%s/Mixed_InvXSectionNLO_OnlySpectrum_Paper.%s",outputDir.Data(),suffix.Data()));
	
// //**************************************************************************************************************************************
// //************************************ Inv Crosssection + NLO only ratio Mixed ***************************************************************
// //**************************************************************************************************************************************
// 	
	TCanvas* canvasInvXSectionNLOOnlyRatioMixed = new TCanvas("canvasInvXSectionNLOOnlyRatioMixed","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlyRatioMixed,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padXSectionNLOLegendOnlyRatioMixed = new TPad("padXSectionNLOLegendOnlyRatioMixed", "",  0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioMixed, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlyRatioMixed->Draw();

	TPad* padXSectionNLOOnlyRatioMixedPi07TeV = new TPad("padXSectionNLOOnlyRatioMixedPi07TeV", "", 0., 0.565, 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioMixedPi07TeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOOnlyRatioMixedPi07TeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioMixedPi02760GeV = new TPad("padXSectionNLOOnlyRatioMixedPi02760GeV", "", 0., 0.33, 1., 0.565,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioMixedPi02760GeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOOnlyRatioMixedPi02760GeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioMixedPi0900GeV = new TPad("padXSectionNLOOnlyRatioMixedPi0900GeV", "", 0., 0., 1., 0.33,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioMixedPi0900GeV, 0.15, 0.02, 0., 0.28);
	padXSectionNLOOnlyRatioMixedPi0900GeV->Draw();
	
	padXSectionNLOOnlyRatioMixedPi07TeV->cd();
	padXSectionNLOOnlyRatioMixedPi07TeV->SetLogx();
	
	ratio2DInvXSectionNLOVar2SepPi0->DrawCopy(); 
	graphRatioCombCombFitSys->Draw("2,same");
	graphRatioCombCombFitStat->Draw("p,same");
	graphRatioCombNLOPi07TeVMuHalf->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOPi07TeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi07TeVMuOne->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOPi07TeVMuOne->Draw("same,c");
	graphRatioCombNLOPi07TeVMuTwo->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOPi07TeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi07TeVMuTwo->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOBKKPi07TeVMuTwo->Draw("same,c");
// 	graphRatioCombNLODSSPi07TeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi07TeVRatio->Draw();
	TLatex *labelLegendBVar2_2 = new TLatex(0.18,0.86,"b)");
	SetStyleTLatex( labelLegendBVar2_2, 0.12,4);
	labelLegendBVar2_2->Draw();
		
	DrawGammaLines(0., 30.,1., 1.,0.2);
	
	
	padXSectionNLOOnlyRatioMixedPi02760GeV->cd();
	padXSectionNLOOnlyRatioMixedPi02760GeV->SetLogx();
	
	ratio2DInvXSectionNLOVar2SepPi0900GeV->DrawCopy(); 
	graphRatioCombCombFit900GeVSys->Draw("2,same");
	graphRatioCombCombFit900GeVStat->Draw("p,same");
	graphRatioCombNLOPi0900GeVMuHalf->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOPi0900GeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuOne->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOPi0900GeVMuOne->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuTwo->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOPi0900GeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi0900GeVMuTwo->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOBKKPi0900GeVMuTwo->Draw("same,c");
	boxErrorSigmaPi0900GeVRatio->Draw();	
	TLatex *labelLegendCVar2_2 = new TLatex(0.18,0.86,"c)");
	SetStyleTLatex( labelLegendCVar2_2, 0.12,4);
	labelLegendCVar2_2->Draw();
	
		
	DrawGammaLines(0., 30.,1., 1.,0.2);
	
	padXSectionNLOOnlyRatioMixedPi0900GeV->cd();
	padXSectionNLOOnlyRatioMixedPi0900GeV->SetLogx();
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOVar2SepEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.081,0.086, 0.081,0.086, 0.95,0.7, 510, 512);
	ratio2DInvXSectionNLOVar2SepEta->GetXaxis()->SetMoreLogLabels(kFALSE);
	ratio2DInvXSectionNLOVar2SepEta->DrawCopy(); 

	boxErrorSigmaEta7TeVRatio->Draw();
	
	
	graphRatioCombCombFitEta7TeVSys->Draw("2,same");
	graphRatioCombCombFitEta7TeVStat->Draw("p,same");
	
	graphRatioCombNLOEta7TeVMuHalf->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOEta7TeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta7TeVMuOne->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOEta7TeVMuOne->Draw("same,c");
	graphRatioCombNLOEta7TeVMuTwo->SetLineWidth(widthCommonFit*1.5);
	graphRatioCombNLOEta7TeVMuTwo->Draw("same,c");
	
	TLatex *labelLegendDVar2_2 = new TLatex(0.18,0.88,"d)");
	SetStyleTLatex( labelLegendDVar2_2, 0.09,4);
	labelLegendDVar2_2->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.2);
		
	padXSectionNLOLegendOnlyRatioMixed->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioMixed, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlyRatioMixed->SetBorderMode(-1);
	padXSectionNLOLegendOnlyRatioMixed->SetBorderSize(3);
	padXSectionNLOLegendOnlyRatioMixed->Draw();
	padXSectionNLOLegendOnlyRatioMixed->cd();

			textSpectrumVar->Draw();
		textFitCombVar2->Draw();
		textNLOMuHalfVar2->Draw();
		textNLOMuOneVar2->Draw();
		textNLOMuTwoVar2->Draw();
		textNLOMuTwoBKKVar2->Draw();
		
		
		//*************** second Column **********************************************************
		textPi07TeVNLOVar2->Draw();
		textPi07TeVNLOsysVar2Sep->Draw();
		boxCombinedPi07TeVVar2Sep->Draw("f");
		markerCombinedPi07TeVVar2->DrawMarker(columnsNLOLegendVar2[1]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		lineFit7TeVNLOVar2->Draw("same");
		lineNLOPi07TeVMuHalfVar2->Draw("same");
		lineNLOPi07TeVMuOneVar2->Draw("same");
		lineNLOPi07TeVMuTwoVar2->Draw("same");
		lineNLOPi07TeVMuTwoBKKVar2->Draw("same");
// 		lineNLOPi07TeVMuTwoDSSVar2->Draw("same");		
		//*************** third Column **********************************************************
		textPi0900GeVNLOVar2->Draw();
		textPi0900GeVNLOsysVar2Sep->Draw();
		boxCombinedPi0900GeVVar2Sep->Draw("f");
		markerCombinedPi0900GeVVar2->DrawMarker(columnsNLOLegendVar2[2]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		lineFitPi0900GeVNLOVar2->Draw("same");
		lineNLOPi0900GeVMuHalfVar2->Draw("same");
		lineNLOPi0900GeVMuOneVar2->Draw("same");
		lineNLOPi0900GeVMuTwoVar2->Draw("same");
		lineNLOPi0900GeVMuTwoBKKVar2->Draw("same");
		
		//**************** forth Column **********************************************************
		textEta7TeVNLOVar2->Draw();
		textEta7TeVNLOsysVar2Sep->Draw();
		boxCombinedEta7TeVVar2Sep->Draw("f");
		markerCombinedEta7TeVVar2->DrawMarker(columnsNLOLegendVar2[3]+offsetNLOLegendMarkerXVar2,rowsNLOLegendVar2[2]+offsetNLOLegendMarkerYVar2);
		lineFitEta7TeVNLOVar2->Draw("same");
		lineNLOEta7TeVMuHalfVar2->Draw("same");
		lineNLOEta7TeVMuOneVar2->Draw("same");
		lineNLOEta7TeVMuTwoVar2->Draw("same");

	
	canvasInvXSectionNLOOnlyRatioMixed->Update();
	canvasInvXSectionNLOOnlyRatioMixed->Print(Form("%s/Mixed_InvXSectionNLO_OnlyRatio_Paper.%s",outputDir.Data(),suffix.Data()));	
	
	//***************************************************************************************************************
	//************************ Pi0 Combined + NLO *******************************************************************
	//***************************************************************************************************************
	
	TCanvas* canvasInvXSectionALLEnergies = new TCanvas("canvasInvXSectionALLEnergies","",200,10,1200,2000);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionALLEnergies,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionALLEnergies = new TPad("padComparisonXSectionALLEnergies", "", 0., 0.42, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionALLEnergies, 0.15, 0.02, 0.03, 0.);
	padComparisonXSectionALLEnergies->Draw();
	
	TPad* padXSectionALLEnergiesRatioPi07TeV = new TPad("padXSectionALLEnergiesRatioPi07TeV", "", 0., 0.30, 1., 0.42,-1, -1, -2);
	DrawGammaPadSettings( padXSectionALLEnergiesRatioPi07TeV, 0.15, 0.02, 0., 0.);
	padXSectionALLEnergiesRatioPi07TeV->Draw();
	
	TPad* padXSectionALLEnergiesRatioPi02760GeV = new TPad("padXSectionALLEnergiesRatioPi02760GeV", "", 0., 0.18, 1., 0.30,-1, -1, -2);
	DrawGammaPadSettings( padXSectionALLEnergiesRatioPi02760GeV, 0.15, 0.02, 0., 0.);
	padXSectionALLEnergiesRatioPi02760GeV->Draw();
	
	TPad* padXSectionALLEnergiesRatioPi0900GeV = new TPad("padXSectionALLEnergiesRatioPi0900GeV", "", 0., 0., 1., 0.18,-1, -1, -2);
	DrawGammaPadSettings( padXSectionALLEnergiesRatioPi0900GeV, 0.15, 0.02, 0., 0.28);
	padXSectionALLEnergiesRatioPi0900GeV->Draw();
	
	
	padComparisonXSectionALLEnergies->cd();
	padComparisonXSectionALLEnergies->SetLogy();		
	padComparisonXSectionALLEnergies->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionALLEnergies = new TH2F("histo2DInvXSectionALLEnergies","histo2DInvXSectionALLEnergies",1000,0.23,30.,1000,2e-3,10e11);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionALLEnergies, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionALLEnergies->GetXaxis()->SetLabelFont(42);
	histo2DInvXSectionALLEnergies->GetYaxis()->SetLabelFont(42); 
	histo2DInvXSectionALLEnergies->GetXaxis()->SetTitleFont(62);
	histo2DInvXSectionALLEnergies->GetYaxis()->SetTitleFont(62);
	histo2DInvXSectionALLEnergies->DrawCopy(); 
	
	graphInvCrossSectionPi0Comb900GeV = ScaleGraph(graphInvCrossSectionPi0Comb900GeV,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb900GeV, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionPi0Comb900GeV->SetLineWidth(widthCommonErrors);
	//graphInvCrossSectionPi0Comb900GeV->Draw("same,p");
	graphInvCrossSectionPi0Comb900GeV->Draw("p,E2same");
	graphInvCrossSectionPi0Comb2760GeV = ScaleGraph(graphInvCrossSectionPi0Comb2760GeV,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb2760GeV, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionPi0Comb2760GeV->SetLineWidth(widthCommonErrors);
	//graphInvCrossSectionPi0Comb2760GeV->Draw("same,p");
	graphInvCrossSectionPi0Comb2760GeV->Draw("p,E2same");
	
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb7TeV, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionPi0Comb7TeV->SetLineWidth(widthCommonErrors);
	//graphInvCrossSectionPi0Comb7TeV->Draw("p,same");
	graphInvCrossSectionPi0Comb7TeV->Draw("p,E2same");
		
	DrawGammaSetMarkerTF1( fitInvCrossSectionPi0, styleFitCommonSpectrum, widthCommonFit, colorCommonSpectrumPi07TeV);
	fitInvCrossSectionPi0->Draw("same");
	
	histoFitInvCrossSectionPi0900GeV->Scale(1e-1);
	SetStyleHisto(histoFitInvCrossSectionPi0900GeV, widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi0900GeV);
	histoFitInvCrossSectionPi0900GeV->Draw("same,c");
	
	histoFitInvCrossSectionPi02760GeV->Scale(1e-1);
	SetStyleHisto(histoFitInvCrossSectionPi02760GeV, widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi02760GeV);
	histoFitInvCrossSectionPi02760GeV->Draw("same,c");

	DrawGammaNLOTGraph( graphNLOMuHalfPi07TeV, widthCommonFit, styleLineNLOMuHalf, colorNLOPi07TeVMuHalf);
	graphNLOMuHalfPi07TeV->Draw("same,c");
	DrawGammaNLOTGraph( graphNLOMuOnePi07TeV, widthCommonFit, styleLineNLOMuOne, colorNLOPi07TeVMuOne);
	graphNLOMuOnePi07TeV->Draw("same,c");
	DrawGammaNLOTGraph( graphNLOMuTwoPi07TeV, widthCommonFit, styleLineNLOMuTwo, colorNLOPi07TeVMuTwo);
	graphNLOMuTwoPi07TeV->Draw("same,c");
	
	graphNLOMuHalfPi02760GeV= ScaleGraph(graphNLOMuHalfPi02760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuHalfPi02760GeV, widthCommonFit, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	graphNLOMuHalfPi02760GeV->Draw("same,c");
	graphNLOMuOnePi02760GeV= ScaleGraph(graphNLOMuOnePi02760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuOnePi02760GeV, widthCommonFit, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	graphNLOMuOnePi02760GeV->Draw("same,c");
	graphNLOMuTwoPi02760GeV = ScaleGraph(graphNLOMuTwoPi02760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuTwoPi02760GeV, widthCommonFit, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);
	graphNLOMuTwoPi02760GeV->Draw("same,c");
	
	
	graphNLOMuHalfPi0900GeV= ScaleGraph(graphNLOMuHalfPi0900GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuHalfPi0900GeV, widthCommonFit, styleLineNLOMuHalf, colorNLOPi0900GeVMuHalf);
	graphNLOMuHalfPi0900GeV->Draw("same,c");
	graphNLOMuOnePi0900GeV= ScaleGraph(graphNLOMuOnePi0900GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuOnePi0900GeV, widthCommonFit, styleLineNLOMuOne, colorNLOPi0900GeVMuOne);
	graphNLOMuOnePi0900GeV->Draw("same,c");
	graphNLOMuTwoPi0900GeV = ScaleGraph(graphNLOMuTwoPi0900GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuTwoPi0900GeV, widthCommonFit, styleLineNLOMuTwo, colorNLOPi0900GeVMuTwo);
	graphNLOMuTwoPi0900GeV->Draw("same,c");
	
	labelScalingPi07TeV->Draw();
	TLatex *labelScalingPi02760GeV = new TLatex(0.27,5E9,"x 10^{-1}");
	SetStyleTLatex( labelScalingPi02760GeV, 0.025,4,histoFitInvCrossSectionPi02760GeV->GetLineColor(),62,kFALSE);
	labelScalingPi02760GeV->Draw();
	TLatex *labelScalingPi0900GeVALLEnergies = new TLatex(0.27,1E8,"x 10^{-2}");
	SetStyleTLatex( labelScalingPi0900GeVALLEnergies, 0.025,4,histoFitInvCrossSectionPi0900GeV->GetLineColor(),62,kFALSE);
	labelScalingPi0900GeVALLEnergies->Draw();	
	
	DrawNormalizationErrorText(normalizationInvX1MesonALLEn[0],normalizationInvX1MesonALLEn[1],normalizationInvX1MesonALLEn[2],
						  normalizationInvX1MesonALLEn[3],normalizationInvX1MesonALLEn[4],"all"); 

	TPad* padXSectionALLEnergiesLegend = new TPad("padXSectionALLEnergiesLegend", "", 0.18, 0.02, 0.98, 0.22,-1, -1, -2); 
	DrawGammaPadSettings( padXSectionALLEnergiesLegend, 0., 0., 0., 0.);
	padXSectionALLEnergiesLegend->SetFillStyle(0);	
	padXSectionALLEnergiesLegend->Draw();
	padXSectionALLEnergiesLegend->cd();
	
	//**************** Row def ************************
		Double_t rowsNLOLegendPrelAndFinal[8] = {0.9,0.81,0.70,0.59,0.48,0.36,0.23,0.11};
		//*************** Label sizes *********************
 		Double_t textSizeLeftLabelsPrelAndFinal = 0.11;
		Double_t textSizeTopLablesPrelAndFinal = 0.115;
 		Double_t textSizeTopLowerLablesPrelAndFinal = 0.11;
		//*************** Column def ***********************
		Double_t columnsNLOLegendPrelAndFinal[4] = {0.,0.23,0.46,0.74};
		//*************** Size factors ********************
// 		Double_t scaleMarkerNLO = 1.5;
// 		Double_t scaleLineWidthNLO = 1.;
		//*************** Offsets *************************
		Double_t offsetNLOLegendPrelAndFinalMarkerX = 0.09;
		Double_t offsetNLOLegendPrelAndFinalMarkerY = 0.03;
		Double_t offsetNLOLegendPrelAndFinalBox = 0.05;
		Double_t offsetNLOLegendPrelAndFinalLine = 0.06;

	
	//*************** first Column **********************************************************
	TLatex *textSpectrumALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[2],"combined Spec.");
	SetStyleTLatex( textSpectrumALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	textSpectrumALLEnergies->Draw();
	TLatex *textFitCombALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[3],"fit combined");
	SetStyleTLatex( textFitCombALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	textFitCombALLEnergies->Draw();
	TLatex *textNLOMuHalfALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[4],"NLO #mu= 0.5 #it{p}_{T}");
	SetStyleTLatex( textNLOMuHalfALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	textNLOMuHalfALLEnergies->Draw();
	TLatex *textNLOMuOneALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[5],"NLO #mu= #it{p}_{T}");
	SetStyleTLatex( textNLOMuOneALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	textNLOMuOneALLEnergies->Draw();
	TLatex *textNLOMuTwoALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[6],"NLO #mu= 2 #it{p}_{T}");
	SetStyleTLatex( textNLOMuTwoALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	textNLOMuTwoALLEnergies->Draw();
	TLatex *textNLOMuTwoBKKALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[7],"NLO #mu= 2 #it{p}_{T} (BKK)");
	SetStyleTLatex( textNLOMuTwoBKKALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	textNLOMuTwoBKKALLEnergies->Draw();
	
	//*************** second Column **********************************************************
	TLatex *textPi07TeVNLOALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[1],rowsNLOLegendPrelAndFinal[0],"#pi^{0}, #sqrt{#it{s}} = 7 TeV (*)");
	SetStyleTLatex( textPi07TeVNLOALLEnergies, textSizeTopLablesPrelAndFinal,4);
	textPi07TeVNLOALLEnergies->Draw();
	TLatex *textPi07TeVNLOsysALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[1],rowsNLOLegendPrelAndFinal[1],"syst. + stat.");
	SetStyleTLatex( textPi07TeVNLOsysALLEnergies, textSizeTopLowerLablesPrelAndFinal,4);
	textPi07TeVNLOsysALLEnergies->Draw();
	TBox* boxCombinedPi07TeVALLEnergies = CreateBoxFromGraph(graphInvCrossSectionPi0Comb7TeV,columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi07TeVALLEnergies->Draw("l");
	TMarker* markerCombinedPi07TeVALLEnergies = CreateMarkerFromGraph(graphInvCrossSectionPi0Comb7TeV,columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	TLine * lineFit7TeVNLOALLEnergies = CreateLineFromFit(fitInvCrossSectionPi0, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO);
	lineFit7TeVNLOALLEnergies->Draw("same");
	TLine * lineNLOPi07TeVMuHalfALLEnergies = CreateLineFromGraph(graphNLOMuHalfPi07TeV, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[4] + offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
	TLine * lineNLOPi07TeVMuOneALLEnergies = CreateLineFromGraph(graphNLOMuOnePi07TeV,columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
	TLine * lineNLOPi07TeVMuTwoALLEnergies = CreateLineFromGraph(graphNLOMuTwoPi07TeV, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[6]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[6]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
	TLine * lineNLOPi07TeVMuTwoBKKALLEnergies = CreateLineFromGraph(graphRatioCombNLOBKKPi07TeVMuTwo, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[7]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[7]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi07TeVMuTwoBKKALLEnergies->Draw("same");
	
	//*************** third Column **********************************************************
	TLatex *textPi02760GeVNLOALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[2],rowsNLOLegendPrelAndFinal[0],"#pi^{0}, #sqrt{#it{s}} = 2.76 TeV (**)");
	SetStyleTLatex( textPi02760GeVNLOALLEnergies, textSizeTopLablesPrelAndFinal,4);
	textPi02760GeVNLOALLEnergies->Draw();
	TLatex *textPi02760GeVNLOsysALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[2],rowsNLOLegendPrelAndFinal[1],"syst. + stat.");
	SetStyleTLatex( textPi02760GeVNLOsysALLEnergies, textSizeTopLowerLablesPrelAndFinal,4);
	textPi02760GeVNLOsysALLEnergies->Draw();
	TBox* boxCombinedPi02760GeVALLEnergies = CreateBoxFromGraph(graphInvCrossSectionPi0Comb2760GeV,columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi02760GeVALLEnergies->Draw("l");
	TMarker* markerCombinedPi02760GeVALLEnergies = CreateMarkerFromGraph(graphInvCrossSectionPi0Comb2760GeV,columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	TLine * lineFitPi02760GeVNLOALLEnergies =CreateLineFromHisto(histoFitInvCrossSectionPi02760GeV, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO);
	lineFitPi02760GeVNLOALLEnergies->Draw("same");
	TLine * lineNLOPi02760GeVMuHalfALLEnergies = CreateLineFromGraph(graphNLOMuHalfPi02760GeV,  columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[4] + offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
	TLine * lineNLOPi02760GeVMuOneALLEnergies = CreateLineFromGraph(graphNLOMuOnePi02760GeV, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
	TLine * lineNLOPi02760GeVMuTwoALLEnergies = CreateLineFromGraph(graphNLOMuTwoPi02760GeV, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[6]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[6]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
	
	
	//**************** forth Column **********************************************************
	TLatex *textPi0900GeVNLOALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[3],rowsNLOLegendPrelAndFinal[0],"#pi^{0}, #sqrt{#it{s}} = 0.9 TeV (*)");
	SetStyleTLatex( textPi0900GeVNLOALLEnergies, textSizeTopLablesPrelAndFinal,4);
	textPi0900GeVNLOALLEnergies->Draw();
	TLatex *textPi0900GeVNLOsysALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[3],rowsNLOLegendPrelAndFinal[1],"syst. + stat.");
	SetStyleTLatex( textPi0900GeVNLOsysALLEnergies, textSizeTopLowerLablesPrelAndFinal,4);
	textPi0900GeVNLOsysALLEnergies->Draw();
	TBox* boxCombinedPi0900GeVALLEnergies = CreateBoxFromGraph(graphInvCrossSectionPi0Comb900GeV,columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi0900GeVALLEnergies->Draw("l");
	TMarker* markerCombinedPi0900GeVALLEnergies = CreateMarkerFromGraph(graphInvCrossSectionPi0Comb900GeV ,columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	markerCombinedPi0900GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	TLine * lineFitPi0900GeVNLOALLEnergies = CreateLineFromHisto(histoFitInvCrossSectionPi0900GeV,columnsNLOLegendPrelAndFinal[3]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[3]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO);
	lineFitPi0900GeVNLOALLEnergies->Draw("same");
	TLine * lineNLOPi0900GeVMuHalfALLEnergies =  CreateLineFromGraph(graphNLOMuHalfPi0900GeV, columnsNLOLegendPrelAndFinal[3]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[4] + offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[3]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi0900GeVMuHalfALLEnergies->Draw("same");
	TLine * lineNLOPi0900GeVMuOneALLEnergies =  CreateLineFromGraph(graphNLOMuOnePi0900GeV,columnsNLOLegendPrelAndFinal[3]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[3]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi0900GeVMuOneALLEnergies->Draw("same");
	TLine * lineNLOPi0900GeVMuTwoALLEnergies =  CreateLineFromGraph(graphNLOMuTwoPi0900GeV, columnsNLOLegendPrelAndFinal[3]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[6]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[3]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[6]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi0900GeVMuTwoALLEnergies->Draw("same");
	TLine * lineNLOPi0900GeVMuTwoBKKALLEnergies = CreateLineFromGraph(graphRatioCombNLOBKKPi0900GeVMuTwo, columnsNLOLegendPrelAndFinal[3]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[7]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[3]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[7]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	lineNLOPi0900GeVMuTwoBKKALLEnergies->Draw("same");
	TLatex *textArxivALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[1],0.015,"(*) PLB 717 (2012) 162-172");
	SetStyleTLatex( textArxivALLEnergies, textSizeTopLablesPrelAndFinal*0.85,4);
	textArxivALLEnergies->Draw();
	TLatex *textArxivALLEnergies2 = new TLatex(columnsNLOLegendPrelAndFinal[2]+0.1,0.015,"(**) arXiv:1405.3794");
	SetStyleTLatex( textArxivALLEnergies2, textSizeTopLablesPrelAndFinal*0.85,4);
	textArxivALLEnergies2->Draw();

	
	//************************************************* End Legend ***************************************************
	
	padXSectionALLEnergiesRatioPi07TeV->cd();
	padXSectionALLEnergiesRatioPi07TeV->SetLogx();
	
	TH2F * ratio2DInvXSectionALLEnergiesPi0 = new TH2F("ratio2DInvXSectionALLEnergiesPi0","ratio2DInvXSectionALLEnergiesPi0",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionALLEnergiesPi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionALLEnergiesPi0->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFit->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFit->Draw("p,E2same");
	
	DrawGammaNLOTGraph( graphRatioCombNLOPi07TeVMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLOPi07TeVMuHalf);
	graphRatioCombNLOPi07TeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOPi07TeVMuOne, widthCommonFit, styleLineNLOMuOne, colorNLOPi07TeVMuOne);
	graphRatioCombNLOPi07TeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOPi07TeVMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLOPi07TeVMuTwo);
	graphRatioCombNLOPi07TeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi07TeVMuTwo->Draw("same,c");
// 	graphRatioCombNLODSSPi07TeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi07TeVRatio->Draw();
	labelRatioNLOPi07TeV->Draw();
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionALLEnergiesRatioPi02760GeV->cd();
	padXSectionALLEnergiesRatioPi02760GeV->SetLogx();
	
	TH2F * ratio2DInvXSectionALLEnergiesPi02760GeV;
	ratio2DInvXSectionALLEnergiesPi02760GeV = new TH2F("ratio2DInvXSectionALLEnergiesPi02760GeV","ratio2DInvXSectionALLEnergiesPi02760GeV",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionALLEnergiesPi02760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionALLEnergiesPi02760GeV->DrawCopy(); 
	
	TBox* boxErrorSigmaPi02760GeVRatio = CreateBoxConv(colorCommonSpectrumPi02760GeVBox, 0.25, 1.-(xSection2760GeVErr/(xSection2760GeV*1e3) ), 0.29, 1.+(xSection2760GeVErr/(xSection2760GeV*1e3)));
	boxErrorSigmaPi02760GeVRatio->Draw();
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit2760GeV, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFit2760GeV->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFit2760GeV->Draw("p,E2,same");
	
	DrawGammaNLOTGraph( graphRatioCombNLOPi02760GeVMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	graphRatioCombNLOPi02760GeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOPi02760GeVMuOne, widthCommonFit, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	graphRatioCombNLOPi02760GeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOPi02760GeVMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);
	graphRatioCombNLOPi02760GeVMuTwo->Draw("same,c");
	
	TLatex *labelRatioNLOPi02760GeV = new TLatex(0.18,0.75,"#pi^{0}, #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( labelRatioNLOPi02760GeV, 0.17,4);
	labelRatioNLOPi02760GeV->Draw();
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionALLEnergiesRatioPi0900GeV->cd();
	padXSectionALLEnergiesRatioPi0900GeV->SetLogx();
	
	TH2F * ratio2DInvXSectionALLEnergiesPi0900GeV=  new TH2F("ratio2DInvXSectionALLEnergiesPi0900GeV","ratio2DInvXSectionALLEnergiesPi0900GeV",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionALLEnergiesPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.12,0.13, 0.125,0.128, 0.9,0.45, 510, 505);
	ratio2DInvXSectionALLEnergiesPi0900GeV->GetXaxis()->SetLabelOffset(-0.03);
	ratio2DInvXSectionALLEnergiesPi0900GeV->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit900GeV, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFit900GeV->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFit900GeV->Draw("p,E2same");
	
	DrawGammaNLOTGraph( graphRatioCombNLOPi0900GeVMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLOPi0900GeVMuHalf);
	graphRatioCombNLOPi0900GeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOPi0900GeVMuOne, widthCommonFit, styleLineNLOMuOne, colorNLOPi0900GeVMuOne);
	graphRatioCombNLOPi0900GeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOPi0900GeVMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLOPi0900GeVMuTwo);
	graphRatioCombNLOPi0900GeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi0900GeVMuTwo->Draw("same,c");

	boxErrorSigmaPi0900GeVRatio->Draw();
	TLatex *labelRatioNLOPi0900GeVALLEnergies = new TLatex(0.18,0.85,"#pi^{0}, #sqrt{#it{s}} = 0.9 TeV");
	SetStyleTLatex( labelRatioNLOPi0900GeVALLEnergies, 0.115,4);
	labelRatioNLOPi0900GeVALLEnergies->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	canvasInvXSectionALLEnergies->Update();
	canvasInvXSectionALLEnergies->Print(Form("%s/Pi0_InvXSectionALLEnergies_Paper.%s",outputDir.Data(),suffix.Data()));
	padComparisonXSectionALLEnergies->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSFinal(pictDrawingCoordinatesCombineMeasCross[0]-0.05, pictDrawingCoordinatesCombineMeasCross[1], 0.032,collisionSystemCombinedReallyAll, kTRUE);
	
	canvasInvXSectionALLEnergies->Update();
	canvasInvXSectionALLEnergies->Print(Form("%s/Pi0_InvXSectionALLEnergies.%s",outputDir.Data(),suffix.Data()));

	
	//*****************************************************************************************************************
	//************************ Pi0 Meson Combined + NLO just spectrum******************************************************
	//*****************************************************************************************************************

	TCanvas* canvasInvXSectionNLOOnlySpectraPi0 = new TCanvas("canvasInvXSectionNLOOnlySpectraPi0","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlySpectraPi0,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionNLOOnlySpectraPi0 = new TPad("padComparisonXSectionNLOOnlySpectraPi0", "", 0., 0., 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionNLOOnlySpectraPi0, 0.15, 0.02, 0.02, 0.09);
	padComparisonXSectionNLOOnlySpectraPi0->Draw();
	
	TPad* padXSectionNLOLegendOnlySpectraPi0 = new TPad("padXSectionNLOLegendOnlySpectraPi0", "", 0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraPi0, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlySpectraPi0->Draw();
	
	padComparisonXSectionNLOOnlySpectraPi0->cd();
	padComparisonXSectionNLOOnlySpectraPi0->SetLogy();		
	padComparisonXSectionNLOOnlySpectraPi0->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionNLOOnlySpectraPi0 = new TH2F("histo2DInvXSectionNLOOnlySpectraPi0","histo2DInvXSectionNLOOnlySpectraPi0",1000,0.23,30.,1000,2e0,10e12);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLOOnlySpectraPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionNLOOnlySpectraPi0->DrawCopy(); 
	
	graphInvCrossSectionPi0Comb900GeV->Draw("p,E2same");
	graphInvCrossSectionPi0Comb7TeV->Draw("p,E2same");	
	graphInvCrossSectionPi0Comb2760GeV->Draw("p,E2same");
	
	fitInvCrossSectionPi0->Draw("same");
	histoFitInvCrossSectionPi0900GeV->Draw("same,c");
	histoFitInvCrossSectionPi02760GeV ->Draw("same,c");
	
	graphNLOMuHalfPi07TeV->Draw("same,c");
	graphNLOMuOnePi07TeV->Draw("same,c");
	graphNLOMuTwoPi07TeV->Draw("same,c");
	
	graphNLOMuHalfPi0900GeV->Draw("same,c");
	graphNLOMuOnePi0900GeV->Draw("same,c");
	graphNLOMuTwoPi0900GeV->Draw("same,c");
	
	graphNLOMuHalfPi02760GeV->Draw("same,c");
	graphNLOMuOnePi02760GeV->Draw("same,c");
	graphNLOMuTwoPi02760GeV->Draw("same,c");
	
	labelScalingPi07TeV->Draw();
	labelScalingPi02760GeV->Draw();
	labelScalingPi0900GeVALLEnergies->Draw();
	
	DrawNormalizationErrorText(normalizationInvXOnlySpec[0],normalizationInvXOnlySpec[1],normalizationInvXOnlySpec[2],
						  normalizationInvXOnlySpec[3],normalizationInvXOnlySpec[4],"all"); 
	
	//************************************************* Begin Legend ***************************************************
	
	padXSectionNLOLegendOnlySpectraPi0->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraPi0, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlySpectraPi0->SetBorderMode(-1);
	padXSectionNLOLegendOnlySpectraPi0->SetBorderSize(3);
	padXSectionNLOLegendOnlySpectraPi0->Draw();
	padXSectionNLOLegendOnlySpectraPi0->cd();

	//*************** first Column **********************************************************	
	textSpectrumALLEnergies->Draw();
	textFitCombALLEnergies->Draw();
	textNLOMuHalfALLEnergies->Draw();
	textNLOMuOneALLEnergies->Draw();
	textNLOMuTwoALLEnergies->Draw();
	//textNLOMuTwoBKKALLEnergies->Draw();
	
	//*************** second Column **********************************************************
	textPi07TeVNLOALLEnergies->Draw();
	textPi07TeVNLOsysALLEnergies->Draw();
	boxCombinedPi07TeVALLEnergies->Draw("l");
	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFit7TeVNLOALLEnergies->Draw("same");
	lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
	lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
	//lineNLOPi07TeVMuTwoBKKALLEnergies->Draw("same");
	
	//*************** third Column **********************************************************
	textPi02760GeVNLOALLEnergies->Draw();
	textPi02760GeVNLOsysALLEnergies->Draw();
	boxCombinedPi02760GeVALLEnergies->Draw("l");
	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi02760GeVNLOALLEnergies->Draw("same");
	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
		
	//**************** forth Column **********************************************************
	textPi0900GeVNLOALLEnergies->Draw();
	textPi0900GeVNLOsysALLEnergies->Draw();
	boxCombinedPi0900GeVALLEnergies->Draw("l");
	markerCombinedPi0900GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi0900GeVNLOALLEnergies->Draw("same");
	lineNLOPi0900GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi0900GeVMuOneALLEnergies->Draw("same");
	lineNLOPi0900GeVMuTwoALLEnergies->Draw("same");
	//lineNLOPi0900GeVMuTwoBKKALLEnergies->Draw("same");
	textArxivALLEnergies->Draw();
	canvasInvXSectionNLOOnlySpectraPi0->Update();
	canvasInvXSectionNLOOnlySpectraPi0->Print(Form("%s/Pi0_InvXSectionNLO_OnlySpectrum_Paper.%s",outputDir.Data(),suffix.Data()));
	
	padComparisonXSectionNLOOnlySpectraPi0->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSPrelim(pictDrawingCoordinatesOnlySpectrum[0]-0.05, pictDrawingCoordinatesOnlySpectrum[1], pictDrawingCoordinatesOnlySpectrum[2], pictDrawingCoordinatesOnlySpectrum[3], pictDrawingCoordinatesOnlySpectrum[4], pictDrawingCoordinatesOnlySpectrum[5], pictDrawingCoordinatesOnlySpectrum[6], pictDrawingCoordinatesOnlySpectrum[7], pictDrawingCoordinates[8],collisionSystemCombinedReallyAll, pictDrawingOptions[1], pictDrawingOptions[2], kTRUE,1200,960);
	canvasInvXSectionNLOOnlySpectraPi0->Update();
	canvasInvXSectionNLOOnlySpectraPi0->Print(Form("%s/Pi0_InvXSectionNLO_OnlySpectrum.%s",outputDir.Data(),suffix.Data()));
	
//**************************************************************************************************************************************
//************************************ Inv Crosssection + NLO only ratio Pi0 ***************************************************************
//**************************************************************************************************************************************
	
	TCanvas* canvasInvXSectionNLOOnlyRatioPi0 = new TCanvas("canvasInvXSectionNLOOnlyRatioPi0","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlyRatioPi0,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padXSectionNLOLegendOnlyRatioPi0 = new TPad("padXSectionNLOLegendOnlyRatioPi0", "",  0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioPi0, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlyRatioPi0->Draw();

	TPad* padXSectionNLOOnlyRatioPi0Pi07TeV = new TPad("padXSectionNLOOnlyRatioPi0Pi07TeV", "", 0., 0.565, 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioPi0Pi07TeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOOnlyRatioPi0Pi07TeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioPi0Pi02760GeV = new TPad("padXSectionNLOOnlyRatioPi0Pi02760GeV", "", 0., 0.33, 1., 0.565,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioPi0Pi02760GeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOOnlyRatioPi0Pi02760GeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioPi0Pi0900GeV = new TPad("padXSectionNLOOnlyRatioPi0Pi0900GeV", "", 0., 0., 1., 0.33,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioPi0Pi0900GeV, 0.15, 0.02, 0., 0.28);
	padXSectionNLOOnlyRatioPi0Pi0900GeV->Draw();
	
	padXSectionNLOOnlyRatioPi0Pi07TeV->cd();
	padXSectionNLOOnlyRatioPi0Pi07TeV->SetLogx();
	padXSectionNLOOnlyRatioPi0Pi07TeV->SetLogy();

	TH2F * ratio2DInvXSectionNLOPi0 = new TH2F("ratio2DInvXSectionNLOPi0","ratio2DInvXSectionNLOPi0",1000,0.23,30.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOPi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionNLOPi0->DrawCopy(); 
	graphRatioCombNLOPi07TeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi07TeVMuOne->Draw("same,c");
	graphRatioCombNLOPi07TeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi07TeVMuTwo->Draw("same,c");
	boxErrorSigmaPi07TeVRatio->Draw();
	
	graphRatioCombCombFit->Draw("p,E2same");
	labelRatioNLOPi07TeV->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOOnlyRatioPi0Pi02760GeV->cd();
	padXSectionNLOOnlyRatioPi0Pi02760GeV->SetLogx();
	padXSectionNLOOnlyRatioPi0Pi02760GeV->SetLogy();
	
	TH2F * ratio2DInvXSectionNLOPi0900GeV;
	ratio2DInvXSectionNLOPi0900GeV = new TH2F("ratio2DInvXSectionNLOPi0900GeV","ratio2DInvXSectionNLOPi0900GeV",1000,0.23,30.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionNLOPi0900GeV->DrawCopy(); 
	
	graphRatioCombNLOPi02760GeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi02760GeVMuOne->Draw("same,c");
	graphRatioCombNLOPi02760GeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi02760GeVRatio->Draw();
	graphRatioCombCombFit2760GeV->Draw("p,E2same");
	labelRatioNLOPi02760GeV->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOOnlyRatioPi0Pi0900GeV->cd();
	padXSectionNLOOnlyRatioPi0Pi0900GeV->SetLogx();
	padXSectionNLOOnlyRatioPi0Pi0900GeV->SetLogy();
	TH2F * ratio2DInvXSectionNLOEta=  new TH2F("ratio2DInvXSectionNLOEta","ratio2DInvXSectionNLOEta",1000,0.23,30.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.12,0.13, 0.125,0.128, 0.9,0.45, 510, 505);
	ratio2DInvXSectionNLOEta->GetXaxis()->SetLabelOffset(-0.03);
	ratio2DInvXSectionNLOEta->DrawCopy(); 
	graphRatioCombNLOBKKPi0900GeVMuTwo->Draw("same,c");
	
	graphRatioCombNLOPi0900GeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuOne->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuTwo->Draw("same,c");
	
	graphRatioCombCombFit900GeV->Draw("p,E2same");
	boxErrorSigmaPi0900GeVRatio->Draw();
	labelRatioNLOPi0900GeVALLEnergies->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
		
	padXSectionNLOLegendOnlyRatioPi0->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioPi0, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlyRatioPi0->SetBorderMode(-1);
	padXSectionNLOLegendOnlyRatioPi0->SetBorderSize(3);
	padXSectionNLOLegendOnlyRatioPi0->Draw();
	padXSectionNLOLegendOnlyRatioPi0->cd();
	
	//****************************** Definition of the Legend ******************************************
	//**************** Row def ************************
	Double_t rowsNLOLegendPrelAndFinalOnlyRatioPi0[7] = {0.9,0.79,0.69,0.54,0.41,0.26,0.12};
	Double_t columnsNLOLegendPrelAndFinalOnlyRatio[4] = {0.,0.23,0.46,0.74};
	Double_t textSizeLeftLabelsPrelAndFinalRatio = 0.105;
	Double_t textSizeTopLablesPrelAndFinalRatio = 0.11;
 	Double_t textSizeTopLowerLablesPrelAndFinalRatio = 0.105;
		
	//*************** first Column **********************************************************
	TLatex *textSpectrumOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[2],"combined Spec.");
	SetStyleTLatex( textSpectrumOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	textSpectrumOnlyRatioPi0->Draw();
	TLatex *textNLOMuHalfOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[3],"NLO #mu= 0.5 #it{p}_{T}");
	SetStyleTLatex( textNLOMuHalfOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	textNLOMuHalfOnlyRatioPi0->Draw();
	TLatex *textNLOMuOneOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[4],"NLO #mu= #it{p}_{T}");
	SetStyleTLatex( textNLOMuOneOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	textNLOMuOneOnlyRatioPi0->Draw();
	TLatex *textNLOMuTwoOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[5],"NLO #mu= 2 #it{p}_{T}");
	SetStyleTLatex( textNLOMuTwoOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	textNLOMuTwoOnlyRatioPi0->Draw();
	TLatex *textNLOBKKMuTwoOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[6],"NLO #mu= 2 #it{p}_{T} (BKK)");
	SetStyleTLatex( textNLOBKKMuTwoOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	textNLOBKKMuTwoOnlyRatioPi0->Draw();
	
	
	//*************** second Column **********************************************************
	TLatex *textPi07TeVNLOOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[1],rowsNLOLegendPrelAndFinalOnlyRatioPi0[0],"#pi^{0}, #sqrt{#it{s}} = 7 TeV (*)");
	SetStyleTLatex( textPi07TeVNLOOnlyRatioPi0, textSizeTopLablesPrelAndFinalRatio,4);
	textPi07TeVNLOOnlyRatioPi0->Draw();
	TLatex *textPi07TeVNLOsysOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[1],rowsNLOLegendPrelAndFinalOnlyRatioPi0[1],"syst. + stat.");
	SetStyleTLatex( textPi07TeVNLOsysOnlyRatioPi0, textSizeTopLowerLablesPrelAndFinalRatio,4);
	textPi07TeVNLOsysOnlyRatioPi0->Draw();
	TBox* boxCombinedPi07TeVOnlyRatioPi0 = CreateBoxFromGraph(graphInvCrossSectionPi0Comb7TeV,columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi07TeVOnlyRatioPi0->Draw("l");
	TMarker* markerCombinedPi07TeVOnlyRatioPi0 = CreateMarkerFromGraph(graphInvCrossSectionPi0Comb7TeV,columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	markerCombinedPi07TeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	TLine * lineNLOPi07TeVMuHalfOnlyRatioPi0 = CreateLineFromGraph(graphNLOMuHalfPi07TeV, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3] + offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOPi07TeVMuHalfOnlyRatioPi0->Draw("same");
	TLine * lineNLOPi07TeVMuOneOnlyRatioPi0 = CreateLineFromGraph(graphNLOMuOnePi07TeV,columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOPi07TeVMuOneOnlyRatioPi0->Draw("same");
	TLine * lineNLOPi07TeVMuTwoOnlyRatioPi0 = CreateLineFromGraph(graphNLOMuTwoPi07TeV, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOPi07TeVMuTwoOnlyRatioPi0->Draw("same");
	TLine * lineNLOBKKPi07TeVMuTwoOnlyRatioPi0 = CreateLineFromGraph(graphRatioCombNLOBKKPi07TeVMuTwo, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[6]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[6]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOBKKPi07TeVMuTwoOnlyRatioPi0->Draw("same");
	
	
	//*************** third Column **********************************************************
	TLatex *textPi02760GeVNLOOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[2],rowsNLOLegendPrelAndFinalOnlyRatioPi0[0],"#pi^{0}, #sqrt{#it{s}} = 2.76 TeV (**)");
	SetStyleTLatex( textPi02760GeVNLOOnlyRatioPi0, textSizeTopLablesPrelAndFinalRatio,4);
	textPi02760GeVNLOOnlyRatioPi0->Draw();
	TLatex *textPi02760GeVNLOsysOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[2],rowsNLOLegendPrelAndFinalOnlyRatioPi0[1],"syst. + stat.");
	SetStyleTLatex( textPi02760GeVNLOsysOnlyRatioPi0, textSizeTopLowerLablesPrelAndFinalRatio,4);
	textPi02760GeVNLOsysOnlyRatioPi0->Draw();
	TBox* boxCombinedPi02760GeVOnlyRatioPi0 = CreateBoxFromGraph(graphInvCrossSectionPi0Comb2760GeV,columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi02760GeVOnlyRatioPi0->Draw("l");
	TMarker* markerCombinedPi02760GeVOnlyRatioPi0 = CreateMarkerFromGraph(graphInvCrossSectionPi0Comb2760GeV,columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	markerCombinedPi02760GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	TLine * lineNLOPi02760GeVMuHalfOnlyRatioPi0 = CreateLineFromGraph(graphNLOMuHalfPi02760GeV,  columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3] + offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOPi02760GeVMuHalfOnlyRatioPi0->Draw("same");
	TLine * lineNLOPi02760GeVMuOneOnlyRatioPi0 = CreateLineFromGraph(graphNLOMuOnePi02760GeV, columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOPi02760GeVMuOneOnlyRatioPi0->Draw("same");
	TLine * lineNLOPi02760GeVMuTwoOnlyRatioPi0 = CreateLineFromGraph(graphNLOMuTwoPi02760GeV, columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOPi02760GeVMuTwoOnlyRatioPi0->Draw("same");
	
	//**************** forth Column **********************************************************
	TLatex *textPi0900GeVNLOOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[3],rowsNLOLegendPrelAndFinalOnlyRatioPi0[0],"#pi^{0}, #sqrt{#it{s}} = 0.9 TeV (*)");
	SetStyleTLatex( textPi0900GeVNLOOnlyRatioPi0, textSizeTopLablesPrelAndFinalRatio,4);
	textPi0900GeVNLOOnlyRatioPi0->Draw();
	TLatex *textPi0900GeVNLOsysOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[3],rowsNLOLegendPrelAndFinalOnlyRatioPi0[1],"syst. + stat.");
	SetStyleTLatex( textPi0900GeVNLOsysOnlyRatioPi0, textSizeTopLowerLablesPrelAndFinalRatio,4);
	textPi0900GeVNLOsysOnlyRatioPi0->Draw();
	TBox* boxCombinedPi0900GeVOnlyRatioPi0 = CreateBoxFromGraph(graphInvCrossSectionPi0Comb900GeV,columnsNLOLegendPrelAndFinalOnlyRatio[3]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[3]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi0900GeVOnlyRatioPi0->Draw("l");
	TMarker* markerCombinedPi0900GeVOnlyRatioPi0 = CreateMarkerFromGraph(graphInvCrossSectionPi0Comb900GeV,columnsNLOLegendPrelAndFinalOnlyRatio[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	markerCombinedPi0900GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	TLine * lineNLOPi0900GeVMuHalfOnlyRatioPi0 =  CreateLineFromGraph(graphNLOMuHalfPi0900GeV, columnsNLOLegendPrelAndFinalOnlyRatio[3]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3] + offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[3]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOPi0900GeVMuHalfOnlyRatioPi0->Draw("same");
	TLine * lineNLOPi0900GeVMuOneOnlyRatioPi0 =  CreateLineFromGraph(graphNLOMuOnePi0900GeV,columnsNLOLegendPrelAndFinalOnlyRatio[3]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[3]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOPi0900GeVMuOneOnlyRatioPi0->Draw("same");
	TLine * lineNLOPi0900GeVMuTwoOnlyRatioPi0 =  CreateLineFromGraph(graphNLOMuTwoPi0900GeV, columnsNLOLegendPrelAndFinalOnlyRatio[3]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[3]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOPi0900GeVMuTwoOnlyRatioPi0->Draw("same");
	TLine * lineNLOBKKPi0900GeVMuTwoOnlyRatioPi0 = CreateLineFromGraph(graphRatioCombNLOBKKPi0900GeVMuTwo, columnsNLOLegendPrelAndFinalOnlyRatio[3]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[6]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[3]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[6]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	lineNLOBKKPi0900GeVMuTwoOnlyRatioPi0->Draw("same");
	TLatex *textArxivALLEnergiesRatio = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[1],0.02,"(*) PLB 717 (2012) 162-172");
	SetStyleTLatex( textArxivALLEnergiesRatio, textSizeTopLablesPrelAndFinalRatio*0.85,4);
	textArxivALLEnergiesRatio->Draw();
	TLatex *textArxivALLEnergiesRatio2 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[2]+0.1,0.02,"(**) arXiv:1405.3794");
	SetStyleTLatex( textArxivALLEnergiesRatio2, textSizeTopLablesPrelAndFinalRatio*0.85,4);
	textArxivALLEnergiesRatio2->Draw();

	//************************************************* End Legend ***************************************************
	
	
	canvasInvXSectionNLOOnlyRatioPi0->Update();
	canvasInvXSectionNLOOnlyRatioPi0->Print(Form("%s/Pi0_InvXSectionNLO_OnlyRatio_Paper.%s",outputDir.Data(),suffix.Data()));
	
	
	//***************************************************************************************************************
	//************************ Pi0 Combined + NLO separated systematic & statistical errors *************************
	//***************************************************************************************************************
	
	TCanvas* canvasInvXSectionALLEnergiesSep = new TCanvas("canvasInvXSectionALLEnergiesSep","",200,10,1200,2000);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionALLEnergiesSep,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionALLEnergiesSep = new TPad("padComparisonXSectionALLEnergiesSep", "", 0., 0.42, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionALLEnergiesSep, 0.15, 0.02, 0.02, 0.);
	padComparisonXSectionALLEnergiesSep->Draw();
	
	TPad* padXSectionALLEnergiesSepRatioPi07TeV = new TPad("padXSectionALLEnergiesSepRatioPi07TeV", "", 0., 0.30, 1., 0.42,-1, -1, -2);
	DrawGammaPadSettings( padXSectionALLEnergiesSepRatioPi07TeV, 0.15, 0.02, 0., 0.);
	padXSectionALLEnergiesSepRatioPi07TeV->Draw();
	
	TPad* padXSectionALLEnergiesSepRatioPi02760GeV = new TPad("padXSectionALLEnergiesSepRatioPi02760GeV", "", 0., 0.18, 1., 0.30,-1, -1, -2);
	DrawGammaPadSettings( padXSectionALLEnergiesSepRatioPi02760GeV, 0.15, 0.02, 0., 0.);
	padXSectionALLEnergiesSepRatioPi02760GeV->Draw();
	
	TPad* padXSectionALLEnergiesSepRatioPi0900GeV = new TPad("padXSectionALLEnergiesSepRatioPi0900GeV", "", 0., 0., 1., 0.18,-1, -1, -2);
	DrawGammaPadSettings( padXSectionALLEnergiesSepRatioPi0900GeV, 0.15, 0.02, 0., 0.28);
	padXSectionALLEnergiesSepRatioPi0900GeV->Draw();
	
	
	padComparisonXSectionALLEnergiesSep->cd();
	padComparisonXSectionALLEnergiesSep->SetLogy();		
	padComparisonXSectionALLEnergiesSep->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionALLEnergiesSep = new TH2F("histo2DInvXSectionALLEnergiesSep","histo2DInvXSectionALLEnergiesSep",1000,0.23,30.,1000,2e-3,10e11);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionALLEnergiesSep, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionALLEnergiesSep->GetXaxis()->SetTickLength(0.02);
	histo2DInvXSectionALLEnergiesSep->GetYaxis()->SetTickLength(0.02);
	histo2DInvXSectionALLEnergiesSep->DrawCopy(); 
	
	graphInvCrossSectionPi0Comb900GeVStatErr = ScaleGraph(graphInvCrossSectionPi0Comb900GeVStatErr,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb900GeVStatErr, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionPi0Comb900GeVStatErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionPi0Comb900GeVSysErr = ScaleGraph(graphInvCrossSectionPi0Comb900GeVSysErr,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb900GeVSysErr, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE,colorCommonSpectrumPi0900GeVBox);
	graphInvCrossSectionPi0Comb900GeVSysErr->SetLineWidth(0);
	graphInvCrossSectionPi0Comb900GeVSysErr->Draw("2,same");
	graphInvCrossSectionPi0Comb900GeVStatErr->Draw("p,same");
	
	graphInvCrossSectionPi0Comb2760GeVStatErr = ScaleGraph(graphInvCrossSectionPi0Comb2760GeVStatErr,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb2760GeVStatErr, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionPi0Comb2760GeVStatErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionPi0Comb2760GeVSysErr = ScaleGraph(graphInvCrossSectionPi0Comb2760GeVSysErr,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb2760GeVSysErr, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi02760GeVBox);
	graphInvCrossSectionPi0Comb2760GeVSysErr->SetLineWidth(0);
	graphInvCrossSectionPi0Comb2760GeVSysErr->Draw("2,same");
	graphInvCrossSectionPi0Comb2760GeVStatErr->Draw("p,same");
	
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb7TeVStatErr, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionPi0Comb7TeVStatErr->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb7TeVSysErr, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi07TeVBox);
	graphInvCrossSectionPi0Comb7TeVSysErr->SetLineWidth(0);
	graphInvCrossSectionPi0Comb7TeVSysErr->Draw("2,same");
	graphInvCrossSectionPi0Comb7TeVStatErr->Draw("p,same");
		
	fitInvCrossSectionPi0->Draw("same");
	histoFitInvCrossSectionPi0900GeV->Draw("same,c");
	histoFitInvCrossSectionPi02760GeV->Draw("same,c");

	graphNLOMuHalfPi07TeV->Draw("same,c");
	graphNLOMuOnePi07TeV->Draw("same,c");
	graphNLOMuTwoPi07TeV->Draw("same,c");
	
	graphNLOMuHalfPi02760GeV->Draw("same,c");
	graphNLOMuOnePi02760GeV->Draw("same,c");
	graphNLOMuTwoPi02760GeV->Draw("same,c");
	
	
	graphNLOMuHalfPi0900GeV->Draw("same,c");
	graphNLOMuOnePi0900GeV->Draw("same,c");
	graphNLOMuTwoPi0900GeV->Draw("same,c");
	
	labelScalingPi07TeV->Draw();
	labelScalingPi02760GeV->Draw();
	labelScalingPi0900GeVALLEnergies->Draw();	
	
	DrawNormalizationErrorText(normalizationInvX1MesonALLEn[0],normalizationInvX1MesonALLEn[1],normalizationInvX1MesonALLEn[2],
						  normalizationInvX1MesonALLEn[3],normalizationInvX1MesonALLEn[4],"all"); 

	TPad* padXSectionALLEnergiesLegendSep = new TPad("padXSectionALLEnergiesLegendSep", "", 0.18, 0.02, 0.98, 0.22,-1, -1, -2); 
	DrawGammaPadSettings( padXSectionALLEnergiesLegendSep, 0., 0., 0., 0.);
	padXSectionALLEnergiesLegendSep->SetFillStyle(0);	
	padXSectionALLEnergiesLegendSep->Draw();	
	padXSectionALLEnergiesLegendSep->cd();
	//*************** first Column **********************************************************
	textSpectrumALLEnergies->Draw();
	textFitCombALLEnergies->Draw();
	textNLOMuHalfALLEnergies->Draw();
	textNLOMuOneALLEnergies->Draw();
	textNLOMuTwoALLEnergies->Draw();
	textNLOMuTwoBKKALLEnergies->Draw();
	
	//*************** second Column **********************************************************
	textPi07TeVNLOALLEnergies->Draw();
	TLatex *textPi07TeVNLOsysALLEnergiesSep = new TLatex(columnsNLOLegendPrelAndFinal[1],rowsNLOLegendPrelAndFinal[1],"syst., stat.");
	SetStyleTLatex( textPi07TeVNLOsysALLEnergiesSep, textSizeTopLowerLablesPrelAndFinal,4);
	textPi07TeVNLOsysALLEnergiesSep->Draw();
	TBox* boxCombinedPi07TeVALLEnergiesSep = CreateBoxFromGraphWithFill(graphInvCrossSectionPi0Comb7TeVSysErr,columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi07TeVALLEnergiesSep->Draw("f");
	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFit7TeVNLOALLEnergies->Draw("same");
	lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
	lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
	lineNLOPi07TeVMuTwoBKKALLEnergies->Draw("same");
	
	//*************** third Column **********************************************************
	textPi02760GeVNLOALLEnergies->Draw();
	TLatex *textPi02760GeVNLOsysALLEnergiesSep = new TLatex(columnsNLOLegendPrelAndFinal[2],rowsNLOLegendPrelAndFinal[1],"syst., stat.");
	SetStyleTLatex( textPi02760GeVNLOsysALLEnergiesSep, textSizeTopLowerLablesPrelAndFinal,4);
	textPi02760GeVNLOsysALLEnergiesSep->Draw();
	TBox* boxCombinedPi02760GeVALLEnergiesSep = CreateBoxFromGraphWithFill(graphInvCrossSectionPi0Comb2760GeVSysErr,columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi02760GeVALLEnergiesSep->Draw("f");
	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi02760GeVNLOALLEnergies->Draw("same");
	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
	
	
	//**************** forth Column **********************************************************
	textPi0900GeVNLOALLEnergies->Draw();
	TLatex *textPi0900GeVNLOsysALLEnergiesSep = new TLatex(columnsNLOLegendPrelAndFinal[3],rowsNLOLegendPrelAndFinal[1],"syst., stat.");
	SetStyleTLatex( textPi0900GeVNLOsysALLEnergiesSep, textSizeTopLowerLablesPrelAndFinal,4);
	textPi0900GeVNLOsysALLEnergiesSep->Draw();
	TBox* boxCombinedPi0900GeVALLEnergiesSep = CreateBoxFromGraphWithFill(graphInvCrossSectionPi0Comb900GeVSysErr,columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi0900GeVALLEnergiesSep->Draw("f");
	markerCombinedPi0900GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi0900GeVNLOALLEnergies->Draw("same");
	lineNLOPi0900GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi0900GeVMuOneALLEnergies->Draw("same");
	lineNLOPi0900GeVMuTwoALLEnergies->Draw("same");
	lineNLOPi0900GeVMuTwoBKKALLEnergies->Draw("same");
	textArxivALLEnergies->Draw();
	textArxivALLEnergies2->Draw();
	
	padXSectionALLEnergiesSepRatioPi07TeV->cd();
	padXSectionALLEnergiesSepRatioPi07TeV->SetLogx();
	
	TH2F * ratio2DInvXSectionALLEnergiesSepPi0 = new TH2F("ratio2DInvXSectionALLEnergiesSepPi0","ratio2DInvXSectionALLEnergiesSepPi0",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionALLEnergiesSepPi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionALLEnergiesSepPi0->GetXaxis()->SetTickLength(0.07);
	ratio2DInvXSectionALLEnergiesSepPi0->GetYaxis()->SetTickLength(0.02);
	ratio2DInvXSectionALLEnergiesSepPi0->GetYaxis()->SetLabelOffset(0.01);
	ratio2DInvXSectionALLEnergiesSepPi0->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kFALSE);
	graphRatioCombCombFitStat->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi07TeVBox);
	graphRatioCombCombFitSys->SetLineWidth(0);
	graphRatioCombCombFitSys->Draw("2,same");
	graphRatioCombCombFit->Draw("p,same");
	
	graphRatioCombNLOPi07TeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi07TeVMuOne->Draw("same,c");
	graphRatioCombNLOPi07TeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi07TeVMuTwo->Draw("same,c");
// 	graphRatioCombNLODSSPi07TeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi07TeVRatio->Draw();
	labelRatioNLOPi07TeV->Draw();
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionALLEnergiesSepRatioPi02760GeV->cd();
	padXSectionALLEnergiesSepRatioPi02760GeV->SetLogx();
	TH2F * ratio2DInvXSectionALLEnergiesSepPi02760GeV;
	ratio2DInvXSectionALLEnergiesSepPi02760GeV = new TH2F("ratio2DInvXSectionALLEnergiesSepPi02760GeV","ratio2DInvXSectionALLEnergiesSepPi02760GeV",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionALLEnergiesSepPi02760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionALLEnergiesSepPi02760GeV->GetXaxis()->SetTickLength(0.07);
	ratio2DInvXSectionALLEnergiesSepPi02760GeV->GetYaxis()->SetTickLength(0.02);
	ratio2DInvXSectionALLEnergiesSepPi02760GeV->GetYaxis()->SetLabelOffset(0.01);
	ratio2DInvXSectionALLEnergiesSepPi02760GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit2760GeVStat, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kFALSE);
	graphRatioCombCombFit2760GeVStat->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit2760GeVSys, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi02760GeVBox);
	graphRatioCombCombFit2760GeVSys->SetLineWidth(0);
	graphRatioCombCombFit2760GeVSys->Draw("2,same");
	graphRatioCombCombFit2760GeVStat->Draw("p,same");
	
	graphRatioCombNLOPi02760GeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi02760GeVMuOne->Draw("same,c");
	graphRatioCombNLOPi02760GeVMuTwo->Draw("same,c");
	
	labelRatioNLOPi02760GeV->Draw();
	boxErrorSigmaPi02760GeVRatio->Draw();
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionALLEnergiesSepRatioPi0900GeV->cd();
	padXSectionALLEnergiesSepRatioPi0900GeV->SetLogx();
	TH2F * ratio2DInvXSectionALLEnergiesSepPi0900GeV=  new TH2F("ratio2DInvXSectionALLEnergiesSepPi0900GeV","ratio2DInvXSectionALLEnergiesSepPi0900GeV",1000,0.23,30.,1000,0.3,3.62);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionALLEnergiesSepPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.12,0.13, 0.125,0.128, 0.9,0.45, 510, 505);
	ratio2DInvXSectionALLEnergiesSepPi0900GeV->GetXaxis()->SetLabelOffset(-0.03);
	ratio2DInvXSectionALLEnergiesSepPi0900GeV->GetYaxis()->SetLabelOffset(0.01);
	ratio2DInvXSectionALLEnergiesSepPi0900GeV->GetYaxis()->SetTickLength(0.03);
	ratio2DInvXSectionALLEnergiesSepPi0900GeV->GetXaxis()->SetTickLength(0.05);
// 	ratio2DInvXSectionALLEnergiesSepPi0900GeV->GetYaxis()->SetNdivisions(512);
	ratio2DInvXSectionALLEnergiesSepPi0900GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit900GeVStat, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kFALSE);
	graphRatioCombCombFit900GeVStat->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit900GeVSys, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi0900GeVBox);
	graphRatioCombCombFit900GeVSys->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFit900GeVSys->Draw("2same");
	graphRatioCombCombFit900GeVStat->Draw("psame");
	
	graphRatioCombNLOPi0900GeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuOne->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi0900GeVMuTwo->Draw("same,c");

	boxErrorSigmaPi0900GeVRatio->Draw();
	TLatex *labelRatioNLOPi0900GeVALLEnergiesSep = new TLatex(0.18,0.85,"#pi^{0}, #sqrt{#it{s}} = 0.9 TeV");
	SetStyleTLatex( labelRatioNLOPi0900GeVALLEnergiesSep, 0.115,4);
	labelRatioNLOPi0900GeVALLEnergiesSep->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	canvasInvXSectionALLEnergiesSep->Update();
	canvasInvXSectionALLEnergiesSep->Print(Form("%s/Pi0_InvXSectionALLEnergiesSep_Paper.%s",outputDir.Data(),suffix.Data()));
	padComparisonXSectionALLEnergiesSep->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSFinal(pictDrawingCoordinatesCombineMeasCross[0]-0.05, pictDrawingCoordinatesCombineMeasCross[1], 0.032, collisionSystemCombinedReallyAll, kTRUE);
	
	canvasInvXSectionALLEnergiesSep->Update();
	canvasInvXSectionALLEnergiesSep->Print(Form("%s/Pi0_InvXSectionALLEnergiesSep.%s",outputDir.Data(),suffix.Data()));

	
	//*****************************************************************************************************************
	//************************ Pi0 Meson Combined + NLO just spectrum + separated sys & stat errors *******************
	//*****************************************************************************************************************

	TCanvas* canvasInvXSectionNLOOnlySpectraPi0Sep = new TCanvas("canvasInvXSectionNLOOnlySpectraPi0Sep","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlySpectraPi0Sep,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionNLOOnlySpectraPi0Sep = new TPad("padComparisonXSectionNLOOnlySpectraPi0Sep", "", 0., 0., 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionNLOOnlySpectraPi0Sep, 0.15, 0.02, 0.02, 0.09);
	padComparisonXSectionNLOOnlySpectraPi0Sep->Draw();
	
	TPad* padXSectionNLOLegendOnlySpectraPi0Sep = new TPad("padXSectionNLOLegendOnlySpectraPi0Sep", "", 0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraPi0Sep, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlySpectraPi0Sep->Draw();

	padComparisonXSectionNLOOnlySpectraPi0Sep->cd();
	padComparisonXSectionNLOOnlySpectraPi0Sep->SetLogy();		
	padComparisonXSectionNLOOnlySpectraPi0Sep->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionNLOOnlySpectraPi0Sep = new TH2F("histo2DInvXSectionNLOOnlySpectraPi0Sep","histo2DInvXSectionNLOOnlySpectraPi0Sep",1000,0.23,30.,1000,2e0,10e11);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLOOnlySpectraPi0Sep, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionNLOOnlySpectraPi0Sep->DrawCopy(); 
	graphInvCrossSectionPi0Comb900GeVSysErr->Draw("2,same");
	graphInvCrossSectionPi0Comb900GeVStatErr->Draw("p,same");
	graphInvCrossSectionPi0Comb2760GeVSysErr->Draw("2,same");
	graphInvCrossSectionPi0Comb2760GeVStatErr->Draw("p,same");
	graphInvCrossSectionPi0Comb7TeVSysErr->Draw("2,same");
	graphInvCrossSectionPi0Comb7TeVStatErr->Draw("p,same");
		
	fitInvCrossSectionPi0->Draw("same");
	histoFitInvCrossSectionPi0900GeV->Draw("same,c");
	histoFitInvCrossSectionPi02760GeV->Draw("same,c");

	graphNLOMuHalfPi07TeV->Draw("same,c");
	graphNLOMuOnePi07TeV->Draw("same,c");
	graphNLOMuTwoPi07TeV->Draw("same,c");
	
	graphNLOMuHalfPi02760GeV->Draw("same,c");
	graphNLOMuOnePi02760GeV->Draw("same,c");
	graphNLOMuTwoPi02760GeV->Draw("same,c");
	
	graphNLOMuHalfPi0900GeV->Draw("same,c");
	graphNLOMuOnePi0900GeV->Draw("same,c");
	graphNLOMuTwoPi0900GeV->Draw("same,c");
	
	labelScalingPi07TeV->Draw();
	labelScalingPi02760GeV->Draw();
	labelScalingPi0900GeVALLEnergies->Draw();	
	
	DrawNormalizationErrorText(normalizationInvX1MesonALLEn[0],normalizationInvX1MesonALLEn[1],0.14,
						  0.04*1.1,0.04,"all"); 
	
	padXSectionNLOLegendOnlySpectraPi0Sep->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraPi0Sep, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlySpectraPi0Sep->SetBorderMode(-1);
	padXSectionNLOLegendOnlySpectraPi0Sep->SetBorderSize(3);
	padXSectionNLOLegendOnlySpectraPi0Sep->Draw();
	padXSectionNLOLegendOnlySpectraPi0Sep->cd();

	//*************** first Column **********************************************************	
	textSpectrumALLEnergies->Draw();
	textFitCombALLEnergies->Draw();
	textNLOMuHalfALLEnergies->Draw();
	textNLOMuOneALLEnergies->Draw();
	textNLOMuTwoALLEnergies->Draw();
	//textNLOMuTwoBKKALLEnergies->Draw();
	
	//*************** second Column **********************************************************
	textPi07TeVNLOALLEnergies->Draw();
	textPi07TeVNLOsysALLEnergiesSep->Draw();
	boxCombinedPi07TeVALLEnergiesSep->Draw("f");
	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFit7TeVNLOALLEnergies->Draw("same");
	lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
	lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
	//lineNLOPi07TeVMuTwoBKKALLEnergies->Draw("same");
	
	//*************** third Column **********************************************************
	textPi02760GeVNLOALLEnergies->Draw();
	textPi02760GeVNLOsysALLEnergiesSep->Draw();
	boxCombinedPi02760GeVALLEnergiesSep->Draw("f");
	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi02760GeVNLOALLEnergies->Draw("same");
	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
		
	//**************** forth Column **********************************************************
	textPi0900GeVNLOALLEnergies->Draw();
	textPi0900GeVNLOsysALLEnergiesSep->Draw();
	boxCombinedPi0900GeVALLEnergiesSep->Draw("f");
	markerCombinedPi0900GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi0900GeVNLOALLEnergies->Draw("same");
	lineNLOPi0900GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi0900GeVMuOneALLEnergies->Draw("same");
	lineNLOPi0900GeVMuTwoALLEnergies->Draw("same");
	//lineNLOPi0900GeVMuTwoBKKALLEnergies->Draw("same");
	textArxivALLEnergies->Draw();
	textArxivALLEnergies2->Draw();
	
	canvasInvXSectionNLOOnlySpectraPi0Sep->Update();
	canvasInvXSectionNLOOnlySpectraPi0Sep->Print(Form("%s/Pi0_InvXSectionNLO_OnlySpectrumSep_Paper.%s",outputDir.Data(),suffix.Data()));
	
	padComparisonXSectionNLOOnlySpectraPi0Sep->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSFinal(pictDrawingCoordinatesOnlySpectrum[0]-0.15, pictDrawingCoordinatesOnlySpectrum[1],0.04,collisionSystemCombinedReallyAll,kTRUE);
	
	canvasInvXSectionNLOOnlySpectraPi0Sep->Update();
	canvasInvXSectionNLOOnlySpectraPi0Sep->Print(Form("%s/Pi0_InvXSectionNLO_OnlySpectrumSep.%s",outputDir.Data(),suffix.Data()));
	
//**************************************************************************************************************************************
//************************************ Inv Crosssection + NLO only ratio Pi0 + separated sys & stat errrors ****************************
//**************************************************************************************************************************************
	
	TCanvas* canvasInvXSectionNLOOnlyRatioPi0Sep = new TCanvas("canvasInvXSectionNLOOnlyRatioPi0Sep","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlyRatioPi0Sep,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padXSectionNLOLegendOnlyRatioPi0Sep = new TPad("padXSectionNLOLegendOnlyRatioPi0Sep", "",  0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioPi0Sep, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlyRatioPi0Sep->Draw();

	TPad* padXSectionNLOOnlyRatioPi0SepPi07TeV = new TPad("padXSectionNLOOnlyRatioPi0SepPi07TeV", "", 0., 0.565, 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioPi0SepPi07TeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOOnlyRatioPi0SepPi07TeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioPi0SepPi02760GeV = new TPad("padXSectionNLOOnlyRatioPi0SepPi02760GeV", "", 0., 0.33, 1., 0.565,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioPi0SepPi02760GeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOOnlyRatioPi0SepPi02760GeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioPi0SepPi0900GeV = new TPad("padXSectionNLOOnlyRatioPi0SepPi0900GeV", "", 0., 0., 1., 0.33,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioPi0SepPi0900GeV, 0.15, 0.02, 0., 0.28);
	padXSectionNLOOnlyRatioPi0SepPi0900GeV->Draw();
	
	padXSectionNLOOnlyRatioPi0SepPi07TeV->cd();
	padXSectionNLOOnlyRatioPi0SepPi07TeV->SetLogx();
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOPi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionNLOPi0->GetYaxis()->SetLabelOffset(0.01);
	ratio2DInvXSectionNLOPi0->GetYaxis()->SetTickLength(0.02);
	ratio2DInvXSectionNLOPi0->GetXaxis()->SetTickLength(0.07);
	ratio2DInvXSectionNLOPi0->DrawCopy(); 
	
	graphRatioCombCombFitSys->Draw("2,same");
	graphRatioCombCombFit->Draw("p,same");
	graphRatioCombNLOPi07TeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi07TeVMuOne->Draw("same,c");
	graphRatioCombNLOPi07TeVMuTwo->Draw("same,c");
	graphRatioCombNLOBKKPi07TeVMuTwo->Draw("same,c");
	boxErrorSigmaPi07TeVRatio->Draw();
	labelRatioNLOPi07TeV->Draw();
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOOnlyRatioPi0SepPi02760GeV->cd();
	padXSectionNLOOnlyRatioPi0SepPi02760GeV->SetLogx();
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionNLOPi0900GeV->GetYaxis()->SetLabelOffset(0.01);
	ratio2DInvXSectionNLOPi0900GeV->GetYaxis()->SetTickLength(0.02);
	ratio2DInvXSectionNLOPi0900GeV->GetXaxis()->SetTickLength(0.07);
	ratio2DInvXSectionNLOPi0900GeV->DrawCopy(); 

	graphRatioCombCombFit2760GeVSys->Draw("2,same");
	graphRatioCombCombFit2760GeVStat->Draw("p,same");
	graphRatioCombNLOPi02760GeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi02760GeVMuOne->Draw("same,c");
	graphRatioCombNLOPi02760GeVMuTwo->Draw("same,c");
	boxErrorSigmaPi02760GeVRatio->Draw();
	labelRatioNLOPi02760GeV->Draw();
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOOnlyRatioPi0SepPi0900GeV->cd();
	padXSectionNLOOnlyRatioPi0SepPi0900GeV->SetLogx();
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.12,0.13, 0.125,0.128, 0.9,0.45, 510, 505);
	ratio2DInvXSectionNLOEta->GetYaxis()->SetLabelOffset(0.01);
	ratio2DInvXSectionNLOEta->GetXaxis()->SetLabelOffset(-0.03);
	ratio2DInvXSectionNLOEta->GetYaxis()->SetTickLength(0.03);
	ratio2DInvXSectionNLOEta->GetXaxis()->SetTickLength(0.06);
	ratio2DInvXSectionNLOEta->DrawCopy(); 
	
	graphRatioCombCombFit900GeVSys->Draw("2,same");
	graphRatioCombCombFit900GeVStat->Draw("p,same");
	graphRatioCombNLOBKKPi0900GeVMuTwo->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuHalf->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuOne->Draw("same,c");
	graphRatioCombNLOPi0900GeVMuTwo->Draw("same,c");
	boxErrorSigmaPi0900GeVRatio->Draw();
	labelRatioNLOPi0900GeVALLEnergiesSep->Draw();
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	
	padXSectionNLOLegendOnlyRatioPi0Sep->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioPi0Sep, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlyRatioPi0Sep->SetBorderMode(-1);
	padXSectionNLOLegendOnlyRatioPi0Sep->SetBorderSize(3);
	padXSectionNLOLegendOnlyRatioPi0Sep->Draw();
	padXSectionNLOLegendOnlyRatioPi0Sep->cd();
		
	//*************** first Column **********************************************************
	textSpectrumOnlyRatioPi0->Draw();
	textNLOMuHalfOnlyRatioPi0->Draw();
	textNLOMuOneOnlyRatioPi0->Draw();
	textNLOMuTwoOnlyRatioPi0->Draw();
	textNLOBKKMuTwoOnlyRatioPi0->Draw();
	
	//*************** second Column **********************************************************
	textPi07TeVNLOOnlyRatioPi0->Draw();
	TLatex *textPi07TeVNLOsysOnlyRatioPi0Sep = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[1],rowsNLOLegendPrelAndFinalOnlyRatioPi0[1],"syst., stat.");
	SetStyleTLatex( textPi07TeVNLOsysOnlyRatioPi0Sep, textSizeTopLowerLablesPrelAndFinalRatio,4);
	textPi07TeVNLOsysOnlyRatioPi0Sep->Draw();
	TBox* boxCombinedPi07TeVOnlyRatioPi0Sep = CreateBoxFromGraphWithFill(graphInvCrossSectionPi0Comb7TeVSysErr,columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi07TeVOnlyRatioPi0Sep->Draw("f");
	markerCombinedPi07TeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineNLOPi07TeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi07TeVMuOneOnlyRatioPi0->Draw("same");
	lineNLOPi07TeVMuTwoOnlyRatioPi0->Draw("same");
	lineNLOBKKPi07TeVMuTwoOnlyRatioPi0->Draw("same");
	
	//*************** third Column **********************************************************
	textPi02760GeVNLOOnlyRatioPi0->Draw();
	TLatex *textPi02760GeVNLOsysOnlyRatioPi0Sep = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[2],rowsNLOLegendPrelAndFinalOnlyRatioPi0[1],"syst., stat.");
	SetStyleTLatex( textPi02760GeVNLOsysOnlyRatioPi0Sep, textSizeTopLowerLablesPrelAndFinalRatio,4);
	textPi02760GeVNLOsysOnlyRatioPi0Sep->Draw();
	TBox* boxCombinedPi02760GeVOnlyRatioPi0Sep = CreateBoxFromGraphWithFill(graphInvCrossSectionPi0Comb2760GeVSysErr,columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi02760GeVOnlyRatioPi0Sep->Draw("f");
	markerCombinedPi02760GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineNLOPi02760GeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi02760GeVMuOneOnlyRatioPi0->Draw("same");
	lineNLOPi02760GeVMuTwoOnlyRatioPi0->Draw("same");
	
	//**************** forth Column **********************************************************
	textPi0900GeVNLOOnlyRatioPi0->Draw();
	TLatex *textPi0900GeVNLOsysOnlyRatioPi0Sep = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[3],rowsNLOLegendPrelAndFinalOnlyRatioPi0[1],"syst., stat.");
	SetStyleTLatex( textPi0900GeVNLOsysOnlyRatioPi0Sep, textSizeTopLowerLablesPrelAndFinalRatio,4);
	textPi0900GeVNLOsysOnlyRatioPi0Sep->Draw();
	TBox* boxCombinedPi0900GeVOnlyRatioPi0Sep = CreateBoxFromGraphWithFill(graphInvCrossSectionPi0Comb900GeVSysErr,columnsNLOLegendPrelAndFinalOnlyRatio[3]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[3]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	boxCombinedPi0900GeVOnlyRatioPi0Sep->Draw("f");
	markerCombinedPi0900GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineNLOPi0900GeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi0900GeVMuOneOnlyRatioPi0->Draw("same");
	lineNLOPi0900GeVMuTwoOnlyRatioPi0->Draw("same");
	lineNLOBKKPi0900GeVMuTwoOnlyRatioPi0->Draw("same");
	textArxivALLEnergiesRatio->Draw();
	textArxivALLEnergiesRatio2->Draw();
	
	canvasInvXSectionNLOOnlyRatioPi0Sep->Update();
	canvasInvXSectionNLOOnlyRatioPi0Sep->Print(Form("%s/Pi0_InvXSectionNLO_OnlyRatioSep_Paper.%s",outputDir.Data(),suffix.Data()));	

	
	//***************************************************************************************************************
	//************************ Eta Combined + NLO *******************************************************************
	//***************************************************************************************************************
	
	TCanvas* canvasInvXSectionEtaALLEnergies = new TCanvas("canvasInvXSectionEtaALLEnergies","",200,10,1200,2000);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionEtaALLEnergies,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionEtaALLEnergies = new TPad("padComparisonXSectionEtaALLEnergies", "", 0., 0.42, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionEtaALLEnergies, 0.15, 0.02, 0.03, 0.);
	padComparisonXSectionEtaALLEnergies->Draw();
	
	TPad* padXSectionEtaALLEnergiesRatioEta7TeV = new TPad("padXSectionEtaALLEnergiesRatioEta7TeV", "", 0., 0.30, 1., 0.42,-1, -1, -2);
	DrawGammaPadSettings( padXSectionEtaALLEnergiesRatioEta7TeV, 0.15, 0.02, 0., 0.);
	padXSectionEtaALLEnergiesRatioEta7TeV->Draw();
	
	TPad* padXSectionEtaALLEnergiesRatioEta2760GeV = new TPad("padXSectionEtaALLEnergiesRatioEta2760GeV", "", 0., 0.18, 1., 0.30,-1, -1, -2);
	DrawGammaPadSettings( padXSectionEtaALLEnergiesRatioEta2760GeV, 0.15, 0.02, 0., 0.);
	padXSectionEtaALLEnergiesRatioEta2760GeV->Draw();
	
	TPad* padXSectionEtaALLEnergiesRatioEta900GeV = new TPad("padXSectionEtaALLEnergiesRatioEta900GeV", "", 0., 0., 1., 0.18,-1, -1, -2);
	DrawGammaPadSettings( padXSectionEtaALLEnergiesRatioEta900GeV, 0.15, 0.02, 0., 0.28);
	padXSectionEtaALLEnergiesRatioEta900GeV->Draw();
	
	
	padComparisonXSectionEtaALLEnergies->cd();
	padComparisonXSectionEtaALLEnergies->SetLogy();		
	padComparisonXSectionEtaALLEnergies->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionEtaALLEnergies = new TH2F("histo2DInvXSectionEtaALLEnergies","histo2DInvXSectionEtaALLEnergies",1000,0.23,20.,1000,2,10e10);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionEtaALLEnergies, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionEtaALLEnergies->DrawCopy(); 
	
	graphInvCrossSectionEtaComb900GeV= ScaleGraph(graphInvCrossSectionEtaComb900GeV,1e-2);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb900GeV, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionEtaComb900GeV->SetLineWidth(widthCommonErrors);
	//graphInvCrossSectionEtaComb900GeV->Draw("same,p");
	graphInvCrossSectionEtaComb900GeV->Draw("p,E2same");
	
	graphInvCrossSectionEtaComb2760GeV = ScaleGraph(graphInvCrossSectionEtaComb2760GeV,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb2760GeV, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionEtaComb2760GeV->SetLineWidth(widthCommonErrors);
	//graphInvCrossSectionEtaComb2760GeV->Draw("same,p");
	graphInvCrossSectionEtaComb2760GeV->Draw("p,E2same");
	
	graphInvCrossSectionEtaComb7TeV = ScaleGraph(graphInvCrossSectionEtaComb7TeV,1e+3);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb7TeV, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionEtaComb7TeV->SetLineWidth(widthCommonErrors);
	//graphInvCrossSectionEtaComb7TeV->Draw("p,same");
	graphInvCrossSectionEtaComb7TeV->Draw("p,E2same");
	
	DrawGammaSetMarkerTF1( fitInvCrossSectionEta7TeV, styleFitCommonSpectrum, widthCommonFit, colorCommonSpectrumPi07TeV);
	fitInvCrossSectionEta7TeV->Draw("same");
	
	histoFitInvCrossSectionEta900GeV->Scale(1e-2);
	SetStyleHisto(histoFitInvCrossSectionEta900GeV, widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi0900GeV);
	histoFitInvCrossSectionEta900GeV->Draw("same,c");
	
	histoFitInvCrossSectionEta2760GeV->Scale(1e-1);
	SetStyleHisto(histoFitInvCrossSectionEta2760GeV, widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi02760GeV);
	histoFitInvCrossSectionEta2760GeV->Draw("same,c");
	
	graphNLOMuHalfEta7TeV= ScaleGraph(graphNLOMuHalfEta7TeV,1e+3);
	DrawGammaNLOTGraph( graphNLOMuHalfEta7TeV, widthCommonFit, styleLineNLOMuHalf, colorNLOPi07TeVMuHalf);
	graphNLOMuHalfEta7TeV->Draw("same,c");
	graphNLOMuOneEta7TeV= ScaleGraph(graphNLOMuOneEta7TeV,1e+3);
	DrawGammaNLOTGraph( graphNLOMuOneEta7TeV, widthCommonFit, styleLineNLOMuOne, colorNLOPi07TeVMuOne);
	graphNLOMuOneEta7TeV->Draw("same,c");
	graphNLOMuTwoEta7TeV = ScaleGraph(graphNLOMuTwoEta7TeV,1e+3);
	DrawGammaNLOTGraph( graphNLOMuTwoEta7TeV, widthCommonFit, styleLineNLOMuTwo, colorNLOPi07TeVMuTwo);
	graphNLOMuTwoEta7TeV->Draw("same,c");
	
	graphNLOMuHalfEta2760GeV= ScaleGraph(graphNLOMuHalfEta2760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuHalfEta2760GeV, widthCommonFit, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	graphNLOMuHalfEta2760GeV->Draw("same,c");
	graphNLOMuOneEta2760GeV= ScaleGraph(graphNLOMuOneEta2760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuOneEta2760GeV, widthCommonFit, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	graphNLOMuOneEta2760GeV->Draw("same,c");
	graphNLOMuTwoEta2760GeV = ScaleGraph(graphNLOMuTwoEta2760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuTwoEta2760GeV, widthCommonFit, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);
	graphNLOMuTwoEta2760GeV->Draw("same,c");
	
	graphNLOMuHalfEta900GeV= ScaleGraph(graphNLOMuHalfEta900GeV,1e-2);
	DrawGammaNLOTGraph( graphNLOMuHalfEta900GeV, widthCommonFit, styleLineNLOMuHalf, colorNLOEta900GeVMuHalf);
	graphNLOMuHalfEta900GeV->Draw("same,c");
	graphNLOMuOneEta900GeV= ScaleGraph(graphNLOMuOneEta900GeV,1e-2);
	DrawGammaNLOTGraph( graphNLOMuOneEta900GeV, widthCommonFit, styleLineNLOMuOne, colorNLOEta900GeVMuOne);
	graphNLOMuOneEta900GeV->Draw("same,c");
	graphNLOMuTwoEta900GeV = ScaleGraph(graphNLOMuTwoEta900GeV,1e-2);
	DrawGammaNLOTGraph( graphNLOMuTwoEta900GeV, widthCommonFit, styleLineNLOMuTwo, colorNLOEta900GeVMuTwo);
	graphNLOMuTwoEta900GeV->Draw("same,c");
	
	// 		TLatex *labelScalingEta7TeVAllEnergiesSp = new TLatex(0.33,5E9,"x 1");
// 	SetStyleTLatex( labelScalingEta7TeVAllEnergiesSp, 0.028,4,colorCommonSpectrumPi07TeV,kFALSE);
// 	labelScalingEta7TeVAllEnergiesSp->Draw();
// 	TLatex *labelScalingEta2760GeVAllEnergiesSp = new TLatex(0.46,3E8,"x 10^{-1}");
// 	SetStyleTLatex( labelScalingEta2760GeVAllEnergiesSp, 0.028,4,colorCommonSpectrumPi02760GeV,kFALSE);
// 	labelScalingEta2760GeVAllEnergiesSp->Draw();
// 	TLatex *labelScalingEta900GeVAllEnergiesSp = new TLatex(0.71,2.5E6,"x 10^{-2}");
// 	SetStyleTLatex( labelScalingEta900GeVAllEnergiesSp, 0.028,4,colorCommonSpectrumPi0900GeV,kFALSE);
// 	labelScalingEta900GeVAllEnergiesSp->Draw();

	TLatex *labelScalingEta7TeVALLEnergies = new TLatex(0.33,5E9,"x 1");
	SetStyleTLatex( labelScalingEta7TeVALLEnergies, 0.025,4,fitInvCrossSectionPi0->GetLineColor(),62,kFALSE);
	labelScalingEta7TeVALLEnergies->Draw();

	TLatex *labelScalingEta2760GeVALLEnergies = new TLatex(0.46,3E8,"x 10^{-1}");
	SetStyleTLatex( labelScalingEta2760GeVALLEnergies, 0.025,4,histoFitInvCrossSectionPi02760GeV->GetLineColor(),62,kFALSE);
	labelScalingEta2760GeVALLEnergies->Draw();

	TLatex *labelScalingEta900GeVALLEnergies = new TLatex(0.70,2.5E7,"x 10^{-2}");
	SetStyleTLatex( labelScalingEta900GeVALLEnergies, 0.025,4,histoFitInvCrossSectionPi0900GeV->GetLineColor(),62,kFALSE);
	labelScalingEta900GeVALLEnergies->Draw();

	DrawNormalizationErrorText(normalizationInvX1MesonALLEn[0],normalizationInvX1MesonALLEn[1],normalizationInvX1MesonALLEn[2],
						  normalizationInvX1MesonALLEn[3],normalizationInvX1MesonALLEn[4],"all"); 
			
	TPad* padXSectionEtaALLEnergiesLegend = new TPad("padXSectionEtaALLEnergiesLegend", "", 0.17, 0.005, 0.95, 0.21,-1, -1, -2); 
	DrawGammaPadSettings( padXSectionEtaALLEnergiesLegend, 0., 0., 0., 0.);
	padXSectionEtaALLEnergiesLegend->Draw();
	padXSectionEtaALLEnergiesLegend->cd();
	
	//*************** first Column **********************************************************
	textSpectrumALLEnergies->Draw();
	textFitCombALLEnergies->Draw();
	textNLOMuHalfALLEnergies->Draw();
	textNLOMuOneALLEnergies->Draw();
	textNLOMuTwoALLEnergies->Draw();
	
	//*************** second Column **********************************************************
	TLatex *textEta7TeVNLOEtaALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[1],rowsNLOLegendPrelAndFinal[0],"#eta, #sqrt{#it{s}} = 7 TeV (*)");
	SetStyleTLatex( textEta7TeVNLOEtaALLEnergies, textSizeTopLablesPrelAndFinal,4);
	textEta7TeVNLOEtaALLEnergies->Draw();
	textPi07TeVNLOsysALLEnergies->Draw();
	boxCombinedPi07TeVALLEnergies->Draw("l");
	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFit7TeVNLOALLEnergies->Draw("same");
	lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
	lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
	//lineNLOPi07TeVMuTwoBKKALLEnergies->Draw("same");
	
	//*************** third Column **********************************************************
	TLatex *textEta2760GeVNLOEtaALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[2],rowsNLOLegendPrelAndFinal[0],"#eta, #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( textEta2760GeVNLOEtaALLEnergies, textSizeTopLablesPrelAndFinal,4);
	textEta2760GeVNLOEtaALLEnergies->Draw();
	textPi02760GeVNLOsysALLEnergies->Draw();
	boxCombinedPi02760GeVALLEnergies->Draw("l");
	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi02760GeVNLOALLEnergies->Draw("same");
	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
		
	//**************** forth Column **********************************************************
	TLatex *textEta900GeVNLOEtaALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[3],rowsNLOLegendPrelAndFinal[0],"#eta, #sqrt{#it{s}} = 0.9 TeV");
	SetStyleTLatex( textEta900GeVNLOEtaALLEnergies, textSizeTopLablesPrelAndFinal,4);
	textEta900GeVNLOEtaALLEnergies->Draw();
	textPi0900GeVNLOsysALLEnergies->Draw();
	boxCombinedPi0900GeVALLEnergies->Draw("l");
	markerCombinedPi0900GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi0900GeVNLOALLEnergies->Draw("same");
	lineNLOPi0900GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi0900GeVMuOneALLEnergies->Draw("same");
	lineNLOPi0900GeVMuTwoALLEnergies->Draw("same");
	//lineNLOPi0900GeVMuTwoBKKALLEnergies->Draw("same");
	textArxivALLEnergies->Draw();
	
	//************************************************* End Legend ***************************************************
	
	padXSectionEtaALLEnergiesRatioEta7TeV->cd();
	padXSectionEtaALLEnergiesRatioEta7TeV->SetLogx();
	
	TH2F * ratio2DInvXSectionEtaALLEnergiesEta = new TH2F("ratio2DInvXSectionEtaALLEnergiesEta","ratio2DInvXSectionEtaALLEnergiesEta",1000,0.23,20.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtaALLEnergiesEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionEtaALLEnergiesEta->DrawCopy(); 
	DrawGammaNLOTGraph( graphRatioCombNLOEta7TeVMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLOPi07TeVMuHalf);
	graphRatioCombNLOEta7TeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta7TeVMuOne, widthCommonFit, styleLineNLOMuOne, colorNLOPi07TeVMuOne);
	graphRatioCombNLOEta7TeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta7TeVMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLOPi07TeVMuTwo);
	graphRatioCombNLOEta7TeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi07TeVRatio->Draw();

	TLatex *labelRatioNLOEta7teVALLEnergies = new TLatex(0.18,0.75,"#eta, #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( labelRatioNLOEta7teVALLEnergies, 0.17,4);
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta7TeV, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFitEta7TeV->SetLineWidth(widthCommonErrors);
	//graphRatioCombCombFitEta7TeV->Draw("p,same");
	graphRatioCombCombFitEta7TeV->Draw("p,E2same");
	labelRatioNLOEta7teVALLEnergies->Draw();

	DrawGammaLines(0., 20.,1., 1.,0.1);
	
	padXSectionEtaALLEnergiesRatioEta2760GeV->cd();
	padXSectionEtaALLEnergiesRatioEta2760GeV->SetLogx();
	TH2F * ratio2DInvXSectionEtaALLEnergiesEta2760GeV;
	ratio2DInvXSectionEtaALLEnergiesEta2760GeV = new TH2F("ratio2DInvXSectionEtaALLEnergiesEta2760GeV","ratio2DInvXSectionEtaALLEnergiesEta2760GeV",1000,0.23,20.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtaALLEnergiesEta2760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionEtaALLEnergiesEta2760GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta2760GeV, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFitEta2760GeV->SetLineWidth(widthCommonErrors);
	//graphRatioCombCombFitEta2760GeV->Draw("same,p");
	graphRatioCombCombFitEta2760GeV->Draw("p,E2same");
	
	DrawGammaNLOTGraph( graphRatioCombNLOEta2760GeVMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	graphRatioCombNLOEta2760GeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta2760GeVMuOne, widthCommonFit, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	graphRatioCombNLOEta2760GeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta2760GeVMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);
	graphRatioCombNLOEta2760GeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi02760GeVRatio->Draw();
	
	TLatex *labelRatioNLOEta2760GeV = new TLatex(0.18,0.75,"#eta, #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( labelRatioNLOEta2760GeV, 0.17,4);
	labelRatioNLOEta2760GeV->Draw();
	
	DrawGammaLines(0., 20.,1., 1.,0.1);
	
	padXSectionEtaALLEnergiesRatioEta900GeV->cd();
	padXSectionEtaALLEnergiesRatioEta900GeV->SetLogx();
	TH2F * ratio2DInvXSectionEtaALLEnergiesEta900GeV=  new TH2F("ratio2DInvXSectionEtaALLEnergiesEta900GeV","ratio2DInvXSectionEtaALLEnergiesEta900GeV",1000,0.23,20.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtaALLEnergiesEta900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.12,0.13, 0.125,0.128, 0.9,0.45, 510, 505);
	ratio2DInvXSectionEtaALLEnergiesEta900GeV->GetXaxis()->SetLabelOffset(-0.03);
	ratio2DInvXSectionEtaALLEnergiesEta900GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta900GeV, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioCombCombFitEta900GeV->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFitEta900GeV->Draw("p,E2same");
	
	DrawGammaNLOTGraph( graphRatioCombNLOEta900GeVMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLOPi0900GeVMuHalf);
	graphRatioCombNLOEta900GeVMuHalf->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta900GeVMuOne, widthCommonFit, styleLineNLOMuOne, colorNLOPi0900GeVMuOne);
	graphRatioCombNLOEta900GeVMuOne->Draw("same,c");
	DrawGammaNLOTGraph( graphRatioCombNLOEta900GeVMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLOPi0900GeVMuTwo);
	graphRatioCombNLOEta900GeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi0900GeVRatio->Draw();
	TLatex *labelRatioNLOEta900GeVEtaALLEnergies = new TLatex(0.18,0.85,"#eta, #sqrt{#it{s}} = 0.9 TeV");
	SetStyleTLatex( labelRatioNLOEta900GeVEtaALLEnergies, 0.115,4);
	labelRatioNLOEta900GeVEtaALLEnergies->Draw();
	
	DrawGammaLines(0., 20.,1., 1.,0.1);
	
	canvasInvXSectionEtaALLEnergies->Update();
	canvasInvXSectionEtaALLEnergies->Print(Form("%s/Eta_InvXSectionALLEnergies_Paper.%s",outputDir.Data(),suffix.Data()));
	padComparisonXSectionEtaALLEnergies->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSPrelim(pictDrawingCoordinatesCombineMeasCross[0], pictDrawingCoordinatesCombineMeasCross[1], pictDrawingCoordinatesCombineMeasCross[2], pictDrawingCoordinatesCombineMeasCross[3], pictDrawingCoordinatesCombineMeasCross[4], pictDrawingCoordinatesCombineMeasCross[5], pictDrawingCoordinatesCombineMeasCross[6], pictDrawingCoordinatesCombineMeasCross[7], pictDrawingCoordinates[8],collisionSystemCombinedReallyAll, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1200,1220);

	canvasInvXSectionEtaALLEnergies->Update();
	canvasInvXSectionEtaALLEnergies->Print(Form("%s/Eta_InvXSectionALLEnergies.%s",outputDir.Data(),suffix.Data()));

	//*****************************************************************************************************************
	//************************ Eta Meson Combined + NLO just spectrum******************************************************
	//*****************************************************************************************************************

	TCanvas* canvasInvXSectionNLOOnlySpectraEta = new TCanvas("canvasInvXSectionNLOOnlySpectraEta","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlySpectraEta,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionNLOOnlySpectraEta = new TPad("padComparisonXSectionNLOOnlySpectraEta", "", 0., 0., 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionNLOOnlySpectraEta, 0.15, 0.02, 0.02, 0.09);
	padComparisonXSectionNLOOnlySpectraEta->Draw();
	
	TPad* padXSectionNLOLegendOnlySpectraEta = new TPad("padXSectionNLOLegendOnlySpectraEta", "", 0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraEta, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlySpectraEta->Draw();
	
	padComparisonXSectionNLOOnlySpectraEta->cd();
	padComparisonXSectionNLOOnlySpectraEta->SetLogy();		
	padComparisonXSectionNLOOnlySpectraEta->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionNLOOnlySpectraEta = new TH2F("histo2DInvXSectionNLOOnlySpectraEta","histo2DInvXSectionNLOOnlySpectraEta",1000,0.23,30.,1000,1e2,1e11);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLOOnlySpectraEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionNLOOnlySpectraEta->DrawCopy(); 
	
	graphInvCrossSectionEtaComb900GeV->Draw("p,E2same");
	graphInvCrossSectionEtaComb7TeV->Draw("p,E2same");
	fitInvCrossSectionEta7TeV->Draw("same");
	graphInvCrossSectionEtaComb2760GeV->Draw("p,E2same");
	histoFitInvCrossSectionEta900GeV->Draw("same,c");
	histoFitInvCrossSectionEta2760GeV ->Draw("same,c");
	
	graphNLOMuHalfEta7TeV->Draw("same,c");
	graphNLOMuOneEta7TeV->Draw("same,c");
	graphNLOMuTwoEta7TeV->Draw("same,c");
	
	graphNLOMuHalfEta900GeV->Draw("same,c");
	graphNLOMuOneEta900GeV->Draw("same,c");
	graphNLOMuTwoEta900GeV->Draw("same,c");
	
	graphNLOMuHalfEta2760GeV->Draw("same,c");
	graphNLOMuOneEta2760GeV->Draw("same,c");
	graphNLOMuTwoEta2760GeV->Draw("same,c");

	labelScalingEta7TeVALLEnergies->Draw();
	labelScalingEta900GeVALLEnergies->Draw();
	labelScalingEta2760GeVALLEnergies->Draw();

	DrawNormalizationErrorText(normalizationInvXOnlySpec[0]+0.1,normalizationInvXOnlySpec[1],normalizationInvXOnlySpec[2],
						  normalizationInvXOnlySpec[3],normalizationInvXOnlySpec[4],"all"); 
	
	//************************************************* Begin Legend ***************************************************

	padXSectionNLOLegendOnlySpectraEta->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraEta, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlySpectraEta->SetBorderMode(-1);
	padXSectionNLOLegendOnlySpectraEta->SetBorderSize(3);
	padXSectionNLOLegendOnlySpectraEta->Draw();
	padXSectionNLOLegendOnlySpectraEta->cd();

	//*************** first Column **********************************************************
	textSpectrumALLEnergies->Draw();
	textFitCombALLEnergies->Draw();
	textNLOMuHalfALLEnergies->Draw();
	textNLOMuOneALLEnergies->Draw();
	textNLOMuTwoALLEnergies->Draw();
	
	//*************** second Column **********************************************************
	textEta7TeVNLOEtaALLEnergies->Draw();
	textPi07TeVNLOsysALLEnergies->Draw();
	boxCombinedPi07TeVALLEnergies->Draw("l");
	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFit7TeVNLOALLEnergies->Draw("same");
	lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
	lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
	
	//*************** third Column **********************************************************
	textEta2760GeVNLOEtaALLEnergies->Draw();
	textPi02760GeVNLOsysALLEnergies->Draw();
	boxCombinedPi02760GeVALLEnergies->Draw("l");
	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi02760GeVNLOALLEnergies->Draw("same");
	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
		
	//**************** forth Column **********************************************************
	textEta900GeVNLOEtaALLEnergies->Draw();
	textPi0900GeVNLOsysALLEnergies->Draw();
	boxCombinedPi0900GeVALLEnergies->Draw("l");
	markerCombinedPi0900GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi0900GeVNLOALLEnergies->Draw("same");
	lineNLOPi0900GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi0900GeVMuOneALLEnergies->Draw("same");
	lineNLOPi0900GeVMuTwoALLEnergies->Draw("same");
	textArxivALLEnergies->Draw();
	
	
	canvasInvXSectionNLOOnlySpectraEta->Update();
	canvasInvXSectionNLOOnlySpectraEta->Print(Form("%s/Eta_InvXSectionNLO_OnlySpectrum_Paper.%s",outputDir.Data(),suffix.Data()));
	
	padComparisonXSectionNLOOnlySpectraEta->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSPrelim(pictDrawingCoordinatesOnlySpectrum[0]-0.05, pictDrawingCoordinatesOnlySpectrum[1], pictDrawingCoordinatesOnlySpectrum[2], pictDrawingCoordinatesOnlySpectrum[3], pictDrawingCoordinatesOnlySpectrum[4], pictDrawingCoordinatesOnlySpectrum[5], pictDrawingCoordinatesOnlySpectrum[6], pictDrawingCoordinatesOnlySpectrum[7], pictDrawingCoordinates[8],collisionSystemCombinedReallyAll, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1200,960);
	
	canvasInvXSectionNLOOnlySpectraEta->Update();
	canvasInvXSectionNLOOnlySpectraEta->Print(Form("%s/Eta_InvXSectionNLO_OnlySpectrum.%s",outputDir.Data(),suffix.Data()));
	
//**************************************************************************************************************************************
//************************************ Inv Crosssection + NLO only ratio Pi0 ***************************************************************
//**************************************************************************************************************************************
	
	TCanvas* canvasInvXSectionNLOOnlyRatioEta = new TCanvas("canvasInvXSectionNLOOnlyRatioEta","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlyRatioEta,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padXSectionNLOLegendOnlyRatioEta = new TPad("padXSectionNLOLegendOnlyRatioEta", "", 0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioEta, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlyRatioEta->Draw();

	TPad* padXSectionNLOOnlyRatioEtaEta7TeV = new TPad("padXSectionNLOOnlyRatioEtaEta7TeV", "", 0., 0.565, 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioEtaEta7TeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOOnlyRatioEtaEta7TeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioEtaEta2760GeV = new TPad("padXSectionNLOOnlyRatioEtaEta2760GeV", "", 0., 0.33, 1., 0.565,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioEtaEta2760GeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOOnlyRatioEtaEta2760GeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioEtaEta900GeV = new TPad("padXSectionNLOOnlyRatioEtaEta900GeV", "", 0., 0., 1., 0.33,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioEtaEta900GeV, 0.15, 0.02, 0., 0.28);
	padXSectionNLOOnlyRatioEtaEta900GeV->Draw();
	
	padXSectionNLOOnlyRatioEtaEta7TeV->cd();
	padXSectionNLOOnlyRatioEtaEta7TeV->SetLogx();
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOPi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionNLOPi0->DrawCopy(); 
	graphRatioCombNLOEta7TeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta7TeVMuOne->Draw("same,c");
	graphRatioCombNLOEta7TeVMuTwo->Draw("same,c");
	boxErrorSigmaPi07TeVRatio->Draw();
	graphRatioCombCombFitEta7TeV ->Draw("p,E2same");
	labelRatioNLOEta7teVALLEnergies->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOOnlyRatioEtaEta2760GeV->cd();
	padXSectionNLOOnlyRatioEtaEta2760GeV->SetLogx();
	
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionNLOPi0900GeV->DrawCopy(); 
	
	graphRatioCombNLOEta2760GeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta2760GeVMuOne->Draw("same,c");
	graphRatioCombNLOEta2760GeVMuTwo->Draw("same,c");
	boxErrorSigmaPi02760GeVRatio->Draw();
	graphRatioCombCombFitEta2760GeV->Draw("p,E2same");
	labelRatioNLOEta2760GeV->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOOnlyRatioEtaEta900GeV->cd();
	padXSectionNLOOnlyRatioEtaEta900GeV->SetLogx();
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.12,0.13, 0.125,0.128, 0.9,0.45, 510, 505);
	ratio2DInvXSectionNLOEta->GetXaxis()->SetLabelOffset(-0.03);
	ratio2DInvXSectionNLOEta->DrawCopy(); 
	
	
	graphRatioCombNLOEta900GeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta900GeVMuOne->Draw("same,c");
	graphRatioCombNLOEta900GeVMuTwo->Draw("same,c");
	
	//graphRatioCombCombFitPi0900GeV->Draw("same,p");
	
	graphRatioCombCombFitEta900GeV->Draw("p,E2same");
	boxErrorSigmaPi0900GeVRatio->Draw();
	labelRatioNLOEta900GeVEtaALLEnergies->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	
	padXSectionNLOLegendOnlyRatioEta->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioEta, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlyRatioEta->SetBorderMode(-1);
	padXSectionNLOLegendOnlyRatioEta->SetBorderSize(3);
	padXSectionNLOLegendOnlyRatioEta->Draw();
	padXSectionNLOLegendOnlyRatioEta->cd();
	
	//****************************** Definition of the Legend ******************************************
	//**************** Row def ************************
	Double_t rowsNLOLegendPrelAndFinalOnlyRatioEta[7] = {0.88,0.75,0.65,0.5,0.37,0.22,0.07};
	
	//*************** first Column **********************************************************
	textSpectrumOnlyRatioPi0->Draw();
	textNLOMuHalfOnlyRatioPi0->Draw();
	textNLOMuOneOnlyRatioPi0->Draw();
	textNLOMuTwoOnlyRatioPi0->Draw();	
	
	//*************** second Column **********************************************************
	TLatex *textEta7TeVNLOOnlyRatioEta = new TLatex(columnsNLOLegendPrelAndFinal[1],rowsNLOLegendPrelAndFinalOnlyRatioEta[0],"#eta, #sqrt{#it{s}} = 7 TeV (*)");
	SetStyleTLatex( textEta7TeVNLOOnlyRatioEta, textSizeTopLablesPrelAndFinal,4);
	textEta7TeVNLOOnlyRatioEta->Draw();
	textPi07TeVNLOsysOnlyRatioPi0->Draw();
	boxCombinedPi07TeVOnlyRatioPi0->Draw("l");
	markerCombinedPi07TeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineNLOPi07TeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi07TeVMuOneOnlyRatioPi0->Draw("same");
	lineNLOPi07TeVMuTwoOnlyRatioPi0->Draw("same");
		
	
	//*************** third Column **********************************************************
	TLatex *textEta2760GeVNLOOnlyRatioEta = new TLatex(columnsNLOLegendPrelAndFinal[2],rowsNLOLegendPrelAndFinalOnlyRatioEta[0],"#eta, #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( textEta2760GeVNLOOnlyRatioEta, textSizeTopLablesPrelAndFinal,4);
	textEta2760GeVNLOOnlyRatioEta->Draw();
	textPi02760GeVNLOsysOnlyRatioPi0->Draw();
	boxCombinedPi02760GeVOnlyRatioPi0->Draw("l");
	markerCombinedPi02760GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineNLOPi02760GeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi02760GeVMuOneOnlyRatioPi0->Draw("same");
	lineNLOPi02760GeVMuTwoOnlyRatioPi0->Draw("same");
	
	//**************** forth Column **********************************************************
	TLatex *textEta900GeVNLOOnlyRatioEta = new TLatex(columnsNLOLegendPrelAndFinal[3],rowsNLOLegendPrelAndFinalOnlyRatioEta[0],"#eta, #sqrt{#it{s}} = 0.9 TeV");
	SetStyleTLatex( textEta900GeVNLOOnlyRatioEta, textSizeTopLablesPrelAndFinal,4);
	textEta900GeVNLOOnlyRatioEta->Draw();
	textPi0900GeVNLOsysOnlyRatioPi0->Draw();
	boxCombinedPi0900GeVOnlyRatioPi0->Draw("l");
	markerCombinedPi0900GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineNLOPi0900GeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi0900GeVMuOneOnlyRatioPi0->Draw("same");
	lineNLOPi0900GeVMuTwoOnlyRatioPi0->Draw("same");
	textArxivALLEnergiesRatio->Draw();

	//************************************************* End Legend ***************************************************
	canvasInvXSectionNLOOnlyRatioEta->Update();
	canvasInvXSectionNLOOnlyRatioEta->Print(Form("%s/Eta_InvXSectionNLO_OnlyRatio_Paper.%s",outputDir.Data(),suffix.Data()));
	padXSectionNLOOnlyRatioEtaEta900GeV->cd();
	DrawAliceLogoSimplePreliminary(pictDrawingCoordinatesOnlyRatio[4],pictDrawingCoordinatesOnlyRatio[5], pictDrawingCoordinatesOnlyRatio[6], pictDrawingCoordinatesOnlyRatio[7], 1200, 396);
	
	canvasInvXSectionNLOOnlyRatioEta->Update();
	canvasInvXSectionNLOOnlyRatioEta->Print(Form("%s/Eta_InvXSectionNLO_OnlyRatio.%s",outputDir.Data(),suffix.Data()));
	
	
	//***************************************************************************************************************
	//************************ Eta Combined + NLO + separated sys & stat errors *************************************
	//***************************************************************************************************************
	
	TCanvas* canvasInvXSectionEtaALLEnergiesSep = new TCanvas("canvasInvXSectionEtaALLEnergiesSep","",200,10,1200,2000);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionEtaALLEnergiesSep,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padComparisonXSectionEtaALLEnergiesSep = new TPad("padComparisonXSectionEtaALLEnergiesSep", "", 0., 0.42, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionEtaALLEnergiesSep, 0.15, 0.02, 0.03, 0.);
	padComparisonXSectionEtaALLEnergiesSep->Draw();
	
	TPad* padXSectionEtaALLEnergiesSepRatioEta7TeV = new TPad("padXSectionEtaALLEnergiesSepRatioEta7TeV", "", 0., 0.30, 1., 0.42,-1, -1, -2);
	DrawGammaPadSettings( padXSectionEtaALLEnergiesSepRatioEta7TeV, 0.15, 0.02, 0., 0.);
	padXSectionEtaALLEnergiesSepRatioEta7TeV->Draw();
	
	TPad* padXSectionEtaALLEnergiesSepRatioEta2760GeV = new TPad("padXSectionEtaALLEnergiesSepRatioEta2760GeV", "", 0., 0.18, 1., 0.30,-1, -1, -2);
	DrawGammaPadSettings( padXSectionEtaALLEnergiesSepRatioEta2760GeV, 0.15, 0.02, 0., 0.);
	padXSectionEtaALLEnergiesSepRatioEta2760GeV->Draw();
	
	TPad* padXSectionEtaALLEnergiesSepRatioEta900GeV = new TPad("padXSectionEtaALLEnergiesSepRatioEta900GeV", "", 0., 0., 1., 0.18,-1, -1, -2);
	DrawGammaPadSettings( padXSectionEtaALLEnergiesSepRatioEta900GeV, 0.15, 0.02, 0., 0.28);
	padXSectionEtaALLEnergiesSepRatioEta900GeV->Draw();
	
	
	padComparisonXSectionEtaALLEnergiesSep->cd();
	padComparisonXSectionEtaALLEnergiesSep->SetLogy();		
	padComparisonXSectionEtaALLEnergiesSep->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionEtaALLEnergiesSep = new TH2F("histo2DInvXSectionEtaALLEnergiesSep","histo2DInvXSectionEtaALLEnergiesSep",1000,0.23,20.,1000,2,10e10);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionEtaALLEnergiesSep, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionEtaALLEnergiesSep->DrawCopy(); 
	
	graphInvCrossSectionEtaComb900GeVStatErr= ScaleGraph(graphInvCrossSectionEtaComb900GeVStatErr,1e-2);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb900GeVStatErr, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionEtaComb900GeVStatErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb900GeVSysErr= ScaleGraph(graphInvCrossSectionEtaComb900GeVSysErr,1e-2);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb900GeVSysErr, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi0900GeVBox );
	graphInvCrossSectionEtaComb900GeVSysErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb900GeVSysErr->Draw("2same");
	graphInvCrossSectionEtaComb900GeVStatErr->Draw("psame");
	
	graphInvCrossSectionEtaComb2760GeVStatErr = ScaleGraph(graphInvCrossSectionEtaComb2760GeVStatErr,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb2760GeVStatErr, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionEtaComb2760GeVStatErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb2760GeVSysErr = ScaleGraph(graphInvCrossSectionEtaComb2760GeVSysErr,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb2760GeVSysErr, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE,colorCommonSpectrumPi02760GeVBox );
	graphInvCrossSectionEtaComb2760GeVSysErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb2760GeVSysErr->Draw("2same");
	graphInvCrossSectionEtaComb2760GeVStatErr->Draw("p,same");
	
	graphInvCrossSectionEtaComb7TeVStatErr = ScaleGraph(graphInvCrossSectionEtaComb7TeVStatErr,1e+3);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb7TeVStatErr, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionEtaComb7TeVStatErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb7TeVSysErr = ScaleGraph(graphInvCrossSectionEtaComb7TeVSysErr,1e+3);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb7TeVSysErr, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi07TeVBox);
	graphInvCrossSectionEtaComb7TeVSysErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb7TeVSysErr->Draw("2same");
	graphInvCrossSectionEtaComb7TeVStatErr->Draw("psame");
	
	fitInvCrossSectionEta7TeV->Draw("same");
	histoFitInvCrossSectionEta900GeV->Draw("same,c");
	histoFitInvCrossSectionEta2760GeV->Draw("same,c");
	
	graphNLOMuHalfEta7TeV->Draw("same,c");
	graphNLOMuOneEta7TeV->Draw("same,c");
	graphNLOMuTwoEta7TeV->Draw("same,c");
	
	graphNLOMuHalfEta2760GeV->Draw("same,c");
	graphNLOMuOneEta2760GeV->Draw("same,c");
	graphNLOMuTwoEta2760GeV->Draw("same,c");
	
	graphNLOMuHalfEta900GeV->Draw("same,c");
	graphNLOMuOneEta900GeV->Draw("same,c");
	graphNLOMuTwoEta900GeV->Draw("same,c");
	
	labelScalingEta7TeVALLEnergies->Draw();
	labelScalingEta2760GeVALLEnergies->Draw();
	labelScalingEta900GeVALLEnergies->Draw();

	DrawNormalizationErrorText(normalizationInvX1MesonALLEn[0],normalizationInvX1MesonALLEn[1],normalizationInvX1MesonALLEn[2],
						  normalizationInvX1MesonALLEn[3],normalizationInvX1MesonALLEn[4],"all"); 
			
	TPad* padXSectionEtaALLEnergiesLegendSep = new TPad("padXSectionEtaALLEnergiesLegendSep", "", 0.17, 0.005, 0.95, 0.21,-1, -1, -2); 
	DrawGammaPadSettings( padXSectionEtaALLEnergiesLegendSep, 0., 0., 0., 0.);
	padXSectionEtaALLEnergiesLegendSep->Draw();
	padXSectionEtaALLEnergiesLegendSep->cd();
	
	//*************** first Column **********************************************************
	textSpectrumALLEnergies->Draw();
	textFitCombALLEnergies->Draw();
	textNLOMuHalfALLEnergies->Draw();
	textNLOMuOneALLEnergies->Draw();
	textNLOMuTwoALLEnergies->Draw();
	
	//*************** second Column **********************************************************
	textEta7TeVNLOEtaALLEnergies->Draw();
	textPi07TeVNLOsysALLEnergiesSep->Draw();
	boxCombinedPi07TeVALLEnergiesSep->Draw("f");
	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFit7TeVNLOALLEnergies->Draw("same");
	lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
	lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
	//lineNLOPi07TeVMuTwoBKKALLEnergies->Draw("same");
	
	//*************** third Column **********************************************************
	textEta2760GeVNLOEtaALLEnergies->Draw();
	textPi02760GeVNLOsysALLEnergiesSep->Draw();
	boxCombinedPi02760GeVALLEnergiesSep->Draw("f");
	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi02760GeVNLOALLEnergies->Draw("same");
	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
		
	//**************** forth Column **********************************************************
	textEta900GeVNLOEtaALLEnergies->Draw();
	textPi0900GeVNLOsysALLEnergiesSep->Draw();
	boxCombinedPi0900GeVALLEnergiesSep->Draw("f");
	markerCombinedPi0900GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi0900GeVNLOALLEnergies->Draw("same");
	lineNLOPi0900GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi0900GeVMuOneALLEnergies->Draw("same");
	lineNLOPi0900GeVMuTwoALLEnergies->Draw("same");
	//lineNLOPi0900GeVMuTwoBKKALLEnergies->Draw("same");
	textArxivALLEnergies->Draw();
	
	
	padXSectionEtaALLEnergiesSepRatioEta7TeV->cd();
	padXSectionEtaALLEnergiesSepRatioEta7TeV->SetLogx();
	
	TH2F * ratio2DInvXSectionEtaALLEnergiesSepEta = new TH2F("ratio2DInvXSectionEtaALLEnergiesSepEta","ratio2DInvXSectionEtaALLEnergiesSepEta",1000,0.23,20.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtaALLEnergiesSepEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionEtaALLEnergiesSepEta->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta7TeVStat, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kFALSE);
	graphRatioCombCombFitEta7TeVStat->SetLineWidth(widthCommonErrors);	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta7TeVSys, markerStyleCommmonSpectrum7TeV,markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi07TeVBox);
	graphRatioCombCombFitEta7TeVSys->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFitEta7TeVSys->Draw("2same");
	graphRatioCombCombFitEta7TeVStat->Draw("p,same");
	
	graphRatioCombNLOEta7TeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta7TeVMuOne->Draw("same,c");
	graphRatioCombNLOEta7TeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi07TeVRatio->Draw();
	labelRatioNLOEta7teVALLEnergies->Draw();

	DrawGammaLines(0., 20.,1., 1.,0.1);
	
	padXSectionEtaALLEnergiesSepRatioEta2760GeV->cd();
	padXSectionEtaALLEnergiesSepRatioEta2760GeV->SetLogx();
	TH2F * ratio2DInvXSectionEtaALLEnergiesSepEta2760GeV;
	ratio2DInvXSectionEtaALLEnergiesSepEta2760GeV = new TH2F("ratio2DInvXSectionEtaALLEnergiesSepEta2760GeV","ratio2DInvXSectionEtaALLEnergiesSepEta2760GeV",1000,0.23,20.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtaALLEnergiesSepEta2760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionEtaALLEnergiesSepEta2760GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta2760GeVStat, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kFALSE);
	graphRatioCombCombFitEta2760GeVStat->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta2760GeVSys, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi02760GeVBox);
	graphRatioCombCombFitEta2760GeVSys->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFitEta2760GeVSys->Draw("2same");
	graphRatioCombCombFitEta2760GeVStat->Draw("p,same");
	
	graphRatioCombNLOEta2760GeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta2760GeVMuOne->Draw("same,c");
	graphRatioCombNLOEta2760GeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi02760GeVRatio->Draw();
	labelRatioNLOEta2760GeV->Draw();
	
	DrawGammaLines(0., 20.,1., 1.,0.1);
	
	padXSectionEtaALLEnergiesSepRatioEta900GeV->cd();
	padXSectionEtaALLEnergiesSepRatioEta900GeV->SetLogx();
	TH2F * ratio2DInvXSectionEtaALLEnergiesSepEta900GeV=  new TH2F("ratio2DInvXSectionEtaALLEnergiesSepEta900GeV","ratio2DInvXSectionEtaALLEnergiesSepEta900GeV",1000,0.23,20.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtaALLEnergiesSepEta900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.12,0.13, 0.125,0.128, 0.9,0.45, 510, 505);
	ratio2DInvXSectionEtaALLEnergiesSepEta900GeV->GetXaxis()->SetLabelOffset(-0.03);
	ratio2DInvXSectionEtaALLEnergiesSepEta900GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta900GeVStat, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kFALSE);
	graphRatioCombCombFitEta900GeVStat->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitEta900GeVSys, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi0900GeVBox);
	graphRatioCombCombFitEta900GeVSys->SetLineWidth(widthCommonErrors);
	graphRatioCombCombFitEta900GeVSys->Draw("2same");
	graphRatioCombCombFitEta900GeVStat->Draw("p,same");
	
	graphRatioCombNLOEta900GeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta900GeVMuOne->Draw("same,c");
	graphRatioCombNLOEta900GeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi0900GeVRatio->Draw();
	labelRatioNLOEta900GeVEtaALLEnergies->Draw();
	
	DrawGammaLines(0., 20.,1., 1.,0.1);
	
	canvasInvXSectionEtaALLEnergiesSep->Update();
	canvasInvXSectionEtaALLEnergiesSep->Print(Form("%s/Eta_InvXSectionALLEnergiesSep_Paper.%s",outputDir.Data(),suffix.Data()));
	padComparisonXSectionEtaALLEnergiesSep->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSPrelim(pictDrawingCoordinatesCombineMeasCross[0], pictDrawingCoordinatesCombineMeasCross[1], pictDrawingCoordinatesCombineMeasCross[2], pictDrawingCoordinatesCombineMeasCross[3], pictDrawingCoordinatesCombineMeasCross[4], pictDrawingCoordinatesCombineMeasCross[5], pictDrawingCoordinatesCombineMeasCross[6], pictDrawingCoordinatesCombineMeasCross[7], pictDrawingCoordinates[8],collisionSystemCombinedReallyAll, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1200,1220);

	canvasInvXSectionEtaALLEnergiesSep->Update();
	canvasInvXSectionEtaALLEnergiesSep->Print(Form("%s/Eta_InvXSectionALLEnergiesSep.%s",outputDir.Data(),suffix.Data()));

	//*****************************************************************************************************************
	//************************ Eta Meson Combined + NLO just spectrum******************************************************
	//*****************************************************************************************************************

	TCanvas* canvasInvXSectionNLOOnlySpectraEtaSep = new TCanvas("canvasInvXSectionNLOOnlySpectraEtaSep","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlySpectraEtaSep,  0.15, 0.02, 0.03, 0.06);
	
	TPad* padXSectionNLOLegendOnlySpectraEtaSep = new TPad("padXSectionNLOLegendOnlySpectraEtaSep", "", 0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraEtaSep, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlySpectraEtaSep->Draw();

	TPad* padComparisonXSectionNLOOnlySpectraEtaSep = new TPad("padComparisonXSectionNLOOnlySpectraEtaSep", "", 0., 0., 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionNLOOnlySpectraEtaSep, 0.15, 0.02, 0.02, 0.09);
	padComparisonXSectionNLOOnlySpectraEtaSep->Draw();
		
	padComparisonXSectionNLOOnlySpectraEtaSep->cd();
	padComparisonXSectionNLOOnlySpectraEtaSep->SetLogy();		
	padComparisonXSectionNLOOnlySpectraEtaSep->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionNLOOnlySpectraEtaSep = new TH2F("histo2DInvXSectionNLOOnlySpectraEtaSep","histo2DInvXSectionNLOOnlySpectraEtaSep",1000,0.23,30.,1000,1e2,1e11);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLOOnlySpectraEtaSep, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionNLOOnlySpectraEtaSep->DrawCopy(); 
	
	graphInvCrossSectionEtaComb900GeVSysErr->Draw("2same");
	graphInvCrossSectionEtaComb900GeVStatErr->Draw("psame");
	graphInvCrossSectionEtaComb2760GeVSysErr->Draw("2same");
	graphInvCrossSectionEtaComb2760GeVStatErr->Draw("p,same");
	graphInvCrossSectionEtaComb7TeVSysErr->Draw("2same");
	graphInvCrossSectionEtaComb7TeVStatErr->Draw("psame");
	fitInvCrossSectionEta7TeV->Draw("same");
	histoFitInvCrossSectionEta900GeV->Draw("same,c");
	histoFitInvCrossSectionEta2760GeV ->Draw("same,c");
	
	graphNLOMuHalfEta7TeV->Draw("same,c");
	graphNLOMuOneEta7TeV->Draw("same,c");
	graphNLOMuTwoEta7TeV->Draw("same,c");
	
	graphNLOMuHalfEta900GeV->Draw("same,c");
	graphNLOMuOneEta900GeV->Draw("same,c");
	graphNLOMuTwoEta900GeV->Draw("same,c");
	
	graphNLOMuHalfEta2760GeV->Draw("same,c");
	graphNLOMuOneEta2760GeV->Draw("same,c");
	graphNLOMuTwoEta2760GeV->Draw("same,c");

	labelScalingEta7TeVALLEnergies->Draw();
	labelScalingEta900GeVALLEnergies->Draw();
	labelScalingEta2760GeVALLEnergies->Draw();
	DrawNormalizationErrorText(normalizationInvXOnlySpec[0]+0.1,normalizationInvXOnlySpec[1],normalizationInvXOnlySpec[2],
						  normalizationInvXOnlySpec[3],normalizationInvXOnlySpec[4],"all"); 

	//************************************************* Begin Legend ***************************************************

	padXSectionNLOLegendOnlySpectraEtaSep->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraEtaSep, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlySpectraEtaSep->SetBorderMode(-1);
	padXSectionNLOLegendOnlySpectraEtaSep->SetBorderSize(3);
	padXSectionNLOLegendOnlySpectraEtaSep->Draw();
	padXSectionNLOLegendOnlySpectraEtaSep->cd();

	//*************** first Column **********************************************************
	textSpectrumALLEnergies->Draw();
	textFitCombALLEnergies->Draw();
	textNLOMuHalfALLEnergies->Draw();
	textNLOMuOneALLEnergies->Draw();
	textNLOMuTwoALLEnergies->Draw();
	
	//*************** second Column **********************************************************
	textEta7TeVNLOEtaALLEnergies->Draw();
	textPi07TeVNLOsysALLEnergiesSep->Draw();
	boxCombinedPi07TeVALLEnergiesSep->Draw("f");
	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFit7TeVNLOALLEnergies->Draw("same");
	lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
	lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
	
	//*************** third Column **********************************************************
	textEta2760GeVNLOEtaALLEnergies->Draw();
	textPi02760GeVNLOsysALLEnergiesSep->Draw();
	boxCombinedPi02760GeVALLEnergiesSep->Draw("f");
	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi02760GeVNLOALLEnergies->Draw("same");
	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
		
	//**************** forth Column **********************************************************
	textEta900GeVNLOEtaALLEnergies->Draw();
	textPi0900GeVNLOsysALLEnergiesSep->Draw();
	boxCombinedPi0900GeVALLEnergiesSep->Draw("f");
	markerCombinedPi0900GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineFitPi0900GeVNLOALLEnergies->Draw("same");
	lineNLOPi0900GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi0900GeVMuOneALLEnergies->Draw("same");
	lineNLOPi0900GeVMuTwoALLEnergies->Draw("same");
	textArxivALLEnergies->Draw();
	
	canvasInvXSectionNLOOnlySpectraEtaSep->Update();
	canvasInvXSectionNLOOnlySpectraEtaSep->Print(Form("%s/Eta_InvXSectionNLO_OnlySpectrumSep_Paper.%s",outputDir.Data(),suffix.Data()));
	
	padComparisonXSectionNLOOnlySpectraEtaSep->cd();
	if(!thesis)DrawAliceLogoPi0WithPHOSPrelim(pictDrawingCoordinatesOnlySpectrum[0]-0.05, pictDrawingCoordinatesOnlySpectrum[1], pictDrawingCoordinatesOnlySpectrum[2], pictDrawingCoordinatesOnlySpectrum[3], pictDrawingCoordinatesOnlySpectrum[4], pictDrawingCoordinatesOnlySpectrum[5], pictDrawingCoordinatesOnlySpectrum[6], pictDrawingCoordinatesOnlySpectrum[7], pictDrawingCoordinates[8],collisionSystemCombinedReallyAll, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1200,960);
	
	canvasInvXSectionNLOOnlySpectraEtaSep->Update();
	canvasInvXSectionNLOOnlySpectraEtaSep->Print(Form("%s/Eta_InvXSectionNLO_OnlySpectrumSep.%s",outputDir.Data(),suffix.Data()));
	
//**************************************************************************************************************************************
//************************************ Inv Crosssection + NLO only ratio Eta ***************************************************************
//**************************************************************************************************************************************
	
	TCanvas* canvasInvXSectionNLOOnlyRatioEtaSep = new TCanvas("canvasInvXSectionNLOOnlyRatioEtaSep","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlyRatioEtaSep,  0.15, 0.02, 0.03, 0.06);
	canvasInvXSectionNLOOnlyRatioEtaSep->Draw();
	
	TPad* padXSectionNLOLegendOnlyRatioEtaSep = new TPad("padXSectionNLOLegendOnlyRatioEtaSep", "", 0.15, 0.8, 0.98, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioEtaSep, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlyRatioEtaSep->Draw();

	TPad* padXSectionNLOOnlyRatioEtaSepEta7TeV = new TPad("padXSectionNLOOnlyRatioEtaSepEta7TeV", "", 0., 0.565, 1., 0.8,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioEtaSepEta7TeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOOnlyRatioEtaSepEta7TeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioEtaSepEta2760GeV = new TPad("padXSectionNLOOnlyRatioEtaSepEta2760GeV", "", 0., 0.33, 1., 0.565,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioEtaSepEta2760GeV, 0.15, 0.02, 0., 0.);
	padXSectionNLOOnlyRatioEtaSepEta2760GeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioEtaSepEta900GeV = new TPad("padXSectionNLOOnlyRatioEtaSepEta900GeV", "", 0., 0., 1., 0.33,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioEtaSepEta900GeV, 0.15, 0.02, 0., 0.28);
	padXSectionNLOOnlyRatioEtaSepEta900GeV->Draw();
	
	padXSectionNLOOnlyRatioEtaSepEta7TeV->cd();
	padXSectionNLOOnlyRatioEtaSepEta7TeV->SetLogx();
	
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOPi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionNLOPi0->DrawCopy(); 
	graphRatioCombCombFitEta7TeVSys->Draw("2same");
	graphRatioCombCombFitEta7TeVStat->Draw("p,same");
	
	graphRatioCombNLOEta7TeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta7TeVMuOne->Draw("same,c");
	graphRatioCombNLOEta7TeVMuTwo->Draw("same,c");
	boxErrorSigmaPi07TeVRatio->Draw();
	labelRatioNLOEta7teVALLEnergies->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOOnlyRatioEtaSepEta2760GeV->cd();
	padXSectionNLOOnlyRatioEtaSepEta2760GeV->SetLogx();
	
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOPi0900GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.13,0.13, 0.175,0.18, 1,0.32, 512, 505);
	ratio2DInvXSectionNLOPi0900GeV->DrawCopy(); 
	
	graphRatioCombCombFitEta2760GeVSys->Draw("2same");
	graphRatioCombCombFitEta2760GeVStat->Draw("p,same");
	
	graphRatioCombNLOEta2760GeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta2760GeVMuOne->Draw("same,c");
	graphRatioCombNLOEta2760GeVMuTwo->Draw("same,c");
	boxErrorSigmaPi02760GeVRatio->Draw();
	labelRatioNLOEta2760GeV->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOOnlyRatioEtaSepEta900GeV->cd();
	padXSectionNLOOnlyRatioEtaSepEta900GeV->SetLogx();
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 0.12,0.13, 0.125,0.128, 0.9,0.45, 510, 505);
	ratio2DInvXSectionNLOEta->GetXaxis()->SetLabelOffset(-0.03);
	ratio2DInvXSectionNLOEta->DrawCopy(); 
	
	graphRatioCombCombFitEta900GeVSys->Draw("2same");
	graphRatioCombCombFitEta900GeVStat->Draw("p,same");
	
	graphRatioCombNLOEta900GeVMuHalf->Draw("same,c");
	graphRatioCombNLOEta900GeVMuOne->Draw("same,c");
	graphRatioCombNLOEta900GeVMuTwo->Draw("same,c");
	boxErrorSigmaPi0900GeVRatio->Draw();
	labelRatioNLOEta900GeVEtaALLEnergies->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);

		padXSectionNLOLegendOnlyRatioEtaSep->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioEtaSep, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlyRatioEtaSep->SetBorderMode(-1);
	padXSectionNLOLegendOnlyRatioEtaSep->SetBorderSize(3);
	padXSectionNLOLegendOnlyRatioEtaSep->Draw();
	padXSectionNLOLegendOnlyRatioEtaSep->cd();
	
	//*************** first Column **********************************************************
	textSpectrumOnlyRatioPi0->Draw();
	textNLOMuHalfOnlyRatioPi0->Draw();
	textNLOMuOneOnlyRatioPi0->Draw();
	textNLOMuTwoOnlyRatioPi0->Draw();	
	
	//*************** second Column **********************************************************
	textEta7TeVNLOOnlyRatioEta->Draw();
	textPi07TeVNLOsysOnlyRatioPi0Sep->Draw();
	boxCombinedPi07TeVOnlyRatioPi0Sep->Draw("f");
	markerCombinedPi07TeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineNLOPi07TeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi07TeVMuOneOnlyRatioPi0->Draw("same");
	lineNLOPi07TeVMuTwoOnlyRatioPi0->Draw("same");
		
	
	//*************** third Column **********************************************************
	textEta2760GeVNLOOnlyRatioEta->Draw();
	textPi02760GeVNLOsysOnlyRatioPi0Sep->Draw();
	boxCombinedPi02760GeVOnlyRatioPi0Sep->Draw("f");
	markerCombinedPi02760GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineNLOPi02760GeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi02760GeVMuOneOnlyRatioPi0->Draw("same");
	lineNLOPi02760GeVMuTwoOnlyRatioPi0->Draw("same");
	
	//**************** forth Column **********************************************************
	textEta900GeVNLOOnlyRatioEta->Draw();
	textPi0900GeVNLOsysOnlyRatioPi0Sep->Draw();
	boxCombinedPi0900GeVOnlyRatioPi0Sep->Draw("f");
	markerCombinedPi0900GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[3]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineNLOPi0900GeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi0900GeVMuOneOnlyRatioPi0->Draw("same");
	lineNLOPi0900GeVMuTwoOnlyRatioPi0->Draw("same");
	textArxivALLEnergiesRatio->Draw();

	
	canvasInvXSectionNLOOnlyRatioEtaSep->Update();
	canvasInvXSectionNLOOnlyRatioEtaSep->Print(Form("%s/Eta_InvXSectionNLO_OnlyRatioSep_Paper.%s",outputDir.Data(),suffix.Data()));
	padXSectionNLOOnlyRatioEtaSepEta900GeV->cd();
	DrawAliceLogoSimplePreliminary(pictDrawingCoordinatesOnlyRatio[4],pictDrawingCoordinatesOnlyRatio[5], pictDrawingCoordinatesOnlyRatio[6], pictDrawingCoordinatesOnlyRatio[7], 1200, 396);
	
	canvasInvXSectionNLOOnlyRatioEtaSep->Update();
	canvasInvXSectionNLOOnlyRatioEtaSep->Print(Form("%s/Eta_InvXSectionNLO_OnlyRatioSep.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	TF1* fitTsallisPi02760GeVPtYShift = ApplyYShift(graphInvCrossSectionPi0Comb2760GeVUnShifted,&graphInvCrossSectionPi0Comb2760GeVYShifted, "l", "InvXsectionY",0.3, paramGraph,0.00001);
	graphInvCrossSectionPi0Comb2760GeVYShiftedSys = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVSysErrUnShifted->Clone("YShiftedCombSys");
	graphInvCrossSectionPi0Comb2760GeVYShiftedSys= ApplyYshiftIndividualSpectra( graphInvCrossSectionPi0Comb2760GeVYShiftedSys, fitTsallisPi02760GeVPtYShift);
	graphInvCrossSectionPi0Comb2760GeVYShiftedStat = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVStatErrUnShifted->Clone("YShiftedCombStat");
	graphInvCrossSectionPi0Comb2760GeVYShiftedStat= ApplyYshiftIndividualSpectra( graphInvCrossSectionPi0Comb2760GeVYShiftedStat, fitTsallisPi02760GeVPtYShift);
		
	graphConversionXSectionPi02760GeVYShifted = (TGraphAsymmErrors*)graphConversionXSectionPi02760GeVUnShifted->Clone("YShiftedConversion");
	graphConversionXSectionPi02760GeVYShifted= ApplyYshiftIndividualSpectra( graphConversionXSectionPi02760GeVYShifted, fitTsallisPi02760GeVPtYShift);
	graphConversionXSectionPi02760GeVYShiftedSys = (TGraphAsymmErrors*)graphConversionXSectionSysPi02760GeVUnShifted->Clone("YShiftedConversionSys");
	graphConversionXSectionPi02760GeVYShiftedSys = ApplyYshiftIndividualSpectra( graphConversionXSectionPi02760GeVYShiftedSys, fitTsallisPi02760GeVPtYShift);
	graphConversionXSectionPi02760GeVYShiftedSysRAA = (TGraphAsymmErrors*)graphInvCrossSectionSysAPi02760GeV->Clone("YShiftedConversionSysRAA");
	graphConversionXSectionPi02760GeVYShiftedSysRAA = ApplyYshiftIndividualSpectra( graphConversionXSectionPi02760GeVYShiftedSysRAA, fitTsallisPi02760GeVPtYShift);
	
	graphInvCrossSectionSysAPi02760GeV= ApplyYshiftIndividualSpectra( graphInvCrossSectionSysAPi02760GeV, fitTsallisPi02760GeVPtYShift);
	graphPhosXSectionPi02760GeVYShifted = (TGraphAsymmErrors*)graphPHOSXSectionPi02760GeVUnshifted->Clone("YShiftedPHOS");
	graphPhosXSectionPi02760GeVYShifted= ApplyYshiftIndividualSpectra( graphPhosXSectionPi02760GeVYShifted, fitTsallisPi02760GeVPtYShift);
	graphPhosXSectionPi02760GeVYShiftedSys = (TGraphAsymmErrors*)graphPHOSXSectionSysPi02760GeVUnshifted->Clone("YShiftedPHOSSys");
	graphPhosXSectionPi02760GeVYShiftedSys= ApplyYshiftIndividualSpectra( graphPhosXSectionPi02760GeVYShiftedSys, fitTsallisPi02760GeVPtYShift);
	if (!conference){
		graphPhosXSectionPi02760GeVYShiftedSysRAA = (TGraphAsymmErrors*)graphSysErrRAAInvCrossSectionPi0PHOS2760GeV->Clone("YShiftedPHOSSysRAA");
		graphPhosXSectionPi02760GeVYShiftedSysRAA= ApplyYshiftIndividualSpectra( graphPhosXSectionPi02760GeVYShiftedSysRAA, fitTsallisPi02760GeVPtYShift);
	}
	cout << WriteParameterToFile(fitTsallisPi02760GeVPtYShift)<< endl;

	
	TF1*	fitInvCrossSectionOmega7TeV = FitObject("l","fitInvCrossSectionOmega7TeV","Omega",graphOmegaPhosComb7TeV,2.,17 ,paramGraph,"QNRMEX0+");
	cout << "Omega 7TeV ____________________________________________" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionOmega7TeV)<< endl;	
	

	
	TCanvas* canvasInvXSectionNLOOnlySpectra7TeVAllMesons = new TCanvas("canvasInvXSectionNLOOnlySpectra7TeVAllMesons","",200,10,1200,960);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlySpectra7TeVAllMesons,  0.15, 0.02, 0.03, 0.09);
	
	canvasInvXSectionNLOOnlySpectra7TeVAllMesons->cd();
	canvasInvXSectionNLOOnlySpectra7TeVAllMesons->SetLogy();		
	canvasInvXSectionNLOOnlySpectra7TeVAllMesons->SetLogx();		
	
	//-------------- Plotting ------------------------------------------------------
	TH2F * histo2DInvXSectionNLOOnlySpectra7TeVAllMesons = new TH2F("histo2DInvXSectionNLOOnlySpectra7TeVAllMesons","histo2DInvXSectionNLOOnlySpectra7TeVAllMesons",1000,0.23,30.,1000,4e0,5e12);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLOOnlySpectra7TeVAllMesons, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.032,0.04, 1,1.55);
	histo2DInvXSectionNLOOnlySpectra7TeVAllMesons->DrawCopy(); 
	
	graphInvCrossSectionPi0Comb7TeVSysErr->SetMarkerSize(markerSizeCommonSpectrum*1.5);
	graphInvCrossSectionPi0Comb7TeVSysErr->Draw("2same");
	graphInvCrossSectionPi0Comb7TeVStatErr->SetMarkerSize(markerSizeCommonSpectrum*1.5);
	graphInvCrossSectionPi0Comb7TeVStatErr->Draw("psame");
	fitInvCrossSectionPi0->SetLineWidth(1.5*widthCommonFit);
	fitInvCrossSectionPi0->Draw("same");
	
	graphInvCrossSectionEtaComb7TeVStatErr = ScaleGraph(graphInvCrossSectionEtaComb7TeVStatErr,1e-2);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb7TeVStatErr, markerStyleCommmonSpectrum2760GeV,markerSizeCommonSpectrum*1.5, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kFALSE);
	graphInvCrossSectionEtaComb7TeVStatErr->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionEtaComb7TeVSysErr = ScaleGraph(graphInvCrossSectionEtaComb7TeVSysErr,1e-2);
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionEtaComb7TeVSysErr, markerStyleCommmonSpectrum2760GeV,markerSizeCommonSpectrum*1.5, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi02760GeVBox);
	graphInvCrossSectionEtaComb7TeVSysErr->SetLineWidth(0);
	graphInvCrossSectionEtaComb7TeVSysErr->Draw("2same");
	graphInvCrossSectionEtaComb7TeVStatErr->Draw("psame");
	cout << "here" << endl;
	TH1D* histoFitInvCrossSectionEta7TeV_2 = (TH1D*)fitInvCrossSectionEta7TeV->GetHistogram();
	histoFitInvCrossSectionEta7TeV_2->Scale(1e-2);
	SetStyleHisto(histoFitInvCrossSectionEta7TeV_2, 1.5*widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi02760GeV);
	histoFitInvCrossSectionEta7TeV_2->Draw("same,c");
	cout << "here" << endl;
	
	graphOmegaPhos7TeV = ScaleGraph(graphOmegaPhos7TeV,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphOmegaPhos7TeV, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum*1.5, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kFALSE);
	graphOmegaPhos7TeV->SetLineWidth(widthCommonErrors);
	graphOmegaPhosSys7TeV = ScaleGraph(graphOmegaPhosSys7TeV,1e-1);
	DrawGammaSetMarkerTGraphAsym(graphOmegaPhosSys7TeV, markerStyleCommmonSpectrum900GeV,markerSizeCommonSpectrum*1.5, colorCommonSpectrumPi0900GeV, colorCommonSpectrumPi0900GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi0900GeVBox);
	graphOmegaPhosSys7TeV->SetLineWidth(0);
	graphOmegaPhosSys7TeV->Draw("2same");
	graphOmegaPhos7TeV->Draw("p,same");
	cout << "here" << endl;
	TH1D* histoFitInvCrossSectionOmega7TeV = (TH1D*)fitInvCrossSectionOmega7TeV->GetHistogram();
	histoFitInvCrossSectionOmega7TeV->Scale(1e-1);
	SetStyleHisto(histoFitInvCrossSectionOmega7TeV, 1.5*widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi0900GeV);
	histoFitInvCrossSectionOmega7TeV->Draw("same,c");
	cout << "here" << endl;
	
	
	graphNLOMuHalfPi07TeV->Draw("same,c");
	graphNLOMuOnePi07TeV->Draw("same,c");
	graphNLOMuTwoPi07TeV->Draw("same,c");
	
	graphNLOMuHalfEta7TeV= ScaleGraph(graphNLOMuHalfEta7TeV,1e-2);
	DrawGammaNLOTGraph( graphNLOMuHalfEta7TeV, widthCommonFit, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	graphNLOMuHalfEta7TeV->Draw("same,c");
	graphNLOMuOneEta7TeV= ScaleGraph(graphNLOMuOneEta7TeV,1e-2);
	DrawGammaNLOTGraph( graphNLOMuOneEta7TeV, widthCommonFit, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	graphNLOMuOneEta7TeV->Draw("same,c");
	graphNLOMuTwoEta7TeV= ScaleGraph(graphNLOMuTwoEta7TeV,1e-2);
	DrawGammaNLOTGraph( graphNLOMuTwoEta7TeV, widthCommonFit, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);
	graphNLOMuTwoEta7TeV->Draw("same,c");
	
	TLatex *labelScalingPi07TeVALLMesons = new TLatex(0.27,3E11,"#pi^{0}	");
	SetStyleTLatex( labelScalingPi07TeVALLMesons, 0.03,4,fitInvCrossSectionPi0->GetLineColor(),62,kFALSE);
	labelScalingPi07TeVALLMesons->Draw();
	TLatex *labelScalingEta7TeVALLMesons = new TLatex(0.33,8E7,"#eta #times 10^{-2}");
	SetStyleTLatex( labelScalingEta7TeVALLMesons, 0.03,4,histoFitInvCrossSectionEta7TeV_2->GetLineColor(),62,kFALSE);
	labelScalingEta7TeVALLMesons->Draw();
	TLatex *labelScalingOmega7TeVALLMesons = new TLatex(1.5,2E7,"#omega #times 10^{-1}");
	SetStyleTLatex( labelScalingOmega7TeVALLMesons, 0.03,4,histoFitInvCrossSectionOmega7TeV->GetLineColor(),62,kFALSE);
	labelScalingOmega7TeVALLMesons->Draw();

	TLegend* legendAllMesons = new TLegend(0.2,0.16,0.5,0.28);
	legendAllMesons->SetTextSize(0.03);			
	legendAllMesons->SetFillColor(0);
	legendAllMesons->SetBorderSize(0);
	legendAllMesons->AddEntry(graphInvCrossSectionPi0Comb7TeVSysErr,"#pi^{0}, PLB 717 (2012) 162-172","fp");
	legendAllMesons->AddEntry(graphOmegaPhosSys7TeV,"#omega, preliminary","pf");
	legendAllMesons->AddEntry(graphInvCrossSectionEtaComb7TeVSysErr,"#eta, PLB 717 (2012) 162-172","fp");
	legendAllMesons->Draw();	
	
	DrawNormalizationErrorText(normalizationInvXOnlySpec[0]+0.1,normalizationInvXOnlySpec[1]-0.1,normalizationInvXOnlySpec[2],
						  normalizationInvXOnlySpec[3],normalizationInvXOnlySpec[4],"7TeV"); 

	
	canvasInvXSectionNLOOnlySpectra7TeVAllMesons->Update();
	canvasInvXSectionNLOOnlySpectra7TeVAllMesons->Print(Form("%s/AllMesons7TeV_InvXSectionNLO_OnlySpectrumSep_Paper.%s",outputDir.Data(),suffix.Data()));
	
	canvasInvXSectionNLOOnlySpectra7TeVAllMesons->cd();
	if(!thesis)DrawAliceLogoAllMesonsWithPHOSPrelim(normalizationInvXOnlySpec[0]+0.1, pictDrawingCoordinatesOnlySpectrum[1], pictDrawingCoordinatesOnlySpectrum[2]+0.1, pictDrawingCoordinatesOnlySpectrum[3], normalizationInvXOnlySpec[0]+0.15, normalizationInvXOnlySpec[1]-0.05, pictDrawingCoordinatesOnlySpectrum[6]-0.02, pictDrawingCoordinatesOnlySpectrum[7], pictDrawingCoordinates[8],collisionSystem7TeV, pictDrawingOptions[1], pictDrawingOptions[2], kFALSE,1200,960);
	
	canvasInvXSectionNLOOnlySpectra7TeVAllMesons->Update();
	canvasInvXSectionNLOOnlySpectra7TeVAllMesons->Print(Form("%s/AllMesons7TeV_InvXSectionNLO_OnlySpectrumSep_Paper.%s",outputDir.Data(),suffix.Data()));

	
   TF1* fitInvCrossSectionPi02760GeVPCMStat = FitObject("l","fitInvCrossSectionPi02760GeVPCMStat","Pi0",graphConversionXSectionPi02760GeV,0.4,8.,paramGraph,"QNRMEX0+");
   cout << WriteParameterToFile(fitInvCrossSectionPi02760GeVPCMStat)<< endl;
   forOutput= WriteParameterToFile(fitInvCrossSectionPi02760GeVPCMStat);
   fileFinalResults<< forOutput.Data()<< endl;
   
   TGraphAsymmErrors* graphConversionXSectionPi02760GeVSysAndStat = CalculateCombinedSysAndStatError( graphConversionXSectionPi02760GeV, graphInvCrossSectionSysPi02760GeV);
   
   cout << "stat errors" << endl;
   graphConversionXSectionPi02760GeV->Print();
   cout << "sys errors" << endl;
   graphInvCrossSectionSysPi02760GeV->Print();
   
   cout << "combined errors" << endl;
   graphConversionXSectionPi02760GeVSysAndStat->Print();
      
      
   TF1* fitInvCrossSectionPi02760GeVPCMSysAndStat = FitObject("l","fitInvCrossSectionPi02760GeVPCMSysAndStat","Pi0",graphConversionXSectionPi02760GeVSysAndStat,0.4,8.,paramGraph,"QNRMEX0+");
   cout << WriteParameterToFile(fitInvCrossSectionPi02760GeVPCMSysAndStat)<< endl;
   forOutput= WriteParameterToFile(fitInvCrossSectionPi02760GeVPCMSysAndStat);
   fileFinalResults<< forOutput.Data()<< endl;
   
   Double_t fitresultsPCM[9];
   CalculateFitResults(fitInvCrossSectionPi02760GeVPCMStat,fitInvCrossSectionPi02760GeVPCMSysAndStat, fitresultsPCM ,"Levy",xSection2760GeV*recalcBarn);

   
   
   TF1* fitInvCrossSectionPi02760GeVPHOS = FitObject("l","fitInvCrossSectionPi02760GeVPHOS","Pi0",graphPHOSXSectionPi02760GeV,0.4,10.,paramGraph,"QNRMEX0+");
   cout << WriteParameterToFile(fitInvCrossSectionPi02760GeVPHOS)<< endl;
   forOutput= WriteParameterToFile(fitInvCrossSectionPi02760GeVPHOS);
   fileFinalResults<< forOutput.Data()<< endl;
   
   fileFinalResults.close();        
   
	//Comparison between Y and x Shift
	TCanvas* canvasComparisonXAndYShift = new TCanvas("canvasComparisonXAndYShift","",1200,900);  // gives the page size
	DrawGammaCanvasSettings( canvasComparisonXAndYShift, 0.13, 0.02, 0.02, 0.09);
	canvasComparisonXAndYShift->SetLogy();
	canvasComparisonXAndYShift->SetLogx();
	
	TH2F * histo2DInvXSectionComparisonShifts = new TH2F("histo2DInvXSectionComparisonShifts","histo2DInvXSectionComparisonShifts",1000,0.23,30.,1000,1e2,4e11);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionComparisonShifts, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionComparisonShifts->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb2760GeVUnscaled, markerStyleCommmonSpectrum,markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionPi0Comb2760GeVUnscaled->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionPi0Comb2760GeVUnscaled->Draw("p,E2same");
	
	DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb2760GeVYShifted, markerStyleCommmonSpectrum,markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeVBox, colorCommonSpectrumPi02760GeVBox, widthCommonSpectrumBoxes, kTRUE);
	graphInvCrossSectionPi0Comb2760GeVYShifted->SetLineWidth(widthCommonErrors);
	graphInvCrossSectionPi0Comb2760GeVYShifted->Draw("p,E2same");

	DrawGammaSetMarkerTGraphAsym(graphConversionXSectionPi02760GeV, markerStyleCommmonSpectrum, markerSizeCommonSpectrum,kYellow , kYellow, widthCommonSpectrumBoxes, kTRUE);
	graphConversionXSectionPi02760GeV->SetLineWidth(widthCommonErrors);
	graphConversionXSectionPi02760GeV->Draw("p,E2same");
	
	DrawGammaSetMarkerTGraphAsym(graphPHOSXSectionPi02760GeV, markerStyleCommmonSpectrum, markerSizeCommonSpectrum, kBlue, kBlue, widthCommonSpectrumBoxes, kTRUE);
	graphPHOSXSectionPi02760GeV->SetLineWidth(widthCommonErrors);
	graphPHOSXSectionPi02760GeV->Draw("p,E2same");
		
	DrawGammaSetMarkerTF1( fitInvCrossSectionPi02760GeV, styleFitCommonSpectrum, widthCommonFit, colorCommonSpectrumPi02760GeV);
	fitInvCrossSectionPi02760GeV->Draw("same");
	
	fitTsallisPi02760GeVPtYShift->SetRange(0.2,20);
	DrawGammaSetMarkerTF1( fitTsallisPi02760GeVPtYShift, styleFitCommonSpectrum+2, widthCommonFit, colorCommonSpectrumPi02760GeVBox);
	fitTsallisPi02760GeVPtYShift->Draw("same");
	
	canvasComparisonXAndYShift->Update();	
	canvasComparisonXAndYShift->SaveAs(Form("%s/ComparisonShiftingPi02760GeV.%s",outputDir.Data(),suffix.Data()));

	
	
// 	TGraphAsymmErrors* graphInvCrossSectionINEL2760GeV = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVUnscaled->Clone("graphInvCrossSectionPi0Comb2760GeVINEL");
// 	graphInvCrossSectionINEL2760GeV = ScaleGraph(graphInvCrossSectionINEL2760GeV,1/1.158);
	TFile fCombResults(Form("CombinedResultsPaper%s_%s.root",bWCorrection.Data(),dateForOutput.Data()),"RECREATE");
		graphInvCrossSectionPi0Comb7TeVUnscaled->Write("graphInvCrossSectionPi0Comb7TeV");
		graphInvCrossSectionPi0Comb7TeVStatErrUnscaled->Write("graphInvCrossSectionPi0Comb7TeVStatErr");
		graphInvCrossSectionPi0Comb7TeVSysErrUnscaled->Write("graphInvCrossSectionPi0Comb7TeVSysErr");
		graphInvCrossSectionPi0Comb7TeVUnscaled->Write("graphInvCrossSectionPi0Comb7TeV");
		graphConversionXSectionPi07TeVUnscaled->Write("graphInvCrossSectionPi0PCMStat7TeV");
		graphInvCrossSectionSysPi0Unscaled->Write("graphInvCrossSectionPi0PCMSys7TeV");
		graphPHOSXSectionPi07TeVUnscaled->Write("graphInvCrossSectionPi0PHOSStat7TeV");
		graphSysErrInvCrossSectionPi0PHOSUnscaled->Write("graphInvCrossSectionPi0PHOSSys7TeV");
		
// 		if (conference){
			cout << "writing 2.76 TeV data" << endl;
			graphInvCrossSectionPi0Comb2760GeVUnscaled->Write("graphInvCrossSectionPi0Comb2760GeV");
			graphInvCrossSectionPi0Comb2760GeVStatErrUnscaled->Write("graphInvCrossSectionPi0Comb2760GeVStatErr");
			graphInvCrossSectionPi0Comb2760GeVSysErrUnscaled->Write("graphInvCrossSectionPi0Comb2760GeVSysErr");
			graphConversionXSectionPi02760GeV->Write("graphInvCrossSectionPi0PCM2760GeVStatErr");
			graphInvCrossSectionSysPi02760GeV->Write("graphInvCrossSectionPi0PCM2760GeVSysErr");
			graphPHOSXSectionPi02760GeV->Write("graphInvCrossSectionPi0PHOS2760GeVStatErr");
			graphSysErrInvCrossSectionPi0PHOS2760GeV->Write("graphInvCrossSectionPi0PHOS2760GeVSysErr");
			fitTsallisPi02760GeVPtYShift->Write("fitInvCrossSectionPi0Comb2760GeV_YShift");
			graphInvCrossSectionPi0Comb2760GeVYShifted->Write("graphInvCrossSectionPi0Comb2760GeV_YShifted");
			graphInvCrossSectionPi0Comb2760GeVYShiftedSys->Write("graphInvCrossSectionPi0Comb2760GeVSysErr_YShifted");
			graphInvCrossSectionPi0Comb2760GeVYShiftedStat->Write("graphInvCrossSectionPi0Comb2760GeVStatErr_YShifted");
			graphConversionXSectionPi02760GeVYShifted->Write("graphInvCrossSectionPi0PCMStat2760GeV_YShifted");
			graphConversionXSectionPi02760GeVYShiftedSys->Write("graphInvCrossSectionPi0PCMSys2760GeV_YShifted");
			graphConversionXSectionPi02760GeVYShiftedSysRAA->Write("graphInvCrossSectionPi0PCMSysForRAA2760GeV_YShifted");
			graphPhosXSectionPi02760GeVYShifted->Write("graphInvCrossSectionPi0PHOSStat2760GeV_YShifted");
			graphPhosXSectionPi02760GeVYShiftedSys->Write("graphInvCrossSectionPi0PHOSSys2760GeV_YShifted");
			if (!conference){
				graphPhosXSectionPi02760GeVYShiftedSysRAA->Write("graphInvCrossSectionPi0PHOSSysForRAA2760GeV_YShifted");
			}
			graphRatioCombCombFit2760GeVSys->Write("graphRatioCombCombFit2760GeVSys");
         graphRatioCombCombFit2760GeVStat->Write("graphRatioCombCombFit2760GeVStat");
			graphRatioCombNLOPi02760GeVMuHalf->Write("graphRatioCombNLOPi02760GeVMuHalf");
			graphRatioCombNLOPi02760GeVMuOne->Write("graphRatioCombNLOPi02760GeVMuOne");
         graphRatioCombNLOPi02760GeVMuTwo->Write("graphRatioCombNLOPi02760GeVMuTwo");
         histoRatioPythia8ToFit2760GeV->Write("histoRatioPythia8ToFit2760GeV");
			histoRatioPythia8VarBinningToFit2760GeV->Write("histoRatioPythia8VarBinningToFit2760GeV");
	 		graphInvCrossSectionEtaComb2760GeVUnscaled->Write("graphInvCrossSectionEtaComb2760GeV");
			
			graphInvCrossSectionEtaComb2760GeVStatErrUnscaled->Write("graphInvCrossSectionEtaComb2760GeVStatErr_PrelimQM2011");
			graphInvCrossSectionEtaComb2760GeVSysErrUnscaled->Write("graphInvCrossSectionEtaComb2760GeVSysErr_PrelimQM2011");
			graphInvCrossSectionEtaComb900GeVUnscaled->Write("graphInvCrossSectionEtaComb900GeV_PrelimQM2011");
			graphInvCrossSectionEtaComb900GeVStatErrUnscaled->Write("graphInvCrossSectionEtaComb900GeVStatErr_PrelimQM2011");
			graphInvCrossSectionEtaComb900GeVSysErrUnscaled->Write("graphInvCrossSectionEtaComb900GeVSysErr_PrelimQM2011");
			graphRatioEtaPi0ComplErr2760GeV->Write("graphEtaToPi0CTS2760GeV");
			graphRatioEtaPi0ComplErr900GeV->Write("graphEtaToPi0CTS900GeV");

// 		}
		graphInvCrossSectionPi0Comb900GeVUnscaled->Write("graphInvCrossSectionPi0Comb900GeV");
		graphInvCrossSectionPi0Comb900GeVStatErrUnscaled->Write("graphInvCrossSectionPi0Comb900GeVStatErr");
		graphInvCrossSectionPi0Comb900GeVSysErrUnscaled->Write("graphInvCrossSectionPi0Comb900GeVSysErr");
		graphConversionXSectionPi0900GeVUnscaled->Write("graphInvCrossSectionPi0PCMStat900GeV");
		graphInvCrossSectionSysPi0900GeVUnscaled->Write("graphInvCrossSectionPi0PCMSys900GeV");
		graphPHOSXSectionPi0900GeVUnscaled->Write("graphInvCrossSectionPi0PHOSStat900GeV");
		graphSysErrInvCrossSectionPi0PHOS900GeVUnscaled->Write("graphInvCrossSectionPi0PHOSSys900GeV");

		graphInvCrossSectionEtaComb7TeVUnscaled->Write("graphInvCrossSectionEtaComb7TeV");
		graphInvCrossSectionEtaComb7TeVStatErrUnscaled->Write("graphInvCrossSectionEtaComb7TeVStatErr");
		graphInvCrossSectionEtaComb7TeVSysErrUnscaled->Write("graphInvCrossSectionEtaComb7TeVSysErr");
		graphConversionXSectionEta7TeVUnscaled->Write("graphInvCrossSectionEtaPCMStat7TeV");
		graphInvCrossSectionSysEtaUnscaled->Write("graphInvCrossSectionEtaPCMSys7TeV");
		graphPHOSXSectionEta7TeVUnscaled->Write("graphInvCrossSectionEtaPHOSStat7TeV");
		graphSysErrInvCrossSectionEtaPHOSUnscaled->Write("graphInvCrossSectionEtaPHOSSys7TeV");
	
		graphCombinedEtaToPi0->Write("graphEtaToPi0Comb7TeV");
		graphStatErrCombinedEtaToPi0->Write("graphEtaToPi0Comb7TeVStat");
		graphSysErrCombinedEtaToPi0->Write("graphEtaToPi0Comb7TeVSys");

      graphRatioCombPHOSPi0Sys->Write("graphRatioCombFitPHOSPi0Sys7TeV");
      graphRatioCombConvPi0Sys->Write("graphRatioCombFitPCMPi0Sys7TeV");
      graphRatioCombPHOSPi0->Write("graphRatioCombFitPHOSPi0Stat7TeV");
      graphRatioCombConvPi0->Write("graphRatioCombFitPCMPi0Stat7TeV");
      graphRatioCombPHOSPi0900GeVSys->Write("graphRatioCombFitPHOSPi0Sys900GeV");
      graphRatioCombConvPi0900GeVSys->Write("graphRatioCombFitPCMPi0Sys900GeV");
      graphRatioCombPHOSPi0900GeV->Write("graphRatioCombFitPHOSPi0Stat900GeV");
      graphRatioCombConvPi0900GeV->Write("graphRatioCombFitPCMPi0Stat900GeV");
      graphSysErrRatioCombPHOSPi02760GeV->Write("graphRatioCombFitPHOSPi0Sys2760GeV");
      graphSysErrRatioCombConvPi02760GeV->Write("graphRatioCombFitPCMPi0Sys2760GeV");
      histoRatioCombPHOSPi02760GeV->Write("histoRatioCombFitPHOSPi0Stat2760GeV");
      histoRatioCombConvPi02760GeV->Write("histoRatioCombFitPCMPi0Stat2760GeV");
      graphRatioCombPHOSEta7TeVSys->Write("graphRatioCombFitPHOSEtaSys7TeV");
      graphRatioCombConvEta7TeVSys->Write("graphRatioCombFitPCMEtaSys7TeV");
      graphRatioCombPHOSEta7TeV->Write("graphRatioCombFitPHOSEtaStat7TeV");
      graphRatioCombConvEta7TeV->Write("graphRatioCombFitPCMEtaStat7TeV");

      
	fCombResults.Close();
}
	
