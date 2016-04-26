/****************************************************************************************************************************
********    provided by Gamma Conversion Group, PWGGA,                                                                  *****
********    Lucia Leardini, leardini@cern.ch                                                                            *****
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

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;	

struct SysErrorConversion {
	Double_t value;
	Double_t error;
	//	TString name;
};

void CombineMesonMeasurementsPbPbLHC11h(TString meson = "Eta", 	
										TString fileNamePCM = "data_PCMResults_PbPb_2.76TeV_30June.root", 
										TString fileNameEMCalFull = "/home/admin1/leardini/alicepietapaper/pbpbrootfiles/EMCALTH1Results_PbPbJuly12015.root",  
										TString centrality = "0010",
										TString suffix = "pdf", 
										TString isMC= "", 
										TString thesisPlots = "", 
										TString bWCorrection=""){	
//root -b -l -q -x 'CombineMesonMeasurements2760GeV.C++("FinalResults/InputPCMPP_18_Feb_2014.root","FinalResults/InputPCMEMCalPP_11_Jun_2015.root","ExternalInput/EMCal/2.76TeV/EMCalResults_11Sept.root","ExternalInput/EMCal/2.76TeV/FinalCombinedXsec_pi0276TeV_25Apr2015_11a13g.root","eps","kFALSE","","X")'
//root -b -l -q -x 'CombineMesonMeasurementsPbPb2760GeV.C++("data_PCMResults_PbPb_2.76TeV_22May.root","EMCalPion010_EMCal.root")'
// root -b -l -q -x 'CombineMesonMeasurementsPbPb2760GeV.C++("data_PCMResults_PbPb_2.76TeV_15June","EMCalPion010_EMCalokbin.root")'

// root -b -l -q -x 'CombineMesonMeasurementsPbPbLHC11h.C++("Eta","data_PCMResults_PbPb_2.76TeV_2July.root","/home/admin1/leardini/alicepietapaper/pbpbrootfiles/EMCALTH1Results_PbPbJuly12015.root")'

	
	TString date = ReturnDateString();
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TString dateForOutput 			= ReturnDateStringForOutput();
	cout << dateForOutput.Data() << endl;
	//___________________________________ Declaration of files _____________________________________________
	TString collisionSystem2760GeV 			= "Pb-Pb, #sqrt{#it{s}_{NN}} = 2.76 TeV";	
	if(centrality.Contains("0010")){
		collisionSystem2760GeV = "0-10% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV"	;      
	} else if(centrality.Contains("2050")){
		collisionSystem2760GeV = "20-50% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";      
	} else if(centrality.Contains("2040")){
		collisionSystem2760GeV = "20-40% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";      
	}

	
	TString collisionSystemPbPb0010 = "0-10% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV"	;      
	TString collisionSystemPbPb2050 = "20-50% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";      
	
	TString collisionSystemPP2760GeV = "pp #sqrt{#it{s}} = 2.76 TeV";		
	TString collisionSystemPP7TeV = "pp #sqrt{#it{s}} = 7 TeV";		
	TString collisionSystemPP900GeV = "pp #sqrt{#it{s}} = 0.9 TeV";		
	TString collisionSystempPb = "p-Pb  #sqrt{s_{_{NN}}} = 5.02 TeV";
		
	
	
	TString fileNameTheory					= "ExternalInput/TheoryCompilationPP.root";	
	TFile* fileDataALICEChargedHadrons = new TFile("ExternalInputPbPb/IdentifiedCharged/PbPb_RAA_sigma_2760GeV_20120809.root");
	TFile* fileChargedPionInputPrelim2012 = new TFile("ExternalInputPbPb/IdentifiedCharged/ChargedPionSpectraPbPb_8_Nov_2013.root");
		
	TString fineNamePCMPublished 			= "CombinedResultsPbPb_18_Feb_2014.root";
// 	TString fileNamePHOSPP = "ExternalInput/PHOS/2.76TeV/LHC11a_PHOS_pi0_pp2760_noBWCorr_FDcorr_20131112_v2.root";
// 	TString fileNameChargedPionPP 			= "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_20_May_2015.root";
// 	TString fileNameChargedHadronPP			= "ExternalInput/UnidentifiedCharged/ChargedHadrinSpectraPP_20_May_2015.root";
// 	TString fileNameChargedKaons			= "";
	TString outputDir 						= Form("%s/%s/CombineMesonMeasurements2760GeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
	TString nameFinalResDat 				= Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
	cout << outputDir.Data() << endl;
	cout << fileNamePCM.Data() << endl;	
	cout << fileNameEMCalFull.Data() << endl;

	gSystem->Exec("mkdir -p "+outputDir);
 	gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPCMPublished.root", fineNamePCMPublished.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputEMCalFull.root", fileNameEMCalFull.Data(), outputDir.Data()));
// 	gSystem->Exec(Form("cp %s %s/Theory.root", fileNameTheory.Data(), outputDir.Data()));
// 	gSystem->Exec(Form("cp %s %s/ChargedPionsPrel2012.root", fileNameChargedPionPP.Data(), outputDir.Data()));
// 	gSystem->Exec(Form("cp %s %s/ChargedHadrons.root", fileNameChargedHadronPP.Data(), outputDir.Data()));
// 	gSystem->Exec(Form("cp %s %s/ChargedKaons.root", fileNameChargedKaons.Data(), outputDir.Data()));
	
	Bool_t thesis 							= kFALSE;
	if(thesisPlots.CompareTo("thesis") == 0){
		thesis 								= kTRUE;
	}
	
	TString prefix2							= "";	
	if (isMC.CompareTo("kTRUE")==0){ 
		prefix2 = "MC";
	} else {	
		prefix2 = "Data";
	}
	
	Double_t mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	Double_t mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();
	
	Double_t xSection2760GeV 		= 55.416*1e-3;
	Double_t xSection2760GeVV0AND 	= 47.73*1e-3;	
	Double_t xSection2760GeVErr 	= 3.9;
	Double_t xSection2760GeVppINEL 	= 62.8*1e9;
	Double_t recalcBarn 			= 1e12; //NLO in pbarn!!!!

	Width_t		widthLinesBoxes		= 1.4;
	Width_t		widthCommonFit		= 2;
	
	// Definition of colors, styles and markers sizes
	Color_t		colorComb			= kMagenta+2;
	Style_t		markerStyleComb		= 20;
	Size_t		markerSizeComb		= 2;
	
	Color_t 	colorCombLowPt 			= GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kFALSE);
	Color_t 	colorCombHighPt 		= GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kTRUE);
	Style_t 	markerStyleCombLowPt	= 20;
	Style_t 	markerStyleCombHighPt	= 20;	
	Size_t 		markerSizeComparison = 2;
	
	TString 	nameMeasGlobal[11] 	= {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCMOtherDataset"};
	Color_t 	colorDet[11];
	Color_t 	colorDetMC[11];
	Style_t 	markerStyleDet[11];
	Style_t 	markerStyleDetMC[11];
	Size_t 		markerSizeDet[11];
	Size_t 		markerSizeDetMC[11];

	Style_t 	styleMarkerNLOMuHalf	= 24;
	Style_t 	styleMarkerNLOMuOne		= 27;
	Style_t 	styleMarkerNLOMuTwo		= 30;
	Style_t 	styleLineNLOMuHalf		= 8;
	Style_t 	styleLineNLOMuOne		= 7;
	Style_t 	styleLineNLOMuTwo		= 4;
	Style_t 	styleLineNLOMuTwoBKK	= 3;
	Style_t 	styleLineNLOMuTwoDSS	= 6;
	Size_t		sizeMarkerNLO			= 1;
	Width_t		widthLineNLO			= 2.;

	for (Int_t i = 0; i < 11; i++){
		colorDet[i]					= GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
		colorDetMC[i]				= GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
		markerStyleDet[i]			= GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
		markerStyleDetMC[i]			= GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
		markerSizeDet[i]			= GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
		markerSizeDetMC[i]			= GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
	}	
	

	// ************************** Read data for PCM LHC11h **************************************************
	cout << "Loading PCM histos published 2010 data" << endl;
	TFile* filePCMPublished	= new TFile(fineNamePCMPublished.Data());
	
	if(centrality.Contains("0010")){
		collisionSystem2760GeV = "0-10% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV"	;      
	} else if(centrality.Contains("2050")){
		collisionSystem2760GeV = "20-50% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";      
	} else if(centrality.Contains("2040")){
		collisionSystem2760GeV = "20-40% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";      
	}

	
	// ************************** Read data for PCM LHC11h **************************************************
	cout << "Loading PCM histos" << endl;
	TFile* filePCM 									= new TFile(fileNamePCM.Data());
	
 	TH1D* histoPCMNumberOfEventsPbPb2760GeV_0010 			= (TH1D*)filePCM->Get("histoNumberOfEventsPbPb_2.76TeV0-10%");
	TH1D* histoPCMNumberOfEventsPbPb2760GeV_2050 			= (TH1D*)filePCM->Get("histoNumberOfEventsPbPb_2.76TeV20-50%");

	cout << "For the Pi0 in 0-10% " << endl;
	TDirectoryFile* directoryPCMPi0PbPb2760GeV_0010 				= (TDirectoryFile*)filePCM->Get("Pi0_PbPb_2.76TeV_0-10%"); 
		TH1D* histoPCMPi0MassPbPb2760GeV_0010 						= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("MassPi0");
		TH1D* histoPCMPi0FWHMMeVPbPb2760GeV_0010 					= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0TrueMassPbPb2760GeV_0010 					= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("TrueMassPi0");
		TH1D* histoPCMPi0TrueFWHMMeVPbPb2760GeV_0010 				= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("TrueFWHMPi0MeV");
		TH1D* histoPCMPi0AccPbPb2760GeV_0010 						= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0_Acceptance");
		TH1D* histoPCMPi0TrueEffPtPbPb2760GeV_0010 					= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0_Efficiency");
		TH1D* histoPCMPi0InvYieldPbPb2760GeV_0010 					= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("CorrectedYieldPi0");   

		TGraphAsymmErrors* graphPCMPi0RAAStat2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("RAA");
		TGraphAsymmErrors* graphPCMPi0RAASys2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("RAASys");

		TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeV_0010);
			graphPCMPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(graphPCMPi0InvYieldStatPbPb2760GeV_0010->GetN()-1);
			graphPCMPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(0);

		TGraphAsymmErrors* graphPCMPi0InvYieldSysA2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0SystErrorA");
		TGraphAsymmErrors* graphPCMPi0CorrYieldSysErrPbPb2760GeV_0010 	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0SystError"); 
		TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_0010		= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0SystError");
			graphPCMPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(graphPCMPi0InvYieldSysPbPb2760GeV_0010->GetN()-1);

		TH1D* ratioPCMPi0MassMCDiffData_0010						= (TH1D*)histoPCMPi0TrueMassPbPb2760GeV_0010->Clone("ratioPCMPi0MassMCDiffData_0010");
			ratioPCMPi0MassMCDiffData_0010->Sumw2();
			ratioPCMPi0MassMCDiffData_0010->Add(histoPCMPi0MassPbPb2760GeV_0010,-1);
			ratioPCMPi0MassMCDiffData_0010->Divide(ratioPCMPi0MassMCDiffData_0010, histoPCMPi0MassPbPb2760GeV_0010,1.,1.,"");
			ratioPCMPi0MassMCDiffData_0010->Scale(1./mesonMassExpectPi0);

 		TH1D *histoPCMPi0AccTimesEffPbPb2760GeV_0010				= (TH1D*)histoPCMPi0TrueEffPtPbPb2760GeV_0010->Clone("histoPCMPi0AccTimesEffPbPb2760GeV_0010");
			histoPCMPi0AccTimesEffPbPb2760GeV_0010->Sumw2();
			histoPCMPi0AccTimesEffPbPb2760GeV_0010->Multiply(histoPCMPi0AccPbPb2760GeV_0010);
		
		
	cout << "For the Pi0 in 20-50% " << endl;
	TDirectory* directoryPCMPi0PbPb2760GeV_2050 					= (TDirectory*)filePCM->Get("Pi0_PbPb_2.76TeV_20-50%"); 
		TH1D* histoPCMPi0MassPbPb2760GeV_2050 						= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("MassPi0");
		TH1D* histoPCMPi0FWHMMeVPbPb2760GeV_2050 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0TrueMassPbPb2760GeV_2050 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("TrueMassPi0");
		TH1D* histoPCMPi0TrueFWHMMeVPbPb2760GeV_2050 				= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("TrueFWHMPi0MeV");
		TH1D* histoPCMPi0AccPbPb2760GeV_2050 						= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0_Acceptance");
		TH1D* histoPCMPi0TrueEffPtPbPb2760GeV_2050 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0_Efficiency");
	    TH1D* histoPCMPi0InvYieldPbPb2760GeV_2050 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("CorrectedYieldPi0");   

		TGraphAsymmErrors* graphPCMPi0RAAStat2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("RAA");
		TGraphAsymmErrors* graphPCMPi0RAASys2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("RAASys");

		TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeV_2050);
			graphPCMPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(graphPCMPi0InvYieldStatPbPb2760GeV_2050->GetN()-1);
			graphPCMPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(0);
				
		TGraphAsymmErrors* graphPCMPi0InvYieldSysA2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0SystErrorA");
        TGraphAsymmErrors* graphPCMPi0CorrYieldSysErrPbPb2760GeV_2050 	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0SystError"); 
		TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_2050		= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0SystError");
				graphPCMPi0InvYieldSysPbPb2760GeV_2050->RemovePoint(graphPCMPi0InvYieldSysPbPb2760GeV_2050->GetN()-1);
		
		TH1D* ratioPCMPi0MassMCDiffData_2050						= (TH1D*)histoPCMPi0TrueMassPbPb2760GeV_2050->Clone("ratioPCMPi0MassMCDiffData_2050");
			ratioPCMPi0MassMCDiffData_2050->Sumw2();
			ratioPCMPi0MassMCDiffData_2050->Add(histoPCMPi0MassPbPb2760GeV_2050,-1);
			ratioPCMPi0MassMCDiffData_2050->Divide(ratioPCMPi0MassMCDiffData_2050, histoPCMPi0MassPbPb2760GeV_2050,1.,1.,"");
			ratioPCMPi0MassMCDiffData_2050->Scale(1./mesonMassExpectPi0);
		
 		TH1D* histoPCMPi0AccTimesEffPbPb2760GeV_2050				= (TH1D*)histoPCMPi0TrueEffPtPbPb2760GeV_2050->Clone("histoPCMPi0AccTimesEffPbPb2760GeV_2050");
			histoPCMPi0AccTimesEffPbPb2760GeV_2050->Sumw2();
			histoPCMPi0AccTimesEffPbPb2760GeV_2050->Multiply(histoPCMPi0AccPbPb2760GeV_2050);

			
		TGraphAsymmErrors* graphPCMPi0RCPStat2760GeV	= (TGraphAsymmErrors*)filePCM->Get("Pi0RCP");
		TGraphAsymmErrors* graphPCMPi0RCPSys2760GeV	= (TGraphAsymmErrors*)filePCM->Get("Pi0RCPsys");
		

	cout << "For the Eta in 0-10% " << endl;
	TDirectory* directoryPCMEtaPbPb2760GeV_0010 					= (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_0-10%"); 
		TH1D* histoPCMEtaMassPbPb2760GeV_0010 						= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("MassEta");
		TH1D* histoPCMEtaFWHMMeVPbPb2760GeV_0011 					= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("FWHMEtaMeV");
		TH1D* histoPCMEtaTrueMassPbPb2760GeV_0010 					= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("TrueMassEta");
		TH1D* histoPCMEtaTrueFWHMMeVPbPb2760GeV_0010 				= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("TrueFWHMEtaMeV");
		TH1D* histoPCMEtaAccPbPb2760GeV_0010 						= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("Eta_Acceptance");
		TH1D* histoPCMEtaTrueEffPtPbPb2760GeV_0010 					= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("Eta_Efficiency");
	    TH1D* histoPCMEtaInvYieldPbPb2760GeV_0010 					= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("CorrectedYieldEta");   
		
		TGraphAsymmErrors* graphPCMEtaRAAStat2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("RAA");
		TGraphAsymmErrors* graphPCMEtaRAASys2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("RAASys");

		TH1D* histoPCMEtatoPi0Stat2760GeV_0010 					= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("EtatoPi0Ratio"); 
		TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_0010	= new TGraphAsymmErrors(histoPCMEtatoPi0Stat2760GeV_0010);
		TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtatoPi0RatioSys");

		TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeV_0010);
// 			graphPCMEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(graphPCMEtaInvYieldStatPbPb2760GeV_0010->GetN()-1);
// 			graphPCMEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(0);

		TGraphAsymmErrors* graphPCMEtaInvYieldSysA2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaSystErrorA");
        TGraphAsymmErrors* graphPCMEtaCorrYieldSysErrPbPb2760GeV_0010 	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaSystError"); 
// 		graphPCMEtaCorrYieldSysErrPbPb2760GeV_0010->RemovePoint(0);
		TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_0010		= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaSystError");
// 			graphPCMEtaInvYieldSysPbPb2760GeV_0010->RemovePoint(graphPCMEtaInvYieldSysPbPb2760GeV_0010->GetN()-1);

		TH1D* ratioPCMEtaMassMCDiffData_0010						= (TH1D*)histoPCMEtaTrueMassPbPb2760GeV_0010->Clone("ratioPCMEtaMassMCDiffData_0010");
			ratioPCMEtaMassMCDiffData_0010->Sumw2();
			ratioPCMEtaMassMCDiffData_0010->Add(histoPCMEtaMassPbPb2760GeV_0010,-1);
			ratioPCMEtaMassMCDiffData_0010->Divide(ratioPCMEtaMassMCDiffData_0010, histoPCMEtaMassPbPb2760GeV_0010,1.,1.,"");
			ratioPCMEtaMassMCDiffData_0010->Scale(1./mesonMassExpectEta);
		
 		TH1D* histoPCMEtaAccTimesEffPbPb2760GeV_0010				= (TH1D*)histoPCMEtaTrueEffPtPbPb2760GeV_0010->Clone("histoPCMEtaAccTimesEffPbPb2760GeV_0010");
			histoPCMEtaAccTimesEffPbPb2760GeV_0010->Sumw2();
			histoPCMEtaAccTimesEffPbPb2760GeV_0010->Multiply(histoPCMEtaAccPbPb2760GeV_0010);
		
		
	cout << "For the Eta in 20-50% " << endl;
	TDirectory* directoryPCMEtaPbPb2760GeV_2050 					= (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_20-50%"); 
		TH1D* histoPCMEtaMassPbPb2760GeV_2050 						= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("MassEta");
		TH1D* histoPCMEtaFWHMMeVPbPb2760GeV_2050 					= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("FWHMEtaMeV");
		TH1D* histoPCMEtaTrueMassPbPb2760GeV_2050 					= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("TrueMassEta");
		TH1D* histoPCMEtaTrueFWHMMeVPbPb2760GeV_2050 				= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("TrueFWHMEtaMeV");
		TH1D* histoPCMEtaAccPbPb2760GeV_2050 						= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("Eta_Acceptance");
		TH1D* histoPCMEtaTrueEffPtPbPb2760GeV_2050 					= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("Eta_Efficiency");
	    TH1D* histoPCMEtaInvYieldPbPb2760GeV_2050 					= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("CorrectedYieldEta"); 

		TGraphAsymmErrors* graphPCMEtaRAAStat2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("RAA");
		TGraphAsymmErrors* graphPCMEtaRAASys2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("RAASys");

		TH1D* histoPCMEtatoPi0Stat2760GeV_2050 					= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("EtatoPi0Ratio"); 
		TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_2050	= new TGraphAsymmErrors(histoPCMEtatoPi0Stat2760GeV_2050);
		TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtatoPi0RatioSys");

		TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeV_2050);
// 			graphPCMEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(graphPCMEtaInvYieldStatPbPb2760GeV_2050->GetN()-1);
// 			graphPCMEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(0);
		TGraphAsymmErrors* graphPCMEtaInvYieldSysA2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaSystErrorA");
        TGraphAsymmErrors* graphPCMEtaCorrYieldSysErrPbPb2760GeV_2050 	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaSystError"); 
// 		graphPCMEtaCorrYieldSysErrPbPb2760GeV_2050->RemovePoint(0);
		TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_2050		= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaSystError");
// 			graphPCMEtaInvYieldSysPbPb2760GeV_2050->RemovePoint(graphPCMEtaInvYieldSysPbPb2760GeV_2050->GetN()-1);
		
		TH1D* ratioPCMEtaMassMCDiffData_2050						= (TH1D*)histoPCMEtaTrueMassPbPb2760GeV_2050->Clone("ratioPCMEtaMassMCDiffData_2050");
			ratioPCMEtaMassMCDiffData_2050->Sumw2();
			ratioPCMEtaMassMCDiffData_2050->Add(histoPCMEtaMassPbPb2760GeV_2050,-1);
			ratioPCMEtaMassMCDiffData_2050->Divide(ratioPCMEtaMassMCDiffData_2050, histoPCMEtaMassPbPb2760GeV_2050,1.,1.,"");
			ratioPCMEtaMassMCDiffData_2050->Scale(1./mesonMassExpectEta);
		
 		TH1D* histoPCMEtaAccTimesEffPbPb2760GeV_2050				= (TH1D*)histoPCMEtaTrueEffPtPbPb2760GeV_2050->Clone("histoPCMEtaAccTimesEffPbPb2760GeV_2050");
			histoPCMEtaAccTimesEffPbPb2760GeV_2050->Sumw2();
			histoPCMEtaAccTimesEffPbPb2760GeV_2050->Multiply(histoPCMEtaAccPbPb2760GeV_2050);	
			
			
		TGraphAsymmErrors* graphPCMEtaRCPStat2760GeV	= (TGraphAsymmErrors*)filePCM->Get("EtaRCP");
		TGraphAsymmErrors* graphPCMEtaRCPSys2760GeV	= (TGraphAsymmErrors*)filePCM->Get("EtaRCPsys");
			
		
// 	Int_t nEvtPCM2760GeV 							= histoPCMNumberOfEventsPbPb2760GeV_0010->GetBinContent(1);
// 	cout << "here" << endl;
	histoPCMPi0MassPbPb2760GeV_0010->Scale(1000.);
	histoPCMPi0TrueMassPbPb2760GeV_0010->Scale(1000.);
	histoPCMEtaMassPbPb2760GeV_0010->Scale(1000.);
	histoPCMEtaTrueMassPbPb2760GeV_0010->Scale(1000.);
	
	histoPCMPi0MassPbPb2760GeV_2050->Scale(1000.);
	histoPCMPi0TrueMassPbPb2760GeV_2050->Scale(1000.);
	histoPCMEtaMassPbPb2760GeV_2050->Scale(1000.);
	histoPCMEtaTrueMassPbPb2760GeV_2050->Scale(1000.);

	cout << "here" << endl;
	

	// ************************** Read data for EMCal ****************************************************		
	cout << "Loading EMCal histos" << endl;	
	TFile* fileEMCal								= new TFile(fileNameEMCalFull.Data());
	TFile* fileEMCalbis								= new TFile("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/EMCALyields_1July.root");
//for EMCal file with Graph only
// 	cout << "Pi0 0-10%" << endl;
// 		TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeV_0010 		= (TGraphAsymmErrors*)fileEMCal->Get("EMCalStatPions010");
// 		TGraphAsymmErrors* graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010 	= (TGraphAsymmErrors*)fileEMCal->Get("EMCalSysPions010");
// 		TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeV_0010 	= (TGraphAsymmErrors*)fileEMCal->Get("EMCalSysPions010");

// 		TGraphAsymmErrors* graphEMCalEtatoPi0Stat2760GeV_0010 	= (TGraphAsymmErrors*)fileEMCalbis->Get("EMCalStatRatio010");
		TGraphAsymmErrors* graphEMCalEtatoPi0Sys2760GeV_0010 	= (TGraphAsymmErrors*)fileEMCalbis->Get("EMCalSysRatio010");

// 	cout << "Pi0 20-50%" << endl;		
// 		TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeV_2050 		= (TGraphAsymmErrors*)fileEMCal->Get("EMCalStatPions2050");
// 		TGraphAsymmErrors* graphEMCalPi0CorrYieldSysErrPbPb2760GeV_2050 	= (TGraphAsymmErrors*)fileEMCal->Get("EMCalSysPions2050");
// 		TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeV_2050 	= (TGraphAsymmErrors*)fileEMCal->Get("EMCalSysPions2050");

// 		TGraphAsymmErrors* graphEMCalEtatoPi0Stat2760GeV_2050 	= (TGraphAsymmErrors*)fileEMCalbis->Get("EMCalStatRatio2050");
		TGraphAsymmErrors* graphEMCalEtatoPi0Sys2760GeV_2050 	= (TGraphAsymmErrors*)fileEMCalbis->Get("EMCalSysRatio2050");
		
// 	cout << "Eta 0-10%" << endl;
// 		TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeV_0010 		= (TGraphAsymmErrors*)fileEMCal->Get("EMCalStatEta010");
// 		TGraphAsymmErrors* graphEMCalEtaCorrYieldSysErrPbPb2760GeV_0010 	= (TGraphAsymmErrors*)fileEMCal->Get("EMCalSysEta010");
// 		TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeV_0010 	= (TGraphAsymmErrors*)fileEMCal->Get("EMCalSysEta010");
// 
// 	cout << "Eta 20-50%" << endl;		
// 		TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeV_2050 		= (TGraphAsymmErrors*)fileEMCal->Get("EMCalStatEta2050");
// 		TGraphAsymmErrors* graphEMCalEtaCorrYieldSysErrPbPb2760GeV_2050 	= (TGraphAsymmErrors*)fileEMCal->Get("EMCalSysEta2050");
// 		TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeV_2050 	= (TGraphAsymmErrors*)fileEMCal->Get("EMCalSysEta2050");

	
//for EMCal file with histos, used naming below
	TDirectory* directoryEMCalPi0PbPb2760GeV 				= (TDirectory*)fileEMCal->Get("Pi02.76TeV_PbPb");
	cout << "Pi0 0-10%" << endl;
		TH1D* histoEMCalPi0InvYieldPbPb2760GeV_0010			= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("InvYieldPbPbStatErrPi_0010");
		TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeV_0010 		= new TGraphAsymmErrors(histoEMCalPi0InvYieldPbPb2760GeV_0010);
// 			graphEMCalPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(0);

		TH1D* histoEMCalPi0InvYieldSysPbPb2760GeV_0010 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("InvYieldPbPbSysErrPi_0010");
// 		TGraphAsymmErrors* graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010 	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("InvYieldPbPbSysErrPi_0010"); 
		TGraphAsymmErrors* graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_0010); 
		TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_0010);
// 			graphEMCalPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(0);

		TH1D*	histoEMCalEtatoPi0StatPbPb2760GeV_0010 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("StatErrEtatoPi0Ratio_0010");
		TGraphAsymmErrors* graphEMCalEtatoPi0Stat2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalEtatoPi0StatPbPb2760GeV_0010);
// 		TH1D*	histoEMCalEtatoPi0SysPbPb2760GeV_0010 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("SysErrEtatoPi0Ratio_0010");
// 		TGraphAsymmErrors* graphEMCalEtatoPi0Sys2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalEtatoPi0SysPbPb2760GeV_0010);
		
	cout << "Pi0 20-50%" << endl;		
		TH1D* histoEMCalPi0InvYieldPbPb2760GeV_2050			= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("InvYieldPbPbStatErrPi_2050");
		TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeV_2050 		= new TGraphAsymmErrors(histoEMCalPi0InvYieldPbPb2760GeV_2050);
// 			graphEMCalPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(0);

		TH1D*	histoEMCalPi0InvYieldSysPbPb2760GeV_2050 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("InvYieldPbPbSysErrPi_2050");
// 		TGraphAsymmErrors* graphEMCalPi0CorrYieldSysErrPbPb2760GeV_2050 	= (TGraphAsymmErrors*)directoryEMCalPi0PbPb2760GeV->Get("InvYieldPbPbStatErrPi_0010"); 
		TGraphAsymmErrors* graphEMCalPi0CorrYieldSysErrPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_2050);
		TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_2050);
// 			graphEMCalPi0InvYieldSysPbPb2760GeV_2050->RemovePoint(0);

		TH1D*	histoEMCalEtatoPi0StatPbPb2760GeV_2050 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("StatErrEtatoPi0Ratio_2050");
		TGraphAsymmErrors* graphEMCalEtatoPi0Stat2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalEtatoPi0StatPbPb2760GeV_2050);
// 		TH1D*	histoEMCalEtatoPi0SysPbPb2760GeV_2050 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("SysErrEtatoPi0Ratio_2050");
// 		TGraphAsymmErrors* graphEMCalEtatoPi0Sys2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalEtatoPi0SysPbPb2760GeV_2050);

		TH1D*	histoEMCalPi0RCPStatPbPb2760GeV 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("StatErrRcpPi");
		TGraphAsymmErrors* graphEMCalPi0RCPStat2760GeV 	= new TGraphAsymmErrors(histoEMCalPi0RCPStatPbPb2760GeV);
		TH1D*	histoEMCalPi0RCPSysPbPb2760GeV 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("SysErrRcpPi");
		TGraphAsymmErrors* graphEMCalPi0RCPSys2760GeV 	= new TGraphAsymmErrors(histoEMCalPi0RCPSysPbPb2760GeV);

		
		
		
	TDirectory* directoryEMCalEtaPbPb2760GeV 				= (TDirectory*)fileEMCal->Get("Eta2.76TeV_PbPb");
	cout << "Eta 0-10%" << endl;
		TH1D* histoEMCalEtaInvYieldPbPb2760GeV_0010			= (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbStatErrEta_0010");
		TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeV_0010 		= new TGraphAsymmErrors(histoEMCalEtaInvYieldPbPb2760GeV_0010);
// 			graphEMCalEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(0);

		TH1D*	histoEMCalEtaInvYieldSysPbPb2760GeV_0010 	= (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbSysErrEta_0010");
// 		TGraphAsymmErrors* graphEMCalEtaCorrYieldSysErrPbPb2760GeV_0010 	= (TGraphAsymmErrors*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbStatErrEta_0010"); 
		TGraphAsymmErrors* graphEMCalEtaCorrYieldSysErrPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_0010);
		TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_0010);
// 			graphEMCalEtaInvYieldSysPbPb2760GeV_0010->RemovePoint(0);

	cout << "Eta 20-50%" << endl;		
		TH1D* histoEMCalEtaInvYieldPbPb2760GeV_2050			= (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbStatErrEta_2050");
		TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeV_2050 		= new TGraphAsymmErrors(histoEMCalEtaInvYieldPbPb2760GeV_2050);
// 			graphEMCalEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(0);

		TH1D*	histoEMCalEtaInvYieldSysPbPb2760GeV_2050 	= (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbSysErrEta_2050");
// 		TGraphAsymmErrors* graphEMCalEtaCorrYieldSysErrPbPb2760GeV_2050 	= (TGraphAsymmErrors*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbSysErrEta_2050"); 
		TGraphAsymmErrors* graphEMCalEtaCorrYieldSysErrPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_2050);
		TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_2050);
// 			graphEMCalEtaInvYieldSysPbPb2760GeV_2050->RemovePoint(0);
		
		TH1D*	histoEMCalEtaRCPStatPbPb2760GeV 	= (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("StatErrRcpEta");
		TGraphAsymmErrors* graphEMCalEtaRCPStat2760GeV 	= new TGraphAsymmErrors(histoEMCalEtaRCPStatPbPb2760GeV);
		TH1D*	histoEMCalEtaRCPSysPbPb2760GeV 	= (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("SysErrRcpEta");
		TGraphAsymmErrors* graphEMCalEtaRCPSys2760GeV 	= new TGraphAsymmErrors(histoEMCalEtaRCPSysPbPb2760GeV);

		
	TH1D *histoPCMPi0InvYieldPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphPCMPi0CorrYieldSysErrPbPb2760GeV = NULL;
	TH1D *histoEMCalPi0InvYieldPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphEMCalPi0CorrYieldSysErrPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphPCMPi0InvYieldStatPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphPCMPi0InvYieldSysPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphEMCalPi0InvYieldStatPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphEMCalPi0InvYieldSysPbPb2760GeV = NULL;
	
	TH1D *histoPCMEtaInvYieldPbPb2760GeV = NULL;
	TH1D *histoPCMEtatoPi0Stat2760GeV = NULL;
	TGraphAsymmErrors *graphPCMEtaCorrYieldSysErrPbPb2760GeV = NULL;
	TH1D *histoEMCalEtaInvYieldPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphEMCalEtaCorrYieldSysErrPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphPCMEtaInvYieldStatPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphPCMEtaInvYieldSysPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphEMCalEtaInvYieldStatPbPb2760GeV = NULL;
	TGraphAsymmErrors *graphEMCalEtaInvYieldSysPbPb2760GeV = NULL;
	
	TH1D *histoEMCalEtatoPi0StatPbPb2760GeV = NULL;
	TGraphAsymmErrors *	graphPCMPi0RAAStat2760GeV = NULL;
	TGraphAsymmErrors *	graphPCMPi0RAASys2760GeV = NULL;
	TGraphAsymmErrors *	graphEMCalPi0RAAStat2760GeV = NULL;
	TGraphAsymmErrors *	graphEMCalPi0RAASys2760GeV  = NULL;
	
	TGraphAsymmErrors *	graphEMCalEtatoPi0Stat2760GeV  = NULL;
	TGraphAsymmErrors *	graphEMCalEtatoPi0Sys2760GeV  = NULL;
	TGraphAsymmErrors *	graphPCMEtaRAAStat2760GeV  = NULL;
	TGraphAsymmErrors *	graphPCMEtaRAASys2760GeV  = NULL;

	TGraphAsymmErrors *	graphPCMEtatoPi0Stat2760GeV  = NULL;
	TGraphAsymmErrors *	graphPCMEtatoPi0Sys2760GeV  = NULL;
	TGraphAsymmErrors *	graphEMCalEtaRAAStat2760GeV  = NULL;
	TGraphAsymmErrors *	graphEMCalEtaRAASys2760GeV  = NULL;

	
	if (centrality.CompareTo("0010")==0){
		
		//=============================================================== PCM 
		histoPCMPi0InvYieldPbPb2760GeV = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_0010->Clone("histoPCMPi0InvYieldPbPb2760GeV");
		graphPCMPi0CorrYieldSysErrPbPb2760GeV = (TGraphAsymmErrors*)graphPCMPi0CorrYieldSysErrPbPb2760GeV_0010->Clone("graphPCMPi0CorrYieldSysErrPbPb2760GeV");
		graphPCMPi0InvYieldStatPbPb2760GeV = (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_0010->Clone("graphPCMPi0InvYieldStatPbPb2760GeV");
		graphPCMPi0InvYieldSysPbPb2760GeV = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0010->Clone("graphPCMPi0InvYieldSysPbPb2760GeV");

		histoPCMEtaInvYieldPbPb2760GeV = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_0010->Clone("histoPCMEtaInvYieldPbPb2760GeV");
		graphPCMEtaCorrYieldSysErrPbPb2760GeV = (TGraphAsymmErrors*)graphPCMEtaCorrYieldSysErrPbPb2760GeV_0010->Clone("graphPCMEtaCorrYieldSysErrPbPb2760GeV");
		graphPCMEtaInvYieldStatPbPb2760GeV = (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_0010->Clone("graphPCMEtaInvYieldStatPbPb2760GeV");
		graphPCMEtaInvYieldSysPbPb2760GeV = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0010->Clone("graphPCMEtaInvYieldSysPbPb2760GeV");

		graphPCMPi0RAAStat2760GeV = (TGraphAsymmErrors*)graphPCMPi0RAAStat2760GeV_0010->Clone("graphPCMPi0RAAStat2760GeV");
		graphPCMPi0RAASys2760GeV = (TGraphAsymmErrors*)graphPCMPi0RAASys2760GeV_0010->Clone("graphPCMPi0RAASys2760GeV");
		graphPCMEtaRAAStat2760GeV = (TGraphAsymmErrors*)graphPCMEtaRAAStat2760GeV_0010->Clone("graphPCMEtaRAAStat2760GeV");
		graphPCMEtaRAASys2760GeV = (TGraphAsymmErrors*)graphPCMEtaRAAStat2760GeV_0010->Clone("graphPCMEtaRAASys2760GeV");

		histoPCMEtatoPi0Stat2760GeV = (TH1D*)histoPCMEtatoPi0Stat2760GeV_0010->Clone("histoPCMEtatoPi0Stat2760GeV");

		graphPCMEtatoPi0Stat2760GeV = (TGraphAsymmErrors*)graphPCMEtatoPi0Stat2760GeV_0010->Clone("graphPCMEtatoPi0Stat2760GeV");
		graphPCMEtatoPi0Sys2760GeV = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_0010->Clone("graphPCMEtatoPi0Sys2760GeV");



		//=============================================================== EMCal
		histoEMCalPi0InvYieldPbPb2760GeV = (TH1D*)histoEMCalPi0InvYieldPbPb2760GeV_0010->Clone("histoEMCalPi0InvYieldPbPb2760GeV");
		graphEMCalPi0CorrYieldSysErrPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010->Clone("graphEMCalPi0CorrYieldSysErrPbPb2760GeV");
		graphEMCalPi0InvYieldStatPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Clone("graphEMCalPi0InvYieldStatPbPb2760GeV");
		graphEMCalPi0InvYieldSysPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Clone("graphEMCalPi0InvYieldSysPbPb2760GeV");

		histoEMCalEtaInvYieldPbPb2760GeV = (TH1D*)histoEMCalEtaInvYieldPbPb2760GeV_0010->Clone("histoEMCalEtaInvYieldPbPb2760GeV");
		graphEMCalEtaCorrYieldSysErrPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalEtaCorrYieldSysErrPbPb2760GeV_0010->Clone("graphEMCalEtaCorrYieldSysErrPbPb2760GeV");
		graphEMCalEtaInvYieldStatPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_0010->Clone("graphEMCalEtaInvYieldStatPbPb2760GeV");
		graphEMCalEtaInvYieldSysPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Clone("graphEMCalEtaInvYieldSysPbPb2760GeV");

// 		graphEMCalPi0RAAStat2760GeV = (TGraphAsymmErrors*)graphEMCalPi0RAAStat2760GeV_0010->Clone("graphEMCalPi0RAAStat2760GeV");
// 		graphEMCalPi0RAASys2760GeV = (TGraphAsymmErrors*)graphEMCalPi0RAASys2760GeV_0010->Clone("graphEMCalPi0RAASys2760GeV");
// 		graphEMCalEtaRAAStat2760GeV = (TGraphAsymmErrors*)graphEMCalPi0RAAStat2760GeV_0010->Clone("graphEMCalPi0RAAStat2760GeV");
// 		graphEMCalEtaRAASys2760GeV = (TGraphAsymmErrors*)graphEMCalPi0RAASys2760GeV_0010->Clone("graphEMCalPi0RAASys2760GeV");

		histoEMCalEtatoPi0StatPbPb2760GeV = (TH1D*)histoEMCalEtatoPi0StatPbPb2760GeV_0010->Clone("histoEMCalEtatoPi0StatPbPb2760GeV");
		graphEMCalEtatoPi0Stat2760GeV = new TGraphAsymmErrors(histoEMCalEtatoPi0StatPbPb2760GeV);
//		graphEMCalEtatoPi0Stat2760GeV = (TGraphAsymmErrors*)graphEMCalEtatoPi0Stat2760GeV_0010->Clone("graphEMCalEtatoPi0Stat2760GeV");
		graphEMCalEtatoPi0Sys2760GeV = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV_0010->Clone("graphEMCalEtatoPi0Sys2760GeV");		
		
	} else if (centrality.CompareTo("2050")==0){
		
		
		//=============================================================== PCM
		histoPCMPi0InvYieldPbPb2760GeV = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_2050->Clone("histoPCMPi0InvYieldPbPb2760GeV");
		graphPCMPi0CorrYieldSysErrPbPb2760GeV = (TGraphAsymmErrors*)graphPCMPi0CorrYieldSysErrPbPb2760GeV_2050->Clone("graphPCMPi0CorrYieldSysErrPbPb2760GeV");
		graphPCMPi0InvYieldStatPbPb2760GeV = (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_2050->Clone("graphPCMPi0InvYieldStatPbPb2760GeV");
		graphPCMPi0InvYieldSysPbPb2760GeV = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2050->Clone("graphPCMPi0InvYieldSysPbPb2760GeV");
		
		histoPCMEtaInvYieldPbPb2760GeV = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_2050->Clone("histoPCMEtaInvYieldPbPb2760GeV");
		graphPCMEtaCorrYieldSysErrPbPb2760GeV = (TGraphAsymmErrors*)graphPCMEtaCorrYieldSysErrPbPb2760GeV_2050->Clone("graphPCMEtaCorrYieldSysErrPbPb2760GeV");
		graphPCMEtaInvYieldStatPbPb2760GeV = (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_2050->Clone("graphPCMEtaInvYieldStatPbPb2760GeV");
		graphPCMEtaInvYieldSysPbPb2760GeV = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2050->Clone("graphPCMEtaInvYieldSysPbPb2760GeV");

		graphPCMPi0RAAStat2760GeV = (TGraphAsymmErrors*)graphPCMPi0RAAStat2760GeV_2050->Clone("graphPCMPi0RAAStat2760GeV");
		graphPCMPi0RAASys2760GeV = (TGraphAsymmErrors*)graphPCMPi0RAASys2760GeV_2050->Clone("graphPCMPi0RAASys2760GeV");		
		graphPCMEtaRAAStat2760GeV = (TGraphAsymmErrors*)graphPCMEtaRAAStat2760GeV_2050->Clone("graphPCMEtaRAAStat2760GeV");
		graphPCMEtaRAASys2760GeV = (TGraphAsymmErrors*)graphPCMEtaRAAStat2760GeV_2050->Clone("graphPCMEtaRAASys2760GeV");
		
		histoPCMEtatoPi0Stat2760GeV = (TH1D*)histoPCMEtatoPi0Stat2760GeV_2050->Clone("histoPCMEtatoPi0Stat2760GeV");

 		graphPCMEtatoPi0Stat2760GeV = (TGraphAsymmErrors*)graphPCMEtatoPi0Stat2760GeV_2050->Clone("graphPCMEtatoPi0Stat2760GeV");
		graphPCMEtatoPi0Sys2760GeV = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_2050->Clone("graphPCMEtatoPi0Sys2760GeV");


		//=============================================================== EMCal
		histoEMCalPi0InvYieldPbPb2760GeV = (TH1D*)histoEMCalPi0InvYieldPbPb2760GeV_2050->Clone("histoEMCalPi0InvYieldPbPb2760GeV");
		graphEMCalPi0CorrYieldSysErrPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalPi0CorrYieldSysErrPbPb2760GeV_2050->Clone("graphEMCalPi0CorrYieldSysErrPbPb2760GeV");
		graphEMCalPi0InvYieldStatPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_2050->Clone("graphEMCalPi0InvYieldStatPbPb2760GeV");
		graphEMCalPi0InvYieldSysPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Clone("graphEMCalPi0InvYieldSysPbPb2760GeV");
		
		histoEMCalEtaInvYieldPbPb2760GeV = (TH1D*)histoEMCalEtaInvYieldPbPb2760GeV_2050->Clone("histoEMCalEtaInvYieldPbPb2760GeV");
		graphEMCalEtaCorrYieldSysErrPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalEtaCorrYieldSysErrPbPb2760GeV_2050->Clone("graphEMCalEtaCorrYieldSysErrPbPb2760GeV");
		graphEMCalEtaInvYieldStatPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Clone("graphEMCalEtaInvYieldStatPbPb2760GeV");
		graphEMCalEtaInvYieldSysPbPb2760GeV = (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Clone("graphEMCalEtaInvYieldSysPbPb2760GeV");

// 		graphEMCalPi0RAAStat2760GeV = (TGraphAsymmErrors*)graphEMCalPi0RAAStat2760GeVV_2050->Clone("graphEMCalPi0RAAStat2760GeV");
// 		graphEMCalPi0RAASys2760GeV = (TGraphAsymmErrors*)graphEMCalPi0RAASys2760GeV_2050->Clone("graphEMCalPi0RAASys2760GeV");
// 		graphEMCalEtaRAAStat2760GeV = (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Clone("graphEMCalPi0RAAStat2760GeV");
// 		graphEMCalEtaRAASys2760GeV = (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Clone("graphEMCalPi0RAASys2760GeV");

		histoEMCalEtatoPi0StatPbPb2760GeV = (TH1D*)histoEMCalEtatoPi0StatPbPb2760GeV_2050->Clone("histoEMCalEtatoPi0StatPbPb2760GeV");
		graphEMCalEtatoPi0Stat2760GeV = new TGraphAsymmErrors(histoEMCalEtatoPi0StatPbPb2760GeV);
//		graphEMCalEtatoPi0Stat2760GeV = (TGraphAsymmErrors*)graphEMCalEtatoPi0Stat2760GeV_2050->Clone("graphEMCalEtatoPi0Stat2760GeV");
		graphEMCalEtatoPi0Sys2760GeV = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV_2050->Clone("graphEMCalEtatoPi0Sys2760GeV");

	}
			cout << "Eta 20-50%" << endl;		

// 	cout << "Spectra print: " << endl;
// 	if(meson.CompareTo("Pi0")==0){
// 		cout << "\n graphPCMPi0CorrYieldSysErrPbPb2760GeV" << endl;
// 			graphPCMPi0CorrYieldSysErrPbPb2760GeV->Print();
// 			
// 		cout << "\n graphPCMPi0InvYieldStatPbPb2760GeV" << endl;
// 			graphPCMPi0InvYieldStatPbPb2760GeV->Print();
// 
// 		cout << "\n graphPCMPi0InvYieldSysPbPb2760GeV" << endl;
// 			graphPCMPi0InvYieldSysPbPb2760GeV->Print();
// 			
// 		cout << "\n graphEMCalPi0CorrYieldSysErrPbPb2760GeV" << endl;
// 			graphEMCalPi0CorrYieldSysErrPbPb2760GeV->Print();
// 			
// 		cout << "\n graphEMCalPi0InvYieldStatPbPb2760GeV" << endl;
// 			graphEMCalPi0InvYieldStatPbPb2760GeV->Print();
// 			
// 		cout << "\n graphEMCalPi0InvYieldSysPbPb2760GeV" << endl;
// 			graphEMCalPi0InvYieldSysPbPb2760GeV->Print();
// 	} else {
// 		cout << "\n graphPCMEtaCorrYieldSysErrPbPb2760GeV" << endl;
// 			graphPCMEtaCorrYieldSysErrPbPb2760GeV->Print();
// 			
// 		cout << "\n graphPCMEtaInvYieldStatPbPb2760GeV" << endl;
// 			graphPCMEtaInvYieldStatPbPb2760GeV->Print();
// 			
// 		cout << "\n graphPCMEtaInvYieldSysPbPb2760GeV" << endl;
// 			graphPCMEtaInvYieldSysPbPb2760GeV->Print();
// 			
// 		cout << "\n graphEMCalEtaCorrYieldSysErrPbPb2760GeV" << endl;
// 			graphEMCalEtaCorrYieldSysErrPbPb2760GeV->Print();
// 			
// 		cout << "\n graphEMCalEtaInvYieldStatPbPb2760GeV" << endl;
// 			graphEMCalEtaInvYieldStatPbPb2760GeV->Print();
// 			
// 		cout << "\n graphEMCalEtaInvYieldSysPbPb2760GeV" << endl;
// 			graphEMCalEtaInvYieldSysPbPb2760GeV->Print();
// 	}
	
	// *******************************************************************************************************
	// ************************** Combination of different measurements **************************************
	// *******************************************************************************************************
	// REMARKS: 
	// 		- order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
	//		- correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
	// 		- currently only PCM-EMCal vs others fully implemeted energy independent
	// 		- extendable to other energies
	//		- offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented
		
	
	// definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
	// for statistical error and final value from respective method
	TH1D* statErrorCollection[11];
	TH1D* statErrorCollectionEtatoPi0[11];
	for (Int_t i = 0; i< 11; i++){
		statErrorCollection[i] = NULL;
	}	
	if(meson.CompareTo("Pi0")==0){
		statErrorCollection[0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV->Clone("statErrPCMPi0");
		statErrorCollection[2] = (TH1D*)histoEMCalPi0InvYieldPbPb2760GeV->Clone("statErrEMCalPi0");
	} else if(meson.CompareTo("Eta")==0) {
		statErrorCollection[0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV->Clone("statErrPCMEta");
		statErrorCollection[2] = (TH1D*)histoEMCalEtaInvYieldPbPb2760GeV->Clone("statErrEMCalEta");		

		statErrorCollectionEtatoPi0[0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV->Clone("statErrPCMEtatoPi0");
		statErrorCollectionEtatoPi0[2] = (TH1D*)histoEMCalEtatoPi0StatPbPb2760GeV->Clone("statErrEMCalEtatoPi0");		
		
	}
	
	// definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)	
	// for systematic error from respective method
	TGraphAsymmErrors* sysErrorCollection[11];
	TGraphAsymmErrors* sysErrorCollectionEtatoPi0[11];
	for (Int_t i = 0; i< 11; i++){
		sysErrorCollection[i] = NULL;
	}	
	if(meson.CompareTo("Pi0")==0){
		sysErrorCollection[0] = (TGraphAsymmErrors*)graphPCMPi0CorrYieldSysErrPbPb2760GeV->Clone("sysErrPCMPi0");
		sysErrorCollection[2] = (TGraphAsymmErrors*)graphEMCalPi0CorrYieldSysErrPbPb2760GeV->Clone("sysErrEMCalPi0");
	} else if(meson.CompareTo("Eta")==0) {
		sysErrorCollection[0] = (TGraphAsymmErrors*)graphPCMEtaCorrYieldSysErrPbPb2760GeV->Clone("sysErrPCMEta");
		sysErrorCollection[2] = (TGraphAsymmErrors*)graphEMCalEtaCorrYieldSysErrPbPb2760GeV->Clone("sysErrEMCalEta");	

		sysErrorCollectionEtatoPi0[0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV->Clone("sysErrPCMEtatoPi0");
		sysErrorCollectionEtatoPi0[2] = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV->Clone("sysErrEMCalEtatoPi0");	

	}
	
	// Definition of final pt binning (has to be set manually)
	Double_t xPtLimitsPi0[24] =  {0.0, 0.4, 0.6, 0.8, 1.0,
								  1.2, 1.4, 1.6, 1.8, 2.0,
								  2.2, 2.4, 2.6, 3.0, 3.5,
								  4.0, 5.0, 6.0, 8.0, 10.0, 
								  12.0, 14.,20.,30.};

	Double_t xPtLimitsEta[14] =  {0.0, 0.5, 1.0, 1.5, 2,
								  3.0, 4.0, 6.0, 8.0, 10,
								  12., 14., 20., 30.};
// 	Double_t xPtLimitsEta[18] =  {0.0, 0.6, 1.0, 1.4, 1.8,
// 								  2.2, 2.6, 3.0 ,3.5, 4., 5.,
// 				  				  6.0, 8.0, 10, 12., 14.,
// 								  20., 30.};

	// matrix order: 	= {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", 
	//					  "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged"};
	// Definition of offsets for stat & sys see output of function in shell, make sure pt bins match							
	Int_t offSetsPi0[11]	= 	{ 0, 0, 15, 0, 0,
								  0, 0, 0, 0, 0, 0};
	Int_t offSetsPi0Sys[11]= 	{ 1, 0, 15, 0, 0,
								  0, 0, 0, 0, 0, 0};
								
	Int_t offSetsEta[11]	= 	{ 0, 0, 6, 0, 0,
								  0, 0, 0, 0, 0, 0};//qui
	Int_t offSetsEtaSys[11]= 	{ 1, 0, 6, 0, 0,
								  0, 0, 0, 0, 0, 0};
// 	Int_t offSetsEta[11]	= 	{ 0, 0, 10, 0, 0,
// 								  0, 0, 0, 0, 0, 0};//qui
// 	Int_t offSetsEtaSys[11]= 	{ 1, 0, 10, 0, 0,
// 								  0, 0, 0, 0, 0, 0};

	//	**********************************************************************************************************************
	//	******************************************* Recalculating published spectrum *****************************************
	//	**********************************************************************************************************************
// 	TGraph* graphWeightsOld[11];
// 	for (Int_t i = 0; i< 11; i++){
// 		graphWeightsOld[i] = NULL;
// 	}	
						
	// Declaration & calculation of combined spectrum			
	if(meson.CompareTo("Pi0")==0){
// 		TString fileNameOutputWeightingOld							= Form("%s/%s_WeightingOldPi0.dat",outputDir.Data(),centrality.Data());
// 		TGraphAsymmErrors* graphCombPi0InvYieldStat2760GeVOld= NULL;
// 		TGraphAsymmErrors* graphCombPi0InvYieldSys2760GeVOld = NULL;
// 		TGraphAsymmErrors* graphCombPi0InvYieldTot2760GeVOld = CombinePtPointsSpectraFullCorrMat( 	statErrorCollection,	sysErrorCollection, 	
// 																											xPtLimitsPi0, 24, //combo Pi0 - 24 bins
// 																											offSetsPi0, offSetsPi0Sys,
// 																											graphCombPi0InvYieldStat2760GeVOld, graphCombPi0InvYieldSys2760GeVOld,
// 																											fileNameOutputWeightingOld,1
// 																										);
// 		cout << "Printing the graph:    " << endl;
// 		graphCombPi0InvYieldStat2760GeVOld->Print();
// 
// 		// Reading weights from output file for plotting
// 		ifstream fileWeightsReadOld;
// 		fileWeightsReadOld.open(fileNameOutputWeightingOld,ios_base::in);
// 		cout << "reading" << fileNameOutputWeightingOld << endl;
// 		Double_t xValuesReadOld[50];
// 		Double_t weightsReadOld[11][50];
// 		Int_t availableMeasOld[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
// 		Int_t nMeasSetOld 			= 2;
// 		Int_t nPtBinsReadOld 		= 0;
// 		while(!fileWeightsReadOld.eof() && nPtBinsReadOld < 50){
// 			TString garbage = "";
// 			if (nPtBinsReadOld == 0){
// 				fileWeightsReadOld >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
// 				for (Int_t i = 0; i < nMeasSetOld; i++){
// 					fileWeightsReadOld >> availableMeasOld[i] ;
// 				}	
// 				cout << "read following measurements: "; 
// 				for (Int_t i = 0; i < 11; i++){
// 					cout << availableMeasOld[i] << "\t" ;
// 				}	
// 				cout << endl;
// 			} else {
// 				fileWeightsReadOld >> xValuesReadOld[nPtBinsReadOld-1];
// 				for (Int_t i = 0; i < nMeasSetOld; i++){
// 					fileWeightsReadOld >> weightsReadOld[availableMeasOld[i]][nPtBinsReadOld-1] ;
// 				}	
// 				cout << "read: "<<  nPtBinsReadOld << "\t"<< xValuesReadOld[nPtBinsReadOld-1] << "\t" ;
// 				for (Int_t i = 0; i < nMeasSetOld; i++){
// 					cout << weightsReadOld[availableMeasOld[i]][nPtBinsReadOld-1] << "\t";
// 				}
// 				cout << endl;
// 			}
// 			nPtBinsReadOld++;
// 		}
// 		nPtBinsReadOld = nPtBinsReadOld-2 ;
// 		fileWeightsReadOld.close();
// 
// 		for (Int_t i = 0; i < nMeasSetOld; i++){
// 			graphWeightsOld[availableMeasOld[i]] = new TGraph(nPtBinsReadOld,xValuesReadOld,weightsReadOld[availableMeasOld[i]]);
// 			Int_t bin = 0;
// 			for (Int_t n = 0; n< nPtBinsReadOld; n++){
// 				if (graphWeightsOld[availableMeasOld[i]]->GetY()[bin] == 0) graphWeightsOld[availableMeasOld[i]]->RemovePoint(bin);
// 				else bin++;
// 			}	
// 		}	
	
	
		//	**********************************************************************************************************************
		//	******************************************* Plotting weights method only EMC *****************************************
		//	**********************************************************************************************************************
		Int_t textSizeLabelsPixel = 900*0.04;

		TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
		canvasWeights->SetLogx();
	
		TH2F * histo2DWeights;
		histo2DWeights = new TH2F("histo2DWeights","histo2DWeights",11000,0.23,70.,1000,-0.5,1.1);
		SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DWeights->GetXaxis()->SetMoreLogLabels();
		histo2DWeights->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DWeights->GetYaxis()->SetRangeUser(-10,10);
		histo2DWeights->Draw("copy");
// 		
// 		TLegend* legendAccWeightsOld = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetOld*1.35), 32);
// 		for (Int_t i = 0; i < nMeasSetOld; i++){
// 			DrawGammaSetMarkerTGraph(graphWeightsOld[availableMeasOld[i]],
// 									markerStyleDet[availableMeasOld[i]], 
// 									markerSizeDet[availableMeasOld[i]]*0.5, 
// 									colorDet[availableMeasOld[i]] , 
// 									colorDet[availableMeasOld[i]]);
// 			graphWeightsOld[availableMeasOld[i]]->Draw("p,same,e1");
// 			legendAccWeightsOld->AddEntry(graphWeightsOld[availableMeasOld[i]],nameMeasGlobal[availableMeasOld[i]],"p");
// 		}	
// 		legendAccWeightsOld->Draw();
// 
		TLatex *labelWeightsEnergy = new TLatex(0.7,0.20,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
		labelWeightsEnergy->SetTextFont(43);
// 		labelWeightsEnergy->Draw();
		TLatex *labelWeightsPi0 = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
		labelWeightsPi0->SetTextFont(43);
// 		labelWeightsPi0->Draw();
// 
// 	//	DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
// 		DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
// 		DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
// 		DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
// 		DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
// 			
// 		canvasWeights->SaveAs(Form("%s/%s_%s_WeightsOldPublished.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		
		//	**********************************************************************************************************************
		//	******************************************* Calculation of spectrum including EMCal only *****************************
		//	**********************************************************************************************************************
							
		TGraph* graphWeightsOEMC[11];
		for (Int_t i = 0; i< 11; i++){
			graphWeightsOEMC[i] = NULL;
		}	

		// Declaration & calculation of combined spectrum							
		TString fileNameOutputWeightingOEMC							= Form("%s/%s_WeightingEMCalPi0.dat",outputDir.Data(),centrality.Data());
		TGraphAsymmErrors* graphCombPi0InvYieldStat2760GeVOEMC= NULL;
		TGraphAsymmErrors* graphCombPi0InvYieldSys2760GeVOEMC = NULL;
		TGraphAsymmErrors* graphCombPi0InvYieldTot2760GeVOEMC = CombinePtPointsSpectraFullCorrMat( 	statErrorCollection,	sysErrorCollection, 	
																											xPtLimitsPi0, 23,
																											offSetsPi0, offSetsPi0Sys,
																											graphCombPi0InvYieldStat2760GeVOEMC, graphCombPi0InvYieldSys2760GeVOEMC,
																											fileNameOutputWeightingOEMC,1
																										);
		graphCombPi0InvYieldStat2760GeVOEMC->RemovePoint(0);
		graphCombPi0InvYieldSys2760GeVOEMC->RemovePoint(0);
		cout << "\n\n\n\n Printing of EMCal only spectra" << endl;
		graphCombPi0InvYieldStat2760GeVOEMC->Print();
		
		// Reading weights from output file for plotting
		ifstream fileWeightsReadOEMC;
		fileWeightsReadOEMC.open(fileNameOutputWeightingOEMC,ios_base::in);
		cout << "reading" << fileNameOutputWeightingOEMC << endl;
		Double_t xValuesReadOEMC[50];
		Double_t weightsReadOEMC[11][50];
		Int_t availableMeasOEMC[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
		Int_t nMeasSetOEMC 			= 2;
		Int_t nPtBinsReadOEMC 		= 0;
		while(!fileWeightsReadOEMC.eof() && nPtBinsReadOEMC < 50){
			TString garbage = "";
			if (nPtBinsReadOEMC == 0){
				fileWeightsReadOEMC >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
				for (Int_t i = 0; i < nMeasSetOEMC; i++){
					fileWeightsReadOEMC >> availableMeasOEMC[i] ;
				}	
				cout << "read following measurements: "; 
				for (Int_t i = 0; i < 11; i++){
					cout << availableMeasOEMC[i] << "\t" ;
				}	
				cout << endl;
			} else {
				fileWeightsReadOEMC >> xValuesReadOEMC[nPtBinsReadOEMC-1];
				for (Int_t i = 0; i < nMeasSetOEMC; i++){
					fileWeightsReadOEMC >> weightsReadOEMC[availableMeasOEMC[i]][nPtBinsReadOEMC-1] ;
				}	
				cout << "read: "<<  nPtBinsReadOEMC << "\t"<< xValuesReadOEMC[nPtBinsReadOEMC-1] << "\t" ;
				for (Int_t i = 0; i < nMeasSetOEMC; i++){
					cout << weightsReadOEMC[availableMeasOEMC[i]][nPtBinsReadOEMC-1] << "\t";
				}
				cout << endl;
			}
			nPtBinsReadOEMC++;
		}
		nPtBinsReadOEMC = nPtBinsReadOEMC-2 ;
		fileWeightsReadOEMC.close();

		for (Int_t i = 0; i < nMeasSetOEMC; i++){
			graphWeightsOEMC[availableMeasOEMC[i]] = new TGraph(nPtBinsReadOEMC,xValuesReadOEMC,weightsReadOEMC[availableMeasOEMC[i]]);
			Int_t bin = 0;
			for (Int_t n = 0; n< nPtBinsReadOEMC; n++){
				if (graphWeightsOEMC[availableMeasOEMC[i]]->GetY()[bin] == 0) graphWeightsOEMC[availableMeasOEMC[i]]->RemovePoint(bin);
				else bin++;
			}	
		}	

		//	**********************************************************************************************************************
		//	******************************************* Plotting weights method only EMC *****************************************
		//	**********************************************************************************************************************
		canvasWeights->cd();
		histo2DWeights->Draw("copy");

		TLegend* legendAccWeightsOEMC = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetOEMC*1.35), 32);
		for (Int_t i = 0; i < nMeasSetOEMC; i++){
			DrawGammaSetMarkerTGraph(graphWeightsOEMC[availableMeasOEMC[i]],
									markerStyleDet[availableMeasOEMC[i]], 
									markerSizeDet[availableMeasOEMC[i]]*0.5, 
									colorDet[availableMeasOEMC[i]] , 
									colorDet[availableMeasOEMC[i]]);
			graphWeightsOEMC[availableMeasOEMC[i]]->Draw("p,same,e1");
			legendAccWeightsOEMC->AddEntry(graphWeightsOEMC[availableMeasOEMC[i]],nameMeasGlobal[availableMeasOEMC[i]],"p");
		}	
		legendAccWeightsOEMC->Draw();
		labelWeightsEnergy->Draw();
		labelWeightsPi0->Draw();

	//	DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
		DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
		DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
			
		canvasWeights->SaveAs(Form("%s/%s_%s_WeightsOnlyEMCal.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		
		TH1D* statErrorRelCollection[11];
		for (Int_t i = 0; i< 11; i++){
			statErrorRelCollection[i] = NULL;
		}	
		for (Int_t i = 0; i < 11; i++){
			if (statErrorCollection[i]) statErrorRelCollection[i] = CalculateRelErrUpTH1D( statErrorCollection[i], Form("relativeStatError_%s", nameMeasGlobal[i].Data()));
		}
		
		TGraphAsymmErrors* sysErrorRelCollection[11];
		for (Int_t i = 0; i< 11; i++){
			sysErrorRelCollection[i] = NULL;
		}	
		for (Int_t i = 0; i < 11; i++){
			if (sysErrorCollection[i]) sysErrorRelCollection[i] = CalculateRelErrUpAsymmGraph( sysErrorCollection[i], Form("relativeSysError_%s", nameMeasGlobal[i].Data()));
		}
		
		
		//	**********************************************************************************************************************
		//	******************************************* Assuming maximal correlation *********************************************
		//	**********************************************************************************************************************
									
		TGraph* graphWeightsA[11];
		for (Int_t i = 0; i< 11; i++){
			graphWeightsA[i] = NULL;
		}	
									
		// Declaration & calculation of combined spectrum							
		TString fileNameOutputWeightingA 							= Form("%s/%s_WeightingMethodAPi0.dat",outputDir.Data(),centrality.Data());
		TGraphAsymmErrors* graphCombPi0InvYieldStat2760GeVA 	= NULL;
		TGraphAsymmErrors* graphCombPi0InvYieldSys2760GeVA 	= NULL;
		TGraphAsymmErrors* graphCombPi0InvYieldTot2760GeVA 	= CombinePtPointsSpectraFullCorrMat( 	statErrorCollection,	sysErrorCollection, 	
																											xPtLimitsPi0, 23,
																											offSetsPi0, offSetsPi0Sys,
																											graphCombPi0InvYieldStat2760GeVA, graphCombPi0InvYieldSys2760GeVA,
																											fileNameOutputWeightingA,1
																										);
		graphCombPi0InvYieldStat2760GeVA->RemovePoint(0);
        graphCombPi0InvYieldSys2760GeVA->RemovePoint(0);
        graphCombPi0InvYieldTot2760GeVA->RemovePoint(0);

// 		cout << "\n\n\n\n Printing spectra method A" << endl;
// 		graphCombPi0InvYieldStat2760GeVA->Print();
//         graphCombPi0InvYieldSys2760GeVA->Print();
//         graphCombPi0InvYieldTot2760GeVA->Print();

		
		// Reading weights from output file for plotting
		ifstream fileWeightsReadA;
		fileWeightsReadA.open(fileNameOutputWeightingA,ios_base::in);
		cout << "reading" << fileNameOutputWeightingA << endl;
		Double_t xValuesReadA[50];
		Double_t weightsReadA[11][50];
		Int_t availableMeasA[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
		Int_t nMeasSetA 			= 2;
		Int_t nPtBinsReadA 		= 0;
		while(!fileWeightsReadA.eof() && nPtBinsReadA < 50){
			TString garbage = "";
			if (nPtBinsReadA == 0){
				fileWeightsReadA >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
				for (Int_t i = 0; i < nMeasSetA; i++){
					fileWeightsReadA >> availableMeasA[i] ;
				}	
				cout << "read following measurements: "; 
				for (Int_t i = 0; i < 11; i++){
					cout << availableMeasA[i] << "\t" ;
				}	
				cout << endl;
			} else {
				fileWeightsReadA >> xValuesReadA[nPtBinsReadA-1];
				for (Int_t i = 0; i < nMeasSetA; i++){
					fileWeightsReadA >> weightsReadA[availableMeasA[i]][nPtBinsReadA-1] ;
				}	
				cout << "read: "<<  nPtBinsReadA << "\t"<< xValuesReadA[nPtBinsReadA-1] << "\t" ;
				for (Int_t i = 0; i < nMeasSetA; i++){
					cout << weightsReadA[availableMeasA[i]][nPtBinsReadA-1] << "\t";
				}
				cout << endl;
			}
			nPtBinsReadA++;
		}
		nPtBinsReadA = nPtBinsReadA-2 ;
		fileWeightsReadA.close();

		for (Int_t i = 0; i < nMeasSetA; i++){
			graphWeightsA[availableMeasA[i]] = new TGraph(nPtBinsReadA,xValuesReadA,weightsReadA[availableMeasA[i]]);
			Int_t bin = 0;
			for (Int_t n = 0; n< nPtBinsReadA; n++){
				if (graphWeightsA[availableMeasA[i]]->GetY()[bin] == 0) graphWeightsA[availableMeasA[i]]->RemovePoint(bin);
				else bin++;
			}	
		}	

		//	**********************************************************************************************************************
		//	******************************************* Plotting weights Method A ************************************************
		//	**********************************************************************************************************************
		canvasWeights->cd();
	
		histo2DWeights->Draw("copy");

		TLegend* legendAccWeights = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetA*1.35), 32);
		for (Int_t i = 0; i < nMeasSetA; i++){
			DrawGammaSetMarkerTGraph(graphWeightsA[availableMeasA[i]], markerStyleDet[availableMeasA[i]], markerSizeDet[availableMeasA[i]]*0.5, colorDet[availableMeasA[i]] , colorDet[availableMeasA[i]]);
			graphWeightsA[availableMeasA[i]]->Draw("p,same,e1");
			legendAccWeights->AddEntry(graphWeightsA[availableMeasA[i]],nameMeasGlobal[availableMeasA[i]],"p");
		}	
		legendAccWeights->Draw();

		labelWeightsEnergy->Draw();
		labelWeightsPi0->Draw();

// 		DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
		DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
		DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
		
		canvasWeights->SaveAs(Form("%s/%s_%s_WeightsA.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		//	**********************************************************************************************************************
		//	******************************************* Assuming maximal correlation *********************************************
		//	**********************************************************************************************************************
		
		TGraph* graphWeightsB[11];
		for (Int_t i = 0; i< 11; i++){
			graphWeightsB[i] = NULL;
		}	
									
		// Declaration & calculation of combined spectrum							
		TString fileNameOutputWeightingB 							= Form("%s/%s_WeightingMethodBPi0.dat",outputDir.Data(),centrality.Data());
		TGraphAsymmErrors* graphCombPi0InvYieldStat2760GeVB 	= NULL;
		TGraphAsymmErrors* graphCombPi0InvYieldSys2760GeVB 	= NULL;
		TGraphAsymmErrors* graphCombPi0InvYieldTot2760GeVB 	= CombinePtPointsSpectraFullCorrMat( 	statErrorCollection,	sysErrorCollection, 	
																											xPtLimitsPi0, 23,
																											offSetsPi0, offSetsPi0Sys,
																											graphCombPi0InvYieldStat2760GeVB, graphCombPi0InvYieldSys2760GeVB,
																											fileNameOutputWeightingB,2
																										);
		graphCombPi0InvYieldStat2760GeVB->RemovePoint(0);
		graphCombPi0InvYieldTot2760GeVB->RemovePoint(0);
		graphCombPi0InvYieldSys2760GeVB->RemovePoint(0);
// 		cout << "\n\n\n\n Printing spectra method B" << endl;
// 		graphCombPi0InvYieldStat2760GeVB->Print();
//         graphCombPi0InvYieldSys2760GeVB->Print();
//         graphCombPi0InvYieldTot2760GeVB->Print();

		
		// Reading weights from output file for plotting
		ifstream fileWeightsReadB;
		fileWeightsReadB.open(fileNameOutputWeightingB,ios_base::in);
		cout << "reading" << fileNameOutputWeightingB << endl;
		Double_t xValuesReadB[50];
		Double_t weightsReadB[11][50];
		Int_t availableMeasB[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
		Int_t nMeasSetB 			= 2;
		Int_t nPtBinsReadB 		= 0;
		while(!fileWeightsReadB.eof() && nPtBinsReadB < 50){
			TString garbage = "";
			if (nPtBinsReadB == 0){
				fileWeightsReadB >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
				for (Int_t i = 0; i < nMeasSetB; i++){
					fileWeightsReadB >> availableMeasB[i] ;
				}	
				cout << "read following measurements: "; 
				for (Int_t i = 0; i < 11; i++){
					cout << availableMeasB[i] << "\t" ;
				}	
				cout << endl;
			} else {
				fileWeightsReadB >> xValuesReadB[nPtBinsReadB-1];
				for (Int_t i = 0; i < nMeasSetB; i++){
					fileWeightsReadB >> weightsReadB[availableMeasB[i]][nPtBinsReadB-1] ;
				}	
				cout << "read: "<<  nPtBinsReadB << "\t"<< xValuesReadB[nPtBinsReadB-1] << "\t" ;
				for (Int_t i = 0; i < nMeasSetB; i++){
					cout << weightsReadB[availableMeasB[i]][nPtBinsReadB-1] << "\t";
				}
				cout << endl;
			}
			nPtBinsReadB++;
		}
		nPtBinsReadB = nPtBinsReadB-2 ;
		fileWeightsReadB.close();

		for (Int_t i = 0; i < nMeasSetB; i++){
			graphWeightsB[availableMeasB[i]] = new TGraph(nPtBinsReadB,xValuesReadB,weightsReadB[availableMeasB[i]]);
			Int_t bin = 0;
			for (Int_t n = 0; n< nPtBinsReadB; n++){
				if (graphWeightsB[availableMeasB[i]]->GetY()[bin] == 0) graphWeightsB[availableMeasB[i]]->RemovePoint(bin);
				else bin++;
			}	
		}	

		//	**********************************************************************************************************************
		//	******************************************* Plotting weights method B ************************************************
		//	**********************************************************************************************************************
		
		canvasWeights->cd();
	
		histo2DWeights->Draw("copy");

		for (Int_t i = 0; i < nMeasSetB; i++){
			DrawGammaSetMarkerTGraph(graphWeightsB[availableMeasB[i]], markerStyleDet[availableMeasB[i]], markerSizeDet[availableMeasB[i]]*0.5, colorDet[availableMeasB[i]] , colorDet[availableMeasB[i]]);
			graphWeightsB[availableMeasB[i]]->Draw("p,same,e1");
		}	
		legendAccWeights->Draw();

		labelWeightsEnergy->Draw();
		labelWeightsPi0->Draw();

// 		DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
		DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
		DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
		
		canvasWeights->SaveAs(Form("%s/%s_%s_WeightsB.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		

		//	**********************************************************************************************************************
		//	******************************************* Compare new spectrum to old average **************************************
		//	**********************************************************************************************************************
		
// 		TGraphAsymmErrors* graphCombPi0InvYieldStat2760GeVAForCompToOld	 = (TGraphAsymmErrors*)graphCombPi0InvYieldStat2760GeVA->Clone("graphCombPi0InvYieldStat2760GeVAForCompToOld");
// 		TGraphAsymmErrors* graphCombPi0InvYieldSys2760GeVAForCompToOld	 = (TGraphAsymmErrors*)graphCombPi0InvYieldSys2760GeVA->Clone("graphCombPi0InvYieldSys2760GeVAForCompToOld");
// 		
// 		TGraphAsymmErrors* graphCombPi0InvYieldStat2760GeVBForCompToOld	 = (TGraphAsymmErrors*)graphCombPi0InvYieldStat2760GeVB->Clone("graphCombPi0InvYieldStat2760GeVBForCompToOld");
// 		TGraphAsymmErrors* graphCombPi0InvYieldSys2760GeVBForCompToOld	 = (TGraphAsymmErrors*)graphCombPi0InvYieldSys2760GeVB->Clone("graphCombPi0InvYieldSys2760GeVBForCompToOld");
// 
// 		TGraphAsymmErrors* graphCombPi0InvYieldStat2760GeVOEMCBForCompToOld= (TGraphAsymmErrors*)graphCombPi0InvYieldStat2760GeVOEMC->Clone("graphCombPi0InvYieldStat2760GeVOEMCBForCompToOld");
// 		TGraphAsymmErrors* graphCombPi0InvYieldSys2760GeVOEMCForCompToOld = (TGraphAsymmErrors*)graphCombPi0InvYieldSys2760GeVOEMC->Clone("graphCombPi0InvYieldSys2760GeVOEMCForCompToOld");
// 		
// 		for (Int_t n = graphCombPi0InvYieldStat2760GeVAForCompToOld->GetN()-1; graphCombPi0InvYieldStat2760GeVAForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombPi0InvYieldStat2760GeVAForCompToOld->RemovePoint(n);
// 		}	
// 		for (Int_t n = graphCombPi0InvYieldSys2760GeVAForCompToOld->GetN()-1; graphCombPi0InvYieldSys2760GeVAForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombPi0InvYieldSys2760GeVAForCompToOld->RemovePoint(n);
// 		}	
// 		for (Int_t n = graphCombPi0InvYieldStat2760GeVBForCompToOld->GetN()-1; graphCombPi0InvYieldStat2760GeVBForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombPi0InvYieldStat2760GeVBForCompToOld->RemovePoint(n);
// 		}	
// 		for (Int_t n = graphCombPi0InvYieldSys2760GeVBForCompToOld->GetN()-1; graphCombPi0InvYieldSys2760GeVBForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombPi0InvYieldSys2760GeVBForCompToOld->RemovePoint(n);
// 		}	
// 		for (Int_t n = graphCombPi0InvYieldStat2760GeVOEMCBForCompToOld->GetN()-1; graphCombPi0InvYieldStat2760GeVOEMCBForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombPi0InvYieldStat2760GeVOEMCBForCompToOld->RemovePoint(n);
// 		}	
// 		for (Int_t n = graphCombPi0InvYieldSys2760GeVOEMCForCompToOld->GetN()-1; graphCombPi0InvYieldSys2760GeVOEMCForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombPi0InvYieldSys2760GeVOEMCForCompToOld->RemovePoint(n);
// 		}	
// 
// 		graphCombPi0InvYieldStat2760GeVAForCompToOld->Print();
// 		graphCombPi0InvYieldSys2760GeVAForCompToOld->Print();
// 		graphCombPi0InvYieldStat2760GeVBForCompToOld->Print();
// 		graphCombPi0InvYieldSys2760GeVBForCompToOld->Print();
// 		
// 		TGraphAsymmErrors*  graphRatioCombNewADivCombOldStat 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombPi0InvYieldStat2760GeVAForCompToOld, graphCombPi0InvYieldStat2760GeVOld, kTRUE);
// 		TGraphAsymmErrors*  graphRatioCombNewADivCombOldSys 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombPi0InvYieldSys2760GeVAForCompToOld, graphCombPi0InvYieldSys2760GeVOld, kTRUE);
// 		TGraphAsymmErrors*  graphRatioCombNewBDivCombOldStat 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombPi0InvYieldStat2760GeVBForCompToOld, graphCombPi0InvYieldStat2760GeVOld, kTRUE);
// 		TGraphAsymmErrors*  graphRatioCombNewBDivCombOldSys 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombPi0InvYieldSys2760GeVBForCompToOld, graphCombPi0InvYieldSys2760GeVOld, kTRUE);
// 		TGraphAsymmErrors*  graphRatioCombNewOEMCDivCombOldStat = CalculateGraphErrRatioToOtherTGraphErr(graphCombPi0InvYieldStat2760GeVOEMCBForCompToOld,
// 																										graphCombPi0InvYieldStat2760GeVOld, kTRUE);
// 		TGraphAsymmErrors*  graphRatioCombNewOECMDivCombOldSys 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombPi0InvYieldSys2760GeVOEMCForCompToOld, 
// 																										graphCombPi0InvYieldSys2760GeVOld, kTRUE);
		

		
		

		//	**********************************************************************************************************************
		//	******************************************* Ratio of Comb to Fit ****************************************
		//	**********************************************************************************************************************
		
// 		TCanvas* canvasRatioToOldCombined = new TCanvas("canvasRatioToOldCombined","",200,10,1350,900);  // gives the page size
// 		DrawGammaCanvasSettings( canvasRatioToOldCombined, 0.1, 0.02, 0.035, 0.09);
// 		canvasRatioToOldCombined->SetLogx();
// 	
// 		TH2F * histo2DRatioToOldCombined;
// 		histo2DRatioToOldCombined = new TH2F("histo2DRatioToOldCombined","histo2DRatioToOldCombined",11000,0.23,10.,1000,0.75,1.25);
// 		SetStyleHistoTH2ForGraphs(histo2DRatioToOldCombined, "#it{p}_{T} (GeV/#it{c})","#frac{Comb}{Comb Old}",0.035,0.04, 0.035,0.04, 1.,1.,510,505);
// 		histo2DRatioToOldCombined->GetXaxis()->SetMoreLogLabels();
// 		histo2DRatioToOldCombined->GetXaxis()->SetLabelOffset(-0.01);
// 	// 	histo2DRatioToOldCombined->GetYaxis()->SetRangeUser(-10,10);
// 		histo2DRatioToOldCombined->Draw("copy");
// 
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewBDivCombOldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
// 			graphRatioCombNewBDivCombOldSys->Draw("E2same");
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewBDivCombOldStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
// 			graphRatioCombNewBDivCombOldStat->Draw("p,same,e1");
// 
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewADivCombOldSys, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2, widthLinesBoxes, kTRUE);
// 			graphRatioCombNewADivCombOldSys->Draw("E2same");
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewADivCombOldStat, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
// 			graphRatioCombNewADivCombOldStat->Draw("p,same,e1");
// 
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewOECMDivCombOldSys, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2, widthLinesBoxes, kTRUE);
// 			graphRatioCombNewOECMDivCombOldSys->Draw("E2same");
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewOEMCDivCombOldStat, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
// 			graphRatioCombNewOEMCDivCombOldStat->Draw("p,same,e1");
// 			
// 			DrawGammaLines(0.23, 10. , 1., 1.,0.1, kGray+2);
// 			DrawGammaLines(0.23, 10. , 1.1, 1.1,0.1, kGray, 7);
// 			DrawGammaLines(0.23, 10. , 0.9, 0.9,0.1, kGray, 7);
// 			DrawGammaLines(0.23, 10. , 1.05, 1.05,0.1, kGray, 3);
// 			DrawGammaLines(0.23, 10. , 0.95, 0.95,0.1, kGray, 3);
// 			
// 			TLegend* legendRatioToOld = GetAndSetLegend2(0.67, 0.96-(0.035*3*1.35), 0.93, 0.96, 32);
// 			legendRatioToOld->AddEntry(graphRatioCombNewOECMDivCombOldSys,"PCM, EMCal");
// 			legendRatioToOld->AddEntry(graphRatioCombNewADivCombOldSys,"Method A");
// 			legendRatioToOld->AddEntry(graphRatioCombNewBDivCombOldSys,"Method B");
// 			legendRatioToOld->Draw();
// 
			TLatex *labelRatioToOldEnergy = new TLatex(0.15,0.89,collisionSystem2760GeV.Data());
			SetStyleTLatex( labelRatioToOldEnergy, 0.85*textSizeLabelsPixel,4);
			labelRatioToOldEnergy->SetTextFont(43);
// 			labelRatioToOldEnergy->Draw();
			TLatex *labelRatioToOldPi0 = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
			SetStyleTLatex( labelRatioToOldPi0, 0.85*textSizeLabelsPixel,4);
			labelRatioToOldPi0->SetTextFont(43);
// 			labelRatioToOldPi0->Draw();
// 			
// 		canvasRatioToOldCombined->SaveAs(Form("%s/%s_%s_RatioOfCombToCombOld_PbPb2760GeV.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
// 
		//	**********************************************************************************************************************
		//	************************************** Ratio of Comb to Fit without Method A *****************************************
		//	**********************************************************************************************************************

// 		canvasRatioToOldCombined->cd();
// 		histo2DRatioToOldCombined->Draw("copy");
// 
// 			graphRatioCombNewBDivCombOldSys->Draw("E2same");
// 			graphRatioCombNewBDivCombOldStat->Draw("p,same,e1");
// 
// 			graphRatioCombNewOECMDivCombOldSys->Draw("E2same");
// 			graphRatioCombNewOEMCDivCombOldStat->Draw("p,same,e1");
// 			
// 			DrawGammaLines(0.23, 10. , 1., 1.,0.1, kGray+2);
// 			DrawGammaLines(0.23, 10. , 1.1, 1.1,0.1, kGray, 7);
// 			DrawGammaLines(0.23, 10. , 0.9, 0.9,0.1, kGray, 7);
// 			DrawGammaLines(0.23, 10. , 1.05, 1.05,0.1, kGray, 3);
// 			DrawGammaLines(0.23, 10. , 0.95, 0.95,0.1, kGray, 3);
// 			
// 			TLegend* legendRatioToOld2 = GetAndSetLegend2(0.67, 0.96-(0.035*3*1.35), 0.93, 0.96, 32);
// 			legendRatioToOld2->AddEntry(graphRatioCombNewOECMDivCombOldSys,"PCM, EMCal");
// 			legendRatioToOld2->AddEntry(graphRatioCombNewBDivCombOldSys,"All");
// 			legendRatioToOld2->Draw();
// 
// 			labelRatioToOldEnergy->Draw();
// 			labelRatioToOldPi0->Draw();
// 			
// 		canvasRatioToOldCombined->SaveAs(Form("%s/%s_%s_RatioOfCombToCombOld_PbPb2760GeV_withoutMethodA.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
// 
// 
		
		// 	*********************************************************************************************************************
		// 	************************************ Visualize relative errors ******************************************************
		// 	*********************************************************************************************************************
		
		TCanvas* canvasRelSysErr = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
		canvasRelSysErr->SetLogx();
	
		TH2F * histo2DRelSysErr;
		histo2DRelSysErr = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23,70.,1000,0,80.5);
		SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
		histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DRelSysErr->GetYaxis()->SetRangeUser(-10,10);
		histo2DRelSysErr->Draw("copy");	
			TLegend* legendRelSysErr = GetAndSetLegend2(0.62, 0.94-(0.035*nMeasSetA*1.35), 0.95, 0.94, 32);
			for (Int_t i = 0; i < nMeasSetA; i++){
				DrawGammaSetMarkerTGraph(sysErrorRelCollection[availableMeasA[i]], markerStyleDet[availableMeasA[i]], markerSizeDet[availableMeasA[i]]*0.5, colorDet[availableMeasA[i]],
										colorDet[availableMeasA[i]]);
				sysErrorRelCollection[availableMeasA[i]]->Draw("p,same,e1");
				legendRelSysErr->AddEntry(sysErrorRelCollection[availableMeasA[i]],nameMeasGlobal[availableMeasA[i]],"p");
			}	
			legendRelSysErr->Draw();

			TLatex *labelRelSysErrEnergy = new TLatex(0.15,0.89,collisionSystem2760GeV.Data());
			SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel,4);
			labelRelSysErrEnergy->SetTextFont(43);
			labelRelSysErrEnergy->Draw();
			TLatex *labelRelSysErrPi0 = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
			SetStyleTLatex( labelRelSysErrPi0, 0.85*textSizeLabelsPixel,4);
			labelRelSysErrPi0->SetTextFont(43);
			labelRelSysErrPi0->Draw();
			
		canvasRelSysErr->SaveAs(Form("%s/%s_%s_RelSysErr.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		histo2DRelSysErr->GetYaxis()->SetRangeUser(0,30.5);
		histo2DRelSysErr->Draw("copy");	

			for (Int_t i = 0; i < nMeasSetA; i++){
				sysErrorRelCollection[availableMeasA[i]]->Draw("p,same,e1");
			}	
			legendRelSysErr->Draw();

			labelRelSysErrEnergy->Draw();
			labelRelSysErrPi0->Draw();
			
		canvasRelSysErr->SaveAs(Form("%s/%s_%s_RelSysErrZoomed.%s",outputDir.Data(),meson.Data(), centrality.Data(),suffix.Data()));
		
		
		
		// 	*********************************************************************************************************************
		// 	************************************ Visualize relative errors ******************************************************
		// 	*********************************************************************************************************************
		
		TCanvas* canvasRelStatErr = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
		canvasRelStatErr->SetLogx();
	
		TH2F * histo2DRelStatErr;
		histo2DRelStatErr = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23,70.,1000,0,80.5);
		SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
		histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DRelStatErr->GetYaxis()->SetRangeUser(-10,10);
		histo2DRelStatErr->Draw("copy");	
//  			statErrorRelCollection[0]->SetBinContent(statErrorRelCollection[0]->GetNbinsX(), -1);
		
			TLegend* legendRelStatErr = GetAndSetLegend2(0.14, 0.94-(0.035*nMeasSetA*1.35), 0.45, 0.94, 32);
			for (Int_t i = 0; i < nMeasSetA; i++){
				DrawGammaSetMarker(statErrorRelCollection[availableMeasA[i]], markerStyleDet[availableMeasA[i]], markerSizeDet[availableMeasA[i]]*0.5, colorDet[availableMeasA[i]] , 
								colorDet[availableMeasA[i]]);
				statErrorRelCollection[availableMeasA[i]]->Draw("p,same,e1");
				legendRelStatErr->AddEntry(statErrorRelCollection[availableMeasA[i]],nameMeasGlobal[availableMeasA[i]],"p");
			}	
			legendRelStatErr->Draw();

			TLatex *labelRelStatErrEnergy = new TLatex(0.75,0.89,collisionSystem2760GeV.Data());
			SetStyleTLatex( labelRelStatErrEnergy, 0.85*textSizeLabelsPixel,4);
			labelRelStatErrEnergy->SetTextFont(43);
			labelRelStatErrEnergy->Draw();
			TLatex *labelRelStatErrPi0 = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
			SetStyleTLatex( labelRelStatErrPi0, 0.85*textSizeLabelsPixel,4);
			labelRelStatErrPi0->SetTextFont(43);
			labelRelStatErrPi0->Draw();
			
		canvasRelStatErr->SaveAs(Form("%s/%s_%s_RelStatErr.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

// 		histo2DRelStatErr->GetYaxis()->SetRangeUser(0,30.5);
// 		histo2DRelStatErr->Draw("copy");	
// 			statErrorRelCollection[1]->SetBinContent(statErrorRelCollection[1]->GetNbinsX(), -1);
// 			statErrorRelCollection[2]->SetBinContent(statErrorRelCollection[2]->GetNbinsX(), -1);
// 			statErrorRelCollection[2]->SetBinContent(statErrorRelCollection[2]->GetNbinsX()-1, -1);
// 			for (Int_t i = 0; i < nMeasSetA; i++){
// 				statErrorRelCollection[availableMeasA[i]]->Draw("p,same,e1");
// 			}	
// 			legendRelStatErr->Draw();
// 
// 			labelRelStatErrEnergy->Draw();
// 			labelRelStatErrPi0->Draw();
// 			
// 		canvasRelStatErr->SaveAs(Form("%s/%s_RelStatErrZoomed.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		
		// 	*********************************************************************************************************************
		// 	************************ Visualize relative total errors of different combination methods ***************************
		// 	*********************************************************************************************************************
		TGraphAsymmErrors* graphCombPi0InvYieldRelTot2760GeVB 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldTot2760GeVB, "relativeTotalError_MethodB");
		TGraphAsymmErrors* graphCombPi0InvYieldRelStat2760GeVB 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldStat2760GeVB, "relativeStatError_MethodB");
		TGraphAsymmErrors* graphCombPi0InvYieldRelSys2760GeVB 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldSys2760GeVB, "relativeSysError_MethodB");
		TGraphAsymmErrors* graphCombPi0InvYieldRelStat2760GeVA 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldStat2760GeVA, "relativeStatError_MethodA");
		TGraphAsymmErrors* graphCombPi0InvYieldRelSys2760GeVA 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldSys2760GeVA, "relativeSysError_MethodA");

		TGraphAsymmErrors* graphCombPi0InvYieldRelTot2760GeVA 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldTot2760GeVA, "relativeTotalError_MethodA");
		TGraphAsymmErrors* graphCombPi0InvYieldRelTot2760GeVOEMC = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldTot2760GeVOEMC, "relativeTotalError_OEMC");
// 		TGraphAsymmErrors* graphCombPi0InvYieldRelTot2760GeVOld 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldTot2760GeVOld, "relativeTotalError_Old");

		
		TCanvas* canvasRelTotErr = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
		canvasRelTotErr->SetLogx();
	
		TH2F * histo2DRelTotErr;
		histo2DRelTotErr = new TH2F("histo2DRelTotErr","histo2DRelTotErr",11000,0.23,70.,1000,0,80.5);
		SetStyleHistoTH2ForGraphs(histo2DRelTotErr, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DRelTotErr->GetXaxis()->SetMoreLogLabels();
		histo2DRelTotErr->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DRelTotErr->GetYaxis()->SetRangeUser(-10,10);
		histo2DRelTotErr->Draw("copy");	

// 			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot2760GeVOld, markerStyleComb+1, markerSizeComb, kBlack , kBlack);
// 			graphCombPi0InvYieldRelTot2760GeVOld->Draw("p,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot2760GeVB, markerStyleComb, markerSizeComb, colorComb , colorComb);
			graphCombPi0InvYieldRelTot2760GeVB->Draw("p,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot2760GeVA, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
			graphCombPi0InvYieldRelTot2760GeVA->Draw("p,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot2760GeVOEMC, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
			graphCombPi0InvYieldRelTot2760GeVOEMC->Draw("p,same,e1");

			TLegend* legendRelTotErr = GetAndSetLegend2(0.14, 0.94-(0.035*4*1.35), 0.45, 0.94, 32);
// 			legendRelTotErr->AddEntry(graphCombPi0InvYieldRelTot2760GeVOld,"PCMold, EMCalold","p");
			legendRelTotErr->AddEntry(graphCombPi0InvYieldRelTot2760GeVOEMC,"PCM, EMCal","p");
			legendRelTotErr->AddEntry(graphCombPi0InvYieldRelTot2760GeVA,"Method A","p");
			legendRelTotErr->AddEntry(graphCombPi0InvYieldRelTot2760GeVB,"Method B","p");
			legendRelTotErr->Draw();

			TLatex *labelRelTotErrEnergy = new TLatex(0.75,0.89,collisionSystem2760GeV.Data());
			SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel,4);
			labelRelTotErrEnergy->SetTextFont(43);
			labelRelTotErrEnergy->Draw();
			TLatex *labelRelTotErrPi0 = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
			SetStyleTLatex( labelRelTotErrPi0, 0.85*textSizeLabelsPixel,4);
			labelRelTotErrPi0->SetTextFont(43);
			labelRelTotErrPi0->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelTotErr.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		histo2DRelTotErr->Draw("copy");	
// 			graphCombPi0InvYieldRelTot2760GeVOld->Draw("p,same,e1");
			graphCombPi0InvYieldRelTot2760GeVB->Draw("p,same,e1");
			graphCombPi0InvYieldRelTot2760GeVA->Draw("p,same,e1");
			graphCombPi0InvYieldRelTot2760GeVOEMC->Draw("p,same,e1");

			TLegend* legendRelTotErrWOA = GetAndSetLegend2(0.14, 0.94-(0.035*4*1.35), 0.45, 0.94, 32);
// 			legendRelTotErrWOA->AddEntry(graphCombPi0InvYieldRelTot2760GeVOld,"PCMold, EMCalold","p");
			legendRelTotErrWOA->AddEntry(graphCombPi0InvYieldRelTot2760GeVOEMC,"PCM, EMCal","p");
			legendRelTotErrWOA->AddEntry(graphCombPi0InvYieldRelTot2760GeVB,"All","p");
			legendRelTotErrWOA->Draw();

			labelRelTotErrEnergy->Draw();
			labelRelTotErrPi0->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelTotErr_withoutMethodA.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		
		histo2DRelTotErr->GetYaxis()->SetRangeUser(0,30.5);
		histo2DRelTotErr->Draw("copy");	
// 			graphCombPi0InvYieldRelTot2760GeVOld->Draw("p,same,e1");
			graphCombPi0InvYieldRelTot2760GeVB->Draw("p,same,e1");
			graphCombPi0InvYieldRelTot2760GeVA->Draw("p,same,e1");
			graphCombPi0InvYieldRelTot2760GeVOEMC->Draw("p,same,e1");

			legendRelTotErr->Draw();

			labelRelTotErrEnergy->Draw();
			labelRelTotErrPi0->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelTotErrZoomed.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		histo2DRelTotErr->Draw("copy");	
// 			graphCombPi0InvYieldRelTot2760GeVOld->Draw("p,same,e1");
			graphCombPi0InvYieldRelTot2760GeVB->Draw("p,same,e1");
			graphCombPi0InvYieldRelTot2760GeVOEMC->Draw("p,same,e1");

			legendRelTotErrWOA->Draw();

			labelRelTotErrEnergy->Draw();
			labelRelTotErrPi0->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelTotErrZoomed_withoutMethodA.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		histo2DRelTotErr->GetYaxis()->SetRangeUser(0,80.5);
		histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
			histo2DRelTotErr->Draw("copy");	

			graphCombPi0InvYieldRelTot2760GeVB->Draw("p,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelStat2760GeVB, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
			graphCombPi0InvYieldRelStat2760GeVB->Draw("l,x0,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelSys2760GeVB, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
			graphCombPi0InvYieldRelSys2760GeVB->SetLineStyle(7);
			graphCombPi0InvYieldRelSys2760GeVB->Draw("l,x0,same,e1");

			TLegend* legendRelTotErr2 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
			legendRelTotErr2->AddEntry(graphCombPi0InvYieldRelTot2760GeVB,"tot","p");
			legendRelTotErr2->AddEntry(graphCombPi0InvYieldRelStat2760GeVB,"stat","l");
			legendRelTotErr2->AddEntry(graphCombPi0InvYieldRelSys2760GeVB,"sys","l");
			legendRelTotErr2->Draw();

			labelRelTotErrEnergy->Draw();
			labelRelTotErrPi0->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelMethodBdecomp.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		histo2DRelTotErr->GetYaxis()->SetRangeUser(0,80.5);
		histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
			histo2DRelTotErr->Draw("copy");	

			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot2760GeVA, markerStyleComb, markerSizeComb, colorComb , colorComb);		
			graphCombPi0InvYieldRelTot2760GeVA->Draw("p,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelStat2760GeVA, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
			graphCombPi0InvYieldRelStat2760GeVA->Draw("l,x0,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelSys2760GeVA, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
			graphCombPi0InvYieldRelSys2760GeVA->SetLineStyle(7);
			graphCombPi0InvYieldRelSys2760GeVA->Draw("l,x0,same,e1");

			TLegend* legendRelTotErr3 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
			legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelTot2760GeVA,"tot","p");
			legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelStat2760GeVA,"stat","l");
			legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelSys2760GeVA,"sys","l");
			legendRelTotErr3->Draw();

			labelRelTotErrEnergy->Draw();
			labelRelTotErrPi0->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelMethodAdecomp.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		
		
		//	**********************************************************************************************************************
		//	************************************* Calculating bin shifted spectra & fitting **************************************
		//	**********************************************************************************************************************
		
		// Cloning spectra
		TGraphAsymmErrors* graphCombPi0InvYieldTot2760GeVAUnShifted 		= (TGraphAsymmErrors*)graphCombPi0InvYieldTot2760GeVA->Clone("Unshifted"); 
		TGraphAsymmErrors* graphCombPi0InvYieldStat2760GeVAUnShifted 		= (TGraphAsymmErrors*)graphCombPi0InvYieldStat2760GeVA->Clone("UnshiftedStat"); 
		TGraphAsymmErrors* graphCombPi0InvYieldSys2760GeVAUnShifted 		= (TGraphAsymmErrors*)graphCombPi0InvYieldSys2760GeVA->Clone("UnshiftedSys"); 

		TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeVUnShifted 		= (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV->Clone("UnshiftedStatPCM"); 
		TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeVUnShifted 		= (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV->Clone("UnshiftedSysPCM"); 

		TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeVUnshifted 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV->Clone("UnshiftedStatEMCal"); 
		TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeVUnshifted 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV->Clone("UnshiftedSysEMCal"); 
	
		// Calculating binshifts
		TF1* fitBylinkinPi0PbPb2760GeVPt = FitObject("tcm","BylinkinFitPi0","Pi0",graphCombPi0InvYieldStat2760GeVA,0.4,30.);
		graphCombPi0InvYieldStat2760GeVA->Fit(fitBylinkinPi0PbPb2760GeVPt,"QNRMEX0+","",0.4,30.);
	//	graphPCMPi0InvYieldStatPbPb2760GeV->Fit(fitBylinkinPi0PbPb2760GeVPt,"QNRMEX0+","",0.4,30.);

		if(bWCorrection.CompareTo("X")==0 ){
			TF1* fitBylinkinPi0PbPb2760GeVPtMult = FitObject("tcm","BylinkinFitPi0","Pi0");		
			graphCombPi0InvYieldTot2760GeVA			= ApplyXshift(graphCombPi0InvYieldTot2760GeVA, fitBylinkinPi0PbPb2760GeVPtMult,"Pi0");
			
			graphCombPi0InvYieldStat2760GeVA 		= ApplyXshiftIndividualSpectra (graphCombPi0InvYieldTot2760GeVA, 
																							graphCombPi0InvYieldStat2760GeVA, 
																							fitBylinkinPi0PbPb2760GeVPtMult,
																							0, 5,"Pi0");
			
			graphCombPi0InvYieldSys2760GeVA 			= ApplyXshiftIndividualSpectra (graphCombPi0InvYieldTot2760GeVA, 
																							graphCombPi0InvYieldSys2760GeVA, 
																							fitBylinkinPi0PbPb2760GeVPtMult, 
																							0, 5,"Pi0");
			
			graphPCMPi0InvYieldStatPbPb2760GeV			= ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot2760GeVA,
																							graphPCMPi0InvYieldStatPbPb2760GeV,
																							fitBylinkinPi0PbPb2760GeVPtMult, 
																							0, 5,"Pi0");
			
			graphPCMPi0InvYieldSysPbPb2760GeV			= ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot2760GeVA, 
																							graphPCMPi0InvYieldSysPbPb2760GeV, 
																							fitBylinkinPi0PbPb2760GeVPtMult, 
																							0, 5,"Pi0");
			
			graphEMCalPi0InvYieldStatPbPb2760GeV 		= ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot2760GeVA, 
																							graphEMCalPi0InvYieldStatPbPb2760GeV, 
																							fitBylinkinPi0PbPb2760GeVPtMult,
																							15, 5,"Pi0");	
			
			graphEMCalPi0InvYieldSysPbPb2760GeV 		= ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot2760GeVA, 
																							graphEMCalPi0InvYieldSysPbPb2760GeV, 
																							fitBylinkinPi0PbPb2760GeVPtMult,
																							15, 5,"Pi0");	

					
			TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
			DrawGammaCanvasSettings( canvasDummy2,  0.1, 0.01, 0.015, 0.08);
			canvasDummy2->SetLogy();
			canvasDummy2->SetLogx();
			TH2F * histo2DDummy2;
			histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,40.,1000,1e-7,1e4);
			SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
			histo2DDummy2->DrawCopy(); 
		
			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat2760GeVAUnShifted, 20, 3, kRed, kRed, widthLinesBoxes, kTRUE);
			graphCombPi0InvYieldStat2760GeVAUnShifted->Draw("pEsame");
			DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldTot2760GeVA, 24, 3, kBlack, kBlack, widthLinesBoxes, kTRUE);
			graphCombPi0InvYieldTot2760GeVA->Draw("pEsame");
			DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldSysPbPb2760GeV, markerStyleDet[0] ,/*markerSizeDet[0]*/2, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
			graphPCMPi0InvYieldSysPbPb2760GeV->Draw("pEsame");
			DrawGammaSetMarkerTGraphAsym(graphEMCalPi0InvYieldSysPbPb2760GeV, markerStyleDet[2] ,/*markerSizeDet[2]*/2.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
			graphEMCalPi0InvYieldSysPbPb2760GeV->Draw("pEsame");
						
			TLegend* legendXdummy = new TLegend(0.6,0.8,0.9,0.95);
			legendXdummy->SetFillColor(0);
			legendXdummy->SetLineColor(0);
			legendXdummy->SetTextFont(42);
			legendXdummy->SetTextSize(0.035);
			legendXdummy->AddEntry(graphCombPi0InvYieldStat2760GeVAUnShifted,"combined unshifted","lp");
// 			legendXdummy->AddEntry(graphCombPi0InvYieldTot2760GeVA,"combined unshifted (stat. err.)","lp");
			legendXdummy->AddEntry(graphPCMPi0InvYieldSysPbPb2760GeV,"PCM","fp");
			legendXdummy->AddEntry(graphEMCalPi0InvYieldSysPbPb2760GeV,"EMCal","fp");
			legendXdummy->Draw();
			
// 			TLatex *labelRelTotErrEnergy = new TLatex(0.75,0.89,collisionSystem2760GeV.Data());
// 			SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel,4);
// 			labelRelTotErrEnergy->SetTextFont(43);
// 			labelRelTotErrEnergy->Draw();
// 			TLatex *labelRelTotErrEta = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
// 			SetStyleTLatex( labelRelTotErrEta, 0.85*textSizeLabelsPixel,4);
// 			labelRelTotErrEta->SetTextFont(43);
// 			labelRelTotErrEta->Draw();
// 			labelRelTotErrEnergy->Draw();
// 			labelRelTotErrPi0->Draw();

			canvasDummy2->Update();
			canvasDummy2->Print(Form("%s/%s_%s_ComparisonShiftedPi0_PbPb2760GeV.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		}
			
		
		Double_t paramTCM[5] = {graphCombPi0InvYieldTot2760GeVA->GetY()[0],0.3,graphCombPi0InvYieldTot2760GeVA->GetY()[0]/10000,0.8,3};
		TF1* fitTCMInvYieldPi02760GeV = FitObject("tcm","fitTCMInvYieldPi02760GeV","Pi0",graphCombPi0InvYieldTot2760GeVA,0.4,50.,paramTCM,"QNRMEX0+");
		fitTCMInvYieldPi02760GeV = FitObject("tcm","fitTCMInvYieldPi02760GeV","Pi0",graphCombPi0InvYieldTot2760GeVA,0.4,50. ,paramTCM,"QNRMEX0+");

		cout << WriteParameterToFile(fitTCMInvYieldPi02760GeV)<< endl;
		
		TGraphAsymmErrors* graphRatioCombCombFitTot2760GeVA 	= (TGraphAsymmErrors*)graphCombPi0InvYieldTot2760GeVA->Clone();
		graphRatioCombCombFitTot2760GeVA 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitTot2760GeVA, fitTCMInvYieldPi02760GeV); 
		TGraphAsymmErrors* graphRatioCombCombFitStat2760GeVA 	= (TGraphAsymmErrors*)graphCombPi0InvYieldStat2760GeVA->Clone();
		graphRatioCombCombFitStat2760GeVA 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitStat2760GeVA, fitTCMInvYieldPi02760GeV); 
		TGraphAsymmErrors* graphRatioCombCombFitSys2760GeVA 	= (TGraphAsymmErrors*)graphCombPi0InvYieldSys2760GeVA->Clone();
		graphRatioCombCombFitSys2760GeVA 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitSys2760GeVA, fitTCMInvYieldPi02760GeV); 

		TGraphAsymmErrors* graphRatioPCMCombFitStat2760GeV 		= (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV->Clone();
		graphRatioPCMCombFitStat2760GeV 						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat2760GeV, fitTCMInvYieldPi02760GeV); 
		TGraphAsymmErrors* graphRatioPCMCombFitSys2760GeV 		= (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV->Clone();
		graphRatioPCMCombFitSys2760GeV 							= CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys2760GeV, fitTCMInvYieldPi02760GeV); 
		TGraphAsymmErrors* graphRatioEMCalCombFitStat2760GeV 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV->Clone();
		graphRatioEMCalCombFitStat2760GeV 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitStat2760GeV, fitTCMInvYieldPi02760GeV); 
		TGraphAsymmErrors* graphRatioEMCalCombFitSys2760GeV 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV->Clone();
		graphRatioEMCalCombFitSys2760GeV 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitSys2760GeV, fitTCMInvYieldPi02760GeV); 
		

		if(centrality.Contains("2050")){

			graphPCMPi0InvYieldSysPbPb2760GeV->RemovePoint(0);
			graphPCMPi0InvYieldStatPbPb2760GeV->RemovePoint(0);
			graphCombPi0InvYieldSys2760GeVA->RemovePoint(0);
			graphCombPi0InvYieldStat2760GeVA->RemovePoint(0);
			
		}

		//	**********************************************************************************************************************
		//	******************************************* Ratio of Comb to Fit ****************************************
		//	**********************************************************************************************************************
		textSizeLabelsPixel = 48;
		TCanvas* canvasRatioToCombFit = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRatioToCombFit, 0.12, 0.01, 0.01, 0.11);
		canvasRatioToCombFit->SetLogx();

			Double_t textsizeLabelsPP = 0;
			Double_t textsizeFacPP= 0;
			if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
				textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
				textsizeFacPP = (Double_t)1./canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
			} else {
				textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
				textsizeFacPP = (Double_t)1./canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
			}
			cout << textsizeLabelsPP << endl;
		
		TH2F * histo2DRatioToCombFit;
		histo2DRatioToCombFit = new TH2F("histo2DRatioToCombFit","histo2DRatioToCombFit",1000,0.23,70.,1000,0.2,4.	);
		SetStyleHistoTH2ForGraphs(histo2DRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
								0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
		histo2DRatioToCombFit->GetXaxis()->SetMoreLogLabels();
		histo2DRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
		histo2DRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
		histo2DRatioToCombFit->Draw("copy");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys2760GeVA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
		graphRatioCombCombFitSys2760GeVA->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat2760GeVA, markerStyleComb, markerSizeComb, colorComb , colorComb);
		graphRatioCombCombFitStat2760GeVA->Draw("p,same,e1");

		DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
		DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);

		TLatex *labelRatioToFitEnergy = new TLatex(0.73,0.92,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelRatioToFitEnergy, 0.85*textSizeLabelsPixel,4);
		labelRatioToFitEnergy->SetTextFont(43);
		labelRatioToFitEnergy->Draw();
		TLatex *labelRatioToFitPi0 = new TLatex(0.73,0.87,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRatioToFitPi0, 0.85*textSizeLabelsPixel,4);
		labelRatioToFitPi0->SetTextFont(43);
		labelRatioToFitPi0->Draw();

			
		canvasRatioToCombFit->SaveAs(Form("%s/%s_%s_RatioOfCombToCombFit_PbPb2760GeV.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		//	**********************************************************************************************************************
		//	******************************************* Ratio of Individual meas to Fit ******************************************
		//	**********************************************************************************************************************
		
		canvasRatioToCombFit->cd();
		histo2DRatioToCombFit->Draw("copy");
	
		DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitSys2760GeV, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitStat2760GeV, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
		DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitSys2760GeV, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitStat2760GeV, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
		
		graphRatioPCMCombFitSys2760GeV->Draw("E2same");
		graphRatioEMCalCombFitSys2760GeV->Draw("E2same");
		
		graphRatioPCMCombFitStat2760GeV->Draw("p,same,e");
		graphRatioEMCalCombFitStat2760GeV->Draw("p,same,e");

		DrawGammaLines(0.23, 70. , 1., 1.,0.5, kGray+2);
		DrawGammaLines(0.23, 70. , 1.1, 1.1,0.5, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.9, 0.9,0.5, kGray, 7);
		
		labelRatioToFitEnergy->Draw();
		labelRatioToFitPi0->Draw();
	
		// ****************************** Definition of the Legend ******************************************
		// **************** Row def ************************
		Double_t rowsLegendOnlyPi0Ratio[5] 		= {0.92,0.88,0.84,0.80,0.76};
		Double_t rowsLegendOnlyPi0RatioAbs[5] 	= {0.91,2.2,2.1,2.0,1.9};
		Double_t columnsLegendOnlyPi0Ratio[3] 	= {0.15,0.32, 0.38};
		Double_t columnsLegendOnlyPi0RatioAbs[3]= {0.15,1.04, 1.37};
		Double_t lengthBox						= 0.2/2;
		Double_t heightBox						= 0.08/2;
		// ****************** first Column **************************************************
		TLatex *textPCMOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],"PCM");
		SetStyleTLatex( textPCMOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
		textPCMOnlyRatioPi0->SetTextFont(43);
		textPCMOnlyRatioPi0->Draw();
		TLatex *textEMCalOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],"EMCal");
		SetStyleTLatex( textEMCalOnlyRatioPi0,  0.85*textSizeLabelsPixel,4);
		textEMCalOnlyRatioPi0->SetTextFont(43);
		textEMCalOnlyRatioPi0->Draw();
		
		// ****************** second Column *************************************************
		TLatex *textStatOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
		SetStyleTLatex( textStatOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
		textStatOnlyRatioPi0->SetTextFont(43);
		textStatOnlyRatioPi0->Draw();
		TLatex *textSysOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
		SetStyleTLatex( textSysOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
		textSysOnlyRatioPi0->SetTextFont(43);
		textSysOnlyRatioPi0->Draw();
		TMarker* markerPCMPi0OnlyRatioPi0 = CreateMarkerFromGraph(graphRatioPCMCombFitSys2760GeV,columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
		markerPCMPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
		TMarker* markerEMCalPi0OnlyRatioPi0 = CreateMarkerFromGraph(graphRatioEMCalCombFitSys2760GeV, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
		markerEMCalPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);

		TBox* boxPCMPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioPCMCombFitSys2760GeV, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
														 columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
		boxPCMPi0OnlyRatioPi0->Draw("l");
		TBox* boxEMCalPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioEMCalCombFitSys2760GeV, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
														   columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
		boxEMCalPi0OnlyRatioPi0->Draw("l");
		
		canvasRatioToCombFit->SaveAs(Form("%s/%s_%s_RatioOfIndividualMeasToCombFit_PbPb2760GeV.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

	
	
	
	

		//	**********************************************************************************************************************
		//	******************************** Cross section for pi0 single measurement 2.76TeV ************************************
		//	**********************************************************************************************************************
		
		TCanvas* canvasXSectionPi0 = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
		DrawGammaCanvasSettings( canvasXSectionPi0, 0.14, 0.02, 0.02, 0.09);
		canvasXSectionPi0->SetLogx();
		canvasXSectionPi0->SetLogy();
		
		TH2F * histo2DXSectionPi0;
		histo2DXSectionPi0 = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,70.,1000,2e-7,1e3);
		SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 1.,1.45);
		histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
		histo2DXSectionPi0->GetXaxis()->SetLabelOffset(-0.01);
		histo2DXSectionPi0->Draw("copy");


		DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldSysPbPb2760GeV, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
		graphPCMPi0InvYieldSysPbPb2760GeV->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphEMCalPi0InvYieldSysPbPb2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
		graphEMCalPi0InvYieldSysPbPb2760GeV->Draw("E2same");
	
		DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldStatPbPb2760GeV,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
		graphPCMPi0InvYieldStatPbPb2760GeV->Draw("p,same,e1");
			
		DrawGammaSetMarkerTGraphAsym(graphEMCalPi0InvYieldStatPbPb2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
		graphEMCalPi0InvYieldStatPbPb2760GeV->Draw("p,same,e1");
		
// 		fitBylinkinPi0PbPb2760GeVPt->SetRange(0.4,30.);
		fitBylinkinPi0PbPb2760GeVPt->Draw("same");

		TLatex *labelEnergyXSectionPi0 = new TLatex(0.5,0.92,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4);
		labelEnergyXSectionPi0->Draw();
		TLatex *labelDetSysXSectionPi0 = new TLatex(0.5,0.88,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4);
		labelDetSysXSectionPi0->Draw();

		TLegend* legendXSectionPi0 = new TLegend(0.62,0.66,0.9,0.86);
		legendXSectionPi0->SetFillColor(0);
		legendXSectionPi0->SetLineColor(0);
		legendXSectionPi0->SetTextFont(42);
		legendXSectionPi0->SetTextSize(0.035);
		legendXSectionPi0->AddEntry(graphPCMPi0InvYieldSysPbPb2760GeV,"PCM","fp");
		legendXSectionPi0->AddEntry(graphEMCalPi0InvYieldStatPbPb2760GeV,"EMCal","fp");
// 		legendXSectionPi0->AddEntry(fitBylinkinPi0PbPb2760GeVPt, "Bylinkin Fit","l");
		legendXSectionPi0->Draw();
   
		canvasXSectionPi0->SaveAs(Form("%s/%s_%s_YieldsPi0CompAllSystems.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		canvasXSectionPi0->cd();

		histo2DXSectionPi0->Draw("copy");

		graphPCMPi0InvYieldSysPbPb2760GeV->Draw("E2same");
		graphPCMPi0InvYieldStatPbPb2760GeV->Draw("p,same,e1");
		
		graphEMCalPi0InvYieldSysPbPb2760GeV->Draw("E2same");
		graphEMCalPi0InvYieldStatPbPb2760GeV->Draw("p,same,e1");
		
		
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys2760GeVA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat2760GeVA, markerStyleComb, markerSizeComb, colorComb , colorComb);
		graphCombPi0InvYieldSys2760GeVA->Draw("E2same");
		graphCombPi0InvYieldStat2760GeVA->Draw("p,same,e1");
	
		
		fitBylinkinPi0PbPb2760GeVPt->Draw("same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();

		legendXSectionPi0->AddEntry(graphCombPi0InvYieldSys2760GeVA,"comb","fp");
		legendXSectionPi0->AddEntry(fitBylinkinPi0PbPb2760GeVPt, "Bylinkin Fit","l");
		legendXSectionPi0->Draw();

// 		DrawGammaSetMarkerTGraphAsym(graphChargedHadronsStatPP2760GeV, markerStyleDet[1], markerSizeDet[1]*0.75, kBlack , kBlack, widthLinesBoxes);
// 		graphChargedHadronsStatPP2760GeV->Draw("E2same");
		
		canvasXSectionPi0->SaveAs(Form("%s/%s_%s_YieldPi0CompAllSystems_Comb.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
  		
		
		//	**********************************************************************************************************************
		//	*************************	 Combination of RAA 	*******************************************************************
	
// 		Int_t bin0PCM =0;
// 		Int_t bin0EMCal =0;
// 		
// 		TGraphAsymmErrors* graphCombPi0RAAStat2760GeVOEMC= NULL;
// 		TGraphAsymmErrors* graphCombPi0RAASys2760GeVOEMC = NULL;
// 		TGraphAsymmErrors* graphCombPi0RAATot2760GeVOEMC = CombinePtPointsRAA( 	graphPCMPi0RAAStat2760GeV,	graphPCMPi0RAASys2760GeV, 	
// 																				graphEMCalPi0RAAStat2760GeV, graphEMCalPi0RAASys2760GeV,																graphCombPi0RAAStat2760GeVOEMC, graphCombPi0RAASys2760GeVOEMC,
// 																				xPtLimitsPi0, 23, offSetsPi0, bin0PCM, bin0EMCal, 
// 																				kFALSE,1);
// 		
// // TGraphAsymmErrors* CombinePtPointsRAA( TGraphAsymmErrors* graphStatErrPCM,		TGraphAsymmErrors* graphSystPCM,
// // 										  TGraphAsymmErrors* graphStatErrPHOS,		TGraphAsymmErrors* graphSystPHOS,
// //  									  TGraphAsymmErrors* &graphStatComb, 	TGraphAsymmErrors* &graphSystComb,  
// // 													Double_t* xPtLimits,	Int_t nPtLimits,
// // 													Int_t offset,			Int_t bin0PCM,	 Int_t bin0PHOS, Bool_t kRemoveLastPCMPoint=kFALSE, Int_t nBinsPCMRem = 1){
// 		
// // 		cout << "Printing of PCM" << endl;
// // 		graphPCMPi0RAAStat2760GeV->Print();
// // 		graphPCMPi0RAASys2760GeV->Print();
// // 		
// // 		cout << "Printing of EMCal" << endl;
// // 		graphEMCalPi0RAAStat2760GeV->Print();
// // 		graphEMCalPi0RAASys2760GeV->Print();
// 
// 		cout << "Printing of combined RAA" << endl;
// // 		graphCombPi0RAAStat2760GeVOEMC->RemovePoint(0);
// //         graphCombPi0RAASys2760GeVOEMC->RemovePoint(0);
// 		graphCombPi0RAAStat2760GeVOEMC->Print();
// 		graphCombPi0RAASys2760GeVOEMC->Print();
// 

		TCanvas* canvasRAAcombo = new TCanvas("canvasRAAcombo","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRAAcombo, 0.08, 0.02, 0.035, 0.09);
// 		canvasRAAcombo->SetLogx();
	
		TH2F * histo2DRAAcombo;
		histo2DRAAcombo = new TH2F("histo2DRAAcombo","histo2DRAAcombo",11000,0.23,70.,1000,-0.5,1.1);
		SetStyleHistoTH2ForGraphs(histo2DRAAcombo, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,1.);
// 		histo2DRAAcombo->GetXaxis()->SetMoreLogLabels();	
// 		histo2DRAAcombo->GetXaxis()->SetLabelOffset(-0.01);
		histo2DRAAcombo->GetYaxis()->SetRangeUser(0.,1.);
		histo2DRAAcombo->GetXaxis()->SetRangeUser(0.,30.);
		histo2DRAAcombo->Draw("copy");
		
		DrawGammaSetMarkerTGraphAsym(graphPCMPi0RAASys2760GeV, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
		graphPCMPi0RAASys2760GeV->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPCMPi0RAAStat2760GeV,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
		graphPCMPi0RAAStat2760GeV->Draw("p,same,e1");		
		
// 		DrawGammaSetMarkerTGraphAsym(graphEMCalPi0RAASys2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
// 		graphEMCalPi0RAASys2760GeV->Draw("E2same");			
// 		DrawGammaSetMarkerTGraphAsym(graphEMCalPi0RAAStat2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
// 		graphEMCalPi0RAAStat2760GeV->Draw("p,same,e1");
		
// 		DrawGammaSetMarkerTGraphAsym(graphCombPi0RAASys2760GeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
// 		graphCombPi0RAASys2760GeVOEMC->Draw("E2same");
// 		DrawGammaSetMarkerTGraphAsym(graphCombPi0RAAStat2760GeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb);
// 		graphCombPi0RAAStat2760GeVOEMC->Draw("p,same,e1");

		TLegend* legendRAAcombo = new TLegend(0.62,0.66,0.9,0.86);
		legendRAAcombo->SetFillColor(0);
		legendRAAcombo->SetLineColor(0);
		legendRAAcombo->SetTextFont(42);
		legendRAAcombo->SetTextSize(0.035);
// 		legendRAAcombo->AddEntry(graphCombPi0RAAStat2760GeVOEMC,"combined","p");		
		legendRAAcombo->AddEntry(graphPCMPi0RAASys2760GeV,"PCM","p");
// 		legendRAAcombo->AddEntry(graphEMCalPi0RAASys2760GeV,"EMCal","p");
		legendRAAcombo->Draw();

		TLatex *labelRAAEnergy = new TLatex(0.7,0.20,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelRAAEnergy, 0.85*textSizeLabelsPixel,4);
		labelRAAEnergy->SetTextFont(43);
// 		labelRAAEnergy->Draw();
		TLatex *labelRAAPi0 = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRAAPi0, 0.85*textSizeLabelsPixel,4);
		labelRAAPi0->SetTextFont(43);
// 		labelRAAPi0->Draw();

		canvasRAAcombo->SaveAs(Form("%s/%s_%s_RAAPi0_combined.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		//	**********************************************************************************************************************
		//	*************************	 Combination of RCP 	*******************************************************************
		cout << "\n\n\n *************************	 Combination of RCP 	************************* \n\n\n" << endl;
		Int_t bin0PCMRCP =0;
		Int_t bin0EMCalRCP =16;
		graphPCMPi0RCPStat2760GeV->SetName("Pi0RCP_pcm");
		graphEMCalPi0RCPStat2760GeV->RemovePoint(0);
		graphEMCalPi0RCPSys2760GeV->RemovePoint(0);
		TGraphAsymmErrors* graphCombPi0RCPStat2760GeVOEMC= NULL;
		TGraphAsymmErrors* graphCombPi0RCPSys2760GeVOEMC = NULL;
		TGraphAsymmErrors* graphCombPi0RCPTot2760GeVOEMC = CombinePtPointsRAA( 	graphPCMPi0RCPStat2760GeV,	graphPCMPi0RCPSys2760GeV, 	
																				graphEMCalPi0RCPStat2760GeV, graphEMCalPi0RCPSys2760GeV,																graphCombPi0RCPStat2760GeVOEMC, graphCombPi0RCPSys2760GeVOEMC,
																				xPtLimitsPi0, 24, /*offSetsPi0*/0, bin0PCMRCP, bin0EMCalRCP);
		
// TGraphAsymmErrors* CombinePtPointsRCP( TGraphAsymmErrors* graphStatErrPCM,		TGraphAsymmErrors* graphSystPCM,
// 										  TGraphAsymmErrors* graphStatErrPHOS,		TGraphAsymmErrors* graphSystPHOS,
//  									  TGraphAsymmErrors* &graphStatComb, 	TGraphAsymmErrors* &graphSystComb,  
// 													Double_t* xPtLimits,	Int_t nPtLimits,
// 													Int_t offset,			Int_t bin0PCM,	 Int_t bin0PHOS, Bool_t kRemoveLastPCMPoint=kFALSE, Int_t nBinsPCMRem = 1){
		
// 		cout << "Printing of PCM" << endl;
// 		graphPCMPi0RCPStat2760GeV->Print();
// 		graphPCMPi0RCPSys2760GeV->Print();
// 		
		cout << "Printing of EMCal" << endl;
		graphEMCalPi0RCPStat2760GeV->Print();
		graphEMCalPi0RCPSys2760GeV->Print();

		cout << "Printing of combined RCP" << endl;
// 		graphCombPi0RCPStat2760GeVOEMC->RemovePoint(0);
//         graphCombPi0RCPSys2760GeVOEMC->RemovePoint(0);
		graphCombPi0RCPStat2760GeVOEMC->Print();
		graphCombPi0RCPSys2760GeVOEMC->Print();


		TCanvas* canvasRCPcombo = new TCanvas("canvasRCPcombo","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRCPcombo, 0.08, 0.02, 0.035, 0.09);
// 		canvasRCPcombo->SetLogx();
	
		TH2F * histo2DRCPcombo;
		histo2DRCPcombo = new TH2F("histo2DRCPcombo","histo2DRCPcombo",11000,0.23,70.,1000,-0.5,1.1);
		SetStyleHistoTH2ForGraphs(histo2DRCPcombo, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,1.);
// 		histo2DRCPcombo->GetXaxis()->SetMoreLogLabels();	
// 		histo2DRCPcombo->GetXaxis()->SetLabelOffset(-0.01);
		histo2DRCPcombo->GetYaxis()->SetRangeUser(0.,1.);
		histo2DRCPcombo->GetXaxis()->SetRangeUser(0.,30.);
		histo2DRCPcombo->Draw("copy");
		
		DrawGammaSetMarkerTGraphAsym(graphPCMPi0RCPSys2760GeV, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
		graphPCMPi0RCPSys2760GeV->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPCMPi0RCPStat2760GeV,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
		graphPCMPi0RCPStat2760GeV->Draw("p,same,e1");		
		
		DrawGammaSetMarkerTGraphAsym(graphEMCalPi0RCPSys2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
		graphEMCalPi0RCPSys2760GeV->Draw("E2same");			
		DrawGammaSetMarkerTGraphAsym(graphEMCalPi0RCPStat2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
		graphEMCalPi0RCPStat2760GeV->Draw("p,same,e1");
		
		DrawGammaSetMarkerTGraphAsym(graphCombPi0RCPSys2760GeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
		graphCombPi0RCPSys2760GeVOEMC->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombPi0RCPStat2760GeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb);
		graphCombPi0RCPStat2760GeVOEMC->Draw("p,same,e1");

		TLegend* legendRCPcombo = new TLegend(0.62,0.66,0.9,0.86);
		legendRCPcombo->SetFillColor(0);
		legendRCPcombo->SetLineColor(0);
		legendRCPcombo->SetTextFont(42);
		legendRCPcombo->SetTextSize(0.035);
		legendRCPcombo->AddEntry(graphCombPi0RCPStat2760GeVOEMC,"combined","p");		
		legendRCPcombo->AddEntry(graphPCMPi0RCPSys2760GeV,"PCM","p");
		legendRCPcombo->AddEntry(graphEMCalPi0RCPSys2760GeV,"EMCal","p");
		legendRCPcombo->Draw();

		TLatex *labelRCPEnergy = new TLatex(0.7,0.20,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelRCPEnergy, 0.85*textSizeLabelsPixel,4);
		labelRCPEnergy->SetTextFont(43);
// 		labelRCPEnergy->Draw();
		TLatex *labelRCPPi0 = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRCPPi0, 0.85*textSizeLabelsPixel,4);
		labelRCPPi0->SetTextFont(43);
// 		labelRCPPi0->Draw();

		canvasRCPcombo->SaveAs(Form("%s/%s_%s_RCPPi0_combined.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		

		
		//	**********************************************************************************************************************
		//	************************* Saving of final results ********************************************************************
		//	**********************************************************************************************************************
		
		cout <<"method A" << endl;
		graphCombPi0InvYieldStat2760GeVA->Print();
		cout <<"method B" << endl;
		graphCombPi0InvYieldStat2760GeVB->Print();
		TFile fCombResults(Form("%s/%s_%s_CombinedResultsPaperPbPb2760GeV_%s.root", outputDir.Data(),meson.Data(),centrality.Data(),dateForOutput.Data()), "RECREATE");
			
			graphPCMPi0CorrYieldSysErrPbPb2760GeV_0010->Write("graphPCMPi0CorrYieldSysErrPbPb2760GeV_0010");
			graphPCMPi0InvYieldStatPbPb2760GeV_0010->Write("graphPCMPi0InvYieldStatPbPb2760GeV_0010");
			graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010->Write("graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010");
			graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Write("graphEMCalPi0InvYieldStatPbPb2760GeV_0010");

			graphCombPi0InvYieldTot2760GeVA->Write("graphInvYieldPi0Comb2760GeVA");
			graphCombPi0InvYieldStat2760GeVA->Write("graphInvYieldPi0Comb2760GeVAStatErr");
			graphCombPi0InvYieldSys2760GeVA->Write("graphInvYieldPi0Comb2760GeVASysErr");      
			graphCombPi0InvYieldTot2760GeVB->Write("graphInvYieldPi0Comb2760GeVB");
			graphCombPi0InvYieldStat2760GeVB->Write("graphInvYieldPi0Comb2760GeVBStatErr");
			graphCombPi0InvYieldSys2760GeVB->Write("graphInvYieldPi0Comb2760GeVBSysErr");      
		fCombResults.Close();
		
	} else if(meson.CompareTo("Eta")==0) {
				
// 		TString fileNameOutputWeightingOld							= Form("%s/%s_WeightingOldEta.dat",outputDir.Data(),centrality.Data());
// 		TGraphAsymmErrors* graphCombEtaInvYieldStat2760GeVOld= NULL;
// 		TGraphAsymmErrors* graphCombEtaInvYieldSys2760GeVOld = NULL;
// 		TGraphAsymmErrors* graphCombEtaInvYieldTot2760GeVOld = CombinePtPointsSpectraFullCorrMat( 	statErrorCollection,	sysErrorCollection, 	
// 																											xPtLimitsEta, 16,
// 																											offSetsEta, offSetsEtaSys,
// 																											graphCombEtaInvYieldStat2760GeVOld, graphCombEtaInvYieldSys2760GeVOld,
// 																											fileNameOutputWeightingOld,1
// 																										);
// 		cout << "Printing the graph:    " << endl;
// 		graphCombEtaInvYieldStat2760GeVOld->Print();
// 
// 		// Reading weights from output file for plotting
// 		ifstream fileWeightsReadOld;
// 		fileWeightsReadOld.open(fileNameOutputWeightingOld,ios_base::in);
// 		cout << "reading" << fileNameOutputWeightingOld << endl;
// 		Double_t xValuesReadOld[50];
// 		Double_t weightsReadOld[11][50];
// 		Int_t availableMeasOld[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
// 		Int_t nMeasSetOld 			= 2;
// 		Int_t nPtBinsReadOld 		= 0;
// 		while(!fileWeightsReadOld.eof() && nPtBinsReadOld < 50){
// 			TString garbage = "";
// 			if (nPtBinsReadOld == 0){
// 				fileWeightsReadOld >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
// 				for (Int_t i = 0; i < nMeasSetOld; i++){
// 					fileWeightsReadOld >> availableMeasOld[i] ;
// 				}	
// 				cout << "read following measurements: "; 
// 				for (Int_t i = 0; i < 11; i++){
// 					cout << availableMeasOld[i] << "\t" ;
// 				}	
// 				cout << endl;
// 			} else {
// 				fileWeightsReadOld >> xValuesReadOld[nPtBinsReadOld-1];
// 				for (Int_t i = 0; i < nMeasSetOld; i++){
// 					fileWeightsReadOld >> weightsReadOld[availableMeasOld[i]][nPtBinsReadOld-1] ;
// 				}	
// 				cout << "read: "<<  nPtBinsReadOld << "\t"<< xValuesReadOld[nPtBinsReadOld-1] << "\t" ;
// 				for (Int_t i = 0; i < nMeasSetOld; i++){
// 					cout << weightsReadOld[availableMeasOld[i]][nPtBinsReadOld-1] << "\t";
// 				}
// 				cout << endl;
// 			}
// 			nPtBinsReadOld++;
// 		}
// 		nPtBinsReadOld = nPtBinsReadOld-2 ;
// 		fileWeightsReadOld.close();
// 
// 		for (Int_t i = 0; i < nMeasSetOld; i++){
// 			graphWeightsOld[availableMeasOld[i]] = new TGraph(nPtBinsReadOld,xValuesReadOld,weightsReadOld[availableMeasOld[i]]);
// 			Int_t bin = 0;
// 			for (Int_t n = 0; n< nPtBinsReadOld; n++){
// 				if (graphWeightsOld[availableMeasOld[i]]->GetY()[bin] == 0) graphWeightsOld[availableMeasOld[i]]->RemovePoint(bin);
// 				else bin++;
// 			}	
// 		}	
	
	
		//	**********************************************************************************************************************
		//	******************************************* Plotting weights method only EMC *****************************************
		//	**********************************************************************************************************************
		Int_t textSizeLabelsPixel = 900*0.04;

		TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
		canvasWeights->SetLogx();
	
		TH2F * histo2DWeights;
		histo2DWeights = new TH2F("histo2DWeights","histo2DWeights",11000,0.23,70.,1000,-0.5,1.1);
		SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DWeights->GetXaxis()->SetMoreLogLabels();
		histo2DWeights->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DWeights->GetYaxis()->SetRangeUser(-10,10);
		histo2DWeights->Draw("copy");
// 		
// 		TLegend* legendAccWeightsOld = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetOld*1.35), 32);
// 		for (Int_t i = 0; i < nMeasSetOld; i++){
// 			DrawGammaSetMarkerTGraph(graphWeightsOld[availableMeasOld[i]],
// 									markerStyleDet[availableMeasOld[i]], 
// 									markerSizeDet[availableMeasOld[i]]*0.5, 
// 									colorDet[availableMeasOld[i]] , 
// 									colorDet[availableMeasOld[i]]);
// 			graphWeightsOld[availableMeasOld[i]]->Draw("p,same,e1");
// 			legendAccWeightsOld->AddEntry(graphWeightsOld[availableMeasOld[i]],nameMeasGlobal[availableMeasOld[i]],"p");
// 		}	
// 		legendAccWeightsOld->Draw();

		TLatex *labelWeightsEnergy = new TLatex(0.7,0.20,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
		labelWeightsEnergy->SetTextFont(43);
		labelWeightsEnergy->Draw();
		TLatex *labelWeightsEta = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
		SetStyleTLatex( labelWeightsEta, 0.85*textSizeLabelsPixel,4);
		labelWeightsEta->SetTextFont(43);
		labelWeightsEta->Draw();

// 	//	DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
// 		DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
// 		DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
// 		DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
// 		DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
// 			
// 		canvasWeights->SaveAs(Form("%s/%s_%s_WeightsOldPublished.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
// 
// 		
		//	**********************************************************************************************************************
		//	******************************************* Calculation of spectrum including EMCal only *****************************
		//	**********************************************************************************************************************
							
		TGraph* graphWeightsOEMC[11];
		for (Int_t i = 0; i< 11; i++){
			graphWeightsOEMC[i] = NULL;
		}	
		// Declaration & calculation of combined spectrum							
		TString fileNameOutputWeightingOEMC							= Form("%s/%s_WeightingEMCalEta.dat",outputDir.Data(),centrality.Data());
		TGraphAsymmErrors* graphCombEtaInvYieldStat2760GeVOEMC= NULL;
		TGraphAsymmErrors* graphCombEtaInvYieldSys2760GeVOEMC = NULL;
		TGraphAsymmErrors* graphCombEtaInvYieldTot2760GeVOEMC = CombinePtPointsSpectraFullCorrMat( 	statErrorCollection,	sysErrorCollection, 	
																											xPtLimitsEta, 13,
																											offSetsEta, offSetsEtaSys,
																											graphCombEtaInvYieldStat2760GeVOEMC, graphCombEtaInvYieldSys2760GeVOEMC,
																											fileNameOutputWeightingOEMC,1
																										);
		cout << "Printing of EMCal only spectra" << endl;
		graphCombEtaInvYieldStat2760GeVOEMC->RemovePoint(0);
        graphCombEtaInvYieldSys2760GeVOEMC->RemovePoint(0);
		graphCombEtaInvYieldStat2760GeVOEMC->Print();
		graphCombEtaInvYieldSys2760GeVOEMC->Print();

		// Reading weights from output file for plotting
		ifstream fileWeightsReadOEMC;
		fileWeightsReadOEMC.open(fileNameOutputWeightingOEMC,ios_base::in);
		cout << "reading" << fileNameOutputWeightingOEMC << endl;
		Double_t xValuesReadOEMC[50];
		Double_t weightsReadOEMC[11][50];
		Int_t availableMeasOEMC[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
		Int_t nMeasSetOEMC 			= 2;
		Int_t nPtBinsReadOEMC 		= 0;
		while(!fileWeightsReadOEMC.eof() && nPtBinsReadOEMC < 50){
			TString garbage = "";
			if (nPtBinsReadOEMC == 0){
				fileWeightsReadOEMC >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
				for (Int_t i = 0; i < nMeasSetOEMC; i++){
					fileWeightsReadOEMC >> availableMeasOEMC[i] ;
				}	
				cout << "read following measurements: "; 
				for (Int_t i = 0; i < 11; i++){
					cout << availableMeasOEMC[i] << "\t" ;
				}	
				cout << endl;
			} else {
				fileWeightsReadOEMC >> xValuesReadOEMC[nPtBinsReadOEMC-1];
				for (Int_t i = 0; i < nMeasSetOEMC; i++){
					fileWeightsReadOEMC >> weightsReadOEMC[availableMeasOEMC[i]][nPtBinsReadOEMC-1] ;
				}	
				cout << "read: "<<  nPtBinsReadOEMC << "\t"<< xValuesReadOEMC[nPtBinsReadOEMC-1] << "\t" ;
				for (Int_t i = 0; i < nMeasSetOEMC; i++){
					cout << weightsReadOEMC[availableMeasOEMC[i]][nPtBinsReadOEMC-1] << "\t";
				}
				cout << endl;
			}
			nPtBinsReadOEMC++;
		}
		nPtBinsReadOEMC = nPtBinsReadOEMC-2 ;
		fileWeightsReadOEMC.close();

		for (Int_t i = 0; i < nMeasSetOEMC; i++){
			graphWeightsOEMC[availableMeasOEMC[i]] = new TGraph(nPtBinsReadOEMC,xValuesReadOEMC,weightsReadOEMC[availableMeasOEMC[i]]);
			Int_t bin = 0;
			for (Int_t n = 0; n< nPtBinsReadOEMC; n++){
				if (graphWeightsOEMC[availableMeasOEMC[i]]->GetY()[bin] == 0) graphWeightsOEMC[availableMeasOEMC[i]]->RemovePoint(bin);
				else bin++;
			}	
		}	

		//	**********************************************************************************************************************
		//	******************************************* Plotting weights method only EMC *****************************************
		//	**********************************************************************************************************************
		canvasWeights->cd();
		histo2DWeights->Draw("copy");
		
		TLegend* legendAccWeightsOEMC = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetOEMC*1.35), 32);
		for (Int_t i = 0; i < nMeasSetOEMC; i++){
			DrawGammaSetMarkerTGraph(graphWeightsOEMC[availableMeasOEMC[i]],
									markerStyleDet[availableMeasOEMC[i]], 
									markerSizeDet[availableMeasOEMC[i]]*0.5, 
									colorDet[availableMeasOEMC[i]] , 
									colorDet[availableMeasOEMC[i]]);
			graphWeightsOEMC[availableMeasOEMC[i]]->Draw("p,same,e1");
			legendAccWeightsOEMC->AddEntry(graphWeightsOEMC[availableMeasOEMC[i]],nameMeasGlobal[availableMeasOEMC[i]],"p");
		}	
		legendAccWeightsOEMC->Draw();
		labelWeightsEnergy->Draw();
		labelWeightsEta->Draw();

	//	DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
		DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
		DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
			
		canvasWeights->SaveAs(Form("%s/%s_%s_WeightsOnlyEMCal.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		
		TH1D* statErrorRelCollection[11];
		for (Int_t i = 0; i< 11; i++){
			statErrorRelCollection[i] = NULL;
		}	
		for (Int_t i = 0; i < 11; i++){
			if (statErrorCollection[i]) statErrorRelCollection[i] = CalculateRelErrUpTH1D( statErrorCollection[i], Form("relativeStatError_%s", nameMeasGlobal[i].Data()));
		}
		
		TGraphAsymmErrors* sysErrorRelCollection[11];
		for (Int_t i = 0; i< 11; i++){
			sysErrorRelCollection[i] = NULL;
		}	
		for (Int_t i = 0; i < 11; i++){
			if (sysErrorCollection[i]) sysErrorRelCollection[i] = CalculateRelErrUpAsymmGraph( sysErrorCollection[i], Form("relativeSysError_%s", nameMeasGlobal[i].Data()));
		}

		
		//	**********************************************************************************************************************
		//	******************************************* Assuming maximal correlation *********************************************
		//	**********************************************************************************************************************

									
		TGraph* graphWeightsA[11];
		for (Int_t i = 0; i< 11; i++){
			graphWeightsA[i] = NULL;
		}	
									
		// Declaration & calculation of combined spectrum							
		TString fileNameOutputWeightingA 							= Form("%s/%s_WeightingMethodAEta.dat",outputDir.Data(),centrality.Data());
		TGraphAsymmErrors* graphCombEtaInvYieldStat2760GeVA 	= NULL;
		TGraphAsymmErrors* graphCombEtaInvYieldSys2760GeVA 	= NULL;
		TGraphAsymmErrors* graphCombEtaInvYieldTot2760GeVA 	= CombinePtPointsSpectraFullCorrMat( 	statErrorCollection,	sysErrorCollection, 	
																											xPtLimitsEta, 13,
																											offSetsEta, offSetsEtaSys,
																											graphCombEtaInvYieldStat2760GeVA, graphCombEtaInvYieldSys2760GeVA,
																											fileNameOutputWeightingA,1
																										);
		graphCombEtaInvYieldStat2760GeVA->RemovePoint(0);
		graphCombEtaInvYieldSys2760GeVA->RemovePoint(0);
		graphCombEtaInvYieldTot2760GeVA->RemovePoint(0);
		graphCombEtaInvYieldStat2760GeVA->RemovePoint(0);
		graphCombEtaInvYieldSys2760GeVA->RemovePoint(0);
		graphCombEtaInvYieldTot2760GeVA->RemovePoint(0);
		cout << "\n\n\n\n Printing spectra method A" << endl;
// 		graphCombEtaInvYieldStat2760GeVA->Print();
//         graphCombEtaInvYieldSys2760GeVA->Print();
//         graphCombEtaInvYieldTot2760GeVA->Print();
				
		// Reading weights from output file for plotting
		ifstream fileWeightsReadA;
		fileWeightsReadA.open(fileNameOutputWeightingA,ios_base::in);
		cout << "reading" << fileNameOutputWeightingA << endl;
		Double_t xValuesReadA[50];
		Double_t weightsReadA[11][50];
		Int_t availableMeasA[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
		Int_t nMeasSetA 			= 2;
		Int_t nPtBinsReadA 		= 0;
		while(!fileWeightsReadA.eof() && nPtBinsReadA < 50){
			TString garbage = "";
			if (nPtBinsReadA == 0){
				fileWeightsReadA >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
				for (Int_t i = 0; i < nMeasSetA; i++){
					fileWeightsReadA >> availableMeasA[i] ;
				}	
				cout << "read following measurements: "; 
				for (Int_t i = 0; i < 11; i++){
					cout << availableMeasA[i] << "\t" ;
				}	
				cout << endl;
			} else {
				fileWeightsReadA >> xValuesReadA[nPtBinsReadA-1];
				for (Int_t i = 0; i < nMeasSetA; i++){
					fileWeightsReadA >> weightsReadA[availableMeasA[i]][nPtBinsReadA-1] ;
				}	
				cout << "read: "<<  nPtBinsReadA << "\t"<< xValuesReadA[nPtBinsReadA-1] << "\t" ;
				for (Int_t i = 0; i < nMeasSetA; i++){
					cout << weightsReadA[availableMeasA[i]][nPtBinsReadA-1] << "\t";
				}
				cout << endl;
			}
			nPtBinsReadA++;
		}
		nPtBinsReadA = nPtBinsReadA-2 ;
		fileWeightsReadA.close();

		for (Int_t i = 0; i < nMeasSetA; i++){
			graphWeightsA[availableMeasA[i]] = new TGraph(nPtBinsReadA,xValuesReadA,weightsReadA[availableMeasA[i]]);
			Int_t bin = 0;
			for (Int_t n = 0; n< nPtBinsReadA; n++){
				if (graphWeightsA[availableMeasA[i]]->GetY()[bin] == 0) graphWeightsA[availableMeasA[i]]->RemovePoint(bin);
				else bin++;
			}	
		}	

		//	**********************************************************************************************************************
		//	******************************************* Plotting weights Method A ************************************************
		//	**********************************************************************************************************************
		canvasWeights->cd();
	
		histo2DWeights->Draw("copy");

		TLegend* legendAccWeights = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetA*1.35), 32);
		for (Int_t i = 0; i < nMeasSetA; i++){
			DrawGammaSetMarkerTGraph(graphWeightsA[availableMeasA[i]], markerStyleDet[availableMeasA[i]], markerSizeDet[availableMeasA[i]]*0.5, colorDet[availableMeasA[i]] , colorDet[availableMeasA[i]]);
			graphWeightsA[availableMeasA[i]]->Draw("p,same,e1");
			legendAccWeights->AddEntry(graphWeightsA[availableMeasA[i]],nameMeasGlobal[availableMeasA[i]],"p");
		}	
		legendAccWeights->Draw();

		labelWeightsEnergy->Draw();
		labelWeightsEta->Draw();

// 		DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
		DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
		DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
		
		canvasWeights->SaveAs(Form("%s/%s_%s_WeightsA.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		//	**********************************************************************************************************************
		//	******************************************* Assuming maximal correlation *********************************************
		//	**********************************************************************************************************************
		
		TGraph* graphWeightsB[11];
		for (Int_t i = 0; i< 11; i++){
			graphWeightsB[i] = NULL;
		}	
									
		// Declaration & calculation of combined spectrum							
		TString fileNameOutputWeightingB 							= Form("%s/%s_WeightingMethodBEta.dat",outputDir.Data(),centrality.Data());
		TGraphAsymmErrors* graphCombEtaInvYieldStat2760GeVB 	= NULL;
		TGraphAsymmErrors* graphCombEtaInvYieldSys2760GeVB 	= NULL;
		TGraphAsymmErrors* graphCombEtaInvYieldTot2760GeVB 	= CombinePtPointsSpectraFullCorrMat( 	statErrorCollection,	sysErrorCollection, 	
																											xPtLimitsEta, 13,
																											offSetsEta, offSetsEtaSys,
																											graphCombEtaInvYieldStat2760GeVB, graphCombEtaInvYieldSys2760GeVB,
																											fileNameOutputWeightingB,2
																										);
		graphCombEtaInvYieldStat2760GeVB->RemovePoint(0);
		graphCombEtaInvYieldTot2760GeVB->RemovePoint(0);
		graphCombEtaInvYieldSys2760GeVB->RemovePoint(0);
		graphCombEtaInvYieldStat2760GeVB->RemovePoint(0);
		graphCombEtaInvYieldTot2760GeVB->RemovePoint(0);
		graphCombEtaInvYieldSys2760GeVB->RemovePoint(0);
// 		cout << "\n\n\n\n Printing spectra method B" << endl;
// 		graphCombEtaInvYieldStat2760GeVB->Print();
//         graphCombEtaInvYieldSys2760GeVB->Print();
//         graphCombEtaInvYieldTot2760GeVB->Print();

		// Reading weights from output file for plotting
		ifstream fileWeightsReadB;
		fileWeightsReadB.open(fileNameOutputWeightingB,ios_base::in);
		cout << "reading" << fileNameOutputWeightingB << endl;
		Double_t xValuesReadB[50];
		Double_t weightsReadB[11][50];
		Int_t availableMeasB[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
		Int_t nMeasSetB 			= 2;
		Int_t nPtBinsReadB 		= 0;
		while(!fileWeightsReadB.eof() && nPtBinsReadB < 50){
			TString garbage = "";
			if (nPtBinsReadB == 0){
				fileWeightsReadB >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
				for (Int_t i = 0; i < nMeasSetB; i++){
					fileWeightsReadB >> availableMeasB[i] ;
				}	
				cout << "read following measurements: "; 
				for (Int_t i = 0; i < 11; i++){
					cout << availableMeasB[i] << "\t" ;
				}	
				cout << endl;
			} else {
				fileWeightsReadB >> xValuesReadB[nPtBinsReadB-1];
				for (Int_t i = 0; i < nMeasSetB; i++){
					fileWeightsReadB >> weightsReadB[availableMeasB[i]][nPtBinsReadB-1] ;
				}	
				cout << "read: "<<  nPtBinsReadB << "\t"<< xValuesReadB[nPtBinsReadB-1] << "\t" ;
				for (Int_t i = 0; i < nMeasSetB; i++){
					cout << weightsReadB[availableMeasB[i]][nPtBinsReadB-1] << "\t";
				}
				cout << endl;
			}
			nPtBinsReadB++;
		}
		nPtBinsReadB = nPtBinsReadB-2 ;
		fileWeightsReadB.close();

		for (Int_t i = 0; i < nMeasSetB; i++){
			graphWeightsB[availableMeasB[i]] = new TGraph(nPtBinsReadB,xValuesReadB,weightsReadB[availableMeasB[i]]);
			Int_t bin = 0;
			for (Int_t n = 0; n< nPtBinsReadB; n++){
				if (graphWeightsB[availableMeasB[i]]->GetY()[bin] == 0) graphWeightsB[availableMeasB[i]]->RemovePoint(bin);
				else bin++;
			}	
		}	

		//	**********************************************************************************************************************
		//	******************************************* Plotting weights method B ************************************************
		//	**********************************************************************************************************************
		
		canvasWeights->cd();
	
		histo2DWeights->Draw("copy");

		for (Int_t i = 0; i < nMeasSetB; i++){
			DrawGammaSetMarkerTGraph(graphWeightsB[availableMeasB[i]], markerStyleDet[availableMeasB[i]], markerSizeDet[availableMeasB[i]]*0.5, colorDet[availableMeasB[i]] , colorDet[availableMeasB[i]]);
			graphWeightsB[availableMeasB[i]]->Draw("p,same,e1");
		}	
		legendAccWeights->Draw();

		labelWeightsEnergy->Draw();
		labelWeightsEta->Draw();

// 		DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
		DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
		DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
		
		canvasWeights->SaveAs(Form("%s/%s_%s_WeightsB.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		

		//	**********************************************************************************************************************
		//	******************************************* Compare new spectrum to old average **************************************
		//	**********************************************************************************************************************
		
// 		TGraphAsymmErrors* graphCombEtaInvYieldStat2760GeVAForCompToOld	 = (TGraphAsymmErrors*)graphCombEtaInvYieldStat2760GeVA->Clone("graphCombEtaInvYieldStat2760GeVAForCompToOld");
// 		TGraphAsymmErrors* graphCombEtaInvYieldSys2760GeVAForCompToOld	 = (TGraphAsymmErrors*)graphCombEtaInvYieldSys2760GeVA->Clone("graphCombEtaInvYieldSys2760GeVAForCompToOld");
// 		
// 		TGraphAsymmErrors* graphCombEtaInvYieldStat2760GeVBForCompToOld	 = (TGraphAsymmErrors*)graphCombEtaInvYieldStat2760GeVB->Clone("graphCombEtaInvYieldStat2760GeVBForCompToOld");
// 		TGraphAsymmErrors* graphCombEtaInvYieldSys2760GeVBForCompToOld	 = (TGraphAsymmErrors*)graphCombEtaInvYieldSys2760GeVB->Clone("graphCombEtaInvYieldSys2760GeVBForCompToOld");
// 
// 		TGraphAsymmErrors* graphCombEtaInvYieldStat2760GeVOEMCBForCompToOld= (TGraphAsymmErrors*)graphCombEtaInvYieldStat2760GeVOEMC->Clone("graphCombEtaInvYieldStat2760GeVOEMCBForCompToOld");
// 		TGraphAsymmErrors* graphCombEtaInvYieldSys2760GeVOEMCForCompToOld = (TGraphAsymmErrors*)graphCombEtaInvYieldSys2760GeVOEMC->Clone("graphCombEtaInvYieldSys2760GeVOEMCForCompToOld");
// 		
// 		for (Int_t n = graphCombEtaInvYieldStat2760GeVAForCompToOld->GetN()-1; graphCombEtaInvYieldStat2760GeVAForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombEtaInvYieldStat2760GeVAForCompToOld->RemovePoint(n);
// 		}	
// 		for (Int_t n = graphCombEtaInvYieldSys2760GeVAForCompToOld->GetN()-1; graphCombEtaInvYieldSys2760GeVAForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombEtaInvYieldSys2760GeVAForCompToOld->RemovePoint(n);
// 		}	
// 		for (Int_t n = graphCombEtaInvYieldStat2760GeVBForCompToOld->GetN()-1; graphCombEtaInvYieldStat2760GeVBForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombEtaInvYieldStat2760GeVBForCompToOld->RemovePoint(n);
// 		}	
// 		for (Int_t n = graphCombEtaInvYieldSys2760GeVBForCompToOld->GetN()-1; graphCombEtaInvYieldSys2760GeVBForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombEtaInvYieldSys2760GeVBForCompToOld->RemovePoint(n);
// 		}	
// 		for (Int_t n = graphCombEtaInvYieldStat2760GeVOEMCBForCompToOld->GetN()-1; graphCombEtaInvYieldStat2760GeVOEMCBForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombEtaInvYieldStat2760GeVOEMCBForCompToOld->RemovePoint(n);
// 		}	
// 		for (Int_t n = graphCombEtaInvYieldSys2760GeVOEMCForCompToOld->GetN()-1; graphCombEtaInvYieldSys2760GeVOEMCForCompToOld->GetX()[n] > 6; n-- ){
// 			graphCombEtaInvYieldSys2760GeVOEMCForCompToOld->RemovePoint(n);
// 		}	
// 
// 		graphCombEtaInvYieldStat2760GeVAForCompToOld->Print();
// 		graphCombEtaInvYieldSys2760GeVAForCompToOld->Print();
// 		graphCombEtaInvYieldStat2760GeVBForCompToOld->Print();
// 		graphCombEtaInvYieldSys2760GeVBForCompToOld->Print();
// 		
// 		TGraphAsymmErrors*  graphRatioCombNewADivCombOldStat 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombEtaInvYieldStat2760GeVAForCompToOld, graphCombEtaInvYieldStat2760GeVOld, kTRUE);
// 		TGraphAsymmErrors*  graphRatioCombNewADivCombOldSys 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombEtaInvYieldSys2760GeVAForCompToOld, graphCombEtaInvYieldSys2760GeVOld, kTRUE);
// 		TGraphAsymmErrors*  graphRatioCombNewBDivCombOldStat 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombEtaInvYieldStat2760GeVBForCompToOld, graphCombEtaInvYieldStat2760GeVOld, kTRUE);
// 		TGraphAsymmErrors*  graphRatioCombNewBDivCombOldSys 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombEtaInvYieldSys2760GeVBForCompToOld, graphCombEtaInvYieldSys2760GeVOld, kTRUE);
// 		TGraphAsymmErrors*  graphRatioCombNewOEMCDivCombOldStat = CalculateGraphErrRatioToOtherTGraphErr(graphCombEtaInvYieldStat2760GeVOEMCBForCompToOld,
// 																										graphCombEtaInvYieldStat2760GeVOld, kTRUE);
// 		TGraphAsymmErrors*  graphRatioCombNewOECMDivCombOldSys 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombEtaInvYieldSys2760GeVOEMCForCompToOld, 
// 																										graphCombEtaInvYieldSys2760GeVOld, kTRUE);
		

		
		

		//	**********************************************************************************************************************
		//	******************************************* Ratio of Comb to Fit ****************************************
		//	**********************************************************************************************************************
		
// 		TCanvas* canvasRatioToOldCombined = new TCanvas("canvasRatioToOldCombined","",200,10,1350,900);  // gives the page size
// 		DrawGammaCanvasSettings( canvasRatioToOldCombined, 0.1, 0.02, 0.035, 0.09);
// 		canvasRatioToOldCombined->SetLogx();
// 	
// 		TH2F * histo2DRatioToOldCombined;
// 		histo2DRatioToOldCombined = new TH2F("histo2DRatioToOldCombined","histo2DRatioToOldCombined",11000,0.23,10.,1000,0.75,1.25);
// 		SetStyleHistoTH2ForGraphs(histo2DRatioToOldCombined, "#it{p}_{T} (GeV/#it{c})","#frac{Comb}{Comb Old}",0.035,0.04, 0.035,0.04, 1.,1.,510,505);
// 		histo2DRatioToOldCombined->GetXaxis()->SetMoreLogLabels();
// 		histo2DRatioToOldCombined->GetXaxis()->SetLabelOffset(-0.01);
// 	// 	histo2DRatioToOldCombined->GetYaxis()->SetRangeUser(-10,10);
// 		histo2DRatioToOldCombined->Draw("copy");
// 
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewBDivCombOldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
// 			graphRatioCombNewBDivCombOldSys->Draw("E2same");
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewBDivCombOldStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
// 			graphRatioCombNewBDivCombOldStat->Draw("p,same,e1");
// 
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewADivCombOldSys, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2, widthLinesBoxes, kTRUE);
// 			graphRatioCombNewADivCombOldSys->Draw("E2same");
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewADivCombOldStat, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
// 			graphRatioCombNewADivCombOldStat->Draw("p,same,e1");
// 
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewOECMDivCombOldSys, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2, widthLinesBoxes, kTRUE);
// 			graphRatioCombNewOECMDivCombOldSys->Draw("E2same");
// 			DrawGammaSetMarkerTGraphAsym(graphRatioCombNewOEMCDivCombOldStat, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
// 			graphRatioCombNewOEMCDivCombOldStat->Draw("p,same,e1");
// 			
// 			DrawGammaLines(0.23, 10. , 1., 1.,0.1, kGray+2);
// 			DrawGammaLines(0.23, 10. , 1.1, 1.1,0.1, kGray, 7);
// 			DrawGammaLines(0.23, 10. , 0.9, 0.9,0.1, kGray, 7);
// 			DrawGammaLines(0.23, 10. , 1.05, 1.05,0.1, kGray, 3);
// 			DrawGammaLines(0.23, 10. , 0.95, 0.95,0.1, kGray, 3);
// 			
// 			TLegend* legendRatioToOld = GetAndSetLegend2(0.67, 0.96-(0.035*3*1.35), 0.93, 0.96, 32);
// 			legendRatioToOld->AddEntry(graphRatioCombNewOECMDivCombOldSys,"PCM, EMCal");
// 			legendRatioToOld->AddEntry(graphRatioCombNewADivCombOldSys,"Method A");
// 			legendRatioToOld->AddEntry(graphRatioCombNewBDivCombOldSys,"Method B");
// 			legendRatioToOld->Draw();
// 
			TLatex *labelRatioToOldEnergy = new TLatex(0.15,0.89,collisionSystem2760GeV.Data());
			SetStyleTLatex( labelRatioToOldEnergy, 0.85*textSizeLabelsPixel,4);
			labelRatioToOldEnergy->SetTextFont(43);
// 			labelRatioToOldEnergy->Draw();
			TLatex *labelRatioToOldEta = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
			SetStyleTLatex( labelRatioToOldEta, 0.85*textSizeLabelsPixel,4);
			labelRatioToOldEta->SetTextFont(43);
// 			labelRatioToOldEta->Draw();
// 			
// 		canvasRatioToOldCombined->SaveAs(Form("%s/%s_%s_RatioOfCombToCombOld_PbPb2760GeV.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		//	**********************************************************************************************************************
		//	************************************** Ratio of Comb to Fit without Method A *****************************************
		//	**********************************************************************************************************************

// 		canvasRatioToOldCombined->cd();
// 		histo2DRatioToOldCombined->Draw("copy");
// 
// 			graphRatioCombNewBDivCombOldSys->Draw("E2same");
// 			graphRatioCombNewBDivCombOldStat->Draw("p,same,e1");
// 
// 			graphRatioCombNewOECMDivCombOldSys->Draw("E2same");
// 			graphRatioCombNewOEMCDivCombOldStat->Draw("p,same,e1");
// 			
// 			DrawGammaLines(0.23, 10. , 1., 1.,0.1, kGray+2);
// 			DrawGammaLines(0.23, 10. , 1.1, 1.1,0.1, kGray, 7);
// 			DrawGammaLines(0.23, 10. , 0.9, 0.9,0.1, kGray, 7);
// 			DrawGammaLines(0.23, 10. , 1.05, 1.05,0.1, kGray, 3);
// 			DrawGammaLines(0.23, 10. , 0.95, 0.95,0.1, kGray, 3);
// 			
// 			TLegend* legendRatioToOld2 = GetAndSetLegend2(0.67, 0.96-(0.035*3*1.35), 0.93, 0.96, 32);
// 			legendRatioToOld2->AddEntry(graphRatioCombNewOECMDivCombOldSys,"PCM, EMCal");
// 			legendRatioToOld2->AddEntry(graphRatioCombNewBDivCombOldSys,"All");
// 			legendRatioToOld2->Draw();
// 
// 			labelRatioToOldEnergy->Draw();
// 			labelRatioToOldEta->Draw();
// 			
// 		canvasRatioToOldCombined->SaveAs(Form("%s/%s_%s_RatioOfCombToCombOld_PbPb2760GeV_withoutMethodA.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));


		
		// 	*********************************************************************************************************************
		// 	************************************ Visualize relative errors ******************************************************
		// 	*********************************************************************************************************************
		
		TCanvas* canvasRelSysErr = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
		canvasRelSysErr->SetLogx();
	
		TH2F * histo2DRelSysErr;
		histo2DRelSysErr = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23,70.,1000,0,80.5);
		SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
		histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DRelSysErr->GetYaxis()->SetRangeUser(-10,10);
		histo2DRelSysErr->Draw("copy");	

			TLegend* legendRelSysErr = GetAndSetLegend2(0.62, 0.94-(0.035*nMeasSetA*1.35), 0.95, 0.94, 32);
			for (Int_t i = 0; i < nMeasSetA; i++){
				DrawGammaSetMarkerTGraph(sysErrorRelCollection[availableMeasA[i]], markerStyleDet[availableMeasA[i]], markerSizeDet[availableMeasA[i]]*0.5, colorDet[availableMeasA[i]],
										colorDet[availableMeasA[i]]);
				sysErrorRelCollection[availableMeasA[i]]->Draw("p,same,e1");
				legendRelSysErr->AddEntry(sysErrorRelCollection[availableMeasA[i]],nameMeasGlobal[availableMeasA[i]],"p");
			}	
			legendRelSysErr->Draw();

			TLatex *labelRelSysErrEnergy = new TLatex(0.15,0.89,collisionSystem2760GeV.Data());
			SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel,4);
			labelRelSysErrEnergy->SetTextFont(43);
			labelRelSysErrEnergy->Draw();
			TLatex *labelRelSysErrEta = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
			SetStyleTLatex( labelRelSysErrEta, 0.85*textSizeLabelsPixel,4);
			labelRelSysErrEta->SetTextFont(43);
			labelRelSysErrEta->Draw();
			
		canvasRelSysErr->SaveAs(Form("%s/%s_%s_RelSysErr.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		histo2DRelSysErr->GetYaxis()->SetRangeUser(0,30.5);
		histo2DRelSysErr->Draw("copy");	

			for (Int_t i = 0; i < nMeasSetA; i++){
				sysErrorRelCollection[availableMeasA[i]]->Draw("p,same,e1");
			}	
			legendRelSysErr->Draw();

			labelRelSysErrEnergy->Draw();
			labelRelSysErrEta->Draw();
			
		canvasRelSysErr->SaveAs(Form("%s/%s_%s_RelSysErrZoomed.%s",outputDir.Data(),meson.Data(), centrality.Data(),suffix.Data()));
		
		
		
		// 	*********************************************************************************************************************
		// 	************************************ Visualize relative errors ******************************************************
		// 	*********************************************************************************************************************
		
		TCanvas* canvasRelStatErr = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
		canvasRelStatErr->SetLogx();
	
		TH2F * histo2DRelStatErr;
		histo2DRelStatErr = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23,70.,1000,0,80.5);
		SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
		histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DRelStatErr->GetYaxis()->SetRangeUser(-10,10);
		histo2DRelStatErr->Draw("copy");	
// 			statErrorRelCollection[0]->SetBinContent(statErrorRelCollection[0]->GetNbinsX(), -1);
		
			TLegend* legendRelStatErr = GetAndSetLegend2(0.14, 0.94-(0.035*nMeasSetA*1.35), 0.45, 0.94, 32);
			for (Int_t i = 0; i < nMeasSetA; i++){
				DrawGammaSetMarker(statErrorRelCollection[availableMeasA[i]], markerStyleDet[availableMeasA[i]], markerSizeDet[availableMeasA[i]]*0.5, colorDet[availableMeasA[i]] , 
								colorDet[availableMeasA[i]]);
				statErrorRelCollection[availableMeasA[i]]->Draw("p,same,e1");
				legendRelStatErr->AddEntry(statErrorRelCollection[availableMeasA[i]],nameMeasGlobal[availableMeasA[i]],"p");
			}	
			legendRelStatErr->Draw();

			TLatex *labelRelStatErrEnergy = new TLatex(0.75,0.89,collisionSystem2760GeV.Data());
			SetStyleTLatex( labelRelStatErrEnergy, 0.85*textSizeLabelsPixel,4);
			labelRelStatErrEnergy->SetTextFont(43);
			labelRelStatErrEnergy->Draw();
			TLatex *labelRelStatErrEta = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
			SetStyleTLatex( labelRelStatErrEta, 0.85*textSizeLabelsPixel,4);
			labelRelStatErrEta->SetTextFont(43);
			labelRelStatErrEta->Draw();
			
		canvasRelStatErr->SaveAs(Form("%s/%s_%s_RelStatErr.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

// 		histo2DRelStatErr->GetYaxis()->SetRangeUser(0,30.5);
// 		histo2DRelStatErr->Draw("copy");	
// 			statErrorRelCollection[1]->SetBinContent(statErrorRelCollection[1]->GetNbinsX(), -1);
// 			statErrorRelCollection[2]->SetBinContent(statErrorRelCollection[2]->GetNbinsX(), -1);
// 			statErrorRelCollection[2]->SetBinContent(statErrorRelCollection[2]->GetNbinsX()-1, -1);
// 			for (Int_t i = 0; i < nMeasSetA; i++){
// 				statErrorRelCollection[availableMeasA[i]]->Draw("p,same,e1");
// 			}	
// 			legendRelStatErr->Draw();
// 
// 			labelRelStatErrEnergy->Draw();
// 			labelRelStatErrEta->Draw();
// 			
// 		canvasRelStatErr->SaveAs(Form("%s/%s_RelStatErrZoomed.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		
		
		// 	*********************************************************************************************************************
		// 	************************ Visualize relative total errors of different combination methods ***************************
		// 	*********************************************************************************************************************
		TGraphAsymmErrors* graphCombEtaInvYieldRelTot2760GeVB 	= CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldTot2760GeVB, "relativeTotalError_MethodB");
		TGraphAsymmErrors* graphCombEtaInvYieldRelStat2760GeVB 	= CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldStat2760GeVB, "relativeStatError_MethodB");
		TGraphAsymmErrors* graphCombEtaInvYieldRelSys2760GeVB 	= CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldSys2760GeVB, "relativeSysError_MethodB");
		
		TGraphAsymmErrors* graphCombEtaInvYieldRelStat2760GeVA 	= CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldStat2760GeVA, "relativeStatError_MethodA");
		TGraphAsymmErrors* graphCombEtaInvYieldRelSys2760GeVA 	= CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldSys2760GeVA, "relativeSysError_MethodA");
		TGraphAsymmErrors* graphCombEtaInvYieldRelTot2760GeVA 	= CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldTot2760GeVA, "relativeTotalError_MethodA");
		
		TGraphAsymmErrors* graphCombEtaInvYieldRelTot2760GeVOEMC = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldTot2760GeVOEMC, "relativeTotalError_OEMC");
// 		TGraphAsymmErrors* graphCombEtaInvYieldRelTot2760GeVOld 	= CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldTot2760GeVOld, "relativeTotalError_Old");

		
		TCanvas* canvasRelTotErr = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
		canvasRelTotErr->SetLogx();
	
		TH2F * histo2DRelTotErr;
		histo2DRelTotErr = new TH2F("histo2DRelTotErr","histo2DRelTotErr",11000,0.23,70.,1000,0,80.5);
		SetStyleHistoTH2ForGraphs(histo2DRelTotErr, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DRelTotErr->GetXaxis()->SetMoreLogLabels();
		histo2DRelTotErr->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DRelTotErr->GetYaxis()->SetRangeUser(-10,10);
		histo2DRelTotErr->Draw("copy");	

// 			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot2760GeVOld, markerStyleComb+1, markerSizeComb, kBlack , kBlack);
// 			graphCombEtaInvYieldRelTot2760GeVOld->Draw("p,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot2760GeVB, markerStyleComb, markerSizeComb, colorComb , colorComb);
			graphCombEtaInvYieldRelTot2760GeVB->Draw("p,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot2760GeVA, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
			graphCombEtaInvYieldRelTot2760GeVA->Draw("p,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot2760GeVOEMC, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
			graphCombEtaInvYieldRelTot2760GeVOEMC->Draw("p,same,e1");

			TLegend* legendRelTotErr = GetAndSetLegend2(0.14, 0.94-(0.035*4*1.35), 0.45, 0.94, 32);
// 			legendRelTotErr->AddEntry(graphCombEtaInvYieldRelTot2760GeVOld,"PCMold, EMCalold","p");
			legendRelTotErr->AddEntry(graphCombEtaInvYieldRelTot2760GeVOEMC,"PCM, EMCal","p");
			legendRelTotErr->AddEntry(graphCombEtaInvYieldRelTot2760GeVA,"Method A","p");
			legendRelTotErr->AddEntry(graphCombEtaInvYieldRelTot2760GeVB,"Method B","p");
			legendRelTotErr->Draw();

			TLatex *labelRelTotErrEnergy = new TLatex(0.75,0.89,collisionSystem2760GeV.Data());
			SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel,4);
			labelRelTotErrEnergy->SetTextFont(43);
			labelRelTotErrEnergy->Draw();
			TLatex *labelRelTotErrEta = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
			SetStyleTLatex( labelRelTotErrEta, 0.85*textSizeLabelsPixel,4);
			labelRelTotErrEta->SetTextFont(43);
			labelRelTotErrEta->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelTotErr.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		histo2DRelTotErr->Draw("copy");	
// 			graphCombEtaInvYieldRelTot2760GeVOld->Draw("p,same,e1");
			graphCombEtaInvYieldRelTot2760GeVB->Draw("p,same,e1");
			graphCombEtaInvYieldRelTot2760GeVA->Draw("p,same,e1");
			graphCombEtaInvYieldRelTot2760GeVOEMC->Draw("p,same,e1");

			TLegend* legendRelTotErrWOA = GetAndSetLegend2(0.14, 0.94-(0.035*4*1.35), 0.45, 0.94, 32);
// 			legendRelTotErrWOA->AddEntry(graphCombEtaInvYieldRelTot2760GeVOld,"PCMold, EMCalold","p");
			legendRelTotErrWOA->AddEntry(graphCombEtaInvYieldRelTot2760GeVOEMC,"PCM, EMCal","p");
			legendRelTotErrWOA->AddEntry(graphCombEtaInvYieldRelTot2760GeVB,"All","p");
			legendRelTotErrWOA->Draw();

			labelRelTotErrEnergy->Draw();
			labelRelTotErrEta->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelTotErr_withoutMethodA.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		
		histo2DRelTotErr->GetYaxis()->SetRangeUser(0,30.5);
		histo2DRelTotErr->Draw("copy");	
// 			graphCombEtaInvYieldRelTot2760GeVOld->Draw("p,same,e1");
			graphCombEtaInvYieldRelTot2760GeVB->Draw("p,same,e1");
			graphCombEtaInvYieldRelTot2760GeVA->Draw("p,same,e1");
			graphCombEtaInvYieldRelTot2760GeVOEMC->Draw("p,same,e1");

			legendRelTotErr->Draw();

			labelRelTotErrEnergy->Draw();
			labelRelTotErrEta->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelTotErrZoomed.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		histo2DRelTotErr->Draw("copy");	
// 			graphCombEtaInvYieldRelTot2760GeVOld->Draw("p,same,e1");
			graphCombEtaInvYieldRelTot2760GeVB->Draw("p,same,e1");
			graphCombEtaInvYieldRelTot2760GeVOEMC->Draw("p,same,e1");

			legendRelTotErrWOA->Draw();

			labelRelTotErrEnergy->Draw();
			labelRelTotErrEta->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelTotErrZoomed_withoutMethodA.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		histo2DRelTotErr->GetYaxis()->SetRangeUser(0,80.5);
		histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
		histo2DRelTotErr->Draw("copy");	

			graphCombEtaInvYieldRelTot2760GeVB->Draw("p,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelStat2760GeVB, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
			graphCombEtaInvYieldRelStat2760GeVB->Draw("l,x0,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelSys2760GeVB, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
			graphCombEtaInvYieldRelSys2760GeVB->SetLineStyle(7);
			graphCombEtaInvYieldRelSys2760GeVB->Draw("l,x0,same,e1");

			TLegend* legendRelTotErr2 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
			legendRelTotErr2->AddEntry(graphCombEtaInvYieldRelTot2760GeVB,"tot","p");
			legendRelTotErr2->AddEntry(graphCombEtaInvYieldRelStat2760GeVB,"stat","l");
			legendRelTotErr2->AddEntry(graphCombEtaInvYieldRelSys2760GeVB,"sys","l");
			legendRelTotErr2->Draw();

			labelRelTotErrEnergy->Draw();
			labelRelTotErrEta->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelMethodBdecomp.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		histo2DRelTotErr->GetYaxis()->SetRangeUser(0,80.5);
		histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
			histo2DRelTotErr->Draw("copy");	

			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot2760GeVA, markerStyleComb, markerSizeComb, colorComb , colorComb);		
			graphCombEtaInvYieldRelTot2760GeVA->Draw("p,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelStat2760GeVA, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
			graphCombEtaInvYieldRelStat2760GeVA->Draw("l,x0,same,e1");
			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelSys2760GeVA, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
			graphCombEtaInvYieldRelSys2760GeVA->SetLineStyle(7);
			graphCombEtaInvYieldRelSys2760GeVA->Draw("l,x0,same,e1");

			TLegend* legendRelTotErr3 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
			legendRelTotErr3->AddEntry(graphCombEtaInvYieldRelTot2760GeVA,"tot","p");
			legendRelTotErr3->AddEntry(graphCombEtaInvYieldRelStat2760GeVA,"stat","l");
			legendRelTotErr3->AddEntry(graphCombEtaInvYieldRelSys2760GeVA,"sys","l");
			legendRelTotErr3->Draw();

			labelRelTotErrEnergy->Draw();
			labelRelTotErrEta->Draw();
			
		canvasRelTotErr->SaveAs(Form("%s/%s_%s_RelMethodAdecomp.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		
		
		//	**********************************************************************************************************************
		//	************************************* Calculating bin shifted spectra & fitting **************************************
		//	**********************************************************************************************************************
		
		// Cloning spectra
		TGraphAsymmErrors* graphCombEtaInvYieldTot2760GeVAUnShifted 		= (TGraphAsymmErrors*)graphCombEtaInvYieldTot2760GeVA->Clone("Unshifted"); 
		TGraphAsymmErrors* graphCombEtaInvYieldStat2760GeVAUnShifted 	= (TGraphAsymmErrors*)graphCombEtaInvYieldStat2760GeVA->Clone("UnshiftedStat"); 
		TGraphAsymmErrors* graphCombEtaInvYieldSys2760GeVAUnShifted 		= (TGraphAsymmErrors*)graphCombEtaInvYieldSys2760GeVA->Clone("UnshiftedSys"); 

		TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeVUnShifted 		= (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV->Clone("UnshiftedStatPCM"); 
		TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeVUnShifted 		= (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV->Clone("UnshiftedSysPCM"); 

		TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeVUnshifted 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV->Clone("UnshiftedStatEMCal"); 
		TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeVUnshifted 		= (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV->Clone("UnshiftedSysEMCal"); 
		
		// Calculating binshifts
		TF1* fitBylinkinEtaPbPb2760GeVPt 				= FitObject("tcm","BylinkinFitEta","Eta",graphCombEtaInvYieldStat2760GeVA,1.,30.);
		graphCombEtaInvYieldStat2760GeVA->Fit(fitBylinkinEtaPbPb2760GeVPt,"QNRMEX0+","",1.0,30.);
// 		graphPCMEtaInvYieldStatPbPb2760GeV->Fit(fitBylinkinEtaPbPb2760GeVPt,"QNRMEX0+","",1.0,30.);
		

		if(bWCorrection.CompareTo("X")==0 ){
			TF1* fitBylinkinEtaPbPb2760GeVPtMult 				= FitObject("tcm","BylinkinFitEta","Eta");		
			graphCombEtaInvYieldTot2760GeVA			= ApplyXshift(graphCombEtaInvYieldTot2760GeVA, fitBylinkinEtaPbPb2760GeVPtMult,"Eta");
			
			graphCombEtaInvYieldStat2760GeVA 		= ApplyXshiftIndividualSpectra (graphCombEtaInvYieldTot2760GeVA, 
																							graphCombEtaInvYieldStat2760GeVA, 
																							fitBylinkinEtaPbPb2760GeVPtMult,
																							2, 3,"Eta");
			
			graphCombEtaInvYieldSys2760GeVA 			= ApplyXshiftIndividualSpectra (graphCombEtaInvYieldTot2760GeVA, 
																							graphCombEtaInvYieldSys2760GeVA, 
																							fitBylinkinEtaPbPb2760GeVPtMult, 
																							2, 3,"Eta");
			
			graphPCMEtaInvYieldStatPbPb2760GeV			= ApplyXshiftIndividualSpectra( graphCombEtaInvYieldTot2760GeVA,
																							graphPCMEtaInvYieldStatPbPb2760GeV,
																							fitBylinkinEtaPbPb2760GeVPtMult, 
																							2, 3,"Eta");
			
			graphPCMEtaInvYieldSysPbPb2760GeV			= ApplyXshiftIndividualSpectra( graphCombEtaInvYieldTot2760GeVA, 
																							graphPCMEtaInvYieldSysPbPb2760GeV, 
																							fitBylinkinEtaPbPb2760GeVPtMult, 
																							2, 3,"Eta");
			
			graphEMCalEtaInvYieldStatPbPb2760GeV 		= ApplyXshiftIndividualSpectra( graphCombEtaInvYieldTot2760GeVA, 
																							graphEMCalEtaInvYieldStatPbPb2760GeV, 
																							fitBylinkinEtaPbPb2760GeVPtMult,
																							9, 3,"Eta");	
			
			graphEMCalEtaInvYieldSysPbPb2760GeV 			= ApplyXshiftIndividualSpectra( graphCombEtaInvYieldTot2760GeVA, 
																							graphEMCalEtaInvYieldSysPbPb2760GeV, 
																							fitBylinkinEtaPbPb2760GeVPtMult,
																							9, 3,"Eta");		

					
			TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
			DrawGammaCanvasSettings( canvasDummy2,  0.1, 0.01, 0.015, 0.08);
			canvasDummy2->SetLogy();
			canvasDummy2->SetLogx();
			TH2F * histo2DDummy2;
			histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.8,30.,1000,1e-8,1e1);
			SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
			histo2DDummy2->DrawCopy(); 
		
			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat2760GeVAUnShifted, 20,3, kRed, kRed, widthLinesBoxes, kTRUE);
			graphCombEtaInvYieldStat2760GeVAUnShifted->Draw("pEsame");
			DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldTot2760GeVA, 24, 3, kBlack, kBlack, widthLinesBoxes, kTRUE);
			graphCombEtaInvYieldTot2760GeVA->Draw("pEsame");
			DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldSysPbPb2760GeV, markerStyleDet[0] ,/*markerSizeDet[0]*/2, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
			graphPCMEtaInvYieldSysPbPb2760GeV->Draw("pEsame");
			DrawGammaSetMarkerTGraphAsym(graphEMCalEtaInvYieldSysPbPb2760GeV, markerStyleDet[2] ,/*markerSizeDet[2]*/2.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
			graphEMCalEtaInvYieldSysPbPb2760GeV->Draw("pEsame");

			TLegend* legendXdummy = new TLegend(0.6,0.8,0.9,0.95);
			legendXdummy->SetFillColor(0);
			legendXdummy->SetLineColor(0);
			legendXdummy->SetTextFont(42);
			legendXdummy->SetTextSize(0.035);
			legendXdummy->AddEntry(graphCombEtaInvYieldStat2760GeVAUnShifted,"combined unshifted","lp");
// 			legendXdummy->AddEntry(graphCombEtaInvYieldTot2760GeVA,"combined unshifted (stat. err.)","lp");
			legendXdummy->AddEntry(graphPCMEtaInvYieldSysPbPb2760GeV,"PCM","fp");
			legendXdummy->AddEntry(graphEMCalEtaInvYieldSysPbPb2760GeV,"EMCal","fp");
			legendXdummy->Draw();

			canvasDummy2->Update();
			canvasDummy2->Print(Form("%s/%s_%s_ComparisonShiftedEta_PbPb2760GeV.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		}
			
		
		Double_t paramTCM[5] = {graphCombEtaInvYieldTot2760GeVA->GetY()[0],0.3,graphCombEtaInvYieldTot2760GeVA->GetY()[0]/10000,0.8,3};
		TF1* fitTCMInvYieldEta2760GeV = FitObject("tcm","fitTCMInvYieldEta2760GeV","Eta",graphCombEtaInvYieldTot2760GeVA,0.8,30.,paramTCM,"QNRMEX0+");
		fitTCMInvYieldEta2760GeV = FitObject("tcm","fitTCMInvYieldEta2760GeV","Eta",graphCombEtaInvYieldTot2760GeVA,0.8,30. ,paramTCM,"QNRMEX0+");

		cout << WriteParameterToFile(fitTCMInvYieldEta2760GeV)<< endl;
		
		TGraphAsymmErrors* graphRatioCombCombFitTot2760GeVA 	= (TGraphAsymmErrors*)graphCombEtaInvYieldTot2760GeVA->Clone();
		graphRatioCombCombFitTot2760GeVA 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitTot2760GeVA, fitTCMInvYieldEta2760GeV); 
		TGraphAsymmErrors* graphRatioCombCombFitStat2760GeVA 	= (TGraphAsymmErrors*)graphCombEtaInvYieldStat2760GeVA->Clone();
		graphRatioCombCombFitStat2760GeVA 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitStat2760GeVA, fitTCMInvYieldEta2760GeV); 
		TGraphAsymmErrors* graphRatioCombCombFitSys2760GeVA 	= (TGraphAsymmErrors*)graphCombEtaInvYieldSys2760GeVA->Clone();
		graphRatioCombCombFitSys2760GeVA 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitSys2760GeVA, fitTCMInvYieldEta2760GeV); 

		TGraphAsymmErrors* graphRatioPCMCombFitStat2760GeV 		= (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV->Clone();
		graphRatioPCMCombFitStat2760GeV 						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat2760GeV, fitTCMInvYieldEta2760GeV); 
		TGraphAsymmErrors* graphRatioPCMCombFitSys2760GeV 		= (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV->Clone();
		graphRatioPCMCombFitSys2760GeV 							= CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys2760GeV, fitTCMInvYieldEta2760GeV); 
		TGraphAsymmErrors* graphRatioEMCalCombFitStat2760GeV 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV->Clone();
		graphRatioEMCalCombFitStat2760GeV 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitStat2760GeV, fitTCMInvYieldEta2760GeV); 
		TGraphAsymmErrors* graphRatioEMCalCombFitSys2760GeV 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV->Clone();
		graphRatioEMCalCombFitSys2760GeV 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitSys2760GeV, fitTCMInvYieldEta2760GeV); 
		

		//	**********************************************************************************************************************
		//	******************************************* Ratio of Comb to Fit ****************************************
		//	**********************************************************************************************************************
		textSizeLabelsPixel = 48;
		TCanvas* canvasRatioToCombFit = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRatioToCombFit, 0.12, 0.01, 0.01, 0.11);
		canvasRatioToCombFit->SetLogx();

			Double_t textsizeLabelsPP = 0;
			Double_t textsizeFacPP= 0;
			if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
				textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
				textsizeFacPP = (Double_t)1./canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
			} else {
				textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
				textsizeFacPP = (Double_t)1./canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
			}
			cout << textsizeLabelsPP << endl;
		
		TH2F * histo2DRatioToCombFit;
		histo2DRatioToCombFit = new TH2F("histo2DRatioToCombFit","histo2DRatioToCombFit",1000,0.23,70.,1000,0.2,4.	);
		SetStyleHistoTH2ForGraphs(histo2DRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
								0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
		histo2DRatioToCombFit->GetXaxis()->SetMoreLogLabels();
		histo2DRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
		histo2DRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
		histo2DRatioToCombFit->Draw("copy");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys2760GeVA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
		graphRatioCombCombFitSys2760GeVA->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat2760GeVA, markerStyleComb, markerSizeComb, colorComb , colorComb);
		graphRatioCombCombFitStat2760GeVA->Draw("p,same,e1");

		DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
		DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);

		TLatex *labelRatioToFitEnergy = new TLatex(0.67,0.92,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelRatioToFitEnergy, 0.85*textSizeLabelsPixel,4);
		labelRatioToFitEnergy->SetTextFont(43);
		labelRatioToFitEnergy->Draw();
		TLatex *labelRatioToFitEta = new TLatex(0.73,0.87,"#eta #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRatioToFitEta, 0.85*textSizeLabelsPixel,4);
		labelRatioToFitEta->SetTextFont(43);
		labelRatioToFitEta->Draw();

			
		canvasRatioToCombFit->SaveAs(Form("%s/%s_%s_RatioOfCombToCombFit_PbPb2760GeV.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		//	**********************************************************************************************************************
		//	******************************************* Ratio of Individual meas to Fit ******************************************
		//	**********************************************************************************************************************
		
		canvasRatioToCombFit->cd();
		histo2DRatioToCombFit->Draw("copy");
	
		DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitSys2760GeV, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitStat2760GeV, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
		DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitSys2760GeV, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitStat2760GeV, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
		
		graphRatioPCMCombFitSys2760GeV->Draw("E2same");
		graphRatioEMCalCombFitSys2760GeV->Draw("E2same");
		
		graphRatioPCMCombFitStat2760GeV->Draw("p,same,e");
		graphRatioEMCalCombFitStat2760GeV->Draw("p,same,e");

		DrawGammaLines(0.23, 70. , 1., 1.,0.5, kGray+2);
		DrawGammaLines(0.23, 70. , 1.1, 1.1,0.5, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.9, 0.9,0.5, kGray, 7);
		
		labelRatioToFitEnergy->Draw();
		labelRatioToFitEta->Draw();
	
		// ****************************** Definition of the Legend ******************************************
		// **************** Row def ************************
		Double_t rowsLegendOnlyEtaRatio[5] 		= {0.92,0.88,0.84,0.80,0.76};
		Double_t rowsLegendOnlyEtaRatioAbs[5] 	= {0.91,2.2,2.1,2.0,1.9};
		Double_t columnsLegendOnlyEtaRatio[3] 	= {0.15,0.32, 0.38};
		Double_t columnsLegendOnlyEtaRatioAbs[3]= {0.15,1.04, 1.37};
		Double_t lengthBox						= 0.2/2;
		Double_t heightBox						= 0.08/2;
		// ****************** first Column **************************************************
		TLatex *textPCMOnlyRatioEta = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[1],"PCM");
		SetStyleTLatex( textPCMOnlyRatioEta, 0.85*textSizeLabelsPixel,4);
		textPCMOnlyRatioEta->SetTextFont(43);
		textPCMOnlyRatioEta->Draw();
		TLatex *textEMCalOnlyRatioEta = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[3],"EMCal");
		SetStyleTLatex( textEMCalOnlyRatioEta,  0.85*textSizeLabelsPixel,4);
		textEMCalOnlyRatioEta->SetTextFont(43);
		textEMCalOnlyRatioEta->Draw();
		
		// ****************** second Column *************************************************
		TLatex *textStatOnlyRatioEta = new TLatex(columnsLegendOnlyEtaRatio[1],rowsLegendOnlyEtaRatio[0] ,"stat");
		SetStyleTLatex( textStatOnlyRatioEta, 0.85*textSizeLabelsPixel,4);
		textStatOnlyRatioEta->SetTextFont(43);
		textStatOnlyRatioEta->Draw();
		TLatex *textSysOnlyRatioEta = new TLatex(columnsLegendOnlyEtaRatio[2] ,rowsLegendOnlyEtaRatio[0],"syst");
		SetStyleTLatex( textSysOnlyRatioEta, 0.85*textSizeLabelsPixel,4);
		textSysOnlyRatioEta->SetTextFont(43);
		textSysOnlyRatioEta->Draw();
		TMarker* markerPCMEtaOnlyRatioEta = CreateMarkerFromGraph(graphRatioPCMCombFitSys2760GeV,columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[1],1);
		markerPCMEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[1]);
		TMarker* markerEMCalEtaOnlyRatioEta = CreateMarkerFromGraph(graphRatioEMCalCombFitSys2760GeV, columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[3],1);
		markerEMCalEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[3]);

		TBox* boxPCMEtaOnlyRatioEta = CreateBoxFromGraph(graphRatioPCMCombFitSys2760GeV, columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[1]- heightBox,
														 columnsLegendOnlyEtaRatioAbs[2]+ 3*lengthBox, rowsLegendOnlyEtaRatioAbs[1]+ heightBox);
		boxPCMEtaOnlyRatioEta->Draw("l");
		TBox* boxEMCalEtaOnlyRatioEta = CreateBoxFromGraph(graphRatioEMCalCombFitSys2760GeV, columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[3]- heightBox,
														   columnsLegendOnlyEtaRatioAbs[2]+ 3*lengthBox, rowsLegendOnlyEtaRatioAbs[3]+ heightBox);
		boxEMCalEtaOnlyRatioEta->Draw("l");
		
		canvasRatioToCombFit->SaveAs(Form("%s/%s_%s_RatioOfIndividualMeasToCombFit_PbPb2760GeV.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

	
	
	
	

		//	**********************************************************************************************************************
		//	******************************** Cross section for eta single measurement 2.76TeV ************************************
		//	**********************************************************************************************************************
		
		TCanvas* canvasXSectionEta = new TCanvas("canvasXSectionEta","",200,10,1350,1350*1.15);  // gives the page size
		DrawGammaCanvasSettings( canvasXSectionEta, 0.14, 0.02, 0.02, 0.09);
		canvasXSectionEta->SetLogx();
		canvasXSectionEta->SetLogy();
		
		TH2F * histo2DXSectionEta;
		histo2DXSectionEta = new TH2F("histo2DXSectionEta","histo2DXSectionEta",11000,0.23,70.,1000,2e-8,1e3);
		SetStyleHistoTH2ForGraphs(histo2DXSectionEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 1.,1.45);
		histo2DXSectionEta->GetXaxis()->SetMoreLogLabels();
		histo2DXSectionEta->GetXaxis()->SetLabelOffset(-0.01);
		histo2DXSectionEta->Draw("copy");
		histo2DXSectionEta->GetXaxis()->SetRangeUser(0.8,30.);

		DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldSysPbPb2760GeV, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
		graphPCMEtaInvYieldSysPbPb2760GeV->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphEMCalEtaInvYieldSysPbPb2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
		graphEMCalEtaInvYieldSysPbPb2760GeV->Draw("E2same");
	
		DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldStatPbPb2760GeV,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
		graphPCMEtaInvYieldStatPbPb2760GeV->Draw("p,same,e1");
			
		DrawGammaSetMarkerTGraphAsym(graphEMCalEtaInvYieldStatPbPb2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
// 		graphEMCalEtaInvYieldStatPbPb2760GeV->Draw("p,same,e1");
		
// 		fitBylinkinEtaPbPb2760GeVPt->SetRange(1.0,30.);
		fitBylinkinEtaPbPb2760GeVPt->Draw("same");

		TLatex *labelEnergyXSectionEta = new TLatex(0.5,0.92,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelEnergyXSectionEta, 0.035,4);
		labelEnergyXSectionEta->Draw();
		TLatex *labelDetSysXSectionEta = new TLatex(0.5,0.88,"#eta #rightarrow #gamma#gamma");
		SetStyleTLatex( labelDetSysXSectionEta, 0.035,4);
		labelDetSysXSectionEta->Draw();

		TLegend* legendXSectionEta = new TLegend(0.62,0.66,0.9,0.86);
		legendXSectionEta->SetFillColor(0);
		legendXSectionEta->SetLineColor(0);
		legendXSectionEta->SetTextFont(42);
		legendXSectionEta->SetTextSize(0.035);
		legendXSectionEta->AddEntry(graphPCMEtaInvYieldSysPbPb2760GeV,"PCM","fp");
		legendXSectionEta->AddEntry(graphEMCalEtaInvYieldStatPbPb2760GeV,"EMCal","fp");
// 		legendXSectionEta->AddEntry(fitBylinkinEtaPbPb2760GeVPt, "Bylinkin Fit","l");
		legendXSectionEta->Draw();
   
		canvasXSectionEta->SaveAs(Form("%s/%s_%s_YieldsEtaCompAllSystems.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));

		canvasXSectionEta->cd();

		histo2DXSectionEta->Draw("copy");

		graphPCMEtaInvYieldSysPbPb2760GeV->Draw("E2same");
		graphEMCalEtaInvYieldSysPbPb2760GeV->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys2760GeVA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
		graphCombEtaInvYieldSys2760GeVA->Draw("E2same");
	
		graphPCMEtaInvYieldStatPbPb2760GeV->Draw("p,same,e1");
		graphEMCalEtaInvYieldStatPbPb2760GeV->Draw("p,same,e1");
		DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat2760GeVA, markerStyleComb, markerSizeComb, colorComb , colorComb);
		graphCombEtaInvYieldStat2760GeVA->Draw("p,same,e1");
		
		fitBylinkinEtaPbPb2760GeVPt->Draw("same");

		labelEnergyXSectionEta->Draw();
		labelDetSysXSectionEta->Draw();

		legendXSectionEta->AddEntry(graphCombEtaInvYieldSys2760GeVA,"comb","fp");
		legendXSectionEta->AddEntry(fitBylinkinEtaPbPb2760GeVPt, "Bylinkin Fit","l");
		legendXSectionEta->Draw();

// 		DrawGammaSetMarkerTGraphAsym(graphChargedHadronsStatPP2760GeV, markerStyleDet[1], markerSizeDet[1]*0.75, kBlack , kBlack, widthLinesBoxes);
// 		graphChargedHadronsStatPP2760GeV->Draw("E2same");
		
		canvasXSectionEta->SaveAs(Form("%s/%s_%s_YieldEtaCompAllSystems_Comb.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
  		
		
		//	**********************************************************************************************************************
		//	*************************	 Combination of RAA 	*******************************************************************

// 		Int_t bin0PCM =0;
// 		Int_t bin0EMCal =0;
// 		
// 		TGraphAsymmErrors* graphCombEtaRAAStat2760GeVOEMC= NULL;
// 		TGraphAsymmErrors* graphCombEtaRAASys2760GeVOEMC = NULL;
// 		TGraphAsymmErrors* graphCombEtaRAATot2760GeVOEMC = CombinePtPointsRAA( 	graphPCMEtaRAAStat2760GeV,	graphPCMEtaRAASys2760GeV, 	
// 																				graphEMCalEtaRAAStat2760GeV, graphEMCalEtaRAASys2760GeV,																graphCombEtaRAAStat2760GeVOEMC, graphCombEtaRAASys2760GeVOEMC,
// 																				xPtLimitsEta, 23, offSetsEta, bin0PCM, bin0EMCal, 
// 																				kFALSE,1);
// 		
// // TGraphAsymmErrors* CombinePtPointsRAA( TGraphAsymmErrors* graphStatErrPCM,		TGraphAsymmErrors* graphSystPCM,
// // 										  TGraphAsymmErrors* graphStatErrPHOS,		TGraphAsymmErrors* graphSystPHOS,
// //  									  TGraphAsymmErrors* &graphStatComb, 	TGraphAsymmErrors* &graphSystComb,  
// // 													Double_t* xPtLimits,	Int_t nPtLimits,
// // 													Int_t offset,			Int_t bin0PCM,	 Int_t bin0PHOS, Bool_t kRemoveLastPCMPoint=kFALSE, Int_t nBinsPCMRem = 1){
// 		
// // 		cout << "Printing of PCM" << endl;
// // 		graphPCMEtaRAAStat2760GeV->Print();
// // 		graphPCMEtaRAASys2760GeV->Print();
// // 		
// // 		cout << "Printing of EMCal" << endl;
// // 		graphEMCalEtaRAAStat2760GeV->Print();
// // 		graphEMCalEtaRAASys2760GeV->Print();
// 
// 		cout << "Printing of combined RAA" << endl;
// // 		graphCombEtaRAAStat2760GeVOEMC->RemovePoint(0);
// //         graphCombEtaRAASys2760GeVOEMC->RemovePoint(0);
// 		graphCombEtaRAAStat2760GeVOEMC->Print();
// 		graphCombEtaRAASys2760GeVOEMC->Print();
// 

		TCanvas* canvasRAAcombo = new TCanvas("canvasRAAcombo","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRAAcombo, 0.08, 0.02, 0.035, 0.09);
// 		canvasRAAcombo->SetLogx();
	
		TH2F * histo2DRAAcombo;
		histo2DRAAcombo = new TH2F("histo2DRAAcombo","histo2DRAAcombo",11000,0.23,70.,1000,-0.5,1.1);
		SetStyleHistoTH2ForGraphs(histo2DRAAcombo, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,1.);
// 		histo2DRAAcombo->GetXaxis()->SetMoreLogLabels();	
// 		histo2DRAAcombo->GetXaxis()->SetLabelOffset(-0.01);
		histo2DRAAcombo->GetYaxis()->SetRangeUser(0.,1.);
		histo2DRAAcombo->GetXaxis()->SetRangeUser(0.,30.);
		histo2DRAAcombo->Draw("copy");
		
		DrawGammaSetMarkerTGraphAsym(graphPCMEtaRAASys2760GeV, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
		graphPCMEtaRAASys2760GeV->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPCMEtaRAAStat2760GeV,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
		graphPCMEtaRAAStat2760GeV->Draw("p,same,e1");		
		
// 		DrawGammaSetMarkerTGraphAsym(graphEMCalEtaRAASys2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
// 		graphEMCalEtaRAASys2760GeV->Draw("E2same");			
// 		DrawGammaSetMarkerTGraphAsym(graphEMCalEtaRAAStat2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
// 		graphEMCalEtaRAAStat2760GeV->Draw("p,same,e1");
		
// 		DrawGammaSetMarkerTGraphAsym(graphCombEtaRAASys2760GeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
// 		graphCombEtaRAASys2760GeVOEMC->Draw("E2same");
// 		DrawGammaSetMarkerTGraphAsym(graphCombEtaRAAStat2760GeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb);
// 		graphCombEtaRAAStat2760GeVOEMC->Draw("p,same,e1");

		TLegend* legendRAAcombo = new TLegend(0.62,0.66,0.9,0.86);
		legendRAAcombo->SetFillColor(0);
		legendRAAcombo->SetLineColor(0);
		legendRAAcombo->SetTextFont(42);
		legendRAAcombo->SetTextSize(0.035);
// 		legendRAAcombo->AddEntry(graphCombEtaRAAStat2760GeVOEMC,"combined","p");		
		legendRAAcombo->AddEntry(graphPCMEtaRAASys2760GeV,"PCM","p");
// 		legendRAAcombo->AddEntry(graphEMCalEtaRAASys2760GeV,"EMCal","p");
		legendRAAcombo->Draw();

		TLatex *labelRAAEnergy = new TLatex(0.7,0.20,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelRAAEnergy, 0.85*textSizeLabelsPixel,4);
		labelRAAEnergy->SetTextFont(43);
// 		labelRAAEnergy->Draw();
		TLatex *labelRAAEta = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRAAEta, 0.85*textSizeLabelsPixel,4);
		labelRAAEta->SetTextFont(43);
// 		labelRAAEta->Draw();

		canvasRAAcombo->SaveAs(Form("%s/%s_%s_RAAEta_combined.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
		
		//	**********************************************************************************************************************
		//	*************************	 Combination of RCP 	*******************************************************************
		cout << "\n\n\n *************************	 Combination of RCP 	************************* \n\n\n" << endl;
		Int_t bin0PCMRCP =0;
		Int_t bin0EMCalRCP =3;
// 		Int_t bin0EMCalRCP =8;
		graphPCMEtaRCPStat2760GeV->SetName("EtaRCP_pcm");
		TGraphAsymmErrors* graphCombEtaRCPStat2760GeVOEMC= NULL;
		TGraphAsymmErrors* graphCombEtaRCPSys2760GeVOEMC = NULL;
		TGraphAsymmErrors* graphCombEtaRCPTot2760GeVOEMC = CombinePtPointsRAA( 	graphPCMEtaRCPStat2760GeV,	graphPCMEtaRCPSys2760GeV, 	
																				graphEMCalEtaRCPStat2760GeV, graphEMCalEtaRCPSys2760GeV,																graphCombEtaRCPStat2760GeVOEMC, graphCombEtaRCPSys2760GeVOEMC,
																				xPtLimitsEta, 14, /*offSetsEta*/0, bin0PCMRCP, bin0EMCalRCP);
		
// TGraphAsymmErrors* CombinePtPointsRCP( TGraphAsymmErrors* graphStatErrPCM,		TGraphAsymmErrors* graphSystPCM,
// 										  TGraphAsymmErrors* graphStatErrPHOS,		TGraphAsymmErrors* graphSystPHOS,
//  									  TGraphAsymmErrors* &graphStatComb, 	TGraphAsymmErrors* &graphSystComb,  
// 													Double_t* xPtLimits,	Int_t nPtLimits,
// 													Int_t offset,			Int_t bin0PCM,	 Int_t bin0PHOS, Bool_t kRemoveLastPCMPoint=kFALSE, Int_t nBinsPCMRem = 1){

		
// 		cout << "Printing of PCM" << endl;
// 		graphPCMEtaRCPStat2760GeV->Print();
// 		graphPCMEtaRCPSys2760GeV->Print();

// 		cout << "Printing of EMCal" << endl;
// 		graphEMCalEtaRCPStat2760GeV->Print();
// 		graphEMCalEtaRCPSys2760GeV->Print();

		cout << "Printing of combined RCP" << endl;
		graphCombEtaRCPStat2760GeVOEMC->Print();
		graphCombEtaRCPSys2760GeVOEMC->Print();


		TCanvas* canvasRCPcombo = new TCanvas("canvasRCPcombo","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRCPcombo, 0.08, 0.02, 0.035, 0.09);
// 		canvasRCPcombo->SetLogx();
	
		TH2F * histo2DRCPcombo;
		histo2DRCPcombo = new TH2F("histo2DRCPcombo","histo2DRCPcombo",11000,0.23,70.,1000,-0.5,1.1);
		SetStyleHistoTH2ForGraphs(histo2DRCPcombo, "#it{p}_{T} (GeV/#it{c})","#it{R}_{CP}",0.035,0.04, 0.035,0.04, 1.,1.);
// 		histo2DRCPcombo->GetXaxis()->SetMoreLogLabels();	
// 		histo2DRCPcombo->GetXaxis()->SetLabelOffset(-0.01);
		histo2DRCPcombo->GetYaxis()->SetRangeUser(0.,1.);
		histo2DRCPcombo->GetXaxis()->SetRangeUser(0.,30.);
		histo2DRCPcombo->Draw("copy");
		
		DrawGammaSetMarkerTGraphAsym(graphPCMEtaRCPSys2760GeV, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
		graphPCMEtaRCPSys2760GeV->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPCMEtaRCPStat2760GeV,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
		graphPCMEtaRCPStat2760GeV->Draw("p,same,e1");		
		
		DrawGammaSetMarkerTGraphAsym(graphEMCalEtaRCPSys2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
		graphEMCalEtaRCPSys2760GeV->Draw("E2same");			
		DrawGammaSetMarkerTGraphAsym(graphEMCalEtaRCPStat2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
		graphEMCalEtaRCPStat2760GeV->Draw("p,same,e1");
		
		DrawGammaSetMarkerTGraphAsym(graphCombEtaRCPSys2760GeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
		graphCombEtaRCPSys2760GeVOEMC->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombEtaRCPStat2760GeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb);
		graphCombEtaRCPStat2760GeVOEMC->Draw("p,same,e1");

		TLegend* legendRCPcombo = new TLegend(0.62,0.66,0.9,0.86);
		legendRCPcombo->SetFillColor(0);
		legendRCPcombo->SetLineColor(0);
		legendRCPcombo->SetTextFont(42);
		legendRCPcombo->SetTextSize(0.035);
		legendRCPcombo->AddEntry(graphCombEtaRCPStat2760GeVOEMC,"combined","p");		
		legendRCPcombo->AddEntry(graphPCMEtaRCPSys2760GeV,"PCM","p");
		legendRCPcombo->AddEntry(graphEMCalEtaRCPSys2760GeV,"EMCal","p");
		legendRCPcombo->Draw();

		TLatex *labelRCPEnergy = new TLatex(0.7,0.20,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelRCPEnergy, 0.85*textSizeLabelsPixel,4);
		labelRCPEnergy->SetTextFont(43);
// 		labelRCPEnergy->Draw();
		TLatex *labelRCPEta = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRCPEta, 0.85*textSizeLabelsPixel,4);
		labelRCPEta->SetTextFont(43);
// 		labelRCPEta->Draw();

		canvasRCPcombo->SaveAs(Form("%s/%s_%s_RCPEta_combined.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		

		//	**********************************************************************************************************************
		//	*************************	 Combination of EtatoPi0 	**************************************************************
		
		
		TFile* fEtatoPi0input = new TFile("EtaToPi0InputsForCombination.root");
		
		TH1D *histoPCMEtaToPi0RatioPbPb0010 = (TH1D*)fEtatoPi0input->Get("histoPCMEtaToPi0RatioPbPb0010");
		TH1D *histoPCMEtaToPi0RatioPbPb2050 = (TH1D*)fEtatoPi0input->Get("histoPCMEtaToPi0RatioPbPb2050");
		TH1D *histoPCMEtaToPi0RatiopPb = (TH1D*)fEtatoPi0input->Get("histoPCMEtaToPi0RatiopPb");
		TGraphAsymmErrors *graphPCMEtaToPi0RatioSysErrpPb = (TGraphAsymmErrors*)fEtatoPi0input->Get("graphPCMEtaToPi0RatioSysErrpPb");
		
		TGraphAsymmErrors *graphCombEtaToPi0RatioSysErrpp7TeV = (TGraphAsymmErrors*)fEtatoPi0input->Get("graphCombEtaToPi0RatioSysErrpp7TeV");
		TGraphAsymmErrors *graphCombEtaToPi0Ratiopp7TeVNoXErrors = (TGraphAsymmErrors*)fEtatoPi0input->Get("graphCombEtaToPi0Ratiopp7TeVNoXErrors");

		TH1D *cocktailEtaToPi0Ratio7TeVRebined = (TH1D*)fEtatoPi0input->Get("cocktailEtaToPi0Ratio7TeVRebined");
		TH1D *cocktailEtaToPi0Ratio_MtScaledRebinned = (TH1D*)fEtatoPi0input->Get("cocktailEtaToPi0Ratio_MtScaledRebinned");
		TH1D *cocktailEtaToPi0Ratio_K0ScaledRebinned = (TH1D*)fEtatoPi0input->Get("cocktailEtaToPi0Ratio_K0ScaledRebinned");

		cout << "\n\n\n *************************	 Combination of EtatoPi0 	************************* \n\n\n" << endl;
		// Declaration & calculation of combined spectrum				
// 		Int_t bin0PCM =1;
// 		Int_t bin0EMCal =3;
// 		Int_t bin0EMCal =10;
// 		
// 		graphPCMEtatoPi0Stat2760GeV->RemovePoint(0);
// 		graphPCMEtatoPi0Stat2760GeV->RemovePoint(0);
// 		graphPCMEtatoPi0Sys2760GeV->RemovePoint(0);
		cout << "here" << endl;
// 		graphEMCalEtatoPi0Stat2760GeV->Print();
// 		graphEMCalEtatoPi0Sys2760GeV->Print();
		TGraphAsymmErrors* graphCombEtatoPi0Stat2760GeVOEMC= NULL;
		TGraphAsymmErrors* graphCombEtatoPi0Sys2760GeVOEMC = NULL;
		TGraphAsymmErrors* graphCombEtatoPi0Tot2760GeVOEMC = CombinePtPointsSpectraFullCorrMat( 	statErrorCollectionEtatoPi0,	sysErrorCollectionEtatoPi0, 	
																											xPtLimitsEta, 13,
																											offSetsEta, offSetsEtaSys,
																											graphCombEtatoPi0Stat2760GeVOEMC, graphCombEtatoPi0Sys2760GeVOEMC,
																											"weightEtatoPi0.dat",1 );
// 		CombinePtPointsRAA( 	graphPCMEtatoPi0Stat2760GeV,	graphPCMEtatoPi0Sys2760GeV, 	
// 																					graphEMCalEtatoPi0Stat2760GeV, graphEMCalEtatoPi0Sys2760GeV,														graphCombEtatoPi0Stat2760GeVOEMC, graphCombEtatoPi0Sys2760GeVOEMC,
// 																					xPtLimitsEta, 14, /*offSetsEta*/0, bin0PCM, bin0EMCal);
		

		


// 		graphPCMEtatoPi0Stat2760GeV->RemovePoint(0);
// 		graphPCMEtatoPi0Stat2760GeV->RemovePoint(0);
		graphPCMEtatoPi0Sys2760GeV->RemovePoint(0);
		graphPCMEtatoPi0Sys2760GeV->RemovePoint(graphPCMEtatoPi0Sys2760GeV->GetN()-1);
// 		graphPCMEtatoPi0Sys2760GeV->RemovePoint(graphPCMEtatoPi0Sys2760GeV->GetN()-1);
// 		graphPCMEtatoPi0Sys2760GeV->RemovePoint(graphPCMEtatoPi0Sys2760GeV->GetN()-1);
// 		graphPCMEtatoPi0Sys2760GeV->RemovePoint(graphPCMEtatoPi0Sys2760GeV->GetN()-1);
		cout << "Printing of PCM" << endl;
		graphPCMEtatoPi0Stat2760GeV->Print();
		graphPCMEtatoPi0Sys2760GeV->Print();
// 
// 		cout << "Printing of EMCal" << endl;
// 		graphEMCalEtatoPi0Stat2760GeV->Print();
// 		graphEMCalEtatoPi0Sys2760GeV->Print();

		cout << "Printing of combined EtatoPi0" << endl;
		graphCombEtatoPi0Tot2760GeVOEMC->RemovePoint(0);
        graphCombEtatoPi0Stat2760GeVOEMC->RemovePoint(0);
		graphCombEtatoPi0Sys2760GeVOEMC->RemovePoint(0);
//         graphCombEtatoPi0Sys2760GeVOEMC->RemovePoint(0);

		graphCombEtatoPi0Tot2760GeVOEMC->Print();
		graphCombEtatoPi0Stat2760GeVOEMC->Print();
		graphCombEtatoPi0Sys2760GeVOEMC->Print();


		TCanvas* canvasEtatoPi0combo = new TCanvas("canvasEtatoPi0combo","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.08, 0.02, 0.035, 0.09);
		canvasEtatoPi0combo->SetLogx();
	
		TH2F * histo2DEtatoPi0combo;
		histo2DEtatoPi0combo = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",11000,0.23,70.,1000,-0.5,1.1);
		SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();	
// 		histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(-10,10);
		histo2DEtatoPi0combo->Draw("copy");
		histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.,1.2);
		histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(1.,30);
		histo2DEtatoPi0combo->GetYaxis()->SetTitle("#eta/#pi^{0}");

		DrawGammaSetMarkerTGraphAsym(graphPCMEtatoPi0Sys2760GeV, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
		graphPCMEtatoPi0Sys2760GeV->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPCMEtatoPi0Stat2760GeV,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
		graphPCMEtatoPi0Stat2760GeV->Draw("p,same,e1");		
		
		DrawGammaSetMarkerTGraphAsym(graphEMCalEtatoPi0Sys2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
		graphEMCalEtatoPi0Sys2760GeV->Draw("E2same");			
		DrawGammaSetMarkerTGraphAsym(graphEMCalEtatoPi0Stat2760GeV, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
		graphEMCalEtatoPi0Stat2760GeV->Draw("p,same,e1");

		
		DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0Sys2760GeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
		graphCombEtatoPi0Sys2760GeVOEMC->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0Stat2760GeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb);
		graphCombEtatoPi0Stat2760GeVOEMC->Draw("p,same,e1");

		
		TLegend* legendEtatoPi0combo = new TLegend(0.75,0.67,0.95,0.81);
		legendEtatoPi0combo->SetFillColor(0);
		legendEtatoPi0combo->SetLineColor(0);
		legendEtatoPi0combo->SetTextFont(42);
		legendEtatoPi0combo->SetTextSize(0.04);
		legendEtatoPi0combo->AddEntry(graphCombEtatoPi0Sys2760GeVOEMC,"combined","p");
		legendEtatoPi0combo->AddEntry(graphPCMEtatoPi0Sys2760GeV,"PCM","p");
		legendEtatoPi0combo->AddEntry(graphEMCalEtatoPi0Sys2760GeV,"EMCal","p");
		legendEtatoPi0combo->Draw();
		
// 		cocktailEtaToPi0Ratio7TeVRebined->GetXaxis()->SetRangeUser(0.2,20);
// 		cocktailEtaToPi0Ratio_K0ScaledRebinned->GetXaxis()->SetRangeUser(0.2,20);
// 		cocktailEtaToPi0Ratio_MtScaledRebinned->GetXaxis()->SetRangeUser(0.2,20);
		cocktailEtaToPi0Ratio7TeVRebined->SetLineStyle(5);
		cocktailEtaToPi0Ratio7TeVRebined->SetLineWidth(2.5);
		cocktailEtaToPi0Ratio_K0ScaledRebinned->SetLineStyle(6);
		cocktailEtaToPi0Ratio_K0ScaledRebinned->SetLineWidth(2.5);
		cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineStyle(7);
		cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineWidth(2.5);
		cocktailEtaToPi0Ratio7TeVRebined->Draw("same,hist,c");
		cocktailEtaToPi0Ratio_K0ScaledRebinned->Draw("same,hist,c");
		cocktailEtaToPi0Ratio_MtScaledRebinned->Draw("same,hist,c");

		TLegend* legendRatioALICE2 = new TLegend(0.3,0.16,0.93,0.3);
		legendRatioALICE2->SetTextSize(0.04);			
		legendRatioALICE2->SetFillColor(0);
		legendRatioALICE2->SetFillStyle(0);
		legendRatioALICE2->SetBorderSize(0);
		legendRatioALICE2->SetMargin(0.1);
		legendRatioALICE2->AddEntry(cocktailEtaToPi0Ratio7TeVRebined,Form("#eta from %s as input",collisionSystemPP7TeV.Data()),"pl");
		legendRatioALICE2->AddEntry(cocktailEtaToPi0Ratio_MtScaledRebinned,"#eta from m_{T} scaled #pi^{0}","pl");
		legendRatioALICE2->AddEntry(cocktailEtaToPi0Ratio_K0ScaledRebinned,"0-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV, #eta from K^{0}_{s} scaled ","pl");
		legendRatioALICE2->Draw();


		TLatex *labeltoPi0Energy = new TLatex(0.55,0.87,collisionSystem2760GeV.Data());
		SetStyleTLatex( labeltoPi0Energy, 0.85*textSizeLabelsPixel,4);
		labeltoPi0Energy->SetTextFont(43);
		labeltoPi0Energy->Draw();
		TLatex *labeltoPi0Eta = new TLatex(0.55,0.81,"#eta #rightarrow #gamma#gamma");
		SetStyleTLatex( labeltoPi0Eta, 0.85*textSizeLabelsPixel,4);
		labeltoPi0Eta->SetTextFont(43);
		labeltoPi0Eta->Draw();


		canvasEtatoPi0combo->SaveAs(Form("%s/%s_%s_EtatoPi0Ratio_combined.%s",outputDir.Data(),meson.Data(),centrality.Data(),suffix.Data()));
		
						
		TCanvas* canvasRatioEtaToPi0ALICEPbPb = new TCanvas("canvasRatioEtaToPi0ALICEPbPb","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRatioEtaToPi0ALICEPbPb, 0.09, 0.01, 0.015, 0.115);
		
		TH2D *histo2DRatioEtaToPi0ALICEPbPb = new TH2D("histo2DRatioEtaToPi0ALICEPbPb", "histo2DRatioEtaToPi0ALICEPbPb", 20,0.,20.,1000.,-0.4,1.3);
		SetStyleHistoTH2ForGraphs(histo2DRatioEtaToPi0ALICEPbPb, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
		histo2DRatioEtaToPi0ALICEPbPb->GetYaxis()->SetRangeUser(0.,1.32);
		histo2DRatioEtaToPi0ALICEPbPb->Draw("copy");

	// 	histo2DRatioEtaToPi0ALICE->Draw("copy");
		
		graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,pE2");
		graphCombEtaToPi0Ratiopp7TeVNoXErrors->Draw("same,pe");
	// 	cocktailEtaToPi0Ratio7TeVRebined->GetXaxis()->SetRangeUser(0.2,20);
	// 	cocktailEtaToPi0Ratio_K0ScaledRebinned->GetXaxis()->SetRangeUser(0.2,20);
	// 	cocktailEtaToPi0Ratio_MtScaledRebinned->GetXaxis()->SetRangeUser(0.2,20);
		cocktailEtaToPi0Ratio7TeVRebined->SetLineStyle(5);
		cocktailEtaToPi0Ratio7TeVRebined->SetLineWidth(2.5);
		cocktailEtaToPi0Ratio_K0ScaledRebinned->SetLineStyle(6);
		cocktailEtaToPi0Ratio_K0ScaledRebinned->SetLineWidth(2.5);
		cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineStyle(7);
		cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineWidth(2.5);
		cocktailEtaToPi0Ratio7TeVRebined->Draw("same,hist,c");
		cocktailEtaToPi0Ratio_K0ScaledRebinned->Draw("same,hist,c");
		cocktailEtaToPi0Ratio_MtScaledRebinned->Draw("same,hist,c");
		
	// 	graphPCMEtaToPi0RatioSysErrPbPb2050->Draw("same,pE2");  
	// 	graphPCMEtaToPi0RatioSysErrPbPb0010->Draw("same,pE2");
	// 	histoPCMEtaToPi0RatioPbPb2050->Draw("same,peX0");
	// 	histoPCMEtaToPi0RatioPbPb0010->Draw("same,peX0");
		if (centrality.CompareTo("0010")==0){
			DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0Sys2760GeVOEMC, markerStyleComb, markerSizeComb, kRed+1 , kRed+1, widthLinesBoxes, kTRUE);
			DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0Stat2760GeVOEMC, markerStyleComb, markerSizeComb, kRed+1 , kRed+1);
		} else if (centrality.CompareTo("2050")==0){
			DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0Sys2760GeVOEMC, markerStyleComb, markerSizeComb, kAzure+7 , kAzure+7, widthLinesBoxes, kTRUE);
			DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0Stat2760GeVOEMC, markerStyleComb, markerSizeComb, kAzure+7 , kAzure+7);			
		}
		graphCombEtatoPi0Stat2760GeVOEMC->Draw("same,peX0");
		graphCombEtatoPi0Sys2760GeVOEMC->Draw("same,pE2");
		histoPCMEtaToPi0RatiopPb->Draw("same,peX0");
		graphPCMEtaToPi0RatioSysErrpPb->Draw("same,pE2");

		TLegend* legendRatioALICEdata2 = new TLegend(0.12,0.72,0.32,0.96);
		legendRatioALICEdata2->SetTextSize(0.04);			
		legendRatioALICEdata2->SetFillColor(0);
		legendRatioALICEdata2->SetFillStyle(0);
		legendRatioALICEdata2->SetBorderSize(0);
	// 	legendRatioALICEdata2->SetMargin(0.1);
		legendRatioALICEdata2->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,Form("ALICE, %s",collisionSystemPP7TeV.Data()),"pe");
		legendRatioALICEdata2->AddEntry((TObject*)0,"Phys.Lett. B717 (2012) 162-172","");
	// 	legendRatioALICEdata2->AddEntry(histoPCMEtaToPi0RatioPbPb0010,Form("%s, LHC11h",collisionSystemPbPb0010.Data()),"pe");
	// 	legendRatioALICEdata2->AddEntry(histoPCMEtaToPi0RatioPbPb2050,Form("%s, LHC11h",collisionSystemPbPb2050.Data()),"pe");
		if (centrality.CompareTo("0010")==0){
			legendRatioALICEdata2->AddEntry(graphCombEtatoPi0Sys2760GeVOEMC,Form("%s, LHC11h",collisionSystemPbPb0010.Data()),"pe");
		} else if (centrality.CompareTo("2050")==0){
			legendRatioALICEdata2->AddEntry(graphCombEtatoPi0Sys2760GeVOEMC,Form("%s, LHC11h",collisionSystemPbPb2050.Data()),"pe");
		}
		legendRatioALICEdata2->AddEntry(histoPCMEtaToPi0RatiopPb,Form("%s, LHC13b+c",collisionSystempPb.Data()),"pe");
		legendRatioALICEdata2->Draw();
		
// 		TLegend* legendRatioALICE2 = new TLegend(0.3,0.16,0.977,0.3);
// 		legendRatioALICE2->SetTextSize(0.04);			
// 		legendRatioALICE2->SetFillColor(0);
// 		legendRatioALICE2->SetFillStyle(0);
// 		legendRatioALICE2->SetBorderSize(0);
// 		legendRatioALICE2->SetMargin(0.1);
// 		legendRatioALICE2->AddEntry(cocktailEtaToPi0Ratio7TeVRebined,Form("#eta from %s as input",collisionSystemPP7TeV.Data()),"pl");
// 		legendRatioALICE2->AddEntry(cocktailEtaToPi0Ratio_MtScaledRebinned,"#eta from m_{T} scaled #pi^{0}","pl");
// 		legendRatioALICE2->AddEntry(cocktailEtaToPi0Ratio_K0ScaledRebinned,"0-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV, #eta from K^{0}_{s} scaled ","pl");
		legendRatioALICE2->Draw();
		
		canvasRatioEtaToPi0ALICEPbPb->Update();  //qui
		canvasRatioEtaToPi0ALICEPbPb->SaveAs(Form("%s/EtaToPi0Ratio_DiffSystems_PbPb%s.%s",outputDir.Data(), centrality.Data(), suffix.Data()));
		
		
		
		//	**********************************************************************************************************************
		//	************************* Saving of final results ********************************************************************
		//	**********************************************************************************************************************
		
// 		cout <<"method A" << endl;
// 		graphCombEtaInvYieldStat2760GeVA->Print();
// 		cout <<"method B" << endl;
// 		graphCombEtaInvYieldStat2760GeVB->Print();
		TFile fCombResults(Form("%s/%s_%s_CombinedResultsPaperPbPb2760GeV_%s.root", outputDir.Data(),meson.Data(),centrality.Data(),dateForOutput.Data()), "RECREATE");
			graphCombEtaInvYieldTot2760GeVA->Write("graphInvYieldEtaComb2760GeVA");
			graphCombEtaInvYieldStat2760GeVA->Write("graphInvYieldEtaComb2760GeVAStatErr");
			graphCombEtaInvYieldSys2760GeVA->Write("graphInvYieldEtaComb2760GeVASysErr");      
			graphCombEtaInvYieldTot2760GeVB->Write("graphInvYieldEtaComb2760GeVB");
			graphCombEtaInvYieldStat2760GeVB->Write("graphInvYieldEtaComb2760GeVBStatErr");
			graphCombEtaInvYieldSys2760GeVB->Write("graphInvYieldEtaComb2760GeVBSysErr");      
		fCombResults.Close();

		
		
		
	}
	

}
	
