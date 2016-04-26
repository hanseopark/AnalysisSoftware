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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"
#include "CombineMesonMeasurementsPbPbVer2.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

void CombineMesonMeasurementsPbPbVer2(TString fileNameConversionPbPb = "", TString fileNameCombPP2760GeV = "", TString suffix = "eps", TString conferencePlots ="", TString bWCorrection="X"){	

	date = ReturnDateString();
	gROOT->Reset();	
	gROOT->SetStyle("Plain");

	StyleSettingsThesis();	
	SetPlotStyle();

	//******************************** Declaration of files *****************************************
	collisionSystemPP = "pp #sqrt{#it{s}} = 2.76 TeV";		
	collisionSystemPbPb = "Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";     
	collisionSystemCent0005 = "0-5% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemCent0005Legend = "  0-  5% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";    
	collisionSystemCent0510 = "5-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	TString collisionSystemCent0510Legend = "  5-10% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";     
	collisionSystemCent1020 = "10-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	collisionSystemCent0020 = "0-20% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	collisionSystemCent2040 = "20-40% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	collisionSystemCent4060 = "40-60% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		
	collisionSystemCent6080 = "60-80% Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";		

	fileNameCaloPHOSPbPb = "ExternalInputPbPb/PHOS/LHC10h_PHOS_pi0_PbPb_06022014.root";
	fileNameCaloPHOSPP = "ExternalInput/PHOS/2.76TeV/LHC11a_PHOS_pi0_pp2760_noBWCorr_FDcorr_20131112_v2.root";
	cout << fileNameCaloPHOSPbPb << endl;

	fileNameCaloEMCALPbPb = "ExternalInputPbPb/EMCAL/EMCALResults_PbPb.root";

	fileNameTheoryInput = "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";
	fileNameDataOtherEnergyInput = "ExternalInputPbPb/OtherExperiments/DataCompilationFromOtherEnergiesPbPb.root";

	fileNamePCMPP2760GeV =  "data_PCMResultsFullCorrection_PP_NoBinShifting_4_Aug_2013.root";
	TString dateForOutput = ReturnDateStringForOutput();
	TString outputDir = Form("%s/%s/CombineMesonMeasurementsPbPb%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
	nameFinalResDat = Form("%s/CombinedResultsPbPb%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	gSystem->Exec(Form("cp %s %s/InputPCMPbPb.root", fileNameConversionPbPb.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPCMPP.root", fileNamePCMPP2760GeV.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputCombPP.root", fileNameCombPP2760GeV.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPHOSPbPb.root", fileNameCaloPHOSPbPb.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPHOSPP.root", fileNameCaloPHOSPP.Data(), outputDir.Data()));

	TString nameHistoConv = "CorrectedYieldPi0";
	TString nameGraphConv = "Pi0SystError";
	TString namePPComb;
	TString namePPCombSys;

	fileFinalResultsPP = 		new TFile(fileNameCombPP2760GeV);
	graphInvSectionCombStatPi02760GeVPlot = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0Comb2760GeVStatErr");
	graphInvSectionCombStatPi02760GeVPlot = ScaleGraph(graphInvSectionCombStatPi02760GeVPlot,1./xSection2769GeVppINEL);
	graphInvSectionCombStatPi02760GeVPlot->RemovePoint(graphInvSectionCombStatPi02760GeVPlot->GetN()-1);
	graphInvSectionCombSysPi02760GeVPlot = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0Comb2760GeVSysErr");
	graphInvSectionCombSysPi02760GeVPlot = ScaleGraph(graphInvSectionCombSysPi02760GeVPlot,1./xSection2769GeVppINEL);
	graphInvSectionCombSysPi02760GeVPlot->RemovePoint(graphInvSectionCombSysPi02760GeVPlot->GetN()-1);

	graphInvSectionPCMStatPi02760GeVPlot = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
	graphInvSectionPCMStatPi02760GeVPlot = ScaleGraph(graphInvSectionPCMStatPi02760GeVPlot,1./xSection2769GeVppINEL);
	graphInvSectionPCMSysPi02760GeVPlot = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0PCM2760GeVSysErr");
	graphInvSectionPCMSysPi02760GeVPlot = ScaleGraph(graphInvSectionPCMSysPi02760GeVPlot,1./xSection2769GeVppINEL);

	graphInvSectionPHOSStatPi02760GeVPlot = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0PHOS2760GeVStatErr");
	graphInvSectionPHOSStatPi02760GeVPlot = ScaleGraph(graphInvSectionPHOSStatPi02760GeVPlot,1./xSection2769GeVppINEL);
	graphInvSectionPHOSSysPi02760GeVPlot = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0PHOS2760GeVSysErr");
	graphInvSectionPHOSSysPi02760GeVPlot = ScaleGraph(graphInvSectionPHOSSysPi02760GeVPlot,1./xSection2769GeVppINEL);

	graphInvSectionCombStatPi02760GeV = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0Comb2760GeVStatErr_YShifted");
	graphInvSectionCombStatPi02760GeV = ScaleGraph(graphInvSectionCombStatPi02760GeV,1./xSection2769GeVppINEL);
	graphInvSectionCombStatPi02760GeV->RemovePoint(graphInvSectionCombStatPi02760GeV->GetN()-1);
	graphInvSectionCombSysPi02760GeV = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0Comb2760GeVSysErr_YShifted");
	graphInvSectionCombSysPi02760GeV = ScaleGraph(graphInvSectionCombSysPi02760GeV,1./xSection2769GeVppINEL);
	graphInvSectionCombSysPi02760GeV->RemovePoint(graphInvSectionCombStatPi02760GeV->GetN()-1);
	graphInvSectionPHOSSysPi02760GeV = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0PHOSSysForRAA2760GeV_YShifted");
	graphInvSectionPHOSSysPi02760GeV = ScaleGraph(graphInvSectionPHOSSysPi02760GeV,1./xSection2769GeVppINEL);
	graphInvSectionPHOSStatPi02760GeV = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0PHOSStat2760GeV_YShifted");
	graphInvSectionPHOSStatPi02760GeV = ScaleGraph(graphInvSectionPHOSStatPi02760GeV,1./xSection2769GeVppINEL);
	graphInvSectionCombPi02760GeV = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0Comb2760GeV_YShifted");
	graphInvSectionCombPi02760GeV = ScaleGraph(graphInvSectionCombPi02760GeV,1./xSection2769GeVppINEL);
	TGraphAsymmErrors* graphInvSectionCombPi02760GeVOnlyStat = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0Comb2760GeVStatErr_YShifted");
	graphInvSectionCombPi02760GeVOnlyStat = ScaleGraph(graphInvSectionCombPi02760GeVOnlyStat,1./xSection2769GeVppINEL);

	graphInvSectionPCMSysPi02760GeV = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0PCMSysForRAA2760GeV_YShifted");
	graphInvSectionPCMSysPi02760GeV = ScaleGraph(graphInvSectionPCMSysPi02760GeV,1./xSection2769GeVppINEL);
	graphInvSectionPCMStatPi02760GeV = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphInvCrossSectionPi0PCMStat2760GeV_YShifted");
	graphInvSectionPCMStatPi02760GeV = ScaleGraph(graphInvSectionPCMStatPi02760GeV,1./xSection2769GeVppINEL);

	TGraphAsymmErrors* graphRatioCombCombFit2760GeVSys = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphRatioCombCombFit2760GeVSys");
	TGraphAsymmErrors* graphRatioCombCombFit2760GeVStat = (TGraphAsymmErrors*)fileFinalResultsPP->Get("graphRatioCombCombFit2760GeVStat");
	TGraph* graphRatioCombNLOPi02760GeVMuHalf = (TGraph*)fileFinalResultsPP->Get("graphRatioCombNLOPi02760GeVMuHalf");
	TGraph* graphRatioCombNLOPi02760GeVMuOne = (TGraph*)fileFinalResultsPP->Get("graphRatioCombNLOPi02760GeVMuOne");
	TGraph* graphRatioCombNLOPi02760GeVMuTwo = (TGraph*)fileFinalResultsPP->Get("graphRatioCombNLOPi02760GeVMuTwo");
	TH1F* histoRatioPythia8ToFit2760GeV = (TH1F*)fileFinalResultsPP->Get("histoRatioPythia8ToFit2760GeV");

	graphInvSectionPCMStatPi02760GeVRed = new TGraphAsymmErrors(17);
	graphInvSectionPCMSysPi02760GeVRed = new TGraphAsymmErrors(17);
	Double_t* xValPCMStat = graphInvSectionPCMStatPi02760GeV->GetX();
	Double_t* yValPCMStat = graphInvSectionPCMStatPi02760GeV->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphInvSectionPCMStatPi02760GeVRed->SetPoint(i,xValPCMStat[i], yValPCMStat[i]);
		graphInvSectionPCMStatPi02760GeVRed->SetPointError(i, graphInvSectionPCMStatPi02760GeV->GetErrorXlow(i),graphInvSectionPCMStatPi02760GeV->GetErrorXhigh(i), graphInvSectionPCMStatPi02760GeV->GetErrorYlow(i), graphInvSectionPCMStatPi02760GeV->GetErrorYhigh(i));
		graphInvSectionPCMSysPi02760GeVRed->SetPoint(i,xValPCMStat[i], yValPCMStat[i]);
		graphInvSectionPCMSysPi02760GeVRed->SetPointError(i, graphInvSectionPCMSysPi02760GeV->GetErrorXlow(i),graphInvSectionPCMSysPi02760GeV->GetErrorXhigh(i), graphInvSectionPCMSysPi02760GeV->GetErrorYlow(i), graphInvSectionPCMSysPi02760GeV->GetErrorYhigh(i));
		cout << "PCM: \t" << i << "\t" << xValPCMStat[i] << "\t" << yValPCMStat[i] << endl;
	}

	graphInvSectionPHOSStatPi02760GeVRed = new TGraphAsymmErrors(17);
	graphInvSectionPHOSSysPi02760GeVRed = new TGraphAsymmErrors(17);
	Double_t* xValPHOSStat = graphInvSectionPHOSStatPi02760GeV->GetX();
	Double_t* yValPHOSStat = graphInvSectionPHOSStatPi02760GeV->GetY();
	// 	graphInvSectionPHOSStatPi02760GeV->Print();
	graphInvSectionPHOSSysPi02760GeV->RemovePoint(0);
	// 	graphInvSectionPHOSSysPi02760GeV->Print();
	Int_t nOffsetPhosPP = 0;
	for (Int_t i = 0; i < 17; i++){
		graphInvSectionPHOSStatPi02760GeVRed->SetPoint(i,xValPHOSStat[i+nOffsetPhosPP], yValPHOSStat[i+nOffsetPhosPP]);
		graphInvSectionPHOSStatPi02760GeVRed->SetPointError(i, graphInvSectionPHOSStatPi02760GeV->GetErrorXlow(i+nOffsetPhosPP),graphInvSectionPHOSStatPi02760GeV->GetErrorXhigh(i+nOffsetPhosPP), graphInvSectionPHOSStatPi02760GeV->GetErrorYlow(i+nOffsetPhosPP), graphInvSectionPHOSStatPi02760GeV->GetErrorYhigh(i+nOffsetPhosPP));
		graphInvSectionPHOSSysPi02760GeVRed->SetPoint(i,xValPHOSStat[i+nOffsetPhosPP], yValPHOSStat[i+nOffsetPhosPP]);
		graphInvSectionPHOSSysPi02760GeVRed->SetPointError(i, graphInvSectionPHOSSysPi02760GeV->GetErrorXlow(i+nOffsetPhosPP),graphInvSectionPHOSSysPi02760GeV->GetErrorXhigh(i+nOffsetPhosPP), graphInvSectionPHOSSysPi02760GeV->GetErrorYlow(i+nOffsetPhosPP), graphInvSectionPHOSSysPi02760GeV->GetErrorYhigh(i+nOffsetPhosPP));
		cout << "PHOS: \t" << i << "\t" << xValPHOSStat[i+nOffsetPhosPP] << "\t" << yValPHOSStat[i+nOffsetPhosPP] << "\t" << graphInvSectionPHOSSysPi02760GeV->GetErrorYlow(i+nOffsetPhosPP)<< endl;
	}

	filePCMPP = 		new TFile(fileNamePCMPP2760GeV);
	directoryPi0PP = 			(TDirectory*)filePCMPP->Get("Pi02.76TeV"); 
	histoPCMMassDataPP = 				(TH1D*)directoryPi0PP->Get("MassPi0");
	histoPCMWidthDataPP = 				(TH1D*)directoryPi0PP->Get("FWHMPi0MeV");
	histoPCMMassMCPP = 			(TH1D*)directoryPi0PP->Get("TrueMassPi0");
	histoPCMWidthMCPP = 			(TH1D*)directoryPi0PP->Get("TrueFWHMPi0MeV");
	histoPCMMassDataPP->Scale(1000.);
	histoPCMMassMCPP->Scale(1000.);
	histoPCMMassDataPP->SetBinContent(histoPCMMassDataPP->GetNbinsX(),0);
	histoPCMMassMCPP->SetBinContent(histoPCMMassMCPP->GetNbinsX(),0);
	histoPCMWidthDataPP->SetBinContent(histoPCMWidthDataPP->GetNbinsX(),10000.);
	histoPCMWidthMCPP->SetBinContent(histoPCMWidthMCPP->GetNbinsX(),10000.);

	minPtForFits = 0.4;
	if(conferencePlots.CompareTo("conference") == 0){
		conference = kTRUE;
	}	
	fstream  fileFinalResultsFits;
	fileFinalResultsFits.open(nameFinalResDat.Data(), ios::out);  

	//declaration for printing logo 	
	prefix2 = "data";
	pictDrawingOptions[1] = kFALSE;

	mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	cout << mesonMassExpectPi0 << endl;
	//************************** Read data for conversions **************************************************
	// File definitions
	filePCMPbPb = 					new TFile(fileNameConversionPbPb.Data());
	directoryPi0PbPb0005 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_0-5%"); 
	cout << "PCM: 0-5%" << endl;
	histoNumberOfEvents0005= 			(TH1D*)filePCMPbPb->Get("histoNumberOfEventsPbPb_2.76TeV0-5%");
	histoYieldPbPb0005 = 				(TH1D*)directoryPi0PbPb0005->Get(nameHistoConv.Data());
	graphPCMYieldPi0SysErrPbPb0005= 	(TGraphAsymmErrors*)directoryPi0PbPb0005->Get(nameGraphConv.Data());	
	graphPCMYieldPi0SysErrRAAPbPb0005= 	(TGraphAsymmErrors*)directoryPi0PbPb0005->Get("Pi0SystErrorA");	
	histoPCMMassData0005 = 			(TH1D*)directoryPi0PbPb0005->Get("MassPi0");
	histoPCMMassMC0005 = 			(TH1D*)directoryPi0PbPb0005->Get("TrueMassPi0");
	histoPCMWidthData0005 = 		(TH1D*)directoryPi0PbPb0005->Get("FWHMPi0MeV");
	histoPCMWidthMC0005 = 			(TH1D*)directoryPi0PbPb0005->Get("TrueFWHMPi0MeV");
	histoPCMMassData0005->Scale(1000.);
	histoPCMMassMC0005->Scale(1000.);
	histoPCMMassData0005->SetBinContent(histoPCMMassData0005->GetNbinsX(),0);
	histoPCMMassData0005->SetBinContent(histoPCMMassData0005->GetNbinsX()-1,0);
	histoPCMMassMC0005->SetBinContent(histoPCMMassMC0005->GetNbinsX(),0);
	histoPCMMassMC0005->SetBinContent(histoPCMMassMC0005->GetNbinsX()-1,0);
	histoPCMWidthData0005->SetBinContent(histoPCMWidthData0005->GetNbinsX(),10000.);
	histoPCMWidthData0005->SetBinContent(histoPCMWidthData0005->GetNbinsX()-1,10000.);
	histoPCMWidthMC0005->SetBinContent(histoPCMWidthMC0005->GetNbinsX(),10000.);
	histoPCMWidthMC0005->SetBinContent(histoPCMWidthMC0005->GetNbinsX()-1,10000.);
	for (Int_t  i = 0; i < histoPCMMassData0005->GetXaxis()->FindBin(0.6); i++){
		histoPCMMassData0005->SetBinContent(i,0);
		histoPCMMassMC0005->SetBinContent(i,0);
		histoPCMWidthData0005->SetBinContent(i,10000.);  
		histoPCMWidthMC0005->SetBinContent(i,10000.);  
	}

	nEvent0005 = 							histoNumberOfEvents0005->GetBinContent(1);

	graphPCMYieldPi0SysErrPbPb0005Red = new TGraphAsymmErrors(17);
	graphPCMYieldPi0SysErrRAAPbPb0005Red = new TGraphAsymmErrors(17);
	Double_t* xValPCMStat0005 = graphPCMYieldPi0SysErrPbPb0005->GetX();
	Double_t* yValPCMStat0005 = graphPCMYieldPi0SysErrPbPb0005->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPCMYieldPi0SysErrPbPb0005Red->SetPoint(i,xValPCMStat0005[i], yValPCMStat0005[i]);
		graphPCMYieldPi0SysErrPbPb0005Red->SetPointError(i, graphPCMYieldPi0SysErrPbPb0005->GetErrorXlow(i),graphPCMYieldPi0SysErrPbPb0005->GetErrorXhigh(i), graphPCMYieldPi0SysErrPbPb0005->GetErrorYlow(i), graphPCMYieldPi0SysErrPbPb0005->GetErrorYhigh(i));
		graphPCMYieldPi0SysErrRAAPbPb0005Red->SetPoint(i,xValPCMStat0005[i], yValPCMStat0005[i]);
		graphPCMYieldPi0SysErrRAAPbPb0005Red->SetPointError(i, graphPCMYieldPi0SysErrRAAPbPb0005->GetErrorXlow(i),graphPCMYieldPi0SysErrRAAPbPb0005->GetErrorXhigh(i), graphPCMYieldPi0SysErrRAAPbPb0005->GetErrorYlow(i), graphPCMYieldPi0SysErrRAAPbPb0005->GetErrorYhigh(i));
	}
	graphPCMYieldPi0SysErrPbPb0005Red->SetPoint(0,0.5, 0);
	graphPCMYieldPi0SysErrRAAPbPb0005Red->SetPoint(0,0.5, 0);
	histoYieldPbPb0005Red = 		new TH1D("histoYieldPbPb0005Red","",18,fBinsPi0HIPtPCM);
	for (Int_t i =1; i < 19; i++){
		cout << i<< "\t"<<histoYieldPbPb0005->GetBinCenter(i)<< "\t"<<histoYieldPbPb0005->GetBinContent(i) << endl;
		histoYieldPbPb0005Red->SetBinContent(i,histoYieldPbPb0005->GetBinContent(i));
		histoYieldPbPb0005Red->SetBinError(i,histoYieldPbPb0005->GetBinError(i));
	}
	histoYieldPbPb0005Red->SetBinContent(1,0);
	histoYieldPbPb0005Red->SetBinContent(2,0);
	nPointsPi00005 = graphPCMYieldPi0SysErrPbPb0005Red->GetN();

	directoryPi0PbPb0010 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_0-10%"); 
	cout << "PCM: 0-10%" << endl;
	histoNumberOfEvents0010= 			(TH1D*)filePCMPbPb->Get("histoNumberOfEventsPbPb_2.76TeV0-10%");
	histoYieldPbPb0010 = 				(TH1D*)directoryPi0PbPb0010->Get(nameHistoConv.Data());
	graphPCMYieldPi0SysErrPbPb0010= 	(TGraphAsymmErrors*)directoryPi0PbPb0010->Get(nameGraphConv.Data());	
	graphPCMYieldPi0SysErrRAAPbPb0010= 	(TGraphAsymmErrors*)directoryPi0PbPb0010->Get("Pi0SystErrorA");	
	histoPCMMassData0010 =        (TH1D*)directoryPi0PbPb0010->Get("MassPi0");
	histoPCMMassMC0010 =          (TH1D*)directoryPi0PbPb0010->Get("TrueMassPi0");
	histoPCMWidthData0010 =       (TH1D*)directoryPi0PbPb0010->Get("FWHMPi0MeV");
	histoPCMWidthMC0010 =         (TH1D*)directoryPi0PbPb0010->Get("TrueFWHMPi0MeV");
	histoPCMMassData0010->Scale(1000.);
	histoPCMMassMC0010->Scale(1000.);
	for (Int_t  i = 0; i < histoPCMMassMC0010->GetXaxis()->FindBin(0.6); i++){
		histoPCMMassData0010->SetBinContent(i,0);
		histoPCMMassMC0010->SetBinContent(i,0);
		histoPCMWidthData0010->SetBinContent(i,10000.);  
		histoPCMWidthMC0010->SetBinContent(i,10000.);  
	}

	histoPCMMassData0010->SetBinContent(histoPCMMassData0010->GetNbinsX(),0);
	histoPCMMassData0010->SetBinContent(histoPCMMassData0010->GetNbinsX()-1,0);
	histoPCMMassMC0010->SetBinContent(histoPCMMassMC0010->GetNbinsX(),0);
	histoPCMMassMC0010->SetBinContent(histoPCMMassMC0010->GetNbinsX()-1,0);
	histoPCMWidthData0010->SetBinContent(histoPCMWidthData0010->GetNbinsX(),10000.);
	histoPCMWidthData0010->SetBinContent(histoPCMWidthData0010->GetNbinsX()-1,10000.);
	histoPCMWidthMC0010->SetBinContent(histoPCMWidthMC0010->GetNbinsX(),10000.);
	histoPCMWidthMC0010->SetBinContent(histoPCMWidthMC0010->GetNbinsX()-1,10000.);
	nEvent0010 = 							histoNumberOfEvents0010->GetBinContent(1);

	graphPCMYieldPi0SysErrPbPb0010Red = new TGraphAsymmErrors(17);
	graphPCMYieldPi0SysErrRAAPbPb0010Red = new TGraphAsymmErrors(17);
	Double_t* xValPCMStat0010 = graphPCMYieldPi0SysErrPbPb0010->GetX();
	Double_t* yValPCMStat0010 = graphPCMYieldPi0SysErrPbPb0010->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPCMYieldPi0SysErrPbPb0010Red->SetPoint(i,xValPCMStat0010[i], yValPCMStat0010[i]);
		graphPCMYieldPi0SysErrPbPb0010Red->SetPointError(i, graphPCMYieldPi0SysErrPbPb0010->GetErrorXlow(i),graphPCMYieldPi0SysErrPbPb0010->GetErrorXhigh(i), graphPCMYieldPi0SysErrPbPb0010->GetErrorYlow(i), graphPCMYieldPi0SysErrPbPb0010->GetErrorYhigh(i));
		graphPCMYieldPi0SysErrRAAPbPb0010Red->SetPoint(i,xValPCMStat0010[i], yValPCMStat0010[i]);
		graphPCMYieldPi0SysErrRAAPbPb0010Red->SetPointError(i, graphPCMYieldPi0SysErrRAAPbPb0010->GetErrorXlow(i),graphPCMYieldPi0SysErrRAAPbPb0010->GetErrorXhigh(i), graphPCMYieldPi0SysErrRAAPbPb0010->GetErrorYlow(i), graphPCMYieldPi0SysErrRAAPbPb0010->GetErrorYhigh(i));
	}
	graphPCMYieldPi0SysErrPbPb0010Red->SetPoint(0,0.5, 0);
	graphPCMYieldPi0SysErrRAAPbPb0010Red->SetPoint(0,0.5, 0);
	histoYieldPbPb0010Red = 		new TH1D("histoYieldPbPb0010Red","",18,fBinsPi0HIPtPCM);
	for (Int_t i =1; i < 19; i++){
		histoYieldPbPb0010Red->SetBinContent(i,histoYieldPbPb0010->GetBinContent(i));
		histoYieldPbPb0010Red->SetBinError(i,histoYieldPbPb0010->GetBinError(i));
	}
	histoYieldPbPb0010Red->SetBinContent(1,0);
	histoYieldPbPb0010Red->SetBinContent(2,0);
	nPointsPi00010 = graphPCMYieldPi0SysErrPbPb0010Red->GetN();


	cout << "PCM: 5-10%" << endl;
	directoryPi0PbPb0510 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_5-10%"); 
	histoNumberOfEvents0510= 			(TH1D*)filePCMPbPb->Get("histoNumberOfEventsPbPb_2.76TeV5-10%");
	histoYieldPbPb0510 = 				(TH1D*)directoryPi0PbPb0510->Get(nameHistoConv.Data());
	graphPCMYieldPi0SysErrPbPb0510= 	(TGraphAsymmErrors*)directoryPi0PbPb0510->Get(nameGraphConv.Data());	
	graphPCMYieldPi0SysErrRAAPbPb0510= 	(TGraphAsymmErrors*)directoryPi0PbPb0510->Get("Pi0SystErrorA");	
	histoPCMMassData0510 = 			(TH1D*)directoryPi0PbPb0510->Get("MassPi0");
	histoPCMMassMC0510 = 			(TH1D*)directoryPi0PbPb0510->Get("TrueMassPi0");
	histoPCMWidthData0510 = 		(TH1D*)directoryPi0PbPb0510->Get("FWHMPi0MeV");
	histoPCMWidthMC0510 = 			(TH1D*)directoryPi0PbPb0510->Get("TrueFWHMPi0MeV");
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

	nEvent0510 = 							histoNumberOfEvents0510->GetBinContent(1);

	//  	graphPCMYieldPi0SysErrRAAPbPb0510->Print();
	graphPCMYieldPi0SysErrPbPb0510Red = new TGraphAsymmErrors(17);
	graphPCMYieldPi0SysErrRAAPbPb0510Red = new TGraphAsymmErrors(17);
	Double_t* xValPCMStat0510 = graphPCMYieldPi0SysErrPbPb0510->GetX();
	Double_t* yValPCMStat0510 = graphPCMYieldPi0SysErrPbPb0510->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPCMYieldPi0SysErrPbPb0510Red->SetPoint(i,xValPCMStat0510[i], yValPCMStat0510[i]);
		graphPCMYieldPi0SysErrPbPb0510Red->SetPointError(i, graphPCMYieldPi0SysErrPbPb0510->GetErrorXlow(i),graphPCMYieldPi0SysErrPbPb0510->GetErrorXhigh(i), graphPCMYieldPi0SysErrPbPb0510->GetErrorYlow(i), graphPCMYieldPi0SysErrPbPb0510->GetErrorYhigh(i));
		graphPCMYieldPi0SysErrRAAPbPb0510Red->SetPoint(i,xValPCMStat0510[i], yValPCMStat0510[i]);
		graphPCMYieldPi0SysErrRAAPbPb0510Red->SetPointError(i, graphPCMYieldPi0SysErrRAAPbPb0510->GetErrorXlow(i),graphPCMYieldPi0SysErrRAAPbPb0510->GetErrorXhigh(i), graphPCMYieldPi0SysErrRAAPbPb0510->GetErrorYlow(i), graphPCMYieldPi0SysErrRAAPbPb0510->GetErrorYhigh(i));
	}
	graphPCMYieldPi0SysErrPbPb0510Red->SetPoint(0,0.5, 0);
	graphPCMYieldPi0SysErrRAAPbPb0510Red->SetPoint(0,0.5, 0);
	// 	graphPCMYieldPi0SysErrPbPb0510Red->SetPoint(16,7, 0);
	// 	graphPCMYieldPi0SysErrRAAPbPb0510Red->SetPoint(16,7, 0);

	// 	cout << "" << endl;
	// 	graphPCMYieldPi0SysErrRAAPbPb0510Red->Print();
	histoYieldPbPb0510Red = 		new TH1D("histoYieldPbPb0510Red","",18,fBinsPi0HIPtPCM);
	for (Int_t i =1; i < 19; i++){
	// 		cout << histoYieldPbPb0510Red->GetBinCenter(i) << "\t"<< histoYieldPbPb0510->GetBinCenter(i) << "\t"<< histoYieldPbPb0510->GetBinContent(i) << "\t"<< histoYieldPbPb0510->GetBinError(i) << endl; 
		histoYieldPbPb0510Red->SetBinContent(i,histoYieldPbPb0510->GetBinContent(i));
		histoYieldPbPb0510Red->SetBinError(i,histoYieldPbPb0510->GetBinError(i));
	}
	histoYieldPbPb0510Red->SetBinContent(1,0);
	histoYieldPbPb0510Red->SetBinContent(2,0);
	// 	histoYieldPbPb0510Red->SetBinContent(18,0);
	nPointsPi00510 = graphPCMYieldPi0SysErrPbPb0510Red->GetN();

	cout << "PCM: 10-20%" << endl;
	directoryPi0PbPb1020 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_10-20%"); 
	histoNumberOfEvents1020= 			(TH1D*)filePCMPbPb->Get("histoNumberOfEventsPbPb_2.76TeV10-20%");
	histoYieldPbPb1020 = 				(TH1D*)directoryPi0PbPb1020->Get(nameHistoConv.Data());
	graphPCMYieldPi0SysErrPbPb1020= 	(TGraphAsymmErrors*)directoryPi0PbPb1020->Get(nameGraphConv.Data());	
	graphPCMYieldPi0SysErrRAAPbPb1020= 	(TGraphAsymmErrors*)directoryPi0PbPb1020->Get("Pi0SystErrorA");	
	histoPCMMassData1020 = 			(TH1D*)directoryPi0PbPb1020->Get("MassPi0");
	histoPCMMassMC1020 = 			(TH1D*)directoryPi0PbPb1020->Get("TrueMassPi0");
	histoPCMWidthData1020 = 		(TH1D*)directoryPi0PbPb1020->Get("FWHMPi0MeV");
	histoPCMWidthMC1020 = 			(TH1D*)directoryPi0PbPb1020->Get("TrueFWHMPi0MeV");
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

	nEvent1020 = 							histoNumberOfEvents1020->GetBinContent(1);

	//  	graphPCMYieldPi0SysErrRAAPbPb1020->Print();
	graphPCMYieldPi0SysErrPbPb1020Red = new TGraphAsymmErrors(17);
	graphPCMYieldPi0SysErrRAAPbPb1020Red = new TGraphAsymmErrors(17);
	Double_t* xValPCMStat1020 = graphPCMYieldPi0SysErrPbPb1020->GetX();
	Double_t* yValPCMStat1020 = graphPCMYieldPi0SysErrPbPb1020->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPCMYieldPi0SysErrPbPb1020Red->SetPoint(i,xValPCMStat1020[i], yValPCMStat1020[i]);
		graphPCMYieldPi0SysErrPbPb1020Red->SetPointError(i, graphPCMYieldPi0SysErrPbPb1020->GetErrorXlow(i),graphPCMYieldPi0SysErrPbPb1020->GetErrorXhigh(i), graphPCMYieldPi0SysErrPbPb1020->GetErrorYlow(i), graphPCMYieldPi0SysErrPbPb1020->GetErrorYhigh(i));
		graphPCMYieldPi0SysErrRAAPbPb1020Red->SetPoint(i,xValPCMStat1020[i], yValPCMStat1020[i]);
		graphPCMYieldPi0SysErrRAAPbPb1020Red->SetPointError(i, graphPCMYieldPi0SysErrRAAPbPb1020->GetErrorXlow(i),graphPCMYieldPi0SysErrRAAPbPb1020->GetErrorXhigh(i), graphPCMYieldPi0SysErrRAAPbPb1020->GetErrorYlow(i), graphPCMYieldPi0SysErrRAAPbPb1020->GetErrorYhigh(i));
	}
	graphPCMYieldPi0SysErrPbPb1020Red->SetPoint(0,0.5, 0);
	graphPCMYieldPi0SysErrRAAPbPb1020Red->SetPoint(0,0.5, 0);
	// 	cout << "" << endl;
	// 	graphPCMYieldPi0SysErrRAAPbPb1020Red->Print();
	histoYieldPbPb1020Red = 		new TH1D("histoYieldPbPb1020Red","",18,fBinsPi0HIPtPCM);
	for (Int_t i =1; i < 19; i++){
	// 		cout << histoYieldPbPb1020Red->GetBinCenter(i) << "\t"<< histoYieldPbPb1020->GetBinCenter(i) << "\t"<< histoYieldPbPb1020->GetBinContent(i) << "\t"<< histoYieldPbPb1020->GetBinError(i) << endl; 
		histoYieldPbPb1020Red->SetBinContent(i,histoYieldPbPb1020->GetBinContent(i));
		histoYieldPbPb1020Red->SetBinError(i,histoYieldPbPb1020->GetBinError(i));
	}
	histoYieldPbPb1020Red->SetBinContent(1,0);
	histoYieldPbPb1020Red->SetBinContent(2,0);
	nPointsPi01020 = graphPCMYieldPi0SysErrPbPb1020Red->GetN();

	cout << "PCM: 20-40%" << endl;
	directoryPi0PbPb2040 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_20-40%"); 
	histoNumberOfEvents2040= 			(TH1D*)filePCMPbPb->Get("histoNumberOfEventsPbPb_2.76TeV20-40%");
	histoYieldPbPb2040 = 				(TH1D*)directoryPi0PbPb2040->Get(nameHistoConv.Data());
	graphPCMYieldPi0SysErrPbPb2040= 	(TGraphAsymmErrors*)directoryPi0PbPb2040->Get(nameGraphConv.Data());	
	graphPCMYieldPi0SysErrRAAPbPb2040= 	(TGraphAsymmErrors*)directoryPi0PbPb2040->Get("Pi0SystErrorA");		
	nEvent2040 = 							histoNumberOfEvents2040->GetBinContent(1);

	// 	graphPCMYieldPi0SysErrPbPb2040->Print();
	graphPCMYieldPi0SysErrRAAPbPb2040Red = new TGraphAsymmErrors(17);
	graphPCMYieldPi0SysErrPbPb2040Red = new TGraphAsymmErrors(17);
	Double_t* xValPCMStat2040 = graphPCMYieldPi0SysErrPbPb2040->GetX();
	Double_t* yValPCMStat2040 = graphPCMYieldPi0SysErrPbPb2040->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPCMYieldPi0SysErrPbPb2040Red->SetPoint(i,xValPCMStat2040[i], yValPCMStat2040[i]);
		graphPCMYieldPi0SysErrPbPb2040Red->SetPointError(i, graphPCMYieldPi0SysErrPbPb2040->GetErrorXlow(i),graphPCMYieldPi0SysErrPbPb2040->GetErrorXhigh(i), graphPCMYieldPi0SysErrPbPb2040->GetErrorYlow(i), graphPCMYieldPi0SysErrPbPb2040->GetErrorYhigh(i));
		graphPCMYieldPi0SysErrRAAPbPb2040Red->SetPoint(i,xValPCMStat2040[i], yValPCMStat2040[i]);
		graphPCMYieldPi0SysErrRAAPbPb2040Red->SetPointError(i, graphPCMYieldPi0SysErrRAAPbPb2040->GetErrorXlow(i),graphPCMYieldPi0SysErrRAAPbPb2040->GetErrorXhigh(i), graphPCMYieldPi0SysErrRAAPbPb2040->GetErrorYlow(i), graphPCMYieldPi0SysErrRAAPbPb2040->GetErrorYhigh(i));
	}
	graphPCMYieldPi0SysErrPbPb2040Red->SetPoint(0,0.5, 0);
	graphPCMYieldPi0SysErrRAAPbPb2040Red->SetPoint(0,0.5, 0);
	// 	cout << "" << endl;
	// 	graphPCMYieldPi0SysErrPbPb2040Red->Print();
	histoYieldPbPb2040Red = 		new TH1D("histoYieldPbPb2040Red","",18,fBinsPi0HIPtPCM);
	for (Int_t i =1; i < 19; i++){
	// 		cout << histoYieldPbPb2040Red->GetBinCenter(i) << "\t"<< histoYieldPbPb2040->GetBinCenter(i) << "\t"<< histoYieldPbPb2040->GetBinContent(i) << "\t"<< histoYieldPbPb2040->GetBinError(i) << endl; 
		histoYieldPbPb2040Red->SetBinContent(i,histoYieldPbPb2040->GetBinContent(i));
		histoYieldPbPb2040Red->SetBinError(i,histoYieldPbPb2040->GetBinError(i));
	}
	histoYieldPbPb2040Red->SetBinContent(1,0);
	histoYieldPbPb2040Red->SetBinContent(2,0);
	nPointsPi02040 = graphPCMYieldPi0SysErrPbPb2040Red->GetN();


	cout << "PCM: 40-60%" << endl;
	directoryPi0PbPb4060 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_40-60%"); 
	histoNumberOfEvents4060= 			(TH1D*)filePCMPbPb->Get("histoNumberOfEventsPbPb_2.76TeV40-60%");
	histoYieldPbPb4060 = 				(TH1D*)directoryPi0PbPb4060->Get(nameHistoConv.Data());
	graphPCMYieldPi0SysErrPbPb4060= 	(TGraphAsymmErrors*)directoryPi0PbPb4060->Get(nameGraphConv.Data());	
	graphPCMYieldPi0SysErrRAAPbPb4060= 	(TGraphAsymmErrors*)directoryPi0PbPb4060->Get("Pi0SystErrorA");	
	nEvent4060 = 							histoNumberOfEvents4060->GetBinContent(1);

	// 	graphPCMYieldPi0SysErrPbPb4060->Print();
	graphPCMYieldPi0SysErrRAAPbPb4060Red = new TGraphAsymmErrors(17);
	graphPCMYieldPi0SysErrPbPb4060Red = new TGraphAsymmErrors(17);
	Double_t* xValPCMStat4060 = graphPCMYieldPi0SysErrPbPb4060->GetX();
	Double_t* yValPCMStat4060 = graphPCMYieldPi0SysErrPbPb4060->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPCMYieldPi0SysErrPbPb4060Red->SetPoint(i,xValPCMStat4060[i], yValPCMStat4060[i]);
		graphPCMYieldPi0SysErrPbPb4060Red->SetPointError(i, graphPCMYieldPi0SysErrPbPb4060->GetErrorXlow(i),graphPCMYieldPi0SysErrPbPb4060->GetErrorXhigh(i), graphPCMYieldPi0SysErrPbPb4060->GetErrorYlow(i), graphPCMYieldPi0SysErrPbPb4060->GetErrorYhigh(i));
		graphPCMYieldPi0SysErrRAAPbPb4060Red->SetPoint(i,xValPCMStat4060[i], yValPCMStat4060[i]);
		graphPCMYieldPi0SysErrRAAPbPb4060Red->SetPointError(i, graphPCMYieldPi0SysErrRAAPbPb4060->GetErrorXlow(i),graphPCMYieldPi0SysErrRAAPbPb4060->GetErrorXhigh(i), graphPCMYieldPi0SysErrRAAPbPb4060->GetErrorYlow(i), graphPCMYieldPi0SysErrRAAPbPb4060->GetErrorYhigh(i));
	}
	graphPCMYieldPi0SysErrPbPb4060Red->SetPoint(0,0.5, 0);
	graphPCMYieldPi0SysErrRAAPbPb4060Red->SetPoint(0,0.5, 0);

	// 	cout << "" << endl;
	// 	graphPCMYieldPi0SysErrPbPb4060Red->Print();
	histoYieldPbPb4060Red = 		new TH1D("histoYieldPbPb4060Red","",18,fBinsPi0HIPtPCM);
	for (Int_t i =1; i < 19; i++){
	// 		cout << histoYieldPbPb4060Red->GetBinCenter(i) << "\t"<< histoYieldPbPb4060->GetBinCenter(i) << "\t"<< histoYieldPbPb4060->GetBinContent(i) << "\t"<< histoYieldPbPb4060->GetBinError(i) << endl; 
		histoYieldPbPb4060Red->SetBinContent(i,histoYieldPbPb4060->GetBinContent(i));
		histoYieldPbPb4060Red->SetBinError(i,histoYieldPbPb4060->GetBinError(i));
	}
	histoYieldPbPb4060Red->SetBinContent(1,0);
	histoYieldPbPb4060Red->SetBinContent(2,0);	
	nPointsPi04060 = graphPCMYieldPi0SysErrPbPb4060Red->GetN();

	cout << "PCM: 60-80%" << endl;
	directoryPi0PbPb6080 = 				(TDirectory*)filePCMPbPb->Get("Pi0_PbPb_2.76TeV_60-80%"); 
	histoNumberOfEvents6080= 			(TH1D*)filePCMPbPb->Get("histoNumberOfEventsPbPb_2.76TeV60-80%");
	histoYieldPbPb6080 = 				(TH1D*)directoryPi0PbPb6080->Get(nameHistoConv.Data());
	graphPCMYieldPi0SysErrPbPb6080= 	(TGraphAsymmErrors*)directoryPi0PbPb6080->Get(nameGraphConv.Data());	
	graphPCMYieldPi0SysErrRAAPbPb6080= 	(TGraphAsymmErrors*)directoryPi0PbPb6080->Get("Pi0SystErrorA");	
	histoPCMMassData6080 = 			(TH1D*)directoryPi0PbPb6080->Get("MassPi0");
	histoPCMMassMC6080 = 			(TH1D*)directoryPi0PbPb6080->Get("TrueMassPi0");
	histoPCMWidthData6080 = 		(TH1D*)directoryPi0PbPb6080->Get("FWHMPi0MeV");
	histoPCMWidthMC6080 = 			(TH1D*)directoryPi0PbPb6080->Get("TrueFWHMPi0MeV");
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
	for (Int_t  i = 0; i < histoPCMMassData6080->GetXaxis()->FindBin(0.6); i++){
		histoPCMMassData6080->SetBinContent(i,0);
		histoPCMMassMC6080->SetBinContent(i,0);
		histoPCMWidthData6080->SetBinContent(i,10000.);  
		histoPCMWidthMC6080->SetBinContent(i,10000.);  
	}



	nEvent6080 = 							histoNumberOfEvents6080->GetBinContent(1);

	// 	graphPCMYieldPi0SysErrPbPb6080->Print();
	graphPCMYieldPi0SysErrRAAPbPb6080Red = new TGraphAsymmErrors(17);
	graphPCMYieldPi0SysErrPbPb6080Red = new TGraphAsymmErrors(17);
	Double_t* xValPCMStat6080 = graphPCMYieldPi0SysErrPbPb6080->GetX();
	Double_t* yValPCMStat6080 = graphPCMYieldPi0SysErrPbPb6080->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPCMYieldPi0SysErrPbPb6080Red->SetPoint(i,xValPCMStat6080[i], yValPCMStat6080[i]);
		graphPCMYieldPi0SysErrPbPb6080Red->SetPointError(i, graphPCMYieldPi0SysErrPbPb6080->GetErrorXlow(i),graphPCMYieldPi0SysErrPbPb6080->GetErrorXhigh(i), graphPCMYieldPi0SysErrPbPb6080->GetErrorYlow(i), graphPCMYieldPi0SysErrPbPb6080->GetErrorYhigh(i));
		graphPCMYieldPi0SysErrRAAPbPb6080Red->SetPoint(i,xValPCMStat6080[i], yValPCMStat6080[i]);
		graphPCMYieldPi0SysErrRAAPbPb6080Red->SetPointError(i, graphPCMYieldPi0SysErrRAAPbPb6080->GetErrorXlow(i),graphPCMYieldPi0SysErrRAAPbPb6080->GetErrorXhigh(i), graphPCMYieldPi0SysErrRAAPbPb6080->GetErrorYlow(i), graphPCMYieldPi0SysErrRAAPbPb6080->GetErrorYhigh(i));
	}
	graphPCMYieldPi0SysErrPbPb6080Red->SetPoint(0,0.5, 0);
	graphPCMYieldPi0SysErrRAAPbPb6080Red->SetPoint(0,0.5, 0);

	// 	cout << "" << endl;
	// 	graphPCMYieldPi0SysErrPbPb6080Red->Print();
	histoYieldPbPb6080Red = 		new TH1D("histoYieldPbPb6080Red","",18,fBinsPi0HIPtPCM);
	for (Int_t i =1; i < 19; i++){
	// 		cout << histoYieldPbPb6080Red->GetBinCenter(i) << "\t"<< histoYieldPbPb6080->GetBinCenter(i) << "\t"<< histoYieldPbPb6080->GetBinContent(i) << "\t"<< histoYieldPbPb6080->GetBinError(i) << endl; 
		histoYieldPbPb6080Red->SetBinContent(i,histoYieldPbPb6080->GetBinContent(i));
		histoYieldPbPb6080Red->SetBinError(i,histoYieldPbPb6080->GetBinError(i));
	}
	histoYieldPbPb6080Red->SetBinContent(1,0);
	histoYieldPbPb6080Red->SetBinContent(2,0);   

	nPointsPi06080 = graphPCMYieldPi0SysErrPbPb6080Red->GetN();

	//************************** Read data for PHOS **************************************************
	filePHOSPbPb = 		new TFile(fileNameCaloPHOSPbPb);
	directoryPHOSPi0PbPb0005 = (TDirectory*)filePHOSPbPb->Get("pi0_PbPb_2760_Centrality_0-5%");
	directoryPHOSPi0PbPb0510 = (TDirectory*)filePHOSPbPb->Get("pi0_PbPb_2760_Centrality_5-10%");
	directoryPHOSPi0PbPb0010 = (TDirectory*)filePHOSPbPb->Get("pi0_PbPb_2760_Centrality_0-10%");
	directoryPHOSPi0PbPb1020 = (TDirectory*)filePHOSPbPb->Get("pi0_PbPb_2760_Centrality_10-20%");
	directoryPHOSPi0PbPb2040 = (TDirectory*)filePHOSPbPb->Get("pi0_PbPb_2760_Centrality_20-40%");
	directoryPHOSPi0PbPb4060 = (TDirectory*)filePHOSPbPb->Get("pi0_PbPb_2760_Centrality_40-60%");
	directoryPHOSPi0PbPb6080 = (TDirectory*)filePHOSPbPb->Get("pi0_PbPb_2760_Centrality_60-80%");

	cout << "PHOS: 0-5%" << endl;
	histoPi0PHOSPbPb0005 = 		(TH1D*)directoryPHOSPi0PbPb0005->Get("hPi0_PbPb_cen0_NoBW_Stat");
	histoPi0PHOSSysPbPb0005 = 	(TH1D*)directoryPHOSPi0PbPb0005->Get("hPi0_PbPb_cen0_NoBW_Syst");
	graphPHOSYieldPi0SysErrPbPb0005 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb0005);	
	histoPi0PHOSSysRAAPbPb0005 = 	(TH1D*)directoryPHOSPi0PbPb0005->Get("hPi0_PbPb_cen0_SystRaa");
	graphSysErrRAAYieldPi0PHOSPbPb0005 = new TGraphAsymmErrors(histoPi0PHOSSysRAAPbPb0005);	
	histoPHOSMassData0005 =       (TH1D*)directoryPHOSPi0PbPb0005->Get("mass1_GS_Both2core_cen0");
	histoPHOSMassMC0005 =         (TH1D*)directoryPHOSPi0PbPb0005->Get("MC_mass1_GS_Both2core_cen0");
	histoPHOSWidthData0005 =      (TH1D*)directoryPHOSPi0PbPb0005->Get("width1_GS_Both2core_cen0");
	histoPHOSWidthMC0005 =        (TH1D*)directoryPHOSPi0PbPb0005->Get("MC_width1_GS_Both2core_cen0");
	histoPHOSWidthData0005->Scale(1000.);
	histoPHOSWidthMC0005->Scale(1000.);
	histoPHOSMassData0005->Scale(1000.);
	histoPHOSMassMC0005->Scale(1000.);
	for (Int_t  i = 0; i < histoPHOSWidthData0005->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthData0005->SetBinContent(i,10000.);
		histoPHOSMassData0005->SetBinContent(i,10000.);
	}
	for (Int_t  i = histoPHOSWidthData0005->GetXaxis()->FindBin(12.1); i < histoPHOSWidthData0005->GetNbinsX()+1; i++){
		histoPHOSWidthData0005->SetBinContent(i,10000.);
		histoPHOSMassData0005->SetBinContent(i,10000.);
	}
	for (Int_t  i = 0; i < histoPHOSWidthMC0005->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthMC0005->SetBinContent(i,10000.);
		histoPHOSMassMC0005->SetBinContent(i,10000.);
	}
	for (Int_t  i = histoPHOSWidthMC0005->GetXaxis()->FindBin(12.1); i < histoPHOSWidthMC0005->GetNbinsX()+1; i++){
		histoPHOSWidthMC0005->SetBinContent(i,10000.);
		histoPHOSMassMC0005->SetBinContent(i,10000.);
	}

	// 	graphPHOSYieldPi0SysErrPbPb0005->Print();
	graphPHOSYieldPi0SysErrPbPb0005Red = new TGraphAsymmErrors(17);
	graphSysErrRAAYieldPi0PHOSPbPb0005Red = new TGraphAsymmErrors(17);
	Double_t* xValPHOStat0005 = graphPHOSYieldPi0SysErrPbPb0005->GetX();
	Double_t* yValPHOSStat0005 = graphPHOSYieldPi0SysErrPbPb0005->GetY();
	Int_t phosOffset = 1;
	for (Int_t i = 0; i < 17; i++){
		graphPHOSYieldPi0SysErrPbPb0005Red->SetPoint(i,xValPHOStat0005[i+phosOffset], yValPHOSStat0005[i+phosOffset]);
		graphPHOSYieldPi0SysErrPbPb0005Red->SetPointError(i, graphPHOSYieldPi0SysErrPbPb0005->GetErrorXlow(i+phosOffset),graphPHOSYieldPi0SysErrPbPb0005->GetErrorXhigh(i+phosOffset), graphPHOSYieldPi0SysErrPbPb0005->GetErrorYlow(i+phosOffset), graphPHOSYieldPi0SysErrPbPb0005->GetErrorYhigh(i+phosOffset));
		graphSysErrRAAYieldPi0PHOSPbPb0005Red->SetPoint(i,xValPHOStat0005[i+phosOffset], yValPHOSStat0005[i+phosOffset]);
		graphSysErrRAAYieldPi0PHOSPbPb0005Red->SetPointError(i, graphSysErrRAAYieldPi0PHOSPbPb0005->GetErrorXlow(i+phosOffset),graphSysErrRAAYieldPi0PHOSPbPb0005->GetErrorXhigh(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb0005->GetErrorYlow(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb0005->GetErrorYhigh(i+phosOffset));
	}
	cout << "" << endl;
	graphPHOSYieldPi0SysErrPbPb0005Red->SetPoint(0,0.9, 0);
	graphSysErrRAAYieldPi0PHOSPbPb0005Red->SetPoint(0,0.9, 0);
	//  	graphPHOSYieldPi0SysErrPbPb0005Red->Print();

	histoPi0PHOSPbPb0005Red = 		new TH1D("histoPi0PHOSPbPb0005Red","",18,fBinsPi0HIPtPHOS);
	for (Int_t i =1; i < 18; i++){
		cout << histoPi0PHOSPbPb0005Red->GetBinCenter(i) << "\t"<< histoPi0PHOSPbPb0005->GetBinCenter(i+phosOffset) << "\t"<< histoPi0PHOSPbPb0005->GetBinContent(i+phosOffset) << "\t"<< histoPi0PHOSPbPb0005->GetBinError(i+phosOffset) << endl; 
		histoPi0PHOSPbPb0005Red->SetBinContent(i,histoPi0PHOSPbPb0005->GetBinContent(i+phosOffset));
		histoPi0PHOSPbPb0005Red->SetBinError(i,histoPi0PHOSPbPb0005->GetBinError(i+phosOffset));
	}
	histoPi0PHOSPbPb0005Red->SetBinContent(1,0);

	cout << "PHOS: 0-10%" << endl;
	histoPi0PHOSPbPb0010 = 		(TH1D*)directoryPHOSPi0PbPb0010->Get("hPi0_PbPb_cen7_NoBW_Stat");
	histoPi0PHOSSysPbPb0010 = 	(TH1D*)directoryPHOSPi0PbPb0010->Get("hPi0_PbPb_cen7_NoBW_Syst");
	graphPHOSYieldPi0SysErrPbPb0010 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb0010);	
	histoPi0PHOSSysRAAPbPb0010 = 	(TH1D*)directoryPHOSPi0PbPb0010->Get("hPi0_PbPb_cen7_SystRaa");
	graphSysErrRAAYieldPi0PHOSPbPb0010 = new TGraphAsymmErrors(histoPi0PHOSSysRAAPbPb0010);	
	histoPHOSMassData0010 =       (TH1D*)directoryPHOSPi0PbPb0010->Get("mass1_GS_Both2core_cen7");
	histoPHOSMassMC0010 =         (TH1D*)directoryPHOSPi0PbPb0010->Get("MC_mass1_GS_Both2core_cen7");
	histoPHOSWidthData0010 =      (TH1D*)directoryPHOSPi0PbPb0010->Get("width1_GS_Both2core_cen7");
	histoPHOSWidthMC0010 =        (TH1D*)directoryPHOSPi0PbPb0010->Get("MC_width1_GS_Both2core_cen7");
	histoPHOSWidthData0010->Scale(1000.);
	histoPHOSWidthMC0010->Scale(1000.);
	histoPHOSMassData0010->Scale(1000.);
	histoPHOSMassMC0010->Scale(1000.);
	for (Int_t  i = 0; i < histoPHOSWidthData0010->GetXaxis()->FindBin(1.0); i++){
		histoPHOSMassData0010->SetBinContent(i,0.);
		histoPHOSWidthData0010->SetBinContent(i,10000.);
	}   
	for (Int_t  i = 0; i < histoPHOSWidthMC0010->GetXaxis()->FindBin(1.0); i++){  
		
		histoPHOSMassMC0010->SetBinContent(i,0.);
		histoPHOSWidthMC0010->SetBinContent(i,10000.);
		
	}
	cout << histoPHOSWidthData0010->GetXaxis()->FindBin(15.0) << endl;
	for (Int_t  i = histoPHOSWidthData0010->GetXaxis()->FindBin(15.0); i < histoPHOSWidthData0010->GetNbinsX()+1; i++){
		histoPHOSMassData0010->SetBinContent(i,0);
		histoPHOSWidthData0010->SetBinContent(i,10000.);
		
	}
	for (Int_t  i = histoPHOSWidthMC0010->GetXaxis()->FindBin(15.0); i < histoPHOSWidthMC0010->GetNbinsX()+1; i++){   
		histoPHOSMassMC0010->SetBinContent(i,0);;
		histoPHOSWidthMC0010->SetBinContent(i,10000.);   
	}      

	// 	graphPHOSYieldPi0SysErrPbPb0010->Print();
	graphPHOSYieldPi0SysErrPbPb0010Red = new TGraphAsymmErrors(17);
	graphSysErrRAAYieldPi0PHOSPbPb0010Red = new TGraphAsymmErrors(17);
	Double_t* xValPHOStat0010 = graphPHOSYieldPi0SysErrPbPb0010->GetX();
	Double_t* yValPHOSStat0010 = graphPHOSYieldPi0SysErrPbPb0010->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPHOSYieldPi0SysErrPbPb0010Red->SetPoint(i,xValPHOStat0010[i+phosOffset], yValPHOSStat0010[i+phosOffset]);
		graphPHOSYieldPi0SysErrPbPb0010Red->SetPointError(i, graphPHOSYieldPi0SysErrPbPb0010->GetErrorXlow(i+phosOffset),graphPHOSYieldPi0SysErrPbPb0010->GetErrorXhigh(i+phosOffset), graphPHOSYieldPi0SysErrPbPb0010->GetErrorYlow(i+phosOffset), graphPHOSYieldPi0SysErrPbPb0010->GetErrorYhigh(i+phosOffset));
		graphSysErrRAAYieldPi0PHOSPbPb0010Red->SetPoint(i,xValPHOStat0010[i+phosOffset], yValPHOSStat0010[i+phosOffset]);
		graphSysErrRAAYieldPi0PHOSPbPb0010Red->SetPointError(i, graphSysErrRAAYieldPi0PHOSPbPb0010->GetErrorXlow(i+phosOffset),graphSysErrRAAYieldPi0PHOSPbPb0010->GetErrorXhigh(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb0010->GetErrorYlow(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb0010->GetErrorYhigh(i+phosOffset));
	}
	cout << "" << endl;
	graphPHOSYieldPi0SysErrPbPb0010Red->SetPoint(0,0.9, 0);
	graphSysErrRAAYieldPi0PHOSPbPb0010Red->SetPoint(0,0.9, 0);
	//  	graphPHOSYieldPi0SysErrPbPb0010Red->Print();

	histoPi0PHOSPbPb0010Red = 		new TH1D("histoPi0PHOSPbPb0010Red","",18,fBinsPi0HIPtPHOS);
	for (Int_t i =1; i < 18; i++){
		cout << histoPi0PHOSPbPb0010Red->GetBinCenter(i) << "\t"<< histoPi0PHOSPbPb0010->GetBinCenter(i+phosOffset) << "\t"<< histoPi0PHOSPbPb0010->GetBinContent(i+phosOffset) << "\t"<< histoPi0PHOSPbPb0010->GetBinError(i+phosOffset) << endl; 
		histoPi0PHOSPbPb0010Red->SetBinContent(i,histoPi0PHOSPbPb0010->GetBinContent(i+phosOffset));
		histoPi0PHOSPbPb0010Red->SetBinError(i,histoPi0PHOSPbPb0010->GetBinError(i+phosOffset));
	}
	histoPi0PHOSPbPb0010Red->SetBinContent(1,0);

	cout << "PHOS: 5-10%" << endl;
	histoPi0PHOSPbPb0510 = 		(TH1D*)directoryPHOSPi0PbPb0510->Get("hPi0_PbPb_cen1_NoBW_Stat");
	histoPi0PHOSSysPbPb0510 = 	(TH1D*)directoryPHOSPi0PbPb0510->Get("hPi0_PbPb_cen1_NoBW_Syst");
	graphPHOSYieldPi0SysErrPbPb0510 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb0510);	
	histoPi0PHOSSysRAAPbPb0510 = 	(TH1D*)directoryPHOSPi0PbPb0510->Get("hPi0_PbPb_cen1_SystRaa");
	graphSysErrRAAYieldPi0PHOSPbPb0510 = new TGraphAsymmErrors(histoPi0PHOSSysRAAPbPb0510);	
	histoPHOSMassData0510 =       (TH1D*)directoryPHOSPi0PbPb0510->Get("mass1_GS_Both2core_cen1");
	histoPHOSMassMC0510 =         (TH1D*)directoryPHOSPi0PbPb0510->Get("MC_mass1_GS_Both2core_cen1");
	histoPHOSWidthData0510 =      (TH1D*)directoryPHOSPi0PbPb0510->Get("width1_GS_Both2core_cen1");
	histoPHOSWidthMC0510 =        (TH1D*)directoryPHOSPi0PbPb0510->Get("MC_width1_GS_Both2core_cen1");
	histoPHOSWidthData0510->Scale(1000.);
	histoPHOSWidthMC0510->Scale(1000.);
	histoPHOSMassData0510->Scale(1000.);
	histoPHOSMassMC0510->Scale(1000.);
	for (Int_t  i = 0; i < histoPHOSWidthData0510->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthData0510->SetBinContent(i,10000.);
		histoPHOSMassData0510->SetBinContent(i,10000.);
	}
	for (Int_t  i = 0; i < histoPHOSWidthMC0510->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthMC0510->SetBinContent(i,10000.);
		histoPHOSMassMC0510->SetBinContent(i,10000.);
	}
	histoPHOSMassData0510->SetBinContent(histoPHOSMassData0510->GetNbinsX(),0);
	histoPHOSMassMC0510->SetBinContent(histoPHOSMassMC0510->GetNbinsX(),0);;
	histoPHOSWidthData0510->SetBinContent(histoPHOSWidthData0510->GetNbinsX(),10000.);
	histoPHOSWidthMC0510->SetBinContent(histoPHOSWidthMC0510->GetNbinsX(),10000.);

	// 	graphPHOSYieldPi0SysErrPbPb0510->Print();
	graphPHOSYieldPi0SysErrPbPb0510Red = new TGraphAsymmErrors(17);
	graphSysErrRAAYieldPi0PHOSPbPb0510Red = new TGraphAsymmErrors(17);
	Double_t* xValPHOStat0510 = graphPHOSYieldPi0SysErrPbPb0510->GetX();
	Double_t* yValPHOSStat0510 = graphPHOSYieldPi0SysErrPbPb0510->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPHOSYieldPi0SysErrPbPb0510Red->SetPoint(i,xValPHOStat0510[i+phosOffset], yValPHOSStat0510[i+phosOffset]);
		graphPHOSYieldPi0SysErrPbPb0510Red->SetPointError(i, graphPHOSYieldPi0SysErrPbPb0510->GetErrorXlow(i+phosOffset),graphPHOSYieldPi0SysErrPbPb0510->GetErrorXhigh(i+phosOffset), graphPHOSYieldPi0SysErrPbPb0510->GetErrorYlow(i+phosOffset), graphPHOSYieldPi0SysErrPbPb0510->GetErrorYhigh(i+phosOffset));
		graphSysErrRAAYieldPi0PHOSPbPb0510Red->SetPoint(i,xValPHOStat0510[i+phosOffset], yValPHOSStat0510[i+phosOffset]);
		graphSysErrRAAYieldPi0PHOSPbPb0510Red->SetPointError(i, graphSysErrRAAYieldPi0PHOSPbPb0510->GetErrorXlow(i+phosOffset),graphSysErrRAAYieldPi0PHOSPbPb0510->GetErrorXhigh(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb0510->GetErrorYlow(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb0510->GetErrorYhigh(i+phosOffset));
	}
	graphPHOSYieldPi0SysErrPbPb0510Red->SetPoint(0,0.9, 0);
	graphSysErrRAAYieldPi0PHOSPbPb0510Red->SetPoint(0,0.9, 0);

	histoPi0PHOSPbPb0510Red = 		new TH1D("histoPi0PHOSPbPb0510Red","",18,fBinsPi0HIPtPHOS);
	for (Int_t i =1; i < 18; i++){
	//  		cout << histoPi0PHOSPbPb0510Red->GetBinCenter(i) << "\t"<< histoPi0PHOSPbPb0510->GetBinCenter(i+phosOffset) << "\t"<< histoPi0PHOSPbPb0510->GetBinContent(i+phosOffset) << "\t"<< histoPi0PHOSPbPb0510->GetBinError(i+phosOffset) << endl; 
		histoPi0PHOSPbPb0510Red->SetBinContent(i,histoPi0PHOSPbPb0510->GetBinContent(i+phosOffset));
		histoPi0PHOSPbPb0510Red->SetBinError(i,histoPi0PHOSPbPb0510->GetBinError(i+phosOffset));
	}
	histoPi0PHOSPbPb0510Red->SetBinContent(1,0);

	cout << "PHOS: 10-20%" << endl;
	histoPi0PHOSPbPb1020 = 		(TH1D*)directoryPHOSPi0PbPb1020->Get("hPi0_PbPb_cen2_NoBW_Stat");
	histoPi0PHOSSysPbPb1020 = 	(TH1D*)directoryPHOSPi0PbPb1020->Get("hPi0_PbPb_cen2_NoBW_Syst");
	graphPHOSYieldPi0SysErrPbPb1020 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb1020);	
	histoPi0PHOSSysRAAPbPb1020 = 	(TH1D*)directoryPHOSPi0PbPb1020->Get("hPi0_PbPb_cen2_SystRaa");
	graphSysErrRAAYieldPi0PHOSPbPb1020 = new TGraphAsymmErrors(histoPi0PHOSSysRAAPbPb1020);	
	histoPHOSMassData1020 =       (TH1D*)directoryPHOSPi0PbPb1020->Get("mass1_GS_Both2core_cen2");
	histoPHOSMassMC1020 =         (TH1D*)directoryPHOSPi0PbPb1020->Get("MC_mass1_GS_Both2core_cen2");
	histoPHOSWidthData1020 =      (TH1D*)directoryPHOSPi0PbPb1020->Get("width1_GS_Both2core_cen2");
	histoPHOSWidthMC1020 =        (TH1D*)directoryPHOSPi0PbPb1020->Get("MC_width1_GS_Both2core_cen2");
	histoPHOSWidthData1020->Scale(1000.);
	histoPHOSWidthMC1020->Scale(1000.);
	histoPHOSMassData1020->Scale(1000.);
	histoPHOSMassMC1020->Scale(1000.);
	for (Int_t  i = 0; i < histoPHOSWidthData1020->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthData1020->SetBinContent(i,10000.);
		histoPHOSMassData1020->SetBinContent(i,10000.);
	}
	for (Int_t  i = 0; i < histoPHOSWidthMC1020->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthMC1020->SetBinContent(i,10000.);
		histoPHOSMassMC1020->SetBinContent(i,10000.);
	}
	histoPHOSMassData1020->SetBinContent(histoPHOSMassData1020->GetNbinsX(),0);
	histoPHOSMassMC1020->SetBinContent(histoPHOSMassMC1020->GetNbinsX(),0);;
	histoPHOSWidthData1020->SetBinContent(histoPHOSWidthData1020->GetNbinsX(),10000.);
	histoPHOSWidthMC1020->SetBinContent(histoPHOSWidthMC1020->GetNbinsX(),10000.);

	// 	graphPHOSYieldPi0SysErrPbPb1020->Print();
	graphPHOSYieldPi0SysErrPbPb1020Red = new TGraphAsymmErrors(17);
	graphSysErrRAAYieldPi0PHOSPbPb1020Red = new TGraphAsymmErrors(17);
	Double_t* xValPHOStat1020 = graphPHOSYieldPi0SysErrPbPb1020->GetX();
	Double_t* yValPHOSStat1020 = graphPHOSYieldPi0SysErrPbPb1020->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPHOSYieldPi0SysErrPbPb1020Red->SetPoint(i,xValPHOStat1020[i+phosOffset], yValPHOSStat1020[i+phosOffset]);
		graphPHOSYieldPi0SysErrPbPb1020Red->SetPointError(i, graphPHOSYieldPi0SysErrPbPb1020->GetErrorXlow(i+phosOffset),graphPHOSYieldPi0SysErrPbPb1020->GetErrorXhigh(i+phosOffset), graphPHOSYieldPi0SysErrPbPb1020->GetErrorYlow(i+phosOffset), graphPHOSYieldPi0SysErrPbPb1020->GetErrorYhigh(i+phosOffset));
		graphSysErrRAAYieldPi0PHOSPbPb1020Red->SetPoint(i,xValPHOStat1020[i+phosOffset], yValPHOSStat1020[i+phosOffset]);
		graphSysErrRAAYieldPi0PHOSPbPb1020Red->SetPointError(i, graphSysErrRAAYieldPi0PHOSPbPb1020->GetErrorXlow(i+phosOffset),graphSysErrRAAYieldPi0PHOSPbPb1020->GetErrorXhigh(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb1020->GetErrorYlow(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb1020->GetErrorYhigh(i+phosOffset));
	}
	graphPHOSYieldPi0SysErrPbPb1020Red->SetPoint(0,0.9, 0);
	graphSysErrRAAYieldPi0PHOSPbPb1020Red->SetPoint(0,0.9, 0);

	histoPi0PHOSPbPb1020Red = 		new TH1D("histoPi0PHOSPbPb1020Red","",18,fBinsPi0HIPtPHOS);
	for (Int_t i =1; i < 18; i++){
	//  		cout << histoPi0PHOSPbPb1020Red->GetBinCenter(i) << "\t"<< histoPi0PHOSPbPb1020->GetBinCenter(i+phosOffset) << "\t"<< histoPi0PHOSPbPb1020->GetBinContent(i+phosOffset) << "\t"<< histoPi0PHOSPbPb1020->GetBinError(i+phosOffset) << endl; 
		histoPi0PHOSPbPb1020Red->SetBinContent(i,histoPi0PHOSPbPb1020->GetBinContent(i+phosOffset));
		histoPi0PHOSPbPb1020Red->SetBinError(i,histoPi0PHOSPbPb1020->GetBinError(i+phosOffset));
	}
	histoPi0PHOSPbPb1020Red->SetBinContent(1,0);

	cout << "PHOS: 20-40%" << endl;
	histoPi0PHOSPbPb2040 = 		(TH1D*)directoryPHOSPi0PbPb2040->Get("hPi0_PbPb_cen3_NoBW_Stat");
	histoPi0PHOSSysPbPb2040 = 	(TH1D*)directoryPHOSPi0PbPb2040->Get("hPi0_PbPb_cen3_NoBW_Syst");
	graphPHOSYieldPi0SysErrPbPb2040 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb2040);	
	histoPi0PHOSSysRAAPbPb2040 = 	(TH1D*)directoryPHOSPi0PbPb2040->Get("hPi0_PbPb_cen3_SystRaa");
	graphSysErrRAAYieldPi0PHOSPbPb2040 = new TGraphAsymmErrors(histoPi0PHOSSysRAAPbPb2040);	
	histoPHOSMassData2040 =       (TH1D*)directoryPHOSPi0PbPb2040->Get("mass1_GS_Both2core_cen3");
	histoPHOSMassMC2040 =         (TH1D*)directoryPHOSPi0PbPb2040->Get("MC_mass1_GS_Both2core_cen3");
	histoPHOSWidthData2040 =      (TH1D*)directoryPHOSPi0PbPb2040->Get("width1_GS_Both2core_cen3");
	histoPHOSWidthMC2040 =        (TH1D*)directoryPHOSPi0PbPb2040->Get("MC_width1_GS_Both2core_cen3");
	histoPHOSWidthData2040->Scale(1000.);
	histoPHOSWidthMC2040->Scale(1000.);
	histoPHOSMassData2040->Scale(1000.);
	histoPHOSMassMC2040->Scale(1000.);
	for (Int_t  i = 0; i < histoPHOSWidthData2040->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthData2040->SetBinContent(i,10000.);
		histoPHOSMassData2040->SetBinContent(i,10000.);
	}
	for (Int_t  i = 0; i < histoPHOSWidthMC2040->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthMC2040->SetBinContent(i,10000.);
		histoPHOSMassMC2040->SetBinContent(i,10000.);
	}
	histoPHOSMassData2040->SetBinContent(histoPHOSMassData2040->GetNbinsX(),0);
	histoPHOSMassMC2040->SetBinContent(histoPHOSMassMC2040->GetNbinsX(),0);;
	histoPHOSWidthData2040->SetBinContent(histoPHOSWidthData2040->GetNbinsX(),10000.);
	histoPHOSWidthMC2040->SetBinContent(histoPHOSWidthMC2040->GetNbinsX(),10000.);

	graphPHOSYieldPi0SysErrPbPb2040Red = new TGraphAsymmErrors(17);
	graphSysErrRAAYieldPi0PHOSPbPb2040Red = new TGraphAsymmErrors(17);
	Double_t* xValPHOStat2040 = graphPHOSYieldPi0SysErrPbPb2040->GetX();
	Double_t* yValPHOSStat2040 = graphPHOSYieldPi0SysErrPbPb2040->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPHOSYieldPi0SysErrPbPb2040Red->SetPoint(i,xValPHOStat2040[i+phosOffset], yValPHOSStat2040[i+phosOffset]);
		graphPHOSYieldPi0SysErrPbPb2040Red->SetPointError(i, graphPHOSYieldPi0SysErrPbPb2040->GetErrorXlow(i+phosOffset),graphPHOSYieldPi0SysErrPbPb2040->GetErrorXhigh(i+phosOffset), graphPHOSYieldPi0SysErrPbPb2040->GetErrorYlow(i+phosOffset), graphPHOSYieldPi0SysErrPbPb2040->GetErrorYhigh(i+phosOffset));
		graphSysErrRAAYieldPi0PHOSPbPb2040Red->SetPoint(i,xValPHOStat2040[i+phosOffset], yValPHOSStat2040[i+phosOffset]);
		graphSysErrRAAYieldPi0PHOSPbPb2040Red->SetPointError(i, graphSysErrRAAYieldPi0PHOSPbPb2040->GetErrorXlow(i+phosOffset),graphSysErrRAAYieldPi0PHOSPbPb2040->GetErrorXhigh(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb2040->GetErrorYlow(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb2040->GetErrorYhigh(i+phosOffset));
	}
	graphPHOSYieldPi0SysErrPbPb2040Red->SetPoint(0,0.9, 0);
	graphPHOSYieldPi0SysErrPbPb2040Red->SetPoint(0,0.9, 0);

	// 	cout << "" << endl;
	// 	graphPHOSYieldPi0SysErrPbPb2040Red->Print();

	histoPi0PHOSPbPb2040Red = 		new TH1D("histoPi0PHOSPbPb2040Red","",18,fBinsPi0HIPtPHOS);
	for (Int_t i =1; i < 18; i++){
	//  		cout << histoPi0PHOSPbPb2040Red->GetBinCenter(i) << "\t"<< histoPi0PHOSPbPb2040->GetBinCenter(i+phosOffset) << "\t"<< histoPi0PHOSPbPb2040->GetBinContent(i+phosOffset) << "\t"<< histoPi0PHOSPbPb2040->GetBinError(i+phosOffset) << endl; 
		histoPi0PHOSPbPb2040Red->SetBinContent(i,histoPi0PHOSPbPb2040->GetBinContent(i+phosOffset));
		histoPi0PHOSPbPb2040Red->SetBinError(i,histoPi0PHOSPbPb2040->GetBinError(i+phosOffset));
	}
	histoPi0PHOSPbPb2040Red->SetBinContent(1,0);

	cout << "PHOS: 40-60%" << endl;
	histoPi0PHOSPbPb4060 = 		(TH1D*)directoryPHOSPi0PbPb4060->Get("hPi0_PbPb_cen4_NoBW_Stat");
	histoPi0PHOSSysPbPb4060 = 	(TH1D*)directoryPHOSPi0PbPb4060->Get("hPi0_PbPb_cen4_NoBW_Syst");
	graphPHOSYieldPi0SysErrPbPb4060 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb4060);	
	histoPi0PHOSSysRAAPbPb4060 = 	(TH1D*)directoryPHOSPi0PbPb4060->Get("hPi0_PbPb_cen4_SystRaa");
	graphSysErrRAAYieldPi0PHOSPbPb4060 = new TGraphAsymmErrors(histoPi0PHOSSysRAAPbPb4060);	
	histoPHOSMassData4060 =       (TH1D*)directoryPHOSPi0PbPb4060->Get("mass1_GS_Both2core_cen4");
	histoPHOSMassMC4060 =         (TH1D*)directoryPHOSPi0PbPb4060->Get("MC_mass1_GS_Both2core_cen4");
	histoPHOSWidthData4060 =      (TH1D*)directoryPHOSPi0PbPb4060->Get("width1_GS_Both2core_cen4");
	histoPHOSWidthMC4060 =        (TH1D*)directoryPHOSPi0PbPb4060->Get("MC_width1_GS_Both2core_cen4");
	histoPHOSWidthData4060->Scale(1000.);
	histoPHOSWidthMC4060->Scale(1000.);
	histoPHOSMassData4060->Scale(1000.);
	histoPHOSMassMC4060->Scale(1000.);
	for (Int_t  i = 0; i < histoPHOSWidthData4060->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthData4060->SetBinContent(i,10000.);
		histoPHOSMassData4060->SetBinContent(i,10000.);
	}
	for (Int_t  i = 0; i < histoPHOSWidthMC4060->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthMC4060->SetBinContent(i,10000.);
		histoPHOSMassMC4060->SetBinContent(i,10000.);
	}
	histoPHOSMassData4060->SetBinContent(histoPHOSMassData4060->GetNbinsX(),0);
	histoPHOSMassMC4060->SetBinContent(histoPHOSMassMC4060->GetNbinsX(),0);;
	histoPHOSWidthData4060->SetBinContent(histoPHOSWidthData4060->GetNbinsX(),10000.);
	histoPHOSWidthMC4060->SetBinContent(histoPHOSWidthMC4060->GetNbinsX(),10000.);

	graphPHOSYieldPi0SysErrPbPb4060Red = new TGraphAsymmErrors(17);
	graphSysErrRAAYieldPi0PHOSPbPb4060Red = new TGraphAsymmErrors(17);
	Double_t* xValPHOStat4060 = graphPHOSYieldPi0SysErrPbPb4060->GetX();
	Double_t* yValPHOSStat4060 = graphPHOSYieldPi0SysErrPbPb4060->GetY();
	for (Int_t i = 0; i < 17; i++){
		graphPHOSYieldPi0SysErrPbPb4060Red->SetPoint(i,xValPHOStat4060[i+phosOffset], yValPHOSStat4060[i+phosOffset]);
		graphPHOSYieldPi0SysErrPbPb4060Red->SetPointError(i, graphPHOSYieldPi0SysErrPbPb4060->GetErrorXlow(i+phosOffset),graphPHOSYieldPi0SysErrPbPb4060->GetErrorXhigh(i+phosOffset), graphPHOSYieldPi0SysErrPbPb4060->GetErrorYlow(i+phosOffset), graphPHOSYieldPi0SysErrPbPb4060->GetErrorYhigh(i+phosOffset));
		graphSysErrRAAYieldPi0PHOSPbPb4060Red->SetPoint(i,xValPHOStat4060[i+phosOffset], yValPHOSStat4060[i+phosOffset]);
		graphSysErrRAAYieldPi0PHOSPbPb4060Red->SetPointError(i, graphSysErrRAAYieldPi0PHOSPbPb4060->GetErrorXlow(i+phosOffset),graphSysErrRAAYieldPi0PHOSPbPb4060->GetErrorXhigh(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb4060->GetErrorYlow(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb4060->GetErrorYhigh(i+phosOffset));
	}
	// 	cout << "" << endl;
	// 	graphPHOSYieldPi0SysErrPbPb4060Red->Print();

	histoPi0PHOSPbPb4060Red = 		new TH1D("histoPi0PHOSPbPb4060Red","",18,fBinsPi0HIPtPHOS);
	for (Int_t i =1; i < 18; i++){
	//  		cout << histoPi0PHOSPbPb4060Red->GetBinCenter(i) << "\t"<< histoPi0PHOSPbPb4060->GetBinCenter(i+phosOffset) << "\t"<< histoPi0PHOSPbPb4060->GetBinContent(i+phosOffset) << "\t"<< histoPi0PHOSPbPb4060->GetBinError(i+phosOffset) << endl; 
		histoPi0PHOSPbPb4060Red->SetBinContent(i,histoPi0PHOSPbPb4060->GetBinContent(i+phosOffset));
		histoPi0PHOSPbPb4060Red->SetBinError(i,histoPi0PHOSPbPb4060->GetBinError(i+phosOffset));
	}

	cout << "PHOS: 60-80%" << endl;
	histoPi0PHOSPbPb6080 = 		(TH1D*)directoryPHOSPi0PbPb6080->Get("hPi0_PbPb_cen5_NoBW_Stat");
	histoPi0PHOSSysPbPb6080 = 	(TH1D*)directoryPHOSPi0PbPb6080->Get("hPi0_PbPb_cen5_NoBW_Syst");
	graphPHOSYieldPi0SysErrPbPb6080 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb6080);	
	histoPi0PHOSSysRAAPbPb6080 = 	(TH1D*)directoryPHOSPi0PbPb6080->Get("hPi0_PbPb_cen5_SystRaa");
	graphSysErrRAAYieldPi0PHOSPbPb6080 = new TGraphAsymmErrors(histoPi0PHOSSysRAAPbPb6080);	
	histoPHOSMassData6080 =       (TH1D*)directoryPHOSPi0PbPb6080->Get("mass1_GS_Both2core_cen5");
	histoPHOSMassMC6080 =         (TH1D*)directoryPHOSPi0PbPb6080->Get("MC_mass1_GS_Both2core_cen5");
	histoPHOSWidthData6080 =      (TH1D*)directoryPHOSPi0PbPb6080->Get("width1_GS_Both2core_cen5");
	histoPHOSWidthMC6080 =        (TH1D*)directoryPHOSPi0PbPb6080->Get("MC_width1_GS_Both2core_cen5");
	histoPHOSWidthData6080->Scale(1000.);
	histoPHOSWidthMC6080->Scale(1000.);
	histoPHOSMassData6080->Scale(1000.);
	histoPHOSMassMC6080->Scale(1000.);
	for (Int_t  i = 0; i < histoPHOSWidthData6080->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthData6080->SetBinContent(i,10000.);
		histoPHOSMassData6080->SetBinContent(i,10000.);
	}
	for (Int_t  i = histoPHOSWidthMC6080->GetXaxis()->FindBin(10.1); i < histoPHOSWidthMC6080->GetNbinsX()+1; i++){
		histoPHOSWidthMC6080->SetBinContent(i,10000.);
		histoPHOSMassMC6080->SetBinContent(i,10000.);
	}
	for (Int_t  i = 0; i < histoPHOSWidthMC6080->GetXaxis()->FindBin(1.0); i++){
		histoPHOSWidthMC6080->SetBinContent(i,10000.);
		histoPHOSMassMC6080->SetBinContent(i,10000.);
	}
	for (Int_t  i = histoPHOSWidthData6080->GetXaxis()->FindBin(10.1); i < histoPHOSWidthData6080->GetNbinsX()+1; i++){
		histoPHOSWidthData6080->SetBinContent(i,10000.);
		histoPHOSMassData6080->SetBinContent(i,10000.);
	}

	graphPHOSYieldPi0SysErrPbPb6080Red = new TGraphAsymmErrors(16);
	graphSysErrRAAYieldPi0PHOSPbPb6080Red = new TGraphAsymmErrors(16);
	Double_t* xValPHOStat6080 = graphPHOSYieldPi0SysErrPbPb6080->GetX();
	Double_t* yValPHOSStat6080 = graphPHOSYieldPi0SysErrPbPb6080->GetY();
	for (Int_t i = 0; i < 16; i++){
		graphPHOSYieldPi0SysErrPbPb6080Red->SetPoint(i,xValPHOStat6080[i+phosOffset], yValPHOSStat6080[i+phosOffset]);
		graphPHOSYieldPi0SysErrPbPb6080Red->SetPointError(i, graphPHOSYieldPi0SysErrPbPb6080->GetErrorXlow(i+phosOffset),graphPHOSYieldPi0SysErrPbPb6080->GetErrorXhigh(i+phosOffset), graphPHOSYieldPi0SysErrPbPb6080->GetErrorYlow(i+phosOffset), graphPHOSYieldPi0SysErrPbPb6080->GetErrorYhigh(i+phosOffset));
		graphSysErrRAAYieldPi0PHOSPbPb6080Red->SetPoint(i,xValPHOStat6080[i+phosOffset], yValPHOSStat6080[i+phosOffset]);
		graphSysErrRAAYieldPi0PHOSPbPb6080Red->SetPointError(i, graphSysErrRAAYieldPi0PHOSPbPb6080->GetErrorXlow(i+phosOffset),graphSysErrRAAYieldPi0PHOSPbPb6080->GetErrorXhigh(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb6080->GetErrorYlow(i+phosOffset), graphSysErrRAAYieldPi0PHOSPbPb6080->GetErrorYhigh(i+phosOffset));
	}
	histoPi0PHOSPbPb6080Red = 		new TH1D("histoPi0PHOSPbPb6080Red","",17,fBinsPi0HIPtPHOS);
	for (Int_t i =1; i < 17; i++){
		histoPi0PHOSPbPb6080Red->SetBinContent(i,histoPi0PHOSPbPb6080->GetBinContent(i+phosOffset));
		histoPi0PHOSPbPb6080Red->SetBinError(i,histoPi0PHOSPbPb6080->GetBinError(i+phosOffset));
	}


	filePHOSPP =       new TFile(fileNameCaloPHOSPP);
	directoryPHOSPi0PP =   (TDirectory*)filePHOSPP->Get("pp2760"); 
	histoPHOSMassDataPP = 	(TH1D*)directoryPHOSPi0PP->Get("mass2_GS_Data");
	histoPHOSWidthDataPP = 	(TH1D*)directoryPHOSPi0PP->Get("width2_GS_Data");
	histoPHOSMassMCPP = 		(TH1D*)directoryPHOSPi0PP->Get("mass2_GS_MC");
	histoPHOSWidthMCPP = 	(TH1D*)directoryPHOSPi0PP->Get("width2_GS_MC");
	histoPHOSMassDataPP->Scale(1000.);
	histoPHOSMassMCPP->Scale(1000.);
	histoPHOSWidthDataPP->Scale(1000.);
	histoPHOSWidthMCPP->Scale(1000.);
	histoPHOSMassDataPP->SetBinContent(histoPHOSMassDataPP->GetNbinsX(),0);
	histoPHOSMassDataPP->SetBinContent(1,0);
	histoPHOSMassMCPP->SetBinContent(histoPHOSMassMCPP->GetNbinsX(),0);
	histoPHOSMassMCPP->SetBinContent(1,0);
	histoPHOSWidthDataPP->SetBinContent(histoPHOSWidthDataPP->GetNbinsX(),10000.);
	histoPHOSWidthDataPP->SetBinContent(1,10000.);
	histoPHOSWidthMCPP->SetBinContent(histoPHOSWidthMCPP->GetNbinsX(),10000.);
	histoPHOSWidthMCPP->SetBinContent(1,10000.);

	maxPtPi0PbPb0005 = histoPi0PHOSPbPb0005Red->GetXaxis()->GetBinUpEdge(histoPi0PHOSPbPb0005Red->GetNbinsX()-1);
	cout << maxPtPi0PbPb0005 << endl;
	maxPtPi0PbPb0010 = histoPi0PHOSPbPb0010Red->GetXaxis()->GetBinUpEdge(histoPi0PHOSPbPb0010Red->GetNbinsX()-1);
	cout << maxPtPi0PbPb0010 << endl;

	maxPtPi0PbPb0510 = histoPi0PHOSPbPb0510Red->GetXaxis()->GetBinUpEdge(histoPi0PHOSPbPb0510Red->GetNbinsX()-1);
	cout << maxPtPi0PbPb0510 << endl;
	maxPtPi0PbPb1020 = histoPi0PHOSPbPb1020Red->GetXaxis()->GetBinUpEdge(histoPi0PHOSPbPb1020Red->GetNbinsX()-1);
	cout << maxPtPi0PbPb4060 << endl;
	cout << maxPtPi0PbPb1020 << endl;
	maxPtPi0PbPb2040 = histoPi0PHOSPbPb2040Red->GetXaxis()->GetBinUpEdge(histoPi0PHOSPbPb2040Red->GetNbinsX()-1);
	cout << maxPtPi0PbPb2040 << endl;
	maxPtPi0PbPb4060 = histoPi0PHOSPbPb4060Red->GetXaxis()->GetBinUpEdge(histoPi0PHOSPbPb4060Red->GetNbinsX()-1);
	maxPtPi0PbPb6080 = histoPi0PHOSPbPb6080Red->GetXaxis()->GetBinUpEdge(histoPi0PHOSPbPb6080Red->GetNbinsX()-1);
	cout << maxPtPi0PbPb6080 << endl;

	fileEMCALPbPb =       new TFile(fileNameCaloEMCALPbPb);
	directoryEMCALPi0PbPb0010 = (TDirectory*)fileEMCALPbPb->Get("Pi02.76TeV_PbPb_0-10");

	cout << "**************************************************************************" << endl;
	cout << "*****************************EMCAL mass data******************************" << endl;
	cout << "**************************************************************************" << endl;
	graphEMCALMassData0010 =       (TGraphErrors*)directoryEMCALPi0PbPb0010->Get("MassPi0");
	for (Int_t i = graphEMCALMassData0010->GetN()-1; i > 8; i--){  
		graphEMCALMassData0010->RemovePoint(graphEMCALMassData0010->GetN()-1);
	}
	graphEMCALMassData0010 = ScaleGraph(graphEMCALMassData0010,1000);


	cout << "**************************************************************************" << endl;
	cout << "*****************************EMCAL mass MC********************************" << endl;
	cout << "**************************************************************************" << endl;
	graphEMCALMassMC0010 =         (TGraphErrors*)directoryEMCALPi0PbPb0010->Get("TrueMassPi0");
	graphEMCALMassMC0010->Print();
	for (Int_t i = graphEMCALMassMC0010->GetN()-1; i > 8; i--){  
		graphEMCALMassMC0010->RemovePoint(graphEMCALMassMC0010->GetN()-1);
	}
	graphEMCALMassMC0010->Print();
	graphEMCALMassMC0010 = ScaleGraph(graphEMCALMassMC0010,1000);

	graphEMCALWidthData0010 =      (TGraphErrors*)directoryEMCALPi0PbPb0010->Get("FWHMPi0MeV");
	for (Int_t i = graphEMCALWidthData0010->GetN()-1; i > 8; i--){  
		graphEMCALWidthData0010->RemovePoint(graphEMCALWidthData0010->GetN()-1);
	}
	graphEMCALWidthData0010->Print();

	graphEMCALWidthMC0010 =        (TGraphErrors*)directoryEMCALPi0PbPb0010->Get("TrueFWHMPi0MeV");
	for (Int_t i = graphEMCALWidthMC0010->GetN()-1; i > 8; i--){  
		graphEMCALWidthMC0010->RemovePoint(graphEMCALWidthMC0010->GetN()-1);
	}
	graphEMCALWidthMC0010->Print();


	//    
	//    histoEMCALWidthData0005->Scale(1000.);
	//    histoEMCALWidthMC0005->Scale(1000.);
	//    histoEMCALMassData0010->Scale(1000.);
	//    histoEMCALMassMC0010->Scale(1000.);
	//    for (Int_t  i = 0; i < histoEMCALWidthData0005->GetXaxis()->FindBin(2.0); i++){
	//       histoEMCALWidthData0005->SetBinContent(i,10000.);
	//       histoEMCALMassData0005->SetBinContent(i,10000.);
	//    }
	//    for (Int_t  i = 0; i < histoEMCALWidthMC0005->GetXaxis()->FindBin(2.0); i++){
	//       histoEMCALWidthMC0005->SetBinContent(i,10000.);
	//       histoEMCALMassMC0005->SetBinContent(i,10000.);
	//    }
	//    histoEMCALMassData0005->SetBinContent(histoEMCALMassData0005->GetNbinsX(),0);
	//    histoEMCALMassMC0005->SetBinContent(histoEMCALMassMC0005->GetNbinsX(),0);;
	//    histoEMCALWidthData0005->SetBinContent(histoEMCALWidthData0005->GetNbinsX(),10000.);
	//    histoEMCALWidthMC0005->SetBinContent(histoEMCALWidthMC0005->GetNbinsX(),10000.);



	nColl0005 = GetNCollFromName("0005");
	nColl0010 = GetNCollFromName("0010");
	nColl0510 = GetNCollFromName("0510");
	nColl1020 = GetNCollFromName("1020");
	nColl2040 = GetNCollFromName("2040");
	nColl4060 = GetNCollFromName("4060");
	nColl6080 = GetNCollFromName("6080");
	nCollErr0005 = GetNCollErrFromName("0005");
	nCollErr0010 = GetNCollErrFromName("0010");
	nCollErr0510 = GetNCollErrFromName("0510");
	nCollErr1020 = GetNCollErrFromName("1020");
	nCollErr2040 = GetNCollErrFromName("2040");
	nCollErr4060 = GetNCollErrFromName("4060");
	nCollErr6080 = GetNCollErrFromName("6080");


	TFile* fileTheoryGraphs = new TFile(fileNameTheoryInput);
	Vitev_Bas_Raa_0020 = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA0020");
	cout << "Vitev_Bas_Raa_0020 " << Vitev_Bas_Raa_0020 << endl;
	Vitev_Bas_Raa_0005 = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA0005");
	cout << "Vitev_Bas_Raa_0005 " << Vitev_Bas_Raa_0005 << endl;
	Vitev_Bas_Raa_0510 = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA0510");
	cout << "Vitev_Bas_Raa_0510 " << Vitev_Bas_Raa_0510 << endl;
	Vitev_Bas_Raa_1020 = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA1020");
	cout << "Vitev_Bas_Raa_1020 " << Vitev_Bas_Raa_1020 << endl;
	Vitev_Bas_Raa_2040 = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA2040");
	cout << "Vitev_Bas_Raa_2040 " << Vitev_Bas_Raa_2040 << endl;
	Vitev_Bas_Raa_4060 = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA4060");
	cout << "Vitev_Bas_Raa_4060 " << Vitev_Bas_Raa_4060 << endl;
	Vitev_Bas_Raa_6080 = (TGraphErrors*)fileTheoryGraphs->Get("graphVitevBasRAA6080");
	cout << "Vitev_Bas_Raa_6080 " << Vitev_Bas_Raa_6080 << endl;

	// 	Vitev_ShlSel_Raa_0005 = (TGraphErrors*)Vitev_ShlSel_Raa_0020->Clone("graphVitevShlSelRAA0005");
	//    cout << "Vitev_ShlSel_Raa_0005 " << Vitev_ShlSel_Raa_0005 << endl;
	Xiao_Raa_0020 = (TGraphErrors*)fileTheoryGraphs->Get("graphXiaoRAA0020");
	Xiao_Raa_0510 = (TGraphErrors*)Xiao_Raa_0020->Clone("Xiao_Raa_0510");
	Xiao_Raa_0005 = (TGraphErrors*)Xiao_Raa_0020->Clone("Xiao_Raa_0005");
	Xiao_Raa_1020 = (TGraphErrors*)Xiao_Raa_0020->Clone("Xiao_Raa_1020");
	Xiao_Raa_2040 = (TGraphErrors*)fileTheoryGraphs->Get("graphXiaoRAA2040");
	Xiao_Raa_4060 = (TGraphErrors*)fileTheoryGraphs->Get("graphXiaoRAA4060");
	Xiao_Raa_6080 = (TGraphErrors*)fileTheoryGraphs->Get("graphXiaoRAA6080");
	gWHDG_Raa_0020 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0020");
	gWHDG_Raa_1020 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA1020");
	gWHDG_Raa_0510 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0510");
	gWHDG_Raa_0005 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0005");
	gWHDG_Raa_2040 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA2040");
	gWHDG_Raa_4060 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA4060");
	gWHDG_Raa_6080 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA6080");

	gEPOS_Spec_0005 = (TGraph*)fileTheoryGraphs->Get("graphEPOSSpecWOErr0005");
	gEPOS_Spec_0510 = (TGraph*)fileTheoryGraphs->Get("graphEPOSSpecWOErr0510");
	gEPOS_Spec_1020 = (TGraph*)fileTheoryGraphs->Get("graphEPOSSpecWOErr1020");
	gEPOS_Spec_0020 = (TGraph*)gEPOS_Spec_0005->Clone("graphEPOSSpecWOErr0020");
	gEPOS_Spec_2040 = (TGraph*)fileTheoryGraphs->Get("graphEPOSSpecWOErr2040");
	gEPOS_Spec_4060 = (TGraph*)fileTheoryGraphs->Get("graphEPOSSpecWOErr4060");
	gEPOS_Spec_6080 = (TGraph*)fileTheoryGraphs->Get("graphEPOSSpecWOErr6080");

	gKopeliovichELoss_Spec_0005 = (TGraph*)fileTheoryGraphs->Get("graphKopeliovichELossYield0005");
	gKopeliovichHydro_Spec_0005 = (TGraph*)fileTheoryGraphs->Get("graphKopeliovichHydroYield0005");
	gKopeliovichTotal_Spec_0005 = (TGraph*)fileTheoryGraphs->Get("graphKopeliovichTotalYield0005");
	gKopeliovichELoss_Spec_2040 = (TGraph*)fileTheoryGraphs->Get("graphKopeliovichELossYield2040");
	gKopeliovichHydro_Spec_2040 = (TGraph*)fileTheoryGraphs->Get("graphKopeliovichHydroYield2040");
	gKopeliovichTotal_Spec_2040 = (TGraph*)fileTheoryGraphs->Get("graphKopeliovichTotalYield2040");
	gKopeliovichELoss_Spec_6080 = (TGraph*)fileTheoryGraphs->Get("graphKopeliovichELossYield6080");
	gKopeliovichHydro_Spec_6080 = (TGraph*)fileTheoryGraphs->Get("graphKopeliovichHydroYield6080");
	gKopeliovichTotal_Spec_6080 = (TGraph*)fileTheoryGraphs->Get("graphKopeliovichTotalYield6080");
	TGraph* gKopeliovichTotal_RAA_6080 = (TGraph*)fileTheoryGraphs->Get("graphKopeliovichRAA6080_finerBinning");

	// 	cout << "WHDG graphs 0-5%" << endl;
	// 	gWHDG_Raa_0005->Print();
	// 	cout << "WHDG graphs 5-10%" << endl;
	// 	gWHDG_Raa_0510->Print();
	// 	cout << "WHDG graphs 10-20%" << endl;
	// 	gWHDG_Raa_1020->Print();


	TFile* fileDataOtherEnergies = new TFile(fileNameDataOtherEnergyInput);
	graphWA98_17_3GeVRAA_0013= (TGraphErrors*)fileDataOtherEnergies->Get("graphWA98RAA_0013");
	graphPHENIX200GeVRAA_0010= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_0010");
	graphPHENIX200GeVRAA_1020= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_1020");
	graphPHENIX200GeVRAA_0020= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_0020");
	graphPHENIX200GeVRAA_2040= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_2040");
	graphPHENIX200GeVRAA_4060= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_4060");
	graphPHENIX200GeVRAA_6080= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_6080");
	graphPHENIX39GeVRAA_0010= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_0010");
	graphPHENIX39GeVRAA_0020= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_0020");
	graphPHENIX39GeVRAA_2040= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_2040");
	graphPHENIX39GeVRAA_4060= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_4060");
	graphPHENIX39GeVRAA_6080= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_6086");
	graphPHENIX62GeVRAA_0010= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_0010");
	graphPHENIX62GeVRAA_0020= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_0020");
	graphPHENIX62GeVRAA_2040= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_2040");
	graphPHENIX62GeVRAA_4060= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_4060");
	graphPHENIX62GeVRAA_6080= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_6086");
	TGraphErrors* graphPHENIX200GeVRAAvsNPartAt7GeVc = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAAvsNPartAt7GeVc");
	TGraphErrors* graphPHENIX62GeVRAAvsNPartAt7GeVc = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAAvsNPartAt7GeVc");
	TGraphErrors* graphPHENIX39GeVRAAvsNPartAt7GeVc = (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAAvsNPartAt7GeVc");

	TFile* fileDataALICEChargedHadrons = new TFile("ExternalInputPbPb/IdentifiedCharged/PbPb_RAA_sigma_2760GeV_20120809.root");
	TGraphErrors*	graphChargedRAAStat_0005 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c0_5_stat");
	TGraphErrors*	graphChargedRAASys_0005 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c0_5_syst");
	TGraphErrors*	graphChargedRAAStat_0510 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c5_10_stat");
	TGraphErrors*	graphChargedRAASys_0510 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c5_10_syst");
	TGraphErrors*	graphChargedRAAStat_1020 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c10_20_stat");
	TGraphErrors*	graphChargedRAASys_1020 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c10_20_syst");
	TGraphErrors*	graphChargedRAAStat_2040 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c20_40_stat");
	TGraphErrors*	graphChargedRAASys_2040 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c20_40_syst");
	TGraphErrors*	graphChargedRAAStat_4060 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c40_60_stat");
	TGraphErrors*	graphChargedRAASys_4060 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c40_60_syst");
	TGraphErrors*	graphChargedRAAStat_6080 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c60_80_stat");
	TGraphErrors*	graphChargedRAASys_6080 = (TGraphErrors*)fileDataALICEChargedHadrons->Get("raa_c60_80_syst");

	TFile* fileChargedPionInputPrelim2012 = new TFile("ExternalInputPbPb/IdentifiedCharged/ChargedPionSpectraPbPb_8_Nov_2013.root");
	TH1D*	histoChargedPionRAAStat2040 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAAStat2040");
	TH1D*	histoChargedPionRAASyst2040 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAASyst2040");
	TH1D*	histoChargedPionRAAStat4060 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAAStat4060");
	TH1D*	histoChargedPionRAASyst4060 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAASyst4060");
	TH1D*	histoChargedPionRAAStat6080 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAAStat6080");
	TH1D*	histoChargedPionRAASyst6080 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAASyst6080");
	TH1D*	histoChargedPionRAAStat1020 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAAStat1020");
	TH1D*	histoChargedPionRAASyst1020 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAASyst1020");
	TH1D*	histoChargedPionRAAStat0005 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAAStat0005");
	TH1D*	histoChargedPionRAASyst0005 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAASyst0005");
	TH1D*	histoChargedPionRAAStat0510 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAAStat0510");
	TH1D*	histoChargedPionRAASyst0510 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionRAASyst0510");

	// 	TH1D*	histoChargedPionSpecHighPtStatPP = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtStatPP");
	// 	TH1D*	histoChargedPionSpecHighPtSystPP = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtSystPP");
	// 	TH1D*	histoChargedPionSpecHighPtStat0005 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtStat0005");
	// 	TH1D*	histoChargedPionSpecHighPtSyst0005 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtSyst0005");
	// 	TH1D*	histoChargedPionSpecHighPtStat0510 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtStat0510");
	// 	TH1D*	histoChargedPionSpecHighPtSyst0510 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtSyst0510");
	TH1D*	histoChargedPionSpecHighPtStat1020 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtStat1020");
	TH1D*	histoChargedPionSpecHighPtSyst1020 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtSyst1020");
	// 	TH1D*	histoChargedPionSpecHighPtSyst2040 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtSyst2040");
	// 	TH1D*	histoChargedPionSpecHighPtStat2040 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtStat2040");
	// 	TH1D*	histoChargedPionSpecHighPtSyst4060 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtSyst4060");
	// 	TH1D*	histoChargedPionSpecHighPtStat4060 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtStat4060");
	// 	TH1D*	histoChargedPionSpecHighPtSyst6080 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtSyst6080");
	// 	TH1D*	histoChargedPionSpecHighPtStat6080 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecHighPtStat6080");

	// 	TH1D*	histoChargedPionSpecLowPtStat0005 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtStat0005");
	// 	TH1D*	histoChargedPionSpecLowPtSyst0005 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtSyst0005");
	// 	TH1D*	histoChargedPionSpecLowPtStat0510 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtStat0510");
	// 	TH1D*	histoChargedPionSpecLowPtSyst0510 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtSyst0510");
	TH1D*	histoChargedPionSpecLowPtStat1020 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtStat1020");
	TH1D*	histoChargedPionSpecLowPtSyst1020 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtSyst1020");
	// 	TH1D*	histoChargedPionSpecLowPtStat2040 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtStat2040");
	// 	TH1D*	histoChargedPionSpecLowPtSyst2040 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtSyst2040");
	// 	TH1D*	histoChargedPionSpecLowPtStat4060 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtStat4060");
	// 	TH1D*	histoChargedPionSpecLowPtSyst4060 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtSyst4060");
	// 	TH1D*	histoChargedPionSpecLowPtStat6080 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtStat6080");
	// 	TH1D*	histoChargedPionSpecLowPtSyst6080 = (TH1D*)fileChargedPionInputPrelim2012->Get("histoChargedPionSpecLowPtSyst6080");

	//******************* CombinePoints Pi0 0-5 % PbPb **********************
	cout << endl << "***********************************************" << endl;
	cout << "Pi0 0-5% Spectra*********************************" << endl;
	cout << "***********************************************" << endl<< endl ;

	graphYieldPi0CombPbPb0005 = CombinePtPointsSpectra(	histoYieldPbPb0005Red, 		graphPCMYieldPi0SysErrPbPb0005Red,
												histoPi0PHOSPbPb0005Red, 	graphPHOSYieldPi0SysErrPbPb0005Red,
												graphYieldPi0CombPbPb0005StatErr, graphYieldPi0CombPbPb0005SysErr,
												xPtLimitsPbPbPHOS, 20,  1, 0,2);

	//    graphYieldPi0CombPbPb0005->RemovePoint(graphYieldPi0CombPbPb0005->GetN()-1);
	graphYieldPi0CombPbPb0005->Print();
	graphYieldPi0CombPbPb0005StatErr->Print();   
	graphYieldPi0CombPbPb0005Unshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005->Clone("graphYieldPi0CombPbPb0005Unshifted");
	graphYieldPi0CombPbPb0005StatErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005StatErr->Clone("graphYieldPi0CombPbPb0005StatErrUnshifted");

	graphYieldPi0CombPbPb0005SysErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005SysErr->Clone("graphYieldPi0CombPbPb0005SysErrUnshifted");

	Double_t parameters0005[15];
	ReturnParameterSetFittingPbPbFromString("0005",parameters0005);
	graphPCMYieldPi0PbPb0005 = new TGraphAsymmErrors(histoYieldPbPb0005Red);
	graphPCMYieldPi0PbPb0005->RemovePoint(0);
	graphPCMYieldPi0PbPb0005->RemovePoint(0);
	graphPCMYieldPi0SysErrPbPb0005Red->RemovePoint(0);

	graphPHOSYieldPi0PbPb0005 = new TGraphAsymmErrors(histoPi0PHOSPbPb0005Red);
	graphPHOSYieldPi0PbPb0005->Print();
	graphPCMYieldPi0PbPb0005Unshifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb0005->Clone("graphPCMYieldPi0PbPb0005Unshifted");
	graphPCMYieldPi0SysErrPbPb0005Unshifted= (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb0005Red->Clone("graphPCMYieldPi0SysErrPbPb0005Unshifted");
	graphPHOSYieldPi0PbPb0005Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb0005->Clone("graphPHOSYieldPi0PbPb0005Unshifted");
	graphPHOSYieldPi0SysErrPbPb0005Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb0005Red->Clone("graphPHOSYieldPi0SysErrPbPb0005Unshifted");
		
	//Shifting in X 
		fitModTsallisPi0PbPb0005 = FitObject("rad","ModTsallisPi0PbPb0005","Pi0");		
		SetParametersLimitsForFit (fitModTsallisPi0PbPb0005, 5, parameters0005);
		graphYieldPi0CombPbPb0005   = ApplyXshift(graphYieldPi0CombPbPb0005 ,fitModTsallisPi0PbPb0005);
		
		graphYieldPi0CombPbPb0005Unscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005->Clone("graphYieldPi0CombPbPb0005Unscaled");
		fitYieldPi0PbPb0005 = FitObject("rad","fitYieldPi0PbPb0005","Pi0",graphYieldPi0CombPbPb0005,0.6,maxPtPi0PbPb0005,parameters0005,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldPi0PbPb0005)<< endl;	
		TString forOutputToFile =  WriteParameterToFileLatexTable(fitYieldPi0PbPb0005,kTRUE);
		fileFinalResultsFits << forOutputToFile << endl;
		
		TF1* fitYieldDataQCDPi0PbPb0005 = FitObject("qcd","fitYieldDataQCDPi0PbPb0005","Pi0",graphYieldPi0CombPbPb0005,0.6,maxPtPi0PbPb0005,NULL,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldDataQCDPi0PbPb0005)<< endl;   
		forOutputToFile =  WriteParameterToFile(fitYieldDataQCDPi0PbPb0005);
		fileFinalResultsFits << forOutputToFile << endl;
		
		graphPCMYieldPi0PbPb0005= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0005, graphPCMYieldPi0PbPb0005, fitModTsallisPi0PbPb0005, 0, 16);
		graphPCMYieldPi0SysErrPbPb0005Red= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0005,graphPCMYieldPi0SysErrPbPb0005Red, fitModTsallisPi0PbPb0005,  0, 16);
		graphPHOSYieldPi0PbPb0005 = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0005, graphPHOSYieldPi0PbPb0005, fitModTsallisPi0PbPb0005,1,18);
		graphPHOSYieldPi0SysErrPbPb0005Red = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0005, graphPHOSYieldPi0SysErrPbPb0005Red, fitModTsallisPi0PbPb0005,1,18);
		graphYieldPi0CombPbPb0005StatErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0005, graphYieldPi0CombPbPb0005StatErr, fitModTsallisPi0PbPb0005, 0, graphYieldPi0CombPbPb0005StatErr->GetN());
		graphYieldPi0CombPbPb0005SysErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0005, graphYieldPi0CombPbPb0005SysErr, fitModTsallisPi0PbPb0005, 0, graphYieldPi0CombPbPb0005SysErr->GetN());
		
		graphYieldPi0CombPbPb0005StatErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005StatErr->Clone("graphYieldPi0CombPbPb0005StatErrUnscaled");
		graphYieldPi0CombPbPb0005SysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005SysErr->Clone("graphYieldPi0CombPbPb0005SysErrUnscaled");
		
	//Shifting in Y

	TF1* fitModTsallisPi0PbPb0005YShift = ApplyYShift(graphYieldPi0CombPbPb0005Unshifted,&graphYieldPi0CombPbPb0005YShifted, "rad", "",0.6, parameters0005,0.00001,kTRUE);
	cout << "combined binshift Y" << endl;
	graphYieldPi0CombPbPb0005SysErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005SysErrUnshifted->Clone("YShiftedCombSys0005");
	graphYieldPi0CombPbPb0005SysErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb0005SysErrYShifted, fitModTsallisPi0PbPb0005YShift);
	graphYieldPi0CombPbPb0005StatErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005StatErrUnshifted->Clone("YShiftedCombStat0005");
	graphYieldPi0CombPbPb0005StatErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb0005StatErrYShifted, fitModTsallisPi0PbPb0005YShift);		

	cout << "PCM binshift Y" << endl;
	graphYieldPi0PCMPbPb0005SysErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb0005Unshifted->Clone("YShiftedPCMSys0005");
	graphYieldPi0PCMPbPb0005SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb0005SysErrYShifted, fitModTsallisPi0PbPb0005YShift);
	graphPCMYieldPi0SysErrRAAPbPb0005Red->RemovePoint(0);
	graphYieldPi0PCMPbPb0005SysRAAErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrRAAPbPb0005Red->Clone("YShiftedPCMSysRAA0005");
	graphYieldPi0PCMPbPb0005SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb0005SysRAAErrYShifted, fitModTsallisPi0PbPb0005YShift);
	graphYieldPi0PCMPbPb0005StatErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb0005Unshifted->Clone("YShiftedPCMStat0005");
	graphYieldPi0PCMPbPb0005StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb0005StatErrYShifted, fitModTsallisPi0PbPb0005YShift);

	cout << "PHOS binshift Y" << endl;
	graphYieldPi0PHOSPbPb0005SysErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb0005Unshifted->Clone("YShiftedPHOSSys0005");
	graphYieldPi0PHOSPbPb0005SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb0005SysErrYShifted, fitModTsallisPi0PbPb0005YShift);
	graphYieldPi0PHOSPbPb0005SysRAAErrYShifted = (TGraphAsymmErrors*) graphSysErrRAAYieldPi0PHOSPbPb0005Red->Clone("YShiftedPHOSSysRAA0005");
	graphYieldPi0PHOSPbPb0005SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb0005SysRAAErrYShifted, fitModTsallisPi0PbPb0005YShift);
	graphYieldPi0PHOSPbPb0005StatErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb0005Unshifted->Clone("YShiftedPHOSStat0005");
	graphYieldPi0PHOSPbPb0005StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb0005StatErrYShifted, fitModTsallisPi0PbPb0005YShift);

	cout << endl << "***********************************************" << endl;
	cout << "Pi0 0-5% RAA*************************************" << endl;
	cout << "***********************************************" << endl<< endl;
	cout << "PCM*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;
	TF1* fitInvCrossSectionPi0Comb2760GeV = FitObject("l","fitInvCrossSectionPi0Comb2760GeV","Pi0");
	fitInvCrossSectionPi0Comb2760GeV->SetParameter(1,7.5);
	fitInvCrossSectionPi0Comb2760GeV->SetParameter(2,0.152);
	fitInvCrossSectionPi0Comb2760GeV->SetRange(0,15);

	graphInvSectionCombPi02760GeV->Fit(fitInvCrossSectionPi0Comb2760GeV,"QNRMEX0+","",0.4,8.);
	cout << "Initial Fit result" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionPi0Comb2760GeV)<< endl;

	TGraphAsymmErrors* graphRAAPCM0005;
	TGraphAsymmErrors* graphRAASysPCM0005;
	// 	graphInvSectionPCMSysPi02760GeVRed->Print();

	TCanvas* canvasDummy5 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasDummy5,  0.1, 0.01, 0.015, 0.08);
	canvasDummy5->SetLogy();
	canvasDummy5->SetLogx();
	TH2F * histo2DDummy5;
	histo2DDummy5 = new TH2F("histo2DDummy5","histo2DDummy5",1000,0.23,30.,1000,1e-8,10);
	SetStyleHistoTH2ForGraphs(histo2DDummy5, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DDummy5->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphInvSectionCombPi02760GeV, 20,markerSizeCommonSpectrumPi0900GeV, kRed, kRed, widthLinesBoxes, kTRUE);
	graphInvSectionCombPi02760GeV->Draw("pEsame");

	fitInvCrossSectionPi0Comb2760GeV->SetLineColor(kBlue+2);
	fitInvCrossSectionPi0Comb2760GeV->Draw("same");

	canvasDummy5->Update();
	canvasDummy5->Print(Form("%s/PPDataPlusFit.%s",outputDir.Data(),suffix.Data()));


	cout << "PCM*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	CalcRaa( 	graphInvSectionPCMStatPi02760GeVRed, graphInvSectionPCMSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PCMPbPb0005StatErrYShifted, graphYieldPi0PCMPbPb0005SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPCM0005, &graphRAASysPCM0005,
					nColl0005, nCollErr0005,8.,0);

	// 	return;
	cout << "PHOS*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAAPHOS0005;
	TGraphAsymmErrors* graphRAASysPHOS0005;
	CalcRaa( 	graphInvSectionPHOSStatPi02760GeVRed, graphInvSectionPHOSSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PHOSPbPb0005StatErrYShifted, graphYieldPi0PHOSPbPb0005SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPHOS0005, &graphRAASysPHOS0005,
					nColl0005, nCollErr0005, 8.,1);

	// 	return;	
	cout << "Combined*****************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAACombInd0005;
	TGraphAsymmErrors* graphRAASysCombInd0005;
	TGraphAsymmErrors* graphRAAPi0CombPbPb0005_2 = CombinePtPointsRAA(	graphRAAPCM0005, 		graphRAASysPCM0005,
												graphRAAPHOS0005, 	graphRAASysPHOS0005,
												graphRAACombInd0005, graphRAASysCombInd0005,
												xPtLimitsPbPbPHOS2, 19, 0, 0,2);

	TH1F *histoYieldPi0CombPbPb0005CombErrYShifted = new TH1F("histoYieldPi0CombPbPb0005CombErrYShifted", "", 19, xPtLimitsPbPbPHOS2) ; 
	Double_t* yvalues0005 = graphYieldPi0CombPbPb0005SysErrYShifted->GetY();
	//    graphYieldPi0CombPbPb0005SysErrYShifted->Print();
	for (Int_t i = 1; i < histoYieldPi0CombPbPb0005CombErrYShifted->GetNbinsX()+1; i++){
		histoYieldPi0CombPbPb0005CombErrYShifted->SetBinContent(i,yvalues0005[i-1]*2*TMath::Pi()* histoYieldPi0CombPbPb0005CombErrYShifted->GetBinCenter(i) );
		histoYieldPi0CombPbPb0005CombErrYShifted->SetBinError(i,TMath::Sqrt(pow(graphYieldPi0CombPbPb0005SysErrYShifted->GetErrorY(i-1),2)+pow(graphYieldPi0CombPbPb0005StatErrYShifted->GetErrorY(i-1),2)) *2*TMath::Pi()* histoYieldPi0CombPbPb0005CombErrYShifted->GetBinCenter(i) );
	}

	TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasDummy2,  0.1, 0.01, 0.015, 0.08);

	// 	canvasDummy2->SetLogy();
	canvasDummy2->SetLogx();
	TH2F * histo2DDummy2;
	histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.3,20.,1000,0,1.5	);
	SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.03,0.04, 0.03,0.04, 0.8,1., 512, 505);
	histo2DDummy2->DrawCopy(); 

	graphRAAPCM0005->Print();
	graphRAASysPCM0005->Print();

	DrawGammaSetMarkerTGraphAsym(graphRAAPCM0005, 20,markerSizeCommonSpectrum0005, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAAPCM0005->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPCM0005 , 20,markerSizeCommonSpectrum0005, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAASysPCM0005->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAAPHOS0005, 20,markerSizeCommonSpectrum0005, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAAPHOS0005->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPHOS0005 , 20,markerSizeCommonSpectrum0005, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAASysPHOS0005->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd0005, 24,markerSizeCommonSpectrum0005,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAACombInd0005->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd0005 , 20,markerSizeCommonSpectrum0005, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAASysCombInd0005->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_RAA_0005.%s",outputDir.Data(),suffix.Data()));

	canvasDummy2->SetLogy();
	canvasDummy2->SetLogx();
	TH2F * histo2DDummy3;
	histo2DDummy3 = new TH2F("histo2DDummy3","histo2DDummy3",1000,0.28,20.,1000,1e-8,7e3 );
	SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
	histo2DDummy3->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0PbPb0005, 20,markerSizeCommonSpectrum0005, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0PbPb0005->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0SysErrPbPb0005Red , 20,markerSizeCommonSpectrum0005, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0SysErrPbPb0005Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0PbPb0005, 20,markerSizeCommonSpectrum0005, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0PbPb0005->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrPbPb0005Red , 20,markerSizeCommonSpectrum0005, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0SysErrPbPb0005Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0005StatErr, 24,markerSizeCommonSpectrum0005,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb0005StatErr->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0005SysErr , 20,markerSizeCommonSpectrum0005, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb0005SysErr->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_Spectra_0005.%s",outputDir.Data(),suffix.Data()));



	graphPCMYieldPi0PbPb0005Unscaled = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb0005->Clone();	
	graphRatioCombPHOSPi0PbPb0005 = (TGraphAsymmErrors*) graphPHOSYieldPi0PbPb0005->Clone();	
	graphRatioCombPHOSPi0PbPb0005Sys = (TGraphAsymmErrors*) graphPHOSYieldPi0SysErrPbPb0005Red->Clone();	
	graphRatioCombConvPi0PbPb0005 = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb0005->Clone();	
	graphRatioCombConvPi0PbPb0005Sys = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrPbPb0005Red->Clone();	
	graphRatioCombPHOSPi0PbPb0005 = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb0005, fitYieldDataQCDPi0PbPb0005); 
	graphRatioCombPHOSPi0PbPb0005Sys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb0005Sys, fitYieldDataQCDPi0PbPb0005); 
	graphRatioCombConvPi0PbPb0005 = CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb0005, fitYieldDataQCDPi0PbPb0005); 
	graphRatioCombConvPi0PbPb0005Sys= CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb0005Sys, fitYieldDataQCDPi0PbPb0005); 
	graphRatioCombCombFitPbPb0005 = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005->Clone();
	graphRatioCombCombFitPbPb0005 = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb0005, fitYieldDataQCDPi0PbPb0005); 
	graphRatioCombCombFitPbPb0005Stat = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005StatErr->Clone();
	graphRatioCombCombFitPbPb0005Stat = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb0005Stat, fitYieldDataQCDPi0PbPb0005); 
	graphRatioCombCombFitPbPb0005Sys = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005SysErr->Clone();
	graphRatioCombCombFitPbPb0005Sys = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb0005Sys, fitYieldDataQCDPi0PbPb0005); 



	cout << endl << "***********************************************" << endl;
	cout << "Pi0 0-10% Spectra*********************************" << endl;
	cout << "***********************************************" << endl<< endl ;
	// 	graphPCMYieldPi0SysErrPbPb1020->Print();
	// 	graphPCMYieldPi0SysErrPbPb1020->RemovePoint(graphPCMYieldPi0SysErrPbPb1020->GetN()-1);
	graphYieldPi0CombPbPb0010 = CombinePtPointsSpectra(	histoYieldPbPb0010Red, 		graphPCMYieldPi0SysErrPbPb0010Red,
												histoPi0PHOSPbPb0010Red, 	graphPHOSYieldPi0SysErrPbPb0010Red,
												graphYieldPi0CombPbPb0010StatErr, graphYieldPi0CombPbPb0010SysErr,
												xPtLimitsPbPbPHOS, 20,  1, 0,2);

	graphYieldPi0CombPbPb0010Unshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010->Clone("graphYieldPi0CombPbPb0010Unshifted");
	graphYieldPi0CombPbPb0010StatErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010StatErr->Clone("graphYieldPi0CombPbPb0010StatErrUnshifted");
	graphYieldPi0CombPbPb0010SysErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010SysErr->Clone("graphYieldPi0CombPbPb0010SysErrUnshifted");

	Double_t parameters0010[15];
	ReturnParameterSetFittingPbPbFromString("0010",parameters0010);
	graphPCMYieldPi0PbPb0010 = new TGraphAsymmErrors(histoYieldPbPb0010Red);
	graphPCMYieldPi0PbPb0010->RemovePoint(0);
	graphPCMYieldPi0PbPb0010->RemovePoint(0);
	graphPCMYieldPi0SysErrPbPb0010Red->RemovePoint(0);
	graphPHOSYieldPi0PbPb0010 = new TGraphAsymmErrors(histoPi0PHOSPbPb0010Red);
	graphPCMYieldPi0PbPb0010Unshifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb0010->Clone("graphPCMYieldPi0PbPb0010Unshifted");
	graphPCMYieldPi0SysErrPbPb0010Unshifted= (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb0010Red->Clone("graphPCMYieldPi0SysErrPbPb0010Unshifted");
	graphPHOSYieldPi0PbPb0010Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb0010->Clone("graphPHOSYieldPi0PbPb0010Unshifted");
	graphPHOSYieldPi0SysErrPbPb0010Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb0010Red->Clone("graphPHOSYieldPi0SysErrPbPb0010Unshifted");
		
	//Shifting in X 
		fitModTsallisPi0PbPb0010 = FitObject("rad","ModTsallisPi0PbPb0010","Pi0");		
		SetParametersLimitsForFit (fitModTsallisPi0PbPb0010, 5, parameters0010);
		graphYieldPi0CombPbPb0010   = ApplyXshift(graphYieldPi0CombPbPb0010 ,fitModTsallisPi0PbPb0010);
		
		graphYieldPi0CombPbPb0010Unscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010->Clone("graphYieldPi0CombPbPb0010Unscaled");
		fitYieldPi0PbPb0010 = FitObject("rad","fitYieldPi0PbPb0010","Pi0",graphYieldPi0CombPbPb0010,0.6,maxPtPi0PbPb0010,parameters0010,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldPi0PbPb0010)<< endl;	
		forOutputToFile =  WriteParameterToFileLatexTable(fitYieldPi0PbPb0010,kTRUE);
		fileFinalResultsFits << forOutputToFile << endl;

		TF1* fitYieldDataQCDPi0PbPb0010 = FitObject("qcd","fitYieldDataQCDPi0PbPb0010","Pi0",graphYieldPi0CombPbPb0010,0.6,maxPtPi0PbPb0010,NULL,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldDataQCDPi0PbPb0010)<< endl;   
		forOutputToFile =  WriteParameterToFile(fitYieldDataQCDPi0PbPb0010);
		fileFinalResultsFits << forOutputToFile << endl;
		
		
		graphPCMYieldPi0PbPb0010= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0010, graphPCMYieldPi0PbPb0010, fitModTsallisPi0PbPb0010, 0, 16);		
		graphPCMYieldPi0SysErrPbPb0010Red= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0010,graphPCMYieldPi0SysErrPbPb0010Red, fitModTsallisPi0PbPb0010,  0, 16);
		graphPHOSYieldPi0PbPb0010 = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0010, graphPHOSYieldPi0PbPb0010, fitModTsallisPi0PbPb0010,1,18);
		graphPHOSYieldPi0SysErrPbPb0010Red = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0010, graphPHOSYieldPi0SysErrPbPb0010Red, fitModTsallisPi0PbPb0010,1,18);
		graphYieldPi0CombPbPb0010StatErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0010, graphYieldPi0CombPbPb0010StatErr, fitModTsallisPi0PbPb0010, 0, graphYieldPi0CombPbPb0010StatErr->GetN());
		graphYieldPi0CombPbPb0010SysErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0010, graphYieldPi0CombPbPb0010SysErr, fitModTsallisPi0PbPb0010, 0, graphYieldPi0CombPbPb0010SysErr->GetN());
		
		graphYieldPi0CombPbPb0010StatErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010StatErr->Clone("graphYieldPi0CombPbPb0010StatErrUnscaled");
		graphYieldPi0CombPbPb0010SysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010SysErr->Clone("graphYieldPi0CombPbPb0010SysErrUnscaled");
		
	//Shifting in Y

	TF1* fitModTsallisPi0PbPb0010YShift = ApplyYShift(graphYieldPi0CombPbPb0010Unshifted,&graphYieldPi0CombPbPb0010YShifted, "rad", "",0.6, parameters0010,0.00001,kTRUE);
	cout << "combined binshift Y" << endl;
	graphYieldPi0CombPbPb0010SysErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010SysErrUnshifted->Clone("YShiftedCombSys0010");
	graphYieldPi0CombPbPb0010SysErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb0010SysErrYShifted, fitModTsallisPi0PbPb0010YShift);
	graphYieldPi0CombPbPb0010StatErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010StatErrUnshifted->Clone("YShiftedCombStat0010");
	graphYieldPi0CombPbPb0010StatErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb0010StatErrYShifted, fitModTsallisPi0PbPb0010YShift);		

	cout << "PCM binshift Y" << endl;
	graphYieldPi0PCMPbPb0010SysErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb0010Unshifted->Clone("YShiftedPCMSys0010");
	graphYieldPi0PCMPbPb0010SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb0010SysErrYShifted, fitModTsallisPi0PbPb0010YShift);

	graphYieldPi0PCMPbPb0010SysRAAErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrRAAPbPb0010Red->Clone("YShiftedPCMSysRAA0010");
	graphYieldPi0PCMPbPb0010SysRAAErrYShifted->RemovePoint(0);
	graphYieldPi0PCMPbPb0010SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb0010SysRAAErrYShifted, fitModTsallisPi0PbPb0010YShift);
	graphYieldPi0PCMPbPb0010StatErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb0010Unshifted->Clone("YShiftedPCMStat0010");
	graphYieldPi0PCMPbPb0010StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb0010StatErrYShifted, fitModTsallisPi0PbPb0010YShift);

	cout << "PHOS binshift Y" << endl;
	graphYieldPi0PHOSPbPb0010SysErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb0010Unshifted->Clone("YShiftedPHOSSys0010");
	graphYieldPi0PHOSPbPb0010SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb0010SysErrYShifted, fitModTsallisPi0PbPb0010YShift);
	graphYieldPi0PHOSPbPb0010SysRAAErrYShifted = (TGraphAsymmErrors*) graphSysErrRAAYieldPi0PHOSPbPb0010Red->Clone("YShiftedPHOSSysRAA0010");
	graphYieldPi0PHOSPbPb0010SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb0010SysRAAErrYShifted, fitModTsallisPi0PbPb0010YShift);
	graphYieldPi0PHOSPbPb0010StatErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb0010Unshifted->Clone("YShiftedPHOSStat0010");
	graphYieldPi0PHOSPbPb0010StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb0010StatErrYShifted, fitModTsallisPi0PbPb0010YShift);
	cout << endl << "***********************************************" << endl;
	cout << "Pi0 0-10% RAA*************************************" << endl;
	cout << "***********************************************" << endl<< endl;
	cout << "PCM*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	graphInvSectionCombPi02760GeV->Fit(fitInvCrossSectionPi0Comb2760GeV,"QNRMEX0+","",0.4,8.);
	cout << "Initial Fit result" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionPi0Comb2760GeV)<< endl;

	TGraphAsymmErrors* graphRAAPCM0010;
	TGraphAsymmErrors* graphRAASysPCM0010;
	CalcRaa( 	graphInvSectionPCMStatPi02760GeVRed, graphInvSectionPCMSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PCMPbPb0010StatErrYShifted, graphYieldPi0PCMPbPb0010SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPCM0010, &graphRAASysPCM0010,
					nColl0010, nCollErr0010,8.,0);

	cout << "PHOS*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAAPHOS0010;
	TGraphAsymmErrors* graphRAASysPHOS0010;
	CalcRaa( 	graphInvSectionPHOSStatPi02760GeVRed, graphInvSectionPHOSSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PHOSPbPb0010StatErrYShifted, graphYieldPi0PHOSPbPb0010SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPHOS0010, &graphRAASysPHOS0010,
					nColl0010, nCollErr0010, 8.,1);

	cout << "Combined*****************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAACombInd0010;
	TGraphAsymmErrors* graphRAASysCombInd0010;
	TGraphAsymmErrors* graphRAAPi0CombPbPb0010_2 = CombinePtPointsRAA(	graphRAAPCM0010, 		graphRAASysPCM0010,
												graphRAAPHOS0010, 	graphRAASysPHOS0010,
												graphRAACombInd0010, graphRAASysCombInd0010,
												xPtLimitsPbPbPHOS2, 19, 0, 0,2);

	TH1F *histoYieldPi0CombPbPb0010CombErrYShifted = new TH1F("histoYieldPi0CombPbPb0010CombErrYShifted", "", 19, xPtLimitsPbPbPHOS2) ; 
	Double_t* yvalues0010 = graphYieldPi0CombPbPb0010SysErrYShifted->GetY();
	//    graphYieldPi0CombPbPb0010SysErrYShifted->Print();
	for (Int_t i = 1; i < histoYieldPi0CombPbPb0010CombErrYShifted->GetNbinsX()+1; i++){
		histoYieldPi0CombPbPb0010CombErrYShifted->SetBinContent(i,yvalues0010[i-1]*2*TMath::Pi()* histoYieldPi0CombPbPb0010CombErrYShifted->GetBinCenter(i) );
		histoYieldPi0CombPbPb0010CombErrYShifted->SetBinError(i,TMath::Sqrt(pow(graphYieldPi0CombPbPb0010SysErrYShifted->GetErrorY(i-1),2)+pow(graphYieldPi0CombPbPb0010StatErrYShifted->GetErrorY(i-1),2)) *2*TMath::Pi()* histoYieldPi0CombPbPb0010CombErrYShifted->GetBinCenter(i) );
	}

	canvasDummy2->cd();
	histo2DDummy2->DrawCopy();
	canvasDummy2->SetLogy(0);

	DrawGammaSetMarkerTGraphAsym(graphRAAPCM0010, 20,markerSizeCommonSpectrum0010, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAAPCM0010->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPCM0010 , 20,markerSizeCommonSpectrum0010, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAASysPCM0010->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAAPHOS0010, 20,markerSizeCommonSpectrum0010, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAAPHOS0010->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPHOS0010 , 20,markerSizeCommonSpectrum0010, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAASysPHOS0010->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd0010, 24,markerSizeCommonSpectrum0010,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAACombInd0010->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd0010 , 20,markerSizeCommonSpectrum0010, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAASysCombInd0010->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_RAA_0010.%s",outputDir.Data(),suffix.Data()));

	canvasDummy2->SetLogy();
	histo2DDummy3->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0PbPb0010, 20,markerSizeCommonSpectrum0010, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0PbPb0010->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0SysErrPbPb0010Red , 20,markerSizeCommonSpectrum0010, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0SysErrPbPb0010Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0PbPb0010, 20,markerSizeCommonSpectrum0010, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0PbPb0010->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrPbPb0010Red , 20,markerSizeCommonSpectrum0010, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0SysErrPbPb0010Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0010StatErr, 24,markerSizeCommonSpectrum0010,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb0010StatErr->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0010SysErr , 20,markerSizeCommonSpectrum0010, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb0010SysErr->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_Spectra_0010.%s",outputDir.Data(),suffix.Data()));


	graphPCMYieldPi0PbPb0010Unscaled = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb0010->Clone();	
	graphRatioCombPHOSPi0PbPb0010 = (TGraphAsymmErrors*) graphPHOSYieldPi0PbPb0010->Clone();	
	graphRatioCombPHOSPi0PbPb0010Sys = (TGraphAsymmErrors*) graphPHOSYieldPi0SysErrPbPb0010Red->Clone();	
	graphRatioCombConvPi0PbPb0010 = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb0010->Clone();	
	graphRatioCombConvPi0PbPb0010Sys = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrPbPb0010Red->Clone();	
	graphRatioCombPHOSPi0PbPb0010 = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb0010, fitYieldDataQCDPi0PbPb0010); 
	graphRatioCombPHOSPi0PbPb0010Sys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb0010Sys, fitYieldDataQCDPi0PbPb0010); 
	graphRatioCombConvPi0PbPb0010 = CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb0010, fitYieldDataQCDPi0PbPb0010); 
	graphRatioCombConvPi0PbPb0010Sys= CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb0010Sys, fitYieldDataQCDPi0PbPb0010); 
	graphRatioCombCombFitPbPb0010 = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010->Clone();
	graphRatioCombCombFitPbPb0010 = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb0010, fitYieldDataQCDPi0PbPb0010); 
	graphRatioCombCombFitPbPb0010Stat = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010StatErr->Clone();
	graphRatioCombCombFitPbPb0010Stat = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb0010Stat, fitYieldDataQCDPi0PbPb0010); 
	graphRatioCombCombFitPbPb0010Sys = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010SysErr->Clone();
	graphRatioCombCombFitPbPb0010Sys = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb0010Sys, fitYieldDataQCDPi0PbPb0010); 



	//******************* CombinePoints Pi0 5-10 % PbPb **********************
	cout << endl << "***********************************************" << endl;
	cout << "Pi0 5-10% Spectra*********************************" << endl;
	cout << "***********************************************" << endl<< endl ;
	// 	graphPCMYieldPi0SysErrPbPb0510->Print();
	// 	graphPCMYieldPi0SysErrPbPb0510->RemovePoint(graphPCMYieldPi0SysErrPbPb0510->GetN()-1);
	graphYieldPi0CombPbPb0510 = CombinePtPointsSpectra(	histoYieldPbPb0510Red, 		graphPCMYieldPi0SysErrPbPb0510Red,
												histoPi0PHOSPbPb0510Red, 	graphPHOSYieldPi0SysErrPbPb0510Red,
												graphYieldPi0CombPbPb0510StatErr, graphYieldPi0CombPbPb0510SysErr,
												xPtLimitsPbPbPHOS, 20,  1, 0,2);


	graphYieldPi0CombPbPb0510Unshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510->Clone("graphYieldPi0CombPbPb0510Unshifted");
	graphYieldPi0CombPbPb0510StatErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510StatErr->Clone("graphYieldPi0CombPbPb0510StatErrUnshifted");
	graphYieldPi0CombPbPb0510SysErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510SysErr->Clone("graphYieldPi0CombPbPb0510SysErrUnshifted");

	Double_t parameters0510[15];
	ReturnParameterSetFittingPbPbFromString("0020",parameters0510);
	graphPCMYieldPi0PbPb0510 = new TGraphAsymmErrors(histoYieldPbPb0510Red);
	graphPCMYieldPi0PbPb0510->RemovePoint(0);
	graphPCMYieldPi0PbPb0510->RemovePoint(0);
	graphPCMYieldPi0SysErrPbPb0510Red->RemovePoint(0);
	graphPHOSYieldPi0PbPb0510 = new TGraphAsymmErrors(histoPi0PHOSPbPb0510Red);
	graphPCMYieldPi0PbPb0510Unshifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb0510->Clone("graphPCMYieldPi0PbPb0510Unshifted");
	graphPCMYieldPi0SysErrPbPb0510Unshifted= (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb0510Red->Clone("graphPCMYieldPi0SysErrPbPb0510Unshifted");
	graphPHOSYieldPi0PbPb0510Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb0510->Clone("graphPHOSYieldPi0PbPb0510Unshifted");
	graphPHOSYieldPi0SysErrPbPb0510Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb0510Red->Clone("graphPHOSYieldPi0SysErrPbPb0510Unshifted");
		
	//Shifting in X 
		fitModTsallisPi0PbPb0510 = FitObject("rad","ModTsallisPi0PbPb0510","Pi0");		
		SetParametersLimitsForFit (fitModTsallisPi0PbPb0510, 5, parameters0510);
		graphYieldPi0CombPbPb0510   = ApplyXshift(graphYieldPi0CombPbPb0510 ,fitModTsallisPi0PbPb0510);
		
		graphYieldPi0CombPbPb0510Unscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510->Clone("graphYieldPi0CombPbPb0510Unscaled");
		fitYieldPi0PbPb0510 = FitObject("rad","fitYieldPi0PbPb0510","Pi0",graphYieldPi0CombPbPb0510,0.6,maxPtPi0PbPb0510,parameters0510,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldPi0PbPb0510)<< endl;	
		forOutputToFile =  WriteParameterToFileLatexTable(fitYieldPi0PbPb0510,kTRUE);
		fileFinalResultsFits << forOutputToFile << endl;

		TF1* fitYieldDataQCDPi0PbPb0510 = FitObject("qcd","fitYieldDataQCDPi0PbPb0510","Pi0",graphYieldPi0CombPbPb0510,0.6,maxPtPi0PbPb0510,NULL,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldDataQCDPi0PbPb0510)<< endl;   
		forOutputToFile =  WriteParameterToFile(fitYieldDataQCDPi0PbPb0510);
		fileFinalResultsFits << forOutputToFile << endl;
		
		graphPCMYieldPi0PbPb0510= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0510, graphPCMYieldPi0PbPb0510, fitModTsallisPi0PbPb0510, 0, 15);
		graphPCMYieldPi0SysErrPbPb0510Red= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0510,graphPCMYieldPi0SysErrPbPb0510Red, fitModTsallisPi0PbPb0510,  0, 15);
		graphPHOSYieldPi0PbPb0510 = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0510, graphPHOSYieldPi0PbPb0510, fitModTsallisPi0PbPb0510,1,18);
		graphPHOSYieldPi0SysErrPbPb0510Red = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0510, graphPHOSYieldPi0SysErrPbPb0510Red, fitModTsallisPi0PbPb0510,1,18);
		
		graphYieldPi0CombPbPb0510StatErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0510, graphYieldPi0CombPbPb0510StatErr, fitModTsallisPi0PbPb0510, 0, graphYieldPi0CombPbPb0510StatErr->GetN());
		graphYieldPi0CombPbPb0510SysErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb0510, graphYieldPi0CombPbPb0510SysErr, fitModTsallisPi0PbPb0510, 0, graphYieldPi0CombPbPb0510SysErr->GetN());
		
		graphYieldPi0CombPbPb0510StatErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510StatErr->Clone("graphYieldPi0CombPbPb0510StatErrUnscaled");
		graphYieldPi0CombPbPb0510SysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510SysErr->Clone("graphYieldPi0CombPbPb0510SysErrUnscaled");
		
	//Shifting in Y

	TF1* fitModTsallisPi0PbPb0510YShift = ApplyYShift(graphYieldPi0CombPbPb0510Unshifted,&graphYieldPi0CombPbPb0510YShifted, "rad", "",0.3, parameters0510,0.00001,kTRUE);
	graphYieldPi0CombPbPb0510SysErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510SysErrUnshifted->Clone("YShiftedCombSys0510");
	graphYieldPi0CombPbPb0510SysErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb0510SysErrYShifted, fitModTsallisPi0PbPb0510YShift);
	graphYieldPi0CombPbPb0510StatErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510StatErrUnshifted->Clone("YShiftedCombStat0510");
	graphYieldPi0CombPbPb0510StatErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb0510StatErrYShifted, fitModTsallisPi0PbPb0510YShift);		

	graphYieldPi0PCMPbPb0510SysErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb0510Unshifted->Clone("YShiftedPCMSys0510");
	graphYieldPi0PCMPbPb0510SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb0510SysErrYShifted, fitModTsallisPi0PbPb0510YShift);
	graphYieldPi0PCMPbPb0510SysRAAErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrRAAPbPb0510Red->Clone("YShiftedPCMSysRAA0510");
	graphYieldPi0PCMPbPb0510SysRAAErrYShifted->RemovePoint(0);
	graphYieldPi0PCMPbPb0510SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb0510SysRAAErrYShifted, fitModTsallisPi0PbPb0510YShift);
	graphYieldPi0PCMPbPb0510StatErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb0510Unshifted->Clone("YShiftedPCMStat0510");
	graphYieldPi0PCMPbPb0510StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb0510StatErrYShifted, fitModTsallisPi0PbPb0510YShift);

	graphYieldPi0PHOSPbPb0510SysErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb0510Unshifted->Clone("YShiftedPHOSSys0510");
	graphYieldPi0PHOSPbPb0510SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb0510SysErrYShifted, fitModTsallisPi0PbPb0510YShift);
	graphYieldPi0PHOSPbPb0510SysRAAErrYShifted = (TGraphAsymmErrors*) graphSysErrRAAYieldPi0PHOSPbPb0510Red->Clone("YShiftedPHOSSysRAA0510");
	graphYieldPi0PHOSPbPb0510SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb0510SysRAAErrYShifted, fitModTsallisPi0PbPb0510YShift);
	graphYieldPi0PHOSPbPb0510StatErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb0510Unshifted->Clone("YShiftedPHOSStat0510");
	graphYieldPi0PHOSPbPb0510StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb0510StatErrYShifted, fitModTsallisPi0PbPb0510YShift);

	cout << endl << "***********************************************" << endl;
	cout << "Pi0 5-10% RAA*************************************" << endl;
	cout << "***********************************************" << endl<< endl;
	cout << "PCM*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	graphInvSectionCombPi02760GeV->Fit(fitInvCrossSectionPi0Comb2760GeV,"QNRMEX0+","",0.4,8.);
	cout << "Initial Fit result" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionPi0Comb2760GeV)<< endl;


	TGraphAsymmErrors* graphRAAPCM0510;
	TGraphAsymmErrors* graphRAASysPCM0510;
	CalcRaa( 	graphInvSectionPCMStatPi02760GeVRed, graphInvSectionPCMSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PCMPbPb0510StatErrYShifted, graphYieldPi0PCMPbPb0510SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPCM0510, &graphRAASysPCM0510,
					nColl0510, nCollErr0510,8.,0);

	cout << "PHOS*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;


	TGraphAsymmErrors* graphRAAPHOS0510;
	TGraphAsymmErrors* graphRAASysPHOS0510;
	CalcRaa( 	graphInvSectionPHOSStatPi02760GeVRed, graphInvSectionPHOSSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PHOSPbPb0510StatErrYShifted, graphYieldPi0PHOSPbPb0510SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPHOS0510, &graphRAASysPHOS0510,
					nColl0510, nCollErr0510, 8.,1);

	cout << "Combined*****************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAACombInd0510;
	TGraphAsymmErrors* graphRAASysCombInd0510;
	TGraphAsymmErrors* graphRAAPi0CombPbPb0510_2 = CombinePtPointsRAA(	graphRAAPCM0510, 		graphRAASysPCM0510,
												graphRAAPHOS0510, 	graphRAASysPHOS0510,
												graphRAACombInd0510, graphRAASysCombInd0510,
												xPtLimitsPbPbPHOS2, 19, 0, 0,2);

	TH1F *histoYieldPi0CombPbPb0510CombErrYShifted = new TH1F("histoYieldPi0CombPbPb0510CombErrYShifted", "", 19, xPtLimitsPbPbPHOS2) ; 
	Double_t* yvalues0510 = graphYieldPi0CombPbPb0510SysErrYShifted->GetY();
	graphYieldPi0CombPbPb0510SysErrYShifted->Print();
	for (Int_t i = 1; i < histoYieldPi0CombPbPb0510CombErrYShifted->GetNbinsX()+1; i++){
		histoYieldPi0CombPbPb0510CombErrYShifted->SetBinContent(i,yvalues0510[i-1]*2*TMath::Pi()* histoYieldPi0CombPbPb0510CombErrYShifted->GetBinCenter(i) );
		histoYieldPi0CombPbPb0510CombErrYShifted->SetBinError(i,TMath::Sqrt(pow(graphYieldPi0CombPbPb0510SysErrYShifted->GetErrorY(i-1),2)+pow(graphYieldPi0CombPbPb0510StatErrYShifted->GetErrorY(i-1),2)) *2*TMath::Pi()* histoYieldPi0CombPbPb0510CombErrYShifted->GetBinCenter(i) );
	}
	canvasDummy2->cd();
	canvasDummy2->SetLogy(0);
	histo2DDummy2->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphRAAPCM0510, 20,markerSizeCommonSpectrum0510, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAAPCM0510->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPCM0510 , 20,markerSizeCommonSpectrum0510, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAASysPCM0510->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAAPHOS0510, 20,markerSizeCommonSpectrum0510, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAAPHOS0510->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPHOS0510 , 20,markerSizeCommonSpectrum0510, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAASysPHOS0510->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd0510, 24,markerSizeCommonSpectrum0510,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAACombInd0510->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd0510 , 20,markerSizeCommonSpectrum0510, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAASysCombInd0510->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_RAA_0510.%s",outputDir.Data(),suffix.Data()));

	canvasDummy2->SetLogy();
	histo2DDummy3->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0PbPb0510, 20,markerSizeCommonSpectrum0510, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0PbPb0510->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0SysErrPbPb0510Red , 20,markerSizeCommonSpectrum0510, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0SysErrPbPb0510Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0PbPb0510, 20,markerSizeCommonSpectrum0510, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0PbPb0510->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrPbPb0510Red , 20,markerSizeCommonSpectrum0510, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0SysErrPbPb0510Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0510StatErr, 24,markerSizeCommonSpectrum0510,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb0510StatErr->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0510SysErr , 20,markerSizeCommonSpectrum0510, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb0510SysErr->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_Spectra_0510.%s",outputDir.Data(),suffix.Data()));

	cout << WriteParameterToFile(fitModTsallisPi0PbPb0510YShift)<< endl;	
	forOutputToFile =  WriteParameterToFileLatexTable(fitModTsallisPi0PbPb0510YShift,kTRUE);
	fileFinalResultsFits << forOutputToFile << endl;


	graphPCMYieldPi0PbPb0510Unscaled = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb0510->Clone();	
	graphRatioCombPHOSPi0PbPb0510 = (TGraphAsymmErrors*) graphPHOSYieldPi0PbPb0510->Clone();	
	graphRatioCombPHOSPi0PbPb0510Sys = (TGraphAsymmErrors*) graphPHOSYieldPi0SysErrPbPb0510Red->Clone();	
	graphRatioCombConvPi0PbPb0510 = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb0510->Clone();	
	graphRatioCombConvPi0PbPb0510Sys = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrPbPb0510Red->Clone();	
	graphRatioCombPHOSPi0PbPb0510 = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb0510, fitYieldDataQCDPi0PbPb0510); 
	graphRatioCombPHOSPi0PbPb0510Sys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb0510Sys, fitYieldDataQCDPi0PbPb0510); 
	graphRatioCombConvPi0PbPb0510 = CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb0510, fitYieldDataQCDPi0PbPb0510); 
	// 	return;
	graphRatioCombConvPi0PbPb0510Sys= CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb0510Sys, fitYieldDataQCDPi0PbPb0510); 
	graphRatioCombCombFitPbPb0510 = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510->Clone();
	graphRatioCombCombFitPbPb0510 = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb0510, fitYieldDataQCDPi0PbPb0510); 
	graphRatioCombCombFitPbPb0510Stat = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510StatErr->Clone();
	graphRatioCombCombFitPbPb0510Stat = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb0510Stat, fitYieldDataQCDPi0PbPb0510); 
	graphRatioCombCombFitPbPb0510Sys = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510SysErr->Clone();
	graphRatioCombCombFitPbPb0510Sys = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb0510Sys, fitYieldDataQCDPi0PbPb0510); 


	//******************* CombinePoints Pi0 10-20 % PbPb **********************
	cout << endl << "***********************************************" << endl;
	cout << "Pi0 10-20% Spectra*********************************" << endl;
	cout << "***********************************************" << endl<< endl ;
	graphYieldPi0CombPbPb1020 = CombinePtPointsSpectra(	histoYieldPbPb1020Red, 		graphPCMYieldPi0SysErrPbPb1020Red,
												histoPi0PHOSPbPb1020Red, 	graphPHOSYieldPi0SysErrPbPb1020Red,
												graphYieldPi0CombPbPb1020StatErr, graphYieldPi0CombPbPb1020SysErr,
												xPtLimitsPbPbPHOS, 20,  1, 0,2);

	graphYieldPi0CombPbPb1020Unshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020->Clone("graphYieldPi0CombPbPb1020Unshifted");
	graphYieldPi0CombPbPb1020StatErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020StatErr->Clone("graphYieldPi0CombPbPb1020StatErrUnshifted");
	graphYieldPi0CombPbPb1020SysErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020SysErr->Clone("graphYieldPi0CombPbPb1020SysErrUnshifted");

	Double_t parameters1020[15];
	ReturnParameterSetFittingPbPbFromString("0020",parameters1020);
	graphPCMYieldPi0PbPb1020 = new TGraphAsymmErrors(histoYieldPbPb1020Red);
	graphPCMYieldPi0PbPb1020->RemovePoint(0);
	graphPCMYieldPi0PbPb1020->RemovePoint(0);
	graphPCMYieldPi0SysErrPbPb1020Red->RemovePoint(0);
	graphPHOSYieldPi0PbPb1020 = new TGraphAsymmErrors(histoPi0PHOSPbPb1020Red);
	graphPCMYieldPi0PbPb1020Unshifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb1020->Clone("graphPCMYieldPi0PbPb1020Unshifted");
	graphPCMYieldPi0SysErrPbPb1020Unshifted= (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb1020Red->Clone("graphPCMYieldPi0SysErrPbPb1020Unshifted");
	graphPHOSYieldPi0PbPb1020Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb1020->Clone("graphPHOSYieldPi0PbPb1020Unshifted");
	graphPHOSYieldPi0SysErrPbPb1020Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb1020Red->Clone("graphPHOSYieldPi0SysErrPbPb1020Unshifted");
		
	//Shifting in X 
		fitModTsallisPi0PbPb1020 = FitObject("rad","ModTsallisPi0PbPb1020","Pi0");		
		SetParametersLimitsForFit (fitModTsallisPi0PbPb1020, 5, parameters1020);
		graphYieldPi0CombPbPb1020   = ApplyXshift(graphYieldPi0CombPbPb1020 ,fitModTsallisPi0PbPb1020);
		
		graphYieldPi0CombPbPb1020Unscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020->Clone("graphYieldPi0CombPbPb1020Unscaled");
		fitYieldPi0PbPb1020 = FitObject("rad","fitYieldPi0PbPb1020","Pi0",graphYieldPi0CombPbPb1020,0.6,maxPtPi0PbPb1020,parameters1020,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldPi0PbPb1020)<< endl;	
		forOutputToFile =  WriteParameterToFileLatexTable(fitYieldPi0PbPb1020,kTRUE);
		fileFinalResultsFits << forOutputToFile << endl;

		TF1* fitYieldDataQCDPi0PbPb1020 = FitObject("qcd","fitYieldDataQCDPi0PbPb1020","Pi0",graphYieldPi0CombPbPb1020,0.6,maxPtPi0PbPb1020,NULL,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldDataQCDPi0PbPb1020)<< endl;   
		forOutputToFile =  WriteParameterToFile(fitYieldDataQCDPi0PbPb1020);
		fileFinalResultsFits << forOutputToFile << endl;
		
				
		graphPCMYieldPi0PbPb1020= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb1020, graphPCMYieldPi0PbPb1020, fitModTsallisPi0PbPb1020, 0, 15);
		graphPCMYieldPi0SysErrPbPb1020Red= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb1020,graphPCMYieldPi0SysErrPbPb1020Red, fitModTsallisPi0PbPb1020,  0, 15);
		

		graphPHOSYieldPi0PbPb1020 = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb1020, graphPHOSYieldPi0PbPb1020, fitModTsallisPi0PbPb1020,1,18);
		graphPHOSYieldPi0SysErrPbPb1020Red = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb1020, graphPHOSYieldPi0SysErrPbPb1020Red, fitModTsallisPi0PbPb1020,1,18);
		
		graphYieldPi0CombPbPb1020StatErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb1020, graphYieldPi0CombPbPb1020StatErr, fitModTsallisPi0PbPb1020, 0, graphYieldPi0CombPbPb1020StatErr->GetN());
		graphYieldPi0CombPbPb1020SysErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb1020, graphYieldPi0CombPbPb1020SysErr, fitModTsallisPi0PbPb1020, 0, graphYieldPi0CombPbPb1020SysErr->GetN());
		
		graphYieldPi0CombPbPb1020StatErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020StatErr->Clone("graphYieldPi0CombPbPb1020StatErrUnscaled");
		graphYieldPi0CombPbPb1020SysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020SysErr->Clone("graphYieldPi0CombPbPb1020SysErrUnscaled");

		
	//Shifting in Y

	TF1* fitModTsallisPi0PbPb1020YShift = ApplyYShift(graphYieldPi0CombPbPb1020Unshifted,&graphYieldPi0CombPbPb1020YShifted, "rad", "",0.3, parameters1020,0.00001,kTRUE);
	graphYieldPi0CombPbPb1020SysErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020SysErrUnshifted->Clone("YShiftedCombSys1020");
	graphYieldPi0CombPbPb1020SysErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb1020SysErrYShifted, fitModTsallisPi0PbPb1020YShift);
	graphYieldPi0CombPbPb1020StatErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020StatErrUnshifted->Clone("YShiftedCombStat1020");
	graphYieldPi0CombPbPb1020StatErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb1020StatErrYShifted, fitModTsallisPi0PbPb1020YShift);		

	graphYieldPi0PCMPbPb1020SysErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb1020Unshifted->Clone("YShiftedPCMSys1020");
	graphYieldPi0PCMPbPb1020SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb1020SysErrYShifted, fitModTsallisPi0PbPb1020YShift);
	graphYieldPi0PCMPbPb1020SysRAAErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrRAAPbPb1020Red->Clone("YShiftedPCMSysRAA1020");
	graphYieldPi0PCMPbPb1020SysRAAErrYShifted->RemovePoint(0);
	graphYieldPi0PCMPbPb1020SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb1020SysRAAErrYShifted, fitModTsallisPi0PbPb1020YShift);
	graphYieldPi0PCMPbPb1020StatErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb1020Unshifted->Clone("YShiftedPCMStat1020");
	graphYieldPi0PCMPbPb1020StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb1020StatErrYShifted, fitModTsallisPi0PbPb1020YShift);

	graphYieldPi0PHOSPbPb1020SysErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb1020Unshifted->Clone("YShiftedPHOSSys1020");
	graphYieldPi0PHOSPbPb1020SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb1020SysErrYShifted, fitModTsallisPi0PbPb1020YShift);
	graphYieldPi0PHOSPbPb1020SysRAAErrYShifted = (TGraphAsymmErrors*) graphSysErrRAAYieldPi0PHOSPbPb1020Red->Clone("YShiftedPHOSSysRAA1020");
	graphYieldPi0PHOSPbPb1020SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb1020SysRAAErrYShifted, fitModTsallisPi0PbPb1020YShift);
	graphYieldPi0PHOSPbPb1020StatErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb1020Unshifted->Clone("YShiftedPHOSStat1020");
	graphYieldPi0PHOSPbPb1020StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb1020StatErrYShifted, fitModTsallisPi0PbPb1020YShift);

	cout << endl << "***********************************************" << endl;
	cout << "Pi0 10-20% RAA*************************************" << endl;
	cout << "***********************************************" << endl<< endl;
	cout << "PCM*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	graphInvSectionCombPi02760GeV->Fit(fitInvCrossSectionPi0Comb2760GeV,"QNRMEX0+","",0.4,8.);
	cout << "Initial Fit result" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionPi0Comb2760GeV)<< endl;

	TGraphAsymmErrors* graphRAAPCM1020;
	TGraphAsymmErrors* graphRAASysPCM1020;
	CalcRaa( 	graphInvSectionPCMStatPi02760GeVRed, graphInvSectionPCMSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PCMPbPb1020StatErrYShifted, graphYieldPi0PCMPbPb1020SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPCM1020, &graphRAASysPCM1020,
					nColl1020, nCollErr1020,8.,0);

	cout << "PHOS*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAAPHOS1020;
	TGraphAsymmErrors* graphRAASysPHOS1020;
	CalcRaa( 	graphInvSectionPHOSStatPi02760GeVRed, graphInvSectionPHOSSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PHOSPbPb1020StatErrYShifted, graphYieldPi0PHOSPbPb1020SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPHOS1020, &graphRAASysPHOS1020,
					nColl1020, nCollErr1020, 8.,1);

	cout << "Combined*****************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAACombInd1020;
	TGraphAsymmErrors* graphRAASysCombInd1020;
	TGraphAsymmErrors* graphRAAPi0CombPbPb1020_2 = CombinePtPointsRAA(	graphRAAPCM1020, 		graphRAASysPCM1020,
												graphRAAPHOS1020, 	graphRAASysPHOS1020,
												graphRAACombInd1020, graphRAASysCombInd1020,
												xPtLimitsPbPbPHOS2, 19, 0, 0,2);

	TH1F *histoYieldPi0CombPbPb1020CombErrYShifted = new TH1F("histoYieldPi0CombPbPb1020CombErrYShifted", "", 19, xPtLimitsPbPbPHOS2) ; 
	Double_t* yvalues1020 = graphYieldPi0CombPbPb1020SysErrYShifted->GetY();
	graphYieldPi0CombPbPb1020SysErrYShifted->Print();
	for (Int_t i = 1; i < histoYieldPi0CombPbPb1020CombErrYShifted->GetNbinsX()+1; i++){
		histoYieldPi0CombPbPb1020CombErrYShifted->SetBinContent(i,yvalues1020[i-1]*2*TMath::Pi()* histoYieldPi0CombPbPb1020CombErrYShifted->GetBinCenter(i) );
		histoYieldPi0CombPbPb1020CombErrYShifted->SetBinError(i,TMath::Sqrt(pow(graphYieldPi0CombPbPb1020SysErrYShifted->GetErrorY(i-1),2)+pow(graphYieldPi0CombPbPb1020StatErrYShifted->GetErrorY(i-1),2)) *2*TMath::Pi()* histoYieldPi0CombPbPb1020CombErrYShifted->GetBinCenter(i) );
	}

	canvasDummy2->cd();
	canvasDummy2->SetLogy(0);
	histo2DDummy2->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphRAAPCM1020, 20,markerSizeCommonSpectrum1020, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAAPCM1020->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPCM1020 , 20,markerSizeCommonSpectrum1020, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAASysPCM1020->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAAPHOS1020, 20,markerSizeCommonSpectrum1020, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAAPHOS1020->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPHOS1020 , 20,markerSizeCommonSpectrum1020, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAASysPHOS1020->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd1020, 24,markerSizeCommonSpectrum1020,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAACombInd1020->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd1020 , 20,markerSizeCommonSpectrum1020, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAASysCombInd1020->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_RAA_1020.%s",outputDir.Data(),suffix.Data()));

	canvasDummy2->SetLogy();
	histo2DDummy3->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0PbPb1020, 20,markerSizeCommonSpectrum1020, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0PbPb1020->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0SysErrPbPb1020Red , 20,markerSizeCommonSpectrum1020, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0SysErrPbPb1020Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0PbPb1020, 20,markerSizeCommonSpectrum1020, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0PbPb1020->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrPbPb1020Red , 20,markerSizeCommonSpectrum1020, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0SysErrPbPb1020Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb1020StatErr, 24,markerSizeCommonSpectrum1020,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb1020StatErr->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb1020SysErr , 20,markerSizeCommonSpectrum1020, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb1020SysErr->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_Spectra_1020.%s",outputDir.Data(),suffix.Data()));


	// 	TGraphAsymmErrors* graphRAACombComb1020;
	// 	TGraphAsymmErrors* graphRAASysCombComb1020;
	// 	CalcRaa( 	graphInvSectionCombStatPi02760GeV, graphInvSectionCombSysPi02760GeV,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
	// 					graphYieldPi0CombPbPb1020StatErrYShifted, graphYieldPi0CombPbPb1020SysErrYShifted,  //PbPb Yields
	// 					&graphRAACombComb1020, &graphRAASysCombComb1020,
	// 					nColl1020, nCollErr1020);
	cout << "combined RAA with individual spectra calculated in this macro" << endl;
	graphRAACombInd1020->Print();
	// 	cout << "combined RAA based on combined spectra" << endl;
	// 	graphRAACombComb1020->Print();
	cout << WriteParameterToFile(fitModTsallisPi0PbPb1020YShift)<< endl;	
	forOutputToFile =  WriteParameterToFileLatexTable(fitModTsallisPi0PbPb1020YShift,1);
	fileFinalResultsFits << forOutputToFile << endl;

	graphPCMYieldPi0PbPb1020Unscaled = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb1020->Clone();	
	graphRatioCombPHOSPi0PbPb1020 = (TGraphAsymmErrors*) graphPHOSYieldPi0PbPb1020->Clone();	
	graphRatioCombPHOSPi0PbPb1020Sys = (TGraphAsymmErrors*) graphPHOSYieldPi0SysErrPbPb1020Red->Clone();	
	graphRatioCombConvPi0PbPb1020 = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb1020->Clone();	
	graphRatioCombConvPi0PbPb1020Sys = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrPbPb1020Red->Clone();	
	graphRatioCombPHOSPi0PbPb1020 = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb1020, fitYieldDataQCDPi0PbPb1020); 
	graphRatioCombPHOSPi0PbPb1020Sys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb1020Sys, fitYieldDataQCDPi0PbPb1020); 
	graphRatioCombConvPi0PbPb1020 = CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb1020, fitYieldDataQCDPi0PbPb1020); 
	graphRatioCombConvPi0PbPb1020Sys= CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb1020Sys, fitYieldDataQCDPi0PbPb1020); 
	graphRatioCombCombFitPbPb1020 = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020->Clone();
	graphRatioCombCombFitPbPb1020 = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb1020, fitYieldDataQCDPi0PbPb1020); 
	graphRatioCombCombFitPbPb1020Stat = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020StatErr->Clone();
	graphRatioCombCombFitPbPb1020Stat = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb1020Stat, fitYieldDataQCDPi0PbPb1020); 
	graphRatioCombCombFitPbPb1020Sys = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020SysErr->Clone();
	graphRatioCombCombFitPbPb1020Sys = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb1020Sys, fitYieldDataQCDPi0PbPb1020); 
		


	//******************* CombinePoints Pi0 20-40 % PbPb **********************
	cout << endl << "***********************************************" << endl;
	cout << "Pi0 20-40% Spectra*********************************" << endl;
	cout << "***********************************************" << endl<< endl ;
	graphYieldPi0CombPbPb2040 = CombinePtPointsSpectra(	histoYieldPbPb2040Red, 		graphPCMYieldPi0SysErrPbPb2040Red,
												histoPi0PHOSPbPb2040Red, 	graphPHOSYieldPi0SysErrPbPb2040Red,
												graphYieldPi0CombPbPb2040StatErr, graphYieldPi0CombPbPb2040SysErr,
												xPtLimitsPbPbPHOS, 20,   1, 0,2);
	graphYieldPi0CombPbPb2040Unshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040->Clone("graphYieldPi0CombPbPb2040Unshifted");
	graphYieldPi0CombPbPb2040StatErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040StatErr->Clone("graphYieldPi0CombPbPb2040StatErrUnshifted");
	graphYieldPi0CombPbPb2040SysErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040SysErr->Clone("graphYieldPi0CombPbPb2040SysErrUnshifted");

	Double_t parameters2040[15];
	ReturnParameterSetFittingPbPbFromString("2040",parameters2040);
	graphPCMYieldPi0PbPb2040 = new TGraphAsymmErrors(histoYieldPbPb2040Red);
	graphPCMYieldPi0PbPb2040->RemovePoint(0);
	graphPCMYieldPi0PbPb2040->RemovePoint(0);
	graphPCMYieldPi0SysErrPbPb2040Red->RemovePoint(0);
	graphPHOSYieldPi0PbPb2040 = new TGraphAsymmErrors(histoPi0PHOSPbPb2040Red);
	graphPCMYieldPi0PbPb2040Unshifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb2040->Clone("graphPCMYieldPi0PbPb2040Unshifted");
	graphPCMYieldPi0SysErrPbPb2040Unshifted= (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb2040Red->Clone("graphPCMYieldPi0SysErrPbPb2040Unshifted");
	graphPHOSYieldPi0PbPb2040Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb2040->Clone("graphPHOSYieldPi0PbPb2040Unshifted");
	graphPHOSYieldPi0SysErrPbPb2040Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb2040Red->Clone("graphPHOSYieldPi0SysErrPbPb2040Unshifted");

	//Shifting in X 
		fitModTsallisPi0PbPb2040 = FitObject("rad","ModTsallisPi0PbPb2040","Pi0");		
		SetParametersLimitsForFit (fitModTsallisPi0PbPb2040, 5, parameters2040);
		graphYieldPi0CombPbPb2040   = ApplyXshift(graphYieldPi0CombPbPb2040 ,fitModTsallisPi0PbPb2040);
		
		graphYieldPi0CombPbPb2040Unscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040->Clone("graphYieldPi0CombPbPb2040Unscaled");
		fitYieldPi0PbPb2040 = FitObject("rad","fitYieldPi0PbPb2040","Pi0",graphYieldPi0CombPbPb2040,0.6,maxPtPi0PbPb2040,parameters2040,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldPi0PbPb2040)<< endl;	
		forOutputToFile =  WriteParameterToFileLatexTable(fitYieldPi0PbPb2040,1);
		fileFinalResultsFits << forOutputToFile << endl;

		TF1* fitYieldDataQCDPi0PbPb2040 = FitObject("qcd","fitYieldDataQCDPi0PbPb2040","Pi0",graphYieldPi0CombPbPb2040,0.6,maxPtPi0PbPb2040,NULL,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldDataQCDPi0PbPb2040)<< endl;   
		forOutputToFile =  WriteParameterToFile(fitYieldDataQCDPi0PbPb2040);
		fileFinalResultsFits << forOutputToFile << endl;
		
		
		graphPCMYieldPi0PbPb2040= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb2040, graphPCMYieldPi0PbPb2040, fitModTsallisPi0PbPb2040, 0, 15);
		graphPCMYieldPi0SysErrPbPb2040Red= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb2040,graphPCMYieldPi0SysErrPbPb2040Red, fitModTsallisPi0PbPb2040,  0, 15);
		graphPHOSYieldPi0PbPb2040 = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb2040, graphPHOSYieldPi0PbPb2040, fitModTsallisPi0PbPb2040,1,18);
		graphPHOSYieldPi0SysErrPbPb2040Red = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb2040, graphPHOSYieldPi0SysErrPbPb2040Red, fitModTsallisPi0PbPb2040,1,18);
		
		graphYieldPi0CombPbPb2040StatErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb2040, graphYieldPi0CombPbPb2040StatErr, fitModTsallisPi0PbPb2040, 0, graphYieldPi0CombPbPb2040StatErr->GetN());
		graphYieldPi0CombPbPb2040SysErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb2040, graphYieldPi0CombPbPb2040SysErr, fitModTsallisPi0PbPb2040, 0, graphYieldPi0CombPbPb2040SysErr->GetN());
		
		graphYieldPi0CombPbPb2040StatErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040StatErr->Clone("graphYieldPi0CombPbPb2040StatErrUnscaled");
		graphYieldPi0CombPbPb2040SysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040SysErr->Clone("graphYieldPi0CombPbPb2040SysErrUnscaled");
		
	//Shifting in Y

	TF1* fitModTsallisPi0PbPb2040YShift = ApplyYShift(graphYieldPi0CombPbPb2040Unshifted,&graphYieldPi0CombPbPb2040YShifted, "rad", "",0.3, parameters2040,0.00001,kTRUE);
	graphYieldPi0CombPbPb2040SysErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040SysErrUnshifted->Clone("YShiftedCombSys2040");
	graphYieldPi0CombPbPb2040SysErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb2040SysErrYShifted, fitModTsallisPi0PbPb2040YShift);
	graphYieldPi0CombPbPb2040StatErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040StatErrUnshifted->Clone("YShiftedCombStat2040");
	graphYieldPi0CombPbPb2040StatErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb2040StatErrYShifted, fitModTsallisPi0PbPb2040YShift);		

	graphYieldPi0PCMPbPb2040SysErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb2040Unshifted->Clone("YShiftedPCMSys2040");
	graphYieldPi0PCMPbPb2040SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb2040SysErrYShifted, fitModTsallisPi0PbPb2040YShift);
	graphYieldPi0PCMPbPb2040SysRAAErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrRAAPbPb2040Red->Clone("YShiftedPCMSysRAA2040");
	graphYieldPi0PCMPbPb2040SysRAAErrYShifted->RemovePoint(0);
	graphYieldPi0PCMPbPb2040SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb2040SysRAAErrYShifted, fitModTsallisPi0PbPb2040YShift);
	graphYieldPi0PCMPbPb2040StatErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb2040Unshifted->Clone("YShiftedPCMStat2040");
	graphYieldPi0PCMPbPb2040StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb2040StatErrYShifted, fitModTsallisPi0PbPb2040YShift);

	graphYieldPi0PHOSPbPb2040SysErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb2040Unshifted->Clone("YShiftedPHOSSys2040");
	graphYieldPi0PHOSPbPb2040SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb2040SysErrYShifted, fitModTsallisPi0PbPb2040YShift);
	graphYieldPi0PHOSPbPb2040SysRAAErrYShifted = (TGraphAsymmErrors*) graphSysErrRAAYieldPi0PHOSPbPb2040Red->Clone("YShiftedPHOSSysRAA2040");
	graphYieldPi0PHOSPbPb2040SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb2040SysRAAErrYShifted, fitModTsallisPi0PbPb2040YShift);
	graphYieldPi0PHOSPbPb2040StatErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb2040Unshifted->Clone("YShiftedPHOSStat2040");
	graphYieldPi0PHOSPbPb2040StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb2040StatErrYShifted, fitModTsallisPi0PbPb2040YShift);

	cout << endl << "***********************************************" << endl;
	cout << "Pi0 20-40% RAA*************************************" << endl;
	cout << "***********************************************" << endl<< endl;
	cout << "PCM*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	graphInvSectionCombPi02760GeV->Fit(fitInvCrossSectionPi0Comb2760GeV,"QNRMEX0+","",0.4,8.);
	cout << "Initial Fit result" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionPi0Comb2760GeV)<< endl;

	TGraphAsymmErrors* graphRAAPCM2040;
	TGraphAsymmErrors* graphRAASysPCM2040;
	CalcRaa( 	graphInvSectionPCMStatPi02760GeVRed, graphInvSectionPCMSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PCMPbPb2040StatErrYShifted, graphYieldPi0PCMPbPb2040SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPCM2040, &graphRAASysPCM2040,
					nColl2040, nCollErr2040,8.,0);

	cout << "PHOS*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAAPHOS2040;
	TGraphAsymmErrors* graphRAASysPHOS2040;
	CalcRaa( 	graphInvSectionPHOSStatPi02760GeVRed, graphInvSectionPHOSSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PHOSPbPb2040StatErrYShifted, graphYieldPi0PHOSPbPb2040SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPHOS2040, &graphRAASysPHOS2040,
					nColl2040, nCollErr2040,8.,1);

	cout << "Combined*****************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAACombInd2040;
	TGraphAsymmErrors* graphRAASysCombInd2040;

	TGraphAsymmErrors* graphRAAPi0CombPbPb2040_2 = CombinePtPointsRAA(	graphRAAPCM2040, 		graphRAASysPCM2040,
												graphRAAPHOS2040, 	graphRAASysPHOS2040,
												graphRAACombInd2040, graphRAASysCombInd2040,
												xPtLimitsPbPbPHOS2, 19, 0,0,2);

	// 	TGraphAsymmErrors* graphRAACombComb2040;
	// 	TGraphAsymmErrors* graphRAASysCombComb2040;
	// 	CalcRaa( 	graphInvSectionCombStatPi02760GeV, graphInvSectionCombSysPi02760GeV,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
	// 					graphYieldPi0CombPbPb2040StatErrYShifted, graphYieldPi0CombPbPb2040SysErrYShifted,  //PbPb Yields
	// 					&graphRAACombComb2040, &graphRAASysCombComb2040,
	// 					nColl2040, nCollErr2040);

	TH1F *histoYieldPi0CombPbPb2040CombErrYShifted = new TH1F("histoYieldPi0CombPbPb2040CombErrYShifted", "", 19, xPtLimitsPbPbPHOS2) ; 
	Double_t* yvalues2040 = graphYieldPi0CombPbPb2040SysErrYShifted->GetY();
	graphYieldPi0CombPbPb2040SysErrYShifted->Print();
	for (Int_t i = 1; i < histoYieldPi0CombPbPb2040CombErrYShifted->GetNbinsX()+1; i++){
		histoYieldPi0CombPbPb2040CombErrYShifted->SetBinContent(i,yvalues2040[i-1]*2*TMath::Pi()* histoYieldPi0CombPbPb2040CombErrYShifted->GetBinCenter(i) );
		histoYieldPi0CombPbPb2040CombErrYShifted->SetBinError(i,TMath::Sqrt(pow(graphYieldPi0CombPbPb2040SysErrYShifted->GetErrorY(i-1),2)+pow(graphYieldPi0CombPbPb2040StatErrYShifted->GetErrorY(i-1),2)) *2*TMath::Pi()* histoYieldPi0CombPbPb2040CombErrYShifted->GetBinCenter(i) );
	}

	canvasDummy2->cd();
	canvasDummy2->SetLogy(0);
	histo2DDummy2->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphRAAPCM2040, 20,markerSizeCommonSpectrum2040, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAAPCM2040->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPCM2040 , 20,markerSizeCommonSpectrum2040, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAASysPCM2040->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAAPHOS2040, 20,markerSizeCommonSpectrum2040, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAAPHOS2040->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPHOS2040 , 20,markerSizeCommonSpectrum2040, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAASysPHOS2040->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd2040, 24,markerSizeCommonSpectrum2040,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAACombInd2040->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd2040 , 20,markerSizeCommonSpectrum2040, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAASysCombInd2040->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_RAA_2040.%s",outputDir.Data(),suffix.Data()));

	canvasDummy2->SetLogy();
	histo2DDummy3->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0PbPb2040, 20,markerSizeCommonSpectrum2040, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0PbPb2040->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0SysErrPbPb2040Red , 20,markerSizeCommonSpectrum2040, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0SysErrPbPb2040Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0PbPb2040, 20,markerSizeCommonSpectrum2040, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0PbPb2040->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrPbPb2040Red , 20,markerSizeCommonSpectrum2040, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0SysErrPbPb2040Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb2040StatErr, 24,markerSizeCommonSpectrum2040,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb2040StatErr->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb2040SysErr , 20,markerSizeCommonSpectrum2040, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb2040SysErr->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_Spectra_2040.%s",outputDir.Data(),suffix.Data()));


	cout << "combined RAA with individual spectra calculated in this macro" << endl;
	graphRAACombInd2040->Print();
	// 	cout << "combined RAA based on combined spectra" << endl;
	// 	graphRAACombComb2040->Print();
	cout << WriteParameterToFile(fitModTsallisPi0PbPb2040YShift)<< endl;	
	forOutputToFile =  WriteParameterToFileLatexTable(fitModTsallisPi0PbPb2040YShift,1);
	fileFinalResultsFits << forOutputToFile << endl;


	graphYieldPi0CombPbPb2040StatErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040StatErr->Clone("graphYieldPi0CombPbPb2040StatErrUnscaled");
	graphYieldPi0CombPbPb2040SysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040SysErr->Clone("graphYieldPi0CombPbPb2040SysErrUnscaled");

	graphPCMYieldPi0PbPb2040Unscaled = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb2040->Clone();	
	graphRatioCombPHOSPi0PbPb2040 = (TGraphAsymmErrors*) graphPHOSYieldPi0PbPb2040->Clone();	
	graphRatioCombPHOSPi0PbPb2040Sys = (TGraphAsymmErrors*) graphPHOSYieldPi0SysErrPbPb2040Red->Clone();	
	graphRatioCombConvPi0PbPb2040 = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb2040->Clone();	
	graphRatioCombConvPi0PbPb2040Sys = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrPbPb2040Red->Clone();	
	graphRatioCombPHOSPi0PbPb2040 = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb2040, fitYieldDataQCDPi0PbPb2040); 
	graphRatioCombPHOSPi0PbPb2040Sys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb2040Sys, fitYieldDataQCDPi0PbPb2040); 
	graphRatioCombConvPi0PbPb2040 = CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb2040, fitYieldDataQCDPi0PbPb2040); 
	graphRatioCombConvPi0PbPb2040Sys= CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb2040Sys, fitYieldDataQCDPi0PbPb2040); 
	graphRatioCombCombFitPbPb2040 = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040->Clone();
	graphRatioCombCombFitPbPb2040 = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb2040, fitYieldDataQCDPi0PbPb2040); 
	graphRatioCombCombFitPbPb2040Stat = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040StatErr->Clone();
	graphRatioCombCombFitPbPb2040Stat = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb2040Stat, fitYieldDataQCDPi0PbPb2040); 
	graphRatioCombCombFitPbPb2040Sys = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040SysErr->Clone();
	graphRatioCombCombFitPbPb2040Sys = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb2040Sys, fitYieldDataQCDPi0PbPb2040); 

	//******************* CombinePoints Pi0 40-60 % PbPb **********************
	cout << endl << "***********************************************" << endl;
	cout << "Pi0 40-60% Spectra*********************************" << endl;
	cout << "***********************************************" << endl<< endl ;
	graphYieldPi0CombPbPb4060 = CombinePtPointsSpectra(	histoYieldPbPb4060Red, 		graphPCMYieldPi0SysErrPbPb4060Red,
												histoPi0PHOSPbPb4060Red, 	graphPHOSYieldPi0SysErrPbPb4060Red,
												graphYieldPi0CombPbPb4060StatErr, graphYieldPi0CombPbPb4060SysErr,
												xPtLimitsPbPbPHOS, 20,   1, 0,2);
	graphYieldPi0CombPbPb4060Unshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060->Clone("graphYieldPi0CombPbPb4060Unshifted");
	graphYieldPi0CombPbPb4060StatErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060StatErr->Clone("graphYieldPi0CombPbPb4060StatErrUnshifted");
	graphYieldPi0CombPbPb4060SysErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060SysErr->Clone("graphYieldPi0CombPbPb4060SysErrUnshifted");

	Double_t parameters4060[15];
	ReturnParameterSetFittingPbPbFromString("4060",parameters4060);
	graphPCMYieldPi0PbPb4060 = new TGraphAsymmErrors(histoYieldPbPb4060Red);
	graphPCMYieldPi0PbPb4060->RemovePoint(0);
	graphPCMYieldPi0PbPb4060->RemovePoint(0);
	graphPCMYieldPi0SysErrPbPb4060Red->RemovePoint(0);
	graphPHOSYieldPi0PbPb4060 = new TGraphAsymmErrors(histoPi0PHOSPbPb4060Red);
	graphPCMYieldPi0PbPb4060Unshifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb4060->Clone("graphPCMYieldPi0PbPb4060Unshifted");
	graphPCMYieldPi0SysErrPbPb4060Unshifted= (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb4060Red->Clone("graphPCMYieldPi0SysErrPbPb4060Unshifted");
	graphPHOSYieldPi0PbPb4060Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb4060->Clone("graphPHOSYieldPi0PbPb4060Unshifted");
	graphPHOSYieldPi0SysErrPbPb4060Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb4060Red->Clone("graphPHOSYieldPi0SysErrPbPb4060Unshifted");

	//Shifting in X 
		fitModTsallisPi0PbPb4060 = FitObject("rad","ModTsallisPi0PbPb4060","Pi0");		
		SetParametersLimitsForFit (fitModTsallisPi0PbPb4060, 5, parameters4060);
		graphYieldPi0CombPbPb4060   = ApplyXshift(graphYieldPi0CombPbPb4060 ,fitModTsallisPi0PbPb4060);
		
		graphYieldPi0CombPbPb4060Unscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060->Clone("graphYieldPi0CombPbPb4060Unscaled");
		fitYieldPi0PbPb4060 = FitObject("rad","fitYieldPi0PbPb4060","Pi0",graphYieldPi0CombPbPb4060,0.6,maxPtPi0PbPb4060,parameters4060,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldPi0PbPb4060)<< endl;	
		forOutputToFile =  WriteParameterToFileLatexTable(fitYieldPi0PbPb4060,1);
		fileFinalResultsFits << forOutputToFile << endl;

		TF1* fitYieldDataQCDPi0PbPb4060 = FitObject("qcd","fitYieldDataQCDPi0PbPb4060","Pi0",graphYieldPi0CombPbPb4060,0.6,maxPtPi0PbPb4060,NULL,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldDataQCDPi0PbPb4060)<< endl;   
		forOutputToFile =  WriteParameterToFile(fitYieldDataQCDPi0PbPb4060);
		fileFinalResultsFits << forOutputToFile << endl;

				
		graphPCMYieldPi0PbPb4060= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb4060, graphPCMYieldPi0PbPb4060, fitModTsallisPi0PbPb4060, 0, 15);
		graphPCMYieldPi0SysErrPbPb4060Red= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb4060,graphPCMYieldPi0SysErrPbPb4060Red, fitModTsallisPi0PbPb4060,  0, 15);
		graphPHOSYieldPi0PbPb4060 = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb4060, graphPHOSYieldPi0PbPb4060, fitModTsallisPi0PbPb4060,1,18);
		graphPHOSYieldPi0SysErrPbPb4060Red = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb4060, graphPHOSYieldPi0SysErrPbPb4060Red, fitModTsallisPi0PbPb4060,1,18);
		
		graphYieldPi0CombPbPb4060StatErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb4060, graphYieldPi0CombPbPb4060StatErr, fitModTsallisPi0PbPb4060, 0, graphYieldPi0CombPbPb4060StatErr->GetN());
		graphYieldPi0CombPbPb4060SysErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb4060, graphYieldPi0CombPbPb4060SysErr, fitModTsallisPi0PbPb4060, 0, graphYieldPi0CombPbPb4060SysErr->GetN());
		
		graphYieldPi0CombPbPb4060StatErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060StatErr->Clone("graphYieldPi0CombPbPb4060StatErrUnscaled");
		graphYieldPi0CombPbPb4060SysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060SysErr->Clone("graphYieldPi0CombPbPb4060SysErrUnscaled");
		
	//Shifting in Y

	TF1* fitModTsallisPi0PbPb4060YShift = ApplyYShift(graphYieldPi0CombPbPb4060Unshifted,&graphYieldPi0CombPbPb4060YShifted, "rad", "",0.3, parameters4060,0.00001,kTRUE);
	graphYieldPi0CombPbPb4060SysErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060SysErrUnshifted->Clone("YShiftedCombSys4060");
	graphYieldPi0CombPbPb4060SysErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb4060SysErrYShifted, fitModTsallisPi0PbPb4060YShift);
	graphYieldPi0CombPbPb4060StatErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060StatErrUnshifted->Clone("YShiftedCombStat4060");
	graphYieldPi0CombPbPb4060StatErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb4060StatErrYShifted, fitModTsallisPi0PbPb4060YShift);		

	graphYieldPi0PCMPbPb4060SysErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb4060Unshifted->Clone("YShiftedPCMSys4060");
	graphYieldPi0PCMPbPb4060SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb4060SysErrYShifted, fitModTsallisPi0PbPb4060YShift);
	graphYieldPi0PCMPbPb4060SysRAAErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrRAAPbPb4060Red->Clone("YShiftedPCMSysRAA4060");
	graphYieldPi0PCMPbPb4060SysRAAErrYShifted->RemovePoint(0);
	graphYieldPi0PCMPbPb4060SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb4060SysRAAErrYShifted, fitModTsallisPi0PbPb4060YShift);
	graphYieldPi0PCMPbPb4060StatErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb4060Unshifted->Clone("YShiftedPCMStat4060");
	graphYieldPi0PCMPbPb4060StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb4060StatErrYShifted, fitModTsallisPi0PbPb4060YShift);

	graphYieldPi0PHOSPbPb4060SysErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb4060Unshifted->Clone("YShiftedPHOSSys4060");
	graphYieldPi0PHOSPbPb4060SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb4060SysErrYShifted, fitModTsallisPi0PbPb4060YShift);
	graphYieldPi0PHOSPbPb4060SysRAAErrYShifted = (TGraphAsymmErrors*) graphSysErrRAAYieldPi0PHOSPbPb4060Red->Clone("YShiftedPHOSSysRAA4060");
	graphYieldPi0PHOSPbPb4060SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb4060SysRAAErrYShifted, fitModTsallisPi0PbPb4060YShift);
	graphYieldPi0PHOSPbPb4060StatErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb4060Unshifted->Clone("YShiftedPHOSStat4060");
	graphYieldPi0PHOSPbPb4060StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb4060StatErrYShifted, fitModTsallisPi0PbPb4060YShift);

	cout << endl << "***********************************************" << endl;
	cout << "Pi0 40-60% RAA*************************************" << endl;
	cout << "***********************************************" << endl<< endl;
	cout << "PCM*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	graphInvSectionCombPi02760GeV->Fit(fitInvCrossSectionPi0Comb2760GeV,"QNRMEX0+","",0.4,8.);
	cout << "Initial Fit result" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionPi0Comb2760GeV)<< endl;

	TGraphAsymmErrors* graphRAAPCM4060;
	TGraphAsymmErrors* graphRAASysPCM4060;
	CalcRaa( 	graphInvSectionPCMStatPi02760GeVRed, graphInvSectionPCMSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PCMPbPb4060StatErrYShifted, graphYieldPi0PCMPbPb4060SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPCM4060, &graphRAASysPCM4060,
					nColl4060, nCollErr4060,8.,0);

	cout << "PHOS*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAAPHOS4060;
	TGraphAsymmErrors* graphRAASysPHOS4060;
	CalcRaa( 	graphInvSectionPHOSStatPi02760GeVRed, graphInvSectionPHOSSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PHOSPbPb4060StatErrYShifted, graphYieldPi0PHOSPbPb4060SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPHOS4060, &graphRAASysPHOS4060,
					nColl4060, nCollErr4060,8.,1);

	cout << "Combined*****************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAACombInd4060;
	TGraphAsymmErrors* graphRAASysCombInd4060;

	TGraphAsymmErrors* graphRAAPi0CombPbPb4060_2 = CombinePtPointsRAA(	graphRAAPCM4060, 		graphRAASysPCM4060,
												graphRAAPHOS4060, 	graphRAASysPHOS4060,
												graphRAACombInd4060, graphRAASysCombInd4060,
												xPtLimitsPbPbPHOS2, 19, 0,0,2);

	TH1F *histoYieldPi0CombPbPb4060CombErrYShifted = new TH1F("histoYieldPi0CombPbPb4060CombErrYShifted", "", 19, xPtLimitsPbPbPHOS2) ; 
	Double_t* yvalues4060 = graphYieldPi0CombPbPb4060SysErrYShifted->GetY();
	graphYieldPi0CombPbPb4060SysErrYShifted->Print();
	for (Int_t i = 1; i < histoYieldPi0CombPbPb4060CombErrYShifted->GetNbinsX()+1; i++){
		histoYieldPi0CombPbPb4060CombErrYShifted->SetBinContent(i,yvalues4060[i-1]*2*TMath::Pi()* histoYieldPi0CombPbPb4060CombErrYShifted->GetBinCenter(i) );
		histoYieldPi0CombPbPb4060CombErrYShifted->SetBinError(i,TMath::Sqrt(pow(graphYieldPi0CombPbPb4060SysErrYShifted->GetErrorY(i-1),2)+pow(graphYieldPi0CombPbPb4060StatErrYShifted->GetErrorY(i-1),2)) *2*TMath::Pi()* histoYieldPi0CombPbPb4060CombErrYShifted->GetBinCenter(i) );
	}


	canvasDummy2->cd();
	canvasDummy2->SetLogy(0);
	histo2DDummy2->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphRAAPCM4060, 20,markerSizeCommonSpectrum4060, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAAPCM4060->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPCM4060 , 20,markerSizeCommonSpectrum4060, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAASysPCM4060->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAAPHOS4060, 20,markerSizeCommonSpectrum4060, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAAPHOS4060->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPHOS4060 , 20,markerSizeCommonSpectrum4060, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAASysPHOS4060->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd4060, 24,markerSizeCommonSpectrum4060,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAACombInd4060->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd4060 , 20,markerSizeCommonSpectrum4060, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAASysCombInd4060->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_RAA_4060.%s",outputDir.Data(),suffix.Data()));

	canvasDummy2->SetLogy();
	histo2DDummy3->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0PbPb4060, 20,markerSizeCommonSpectrum4060, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0PbPb4060->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrPbPb4060Red , 20,markerSizeCommonSpectrum4060, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0SysErrPbPb4060Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0PbPb4060, 20,markerSizeCommonSpectrum4060, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0PbPb4060->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0SysErrPbPb4060Red , 20,markerSizeCommonSpectrum4060, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0SysErrPbPb4060Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb4060StatErr, 24,markerSizeCommonSpectrum4060,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb4060StatErr->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb4060SysErr , 20,markerSizeCommonSpectrum4060, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb4060SysErr->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_Spectra_4060.%s",outputDir.Data(),suffix.Data()));

	cout << "comparison Graph 40-60% PCM XXXXXX" << endl,
	graphPCMYieldPi0PbPb4060->Print();


	cout << "comparison Graph 40-60% PHOS XXXXXX" << endl,
	graphPHOSYieldPi0PbPb4060->Print();


	cout << WriteParameterToFile(fitModTsallisPi0PbPb4060YShift)<< endl;	
	forOutputToFile =  WriteParameterToFileLatexTable(fitModTsallisPi0PbPb4060YShift,1);
	fileFinalResultsFits << forOutputToFile << endl;


	graphPCMYieldPi0PbPb4060Unscaled = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb4060->Clone();	
	graphRatioCombPHOSPi0PbPb4060 = (TGraphAsymmErrors*) graphPHOSYieldPi0PbPb4060->Clone();	
	graphRatioCombPHOSPi0PbPb4060Sys = (TGraphAsymmErrors*) graphPHOSYieldPi0SysErrPbPb4060Red->Clone();	
	graphRatioCombConvPi0PbPb4060 = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb4060->Clone();	
	graphRatioCombConvPi0PbPb4060Sys = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrPbPb4060Red->Clone();	
	graphRatioCombCombFitPbPb4060 = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060->Clone();
	graphRatioCombCombFitPbPb4060Stat = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060StatErr->Clone();
	graphRatioCombCombFitPbPb4060Sys = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060SysErr->Clone();

	// 	graphRatioCombPHOSPi0PbPb4060 = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb4060, fitYieldPi0PbPb4060); 
	// 	graphRatioCombPHOSPi0PbPb4060Sys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb4060Sys, fitYieldPi0PbPb4060); 
	// 	graphRatioCombConvPi0PbPb4060 = CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb4060, fitYieldPi0PbPb4060); 
	// 	graphRatioCombConvPi0PbPb4060Sys= CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb4060Sys, fitYieldPi0PbPb4060); 
	// 	graphRatioCombCombFitPbPb4060 = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb4060, fitYieldPi0PbPb4060); 
	// 	graphRatioCombCombFitPbPb4060Stat = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb4060Stat, fitYieldPi0PbPb4060); 
	// 	graphRatioCombCombFitPbPb4060Sys = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb4060Sys, fitYieldPi0PbPb4060); 

	graphRatioCombPHOSPi0PbPb4060 = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb4060, fitYieldDataQCDPi0PbPb4060); 
	graphRatioCombPHOSPi0PbPb4060Sys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb4060Sys, fitYieldDataQCDPi0PbPb4060); 
	graphRatioCombConvPi0PbPb4060 = CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb4060, fitYieldDataQCDPi0PbPb4060); 
	graphRatioCombConvPi0PbPb4060Sys= CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb4060Sys, fitYieldDataQCDPi0PbPb4060); 
	graphRatioCombCombFitPbPb4060 = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb4060, fitYieldDataQCDPi0PbPb4060); 
	graphRatioCombCombFitPbPb4060Stat = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb4060Stat, fitYieldDataQCDPi0PbPb4060); 
	graphRatioCombCombFitPbPb4060Sys = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb4060Sys, fitYieldDataQCDPi0PbPb4060); 

	//******************* CombinePoints Pi0 60-80 % PbPb **********************
	cout << endl << "***********************************************" << endl;
	cout << "Pi0 60-80% Spectra*********************************" << endl;
	cout << "***********************************************" << endl<< endl ;
	graphYieldPi0CombPbPb6080 = CombinePtPointsSpectra(	histoYieldPbPb6080Red, 		graphPCMYieldPi0SysErrPbPb6080Red,
												histoPi0PHOSPbPb6080Red, 	graphPHOSYieldPi0SysErrPbPb6080Red,
												graphYieldPi0CombPbPb6080StatErr, graphYieldPi0CombPbPb6080SysErr,
												xPtLimitsPbPbPHOS, 19,   1, 0,2);
	graphYieldPi0CombPbPb6080Unshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080->Clone("graphYieldPi0CombPbPb6080Unshifted");
	graphYieldPi0CombPbPb6080StatErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080StatErr->Clone("graphYieldPi0CombPbPb6080StatErrUnshifted");
	graphYieldPi0CombPbPb6080SysErrUnshifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080SysErr->Clone("graphYieldPi0CombPbPb6080SysErrUnshifted");

	Double_t parameters6080[15];
	ReturnParameterSetFittingPbPbFromString("6080",parameters6080);
	graphPCMYieldPi0PbPb6080 = new TGraphAsymmErrors(histoYieldPbPb6080Red);
	graphPCMYieldPi0PbPb6080->RemovePoint(0);
	graphPCMYieldPi0PbPb6080->RemovePoint(0);
	graphPCMYieldPi0SysErrPbPb6080Red->RemovePoint(0);
	graphPHOSYieldPi0PbPb6080 = new TGraphAsymmErrors(histoPi0PHOSPbPb6080Red);
	graphPCMYieldPi0PbPb6080Unshifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb6080->Clone("graphPCMYieldPi0PbPb6080Unshifted");
	graphPCMYieldPi0SysErrPbPb6080Unshifted= (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb6080Red->Clone("graphPCMYieldPi0SysErrPbPb6080Unshifted");
	graphPHOSYieldPi0PbPb6080Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb6080->Clone("graphPHOSYieldPi0PbPb6080Unshifted");
	graphPHOSYieldPi0SysErrPbPb6080Unshifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb6080Red->Clone("graphPHOSYieldPi0SysErrPbPb6080Unshifted");

	//Shifting in X 
		fitModTsallisPi0PbPb6080 = FitObject("l","fitModTsallisPi0PbPb6080","Pi0"); //fitModTsallisPi0PbPb6080 FitObject("rad","ModTsallisPi0PbPb6080","Pi0");		
	// 		SetParametersLimitsForFit (fitModTsallisPi0PbPb6080, 5, parameters6080);
		graphYieldPi0CombPbPb6080   = ApplyXshift(graphYieldPi0CombPbPb6080 ,fitModTsallisPi0PbPb6080);
		
		graphYieldPi0CombPbPb6080Unscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080->Clone("graphYieldPi0CombPbPb6080Unscaled");
		fitYieldPi0PbPb6080 = FitObject("l","fitYieldPi0PbPb6080","Pi0",graphYieldPi0CombPbPb6080,0.6,maxPtPi0PbPb6080,NULL,"QNRMEX0+");
	// 		FitObject("rad","fitYieldPi0PbPb6080","Pi0",graphYieldPi0CombPbPb6080,0.6,maxPtPi0PbPb6080,parameters6080);
		cout << WriteParameterToFile(fitYieldPi0PbPb6080)<< endl;	
		forOutputToFile =  WriteParameterToFileLatexTable(fitYieldPi0PbPb6080,1);
		fileFinalResultsFits << forOutputToFile << endl;

		TF1* fitYieldDataQCDPi0PbPb6080 = FitObject("l","fitYieldDataQCDPi0PbPb6080","Pi0",graphYieldPi0CombPbPb6080,0.6,maxPtPi0PbPb6080,NULL,"QNRMEX0+");
		cout << WriteParameterToFile(fitYieldDataQCDPi0PbPb6080)<< endl;   
		forOutputToFile =  WriteParameterToFile(fitYieldDataQCDPi0PbPb6080);
		fileFinalResultsFits << forOutputToFile << endl;
		
				
		graphPCMYieldPi0PbPb6080= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb6080, graphPCMYieldPi0PbPb6080, fitModTsallisPi0PbPb6080, 0, 15);
		graphPCMYieldPi0SysErrPbPb6080Red= ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb6080,graphPCMYieldPi0SysErrPbPb6080Red, fitModTsallisPi0PbPb6080,  0, 15);
		graphPHOSYieldPi0PbPb6080 = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb6080, graphPHOSYieldPi0PbPb6080, fitModTsallisPi0PbPb6080,1,18);
		graphPHOSYieldPi0SysErrPbPb6080Red = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb6080, graphPHOSYieldPi0SysErrPbPb6080Red, fitModTsallisPi0PbPb6080,1,18);
		
		graphYieldPi0CombPbPb6080StatErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb6080, graphYieldPi0CombPbPb6080StatErr, fitModTsallisPi0PbPb6080, 0, graphYieldPi0CombPbPb6080StatErr->GetN());
		graphYieldPi0CombPbPb6080SysErr = ApplyXshiftIndividualSpectra(graphYieldPi0CombPbPb6080, graphYieldPi0CombPbPb6080SysErr, fitModTsallisPi0PbPb6080, 0, graphYieldPi0CombPbPb6080SysErr->GetN());
		
		graphYieldPi0CombPbPb6080StatErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080StatErr->Clone("graphYieldPi0CombPbPb6080StatErrUnscaled");
		graphYieldPi0CombPbPb6080SysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080SysErr->Clone("graphYieldPi0CombPbPb6080SysErrUnscaled");
		
	//Shifting in Y

	TF1* fitModTsallisPi0PbPb6080YShift = ApplyYShift(graphYieldPi0CombPbPb6080Unshifted,&graphYieldPi0CombPbPb6080YShifted,"l", "",0.3);
																		//"rad", "",0.3, parameters6080,0.00001,kTRUE);
	graphYieldPi0CombPbPb6080SysErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080SysErrUnshifted->Clone("YShiftedCombSys6080");
	graphYieldPi0CombPbPb6080SysErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb6080SysErrYShifted, fitModTsallisPi0PbPb6080YShift);
	graphYieldPi0CombPbPb6080StatErrYShifted = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080StatErrUnshifted->Clone("YShiftedCombStat6080");
	graphYieldPi0CombPbPb6080StatErrYShifted= ApplyYshiftIndividualSpectra( graphYieldPi0CombPbPb6080StatErrYShifted, fitModTsallisPi0PbPb6080YShift);		

	graphYieldPi0PCMPbPb6080SysErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrPbPb6080Unshifted->Clone("YShiftedPCMSys6080");
	graphYieldPi0PCMPbPb6080SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb6080SysErrYShifted, fitModTsallisPi0PbPb6080YShift);
	graphYieldPi0PCMPbPb6080SysRAAErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrRAAPbPb6080Red->Clone("YShiftedPCMSysRAA6080");
	graphYieldPi0PCMPbPb6080SysRAAErrYShifted->RemovePoint(0);
	graphYieldPi0PCMPbPb6080SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb6080SysRAAErrYShifted, fitModTsallisPi0PbPb6080YShift);
	graphYieldPi0PCMPbPb6080StatErrYShifted = (TGraphAsymmErrors*)graphPCMYieldPi0PbPb6080Unshifted->Clone("YShiftedPCMStat6080");
	graphYieldPi0PCMPbPb6080StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PCMPbPb6080StatErrYShifted, fitModTsallisPi0PbPb6080YShift);

	graphYieldPi0PHOSPbPb6080SysErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0SysErrPbPb6080Unshifted->Clone("YShiftedPHOSSys6080");
	graphYieldPi0PHOSPbPb6080SysErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb6080SysErrYShifted, fitModTsallisPi0PbPb6080YShift);
	graphYieldPi0PHOSPbPb6080SysRAAErrYShifted = (TGraphAsymmErrors*) graphSysErrRAAYieldPi0PHOSPbPb6080Red->Clone("YShiftedPHOSSysRAA6080");
	graphYieldPi0PHOSPbPb6080SysRAAErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb6080SysRAAErrYShifted, fitModTsallisPi0PbPb6080YShift);
	graphYieldPi0PHOSPbPb6080StatErrYShifted = (TGraphAsymmErrors*)graphPHOSYieldPi0PbPb6080Unshifted->Clone("YShiftedPHOSStat6080");
	graphYieldPi0PHOSPbPb6080StatErrYShifted = ApplyYshiftIndividualSpectra( graphYieldPi0PHOSPbPb6080StatErrYShifted, fitModTsallisPi0PbPb6080YShift);

	cout << endl << "***********************************************" << endl;
	cout << "Pi0 60-80% RAA*************************************" << endl;
	cout << "***********************************************" << endl<< endl;
	cout << "PCM*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	graphInvSectionCombPi02760GeV->Fit(fitInvCrossSectionPi0Comb2760GeV,"QNRMEX0+","",0.4,8.);
	cout << "Initial Fit result" << endl;
	cout << WriteParameterToFile(fitInvCrossSectionPi0Comb2760GeV)<< endl;

	TGraphAsymmErrors* graphRAAPCM6080;
	TGraphAsymmErrors* graphRAASysPCM6080;
	CalcRaa( 	graphInvSectionPCMStatPi02760GeVRed, graphInvSectionPCMSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PCMPbPb6080StatErrYShifted, graphYieldPi0PCMPbPb6080SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPCM6080, &graphRAASysPCM6080,
					nColl6080, nCollErr6080,8.,0);

	cout << "PHOS*********************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAAPHOS6080;
	TGraphAsymmErrors* graphRAASysPHOS6080;
	CalcRaa( 	graphInvSectionPHOSStatPi02760GeVRed, graphInvSectionPHOSSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
					graphYieldPi0PHOSPbPb6080StatErrYShifted, graphYieldPi0PHOSPbPb6080SysRAAErrYShifted,  //PbPb Yields
					&graphRAAPHOS6080, &graphRAASysPHOS6080,
					nColl6080, nCollErr6080,8.,1);

	cout << "Combined*****************************************" << endl ;
	cout << "***********************************************" << endl<< endl ;

	TGraphAsymmErrors* graphRAACombInd6080;
	TGraphAsymmErrors* graphRAASysCombInd6080;

	TGraphAsymmErrors* graphRAAPi0CombPbPb6080_2 = CombinePtPointsRAA(	graphRAAPCM6080, 		graphRAASysPCM6080,
												graphRAAPHOS6080, 	graphRAASysPHOS6080,
												graphRAACombInd6080, graphRAASysCombInd6080,
												xPtLimitsPbPbPHOS2, 18, 0,0,2);

	TH1F *histoYieldPi0CombPbPb6080CombErrYShifted = new TH1F("histoYieldPi0CombPbPb6080CombErrYShifted", "", 19, xPtLimitsPbPbPHOS2) ; 
	Double_t* yvalues6080 = graphYieldPi0CombPbPb6080SysErrYShifted->GetY();
	graphYieldPi0CombPbPb6080SysErrYShifted->Print();
	for (Int_t i = 1; i < histoYieldPi0CombPbPb6080CombErrYShifted->GetNbinsX()+1; i++){
		histoYieldPi0CombPbPb6080CombErrYShifted->SetBinContent(i,yvalues6080[i-1]*2*TMath::Pi()* histoYieldPi0CombPbPb6080CombErrYShifted->GetBinCenter(i) );
		histoYieldPi0CombPbPb6080CombErrYShifted->SetBinError(i,TMath::Sqrt(pow(graphYieldPi0CombPbPb6080SysErrYShifted->GetErrorY(i-1),2)+pow(graphYieldPi0CombPbPb6080StatErrYShifted->GetErrorY(i-1),2)) *2*TMath::Pi()* histoYieldPi0CombPbPb6080CombErrYShifted->GetBinCenter(i) );
	}

	canvasDummy2->cd();
	canvasDummy2->SetLogy(0);
	histo2DDummy2->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphRAAPCM6080, 20,markerSizeCommonSpectrum6080, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAAPCM6080->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPCM6080 , 20,markerSizeCommonSpectrum6080, kRed, kRed, widthLinesBoxes, kTRUE);
	graphRAASysPCM6080->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAAPHOS6080, 20,markerSizeCommonSpectrum6080, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAAPHOS6080->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysPHOS6080 , 20,markerSizeCommonSpectrum6080, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphRAASysPHOS6080->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd6080, 24,markerSizeCommonSpectrum6080,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAACombInd6080->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd6080 , 20,markerSizeCommonSpectrum6080, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphRAASysCombInd6080->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_RAA_6080.%s",outputDir.Data(),suffix.Data()));

	canvasDummy2->SetLogy();
	histo2DDummy3->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0PbPb6080, 20,markerSizeCommonSpectrum6080, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0PbPb6080->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0SysErrPbPb6080Red , 20,markerSizeCommonSpectrum6080, kRed, kRed, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0SysErrPbPb6080Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0PbPb6080, 20,markerSizeCommonSpectrum6080, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPHOSYieldPi0PbPb6080->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0SysErrPbPb6080Red , 20,markerSizeCommonSpectrum6080, kBlue, kBlue, widthLinesBoxes, kTRUE);
	graphPCMYieldPi0SysErrPbPb6080Red->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb6080StatErr, 24,markerSizeCommonSpectrum6080,  kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb6080StatErr->Draw("pEsame");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb6080SysErr , 20,markerSizeCommonSpectrum6080, kBlack , kBlack, widthLinesBoxes, kTRUE);
	graphYieldPi0CombPbPb6080SysErr->Draw("E2same");
		
	canvasDummy2->Update();
	canvasDummy2->Print(Form("%s/Comparison_Spectra_6080.%s",outputDir.Data(),suffix.Data()));


	// 	TGraphAsymmErrors* graphRAACombComb6080;
	// 	TGraphAsymmErrors* graphRAASysCombComb6080;
	// 	CalcRaa( 	graphInvSectionCombStatPi02760GeV, graphInvSectionCombSysPi02760GeV,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
	// 					graphYieldPi0CombPbPb6080StatErrYShifted, graphYieldPi0CombPbPb6080SysErrYShifted,  //PbPb Yields
	// 					&graphRAACombComb6080, &graphRAASysCombComb6080,
	// 					nColl6080, nCollErr6080);

	cout << "combined RAA with individual spectra calculated in this macro" << endl;
	graphRAACombInd6080->Print();
	// 	cout << "combined RAA based on combined spectra" << endl;
	// 	graphRAACombComb6080->Print();
	cout << WriteParameterToFile(fitModTsallisPi0PbPb6080YShift)<< endl;	
	forOutputToFile =  WriteParameterToFileLatexTable(fitModTsallisPi0PbPb6080YShift,1);
	fileFinalResultsFits << forOutputToFile << endl;


	graphPCMYieldPi0PbPb6080Unscaled = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb6080->Clone();	
	graphRatioCombPHOSPi0PbPb6080 = (TGraphAsymmErrors*) graphPHOSYieldPi0PbPb6080->Clone();	
	graphRatioCombPHOSPi0PbPb6080Sys = (TGraphAsymmErrors*) graphPHOSYieldPi0SysErrPbPb6080Red->Clone();	
	graphRatioCombConvPi0PbPb6080 = (TGraphAsymmErrors*) graphPCMYieldPi0PbPb6080->Clone();	
	graphRatioCombConvPi0PbPb6080Sys = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrPbPb6080Red->Clone();	
	graphRatioCombCombFitPbPb6080 = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080->Clone();
	graphRatioCombCombFitPbPb6080Stat = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080StatErr->Clone();
	graphRatioCombCombFitPbPb6080Sys = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080SysErr->Clone();
	// 	graphRatioCombPHOSPi0PbPb6080 = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb6080, fitYieldPi0PbPb6080); 
	// 	graphRatioCombPHOSPi0PbPb6080Sys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb6080Sys, fitYieldPi0PbPb6080); 
	// 	graphRatioCombConvPi0PbPb6080 = CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb6080, fitYieldPi0PbPb6080); 
	// 	graphRatioCombConvPi0PbPb6080Sys= CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb6080Sys, fitYieldPi0PbPb6080); 
	// 	graphRatioCombCombFitPbPb6080 = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb6080, fitYieldPi0PbPb6080); 
	// 	graphRatioCombCombFitPbPb6080Stat = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb6080Stat, fitYieldPi0PbPb6080); 
	// 	graphRatioCombCombFitPbPb6080Sys = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb6080Sys, fitYieldPi0PbPb6080); 
	graphRatioCombPHOSPi0PbPb6080 = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb6080, fitYieldDataQCDPi0PbPb6080); 
	graphRatioCombPHOSPi0PbPb6080Sys = CalculateGraphErrRatioToFit (graphRatioCombPHOSPi0PbPb6080Sys, fitYieldDataQCDPi0PbPb6080); 
	graphRatioCombConvPi0PbPb6080 = CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb6080, fitYieldDataQCDPi0PbPb6080); 
	graphRatioCombConvPi0PbPb6080Sys= CalculateGraphErrRatioToFit (graphRatioCombConvPi0PbPb6080Sys, fitYieldDataQCDPi0PbPb6080); 
	graphRatioCombCombFitPbPb6080 = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb6080, fitYieldDataQCDPi0PbPb6080); 
	graphRatioCombCombFitPbPb6080Stat = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb6080Stat, fitYieldDataQCDPi0PbPb6080); 
	graphRatioCombCombFitPbPb6080Sys = CalculateGraphErrRatioToFit(graphRatioCombCombFitPbPb6080Sys, fitYieldDataQCDPi0PbPb6080); 


	fitInvCrossSectionPi0Comb2760GeV = FitObject("l","fitInvCrossSectionPi0Comb2760GeV","Pi0",graphInvSectionCombPi02760GeV,0.6,maxPtPi0PbPb6080,NULL,"QNRMEX0+");
	cout << WriteParameterToFile(fitInvCrossSectionPi0Comb2760GeV)<< endl;
	forOutputToFile =  WriteParameterToFileLatexTable(fitInvCrossSectionPi0Comb2760GeV);
	fileFinalResultsFits << forOutputToFile << endl;

	TF1* fitInvCrossSectionPi0Comb2760GeVPow = FitObject("powPure","fitInvCrossSectionPi0Comb2760GeVPow","Pi0",graphInvSectionCombPi02760GeV,3.,maxPtPi0PbPb6080,NULL,"QNRMEX0+");
	cout << WriteParameterToFile(fitInvCrossSectionPi0Comb2760GeVPow)<< endl;	
	forOutputToFile =  WriteParameterToFileLatexTable(fitInvCrossSectionPi0Comb2760GeVPow);
	fileFinalResultsFits << forOutputToFile << endl;


	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0005SysErr, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005, widthLinesBoxes, kTRUE);//, colorComb0005-5);
	graphYieldPi0CombPbPb0005SysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0510SysErr, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510, widthLinesBoxes, kTRUE);//, colorComb0510-5);
	graphYieldPi0CombPbPb0510SysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb1020SysErr, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020, widthLinesBoxes, kTRUE);//, colorComb1020-5);
	graphYieldPi0CombPbPb1020SysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb2040SysErr, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);//, colorComb2040-5);
	graphYieldPi0CombPbPb2040SysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb4060SysErr, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060, widthLinesBoxes, kTRUE);//, colorComb4060-5);
	graphYieldPi0CombPbPb4060SysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb6080SysErr, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080, widthLinesBoxes, kTRUE);//, colorComb6080-5);
	graphYieldPi0CombPbPb6080SysErr->Draw("E2same");


	TGraphAsymmErrors* graphYieldPi0CombPbPb0005FullSysErr = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005SysErr->Clone("graphYieldPi0CombPbPb0005FullSysErr");
	Double_t* Yield0005 = graphYieldPi0CombPbPb0005FullSysErr->GetY();
	Double_t* SysErrorHigh0005 = graphYieldPi0CombPbPb0005FullSysErr->GetEYhigh();
	Double_t* SysErrorLow0005 = graphYieldPi0CombPbPb0005FullSysErr->GetEYlow();
	for (Int_t i = 0; i < graphYieldPi0CombPbPb0005FullSysErr->GetN(); i++){
		Double_t relativeErr = SysErrorHigh0005[i]/Yield0005[i]*100;
		Double_t relativeErrNew = TMath::Sqrt(TMath::Power(SysErrorHigh0005[i]/Yield0005[i]*100,2)+  commonCentralityErr0005*commonCentralityErr0005);
		Double_t absoluteErr = SysErrorHigh0005[i];
		Double_t absoluteErrNew = relativeErrNew/100*Yield0005[i];
		cout << i << "\t" <<relativeErr << "\t" << relativeErrNew << "\t" << absoluteErr << "\t" << absoluteErrNew <<endl;
		SysErrorHigh0005[i] = absoluteErrNew;
		SysErrorLow0005[i] = absoluteErrNew;
	}
	TGraphAsymmErrors* graphYieldPi0CombPbPb0005FullSysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0005FullSysErr->Clone("graphYieldPi0CombPbPb0005FullSysErrUnscaled");

	TGraphAsymmErrors* graphYieldPi0CombPbPb0510FullSysErr = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510SysErr->Clone("graphYieldPi0CombPbPb0510FullSysErr");
	Double_t* Yield0510 = graphYieldPi0CombPbPb0510FullSysErr->GetY();
	Double_t* SysErrorHigh0510 = graphYieldPi0CombPbPb0510FullSysErr->GetEYhigh();
	Double_t* SysErrorLow0510 = graphYieldPi0CombPbPb0510FullSysErr->GetEYlow();
	for (Int_t i = 0; i < graphYieldPi0CombPbPb0510FullSysErr->GetN(); i++){
		Double_t relativeErr = SysErrorHigh0510[i]/Yield0510[i]*100;
		Double_t relativeErrNew = TMath::Sqrt(TMath::Power(SysErrorHigh0510[i]/Yield0510[i]*100,2)+  commonCentralityErr0510*commonCentralityErr0510);
		Double_t absoluteErr = SysErrorHigh0510[i];
		Double_t absoluteErrNew = relativeErrNew/100*Yield0510[i];
		cout << i << "\t" <<relativeErr << "\t" << relativeErrNew << "\t" << absoluteErr << "\t" << absoluteErrNew <<endl;
		SysErrorHigh0510[i] = absoluteErrNew;
		SysErrorLow0510[i] = absoluteErrNew;
	}
	TGraphAsymmErrors* graphYieldPi0CombPbPb0510FullSysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0510FullSysErr->Clone("graphYieldPi0CombPbPb0510FullSysErrUnscaled");

	TGraphAsymmErrors* graphYieldPi0CombPbPb0010FullSysErr = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010SysErr->Clone("graphYieldPi0CombPbPb0010FullSysErr");
	Double_t* Yield0010 = graphYieldPi0CombPbPb0010FullSysErr->GetY();
	Double_t* SysErrorHigh0010 = graphYieldPi0CombPbPb0010FullSysErr->GetEYhigh();
	Double_t* SysErrorLow0010 = graphYieldPi0CombPbPb0010FullSysErr->GetEYlow();
	for (Int_t i = 0; i < graphYieldPi0CombPbPb0010FullSysErr->GetN(); i++){
		Double_t relativeErr = SysErrorHigh0010[i]/Yield0010[i]*100;
		Double_t relativeErrNew = TMath::Sqrt(TMath::Power(SysErrorHigh0010[i]/Yield0010[i]*100,2)+  commonCentralityErr0010*commonCentralityErr0010);
		Double_t absoluteErr = SysErrorHigh0010[i];
		Double_t absoluteErrNew = relativeErrNew/100*Yield0010[i];
		cout << i << "\t" <<relativeErr << "\t" << relativeErrNew << "\t" << absoluteErr << "\t" << absoluteErrNew <<endl;
		SysErrorHigh0010[i] = absoluteErrNew;
		SysErrorLow0010[i] = absoluteErrNew;
	}
	TGraphAsymmErrors* graphYieldPi0CombPbPb0010FullSysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb0010FullSysErr->Clone("graphYieldPi0CombPbPb0010FullSysErrUnscaled");

	TGraphAsymmErrors* graphYieldPi0CombPbPb1020FullSysErr = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020SysErr->Clone("graphYieldPi0CombPbPb1020FullSysErr");
	Double_t* Yield1020 = graphYieldPi0CombPbPb1020FullSysErr->GetY();
	Double_t* SysErrorHigh1020 = graphYieldPi0CombPbPb1020FullSysErr->GetEYhigh();
	Double_t* SysErrorLow1020 = graphYieldPi0CombPbPb1020FullSysErr->GetEYlow();
	for (Int_t i = 0; i < graphYieldPi0CombPbPb1020FullSysErr->GetN(); i++){
		Double_t relativeErr = SysErrorHigh1020[i]/Yield1020[i]*100;
		Double_t relativeErrNew = TMath::Sqrt(TMath::Power(SysErrorHigh1020[i]/Yield1020[i]*100,2)+  commonCentralityErr1020*commonCentralityErr1020);
		Double_t absoluteErr = SysErrorHigh1020[i];
		Double_t absoluteErrNew = relativeErrNew/100*Yield1020[i];
		cout << i << "\t" <<relativeErr << "\t" << relativeErrNew << "\t" << absoluteErr << "\t" << absoluteErrNew <<endl;
		SysErrorHigh1020[i] = absoluteErrNew;
		SysErrorLow1020[i] = absoluteErrNew;
	}
	TGraphAsymmErrors* graphYieldPi0CombPbPb1020FullSysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb1020FullSysErr->Clone("graphYieldPi0CombPbPb1020FullSysErrUnscaled");

	TGraphAsymmErrors* graphYieldPi0CombPbPb2040FullSysErr = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040SysErr->Clone("graphYieldPi0CombPbPb2040FullSysErr");
	Double_t* Yield2040 = graphYieldPi0CombPbPb2040FullSysErr->GetY();
	Double_t* SysErrorHigh2040 = graphYieldPi0CombPbPb2040FullSysErr->GetEYhigh();
	Double_t* SysErrorLow2040 = graphYieldPi0CombPbPb2040FullSysErr->GetEYlow();
	for (Int_t i = 0; i < graphYieldPi0CombPbPb2040FullSysErr->GetN(); i++){
		Double_t relativeErr = SysErrorHigh2040[i]/Yield2040[i]*100;
		Double_t relativeErrNew = TMath::Sqrt(TMath::Power(SysErrorHigh2040[i]/Yield2040[i]*100,2)+  commonCentralityErr2040*commonCentralityErr2040);
		Double_t absoluteErr = SysErrorHigh2040[i];
		Double_t absoluteErrNew = relativeErrNew/100*Yield2040[i];
		cout << i << "\t" <<relativeErr << "\t" << relativeErrNew << "\t" << absoluteErr << "\t" << absoluteErrNew <<endl;
		SysErrorHigh2040[i] = absoluteErrNew;
		SysErrorLow2040[i] = absoluteErrNew;
	}
	TGraphAsymmErrors* graphYieldPi0CombPbPb2040FullSysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb2040FullSysErr->Clone("graphYieldPi0CombPbPb2040FullSysErrUnscaled");

	TGraphAsymmErrors* graphYieldPi0CombPbPb4060FullSysErr = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060SysErr->Clone("graphYieldPi0CombPbPb4060FullSysErr");
	Double_t* Yield4060 = graphYieldPi0CombPbPb4060FullSysErr->GetY();
	Double_t* SysErrorHigh4060 = graphYieldPi0CombPbPb4060FullSysErr->GetEYhigh();
	Double_t* SysErrorLow4060 = graphYieldPi0CombPbPb4060FullSysErr->GetEYlow();
	for (Int_t i = 0; i < graphYieldPi0CombPbPb4060FullSysErr->GetN(); i++){
		Double_t relativeErr = SysErrorHigh4060[i]/Yield4060[i]*100;
		Double_t relativeErrNew = TMath::Sqrt(TMath::Power(SysErrorHigh4060[i]/Yield4060[i]*100,2)+  commonCentralityErr4060*commonCentralityErr4060);
		Double_t absoluteErr = SysErrorHigh4060[i];
		Double_t absoluteErrNew = relativeErrNew/100*Yield4060[i];
		cout << i << "\t" <<relativeErr << "\t" << relativeErrNew << "\t" << absoluteErr << "\t" << absoluteErrNew <<endl;
		SysErrorHigh4060[i] = absoluteErrNew;
		SysErrorLow4060[i] = absoluteErrNew;
	}
	TGraphAsymmErrors* graphYieldPi0CombPbPb4060FullSysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb4060FullSysErr->Clone("graphYieldPi0CombPbPb4060FullSysErrUnscaled");

	TGraphAsymmErrors* graphYieldPi0CombPbPb6080FullSysErr = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080SysErr->Clone("graphYieldPi0CombPbPb6080FullSysErr");
	Double_t* Yield6080 = graphYieldPi0CombPbPb6080FullSysErr->GetY();
	Double_t* SysErrorHigh6080 = graphYieldPi0CombPbPb6080FullSysErr->GetEYhigh();
	Double_t* SysErrorLow6080 = graphYieldPi0CombPbPb6080FullSysErr->GetEYlow();
	for (Int_t i = 0; i < graphYieldPi0CombPbPb6080FullSysErr->GetN(); i++){
		Double_t relativeErr = SysErrorHigh6080[i]/Yield6080[i]*100;
		Double_t relativeErrNew = TMath::Sqrt(TMath::Power(SysErrorHigh6080[i]/Yield6080[i]*100,2)+  commonCentralityErr6080*commonCentralityErr6080);
		Double_t absoluteErr = SysErrorHigh6080[i];
		Double_t absoluteErrNew = relativeErrNew/100*Yield6080[i];
		cout << i << "\t" <<relativeErr << "\t" << relativeErrNew << "\t" << absoluteErr << "\t" << absoluteErrNew <<endl;
		SysErrorHigh6080[i] = absoluteErrNew;
		SysErrorLow6080[i] = absoluteErrNew;
	}
	TGraphAsymmErrors* graphYieldPi0CombPbPb6080FullSysErrUnscaled = (TGraphAsymmErrors*)graphYieldPi0CombPbPb6080FullSysErr->Clone("graphYieldPi0CombPbPb6080FullSysErrUnscaled");

	Double_t normErr0005 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0005/tAA0005),2)+pow((commonCentralityErr0005/100),2),0.5);
	Double_t normErr0010 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0010/tAA0010),2)+pow((commonCentralityErr0010/100),2),0.5);
	Double_t normErr0510 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0510/tAA0510),2)+pow((commonCentralityErr0510/100),2),0.5);
	Double_t normErr1020 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr1020/tAA1020),2)+pow((commonCentralityErr1020/100),2),0.5);
	Double_t normErr2040 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2040/tAA2040),2)+pow((commonCentralityErr2040/100),2),0.5);
	Double_t normErr4060 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr4060/tAA4060),2)+pow((commonCentralityErr4060/100),2),0.5);
	Double_t normErr6080 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr6080/tAA6080),2)+pow((commonCentralityErr6080/100),2),0.5);

	cout << "norm Err \t" << normErr0005 << "\t" << normErr0010 << "\t" << normErr0510 << "\t" << normErr1020 << "\t" << normErr2040 << "\t" << normErr4060 << "\t" << normErr6080 << endl;
	//    return;

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


	TCanvas* canvasMassPlusFWHMDifPads = new TCanvas("canvasMassPlusFWHMDifPads","",200,10,2700,2500);  // gives the page size
	DrawGammaCanvasSettings( canvasMassPlusFWHMDifPads, 0.00, 0.00, 0.00, 0.0);   
	canvasMassPlusFWHMDifPads->cd();

	TPad* padFWHM = new TPad("padFWHM", "", 0., 0.56, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padFWHM, 0.11, 0.005, 0.02, 0.);
	padFWHM->SetLogx();
	padFWHM->Draw();

	TPad* padMass = new TPad("padMass", "", 0., 0., 1., 0.56,-1, -1, -2);
	DrawGammaPadSettings( padMass, 0.11, 0.005, 0., 0.15);
	padMass->Draw();

	TPad* padMassLegendSingle1 = new TPad("padMassLegendSingle1", "", 0.15, 0.41, 0.42, 0.53,-1, -1, -2);
	DrawGammaPadSettings( padMassLegendSingle1, 0., 0., 0., 0.);
	padMassLegendSingle1->Draw();

	//    TPad* padFWHMLegend1 = new TPad("padFWHMLegend1", "", 0.62, 0.54, 0.97, 0.635,-1, -1, -2);
	//    DrawGammaPadSettings( padFWHMLegend1, 0., 0., 0., 0.);
	//    padFWHMLegend1->Draw();

	padFWHM->cd();
	padFWHM->SetLogx();

	TH2D *histo2DAllPi0FWHM;
	histo2DAllPi0FWHM = new TH2D("histo2DAllPi0FWHM", "histo2DAllPi0FWHM", 20,0.25,35. ,1000.,-30,60);
	histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,25.5);
	histo2DAllPi0FWHM->SetYTitle("peak width (MeV/#it{c}^{2})");
	histo2DAllPi0FWHM->SetXTitle("#it{p}_{T} (GeV/#it{c})");
	histo2DAllPi0FWHM->SetLineWidth(2);
	histo2DAllPi0FWHM->SetTitleFont(42,"Y");
	histo2DAllPi0FWHM->SetLabelFont(42,"Y");
	histo2DAllPi0FWHM->GetXaxis()->SetNdivisions(515,kTRUE);
	histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DAllPi0FWHM->GetYaxis()->SetLabelSize(0.074);
	histo2DAllPi0FWHM->GetYaxis()->SetTitleSize(0.085);  
	histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.015);  
	histo2DAllPi0FWHM->GetYaxis()->SetDecimals();
	histo2DAllPi0FWHM->GetYaxis()->SetLabelOffset(0.01);
	histo2DAllPi0FWHM->GetYaxis()->SetTitleOffset(0.67);
	histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.04);  
	histo2DAllPi0FWHM->SetTitle("");
	histo2DAllPi0FWHM->DrawCopy();

	DrawGammaSetMarker(histoPHOSWidthData0010, markerStylePHOS, markerSizeMass*1.5, colorPHOSMass, colorPHOSMass);
	histoPHOSWidthData0010->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoPHOSWidthMC0010, markerStylePHOSMC, markerSizeMass*1.5, colorPHOSMCMass , colorPHOSMCMass);
	histoPHOSWidthMC0010->DrawCopy("same,e1,p,x0"); 

	DrawGammaSetMarkerTGraphErr( graphEMCALWidthData0010, markerStyleEMCAL, markerSizeMass*1.5*1.5, colorEMCALMass, colorEMCALMass);
	graphEMCALWidthData0010->Draw("same,e1,p"); 
	DrawGammaSetMarkerTGraphErr(graphEMCALWidthMC0010, markerStyleEMCALMC, markerSizeMass*1.5*1.5, colorEMCALMCMass , colorEMCALMCMass);
	graphEMCALWidthMC0010->Draw("same,e1,p"); 
	//    
	DrawGammaSetMarker(histoPCMWidthData0010, markerStyleConv, markerSizeMass*1.5, colorConv, colorConv);
	histoPCMWidthData0010->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoPCMWidthMC0010, markerStyleConvMC, markerSizeMass*1.5, colorConvMC, colorConvMC);
	histoPCMWidthMC0010->DrawCopy("same,e1,p,x0"); 

	TLatex *labelMassPi0PbPb0010_2 = new TLatex(0.7,0.89,collisionSystemPbPb.Data());
	SetStyleTLatex( labelMassPi0PbPb0010_2, 0.07,4);
	labelMassPi0PbPb0010_2->SetTextFont(42);
	labelMassPi0PbPb0010_2->Draw();
	TLatex *labelMassPi0PbPb0010_3 = new TLatex(0.77,0.8,"centrality 0-10%");
	SetStyleTLatex( labelMassPi0PbPb0010_3, 0.07,4);
	labelMassPi0PbPb0010_3->SetTextFont(42);
	labelMassPi0PbPb0010_3->Draw();
		
	TLatex *labelLegendAMass0010 = new TLatex(0.15,0.04,"a)");
	SetStyleTLatex( labelLegendAMass0010, 0.07,4);
	labelLegendAMass0010->SetTextFont(42);
	labelLegendAMass0010->Draw();


	//********************************** Defintion of the Legend **************************************************   
	Double_t columnsLegendFWHM2[4]    = {0.15,0.4,0.48,0.39};
	Double_t columnsLegendFWHM2Abs[4]    = {4,1.5,2.3,12};
	//    Double_t rowsLegendFWHM2[3]       = {0.66,0.33,0.0};
	Double_t rowsLegendFWHM2[4]       = {0.9,0.82,0.76,0.7};
	Double_t rowsLegendFWHM2Abs[4]       = {0.2,22,20.,18.};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnFWHM2 = 0.07;
	Size_t textSizeTopRowFWHM2  = 0.07; 
	//****************** Scale factors ******************
	Double_t scaleMarkerFWHM2      = 1.2;

	//    padFWHMLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textFWHM2CTS = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[1],"PCM (FWHM/2.35)");
	SetStyleTLatex( textFWHM2CTS, textSizeLeftColumnFWHM2,4);
	textFWHM2CTS->SetTextFont(42);
	textFWHM2CTS->Draw();
	TLatex *textFWHM2PHOS = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[2],"PHOS (#sigma)");
	SetStyleTLatex( textFWHM2PHOS, textSizeLeftColumnFWHM2,4);
	textFWHM2PHOS->SetTextFont(42);
	textFWHM2PHOS->Draw();
	TLatex *textFWHM2EMCAL = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[3],"EMCal (#sigma)");
	SetStyleTLatex( textFWHM2EMCAL, textSizeLeftColumnFWHM2,4);
	textFWHM2EMCAL->SetTextFont(42);
	textFWHM2EMCAL->Draw();

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

	TMarker* markerCTSPi0FWHM2 = CreateMarkerFromHisto(histoPCMWidthData0010,columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[1] ,scaleMarkerFWHM2);
	markerCTSPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[1]);
	TMarker* markerPHOSPi0FWHM2 = CreateMarkerFromHisto(histoPHOSWidthData0010,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[2],scaleMarkerFWHM2);
	markerPHOSPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[2]);
	TMarker* markerEMCALPi0FWHM2 = CreateMarkerFromGraph(graphEMCALWidthData0010,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[3],scaleMarkerFWHM2);
	markerEMCALPi0FWHM2->DrawMarker(columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[3]);

	TMarker* markerCTSPi0FWHM2MC = CreateMarkerFromHisto(histoPCMWidthMC0010,columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[1],scaleMarkerFWHM2);
	markerCTSPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[1]);
	TMarker* markerPHOSPi0FWHM2MC = CreateMarkerFromHisto(histoPHOSWidthMC0010,columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[2] ,scaleMarkerFWHM2);
	markerPHOSPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[2]);
	TMarker* markerEMCALPi0FWHM2MC = CreateMarkerFromGraph(graphEMCALWidthMC0010,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[3] ,scaleMarkerFWHM2);
	markerEMCALPi0FWHM2MC->DrawMarker(columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[3]);

	padMass->cd();
	padMass->SetLogx();

	TH2D *histo2DAllPi0Mass;
	histo2DAllPi0Mass = new TH2D("histo2DAllPi0Mass", "histo2DAllPi0Mass", 20,0.25,35. ,1000.,125.,170);
	histo2DAllPi0Mass->GetYaxis()->SetRangeUser(125.,156);
	histo2DAllPi0Mass->SetYTitle("peak position (MeV/#it{c}^{2})");
	histo2DAllPi0Mass->SetXTitle("#it{p}_{T} (GeV/#it{c})");
	histo2DAllPi0Mass->SetTitleFont(42,"X");
	histo2DAllPi0Mass->SetTitleFont(42,"Y");
	histo2DAllPi0Mass->SetLabelFont(42,"X");
	histo2DAllPi0Mass->SetLabelFont(42,"Y");
	histo2DAllPi0Mass->SetLineWidth(2);
	histo2DAllPi0Mass->GetXaxis()->SetNdivisions(515,kTRUE);
	histo2DAllPi0Mass->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DAllPi0Mass->GetYaxis()->SetLabelSize(0.063);
	histo2DAllPi0Mass->GetYaxis()->SetTitleSize(0.068);   
	histo2DAllPi0Mass->GetYaxis()->SetDecimals();
	histo2DAllPi0Mass->GetYaxis()->SetLabelOffset(0.01);
	histo2DAllPi0Mass->GetYaxis()->SetTitleOffset(0.82);
	histo2DAllPi0Mass->GetXaxis()->SetTitleSize(0.068);   
	histo2DAllPi0Mass->GetXaxis()->SetLabelSize(0.063);
	histo2DAllPi0Mass->GetYaxis()->SetTickLength(0.02);  
	histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.04);  
	histo2DAllPi0Mass->SetTitle("");
	histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
	histo2DAllPi0Mass->GetXaxis()->SetTitleOffset(0.85);
	histo2DAllPi0Mass->DrawCopy();

	DrawGammaSetMarker(histoPHOSMassData0010 , markerStylePHOS, markerSizeMass*1.5, colorPHOSMass,colorPHOSMass);
	histoPHOSMassData0010->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoPHOSMassMC0010 , markerStylePHOSMC, markerSizeMass*1.5, colorPHOSMCMass, colorPHOSMCMass);
	histoPHOSMassMC0010->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarkerTGraphErr( graphEMCALMassData0010, markerStyleEMCAL, markerSizeMass*1.5*1.5, colorEMCALMass, colorEMCALMass);
	graphEMCALMassData0010->Draw("same,e1,p"); 
	DrawGammaSetMarkerTGraphErr(graphEMCALMassMC0010, markerStyleEMCALMC, markerSizeMass*1.5*1.5, colorEMCALMCMass , colorEMCALMCMass);
	graphEMCALMassMC0010->Draw("same,e1,p"); 

	DrawGammaSetMarker(histoPCMMassData0010 , markerStyleConv, markerSizeMass*1.5, colorConv, colorConv);               
	histoPCMMassData0010->DrawCopy("same,e1,p,x0"); 
	DrawGammaSetMarker(histoPCMMassMC0010 , markerStyleConvMC , markerSizeMass*1.5, colorConvMC, colorConvMC);                
	histoPCMMassMC0010->DrawCopy("same,e1,p,x0"); 
	//    DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,colorConv);



	DrawGammaLines(0.25, 35. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1,colorConv);



	TLatex *labelLegendBMass0010 = new TLatex(0.15,0.18,"b)");
	SetStyleTLatex( labelLegendBMass0010, 0.06,4);
	labelLegendBMass0010->SetTextFont(42);
	labelLegendBMass0010->Draw();

	//********************************** Defintion of the Legend **************************************************   
	Double_t columnsLegendMass2[3]   = {0.,0.45,0.75};
	Double_t rowsLegendMass2[4]      = {0.75,0.5,0.25,0.0};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnMass2   = 0.26;
	Size_t textSizeSecondRowMass2    = 0.25;
	//******************* Offsets ***********************
	Double_t offsetMarkerXMass2   = 0.1;
	Double_t offsetMarkerYMass2   = 0.1;
	//****************** Scale factors ******************
	Double_t scaleMarkerMass2     = 1.2;

	padMassLegendSingle1->cd();
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

	TMarker* markerCTSPi0Mass2 = CreateMarkerFromHisto(histoPCMMassData0010,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerCTSPi0Mass2->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
	TMarker* markerPHOSPi0Mass2 = CreateMarkerFromHisto(histoPHOSMassData0010,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerPHOSPi0Mass2->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
	TMarker* markerEMCALPi0Mass2 = CreateMarkerFromGraph(graphEMCALMassData0010,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerEMCALPi0Mass2->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
	//    
	TMarker* markerCTSPi0MassMC2 = CreateMarkerFromHisto(histoPCMMassMC0010,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerCTSPi0MassMC2->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
	TMarker* markerPHOSPi0MassMC2 = CreateMarkerFromHisto(histoPHOSMassMC0010,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerPHOSPi0MassMC2->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
	TMarker* markerEMCALPi0MassMC2 = CreateMarkerFromGraph(graphEMCALMassMC0010,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerEMCALPi0MassMC2->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
	//    
	canvasMassPlusFWHMDifPads->Update();
	canvasMassPlusFWHMDifPads->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurementsDP_0010_Paper.%s",outputDir.Data(), prefix2.Data(),suffix.Data()));
	canvasMassPlusFWHMDifPads->SaveAs(Form("%s/%s_Pi0_CombinedMassAndWidthALLMeasurementsDP_0010_Paper.C",outputDir.Data(), prefix2.Data()));
	canvasMassPlusFWHMDifPads->cd();

		
	TCanvas * canvas6PartMassWidth = new TCanvas("canvas6PartMassWidth","",10,10,2400,1300);  // gives the page size		
	canvas6PartMassWidth->cd();
	DrawGammaCanvasSettings( canvas6PartMassWidth, 0.13, 0.0, 0.02, 0.09);

	TPad* pad6PartMassWidth1 = new TPad("pad6PartMassWidth1", "", 0., 0.54, 0.37, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth1, 0.16, 0.0, 0.02, 0.);
	pad6PartMassWidth1->Draw();
	TPad* pad6PartMassWidth2 = new TPad("pad6PartMassWidth2", "", 0., 0., 0.37, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth2, 0.16, 0.0, 0., 0.15);
	pad6PartMassWidth2->Draw();

	TPad* pad6PartMassWidth3 = new TPad("pad6PartMassWidth3", "", 0.37, 0.54, 0.685, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth3, 0.0, 0.0, 0.02, 0.);
	pad6PartMassWidth3->Draw();
	TPad* pad6PartMassWidth4 = new TPad("pad6PartMassWidth4", "", 0.37, 0., 0.685, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth4, 0.0, 0.0, 0., 0.15);
	pad6PartMassWidth4->Draw();

	TPad* pad6PartMassWidth5 = new TPad("pad6PartMassWidth5", "", 0.685, 0.54, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth5, 0.0, 0.02, 0.02, 0.);
	pad6PartMassWidth5->Draw();
	TPad* pad6PartMassWidth6 = new TPad("pad6PartMassWidth6", "", 0.685, 0., 1., 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth6, 0.0, 0.02, 0., 0.15);
	pad6PartMassWidth6->Draw();

	Int_t textSizeLabelsPixelMass = 50;
	Double_t marginMass = 0.13*2400;
	Double_t textsizeLabelsMass = 0;
	Double_t textsizeFacMass = 0;
	Double_t textsizeLabelsWidth = 0;
	Double_t textsizeFacWidth = 0;

	if (pad6PartMassWidth1->XtoPixel(pad6PartMassWidth1->GetX2()) < pad6PartMassWidth1->YtoPixel(pad6PartMassWidth1->GetY1())){
		textsizeLabelsWidth = (Double_t)textSizeLabelsPixelMass/pad6PartMassWidth1->XtoPixel(pad6PartMassWidth1->GetX2()) ;
		textsizeFacWidth = (Double_t)1./pad6PartMassWidth1->XtoPixel(pad6PartMassWidth1->GetX2()) ;
	} else {
		textsizeLabelsWidth = (Double_t)textSizeLabelsPixelMass/pad6PartMassWidth1->YtoPixel(pad6PartMassWidth1->GetY1());
		textsizeFacWidth = (Double_t)1./pad6PartMassWidth1->YtoPixel(pad6PartMassWidth1->GetY1());
	}
	if (pad6PartMassWidth2->XtoPixel(pad6PartMassWidth2->GetX2()) < pad6PartMassWidth2->YtoPixel(pad6PartMassWidth2->GetY1())){
		textsizeLabelsMass = (Double_t)textSizeLabelsPixelMass/pad6PartMassWidth2->XtoPixel(pad6PartMassWidth2->GetX2()) ;
		textsizeFacMass = (Double_t)1./pad6PartMassWidth2->XtoPixel(pad6PartMassWidth2->GetX2()) ;
	} else {
		textsizeLabelsMass = (Double_t)textSizeLabelsPixelMass/pad6PartMassWidth2->YtoPixel(pad6PartMassWidth2->GetY1());
		textsizeFacMass = (Double_t)1./pad6PartMassWidth2->YtoPixel(pad6PartMassWidth2->GetY1());
	}

	cout << textsizeLabelsMass << endl;

	TPad* padFWHMLegend1 = new TPad("padFWHMLegend1", "", 0.07, 0.815, 0.35, 0.93,-1, -1, -2);
	DrawGammaPadSettings( padFWHMLegend1, 0., 0., 0., 0.);
	padFWHMLegend1->Draw();


	TH2D *histo2DPi0FWHM;
	histo2DPi0FWHM = new TH2D("histo2DPi0FWHM", "histo2DPi0FWHM", 20,0.35,20. ,1000.,-30,40);
	SetStyleHistoTH2ForGraphs(histo2DPi0FWHM, "#it{p}_{T} (GeV/#it{c})","peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth,textsizeLabelsWidth, 0.85*textsizeLabelsWidth,textsizeLabelsWidth, 1,0.5/(textsizeFacWidth*marginMass), 515, 504); 
	histo2DPi0FWHM->GetYaxis()->SetRangeUser(-1.,18);
	histo2DPi0FWHM->GetYaxis()->SetLabelOffset(0.01);
	histo2DPi0FWHM->GetYaxis()->SetLabelFont(42);
	histo2DPi0FWHM->GetXaxis()->SetLabelFont(42);

	TH2D *histo2DPi0Mass;
	histo2DPi0Mass = new TH2D("histo2DPi0Mass", "histo2DPi0Mass", 20,0.35,20. ,1000.,125.,150);
	SetStyleHistoTH2ForGraphs(histo2DPi0Mass, "#it{p}_{T} (GeV/#it{c})","peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass,textsizeLabelsMass, 0.85*textsizeLabelsMass,textsizeLabelsMass, 0.9,0.5/(textsizeFacMass*marginMass), 515, 510); 
	histo2DPi0Mass->GetYaxis()->SetRangeUser(128.,143.5);
	histo2DPi0Mass->GetXaxis()->SetLabelOffset(-0.02);
	histo2DPi0Mass->GetYaxis()->SetLabelOffset(0.01);
	histo2DPi0Mass->GetYaxis()->SetLabelFont(42);
	histo2DPi0Mass->GetXaxis()->SetLabelFont(42);

	pad6PartMassWidth1->cd();
	pad6PartMassWidth1->SetLogx();
	histo2DPi0FWHM->DrawCopy();
		
	DrawGammaSetMarker(histoPHOSWidthDataPP, markerStylePHOS, markerSizeMass, colorPHOSMass, colorPHOSMass);
	histoPHOSWidthDataPP->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPHOSWidthMCPP, markerStylePHOSMC, markerSizeMass, colorPHOSMCMass , colorPHOSMCMass);
	histoPHOSWidthMCPP->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoPCMWidthDataPP, markerStyleConv, markerSizeMass, colorConv, colorConv);
	histoPCMWidthDataPP->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMWidthMCPP, markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);
	histoPCMWidthMCPP->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoPHOSWidthData0010, markerStylePHOS, markerSizeMass, colorPHOSMass, colorPHOSMass);
	DrawGammaSetMarker(histoPHOSWidthMC0010, markerStylePHOSMC, markerSizeMass, colorPHOSMCMass , colorPHOSMCMass);

	TLatex *labelMassPi0PP = new TLatex(0.2,0.88,collisionSystemPP.Data());
	SetStyleTLatex( labelMassPi0PP, 0.85*textsizeLabelsWidth,4);
	labelMassPi0PP->Draw();
	TLatex *labelLegendAMass = new TLatex(0.92,0.88,"a)");
	SetStyleTLatex( labelLegendAMass,0.85*textsizeLabelsWidth,4);
	labelLegendAMass->Draw();

	//********************************** Defintion of the Legend **************************************************	
	Double_t columnsLegendFWHM[4] 	= {0.,0.2,0.37,0.55};
	Double_t rowsLegendFWHM[3] 		= {0.66,0.33,0.0};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnFWHM	= 0.301;
	Size_t textSizeTopRowFWHM	= 0.301; 
	Size_t textSizeSecondRowFWHM 	= 0.301;
	//******************* Offsets ***********************
	Double_t offsetMarkerXFWHM	= 0.07;
	Double_t offsetMarkerYFWHM	= 0.07;
	//****************** Scale factors ******************
	Double_t scaleMarkerFWHM		= 1.;

	padFWHMLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textFWHMCTS = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM[1],"PCM");
	SetStyleTLatex( textFWHMCTS, textSizeLeftColumnFWHM,4);
	textFWHMCTS->Draw();
	TLatex *textFWHMPHOS = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM[2],"PHOS");
	SetStyleTLatex( textFWHMPHOS, textSizeLeftColumnFWHM,4);
	textFWHMPHOS->Draw();

	//****************** second Column *************************************************
	TLatex *textFWHMData2 = new TLatex(columnsLegendFWHM[1],rowsLegendFWHM[0] ,"Data");
	SetStyleTLatex( textFWHMData2, textSizeTopRowFWHM ,4);
	textFWHMData2->Draw();
	TLatex *textFWHMMC2 = new TLatex(columnsLegendFWHM[2] ,rowsLegendFWHM[0],"MC");
	SetStyleTLatex( textFWHMMC2, textSizeTopRowFWHM,4);
	textFWHMMC2->Draw();

	TMarker* markerCTSPi0FWHM = CreateMarkerFromHisto(histoPCMWidthDataPP,columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
	markerCTSPi0FWHM->DrawMarker(columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM);
	TMarker* markerPHOSPi0FWHM = CreateMarkerFromHisto(histoPHOSWidthData0010,columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
	markerPHOSPi0FWHM->DrawMarker(columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM);

	TMarker* markerCTSPi0FWHMMC = CreateMarkerFromHisto(histoPCMWidthMCPP,columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
	markerCTSPi0FWHMMC->DrawMarker(columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM);
	TMarker* markerPHOSPi0FWHMMC = CreateMarkerFromHisto(histoPHOSWidthMC0010,columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
	markerPHOSPi0FWHMMC->DrawMarker(columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM);

	TLatex *textWidthConv2 = new TLatex(columnsLegendFWHM[3],rowsLegendFWHM[1] ,"FWHM/2.35");
	SetStyleTLatex( textWidthConv2, textSizeSecondRowFWHM,4);
	textWidthConv2->Draw();
	TLatex *textWidthPHOS2 = new TLatex(columnsLegendFWHM[3] ,rowsLegendFWHM[2],"#sigma");
	SetStyleTLatex( textWidthPHOS2, textSizeSecondRowFWHM,4);
	textWidthPHOS2->Draw();


	pad6PartMassWidth1->Update();
	pad6PartMassWidth2->cd();
	pad6PartMassWidth2->SetLogx();
	histo2DPi0Mass->DrawCopy();

	DrawGammaSetMarker(histoPHOSMassDataPP , markerStylePHOS, markerSizeMass, colorPHOSMass,colorPHOSMass);
	histoPHOSMassDataPP->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPHOSMassMCPP , markerStylePHOSMC, markerSizeMass, colorPHOSMCMass, colorPHOSMCMass);
	histoPHOSMassMCPP->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoPHOSMassData0005 , markerStylePHOS, markerSizeMass, colorPHOSMass,colorPHOSMass);
	DrawGammaSetMarker(histoPHOSMassMC0005 , markerStylePHOSMC, markerSizeMass, colorPHOSMCMass, colorPHOSMCMass);

	DrawGammaSetMarker(histoPCMMassDataPP , markerStyleConv, markerSizeMass, colorConv, colorConv);					 
	histoPCMMassDataPP->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMMassMCPP , markerStyleConvMC , markerSizeMass, colorConvMC, colorConvMC);					 
	histoPCMMassMCPP->DrawCopy("same,e1,p"); 

	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,colorConv);
	TLatex *labelLegendDMass = new TLatex(0.92,0.9,"d)");
	SetStyleTLatex( labelLegendDMass, 0.85*textsizeLabelsMass,4);
	labelLegendDMass->Draw();

	pad6PartMassWidth2->Update();

	pad6PartMassWidth5->cd();
	pad6PartMassWidth5->SetLogx();
	histo2DPi0FWHM->DrawCopy();
		
	DrawGammaSetMarker(histoPHOSWidthData0005, markerStylePHOS, markerSizeMass, colorPHOSMass, colorPHOSMass);
	histoPHOSWidthData0005->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPHOSWidthMC0005, markerStylePHOSMC, markerSizeMass, colorPHOSMCMass , colorPHOSMCMass);
	histoPHOSWidthMC0005->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoPCMWidthData0005, markerStyleConv, markerSizeMass, colorConv, colorConv);
	histoPCMWidthData0005->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMWidthMC0005, markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);
	histoPCMWidthMC0005->DrawCopy("same,e1,p"); 

	TLatex *labelMassPi0PbPb0005 = new TLatex(0.05,0.88,collisionSystemCent0005.Data());
	SetStyleTLatex( labelMassPi0PbPb0005, 0.85*textsizeLabelsWidth,4);
	labelMassPi0PbPb0005->Draw();
	TLatex *labelLegendBMass = new TLatex(0.89,0.88,"c)");
	SetStyleTLatex( labelLegendBMass, 0.85*textsizeLabelsWidth,4);
	labelLegendBMass->Draw();

	pad6PartMassWidth5->Update();
	pad6PartMassWidth6->cd();
	pad6PartMassWidth6->SetLogx();
	histo2DPi0Mass->DrawCopy();

	DrawGammaSetMarker(histoPHOSMassData0005 , markerStylePHOS, markerSizeMass, colorPHOSMass,colorPHOSMass);
	histoPHOSMassData0005->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPHOSMassMC0005 , markerStylePHOSMC, markerSizeMass, colorPHOSMCMass, colorPHOSMCMass);
	histoPHOSMassMC0005->DrawCopy("same,e1,p"); 
		
	DrawGammaSetMarker(histoPCMMassData0005 , markerStyleConv, markerSizeMass, colorConv, colorConv);					 
	histoPCMMassData0005->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMMassMC0005 , markerStyleConvMC , markerSizeMass, colorConvMC, colorConvMC);					 
	histoPCMMassMC0005->DrawCopy("same,e1,p"); 
	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,colorConv);
	TLatex *labelLegendFMass = new TLatex(0.89,0.9,"f)");
	SetStyleTLatex( labelLegendFMass, 0.85*textsizeLabelsMass,4);
	labelLegendFMass->Draw();

	pad6PartMassWidth6->Update();

	pad6PartMassWidth3->cd();
	pad6PartMassWidth3->SetLogx();
	histo2DPi0FWHM->DrawCopy();
		
	DrawGammaSetMarker(histoPHOSWidthData6080, markerStylePHOS, markerSizeMass, colorPHOSMass, colorPHOSMass);
	histoPHOSWidthData6080->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPHOSWidthMC6080, markerStylePHOSMC, markerSizeMass, colorPHOSMCMass , colorPHOSMCMass);
	histoPHOSWidthMC6080->DrawCopy("same,e1,p"); 

	DrawGammaSetMarker(histoPCMWidthData6080, markerStyleConv, markerSizeMass, colorConv, colorConv);
	histoPCMWidthData6080->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMWidthMC6080, markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);
	histoPCMWidthMC6080->DrawCopy("same,e1,p"); 

	TLatex *labelMassPi0PbPb6080 = new TLatex(0.05,0.88,collisionSystemCent6080.Data());
	SetStyleTLatex( labelMassPi0PbPb6080, 0.85*textsizeLabelsWidth,4);
	labelMassPi0PbPb6080->Draw();
	TLatex *labelLegendCMass = new TLatex(0.91,0.88,"b)");
	SetStyleTLatex( labelLegendCMass, 0.85*textsizeLabelsWidth,4);
	labelLegendCMass->Draw();

	// 	TLatex *labelMassPi0LabelPbPb6080 = new TLatex(0.05,0.83,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelMassPi0LabelPbPb6080, 0.052,4);
	// 	labelMassPi0LabelPbPb6080->Draw();

	pad6PartMassWidth3->Update();
	pad6PartMassWidth4->cd();
	pad6PartMassWidth4->SetLogx();
	histo2DPi0Mass->DrawCopy();

	DrawGammaSetMarker(histoPHOSMassData6080 , markerStylePHOS, markerSizeMass, colorPHOSMass,colorPHOSMass);
	histoPHOSMassData6080->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPHOSMassMC6080 , markerStylePHOSMC, markerSizeMass, colorPHOSMCMass, colorPHOSMCMass);
	histoPHOSMassMC6080->DrawCopy("same,e1,p"); 
		
	DrawGammaSetMarker(histoPCMMassData6080 , markerStyleConv, markerSizeMass, colorConv, colorConv);					 
	histoPCMMassData6080->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMMassMC6080 , markerStyleConvMC , markerSizeMass, colorConvMC, colorConvMC);					 
	histoPCMMassMC6080->DrawCopy("same,e1,p"); 
	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,colorConv);
	TLatex *labelLegendEMass = new TLatex(0.91,0.9,"e)");
	SetStyleTLatex( labelLegendEMass, 0.85*textsizeLabelsMass,4);
	labelLegendEMass->Draw();

	pad6PartMassWidth4->Update();

	canvas6PartMassWidth->Update();	
	canvas6PartMassWidth->SaveAs(Form("%s/MassWidth_6Parted_Paper.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartMassWidth1;	
	delete pad6PartMassWidth2;	
	delete pad6PartMassWidth3;	
	delete pad6PartMassWidth4;	
	delete canvas6PartMassWidth;	

	cout << "here ............................." << endl;
		
	//****************************************************************************************************************************
	//************************************ Inv Yield only ratio Pi0 **************************************************************
	//****************************************************************************************************************************

	Double_t arrayBoundsXIndMeasRatio[4];
	Double_t arrayBoundsYIndMeasRatio[3];
	Double_t relativeMarginsIndMeasRatioX[3];
	Double_t relativeMarginsIndMeasRatioY[3];
	ReturnCorrectValuesForCanvasScaling(2400,1300, 3, 2,0.05, 0.005, 0.005,0.09,arrayBoundsXIndMeasRatio,arrayBoundsYIndMeasRatio,relativeMarginsIndMeasRatioX,relativeMarginsIndMeasRatioY);


	TCanvas * canvas6PartRatioIndMeas = new TCanvas("canvas6PartRatioIndMeas","",10,10,2400,1300);  // gives the page size		
	canvas6PartRatioIndMeas->cd();
	// 	DrawGammaCanvasSettings( canvas6PartRatioIndMeas, 0.13, 0.0, 0.02, 0.09);

	TPad* pad6PartRatioIndMeas1 = new TPad("pad6PartRatioIndMeas1", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas1, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRatioIndMeas1->Draw();
	TPad* pad6PartRatioIndMeas2 = new TPad("pad6PartRatioIndMeas2", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas2, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRatioIndMeas2->Draw();

	TPad* pad6PartRatioIndMeas3 = new TPad("pad6PartRatioIndMeas3", "", arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas3, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRatioIndMeas3->Draw();
	TPad* pad6PartRatioIndMeas4 = new TPad("pad6PartRatioIndMeas4", "", arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas4, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRatioIndMeas4->Draw();

	TPad* pad6PartRatioIndMeas5 = new TPad("pad6PartRatioIndMeas5", "", arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas5, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRatioIndMeas5->Draw();
	TPad* pad6PartRatioIndMeas6 = new TPad("pad6PartRatioIndMeas6", "", arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas6, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRatioIndMeas6->Draw();

	Int_t textSizeLabelsPixelRatio = 50;
	Double_t marginRatio = 0.16*2400;
	Double_t textsizeLabelsRatioUp = 0;
	Double_t textsizeFacRatioUp = 0;
	Double_t textsizeLabelsRatioDown = 0;
	Double_t textsizeFacRatioDown = 0;

	if (pad6PartRatioIndMeas1->XtoPixel(pad6PartRatioIndMeas1->GetX2()) < pad6PartRatioIndMeas1->YtoPixel(pad6PartRatioIndMeas1->GetY1())){
		textsizeLabelsRatioUp = (Double_t)textSizeLabelsPixelRatio/pad6PartRatioIndMeas1->XtoPixel(pad6PartRatioIndMeas1->GetX2()) ;
		textsizeFacRatioUp = (Double_t)1./pad6PartRatioIndMeas1->XtoPixel(pad6PartRatioIndMeas1->GetX2()) ;
	} else {
		textsizeLabelsRatioUp = (Double_t)textSizeLabelsPixelRatio/pad6PartRatioIndMeas1->YtoPixel(pad6PartRatioIndMeas1->GetY1());
		textsizeFacRatioUp = (Double_t)1./pad6PartRatioIndMeas1->YtoPixel(pad6PartRatioIndMeas1->GetY1());
	}
	if (pad6PartRatioIndMeas2->XtoPixel(pad6PartRatioIndMeas2->GetX2()) < pad6PartRatioIndMeas2->YtoPixel(pad6PartRatioIndMeas2->GetY1())){
		textsizeLabelsRatioDown = (Double_t)textSizeLabelsPixelRatio/pad6PartRatioIndMeas2->XtoPixel(pad6PartRatioIndMeas2->GetX2()) ;
		textsizeFacRatioDown = (Double_t)1./pad6PartRatioIndMeas2->XtoPixel(pad6PartRatioIndMeas2->GetX2()) ;
	} else {
		textsizeLabelsRatioDown = (Double_t)textSizeLabelsPixelRatio/pad6PartRatioIndMeas2->YtoPixel(pad6PartRatioIndMeas2->GetY1());
		textsizeFacRatioDown = (Double_t)1./pad6PartRatioIndMeas2->YtoPixel(pad6PartRatioIndMeas2->GetY1());
	}

	pad6PartRatioIndMeas1->cd();
	pad6PartRatioIndMeas1->SetLogx();
	TH2F * ratio2DInvXSectionOnlyPi0Up;
	ratio2DInvXSectionOnlyPi0Up = new TH2F("ratio2DInvXSectionOnlyPi0Up","ratio2DInvXSectionOnlyPi0Up",1000,0.33,20.,1000,0.05,2.4	);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionOnlyPi0Up, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsRatioUp,textsizeLabelsRatioUp,  0.85*textsizeLabelsRatioUp,textsizeLabelsRatioUp, 1,0.5/(textsizeFacRatioUp*marginRatio), 512, 505);
	ratio2DInvXSectionOnlyPi0Up->GetXaxis()->SetLabelOffset(-0.015);
	ratio2DInvXSectionOnlyPi0Up->GetYaxis()->SetLabelOffset(0.01);
	ratio2DInvXSectionOnlyPi0Up->GetYaxis()->SetLabelFont(42);
	ratio2DInvXSectionOnlyPi0Up->GetXaxis()->SetLabelFont(42);
	ratio2DInvXSectionOnlyPi0Up->DrawCopy(); 

	TH2F * ratio2DInvXSectionOnlyPi0Down;
	ratio2DInvXSectionOnlyPi0Down = new TH2F("ratio2DInvXSectionOnlyPi0Down","ratio2DInvXSectionOnlyPi0Down",1000,0.33,20.,1000,0.05,2.4	);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionOnlyPi0Down, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsRatioDown,textsizeLabelsRatioDown,  0.85*textsizeLabelsRatioDown,textsizeLabelsRatioDown, 1,0.5/(textsizeFacRatioDown*marginRatio), 512, 505);
	ratio2DInvXSectionOnlyPi0Down->GetXaxis()->SetLabelOffset(-0.015);
	ratio2DInvXSectionOnlyPi0Down->GetXaxis()->SetLabelFont(42);
	ratio2DInvXSectionOnlyPi0Down->GetYaxis()->SetLabelFont(42);
	ratio2DInvXSectionOnlyPi0Down->GetYaxis()->SetLabelOffset(0.01);


	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb0005Sys, markerStylePHOS,markerSizeSpectrum*0.6, colorPHOS , colorPHOS, widthLinesBoxes, kTRUE);
	graphRatioCombPHOSPi0PbPb0005Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb0005Sys, markerStyleConv,markerSizeSpectrum*0.6, colorConv , colorConv, widthLinesBoxes, kTRUE);
	graphRatioCombConvPi0PbPb0005Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb0005, markerStylePHOS,markerSizeSpectrum*2.*0.6, colorPHOS , colorPHOS);
	graphRatioCombPHOSPi0PbPb0005->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb0005 , markerStyleConv,markerSizeSpectrum*2.*0.6, colorConv , colorConv);
	graphRatioCombConvPi0PbPb0005->Draw("p,same,e1");
	DrawGammaLines(0., 20.,1., 1.,1.,kBlack,2);

	TLatex *labelRatioIndMeasPbPb0005 = new TLatex(0.18,0.88,collisionSystemCent0005.Data());
	SetStyleTLatex( labelRatioIndMeasPbPb0005, 0.85*textsizeLabelsRatioUp,4);
	labelRatioIndMeasPbPb0005->Draw();


	pad6PartRatioIndMeas1->Update();
	pad6PartRatioIndMeas2->cd();
	pad6PartRatioIndMeas2->SetLogx();
	ratio2DInvXSectionOnlyPi0Down->DrawCopy(); 
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb2040Sys, markerStylePHOS,markerSizeSpectrum*0.6, colorPHOS , colorPHOS, widthLinesBoxes, kTRUE);
	graphRatioCombPHOSPi0PbPb2040Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb2040Sys, markerStyleConv,markerSizeSpectrum*0.6, colorConv , colorConv, widthLinesBoxes, kTRUE);
	graphRatioCombConvPi0PbPb2040Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb2040, markerStylePHOS,markerSizeSpectrum*2.*0.6, colorPHOS , colorPHOS);
	graphRatioCombPHOSPi0PbPb2040->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb2040 , markerStyleConv,markerSizeSpectrum*2.*0.6, colorConv , colorConv);
	graphRatioCombConvPi0PbPb2040->Draw("p,same,e1");
	DrawGammaLines(0., 20.,1., 1.,1.,kBlack,2);

	TLegend* legendRatioInvSec = new TLegend(0.15,0.65,0.6,0.78);
	legendRatioInvSec->SetFillColor(0);
	legendRatioInvSec->SetLineColor(0);
	legendRatioInvSec->SetTextSize(0.85*textsizeLabelsRatioDown);
	legendRatioInvSec->SetNColumns(3);
	legendRatioInvSec->SetMargin(0.4);
	legendRatioInvSec->SetTextFont(42);
	TLatex *labelStat = new TLatex(0.37,0.8,"stat.");
	SetStyleTLatex( labelStat, 0.85*textsizeLabelsRatioDown,4);
	labelStat->Draw();
	TLatex *labelSyst = new TLatex(0.5,0.8,"syst.");
	SetStyleTLatex( labelSyst, 0.85*textsizeLabelsRatioDown,4);
	labelSyst->Draw();
	legendRatioInvSec->AddEntry((TObject*)0, "PCM","");
	legendRatioInvSec->AddEntry(graphRatioCombConvPi0PbPb0005,"    ","p");
	legendRatioInvSec->AddEntry(graphRatioCombConvPi0PbPb0005Sys," ","f");
	legendRatioInvSec->AddEntry((TObject*)0, "PHOS","");
	legendRatioInvSec->AddEntry(graphRatioCombPHOSPi0PbPb0005," ","p");
	legendRatioInvSec->AddEntry(graphRatioCombPHOSPi0PbPb0005Sys," ","f");
	legendRatioInvSec->Draw();

	TLatex *labelRatioIndMeasPbPb2040 = new TLatex(0.18,0.9,collisionSystemCent2040.Data());
	SetStyleTLatex( labelRatioIndMeasPbPb2040, 0.85*textsizeLabelsRatioDown,4);
	labelRatioIndMeasPbPb2040->Draw();

	pad6PartRatioIndMeas2->Update();

	pad6PartRatioIndMeas3->cd();
	pad6PartRatioIndMeas3->SetLogx();
	ratio2DInvXSectionOnlyPi0Up->DrawCopy(); 
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb0510Sys, markerStylePHOS,markerSizeSpectrum*0.6, colorPHOS , colorPHOS, widthLinesBoxes, kTRUE);
	graphRatioCombPHOSPi0PbPb0510Sys->Draw("E2same");

	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb0510Sys, markerStyleConv,markerSizeSpectrum*0.6, colorConv , colorConv, widthLinesBoxes, kTRUE);
	graphRatioCombConvPi0PbPb0510Sys->Draw("E2same");

	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb0510, markerStylePHOS,markerSizeSpectrum*2.*0.6, colorPHOS , colorPHOS);
	graphRatioCombPHOSPi0PbPb0510->Draw("p,same,e1");

	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb0510 , markerStyleConv,markerSizeSpectrum*2.*0.6, colorConv , colorConv);
	graphRatioCombConvPi0PbPb0510->Draw("p,same,e1");
	DrawGammaLines(0., 20.,1., 1.,1.,kBlack,2);

	TLatex *labelRatioIndMeasPbPb0510 = new TLatex(0.05,0.88,collisionSystemCent0510.Data());
	SetStyleTLatex( labelRatioIndMeasPbPb0510, 0.85*textsizeLabelsRatioUp,4);
	labelRatioIndMeasPbPb0510->Draw();


	pad6PartRatioIndMeas3->Update();
	pad6PartRatioIndMeas4->cd();
	pad6PartRatioIndMeas4->SetLogx();
	ratio2DInvXSectionOnlyPi0Down->DrawCopy(); 
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb4060Sys, markerStylePHOS,markerSizeSpectrum*0.6, colorPHOS , colorPHOS, widthLinesBoxes, kTRUE);
	graphRatioCombPHOSPi0PbPb4060Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb4060Sys, markerStyleConv,markerSizeSpectrum*0.6, colorConv , colorConv, widthLinesBoxes, kTRUE);
	graphRatioCombConvPi0PbPb4060Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb4060, markerStylePHOS,markerSizeSpectrum*2.*0.6, colorPHOS , colorPHOS);
	graphRatioCombPHOSPi0PbPb4060->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb4060 , markerStyleConv,markerSizeSpectrum*2.*0.6, colorConv , colorConv);
	graphRatioCombConvPi0PbPb4060->Draw("p,same,e1");
	DrawGammaLines(0., 20.,1., 1.,1.,kBlack,2);

	TLatex *labelRatioIndMeasPbPb04060 = new TLatex(0.05,0.9,collisionSystemCent4060.Data());
	SetStyleTLatex( labelRatioIndMeasPbPb04060, 0.85*textsizeLabelsRatioDown,4);	
	labelRatioIndMeasPbPb04060->Draw();

	// 	TLatex *labelRatioIndMeasLabelPbPb4060 = new TLatex(0.15,0.16,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelRatioIndMeasLabelPbPb4060, 0.038,4);
	// 	labelRatioIndMeasLabelPbPb4060->Draw();


	pad6PartRatioIndMeas5->cd();
	pad6PartRatioIndMeas5->SetLogx();
	ratio2DInvXSectionOnlyPi0Up->DrawCopy(); 
		
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb1020Sys, markerStylePHOS,markerSizeSpectrum*0.6, colorPHOS , colorPHOS, widthLinesBoxes, kTRUE);
	graphRatioCombPHOSPi0PbPb1020Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb1020Sys, markerStyleConv,markerSizeSpectrum*0.6, colorConv , colorConv, widthLinesBoxes, kTRUE);
	graphRatioCombConvPi0PbPb1020Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb1020, markerStylePHOS,markerSizeSpectrum*2.*0.6, colorPHOS , colorPHOS);
	graphRatioCombPHOSPi0PbPb1020->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb1020 , markerStyleConv,markerSizeSpectrum*2.*0.6, colorConv , colorConv);
	graphRatioCombConvPi0PbPb1020->Draw("p,same,e1");

	DrawGammaLines(0., 20.,1., 1.,1.,kBlack,2);
	TLatex *labelRatioIndMeasPbPb1020 = new TLatex(0.05,0.88,collisionSystemCent1020.Data());
	SetStyleTLatex( labelRatioIndMeasPbPb1020, 0.85*textsizeLabelsRatioUp,4);
	labelRatioIndMeasPbPb1020->Draw();
	// 	TLatex *labelRatioIndMeasLabelPbPb1020 = new TLatex(0.15,0.87,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelRatioIndMeasLabelPbPb1020, 0.04,4);
	// 	labelRatioIndMeasLabelPbPb1020->Draw();


	pad6PartRatioIndMeas5->Update();

	pad6PartRatioIndMeas6->Update();

	pad6PartRatioIndMeas6->cd();
	pad6PartRatioIndMeas6->SetLogx();
	ratio2DInvXSectionOnlyPi0Down->DrawCopy(); 
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb6080Sys, markerStylePHOS,markerSizeSpectrum*0.6, colorPHOS , colorPHOS, widthLinesBoxes, kTRUE);
	graphRatioCombPHOSPi0PbPb6080Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb6080Sys, markerStyleConv,markerSizeSpectrum*0.6, colorConv , colorConv, widthLinesBoxes, kTRUE);
	graphRatioCombConvPi0PbPb6080Sys->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSPi0PbPb6080, markerStylePHOS,markerSizeSpectrum*2.*0.6, colorPHOS , colorPHOS);
	graphRatioCombPHOSPi0PbPb6080->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRatioCombConvPi0PbPb6080 , markerStyleConv,markerSizeSpectrum*2.*0.6, colorConv , colorConv);
	graphRatioCombConvPi0PbPb6080->Draw("p,same,e1");
	DrawGammaLines(0., 20.,1., 1.,1.,kBlack,2);

	TLatex *labelRatioIndMeasPbPb6080 = new TLatex(0.05,0.9,collisionSystemCent6080.Data());
	SetStyleTLatex( labelRatioIndMeasPbPb6080, 0.85*textsizeLabelsRatioDown,4);
	labelRatioIndMeasPbPb6080->Draw();
	// 	TLatex *labelRatioIndMeasLabelPbPb6080 = new TLatex(0.03,0.16,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelRatioIndMeasLabelPbPb6080, 0.038,4);
	// 	labelRatioIndMeasLabelPbPb6080->Draw();

	pad6PartRatioIndMeas6->Update();

	canvas6PartRatioIndMeas->Update();	
	canvas6PartRatioIndMeas->SaveAs(Form("%s/InvYieldOnlyRatio_All6Parted_Paper_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
	delete pad6PartRatioIndMeas1;	
	delete pad6PartRatioIndMeas2;	
	delete pad6PartRatioIndMeas3;	
	delete pad6PartRatioIndMeas4;	
	delete pad6PartRatioIndMeas5;	
	delete pad6PartRatioIndMeas6;	
	delete canvas6PartRatioIndMeas;	

	TCanvas* canvasRAAAllTogether = new TCanvas("canvasRAAAllTogether","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasRAAAllTogether,  0.1, 0.01, 0.015, 0.1);

	Int_t textSizeLabelsPixelRAA = 50;
	Double_t marginRAA = 0.14*1200;
	Double_t textsizeLabelsRAA = 0;
	Double_t textsizeFacRAA = 0;

	if (canvasRAAAllTogether->XtoPixel(canvasRAAAllTogether->GetX2()) < canvasRAAAllTogether->YtoPixel(canvasRAAAllTogether->GetY1())){
		textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAAAllTogether->XtoPixel(canvasRAAAllTogether->GetX2()) ;
		textsizeFacRAA = (Double_t)1./canvasRAAAllTogether->XtoPixel(canvasRAAAllTogether->GetX2()) ;
	} else {
		textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAAAllTogether->YtoPixel(canvasRAAAllTogether->GetY1());
		textsizeFacRAA = (Double_t)1./canvasRAAAllTogether->YtoPixel(canvasRAAAllTogether->GetY1());
	}


	TH2F * histo2DRAAAll;
	histo2DRAAAll = new TH2F("histo2DRAAAll","histo2DRAAAll",1000,0.,15.,1000,0.0,1.3);
	SetStyleHistoTH2ForGraphs(histo2DRAAAll, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.9,0.9, 512, 510); 
	histo2DRAAAll->GetXaxis()->SetRangeUser(0.,12.5);
	histo2DRAAAll->GetXaxis()->SetLabelFont(42);
	histo2DRAAAll->GetYaxis()->SetLabelFont(42);
	histo2DRAAAll->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005, widthLinesBoxes, kTRUE);//, colorComb0005-5);
	graphRAASysCombInd0005->Draw("E2same");

	// 	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd0510, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510, widthLinesBoxes, kTRUE);//, colorComb0510-5);
	// 	graphRAASysCombInd0510->Draw("E2same");

	// 	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd1020, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020, widthLinesBoxes, kTRUE);//, colorComb1020-5);
	// 	graphRAASysCombInd1020->Draw("E2same");

	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);//, colorComb2040-5);
	graphRAASysCombInd2040->Draw("E2same");
	// 	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060, widthLinesBoxes, kTRUE);//, colorComb4060-5);
	// 	graphRAASysCombInd4060->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080, widthLinesBoxes, kTRUE);//, colorComb6080-5);
	graphRAASysCombInd6080->Draw("E2same");

	DrawGammaSetMarkerTGraphAsym(graphRAACombInd0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005);
	graphRAACombInd0005->Draw("p,same,e1");
	// 	DrawGammaSetMarkerTGraphAsym(graphRAACombInd0510, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510);
	// 	graphRAACombInd0510->Draw("p,same,e1");
	// 	DrawGammaSetMarkerTGraphAsym(graphRAACombInd1020, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020);
	// 	graphRAACombInd1020->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
	graphRAACombInd2040->Draw("p,same,e1");
	// 	DrawGammaSetMarkerTGraphAsym(graphRAACombInd4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
	// 	graphRAACombInd4060->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
	graphRAACombInd6080->Draw("p,same,e1");

	DrawGammaLines(0., 12.5,1., 1.,0.1,kGray+2);

	TLegend* legendRAAAll = new TLegend(0.61,0.8,0.96,0.96);
	legendRAAAll->SetFillColor(0);
	legendRAAAll->SetLineColor(0);
	legendRAAAll->SetTextFont(42);
	legendRAAAll->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAAAll->AddEntry(graphRAASysCombInd0005,"  0 -   5% #pi^{0}","pf");
	legendRAAAll->AddEntry(graphRAASysCombInd2040,"20 - 40% #pi^{0}","pf");
	legendRAAAll->AddEntry(graphRAASysCombInd6080,"60 - 80% #pi^{0}","pf");
	legendRAAAll->Draw();

	TBox* boxErrorNorm0005All = CreateBoxConv(colorComb0005Box, 0.4, 1.-normErr0005 , 0.6, 1.+normErr0005);
	boxErrorNorm0005All->Draw();

	TBox* boxErrorNorm2040All = CreateBoxConv(colorComb2040Box, 0.65, 1.-normErr2040 , 0.85, 1.+normErr2040);
	boxErrorNorm2040All->Draw();

	TBox* boxErrorNorm6080All = CreateBoxConv(colorComb6080Box, 0.9, 1.-normErr6080 , 1.1, 1.+normErr6080);
	boxErrorNorm6080All->Draw();
		

	TLatex *labelRAAPi0All = new TLatex(0.15,0.92,"Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV");
	SetStyleTLatex( labelRAAPi0All, 0.85*textsizeLabelsRAA,4);
	labelRAAPi0All->Draw();


	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);


	canvasRAAAllTogether->Update();
	canvasRAAAllTogether->Print(Form("%s/RAA_All_Paper.%s",outputDir.Data(),suffix.Data()));


	cout << "here ............................." << endl;


	TCanvas* canvasSpectraAllChargesTogether = new TCanvas("canvasSpectraAllChargesTogether","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasSpectraAllChargesTogether,  0.1, 0.01, 0.015, 0.08);

	canvasSpectraAllChargesTogether->SetLogy();
	canvasSpectraAllChargesTogether->SetLogx();
	TH2F * histo2DSpectraAllCharges;
	histo2DSpectraAllCharges = new TH2F("histo2DSpectraAllCharges","histo2DSpectraAllCharges",1000,0.1,20.,1000,1e-8,7e3	);
	SetStyleHistoTH2ForGraphs(histo2DSpectraAllCharges, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{#it{N}_{evt}^{AA}}#frac{dN^{AA}}{d#it{p}_{T} d#it{y}}}{ #frac{1}{#it{N}_{evt}^{pp}}#frac{dN^{pp}}{d#it{p}_{T} d#it{y}}}
	histo2DSpectraAllCharges->GetXaxis()->SetLabelOffset(-0.01);
	histo2DSpectraAllCharges->GetYaxis()->SetLabelOffset(0.01);
	histo2DSpectraAllCharges->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb1020FullSysErr, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, kRed ,kRed, widthLinesBoxes, kTRUE);//, colorComb1020-5);
	graphYieldPi0CombPbPb1020FullSysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb1020StatErr, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, kRed ,kRed);
	graphYieldPi0CombPbPb1020StatErr->Draw("p,same,e1");
	// 	DrawGammaSetMarkerTF1( fitYieldPi0PbPb1020, styleFitCommonSpectrum, widthCommonFit, colorComb1020);
	// 	fitYieldPi0PbPb1020->Draw("same");

	DrawGammaSetMarker(histoChargedPionSpecHighPtStat1020 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	DrawGammaSetMarker(histoChargedPionSpecHighPtSyst1020 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	histoChargedPionSpecHighPtStat1020->Draw("p,same,e1");
	histoChargedPionSpecHighPtSyst1020->Draw("E[]same");

	DrawGammaSetMarker(histoChargedPionSpecLowPtStat1020 ,24,markerSizeChargedHadronSpectrum, kBlue , kBlue);					 
	DrawGammaSetMarker(histoChargedPionSpecLowPtStat1020 ,24,markerSizeChargedHadronSpectrum, kBlue , kBlue);					 
	histoChargedPionSpecLowPtStat1020->Draw("p,same,e1");
	histoChargedPionSpecLowPtStat1020->Draw("E[]same");



	canvasSpectraAllChargesTogether->Update();
	canvasSpectraAllChargesTogether->Print(Form("%s/Spectra_AllCharges_Paper.%s",outputDir.Data(),suffix.Data()));

	TH1D* histoFitYieldPi0PbPb6080 = (TH1D*)fitYieldPi0PbPb6080->GetHistogram();
	TH1D* histoFitQCDYieldPi0PbPb6080 = (TH1D*)fitYieldDataQCDPi0PbPb6080->GetHistogram();
	TH1D* histoFitQCDYieldPi0PbPb6080Draw = (TH1D*)histoFitQCDYieldPi0PbPb6080->Clone("histoFitQCDYieldPi0PbPb6080Draw");

	graphYieldPi0CombPbPb4060FullSysErr = ScaleGraph(graphYieldPi0CombPbPb4060FullSysErr,2);   
	graphYieldPi0CombPbPb4060StatErr = ScaleGraph(graphYieldPi0CombPbPb4060StatErr,2);
	TH1D* histoFitYieldPi0PbPb4060 = (TH1D*)fitYieldPi0PbPb4060->GetHistogram();
	histoFitYieldPi0PbPb4060->Scale(2);
	TH1D* histoFitQCDYieldPi0PbPb4060 = (TH1D*)fitYieldDataQCDPi0PbPb4060->GetHistogram();
	histoFitQCDYieldPi0PbPb4060->Scale(2);

	graphYieldPi0CombPbPb2040FullSysErr = ScaleGraph(graphYieldPi0CombPbPb2040FullSysErr,4);   
	graphYieldPi0CombPbPb2040StatErr = ScaleGraph(graphYieldPi0CombPbPb2040StatErr,4);
	TH1D* histoFitYieldPi0PbPb2040 = (TH1D*)fitYieldPi0PbPb2040->GetHistogram();
	histoFitYieldPi0PbPb2040->Scale(4);
	TH1D* histoFitQCDYieldPi0PbPb2040 = (TH1D*)fitYieldDataQCDPi0PbPb2040->GetHistogram();
	histoFitQCDYieldPi0PbPb2040->Scale(4);

	graphYieldPi0CombPbPb1020FullSysErr = ScaleGraph(graphYieldPi0CombPbPb1020FullSysErr,8);   
	graphYieldPi0CombPbPb1020StatErr = ScaleGraph(graphYieldPi0CombPbPb1020StatErr,8);
	TH1D* histoFitYieldPi0PbPb1020 = (TH1D*)fitYieldPi0PbPb1020->GetHistogram();
	histoFitYieldPi0PbPb1020->Scale(8);
	TH1D* histoFitQCDYieldPi0PbPb1020 = (TH1D*)fitYieldDataQCDPi0PbPb1020->GetHistogram();
	histoFitQCDYieldPi0PbPb1020->Scale(8);

	graphYieldPi0CombPbPb0510FullSysErr = ScaleGraph(graphYieldPi0CombPbPb0510FullSysErr,32);   
	graphYieldPi0CombPbPb0510StatErr = ScaleGraph(graphYieldPi0CombPbPb0510StatErr,32);
	TH1D* histoFitYieldPi0PbPb0510 = (TH1D*)fitYieldPi0PbPb0510->GetHistogram();
	histoFitYieldPi0PbPb0510->Scale(32);
	TH1D* histoFitQCDYieldPi0PbPb0510 = (TH1D*)fitYieldDataQCDPi0PbPb0510->GetHistogram();
	histoFitQCDYieldPi0PbPb0510->Scale(32);

	graphYieldPi0CombPbPb0005FullSysErr = ScaleGraph(graphYieldPi0CombPbPb0005FullSysErr,128);
	graphYieldPi0CombPbPb0005StatErr = ScaleGraph(graphYieldPi0CombPbPb0005StatErr,128);
	TH1D* histoFitYieldPi0PbPb0005 = (TH1D*)fitYieldPi0PbPb0005->GetHistogram();
	histoFitYieldPi0PbPb0005->Scale(128);
	TH1D* histoFitQCDYieldPi0PbPb0005 = (TH1D*)fitYieldDataQCDPi0PbPb0005->GetHistogram();
	histoFitQCDYieldPi0PbPb0005->Scale(128);




	TCanvas* canvasSpectraAllTogether = new TCanvas("canvasSpectraAllTogether","",200,10,1200,1700);  // gives the page size
	DrawGammaCanvasSettings( canvasSpectraAllTogether,  0.17, 0.005, 0.005, 0.08);

	Int_t textSizeLabelsPixelSpectra = 50;
	Double_t marginSpectra = 0.14*1200;
	Double_t textsizeLabelsSpectra = 0;
	Double_t textsizeFacSpectra = 0;

	if (canvasSpectraAllTogether->XtoPixel(canvasSpectraAllTogether->GetX2()) < canvasSpectraAllTogether->YtoPixel(canvasSpectraAllTogether->GetY1())){
		textsizeLabelsSpectra = (Double_t)textSizeLabelsPixelSpectra/canvasSpectraAllTogether->XtoPixel(canvasSpectraAllTogether->GetX2()) ;
		textsizeFacSpectra = (Double_t)1./canvasSpectraAllTogether->XtoPixel(canvasSpectraAllTogether->GetX2()) ;
	} else {
		textsizeLabelsSpectra = (Double_t)textSizeLabelsPixelSpectra/canvasSpectraAllTogether->YtoPixel(canvasSpectraAllTogether->GetY1());
		textsizeFacSpectra = (Double_t)1./canvasSpectraAllTogether->YtoPixel(canvasSpectraAllTogether->GetY1());
	}


	canvasSpectraAllTogether->SetLogy();
	canvasSpectraAllTogether->SetLogx();
	TH2F * histo2DSpectraAll;
	histo2DSpectraAll = new TH2F("histo2DSpectraAll","histo2DSpectraAll",1000,0.09,20.,1000,1e-9,7e4	);
	SetStyleHistoTH2ForGraphs(histo2DSpectraAll, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.85*textsizeLabelsSpectra,textsizeLabelsSpectra, 0.85*textsizeLabelsSpectra,textsizeLabelsSpectra, 0.87,1.8);// 512, 505); //#frac{#frac{1}{#it{N}_{evt}^{AA}}#frac{dN^{AA}}{d#it{p}_{T} d#it{y}}}{ #frac{1}{#it{N}_{evt}^{pp}}#frac{dN^{pp}}{d#it{p}_{T} d#it{y}}}
	histo2DSpectraAll->GetXaxis()->SetLabelFont(42);
	histo2DSpectraAll->GetYaxis()->SetLabelFont(42); 
	histo2DSpectraAll->GetXaxis()->SetTitleFont(62);
	histo2DSpectraAll->GetYaxis()->SetTitleFont(62);
	histo2DSpectraAll->GetXaxis()->SetLabelOffset(-0.015);
	histo2DSpectraAll->GetYaxis()->SetLabelOffset(0.01);
	histo2DSpectraAll->GetXaxis()->SetRangeUser(0.33,15);
	histo2DSpectraAll->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysPi02760GeVPlot, markerStyleCommmonSpectrumpp,markerSizeSpectrum, kBlack , kBlack, widthLinesBoxes, kTRUE);//, colorComb1020-5);
	graphInvSectionCombSysPi02760GeVPlot->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0005FullSysErr, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005, widthLinesBoxes, kTRUE);//, colorComb0005-5);
	graphYieldPi0CombPbPb0005FullSysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0510FullSysErr, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510, widthLinesBoxes, kTRUE);//, colorComb0510-5);
	graphYieldPi0CombPbPb0510FullSysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb1020FullSysErr, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020, widthLinesBoxes, kTRUE);//, colorComb1020-5);
	graphYieldPi0CombPbPb1020FullSysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb2040FullSysErr, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);//, colorComb2040-5);
	graphYieldPi0CombPbPb2040FullSysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb4060FullSysErr, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060, widthLinesBoxes, kTRUE);//, colorComb4060-5);
	graphYieldPi0CombPbPb4060FullSysErr->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb6080FullSysErr, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080, widthLinesBoxes, kTRUE);//, colorComb6080-5);
	graphYieldPi0CombPbPb6080FullSysErr->Draw("E2same");

	DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatPi02760GeVPlot, markerStyleCommmonSpectrumpp,markerSizeSpectrum, kBlack , kBlack);
	graphInvSectionCombStatPi02760GeVPlot->Draw("p,same,e1");
	fitInvCrossSectionPi0Comb2760GeVPow->SetRange(3., maxPtPi0PbPb6080);
	DrawGammaSetMarkerTF1( fitInvCrossSectionPi0Comb2760GeVPow, styleFitCommonSpectrum, widthCommonFit, kBlack);
	fitInvCrossSectionPi0Comb2760GeVPow->Draw("same");
	fitInvCrossSectionPi0Comb2760GeV->SetRange(0.4, 15.);
	DrawGammaSetMarkerTF1( fitInvCrossSectionPi0Comb2760GeV, 7, widthCommonFit, kGray+2);
	fitInvCrossSectionPi0Comb2760GeV->Draw("same");

	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0005StatErr, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005);
	graphYieldPi0CombPbPb0005StatErr->Draw("p,same,e1");
	// 	SetStyleHisto(histoFitYieldPi0PbPb0005, widthCommonFit, styleFitCommonSpectrum, colorComb0005);
	// 	histoFitYieldPi0PbPb0005->Draw("same,c");
	SetStyleHisto(histoFitQCDYieldPi0PbPb0005, widthCommonFit, 5, colorComb0005);
	histoFitQCDYieldPi0PbPb0005->Draw("same,c");

	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb0510StatErr, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510);
	graphYieldPi0CombPbPb0510StatErr->Draw("p,same,e1");
	// 	SetStyleHisto(histoFitYieldPi0PbPb0510, widthCommonFit, styleFitCommonSpectrum, colorComb0510);
	// 	histoFitYieldPi0PbPb0510->Draw("same,c");
	SetStyleHisto(histoFitQCDYieldPi0PbPb0510, widthCommonFit, 5, colorComb0510);
	histoFitQCDYieldPi0PbPb0510->Draw("same,c");


	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb1020StatErr, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020);
	graphYieldPi0CombPbPb1020StatErr->Draw("p,same,e1");
	// 	DrawGammaSetMarkerTF1( fitYieldPi0PbPb1020, styleFitCommonSpectrum, widthCommonFit, colorComb1020);
	// 	fitYieldPi0PbPb1020->Draw("same");
	//    DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb1020, 1, widthCommonFit, colorComb1020);
	//    fitYieldDataQCDPi0PbPb1020->Draw("same");
	SetStyleHisto(histoFitQCDYieldPi0PbPb1020, widthCommonFit, 5, colorComb1020);
	histoFitQCDYieldPi0PbPb1020->Draw("same,c");


	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb2040StatErr, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
	graphYieldPi0CombPbPb2040StatErr->Draw("p,same,e1");
	// 	DrawGammaSetMarkerTF1( fitYieldPi0PbPb2040, styleFitCommonSpectrum, widthCommonFit, colorComb2040);
	// 	fitYieldPi0PbPb2040->Draw("same");
	//    DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb2040, 1, widthCommonFit, colorComb2040);
	//    fitYieldDataQCDPi0PbPb2040->Draw("same");
	SetStyleHisto(histoFitQCDYieldPi0PbPb2040, widthCommonFit, 5, colorComb2040);
	histoFitQCDYieldPi0PbPb2040->Draw("same,c");   

	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb4060StatErr, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
	graphYieldPi0CombPbPb4060StatErr->Draw("p,same,e1");
	// 	DrawGammaSetMarkerTF1( fitYieldPi0PbPb4060, styleFitCommonSpectrum, widthCommonFit, colorComb4060);
	// 	fitYieldPi0PbPb4060->Draw("same");
	SetStyleHisto(histoFitQCDYieldPi0PbPb4060, widthCommonFit, 5, colorComb4060);
	histoFitQCDYieldPi0PbPb4060->Draw("same,c");   
	//    DrawGammaSetMarkerTF1(fitYieldPi0PbPb4060, 1, widthCommonFit, colorComb4060);
	//    fitYieldPi0PbPb4060->Draw("same");

	DrawGammaSetMarkerTGraphAsym(graphYieldPi0CombPbPb6080StatErr, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
	graphYieldPi0CombPbPb6080StatErr->Draw("p,same,e1");
	// 	DrawGammaSetMarkerTF1( fitYieldPi0PbPb6080, styleFitCommonSpectrum, widthCommonFit, colorComb6080);
	// 	fitYieldPi0PbPb6080->Draw("same");
	DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPb6080, 5, widthCommonFit, colorComb6080);
	fitYieldDataQCDPi0PbPb6080->Draw("same");

	SetStyleHisto(histoFitQCDYieldPi0PbPb6080Draw, widthCommonFit, 5, kGray+2);
	histoFitQCDYieldPi0PbPb6080Draw->SetLineColor(kGray+2);
	histoFitQCDYieldPi0PbPb6080Draw->SetMarkerColor(kGray+2);
		
		
	// 	TLatex *labelSpectraPi0LabelPbPb = new TLatex(0.175,0.34,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelSpectraPi0LabelPbPb, 0.03,4);
	// 	labelSpectraPi0LabelPbPb->Draw();

	//    TLatex *labelSpectraPbPb = new TLatex(0.16,0.34,"Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV");
	//    SetStyleTLatex( labelSpectraPbPb, 0.03,4);
	//    labelSpectraPbPb->Draw();

	TLegend* legendSpectra = new TLegend(0.2,0.09,0.54,0.31);
	legendSpectra->SetFillColor(0);
	legendSpectra->SetLineColor(0);
	legendSpectra->SetTextSize(0.85*textsizeLabelsSpectra);
	legendSpectra->SetTextFont(42);
	legendSpectra->SetMargin(0.14);
	legendSpectra->AddEntry(graphYieldPi0CombPbPb0005FullSysErr,Form("%s #times 2^{7}",collisionSystemCent0005Legend.Data()),"pf");
	legendSpectra->AddEntry(graphYieldPi0CombPbPb0510FullSysErr,Form("%s #times 2^{5}",collisionSystemCent0510Legend.Data()),"pf");
	legendSpectra->AddEntry(graphYieldPi0CombPbPb1020FullSysErr,Form("%s #times 2^{3}",collisionSystemCent1020.Data()),"pf");
	legendSpectra->AddEntry(graphYieldPi0CombPbPb2040FullSysErr,Form("%s #times 2^{2}",collisionSystemCent2040.Data()),"pf");
	legendSpectra->AddEntry(graphYieldPi0CombPbPb4060FullSysErr,Form("%s #times 2^{1}",collisionSystemCent4060.Data()),"pf");
	legendSpectra->AddEntry(graphYieldPi0CombPbPb6080FullSysErr,Form("%s #times 2^{0}",collisionSystemCent6080.Data()),"pf");
	legendSpectra->AddEntry(histoFitQCDYieldPi0PbPb6080Draw,"fits to Pb-Pb","l");
	legendSpectra->Draw();

	TLegend* legendSpectraPP = new TLegend(0.65,0.85,0.95,0.95);
	legendSpectraPP->SetFillColor(0);
	legendSpectraPP->SetLineColor(0);
	legendSpectraPP->SetTextSize(0.85*textsizeLabelsSpectra);
	legendSpectraPP->SetMargin(0.16);
	legendSpectraPP->SetTextFont(42);
	legendSpectraPP->AddEntry(graphInvSectionCombSysPi02760GeVPlot,collisionSystemPP.Data(),"pf");
	legendSpectraPP->AddEntry(fitInvCrossSectionPi0Comb2760GeV,"Tsallis fit","l");
	legendSpectraPP->AddEntry(fitInvCrossSectionPi0Comb2760GeVPow,"power law fit","l");
	legendSpectraPP->Draw();

	canvasSpectraAllTogether->Update();
	canvasSpectraAllTogether->Print(Form("%s/Spectra_All_Paper.%s",outputDir.Data(),suffix.Data()));

	canvasSpectraAllTogether->cd();
	canvasSpectraAllTogether->SetLogy();
	canvasSpectraAllTogether->SetLogx(0);
	TH2F * histo2DSpectraAll_linX;
	histo2DSpectraAll_linX = new TH2F("histo2DSpectraAll_linX","histo2DSpectraAll_linX",1000,0.0,12.5,1000,1e-8,7e4	);
	SetStyleHistoTH2ForGraphs(histo2DSpectraAll_linX, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.85*textsizeLabelsSpectra,textsizeLabelsSpectra, 0.85*textsizeLabelsSpectra,textsizeLabelsSpectra, 0.87,1.8);// 512, 505); //#frac{#frac{1}{#it{N}_{evt}^{AA}}#frac{dN^{AA}}{d#it{p}_{T} d#it{y}}}{ #frac{1}{#it{N}_{evt}^{pp}}#frac{dN^{pp}}{d#it{p}_{T} d#it{y}}}
	histo2DSpectraAll_linX->GetXaxis()->SetLabelFont(42);
	histo2DSpectraAll_linX->GetYaxis()->SetLabelFont(42); 
	histo2DSpectraAll_linX->GetXaxis()->SetTitleFont(62);
	histo2DSpectraAll_linX->GetYaxis()->SetTitleFont(62);
	//  	histo2DSpectraAll_linX->GetXaxis()->SetLabelOffset(-0.015);
	histo2DSpectraAll_linX->GetYaxis()->SetLabelOffset(0.01);
	// 	histo2DSpectraAll_linX->GetXaxis()->SetRangeUser(0.33,15);
	histo2DSpectraAll_linX->DrawCopy(); 

	graphInvSectionCombSysPi02760GeVPlot->Draw("E2same");
	graphYieldPi0CombPbPb0005FullSysErr->Draw("E2same");
	graphYieldPi0CombPbPb0510FullSysErr->Draw("E2same");
	graphYieldPi0CombPbPb1020FullSysErr->Draw("E2same");
	graphYieldPi0CombPbPb2040FullSysErr->Draw("E2same");
	graphYieldPi0CombPbPb4060FullSysErr->Draw("E2same");
	graphYieldPi0CombPbPb6080FullSysErr->Draw("E2same");

	graphInvSectionCombStatPi02760GeVPlot->Draw("p,same,e1");
	fitInvCrossSectionPi0Comb2760GeVPow->Draw("same");
	fitInvCrossSectionPi0Comb2760GeV->Draw("same");

	graphYieldPi0CombPbPb0005StatErr->Draw("p,same,e1");
	histoFitQCDYieldPi0PbPb0005->Draw("same,c");	
	graphYieldPi0CombPbPb0510StatErr->Draw("p,same,e1");
	histoFitQCDYieldPi0PbPb0510->Draw("same,c");
	graphYieldPi0CombPbPb1020StatErr->Draw("p,same,e1");
	histoFitQCDYieldPi0PbPb1020->Draw("same,c");
	graphYieldPi0CombPbPb2040StatErr->Draw("p,same,e1");
	histoFitQCDYieldPi0PbPb2040->Draw("same,c");   
	graphYieldPi0CombPbPb4060StatErr->Draw("p,same,e1");
	histoFitYieldPi0PbPb4060->Draw("same,c");   
	graphYieldPi0CombPbPb6080StatErr->Draw("p,same,e1");
	fitYieldPi0PbPb6080->Draw("same");

	TLegend* legendSpectra_lin = new TLegend(0.38,0.75,0.73,0.97);
	legendSpectra_lin->SetFillColor(0);
	legendSpectra_lin->SetLineColor(0);
	legendSpectra_lin->SetTextSize(0.85*textsizeLabelsSpectra);
	legendSpectra_lin->SetTextFont(42);
	legendSpectra_lin->SetMargin(0.14);
	legendSpectra_lin->AddEntry(graphYieldPi0CombPbPb0005FullSysErr,Form("%s #times 2^{7}",collisionSystemCent0005Legend.Data()),"pf");
	legendSpectra_lin->AddEntry(graphYieldPi0CombPbPb0510FullSysErr,Form("%s #times 2^{5}",collisionSystemCent0510Legend.Data()),"pf");
	legendSpectra_lin->AddEntry(graphYieldPi0CombPbPb1020FullSysErr,Form("%s #times 2^{3}",collisionSystemCent1020.Data()),"pf");
	legendSpectra_lin->AddEntry(graphYieldPi0CombPbPb2040FullSysErr,Form("%s #times 2^{2}",collisionSystemCent2040.Data()),"pf");
	legendSpectra_lin->AddEntry(graphYieldPi0CombPbPb4060FullSysErr,Form("%s #times 2^{1}",collisionSystemCent4060.Data()),"pf");
	legendSpectra_lin->AddEntry(graphYieldPi0CombPbPb6080FullSysErr,Form("%s #times 2^{0}",collisionSystemCent6080.Data()),"pf");
	legendSpectra_lin->AddEntry(histoFitYieldPi0PbPb6080,"fits to Pb-Pb","l");
	legendSpectra_lin->Draw();

	TLegend* legendSpectraPP_lin = new TLegend(0.21,0.11,0.56,0.21);
	legendSpectraPP_lin->SetFillColor(0);
	legendSpectraPP_lin->SetLineColor(0);
	legendSpectraPP_lin->SetTextSize(0.85*textsizeLabelsSpectra);
	legendSpectraPP_lin->SetMargin(0.14);
	legendSpectraPP_lin->SetTextFont(42);
	legendSpectraPP_lin->AddEntry(graphInvSectionCombSysPi02760GeVPlot,collisionSystemPP.Data(),"pf");
	legendSpectraPP_lin->AddEntry(fitInvCrossSectionPi0Comb2760GeV,"Tsallis fit","l");
	legendSpectraPP_lin->AddEntry(fitInvCrossSectionPi0Comb2760GeVPow,"power law fit","l");
	legendSpectraPP_lin->Draw();

	canvasSpectraAllTogether->Update();
	canvasSpectraAllTogether->Print(Form("%s/Spectra_All_Paper_linX.%s",outputDir.Data(),suffix.Data()));


	gEPOS_RAA_6080 =   ScaleGraph(gEPOS_Spec_6080,1./nColl6080);
	gEPOS_RAA_6080 =  CalculateGraphRatioToFit (gEPOS_RAA_6080, fitInvCrossSectionPi0Comb2760GeV);
	gEPOS_RAA_4060 =   ScaleGraph(gEPOS_Spec_4060,1./nColl4060);
	gEPOS_RAA_4060 =  CalculateGraphRatioToFit (gEPOS_RAA_4060, fitInvCrossSectionPi0Comb2760GeV);
	gEPOS_RAA_2040 =   ScaleGraph(gEPOS_Spec_2040,1./nColl2040);
	gEPOS_RAA_2040 =  CalculateGraphRatioToFit (gEPOS_RAA_2040, fitInvCrossSectionPi0Comb2760GeV);
	gEPOS_RAA_1020 =   ScaleGraph(gEPOS_Spec_1020,1./nColl1020);
	gEPOS_RAA_1020 =  CalculateGraphRatioToFit (gEPOS_RAA_1020, fitInvCrossSectionPi0Comb2760GeV);
	gEPOS_RAA_0510 =   ScaleGraph(gEPOS_Spec_0510,1./nColl0510);
	gEPOS_RAA_0510 =  CalculateGraphRatioToFit (gEPOS_RAA_0510, fitInvCrossSectionPi0Comb2760GeV);
	gEPOS_RAA_0005 =   ScaleGraph(gEPOS_Spec_0005,1./nColl0005);
	gEPOS_RAA_0005 =  CalculateGraphRatioToFit (gEPOS_RAA_0005, fitInvCrossSectionPi0Comb2760GeV);
	gEPOS_RAA_0020 = (TGraphErrors*)gEPOS_RAA_0005->Clone("gEPOS_RAA_0020");


	gRatioEPOSToFit6080 =  CalculateGraphRatioToFit (gEPOS_Spec_6080, fitYieldDataQCDPi0PbPb6080);
	gRatioEPOSToFit4060 =  CalculateGraphRatioToFit (gEPOS_Spec_4060, fitYieldDataQCDPi0PbPb4060);
	gRatioEPOSToFit2040 =  CalculateGraphRatioToFit (gEPOS_Spec_2040, fitYieldDataQCDPi0PbPb2040);
	gRatioEPOSToFit1020 =  CalculateGraphRatioToFit (gEPOS_Spec_1020, fitYieldDataQCDPi0PbPb1020);
	gRatioEPOSToFit0510 =  CalculateGraphRatioToFit (gEPOS_Spec_0510, fitYieldDataQCDPi0PbPb0510);
	gRatioEPOSToFit0005 =  CalculateGraphRatioToFit (gEPOS_Spec_0005, fitYieldDataQCDPi0PbPb0005);


	gRatioKopeliovichELossToFit0005 = CalculateGraphRatioToFit( gKopeliovichELoss_Spec_0005, fitYieldDataQCDPi0PbPb0005);
	gRatioKopeliovichHydroToFit0005 = CalculateGraphRatioToFit( gKopeliovichHydro_Spec_0005, fitYieldDataQCDPi0PbPb0005);
	gRatioKopeliovichTotalToFit0005 = CalculateGraphRatioToFit( gKopeliovichTotal_Spec_0005, fitYieldDataQCDPi0PbPb0005);

	gRatioKopeliovichELossToFit2040 = CalculateGraphRatioToFit( gKopeliovichELoss_Spec_2040, fitYieldDataQCDPi0PbPb2040);
	gRatioKopeliovichHydroToFit2040 = CalculateGraphRatioToFit( gKopeliovichHydro_Spec_2040, fitYieldDataQCDPi0PbPb2040);
	gRatioKopeliovichTotalToFit2040 = CalculateGraphRatioToFit( gKopeliovichTotal_Spec_2040, fitYieldDataQCDPi0PbPb2040);

	cout << fitYieldDataQCDPi0PbPb6080->Eval(7.) << endl;
	cout << fitYieldDataQCDPi0PbPb6080->Eval(9.) << endl;
	gRatioKopeliovichELossToFit6080 = CalculateGraphRatioToFit( gKopeliovichELoss_Spec_6080, fitYieldDataQCDPi0PbPb6080);
	gRatioKopeliovichHydroToFit6080 = CalculateGraphRatioToFit( gKopeliovichHydro_Spec_6080, fitYieldDataQCDPi0PbPb6080);
	gRatioKopeliovichTotalToFit6080 = CalculateGraphRatioToFit( gKopeliovichTotal_Spec_6080, fitYieldDataQCDPi0PbPb6080);


	canvasSpectraAllTogether->cd();
	canvasSpectraAllTogether->SetLogy();
	canvasSpectraAllTogether->SetLogx();
		histo2DSpectraAll->GetXaxis()->SetRangeUser(0.29,20);
		histo2DSpectraAll->DrawCopy(); 
		
		graphYieldPi0CombPbPb0005FullSysErr->Draw("E2same");
		graphYieldPi0CombPbPb0510FullSysErr->Draw("E2same");
		graphYieldPi0CombPbPb1020FullSysErr->Draw("E2same");
		graphYieldPi0CombPbPb2040FullSysErr->Draw("E2same");
		graphYieldPi0CombPbPb4060FullSysErr->Draw("E2same");
		graphYieldPi0CombPbPb6080FullSysErr->Draw("E2same");

		graphYieldPi0CombPbPb0005StatErr->Draw("p,same,e1");   
		graphYieldPi0CombPbPb0510StatErr->Draw("p,same,e1");
		graphYieldPi0CombPbPb1020StatErr->Draw("p,same,e1");   
		graphYieldPi0CombPbPb2040StatErr->Draw("p,same,e1");
		graphYieldPi0CombPbPb4060StatErr->Draw("p,same,e1");
		graphYieldPi0CombPbPb6080StatErr->Draw("p,same,e1");

		
		for (Int_t i = 0; i < 1; i++){
			gEPOS_Spec_6080->RemovePoint(0);
			gEPOS_Spec_4060->RemovePoint(0);
			gEPOS_Spec_2040->RemovePoint(0);
			gEPOS_Spec_1020->RemovePoint(0);
			gEPOS_Spec_0510->RemovePoint(0);
			gEPOS_Spec_0005->RemovePoint(0);
			
		}   
		
		DrawGammaSetMarkerTGraph(gEPOS_Spec_6080, markerStyleCommmonSpectrum6080,markerSizeSpectrum, colorEPOS6080, colorEPOS6080,10,kTRUE);
		gEPOS_Spec_6080->Draw("3 same");
		gEPOS_Spec_4060 = ScaleGraph(gEPOS_Spec_4060,2);
		DrawGammaSetMarkerTGraph(gEPOS_Spec_4060, markerStyleCommmonSpectrum4060,markerSizeSpectrum, colorEPOS4060, colorEPOS4060,2);
		gEPOS_Spec_4060->Draw("3 same");
		gEPOS_Spec_2040 = ScaleGraph(gEPOS_Spec_2040,4);
		DrawGammaSetMarkerTGraph(gEPOS_Spec_2040, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorEPOS2040, colorEPOS2040,2);
		gEPOS_Spec_2040->Draw("3 same");
		gEPOS_Spec_1020 = ScaleGraph(gEPOS_Spec_1020,8);
		DrawGammaSetMarkerTGraph(gEPOS_Spec_1020, markerStyleCommmonSpectrum1020,markerSizeSpectrum, colorEPOS1020, colorEPOS1020,2);
		gEPOS_Spec_1020->Draw("3 same");
		gEPOS_Spec_0510 = ScaleGraph(gEPOS_Spec_0510,32);
		DrawGammaSetMarkerTGraph(gEPOS_Spec_0510, markerStyleCommmonSpectrum0510,markerSizeSpectrum, colorEPOS0510, colorEPOS0510,2);
		gEPOS_Spec_0510->Draw("3 same");
		gEPOS_Spec_0005 = ScaleGraph(gEPOS_Spec_0005,128);
		DrawGammaSetMarkerTGraph(gEPOS_Spec_0005, markerStyleCommmonSpectrum0005,markerSizeSpectrum, colorEPOS0005, colorEPOS0005,2);
		gEPOS_Spec_0005->Draw("3 same");

	//    legendSpectra->Draw();

	canvasSpectraAllTogether->Update();
	canvasSpectraAllTogether->Print(Form("%s/Spectra_All_TheoryEPOS_Paper.%s",outputDir.Data(),suffix.Data()));




	TCanvas* canvasRatioEPOSToData = new TCanvas("canvasRatioEPOSToData","",200,10,1300,3000);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioEPOSToData,  0.15, 0.02, 0.03, 0.06);

	TPad* padRatioEPOSToData0005 = new TPad("padRatioEPOSToData0005", "", 0., 0.84, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padRatioEPOSToData0005, 0.15, 0.02, 0.03, 0.);
	padRatioEPOSToData0005->Draw();

	TPad* padRatioEPOSToData0510 = new TPad("padRatioEPOSToData0510", "", 0., 0.68, 1., 0.84,-1, -1, -2);
	DrawGammaPadSettings( padRatioEPOSToData0510, 0.15, 0.02, 0.0, 0.);
	padRatioEPOSToData0510->Draw();

	TPad* padRatioEPOSToData1020 = new TPad("padRatioEPOSToData1020", "", 0., 0.53, 1., 0.68,-1, -1, -2);
	DrawGammaPadSettings( padRatioEPOSToData1020, 0.15, 0.02, 0., 0.);
	padRatioEPOSToData1020->Draw();

	TPad* padRatioEPOSToData2040 = new TPad("padRatioEPOSToData2040", "", 0., 0.37, 1., 0.53,-1, -1, -2);
	DrawGammaPadSettings( padRatioEPOSToData2040, 0.15, 0.02, 0., 0.);
	padRatioEPOSToData2040->Draw();

	TPad* padRatioEPOSToData4060 = new TPad("padRatioEPOSToData4060", "", 0., 0.21, 1., 0.37,-1, -1, -2);
	DrawGammaPadSettings( padRatioEPOSToData4060, 0.15, 0.02, 0., 0.);
	padRatioEPOSToData4060->Draw();

	TPad* padRatioEPOSToData6080 = new TPad("padRatioEPOSToData6080", "", 0., 0., 1., 0.21,-1, -1, -2);
	DrawGammaPadSettings( padRatioEPOSToData6080, 0.15, 0.02, 0., 0.2);
	padRatioEPOSToData6080->Draw();

	for (Int_t i = 0; i < 2; i++){
			gRatioEPOSToFit6080->RemovePoint(0);
			gRatioEPOSToFit4060->RemovePoint(0);
			gRatioEPOSToFit2040->RemovePoint(0);
			gRatioEPOSToFit1020->RemovePoint(0);
			gRatioEPOSToFit0510->RemovePoint(0);
			gRatioEPOSToFit0005->RemovePoint(0);
			
	}

	TGraphAsymmErrors*  graphPi0SysDataToFit0005 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb0005FullSysErrUnscaled->Clone("graphPi0SysDataToFit0005");
	graphPi0SysDataToFit0005 = CalculateGraphErrRatioToFit(graphPi0SysDataToFit0005, fitYieldDataQCDPi0PbPb0005);
	TGraphAsymmErrors*  graphPi0StatDataToFit0005 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb0005StatErrUnscaled->Clone("graphPi0StatDataToFit0005");
	graphPi0StatDataToFit0005 = CalculateGraphErrRatioToFit(graphPi0StatDataToFit0005, fitYieldDataQCDPi0PbPb0005);
	TGraphAsymmErrors*  graphPi0SysDataToFit0510 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb0510FullSysErrUnscaled->Clone("graphPi0SysDataToFit0510");
	graphPi0SysDataToFit0510 = CalculateGraphErrRatioToFit(graphPi0SysDataToFit0510, fitYieldDataQCDPi0PbPb0510);
	TGraphAsymmErrors*  graphPi0StatDataToFit0510 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb0510StatErrUnscaled->Clone("graphPi0StatDataToFit0510");
	graphPi0StatDataToFit0510 = CalculateGraphErrRatioToFit(graphPi0StatDataToFit0510, fitYieldDataQCDPi0PbPb0510);
	TGraphAsymmErrors*  graphPi0SysDataToFit1020 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb1020FullSysErrUnscaled->Clone("graphPi0SysDataToFit1020");
	graphPi0SysDataToFit1020 = CalculateGraphErrRatioToFit(graphPi0SysDataToFit1020, fitYieldDataQCDPi0PbPb1020);
	TGraphAsymmErrors*  graphPi0StatDataToFit1020 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb1020StatErrUnscaled->Clone("graphPi0StatDataToFit1020");
	graphPi0StatDataToFit1020 = CalculateGraphErrRatioToFit(graphPi0StatDataToFit1020, fitYieldDataQCDPi0PbPb1020);
	TGraphAsymmErrors*  graphPi0SysDataToFit2040 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb2040FullSysErrUnscaled->Clone("graphPi0SysDataToFit2040");
	graphPi0SysDataToFit2040 = CalculateGraphErrRatioToFit(graphPi0SysDataToFit2040, fitYieldDataQCDPi0PbPb2040);
	TGraphAsymmErrors*  graphPi0StatDataToFit2040 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb2040StatErrUnscaled->Clone("graphPi0StatDataToFit2040");
	graphPi0StatDataToFit2040 = CalculateGraphErrRatioToFit(graphPi0StatDataToFit2040, fitYieldDataQCDPi0PbPb2040);
	TGraphAsymmErrors*  graphPi0SysDataToFit4060 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb4060FullSysErrUnscaled->Clone("graphPi0SysDataToFit4060");
	graphPi0SysDataToFit4060 = CalculateGraphErrRatioToFit(graphPi0SysDataToFit4060, fitYieldPi0PbPb4060);
	TGraphAsymmErrors*  graphPi0StatDataToFit4060 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb4060StatErrUnscaled->Clone("graphPi0StatDataToFit4060");
	graphPi0StatDataToFit4060 = CalculateGraphErrRatioToFit(graphPi0StatDataToFit4060, fitYieldPi0PbPb4060);
	TGraphAsymmErrors*  graphPi0SysDataToFit6080 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb6080FullSysErrUnscaled->Clone("graphPi0SysDataToFit6080");
	graphPi0SysDataToFit6080 = CalculateGraphErrRatioToFit(graphPi0SysDataToFit6080, fitYieldPi0PbPb6080);
	TGraphAsymmErrors*  graphPi0StatDataToFit6080 = (TGraphAsymmErrors*) graphYieldPi0CombPbPb6080StatErrUnscaled->Clone("graphPi0StatDataToFit6080");
	graphPi0StatDataToFit6080 = CalculateGraphErrRatioToFit(graphPi0StatDataToFit6080, fitYieldPi0PbPb6080);
		
		
		
		

	padRatioEPOSToData0005->cd();
	padRatioEPOSToData0005->SetLogx();     

		TH2F * ratio2DEPOSToDataFirst = new TH2F("ratio2DEPOSToDataFirst","ratio2DEPOSToDataFirst",1000,0.35,15.,1000,0.1,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToDataFirst, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.13,0.13, 0.11,0.115, 1,0.53, 512, 505);
		ratio2DEPOSToDataFirst->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataFirst->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToDataFirst->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataFirst->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToDataFirst->DrawCopy(); 

		DrawGammaSetMarkerTGraph(gRatioEPOSToFit0005, markerStyleCommmonSpectrum0005,markerSizeSpectrum, colorEPOS0005, colorEPOS0005,2,widthLinesBoxes, kTRUE);
		gRatioEPOSToFit0005->Draw("3 same");
		
		DrawGammaSetMarkerTGraphAsym(graphPi0SysDataToFit0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005, widthLinesBoxes, kTRUE);//, colorComb0005-5);
		graphPi0SysDataToFit0005->Draw("E2same");
		
		DrawGammaSetMarkerTGraphAsym(graphPi0StatDataToFit0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005);
		graphPi0StatDataToFit0005->Draw("p,same,e1");
		
		
		TGraphErrors* gRatioEPOSToFit0005_Legend = (TGraphErrors*)gRatioEPOSToFit0005->Clone("gRatioEPOSToFit0005_Legend");
		TGraphAsymmErrors* graphPi0SysDataToFit0005_Legend = (TGraphAsymmErrors*)graphPi0SysDataToFit0005->Clone("graphPi0SysDataToFit0005_Legend");
		TGraphAsymmErrors* graphPi0StatDataToFit0005_Legend = (TGraphAsymmErrors*)graphPi0StatDataToFit0005->Clone("graphPi0StatDataToFit0005_Legend");
		DrawGammaSetMarkerTGraph(gRatioEPOSToFit0005_Legend, markerStyleCommmonSpectrum0005,markerSizeSpectrum, kGray+1, kGray+1,2);
		
		DrawGammaSetMarkerTGraphAsym(graphPi0SysDataToFit0005_Legend, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, kBlack , kBlack, widthLinesBoxes, kTRUE);      
		DrawGammaSetMarkerTGraphAsym(graphPi0StatDataToFit0005_Legend, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, kBlack , kBlack);
		
		TLegend* legendRatioTheory = new TLegend(0.18,0.8,0.7,0.9);
		legendRatioTheory->SetFillColor(0);
		legendRatioTheory->SetLineColor(0);
		legendRatioTheory->SetTextSize(0.1);
		legendRatioTheory->SetNColumns(2);
		legendRatioTheory->AddEntry(graphPi0SysDataToFit0005_Legend,"#pi^{0} ALICE","pf");
		legendRatioTheory->AddEntry(gRatioEPOSToFit0005_Legend, "EPOS", "f");  
		legendRatioTheory->Draw();

		DrawGammaLines(0., 30.,1., 1.,0.1);

		TLatex *labelRatioEPOSToFit0005 = new TLatex(0.18,0.07,collisionSystemCent0005.Data());
		SetStyleTLatex( labelRatioEPOSToFit0005, 0.1,4);
		labelRatioEPOSToFit0005->Draw();

		
	padRatioEPOSToData0510->cd();
	padRatioEPOSToData0510->SetLogx();     
		TH2F * ratio2DEPOSToDataMiddle = new TH2F("ratio2DEPOSToDataMiddle","ratio2DEPOSToDataMiddle",1000,0.35,15.,1000,0.1,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToDataMiddle, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.13,0.13, 0.11,0.115, 1,0.53, 512, 505);
		ratio2DEPOSToDataMiddle->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataMiddle->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToDataMiddle->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataMiddle->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToDataMiddle->DrawCopy(); 

		DrawGammaSetMarkerTGraph(gRatioEPOSToFit0510, markerStyleCommmonSpectrum0510,markerSizeSpectrum, colorEPOS0510, colorEPOS0510,2);
		gRatioEPOSToFit0510->Draw("3 same");

		DrawGammaSetMarkerTGraphAsym(graphPi0SysDataToFit0510, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510, widthLinesBoxes, kTRUE);//, colorComb0510-5);
		graphPi0SysDataToFit0510->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPi0StatDataToFit0510, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510);
		graphPi0StatDataToFit0510->Draw("p,same,e1");
			
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);
		
		TLatex *labelRatioEPOSToFit0510 = new TLatex(0.18,0.07,collisionSystemCent0510.Data());
		SetStyleTLatex( labelRatioEPOSToFit0510, 0.1,4);
		labelRatioEPOSToFit0510->Draw();

	padRatioEPOSToData1020->cd();
	padRatioEPOSToData1020->SetLogx();     
		
		ratio2DEPOSToDataMiddle->DrawCopy(); 
		
		DrawGammaSetMarkerTGraph(gRatioEPOSToFit1020, markerStyleCommmonSpectrum1020,markerSizeSpectrum, colorEPOS1020, colorEPOS1020,2);
		gRatioEPOSToFit1020->Draw("3 same");
		
		DrawGammaSetMarkerTGraphAsym(graphPi0SysDataToFit1020, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020, widthLinesBoxes, kTRUE);//, colorComb1020-5);
		graphPi0SysDataToFit1020->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPi0StatDataToFit1020, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020);
		graphPi0StatDataToFit1020->Draw("p,same,e1");
		
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);

		TLatex *labelRatioEPOSToFit1020 = new TLatex(0.18,0.07,collisionSystemCent1020.Data());
		SetStyleTLatex( labelRatioEPOSToFit1020, 0.1,4);
		labelRatioEPOSToFit1020->Draw();
		
	padRatioEPOSToData2040->cd();
	padRatioEPOSToData2040->SetLogx();

		ratio2DEPOSToDataMiddle->DrawCopy(); 
		
		DrawGammaSetMarkerTGraph(gRatioEPOSToFit2040, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorEPOS2040, colorEPOS2040,2);
		gRatioEPOSToFit2040->Draw("3 same");
		
		DrawGammaSetMarkerTGraphAsym(graphPi0SysDataToFit2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);//, colorComb2040-5);
		graphPi0SysDataToFit2040->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPi0StatDataToFit2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
		graphPi0StatDataToFit2040->Draw("p,same,e1");
			
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);

		TLatex *labelRatioEPOSToFit2040 = new TLatex(0.18,0.07,collisionSystemCent2040.Data());
		SetStyleTLatex( labelRatioEPOSToFit2040, 0.1,4);
		labelRatioEPOSToFit2040->Draw();

	padRatioEPOSToData4060->cd();
	padRatioEPOSToData4060->SetLogx();
		ratio2DEPOSToDataMiddle->DrawCopy(); 
		
		DrawGammaSetMarkerTGraph(gRatioEPOSToFit4060, markerStyleCommmonSpectrum4060,markerSizeSpectrum, colorEPOS4060, colorEPOS4060,2);
		gRatioEPOSToFit4060->Draw("3 same");
		
		DrawGammaSetMarkerTGraphAsym(graphPi0SysDataToFit4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060, widthLinesBoxes, kTRUE);//, colorComb4060-5);
		graphPi0SysDataToFit4060->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPi0StatDataToFit4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
		graphPi0StatDataToFit4060->Draw("p,same,e1");
		
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);
		
		TLatex *labelRatioEPOSToFit4060 = new TLatex(0.18,0.07,collisionSystemCent4060.Data());
		SetStyleTLatex( labelRatioEPOSToFit4060, 0.1,4);
		labelRatioEPOSToFit4060->Draw();

	padRatioEPOSToData6080->cd();
	padRatioEPOSToData6080->SetLogx();

		TH2F * ratio2DEPOSToData6080=  new TH2F("ratio2DEPOSToData6080","ratio2DEPOSToData6080",1000,0.35,15.,1000,0.1,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToData6080, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.09,0.095, 0.09,0.095, 0.95,0.63, 510, 505);
		ratio2DEPOSToData6080->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToData6080->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToData6080->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToData6080->GetXaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToData6080->GetXaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToData6080->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToData6080->DrawCopy(); 

		DrawGammaSetMarkerTGraph(gRatioEPOSToFit6080, markerStyleCommmonSpectrum6080,markerSizeSpectrum, colorEPOS6080, colorEPOS6080,2);
		gRatioEPOSToFit6080->Draw("3 same");
		
		DrawGammaSetMarkerTGraphAsym(graphPi0SysDataToFit6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080, widthLinesBoxes, kTRUE);//, colorComb6080-5);
		graphPi0SysDataToFit6080->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPi0StatDataToFit6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
		graphPi0StatDataToFit6080->Draw("p,same,e1");
		
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);

		TLatex *labelRatioEPOSToFit6080 = new TLatex(0.18,0.25,collisionSystemCent6080.Data());
		SetStyleTLatex( labelRatioEPOSToFit6080, 0.075,4);
		labelRatioEPOSToFit6080->Draw();

	canvasRatioEPOSToData->Update();
	canvasRatioEPOSToData->Print(Form("%s/RatioEPOSToData.%s",outputDir.Data(),suffix.Data()));


	Double_t arrayBoundariesX1_4[2];
	Double_t arrayBoundariesY1_4[5];
	Double_t relativeMarginsX[3];
	Double_t relativeMarginsY[3];
	ReturnCorrectValuesForCanvasScaling(1200,1570, 1, 3,0.13, 0.005, 0.005,0.07,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);
	cout << arrayBoundariesX1_4[0] << "\t" << arrayBoundariesX1_4[1] <<endl; 
	cout << arrayBoundariesY1_4[0] << "\t" << arrayBoundariesY1_4[1] << "\t" << arrayBoundariesY1_4[2] << "\t" << arrayBoundariesY1_4[3]<< "\t" << arrayBoundariesY1_4[4]<<endl; 
	cout << relativeMarginsX[0] << "\t" << relativeMarginsX[1] << "\t" << relativeMarginsX[2] << endl;
	cout << relativeMarginsY[0] << "\t" << relativeMarginsY[1] << "\t" << relativeMarginsY[2] << endl;

	TCanvas* canvasRatioEPOSToData3Parted = new TCanvas("canvasRatioEPOSToData3Parted","",200,10,1200,1570);  // gives the page size
	//    DrawGammaCanvasSettings( canvasRatioEPOSToData3Parted,  0.13, 0.02, 0.03, 0.06);

	TPad* padRatioEPOSToData0005_3Parted = new TPad("padRatioEPOSToData0005_3Parted", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
	DrawGammaPadSettings( padRatioEPOSToData0005_3Parted, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
	padRatioEPOSToData0005_3Parted->Draw();

	TPad* padRatioEPOSToData2040_3Parted = new TPad("padRatioEPOSToData2040_3Parted", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
	DrawGammaPadSettings( padRatioEPOSToData2040_3Parted, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[1]);
	padRatioEPOSToData2040_3Parted->Draw();

	TPad* padRatioEPOSToData6080_3Parted = new TPad("padRatioEPOSToData6080_3Parted", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[3], arrayBoundariesX1_4[1], arrayBoundariesY1_4[2],-1, -1, -2);
	DrawGammaPadSettings( padRatioEPOSToData6080_3Parted, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
	padRatioEPOSToData6080_3Parted->Draw();
			
	padRatioEPOSToData0005_3Parted->cd();
	padRatioEPOSToData0005_3Parted->SetLogx();     

		Int_t textSizeLabelsPixel = 50;
		Double_t margin = relativeMarginsX[0]*1700;
		Double_t textsizeLabels0005 = 0;
		Double_t textsizeFac0005 = 0;
		if (padRatioEPOSToData0005_3Parted->XtoPixel(padRatioEPOSToData0005_3Parted->GetX2()) < padRatioEPOSToData0005_3Parted->YtoPixel(padRatioEPOSToData0005_3Parted->GetY1())){
			textsizeLabels0005 = (Double_t)textSizeLabelsPixel/padRatioEPOSToData0005_3Parted->XtoPixel(padRatioEPOSToData0005_3Parted->GetX2()) ;
			textsizeFac0005 = (Double_t)1./padRatioEPOSToData0005_3Parted->XtoPixel(padRatioEPOSToData0005_3Parted->GetX2()) ;
		} else {
		textsizeLabels0005 = (Double_t)textSizeLabelsPixel/padRatioEPOSToData0005_3Parted->YtoPixel(padRatioEPOSToData0005_3Parted->GetY1());
		textsizeFac0005 = (Double_t)1./padRatioEPOSToData0005_3Parted->YtoPixel(padRatioEPOSToData0005_3Parted->GetY1());
		}
		cout << textsizeLabels0005 << endl;

		TH2F * ratio2DEPOSToDataFirst_3Parted = new TH2F("ratio2DEPOSToDataFirst_3Parted","ratio2DEPOSToDataFirst_3Parted",1000,0.35,15.,1000,-0.1,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToDataFirst_3Parted, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabels0005,textsizeLabels0005, 0.85*textsizeLabels0005,textsizeLabels0005, 0.8,0.25/(textsizeFac0005*margin), 512, 505);
		cout << "bla "<<  0.22/(textsizeFac0005*margin) << endl;
		ratio2DEPOSToDataFirst_3Parted->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataFirst_3Parted->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToDataFirst_3Parted->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataFirst_3Parted->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToDataFirst_3Parted->GetXaxis()->SetLabelFont(42);
		ratio2DEPOSToDataFirst_3Parted->GetYaxis()->SetLabelFont(42);

		ratio2DEPOSToDataFirst_3Parted->DrawCopy(); 

		DrawGammaNLOTGraph( gRatioEPOSToFit0005, widthCommonFit, 1, colorEPOS);
		gRatioEPOSToFit0005->Draw("l same");

		DrawGammaNLOTGraph( gRatioKopeliovichHydroToFit0005, widthCommonFit, styleKopeliovichHydro, colorKopeliovichHydro);
		gRatioKopeliovichHydroToFit0005->Draw("l same");

		DrawGammaNLOTGraph( gRatioKopeliovichELossToFit0005, widthCommonFit, styleKopeliovichELoss, colorKopeliovichELoss);
		gRatioKopeliovichELossToFit0005->Draw("l same");
		
		DrawGammaNLOTGraph( gRatioKopeliovichTotalToFit0005, widthCommonFit, styleKopeliovichComb, colorKopeliovichComb);
		gRatioKopeliovichTotalToFit0005->Draw("l same");

		DrawGammaSetMarkerTGraphAsym(graphPi0SysDataToFit0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005, widthLinesBoxes, kTRUE);//, colorComb0005-5);
		graphPi0SysDataToFit0005->Draw("E2same");
		
		DrawGammaSetMarkerTGraphAsym(graphPi0StatDataToFit0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005);
		graphPi0StatDataToFit0005->Draw("p,same,e1");
		
		TLegend* legendRatioTheory0005_3Parted = new TLegend(0.155,0.65,0.9,0.92);
		legendRatioTheory0005_3Parted->SetFillColor(0);
		legendRatioTheory0005_3Parted->SetLineColor(0);
		legendRatioTheory0005_3Parted->SetTextFont(42);
		legendRatioTheory0005_3Parted->SetTextSize(0.85*textsizeLabels0005);
		legendRatioTheory0005_3Parted->SetNColumns(3);
		legendRatioTheory0005_3Parted->SetMargin(0.24);
		legendRatioTheory0005_3Parted->AddEntry(graphPi0SysDataToFit0005,"#pi^{0} ALICE","pf");
		legendRatioTheory0005_3Parted->AddEntry(gRatioEPOSToFit0005, "EPOS", "l");  
		legendRatioTheory0005_3Parted->AddEntry(gRatioKopeliovichHydroToFit0005, "Nemchik hydro", "l");  
		legendRatioTheory0005_3Parted->AddEntry((TObject*)0, "","");
		legendRatioTheory0005_3Parted->AddEntry((TObject*)0, "","");  
		legendRatioTheory0005_3Parted->AddEntry(gRatioKopeliovichELossToFit0005, "Nemchik dipole abs.", "l");  
		legendRatioTheory0005_3Parted->AddEntry((TObject*)0, "","");
		legendRatioTheory0005_3Parted->AddEntry((TObject*)0, "","");  
		legendRatioTheory0005_3Parted->AddEntry(gRatioKopeliovichTotalToFit0005, "Nemchik sum", "l");  
		
		legendRatioTheory0005_3Parted->Draw();
		
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);

		TLatex *labelRatioEPOSToFit0005_3Parted = new TLatex(0.16,0.06,collisionSystemCent0005.Data());
		SetStyleTLatex( labelRatioEPOSToFit0005_3Parted, 0.85*textsizeLabels0005,4);
		labelRatioEPOSToFit0005_3Parted->Draw();
		
	padRatioEPOSToData2040_3Parted->cd();
	padRatioEPOSToData2040_3Parted->SetLogx();

		Double_t textsizeLabels2040 = 0;
		Double_t textsizeFac2040 = 0;
		if (padRatioEPOSToData2040_3Parted->XtoPixel(padRatioEPOSToData2040_3Parted->GetX2()) <padRatioEPOSToData2040_3Parted->YtoPixel(padRatioEPOSToData2040_3Parted->GetY1()) ){
			textsizeLabels2040 = (Double_t)textSizeLabelsPixel/padRatioEPOSToData2040_3Parted->XtoPixel(padRatioEPOSToData2040_3Parted->GetX2()) ;
			textsizeFac2040 = (Double_t)1./padRatioEPOSToData2040_3Parted->XtoPixel(padRatioEPOSToData2040_3Parted->GetX2()) ;
		} else {
		textsizeLabels2040 = (Double_t)textSizeLabelsPixel/padRatioEPOSToData2040_3Parted->YtoPixel(padRatioEPOSToData2040_3Parted->GetY1());
		textsizeFac2040 = (Double_t)1./padRatioEPOSToData2040_3Parted->YtoPixel(padRatioEPOSToData2040_3Parted->GetY1());
		}
		cout << textsizeLabels2040 << endl;

		TH2F * ratio2DEPOSToDataMiddle_3Parted = new TH2F("ratio2DEPOSToDataMiddle_3Parted","ratio2DEPOSToDataMiddle_3Parted",1000,0.35,15.,1000,-0.1,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToDataMiddle_3Parted, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabels2040,textsizeLabels2040, 0.85*textsizeLabels2040,textsizeLabels2040, 0.8,0.25/(textsizeFac2040*margin), 512, 505);
		ratio2DEPOSToDataMiddle_3Parted->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataMiddle_3Parted->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToDataMiddle_3Parted->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataMiddle_3Parted->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToDataMiddle_3Parted->GetXaxis()->SetLabelFont(42);
		ratio2DEPOSToDataMiddle_3Parted->GetYaxis()->SetLabelFont(42);

		ratio2DEPOSToDataMiddle_3Parted->DrawCopy(); 
		
		DrawGammaNLOTGraph( gRatioEPOSToFit2040, widthCommonFit, 1, colorEPOS);
		gRatioEPOSToFit2040->Draw("l same");

		DrawGammaNLOTGraph( gRatioKopeliovichHydroToFit2040, widthCommonFit, styleKopeliovichHydro, colorKopeliovichHydro);
		gRatioKopeliovichHydroToFit2040->Draw("l same");

		DrawGammaNLOTGraph( gRatioKopeliovichELossToFit2040, widthCommonFit, styleKopeliovichELoss, colorKopeliovichELoss);
		gRatioKopeliovichELossToFit2040->Draw("l same");
		
		DrawGammaNLOTGraph( gRatioKopeliovichTotalToFit2040, widthCommonFit, styleKopeliovichComb, colorKopeliovichComb);
		gRatioKopeliovichTotalToFit2040->Draw("l same");
		
		DrawGammaSetMarkerTGraphAsym(graphPi0SysDataToFit2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);//, colorComb2040-5);
		graphPi0SysDataToFit2040->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPi0StatDataToFit2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
		graphPi0StatDataToFit2040->Draw("p,same,e1");
			
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);

	/*     TLegend* legendRatioTheory2040_3Parted = new TLegend(0.16,0.7,0.7,0.8);
		legendRatioTheory2040_3Parted->SetFillColor(0);
		legendRatioTheory2040_3Parted->SetLineColor(0);
		legendRatioTheory2040_3Parted->SetTextSize(textsizeLabels2040);
		legendRatioTheory2040_3Parted->SetNColumns(2);
		legendRatioTheory2040_3Parted->AddEntry(graphPi0SysDataToFit2040,"#pi^{0} ALICE","pf");
		legendRatioTheory2040_3Parted->AddEntry(gRatioEPOSToFit2040, "EPOS", "f");  
		legendRatioTheory2040_3Parted->Draw();

	*/     
		TLatex *labelRatioEPOSToFit2040_3Parted = new TLatex(0.16,0.06,collisionSystemCent2040.Data());
		SetStyleTLatex( labelRatioEPOSToFit2040_3Parted, 0.85*textsizeLabels2040,4);
		labelRatioEPOSToFit2040_3Parted->Draw();

		
	padRatioEPOSToData6080_3Parted->cd();
	padRatioEPOSToData6080_3Parted->SetLogx();

		Double_t textsizeLabels6080 = 0;
		Double_t textsizeFac6080= 0;
		if (padRatioEPOSToData6080_3Parted->XtoPixel(padRatioEPOSToData6080_3Parted->GetX2()) <padRatioEPOSToData6080_3Parted->YtoPixel(padRatioEPOSToData6080_3Parted->GetY1()) ){
			textsizeLabels6080 = (Double_t)textSizeLabelsPixel/padRatioEPOSToData6080_3Parted->XtoPixel(padRatioEPOSToData6080_3Parted->GetX2()) ;
			textsizeFac6080 = (Double_t)1./padRatioEPOSToData6080_3Parted->XtoPixel(padRatioEPOSToData6080_3Parted->GetX2()) ;
		} else {
		textsizeLabels6080 = (Double_t)textSizeLabelsPixel/padRatioEPOSToData6080_3Parted->YtoPixel(padRatioEPOSToData6080_3Parted->GetY1());
		textsizeFac6080 = (Double_t)1./padRatioEPOSToData6080_3Parted->YtoPixel(padRatioEPOSToData6080_3Parted->GetY1());
		}
		cout << textsizeLabels6080 << endl;

		TH2F * ratio2DEPOSToDataBottom=  new TH2F("ratio2DEPOSToDataBottom","ratio2DEPOSToDataBottom",1000,0.35,15.,1000,-0.1,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToDataBottom, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabels6080,textsizeLabels6080, 0.85*textsizeLabels6080,textsizeLabels6080, 0.9,0.25/(textsizeFac6080*margin), 510, 505);
		
		ratio2DEPOSToDataBottom->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataBottom->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToDataBottom->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataBottom->GetXaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataBottom->GetXaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataBottom->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToDataBottom->GetXaxis()->SetLabelFont(42);
		ratio2DEPOSToDataBottom->GetYaxis()->SetLabelFont(42);
		ratio2DEPOSToDataBottom->DrawCopy(); 

		DrawGammaNLOTGraph( gRatioEPOSToFit6080, widthCommonFit, 1, colorEPOS);
		gRatioEPOSToFit6080->Draw("l same");
		
		DrawGammaNLOTGraph( gRatioKopeliovichHydroToFit6080, widthCommonFit, styleKopeliovichHydro, colorKopeliovichHydro);
		gRatioKopeliovichHydroToFit6080->Draw("l same");

		DrawGammaNLOTGraph( gRatioKopeliovichELossToFit6080, widthCommonFit, styleKopeliovichELoss, colorKopeliovichELoss);
		gRatioKopeliovichELossToFit6080->Draw("l same");
		
		DrawGammaNLOTGraph( gRatioKopeliovichTotalToFit6080, widthCommonFit, styleKopeliovichComb, colorKopeliovichComb);
		gRatioKopeliovichTotalToFit6080->Draw("l same");

		
		DrawGammaSetMarkerTGraphAsym(graphPi0SysDataToFit6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080, widthLinesBoxes, kTRUE);//, colorComb6080-5);
		graphPi0SysDataToFit6080->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphPi0StatDataToFit6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
		graphPi0StatDataToFit6080->Draw("p,same,e1");
		
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);

	//       TLegend* legendRatioTheory6080_3Parted = new TLegend(0.16,0.7,0.7,0.8);
	//       legendRatioTheory6080_3Parted->SetFillColor(0);
	//       legendRatioTheory6080_3Parted->SetLineColor(0);
	//       legendRatioTheory6080_3Parted->SetTextSize(textsizeLabels6080);
	//       legendRatioTheory6080_3Parted->SetNColumns(2);
	//       legendRatioTheory6080_3Parted->AddEntry(graphPi0SysDataToFit6080,"#pi^{0} ALICE","pf");
	//       legendRatioTheory6080_3Parted->AddEntry(gRatioEPOSToFit6080, "EPOS", "f");  
	//       legendRatioTheory6080_3Parted->Draw();

		TLatex *labelRatioEPOSToFit6080_3Parted = new TLatex(0.16,0.21,collisionSystemCent6080.Data());
		SetStyleTLatex( labelRatioEPOSToFit6080_3Parted, 0.85*textsizeLabels6080,4);
		labelRatioEPOSToFit6080_3Parted->Draw();

	canvasRatioEPOSToData3Parted->Update();
	canvasRatioEPOSToData3Parted->Print(Form("%s/RatioEPOSToData_3Parted.%s",outputDir.Data(),suffix.Data()));


	canvasRatioEPOSToData3Parted->cd();
	padRatioEPOSToData0005_3Parted->Draw();

	padRatioEPOSToData2040_3Parted->Draw();

	padRatioEPOSToData6080_3Parted->Draw();
			
	padRatioEPOSToData0005_3Parted->cd();
	padRatioEPOSToData0005_3Parted->SetLogx(0);     

		TH2F * ratio2DEPOSToDataFirst_3Parted_lin = new TH2F("ratio2DEPOSToDataFirst_3Parted_lin","ratio2DEPOSToDataFirst_3Parted_lin",1000,0.,12.5,1000,-0.1,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToDataFirst_3Parted_lin, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabels0005,textsizeLabels0005, 0.85*textsizeLabels0005,textsizeLabels0005, 0.8,0.25/(textsizeFac0005*margin), 512, 505);
		ratio2DEPOSToDataFirst_3Parted_lin->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataFirst_3Parted_lin->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToDataFirst_3Parted_lin->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataFirst_3Parted_lin->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToDataFirst_3Parted_lin->GetXaxis()->SetLabelFont(42);
		ratio2DEPOSToDataFirst_3Parted_lin->GetYaxis()->SetLabelFont(42);

		ratio2DEPOSToDataFirst_3Parted_lin->DrawCopy(); 

		gRatioEPOSToFit0005->Draw("l same");
		gRatioKopeliovichHydroToFit0005->Draw("l same");
		gRatioKopeliovichELossToFit0005->Draw("l same");
		gRatioKopeliovichTotalToFit0005->Draw("l same");
		graphPi0SysDataToFit0005->Draw("E2same");
		graphPi0StatDataToFit0005->Draw("p,same,e1");
		
	//       TLegend* legendRatioTheory0005_3Parted = new TLegend(0.155,0.65,0.9,0.92);
	//       legendRatioTheory0005_3Parted->SetFillColor(0);
	//       legendRatioTheory0005_3Parted->SetLineColor(0);
	//       legendRatioTheory0005_3Parted->SetTextFont(42);
	//       legendRatioTheory0005_3Parted->SetTextSize(0.85*textsizeLabels0005);
	//       legendRatioTheory0005_3Parted->SetNColumns(3);
	//       legendRatioTheory0005_3Parted->SetMargin(0.24);
	//       legendRatioTheory0005_3Parted->AddEntry(graphPi0SysDataToFit0005,"#pi^{0} ALICE","pf");
	//       legendRatioTheory0005_3Parted->AddEntry(gRatioEPOSToFit0005, "EPOS", "l");  
	//       legendRatioTheory0005_3Parted->AddEntry(gRatioKopeliovichHydroToFit0005, "Nemchik hydro", "l");  
	//       legendRatioTheory0005_3Parted->AddEntry((TObject*)0, "","");
	//       legendRatioTheory0005_3Parted->AddEntry((TObject*)0, "","");  
	//       legendRatioTheory0005_3Parted->AddEntry(gRatioKopeliovichELossToFit0005, "Nemchik dipole abs.", "l");  
	//       legendRatioTheory0005_3Parted->AddEntry((TObject*)0, "","");
	//       legendRatioTheory0005_3Parted->AddEntry((TObject*)0, "","");  
	//       legendRatioTheory0005_3Parted->AddEntry(gRatioKopeliovichTotalToFit0005, "Nemchik sum", "l");  
	//       
		legendRatioTheory0005_3Parted->Draw();
		
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);

		labelRatioEPOSToFit0005_3Parted->Draw();
		
	padRatioEPOSToData2040_3Parted->cd();
	padRatioEPOSToData2040_3Parted->SetLogx(0);

		TH2F * ratio2DEPOSToDataMiddle_3Parted_lin = new TH2F("ratio2DEPOSToDataMiddle_3Parted_lin","ratio2DEPOSToDataMiddle_3Parted_lin",1000,0.,12.5,1000,-0.1,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToDataMiddle_3Parted_lin, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabels2040,textsizeLabels2040, 0.85*textsizeLabels2040,textsizeLabels2040, 0.8,0.25/(textsizeFac2040*margin), 512, 505);
		ratio2DEPOSToDataMiddle_3Parted_lin->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataMiddle_3Parted_lin->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToDataMiddle_3Parted_lin->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataMiddle_3Parted_lin->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToDataMiddle_3Parted_lin->GetXaxis()->SetLabelFont(42);
		ratio2DEPOSToDataMiddle_3Parted_lin->GetYaxis()->SetLabelFont(42);
		ratio2DEPOSToDataMiddle_3Parted_lin->DrawCopy(); 
		
		gRatioEPOSToFit2040->Draw("l same");
		gRatioKopeliovichHydroToFit2040->Draw("l same");
		gRatioKopeliovichELossToFit2040->Draw("l same");
		gRatioKopeliovichTotalToFit2040->Draw("l same");
		graphPi0SysDataToFit2040->Draw("E2same");
		graphPi0StatDataToFit2040->Draw("p,same,e1");
			
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);
		labelRatioEPOSToFit2040_3Parted->Draw();

		
	padRatioEPOSToData6080_3Parted->cd();
	padRatioEPOSToData6080_3Parted->SetLogx(0);

		TH2F * ratio2DEPOSToDataBottom_lin=  new TH2F("ratio2DEPOSToDataBottom_lin","ratio2DEPOSToDataBottom_lin",1000,0.,12.5,1000,-0.1,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToDataBottom_lin, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabels6080,textsizeLabels6080, 0.85*textsizeLabels6080,textsizeLabels6080, 0.9,0.25/(textsizeFac6080*margin), 510, 505);
		
		ratio2DEPOSToDataBottom_lin->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataBottom_lin->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToDataBottom_lin->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataBottom_lin->GetXaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataBottom_lin->GetXaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataBottom_lin->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToDataBottom_lin->GetXaxis()->SetLabelFont(42);
		ratio2DEPOSToDataBottom_lin->GetYaxis()->SetLabelFont(42);
		ratio2DEPOSToDataBottom_lin->DrawCopy(); 

		gRatioEPOSToFit6080->Draw("l same");
		gRatioKopeliovichHydroToFit6080->Draw("l same");
		gRatioKopeliovichELossToFit6080->Draw("l same");
		gRatioKopeliovichTotalToFit6080->Draw("l same");
		graphPi0SysDataToFit6080->Draw("E2same");
		graphPi0StatDataToFit6080->Draw("p,same,e1");
		
		DrawGammaLines(0., 21.,1., 1.,0.1,kGray);
		labelRatioEPOSToFit6080_3Parted->Draw();

	canvasRatioEPOSToData3Parted->Update();
	canvasRatioEPOSToData3Parted->Print(Form("%s/RatioEPOSToData_3Parted_linX.%s",outputDir.Data(),suffix.Data()));

		cout << "here ............................." << endl;
		
	TCanvas* canvasRatioPP = new TCanvas("canvasRatioPP","",200,10,1200,700);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioPP,  0.105, 0.01, 0.01, 0.15);
	canvasRatioPP->cd();
	canvasRatioPP->SetLogx();

		Double_t textsizeLabelsPP = 0;
		Double_t textsizeFacPP= 0;
		if (canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) <canvasRatioPP->YtoPixel(canvasRatioPP->GetY1()) ){
			textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
			textsizeFacPP = (Double_t)1./canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
		} else {
		textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
		textsizeFacPP = (Double_t)1./canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
		}
		cout << textsizeLabelsPP << endl;

		TH2F * ratio2DEPOSToDataPP=  new TH2F("ratio2DEPOSToDataPP","ratio2DEPOSToDataPP",1000,0.35,15.,1000,0.1,3.6);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToDataPP, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsPP,textsizeLabelsPP, 0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9,0.17/(textsizeFacPP*margin), 510, 505);
		ratio2DEPOSToDataPP->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataPP->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToDataPP->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataPP->GetXaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataPP->GetXaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataPP->GetYaxis()->SetLabelOffset(0.01);
		ratio2DEPOSToDataPP->GetXaxis()->SetLabelFont(42);
		ratio2DEPOSToDataPP->GetYaxis()->SetLabelFont(42);
	//       ratio2DEPOSToDataPP->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToDataPP->DrawCopy(); 
		
		DrawGammaSetMarker(histoRatioPythia8ToFit2760GeV,24,markerSizeChargedHadronSpectrum, kRed+2 , kRed+2);  
		histoRatioPythia8ToFit2760GeV->SetLineWidth(widthCommonFit);
		histoRatioPythia8ToFit2760GeV->GetXaxis()->SetRangeUser(0.5,14);
		histoRatioPythia8ToFit2760GeV->Draw("same,hist,c");
		DrawGammaNLOTGraph( graphRatioCombNLOPi02760GeVMuHalf, widthCommonFit, styleLineNLOMuHalf, kGray+2);
		graphRatioCombNLOPi02760GeVMuHalf->Draw("same,c");
		DrawGammaNLOTGraph( graphRatioCombNLOPi02760GeVMuOne, widthCommonFit, styleLineNLOMuOne, kGray+1);
		graphRatioCombNLOPi02760GeVMuOne->Draw("same,c");
		DrawGammaNLOTGraph( graphRatioCombNLOPi02760GeVMuTwo, widthCommonFit, styleLineNLOMuTwo, kGray+1);
		graphRatioCombNLOPi02760GeVMuTwo->Draw("same,c");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit2760GeVStat, markerStyleCommmonSpectrumpp,markerSizeSpectrum, kBlack, kBlack, widthCommonSpectrumBoxes, kFALSE);
		graphRatioCombCombFit2760GeVStat->SetLineWidth(widthCommonErrors);
		graphRatioCombCombFit2760GeVStat->RemovePoint(graphRatioCombCombFit2760GeVStat->GetN()-1);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit2760GeVSys,markerStyleCommmonSpectrumpp,markerSizeSpectrum, kBlack, kBlack, widthCommonSpectrumBoxes, kTRUE, 0);
		graphRatioCombCombFit2760GeVSys->SetLineWidth(0);
		graphRatioCombCombFit2760GeVSys->RemovePoint(graphRatioCombCombFit2760GeVSys->GetN()-1);
		graphRatioCombCombFit2760GeVSys->Draw("2,same");
		graphRatioCombCombFit2760GeVStat->Draw("p,same");

		
		DrawGammaLines(0., 15.,1., 1.,0.1,kGray);

		TLegend* legendRatioTheorypp_3Parted = new TLegend(0.13,0.55,0.4,0.87);
		legendRatioTheorypp_3Parted->SetFillColor(0);
		legendRatioTheorypp_3Parted->SetLineColor(0);
		legendRatioTheorypp_3Parted->SetTextFont(42);
		legendRatioTheorypp_3Parted->SetTextSize(0.85*textsizeLabelsPP);
	//       legendRatioTheorypp_3Parted->SetNColumns(2);
		legendRatioTheorypp_3Parted->AddEntry(graphRatioCombCombFit2760GeVSys,"#pi^{0} ALICE","pf");
		legendRatioTheorypp_3Parted->AddEntry(graphRatioCombNLOPi02760GeVMuHalf, "NLO  #mu = 0.5 #it{p}_{T}", "l");  
		legendRatioTheorypp_3Parted->AddEntry(graphRatioCombNLOPi02760GeVMuOne,  "NLO  #mu = #it{p}_{T}", "l");  
		legendRatioTheorypp_3Parted->AddEntry(graphRatioCombNLOPi02760GeVMuTwo,  "NLO  #mu = 2 #it{p}_{T}", "l");  
		legendRatioTheorypp_3Parted->AddEntry(histoRatioPythia8ToFit2760GeV,  "Pythia 8, Tune 4C", "l");  
		legendRatioTheorypp_3Parted->Draw();
		
		TLatex *labelRatioEPOSToFitpp2760GeV = new TLatex(0.13,0.9,collisionSystemPP.Data());
		SetStyleTLatex( labelRatioEPOSToFitpp2760GeV, 0.85*textsizeLabelsPP,4);
		labelRatioEPOSToFitpp2760GeV->Draw();

		
	canvasRatioPP->Update();
	canvasRatioPP->Print(Form("%s/RatioTheoryToData_PP.%s",outputDir.Data(),suffix.Data()));

	canvasRatioPP->cd();
	canvasRatioPP->SetLogx(0);

		TH2F * ratio2DEPOSToDataPP_linX=  new TH2F("ratio2DEPOSToDataPP_linX","ratio2DEPOSToDataPP_linX",1000,0.,12.5,1000,0.1,3.6);
		SetStyleHistoTH2ForGraphs(ratio2DEPOSToDataPP_linX, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsPP,textsizeLabelsPP, 0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9,0.17/(textsizeFacPP*margin), 510, 505);
		ratio2DEPOSToDataPP_linX->GetYaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataPP_linX->GetYaxis()->SetNdivisions(505);
		ratio2DEPOSToDataPP_linX->GetYaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataPP_linX->GetXaxis()->SetMoreLogLabels(kTRUE);
		ratio2DEPOSToDataPP_linX->GetXaxis()->SetNoExponent(kTRUE);
		ratio2DEPOSToDataPP_linX->GetYaxis()->SetLabelOffset(0.01);
		ratio2DEPOSToDataPP_linX->GetXaxis()->SetLabelFont(42);
		ratio2DEPOSToDataPP_linX->GetYaxis()->SetLabelFont(42);
	//       ratio2DEPOSToDataPP_linX->GetXaxis()->SetTickLength(0.07);
		ratio2DEPOSToDataPP_linX->DrawCopy(); 
		
		histoRatioPythia8ToFit2760GeV->Draw("same,hist,c");
		graphRatioCombNLOPi02760GeVMuHalf->Draw("same,c");
		graphRatioCombNLOPi02760GeVMuOne->Draw("same,c");
		graphRatioCombNLOPi02760GeVMuTwo->Draw("same,c");

		graphRatioCombCombFit2760GeVSys->Draw("2,same");
		graphRatioCombCombFit2760GeVStat->Draw("p,same");

		
		DrawGammaLines(0., 15.,1., 1.,0.1,kGray);
		legendRatioTheorypp_3Parted->Draw();      
		labelRatioEPOSToFitpp2760GeV->Draw();

		
	canvasRatioPP->Update();
	canvasRatioPP->Print(Form("%s/RatioTheoryToData_PP_linX.%s",outputDir.Data(),suffix.Data()));





	cout << "here ............................." << endl;


		
	TCanvas * canvas6PartRAA = new TCanvas("canvas6PartRAA","",10,10,2400,1300);  // gives the page size		
	canvas6PartRAA->cd();
	DrawGammaCanvasSettings( canvas6PartRAA, 0.13, 0.0, 0.02, 0.09);

	TPad* pad6PartPAA1 = new TPad("pad6PartPAA1", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartPAA1, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartPAA1->Draw();
	TPad* pad6PartPAA2 = new TPad("pad6PartPAA2", "",arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartPAA2, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartPAA2->Draw();

	TPad* pad6PartPAA3 = new TPad("pad6PartPAA3", "",arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartPAA3, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartPAA3->Draw();
	TPad* pad6PartPAA4 = new TPad("pad6PartPAA4", "",arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartPAA4, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartPAA4->Draw();

	TPad* pad6PartPAA5 = new TPad("pad6PartPAA5", "",arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartPAA5, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartPAA5->Draw();
	TPad* pad6PartPAA6 = new TPad("pad6PartPAA6", "", arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartPAA6, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartPAA6->Draw();

	TH2F * histo2DRAAAll3Up;
	histo2DRAAAll3Up = new TH2F("histo2DRAAAll3Up","histo2DRAAAll3Up",1000,0.,20.,1000,-0.05,10.	);
	histo2DRAAAll3Up->GetYaxis()->SetRangeUser(-0.05,1.4);
	histo2DRAAAll3Up->GetXaxis()->SetRangeUser(0.,15.);
	SetStyleHistoTH2ForGraphs(histo2DRAAAll3Up, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRatioUp,textsizeLabelsRatioUp,  0.85*textsizeLabelsRatioUp,textsizeLabelsRatioUp, 1,0.5/(textsizeFacRatioUp*marginRatio), 512, 505); //#frac{#frac{1}{#it{N}_{evt}^{AA}}#frac{dN^{AA}}{d#it{p}_{T} d#it{y}}}{<#it{N}_{coll}> #frac{1}{#it{N}_{evt}^{pp}}#frac{dN^{pp}}{d#it{p}_{T} d#it{y}}}
	histo2DRAAAll3Up->GetXaxis()->SetLabelFont(42);
	histo2DRAAAll3Up->GetYaxis()->SetLabelFont(42);

	TH2F * histo2DRAAAll3Down;
	histo2DRAAAll3Down = new TH2F("histo2DRAAAll3Down","histo2DRAAAll3Down",1000,0.,20.,1000,-0.05,10.	);
	histo2DRAAAll3Down->GetYaxis()->SetRangeUser(-0.05,1.4);
	histo2DRAAAll3Down->GetXaxis()->SetRangeUser(0.,15.);
	SetStyleHistoTH2ForGraphs(histo2DRAAAll3Down, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.85*textsizeLabelsRatioDown,textsizeLabelsRatioDown,  0.85*textsizeLabelsRatioDown,textsizeLabelsRatioDown, 1,0.5/(textsizeFacRatioDown*marginRatio), 512, 505); //#frac{#frac{1}{#it{N}_{evt}^{AA}}#frac{dN^{AA}}{d#it{p}_{T} d#it{y}}}{<#it{N}_{coll}> #frac{1}{#it{N}_{evt}^{pp}}#frac{dN^{pp}}{d#it{p}_{T} d#it{y}}}
	histo2DRAAAll3Down->GetXaxis()->SetLabelFont(42);
	histo2DRAAAll3Down->GetYaxis()->SetLabelFont(42);


	TH2F* histo2DRAAAll2Up;
	histo2DRAAAll2Up = new TH2F("histo2DRAAAll2Up","histo2DRAAAll2Up",1000,-0.25,20.,1000,-0.05,10.);
	histo2DRAAAll2Up->GetYaxis()->SetRangeUser(-0.05,1.4);
	histo2DRAAAll2Up->GetXaxis()->SetRangeUser(-0.18,15.);
	SetStyleHistoTH2ForGraphs(histo2DRAAAll2Up, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRatioUp,textsizeLabelsRatioUp,  0.85*textsizeLabelsRatioUp,textsizeLabelsRatioUp, 1,0.5/(textsizeFacRatioUp*marginRatio), 512, 505); //#frac{#frac{1}{#it{N}_{evt}^{AA}}#frac{dN^{AA}}{d#it{p}_{T} d#it{y}}}{<#it{N}_{coll}> #frac{1}{#it{N}_{evt}^{pp}}#frac{dN^{pp}}{d#it{p}_{T} d#it{y}}}
	histo2DRAAAll2Up->GetXaxis()->SetLabelFont(42);
	histo2DRAAAll2Up->GetYaxis()->SetLabelFont(42);

	TH2F* histo2DRAAAll2Down;
	histo2DRAAAll2Down = new TH2F("histo2DRAAAll2Down","histo2DRAAAll2Down",1000,-0.25,20.,1000,-0.05,10.);
	histo2DRAAAll2Down->GetYaxis()->SetRangeUser(-0.05,1.4);
	histo2DRAAAll2Down->GetXaxis()->SetRangeUser(-0.18,15.);
	SetStyleHistoTH2ForGraphs(histo2DRAAAll2Down, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRatioDown,textsizeLabelsRatioDown,  0.85*textsizeLabelsRatioDown,textsizeLabelsRatioDown, 1,0.5/(textsizeFacRatioDown*marginRatio), 512, 505); //#frac{#frac{1}{#it{N}_{evt}^{AA}}#frac{dN^{AA}}{d#it{p}_{T} d#it{y}}}{<#it{N}_{coll}> #frac{1}{#it{N}_{evt}^{pp}}#frac{dN^{pp}}{d#it{p}_{T} d#it{y}}}
	histo2DRAAAll2Down->GetXaxis()->SetLabelFont(42);
	histo2DRAAAll2Down->GetYaxis()->SetLabelFont(42);


	pad6PartPAA1->cd();
	// 	pad6PartPAA1->SetLogy();
	histo2DRAAAll3Up->DrawCopy();
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005, widthLinesBoxes, kTRUE, kGray+1);
	graphRAASysCombInd0005->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005);
	graphRAACombInd0005->Draw("p,same,e1");
	TBox* boxErrorNorm0005 = CreateBoxConv(colorComb0005Box, 14., 1.-normErr0005 , 14.3, 1.+normErr0005);
	boxErrorNorm0005->Draw();
	DrawGammaLines(0., 15. , 1, 1 ,1,kGray);

	TLatex *labelRAAPi0PbPb0005 = new TLatex(0.18,0.88,collisionSystemCent0005.Data());
	SetStyleTLatex( labelRAAPi0PbPb0005, 0.85*textsizeLabelsRatioUp,4);
	labelRAAPi0PbPb0005->Draw();
	// 	TLatex *labelRAAPi0LabelPbPb0005 = new TLatex(0.15,0.83,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelRAAPi0LabelPbPb0005, 0.04,4);
	// 	labelRAAPi0LabelPbPb0005->Draw();


	histo2DRAAAll3Up->Draw("axis,same");
	pad6PartPAA1->Update();
	pad6PartPAA2->cd();
	// 	pad6PartPAA2->SetLogy();
	histo2DRAAAll3Down->DrawCopy();
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE, kGray+1);//, colorComb2040-5);
	graphRAASysCombInd2040->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
	graphRAACombInd2040->Draw("p,same,e1");	
	TBox* boxErrorNorm2040 = CreateBoxConv(colorComb2040Box, 14., 1.-normErr2040 , 14.3, 1.+normErr2040);
	boxErrorNorm2040->Draw();
	DrawGammaLines(0., 15. , 1, 1 ,1,kGray);
	TLatex *labelRAAPi0PbPb2040 = new TLatex(0.18,0.92,collisionSystemCent2040.Data());
	SetStyleTLatex( labelRAAPi0PbPb2040, 0.85*textsizeLabelsRatioDown,4);
	labelRAAPi0PbPb2040->Draw();
	// 	TLatex *labelRAAPi0LabelPbPb2040 = new TLatex(0.03,0.83,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelRAAPi0LabelPbPb2040, 0.04,4);
	// 	labelRAAPi0LabelPbPb2040->Draw();

	histo2DRAAAll3Down->Draw("axis,same");
	pad6PartPAA2->Update();

	pad6PartPAA3->cd();
	// 	pad6PartPAA3->SetLogy();
	histo2DRAAAll2Up->DrawCopy();
	graphRAASysCombInd0510->RemovePoint(19);
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd0510, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510, widthLinesBoxes, kTRUE, kGray+1);
	graphRAASysCombInd0510->Draw("E2same");
	graphRAACombInd0510->RemovePoint(19);
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd0510, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510);
	graphRAACombInd0510->Draw("p,same,e1");
	TBox* boxErrorNorm0510 = CreateBoxConv(colorComb0510Box, 14., 1.-normErr0510 , 14.3, 1.+normErr0510);
	boxErrorNorm0510->Draw();
	DrawGammaLines(0., 15. , 1, 1 ,1,kGray);

	TLatex *labelRAAPi0PbPb0510 = new TLatex(0.05,0.88,collisionSystemCent0510.Data());
	SetStyleTLatex( labelRAAPi0PbPb0510, 0.85*textsizeLabelsRatioUp,4);
	labelRAAPi0PbPb0510->Draw();
	// 	TLatex *labelRAAPi0LabelPbPb0510 = new TLatex(0.15,0.83,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelRAAPi0LabelPbPb0510, 0.04,4);
	// 	labelRAAPi0LabelPbPb0510->Draw();

	histo2DRAAAll2Up->Draw("axis,same");
	pad6PartPAA3->Update();
	pad6PartPAA4->cd();
	// 	pad6PartPAA4->SetLogy();
	histo2DRAAAll2Down->DrawCopy();
		DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060, widthLinesBoxes, kTRUE, kGray+1);//, colorComb4060-5);
	graphRAASysCombInd4060->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
	graphRAACombInd4060->Draw("p,same,e1");
	TBox* boxErrorNorm4060 = CreateBoxConv(colorComb4060Box, 14., 1.-normErr4060 , 14.3, 1.+normErr4060);
	boxErrorNorm4060->Draw();
	DrawGammaLines(0., 15. , 1, 1 ,1,kGray);
	TLatex *labelRAAPi0PbPb04060 = new TLatex(0.05,0.92,collisionSystemCent4060.Data());
	SetStyleTLatex( labelRAAPi0PbPb04060, 0.85*textsizeLabelsRatioDown,4);
	labelRAAPi0PbPb04060->Draw();
	// 	TLatex *labelRAAPi0LabelPbPb4060 = new TLatex(0.15,0.87,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelRAAPi0LabelPbPb4060, 0.038,4);
	// 	labelRAAPi0LabelPbPb4060->Draw();
	histo2DRAAAll2Down->Draw("axis,same");
	pad6PartPAA4->Update();

	pad6PartPAA5->cd();
	// 	pad6PartPAA5->SetLogy();

	histo2DRAAAll2Up->DrawCopy();
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd1020, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020, widthLinesBoxes, kTRUE, kGray+1);
	graphRAASysCombInd1020->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd1020, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020);
	graphRAACombInd1020->Draw("p,same,e1");
	TBox* boxErrorNorm1020 = CreateBoxConv(colorComb1020Box, 14., 1.-normErr1020 , 14.3, 1.+normErr1020);
	boxErrorNorm1020->Draw();
	DrawGammaLines(0., 15. , 1, 1 ,1,kGray);

	TLatex *labelRAAPi0PbPb1020 = new TLatex(0.05,0.88,collisionSystemCent1020.Data());
	SetStyleTLatex( labelRAAPi0PbPb1020, 0.85*textsizeLabelsRatioUp,4);
	labelRAAPi0PbPb1020->Draw();
	// 	TLatex *labelRAAPi0LabelPbPb1020 = new TLatex(0.15,0.83,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelRAAPi0LabelPbPb1020, 0.04,4);
	// 	labelRAAPi0LabelPbPb1020->Draw();
	histo2DRAAAll2Up->Draw("axis,same");

	pad6PartPAA5->Update();

		pad6PartPAA6->cd();
	// 	pad6PartPAA6->SetLogy();
	histo2DRAAAll2Down->DrawCopy();
	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080, widthLinesBoxes, kTRUE, kGray+1);//, colorComb6080-5);
	graphRAASysCombInd6080->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
	graphRAACombInd6080->Draw("p,same,e1");
	TBox* boxErrorNorm6080 = CreateBoxConv(colorComb6080Box, 14., 1.-normErr6080 , 14.3, 1.+normErr6080);
	boxErrorNorm6080->Draw();
	DrawGammaLines(0., 15. , 1, 1 ,1,kGray);
	TLatex *labelRAAPi0PbPb6080 = new TLatex(0.05,0.92,collisionSystemCent6080.Data());
	SetStyleTLatex( labelRAAPi0PbPb6080, 0.85*textsizeLabelsRatioDown,4);
	labelRAAPi0PbPb6080->Draw();
	// 	TLatex *labelRAAPi0LabelPbPb6080 = new TLatex(0.03,0.87,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	// 	SetStyleTLatex( labelRAAPi0LabelPbPb6080, 0.038,4);
	// 	labelRAAPi0LabelPbPb6080->Draw();
	histo2DRAAAll2Down->Draw("axis,same");
	pad6PartPAA6->Update();

	canvas6PartRAA->Update();	
	canvas6PartRAA->SaveAs(Form("%s/RAA_All6Parted_Paper.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartPAA1;	
	delete pad6PartPAA2;	
	delete pad6PartPAA3;	
	delete pad6PartPAA4;	
	delete pad6PartPAA5;	
	delete pad6PartPAA6;	
	delete canvas6PartRAA;	
	cout << "here ............................." << endl;

		//**************************************************************************************************************
	//*********************************** Comparison to Theory *****************************************************
	//**************************************************************************************************************

		
	TCanvas * canvas6PartRAATheory = new TCanvas("canvas6PartRAATheory","",10,10,2400,1300);  // gives the page size		
	canvas6PartRAATheory->cd();
	DrawGammaCanvasSettings( canvas6PartRAATheory, 0.13, 0.0, 0.02, 0.09);

	TPad* pad6PartRAATheory1 = new TPad("pad6PartRAATheory1", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAATheory1, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRAATheory1->Draw();
	TPad* pad6PartRAATheory2 = new TPad("pad6PartRAATheory2", "",arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAATheory2, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRAATheory2->Draw();

	TPad* pad6PartRAATheory3 = new TPad("pad6PartRAATheory3", "",arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAATheory3, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRAATheory3->Draw();
	TPad* pad6PartRAATheory4 = new TPad("pad6PartRAATheory4", "",arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAATheory4, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRAATheory4->Draw();

	TPad* pad6PartRAATheory5 = new TPad("pad6PartRAATheory5", "",arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAATheory5, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRAATheory5->Draw();
	TPad* pad6PartRAATheory6 = new TPad("pad6PartRAATheory6", "", arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAATheory6, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRAATheory6->Draw();


	pad6PartRAATheory1->cd();
	histo2DRAAAll3Up->GetXaxis()->SetRangeUser(0.,19.5);
	histo2DRAAAll3Up->DrawCopy();

		graphRAASysCombInd0005->Draw("E2same");
		graphRAACombInd0005->Draw("p,same,e1");
		
	//       DrawGammaSetMarkerTGraphErr(gEPOS_RAA_0005, markerStyleCommmonSpectrum0005,markerSizeSpectrum, colorEPOS0005, colorEPOS0005,2);
	//       gEPOS_RAA_0005->SetFillStyle(fillStyleEPOS);
	//       gEPOS_RAA_0005->SetFillColor(colorEPOS0005);
	//       gEPOS_RAA_0005->Draw("3 same");
	//       DrawGammaSetMarkerTGraphErr(Xiao_Raa_0005, markerStyleCommmonSpectrum0005,markerSizeSpectrum, colorXiao0005 ,colorXiao0005,2);
	//       Xiao_Raa_0005->SetFillStyle(fillStyleXiao); 
	//       Xiao_Raa_0005->SetFillColor(colorXiao0005);
	//       Xiao_Raa_0005->Draw("3 same");
		DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_0005, markerStyleCommmonSpectrum0005,markerSizeSpectrum, colorVitevBas0005 ,colorVitevBas0005,2);
		Vitev_Bas_Raa_0005->SetFillStyle(fillStyleVitev);
		Vitev_Bas_Raa_0005->SetFillColor(colorVitevBas0005);
		Vitev_Bas_Raa_0005->Draw("3 same");
		DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_0005, markerStyleCommmonSpectrum0005,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,0.8);
		gWHDG_Raa_0005->SetFillStyle(fillStyleWHDG);
		gWHDG_Raa_0005->SetFillColor(colorWHDG0005);
		gWHDG_Raa_0005->Draw("3 same");
		
		graphRAACombInd0005->Draw("p,same,e1");
	//       gEPOS_RAA_0005->Draw("3 same");
		TLatex *labelRAAPi0TheoryPbPb0005 = new TLatex(0.18,0.9,collisionSystemCent0005.Data());
		SetStyleTLatex( labelRAAPi0TheoryPbPb0005, 0.85*textsizeLabelsRatioUp,4);
		labelRAAPi0TheoryPbPb0005->Draw();
		
		
		TLegend* legendRAAPi0TheoryPbPb0005 = new TLegend(0.17,0.79,0.9,0.86);
		legendRAAPi0TheoryPbPb0005->SetFillColor(0);
		legendRAAPi0TheoryPbPb0005->SetLineColor(0);
		legendRAAPi0TheoryPbPb0005->SetNColumns(2);
		legendRAAPi0TheoryPbPb0005->SetMargin(0.23);
		legendRAAPi0TheoryPbPb0005->SetTextSize(0.85*textsizeLabelsRatioUp);
		legendRAAPi0TheoryPbPb0005->SetTextFont(42);
		legendRAAPi0TheoryPbPb0005->AddEntry(graphRAASysCombInd0005,"#pi^{0} ALICE","pf");
		legendRAAPi0TheoryPbPb0005->AddEntry((TObject*)0, "", "");  
		legendRAAPi0TheoryPbPb0005->Draw();

		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		TBox* boxErrorNorm0005_Theory = CreateBoxConv(colorComb0005Box, 18., 1.-normErr0005 , 18.5, 1.+normErr0005);
		boxErrorNorm0005_Theory->Draw();


	//       TGraphErrors* Xiao_Raa_0020_Legend = (TGraphErrors*)Xiao_Raa_0020->Clone();
		TGraphErrors* Vitev_Bas_Raa_0020_Legend = (TGraphErrors*)Vitev_Bas_Raa_0020->Clone();
		TGraphAsymmErrors* gWHDG_Raa_0020_Legend = (TGraphAsymmErrors*)gWHDG_Raa_0020->Clone();
		
	//       DrawGammaSetMarkerTGraphErr(gEPOS_RAA_0020, markerStyleCommmonSpectrum0005,markerSizeSpectrum, kBlack, kBlack);
	//       gEPOS_RAA_0020->SetFillStyle(fillStyleEPOS);
	//       gEPOS_RAA_0020->SetFillColor(kBlack);
	//       DrawGammaSetMarkerTGraphErr(Xiao_Raa_0020_Legend, markerStyleCommmonSpectrum0005,markerSizeSpectrum, kGray+2 ,kGray+2);
	//       Xiao_Raa_0020_Legend->SetFillStyle(fillStyleXiao); 
	//       Xiao_Raa_0020_Legend->SetFillColor(kGray+2);
		DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_0020_Legend, markerStyleCommmonSpectrum0005,markerSizeSpectrum, kBlack ,kBlack);
		Vitev_Bas_Raa_0020_Legend->SetFillStyle(fillStyleVitev);
		Vitev_Bas_Raa_0020_Legend->SetFillColor(kBlack);
		DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_0020_Legend, markerStyleCommmonSpectrum0005,markerSizeSpectrum, kBlack, kBlack,0.8);
		gWHDG_Raa_0020_Legend->SetFillStyle(fillStyleWHDG);
		gWHDG_Raa_0020_Legend->SetFillColor(kBlack);
		
		TLegend* legendRAAPi0TheoryPbPb = new TLegend(0.5,0.79,0.98,0.86);
		legendRAAPi0TheoryPbPb->SetFillColor(0);
		legendRAAPi0TheoryPbPb->SetLineColor(0);
		legendRAAPi0TheoryPbPb->SetMargin(0.35);
		legendRAAPi0TheoryPbPb->SetNColumns(2);
		legendRAAPi0TheoryPbPb->SetTextSize(0.85*textsizeLabelsRatioUp);
		legendRAAPi0TheoryPbPb->SetTextFont(42);
	//       legendRAAPi0TheoryPbPb->AddEntry(Xiao_Raa_0020_Legend,"Chen (HT)","f");
		legendRAAPi0TheoryPbPb->AddEntry(Vitev_Bas_Raa_0020_Legend,"GLV","f");
		legendRAAPi0TheoryPbPb->AddEntry(gWHDG_Raa_0020_Legend,"WHDG","f");
	//       legendRAAPi0TheoryPbPb->AddEntry(gEPOS_RAA_0020,"EPOS","f");
		legendRAAPi0TheoryPbPb->Draw();

		
		histo2DRAAAll3Up->Draw("axis,same");

	pad6PartRAATheory1->Update();

	pad6PartRAATheory2->cd();
		histo2DRAAAll3Down->GetXaxis()->SetRangeUser(0.,19.5);
		histo2DRAAAll3Down->DrawCopy();
		
		graphRAASysCombInd2040->Draw("E2same");
		graphRAACombInd2040->Draw("p,same,e1");
			
	//       DrawGammaSetMarkerTGraphErr(gEPOS_RAA_2040, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorEPOS2040, colorEPOS2040,2);
	//       gEPOS_RAA_2040->SetFillStyle(fillStyleEPOS);
	//       gEPOS_RAA_2040->SetFillColor(colorEPOS2040);
	//       gEPOS_RAA_2040->Draw("3 same");

	//       DrawGammaSetMarkerTGraphErr(Xiao_Raa_2040, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorXiao2040 ,colorXiao2040,2);
	//       Xiao_Raa_2040->SetFillStyle(fillStyleXiao); 
	//       Xiao_Raa_2040->SetFillColor(colorXiao2040);
	//       Xiao_Raa_2040->Draw("3 same");   
		DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_2040, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorVitevBas2040 ,colorVitevBas2040,2);
		Vitev_Bas_Raa_2040->SetFillStyle(fillStyleVitev);
		Vitev_Bas_Raa_2040->SetFillColor(colorVitevBas2040);
		Vitev_Bas_Raa_2040->Draw("3 same");
		DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_2040, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG2040, colorWHDG2040,0.8);
		gWHDG_Raa_2040->SetFillStyle(fillStyleWHDG);
		gWHDG_Raa_2040->SetFillColor(colorWHDG2040);
		gWHDG_Raa_2040->Draw("3 same");
		graphRAACombInd2040->Draw("p,same,e1");
	//       gEPOS_RAA_2040->Draw("3 same");
		
		TLatex *labelRAAPi0TheoryPbPb2040 = new TLatex(0.18,0.93,collisionSystemCent2040.Data());
		SetStyleTLatex( labelRAAPi0TheoryPbPb2040, 0.85*textsizeLabelsRatioDown,4);
		labelRAAPi0TheoryPbPb2040->Draw();

		TLegend* legendRAAPi0TheoryPbPb2040 = new TLegend(0.17,0.85,0.9,0.91);
		legendRAAPi0TheoryPbPb2040->SetFillColor(0);
		legendRAAPi0TheoryPbPb2040->SetLineColor(0);
		legendRAAPi0TheoryPbPb2040->SetNColumns(2);
		legendRAAPi0TheoryPbPb2040->SetMargin(0.23);
		legendRAAPi0TheoryPbPb2040->SetTextSize(0.85*textsizeLabelsRatioDown);
		legendRAAPi0TheoryPbPb2040->SetTextFont(42);
		legendRAAPi0TheoryPbPb2040->AddEntry(graphRAASysCombInd2040,"#pi^{0} ALICE","pf");
		legendRAAPi0TheoryPbPb2040->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		TBox* boxErrorNorm2040_Theory = CreateBoxConv(colorComb2040Box, 18., 1.-normErr2040 , 18.5, 1.+normErr2040);
		boxErrorNorm2040_Theory->Draw();

		histo2DRAAAll3Down->Draw("axis,same");
		

	pad6PartRAATheory2->Update();

	pad6PartRAATheory3->cd();
		histo2DRAAAll2Up->GetXaxis()->SetRangeUser(-0.25,19.5);
		histo2DRAAAll2Up->DrawCopy();
		
		graphRAASysCombInd0510->Draw("E2same");
		graphRAACombInd0510->Draw("p,same,e1");
		
	//       DrawGammaSetMarkerTGraphErr(gEPOS_RAA_0510, markerStyleCommmonSpectrum0510,markerSizeSpectrum, colorEPOS0510, colorEPOS0510,2);
	//       gEPOS_RAA_0510->SetFillStyle(fillStyleEPOS);
	//       gEPOS_RAA_0510->SetFillColor(colorEPOS0510);
	//       gEPOS_RAA_0510->Draw("3 same");  
	//       DrawGammaSetMarkerTGraphErr(Xiao_Raa_0510, markerStyleCommmonSpectrum0510,markerSizeSpectrum, colorXiao0510 ,colorXiao0510,2);
	//       Xiao_Raa_0510->SetFillStyle(fillStyleXiao); 
	//       Xiao_Raa_0510->SetFillColor(colorXiao0510);
	//       Xiao_Raa_0510->Draw("3 same");
		DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_0510, markerStyleCommmonSpectrum0510,markerSizeSpectrum, colorVitevBas0510 ,colorVitevBas0510,2);
		Vitev_Bas_Raa_0510->SetFillStyle(fillStyleVitev);
		Vitev_Bas_Raa_0510->SetFillColor(colorVitevBas0510);
		Vitev_Bas_Raa_0510->Draw("3 same");
		DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_0510, markerStyleCommmonSpectrum0510,markerSizeSpectrum, colorWHDG0510, colorWHDG0510,0.8);
		gWHDG_Raa_0510->SetFillStyle(fillStyleWHDG);
		gWHDG_Raa_0510->SetFillColor(colorWHDG0510);
		gWHDG_Raa_0510->Draw("3 same");
		graphRAACombInd0510->Draw("p,same,e1");
	//       gEPOS_RAA_0510->Draw("3 same");  
		TLatex *labelRAAPi0TheoryPbPb0510 = new TLatex(0.05,0.9,collisionSystemCent0510.Data());
		SetStyleTLatex( labelRAAPi0TheoryPbPb0510, 0.85*textsizeLabelsRatioUp,4);
		labelRAAPi0TheoryPbPb0510->Draw();

	TLegend* legendRAAPi0TheoryPbPb0510 = new TLegend(0.05,0.79,0.76,0.86);
	legendRAAPi0TheoryPbPb0510->SetFillColor(0);
	legendRAAPi0TheoryPbPb0510->SetLineColor(0);
	legendRAAPi0TheoryPbPb0510->SetNColumns(2);
	legendRAAPi0TheoryPbPb0510->SetMargin(0.25);
	legendRAAPi0TheoryPbPb0510->SetTextSize(0.85*textsizeLabelsRatioUp);
	legendRAAPi0TheoryPbPb0510->SetTextFont(42);
	legendRAAPi0TheoryPbPb0510->AddEntry(graphRAASysCombInd0510,"#pi^{0} ALICE","pf");
	legendRAAPi0TheoryPbPb0510->AddEntry((TObject*)0, "", "");  
	legendRAAPi0TheoryPbPb0510->Draw();
	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm0510_Theory = CreateBoxConv(colorComb0510Box, 18., 1.-normErr0510 , 18.5, 1.+normErr0510);
	boxErrorNorm0510_Theory->Draw();

	histo2DRAAAll2Up->Draw("axis,same");

	pad6PartRAATheory3->Update();
	pad6PartRAATheory4->cd();
		histo2DRAAAll2Down->GetXaxis()->SetRangeUser(-0.25,19.5);
		histo2DRAAAll2Down->DrawCopy();

		graphRAASysCombInd4060->Draw("E2same");
		graphRAACombInd4060->Draw("p,same,e1");

	//       DrawGammaSetMarkerTGraphErr(gEPOS_RAA_4060, markerStyleCommmonSpectrum4060,markerSizeSpectrum, colorEPOS4060, colorEPOS4060,2);
	//       gEPOS_RAA_4060->SetFillStyle(fillStyleEPOS);
	//       gEPOS_RAA_4060->SetFillColor(colorEPOS4060);
	//       gEPOS_RAA_4060->Draw("3 same");  
	//       DrawGammaSetMarkerTGraphErr(Xiao_Raa_4060, markerStyleCommmonSpectrum4060,markerSizeSpectrum, colorXiao4060 ,colorXiao4060,2);
	//       Xiao_Raa_4060->SetFillStyle(fillStyleXiao); 
	//       Xiao_Raa_4060->SetFillColor(colorXiao4060);
	//       Xiao_Raa_4060->Draw("3 same");   
		DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_4060, markerStyleCommmonSpectrum4060,markerSizeSpectrum, colorVitevBas4060 ,colorVitevBas4060,2);
		Vitev_Bas_Raa_4060->SetFillStyle(fillStyleVitev);
		Vitev_Bas_Raa_4060->SetFillColor(colorVitevBas4060);
		Vitev_Bas_Raa_4060->Draw("3 same");
		DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_4060, markerStyleCommmonSpectrum4060,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,0.8);
		gWHDG_Raa_4060->SetFillStyle(fillStyleWHDG);
		gWHDG_Raa_4060->SetFillColor(colorWHDG4060);
		gWHDG_Raa_4060->Draw("3 same");
		graphRAACombInd4060->Draw("p,same,e1");
	//       gEPOS_RAA_4060->Draw("3 same");  
		TLatex *labelRAAPi0TheoryPbPb4060 = new TLatex(0.05,0.93,collisionSystemCent4060.Data());
		SetStyleTLatex( labelRAAPi0TheoryPbPb4060, 0.85*textsizeLabelsRatioDown,4);
		labelRAAPi0TheoryPbPb4060->Draw();
		TLegend* legendRAAPi0TheoryPbPb4060 = new TLegend(0.05,0.85,0.76,0.91);
		legendRAAPi0TheoryPbPb4060->SetFillColor(0);
		legendRAAPi0TheoryPbPb4060->SetLineColor(0);
		legendRAAPi0TheoryPbPb4060->SetTextSize(0.85*textsizeLabelsRatioDown);
		legendRAAPi0TheoryPbPb4060->SetNColumns(2);
		legendRAAPi0TheoryPbPb4060->SetTextFont(42);
		legendRAAPi0TheoryPbPb4060->SetMargin(0.25);
		legendRAAPi0TheoryPbPb4060->AddEntry(graphRAASysCombInd4060,"#pi^{0} ALICE","pf");
		legendRAAPi0TheoryPbPb4060->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		
		TBox* boxErrorNorm4060_Theory = CreateBoxConv(colorComb4060Box, 18., 1.-normErr4060 , 18.5, 1.+normErr4060);
		boxErrorNorm4060_Theory->Draw();
		histo2DRAAAll2Down->Draw("axis,same");
		
	pad6PartRAATheory4->Update();

	pad6PartRAATheory5->cd();
	histo2DRAAAll2Up->DrawCopy();

		graphRAASysCombInd1020->Draw("E2same");
		graphRAACombInd1020->Draw("p,same,e1");

	//       DrawGammaSetMarkerTGraphErr(gEPOS_RAA_1020, markerStyleCommmonSpectrum1020,markerSizeSpectrum, colorEPOS1020, colorEPOS1020,2);
	//       gEPOS_RAA_1020->SetFillStyle(fillStyleEPOS);         
	//       gEPOS_RAA_1020->SetFillColor(colorEPOS1020);
	//       gEPOS_RAA_1020->Draw("3 same");  
	//       DrawGammaSetMarkerTGraphErr(Xiao_Raa_1020, markerStyleCommmonSpectrum1020,markerSizeSpectrum, colorXiao1020 ,colorXiao1020,2);
	//       Xiao_Raa_1020->SetFillStyle(fillStyleXiao); 
	//       Xiao_Raa_1020->SetFillColor(colorXiao1020);
	//       Xiao_Raa_1020->Draw("3 same");
		DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_1020, markerStyleCommmonSpectrum1020,markerSizeSpectrum, colorVitevBas1020 ,colorVitevBas1020,2);
		Vitev_Bas_Raa_1020->SetFillStyle(fillStyleVitev);
		Vitev_Bas_Raa_1020->SetFillColor(colorVitevBas1020);
		Vitev_Bas_Raa_1020->Draw("3 same");
		DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_1020, markerStyleCommmonSpectrum1020,markerSizeSpectrum, colorWHDG1020, colorWHDG1020,0.8);
		gWHDG_Raa_1020->SetFillStyle(fillStyleWHDG);
		gWHDG_Raa_1020->SetFillColor(colorWHDG1020);
		gWHDG_Raa_1020->Draw("3 same");
		graphRAACombInd1020->Draw("p,same,e1");
	//       gEPOS_RAA_1020->Draw("3 same");  
		
		TLatex *labelRAAPi0TheoryPbPb1020 = new TLatex(0.05,0.9,collisionSystemCent1020.Data());
		SetStyleTLatex( labelRAAPi0TheoryPbPb1020, 0.85*textsizeLabelsRatioUp,4);
		labelRAAPi0TheoryPbPb1020->Draw();
		
		TLegend* legendRAAPi0TheoryPbPb1020 = new TLegend(0.05,0.79,0.76,0.86);
		legendRAAPi0TheoryPbPb1020->SetFillColor(0);
		legendRAAPi0TheoryPbPb1020->SetLineColor(0);
		legendRAAPi0TheoryPbPb1020->SetNColumns(2);
		legendRAAPi0TheoryPbPb1020->SetTextFont(42);
		legendRAAPi0TheoryPbPb1020->SetMargin(0.25);
		legendRAAPi0TheoryPbPb1020->SetTextSize(0.85*textsizeLabelsRatioUp);
		legendRAAPi0TheoryPbPb1020->AddEntry(graphRAASysCombInd1020,"#pi^{0} ALICE","pf");
		legendRAAPi0TheoryPbPb1020->AddEntry((TObject*)0, "", "");  
		legendRAAPi0TheoryPbPb1020->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		TBox* boxErrorNorm1020_Theory = CreateBoxConv(colorComb1020Box, 18., 1.-normErr1020 , 18.5, 1.+normErr1020);
		boxErrorNorm1020_Theory->Draw();

		histo2DRAAAll2Up->Draw("axis,same");

	pad6PartRAATheory5->Update();

	pad6PartRAATheory6->cd();
		histo2DRAAAll2Down->DrawCopy();

		
		TLegend* legendRAAPi0TheoryPbPb6080 = new TLegend(0.05,0.85,0.76,0.91);
		legendRAAPi0TheoryPbPb6080->SetFillColor(0);
		legendRAAPi0TheoryPbPb6080->SetLineColor(0);
		legendRAAPi0TheoryPbPb6080->SetTextSize(0.85*textsizeLabelsRatioDown);
		legendRAAPi0TheoryPbPb6080->SetNColumns(2);
		legendRAAPi0TheoryPbPb6080->SetTextFont(42);
		legendRAAPi0TheoryPbPb6080->SetMargin(0.25);
		legendRAAPi0TheoryPbPb6080->AddEntry(graphRAASysCombInd6080,"#pi^{0} ALICE","pf");
		legendRAAPi0TheoryPbPb6080->Draw();
		
		graphRAASysCombInd6080->Draw("E2same");
		graphRAACombInd6080->Draw("p,same,e1");
		
	//       DrawGammaSetMarkerTGraphErr(gEPOS_RAA_6080, markerStyleCommmonSpectrum6080,markerSizeSpectrum, colorEPOS6080, colorEPOS6080,2);
	//       gEPOS_RAA_6080->SetFillStyle(fillStyleEPOS);
	//       gEPOS_RAA_6080->SetFillColor(colorEPOS6080);
	//       gEPOS_RAA_6080->Draw("3 same");     
	//       DrawGammaSetMarkerTGraphErr(Xiao_Raa_6080, markerStyleCommmonSpectrum6080,markerSizeSpectrum, colorXiao6080 ,colorXiao6080,2);
	//       Xiao_Raa_6080->SetFillStyle(fillStyleXiao); 
	//       Xiao_Raa_6080->SetFillColor(colorXiao6080);
	//       Xiao_Raa_6080->Draw("3 same");   
		DrawGammaSetMarkerTGraphErr(Vitev_Bas_Raa_6080, markerStyleCommmonSpectrum6080,markerSizeSpectrum, colorVitevBas6080 ,colorVitevBas6080,2);
		Vitev_Bas_Raa_6080->SetFillStyle(fillStyleVitev);
		Vitev_Bas_Raa_6080->SetFillColor(colorVitevBas6080);
		Vitev_Bas_Raa_6080->Draw("3 same");
		DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_6080, markerStyleCommmonSpectrum6080,markerSizeSpectrum, colorWHDG6080, colorWHDG6080,0.8);
		gWHDG_Raa_6080->SetFillStyle(fillStyleWHDG);
		gWHDG_Raa_6080->SetFillColor(colorWHDG6080);
		gWHDG_Raa_6080->Draw("3 same");
		graphRAACombInd6080->Draw("p,same,e1");

	//       DrawGammaSetMarkerTGraph(gKopeliovichTotal_RAA_6080, markerStyleCommmonSpectrum6080,markerSizeSpectrum, colorEPOS6080, colorEPOS6080,2);
	//       gKopeliovichTotal_RAA_6080->SetFillStyle(fillStyleEPOS);
	//       gKopeliovichTotal_RAA_6080->SetFillColor(colorEPOS6080);
	//       gKopeliovichTotal_RAA_6080->Draw("same");     
		

		TLatex *labelRAAPi0TheoryPbPb6080 = new TLatex(0.05,0.93,collisionSystemCent6080.Data());
		SetStyleTLatex( labelRAAPi0TheoryPbPb6080, 0.85*textsizeLabelsRatioDown,4);
		labelRAAPi0TheoryPbPb6080->Draw();
		
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		TBox* boxErrorNorm6080_Theory = CreateBoxConv(colorComb6080Box, 18., 1.-normErr6080 , 18.5, 1.+normErr6080);
		boxErrorNorm6080_Theory->Draw();
		histo2DRAAAll2Down->Draw("axis,same");
		
	pad6PartRAATheory6->Update();


	canvas6PartRAATheory->Update();  
	canvas6PartRAATheory->SaveAs(Form("%s/RAA_All6PartedTheory_Paper_LinY.%s",outputDir.Data(),suffix.Data()));


	//**************************************************************************************************************
	//*********************************** Comparison to Theory *****************************************************
	//**************************************************************************************************************
	canvas6PartRAATheory->cd();
	pad6PartRAATheory1->Draw();
	pad6PartRAATheory2->Draw();

	pad6PartRAATheory3->Draw();
	pad6PartRAATheory4->Draw();

	pad6PartRAATheory5->Draw();
	pad6PartRAATheory6->Draw();


	pad6PartRAATheory1->cd();
	pad6PartRAATheory1->SetLogy();
		histo2DRAAAll3Up->GetXaxis()->SetRangeUser(0.,19.5);
		histo2DRAAAll3Up->GetYaxis()->SetRangeUser(0.05,3.);
		histo2DRAAAll3Up->DrawCopy();
		graphRAASysCombInd0005->Draw("E2same");
		graphRAACombInd0005->Draw("p,same,e1");
	//       Xiao_Raa_0005->Draw("3 same");
		Vitev_Bas_Raa_0005->Draw("3 same");
		gWHDG_Raa_0005->Draw("3 same");
		graphRAACombInd0005->Draw("p,same,e1");
	//       gEPOS_RAA_0005->Draw("3 same");
		labelRAAPi0TheoryPbPb0005->Draw();
		legendRAAPi0TheoryPbPb0005->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		boxErrorNorm0005_Theory->Draw();
	//       legendRAAPi0TheoryPbPb->SetY1NDC(0.15);
	//       legendRAAPi0TheoryPbPb->SetY2NDC(0.27);
		legendRAAPi0TheoryPbPb->Draw();

		histo2DRAAAll3Up->Draw("axis,same");

	pad6PartRAATheory1->Update();

	pad6PartRAATheory2->cd();
	pad6PartRAATheory2->SetLogy();
		histo2DRAAAll3Down->GetXaxis()->SetRangeUser(0.,19.5);
		histo2DRAAAll3Down->GetYaxis()->SetRangeUser(0.05,3.);
		histo2DRAAAll3Down->DrawCopy();
		graphRAASysCombInd2040->Draw("E2same");
		graphRAACombInd2040->Draw("p,same,e1");
		
	//       Xiao_Raa_2040->Draw("3 same");	
		Vitev_Bas_Raa_2040->Draw("3 same");
		gWHDG_Raa_2040->Draw("3 same");
		graphRAACombInd2040->Draw("p,same,e1");
	//       gEPOS_RAA_2040->Draw("3 same");
		
		labelRAAPi0TheoryPbPb2040->Draw();
		legendRAAPi0TheoryPbPb2040->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		boxErrorNorm2040_Theory->Draw();
		histo2DRAAAll3Down->Draw("axis,same");

		
	pad6PartRAATheory2->Update();

	pad6PartRAATheory3->cd();
	pad6PartRAATheory3->SetLogy();
		histo2DRAAAll2Up->GetXaxis()->SetRangeUser(-0.25,19.5);
		histo2DRAAAll2Up->GetYaxis()->SetRangeUser(0.05,3.);
		histo2DRAAAll2Up->DrawCopy();
		
		graphRAASysCombInd0510->Draw("E2same");
		graphRAACombInd0510->Draw("p,same,e1");
		
	//       Xiao_Raa_0510->Draw("3 same");
		Vitev_Bas_Raa_0510->Draw("3 same");
		gWHDG_Raa_0510->Draw("3 same");
		graphRAACombInd0510->Draw("p,same,e1");
	//       gEPOS_RAA_0510->Draw("3 same");
		
		labelRAAPi0TheoryPbPb0510->Draw();
		legendRAAPi0TheoryPbPb0510->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		boxErrorNorm0510_Theory->Draw();
		histo2DRAAAll2Up->Draw("axis,same");

	pad6PartRAATheory3->Update();
	pad6PartRAATheory4->cd();
	pad6PartRAATheory4->SetLogy();

		histo2DRAAAll2Down->GetYaxis()->SetRangeUser(0.05,3.);
		histo2DRAAAll2Down->DrawCopy();
		graphRAASysCombInd4060->Draw("E2same");
		graphRAACombInd4060->Draw("p,same,e1");
		
	//       Xiao_Raa_4060->Draw("3 same");	
		Vitev_Bas_Raa_4060->Draw("3 same");
		gWHDG_Raa_4060->Draw("3 same");
		graphRAACombInd4060->Draw("p,same,e1");
	//       gEPOS_RAA_4060->Draw("3 same");
		
		labelRAAPi0TheoryPbPb4060->Draw();
		legendRAAPi0TheoryPbPb4060->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		histo2DRAAAll2Down->Draw("axis,same");

	pad6PartRAATheory4->Update();

	pad6PartRAATheory5->cd();
	pad6PartRAATheory5->SetLogy();

		histo2DRAAAll2Up->GetYaxis()->SetRangeUser(0.05,3.);
		histo2DRAAAll2Up->DrawCopy();
		graphRAASysCombInd1020->Draw("E2same");
		graphRAACombInd1020->Draw("p,same,e1");
		
	//       Xiao_Raa_1020->Draw("3 same");
		Vitev_Bas_Raa_1020->Draw("3 same");
		gWHDG_Raa_1020->Draw("3 same");
		graphRAACombInd1020->Draw("p,same,e1");
	//       gEPOS_RAA_1020->Draw("3 same");
		
		labelRAAPi0TheoryPbPb1020->Draw();
		legendRAAPi0TheoryPbPb1020->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		boxErrorNorm1020_Theory->Draw();
		histo2DRAAAll2Up->Draw("axis,same");

	pad6PartRAATheory5->Update();

	pad6PartRAATheory6->cd();
	pad6PartRAATheory6->SetLogy();
		histo2DRAAAll2Down->GetYaxis()->SetRangeUser(0.05,3.);
		histo2DRAAAll2Down->DrawCopy();

		legendRAAPi0TheoryPbPb6080->Draw();
		graphRAASysCombInd6080->Draw("E2same");
		graphRAACombInd6080->Draw("p,same,e1");
		
	//       Xiao_Raa_6080->Draw("3 same");	
		Vitev_Bas_Raa_6080->Draw("3 same");
		gWHDG_Raa_6080->Draw("3 same");
		graphRAACombInd6080->Draw("p,same,e1");
	//       gEPOS_RAA_6080->Draw("3 same");
		
		labelRAAPi0TheoryPbPb6080->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		histo2DRAAAll2Down->Draw("axis,same");

	pad6PartRAATheory6->Update();

	canvas6PartRAATheory->Update();	
	canvas6PartRAATheory->SaveAs(Form("%s/RAA_All6PartedTheory_Paper_LogY.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartRAATheory1;	
	delete pad6PartRAATheory2;	
	delete pad6PartRAATheory3;	
	delete pad6PartRAATheory4;	
	delete pad6PartRAATheory5;	
	delete pad6PartRAATheory6;	
	delete canvas6PartRAATheory;

	TCanvas * canvas6PartRAACompChargedHadrons = new TCanvas("canvas6PartRAATheory","",10,10,2400,1300);  // gives the page size		
	canvas6PartRAACompChargedHadrons->cd();
	DrawGammaCanvasSettings( canvas6PartRAACompChargedHadrons, 0.13, 0.0, 0.02, 0.09);

	TPad* pad6PartRAACompChargedHadrons1 = new TPad("pad6PartRAACompChargedHadrons1", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedHadrons1, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRAACompChargedHadrons1->Draw();
	TPad* pad6PartRAACompChargedHadrons2 = new TPad("pad6PartRAACompChargedHadrons2", "",arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedHadrons2, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRAACompChargedHadrons2->Draw();

	TPad* pad6PartRAACompChargedHadrons3 = new TPad("pad6PartRAACompChargedHadrons3", "",arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedHadrons3, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRAACompChargedHadrons3->Draw();
	TPad* pad6PartRAACompChargedHadrons4 = new TPad("pad6PartRAACompChargedHadrons4", "",arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedHadrons4, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRAACompChargedHadrons4->Draw();

	TPad* pad6PartRAACompChargedHadrons5 = new TPad("pad6PartRAACompChargedHadrons5", "",arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedHadrons5, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRAACompChargedHadrons5->Draw();
	TPad* pad6PartRAACompChargedHadrons6 = new TPad("pad6PartRAACompChargedHadrons6", "", arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedHadrons6, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRAACompChargedHadrons6->Draw();

	pad6PartRAACompChargedHadrons1->cd();
	//  	pad6PartRAACompChargedHadrons1->SetLogy();
	histo2DRAAAll3Up->GetXaxis()->SetRangeUser(0.,19.5);
	histo2DRAAAll3Up->GetYaxis()->SetRangeUser(0.05,1.4);
	// 	histo2DRAAAll3Up->GetYaxis()->SetRangeUser(0.05,3.);
	histo2DRAAAll3Up->DrawCopy();

		
	graphRAASysCombInd0005->Draw("E2same");
	DrawGammaSetMarkerTGraphErr(graphChargedRAAStat_0005, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphChargedRAASys_0005, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	graphChargedRAAStat_0005->Draw("p,same,e1");
	graphChargedRAASys_0005->Draw("E[]same");
	graphRAACombInd0005->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedHadronsPbPb0005 = new TLatex(0.18,0.9,collisionSystemCent0005.Data());
	SetStyleTLatex( labelRAAPi0CompChargedHadronsPbPb0005, 0.85*textsizeLabelsRatioUp,4);
	labelRAAPi0CompChargedHadronsPbPb0005->Draw();


	TLegend* legendRAAPi0CompChargedHadronsPbPb0005 = new TLegend(0.18,0.82,0.9,0.88);
	legendRAAPi0CompChargedHadronsPbPb0005->SetFillColor(0);
	legendRAAPi0CompChargedHadronsPbPb0005->SetLineColor(0);
	legendRAAPi0CompChargedHadronsPbPb0005->SetNColumns(2);
	legendRAAPi0CompChargedHadronsPbPb0005->SetTextFont(42);
	legendRAAPi0CompChargedHadronsPbPb0005->SetTextSize(0.85*textsizeLabelsRatioUp);
	legendRAAPi0CompChargedHadronsPbPb0005->AddEntry(graphRAASysCombInd0005,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedHadronsPbPb0005->AddEntry(graphChargedRAAStat_0005,"h^{#pm} ALICE","p");
	legendRAAPi0CompChargedHadronsPbPb0005->Draw();

	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm0005_CompChargedHadrons = CreateBoxConv(colorComb0005Box, 18., 1.-normErr0005 , 18.5, 1.+normErr0005);
	boxErrorNorm0005_CompChargedHadrons->Draw();

	histo2DRAAAll3Up->Draw("axis,same");

	pad6PartRAACompChargedHadrons1->Update();
	pad6PartRAACompChargedHadrons2->cd();
	//  	pad6PartRAACompChargedHadrons2->SetLogy();
	histo2DRAAAll3Down->GetYaxis()->SetRangeUser(0.05,1.4);
	// 	histo2DRAAAll3Down->GetYaxis()->SetRangeUser(0.05,3.);
	histo2DRAAAll3Down->DrawCopy();
	graphRAASysCombInd2040->Draw("E2same");
	DrawGammaSetMarkerTGraphErr(graphChargedRAAStat_2040, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphChargedRAASys_2040, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	graphChargedRAAStat_2040->Draw("p,same,e1");
	graphChargedRAASys_2040->Draw("E[]same");
	graphRAACombInd2040->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedHadronsPbPb2040 = new TLatex(0.18,0.93,collisionSystemCent2040.Data());
	SetStyleTLatex( labelRAAPi0CompChargedHadronsPbPb2040, 0.85*textsizeLabelsRatioDown,4);
	labelRAAPi0CompChargedHadronsPbPb2040->Draw();

	TLegend* legendRAAPi0CompChargedHadronsPbPb2040 = new TLegend(0.18,0.85,0.76,0.91);
	legendRAAPi0CompChargedHadronsPbPb2040->SetFillColor(0);
	legendRAAPi0CompChargedHadronsPbPb2040->SetLineColor(0);
	legendRAAPi0CompChargedHadronsPbPb2040->SetNColumns(2);
	legendRAAPi0CompChargedHadronsPbPb2040->SetTextFont(42);
	legendRAAPi0CompChargedHadronsPbPb2040->SetTextSize(0.85*textsizeLabelsRatioDown);
	legendRAAPi0CompChargedHadronsPbPb2040->AddEntry(graphRAASysCombInd2040,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedHadronsPbPb2040->Draw();
	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm2040_CompChargedHadrons = CreateBoxConv(colorComb2040Box, 18., 1.-normErr2040 , 18.5, 1.+normErr2040);
	boxErrorNorm2040_CompChargedHadrons->Draw();


	histo2DRAAAll3Down->Draw("axis,same");
	pad6PartRAACompChargedHadrons2->Update();

	pad6PartRAACompChargedHadrons3->cd();
	// 	pad6PartRAACompChargedHadrons3->SetLogy();
	histo2DRAAAll2Up->GetXaxis()->SetRangeUser(-0.25,19.5);
	histo2DRAAAll2Up->GetYaxis()->SetRangeUser(0.05,1.4);
	// 	histo2DRAAAll2->GetYaxis()->SetRangeUser(0.05,3.);
	histo2DRAAAll2Up->DrawCopy();
			
	graphRAASysCombInd0510->Draw("E2same");
	DrawGammaSetMarkerTGraphErr(graphChargedRAAStat_0510, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphChargedRAASys_0510, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	graphChargedRAAStat_0510->Draw("p,same,e1");
	graphChargedRAASys_0510->Draw("E[]same");
	graphRAACombInd0510->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedHadronsPbPb0510 = new TLatex(0.05,0.9,collisionSystemCent0510.Data());
	SetStyleTLatex( labelRAAPi0CompChargedHadronsPbPb0510, 0.85*textsizeLabelsRatioUp,4);
	labelRAAPi0CompChargedHadronsPbPb0510->Draw();


	TLegend* legendRAAPi0CompChargedHadronsPbPb0510 = new TLegend(0.06,0.82,0.76,0.88);
	legendRAAPi0CompChargedHadronsPbPb0510->SetFillColor(0);
	legendRAAPi0CompChargedHadronsPbPb0510->SetLineColor(0);
	legendRAAPi0CompChargedHadronsPbPb0510->SetNColumns(2);
	legendRAAPi0CompChargedHadronsPbPb0510->SetTextFont(42);
	legendRAAPi0CompChargedHadronsPbPb0510->SetTextSize(0.85*textsizeLabelsRatioUp);
	legendRAAPi0CompChargedHadronsPbPb0510->AddEntry(graphRAASysCombInd0510,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedHadronsPbPb0510->Draw();

	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm0510_CompChargedHadrons = CreateBoxConv(colorComb0510Box, 18., 1.-normErr0510 , 18.5, 1.+normErr0510);
	boxErrorNorm0510_CompChargedHadrons->Draw();
	histo2DRAAAll2Up->Draw("axis,same");

	pad6PartRAACompChargedHadrons3->Update();
	pad6PartRAACompChargedHadrons4->cd();
	//  	pad6PartRAACompChargedHadrons4->SetLogy();
	histo2DRAAAll2Down->GetYaxis()->SetRangeUser(0.05,1.4);
	// 	histo2DRAAAll2Down->GetYaxis()->SetRangeUser(0.05,3.);
	histo2DRAAAll2Down->DrawCopy();

	graphRAASysCombInd4060->Draw("E2same");
	DrawGammaSetMarkerTGraphErr(graphChargedRAAStat_4060, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphChargedRAASys_4060, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	graphChargedRAAStat_4060->Draw("p,same,e1");
	graphChargedRAASys_4060->Draw("E[]same");
	graphRAACombInd4060->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedHadronsPbPb4060 = new TLatex(0.05,0.93,collisionSystemCent4060.Data());
	SetStyleTLatex( labelRAAPi0CompChargedHadronsPbPb4060, 0.85*textsizeLabelsRatioDown,4);
	labelRAAPi0CompChargedHadronsPbPb4060->Draw();

	TLegend* legendRAAPi0CompChargedHadronsPbPb4060 = new TLegend(0.06,0.85,0.9,0.91);
	legendRAAPi0CompChargedHadronsPbPb4060->SetFillColor(0);
	legendRAAPi0CompChargedHadronsPbPb4060->SetLineColor(0);
	legendRAAPi0CompChargedHadronsPbPb4060->SetTextSize(0.85*textsizeLabelsRatioDown);
	legendRAAPi0CompChargedHadronsPbPb4060->SetNColumns(2);
	legendRAAPi0CompChargedHadronsPbPb4060->SetTextFont(42);
	legendRAAPi0CompChargedHadronsPbPb4060->AddEntry(graphRAASysCombInd4060,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedHadronsPbPb4060->Draw();
	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

	TBox* boxErrorNorm4060_CompChargedHadrons = CreateBoxConv(colorComb4060Box, 18., 1.-normErr4060 , 18.5, 1.+normErr4060);
	boxErrorNorm4060_CompChargedHadrons->Draw();
	histo2DRAAAll2Down->Draw("axis,same");

	pad6PartRAACompChargedHadrons4->Update();

	pad6PartRAACompChargedHadrons5->cd();
	histo2DRAAAll2Up->GetYaxis()->SetRangeUser(0.05,1.4);
	histo2DRAAAll2Up->DrawCopy();

	graphRAASysCombInd1020->Draw("E2same");
	DrawGammaSetMarkerTGraphErr(graphChargedRAAStat_1020, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphChargedRAASys_1020, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	graphChargedRAAStat_1020->Draw("p,same,e1");
	graphChargedRAASys_1020->Draw("E[]same");
	graphRAACombInd1020->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedHadronsPbPb1020 = new TLatex(0.05,0.9,collisionSystemCent1020.Data());
	SetStyleTLatex( labelRAAPi0CompChargedHadronsPbPb1020, 0.85*textsizeLabelsRatioUp,4);
	labelRAAPi0CompChargedHadronsPbPb1020->Draw();


	TLegend* legendRAAPi0CompChargedHadronsPbPb1020 = new TLegend(0.06,0.82,0.76,0.88);
	legendRAAPi0CompChargedHadronsPbPb1020->SetFillColor(0);
	legendRAAPi0CompChargedHadronsPbPb1020->SetLineColor(0);
	legendRAAPi0CompChargedHadronsPbPb1020->SetNColumns(2);
	legendRAAPi0CompChargedHadronsPbPb1020->SetTextFont(42);
	legendRAAPi0CompChargedHadronsPbPb1020->SetTextSize(0.85*textsizeLabelsRatioUp);
	legendRAAPi0CompChargedHadronsPbPb1020->AddEntry(graphRAASysCombInd1020,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedHadronsPbPb1020->Draw();

	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm1020_CompChargedHadrons = CreateBoxConv(colorComb1020Box, 18., 1.-normErr1020 , 18.5, 1.+normErr1020);
	boxErrorNorm1020_CompChargedHadrons->Draw();
	histo2DRAAAll2Up->Draw("axis,same");

	pad6PartRAACompChargedHadrons5->Update();

	pad6PartRAACompChargedHadrons6->cd();
	//  	pad6PartRAACompChargedHadrons4->SetLogy();
	histo2DRAAAll2Down->GetYaxis()->SetRangeUser(0.05,1.4);
	// 	histo2DRAAAll2Down->GetYaxis()->SetRangeUser(0.05,3.);
	histo2DRAAAll2Down->DrawCopy();

	graphRAASysCombInd6080->Draw("E2same");
	DrawGammaSetMarkerTGraphErr(graphChargedRAAStat_6080, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphChargedRAASys_6080, 24,markerSizeChargedHadronSpectrum, kBlack , kBlack);
	graphChargedRAAStat_6080->Draw("p,same,e1");
	graphChargedRAASys_6080->Draw("E[]same");
	graphRAACombInd6080->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedHadronsPbPb6080 = new TLatex(0.05,0.93,collisionSystemCent6080.Data());
	SetStyleTLatex( labelRAAPi0CompChargedHadronsPbPb6080, 0.85*textsizeLabelsRatioDown,4);
	labelRAAPi0CompChargedHadronsPbPb6080->Draw();

	TLegend* legendRAAPi0CompChargedHadronsPbPb6080 = new TLegend(0.06,0.85,0.78,0.91);
	legendRAAPi0CompChargedHadronsPbPb6080->SetFillColor(0);
	legendRAAPi0CompChargedHadronsPbPb6080->SetLineColor(0);
	legendRAAPi0CompChargedHadronsPbPb6080->SetTextSize(0.85*textsizeLabelsRatioDown);
	legendRAAPi0CompChargedHadronsPbPb6080->SetNColumns(2);
	legendRAAPi0CompChargedHadronsPbPb6080->SetTextFont(42);
	legendRAAPi0CompChargedHadronsPbPb6080->AddEntry(graphRAASysCombInd6080,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedHadronsPbPb6080->Draw();
	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm6080_CompChargedHadrons = CreateBoxConv(colorComb6080Box, 18., 1.-normErr6080 , 18.5, 1.+normErr6080);
	boxErrorNorm6080_CompChargedHadrons->Draw();
	histo2DRAAAll2Down->Draw("axis,same");

	pad6PartRAACompChargedHadrons6->Update();

	canvas6PartRAACompChargedHadrons->Update();	
	canvas6PartRAACompChargedHadrons->SaveAs(Form("%s/RAA_All6PartedCompChargedHadrons_Paper.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartRAACompChargedHadrons1;	
	delete pad6PartRAACompChargedHadrons2;	
	delete pad6PartRAACompChargedHadrons3;	
	delete pad6PartRAACompChargedHadrons4;	
	delete canvas6PartRAACompChargedHadrons;

	cout << "here ............................." << endl;
		
	// ***************************************************************************************************************
	// ************************************ Charged Pion RAA *********************************************************
	// ***************************************************************************************************************

	TCanvas * canvas6PartRAACompChargedPions = new TCanvas("canvas6PartRAATheory","",10,10,2400,1300);  // gives the page size		
	canvas6PartRAACompChargedPions->cd();
	DrawGammaCanvasSettings( canvas6PartRAACompChargedPions, 0.13, 0.0, 0.02, 0.09);

	TPad* pad6PartRAACompChargedPions1 = new TPad("pad6PartRAACompChargedPions1", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedPions1, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRAACompChargedPions1->Draw();
	TPad* pad6PartRAACompChargedPions2 = new TPad("pad6PartRAACompChargedPions2", "",arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedPions2, relativeMarginsIndMeasRatioX[0], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRAACompChargedPions2->Draw();

	TPad* pad6PartRAACompChargedPions3 = new TPad("pad6PartRAACompChargedPions3", "",arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedPions3, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRAACompChargedPions3->Draw();
	TPad* pad6PartRAACompChargedPions4 = new TPad("pad6PartRAACompChargedPions4", "",arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedPions4, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRAACompChargedPions4->Draw();

	TPad* pad6PartRAACompChargedPions5 = new TPad("pad6PartRAACompChargedPions5", "",arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedPions5, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[0], relativeMarginsIndMeasRatioY[1]);
	pad6PartRAACompChargedPions5->Draw();
	TPad* pad6PartRAACompChargedPions6 = new TPad("pad6PartRAACompChargedPions6", "", arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRAACompChargedPions6, relativeMarginsIndMeasRatioX[1], relativeMarginsIndMeasRatioX[2], relativeMarginsIndMeasRatioY[1], relativeMarginsIndMeasRatioY[2]);
	pad6PartRAACompChargedPions6->Draw();

	pad6PartRAACompChargedPions1->cd();
	//  	pad6PartRAACompChargedPions1->SetLogy();
	histo2DRAAAll3Up->GetXaxis()->SetRangeUser(0.,19.5);
	histo2DRAAAll3Up->GetYaxis()->SetRangeUser(0.05,1.4);
	// 	histo2DRAAAll3Up->GetYaxis()->SetRangeUser(0.05,3.);
	histo2DRAAAll3Up->DrawCopy();

		
	graphRAASysCombInd0005->Draw("E2same");
	DrawGammaSetMarker(histoChargedPionRAAStat0005 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	DrawGammaSetMarker(histoChargedPionRAASyst0005 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	histoChargedPionRAAStat0005->Draw("p,same,e1");
	histoChargedPionRAASyst0005->Draw("E[]same");
	graphRAACombInd0005->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedPionsPbPb0005 = new TLatex(0.15,0.9,collisionSystemCent0005.Data());
	SetStyleTLatex( labelRAAPi0CompChargedPionsPbPb0005, 0.85*textsizeLabelsRatioUp,4);
	labelRAAPi0CompChargedPionsPbPb0005->Draw();


	TLegend* legendRAAPi0CompChargedPionsPbPb0005 = new TLegend(0.18,0.82,0.9,0.88);
	legendRAAPi0CompChargedPionsPbPb0005->SetFillColor(0);
	legendRAAPi0CompChargedPionsPbPb0005->SetLineColor(0);
	legendRAAPi0CompChargedPionsPbPb0005->SetNColumns(2);
	legendRAAPi0CompChargedPionsPbPb0005->SetTextFont(42);
	legendRAAPi0CompChargedPionsPbPb0005->SetTextSize(0.85*textsizeLabelsRatioUp);
	legendRAAPi0CompChargedPionsPbPb0005->AddEntry(graphRAASysCombInd0005,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedPionsPbPb0005->AddEntry(histoChargedPionRAAStat0005,"#pi^{#pm} ALICE","p");
	legendRAAPi0CompChargedPionsPbPb0005->Draw();

	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm0005_CompChargedPions = CreateBoxConv(colorComb0005Box, 18., 1.-normErr0005 , 18.5, 1.+normErr0005);
	boxErrorNorm0005_CompChargedPions->Draw();

	histo2DRAAAll3Up->Draw("axis,same");

	pad6PartRAACompChargedPions1->Update();
	pad6PartRAACompChargedPions2->cd();
	//  	pad6PartRAACompChargedPions2->SetLogy();
	histo2DRAAAll3Down->GetYaxis()->SetRangeUser(0.05,1.4);
	// 	histo2DRAAAll3Down->GetYaxis()->SetRangeUser(0.05,3.);
	histo2DRAAAll3Down->DrawCopy();
	graphRAASysCombInd2040->Draw("E2same");
	DrawGammaSetMarker(histoChargedPionRAAStat2040 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	DrawGammaSetMarker(histoChargedPionRAASyst2040 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	histoChargedPionRAAStat2040->Draw("p,same,e1");
	histoChargedPionRAASyst2040->Draw("E[]same");
	graphRAACombInd2040->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedPionsPbPb2040 = new TLatex(0.15,0.93,collisionSystemCent2040.Data());
	SetStyleTLatex( labelRAAPi0CompChargedPionsPbPb2040, 0.85*textsizeLabelsRatioDown,4);
	labelRAAPi0CompChargedPionsPbPb2040->Draw();

	TLegend* legendRAAPi0CompChargedPionsPbPb2040 = new TLegend(0.18,0.85,0.76,0.91);
	legendRAAPi0CompChargedPionsPbPb2040->SetFillColor(0);
	legendRAAPi0CompChargedPionsPbPb2040->SetLineColor(0);
	legendRAAPi0CompChargedPionsPbPb2040->SetNColumns(2);
	legendRAAPi0CompChargedPionsPbPb2040->SetTextFont(42);
	legendRAAPi0CompChargedPionsPbPb2040->SetTextSize(0.85*textsizeLabelsRatioDown);
	legendRAAPi0CompChargedPionsPbPb2040->AddEntry(graphRAASysCombInd2040,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedPionsPbPb2040->Draw();
	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm2040_CompChargedPions = CreateBoxConv(colorComb2040Box, 18., 1.-normErr2040 , 18.5, 1.+normErr2040);
	boxErrorNorm2040_CompChargedPions->Draw();


	histo2DRAAAll3Down->Draw("axis,same");
	pad6PartRAACompChargedPions2->Update();

	pad6PartRAACompChargedPions3->cd();
	// 	pad6PartRAACompChargedPions3->SetLogy();
	histo2DRAAAll2Up->GetXaxis()->SetRangeUser(-0.25,19.5);
	histo2DRAAAll2Up->GetYaxis()->SetRangeUser(0.05,1.4);
	// 	histo2DRAAAll2Up->GetYaxis()->SetRangeUser(0.05,3.);
	histo2DRAAAll2Up->DrawCopy();
			
	graphRAASysCombInd0510->Draw("E2same");
	DrawGammaSetMarker(histoChargedPionRAAStat0510 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	DrawGammaSetMarker(histoChargedPionRAASyst0510 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	histoChargedPionRAAStat0510->Draw("p,same,e1");
	histoChargedPionRAASyst0510->Draw("E[]same");
	graphRAACombInd0510->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedPionsPbPb0510 = new TLatex(0.03,0.9,collisionSystemCent0510.Data());
	SetStyleTLatex( labelRAAPi0CompChargedPionsPbPb0510, 0.85*textsizeLabelsRatioUp,4);
	labelRAAPi0CompChargedPionsPbPb0510->Draw();


	TLegend* legendRAAPi0CompChargedPionsPbPb0510 = new TLegend(0.06,0.82,0.76,0.88);
	legendRAAPi0CompChargedPionsPbPb0510->SetFillColor(0);
	legendRAAPi0CompChargedPionsPbPb0510->SetLineColor(0);
	legendRAAPi0CompChargedPionsPbPb0510->SetNColumns(2);
	legendRAAPi0CompChargedPionsPbPb0510->SetTextFont(42);
	legendRAAPi0CompChargedPionsPbPb0510->SetTextSize(0.85*textsizeLabelsRatioUp);
	legendRAAPi0CompChargedPionsPbPb0510->AddEntry(graphRAASysCombInd0510,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedPionsPbPb0510->Draw();

	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm0510_CompChargedPions = CreateBoxConv(colorComb0510Box, 18., 1.-normErr0510 , 18.5, 1.+normErr0510);
	boxErrorNorm0510_CompChargedPions->Draw();
	histo2DRAAAll2Up->Draw("axis,same");

	pad6PartRAACompChargedPions3->Update();
	pad6PartRAACompChargedPions4->cd();
	//  	pad6PartRAACompChargedPions4->SetLogy();
	histo2DRAAAll2Down->GetYaxis()->SetRangeUser(0.05,1.4);
	// 	histo2DRAAAll2Down->GetYaxis()->SetRangeUser(0.05,3.);
	histo2DRAAAll2Down->DrawCopy();

	graphRAASysCombInd4060->Draw("E2same");
	DrawGammaSetMarker(histoChargedPionRAAStat4060 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	DrawGammaSetMarker(histoChargedPionRAASyst4060 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	histoChargedPionRAAStat4060->Draw("p,same,e1");
	histoChargedPionRAASyst4060->Draw("E[]same");
	graphRAACombInd4060->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedPionsPbPb4060 = new TLatex(0.03,0.93,collisionSystemCent4060.Data());
	SetStyleTLatex( labelRAAPi0CompChargedPionsPbPb4060, 0.85*textsizeLabelsRatioDown,4);
	labelRAAPi0CompChargedPionsPbPb4060->Draw();

	TLegend* legendRAAPi0CompChargedPionsPbPb4060 = new TLegend(0.06,0.85,0.9,0.91);
	legendRAAPi0CompChargedPionsPbPb4060->SetFillColor(0);
	legendRAAPi0CompChargedPionsPbPb4060->SetLineColor(0);
	legendRAAPi0CompChargedPionsPbPb4060->SetTextSize(0.85*textsizeLabelsRatioDown);
	legendRAAPi0CompChargedPionsPbPb4060->SetNColumns(2);
	legendRAAPi0CompChargedPionsPbPb4060->SetTextFont(42);
	legendRAAPi0CompChargedPionsPbPb4060->AddEntry(graphRAASysCombInd4060,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedPionsPbPb4060->Draw();
	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

	TBox* boxErrorNorm4060_CompChargedPions = CreateBoxConv(colorComb4060Box, 18., 1.-normErr4060 , 18.5, 1.+normErr4060);
	boxErrorNorm4060_CompChargedPions->Draw();
	histo2DRAAAll2Down->Draw("axis,same");

	pad6PartRAACompChargedPions4->Update();

	pad6PartRAACompChargedPions5->cd();
	histo2DRAAAll2Up->GetYaxis()->SetRangeUser(0.05,1.4);
	histo2DRAAAll2Up->DrawCopy();

	graphRAASysCombInd1020->Draw("E2same");
	DrawGammaSetMarker(histoChargedPionRAAStat1020 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	DrawGammaSetMarker(histoChargedPionRAASyst1020 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	histoChargedPionRAAStat1020->Draw("p,same,e1");
	histoChargedPionRAASyst1020->Draw("E[]same");
	graphRAACombInd1020->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedPionsPbPb1020 = new TLatex(0.03,0.9,collisionSystemCent1020.Data());
	SetStyleTLatex( labelRAAPi0CompChargedPionsPbPb1020, 0.85*textsizeLabelsRatioUp,4);
	labelRAAPi0CompChargedPionsPbPb1020->Draw();


	TLegend* legendRAAPi0CompChargedPionsPbPb1020 = new TLegend(0.06,0.82,0.76,0.88);
	legendRAAPi0CompChargedPionsPbPb1020->SetFillColor(0);
	legendRAAPi0CompChargedPionsPbPb1020->SetLineColor(0);
	legendRAAPi0CompChargedPionsPbPb1020->SetNColumns(2);
	legendRAAPi0CompChargedPionsPbPb1020->SetTextFont(42);
	legendRAAPi0CompChargedPionsPbPb1020->SetTextSize(0.85*textsizeLabelsRatioUp);
	legendRAAPi0CompChargedPionsPbPb1020->AddEntry(graphRAASysCombInd1020,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedPionsPbPb1020->Draw();

	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm1020_CompChargedPions = CreateBoxConv(colorComb1020Box, 18., 1.-normErr1020 , 18.5, 1.+normErr1020);
	boxErrorNorm1020_CompChargedPions->Draw();
	histo2DRAAAll2Up->Draw("axis,same");

	pad6PartRAACompChargedPions5->Update();

	pad6PartRAACompChargedPions6->cd();
	//  	pad6PartRAACompChargedPions4->SetLogy();
	histo2DRAAAll2Down->GetYaxis()->SetRangeUser(0.05,1.4);
	// 	histo2DRAAAll2Down->GetYaxis()->SetRangeUser(0.05,3.);
	histo2DRAAAll2Down->DrawCopy();

	graphRAASysCombInd6080->Draw("E2same");
	DrawGammaSetMarker(histoChargedPionRAAStat6080 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	DrawGammaSetMarker(histoChargedPionRAASyst6080 ,24,markerSizeChargedHadronSpectrum, kBlack , kBlack);					 
	histoChargedPionRAAStat6080->Draw("p,same,e1");
	histoChargedPionRAASyst6080->Draw("E[]same");
		graphRAACombInd6080->Draw("p,same,e1");

	TLatex *labelRAAPi0CompChargedPionsPbPb6080 = new TLatex(0.04,0.93,collisionSystemCent6080.Data());
	SetStyleTLatex( labelRAAPi0CompChargedPionsPbPb6080, 0.85*textsizeLabelsRatioDown,4);
	labelRAAPi0CompChargedPionsPbPb6080->Draw();

	TLegend* legendRAAPi0CompChargedPionsPbPb6080 = new TLegend(0.04,0.85,0.78,0.91);
	legendRAAPi0CompChargedPionsPbPb6080->SetFillColor(0);
	legendRAAPi0CompChargedPionsPbPb6080->SetLineColor(0);
	legendRAAPi0CompChargedPionsPbPb6080->SetTextSize(0.85*textsizeLabelsRatioDown);
	legendRAAPi0CompChargedPionsPbPb6080->SetNColumns(2);
	legendRAAPi0CompChargedPionsPbPb6080->SetTextFont(42);
	legendRAAPi0CompChargedPionsPbPb6080->AddEntry(graphRAASysCombInd6080,"#pi^{0} ALICE","pf");
	legendRAAPi0CompChargedPionsPbPb6080->Draw();
	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	TBox* boxErrorNorm6080_CompChargedPions = CreateBoxConv(colorComb6080Box, 18., 1.-normErr6080 , 18.5, 1.+normErr6080);
	boxErrorNorm6080_CompChargedPions->Draw();
	histo2DRAAAll2Down->Draw("axis,same");

	pad6PartRAACompChargedPions6->Update();

	canvas6PartRAACompChargedPions->Update();	
	canvas6PartRAACompChargedPions->SaveAs(Form("%s/RAA_All6PartedCompChargedPions_Paper.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartRAACompChargedPions1;	
	delete pad6PartRAACompChargedPions2;	
	delete pad6PartRAACompChargedPions3;	
	delete pad6PartRAACompChargedPions4;	
	delete canvas6PartRAACompChargedPions;
	cout << "here ............................." << endl;



	TCanvas* canvasRAA_0010 = new TCanvas("canvasRAA_0010","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasRAA_0010,  0.1, 0.01, 0.015, 0.1);

	// 	canvasRAA_0010->SetLogy();
	// 	canvasRAA_0010->SetLogx();
	SetStyleHistoTH2ForGraphs(histo2DRAAAll3Up, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA,0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.95,1., 512, 505); //#frac{#frac{1
	// 	histo2DRAAAll3->GetYaxis()->SetRangeUser(0.05,8.);
	histo2DRAAAll3Up->GetYaxis()->SetLabelOffset(0.005);
	histo2DRAAAll3Up->GetXaxis()->SetRangeUser(0.,14.2);
	histo2DRAAAll3Up->GetYaxis()->SetRangeUser(0.00,2.1);
	histo2DRAAAll3Up->DrawCopy("");

	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd0010, markerStyleCommmonSpectrum0010,markerSizeCommonSpectrum0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE, kGray+1);
	graphRAASysCombInd0010->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd0010, markerStyleCommmonSpectrum0010,markerSizeCommonSpectrum0010, colorComb0010 , colorComb0010);
	graphRAACombInd0010->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAA_0010, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAA_0010, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
	DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAA_0010, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphWA98_17_3GeVRAA_0013, markerStyleWA98,markerSizeWA98, kGray+2 , kGray+2);

	graphPHENIX200GeVRAA_0010->Draw("p,same,e1");	
	graphPHENIX39GeVRAA_0010->Draw("p,same,e1");	
	graphPHENIX62GeVRAA_0010->Draw("p,same,e1");	
	graphWA98_17_3GeVRAA_0013->Draw("p,same,e1");	


	histo2DRAAAll3Up->Draw("axis,same");

	TLatex *labelRAAALICEPbPb0010 = new TLatex(0.35,0.93,"#pi^{0} ALICE    0-10% Pb-Pb");
	SetStyleTLatex( labelRAAALICEPbPb0010, 0.85*textsizeLabelsRAA,4);
	labelRAAALICEPbPb0010->Draw();
	TLegend* legendRAASinglePbPb0010 = new TLegend(0.35,0.88,0.65,0.92);
	legendRAASinglePbPb0010->SetFillColor(0);
	legendRAASinglePbPb0010->SetLineColor(0);
	legendRAASinglePbPb0010->SetNColumns(1);
	legendRAASinglePbPb0010->SetTextFont(42);
	legendRAASinglePbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAASinglePbPb0010->AddEntry(graphRAASysCombInd0010,"#sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
	legendRAASinglePbPb0010->Draw();

	TLatex *labelRAAPHENIXPbPb0010 = new TLatex(0.35,0.83,"#pi^{0} PHENIX 0-10% Au-Au");
	SetStyleTLatex( labelRAAPHENIXPbPb0010, 0.85*textsizeLabelsRAA,4);
	labelRAAPHENIXPbPb0010->Draw();

	TLegend* legendRAARHICPbPb0010 = new TLegend(0.35,0.73,0.95,0.82);
	legendRAARHICPbPb0010->SetFillColor(0);
	legendRAARHICPbPb0010->SetLineColor(0);
	legendRAARHICPbPb0010->SetNColumns(2);
	legendRAARHICPbPb0010->SetTextFont(42);
	legendRAARHICPbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAARHICPbPb0010->AddEntry(graphPHENIX200GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
	legendRAARHICPbPb0010->AddEntry(graphPHENIX62GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
	legendRAARHICPbPb0010->AddEntry(graphPHENIX39GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
	legendRAARHICPbPb0010->Draw();

	TLatex *labelRAAWA98PbPb0010 = new TLatex(0.35,0.68,"#pi^{0} WA98     0-13% Pb-Pb");
	SetStyleTLatex( labelRAAWA98PbPb0010, 0.85*textsizeLabelsRAA,4);
	labelRAAWA98PbPb0010->Draw();

	TLegend* legendRAASPSPbPb0010 = new TLegend(0.35,0.63,0.95,0.67);
	legendRAASPSPbPb0010->SetFillColor(0);
	legendRAASPSPbPb0010->SetLineColor(0);
	legendRAASPSPbPb0010->SetNColumns(2);
	legendRAASPSPbPb0010->SetTextFont(42);
	legendRAASPSPbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAASPSPbPb0010->AddEntry(graphWA98_17_3GeVRAA_0013,"#sqrt{#it{s}_{_{NN}}} = 17.3 GeV","p");
	legendRAASPSPbPb0010->Draw();

	TBox* boxErrorNorm0010_Single = CreateBoxConv(colorComb0005Box, 0.2, 1.-normErr0010 , 0.5, 1.+normErr0010);
	boxErrorNorm0010_Single->Draw();


	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);


	canvasRAA_0010->Update();
	canvasRAA_0010->Print(Form("%s/RAA_0010_Paper.%s",outputDir.Data(),suffix.Data()));



		TCanvas* canvasRAA_6080 = new TCanvas("canvasRAA_6080","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasRAA_6080,  0.11, 0.01, 0.015, 0.11);

	histo2DRAAAll3Up->DrawCopy("");

	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080, widthLinesBoxes, kTRUE, kGray+1);
	graphRAASysCombInd6080->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
	graphRAACombInd6080->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAA_6080, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAA_6080, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
	DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAA_6080, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);

	graphPHENIX200GeVRAA_6080->Draw("p,same,e1");   
	graphPHENIX39GeVRAA_6080->Draw("p,same,e1"); 
	graphPHENIX62GeVRAA_6080->Draw("p,same,e1"); 


	histo2DRAAAll3Up->Draw("axis,same");

	TLatex *labelRAAALICEPbPb6080 = new TLatex(0.4,0.93,"#pi^{0} ALICE   60-80% Pb-Pb");
	SetStyleTLatex( labelRAAALICEPbPb6080, 0.85*textsizeLabelsRAA,4);
	labelRAAALICEPbPb6080->Draw();
	TLegend* legendRAASinglePbPb6080 = new TLegend(0.4,0.88,0.65,0.92);
	legendRAASinglePbPb6080->SetFillColor(0);
	legendRAASinglePbPb6080->SetLineColor(0);
	legendRAASinglePbPb6080->SetNColumns(1);
	legendRAASinglePbPb6080->SetTextFont(42);
	legendRAASinglePbPb6080->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAASinglePbPb6080->AddEntry(graphRAASysCombInd6080,"#sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
	legendRAASinglePbPb6080->Draw();

	TLatex *labelRAAPHENIXPbPb6080 = new TLatex(0.4,0.82,"#pi^{0} PHENIX 60-80% Au-Au");
	SetStyleTLatex( labelRAAPHENIXPbPb6080, 0.85*textsizeLabelsRAA,4);
	labelRAAPHENIXPbPb6080->Draw();

	TLegend* legendRAARHICPbPb6080 = new TLegend(0.4,0.72,0.95,0.81);
	legendRAARHICPbPb6080->SetFillColor(0);
	legendRAARHICPbPb6080->SetLineColor(0);
	legendRAARHICPbPb6080->SetNColumns(2);
	legendRAARHICPbPb6080->SetTextFont(42);
	legendRAARHICPbPb6080->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAARHICPbPb6080->AddEntry(graphPHENIX200GeVRAA_6080,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
	legendRAARHICPbPb6080->AddEntry(graphPHENIX62GeVRAA_6080,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
	legendRAARHICPbPb6080->AddEntry(graphPHENIX39GeVRAA_6080,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
	legendRAARHICPbPb6080->Draw();

	TBox* boxErrorNorm6080_Single = CreateBoxConv(colorComb6080Box, 0.2, 1.-normErr6080 , 0.5, 1.+normErr6080);
	boxErrorNorm6080_Single->Draw();


	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);


	canvasRAA_6080->Update();
	canvasRAA_6080->Print(Form("%s/RAA_6080_Paper.%s",outputDir.Data(),suffix.Data()));


	//    cout << "here ............................." << endl;
	// 
	// 	Double_t DMesonRAAPt_0020[7] = {2.5, 3.5, 4.5, 5.5,  7., 10., 14.};
	// 	Double_t DMesonRAAPtErr_0020[7] = {0.5, 0.5, 0.5, 0.5,  1., 2., 2.};
	// 	Double_t DMesonRAA_0020[7] =  {0.51, 0.37, 0.33, 0.27, 0.28, 0.26, 0.35};
	// 	Double_t DMesonRAAStatErr_0020[7] = {0.10,0.06,0.05,0.07,0.04,0.03,0.06};
	// 	Double_t DMesonRAASysErrUp_0020[7] = {0.18,0.11,0.10,0.08,0.07,0.06,0.10};
	// 	Double_t DMesonRAASysErrDown_0020[7] = {0.22,0.13,0.11,0.09,0.08,0.07,0.12};
	// 
	// 	TGraphAsymmErrors* graphDMesonRAA_0020 = new TGraphAsymmErrors(7);
	// 	TGraphAsymmErrors* graphDMesonRAASys_0020 = new TGraphAsymmErrors(7);
	// 	for(Int_t i=0; i<7; i++){
	// // 		Double_t quadErr = TMath::Sqrt(DMesonRAAStatErr_0020[i]*DMesonRAAStatErr_0020[i] + DMesonRAAErrPt_0020[i] * DMesonRAAErrPt_0020[i]);
	// 		graphDMesonRAA_0020->SetPoint(i, DMesonRAAPt_0020[i],  DMesonRAA_0020[i]);
	// 		graphDMesonRAA_0020->SetPointError(i, 0. ,0., DMesonRAAStatErr_0020[i],DMesonRAAStatErr_0020[i]);
	// 		graphDMesonRAASys_0020->SetPoint(i, DMesonRAAPt_0020[i],  DMesonRAA_0020[i]);
	// 		graphDMesonRAASys_0020->SetPointError(i, DMesonRAAPtErr_0020[i] ,DMesonRAAPtErr_0020[i], DMesonRAASysErrDown_0020[i],DMesonRAASysErrUp_0020[i]);
	// 	}
	// 
	// 	
	TCanvas* canvasRAA_0005 = new TCanvas("canvasRAA_0005","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasRAA_0005,  0.11, 0.01, 0.015, 0.11);

	// 	canvasRAA_0005->SetLogy();
	// 	canvasRAA_0005->SetLogx();
	SetStyleHistoTH2ForGraphs(histo2DRAAAll3Up, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.8,0.8, 512, 505); //#frac{#frac{1
	// 	histo2DRAAAll3Up->GetYaxis()->SetRangeUser(0.05,8.);
	histo2DRAAAll3Up->GetYaxis()->SetRangeUser(0.05,1.45);
	histo2DRAAAll3Up->DrawCopy("");

	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005, widthLinesBoxes, kTRUE, kGray+1);
	graphRAASysCombInd0005->Draw("E2same");
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005);
	graphRAACombInd0005->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphErr(graphChargedRAAStat_0005 ,26,markerSizeChargedHadronSpectrum, kGray+2 , kGray+2);					 
	DrawGammaSetMarkerTGraphErr(graphChargedRAASys_0005 ,26,markerSizeChargedHadronSpectrum, kGray+2 , kGray+2);					 

	graphChargedRAAStat_0005->Draw("p,same,e1");
	graphChargedRAASys_0005->Draw("E[]same");
	DrawGammaSetMarker(histoChargedPionRAAStat0005 ,24,markerSizeChargedHadronSpectrum, colorComb0005-2 , colorComb0005-2);					 
	TGraphAsymmErrors* graphChargedPionRAASyst0005 = new TGraphAsymmErrors(histoChargedPionRAASyst0005);   
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASyst0005, 25,markerSizeChargedHadronSpectrum, colorComb0005-2 , colorComb0005-2, widthLinesBoxes, kTRUE);//, colorComb0005-5);
	graphChargedPionRAASyst0005->Draw("E2same");
	histoChargedPionRAAStat0005->Draw("p,same,e1");


	DrawGammaSetMarkerTGraphAsym(graphRAASysCombInd6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080, widthLinesBoxes, kTRUE, kGray+1);
	DrawGammaSetMarkerTGraphAsym(graphRAACombInd6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
	DrawGammaSetMarkerTGraphErr(graphChargedRAAStat_6080 ,27,markerSizeChargedHadronSpectrum, kGray+2 , kGray+2);					 
	DrawGammaSetMarkerTGraphErr(graphChargedRAASys_6080 ,27,markerSizeChargedHadronSpectrum, kGray+2 , kGray+2);					 


	TGraphAsymmErrors* graphChargedPionRAASyst6080 = new TGraphAsymmErrors(histoChargedPionRAASyst6080);   
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASyst6080, 25,markerSizeChargedHadronSpectrum, colorComb6080-2 , colorComb6080-2, widthLinesBoxes, kTRUE);//, colorComb6080-5);
	DrawGammaSetMarker(histoChargedPionRAAStat6080 ,25,markerSizeChargedHadronSpectrum, colorComb6080-2 , colorComb6080-2);					 


	graphRAASysCombInd6080->Draw("E2same");
	graphChargedPionRAASyst6080->Draw("E2same");
	graphRAACombInd6080->Draw("p,same,e1");

	histoChargedPionRAAStat6080->Draw("p,same,e1");

	graphChargedRAAStat_6080->Draw("p,same,e1");
	graphChargedRAASys_6080->Draw("E[]same");
	// 	DrawGammaSetMarkerTGraphAsym(graphDMesonRAA_0020, 21,2., kViolet+2 , kViolet+2);
	// 	DrawGammaSetMarkerTGraphAsym(graphDMesonRAASys_0020,  21,2., kViolet+2 , kViolet+2,widthLinesBoxes, kTRUE,kWhite);
	// 	graphDMesonRAASys_0020->Draw("E2same");
	// 	graphDMesonRAA_0020->Draw("p,same,e1");

	histo2DRAAAll3Up->Draw("axis,same");

	// 	TLatex *labelRAAALICEPbPb0005 = new TLatex(0.4,0.93,"");
	// 	SetStyleTLatex( labelRAAALICEPbPb0005, 0.025,4);
	// 	labelRAAALICEPbPb0005->Draw();
	TLegend* legendRAASinglePbPb0005 = new TLegend(0.61,0.825,0.96,0.96);
	legendRAASinglePbPb0005->SetFillColor(0);
	legendRAASinglePbPb0005->SetLineColor(0);
	legendRAASinglePbPb0005->SetNColumns(2);
	legendRAASinglePbPb0005->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAASinglePbPb0005->AddEntry(graphRAASysCombInd0005,"0-5% #pi^{0}","pf");
	legendRAASinglePbPb0005->AddEntry(graphRAASysCombInd6080,"60-80% #pi^{0}","pf");
	legendRAASinglePbPb0005->AddEntry(graphChargedPionRAASyst0005,"0-5% #pi^{#pm}","pf");
	legendRAASinglePbPb0005->AddEntry(graphChargedPionRAASyst6080,"60-80% #pi^{#pm}","pf");
	legendRAASinglePbPb0005->AddEntry(graphChargedRAASys_0005,"0-5% h^{#pm}","pe");
	legendRAASinglePbPb0005->AddEntry(graphChargedRAASys_6080,"60-80% h^{#pm}","pe");

	legendRAASinglePbPb0005->Draw();


	// 	TBox* boxErrorNorm0005_Single = CreateBoxConv(colorComb0005Box, 0.2, 1.-normErr0005 , 0.5, 1.+normErr0005);
	// 	boxErrorNorm0005_Single->Draw();

	TLatex *labelRAAPi0Comb = new TLatex(0.15,0.9,"Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV");
	SetStyleTLatex( labelRAAPi0Comb, 0.85*textsizeLabelsRAA,4);
	labelRAAPi0Comb->Draw();


	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);


	canvasRAA_0005->Update();
	canvasRAA_0005->Print(Form("%s/RAA_0005_DiffPart.%s",outputDir.Data(),suffix.Data()));

	//***************************************************************************************************
	//***************************** Creating new Delta Pt/ Pt plot **************************************
	//***************************************************************************************************

	//************************************ 0-5% *********************************************************   
	cout << "*********************** 0-5% **********************" << endl;
	TGraphAsymmErrors* graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0005 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombStatPi02760GeVPlot->Clone());
	TGraphAsymmErrors* graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0005 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombSysPi02760GeVPlot->Clone());

	for (Int_t i=0; i < graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0005->GetN(); i++) {
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0005->GetY()[i] = nColl0005 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0005->GetY()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0005->GetEYhigh()[i] = nColl0005 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0005->GetEYhigh()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0005->GetEYlow()[i] = nColl0005 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0005->GetEYlow()[i];

		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0005->GetY()[i] = nColl0005 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0005->GetY()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0005->GetEYhigh()[i] = nColl0005 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0005->GetEYhigh()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0005->GetEYlow()[i] = nColl0005 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0005->GetEYlow()[i];
	}

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_StatErr0005 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0005StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_SysErr0005 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0005FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0005, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0005,
			graphYieldPi0CombPbPb0005StatErrUnscaled, graphYieldPi0CombPbPb0005FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPP_StatErr0005, graphDeltaPtOverPtComb_PtPP_SysErr0005);

	graphDeltaPtOverPtComb_PtPP_StatErr0005->RemovePoint(graphDeltaPtOverPtComb_PtPP_StatErr0005->GetN()-1);
	graphDeltaPtOverPtComb_PtPP_SysErr0005->RemovePoint(graphDeltaPtOverPtComb_PtPP_SysErr0005->GetN()-1);


	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_StatErr0005 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0005StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_SysErr0005 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0005FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0005, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0005,
			graphYieldPi0CombPbPb0005StatErrUnscaled, graphYieldPi0CombPbPb0005FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPbPb_StatErr0005, graphDeltaPtOverPtComb_PtPbPb_SysErr0005,0);

	graphDeltaPtOverPtComb_PtPbPb_StatErr0005->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_StatErr0005->GetN()-1);
	graphDeltaPtOverPtComb_PtPbPb_SysErr0005->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_SysErr0005->GetN()-1);

	//************************************ 5-10% ********************************************************   
	cout << "*********************** 5-10% **********************" << endl;
	TGraphAsymmErrors* graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0510 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombStatPi02760GeVPlot->Clone());
	TGraphAsymmErrors* graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0510 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombSysPi02760GeVPlot->Clone());

	for (Int_t i=0; i < graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0510->GetN(); i++) {
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0510->GetY()[i] = nColl0510 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0510->GetY()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0510->GetEYhigh()[i] = nColl0510 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0510->GetEYhigh()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0510->GetEYlow()[i] = nColl0510 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0510->GetEYlow()[i];

		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0510->GetY()[i] = nColl0510 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0510->GetY()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0510->GetEYhigh()[i] = nColl0510 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0510->GetEYhigh()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0510->GetEYlow()[i] = nColl0510 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0510->GetEYlow()[i];
	}

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_StatErr0510 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0510StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_SysErr0510 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0510FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0510, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0510,
			graphYieldPi0CombPbPb0510StatErrUnscaled, graphYieldPi0CombPbPb0510FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPP_StatErr0510, graphDeltaPtOverPtComb_PtPP_SysErr0510);

	graphDeltaPtOverPtComb_PtPP_StatErr0510->RemovePoint(graphDeltaPtOverPtComb_PtPP_StatErr0510->GetN()-1);
	graphDeltaPtOverPtComb_PtPP_SysErr0510->RemovePoint(graphDeltaPtOverPtComb_PtPP_SysErr0510->GetN()-1);

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_StatErr0510 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0510StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_SysErr0510 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0510FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0510, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0510,
			graphYieldPi0CombPbPb0510StatErrUnscaled, graphYieldPi0CombPbPb0510FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPbPb_StatErr0510, graphDeltaPtOverPtComb_PtPbPb_SysErr0510,0);

	graphDeltaPtOverPtComb_PtPbPb_StatErr0510->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_StatErr0510->GetN()-1);
	graphDeltaPtOverPtComb_PtPbPb_SysErr0510->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_SysErr0510->GetN()-1);

	//************************************ 0-10% ********************************************************   
	cout << "*********************** 0-10% **********************" << endl;
	TGraphAsymmErrors* graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0010 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombStatPi02760GeVPlot->Clone());
	TGraphAsymmErrors* graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0010 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombSysPi02760GeVPlot->Clone());

	for (Int_t i=0; i < graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0010->GetN(); i++) {
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0010->GetY()[i] = nColl0010 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0010->GetY()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0010->GetEYhigh()[i] = nColl0010 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0010->GetEYhigh()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0010->GetEYlow()[i] = nColl0010 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0010->GetEYlow()[i];

		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0010->GetY()[i] = nColl0010 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0010->GetY()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0010->GetEYhigh()[i] = nColl0010 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0010->GetEYhigh()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0010->GetEYlow()[i] = nColl0010 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0010->GetEYlow()[i];
	}

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_StatErr0010 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0010StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_SysErr0010 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0010FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0010, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0010,
			graphYieldPi0CombPbPb0010StatErrUnscaled, graphYieldPi0CombPbPb0010FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPP_StatErr0010, graphDeltaPtOverPtComb_PtPP_SysErr0010);

	graphDeltaPtOverPtComb_PtPP_StatErr0010->RemovePoint(graphDeltaPtOverPtComb_PtPP_StatErr0010->GetN()-1);
	graphDeltaPtOverPtComb_PtPP_SysErr0010->RemovePoint(graphDeltaPtOverPtComb_PtPP_SysErr0010->GetN()-1);

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_StatErr0010 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0010StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_SysErr0010 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb0010FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled0010, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled0010,
			graphYieldPi0CombPbPb0010StatErrUnscaled, graphYieldPi0CombPbPb0010FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPbPb_StatErr0010, graphDeltaPtOverPtComb_PtPbPb_SysErr0010,0);

	graphDeltaPtOverPtComb_PtPbPb_StatErr0010->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_StatErr0010->GetN()-1);
	graphDeltaPtOverPtComb_PtPbPb_SysErr0010->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_SysErr0010->GetN()-1);


	//************************************ 10-20% *******************************************************  
	cout << "*********************** 10-20% **********************" << endl;
	TGraphAsymmErrors* graphInvYieldComb_StatErr_PP2760GeV_NCollScaled1020 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombStatPi02760GeVPlot->Clone());
	TGraphAsymmErrors* graphInvYieldComb_SysErr_PP2760GeV_NCollScaled1020 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombSysPi02760GeVPlot->Clone());

	for (Int_t i=0; i < graphInvYieldComb_StatErr_PP2760GeV_NCollScaled1020->GetN(); i++) {
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled1020->GetY()[i] = nColl1020 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled1020->GetY()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled1020->GetEYhigh()[i] = nColl1020 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled1020->GetEYhigh()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled1020->GetEYlow()[i] = nColl1020 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled1020->GetEYlow()[i];

		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled1020->GetY()[i] = nColl1020 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled1020->GetY()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled1020->GetEYhigh()[i] = nColl1020 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled1020->GetEYhigh()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled1020->GetEYlow()[i] = nColl1020 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled1020->GetEYlow()[i];
	}

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_StatErr1020 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb1020StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_SysErr1020 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb1020FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled1020, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled1020,
			graphYieldPi0CombPbPb1020StatErrUnscaled, graphYieldPi0CombPbPb1020FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPP_StatErr1020, graphDeltaPtOverPtComb_PtPP_SysErr1020);

	graphDeltaPtOverPtComb_PtPP_StatErr1020->RemovePoint(graphDeltaPtOverPtComb_PtPP_StatErr1020->GetN()-1);
	graphDeltaPtOverPtComb_PtPP_SysErr1020->RemovePoint(graphDeltaPtOverPtComb_PtPP_SysErr1020->GetN()-1);

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_StatErr1020 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb1020StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_SysErr1020 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb1020FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled1020, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled1020,
			graphYieldPi0CombPbPb1020StatErrUnscaled, graphYieldPi0CombPbPb1020FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPbPb_StatErr1020, graphDeltaPtOverPtComb_PtPbPb_SysErr1020,0);

	graphDeltaPtOverPtComb_PtPbPb_StatErr1020->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_StatErr1020->GetN()-1);
	graphDeltaPtOverPtComb_PtPbPb_SysErr1020->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_SysErr1020->GetN()-1);

	//************************************ 20-40% *******************************************************   
	cout << "*********************** 20-40% **********************" << endl;
	TGraphAsymmErrors* graphInvYieldComb_StatErr_PP2760GeV_NCollScaled2040 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombStatPi02760GeVPlot->Clone());
	TGraphAsymmErrors* graphInvYieldComb_SysErr_PP2760GeV_NCollScaled2040 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombSysPi02760GeVPlot->Clone());

	for (Int_t i=0; i < graphInvYieldComb_StatErr_PP2760GeV_NCollScaled2040->GetN(); i++) {
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled2040->GetY()[i] = nColl2040 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled2040->GetY()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled2040->GetEYhigh()[i] = nColl2040 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled2040->GetEYhigh()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled2040->GetEYlow()[i] = nColl2040 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled2040->GetEYlow()[i];

		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled2040->GetY()[i] = nColl2040 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled2040->GetY()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled2040->GetEYhigh()[i] = nColl2040 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled2040->GetEYhigh()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled2040->GetEYlow()[i] = nColl2040 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled2040->GetEYlow()[i];
	}

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_StatErr2040 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb2040StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_SysErr2040 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb2040FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled2040, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled2040,
			graphYieldPi0CombPbPb2040StatErrUnscaled, graphYieldPi0CombPbPb2040FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPP_StatErr2040, graphDeltaPtOverPtComb_PtPP_SysErr2040);

	graphDeltaPtOverPtComb_PtPP_StatErr2040->RemovePoint(graphDeltaPtOverPtComb_PtPP_StatErr2040->GetN()-1);
	graphDeltaPtOverPtComb_PtPP_SysErr2040->RemovePoint(graphDeltaPtOverPtComb_PtPP_SysErr2040->GetN()-1);

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_StatErr2040 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb2040StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_SysErr2040 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb2040FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled2040, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled2040,
			graphYieldPi0CombPbPb2040StatErrUnscaled, graphYieldPi0CombPbPb2040FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPbPb_StatErr2040, graphDeltaPtOverPtComb_PtPbPb_SysErr2040,0);

	graphDeltaPtOverPtComb_PtPbPb_StatErr2040->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_StatErr2040->GetN()-1);
	graphDeltaPtOverPtComb_PtPbPb_SysErr2040->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_SysErr2040->GetN()-1);


	//************************************ 40-60% *******************************************************   
	cout << "*********************** 40-60% **********************" << endl;
	TGraphAsymmErrors* graphInvYieldComb_StatErr_PP2760GeV_NCollScaled4060 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombStatPi02760GeVPlot->Clone());
	TGraphAsymmErrors* graphInvYieldComb_SysErr_PP2760GeV_NCollScaled4060 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombSysPi02760GeVPlot->Clone());

	for (Int_t i=0; i < graphInvYieldComb_StatErr_PP2760GeV_NCollScaled4060->GetN(); i++) {
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled4060->GetY()[i] = nColl4060 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled4060->GetY()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled4060->GetEYhigh()[i] = nColl4060 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled4060->GetEYhigh()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled4060->GetEYlow()[i] = nColl4060 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled4060->GetEYlow()[i];

		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled4060->GetY()[i] = nColl4060 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled4060->GetY()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled4060->GetEYhigh()[i] = nColl4060 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled4060->GetEYhigh()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled4060->GetEYlow()[i] = nColl4060 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled4060->GetEYlow()[i];
	}

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_StatErr4060 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb4060StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_SysErr4060 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb4060FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled4060, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled4060,
			graphYieldPi0CombPbPb4060StatErrUnscaled, graphYieldPi0CombPbPb4060FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPP_StatErr4060, graphDeltaPtOverPtComb_PtPP_SysErr4060);

	graphDeltaPtOverPtComb_PtPP_StatErr4060->RemovePoint(graphDeltaPtOverPtComb_PtPP_StatErr4060->GetN()-1);
	graphDeltaPtOverPtComb_PtPP_SysErr4060->RemovePoint(graphDeltaPtOverPtComb_PtPP_SysErr4060->GetN()-1);

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_StatErr4060 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb4060StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_SysErr4060 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb4060FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled4060, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled4060,
			graphYieldPi0CombPbPb4060StatErrUnscaled, graphYieldPi0CombPbPb4060FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPbPb_StatErr4060, graphDeltaPtOverPtComb_PtPbPb_SysErr4060,0);

	graphDeltaPtOverPtComb_PtPbPb_StatErr4060->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_StatErr4060->GetN()-1);
	graphDeltaPtOverPtComb_PtPbPb_SysErr4060->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_SysErr4060->GetN()-1);


	//************************************ 60-80% *******************************************************   
	cout << "*********************** 60-80% **********************" << endl;
	TGraphAsymmErrors* graphInvYieldComb_StatErr_PP2760GeV_NCollScaled6080 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombStatPi02760GeVPlot->Clone());
	TGraphAsymmErrors* graphInvYieldComb_SysErr_PP2760GeV_NCollScaled6080 = dynamic_cast<TGraphAsymmErrors*>(graphInvSectionCombSysPi02760GeVPlot->Clone());

	for (Int_t i=0; i < graphInvYieldComb_StatErr_PP2760GeV_NCollScaled6080->GetN(); i++) {
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled6080->GetY()[i] = nColl6080 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled6080->GetY()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled6080->GetEYhigh()[i] = nColl6080 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled6080->GetEYhigh()[i];
		graphInvYieldComb_StatErr_PP2760GeV_NCollScaled6080->GetEYlow()[i] = nColl6080 * graphInvYieldComb_StatErr_PP2760GeV_NCollScaled6080->GetEYlow()[i];

		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled6080->GetY()[i] = nColl6080 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled6080->GetY()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled6080->GetEYhigh()[i] = nColl6080 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled6080->GetEYhigh()[i];
		graphInvYieldComb_SysErr_PP2760GeV_NCollScaled6080->GetEYlow()[i] = nColl6080 * graphInvYieldComb_SysErr_PP2760GeV_NCollScaled6080->GetEYlow()[i];
	}

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_StatErr6080 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb6080StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPP_SysErr6080 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb6080FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled6080, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled6080,
			graphYieldPi0CombPbPb6080StatErrUnscaled, graphYieldPi0CombPbPb6080FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPP_StatErr6080, graphDeltaPtOverPtComb_PtPP_SysErr6080);

	graphDeltaPtOverPtComb_PtPP_StatErr6080->RemovePoint(graphDeltaPtOverPtComb_PtPP_StatErr6080->GetN()-1);
	graphDeltaPtOverPtComb_PtPP_SysErr6080->RemovePoint(graphDeltaPtOverPtComb_PtPP_SysErr6080->GetN()-1);

	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_StatErr6080 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb6080StatErrUnscaled->Clone());
	TGraphAsymmErrors* graphDeltaPtOverPtComb_PtPbPb_SysErr6080 = dynamic_cast<TGraphAsymmErrors*>(graphYieldPi0CombPbPb6080FullSysErrUnscaled->Clone());

	CalcDeltaPtOverPt(graphInvYieldComb_StatErr_PP2760GeV_NCollScaled6080, graphInvYieldComb_SysErr_PP2760GeV_NCollScaled6080,
			graphYieldPi0CombPbPb6080StatErrUnscaled, graphYieldPi0CombPbPb6080FullSysErrUnscaled, 
			graphDeltaPtOverPtComb_PtPbPb_StatErr6080, graphDeltaPtOverPtComb_PtPbPb_SysErr6080,0);

	graphDeltaPtOverPtComb_PtPbPb_StatErr6080->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_StatErr6080->GetN()-1);
	graphDeltaPtOverPtComb_PtPbPb_SysErr6080->RemovePoint(graphDeltaPtOverPtComb_PtPbPb_SysErr6080->GetN()-1);

	TCanvas* canvasDeltaPtOverPt_PtPP = new TCanvas("canvasDeltaPtOverPt_PtPP","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasDeltaPtOverPt_PtPP, 0.1, 0.02, 0.02, 0.09);
	canvasDeltaPtOverPt_PtPP->cd();

	TH2F* histo2DDeltaPtOverPtAll_PtPP = new TH2F("histo2DDeltaPtOverPtAll_PtPP", "", 1, 0., 18., 1, -0.05, 0.7); 
	SetStyleHistoTH2ForGraphs(histo2DDeltaPtOverPtAll_PtPP, "#it{p}_{T,pp} (GeV/#it{c})","#it{S}_{loss}( #equiv #delta #it{p}_{T}/#it{p}_{T,pp})", 0.035,0.045, 0.035,0.045, 0.82,1., 512, 510); 
	histo2DDeltaPtOverPtAll_PtPP->DrawCopy();

			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_SysErr0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005, widthLinesBoxes, kTRUE);//, colorComb0005-5);
			graphDeltaPtOverPtComb_PtPP_SysErr0005->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_StatErr0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005);
			graphDeltaPtOverPtComb_PtPP_StatErr0005->Draw("p,same,e1");
			
	/*        DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_SysErr0510, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510, widthLinesBoxes, kTRUE);//, colorComb0510-5);
			graphDeltaPtOverPtComb_PtPP_SysErr0510->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_StatErr0510, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510);
			graphDeltaPtOverPtComb_PtPP_StatErr0510->Draw("p,same,e1");
	*/        
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_SysErr1020, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020, widthLinesBoxes, kTRUE);//, colorComb1020-5);
			graphDeltaPtOverPtComb_PtPP_SysErr1020->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_StatErr1020, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020);
			graphDeltaPtOverPtComb_PtPP_StatErr1020->Draw("p,same,e1");
			
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_SysErr2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);//, colorComb2040-5);
			graphDeltaPtOverPtComb_PtPP_SysErr2040->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_StatErr2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
			graphDeltaPtOverPtComb_PtPP_StatErr2040->Draw("p,same,e1");
			
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_SysErr4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060, widthLinesBoxes, kTRUE);//, colorComb4060-5);
			graphDeltaPtOverPtComb_PtPP_SysErr4060->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_StatErr4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
			graphDeltaPtOverPtComb_PtPP_StatErr4060->Draw("p,same,e1");
			
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_SysErr6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080, widthLinesBoxes, kTRUE);//, colorComb6080-5);
			graphDeltaPtOverPtComb_PtPP_SysErr6080->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPP_StatErr6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
			graphDeltaPtOverPtComb_PtPP_StatErr6080->Draw("p,same,e1");
			
			TLegend* legendDeltPtOverPt_PtPP = new TLegend(0.45,0.7,0.95,0.95);
			legendDeltPtOverPt_PtPP->SetFillColor(0);
			legendDeltPtOverPt_PtPP->SetLineColor(0);
			legendDeltPtOverPt_PtPP->SetTextSize(0.035);
			legendDeltPtOverPt_PtPP->SetMargin(0.12);
			legendDeltPtOverPt_PtPP->AddEntry(graphYieldPi0CombPbPb0005FullSysErr,collisionSystemCent0005Legend.Data(),"pf");
	//          legendDeltPtOverPt_PtPP->AddEntry(graphYieldPi0CombPbPb0510SysErr,collisionSystemCent0510Legend.Data(),"pf");
			legendDeltPtOverPt_PtPP->AddEntry(graphYieldPi0CombPbPb1020FullSysErr,collisionSystemCent1020.Data(),"pf");
			legendDeltPtOverPt_PtPP->AddEntry(graphYieldPi0CombPbPb2040FullSysErr,collisionSystemCent2040.Data(),"pf");
			legendDeltPtOverPt_PtPP->AddEntry(graphYieldPi0CombPbPb4060FullSysErr,collisionSystemCent4060.Data(),"pf");
			legendDeltPtOverPt_PtPP->AddEntry(graphYieldPi0CombPbPb6080FullSysErr,collisionSystemCent6080.Data(),"pf");
			legendDeltPtOverPt_PtPP->Draw();
			
			DrawGammaLines(0., 18. , 0, 0.,1.,kGray);
			
	canvasDeltaPtOverPt_PtPP->SaveAs(Form("%s/DeltaPtOverPt_PtPP_All.%s",outputDir.Data(),suffix.Data()));


	TCanvas* canvasDeltaPtOverPt_PtPbPb = new TCanvas("canvasDeltaPtOverPt_PtPbPb","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasDeltaPtOverPt_PtPbPb, 0.1, 0.02, 0.02, 0.09);
	canvasDeltaPtOverPt_PtPbPb->cd();

	TH2F* histo2DDeltaPtOverPtAll_PtPbPb = new TH2F("histo2DDeltaPtOverPtAll_PtPbPb", "", 1, 0., 15., 1, -0.05, 0.7); 
	SetStyleHistoTH2ForGraphs(histo2DDeltaPtOverPtAll_PtPbPb, "#it{p}_{T,Pb-Pb} (GeV/#it{c})","#it{S}_{loss}( #equiv #delta #it{p}_{T}/#it{p}_{T,pp})", 0.035,0.045, 0.035,0.045, 0.82,1., 512, 510); 
	histo2DDeltaPtOverPtAll_PtPbPb->Draw();

			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_SysErr0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005, widthLinesBoxes, kTRUE);//, colorComb0005-5);
			graphDeltaPtOverPtComb_PtPbPb_SysErr0005->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_StatErr0005, markerStyleCommmonSpectrum0005,markerSizeCommonSpectrum0005, colorComb0005 , colorComb0005);
			graphDeltaPtOverPtComb_PtPbPb_StatErr0005->Draw("p,same,e1");
			
	/*        DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_SysErr0510, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510, widthLinesBoxes, kTRUE);//, colorComb0510-5);
			graphDeltaPtOverPtComb_PtPbPb_SysErr0510->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_StatErr0510, markerStyleCommmonSpectrum0510,markerSizeCommonSpectrum0510, colorComb0510 , colorComb0510);
			graphDeltaPtOverPtComb_PtPbPb_StatErr0510->Draw("p,same,e1");
	*/        
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_SysErr1020, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020, widthLinesBoxes, kTRUE);//, colorComb1020-5);
			graphDeltaPtOverPtComb_PtPbPb_SysErr1020->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_StatErr1020, markerStyleCommmonSpectrum1020,markerSizeCommonSpectrum1020, colorComb1020 , colorComb1020);
			graphDeltaPtOverPtComb_PtPbPb_StatErr1020->Draw("p,same,e1");
			
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_SysErr2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE);//, colorComb2040-5);
			graphDeltaPtOverPtComb_PtPbPb_SysErr2040->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_StatErr2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
			graphDeltaPtOverPtComb_PtPbPb_StatErr2040->Draw("p,same,e1");
			
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_SysErr4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060, widthLinesBoxes, kTRUE);//, colorComb4060-5);
			graphDeltaPtOverPtComb_PtPbPb_SysErr4060->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_StatErr4060, markerStyleCommmonSpectrum4060,markerSizeCommonSpectrum4060, colorComb4060 , colorComb4060);
			graphDeltaPtOverPtComb_PtPbPb_StatErr4060->Draw("p,same,e1");
			
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_SysErr6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080, widthLinesBoxes, kTRUE);//, colorComb6080-5);
			graphDeltaPtOverPtComb_PtPbPb_SysErr6080->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphDeltaPtOverPtComb_PtPbPb_StatErr6080, markerStyleCommmonSpectrum6080,markerSizeCommonSpectrum6080, colorComb6080 , colorComb6080);
			graphDeltaPtOverPtComb_PtPbPb_StatErr6080->Draw("p,same,e1");
			
			legendDeltPtOverPt_PtPP->Draw();
			
			DrawGammaLines(0., 15. , 0, 0.,1.,kGray);
			
	canvasDeltaPtOverPt_PtPbPb->SaveAs(Form("%s/DeltaPtOverPt_PtPbPb_All.%s",outputDir.Data(),suffix.Data()));


	TCanvas* canvasRAAvsNPart = new TCanvas("canvasRAAvsNPart","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasRAAvsNPart, 0.1, 0.01, 0.015, 0.1);
	canvasRAAvsNPart->cd();

	TH2F* histo2DRAAVsNPart = new TH2F("histo2DRAAVsNPart", "", 100, 0., 420., 100, -0.05, 1.3); 
	SetStyleHistoTH2ForGraphs(histo2DRAAVsNPart, "#it{N}_{part}","#it{R}_{AA} (#it{p}_{T} = 7 GeV/#it{c})",0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.9,0.9, 512, 505);
	histo2DRAAVsNPart->GetXaxis()->SetLabelFont(42);
	histo2DRAAVsNPart->GetYaxis()->SetLabelFont(42);
	histo2DRAAVsNPart->Draw();

	DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAAvsNPartAt7GeVc, markerStylePHENIX200GeV,markerSizePHENIX200GeV*1.2, kGray+2 , kGray+2);
	DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAAvsNPartAt7GeVc, markerStylePHENIX62GeV,markerSizePHENIX62GeV*1.2, kBlack, kBlack);
	DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAAvsNPartAt7GeVc, markerStylePHENIX39GeV,markerSizePHENIX39GeV*1.2, kGray+1 , kGray+1);

	//    DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAAvsNPartAt7GeVc, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kGray+2,kGray+2); //kGreen+2 , kGreen+2);
	//    DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAAvsNPartAt7GeVc, markerStylePHENIX62GeV,markerSizePHENIX62GeV*2,kGray+2, kGray+2); //kCyan+2, kCyan+2);
	//    DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAAvsNPartAt7GeVc, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kGray+2, kGray+2);// kBlue+2 , kBlue+2);
	graphPHENIX200GeVRAAvsNPartAt7GeVc->Draw("p,same,e1");   
	graphPHENIX62GeVRAAvsNPartAt7GeVc->Draw("p,same,e1"); 
	graphPHENIX39GeVRAAvsNPartAt7GeVc->Draw("p,same,e1"); 

	TGraphErrors* graphALICE2760GeVStat = new TGraphErrors(0);
	TGraphErrors* graphALICE2760GeVSys = new TGraphErrors(0);

	Double_t x6080 = 0;
	Double_t y6080 = 0;
	Double_t yErrStat6080 = graphRAACombInd6080->GetErrorYhigh(15);
	graphRAACombInd6080->GetPoint(15, x6080,y6080);
	Double_t yErrSys6080Rel = graphRAASysCombInd6080->GetErrorYhigh(15)/y6080;
	Double_t yErrSys6080 = TMath::Sqrt(yErrSys6080Rel*yErrSys6080Rel+normErr6080*normErr6080)*y6080;
	graphALICE2760GeVStat->SetPoint(graphALICE2760GeVStat->GetN(), 22.52, y6080); 
	graphALICE2760GeVStat->SetPointError(graphALICE2760GeVStat->GetN()-1, 0, yErrStat6080);
	graphALICE2760GeVSys->SetPoint(graphALICE2760GeVSys->GetN(), 22.52, y6080); 
	graphALICE2760GeVSys->SetPointError(graphALICE2760GeVSys->GetN()-1, 2.5, yErrSys6080);

	Double_t x4060 = 0;
	Double_t y4060 = 0;
	Double_t yErrStat4060 = graphRAACombInd4060->GetErrorYhigh(15);
	graphRAACombInd4060->GetPoint(15, x4060,y4060);
	Double_t yErrSys4060Rel = graphRAASysCombInd4060->GetErrorYhigh(15)/y4060;
	Double_t yErrSys4060 = TMath::Sqrt(yErrSys4060Rel*yErrSys4060Rel+normErr4060*normErr4060)*y4060;
	graphALICE2760GeVStat->SetPoint(graphALICE2760GeVStat->GetN(), 68.56, y4060); 
	graphALICE2760GeVStat->SetPointError(graphALICE2760GeVStat->GetN()-1, 0, yErrStat4060);
	graphALICE2760GeVSys->SetPoint(graphALICE2760GeVSys->GetN(), 68.56, y4060); 
	graphALICE2760GeVSys->SetPointError(graphALICE2760GeVSys->GetN()-1, 2.5, yErrSys4060);

	Double_t x2040 = 0;
	Double_t y2040 = 0;
	Double_t yErrStat2040 = graphRAACombInd2040->GetErrorYhigh(15);
	graphRAACombInd2040->GetPoint(15, x2040,y2040);
	Double_t yErrSys2040Rel = graphRAASysCombInd2040->GetErrorYhigh(15)/y2040;
	Double_t yErrSys2040 = TMath::Sqrt(yErrSys2040Rel*yErrSys2040Rel+normErr2040*normErr2040)*y2040;
	graphALICE2760GeVStat->SetPoint(graphALICE2760GeVStat->GetN(), 157.2, y2040); 
	graphALICE2760GeVStat->SetPointError(graphALICE2760GeVStat->GetN()-1, 0, yErrStat2040);
	graphALICE2760GeVSys->SetPoint(graphALICE2760GeVSys->GetN(), 157.2, y2040); 
	graphALICE2760GeVSys->SetPointError(graphALICE2760GeVSys->GetN()-1, 2.5, yErrSys2040);

	Double_t x1020 = 0;
	Double_t y1020 = 0;
	Double_t yErrStat1020 = graphRAACombInd1020->GetErrorYhigh(15);
	graphRAACombInd1020->GetPoint(15, x1020,y1020);
	Double_t yErrSys1020Rel = graphRAASysCombInd1020->GetErrorYhigh(15)/y1020;
	Double_t yErrSys1020 = TMath::Sqrt(yErrSys1020Rel*yErrSys1020Rel+normErr1020*normErr1020)*y1020;
	graphALICE2760GeVStat->SetPoint(graphALICE2760GeVStat->GetN(), 260.1, y1020); 
	graphALICE2760GeVStat->SetPointError(graphALICE2760GeVStat->GetN()-1, 0, yErrStat1020);
	graphALICE2760GeVSys->SetPoint(graphALICE2760GeVSys->GetN(), 260.1, y1020); 
	graphALICE2760GeVSys->SetPointError(graphALICE2760GeVSys->GetN()-1, 2.5, yErrSys1020);

	Double_t x0510 = 0;
	Double_t y0510 = 0;
	Double_t yErrStat0510 = graphRAACombInd0510->GetErrorYhigh(15);
	graphRAACombInd0510->GetPoint(15, x0510,y0510);
	Double_t yErrSys0510Rel = graphRAASysCombInd0510->GetErrorYhigh(15)/y0510;
	Double_t yErrSys0510 = TMath::Sqrt(yErrSys0510Rel*yErrSys0510Rel+normErr0510*normErr0510)*y0510;
	graphALICE2760GeVStat->SetPoint(graphALICE2760GeVStat->GetN(), 329.4, y0510); 
	graphALICE2760GeVStat->SetPointError(graphALICE2760GeVStat->GetN()-1, 0., yErrStat0510);
	graphALICE2760GeVSys->SetPoint(graphALICE2760GeVSys->GetN(), 329.4, y0510); 
	graphALICE2760GeVSys->SetPointError(graphALICE2760GeVSys->GetN()-1, 2.5, yErrSys0510);

	Double_t x0005 = 0;
	Double_t y0005 = 0;
	Double_t yErrStat0005 = graphRAACombInd0005->GetErrorYhigh(15);
	graphRAACombInd0005->GetPoint(15, x0005,y0005);
	Double_t yErrSys0005Rel = graphRAASysCombInd0005->GetErrorYhigh(15)/y0005;
	Double_t yErrSys0005 = TMath::Sqrt(yErrSys0005Rel*yErrSys0005Rel+normErr0005*normErr0005)*y0005;
	graphALICE2760GeVStat->SetPoint(graphALICE2760GeVStat->GetN(), 382.7, y0005); 
	graphALICE2760GeVStat->SetPointError(graphALICE2760GeVStat->GetN()-1, 0., yErrStat0005);
	graphALICE2760GeVSys->SetPoint(graphALICE2760GeVSys->GetN(), 382.7, y0005); 
	graphALICE2760GeVSys->SetPointError(graphALICE2760GeVSys->GetN()-1, 2.5, yErrSys0005);

	DrawGammaSetMarkerTGraphErr(graphALICE2760GeVSys,  markerStyleCommmonSpectrum0010,markerSizeCommonSpectrum0010*1.3, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE, kGray+1);
	graphALICE2760GeVSys->Draw("E2same");

	DrawGammaSetMarkerTGraphErr(graphALICE2760GeVStat, markerStyleCommmonSpectrum0010,markerSizeCommonSpectrum0010*1.3, colorComb0010 , colorComb0010);
	graphALICE2760GeVStat->SetLineWidth(2.);
	graphALICE2760GeVStat->Draw("p,same,e1");
		

	TLatex *labelRAAALICEPbPbNPart = new TLatex(0.55,0.93,"#pi^{0} ALICE   Pb-Pb");
	SetStyleTLatex( labelRAAALICEPbPbNPart, 0.85*textsizeLabelsRAA,4);
	labelRAAALICEPbPbNPart->Draw();
	TLegend* legendRAASinglePbPbNPart = new TLegend(0.55,0.88,0.85,0.92);
	legendRAASinglePbPbNPart->SetFillColor(0);
	legendRAASinglePbPbNPart->SetLineColor(0);
	legendRAASinglePbPbNPart->SetNColumns(1);
	legendRAASinglePbPbNPart->SetTextFont(42);
	legendRAASinglePbPbNPart->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAASinglePbPbNPart->AddEntry(graphALICE2760GeVSys,"#sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
	legendRAASinglePbPbNPart->Draw();

	TLatex *labelRAAPHENIXPbPbNPart = new TLatex(0.55,0.82,"#pi^{0} PHENIX Au-Au");
	SetStyleTLatex( labelRAAPHENIXPbPbNPart, 0.85*textsizeLabelsRAA,4);
	labelRAAPHENIXPbPbNPart->Draw();

	TLegend* legendRAARHICPbPbNPart = new TLegend(0.55,0.66,0.85,0.8);
	legendRAARHICPbPbNPart->SetFillColor(0);
	legendRAARHICPbPbNPart->SetLineColor(0);
	legendRAARHICPbPbNPart->SetNColumns(1);
	legendRAARHICPbPbNPart->SetTextFont(42);
	legendRAARHICPbPbNPart->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAARHICPbPbNPart->AddEntry(graphPHENIX200GeVRAAvsNPartAt7GeVc,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
	legendRAARHICPbPbNPart->AddEntry(graphPHENIX62GeVRAAvsNPartAt7GeVc,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
	legendRAARHICPbPbNPart->AddEntry(graphPHENIX39GeVRAAvsNPartAt7GeVc,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
	legendRAARHICPbPbNPart->Draw();


	canvasRAAvsNPart->SaveAs(Form("%s/RAA_vs_NPart_Paper.%s",outputDir.Data(),suffix.Data()));

	TFile fCombResults(Form("CombinedResultsPbPb_%s.root",dateForOutput.Data()),"RECREATE");
		graphInvSectionCombStatPi02760GeVPlot->Write("InvYieldPPComb2760GeV_StatErr");
		graphInvSectionCombSysPi02760GeVPlot->Write("InvYieldPPComb2760GeV_SysErr");
		graphInvSectionPCMStatPi02760GeVPlot->Write("InvYieldPPPCM2760GeV_StatErr");
		graphInvSectionPCMSysPi02760GeVPlot->Write("InvYieldPPPCM2760GeV_SysErr");
		graphInvSectionPHOSStatPi02760GeVPlot->Write("InvYieldPPPHOS2760GeV_StatErr");
		graphInvSectionPHOSSysPi02760GeVPlot->Write("InvYieldPPPHOS2760GeV_SysErr");
		
		graphRAASysCombInd0005->Write("graphRAASysErr_0005");
		graphRAACombInd0005->Write("graphRAAStatErr_0005");
		graphRAASysPCM0005->Write("graphRAAPCMSysErr_0005");
		graphRAAPCM0005->Write("graphRAAPCMStatErr_0005");
		graphRAASysPHOS0005->Write("graphRAAPHOSSysErr_0005");
		graphRAAPHOS0005->Write("graphRAAPHOSStatErr_0005");
		graphYieldPi0CombPbPb0005StatErrUnscaled->Write("InvYieldPbPbStatErr_0005");
		graphYieldPi0CombPbPb0005SysErrUnscaled->Write("InvYieldPbPbSysErr_0005");
		graphYieldPi0CombPbPb0005FullSysErrUnscaled->Write("InvYieldPbPbFullSysErr_0005");
		graphPCMYieldPi0PbPb0005->Write("InvYieldPbPbPCMStatErr_0005");
		graphPCMYieldPi0SysErrPbPb0005Red->Write("InvYieldPbPbPCMSysErr_0005");
		graphPHOSYieldPi0PbPb0005->Write("InvYieldPbPbPHOSStatErr_0005");
		graphPHOSYieldPi0SysErrPbPb0005Red->Write("InvYieldPbPbPHOSSysErr_0005");

		graphRAASysCombInd0010->Write("graphRAASysErr_0010");
		graphRAACombInd0010->Write("graphRAAStatErr_0010");
		graphRAASysPCM0010->Write("graphRAAPCMSysErr_0010");
		graphRAAPCM0010->Write("graphRAAPCMStatErr_0010");
		graphRAASysPHOS0010->Write("graphRAAPHOSSysErr_0010");
		graphRAAPHOS0010->Write("graphRAAPHOSStatErr_0010");
		graphYieldPi0CombPbPb0010StatErrUnscaled->Write("InvYieldPbPbStatErr_0010");
		graphYieldPi0CombPbPb0010SysErrUnscaled->Write("InvYieldPbPbSysErr_0010");
		graphYieldPi0CombPbPb0010FullSysErrUnscaled->Write("InvYieldPbPbFullSysErr_0010");
		graphPCMYieldPi0PbPb0010->Write("InvYieldPbPbPCMStatErr_0010");
		graphPCMYieldPi0SysErrPbPb0010Red->Write("InvYieldPbPbPCMSysErr_0010");
		graphPHOSYieldPi0PbPb0010->Write("InvYieldPbPbPHOSStatErr_0010");
		graphPHOSYieldPi0SysErrPbPb0010Red->Write("InvYieldPbPbPHOSSysErr_0010");
		
		graphRAASysCombInd0510->Write("graphRAASysErr_0510");
		graphRAACombInd0510->Write("graphRAAStatErr_0510");
		graphRAASysPCM0510->Write("graphRAAPCMSysErr_0510");
		graphRAAPCM0510->Write("graphRAAPCMStatErr_0510");
		graphRAASysPHOS0510->Write("graphRAAPHOSSysErr_0510");
		graphRAAPHOS0510->Write("graphRAAPHOSStatErr_0510");
		graphYieldPi0CombPbPb0510StatErrUnscaled->Write("InvYieldPbPbStatErr_0510");
		graphYieldPi0CombPbPb0510SysErrUnscaled->Write("InvYieldPbPbSysErr_0510");
		graphYieldPi0CombPbPb0510FullSysErrUnscaled->Write("InvYieldPbPbFullSysErr_0510");
		graphPCMYieldPi0PbPb0510->Write("InvYieldPbPbPCMStatErr_0510");
		graphPCMYieldPi0SysErrPbPb0510Red->Write("InvYieldPbPbPCMSysErr_0510");
		graphPHOSYieldPi0PbPb0510->Write("InvYieldPbPbPHOSStatErr_0510");
		graphPHOSYieldPi0SysErrPbPb0510Red->Write("InvYieldPbPbPHOSSysErr_0510");

		graphRAASysCombInd1020->Write("graphRAASysErr_1020");
		graphRAACombInd1020->Write("graphRAAStatErr_1020");
		graphRAASysPCM1020->Write("graphRAAPCMSysErr_1020");
		graphRAAPCM1020->Write("graphRAAPCMStatErr_1020");
		graphRAASysPHOS1020->Write("graphRAAPHOSSysErr_1020");
		graphRAAPHOS1020->Write("graphRAAPHOSStatErr_1020");
		graphYieldPi0CombPbPb1020StatErrUnscaled->Write("InvYieldPbPbStatErr_1020");
		graphYieldPi0CombPbPb1020SysErrUnscaled->Write("InvYieldPbPbSysErr_1020");
		graphYieldPi0CombPbPb1020FullSysErrUnscaled->Write("InvYieldPbPbFullSysErr_1020");
		graphPCMYieldPi0PbPb1020->Write("InvYieldPbPbPCMStatErr_1020");
		graphPCMYieldPi0SysErrPbPb1020Red->Write("InvYieldPbPbPCMSysErr_1020");
		graphPHOSYieldPi0PbPb1020->Write("InvYieldPbPbPHOSStatErr_1020");
		graphPHOSYieldPi0SysErrPbPb1020Red->Write("InvYieldPbPbPHOSSysErr_1020");
		

		graphRAASysCombInd2040->Write("graphRAASysErr_2040");
		graphRAACombInd2040->Write("graphRAAStatErr_2040");
		graphRAASysPCM2040->Write("graphRAAPCMSysErr_2040");
		graphRAAPCM2040->Write("graphRAAPCMStatErr_2040");
		graphRAASysPHOS2040->Write("graphRAAPHOSSysErr_2040");
		graphRAAPHOS2040->Write("graphRAAPHOSStatErr_2040");
		graphYieldPi0CombPbPb2040StatErrUnscaled->Write("InvYieldPbPbStatErr_2040");
		graphYieldPi0CombPbPb2040SysErrUnscaled->Write("InvYieldPbPbSysErr_2040");
		graphYieldPi0CombPbPb2040FullSysErrUnscaled->Write("InvYieldPbPbFullSysErr_2040");
		graphPCMYieldPi0PbPb2040->Write("InvYieldPbPbPCMStatErr_2040");
		graphPCMYieldPi0SysErrPbPb2040Red->Write("InvYieldPbPbPCMSysErr_2040");
		graphPHOSYieldPi0PbPb2040->Write("InvYieldPbPbPHOSStatErr_2040");
		graphPHOSYieldPi0SysErrPbPb2040Red->Write("InvYieldPbPbPHOSSysErr_2040");

		graphRAASysCombInd4060->Write("graphRAASysErr_4060");
		graphRAACombInd4060->Write("graphRAAStatErr_4060");
		graphRAASysPCM4060->Write("graphRAAPCMSysErr_4060");
		graphRAAPCM4060->Write("graphRAAPCMStatErr_4060");
		graphRAASysPHOS4060->Write("graphRAAPHOSSysErr_4060");
		graphRAAPHOS4060->Write("graphRAAPHOSStatErr_4060");
		graphYieldPi0CombPbPb4060StatErrUnscaled->Write("InvYieldPbPbStatErr_4060");
		graphYieldPi0CombPbPb4060SysErrUnscaled->Write("InvYieldPbPbSysErr_4060");
		graphYieldPi0CombPbPb4060SysErrUnscaled->Write("InvYieldPbPbFullSysErr_4060");
		graphPCMYieldPi0PbPb4060->Write("InvYieldPbPbPCMStatErr_4060");
		graphPCMYieldPi0SysErrPbPb4060Red->Write("InvYieldPbPbPCMSysErr_4060");
		
		graphPHOSYieldPi0PbPb4060->Write("InvYieldPbPbPHOSStatErr_4060");
		graphPHOSYieldPi0SysErrPbPb4060Red->Write("InvYieldPbPbPHOSSysErr_4060");

		graphRAASysCombInd6080->Write("graphRAASysErr_6080");
		graphRAACombInd6080->Write("graphRAAStatErr_6080");
		graphRAASysPCM6080->Write("graphRAAPCMSysErr_6080");
		graphRAAPCM6080->Write("graphRAAPCMStatErr_6080");
		graphRAASysPHOS6080->Write("graphRAAPHOSSysErr_6080");
		graphRAAPHOS6080->Write("graphRAAPHOSStatErr_6080");
		graphYieldPi0CombPbPb6080StatErrUnscaled->Write("InvYieldPbPbStatErr_6080");
		graphYieldPi0CombPbPb6080SysErrUnscaled->Write("InvYieldPbPbSysErr_6080");
		graphYieldPi0CombPbPb6080SysErrUnscaled->Write("InvYieldPbPbFullSysErr_6080");
		graphPCMYieldPi0PbPb6080->Write("InvYieldPbPbPCMStatErr_6080");
		graphPCMYieldPi0SysErrPbPb6080Red->Write("InvYieldPbPbPCMSysErr_6080");
		graphPHOSYieldPi0PbPb6080->Write("InvYieldPbPbPHOSStatErr_6080");
		graphPHOSYieldPi0SysErrPbPb6080Red->Write("InvYieldPbPbPHOSSysErr_6080");

		
		graphRatioCombPHOSPi0PbPb0005Sys->Write("graphRatioCombFitPHOSPi0SysPbPb0005");
		graphRatioCombConvPi0PbPb0005Sys->Write("graphRatioCombFitPCMPi0SysPbPb0005");
		graphRatioCombPHOSPi0PbPb0005->Write("graphRatioCombFitPHOSPi0StatPbPb0005");
		graphRatioCombConvPi0PbPb0005->Write("graphRatioCombFitPCMPi0StatPbPb0005");
		graphRatioCombPHOSPi0PbPb0510Sys->Write("graphRatioCombFitPHOSPi0SysPbPb0510");
		graphRatioCombConvPi0PbPb0510Sys->Write("graphRatioCombFitPCMPi0SysPbPb0510");
		graphRatioCombPHOSPi0PbPb0510->Write("graphRatioCombFitPHOSPi0StatPbPb0510");
		graphRatioCombConvPi0PbPb0510->Write("graphRatioCombFitPCMPi0StatPbPb0510");
		graphRatioCombPHOSPi0PbPb1020Sys->Write("graphRatioCombFitPHOSPi0SysPbPb1020");
		graphRatioCombConvPi0PbPb1020Sys->Write("graphRatioCombFitPCMPi0SysPbPb1020");
		graphRatioCombPHOSPi0PbPb1020->Write("graphRatioCombFitPHOSPi0StatPbPb1020");
		graphRatioCombConvPi0PbPb1020->Write("graphRatioCombFitPCMPi0StatPbPb1020");
		graphRatioCombPHOSPi0PbPb2040Sys->Write("graphRatioCombFitPHOSPi0SysPbPb2040");
		graphRatioCombConvPi0PbPb2040Sys->Write("graphRatioCombFitPCMPi0SysPbPb2040");
		graphRatioCombPHOSPi0PbPb2040->Write("graphRatioCombFitPHOSPi0StatPbPb2040");
		graphRatioCombConvPi0PbPb2040->Write("graphRatioCombFitPCMPi0StatPbPb2040");
		graphRatioCombPHOSPi0PbPb4060Sys->Write("graphRatioCombFitPHOSPi0SysPbPb4060");
		graphRatioCombConvPi0PbPb4060Sys->Write("graphRatioCombFitPCMPi0SysPbPb4060");
		graphRatioCombPHOSPi0PbPb4060->Write("graphRatioCombFitPHOSPi0StatPbPb4060");
		graphRatioCombConvPi0PbPb4060->Write("graphRatioCombFitPCMPi0StatPbPb4060");
		graphRatioCombPHOSPi0PbPb6080Sys->Write("graphRatioCombFitPHOSPi0SysPbPb6080");
		graphRatioCombConvPi0PbPb6080Sys->Write("graphRatioCombFitPCMPi0SysPbPb6080");
		graphRatioCombPHOSPi0PbPb6080->Write("graphRatioCombFitPHOSPi0StatPbPb6080");
		graphRatioCombConvPi0PbPb6080->Write("graphRatioCombFitPCMPi0StatPbPb6080");
		
		graphDeltaPtOverPtComb_PtPP_StatErr0005->Write("DeltaPtOverPt_PtPP_CombStatErr_0005");
		graphDeltaPtOverPtComb_PtPP_SysErr0005->Write("DeltaPtOverPt_PtPP_SysErr_0005");
		graphDeltaPtOverPtComb_PtPP_StatErr0510->Write("DeltaPtOverPt_PtPP_StatErr_0510");
		graphDeltaPtOverPtComb_PtPP_SysErr0510->Write("DeltaPtOverPt_PtPP_SysErr_0510");
		graphDeltaPtOverPtComb_PtPP_StatErr0010->Write("DeltaPtOverPt_PtPP_StatErr_0010");
		graphDeltaPtOverPtComb_PtPP_SysErr0010->Write("DeltaPtOverPt_PtPP_SysErr_0010");
		graphDeltaPtOverPtComb_PtPP_StatErr1020->Write("DeltaPtOverPt_PtPP_StatErr_1020");
		graphDeltaPtOverPtComb_PtPP_SysErr1020->Write("DeltaPtOverPt_PtPP_SysErr_1020");
		graphDeltaPtOverPtComb_PtPP_StatErr2040->Write("DeltaPtOverPt_PtPP_StatErr_2040");
		graphDeltaPtOverPtComb_PtPP_SysErr2040->Write("DeltaPtOverPt_PtPP_SysErr_2040");
		graphDeltaPtOverPtComb_PtPP_StatErr4060->Write("DeltaPtOverPt_PtPP_StatErr_4060");
		graphDeltaPtOverPtComb_PtPP_SysErr4060->Write("DeltaPtOverPt_PtPP_SysErr_4060");
		graphDeltaPtOverPtComb_PtPP_StatErr6080->Write("DeltaPtOverPt_PtPP_StatErr_6080");
		graphDeltaPtOverPtComb_PtPP_SysErr6080->Write("DeltaPtOverPt_PtPP_SysErr_6080");
		
		graphDeltaPtOverPtComb_PtPbPb_StatErr0005->Write("DeltaPtOverPt_PtPbPb_StatErr_0005");
		graphDeltaPtOverPtComb_PtPbPb_SysErr0005->Write("DeltaPtOverPt_PtPbPb_SysErr_0005");
		graphDeltaPtOverPtComb_PtPbPb_StatErr0510->Write("DeltaPtOverPt_PtPbPb_StatErr_0510");
		graphDeltaPtOverPtComb_PtPbPb_SysErr0510->Write("DeltaPtOverPt_PtPbPb_SysErr_0510");
		graphDeltaPtOverPtComb_PtPbPb_StatErr0010->Write("DeltaPtOverPt_PtPbPb_StatErr_0010");
		graphDeltaPtOverPtComb_PtPbPb_SysErr0010->Write("DeltaPtOverPt_PtPbPb_SysErr_0010");
		graphDeltaPtOverPtComb_PtPbPb_StatErr1020->Write("DeltaPtOverPt_PtPbPb_StatErr_1020");
		graphDeltaPtOverPtComb_PtPbPb_SysErr1020->Write("DeltaPtOverPt_PtPbPb_SysErr_1020");
		graphDeltaPtOverPtComb_PtPbPb_StatErr2040->Write("DeltaPtOverPt_PtPbPb_StatErr_2040");
		graphDeltaPtOverPtComb_PtPbPb_SysErr2040->Write("DeltaPtOverPt_PtPbPb_SysErr_2040");
		graphDeltaPtOverPtComb_PtPbPb_StatErr4060->Write("DeltaPtOverPt_PtPbPb_StatErr_4060");
		graphDeltaPtOverPtComb_PtPbPb_SysErr4060->Write("DeltaPtOverPt_PtPbPb_SysErr_4060");
		graphDeltaPtOverPtComb_PtPbPb_StatErr6080->Write("DeltaPtOverPt_PtPbPb_StatErr_6080");
		graphDeltaPtOverPtComb_PtPbPb_SysErr6080->Write("DeltaPtOverPt_PtPbPb_SysErr_6080");
		
		graphALICE2760GeVStat->Write("graphRAAvsNPart_StatErr");
		graphALICE2760GeVSys->Write("graphRAAvsNPart_SysErr");
		
		
	fCombResults.Close();


}
	
