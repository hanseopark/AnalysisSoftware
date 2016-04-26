#include <stdlib.h>
#include <iostream>
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
#include "TProfile2D.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
//#include "../CommonHeaders/FittingGammaConversion.h"
//#include "../CommonHeaders/PlottingMeson.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
//#include "../CommonHeaders/ConversionFunctions.h"
//#include "../CommonHeaders/ExtractSignalBinning.h"
//#include "../CommonHeaders/ExtractSignalPlotting.h"

TF1* FitRecursiveGaussian (TH1* histo, Double_t precision, Double_t correctRange, Double_t fitRangeMin, Double_t fitRangeMax);
TF1* FitBckg(TH1* fHisto, Double_t minFit, Double_t maxFit);
TF1* FitDataMC(TH1* fHisto, Double_t minFit, Double_t maxFit, Double_t constPar = -1);
Float_t	FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2);
template<class ForwardIt>
void SetLogBinningXTH(ForwardIt* histoRebin){
	TAxis *axisafter = histoRebin->GetXaxis();
	Int_t bins = axisafter->GetNbins();
	Double_t from = axisafter->GetXmin();
	Double_t to = axisafter->GetXmax();
	Double_t *newbins = new Double_t[bins+1];
	newbins[0] = from;
	Double_t factor = TMath::Power(to/from, 1./bins);
	for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
	axisafter->Set(bins, newbins);
	delete [] newbins;
	return;
}

//****************************************************************************
//************** Function to Correct CaloNonLinearity2 ***********************
//****************************************************************************
void CorrectCaloNonLinearity2(TString select = "LHC12-ConvCalo")
{
	gROOT->Reset();

	StyleSettingsThesis();
	SetPlotStyle();

	TString outputDir = "CorrectCaloNonLinearity2";
	gSystem->Exec("mkdir -p "+outputDir);

	TString suffix = "eps";
	TString optionEnergy = "";
	TString fPlot = "";
	TString fPlotTrigger = "";

	TString strDataFile = "";
	TString strTriggerDataFile = "";
	TString strMCFile = "";
	TString strTriggerMCFile = "";
	TString dataCut = "";
	TString mcCut = "";
	TString triggerDataCut = "";
	TString triggerMCCut = "";
	TString dataMainDir = "";
	TString mcMainDir = "";
	TString dataTriggerMainDir = "";
	TString mcTriggerMainDir = "";
	Int_t mode = 2;

	// pT range for mass fitting
	Int_t startPt = 0;
	Int_t endPt = 17;
	Int_t firstTriggerBin = -1;
	Int_t exampleBin1 = -1;
	Int_t exampleBin2 = -1;
//*******************************************************************************
// Choosing data set
	if(select.CompareTo("LHC11a-Pythia-Calo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaCalo_70.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaCalo_Pythia_70.root";
		dataCut = "00003113_1111100050032220000_0163103100000000";
		mcCut = "00003113_1111100050032220000_0163103100000000";
		dataMainDir = "GammaCalo_70";
		mcMainDir = "GammaCalo_70";

		firstTriggerBin = 16;
		strTriggerDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaCalo_6.root";
		strTriggerMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaCalo_JJ_70.root";
		triggerDataCut = "00051113_1111100050032230000_0163103100000050";
		triggerMCCut = "00003113_1111100050032220000_0163103100000000";
		dataTriggerMainDir = "GammaCalo_6";
		mcTriggerMainDir = "GammaCalo_70";

		optionEnergy = "2.76TeV";
		fPlot = "#frac{LHC12f1a & LHC12i3}{LHC11a - V0OR}";
		fPlotTrigger = "#frac{LHC15g1a}{LHC11a - EMC L0, INT1}";
		mode = 4;
		startPt = 4;
		endPt = 20;
	}else if(select.CompareTo("LHC11a-Pythia-ConvCalo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_1.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_Pythia_1.root";
		dataCut = "00003113_00200009327000008250400000_1111100053032230000_0163103100000000";
		mcCut = "00003113_00200009327000008250400000_1111100053032230000_0163103100000000";
		dataMainDir = "GammaConvCalo_1";
		mcMainDir = "GammaConvCalo_1";

		firstTriggerBin = 16;
		strTriggerDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_1.root";
		strTriggerMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_JJ_1.root";
		triggerDataCut = "00051113_00200009327000008250400000_1111100053032230000_0163103100000000";
		triggerMCCut = "00051113_00200009327000008250400000_1111100053032230000_0163103100000000";
		dataTriggerMainDir = "GammaConvCalo_1";
		mcTriggerMainDir = "GammaConvCalo_1";

		optionEnergy = "2.76TeV";
		fPlot = "#frac{LHC12f1a & LHC12i3}{LHC11a - V0OR}";
		fPlotTrigger = "#frac{LHC15g1a}{LHC11a - EMC L0, INT1}";
		mode = 2;
		startPt = 0;
		endPt = 21;

		exampleBin1 = 4;
		exampleBin2 = 19;
	}else if(select.CompareTo("LHC11a-Phojet-Calo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaCalo_70.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaCalo_Phojet_70.root";
		dataCut = "00003113_1111100050032220000_0163103100000000";
		mcCut = "00003113_1111100050032220000_0163103100000000";
		dataMainDir = "GammaCalo_70";
		mcMainDir = "GammaCalo_70";

		firstTriggerBin = 16;
		strTriggerDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaCalo_6.root";
		strTriggerMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaCalo_JJ_70.root";
		triggerDataCut = "00051113_1111100050032230000_0163103100000050";
		triggerMCCut = "00003113_1111100050032220000_0163103100000000";
		dataTriggerMainDir = "GammaCalo_6";
		mcTriggerMainDir = "GammaCalo_70";

		optionEnergy = "2.76TeV";
		fPlot = "#frac{LHC12f1b}{LHC11a - V0OR}";
		fPlotTrigger = "#frac{LHC15g1a}{LHC11a - EMC L0, INT1}";
		mode = 4;
		startPt = 4;
		endPt = 20;
	}else if(select.CompareTo("LHC11a-Phojet-ConvCalo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_1.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_Phojet_1.root";
		dataCut = "00003113_00200009327000008250400000_1111100053032230000_0163103100000000";
		mcCut = "00003113_00200009327000008250400000_1111100053032230000_0163103100000000";
		dataMainDir = "GammaConvCalo_1";
		mcMainDir = "GammaConvCalo_1";

		firstTriggerBin = 16;
		strTriggerDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_1.root";
		strTriggerMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_JJ_1.root";
		triggerDataCut = "00051113_00200009327000008250400000_1111100053032230000_0163103100000000";
		triggerMCCut = "00051113_00200009327000008250400000_1111100053032230000_0163103100000000";
		dataTriggerMainDir = "GammaConvCalo_1";
		mcTriggerMainDir = "GammaConvCalo_1";

		optionEnergy = "2.76TeV";
		fPlot = "#frac{LHC12f1b}{LHC11a - V0OR}";
        fPlotTrigger = "#frac{LHC15g1a}{LHC11a - EMC L0, INT1}";
		mode = 2;
		startPt = 0;
		endPt = 21;
	}else if(select.CompareTo("LHC12-Pythia-Calo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC12/GammaCalo_LHC12_p1_110.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC12/GammaCalo_Pythia_110.root";
		dataCut = "00000113_1111100063032230000_0163103100000050";
		mcCut = "00000113_1111100063032230000_0163103100000050";
		dataMainDir = "GammaCalo_110";
		mcMainDir = "GammaCalo_110";

		optionEnergy = "8TeV";
		fPlot = "#frac{LHC14e2a & LHC14e2b}{LHC12 - V0AND}";
		mode = 4;
		startPt = 4;
		endPt = 18;
	}else if(select.CompareTo("LHC12-Pythia-ConvCalo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC12/GammaConvCalo_LHC12_p1_111.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC12/GammaConvCalo_Pythia_101.root";
		dataCut = "00000113_00200009327000008250400000_1111100063032230000_0163103100000000";
		mcCut = "00000113_00200009327000008250400000_1111100063032230000_0163103100000000";
		dataMainDir = "GammaConvCalo_111";
		mcMainDir = "GammaConvCalo_101";

		optionEnergy = "8TeV";
		fPlot = "#frac{LHC14e2a & LHC14e2b}{LHC12 - V0AND}";
		mode = 2;
		startPt = 0;
		endPt = 19;
	}else if(select.CompareTo("LHC12-Phojet-Calo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC12/GammaCalo_LHC12_p1_110.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC12/GammaCalo_Phojet_110.root";
		dataCut = "00000113_1111100063032230000_0163103100000050";
		mcCut = "00000113_1111100063032230000_0163103100000050";
		dataMainDir = "GammaCalo_110";
		mcMainDir = "GammaCalo_110";

		optionEnergy = "8TeV";
		fPlot = "#frac{LHC14e2c}{LHC12 - V0AND}";
		mode = 4;
		startPt = 4;
		endPt = 18;
	}else if(select.CompareTo("LHC12-Phojet-ConvCalo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC12/GammaConvCalo_LHC12_p1_111.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC12/GammaConvCalo_Phojet_101.root";
		dataCut = "00000113_00200009327000008250400000_1111100063032230000_0163103100000000";
		mcCut = "00000113_00200009327000008250400000_1111100063032230000_0163103100000000";
		dataMainDir = "GammaConvCalo_111";
		mcMainDir = "GammaConvCalo_101";

		optionEnergy = "8TeV";
		fPlot = "#frac{LHC14e2c}{LHC12 - V0AND}";
		mode = 2;
		startPt = 0;
		endPt = 19;
	}else if(select.CompareTo("LHC13bc-Calo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaCalo_LHC13bc_1.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaCalo_LHC13b2_1.root";
		dataCut = "80000013_1111100050022230000_0163103100000000";
		mcCut = "80000013_1111100050022230000_0163103100000000";
		dataMainDir = "GammaCalo_1";
		mcMainDir = "GammaCalo_1";

		optionEnergy = "pPb_5.023TeV";
		fPlot = "LHC13b&c, V0AND";
		mode = 4;
		startPt = 4;
		endPt = 20;
//	}else if(select.CompareTo("LHC13bc-Calo")==0){
//		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaCalo_LHC13bc_1.root";
//		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaCalo_LHC13b2_1.root";
//		dataCut = "80000013_1111100050022230000_0163103100000000";
//		mcCut = "80000013_1111100050022230000_0163103100000000";
//		dataMainDir = "GammaCalo_1";
//		mcMainDir = "GammaCalo_1";

//		firstTriggerBin = 17;
//		strTriggerDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaCalo_LHC13_2.root";
//		strTriggerMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaCalo_LHC13b2_1.root";
//		triggerDataCut = "80083013_1111100050022230000_0163103100000000";
//		triggerMCCut = "80000013_1111100050022230000_0163103100000000";
//		dataTriggerMainDir = "GammaCalo_2";
//		mcTriggerMainDir = "GammaCalo_1";

//		optionEnergy = "pPb_5.023TeV";
//		fPlot = "#frac{LHC13b2_efix}{LHC13 - V0AND}";
//		fPlotTrigger = "#frac{LHC13b2_efix}{LHC13 - EMC L1, EG1}";
//		mode = 4;
//		startPt = 4;
//		endPt = 20;
//	}else if(select.CompareTo("LHC13bc-ConvCalo")==0){
//		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaConvCalo_LHC13bc_1.root";
//		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaConvCalo_LHC13b2_1.root";
//		dataCut = "80000013_00200009327002008250400000_1111100053022230000_0163103100000000";
//		mcCut = "80000013_00200009327002008250400000_1111100053022230000_0163103100000000";
//		dataMainDir = "GammaConvCalo_1";
//		mcMainDir = "GammaConvCalo_1";

//		optionEnergy = "pPb_5.023TeV";
//		fPlot = "#frac{LHC13b2_efix}{LHC13 - V0AND}";
//		mode = 2;
//		startPt = 0;
//		endPt = 21;
	}else if(select.CompareTo("LHC13bc-ConvCalo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaConvCalo_LHC13bc_1.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaConvCalo_LHC13b2_1.root";
		dataCut = "80000013_00200009327002008250400000_1111100053022230000_0163103100000000";
		mcCut = "80000013_00200009327002008250400000_1111100053022230000_0163103100000000";
		dataMainDir = "GammaConvCalo_1";
		mcMainDir = "GammaConvCalo_1";

		firstTriggerBin = 17;
		strTriggerDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaConvCalo_LHC13_2.root";
		strTriggerMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC13/GammaConvCalo_LHC13b2_1.root";
		triggerDataCut = "80083013_00200009327002008250400000_1111100053022230000_0163103100000000";
		triggerMCCut = "80000013_00200009327002008250400000_1111100053022230000_0163103100000000";
		dataTriggerMainDir = "GammaConvCalo_2";
		mcTriggerMainDir = "GammaConvCalo_1";

		optionEnergy = "pPb_5.023TeV";
		fPlot = "#frac{LHC13b2_efix}{LHC13 - V0AND}";
		fPlotTrigger = "#frac{LHC13b2_efix}{LHC13 - EMC L1, EG1}";
		mode = 2;
		startPt = 0;
		endPt = 21;
	}else if(select.CompareTo("LHC11a-JJ-Calo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaCalo_70.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaCalo_JJ_70.root";
		dataCut = "00003113_1111100050032220000_0163103100000000";
		mcCut = "00003113_1111100050032220000_0163103100000000";
		dataMainDir = "GammaCalo_70";
		mcMainDir = "GammaCalo_70";

		optionEnergy = "2.76TeV";
		fPlot = "#frac{LHC13b2_efix}{LHC13 - V0AND}";
		mode = 4;
		startPt = 4;
		endPt = 18;
	}else if(select.CompareTo("LHC11a-JJ-ConvCalo")==0){
		strDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_1.root";
		strMCFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_JJ_1.root";
		dataCut = "00003113_00200009327000008250400000_1111100053032230000_0163103100000000";
		mcCut = "00003113_00200009327000008250400000_1111100053032230000_0163103100000000";
		dataMainDir = "GammaConvCalo_1";
		mcMainDir = "GammaConvCalo_1";

		firstTriggerBin = 16;
		strTriggerDataFile = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150807-NonLinearity/LHC11a/GammaConvCalo_1.root";
		strTriggerMCFile = "";
		triggerDataCut = "00051113_00200009327000008250400000_1111100053032230000_0163103100000000";
		triggerMCCut = "";
		dataTriggerMainDir = "GammaConvCalo_1";
		mcTriggerMainDir = "";

		optionEnergy = "2.76TeV";
		fPlot = "#frac{LHC15g1a}{LHC11a - V0OR}";
		fPlotTrigger = "#frac{LHC15g1a}{LHC11a - EMC L0, INT1}";
		mode = 2;
		startPt = 0;
		endPt = 21;
	}else{
		cout << "No valid selection '" << select.Data() << "'' given, returning..." << endl;
		return;
	}
//*******************************************************************************

	// binning for fits Data vs MC mean mass position
	const Int_t fNBins = 24;
	Double_t fBins[fNBins+1] = {	0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6,
									1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0,
									5.0, 6.0, 8.0, 10.0, 12.0, 16.0,
									20.0, 25.0, 30.0};

	if(firstTriggerBin>-1){
		if(mode==2||mode==3){
			fPlot = Form("E_{Cluster} < %0.1f GeV : %s",fBins[firstTriggerBin],fPlot.Data());
			fPlotTrigger = Form("E_{Cluster} #geq %0.1f GeV : %s",fBins[firstTriggerBin],fPlotTrigger.Data());
		}else{
			fPlot = Form("#it{p}_{T} < %0.1f GeV/c : %s",fBins[firstTriggerBin],fPlot.Data());
			fPlotTrigger = Form("#it{p}_{T} #geq %0.1f GeV/c : %s",fBins[firstTriggerBin],fPlotTrigger.Data());
		}
	}

	TString recGamma = "";
	if(select.Contains("-Calo")) recGamma="#gamma's rec. with EMCal";
	else if(select.Contains("-ConvCalo")) recGamma="#gamma's rec. with PCM, EMCal";

	TString fTextMeasurement = Form("#pi^{0} #rightarrow #gamma#gamma");
	TString fCollisionSystem = ReturnFullCollisionsSystem(optionEnergy);
	if (fCollisionSystem.CompareTo("") == 0){
		cout << "No correct collision system specification, has been given" << endl;
		return;
	}
//*******************************************************************************
// Input
	TFile* dataFile = new TFile(strDataFile.Data(),"READ"); if(dataFile->IsZombie()) {cout << "Info: ROOT file '" << strDataFile.Data() << "' could not be openend, return!" << endl; return;}
	TList* dataTopDir = (TList*) dataFile->Get(dataMainDir.Data()); if(dataTopDir == NULL) {cout << "ERROR: dataTopDir not Found"<<endl; return;}
	TList* dataTopContainer = (TList*) dataTopDir->FindObject(Form("Cut Number %s",dataCut.Data())); if(dataTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",dataCut.Data()) << " not found in Data-File" << endl; return;}
	TList* dataESDContainer = (TList*) dataTopContainer->FindObject(Form("%s ESD histograms",dataCut.Data())); if(dataESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",dataCut.Data()) << " not found in Data-File" << endl; return;}

	TFile* mcFile = new TFile(strMCFile.Data(),"READ"); if(mcFile->IsZombie()) {cout << "Info: ROOT file '" << strMCFile.Data() << "' could not be openend, return!" << endl; return;}
	TList* mcTopDir = (TList*) mcFile->Get(mcMainDir.Data()); if(mcTopDir == NULL) {cout << "ERROR: mcTopDir not Found"<<endl; return;}
	TList* mcTopContainer = (TList*) mcTopDir->FindObject(Form("Cut Number %s",mcCut.Data())); if(mcTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",mcCut.Data()) << " not found in MC-File" << endl; return;}
	TList* mcESDContainer = (TList*) mcTopContainer->FindObject(Form("%s ESD histograms",mcCut.Data())); if(mcESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",mcCut.Data()) << " not found in MC-File" << endl; return;}
//*******************************************************************************
// Output
	TString nameOutput = "";
	nameOutput = Form("%s/CorrectCaloNonLinearity2_%s.root",outputDir.Data(),select.Data());
	TFile* fOutput = new TFile(nameOutput,"RECREATE");
//*******************************************************************************
// Fitting Data+MC
	TH2F* dataInvMassPtAlpha = NULL;
	TH2F* mcInvMassPtAlpha = NULL;
	if(mode==2||mode==3){
		dataInvMassPtAlpha = (TH2F*) dataESDContainer->FindObject("ESD_Mother_InvMass_E_Calib");
		mcInvMassPtAlpha = (TH2F*) mcESDContainer->FindObject("ESD_Mother_InvMass_E_Calib");
	}else{
		dataInvMassPtAlpha = (TH2F*) dataESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt_Alpha");
		mcInvMassPtAlpha = (TH2F*) mcESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt_Alpha");
	}
	if(!dataInvMassPtAlpha){cout << "did not find ESD_Mother_InvMass_E_Calib in data" << endl; return;}
	if(!mcInvMassPtAlpha){cout << "did not find ESD_Mother_InvMass_E_Calib in mc" << endl; return;}

	dataInvMassPtAlpha->Write("Data - ESD_Mother_InvMass");
	mcInvMassPtAlpha->Write("MC - ESD_Mother_InvMass");

	TH1D* histMCResults = NULL;
	TH1D* histDataResults = NULL;
	TH1D* histDataMCResults = NULL;
	if(mode==2||mode==3){
		histMCResults = new TH1D("Mean mass MC","; E_{Cluster} (GeV); mean mass MC",fNBins,fBins);
		histDataResults = new TH1D("Mean mass Data","; E_{Cluster} (GeV); mean mass Data",fNBins,fBins);
		histDataMCResults = new TH1D("Mean mass ratio MC/Data","; E_{Cluster} (GeV); mean mass ratio (MC/Data)",fNBins,fBins);
	}else{
		histMCResults = new TH1D("Mean mass MC","; #it{p}_{T} (GeV/c); mean #pi^{0} mass MC",fNBins,fBins);
		histDataResults = new TH1D("Mean mass Data","; #it{p}_{T} (GeV/c); mean #pi^{0} mass Data",fNBins,fBins);
		histDataMCResults = new TH1D("Mean mass ratio MC/Data","; #it{p}_{T} (GeV/c); mean #pi^{0} mass ratio (MC/Data)",fNBins,fBins);
	}
	histMCResults->SetDirectory(0);
	histDataResults->SetDirectory(0);
	histMCResults->GetXaxis()->SetRangeUser(fBins[startPt],fBins[endPt]);
	histDataResults->GetXaxis()->SetRangeUser(fBins[startPt],fBins[endPt]);
	histMCResults->GetYaxis()->SetRangeUser(0.12,0.14);
	histDataResults->GetYaxis()->SetRangeUser(0.12,0.14);

	TF1* fFitReco;
	TF1* fFitMassPos;

	TCanvas *canvas = new TCanvas("canvas","",200,10,1350,900);  // gives the page size
	Double_t leftMargin = 0.1; Double_t rightMargin = 0.02; Double_t topMargin = 0.06; Double_t bottomMargin = 0.1;
	DrawGammaCanvasSettings(canvas, leftMargin, rightMargin, topMargin, bottomMargin);

	TString dataMC[2]={"Data","MC"};

	for(Int_t iClusterPt=startPt; iClusterPt<endPt; iClusterPt++)
	{
		if(iClusterPt==firstTriggerBin){
			// switching to trigger
			cout << endl;
			cout << "-----------------------------------------------------" << endl;
			cout << "\t Closing open files, switching to Trigger!" << endl;
			cout << "bin: " << firstTriggerBin << endl;
			cout << "-----------------------------------------------------" << endl;

			if(!strTriggerDataFile.IsNull()){
				dataESDContainer->Clear(); dataTopContainer->Clear(); dataTopDir->Clear(); dataFile->Delete();
				dataFile = new TFile(strTriggerDataFile.Data(),"READ"); if(dataFile->IsZombie()) {cout << "Info: ROOT file '" << strTriggerDataFile.Data() << "' could not be openend, return!" << endl; return;}
				dataTopDir = (TList*) dataFile->Get(dataTriggerMainDir.Data()); if(dataTopDir == NULL) {cout << "ERROR: dataTopDir not Found"<<endl; return;}
				dataTopContainer = (TList*) dataTopDir->FindObject(Form("Cut Number %s",triggerDataCut.Data())); if(dataTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",triggerDataCut.Data()) << " not found in File" << endl; return;}
				dataESDContainer = (TList*) dataTopContainer->FindObject(Form("%s ESD histograms",triggerDataCut.Data())); if(dataESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",triggerDataCut.Data()) << " not found in File" << endl; return;}
			}
			if(!strTriggerMCFile.IsNull()){
				mcESDContainer->Clear(); mcTopContainer->Clear(); mcTopDir->Clear(); mcFile->Delete();
				mcFile = new TFile(strTriggerMCFile.Data(),"READ"); if(mcFile->IsZombie()) {cout << "Info: ROOT file '" << strTriggerMCFile.Data() << "' could not be openend, return!" << endl; return;}
				mcTopDir = (TList*) mcFile->Get(mcTriggerMainDir.Data()); if(mcTopDir == NULL) {cout << "ERROR: mcTopDir not Found"<<endl; return;}
				mcTopContainer = (TList*) mcTopDir->FindObject(Form("Cut Number %s",triggerMCCut.Data())); if(mcTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",triggerMCCut.Data()) << " not found in File" << endl; return;}
				mcESDContainer = (TList*) mcTopContainer->FindObject(Form("%s ESD histograms",triggerMCCut.Data())); if(mcESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",triggerMCCut.Data()) << " not found in File" << endl; return;}
			}
			if(mode==2||mode==3){
				dataInvMassPtAlpha = (TH2F*) dataESDContainer->FindObject("ESD_Mother_InvMass_E_Calib");
				mcInvMassPtAlpha = (TH2F*) mcESDContainer->FindObject("ESD_Mother_InvMass_E_Calib");
			}else{
				dataInvMassPtAlpha = (TH2F*) dataESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt_Alpha");
				mcInvMassPtAlpha = (TH2F*) mcESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt_Alpha");
			}
			if(!dataInvMassPtAlpha){cout << "did not find ESD_Mother_InvMass_E_Calib in data" << endl; return;}
			if(!mcInvMassPtAlpha){cout << "did not find ESD_Mother_InvMass_E_Calib in mc" << endl; return;}

			fOutput->cd();
		}

		cout << endl;
		cout << "-----------------------------------------------------" << endl;
		cout << "\t MC/Data Fitting mass positions" << endl;
		cout << "loop: " << iClusterPt << ", " << fBins[iClusterPt] << " - " << fBins[iClusterPt+1] << " GeV" << endl;
		cout << "-----------------------------------------------------" << endl;

		TH2* Hist2D;
		for(Int_t iDataMC=0;iDataMC<2;iDataMC++)
		{
			if(iDataMC==0) Hist2D = dataInvMassPtAlpha;
			else if(iDataMC==1) Hist2D = mcInvMassPtAlpha;
			else{cout << "ERROR: data/mc loop, returning..." << endl; return;}

			Double_t projectMin = Hist2D->GetYaxis()->FindBin(fBins[iClusterPt]+0.001);
			Double_t projectMax = Hist2D->GetYaxis()->FindBin(fBins[iClusterPt+1]-0.001);
			TH1D* sliceHist = (TH1D*) Hist2D->ProjectionX(Form("slice%sAlpha_%f-%f",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]),projectMin,projectMax);
			sliceHist->SetDirectory(0);
            if(select.Contains("ConvCalo")) sliceHist->SetTitle(Form("%s - %.01f < E_{Cluster} < %.01f (GeV)",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
            else sliceHist->SetTitle(Form("%s - %.01f < #it{p}_{T} < %.01f (GeV/#it{c})",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
			sliceHist->GetYaxis()->SetTitle("#frac{d#it{M}_{inv}}{dN}");
//*******************************************************************************
// Rebin
			 if( select.Contains("LHC11a-Pythia-ConvCalo") || select.Contains("LHC11a-Phojet-ConvCalo") || select.Contains("LHC11a-JJ-ConvCalo")){
				if(fBins[iClusterPt]>=12) sliceHist->Rebin(8);
				else if(fBins[iClusterPt]>=10) sliceHist->Rebin(4);
				else sliceHist->Rebin(2);
			}else if(select.Contains("LHC11a-Pythia-Calo") || select.Contains("LHC11a-Phojet-Calo") || select.Contains("LHC11a-JJ-Calo")){
				if(fBins[iClusterPt]>=5) sliceHist->Rebin(5);
				else if(fBins[iClusterPt]>=3) sliceHist->Rebin(4);
				else if(fBins[iClusterPt]<1) sliceHist->Rebin(4);
				else sliceHist->Rebin(2);
			}else if(select.Contains("LHC13bc") && iDataMC == 0) sliceHist->Rebin(2);
			 else if(select.Contains("LHC13bc") && iDataMC == 1 && fBins[iClusterPt]<10) sliceHist->Rebin(2);
			 else if(select.Contains("LHC13bc") && iDataMC == 1 && fBins[iClusterPt]>=10) sliceHist->Rebin(4);
/*			 else if(select.Contains("LHC13bc") && fBins[iClusterPt]>=12){
				 if(select.Contains("LHC13bc-Calo")) sliceHist->Rebin(8);
				 if(select.Contains("LHC13bc-ConvCalo")) sliceHist->Rebin(5);
			 }*/else if(fBins[iClusterPt]>=8) sliceHist->Rebin(5);
			 else sliceHist->Rebin(4);
//*******************************************************************************
// Background subtraction ranges
			if( ( ( select.Contains("LHC11a-Pythia-ConvCalo") || select.Contains("LHC11a-Phojet-ConvCalo") || select.Contains("LHC11a-JJ-ConvCalo")) && iDataMC==0 && fBins[iClusterPt]<3.2 )
				|| ( (select.Contains("LHC11a-Pythia-ConvCalo") || select.Contains("LHC11a-Phojet-ConvCalo")) && iDataMC==1 && fBins[iClusterPt]<2.8 )
				|| ( (select.Contains("LHC11a-Pythia-ConvCalo") || select.Contains("LHC11a-Phojet-ConvCalo")) && iDataMC==1 && fBins[iClusterPt]>=5 )
				|| ( select.Contains("LHC11a-JJ-ConvCalo") && iDataMC==1 )

				|| ( (select.Contains("LHC11a-Pythia-Calo") || select.Contains("LHC11a-Phojet-Calo")) && fBins[iClusterPt]<8 )
				|| ( select.Contains("LHC11a-JJ-Calo") && iDataMC==0 && fBins[iClusterPt]<6 )
				|| ( select.Contains("LHC11a-JJ-Calo") && iDataMC==1 )

				|| ( select.Contains("LHC12") && select.Contains("ConvCalo") && fBins[iClusterPt]<2.8 )
				|| ( select.Contains("LHC12") && !select.Contains("ConvCalo") && fBins[iClusterPt]<6 )
				|| ( select.Contains("LHC13bc-ConvCalo") && iDataMC==0 )
				|| ( select.Contains("LHC13bc-Calo") && fBins[iClusterPt]<5 && iDataMC==0 )
				|| ( select.Contains("LHC13") && fBins[iClusterPt]<5 && iDataMC==1 )
				){
				Double_t range = 0.2;
				if(select.Contains("ConvCalo")) range = 0.3;
				TF1* fBckFit = FitBckg(sliceHist,0.03,range);
				sliceHist->GetListOfFunctions()->Add(fBckFit);
				sliceHist->Write(Form("slice%sAlpha_%f-%f-withBckgAndFit",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
				if(iClusterPt==exampleBin1 || iClusterPt==exampleBin2){
					sliceHist->GetXaxis()->SetRangeUser(0.05,0.3);
					sliceHist->DrawCopy();
					canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
                    canvas->SaveAs(Form("%s/slice%sAlpha_%.01f-%.01f-withBckgAndFit.%s",outputDir.Data(),dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1],suffix.Data()));
					canvas->Clear();
				}
				sliceHist->GetXaxis()->SetRangeUser(0.05,0.2);

				TF1* fBckg = new TF1(Form("fBckg%s",dataMC[iDataMC].Data()),"[0]+[1]*x+[2]*x*x",0.05,0.2);
				fBckg->SetParameter(0,fBckFit->GetParameter(0));
				fBckg->SetParameter(1,fBckFit->GetParameter(1));
				fBckg->SetParameter(2,fBckFit->GetParameter(2));

				sliceHist->Add(fBckg,-1);
			}else sliceHist->GetXaxis()->SetRangeUser(0.05,0.2);

			Double_t sigmaRangeAdjust = 1;
			Double_t precision = 0.1;
			Double_t minMax[2]={0.1,0.15};
//*******************************************************************************
// Adjusting sigma fitting
			if(select.Contains("LHC11a-Pythia-ConvCalo") || select.Contains("LHC11a-Phojet-ConvCalo") || select.Contains("LHC11a-JJ-ConvCalo")){
				if(fBins[iClusterPt]>=2 && fBins[iClusterPt]<3.2) sigmaRangeAdjust = 1.5;
				else if(fBins[iClusterPt]>=3.2 && fBins[iClusterPt]<5) sigmaRangeAdjust = 2;
				else if(fBins[iClusterPt]>=5 && fBins[iClusterPt]<10) sigmaRangeAdjust = 1;
				else if(fBins[iClusterPt]>=10) sigmaRangeAdjust = 2;

			}else if(select.Contains("LHC11a-Pythia-Calo") || select.Contains("LHC11a-Phojet-Calo") || select.Contains("LHC11a-JJ-Calo")){
				sigmaRangeAdjust = 1.5;
				if(fBins[iClusterPt]>=6) sigmaRangeAdjust = 2;

			}else if(select.Contains("LHC12")){
				if(fBins[iClusterPt]<8) sigmaRangeAdjust = 1.5;
				else if(fBins[iClusterPt]>=8) sigmaRangeAdjust = 2.5;

				if(fBins[iClusterPt]==6 && select.Contains("Phojet-ConvCalo") && iDataMC==1) sigmaRangeAdjust = 2;

			}else if(select.Contains("LHC13bc-ConvCalo")){
				sigmaRangeAdjust = 1;
				if(fBins[iClusterPt]>=1) sigmaRangeAdjust = 1.5;
			}else if(select.Contains("LHC13bc-Calo")){
				sigmaRangeAdjust = 1;
//				if(iDataMC==0){
//					if(fBins[iClusterPt]==6) sigmaRangeAdjust = 0.5;
//					else if(fBins[iClusterPt]==8) sigmaRangeAdjust = 1;
//					else if(fBins[iClusterPt]==10) sigmaRangeAdjust = 1;
//				}
				if(fBins[iClusterPt]>=1) sigmaRangeAdjust = 1.5;
			}
//*******************************************************************************
// Fit
			if((select.Contains("LHC11a-Pythia-ConvCalo") || select.Contains("LHC11a-Phojet-ConvCalo") || select.Contains("LHC11a-JJ-ConvCalo")) && fBins[iClusterPt]==8){minMax[0]=0.13; minMax[1]=0.15;}
			if(select.Contains("LHC12-Phojet-ConvCalo") && fBins[iClusterPt]==8){minMax[0]=0.12; minMax[1]=0.15;}
			if(select.Contains("LHC13bc-ConvCalo") && fBins[iClusterPt]==12){minMax[0]=0.12; minMax[1]=0.15;}
			if(select.Contains("LHC13bc-ConvCalo") && iDataMC == 1 && fBins[iClusterPt]==16){minMax[0]=0.1; minMax[1]=0.17;}
			if(select.Contains("LHC13bc-Calo") && fBins[iClusterPt]>=12){minMax[0]=0.11; minMax[1]=0.17;}

			fFitReco = FitRecursiveGaussian (sliceHist, precision, sigmaRangeAdjust, minMax[0], minMax[1]);

			if(iDataMC==0) {
				histDataResults->SetBinContent(iClusterPt+1,fFitReco->GetParameter(1));
				histDataResults->SetBinError(iClusterPt+1,fFitReco->GetParError(1));
			}
			else if(iDataMC==1) {
				histMCResults->SetBinContent(iClusterPt+1,fFitReco->GetParameter(1));
				histMCResults->SetBinError(iClusterPt+1,fFitReco->GetParError(1));
			}

			sliceHist->GetListOfFunctions()->Add(fFitReco);
			sliceHist->Write();
			if(iClusterPt==exampleBin1 || iClusterPt==exampleBin2){
				sliceHist->DrawCopy();
				canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
                canvas->SaveAs(Form("%s/slice%sAlpha_%.01f-%.01f.%s",outputDir.Data(),dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1],suffix.Data()));
				canvas->Clear();
			}
		}
	}

	cout << endl;
	cout << "-----------------------------------------------------" << endl;
	cout << "-----------------------------------------------------" << endl;

	histMCResults->Write();
	histDataResults->Write();

	histDataMCResults->Divide(histMCResults,histDataResults,1,1);
	DrawGammaSetMarker(histDataMCResults, 24, 2, kBlack, kBlack);

	Double_t startFit = 6;
	if(select.Contains("LHC12") && select.Contains("ConvCalo")) startFit = 3;
	else if(select.Contains("LHC12")) startFit = 4;

	TF1* fFitConst = new TF1("DataMCConst", "[0]" ,startFit,fBins[endPt]);
	histDataMCResults->Fit(fFitConst,"QRME0");

	fFitMassPos = FitDataMC(histDataMCResults, fBins[startPt], fBins[endPt],fFitConst->GetParameter(0));

	histDataMCResults->GetYaxis()->SetRangeUser(0.95,1.05);
	histDataMCResults->GetXaxis()->SetRangeUser(fBins[startPt],fBins[endPt]);

	histDataMCResults->Write("MeanMassRatioMCData-noFit");
	fFitMassPos->Write("MeanMassRatioMCData-Fit");
	histDataMCResults->GetListOfFunctions()->Add(fFitMassPos);

	fstream fLog;
	fLog.open(Form("%s/CorrectCaloNonLinearity2_%s.log",outputDir.Data(),select.Data()), ios::out);
	fLog << "FitDataMC results:" << endl;
	for(Int_t i=0;i<=2;i++) fLog << "Par " << i << ": " << fFitMassPos->GetParameter(i) << " +- " << fFitMassPos->GetParError(i) << endl;
	fLog.close();
//*******************************************************************************
// plotting mass ratios
	histDataMCResults->GetYaxis()->SetLabelFont(42);
	histDataMCResults->GetXaxis()->SetLabelFont(42);
	histDataMCResults->GetYaxis()->SetTitleFont(62);
	histDataMCResults->GetXaxis()->SetTitleFont(62);
	histDataMCResults->GetYaxis()->SetLabelSize(0.035);
	histDataMCResults->GetYaxis()->SetTitleSize(0.043);
	histDataMCResults->GetYaxis()->SetDecimals();
	histDataMCResults->GetYaxis()->SetTitleOffset(1);
	histDataMCResults->GetXaxis()->SetTitleSize(0.043);
	histDataMCResults->GetXaxis()->SetLabelSize(0.035);
	histDataMCResults->GetXaxis()->SetTitleOffset(0.9);
	histDataMCResults->Draw();
	fFitMassPos->Draw("same");
	PutProcessLabelAndEnergyOnPlot(0.7, 0.89, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data());
	PutProcessLabelAndEnergyOnPlot(0.2, 0.89, 0.03, fPlot.Data(),"", "");
	PutProcessLabelAndEnergyOnPlot(0.2, 0.82, 0.03, fPlotTrigger.Data(),"", "");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/MeanMassRatio_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
//*******************************************************************************
// plotting total correction
	TH1D* totalCorrection = new TH1D("Total Correction","; E_{Cluster} (GeV); correction factor",1000,0.4,50);
	totalCorrection->GetYaxis()->SetLabelFont(42);
	totalCorrection->GetXaxis()->SetLabelFont(42);
	totalCorrection->GetYaxis()->SetTitleFont(62);
	totalCorrection->GetXaxis()->SetTitleFont(62);
	totalCorrection->GetYaxis()->SetLabelSize(0.035);
	totalCorrection->GetYaxis()->SetTitleSize(0.043);
	totalCorrection->GetYaxis()->SetDecimals();
	totalCorrection->GetYaxis()->SetTitleOffset(1);
	totalCorrection->GetXaxis()->SetTitleSize(0.043);
	totalCorrection->GetXaxis()->SetLabelSize(0.035);
	totalCorrection->GetXaxis()->SetTitleOffset(0.9);
	totalCorrection->GetYaxis()->SetRangeUser(0.95,1.1);
	totalCorrection->GetXaxis()->SetMoreLogLabels();
	totalCorrection->GetXaxis()->SetNoExponent();
	SetLogBinningXTH(totalCorrection);
	DrawGammaSetMarker(totalCorrection, 8, 0.8, kBlack, kBlack);

	for(Int_t iBin=1; iBin<totalCorrection->GetNbinsX()+1;iBin++){
		Float_t e = totalCorrection->GetXaxis()->GetBinCenter(iBin);
		Float_t factor = 1;
		if(mode==2||mode==3) factor /= FunctionNL_kSDM(e,fFitMassPos->GetParameter(0),fFitMassPos->GetParameter(1),fFitMassPos->GetParameter(2));
		else factor /= FunctionNL_kSDM(2.0*e,fFitMassPos->GetParameter(0),fFitMassPos->GetParameter(1),fFitMassPos->GetParameter(2));
		totalCorrection->SetBinContent(iBin,factor);
	}
	totalCorrection->DrawCopy("p");
	TLegend *legend = new TLegend(0.2,0.65,0.7,0.78);
	legend->SetNColumns(1);
	legend->SetFillColor(0);
	legend->SetLineColor(0);
	legend->SetTextSize(0.03);
	legend->SetTextFont(42);
	legend->AddEntry(totalCorrection,"Correction factor for MC");
	legend->Draw("same");

	PutProcessLabelAndEnergyOnPlot(0.7, 0.89, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data());
	PutProcessLabelAndEnergyOnPlot(0.2, 0.89, 0.03, fPlot.Data(),"", "");
	PutProcessLabelAndEnergyOnPlot(0.2, 0.82, 0.03, fPlotTrigger.Data(),"", "");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/TotalCorrection_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();

	cout << "-----------------------------------------------------" << endl;
	cout << "-----------------------------------------------------" << endl;

	fOutput->Write();
	fOutput->Close();

	delete canvas;

	return;
}

//*******************************************************************************
//*******************************************************************************
//*******************************************************************************

Float_t FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2){
	return ( p0 + exp( p1 + ( p2 * e ) ) );
}

Double_t fitExcludeSignal(Double_t *x, Double_t *par)
{
	if (x[0] > 0.08 && x[0] < 0.17) {
	  TF1::RejectPoint();
	  return 0;
   }
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}


TF1* FitBckg(TH1* fHisto, Double_t minFit, Double_t maxFit){
	TF1* fFitBckg = new TF1("fFitBckg",fitExcludeSignal,minFit,maxFit,3);
	fFitBckg->SetLineColor(kBlue);
	fFitBckg->SetLineWidth(2);
	fFitBckg->SetLineStyle(1);
	fHisto->Fit(fFitBckg,"QRME0");
	return fFitBckg;
}

TF1* FitDataMC(TH1* fHisto, Double_t minFit, Double_t maxFit, Double_t constPar){

	TF1* fFitReco = new TF1("DataMC", "[0]+exp([1]+([2]*x))" ,minFit,maxFit);

	fFitReco->SetParameter(0,1.);
	fFitReco->SetParameter(1,-1.);
	fFitReco->SetParameter(2,-0.5);
	if(constPar!=-1) fFitReco->FixParameter(0,constPar);

	fHisto->Fit(fFitReco,"QRME0");

	fFitReco->SetLineColor(kRed);
	fFitReco->SetLineWidth(2);
	fFitReco->SetLineStyle(1);

	if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
		cout << "Parameters for DataMC: " << endl;
		for(Int_t i=0;i<=2;i++) cout << "Par " << i << ": " << fFitReco->GetParameter(i) << " +- " << fFitReco->GetParError(i) << endl;
	} else {
		cout << "DataMC fitting failed with status " << gMinuit->fCstatu.Data() <<endl << endl;
	}

	return fFitReco;
}

TF1* FitRecursiveGaussian (TH1* histo, Double_t precision, Double_t correctRange, Double_t fitRangeMin, Double_t fitRangeMax ){
	TF1 *f0 = new TF1("f0", "gaus", fitRangeMin,fitRangeMax);
	histo->Fit(f0,"0RMEQ");
	Double_t rp = f0->GetParameter(2);
	Double_t mp = f0->GetParameter(1);
	Double_t ymin = mp -(rp * correctRange);
	Double_t ymax = mp + (rp * correctRange);
	Double_t deviation = 100;
	Int_t counter = 0;
	TF1* f1 = new TF1 ("f1", "gaus", ymin, ymax);
	while(deviation > precision && counter < 100){
		f1->SetRange(ymin,ymax);
		histo->Fit(f1,"0RMEQ");
		Double_t rp2 = f1->GetParameter(2);
		if (rp2>rp){ deviation = rp2-rp;}
			else {deviation = rp -rp2 ;}
		rp = rp2 ;
		mp = f1->GetParameter(1);
		ymin = mp -(rp * correctRange);
		ymax = mp +(rp * correctRange);
		counter++;
	}
	delete f0;

	if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
		cout << "Parameters for FitRecursiveGaussian: " << endl;
		for(Int_t i=0;i<=2;i++) cout << "Par " << i << ": " << f1->GetParameter(i) << " +- " << f1->GetParError(i) << endl;
	} else {
		cout << "FitRecursiveGaussian fitting failed with status " << gMinuit->fCstatu.Data() <<endl << endl;
	}

	return f1;
}
