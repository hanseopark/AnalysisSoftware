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

//****************************************************************************
//************** Function to compare different CaloNonLinearities ************
//****************************************************************************
void CorrectCaloNonLinearity_Compare(TString select = "LHC12")
{
	gROOT->Reset();

	StyleSettingsThesis();
	SetPlotStyle();

	TString inputDir = "CorrectCaloNonLinearity";
	TString outputDir = Form("%s/Compare",inputDir.Data());
	gSystem->Exec("mkdir -p "+outputDir);

	TString suffix = "eps";

	Int_t nNL;
	TString *inputFileNames;
	TString *inputFilePaths;
	TFile **inputFiles;

	Color_t color[13] = {kBlack,633,807,800,418,kGreen+4,435,601,879,806,852,kCyan+3,426};

//*******************************************************************************
// Choosing data set
	if(select.CompareTo("LHC12")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC12-ConvCalo";
		inputFileNames[1] = "LHC12-Calo";
	}else if(select.CompareTo("LHC11a")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC11a-ConvCalo";
		inputFileNames[1] = "LHC11a-Calo";
	}else if(select.CompareTo("LHC13bc")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC13bc-ConvCalo";
		inputFileNames[1] = "LHC13bc-Calo";
	}else if(select.CompareTo("ConvCalo")==0){
		nNL = 3;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC11a-ConvCalo";
		inputFileNames[1] = "LHC12-ConvCalo";
		inputFileNames[2] = "LHC13bc-ConvCalo";
	}else if(select.CompareTo("Calo")==0){
		nNL = 3;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC11a-Calo";
		inputFileNames[1] = "LHC12-Calo";
		inputFileNames[2] = "LHC13bc-Calo";
	}else{
		cout << "No valid selection '" << select.Data() << "'' given, returning..." << endl;
		return;
	}

//*******************************************************************************
// Input
	inputFilePaths = new TString[nNL];
	inputFiles = new TFile*[nNL];
	for(Int_t i=0; i<nNL; i++){
		inputFilePaths[i] = Form("%s/CorrectCaloNonLinearity_%s.root",inputDir.Data(),inputFileNames[i].Data());
		inputFiles[i] = new TFile(inputFilePaths[i].Data(),"READ");
			if(inputFiles[i]->IsZombie()) {cout << "Info: ROOT file '" << inputFilePaths[i].Data() << "' could not be openend, return!" << endl; return;}
	}

//*******************************************************************************
// Output
	TString nameOutput = "";
	nameOutput = Form("%s/CorrectCaloNonLinearity_%s.root",outputDir.Data(),select.Data());
	TFile* fOutput = new TFile(nameOutput,"RECREATE");

	TCanvas *canvas = new TCanvas("canvas","",200,10,1350,900);  // gives the page size
	Double_t leftMargin = 0.1; Double_t rightMargin = 0.02; Double_t topMargin = 0.06; Double_t bottomMargin = 0.1;
	DrawGammaCanvasSettings(canvas, leftMargin, rightMargin, topMargin, bottomMargin);

	TLegend *legend = new TLegend(0.3,0.8,0.5,0.88);
	legend->SetNColumns(1);
	legend->SetFillColor(0);
	legend->SetLineColor(0);
	legend->SetTextSize(0.03);
	legend->SetTextFont(42);
//*******************************************************************************
// plotting masses MC+Data
	TH1D* histoMassMC[nNL];
	for(Int_t i=0; i<nNL; i++){
		TString drawOption = (i==0)?"p":"p, same";
		histoMassMC[i] = (TH1D*)inputFiles[i]->Get("Mean mass MC");
		histoMassMC[i]->SetMarkerColor(color[i]);
		histoMassMC[i]->SetMarkerStyle(24+i);
		histoMassMC[i]->SetLineColor(color[i]);
		histoMassMC[i]->Draw(drawOption.Data());
		histoMassMC[i]->GetXaxis()->SetTitleOffset(1.1);
		histoMassMC[i]->GetXaxis()->SetTitle("ConvCalo: E_{Cluster} (GeV) - Calo: #it{p}^{#pi^{0}}_{T} (GeV/c)");
		legend->AddEntry(histoMassMC[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/Mass_MC_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

	TH1D* histoMassData[nNL];
	for(Int_t i=0; i<nNL; i++){
		TString drawOption = (i==0)?"p":"p, same";
		histoMassData[i] = (TH1D*)inputFiles[i]->Get("Mean mass Data");
		histoMassData[i]->SetMarkerColor(color[i]);
		histoMassData[i]->SetMarkerStyle(24+i);
		histoMassData[i]->SetLineColor(color[i]);
		histoMassData[i]->Draw(drawOption.Data());
		histoMassData[i]->GetXaxis()->SetTitleOffset(1.1);
		histoMassData[i]->GetXaxis()->SetTitle("ConvCalo: E_{Cluster} (GeV) - Calo: #it{p}^{#pi^{0}}_{T} (GeV/c)");
		legend->AddEntry(histoMassData[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/Mass_Data_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();
	return;
//*******************************************************************************
// plotting E_truth/E_rec
	TH1D* histoEtruthErec[nNL];
	for(Int_t i=0; i<nNL; i++){
		TString drawOption = (i==0)?"p":"p, same";
		histoEtruthErec[i] = (TH1D*)inputFiles[i]->Get("MeanEtruthErec-noFit");
		histoEtruthErec[i]->SetMarkerColor(color[i]);
		histoEtruthErec[i]->SetMarkerStyle(24+i);
		histoEtruthErec[i]->SetLineColor(color[i]);
		histoEtruthErec[i]->Draw(drawOption.Data());
		legend->AddEntry(histoEtruthErec[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/E_truth_recVsE_rec_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

	TF1* histoEtruthErecFits[nNL];
	histoEtruthErec[0]->Draw("axis");
	for(Int_t i=0; i<nNL; i++){
		histoEtruthErecFits[i] = (TF1*)inputFiles[i]->Get("MeanEtruthErec-Fit");
		histoEtruthErecFits[i]->SetLineColor(color[i]);
		histoEtruthErecFits[i]->Draw("same");
		legend->AddEntry(histoEtruthErecFits[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/E_truth_recVsE_rec_Fits_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

	cout << "-----------------------------------------------------" << endl;
	cout << "-----------------------------------------------------" << endl;
	cout << endl;

//*******************************************************************************
// plotting mass ratios
	TH1D* histoMeanMassRatio[nNL];
	for(Int_t i=0; i<nNL; i++){
		TString drawOption = (i==0)?"p":"p, same";
		histoMeanMassRatio[i] = (TH1D*)inputFiles[i]->Get("MeanMassRatioMCData-noFit");
		histoMeanMassRatio[i]->SetMarkerColor(color[i]);
		histoMeanMassRatio[i]->SetMarkerStyle(24+i);
		histoMeanMassRatio[i]->SetLineColor(color[i]);
		histoMeanMassRatio[i]->Draw(drawOption.Data());
		histoMeanMassRatio[i]->GetXaxis()->SetTitleOffset(1.1);
		histoMeanMassRatio[i]->GetXaxis()->SetTitle("ConvCalo: E_{Cluster} (GeV) - Calo: #it{p}^{#pi^{0}}_{T} (GeV/c)");
		legend->AddEntry(histoMeanMassRatio[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/MeanMassRatio_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

	TF1* histoMeanMassRatioFits[nNL];
	histoMeanMassRatio[0]->Draw("axis");
	for(Int_t i=0; i<nNL; i++){
		histoMeanMassRatioFits[i] = (TF1*)inputFiles[i]->Get("MeanMassRatioMCData-Fit");
		histoMeanMassRatioFits[i]->SetLineColor(color[i]);
		histoMeanMassRatioFits[i]->GetXaxis()->SetTitle("ConvCalo: E_{Cluster} (GeV) - Calo: #it{p}^{#pi^{0}}_{T} (GeV/c)");
		histoMeanMassRatioFits[i]->Draw("same");
		legend->AddEntry(histoMeanMassRatioFits[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/MeanMassRatio_Fits_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

//*******************************************************************************
// plotting total correction
	TH1D* histoTotalCorrection[nNL];
	for(Int_t i=0; i<nNL; i++){
		TString drawOption = (i==0)?"p":"p, same";
		histoTotalCorrection[i] = (TH1D*)inputFiles[i]->Get("Total Correction");
		histoTotalCorrection[i]->SetMarkerColor(color[i]);
		histoTotalCorrection[i]->SetMarkerStyle(24+i);
		histoTotalCorrection[i]->SetMarkerSize(0.4);
		histoTotalCorrection[i]->SetLineColor(color[i]);
		histoTotalCorrection[i]->Draw(drawOption.Data());
		legend->AddEntry(histoTotalCorrection[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/TotalCorrection_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

	TH1D* histoTotalCorrectionMC[nNL];
	for(Int_t i=0; i<nNL; i++){
		TString drawOption = (i==0)?"p":"p, same";
		histoTotalCorrectionMC[i] = (TH1D*)inputFiles[i]->Get("TotalCorrectionMC");
		histoTotalCorrectionMC[i]->SetMarkerColor(color[i]);
		histoTotalCorrectionMC[i]->SetMarkerStyle(24+i);
		histoTotalCorrectionMC[i]->SetMarkerSize(0.4);
		histoTotalCorrectionMC[i]->SetLineColor(color[i]);
		histoTotalCorrectionMC[i]->Draw(drawOption.Data());
		legend->AddEntry(histoTotalCorrectionMC[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/TotalCorrectionMC_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

	cout << "-----------------------------------------------------" << endl;
	cout << "-----------------------------------------------------" << endl;

	fOutput->Write();
	fOutput->Close();

	delete canvas;

	return;
}
