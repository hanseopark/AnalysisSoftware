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
void CorrectCaloNonLinearity2_Compare(TString select = "LHC11a-Pythia")
{
	gROOT->Reset();

	StyleSettingsThesis();
	SetPlotStyle();

	TString inputDir = "CorrectCaloNonLinearity2";
	TString outputDir = Form("%s/Compare",inputDir.Data());
	gSystem->Exec("mkdir -p "+outputDir);

	TString suffix = "eps";

	Int_t nNL = 0;
	TString *inputFileNames = 0x0;
	TString *inputFilePaths = 0x0;
	Bool_t *plotMassData = 0x0;
	TString xTitle = "ConvCalo: E_{Cluster} (GeV) - Calo: #it{p}^{#pi^{0}}_{T} (GeV/c)";
	TFile **inputFiles = 0x0;

	const Int_t nColor = 13;
	const Int_t nMarkerStyle = 7;
	Color_t color[nColor] = {kBlack,633,807,/*800,*/418,/*kGreen+4,*/435,601,879,806,852,kCyan+3,426};
	Int_t markerStyle[nMarkerStyle] = {24,25,27,28,29,30,31};

//*******************************************************************************
// Choosing data set
	if(select.CompareTo("LHC11a-Pythia")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC11a-Pythia-ConvCalo";
		inputFileNames[1] = "LHC11a-Pythia-Calo";
	}else if(select.CompareTo("LHC11a-Phojet")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC11a-Phojet-ConvCalo";
		inputFileNames[1] = "LHC11a-Phojet-Calo";
	}else if(select.CompareTo("LHC11a-JJ")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC11a-JJ-ConvCalo";
		inputFileNames[1] = "LHC11a-JJ-Calo";
	}else if(select.CompareTo("LHC12-Pythia")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC12-Pythia-ConvCalo";
		inputFileNames[1] = "LHC12-Pythia-Calo";
	}else if(select.CompareTo("LHC12-Phojet")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC12-Phojet-ConvCalo";
		inputFileNames[1] = "LHC12-Phojet-Calo";
	}else if(select.CompareTo("LHC13bc")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		inputFileNames[0] = "LHC13bc-ConvCalo";
		inputFileNames[1] = "LHC13bc-Calo";
	}else if(select.CompareTo("ConvCalo")==0){
		nNL = 6;
		inputFileNames = new TString[nNL];
		plotMassData = new Bool_t[nNL];
		inputFileNames[0] = "LHC11a-JJ-ConvCalo";		plotMassData[0] = kTRUE;
		inputFileNames[1] = "LHC11a-Pythia-ConvCalo";	plotMassData[1] = kFALSE;
		inputFileNames[2] = "LHC11a-Phojet-ConvCalo";	plotMassData[2] = kFALSE;
		inputFileNames[3] = "LHC12-Pythia-ConvCalo";	plotMassData[3] = kTRUE;
		inputFileNames[4] = "LHC12-Phojet-ConvCalo";	plotMassData[4] = kFALSE;
		inputFileNames[5] = "LHC13bc-ConvCalo";			plotMassData[5] = kTRUE;
		xTitle = "E_{Cluster} (GeV)";
	}else if(select.CompareTo("Calo")==0){
		nNL = 6;
		inputFileNames = new TString[nNL];
		plotMassData = new Bool_t[nNL];
		inputFileNames[0] = "LHC11a-JJ-Calo";		plotMassData[0] = kTRUE;
		inputFileNames[1] = "LHC11a-Pythia-Calo";	plotMassData[1] = kFALSE;
		inputFileNames[2] = "LHC11a-Phojet-Calo";	plotMassData[2] = kFALSE;
		inputFileNames[3] = "LHC12-Pythia-Calo";	plotMassData[3] = kTRUE;
		inputFileNames[4] = "LHC12-Phojet-Calo";	plotMassData[4] = kFALSE;
		inputFileNames[5] = "LHC13bc-Calo";			plotMassData[5] = kTRUE;
		xTitle = "#it{p}^{#pi^{0}}_{T} (GeV/c)";
	}else if(select.CompareTo("LHC11a-ConvCalo")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		plotMassData = new Bool_t[nNL];
		inputFileNames[0] = "LHC11a-Pythia-ConvCalo";	plotMassData[0] = kTRUE;
		inputFileNames[1] = "LHC11a-Phojet-ConvCalo";	plotMassData[1] = kFALSE;
		xTitle = "E_{Cluster} (GeV)";
	}else if(select.CompareTo("LHC11a-Calo")==0){
		nNL = 2;
		inputFileNames = new TString[nNL];
		plotMassData = new Bool_t[nNL];
		inputFileNames[0] = "LHC11a-Pythia-Calo";	plotMassData[0] = kTRUE;
		inputFileNames[1] = "LHC11a-Phojet-Calo";	plotMassData[1] = kFALSE;
		xTitle = "#it{p}^{#pi^{0}}_{T} (GeV/c)";
	}else if(select.CompareTo("LHC11a-ConvCalo-Calo")==0){
		nNL = 4;
		inputFileNames = new TString[nNL];
		plotMassData = new Bool_t[nNL];
		inputFileNames[0] = "LHC11a-Pythia-ConvCalo";	plotMassData[0] = kTRUE;
		inputFileNames[1] = "LHC11a-Phojet-ConvCalo";	plotMassData[1] = kFALSE;
		inputFileNames[2] = "LHC11a-Pythia-Calo";	plotMassData[2] = kFALSE;
		inputFileNames[3] = "LHC11a-Phojet-Calo";	plotMassData[3] = kFALSE;
	}else{
		cout << "No valid selection '" << select.Data() << "'' given, returning..." << endl;
		return;
	}

	if(nNL>nColor||nNL>nMarkerStyle){
		cout << "nNL: " << nNL << ", but defined nColor " << nColor << " and defined nMarkerStyle " << nMarkerStyle << endl;
		return;
	}

//*******************************************************************************
// Input
	inputFilePaths = new TString[nNL];
	inputFiles = new TFile*[nNL];
	for(Int_t i=0; i<nNL; i++){
		inputFilePaths[i] = Form("%s/CorrectCaloNonLinearity2_%s.root",inputDir.Data(),inputFileNames[i].Data());
		inputFiles[i] = new TFile(inputFilePaths[i].Data(),"READ");
			if(inputFiles[i]->IsZombie()) {cout << "Info: ROOT file '" << inputFilePaths[i].Data() << "' could not be openend, return!" << endl; return;}
	}

//*******************************************************************************
// Output
//	TString nameOutput = "";
//	nameOutput = Form("%s/CorrectCaloNonLinearity2_%s.root",outputDir.Data(),select.Data());
//	TFile* fOutput = new TFile(nameOutput,"RECREATE");

	TCanvas *canvas = new TCanvas("canvas","",200,10,1350,900);  // gives the page size
	Double_t leftMargin = 0.1; Double_t rightMargin = 0.02; Double_t topMargin = 0.06; Double_t bottomMargin = 0.1;
	DrawGammaCanvasSettings(canvas, leftMargin, rightMargin, topMargin, bottomMargin);

	TLegend *legend = new TLegend(0.2,0.75,0.4,0.9);
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
		histoMassMC[i]->SetMarkerStyle(markerStyle[i]);
		histoMassMC[i]->SetLineColor(color[i]);
		histoMassMC[i]->Draw(drawOption.Data());
		histoMassMC[i]->GetXaxis()->SetTitleOffset(1.1);
		histoMassMC[i]->GetXaxis()->SetTitle(xTitle.Data());
		legend->AddEntry(histoMassMC[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/Mass_MC_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

	TH1D* histoMassData[nNL];
	for(Int_t i=0; i<nNL; i++){
		TString plotNameData = inputFileNames[i];
		if(plotMassData && !plotMassData[i]) continue;
		else{
			TObjArray *d = plotNameData.Tokenize("-");
			TObjString* rString = (TObjString*)d->At(0);
			plotNameData = rString->GetString();
			for(Int_t iObjArr=1; iObjArr<d->GetEntries(); iObjArr++){
				TObjString* tmpObjString = (TObjString*)d->At(iObjArr);
				TString tmpString = tmpObjString->GetString();
				if(tmpString.Contains("Calo")){
					plotNameData+="-";
					plotNameData+=tmpString;
				}
			}
		}
		TString drawOption = (i==0)?"p":"p, same";
		histoMassData[i] = (TH1D*)inputFiles[i]->Get("Mean mass Data");
		histoMassData[i]->SetMarkerColor(color[i]);
		histoMassData[i]->SetMarkerStyle(markerStyle[i]);
		histoMassData[i]->SetLineColor(color[i]);
		histoMassData[i]->Draw(drawOption.Data());
		histoMassData[i]->GetXaxis()->SetTitleOffset(1.1);
		histoMassData[i]->GetXaxis()->SetTitle(xTitle.Data());
		legend->AddEntry(histoMassData[i],plotNameData.Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/Mass_Data_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

//*******************************************************************************
// plotting mass ratios
	TH1D* histoMeanMassRatio[nNL];
	for(Int_t i=0; i<nNL; i++){
		TString drawOption = (i==0)?"p":"p, same";
		histoMeanMassRatio[i] = (TH1D*)inputFiles[i]->Get("MeanMassRatioMCData-noFit");
		histoMeanMassRatio[i]->SetTitle(" ");
		histoMeanMassRatio[i]->SetMarkerColor(color[i]);
		histoMeanMassRatio[i]->SetMarkerStyle(markerStyle[i]);
		histoMeanMassRatio[i]->SetLineColor(color[i]);
		histoMeanMassRatio[i]->Draw(drawOption.Data());
		histoMeanMassRatio[i]->GetXaxis()->SetTitleOffset(1.1);
		histoMeanMassRatio[i]->GetXaxis()->SetTitle(xTitle.Data());
		legend->AddEntry(histoMeanMassRatio[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/MeanMassRatio_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

	TF1* histoMeanMassRatioFits[nNL];
	histoMeanMassRatio[0]->Draw("axis");
	for(Int_t i=0; i<nNL; i++){
		histoMeanMassRatioFits[i] = (TF1*)inputFiles[i]->Get("MeanMassRatioMCData-Fit");
		histoMeanMassRatioFits[i]->SetTitle(" ");
		histoMeanMassRatioFits[i]->SetLineColor(color[i]);
		histoMeanMassRatioFits[i]->GetXaxis()->SetTitle(xTitle.Data());
		histoMeanMassRatioFits[i]->Draw("same");
		legend->AddEntry(histoMeanMassRatioFits[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/MeanMassRatio_Fits_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();


	for(Int_t i=0; i<nNL; i++){
		TString drawOption = (i==0)?"p":"p, same";
		histoMeanMassRatio[i]->Draw(drawOption.Data());
		histoMeanMassRatioFits[i]->Draw("same");
		legend->AddEntry(histoMeanMassRatio[i],inputFileNames[i].Data());
	}

	legend->Draw("same");
	canvas->SetLogx(1); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
	canvas->SaveAs(Form("%s/MeanMassRatioAndFits_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
	canvas->Clear();
	legend->Clear();

//*******************************************************************************
// plotting total correction
	TH1D* histoTotalCorrection[nNL];
	for(Int_t i=0; i<nNL; i++){
		TString drawOption = (i==0)?"p":"p, same";
		histoTotalCorrection[i] = (TH1D*)inputFiles[i]->Get("Total Correction");
		histoTotalCorrection[i]->SetMarkerColor(color[i]);
		histoTotalCorrection[i]->SetMarkerStyle(markerStyle[i]);
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

	cout << "-----------------------------------------------------" << endl;
	cout << "-----------------------------------------------------" << endl;

//	fOutput->Write();
//	fOutput->Close();

	delete canvas;

	return;
}
