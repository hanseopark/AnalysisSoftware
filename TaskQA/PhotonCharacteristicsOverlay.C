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
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TMarker.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"

void  PhotonCharacteristicsOverlay(TString nameFileInputPythia6 = "myOutput", TString nameFileInputPythia8 = "myOutput",TString nameFileInputPhojet = "myOutput", TString cutSelection = "0001"){

	
	gROOT->Reset();
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();
	SetPlotStyle();
	
	
	TString etaCutNumber = cutSelection(1,1);
	cout << "eta cutnumber = " << etaCutNumber.Data() << endl;
	TLatex *latexEtaRange;
	
	Double_t maxEta = 0.;
	TString nameNTrackHisto = "ESD_NumberOfGoodESDTracks09Vtx";
 	if(etaCutNumber.CompareTo("0") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"|#eta| < 0.9 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);      
			nameNTrackHisto = "ESD_NumberOfGoodESDTracks09Vtx";
	} else if(etaCutNumber.CompareTo("2") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"0.9 < |#eta| < 1.4 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);		
			nameNTrackHisto = "ESD_NumberOfGoodESDTracks0914Vtx";
	} else if(etaCutNumber.CompareTo("4") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92," |#eta| < 1.4 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);
			nameNTrackHisto = "ESD_NumberOfGoodESDTracks14Vtx";
	} else if(etaCutNumber.CompareTo("6") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92," |#eta| < 10.0 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);      
			nameNTrackHisto = "ESD_NumberOfGoodESDTracks09Vtx";
	}    

	
	// the outputfile
	TFile* fileDataPythia6 = new TFile(nameFileInputPythia6); 
	TFile* fileDataPythia8 = new TFile(nameFileInputPythia8); 
	TFile* fileDataPhojet = new TFile(nameFileInputPhojet); 
	
	TH1D *histoInputSpecGammaPtPhojet;
	TH1D *histoInputSpecGammaPtPythia6;
	TH1D *histoInputSpecGammaPtPythia8;
	TH1D *histoInputSpecGammaPtRatioPhojetPythia6;
	TH1D *histoInputSpecGammaPtRatioPythia8Pythia6;
	TH1D *histoInputSpecGammaPtRatioPythia6Pythia6;
		
	TDirectory* dirPhojet = (TDirectory*)fileDataPhojet->Get(Form("GammaConv_%s",cutSelection.Data())); 
	TH1F* 	histoEventQualityPhojet=	(TH1F*)dirPhojet->Get("ESD_EventQuality");
	TH1F* 	histoNumberOfGoodESDTracksPhojet=	(TH1F*)dirPhojet->Get(nameNTrackHisto.Data());                       
	Float_t numberGoodEventsPhojet = 					histoEventQualityPhojet->GetBinContent(1);
	Float_t numberGoodTriggerPhojet = 					histoEventQualityPhojet->GetEntries();
	Double_t meanMultiplitcityPhojet = 					histoNumberOfGoodESDTracksPhojet->GetMean();
	Float_t normFactorReconstPhojet=						1./numberGoodEventsPhojet *1/1.6* 1./meanMultiplitcityPhojet;
	histoInputSpecGammaPtPhojet= (TH1D*)dirPhojet->Get("MC_ConvQualAftCuts_Pt");        
	cout << "Phojet :" << numberGoodEventsPhojet << "\t" << meanMultiplitcityPhojet << "\t" << histoInputSpecGammaPtPhojet->GetEntries()<<endl;
	
	GammaScalingHistogramm(histoInputSpecGammaPtPhojet,normFactorReconstPhojet);
	

	TDirectory* dirPythia6 = (TDirectory*)fileDataPythia6->Get(Form("GammaConv_%s",cutSelection.Data())); 
	TH1F* 	histoEventQualityPythia6=	(TH1F*)dirPythia6->Get("ESD_EventQuality");
	TH1F* 	histoNumberOfGoodESDTracksPythia6=	(TH1F*)dirPythia6->Get(nameNTrackHisto.Data());                       
	Float_t numberGoodEventsPythia6 = 					histoEventQualityPythia6->GetBinContent(1);
	Float_t numberGoodTriggerPythia6 = 					histoEventQualityPythia6->GetEntries();
	Double_t meanMultiplitcityPythia6 = 					histoNumberOfGoodESDTracksPythia6->GetMean();
	Float_t normFactorReconstPythia6=						1./numberGoodEventsPythia6 *1/1.6* 1./meanMultiplitcityPythia6;
	histoInputSpecGammaPtPythia6=	(TH1D*)dirPythia6->Get("MC_ConvQualAftCuts_Pt");        
	cout << "Pythia6 :" << numberGoodEventsPythia6 << "\t" << meanMultiplitcityPythia6 << "\t" << histoInputSpecGammaPtPythia6->GetEntries()<<endl;

	GammaScalingHistogramm(histoInputSpecGammaPtPythia6,normFactorReconstPythia6);
	
	TDirectory* dirPythia8 = (TDirectory*)fileDataPythia8->Get(Form("GammaConv_%s",cutSelection.Data())); 
	TH1F* 	histoEventQualityPythia8=	(TH1F*)dirPythia8->Get("ESD_EventQuality");
	TH1F* 	histoNumberOfGoodESDTracksPythia8=	(TH1F*)dirPythia8->Get(nameNTrackHisto.Data());                       
	Float_t numberGoodEventsPythia8 = 					histoEventQualityPythia8->GetBinContent(1);
	Float_t numberGoodTriggerPythia8 = 					histoEventQualityPythia8->GetEntries();
	Double_t meanMultiplitcityPythia8 = 					histoNumberOfGoodESDTracksPythia8->GetMean();
	Float_t normFactorReconstPythia8=						1./numberGoodEventsPythia8 *1/1.6* 1./meanMultiplitcityPythia8;
	histoInputSpecGammaPtPythia8=	(TH1D*)dirPythia8->Get("MC_ConvQualAftCuts_Pt");        
	cout << "Pythia8 :" << numberGoodEventsPythia8 << "\t" << meanMultiplitcityPythia8 << "\t" << histoInputSpecGammaPtPythia8->GetEntries()<<endl;
	
	GammaScalingHistogramm(histoInputSpecGammaPtPythia8,normFactorReconstPythia8);
	
	cout << "here" << endl;
	histoInputSpecGammaPtRatioPhojetPythia6 = (TH1D*)histoInputSpecGammaPtPhojet->Clone();
	histoInputSpecGammaPtRatioPhojetPythia6->Sumw2();
	cout << "here" << endl;
	histoInputSpecGammaPtRatioPhojetPythia6->Divide(histoInputSpecGammaPtPhojet,histoInputSpecGammaPtPythia6,1,1,"B");
	histoInputSpecGammaPtRatioPythia8Pythia6 = (TH1D*)histoInputSpecGammaPtPythia8->Clone();
	histoInputSpecGammaPtRatioPythia8Pythia6->Sumw2();
	cout << "here" << endl;
	histoInputSpecGammaPtRatioPythia8Pythia6->Divide(histoInputSpecGammaPtPythia8,histoInputSpecGammaPtPythia6,1,1,"B");
	histoInputSpecGammaPtRatioPythia6Pythia6 = (TH1D*)histoInputSpecGammaPtPythia6->Clone();
	histoInputSpecGammaPtRatioPythia6Pythia6->Sumw2();
	cout << "here" << endl;
	histoInputSpecGammaPtRatioPythia6Pythia6->Divide(histoInputSpecGammaPtPythia6,histoInputSpecGammaPtPythia6,1,1,"B");
	
	
	TCanvas* canvasInputSpectrum = new TCanvas("canvasInputSpectrum","",1350,1500);  // gives the page size
	DrawGammaCanvasSettings( canvasInputSpectrum, 0.15, 0.02, 0.02, 0.09);
	
	TPad* padInputSpecHistos = new TPad("padInputSpecHistos", "", 0., 0.33, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padInputSpecHistos, 0.15, 0.02, 0.02, 0.00);
	padInputSpecHistos->Draw();
	
	TPad* padInputRatio = new TPad("padInputRatio", "", 0., 0., 1., 0.33,-1, -1, -2);
	DrawGammaPadSettings( padInputRatio, 0.15, 0.02, 0.00, 0.2);
	padInputRatio->Draw();
	
	padInputSpecHistos->cd();
	padInputSpecHistos->SetLogy();		
	padInputSpecHistos->SetLogx();		

	DrawGammaSetMarker(histoInputSpecGammaPtPhojet, 20,1., kBlue+2, kBlue+2);	
	DrawAutoGammaMesonHistos( histoInputSpecGammaPtPhojet,
				"", "p_{t} (GeV/c)", "#frac{1}{N_{ch}} #frac{d^{2}N_{#gamma}}{p_{t}dp_{t}dy} P (c/GeV)^{2}",
				kTRUE, 3., 6e-7,kFALSE,
				kFALSE, 1e-10, 2.,
				kTRUE, 0.1, 20.);
	DrawGammaSetMarker(histoInputSpecGammaPtPythia6, 20,1., kGreen+2, kGreen+2);	
	histoInputSpecGammaPtPythia6->Draw("same,e1,p");
	DrawGammaSetMarker(histoInputSpecGammaPtPythia8, 20,1., kRed+2, kRed+2);	
	histoInputSpecGammaPtPythia8->Draw("same,e1,p");

	
	TLegend* legendInputSpec;
	legendInputSpec = new TLegend( 0.55,0.84,0.97,0.97);
	legendInputSpec->SetTextSize(0.03);                        
	legendInputSpec->SetFillColor(0);
	legendInputSpec->AddEntry(histoInputSpecGammaPtPhojet,("Phojet 2.76 TeV"),"p");
	legendInputSpec->AddEntry(histoInputSpecGammaPtPythia6,("Pythia 6 2.76 TeV"),"p");
	legendInputSpec->AddEntry(histoInputSpecGammaPtPythia8,("Pythia 8 2.76 TeV"),"p");
	legendInputSpec->Draw();

	padInputRatio->cd();
	padInputRatio->SetLogx(1);
// 	padInputRatio->SetLogy(1);		

	histoInputSpecGammaPtRatioPhojetPythia6->SetYTitle("#frac{Generator}{Pythia6}");	
	histoInputSpecGammaPtRatioPhojetPythia6->SetXTitle("p_{T} (GeV/c)");	
	histoInputSpecGammaPtRatioPhojetPythia6->GetYaxis()->SetRangeUser(0,2.1);
	histoInputSpecGammaPtRatioPhojetPythia6->GetXaxis()->SetRangeUser(1e-1,20.);
	histoInputSpecGammaPtRatioPhojetPythia6->GetYaxis()->SetNdivisions(505);
	histoInputSpecGammaPtRatioPhojetPythia6->GetYaxis()->SetLabelSize(0.07);
	histoInputSpecGammaPtRatioPhojetPythia6->GetYaxis()->SetTitleSize(0.09);
	histoInputSpecGammaPtRatioPhojetPythia6->GetYaxis()->SetDecimals();
	histoInputSpecGammaPtRatioPhojetPythia6->GetYaxis()->SetTitleOffset(0.5);
	histoInputSpecGammaPtRatioPhojetPythia6->GetXaxis()->SetTitleSize(0.09);
	histoInputSpecGammaPtRatioPhojetPythia6->GetXaxis()->SetTitleOffset(1);
	histoInputSpecGammaPtRatioPhojetPythia6->GetXaxis()->SetLabelSize(0.075);

	DrawGammaSetMarker(histoInputSpecGammaPtRatioPhojetPythia6,20, 1., kBlue+2, kBlue+2);
	histoInputSpecGammaPtRatioPhojetPythia6->DrawCopy("p,e1"); 	
	DrawGammaSetMarker(histoInputSpecGammaPtRatioPythia8Pythia6,20, 1., kRed+2, kRed+2);
	histoInputSpecGammaPtRatioPythia8Pythia6->DrawCopy("same,p,e1"); 	
	DrawGammaSetMarker(histoInputSpecGammaPtRatioPythia6Pythia6,20, 1., kGreen+2, kGreen+2);
	histoInputSpecGammaPtRatioPythia6Pythia6->DrawCopy("same,p,e1"); 	

	canvasInputSpectrum->Update();

	canvasInputSpectrum->SaveAs(Form("InputSpecGammaPt_2.76TeV_%s.eps",cutSelection.Data()));
	delete legendInputSpec;
	delete padInputRatio;
	delete padInputSpecHistos;
	delete canvasInputSpectrum;

}        
