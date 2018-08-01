// provided by Gamma Conversion Group, PWGGA, Lucia Leardini lucia.leardini@cern.ch
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
#include "TString.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TH1D.h"
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
#include "TMarker.h"
#include "TGraphAsymmErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"



void PPRstyle()
{
//   gStyle->SetPalette(1);
//   gStyle->SetCanvasBorderMode(-1);
//   gStyle->SetCanvasBorderSize(1);
//   gStyle->SetCanvasColor(10);
//
//   gStyle->SetFrameFillColor(10);
//   gStyle->SetFrameBorderSize(1);
//   gStyle->SetFrameBorderMode(-1);
//   gStyle->SetFrameLineWidth(1);
//   gStyle->SetFrameLineColor(1);
//
//   gStyle->SetHistFillColor(0);
//   gStyle->SetHistLineWidth(0.5);
//   gStyle->SetHistLineColor(1);

  gStyle->SetPadColor(10);
  gStyle->SetPadTickY(1);  	// TICKS ON BOTH SIDES OF THE PAD
  gStyle->SetPadTickX(1);	// TICKS ON BOTH SIDES OF THE PAD
//   gStyle->SetPadBorderSize(1);
//   gStyle->SetPadBorderMode(-1);

//   gStyle->SetStatColor(10);
//   gStyle->SetTitleColor(kBlack,"X");
//   gStyle->SetTitleColor(kBlack,"Y");

  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetLabelSize(0.04,"Z");
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleSize(0.04,"Z");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelFont(42,"Z");
  gStyle->SetStatFont(42);

  gStyle->SetTitleOffset(1.9,"X");
  gStyle->SetTitleOffset(1.4,"Y");

  gStyle->SetFillColor(kWhite);
  gStyle->SetTitleFillColor(kWhite);

  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

}


void DrawHisto(TH1* histo1,TH1* histo2,TString Output,TString Name,TString Title, TString XTitle, TString YTitle,Float_t xMin, Float_t xMax, Float_t yMin, Float_t yMax, Bool_t line = 0) {
	PPRstyle();
	TH1D* histodummy= new TH1D("","",100,xMin,xMax);
	histodummy->SetAxisRange(yMin,yMax,"Y");

	if(XTitle.Length() > 0){
		histodummy->SetXTitle(XTitle.Data());
	}
	if(YTitle.Length() > 0){
		histodummy->SetYTitle(YTitle.Data());
	}
	TLatex *alice;
	if(Title.Length() > 0){
		alice = new TLatex(0.4,0.95,Form("%s",Title.Data()));
		alice->SetTextSize(0.055);
		alice->SetNDC();
		alice->SetTextColor(1);

	}


	Double_t markersize=1.;

// 	histodummy->GetYaxis()->SetDecimals();
// 	histodummy->GetYaxis()->SetLabelOffset(0.005);
	histodummy->SetTitleOffset(1.3,"xy");
	histodummy->SetTitleSize(0.04,"xy");
	histodummy->GetYaxis()->SetLabelSize(0.035);
	histodummy->GetXaxis()->SetLabelSize(0.035);
// 	histodummy->GetXaxis()->SetNdivisions(507,kTRUE);
// 	histodummy->GetYaxis()->SetTickLength(0.02);
// 	histodummy->GetXaxis()->SetTickLength(0.02);

	histo1->SetMarkerStyle(24);
	histo2->SetMarkerStyle(24);
	histo2->SetMarkerColor(kBlue);
	histo2->SetLineColor(kBlue);
	histo1->SetMarkerColor(kRed);
	histo1->SetLineColor(kRed);
// 	histo->SetLineWidth(1);
	histo1->SetMarkerSize(markersize);
	histo2->SetMarkerSize(markersize);

	TCanvas * c = new TCanvas("c","",2800,1800);  // gives the page size
	c->SetTopMargin(0.00);
	c->SetBottomMargin(0.00);
	c->SetRightMargin(0.00);
	c->SetLeftMargin(0.00);
	if(Name.Contains("histoRatioRecConvGamma_Pt_OnFly")) c->SetLogy();
	c->cd();

	TPad *p = new TPad("p","",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	p->SetFillColor(0);
	p->GetFrame()->SetFillColor(0);
	p->SetLeftMargin(0.1);
	p->SetRightMargin(0.02);
	p->SetTopMargin(0.08);
	p->SetBottomMargin(0.1);
	//p->SetBorderMode(0);
	if(Name.Contains("histoRatioRecConvGamma_Pt_OnFly")) p->SetLogy();
	//p->Divide(fColumn,fRow,0.0,0.0);
	p->Draw();
	p->cd();


	histodummy->DrawCopy("hist");

	if(line){
		histo1->DrawCopy("hist,same");
		histo2->DrawCopy("hist,same");
	}else{
		histo1->DrawCopy("same,x0");
		histo2->DrawCopy("same,x0");
	}
	if(Title.Length() > 0) alice->Draw();
	TString name = histo1->GetName();

	if(name.Contains("mid")){

		TLegend* leg = new TLegend(0.8,0.15,0.9,0.25);
		leg->SetTextSize(0.04);
		leg->SetLineColor(0);
		leg->SetFillColor(0);
		leg->SetFillStyle(0);
		leg->AddEntry(histo1,("On-the-Fly"));
		leg->AddEntry(histo2,("Offline"));
		leg->Draw();

		TLatex *labelEnergy = new TLatex(0.15,0.18,"pp, #sqrt{s_{NN}} = 5 TeV (2017)");
		SetStyleTLatex( labelEnergy, 0.04,4);
		labelEnergy->Draw();

	} else {

		TLegend* leg = new TLegend(0.8,0.76,0.9,0.86);
		leg->SetTextSize(0.04);
		leg->SetLineColor(0);
		leg->SetFillColor(0);
		leg->SetFillStyle(0);
		leg->AddEntry(histo1,("On-the-Fly"));
		leg->AddEntry(histo2,("Offline"));
		leg->Draw();

		TLatex *labelEnergy = new TLatex(0.15,0.82,"pp, #sqrt{s_{NN}} = 5 TeV (2017)");
		SetStyleTLatex( labelEnergy, 0.04,4);
		labelEnergy->Draw();
	}

	c->SaveAs(Form("%s/%s.pdf",Output.Data(),Name.Data()));

	delete p;
	delete c;



}


void ProduceV0FindingEffPlotting(TString filename_OnFly = "",
                                 TString filename_Offline = "",
                                 TString nameFolder ="",
                                 TString cutString_OnFly = "",
                                 TString cutString_Offline = "",
                                 TString outputDir = "V0EfficiencyStudies"
){

    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetEndErrorSize(0);
    StyleSettingsThesis();
    SetPlotStyle();

    gSystem->Exec(Form("mkdir -p %s",outputDir.Data()));

	TFile *file_OnFly              = new TFile(filename_OnFly.Data());
	TList *list_OnFly              = (TList*) file_OnFly->Get(nameFolder.Data());
	TList *HistosV0Finding_OnFly   = (TList*)list_OnFly->FindObject(Form("V0FindingEfficiencyInput_%s", cutString_OnFly.Data()));

	//loading histos
	TH2F* histoConvGamma_Pt_R_OnFly        = (TH2F*)HistosV0Finding_OnFly->FindObject("MCconvGamma_Pt_R");
	TH2F* histoConvGamma_Pt_Eta_OnFly      = (TH2F*)HistosV0Finding_OnFly->FindObject("MCconvGamma_Pt_Eta");
	TH2F* histoConvGamma_Pt_Phi_OnFly      = (TH2F*)HistosV0Finding_OnFly->FindObject("MCconvGamma_Pt_Phi");

	TH2F* histoRecConvGamma_Pt_R_OnFly 	   = (TH2F*)HistosV0Finding_OnFly->FindObject("RecMCconvGamma_Pt_R");
	TH2F* histoRecConvGamma_Pt_Eta_OnFly   = (TH2F*)HistosV0Finding_OnFly->FindObject("RecMCconvGamma_Pt_Eta");
	TH2F* histoRecConvGamma_Pt_Phi_OnFly   = (TH2F*)HistosV0Finding_OnFly->FindObject("RecMCconvGamma_Pt_Phi");

	TH1D* histoRecConvGammaMult_Pt_OnFly   = (TH1D*)HistosV0Finding_OnFly->FindObject("RecMCconvGammaMulti_Pt");
	TH1D* histoRecConvGammaMult_R_OnFly    = (TH1D*)HistosV0Finding_OnFly->FindObject("RecMCconvGammaMulti_R");
	TH1D* histoRecConvGammaMult_Phi_OnFly  = (TH1D*)HistosV0Finding_OnFly->FindObject("RecMCconvGammaMulti_Phi");

	//projecting histos
	TH1D* histoConvGamma_Pt_OnFly 			= (TH1D*)histoConvGamma_Pt_R_OnFly->ProjectionX("histoConvGamma_Pt_OnFly"); //restricted to eta cut range
	TH1D* histoConvGamma_PtlargerEta_OnFly	= (TH1D*)histoConvGamma_Pt_Eta_OnFly->ProjectionX("histoConvGamma_PtlargerEta_OnFly"); //full eta range

	TH1D* histoConvGamma_R_OnFly			= (TH1D*)histoConvGamma_Pt_R_OnFly->ProjectionY("histoConvGamma_R_OnFly"); //full pT range
	TH1D* histoConvGamma_R_lowpT_OnFly 		= (TH1D*)histoConvGamma_Pt_R_OnFly->ProjectionY("histoConvGamma_R_lowpT_OnFly",
                                                                                            histoConvGamma_Pt_R_OnFly->GetXaxis()->FindBin(0.),
                                                                                            histoConvGamma_Pt_R_OnFly->GetXaxis()->FindBin(1.0)); // 0. < pT < 1. GeV/c
	TH1D* histoConvGamma_R_midpT_OnFly 		= (TH1D*)histoConvGamma_Pt_R_OnFly->ProjectionY("histoConvGamma_R_midpT_OnFly",
                                                                                            histoConvGamma_Pt_R_OnFly->GetXaxis()->FindBin(1.0),
                                                                                            histoConvGamma_Pt_R_OnFly->GetXaxis()->FindBin(4.0)); // 1.0 < pT < 4. GeV/c
	TH1D* histoConvGamma_R_highpT_OnFly		= (TH1D*)histoConvGamma_Pt_R_OnFly->ProjectionY("histoConvGamma_R_highpT_OnFly",
                                                                                            histoConvGamma_Pt_R_OnFly->GetXaxis()->FindBin(4.0),
                                                                                            histoConvGamma_Pt_R_OnFly->GetNbinsX()); // pT > 4. GeV/c

	TH1D* histoConvGamma_Eta_OnFly			= (TH1D*)histoConvGamma_Pt_Eta_OnFly->ProjectionY("histoConvGamma_Eta_OnFly");
	TH1D* histoConvGamma_Phi_OnFly 			= (TH1D*)histoConvGamma_Pt_Phi_OnFly->ProjectionY("histoConvGamma_Phi_OnFly");

	TH1D* histoRecConvGamma_Pt_OnFly 			= (TH1D*)histoRecConvGamma_Pt_R_OnFly->ProjectionX("histoRecConvGamma_Pt_OnFly"); //restricted to eta cut range
	TH1D* histoRecConvGamma_PtlargerEta_OnFly 	= (TH1D*)histoRecConvGamma_Pt_Eta_OnFly->ProjectionX("histoRecConvGamma_PtlargerEta_OnFly"); //full eta range
	TH1D* histoRecConvGamma_R_OnFly 			= (TH1D*)histoRecConvGamma_Pt_R_OnFly->ProjectionY("histoRecConvGamma_R_OnFly"); //full pT range
	TH1D* histoRecConvGamma_R_lowpT_OnFly 		= (TH1D*)histoRecConvGamma_Pt_R_OnFly->ProjectionY("histoRecConvGamma_R_lowpT_OnFly",
                                                                                                   histoRecConvGamma_Pt_R_OnFly->GetXaxis()->FindBin(0.),
                                                                                                   histoRecConvGamma_Pt_R_OnFly->GetXaxis()->FindBin(1.0)); // 0. < pT < 1. GeV/c
	TH1D* histoRecConvGamma_R_midpT_OnFly 		= (TH1D*)histoRecConvGamma_Pt_R_OnFly->ProjectionY("histoRecConvGamma_R_midpT_OnFly",
                                                                                                   histoRecConvGamma_Pt_R_OnFly->GetXaxis()->FindBin(1.0),
                                                                                                   histoRecConvGamma_Pt_R_OnFly->GetXaxis()->FindBin(4.0)); // 1.0 < pT < 4. GeV/c
	TH1D* histoRecConvGamma_R_highpT_OnFly		= (TH1D*)histoRecConvGamma_Pt_R_OnFly->ProjectionY("histoRecConvGamma_R_highpT_OnFly",
                                                                                                   histoRecConvGamma_Pt_R_OnFly->GetXaxis()->FindBin(4.0),
                                                                                                   histoRecConvGamma_Pt_R_OnFly->GetNbinsX()); // pT > 4. GeV/c

	TH1D* histoRecConvGamma_Eta_OnFly 		= (TH1D*)histoRecConvGamma_Pt_Eta_OnFly->ProjectionY("histoRecConvGamma_Eta_OnFly");
	TH1D* histoRecConvGamma_Phi_OnFly 		= (TH1D*)histoRecConvGamma_Pt_Phi_OnFly->ProjectionY("histoRecConvGamma_Phi_OnFly");

	//pt rebinning array
	Float_t ptcalcbinning1[5] = {0.1,0.5 ,1 ,2 ,5};
// 	Float_t ptcalcbinning2[4] = {14,40,100,200};
	Float_t ptcalcbinning3[6] = {1,11,17,19,21,23};
	Double_t ptbinning[23];
	ptbinning[0] = 0;
	for ( int j = 0; j < 5 ; j ++ ){
		for ( int i = ptcalcbinning3[j]; i<ptcalcbinning3[j+1]; i++){
			ptbinning[i] = ptbinning[i-1] +	 ptcalcbinning1[j];
		}
	}

		//R rebinning array
	Float_t Rcalcbinning1[18] =  {5, 0.5, 2, 0.5, 4,0.5,1,0.5,3.5,0.5,4,0.5,5,0.5,4,1,11,15};
	Float_t Rcalcbinning3[19] = {1,2,16,17,27,28,38,39,41,43,65,68,71,74,79,80,81,82,87};
	Double_t Rbinning[87];
	Rbinning[0] = 0;
	for ( int j = 0; j < 18 ; j ++ ){
		for ( int i = Rcalcbinning3[j]; i<Rcalcbinning3[j+1]; i++){
			Rbinning[i] = Rbinning[i-1] +	 Rcalcbinning1[j];
		}
	}

	//making ratios = efficiency
	TH1D* histoRecConvGamma_Pt_rebin_OnFly 	= (TH1D*)histoRecConvGamma_Pt_OnFly->Clone(); //Rebin(22,"recPteff",ptbinning);
	TH1D* histoConvGamma_Pt_rebin_OnFly		= (TH1D*)histoConvGamma_Pt_OnFly->Clone(); //->Rebin(22,"totPteff",ptbinning);
	TH1D* histoRatioRecConvGamma_Pt_OnFly 	= (TH1D*)histoRecConvGamma_Pt_rebin_OnFly->Clone("histoRatioRecConvGamma_Pt_OnFly");
	histoRatioRecConvGamma_Pt_OnFly->Sumw2();
	histoRatioRecConvGamma_Pt_OnFly->Divide(histoRecConvGamma_Pt_rebin_OnFly,histoConvGamma_Pt_rebin_OnFly,1.,1.,"B");

	TH1D* histoRecConvGamma_PtlargerEta_rebin_OnFly 	= (TH1D*)histoRecConvGamma_PtlargerEta_OnFly->Clone(); //->Rebin(22,"recPtlargeEta",ptbinning);
	TH1D* histoConvGamma_PtlargerEta_rebin_OnFly		= (TH1D*)histoConvGamma_PtlargerEta_OnFly->Clone(); //->Rebin(22,"totPtlargeEta",ptbinning);
	TH1D* histoRatioRecConvGamma_PtlargerEta_OnFly 		= (TH1D*)histoRecConvGamma_PtlargerEta_rebin_OnFly->Clone("histoRatioRecConvGamma_PtlargerEta_OnFly");
	histoRatioRecConvGamma_PtlargerEta_OnFly->Sumw2();
	histoRatioRecConvGamma_PtlargerEta_OnFly->Divide(histoRecConvGamma_PtlargerEta_rebin_OnFly,histoConvGamma_PtlargerEta_rebin_OnFly,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_OnFly = (TH1D*)histoRecConvGamma_R_OnFly->Clone("histoRatioRecConvGamma_R_OnFly");
	histoRatioRecConvGamma_R_OnFly->Sumw2();
	histoRatioRecConvGamma_R_OnFly->Divide(histoRatioRecConvGamma_R_OnFly,histoConvGamma_R_OnFly,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_lowpT_OnFly = (TH1D*)histoRecConvGamma_R_lowpT_OnFly->Clone("histoRatioRecConvGamma_R_lowpT_OnFly");
	histoRatioRecConvGamma_R_lowpT_OnFly->Sumw2();
	histoRatioRecConvGamma_R_lowpT_OnFly->Divide(histoRatioRecConvGamma_R_lowpT_OnFly,histoConvGamma_R_lowpT_OnFly,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_midpT_OnFly = (TH1D*)histoRecConvGamma_R_midpT_OnFly->Clone("histoRatioRecConvGamma_R_midpT_OnFly");
	histoRatioRecConvGamma_R_midpT_OnFly->Sumw2();
	histoRatioRecConvGamma_R_midpT_OnFly->Divide(histoRatioRecConvGamma_R_midpT_OnFly,histoConvGamma_R_midpT_OnFly,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_highpT_OnFly = (TH1D*)histoRecConvGamma_R_highpT_OnFly->Clone("histoRatioRecConvGamma_R_highpT_OnFly");
	histoRatioRecConvGamma_R_highpT_OnFly->Sumw2();
	histoRatioRecConvGamma_R_highpT_OnFly->Divide(histoRatioRecConvGamma_R_highpT_OnFly,histoConvGamma_R_highpT_OnFly,1.,1.,"B");



	TH1D* histoRatioRecConvGamma_Eta_OnFly = (TH1D*)histoRecConvGamma_Eta_OnFly->Clone("histoRatioRecConvGamma_Eta_OnFly");
	histoRatioRecConvGamma_Eta_OnFly->Sumw2();
	histoRatioRecConvGamma_Eta_OnFly->Divide(histoRatioRecConvGamma_Eta_OnFly,histoConvGamma_Eta_OnFly,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_Phi_OnFly = (TH1D*)histoRecConvGamma_Phi_OnFly->Clone("histoRatioRecConvGamma_Phi_OnFly");
	histoRatioRecConvGamma_Phi_OnFly->Sumw2();
	histoRatioRecConvGamma_Phi_OnFly->Divide(histoRatioRecConvGamma_Phi_OnFly,histoConvGamma_Phi_OnFly,1.,1.,"B");

	//calculating the tot number of photons = 1st photon counted + 2nd time the same photon has been counted
	TH1D* histoTotalRecConvGamma_Pt_OnFly = (TH1D*)histoRecConvGamma_Pt_OnFly->Clone("histoRatioRecConvGamma_Pt_OnFly");
	histoRecConvGamma_Pt_OnFly->Add(histoRecConvGammaMult_Pt_OnFly);
	TH1D* histoTotalRecConvGamma_R_OnFly = (TH1D*)histoRecConvGamma_R_OnFly->Clone("histoRatioRecConvGamma_R_OnFly");
	histoTotalRecConvGamma_R_OnFly->Add(histoRecConvGammaMult_R_OnFly);
	TH1D* histoTotalRecConvGamma_Phi_OnFly = (TH1D*)histoRecConvGamma_Phi_OnFly->Clone("histoRatioRecConvGamma_Phi_OnFly");
	histoTotalRecConvGamma_Phi_OnFly->Add(histoRecConvGammaMult_Phi_OnFly);



	//making ratio with the Mult histos = double counting ratio
	TH1D* histoRecConvGammaMult_Pt_rebin_OnFly  = (TH1D*)histoRecConvGammaMult_Pt_OnFly->Clone(); //->Rebin(22,"recPt",ptbinning);
	TH1D* histoTotalRecConvGamma_Pt_rebin_OnFly = (TH1D*)histoTotalRecConvGamma_Pt_OnFly->Clone(); //->Rebin(22,"totPt",ptbinning);
	TH1D* histoRatioRecConvGammaMult_Pt_OnFly   = (TH1D*)histoRecConvGammaMult_Pt_rebin_OnFly->Clone("histoRatioRecConvGammaMult_Pt_OnFly");
	histoRatioRecConvGammaMult_Pt_OnFly->Sumw2();
	histoRatioRecConvGammaMult_Pt_OnFly->Divide(histoRecConvGammaMult_Pt_rebin_OnFly,histoTotalRecConvGamma_Pt_rebin_OnFly,1.,1.,"B");

	TH1D* histoRecConvGammaMult_R_rebin_OnFly  = (TH1D*)histoRecConvGammaMult_R_OnFly->Clone(); //->Rebin(86,"recR",Rbinning);
	TH1D* histoTotalRecConvGamma_R_rebin_OnFly = (TH1D*)histoTotalRecConvGamma_R_OnFly->Clone(); //->Rebin(86,"recR",Rbinning);
	TH1D* histoRatioRecConvGammaMult_R_OnFly   = (TH1D*)histoRecConvGammaMult_R_rebin_OnFly->Clone("histoRatioRecConvGammaMult_R_OnFly");
	histoRatioRecConvGammaMult_R_OnFly->Sumw2();
	histoRatioRecConvGammaMult_R_OnFly->Divide(histoRecConvGammaMult_R_rebin_OnFly,histoTotalRecConvGamma_R_rebin_OnFly,1.,1.,"B");

	TH1D* histoRatioRecConvGammaMult_Phi_OnFly 	= (TH1D*)histoRecConvGammaMult_Phi_OnFly->Clone("histoRatioRecConvGammaMult_Phi_OnFly");
	histoRatioRecConvGammaMult_Phi_OnFly->Sumw2();
	histoRatioRecConvGammaMult_Phi_OnFly->Divide(histoRatioRecConvGammaMult_Phi_OnFly,histoTotalRecConvGamma_Phi_OnFly,1.,1.,"B");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	TFile *file_Offline = new TFile(filename_Offline.Data());
	TList *list_Offline = (TList*) file_Offline->Get(nameFolder.Data());
	TList *HistosV0Finding_Offline = (TList*)list_Offline->FindObject(Form("V0FindingEfficiencyInput_%s", cutString_Offline.Data()));

	//loading histos
	TH2F* histoConvGamma_Pt_R_Offline 	= (TH2F*)HistosV0Finding_Offline->FindObject("MCconvGamma_Pt_R");
	TH2F* histoConvGamma_Pt_Eta_Offline = (TH2F*)HistosV0Finding_Offline->FindObject("MCconvGamma_Pt_Eta");
	TH2F* histoConvGamma_Pt_Phi_Offline = (TH2F*)HistosV0Finding_Offline->FindObject("MCconvGamma_Pt_Phi");

	TH2F* histoRecConvGamma_Pt_R_Offline 	= (TH2F*)HistosV0Finding_Offline->FindObject("RecMCconvGamma_Pt_R");
	TH2F* histoRecConvGamma_Pt_Eta_Offline 	= (TH2F*)HistosV0Finding_Offline->FindObject("RecMCconvGamma_Pt_Eta");
	TH2F* histoRecConvGamma_Pt_Phi_Offline 	= (TH2F*)HistosV0Finding_Offline->FindObject("RecMCconvGamma_Pt_Phi");

	TH1D* histoRecConvGammaMult_Pt_Offline 	= (TH1D*)HistosV0Finding_Offline->FindObject("RecMCconvGammaMulti_Pt");
	TH1D* histoRecConvGammaMult_R_Offline 	= (TH1D*)HistosV0Finding_Offline->FindObject("RecMCconvGammaMulti_R");
	TH1D* histoRecConvGammaMult_Phi_Offline = (TH1D*)HistosV0Finding_Offline->FindObject("RecMCconvGammaMulti_Phi");

	//projecting histos
	TH1D* histoConvGamma_Pt_Offline 			= (TH1D*)histoConvGamma_Pt_R_Offline->ProjectionX("histoConvGamma_Pt_Offline"); //restricted to eta cut range
	TH1D* histoConvGamma_PtlargerEta_Offline	= (TH1D*)histoConvGamma_Pt_Eta_Offline->ProjectionX("histoConvGamma_PtlargerEta_Offline"); //full eta range

	TH1D* histoConvGamma_R_Offline				= (TH1D*)histoConvGamma_Pt_R_Offline->ProjectionY("histoConvGamma_R_Offline"); //full pT range
	TH1D* histoConvGamma_R_lowpT_Offline 		= (TH1D*)histoConvGamma_Pt_R_Offline->ProjectionY("histoConvGamma_R_lowpT_Offline",
                                                                                                  histoConvGamma_Pt_R_Offline->GetXaxis()->FindBin(0.),
                                                                                                  histoConvGamma_Pt_R_Offline->GetXaxis()->FindBin(1.0)); // 0. < pT < 1. GeV/c
	TH1D* histoConvGamma_R_midpT_Offline 		= (TH1D*)histoConvGamma_Pt_R_Offline->ProjectionY("histoConvGamma_R_midpT_Offline",
                                                                                                  histoConvGamma_Pt_R_Offline->GetXaxis()->FindBin(1.0),
                                                                                                  histoConvGamma_Pt_R_Offline->GetXaxis()->FindBin(4.0)); // 1.0 < pT < 4. GeV/c
	TH1D* histoConvGamma_R_highpT_Offline		= (TH1D*)histoConvGamma_Pt_R_Offline->ProjectionY("histoConvGamma_R_highpT_Offline",
                                                                                                  histoConvGamma_Pt_R_Offline->GetXaxis()->FindBin(4.0),
                                                                                                  histoConvGamma_Pt_R_Offline->GetNbinsX()); // pT > 4. GeV/c

	TH1D* histoConvGamma_Eta_Offline			= (TH1D*)histoConvGamma_Pt_Eta_Offline->ProjectionY("histoConvGamma_Eta_Offline");
	TH1D* histoConvGamma_Phi_Offline 			= (TH1D*)histoConvGamma_Pt_Phi_Offline->ProjectionY("histoConvGamma_Phi_Offline");

	TH1D* histoRecConvGamma_Pt_Offline 			= (TH1D*)histoRecConvGamma_Pt_R_Offline->ProjectionX("histoRecConvGamma_Pt_Offline"); //restricted to eta cut range
	TH1D* histoRecConvGamma_PtlargerEta_Offline = (TH1D*)histoRecConvGamma_Pt_Eta_Offline->ProjectionX("histoRecConvGamma_PtlargerEta_Offline"); //full eta range
	TH1D* histoRecConvGamma_R_Offline 			= (TH1D*)histoRecConvGamma_Pt_R_Offline->ProjectionY("histoRecConvGamma_R_Offline"); //full pT range
	TH1D* histoRecConvGamma_R_lowpT_Offline 	= (TH1D*)histoRecConvGamma_Pt_R_Offline->ProjectionY("histoRecConvGamma_R_lowpT_Offline",
                                                                                                         histoRecConvGamma_Pt_R_Offline->GetXaxis()->FindBin(0.),
                                                                                                         histoRecConvGamma_Pt_R_Offline->GetXaxis()->FindBin(1.0)); // 0. < pT < 1. GeV/c
	TH1D* histoRecConvGamma_R_midpT_Offline 	= (TH1D*)histoRecConvGamma_Pt_R_Offline->ProjectionY("histoRecConvGamma_R_midpT_Offline",
                                                                                                         histoRecConvGamma_Pt_R_Offline->GetXaxis()->FindBin(1.0),
                                                                                                         histoRecConvGamma_Pt_R_Offline->GetXaxis()->FindBin(4.0)); // 1.0 < pT < 4. GeV/c
	TH1D* histoRecConvGamma_R_highpT_Offline	= (TH1D*)histoRecConvGamma_Pt_R_Offline->ProjectionY("histoRecConvGamma_R_highpT_Offline",
                                                                                                         histoRecConvGamma_Pt_R_Offline->GetXaxis()->FindBin(4.0),
                                                                                                         histoRecConvGamma_Pt_R_Offline->GetNbinsX()); // pT > 4. GeV/c

	TH1D* histoRecConvGamma_Eta_Offline 		= (TH1D*)histoRecConvGamma_Pt_Eta_Offline->ProjectionY("histoRecConvGamma_Eta_Offline");
	TH1D* histoRecConvGamma_Phi_Offline 		= (TH1D*)histoRecConvGamma_Pt_Phi_Offline->ProjectionY("histoRecConvGamma_Phi_Offline");

// 	//pt rebinning array
// 	Float_t ptcalcbinning1[5] = {0.1,0.5 ,1 ,2 ,5};
// // 	Float_t ptcalcbinning2[4] = {14,40,100,200};
// 	Float_t ptcalcbinning3[6] = {1,11,17,19,21,23};
// 	Double_t ptbinning[23];
// 	ptbinning[0] = 0;
// 	for ( int j = 0; j < 5 ; j ++ ){
// 		for ( int i = ptcalcbinning3[j]; i<ptcalcbinning3[j+1]; i++){
// 			ptbinning[i] = ptbinning[i-1] +	 ptcalcbinning1[j];
// 		}
// 	}
//
// 		//R rebinning array
// 	Float_t Rcalcbinning1[18] =  {5, 0.5, 2, 0.5, 4,0.5,1,0.5,3.5,0.5,4,0.5,5,0.5,4,1,11,15};
// 	Float_t Rcalcbinning3[19] = {1,2,16,17,27,28,38,39,41,43,65,68,71,74,79,80,81,82,87};
// 	Double_t Rbinning[87];
// 	Rbinning[0] = 0;
// 	for ( int j = 0; j < 18 ; j ++ ){
// 		for ( int i = Rcalcbinning3[j]; i<Rcalcbinning3[j+1]; i++){
// 			Rbinning[i] = Rbinning[i-1] +	 Rcalcbinning1[j];
// 		}
// 	}
//
	//making ratios = efficiency
	TH1D* histoRecConvGamma_Pt_rebin_Offline 	= (TH1D*)histoRecConvGamma_Pt_Offline->Clone(); //->Rebin(22,"recPteff",ptbinning);
	TH1D* histoConvGamma_Pt_rebin_Offline		= (TH1D*)histoConvGamma_Pt_Offline->Clone(); //->Rebin(22,"totPteff",ptbinning);
	TH1D* histoRatioRecConvGamma_Pt_Offline 	= (TH1D*)histoRecConvGamma_Pt_rebin_Offline->Clone("histoRatioRecConvGamma_Pt_Offline");
	histoRatioRecConvGamma_Pt_Offline->Sumw2();
	histoRatioRecConvGamma_Pt_Offline->Divide(histoRecConvGamma_Pt_rebin_Offline,histoConvGamma_Pt_rebin_Offline,1.,1.,"B");

	TH1D* histoRecConvGamma_PtlargerEta_rebin_Offline 	= (TH1D*)histoRecConvGamma_PtlargerEta_Offline->Clone(); //->Rebin(22,"recPtlargeEta",ptbinning);
	TH1D* histoConvGamma_PtlargerEta_rebin_Offline		= (TH1D*)histoConvGamma_PtlargerEta_Offline->Clone(); //->Rebin(22,"totPtlargeEta",ptbinning);
	TH1D* histoRatioRecConvGamma_PtlargerEta_Offline 	= (TH1D*)histoRecConvGamma_PtlargerEta_rebin_Offline->Clone("histoRatioRecConvGamma_PtlargerEta_Offline");
	histoRatioRecConvGamma_PtlargerEta_Offline->Sumw2();
	histoRatioRecConvGamma_PtlargerEta_Offline->Divide(histoRecConvGamma_PtlargerEta_rebin_Offline,histoConvGamma_PtlargerEta_rebin_Offline,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_Offline = (TH1D*)histoRecConvGamma_R_Offline->Clone("histoRatioRecConvGamma_R_Offline");
	histoRatioRecConvGamma_R_Offline->Sumw2();
	histoRatioRecConvGamma_R_Offline->Divide(histoRatioRecConvGamma_R_Offline,histoConvGamma_R_Offline,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_lowpT_Offline = (TH1D*)histoRecConvGamma_R_lowpT_Offline->Clone("histoRatioRecConvGamma_R_lowpT_Offline");
	histoRatioRecConvGamma_R_lowpT_Offline->Sumw2();
	histoRatioRecConvGamma_R_lowpT_Offline->Divide(histoRatioRecConvGamma_R_lowpT_Offline,histoConvGamma_R_lowpT_Offline,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_midpT_Offline = (TH1D*)histoRecConvGamma_R_midpT_Offline->Clone("histoRatioRecConvGamma_R_midpT_Offline");
	histoRatioRecConvGamma_R_midpT_Offline->Sumw2();
	histoRatioRecConvGamma_R_midpT_Offline->Divide(histoRatioRecConvGamma_R_midpT_Offline,histoConvGamma_R_midpT_Offline,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_highpT_Offline = (TH1D*)histoRecConvGamma_R_highpT_Offline->Clone("histoRatioRecConvGamma_R_highpT_Offline");
	histoRatioRecConvGamma_R_highpT_Offline->Sumw2();
	histoRatioRecConvGamma_R_highpT_Offline->Divide(histoRatioRecConvGamma_R_highpT_Offline,histoConvGamma_R_highpT_Offline,1.,1.,"B");



	TH1D* histoRatioRecConvGamma_Eta_Offline 	= (TH1D*)histoRecConvGamma_Eta_Offline->Clone("histoRatioRecConvGamma_Eta_Offline");
	histoRatioRecConvGamma_Eta_Offline->Sumw2();
	histoRatioRecConvGamma_Eta_Offline->Divide(histoRatioRecConvGamma_Eta_Offline,histoConvGamma_Eta_Offline,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_Phi_Offline	= (TH1D*)histoRecConvGamma_Phi_Offline->Clone("histoRatioRecConvGamma_Phi_Offline");
	histoRatioRecConvGamma_Phi_Offline->Sumw2();
	histoRatioRecConvGamma_Phi_Offline->Divide(histoRatioRecConvGamma_Phi_Offline,histoConvGamma_Phi_Offline,1.,1.,"B");

	//calculating the tot number of photons = 1st photon counted + 2nd time the same photon has been counted
	TH1D* histoTotalRecConvGamma_Pt_Offline 	= (TH1D*)histoRecConvGamma_Pt_Offline->Clone("histoRatioRecConvGamma_Pt_Offline");
	histoRecConvGamma_Pt_Offline->Add(histoRecConvGammaMult_Pt_Offline);
	TH1D* histoTotalRecConvGamma_R_Offline 		= (TH1D*)histoRecConvGamma_R_Offline->Clone("histoRatioRecConvGamma_R_Offline");
	histoTotalRecConvGamma_R_Offline->Add(histoRecConvGammaMult_R_Offline);
	TH1D* histoTotalRecConvGamma_Phi_Offline 	= (TH1D*)histoRecConvGamma_Phi_Offline->Clone("histoRatioRecConvGamma_Phi_Offline");
	histoTotalRecConvGamma_Phi_Offline->Add(histoRecConvGammaMult_Phi_Offline);

	//making ratio with the Mult histos = double counting ratio
	TH1D* histoRecConvGammaMult_Pt_rebin_Offline 	  = (TH1D*)histoRecConvGammaMult_Pt_Offline->Clone(); //->Rebin(22,"recPt",ptbinning);
	TH1D* histoTotalRecConvGamma_Pt_rebin_Offline	  = (TH1D*)histoTotalRecConvGamma_Pt_Offline->Clone(); //->Rebin(22,"totPt",ptbinning);
	TH1D* histoRatioRecConvGammaMult_Pt_Offline 	  = (TH1D*)histoRecConvGammaMult_Pt_rebin_Offline->Clone("histoRatioRecConvGammaMult_Pt_Offline");
	histoRatioRecConvGammaMult_Pt_Offline->Sumw2();
	histoRatioRecConvGammaMult_Pt_Offline->Divide(histoRecConvGammaMult_Pt_rebin_Offline,histoTotalRecConvGamma_Pt_rebin_Offline,1.,1.,"B");

	TH1D* histoRecConvGammaMult_R_rebin_Offline 	= (TH1D*)histoRecConvGammaMult_R_Offline->Clone(); //->Rebin(86,"recR",Rbinning);
	TH1D* histoTotalRecConvGamma_R_rebin_Offline	= (TH1D*)histoTotalRecConvGamma_R_Offline->Clone(); //->Rebin(86,"recR",Rbinning);
	TH1D* histoRatioRecConvGammaMult_R_Offline 	= (TH1D*)histoRecConvGammaMult_R_rebin_Offline->Clone("histoRatioRecConvGammaMult_R_Offline");
	histoRatioRecConvGammaMult_R_Offline->Sumw2();
	histoRatioRecConvGammaMult_R_Offline->Divide(histoRecConvGammaMult_R_rebin_Offline,histoTotalRecConvGamma_R_rebin_Offline,1.,1.,"B");

	TH1D* histoRatioRecConvGammaMult_Phi_Offline 	= (TH1D*)histoRecConvGammaMult_Phi_Offline->Clone("histoRatioRecConvGammaMult_Phi_Offline");
	histoRatioRecConvGammaMult_Phi_Offline->Sumw2();
	histoRatioRecConvGammaMult_Phi_Offline->Divide(histoRatioRecConvGammaMult_Phi_Offline,histoTotalRecConvGamma_Phi_Offline,1.,1.,"B");


	TH1D *ratioV0Finders = (TH1D*)histoRatioRecConvGamma_Pt_Offline->Clone();
	ratioV0Finders->Divide(histoRatioRecConvGamma_Pt_Offline,histoRatioRecConvGamma_Pt_OnFly,1.,1.,"");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    cout << "Plotting..." << endl;
    Double_t EtaRange[2] = {-1.5,1.5};
    Double_t PhiRange[2] = {0.,6.28};
    Double_t PtRange[2]  = {0.,15.};
    Double_t RRange[2]   = {0.,180.};

    Double_t minYRatio = 0.9;
    Double_t maxYRatio = 1.5;

    Color_t colorOnFly = kRed+1;
    Color_t colorOffline = kBlue+1;

	TString OutputNames[11] = {"histoRatioRecConvGamma_Pt_OnFly", "histoRatioRecConvGamma_PtlargerEta_OnFly", "histoRatioRecConvGamma_R_OnFly","histoRatioRecConvGamma_R_lowpT_OnFly", "histoRatioRecConvGamma_R_midpT_OnFly", "histoRatioRecConvGamma_R_highpT_OnFly", "histoRatioRecConvGamma_Eta_OnFly", "histoRatioRecConvGamma_Phi_OnFly", "histoRatioRecConvGammaMult_Pt_OnFly", "histoRatioRecConvGammaMult_R_OnFly", "histoRatioRecConvGammaMult_Phi_OnFly"};

    TCanvas* canvasV0Finder = new TCanvas("canvasV0Finder","",1300,1000);
    DrawGammaCanvasSettings( canvasV0Finder,  0.13, 0.02, 0.02, 0.09);
    TH2F * histo2DDummy = new TH2F("","",1000,0.,200.,1000,0.,2.);

    // draw efficiency vs Pt, smaller eta range
    canvasV0Finder->SetLogy(1);
    SetStyleHistoTH2ForGraphs(histo2DDummy,"#it{p}_{T} (GeV/#it{c})","#epsilon = #frac{rec. #gamma}{true #gamma} (|#eta| < 0.9)",0.04,0.04, 0.04,0.04, 1.,1.);
    histo2DDummy->GetXaxis()->SetRangeUser(PtRange[0],PtRange[1]);
    histo2DDummy->GetYaxis()->SetRangeUser(1e-2,5.);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoRatioRecConvGamma_Pt_OnFly, 20, 1.5,colorOnFly,colorOnFly);
        histoRatioRecConvGamma_Pt_OnFly->Draw("same,l");
        DrawGammaSetMarker(histoRatioRecConvGamma_Pt_Offline, 20, 1.5,colorOffline,colorOffline);
        histoRatioRecConvGamma_Pt_Offline->Draw("same,l");

        TLegend* legendV0Finder = GetAndSetLegend2(0.2,0.14,0.5,0.14+0.04*3.5,36);
        legendV0Finder->SetHeader("pp, #sqrt{s} = 5 TeV (2017)");
        legendV0Finder->AddEntry(histoRatioRecConvGamma_Pt_OnFly,"On-the-Fly","pl");
        legendV0Finder->AddEntry(histoRatioRecConvGamma_Pt_Offline,"Offline","pl");
        legendV0Finder->Draw();


    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/RecEffConvGamma_Pt.pdf",outputDir.Data()));

    // draw efficiency vs Pt, full eta range
    canvasV0Finder->SetLogy(1);
    SetStyleHistoTH2ForGraphs(histo2DDummy,"#it{p}_{T} (GeV/#it{c})","#epsilon = #frac{rec. #gamma}{true #gamma} (full #eta range)",0.04,0.04, 0.04,0.04, 1.,1.);
    histo2DDummy->GetXaxis()->SetRangeUser(PtRange[0],PtRange[1]);
    histo2DDummy->GetYaxis()->SetRangeUser(1e-2,5.);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoRatioRecConvGamma_PtlargerEta_OnFly, 20, 1.5,colorOnFly,colorOnFly);
        histoRatioRecConvGamma_PtlargerEta_OnFly->Draw("same,l");
        DrawGammaSetMarker(histoRatioRecConvGamma_PtlargerEta_Offline, 20, 1.5,colorOffline,colorOffline);
        histoRatioRecConvGamma_PtlargerEta_Offline->Draw("same,l");

        legendV0Finder->Draw();

    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/RecEffConvGamma_Pt_largerEta.pdf",outputDir.Data()));

    // draw efficiency vs R
    canvasV0Finder->SetLogy(0);
    SetStyleHistoTH2ForGraphs(histo2DDummy,"R (cm)","#epsilon = #frac{rec. #gamma}{true #gamma}",0.04,0.04, 0.04,0.04, 1.,1.);
    histo2DDummy->GetXaxis()->SetRangeUser(RRange[0],RRange[1]+5.);
    histo2DDummy->GetYaxis()->SetRangeUser(0.,.55);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoRatioRecConvGamma_R_OnFly, 20, 1.5,colorOnFly,colorOnFly);
        histoRatioRecConvGamma_R_OnFly->Draw("same,l");
        DrawGammaSetMarker(histoRatioRecConvGamma_R_Offline, 20, 1.5,colorOffline,colorOffline);
        histoRatioRecConvGamma_R_Offline->Draw("same,l");

        TLegend* legendV0Finder_R = GetAndSetLegend2(0.17,0.93-0.04*3.5,0.5,0.93,36);
        legendV0Finder_R->SetHeader("pp, #sqrt{s} = 5 TeV (2017)");
        legendV0Finder_R->AddEntry(histoRatioRecConvGamma_R_OnFly,"On-the-Fly","pl");
        legendV0Finder_R->AddEntry(histoRatioRecConvGamma_R_Offline,"Offline","pl");
        legendV0Finder_R->Draw();

    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/RecEffConvGamma_R.pdf",outputDir.Data()));

    // draw efficiency vs R, low pT
    histo2DDummy->GetYaxis()->SetRangeUser(0.,.55);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoRatioRecConvGamma_R_lowpT_OnFly, 20, 1.5,colorOnFly,colorOnFly);
        histoRatioRecConvGamma_R_lowpT_OnFly->Draw("same,l");
        DrawGammaSetMarker(histoRatioRecConvGamma_R_lowpT_Offline, 20, 1.5,colorOffline,colorOffline);
        histoRatioRecConvGamma_R_lowpT_Offline->Draw("same,l");

        TLegend* legendV0Finder_lowPt = GetAndSetLegend2(0.17,0.93-0.04*3.5,0.5,0.93,36);
        legendV0Finder_lowPt->SetHeader("#it{p}_{T} < 1.0 GeV/#it{c}");
        legendV0Finder_lowPt->AddEntry(histoRatioRecConvGamma_R_lowpT_OnFly,"On-the-Fly","pl");
        legendV0Finder_lowPt->AddEntry(histoRatioRecConvGamma_R_lowpT_Offline,"Offline","pl");
        legendV0Finder_lowPt->Draw();

    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/RecEffConvGamma_R_lowPt.pdf",outputDir.Data()));

    // draw efficiency vs R, mid pT
    histo2DDummy->GetYaxis()->SetRangeUser(0.,1.3);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoRatioRecConvGamma_R_midpT_OnFly, 20, 1.5,colorOnFly,colorOnFly);
        histoRatioRecConvGamma_R_midpT_OnFly->Draw("same,l");
        DrawGammaSetMarker(histoRatioRecConvGamma_R_midpT_Offline, 20, 1.5,colorOffline,colorOffline);
        histoRatioRecConvGamma_R_midpT_Offline->Draw("same,l");

        TLegend* legendV0Finder_midPt = GetAndSetLegend2(0.17,0.93-0.04*3.5,0.5,0.93,36);
        legendV0Finder_midPt->SetHeader("1.0 < #it{p}_{T} < 4.0 GeV/#it{c}");
        legendV0Finder_midPt->AddEntry(histoRatioRecConvGamma_R_midpT_OnFly,"On-the-Fly","pl");
        legendV0Finder_midPt->AddEntry(histoRatioRecConvGamma_R_midpT_Offline,"Offline","pl");
        legendV0Finder_midPt->Draw();

    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/RecEffConvGamma_R_midPt.pdf",outputDir.Data()));

    // draw efficiency vs R, high pT
//     canvasV0Finder->SetLogy(1);
    histo2DDummy->GetYaxis()->SetRangeUser(0.,1.3);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoRatioRecConvGamma_R_highpT_OnFly, 20, 1.5,colorOnFly,colorOnFly);
        histoRatioRecConvGamma_R_highpT_OnFly->Draw("same,l");
        DrawGammaSetMarker(histoRatioRecConvGamma_R_highpT_Offline, 20, 1.5,colorOffline,colorOffline);
        histoRatioRecConvGamma_R_highpT_Offline->Draw("same,l");

        TLegend* legendV0Finder_highPt = GetAndSetLegend2(0.17,0.93-0.04*3.5,0.5,0.93,36);
        legendV0Finder_highPt->SetHeader("#it{p}_{T} > 4.0 GeV/#it{c}");
        legendV0Finder_highPt->AddEntry(histoRatioRecConvGamma_R_highpT_OnFly,"On-the-Fly","pl");
        legendV0Finder_highPt->AddEntry(histoRatioRecConvGamma_R_highpT_Offline,"Offline","pl");
        legendV0Finder_highPt->Draw();

    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/RecEffConvGamma_R_highPt.pdf",outputDir.Data()));

    // draw efficiency vs phi
//     canvasV0Finder->SetLogy(1);
    SetStyleHistoTH2ForGraphs(histo2DDummy,"#phi (rad)","#epsilon = #frac{rec. #gamma}{true #gamma}",0.04,0.04, 0.04,0.04, 1.,1.3);
    histo2DDummy->GetXaxis()->SetRangeUser(PhiRange[0],PhiRange[1]);
    histo2DDummy->GetYaxis()->SetRangeUser(0.,.4);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoRatioRecConvGamma_Phi_OnFly, 20, 1.5,colorOnFly,colorOnFly);
        histoRatioRecConvGamma_Phi_OnFly->Draw("same,l");
        DrawGammaSetMarker(histoRatioRecConvGamma_Phi_Offline, 20, 1.5,colorOffline,colorOffline);
        histoRatioRecConvGamma_Phi_Offline->Draw("same,l");

        TLegend* legendV0Finder_Phi = GetAndSetLegend2(0.17,0.93-0.04*3.5,0.5,0.93,36);
        legendV0Finder_Phi->SetHeader("pp, #sqrt{s} = 5 TeV (2017)");
        legendV0Finder_Phi->AddEntry(histoRatioRecConvGamma_Phi_OnFly,"On-the-Fly","pl");
        legendV0Finder_Phi->AddEntry(histoRatioRecConvGamma_Phi_Offline,"Offline","pl");
        legendV0Finder_Phi->Draw();

    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/RecEffConvGamma_Phi.pdf",outputDir.Data()));

    // draw efficiency vs eta
//     canvasV0Finder->SetLogy(1);
    histo2DDummy = new TH2F("","",1000,-2.,2.,1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histo2DDummy,"#eta","#epsilon = #frac{rec. #gamma}{true #gamma}",0.04,0.04, 0.04,0.04, 1.,1.3);
    histo2DDummy->GetXaxis()->SetRangeUser(EtaRange[0],EtaRange[1]);
    histo2DDummy->GetYaxis()->SetRangeUser(0.,.4);
    histo2DDummy->Draw("copy");

        DrawGammaSetMarker(histoRatioRecConvGamma_Eta_OnFly, 20, 1.5,colorOnFly,colorOnFly);
        histoRatioRecConvGamma_Eta_OnFly->Draw("same,l");
        DrawGammaSetMarker(histoRatioRecConvGamma_Eta_Offline, 20, 1.5,colorOffline,colorOffline);
        histoRatioRecConvGamma_Eta_Offline->Draw("same,l");

        TLegend* legendV0Finder_Eta = GetAndSetLegend2(0.17,0.93-0.04*3.5,0.5,0.93,36);
        legendV0Finder_Eta->SetHeader("pp, #sqrt{s} = 5 TeV (2017)");
        legendV0Finder_Eta->AddEntry(histoRatioRecConvGamma_Eta_OnFly,"On-the-Fly","pl");
        legendV0Finder_Eta->AddEntry(histoRatioRecConvGamma_Eta_Offline,"Offline","pl");
        legendV0Finder_Eta->Draw();

    histo2DDummy->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/RecEffConvGamma_Eta.pdf",outputDir.Data()));


    canvasV0Finder->cd();

    TH2F * histo2DDummyRatio = new TH2F("","",1000,0.,20.,1000,0.,2.);
    SetStyleHistoTH2ForGraphs(histo2DDummyRatio,"#it{p}_{T} (GeV/#it{c})","efficiency ratio #frac{Offline}{On-the-Fly}",0.04,0.04, 0.04,0.04, 1.,1.2);
    histo2DDummyRatio->GetXaxis()->SetRangeUser(PtRange[0],PtRange[1]);
    histo2DDummyRatio->GetYaxis()->SetRangeUser(1e-2,5.);
    histo2DDummyRatio->Draw("copy");

        DrawGammaSetMarker(ratioV0Finders, 20, 1.5,colorOnFly,colorOnFly);
        ratioV0Finders->Draw("same,l");

        TLegend* legendV0FinderRatio = GetAndSetLegend2(0.17,0.93-0.04*3.5,0.5,0.93,36);
        legendV0FinderRatio->SetHeader("pp, #sqrt{s} = 5 TeV (2017)");
        legendV0FinderRatio->AddEntry(ratioV0Finders,"#frac{Offline}{On-the-Fly}","pl");
        legendV0FinderRatio->Draw();

    histo2DDummyRatio->Draw("same,axis");
    canvasV0Finder->Update();
    canvasV0Finder->SaveAs(Form("%s/V0FinderRatio.pdf",outputDir.Data()));
    delete canvasV0Finder;

// 	histoRatioRecConvGammaMult_Pt_OnFly->GetYaxis()->SetDecimals(1);
// 	histoRatioRecConvGammaMult_Pt_Offline->GetYaxis()->SetDecimals(1);
// 	DrawHisto(histoRatioRecConvGammaMult_Pt_OnFly,histoRatioRecConvGammaMult_Pt_Offline,outputDir,OutputNames[8],"Double counting ratio","#it{p}_{T} (GeV/#it{c})","#frac{2nd time same #gamma counted}{1st #gamma counted + 2nd time same #gamma counted}",0.,20., 0,0.03);
// 	DrawHisto(histoRatioRecConvGammaMult_R_OnFly,histoRatioRecConvGammaMult_R_Offline,outputDir,OutputNames[9],"Double counting ratio","R (cm)","#frac{2nd time same #gamma counted}{1st #gamma counted + 2nd time same #gamma counted}",0.,180, 0,0.2);
// 	DrawHisto(histoRatioRecConvGammaMult_Phi_OnFly,histoRatioRecConvGammaMult_Phi_Offline,outputDir,OutputNames[10],"Double counting ratio","#phi","#frac{2nd time same #gamma counted}{1st #gamma counted + 2nd time same #gamma counted}",0.,6.20, 0,0.05);

    cout << "Saving the histos" << endl;
	TFile* output = new TFile(Form("%s/HistosV0FindingQA_%s.root",outputDir.Data(), nameFolder.Data()),"RECREATE");
	output->cd();

        histoConvGamma_Pt_OnFly->Write();
        histoRecConvGamma_Pt_OnFly->Write();
        histoRatioRecConvGamma_Pt_OnFly->Write();
        histoRatioRecConvGamma_PtlargerEta_OnFly->Write();
        histoRatioRecConvGamma_R_OnFly->Write();
        histoConvGamma_R_highpT_OnFly->Write();
        histoConvGamma_R_lowpT_OnFly->Write();
        histoConvGamma_R_midpT_OnFly->Write();
        histoRecConvGamma_R_highpT_OnFly->Write();
        histoRecConvGamma_R_lowpT_OnFly->Write();
        histoRecConvGamma_R_midpT_OnFly->Write();
        histoRatioRecConvGamma_R_lowpT_OnFly->Write();
        histoRatioRecConvGamma_R_midpT_OnFly->Write();
        histoRatioRecConvGamma_R_highpT_OnFly->Write();
        histoRatioRecConvGamma_Eta_OnFly->Write();
        histoRatioRecConvGamma_Phi_OnFly->Write();
        histoRatioRecConvGammaMult_Pt_OnFly->Write();
        histoRatioRecConvGammaMult_R_OnFly->Write();
        histoRatioRecConvGammaMult_Phi_OnFly->Write();

        histoConvGamma_Pt_Offline->Write();
        histoRecConvGamma_Pt_Offline->Write();
        histoRatioRecConvGamma_Pt_Offline->Write();
        histoRatioRecConvGamma_PtlargerEta_Offline->Write();
        histoRatioRecConvGamma_R_Offline->Write();
        histoConvGamma_R_highpT_Offline->Write();
        histoConvGamma_R_lowpT_Offline->Write();
        histoConvGamma_R_midpT_Offline->Write();
        histoRecConvGamma_R_highpT_Offline->Write();
        histoRecConvGamma_R_lowpT_Offline->Write();
        histoRecConvGamma_R_midpT_Offline->Write();
        histoRatioRecConvGamma_R_lowpT_Offline->Write();
        histoRatioRecConvGamma_R_midpT_Offline->Write();
        histoRatioRecConvGamma_R_highpT_Offline->Write();
        histoRatioRecConvGamma_Eta_Offline->Write();
        histoRatioRecConvGamma_Phi_Offline->Write();
        histoRatioRecConvGammaMult_Pt_Offline->Write();
        histoRatioRecConvGammaMult_R_Offline->Write();
        histoRatioRecConvGammaMult_Phi_Offline->Write();

	output->Write();
	output->Close();


}


