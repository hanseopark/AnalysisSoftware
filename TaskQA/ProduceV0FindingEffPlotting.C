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
	histodummy->GetYaxis()->SetLabelOffset(0.005);
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
	if(Name.Contains("histoRatioRecConvGamma_Pt_0010")) c->SetLogy();	
	c->cd();
	
	TPad *p = new TPad("p","",0.0,0.0,1.,1.,0);   // gives the size of the histo areas 
	p->SetFillColor(0);
	p->GetFrame()->SetFillColor(0);
	p->SetLeftMargin(0.12);
	p->SetRightMargin(0.05);
	p->SetTopMargin(0.1);
	p->SetBottomMargin(0.12);
	//p->SetBorderMode(0);
	if(Name.Contains("histoRatioRecConvGamma_Pt_0010")) p->SetLogy();
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
		leg->AddEntry(histo1,("0-10%"));
		leg->AddEntry(histo2,("20-50%"));
		leg->Draw();

		TLatex *labelEnergy = new TLatex(0.15,0.18,"Pb - Pb, #sqrt{s_{NN}} = 2.76 TeV");
		SetStyleTLatex( labelEnergy, 0.04,4);
		labelEnergy->Draw();

	} else {
	
		TLegend* leg = new TLegend(0.8,0.76,0.9,0.86);
		leg->SetTextSize(0.04);
		leg->SetLineColor(0);
		leg->SetFillColor(0);
		leg->SetFillStyle(0);
		leg->AddEntry(histo1,("0-10%"));
		leg->AddEntry(histo2,("20-50%"));
		leg->Draw();

		TLatex *labelEnergy = new TLatex(0.15,0.82,"Pb - Pb, #sqrt{s_{NN}} = 2.76 TeV");
		SetStyleTLatex( labelEnergy, 0.04,4);
		labelEnergy->Draw();
	}

	c->SaveAs(Form("%s/%s.pdf",Output.Data(),Name.Data()));
	
	delete p;
	delete c;
	
	
		
}



void ProduceV0FindingEffPlotting(TString filename_0010 = "/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20150803/GA_PbPb_MC-222_20150804-0801/merge/GammaConvV1_162.root", TString filename_2050 = "/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20150803/GA_PbPb_MC-223_20150803-1855/merge/GammaConvV1_162.root", TString nameFolder ="GammaConvV1_162", TString EventCutString_0010 = "5010001", TString EventCutString_2050 = "5250001", TString outputDir = "V0EfficiencyStudies"){
//root -l -b -x -q 'ProduceV0FindingEffPlotting.C++("/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20150413/GA_PbPb_MC-202_20150415-0951/merge/GammaConvV1_GA_PbPb_MC_LHC14a1a_162.root","GammaConvV1","5010001","V0EfficiencyStudies")'
//root -l -b -x -q 'ProduceV0FindingEffPlotting.C++("/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20150413/GA_PbPb_MC-205_20150415-1147/merge/GammaConvV1_GA_PbPb_MC_LHC14a1b_162.root","GammaConvV1","5250001","V0EfficiencyStudies2050")'	
	
//root -l -b -x -q 'ProduceV0FindingEffPlotting.C++("/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20150413/GA_PbPb_MC-202_20150415-0951/merge/GammaConvV1_GA_PbPb_MC_LHC14a1a_162.root","/home/admin1/leardini/GridOutput/PbPb/Legotrain-vAN-20150413/GA_PbPb_MC-205_20150415-1147/merge/GammaConvV1_GA_PbPb_MC_LHC14a1b_162.root","GammaConvV1","5010001","5250001","V0EfficiencyStudiesTogether2")'
	
	gSystem->Exec("mkdir -p V0EfficiencyStudies");
	
	TFile *file_0010 = new TFile(filename_0010.Data()); 
	TList *list_0010 = (TList*) file_0010->Get(nameFolder.Data());
	TList *HistosV0Finding_0010 = (TList*)list_0010->FindObject(Form("V0FindingEfficiencyInput_%s", EventCutString_0010.Data()));

	//loading histos
	TH2F* histoConvGamma_Pt_R_0010 	= (TH2F*)HistosV0Finding_0010->FindObject("MCconvGamma_Pt_R");   
	TH2F* histoConvGamma_Pt_Eta_0010 = (TH2F*)HistosV0Finding_0010->FindObject("MCconvGamma_Pt_Eta");  
	TH2F* histoConvGamma_Pt_Phi_0010 = (TH2F*)HistosV0Finding_0010->FindObject("MCconvGamma_Pt_Phi"); 

	TH2F* histoRecConvGamma_Pt_R_0010 	= (TH2F*)HistosV0Finding_0010->FindObject("RecMCconvGamma_Pt_R");   
	TH2F* histoRecConvGamma_Pt_Eta_0010 	= (TH2F*)HistosV0Finding_0010->FindObject("RecMCconvGamma_Pt_Eta");  
	TH2F* histoRecConvGamma_Pt_Phi_0010 	= (TH2F*)HistosV0Finding_0010->FindObject("RecMCconvGamma_Pt_Phi");

	TH1D* histoRecConvGammaMult_Pt_0010 	= (TH1D*)HistosV0Finding_0010->FindObject("RecMCconvGammaMulti_Pt");   
	TH1D* histoRecConvGammaMult_R_0010 	= (TH1D*)HistosV0Finding_0010->FindObject("RecMCconvGammaMulti_R");  
	TH1D* histoRecConvGammaMult_Phi_0010 = (TH1D*)HistosV0Finding_0010->FindObject("RecMCconvGammaMulti_Phi");
	
	//projecting histos
	TH1D* histoConvGamma_Pt_0010 			= (TH1D*)histoConvGamma_Pt_R_0010->ProjectionX("histoConvGamma_Pt_0010"); //restricted to eta cut range
	TH1D* histoConvGamma_PtlargerEta_0010	= (TH1D*)histoConvGamma_Pt_Eta_0010->ProjectionX("histoConvGamma_PtlargerEta_0010"); //full eta range
		
	TH1D* histoConvGamma_R_0010				= (TH1D*)histoConvGamma_Pt_R_0010->ProjectionY("histoConvGamma_R_0010"); //full pT range
	TH1D* histoConvGamma_R_lowpT_0010 		= (TH1D*)histoConvGamma_Pt_R_0010->ProjectionY("histoConvGamma_R_lowpT_0010",histoConvGamma_Pt_R_0010->GetXaxis()->FindBin(0.),histoConvGamma_Pt_R_0010->GetXaxis()->FindBin(1.0)); // 0. < pT < 1. GeV/c 
	TH1D* histoConvGamma_R_midpT_0010 		= (TH1D*)histoConvGamma_Pt_R_0010->ProjectionY("histoConvGamma_R_midpT_0010",histoConvGamma_Pt_R_0010->GetXaxis()->FindBin(1.0),histoConvGamma_Pt_R_0010->GetXaxis()->FindBin(4.0)); // 1.0 < pT < 4. GeV/c 
	TH1D* histoConvGamma_R_highpT_0010		= (TH1D*)histoConvGamma_Pt_R_0010->ProjectionY("histoConvGamma_R_highpT_0010",histoConvGamma_Pt_R_0010->GetXaxis()->FindBin(4.0),histoConvGamma_Pt_R_0010->GetNbinsX()); // pT > 4. GeV/c
	
	TH1D* histoConvGamma_Eta_0010			= (TH1D*)histoConvGamma_Pt_Eta_0010->ProjectionY("histoConvGamma_Eta_0010");
	TH1D* histoConvGamma_Phi_0010 			= (TH1D*)histoConvGamma_Pt_Phi_0010->ProjectionY("histoConvGamma_Phi_0010");

	TH1D* histoRecConvGamma_Pt_0010 			= (TH1D*)histoRecConvGamma_Pt_R_0010->ProjectionX("histoRecConvGamma_Pt_0010"); //restricted to eta cut range
	TH1D* histoRecConvGamma_PtlargerEta_0010 	= (TH1D*)histoRecConvGamma_Pt_Eta_0010->ProjectionX("histoRecConvGamma_PtlargerEta_0010"); //full eta range
	TH1D* histoRecConvGamma_R_0010 				= (TH1D*)histoRecConvGamma_Pt_R_0010->ProjectionY("histoRecConvGamma_R_0010"); //full pT range
	TH1D* histoRecConvGamma_R_lowpT_0010 		= (TH1D*)histoRecConvGamma_Pt_R_0010->ProjectionY("histoRecConvGamma_R_lowpT_0010",histoRecConvGamma_Pt_R_0010->GetXaxis()->FindBin(0.),histoRecConvGamma_Pt_R_0010->GetXaxis()->FindBin(1.0)); // 0. < pT < 1. GeV/c 
	TH1D* histoRecConvGamma_R_midpT_0010 		= (TH1D*)histoRecConvGamma_Pt_R_0010->ProjectionY("histoRecConvGamma_R_midpT_0010",histoRecConvGamma_Pt_R_0010->GetXaxis()->FindBin(1.0),histoRecConvGamma_Pt_R_0010->GetXaxis()->FindBin(4.0)); // 1.0 < pT < 4. GeV/c 
	TH1D* histoRecConvGamma_R_highpT_0010		= (TH1D*)histoRecConvGamma_Pt_R_0010->ProjectionY("histoRecConvGamma_R_highpT_0010",histoRecConvGamma_Pt_R_0010->GetXaxis()->FindBin(4.0),histoRecConvGamma_Pt_R_0010->GetNbinsX()); // pT > 4. GeV/c
	
	TH1D* histoRecConvGamma_Eta_0010 		= (TH1D*)histoRecConvGamma_Pt_Eta_0010->ProjectionY("histoRecConvGamma_Eta_0010");
	TH1D* histoRecConvGamma_Phi_0010 		= (TH1D*)histoRecConvGamma_Pt_Phi_0010->ProjectionY("histoRecConvGamma_Phi_0010");

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
	TH1D* histoRecConvGamma_Pt_rebin_0010 	= (TH1D*)histoRecConvGamma_Pt_0010->Clone(); //Rebin(22,"recPteff",ptbinning);
	TH1D* histoConvGamma_Pt_rebin_0010		= (TH1D*)histoConvGamma_Pt_0010->Clone(); //->Rebin(22,"totPteff",ptbinning);
	TH1D* histoRatioRecConvGamma_Pt_0010 	= (TH1D*)histoRecConvGamma_Pt_rebin_0010->Clone("histoRatioRecConvGamma_Pt_0010");
	histoRatioRecConvGamma_Pt_0010->Sumw2();
	histoRatioRecConvGamma_Pt_0010->Divide(histoRecConvGamma_Pt_rebin_0010,histoConvGamma_Pt_rebin_0010,1.,1.,"B");
	
	TH1D* histoRecConvGamma_PtlargerEta_rebin_0010 		= (TH1D*)histoRecConvGamma_PtlargerEta_0010->Clone(); //->Rebin(22,"recPtlargeEta",ptbinning);
	TH1D* histoConvGamma_PtlargerEta_rebin_0010			= (TH1D*)histoConvGamma_PtlargerEta_0010->Clone(); //->Rebin(22,"totPtlargeEta",ptbinning);
	TH1D* histoRatioRecConvGamma_PtlargerEta_0010 		= (TH1D*)histoRecConvGamma_PtlargerEta_rebin_0010->Clone("histoRatioRecConvGamma_PtlargerEta_0010");
	histoRatioRecConvGamma_PtlargerEta_0010->Sumw2();
	histoRatioRecConvGamma_PtlargerEta_0010->Divide(histoRecConvGamma_PtlargerEta_rebin_0010,histoConvGamma_PtlargerEta_rebin_0010,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_0010 = (TH1D*)histoRecConvGamma_R_0010->Clone("histoRatioRecConvGamma_R_0010");
	histoRatioRecConvGamma_R_0010->Sumw2();
	histoRatioRecConvGamma_R_0010->Divide(histoRatioRecConvGamma_R_0010,histoConvGamma_R_0010,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_lowpT_0010 = (TH1D*)histoRecConvGamma_R_lowpT_0010->Clone("histoRatioRecConvGamma_R_lowpT_0010");
	histoRatioRecConvGamma_R_lowpT_0010->Sumw2();
	histoRatioRecConvGamma_R_lowpT_0010->Divide(histoRatioRecConvGamma_R_lowpT_0010,histoConvGamma_R_lowpT_0010,1.,1.,"B");	
	
	TH1D* histoRatioRecConvGamma_R_midpT_0010 = (TH1D*)histoRecConvGamma_R_midpT_0010->Clone("histoRatioRecConvGamma_R_midpT_0010");
	histoRatioRecConvGamma_R_midpT_0010->Sumw2();
	histoRatioRecConvGamma_R_midpT_0010->Divide(histoRatioRecConvGamma_R_midpT_0010,histoConvGamma_R_midpT_0010,1.,1.,"B");
	
	TH1D* histoRatioRecConvGamma_R_highpT_0010 = (TH1D*)histoRecConvGamma_R_highpT_0010->Clone("histoRatioRecConvGamma_R_highpT_0010");
	histoRatioRecConvGamma_R_highpT_0010->Sumw2();
	histoRatioRecConvGamma_R_highpT_0010->Divide(histoRatioRecConvGamma_R_highpT_0010,histoConvGamma_R_highpT_0010,1.,1.,"B");
	
	
	
	TH1D* histoRatioRecConvGamma_Eta_0010 		 = (TH1D*)histoRecConvGamma_Eta_0010->Clone("histoRatioRecConvGamma_Eta_0010");
	histoRatioRecConvGamma_Eta_0010->Sumw2();
	histoRatioRecConvGamma_Eta_0010->Divide(histoRatioRecConvGamma_Eta_0010,histoConvGamma_Eta_0010,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_Phi_0010		 = (TH1D*)histoRecConvGamma_Phi_0010->Clone("histoRatioRecConvGamma_Phi_0010");
	histoRatioRecConvGamma_Phi_0010->Sumw2();
	histoRatioRecConvGamma_Phi_0010->Divide(histoRatioRecConvGamma_Phi_0010,histoConvGamma_Phi_0010,1.,1.,"B");
	
	//calculating the tot number of photons = 1st photon counted + 2nd time the same photon has been counted
	TH1D* histoTotalRecConvGamma_Pt_0010 		= (TH1D*)histoRecConvGamma_Pt_0010->Clone("histoRatioRecConvGamma_Pt_0010");
	histoRecConvGamma_Pt_0010->Add(histoRecConvGammaMult_Pt_0010);
	TH1D* histoTotalRecConvGamma_R_0010 		= (TH1D*)histoRecConvGamma_R_0010->Clone("histoRatioRecConvGamma_R_0010");
	histoTotalRecConvGamma_R_0010->Add(histoRecConvGammaMult_R_0010);
	TH1D* histoTotalRecConvGamma_Phi_0010 		= (TH1D*)histoRecConvGamma_Phi_0010->Clone("histoRatioRecConvGamma_Phi_0010");
	histoTotalRecConvGamma_Phi_0010->Add(histoRecConvGammaMult_Phi_0010);
	

	
	//making ratio with the Mult histos = double counting ratio
	TH1D* histoRecConvGammaMult_Pt_rebin_0010 	  = (TH1D*)histoRecConvGammaMult_Pt_0010->Clone(); //->Rebin(22,"recPt",ptbinning);
	TH1D* histoTotalRecConvGamma_Pt_rebin_0010	  = (TH1D*)histoTotalRecConvGamma_Pt_0010->Clone(); //->Rebin(22,"totPt",ptbinning);
	TH1D* histoRatioRecConvGammaMult_Pt_0010 	  = (TH1D*)histoRecConvGammaMult_Pt_rebin_0010->Clone("histoRatioRecConvGammaMult_Pt_0010");
	histoRatioRecConvGammaMult_Pt_0010->Sumw2();
	histoRatioRecConvGammaMult_Pt_0010->Divide(histoRecConvGammaMult_Pt_rebin_0010,histoTotalRecConvGamma_Pt_rebin_0010,1.,1.,"B");

	TH1D* histoRecConvGammaMult_R_rebin_0010 	= (TH1D*)histoRecConvGammaMult_R_0010->Clone(); //->Rebin(86,"recR",Rbinning);
	TH1D* histoTotalRecConvGamma_R_rebin_0010	= (TH1D*)histoTotalRecConvGamma_R_0010->Clone(); //->Rebin(86,"recR",Rbinning);	
	TH1D* histoRatioRecConvGammaMult_R_0010 	= (TH1D*)histoRecConvGammaMult_R_rebin_0010->Clone("histoRatioRecConvGammaMult_R_0010");
	histoRatioRecConvGammaMult_R_0010->Sumw2();
	histoRatioRecConvGammaMult_R_0010->Divide(histoRecConvGammaMult_R_rebin_0010,histoTotalRecConvGamma_R_rebin_0010,1.,1.,"B");
	
	TH1D* histoRatioRecConvGammaMult_Phi_0010 	= (TH1D*)histoRecConvGammaMult_Phi_0010->Clone("histoRatioRecConvGammaMult_Phi_0010");
	histoRatioRecConvGammaMult_Phi_0010->Sumw2();
	histoRatioRecConvGammaMult_Phi_0010->Divide(histoRatioRecConvGammaMult_Phi_0010,histoTotalRecConvGamma_Phi_0010,1.,1.,"B");
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		
	TFile *file_2050 = new TFile(filename_2050.Data()); 
	TList *list_2050 = (TList*) file_2050->Get(nameFolder.Data());
	TList *HistosV0Finding_2050 = (TList*)list_2050->FindObject(Form("V0FindingEfficiencyInput_%s", EventCutString_2050.Data()));

	//loading histos
	TH2F* histoConvGamma_Pt_R_2050 	= (TH2F*)HistosV0Finding_2050->FindObject("MCconvGamma_Pt_R");   
	TH2F* histoConvGamma_Pt_Eta_2050 = (TH2F*)HistosV0Finding_2050->FindObject("MCconvGamma_Pt_Eta");  
	TH2F* histoConvGamma_Pt_Phi_2050 = (TH2F*)HistosV0Finding_2050->FindObject("MCconvGamma_Pt_Phi"); 

	TH2F* histoRecConvGamma_Pt_R_2050 	= (TH2F*)HistosV0Finding_2050->FindObject("RecMCconvGamma_Pt_R");   
	TH2F* histoRecConvGamma_Pt_Eta_2050 	= (TH2F*)HistosV0Finding_2050->FindObject("RecMCconvGamma_Pt_Eta");  
	TH2F* histoRecConvGamma_Pt_Phi_2050 	= (TH2F*)HistosV0Finding_2050->FindObject("RecMCconvGamma_Pt_Phi");

	TH1D* histoRecConvGammaMult_Pt_2050 	= (TH1D*)HistosV0Finding_2050->FindObject("RecMCconvGammaMulti_Pt");   
	TH1D* histoRecConvGammaMult_R_2050 	= (TH1D*)HistosV0Finding_2050->FindObject("RecMCconvGammaMulti_R");  
	TH1D* histoRecConvGammaMult_Phi_2050 = (TH1D*)HistosV0Finding_2050->FindObject("RecMCconvGammaMulti_Phi");
	
	//projecting histos
	TH1D* histoConvGamma_Pt_2050 			= (TH1D*)histoConvGamma_Pt_R_2050->ProjectionX("histoConvGamma_Pt_2050"); //restricted to eta cut range
	TH1D* histoConvGamma_PtlargerEta_2050	= (TH1D*)histoConvGamma_Pt_Eta_2050->ProjectionX("histoConvGamma_PtlargerEta_2050"); //full eta range
		
	TH1D* histoConvGamma_R_2050				= (TH1D*)histoConvGamma_Pt_R_2050->ProjectionY("histoConvGamma_R_2050"); //full pT range
	TH1D* histoConvGamma_R_lowpT_2050 		= (TH1D*)histoConvGamma_Pt_R_2050->ProjectionY("histoConvGamma_R_lowpT_2050",histoConvGamma_Pt_R_2050->GetXaxis()->FindBin(0.),histoConvGamma_Pt_R_2050->GetXaxis()->FindBin(1.0)); // 0. < pT < 1. GeV/c 
	TH1D* histoConvGamma_R_midpT_2050 		= (TH1D*)histoConvGamma_Pt_R_2050->ProjectionY("histoConvGamma_R_midpT_2050",histoConvGamma_Pt_R_2050->GetXaxis()->FindBin(1.0),histoConvGamma_Pt_R_2050->GetXaxis()->FindBin(4.0)); // 1.0 < pT < 4. GeV/c 
	TH1D* histoConvGamma_R_highpT_2050		= (TH1D*)histoConvGamma_Pt_R_2050->ProjectionY("histoConvGamma_R_highpT_2050",histoConvGamma_Pt_R_2050->GetXaxis()->FindBin(4.0),histoConvGamma_Pt_R_2050->GetNbinsX()); // pT > 4. GeV/c
	
	TH1D* histoConvGamma_Eta_2050			= (TH1D*)histoConvGamma_Pt_Eta_2050->ProjectionY("histoConvGamma_Eta_2050");
	TH1D* histoConvGamma_Phi_2050 			= (TH1D*)histoConvGamma_Pt_Phi_2050->ProjectionY("histoConvGamma_Phi_2050");

	TH1D* histoRecConvGamma_Pt_2050 			= (TH1D*)histoRecConvGamma_Pt_R_2050->ProjectionX("histoRecConvGamma_Pt_2050"); //restricted to eta cut range
	TH1D* histoRecConvGamma_PtlargerEta_2050 	= (TH1D*)histoRecConvGamma_Pt_Eta_2050->ProjectionX("histoRecConvGamma_PtlargerEta_2050"); //full eta range
	TH1D* histoRecConvGamma_R_2050 				= (TH1D*)histoRecConvGamma_Pt_R_2050->ProjectionY("histoRecConvGamma_R_2050"); //full pT range
	TH1D* histoRecConvGamma_R_lowpT_2050 		= (TH1D*)histoRecConvGamma_Pt_R_2050->ProjectionY("histoRecConvGamma_R_lowpT_2050",histoRecConvGamma_Pt_R_2050->GetXaxis()->FindBin(0.),histoRecConvGamma_Pt_R_2050->GetXaxis()->FindBin(1.0)); // 0. < pT < 1. GeV/c 
	TH1D* histoRecConvGamma_R_midpT_2050 		= (TH1D*)histoRecConvGamma_Pt_R_2050->ProjectionY("histoRecConvGamma_R_midpT_2050",histoRecConvGamma_Pt_R_2050->GetXaxis()->FindBin(1.0),histoRecConvGamma_Pt_R_2050->GetXaxis()->FindBin(4.0)); // 1.0 < pT < 4. GeV/c 
	TH1D* histoRecConvGamma_R_highpT_2050		= (TH1D*)histoRecConvGamma_Pt_R_2050->ProjectionY("histoRecConvGamma_R_highpT_2050",histoRecConvGamma_Pt_R_2050->GetXaxis()->FindBin(4.0),histoRecConvGamma_Pt_R_2050->GetNbinsX()); // pT > 4. GeV/c
	
	TH1D* histoRecConvGamma_Eta_2050 		= (TH1D*)histoRecConvGamma_Pt_Eta_2050->ProjectionY("histoRecConvGamma_Eta_2050");
	TH1D* histoRecConvGamma_Phi_2050 		= (TH1D*)histoRecConvGamma_Pt_Phi_2050->ProjectionY("histoRecConvGamma_Phi_2050");

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
	TH1D* histoRecConvGamma_Pt_rebin_2050 	= (TH1D*)histoRecConvGamma_Pt_2050->Clone(); //->Rebin(22,"recPteff",ptbinning);
	TH1D* histoConvGamma_Pt_rebin_2050		= (TH1D*)histoConvGamma_Pt_2050->Clone(); //->Rebin(22,"totPteff",ptbinning);
	TH1D* histoRatioRecConvGamma_Pt_2050 	= (TH1D*)histoRecConvGamma_Pt_rebin_2050->Clone("histoRatioRecConvGamma_Pt_2050");
	histoRatioRecConvGamma_Pt_2050->Sumw2();
	histoRatioRecConvGamma_Pt_2050->Divide(histoRecConvGamma_Pt_rebin_2050,histoConvGamma_Pt_rebin_2050,1.,1.,"B");
	
	TH1D* histoRecConvGamma_PtlargerEta_rebin_2050 		= (TH1D*)histoRecConvGamma_PtlargerEta_2050->Clone(); //->Rebin(22,"recPtlargeEta",ptbinning);
	TH1D* histoConvGamma_PtlargerEta_rebin_2050			= (TH1D*)histoConvGamma_PtlargerEta_2050->Clone(); //->Rebin(22,"totPtlargeEta",ptbinning);
	TH1D* histoRatioRecConvGamma_PtlargerEta_2050 		= (TH1D*)histoRecConvGamma_PtlargerEta_rebin_2050->Clone("histoRatioRecConvGamma_PtlargerEta_2050");
	histoRatioRecConvGamma_PtlargerEta_2050->Sumw2();
	histoRatioRecConvGamma_PtlargerEta_2050->Divide(histoRecConvGamma_PtlargerEta_rebin_2050,histoConvGamma_PtlargerEta_rebin_2050,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_2050 = (TH1D*)histoRecConvGamma_R_2050->Clone("histoRatioRecConvGamma_R_2050");
	histoRatioRecConvGamma_R_2050->Sumw2();
	histoRatioRecConvGamma_R_2050->Divide(histoRatioRecConvGamma_R_2050,histoConvGamma_R_2050,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_lowpT_2050 = (TH1D*)histoRecConvGamma_R_lowpT_2050->Clone("histoRatioRecConvGamma_R_lowpT_2050");
	histoRatioRecConvGamma_R_lowpT_2050->Sumw2();
	histoRatioRecConvGamma_R_lowpT_2050->Divide(histoRatioRecConvGamma_R_lowpT_2050,histoConvGamma_R_lowpT_2050,1.,1.,"B");	
	
	TH1D* histoRatioRecConvGamma_R_midpT_2050 = (TH1D*)histoRecConvGamma_R_midpT_2050->Clone("histoRatioRecConvGamma_R_midpT_2050");
	histoRatioRecConvGamma_R_midpT_2050->Sumw2();
	histoRatioRecConvGamma_R_midpT_2050->Divide(histoRatioRecConvGamma_R_midpT_2050,histoConvGamma_R_midpT_2050,1.,1.,"B");
	
	TH1D* histoRatioRecConvGamma_R_highpT_2050 = (TH1D*)histoRecConvGamma_R_highpT_2050->Clone("histoRatioRecConvGamma_R_highpT_2050");
	histoRatioRecConvGamma_R_highpT_2050->Sumw2();
	histoRatioRecConvGamma_R_highpT_2050->Divide(histoRatioRecConvGamma_R_highpT_2050,histoConvGamma_R_highpT_2050,1.,1.,"B");
	
	
	
	TH1D* histoRatioRecConvGamma_Eta_2050 		 = (TH1D*)histoRecConvGamma_Eta_2050->Clone("histoRatioRecConvGamma_Eta_2050");
	histoRatioRecConvGamma_Eta_2050->Sumw2();
	histoRatioRecConvGamma_Eta_2050->Divide(histoRatioRecConvGamma_Eta_2050,histoConvGamma_Eta_2050,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_Phi_2050		 = (TH1D*)histoRecConvGamma_Phi_2050->Clone("histoRatioRecConvGamma_Phi_2050");
	histoRatioRecConvGamma_Phi_2050->Sumw2();
	histoRatioRecConvGamma_Phi_2050->Divide(histoRatioRecConvGamma_Phi_2050,histoConvGamma_Phi_2050,1.,1.,"B");
	
	//calculating the tot number of photons = 1st photon counted + 2nd time the same photon has been counted
	TH1D* histoTotalRecConvGamma_Pt_2050 		= (TH1D*)histoRecConvGamma_Pt_2050->Clone("histoRatioRecConvGamma_Pt_2050");
	histoRecConvGamma_Pt_2050->Add(histoRecConvGammaMult_Pt_2050);
	TH1D* histoTotalRecConvGamma_R_2050 		= (TH1D*)histoRecConvGamma_R_2050->Clone("histoRatioRecConvGamma_R_2050");
	histoTotalRecConvGamma_R_2050->Add(histoRecConvGammaMult_R_2050);
	TH1D* histoTotalRecConvGamma_Phi_2050 		= (TH1D*)histoRecConvGamma_Phi_2050->Clone("histoRatioRecConvGamma_Phi_2050");
	histoTotalRecConvGamma_Phi_2050->Add(histoRecConvGammaMult_Phi_2050);
	

	
	//making ratio with the Mult histos = double counting ratio
	TH1D* histoRecConvGammaMult_Pt_rebin_2050 	  = (TH1D*)histoRecConvGammaMult_Pt_2050->Clone(); //->Rebin(22,"recPt",ptbinning);
	TH1D* histoTotalRecConvGamma_Pt_rebin_2050	  = (TH1D*)histoTotalRecConvGamma_Pt_2050->Clone(); //->Rebin(22,"totPt",ptbinning);
	TH1D* histoRatioRecConvGammaMult_Pt_2050 	  = (TH1D*)histoRecConvGammaMult_Pt_rebin_2050->Clone("histoRatioRecConvGammaMult_Pt_2050");
	histoRatioRecConvGammaMult_Pt_2050->Sumw2();
	histoRatioRecConvGammaMult_Pt_2050->Divide(histoRecConvGammaMult_Pt_rebin_2050,histoTotalRecConvGamma_Pt_rebin_2050,1.,1.,"B");

	TH1D* histoRecConvGammaMult_R_rebin_2050 	= (TH1D*)histoRecConvGammaMult_R_2050->Clone(); //->Rebin(86,"recR",Rbinning);
	TH1D* histoTotalRecConvGamma_R_rebin_2050	= (TH1D*)histoTotalRecConvGamma_R_2050->Clone(); //->Rebin(86,"recR",Rbinning);	
	TH1D* histoRatioRecConvGammaMult_R_2050 	= (TH1D*)histoRecConvGammaMult_R_rebin_2050->Clone("histoRatioRecConvGammaMult_R_2050");
	histoRatioRecConvGammaMult_R_2050->Sumw2();
	histoRatioRecConvGammaMult_R_2050->Divide(histoRecConvGammaMult_R_rebin_2050,histoTotalRecConvGamma_R_rebin_2050,1.,1.,"B");
	
	TH1D* histoRatioRecConvGammaMult_Phi_2050 	= (TH1D*)histoRecConvGammaMult_Phi_2050->Clone("histoRatioRecConvGammaMult_Phi_2050");
	histoRatioRecConvGammaMult_Phi_2050->Sumw2();
	histoRatioRecConvGammaMult_Phi_2050->Divide(histoRatioRecConvGammaMult_Phi_2050,histoTotalRecConvGamma_Phi_2050,1.,1.,"B");
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	TFile* output = new TFile(Form("%s/HistosV0FindingQA_%s.root",outputDir.Data(), nameFolder.Data()),"RECREATE");
	output->cd();
	
	histoConvGamma_Pt_0010->Write();
	histoRecConvGamma_Pt_0010->Write();
	histoRatioRecConvGamma_Pt_0010->Write();
	histoRatioRecConvGamma_PtlargerEta_0010->Write();
	histoRatioRecConvGamma_R_0010->Write();
	histoConvGamma_R_highpT_0010->Write();
	histoConvGamma_R_lowpT_0010->Write();
	histoConvGamma_R_midpT_0010->Write();
	histoRecConvGamma_R_highpT_0010->Write();
	histoRecConvGamma_R_lowpT_0010->Write();
	histoRecConvGamma_R_midpT_0010->Write();
	histoRatioRecConvGamma_R_lowpT_0010->Write();
	histoRatioRecConvGamma_R_midpT_0010->Write();
	histoRatioRecConvGamma_R_highpT_0010->Write();
	histoRatioRecConvGamma_Eta_0010->Write();
	histoRatioRecConvGamma_Phi_0010->Write();
	histoRatioRecConvGammaMult_Pt_0010->Write();
	histoRatioRecConvGammaMult_R_0010->Write();
	histoRatioRecConvGammaMult_Phi_0010->Write();
	
	histoConvGamma_Pt_2050->Write();
	histoRecConvGamma_Pt_2050->Write();
	histoRatioRecConvGamma_Pt_2050->Write();
	histoRatioRecConvGamma_PtlargerEta_2050->Write();
	histoRatioRecConvGamma_R_2050->Write();
	histoConvGamma_R_highpT_2050->Write();
	histoConvGamma_R_lowpT_2050->Write();
	histoConvGamma_R_midpT_2050->Write();
	histoRecConvGamma_R_highpT_2050->Write();
	histoRecConvGamma_R_lowpT_2050->Write();
	histoRecConvGamma_R_midpT_2050->Write();
	histoRatioRecConvGamma_R_lowpT_2050->Write();
	histoRatioRecConvGamma_R_midpT_2050->Write();
	histoRatioRecConvGamma_R_highpT_2050->Write();
	histoRatioRecConvGamma_Eta_2050->Write();
	histoRatioRecConvGamma_Phi_2050->Write();
	histoRatioRecConvGammaMult_Pt_2050->Write();
	histoRatioRecConvGammaMult_R_2050->Write();
	histoRatioRecConvGammaMult_Phi_2050->Write();

	
	output->Write();
	output->Close();
	
	//Plotting

	TString OutputNames[11]={"histoRatioRecConvGamma_Pt_0010","histoRatioRecConvGamma_PtlargerEta_0010","histoRatioRecConvGamma_R_0010","histoRatioRecConvGamma_R_lowpT_0010","histoRatioRecConvGamma_R_midpT_0010","histoRatioRecConvGamma_R_highpT_0010","histoRatioRecConvGamma_Eta_0010","histoRatioRecConvGamma_Phi_0010","histoRatioRecConvGammaMult_Pt_0010","histoRatioRecConvGammaMult_R_0010","histoRatioRecConvGammaMult_Phi_0010"};//,"histoRatioRecConvGamma_Pt_2050","histoRatioRecConvGamma_PtlargerEta_2050","histoRatioRecConvGamma_R_2050","histoRatioRecConvGamma_R_lowpT_2050","histoRatioRecConvGamma_R_midpT_2050","histoRatioRecConvGamma_R_highpT_2050","histoRatioRecConvGamma_Eta_2050","histoRatioRecConvGamma_Phi_2050","histoRatioRecConvGammaMult_Pt_2050","histoRatioRecConvGammaMult_R_2050","histoRatioRecConvGammaMult_Phi_2050"};
	//	TString OutputNames[11]={"histoRatioRecConvGamma_Pt_2050","histoRatioRecConvGamma_PtlargerEta_2050","histoRatioRecConvGamma_R_2050","histoRatioRecConvGamma_R_lowpT_2050","histoRatioRecConvGamma_R_midpT_2050","histoRatioRecConvGamma_R_highpT_2050","histoRatioRecConvGamma_Eta_2050","histoRatioRecConvGamma_Phi_2050","histoRatioRecConvGammaMult_Pt_2050","histoRatioRecConvGammaMult_R_2050","histoRatioRecConvGammaMult_Phi_2050"};

	/*
	 * DrawHisto input
	 * 1. TH1D Pointer to Histogram
	 * 2. TString OutputDir
	 * 3. TString Name for saving file
	 * 4. TString Title to print in Histogram
	 * 5. TString X-Axis Title
	 * 6. TString Y-Axis Title
	 * 7. Float_t x min
	 * 8. Float_t x max
	 * 9. Float_t y min
	 * 10.Float_t ymax
	 * 11.Int_t color
	 * 12.Bool_t 0 for bins with errors and 1 for line
	 */
	
	
	DrawHisto(histoRatioRecConvGamma_Pt_0010,histoRatioRecConvGamma_Pt_2050,outputDir,OutputNames[0],"","#it{p}_{T} (GeV/#it{c})","efficiency = #frac{reconstructed #gamma}{true #gamma} (|#eta| < 0.9)",0.,20, 1e-2,5.);	

	TH1D *RatioCent = (TH1D*)histoRatioRecConvGamma_Pt_2050->Clone();
	RatioCent->Divide(histoRatioRecConvGamma_Pt_2050,histoRatioRecConvGamma_Pt_0010,1.,1.,"");
	DrawHisto(RatioCent,RatioCent,outputDir,"ratioCent","","#it{p}_{T} (GeV/#it{c})","efficiency ratio #frac{20-50% }{0-10% } (|#eta| < 0.9)",0.,20, 0.5,1.5);	


	DrawHisto(histoRatioRecConvGamma_PtlargerEta_0010,histoRatioRecConvGamma_PtlargerEta_2050,outputDir,OutputNames[1],"","#it{p}_{T} (GeV/#it{c})","efficiency = #frac{reconstructed #gamma}{true #gamma} (full #eta range)",0.0001,20., 0,1);
	DrawHisto(histoRatioRecConvGamma_R_0010,histoRatioRecConvGamma_R_2050,outputDir,OutputNames[2],"","R (cm)","#frac{reconstructed #gamma}{true #gamma}",0.0001,180., 0,1);
	DrawHisto(histoRatioRecConvGamma_R_lowpT_0010,histoRatioRecConvGamma_R_lowpT_2050,outputDir,OutputNames[3],"#it{p}_{T} < 1.0 GeV/#it{c}","R (cm)","#frac{reconstructed #gamma}{true #gamma}",0.,180, 0,1);
	DrawHisto(histoRatioRecConvGamma_R_midpT_0010,histoRatioRecConvGamma_R_midpT_2050,outputDir,OutputNames[4],"1.0 < #it{p}_{T} < 4.0 GeV/#it{c}","R (cm)","#frac{reconstructed #gamma}{true #gamma}",0.,180, 0,1);	
	DrawHisto(histoRatioRecConvGamma_R_highpT_0010,histoRatioRecConvGamma_R_highpT_2050,outputDir,OutputNames[5],"#it{p}_{T} > 4.0 GeV/#it{c}","R (cm)","#frac{reconstructed #gamma}{true #gamma}",0.,180., 0,1);
	DrawHisto(histoRatioRecConvGamma_Eta_0010,histoRatioRecConvGamma_Eta_2050,outputDir,OutputNames[6],"","#eta","#frac{reconstructed #gamma}{true #gamma}",-1.5,1.5, 0,0.6);
	DrawHisto(histoRatioRecConvGamma_Phi_0010,histoRatioRecConvGamma_Phi_2050,outputDir,OutputNames[7],"","#phi","#frac{reconstructed #gamma}{true #gamma}",0.,6.20, 0.,0.7,1);
	histoRatioRecConvGammaMult_Pt_0010->GetYaxis()->SetDecimals(1);
	histoRatioRecConvGammaMult_Pt_2050->GetYaxis()->SetDecimals(1);
	DrawHisto(histoRatioRecConvGammaMult_Pt_0010,histoRatioRecConvGammaMult_Pt_2050,outputDir,OutputNames[8],"Double counting ratio","#it{p}_{T} (GeV/#it{c})","#frac{2nd time same #gamma counted}{1st #gamma counted + 2nd time same #gamma counted}",0.,20., 0,0.03);
	DrawHisto(histoRatioRecConvGammaMult_R_0010,histoRatioRecConvGammaMult_R_2050,outputDir,OutputNames[9],"Double counting ratio","R (cm)","#frac{2nd time same #gamma counted}{1st #gamma counted + 2nd time same #gamma counted}",0.,180, 0,0.2);
	DrawHisto(histoRatioRecConvGammaMult_Phi_0010,histoRatioRecConvGammaMult_Phi_2050,outputDir,OutputNames[10],"Double counting ratio","#phi","#frac{2nd time same #gamma counted}{1st #gamma counted + 2nd time same #gamma counted}",0.,6.20, 0,0.05);

// 	DrawHisto(histoRatioRecConvGamma_Pt_2050,outputDir,,"","#it{p}_{T} (GeV/#it{c})","efficiency = #frac{reconstructed #gamma}{true #gamma} (|#eta| < 0.9)",0,20, 0,1);
// 	DrawHisto(histoRatioRecConvGamma_PtlargerEta_2050,outputDir,OutputNames[12],"","#it{p}_{T} (GeV/#it{c})","efficiency = #frac{reconstructed #gamma}{true #gamma} (full #eta range)",0.,20., 0,1);
// 	DrawHisto(histoRatioRecConvGamma_R_2050,outputDir,OutputNames[13],"","R (cm)","#frac{reconstructed #gamma}{true #gamma}",0.,180., 0,1);
// 	DrawHisto(histoRatioRecConvGamma_R_lowpT_2050,outputDir,OutputNames[14],"#it{p}_{T} < 1.0 GeV/#it{c}","R (cm)","#frac{reconstructed #gamma}{true #gamma}",0.,180, 0,1);
// 	DrawHisto(histoRatioRecConvGamma_R_midpT_2050,outputDir,OutputNames[15],"1.0 < #it{p}_{T} < 4.0 GeV/#it{c}","R (cm)","#frac{reconstructed #gamma}{true #gamma}",0.,180, 0,1);
// 	DrawHisto(histoRatioRecConvGamma_R_highpT_2050,outputDir,OutputNames[16],"#it{p}_{T} > 4.0 GeV/#it{c}","R (cm)","#frac{reconstructed #gamma}{true #gamma}",0.,180., 0,1);
// 	DrawHisto(histoRatioRecConvGamma_Eta_2050,outputDir,OutputNames[17],"","#eta","#frac{reconstructed #gamma}{true #gamma}",-1.5,1.5, 0,0.5);
// 	DrawHisto(histoRatioRecConvGamma_Phi_2050,outputDir,OutputNames[18],"","#phi","#frac{reconstructed #gamma}{true #gamma}",0.,6.20, 0.,0.7,1);
// 	DrawHisto(histoRatioRecConvGammaMult_Pt_2050,outputDir,OutputNames[19],"Double counting ratio","#it{p}_{T} (GeV/#it{c})","#frac{2nd time same #gamma counted}{1st #gamma counted + 2nd time same #gamma counted}",0.,20., 0,0.03);
// 	DrawHisto(histoRatioRecConvGammaMult_R_2050,outputDir,OutputNames[20],"Double counting ratio","R (cm)","#frac{2nd time same #gamma counted}{1st #gamma counted + 2nd time same #gamma counted}",0.,180, 0,0.2);
// 	DrawHisto(histoRatioRecConvGammaMult_Phi_2050,outputDir,OutputNames[21],"Double counting ratio","#phi","#frac{2nd time same #gamma counted}{1st #gamma counted + 2nd time same #gamma counted}",0.,6.20, 0,0.05);
	
	
}
	
	
