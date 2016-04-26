/*********************************************************************************************************
 ******  provided by Gamma Conversion Group, PWG4:                                                 ********
 ******                 Kathrin Koch, kkoch@physi.uni-heidelberg.de                                                                ********
 ******         Friederike Bock, fbock@physi.uni-heidelberg.de                                                     ********
 **********************************************************************************************************
 ***     This Macro can be used for the material analysis of the ALICE ITS and TPC up to an R of 180 cm   ***
 ***     To use this macro the GammaConversionTaks needs to be run first.                                                  ***
 ***     The macro can be run with 2 input files, the second one should be a Montecarlo file, if it is a  ***
 ***     data-file, modifications are necessary, taking out the histograms, which are only produced if    *** 
 ***     the ConversionTask is run on MC. They usually have MCtruth somewhere in its name.                ***
 **********************************************************************************************************/

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
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
// #include "../CommonHeaders/PlottingGammaConversionHistos.h"
// #include "../CommonHeaders/PlottingGammaConversionAdditional.h"
// #include "../CommonHeaders/FittingGammaConversion.h"
// #include "../CommonHeaders/ConversionFunctions.h"

TString textGenerator;
TString collisionSystem;
TString textPeriod;
TString textDate;
Bool_t thesis = kFALSE;

using namespace std;

/***************************************************************************************************** 
DrawStructure() draws the structure of the Inner Alice Detectors labeled for a xy - Plot
******************************************************************************************************
******************************************************************************************************/

void DrawStructureNew(){
	TLatex *ssdText = new TLatex(0.16,0.87,"SSD");
	ssdText->SetTextColor(kWhite);
	ssdText->SetNDC();
	ssdText->SetTextFont(72);
	ssdText->SetTextSize(0.03);
	ssdText->SetLineWidth(4);
	ssdText->Draw();
	
	TLatex *sddText = new TLatex(0.14,0.78,"SDD");
	sddText->SetTextColor(kWhite);
	sddText->SetNDC();
	sddText->SetTextFont(72);
	sddText->SetTextSize(0.03);
	sddText->SetLineWidth(4);
	sddText->Draw();
	
	TLatex *spdText = new TLatex(0.14,0.27,"SPD");
	spdText->SetTextColor(kWhite);
	spdText->SetNDC();
	spdText->SetTextFont(72);
	spdText->SetTextSize(0.03);
	spdText->SetLineWidth(4);
	spdText->Draw();
	
	TLatex *tpcRText = new TLatex(0.52,0.095,"TPC Rods");
	tpcRText->SetTextColor(kWhite);
	tpcRText->SetNDC();
	tpcRText->SetTextFont(72);
	tpcRText->SetTextSize(0.03);
	tpcRText->SetLineWidth(4);
	tpcRText->Draw();
	
	TLatex *tpcIText = new TLatex(0.14,0.18,"TPC inner");
	tpcIText->SetTextColor(kWhite);
	tpcIText->SetNDC();
	tpcIText->SetTextFont(72);
	tpcIText->SetTextSize(0.03);
	tpcIText->SetLineWidth(4);
	tpcIText->Draw();
	
	TLatex *tpcIFText = new TLatex(0.14,0.15,"field cage");
	tpcIFText->SetTextColor(kWhite);
	tpcIFText->SetNDC();
	tpcIFText->SetTextFont(72);
	tpcIFText->SetTextSize(0.03);
	tpcIFText->SetLineWidth(4);
	tpcIFText->Draw();
	
	TLatex *tpcIFVText = new TLatex(0.16,0.12,"vessel");
	tpcIFVText->SetTextColor(kWhite);
	tpcIFVText->SetNDC();
	tpcIFVText->SetTextFont(72);
	tpcIFVText->SetTextSize(0.03);
	tpcIFVText->SetLineWidth(4);
	tpcIFVText->Draw();
	
	TLatex *tpc1IText = new TLatex(0.705,0.18,"TPC inner");
	tpc1IText->SetTextColor(kWhite);
	tpc1IText->SetNDC();
	tpc1IText->SetTextFont(72);
	tpc1IText->SetTextSize(0.03);
	tpc1IText->SetLineWidth(4);
	tpc1IText->Draw();
	
	TLatex *tpcICText = new TLatex(0.69,0.15,"containment");
	tpcICText->SetTextColor(kWhite);
	tpcICText->SetNDC();
	tpcICText->SetTextFont(72);
	tpcICText->SetTextSize(0.03);
	tpcICText->SetLineWidth(4);
	tpcICText->Draw();
	
	TLatex *tpcICVText = new TLatex(0.72,0.12,"vessel");
	tpcICVText->SetTextColor(kWhite);
	tpcICVText->SetNDC();
	tpcICVText->SetTextFont(72);
	tpcICVText->SetTextSize(0.03);
	tpcICVText->SetLineWidth(10);
	tpcICVText->Draw();
	
	TLatex *tpcGasText = new TLatex(0.79,0.30,"TPC");
	tpcGasText->SetTextColor(kWhite);
	tpcGasText->SetNDC();
	tpcGasText->SetTextFont(72);
	tpcGasText->SetTextSize(0.03);
	tpcGasText->SetLineWidth(10);
	tpcGasText->Draw();
	
	TLatex *tpcGas2Text = new TLatex(0.77,0.27,"drift gas");
	tpcGas2Text->SetTextColor(kWhite);
	tpcGas2Text->SetNDC();
	tpcGas2Text->SetTextFont(72);
	tpcGas2Text->SetTextSize(0.03);
	tpcGas2Text->SetLineWidth(10);
	tpcGas2Text->Draw();
	
	
	TArrow *arrow = new TArrow(-82.,90.,-11.99843,38.629599,0.02,">"); //SSD arrow
	arrow->SetLineColor(kWhite);
	arrow->SetFillColor(kWhite);
	arrow->SetFillStyle(1001);
	arrow->SetLineWidth(2.);
	arrow->Draw();
	
	TArrow *arrow1 = new TArrow(-90.,65.,-11.99843,25.,0.02,">"); //SDD arrow
	arrow1->SetLineColor(kWhite);
	arrow1->SetFillColor(kWhite);
	arrow1->SetFillStyle(1001);
	arrow1->SetLineWidth(2.);
	arrow1->Draw();
	
	
	TArrow *arrow2 = new TArrow(-90.,-65.,-7.,2.,0.02,">");  //SPD arrow
	arrow2->SetLineColor(kWhite);
	arrow2->SetFillColor(kWhite);
	arrow2->SetFillStyle(1001);
	arrow2->SetLineWidth(2.);
	arrow2->Draw();
	
	TArrow *arrow3 = new TArrow(-70.,-100.,-30.,-75.,0.02,">"); //TPC field cage vessel arrow
	arrow3->SetLineColor(kWhite);
	arrow3->SetFillColor(kWhite);
	arrow3->SetFillStyle(1001);
	arrow3->SetLineWidth(2.);
	arrow3->Draw();
	
	
	TArrow *arrow4 = new TArrow(20.,-110.,15.,-83.,0.02,">"); //TPC rods arrow
	arrow4->SetLineColor(kWhite);
	arrow4->SetFillColor(kWhite);
	arrow4->SetFillStyle(1001);
	arrow4->SetLineWidth(2.);
	arrow4->Draw();
	
	TArrow *arrow5 = new TArrow(80.,-85.,50.,-38.,0.02,">");// TPC inner constainment vessel arrow
	arrow5->SetLineColor(kWhite);
	arrow5->SetFillColor(kWhite);
	arrow5->SetFillStyle(1001);
	arrow5->SetLineWidth(2.);
	arrow5->Draw();
	
	
	TArrow *arrow6 = new TArrow(100.,-50.,90.,-10.,0.02,">");// TPC gas
	arrow6->SetLineColor(kWhite);
	arrow6->SetFillColor(kWhite);
	arrow6->SetFillStyle(1001);
	arrow6->SetLineWidth(2.);
	arrow6->Draw();
	
}

void DrawAutoGammaHisto2D(	TH2 *histo,  
						TString Title, TString XTitle, TString YTitle, TString Input,
						Bool_t YRange, Float_t YMin ,Float_t YMax, 
						Bool_t XRange, Float_t XMin, Float_t XMax,Float_t titleOffsetX=1.4, Float_t titleOffsetY=1.2) {
	
	
	if (YRange && XRange){
		histo->GetYaxis()->SetRangeUser(YMin, YMax);	
		histo->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	if ( !YRange && XRange){
		histo->GetXaxis()->SetRangeUser(XMin, XMax);	
	}
	
	if (YRange && !XRange){
		histo->GetYaxis()->SetRangeUser(YMin, YMax);
	}
	
	if(Title.CompareTo("") != 0){
		histo->SetTitle(Title.Data());
	}
	if(XTitle.CompareTo("") != 0){
		histo->SetXTitle(XTitle.Data());
	}
	if(YTitle.CompareTo("") != 0){
		histo->SetYTitle(YTitle.Data());
	}
	histo->GetYaxis()->SetTitleSize(0.043);	
	histo->GetYaxis()->SetLabelSize(0.035);
	histo->GetXaxis()->SetLabelSize(0.035);
	histo->GetYaxis()->SetDecimals();
	histo->GetYaxis()->SetTitleOffset(titleOffsetY);
	histo->GetXaxis()->SetTitleOffset(titleOffsetX);
	histo->GetXaxis()->SetTitleSize(0.043);	
	histo->GetXaxis()->SetTitleOffset(0.88);
	histo->DrawCopy("colz");
	if(Input.CompareTo("") != 0){
		TLegend* leg2 = new TLegend(0.6,0.82,0.83,0.9);
		leg2->SetTextSize(0.04);			
		leg2->SetFillColor(0);
		leg2->AddEntry(histo,(Input.Data()));
		leg2->Draw("same");
	}
}



void StyleSettingsThesis(){
	//gStyle->SetOptTitle(kFALSE);
	gStyle->SetOptDate(0);   //show day and time
	gStyle->SetOptStat(0);  //show statistic
	gStyle->SetPalette(1,0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetTextSize(0.5);
	gStyle->SetLabelSize(0.03,"xyz");
	gStyle->SetLabelOffset(0.002,"xyz");
	gStyle->SetTitleFontSize(0.04);
	gStyle->SetTitleOffset(1,"y");
	gStyle->SetTitleOffset(0.7,"x");		
	gStyle->SetCanvasColor(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetLineWidth(0.01);
	
	gStyle->SetPadTopMargin(0.03);
	gStyle->SetPadBottomMargin(0.09);
	gStyle->SetPadRightMargin(0.03);
	gStyle->SetPadLeftMargin(0.13);
	
	
	TGaxis::SetMaxDigits(3);
	gErrorIgnoreLevel=kError;
}

void SetPlotStyle() {
// 	const Int_t nRGBs = 7;
	const Int_t nRGBs = 5;
	const Int_t nCont = 255;
	
	Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[nRGBs] = { 0.31, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[nRGBs]  = { 0.51, 1., 0.12, 0.00, 0.00};

// 	Double_t stops[nRGBs] = {  0.34, 0.61, 0.84, 1.00 };
// 	Double_t red[nRGBs]   = {  0.00, 0.87, 1.00, 0.51 };
// 	Double_t green[nRGBs] = {  0.81, 1.00, 0.20, 0.00 };
// // 	Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
// 	Double_t blue[nRGBs]  = {  1., 0.00, 0.00, 0.00};

// 	Double_t blue[nRGBs]  = { 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
	TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
	gStyle->SetNumberContours(nCont);
}

TString ReturnFullCollisionsSystem(TString fEnergyFlagOpt){ 
   if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
      return  "pp, #sqrt{#it{s}} = 7 TeV";
   } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
      return  "pp, #sqrt{#it{s}} = 900 GeV";
   } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
      return  "pp, #sqrt{#it{s}} = 2.76 TeV";
   } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
      return "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
   } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
      return "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
   } else {
      cout << "No correct collision system specification, has been given" << endl;
      return "";
   }
}

TString ReturnDateString(){
	TDatime 	today;
	int 		iDate = today.GetDate();
	int 		iYear=iDate/10000;
	int 		iMonth=(iDate%10000)/100;
	int 		iDay=iDate%100;
	TString 	cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
		    "Jul","Aug","Sep","Oct","Nov","Dec"};
	TString 	textDayth;
	if (iDay== 11){
		textDayth = "st";
	} else if  (iDay== 12){
		textDayth = "th";
	} else if  (iDay== 13){
		textDayth = "th";
	} else if  (iDay%10 == 1){
		textDayth = "st";
	} else if (iDay%10 == 2){
		textDayth = "nd";
	} else if (iDay%10 == 3){
		textDayth = "rd";
	} else {
		textDayth = "th";
	}
	return Form("%i^{%s} %s %i",iDay, textDayth.Data(),cMonth[iMonth-1].Data(), iYear);
}

void  Material2DXY(TString data = "myOutput", TString cutSel = "", TString suffix = "gif", TString optEnergy=""){        
        
	gROOT->Reset();
	gROOT->SetStyle("Plain");
 	
	collisionSystem = ReturnFullCollisionsSystem(optEnergy);
	if (collisionSystem.CompareTo("") == 0){
		cout << "No correct collision system specification, has been given" << endl;
		return;
	}
   
	TString outputDirectory=	 		"MaterialPlots";		
	gSystem->Exec("mkdir -p "+outputDirectory);
	
	     
	StyleSettingsThesis();
	SetPlotStyle();

	// Which file do you want to analyse
	TString 	nameDataFile = 		(Form("%s",data.Data()));
			
	//How big should the right margin in 2D plots be? 
	Float_t 	rightMargin = 			0.12;
	Float_t 	leftMargin = 			0.11;
	
	// Get the histos
	TFile 		fileData(nameDataFile);  
	
	// choice of dateset
	TString nameGammaDirectory;
	TString cutSelMat = cutSel(0,4);
	cout << cutSelMat.Data() << endl;
	if(cutSelMat.CompareTo("") != 0){
		nameGammaDirectory = 		Form("GammaConv_%s",  cutSelMat.Data());
		cout << nameGammaDirectory.Data() << endl;
	}       
	
	textDate = ReturnDateString();
	
	TLatex *latexEtaRange= 	new TLatex(0.15,0.92,"|#eta| < 0.9 "); 
	latexEtaRange->SetNDC();
	latexEtaRange->SetTextFont(62);
	latexEtaRange->SetTextSize(0.04);
	latexEtaRange->SetLineWidth(6);      
	
	//------------------------------- Reading FILES ----------------------------------------------------------------------
	TDirectory 	*directoryGammaConvData = 	new TDirectory(); 

	directoryGammaConvData = 			(TDirectory*)fileData.Get(nameGammaDirectory); 
	TH2F *histoConversionPointVsXYData=				(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_XY");
	
	TH2F *histoConversionPointVsXYPlotData =  new TH2F("ESD_ConversionMapping_XYPlot","",2400,-120.,120.,2400,-120.,120.);
	for (Int_t iX = 1 ; iX < histoConversionPointVsXYData->GetNbinsX(); iX ++){
		for (Int_t iY = 1; iY < histoConversionPointVsXYData->GetNbinsY(); iY ++) {
// 			Double_t radius = TMath::Sqrt(histoConversionPointVsXYData->GetXaxis()->GetBinCenter(iX)* histoConversionPointVsXYData->GetXaxis()->GetBinCenter(iX) + histoConversionPointVsXYData->GetYaxis()->GetBinCenter(iY)* histoConversionPointVsXYData->GetYaxis()->GetBinCenter(iY));
			histoConversionPointVsXYPlotData->Fill(histoConversionPointVsXYData->GetXaxis()->GetBinCenter(iX), histoConversionPointVsXYData->GetYaxis()->GetBinCenter(iY),  histoConversionPointVsXYData->GetBinContent(iX,iY));
		}
	}
	
	Int_t sizeCanvasX = 500;
	Int_t sizeCanvasY = 440;
	if (suffix.CompareTo("eps")!=0){
		sizeCanvasX = 1000;
		sizeCanvasY = 880;
	}
	
	TCanvas * canvas2DimensionXY = new TCanvas("canvas2DimensionXY","",sizeCanvasX,sizeCanvasY);  // gives the page size
	canvas2DimensionXY->SetLogz(1);
	canvas2DimensionXY->SetTopMargin(0.02);                                              
	canvas2DimensionXY->SetRightMargin(rightMargin);                                            
	canvas2DimensionXY->SetLeftMargin(leftMargin);
	canvas2DimensionXY->SetBottomMargin(0.08);
	canvas2DimensionXY->cd();
	histoConversionPointVsXYPlotData->GetZaxis()->SetRangeUser(1,histoConversionPointVsXYPlotData->GetMaximum());
	DrawAutoGammaHisto2D(   histoConversionPointVsXYPlotData,
									"", "X (cm)", "Y (cm)", "",
									kTRUE, -250., 250.,
									kTRUE, -250., 250.);
	DrawStructureNew();
	latexEtaRange->SetTextColor(kWhite);
	latexEtaRange->Draw();

	canvas2DimensionXY->Update();
	canvas2DimensionXY->SaveAs(Form("%s/XY_distribution.%s",outputDirectory.Data(),suffix.Data()));
	delete canvas2DimensionXY;
						    
}
