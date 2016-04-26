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

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;
const Double_t kMean=0.136 ; //Approximate peak position to facilitate error estimate

//-----------------------------------------------------------------------------
Double_t CB(Double_t * x, Double_t * par){
  //Parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean) ;
}
//-----------------------------------------------------------------------------
Double_t BG1(Double_t * x, Double_t * par){
  //Normalizatino of Mixed
  return par[0]+par[1]*(x[0]-kMean) ;
}


void PlotInvMassPHOSAndPCMTogether(){

	TString	date = ReturnDateString();
	TString dateForOutput = ReturnDateStringForOutput();
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();

	TString outputDir = Form("eps/%s/CombineMesonMeasurementsPbPbX",dateForOutput.Data());
	
	gSystem->Exec("mkdir -p "+outputDir);

	// Loading & preparing PCM histograms
	Int_t binPbPb = 3;
	Int_t binpp =3;
	Double_t ptMinPCM = 0.8;
	Double_t ptMaxPCM = 1.;
	
	TFile * filePCMPbPb0010 = new TFile("501000103209297002322000000_01523045009000/PbPb_2.76TeV/Pi0_data_GammaConvV1WithoutCorrection_501000103209297002322000000_01523045009000.root") ;
	TH1D* histoPCMSignalPlusBG0010 = (TH1D*)filePCMPbPb0010->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPCMSignal0010 = (TH1D*)filePCMPbPb0010->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPCMRemainingBG0010 = (TH1D*)histoPCMSignal0010->Clone("histoPCMRemainingBG0010");
	TF1* fitPCMSignal0010 = (TF1*)filePCMPbPb0010->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",binPbPb));
	histoPCMSignal0010->Fit(fitPCMSignal0010,"QRME0");
	for (Int_t i=0; i < 6; i++){
		cout << fitPCMSignal0010->GetParameter(i) << "\t +- " << fitPCMSignal0010->GetParError(i) << endl;
	} 	
	TF1*  fitLinearBck0010 = new TF1("Linear0010","[0]+[1]*x",0.0,0.3);
	fitLinearBck0010->SetParameter(0, fitPCMSignal0010->GetParameter(4));
	fitLinearBck0010->SetParameter(1, fitPCMSignal0010->GetParameter(5));
	TVirtualFitter * fitter0010 = TVirtualFitter::GetFitter();
	Int_t nFreePar0010 = fitPCMSignal0010->GetNumberFreeParameters();
	double * covMatrix0010 = fitter0010->GetCovarianceMatrix();
	for (Int_t i = 1; i < histoPCMSignal0010->GetXaxis()->FindBin(0.3); i++){
		Double_t startBinEdge = histoPCMSignal0010->GetXaxis()->GetBinLowEdge(i);
		Double_t endBinEdge = histoPCMSignal0010->GetXaxis()->GetBinUpEdge(i);
		Double_t intLinearBack = fitLinearBck0010->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
		Double_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fitPCMSignal0010->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPCMSignal0010->GetParError(5),2)+2*covMatrix0010[nFreePar0010*nFreePar0010-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
		histoPCMRemainingBG0010->SetBinContent(i,intLinearBack);
		histoPCMRemainingBG0010->SetBinError(i,errorLinearBck);
		cout << fitLinearBck0010->Eval(startBinEdge) << "\t" <<fitLinearBck0010->Eval(endBinEdge) << "\t" <<histoPCMRemainingBG0010->GetBinContent(i) << "\t" <<histoPCMSignal0010->GetBinContent(i) << endl;
	}
	histoPCMSignal0010->Add(histoPCMRemainingBG0010,-1.);
	fitPCMSignal0010->SetParameter(4,0.);
	fitPCMSignal0010->SetParameter(5,0.);
	
	TFile * filePCMPbPb6080 = new TFile("568000103209297002322000000_01523065009000/PbPb_2.76TeV/Pi0_data_GammaConvV1WithoutCorrection_568000103209297002322000000_01523065009000.root") ;
	TH1D* histoPCMSignalPlusBG6080 = (TH1D*)filePCMPbPb6080->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPCMSignal6080 = (TH1D*)filePCMPbPb6080->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPCMRemainingBG6080 = (TH1D*)histoPCMSignal6080->Clone("histoPCMRemainingBG0010");
	TF1* fitPCMSignal6080 = (TF1*)filePCMPbPb6080->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",binPbPb));
	histoPCMSignal6080->Fit(fitPCMSignal6080,"QRME0");
	for (Int_t i=0; i < 6; i++){
		cout << fitPCMSignal6080->GetParameter(i) << "\t +- " << fitPCMSignal6080->GetParError(i) << endl;
	} 	
	TF1*  fitLinearBck6080 = new TF1("Linear6080","[0]+[1]*x",0.0,0.3);
	fitLinearBck6080->SetParameter(0, fitPCMSignal6080->GetParameter(4));
	fitLinearBck6080->SetParameter(1, fitPCMSignal6080->GetParameter(5));
	TVirtualFitter * fitter6080 = TVirtualFitter::GetFitter();
	Int_t nFreePar6080 = fitPCMSignal6080->GetNumberFreeParameters();
	double * covMatrix6080 = fitter6080->GetCovarianceMatrix();
	for (Int_t i = 1; i < histoPCMSignal6080->GetXaxis()->FindBin(0.3); i++){
		
		Double_t startBinEdge = histoPCMSignal6080->GetXaxis()->GetBinLowEdge(i);
		Double_t endBinEdge = histoPCMSignal6080->GetXaxis()->GetBinUpEdge(i);
		Double_t intLinearBack = fitLinearBck6080->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
		Double_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fitPCMSignal6080->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPCMSignal6080->GetParError(5),2)+2*covMatrix6080[nFreePar6080*nFreePar6080-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
		histoPCMRemainingBG6080->SetBinContent(i,intLinearBack);
		histoPCMRemainingBG6080->SetBinError(i,errorLinearBck);
		cout << fitLinearBck6080->Eval(startBinEdge) << "\t" <<fitLinearBck6080->Eval(endBinEdge) << "\t" <<histoPCMRemainingBG6080->GetBinContent(i) << "\t" <<histoPCMSignal6080->GetBinContent(i) << endl;
	}
	histoPCMSignal6080->Add(histoPCMRemainingBG6080,-1.);
	fitPCMSignal6080->SetParameter(4,0.);
	fitPCMSignal6080->SetParameter(5,0.);
	
	TFile * filePCMpp = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Pi0_data_GammaConvV1WithoutCorrection_000001100209366300380000000_01631031009000.root") ;
	TH1D* histoPCMSignalPlusBGpp = (TH1D*)filePCMpp->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",binpp));
	TH1D* histoPCMSignalpp = (TH1D*)filePCMpp->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",binpp));
	TH1D* histoPCMRemainingBGpp = (TH1D*)histoPCMSignalpp->Clone("histoPCMRemainingBG0010");
	TF1* fitPCMSignalpp = (TF1*)filePCMpp->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",binpp));
	histoPCMSignalpp->Fit(fitPCMSignalpp,"QRME0");
	for (Int_t i=0; i < 6; i++){
		cout << fitPCMSignalpp->GetParameter(i) << "\t +- " << fitPCMSignalpp->GetParError(i) << endl;
	} 	
	TF1*  fitLinearBckpp = new TF1("Linearpp","[0]+[1]*x",0.0,0.3);
	fitLinearBckpp->SetParameter(0, fitPCMSignalpp->GetParameter(4));
	fitLinearBckpp->SetParameter(1, fitPCMSignalpp->GetParameter(5));
	TVirtualFitter * fitterpp = TVirtualFitter::GetFitter();
	Int_t nFreeParpp = fitPCMSignalpp->GetNumberFreeParameters();
	double * covMatrixpp = fitterpp->GetCovarianceMatrix();
	for (Int_t i = 1; i < histoPCMSignalpp->GetXaxis()->FindBin(0.3); i++){
		
		Double_t startBinEdge = histoPCMSignalpp->GetXaxis()->GetBinLowEdge(i);
		Double_t endBinEdge = histoPCMSignalpp->GetXaxis()->GetBinUpEdge(i);
		Double_t intLinearBack = fitLinearBckpp->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
		Double_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fitPCMSignalpp->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPCMSignalpp->GetParError(5),2)+2*covMatrixpp[nFreeParpp*nFreeParpp-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
		histoPCMRemainingBGpp->SetBinContent(i,intLinearBack);
		histoPCMRemainingBGpp->SetBinError(i,errorLinearBck);
		cout << fitLinearBckpp->Eval(startBinEdge) << "\t" <<fitLinearBckpp->Eval(endBinEdge) << "\t" <<histoPCMRemainingBGpp->GetBinContent(i) << "\t" <<histoPCMSignalpp->GetBinContent(i) << endl;
	}
	histoPCMSignalpp->Add(histoPCMRemainingBGpp,-1.);
	fitPCMSignalpp->SetParameter(4,0.);
	fitPCMSignalpp->SetParameter(5,0.);
	
	// Loading & preparing PHOS histograms
	TFile * filePHOSPbPb = new TFile("ExternalInputPbPb/PHOS/PHOS_PbPb_mgg.root") ;
	TH2F *histo2DMass0010 = (TH2F*)filePHOSPbPb->Get("hPi0Both2core_cen0") ;
	histo2DMass0010->Add((TH2F*)filePHOSPbPb->Get("hPi0Both2core_cen1")) ;
	TH2F *histo2DMixedEvent0010 = (TH2F*)filePHOSPbPb->Get("hMiPi0Both2core_cen0") ;
	histo2DMixedEvent0010->Add((TH2F*)filePHOSPbPb->Get("hMiPi0Both2core_cen1")) ;
	TH2F *histo2DMass6080 = (TH2F*)filePHOSPbPb->Get("hPi0Both2core_cen5") ;
	TH2F *histo2DMixedEvent6080= (TH2F*)filePHOSPbPb->Get("hMiPi0Both2core_cen5") ;

	Double_t ptMinPHOS = 2.;
	Double_t ptMaxPHOS = 2.5;
	
	TAxis * ptaxis = histo2DMass0010->GetYaxis() ;
	Int_t iminPHOS=ptaxis->FindBin(ptMinPHOS+0.0001);
	Int_t imaxPHOS=ptaxis->FindBin(ptMaxPHOS-0.0001) ;
	
	TH1D * histoPHOSSignalPlusBG0010 = (TH1D*)histo2DMass0010->ProjectionX("histoPHOSSignalPlusBG0010",iminPHOS,imaxPHOS) ;
	histoPHOSSignalPlusBG0010->Sumw2() ;
	TH1D * histoPHOSSignalPlusBG6080 = (TH1D*)histo2DMass6080->ProjectionX("histoPHOSSignalPlusBG6080",iminPHOS,imaxPHOS) ;
	histoPHOSSignalPlusBG6080->Sumw2() ;
	
	TH1D * histoBGMixedEvent0010= (TH1D*)histo2DMixedEvent0010->ProjectionX("histoBGMixedEvent0010",iminPHOS,imaxPHOS) ;
	TH1D * histoBGMixedEvent6080= (TH1D*)histo2DMixedEvent6080->ProjectionX("histoBGMixedEvent6080",iminPHOS,imaxPHOS) ;


	TF1 * fit1_0010 = new TF1("fit0010",CB,0.,1.,6) ;
	TF1 * fbg1_0010 = new TF1("bg0010",BG1,0.,1.,3) ;
	
	TH1D * remi_0010 = (TH1D*)histoPHOSSignalPlusBG0010->Clone("remi_0010") ;
	remi_0010->Divide(histoBGMixedEvent0010) ;
	remi_0010->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})") ;
	remi_0010->SetYTitle("Real/Mixed") ;
	remi_0010->GetXaxis()->SetRangeUser(0.,0.3) ;
	
	Double_t rangeMin=0.1 ;
	Double_t rangeMax=0.2 ;
	fit1_0010->SetParameters(0.001,0.136,0.0055,0.0002,-0.002,0.0) ;
	fit1_0010->SetParLimits(0,0.000,1.000) ;
	fit1_0010->SetParLimits(1,0.120,0.145) ;
	fit1_0010->SetParLimits(2,0.005,0.012) ;
	remi_0010->Fit(fit1_0010,"Q" ,"",rangeMin,rangeMax) ;
	remi_0010->Fit(fit1_0010,"MQ","",rangeMin,rangeMax) ;  
	
	fbg1_0010->SetParameters(fit1_0010->GetParameter(3),fit1_0010->GetParameter(4),fit1_0010->GetParameter(5)); 
	histoBGMixedEvent0010 ->Multiply(fbg1_0010) ;
	TH1D * histoPHOSSignal0010 = (TH1D*)histoPHOSSignalPlusBG0010->Clone("Signal") ;
	histoPHOSSignal0010->Add(histoBGMixedEvent0010,-1.);

	TF1 * fitPHOSSignal0010  = new TF1("gaus","gaus",0.,1.) ;
	fitPHOSSignal0010->SetLineColor(kBlue);
	fitPHOSSignal0010->SetLineWidth(2);
	histoPHOSSignal0010->Fit(fitPHOSSignal0010,"0","",0.05,0.3) ;

	
	TF1 * fit1_6080 = new TF1("fit6080",CB,0.,1.,6) ;
	TF1 * fbg1_6080 = new TF1("bg6080",BG1,0.,1.,3) ;
	
	TH1D * remi_6080 = (TH1D*)histoPHOSSignalPlusBG6080->Clone("remi_6080") ;
	remi_6080->Divide(histoBGMixedEvent6080) ;
	remi_6080->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})") ;
	remi_6080->SetYTitle("Real/Mixed") ;
	remi_6080->GetXaxis()->SetRangeUser(0.,0.3) ;
	
	fit1_6080->SetParameters(0.001,0.136,0.0055,0.0002,-0.002,0.0) ;
	fit1_6080->SetParLimits(0,0.000,1.000) ;
	fit1_6080->SetParLimits(1,0.120,0.145) ;
	fit1_6080->SetParLimits(2,0.005,0.012) ;
	remi_6080->Fit(fit1_6080,"Q" ,"",rangeMin,rangeMax) ;
	remi_6080->Fit(fit1_6080,"MQ","",rangeMin,rangeMax) ;  
	
	fbg1_6080->SetParameters(fit1_6080->GetParameter(3),fit1_6080->GetParameter(4),fit1_6080->GetParameter(5)); 
	histoBGMixedEvent6080 ->Multiply(fbg1_6080) ;
	TH1D * histoPHOSSignal6080 = (TH1D*)histoPHOSSignalPlusBG6080->Clone("Signal") ;
	histoPHOSSignal6080->Add(histoBGMixedEvent6080,-1.);
	
	TF1 * fitPHOSSignal6080  = new TF1("gaus","gaus",0.,1.);
	fitPHOSSignal6080->SetLineColor(kBlue);
	fitPHOSSignal6080->SetLineWidth(2);
	histoPHOSSignal6080->Fit(fitPHOSSignal6080,"0","",0.05,0.3) ;
	
	TFile* filePHOSpp = new TFile("ExternalInput/PHOS/2.76TeV/PHOS_pp2760_pi0_FitResult.root") ;
	TH1F* histoPHOSSignalPlusBGpp = (TH1F*)filePHOSpp->Get("Real");
	TH1F* histoPHOSSignalpp = (TH1F*)filePHOSpp->Get("hp2");
  
	
	
	Double_t arrayBoundsXIndMeasRatio[4];
	Double_t arrayBoundsYIndMeasRatio[3];
	Double_t relativeMarginsIndMeasRatioX[3];
	Double_t relativeMarginsIndMeasRatioY[3];
	ReturnCorrectValuesForCanvasScaling(1200,800, 3, 2,0.0, 0.00, 0.00,0.0,arrayBoundsXIndMeasRatio,arrayBoundsYIndMeasRatio,relativeMarginsIndMeasRatioX,relativeMarginsIndMeasRatioY);
	
	Double_t marginUp = 0.06;
	Double_t marginDown = 0.14;
	Double_t marginLeft = 0.14;
	Double_t marginRight = 0.01;
	Double_t scaleFactorPi0PCM = 15; 
	Double_t scaleFactorPi0PHOS = 5; 
	
	TCanvas * canvas6PartRatioIndMeas = new TCanvas("canvas6PartRatioIndMeas","",10,10,1200,800);  // gives the page size		
	canvas6PartRatioIndMeas->cd();

	TPad* pad6PartRatioIndMeas1 = new TPad("pad6PartRatioIndMeas1", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[1],arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas1, marginLeft,marginRight,marginUp,marginDown);
	pad6PartRatioIndMeas1->Draw();
	TPad* pad6PartRatioIndMeas2 = new TPad("pad6PartRatioIndMeas2", "", arrayBoundsXIndMeasRatio[0], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas2, marginLeft,marginRight,marginUp,marginDown);
	pad6PartRatioIndMeas2->Draw();
	
	TPad* pad6PartRatioIndMeas3 = new TPad("pad6PartRatioIndMeas3", "", arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas3, marginLeft,marginRight,marginUp,marginDown);
	pad6PartRatioIndMeas3->Draw();
	TPad* pad6PartRatioIndMeas4 = new TPad("pad6PartRatioIndMeas4", "", arrayBoundsXIndMeasRatio[1], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas4, marginLeft,marginRight,marginUp,marginDown);
	pad6PartRatioIndMeas4->Draw();
	
	TPad* pad6PartRatioIndMeas5 = new TPad("pad6PartRatioIndMeas5", "", arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas5, marginLeft,marginRight,marginUp,marginDown);
	pad6PartRatioIndMeas5->Draw();
	TPad* pad6PartRatioIndMeas6 = new TPad("pad6PartRatioIndMeas6", "", arrayBoundsXIndMeasRatio[2], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas6, marginLeft,marginRight,marginUp,marginDown); 
	pad6PartRatioIndMeas6->Draw();

	Int_t textSizeLabelsPixelRatio = 25;
	Double_t marginRatio = 0.145*1200;
	Double_t textSizeLabels = 0;
	Double_t textSizeFac = 0;
	
	if (pad6PartRatioIndMeas1->XtoPixel(pad6PartRatioIndMeas1->GetX2()) < pad6PartRatioIndMeas1->YtoPixel(pad6PartRatioIndMeas1->GetY1())){
	    textSizeLabels = (Double_t)textSizeLabelsPixelRatio/pad6PartRatioIndMeas1->XtoPixel(pad6PartRatioIndMeas1->GetX2()) ;
	    textSizeFac = (Double_t)1./pad6PartRatioIndMeas1->XtoPixel(pad6PartRatioIndMeas1->GetX2()) ;
	} else {
	  textSizeLabels = (Double_t)textSizeLabelsPixelRatio/pad6PartRatioIndMeas1->YtoPixel(pad6PartRatioIndMeas1->GetY1());
	  textSizeFac = (Double_t)1./pad6PartRatioIndMeas1->YtoPixel(pad6PartRatioIndMeas1->GetY1());
	}
	
	Double_t minInvMass = 0.095;
	Double_t maxInvMass = 0.205;
	Size_t markerSizeInvMass = 1;
	Width_t widthCommonFit = 1;
	
	pad6PartRatioIndMeas1->cd();
	
	SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBGpp, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.85*textSizeLabels,textSizeLabels,  0.85*textSizeLabels,textSizeLabels, 0.95,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPCMSignalPlusBGpp->GetYaxis()->SetTitleFont(62);
	histoPCMSignalPlusBGpp->GetXaxis()->SetTitleFont(62);
	histoPCMSignalPlusBGpp->GetYaxis()->SetLabelFont(42);
	histoPCMSignalPlusBGpp->GetXaxis()->SetLabelFont(42);
	histoPCMSignalPlusBGpp->SetLineWidth(0.8);
	histoPCMSignalPlusBGpp->GetYaxis()->SetRangeUser(-20,histoPCMSignalPlusBGpp->GetMaximum());
	histoPCMSignalPlusBGpp->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	histoPCMSignalPlusBGpp->Draw("hist,e");
	
	DrawGammaSetMarker(histoPCMSignalpp,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoPCMSignalpp->Draw("same,pe");
	DrawGammaSetMarkerTF1( fitPCMSignalpp, 1, widthCommonFit, kBlue+2);
	fitPCMSignalpp->Draw("same");
	
	TLatex *labelPP = new TLatex(0.5,0.80,Form("pp #sqrt{#it{s}} = 2.76 TeV"));
	SetStyleTLatex( labelPP, 0.85*textSizeLabels,4,kBlack);
	labelPP->Draw();
	TLatex *labelMethod1 = new TLatex(0.18,0.87,Form("PCM"));
	SetStyleTLatex( labelMethod1, 0.85*textSizeLabels,4,kBlack);
	labelMethod1->Draw();
	TLatex *labelRangePPPCM = new TLatex(0.5,0.87,Form("0.8 < #it{p}^{#gamma#gamma}_{T}< 1.0 GeV/c"));
	SetStyleTLatex( labelRangePPPCM, 0.85*textSizeLabels,4,kBlack);
	labelRangePPPCM->Draw();
	
	TLegend* legendInvMass = new TLegend(0.6,0.4,0.9,0.61);
	legendInvMass->SetFillColor(0);
	legendInvMass->SetLineColor(0);
	legendInvMass->SetTextSize(0.85*textSizeLabels);
	legendInvMass->SetTextFont(42);
// 	legendInvMass->SetMargin(0.14);
	legendInvMass->AddEntry(histoPCMSignalPlusBGpp,"signal + bkg","l,e");
	legendInvMass->AddEntry(histoPCMSignalpp,"signal","p");
	legendInvMass->AddEntry(fitPCMSignalpp,"fit","l");
	legendInvMass->Draw();

	
	
	pad6PartRatioIndMeas2->cd();
	SetStyleHistoTH1ForGraphs(histoPHOSSignalPlusBGpp, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.85*textSizeLabels,textSizeLabels,  0.85*textSizeLabels,textSizeLabels, 0.95,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPHOSSignalPlusBGpp->GetYaxis()->SetTitleFont(62);
	histoPHOSSignalPlusBGpp->GetXaxis()->SetTitleFont(62);	
	histoPHOSSignalPlusBGpp->GetYaxis()->SetLabelFont(42);
	histoPHOSSignalPlusBGpp->GetXaxis()->SetLabelFont(42);
	DrawGammaSetMarker(histoPHOSSignalPlusBGpp,20,0, kBlack , kBlack);  
	histoPHOSSignalPlusBGpp->SetLineWidth(0.8);
	histoPHOSSignalPlusBGpp->GetYaxis()->SetRangeUser(-20,0.55*histoPHOSSignalPlusBGpp->GetMaximum());
	histoPHOSSignalPlusBGpp->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	histoPHOSSignalPlusBGpp->Draw("hist,e");
	DrawGammaSetMarker(histoPHOSSignalpp,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoPHOSSignalpp->SetLineWidth(0.8);
	histoPHOSSignalpp->Draw("same,pe");
	histoPHOSSignalpp->GetFunction("gs")->SetLineColor(kBlue+2);
	histoPHOSSignalpp->GetFunction("gs")->SetLineWidth(widthCommonFit);
	histoPHOSSignalpp->GetFunction("gs")->Draw("same");
	
	labelPP->Draw();
	TLatex *labelMethod2 = new TLatex(0.18,0.87,Form("PHOS"));
	SetStyleTLatex( labelMethod2, 0.85*textSizeLabels,4,kBlack);
	labelMethod2->Draw();
	TLatex *labelRangePPPHOS = new TLatex(0.5,0.87,Form("2.0 < #it{p}^{#gamma#gamma}_{T}< 2.5 GeV/c"));
	SetStyleTLatex( labelRangePPPHOS, 0.85*textSizeLabels,4,kBlack);
	labelRangePPPHOS->Draw();

	
	
	pad6PartRatioIndMeas3->cd();
	SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBG6080, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.85*textSizeLabels,textSizeLabels,  0.85*textSizeLabels,textSizeLabels, 0.95,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPCMSignalPlusBG6080->GetYaxis()->SetTitleFont(62);
	histoPCMSignalPlusBG6080->GetXaxis()->SetTitleFont(62);
	histoPCMSignalPlusBG6080->GetYaxis()->SetLabelFont(42);
	histoPCMSignalPlusBG6080->GetXaxis()->SetLabelFont(42);
	histoPCMSignalPlusBG6080->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	histoPCMSignalPlusBG6080->GetYaxis()->SetRangeUser(-100,1.2*histoPCMSignalPlusBG6080->GetMaximum());
	histoPCMSignalPlusBG6080->Draw("hist,e");
	DrawGammaSetMarker(histoPCMSignal6080,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoPCMSignal6080->Draw("same,pe");
	DrawGammaSetMarkerTF1( fitPCMSignal6080, 1, widthCommonFit, kBlue+2);
	fitPCMSignal6080->Draw("same");
	TLatex *labelPbPb = new TLatex(0.45,0.80,Form("Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV"));
	SetStyleTLatex( labelPbPb, 0.85*textSizeLabels,4,kBlack);
	labelPbPb->Draw();
	TLatex *labelPbPb6080 = new TLatex(0.45,0.74,"60-80%");
	SetStyleTLatex( labelPbPb6080, 0.85*textSizeLabels,4,kBlack);
	labelPbPb6080->Draw();
	labelMethod1->Draw();
	TLatex *labelRangePbPbPCM = new TLatex(0.45,0.87,Form("0.8 < #it{p}^{#gamma#gamma}_{T}< 1.0 GeV/c"));
	SetStyleTLatex( labelRangePbPbPCM, 0.85*textSizeLabels,4,kBlack);
	labelRangePbPbPCM->Draw();

	pad6PartRatioIndMeas4->cd();
	SetStyleHistoTH1ForGraphs(histoPHOSSignalPlusBG6080, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.85*textSizeLabels,textSizeLabels,  0.85*textSizeLabels,textSizeLabels, 0.95,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPHOSSignalPlusBG6080->GetYaxis()->SetTitleFont(62);
	histoPHOSSignalPlusBG6080->GetXaxis()->SetTitleFont(62);
	histoPHOSSignalPlusBG6080->GetYaxis()->SetLabelFont(42);
	histoPHOSSignalPlusBG6080->GetXaxis()->SetLabelFont(42);
	DrawGammaSetMarker(histoPHOSSignalPlusBG6080,20,0, kBlack , kBlack);  
	histoPHOSSignalPlusBG6080->SetLineWidth(0.8);
	histoPHOSSignalPlusBG6080->GetYaxis()->SetRangeUser(-40,1.4*histoPHOSSignalPlusBG6080->GetMaximum());
	histoPHOSSignalPlusBG6080->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	histoPHOSSignalPlusBG6080->Draw("hist,e");
	DrawGammaSetMarker(histoPHOSSignal6080,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoPHOSSignal6080->SetLineWidth(0.8);
	histoPHOSSignal6080->Draw("same,pe");
	DrawGammaSetMarkerTF1( fitPHOSSignal6080, 1, widthCommonFit, kBlue+2);
	fitPHOSSignal6080->Draw("same");

	labelPbPb->Draw();
	labelPbPb6080->Draw();
	labelMethod2->Draw();
	TLatex *labelRangePbPbPHOS = new TLatex(0.45,0.87,Form("2.0 < #it{p}^{#gamma#gamma}_{T}< 2.5 GeV/c"));
	SetStyleTLatex( labelRangePbPbPHOS, 0.85*textSizeLabels,4,kBlack);
	labelRangePbPbPHOS->Draw();

	
	pad6PartRatioIndMeas5->cd();
	
	SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBG0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.85*textSizeLabels,textSizeLabels,  0.85*textSizeLabels,textSizeLabels, 0.95,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPCMSignalPlusBG0010->GetYaxis()->SetTitleFont(62);
	histoPCMSignalPlusBG0010->GetXaxis()->SetTitleFont(62);
	histoPCMSignalPlusBG0010->GetYaxis()->SetLabelFont(42);
	histoPCMSignalPlusBG0010->GetXaxis()->SetLabelFont(42);
	histoPCMSignalPlusBG0010->GetYaxis()->SetRangeUser(-10000,1.*histoPCMSignalPlusBG0010->GetMaximum());
	histoPCMSignalPlusBG0010->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	histoPCMSignalPlusBG0010->Draw("hist,e");
	DrawGammaSetMarker(histoPCMSignal0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoPCMSignal0010->Draw("same,pe");
	histoPCMSignal0010->Scale(scaleFactorPi0PCM);
	TH1D* histoFitPCMSignal0010 = (TH1D*)fitPCMSignal0010->GetHistogram();

	DrawGammaSetMarker(histoFitPCMSignal0010,20,0, kBlue+2 , kBlue+2);
	histoFitPCMSignal0010->SetLineWidth(widthCommonFit);
	histoFitPCMSignal0010->Scale(scaleFactorPi0PCM);
	histoFitPCMSignal0010->Draw("same,hist,c");

// 	TH1D* histoPCMRemainingBGFit0010 = (TH1D*)fitLinearBck0010->GetHistogram();
// 	histoPCMRemainingBGFit0010->Fit(fitLinearBck0010,"SQRME0");
	
// 	DrawGammaSetMarker(histoPCMRemainingBG0010,20,0, kGreen+2 , kGreen+2);
// 	histoPCMRemainingBG0010->Scale(scaleFactorPi0PCM);
// 	histoPCMRemainingBG0010->Draw("same,pe");
	
	TLatex *labelScalingPi0PCM = new TLatex(0.57,0.25,Form("Signal #times %1.0f",scaleFactorPi0PCM));
	SetStyleTLatex( labelScalingPi0PCM, 0.85*textSizeLabels,4,histoPCMSignal0010->GetLineColor());
	labelScalingPi0PCM->Draw();

	labelPbPb->Draw();
	TLatex *labelPbPb0010 = new TLatex(0.45,0.74,"0-10%");
	SetStyleTLatex( labelPbPb0010, 0.85*textSizeLabels,4,kBlack);
	labelPbPb0010->Draw();
	labelMethod1->Draw();
	labelRangePbPbPCM->Draw();

	pad6PartRatioIndMeas6->cd();
	SetStyleHistoTH1ForGraphs(histoPHOSSignalPlusBG0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.85*textSizeLabels,textSizeLabels,  0.85*textSizeLabels,textSizeLabels, 0.95,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPHOSSignalPlusBG0010->GetYaxis()->SetTitleFont(62);
	histoPHOSSignalPlusBG0010->GetXaxis()->SetTitleFont(62);
	histoPHOSSignalPlusBG0010->GetYaxis()->SetLabelFont(42);
	histoPHOSSignalPlusBG0010->GetXaxis()->SetLabelFont(42);
	DrawGammaSetMarker(histoPHOSSignalPlusBG0010,20,0, kBlack , kBlack);  
	histoPHOSSignalPlusBG0010->SetLineWidth(0.8);
	histoPHOSSignalPlusBG0010->GetYaxis()->SetRangeUser(-2000,1.4*histoPHOSSignalPlusBG0010->GetMaximum());
	histoPHOSSignalPlusBG0010->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	histoPHOSSignalPlusBG0010->Draw("hist,e");
	histoPHOSSignal0010->SetLineWidth(0.8);
	DrawGammaSetMarker(histoPHOSSignal0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoPHOSSignal0010->Draw("same,pe");
	histoPHOSSignal0010->Scale(scaleFactorPi0PHOS);
	TH1D* histoFitPHOSSignal0010 = (TH1D*)fitPHOSSignal0010->GetHistogram();
	DrawGammaSetMarker(histoFitPHOSSignal0010,20,0, kBlue+2 , kBlue+2);
	histoFitPHOSSignal0010->SetLineWidth(widthCommonFit);
	histoFitPHOSSignal0010->Scale(scaleFactorPi0PHOS);
	histoFitPHOSSignal0010->Draw("same,hist,c");
	TLatex *labelScalingPi0PHOS = new TLatex(0.57,0.35,Form("Signal #times %1.0f",scaleFactorPi0PHOS));
	SetStyleTLatex( labelScalingPi0PHOS, 0.85*textSizeLabels,4,histoPCMSignal0010->GetLineColor());
	labelScalingPi0PHOS->Draw();
	
	
	labelPbPb->Draw();
	labelPbPb0010->Draw();
	labelMethod2->Draw();
	labelRangePbPbPHOS->Draw();

	
	canvas6PartRatioIndMeas->SaveAs(Form("%s/InvMassPaper_%s.eps",outputDir.Data(),dateForOutput.Data()));	
}

