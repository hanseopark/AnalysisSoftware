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


void PlotInvMassEMCalAndPCMTogether(){

	TString	date = ReturnDateString();
	TString dateForOutput = ReturnDateStringForOutput();
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();

	TString outputDir = Form("pdf/%s/CombineMesonMeasurementsPbPbX",dateForOutput.Data());
	gSystem->Exec("mkdir -p "+outputDir);

	// Loading & preparing PCM histograms
	Int_t binPbPb = 4;
	Int_t binPbPbEMCal;
// 	if(meson.Contains("Pi0")){
// 		binPbPbPCM = 4;
// 		
// 	} else if(meson.Contains("Eta")){
// 		binPbPbPCM = 4;
// 	}

//===============================================================================================
	TFile * filePCMPbPb0010 = new TFile("50100013_00200009247602008250404000_0152501500000000/PbPb_2.76TeV/Pi0_data_GammaConvV1WithoutCorrection_50100013_00200009247602008250404000_0152501500000000.root");
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

//===============================================================================================	
	TFile * filePCMPbPb2050 = new TFile("52500013_00200009247602008250404000_0152501500000000/PbPb_2.76TeV/Pi0_data_GammaConvV1WithoutCorrection_52500013_00200009247602008250404000_0152501500000000.root");
	TH1D* histoPCMSignalPlusBG2050 = (TH1D*)filePCMPbPb2050->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPCMSignal2050 = (TH1D*)filePCMPbPb2050->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPCMRemainingBG2050 = (TH1D*)histoPCMSignal2050->Clone("histoPCMRemainingBG0010");
	TF1* fitPCMSignal2050 = (TF1*)filePCMPbPb2050->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",binPbPb));
	histoPCMSignal2050->Fit(fitPCMSignal2050,"QRME0");
	for (Int_t i=0; i < 6; i++){
		cout << fitPCMSignal2050->GetParameter(i) << "\t +- " << fitPCMSignal2050->GetParError(i) << endl;
	} 	
	TF1*  fitLinearBck2050 = new TF1("Linear2050","[0]+[1]*x",0.0,0.3);
	fitLinearBck2050->SetParameter(0, fitPCMSignal2050->GetParameter(4));
	fitLinearBck2050->SetParameter(1, fitPCMSignal2050->GetParameter(5));
	TVirtualFitter * fitter2050 = TVirtualFitter::GetFitter();
	Int_t nFreePar2050 = fitPCMSignal2050->GetNumberFreeParameters();
	double * covMatrix2050 = fitter2050->GetCovarianceMatrix();
	for (Int_t i = 1; i < histoPCMSignal2050->GetXaxis()->FindBin(0.3); i++){
		
		Double_t startBinEdge = histoPCMSignal2050->GetXaxis()->GetBinLowEdge(i);
		Double_t endBinEdge = histoPCMSignal2050->GetXaxis()->GetBinUpEdge(i);
		Double_t intLinearBack = fitLinearBck2050->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
		Double_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fitPCMSignal2050->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPCMSignal2050->GetParError(5),2)+2*covMatrix2050[nFreePar2050*nFreePar2050-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
		histoPCMRemainingBG2050->SetBinContent(i,intLinearBack);
		histoPCMRemainingBG2050->SetBinError(i,errorLinearBck);
		cout << fitLinearBck2050->Eval(startBinEdge) << "\t" <<fitLinearBck2050->Eval(endBinEdge) << "\t" <<histoPCMRemainingBG2050->GetBinContent(i) << "\t" <<histoPCMSignal2050->GetBinContent(i) << endl;
	}
	histoPCMSignal2050->Add(histoPCMRemainingBG2050,-1.);
	fitPCMSignal2050->SetParameter(4,0.);
	fitPCMSignal2050->SetParameter(5,0.);
	
//===============================================================================================
	TFile * filePCMPbPbEta0010 = new TFile("50100013_00200009247602008250404000_0152501500000000/PbPb_2.76TeV/Eta_data_GammaConvV1WithoutCorrection_50100013_00200009247602008250404000_0152501500000000.root");
	TH1D* histoPCMSignalPlusBGEta0010 = (TH1D*)filePCMPbPbEta0010->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPCMSignalEta0010 = (TH1D*)filePCMPbPbEta0010->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPCMRemainingBGEta0010 = (TH1D*)histoPCMSignalEta0010->Clone("histoPCMRemainingBGEta0010");
	TF1* fitPCMSignalEta0010 = (TF1*)filePCMPbPbEta0010->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",binPbPb));
	histoPCMSignalEta0010->Fit(fitPCMSignalEta0010,"QRME0");
	for (Int_t i=0; i < 6; i++){
		cout << fitPCMSignalEta0010->GetParameter(i) << "\t +- " << fitPCMSignalEta0010->GetParError(i) << endl;
	} 	
	TF1*  fitLinearBckEta0010 = new TF1("LinearEta0010","[0]+[1]*x",0.4,0.7);
	fitLinearBckEta0010->SetParameter(0, fitPCMSignalEta0010->GetParameter(4));
	fitLinearBckEta0010->SetParameter(1, fitPCMSignalEta0010->GetParameter(5));
	TVirtualFitter * fitterEta0010 = TVirtualFitter::GetFitter();
	Int_t nFreeParEta0010 = fitPCMSignalEta0010->GetNumberFreeParameters();
	double * covMatrixEta0010 = fitterEta0010->GetCovarianceMatrix();
	for (Int_t i = 1; i < histoPCMSignalEta0010->GetXaxis()->FindBin(0.7); i++){
		Double_t startBinEdge = histoPCMSignalEta0010->GetXaxis()->GetBinLowEdge(i);
		Double_t endBinEdge = histoPCMSignalEta0010->GetXaxis()->GetBinUpEdge(i);
		Double_t intLinearBack = fitLinearBckEta0010->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
		Double_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fitPCMSignalEta0010->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPCMSignalEta0010->GetParError(5),2)+2*covMatrixEta0010[nFreeParEta0010*nFreeParEta0010-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
		histoPCMRemainingBGEta0010->SetBinContent(i,intLinearBack);
		histoPCMRemainingBGEta0010->SetBinError(i,errorLinearBck);
		cout << fitLinearBckEta0010->Eval(startBinEdge) << "\t" <<fitLinearBckEta0010->Eval(endBinEdge) << "\t" <<histoPCMRemainingBGEta0010->GetBinContent(i) << "\t" <<histoPCMSignalEta0010->GetBinContent(i) << endl;
	}
	histoPCMSignalEta0010->Add(histoPCMRemainingBGEta0010,-1.);
	fitPCMSignalEta0010->SetParameter(4,0.);
	fitPCMSignalEta0010->SetParameter(5,0.);
	
//===============================================================================================
	TFile * filePCMPbPbEta2050 = new TFile("52500013_00200009247602008250404000_0152501500000000/PbPb_2.76TeV/Eta_data_GammaConvV1WithoutCorrection_52500013_00200009247602008250404000_0152501500000000.root");
	TH1D* histoPCMSignalPlusBGEta2050 = (TH1D*)filePCMPbPbEta2050->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPCMSignalEta2050 = (TH1D*)filePCMPbPbEta2050->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPCMRemainingBGEta2050 = (TH1D*)histoPCMSignalEta2050->Clone("histoPCMRemainingBGEta0010");
	TF1* fitPCMSignalEta2050 = (TF1*)filePCMPbPbEta2050->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",binPbPb));
	histoPCMSignalEta2050->Fit(fitPCMSignalEta2050,"QRME0");
	for (Int_t i=0; i < 6; i++){
		cout << fitPCMSignalEta2050->GetParameter(i) << "\t +- " << fitPCMSignalEta2050->GetParError(i) << endl;
	} 	
	TF1*  fitLinearBckEta2050 = new TF1("LinearEta2050","[0]+[1]*x",0.4,0.7);
	fitLinearBckEta2050->SetParameter(0, fitPCMSignalEta2050->GetParameter(4));
	fitLinearBckEta2050->SetParameter(1, fitPCMSignalEta2050->GetParameter(5));
	TVirtualFitter * fitterEta2050 = TVirtualFitter::GetFitter();
	Int_t nFreeParEta2050 = fitPCMSignalEta2050->GetNumberFreeParameters();
	double * covMatrixEta2050 = fitterEta2050->GetCovarianceMatrix();
	for (Int_t i = 1; i < histoPCMSignalEta2050->GetXaxis()->FindBin(0.7); i++){
		
		Double_t startBinEdge = histoPCMSignalEta2050->GetXaxis()->GetBinLowEdge(i);
		Double_t endBinEdge = histoPCMSignalEta2050->GetXaxis()->GetBinUpEdge(i);
		Double_t intLinearBack = fitLinearBckEta2050->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
		Double_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fitPCMSignalEta2050->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPCMSignalEta2050->GetParError(5),2)+2*covMatrixEta2050[nFreeParEta2050*nFreeParEta2050-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
		histoPCMRemainingBGEta2050->SetBinContent(i,intLinearBack);
		histoPCMRemainingBGEta2050->SetBinError(i,errorLinearBck);
		cout << fitLinearBckEta2050->Eval(startBinEdge) << "\t" <<fitLinearBckEta2050->Eval(endBinEdge) << "\t" <<histoPCMRemainingBGEta2050->GetBinContent(i) << "\t" <<histoPCMSignalEta2050->GetBinContent(i) << endl;
	}
	histoPCMSignalEta2050->Add(histoPCMRemainingBGEta2050,-1.);
	fitPCMSignalEta2050->SetParameter(4,0.);
	fitPCMSignalEta2050->SetParameter(5,0.);

//===============================================================================================
//===============================================================================================
	
	TFile *fEMC0      = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassEta0.root");
	TFile *fEMC2 = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassEta2.root");
    TFile *fEMCPI0 = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassPi0.root");
    TFile *fEMCPI2 = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassPi2.root");

	TH1F *histoEMCalSignalPlusBG0010;
	TF1 *fitEMCalSignal0010;
	TF1 *FitBack0010;
	TH1F *histoEMCalSignal0010;

	TF1 * fitpeak0010  = new TF1("gaus","gaus",0.08,0.22);
	DrawGammaSetMarkerTF1(fitpeak0010,1,1,kBlue+1);
	
	histoEMCalSignalPlusBG0010     = (TH1F*)fEMCPI0->Get("hDiff1D"); //hMassPi0
	histoEMCalSignalPlusBG0010->SetStats(0);
	fitEMCalSignal0010  	  = (TF1*)fEMCPI0->Get("fitPolyPi");	//fFitMassPi0
	FitBack0010  	  = (TF1*)fEMCPI0->Get("fitFullBackpi0");		//fFitBackPi0
	histoEMCalSignal0010 = (TH1F*)histoEMCalSignalPlusBG0010->Clone("histoEMCalSignal0010");
	histoEMCalSignal0010->Add(FitBack0010,-1.);
	histoEMCalSignal0010->Fit(fitpeak0010,"","",0.08,0.22);
	histoEMCalSignalPlusBG0010->SetMinimum(-20.);
//     fitpeak0010 = (TF1*)histoEMCalSignal0010->GetFunction("gaus");
	


//===============================================================================================

	TH1F *histoEMCalSignal2050;
	TH1F *histoEMCalSignalPlusBG2050;
	TF1 *fitEMCalSignal2050;
	TF1 *FitBack2050;

	TF1 * fitpeak2050  = new TF1("gaus","gaus",0.08,0.22);
	DrawGammaSetMarkerTF1(fitpeak2050,1,1,kBlue+1);

	histoEMCalSignalPlusBG2050     = (TH1F*)fEMCPI2->Get("hDiff1D");
	histoEMCalSignalPlusBG2050->SetStats(0);
	fitEMCalSignal2050  	  = (TF1*)fEMCPI2->Get("fitPolyPi");
	FitBack2050    	  = (TF1*)fEMCPI2->Get("fitFullBackpi0");
	histoEMCalSignal2050 = (TH1F*)histoEMCalSignalPlusBG2050->Clone("histoEMCalSignal2050");
	histoEMCalSignal2050->Add(FitBack2050,-1.);
	histoEMCalSignal2050->Fit(fitpeak2050,"","",0.08,0.22);
	histoEMCalSignalPlusBG2050->SetMinimum(-20.);
//     TF1 *fitpeak2050 = (TF1*)histoEMCalSignal2050->GetFunction("gaus");

//===============================================================================================
	
	TH1F *histoEMCalSignalPlusBGEta0010;
	TF1 *fitEMCalSignalEta0010;
	TF1 *FitBackEta0010;	
	TH1F *histoEMCalSignalEta0010;
	
	TF1 * fitpeakEta0010  = new TF1("gaus","gaus",0.45,0.7);
	DrawGammaSetMarkerTF1(fitpeakEta0010,1,1,kBlue+1);

	histoEMCalSignalPlusBGEta0010     = (TH1F*)fEMC0->Get("hDiff1DEta__1");
	histoEMCalSignalPlusBGEta0010->SetStats(0);

	fitEMCalSignalEta0010  	  = (TF1*)fEMC0->Get("fitPolyEta");
	FitBackEta0010  	  = (TF1*)fEMC0->Get("fitFullBackEta");
	histoEMCalSignalEta0010 = (TH1F*)histoEMCalSignalPlusBGEta0010->Clone("histoEMCalSignalEta0010");
	histoEMCalSignalEta0010->Add(FitBackEta0010,-1.);
	histoEMCalSignalEta0010->Fit(fitpeakEta0010,"","",0.4,0.75);
	histoEMCalSignalPlusBGEta0010->SetMinimum(-20.);
//     TF1 *fitpeakEta0010 = (TF1*)histoEMCalSignalEta0010->GetFunction("gaus");

//===============================================================================================

	TH1F *histoEMCalSignalPlusBGEta2050;
	TF1 *fitEMCalSignalEta2050;
	TF1 *FitBackEta2050;
	TH1F *histoEMCalSignalEta2050;

	TF1 * fitpeakEta2050  = new TF1("gaus","gaus",0.45,0.7);
	DrawGammaSetMarkerTF1(fitpeakEta2050,1,1,kBlue+1);
	
	histoEMCalSignalPlusBGEta2050     = (TH1F*)fEMC2->Get("hDiff1DEta__1");
	histoEMCalSignalPlusBGEta2050->SetStats(0);

	fitEMCalSignalEta2050  	  = (TF1*)fEMC2->Get("fitPolyEta");
	FitBackEta2050    	  = (TF1*)fEMC2->Get("fitFullBackEta");
	histoEMCalSignalEta2050 = (TH1F*)histoEMCalSignalPlusBGEta2050->Clone("histoEMCalSignalEta2050");
	histoEMCalSignalEta2050->Add(FitBackEta2050,-1.);
	histoEMCalSignalEta2050->Fit(fitpeakEta2050,"","",0.4,0.75);
	histoEMCalSignalPlusBGEta2050->SetMinimum(-20.);
//     TF1 *fitpeakEta2050 = (TF1*)histoEMCalSignalEta2050->GetFunction("gaus");


//===============================================================================================
//===============================================================================================


	Double_t ptMinEMCal = 12;
	Double_t ptMaxEMCal = 14;	
	
	
	Double_t marginUp = 0.06;
	Double_t marginDown = 0.14;
	Double_t marginLeft = 0.14;
	Double_t marginRight = 0.01;
	Double_t scaleFactorPi0PCM = 10; 
	Double_t scaleFactorEtaPCM = 40; 
	Double_t scaleFactorPi0EMCal = 1; 

	Double_t arrayBoundsXIndMeasRatio[5];
	Double_t arrayBoundsYIndMeasRatio[3];
	Double_t relativeMarginsIndMeasRatioX[4];
	Double_t relativeMarginsIndMeasRatioY[3];
	ReturnCorrectValuesForCanvasScaling(1200,800, 4, 2,0.0, 0.00, 0.00,0.0,arrayBoundsXIndMeasRatio,arrayBoundsYIndMeasRatio,relativeMarginsIndMeasRatioX,relativeMarginsIndMeasRatioY);
	
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

	TPad* pad6PartRatioIndMeas7 = new TPad("pad6PartRatioIndMeas7", "", arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[1], arrayBoundsXIndMeasRatio[4], arrayBoundsYIndMeasRatio[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas7, marginLeft,marginRight,marginUp,marginDown);
	pad6PartRatioIndMeas7->Draw();
	TPad* pad6PartRatioIndMeas8 = new TPad("pad6PartRatioIndMeas8", "", arrayBoundsXIndMeasRatio[3], arrayBoundsYIndMeasRatio[2], arrayBoundsXIndMeasRatio[4], arrayBoundsYIndMeasRatio[1],-1, -1, -2);
	DrawGammaPadSettings( pad6PartRatioIndMeas8, marginLeft,marginRight,marginUp,marginDown); 
	pad6PartRatioIndMeas8->Draw();

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
	
	
	
    TLatex *labelMethod1 = new TLatex(0.2,0.85,Form("PCM"));
    SetStyleTLatex( labelMethod1, 0.7*textSizeLabels,4,kBlack);
    TLatex *labelMethod2 = new TLatex(0.2,0.85,Form("EMCal"));
    SetStyleTLatex( labelMethod2, 0.7*textSizeLabels,4,kBlack);


	TLatex *labelPbPb = new TLatex(0.2,0.78,Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV"));
	SetStyleTLatex( labelPbPb, 0.7*textSizeLabels,4,kBlack);
	
	TLatex *labelScalingPi0PCM = new TLatex(0.57,0.29,Form("Signal #times %1.0f",scaleFactorPi0PCM));
	SetStyleTLatex( labelScalingPi0PCM, 0.7*textSizeLabels,4,histoPCMSignal0010->GetLineColor());
	TLatex *labelScalingPi0PCM2050 = new TLatex(0.57,0.29,Form("Signal #times %1.0f",scaleFactorPi0PCM-8));
	SetStyleTLatex( labelScalingPi0PCM2050, 0.7*textSizeLabels,4,histoPCMSignal2050->GetLineColor());

	TLatex *labelScalingEtaPCM = new TLatex(0.57,0.29,Form("Signal #times %1.0f",scaleFactorEtaPCM));
	SetStyleTLatex( labelScalingEtaPCM, 0.7*textSizeLabels,4,histoPCMSignalEta0010->GetLineColor());
	TLatex *labelScalingEtaPCM2050 = new TLatex(0.57,0.29,Form("Signal #times %1.0f",scaleFactorEtaPCM-20));
	SetStyleTLatex( labelScalingEtaPCM2050, 0.7*textSizeLabels,4,histoPCMSignalEta2050->GetLineColor());

	TLatex *labelPbPb0010 = new TLatex(0.2,0.72,"0#font[122]{-}10%");
	SetStyleTLatex( labelPbPb0010, 0.7*textSizeLabels,4,kBlack);
	TLatex *labelPbPb2050 = new TLatex(0.2,0.72,"20#font[122]{-}50%");
	SetStyleTLatex( labelPbPb2050, 0.7*textSizeLabels,4,kBlack);
	TLatex *labelRangePbPbPCM = new TLatex(0.2,0.66,Form("1.0 < #it{p}^{#gamma#gamma}_{T}< 1.2 GeV/#it{c}^{2}"));
	SetStyleTLatex( labelRangePbPbPCM, 0.7*textSizeLabels,4,kBlack);

	TLatex *labelRangePbPbPCMEta = new TLatex(0.2,0.66,Form("2.0 < #it{p}^{#gamma#gamma}_{T}< 3.0 GeV/#it{c}^{2}"));
	SetStyleTLatex( labelRangePbPbPCMEta, 0.7*textSizeLabels,4,kBlack);

	TLatex *labelRangePbPbEMCal = new TLatex(0.2,0.66,Form("12 < #it{p}^{#gamma#gamma}_{T}< 14 GeV/#it{c}^{2}"));
	SetStyleTLatex( labelRangePbPbEMCal, 0.7*textSizeLabels,4,kBlack);

	TLatex *labelPi0 = new TLatex(0.43,0.85,"#pi^{0} #rightarrow #gamma#gamma");
	SetStyleTLatex( labelPi0, 0.7*textSizeLabels,4);
	TLatex *labelEta= new TLatex(0.43,0.85,"#eta #rightarrow #gamma#gamma");
	SetStyleTLatex( labelEta, 0.7*textSizeLabels,4);

	Double_t minInvMass = 0.095;
	Double_t maxInvMass = 0.205;
	Double_t minInvMassEta = 0.400;
	Double_t maxInvMassEta = 0.750;
	Size_t markerSizeInvMass = 1;
	Width_t widthCommonFit = 1;
	
	
	
	pad6PartRatioIndMeas1->cd();
	
	SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBG0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPCMSignalPlusBG0010->GetYaxis()->SetLabelOffset(0.005);
	histoPCMSignalPlusBG0010->GetYaxis()->SetTitleFont(62);
	histoPCMSignalPlusBG0010->GetXaxis()->SetTitleFont(62);
	histoPCMSignalPlusBG0010->GetYaxis()->SetLabelFont(42);
	histoPCMSignalPlusBG0010->GetXaxis()->SetLabelFont(42);
	histoPCMSignalPlusBG0010->GetYaxis()->SetRangeUser(-10000,.95*histoPCMSignalPlusBG0010->GetMaximum());
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

	labelScalingPi0PCM->Draw();
	labelPbPb0010->Draw();
	labelMethod1->Draw();
	labelRangePbPbPCM->Draw();
	labelPi0->Draw();
	labelPbPb->Draw();

	
	
	pad6PartRatioIndMeas2->cd();

	SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBG2050, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts",  0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPCMSignalPlusBG2050->GetYaxis()->SetLabelOffset(0.005);
	histoPCMSignalPlusBG2050->GetYaxis()->SetTitleFont(62);
	histoPCMSignalPlusBG2050->GetXaxis()->SetTitleFont(62);
	histoPCMSignalPlusBG2050->GetYaxis()->SetLabelFont(42);
	histoPCMSignalPlusBG2050->GetXaxis()->SetLabelFont(42);
	histoPCMSignalPlusBG2050->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	histoPCMSignalPlusBG2050->GetYaxis()->SetRangeUser(-100,1.*histoPCMSignalPlusBG2050->GetMaximum());
	histoPCMSignalPlusBG2050->Draw("hist,e");
	DrawGammaSetMarker(histoPCMSignal2050,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoPCMSignal2050->Draw("same,pe");
	histoPCMSignal2050->Scale(scaleFactorPi0PCM-8);
	
	TH1D* histoFitPCMSignal2050 = (TH1D*)fitPCMSignal2050->GetHistogram();
	DrawGammaSetMarker(histoFitPCMSignal2050,20,0, kBlue+2 , kBlue+2);
	histoFitPCMSignal2050->SetLineWidth(widthCommonFit);
	histoFitPCMSignal2050->Scale(scaleFactorPi0PCM-8);
	histoFitPCMSignal2050->Draw("same,hist,c");

// 	DrawGammaSetMarkerTF1( fitPCMSignal2050, 1, widthCommonFit, kBlue+2);
// 	fitPCMSignal2050->Draw("same");
	labelPbPb->Draw();
	labelPbPb2050->Draw();
	labelMethod1->Draw();
	labelPi0->Draw();
	labelScalingPi0PCM2050->Draw();
	labelRangePbPbPCM->Draw();
	
	
	
	pad6PartRatioIndMeas3->cd();

	SetStyleHistoTH1ForGraphs(histoEMCalSignalPlusBG0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoEMCalSignalPlusBG0010->GetYaxis()->SetLabelOffset(0.005);
	histoEMCalSignalPlusBG0010->GetYaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBG0010->GetXaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBG0010->GetYaxis()->SetLabelFont(42);
	histoEMCalSignalPlusBG0010->GetXaxis()->SetLabelFont(42);
// 	DrawGammaSetMarker(histoEMCalSignalPlusBG0010,20,0, kBlack , kBlack);  
// 	histoEMCalSignalPlusBG0010->SetLineWidth(0.8);
	histoEMCalSignalPlusBG0010->GetYaxis()->SetRangeUser(-10,1.4*histoEMCalSignalPlusBG0010->GetMaximum());
	histoEMCalSignalPlusBG0010->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	histoEMCalSignalPlusBG0010->Draw("hist,e");
// 	histoEMCalSignal0010->SetLineWidth(0.8);
	DrawGammaSetMarker(histoEMCalSignal0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoEMCalSignal0010->Draw("same,pe");
	histoEMCalSignal0010->Scale(scaleFactorPi0EMCal);

	TH1D* histoFitEMCalSignal0010 = (TH1D*)fitpeak0010->GetHistogram();
// 	TH1D* histoFitEMCalSignal0010 = (TH1D*)fitEMCalSignal0010->GetHistogram();
	DrawGammaSetMarker(histoFitEMCalSignal0010,20,0, kBlue+2 , kBlue+2);
	histoFitEMCalSignal0010->SetLineWidth(widthCommonFit);
	histoFitEMCalSignal0010->Scale(scaleFactorPi0EMCal);
	histoFitEMCalSignal0010->Draw("same,hist,c");
	TLatex *labelScalingPi0EMCal = new TLatex(0.57,0.35,Form("Signal #times %1.0f",scaleFactorPi0EMCal));
	SetStyleTLatex( labelScalingPi0EMCal, 0.78*textSizeLabels,4,histoPCMSignal0010->GetLineColor());
	
	labelPbPb->Draw();
	labelPbPb0010->Draw();
	labelMethod2->Draw();
	labelPi0->Draw();
	labelRangePbPbEMCal->Draw();
// 	labelScalingPi0EMCal->Draw();

	
	pad6PartRatioIndMeas4->cd();

	SetStyleHistoTH1ForGraphs(histoEMCalSignalPlusBG2050, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoEMCalSignalPlusBG2050->GetYaxis()->SetLabelOffset(0.005);
	histoEMCalSignalPlusBG2050->GetYaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBG2050->GetXaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBG2050->GetYaxis()->SetLabelFont(42);
	histoEMCalSignalPlusBG2050->GetXaxis()->SetLabelFont(42);
	DrawGammaSetMarker(histoEMCalSignalPlusBG2050,20,0, kBlack , kBlack);  
	histoEMCalSignalPlusBG2050->SetLineWidth(0.8);
	histoEMCalSignalPlusBG2050->GetYaxis()->SetRangeUser(-8,1.4*histoEMCalSignalPlusBG2050->GetMaximum());
	histoEMCalSignalPlusBG2050->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	histoEMCalSignalPlusBG2050->Draw("hist,e");
	DrawGammaSetMarker(histoEMCalSignal2050,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoEMCalSignal2050->SetLineWidth(0.8);
	histoEMCalSignal2050->Draw("same,pe");
	DrawGammaSetMarkerTF1( fitpeak2050, 1, widthCommonFit, kBlue+2);
	fitpeak2050->Draw("same");
// 	DrawGammaSetMarkerTF1( fitEMCalSignal2050, 1, widthCommonFit, kBlue+2);
// 	fitEMCalSignal2050->Draw("same");

	labelPbPb->Draw();
	labelPbPb2050->Draw();
	labelMethod2->Draw();
	labelPi0->Draw();
	labelRangePbPbEMCal->Draw();

	
	pad6PartRatioIndMeas5->cd();

	SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBGEta0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPCMSignalPlusBGEta0010->GetYaxis()->SetLabelOffset(0.005);
	histoPCMSignalPlusBGEta0010->GetYaxis()->SetTitleFont(62);
	histoPCMSignalPlusBGEta0010->GetXaxis()->SetTitleFont(62);
	histoPCMSignalPlusBGEta0010->GetYaxis()->SetLabelFont(42);
	histoPCMSignalPlusBGEta0010->GetXaxis()->SetLabelFont(42);
	histoPCMSignalPlusBGEta0010->GetYaxis()->SetRangeUser(-10000,1.05*histoPCMSignalPlusBGEta0010->GetMaximum());
	histoPCMSignalPlusBGEta0010->GetXaxis()->SetRangeUser(minInvMassEta,maxInvMassEta);
	histoPCMSignalPlusBGEta0010->Draw("hist,e");
	DrawGammaSetMarker(histoPCMSignalEta0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoPCMSignalEta0010->Draw("same,pe");
	histoPCMSignalEta0010->Scale(scaleFactorEtaPCM);
	
	TH1D* histoFitPCMSignalEta0010 = (TH1D*)fitPCMSignalEta0010->GetHistogram();
	DrawGammaSetMarker(histoFitPCMSignalEta0010,20,0, kBlue+2 , kBlue+2);
	histoFitPCMSignalEta0010->SetLineWidth(widthCommonFit);
	histoFitPCMSignalEta0010->Scale(scaleFactorEtaPCM);
	histoFitPCMSignalEta0010->Draw("same,hist,c");

	labelScalingEtaPCM->Draw();
	labelPbPb0010->Draw();
	labelMethod1->Draw();
	labelEta->Draw();
	labelRangePbPbPCMEta->Draw();
	labelPbPb->Draw();


	pad6PartRatioIndMeas6->cd();
	
	SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBGEta2050, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPCMSignalPlusBGEta2050->GetYaxis()->SetLabelOffset(0.005);
	histoPCMSignalPlusBGEta2050->GetYaxis()->SetTitleFont(62);
	histoPCMSignalPlusBGEta2050->GetXaxis()->SetTitleFont(62);
	histoPCMSignalPlusBGEta2050->GetYaxis()->SetLabelFont(42);
	histoPCMSignalPlusBGEta2050->GetXaxis()->SetLabelFont(42);
	histoPCMSignalPlusBGEta2050->GetXaxis()->SetRangeUser(minInvMassEta,maxInvMassEta);
	histoPCMSignalPlusBGEta2050->GetYaxis()->SetRangeUser(-100,1.*histoPCMSignalPlusBGEta2050->GetMaximum());
	histoPCMSignalPlusBGEta2050->Draw("hist,e");
	DrawGammaSetMarker(histoPCMSignalEta2050,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoPCMSignalEta2050->Draw("same,pe");
	histoPCMSignalEta2050->Scale(scaleFactorEtaPCM-20);
	
	TH1D* histoFitPCMSignalEta2050 = (TH1D*)fitPCMSignalEta2050->GetHistogram();
	DrawGammaSetMarker(histoFitPCMSignalEta2050,20,0, kBlue+2 , kBlue+2);
	histoFitPCMSignalEta2050->SetLineWidth(widthCommonFit);
	histoFitPCMSignalEta2050->Scale(scaleFactorEtaPCM-20);
	histoFitPCMSignalEta2050->Draw("same,hist,c");

// 	DrawGammaSetMarkerTF1( fitPCMSignalEta2050, 1, widthCommonFit, kBlue+2);
// 	fitPCMSignalEta2050->Draw("same");
	labelPbPb->Draw();
	labelPbPb2050->Draw();
	labelMethod1->Draw();
	labelScalingEtaPCM2050->Draw();
	labelEta->Draw();
	labelRangePbPbPCMEta->Draw();


	pad6PartRatioIndMeas7->cd();

	SetStyleHistoTH1ForGraphs(histoEMCalSignalPlusBGEta0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoEMCalSignalPlusBGEta0010->GetYaxis()->SetLabelOffset(0.005);
	histoEMCalSignalPlusBGEta0010->GetYaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBGEta0010->GetXaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBGEta0010->GetYaxis()->SetLabelFont(42);
	histoEMCalSignalPlusBGEta0010->GetXaxis()->SetLabelFont(42);
	DrawGammaSetMarker(histoEMCalSignalPlusBGEta0010,20,0, kBlack , kBlack);  
	histoEMCalSignalPlusBGEta0010->SetLineWidth(0.8);
	histoEMCalSignalPlusBGEta0010->GetYaxis()->SetRangeUser(-1,0.8*histoEMCalSignalPlusBGEta0010->GetMaximum());
	histoEMCalSignalPlusBGEta0010->GetXaxis()->SetRangeUser(minInvMassEta,maxInvMassEta);
	histoEMCalSignalPlusBGEta0010->Draw("hist,e");
	histoEMCalSignalEta0010->SetLineWidth(0.8);
	DrawGammaSetMarker(histoEMCalSignalEta0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoEMCalSignalEta0010->Draw("same,pe");
	histoEMCalSignalEta0010->Scale(scaleFactorPi0EMCal);
	TH1D* histoFitEMCalSignalEta0010 = (TH1D*)fitpeakEta0010->GetHistogram();
// 	TH1D* histoFitEMCalSignalEta0010 = (TH1D*)fitEMCalSignalEta0010->GetHistogram();
	DrawGammaSetMarker(histoFitEMCalSignalEta0010,20,0, kBlue+2 , kBlue+2);
	histoFitEMCalSignalEta0010->SetLineWidth(widthCommonFit);
// 	histoFitEMCalSignalEta0010->Scale(scaleFactorPi0EMCal);
	histoFitEMCalSignalEta0010->Draw("same,hist,c");
// 	TLatex *labelScalingPi0EMCal = new TLatex(0.57,0.35,Form("Signal #times %1.0f",scaleFactorPi0EMCal));
// 	SetStyleTLatex( labelScalingPi0EMCal, 0.8*textSizeLabels,4,histoPCMSignalEta0010->GetLineColor());
	
	labelPbPb->Draw();
	labelPbPb0010->Draw();
	labelMethod2->Draw();
	labelEta->Draw();
	labelRangePbPbEMCal->Draw();
// 	labelScalingPi0EMCal->Draw();



	pad6PartRatioIndMeas8->cd();

	SetStyleHistoTH1ForGraphs(histoEMCalSignalPlusBGEta2050, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoEMCalSignalPlusBGEta2050->GetYaxis()->SetLabelOffset(0.005);
	histoEMCalSignalPlusBGEta2050->GetYaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBGEta2050->GetXaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBGEta2050->GetYaxis()->SetLabelFont(42);
	histoEMCalSignalPlusBGEta2050->GetXaxis()->SetLabelFont(42);
	DrawGammaSetMarker(histoEMCalSignalPlusBGEta2050,20,0, kBlack , kBlack);  
	histoEMCalSignalPlusBGEta2050->SetLineWidth(0.8);
	histoEMCalSignalPlusBGEta2050->GetYaxis()->SetRangeUser(-1,0.8*histoEMCalSignalPlusBGEta2050->GetMaximum());
	histoEMCalSignalPlusBGEta2050->GetXaxis()->SetRangeUser(minInvMassEta,maxInvMassEta);
	histoEMCalSignalPlusBGEta2050->Draw("hist,e");
	DrawGammaSetMarker(histoEMCalSignalEta2050,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoEMCalSignalEta2050->SetLineWidth(0.8);
	histoEMCalSignalEta2050->Draw("same,pe");
	DrawGammaSetMarkerTF1( fitpeakEta2050, 1, widthCommonFit, kBlue+2);
	fitpeakEta2050->Draw("same");
// 	DrawGammaSetMarkerTF1( fitEMCalSignalEta2050, 1, widthCommonFit, kBlue+2);
// 	fitEMCalSignalEta2050->Draw("same");

	labelPbPb->Draw();
	labelPbPb2050->Draw();
	labelMethod2->Draw();
	labelEta->Draw();
	labelRangePbPbEMCal->Draw();

	
	canvas6PartRatioIndMeas->SaveAs(Form("%s/InvMassPaper_%s.pdf",outputDir.Data(),dateForOutput.Data()));	

  
    
    
    
    Double_t arrayBoundsXIndMeasRatio2[3];
	Double_t arrayBoundsYIndMeasRatio2[3];
	Double_t relativeMarginsIndMeasRatioX2[3];
	Double_t relativeMarginsIndMeasRatioY2[3];
	ReturnCorrectValuesForCanvasScaling(800,800, 2, 2,0.0, 0.00, 0.00,0.0,arrayBoundsXIndMeasRatio2,arrayBoundsYIndMeasRatio2,relativeMarginsIndMeasRatioX2,relativeMarginsIndMeasRatioY2);

	TCanvas * canvas2PartRatioIndMeas = new TCanvas("canvas2PartRatioIndMeas","",0,0,800,800);  // gives the page size		
	canvas2PartRatioIndMeas->cd();

	TPad* pad2PartRatioIndMeas1 = new TPad("pad2PartRatioIndMeas1", "", arrayBoundsXIndMeasRatio2[0], arrayBoundsYIndMeasRatio2[1],arrayBoundsXIndMeasRatio2[1], arrayBoundsYIndMeasRatio2[0],-1, -1, -2);
	DrawGammaPadSettings( pad2PartRatioIndMeas1, marginLeft,marginRight,marginUp,marginDown);
	pad2PartRatioIndMeas1->Draw();
	TPad* pad2PartRatioIndMeas2 = new TPad("pad2PartRatioIndMeas2", "", arrayBoundsXIndMeasRatio2[0], arrayBoundsYIndMeasRatio2[2], arrayBoundsXIndMeasRatio2[1], arrayBoundsYIndMeasRatio2[1],-1, -1, -2);
	DrawGammaPadSettings( pad2PartRatioIndMeas2, marginLeft,marginRight,marginUp,marginDown);
	pad2PartRatioIndMeas2->Draw();
	
	TPad* pad2PartRatioIndMeas3 = new TPad("pad2PartRatioIndMeas3", "", arrayBoundsXIndMeasRatio2[1], arrayBoundsYIndMeasRatio2[1], arrayBoundsXIndMeasRatio2[2], arrayBoundsYIndMeasRatio2[0],-1, -1, -2);
	DrawGammaPadSettings( pad2PartRatioIndMeas3, marginLeft,marginRight,marginUp,marginDown);
	pad2PartRatioIndMeas3->Draw();
	TPad* pad2PartRatioIndMeas4 = new TPad("pad2PartRatioIndMeas4", "", arrayBoundsXIndMeasRatio2[1], arrayBoundsYIndMeasRatio2[2], arrayBoundsXIndMeasRatio2[2], arrayBoundsYIndMeasRatio2[1],-1, -1, -2);
	DrawGammaPadSettings( pad2PartRatioIndMeas4, marginLeft,marginRight,marginUp,marginDown);
	pad2PartRatioIndMeas4->Draw();
	
	if (pad2PartRatioIndMeas1->XtoPixel(pad2PartRatioIndMeas1->GetX2()) < pad2PartRatioIndMeas1->YtoPixel(pad2PartRatioIndMeas1->GetY1())){
	    textSizeLabels = (Double_t)textSizeLabelsPixelRatio/pad2PartRatioIndMeas1->XtoPixel(pad2PartRatioIndMeas1->GetX2()) ;
	    textSizeFac = (Double_t)1./pad2PartRatioIndMeas1->XtoPixel(pad2PartRatioIndMeas1->GetX2()) ;
	} else {
	  textSizeLabels = (Double_t)textSizeLabelsPixelRatio/pad2PartRatioIndMeas1->YtoPixel(pad2PartRatioIndMeas1->GetY1());
	  textSizeFac = (Double_t)1./pad2PartRatioIndMeas1->YtoPixel(pad2PartRatioIndMeas1->GetY1());
	}
	
	pad2PartRatioIndMeas1->cd();
	
	SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBG0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPCMSignalPlusBG0010->GetYaxis()->SetLabelOffset(0.005);
	histoPCMSignalPlusBG0010->GetYaxis()->SetTitleFont(62);
	histoPCMSignalPlusBG0010->GetXaxis()->SetTitleFont(62);
	histoPCMSignalPlusBG0010->GetYaxis()->SetLabelFont(42);
	histoPCMSignalPlusBG0010->GetXaxis()->SetLabelFont(42);
	histoPCMSignalPlusBG0010->GetYaxis()->SetRangeUser(-10000,histoPCMSignalPlusBG0010->GetMaximum());
	histoPCMSignalPlusBG0010->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	DrawGammaSetMarker(histoPCMSignalPlusBG0010,1,0, kBlack , kBlack);  
	histoPCMSignalPlusBG0010->Draw("hist,e");
// 	DrawGammaSetMarker(histoPCMSignal0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
// 	histoPCMSignal0010->Scale(scaleFactorPi0PCM);
	histoPCMSignal0010->Draw("same,pe");
	
// 	TH1D* histoFitPCMSignal0010 = (TH1D*)fitPCMSignal0010->GetHistogram();
	DrawGammaSetMarker(histoFitPCMSignal0010,20,0, kBlue+2 , kBlue+2);
	histoFitPCMSignal0010->SetLineWidth(widthCommonFit);
// 	histoFitPCMSignal0010->Scale(scaleFactorPi0PCM);
	histoFitPCMSignal0010->Draw("same,hist,c");

	labelScalingPi0PCM->Draw();
	labelPbPb0010->Draw();
	labelMethod1->Draw();
	labelRangePbPbPCM->Draw();
	labelPi0->Draw();
	labelPbPb->Draw();
	
	
	pad2PartRatioIndMeas2->cd();

	SetStyleHistoTH1ForGraphs(histoEMCalSignalPlusBG0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoEMCalSignalPlusBG0010->GetYaxis()->SetLabelOffset(0.005);
	histoEMCalSignalPlusBG0010->GetYaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBG0010->GetXaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBG0010->GetYaxis()->SetLabelFont(42);
	histoEMCalSignalPlusBG0010->GetXaxis()->SetLabelFont(42);
	histoEMCalSignalPlusBG0010->GetYaxis()->SetRangeUser(-10,histoEMCalSignalPlusBG0010->GetMaximum());
	histoEMCalSignalPlusBG0010->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
	DrawGammaSetMarker(histoEMCalSignalPlusBG0010,1,0, kBlack , kBlack);  
	histoEMCalSignalPlusBG0010->Draw("hist,e");
// 	histoEMCalSignalPlusBG0010->SetLineWidth(0.8);
// 	histoEMCalSignal0010->SetLineWidth(0.8);
	DrawGammaSetMarker(histoEMCalSignal0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoEMCalSignal0010->Draw("same,pe");
// 	histoEMCalSignal0010->Scale(scaleFactorPi0EMCal);

// 	TH1D* histoFitEMCalSignal0010 = (TH1D*)fitpeak0010->GetHistogram();
// // 	TH1D* histoFitEMCalSignal0010 = (TH1D*)fitEMCalSignal0010->GetHistogram();
	DrawGammaSetMarker(histoFitEMCalSignal0010,20,0, kBlue+2 , kBlue+2);
	histoFitEMCalSignal0010->SetLineWidth(widthCommonFit);
// 	histoFitEMCalSignal0010->Scale(scaleFactorPi0EMCal);
	histoFitEMCalSignal0010->Draw("same,hist,c");
	
	labelPbPb->Draw();
	labelPbPb0010->Draw();
	labelMethod2->Draw();
	labelPi0->Draw();
	labelRangePbPbEMCal->Draw();

	
	pad2PartRatioIndMeas3->cd();

	SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBGEta0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoPCMSignalPlusBGEta0010->GetYaxis()->SetLabelOffset(0.005);
	histoPCMSignalPlusBGEta0010->GetYaxis()->SetTitleFont(62);
	histoPCMSignalPlusBGEta0010->GetXaxis()->SetTitleFont(62);
	histoPCMSignalPlusBGEta0010->GetYaxis()->SetLabelFont(42);
	histoPCMSignalPlusBGEta0010->GetXaxis()->SetLabelFont(42);
	histoPCMSignalPlusBGEta0010->GetYaxis()->SetRangeUser(-10000,histoPCMSignalPlusBGEta0010->GetMaximum());
	histoPCMSignalPlusBGEta0010->GetXaxis()->SetRangeUser(minInvMassEta,maxInvMassEta);
	DrawGammaSetMarker(histoPCMSignalPlusBGEta0010,1,0, kBlack , kBlack);  
	histoPCMSignalPlusBGEta0010->Draw("hist,e");
// 	DrawGammaSetMarker(histoPCMSignalEta0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoPCMSignalEta0010->Draw("same,pe");
// 	histoPCMSignalEta0010->Scale(scaleFactorEtaPCM);
	
// 	TH1D* histoFitPCMSignalEta0010 = (TH1D*)fitPCMSignalEta0010->GetHistogram();
	DrawGammaSetMarker(histoFitPCMSignalEta0010,20,0, kBlue+2 , kBlue+2);
	histoFitPCMSignalEta0010->SetLineWidth(widthCommonFit);
// 	histoFitPCMSignalEta0010->Scale(scaleFactorEtaPCM);
	histoFitPCMSignalEta0010->Draw("same,hist,c");

	labelScalingEtaPCM->Draw();
	labelPbPb0010->Draw();
	labelMethod1->Draw();
	labelEta->Draw();
	labelRangePbPbPCMEta->Draw();
	labelPbPb->Draw();

	
	pad2PartRatioIndMeas4->cd();

	SetStyleHistoTH1ForGraphs(histoEMCalSignalPlusBGEta0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
	histoEMCalSignalPlusBGEta0010->GetYaxis()->SetLabelOffset(0.005);
	histoEMCalSignalPlusBGEta0010->GetYaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBGEta0010->GetXaxis()->SetTitleFont(62);
	histoEMCalSignalPlusBGEta0010->GetYaxis()->SetLabelFont(42);
	histoEMCalSignalPlusBGEta0010->GetXaxis()->SetLabelFont(42);
// 	DrawGammaSetMarker(histoEMCalSignalPlusBGEta0010,20,0, kBlack , kBlack);  
// 	histoEMCalSignalPlusBGEta0010->SetLineWidth(0.8);
// 	histoEMCalSignalPlusBGEta0010->GetYaxis()->SetRangeUser(-1,0.8*histoEMCalSignalPlusBGEta0010->GetMaximum());
// 	histoEMCalSignalPlusBGEta0010->GetXaxis()->SetRangeUser(minInvMassEta,maxInvMassEta);
	histoEMCalSignalPlusBGEta0010->Draw("hist,e");
// 	histoEMCalSignalEta0010->SetLineWidth(0.8);
	DrawGammaSetMarker(histoEMCalSignalEta0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
	histoEMCalSignalEta0010->Draw("same,pe");
// 	histoEMCalSignalEta0010->Scale(scaleFactorPi0EMCal);
// 	TH1D* histoFitEMCalSignalEta0010 = (TH1D*)fitpeakEta0010->GetHistogram();
// 	TH1D* histoFitEMCalSignalEta0010 = (TH1D*)fitEMCalSignalEta0010->GetHistogram();
	DrawGammaSetMarker(histoFitEMCalSignalEta0010,20,0, kBlue+2 , kBlue+2);
	histoFitEMCalSignalEta0010->SetLineWidth(widthCommonFit);
// 	histoFitEMCalSignalEta0010->Scale(scaleFactorPi0EMCal);
	histoFitEMCalSignalEta0010->Draw("same,hist,c");
	
	labelPbPb->Draw();
	labelPbPb0010->Draw();
	labelMethod2->Draw();
	labelEta->Draw();
	labelRangePbPbEMCal->Draw();
	
	canvas2PartRatioIndMeas->SaveAs(Form("%s/InvMassPaper_4pad_%s.pdf",outputDir.Data(),dateForOutput.Data()));	    
    
    
    
    TCanvas * canvasSingleMeasInvMass = new TCanvas("canvasSingleMeasInvMass","",0,0,1000,1000);  // gives the page size      
    DrawGammaCanvasSettings( canvasSingleMeasInvMass,  0.1, 0.04, 0.04, 0.1);
    canvasSingleMeasInvMass->cd(); //0.035,0.04, 0.035,0.04, 1.2,1.);
    
//     SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBG0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
//     histoPCMSignalPlusBG0010->GetYaxis()->SetLabelOffset(0.005);
//     histoPCMSignalPlusBG0010->GetYaxis()->SetTitleFont(62);
//     histoPCMSignalPlusBG0010->GetXaxis()->SetTitleFont(62);
//     histoPCMSignalPlusBG0010->GetYaxis()->SetLabelFont(42);
//     histoPCMSignalPlusBG0010->GetXaxis()->SetLabelFont(42);
//     histoPCMSignalPlusBG0010->GetYaxis()->SetRangeUser(-10000,histoPCMSignalPlusBG0010->GetMaximum());
//     histoPCMSignalPlusBG0010->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
//     DrawGammaSetMarker(histoPCMSignalPlusBG0010,1,0, kBlack , kBlack);  
//     histoPCMSignalPlusBG0010->Draw("hist,e");
// //  DrawGammaSetMarker(histoPCMSignal0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
// //  histoPCMSignal0010->Scale(scaleFactorPi0PCM);
//     histoPCMSignal0010->Draw("same,pe");
//     
// //  TH1D* histoFitPCMSignal0010 = (TH1D*)fitPCMSignal0010->GetHistogram();
//     DrawGammaSetMarker(histoFitPCMSignal0010,20,0, kBlue+2 , kBlue+2);
//     histoFitPCMSignal0010->SetLineWidth(widthCommonFit);
// //  histoFitPCMSignal0010->Scale(scaleFactorPi0PCM);
//     histoFitPCMSignal0010->Draw("same,hist,c");
// 
//     labelScalingPi0PCM->Draw();
//     labelPbPb0010->Draw();
//     labelMethod1->Draw();
//     labelRangePbPbPCM->Draw();
//     labelPi0->Draw();
//     labelPbPb->Draw();
    
    Double_t mBinPi0PCM = histoPCMSignalPlusBGEta0010->GetBinWidth(1);
    
    SetStyleHistoTH1ForGraphs(histoPCMSignalPlusBGEta0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})",Form("Events / %.0f MeV/#it{c}^{2}",mBinPi0PCM*1000), 0.035,0.04,0.035,0.04,1.,1.2,505,505);
//     histoPCMSignalPlusBGEta0010->GetYaxis()->SetLabelOffset(0.005);
    histoPCMSignalPlusBGEta0010->GetYaxis()->SetTitleFont(62);
    histoPCMSignalPlusBGEta0010->GetXaxis()->SetTitleFont(62);
    histoPCMSignalPlusBGEta0010->GetYaxis()->SetLabelFont(42);
    histoPCMSignalPlusBGEta0010->GetXaxis()->SetLabelFont(42);
    histoPCMSignalPlusBGEta0010->GetYaxis()->SetRangeUser(1.45*histoPCMSignalEta0010->GetMinimum(),1.2*histoPCMSignalPlusBGEta0010->GetMaximum());
    histoPCMSignalPlusBGEta0010->GetXaxis()->SetRangeUser(minInvMassEta,0.7);//maxInvMassEta);
    DrawGammaSetMarker(histoPCMSignalPlusBGEta0010,1,0, kBlack , kBlack);  
    histoPCMSignalPlusBGEta0010->Draw("hist,e");
 DrawGammaSetMarker(histoPCMSignalEta0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
    histoPCMSignalEta0010->Draw("same,pe");
//  histoPCMSignalEta0010->Scale(scaleFactorEtaPCM);
    
    histoPCMSignalPlusBGEta0010->SetLineWidth(2);
    histoPCMSignalEta0010->SetLineWidth(2);
    histoFitPCMSignalEta0010->SetLineWidth(2);

//  TH1D* histoFitPCMSignalEta0010 = (TH1D*)fitPCMSignalEta0010->GetHistogram();
    DrawGammaSetMarker(histoFitPCMSignalEta0010,20,0, kBlue+2 , kBlue+2);
    histoFitPCMSignalEta0010->SetLineWidth(widthCommonFit);
//  histoFitPCMSignalEta0010->Scale(scaleFactorEtaPCM);
    histoFitPCMSignalEta0010->Draw("same,hist,c");


    TLatex *labelMeth1 = new TLatex(0.15,0.88,Form("PCM"));
    SetStyleTLatex( labelMeth1,0.04,4,kBlack);
    labelMeth1->Draw();
    TLatex *labelE= new TLatex(0.3,0.88,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelE,0.04,4);
    labelE->Draw();
    TLatex *labelPb = new TLatex(0.15,0.83,Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV"));
    SetStyleTLatex( labelPb,0.04,4,kBlack);    
    labelPb->Draw();
    TLatex *labelPbPb010 = new TLatex(0.15,0.78,"0#font[122]{-}10%");
    SetStyleTLatex( labelPbPb010,0.04,4,kBlack);
    labelPbPb010->Draw();
    TLatex *labelRangePbPCMEta = new TLatex(0.15,0.73,Form("2 < #it{p}^{#gamma#gamma}_{T}< 3 GeV/#it{c}^{2}"));
    SetStyleTLatex( labelRangePbPCMEta,0.04,4,kBlack);
    labelRangePbPCMEta->Draw();
    TLatex *labelEtaPCM = new TLatex(0.6,0.29,Form("Signal #times %1.0f",scaleFactorEtaPCM));
    SetStyleTLatex( labelEtaPCM,0.04,4,histoPCMSignalEta0010->GetLineColor());
    labelEtaPCM->Draw();
    
    TLatex *txtALICE = new TLatex(0.6,0.88,"ALICE performance");
    SetStyleTLatex( txtALICE,0.04,4);
    TLatex *txtdate = new TLatex(0.735,0.83,"19.01.2016");
    SetStyleTLatex( txtdate,0.04,4);
    txtALICE->Draw();
    txtdate->Draw();
    
    canvasSingleMeasInvMass->SaveAs(Form("%s/InvMassPaper_PCM_%s.pdf",outputDir.Data(),dateForOutput.Data()));     
    
    
    canvasSingleMeasInvMass->cd();

//     SetStyleHistoTH1ForGraphs(histoEMCalSignalPlusBG0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts", 0.78*textSizeLabels,textSizeLabels,  0.78*textSizeLabels,textSizeLabels, 0.8,0.5/(textSizeFac*marginRatio), 505, 505);
//     histoEMCalSignalPlusBG0010->GetYaxis()->SetLabelOffset(0.005);
//     histoEMCalSignalPlusBG0010->GetYaxis()->SetTitleFont(62);
//     histoEMCalSignalPlusBG0010->GetXaxis()->SetTitleFont(62);
//     histoEMCalSignalPlusBG0010->GetYaxis()->SetLabelFont(42);
//     histoEMCalSignalPlusBG0010->GetXaxis()->SetLabelFont(42);
//     histoEMCalSignalPlusBG0010->GetYaxis()->SetRangeUser(-10,histoEMCalSignalPlusBG0010->GetMaximum());
//     histoEMCalSignalPlusBG0010->GetXaxis()->SetRangeUser(minInvMass,maxInvMass);
//     DrawGammaSetMarker(histoEMCalSignalPlusBG0010,1,0, kBlack , kBlack);  
//     histoEMCalSignalPlusBG0010->Draw("hist,e");
// //  histoEMCalSignalPlusBG0010->SetLineWidth(0.8);
// //  histoEMCalSignal0010->SetLineWidth(0.8);
//     DrawGammaSetMarker(histoEMCalSignal0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
//     histoEMCalSignal0010->Draw("same,pe");
// //  histoEMCalSignal0010->Scale(scaleFactorPi0EMCal);
// 
// //  TH1D* histoFitEMCalSignal0010 = (TH1D*)fitpeak0010->GetHistogram();
// // //   TH1D* histoFitEMCalSignal0010 = (TH1D*)fitEMCalSignal0010->GetHistogram();
//     DrawGammaSetMarker(histoFitEMCalSignal0010,20,0, kBlue+2 , kBlue+2);
//     histoFitEMCalSignal0010->SetLineWidth(widthCommonFit);
// //  histoFitEMCalSignal0010->Scale(scaleFactorPi0EMCal);
//     histoFitEMCalSignal0010->Draw("same,hist,c");
//     
//     labelPbPb->Draw();
//     labelPbPb0010->Draw();
//     labelMethod2->Draw();
//     labelPi0->Draw();
//     labelRangePbPbEMCal->Draw();

    Double_t mBinPi0EMCal = histoEMCalSignalPlusBGEta0010->GetBinWidth(1);

    SetStyleHistoTH1ForGraphs(histoEMCalSignalPlusBGEta0010, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})",Form("Events / %.0f MeV/#it{c}^{2}",mBinPi0EMCal*1000), 0.035,0.04,0.035,0.04, 1.,1.2, 505, 505);
//     histoEMCalSignalPlusBGEta0010->GetYaxis()->SetLabelOffset(0.005);
    histoEMCalSignalPlusBGEta0010->GetYaxis()->SetTitleFont(62);
    histoEMCalSignalPlusBGEta0010->GetXaxis()->SetTitleFont(62);
    histoEMCalSignalPlusBGEta0010->GetYaxis()->SetLabelFont(42);
    histoEMCalSignalPlusBGEta0010->GetXaxis()->SetLabelFont(42);
//  DrawGammaSetMarker(histoEMCalSignalPlusBGEta0010,20,0, kBlack , kBlack);  
//  histoEMCalSignalPlusBGEta0010->SetLineWidth(0.8);
    histoEMCalSignalPlusBGEta0010->GetYaxis()->SetRangeUser(1.32*histoEMCalSignalEta0010->GetMinimum(),0.7*histoEMCalSignalPlusBGEta0010->GetMaximum());
//  histoEMCalSignalPlusBGEta0010->GetXaxis()->SetRangeUser(minInvMassEta,0.7);//maxInvMassEta);
    histoEMCalSignalPlusBGEta0010->Draw("hist,e");
//  histoEMCalSignalEta0010->SetLineWidth(0.8);
    DrawGammaSetMarker(histoEMCalSignalEta0010,20,markerSizeInvMass, kRed+1 , kRed+1);  
    histoEMCalSignalEta0010->Draw("same,pe");
//  histoEMCalSignalEta0010->Scale(scaleFactorPi0EMCal);
//  TH1D* histoFitEMCalSignalEta0010 = (TH1D*)fitpeakEta0010->GetHistogram();
//  TH1D* histoFitEMCalSignalEta0010 = (TH1D*)fitEMCalSignalEta0010->GetHistogram();
    DrawGammaSetMarker(histoFitEMCalSignalEta0010,20,0, kBlue+2 , kBlue+2);
    histoFitEMCalSignalEta0010->SetLineWidth(widthCommonFit);
//  histoFitEMCalSignalEta0010->Scale(scaleFactorPi0EMCal);
    histoFitEMCalSignalEta0010->Draw("same,hist,c");
    
    histoEMCalSignalPlusBGEta0010->SetLineWidth(2);
    histoEMCalSignalEta0010->SetLineWidth(2);
    histoFitEMCalSignalEta0010->SetLineWidth(2);

    labelPb->Draw();
    labelPbPb010->Draw();
    TLatex *labelMeth2 = new TLatex(0.15,0.88,Form("EMCal"));
    SetStyleTLatex( labelMeth2,0.04,4,kBlack);
    labelMeth2->Draw();
    labelE->Draw();
    TLatex *labelRangePbEMCal = new TLatex(0.15,0.73,Form("12 < #it{p}^{#gamma#gamma}_{T}< 14 GeV/#it{c}^{2}"));
    SetStyleTLatex( labelRangePbEMCal, 0.04,4,kBlack);
    labelRangePbEMCal->Draw();

    txtALICE->Draw();
    txtdate->Draw();

    
    canvasSingleMeasInvMass->SaveAs(Form("%s/InvMassPaper_EMCal_%s.pdf",outputDir.Data(),dateForOutput.Data()));     
    
    
}

