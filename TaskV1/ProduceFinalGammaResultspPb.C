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
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
//#include "AliHEPDataParser.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "CalculateGammaToPi0.h"


//-----------------------------------------------------------------------------
Double_t FitSpectrum(Double_t *x, Double_t *p)
{
   // Double_t hagd = p[0]*(p[2]-1)*(p[2]-2)/p[1]/p[2]/TMath::TwoPi()/TMath::Power((1.+x[0]/p[1]/p[2]),p[2]);
   Double_t hagd = p[0]/TMath::TwoPi()/TMath::Power((1.+x[0]/p[1]/p[2]),p[2]);
   Double_t expo = p[3]*TMath::Exp(-x[0]/p[4]);
   return hagd+expo;
}
//-----------------------------------------------------------------------------
Double_t Hagedorn(Double_t *x, Double_t *p)
{
   Double_t hagd = p[0]/TMath::TwoPi()/TMath::Power((1.+x[0]/p[1]/p[2]),p[2]);
   return hagd;
}
//-----------------------------------------------------------------------------
Double_t Exponent(Double_t *x, Double_t *p)
{
   Double_t expo = p[0]*TMath::Exp(-x[0]/p[1]);
   return expo;
}



void ProduceFinalGammaResultspPb(TString inputFileName = "",TString cutSel = "", TString option = "",TString suffix = "eps"){

	gSystem->Load("libCore.so");
	gSystem->Load("libTree.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libMinuit");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");

	StyleSettingsThesis();
	SetPlotStyle();
	TColor *white = gROOT->GetColor(0);
	white->SetAlpha(0.);

	TString outputDir = Form("FinalGammaResults/%s",option.Data());
	gSystem->Exec("mkdir -p "+outputDir);
  
	TString date = ReturnDateString();
	TString fileNameSysErrPi0 ="SystematicErrorsInput/SystematicErrorAveraged_Pi0_7TeV_24_Apr_2012.dat"; // default
	TString collisionSystem= ReturnFullCollisionsSystem(option);
	
	TFile *fileInput = new TFile(inputFileName);
	cout << "here" << endl;
	TH1D *DoubleRatio = (TH1D*) fileInput->Get("DoubleRatioConversionTrueEffPurity");
	cout << "here" << endl;
	TH1D *DoubleRatioPi0Fit = (TH1D*) fileInput->Get("DoubleRatioConversionFitPurity");
	cout << "here" << endl;
	TH1D *IncRatio = (TH1D*) fileInput->Get("IncRatioPurity_trueEff");
	cout << "here" << endl;
	TH1D *IncRatioPi0Fit = (TH1D*) fileInput->Get("histoIncRatioFitPurity");
	cout << "here" << endl;
	TH1D *Gamma = (TH1D*) fileInput->Get("histoGammaSpecCorrPurity");
	cout << "here" << endl;
	TH1D *Pi0 = (TH1D*) fileInput->Get("CorrectedYieldTrueEff");
	cout << "here" << endl;
	TH1D *Pi0Fit = (TH1D*) fileInput->Get("CorrectedYieldTrueEffPi0Fit");
	cout << "here" << endl;
// 	TH1D *CocktailRatioPi0 = (TH1D*) fileInput->Get("cocktailAllGammaPi0");
// 	cout << "here" << endl;
// 	CocktailRatioPi0->SetName("cocktailAllGammaPi0");
// 	cout << "here" << endl;
	
	TString fileNameSysErrDoubleRatio ="SystematicErrorsNew/SystematicErrorAveraged_DoubleRatio_pPb_5.023TeVMB_4_Apr_2014.dat";
	TString fileNameSysErrIncGamma ="SystematicErrorsNew/SystematicErrorAveraged_GammaInc_pPb_5.023TeVMB_4_Apr_2014.dat"; 
	TString fileNameSysErrIncRatio ="SystematicErrorsNew/SystematicErrorAveraged_IncRatio_pPb_5.023TeVMB_4_Apr_2014.dat"; 
	
	ifstream 		fileSysErrDoubleRatio;
	fileSysErrDoubleRatio.open(fileNameSysErrDoubleRatio,ios_base::in);
	Double_t 		relSystErrorDoubleRatioUp[50];
	Double_t 		relSystErrorDoubleRatioDown[50];
	Double_t 		relSystErrorWOMaterialDoubleRatioUp[50];
	Double_t 		relSystErrorWOMaterialDoubleRatioDown[50];
	Double_t 		relSystErrorADoubleRatioUp[50];
	Double_t 		relSystErrorADoubleRatioDown[50];
	Double_t 		relSystErrorBDoubleRatioUp[50];
	Double_t 		relSystErrorBDoubleRatioDown[50];
	Double_t 		relSystErrorCDoubleRatioUp[50];
	Double_t 		relSystErrorCDoubleRatioDown[50];
	
	cout << fileNameSysErrDoubleRatio << endl;
	Int_t nPointsErrors=0;
	while(!fileSysErrDoubleRatio.eof() && nPointsErrors < 100){
		fileSysErrDoubleRatio >> relSystErrorDoubleRatioDown[nPointsErrors] >> relSystErrorDoubleRatioUp[nPointsErrors]>>	relSystErrorWOMaterialDoubleRatioDown[nPointsErrors] >> relSystErrorWOMaterialDoubleRatioUp[nPointsErrors] >> relSystErrorADoubleRatioDown[nPointsErrors] >> relSystErrorADoubleRatioUp[nPointsErrors]>> relSystErrorBDoubleRatioDown[nPointsErrors] >> relSystErrorBDoubleRatioUp[nPointsErrors]>> relSystErrorCDoubleRatioDown[nPointsErrors] >> relSystErrorCDoubleRatioUp[nPointsErrors];
		cout << nPointsErrors << "\t"  << relSystErrorDoubleRatioDown[nPointsErrors] << "\t"  <<relSystErrorDoubleRatioUp[nPointsErrors] << "\t" << relSystErrorWOMaterialDoubleRatioDown[nPointsErrors] << "\t"  <<relSystErrorWOMaterialDoubleRatioUp[nPointsErrors] << endl;;		
		nPointsErrors++;
	}

	ifstream 		fileSysErrIncRatio;
	fileSysErrIncRatio.open(fileNameSysErrIncRatio,ios_base::in);
	Double_t 		relSystErrorIncRatioUp[50];
	Double_t 		relSystErrorIncRatioDown[50];
	Double_t 		relSystErrorWOMaterialIncRatioUp[50];
	Double_t 		relSystErrorWOMaterialIncRatioDown[50];
	Double_t 		relSystErrorAIncRatioUp[50];
	Double_t 		relSystErrorAIncRatioDown[50];
	Double_t 		relSystErrorBIncRatioUp[50];
	Double_t 		relSystErrorBIncRatioDown[50];
	Double_t 		relSystErrorCIncRatioUp[50];
	Double_t 		relSystErrorCIncRatioDown[50];
	
	cout << fileNameSysErrIncRatio << endl;
	nPointsErrors=0;
	while(!fileSysErrIncRatio.eof() && nPointsErrors < 100){
		fileSysErrIncRatio >> relSystErrorIncRatioDown[nPointsErrors] >> relSystErrorIncRatioUp[nPointsErrors]>>	relSystErrorWOMaterialIncRatioDown[nPointsErrors] >> relSystErrorWOMaterialIncRatioUp[nPointsErrors] >> relSystErrorAIncRatioDown[nPointsErrors] >> relSystErrorAIncRatioUp[nPointsErrors]>> relSystErrorBIncRatioDown[nPointsErrors] >> relSystErrorBIncRatioUp[nPointsErrors]>> relSystErrorCIncRatioDown[nPointsErrors] >> relSystErrorCIncRatioUp[nPointsErrors];
		cout << nPointsErrors << "\t"  << relSystErrorIncRatioDown[nPointsErrors] << "\t"  <<relSystErrorIncRatioUp[nPointsErrors] << "\t" << relSystErrorWOMaterialIncRatioDown[nPointsErrors] << "\t"  <<relSystErrorWOMaterialIncRatioUp[nPointsErrors] << endl;;		
		nPointsErrors++;
	}

	ifstream 		fileSysErrIncGamma;
	fileSysErrIncGamma.open(fileNameSysErrIncGamma,ios_base::in);
	Double_t 		relSystErrorIncGammaUp[50];
	Double_t 		relSystErrorIncGammaDown[50];
	Double_t 		relSystErrorWOMaterialIncGammaUp[50];
	Double_t 		relSystErrorWOMaterialIncGammaDown[50];
	Double_t 		relSystErrorAIncGammaUp[50];
	Double_t 		relSystErrorAIncGammaDown[50];
	Double_t 		relSystErrorBIncGammaUp[50];
	Double_t 		relSystErrorBIncGammaDown[50];
	Double_t 		relSystErrorCIncGammaUp[50];
	Double_t 		relSystErrorCIncGammaDown[50];
	
	cout << fileNameSysErrIncGamma << endl;
	nPointsErrors=0;
	while(!fileSysErrIncGamma.eof() && nPointsErrors < 100){
		fileSysErrIncGamma >> relSystErrorIncGammaDown[nPointsErrors] >> relSystErrorIncGammaUp[nPointsErrors]>>	relSystErrorWOMaterialIncGammaDown[nPointsErrors] >> relSystErrorWOMaterialIncGammaUp[nPointsErrors] >> relSystErrorAIncGammaDown[nPointsErrors] >> relSystErrorAIncGammaUp[nPointsErrors]>> relSystErrorBIncGammaDown[nPointsErrors] >> relSystErrorBIncGammaUp[nPointsErrors]>> relSystErrorCIncGammaDown[nPointsErrors] >> relSystErrorCIncGammaUp[nPointsErrors];
		cout << nPointsErrors << "\t"  << relSystErrorIncGammaDown[nPointsErrors] << "\t"  <<relSystErrorIncGammaUp[nPointsErrors] << "\t" << relSystErrorWOMaterialIncGammaDown[nPointsErrors] << "\t"  <<relSystErrorWOMaterialIncGammaUp[nPointsErrors] << endl;;		
		nPointsErrors++;
	}
	nPointsErrors--;
	
	TGraphAsymmErrors *NLODoubleRatio 	=(TGraphAsymmErrors*) fileInput->Get("graphNLODoubleRatio");
	TGraphErrors *NLO =(TGraphErrors*) fileInput->Get("graphNLODirGamma");
	
	TGraphAsymmErrors* graphDoubleRatioFitPi0SysErr = CalculateSysErrFromRelSysHisto( DoubleRatioPi0Fit, "Pi0SystError",relSystErrorDoubleRatioDown,relSystErrorDoubleRatioUp,0, nPointsErrors);
// 			graphCorrectedYieldPi0SysErrBinShifted = CalculateSysErrFromRelSysHisto( histoPi0CorrYieldBinShifted, "Pi0SystErrorBinShifted",relSystErrorPi0Down , relSystErrorPi0Up, 2, nPointsPi0);		
// 			graphCorrectedYieldPi0SysErrA = CalculateSysErrFromRelSysHisto( histoCorrectedYieldPi0, "Pi0SystErrorA",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);
// 			graphCorrectedYieldPi0SysErrABinShifted = CalculateSysErrFromRelSysHisto( histoPi0CorrYieldBinShifted, "Pi0SystErrorABinShifted",relSystErrorWOMaterialPi0Down , relSystErrorWOMaterialPi0Up, 2, nPointsPi0);		

	TCanvas *canvasConversionFitDoubleRatioSum =  new TCanvas("canvasConversionFitDoubleRatioSum","",0,0,1350,1000);  // gives the page size
	DrawGammaCanvasSettings( canvasConversionFitDoubleRatioSum, 0.09, 0.02, 0.02, 0.09);
	canvasConversionFitDoubleRatioSum->cd();
		
		TH2F * histo2DDoubleRatioPlotting;
		histo2DDoubleRatioPlotting = new TH2F("histo2DDoubleRatioPlotting","histo2DDoubleRatioPlotting",1000,0.0,20.,1000,0.75,1.65);
		SetStyleHistoTH2ForGraphs(histo2DDoubleRatioPlotting, "#it{p}_{T} (GeV/#it{c})","(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.04,0.04, 0.04,0.04, 1.,1.);
		histo2DDoubleRatioPlotting->GetXaxis()->SetRangeUser(0.,DoubleRatioPi0Fit->GetXaxis()->GetBinUpEdge(DoubleRatioPi0Fit->GetNbinsX()));
		histo2DDoubleRatioPlotting->GetXaxis()->SetTitleFont(62);
		histo2DDoubleRatioPlotting->GetYaxis()->SetTitleFont(62);
	//       histo2DDoubleRatioPlotting->GetYaxis()->CenterTitle(kTRUE);
		histo2DDoubleRatioPlotting->DrawCopy(); 

		DrawGammaSetMarkerTGraph(graphDoubleRatioFitPi0SysErr, 20, 0.2, kBlue+2, kBlue+2);   
		graphDoubleRatioFitPi0SysErr->SetFillColor(kGray+1);
		graphDoubleRatioFitPi0SysErr->Draw("same,2,p");

		DrawGammaSetMarker(DoubleRatioPi0Fit, 20, 2.0, kBlue+2, kBlue+2);   
		
		TF1 *One = new TF1("One","1",0,16);
		One->SetLineWidth(1.2);
		One->SetLineColor(kBlack);
		One->Draw("same");

		NLODoubleRatio->SetLineColor(kAzure);
		NLODoubleRatio->SetFillColor(kAzure);
		NLODoubleRatio->SetLineWidth(3.0);
		NLODoubleRatio->SetMarkerSize(0);
		NLODoubleRatio->Draw("lp3");

		DoubleRatioPi0Fit->DrawCopy("same");
		
		TLegend* legendDoubleConversionFit = GetAndSetLegend(0.14,0.8,2,1,"Direct Photon Signal via Conversions");
		legendDoubleConversionFit->AddEntry(DoubleRatioPi0Fit,"Measured Direct Photon Signal fitted #pi^{0}, measured #eta","p");
		legendDoubleConversionFit->AddEntry(NLODoubleRatio,"pp NLO Direct Photon  pp 5.023TeV scaled N_{coll}","l");
		legendDoubleConversionFit->Draw();

	canvasConversionFitDoubleRatioSum->Print(Form("%s/DoubleRatioFit.%s",outputDir.Data(),suffix.Data()));

	


}

