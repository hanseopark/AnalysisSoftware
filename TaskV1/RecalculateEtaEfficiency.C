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
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"


TH1D* RecalcEtaEffi(TH1D* histoEffiEta, TH1D* histoEffiPi0, Double_t startValue){
   TH1D* histoEffiEtaNew = (TH1D*)histoEffiEta->Clone();
   for (Int_t i = histoEffiEta->FindBin(startValue); i < histoEffiEta->GetNbinsX() +1 ; i++ ){
      Double_t valueEffi = (1/TMath::Power(histoEffiPi0->GetBinError(i),2) *histoEffiPi0->GetBinContent(i) + 1/TMath::Power(histoEffiEta->GetBinError(i),2)*histoEffiEta->GetBinContent(i))/(1/TMath::Power(histoEffiPi0->GetBinError(i),2) + 1/TMath::Power(histoEffiEta->GetBinError(i),2));
      Double_t errorEffi = TMath::Power((1/TMath::Power(histoEffiPi0->GetBinError(i),2) + 1/TMath::Power(histoEffiEta->GetBinError(i),2)),-0.5);
      histoEffiEtaNew->SetBinContent(i,valueEffi);
      histoEffiEtaNew->SetBinError(i,errorEffi);
   }
   return histoEffiEtaNew;
}


void RecalculateEtaEfficiency(TString fileNameEta, TString fileNamePi0EtaBinning, TString fCutSelection, TString fSuffix, TString fEnergyFlag ,Double_t startValueChanging){
	
	TFile* filePi0Correction = 			new TFile(fileNamePi0EtaBinning,"READ");
	TH1D *histoEffiPtPi0 =					(TH1D*)filePi0Correction->Get("MesonEffiPt"); //not yet correct MesonEffiPt
	TH1D *histoEffiNarrowPtPi0 = 			(TH1D*)filePi0Correction->Get("MesonNarrowEffiPt");
	TH1D *histoEffiWidePtPi0 = 				(TH1D*)filePi0Correction->Get("MesonWideEffiPt");
	TH1D *histoEffiLeftPtPi0 = 				(TH1D*)filePi0Correction->Get("MesonLeftEffiPt");
	TH1D *histoEffiLeftNarrowPtPi0 = 		(TH1D*)filePi0Correction->Get("MesonLeftNarrowEffiPt");
	TH1D *histoEffiLeftWidePtPi0 = 			(TH1D*)filePi0Correction->Get("MesonLeftWideEffiPt");
	TH1D *histoTrueEffiPtPi0 =				(TH1D*)filePi0Correction->Get("TrueMesonEffiPt"); //not yet correct MesonEffiPt
	TH1D *histoTrueEffiNarrowPtPi0 = 		(TH1D*)filePi0Correction->Get("TrueMesonNarrowEffiPt");
	TH1D *histoTrueEffiWidePtPi0 = 			(TH1D*)filePi0Correction->Get("TrueMesonWideEffiPt");
		
	TFile* fileEtaCorrection = 			new TFile(fileNameEta,"READ");
	TH1D *histoEffiPtEta =					(TH1D*)fileEtaCorrection->Get("MesonEffiPt"); //not yet correct MesonEffiPt
	TH1D *histoEffiNarrowPtEta = 			(TH1D*)fileEtaCorrection->Get("MesonNarrowEffiPt");
	TH1D *histoEffiWidePtEta = 				(TH1D*)fileEtaCorrection->Get("MesonWideEffiPt");
	TH1D *histoEffiLeftPtEta = 				(TH1D*)fileEtaCorrection->Get("MesonLeftEffiPt");
	TH1D *histoEffiLeftNarrowPtEta = 		(TH1D*)fileEtaCorrection->Get("MesonLeftNarrowEffiPt");
	TH1D *histoEffiLeftWidePtEta = 			(TH1D*)fileEtaCorrection->Get("MesonLeftWideEffiPt");
	TH1D *histoTrueEffiPtEta =				(TH1D*)fileEtaCorrection->Get("TrueMesonEffiPt"); //not yet correct MesonEffiPt
	TH1D *histoTrueEffiNarrowPtEta = 		(TH1D*)fileEtaCorrection->Get("TrueMesonNarrowEffiPt");
	TH1D *histoTrueEffiWidePtEta = 			(TH1D*)fileEtaCorrection->Get("TrueMesonWideEffiPt");
	
	TH1D* histoTrueEffiPtEtaOld = (TH1D*)histoTrueEffiPtEta->Clone("histoTrueEffiEtaOld");
	histoEffiPtEta= RecalcEtaEffi(histoEffiPtEta, histoEffiPtPi0 ,startValueChanging);
	histoEffiNarrowPtEta= RecalcEtaEffi(histoEffiNarrowPtEta, histoEffiNarrowPtPi0 ,startValueChanging);
	histoEffiWidePtEta= RecalcEtaEffi(histoEffiWidePtEta, histoEffiWidePtPi0 ,startValueChanging);
	histoEffiLeftPtEta= RecalcEtaEffi(histoEffiLeftPtEta, histoEffiLeftPtPi0 ,startValueChanging);
	histoEffiLeftNarrowPtEta= RecalcEtaEffi(histoEffiLeftNarrowPtEta, histoEffiLeftNarrowPtPi0 ,startValueChanging);
	histoEffiLeftWidePtEta= RecalcEtaEffi(histoEffiLeftWidePtEta, histoEffiLeftWidePtPi0 ,startValueChanging);
	histoTrueEffiPtEta= RecalcEtaEffi(histoTrueEffiPtEta, histoTrueEffiPtPi0 ,startValueChanging);
	histoTrueEffiNarrowPtEta= RecalcEtaEffi(histoTrueEffiNarrowPtEta, histoTrueEffiNarrowPtPi0 ,startValueChanging);
	histoTrueEffiWidePtEta= RecalcEtaEffi(histoTrueEffiWidePtEta, histoTrueEffiWidePtPi0 ,startValueChanging);
	
	TString outputDir = Form("%s/%s/%s/ExtractSignal",fCutSelection.Data(),fEnergyFlag.Data(),fSuffix.Data());
	
		TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasEffSimple, 0.1, 0.02, 0.02, 0.09);
	canvasEffSimple->SetLogy(1);	
	
	DrawGammaSetMarker(histoTrueEffiPtEta, 20, 1.7, kBlack, kBlack);										 					
	DrawAutoGammaMesonHistos( histoTrueEffiPtEta, 
			"", "p_{t} (GeV/c)", "#epsilon_{#pi^{0},#eta}", 
			kTRUE, 2., 1e-5, kFALSE,
			kFALSE, 0., 0.7, 
			kFALSE, 0., 10.);

	DrawGammaSetMarker(histoTrueEffiPtEtaOld, 22, 1.7, kBlue, kBlue);										 
	histoTrueEffiPtEtaOld->DrawCopy("e1,same"); 	
	
	DrawGammaSetMarker(histoTrueEffiPtPi0, 21, 1.7, kRed, kRed);										 
	histoTrueEffiPtPi0->DrawCopy("e1,same"); 		
	
	
	TLegend* legendYield = new TLegend(0.4,0.15,0.95,0.35);
	legendYield->SetTextSize(0.033);			
	legendYield->SetFillColor(0);
	legendYield->SetLineColor(0);
	legendYield->AddEntry(histoTrueEffiPtEta,"recalculated");
	legendYield->AddEntry(histoTrueEffiPtEtaOld,"pure #eta");
	legendYield->AddEntry(histoTrueEffiPtPi0,"pure #pi^{0}");
	legendYield->Draw();
		
		
	//	DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3]);
	canvasEffSimple->Update();
	canvasEffSimple->SaveAs(Form("%s/Eta_ComparisonEffiPi0InEtaBinningToEta_%s.%s",outputDir.Data(),fCutSelection.Data(),fSuffix.Data()));
	delete canvasEffSimple;

	
	TFile* fileEtaCorrectionUpdated = 			new TFile(fileNameEta,"UPDATE");
		histoEffiPtEta->Write("MesonEffiPt",TObject::kOverwrite);
		histoEffiNarrowPtEta->Write("MesonNarrowEffiPt",TObject::kOverwrite);
		histoEffiWidePtEta->Write("MesonWideEffiPt",TObject::kOverwrite);
		histoEffiLeftPtEta->Write("MesonLeftEffiPt",TObject::kOverwrite);
		histoEffiLeftNarrowPtEta->Write("MesonLeftNarrowEffiPt",TObject::kOverwrite);
		histoEffiLeftWidePtEta->Write("MesonLeftWideEffiPt",TObject::kOverwrite);
		histoTrueEffiPtEta->Write("TrueMesonEffiPt",TObject::kOverwrite);
		histoTrueEffiNarrowPtEta->Write("TrueMesonNarrowEffiPt",TObject::kOverwrite);
		histoTrueEffiWidePtEta->Write("TrueMesonWideEffiPt",TObject::kOverwrite);
	fileEtaCorrectionUpdated->Write();
	fileEtaCorrectionUpdated->Close();
	
	
}
