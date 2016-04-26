/****************************************************************************************************************************
****** 		provided by Gamma Conversion Group, PWG4, 													*****
******		Ana Marin, marin@physi.uni-heidelberg.de													*****
******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
*****************************************************************************************************************************/

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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

void CompareNeutralPionGGAndDalitzDataALICE(TString outputDir = "eps", TString suffix = "eps"){

	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	Width_t 	widthLinesBoxes;

	TString collisionSystemPP = "pp #sqrt{#it{s}} = 2.76 TeV";		
	Size_t markerSizeComparison = 1.5;
	
	TFile* fileNeutralPionsDalitz = new TFile("Pi0_data_GammaConvV1DalitzCorrection_900366809010333211361000000900000_9079543454102_1131213600.root");
		TH1D*	histoNeutralPionsDalitzStat = (TH1D*)fileNeutralPionsDalitz->Get("CorrectedYieldTrueEff");
	
	TFile* fCombResults= new TFile("CombinedResultsPbPb.root");
		TGraphAsymmErrors*	graphYieldCombStatPi02760GeV = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPComb2760GeV_StatErr");
		TGraphAsymmErrors*	graphYieldCombSysPi02760GeV = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPComb2760GeV_SysErr");
		TGraphAsymmErrors*	graphYieldPCMStatPi02760GeV = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPCM2760GeV_StatErr");
		TGraphAsymmErrors*	graphYieldPCMSysPi02760GeV = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPCM2760GeV_SysErr");
		TGraphAsymmErrors*	graphYieldPHOSStatPi02760GeV = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPHOS2760GeV_StatErr");
		TGraphAsymmErrors*	graphYieldPHOSSysPi02760GeV = (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPHOS2760GeV_SysErr");
	
		
	
	cout << "*************************************************************************"<< endl;	
	cout << "******************************  PP  *************************************"<< endl;
	cout << "*************************************************************************"<< endl;
	TGraphAsymmErrors* graphYieldCombStatPi02760GeVCopy = (TGraphAsymmErrors*) graphYieldCombStatPi02760GeV->Clone("graphYieldCombStatPi02760GeVCopy");
	TGraphAsymmErrors* graphYieldCombSysPi02760GeVCopy = (TGraphAsymmErrors*) graphYieldCombSysPi02760GeV->Clone("graphYieldCombSysPi02760GeVCopy");
	TGraphAsymmErrors* graphYieldPCMStatPi02760GeVCopy = (TGraphAsymmErrors*) graphYieldPCMStatPi02760GeV->Clone("graphYieldPCMStatPi02760GeVCopy");
	TGraphAsymmErrors* graphYieldPCMSysPi02760GeVCopy = (TGraphAsymmErrors*) graphYieldPCMSysPi02760GeV->Clone("graphYieldPCMSysPi02760GeVCopy");
	TGraphAsymmErrors* graphYieldPHOSStatPi02760GeVCopy = (TGraphAsymmErrors*) graphYieldPHOSStatPi02760GeV->Clone("graphYieldPHOSStatPi02760GeVCopy");
	TGraphAsymmErrors* graphYieldPHOSSysPi02760GeVCopy = (TGraphAsymmErrors*) graphYieldPHOSSysPi02760GeV->Clone("graphYieldPHOSSysPi02760GeVCopy");
	
	cout << "combined Spectrum" << endl;
	TGraphErrors* graphRatioDalitzGGComb = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldCombStatPi02760GeVCopy, graphYieldCombSysPi02760GeVCopy, histoNeutralPionsDalitzStat, histoNeutralPionsDalitzStat,  kTRUE,  kTRUE)	;
	graphRatioDalitzGGComb->Print();
	
	cout << "PCM" << endl;
	TGraphErrors* graphRatioDalitzGGPCM = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPCMStatPi02760GeVCopy, graphYieldPCMSysPi02760GeVCopy, histoNeutralPionsDalitzStat, histoNeutralPionsDalitzStat,  kTRUE,  kTRUE)	;
	graphRatioDalitzGGPCM->Print();
	
	cout << "PHOS" << endl;
	TGraphErrors* graphRatioDalitzGGPHOS = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPHOSStatPi02760GeVCopy, graphYieldPHOSSysPi02760GeVCopy,histoNeutralPionsDalitzStat, histoNeutralPionsDalitzStat,  kTRUE,  kTRUE)	;
	graphRatioDalitzGGPHOS->Print();
		
		

	TCanvas* canvasCompYieldPPComb = new TCanvas("canvasCompYieldPPComb","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvasCompYieldPPComb,  0.12, 0.02, 0.02, 0.12);

	TH2F * histo2DCompCombinedRatio2;
	histo2DCompCombinedRatio2 = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.,6.,1000,0.2,10.	);
	histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,10.);
	histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,6.);
	SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "p_{T} (GeV/c)","#pi^{0}/#pi^{#pm}", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);

// 	canvasCompYieldPPComb->SetLogx();
	histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,6.);
	histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
	histo2DCompCombinedRatio2->DrawCopy();
	
		DrawGammaSetMarkerTGraphErr(graphRatioDalitzGGComb,21,markerSizeComparison, kBlack , kBlack);
		graphRatioDalitzGGComb->Draw("E1psame");
		
		TLatex *labelRatioPi02760GeV = new TLatex(0.16,0.9,"pp #sqrt{#it{s}} = 2.76 TeV");
		SetStyleTLatex( labelRatioPi02760GeV, 0.06,4);
		labelRatioPi02760GeV->Draw();

		TLegend* legendPi0CompChargedPionsPP = new TLegend(0.18,0.15,0.9,0.21);
		legendPi0CompChargedPionsPP->SetFillColor(0);
		legendPi0CompChargedPionsPP->SetLineColor(0);
		legendPi0CompChargedPionsPP->SetNColumns(2);
		legendPi0CompChargedPionsPP->SetTextSize(0.045);
		legendPi0CompChargedPionsPP->AddEntry(graphRatioDalitzGGComb,"#pi^{0} #gamma#gamma/#pi^{0} Dalitz","p");
		legendPi0CompChargedPionsPP->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	
	DrawGammaLines(0., 20.,1., 1.,0.1,kGray,2);
	
	
	canvasCompYieldPPComb->Update();
	canvasCompYieldPPComb->Print(Form("%s/ComparisonPi0GGCombDal.%s",outputDir.Data(),suffix.Data()));

	TCanvas* canvasCompYieldPPInd = new TCanvas("canvasCompYieldPPInd","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvasCompYieldPPInd,  0.12, 0.02, 0.02, 0.12);
	
// 	canvasCompYieldPPInd->SetLogx();
	histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
	histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.6,2.1);
	histo2DCompCombinedRatio2->DrawCopy();
	
		DrawGammaSetMarkerTGraphErr(graphRatioDalitzGGPCM,21,markerSizeComparison, kBlack , kBlack);
		graphRatioDalitzGGPCM->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatioDalitzGGPHOS,21,markerSizeComparison, kRed+2 , kRed+2);
		graphRatioDalitzGGPHOS->Draw("E1psame");
		
		labelRatioPi02760GeV->Draw();

		TLegend* legendPi0CompIndChargedPionsPP = new TLegend(0.18,0.15,0.9,0.21);
		legendPi0CompIndChargedPionsPP->SetFillColor(0);
		legendPi0CompIndChargedPionsPP->SetLineColor(0);
		legendPi0CompIndChargedPionsPP->SetNColumns(2);
		legendPi0CompIndChargedPionsPP->SetTextSize(0.045);
		legendPi0CompIndChargedPionsPP->AddEntry(graphRatioDalitzGGPCM,"#pi^{0} #gamma#gamma/#pi^{0} Dalitz (PCM)","p");
		legendPi0CompIndChargedPionsPP->AddEntry(graphRatioDalitzGGPHOS,"#pi^{0} #gamma#gamma/#pi^{0} Dalitz (PHOS)","p");
		legendPi0CompIndChargedPionsPP->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	
	DrawGammaLines(0., 20.,1., 1.,0.1,kGray,2);
	
	
	canvasCompYieldPPInd->Update();
	canvasCompYieldPPInd->Print(Form("%s/ComparisonPi0GGIndDal.%s",outputDir.Data(),suffix.Data()));

	
}