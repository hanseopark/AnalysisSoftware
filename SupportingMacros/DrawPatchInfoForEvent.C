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

void DrawPatchInfoForEvent(TString fileNameInput = ""){	
	
	TString date = ReturnDateString();
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();

	ifstream fileInput;
	fileInput.open(fileNameInput.Data(),ios_base::in);
	cout << fileNameInput.Data() << endl;

	Double_t 	patchEnergy[200];
	Double_t 	patchADC[200];
	Bool_t 		ktriggerL1G1[200];
	Bool_t 		ktriggerL1G2[200];
	Bool_t 		ktriggerL1J1[200];
	Bool_t 		ktriggerL1J2[200];
	Bool_t 		ktriggerL0[200];
	Double_t 	patchPhiLow[200];
	Double_t 	patchPhiHigh[200];
	Double_t 	patchPhiSize[200];
	Double_t 	patchEtaLow[200];
	Double_t 	patchEtaHigh[200];
	Double_t 	patchEtaSize[200];
	
	Int_t nPoints = 0;
	while(!fileInput.eof()){
		fileInput >> patchEnergy[nPoints] >> patchADC[nPoints]>> ktriggerL1G1[nPoints] >> ktriggerL1G2[nPoints] >> ktriggerL1J1[nPoints] >> ktriggerL1J2[nPoints] >> ktriggerL0[nPoints] >> patchPhiLow[nPoints] >>patchPhiHigh[nPoints] >>patchPhiSize[nPoints] >> patchEtaLow[nPoints] >>patchEtaHigh[nPoints] >>patchEtaSize[nPoints];
		cout << nPoints << "\t"  << patchEnergy[nPoints] << "\t"  <<patchADC[nPoints] << "\t" << patchPhiLow[nPoints] << "\t"  <<patchPhiHigh[nPoints] << endl;;		
		nPoints++;
	}
	fileInput.close();
	nPoints = nPoints-1;

	for (Int_t i = 0; i < nPoints; i++){
		if (patchPhiLow[i]< 0) patchPhiLow[i] = patchPhiLow[nPoints]+ TMath::Pi()*2;
		if (patchPhiHigh[i]< 0) patchPhiHigh[i] = patchPhiHigh[nPoints]+ TMath::Pi()*2;
	}
	
	TCanvas* canvasEMCALtriggerPatches = new TCanvas("canvasEMCALtriggerPatches","",200,10,1000,1000);  // gives the page size
	DrawGammaCanvasSettings( canvasEMCALtriggerPatches, 0.08, 0.01, 0.01, 0.08);   
	canvasEMCALtriggerPatches->cd();

	
	TH2F *histo2DEMCAL;
	histo2DEMCAL = new TH2F("histo2DEMCAL", "", 400,1.362,3.178, 300,-0.728,0.728);
	histo2DEMCAL->GetXaxis()->SetTitle("#phi [rad]");
	histo2DEMCAL->GetYaxis()->SetTitle("#eta");
	histo2DEMCAL->Draw();
	
	Int_t counterL1G = 0; 
	Int_t counterL1J = 0; 
	Int_t counterL0 = 0; 
	for (Int_t i = 0; i< nPoints; i++){
		Color_t boxColor = kBlack;
		if (ktriggerL1G1[i] || ktriggerL1G2[i]){
			boxColor = kGreen+2; //+counterL1G
			counterL1G++;
		}	
		if (ktriggerL1J1[i] || ktriggerL1J2[i]){
			boxColor = kRed+2;//+counterL1J;
			counterL1J++;
		}	
		if (ktriggerL0[i] ){
			boxColor = 800 + counterL0;
			counterL0++;
		}	
		TBox* boxErrorNorm0 = CreateBoxConv(boxColor, patchPhiLow[i], patchEtaLow[i] , patchPhiHigh[i], patchEtaHigh[i]);
		if (!ktriggerL0[i] ){
			boxErrorNorm0->SetFillStyle(0);
		} else {
			boxErrorNorm0->SetFillStyle(3011+counterL0);
		}	
		if (ktriggerL0[i] ){
			boxErrorNorm0->SetLineWidth(1);
			boxErrorNorm0->SetLineColor(kBlack);
			boxErrorNorm0->SetLineStyle(1);
		}	
		boxErrorNorm0->Draw();
	}
	for (Int_t i = 0; i< nPoints; i++){
		Color_t boxColor = kBlack;
		if (ktriggerL1G1[i] || ktriggerL1G2[i]) boxColor = kGreen+2;
		if (ktriggerL1J1[i] || ktriggerL1J2[i]) boxColor = kRed+2;
		if (ktriggerL0[i] ) boxColor = 800;
		TLatex *labelEnergyAdc;
		if (ktriggerL1G1[i]) labelEnergyAdc= new TLatex(patchPhiLow[i],patchEtaHigh[i]+0.004,Form("%2.2f", patchEnergy[i]));
		else  
			labelEnergyAdc= new TLatex(patchPhiLow[i],patchEtaLow[i],Form("%2.2f", patchEnergy[i]));
		SetStyleTLatex( labelEnergyAdc, 0.005,4);
		if (ktriggerL0[i]) boxColor = kBlack;
		labelEnergyAdc->SetTextColor(boxColor);
		labelEnergyAdc->SetNDC(kFALSE);
		labelEnergyAdc->Draw();
	}
	
	cout << "L0: " << counterL0 << "\t L1GA: " << counterL1G << "\t L1JE: " << counterL1J << endl;
	
	canvasEMCALtriggerPatches->SaveAs("TriggerVisualization.eps");
	
	
	
}
	
