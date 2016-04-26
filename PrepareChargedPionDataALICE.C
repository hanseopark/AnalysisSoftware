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

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

void PrepareChargedPionDataALICE(){	
	TDatime now;
	int iDate = now.GetDate();
	int iYear=iDate/10000;
	int iMonth=(iDate%10000)/100;
	int iDay=iDate%100;
	char* cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
		    "Jul","Aug","Sep","Oct","Nov","Dec"};
	char cStamp1[25],cStamp2[25];
	sprintf(cStamp1,"%i_%s_%i",iDay, cMonth[iMonth-1], iYear);
 
	
	// ************************** High Pt  charged pions 7 TeV ALICE ************************************************
	TFile *fPionChargedNew = TFile::Open("ExternalInput/IdentifiedCharged/charged_pion_pectrum_correctedfraction-20120314.root");
	TH1D *histoChargedPionSpecHighPtStat7TeVALICE = (TH1D*)fPionChargedNew->Get("hPionSpectrum_PP");
	histoChargedPionSpecHighPtStat7TeVALICE->Scale(0.50);
	TGraphAsymmErrors *graphChargedPionSpecHighPtSys7TeVALICE = (TGraphAsymmErrors*)fPionChargedNew->Get("hPionGraphAsymUncertainties_PP");
	graphChargedPionSpecHighPtSys7TeVALICE = ScaleGraph(graphChargedPionSpecHighPtSys7TeVALICE, 0.5);
	graphChargedPionSpecHighPtSys7TeVALICE->Print();
	TGraphAsymmErrors *graphChargedPionSpecHighPtSys7TeVALICECopy = (TGraphAsymmErrors*)fPionChargedNew->Get("hPionGraphAsymUncertainties_PP");
	Double_t* yValue = graphChargedPionSpecHighPtSys7TeVALICECopy->GetY();
	Int_t counter = 0;
	while ((abs(yValue[counter] - 0.) )< 1e-10){
		cout << "removed point" << counter << endl;
		counter++;
		graphChargedPionSpecHighPtSys7TeVALICE->RemovePoint(0);	  
	}
	
	
// 	// ************************** High Pt  charged pions 7 TeV ALICE 26Mar2012************************************************
// 	TFile *fPionLowPt = TFile::Open("ExternalInput/IdentifiedCharged/ChargedPionUpdate_26Mar2012.root");
//   
// 	TH1D* histoChargedPionPlusSpecLowPtStat7TeVALICE = (TH1D*)fPionLowPt->Get("hComb_ITSsa0_TPC0_TOF0");
// 	TH1D* histoChargedPionPlusSpecLowPtSys7TeVALICE = (TH1D*)fPionLowPt->Get("Syst_hComb_ITSsa0_TPC0_TOF0");
// 	TH1D* histoChargedPionMinusSpecLowPtStat7TeVALICE = (TH1D*)fPionLowPt->Get("hComb_ITSsa3_TPC3_TOF3");
// 	TH1D* histoChargedPionMinusSpecLowPtSys7TeVALICE = (TH1D*)fPionLowPt->Get("Syst_hComb_ITSsa3_TPC3_TOF3");
// 	
// 	TH1D*	histoChargedPionSpecLowPtStat7TeVALICE = (TH1D*)histoChargedPionMinusSpecLowPtStat7TeVALICE->Clone("histoChargedPionSpecLowPtStat7TeVALICE");
// 	histoChargedPionSpecLowPtStat7TeVALICE->Add(histoChargedPionPlusSpecLowPtStat7TeVALICE);
// 	histoChargedPionSpecLowPtStat7TeVALICE->Scale(0.5);
// 
// 	TH1D*	histoChargedPionSpecLowPtSys7TeVALICE = (TH1D*)histoChargedPionMinusSpecLowPtSys7TeVALICE->Clone("histoChargedPionSpecLowPtSys7TeVALICE");	
// 	
// 	for (Int_t i = 1; i < histoChargedPionSpecLowPtSys7TeVALICE->GetNbinsX()+1; i++){
// 		histoChargedPionSpecLowPtStat7TeVALICE->SetBinContent(i, histoChargedPionSpecLowPtStat7TeVALICE->GetBinContent(i)/histoChargedPionSpecLowPtStat7TeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
// 		histoChargedPionSpecLowPtStat7TeVALICE->SetBinError(i, i, histoChargedPionSpecLowPtStat7TeVALICE->GetBinError(i)/histoChargedPionSpecLowPtStat7TeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
// 		Double_t error = (histoChargedPionMinusSpecLowPtSys7TeVALICE->GetBinContent(i) *histoChargedPionMinusSpecLowPtStat7TeVALICE->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSys7TeVALICE->GetBinContent(i) *histoChargedPionPlusSpecLowPtStat7TeVALICE->GetBinContent(i))/2.;
// 		histoChargedPionSpecLowPtSys7TeVALICE->SetBinContent(i, histoChargedPionSpecLowPtStat7TeVALICE->GetBinContent(i));
// 		histoChargedPionSpecLowPtSys7TeVALICE->SetBinError(i,error/histoChargedPionSpecLowPtSys7TeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
// 		
// 	}
// 	

	// ************************** High Pt  charged pions 7 TeV ALICE 26Mar2012************************************************
	TFile *fPionLowPt = TFile::Open("ExternalInput/IdentifiedCharged/SpectraCombined_13_Marzo_2013.root");
  
	TH1D* histoChargedPionPlusSpecLowPtStat7TeVALICE = (TH1D*)fPionLowPt->Get("hComb_ITSsa0_TPC0_TOF0_HMPID0");
	TH1D* histoChargedPionPlusSpecLowPtSys7TeVALICE = (TH1D*)fPionLowPt->Get("Syst_hComb_ITSsa0_TPC0_TOF0_HMPID0");
	TH1D* histoChargedPionMinusSpecLowPtStat7TeVALICE = (TH1D*)fPionLowPt->Get("hComb_ITSsa3_TPC3_TOF3_HMPID3");
	TH1D* histoChargedPionMinusSpecLowPtSys7TeVALICE = (TH1D*)fPionLowPt->Get("Syst_hComb_ITSsa3_TPC3_TOF3_HMPID3");
	
	TH1D*	histoChargedPionSpecLowPtStat7TeVALICE = (TH1D*)histoChargedPionMinusSpecLowPtStat7TeVALICE->Clone("histoChargedPionSpecLowPtStat7TeVALICE");
	histoChargedPionSpecLowPtStat7TeVALICE->Add(histoChargedPionPlusSpecLowPtStat7TeVALICE);
	histoChargedPionSpecLowPtStat7TeVALICE->Scale(0.5);

	TH1D*	histoChargedPionSpecLowPtSys7TeVALICE = (TH1D*)histoChargedPionMinusSpecLowPtSys7TeVALICE->Clone("histoChargedPionSpecLowPtSys7TeVALICE");	
	
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSys7TeVALICE->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat7TeVALICE->SetBinContent(i, histoChargedPionSpecLowPtStat7TeVALICE->GetBinContent(i)/histoChargedPionSpecLowPtStat7TeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat7TeVALICE->SetBinError(i, i, histoChargedPionSpecLowPtStat7TeVALICE->GetBinError(i)/histoChargedPionSpecLowPtStat7TeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
		Double_t error = (histoChargedPionMinusSpecLowPtSys7TeVALICE->GetBinContent(i) *histoChargedPionMinusSpecLowPtStat7TeVALICE->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSys7TeVALICE->GetBinContent(i) *histoChargedPionPlusSpecLowPtStat7TeVALICE->GetBinContent(i))/2.;
		histoChargedPionSpecLowPtSys7TeVALICE->SetBinContent(i, histoChargedPionSpecLowPtStat7TeVALICE->GetBinContent(i));
		histoChargedPionSpecLowPtSys7TeVALICE->SetBinError(i,error/histoChargedPionSpecLowPtSys7TeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
		
	}
	



	//*********************************** CMS charged pion results 0.9 TeV**********************************************************
	  Double_t chargedPionCMS900GeV_xval[22] = { 0.125, 0.175, 0.225, 0.275, 0.325,
												  0.375, 0.425, 0.475, 0.525, 0.575,
												  0.625, 0.675, 0.725, 0.775, 0.825, 
												  0.875, 0.925, 0.975, 1.025, 1.075,
												  1.125, 1.175 };
	Double_t ptbins900GeV[23] = {	0.1, 0.15, 0.2, 0.25, 0.3,
								0.35, 0.4, 0.45, 0.5, 0.55,
								0.6, 0.65, 0.7, 0.75, 0.8, 
								0.85, 0.9, 0.95, 1.0, 1.05,
								1.1, 1.15, 1.2};
	Double_t chargedPionPlusCMS900GeV_xerrminus[22] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
	0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
	0.025, 0.025, 0.025 };
	Double_t chargedPionPlusCMS900GeV_yval[22] = { 3.094, 3.861, 3.843, 3.494, 3.023,
																2.584, 2.187, 1.836, 1.541, 1.293,
																1.085, 0.9125, 0.77, 0.6521, 0.5514,
																0.4708, 0.4042, 0.347, 0.3008, 0.2591,
																0.2201, 0.1933 };
	Double_t chargedPionPlusCMS900GeV_yerrminus[22] = {0.269, 0.124, 0.092, 0.076, 0.063,
																0.053, 0.044, 0.037, 0.031, 0.026,
																0.022, 0.0188, 0.0163, 0.0145, 0.0126,
																0.0111, 0.0101, 0.0100, 0.0099, 0.0095,
																0.0094, 0.0089};
	Double_t chargedPionPlusCMS900GeV_ystatminus[22] = {0.002, 0.002, 0.002, 0.002, 0.002,
																		0.001, 0.001, 0.001, 0.001, 0.001,
																		0.001, 0.001, 9.0E-4, 9.0E-4, 9.0E-4,
																		9.0E-4, 9.0E-4, 9.0E-4, 0.001, 0.001,
																		0.0012, 0.0017 };
	Double_t chargedPionMinusCMS900GeV_yval[22] = {2.952, 3.704, 3.768, 3.43, 2.989, 
																2.546, 2.157, 1.808, 1.519, 1.271,
																1.068, 0.8984, 0.7608, 0.6457, 0.5527,
																0.4681, 0.3991, 0.3386, 0.2913, 0.2646,
																0.2335, 0.2064};
	Double_t chargedPionMinusCMS900GeV_yerrminus[22] = {0.271, 0.124, 0.092, 0.073, 0.061,
																		0.052, 0.044, 0.037, 0.031, 0.026,
																		0.022, 0.0186, 0.0161, 0.0142, 0.0125,
																		0.0108, 0.0097, 0.0093, 0.0086, 0.0087,
																		0.0091, 0.0081};
	Double_t chargedPionMinusCMS900GeV_ystatminus[22] = {0.002, 0.002, 0.002, 0.002, 0.002, 
																		0.001, 0.001, 0.001, 0.001, 0.001, 
																		0.001, 0.001, 9.0E-4, 9.0E-4, 9.0E-4,
																		9.0E-4, 9.0E-4, 9.0E-4, 9.0E-4, 0.001,
																		0.0011, 0.0016 };
  
	TH1D* histoChargedPionPlusSpecLowPtStat900GeVCMS = new TH1D("histoChargedPionPlusSpecLowPtStat900GeVCMS","histoChargedPionPlusSpecLowPtStat900GeVCMS", 22, ptbins900GeV);
	TH1D* histoChargedPionPlusSpecLowPtSys900GeVCMS = new TH1D("histoChargedPionPlusSpecLowPtSys900GeVCMS","histoChargedPionPlusSpecLowPtSys900GeVCMS", 22, ptbins900GeV);
	TH1D* histoChargedPionMinusSpecLowPtStat900GeVCMS = new TH1D("histoChargedPionMinusSpecLowPtStat900GeVCMS","histoChargedPionMinusSpecLowPtStat900GeVCMS", 22, ptbins900GeV);
	TH1D* histoChargedPionMinusSpecLowPtSys900GeVCMS = new TH1D("histoChargedPionMinusSpecLowPtSys900GeVCMS","histoChargedPionMinusSpecLowPtSys900GeVCMS", 22, ptbins900GeV);
	
	for(Int_t i=1; i<23; i++){
		histoChargedPionPlusSpecLowPtStat900GeVCMS->SetBinContent(i, chargedPionPlusCMS900GeV_yval[i-1]);
		histoChargedPionPlusSpecLowPtStat900GeVCMS->SetBinError(i, chargedPionPlusCMS900GeV_ystatminus[i-1]);
		histoChargedPionPlusSpecLowPtSys900GeVCMS->SetBinContent(i, chargedPionPlusCMS900GeV_yval[i-1]);
		histoChargedPionPlusSpecLowPtSys900GeVCMS->SetBinError(i, chargedPionPlusCMS900GeV_yerrminus[i-1]);
		
		histoChargedPionMinusSpecLowPtStat900GeVCMS->SetBinContent(i, chargedPionMinusCMS900GeV_yval[i-1]);
		histoChargedPionMinusSpecLowPtStat900GeVCMS->SetBinError(i, chargedPionMinusCMS900GeV_ystatminus[i-1]);
		histoChargedPionMinusSpecLowPtSys900GeVCMS->SetBinContent(i, chargedPionMinusCMS900GeV_yval[i-1]);
		histoChargedPionMinusSpecLowPtSys900GeVCMS->SetBinError(i, chargedPionMinusCMS900GeV_yerrminus[i-1]);
	}
	TH1D*	histoChargedPionSpecLowPtSys900GeVCMS = (TH1D*)histoChargedPionMinusSpecLowPtSys900GeVCMS->Clone("histoChargedPionSpecLowPtSys900GeVCMS");
	histoChargedPionSpecLowPtSys900GeVCMS->Add(histoChargedPionPlusSpecLowPtSys900GeVCMS);
	histoChargedPionSpecLowPtSys900GeVCMS->Scale(0.5);
	for (Int_t i = 0; i < histoChargedPionSpecLowPtSys900GeVCMS->GetNbinsX(); i++){
		Double_t fractionalSystematicError =  0;
		if (histoChargedPionPlusSpecLowPtSys900GeVCMS->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSys900GeVCMS->GetBinContent(i) !=0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSys900GeVCMS->GetBinError(i)/histoChargedPionPlusSpecLowPtSys900GeVCMS->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSys900GeVCMS->GetBinError(i)/histoChargedPionMinusSpecLowPtSys900GeVCMS->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSys900GeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys900GeVCMS->GetBinContent(i)*fractionalSystematicError/100.);
		
	}   
	TH1D*	histoChargedPionSpecLowPtStat900GeVCMS = (TH1D*)histoChargedPionMinusSpecLowPtStat900GeVCMS->Clone("histoChargedPionSpecLowPtStat900GeVCMS");
	histoChargedPionSpecLowPtStat900GeVCMS->Add(histoChargedPionPlusSpecLowPtStat900GeVCMS);
	histoChargedPionSpecLowPtStat900GeVCMS->Scale(0.5);

	for (Int_t i = 1; i < histoChargedPionSpecLowPtSys900GeVCMS->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtSys900GeVCMS->SetBinContent(i, histoChargedPionSpecLowPtSys900GeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtSys900GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
		histoChargedPionSpecLowPtSys900GeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys900GeVCMS->GetBinError(i)/histoChargedPionSpecLowPtSys900GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
		histoChargedPionSpecLowPtStat900GeVCMS->SetBinContent(i, histoChargedPionSpecLowPtStat900GeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtStat900GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
		histoChargedPionSpecLowPtStat900GeVCMS->SetBinError(i, i, histoChargedPionSpecLowPtStat900GeVCMS->GetBinError(i)/histoChargedPionSpecLowPtStat900GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
	}
	

	//*********************************** ALICE charged pion results 0.9 TeV**********************************************************
	Double_t chargedPionALICE900GeV_xval[33] = { 0.11, 0.13, 0.15000000000000002, 0.16999999999999998, 0.19, 0.225, 0.275, 0.32499999999999996, 0.375, 
    0.42500000000000004, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 
    0.925, 0.975, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 
    1.9, 2.1, 2.3, 2.5 };
	Double_t ptbinsALICE900GeV[34] = {	0.1, 0.12, 0.14, 0.16, 0.18,
													0.20, 0.25, 0.30, 0.35, 0.40,
												   0.45, 0.50, 0.55, 0.60, 0.65,
													0.70, 0.75, 0.80, 0.85, 0.9,
													0.95, 1., 1.1, 1.2, 1.3, 
													1.4, 1.5, 1.6, 1.7, 1.8, 
													2., 2.2, 2.4, 2.6};

	Double_t p8038_d1x1y1_xerrminus[33] = { 	0.009999999999999995, 0.010000000000000009, 0.010000000000000009, 0.009999999999999981, 0.010000000000000009, 0.024999999999999994, 0.025000000000000022, 0.024999999999999967, 0.025000000000000022, 
    0.025000000000000022, 0.024999999999999967, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 
    0.025000000000000022, 0.025000000000000022, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 
    0.09999999999999987, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009 };
	Double_t chargedPionPlusALICE900GeV_yval[33] = { 2.658354, 2.888623, 3.064501, 2.9925, 3.115569,
															  2.977894, 2.673851, 2.362263, 2.005345, 1.663289,
															  1.412473, 1.235371, 1.00173, 0.82576, 0.697456, 
															  0.616081, 0.504638, 0.40602, 0.350737, 0.299044, 
															  0.258118, 0.206977, 0.143769, 0.109963, 0.083156,
															  0.065611, 0.047197, 0.037742, 0.030838, 0.021651, 
															  0.012858, 0.008694, 0.006367 };
	Double_t chargedPionPlusALICE900GeV_yerrminus[33] = { 0.368345, 0.09928, 0.093292, 0.069068, 0.070011, 
												0.059972, 0.055397, 0.048109, 0.040615, 0.033297,
												0.029966, 0.020852, 0.017897, 0.014441, 0.023033,
												0.020305, 0.0167, 0.013316, 0.011571, 0.00985, 
												0.008538, 0.007001, 0.004977, 0.003897, 0.00304, 
												0.002433, 0.001784, 0.001482, 0.001231,8.88E-4,
												5.57E-4, 3.85E-4, 2.93E-4};
	Double_t chargedPionPlusALICE900GeV_ystatminus[33] = { 0.034355, 0.031949, 0.03253, 0.029677, 0.030266, 0.020327, 0.01895, 0.017561, 0.015952, 
    0.014239, 0.013051, 0.01227, 0.010865, 0.009732, 0.010077, 0.009482, 0.008571, 0.007629, 0.007099, 
    0.006588, 0.006106, 0.00385, 0.003184, 0.002797, 0.002437, 0.002176, 0.001843, 0.001662, 0.001514, 
    8.97E-4, 6.8E-4, 5.69E-4, 4.9E-4 };
   
	Double_t chargedPionMinusALICE900GeV_yval[33] = { 2.664472, 2.790213, 2.976358, 3.066082, 3.026797, 2.950636, 2.668285, 2.343535, 1.998761, 
    1.716659, 1.436463, 1.221818, 1.024817, 0.81446, 0.69469, 0.595425, 0.513817, 0.410521, 0.346321, 
    0.308068, 0.256686, 0.198254, 0.143069, 0.110292, 0.083025, 0.060515, 0.046562, 0.036817, 0.028274, 
    0.020431, 0.012429, 0.008739, 0.005992 };
	Double_t chargedPionMinusALICE900GeV_yerrminus[33] = { 0.361721, 0.094854, 0.087912, 0.070138, 0.072816, 
												0.059373, 0.056005, 0.047364, 0.041263, 0.033914,
												0.029768, 0.021121, 0.018372, 0.014214, 0.023628,
												0.020205, 0.017292, 0.013744, 0.011479, 0.010148, 
												0.008487, 0.006695, 0.004908, 0.003895, 0.002985, 
												0.002261, 0.001767, 0.00144, 0.001115,8.36E-4,
												5.29E-4, 3.85E-4, 2.81E-4};
	Double_t chargedPionMinusALICE900GeV_ystatminus[33] = { 0.034388, 0.031276, 0.031371, 0.029964, 0.029549, 0.020477, 0.018907, 0.01754, 0.015626, 
    0.014267, 0.01317, 0.012069, 0.011089, 0.009505, 0.010239, 0.009452, 0.008735, 0.007768, 0.007097, 
    0.006716, 0.006117, 0.00378, 0.003204, 0.002812, 0.002439, 0.002088, 0.001829, 0.001643, 0.001451, 
    8.66E-4, 6.74E-4, 5.64E-4, 4.75E-4 };
  
	TH1D* histoChargedPionPlusSpecLowPtStat900GeVALICE = new TH1D("histoChargedPionPlusSpecLowPtStat900GeVALICE","histoChargedPionPlusSpecLowPtStat900GeVALICE", 33, ptbinsALICE900GeV);
	TH1D* histoChargedPionPlusSpecLowPtSys900GeVALICE = new TH1D("histoChargedPionPlusSpecLowPtSys900GeVALICE","histoChargedPionPlusSpecLowPtSys900GeVALICE", 33, ptbinsALICE900GeV);
	TH1D* histoChargedPionMinusSpecLowPtStat900GeVALICE = new TH1D("histoChargedPionMinusSpecLowPtStat900GeVALICE","histoChargedPionMinusSpecLowPtStat900GeVALICE", 33, ptbinsALICE900GeV);
	TH1D* histoChargedPionMinusSpecLowPtSys900GeVALICE = new TH1D("histoChargedPionMinusSpecLowPtSys900GeVALICE","histoChargedPionMinusSpecLowPtSys900GeVALICE", 33, ptbinsALICE900GeV);
	
	for(Int_t i=1; i<34; i++){
		histoChargedPionPlusSpecLowPtStat900GeVALICE->SetBinContent(i, chargedPionPlusALICE900GeV_yval[i-1]);
		histoChargedPionPlusSpecLowPtStat900GeVALICE->SetBinError(i, chargedPionPlusALICE900GeV_ystatminus[i-1]);
		histoChargedPionPlusSpecLowPtSys900GeVALICE->SetBinContent(i, chargedPionPlusALICE900GeV_yval[i-1]);
		histoChargedPionPlusSpecLowPtSys900GeVALICE->SetBinError(i, chargedPionPlusALICE900GeV_yerrminus[i-1]);
		
		histoChargedPionMinusSpecLowPtStat900GeVALICE->SetBinContent(i, chargedPionMinusALICE900GeV_yval[i-1]);
		histoChargedPionMinusSpecLowPtStat900GeVALICE->SetBinError(i, chargedPionMinusALICE900GeV_ystatminus[i-1]);
		histoChargedPionMinusSpecLowPtSys900GeVALICE->SetBinContent(i, chargedPionMinusALICE900GeV_yval[i-1]);
		histoChargedPionMinusSpecLowPtSys900GeVALICE->SetBinError(i, chargedPionMinusALICE900GeV_yerrminus[i-1]);
	}
	TH1D*	histoChargedPionSpecLowPtSys900GeVALICE = (TH1D*)histoChargedPionMinusSpecLowPtSys900GeVALICE->Clone("histoChargedPionSpecLowPtSys900GeVALICE");
	histoChargedPionSpecLowPtSys900GeVALICE->Add(histoChargedPionPlusSpecLowPtSys900GeVALICE);
	histoChargedPionSpecLowPtSys900GeVALICE->Scale(0.5);
	TH1D*	histoChargedPionSpecLowPtStat900GeVALICE = (TH1D*)histoChargedPionMinusSpecLowPtStat900GeVALICE->Clone("histoChargedPionSpecLowPtStat900GeVALICE");
	histoChargedPionSpecLowPtStat900GeVALICE->Add(histoChargedPionPlusSpecLowPtStat900GeVALICE);
	histoChargedPionSpecLowPtStat900GeVALICE->Scale(0.5);

	for (Int_t i = 1; i < histoChargedPionSpecLowPtSys900GeVALICE->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtSys900GeVALICE->SetBinContent(i, histoChargedPionSpecLowPtSys900GeVALICE->GetBinContent(i)/histoChargedPionSpecLowPtSys900GeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSys900GeVALICE->SetBinError(i, histoChargedPionSpecLowPtSys900GeVALICE->GetBinError(i)/histoChargedPionSpecLowPtSys900GeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat900GeVALICE->SetBinContent(i, histoChargedPionSpecLowPtStat900GeVALICE->GetBinContent(i)/histoChargedPionSpecLowPtStat900GeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSys900GeVALICE->GetBinContent(i) !=  0 && histoChargedPionMinusSpecLowPtSys900GeVALICE->GetBinContent(i)){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSys900GeVALICE->GetBinError(i)/histoChargedPionPlusSpecLowPtSys900GeVALICE->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSys900GeVALICE->GetBinError(i)/histoChargedPionMinusSpecLowPtSys900GeVALICE->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtStat900GeVALICE->SetBinError(i, histoChargedPionSpecLowPtStat900GeVALICE->GetBinContent(i)*fractionalSystematicError/100.);
	}
	
	
	//*********************************** CMS charged pion results 2.76 TeV**********************************************************
	Double_t chargedPionCMS_xval[22] = { 0.125, 0.175, 0.225, 0.275, 0.325,
													0.375, 0.425, 0.475, 0.525, 0.575,
													0.625, 0.675, 0.725, 0.775, 0.825,
													0.875, 0.925, 0.975, 1.025, 1.075, 
													1.125, 1.175 };
	Double_t chargedPionPlusCMS_xerrminus[22] = { 	0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
		0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
		0.025, 0.025, 0.025 };
	Double_t ptbins[23] = {	0.1, 0.15, 0.2, 0.25, 0.3,
									0.35, 0.4, 0.45, 0.5, 0.55,
									0.6, 0.65, 0.7, 0.75, 0.8, 
									0.85, 0.9, 0.95, 1.0, 1.05,
									1.1, 1.15, 1.2};
	Double_t chargedPionPlusCMS_yval[22] = { 3.828, 4.749, 4.668, 4.201, 3.66, 3.153, 2.69, 2.275, 1.928, 
		1.638, 1.389, 1.185, 1.011, 0.8673, 0.7424, 0.6375, 0.553, 0.4952, 0.4314, 
		0.3836, 0.3383, 0.2812 };
	Double_t chargedPionPlusCMS_yerrminus[22] = { 0.339, 0.147, 0.108, 0.090, 0.077, 0.064, 0.055, 0.046, 0.039, 
		0.033, 0.028, 0.025, 0.021, 0.0192, 0.0170, 0.0151, 0.0140, 0.0141, 0.0135, 
		0.0132, 0.0133, 0.0124 };
	Double_t chargedPionPlusCMS_ystatminus[22] = { 0.003, 0.003, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 
		0.001, 0.001, 0.001, 0.001, 0.0012, 0.0012, 0.0012, 0.0013, 0.0012, 0.0013, 
		0.0014, 0.0016, 0.0022 };
	
	Double_t chargedPionMinusCMS_yval[22] = { 3.596, 4.526, 4.593, 4.159, 3.612, 3.109, 2.648, 2.243, 1.903, 
		1.616, 1.375, 1.168, 0.9971, 0.8497, 0.7288, 0.6328, 0.5609, 0.4865, 0.4147, 
		0.37, 0.321, 0.2809 };
	Double_t chargedPionMinusCMS_yerrminus[22] = { 0.328, 0.145, 0.109, 0.089, 0.074, 0.063, 0.054, 0.045, 0.038, 
		0.033, 0.028, 0.024, 0.0211, 0.0187, 0.0166, 0.0147, 0.0137, 0.0134, 0.0123, 
		0.0122, 0.0120, 0.0112};
	Double_t chargedPionMinusCMS_ystatminus[22] = { 0.003, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 
		0.001, 0.001, 0.001, 0.0012, 0.0012, 0.0012, 0.0012, 0.0012, 0.0012, 0.0013, 
		0.0014, 0.0015, 0.0021 };
  
	TH1D* histoChargedPionPlusSpecLowPtStat2760GeVCMS = new TH1D("histoChargedPionPlusSpecLowPtStat2760GeVCMS","histoChargedPionPlusSpecLowPtStat2760GeVCMS", 22, ptbins);
	TH1D* histoChargedPionPlusSpecLowPtSys2760GeVCMS = new TH1D("histoChargedPionPlusSpecLowPtSys2760GeVCMS","histoChargedPionPlusSpecLowPtSys2760GeVCMS", 22, ptbins);
	TH1D* histoChargedPionMinusSpecLowPtStat2760GeVCMS = new TH1D("histoChargedPionMinusSpecLowPtStat2760GeVCMS","histoChargedPionMinusSpecLowPtStat2760GeVCMS", 22, ptbins);
	TH1D* histoChargedPionMinusSpecLowPtSys2760GeVCMS = new TH1D("histoChargedPionMinusSpecLowPtSys2760GeVCMS","histoChargedPionMinusSpecLowPtSys2760GeVCMS", 22, ptbins);
	
	for(Int_t i=1; i<23; i++){
		histoChargedPionPlusSpecLowPtStat2760GeVCMS->SetBinContent(i, chargedPionPlusCMS_yval[i-1]);
		histoChargedPionPlusSpecLowPtStat2760GeVCMS->SetBinError(i, chargedPionPlusCMS_ystatminus[i-1]);
		histoChargedPionPlusSpecLowPtSys2760GeVCMS->SetBinContent(i, chargedPionPlusCMS_yval[i-1]);
		histoChargedPionPlusSpecLowPtSys2760GeVCMS->SetBinError(i, chargedPionPlusCMS_yerrminus[i-1]);
		
		histoChargedPionMinusSpecLowPtStat2760GeVCMS->SetBinContent(i, chargedPionMinusCMS_yval[i-1]);
		histoChargedPionMinusSpecLowPtStat2760GeVCMS->SetBinError(i, chargedPionMinusCMS_ystatminus[i-1]);
		histoChargedPionMinusSpecLowPtSys2760GeVCMS->SetBinContent(i, chargedPionMinusCMS_yval[i-1]);
		histoChargedPionMinusSpecLowPtSys2760GeVCMS->SetBinError(i, chargedPionMinusCMS_yerrminus[i-1]);
	}
	TH1D*	histoChargedPionSpecLowPtSys2760GeVCMS = (TH1D*)histoChargedPionMinusSpecLowPtSys2760GeVCMS->Clone("histoChargedPionSpecLowPtSys2760GeVCMS");
	histoChargedPionSpecLowPtSys2760GeVCMS->Add(histoChargedPionPlusSpecLowPtSys2760GeVCMS);
	histoChargedPionSpecLowPtSys2760GeVCMS->Scale(0.5);
	TH1D*	histoChargedPionSpecLowPtStat2760GeVCMS = (TH1D*)histoChargedPionMinusSpecLowPtStat2760GeVCMS->Clone("histoChargedPionSpecLowPtStat2760GeVCMS");
	histoChargedPionSpecLowPtStat2760GeVCMS->Add(histoChargedPionPlusSpecLowPtStat2760GeVCMS);
	histoChargedPionSpecLowPtStat2760GeVCMS->Scale(0.5);

	for (Int_t i = 1; i < histoChargedPionSpecLowPtSys2760GeVCMS->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat2760GeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys2760GeVCMS->GetBinError(i)/histoChargedPionSpecLowPtStat2760GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
		histoChargedPionSpecLowPtStat2760GeVCMS->SetBinContent(i, histoChargedPionSpecLowPtStat2760GeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtStat2760GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
		histoChargedPionSpecLowPtSys2760GeVCMS->SetBinContent(i, histoChargedPionSpecLowPtSys2760GeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtSys2760GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSys2760GeVCMS->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtStat2760GeVCMS->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSys2760GeVCMS->GetBinError(i)/histoChargedPionPlusSpecLowPtSys2760GeVCMS->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtStat2760GeVCMS->GetBinError(i)/histoChargedPionMinusSpecLowPtStat2760GeVCMS->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSys2760GeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys2760GeVCMS->GetBinContent(i)*fractionalSystematicError/100.);
	}
	
	
	//*********************************** CMS charged pion results 7 TeV**********************************************************
	Double_t chargedPionCMS7TeV_xval[22] = { 0.125, 0.175, 0.225, 0.275, 0.325,
												  0.375, 0.425, 0.475, 0.525, 0.575,
												  0.625, 0.675, 0.725, 0.775, 0.825, 
												  0.875, 0.925, 0.975, 1.025, 1.075,
												  1.125, 1.175 };
	Double_t ptbins7TeV[23] = {	0.1, 0.15, 0.2, 0.25, 0.3,
									0.35, 0.4, 0.45, 0.5, 0.55,
									0.6, 0.65, 0.7, 0.75, 0.8, 
									0.85, 0.9, 0.95, 1.0, 1.05,
									1.1, 1.15, 1.2};
	Double_t chargedPionPlusCMS7TeV_xerrminus[22] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
		0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
		0.025, 0.025, 0.025 };
	Double_t chargedPionPlusCMS7TeV_yval[22] = {4.705, 5.801, 5.673, 5.107, 4.44,
																3.818, 3.273, 2.798, 2.392, 2.053, 
																1.753, 1.506, 1.304, 1.136, 0.9809,
																0.8585, 0.748, 0.6617, 0.5889, 0.5137,
																0.4518, 0.3356 };
	Double_t chargedPionPlusCMS7TeV_yerrminus[22] = {0.413, 0.179, 0.134, 0.111, 0.094, 
																		 0.078, 0.067, 0.057, 0.048, 0.042, 
																		 0.036, 0.031, 0.028, 0.026, 0.0229,
																		0.0206, 0.0191, 0.0200, 0.196, 0.0193, 
																		0.0206, 0.0212};
	Double_t chargedPionPlusCMS7TeV_ystatminus[22] = {0.003, 0.003, 0.003, 0.003, 0.002,
																		  0.002, 0.002, 0.002, 0.002, 0.002,
																		  0.002, 0.002, 0.002, 0.001, 0.0015,
																		  0.0015, 0.0016, 0.0016, 0.0017, 0.0019,
																		  0.0022, 0.003 };
	Double_t chargedPionMinusCMS7TeV_yval[22] = {4.44, 5.575, 5.588, 5.041, 4.397,
																	3.794, 3.246, 2.769, 2.375, 2.032,
																	1.739, 1.499, 1.29, 1.117, 0.977,
																	0.8601, 0.7628, 0.6559, 0.5884, 0.507,
																	0.4364, 0.3664 };
	Double_t chargedPionMinusCMS7TeV_yerrminus[22] = {0.406, 0.179, 0.134, 0.107, 0.090,
																		0.077, 0.066, 0.056, 0.048, 0.041,
																		0.035, 0.031, 0.028, 0.025, 0.0224,
																		0.0200, 0.0188, 0.0185, 0.0180, 0.0176, 
																		0.0183, 0.0186};
	Double_t chargedPionMinusCMS7TeV_ystatminus[22] = {0.003, 0.003, 0.003, 0.002, 0.002, 
																			0.002, 0.002, 0.002, 0.002, 0.002,
																			0.002, 0.002, 0.002, 0.001, 0.0015,
																			0.0015, 0.0016, 0.0016, 0.0017, 0.0018,
																			0.002, 0.0028 };
  
	TH1D* histoChargedPionPlusSpecLowPtStat7TeVCMS = new TH1D("histoChargedPionPlusSpecLowPtStat7TeVCMS","histoChargedPionPlusSpecLowPtStat7TeVCMS", 22, ptbins7TeV);
	TH1D* histoChargedPionPlusSpecLowPtSys7TeVCMS = new TH1D("histoChargedPionPlusSpecLowPtSys7TeVCMS","histoChargedPionPlusSpecLowPtSys7TeVCMS", 22, ptbins7TeV);
	TH1D* histoChargedPionMinusSpecLowPtStat7TeVCMS = new TH1D("histoChargedPionMinusSpecLowPtStat7TeVCMS","histoChargedPionMinusSpecLowPtStat7TeVCMS", 22, ptbins7TeV);
	TH1D* histoChargedPionMinusSpecLowPtSys7TeVCMS = new TH1D("histoChargedPionMinusSpecLowPtSys7TeVCMS","histoChargedPionMinusSpecLowPtSys7TeVCMS", 22, ptbins7TeV);
	
	for(Int_t i=1; i<23; i++){
		histoChargedPionPlusSpecLowPtStat7TeVCMS->SetBinContent(i, chargedPionPlusCMS7TeV_yval[i-1]);
		histoChargedPionPlusSpecLowPtStat7TeVCMS->SetBinError(i, chargedPionPlusCMS7TeV_ystatminus[i-1]);
		histoChargedPionPlusSpecLowPtSys7TeVCMS->SetBinContent(i, chargedPionPlusCMS7TeV_yval[i-1]);
		histoChargedPionPlusSpecLowPtSys7TeVCMS->SetBinError(i, chargedPionPlusCMS7TeV_yerrminus[i-1]);
		
		histoChargedPionMinusSpecLowPtStat7TeVCMS->SetBinContent(i, chargedPionMinusCMS7TeV_yval[i-1]);
		histoChargedPionMinusSpecLowPtStat7TeVCMS->SetBinError(i, chargedPionMinusCMS7TeV_ystatminus[i-1]);
		histoChargedPionMinusSpecLowPtSys7TeVCMS->SetBinContent(i, chargedPionMinusCMS7TeV_yval[i-1]);
		histoChargedPionMinusSpecLowPtSys7TeVCMS->SetBinError(i, chargedPionMinusCMS7TeV_yerrminus[i-1]);
	}
	TH1D*	histoChargedPionSpecLowPtSys7TeVCMS = (TH1D*)histoChargedPionMinusSpecLowPtSys7TeVCMS->Clone("histoChargedPionSpecLowPtSys7TeVCMS");
	histoChargedPionSpecLowPtSys7TeVCMS->Add(histoChargedPionPlusSpecLowPtSys7TeVCMS);
	histoChargedPionSpecLowPtSys7TeVCMS->Scale(0.5);
	TH1D*	histoChargedPionSpecLowPtStat7TeVCMS = (TH1D*)histoChargedPionMinusSpecLowPtStat7TeVCMS->Clone("histoChargedPionSpecLowPtStat7TeVCMS");
	histoChargedPionSpecLowPtStat7TeVCMS->Add(histoChargedPionPlusSpecLowPtStat7TeVCMS);
	histoChargedPionSpecLowPtStat7TeVCMS->Scale(0.5);

	for (Int_t i = 1; i < histoChargedPionSpecLowPtSys7TeVCMS->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat7TeVCMS->SetBinContent(i, histoChargedPionSpecLowPtStat7TeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtStat7TeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
		histoChargedPionSpecLowPtStat7TeVCMS->SetBinError(i, i, histoChargedPionSpecLowPtStat7TeVCMS->GetBinError(i)/histoChargedPionSpecLowPtStat7TeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
		histoChargedPionSpecLowPtSys7TeVCMS->SetBinContent(i, histoChargedPionSpecLowPtSys7TeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtSys7TeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
		histoChargedPionSpecLowPtSys7TeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys7TeVCMS->GetBinError(i)/histoChargedPionSpecLowPtSys7TeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
      
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionMinusSpecLowPtSys7TeVCMS->GetBinContent(i) != 0 && histoChargedPionPlusSpecLowPtSys7TeVCMS->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionMinusSpecLowPtSys7TeVCMS->GetBinError(i)/histoChargedPionMinusSpecLowPtSys7TeVCMS->GetBinContent(i)*100 + histoChargedPionPlusSpecLowPtSys7TeVCMS->GetBinError(i)/histoChargedPionPlusSpecLowPtSys7TeVCMS->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSys7TeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys7TeVCMS->GetBinContent(i)*fractionalSystematicError/100.);
	}
	
	//***********************************high pt spectra ********************************************************************
	TFile* fileChargedPionSpectraPublishedPP = new TFile("ExternalInput/IdentifiedCharged/pp276.fullpT.INEL.20140504.root");
	TH1D*	histoChargedPionSpecPubStatPP = (TH1D*)fileChargedPionSpectraPublishedPP->Get("hstat_pp276_pion_sum");
	TH1D*	histoChargedPionSpecPubSystPP = (TH1D*)fileChargedPionSpectraPublishedPP->Get("hsys_pp276_pion_sum");
	histoChargedPionSpecPubStatPP->Scale(0.5);
	histoChargedPionSpecPubSystPP->Scale(0.5);
	TH1D*	histoChargedKaonSpecPubStatPP = (TH1D*)fileChargedPionSpectraPublishedPP->Get("hstat_pp276_kaon_sum");
	TH1D*	histoChargedKaonSpecPubSystPP = (TH1D*)fileChargedPionSpectraPublishedPP->Get("hsys_pp276_kaon_sum");
	histoChargedKaonSpecPubStatPP->Scale(0.5);
	histoChargedKaonSpecPubSystPP->Scale(0.5);
	TH1D*	histoProtonSpecPubStatPP = (TH1D*)fileChargedPionSpectraPublishedPP->Get("hstat_pp276_proton_sum");
	TH1D*	histoProtonSpecPubSystPP = (TH1D*)fileChargedPionSpectraPublishedPP->Get("hsys_pp276_proton_sum");
	histoProtonSpecPubStatPP->Scale(0.5);
	histoProtonSpecPubSystPP->Scale(0.5);
	
	
	//***********************************high pt spectra ********************************************************************
	TFile* fileChargedPionSpectraHighPtFinal = new TFile("ExternalInputPbPb/IdentifiedCharged/SpectraHighPtPionFinal_20131101.root");
	
	TH1D*	histoChargedPionSpecHighPtStatPP = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_pp2760");
	TH1D*	histoChargedPionSpecHighPtSystPP = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_pp2760");
	histoChargedPionSpecHighPtStatPP->Scale(0.5);
	histoChargedPionSpecHighPtSystPP->Scale(0.5);
	
	TH1D*	histoChargedPionSpecHighPtStat0005 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_0_5");
	TH1D*	histoChargedPionSpecHighPtSyst0005 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_0_5");
	histoChargedPionSpecHighPtStat0005->Scale(0.5);
	histoChargedPionSpecHighPtSyst0005->Scale(0.5);

	TH1D*	histoChargedPionSpecHighPtStat0510 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_5_10");
	TH1D*	histoChargedPionSpecHighPtSyst0510 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_5_10");
	histoChargedPionSpecHighPtStat0510->Scale(0.5);
	histoChargedPionSpecHighPtSyst0510->Scale(0.5);

	TH1D*	histoChargedPionSpecHighPtStat0010 = (TH1D*)histoChargedPionSpecHighPtStat0510->Clone("histoChargedPionSpecHighPtStat0010");
	TH1D*	histoChargedPionSpecHighPtSyst0010 = (TH1D*)histoChargedPionSpecHighPtSyst0510->Clone("histoChargedPionSpecHighPtSyst0010");
	histoChargedPionSpecHighPtStat0010->Add(histoChargedPionSpecHighPtStat0005);
	histoChargedPionSpecHighPtSyst0010->Add(histoChargedPionSpecHighPtSyst0005);
	histoChargedPionSpecHighPtStat0010->Scale(0.5);
	histoChargedPionSpecHighPtSyst0010->Scale(0.5);

	for (Int_t i = 1; i < histoChargedPionSpecHighPtSyst0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoChargedPionSpecHighPtSyst0005->GetBinContent(i) != 0){
			relErrLowerCent= histoChargedPionSpecHighPtSyst0005->GetBinError(i)/histoChargedPionSpecHighPtSyst0005->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoChargedPionSpecHighPtSyst0510->GetBinContent(i) != 0){
			relErrHigherCent = histoChargedPionSpecHighPtSyst0510->GetBinError(i)/histoChargedPionSpecHighPtSyst0510->GetBinContent(i)*100 ;
		}
		
		if (relErrHigherCent > relErrLowerCent){
			histoChargedPionSpecHighPtSyst0010->SetBinError(i, histoChargedPionSpecHighPtSyst0010->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoChargedPionSpecHighPtSyst0010->SetBinError(i, histoChargedPionSpecHighPtSyst0010->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   

   
	TH1D*	histoChargedPionSpecHighPtStat1020 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_10_20");
	TH1D*	histoChargedPionSpecHighPtSyst1020 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_10_20");
	histoChargedPionSpecHighPtStat1020->Scale(0.5);
	histoChargedPionSpecHighPtSyst1020->Scale(0.5);
	
	TH1D*	histoChargedPionSpecHighPtStat2040 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_20_40");
	TH1D*	histoChargedPionSpecHighPtSyst2040 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_20_40");
	histoChargedPionSpecHighPtStat2040->Scale(0.5);
	histoChargedPionSpecHighPtSyst2040->Scale(0.5);
	
	TH1D*	histoChargedPionSpecHighPtStat4060 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_40_60");
	TH1D*	histoChargedPionSpecHighPtSyst4060 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_40_60");
	histoChargedPionSpecHighPtStat4060->Scale(0.5);
	histoChargedPionSpecHighPtSyst4060->Scale(0.5);
	
	TH1D*	histoChargedPionSpecHighPtStat6080 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_60_80");
	TH1D*	histoChargedPionSpecHighPtSyst6080 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_60_80");
	histoChargedPionSpecHighPtStat6080->Scale(0.5);
	histoChargedPionSpecHighPtSyst6080->Scale(0.5);
	
      //***********************************high pt spectra ********************************************************************
	TFile* fileChargedKaonSpectraHighPtFinal = new TFile("ExternalInputPbPb/IdentifiedCharged/SpectraHighPtKaonFinal_20131108.root");
	
	TH1D* histoChargedKaonSpecHighPtStatPP = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_pp2760");
	TH1D* histoChargedKaonSpecHighPtSystPP = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_pp2760");
	histoChargedKaonSpecHighPtStatPP->Scale(0.5);
	histoChargedKaonSpecHighPtSystPP->Scale(0.5);
	
	TH1D* histoChargedKaonSpecHighPtStat0005 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_0_5");
	TH1D* histoChargedKaonSpecHighPtSyst0005 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_0_5");
	histoChargedKaonSpecHighPtStat0005->Scale(0.5);
	histoChargedKaonSpecHighPtSyst0005->Scale(0.5);

	TH1D* histoChargedKaonSpecHighPtStat0510 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_5_10");
	TH1D* histoChargedKaonSpecHighPtSyst0510 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_5_10");
	histoChargedKaonSpecHighPtStat0510->Scale(0.5);
	histoChargedKaonSpecHighPtSyst0510->Scale(0.5);

	TH1D* histoChargedKaonSpecHighPtStat0010 = (TH1D*)histoChargedKaonSpecHighPtStat0510->Clone("histoChargedKaonSpecHighPtStat0010");
	TH1D* histoChargedKaonSpecHighPtSyst0010 = (TH1D*)histoChargedKaonSpecHighPtSyst0510->Clone("histoChargedKaonSpecHighPtSyst0010");
	histoChargedKaonSpecHighPtStat0010->Add(histoChargedKaonSpecHighPtStat0005);
	histoChargedKaonSpecHighPtSyst0010->Add(histoChargedKaonSpecHighPtSyst0005);
	histoChargedKaonSpecHighPtStat0010->Scale(0.5);
	histoChargedKaonSpecHighPtSyst0010->Scale(0.5);

	for (Int_t i = 1; i < histoChargedKaonSpecHighPtSyst0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoChargedKaonSpecHighPtSyst0005->GetBinContent(i) != 0){
			relErrLowerCent= histoChargedKaonSpecHighPtSyst0005->GetBinError(i)/histoChargedKaonSpecHighPtSyst0005->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoChargedKaonSpecHighPtSyst0510->GetBinContent(i) != 0){
			relErrHigherCent = histoChargedKaonSpecHighPtSyst0510->GetBinError(i)/histoChargedKaonSpecHighPtSyst0510->GetBinContent(i)*100 ;
		}
		
		if (relErrHigherCent > relErrLowerCent){
			histoChargedKaonSpecHighPtSyst0010->SetBinError(i, histoChargedKaonSpecHighPtSyst0010->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoChargedKaonSpecHighPtSyst0010->SetBinError(i, histoChargedKaonSpecHighPtSyst0010->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   

	TH1D* histoChargedKaonSpecHighPtStat1020 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_10_20");
	TH1D* histoChargedKaonSpecHighPtSyst1020 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_10_20");
	histoChargedKaonSpecHighPtStat1020->Scale(0.5);
	histoChargedKaonSpecHighPtSyst1020->Scale(0.5);
	
	TH1D* histoChargedKaonSpecHighPtStat2040 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_20_40");
	TH1D* histoChargedKaonSpecHighPtSyst2040 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_20_40");
	histoChargedKaonSpecHighPtStat2040->Scale(0.5);
	histoChargedKaonSpecHighPtSyst2040->Scale(0.5);
	
	TH1D* histoChargedKaonSpecHighPtStat4060 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_40_60");
	TH1D* histoChargedKaonSpecHighPtSyst4060 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_40_60");
	histoChargedKaonSpecHighPtStat4060->Scale(0.5);
	histoChargedKaonSpecHighPtSyst4060->Scale(0.5);
	
	TH1D* histoChargedKaonSpecHighPtStat6080 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_60_80");
	TH1D* histoChargedKaonSpecHighPtSyst6080 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_60_80");
	histoChargedKaonSpecHighPtStat6080->Scale(0.5);
	histoChargedKaonSpecHighPtSyst6080->Scale(0.5);
   
	// ****************************************** Low Pt pion data ********************************************************************
// 	TFile* fileChargedPionSpectraLowPtPP2760GeVInternal = new TFile("ExternalInput/IdentifiedCharged/spectracombined276TeV_4022013.root");
//    TList* listChargedPionSpectraLowPtPP2760GeVInternal = (TList*)fileChargedPionSpectraLowPtPP2760GeVInternal->Get("output");   
//    TH1D* histoChargedPionPlusSpecLowPtStatPP2760GeV = (TH1D*)listChargedPionSpectraLowPtPP2760GeVInternal->FindObject("CombinedSpectraFinalPionPlus");
//    TH1D* histoChargedPionMinusSpecLowPtStatPP2760GeV = (TH1D*)listChargedPionSpectraLowPtPP2760GeVInternal->FindObject("CombinedSpectraFinalPionMinus");
//    TH1D* histoChargedPionSpecLowPtStatPP2760GeV = (TH1D*)histoChargedPionMinusSpecLowPtStatPP2760GeV->Clone("histoChargedPionSpecLowPtStatPP2760GeV");
//    histoChargedPionSpecLowPtStatPP2760GeV->Add(histoChargedPionPlusSpecLowPtStatPP2760GeV);
//    histoChargedPionSpecLowPtStatPP2760GeV->Scale(0.5);
//    TH1D* histoChargedPionSpecLowPtSysPP2760GeV = (TH1D*)histoChargedPionSpecLowPtStatPP2760GeV->Clone("histoChargedPionSpecLowPtSysPP2760GeV");
//    
   
	TFile* fileChargedPionSpectraLowPtPP2760GeVInternal = new TFile("ExternalInput/IdentifiedCharged/276SpectraCombined_INELstat4072013.root");
	TH1D* histoChargedPionPlusSpecLowPtStatPP2760GeV = (TH1D*)fileChargedPionSpectraLowPtPP2760GeVInternal->Get("pion_plus3");
	TH1D* histoChargedPionMinusSpecLowPtStatPP2760GeV = (TH1D*)fileChargedPionSpectraLowPtPP2760GeVInternal->Get("pion_minus3");
	TH1D* histoChargedPionSpecLowPtStatPP2760GeV = (TH1D*)histoChargedPionMinusSpecLowPtStatPP2760GeV->Clone("histoChargedPionSpecLowPtStatPP2760GeV");
	histoChargedPionSpecLowPtStatPP2760GeV->Add(histoChargedPionPlusSpecLowPtStatPP2760GeV);
	histoChargedPionSpecLowPtStatPP2760GeV->Scale(0.5);
	
	TFile* fileChargedPionSpectraLowPtPP2760GeVInternalSys = new TFile("ExternalInput/IdentifiedCharged/276SpectraCombined_INEL4072013.root");
	TH1D* histoChargedPionPlusSpecLowPtSysPP2760GeV = (TH1D*)fileChargedPionSpectraLowPtPP2760GeVInternalSys->Get("pion_plus3");
	TH1D* histoChargedPionMinusSpecLowPtSysPP2760GeV = (TH1D*)fileChargedPionSpectraLowPtPP2760GeVInternalSys->Get("pion_minus3");
	TH1D* histoChargedPionSpecLowPtSysPP2760GeV = (TH1D*)histoChargedPionPlusSpecLowPtSysPP2760GeV->Clone("histoChargedPionSpecLowPtSysPP2760GeV");
	histoChargedPionSpecLowPtSysPP2760GeV->Add(histoChargedPionMinusSpecLowPtSysPP2760GeV);
	histoChargedPionSpecLowPtSysPP2760GeV->Scale(0.5);
	
	for (Int_t i = 1; i < histoChargedPionSpecLowPtStatPP2760GeV->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStatPP2760GeV->SetBinContent(i, histoChargedPionSpecLowPtStatPP2760GeV->GetBinContent(i)/histoChargedPionSpecLowPtStatPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStatPP2760GeV->SetBinError(i, histoChargedPionSpecLowPtStatPP2760GeV->GetBinError(i)/histoChargedPionSpecLowPtStatPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSysPP2760GeV->SetBinContent(i, histoChargedPionSpecLowPtSysPP2760GeV->GetBinContent(i)/histoChargedPionSpecLowPtSysPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if ( histoChargedPionMinusSpecLowPtSysPP2760GeV->GetBinContent(i) != 0 && histoChargedPionPlusSpecLowPtSysPP2760GeV->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionMinusSpecLowPtSysPP2760GeV->GetBinError(i)/histoChargedPionMinusSpecLowPtSysPP2760GeV->GetBinContent(i)*100 + histoChargedPionPlusSpecLowPtSysPP2760GeV->GetBinError(i)/histoChargedPionPlusSpecLowPtSysPP2760GeV->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSysPP2760GeV->SetBinError(i, histoChargedPionSpecLowPtSysPP2760GeV->GetBinContent(i)*fractionalSystematicError/100);
	}
	
	TFile* fileChargedIndentifiedSpectraLowPtPrelim2012 = new TFile("ExternalInputPbPb/IdentifiedCharged/SPECTRA_COMB_20120809_default.root");
	TH1D*	histoChargedPionMinusSpecLowPtStat0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent0_pion_minus");
	histoChargedPionMinusSpecLowPtStat0005->Sumw2();
	TH1D*	histoChargedPionMinusSpecLowPtSyst0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent0_pion_minus");
	histoChargedPionMinusSpecLowPtSyst0005->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtStat0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent0_pion_plus");
	histoChargedPionPlusSpecLowPtStat0005->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtSyst0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent0_pion_plus");
	histoChargedPionPlusSpecLowPtSyst0005->Sumw2();
	TH1D*	histoChargedPionSpecLowPtStat0005 = (TH1D*)histoChargedPionMinusSpecLowPtStat0005->Clone("histoChargedPionSpecLowPtStat0005");
	TH1D*	histoChargedPionSpecLowPtSyst0005 = (TH1D*)histoChargedPionMinusSpecLowPtSyst0005->Clone("histoChargedPionSpecLowPtSyst0005");
	histoChargedPionSpecLowPtStat0005->Add(histoChargedPionPlusSpecLowPtStat0005);
	histoChargedPionSpecLowPtSyst0005->Add(histoChargedPionPlusSpecLowPtSyst0005);
	histoChargedPionSpecLowPtStat0005->Scale(0.5);
	histoChargedPionSpecLowPtSyst0005->Scale(0.5);
     
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst0005->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat0005->SetBinContent(i, histoChargedPionSpecLowPtStat0005->GetBinContent(i)/histoChargedPionSpecLowPtStat0005->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat0005->SetBinError(i, histoChargedPionSpecLowPtStat0005->GetBinError(i)/histoChargedPionSpecLowPtStat0005->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSyst0005->SetBinContent(i, histoChargedPionSpecLowPtSyst0005->GetBinContent(i)/(histoChargedPionSpecLowPtSyst0005->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSyst0005->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst0005->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst0005->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst0005->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst0005->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst0005->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSyst0005->SetBinError(i, histoChargedPionSpecLowPtSyst0005->GetBinContent(i)*fractionalSystematicError/100.);
	}	

	TH1D*	histoChargedPionMinusSpecLowPtStat0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent1_pion_minus");
	histoChargedPionMinusSpecLowPtStat0510->Sumw2();
	TH1D*	histoChargedPionMinusSpecLowPtSyst0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent1_pion_minus");
	histoChargedPionMinusSpecLowPtSyst0510->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtStat0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent1_pion_plus");
	histoChargedPionPlusSpecLowPtStat0510->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtSyst0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent1_pion_plus");
	histoChargedPionPlusSpecLowPtSyst0510->Sumw2();
	TH1D*	histoChargedPionSpecLowPtStat0510 = (TH1D*)histoChargedPionMinusSpecLowPtStat0510->Clone("histoChargedPionSpecLowPtStat0510");
	TH1D*	histoChargedPionSpecLowPtSyst0510 = (TH1D*)histoChargedPionMinusSpecLowPtSyst0510->Clone("histoChargedPionSpecLowPtSyst0510");
	histoChargedPionSpecLowPtStat0510->Add(histoChargedPionPlusSpecLowPtStat0510);
	histoChargedPionSpecLowPtSyst0510->Add(histoChargedPionPlusSpecLowPtSyst0510);
	histoChargedPionSpecLowPtStat0510->Scale(0.5);
	histoChargedPionSpecLowPtSyst0510->Scale(0.5);
   
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst0510->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat0510->SetBinContent(i, histoChargedPionSpecLowPtStat0510->GetBinContent(i)/histoChargedPionSpecLowPtStat0510->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat0510->SetBinError(i, histoChargedPionSpecLowPtStat0510->GetBinError(i)/histoChargedPionSpecLowPtStat0510->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSyst0510->SetBinContent(i, histoChargedPionSpecLowPtSyst0510->GetBinContent(i)/(histoChargedPionSpecLowPtSyst0510->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSyst0510->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst0510->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst0510->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst0510->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst0510->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst0510->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSyst0510->SetBinError(i, histoChargedPionSpecLowPtSyst0510->GetBinContent(i)*fractionalSystematicError/100.);
	}
	
	TH1D*	histoChargedPionSpecLowPtStat0010 = (TH1D*)histoChargedPionSpecLowPtStat0510->Clone("histoChargedPionSpecLowPtStat0010");
	TH1D*	histoChargedPionSpecLowPtSyst0010 = (TH1D*)histoChargedPionSpecLowPtSyst0510->Clone("histoChargedPionSpecLowPtSyst0010");
	histoChargedPionSpecLowPtStat0010->Add(histoChargedPionSpecLowPtStat0005);
	histoChargedPionSpecLowPtSyst0010->Add(histoChargedPionSpecLowPtSyst0005);
	histoChargedPionSpecLowPtStat0010->Scale(0.5);
	histoChargedPionSpecLowPtSyst0010->Scale(0.5);

	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0. ;
		if (histoChargedPionSpecLowPtSyst0005->GetBinContent(i) != 0){
			relErrLowerCent = histoChargedPionSpecLowPtSyst0005->GetBinError(i)/histoChargedPionSpecLowPtSyst0005->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0. ;
		if (histoChargedPionSpecLowPtSyst0510->GetBinContent(i) != 0){
			relErrHigherCent = histoChargedPionSpecLowPtSyst0510->GetBinError(i)/histoChargedPionSpecLowPtSyst0510->GetBinContent(i)*100 ;
		}
		if (relErrHigherCent > relErrLowerCent){
			histoChargedPionSpecLowPtSyst0010->SetBinError(i, histoChargedPionSpecLowPtSyst0010->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoChargedPionSpecLowPtSyst0010->SetBinError(i, histoChargedPionSpecLowPtSyst0010->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   

	TH1D*	histoChargedPionMinusSpecLowPtStat1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent2_pion_minus");
	histoChargedPionMinusSpecLowPtStat1020->Sumw2();
	TH1D*	histoChargedPionMinusSpecLowPtSyst1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent2_pion_minus");
	histoChargedPionMinusSpecLowPtSyst1020->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtStat1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent2_pion_plus");
	histoChargedPionPlusSpecLowPtStat1020->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtSyst1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent2_pion_plus");
	histoChargedPionPlusSpecLowPtSyst1020->Sumw2();
	TH1D*	histoChargedPionSpecLowPtStat1020 = (TH1D*)histoChargedPionMinusSpecLowPtStat1020->Clone("histoChargedPionSpecLowPtStat1020");
	TH1D*	histoChargedPionSpecLowPtSyst1020 = (TH1D*)histoChargedPionMinusSpecLowPtSyst1020->Clone("histoChargedPionSpecLowPtSyst1020");
	histoChargedPionSpecLowPtStat1020->Add(histoChargedPionPlusSpecLowPtStat1020);
	histoChargedPionSpecLowPtSyst1020->Add(histoChargedPionPlusSpecLowPtSyst1020);
	histoChargedPionSpecLowPtStat1020->Scale(0.5);
	histoChargedPionSpecLowPtSyst1020->Scale(0.5);
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst1020->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat1020->SetBinContent(i, histoChargedPionSpecLowPtStat1020->GetBinContent(i)/histoChargedPionSpecLowPtStat1020->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat1020->SetBinError(i, histoChargedPionSpecLowPtStat1020->GetBinError(i)/histoChargedPionSpecLowPtStat1020->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSyst1020->SetBinContent(i, histoChargedPionSpecLowPtSyst1020->GetBinContent(i)/(histoChargedPionSpecLowPtSyst1020->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSyst1020->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst1020->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst1020->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst1020->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst1020->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst1020->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSyst1020->SetBinError(i, histoChargedPionSpecLowPtSyst1020->GetBinContent(i)*fractionalSystematicError/100.);
	}
	
	TH1D*	histoChargedPionMinusSpecLowPtStat2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent3_pion_minus");
	histoChargedPionMinusSpecLowPtStat2030->Sumw2();
	TH1D*	histoChargedPionMinusSpecLowPtSyst2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent3_pion_minus");
	histoChargedPionMinusSpecLowPtSyst2030->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtStat2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent3_pion_plus");
	histoChargedPionPlusSpecLowPtStat2030->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtSyst2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent3_pion_plus");
	histoChargedPionPlusSpecLowPtSyst2030->Sumw2();
	TH1D*	histoChargedPionSpecLowPtStat2030 = (TH1D*)histoChargedPionMinusSpecLowPtStat2030->Clone("histoChargedPionSpecLowPtStat2030");
	TH1D*	histoChargedPionSpecLowPtSyst2030 = (TH1D*)histoChargedPionMinusSpecLowPtSyst2030->Clone("histoChargedPionSpecLowPtSyst2030");
	histoChargedPionSpecLowPtStat2030->Add(histoChargedPionPlusSpecLowPtStat2030);
	histoChargedPionSpecLowPtSyst2030->Add(histoChargedPionPlusSpecLowPtSyst2030);
	histoChargedPionSpecLowPtStat2030->Scale(0.5);
	histoChargedPionSpecLowPtSyst2030->Scale(0.5);
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst2030->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat2030->SetBinContent(i, histoChargedPionSpecLowPtStat2030->GetBinContent(i)/histoChargedPionSpecLowPtStat2030->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat2030->SetBinError(i, histoChargedPionSpecLowPtStat2030->GetBinError(i)/histoChargedPionSpecLowPtStat2030->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSyst2030->SetBinContent(i, histoChargedPionSpecLowPtSyst2030->GetBinContent(i)/(histoChargedPionSpecLowPtSyst2030->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSyst2030->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst2030->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst2030->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst2030->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst2030->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst2030->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSyst2030->SetBinError(i, histoChargedPionSpecLowPtSyst2030->GetBinContent(i)*fractionalSystematicError/100.);
	}
	TH1D*	histoChargedPionMinusSpecLowPtStat3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent4_pion_minus");
	histoChargedPionMinusSpecLowPtStat3040->Sumw2();
	TH1D*	histoChargedPionMinusSpecLowPtSyst3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent4_pion_minus");
	histoChargedPionMinusSpecLowPtSyst3040->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtStat3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent4_pion_plus");
	histoChargedPionPlusSpecLowPtStat3040->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtSyst3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent4_pion_plus");
	histoChargedPionPlusSpecLowPtSyst3040->Sumw2();
	TH1D*	histoChargedPionSpecLowPtStat3040 = (TH1D*)histoChargedPionMinusSpecLowPtStat3040->Clone("histoChargedPionSpecLowPtStat3040");
	TH1D*	histoChargedPionSpecLowPtSyst3040 = (TH1D*)histoChargedPionMinusSpecLowPtSyst3040->Clone("histoChargedPionSpecLowPtSyst3040");
	histoChargedPionSpecLowPtStat3040->Add(histoChargedPionPlusSpecLowPtStat3040);
	histoChargedPionSpecLowPtSyst3040->Add(histoChargedPionPlusSpecLowPtSyst3040);
	histoChargedPionSpecLowPtStat3040->Scale(0.5);
	histoChargedPionSpecLowPtSyst3040->Scale(0.5);
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst3040->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat3040->SetBinContent(i, histoChargedPionSpecLowPtStat3040->GetBinContent(i)/histoChargedPionSpecLowPtStat3040->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat3040->SetBinError(i, histoChargedPionSpecLowPtStat3040->GetBinError(i)/histoChargedPionSpecLowPtStat3040->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSyst3040->SetBinContent(i, histoChargedPionSpecLowPtSyst3040->GetBinContent(i)/(histoChargedPionSpecLowPtSyst3040->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSyst3040->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst3040->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst3040->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst3040->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst3040->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst3040->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSyst3040->SetBinError(i, histoChargedPionSpecLowPtSyst3040->GetBinContent(i)*fractionalSystematicError/100.);
	}
	TH1D*	histoChargedPionSpecLowPtStat2040 = (TH1D*)histoChargedPionSpecLowPtStat2030->Clone("histoChargedPionSpecLowPtStat2040");
	TH1D*	histoChargedPionSpecLowPtSyst2040 = (TH1D*)histoChargedPionSpecLowPtSyst2030->Clone("histoChargedPionSpecLowPtSyst2040");
	histoChargedPionSpecLowPtStat2040->Add(histoChargedPionSpecLowPtStat3040);
	histoChargedPionSpecLowPtSyst2040->Add(histoChargedPionSpecLowPtSyst3040);
	histoChargedPionSpecLowPtStat2040->Scale(0.5);
	histoChargedPionSpecLowPtSyst2040->Scale(0.5);

   for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst2040->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoChargedPionSpecLowPtSyst2030->GetBinContent(i) != 0){
			relErrLowerCent = histoChargedPionSpecLowPtSyst2030->GetBinError(i)/histoChargedPionSpecLowPtSyst2030->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoChargedPionSpecLowPtSyst3040->GetBinContent(i) != 0){
			relErrHigherCent = histoChargedPionSpecLowPtSyst3040->GetBinError(i)/histoChargedPionSpecLowPtSyst3040->GetBinContent(i)*100 ;
		}
		
		if (relErrHigherCent > relErrLowerCent){
			histoChargedPionSpecLowPtSyst2040->SetBinError(i, histoChargedPionSpecLowPtSyst2040->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoChargedPionSpecLowPtSyst2040->SetBinError(i, histoChargedPionSpecLowPtSyst2040->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   
	
	TH1D*	histoChargedPionMinusSpecLowPtStat4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent5_pion_minus");
	histoChargedPionMinusSpecLowPtStat4050->Sumw2();
	TH1D*	histoChargedPionMinusSpecLowPtSyst4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent5_pion_minus");
	histoChargedPionMinusSpecLowPtSyst4050->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtStat4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent5_pion_plus");
	histoChargedPionPlusSpecLowPtStat4050->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtSyst4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent5_pion_plus");
	histoChargedPionPlusSpecLowPtSyst4050->Sumw2();
	TH1D*	histoChargedPionSpecLowPtStat4050 = (TH1D*)histoChargedPionMinusSpecLowPtStat4050->Clone("histoChargedPionSpecLowPtStat4050");
	TH1D*	histoChargedPionSpecLowPtSyst4050 = (TH1D*)histoChargedPionMinusSpecLowPtSyst4050->Clone("histoChargedPionSpecLowPtSyst4050");
	histoChargedPionSpecLowPtStat4050->Add(histoChargedPionPlusSpecLowPtStat4050);
	histoChargedPionSpecLowPtSyst4050->Add(histoChargedPionPlusSpecLowPtSyst4050);
	histoChargedPionSpecLowPtStat4050->Scale(0.5);
	histoChargedPionSpecLowPtSyst4050->Scale(0.5);
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst4050->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat4050->SetBinContent(i, histoChargedPionSpecLowPtStat4050->GetBinContent(i)/histoChargedPionSpecLowPtStat4050->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat4050->SetBinError(i, histoChargedPionSpecLowPtStat4050->GetBinError(i)/histoChargedPionSpecLowPtStat4050->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSyst4050->SetBinContent(i, histoChargedPionSpecLowPtSyst4050->GetBinContent(i)/(histoChargedPionSpecLowPtSyst4050->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSyst4050->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst4050->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst4050->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst4050->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst4050->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst4050->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSyst4050->SetBinError(i, histoChargedPionSpecLowPtSyst4050->GetBinContent(i)*fractionalSystematicError/100.);
	}

	TH1D*	histoChargedPionMinusSpecLowPtStat5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent6_pion_minus");
	histoChargedPionMinusSpecLowPtStat5060->Sumw2();
	TH1D*	histoChargedPionMinusSpecLowPtSyst5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent6_pion_minus");
	histoChargedPionMinusSpecLowPtSyst5060->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtStat5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent6_pion_plus");
	histoChargedPionPlusSpecLowPtStat5060->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtSyst5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent6_pion_plus");
	histoChargedPionPlusSpecLowPtSyst5060->Sumw2();
	TH1D*	histoChargedPionSpecLowPtStat5060 = (TH1D*)histoChargedPionMinusSpecLowPtStat5060->Clone("histoChargedPionSpecLowPtStat5060");
	TH1D*	histoChargedPionSpecLowPtSyst5060 = (TH1D*)histoChargedPionMinusSpecLowPtSyst5060->Clone("histoChargedPionSpecLowPtSyst5060");
	histoChargedPionSpecLowPtStat5060->Add(histoChargedPionPlusSpecLowPtStat5060);
	histoChargedPionSpecLowPtSyst5060->Add(histoChargedPionPlusSpecLowPtSyst5060);
	histoChargedPionSpecLowPtStat5060->Scale(0.5);
	histoChargedPionSpecLowPtSyst5060->Scale(0.5);
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst5060->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat5060->SetBinContent(i, histoChargedPionSpecLowPtStat5060->GetBinContent(i)/histoChargedPionSpecLowPtStat5060->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat5060->SetBinError(i, histoChargedPionSpecLowPtStat5060->GetBinError(i)/histoChargedPionSpecLowPtStat5060->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSyst5060->SetBinContent(i, histoChargedPionSpecLowPtSyst5060->GetBinContent(i)/(histoChargedPionSpecLowPtSyst5060->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSyst5060->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst5060->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst5060->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst5060->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst5060->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst5060->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSyst5060->SetBinError(i, histoChargedPionSpecLowPtSyst5060->GetBinContent(i)*fractionalSystematicError/100.);
	}
	TH1D*	histoChargedPionSpecLowPtStat4060 = (TH1D*)histoChargedPionSpecLowPtStat4050->Clone("histoChargedPionSpecLowPtStat4060");
	TH1D*	histoChargedPionSpecLowPtSyst4060 = (TH1D*)histoChargedPionSpecLowPtSyst4050->Clone("histoChargedPionSpecLowPtSyst4060");
	histoChargedPionSpecLowPtStat4060->Add(histoChargedPionSpecLowPtStat5060);
	histoChargedPionSpecLowPtSyst4060->Add(histoChargedPionSpecLowPtSyst5060);
	histoChargedPionSpecLowPtStat4060->Scale(0.5);
	histoChargedPionSpecLowPtSyst4060->Scale(0.5);

	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst4060->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0.;
		if (histoChargedPionSpecLowPtSyst4050->GetBinContent(i) != 0){
			relErrLowerCent = histoChargedPionSpecLowPtSyst4050->GetBinError(i)/histoChargedPionSpecLowPtSyst4050->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0.;
		if (histoChargedPionSpecLowPtSyst5060->GetBinContent(i) != 0){
			relErrHigherCent = histoChargedPionSpecLowPtSyst5060->GetBinError(i)/histoChargedPionSpecLowPtSyst5060->GetBinContent(i)*100 ;
		}   
		if (relErrHigherCent > relErrLowerCent){
			histoChargedPionSpecLowPtSyst4060->SetBinError(i, histoChargedPionSpecLowPtSyst4060->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoChargedPionSpecLowPtSyst4060->SetBinError(i, histoChargedPionSpecLowPtSyst4060->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   

   
	TH1D*	histoChargedPionMinusSpecLowPtStat6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent7_pion_minus");
	histoChargedPionMinusSpecLowPtStat6070->Sumw2();
	TH1D*	histoChargedPionMinusSpecLowPtSyst6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent7_pion_minus");
	histoChargedPionMinusSpecLowPtSyst6070->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtStat6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent7_pion_plus");
	histoChargedPionPlusSpecLowPtStat6070->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtSyst6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent7_pion_plus");
	histoChargedPionPlusSpecLowPtSyst6070->Sumw2();
	TH1D*	histoChargedPionSpecLowPtStat6070 = (TH1D*)histoChargedPionMinusSpecLowPtStat6070->Clone("histoChargedPionSpecLowPtStat6070");
	TH1D*	histoChargedPionSpecLowPtSyst6070 = (TH1D*)histoChargedPionMinusSpecLowPtSyst6070->Clone("histoChargedPionSpecLowPtSyst6070");
	histoChargedPionSpecLowPtStat6070->Add(histoChargedPionPlusSpecLowPtStat6070);
	histoChargedPionSpecLowPtSyst6070->Add(histoChargedPionPlusSpecLowPtSyst6070);
	histoChargedPionSpecLowPtStat6070->Scale(0.5);
	histoChargedPionSpecLowPtSyst6070->Scale(0.5);
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst6070->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat6070->SetBinContent(i, histoChargedPionSpecLowPtStat6070->GetBinContent(i)/histoChargedPionSpecLowPtStat6070->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat6070->SetBinError(i, histoChargedPionSpecLowPtStat6070->GetBinError(i)/histoChargedPionSpecLowPtStat6070->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSyst6070->SetBinContent(i, histoChargedPionSpecLowPtSyst6070->GetBinContent(i)/(histoChargedPionSpecLowPtSyst6070->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSyst6070->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst6070->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst6070->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst6070->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst6070->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst6070->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSyst6070->SetBinError(i, histoChargedPionSpecLowPtSyst6070->GetBinContent(i)*fractionalSystematicError/100.);
	}

	TH1D*	histoChargedPionMinusSpecLowPtStat7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent8_pion_minus");
	histoChargedPionMinusSpecLowPtStat7080->Sumw2();
	TH1D*	histoChargedPionMinusSpecLowPtSyst7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent8_pion_minus");
	histoChargedPionMinusSpecLowPtSyst7080->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtStat7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent8_pion_plus");
	histoChargedPionPlusSpecLowPtStat7080->Sumw2();
	TH1D*	histoChargedPionPlusSpecLowPtSyst7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent8_pion_plus");
	histoChargedPionPlusSpecLowPtSyst7080->Sumw2();
	TH1D*	histoChargedPionSpecLowPtStat7080 = (TH1D*)histoChargedPionMinusSpecLowPtStat7080->Clone("histoChargedPionSpecLowPtStat7080");
	TH1D*	histoChargedPionSpecLowPtSyst7080 = (TH1D*)histoChargedPionMinusSpecLowPtSyst7080->Clone("histoChargedPionSpecLowPtSyst7080");
	histoChargedPionSpecLowPtStat7080->Add(histoChargedPionPlusSpecLowPtStat7080);
	histoChargedPionSpecLowPtSyst7080->Add(histoChargedPionPlusSpecLowPtSyst7080);
	histoChargedPionSpecLowPtStat7080->Scale(0.5);
	histoChargedPionSpecLowPtSyst7080->Scale(0.5);
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst7080->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStat7080->SetBinContent(i, histoChargedPionSpecLowPtStat7080->GetBinContent(i)/histoChargedPionSpecLowPtStat7080->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStat7080->SetBinError(i, histoChargedPionSpecLowPtStat7080->GetBinError(i)/histoChargedPionSpecLowPtStat7080->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtSyst7080->SetBinContent(i, histoChargedPionSpecLowPtSyst7080->GetBinContent(i)/(histoChargedPionSpecLowPtSyst7080->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedPionPlusSpecLowPtSyst7080->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst7080->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst7080->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst7080->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst7080->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst7080->GetBinContent(i)*100)/2;
		}
		histoChargedPionSpecLowPtSyst7080->SetBinError(i, histoChargedPionSpecLowPtSyst7080->GetBinContent(i)*fractionalSystematicError/100.);
	}

	TH1D*	histoChargedPionSpecLowPtStat6080 = (TH1D*)histoChargedPionSpecLowPtStat6070->Clone("histoChargedPionSpecLowPtStat6080");
	TH1D*	histoChargedPionSpecLowPtSyst6080 = (TH1D*)histoChargedPionSpecLowPtSyst6070->Clone("histoChargedPionSpecLowPtSyst6080");
	histoChargedPionSpecLowPtStat6080->Add(histoChargedPionSpecLowPtStat7080);
	histoChargedPionSpecLowPtSyst6080->Add(histoChargedPionSpecLowPtSyst7080);
	histoChargedPionSpecLowPtStat6080->Scale(0.5);
	histoChargedPionSpecLowPtSyst6080->Scale(0.5);

	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst6080->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoChargedPionSpecLowPtSyst6070->GetBinContent(i) != 0){
			relErrLowerCent = histoChargedPionSpecLowPtSyst6070->GetBinError(i)/histoChargedPionSpecLowPtSyst6070->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoChargedPionSpecLowPtSyst7080->GetBinContent(i) != 0){
			relErrHigherCent = histoChargedPionSpecLowPtSyst7080->GetBinError(i)/histoChargedPionSpecLowPtSyst7080->GetBinContent(i)*100 ;
		}
		if (relErrHigherCent > relErrLowerCent){
			histoChargedPionSpecLowPtSyst6080->SetBinError(i, histoChargedPionSpecLowPtSyst6080->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoChargedPionSpecLowPtSyst6080->SetBinError(i, histoChargedPionSpecLowPtSyst6080->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   

   
	// ****************************************************************************************************************
	// ******************************************charged Kaons ********************************************************
	// ****************************************************************************************************************
	
	TH1D* histoChargedKaonMinusSpecLowPtStat0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent0_kaon_minus");
	histoChargedKaonMinusSpecLowPtStat0005->Sumw2();
	TH1D* histoChargedKaonMinusSpecLowPtSyst0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent0_kaon_minus");
	histoChargedKaonMinusSpecLowPtSyst0005->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtStat0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent0_kaon_plus");
	histoChargedKaonPlusSpecLowPtStat0005->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtSyst0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent0_kaon_plus");
	histoChargedKaonPlusSpecLowPtSyst0005->Sumw2();
	TH1D* histoChargedKaonSpecLowPtStat0005 = (TH1D*)histoChargedKaonMinusSpecLowPtStat0005->Clone("histoChargedKaonSpecLowPtStat0005");
	TH1D* histoChargedKaonSpecLowPtSyst0005 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst0005->Clone("histoChargedKaonSpecLowPtSyst0005");
	histoChargedKaonSpecLowPtStat0005->Add(histoChargedKaonPlusSpecLowPtStat0005);
	histoChargedKaonSpecLowPtSyst0005->Add(histoChargedKaonPlusSpecLowPtSyst0005);
	histoChargedKaonSpecLowPtStat0005->Scale(0.5);
	histoChargedKaonSpecLowPtSyst0005->Scale(0.5);
		
	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst0005->GetNbinsX()+1; i++){
		histoChargedKaonSpecLowPtStat0005->SetBinContent(i, histoChargedKaonSpecLowPtStat0005->GetBinContent(i)/histoChargedKaonSpecLowPtStat0005->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtStat0005->SetBinError(i, histoChargedKaonSpecLowPtStat0005->GetBinError(i)/histoChargedKaonSpecLowPtStat0005->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtSyst0005->SetBinContent(i, histoChargedKaonSpecLowPtSyst0005->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst0005->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedKaonMinusSpecLowPtSyst0005->GetBinContent(i) != 0 && histoChargedKaonPlusSpecLowPtSyst0005->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst0005->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst0005->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst0005->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst0005->GetBinContent(i)*100)/2;
		}
		histoChargedKaonSpecLowPtSyst0005->SetBinError(i, histoChargedKaonSpecLowPtSyst0005->GetBinContent(i)*fractionalSystematicError/100.);
	}  

	TH1D* histoChargedKaonMinusSpecLowPtStat0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent1_kaon_minus");
	histoChargedKaonMinusSpecLowPtStat0510->Sumw2();
	TH1D* histoChargedKaonMinusSpecLowPtSyst0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent1_kaon_minus");
	histoChargedKaonMinusSpecLowPtSyst0510->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtStat0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent1_kaon_plus");
	histoChargedKaonPlusSpecLowPtStat0510->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtSyst0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent1_kaon_plus");
	histoChargedKaonPlusSpecLowPtSyst0510->Sumw2();
	TH1D* histoChargedKaonSpecLowPtStat0510 = (TH1D*)histoChargedKaonMinusSpecLowPtStat0510->Clone("histoChargedKaonSpecLowPtStat0510");
	TH1D* histoChargedKaonSpecLowPtSyst0510 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst0510->Clone("histoChargedKaonSpecLowPtSyst0510");
	histoChargedKaonSpecLowPtStat0510->Add(histoChargedKaonPlusSpecLowPtStat0510);
	histoChargedKaonSpecLowPtSyst0510->Add(histoChargedKaonPlusSpecLowPtSyst0510);
	histoChargedKaonSpecLowPtStat0510->Scale(0.5);
	histoChargedKaonSpecLowPtSyst0510->Scale(0.5);
	
	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst0510->GetNbinsX()+1; i++){
		histoChargedKaonSpecLowPtStat0510->SetBinContent(i, histoChargedKaonSpecLowPtStat0510->GetBinContent(i)/histoChargedKaonSpecLowPtStat0510->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtStat0510->SetBinError(i, histoChargedKaonSpecLowPtStat0510->GetBinError(i)/histoChargedKaonSpecLowPtStat0510->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtSyst0510->SetBinContent(i, histoChargedKaonSpecLowPtSyst0510->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst0510->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedKaonPlusSpecLowPtSyst0510->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst0510->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst0510->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst0510->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst0510->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst0510->GetBinContent(i)*100)/2;
		}
		histoChargedKaonSpecLowPtSyst0510->SetBinError(i, histoChargedKaonSpecLowPtSyst0510->GetBinContent(i)*fractionalSystematicError/100.);
	}
	
	TH1D* histoChargedKaonSpecLowPtStat0010 = (TH1D*)histoChargedKaonSpecLowPtStat0510->Clone("histoChargedKaonSpecLowPtStat0010");
	TH1D* histoChargedKaonSpecLowPtSyst0010 = (TH1D*)histoChargedKaonSpecLowPtSyst0510->Clone("histoChargedKaonSpecLowPtSyst0010");
	histoChargedKaonSpecLowPtStat0010->Add(histoChargedKaonSpecLowPtStat0005);
	histoChargedKaonSpecLowPtSyst0010->Add(histoChargedKaonSpecLowPtSyst0005);
	histoChargedKaonSpecLowPtStat0010->Scale(0.5);
	histoChargedKaonSpecLowPtSyst0010->Scale(0.5);

	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoChargedKaonSpecLowPtSyst0005->GetBinContent(i) != 0){
			relErrLowerCent = histoChargedKaonSpecLowPtSyst0005->GetBinError(i)/histoChargedKaonSpecLowPtSyst0005->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoChargedKaonSpecLowPtSyst0510->GetBinContent(i) != 0){
			relErrHigherCent = histoChargedKaonSpecLowPtSyst0510->GetBinError(i)/histoChargedKaonSpecLowPtSyst0510->GetBinContent(i)*100 ;
		}
		if (relErrHigherCent > relErrLowerCent){
			histoChargedKaonSpecLowPtSyst0010->SetBinError(i, histoChargedKaonSpecLowPtSyst0010->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoChargedKaonSpecLowPtSyst0010->SetBinError(i, histoChargedKaonSpecLowPtSyst0010->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   

	TH1D* histoChargedKaonMinusSpecLowPtStat1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent2_kaon_minus");
	histoChargedKaonMinusSpecLowPtStat1020->Sumw2();
	TH1D* histoChargedKaonMinusSpecLowPtSyst1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent2_kaon_minus");
	histoChargedKaonMinusSpecLowPtSyst1020->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtStat1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent2_kaon_plus");
	histoChargedKaonPlusSpecLowPtStat1020->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtSyst1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent2_kaon_plus");
	histoChargedKaonPlusSpecLowPtSyst1020->Sumw2();
	TH1D* histoChargedKaonSpecLowPtStat1020 = (TH1D*)histoChargedKaonMinusSpecLowPtStat1020->Clone("histoChargedKaonSpecLowPtStat1020");
	TH1D* histoChargedKaonSpecLowPtSyst1020 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst1020->Clone("histoChargedKaonSpecLowPtSyst1020");
	histoChargedKaonSpecLowPtStat1020->Add(histoChargedKaonPlusSpecLowPtStat1020);
	histoChargedKaonSpecLowPtSyst1020->Add(histoChargedKaonPlusSpecLowPtSyst1020);
	histoChargedKaonSpecLowPtStat1020->Scale(0.5);
	histoChargedKaonSpecLowPtSyst1020->Scale(0.5);
	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst1020->GetNbinsX()+1; i++){
		histoChargedKaonSpecLowPtStat1020->SetBinContent(i, histoChargedKaonSpecLowPtStat1020->GetBinContent(i)/histoChargedKaonSpecLowPtStat1020->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtStat1020->SetBinError(i, histoChargedKaonSpecLowPtStat1020->GetBinError(i)/histoChargedKaonSpecLowPtStat1020->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtSyst1020->SetBinContent(i, histoChargedKaonSpecLowPtSyst1020->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst1020->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedKaonPlusSpecLowPtSyst1020->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst1020->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst1020->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst1020->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst1020->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst1020->GetBinContent(i)*100)/2;
		}
		histoChargedKaonSpecLowPtSyst1020->SetBinError(i, histoChargedKaonSpecLowPtSyst1020->GetBinContent(i)*fractionalSystematicError/100.);
	}
	
	TH1D* histoChargedKaonMinusSpecLowPtStat2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent3_kaon_minus");
	histoChargedKaonMinusSpecLowPtStat2030->Sumw2();
	TH1D* histoChargedKaonMinusSpecLowPtSyst2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent3_kaon_minus");
	histoChargedKaonMinusSpecLowPtSyst2030->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtStat2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent3_kaon_plus");
	histoChargedKaonPlusSpecLowPtStat2030->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtSyst2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent3_kaon_plus");
	histoChargedKaonPlusSpecLowPtSyst2030->Sumw2();
	TH1D* histoChargedKaonSpecLowPtStat2030 = (TH1D*)histoChargedKaonMinusSpecLowPtStat2030->Clone("histoChargedKaonSpecLowPtStat2030");
	TH1D* histoChargedKaonSpecLowPtSyst2030 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst2030->Clone("histoChargedKaonSpecLowPtSyst2030");
	histoChargedKaonSpecLowPtStat2030->Add(histoChargedKaonPlusSpecLowPtStat2030);
	histoChargedKaonSpecLowPtSyst2030->Add(histoChargedKaonPlusSpecLowPtSyst2030);
	histoChargedKaonSpecLowPtStat2030->Scale(0.5);
	histoChargedKaonSpecLowPtSyst2030->Scale(0.5);
	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst2030->GetNbinsX()+1; i++){
		histoChargedKaonSpecLowPtStat2030->SetBinContent(i, histoChargedKaonSpecLowPtStat2030->GetBinContent(i)/histoChargedKaonSpecLowPtStat2030->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtStat2030->SetBinError(i, histoChargedKaonSpecLowPtStat2030->GetBinError(i)/histoChargedKaonSpecLowPtStat2030->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtSyst2030->SetBinContent(i, histoChargedKaonSpecLowPtSyst2030->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst2030->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedKaonPlusSpecLowPtSyst2030->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst2030->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst2030->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst2030->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst2030->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst2030->GetBinContent(i)*100)/2;
		}
		histoChargedKaonSpecLowPtSyst2030->SetBinError(i, histoChargedKaonSpecLowPtSyst2030->GetBinContent(i)*fractionalSystematicError/100.);
	}
	TH1D* histoChargedKaonMinusSpecLowPtStat3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent4_kaon_minus");
	histoChargedKaonMinusSpecLowPtStat3040->Sumw2();
	TH1D* histoChargedKaonMinusSpecLowPtSyst3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent4_kaon_minus");
	histoChargedKaonMinusSpecLowPtSyst3040->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtStat3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent4_kaon_plus");
	histoChargedKaonPlusSpecLowPtStat3040->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtSyst3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent4_kaon_plus");
	histoChargedKaonPlusSpecLowPtSyst3040->Sumw2();
	TH1D* histoChargedKaonSpecLowPtStat3040 = (TH1D*)histoChargedKaonMinusSpecLowPtStat3040->Clone("histoChargedKaonSpecLowPtStat3040");
	TH1D* histoChargedKaonSpecLowPtSyst3040 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst3040->Clone("histoChargedKaonSpecLowPtSyst3040");
	histoChargedKaonSpecLowPtStat3040->Add(histoChargedKaonPlusSpecLowPtStat3040);
	histoChargedKaonSpecLowPtSyst3040->Add(histoChargedKaonPlusSpecLowPtSyst3040);
	histoChargedKaonSpecLowPtStat3040->Scale(0.5);
	histoChargedKaonSpecLowPtSyst3040->Scale(0.5);
	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst3040->GetNbinsX()+1; i++){
		histoChargedKaonSpecLowPtStat3040->SetBinContent(i, histoChargedKaonSpecLowPtStat3040->GetBinContent(i)/histoChargedKaonSpecLowPtStat3040->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtStat3040->SetBinError(i, histoChargedKaonSpecLowPtStat3040->GetBinError(i)/histoChargedKaonSpecLowPtStat3040->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtSyst3040->SetBinContent(i, histoChargedKaonSpecLowPtSyst3040->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst3040->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedKaonPlusSpecLowPtSyst3040->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst3040->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst3040->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst3040->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst3040->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst3040->GetBinContent(i)*100)/2;
		}
		histoChargedKaonSpecLowPtSyst3040->SetBinError(i, histoChargedKaonSpecLowPtSyst3040->GetBinContent(i)*fractionalSystematicError/100.);
	}
	TH1D* histoChargedKaonSpecLowPtStat2040 = (TH1D*)histoChargedKaonSpecLowPtStat2030->Clone("histoChargedKaonSpecLowPtStat2040");
	TH1D* histoChargedKaonSpecLowPtSyst2040 = (TH1D*)histoChargedKaonSpecLowPtSyst2030->Clone("histoChargedKaonSpecLowPtSyst2040");
	histoChargedKaonSpecLowPtStat2040->Add(histoChargedKaonSpecLowPtStat3040);
	histoChargedKaonSpecLowPtSyst2040->Add(histoChargedKaonSpecLowPtSyst3040);
	histoChargedKaonSpecLowPtStat2040->Scale(0.5);
	histoChargedKaonSpecLowPtSyst2040->Scale(0.5);

	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst2040->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoChargedKaonSpecLowPtSyst2030->GetBinContent(i) !=0){
			relErrLowerCent = histoChargedKaonSpecLowPtSyst2030->GetBinError(i)/histoChargedKaonSpecLowPtSyst2030->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoChargedKaonSpecLowPtSyst3040->GetBinContent(i) != 0){
			relErrHigherCent = histoChargedKaonSpecLowPtSyst3040->GetBinError(i)/histoChargedKaonSpecLowPtSyst3040->GetBinContent(i)*100 ;
		}
		if (relErrHigherCent > relErrLowerCent){
			histoChargedKaonSpecLowPtSyst2040->SetBinError(i, histoChargedKaonSpecLowPtSyst2040->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoChargedKaonSpecLowPtSyst2040->SetBinError(i, histoChargedKaonSpecLowPtSyst2040->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   
	
	TH1D* histoChargedKaonMinusSpecLowPtStat4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent5_kaon_minus");
	histoChargedKaonMinusSpecLowPtStat4050->Sumw2();
	TH1D* histoChargedKaonMinusSpecLowPtSyst4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent5_kaon_minus");
	histoChargedKaonMinusSpecLowPtSyst4050->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtStat4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent5_kaon_plus");
	histoChargedKaonPlusSpecLowPtStat4050->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtSyst4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent5_kaon_plus");
	histoChargedKaonPlusSpecLowPtSyst4050->Sumw2();
	TH1D* histoChargedKaonSpecLowPtStat4050 = (TH1D*)histoChargedKaonMinusSpecLowPtStat4050->Clone("histoChargedKaonSpecLowPtStat4050");
	TH1D* histoChargedKaonSpecLowPtSyst4050 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst4050->Clone("histoChargedKaonSpecLowPtSyst4050");
	histoChargedKaonSpecLowPtStat4050->Add(histoChargedKaonPlusSpecLowPtStat4050);
	histoChargedKaonSpecLowPtSyst4050->Add(histoChargedKaonPlusSpecLowPtSyst4050);
	histoChargedKaonSpecLowPtStat4050->Scale(0.5);
	histoChargedKaonSpecLowPtSyst4050->Scale(0.5);
	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst4050->GetNbinsX()+1; i++){
		histoChargedKaonSpecLowPtStat4050->SetBinContent(i, histoChargedKaonSpecLowPtStat4050->GetBinContent(i)/histoChargedKaonSpecLowPtStat4050->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtStat4050->SetBinError(i, histoChargedKaonSpecLowPtStat4050->GetBinError(i)/histoChargedKaonSpecLowPtStat4050->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtSyst4050->SetBinContent(i, histoChargedKaonSpecLowPtSyst4050->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst4050->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedKaonPlusSpecLowPtSyst4050->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst4050->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst4050->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst4050->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst4050->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst4050->GetBinContent(i)*100)/2;
		}
		histoChargedKaonSpecLowPtSyst4050->SetBinError(i, histoChargedKaonSpecLowPtSyst4050->GetBinContent(i)*fractionalSystematicError/100.);
	}

	TH1D* histoChargedKaonMinusSpecLowPtStat5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent6_kaon_minus");
	histoChargedKaonMinusSpecLowPtStat5060->Sumw2();
	TH1D* histoChargedKaonMinusSpecLowPtSyst5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent6_kaon_minus");
	histoChargedKaonMinusSpecLowPtSyst5060->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtStat5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent6_kaon_plus");
	histoChargedKaonPlusSpecLowPtStat5060->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtSyst5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent6_kaon_plus");
	histoChargedKaonPlusSpecLowPtSyst5060->Sumw2();
	TH1D* histoChargedKaonSpecLowPtStat5060 = (TH1D*)histoChargedKaonMinusSpecLowPtStat5060->Clone("histoChargedKaonSpecLowPtStat5060");
	TH1D* histoChargedKaonSpecLowPtSyst5060 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst5060->Clone("histoChargedKaonSpecLowPtSyst5060");
	histoChargedKaonSpecLowPtStat5060->Add(histoChargedKaonPlusSpecLowPtStat5060);
	histoChargedKaonSpecLowPtSyst5060->Add(histoChargedKaonPlusSpecLowPtSyst5060);
	histoChargedKaonSpecLowPtStat5060->Scale(0.5);
	histoChargedKaonSpecLowPtSyst5060->Scale(0.5);
	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst5060->GetNbinsX()+1; i++){
		histoChargedKaonSpecLowPtStat5060->SetBinContent(i, histoChargedKaonSpecLowPtStat5060->GetBinContent(i)/histoChargedKaonSpecLowPtStat5060->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtStat5060->SetBinError(i, histoChargedKaonSpecLowPtStat5060->GetBinError(i)/histoChargedKaonSpecLowPtStat5060->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtSyst5060->SetBinContent(i, histoChargedKaonSpecLowPtSyst5060->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst5060->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedKaonPlusSpecLowPtSyst5060->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst5060->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst5060->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst5060->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst5060->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst5060->GetBinContent(i)*100)/2;
		}
		histoChargedKaonSpecLowPtSyst5060->SetBinError(i, histoChargedKaonSpecLowPtSyst5060->GetBinContent(i)*fractionalSystematicError/100.);
	}
	TH1D* histoChargedKaonSpecLowPtStat4060 = (TH1D*)histoChargedKaonSpecLowPtStat4050->Clone("histoChargedKaonSpecLowPtStat4060");
	TH1D* histoChargedKaonSpecLowPtSyst4060 = (TH1D*)histoChargedKaonSpecLowPtSyst4050->Clone("histoChargedKaonSpecLowPtSyst4060");
	histoChargedKaonSpecLowPtStat4060->Add(histoChargedKaonSpecLowPtStat5060);
	histoChargedKaonSpecLowPtSyst4060->Add(histoChargedKaonSpecLowPtSyst5060);
	histoChargedKaonSpecLowPtStat4060->Scale(0.5);
	histoChargedKaonSpecLowPtSyst4060->Scale(0.5);

	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst4060->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoChargedKaonSpecLowPtSyst4050->GetBinContent(i) !=0 ){
			relErrLowerCent =histoChargedKaonSpecLowPtSyst4050->GetBinError(i)/histoChargedKaonSpecLowPtSyst4050->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent =0;
		if (histoChargedKaonSpecLowPtSyst5060->GetBinContent(i) !=0 ){
			relErrHigherCent = histoChargedKaonSpecLowPtSyst5060->GetBinError(i)/histoChargedKaonSpecLowPtSyst5060->GetBinContent(i)*100 ;
		}
		if (relErrHigherCent > relErrLowerCent){
			histoChargedKaonSpecLowPtSyst4060->SetBinError(i, histoChargedKaonSpecLowPtSyst4060->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoChargedKaonSpecLowPtSyst4060->SetBinError(i, histoChargedKaonSpecLowPtSyst4060->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   

	
	TH1D* histoChargedKaonMinusSpecLowPtStat6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent7_kaon_minus");
	histoChargedKaonMinusSpecLowPtStat6070->Sumw2();
	TH1D* histoChargedKaonMinusSpecLowPtSyst6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent7_kaon_minus");
	histoChargedKaonMinusSpecLowPtSyst6070->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtStat6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent7_kaon_plus");
	histoChargedKaonPlusSpecLowPtStat6070->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtSyst6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent7_kaon_plus");
	histoChargedKaonPlusSpecLowPtSyst6070->Sumw2();
	TH1D* histoChargedKaonSpecLowPtStat6070 = (TH1D*)histoChargedKaonMinusSpecLowPtStat6070->Clone("histoChargedKaonSpecLowPtStat6070");
	TH1D* histoChargedKaonSpecLowPtSyst6070 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst6070->Clone("histoChargedKaonSpecLowPtSyst6070");
	histoChargedKaonSpecLowPtStat6070->Add(histoChargedKaonPlusSpecLowPtStat6070);
	histoChargedKaonSpecLowPtSyst6070->Add(histoChargedKaonPlusSpecLowPtSyst6070);
	histoChargedKaonSpecLowPtStat6070->Scale(0.5);
	histoChargedKaonSpecLowPtSyst6070->Scale(0.5);
	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst6070->GetNbinsX()+1; i++){
		histoChargedKaonSpecLowPtStat6070->SetBinContent(i, histoChargedKaonSpecLowPtStat6070->GetBinContent(i)/histoChargedKaonSpecLowPtStat6070->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtStat6070->SetBinError(i, histoChargedKaonSpecLowPtStat6070->GetBinError(i)/histoChargedKaonSpecLowPtStat6070->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtSyst6070->SetBinContent(i, histoChargedKaonSpecLowPtSyst6070->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst6070->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedKaonPlusSpecLowPtSyst6070->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst6070->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst6070->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst6070->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst6070->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst6070->GetBinContent(i)*100)/2;
		}
		histoChargedKaonSpecLowPtSyst6070->SetBinError(i, histoChargedKaonSpecLowPtSyst6070->GetBinContent(i)*fractionalSystematicError/100.);
	}

	TH1D* histoChargedKaonMinusSpecLowPtStat7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent8_kaon_minus");
	histoChargedKaonMinusSpecLowPtStat7080->Sumw2();
	TH1D* histoChargedKaonMinusSpecLowPtSyst7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent8_kaon_minus");
	histoChargedKaonMinusSpecLowPtSyst7080->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtStat7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent8_kaon_plus");
	histoChargedKaonPlusSpecLowPtStat7080->Sumw2();
	TH1D* histoChargedKaonPlusSpecLowPtSyst7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent8_kaon_plus");
	histoChargedKaonPlusSpecLowPtSyst7080->Sumw2();
	TH1D* histoChargedKaonSpecLowPtStat7080 = (TH1D*)histoChargedKaonMinusSpecLowPtStat7080->Clone("histoChargedKaonSpecLowPtStat7080");
	TH1D* histoChargedKaonSpecLowPtSyst7080 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst7080->Clone("histoChargedKaonSpecLowPtSyst7080");
	histoChargedKaonSpecLowPtStat7080->Add(histoChargedKaonPlusSpecLowPtStat7080);
	histoChargedKaonSpecLowPtSyst7080->Add(histoChargedKaonPlusSpecLowPtSyst7080);
	histoChargedKaonSpecLowPtStat7080->Scale(0.5);
	histoChargedKaonSpecLowPtSyst7080->Scale(0.5);
	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst7080->GetNbinsX()+1; i++){
		histoChargedKaonSpecLowPtStat7080->SetBinContent(i, histoChargedKaonSpecLowPtStat7080->GetBinContent(i)/histoChargedKaonSpecLowPtStat7080->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtStat7080->SetBinError(i, histoChargedKaonSpecLowPtStat7080->GetBinError(i)/histoChargedKaonSpecLowPtStat7080->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedKaonSpecLowPtSyst7080->SetBinContent(i, histoChargedKaonSpecLowPtSyst7080->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst7080->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoChargedKaonPlusSpecLowPtSyst7080->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst7080->GetBinContent(i) != 0){
			fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst7080->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst7080->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst7080->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst7080->GetBinContent(i)*100)/2;
		}
		histoChargedKaonSpecLowPtSyst7080->SetBinError(i, histoChargedKaonSpecLowPtSyst7080->GetBinContent(i)*fractionalSystematicError/100.);
	}

	TH1D* histoChargedKaonSpecLowPtStat6080 = (TH1D*)histoChargedKaonSpecLowPtStat6070->Clone("histoChargedKaonSpecLowPtStat6080");
	TH1D* histoChargedKaonSpecLowPtSyst6080 = (TH1D*)histoChargedKaonSpecLowPtSyst6070->Clone("histoChargedKaonSpecLowPtSyst6080");
	histoChargedKaonSpecLowPtStat6080->Add(histoChargedKaonSpecLowPtStat7080);
	histoChargedKaonSpecLowPtSyst6080->Add(histoChargedKaonSpecLowPtSyst7080);
	histoChargedKaonSpecLowPtStat6080->Scale(0.5);
	histoChargedKaonSpecLowPtSyst6080->Scale(0.5);

	for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst6080->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoChargedKaonSpecLowPtSyst6070->GetBinContent(i) != 0){
			relErrLowerCent = histoChargedKaonSpecLowPtSyst6070->GetBinError(i)/histoChargedKaonSpecLowPtSyst6070->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoChargedKaonSpecLowPtSyst7080->GetBinContent(i) != 0){
			relErrHigherCent = histoChargedKaonSpecLowPtSyst7080->GetBinError(i)/histoChargedKaonSpecLowPtSyst7080->GetBinContent(i)*100 ;
		}
		if (relErrHigherCent > relErrLowerCent){
			histoChargedKaonSpecLowPtSyst6080->SetBinError(i, histoChargedKaonSpecLowPtSyst6080->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoChargedKaonSpecLowPtSyst6080->SetBinError(i, histoChargedKaonSpecLowPtSyst6080->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   

	TFile* fileK0sFinalPbPb = new TFile("ExternalInputPbPb/NeutralKaon/k0s_lambda_final_spectra_12112013.root");
	TH1F* histoNeutralKaonSpecStat0005 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent0005_K0s");
	histoNeutralKaonSpecStat0005->Sumw2();
	TH1F* histoNeutralKaonSpecSyst0005 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent0005_K0s");
	histoNeutralKaonSpecSyst0005->Sumw2();
	for (Int_t i = 1; i < histoNeutralKaonSpecSyst0005->GetNbinsX()+1; i++){
		histoNeutralKaonSpecStat0005->SetBinContent(i, histoNeutralKaonSpecStat0005->GetBinContent(i)/histoNeutralKaonSpecStat0005->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecStat0005->SetBinError(i, histoNeutralKaonSpecStat0005->GetBinError(i)/histoNeutralKaonSpecStat0005->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecSyst0005->SetBinContent(i, histoNeutralKaonSpecSyst0005->GetBinContent(i)/(histoNeutralKaonSpecSyst0005->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoNeutralKaonSpecSyst0005->GetBinContent(i) != 0){
			fractionalSystematicError = histoNeutralKaonSpecSyst0005->GetBinError(i)/histoNeutralKaonSpecSyst0005->GetBinContent(i)*100;
		}
		histoNeutralKaonSpecSyst0005->SetBinError(i, histoNeutralKaonSpecSyst0005->GetBinContent(i)*fractionalSystematicError/100.);
	}
	
	TH1F* histoNeutralKaonSpecStat0510 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent0510_K0s");
	histoNeutralKaonSpecStat0510->Sumw2();
	TH1F* histoNeutralKaonSpecSyst0510 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent0510_K0s");
	histoNeutralKaonSpecSyst0510->Sumw2();
	for (Int_t i = 1; i < histoNeutralKaonSpecSyst0510->GetNbinsX()+1; i++){
		histoNeutralKaonSpecStat0510->SetBinContent(i, histoNeutralKaonSpecStat0510->GetBinContent(i)/histoNeutralKaonSpecStat0510->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecStat0510->SetBinError(i, histoNeutralKaonSpecStat0510->GetBinError(i)/histoNeutralKaonSpecStat0510->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecSyst0510->SetBinContent(i, histoNeutralKaonSpecSyst0510->GetBinContent(i)/(histoNeutralKaonSpecSyst0510->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoNeutralKaonSpecSyst0510->GetBinContent(i) != 0){
			fractionalSystematicError = histoNeutralKaonSpecSyst0510->GetBinError(i)/histoNeutralKaonSpecSyst0510->GetBinContent(i)*100;
		}
		histoNeutralKaonSpecSyst0510->SetBinError(i, histoNeutralKaonSpecSyst0510->GetBinContent(i)*fractionalSystematicError/100.);
	}

	TH1D* histoNeutralKaonSpecStat0010 = (TH1D*)histoNeutralKaonSpecStat0005->Clone("histoNeutralKaonSpecStat0010");
	TH1D* histoNeutralKaonSpecSyst0010 = (TH1D*)histoNeutralKaonSpecSyst0005->Clone("histoChargedKaonSpecLowPtSyst0010");
	histoNeutralKaonSpecStat0010->Add(histoNeutralKaonSpecStat0510);
	histoNeutralKaonSpecSyst0010->Add(histoNeutralKaonSpecSyst0510);
	histoNeutralKaonSpecStat0010->Scale(0.5);
	histoNeutralKaonSpecSyst0010->Scale(0.5);

	for (Int_t i = 1; i < histoNeutralKaonSpecSyst0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoNeutralKaonSpecSyst0005->GetBinContent(i) != 0){
			relErrLowerCent = histoNeutralKaonSpecSyst0005->GetBinError(i)/histoNeutralKaonSpecSyst0005->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoNeutralKaonSpecSyst0510->GetBinContent(i) != 0){
			relErrHigherCent = histoNeutralKaonSpecSyst0510->GetBinError(i)/histoNeutralKaonSpecSyst0510->GetBinContent(i)*100 ;
		}
		if (relErrHigherCent > relErrLowerCent){
			histoNeutralKaonSpecSyst0010->SetBinError(i, histoNeutralKaonSpecSyst0010->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoNeutralKaonSpecSyst0010->SetBinError(i, histoNeutralKaonSpecSyst0010->GetBinContent(i)*relErrLowerCent/100);
		}         
	}
	
	TH1F* histoNeutralKaonSpecStat1020 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent1020_K0s");
	histoNeutralKaonSpecStat1020->Sumw2();
	TH1F* histoNeutralKaonSpecSyst1020 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent1020_K0s");
	histoNeutralKaonSpecSyst1020->Sumw2();
	for (Int_t i = 1; i < histoNeutralKaonSpecSyst1020->GetNbinsX()+1; i++){
		histoNeutralKaonSpecStat1020->SetBinContent(i, histoNeutralKaonSpecStat1020->GetBinContent(i)/histoNeutralKaonSpecStat1020->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecStat1020->SetBinError(i, histoNeutralKaonSpecStat1020->GetBinError(i)/histoNeutralKaonSpecStat1020->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecSyst1020->SetBinContent(i, histoNeutralKaonSpecSyst1020->GetBinContent(i)/(histoNeutralKaonSpecSyst1020->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoNeutralKaonSpecSyst1020->GetBinContent(i) != 0){
			fractionalSystematicError = histoNeutralKaonSpecSyst1020->GetBinError(i)/histoNeutralKaonSpecSyst1020->GetBinContent(i)*100;
		}
		histoNeutralKaonSpecSyst1020->SetBinError(i, histoNeutralKaonSpecSyst1020->GetBinContent(i)*fractionalSystematicError/100.);
	}

	TH1F* histoNeutralKaonSpecStat2040 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent2040_K0s");
	histoNeutralKaonSpecStat2040->Sumw2();
	TH1F* histoNeutralKaonSpecSyst2040 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent2040_K0s");
	histoNeutralKaonSpecSyst2040->Sumw2();
	for (Int_t i = 1; i < histoNeutralKaonSpecSyst2040->GetNbinsX()+1; i++){
		histoNeutralKaonSpecStat2040->SetBinContent(i, histoNeutralKaonSpecStat2040->GetBinContent(i)/histoNeutralKaonSpecStat2040->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecStat2040->SetBinError(i, histoNeutralKaonSpecStat2040->GetBinError(i)/histoNeutralKaonSpecStat2040->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecSyst2040->SetBinContent(i, histoNeutralKaonSpecSyst2040->GetBinContent(i)/(histoNeutralKaonSpecSyst2040->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoNeutralKaonSpecSyst2040->GetBinContent(i) != 0){
			fractionalSystematicError = histoNeutralKaonSpecSyst2040->GetBinError(i)/histoNeutralKaonSpecSyst2040->GetBinContent(i)*100;
		}
		histoNeutralKaonSpecSyst2040->SetBinError(i, histoNeutralKaonSpecSyst2040->GetBinContent(i)*fractionalSystematicError/100.);
	}

	TH1F* histoNeutralKaonSpecStat4060 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent4060_K0s");
	histoNeutralKaonSpecStat4060->Sumw2();
	TH1F* histoNeutralKaonSpecSyst4060 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent4060_K0s");
	histoNeutralKaonSpecSyst4060->Sumw2();
	for (Int_t i = 1; i < histoNeutralKaonSpecSyst4060->GetNbinsX()+1; i++){
		histoNeutralKaonSpecStat4060->SetBinContent(i, histoNeutralKaonSpecStat4060->GetBinContent(i)/histoNeutralKaonSpecStat4060->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecStat4060->SetBinError(i, histoNeutralKaonSpecStat4060->GetBinError(i)/histoNeutralKaonSpecStat4060->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecSyst4060->SetBinContent(i, histoNeutralKaonSpecSyst4060->GetBinContent(i)/(histoNeutralKaonSpecSyst4060->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoNeutralKaonSpecSyst4060->GetBinContent(i) != 0){
			fractionalSystematicError = histoNeutralKaonSpecSyst4060->GetBinError(i)/histoNeutralKaonSpecSyst4060->GetBinContent(i)*100;
		}
		histoNeutralKaonSpecSyst4060->SetBinError(i, histoNeutralKaonSpecSyst4060->GetBinContent(i)*fractionalSystematicError/100.);
	}

	TH1F* histoNeutralKaonSpecStat6080 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent6080_K0s");
	histoNeutralKaonSpecStat6080->Sumw2();
	TH1F* histoNeutralKaonSpecSyst6080 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent6080_K0s");
	histoNeutralKaonSpecSyst6080->Sumw2();  
	for (Int_t i = 1; i < histoNeutralKaonSpecSyst6080->GetNbinsX()+1; i++){
		histoNeutralKaonSpecStat6080->SetBinContent(i, histoNeutralKaonSpecStat6080->GetBinContent(i)/histoNeutralKaonSpecStat6080->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecStat6080->SetBinError(i, histoNeutralKaonSpecStat6080->GetBinError(i)/histoNeutralKaonSpecStat6080->GetBinCenter(i)/(2*TMath::Pi()));
		histoNeutralKaonSpecSyst6080->SetBinContent(i, histoNeutralKaonSpecSyst6080->GetBinContent(i)/(histoNeutralKaonSpecSyst6080->GetBinCenter(i)*2*TMath::Pi()));
		Double_t fractionalSystematicError = 0;
		if (histoNeutralKaonSpecSyst6080->GetBinContent(i) != 0){
			fractionalSystematicError = histoNeutralKaonSpecSyst6080->GetBinError(i)/histoNeutralKaonSpecSyst6080->GetBinContent(i)*100;
		}
		histoNeutralKaonSpecSyst6080->SetBinError(i, histoNeutralKaonSpecSyst6080->GetBinContent(i)*fractionalSystematicError/100.);
	}
   
   
	
	TFile* fileChargedPionRAAFinal2013 = new TFile("ExternalInputPbPb/IdentifiedCharged/RAA_Pion_20131108.root");
	TH1D*	histoChargedPionRAAStat2040 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_20_40");
	TH1D*	histoChargedPionRAASyst2040 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_20_40");
	TH1D*	histoChargedPionRAAStat4060 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_40_60");
	TH1D*	histoChargedPionRAASyst4060 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_40_60");
	TH1D*	histoChargedPionRAAStat6080 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_60_80");
	TH1D*	histoChargedPionRAASyst6080 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_60_80");
	TH1D*	histoChargedPionRAAStat1020 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_10_20");
	TH1D*	histoChargedPionRAASyst1020 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_10_20");
	TH1D*	histoChargedPionRAAStat0005 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_0_5");
	TH1D*	histoChargedPionRAASyst0005 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_0_5");
	TH1D*	histoChargedPionRAAStat0510 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_5_10");
	TH1D*	histoChargedPionRAASyst0510 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_5_10");

   
	TFile* fileChargedKaonRAAFinal2013 = new TFile("ExternalInputPbPb/IdentifiedCharged/RAA_Kaon_20131108.root");
	TH1D* histoChargedKaonRAAStat2040 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_20_40");
	TH1D* histoChargedKaonRAASyst2040 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_20_40");
	TH1D* histoChargedKaonRAAStat4060 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_40_60");
	TH1D* histoChargedKaonRAASyst4060 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_40_60");
	TH1D* histoChargedKaonRAAStat6080 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_60_80");
	TH1D* histoChargedKaonRAASyst6080 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_60_80");
	TH1D* histoChargedKaonRAAStat1020 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_10_20");
	TH1D* histoChargedKaonRAASyst1020 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_10_20");
	TH1D* histoChargedKaonRAAStat0005 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_0_5");
	TH1D* histoChargedKaonRAASyst0005 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_0_5");
	TH1D* histoChargedKaonRAAStat0510 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_5_10");
	TH1D* histoChargedKaonRAASyst0510 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_5_10");

	
	
	cout << "bis hier" << endl;
	//******************************************************************************************************************
	//************************************ Id Charged Xiango ***********************************************************
	//******************************************************************************************************************
	TFile *fIndentifiedChargedXiango = TFile::Open("ExternalInputPbPb/IdentifiedCharged/finalspectra_Xiango_IdentifiedCharged_20130527.root");

	TList* listNegPP2760GeV = (TList*)fIndentifiedChargedXiango->Get("pp2760_NEG");
	TGraphAsymmErrors* graphXiangoChargedPionNegStatPP2760GeV = (TGraphAsymmErrors*)listNegPP2760GeV->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionNegSysPP2760GeV = (TGraphAsymmErrors*)listNegPP2760GeV->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonNegStatPP2760GeV = (TGraphAsymmErrors*)listNegPP2760GeV->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonNegSysPP2760GeV = (TGraphAsymmErrors*)listNegPP2760GeV->FindObject("grsyskaon");
	TList* listPosPP2760GeV = (TList*)fIndentifiedChargedXiango->Get("pp2760_POS");
	TGraphAsymmErrors* graphXiangoChargedPionPosStatPP2760GeV = (TGraphAsymmErrors*)listPosPP2760GeV->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionPosSysPP2760GeV = (TGraphAsymmErrors*)listPosPP2760GeV->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonPosStatPP2760GeV = (TGraphAsymmErrors*)listPosPP2760GeV->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonPosSysPP2760GeV = (TGraphAsymmErrors*)listPosPP2760GeV->FindObject("grsyskaon");   

	TList* listNegPP7TeV = (TList*)fIndentifiedChargedXiango->Get("pp7000_NEG");
	TGraphAsymmErrors* graphXiangoChargedPionNegStatPP7TeV = (TGraphAsymmErrors*)listNegPP7TeV->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionNegSysPP7TeV = (TGraphAsymmErrors*)listNegPP7TeV->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonNegStatPP7TeV = (TGraphAsymmErrors*)listNegPP7TeV->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonNegSysPP7TeV = (TGraphAsymmErrors*)listNegPP7TeV->FindObject("grsyskaon");
	TList* listPosPP7TeV = (TList*)fIndentifiedChargedXiango->Get("pp7000_POS");
	TGraphAsymmErrors* graphXiangoChargedPionPosStatPP7TeV = (TGraphAsymmErrors*)listPosPP7TeV->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionPosSysPP7TeV = (TGraphAsymmErrors*)listPosPP7TeV->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonPosStatPP7TeV = (TGraphAsymmErrors*)listPosPP7TeV->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonPosSysPP7TeV = (TGraphAsymmErrors*)listPosPP7TeV->FindObject("grsyskaon");   

	TList* listNeg0005 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_0_5_NEG");
	TGraphAsymmErrors* graphXiangoChargedPionNegStat0005 = (TGraphAsymmErrors*)listNeg0005->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionNegSys0005 = (TGraphAsymmErrors*)listNeg0005->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonNegStat0005 = (TGraphAsymmErrors*)listNeg0005->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonNegSys0005 = (TGraphAsymmErrors*)listNeg0005->FindObject("grsyskaon");
	TList* listPos0005 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_0_5_POS");
	TGraphAsymmErrors* graphXiangoChargedPionPosStat0005 = (TGraphAsymmErrors*)listPos0005->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionPosSys0005 = (TGraphAsymmErrors*)listPos0005->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonPosStat0005 = (TGraphAsymmErrors*)listPos0005->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonPosSys0005 = (TGraphAsymmErrors*)listPos0005->FindObject("grsyskaon");   
	
	TList* listNeg0510 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_5_10_NEG");
	TGraphAsymmErrors* graphXiangoChargedPionNegStat0510 = (TGraphAsymmErrors*)listNeg0510->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionNegSys0510 = (TGraphAsymmErrors*)listNeg0510->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonNegStat0510 = (TGraphAsymmErrors*)listNeg0510->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonNegSys0510 = (TGraphAsymmErrors*)listNeg0510->FindObject("grsyskaon");
	TList* listPos0510 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_5_10_POS");
	TGraphAsymmErrors* graphXiangoChargedPionPosStat0510 = (TGraphAsymmErrors*)listPos0510->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionPosSys0510 = (TGraphAsymmErrors*)listPos0510->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonPosStat0510 = (TGraphAsymmErrors*)listPos0510->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonPosSys0510 = (TGraphAsymmErrors*)listPos0510->FindObject("grsyskaon");   

	TList* listNeg1020 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_10_20_NEG");
	TGraphAsymmErrors* graphXiangoChargedPionNegStat1020 = (TGraphAsymmErrors*)listNeg1020->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionNegSys1020 = (TGraphAsymmErrors*)listNeg1020->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonNegStat1020 = (TGraphAsymmErrors*)listNeg1020->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonNegSys1020 = (TGraphAsymmErrors*)listNeg1020->FindObject("grsyskaon");
	TList* listPos1020 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_10_20_POS");
	TGraphAsymmErrors* graphXiangoChargedPionPosStat1020 = (TGraphAsymmErrors*)listPos1020->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionPosSys1020 = (TGraphAsymmErrors*)listPos1020->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonPosStat1020 = (TGraphAsymmErrors*)listPos1020->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonPosSys1020 = (TGraphAsymmErrors*)listPos1020->FindObject("grsyskaon");   

	TList* listNeg2040 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_20_40_NEG");
	TGraphAsymmErrors* graphXiangoChargedPionNegStat2040 = (TGraphAsymmErrors*)listNeg2040->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionNegSys2040 = (TGraphAsymmErrors*)listNeg2040->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonNegStat2040 = (TGraphAsymmErrors*)listNeg2040->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonNegSys2040 = (TGraphAsymmErrors*)listNeg2040->FindObject("grsyskaon");
	TList* listPos2040 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_20_40_POS");
	TGraphAsymmErrors* graphXiangoChargedPionPosStat2040 = (TGraphAsymmErrors*)listPos2040->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionPosSys2040 = (TGraphAsymmErrors*)listPos2040->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonPosStat2040 = (TGraphAsymmErrors*)listPos2040->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonPosSys2040 = (TGraphAsymmErrors*)listPos2040->FindObject("grsyskaon");   

	TList* listNeg4060 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_40_60_NEG");
	TGraphAsymmErrors* graphXiangoChargedPionNegStat4060 = (TGraphAsymmErrors*)listNeg4060->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionNegSys4060 = (TGraphAsymmErrors*)listNeg4060->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonNegStat4060 = (TGraphAsymmErrors*)listNeg4060->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonNegSys4060 = (TGraphAsymmErrors*)listNeg4060->FindObject("grsyskaon");
	TList* listPos4060 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_40_60_POS");
	TGraphAsymmErrors* graphXiangoChargedPionPosStat4060 = (TGraphAsymmErrors*)listPos4060->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionPosSys4060 = (TGraphAsymmErrors*)listPos4060->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonPosStat4060 = (TGraphAsymmErrors*)listPos4060->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonPosSys4060 = (TGraphAsymmErrors*)listPos4060->FindObject("grsyskaon");   
	
	TList* listNeg6080 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_60_80_NEG");
	TGraphAsymmErrors* graphXiangoChargedPionNegStat6080 = (TGraphAsymmErrors*)listNeg6080->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionNegSys6080 = (TGraphAsymmErrors*)listNeg6080->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonNegStat6080 = (TGraphAsymmErrors*)listNeg6080->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonNegSys6080 = (TGraphAsymmErrors*)listNeg6080->FindObject("grsyskaon");
	TList* listPos6080 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_60_80_POS");
	TGraphAsymmErrors* graphXiangoChargedPionPosStat6080 = (TGraphAsymmErrors*)listPos6080->FindObject("grstatpion");
	TGraphAsymmErrors* graphXiangoChargedPionPosSys6080 = (TGraphAsymmErrors*)listPos6080->FindObject("grsyspion");
	TGraphAsymmErrors* graphXiangoChargedKaonPosStat6080 = (TGraphAsymmErrors*)listPos6080->FindObject("grstatkaon");
	TGraphAsymmErrors* graphXiangoChargedKaonPosSys6080 = (TGraphAsymmErrors*)listPos6080->FindObject("grsyskaon");   

	
	// ************************** Low Pt  charged pions pPb************************************************
	TFile *fPionLowPtpPbMinBias = TFile::Open("ExternalInputpPb/20130201.CombinedSpectra_pA_pilot_ycms_0005_minbias.root");
	
	TH1D* histoChargedPionPlusSpecLowPtStatpPb = (TH1D*)fPionLowPtpPbMinBias->Get("stat_mb_pion_plus");
	TH1D* histoChargedPionPlusSpecLowPtSyspPb = (TH1D*)fPionLowPtpPbMinBias->Get("sys_mb_pion_plus");
	TH1D* histoChargedPionMinusSpecLowPtStatpPb = (TH1D*)fPionLowPtpPbMinBias->Get("stat_mb_pion_minus");
	TH1D* histoChargedPionMinusSpecLowPtSyspPb = (TH1D*)fPionLowPtpPbMinBias->Get("sys_mb_pion_minus");
	
	TH1D* histoChargedPionSpecLowPtStatpPb = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb->Clone("histoChargedPionSpecLowPtStatpPb");
	histoChargedPionSpecLowPtStatpPb->Add(histoChargedPionPlusSpecLowPtStatpPb);
	histoChargedPionSpecLowPtStatpPb->Scale(0.5);

	TH1D* histoChargedPionSpecLowPtSyspPb = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb->Clone("histoChargedPionSpecLowPtSyspPb");   
	
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStatpPb->SetBinContent(i, histoChargedPionSpecLowPtStatpPb->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStatpPb->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb->GetBinError(i)/histoChargedPionSpecLowPtStatpPb->GetBinCenter(i)/(2*TMath::Pi()));
		Double_t error = 0;
		if (histoChargedPionMinusSpecLowPtSyspPb->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb->GetBinContent(i) != 0){
			error= (histoChargedPionMinusSpecLowPtSyspPb->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb->GetBinContent(i);
		} else {
			error = 0;
		}
		histoChargedPionSpecLowPtSyspPb->SetBinContent(i, histoChargedPionSpecLowPtStatpPb->GetBinContent(i));
		histoChargedPionSpecLowPtSyspPb->SetBinError(i,error);      
	}
	cout << "bis hier" << endl;
	// ************************** Low Pt  charged pions pPb depending on centrality************************************************
	TFile *fPionLowPtpPbCent = TFile::Open("ExternalInputpPb/20130201.CombinedSpectra_pA_pilot_ycms_0005.root");
	
	TH1D* histoChargedPionPlusSpecLowPtStatpPb0020 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent0_pion_plus");
	TH1D* histoChargedPionPlusSpecLowPtSyspPb0020 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent0_pion_plus");
	TH1D* histoChargedPionMinusSpecLowPtStatpPb0020 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent0_pion_minus");
	TH1D* histoChargedPionMinusSpecLowPtSyspPb0020 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent0_pion_minus");
	
	TH1D* histoChargedPionSpecLowPtStatpPb0020 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb0020->Clone("histoChargedPionSpecLowPtStatpPb0020");
	histoChargedPionSpecLowPtStatpPb0020->Add(histoChargedPionPlusSpecLowPtStatpPb0020);
	histoChargedPionSpecLowPtStatpPb0020->Scale(0.5);

	TH1D* histoChargedPionSpecLowPtSyspPb0020 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb0020->Clone("histoChargedPionSpecLowPtSyspPb0020");   
	
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb0020->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStatpPb0020->SetBinContent(i, histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb0020->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStatpPb0020->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb0020->GetBinError(i)/histoChargedPionSpecLowPtStatpPb0020->GetBinCenter(i)/(2*TMath::Pi()));
		Double_t error = 0;
		if (histoChargedPionMinusSpecLowPtSyspPb0020->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb0020->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i) != 0){
			error= (histoChargedPionMinusSpecLowPtSyspPb0020->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb0020->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb0020->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb0020->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i);
		} else {
			error = 0;
		}
		histoChargedPionSpecLowPtSyspPb0020->SetBinContent(i, histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i));
		histoChargedPionSpecLowPtSyspPb0020->SetBinError(i,error);      
	}
	cout << "bis hier" << endl;
	TH1D* histoChargedPionPlusSpecLowPtStatpPb2040 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent1_pion_plus");
	TH1D* histoChargedPionPlusSpecLowPtSyspPb2040 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent1_pion_plus");
	TH1D* histoChargedPionMinusSpecLowPtStatpPb2040 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent1_pion_minus");
	TH1D* histoChargedPionMinusSpecLowPtSyspPb2040 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent1_pion_minus");
	
	TH1D* histoChargedPionSpecLowPtStatpPb2040 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb2040->Clone("histoChargedPionSpecLowPtStatpPb2040");
	histoChargedPionSpecLowPtStatpPb2040->Add(histoChargedPionPlusSpecLowPtStatpPb2040);
	histoChargedPionSpecLowPtStatpPb2040->Scale(0.5);

	TH1D* histoChargedPionSpecLowPtSyspPb2040 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb2040->Clone("histoChargedPionSpecLowPtSyspPb2040");   
	
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb2040->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStatpPb2040->SetBinContent(i, histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb2040->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStatpPb2040->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb2040->GetBinError(i)/histoChargedPionSpecLowPtStatpPb2040->GetBinCenter(i)/(2*TMath::Pi()));
		Double_t error = 0;
		if (histoChargedPionMinusSpecLowPtSyspPb2040->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb2040->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i) != 0){
			error= (histoChargedPionMinusSpecLowPtSyspPb2040->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb2040->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb2040->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb2040->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i);
		} else {
			error = 0;
		}
		histoChargedPionSpecLowPtSyspPb2040->SetBinContent(i, histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i));
		histoChargedPionSpecLowPtSyspPb2040->SetBinError(i,error);      
	}
	TH1D* histoChargedPionPlusSpecLowPtStatpPb4060 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent2_pion_plus");
	TH1D* histoChargedPionPlusSpecLowPtSyspPb4060 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent2_pion_plus");
	TH1D* histoChargedPionMinusSpecLowPtStatpPb4060 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent2_pion_minus");
	TH1D* histoChargedPionMinusSpecLowPtSyspPb4060 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent2_pion_minus");
	
	TH1D* histoChargedPionSpecLowPtStatpPb4060 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb4060->Clone("histoChargedPionSpecLowPtStatpPb4060");
	histoChargedPionSpecLowPtStatpPb4060->Add(histoChargedPionPlusSpecLowPtStatpPb4060);
	histoChargedPionSpecLowPtStatpPb4060->Scale(0.5);

	TH1D* histoChargedPionSpecLowPtSyspPb4060 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb4060->Clone("histoChargedPionSpecLowPtSyspPb4060");   
	
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb4060->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStatpPb4060->SetBinContent(i, histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb4060->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStatpPb4060->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb4060->GetBinError(i)/histoChargedPionSpecLowPtStatpPb4060->GetBinCenter(i)/(2*TMath::Pi()));
		Double_t error = 0;
		if (histoChargedPionMinusSpecLowPtSyspPb4060->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb4060->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i) != 0){
			error= (histoChargedPionMinusSpecLowPtSyspPb4060->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb4060->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb4060->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb4060->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i);
		} else {
			error = 0;
		}
		histoChargedPionSpecLowPtSyspPb4060->SetBinContent(i, histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i));
		histoChargedPionSpecLowPtSyspPb4060->SetBinError(i,error);      
	}
	cout << "bis hier" << endl;
	TH1D* histoChargedPionPlusSpecLowPtStatpPb6080 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent3_pion_plus");
	TH1D* histoChargedPionPlusSpecLowPtSyspPb6080 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent3_pion_plus");
	TH1D* histoChargedPionMinusSpecLowPtStatpPb6080 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent3_pion_minus");
	TH1D* histoChargedPionMinusSpecLowPtSyspPb6080 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent3_pion_minus");
	
	TH1D* histoChargedPionSpecLowPtStatpPb6080 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb6080->Clone("histoChargedPionSpecLowPtStatpPb6080");
	histoChargedPionSpecLowPtStatpPb6080->Add(histoChargedPionPlusSpecLowPtStatpPb6080);
	histoChargedPionSpecLowPtStatpPb6080->Scale(0.5);

	TH1D* histoChargedPionSpecLowPtSyspPb6080 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb6080->Clone("histoChargedPionSpecLowPtSyspPb6080");   
	
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb6080->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStatpPb6080->SetBinContent(i, histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb6080->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStatpPb6080->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb6080->GetBinError(i)/histoChargedPionSpecLowPtStatpPb6080->GetBinCenter(i)/(2*TMath::Pi()));
		Double_t error = 0;
		if (histoChargedPionMinusSpecLowPtSyspPb6080->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb6080->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i) != 0){
			error= (histoChargedPionMinusSpecLowPtSyspPb6080->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb6080->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb6080->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb6080->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i);
		} else {
			error = 0;
		}
		histoChargedPionSpecLowPtSyspPb6080->SetBinContent(i, histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i));
		histoChargedPionSpecLowPtSyspPb6080->SetBinError(i,error);      
	}
	TH1D* histoChargedPionPlusSpecLowPtStatpPb80100 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent4_pion_plus");
	TH1D* histoChargedPionPlusSpecLowPtSyspPb80100 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent4_pion_plus");
	TH1D* histoChargedPionMinusSpecLowPtStatpPb80100 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent4_pion_minus");
	TH1D* histoChargedPionMinusSpecLowPtSyspPb80100 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent4_pion_minus");
	
	TH1D* histoChargedPionSpecLowPtStatpPb80100 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb80100->Clone("histoChargedPionSpecLowPtStatpPb80100");
	histoChargedPionSpecLowPtStatpPb80100->Add(histoChargedPionPlusSpecLowPtStatpPb80100);
	histoChargedPionSpecLowPtStatpPb80100->Scale(0.5);

	TH1D* histoChargedPionSpecLowPtSyspPb80100 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb80100->Clone("histoChargedPionSpecLowPtSyspPb80100");   
	
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb80100->GetNbinsX()+1; i++){
		histoChargedPionSpecLowPtStatpPb80100->SetBinContent(i, histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb80100->GetBinCenter(i)/(2*TMath::Pi()));
		histoChargedPionSpecLowPtStatpPb80100->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb80100->GetBinError(i)/histoChargedPionSpecLowPtStatpPb80100->GetBinCenter(i)/(2*TMath::Pi()));
		Double_t error = 0;
		if (histoChargedPionMinusSpecLowPtSyspPb80100->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb80100->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i) != 0){
			error= (histoChargedPionMinusSpecLowPtSyspPb80100->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb80100->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb80100->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb80100->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i);
		} else {
			error = 0;
		}
		histoChargedPionSpecLowPtSyspPb80100->SetBinContent(i, histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i));
		histoChargedPionSpecLowPtSyspPb80100->SetBinError(i,error);      
	}
   
      //------------------------------------------- pPb Full pT -----------------------------------

	TFile* fileChargedPionsFullpT = new TFile("ExternalInputpPb/FullPreliminaryPionSpectra.root");
	TH1D* histoChargedPionSys0_5= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_0_5"); //0-5% 
	TH1D* histoChargedPionStat0_5= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_0_5");
	TH1D* histoChargedPionSys5_10= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_5_10"); //0-5% 
	TH1D* histoChargedPionStat5_10= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_5_10");
	TH1D* histoChargedPionSys10_20= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_10_20"); //0-5% 
	TH1D* histoChargedPionStat10_20= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_10_20");
	TH1D* histoChargedPionSys20_40= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_20_40"); //0-5% 
	TH1D* histoChargedPionStat20_40= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_20_40");
	TH1D* histoChargedPionSys40_60= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_40_60"); //0-5% 
	TH1D* histoChargedPionStat40_60= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_40_60");
	TH1D* histoChargedPionSys60_80= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_60_80"); //0-5% 
	TH1D* histoChargedPionStat60_80= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_60_80");
	TH1D* histoChargedPionSys80_100= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Syst_80_100"); //0-5% 
	TH1D* histoChargedPionStat80_100= (TH1D*)fileChargedPionsFullpT->Get("hPionSpectrum_Stat_80_100");

	histoChargedPionSys0_5->Scale(0.5);
	histoChargedPionStat0_5->Scale(0.5);
	histoChargedPionSys5_10->Scale(0.5);
	histoChargedPionStat5_10->Scale(0.5);
	histoChargedPionSys10_20->Scale(0.5);
	histoChargedPionStat10_20->Scale(0.5);
	histoChargedPionSys20_40->Scale(0.5);
	histoChargedPionStat20_40->Scale(0.5);
	histoChargedPionSys40_60->Scale(0.5);
	histoChargedPionStat40_60->Scale(0.5);
	histoChargedPionSys60_80->Scale(0.5);
	histoChargedPionStat60_80->Scale(0.5);
	histoChargedPionSys80_100->Scale(0.5);
	histoChargedPionStat80_100->Scale(0.5);


	TH1D *histopPbChargedPionSpecFullpTSys=(TH1D*)histoChargedPionSys0_5->Clone();
	TH1D *histopPbChargedPionSpecFullpTStat=(TH1D*)histoChargedPionStat0_5->Clone();

	histopPbChargedPionSpecFullpTSys->Add(histoChargedPionSys0_5,histoChargedPionSys5_10,0.05,0.05); //weight with the centrality bin width
	histopPbChargedPionSpecFullpTSys->Add(histoChargedPionSys10_20,0.1);
	histopPbChargedPionSpecFullpTSys->Add(histoChargedPionSys20_40,0.2);
	histopPbChargedPionSpecFullpTSys->Add(histoChargedPionSys40_60,0.2);
	histopPbChargedPionSpecFullpTSys->Add(histoChargedPionSys60_80,0.2);
	histopPbChargedPionSpecFullpTSys->Add(histoChargedPionSys80_100,0.2);

	histopPbChargedPionSpecFullpTStat->Add(histoChargedPionStat0_5,histoChargedPionStat5_10,0.05,0.05); //weight with the centrality bin width
	histopPbChargedPionSpecFullpTStat->Add(histoChargedPionStat10_20,0.1);
	histopPbChargedPionSpecFullpTStat->Add(histoChargedPionStat20_40,0.2);
	histopPbChargedPionSpecFullpTStat->Add(histoChargedPionStat40_60,0.2);
	histopPbChargedPionSpecFullpTStat->Add(histoChargedPionStat60_80,0.2);
	histopPbChargedPionSpecFullpTStat->Add(histoChargedPionStat80_100,0.2);

	TH1D *histopPbChargedPionSpecFullpTSys0020=(TH1D*)histoChargedPionSys0_5->Clone();
	TH1D *histopPbChargedPionSpecFullpTStat0020=(TH1D*)histoChargedPionStat0_5->Clone();
	histopPbChargedPionSpecFullpTSys0020->Add(histoChargedPionSys0_5,histoChargedPionSys5_10,0.25,0.25); //weight with the centrality bin width
	histopPbChargedPionSpecFullpTSys0020->Add(histoChargedPionSys10_20,0.5);
	histopPbChargedPionSpecFullpTStat0020->Add(histoChargedPionStat0_5,histoChargedPionStat5_10,0.25,0.25); //weight with the centrality bin width
	histopPbChargedPionSpecFullpTStat0020->Add(histoChargedPionStat10_20,0.5);

	TH1D *histopPbChargedPionSpecFullpTSys2040=(TH1D*)histoChargedPionSys20_40->Clone();
	TH1D *histopPbChargedPionSpecFullpTStat2040=(TH1D*)histoChargedPionStat20_40->Clone();

	TH1D *histopPbChargedPionSpecFullpTSys4060=(TH1D*)histoChargedPionSys40_60->Clone();
	TH1D *histopPbChargedPionSpecFullpTStat4060=(TH1D*)histoChargedPionStat40_60->Clone();

	TH1D *histopPbChargedPionSpecFullpTSys6080=(TH1D*)histoChargedPionSys60_80->Clone();
	TH1D *histopPbChargedPionSpecFullpTStat6080=(TH1D*)histoChargedPionStat60_80->Clone();
	TH1D *histopPbChargedPionSpecFullpTSys60100=(TH1D*)histoChargedPionSys60_80->Clone();
	TH1D *histopPbChargedPionSpecFullpTStat60100=(TH1D*)histoChargedPionStat60_80->Clone();
	histopPbChargedPionSpecFullpTSys60100->Add(histoChargedPionSys60_80,histoChargedPionSys80_100,0.5,0.5);
	histopPbChargedPionSpecFullpTStat60100->Add(histoChargedPionStat60_80,histoChargedPionStat80_100,0.5,0.5);


	Double_t SysErrorRel[7]={0};
	Double_t RelErrorMB=0;
	Double_t RelError0020=0;
	Double_t RelError60100=0;
	for(Int_t i = 1; i < histoChargedPionSys0_5->GetNbinsX()+1; i++){
		if (histoChargedPionSys0_5->GetBinContent(i) != 0){
			SysErrorRel[0]= histoChargedPionSys0_5->GetBinError(i)/histoChargedPionSys0_5->GetBinContent(i)*100 ;
		} else    SysErrorRel[0] = 0;
		if (histoChargedPionSys5_10->GetBinContent(i) != 0){
			SysErrorRel[1]= histoChargedPionSys5_10->GetBinError(i)/histoChargedPionSys5_10->GetBinContent(i)*100 ;
		} else    SysErrorRel[1] = 0;
		if (histoChargedPionSys10_20->GetBinContent(i) != 0){
			SysErrorRel[2]= histoChargedPionSys10_20->GetBinError(i)/histoChargedPionSys10_20->GetBinContent(i)*100 ;
		} else    SysErrorRel[2] = 0;
		if (histoChargedPionSys20_40->GetBinContent(i) != 0){
			SysErrorRel[3]= histoChargedPionSys20_40->GetBinError(i)/histoChargedPionSys20_40->GetBinContent(i)*100 ;
		} else    SysErrorRel[3] = 0;
		if (histoChargedPionSys40_60->GetBinContent(i) != 0){
			SysErrorRel[4]= histoChargedPionSys40_60->GetBinError(i)/histoChargedPionSys40_60->GetBinContent(i)*100 ;
		} else    SysErrorRel[4] = 0;
		if (histoChargedPionSys60_80->GetBinContent(i) != 0){
			SysErrorRel[5]= histoChargedPionSys60_80->GetBinError(i)/histoChargedPionSys60_80->GetBinContent(i)*100 ;
		} else    SysErrorRel[5] = 0;
		if (histoChargedPionSys80_100->GetBinContent(i) != 0){
			SysErrorRel[6]= histoChargedPionSys80_100->GetBinError(i)/histoChargedPionSys80_100->GetBinContent(i)*100 ;
		} else    SysErrorRel[6] = 0;
		RelErrorMB=SysErrorRel[0];
		for(Int_t j= 1; j < 7; j++){
		
		if(SysErrorRel[j-1]<SysErrorRel[j])      RelErrorMB=SysErrorRel[j];
		}
		RelError0020=SysErrorRel[0];
		for(Int_t j= 1; j < 3; j++){
			if(SysErrorRel[j-1]<SysErrorRel[j])      RelError0020=SysErrorRel[j];
		}
		RelError60100=SysErrorRel[5];
		if(SysErrorRel[5]<SysErrorRel[6])      RelError0020=SysErrorRel[6];

		histopPbChargedPionSpecFullpTSys->SetBinError(i, histopPbChargedPionSpecFullpTSys->GetBinContent(i)*RelErrorMB/100);
		histopPbChargedPionSpecFullpTSys0020->SetBinError(i, histopPbChargedPionSpecFullpTSys0020->GetBinContent(i)*RelError0020/100);
		histopPbChargedPionSpecFullpTSys60100->SetBinError(i, histopPbChargedPionSpecFullpTSys60100->GetBinContent(i)*RelError60100/100);
	}

	  gROOT->SetStyle("Plain");

// 	// Plot: p8462_d1x1y3
// 	Double_t p8462_d1x1y1_xval[] = { 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 
// 		0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.05, 1.15, 
// 		1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 
// 		2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.25, 4.75, 
// 		5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 
// 		13.5, 14.5, 15.5, 17.0, 19.0 };
// 	Double_t p8462_d1x1y1_xerrminus[] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
// 								 		  0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 
// 										  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 
// 										  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 
//  									  0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
// 									  	  0.5, 0.5, 0.5, 1.0, 1.0 };
// 	Double_t p8462_d1x1y1_xerrplus[] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
// 								 		  0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 
// 										  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 
// 										  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 
//  									  0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
// 									  	  0.5, 0.5, 0.5, 1.0, 1.0 };
// 	Double_t p8462_d1x1y1_yval[] = { 256.2, 198.9, 154.9, 119.4, 91.68, 70.9, 54.94, 42.71, 33.56, 
// 							 		 26.51, 21.06, 16.9, 13.54, 11.0, 9.002, 7.431, 6.098, 4.615, 3.201, 
// 		 							 2.27, 1.625, 1.162, 0.8488, 0.6375, 0.4721, 0.3581, 0.2728, 0.1874, 0.114, 
// 									 0.07148, 0.04502, 0.02896, 0.01986, 0.01325, 0.008893, 0.006293, 0.004496, 0.0026, 0.001199, 
// 								  	 6.158E-4, 3.394E-4, 1.857E-4, 1.025E-4, 5.206E-5, 2.229E-5, 9.721E-6, 5.3E-6, 3.236E-6, 1.825E-6, 
// 									 5.102E-7, 6.276E-7, 4.19E-7, 1.066E-7, 8.192E-8 };
// 	Double_t p8462_d1x1y1_ysysminus[] = { 20.9, 13.2, 10.3, 7.9, 6.10, 4.72, 3.65, 2.84, 2.23, 
// 										  1.76, 1.40, 1.12, 0.90, 0.73, 0.599, 0.494, 0.406, 0.307, 0.213, 
// 										  0.151, 0.108, 0.077, 0.0565, 0.0424, 0.0314, 0.0238, 0.0181, 0.0125, 0.0076, 
// 										  0.00475, 0.00299, 0.00193, 0.00132, 8.8E-4, 5.91E-4, 4.19E-4, 2.99E-4, 1.73E-4, 8.0E-5, 
// 										  4.1E-5, 2.26E-5, 1.24E-5, 6.8E-6, 3.46E-6, 1.48E-6, 0.647E-6, 3.59E-7, 2.19E-7, 1.24E-7, 
// 								 		  0.346E-7, 0.427E-7, 0.286E-7, 0.73E-8, 0.568E-8 };
// 	Double_t p8462_d1x1y1_ysysplus[] = { 20.9, 13.2, 10.3, 7.9, 6.10, 4.72, 3.65, 2.84, 2.23, 
// 										  1.76, 1.40, 1.12, 0.90, 0.73, 0.599, 0.494, 0.406, 0.307, 0.213, 
// 										  0.151, 0.108, 0.077, 0.0565, 0.0424, 0.0314, 0.0238, 0.0181, 0.0125, 0.0076, 
// 										  0.00475, 0.00299, 0.00193, 0.00132, 8.8E-4, 5.91E-4, 4.19E-4, 2.99E-4, 1.73E-4, 8.0E-5, 
// 										  4.1E-5, 2.26E-5, 1.24E-5, 6.8E-6, 3.46E-6, 1.48E-6, 0.647E-6, 3.59E-7, 2.19E-7, 1.24E-7, 
// 								 		  0.346E-7, 0.427E-7, 0.286E-7, 0.73E-8, 0.568E-8 };
// 	Double_t p8462_d1x1y1_ystatminus[] = { 0.3, 0.2, 0.1, 0.1, 0.09, 0.07, 0.06, 0.05, 0.04, 
// 		0.04, 0.03, 0.03, 0.02, 0.02, 0.018, 0.016, 0.014, 0.008, 0.006, 
// 		0.005, 0.004, 0.003, 0.0028, 0.0024, 0.002, 0.0017, 0.0014, 8.0E-4, 6.0E-4, 
// 		4.6E-4, 3.4E-4, 2.8E-4, 2.1E-4, 1.6E-4, 1.25E-4, 1.02E-4, 8.3E-5, 4.0E-5, 2.5E-5, 
// 		1.68E-5, 1.18E-5, 8.4E-6, 6.0E-6, 2.84E-6, 1.75E-6, 1.093E-6, 7.6E-7, 5.74E-7, 4.09E-7, 
// 		2.09E-7, 2.228E-7, 1.716E-7, 6.16E-8, 4.732E-8 };
// 	Double_t p8462_d1x1y1_ystatplus[] = { 0.3, 0.2, 0.1, 0.1, 0.09, 0.07, 0.06, 0.05, 0.04, 
// 		0.04, 0.03, 0.03, 0.02, 0.02, 0.018, 0.016, 0.014, 0.008, 0.006, 
// 		0.005, 0.004, 0.003, 0.0028, 0.0024, 0.002, 0.0017, 0.0014, 8.0E-4, 6.0E-4, 
// 		4.6E-4, 3.4E-4, 2.8E-4, 2.1E-4, 1.6E-4, 1.25E-4, 1.02E-4, 8.3E-5, 4.0E-5, 2.5E-5, 
// 		1.68E-5, 1.18E-5, 8.4E-6, 6.0E-6, 2.84E-6, 1.75E-6, 1.093E-6, 7.6E-7, 5.74E-7, 4.09E-7, 
// 		2.09E-7, 2.228E-7, 1.716E-7, 6.16E-8, 4.732E-8 };
// 	int p8462_d1x1y1_numpoints = 54;
// 	TGraphAsymmErrors* graphChargedHadronsSysPP900GeV = new TGraphAsymmErrors(p8462_d1x1y1_numpoints, p8462_d1x1y1_xval, p8462_d1x1y1_yval, p8462_d1x1y1_xerrminus, p8462_d1x1y1_xerrplus, 
// 																			   p8462_d1x1y1_ysysminus, p8462_d1x1y1_ysysplus);
// 	TGraphAsymmErrors* graphChargedHadronsStatPP900GeV = new TGraphAsymmErrors(p8462_d1x1y1_numpoints, p8462_d1x1y1_xval, p8462_d1x1y1_yval, p8462_d1x1y1_xerrminus, p8462_d1x1y1_xerrplus, 
// 																				p8462_d1x1y1_ystatminus, p8462_d1x1y1_ystatplus);

	Double_t p8462_d1x1y2_xval[65] = { 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 
		0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.05, 1.15, 
		1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 
		2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.25, 4.75, 
		5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 
		13.5, 14.5, 15.5, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 
		31.0 };
	Double_t p8462_d1x1y2_xerrminus[65] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
										0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 
										0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.10, 0.1, 
										0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 
										0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
										0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
										1.0 };
	Double_t p8462_d1x1y2_xerrplus[65] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
										0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 
										0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.10, 0.1, 
										0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 
										0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
										0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
										1.0 };
	Double_t p8462_d1x1y2_yval[65] = { 357.6, 286.0, 223.9, 172.8, 133.2, 103.3, 80.54, 63.46, 50.38, 
		40.29, 32.43, 26.33, 21.57, 17.65, 14.62, 12.15, 10.13, 7.799, 5.571, 
		4.066, 2.984, 2.224, 1.662, 1.266, 0.9676, 0.7463, 0.5803, 0.4131, 0.2614, 
		0.1667, 0.1104, 0.0749, 0.05148, 0.03607, 0.02516, 0.01807, 0.01326, 0.007935, 0.003978, 
		0.002122, 0.001209, 7.046E-4, 4.346E-4, 2.269E-4, 1.008E-4, 4.915E-5, 2.735E-5, 1.488E-5, 8.995E-6, 
		5.704E-6, 3.021E-6, 2.003E-6, 1.252E-6, 5.819E-7, 3.291E-7, 1.359E-7, 1.163E-7, 5.505E-8, 3.128E-8, 
		2.73E-8 };
	Double_t p8462_d1x1y2_ysysminus[65] = { 27.2, 18.1, 15.1, 12.40, 10.10, 8.3, 6.39, 4.92, 3.81, 
										2.89, 2.33, 1.89, 1.55, 1.27, 1.05, 0.87, 0.73, 0.56, 0.400, 
										0.292, 0.214, 0.160, 0.119, 0.091, 0.0695, 0.0536, 0.0417, 0.0297, 0.0188, 
										0.0120, 0.0079, 0.00538, 0.00370, 0.00259, 0.00181, 0.00130, 9.5E-4, 5.7E-4, 2.86E-4, 
										1.52E-4, 8.7E-5, 5.06E-5, 3.12E-5, 1.63E-5, 7.2E-6, 3.53E-6, 2.07E-6, 1.13E-6, 6.84E-7, 
										4.36E-7, 2.31E-7, 1.54E-7, 0.96E-7, 4.51E-8, 2.55E-8, 1.06E-8, 0.91E-8, 0.431E-8, 0.246E-8, 
										2.16E-9 };
	Double_t p8462_d1x1y2_ysysplus[65] = { 27.2, 18.1, 15.1, 12.40, 10.10, 8.3, 6.39, 4.92, 3.81, 
										2.89, 2.33, 1.89, 1.55, 1.27, 1.05, 0.87, 0.73, 0.56, 0.400, 
										0.292, 0.214, 0.160, 0.119, 0.091, 0.0695, 0.0536, 0.0417, 0.0297, 0.0188, 
										0.0120, 0.0079, 0.00538, 0.00370, 0.00259, 0.00181, 0.00130, 9.5E-4, 5.7E-4, 2.86E-4, 
										1.52E-4, 8.7E-5, 5.06E-5, 3.12E-5, 1.63E-5, 7.2E-6, 3.53E-6, 2.07E-6, 1.13E-6, 6.84E-7, 
										4.36E-7, 2.31E-7, 1.54E-7, 0.96E-7, 4.51E-8, 2.55E-8, 1.06E-8, 0.91E-8, 0.431E-8, 0.246E-8, 
										2.16E-9 };
	Double_t p8462_d1x1y2_ystatminus[65] = { 0.5, 0.2, 0.2, 0.1, 0.1, 0.1, 0.07, 0.06, 0.05, 
		0.04, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.008, 0.006, 
		0.005, 0.004, 0.004, 0.003, 0.003, 0.002, 0.0019, 0.0015, 0.001, 6.0E-4, 
		5.0E-4, 4.0E-4, 4.2E-4, 2.0E-4, 1.5E-4, 1.1E-4, 8.0E-5, 7.0E-5, 4.3E-5, 2.4E-5, 
		1.5E-5, 1.0E-5, 6.6E-6, 4.8E-6, 2.3E-6, 1.4E-6, 9.0E-7, 6.3E-7, 4.4E-7, 3.25E-7, 
		2.48E-7, 1.74E-7, 1.36E-7, 7.3E-8, 4.69E-8, 3.35E-8, 2.05E-8, 1.82E-8, 1.203E-8, 8.68E-9, 
		7.89E-9 };
	Double_t p8462_d1x1y2_ystatplus[65] = { 0.5, 0.2, 0.2, 0.1, 0.1, 0.1, 0.07, 0.06, 0.05, 
		0.04, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.008, 0.006, 
		0.005, 0.004, 0.004, 0.003, 0.003, 0.002, 0.0019, 0.0015, 0.001, 6.0E-4, 
		5.0E-4, 4.0E-4, 4.2E-4, 2.0E-4, 1.5E-4, 1.1E-4, 8.0E-5, 7.0E-5, 4.3E-5, 2.4E-5, 
		1.5E-5, 1.0E-5, 6.6E-6, 4.8E-6, 2.3E-6, 1.4E-6, 9.0E-7, 6.3E-7, 4.4E-7, 3.25E-7, 
		2.48E-7, 1.74E-7, 1.36E-7, 7.3E-8, 4.69E-8, 3.35E-8, 2.05E-8, 1.82E-8, 1.203E-8, 8.68E-9, 
		7.89E-9 };
	int p8462_d1x1y2_numpoints = 60;
	TGraphAsymmErrors* graphChargedHadronsSysPP2760GeV = new TGraphAsymmErrors(p8462_d1x1y2_numpoints, p8462_d1x1y2_xval, p8462_d1x1y2_yval, p8462_d1x1y2_xerrminus, p8462_d1x1y2_xerrplus,
																			   p8462_d1x1y2_ysysminus, p8462_d1x1y2_ysysplus);
	TGraphAsymmErrors* graphChargedHadronsStatPP2760GeV = new TGraphAsymmErrors(p8462_d1x1y2_numpoints, p8462_d1x1y2_xval, p8462_d1x1y2_yval, p8462_d1x1y2_xerrminus, p8462_d1x1y2_xerrplus,
																				p8462_d1x1y2_ystatminus, p8462_d1x1y2_ystatplus);
	graphChargedHadronsStatPP2760GeV->Print();
// 	Double_t p8462_d1x1y3_xval[] = { 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 
// 		0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.05, 1.15, 
// 		1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 
// 		2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.25, 4.75, 
// 		5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 
// 		13.5, 14.5, 15.5, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 
// 		31.0, 33.0, 35.0, 38.0, 42.5, 47.5 };
// 	Double_t p8462_d1x1y3_xerrminus[] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
// 										  0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 
// 										  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.10, 0.1, 
// 										  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 
// 								 		  0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
// 										  0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
// 										  1.0, 1.0, 1.0, 2.0, 2.5, 2.5 };
// 	Double_t p8462_d1x1y3_xerrplus[] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
// 										  0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 
// 										  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.10, 0.1, 
// 										  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 
// 								 		  0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
// 										  0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
// 										  1.0, 1.0, 1.0, 2.0, 2.5, 2.5 };
// 	Double_t p8462_d1x1y3_yval[] = { 490.7, 389.4, 301.1, 232.3, 179.8, 140.3, 110.3, 87.47, 70.0, 
// 									56.52, 46.04, 37.71, 31.06, 25.77, 21.55, 18.13, 15.33, 11.94, 8.722, 
// 									6.463, 4.84, 3.65, 2.775, 2.131, 1.662, 1.315, 1.045, 0.7477, 0.4863, 
// 									0.3226, 0.218, 0.1506, 0.1058, 0.07536, 0.05444, 0.03993, 0.02952, 0.01822, 0.009507, 
// 									0.005201, 0.003045, 0.001835, 0.001159, 6.256E-4, 2.931E-4, 1.483E-4, 8.094E-5, 4.742E-5, 2.804E-5, 
// 									1.786E-5, 1.141E-5, 7.72E-6, 4.512E-6, 2.206E-6, 1.211E-6, 6.881E-7, 4.184E-7, 2.687E-7, 1.721E-7, 
// 									1.071E-7, 8.669E-8, 4.078E-8, 2.977E-8, 1.818E-8, 1.057E-8 };
// 	Double_t p8462_d1x1y3_ysysminus[] = { 37.2, 25.6, 19.8, 15.4, 12.4, 10.2, 8.1, 6.39, 5.11, 
// 									 	  4.13, 3.36, 2.75, 2.27, 1.88, 1.57, 1.32, 1.12, 0.87, 0.637, 
// 									  	  0.472, 0.354, 0.267, 0.203, 0.156, 0.121, 0.096, 0.076, 0.0546, 0.0355, 
// 										  0.0236, 0.0159, 0.0110, 0.0077, 0.00550, 0.00398, 0.00292, 0.00216, 0.00133, 6.94E-4, 
// 										  3.80E-4, 2.22E-4, 1.34E-4, 8.5E-5, 4.57E-5, 2.14E-5, 1.08E-5, 5.91E-6, 3.46E-6, 1.99E-6, 
// 										  1.27E-6, 8.2E-7, 5.55E-7, 3.28E-7, 1.63E-7, 9.0E-8, 5.18E-8, 3.15E-8, 2.02E-8, 1.30E-8, 
// 										  0.81E-8, 0.656E-8, 3.10E-9, 2.77E-9, 1.41E-9, 0.83E-9 };
// 	Double_t p8462_d1x1y3_ysysplus[] = { 37.2, 25.6, 19.8, 15.4, 12.4, 10.2, 8.1, 6.39, 5.11, 
// 									 	  4.13, 3.36, 2.75, 2.27, 1.88, 1.57, 1.32, 1.12, 0.87, 0.637, 
// 									  	  0.472, 0.354, 0.267, 0.203, 0.156, 0.121, 0.096, 0.076, 0.0546, 0.0355, 
// 										  0.0236, 0.0159, 0.0110, 0.0077, 0.00550, 0.00398, 0.00292, 0.00216, 0.00133, 6.94E-4, 
// 										  3.80E-4, 2.22E-4, 1.34E-4, 8.5E-5, 4.57E-5, 2.14E-5, 1.08E-5, 5.91E-6, 3.46E-6, 1.99E-6, 
// 										  1.27E-6, 8.2E-7, 5.55E-7, 3.28E-7, 1.63E-7, 9.0E-8, 5.18E-8, 3.15E-8, 2.02E-8, 1.30E-8, 
// 										  0.81E-8, 0.656E-8, 3.10E-9, 2.77E-9, 1.41E-9, 0.83E-9 };
// 	Double_t p8462_d1x1y3_ystatminus[] = { 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.02, 0.01, 
// 		0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0, 0.002, 
// 		0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 4.0E-4, 3.0E-4, 
// 		2.0E-4, 2.0E-4, 2.0E-4, 1.0E-4, 9.0E-5, 7.0E-5, 6.0E-5, 5.0E-5, 3.0E-5, 1.7E-5, 
// 		1.2E-5, 8.0E-6, 6.0E-6, 5.0E-6, 2.3E-6, 1.5E-6, 1.0E-6, 7.1E-7, 5.2E-7, 3.8E-7, 
// 		2.9E-7, 2.3E-7, 1.79E-7, 9.2E-8, 6.1E-8, 4.3E-8, 3.1E-8, 2.31E-8, 1.78E-8, 1.38E-8, 
// 		1.05E-8, 9.11E-9, 6.1E-9, 3.49E-9, 2.32E-9, 1.65E-9 };
// 	Double_t p8462_d1x1y3_ystatplus[] = { 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.02, 0.01, 
// 		0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0, 0.002, 
// 		0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 4.0E-4, 3.0E-4, 
// 		2.0E-4, 2.0E-4, 2.0E-4, 1.0E-4, 9.0E-5, 7.0E-5, 6.0E-5, 5.0E-5, 3.0E-5, 1.7E-5, 
// 		1.2E-5, 8.0E-6, 6.0E-6, 5.0E-6, 2.3E-6, 1.5E-6, 1.0E-6, 7.1E-7, 5.2E-7, 3.8E-7, 
// 		2.9E-7, 2.3E-7, 1.79E-7, 9.2E-8, 6.1E-8, 4.3E-8, 3.1E-8, 2.31E-8, 1.78E-8, 1.38E-8, 
// 		1.05E-8, 9.11E-9, 6.1E-9, 3.49E-9, 2.32E-9, 1.65E-9 };
// 	int p8462_d1x1y3_numpoints = 65;
// 	TGraphAsymmErrors* graphChargedHadronsSysPP7TeV = new TGraphAsymmErrors(p8462_d1x1y3_numpoints, p8462_d1x1y3_xval, p8462_d1x1y3_yval, p8462_d1x1y3_xerrminus, p8462_d1x1y3_xerrplus, 
// 																			p8462_d1x1y3_ysysminus, p8462_d1x1y3_ysysplus);
// 	TGraphAsymmErrors* graphChargedHadronsStatPP7TeV = new TGraphAsymmErrors(p8462_d1x1y3_numpoints, p8462_d1x1y3_xval, p8462_d1x1y3_yval, p8462_d1x1y3_xerrminus, p8462_d1x1y3_xerrplus, 
// 																			 p8462_d1x1y3_ystatminus, p8462_d1x1y3_ystatplus);

	
	//--------------------- Write Files--------------------------------------------------------------   
	cout << "bis hier" << endl;
	TFile fileChargedPions(Form("ExternalInputPbPb/IdentifiedCharged/ChargedPionSpectraPbPb_%s.root",cStamp1) ,"RECREATE");
		histoChargedPionSpecLowPtStat0005->Write("histoChargedPionSpecLowPtStat0005");
		histoChargedPionSpecLowPtSyst0005->Write("histoChargedPionSpecLowPtSyst0005");
		histoChargedPionSpecLowPtStat0510->Write("histoChargedPionSpecLowPtStat0510");
		histoChargedPionSpecLowPtSyst0510->Write("histoChargedPionSpecLowPtSyst0510");
		histoChargedPionSpecLowPtStat0010->Write("histoChargedPionSpecLowPtStat0010");
		histoChargedPionSpecLowPtSyst0010->Write("histoChargedPionSpecLowPtSyst0010");
		histoChargedPionSpecLowPtStat1020->Write("histoChargedPionSpecLowPtStat1020");
		histoChargedPionSpecLowPtSyst1020->Write("histoChargedPionSpecLowPtSyst1020");
		histoChargedPionSpecLowPtStat2030->Write("histoChargedPionSpecLowPtStat2030");
		histoChargedPionSpecLowPtSyst2030->Write("histoChargedPionSpecLowPtSyst2030");
		histoChargedPionSpecLowPtStat3040->Write("histoChargedPionSpecLowPtStat3040");
		histoChargedPionSpecLowPtSyst3040->Write("histoChargedPionSpecLowPtSyst3040");
		histoChargedPionSpecLowPtStat2040->Write("histoChargedPionSpecLowPtStat2040");
		histoChargedPionSpecLowPtSyst2040->Write("histoChargedPionSpecLowPtSyst2040");
		histoChargedPionSpecLowPtStat4050->Write("histoChargedPionSpecLowPtStat4050");
		histoChargedPionSpecLowPtSyst4050->Write("histoChargedPionSpecLowPtSyst4050");
		histoChargedPionSpecLowPtStat5060->Write("histoChargedPionSpecLowPtStat5060");
		histoChargedPionSpecLowPtSyst5060->Write("histoChargedPionSpecLowPtSyst5060");
		histoChargedPionSpecLowPtStat4060->Write("histoChargedPionSpecLowPtStat4060");
		histoChargedPionSpecLowPtSyst4060->Write("histoChargedPionSpecLowPtSyst4060");
		histoChargedPionSpecLowPtStat6070->Write("histoChargedPionSpecLowPtStat6070");
		histoChargedPionSpecLowPtSyst6070->Write("histoChargedPionSpecLowPtSyst6070");
		histoChargedPionSpecLowPtStat7080->Write("histoChargedPionSpecLowPtStat7080");
		histoChargedPionSpecLowPtSyst7080->Write("histoChargedPionSpecLowPtSyst7080");
		histoChargedPionSpecLowPtStat6080->Write("histoChargedPionSpecLowPtStat6080");
		histoChargedPionSpecLowPtSyst6080->Write("histoChargedPionSpecLowPtSyst6080");
		
		histoChargedPionSpecHighPtStat0005->Write("histoChargedPionSpecHighPtStat0005");
		histoChargedPionSpecHighPtSyst0005->Write("histoChargedPionSpecHighPtSyst0005");
		histoChargedPionSpecHighPtStat0510->Write("histoChargedPionSpecHighPtStat0510");
		histoChargedPionSpecHighPtSyst0510->Write("histoChargedPionSpecHighPtSyst0510");
		histoChargedPionSpecHighPtStat0010->Write("histoChargedPionSpecHighPtStat0010");
		histoChargedPionSpecHighPtSyst0010->Write("histoChargedPionSpecHighPtSyst0010");
		histoChargedPionSpecHighPtStat1020->Write("histoChargedPionSpecHighPtStat1020");
		histoChargedPionSpecHighPtSyst1020->Write("histoChargedPionSpecHighPtSyst1020");
		histoChargedPionSpecHighPtStat2040->Write("histoChargedPionSpecHighPtStat2040");
		histoChargedPionSpecHighPtSyst2040->Write("histoChargedPionSpecHighPtSyst2040");
		histoChargedPionSpecHighPtStat4060->Write("histoChargedPionSpecHighPtStat4060");
		histoChargedPionSpecHighPtSyst4060->Write("histoChargedPionSpecHighPtSyst4060");
		histoChargedPionSpecHighPtStat6080->Write("histoChargedPionSpecHighPtStat6080");
		histoChargedPionSpecHighPtSyst6080->Write("histoChargedPionSpecHighPtSyst6080");
		      
		graphXiangoChargedPionNegStat0005->Write("graphChargedPionNegSpecXiangoStat0005");
		graphXiangoChargedPionNegSys0005->Write("graphChargedPionNegSpecXiangoSyst0005");
		graphXiangoChargedPionPosStat0005->Write("graphChargedPionPosSpecXiangoStat0005");
		graphXiangoChargedPionPosSys0005->Write("graphChargedPionPosSpecXiangoSyst0005");
		graphXiangoChargedPionNegStat0510->Write("graphChargedPionNegSpecXiangoStat0510");
		graphXiangoChargedPionNegSys0510->Write("graphChargedPionNegSpecXiangoSyst0510");
		graphXiangoChargedPionPosStat0510->Write("graphChargedPionPosSpecXiangoStat0510");
		graphXiangoChargedPionPosSys0510->Write("graphChargedPionPosSpecXiangoSyst0510");
		graphXiangoChargedPionNegStat1020->Write("graphChargedPionNegSpecXiangoStat1020");
		graphXiangoChargedPionNegSys1020->Write("graphChargedPionNegSpecXiangoSyst1020");
		graphXiangoChargedPionPosStat1020->Write("graphChargedPionPosSpecXiangoStat1020");
		graphXiangoChargedPionPosSys1020->Write("graphChargedPionPosSpecXiangoSyst1020");
		graphXiangoChargedPionNegStat2040->Write("graphChargedPionNegSpecXiangoStat2040");
		graphXiangoChargedPionNegSys2040->Write("graphChargedPionNegSpecXiangoSyst2040");
		graphXiangoChargedPionPosStat2040->Write("graphChargedPionPosSpecXiangoStat2040");
		graphXiangoChargedPionPosSys2040->Write("graphChargedPionPosSpecXiangoSyst2040");
		graphXiangoChargedPionNegStat4060->Write("graphChargedPionNegSpecXiangoStat4060");
		graphXiangoChargedPionNegSys4060->Write("graphChargedPionNegSpecXiangoSyst4060");
		graphXiangoChargedPionPosStat4060->Write("graphChargedPionPosSpecXiangoStat4060");
		graphXiangoChargedPionPosSys4060->Write("graphChargedPionPosSpecXiangoSyst4060");
		graphXiangoChargedPionNegStat6080->Write("graphChargedPionNegSpecXiangoStat6080");
		graphXiangoChargedPionNegSys6080->Write("graphChargedPionNegSpecXiangoSyst6080");
		graphXiangoChargedPionPosStat6080->Write("graphChargedPionPosSpecXiangoStat6080");
		graphXiangoChargedPionPosSys6080->Write("graphChargedPionPosSpecXiangoSyst6080");
   
		histoChargedPionRAAStat0005->Write("histoChargedPionRAAStat0005");
		histoChargedPionRAASyst0005->Write("histoChargedPionRAASyst0005");
		histoChargedPionRAAStat0510->Write("histoChargedPionRAAStat0510");
		histoChargedPionRAASyst0510->Write("histoChargedPionRAASyst0510");
		histoChargedPionRAAStat1020->Write("histoChargedPionRAAStat1020");
		histoChargedPionRAASyst1020->Write("histoChargedPionRAASyst1020");
		histoChargedPionRAAStat2040->Write("histoChargedPionRAAStat2040");
		histoChargedPionRAASyst2040->Write("histoChargedPionRAASyst2040");
		histoChargedPionRAAStat4060->Write("histoChargedPionRAAStat4060");
		histoChargedPionRAASyst4060->Write("histoChargedPionRAASyst4060");
		histoChargedPionRAAStat6080->Write("histoChargedPionRAAStat6080");
		histoChargedPionRAASyst6080->Write("histoChargedPionRAASyst6080");
		
	fileChargedPions.Close();
   
	TFile fileChargedKaons(Form("ExternalInputPbPb/IdentifiedCharged/KaonSpectraPbPb_%s.root",cStamp1) ,"RECREATE");
		histoChargedKaonSpecLowPtStat0005->Write("histoChargedKaonSpecLowPtStat0005");
		histoChargedKaonSpecLowPtSyst0005->Write("histoChargedKaonSpecLowPtSyst0005");
		histoChargedKaonSpecLowPtStat0510->Write("histoChargedKaonSpecLowPtStat0510");
		histoChargedKaonSpecLowPtSyst0510->Write("histoChargedKaonSpecLowPtSyst0510");
		histoChargedKaonSpecLowPtStat0010->Write("histoChargedKaonSpecLowPtStat0010");
		histoChargedKaonSpecLowPtSyst0010->Write("histoChargedKaonSpecLowPtSyst0010");
		histoChargedKaonSpecLowPtStat1020->Write("histoChargedKaonSpecLowPtStat1020");
		histoChargedKaonSpecLowPtSyst1020->Write("histoChargedKaonSpecLowPtSyst1020");
		histoChargedKaonSpecLowPtStat2030->Write("histoChargedKaonSpecLowPtStat2030");
		histoChargedKaonSpecLowPtSyst2030->Write("histoChargedKaonSpecLowPtSyst2030");
		histoChargedKaonSpecLowPtStat3040->Write("histoChargedKaonSpecLowPtStat3040");
		histoChargedKaonSpecLowPtSyst3040->Write("histoChargedKaonSpecLowPtSyst3040");
		histoChargedKaonSpecLowPtStat2040->Write("histoChargedKaonSpecLowPtStat2040");
		histoChargedKaonSpecLowPtSyst2040->Write("histoChargedKaonSpecLowPtSyst2040");
		histoChargedKaonSpecLowPtStat4050->Write("histoChargedKaonSpecLowPtStat4050");
		histoChargedKaonSpecLowPtSyst4050->Write("histoChargedKaonSpecLowPtSyst4050");
		histoChargedKaonSpecLowPtStat5060->Write("histoChargedKaonSpecLowPtStat5060");
		histoChargedKaonSpecLowPtSyst5060->Write("histoChargedKaonSpecLowPtSyst5060");
		histoChargedKaonSpecLowPtStat4060->Write("histoChargedKaonSpecLowPtStat4060");
		histoChargedKaonSpecLowPtSyst4060->Write("histoChargedKaonSpecLowPtSyst4060");
		histoChargedKaonSpecLowPtStat6070->Write("histoChargedKaonSpecLowPtStat6070");
		histoChargedKaonSpecLowPtSyst6070->Write("histoChargedKaonSpecLowPtSyst6070");
		histoChargedKaonSpecLowPtStat7080->Write("histoChargedKaonSpecLowPtStat7080");
		histoChargedKaonSpecLowPtSyst7080->Write("histoChargedKaonSpecLowPtSyst7080");
		histoChargedKaonSpecLowPtStat6080->Write("histoChargedKaonSpecLowPtStat6080");
		histoChargedKaonSpecLowPtSyst6080->Write("histoChargedKaonSpecLowPtSyst6080");      
		
		graphXiangoChargedKaonNegStat0005->Write("graphChargedKaonNegSpecXiangoStat0005");
		graphXiangoChargedKaonNegSys0005->Write("graphChargedKaonNegSpecXiangoSyst0005");
		graphXiangoChargedKaonPosStat0005->Write("graphChargedKaonPosSpecXiangoStat0005");
		graphXiangoChargedKaonPosSys0005->Write("graphChargedKaonPosSpecXiangoSyst0005");
		graphXiangoChargedKaonNegStat0510->Write("graphChargedKaonNegSpecXiangoStat0510");
		graphXiangoChargedKaonNegSys0510->Write("graphChargedKaonNegSpecXiangoSyst0510");
		graphXiangoChargedKaonPosStat0510->Write("graphChargedKaonPosSpecXiangoStat0510");
		graphXiangoChargedKaonPosSys0510->Write("graphChargedKaonPosSpecXiangoSyst0510");
		graphXiangoChargedKaonNegStat1020->Write("graphChargedKaonNegSpecXiangoStat1020");
		graphXiangoChargedKaonNegSys1020->Write("graphChargedKaonNegSpecXiangoSyst1020");
		graphXiangoChargedKaonPosStat1020->Write("graphChargedKaonPosSpecXiangoStat1020");
		graphXiangoChargedKaonPosSys1020->Write("graphChargedKaonPosSpecXiangoSyst1020");
		graphXiangoChargedKaonNegStat2040->Write("graphChargedKaonNegSpecXiangoStat2040");
		graphXiangoChargedKaonNegSys2040->Write("graphChargedKaonNegSpecXiangoSyst2040");
		graphXiangoChargedKaonPosStat2040->Write("graphChargedKaonPosSpecXiangoStat2040");
		graphXiangoChargedKaonPosSys2040->Write("graphChargedKaonPosSpecXiangoSyst2040");
		graphXiangoChargedKaonNegStat4060->Write("graphChargedKaonNegSpecXiangoStat4060");
		graphXiangoChargedKaonNegSys4060->Write("graphChargedKaonNegSpecXiangoSyst4060");
		graphXiangoChargedKaonPosStat4060->Write("graphChargedKaonPosSpecXiangoStat4060");
		graphXiangoChargedKaonPosSys4060->Write("graphChargedKaonPosSpecXiangoSyst4060");
		graphXiangoChargedKaonNegStat6080->Write("graphChargedKaonNegSpecXiangoStat6080");
		graphXiangoChargedKaonNegSys6080->Write("graphChargedKaonNegSpecXiangoSyst6080");
		graphXiangoChargedKaonPosStat6080->Write("graphChargedKaonPosSpecXiangoStat6080");
		graphXiangoChargedKaonPosSys6080->Write("graphChargedKaonPosSpecXiangoSyst6080");
	
		histoChargedKaonSpecHighPtStat0005->Write("histoChargedKaonSpecHighPtStat0005");
		histoChargedKaonSpecHighPtSyst0005->Write("histoChargedKaonSpecHighPtSyst0005");
		histoChargedKaonSpecHighPtStat0510->Write("histoChargedKaonSpecHighPtStat0510");
		histoChargedKaonSpecHighPtSyst0510->Write("histoChargedKaonSpecHighPtSyst0510");
		histoChargedKaonSpecHighPtStat0010->Write("histoChargedKaonSpecHighPtStat0010");
		histoChargedKaonSpecHighPtSyst0010->Write("histoChargedKaonSpecHighPtSyst0010");
		histoChargedKaonSpecHighPtStat1020->Write("histoChargedKaonSpecHighPtStat1020");
		histoChargedKaonSpecHighPtSyst1020->Write("histoChargedKaonSpecHighPtSyst1020");
		histoChargedKaonSpecHighPtStat2040->Write("histoChargedKaonSpecHighPtStat2040");
		histoChargedKaonSpecHighPtSyst2040->Write("histoChargedKaonSpecHighPtSyst2040");
		histoChargedKaonSpecHighPtStat4060->Write("histoChargedKaonSpecHighPtStat4060");
		histoChargedKaonSpecHighPtSyst4060->Write("histoChargedKaonSpecHighPtSyst4060");
		histoChargedKaonSpecHighPtStat6080->Write("histoChargedKaonSpecHighPtStat6080");
		histoChargedKaonSpecHighPtSyst6080->Write("histoChargedKaonSpecHighPtSyst6080");

		histoNeutralKaonSpecStat0005->Write("histoNeutralKaonSpecStat0005");
		histoNeutralKaonSpecSyst0005->Write("histoNeutralKaonSpecSyst0005");
		histoNeutralKaonSpecStat0010->Write("histoNeutralKaonSpecStat0010");
		histoNeutralKaonSpecSyst0010->Write("histoNeutralKaonSpecSyst0010");
		histoNeutralKaonSpecStat0510->Write("histoNeutralKaonSpecStat0510");
		histoNeutralKaonSpecSyst0510->Write("histoNeutralKaonSpecSyst0510");
		histoNeutralKaonSpecStat1020->Write("histoNeutralKaonSpecStat1020");
		histoNeutralKaonSpecSyst1020->Write("histoNeutralKaonSpecSyst1020");
		histoNeutralKaonSpecStat2040->Write("histoNeutralKaonSpecStat2040");
		histoNeutralKaonSpecSyst2040->Write("histoNeutralKaonSpecSyst2040");
		histoNeutralKaonSpecStat4060->Write("histoNeutralKaonSpecStat4060");
		histoNeutralKaonSpecSyst4060->Write("histoNeutralKaonSpecSyst4060");
		histoNeutralKaonSpecStat6080->Write("histoNeutralKaonSpecStat6080");
		histoNeutralKaonSpecSyst6080->Write("histoNeutralKaonSpecSyst6080");
		
		histoChargedKaonRAAStat0005->Write("histoChargedKaonRAAStat0005");
		histoChargedKaonRAASyst0005->Write("histoChargedKaonRAASyst0005");
		histoChargedKaonRAAStat0510->Write("histoChargedKaonRAAStat0510");
		histoChargedKaonRAASyst0510->Write("histoChargedKaonRAASyst0510");
		histoChargedKaonRAAStat1020->Write("histoChargedKaonRAAStat1020");
		histoChargedKaonRAASyst1020->Write("histoChargedKaonRAASyst1020");
		histoChargedKaonRAAStat2040->Write("histoChargedKaonRAAStat2040");
		histoChargedKaonRAASyst2040->Write("histoChargedKaonRAASyst2040");
		histoChargedKaonRAAStat4060->Write("histoChargedKaonRAAStat4060");
		histoChargedKaonRAASyst4060->Write("histoChargedKaonRAASyst4060");
		histoChargedKaonRAAStat6080->Write("histoChargedKaonRAAStat6080");
		histoChargedKaonRAASyst6080->Write("histoChargedKaonRAASyst6080");
		
	fileChargedKaons.Close();

	TFile fileChargedHadronspp(Form("ExternalInput/UnidentifiedCharged/ChargedHadrinSpectraPP_%s.root",cStamp1) ,"RECREATE");
		graphChargedHadronsStatPP2760GeV->Write("graphChargedHadronsStatPP2760GeV");
		graphChargedHadronsSysPP2760GeV->Write("graphChargedHadronsSysPP2760GeV");

	fileChargedHadronspp.Close();
	cout << "bis hier" << endl;
	
	
	
	cout << "bis hier" << endl;
	TFile fileChargedPionspp(Form("ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_%s.root",cStamp1) ,"RECREATE");
		histoChargedPionSpecPubStatPP->Write("histoChargedPionSpecPubStat2760GeV");
		histoChargedPionSpecPubSystPP->Write("histoChargedPionSpecPubSyst2760GeV");
		histoChargedKaonSpecPubStatPP->Write("histoChargedKaonSpecPubStat2760GeV");
		histoChargedKaonSpecPubSystPP->Write("histoChargedKaonSpecPubSyst2760GeV");
		histoProtonSpecPubStatPP->Write("histoProtonPionSpecPubStat2760GeV");
		histoProtonSpecPubSystPP->Write("histoProtonSpecPubSyst2760GeV");

		histoChargedPionSpecHighPtStatPP->Write("histoChargedPionSpecHighPtStat2760GeV");
		histoChargedPionSpecHighPtSystPP->Write("histoChargedPionSpecHighPtSyst2760GeV");

		graphXiangoChargedPionNegStatPP2760GeV->Write("graphChargedPionNegSpecXiangoStat2760GeV");
		graphXiangoChargedPionNegSysPP2760GeV->Write("graphChargedPionNegSpecXiangoSyst2760GeV");
		graphXiangoChargedPionPosStatPP2760GeV->Write("graphChargedPionPosSpecXiangoStat2760GeV");
		graphXiangoChargedPionPosSysPP2760GeV->Write("graphChargedPionPosSpecXiangoSyst2760GeV");
	
		graphXiangoChargedPionNegStatPP7TeV->Write("graphChargedPionNegSpecXiangoStat7TeV");
		graphXiangoChargedPionNegSysPP7TeV->Write("graphChargedPionNegSpecXiangoSyst7TeV");
		graphXiangoChargedPionPosStatPP7TeV->Write("graphChargedPionPosSpecXiangoStat7TeV");
		graphXiangoChargedPionPosSysPP7TeV->Write("graphChargedPionPosSpecXiangoSyst7TeV");

		graphXiangoChargedKaonNegStatPP2760GeV->Write("graphChargedKaonNegSpecXiangoStat2760GeV");
		graphXiangoChargedKaonNegSysPP2760GeV->Write("graphChargedKaonNegSpecXiangoSyst2760GeV");
		graphXiangoChargedKaonPosStatPP2760GeV->Write("graphChargedKaonPosSpecXiangoStat2760GeV");
		graphXiangoChargedKaonPosSysPP2760GeV->Write("graphChargedKaonPosSpecXiangoSyst2760GeV");
	
		graphXiangoChargedKaonNegStatPP7TeV->Write("graphChargedKaonNegSpecXiangoStat7TeV");
		graphXiangoChargedKaonNegSysPP7TeV->Write("graphChargedKaonNegSpecXiangoSyst7TeV");
		graphXiangoChargedKaonPosStatPP7TeV->Write("graphChargedKaonPosSpecXiangoStat7TeV");
		graphXiangoChargedKaonPosSysPP7TeV->Write("graphChargedKaonPosSpecXiangoSyst7TeV");
		
		histoChargedPionSpecLowPtStat2760GeVCMS->Write("histoChargedPionSpecLowPtStat2760GeVCMS");
		histoChargedPionSpecLowPtSys2760GeVCMS->Write("histoChargedPionSpecLowPtSys2760GeVCMS");
		histoChargedPionSpecLowPtStatPP2760GeV->Write("histoChargedPionSpecLowPtStatPP2760GeV");
		histoChargedPionSpecLowPtSysPP2760GeV->Write("histoChargedPionSpecLowPtSysPP2760GeV");
		
		histoChargedPionSpecLowPtStat900GeVCMS->Write("histoChargedPionSpecLowPtStat900GeVCMS");
		histoChargedPionSpecLowPtSys900GeVCMS->Write("histoChargedPionSpecLowPtSys900GeVCMS");
		histoChargedPionSpecLowPtStat900GeVALICE->Write("histoChargedPionSpecLowPtStat900GeVALICE");
		histoChargedPionSpecLowPtSys900GeVALICE->Write("histoChargedPionSpecLowPtSys900GeVALICE");
		
		histoChargedPionSpecLowPtStat7TeVCMS->Write("histoChargedPionSpecLowPtStat7TeVCMS");
		histoChargedPionSpecLowPtSys7TeVCMS->Write("histoChargedPionSpecLowPtSys7TeVCMS");
		histoChargedPionSpecHighPtStat7TeVALICE->Write("histoChargedPionSpecHighPtStat7TeVALICE");
		graphChargedPionSpecHighPtSys7TeVALICE->Write("graphChargedPionSpecHighPtSys7TeVALICE");

		histoChargedPionSpecLowPtSys7TeVALICE->Write("histoChargedPionSpecLowPtSys7TeVALICE");
		histoChargedPionSpecLowPtStat7TeVALICE->Write("histoChargedPionSpecLowPtStat7TeVALICE");
	fileChargedPionspp.Close();
	cout << "bis hier" << endl;
  
	TFile fileChargedPionspPb(Form("ExternalInputpPb/ChargedPionSpectrapPb_%s.root",cStamp1) ,"RECREATE");
		histoChargedPionSpecLowPtSyspPb->Write("histoChargedPionSpecLowPtSyspPb");
		histoChargedPionSpecLowPtStatpPb->Write("histoChargedPionSpecLowPtStatpPb");
		histoChargedPionSpecLowPtSyspPb0020->Write("histoChargedPionSpecLowPtSyspPb0020");
		histoChargedPionSpecLowPtStatpPb0020->Write("histoChargedPionSpecLowPtStatpPb0020");
		histoChargedPionSpecLowPtSyspPb2040->Write("histoChargedPionSpecLowPtSyspPb2040");
		histoChargedPionSpecLowPtStatpPb2040->Write("histoChargedPionSpecLowPtStatpPb2040");
		histoChargedPionSpecLowPtSyspPb4060->Write("histoChargedPionSpecLowPtSyspPb4060");
		histoChargedPionSpecLowPtStatpPb4060->Write("histoChargedPionSpecLowPtStatpPb4060");
		histoChargedPionSpecLowPtSyspPb6080->Write("histoChargedPionSpecLowPtSyspPb6080");
		histoChargedPionSpecLowPtStatpPb6080->Write("histoChargedPionSpecLowPtStatpPb6080");
		histoChargedPionSpecLowPtSyspPb80100->Write("histoChargedPionSpecLowPtSyspPb80100");
		histoChargedPionSpecLowPtStatpPb80100->Write("histoChargedPionSpecLowPtStatpPb80100");
		histopPbChargedPionSpecFullpTSys->Write("histoChargedPionSpecFullPtSyspPb");
		histopPbChargedPionSpecFullpTStat->Write("histoChargedPionSpecFullPtStatpPb");
		histopPbChargedPionSpecFullpTSys0020->Write("histoChargedPionSpecFullPtSyspPb0020");
		histopPbChargedPionSpecFullpTStat0020->Write("histoChargedPionSpecFullPtStatpPb0020");
		histopPbChargedPionSpecFullpTSys2040->Write("histoChargedPionSpecFullPtSyspPb2040");
		histopPbChargedPionSpecFullpTStat2040->Write("histoChargedPionSpecFullPtStatpPb2040");
		histopPbChargedPionSpecFullpTSys4060->Write("histoChargedPionSpecFullPtSyspPb4060");
		histopPbChargedPionSpecFullpTStat4060->Write("histoChargedPionSpecFullPtStatpPb4060");
		histopPbChargedPionSpecFullpTSys6080->Write("histoChargedPionSpecFullPtSyspPb6080");
		histopPbChargedPionSpecFullpTStat6080->Write("histoChargedPionSpecFullPtStatpPb6080");
		histopPbChargedPionSpecFullpTSys60100->Write("histoChargedPionSpecFullPtSyspPb60100");
		histopPbChargedPionSpecFullpTStat60100->Write("histoChargedPionSpecFullPtStatpPb60100");
	fileChargedPionspPb.Close();
  
}

