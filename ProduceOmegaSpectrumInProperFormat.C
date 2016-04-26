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

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

void ProduceOmegaSpectrumInProperFormat(){	
	
	TFile* fileOmegaSpectrum 				= 	new TFile("ExternalInput/PHOS/7TeV/PHOS_pp_omega_7TeV_07082012.root");
		TGraphErrors* omega_stat_err 		= (TGraphErrors*)fileOmegaSpectrum->Get("omega_stat_err");
		TGraphErrors* omega_syst_err 		= (TGraphErrors*)fileOmegaSpectrum->Get("omega_syst_err");
		TGraphErrors* omega_to_pi0_stat_err = (TGraphErrors*)fileOmegaSpectrum->Get("omega_to_pi0_stat_err");
		TGraphErrors* omega_to_pi0_stat_err = (TGraphErrors*)fileOmegaSpectrum->Get("omega_to_pi0_stat_err");

	Double_t binEdges[8] 					= {2.,4., 6., 8., 10., 12., 14.,17.};
	Double_t* xValuesOmega 					= omega_stat_err->GetX();
	Double_t* yValuesOmega 					= omega_stat_err->GetY();
	Double_t* yErrorsStatOmega 				= omega_stat_err->GetEY();
	Double_t* yErrorsSystOmega 				= omega_syst_err->GetEY();
	Int_t nBins 							= omega_stat_err->GetN();
	Double_t xErrorsLow[10];
	Double_t xErrorsHigh[10];
	Double_t yErrorsComb[10];
	for (Int_t i = 0; i < nBins; i++){
			xErrorsLow[i] 					= xValuesOmega[i] - binEdges[i];
			xErrorsHigh[i] 					= binEdges[i+1] - xValuesOmega[i];
			cout <<i << "\t" << xValuesOmega[i]<< "\t"  << xErrorsLow[i] << "\t" << xErrorsHigh[i] << endl;
			yErrorsComb[i] 					= TMath::Sqrt(TMath::Power(yErrorsStatOmega[i],2.)+TMath::Power(yErrorsSystOmega[i],2.));
	}

	TGraphAsymmErrors* graphOmegaStat 		= new TGraphAsymmErrors(nBins,xValuesOmega,yValuesOmega,xErrorsLow,xErrorsHigh,yErrorsStatOmega,yErrorsStatOmega);
	graphOmegaStat->Print();
	TGraphAsymmErrors* graphOmegaSyst 		= new TGraphAsymmErrors(nBins,xValuesOmega,yValuesOmega,xErrorsLow,xErrorsHigh,yErrorsSystOmega,yErrorsSystOmega);
	graphOmegaSyst->Print();
	TGraphAsymmErrors* graphOmegaComb 		= new TGraphAsymmErrors(nBins,xValuesOmega,yValuesOmega,xErrorsLow,xErrorsHigh,yErrorsComb,yErrorsComb);
	graphOmegaComb->Print();
// 
	TFile* fileOmegaSpectrumNew 			= new TFile("ExternalInput/PHOS/7TeV/PHOS_pp_omega_7TeV_07082012.root","RECREATE");
		omega_stat_err->Write("omega_stat_err");
		graphOmegaStat->Write("graphOmegaStat");
		omega_syst_err->Write("omega_syst_err");
		graphOmegaSyst->Write("graphOmegaSyst");
		graphOmegaComb->Write("graphOmegaComb");
		omega_to_pi0_stat_err->Write("omega_to_pi0_stat_err");
		omega_to_pi0_stat_err->Write("omega_to_pi0_stat_err");
	fileOmegaSpectrumNew->Close();
}
