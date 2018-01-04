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


void ProduceV0FindingEffPlots(TString filename = "GammaConvV1_66.root", TString nameFolder ="GammaConvV1_66", TString EventCutString = "5010001", TString outputDir = ""){
	
	TFile *file = new TFile(filename.Data()); 
	TList *list = (TList*) file->Get(nameFolder.Data());
	TList *HistosV0Finding = (TList*)list->FindObject(Form("V0FindingEfficiencyInput_%s", EventCutString.Data()));

	//loading histos
	TH2F* histoConvGamma_Pt_R 	= (TH2F*)HistosV0Finding->FindObject("MCconvGamma_Pt_R");   
	TH2F* histoConvGamma_Pt_Eta = (TH2F*)HistosV0Finding->FindObject("MCconvGamma_Pt_Eta");  
	TH2F* histoConvGamma_Pt_Phi = (TH2F*)HistosV0Finding->FindObject("MCconvGamma_Pt_Phi"); 

	TH2F* histoRecConvGamma_Pt_R 	= (TH2F*)HistosV0Finding->FindObject("RecMCconvGamma_Pt_R");   
	TH2F* histoRecConvGamma_Pt_Eta 	= (TH2F*)HistosV0Finding->FindObject("RecMCconvGamma_Pt_Eta");  
	TH2F* histoRecConvGamma_Pt_Phi 	= (TH2F*)HistosV0Finding->FindObject("RecMCconvGamma_Pt_Phi");

	TH1D* histoRecConvGammaMult_Pt 	= (TH1D*)HistosV0Finding->FindObject("RecMCconvGammaMulti_Pt");   
	TH1D* histoRecConvGammaMult_R 	= (TH1D*)HistosV0Finding->FindObject("RecMCconvGammaMulti_R");  
	TH1D* histoRecConvGammaMult_Phi = (TH1D*)HistosV0Finding->FindObject("RecMCconvGammaMulti_Phi");
	
	//projecting histos
	TH1D* histoConvGamma_Pt 			= (TH1D*)histoConvGamma_Pt_R->ProjectionX("histoConvGamma_Pt"); //restricted to eta cut range
	TH1D* histoConvGamma_PtlargerEta	= (TH1D*)histoConvGamma_Pt_Eta->ProjectionX("histoConvGamma_PtlargerEta"); //full eta range
		
	TH1D* histoConvGamma_R 				= (TH1D*)histoConvGamma_Pt_R->ProjectionY("histoConvGamma_R"); //full pT range
// 	cout << "proj R in 0. < t <1.0: " <<  histoConvGamma_Pt_R->GetXaxis()->FindBin(0.) << " to " << histoConvGamma_Pt_R->GetXaxis()->FindBin(1.0) << endl;
	TH1D* histoConvGamma_R_lowpT 		= (TH1D*)histoConvGamma_Pt_R->ProjectionY("histoConvGamma_R_lowpT",histoConvGamma_Pt_R->GetXaxis()->FindBin(0.),histoConvGamma_Pt_R->GetXaxis()->FindBin(1.0)); // 0. < pT < 1. GeV/c 
// 	cout << "proj R in 0. < t <1.0: " <<  histoConvGamma_Pt_R->GetXaxis()->FindBin(1.) << " to " << histoConvGamma_Pt_R->GetXaxis()->FindBin(4.0) << endl;
	TH1D* histoConvGamma_R_midpT 		= (TH1D*)histoConvGamma_Pt_R->ProjectionY("histoConvGamma_R_midpT",histoConvGamma_Pt_R->GetXaxis()->FindBin(1.0),histoConvGamma_Pt_R->GetXaxis()->FindBin(4.0)); // 1.0 < pT < 4. GeV/c 
// 	cout << "proj R in 0. < t <1.0: " <<  histoConvGamma_Pt_R->GetXaxis()->FindBin(4.) << " to " << histoConvGamma_Pt_R->GetNbinsX() << endl;
	TH1D* histoConvGamma_R_highpT		= (TH1D*)histoConvGamma_Pt_R->ProjectionY("histoConvGamma_R_highpT",histoConvGamma_Pt_R->GetXaxis()->FindBin(4.0),histoConvGamma_Pt_R->GetNbinsX()); // pT > 4. GeV/c
	
	TH1D* histoConvGamma_Eta			= (TH1D*)histoConvGamma_Pt_Eta->ProjectionY("histoConvGamma_Eta");
	TH1D* histoConvGamma_Phi 			= (TH1D*)histoConvGamma_Pt_Phi->ProjectionY("histoConvGamma_Phi");

	TH1D* histoRecConvGamma_Pt 			= (TH1D*)histoRecConvGamma_Pt_R->ProjectionX("histoRecConvGamma_Pt"); //restricted to eta cut range
	TH1D* histoRecConvGamma_PtlargerEta = (TH1D*)histoRecConvGamma_Pt_Eta->ProjectionX("histoRecConvGamma_PtlargerEta"); //full eta range
	TH1D* histoRecConvGamma_R 			= (TH1D*)histoRecConvGamma_Pt_R->ProjectionY("histoRecConvGamma_R"); //full pT range
	// 	cout << "proj R in 0. < t <1.0: " <<  histoRecConvGamma_Pt_R->GetXaxis()->FindBin(0.) << " to " << histoRecConvGamma_Pt_R->GetXaxis()->FindBin(1.0) << endl;
	TH1D* histoRecConvGamma_R_lowpT 		= (TH1D*)histoRecConvGamma_Pt_R->ProjectionY("histoRecConvGamma_R_lowpT",histoRecConvGamma_Pt_R->GetXaxis()->FindBin(0.),histoRecConvGamma_Pt_R->GetXaxis()->FindBin(1.0)); // 0. < pT < 1. GeV/c 
// 	cout << "proj R in 0. < t <1.0: " <<  histoRecConvGamma_Pt_R->GetXaxis()->FindBin(1.) << " to " << histoRecConvGamma_Pt_R->GetXaxis()->FindBin(4.0) << endl;
	TH1D* histoRecConvGamma_R_midpT 		= (TH1D*)histoRecConvGamma_Pt_R->ProjectionY("histoRecConvGamma_R_midpT",histoRecConvGamma_Pt_R->GetXaxis()->FindBin(1.0),histoRecConvGamma_Pt_R->GetXaxis()->FindBin(4.0)); // 1.0 < pT < 4. GeV/c 
// 	cout << "proj R in 0. < t <1.0: " <<  histoRecConvGamma_Pt_R->GetXaxis()->FindBin(4.) << " to " << histoRecConvGamma_Pt_R->GetNbinsX() << endl;
	TH1D* histoRecConvGamma_R_highpT		= (TH1D*)histoRecConvGamma_Pt_R->ProjectionY("histoRecConvGamma_R_highpT",histoRecConvGamma_Pt_R->GetXaxis()->FindBin(4.0),histoRecConvGamma_Pt_R->GetNbinsX()); // pT > 4. GeV/c
	
	TH1D* histoRecConvGamma_Eta 		= (TH1D*)histoRecConvGamma_Pt_Eta->ProjectionY("histoRecConvGamma_Eta");
	TH1D* histoRecConvGamma_Phi 		= (TH1D*)histoRecConvGamma_Pt_Phi->ProjectionY("histoRecConvGamma_Phi");

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
	TH1D* histoRecConvGamma_Pt_rebin = (TH1D*) histoRecConvGamma_Pt->Rebin(22,"recPteff",ptbinning);
	TH1D* histoConvGamma_Pt_rebin=(TH1D*)  histoConvGamma_Pt->Rebin(22,"totPteff",ptbinning);
	TH1D* histoRatioRecConvGamma_Pt 		 = (TH1D*)histoRecConvGamma_Pt_rebin->Clone("histoRatioRecConvGamma_Pt");
	histoRatioRecConvGamma_Pt->Sumw2();
	histoRatioRecConvGamma_Pt->Divide(histoRecConvGamma_Pt_rebin,histoConvGamma_Pt_rebin,1.,1.,"B");
	
	TH1D* histoRecConvGamma_PtlargerEta_rebin = (TH1D*) histoRecConvGamma_PtlargerEta->Rebin(22,"recPtlargeEta",ptbinning);
	TH1D* histoConvGamma_PtlargerEta_rebin=(TH1D*)  histoConvGamma_PtlargerEta->Rebin(22,"totPtlargeEta",ptbinning);
	TH1D* histoRatioRecConvGamma_PtlargerEta = (TH1D*)histoRecConvGamma_PtlargerEta_rebin->Clone("histoRatioRecConvGamma_PtlargerEta");
	histoRatioRecConvGamma_PtlargerEta->Sumw2();
	histoRatioRecConvGamma_PtlargerEta->Divide(histoRecConvGamma_PtlargerEta_rebin,histoConvGamma_PtlargerEta_rebin,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R = (TH1D*)histoRecConvGamma_R->Clone("histoRatioRecConvGamma_R");
	histoRatioRecConvGamma_R->Sumw2();
	histoRatioRecConvGamma_R->Divide(histoRatioRecConvGamma_R,histoConvGamma_R,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_R_lowpT = (TH1D*)histoRecConvGamma_R_lowpT->Clone("histoRatioRecConvGamma_R_lowpT");
	histoRatioRecConvGamma_R_lowpT->Sumw2();
	histoRatioRecConvGamma_R_lowpT->Divide(histoRatioRecConvGamma_R_lowpT,histoConvGamma_R_lowpT,1.,1.,"B");	
	
	TH1D* histoRatioRecConvGamma_R_midpT = (TH1D*)histoRecConvGamma_R_midpT->Clone("histoRatioRecConvGamma_R_midpT");
	histoRatioRecConvGamma_R_midpT->Sumw2();
	histoRatioRecConvGamma_R_midpT->Divide(histoRatioRecConvGamma_R_midpT,histoConvGamma_R_midpT,1.,1.,"B");
	
	TH1D* histoRatioRecConvGamma_R_highpT = (TH1D*)histoRecConvGamma_R_highpT->Clone("histoRatioRecConvGamma_R_highpT");
	histoRatioRecConvGamma_R_highpT->Sumw2();
	histoRatioRecConvGamma_R_highpT->Divide(histoRatioRecConvGamma_R_highpT,histoConvGamma_R_highpT,1.,1.,"B");
	
	
	
	TH1D* histoRatioRecConvGamma_Eta 		 = (TH1D*)histoRecConvGamma_Eta->Clone("histoRatioRecConvGamma_Eta");
	histoRatioRecConvGamma_Eta->Sumw2();
	histoRatioRecConvGamma_Eta->Divide(histoRatioRecConvGamma_Eta,histoConvGamma_Eta,1.,1.,"B");

	TH1D* histoRatioRecConvGamma_Phi 		 = (TH1D*)histoRecConvGamma_Phi->Clone("histoRatioRecConvGamma_Phi");
	histoRatioRecConvGamma_Phi->Sumw2();
	histoRatioRecConvGamma_Phi->Divide(histoRatioRecConvGamma_Phi,histoConvGamma_Phi,1.,1.,"B");
	
	//calculating the tot number of photons = 1st photon counted + 2nd time the same photon has been counted
	TH1D* histoTotalRecConvGamma_Pt = (TH1D*)histoRecConvGamma_Pt->Clone("histoRatioRecConvGamma_Pt");
	histoRecConvGamma_Pt->Add(histoRecConvGammaMult_Pt);
	TH1D* histoTotalRecConvGamma_R = (TH1D*)histoRecConvGamma_R->Clone("histoRatioRecConvGamma_R");
	histoTotalRecConvGamma_R->Add(histoRecConvGammaMult_R);
	TH1D* histoTotalRecConvGamma_Phi = (TH1D*)histoRecConvGamma_Phi->Clone("histoRatioRecConvGamma_Phi");
	histoTotalRecConvGamma_Phi->Add(histoRecConvGammaMult_Phi);
	

	
	//making ratio with the Mult histos = double counting ratio
	TH1D* histoRecConvGammaMult_Pt_rebin = (TH1D*) histoRecConvGammaMult_Pt->Rebin(22,"recPt",ptbinning);
	TH1D* histoTotalRecConvGamma_Pt_rebin=(TH1D*)  histoTotalRecConvGamma_Pt->Rebin(22,"totPt",ptbinning);
	TH1D* histoRatioRecConvGammaMult_Pt 	= (TH1D*)histoRecConvGammaMult_Pt_rebin->Clone("histoRatioRecConvGammaMult_Pt");
	histoRatioRecConvGammaMult_Pt->Sumw2();
	histoRatioRecConvGammaMult_Pt->Divide(histoRecConvGammaMult_Pt_rebin,histoTotalRecConvGamma_Pt_rebin,1.,1.,"B");

	TH1D* histoRecConvGammaMult_R_rebin = (TH1D*) histoRecConvGammaMult_R->Rebin(86,"recR",Rbinning);
	TH1D* histoTotalRecConvGamma_R_rebin=(TH1D*)  histoTotalRecConvGamma_R->Rebin(86,"recR",Rbinning);	
	TH1D* histoRatioRecConvGammaMult_R 		= (TH1D*)histoRecConvGammaMult_R_rebin->Clone("histoRatioRecConvGammaMult_R");
	histoRatioRecConvGammaMult_R->Sumw2();
	histoRatioRecConvGammaMult_R->Divide(histoRecConvGammaMult_R_rebin,histoTotalRecConvGamma_R_rebin,1.,1.,"B");
	
	TH1D* histoRatioRecConvGammaMult_Phi 	= (TH1D*)histoRecConvGammaMult_Phi->Clone("histoRatioRecConvGammaMult_Phi");
	histoRatioRecConvGammaMult_Phi->Sumw2();
	histoRatioRecConvGammaMult_Phi->Divide(histoRatioRecConvGammaMult_Phi,histoTotalRecConvGamma_Phi,1.,1.,"B");
	
	TFile* output = new TFile(Form("%s/HistosV0FindingQA_%s.root",outputDir.Data(), nameFolder.Data()),"RECREATE");
	output->cd();
	
	histoConvGamma_Pt->Write();
	histoRecConvGamma_Pt->Write();
	histoRatioRecConvGamma_Pt->Write();
	histoRatioRecConvGamma_PtlargerEta->Write();
	histoRatioRecConvGamma_R->Write();
	histoConvGamma_R_highpT->Write();
	histoConvGamma_R_lowpT->Write();
	histoConvGamma_R_midpT->Write();
	histoRecConvGamma_R_highpT->Write();
	histoRecConvGamma_R_lowpT->Write();
	histoRecConvGamma_R_midpT->Write();
	histoRatioRecConvGamma_R_lowpT->Write();
	histoRatioRecConvGamma_R_midpT->Write();
	histoRatioRecConvGamma_R_highpT->Write();
	histoRatioRecConvGamma_Eta->Write();
	histoRatioRecConvGamma_Phi->Write();
	histoRatioRecConvGammaMult_Pt->Write();
	histoRatioRecConvGammaMult_R->Write();
	histoRatioRecConvGammaMult_Phi->Write();
	
	output->Write();
	output->Close();
	
	
// 	TH2F* histoConvGamma_Pt_R = (TH1D*)HistosV0Finding->FindObject("V0 Multiplicity");  

	
}
	
	