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
#include <string>
#include "TGaxis.h"
#include "TFractionFitter.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "THStack.h"
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
#include "TEllipse.h"
#include "DirectPhotonFlowFunctions.h"


void  projectionPt_v9(  Int_t Trainconfig = 55,
                        TString CentralityLow = "20",
                        TString CentralityHigh = "40",
                        TString Cutnumber = "52400013_00200009007000008250400000",
                        Bool_t kMC = 0
                     ){
  
  
    //==================================================================================
    //DEFINING SOME STUFF
    //==================================================================================
  
    Bool_t kAllKappaPlots = 0;
    Bool_t kUseFitConstraints = 1;
    if(kMC){
      kAllKappaPlots = 1;
    }
    
    Double_t nSigmaLow  = -20.0;
    Double_t nSigmaHigh = 20.0;
    
    StyleSettingsThesis();  
    SetPlotStyle();
    
    TString InfoSystem = Form("%s-%s %% PbPb, #sqrt{s}=2.76TeV",CentralityLow.Data(),CentralityHigh.Data());
    TString InfoData = "Data LHC10h";
    TString InfoMC = "MC LHC13d2";
//     TString SigmaStarForm = "#Kappa = #frac{|#kappa^{+}|+|#kappa^{-}|}{2}";
    TString SigmaStarForm = "#Kappa = #frac{|#kappa^{+}|+|#kappa^{-}|}{2}+2(#kappa^{+}+#kappa^{-})";
    //TString SigmaStarForm = "#Kappa = #sqrt{#frac{#kappa^{+}^{2}+#kappa^{+}^{2}}{2}}";
    
//     TFile* fileData = new TFile("PhotonQA_0705415160_redefinedkappalarge_MC.root");
    TFile* fileMC = new TFile(Form("/home/mike/3_PbPb_dirg/1_data/170213_PbPb_v2_final/mc/GammaConvFlow_%i.root",Trainconfig));
    TList* listMC_1   = new TList();
    listMC_1     = (TList*)fileMC->Get(Form("GammaConvV1_%i_v2",Trainconfig));
    cout << listMC_1 << endl;
    TList* listMC_2   = new TList();
    listMC_2     = (TList*)listMC_1->FindObject(Form("Cut Number %s",Cutnumber.Data()));
    cout << listMC_2 << endl;
    TList* listMC_3   = new TList();
    listMC_3     = (TList*)listMC_2->FindObject(Form("%s ESD histograms",Cutnumber.Data()));
    cout << listMC_3 << endl;
    
    TFile* fileData  = new TFile(Form("/home/mike/3_PbPb_dirg/1_data/170213_PbPb_v2_final/data/GammaConvFlow_%i.root",Trainconfig));
    TList* listData_1   = new TList();
    listData_1     = (TList*)fileData->Get(Form("GammaConvV1_%i_v2",Trainconfig));
    
    TList* listData_2   = new TList();
    listData_2     = (TList*)listData_1->FindObject(Form("Cut Number %s",Cutnumber.Data()));
    
    TList* listData_3   = new TList();
    listData_3     = (TList*)listData_2->FindObject(Form("%s ESD histograms",Cutnumber.Data())); 

    TH2F* h2d00 = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp0_ee");
    TH2F* h2d01 = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp9_rem4");
    TH2F* h2d11 = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp1_pipi");
    TH2F* h2d13 = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp2_pie");
    
    TCanvas* c1 = new TCanvas("c1","",1200,800);
    gStyle->SetOptStat(0);
    c1->Divide(3,2,0.001,0.001);
    TCanvas* c2 = new TCanvas("c2","",800,800);
    c2->SetLeftMargin(0.08);
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.08);
    TCanvas* c8 = new TCanvas("c8","",800,800);
    c8->SetLeftMargin(0.08);
    c8->SetTopMargin(0.05);
    c8->SetBottomMargin(0.08);
    TCanvas* c10 = new TCanvas("c10","",800,800);
    c10->SetLeftMargin(0.08);
    c10->SetTopMargin(0.05);
    c10->SetBottomMargin(0.08);
    
    Double_t newBinsComb[19] = {0.0, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.3, 3.7, 4.1, 4.6, 5.4, 6.2, 7.0};
    TString newBinsName[19] = {"00", "09", "11", "13", "15", "17", "19", "21", "23", "25", "27", "30", "33", "37", "41", "46", "54", "62", "70"};
    Double_t centralbinvalues[19];
    Double_t binwidths[19];
    
    for(int i=0;i<18;i++){
      centralbinvalues[i]=newBinsComb[i]+(newBinsComb[i+1]-newBinsComb[i])/2;
      binwidths[i]=(newBinsComb[i+1]-newBinsComb[i])/2;
    }
    
    Double_t values00[18],values01[18],values11[18],values13[18];
    Double_t errors00[18],errors01[18],errors11[18],errors13[18];
    
    Double_t NEntriesData;
    Double_t frac00[18];
    Double_t frac01[18];
    Double_t frac11[18];
    Double_t frac13[18];
    
    Double_t mean00[18];
    Double_t mean01[18];
    Double_t mean11[18];
    Double_t mean13[18];
    
    Double_t h1d00Int;
    Double_t h1d01Int;
    Double_t h1d11Int;
    Double_t h1d13Int;
    Double_t PuritySignal1[18];
    Double_t PuritySignal2[18];
    Double_t PuritySignal3[18];
    Double_t PuritySignal4[18];
    Double_t Signalpurity[18];
    Double_t SignalpurityMC[18];
    Double_t SignalNA[18];
    Double_t SignalNB[18];
    Double_t SignalNC[18];
    Double_t SignalCA[18];
    Double_t SignalCB[18];
    Double_t SignalCC[18];
    Double_t Background1purity[18];
    Double_t Background2NA[18];
    Double_t Background2NB[18];
    Double_t Background2NC[18];    
    
    //==================================================================================
    //MAIN LOOP OVER P_{T}
    //==================================================================================
    

     for(int i=0; i<18; i++){
       
       
	//==================================================================================
	//TEMPLATE PLOTS
	//==================================================================================
   
	Double_t PtLow 	= newBinsComb[i];
	Double_t PtHigh 	= newBinsComb[i+1];
	
// 	Double_t nSigmaLow	= -20.0;
// 	Double_t nSigmaHigh	= 20.0;
	
	TLatex T1;
	T1.SetTextSize(0.06);
	T1.SetTextAlign(12);
	T1.SetNDC();
      
	h2d00->GetYaxis()->SetRangeUser(PtLow,PtHigh);
	h2d01->GetYaxis()->SetRangeUser(PtLow,PtHigh);
	h2d11->GetYaxis()->SetRangeUser(PtLow,PtHigh);
	h2d13->GetYaxis()->SetRangeUser(PtLow,PtHigh);

	TH1D* h1d00 = (TH1D*)h2d00->ProjectionX();
  TH1D* h1d00Ratio = (TH1D*)h1d00->Clone();
	TH1D* h1d01 = (TH1D*)h2d01->ProjectionX();
	TH1D* h1d11 = (TH1D*)h2d11->ProjectionX();
	TH1D* h1d13 = (TH1D*)h2d13->ProjectionX();
	
	h1d00->SetTitle("");
	h1d00->GetYaxis()->SetTitle("N  ");
	h1d00->GetXaxis()->SetTitle("#Kappa");
	h1d00->GetXaxis()->SetTitleSize(0.05);
	h1d00->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	h1d00->GetXaxis()->SetTitleOffset(0.8);
  h1d00->GetYaxis()->SetTitleOffset(1.4);
	h1d00->SetFillColor(12);
	h1d00->SetLineColor(12);
	h1d00->SetFillStyle(3004);
  
  h1d00Ratio->SetTitle("");
  h1d00Ratio->GetYaxis()->SetTitle("Ratio");
  h1d00Ratio->GetXaxis()->SetTitle("#Kappa");
  h1d00Ratio->GetXaxis()->SetTitleSize(0.05);
  h1d00Ratio->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
  h1d00Ratio->GetYaxis()->SetRangeUser(0,1.0);
  h1d00Ratio->GetXaxis()->SetTitleOffset(0.8);
  h1d00Ratio->SetFillColor(12);
  h1d00Ratio->SetLineColor(12);
  h1d00Ratio->SetFillStyle(3004);
	
	h1d01->SetTitle("");
	h1d01->GetYaxis()->SetTitle("counts");
	h1d01->GetXaxis()->SetTitle("#Kappa");
	h1d01->GetXaxis()->SetTitleSize(0.05);
	h1d01->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	h1d01->GetXaxis()->SetTitleOffset(0.8);
	h1d01->SetFillColor(kRed+1);
	h1d01->SetLineColor(kRed+1);
	h1d01->SetFillStyle(3005);
	
	h1d11->SetTitle("");
	h1d11->GetYaxis()->SetTitle("counts");
	h1d11->GetXaxis()->SetTitle("#Kappa");
	h1d11->GetXaxis()->SetTitleSize(0.05);
	h1d11->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	h1d11->GetXaxis()->SetTitleOffset(0.8);
	h1d11->SetFillColor(kBlue-7);
	h1d11->SetLineColor(kBlue-7);
	h1d11->SetFillStyle(3005);
	
	h1d13->SetTitle("");
	h1d13->GetYaxis()->SetTitle("counts");
	h1d13->GetXaxis()->SetTitle("#Kappa");
	h1d13->GetXaxis()->SetTitleSize(0.05);
	h1d13->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	h1d13->GetXaxis()->SetTitleOffset(0.8);
	h1d13->SetFillColor(kMagenta-2);
	h1d13->SetLineColor(kMagenta-2);
	h1d13->SetFillStyle(3004);
	
	TH1D *h1sum=(TH1D*)h1d00->Clone();
	h1sum->GetXaxis()->SetTitleOffset(0.8);
	h1sum->SetFillColor(0);
	h1sum->SetLineColor(kBlack);
	h1sum->SetLineWidth(2.0);
	h1sum->Add(h1d01);
	h1sum->Add(h1d11);
	h1sum->Add(h1d13);
  
  TH1D *h1sumbkg=(TH1D*)h1d01->Clone();
  h1sumbkg->GetXaxis()->SetTitleOffset(0.8);
  h1sumbkg->SetFillColor(0);
  h1sumbkg->SetLineColor(kBlack);
  h1sumbkg->SetLineWidth(2.0);
  h1sumbkg->Add(h1d11);
  h1sumbkg->Add(h1d13);
  
  TH1D *h1ratiobkg1=(TH1D*)h1d01->Clone();
  h1ratiobkg1->GetXaxis()->SetTitleOffset(0.8);
  h1ratiobkg1->GetYaxis()->SetTitle("Ratio");
  h1ratiobkg1->GetYaxis()->SetRangeUser(0,1.8);
  h1ratiobkg1->SetFillColor(0);
  h1ratiobkg1->SetLineColor(kRed+1);
  h1ratiobkg1->SetLineWidth(2.0);
  h1ratiobkg1->Divide(h1sumbkg);
  
  TH1D *h1ratiobkg2=(TH1D*)h1d11->Clone();
  h1ratiobkg2->GetXaxis()->SetTitleOffset(0.8);
  h1ratiobkg2->GetYaxis()->SetTitle("Ratio");
  h1ratiobkg2->GetYaxis()->SetRangeUser(0,1.8);
  h1ratiobkg2->SetFillColor(0);
  h1ratiobkg2->SetLineColor(kBlue-7);
  h1ratiobkg2->SetLineWidth(2.0);
  h1ratiobkg2->Divide(h1sumbkg);
  
  TH1D *h1ratiobkg3=(TH1D*)h1d13->Clone();
  h1ratiobkg3->GetXaxis()->SetTitleOffset(0.8);
  h1ratiobkg3->GetYaxis()->SetTitle("Ratio");
  h1ratiobkg3->GetYaxis()->SetRangeUser(0,1.8);
  h1ratiobkg3->SetFillColor(0);
  h1ratiobkg3->SetLineColor(kMagenta-2);
  h1ratiobkg3->SetLineWidth(2.0);
  h1ratiobkg3->Divide(h1sumbkg);
	
	TH2F* h2d_data;
	if(kMC){
	  h2d_data = (TH2F*)listMC_3->FindObject("KappaTPC_Pt_after");
	}else{
	  h2d_data = (TH2F*)listData_3->FindObject("KappaTPC_Pt_after");
	}
	h2d_data->GetYaxis()->SetRangeUser(PtLow,PtHigh);
	TH1D* data = (TH1D*)h2d_data->ProjectionX();
	data->SetTitle("");
	data->GetXaxis()->SetTitle("#Kappa");
	data->GetYaxis()->SetTitle("counts");
	data->GetXaxis()->SetTitleSize(0.05);
	data->GetXaxis()->SetTitleOffset(0.8);
	data->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	
	Double_t NTotal = h1d00->GetEntries()+h1d01->GetEntries()+h1d11->GetEntries()+h1d13->GetEntries();

	c1->cd(1);
	h1d00->Draw("HIST");
	T1.DrawLatex(0.65, 0.9, "signal");
	T1.DrawLatex(0.65, 0.85, Form("%f",h1d00->GetEntries()/NTotal));
	T1.SetTextSize(0.03);
	T1.DrawLatex(0.18, 0.90, InfoSystem.Data());
	T1.DrawLatex(0.18, 0.85, InfoMC.Data());
	T1.DrawLatex(0.18, 0.75, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
	T1.DrawLatex(0.18, 0.65,  SigmaStarForm.Data());
	
	c1->cd(4);
	h1d01->Draw("HIST");
	T1.SetTextSize(0.06);
	T1.DrawLatex(0.65, 0.9, "remaining");
	T1.DrawLatex(0.65, 0.85, Form("%f",h1d01->GetEntries()/NTotal));
	
	c1->cd(3);
	h1d11->Draw("HIST");
	T1.DrawLatex(0.65, 0.9, "#pi^{#pm} + #pi^{#mp}");
	T1.DrawLatex(0.65, 0.85, Form("%f",h1d11->GetEntries()/NTotal));
	
	c1->cd(2);
	h1d13->Draw("HIST");
	T1.DrawLatex(0.65, 0.9, "#pi^{#pm} + e^{#mp}");
	T1.DrawLatex(0.65, 0.85, Form("%f",h1d13->GetEntries()/NTotal));

	c1->cd(5);
	h1sum->Draw("HIST");
	h1d00->Draw("HIST same");
	h1d01->Draw("HIST same");
	h1d11->Draw("HIST same");
	h1d13->Draw("HIST same");
	
	TLegend* leg = new TLegend(0.7,0.6,0.9,0.9);
	leg->SetTextSize(0.04);
	leg->AddEntry(h1sum,"total","lp");
	leg->AddEntry(h1d00,"real e^{+}e^{-}","f");
	leg->AddEntry(h1d01,"remaining","f");
	leg->AddEntry(h1d11,"#pi^{#pm} + #pi^{#mp}","f");
	leg->AddEntry(h1d13,"#pi^{#pm} + e^{#mp}","f");
	leg->SetBorderSize(0);
	leg->Draw();
  
  TLegend* leg5 = new TLegend(0.7,0.7,0.9,0.9);
  leg5->SetTextSize(0.04);
//   leg5->SetHeader("Ratio to sum");
  leg5->AddEntry(h1ratiobkg1,"remaining","l");
  leg5->AddEntry(h1ratiobkg2,"#pi^{#pm} + #pi^{#mp}","l");
  leg5->AddEntry(h1ratiobkg3,"#pi^{#pm} + e^{#mp}","l");
  leg5->SetBorderSize(0);
//   leg5->Draw();
	
	c2->cd();
	h1sum->Draw("HIST");
	h1d00->Draw("HIST same");
	h1d01->Draw("HIST same");
	h1d11->Draw("HIST same");
	h1d13->Draw("HIST same");
  
	leg->Draw();
	
	T1.SetTextSize(0.03);
	T1.DrawLatex(0.15, 0.97, InfoSystem.Data());
	T1.DrawLatex(0.65, 0.97, InfoMC.Data());
	T1.DrawLatex(0.18, 0.75, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
	T1.DrawLatex(0.18, 0.65,  SigmaStarForm.Data());
//   T1.DrawLatex(0.65, 0.4,  Form("p(-5.0<K<10.0)= %2.3f",Signalpurity[i]));
  
  c8->cd();
  h1ratiobkg1->Draw("HIST");
  h1ratiobkg2->Draw("HIST same");
  h1ratiobkg3->Draw("HIST same");
  
  leg5->Draw();
  
  T1.SetTextSize(0.03);
  T1.DrawLatex(0.15, 0.97, InfoSystem.Data());
  T1.DrawLatex(0.65, 0.97, InfoMC.Data());
  T1.DrawLatex(0.15, 0.8, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
  T1.DrawLatex(0.15, 0.7,  SigmaStarForm.Data());
  T1.DrawLatex(0.15, 0.9, Form("N_{bkg,i} / N_{bkg,tot}"));
  
  c10->cd();
  h1d00Ratio->Divide(h1sum);
  h1d00Ratio->Draw("HIST");
  
  T1.SetTextSize(0.03);
  T1.DrawLatex(0.15, 0.97, InfoSystem.Data());
  T1.DrawLatex(0.65, 0.97, InfoMC.Data());
  T1.DrawLatex(0.15, 0.8, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
  T1.DrawLatex(0.15, 0.7,  SigmaStarForm.Data());
  T1.DrawLatex(0.15, 0.9, Form("N_{signal}/N_{total}"));
	
	
	//==================================================================================
	//TFRACTIONFITTER
	//==================================================================================
	
	c1->cd(6);
	data->Draw("HIST");
	
	Double_t NTotalMC = h1d00->GetEntries()+h1d01->GetEntries()+h1d11->GetEntries()+h1d13->GetEntries();
	Double_t constraint_low;
	Double_t constraint_high;
	Double_t fracwidth;
	
	if(kMC){
	  NEntriesData = data->GetEntries();
	  frac00[i] = h1d00->GetEntries()/NEntriesData;
	  frac01[i] = h1d01->GetEntries()/NEntriesData;
	  frac11[i] = h1d11->GetEntries()/NEntriesData;
	  frac13[i] = h1d13->GetEntries()/NEntriesData;
	  fracwidth = 1.2; //0-20: 1.3, 20-40:1.2, 40-80:2.0
	  constraint_low  = 1/fracwidth;
	  constraint_high = 1*fracwidth;
	}else{
	  NEntriesData = data->GetEntries();
	  frac00[i] = h1d00->GetEntries()/NTotalMC;
	  frac01[i] = h1d01->GetEntries()/NTotalMC;
	  frac11[i] = h1d11->GetEntries()/NTotalMC;
	  frac13[i] = h1d13->GetEntries()/NTotalMC;
	  fracwidth = 1.2;
	  constraint_low  = 1/fracwidth;
	  constraint_high = 1*fracwidth;
	}
	
	//array of all MC that build up the data
	TObjArray *mc = new TObjArray(4);
	mc->Add(h1d00);	//primary + secondary
	mc->Add(h1d01);	//remaining
	mc->Add(h1d11); //pion+pion
	mc->Add(h1d13); //pion+electron
	
	//configure the TFractionFitter
	TFractionFitter* fit = new TFractionFitter(data, mc);
	if(kUseFitConstraints){
	  fit->Constrain(1,frac00[i]*constraint_low,frac00[i]*constraint_high);
	  fit->Constrain(2,frac01[i]*constraint_low,frac01[i]*constraint_high);
	  fit->Constrain(3,frac11[i]*constraint_low,frac11[i]*constraint_high);
	  fit->Constrain(4,frac13[i]*constraint_low,frac13[i]*constraint_high);
	}
	
	TH1D* template1 = NULL;
	TH1D* template2 = NULL;
	TH1D* template3 = NULL;
	TH1D* template4 = NULL;
  TH1D* templatesum = NULL;
  TH1D* hsumMC = NULL;
	//Fit the templates
	Int_t status = fit->Fit();
	//cout << "fit status: " << status << endl;
	if (status == 0) {
	  TH1D* result = (TH1D*)fit->GetPlot();
	  result->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  result->SetLineColor(kRed);
	  result->Draw("HIST same");
	  
	  Double_t value[4],error[4];
	  for(int j=0;j<4;j++){
	    fit->GetResult(j,value[j],error[j]);
	  }
	  
	  values00[i]=value[0];
    errors00[i]=error[0];
    
    cout << " value[0] = " << value[0] << endl;
    
    values01[i]=value[1];
    errors01[i]=error[1];
    
    cout << " value[1] = " << value[1] << endl;
    
    values11[i]=value[2];
    errors11[i]=error[2];
    
    cout << " value[2] = " << value[2] << endl;
    
    values13[i]=value[3];
    errors13[i]=error[3];
    
    cout << " value[3] = " << value[3] << endl;
    
    cout << " frac00[" << i << "] = " << frac00[i] << endl;
    cout << " frac01[" << i << "] = " << frac01[i] << endl;
    cout << " frac11[" << i << "] = " << frac11[i] << endl;
    cout << " frac13[" << i << "] = " << frac13[i] << endl;
    
    cout << "NEntriesData = " << NEntriesData << endl;
    cout << "NTotalMC = " << NTotalMC << endl;
    
    template1 = (TH1D*)h1d00->Clone();
    template1->Scale(value[0]/frac00[i]*(NEntriesData/NTotalMC));
//     template1->Scale((value[0]-frac00[i])/frac00[i]+1);
//     template1->Scale((NEntriesData/NTotalMC));
//     template1->Draw("HIST same");
    
    template2 = (TH1D*)h1d01->Clone();
    template2->Scale(value[1]/frac01[i]*(NEntriesData/NTotalMC));
//     template2->Scale((value[1]-frac01[i])/frac01[i]+1);
//     template2->Scale((NEntriesData/NTotalMC));
//     template2->Draw("HIST same");
    
    template3= (TH1D*)h1d11->Clone();
    template3->Scale(value[2]/frac11[i]*(NEntriesData/NTotalMC));
//     template3->Scale((value[2]-frac11[i])/frac11[i]+1);
//     template3->Scale((NEntriesData/NTotalMC));
//     template3->Draw("HIST same");
    
    template4 = (TH1D*)h1d13->Clone();
    template4->Scale(value[3]/frac13[i]*(NEntriesData/NTotalMC));
//     template4->Scale((value[3]-frac13[i])/frac13[i]+1);
//     template4->Scale((NEntriesData/NTotalMC));
//     template4->Draw("HIST same");
    
    TLegend* leg2 = new TLegend(0.7,0.7,0.9,0.9);
    leg2->AddEntry(data,"data","lp");
    leg2->AddEntry(result,"fit","lp");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->Draw();
    
    mean00[i] = h1d00->GetMean();
    mean01[i] = h1d01->GetMean();
    mean11[i] = h1d11->GetMean();
    mean13[i] = h1d13->GetMean();
    
    templatesum=(TH1D*)template1->Clone();
    templatesum->GetXaxis()->SetTitleOffset(0.8);
    templatesum->SetFillColor(0);
    templatesum->SetLineColor(kBlack);
    templatesum->SetLineWidth(2.0);
    templatesum->Add(template2);
    templatesum->Add(template3);
    templatesum->Add(template4);
    
    hsumMC=(TH1D*)h1d00->Clone();
    hsumMC->Add(h1d01);
    hsumMC->Add(h1d11);
    hsumMC->Add(h1d13);
    
    Int_t IntLowerBound1 = -20.0;
    Int_t IntLowerBound2 = -11.0;
    Int_t IntLowerBound3 = -3.0;
    Int_t IntLowerBound4 = 11.0;
    Int_t IntUpperBound1 = -13.0;
    Int_t IntUpperBound2 = -6.0;
    Int_t IntUpperBound3 = 5.0;
    Int_t IntUpperBound4 = 20.0;
    
    Signalpurity[i] = template1->Integral(template1->GetXaxis()->FindBin(IntLowerBound3),template1->GetXaxis()->FindBin(IntUpperBound3)) / templatesum->Integral(templatesum->GetXaxis()->FindBin(IntLowerBound3),templatesum->GetXaxis()->FindBin(IntUpperBound3));
    SignalpurityMC[i] = h1d00->Integral(h1d00->GetXaxis()->FindBin(IntLowerBound3),h1d00->GetXaxis()->FindBin(IntUpperBound3)) / hsumMC->Integral(hsumMC->GetXaxis()->FindBin(IntLowerBound3),hsumMC->GetXaxis()->FindBin(IntUpperBound3));
    SignalNA[i] = template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound3),template3->GetXaxis()->FindBin(IntUpperBound3)) / (template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound3),template2->GetXaxis()->FindBin(IntUpperBound3))+template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound3),template3->GetXaxis()->FindBin(IntUpperBound3))+template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound3),template4->GetXaxis()->FindBin(IntUpperBound3)));
    SignalNB[i] = template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound3),template4->GetXaxis()->FindBin(IntUpperBound3)) / (template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound3),template2->GetXaxis()->FindBin(IntUpperBound3))+template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound3),template3->GetXaxis()->FindBin(IntUpperBound3))+template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound3),template4->GetXaxis()->FindBin(IntUpperBound3)));
    SignalNC[i] = template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound3),template2->GetXaxis()->FindBin(IntUpperBound3)) / (template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound3),template2->GetXaxis()->FindBin(IntUpperBound3))+template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound3),template3->GetXaxis()->FindBin(IntUpperBound3))+template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound3),template4->GetXaxis()->FindBin(IntUpperBound3)));
    SignalCA[i] = template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound3),template3->GetXaxis()->FindBin(IntUpperBound3)) / (template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound3),template2->GetXaxis()->FindBin(IntUpperBound3))+template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound3),template3->GetXaxis()->FindBin(IntUpperBound3))+template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound3),template4->GetXaxis()->FindBin(IntUpperBound3))+template1->Integral(template1->GetXaxis()->FindBin(IntLowerBound3),template1->GetXaxis()->FindBin(IntUpperBound3)));
    SignalCB[i] = template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound3),template4->GetXaxis()->FindBin(IntUpperBound3)) / (template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound3),template2->GetXaxis()->FindBin(IntUpperBound3))+template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound3),template3->GetXaxis()->FindBin(IntUpperBound3))+template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound3),template4->GetXaxis()->FindBin(IntUpperBound3))+template1->Integral(template1->GetXaxis()->FindBin(IntLowerBound3),template1->GetXaxis()->FindBin(IntUpperBound3)));
    SignalCC[i] = template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound3),template2->GetXaxis()->FindBin(IntUpperBound3)) / (template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound3),template2->GetXaxis()->FindBin(IntUpperBound3))+template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound3),template3->GetXaxis()->FindBin(IntUpperBound3))+template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound3),template4->GetXaxis()->FindBin(IntUpperBound3))+template1->Integral(template1->GetXaxis()->FindBin(IntLowerBound3),template1->GetXaxis()->FindBin(IntUpperBound3)));
    Background1purity[i] = template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound1),template3->GetXaxis()->FindBin(IntUpperBound1)) / (template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound1),template2->GetXaxis()->FindBin(IntUpperBound1))+template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound1),template3->GetXaxis()->FindBin(IntUpperBound1)));
    Background2NA[i] = template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound2),template3->GetXaxis()->FindBin(IntUpperBound2)) / (template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound2),template2->GetXaxis()->FindBin(IntUpperBound2))+template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound2),template3->GetXaxis()->FindBin(IntUpperBound2))+template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound2),template4->GetXaxis()->FindBin(IntUpperBound2)));
    Background2NB[i] = template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound2),template4->GetXaxis()->FindBin(IntUpperBound2)) / (template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound2),template2->GetXaxis()->FindBin(IntUpperBound2))+template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound2),template3->GetXaxis()->FindBin(IntUpperBound2))+template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound2),template4->GetXaxis()->FindBin(IntUpperBound2)));
    Background2NC[i] = template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound2),template2->GetXaxis()->FindBin(IntUpperBound2)) / (template2->Integral(template2->GetXaxis()->FindBin(IntLowerBound2),template2->GetXaxis()->FindBin(IntUpperBound2))+template3->Integral(template3->GetXaxis()->FindBin(IntLowerBound2),template3->GetXaxis()->FindBin(IntUpperBound2))+template4->Integral(template4->GetXaxis()->FindBin(IntLowerBound2),template4->GetXaxis()->FindBin(IntUpperBound2)));
    
  }
  
  cout << " " << endl;
  cout << Signalpurity[i] << endl;
  cout << " " << endl;
  
  TCanvas* c13 = new TCanvas("c12","",800,800);
//   c13->SetLeftMargin(0.08);
  c13->SetTopMargin(0.05);
//   c13->SetBottomMargin(0.08);
  if (status == 0) {
  Double_t maximum = template1->GetMaximum();
  if(template2->GetMaximum()>maximum) maximum = template2->GetMaximum();
  if(template3->GetMaximum()>maximum) maximum = template3->GetMaximum();
  if(template3->GetMaximum()>maximum) maximum = template4->GetMaximum();
  if(templatesum->GetMaximum()>maximum) maximum = templatesum->GetMaximum();
  template1->GetYaxis()->SetRangeUser(0, maximum*1.2);
  templatesum->Draw("HIST same");
  template1->Draw("HIST same");
  template2->Draw("HIST same");
  template3->Draw("HIST same");
  template4->Draw("HIST same");
  }
  
  leg->Draw();
  
  T1.SetTextSize(0.03);
  T1.DrawLatex(0.2, 0.97, InfoSystem.Data());
  if(kMC){
    T1.DrawLatex(0.65, 0.97, InfoMC.Data());
  }else{
    T1.DrawLatex(0.65, 0.97, InfoData.Data());
  }
  T1.DrawLatex(0.18, 0.75, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
	T1.DrawLatex(0.18, 0.65,  SigmaStarForm.Data());
//   T1.DrawLatex(0.65, 0.4,  Form("p(-5.0<K<10.0)= %2.3f",Signalpurity[i]));
  
	gSystem->mkdir(Form("purity_studies_%s",Cutnumber.Data()));
	if(kMC){
	  c1->SaveAs(Form("purity_studies_%s/MC_Kappa_4Templates_pt_%s-%s.eps",Cutnumber.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
    c13->SaveAs(Form("purity_studies_%s/MC_Kappa_4Templates_TFraction_pt_%s-%s.eps",Cutnumber.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
 	  c2->SaveAs(Form("purity_studies_%s/MC_Kappa_4Templates_same_pt_%s-%s.eps",Cutnumber.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
    c8->SaveAs(Form("purity_studies_%s/MC_Kappa_4Templates_ratiobkg_pt_%s-%s.eps",Cutnumber.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
    c10->SaveAs(Form("purity_studies_%s/MC_Kappa_4Templates_ratiosig_pt_%s-%s.eps",Cutnumber.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
	}else{
	  c1->SaveAs(Form("purity_studies_%s/Data_Kappa_4Templates_pt_%s-%s.eps",Cutnumber.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
    c13->SaveAs(Form("purity_studies_%s/Data_Kappa_4Templates_TFraction_pt_%s-%s.eps",Cutnumber.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
	}
	

    }
    //end of pt loop
    
    
    //==================================================================================
    //FRACTION PLOTS
    //==================================================================================

    Double_t ratiofitfraction00[18];
    Double_t ratiofitfraction01[18];
    Double_t ratiofitfraction11[18];
    Double_t ratiofitfraction13[18];
    
    for(int i=0;i<18;i++){
	ratiofitfraction00[i]=values00[i]/frac00[i];
	ratiofitfraction01[i]=values01[i]/frac01[i];
	ratiofitfraction11[i]=values11[i]/frac11[i];
	ratiofitfraction13[i]=values13[i]/frac13[i];
    }

    TCanvas* c3 = new TCanvas("c3","",1350,450);
    c3->Divide(3,1,0.001,0.001);
    c3->cd(1);
    TGraphErrors* gr0 = new TGraphErrors(18,centralbinvalues,values00,binwidths,errors00);
    gr0->RemovePoint(0);
    gr0->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    gr0->GetXaxis()->SetRangeUser(0.0,7);
    gr0->GetXaxis()->SetTitleOffset(1.1);
    gr0->GetYaxis()->SetTitle("fraction");
    gr0->GetYaxis()->SetRangeUser(0,1.49);
    gr0->SetTitle("");
    gr0->SetMarkerColor(12);
    gr0->SetMarkerStyle(20);
    gr0->SetLineColor(12);
    gr0->Draw("AP");
    
    TGraphErrors* gr1 = new TGraphErrors(18,centralbinvalues,values01,binwidths,errors01);
    gr1->RemovePoint(0);
    gr1->SetMarkerColor(kRed+1);
    gr1->SetMarkerStyle(24);
    gr1->SetLineColor(kRed+2);
    gr1->Draw("P");
    
    TGraphErrors* gr2 = new TGraphErrors(18,centralbinvalues,values11,binwidths,errors11);
    gr2->RemovePoint(0);
    gr2->SetMarkerColor(kBlue-7);
    gr2->SetMarkerStyle(24);
    gr2->SetLineColor(kBlue-7);
    gr2->Draw("P");
    
    TGraphErrors* gr3 = new TGraphErrors(18,centralbinvalues,values13,binwidths,errors13);
    gr3->RemovePoint(0);
    gr3->SetMarkerColor(kMagenta-2);
    gr3->SetMarkerStyle(24);
    gr3->SetLineColor(kMagenta-2);
    gr3->Draw("P");
    
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(gr0,"real e^{+}e^{-}","lp");
    leg->AddEntry(gr1,"remaining","lp");
    leg->AddEntry(gr2,"#pi^{#pm} + #pi^{#mp}","lp");
    leg->AddEntry(gr3,"#pi^{#pm} + e^{#mp}","lp");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->Draw();
    
    TLatex T1;
    T1.SetTextSize(0.04);
    T1.SetTextAlign(12);
    T1.SetNDC();
    T1.DrawLatex(0.2, 0.9, "Template fractions");
    T1.SetTextSize(0.03);
    T1.DrawLatex(0.2, 0.8, InfoSystem.Data());
    if(kMC){
      T1.DrawLatex(0.2, 0.75, InfoMC.Data());
    }
    else{
      T1.DrawLatex(0.2, 0.75, InfoData.Data());
    }
    
    c3->cd(2);
    TGraphErrors* gr8 = new TGraphErrors(18,centralbinvalues,ratiofitfraction00,binwidths,0);
    gr8->RemovePoint(0);
    gr8->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    gr8->GetXaxis()->SetRangeUser(0.0,7);
    gr8->GetXaxis()->SetTitleOffset(1.1);
    gr8->GetYaxis()->SetTitle("fit/MC");
    if(kMC){
      gr8->GetYaxis()->SetRangeUser(0.9,1.1);
    }else{
      gr8->GetYaxis()->SetRangeUser(0.0,3.0);
    }
    gr8->SetTitle("");
    gr8->SetMarkerColor(12);
    gr8->SetMarkerStyle(20);
    gr8->SetLineColor(12);
    gr8->Draw("AP");
    
    TGraphErrors* gr9 = new TGraphErrors(18,centralbinvalues,ratiofitfraction01,binwidths,0);
    gr9->RemovePoint(0);
    gr9->SetMarkerColor(kRed+1);
    gr9->SetMarkerStyle(24);
    gr9->SetLineColor(kRed+1);
    gr9->Draw("P");
    
    TGraphErrors* gr10 = new TGraphErrors(18,centralbinvalues,ratiofitfraction11,binwidths,0);
    gr10->RemovePoint(0);
    gr10->SetMarkerColor(kBlue-7);
    gr10->SetMarkerStyle(24);
    gr10->SetLineColor(kBlue-7);
    gr10->Draw("P");
    
    TGraphErrors* gr11 = new TGraphErrors(18,centralbinvalues,ratiofitfraction13,binwidths,0);
    gr11->RemovePoint(0);
    gr11->SetMarkerColor(kMagenta-2);
    gr11->SetMarkerStyle(24);
    gr11->SetLineColor(kMagenta-2);
    gr11->Draw("P");
    
    
    TLegend* leg2 = new TLegend(0.7,0.7,0.9,0.9);
    leg2->AddEntry(gr8,"real e^{+}e^{-}","lp");
    leg2->AddEntry(gr9,"remaining","lp");
    leg2->AddEntry(gr10,"#pi^{#pm} + #pi^{#mp}","lp");
    leg2->AddEntry(gr11,"#pi^{#pm} + e^{#mp}","lp");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->Draw();
    
    T1.SetTextSize(0.04);
    T1.DrawLatex(0.2, 0.9, "Ratio Data/MC");
    
    c3->cd(3);
    TGraphErrors* graphPurity = new TGraphErrors(18,centralbinvalues,Signalpurity,0,0);
    graphPurity->RemovePoint(0);
    graphPurity->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphPurity->GetXaxis()->SetRangeUser(0.0,7);
    graphPurity->GetXaxis()->SetTitleOffset(1.1);
    graphPurity->GetYaxis()->SetTitle("Purity");
    graphPurity->GetYaxis()->SetRangeUser(0.5,1.0);
    graphPurity->GetYaxis()->SetTitleOffset(1.15);
    graphPurity->SetTitle("");
    graphPurity->SetMarkerColor(12);
    graphPurity->SetMarkerStyle(20);
    graphPurity->SetLineColor(12);
    graphPurity->Draw("AP");
    
    TGraphErrors* graphPurityMC = new TGraphErrors(18,centralbinvalues,SignalpurityMC,0,0);
    graphPurityMC->RemovePoint(0);
    graphPurityMC->SetMarkerColor(9);
    graphPurityMC->SetMarkerStyle(21);
    graphPurityMC->SetLineColor(9);
    graphPurityMC->Draw("AP");
    
    TGraphErrors* graphSignalNA = new TGraphErrors(18,centralbinvalues,SignalNA,0,0);
    graphSignalNA->RemovePoint(0);
    graphSignalNA->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphSignalNA->GetXaxis()->SetRangeUser(0.0,7);
    graphSignalNA->GetXaxis()->SetTitleOffset(1.1);
    graphSignalNA->GetYaxis()->SetTitle("Purity");
//     graphSignalNA->GetYaxis()->SetRangeUser(0.0,1.0);
    graphSignalNA->SetTitle("");
    graphSignalNA->SetMarkerColor(kBlue+2);
    graphSignalNA->SetMarkerStyle(20);
    graphSignalNA->SetLineColor(12);
    graphSignalNA->Draw("SAME P");
    
    TGraphErrors* graphSignalNB = new TGraphErrors(18,centralbinvalues,SignalNB,0,0);
    graphSignalNB->RemovePoint(0);
    graphSignalNB->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphSignalNB->GetXaxis()->SetRangeUser(0.0,7);
    graphSignalNB->GetXaxis()->SetTitleOffset(1.1);
    graphSignalNB->GetYaxis()->SetTitle("Purity");
//     graphSignalNB->GetYaxis()->SetRangeUser(0.0,1.0);
    graphSignalNB->SetTitle("");
    graphSignalNB->SetMarkerColor(kMagenta+2);
    graphSignalNB->SetMarkerStyle(20);
    graphSignalNB->SetLineColor(12);
    graphSignalNB->Draw("SAME P");
    
    TGraphErrors* graphSignalNC = new TGraphErrors(18,centralbinvalues,SignalNC,0,0);
    graphSignalNC->RemovePoint(0);
    graphSignalNC->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphSignalNC->GetXaxis()->SetRangeUser(0.0,7);
    graphSignalNC->GetXaxis()->SetTitleOffset(1.1);
    graphSignalNC->GetYaxis()->SetTitle("Purity");
//     graphSignalNC->GetYaxis()->SetRangeUser(0.0,1.0);
    graphSignalNC->SetTitle("");
    graphSignalNC->SetMarkerColor(kRed+2);
    graphSignalNC->SetMarkerStyle(20);
    graphSignalNC->SetLineColor(12);
    graphSignalNC->Draw("SAME P");
    
    TGraphErrors* graphBackground1purity = new TGraphErrors(18,centralbinvalues,Background1purity,0,0);
    graphBackground1purity->RemovePoint(0);
    graphBackground1purity->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphBackground1purity->GetXaxis()->SetRangeUser(0.0,7);
    graphBackground1purity->GetXaxis()->SetTitleOffset(1.1);
    graphBackground1purity->GetYaxis()->SetTitle("purity");
//     graphBackground1purity->GetYaxis()->SetRangeUser(0.0,1.0);
    graphBackground1purity->SetTitle("");
    graphBackground1purity->SetMarkerColor(12);
    graphBackground1purity->SetMarkerStyle(20);
    graphBackground1purity->SetLineColor(12);
    
    TGraphErrors* graphBackground2NA = new TGraphErrors(18,centralbinvalues,Background2NA,0,0);
    graphBackground2NA->RemovePoint(0);
    graphBackground2NA->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphBackground2NA->GetXaxis()->SetRangeUser(0.0,7);
    graphBackground2NA->GetXaxis()->SetTitleOffset(1.1);
    graphBackground2NA->GetYaxis()->SetTitle("N_{A}");
//     graphBackground2NA->GetYaxis()->SetRangeUser(0.0,1.0);
    graphBackground2NA->SetTitle("");
    graphBackground2NA->SetMarkerColor(12);
    graphBackground2NA->SetMarkerStyle(20);
    graphBackground2NA->SetLineColor(12);
    
    TGraphErrors* graphBackground2NB = new TGraphErrors(18,centralbinvalues,Background2NB,0,0);
    graphBackground2NB->RemovePoint(0);
    graphBackground2NB->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphBackground2NB->GetXaxis()->SetRangeUser(0.0,7);
    graphBackground2NB->GetXaxis()->SetTitleOffset(1.1);
    graphBackground2NB->GetYaxis()->SetTitle("N_{B}");
//     graphBackground2NB->GetYaxis()->SetRangeUser(0.0,1.0);
    graphBackground2NB->SetTitle("");
    graphBackground2NB->SetMarkerColor(12);
    graphBackground2NB->SetMarkerStyle(20);
    graphBackground2NB->SetLineColor(12);
    
    TGraphErrors* graphBackground2NC = new TGraphErrors(18,centralbinvalues,Background2NC,0,0);
    graphBackground2NC->RemovePoint(0);
    graphBackground2NC->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphBackground2NC->GetXaxis()->SetRangeUser(0.0,7);
    graphBackground2NC->GetXaxis()->SetTitleOffset(1.1);
    graphBackground2NC->GetYaxis()->SetTitle("N_{C}");
//     graphBackground2NC->GetYaxis()->SetRangeUser(0.0,1.0);
    graphBackground2NC->SetTitle("");
    graphBackground2NC->SetMarkerColor(12);
    graphBackground2NC->SetMarkerStyle(20);
    graphBackground2NC->SetLineColor(12);
    
    TGraphErrors* graphSignalCA = new TGraphErrors(18,centralbinvalues,SignalCA,0,0);
    graphSignalCA->RemovePoint(0);
    graphSignalCA->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphSignalCA->GetXaxis()->SetRangeUser(0.0,7);
    graphSignalCA->GetXaxis()->SetTitleOffset(1.1);
    graphSignalCA->GetYaxis()->SetTitle("c");
    graphSignalCA->GetYaxis()->SetRangeUser(0.0,0.19);
    graphSignalCA->SetTitle("");
    graphSignalCA->SetMarkerColor(kBlue+2);
    graphSignalCA->SetMarkerStyle(20);
    graphSignalCA->SetLineColor(12);
//     graphSignalCA->Draw("SAME P");
    
    TGraphErrors* graphSignalCB = new TGraphErrors(18,centralbinvalues,SignalCB,0,0);
    graphSignalCB->RemovePoint(0);
    graphSignalCB->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphSignalCB->GetXaxis()->SetRangeUser(0.0,7);
    graphSignalCB->GetXaxis()->SetTitleOffset(1.1);
    graphSignalCB->GetYaxis()->SetTitle("Purity");
//     graphSignalCB->GetYaxis()->SetRangeUser(0.0,1.0);
    graphSignalCB->SetTitle("");
    graphSignalCB->SetMarkerColor(kMagenta+2);
    graphSignalCB->SetMarkerStyle(21);
    graphSignalCB->SetLineColor(12);
//     graphSignalCB->Draw("SAME P");
    
    TGraphErrors* graphSignalCC = new TGraphErrors(18,centralbinvalues,SignalCC,0,0);
    graphSignalCC->RemovePoint(0);
    graphSignalCC->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    graphSignalCC->GetXaxis()->SetRangeUser(0.0,7);
    graphSignalCC->GetXaxis()->SetTitleOffset(1.1);
    graphSignalCC->GetYaxis()->SetTitle("Purity");
//     graphSignalCC->GetYaxis()->SetRangeUser(0.0,1.0);
    graphSignalCC->SetTitle("");
    graphSignalCC->SetMarkerColor(kRed+2);
    graphSignalCC->SetMarkerStyle(22);
    graphSignalCC->SetLineColor(12);
//     graphSignalCC->Draw("SAME P");

    
    TLegend* leg4 = new TLegend(0.65,0.3,0.85,0.5);
    leg4->AddEntry(graphPurity,"Data driven","lp");
    leg4->AddEntry(graphPurityMC,"MC driven","lp");
    leg4->SetBorderSize(0);
    leg4->SetTextSize(0.04);
//     leg4->Draw();
    
    TLegend* leg5 = new TLegend(0.65,0.3,0.85,0.5);
    leg5->AddEntry(graphSignalCA,"c_{#pi#pi}","lp");
    leg5->AddEntry(graphSignalCB,"c_{#pie}","lp");
    leg5->AddEntry(graphSignalCC,"c_{rem}","lp");
    leg5->SetBorderSize(0);
    leg5->SetTextSize(0.04);
//     leg5->Draw();
    
    T1.DrawLatex(0.2, 0.15, "Purity #gamma sample (-3<K<5)");
    
    if(kMC){
      c3->SaveAs(Form("purity_studies_%s/MC_fractionfitter_4Templates.eps",Cutnumber.Data()));
    }else{
      c3->SaveAs(Form("purity_studies_%s/Data_fractionfitter_4Templates.eps",Cutnumber.Data()));
    }
    
    TCanvas* c11 = new TCanvas("c4","",800,800);
    c11->SetLeftMargin(0.08);
    c11->SetTopMargin(0.05);
    c11->SetBottomMargin(0.08);
    
    graphPurity->Draw("AP");
    graphPurityMC->Draw("SAME P");
    
    T1.DrawLatex(0.15, 0.97, InfoSystem.Data());
    if(kMC){
      T1.DrawLatex(0.65, 0.97, InfoMC.Data());
    }else{
      T1.DrawLatex(0.65, 0.97, InfoData.Data());
    }
    
    leg4->Draw();
    
    TCanvas* c12 = new TCanvas("c12","",800,800);
    c12->SetLeftMargin(0.08);
    c12->SetTopMargin(0.05);
    c12->SetBottomMargin(0.08);
    
    graphSignalCA->Draw("AP");
    graphSignalCB->Draw("SAME P");
    graphSignalCC->Draw("SAME P");
    
    leg5->Draw();
    
    T1.DrawLatex(0.15, 0.97, InfoSystem.Data());
    if(kMC){
      T1.DrawLatex(0.65, 0.97, InfoMC.Data());
    }else{
      T1.DrawLatex(0.65, 0.97, InfoData.Data());
    }
    T1.DrawLatex(0.65, 0.55, "#gamma_{incl} (-3<K<5)");
    
    if(kMC){
      c11->SaveAs(Form("purity_studies_%s/MC_Purity_4Templates.eps",Cutnumber.Data()));
      c12->SaveAs(Form("purity_studies_%s/MC_Contamination_4Templates.eps",Cutnumber.Data()));
    }else{
      c11->SaveAs(Form("purity_studies_%s/Data_Purity_4Templates.eps",Cutnumber.Data()));
      c12->SaveAs(Form("purity_studies_%s/Data_Contamination_4Templates.eps",Cutnumber.Data()));
    }
	
	TFile *Purity_file;
  if(kMC){
    Purity_file = new TFile(Form("purity_studies_%s/Purity_InclusivePhotonSample_%s%s_%s_%s.root",Cutnumber.Data(),CentralityLow.Data(),CentralityHigh.Data(),Cutnumber.Data(), InfoMC.Data()),"RECREATE");
  }else{
    Purity_file = new TFile(Form("purity_studies_%s/Purity_InclusivePhotonSample_%s%s_%s_%s.root",Cutnumber.Data(),CentralityLow.Data(),CentralityHigh.Data(),Cutnumber.Data(), InfoData.Data()),"RECREATE");
  }
  graphPurity->Write(Form("Gamma_Purity_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  graphSignalNA->Write(Form("Gamma_NA_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  graphSignalNB->Write(Form("Gamma_NB_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  graphSignalNC->Write(Form("Gamma_NC_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  graphBackground1purity->Write(Form("bkg1_Purity_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  graphBackground2NA->Write(Form("bkg2_NA_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  graphBackground2NB->Write(Form("bkg2_NB_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  graphBackground2NC->Write(Form("bkg2_NC_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  Purity_file->Close();
	
	
	
	
	
    //==================================================================================
    //MEAN #Kappa PLOT
    //==================================================================================
    
    
    TCanvas* c4 = new TCanvas("c4","",800,800);
    c4->SetLeftMargin(0.08);
    c4->SetTopMargin(0.05);
    c4->SetBottomMargin(0.08);
    c4->cd();
    TGraphErrors* gr20 = new TGraphErrors(18,centralbinvalues,mean00,binwidths,0);
    gr20->RemovePoint(0);
    gr20->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    gr20->GetXaxis()->SetTitleOffset(1.1);
    gr20->GetXaxis()->SetRangeUser(0,7);
    gr20->GetYaxis()->SetTitle("<#Kappa>");
    gr20->GetYaxis()->SetRangeUser(0,2);
    gr20->SetTitle("");
    gr20->SetMarkerColor(12);
    gr20->SetMarkerStyle(20);
    gr20->SetLineColor(12);
    gr20->Draw("AP");
    
    TGraphErrors* gr21 = new TGraphErrors(18,centralbinvalues,mean01,binwidths,0);
    gr21->RemovePoint(0);
    gr21->SetMarkerColor(kRed+1);
    gr21->SetMarkerStyle(24);
    gr21->SetLineColor(kRed+2);
    gr21->Draw("P");
    
    TGraphErrors* gr22 = new TGraphErrors(18,centralbinvalues,mean11,binwidths,0);
    gr22->RemovePoint(0);
    gr22->SetMarkerColor(kMagenta-2);
    gr22->SetMarkerStyle(24);
    gr22->SetLineColor(kMagenta-2);
    gr22->Draw("P");
    
    TGraphErrors* gr23 = new TGraphErrors(18,centralbinvalues,mean13,binwidths,0);
    gr23->RemovePoint(0);
    gr23->SetMarkerColor(kBlue-7);
    gr23->SetMarkerStyle(24);
    gr23->SetLineColor(kBlue-7);
    gr23->Draw("P");
    
    TLegend* leg3 = new TLegend(0.7,0.7,0.9,0.9);
    leg3->AddEntry(gr20,"real e^{+}e^{-}","lp");
    leg3->AddEntry(gr21,"remaining","lp");
    leg3->AddEntry(gr22,"#pi^{#pm} + #pi^{#mp}","lp");
    leg3->AddEntry(gr23,"#pi^{#pm} + e^{#mp}","lp");
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.04);
    leg3->Draw();
    
    TLatex T2;
    T2.SetTextSize(0.04);
    T2.SetTextAlign(12);
    T2.SetNDC();
    T2.DrawLatex(0.15, 0.85, "mean #Kappa");
    T2.SetTextSize(0.03);
    T2.DrawLatex(0.15, 0.97, InfoSystem.Data());
    if(kMC){
      T2.DrawLatex(0.65, 0.97, InfoMC.Data());
    }
    else{
      T2.DrawLatex(0.45, 0.97, InfoData.Data());
    }
    
    if(kMC){
      c4->SaveAs(Form("purity_studies_%s/MC_Kappa_4Templates_pt_Mean.eps",Cutnumber.Data()));
    }else{
      c4->SaveAs(Form("purity_studies_%s/Data_Kappa_4Templates_pt_Mean.eps",Cutnumber.Data()));
    }
	
	
	
    
    
    //==================================================================================
    //Kappa for all templates
    //==================================================================================
    
    if(kAllKappaPlots){
      
      TH2F* h2d00R = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp0_ee");
      TH2F* h2d01R = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp10_rem10");
      TH2F* h2d02R = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp8_had");
      TH2F* h2d11R = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp1_pipi");
      TH2F* h2d12R = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp4_pip");
      TH2F* h2d13R = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp2_pie");
      TH2F* h2d14R = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp7_KK");
      TH2F* h2d16R = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp6_ep");
      TH2F* h2d17R = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp5_eK");
      TH2F* h2d18R = (TH2F*)listMC_3->FindObject("hKappaTPC_Temp3_piK");
      
//       TH2F* h2dtotal = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPReduced");
//       TH2F* h2d10R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind10Reduced");
    
      TCanvas* c6 = new TCanvas("c6","",1500,600);
      c6->SetLeftMargin(0.08);
      c6->SetTopMargin(0.05);
      c6->SetBottomMargin(0.08);
      c6->Divide(5,2,0.001,0.001);
      gStyle->SetOptStat(0);
      
      TCanvas* c7 = new TCanvas("c7","",800,800);
      c7->SetLeftMargin(0.08);
      c7->SetTopMargin(0.05);
      c7->SetBottomMargin(0.08);
      
      TCanvas* c9 = new TCanvas("c9","",800,800);
      c9->SetLeftMargin(0.08);
      c9->SetTopMargin(0.05);
      c9->SetBottomMargin(0.08);
      
      Double_t SignalElectronComb[18];

      for(int i=0; i<18; i++){
    
        Double_t PtLow 	= newBinsComb[i];
        Double_t PtHigh 	= newBinsComb[i+1];
        
    // 	  Double_t nSigmaLow	= -20.0;
    // 	  Double_t nSigmaHigh	= 20.0;
      
        h2d00R->GetYaxis()->SetRangeUser(PtLow,PtHigh);
        h2d01R->GetYaxis()->SetRangeUser(PtLow,PtHigh);
        h2d02R->GetYaxis()->SetRangeUser(PtLow,PtHigh);
        h2d11R->GetYaxis()->SetRangeUser(PtLow,PtHigh);
        h2d12R->GetYaxis()->SetRangeUser(PtLow,PtHigh);
        h2d13R->GetYaxis()->SetRangeUser(PtLow,PtHigh);
        h2d14R->GetYaxis()->SetRangeUser(PtLow,PtHigh);
        h2d16R->GetYaxis()->SetRangeUser(PtLow,PtHigh);
        h2d17R->GetYaxis()->SetRangeUser(PtLow,PtHigh);
        h2d18R->GetYaxis()->SetRangeUser(PtLow,PtHigh);
//         h2dtotal->GetYaxis()->SetRangeUser(PtLow,PtHigh);
//         h2d10R->GetYaxis()->SetRangeUser(PtLow,PtHigh);

        TH1D* h1d00R = (TH1D*)h2d00R->ProjectionX();
        TH1D* h1d01R = (TH1D*)h2d01R->ProjectionX();
        TH1D* h1d02R = (TH1D*)h2d02R->ProjectionX();
        TH1D* h1d11R = (TH1D*)h2d11R->ProjectionX();
        TH1D* h1d12R = (TH1D*)h2d12R->ProjectionX();
        TH1D* h1d13R = (TH1D*)h2d13R->ProjectionX();
        TH1D* h1d14R = (TH1D*)h2d14R->ProjectionX();
        TH1D* h1d16R = (TH1D*)h2d16R->ProjectionX();
        TH1D* h1d17R = (TH1D*)h2d17R->ProjectionX();
        TH1D* h1d18R = (TH1D*)h2d18R->ProjectionX();
//         TH1D* h1dtotal = (TH1D*)h2dtotal->ProjectionX();
//         TH1D* h1d10R = (TH1D*)h2d10R->ProjectionX();
        
        h1d00R->SetTitle("");
        h1d00R->GetYaxis()->SetTitle("counts");
        h1d00R->GetXaxis()->SetTitle("#Kappa");
        h1d00R->GetXaxis()->SetTitleSize(0.05);
        h1d00R->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
        h1d00R->GetXaxis()->SetTitleOffset(0.8);
        h1d00R->SetFillColor(14);
        h1d00R->SetLineColor(14);
        h1d00R->SetFillStyle(3004);
        
        h1d01R->SetTitle("");
        h1d01R->GetYaxis()->SetTitle("counts");
        h1d01R->GetXaxis()->SetTitle("#Kappa");
        h1d01R->GetXaxis()->SetTitleSize(0.05);
        h1d01R->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
        h1d01R->GetXaxis()->SetTitleOffset(0.8);
        h1d01R->SetFillColor(kRed+1);
        h1d01R->SetLineColor(kRed+1);
        h1d01R->SetFillStyle(3005);
        
        h1d02R->SetTitle("");
        h1d02R->GetYaxis()->SetTitle("counts");
        h1d02R->GetXaxis()->SetTitle("#Kappa");
        h1d02R->GetXaxis()->SetTitleSize(0.05);
        h1d02R->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
        h1d02R->GetXaxis()->SetTitleOffset(0.8);
        h1d02R->SetFillColor(kGreen+2);
        h1d02R->SetLineColor(kGreen+2);
        h1d02R->SetFillStyle(3004);
        
        h1d11R->SetTitle("");
        h1d11R->GetYaxis()->SetTitle("counts");
        h1d11R->GetXaxis()->SetTitle("#Kappa");
        h1d11R->GetXaxis()->SetTitleSize(0.05);
        h1d11R->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
        h1d11R->GetXaxis()->SetTitleOffset(0.8);
        h1d11R->SetFillColor(kBlue-7);
        h1d11R->SetLineColor(kBlue-7);
        h1d11R->SetFillStyle(3005);
        
        h1d12R->SetTitle("");
        h1d12R->GetYaxis()->SetTitle("counts");
        h1d12R->GetXaxis()->SetTitle("#Kappa");
        h1d12R->GetXaxis()->SetTitleSize(0.05);
        h1d12R->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
        h1d12R->GetXaxis()->SetTitleOffset(0.8);
        h1d12R->SetFillColor(807);
        h1d12R->SetLineColor(807);
        h1d12R->SetFillStyle(3005);
        
        h1d13R->SetTitle("");
        h1d13R->GetYaxis()->SetTitle("counts");
        h1d13R->GetXaxis()->SetTitle("#Kappa");
        h1d13R->GetXaxis()->SetTitleSize(0.05);
        h1d13R->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
        h1d13R->GetXaxis()->SetTitleOffset(0.8);
        h1d13R->SetFillColor(kMagenta-2);
        h1d13R->SetLineColor(kMagenta-2);
        h1d13R->SetFillStyle(3004);
        
        h1d14R->SetTitle("");
        h1d14R->GetYaxis()->SetTitle("counts");
        h1d14R->GetXaxis()->SetTitle("#Kappa");
        h1d14R->GetXaxis()->SetTitleSize(0.05);
        h1d14R->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
        h1d14R->GetXaxis()->SetTitleOffset(0.8);
        h1d14R->SetFillColor(kCyan-2);
        h1d14R->SetLineColor(kCyan-2);
        h1d14R->SetFillStyle(3005);
        
        h1d16R->SetTitle("");
        h1d16R->GetYaxis()->SetTitle("counts");
        h1d16R->GetXaxis()->SetTitle("#Kappa");
        h1d16R->GetXaxis()->SetTitleSize(0.05);
        h1d16R->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
        h1d16R->GetXaxis()->SetTitleOffset(0.8);
        h1d16R->SetFillColor(kRed-4);
        h1d16R->SetLineColor(kRed-4);
        h1d16R->SetFillStyle(3005);
        
        h1d17R->SetTitle("");
        h1d17R->GetYaxis()->SetTitle("counts");
        h1d17R->GetXaxis()->SetTitle("#Kappa");
        h1d17R->GetXaxis()->SetTitleSize(0.05);
        h1d17R->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
        h1d17R->GetXaxis()->SetTitleOffset(0.8);
        h1d17R->SetFillColor(kGreen-5);
        h1d17R->SetLineColor(kGreen-5);
        h1d17R->SetFillStyle(3004);
        
        h1d18R->SetTitle("");
        h1d18R->GetYaxis()->SetTitle("counts");
        h1d18R->GetXaxis()->SetTitle("#Kappa");
        h1d18R->GetXaxis()->SetTitleSize(0.05);
        h1d18R->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
        h1d18R->GetXaxis()->SetTitleOffset(0.8);
        h1d18R->SetFillColor(kBlue+2);
        h1d18R->SetLineColor(kBlue+2);
        h1d18R->SetFillStyle(3004);
        
//         h1dtotal->SetTitle("");
//         h1dtotal->GetYaxis()->SetTitle("counts");
//         h1dtotal->GetXaxis()->SetTitle("#Kappa");
//         h1dtotal->GetXaxis()->SetTitleSize(0.05);
//         h1dtotal->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
//         h1dtotal->GetXaxis()->SetTitleOffset(0.8);
//         h1dtotal->SetFillColor(14);
//         h1dtotal->SetLineColor(14);
//         h1dtotal->SetFillStyle(3004);
        
        Double_t NTotal = h1d00R->GetEntries()+h1d01R->GetEntries()+h1d02R->GetEntries()+h1d11R->GetEntries()+h1d12R->GetEntries()+h1d13R->GetEntries()+h1d14R->GetEntries()+h1d16R->GetEntries()+h1d17R->GetEntries()+h1d18R->GetEntries();
        
        TLatex T4;
        T4.SetTextSize(0.06);
        T4.SetTextAlign(12);
        T4.SetNDC();

        c6->cd(1);
        h1d00R->Draw("HIST");
        T4.DrawLatex(0.65, 0.9, "signal");
        T4.DrawLatex(0.65, 0.85, Form("%f",h1d00R->GetEntries()/NTotal));
        T4.SetTextSize(0.03);
        T4.DrawLatex(0.18, 0.90, InfoSystem.Data());
        T4.DrawLatex(0.18, 0.85, InfoMC.Data());
        T4.DrawLatex(0.18, 0.75, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
        T4.DrawLatex(0.18, 0.65, SigmaStarForm.Data());
        
        c6->cd(10);
        h1d01R->Draw("HIST");
        T4.SetTextSize(0.06);
        T4.DrawLatex(0.65, 0.9, "remaining");
        T4.DrawLatex(0.65, 0.85, Form("%f",h1d01R->GetEntries()/NTotal));
        
        c6->cd(9);
        h1d02R->Draw("HIST");
        T4.DrawLatex(0.65, 0.9, "hadronic");
        T4.DrawLatex(0.65, 0.85, Form("%f",h1d02R->GetEntries()/NTotal));
        
        c6->cd(3);
        h1d11R->Draw("HIST");
        T4.DrawLatex(0.65, 0.9, "#pi^{#pm} + #pi^{#mp}");
        T4.DrawLatex(0.65, 0.85, Form("%f",h1d11R->GetEntries()/NTotal));
        
        c6->cd(5);
        h1d12R->Draw("HIST");
        T4.DrawLatex(0.65, 0.9, "#pi^{#pm} + p(#bar{p})");
        T4.DrawLatex(0.65, 0.85, Form("%f",h1d12R->GetEntries()/NTotal));
        
        c6->cd(2);
        h1d13R->Draw("HIST");
        T4.DrawLatex(0.65, 0.9, "#pi^{#pm} + e^{#mp}");
        T4.DrawLatex(0.65, 0.85, Form("%f",h1d13R->GetEntries()/NTotal));
        
        c6->cd(8);
        h1d14R->Draw("HIST");
        T4.DrawLatex(0.65, 0.9, "K^{#pm} + K^{#mp}");
        T4.DrawLatex(0.65, 0.85, Form("%f",h1d14R->GetEntries()/NTotal));
        
        c6->cd(7);
        h1d16R->Draw("HIST");
        T4.DrawLatex(0.65, 0.9,  "e^{#pm} + p(#bar{p})");
        T4.DrawLatex(0.65, 0.85, Form("%f",h1d16R->GetEntries()/NTotal));
        
        c6->cd(6);
        h1d17R->Draw("HIST");
        T4.DrawLatex(0.65, 0.9,  "e^{#pm} + K^{#mp}");
        T4.DrawLatex(0.65, 0.85, Form("%f",h1d17R->GetEntries()/NTotal));
        
        c6->cd(4);
        h1d18R->Draw("HIST");
        T4.DrawLatex(0.65, 0.9,  "#pi^{#pm} + K^{#mp}");
        T4.DrawLatex(0.65, 0.85, Form("%f",h1d18R->GetEntries()/NTotal));
        
        TH1D *h1sum=(TH1D*)h1d00R->Clone();
        h1sum->GetXaxis()->SetTitleOffset(0.8);
        h1sum->SetFillColor(0);
        h1sum->SetLineColor(kBlack);
        h1sum->SetLineWidth(2.0);
        h1sum->Add(h1d01R);
        h1sum->Add(h1d02R);
        h1sum->Add(h1d11R);
        h1sum->Add(h1d12R);
        h1sum->Add(h1d13R);
        h1sum->Add(h1d14R);
        h1sum->Add(h1d16R);
        h1sum->Add(h1d17R);
        h1sum->Add(h1d18R);
        
        c7->cd();
        h1sum->Draw("HIST");
        h1d00R->Draw("HIST same");
        h1d01R->Draw("HIST same");
        h1d02R->Draw("HIST same");
        h1d11R->Draw("HIST same");
        h1d12R->Draw("HIST same");
        h1d13R->Draw("HIST same");
        h1d14R->Draw("HIST same");
        h1d16R->Draw("HIST same");
        h1d17R->Draw("HIST same");
        h1d18R->Draw("HIST same");
        
        TLegend* leg = new TLegend(0.7,0.5,0.9,0.9);
        leg->SetTextSize(0.04);
        leg->AddEntry(h1sum,"total","lp");
        leg->AddEntry(h1d00R,"real e^{+}e^{-}","f");
        leg->AddEntry(h1d01R,"remaining","f");
        leg->AddEntry(h1d02R,"hadronic","f");
        leg->AddEntry(h1d11R,"#pi^{#pm} + #pi^{#mp}","f");
        leg->AddEntry(h1d12R,"#pi^{#pm} + p(#bar{p})","f");
        leg->AddEntry(h1d13R,"#pi^{#pm} + e^{#mp}","f");
        leg->AddEntry(h1d14R,"K^{#pm} + K^{#mp}","f");
        leg->AddEntry(h1d16R,"e^{#pm} + p(#bar{p})","f");
        leg->AddEntry(h1d17R,"e^{#pm} + K^{#mp}","f");
        leg->AddEntry(h1d18R,"#pi^{#pm} + K^{#mp}","f");
        leg->SetBorderSize(0);
        leg->Draw();
        
        T4.SetTextSize(0.03);
        T4.DrawLatex(0.15, 0.97, InfoSystem.Data());
        T4.DrawLatex(0.65, 0.97, InfoMC.Data());
        T1.DrawLatex(0.12, 0.88, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
    // 	  T4.DrawLatex(0.65, 0.3,  SigmaStarForm.Data());
        
//         c9->cd();
//         gPad->SetLogy();
//         Double_t h1dtotalInt, h1d10RInt;
//         h1dtotalInt = h1dtotal->Integral(h1dtotal->GetXaxis()->FindBin(-5.0),h1dtotal->GetXaxis()->FindBin(10.0));
//         h1d10RInt   = h1d10R->Integral(h1d10R->GetXaxis()->FindBin(-5.0),h1d10R->GetXaxis()->FindBin(10.0));
//         SignalElectronComb[i] = h1d10RInt/h1dtotalInt;
//         h1dtotal->Draw("HIST");
//         h1d10R->Draw("HIST same");
//         T4.DrawLatex(0.11, 0.87, Form("#frac{N_{e^{+}e^{-}}}{N_{tot}}(-5<K<10) = %2.3f %%",100*h1d10RInt/h1dtotalInt));
//         T4.DrawLatex(0.15, 0.97, InfoSystem.Data());
//         T4.DrawLatex(0.65, 0.97, InfoMC.Data());
//         T4.DrawLatex(0.11, 0.15, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
        //T4.DrawLatex(0.65, 0.5,  SigmaStarForm.Data());
        
//         TLegend* leg6 = new TLegend(0.7,0.8,0.9,0.9);
//         leg6->SetTextSize(0.04);
//         leg6->AddEntry(h1dtotal,"total","f");
//         leg6->AddEntry(h1d10R,"e^{+}e^{-} bkg","f");
//         leg6->SetBorderSize(0);
//         leg6->Draw();
        
        c6->SaveAs(Form("purity_studies_%s/MC_Kappa_AllTemplates_pt_%s-%s.eps",Cutnumber.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
        c7->SaveAs(Form("purity_studies_%s/MC_Kappa_AllTemplates_same_pt_%s-%s.eps",Cutnumber.Data(),newBinsName[i].Data(),newBinsName[i+1].Data()));
//         c9->SaveAs(Form("purity_studies/MC_Kappa_AllandElectronComb_pt_%s-%s.eps",newBinsName[i].Data(),newBinsName[i+1].Data()));

      }
      
//       TCanvas* c12 = new TCanvas("c12","",800,800);
//       c12->SetLeftMargin(0.08);
//       c12->SetTopMargin(0.05);
//       c12->SetBottomMargin(0.08);
//       
//       TGraphErrors* graphElectronComb = new TGraphErrors(18,centralbinvalues,SignalElectronComb,0,0);
//       graphElectronComb->RemovePoint(0);
//       graphElectronComb->GetXaxis()->SetTitle("p_{T}(GeV/c)");
//       graphElectronComb->GetXaxis()->SetRangeUser(0.0,7);
//       graphElectronComb->GetXaxis()->SetTitleOffset(1.1);
//       graphElectronComb->GetYaxis()->SetTitle("c");
//       graphElectronComb->GetYaxis()->SetRangeUser(0.0,0.025);
//       graphElectronComb->SetTitle("");
//       graphElectronComb->SetMarkerColor(12);
//       graphElectronComb->SetMarkerStyle(20);
//       graphElectronComb->SetLineColor(12);
//       graphElectronComb->Draw("AP");
//       
//       T1.DrawLatex(0.15, 0.97, InfoSystem.Data());
//       T1.DrawLatex(0.65, 0.97, InfoMC.Data());
//       T1.DrawLatex(0.15, 0.9, "e^{+}e^{-} combinatorics (-5<K<10)");
//       
//       if(kMC){
//       c12->SaveAs("purity_studies/MC_ElectronComb.eps");
//       
//       TFile *ElectronComb_file = new TFile(Form("ElectronComb_%s%s.root",CentralityLow.Data(),CentralityHigh.Data()),"RECREATE");
//       graphElectronComb->Write(Form("graphElectronComb_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
//       ElectronComb_file->Close();
//     }
    
    }
    
}
