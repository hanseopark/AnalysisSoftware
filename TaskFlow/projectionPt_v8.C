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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"

void  projectionPt_v8(){
  
  
    //==================================================================================
    //DEFINING SOME STUFF
    //==================================================================================
  
    Bool_t kMC = 1;
    Bool_t kKappavsKappaPlots = 1;
    Bool_t kAllKappaPlots = 1;
    Bool_t kUseFitConstraints = 1;
    
    Double_t nSigmaLow  = 0.0;
    Double_t nSigmaHigh = 8.0;
    
    StyleSettingsThesis();  
    SetPlotStyle();
    
    TString InfoSystem = "20-40% PbPb, #sqrt{s}=2.76TeV";
    TString InfoData = "Data LHC10h";
    TString InfoMC = "MC LHC13d2";
    TString SigmaStarForm = "#Kappa = #frac{|#kappa^{+}|+|#kappa^{-}|}{2}";
//     TString SigmaStarForm = "#Kappa = #frac{|#kappa^{+}|+|#kappa^{-}|}{2}+2(#kappa^{+}+#kappa^{-})";
    //TString SigmaStarForm = "#Kappa = #sqrt{#frac{#kappa^{+}^{2}+#kappa^{+}^{2}}{2}}";
    
//     TFile* fileData = new TFile("PhotonQA_0705415160_redefinedkappalarge_MC.root");
    TFile* fileData = new TFile("PhotonQA_0705415160_MC.root");
    TDirectory* dir 	= new TDirectory();
    dir 		= (TDirectory*)fileData->Get("GammaConv_0705415160");

    TFile* fileData2  = new TFile("PhotonQA_1_Data.root");
    TDirectory* dir2 	= new TDirectory();
    dir2 		= (TDirectory*)fileData2->Get("GammaConv_0705314140");    

    TH2F* h2d00 = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind00Reduced");
    TH2F* h2d01 = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKindBkg2Reduced");
    TH2F* h2d11 = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind11Reduced");
    TH2F* h2d13 = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind13Reduced");
    
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
    
    Double_t newBinsComb[19] = {0.0, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.3, 3.7, 4.1, 4.6, 5.4, 6.0, 7.0};
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
	TH1D* h1d01 = (TH1D*)h2d01->ProjectionX();
	TH1D* h1d11 = (TH1D*)h2d11->ProjectionX();
	TH1D* h1d13 = (TH1D*)h2d13->ProjectionX();
	
	h1d00->SetTitle("");
	h1d00->GetYaxis()->SetTitle("counts");
	h1d00->GetXaxis()->SetTitle("#Kappa");
	h1d00->GetXaxis()->SetTitleSize(0.05);
	h1d00->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	h1d00->GetXaxis()->SetTitleOffset(0.8);
	h1d00->SetFillColor(12);
	h1d00->SetLineColor(12);
	h1d00->SetFillStyle(3004);
	
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
	  h2d_data = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPReduced");
	}else{
	  h2d_data = (TH2F*)dir2->Get("histonSigdEdxElnSigdEdxPosGammaPReduced");
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
	T1.SetTextSize(0.04);
	T1.DrawLatex(0.5, 0.45, InfoSystem.Data());
	T1.DrawLatex(0.5, 0.4, InfoMC.Data());
	T1.DrawLatex(0.5, 0.2, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
	T1.DrawLatex(0.5, 0.3,  SigmaStarForm.Data());
	
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
  
  TLegend* leg5 = new TLegend(0.7,0.6,0.9,0.9);
  leg5->SetTextSize(0.04);
//   leg5->SetHeader("Ratio to sum");
  leg5->AddEntry(h1ratiobkg1,"remaining","l");
  leg5->AddEntry(h1ratiobkg2,"#pi^{#pm} + #pi^{#mp}","l");
  leg5->AddEntry(h1ratiobkg3,"#pi^{#pm} + e^{#mp}","l");
  leg5->SetBorderSize(0);
  leg5->Draw();
	
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
	T1.DrawLatex(0.65, 0.2, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
	T1.DrawLatex(0.65, 0.3,  SigmaStarForm.Data());
  
  c8->cd();
  h1ratiobkg1->Draw("HIST");
  h1ratiobkg2->Draw("HIST same");
  h1ratiobkg3->Draw("HIST same");
  
  leg5->Draw();
  
  T1.SetTextSize(0.03);
  T1.DrawLatex(0.15, 0.97, InfoSystem.Data());
  T1.DrawLatex(0.65, 0.97, InfoMC.Data());
  T1.DrawLatex(0.15, 0.7, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
  T1.DrawLatex(0.15, 0.8,  SigmaStarForm.Data());
  T1.DrawLatex(0.15, 0.9, Form("background/totalbackground"));
	
	
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
	  fracwidth = 0.9;
	  constraint_low  = 1-fracwidth;
	  constraint_high = 1+fracwidth;
	}else{
	  NEntriesData = data->GetEntries();
	  frac00[i] = h1d00->GetEntries()/NTotalMC;
	  frac01[i] = h1d01->GetEntries()/NTotalMC;
	  frac11[i] = h1d11->GetEntries()/NTotalMC;
	  frac13[i] = h1d13->GetEntries()/NTotalMC;
	  fracwidth = 0.9;
	  constraint_low  = 1-fracwidth;
	  constraint_high = 1+fracwidth;
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
	  
	  values01[i]=value[1];
	  errors01[i]=error[1];
	  
	  values11[i]=value[2];
	  errors11[i]=error[2];
	  
	  values13[i]=value[3];
	  errors13[i]=error[3];
	  
	  
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
	}
	
	Int_t IntLowerBound = h1d00->GetXaxis()->FindBin(0.0);
	Int_t IntUpperBound1 = h1d00->GetXaxis()->FindBin(1.0);
	Int_t IntUpperBound2 = h1d00->GetXaxis()->FindBin(2.0);
	Int_t IntUpperBound3 = h1d00->GetXaxis()->FindBin(3.0);
	Int_t IntUpperBound4 = h1d00->GetXaxis()->FindBin(4.0);
	h1d00Int = h1d00->Integral(IntLowerBound,IntUpperBound1)*values00[i]/frac00[i];
	h1d01Int = h1d01->Integral(IntLowerBound,IntUpperBound1)*values01[i]/frac01[i];
	h1d11Int = h1d11->Integral(IntLowerBound,IntUpperBound1)*values11[i]/frac11[i];
	h1d13Int = h1d13->Integral(IntLowerBound,IntUpperBound1)*values13[i]/frac13[i];
	PuritySignal1[i] = h1d00Int/(h1d00Int+h1d01Int+h1d11Int+h1d13Int);
	h1d00Int = h1d00->Integral(IntLowerBound,IntUpperBound2)*values00[i]/frac00[i];
	h1d01Int = h1d01->Integral(IntLowerBound,IntUpperBound2)*values01[i]/frac01[i];
	h1d11Int = h1d11->Integral(IntLowerBound,IntUpperBound2)*values11[i]/frac11[i];
	h1d13Int = h1d13->Integral(IntLowerBound,IntUpperBound2)*values13[i]/frac13[i];
	PuritySignal2[i] = h1d00Int/(h1d00Int+h1d01Int+h1d11Int+h1d13Int);
	h1d00Int = h1d00->Integral(IntLowerBound,IntUpperBound3)*values00[i]/frac00[i];
	h1d01Int = h1d01->Integral(IntLowerBound,IntUpperBound3)*values01[i]/frac01[i];
	h1d11Int = h1d11->Integral(IntLowerBound,IntUpperBound3)*values11[i]/frac11[i];
	h1d13Int = h1d13->Integral(IntLowerBound,IntUpperBound3)*values13[i]/frac13[i];
	PuritySignal3[i] = h1d00Int/(h1d00Int+h1d01Int+h1d11Int+h1d13Int);
	h1d00Int = h1d00->Integral(IntLowerBound,IntUpperBound4)*values00[i]/frac00[i];
	h1d01Int = h1d01->Integral(IntLowerBound,IntUpperBound4)*values01[i]/frac01[i];
	h1d11Int = h1d11->Integral(IntLowerBound,IntUpperBound4)*values11[i]/frac11[i];
	h1d13Int = h1d13->Integral(IntLowerBound,IntUpperBound4)*values13[i]/frac13[i];
	PuritySignal4[i] = h1d00Int/(h1d00Int+h1d01Int+h1d11Int+h1d13Int);

	gSystem->mkdir("purity_studies");
	if(kMC){
	  c1->SaveAs(Form("purity_studies/MC_Kappa_4Templates_pt_%s-%s.eps",newBinsName[i].Data(),newBinsName[i+1].Data()));
 	  c2->SaveAs(Form("purity_studies/MC_Kappa_4Templates_same_pt_%s-%s.eps",newBinsName[i].Data(),newBinsName[i+1].Data()));
    c8->SaveAs(Form("purity_studies/MC_Kappa_4Templates_ratiobkg_pt_%s-%s.eps",newBinsName[i].Data(),newBinsName[i+1].Data()));
	}else{
	  c1->SaveAs(Form("purity_studies/Data_Kappa_4Templates_pt_%s-%s.eps",newBinsName[i].Data(),newBinsName[i+1].Data()));
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
    gr0->GetYaxis()->SetRangeUser(0,1);
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
    TGraphErrors* gr12 = new TGraphErrors(18,centralbinvalues,PuritySignal1,binwidths,0);
    gr12->RemovePoint(0);
    gr12->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    gr12->GetXaxis()->SetRangeUser(0.0,7);
    gr12->GetXaxis()->SetTitleOffset(1.1);
    gr12->GetYaxis()->SetTitle("Purity");
    gr12->GetYaxis()->SetRangeUser(0.0,1.0);
    gr12->SetTitle("");
    gr12->SetMarkerColor(12);
    gr12->SetMarkerStyle(20);
    gr12->SetLineColor(12);
    gr12->Draw("AP");
    
    TGraphErrors* gr13 = new TGraphErrors(18,centralbinvalues,PuritySignal2,binwidths,0);
    gr13->RemovePoint(0);
    gr13->SetMarkerColor(22);
    gr13->SetMarkerStyle(20);
    gr13->SetLineColor(22);
    gr13->Draw("P");
    
    TGraphErrors* gr14 = new TGraphErrors(18,centralbinvalues,PuritySignal3,binwidths,0);
    gr14->RemovePoint(0);
    gr14->SetMarkerColor(32);
    gr14->SetMarkerStyle(20);
    gr14->SetLineColor(32);
    gr14->Draw("P");
    
    TGraphErrors* gr15 = new TGraphErrors(18,centralbinvalues,PuritySignal4,binwidths,0);
    gr15->RemovePoint(0);
    gr15->SetMarkerColor(42);
    gr15->SetMarkerStyle(20);
    gr15->SetLineColor(42);
    gr15->Draw("P");
    
    TLegend* leg4 = new TLegend(0.8,0.1,0.9,0.3);
    leg4->AddEntry(gr12,"#Kappa<1","lp");
    leg4->AddEntry(gr13,"#Kappa<2","lp");
    leg4->AddEntry(gr14,"#Kappa<3","lp");
    leg4->AddEntry(gr15,"#Kappa<4","lp");
    leg4->SetBorderSize(0);
    leg4->SetTextSize(0.04);
    leg4->Draw();
    
    T1.DrawLatex(0.2, 0.15, "Purity #gamma sample");
    
    if(kMC){
      c3->SaveAs("purity_studies/MC_fractionfitter_4Templates.eps");
    }else{
      c3->SaveAs("purity_studies/Data_fractionfitter_4Templates.eps");
    }
	
	
	
	
	
	
	
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
    gr20->GetYaxis()->SetRangeUser(0,10);
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
      c4->SaveAs("purity_studies/MC_Kappa_4Templates_pt_Mean.eps");
    }else{
      c4->SaveAs("purity_studies/Data_Kappa_4Templates_pt_Mean.eps");
    }
	
	
	
    //==================================================================================
    //nSigma vs nSigma plots
    //==================================================================================
    
    if(kKappavsKappaPlots){
    
      TH3F* h3d00 = (TH3F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind00");
      TH3F* h3d01 = (TH3F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKindBkg");
      TH3F* h3d02 = (TH3F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind02");
      TH3F* h3d11 = (TH3F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind11");
      TH3F* h3d12 = (TH3F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind12");
      TH3F* h3d13 = (TH3F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind13");
      TH3F* h3d14 = (TH3F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind14");
      TH3F* h3d16 = (TH3F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind16");
      TH3F* h3d17 = (TH3F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind17");
      TH3F* h3d18 = (TH3F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind18");
    
      TCanvas* c5 = new TCanvas("c5","",1500,600);
      gStyle->SetOptStat(0);
      c5->Divide(5,2,0.001,0.001);
      

      for(int i=0; i<18; i++){
    
	  Double_t PtLow 	= newBinsComb[i];
	  Double_t PtHigh 	= newBinsComb[i+1];
	  
// 	  Double_t nSigmaLow	= -20.0;
// 	  Double_t nSigmaHigh	= 20.0;
	
	  h3d00->GetZaxis()->SetRangeUser(PtLow,PtHigh);
	  h3d01->GetZaxis()->SetRangeUser(PtLow,PtHigh);
	  h3d02->GetZaxis()->SetRangeUser(PtLow,PtHigh);
	  h3d11->GetZaxis()->SetRangeUser(PtLow,PtHigh);
	  h3d12->GetZaxis()->SetRangeUser(PtLow,PtHigh);
	  h3d13->GetZaxis()->SetRangeUser(PtLow,PtHigh);
	  h3d14->GetZaxis()->SetRangeUser(PtLow,PtHigh);
	  h3d16->GetZaxis()->SetRangeUser(PtLow,PtHigh);
	  h3d17->GetZaxis()->SetRangeUser(PtLow,PtHigh);
	  h3d18->GetZaxis()->SetRangeUser(PtLow,PtHigh);
	  
	  TLatex T3;
	  T3.SetTextSize(0.08);
	  T3.SetTextAlign(12);
	  
	  TBox *box1 = new TBox(0,(nSigmaHigh-3),9.5,(nSigmaHigh-1));
	  box1->SetFillColor(19);

	  c5->cd(1);
	  TH2D* h2d00p = (TH2D*)h3d00->Project3D("xyo");
	  h2d00p->SetTitle("");
	  h2d00p->GetXaxis()->SetTitle("#kappa^{-}");
	  h2d00p->GetYaxis()->SetTitle("#kappa^{+}");
	  gPad->SetLogz();
	  h2d00p->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d00p->GetYaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d00p->Draw("colz");
	  box1->Draw();
	  T3.DrawLatex(1, (nSigmaHigh-2), "signal");
	  T3.SetTextSize(0.06);
	  T3.DrawLatex(-8, (nSigmaLow+4), InfoSystem.Data());
	  T3.DrawLatex(-8, (nSigmaLow+2.5), InfoMC.Data());
	  T3.DrawLatex(-8, (nSigmaLow+1), Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
	  
	  Double_t ScaleMin = h2d00p->GetMinimum();
	  Double_t ScaleMax = h2d00p->GetMaximum();
	  
	  c5->cd(2);
	  TH2D* h2d13p = (TH2D*)h3d13->Project3D("xyo");
	  h2d13p->SetTitle("");
	  h2d13p->GetXaxis()->SetTitle("#kappa^{-}");
	  h2d13p->GetYaxis()->SetTitle("#kappa^{+}");
	  gPad->SetLogz();
	  h2d13p->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d13p->GetYaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d13p->SetMinimum(ScaleMin);
	  h2d13p->SetMaximum(ScaleMax);
	  h2d13p->Draw("colz");
	  T3.SetTextSize(0.08);
	  box1->Draw();
	  T3.DrawLatex(1, (nSigmaHigh-2), "#pi^{#pm} + e^{#mp}");
	  
	  c5->cd(3);
	  TH2D* h2d11p = (TH2D*)h3d11->Project3D("xyo");
	  h2d11p->SetTitle("");
	  h2d11p->GetXaxis()->SetTitle("#kappa^{-}");
	  h2d11p->GetYaxis()->SetTitle("#kappa^{+}");
	  gPad->SetLogz();
	  h2d11p->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d11p->GetYaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d11p->SetMinimum(ScaleMin);
	  h2d11p->SetMaximum(ScaleMax);
	  h2d11p->Draw("colz");
	  box1->Draw();
	  T3.DrawLatex(1, (nSigmaHigh-2), "#pi^{#pm} + #pi^{#mp}");
	  
	  c5->cd(4);
	  TH2D* h2d18p = (TH2D*)h3d18->Project3D("xyo");
	  h2d18p->SetTitle("");
	  h2d18p->GetXaxis()->SetTitle("#kappa^{-}");
	  h2d18p->GetYaxis()->SetTitle("#kappa^{+}");
	  gPad->SetLogz();
	  h2d18p->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d18p->GetYaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d18p->SetMinimum(ScaleMin);
	  h2d18p->SetMaximum(ScaleMax);
	  h2d18p->Draw("colz");
	  box1->Draw();
	  T3.DrawLatex(1, (nSigmaHigh-2),  "#pi^{#pm} + K^{#mp}");
	  
	  c5->cd(5);
	  TH2D* h2d12p = (TH2D*)h3d12->Project3D("xyo");
	  h2d12p->SetTitle("");
	  h2d12p->GetXaxis()->SetTitle("#kappa^{-}");
	  h2d12p->GetYaxis()->SetTitle("#kappa^{+}");
	  gPad->SetLogz();
	  h2d12p->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d12p->GetYaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d12p->SetMinimum(ScaleMin);
	  h2d12p->SetMaximum(ScaleMax);
	  h2d12p->Draw("colz");
	  box1->Draw();
	  T3.DrawLatex(1, (nSigmaHigh-2), "#pi^{#pm} + p(#bar{p})");
	  
	  c5->cd(6);
	  TH2D* h2d17p = (TH2D*)h3d17->Project3D("xyo");
	  h2d17p->SetTitle("");
	  h2d17p->GetXaxis()->SetTitle("#kappa^{-}");
	  h2d17p->GetYaxis()->SetTitle("#kappa^{+}");
	  gPad->SetLogz();
	  h2d17p->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d17p->GetYaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d17p->SetMinimum(ScaleMin);
	  h2d17p->SetMaximum(ScaleMax);
	  h2d17p->Draw("colz");
	  box1->Draw();
	  T3.DrawLatex(1, (nSigmaHigh-2),  "e^{#pm} + K^{#mp}");
	  
	  c5->cd(7);
	  TH2D* h2d16p = (TH2D*)h3d16->Project3D("xyo");
	  h2d16p->SetTitle("");
	  h2d16p->GetXaxis()->SetTitle("#kappa^{-}");
	  h2d16p->GetYaxis()->SetTitle("#kappa^{+}");
	  gPad->SetLogz();
	  h2d16p->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d16p->GetYaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d16p->SetMinimum(ScaleMin);
	  h2d16p->SetMaximum(ScaleMax);
	  h2d16p->Draw("colz");
	  box1->Draw();
	  T3.DrawLatex(1, (nSigmaHigh-2),  "e^{#pm} + p(#bar{p})");
	  
	  c5->cd(8);
	  TH2D* h2d14p = (TH2D*)h3d14->Project3D("xyo");
	  h2d14p->SetTitle("");
	  h2d14p->GetXaxis()->SetTitle("#kappa^{-}");
	  h2d14p->GetYaxis()->SetTitle("#kappa^{+}");
	  gPad->SetLogz();
	  h2d14p->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d14p->GetYaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d14p->SetMinimum(ScaleMin);
	  h2d14p->SetMaximum(ScaleMax);
	  h2d14p->Draw("colz");
	  box1->Draw();
	  T3.DrawLatex(1, (nSigmaHigh-2), "K^{#pm} + K^{#mp}");
	  
	  c5->cd(9);
	  TH2D* h2d02p = (TH2D*)h3d02->Project3D("xyo");
	  h2d02p->SetTitle("");
	  h2d02p->GetXaxis()->SetTitle("#kappa^{-}");
	  h2d02p->GetYaxis()->SetTitle("#kappa^{+}");
	  gPad->SetLogz();
	  h2d02p->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d02p->GetYaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d02p->SetMinimum(ScaleMin);
	  h2d02p->SetMaximum(ScaleMax);
	  h2d02p->Draw("colz");
	  box1->Draw();
	  T3.DrawLatex(1, (nSigmaHigh-2), "hadronic");
	  
	  c5->cd(10);
	  TH2D* h2d01p = (TH2D*)h3d01->Project3D("xyo");
	  h2d01p->SetTitle("");
	  h2d01p->GetXaxis()->SetTitle("#kappa^{-}");
	  h2d01p->GetYaxis()->SetTitle("#kappa^{+}");
	  gPad->SetLogz();
	  h2d01p->GetXaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d01p->GetYaxis()->SetRangeUser(nSigmaLow,nSigmaHigh);
	  h2d01p->SetMinimum(ScaleMin);
	  h2d01p->SetMaximum(ScaleMax);
	  h2d01p->Draw("colz");
	  box1->Draw();
	  T3.DrawLatex(1, (nSigmaHigh-2), "remaining");
	  
	  c5->SaveAs(Form("purity_studies/MC_kappa_vs_kappa_AllTemplates_pt_%s-%s.eps",newBinsName[i].Data(),newBinsName[i+1].Data()));
      }
    
    }
    
    //==================================================================================
    //Kappa for all templates
    //==================================================================================
    
    if(kAllKappaPlots){
    
      TH2F* h2d00R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind00Reduced");
      TH2F* h2d01R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKindBkgReduced");
      TH2F* h2d02R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind02Reduced");
      TH2F* h2d11R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind11Reduced");
      TH2F* h2d12R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind12Reduced");
      TH2F* h2d13R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind13Reduced");
      TH2F* h2d14R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind14Reduced");
      TH2F* h2d16R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind16Reduced");
      TH2F* h2d17R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind17Reduced");
      TH2F* h2d18R = (TH2F*)dir->Get("histonSigdEdxElnSigdEdxPosGammaPKind18Reduced");
    
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
	  
	  Double_t NTotal = h1d00R->GetEntries()+h1d01R->GetEntries()+h1d02R->GetEntries()+h1d11R->GetEntries()+h1d12R->GetEntries()+h1d13R->GetEntries()+h1d14R->GetEntries()+h1d16R->GetEntries()+h1d17R->GetEntries()+h1d18R->GetEntries();
	  
	  TLatex T4;
	  T4.SetTextSize(0.06);
	  T4.SetTextAlign(12);
	  T4.SetNDC();

	  c6->cd(1);
	  h1d00R->Draw("HIST");
	  T4.DrawLatex(0.65, 0.9, "signal");
	  T4.DrawLatex(0.65, 0.85, Form("%f",h1d00R->GetEntries()/NTotal));
	  T4.SetTextSize(0.04);
	  T4.DrawLatex(0.5, 0.45, InfoSystem.Data());
	  T4.DrawLatex(0.5, 0.4, InfoMC.Data());
	  T4.DrawLatex(0.5, 0.2, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
	  T4.DrawLatex(0.5, 0.3,  SigmaStarForm.Data());
	  
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
	  T4.DrawLatex(0.65, 0.2, Form("%2.1f < P_{t}(GeV/c) < %2.1f",newBinsComb[i],newBinsComb[i+1]));
	  T4.DrawLatex(0.65, 0.3,  SigmaStarForm.Data());
	  
	  c6->SaveAs(Form("purity_studies/MC_Kappa_AllTemplates_pt_%s-%s.eps",newBinsName[i].Data(),newBinsName[i+1].Data()));
	  c7->SaveAs(Form("purity_studies/MC_Kappa_AllTemplates_same_pt_%s-%s.eps",newBinsName[i].Data(),newBinsName[i+1].Data()));

      }
    }    
}
