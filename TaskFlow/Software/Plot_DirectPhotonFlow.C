#include "DirectPhotonFlowFunctions.h"

void  Plot_DirectPhotonFlow(
                                      TString CentralityLow = "40",
                                      TString CentralityHigh = "80",
                                      TString InputFileName1  = "/home/mike/0_directphoton/13_DirectPhotonAnalysisCode/PCM_DirectPhotonFlow_nocorrection_4080.root",
                                      TString InputFileName2  = "/home/mike/0_directphoton/13_DirectPhotonAnalysisCode/PCM_DirectPhotonFlow_4080.root"
){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  //open datafile with v2 gamma inclusive uncorrected
  TFile* fileInclusive1  = new TFile(InputFileName1.Data());
  TH1F*  histoDirect1 = (TH1F*)fileInclusive1->Get("histoDirect");
  if(!histoDirect1) cout << "histoDirect1 not found in fileInclusive1!!" << endl;
  
  histoDirect1->SetMarkerStyle(21);
  histoDirect1->SetMarkerSize(1.2);
  histoDirect1->SetMarkerColor(kRed+2);
  histoDirect1->SetLineColor(kRed+2);
  
  TFile* fileInclusive2  = new TFile(InputFileName2.Data());
  TH1F*  histoDirect2 = (TH1F*)fileInclusive2->Get("histoDirect");
  if(!histoDirect2) cout << "histoDirect2 not found in fileInclusive2!!" << endl;
  
  histoDirect2->SetMarkerStyle(20);
  histoDirect2->SetMarkerSize(1.2);
  histoDirect2->SetMarkerColor(kGreen+2);
  histoDirect2->SetLineColor(kGreen+2);

  //========================================================
  //histogram cosmetics
  //========================================================
  
  TH1F*  histoDirectEmpty = (TH1F*)Getv2HistStyle();
  histoDirectEmpty->GetYaxis()->SetRangeUser(-0.29,0.49);
  
  //========================================================
  //Plotting and Saving
  //========================================================
  TLatex T1;
  T1.SetTextSize(0.04);
  T1.SetTextAlign(12);
  T1.SetNDC();
  
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  
  TLegend* leg = new TLegend(0.16,0.15,0.6,0.26);
  leg->SetTextSize(0.04);
  leg->AddEntry(histoDirect1,"uncorrected","lp");
  leg->AddEntry(histoDirect2,"background corrected","lp");
  leg->SetBorderSize(0);

  SetProperMargins();
  histoDirectEmpty->Draw();
  histoDirect1->Draw("SAME");
  histoDirect2->Draw("SAME");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  leg->Draw();
  
  //Masking lower pt points
  MaskPoints(histoDirect1,0,10);
  MaskPoints(histoDirect2,0,10);
  
  
  //Masking higher pt points
  MaskPoints(histoDirect1,26,29);
  MaskPoints(histoDirect2,26,29);

  gSystem->mkdir("Results");
  gSystem->mkdir("Results/v2GammaDirectComparison");
  c1->SaveAs(Form("Results/v2GammaDirectComparison/v2_Gamma_Direct_comparison_%s%s.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
}