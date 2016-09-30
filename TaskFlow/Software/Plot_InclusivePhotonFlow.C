#include "DirectPhotonFlowFunctions.h"

void  Plot_InclusivePhotonFlow(       
                                      TString CentralityLow = "20",
                                      TString CentralityHigh = "40",
                                      TString InputFileName1  = "/home/mike/0_directphoton/0_analysis/160303_PbPb_v2_widecents/InclusivePhotonv2_PCM_wideCC.root",
                                      TString InputFileName2  = "/home/mike/0_directphoton/13_DirectPhotonAnalysisCode/InclusivePhotonv2_uncorrected_52.root",
                                      TString InputFileName3  = "/home/mike/0_directphoton/13_DirectPhotonAnalysisCode/InclusivePhotonv2_uncorrected_50.root"
){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  //open datafile with v2 gamma inclusive uncorrected
  TFile* fileInclusive1  = new TFile(InputFileName1.Data());
  TH1F*  histoInclusive1 = (TH1F*)fileInclusive1->Get("2040_tC21");
  if(!histoInclusive1) cout << "histoInclusive1 not found in fileInclusive1!!" << endl;
  
  histoInclusive1->SetMarkerStyle(21);
  histoInclusive1->SetMarkerSize(1.2);
  histoInclusive1->SetMarkerColor(kRed+2);
  histoInclusive1->SetLineColor(kRed+2);
  
  TFile* fileInclusive2  = new TFile(InputFileName2.Data());
  TH1F*  histoInclusive2 = (TH1F*)fileInclusive2->Get("2040_tC52_0");
  if(!histoInclusive2) cout << "histoInclusive2 not found in fileInclusive2!!" << endl;
  
  histoInclusive2->SetMarkerStyle(20);
  histoInclusive2->SetMarkerSize(1.2);
  histoInclusive2->SetMarkerColor(kGreen+2);
  histoInclusive2->SetLineColor(kGreen+2);
  
  TFile* fileInclusive3  = new TFile(InputFileName3.Data());
  TH1F*  histoInclusive3 = (TH1F*)fileInclusive3->Get("2040_tC50_5");
  if(!histoInclusive3) cout << "histoInclusive3 not found in fileInclusive3!!" << endl;
  
  histoInclusive3->SetMarkerStyle(24);
  histoInclusive3->SetMarkerSize(1.2);
  histoInclusive3->SetMarkerColor(kBlue+2);
  histoInclusive3->SetLineColor(kBlue+2);

  //========================================================
  //histogram cosmetics
  //========================================================
  
  TH1F*  histoInclusiveEmpty = (TH1F*)Getv2HistStyle();
  
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
//   leg->SetHeader("v_{2}^{#gamma,inclusive}");
  leg->AddEntry(histoInclusive1,"Closed TPC","lp");
  leg->AddEntry(histoInclusive2,"Open TPC & -5<K<10","lp");
  leg->AddEntry(histoInclusive3,"Open TPC & 0<K<4","lp");
  leg->SetBorderSize(0);

  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoInclusive1->Draw("SAME");
  histoInclusive2->Draw("SAME");
  histoInclusive3->Draw("SAME");
  DrawInfoLabelOnPlot();
  
  leg->Draw();
  
  //Masking lower pt points
  MaskPoints(histoInclusive1,0,10);
  MaskPoints(histoInclusive2,0,10);
  MaskPoints(histoInclusive3,0,10);
  
  
  //Masking higher pt points
  MaskPoints(histoInclusive1,26,29);
  MaskPoints(histoInclusive2,26,29);
  MaskPoints(histoInclusive3,26,29);

  gSystem->mkdir("Results");
  gSystem->mkdir("Results/v2GammaInclusiveComparison");
  c1->SaveAs("Results/v2GammaInclusiveComparison/v2_Gamma_Inclusive_comparison.eps");
  
}