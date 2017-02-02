#include "DirectPhotonFlowFunctions.h"

void  Plot_InclusivePhotonFlow(       
                                      TString CentralityLow = "0",
                                      TString CentralityHigh = "20",
                                      TString InputFileName1  = "/home/mike/3_PbPb_dirg/0_analysis/161017_v2_vs_kappa/Results_50200013_00200009007000008250400000/InclusivePhotonv2_uncorrected_50200013_00200009007000008250400000.root"
){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  //open datafile with v2 gamma inclusive uncorrected
  TFile* fileInclusive1  = new TFile(InputFileName1.Data());
  TH1F*  histoInclusive1 = (TH1F*)fileInclusive1->Get("50200013_00200009007000008250400000_tC58_5");
  if(!histoInclusive1) cout << "histoInclusive1 not found in fileInclusive1!!" << endl;
  
  histoInclusive1->SetMarkerStyle(21);
  histoInclusive1->SetMarkerSize(1.2);
  histoInclusive1->SetMarkerColor(kRed+2);
  histoInclusive1->SetLineColor(kRed+2);
  
  TH1F*  histoInclusive2 = (TH1F*)fileInclusive1->Get("50200013_00200009007000008250400000_tC58_6");
  if(!histoInclusive2) cout << "histoInclusive2 not found in fileInclusive1!!" << endl;
  
  histoInclusive2->SetMarkerStyle(33);
  histoInclusive2->SetMarkerSize(1.2);
  histoInclusive2->SetMarkerColor(kGreen+2);
  histoInclusive2->SetLineColor(kGreen+2);
  
  TH1F*  histoInclusive3 = (TH1F*)fileInclusive1->Get("50200013_00200009007000008250400000_tC58_3");
  if(!histoInclusive3) cout << "histoInclusive3 not found in fileInclusive1!!" << endl;
  
  histoInclusive3->SetMarkerStyle(34);
  histoInclusive3->SetMarkerSize(1.2);
  histoInclusive3->SetMarkerColor(kBlue+2);
  histoInclusive3->SetLineColor(kBlue+2);
  
  TString InputFileName2  = Form("/home/mike/git_afterburner/AnalysisSoftware/TaskFlow/Results/PCM_DirectPhotonFlow_%s%s.root",CentralityLow.Data(),CentralityHigh.Data());
  TFile* fileInclusive2  = new TFile(InputFileName2.Data());
  TH1F*  histoInclusive_cor = (TH1F*)fileInclusive2->Get("histoInclusive");

  //========================================================
  //histogram cosmetics
  //========================================================
  
  TH1F*  histoInclusiveEmpty = (TH1F*)Getv2HistStyle();
  histoInclusiveEmpty->GetYaxis()->SetRangeUser(0,0.49);
  
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
  
  TLegend* leg = new TLegend(0.16,0.75,0.4,0.95);
  leg->SetTextSize(0.04);
//   leg->SetHeader("v_{2}^{#gamma,inclusive}");
  leg->AddEntry(histoInclusive_cor,"#gamma_{inc}","lp");
//   leg->AddEntry(histoInclusive1,"-18<K<-16","lp");
//   leg->AddEntry(histoInclusive2,"-16<K<-14","lp");
  leg->AddEntry(histoInclusive1,"-10<K<-8","lp");
  leg->AddEntry(histoInclusive2,"-8<K<-6","lp");
//   leg->AddEntry(histoInclusive3,"-14<K<-12","lp");
  leg->SetBorderSize(0);

  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoInclusive_cor->Draw("SAME");
  histoInclusive1->Draw("SAME");
  histoInclusive2->Draw("SAME");
//   histoInclusive3->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
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
  c1->SaveAs(Form("Results/v2GammaInclusiveComparison/%s_%s_v2_Gamma_Inclusive_comparison_2.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
}