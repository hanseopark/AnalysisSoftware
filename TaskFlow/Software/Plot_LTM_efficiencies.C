#include "DirectPhotonFlowFunctions.h"

void  Plot_LTM_efficiencies(
){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  TFile* file1  = new TFile("Results/020LTM_efficiencies.root");
  TH1F*  h1 = (TH1F*)file1->Get("hLTM_Eff_pt");
  if(!h1) cout << "h1 not found in file1!!" << endl;
  
  h1->GetYaxis()->SetTitle("|#epsilon_{in} - #epsilon_{out}|");
  h1->GetYaxis()->SetRangeUser(0.0,0.249);
  
  h1->SetMarkerStyle(21);
  h1->SetMarkerSize(1.2);
  h1->SetMarkerColor(kRed+2);
  h1->SetLineColor(kRed+2);
  
  TFile* file2  = new TFile("Results/2040LTM_efficiencies.root");
  TH1F*  h2 = (TH1F*)file2->Get("hLTM_Eff_pt");
  if(!h1) cout << "h2 not found in file2!!" << endl;
  
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(1.2);
  h2->SetMarkerColor(kGreen+2);
  h2->SetLineColor(kGreen+2);
  
  TFile* file3  = new TFile("Results/4080LTM_efficiencies.root");
  TH1F*  h3 = (TH1F*)file3->Get("hLTM_Eff_pt");
  if(!h3) cout << "h3 not found in file3!!" << endl;
  
  h3->SetMarkerStyle(23);
  h3->SetMarkerSize(1.2);
  h3->SetMarkerColor(kBlue+2);
  h3->SetLineColor(kBlue+2);
  
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
  
  TLegend* leg = new TLegend(0.16,0.8,0.5,0.95);
  leg->SetTextSize(0.04);
  leg->AddEntry(h1,"0-20 %","lp");
  leg->AddEntry(h2,"20-40 %","lp");
  leg->AddEntry(h3,"40-80 %","lp");
  leg->SetBorderSize(0);

  SetProperMargins();
  h1->Draw("p");
  h2->Draw("pSAME");
  h3->Draw("pSAME");
  
  T1.DrawLatex(0.62, 0.94, "#bf{PbPb, #sqrt{s_{_{NN}}}=2.76TeV}");
  T1.DrawLatex(0.77, 0.88, "#bf{#gamma_{incl} #rightarrow e^{+}e^{-}}");
  
  leg->Draw();

  gSystem->mkdir("Results");
  c1->SaveAs("Results/LTM_efficiencies_allcents.eps");
  
}