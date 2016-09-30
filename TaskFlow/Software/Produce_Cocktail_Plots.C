#include "DirectPhotonFlowFunctions.h"

void  Produce_Cocktail_Plots(
                              TString CentralityLow = "40",
                              TString CentralityHigh = "80"
){
  
  StyleSettingsThesis();
  SetPlotStyle();
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  //cen3 = 2040
  //cen6 = 020
  //cen8 = 4080
  Int_t centralityInt;
  if(CentralityLow.CompareTo("0")==0) centralityInt = 6;
  if(CentralityLow.CompareTo("20")==0) centralityInt = 3;
  if(CentralityLow.CompareTo("40")==0) centralityInt = 8;
  
  
  //open datafile with v2 gamma inclusive
  TFile* fileData = new TFile("/home/mike/0_directphoton/13_DirectPhotonAnalysisCode/Cocktail.root");
  
  TH1F*  hist1 = (TH1F*)fileData->Get(Form("hPi0Sp_cen%i",centralityInt));
  hist1->SetLineColor(kBlack);
  hist1->SetMarkerColor(kBlack);
  hist1->SetMarkerStyle(20);
  hist1->SetTitle("");
  hist1->GetXaxis()->SetTitle("#it{p}_{T}(GeV/#it{c})");
//   hist1->GetYaxis()->SetTitle("#frac{1}{2#pi p_{T}} #frac{d^{2}N}{dydp_{T}} (GeV/c)^{-2}");
  hist1->GetYaxis()->SetTitle("yield");
  hist1->GetXaxis()->SetTitleSize(0.05);
  hist1->GetYaxis()->SetTitleSize(0.05);
  hist1->GetYaxis()->SetTitleOffset(1.2);
  hist1->GetXaxis()->SetRangeUser(0.9,8);
  hist1->GetYaxis()->SetRangeUser(1E3,2E9);
  
  TH1F*  hist2 = (TH1F*)fileData->Get(Form("hEtaSp_cen%i",centralityInt));
  hist2->SetLineColor(kRed);
  hist2->SetMarkerColor(kRed);
  hist2->SetMarkerStyle(20);
  
  TH1F*  hist3 = (TH1F*)fileData->Get(Form("hK0sSp_cen%i",centralityInt));
  hist3->SetLineColor(kBlue);
  hist3->SetMarkerColor(kBlue);
  hist3->SetMarkerStyle(20);
  
  TH1F*  hist4 = (TH1F*)fileData->Get(Form("hOmegaSp_cen%i",centralityInt));
  hist4->SetLineColor(kMagenta);
  hist4->SetMarkerColor(kMagenta);
  hist4->SetMarkerStyle(20);
  
  TH1F*  hist5 = (TH1F*)fileData->Get(Form("hTotGammaSp_cen%i_Syst",centralityInt));
  hist5->SetLineColor(kGreen+2);
  hist5->SetMarkerColor(kGreen+2);
  hist5->SetMarkerStyle(21);
  hist5->SetTitle("");
  hist5->GetXaxis()->SetTitle("#it{p}_{T}(GeV/#it{c})");
  hist5->GetYaxis()->SetTitle("yield");
  hist5->GetXaxis()->SetTitleSize(0.05);
  hist5->GetYaxis()->SetTitleSize(0.05);
  hist5->GetYaxis()->SetTitleOffset(1.2);
  hist5->GetXaxis()->SetRangeUser(0.9,14.0);
  hist5->GetYaxis()->SetRangeUser(50,2E9);
  
  
  TLatex T1;
  T1.SetTextSize(0.04);
  T1.SetTextAlign(12);
  T1.SetNDC();
  
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  TLegend* leg1 = new TLegend(0.68,0.5,0.9,0.8);
  leg1->SetTextSize(0.05);
  leg1->AddEntry(hist1,"#pi^{0}","lp");
  leg1->AddEntry(hist2,"#eta","lp");
  leg1->AddEntry(hist3,"K^{0}_{S}","lp");
  leg1->AddEntry(hist4,"#omega","lp");
  leg1->AddEntry(hist5,"#gamma^{decay}","lp");
  leg1->SetBorderSize(0);
  
  TLegend* leg2 = new TLegend(0.68,0.5,0.9,0.8);
  leg2->SetTextSize(0.05);
  leg2->AddEntry(hist5,"#gamma^{decay}","lp");
  leg2->SetBorderSize(0);
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  SetProperMargins();
  gPad->SetRightMargin(0.05);
  c1->SetLogy();
  
  hist1->Draw();
  hist2->Draw("SAME");
  hist3->Draw("SAME");
  hist4->Draw("SAME");
  hist5->Draw("SAME");
  leg1->Draw();
  T1.DrawLatex(0.68, 0.92, Form("%s-%s%% PbPb",CentralityLow.Data(),CentralityHigh.Data()));
  T1.DrawLatex(0.68, 0.86, "#sqrt{s_{NN}}=2.76TeV");
  
  TCanvas* c2 = new TCanvas("c2","",800,800);
  SetProperMargins();
  gPad->SetRightMargin(0.05);
  c2->SetLogy();
  hist5->Draw();
  leg2->Draw();
  T1.DrawLatex(0.68, 0.92, Form("%s-%s%% PbPb",CentralityLow.Data(),CentralityHigh.Data()));
  T1.DrawLatex(0.68, 0.86, "#sqrt{s_{NN}}=2.76TeV");
  
  gSystem->mkdir("Results");
  gSystem->mkdir("Results/Cocktail");
  c1->SaveAs(Form("Results/Cocktail/inputspectra_%s_%s.eps",CentralityLow.Data(),CentralityHigh.Data()));
  c2->SaveAs(Form("Results/Cocktail/GammaTotal_%s_%s.eps",CentralityLow.Data(),CentralityHigh.Data()));

  
}






