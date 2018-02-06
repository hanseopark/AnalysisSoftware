#include "DirectPhotonFlowFunctions.h"

void  Plot_InclusivePhotonFlow(       
                                      TString CentralityLow = "20",
                                      TString CentralityHigh = "40"
){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
//   TString InputFileName1  = "/home/mike/3_PbPb_dirg/0_analysis/171030_v2_standardcuts/TaskFlow/Results_50200013_00200009007000008250400000/InclusivePhotonv2_Corrected_55_50200013_00200009007000008250400000.root";
//   TString InputFileName2  = "/home/mike/3_PbPb_dirg/0_analysis/171030_v2_standardcuts/TaskFlow/Results_50200013_00200009007002008750400000/InclusivePhotonv2_Corrected_58_50200013_00200009007002008750400000.root";
//   TString InputFileName3  = "/home/mike/3_PbPb_dirg/0_analysis/171030_v2_standardcuts/TaskFlow/Results_50200013_00200009007002009760600000/InclusivePhotonv2_Corrected_59_50200013_00200009007002009760600000.root";
  TString InputFileName1  = "/home/mike/3_PbPb_dirg/0_analysis/171030_v2_standardcuts/TaskFlow/Results_52400013_00200009007000008250400000/InclusivePhotonv2_Corrected_55_52400013_00200009007000008250400000.root";
  TString InputFileName2  = "/home/mike/3_PbPb_dirg/0_analysis/171030_v2_standardcuts/TaskFlow/Results_52400013_00200009007002008750400000/InclusivePhotonv2_Corrected_58_52400013_00200009007002008750400000.root";
  TString InputFileName3  = "/home/mike/3_PbPb_dirg/0_analysis/171030_v2_standardcuts/TaskFlow/Results_52400013_00200009007002009760600000/InclusivePhotonv2_Corrected_59_52400013_00200009007002009760600000.root";
  //open datafile with v2 gamma inclusive uncorrected
  TFile* fileInclusive1  = new TFile(InputFileName1.Data());
  TFile* fileInclusive2  = new TFile(InputFileName2.Data());
  TFile* fileInclusive3  = new TFile(InputFileName3.Data());
//   TH1F*  histoInclusive1 = (TH1F*)fileInclusive1->Get("v2GammaIncl_corrected_020_tC55");
  TH1F*  histoInclusive1 = (TH1F*)fileInclusive1->Get("v2GammaIncl_corrected_2040_tC55");
  if(!histoInclusive1) cout << "histoInclusive1 not found in fileInclusive1!!" << endl;
  
  histoInclusive1->SetMarkerStyle(21);
  histoInclusive1->SetMarkerSize(2);
  histoInclusive1->SetMarkerColor(kRed+2);
  histoInclusive1->SetLineColor(kRed+2);
  
//   TH1F*  histoInclusive2 = (TH1F*)fileInclusive2->Get("v2GammaIncl_corrected_020_tC58");
  TH1F*  histoInclusive2 = (TH1F*)fileInclusive2->Get("v2GammaIncl_corrected_2040_tC58");
  if(!histoInclusive2) cout << "histoInclusive2 not found in fileInclusive1!!" << endl;
  
  histoInclusive2->SetMarkerStyle(33);
  histoInclusive2->SetMarkerSize(2);
  histoInclusive2->SetMarkerColor(kGreen+2);
  histoInclusive2->SetLineColor(kGreen+2);
  
//   TH1F*  histoInclusive3 = (TH1F*)fileInclusive3->Get("v2GammaIncl_corrected_020_tC59");
  TH1F*  histoInclusive3 = (TH1F*)fileInclusive3->Get("v2GammaIncl_corrected_2040_tC59");
  if(!histoInclusive3) cout << "histoInclusive3 not found in fileInclusive1!!" << endl;
  
  histoInclusive3->SetMarkerStyle(34);
  histoInclusive3->SetMarkerSize(2);
  histoInclusive3->SetMarkerColor(kBlue+2);
  histoInclusive3->SetLineColor(kBlue+2);
  
  TH1F*  histoInclusiveRatio2 = (TH1F*)histoInclusive2->Clone();
  histoInclusiveRatio2->Divide(histoInclusive2,histoInclusive1,1,1,"B");
  
  TH1F*  histoInclusiveRatio3 = (TH1F*)histoInclusive3->Clone();
  histoInclusiveRatio3->Divide(histoInclusive3,histoInclusive1,1,1,"B");

  //========================================================
  //histogram cosmetics
  //========================================================
  
  TH1F*  histoInclusiveEmpty = (TH1F*)Getv2HistStyle();
  histoInclusiveEmpty->GetYaxis()->SetRangeUser(0,0.39);
  
  TH1F*  histoInclusiveRatioEmpty = (TH1F*)Getv2InclRatioHistStyle();
  histoInclusiveRatioEmpty->GetYaxis()->SetRangeUser(0.6,1.4);
  
  MaskPoints(histoInclusiveRatio2,0,10);
  MaskPoints(histoInclusiveRatio3,0,10);
  
  TLine *line = new TLine(0,1,6.2,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  
  TLine *line2 = new TLine(0,1.05,6.2,1.05);
  line2->SetLineStyle(3);
  line2->SetLineWidth(2);
  
  TLine *line3 = new TLine(0,0.95,6.2,0.95);
  line3->SetLineStyle(3);
  line3->SetLineWidth(2);
  
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
  
  TCanvas* c1 = new TCanvas("c1","",800,1600);
  c1->Divide(1,2,1E-10,1E-10);
  TLegend* leg = new TLegend(0.16,0.75,0.4,0.95);
  leg->SetTextSize(0.04);
  leg->AddEntry(histoInclusive1,"set 1","lp");
  leg->AddEntry(histoInclusive2,"set 2","lp");
  leg->AddEntry(histoInclusive3,"set 3","lp");
  leg->SetBorderSize(0);

  c1->cd(1);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoInclusive1->Draw("SAME");
  histoInclusive2->Draw("SAME");
  histoInclusive3->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  leg->Draw();
  
  TLegend* leg2 = new TLegend(0.16,0.75,0.4,0.95);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(histoInclusive2,"set 2 / set 1","lp");
  leg2->AddEntry(histoInclusive3,"set 3 / set 1","lp");
  leg2->SetBorderSize(0);
  
  c1->cd(2);
  SetProperMargins();
  histoInclusiveRatioEmpty->Draw();
  line->Draw();
  line2->Draw();
  line3->Draw();
  histoInclusiveRatio2->Draw("SAME");
  histoInclusiveRatio3->Draw("SAME");
  
  leg2->Draw();
  

  gSystem->mkdir("Results");
  gSystem->mkdir("Results/v2GammaInclusiveComparison");
  c1->SaveAs(Form("Results/v2GammaInclusiveComparison/%s_%s_v2_Gamma_Inclusive_comparison_2.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
}