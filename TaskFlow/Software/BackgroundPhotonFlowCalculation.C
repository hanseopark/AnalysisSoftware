#include "DirectPhotonFlowFunctions.h"

void  BackgroundPhotonFlowCalculation(
                                      Int_t Trainconfig = 55,
                                      TString CentralityLow = "0",
                                      TString CentralityHigh = "20",
                                      TString Cutnumber = "50200013_00200009007000008250400000"
){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  TString InputFileNameBkg = Form("/home/mike/3_PbPb_dirg/0_analysis/170214_v2_final_standardcuts/Results_%s/InclusivePhotonv2_uncorrected_%s.root",Cutnumber.Data(),Cutnumber.Data());
  TString InputFileNamePurity   = Form("/home/mike/3_PbPb_dirg/0_analysis/170214_v2_final_standardcuts/purity_studies_%s/Purity_InclusivePhotonSample_%s%s_%s_Data LHC10h.root",Cutnumber.Data(),CentralityLow.Data(),CentralityHigh.Data(),Cutnumber.Data());
  
  //open datafile with v2 background 1
  TFile* fileBkg  = new TFile(InputFileNameBkg.Data());
  TH1F*  histov2Bkg1 = (TH1F*)fileBkg->Get(Form("%s_tC%i_0",Cutnumber.Data(),Trainconfig));
  if(!histov2Bkg1) cout << "histov2Bkg1 not found in fileBkg1!!" << endl;
  histov2Bkg1->SetMarkerColor(kOrange-3);
  histov2Bkg1->SetLineColor(kOrange-3);
  histov2Bkg1->SetMarkerStyle(21);
  //open datafile with v2 background 2
  TH1F*  histov2Bkg2 = (TH1F*)fileBkg->Get(Form("%s_tC%i_1",Cutnumber.Data(),Trainconfig));
  if(!histov2Bkg2) cout << "histov2Bkg1 not found in fileBkg2!!" << endl;
  histov2Bkg2->SetMarkerColor(kGreen-3);
  histov2Bkg2->SetLineColor(kGreen-3);
  histov2Bkg2->SetMarkerStyle(21);
  //open datafile with v2 background 3
  TH1F*  histov2Bkg3 = (TH1F*)fileBkg->Get(Form("%s_tC%i_3",Cutnumber.Data(),Trainconfig));
  if(!histov2Bkg3) cout << "histov2Bkg1 not found in fileBkg3!!" << endl;
  histov2Bkg3->SetMarkerColor(kAzure);
  histov2Bkg3->SetLineColor(kAzure);
  histov2Bkg3->SetMarkerStyle(21);
  //open datafile with purities
  TFile* filePurity  = new TFile(InputFileNamePurity.Data());
  TGraphErrors*  GraphPurityBkg1NA   = (TGraphErrors*)filePurity->Get(Form("bkg1_Purity_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  TGraphErrors*  GraphPurityBkg2NA = (TGraphErrors*)filePurity->Get(Form("bkg2_NA_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  TGraphErrors*  GraphPurityBkg2NB = (TGraphErrors*)filePurity->Get(Form("bkg2_NB_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  TGraphErrors*  GraphPurityBkg2NC = (TGraphErrors*)filePurity->Get(Form("bkg2_NC_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  TGraphErrors*  GraphPuritySignalNA = (TGraphErrors*)filePurity->Get(Form("Gamma_NA_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  TGraphErrors*  GraphPuritySignalNB = (TGraphErrors*)filePurity->Get(Form("Gamma_NB_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  TGraphErrors*  GraphPuritySignalNC = (TGraphErrors*)filePurity->Get(Form("Gamma_NC_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  if(!GraphPurityBkg1NA) cout << "GraphPurityBkg1NA not found in filePurity!!" << endl;
  if(!GraphPurityBkg2NA) cout << "GraphPurityBkg2NA not found in filePurity!!" << endl;
  if(!GraphPurityBkg2NB) cout << "GraphPurityBkg2NB not found in filePurity!!" << endl;
  if(!GraphPurityBkg2NC) cout << "GraphPurityBkg2NC not found in filePurity!!" << endl;
  if(!GraphPuritySignalNA) cout << "GraphPuritySignalNA not found in filePurity!!" << endl;
  if(!GraphPuritySignalNB) cout << "GraphPuritySignalNB not found in filePurity!!" << endl;
  if(!GraphPuritySignalNC) cout << "GraphPuritySignalNC not found in filePurity!!" << endl;
  
  GraphPurityBkg1NA->SetMarkerColor(kOrange-2);
  GraphPurityBkg1NA->SetLineColor(kOrange-2);
  GraphPurityBkg1NA->SetMarkerStyle(21);
  
  GraphPurityBkg2NA->SetMarkerColor(kGreen);
  GraphPurityBkg2NA->SetLineColor(kGreen);
  GraphPurityBkg2NA->SetMarkerStyle(21);
  
  GraphPurityBkg2NB->SetMarkerColor(kGreen+3);
  GraphPurityBkg2NB->SetLineColor(kGreen+3);
  GraphPurityBkg2NB->SetMarkerStyle(21);
  
  GraphPurityBkg2NC->SetMarkerColor(kGreen-3);
  GraphPurityBkg2NC->SetLineColor(kGreen-3);
  GraphPurityBkg2NC->SetMarkerStyle(21);
  
  GraphPuritySignalNA->SetMarkerColor(kBlue+2);
  GraphPuritySignalNA->SetLineColor(kBlue+2);
  GraphPuritySignalNA->SetMarkerStyle(21);
  
  GraphPuritySignalNB->SetMarkerColor(kMagenta+2);
  GraphPuritySignalNB->SetLineColor(kMagenta+2);
  GraphPuritySignalNB->SetMarkerStyle(21);
  
  GraphPuritySignalNC->SetMarkerColor(kRed+2);
  GraphPuritySignalNC->SetLineColor(kRed+2);
  GraphPuritySignalNC->SetMarkerStyle(21);
  
  
  //========================================================
  //histogram cosmetics
  //========================================================
  
  TH1F*  histov2Empty = (TH1F*)Getv2HistStyle();
  histov2Empty->GetYaxis()->SetRangeUser(0,0.49);
  
  TH1F*  histoPurityEmpty = (TH1F*)GetNHistStyle();
  histoPurityEmpty->GetYaxis()->SetRangeUser(0,1.19);
  
  TH1F*  histov2C = (TH1F*)histov2Bkg3->Clone();
  histov2C->SetLineColor(kRed+2);
  histov2C->SetMarkerColor(kRed+2);
  histov2C->SetMarkerStyle(21);
  
  cout << "//========================================================" << endl;
  cout << "//Calculating v2 of background region 1" << endl;
  cout << "//========================================================" << endl;
  
  TH1F*  histov2A = (TH1F*)histov2C->Clone();
  histov2A->SetLineColor(kBlue+2);
  histov2A->SetMarkerColor(kBlue+2);
  histov2A->SetMarkerStyle(21);
  
  Int_t NBins = histov2Bkg1->GetSize();
  for(Int_t i=0;i<NBins;i++){
    Float_t BinCenter = histov2Bkg1->GetBinCenter(i);
    Float_t v2TotalValue = histov2Bkg1->GetBinContent(i);
    Float_t Purity = GraphPurityBkg1NA->Eval(BinCenter);
    Float_t v2BkgValue = histov2C->Interpolate(BinCenter);
    Float_t v2CorrectedValue = ( v2TotalValue - (1-Purity) * v2BkgValue ) / Purity;
    
    histov2A->SetBinContent(i,v2CorrectedValue);
    
    cout << endl;
    cout << "Correcting pt = " << BinCenter << " GeV/c" << endl;
    cout << "Purity = " << Purity << endl;
    cout << "Uncorrected v2 value = " << v2TotalValue << endl;
    cout << "Backgound   v2 value = " << v2BkgValue << endl;
    cout << "Corrected   v2 value      = " << v2CorrectedValue << endl;
    cout << "Difference -> " << 100*(v2CorrectedValue-v2TotalValue)/v2TotalValue << " % " << endl;
    cout << endl;
  }

  cout << "//========================================================" << endl;
  cout << "//Calculating v2 of background region 2" << endl;
  cout << "//========================================================" << endl;
  
  TH1F*  histov2B = (TH1F*)histov2C->Clone();
  histov2B->SetLineColor(kMagenta+2);
  histov2B->SetMarkerColor(kMagenta+2);
  histov2B->SetMarkerStyle(21);
  
  for(Int_t i=0;i<NBins;i++){
    Float_t BinCenter = histov2Bkg1->GetBinCenter(i);
    Float_t v2TotalValue = histov2Bkg1->GetBinContent(i);
    Float_t NA = GraphPurityBkg2NA->Eval(BinCenter);
    Float_t NB = GraphPurityBkg2NB->Eval(BinCenter);
    Float_t NC = GraphPurityBkg2NC->Eval(BinCenter);
    Float_t v2BkgA = histov2A->Interpolate(BinCenter);
    Float_t v2BkgC = histov2C->Interpolate(BinCenter);
    
    Float_t v2CorrectedValue = (v2TotalValue-NA*v2BkgA-NC*v2BkgC) / NB;
    
    histov2B->SetBinContent(i,v2CorrectedValue);
    
    cout << endl;
    cout << "Correcting pt = " << BinCenter << " GeV/c" << endl;
    cout << "NA = " << NA << endl;
    cout << "NB = " << NB << endl;
    cout << "NC = " << NC << endl;
    cout << "Uncorrected v2 value = " << v2TotalValue << endl;
    cout << "BkgA   v2 value = " << v2BkgA << endl;
    cout << "BkgC   v2 value = " << v2BkgC << endl;
    cout << "BkgB   v2 value = " << v2CorrectedValue << endl;
    cout << "Difference -> " << 100*(v2CorrectedValue-v2TotalValue)/v2TotalValue << " % " << endl;
    cout << endl;
  }
  
  cout << "//========================================================" << endl;
  cout << "//Calculating v2 of background in signal region" << endl;
  cout << "//========================================================" << endl;
  
  TH1F*  histov2BkgRegionSignal = (TH1F*)histov2Bkg1->Clone();
  histov2BkgRegionSignal->SetLineColor(kBlack);
  histov2BkgRegionSignal->SetMarkerColor(kBlack);
  
  for(Int_t i=0;i<NBins;i++){
    Float_t BinCenter = histov2Bkg2->GetBinCenter(i);
    Float_t NA = GraphPuritySignalNA->Eval(BinCenter);
    Float_t NB = GraphPuritySignalNB->Eval(BinCenter);
    Float_t NC = GraphPuritySignalNC->Eval(BinCenter);
    Float_t v2BkgA = histov2A->Interpolate(BinCenter);
    Float_t v2BkgB = histov2B->Interpolate(BinCenter);
    Float_t v2BkgC = histov2C->Interpolate(BinCenter);
    
    Float_t v2BkgTotal = NA*v2BkgA + NB*v2BkgB + NC*v2BkgC;
    
    histov2BkgRegionSignal->SetBinContent(i,v2BkgTotal);
    
  }
  
  TGraphErrors* GraphPurityBkg1NB = (TGraphErrors*)GraphPurityBkg1NA->Clone();
  GraphPurityBkg1NA->SetMarkerColor(kOrange+1);
  GraphPurityBkg1NA->SetLineColor(kOrange+1);
  for(Int_t i=0;i<GraphPurityBkg1NB->GetN();i++){
    GraphPurityBkg1NB->SetPoint(i,GraphPurityBkg1NB->GetX()[i],1-GraphPurityBkg1NB->GetY()[i]);
  }
  
  MaskPoints(histov2A,0,10);
  MaskPoints(histov2B,0,10);
  MaskPoints(histov2C,0,10);
  MaskPoints(histov2Bkg1,0,10);
  MaskPoints(histov2Bkg2,0,10);
  MaskPoints(histov2Bkg3,0,10);
  MaskPoints(histov2BkgRegionSignal,0,10);
  
  //========================================================
  //Plotting and Saving
  //========================================================
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  TLegend* leg1 = new TLegend(0.16,0.8,0.6,0.95);
  leg1->SetTextSize(0.04);
  leg1->AddEntry(histov2Bkg1,"v_{2} region 1","lp");
  leg1->AddEntry(histov2Bkg2,"v_{2} region 2","lp");
  leg1->AddEntry(histov2Bkg3,"v_{2} region 3","lp");
  leg1->SetBorderSize(0);
  
  TLegend* leg2 = new TLegend(0.16,0.7,0.6,0.95);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(histov2A,"v_{2}^{#pi#pi}","lp");
  leg2->AddEntry(histov2B,"v_{2}^{#pie}","lp");
  leg2->AddEntry(histov2C,"v_{2}^{rem}","lp");
  leg2->AddEntry(histov2BkgRegionSignal,"v_{2}^{background}","lp");
  leg2->SetBorderSize(0);
  
  TLegend* leg3 = new TLegend(0.65,0.7,0.95,0.95);
  leg3->SetTextSize(0.04);
  leg3->AddEntry(GraphPurityBkg1NA,"Region 1: n_{#pi#pi}","lp");
  leg3->AddEntry(GraphPurityBkg1NB,"Region 1: n_{rem}","lp");
  leg3->SetBorderSize(0);
  
  TLegend* leg4 = new TLegend(0.65,0.7,0.95,0.95);
  leg4->SetTextSize(0.04);
  leg4->AddEntry(GraphPurityBkg2NA,"Region 2: n_{#pi#pi}","lp");
  leg4->AddEntry(GraphPurityBkg2NB,"Region 2: n_{#pie}","lp");
  leg4->AddEntry(GraphPurityBkg2NC,"Region 2: n_{rem}","lp");
  leg4->SetBorderSize(0);
  
  TLegend* leg5 = new TLegend(0.65,0.7,0.95,0.95);
  leg5->SetTextSize(0.04);
  leg5->AddEntry(GraphPuritySignalNA,"Signal n_{#pi#pi}","lp");
  leg5->AddEntry(GraphPuritySignalNB,"Signal n_{#pie}","lp");
  leg5->AddEntry(GraphPuritySignalNC,"Signal n_{rem}","lp");
  leg5->SetBorderSize(0);
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  SetProperMargins();
  histov2Empty->Draw();
  histov2Bkg1->Draw("same");
  histov2Bkg2->Draw("same");
  histov2Bkg3->Draw("same");
  leg1->Draw();
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  TCanvas* c2 = new TCanvas("c2","",800,800);
  SetProperMargins();
  histov2Empty->Draw();
  histov2Empty->GetYaxis()->SetRangeUser(0,0.599);
  histov2A->Draw("same");
  histov2B->Draw("same");
  histov2C->Draw("same");
  histov2BkgRegionSignal->Draw("same");
  leg2->Draw();
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  TCanvas* c3 = new TCanvas("c3","",800,800);
  SetProperMargins();
  histoPurityEmpty->Draw();
  GraphPurityBkg1NA->Draw("same PL");
  GraphPurityBkg1NB->Draw("same PL");
  leg3->Draw();
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),1);
  
  TCanvas* c4 = new TCanvas("c4","",800,800);
  SetProperMargins();
  histoPurityEmpty->Draw();
  GraphPurityBkg2NA->Draw("same PL");
  GraphPurityBkg2NB->Draw("same PL");
  GraphPurityBkg2NC->Draw("same PL");
  leg4->Draw();
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),1);
  
  TCanvas* c5 = new TCanvas("c5","",800,800);
  SetProperMargins();
  histoPurityEmpty->Draw();
  GraphPuritySignalNA->Draw("same PL");
  GraphPuritySignalNB->Draw("same PL");
  GraphPuritySignalNC->Draw("same PL");
  leg5->Draw();
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),1);
  
  gSystem->mkdir(Form("Results_%s",Cutnumber.Data()));
  gSystem->mkdir(Form("Results_%s/v2GammaBackground",Cutnumber.Data()));
  c1->SaveAs(Form("Results_%s/v2GammaBackground/GammaBackground_v2_region_%s%s.eps",Cutnumber.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  c2->SaveAs(Form("Results_%s/v2GammaBackground/GammaBackground_v2_bkg_%s%s.eps",Cutnumber.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  c3->SaveAs(Form("Results_%s/v2GammaBackground/GammaBackground_N_bkg1_%s%s.eps",Cutnumber.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  c4->SaveAs(Form("Results_%s/v2GammaBackground/GammaBackground_N_bkg2_%s%s.eps",Cutnumber.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  c5->SaveAs(Form("Results_%s/v2GammaBackground/GammaBackground_N_signal_%s%s.eps",Cutnumber.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  
  TFile *BackgroundPhotonFlow_File = new TFile(Form("Results_%s/v2_background_%s%s_%s.root",Cutnumber.Data(),CentralityLow.Data(),CentralityHigh.Data(),Cutnumber.Data()),"RECREATE");
  histov2A->Write(Form("histo_v2_bkg_pipi_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  histov2B->Write(Form("histo_v2_bkg_pie_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  histov2C->Write(Form("histo_v2_bkg_rem_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  histov2BkgRegionSignal->Write(Form("histo_v2_bkg_total_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  BackgroundPhotonFlow_File->Close();
  
}