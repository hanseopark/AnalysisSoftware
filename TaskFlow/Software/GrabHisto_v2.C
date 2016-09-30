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
#include "/home/mike/afterburner_v3/AnalysisSoftware/CommonHeaders/PlottingGammaConversionAdditional.h"
#include "/home/mike/afterburner_v3/AnalysisSoftware/CommonHeaders/PlottingGammaConversionHistos.h"

void MaskPoints(TH1F* histo, const int start, const int end);

void  GrabHisto_v2(){
  
  StyleSettingsThesis();
  SetPlotStyle();
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  //open datafile with v2 gamma inclusive
  TFile* fileData = new TFile("/home/mike/0_directphoton/0_analysis/160203_PbPb_v2_10h/InclusivePhotonv2_SmallCC.root");
  
  TH1F*  histo50 = (TH1F*)fileData->Get("2040_tC50");
  histo50->SetLineColor(kBlack);
  histo50->SetMarkerColor(kBlack);
  histo50->SetMarkerStyle(20);
  histo50->GetXaxis()->SetTitle("#it{p}_{T}(GeV/#it{c})");
  histo50->GetYaxis()->SetTitle("#it{v}_{2}  ");
  histo50->GetXaxis()->SetTitleSize(0.05);
  histo50->GetYaxis()->SetTitleSize(0.05);
  histo50->GetYaxis()->SetTitleOffset(1.2);
  histo50->GetXaxis()->SetRangeUser(0.0,8.0);
  histo50->GetYaxis()->SetRangeUser(0,0.46);
  
  TH1F*  histo51 = (TH1F*)fileData->Get("2040_tC51");
  histo51->SetLineColor(kRed+2);
  histo51->SetMarkerColor(kRed+2);
  histo51->SetMarkerStyle(20);
  TH1F*  histo52 = (TH1F*)fileData->Get("2040_tC52");
  histo52->SetLineColor(kGreen+2);
  histo52->SetMarkerColor(kGreen+2);
  histo52->SetMarkerStyle(24);
  TH1F*  histo53 = (TH1F*)fileData->Get("2040_tC53");
  histo53->SetLineColor(kBlue+2);
  histo53->SetMarkerColor(kBlue+2);
  histo53->SetMarkerStyle(24);
  TH1F*  histo54 = (TH1F*)fileData->Get("2040_tC54");
  histo54->SetLineColor(kYellow+2);
  histo54->SetMarkerColor(kYellow+2);
  histo54->SetMarkerStyle(24);
  TH1F*  histo55 = (TH1F*)fileData->Get("2040_tC55");
  histo55->SetLineColor(kMagenta+2);
  histo55->SetMarkerColor(kMagenta+2);
  histo55->SetMarkerStyle(24);
  
  TH1F*  histo56 = (TH1F*)fileData->Get("2040_tC56");
  histo56->SetLineColor(kBlack);
  histo56->SetMarkerColor(kBlack);
  histo56->SetMarkerStyle(20);
  histo56->GetXaxis()->SetTitle("#it{p}_{T}(GeV/#it{c})");
  histo56->GetYaxis()->SetTitle("#it{v}_{2}  ");
  histo56->GetXaxis()->SetTitleSize(0.05);
  histo56->GetYaxis()->SetTitleSize(0.05);
  histo56->GetYaxis()->SetTitleOffset(1.2);
  histo56->GetXaxis()->SetRangeUser(0.0,8.0);
  histo56->GetYaxis()->SetRangeUser(0,0.46);
  TH1F*  histo57 = (TH1F*)fileData->Get("2040_tC57");
  histo57->SetLineColor(kRed+2);
  histo57->SetMarkerColor(kRed+2);
  histo57->SetMarkerStyle(20);
  
  TLatex T1;
  T1.SetTextSize(0.04);
  T1.SetTextAlign(12);
  T1.SetNDC();
  
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  TLegend* leg = new TLegend(0.15,0.7,0.35,0.95);
  leg->SetTextSize(0.05);
  leg->AddEntry(histo50,"0.0 < #Kappa < 0.5","lp");
  leg->AddEntry(histo51,"0.5 < #Kappa < 2.0","lp");
  leg->AddEntry(histo52,"2.0 < #Kappa < 3.0","lp");
  leg->AddEntry(histo53,"3.0 < #Kappa < 4.0","lp");
  leg->AddEntry(histo54,"4.0 < #Kappa < 5.0","lp");
  leg->AddEntry(histo55,"5.0 < #Kappa < 8.0","lp");
  leg->SetBorderSize(0);
  
  TLegend* leg2 = new TLegend(0.15,0.85,0.35,0.95);
  leg2->SetTextSize(0.05);
  leg2->AddEntry(histo56,"0.0 < #Kappa < 2.0","lp");
  leg2->AddEntry(histo57,"2.0 < #Kappa < 8.0","lp");
  leg2->SetBorderSize(0);
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  c1->SetTopMargin(0.01);
  c1->SetRightMargin(0.01);
  
  histo50->Draw();
  histo51->Draw("SAME");
  histo52->Draw("SAME");
  histo53->Draw("SAME");
  histo54->Draw("SAME");
  histo55->Draw("SAME");
  leg->Draw();
  T1.DrawLatex(0.68, 0.92, "20-40% PbPb");
  T1.DrawLatex(0.68, 0.86, "#sqrt{s_{NN}}=2.76TeV");
  T1.DrawLatex(0.68, 0.80, "#gamma_{incl} #rightarrow e^{+}e^{-}");
  
  TCanvas* c2 = new TCanvas("c2","",800,800);
  c2->SetTopMargin(0.01);
  c2->SetRightMargin(0.01);
  
  histo56->Draw();
  histo57->Draw("SAME");
  leg2->Draw();
  T1.DrawLatex(0.68, 0.92, "20-40% PbPb");
  T1.DrawLatex(0.68, 0.86, "#sqrt{s_{NN}}=2.76TeV");
  T1.DrawLatex(0.68, 0.80, "#gamma_{incl} #rightarrow e^{+}e^{-}");
  
  MaskPoints(histo50,0,10);
  MaskPoints(histo51,0,10);
  MaskPoints(histo52,0,10);
  MaskPoints(histo53,0,10);
  MaskPoints(histo54,0,10);
  MaskPoints(histo55,0,10);
  MaskPoints(histo56,0,10);
  MaskPoints(histo57,0,10);
  
  MaskPoints(histo50,26,29);
  MaskPoints(histo51,26,29);
  MaskPoints(histo52,26,29);
  MaskPoints(histo53,26,29);
  MaskPoints(histo54,26,29);
  MaskPoints(histo55,26,29);
  MaskPoints(histo56,26,29);
  MaskPoints(histo57,26,29);
  
  c1->SaveAs("v2_gamma_inclusive_pt_1.eps");
  c2->SaveAs("v2_gamma_inclusive_pt_2.eps");
  
  TFile *OutputFile = new TFile("InclusivePhotonv2_KappaBinning.root","UPDATE");
  histo50->Write("Kappa_00_05");
  histo51->Write("Kappa_05_20");
  histo52->Write("Kappa_20_30");
  histo53->Write("Kappa_30_40");
  histo54->Write("Kappa_40_50");
  histo55->Write("Kappa_50_80");
  histo56->Write("Kappa_00_20");
  histo57->Write("Kappa_20_80");
  OutputFile->Close();
  
//   //get bin values
//   Double_t bins[] = {0,0.5,2,3,4,5,8};
//   TH1F* histoKappa1 = new TH1F("","",6,bins);
//   histoKappa1->SetLineColor(kBlack);
//   histoKappa1->SetMarkerColor(kBlack);
//   histoKappa1->SetMarkerStyle(20);
//   histoKappa1->GetXaxis()->SetTitle("#Kappa");
//   histoKappa1->GetYaxis()->SetTitle("#it{v}_{2}  ");
//   histoKappa1->GetXaxis()->SetTitleSize(0.05);
//   histoKappa1->GetYaxis()->SetTitleSize(0.05);
//   histoKappa1->GetYaxis()->SetTitleOffset(1.2);
//   histoKappa1->GetXaxis()->SetRangeUser(0.0,8.0);
//   histoKappa1->GetYaxis()->SetRangeUser(0,0.46);
//   
//   TH1F* histoKappa2 = new TH1F("","",6,bins);
//   histoKappa2->SetLineColor(kRed+2);
//   histoKappa2->SetMarkerColor(kRed+2);
//   histoKappa2->SetMarkerStyle(20);
//   TH1F* histoKappa3 = new TH1F("","",6,bins);
//   histoKappa3->SetLineColor(kGreen+2);
//   histoKappa3->SetMarkerColor(kGreen+2);
//   histoKappa3->SetMarkerStyle(20);
//   TH1F* histoKappa4 = new TH1F("","",6,bins);
//   histoKappa4->SetLineColor(kBlue+2);
//   histoKappa4->SetMarkerColor(kBlue+2);
//   histoKappa4->SetMarkerStyle(20);
//   TH1F* histoKappa5 = new TH1F("","",6,bins);
//   histoKappa5->SetLineColor(kYellow+2);
//   histoKappa5->SetMarkerColor(kYellow+2);
//   histoKappa5->SetMarkerStyle(20);
//   TH1F* histoKappa6 = new TH1F("","",6,bins);
//   histoKappa6->SetLineColor(kMagenta+2);
//   histoKappa6->SetMarkerColor(kMagenta+2);
//   histoKappa6->SetMarkerStyle(20);
//   
//   int pt_bin_1 = 10;
//   histoKappa1->SetBinContent(1,histo50->GetBinContent(pt_bin_1));
//   histoKappa1->SetBinContent(2,histo51->GetBinContent(pt_bin_1));
//   histoKappa1->SetBinContent(3,histo52->GetBinContent(pt_bin_1));
//   histoKappa1->SetBinContent(4,histo53->GetBinContent(pt_bin_1));
//   histoKappa1->SetBinContent(5,histo54->GetBinContent(pt_bin_1));
//   histoKappa1->SetBinContent(6,histo55->GetBinContent(pt_bin_1));
//   
//   int pt_bin_2 = 12;
//   histoKappa2->SetBinContent(1,histo50->GetBinContent(pt_bin_2));
//   histoKappa2->SetBinContent(2,histo51->GetBinContent(pt_bin_2));
//   histoKappa2->SetBinContent(3,histo52->GetBinContent(pt_bin_2));
//   histoKappa2->SetBinContent(4,histo53->GetBinContent(pt_bin_2));
//   histoKappa2->SetBinContent(5,histo54->GetBinContent(pt_bin_2));
//   histoKappa2->SetBinContent(6,histo55->GetBinContent(pt_bin_2));
//   
//   int pt_bin_3 = 14;
//   histoKappa3->SetBinContent(1,histo50->GetBinContent(pt_bin_3));
//   histoKappa3->SetBinContent(2,histo51->GetBinContent(pt_bin_3));
//   histoKappa3->SetBinContent(3,histo52->GetBinContent(pt_bin_3));
//   histoKappa3->SetBinContent(4,histo53->GetBinContent(pt_bin_3));
//   histoKappa3->SetBinContent(5,histo54->GetBinContent(pt_bin_3));
//   histoKappa3->SetBinContent(6,histo55->GetBinContent(pt_bin_3));
//   
//   int pt_bin_4 = 16;
//   histoKappa4->SetBinContent(1,histo50->GetBinContent(pt_bin_4));
//   histoKappa4->SetBinContent(2,histo51->GetBinContent(pt_bin_4));
//   histoKappa4->SetBinContent(3,histo52->GetBinContent(pt_bin_4));
//   histoKappa4->SetBinContent(4,histo53->GetBinContent(pt_bin_4));
//   histoKappa4->SetBinContent(5,histo54->GetBinContent(pt_bin_4));
//   histoKappa4->SetBinContent(6,histo55->GetBinContent(pt_bin_4));
//   
//   int pt_bin_5 = 18;
//   histoKappa5->SetBinContent(1,histo50->GetBinContent(pt_bin_5));
//   histoKappa5->SetBinContent(2,histo51->GetBinContent(pt_bin_5));
//   histoKappa5->SetBinContent(3,histo52->GetBinContent(pt_bin_5));
//   histoKappa5->SetBinContent(4,histo53->GetBinContent(pt_bin_5));
//   histoKappa5->SetBinContent(5,histo54->GetBinContent(pt_bin_5));
//   histoKappa5->SetBinContent(6,histo55->GetBinContent(pt_bin_5));
//   
//   int pt_bin_6 = 24;
//   histoKappa6->SetBinContent(1,histo50->GetBinContent(pt_bin_6));
//   histoKappa6->SetBinContent(2,histo51->GetBinContent(pt_bin_6));
//   histoKappa6->SetBinContent(3,histo52->GetBinContent(pt_bin_6));
//   histoKappa6->SetBinContent(4,histo53->GetBinContent(pt_bin_6));
//   histoKappa6->SetBinContent(5,histo54->GetBinContent(pt_bin_6));
//   histoKappa6->SetBinContent(6,histo55->GetBinContent(pt_bin_6));
//   
//   TCanvas* c3 = new TCanvas("c3","",800,800);
//   c3->SetTopMargin(0.01);
//   c3->SetRightMargin(0.01);
//   histoKappa1->Draw("p");
//   histoKappa2->Draw("pSAME");
//   histoKappa3->Draw("pSAME");
//   histoKappa4->Draw("pSAME");
//   histoKappa5->Draw("pSAME");
//   histoKappa6->Draw("pSAME");
  
  
  
}

//Function that masks points from bin# start to bin# end
void MaskPoints(TH1F* histo, const int start, const int end){
  for(int i=start;i<end;i++){
    histo->SetBinContent(i,-10);
  }
}











