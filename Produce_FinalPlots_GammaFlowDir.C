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
#include "TPaveText.h"
#include "TRandom3.h"
#include "Riostream.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"


void SetProperMargins();
TH1F* Getv2HistStyle();
TH1F* GetRatioHistStyle();
void DrawDirectInfoLabelOnPlot(TString cent1, TString cent2, Int_t place);
void EnlargeBinErrorsInX(TGraphAsymmErrors* g1);

void  Produce_FinalPlots_GammaFlowDir(
                                      TString CentralityLow = "20",
                                      TString CentralityHigh = "40",
                                      TString InputFileName1  = "TaskFlow/V2dir_calculation_PCM_PHOS_combined/output/v2dir_pcm_phos_comb_20-40.root"
){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  //open datafile with v2 gamma inclusive uncorrected
  TFile* fileInput  = new TFile(InputFileName1.Data());
  if(!fileInput) cout << "file not found!!" << endl;
  TGraphAsymmErrors*  graph_v2dir_comb_stat = (TGraphAsymmErrors*)fileInput->Get("g_v2_dir_comb_staterr");
  TGraphAsymmErrors*  graph_v2dir_comb_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_dir_comb_toterr");
  
  TGraphAsymmErrors*  graph_v2inc_comb_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_inc_comb_toterr");
  TGraphAsymmErrors*  graph_v2inc_pcm_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_inc_pcm_toterr");
  TGraphAsymmErrors*  graph_v2inc_phos_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_inc_phos_toterr");
  TGraphAsymmErrors*  graph_v2dec_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_dec_comb");
  
  if(!graph_v2dir_comb_stat) cout << "g_v2_dir_comb_toterr not found in fileInput!!" << endl;

  //========================================================
  //cosmetics
  //========================================================
  
  graph_v2inc_comb_tot->SetMarkerStyle(20);
  graph_v2inc_comb_tot->SetMarkerSize(1.2);
  graph_v2inc_comb_tot->SetMarkerColor(kBlack);
  graph_v2inc_comb_tot->SetLineColor(kBlack);
  
  graph_v2inc_pcm_tot->SetMarkerStyle(24);
  graph_v2inc_pcm_tot->SetMarkerSize(1.2);
  graph_v2inc_pcm_tot->SetMarkerColor(kRed+2);
  graph_v2inc_pcm_tot->SetLineColor(kRed+2);
  
  graph_v2inc_phos_tot->SetMarkerStyle(24);
  graph_v2inc_phos_tot->SetMarkerSize(1.2);
  graph_v2inc_phos_tot->SetMarkerColor(kBlue+2);
  graph_v2inc_phos_tot->SetLineColor(kBlue+2);
  
  graph_v2dec_tot->SetMarkerStyle(21);
  graph_v2dec_tot->SetMarkerSize(1.2);
  graph_v2dec_tot->SetMarkerColor(kGreen+2);
  graph_v2dec_tot->SetLineColor(kGreen+2);
  
  graph_v2dir_comb_stat->SetMarkerStyle(21);
  graph_v2dir_comb_stat->SetMarkerSize(1.2);
  graph_v2dir_comb_stat->SetMarkerColor(kBlack);
  graph_v2dir_comb_stat->SetLineColor(kBlack);
  
  graph_v2dir_comb_tot->SetFillColor(kGray);
  graph_v2dir_comb_tot->SetFillStyle(0);
  
  EnlargeBinErrorsInX(graph_v2inc_comb_tot);
  EnlargeBinErrorsInX(graph_v2inc_pcm_tot);
  EnlargeBinErrorsInX(graph_v2inc_phos_tot);
  EnlargeBinErrorsInX(graph_v2dec_tot);
  EnlargeBinErrorsInX(graph_v2dir_comb_tot);
  
  TH1F*  graph_v2dir_comb_statEmpty = (TH1F*)Getv2HistStyle();
  graph_v2dir_comb_statEmpty->GetYaxis()->SetRangeUser(-0.29,0.49);
  
  TH1F*  graph_v2Inc_comb_statEmpty = (TH1F*)Getv2HistStyle();
  graph_v2Inc_comb_statEmpty->GetYaxis()->SetRangeUser(0.01,0.249);
  
  TH1F*  graph_v2Inc_comb_Ratio_statEmpty = (TH1F*)GetRatioHistStyle();
  graph_v2Inc_comb_Ratio_statEmpty->GetYaxis()->SetRangeUser(0.5,1.5);
  
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
  
  gSystem->mkdir("PaperPlots");
  
  //========================================================
  //Plotting and Saving - INCLUSIVE
  //========================================================
  
  //ALL
  TCanvas* cIncComb = new TCanvas("cIncComb","",800,800);
  
  TLegend* legIncComb = new TLegend(0.16,0.15,0.6,0.36);
  legIncComb->SetTextSize(0.05);
  legIncComb->AddEntry(graph_v2inc_comb_tot,"Combined","lp");
  legIncComb->AddEntry(graph_v2dec_tot,"Cocktail","lp");
  legIncComb->AddEntry(graph_v2inc_pcm_tot,"PCM","lp");
  legIncComb->AddEntry(graph_v2inc_phos_tot,"PHOS","lp");
  
  legIncComb->SetBorderSize(0);

  SetProperMargins();
  graph_v2Inc_comb_statEmpty->Draw();
  graph_v2inc_phos_tot->Draw("SAME P");
  graph_v2inc_pcm_tot->Draw("SAME P");
  graph_v2dec_tot->Draw("SAME P");
  graph_v2inc_comb_tot->Draw("SAME P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),4);
  
  legIncComb->Draw();

  cIncComb->SaveAs(Form("PaperPlots/v2inc_combined_all_%s%s.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //ONLY COMB AND DEC
  TCanvas* cIncComb2 = new TCanvas("cIncComb","",800,800);
  
  TLegend* legIncComb2 = new TLegend(0.16,0.15,0.6,0.26);
  legIncComb2->SetTextSize(0.05);
  legIncComb2->AddEntry(graph_v2inc_comb_tot,"Combined","lp");
  legIncComb2->AddEntry(graph_v2dec_tot,"Cocktail","lp");
  
  legIncComb2->SetBorderSize(0);

  SetProperMargins();
  graph_v2Inc_comb_statEmpty->Draw();
  graph_v2dec_tot->Draw("SAME P");
  graph_v2inc_comb_tot->Draw("SAME P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),4);
  
  legIncComb2->Draw();

  cIncComb2->SaveAs(Form("PaperPlots/v2inc_combined_OnlyCombDec_%s%s.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //ONLY COMB AND DEC RATIO
  TCanvas* cIncComb2Ratio = new TCanvas("cIncComb2Ratio","",800,800);
  
  TLegend* legIncComb2Ratio = new TLegend(0.16,0.15,0.6,0.26);
  legIncComb2Ratio->SetTextSize(0.05);
  legIncComb2Ratio->AddEntry(graph_v2inc_comb_tot,"Combined","lp");
  legIncComb2Ratio->AddEntry(graph_v2dec_tot,"Cocktail","lp");
  
  legIncComb2Ratio->SetBorderSize(0);

  SetProperMargins();
  graph_v2Inc_comb_Ratio_statEmpty->Draw();
//   graph_v2dec_tot->Draw("SAME P");
//   graph_v2inc_comb_tot->Draw("SAME P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),4);
  
  legIncComb2Ratio->Draw();

  cIncComb2Ratio->SaveAs(Form("PaperPlots/v2inc_combined_OnlyCombDec_Ratio_%s%s.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //========================================================
  //Plotting and Saving - DIRECT
  //========================================================
  
  TCanvas* cDirComb = new TCanvas("cDirComb","",800,800);
  
  TLegend* legDirComb = new TLegend(0.16,0.15,0.6,0.26);
  legDirComb->SetTextSize(0.05);
  legDirComb->AddEntry(graph_v2dir_comb_stat,"statistical","lp");
  legDirComb->AddEntry(graph_v2dir_comb_tot,"total","f");
  legDirComb->SetBorderSize(0);

  SetProperMargins();
  graph_v2dir_comb_statEmpty->Draw();
  graph_v2dir_comb_stat->Draw("SAME P");
  graph_v2dir_comb_tot->Draw("SAME e2 P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  legDirComb->Draw();

  cDirComb->SaveAs(Form("PaperPlots/v2dir_combined_%s%s.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
}

void SetProperMargins(){
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.01);
  gPad->SetBottomMargin(0.11);
  gPad->SetTopMargin(0.01);
}

TH1F* Getv2HistStyle(){
 
  TH1F*  histoEmpty = new TH1F("","",100,0,10);
  
  for(Int_t i=0;i<100;i++){
    histoEmpty->SetBinContent(i,-100.);
  }
  
  
  histoEmpty->GetXaxis()->SetTitle("#it{p}_{T}(GeV/c)");
  histoEmpty->GetYaxis()->SetTitle("v_{2}");
  histoEmpty->GetXaxis()->SetTitleSize(0.06);
  histoEmpty->GetYaxis()->SetTitleSize(0.07);
  histoEmpty->GetXaxis()->SetTitleOffset(0.8);
  histoEmpty->GetYaxis()->SetTitleOffset(0.85);
  histoEmpty->GetXaxis()->SetRangeUser(0.0,6.99);
  histoEmpty->GetYaxis()->SetRangeUser(0,0.249);
  
  return histoEmpty;
}

TH1F* GetRatioHistStyle(){
 
  TH1F*  histoEmpty = new TH1F("","",100,0,8);
  
  for(Int_t i=0;i<100;i++){
    histoEmpty->SetBinContent(i,-100.);
  }
  
  
  histoEmpty->GetXaxis()->SetTitle("#it{p}_{T}(GeV/c)");
  histoEmpty->GetYaxis()->SetTitle("ratio");
  histoEmpty->GetXaxis()->SetTitleSize(0.06);
  histoEmpty->GetYaxis()->SetTitleSize(0.07);
  histoEmpty->GetXaxis()->SetTitleOffset(0.8);
  histoEmpty->GetYaxis()->SetTitleOffset(0.85);
  histoEmpty->GetXaxis()->SetRangeUser(0.0,6.99);
  histoEmpty->GetYaxis()->SetRangeUser(0,0.249);
  
  return histoEmpty;
}

void DrawDirectInfoLabelOnPlot(TString cent1, TString cent2, Int_t place){
  
  TLatex T1;
  T1.SetTextSize(0.04);
  T1.SetTextAlign(12);
  T1.SetNDC();
  
  if(place==1){
    T1.DrawLatex(0.15, 0.94, Form("#bf{%s-%s %% PbPb, #sqrt{s_{_{NN}}}=2.76TeV}",cent1.Data(),cent2.Data()));
    T1.DrawLatex(0.15, 0.88, "#bf{#gamma_{direct}}");
  }
  
  if(place==2){
    T1.DrawLatex(0.48, 0.94, Form("#bf{%s-%s %% PbPb, #sqrt{s_{_{NN}}}=2.76TeV}",cent1.Data(),cent2.Data()));
    T1.DrawLatex(0.87, 0.88, "#bf{#gamma_{direct}}");
  }
  
  if(place==3){
    T1.DrawLatex(0.15, 0.94, Form("#bf{%s-%s %% PbPb, #sqrt{s_{_{NN}}}=2.76TeV}",cent1.Data(),cent2.Data()));
    T1.DrawLatex(0.15, 0.88, "#bf{#gamma_{inclusive}}");
  }
  
  if(place==4){
    T1.DrawLatex(0.48, 0.94, Form("#bf{%s-%s %% PbPb, #sqrt{s_{_{NN}}}=2.76TeV}",cent1.Data(),cent2.Data()));
    T1.DrawLatex(0.83, 0.88, "#bf{#gamma_{inclusive}}");
  }
  
}

void EnlargeBinErrorsInX(TGraphAsymmErrors* g1){
  
  Int_t NPoints = g1->GetN();
  Double_t* g1_x = g1->GetX();
  Double_t* g1_y = g1->GetY();
 
  Double_t BinEdges[18] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.3, 3.7, 4.1, 4.6, 5.4, 6.2, 7.0};
  
  for(Int_t ibin=0;ibin<NPoints;ibin++){
    g1->SetPointEXlow(ibin,g1_x[ibin]-BinEdges[ibin]); g1->SetPointEXhigh(ibin,BinEdges[ibin+1]-g1_x[ibin]);
  }
  
}