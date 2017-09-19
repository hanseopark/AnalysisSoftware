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
TH1F* GetExampleBinStyle();
void DrawDirectInfoLabelOnPlot(TString cent1, TString cent2, Int_t place);
void EnlargeBinErrorsInX(TGraphAsymmErrors* g1);
TGraphAsymmErrors* calcGraphRatio(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2);
void CalculatepValue(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, Double_t startpT, Double_t endpT);

void  Produce_FinalPlots_GammaFlowDir(
                                      TString CentralityLow = "20",
                                      TString CentralityHigh = "40",
                                      Bool_t IncludeTheory = kTRUE
){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================

  TString centStringInputfile = "";
  if(CentralityLow.CompareTo("0")==0) centStringInputfile  = "00-20";
  if(CentralityLow.CompareTo("20")==0) centStringInputfile = "20-40";
  TString InputFileName1  = Form("TaskFlow/V2dir_calculation_PCM_PHOS_combined/output/v2dir_pcm_phos_comb_%s.root",centStringInputfile.Data());
  
  TString InputFileNameCocktail  = "TaskFlow/Results/CocktailV2_170912.root";
  
  //open datafile with v2 gamma inclusive uncorrected
  TFile* fileInput  = new TFile(InputFileName1.Data());
  if(!fileInput) cout << "file not found!!" << endl;
  TGraphAsymmErrors*  graph_v2dir_comb_stat = (TGraphAsymmErrors*)fileInput->Get("g_v2_dir_comb_staterr");
  TGraphAsymmErrors*  graph_v2dir_comb_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_dir_comb_toterr");
  
  TGraphAsymmErrors*  graph_v2dir_pcm_Rpcm_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_dir_pcm_Rpcm_toterr");
  TGraphAsymmErrors*  graph_v2dir_pcm_Rcomb_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_dir_pcm_Rcomb_toterr");
  TGraphAsymmErrors*  graph_v2dir_phos_Rphos_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_dir_phos_Rphos_toterr");
  TGraphAsymmErrors*  graph_v2dir_phos_Rcomb_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_dir_phos_Rcomb_toterr");
  
  TGraphAsymmErrors*  graph_v2inc_comb_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_inc_comb_toterr");
  TGraphAsymmErrors*  graph_v2inc_pcm_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_inc_pcm_toterr");
  TGraphAsymmErrors*  graph_v2inc_phos_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_inc_phos_toterr");
  TGraphAsymmErrors*  graph_v2dec_tot = (TGraphAsymmErrors*)fileInput->Get("g_v2_dec_comb");
  
  if(!graph_v2dir_comb_stat) cout << "g_v2_dir_comb_toterr not found in fileInput!!" << endl;
  
  //open datafile with v2 cocktail
  TFile* InputFileCocktail  = new TFile(InputFileNameCocktail.Data());
  if(!InputFileCocktail) cout << "file not found!!" << endl;
  
  TString centStringCocktail = "";
  if(CentralityLow.CompareTo("0")==0) centStringCocktail  = "_cen6";
  if(CentralityLow.CompareTo("20")==0) centStringCocktail = "_cen3";
  
  TH1D*  hist_v2gammaCocktail = (TH1D*)InputFileCocktail->Get(Form("v2gammaCocktail%s",centStringCocktail.Data()));
  TH1D*  hist_v2gammaCocktailPi0 = (TH1D*)InputFileCocktail->Get(Form("v2gammaPi0%s",centStringCocktail.Data()));
  TH1D*  hist_v2gammaCocktailK0s = (TH1D*)InputFileCocktail->Get(Form("v2gammaK0s%s",centStringCocktail.Data()));
  TH1D*  hist_v2gammaCocktailEta = (TH1D*)InputFileCocktail->Get(Form("v2gammaEta%s",centStringCocktail.Data()));
  TH1D*  hist_v2gammaCocktailOmega = (TH1D*)InputFileCocktail->Get(Form("v2gammaOmega%s",centStringCocktail.Data()));
  TGraphAsymmErrors*  graph_v2gammaCocktail_sys = (TGraphAsymmErrors*)InputFileCocktail->Get(Form("v2gammaCocktail_sys%s",centStringCocktail.Data()));
  
  TH1D*  exampleBin_pcm = (TH1D*)fileInput->Get("h_vn_dir_pcm_toterr_ptbin_9");
  TH1D*  exampleBin_phos = (TH1D*)fileInput->Get("h_vn_dir_phos_toterr_ptbin_9");
  TH1D*  exampleBin_comb = (TH1D*)fileInput->Get("h_vn_dir_comb_toterr_ptbin_9");
  exampleBin_pcm->Scale(1/exampleBin_pcm->GetEntries());
  exampleBin_phos->Scale(1/exampleBin_pcm->GetEntries());
  exampleBin_comb->Scale(1/exampleBin_pcm->GetEntries());
  

  //========================================================
  //cosmetics
  //========================================================
  
  Double_t globalMarkerSize = 1.5;
  Double_t globalLineWidth = 1.5;
  
  graph_v2inc_comb_tot->SetMarkerStyle(20);
  graph_v2inc_comb_tot->SetMarkerSize(globalMarkerSize);
  graph_v2inc_comb_tot->SetMarkerColor(kGray+3);
  graph_v2inc_comb_tot->SetLineColor(kGray+3);
  graph_v2inc_comb_tot->SetLineWidth(globalLineWidth);
  graph_v2inc_comb_tot->SetFillStyle(0);
  
  graph_v2inc_pcm_tot->SetMarkerStyle(34);
  graph_v2inc_pcm_tot->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2inc_pcm_tot->SetMarkerColor(kRed+2);
  graph_v2inc_pcm_tot->SetLineColor(kRed+2);
  graph_v2inc_pcm_tot->SetLineWidth(globalLineWidth);
  graph_v2inc_pcm_tot->SetFillStyle(0);
  
  graph_v2dir_pcm_Rpcm_tot->SetMarkerStyle(28);
  graph_v2dir_pcm_Rpcm_tot->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2dir_pcm_Rpcm_tot->SetMarkerColor(kRed+2);
  graph_v2dir_pcm_Rpcm_tot->SetLineColor(kRed+2);
  graph_v2dir_pcm_Rpcm_tot->SetLineWidth(globalLineWidth);
  graph_v2dir_pcm_Rpcm_tot->SetFillStyle(0);
  
  graph_v2dir_pcm_Rcomb_tot->SetMarkerStyle(34);
  graph_v2dir_pcm_Rcomb_tot->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2dir_pcm_Rcomb_tot->SetMarkerColor(kRed+2);
  graph_v2dir_pcm_Rcomb_tot->SetLineColor(kRed+2);
  graph_v2dir_pcm_Rcomb_tot->SetLineWidth(globalLineWidth);
  graph_v2dir_pcm_Rcomb_tot->SetFillStyle(0);
  
  graph_v2inc_phos_tot->SetMarkerStyle(33);
  graph_v2inc_phos_tot->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2inc_phos_tot->SetMarkerColor(kBlue+2);
  graph_v2inc_phos_tot->SetLineColor(kBlue+2);
  graph_v2inc_phos_tot->SetLineWidth(globalLineWidth);
  graph_v2inc_phos_tot->SetFillStyle(0);
  
  graph_v2dir_phos_Rphos_tot->SetMarkerStyle(27);
  graph_v2dir_phos_Rphos_tot->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2dir_phos_Rphos_tot->SetMarkerColor(kBlue+2);
  graph_v2dir_phos_Rphos_tot->SetLineColor(kBlue+2);
  graph_v2dir_phos_Rphos_tot->SetLineWidth(globalLineWidth);
  graph_v2dir_phos_Rphos_tot->SetFillStyle(0);
  
  graph_v2dir_phos_Rcomb_tot->SetMarkerStyle(33);
  graph_v2dir_phos_Rcomb_tot->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2dir_phos_Rcomb_tot->SetMarkerColor(kBlue+2);
  graph_v2dir_phos_Rcomb_tot->SetLineColor(kBlue+2);
  graph_v2dir_phos_Rcomb_tot->SetLineWidth(globalLineWidth);
  graph_v2dir_phos_Rcomb_tot->SetFillStyle(0);
  
  graph_v2dec_tot->SetMarkerStyle(1);
  graph_v2dec_tot->SetMarkerSize(0);
  graph_v2dec_tot->SetMarkerColor(kGreen+1);
  graph_v2dec_tot->SetLineColor(kGreen+1);
  graph_v2dec_tot->SetLineWidth(globalLineWidth);
//   graph_v2dec_tot->SetFillStyle(1);
  graph_v2dec_tot->SetFillColor(kGreen+1);
  
  graph_v2dir_comb_stat->SetMarkerStyle(20);
  graph_v2dir_comb_stat->SetMarkerSize(globalMarkerSize);
  graph_v2dir_comb_stat->SetMarkerColor(kGray+3);
  graph_v2dir_comb_stat->SetLineColor(kGray+3);
  graph_v2dir_comb_stat->SetLineWidth(globalLineWidth);
  
  graph_v2dir_comb_tot->SetFillColor(kGray);
  graph_v2dir_comb_tot->SetFillStyle(0);
  
  TGraphAsymmErrors* graph_v2inc_ratio_to_comb_pcm = calcGraphRatio(graph_v2inc_pcm_tot,graph_v2inc_comb_tot);
  graph_v2inc_ratio_to_comb_pcm->SetMarkerStyle(34);
  graph_v2inc_ratio_to_comb_pcm->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2inc_ratio_to_comb_pcm->SetMarkerColor(kRed+2);
  graph_v2inc_ratio_to_comb_pcm->SetLineColor(kRed+2);
  graph_v2inc_ratio_to_comb_pcm->SetLineWidth(globalLineWidth);
  graph_v2inc_ratio_to_comb_pcm->SetFillStyle(0);
  TGraphAsymmErrors* graph_v2inc_ratio_to_comb_phos = calcGraphRatio(graph_v2inc_phos_tot,graph_v2inc_comb_tot);
  graph_v2inc_ratio_to_comb_phos->SetMarkerStyle(33);
  graph_v2inc_ratio_to_comb_phos->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2inc_ratio_to_comb_phos->SetMarkerColor(kBlue+2);
  graph_v2inc_ratio_to_comb_phos->SetLineColor(kBlue+2);
  graph_v2inc_ratio_to_comb_phos->SetLineWidth(globalLineWidth);
  graph_v2inc_ratio_to_comb_phos->SetFillStyle(0);
  
  TGraphAsymmErrors* graph_v2dir_ratio_to_comb_pcm_Rpcm = calcGraphRatio(graph_v2dir_comb_tot,graph_v2dir_pcm_Rpcm_tot);
  graph_v2dir_ratio_to_comb_pcm_Rpcm->SetMarkerStyle(28);
  graph_v2dir_ratio_to_comb_pcm_Rpcm->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2dir_ratio_to_comb_pcm_Rpcm->SetMarkerColor(kRed+2);
  graph_v2dir_ratio_to_comb_pcm_Rpcm->SetLineColor(kRed+2);
  graph_v2dir_ratio_to_comb_pcm_Rpcm->SetLineWidth(globalLineWidth);
  graph_v2dir_ratio_to_comb_pcm_Rpcm->SetFillStyle(0);
  
  TGraphAsymmErrors* graph_v2dir_ratio_to_comb_pcm_Rcomb = calcGraphRatio(graph_v2dir_comb_tot,graph_v2dir_pcm_Rcomb_tot);
  graph_v2dir_ratio_to_comb_pcm_Rcomb->SetMarkerStyle(34);
  graph_v2dir_ratio_to_comb_pcm_Rcomb->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2dir_ratio_to_comb_pcm_Rcomb->SetMarkerColor(kRed+2);
  graph_v2dir_ratio_to_comb_pcm_Rcomb->SetLineColor(kRed+2);
  graph_v2dir_ratio_to_comb_pcm_Rcomb->SetLineWidth(globalLineWidth);
  graph_v2dir_ratio_to_comb_pcm_Rcomb->SetFillStyle(0);
  
  TGraphAsymmErrors* graph_v2dir_ratio_to_comb_phos_Rphos = calcGraphRatio(graph_v2dir_comb_tot,graph_v2dir_phos_Rphos_tot);
  graph_v2dir_ratio_to_comb_phos_Rphos->SetMarkerStyle(27);
  graph_v2dir_ratio_to_comb_phos_Rphos->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2dir_ratio_to_comb_phos_Rphos->SetMarkerColor(kBlue+2);
  graph_v2dir_ratio_to_comb_phos_Rphos->SetLineColor(kBlue+2);
  graph_v2dir_ratio_to_comb_phos_Rphos->SetLineWidth(globalLineWidth);
  graph_v2dir_ratio_to_comb_phos_Rphos->SetFillStyle(0);
  
  TGraphAsymmErrors* graph_v2dir_ratio_to_comb_phos_Rcomb = calcGraphRatio(graph_v2dir_comb_tot,graph_v2dir_phos_Rcomb_tot);
  graph_v2dir_ratio_to_comb_phos_Rcomb->SetMarkerStyle(33);
  graph_v2dir_ratio_to_comb_phos_Rcomb->SetMarkerSize(globalMarkerSize+0.4);
  graph_v2dir_ratio_to_comb_phos_Rcomb->SetMarkerColor(kBlue+2);
  graph_v2dir_ratio_to_comb_phos_Rcomb->SetLineColor(kBlue+2);
  graph_v2dir_ratio_to_comb_phos_Rcomb->SetLineWidth(globalLineWidth);
  graph_v2dir_ratio_to_comb_phos_Rcomb->SetFillStyle(0);
  
  hist_v2gammaCocktail->SetMarkerStyle(20);
  hist_v2gammaCocktail->SetMarkerSize(globalMarkerSize+0.4);
  hist_v2gammaCocktail->SetMarkerColor(kGreen+1);
  hist_v2gammaCocktail->SetLineColor(kGreen+1);
  hist_v2gammaCocktail->SetLineWidth(globalLineWidth);
  hist_v2gammaCocktail->SetFillStyle(0);
  
  graph_v2gammaCocktail_sys->SetMarkerStyle(1);
  graph_v2gammaCocktail_sys->SetMarkerSize(0);
  graph_v2gammaCocktail_sys->SetMarkerColor(kGreen+1);
  graph_v2gammaCocktail_sys->SetLineColor(kGreen+1);
  graph_v2gammaCocktail_sys->SetLineWidth(globalLineWidth);
  graph_v2gammaCocktail_sys->SetFillColor(kGreen+1);
  
  hist_v2gammaCocktailPi0->SetMarkerStyle(31);
  hist_v2gammaCocktailPi0->SetMarkerSize(globalMarkerSize+0.4);
  hist_v2gammaCocktailPi0->SetMarkerColor(kMagenta+2);
  hist_v2gammaCocktailPi0->SetLineColor(kMagenta+2);
  hist_v2gammaCocktailPi0->SetLineWidth(globalLineWidth);
  hist_v2gammaCocktailPi0->SetFillStyle(0);
  
  hist_v2gammaCocktailK0s->SetMarkerStyle(29);
  hist_v2gammaCocktailK0s->SetMarkerSize(globalMarkerSize+0.4);
  hist_v2gammaCocktailK0s->SetMarkerColor(kBlue+2);
  hist_v2gammaCocktailK0s->SetLineColor(kBlue+2);
  hist_v2gammaCocktailK0s->SetLineWidth(globalLineWidth);
  hist_v2gammaCocktailK0s->SetFillStyle(0);
  
  hist_v2gammaCocktailEta->SetMarkerStyle(33);
  hist_v2gammaCocktailEta->SetMarkerSize(globalMarkerSize+0.4);
  hist_v2gammaCocktailEta->SetMarkerColor(kRed+2);
  hist_v2gammaCocktailEta->SetLineColor(kRed+2);
  hist_v2gammaCocktailEta->SetLineWidth(globalLineWidth);
  hist_v2gammaCocktailEta->SetFillStyle(0);
  
  hist_v2gammaCocktailOmega->SetMarkerStyle(34);
  hist_v2gammaCocktailOmega->SetMarkerSize(globalMarkerSize+0.4);
  hist_v2gammaCocktailOmega->SetMarkerColor(kYellow+2);
  hist_v2gammaCocktailOmega->SetLineColor(kYellow+2);
  hist_v2gammaCocktailOmega->SetLineWidth(globalLineWidth);
  hist_v2gammaCocktailOmega->SetFillStyle(0);
  
  exampleBin_pcm->SetMarkerStyle(34);
  exampleBin_pcm->SetMarkerSize(globalMarkerSize+0.4);
  exampleBin_pcm->SetMarkerColor(kRed+2);
  exampleBin_pcm->SetLineColor(kRed+2);
  exampleBin_pcm->SetLineWidth(globalLineWidth);
  exampleBin_pcm->SetFillStyle(0);
  
  exampleBin_phos->SetMarkerStyle(33);
  exampleBin_phos->SetMarkerSize(globalMarkerSize+0.4);
  exampleBin_phos->SetMarkerColor(kBlue+2);
  exampleBin_phos->SetLineColor(kBlue+2);
  exampleBin_phos->SetLineWidth(globalLineWidth);
  exampleBin_phos->SetFillStyle(0);
  
  exampleBin_comb->SetMarkerStyle(20);
  exampleBin_comb->SetMarkerSize(globalMarkerSize);
  exampleBin_comb->SetMarkerColor(kGray+3);
  exampleBin_comb->SetLineColor(kGray+3);
  exampleBin_comb->SetLineWidth(globalLineWidth);
  exampleBin_comb->SetFillStyle(0);
  
  EnlargeBinErrorsInX(graph_v2inc_comb_tot);
  EnlargeBinErrorsInX(graph_v2inc_pcm_tot);
  EnlargeBinErrorsInX(graph_v2inc_phos_tot);
  EnlargeBinErrorsInX(graph_v2dec_tot);
  EnlargeBinErrorsInX(graph_v2dir_comb_tot);
  EnlargeBinErrorsInX(graph_v2inc_ratio_to_comb_pcm);
  EnlargeBinErrorsInX(graph_v2inc_ratio_to_comb_phos);
  EnlargeBinErrorsInX(graph_v2dir_pcm_Rpcm_tot);
  EnlargeBinErrorsInX(graph_v2dir_pcm_Rcomb_tot);
  EnlargeBinErrorsInX(graph_v2dir_phos_Rphos_tot);
  EnlargeBinErrorsInX(graph_v2dir_phos_Rcomb_tot);
  EnlargeBinErrorsInX(graph_v2dir_ratio_to_comb_pcm_Rpcm);
  EnlargeBinErrorsInX(graph_v2dir_ratio_to_comb_pcm_Rcomb);
  EnlargeBinErrorsInX(graph_v2dir_ratio_to_comb_phos_Rphos);
  EnlargeBinErrorsInX(graph_v2dir_ratio_to_comb_phos_Rcomb);
  
  TH1F*  graph_v2dir_comb_statEmpty = (TH1F*)Getv2HistStyle();
  graph_v2dir_comb_statEmpty->GetYaxis()->SetRangeUser(-0.09,0.349);
  graph_v2dir_comb_statEmpty->GetYaxis()->SetTitle("v_{2}^{#gamma,dir}");
  
  TH1F*  graph_v2Inc_comb_statEmpty = (TH1F*)Getv2HistStyle();
  graph_v2Inc_comb_statEmpty->GetYaxis()->SetRangeUser(0.01,0.39);
  graph_v2Inc_comb_statEmpty->GetYaxis()->SetTitle("v_{2}^{#gamma,incl}");
  
  TH1F*  graph_v2Inc_comb_Ratio_statEmpty = (TH1F*)GetRatioHistStyle();
  graph_v2Inc_comb_Ratio_statEmpty->GetYaxis()->SetRangeUser(0.29,1.69);
  graph_v2Inc_comb_Ratio_statEmpty->GetYaxis()->SetTitle("v_{2}^{#gamma,incl} / v_{2}^{#gamma,comb}");
  
  TH1F*  graph_v2dir_comb_Ratio_statEmpty = (TH1F*)GetRatioHistStyle();
  graph_v2dir_comb_Ratio_statEmpty->GetYaxis()->SetRangeUser(-1.99,2.99);
  graph_v2dir_comb_Ratio_statEmpty->GetYaxis()->SetTitle("v_{2}^{#gamma,dir} / v_{2}^{#gamma,comb}");
  
  TH1F*  graph_v2Inc_cocktail_empty = (TH1F*)Getv2HistStyle();
  graph_v2Inc_cocktail_empty->GetYaxis()->SetRangeUser(0.01,0.349);
  graph_v2Inc_cocktail_empty->GetYaxis()->SetTitle("v_{2}^{#gamma,decay}");
  
  TH1F*  examplebin_empty = (TH1F*)GetExampleBinStyle();
//   graph_v2Inc_cocktail_empty->GetYaxis()->SetRangeUser(0.01,0.349);
  
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
  
  Double_t moveLegendDown = 0; Double_t moveLegendDown2 = 0;
  if(CentralityLow.CompareTo("20")==0){ moveLegendDown = 0.6; moveLegendDown2 = 0.7; }
  
  //========================================================
  //Plotting and Saving - INCLUSIVE ALL 
  //========================================================
  TCanvas* cIncComb = new TCanvas("cIncComb","",800,800);
  
  TLegend* legIncComb = new TLegend(0.16,0.75,0.4,0.96);
  legIncComb->SetTextSize(0.05);
  legIncComb->AddEntry(graph_v2dec_tot,"Cocktail","f");
  legIncComb->AddEntry(graph_v2inc_comb_tot,"Combined","p");
  legIncComb->AddEntry(graph_v2inc_pcm_tot,"PCM","p");
  legIncComb->AddEntry(graph_v2inc_phos_tot,"PHOS","p");
  
  legIncComb->SetBorderSize(0);

  SetProperMargins();
  graph_v2Inc_comb_statEmpty->Draw();
  graph_v2dec_tot->Draw("SAME e3 P");
  graph_v2inc_phos_tot->Draw("SAME e2 P");
  graph_v2inc_pcm_tot->Draw("SAME e2 P");
  graph_v2inc_comb_tot->Draw("SAME e2 P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),4);
  
  legIncComb->Draw();

  cIncComb->SaveAs(Form("PaperPlots/%s%s_v2inc_combined_all.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //========================================================
  //Plotting and Saving - INCLUSIVE ONLY COMBINED AND DECAY
  //========================================================
  TCanvas* cIncComb2 = new TCanvas("cIncComb","",800,800);
  
  TLegend* legIncComb2 = new TLegend(0.16,0.85-moveLegendDown2,0.4,0.96-moveLegendDown2);
  legIncComb2->SetTextSize(0.05);
  legIncComb2->AddEntry(graph_v2dec_tot,"Cocktail","f");
  legIncComb2->AddEntry(graph_v2inc_comb_tot,"Combined","p");
  
  legIncComb2->SetBorderSize(0);

  SetProperMargins();
  graph_v2Inc_comb_statEmpty->Draw();
  graph_v2dec_tot->Draw("SAME e3 P");
  graph_v2inc_comb_tot->Draw("SAME e2 P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),4);
  
  legIncComb2->Draw();

  cIncComb2->SaveAs(Form("PaperPlots/%s%s_v2inc_combined_OnlyCombDec.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //========================================================
  //Plotting and Saving - INCLUSIVE RATIO TO COMBINED
  //========================================================
  TCanvas* cIncRatiotoComb = new TCanvas("cIncRatiotoComb","",800,800);
  
  TLegend* legIncRatiotoComb = new TLegend(0.16,0.15,0.4,0.26);
  legIncRatiotoComb->SetTextSize(0.05);
  legIncRatiotoComb->AddEntry(graph_v2inc_ratio_to_comb_pcm,"PCM","p");
  legIncRatiotoComb->AddEntry(graph_v2inc_ratio_to_comb_phos,"PHOS","p");
  
  TLine *line1 = new TLine(0,1,7,1);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  
  TLine *line2 = new TLine(0,1.05,7,1.05);
  line2->SetLineStyle(3);
  line2->SetLineWidth(2);
  
  TLine *line3 = new TLine(0,0.95,7,0.95);
  line3->SetLineStyle(3);
  line3->SetLineWidth(2);
  
  legIncRatiotoComb->SetBorderSize(0);

  SetProperMargins();
  graph_v2Inc_comb_Ratio_statEmpty->Draw();
  graph_v2inc_ratio_to_comb_pcm->Draw("SAME e2 P");
  graph_v2inc_ratio_to_comb_phos->Draw("SAME e2 P");
//   graph_v2dec_tot->Draw("SAME P");
//   graph_v2inc_comb_tot->Draw("SAME P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),4);
  line1->Draw();
//   line2->Draw();
//   line3->Draw();
  legIncRatiotoComb->Draw();

  cIncRatiotoComb->SaveAs(Form("PaperPlots/%s%s_v2inc_Ratio_to_Combined.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //========================================================
  //Plotting and Saving - DIRECT
  //========================================================
  
  TCanvas* cDirComb = new TCanvas("cDirComb","",800,800);
  
  TLegend* legDirComb = new TLegend(0.16,0.85,0.4,0.96);
  legDirComb->SetTextSize(0.05);
  legDirComb->AddEntry(graph_v2dir_comb_stat,"statistical","p");
  legDirComb->AddEntry(graph_v2dir_comb_tot,"total","f");
  legDirComb->SetBorderSize(0);
  
  TLine *line4 = new TLine(0,0,7,0);
  line4->SetLineStyle(2);
  line4->SetLineWidth(2);

  SetProperMargins();
  graph_v2dir_comb_statEmpty->Draw();
  graph_v2dir_comb_stat->Draw("SAME P");
  graph_v2dir_comb_tot->Draw("SAME e2 P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  line4->Draw();
  legDirComb->Draw();

  cDirComb->SaveAs(Form("PaperPlots/%s%s_v2dir_combined.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //========================================================
  //Plotting and Saving - DIRECT ALL
  //========================================================
  
  TCanvas* cDirCombAll = new TCanvas("cDirCombAll","",800,800);
  
  TLegend* legDirCombAll = new TLegend(0.16,0.75,0.4,0.96);
  legDirCombAll->SetTextSize(0.03);
  legDirCombAll->AddEntry(graph_v2dir_comb_stat,"comb","p");
  legDirCombAll->AddEntry(graph_v2dir_pcm_Rpcm_tot,"PCM(R_{#gamma,PCM})","p");
  legDirCombAll->AddEntry(graph_v2dir_pcm_Rcomb_tot,"PCM(R_{#gamma,comb})","p");
  legDirCombAll->AddEntry(graph_v2dir_phos_Rphos_tot,"PHOS(R_{#gamma,PHOS})","p");
  legDirCombAll->AddEntry(graph_v2dir_phos_Rcomb_tot,"PHOS(R_{#gamma,comb})","p");
  legDirCombAll->SetBorderSize(0);

  SetProperMargins();
  graph_v2dir_comb_statEmpty->Draw();
  graph_v2dir_pcm_Rpcm_tot->Draw("SAME e2 P");
  graph_v2dir_pcm_Rcomb_tot->Draw("SAME e2 P");
  graph_v2dir_phos_Rphos_tot->Draw("SAME e2 P");
  graph_v2dir_phos_Rcomb_tot->Draw("SAME e2 P");
  graph_v2dir_comb_stat->Draw("SAME P");
  graph_v2dir_comb_tot->Draw("SAME e2 P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  line4->Draw();
  legDirCombAll->Draw();

  cDirCombAll->SaveAs(Form("PaperPlots/%s%s_v2dir_All.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //========================================================
  //Plotting and Saving - DIRECT ALL RATIO TO COMBINED
  //========================================================
  
  TCanvas* cDirCombAllRatio = new TCanvas("cDirCombAllRatio","",800,800);
  
  TLegend* legDirCombAllRatio = new TLegend(0.16,0.15,0.4,0.36);
  legDirCombAllRatio->SetTextSize(0.03);
  legDirCombAllRatio->AddEntry(graph_v2dir_ratio_to_comb_pcm_Rpcm,"PCM(R_{#gamma,PCM})","p");
  legDirCombAllRatio->AddEntry(graph_v2dir_ratio_to_comb_pcm_Rcomb,"PCM(R_{#gamma,comb})","p");
  legDirCombAllRatio->AddEntry(graph_v2dir_ratio_to_comb_phos_Rphos,"PHOS(R_{#gamma,PHOS})","p");
  legDirCombAllRatio->AddEntry(graph_v2dir_ratio_to_comb_phos_Rcomb,"PHOS(R_{#gamma,comb})","p");
  legDirCombAllRatio->SetBorderSize(0);

  SetProperMargins();
  graph_v2dir_comb_Ratio_statEmpty->Draw();
  
  graph_v2dir_ratio_to_comb_pcm_Rpcm->Draw("SAME e2 P");
  graph_v2dir_ratio_to_comb_pcm_Rcomb->Draw("SAME e2 P");
  graph_v2dir_ratio_to_comb_phos_Rphos->Draw("SAME e2 P");
  graph_v2dir_ratio_to_comb_phos_Rcomb->Draw("SAME e2 P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  line1->Draw();
  legDirCombAllRatio->Draw();

  cDirCombAllRatio->SaveAs(Form("PaperPlots/%s%s_v2dir_All_Ratio.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //========================================================
  //Opening files for theory comparison
  //========================================================
  if(IncludeTheory && (CentralityLow.CompareTo("0")==0 || CentralityLow.CompareTo("20")==0) ){

      
    //========================================================
    //DIRECT THEORY
    //========================================================
    TString centStringTheory = "";
    if(CentralityLow.CompareTo("0")==0) centStringTheory = "0020";
    if(CentralityLow.CompareTo("20")==0) centStringTheory = "2040";
    ifstream inputfile;
    inputfile.open(Form("TaskFlow/Theory/direct_photons_cent%s_Paquet2016_stripped.dat",centStringTheory.Data()));
    
    if(inputfile.is_open()) cout << "file is open..." << endl;
    
    Float_t tpt[11], tyield[11], tyieldstat[11], tv2[11], tv2stat[11], tv3[11], tv3stat[11], tv2stathigh[11], tv2statlow[11];
    Int_t nlinesfile = 0;
    Int_t tn = 11;
    
    while(1){
      cout << nlinesfile << endl;
      inputfile >> tpt[nlinesfile] >> tyield[nlinesfile] >> tyieldstat[nlinesfile] >> tv2[nlinesfile] >> tv2stat[nlinesfile] >> tv3[nlinesfile] >> tv3stat[nlinesfile];
      cout << tpt[nlinesfile] << " " << tyield[nlinesfile] << " " << tyieldstat[nlinesfile] << " " << tv2[nlinesfile] << " " << tv2stat[nlinesfile] << " " << tv3[nlinesfile] << " " << tv3stat[nlinesfile];
      if(!inputfile.good()) break;
      nlinesfile++;
    }
    
    TGraphAsymmErrors* graphDirectTheory = new TGraphAsymmErrors(tn,tpt,tv2,0,0,tv2stat,tv2stat);
    
    TCanvas* cTheory = new TCanvas("cTheory","",800,800);

    SetProperMargins();
    graph_v2dir_comb_statEmpty->Draw();
    graphDirectTheory->SetFillColor(kCyan+1);
//     graphDirectTheory->SetFillStyle(0);
    graphDirectTheory->SetLineColor(kCyan+1);
    graphDirectTheory->Draw("SAME 3");
    graph_v2dec_tot->Draw("SAME e3 P");
    graph_v2dir_comb_stat->Draw("SAME P");
    graph_v2dir_comb_tot->Draw("SAME e2 P");
    DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    line4->Draw();
    TLegend* legDirCombTheory = new TLegend(0.16,0.80,0.4,0.96);
    legDirCombTheory->SetTextSize(0.05);
    legDirCombTheory->AddEntry(graph_v2dir_comb_stat,"Data","p");
    legDirCombTheory->AddEntry(graphDirectTheory,"Theory","f");
    legDirCombTheory->AddEntry(graph_v2dec_tot,"Cocktail","f");
    legDirCombTheory->SetBorderSize(0);
    legDirCombTheory->Draw();

    cTheory->SaveAs(Form("PaperPlots/%s%s_v2dir_combined_theory.eps",CentralityLow.Data(),CentralityHigh.Data()));
    
    
    TCanvas* cTheoryALL = new TCanvas("cTheoryALL","",800,800);

    SetProperMargins();
    graph_v2dir_comb_statEmpty->Draw();
    graphDirectTheory->Draw("SAME 3");
    graph_v2dec_tot->Draw("SAME e3 P");
    graph_v2dir_pcm_Rpcm_tot->Draw("SAME e2 P");
    graph_v2dir_phos_Rphos_tot->Draw("SAME e2 P");
    DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    line4->Draw();
    TLegend* legDirCombTheoryALL = new TLegend(0.16,0.80,0.4,0.96);
    legDirCombTheoryALL->SetTextSize(0.05);
    legDirCombTheoryALL->AddEntry(graph_v2dir_pcm_Rpcm_tot,"PCM","p");
    legDirCombTheoryALL->AddEntry(graph_v2dir_phos_Rphos_tot,"PHOS","p");
    legDirCombTheoryALL->AddEntry(graphDirectTheory,"Theory","f");
    legDirCombTheoryALL->AddEntry(graph_v2dec_tot,"Cocktail","f");
    legDirCombTheoryALL->SetBorderSize(0);
    legDirCombTheoryALL->Draw();

    cTheoryALL->SaveAs(Form("PaperPlots/%s%s_v2dir_ALL_theory.eps",CentralityLow.Data(),CentralityHigh.Data()));
    
    TCanvas* cTheoryPCM = new TCanvas("cTheoryPCM","",800,800);

    SetProperMargins();
    graph_v2dir_comb_statEmpty->Draw();
    graphDirectTheory->Draw("SAME 3");
    graph_v2dec_tot->Draw("SAME e3 P");
    graph_v2dir_pcm_Rpcm_tot->Draw("SAME e2 P");
    DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    line4->Draw();
    TLegend* legDirCombTheoryPCM = new TLegend(0.16,0.80,0.4,0.96);
    legDirCombTheoryPCM->SetTextSize(0.05);
    legDirCombTheoryPCM->AddEntry(graph_v2dir_pcm_Rpcm_tot,"PCM","p");
    legDirCombTheoryPCM->AddEntry(graphDirectTheory,"Theory","f");
    legDirCombTheoryPCM->AddEntry(graph_v2dec_tot,"Cocktail","f");
    legDirCombTheoryPCM->SetBorderSize(0);
    legDirCombTheoryPCM->Draw();

    cTheoryPCM->SaveAs(Form("PaperPlots/%s%s_v2dir_PCM_theory.eps",CentralityLow.Data(),CentralityHigh.Data()));
    
    TCanvas* cTheoryPHOS = new TCanvas("cTheoryPHOS","",800,800);

    SetProperMargins();
    graph_v2dir_comb_statEmpty->Draw();
    graphDirectTheory->Draw("SAME 3");
    graph_v2dec_tot->Draw("SAME e3 P");
    graph_v2dir_phos_Rphos_tot->Draw("SAME e2 P");
    DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    line4->Draw();
    TLegend* legDirCombTheoryPHOS = new TLegend(0.16,0.80,0.4,0.96);
    legDirCombTheoryPHOS->SetTextSize(0.05);
    legDirCombTheoryPHOS->AddEntry(graph_v2dir_phos_Rphos_tot,"PHOS","p");
    legDirCombTheoryPHOS->AddEntry(graphDirectTheory,"Theory","f");
    legDirCombTheoryPHOS->AddEntry(graph_v2dec_tot,"Cocktail","f");
    legDirCombTheoryPHOS->SetBorderSize(0);
    legDirCombTheoryPHOS->Draw();

    cTheoryPHOS->SaveAs(Form("PaperPlots/%s%s_v2dir_PHOS_theory.eps",CentralityLow.Data(),CentralityHigh.Data()));
 
    
    inputfile.close();
    
//       CalculatepValue(graph_v2dir_comb_tot,graph_v2dir_pcm_Rpcm_tot,1,3);
//       CalculatepValue(graph_v2dir_comb_tot,graph_v2dir_phos_Rphos_tot,1,3);
      
//       CalculatepValue(graph_v2dec_tot,graph_v2dir_pcm_Rpcm_tot,1,3);
//       CalculatepValue(graph_v2dec_tot,graph_v2dir_phos_Rphos_tot,1,3);
//       CalculatepValue(graph_v2dec_tot,graph_v2dir_comb_tot,1,3);
      
//       CalculatepValue(graph_v2dir_phos_Rphos_tot,graph_v2dir_pcm_Rpcm_tot,1,3);
//       CalculatepValue(graph_v2dir_pcm_Rpcm_tot,graph_v2dir_pcm_Rpcm_tot,1,3);
    
    
    //========================================================
    //INCLUSIVE THEORY
    //========================================================
    ifstream inputfileInc;
    inputfileInc.open(Form("TaskFlow/Theory/inclusive_photons_cent%s_Paquet2016_stripped.dat",centStringTheory.Data()));
    
    if(inputfileInc.is_open()) cout << "file is open..." << endl;
    nlinesfile = 0;
    Float_t tpt_2[15], tyield_2[15], tyieldstat_2[15], tv2_2[15], tv2stat_2[15], tv3_2[15], tv3stat_2[15];
    tn = 15;
    
    while(1){
      cout << nlinesfile << endl;
      inputfileInc >> tpt_2[nlinesfile] >> tyield_2[nlinesfile] >> tyieldstat_2[nlinesfile] >> tv2_2[nlinesfile] >> tv2stat_2[nlinesfile] >> tv3_2[nlinesfile] >> tv3stat_2[nlinesfile];
      cout << tpt_2[nlinesfile] << " " << tyield_2[nlinesfile] << " " << tyieldstat_2[nlinesfile] << " " << tv2_2[nlinesfile] << " " << tv2stat_2[nlinesfile] << " " << tv3_2[nlinesfile] << " " << tv3stat_2[nlinesfile];
      if(!inputfileInc.good()) break;
      nlinesfile++;
    }
    
    TGraphAsymmErrors* graphIncTheory = new TGraphAsymmErrors(tn,tpt_2,tv2_2,0,0,tv2stat_2,tv2stat_2);
    
    TCanvas* cTheoryInc = new TCanvas("cTheoryInc","",800,800);

    SetProperMargins();
    graph_v2Inc_comb_statEmpty->Draw();
    graphIncTheory->SetFillColor(kCyan+1);
//     graphIncTheory->SetFillColorAlpha(kCyan+1,0.95);
    graphIncTheory->SetLineColor(kCyan+1);
    graph_v2dec_tot->Draw("SAME e3 P");
    graphIncTheory->Draw("SAME 3");
    graph_v2inc_comb_tot->Draw("SAME e2 P");
    DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),4);

    TLegend* legIncCombTheory = new TLegend(0.16,0.80,0.4,0.96);
    legIncCombTheory->SetTextSize(0.05);
    legIncCombTheory->AddEntry(graph_v2inc_comb_tot,"Combined","p");
    legIncCombTheory->AddEntry(graph_v2dec_tot,"Cocktail","f");
    legIncCombTheory->AddEntry(graphIncTheory,"Theory","f");
    legIncCombTheory->SetBorderSize(0);
    legIncCombTheory->Draw();

    cTheoryInc->SaveAs(Form("PaperPlots/%s%s_v2inc_combined_theory.eps",CentralityLow.Data(),CentralityHigh.Data()));
    
    inputfile.close();
  
  }
  
  //========================================================
  //Plotting and Saving - Cocktail
  //========================================================
  
  TCanvas* cIncCocktail = new TCanvas("cIncCocktail","",800,800);
  
  TLegend* legIncCocktail = new TLegend(0.16,0.70,0.4,0.96);
  legIncCocktail->SetTextSize(0.05);
  legIncCocktail->AddEntry(graph_v2gammaCocktail_sys,"Cocktail","f");
  legIncCocktail->AddEntry(hist_v2gammaCocktailPi0,"#pi^{0}#rightarrow2#gamma","p");
  legIncCocktail->AddEntry(hist_v2gammaCocktailEta,"#eta#rightarrow2#gamma","p");
  legIncCocktail->AddEntry(hist_v2gammaCocktailK0s,"K_{s}^{0}#rightarrow2#pi^{0}","p");
  legIncCocktail->AddEntry(hist_v2gammaCocktailOmega,"#omega#rightarrow#gamma#pi^{0}","p");
  
  legIncCocktail->SetBorderSize(0);

  SetProperMargins();
  graph_v2Inc_cocktail_empty->Draw();
  graph_v2gammaCocktail_sys->Draw("SAME e3 P");
  hist_v2gammaCocktailOmega->Draw("SAME e2 P");
  hist_v2gammaCocktailK0s->Draw("SAME e2 P");
  hist_v2gammaCocktailEta->Draw("SAME e2 P");
  hist_v2gammaCocktailPi0->Draw("SAME e2 P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),6);
  
  legIncCocktail->Draw();

  cIncCocktail->SaveAs(Form("PaperPlots/%s%s_v2inc_cocktail.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //========================================================
  //Plotting and Saving - Example bin
  //========================================================
  
  TCanvas* cExampleBin = new TCanvas("cExampleBin","",800,800);
  
  TLegend* legExampleBin = new TLegend(0.16,0.80,0.4,0.96);
  legExampleBin->SetTextSize(0.05);
  legExampleBin->AddEntry(exampleBin_pcm,"PCM","p");
  legExampleBin->AddEntry(exampleBin_phos,"PHOS","p");
  legExampleBin->AddEntry(exampleBin_comb,"Combined","p");
  
  legExampleBin->SetBorderSize(0);

  examplebin_empty->GetYaxis()->SetRangeUser(0,exampleBin_comb->GetMaximum()*1.2);
  SetProperMargins();
  examplebin_empty->Draw();
  exampleBin_pcm->Draw("SAME P");
  exampleBin_phos->Draw("SAME P");
  exampleBin_comb->Draw("SAME P");
  
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  legExampleBin->Draw();

  cExampleBin->SaveAs(Form("PaperPlots/%s%s_examplebin.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  
}

void SetProperMargins(){
  gPad->SetLeftMargin(0.14);
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
  histoEmpty->GetYaxis()->SetTitle("v_{2}^{#gamma,incl}");
  histoEmpty->GetXaxis()->SetTitleSize(0.06);
  histoEmpty->GetYaxis()->SetTitleSize(0.06);
  histoEmpty->GetXaxis()->SetTitleOffset(0.8);
  histoEmpty->GetYaxis()->SetTitleOffset(1.05);
  histoEmpty->GetXaxis()->SetRangeUser(0.0,6.99);
  histoEmpty->GetYaxis()->SetRangeUser(0,0.249);
  
  return histoEmpty;
}

TH1F* GetRatioHistStyle(){
 
  TH1F*  histoEmpty = new TH1F("","",100,0,10);
  
  for(Int_t i=0;i<100;i++){
    histoEmpty->SetBinContent(i,-100.);
  }
  
  
  histoEmpty->GetXaxis()->SetTitle("#it{p}_{T}(GeV/c)");
  histoEmpty->GetYaxis()->SetTitle("v_{2}^{#gamma,meas.} / v_{2}^{#gamma,comb}");
  histoEmpty->GetXaxis()->SetTitleSize(0.06);
  histoEmpty->GetYaxis()->SetTitleSize(0.06);
  histoEmpty->GetXaxis()->SetTitleOffset(0.8);
  histoEmpty->GetYaxis()->SetTitleOffset(1.05);
  histoEmpty->GetXaxis()->SetRangeUser(0.0,6.99);
  histoEmpty->GetYaxis()->SetRangeUser(0,0.249);
  
  return histoEmpty;
}

TH1F* GetExampleBinStyle(){
 
  TH1F*  histoEmpty = new TH1F("","",100,-0.5,0.5);
  
  for(Int_t i=0;i<100;i++){
    histoEmpty->SetBinContent(i,-100.);
  }
  
  
  histoEmpty->GetXaxis()->SetTitle("v_{2}^{#gamma,direct}");
  histoEmpty->GetYaxis()->SetTitle("norm. counts");
  histoEmpty->GetXaxis()->SetTitleSize(0.06);
  histoEmpty->GetYaxis()->SetTitleSize(0.06);
  histoEmpty->GetXaxis()->SetTitleOffset(0.7);
  histoEmpty->GetYaxis()->SetTitleOffset(1.05);
  histoEmpty->GetXaxis()->SetRangeUser(-0.19,0.49);
  histoEmpty->GetYaxis()->SetRangeUser(0,0.179);
  
  return histoEmpty;
}

void DrawDirectInfoLabelOnPlot(TString cent1, TString cent2, Int_t place){
  
  TLatex T1;
  T1.SetTextSize(0.04);
  T1.SetTextAlign(12);
  T1.SetNDC();
  
  if(place==1){
    T1.DrawLatex(0.15, 0.94, Form("#bf{%s-%s %% PbPb, #sqrt{s_{_{NN}}}=2.76TeV}",cent1.Data(),cent2.Data()));
//     T1.DrawLatex(0.15, 0.88, "#bf{#gamma_{direct}}");
  }
  
  if(place==2){
    T1.DrawLatex(0.48, 0.94, Form("#bf{%s-%s %% PbPb, #sqrt{s_{_{NN}}}=2.76TeV}",cent1.Data(),cent2.Data()));
//     T1.DrawLatex(0.87, 0.88, "#bf{#gamma_{direct}}");
  }
  
  if(place==3){
    T1.DrawLatex(0.15, 0.94, Form("#bf{%s-%s %% PbPb, #sqrt{s_{_{NN}}}=2.76TeV}",cent1.Data(),cent2.Data()));
//     T1.DrawLatex(0.15, 0.88, "#bf{#gamma_{inclusive}}");
  }
  
  if(place==4){
    T1.DrawLatex(0.48, 0.94, Form("#bf{%s-%s %% PbPb, #sqrt{s_{_{NN}}}=2.76TeV}",cent1.Data(),cent2.Data()));
//     T1.DrawLatex(0.83, 0.88, "#bf{#gamma_{inclusive}}");
  }
  
  if(place==5){
    T1.DrawLatex(0.15, 0.94, Form("#bf{%s-%s %% PbPb, #sqrt{s_{_{NN}}}=2.76TeV}",cent1.Data(),cent2.Data()));
//     T1.DrawLatex(0.15, 0.88, "#bf{#gamma_{decay}}");
  }
  
  if(place==6){
    T1.DrawLatex(0.48, 0.94, Form("#bf{%s-%s %% PbPb, #sqrt{s_{_{NN}}}=2.76TeV}",cent1.Data(),cent2.Data()));
//     T1.DrawLatex(0.83, 0.88, "#bf{#gamma_{decay}}");
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

TGraphAsymmErrors* calcGraphRatio(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2){
    
  const Int_t NPoints = g1->GetN();
  Double_t x[NPoints], y[NPoints], exl[NPoints], exh[NPoints], eyl[NPoints], eyh[NPoints];
  Double_t* g1_x = g1->GetX();
  Double_t* g1_y = g1->GetY();
  Double_t* g1_EYhigh = g1->GetEYhigh();

  Double_t* g2_x = g2->GetX();
  Double_t* g2_y = g2->GetY();
  Double_t* g2_EYhigh = g2->GetEYhigh();
  
  for(Int_t ibin=0;ibin<NPoints;ibin++){
    x[ibin] = g1_x[ibin];
    y[ibin] = g1_y[ibin]/g2_y[ibin];
    exl[ibin] = 0;
    exh[ibin] = 0;
//     eyl[ibin] = (g1_y[ibin]/g2_y[ibin]) * TMath::Sqrt(pow((g1_EYhigh[ibin]/g1_y[ibin]),2)+pow((g2_EYhigh[ibin]/g2_y[ibin]),2));
    eyl[ibin] = (g1_y[ibin]/g2_y[ibin]) * TMath::Sqrt(pow((g1_EYhigh[ibin]/g1_y[ibin]),2));
    eyh[ibin] = eyl[ibin];
  }
  
  TGraphAsymmErrors* outputGraph = new TGraphAsymmErrors(NPoints,x,y,exl,exh,eyl,eyh);
  return outputGraph;
    
}

void CalculatepValue(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, Double_t startpT, Double_t endpT){
    
    Double_t* g1_x = g1->GetX();
    Double_t* g1_y = g1->GetY();
    Double_t* g1_EYhigh = g1->GetEYhigh();
    
    Double_t* g2_x = g2->GetX();
    Double_t* g2_y = g2->GetY();
    Double_t* g2_EYhigh = g2->GetEYhigh();
    
    Double_t BinEdges[19] = {0, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.3, 3.7, 4.1, 4.6, 5.4, 6.2, 7.0};
    TH1D* h1 = new TH1D("","",18,BinEdges);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(1.6);
    TH1D* h2 = new TH1D("","",18,BinEdges);
    h2->SetMarkerStyle(21);
    h2->SetMarkerSize(1.6);

    for(Int_t i=2 ; i<13; i++){
        h1->Fill(g1_x[i-2],g1_y[i-2]);
        h1->SetBinError(i,g1_EYhigh[i-2]);
        h2->Fill(g2_x[i-2],g2_y[i-2]);
        h2->SetBinError(i,g2_EYhigh[i-2]);
    }
    cout << "CHI2 TESTING ! " << endl;
    h1->Chi2Test(h2,"WWP");
    cout << "DONE WITH CHI2 TESTING ! " << endl;
    
    TCanvas* cTemp = new TCanvas("cTemp","",800,800);
    h1->Draw();
    h2->Draw("SAME");
    g1->Draw("SAME P");
    g2->Draw("SAME P");
    
    
        
    Double_t g1_y_average = 0;
    Double_t g1_y_error_average = 0;
    Double_t g1_y_error_average_corr = 0;
    for(Int_t i=1 ; i<11; i++){
        
        g1_y_average+=(g1_y[i]*(BinEdges[i+2]-BinEdges[i+1]));
        g1_y_error_average+=(g1_EYhigh[i]*g1_EYhigh[i]*(BinEdges[i+2]-BinEdges[i+1]));
        g1_y_error_average_corr+=(((g1_EYhigh[i]-0.01)*(g1_EYhigh[i]-0.01))*(BinEdges[i+2]-BinEdges[i+1]));
        
    }
    
    g1_y_average/=(10*(3.3-1.1));
    cout << "g1_y_average = " << g1_y_average << endl;
    
    g1_y_error_average= TMath::Sqrt(g1_y_error_average);
    g1_y_error_average/=(TMath::Sqrt(10)*(3.3-1.1));
    cout << "g1_y_error_average = " << g1_y_error_average << endl;
    
    cout << "n sigma above 0 = " << g1_y_average/g1_y_error_average << endl;
    
    g1_y_error_average_corr= TMath::Sqrt(g1_y_error_average_corr);
    g1_y_error_average_corr/=(TMath::Sqrt(10)*(3.3-1.1));
    g1_y_error_average_corr = TMath::Sqrt(g1_y_error_average_corr*g1_y_error_average_corr + 0.01*0.01);
    cout << "g1_y_error_average_corr = " << g1_y_error_average_corr << endl;
    
    cout << "n sigma above 0 correlated = " << g1_y_average/g1_y_error_average_corr << endl;
    
}