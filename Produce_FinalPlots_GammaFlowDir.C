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
TGraphAsymmErrors* calcGraphRatio(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2);

void  Produce_FinalPlots_GammaFlowDir(
                                      TString CentralityLow = "0",
                                      TString CentralityHigh = "20",
                                      Bool_t IncludeTheory = kTRUE
){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================

  TString centStringInputfile = "";
  if(CentralityLow.CompareTo("0")==0) centStringInputfile  = "00-20";
  if(CentralityLow.CompareTo("20")==0) centStringInputfile = "20-40";
  TString InputFileName1  = Form("TaskFlow/V2dir_calculation_PCM_PHOS_combined/output/v2dir_pcm_phos_comb_%s.root",centStringInputfile.Data());
  
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
  
  graph_v2dec_tot->SetMarkerStyle(21);
  graph_v2dec_tot->SetMarkerSize(globalMarkerSize);
  graph_v2dec_tot->SetMarkerColor(kGreen+2);
  graph_v2dec_tot->SetLineColor(kGreen+2);
  graph_v2dec_tot->SetLineWidth(globalLineWidth);
  graph_v2dec_tot->SetFillStyle(0);
  
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
  graph_v2dir_comb_statEmpty->GetYaxis()->SetRangeUser(-0.29,0.49);
  
  TH1F*  graph_v2Inc_comb_statEmpty = (TH1F*)Getv2HistStyle();
  graph_v2Inc_comb_statEmpty->GetYaxis()->SetRangeUser(0.01,0.249);
  
  TH1F*  graph_v2Inc_comb_Ratio_statEmpty = (TH1F*)GetRatioHistStyle();
  graph_v2Inc_comb_Ratio_statEmpty->GetYaxis()->SetRangeUser(0.79,1.21);
  
  TH1F*  graph_v2dir_comb_Ratio_statEmpty = (TH1F*)GetRatioHistStyle();
  graph_v2dir_comb_Ratio_statEmpty->GetYaxis()->SetRangeUser(0.01,1.99);
  
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
  
  TLegend* legIncComb = new TLegend(0.16,0.75-moveLegendDown,0.4,0.96-moveLegendDown);
  legIncComb->SetTextSize(0.05);
  legIncComb->AddEntry(graph_v2inc_comb_tot,"Combined","p");
  legIncComb->AddEntry(graph_v2dec_tot,"Cocktail","p");
  legIncComb->AddEntry(graph_v2inc_pcm_tot,"PCM","p");
  legIncComb->AddEntry(graph_v2inc_phos_tot,"PHOS","p");
  
  legIncComb->SetBorderSize(0);

  SetProperMargins();
  graph_v2Inc_comb_statEmpty->Draw();
  graph_v2inc_phos_tot->Draw("SAME e2 P");
  graph_v2inc_pcm_tot->Draw("SAME e2 P");
  graph_v2dec_tot->Draw("SAME e2 P");
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
  legIncComb2->AddEntry(graph_v2inc_comb_tot,"Combined","p");
  legIncComb2->AddEntry(graph_v2dec_tot,"Cocktail","p");
  
  legIncComb2->SetBorderSize(0);

  SetProperMargins();
  graph_v2Inc_comb_statEmpty->Draw();
  graph_v2dec_tot->Draw("SAME e2 P");
  graph_v2inc_comb_tot->Draw("SAME e2 P");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),4);
  
  legIncComb2->Draw();

  cIncComb2->SaveAs(Form("PaperPlots/%s%s_v2inc_combined_OnlyCombDec.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
  //========================================================
  //Plotting and Saving - INCLUSIVE RATIO TO COMBINED
  //========================================================
  TCanvas* cIncRatiotoComb = new TCanvas("cIncRatiotoComb","",800,800);
  
  TLegend* legIncRatiotoComb = new TLegend(0.16,0.15,0.6,0.26);
  legIncRatiotoComb->SetTextSize(0.05);
  legIncRatiotoComb->AddEntry(graph_v2inc_ratio_to_comb_pcm,"PCM/comb","p");
  legIncRatiotoComb->AddEntry(graph_v2inc_ratio_to_comb_phos,"PHOS/comb","p");
  
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
  line2->Draw();
  line3->Draw();
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
  
  TLegend* legDirCombAll = new TLegend(0.16,0.15,0.4,0.36);
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
  
  TLegend* legDirCombAllRatio = new TLegend(0.16,0.15,0.6,0.36);
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
    graphDirectTheory->SetFillColor(kCyan+2);
    graphDirectTheory->SetLineColor(kCyan+2);
//     graphDirectTheory->SetFillStyle(0);
    graphDirectTheory->Draw("SAME 3");
    graph_v2dir_comb_stat->Draw("SAME P");
    graph_v2dir_comb_tot->Draw("SAME e2 P");
    DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    
    line4->Draw();
    TLegend* legDirCombTheory = new TLegend(0.16,0.85,0.4,0.96);
    legDirCombTheory->SetTextSize(0.05);
    legDirCombTheory->AddEntry(graph_v2dir_comb_stat,"Data","p");
    legDirCombTheory->AddEntry(graphDirectTheory,"Theory","f");
    legDirCombTheory->SetBorderSize(0);
    legDirCombTheory->Draw();

    cTheory->SaveAs(Form("PaperPlots/%s%s_v2dir_combined_theory.eps",CentralityLow.Data(),CentralityHigh.Data()));
    
    inputfile.close();
  
  }
  
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
 
  TH1F*  histoEmpty = new TH1F("","",100,0,10);
  
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
    eyl[ibin] = (g1_y[ibin]/g2_y[ibin]) * TMath::Sqrt(pow((g1_EYhigh[ibin]/g1_y[ibin]),2)+pow((g2_EYhigh[ibin]/g2_y[ibin]),2));
    eyh[ibin] = eyl[ibin];
  }
  
  TGraphAsymmErrors* outputGraph = new TGraphAsymmErrors(NPoints,x,y,exl,exh,eyl,eyh);
  return outputGraph;
    
}