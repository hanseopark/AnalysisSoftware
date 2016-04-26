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

void MaskPoints(TH1F* histo, const int start, const int end);

void  DirectPhotonFlowToy(Double_t purity = 0.96, Double_t v2CorrStr = 1.5, TString select = "pion"){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  //open datafile with v2 gamma inclusive
  TFile* fileData = new TFile("/home/mike/0_directphoton/0_analysis/160125_toy_v2_effect/toy_2040.root");
  
  //Histogram settings
  TH1F*  histo2040 = (TH1F*)fileData->Get("2040");
  histo2040->SetLineColor(kRed+1);
  histo2040->SetMarkerColor(kRed+1);
  histo2040->SetMarkerStyle(25);
  histo2040->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  histo2040->GetYaxis()->SetTitle("v_{2}  ");
  histo2040->GetXaxis()->SetTitleSize(0.05);
  histo2040->GetXaxis()->SetTitleOffset(0.9);
  histo2040->GetYaxis()->SetTitleOffset(1.2);
  histo2040->GetXaxis()->SetRangeUser(0.0,8.0);
  histo2040->GetYaxis()->SetRangeUser(0,0.25);
  
  TH1F* histo2040_corrected = (TH1F*)histo2040->Clone();
  histo2040_corrected->SetLineColor(kBlue+1);
  histo2040_corrected->SetMarkerColor(kBlue+1);
  histo2040_corrected->SetMarkerStyle(27);
  
  TH1F* histo2040_v2Dir_uncorrected = (TH1F*)histo2040_corrected->Clone();
  histo2040_v2Dir_uncorrected->SetLineColor(kRed+1);
  histo2040_v2Dir_uncorrected->SetMarkerColor(kRed+1);
  histo2040_v2Dir_uncorrected->SetMarkerStyle(25);
  histo2040_v2Dir_uncorrected->GetXaxis()->SetRangeUser(0.0,8.0);
  histo2040_v2Dir_uncorrected->GetYaxis()->SetRangeUser(-0.3,0.4);
  
  TH1F* histo2040_v2Dir_corrected = (TH1F*)histo2040_corrected->Clone();
  histo2040_v2Dir_corrected->SetMarkerStyle(27);
  
  TH1F* histo2040_cocktail = (TH1F*)histo2040_corrected->Clone();
  histo2040_cocktail->SetLineColor(kGreen+2);
  histo2040_cocktail->SetMarkerColor(kGreen+2);
  histo2040_cocktail->SetMarkerStyle(24);
  
  TH1F* histo2040_RGamma = (TH1F*)histo2040_corrected->Clone();
  histo2040_RGamma->SetLineColor(kBlack);
  histo2040_RGamma->SetMarkerColor(kBlack);
  histo2040_RGamma->SetMarkerStyle(25);
  histo2040_RGamma->GetYaxis()->SetTitle("R_{#gamma}  ");
  histo2040_RGamma->GetYaxis()->SetRangeUser(0.9,1.3);
  
  TLatex T1;
  T1.SetTextSize(0.05);
  T1.SetTextAlign(12);
  T1.SetNDC();
  
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  //========================================================
  //Calculations
  //========================================================
  //Get the v2 values of total sample
  
  //number of signal and background in initial sample
  Double_t Nsig = purity;
  Double_t Nback = 1.0 - Nsig;
  
  //get bin values for v2
  Double_t v2T_arr[28];
  for(int i=0;i<27;i++){
    v2T_arr[i] = histo2040->GetBinContent(i+2);
  }
  
  //arrays with v2 for backgrounds
  //x values 0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.85,3.15,3.45,3.85,4.25,5,5.8,6.6,7.5,9.5
  Double_t R_arr[24] = {0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.986975,1.03893,1.06285,1.05714,1.05956,1.07441,1.06879,1.07893,1.09772,1.10171,1.10937,1.12906,1.12171,1.1481,1.17943,1.18174};
  Double_t v2C_arr[24] = {0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.15,0.166,0.1816,0.185,0.19,0.194,0.201,0.2035,0.205,0.204,0.1935,0.183,0.172,0.158,0.134,0.118};
  Double_t v2_arr_pion[24] = {0.0228068,0.0313851,0.0475462,0.0636896,0.0788676,0.0926556,0.106052,0.118539,0.135,0.155,0.172,0.185,0.195,0.204,0.209,0.211,0.212,0.209874,0.205004,0.198,0.184,0.168,0.147153,0.127347};
  Double_t v2_arr_proton[24] = {0.0001,0.0001,0.00984868,0.0105599,0.0160637,0.0228682,0.0329241,0.0421492,0.059,0.083,0.111,0.138,0.165,0.189,0.207,0.226,0.2395,0.2586,0.268,0.275,0.273,0.271,0.242208,0.22818};
  Double_t v2_arr_kaon[24] = {0.0001,0.0001,0.0166193,0.0326524,0.0438,0.0576794,0.0724384,0.0872708,0.109,0.132,0.152,0.171,0.184,0.197,0.204,0.211,0.2155,0.218099,0.216,0.215,0.209,0.001,0.001,0.001};
  Double_t v2R_arr[24];
  if(select.CompareTo("pion") == 0){
    for(int i=0;i<24;i++){
      v2R_arr[i] = v2_arr_pion[i];
    }
  }else if(select.CompareTo("proton") == 0){
    for(int i=0;i<24;i++){
      v2R_arr[i] = v2_arr_proton[i];
    }
  }else if(select.CompareTo("pion+proton") == 0 || select.CompareTo("proton+pion") == 0){
    for(int i=0;i<24;i++){
      v2R_arr[i] = v2_arr_pion[i]+v2_arr_proton[i];
    }
  }else if(select.CompareTo("kaon") == 0){
    for(int i=0;i<24;i++){
      v2R_arr[i] = v2_arr_kaon[i];
    }
  }else if(select.CompareTo("pion+proton+kaon") == 0){
    for(int i=0;i<24;i++){
      v2R_arr[i] = v2_arr_pion[i]+v2_arr_proton[i]+v2_arr_kaon[i];
    }
  }else{
    cout << "No correct selection" << endl;
    return;
  }
  
  //perform correction on v2 gamma inclusive
  Double_t v2sig_arr[24];
  for(int i=0;i<24;i++){
    v2sig_arr[i] = ( v2T_arr[i] - Nback * v2CorrStr * v2R_arr[i] ) / Nsig;
    histo2040_corrected->SetBinContent(i+2,v2sig_arr[i]);
  }
  
  Double_t v2Dir_arr_old[24];
  Double_t v2Dir_arr_new[24];
//   Double_t v2C_constant = 1.08;
//   Double_t R = 1.07;
  for(int i=0;i<24;i++){
    //cocktail and R constant
//     v2Dir_arr_old[i] = ( R * v2T_arr[i] - v2C_constant * v2T_arr[i] ) / (R-1);
//     v2Dir_arr_new[i] = ( R * v2sig_arr[i] - v2C_constant * v2T_arr[i] ) / (R-1);
    //cocktail constant
//     v2Dir_arr_old[i] = ( R_arr[i] * v2T_arr[i] - v2C_constant * v2T_arr[i] ) / (R_arr[i]-1);
//     v2Dir_arr_new[i] = ( R_arr[i] * v2sig_arr[i] - v2C_constant * v2T_arr[i] ) / (R_arr[i]-1);
    //R constant
//     v2Dir_arr_old[i] = ( R * v2T_arr[i] - v2C_arr[i] ) / (R-1);
//     v2Dir_arr_new[i] = ( R * v2sig_arr[i] - v2C_arr[i] ) / (R-1);
    //none constant
    v2Dir_arr_old[i] = ( R_arr[i] * v2T_arr[i] - v2C_arr[i] ) / (R_arr[i]-1);
    v2Dir_arr_new[i] = ( R_arr[i] * v2sig_arr[i] - v2C_arr[i] ) / (R_arr[i]-1);
    histo2040_v2Dir_uncorrected->SetBinContent(i+2,v2Dir_arr_old[i]);
    histo2040_v2Dir_corrected->SetBinContent(i+2,v2Dir_arr_new[i]);
  }
  
  //setting bincontent for the cocktail and RGamma
  for(int i=0;i<24;i++){
    histo2040_cocktail->SetBinContent(i+2,v2C_arr[i]);
    histo2040_RGamma->SetBinContent(i+2,R_arr[i]);
  }
  
  //========================================================
  //Plotting and Saving
  //========================================================
  
  TLine *line = new TLine(0,1,8,1);
  TLine *line2 = new TLine(0,0,8,0);
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  c1->Divide(2,2,1E-10,1E-10);
  
  TLegend* leg = new TLegend(0.55,0.70,0.8,0.85);
  leg->SetTextSize(0.05);
  leg->SetHeader("v_{2} inclusive #gamma");
  leg->AddEntry(histo2040,"uncorrected","lp");
  leg->AddEntry(histo2040_corrected,"corrected","lp");
  leg->SetBorderSize(0);
  
  TLegend* leg2 = new TLegend(0.6,0.70,0.8,0.85);
  leg2->SetTextSize(0.05);
  leg2->AddEntry(histo2040,"v2_{#gamma}^{incl}","lp");
  leg2->AddEntry(histo2040_cocktail,"v2_{#gamma}^{cocktail}","lp");
  leg2->SetBorderSize(0);
  
  TLegend* leg3 = new TLegend(0.55,0.70,0.8,0.85);
  leg3->SetTextSize(0.05);
  leg3->AddEntry(histo2040_RGamma,"R_{#gamma}","");
  leg3->SetBorderSize(0);
  
  TLegend* leg4 = new TLegend(0.55,0.70,0.8,0.85);
  leg4->SetTextSize(0.05);
  leg4->SetHeader("v_{2} direct #gamma");
  leg4->AddEntry(histo2040,"uncorrected","lp");
  leg4->AddEntry(histo2040_corrected,"corrected","lp");
  leg4->SetBorderSize(0);
  
  c1->cd(1);
  histo2040->Draw();
  histo2040_cocktail->Draw("SAME");
  leg2->Draw();
  T1.DrawLatex(0.15, 0.34, "20-40% PbPb");
  T1.DrawLatex(0.15, 0.26, "#sqrt{s_{NN}}=2.76TeV");
  T1.DrawLatex(0.15, 0.18, "#gamma_{incl} #rightarrow e^{+}e^{-}");
//   T1.DrawLatex(0.15, 0.85, "v_{2} #gamma inclusive and cocktail");
  
  c1->cd(2);
  histo2040_RGamma->Draw();
  line->Draw();
//   leg3->Draw();
//   T1.DrawLatex(0.15, 0.85, "Double ratio R_{#gamma}");
  
  c1->cd(3);
  histo2040->Draw();
  histo2040_corrected->Draw("SAME");
  leg->Draw();
//   T1.DrawLatex(0.15, 0.85, "v_{2} inclusive #gamma");
  T1.DrawLatex(0.15, 0.32, Form("%s corrected",select.Data()));
  T1.DrawLatex(0.15, 0.25, Form("Strength = %2.2f",v2CorrStr));
  T1.DrawLatex(0.15, 0.18, Form("purity = %2.2f",purity));
  
  c1->cd(4);
  histo2040_v2Dir_uncorrected->Draw();
  histo2040_v2Dir_corrected->Draw("SAME");
  leg4->Draw();
  line2->Draw();
//   T1.DrawLatex(0.15, 0.85, "v_{2} direct #gamma");
  
  //Masking lower pt points
  MaskPoints(histo2040,0,11);
  MaskPoints(histo2040_corrected,0,11);
  MaskPoints(histo2040_v2Dir_uncorrected,0,11);
  MaskPoints(histo2040_v2Dir_corrected,0,11);
  MaskPoints(histo2040_cocktail,0,11);
  MaskPoints(histo2040_RGamma,0,11);
  
  //Masking higher pt points
  MaskPoints(histo2040,26,29);
  MaskPoints(histo2040_corrected,26,29);
  MaskPoints(histo2040_v2Dir_uncorrected,26,29);
  MaskPoints(histo2040_v2Dir_corrected,26,29);
  MaskPoints(histo2040_cocktail,26,29);
  MaskPoints(histo2040_RGamma,26,29);
  
  c1->SaveAs("ToyPlots/v2_gamma_inclusive_direct_toy_160302.pdf");
  
}

//Function that masks points from bin# start to bin# end
void MaskPoints(TH1F* histo, const int start, const int end){
  for(int i=start;i<end;i++){
    histo->SetBinContent(i,-10);
  }
}











