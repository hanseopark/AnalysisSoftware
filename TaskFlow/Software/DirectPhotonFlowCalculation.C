#include "DirectPhotonFlowFunctions.h"

void  DirectPhotonFlowCalculation(
                                      Int_t Trainconfig = 55,
                                      TString CentralityLow = "0",
                                      TString CentralityHigh = "20",
                                      TString Cutnumber = "50200013_00200009007000008250400000"
                                 ){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  TFile* fileInclusive  = new TFile(Form("/home/mike/0_directphoton/0_analysis/160714_PbPb_systematics/Results_%s/InclusivePhotonv2_Corrected_%i_%s.root",Cutnumber.Data(),Trainconfig,Cutnumber.Data()));
  TH1F*  histoInclusive = (TH1F*)fileInclusive->Get(Form("v2GammaIncl_corrected_%s%s_tC%i",CentralityLow.Data(),CentralityHigh.Data(),Trainconfig));
//   TFile* fileInclusive  = new TFile("/home/mike/0_directphoton/13_DirectPhotonAnalysisCode/InclusivePhotonv2_uncorrected_53.root");
//   TH1F*  histoInclusive = (TH1F*)fileInclusive->Get(Form("%s%s_tC%i_0",CentralityLow.Data(),CentralityHigh.Data(),Trainconfig));
  if(!histoInclusive) cout << "histoInclusive not found in fileInclusive!!" << endl;
  histoInclusive->SetLineColor(kBlack);
  histoInclusive->SetMarkerColor(kBlack);
  histoInclusive->SetMarkerStyle(20);
  
  TFile* fileRGamma = new TFile("/home/mike/0_directphoton/0_analysis/160203_PbPb_v2_10h/Gamma_CombResults_PbPb_2.76TeV.root");
  TDirectory* dir1 = new TDirectory();
  dir1 = (TDirectory*)fileRGamma->Get(Form("Gamma_PbPb_2.76TeV_%s-%s%%",CentralityLow.Data(),CentralityHigh.Data()));
  TH1F*  histoRGamma = (TH1F*)dir1->Get("hDR_PCM_StatErr");
  if(!histoRGamma) cout << "histoRGamma not found in fileRGamma!!" << endl;
  histoRGamma->SetMarkerStyle(20);
  histoRGamma->SetMarkerSize(histoInclusive->GetMarkerSize());
  histoRGamma->SetMarkerColor(kBlack);
  histoRGamma->SetLineColor(kBlack);
  
  TFile* fileCocktail = new TFile("/home/mike/0_directphoton/13_DirectPhotonAnalysisCode/CocktailV2.root");
  //cen3 = 2040
  //cen6 = 020
  //cen8 = 4080
  TString CocktailString;
  if(CentralityLow.CompareTo("0")==0) CocktailString = "cen6";
  if(CentralityLow.CompareTo("20")==0) CocktailString = "cen3";
  if(CentralityLow.CompareTo("40")==0) CocktailString = "cen8";
  TH1F*  histoCocktail = (TH1F*)fileCocktail->Get(Form("v2gammaCocktail_%s",CocktailString.Data()));
  if(!histoCocktail) cout << "histoCocktail not found in fileCocktail!!" << endl;
  histoCocktail->SetLineColor(kGreen+2);
  histoCocktail->SetMarkerColor(kGreen+2);
  histoCocktail->SetMarkerStyle(21);
  
  //========================================================
  //histogram cosmetics
  //========================================================
  
  TH1F*  histoInclusiveEmpty = (TH1F*)Getv2HistStyle();
  
  TH1F*  histoRGammaEmpty = (TH1F*)GetRGammaHistStyle();
  
  TH1F* histoDirectEmpty = (TH1F*)histoInclusiveEmpty->Clone();
  histoDirectEmpty->GetYaxis()->SetRangeUser(-0.2,0.39);
  
  cout << "//========================================================" << endl;
  cout << "//Calculating v2 Direct" << endl;
  cout << "//========================================================" << endl;
  
  TH1F*  histoDirect = (TH1F*)histoInclusive->Clone();
  histoDirect->SetMarkerStyle(20);
  
  Int_t NBins = histoInclusive->GetSize();
  for(Int_t i=0;i<NBins;i++){
    Float_t BinCenter = histoInclusive->GetBinCenter(i);
    Float_t v2InclValue = 0;
    Float_t RGammaValue = 0;
    Float_t v2CocktailValue = 0;
    
    v2InclValue = histoInclusive->GetBinContent(i);
    RGammaValue = histoRGamma->Interpolate(BinCenter);
    v2CocktailValue = histoCocktail->Interpolate(BinCenter);
    
    Float_t v2DirectValue = (RGammaValue*v2InclValue - v2CocktailValue) / (RGammaValue - 1);
    
    Float_t v2InclError = histoInclusive->GetBinError(i);
    Float_t RGammaError = histoRGamma->GetBinError(histoRGamma->FindBin(BinCenter));
    Float_t v2CocktailError = histoCocktail->GetBinError(histoCocktail->FindBin(BinCenter));
    
    Float_t RGammaFactor = (v2CocktailValue-v2InclValue)/((RGammaValue-1)*(RGammaValue-1));
    Float_t InclFactor = RGammaValue/(RGammaValue-1);
    Float_t CocktailFactor = 1/(RGammaValue-1);
    
    Float_t v2DirectStatError = TMath::Sqrt(RGammaFactor*RGammaFactor*RGammaError*RGammaError+InclFactor*InclFactor*v2InclError*v2InclError+CocktailFactor*CocktailFactor*v2CocktailError*v2CocktailError);
    
    cout << endl;
    cout << "Calculating v2 Direct for pt = " << BinCenter << " GeV/c" << endl;
    cout << "RGamma       = " << RGammaValue << endl;
    cout << "v2 inclusive = " << v2InclValue << endl;
    cout << "v2 Cocktail  = " << v2CocktailValue << endl;
    cout << "v2 Direct    = " << v2DirectValue << " error = " << v2DirectStatError << endl;
    cout << endl;
    histoDirect->SetBinContent(i,v2DirectValue);
    histoDirect->SetBinError(i,v2DirectStatError);
  }
  
  histoDirect->Sumw2();
    
  //========================================================
  //Plotting and Saving
  //========================================================
  
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  TLine *line = new TLine(0,1,6.2,1);
  line->SetLineStyle(2);
  line->SetLineColor(kGray+1);
  TLine *line2 = new TLine(0,0,6.2,0);
  line2->SetLineStyle(2);
  line2->SetLineColor(kGray+1);
  
  TLegend* leg1 = new TLegend(0.16,0.15,0.6,0.26);
  leg1->SetTextSize(0.05);
  leg1->AddEntry(histoInclusive,"v_{2}^{#gamma,inclusive}","p");
  leg1->AddEntry(histoCocktail,"v_{2}^{#gamma,decay}","p");
  leg1->SetBorderSize(0);
  
  TLegend* leg2 = new TLegend(0.16,0.15,0.6,0.26);
  leg2->SetTextSize(0.05);
  leg2->SetHeader("v_{2}^{#gamma,direct}");
  leg2->SetBorderSize(0);
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoInclusive->Draw("SAME");
  histoCocktail->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  leg1->Draw();
  
  TCanvas* c2 = new TCanvas("c2","",800,800);
  SetProperMargins();
  histoRGammaEmpty->Draw();
  histoRGamma->Draw("SAME");
  DrawRGammaInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  line->Draw();
  
  TCanvas* c3 = new TCanvas("c3","",800,800);
  SetProperMargins();
  histoDirectEmpty->Draw();
  histoDirect->Draw("SAME");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  line2->Draw();
  
  //Masking lower pt points
  MaskPoints(histoInclusive,0,histoInclusive->FindBin(1));
  MaskPoints(histoDirect,0,histoDirect->FindBin(1));
  MaskPoints(histoCocktail,0,histoCocktail->FindBin(1));
  //Masking higher pt points
  MaskPoints(histoInclusive,histoInclusive->FindBin(6.3),histoInclusive->FindBin(10));
  MaskPoints(histoDirect,histoDirect->FindBin(6.3),histoDirect->FindBin(10));
  MaskPoints(histoCocktail,histoCocktail->FindBin(6.3),histoCocktail->FindBin(10));
  MaskPoints(histoRGamma,histoRGamma->FindBin(6.3),histoRGamma->FindBin(10));
  
  gSystem->mkdir(Form("Results_%s",Cutnumber.Data()));
  gSystem->mkdir(Form("Results_%s/v2GammaDirect",Cutnumber.Data()));
  c1->SaveAs(Form("Results_%s/v2GammaDirect/v2_Incl_Cocktail_%s%s.eps",Cutnumber.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  c2->SaveAs(Form("Results_%s/v2GammaDirect/RGamma_%s%s.eps",Cutnumber.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  c3->SaveAs(Form("Results_%s/v2GammaDirect/v2_DirectGamma_%s%s.eps",Cutnumber.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  
  TFile *output_File = new TFile(Form("Results_%s/PCM_DirectPhotonFlow_%s%s_%s.root",Cutnumber.Data(), CentralityLow.Data(), CentralityHigh.Data(),Cutnumber.Data()),"RECREATE");
  histoInclusive->Write("histoInclusive");
  histoCocktail->Write("histoCocktail");
  histoDirect->Write("histoDirect");
  histoRGamma->Write("histoRGamma");
  output_File->Close();
  
}