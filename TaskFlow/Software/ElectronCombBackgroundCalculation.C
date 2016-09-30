#include "DirectPhotonFlowFunctions.h"

void  ElectronCombBackgroundCalculation(
                                      TString CentralityLow = "20",
                                      TString CentralityHigh = "40"
){
  
  Double_t MassPion = 0.139570;
  Double_t MassElectron = 0.000501;
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  //open datafile with v2 background 1
  TFile* filebkg  = new TFile(Form("/home/mike/0_directphoton/13_DirectPhotonAnalysisCode/v2_background_%s%s.root",CentralityLow.Data(),CentralityHigh.Data()));
  
  TH1F*  histo_bkg_pipi = (TH1F*)filebkg->Get(Form("histo_v2_bkg_pipi_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  if(!histo_bkg_pipi) cout << "histo_bkg_pipi not found in filebkg!!" << endl;
  histo_bkg_pipi->SetMarkerColor(kOrange-3);
  histo_bkg_pipi->SetLineColor(kOrange-3);
  histo_bkg_pipi->SetMarkerStyle(21);
  
  TH1F*  histo_bkg_pie = (TH1F*)filebkg->Get(Form("histo_v2_bkg_pie_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  if(!histo_bkg_pie) cout << "histo_bkg_pie not found in filebkg!!" << endl;
  histo_bkg_pie->SetMarkerColor(kRed-3);
  histo_bkg_pie->SetLineColor(kRed-3);
  histo_bkg_pie->SetMarkerStyle(21);
  
  
  //========================================================
  //histogram cosmetics
  //========================================================
  
  TH1F*  histov2Empty = (TH1F*)Getv2HistStyle();
  histov2Empty->GetYaxis()->SetRangeUser(0,0.59);
  
  TH1F*  histo_bkg_ee = (TH1F*)histo_bkg_pipi->Clone();
  histo_bkg_ee->SetLineColor(kBlue+2);
  histo_bkg_ee->SetMarkerColor(kBlue+2);
  histo_bkg_ee->SetMarkerStyle(21);
  
  TF1* fitpipi = new TF1("fitpipi","pol3",1.0,6.0);
  TF1* fitpie = new TF1("fitpie","pol3",1.0,6.0);
  histo_bkg_pipi->Fit("fitpipi","R");
  histo_bkg_pie->Fit("fitpie","R");
  
  
  //========================================================
  //KET SCALE
  //========================================================
  
  Bool_t Usepipi = kFALSE;
  
  Int_t NBins = histo_bkg_ee->GetSize();
  for(Int_t i=0;i<NBins;i++){
    
    if(Usepipi){
      Float_t BinCenter = histo_bkg_ee->GetBinCenter(i);    
      Float_t KET = TMath::Sqrt(BinCenter*BinCenter+MassElectron*MassElectron*4.) - MassElectron*2.;
      Float_t ScaledPt = TMath::Sqrt((KET+MassPion*2.)*(KET+MassPion*2.)-MassPion*MassPion*4.);
      Float_t KETScaledValue = fitpipi->GetParameter(0)+fitpipi->GetParameter(1)*ScaledPt+fitpipi->GetParameter(2)*ScaledPt*ScaledPt+fitpipi->GetParameter(3)*ScaledPt*ScaledPt*ScaledPt;
      
      histo_bkg_ee->SetBinContent(i,KETScaledValue);
      
      cout << endl;
      cout << "BinCenter " << BinCenter << endl;
      cout << "KET " << KET << endl;
      cout << "ScaledPt " << ScaledPt << endl;
      cout << "UnscaledValue " << histo_bkg_pipi->GetBinContent(i) << endl;
      cout << "KETScaledValue " << KETScaledValue << endl;
      cout << endl;
    }else{
      Float_t BinCenter = histo_bkg_ee->GetBinCenter(i);    
      Float_t KET = TMath::Sqrt(BinCenter*BinCenter+MassElectron*MassElectron) - MassElectron;
      Float_t ScaledPt = TMath::Sqrt((KET+MassPion)*(KET+MassPion)-MassPion*MassPion);
      Float_t KETScaledValue = fitpie->GetParameter(0)+fitpie->GetParameter(1)*ScaledPt+fitpie->GetParameter(2)*ScaledPt*ScaledPt+fitpie->GetParameter(3)*ScaledPt*ScaledPt*ScaledPt;
      
      histo_bkg_ee->SetBinContent(i,KETScaledValue);
      
      cout << endl;
      cout << "BinCenter " << BinCenter << endl;
      cout << "KET " << KET << endl;
      cout << "ScaledPt " << ScaledPt << endl;
      cout << "UnscaledValue " << histo_bkg_pie->GetBinContent(i) << endl;
      cout << "KETScaledValue " << KETScaledValue << endl;
      cout << endl;
    }
    
    
  }
  
  //========================================================
  //Plotting and Saving
  //========================================================
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  SetProperMargins();
  histov2Empty->Draw();
  histo_bkg_pipi->Draw("SAME");
  histo_bkg_pie->Draw("SAME");
  histo_bkg_ee->Draw("SAME");
  
  
  
}