#include "DirectPhotonFlowFunctions.h"

void  Plot_InclusivePhotonFlow_dRdPhi(       
                                      TString CentralityLow = "20",
                                      TString CentralityHigh = "40",
                                      TString Cutnumber = "52400013_00200009007000008250400000"
                              ){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  //open datafile with v2 gamma inclusive uncorrected
  TString InputFileName1  = Form("/home/mike/3_PbPb_dirg/0_analysis/161205_v2_dRdPhi_rejection_55/Results_%s/InclusivePhotonv2_uncorrected_%s.root",Cutnumber.Data(),Cutnumber.Data());
  TFile* fileInclusive1  = new TFile(InputFileName1.Data());
  TH1F*  histoInclusive1 = (TH1F*)fileInclusive1->Get(Form("%s_tC55_2",Cutnumber.Data()));
  if(!histoInclusive1) cout << "histoInclusive1 not found in fileInclusive1!!" << endl;
  
  histoInclusive1->SetMarkerStyle(21);
  histoInclusive1->SetMarkerSize(1.2);
  histoInclusive1->SetMarkerColor(kRed+2);
  histoInclusive1->SetLineColor(kRed+2);
  
  TString InputFileName2  = Form("/home/mike/Dropbox/msas_ALICE/PbPb_v2dirg/standardcuts/Results_%s/InclusivePhotonv2_uncorrected_%s.root",Cutnumber.Data(),Cutnumber.Data());
  TFile* fileInclusive2  = new TFile(InputFileName2.Data());
  TH1F*  histoInclusive2 = (TH1F*)fileInclusive2->Get(Form("%s_tC55_2",Cutnumber.Data()));
  if(!histoInclusive2) cout << "histoInclusive2 not found in fileInclusive2!!" << endl;
  
  histoInclusive2->SetMarkerStyle(33);
  histoInclusive2->SetMarkerSize(1.2);
  histoInclusive2->SetMarkerColor(kGreen+2);
  histoInclusive2->SetLineColor(kGreen+2);
  
  TH1F* histoRatio = (TH1F*)histoInclusive1->Clone();
  histoRatio->Divide(histoRatio,histoInclusive2,1,1,"b");
  histoRatio->SetMarkerColor(kBlack);
  histoRatio->SetLineColor(kBlack);
  
  TF1* fit = new TF1("fit","pol0",0.9,5.4);
  histoRatio->Fit("fit","R");
  Double_t fitpar = fit->GetParameter(0);
  
  
  TFile* fileMC = new TFile("/home/mike/3_PbPb_dirg/1_data/161103_PbPb_v2_MC/GammaConvFlow_59.root");
  TList* list1   = new TList();
  list1     = (TList*)fileMC->Get("GammaConvV1_59_v2");
  cout << list1 << endl;
  TList* list2   = new TList();
  list2     = (TList*)list1->FindObject(Form("Cut Number %s",Cutnumber.Data()));
  cout << list2 << endl;
  TList* list3   = new TList();
  list3     = (TList*)list2->FindObject(Form("%s ESD histograms",Cutnumber.Data()));
  cout << list3 << endl;

  TH2F* hdRdPhi_all = (TH2F*)list3->FindObject("hdPhidRcandidates");
  TH2F* hdRdPhi_sigsig = (TH2F*)list3->FindObject("hdPhidRcandidates_MCsigsig");
  TH2F* hdRdPhi_sigbkg = (TH2F*)list3->FindObject("hdPhidRcandidates_MCbkgsig");
  TH2F* hdRdPhi_bkgbkg = (TH2F*)list3->FindObject("hdPhidRcandidates_MCbkgbkg");
  
  //========================================================
  //histogram cosmetics
  //========================================================
  
  hdRdPhi_all->SetTitle("");
  hdRdPhi_all->GetXaxis()->SetRangeUser(0,3.15);
  hdRdPhi_all->GetXaxis()->SetTitle("#Delta #phi");
  hdRdPhi_all->GetYaxis()->SetRangeUser(1,299);
  hdRdPhi_all->GetYaxis()->SetTitle("#Delta R");
  hdRdPhi_all->GetXaxis()->SetTitleSize(0.06);
  hdRdPhi_all->GetYaxis()->SetTitleSize(0.06);
  hdRdPhi_all->GetXaxis()->SetTitleOffset(0.8);
  hdRdPhi_all->GetYaxis()->SetTitleOffset(0.85);
  
  hdRdPhi_sigsig->SetTitle("");
  hdRdPhi_sigsig->GetXaxis()->SetRangeUser(0,0.5);
  hdRdPhi_sigsig->GetXaxis()->SetTitle("#Delta #phi");
  hdRdPhi_sigsig->GetYaxis()->SetRangeUser(1,49);
  hdRdPhi_sigsig->GetYaxis()->SetTitle("#Delta R");
  hdRdPhi_sigsig->GetXaxis()->SetTitleSize(0.06);
  hdRdPhi_sigsig->GetYaxis()->SetTitleSize(0.06);
  hdRdPhi_sigsig->GetXaxis()->SetTitleOffset(0.8);
  hdRdPhi_sigsig->GetYaxis()->SetTitleOffset(0.85);
  
  hdRdPhi_sigbkg->SetTitle("");
  hdRdPhi_sigbkg->GetXaxis()->SetRangeUser(0,0.5);
  hdRdPhi_sigbkg->GetXaxis()->SetTitle("#Delta #phi");
  hdRdPhi_sigbkg->GetYaxis()->SetRangeUser(1,49);
  hdRdPhi_sigbkg->GetYaxis()->SetTitle("#Delta R");
  hdRdPhi_sigbkg->GetXaxis()->SetTitleSize(0.06);
  hdRdPhi_sigbkg->GetYaxis()->SetTitleSize(0.06);
  hdRdPhi_sigbkg->GetXaxis()->SetTitleOffset(0.8);
  hdRdPhi_sigbkg->GetYaxis()->SetTitleOffset(0.85);
  
  hdRdPhi_bkgbkg->SetTitle("");
  hdRdPhi_bkgbkg->GetXaxis()->SetRangeUser(0,0.5);
  hdRdPhi_bkgbkg->GetXaxis()->SetTitle("#Delta #phi");
  hdRdPhi_bkgbkg->GetYaxis()->SetRangeUser(1,49);
  hdRdPhi_bkgbkg->GetYaxis()->SetTitle("#Delta R");
  hdRdPhi_bkgbkg->GetXaxis()->SetTitleSize(0.06);
  hdRdPhi_bkgbkg->GetYaxis()->SetTitleSize(0.06);
  hdRdPhi_bkgbkg->GetXaxis()->SetTitleOffset(0.8);
  hdRdPhi_bkgbkg->GetYaxis()->SetTitleOffset(0.85);
  
  TH1F*  histoInclusiveEmpty = (TH1F*)Getv2HistStyle();
  histoInclusiveEmpty->GetYaxis()->SetRangeUser(0,0.49);
  
  TH1F*  histoInclusiveEmptyRatio = (TH1F*)Getv2HistStyle();
  histoInclusiveEmptyRatio->GetYaxis()->SetRangeUser(0.93,1.07);
  histoInclusiveEmptyRatio->GetYaxis()->SetTitle("Ratio");
  
  //========================================================
  //integrals of 2D histograms
  //========================================================
  
  Double_t integral_sigsig_all     = hdRdPhi_sigsig->Integral(hdRdPhi_sigsig->GetXaxis()->FindBin(0.0),hdRdPhi_sigsig->GetXaxis()->FindBin(3.14),hdRdPhi_sigsig->GetYaxis()->FindBin(0.0),hdRdPhi_sigsig->GetYaxis()->FindBin(300));
  Double_t integral_sigsig  = hdRdPhi_sigsig->Integral(hdRdPhi_sigsig->GetXaxis()->FindBin(0.0),hdRdPhi_sigsig->GetXaxis()->FindBin(0.03),hdRdPhi_sigsig->GetYaxis()->FindBin(0.0),hdRdPhi_sigsig->GetYaxis()->FindBin(7.5));
  Double_t integral_sigbkg  = 2*hdRdPhi_sigbkg->Integral(hdRdPhi_sigbkg->GetXaxis()->FindBin(0.0),hdRdPhi_sigbkg->GetXaxis()->FindBin(0.03),hdRdPhi_sigbkg->GetYaxis()->FindBin(0.0),hdRdPhi_sigbkg->GetYaxis()->FindBin(7.5));
  Double_t integral_bkgbkg  = hdRdPhi_bkgbkg->Integral(hdRdPhi_bkgbkg->GetXaxis()->FindBin(0.0),hdRdPhi_bkgbkg->GetXaxis()->FindBin(0.03),hdRdPhi_bkgbkg->GetYaxis()->FindBin(0.0),hdRdPhi_bkgbkg->GetYaxis()->FindBin(7.5));
  Double_t integral_sigbkg_all  = 2*hdRdPhi_sigbkg->Integral(hdRdPhi_sigbkg->GetXaxis()->FindBin(0.0),hdRdPhi_sigbkg->GetXaxis()->FindBin(3.14),hdRdPhi_sigbkg->GetYaxis()->FindBin(0.0),hdRdPhi_sigbkg->GetYaxis()->FindBin(300));
  Double_t integral_bkgbkg_all  = hdRdPhi_bkgbkg->Integral(hdRdPhi_bkgbkg->GetXaxis()->FindBin(0.0),hdRdPhi_bkgbkg->GetXaxis()->FindBin(3.14),hdRdPhi_bkgbkg->GetYaxis()->FindBin(0.0),hdRdPhi_bkgbkg->GetYaxis()->FindBin(300));
  
  cout << endl << endl;
  cout << " integral_sigsig = " << integral_sigsig << endl;
  cout << " integral_sigbkg = " << integral_sigbkg << endl;
  cout << " integral_bkgbkg = " << integral_bkgbkg << endl << endl;
  cout << " all(S+B) = " << integral_sigsig+integral_sigbkg+integral_bkgbkg << endl;
  cout << " S = " << integral_sigsig << endl;
  cout << " B = " << integral_sigbkg+integral_bkgbkg << endl;
  cout << " S/B = " << integral_sigsig/(integral_sigbkg+integral_bkgbkg) << endl;
  cout << " S(cut)/S(all) = " << 100*integral_sigsig/integral_sigsig_all << " %" << endl << endl;
  
  //========================================================
  //Plotting and Saving
  //========================================================
  TLatex T1;
  T1.SetTextSize(0.04);
  T1.SetTextAlign(12);
  T1.SetNDC();
  
  TLatex T2;
  T2.SetTextSize(0.1);
  T2.SetTextAlign(12);
  T2.SetNDC();
  
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  TCanvas* c1 = new TCanvas("c1","",400,800);
  c1->Divide(1,2,0.001,0.001);
  
  TLegend* leg = new TLegend(0.16,0.75,0.4,0.95);
  leg->SetTextSize(0.04);
//   leg->SetHeader("v_{2}^{#gamma,inclusive}");
  leg->AddEntry(histoInclusive1,"v_{2}^{#gamma,incl} unc.","lp");
  leg->AddEntry(histoInclusive2,"v_{2}^{#gamma,incl} unc. dRdPhi cut ","lp");
//   leg->AddEntry(histoInclusive3,"-14<K<-12","lp");
  leg->SetBorderSize(0);

  c1->cd(1);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoInclusive1->Draw("SAME");
  histoInclusive2->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  leg->Draw();
  
  c1->cd(2);
  SetProperMargins();
  histoInclusiveEmptyRatio->Draw();
  histoRatio->Draw("SAME");
  T1.DrawLatex(0.25, 0.85, Form("fit par0 = %2.3f", fitpar));
  
  //Masking lower pt points
  MaskPoints(histoInclusive1,0,10);
  MaskPoints(histoInclusive2,0,10);
  MaskPoints(histoRatio,0,10);
  
  
  //Masking higher pt points
  MaskPoints(histoInclusive1,26,29);
  MaskPoints(histoInclusive2,26,29);
  MaskPoints(histoRatio,26,29);
  
  TLine *l1 = new TLine(0,50,0.5,50);
  l1->SetLineWidth(3);
  TLine *l2 = new TLine(0.5,0,0.5,50);
  l2->SetLineWidth(3);
  
  TCanvas* c2 = new TCanvas("c2","",800,800);
  c2->Divide(2,2,0.001,0.001);
  c2->cd(1);
  c2->SetLogz();
  SetProperMargins();
  hdRdPhi_all->Draw("COLZ");
  T2.DrawLatex(0.7, 0.8, "ALL");
  l1->Draw();
  l2->Draw();
  c2->cd(2);
  c2->SetLogz();
  SetProperMargins();
  hdRdPhi_sigsig->Draw("COLZ");
  T2.DrawLatex(0.7, 0.8, "S + S");
  c2->cd(3);
  c2->SetLogz();
  SetProperMargins();
  hdRdPhi_sigbkg->Draw("COLZ");
  T2.DrawLatex(0.7, 0.8, "S + B");
  c2->cd(4);
  c2->SetLogz();
  SetProperMargins();
  hdRdPhi_bkgbkg->Draw("COLZ");
  T2.DrawLatex(0.7, 0.8, "B + B");
  
  

  gSystem->mkdir("Results");
  gSystem->mkdir("Results/v2GammaInclusiveComparison_dRdPhi");
  c1->SaveAs(Form("Results/v2GammaInclusiveComparison_dRdPhi/%s_%s_v2_Gamma_Inclusive_dRdPhi.eps",CentralityLow.Data(),CentralityHigh.Data()));
  c2->SaveAs(Form("Results/v2GammaInclusiveComparison_dRdPhi/%s_%s_dRdPhi.eps",CentralityLow.Data(),CentralityHigh.Data()));
  
}