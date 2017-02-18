#include "DirectPhotonFlowFunctions.h"

void  InclusivePhotonFlowCorrection(  Int_t Trainconfig = 55,
                                      TString CentralityLow = "0",
                                      TString CentralityHigh = "20",
                                      TString Cutnumber = "50200013_00200009007000008250400000"
){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  TString InputFileNameUncorrected  = Form("/home/mike/3_PbPb_dirg/0_analysis/170216_v2_final_systematics/Results_%s/InclusivePhotonv2_uncorrected_%s.root",Cutnumber.Data(),Cutnumber.Data());
  TString InputFileNameBackground   = Form("/home/mike/3_PbPb_dirg/0_analysis/170216_v2_final_systematics/Results_%s/v2_background_%s%s_%s.root",Cutnumber.Data(),CentralityLow.Data(),CentralityHigh.Data(),Cutnumber.Data());
  TString InputFileNamePurity       = Form("/home/mike/3_PbPb_dirg/0_analysis/170216_v2_final_systematics/purity_studies_%s/Purity_InclusivePhotonSample_%s%s_%s_Data LHC10h.root",Cutnumber.Data(),CentralityLow.Data(),CentralityHigh.Data(),Cutnumber.Data());
  
  //open datafile with v2 gamma inclusive uncorrected
  TFile* fileInclusive  = new TFile(InputFileNameUncorrected.Data());
  TH1F*  histoInclusive = (TH1F*)fileInclusive->Get(Form("%s_tC%i_2",Cutnumber.Data(),Trainconfig));
  if(!histoInclusive) cout << "histoInclusive not found in fileInclusive!!" << endl;
  histoInclusive->SetMarkerColor(kRed+2);
  histoInclusive->SetLineColor(kRed+2);
  histoInclusive->SetMarkerStyle(20);
  //open datafile with v2 gamma inclusive background
  TFile* fileCorrection = new TFile(InputFileNameBackground.Data());
  TH1F*  histoCorrection = (TH1F*)fileCorrection->Get(Form("histo_v2_bkg_total_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  if(!histoCorrection) cout << "histoCorrection not found in fileCorrection!!" << endl;
  histoCorrection->SetMarkerColor(kBlue+2);
  histoCorrection->SetLineColor(kBlue+2);
  //open datafile with purity
  TFile* filePurity     = new TFile(InputFileNamePurity.Data());
  TGraphErrors*  GraphPurity = (TGraphErrors*)filePurity->Get(Form("Gamma_Purity_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  if(!GraphPurity) cout << "GraphPurity not found in filePurity!!" << endl;
  
  //========================================================
  //histogram cosmetics
  //========================================================
  
  TH1F*  histoInclusiveEmpty = (TH1F*)Getv2HistStyle();
  histoInclusiveEmpty->GetYaxis()->SetRangeUser(0,0.34);
  
  TH1F*  GraphPurityEmpty = (TH1F*)GetPurityHistStyle();
  
  TH1F*  histoInclusiveRatioEmpty = (TH1F*)Getv2InclRatioHistStyle();
  
  
  cout << "//========================================================" << endl;
  cout << "//Calculating v2 inclusive correction" << endl;
  cout << "//========================================================" << endl;
  
  TH1F*  histoInclusiveCorrected = (TH1F*)histoInclusive->Clone();
  histoInclusiveCorrected->SetLineColor(kGreen+2);
  histoInclusiveCorrected->SetMarkerColor(kGreen+2);
  histoInclusiveCorrected->SetMarkerStyle(24);
  
  Int_t NBins = histoInclusiveCorrected->GetSize();
  for(Int_t i=0;i<NBins;i++){
    Float_t BinCenter = histoInclusive->GetBinCenter(i);
    Float_t UncorInclValue = histoInclusive->GetBinContent(i);
    Float_t Purity = GraphPurity->Eval(BinCenter);
    Float_t BackgroundValue = histoCorrection->Interpolate(BinCenter);
    Float_t CorInclValue = ( UncorInclValue - (1-Purity) * BackgroundValue ) / Purity;
    
    cout << endl;
    cout << "Correcting pt = " << BinCenter << " GeV/c" << endl;
    cout << "Purity = " << Purity << endl;
    cout << "Uncorrected v2 incl value = " << UncorInclValue << endl;
    cout << "Backgound   v2 incl value = " << BackgroundValue << endl;
    cout << "Corrected value      = " << CorInclValue << endl;
    cout << "Difference -> " << 100*(CorInclValue-UncorInclValue)/UncorInclValue << " % " << endl;
    cout << endl;
    if(CorInclValue>-10 && CorInclValue<10){
      histoInclusiveCorrected->SetBinContent(i,CorInclValue);
    }else{
      histoInclusiveCorrected->SetBinContent(i,-10);
    }
  }
  
  histoInclusive->Sumw2();
  histoInclusiveCorrected->Sumw2();
  TH1F*  histoInclusiveRatio = (TH1F*)histoInclusive->Clone();
  histoInclusiveRatio->Divide(histoInclusiveRatio,histoInclusiveCorrected,1,1,"B");
  histoInclusiveRatio->GetXaxis()->SetTitle("#it{p}_{T}(GeV/c)");
  histoInclusiveRatio->SetMarkerColor(kBlack);
  histoInclusiveRatio->SetLineColor(kBlack);
  
  
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
  
  TLine *line = new TLine(0,1,8,1);
  TLine *line2 = new TLine(0,0,8,0);
  
  TLegend* leg = new TLegend(0.16,0.75,0.6,0.86);
  leg->SetTextSize(0.05);
//   leg->SetHeader("v_{2} inclusive #gamma");
  leg->AddEntry(histoInclusive,"uncorrected","lp");
  leg->AddEntry(histoInclusiveCorrected,"corrected","lp");
  leg->SetBorderSize(0);
  
  TLegend* leg2 = new TLegend(0.16,0.75,0.6,0.86);
  leg2->SetTextSize(0.05);
//   leg2->SetHeader("v_{2} inclusive #gamma");
  leg2->AddEntry(histoInclusive,"Uncorrected","lp");
  leg2->AddEntry(histoCorrection,"Background","lp");
  leg2->SetBorderSize(0);
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  c1->Divide(2,2,1E-10,1E-10);
  
  c1->cd(1);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoInclusive->Draw("SAME");
  histoCorrection->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  leg2->Draw();
  
  c1->cd(2);
  SetProperMargins();
  GraphPurityEmpty->Draw();
  GraphPurity->Draw("SAME P");
  
  c1->cd(3);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoInclusive->Draw("SAME");
  histoInclusiveCorrected->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  leg->Draw();
  
  c1->cd(4);
  SetProperMargins();
  histoInclusiveRatioEmpty->Draw();
  histoInclusiveRatio->Draw("SAME");
  T1.DrawLatex(0.15, 0.92, "v^{#gamma,incl,uncorrected}_{2} / v^{#gamma,incl,corrected}_{2}");
  line->Draw();
  
  //Masking lower pt points
  MaskPoints(histoInclusive,0,10);
  MaskPoints(histoCorrection,0,10);
  MaskPoints(histoInclusiveCorrected,0,10);
  MaskPoints(histoInclusiveRatio,0,10);
//   
//   //Masking higher pt points
  MaskPoints(histoInclusive,26,29);
  MaskPoints(histoCorrection,26,29);
  MaskPoints(histoInclusiveCorrected,26,29);
//   MaskPoints(GraphPurity,18,29);
  MaskPoints(histoInclusiveRatio,26,29);
  
  gSystem->mkdir(Form("Results_%s",Cutnumber.Data()));
  gSystem->mkdir(Form("Results_%s/v2GammaInclusiveCorrected",Cutnumber.Data()));
  c1->SaveAs(Form("Results_%s/v2GammaInclusiveCorrected/v2GammaIncl_Correction_%s%s.eps",Cutnumber.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  
  TFile *InclusivePhotonv2_File = new TFile(Form("Results_%s/InclusivePhotonv2_Corrected_%i_%s.root",Cutnumber.Data(),Trainconfig,Cutnumber.Data()),"RECREATE");
  histoInclusive->Write(Form("v2GammaIncl_uncorrected_%s%s_tC%i",CentralityLow.Data(),CentralityHigh.Data(),Trainconfig));
  histoCorrection->Write(Form("v2GammaIncl_background_%s%s_tC%i",CentralityLow.Data(),CentralityHigh.Data(),Trainconfig));
  GraphPurity->Write(Form("v2GammaIncl_purity_%s%s_tC%i",CentralityLow.Data(),CentralityHigh.Data(),Trainconfig));
  histoInclusiveCorrected->Write(Form("v2GammaIncl_corrected_%s%s_tC%i",CentralityLow.Data(),CentralityHigh.Data(),Trainconfig));
  InclusivePhotonv2_File->Close();
  
}