#include "Photonv2Spectra_V1.h"

void Photonv2Spectra_v2(
                        TString InputFileName = "AnalysisResults53.root",
                        TString CentralityLow = "20",
                        TString CentralityHigh = "40",
                        Int_t Trainconfig = 53,
                        TString SaveName = "53"){
  
  //InputFileName: standard flow output
  //CentralityLow
  //CentralityHigh
  //Centrality list: 05,510,1020,2030,3040,020,040,4060,6080
  //Trainconfig: needed to get the right directory in the InputFile
  
  //======================================================================================================================================
  //======================================================================================================================================
  
  StyleSettingsThesis();
  SetPlotStyle();
  libs();
  
  //======================================================================================================================================
  //======================================================================================================================================
  
  bool reb = kTRUE;
  TString dirdetectSP = "SP_V0";
  
  Color_t colorIncl = GetColorDefaultColor("PbPb_2.76TeV","",Form("%s-%s%%",CentralityLow.Data(),CentralityHigh.Data()));
  Marker_t markerIncl = GetDefaultMarkerStyle("PbPb_2.76TeV","",Form("%s-%s%%",CentralityLow.Data(),CentralityHigh.Data()));
  //Define empty histogram with correct settings
  Double_t xAxisFin[13] = {0, 1, 1.5, 2, 2.5, 3, 4, 5, 8};
  TH1D *emptyhisto = new TH1D("emptyhisto","emptyhisto",8,xAxisFin);
  SetStyleHistoTH1ForGraphs(emptyhisto,"#it{p}_{T} (GeV/#it{c})" , "#it{v}_{2}", 0.035, 0.04, 0.035, 0.04,1,1.2);
  emptyhisto->GetYaxis()->SetRangeUser(-0.07,0.51);
  emptyhisto->GetXaxis()->SetTitleFont(42);
  emptyhisto->GetYaxis()->SetTitleFont(42);
    
  //======================================================================================================================================
  //======================================================================================================================================
    
  TString fileSPv0Redo = InputFileName.Data();
  cout << Form("SPVZEROQa_in_%s%scc_tC%i_v2SPv2Qa",CentralityLow.Data(),CentralityHigh.Data(),Trainconfig) << endl;
  TString directoryQaRedo = Form("SPVZEROQa_in_%s%scc_tC%i_v2_0SPv2Qa",CentralityLow.Data(),CentralityHigh.Data(),Trainconfig);
  TString directoryQbRedo = Form("SPVZEROQb_in_%s%scc_tC%i_v2_0SPv2Qb",CentralityLow.Data(),CentralityHigh.Data(),Trainconfig);

  TH1D *v2InclGamma        = (TH1D*)IncluSPv0(fileSPv0Redo,directoryQaRedo,directoryQbRedo,dirdetectSP,reb);
  TH1D *v2InclGammaSpectra = (TH1D*)Spectra(fileSPv0Redo,directoryQaRedo,directoryQbRedo,dirdetectSP,reb);

  TCanvas *f1 = new TCanvas("c1","c1",200,10,1000,950);
  DrawGammaCanvasSettings(f1,0.1,0.02,0.02,0.1);
  DrawGammaSetMarker(v2InclGamma, markerIncl, 1, colorIncl , colorIncl);
  DrawGammaSetMarker(v2InclGammaSpectra, markerIncl, 1, colorIncl , colorIncl);
  
  emptyhisto->Draw();
  v2InclGamma->Draw("ePSAME");
  
  TLegend *Legend1 = new TLegend(0.65, 0.65, 0.95, 0.8);
  Legend1 -> SetFillColor(0);
  Legend1 -> SetTextSize(0.03);
  Legend1 -> SetTextFont(42);
  Legend1 -> SetTextColor(1);
  Legend1 -> SetBorderSize(1);
  Legend1->SetHeader("v_{2} #gamma inclusive");
  Legend1->AddEntry(v2InclGamma,"Mike, SP-V0","PL");
  Legend1->SetBorderSize(0);
  
  TPaveText* pt1 = new TPaveText(0.4,0.85,0.45,0.95,"brNDC");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetTextFont(42);
  pt1->SetTextSize(0.04);
  TText* text1 = pt1->AddText(Form("%s-%s%% Pb-Pb, #sqrt{#it{s}_{NN}} = 2.76 TeV",CentralityLow.Data(),CentralityHigh.Data()));
  
  pt1->Draw("same");
  Legend1->Draw("same");
    
  //======================================================================================================================================
  //======================================================================================================================================
  
  TCanvas *f2 = new TCanvas();
  v2InclGammaSpectra->SetTitle("");
  v2InclGammaSpectra->Draw("ePSAME");
  
  TLegend *Legend2 = new TLegend(0.7, 0.55, 0.9, 0.85);
  Legend2 -> SetFillColor(0);
  Legend2 -> SetTextSize(0.03);
  Legend2 -> SetTextFont(42);
  Legend2 -> SetTextColor(1);
  Legend2 -> SetBorderSize(1);
  Legend2->SetHeader("inclusive #gamma Spectra");
  Legend2->AddEntry(v2InclGammaSpectra,"0-5%cc","PL");
  Legend2->SetBorderSize(0);
  
  TPaveText* pt2 = new TPaveText(0.15,0.85,0.75,0.95,"brNDC");
  pt2->SetBorderSize(0);
  pt2->SetFillStyle(0);
  pt2->SetTextFont(42);
  pt2->SetTextSize(0.04);
  TText *text2 = pt2->AddText("Pb-Pb, #sqrt{#it{s}_{NN}} = 2.76 TeV");
  
  pt2->Draw();
  Legend2->Draw("same");
  gPad->SetLogy();
  
  //======================================================================================================================================
  //======================================================================================================================================
  
  TFile *InclusivePhotonv2_File = new TFile("InclusivePhotonv2_SmallCC.root","UPDATE");
  v2InclGamma->Write(Form("%s%s_tC%s",CentralityLow.Data(),CentralityHigh.Data(),SaveName.Data()));
  InclusivePhotonv2_File->Close();
  
  TFile *InclusivePhotonv2_spectra_File = new TFile("InclusivePhotonv2_SmallCC_unmerged_spectra.root","UPDATE");
  v2InclGammaSpectra->Write(Form("%s%s_tC%s_Spectra",CentralityLow.Data(),CentralityHigh.Data(),SaveName.Data()));
  InclusivePhotonv2_spectra_File->Close();

  gSystem->mkdir("plots");
  f1->SaveAs(Form("plots/v2_gammaincl_%s%s_%s.eps",CentralityLow.Data(),CentralityHigh.Data(),SaveName.Data()));
  
  //======================================================================================================================================
  //======================================================================================================================================
  
}