#include "DirectPhotonFlowFunctions.h"

TGraphAsymmErrors* CalculateSymmetricSystematics(TH1F* histDefault, TH1F* histvar);
TGraphAsymmErrors* CalculateSymmetricSystematics(TH1F* histDefault, TH1F* histvar1, TH1F* histvar2);
TGraphAsymmErrors* CalculateTotalSystematics(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, TGraphAsymmErrors* g3, TGraphAsymmErrors* g4, TGraphAsymmErrors* g5, TGraphAsymmErrors* g6, TGraphAsymmErrors* g7);
TGraphAsymmErrors* CalculateTotalError(TH1F* h1, TGraphAsymmErrors* g1);
void SmoothSystematicError(TGraphAsymmErrors* g1,TString NameSyst);

void  Systematics(
                                      Int_t Trainconfig = 55,
                                      TString CentralityLow = "0",
                                      TString CentralityHigh = "20",
                                      Bool_t SetSystSmoothing = kTRUE
                                 ){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  TString DefaultCut;
  Int_t SystTrainConfig;
  if(CentralityLow.CompareTo("0")==0)  {DefaultCut = "50200013_00200009007000008250400000"; SystTrainConfig = 60;}
  if(CentralityLow.CompareTo("20")==0) {DefaultCut = "52400013_00200009007000008250400000"; SystTrainConfig = 63;}
  if(CentralityLow.CompareTo("40")==0) {DefaultCut = "54800013_00200009007000008250400000"; SystTrainConfig = 66;}
  
  TString SystNames[7] = {"eta","R","minPt","qt","chi2","psipair","cosp"};
  
  TString SystVariations[12];
  if(CentralityLow.CompareTo("0")==0){
    SystVariations[ 0] = "50200013_04200009007000008250400000"; //eta cut: |eta| <0.75
    SystVariations[ 1] = "50200013_00500009007000008250400000"; //minR= 10 maxR = 180
    SystVariations[ 2] = "50200013_00200049007000008250400000"; //singleptcut = 75MeV
    SystVariations[ 3] = "50200013_00200019007000008250400000"; //singleptcut = 100MeV
    SystVariations[ 4] = "50200013_00200009007000009250400000"; // Qtmax = 0.03
    SystVariations[ 5] = "50200013_00200009007000002250400000"; // Qtmax = 0.06
    SystVariations[ 6] = "50200013_00200009007000008750400000"; //chi2cut = 10
    SystVariations[ 7] = "50200013_00200009007000008150400000"; //chi2cut = 50
    SystVariations[ 8] = "50200013_00200009007000008260400000"; //psipaircut = 0.05
    SystVariations[ 9] = "50200013_00200009007000008280400000"; //psipaircut = 0.2
    SystVariations[10] = "50200013_00200009007000008250600000"; //cos p angle cut = 0.9
    SystVariations[11] = "50200013_00200009007000008250300000"; //cos p angle cut = 0.75
  }else if(CentralityLow.CompareTo("20")==0){
    SystVariations[ 0] = "52400013_04200009007000008250400000"; //eta cut: |eta| <0.75
    SystVariations[ 1] = "52400013_00500009007000008250400000"; //minR= 10 maxR = 180
    SystVariations[ 2] = "52400013_00200049007000008250400000"; //singleptcut = 75MeV
    SystVariations[ 3] = "52400013_00200019007000008250400000"; //singleptcut = 100MeV
    SystVariations[ 4] = "52400013_00200009007000009250400000"; // Qtmax = 0.03
    SystVariations[ 5] = "52400013_00200009007000002250400000"; // Qtmax = 0.06
    SystVariations[ 6] = "52400013_00200009007000008750400000"; //chi2cut = 10
    SystVariations[ 7] = "52400013_00200009007000008150400000"; //chi2cut = 50
    SystVariations[ 8] = "52400013_00200009007000008260400000"; //psipaircut = 0.05
    SystVariations[ 9] = "52400013_00200009007000008280400000"; //psipaircut = 0.2
    SystVariations[10] = "52400013_00200009007000008250600000"; //cos p angle cut = 0.9
    SystVariations[11] = "52400013_00200009007000008250300000"; //cos p angle cut = 0.75
  }else if(CentralityLow.CompareTo("40")==0){
    SystVariations[ 0] = "54800013_04200009007000008250400000"; //eta cut: |eta| <0.75
    SystVariations[ 1] = "54800013_00500009007000008250400000"; //minR= 10 maxR = 180
    SystVariations[ 2] = "54800013_00200049007000008250400000"; //singleptcut = 75MeV
    SystVariations[ 3] = "54800013_00200019007000008250400000"; //singleptcut = 100MeV
    SystVariations[ 4] = "54800013_00200009007000009250400000"; // Qtmax = 0.03
    SystVariations[ 5] = "54800013_00200009007000002250400000"; // Qtmax = 0.06
    SystVariations[ 6] = "54800013_00200009007000008750400000"; //chi2cut = 10
    SystVariations[ 7] = "54800013_00200009007000008150400000"; //chi2cut = 50
    SystVariations[ 8] = "54800013_00200009007000008260400000"; //psipaircut = 0.05
    SystVariations[ 9] = "54800013_00200009007000008280400000"; //psipaircut = 0.2
    SystVariations[10] = "54800013_00200009007000008250600000"; //cos p angle cut = 0.9
    SystVariations[11] = "54800013_00200009007000008250300000"; //cos p angle cut = 0.75
  }
  
  TFile* fileDefaultInclusive  = new TFile(Form("/home/mike/3_PbPb_dirg/0_analysis/170216_v2_final_systematics/Results_%s/InclusivePhotonv2_Corrected_%i_%s.root",DefaultCut.Data(),Trainconfig,DefaultCut.Data()));
  TH1F*  histoDefaultInclusive = (TH1F*)fileDefaultInclusive->Get(Form("v2GammaIncl_corrected_%s%s_tC%i",CentralityLow.Data(),CentralityHigh.Data(),Trainconfig));
  if(!histoDefaultInclusive) cout << "histoDefaultInclusive not found in fileInclusive!!" << endl;
  histoDefaultInclusive->SetLineColor(kBlack);
  histoDefaultInclusive->SetMarkerColor(kBlack);
  histoDefaultInclusive->SetMarkerStyle(20);
  
  TH1F*  histoInclusive[12];
  TH1F*  histoInclusiveRatio[12];
  for(Int_t i=0;i<12;i++){
    if(i==4) SystTrainConfig++;
    if(i==8) SystTrainConfig++;
    TFile* fileInclusive  = new TFile(Form("/home/mike/3_PbPb_dirg/0_analysis/170216_v2_final_systematics/Results_%s/InclusivePhotonv2_Corrected_%i_%s.root",SystVariations[ i].Data(),SystTrainConfig,SystVariations[ i].Data()));
    cout << " fileInclusive " << fileInclusive << endl;
    histoInclusive[i] = (TH1F*)fileInclusive->Get(Form("v2GammaIncl_corrected_%s%s_tC%i",CentralityLow.Data(),CentralityHigh.Data(),SystTrainConfig));
    if(!histoInclusive[i]) cout << i << " histoInclusive not found in fileInclusive!!" << endl;
    if(i==0 || i==1 || i==2 || i==4 || i==6 || i==8 || i==10){
      histoInclusive[i]->SetLineColor(kRed+2);
      histoInclusive[i]->SetMarkerColor(kRed+2);
      histoInclusive[i]->SetMarkerStyle(20);
    }else{
     histoInclusive[i]->SetLineColor(kBlue);
      histoInclusive[i]->SetMarkerColor(kBlue);
      histoInclusive[i]->SetMarkerStyle(24); 
    }
    histoInclusiveRatio[i] = (TH1F*)histoInclusive[i]->Clone();
    histoInclusiveRatio[i]->Divide(histoInclusiveRatio[i],histoDefaultInclusive,1,1,"B");
  }
  
  //========================================================
  //histogram cosmetics
  //========================================================
  
  TH1F*  histoInclusiveEmpty = (TH1F*)Getv2HistStyle();
  histoInclusiveEmpty->GetYaxis()->SetRangeUser(0.0,0.49);
  TH1F*  histoInclusiveEmpty2 = (TH1F*)Getv2HistStyle();
  TH1F*  histoRatioInclusiveEmpty = (TH1F*)Getv2HistStyle();
  histoRatioInclusiveEmpty->GetYaxis()->SetRangeUser(0.5,1.5);
  histoRatioInclusiveEmpty->GetYaxis()->SetTitle("Ratio");
  
  TH1F*  histoSystEmpty = (TH1F*)GetSystHistStyle();
  
  TLatex T1;
  T1.SetTextSize(0.04);
  T1.SetTextAlign(12);
  T1.SetNDC();
  
  //========================================================
  //var1 //eta cut: |eta| <0.75
  //========================================================
  
  TLegend* legVar1 = new TLegend(0.16,0.15,0.6,0.26);
  legVar1->SetTextSize(0.05);
  legVar1->AddEntry(histoDefaultInclusive,"Default:   |#eta|<0.9","lp");
  legVar1->AddEntry(histoInclusive[0],    "Variation: |#eta|<0.75","lp");
  legVar1->SetBorderSize(0);
  
  TLegend* legVar1Ratio = new TLegend(0.16,0.15,0.6,0.26);
  legVar1Ratio->SetTextSize(0.05);
  legVar1Ratio->AddEntry(histoInclusiveRatio[0],"Ratio Variation1/Default","lp");
  legVar1Ratio->SetBorderSize(0);
  
  TCanvas* cVar1 = new TCanvas("cVar1","",400,800);
  cVar1->Divide(1,2);
  cVar1->cd(1);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoDefaultInclusive->Draw("SAME");
  histoInclusive[0]->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  legVar1->Draw();
  cVar1->cd(2);
  SetProperMargins();
  histoRatioInclusiveEmpty->Draw();
  histoInclusiveRatio[0]->Draw("SAME");
//   legVar1Ratio->Draw();
  
  //========================================================
  //var2 //minR= 10 maxR = 180
  //========================================================
  
  TLegend* legVar2 = new TLegend(0.16,0.15,0.6,0.26);
  legVar2->SetTextSize(0.05);
  legVar2->AddEntry(histoDefaultInclusive,"Default:   5<R<180","lp");
  legVar2->AddEntry(histoInclusive[1],    "Variation: 10<R<180","lp");
  legVar2->SetBorderSize(0);
  
  TLegend* legVar2Ratio = new TLegend(0.16,0.15,0.6,0.26);
  legVar2Ratio->SetTextSize(0.05);
  legVar2Ratio->AddEntry(histoInclusiveRatio[1],"Ratio Variation2/Default","lp");
  legVar2Ratio->SetBorderSize(0);
  
  TCanvas* cVar2 = new TCanvas("cVar2","",400,800);
  cVar2->Divide(1,2);
  cVar2->cd(1);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoDefaultInclusive->Draw("SAME");
  histoInclusive[1]->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  legVar2->Draw();
  cVar2->cd(2);
  SetProperMargins();
  histoRatioInclusiveEmpty->Draw();
  histoInclusiveRatio[1]->Draw("SAME");
//   legVar2Ratio->Draw();
  
  //========================================================
  //var3 //singleptcut = 75MeV //singleptcut = 100MeV
  //========================================================
  
  TLegend* legVar3 = new TLegend(0.16,0.65,0.6,0.89);
  legVar3->SetTextSize(0.05);
  legVar3->AddEntry(histoDefaultInclusive,"Default:   p_{T}>50MeV","lp");
  legVar3->AddEntry(histoInclusive[2],    "Variation: p_{T}>75MeV","lp");
  legVar3->AddEntry(histoInclusive[3],    "Variation: p_{T}>100MeV","lp");
  legVar3->SetBorderSize(0);
  
  TLegend* legVar3Ratio = new TLegend(0.16,0.15,0.6,0.26);
  legVar3Ratio->SetTextSize(0.05);
  legVar3Ratio->AddEntry(histoInclusiveRatio[2],"Ratio Variationa/Default","lp");
  legVar3Ratio->AddEntry(histoInclusiveRatio[3],"Ratio Variationb/Default","lp");
  legVar3Ratio->SetBorderSize(0);
  
  TCanvas* cVar3 = new TCanvas("cVar3","",400,800);
  cVar3->Divide(1,2);
  cVar3->cd(1);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoDefaultInclusive->Draw("SAME");
  histoInclusive[2]->Draw("SAME");
  histoInclusive[3]->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  legVar3->Draw();
  cVar3->cd(2);
  SetProperMargins();
  histoRatioInclusiveEmpty->Draw();
  histoInclusiveRatio[2]->Draw("SAME");
  histoInclusiveRatio[3]->Draw("SAME");
//   legVar3Ratio->Draw();
  
  //========================================================
  //var4 // Qtmax = 0.03 // Qtmax = 0.1
  //========================================================
  
  TLegend* legVar4 = new TLegend(0.16,0.65,0.6,0.89);
  legVar4->SetTextSize(0.05);
  legVar4->AddEntry(histoDefaultInclusive,"Default:   q_{t}max=0.05","lp");
  legVar4->AddEntry(histoInclusive[4],    "Variation: q_{t}max=0.03","lp");
  legVar4->AddEntry(histoInclusive[5],    "Variation: q_{t}max=0.10","lp");
  legVar4->SetBorderSize(0);
  
  TLegend* legVar4Ratio = new TLegend(0.16,0.15,0.6,0.26);
  legVar4Ratio->SetTextSize(0.05);
  legVar4Ratio->AddEntry(histoInclusiveRatio[4],"Ratio Variationa/Default","lp");
  legVar4Ratio->AddEntry(histoInclusiveRatio[5],"Ratio Variationb/Default","lp");
  legVar4Ratio->SetBorderSize(0);
  
  TCanvas* cVar4 = new TCanvas("cVar4","",400,800);
  cVar4->Divide(1,2);
  cVar4->cd(1);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoDefaultInclusive->Draw("SAME");
  histoInclusive[4]->Draw("SAME");
  histoInclusive[5]->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  legVar4->Draw();
  cVar4->cd(2);
  SetProperMargins();
  histoRatioInclusiveEmpty->Draw();
  histoInclusiveRatio[4]->Draw("SAME");
  histoInclusiveRatio[5]->Draw("SAME");
//   legVar4Ratio->Draw();
  
  //========================================================
  //var5 //chi2cut = 10 //chi2cut = 50
  //========================================================
  
  TLegend* legVar5 = new TLegend(0.16,0.65,0.6,0.89);
  legVar5->SetTextSize(0.05);
  legVar5->AddEntry(histoDefaultInclusive,"Default:   #chi^{2}<30","lp");
  legVar5->AddEntry(histoInclusive[6],    "Variation: #chi^{2}<10","lp");
  legVar5->AddEntry(histoInclusive[7],    "Variation: #chi^{2}<50","lp");
  legVar5->SetBorderSize(0);
  
  TLegend* legVar5Ratio = new TLegend(0.16,0.15,0.6,0.26);
  legVar5Ratio->SetTextSize(0.05);
  legVar5Ratio->AddEntry(histoInclusiveRatio[6],"Ratio Variationa/Default","lp");
  legVar5Ratio->AddEntry(histoInclusiveRatio[7],"Ratio Variationb/Default","lp");
  legVar5Ratio->SetBorderSize(0);
  
  TCanvas* cVar5 = new TCanvas("cVar5","",400,800);
  cVar5->Divide(1,2);
  cVar5->cd(1);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoDefaultInclusive->Draw("SAME");
  histoInclusive[6]->Draw("SAME");
  histoInclusive[7]->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  legVar5->Draw();
  cVar5->cd(2);
  SetProperMargins();
  histoRatioInclusiveEmpty->Draw();
  histoInclusiveRatio[6]->Draw("SAME");
  histoInclusiveRatio[7]->Draw("SAME");
//   legVar5Ratio->Draw();
  
  //========================================================
  //var6 //psipaircut = 0.05 //psipaircut = 0.2
  //========================================================
  
  TLegend* legVar6 = new TLegend(0.16,0.65,0.6,0.89);
  legVar6->SetTextSize(0.05);
  legVar6->AddEntry(histoDefaultInclusive,"Default:   #psi_{pair}<0.1","lp");
  legVar6->AddEntry(histoInclusive[8],    "Variation: #psi_{pair}<0.05","lp");
  legVar6->AddEntry(histoInclusive[9],    "Variation: #psi_{pair}<0.2","lp");
  legVar6->SetBorderSize(0);
  
  TLegend* legVar6Ratio = new TLegend(0.16,0.15,0.6,0.26);
  legVar6Ratio->SetTextSize(0.05);
  legVar6Ratio->AddEntry(histoInclusiveRatio[8],"Ratio Variationa/Default","lp");
  legVar6Ratio->AddEntry(histoInclusiveRatio[9],"Ratio Variationb/Default","lp");
  legVar6Ratio->SetBorderSize(0);
  
  TCanvas* cVar6 = new TCanvas("cVar6","",400,800);
  cVar6->Divide(1,2);
  cVar6->cd(1);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoDefaultInclusive->Draw("SAME");
  histoInclusive[8]->Draw("SAME");
  histoInclusive[9]->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  legVar6->Draw();
  cVar6->cd(2);
  SetProperMargins();
  histoRatioInclusiveEmpty->Draw();
  histoInclusiveRatio[8]->Draw("SAME");
  histoInclusiveRatio[9]->Draw("SAME");
//   legVar6Ratio->Draw();
  
  //========================================================
  //var7 //cos p angle cut = 0.9 //cos p angle cut = 0.75
  //========================================================
  
  TLegend* legVar7 = new TLegend(0.16,0.65,0.6,0.89);
  legVar7->SetTextSize(0.05);
  legVar7->AddEntry(histoDefaultInclusive,"Default:   cos(#theta)>0.95","lp");
  legVar7->AddEntry(histoInclusive[10],    "Variation: cos(#theta)>0.9","lp");
  legVar7->AddEntry(histoInclusive[11],    "Variation: cos(#theta)>0.75","lp");
  legVar7->SetBorderSize(0);
  
  TLegend* legVar7Ratio = new TLegend(0.16,0.15,0.6,0.26);
  legVar7Ratio->SetTextSize(0.05);
  legVar7Ratio->AddEntry(histoInclusiveRatio[10],"Ratio Variationa/Default","lp");
  legVar7Ratio->AddEntry(histoInclusiveRatio[11],"Ratio Variationb/Default","lp");
  legVar7Ratio->SetBorderSize(0);
  
  TCanvas* cVar7 = new TCanvas("cVar7","",400,800);
  cVar7->Divide(1,2);
  cVar7->cd(1);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  histoDefaultInclusive->Draw("SAME");
  histoInclusive[10]->Draw("SAME");
  histoInclusive[11]->Draw("SAME");
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  legVar7->Draw();
  cVar7->cd(2);
  SetProperMargins();
  histoRatioInclusiveEmpty->Draw();
  histoInclusiveRatio[10]->Draw("SAME");
  histoInclusiveRatio[11]->Draw("SAME");
//   legVar7Ratio->Draw();
  
  //========================================================
  //total
  //========================================================
  
  TGraphAsymmErrors* graphSystErrors1 = (TGraphAsymmErrors*)CalculateSymmetricSystematics(histoDefaultInclusive, histoInclusive[0]);
  TGraphAsymmErrors* graphSystErrors2 = (TGraphAsymmErrors*)CalculateSymmetricSystematics(histoDefaultInclusive, histoInclusive[1]);
  TGraphAsymmErrors* graphSystErrors3 = (TGraphAsymmErrors*)CalculateSymmetricSystematics(histoDefaultInclusive, histoInclusive[2], histoInclusive[3]);
  TGraphAsymmErrors* graphSystErrors4 = (TGraphAsymmErrors*)CalculateSymmetricSystematics(histoDefaultInclusive, histoInclusive[4], histoInclusive[5]);
  TGraphAsymmErrors* graphSystErrors5 = (TGraphAsymmErrors*)CalculateSymmetricSystematics(histoDefaultInclusive, histoInclusive[6], histoInclusive[7]);
  TGraphAsymmErrors* graphSystErrors6 = (TGraphAsymmErrors*)CalculateSymmetricSystematics(histoDefaultInclusive, histoInclusive[8], histoInclusive[9]);
  TGraphAsymmErrors* graphSystErrors7 = (TGraphAsymmErrors*)CalculateSymmetricSystematics(histoDefaultInclusive, histoInclusive[10], histoInclusive[11]);
  TGraphAsymmErrors* graphTotalSystErrors = (TGraphAsymmErrors*)CalculateTotalSystematics(graphSystErrors1,graphSystErrors2,graphSystErrors3,graphSystErrors4,graphSystErrors5,graphSystErrors6,graphSystErrors7);
  graphTotalSystErrors->SetFillColor(kGray);
  graphTotalSystErrors->SetFillStyle(0);
  
  TGraphAsymmErrors* graphTotalErrors = (TGraphAsymmErrors*)CalculateTotalError(histoDefaultInclusive,graphTotalSystErrors);
  graphTotalErrors->SetLineColor(kRed);
  graphTotalErrors->SetFillColor(kRed);
  graphTotalErrors->SetFillStyle(0);
  
  TGraphAsymmErrors* graphSystErrors1_rel = (TGraphAsymmErrors*)graphSystErrors1->Clone();
  TGraphAsymmErrors* graphSystErrors2_rel = (TGraphAsymmErrors*)graphSystErrors2->Clone();
  TGraphAsymmErrors* graphSystErrors3_rel = (TGraphAsymmErrors*)graphSystErrors3->Clone();
  TGraphAsymmErrors* graphSystErrors4_rel = (TGraphAsymmErrors*)graphSystErrors4->Clone();
  TGraphAsymmErrors* graphSystErrors5_rel = (TGraphAsymmErrors*)graphSystErrors5->Clone();
  TGraphAsymmErrors* graphSystErrors6_rel = (TGraphAsymmErrors*)graphSystErrors6->Clone();
  TGraphAsymmErrors* graphSystErrors7_rel = (TGraphAsymmErrors*)graphSystErrors7->Clone();
  graphSystErrors1_rel->SetLineColor(kMagenta+1);   graphSystErrors1_rel->SetMarkerStyle(20); graphSystErrors1_rel->SetMarkerColor(kMagenta+1); graphSystErrors1_rel->SetLineWidth(2);
  graphSystErrors2_rel->SetLineColor(kBlue+2);      graphSystErrors2_rel->SetMarkerStyle(21); graphSystErrors2_rel->SetMarkerColor(kBlue+2); graphSystErrors2_rel->SetLineWidth(2);
  graphSystErrors3_rel->SetLineColor(kBlue-4);      graphSystErrors3_rel->SetMarkerStyle(22); graphSystErrors3_rel->SetMarkerColor(kBlue-4); graphSystErrors3_rel->SetLineWidth(2);
  graphSystErrors4_rel->SetLineColor(kCyan-3);      graphSystErrors4_rel->SetMarkerStyle(23); graphSystErrors4_rel->SetMarkerColor(kCyan-3); graphSystErrors4_rel->SetLineWidth(2);
  graphSystErrors5_rel->SetLineColor(kGreen+3);     graphSystErrors5_rel->SetMarkerStyle(24); graphSystErrors5_rel->SetMarkerColor(kGreen+3); graphSystErrors5_rel->SetLineWidth(2);
  graphSystErrors6_rel->SetLineColor(kGreen-7);     graphSystErrors6_rel->SetMarkerStyle(25); graphSystErrors6_rel->SetMarkerColor(kGreen-7); graphSystErrors6_rel->SetLineWidth(2);
  graphSystErrors7_rel->SetLineColor(kYellow-3);    graphSystErrors7_rel->SetMarkerStyle(26); graphSystErrors7_rel->SetMarkerColor(kYellow-3); graphSystErrors7_rel->SetLineWidth(2);
  
  TGraphAsymmErrors* graphTotalErrors_rel = (TGraphAsymmErrors*)graphTotalSystErrors->Clone();
  graphTotalErrors_rel->SetMarkerStyle(28); graphTotalErrors_rel->SetLineColor(kRed); graphTotalErrors_rel->SetMarkerColor(kRed); graphTotalErrors_rel->SetLineWidth(2);
  
  Int_t NPoints = graphSystErrors1_rel->GetN();
  Double_t* g1_rel_x = graphSystErrors1_rel->GetX();
  Double_t yVAL;
  Double_t g1_rel_y[8];
  
  for(int ibin=0;ibin<NPoints;ibin++){
    if(g1_rel_x[ibin]>0.9){
      yVAL = histoDefaultInclusive->GetBinContent(histoDefaultInclusive->FindBin(g1_rel_x[ibin]));
      
      g1_rel_y[0] = TMath::Abs(100*graphSystErrors1->GetErrorYhigh(ibin)/yVAL); if(g1_rel_y[0]>40.0 || g1_rel_y[0] < 0.0) g1_rel_y[0] = 0.;
      g1_rel_y[1] = TMath::Abs(100*graphSystErrors2->GetErrorYhigh(ibin)/yVAL); if(g1_rel_y[1]>40.0 || g1_rel_y[1] < 0.0) g1_rel_y[1] = 0.;
      g1_rel_y[2] = TMath::Abs(100*graphSystErrors3->GetErrorYhigh(ibin)/yVAL); if(g1_rel_y[2]>40.0 || g1_rel_y[2] < 0.0) g1_rel_y[2] = 0.;
      g1_rel_y[3] = TMath::Abs(100*graphSystErrors4->GetErrorYhigh(ibin)/yVAL); if(g1_rel_y[3]>40.0 || g1_rel_y[3] < 0.0) g1_rel_y[3] = 0.;
      g1_rel_y[4] = TMath::Abs(100*graphSystErrors5->GetErrorYhigh(ibin)/yVAL); if(g1_rel_y[4]>40.0 || g1_rel_y[4] < 0.0) g1_rel_y[4] = 0.;
      g1_rel_y[5] = TMath::Abs(100*graphSystErrors6->GetErrorYhigh(ibin)/yVAL); if(g1_rel_y[5]>40.0 || g1_rel_y[5] < 0.0) g1_rel_y[5] = 0.;
      g1_rel_y[6] = TMath::Abs(100*graphSystErrors7->GetErrorYhigh(ibin)/yVAL); if(g1_rel_y[6]>40.0 || g1_rel_y[6] < 0.0) g1_rel_y[6] = 0.;
      g1_rel_y[7] = TMath::Abs(100*graphTotalSystErrors->GetErrorYhigh(ibin)/yVAL); if(g1_rel_y[7]>40.0 || g1_rel_y[7] < 0.0) g1_rel_y[7] = 0.;
      
      
      graphSystErrors1_rel->SetPoint(ibin,g1_rel_x[ibin],g1_rel_y[0]); graphSystErrors1_rel->SetPointEYlow(ibin,0.0); graphSystErrors1_rel->SetPointEYhigh(ibin,0.0);
      graphSystErrors2_rel->SetPoint(ibin,g1_rel_x[ibin],g1_rel_y[1]); graphSystErrors2_rel->SetPointEYlow(ibin,0.0); graphSystErrors2_rel->SetPointEYhigh(ibin,0.0);
      graphSystErrors3_rel->SetPoint(ibin,g1_rel_x[ibin],g1_rel_y[2]); graphSystErrors3_rel->SetPointEYlow(ibin,0.0); graphSystErrors3_rel->SetPointEYhigh(ibin,0.0);
      graphSystErrors4_rel->SetPoint(ibin,g1_rel_x[ibin],g1_rel_y[3]); graphSystErrors4_rel->SetPointEYlow(ibin,0.0); graphSystErrors4_rel->SetPointEYhigh(ibin,0.0);
      graphSystErrors5_rel->SetPoint(ibin,g1_rel_x[ibin],g1_rel_y[4]); graphSystErrors5_rel->SetPointEYlow(ibin,0.0); graphSystErrors5_rel->SetPointEYhigh(ibin,0.0);
      graphSystErrors6_rel->SetPoint(ibin,g1_rel_x[ibin],g1_rel_y[5]); graphSystErrors6_rel->SetPointEYlow(ibin,0.0); graphSystErrors6_rel->SetPointEYhigh(ibin,0.0);
      graphSystErrors7_rel->SetPoint(ibin,g1_rel_x[ibin],g1_rel_y[6]); graphSystErrors7_rel->SetPointEYlow(ibin,0.0); graphSystErrors7_rel->SetPointEYhigh(ibin,0.0);
      graphTotalErrors_rel->SetPoint(ibin,g1_rel_x[ibin],g1_rel_y[7]); graphTotalErrors_rel->SetPointEYlow(ibin,0.0); graphTotalErrors_rel->SetPointEYhigh(ibin,0.0);
    }
  }
  
  if(SetSystSmoothing){
    //SmoothSystematicError(graphTotalErrors_rel,"");
    SmoothSystematicError(graphSystErrors1_rel,SystNames[0].Data());
    SmoothSystematicError(graphSystErrors2_rel,SystNames[1].Data());
    SmoothSystematicError(graphSystErrors3_rel,SystNames[2].Data());
    SmoothSystematicError(graphSystErrors4_rel,SystNames[3].Data());
    SmoothSystematicError(graphSystErrors5_rel,SystNames[4].Data());
    SmoothSystematicError(graphSystErrors6_rel,SystNames[5].Data());
    SmoothSystematicError(graphSystErrors7_rel,SystNames[6].Data());
  }
  
  
  //========================================================
  //Plotting and Saving
  //========================================================
  
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  TLegend* leg1 = new TLegend(0.16,0.15,0.6,0.29);
  leg1->SetTextSize(0.05);
  leg1->AddEntry(histoDefaultInclusive,"#sigma_{stat}","lp");
  leg1->AddEntry(graphTotalSystErrors,"#sigma_{syst}","f");
  leg1->AddEntry(graphTotalErrors,"#sigma_{total}","f");
  leg1->SetBorderSize(0);
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  SetProperMargins();
  histoInclusiveEmpty2->Draw();
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  T1.DrawLatex(0.17, 0.94, "This thesis");
  leg1->Draw();
  
  graphTotalSystErrors->Draw("same e2 P");
  graphTotalErrors->Draw("same e2 P");
  histoDefaultInclusive->Draw("SAME");
  
  TLegend* leg2 = new TLegend(0.16,0.39,0.4,0.94);
  leg2->SetTextSize(0.05);
  leg2->AddEntry(graphTotalErrors_rel,"total","lp");
  leg2->AddEntry(graphSystErrors1_rel,"|#eta|","lp");
  leg2->AddEntry(graphSystErrors2_rel,"R","lp");
  leg2->AddEntry(graphSystErrors3_rel,"min p_{T}","lp");
  leg2->AddEntry(graphSystErrors4_rel,"q_{T}","lp");
  leg2->AddEntry(graphSystErrors5_rel,"#chi^{2}","lp");
  leg2->AddEntry(graphSystErrors6_rel,"#psi_{pair}","lp");
  leg2->AddEntry(graphSystErrors7_rel,"cos(#theta)","lp");
  leg2->SetBorderSize(0);
  
  TCanvas* c2 = new TCanvas("c2","",800,800);
  SetProperMargins();
  histoSystEmpty->Draw();
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  graphSystErrors1_rel->Draw("same P");
  graphSystErrors2_rel->Draw("same P");
  graphSystErrors3_rel->Draw("same P");
  graphSystErrors4_rel->Draw("same P");
  graphSystErrors5_rel->Draw("same P");
  graphSystErrors6_rel->Draw("same P");
  graphSystErrors7_rel->Draw("same P");
  graphTotalErrors_rel->Draw("same P");
  leg2->Draw();
  
  
  gSystem->mkdir("Systematics");
  c1->SaveAs(Form("Systematics/%s%s_v2GammaInclErrors.eps", CentralityLow.Data(), CentralityHigh.Data()));
  c2->SaveAs(Form("Systematics/%s%s_SystematicContributions.eps", CentralityLow.Data(), CentralityHigh.Data()));
  cVar1->SaveAs(Form("Systematics/%s%s_Var1.eps", CentralityLow.Data(), CentralityHigh.Data()));
  cVar2->SaveAs(Form("Systematics/%s%s_Var2.eps", CentralityLow.Data(), CentralityHigh.Data()));
  cVar3->SaveAs(Form("Systematics/%s%s_Var3.eps", CentralityLow.Data(), CentralityHigh.Data()));
  cVar4->SaveAs(Form("Systematics/%s%s_Var4.eps", CentralityLow.Data(), CentralityHigh.Data()));
  cVar5->SaveAs(Form("Systematics/%s%s_Var5.eps", CentralityLow.Data(), CentralityHigh.Data()));
  cVar6->SaveAs(Form("Systematics/%s%s_Var6.eps", CentralityLow.Data(), CentralityHigh.Data()));
  cVar7->SaveAs(Form("Systematics/%s%s_Var7.eps", CentralityLow.Data(), CentralityHigh.Data()));
  
  TFile *output_File = new TFile(Form("Systematics/PCM_InclusivePhotonFlow_Syst_%s%s.root", CentralityLow.Data(), CentralityHigh.Data()),"RECREATE");
  histoDefaultInclusive->Write(Form("histoInclusive_%s%s", CentralityLow.Data(), CentralityHigh.Data()));
  graphTotalSystErrors->Write(Form("graphTotalSystErrors_%s%s", CentralityLow.Data(), CentralityHigh.Data()));
  graphTotalErrors->Write(Form("graphTotalErrors_%s%s", CentralityLow.Data(), CentralityHigh.Data()));
  output_File->Close();
  
}

TGraphAsymmErrors* CalculateSymmetricSystematics(TH1F* histDefault, TH1F* histvar){
  
  const Int_t NBins = histDefault->GetSize();
  Float_t x[NBins], y[NBins], exl[NBins], exh[NBins], eyl[NBins], eyh[NBins];
  Float_t ErrorUp, ErrorDown;
  Float_t v2InclDefaultValue, v2InclVarValue, BinCenter;
  for(Int_t i=0;i<NBins;i++){    
    v2InclDefaultValue = histDefault->GetBinContent(i);
    v2InclVarValue = histvar->GetBinContent(i);
    ErrorUp = TMath::Abs(v2InclDefaultValue-v2InclVarValue);
    ErrorDown = ErrorUp;
    x[i] = histDefault->GetBinCenter(i);
    y[i] = v2InclDefaultValue;
    exl[i] = histDefault->GetBinWidth(i)/2.0;
    exh[i] = histDefault->GetBinWidth(i)/2.0;
    eyl[i] = ErrorDown;
    eyh[i] = ErrorUp;
  }
  
  TGraphAsymmErrors* outputGraph = new TGraphAsymmErrors(NBins,x,y,exl,exh,eyl,eyh);
  return outputGraph;
  
}
TGraphAsymmErrors* CalculateSymmetricSystematics(TH1F* histDefault, TH1F* histvar1, TH1F* histvar2){
  
  const Int_t NBins = histDefault->GetSize();
  Float_t x[NBins], y[NBins], exl[NBins], exh[NBins], eyl[NBins], eyh[NBins];
  Float_t ErrorUp, ErrorDown;
  Float_t v2InclDefaultValue, v2InclVarValue1, v2InclVarValue2, BinCenter;
  for(Int_t i=0;i<NBins;i++){    
    v2InclDefaultValue = histDefault->GetBinContent(i);
    v2InclVarValue1 = histvar1->GetBinContent(i);
    v2InclVarValue2 = histvar2->GetBinContent(i);
    ErrorUp = TMath::Sqrt(((v2InclDefaultValue-v2InclVarValue1)*(v2InclDefaultValue-v2InclVarValue1)+(v2InclDefaultValue-v2InclVarValue2)*(v2InclDefaultValue-v2InclVarValue2))/2);
    ErrorDown = ErrorUp;
    x[i] = histDefault->GetBinCenter(i);
    y[i] = v2InclDefaultValue;
    exl[i] = histDefault->GetBinWidth(i)/2.0;
    exh[i] = histDefault->GetBinWidth(i)/2.0;
    eyl[i] = ErrorDown;
    eyh[i] = ErrorUp;
  }
  
  TGraphAsymmErrors* outputGraph = new TGraphAsymmErrors(NBins,x,y,exl,exh,eyl,eyh);
  return outputGraph;
  
}
TGraphAsymmErrors* CalculateTotalSystematics(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, TGraphAsymmErrors* g3, TGraphAsymmErrors* g4, TGraphAsymmErrors* g5, TGraphAsymmErrors* g6, TGraphAsymmErrors* g7){
  
  Int_t NPoints = g1->GetN();
  
  TGraphAsymmErrors *gCombined = (TGraphAsymmErrors*)g1->Clone();
  Double_t ErrorLow, ErrorHigh;
  Double_t e1,e2,e3,e4,e5,e6,e7;
  
  for(int ibin=0;ibin<NPoints;ibin++){
    e1 = g1->GetErrorYlow(ibin); if(e1>0.05){e1 = 0; cout << "error e1 " << e1 << " set to zero!" << endl;}
    e2 = g2->GetErrorYlow(ibin); if(e2>0.05){e2 = 0; cout << "error e2 " << e2 << " set to zero!" << endl;}
    e3 = g3->GetErrorYlow(ibin); if(e3>0.05){e3 = 0; cout << "error e3 " << e3 << " set to zero!" << endl;}
    e4 = g4->GetErrorYlow(ibin); if(e4>0.05){e4 = 0; cout << "error e4 " << e4 << " set to zero!" << endl;}
    e5 = g5->GetErrorYlow(ibin); if(e5>0.05){e5 = 0; cout << "error e5 " << e5 << " set to zero!" << endl;}
    e6 = g6->GetErrorYlow(ibin); if(e6>0.05){e6 = 0; cout << "error e6 " << e6 << " set to zero!" << endl;}
    e7 = g7->GetErrorYlow(ibin); if(e7>0.05){e7 = 0; cout << "error e7 " << e7 << " set to zero!" << endl;}
    ErrorLow =  TMath::Sqrt(pow(e1,2)+pow(e2,2)+pow(e3,2)+
                pow(e4,2)+pow(e5,2)+pow(e6,2)+pow(e7,2));
    ErrorHigh=  ErrorLow;
    gCombined->SetPointError(ibin,g1->GetErrorXlow(ibin),g1->GetErrorXhigh(ibin),ErrorLow,ErrorHigh);
  }
  
  return gCombined;
  
}
TGraphAsymmErrors* CalculateTotalError(TH1F* h1, TGraphAsymmErrors* g1){
  
  Int_t NPoints = g1->GetN();
  Double_t* g1_x = g1->GetX();
  
  TGraphAsymmErrors *gCombined = (TGraphAsymmErrors*)g1->Clone();
  Double_t ErrorLow, ErrorHigh;
  
  for(int ibin=0;ibin<NPoints;ibin++){
    ErrorLow =  TMath::Sqrt(pow(g1->GetErrorYlow(ibin),2)+pow(h1->GetBinError(h1->FindBin(g1_x[ibin])),2));
    ErrorHigh=  TMath::Sqrt(pow(g1->GetErrorYhigh(ibin),2)+pow(h1->GetBinError(h1->FindBin(g1_x[ibin])),2));
    gCombined->SetPointError(ibin,g1->GetErrorXlow(ibin),g1->GetErrorXhigh(ibin),ErrorLow,ErrorHigh);
  }
  
  return gCombined;
  
}
void SmoothSystematicError(TGraphAsymmErrors* g1,TString NameSyst){
  
  //TString SystNames[7] = {"eta","R","minPt","qt","chi2","psipair","cosp"};
  
  Int_t NPoints = g1->GetN();
  Double_t* g1_x = g1->GetX();
  
  TF1* fitSystematicPol0;
  TF1* fitSystematicPol1;
  
  fitSystematicPol0 = new TF1("fitSystematicPol0","pol0",0.9,3.0);
  fitSystematicPol0->SetLineColor(g1->GetLineColor());
  
  fitSystematicPol1 = new TF1("fitSystematicPol1","pol1",0.9,5.8);
  fitSystematicPol1->SetLineColor(g1->GetLineColor());
  
  
  cout << "=============" << endl;
  cout << "fitting: " << NameSyst.Data() << endl;
  cout << "=============" << endl;
  g1->Fit("fitSystematicPol0","WR");
  g1->Fit("fitSystematicPol1","WR");
  
  Double_t valPol0 = fitSystematicPol0->GetParameter(0);
  Double_t valPol1;
  
  for(int ibin=0;ibin<NPoints;ibin++){
    valPol1 = 
    g1->SetPoint(ibin,g1_x[ibin],);
  }

  
}