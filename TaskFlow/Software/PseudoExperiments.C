#include "DirectPhotonFlowFunctions.h"

void  PseudoExperiments(
                                      TString CentralityLow = "0",
                                      TString CentralityHigh = "20",
                                      Bool_t IncludeTheory = kTRUE
                                 ){
  
  //========================================================
  //Opening files and creating histograms
  //========================================================
  
  TFile* fileInclusive  = new TFile(Form("/home/mike/git_afterburner/AnalysisSoftware/TaskFlow/Results/PCM_InclusivePhotonFlow_Syst_%s%s.root",CentralityLow.Data(),CentralityHigh.Data()));
  TH1F*  histoInclusive = (TH1F*)fileInclusive->Get(Form("histoInclusive_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  if(!histoInclusive) cout << "histoInclusive not found in fileInclusive!!" << endl;
  
  TGraphAsymmErrors* graphInclusive = (TGraphAsymmErrors*)fileInclusive->Get(Form("graphTotalErrors_%s%s",CentralityLow.Data(),CentralityHigh.Data()));
  if(!graphInclusive) cout << "graphInclusive not found in fileInclusive!!" << endl;
  graphInclusive->SetMarkerStyle(20);
  graphInclusive->SetMarkerColor(kBlack);
  graphInclusive->SetLineColor(kBlack);
  graphInclusive->SetLineWidth(3);

  TFile* fileRGamma = new TFile("/home/mike/git_afterburner/AnalysisSoftware/TaskFlow/Results/Gamma_CombResults_PbPb_2.76TeV.root");
  TDirectory* dir1 = new TDirectory();
  dir1 = (TDirectory*)fileRGamma->Get(Form("Gamma_PbPb_2.76TeV_%s-%s%%",CentralityLow.Data(),CentralityHigh.Data()));
  TGraphAsymmErrors*  graphRGamma = (TGraphAsymmErrors*)dir1->Get("DR_comb_totErr");
  if(!graphRGamma) cout << "graphRGamma not found in fileRGamma!!" << endl;
  graphRGamma->SetMarkerStyle(21);
  graphRGamma->SetMarkerColor(kBlack);
  graphRGamma->SetLineColor(kBlack);
  graphRGamma->SetLineWidth(3);
  
  TFile* fileCocktail = new TFile("/home/mike/git_afterburner/AnalysisSoftware/TaskFlow/Results/CocktailV2_2017.root");
  //cen3 = 2040
  //cen6 = 020
  //cen8 = 4080
  TString CocktailString;
  if(CentralityLow.CompareTo("0")==0) CocktailString = "cen6";
  if(CentralityLow.CompareTo("20")==0) CocktailString = "cen3";
  if(CentralityLow.CompareTo("40")==0) CocktailString = "cen8";
  TH1F*  histoCocktail = (TH1F*)fileCocktail->Get(Form("v2gammaCocktail_%s",CocktailString.Data()));
  if(!histoCocktail) cout << "histoCocktail not found in fileCocktail!!" << endl;
  histoCocktail->SetMarkerStyle(25);
  histoCocktail->SetMarkerColor(kGreen+2);
  histoCocktail->SetLineColor(kGreen+2);
  histoCocktail->SetLineWidth(3);
  MaskPoints(histoCocktail,0,2);

  gSystem->mkdir(Form("PseudoExperiments_%s%s", CentralityLow.Data(), CentralityHigh.Data()));
  //========================================================
  //histogram cosmetics
  //========================================================
  
  TH1F*  histoDirectEmpty = (TH1F*)Getv2HistStyle();
  if(CentralityLow.CompareTo("0")==0) histoDirectEmpty->GetYaxis()->SetRangeUser(-0.49,0.49);
  if(CentralityLow.CompareTo("20")==0) histoDirectEmpty->GetYaxis()->SetRangeUser(-0.49,0.49);
  if(CentralityLow.CompareTo("40")==0) histoDirectEmpty->GetYaxis()->SetRangeUser(-1.09,1.09);
  
  TH1F*  histoRGammaEmpty = (TH1F*)GetRGammaHistStyle();
  
  TH1F*  histoInclusiveEmpty = (TH1F*)Getv2HistStyle();
  
  //========================================================
  //Performing pseudo experiment
  //========================================================
  
  Double_t startPt = 1.0;
  Double_t endPt = 5.8;
  
  Int_t startPtBin = histoInclusive->FindBin(startPt);
  Int_t endPtBin = histoInclusive->FindBin(endPt);
  const int n = endPtBin - startPtBin + 1;
  Double_t x[n],y[n],exl[n],exh[n],eyl[n],eyh[n];
  Double_t systvar1l[n];
  Double_t systvar1h[n];
  
  TRandom3* r = new TRandom3();
  r->SetSeed(0);
  
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  TLatex T1;
  T1.SetTextSize(0.04);
  T1.SetTextAlign(12);
  T1.SetNDC();
  
  TCanvas* c1 = new TCanvas("c1","",400,400);
  SetProperMargins();
  
  TCanvas* c3 = new TCanvas("c3","",1200,400);
  c3->Divide(3,1);
//   SetProperMargins();
  
  Double_t* xRGAMMA = graphRGamma->GetX();
  
  for(Int_t i=0;i<n;i++){
    
    c1->cd();
    
    Double_t BinCenter = histoInclusive->GetBinCenter(startPtBin+i);
    Double_t v2InclValue = histoInclusive->GetBinContent(histoInclusive->FindBin(BinCenter));
    Double_t v2CocktailValue = histoCocktail->GetBinContent(histoCocktail->FindBin(BinCenter));
    Double_t RGammaValue = graphRGamma->Eval(BinCenter);
    
    Double_t v2InclSigma = graphInclusive->GetErrorYhigh(startPtBin+i);
    Double_t v2CocktailSigma = histoCocktail->GetBinError(histoCocktail->FindBin(BinCenter));
    Double_t RGammaSigma = graphRGamma->GetErrorYhigh(i);
    
    Double_t v2InclRand;
    Double_t v2CocktailRand;
    Double_t RGammaRand;
    
    Double_t v2DirectCentralValue = ( v2InclValue * RGammaValue - v2CocktailValue ) / (RGammaValue - 1);
    Double_t v2DirectRandValue;
    
    Float_t RGammaFactor = (v2CocktailValue-v2InclValue)/((RGammaValue-1)*(RGammaValue-1));
    Float_t InclFactor = RGammaValue/(RGammaValue-1);
    Float_t CocktailFactor = 1/(RGammaValue-1);
      
    Float_t v2ErrorApprox = TMath::Sqrt(RGammaFactor*RGammaFactor*RGammaSigma*RGammaSigma+InclFactor*InclFactor*v2InclSigma*v2InclSigma+CocktailFactor*CocktailFactor*v2CocktailSigma*v2CocktailSigma);
    
    
    TH1D* checkhist1 = new TH1D("","",1000,v2InclValue-10*v2InclSigma,v2InclValue+10*v2InclSigma);
    TH1D* checkhist2 = new TH1D("","",1000,v2CocktailValue-10*v2CocktailSigma,v2CocktailValue+10*v2CocktailSigma);
    TH1D* checkhist3 = new TH1D("","",1000,RGammaValue-10*RGammaSigma,RGammaValue+10*RGammaSigma);
    TH1D* checkhist4 = new TH1D("","",1000,v2DirectCentralValue-5*v2ErrorApprox,v2DirectCentralValue+5*v2ErrorApprox);
    checkhist1->GetXaxis()->SetTitle("v_{2}^{#gamma,inc}");
    checkhist2->GetXaxis()->SetTitle("v_{2}^{#gamma,dec}");
    checkhist3->GetXaxis()->SetTitle("R_{#gamma}");
    checkhist4->GetXaxis()->SetTitle("v_{2}^{#gamma,direct}");
    
    checkhist1->GetXaxis()->SetTitleSize(0.05);
    checkhist1->GetXaxis()->SetTitleOffset(0.8);
    
    checkhist2->GetXaxis()->SetTitleSize(0.05);
    checkhist2->GetXaxis()->SetTitleOffset(0.8);
    
    checkhist3->GetXaxis()->SetTitleSize(0.05);
    checkhist3->GetXaxis()->SetTitleOffset(0.8);
    
    checkhist4->GetXaxis()->SetTitleSize(0.05);
    checkhist4->GetXaxis()->SetTitleOffset(0.8);
    
    for(Int_t j=0;j<1000000;j++){
      v2InclRand = v2InclValue + r->Gaus(0,v2InclSigma);
      v2CocktailRand = v2CocktailValue + r->Gaus(0,v2CocktailSigma);
      RGammaRand = RGammaValue + r->Gaus(0,RGammaSigma);
//       if(RGammaRand<1 || RGammaRand == 1) RGammaRand = 1.001;
      if(RGammaRand<1 || RGammaRand == 1) continue;
      
      v2DirectRandValue = ( v2InclRand * RGammaRand - v2CocktailRand ) / (RGammaRand - 1);
      
      checkhist1->Fill(v2InclRand);
      checkhist2->Fill(v2CocktailRand);
      checkhist3->Fill(RGammaRand);
      checkhist4->Fill(v2DirectRandValue);
    }
    checkhist1->Scale(1/checkhist1->GetEntries());
    checkhist2->Scale(1/checkhist2->GetEntries());
    checkhist3->Scale(1/checkhist3->GetEntries());
    checkhist4->Scale(1/checkhist4->GetEntries());
    
    TF1 *myfit = new TF1("myfit","[0]*exp(-0.5*((x-[1])/[2])^2)", v2DirectCentralValue-10*v2ErrorApprox,v2DirectCentralValue+10*v2ErrorApprox);
    myfit->SetParameter(0,checkhist4->GetMaximum());
    myfit->SetParameter(1,v2DirectCentralValue);
    myfit->SetParameter(2,v2ErrorApprox);
    
    Double_t HM = checkhist4->GetMaximum() / 2.0;
    
    Double_t ErrorLow = 0, ErrorHigh = 0, Median = 0, Integral = 0;
    
    for(Int_t j=checkhist4->GetXaxis()->FindBin(v2DirectCentralValue-10*v2ErrorApprox);j<checkhist4->GetXaxis()->FindBin(v2DirectCentralValue+10*v2ErrorApprox);j++){
      Integral = checkhist4->Integral(checkhist4->GetXaxis()->FindBin(v2DirectCentralValue-10*v2ErrorApprox),j);
      if(Integral>0.5){
        Median = checkhist4->GetBinCenter(j);
        break;
      }
    }
    
    for(Int_t j=checkhist4->GetXaxis()->FindBin(v2DirectCentralValue-10*v2ErrorApprox);j<checkhist4->GetXaxis()->FindBin(v2DirectCentralValue+10*v2ErrorApprox);j++){
      Integral = checkhist4->Integral(checkhist4->GetXaxis()->FindBin(v2DirectCentralValue-10*v2ErrorApprox),j);
      if(Integral>0.15865){
        ErrorLow = Median-checkhist4->GetBinCenter(j);
        break;
      }
    }
    
    for(Int_t j=checkhist4->GetXaxis()->FindBin(v2DirectCentralValue-10*v2ErrorApprox);j<checkhist4->GetXaxis()->FindBin(v2DirectCentralValue+10*v2ErrorApprox);j++){
      Integral = checkhist4->Integral(checkhist4->GetXaxis()->FindBin(v2DirectCentralValue-10*v2ErrorApprox),j);
      if(Integral>0.84135){
        ErrorHigh = checkhist4->GetBinCenter(j)-Median;
        break;
      }
    }
    
    TLine *lineCentralValue = new TLine(v2DirectCentralValue,0,v2DirectCentralValue,0.75*checkhist4->GetMaximum());
    TLine *lineErrorApproxLow = new TLine(v2DirectCentralValue-v2ErrorApprox,0,v2DirectCentralValue-v2ErrorApprox,0.75*checkhist4->GetMaximum());
    TLine *lineErrorApproxHigh = new TLine(v2DirectCentralValue+v2ErrorApprox,0,v2DirectCentralValue+v2ErrorApprox,0.75*checkhist4->GetMaximum());
    TLine *lineMeanValue = new TLine(checkhist4->GetMean(),0,checkhist4->GetMean(),0.9*checkhist4->GetMaximum());
    TLine *lineHM = new TLine(v2DirectCentralValue-5*v2ErrorApprox,HM,v2DirectCentralValue+5*v2ErrorApprox,HM);
    
    TLine *lineIntegralMedian = new TLine(Median,0,Median,0.5*checkhist4->GetMaximum());
    TLine *lineIntegralErrorLow = new TLine(v2DirectCentralValue-ErrorLow,0,v2DirectCentralValue-ErrorLow,0.5*checkhist4->GetMaximum());
    TLine *lineIntegralErrorHigh = new TLine(v2DirectCentralValue+ErrorHigh,0,v2DirectCentralValue+ErrorHigh,0.5*checkhist4->GetMaximum());
    
    TLine *linev2InclValue = new TLine(v2InclValue,0,v2InclValue,0.75*checkhist1->GetMaximum());
    TLine *linev2CocktailValue = new TLine(v2CocktailValue,0,v2CocktailValue,0.75*checkhist2->GetMaximum());
    TLine *lineRGammaValue = new TLine(RGammaValue,0,RGammaValue,0.75*checkhist3->GetMaximum());
    
    lineCentralValue->SetLineWidth(2);
    lineErrorApproxLow->SetLineWidth(2);
    lineErrorApproxHigh->SetLineWidth(2);
    lineMeanValue->SetLineWidth(3);
    lineHM->SetLineWidth(2);
    
    lineCentralValue->SetLineStyle(2);
    lineErrorApproxLow->SetLineStyle(2);
    lineErrorApproxHigh->SetLineStyle(2);
    lineHM->SetLineStyle(2);
    
    lineIntegralMedian->SetLineWidth(3);
    lineIntegralErrorLow->SetLineWidth(3);
    lineIntegralErrorHigh->SetLineWidth(3);
    
    lineIntegralMedian->SetLineColor(kRed);
    lineIntegralErrorLow->SetLineColor(kRed);
    lineIntegralErrorHigh->SetLineColor(kRed);
    
    lineMeanValue->SetLineColor(kBlue);
    
    TLegend* leg = new TLegend(0.15,0.75,0.6,0.86);
    leg->SetTextSize(0.04);
    leg->AddEntry(lineCentralValue,"Gaussian error propagation","l");
    leg->AddEntry(lineIntegralMedian,"-1#sigma , Median, +1#sigma","l");
    leg->AddEntry(lineMeanValue,"Mean","l");
    leg->SetBorderSize(0);
    
    checkhist1->GetYaxis()->SetRangeUser(0,checkhist1->GetMaximum()*1.5);
    checkhist2->GetYaxis()->SetRangeUser(0,checkhist2->GetMaximum()*1.5);
    checkhist3->GetYaxis()->SetRangeUser(0,checkhist3->GetMaximum()*1.5);
    checkhist4->GetYaxis()->SetRangeUser(0,checkhist4->GetMaximum()*1.5);
    checkhist4->Draw();
    lineCentralValue->Draw();
    lineErrorApproxLow->Draw();
    lineErrorApproxHigh->Draw();
    lineIntegralMedian->Draw();
    lineIntegralErrorLow->Draw();
    lineIntegralErrorHigh->Draw();
    lineMeanValue->Draw();
//     lineHM->Draw();
    T1.DrawLatex(0.15, 0.9, Form("p_{T} = %2.1f GeV/c",BinCenter));
    DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    leg->Draw();
    c1->SaveAs(Form("PseudoExperiments_%s%s/v2_Gaus_Errors_%s%s_%i.eps", CentralityLow.Data(), CentralityHigh.Data(), CentralityLow.Data(), CentralityHigh.Data(),i));
    
    x[i] = BinCenter;
    y[i] = Median;
    exl[i] = histoInclusive->GetBinWidth(startPtBin+i)/2.0;
    exh[i] = histoInclusive->GetBinWidth(startPtBin+i)/2.0;
    eyl[i] = ErrorLow;
    eyh[i] = ErrorHigh;
    
    c3->cd(1);
    SetProperMargins();
    checkhist1->Draw();
    T1.DrawLatex(0.15, 0.9, Form("p_{T} = %2.1f GeV/c",BinCenter));
    T1.DrawLatex(0.15, 0.8, Form("v_{2}^{#gamma,inc}=%2.3f",v2InclValue));
    T1.DrawLatex(0.15, 0.7, Form("#sigma=%2.3f",v2InclSigma));
    DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    linev2InclValue->Draw();
    c3->cd(2);
    SetProperMargins();
    checkhist2->Draw();
    T1.DrawLatex(0.15, 0.9, Form("p_{T} = %2.1f GeV/c",BinCenter));
    T1.DrawLatex(0.15, 0.8, Form("v_{2}^{#gamma,dec}=%2.3f",v2CocktailValue));
    T1.DrawLatex(0.15, 0.7, Form("#sigma=%2.3f",v2CocktailSigma));
    DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    linev2CocktailValue->Draw();
    c3->cd(3);
    SetProperMargins();
    checkhist3->Draw();
    T1.DrawLatex(0.15, 0.9, Form("p_{T} = %2.1f GeV/c",BinCenter));
    T1.DrawLatex(0.15, 0.8, Form("R_{#gamma}=%2.3f",RGammaValue));
    T1.DrawLatex(0.15, 0.7, Form("#sigma=%2.3f",RGammaSigma));
    DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
    lineRGammaValue->Draw();
    
    c3->SaveAs(Form("PseudoExperiments_%s%s/v2_Gaus_Errors_input_%s%s_%i.eps", CentralityLow.Data(), CentralityHigh.Data(), CentralityLow.Data(), CentralityHigh.Data(),i));
    
  
  }
  
  //========================================================
  //Creating the graph with total error
  //========================================================  
  TGraphAsymmErrors* graphDirect1 = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
  graphDirect1->SetMarkerStyle(21);
  graphDirect1->SetMarkerSize(1.2);
  graphDirect1->SetMarkerColor(kRed+2);
  graphDirect1->SetLineColor(kRed+2);
  graphDirect1->SetFillColor(kGray);
  graphDirect1->SetFillStyle(0);
  
  
  //========================================================
  //Opening files for theory comparison
  //========================================================
  if(IncludeTheory && (CentralityLow.CompareTo("0")==0 || CentralityLow.CompareTo("20")==0) ){

    ifstream inputfile;
//     inputfile.open("/home/mike/3_PbPb_dirg/0_analysis/161025_v2_theory/direct_photons_cent0020_Paquet2016.dat");
    inputfile.open("/home/mike/3_PbPb_dirg/0_analysis/161025_v2_theory/direct_photons_cent2040_Paquet2016.dat");
//     if(CentralityLow.CompareTo("0")==0) inputfile.open("/home/mike/git_afterburner/AnalysisSoftware/TaskFlow/Theory/direct_photons_cent0020_Paquet2016.dat");
//     if(CentralityLow.CompareTo("20")==0) inputfile.open("/home/mike/git_afterburner/AnalysisSoftware/TaskFlow/Theory/direct_photons_cent2040_Paquet2016.dat");
    
    if(inputfile.is_open()) cout << "file is open..." << endl;
    
    Float_t tpt[11], tyield[11], tyieldstat[11], tv2[11], tv2stat[11], tv3[11], tv3stat[11], tv2stathigh[11], tv2statlow[11];
    Int_t nlinesfile = 0;
    Int_t tn = 11;
    
    while(1){
      cout << nlinesfile << endl;
      inputfile >> tpt[nlinesfile] >> tyield[nlinesfile] >> tyieldstat[nlinesfile] >> tv2[nlinesfile] >> tv2stat[nlinesfile] >> tv3[nlinesfile] >> tv3stat[nlinesfile];
//       tv2stathigh[nlinesfile] = tv2[nlinesfile] + tv2stat[nlinesfile];
//       tv2statlow[nlinesfile] = tv2[nlinesfile] - tv2stat[nlinesfile];
      cout << tpt[nlinesfile] << " " << tyield[nlinesfile] << " " << tyieldstat[nlinesfile] << " " << tv2[nlinesfile] << " " << tv2stat[nlinesfile] << " " << tv3[nlinesfile] << " " << tv3stat[nlinesfile];
      if(!inputfile.good()) break;
      nlinesfile++;
    }
    
    TGraphAsymmErrors* graphDirectTheory = new TGraphAsymmErrors(tn,tpt,tv2,0,0,tv2stat,tv2stat);
    
    TCanvas* cTheory = new TCanvas("cTheory","",800,800);
    SetProperMargins();
    histoDirectEmpty->Draw();
    graphDirectTheory->SetFillColor(kBlue+2);
    graphDirectTheory->SetFillStyle(3001);
    graphDirectTheory->Draw("SAME 3");
    graphDirect1->Draw("SAME e1 p");
    
    cTheory->SaveAs(Form("PseudoExperiments_%s%s/v2_DirectGamma_Theory_%s%s.eps", CentralityLow.Data(), CentralityHigh.Data(), CentralityLow.Data(), CentralityHigh.Data()));
    
    inputfile.close();
  
  }
  
  
    
  //========================================================
  //Plotting and Saving
  //========================================================
  
  
  TLegend* leg = new TLegend(0.16,0.15,0.6,0.26);
  leg->SetTextSize(0.04);
  leg->AddEntry(graphDirect1,"v_{2}^{#gamma,direct}","lp");
  leg->SetBorderSize(0);
  
  TCanvas* c2 = new TCanvas("c2","",800,800);
  SetProperMargins();
  histoDirectEmpty->Draw();
  graphDirect1->Draw("SAME e1 p");
  DrawDirectInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  T1.DrawLatex(0.17, 0.94, "work in progress");
  
  TCanvas* c4 = new TCanvas("c4","",800,800);
  SetProperMargins();
  histoRGammaEmpty->Draw();
  graphRGamma->Draw("SAME e1 p");
  DrawRGammaInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  
  TCanvas* c5 = new TCanvas("c5","",800,800);
  SetProperMargins();
  histoInclusiveEmpty->Draw();
  DrawInfoLabelOnPlot(CentralityLow.Data(),CentralityHigh.Data(),2);
  T1.DrawLatex(0.17, 0.94, "work in progress");
  graphInclusive->Draw("SAME e1 p");
  histoCocktail->Draw("SAME");
  
  TLegend* legIncl = new TLegend(0.16,0.15,0.6,0.26);
  legIncl->SetTextSize(0.04);
  legIncl->AddEntry(graphInclusive,"v_{2}^{#gamma,inc}","lp");
  legIncl->AddEntry(histoCocktail,"v_{2}^{#gamma,dec}","lp");
  legIncl->SetBorderSize(0);
  legIncl->Draw();
  
  
  
  
//   leg->Draw();

  c2->SaveAs(Form("PseudoExperiments_%s%s/v2_DirectGamma_Errors_%s%s.eps", CentralityLow.Data(), CentralityHigh.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  c4->SaveAs(Form("PseudoExperiments_%s%s/RGamma_%s%s.eps", CentralityLow.Data(), CentralityHigh.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  c5->SaveAs(Form("PseudoExperiments_%s%s/v2Incl_v2Cocktail_%s%s.eps", CentralityLow.Data(), CentralityHigh.Data(), CentralityLow.Data(), CentralityHigh.Data()));
  
  TFile *output_File = new TFile(Form("PseudoExperiments_%s%s/PCM_DirectPhotonFlow_160724_%s%s.root", CentralityLow.Data(), CentralityHigh.Data(), CentralityLow.Data(), CentralityHigh.Data()),"RECREATE");
  histoInclusive->Write("PCM_v2_Inclusive");
  histoCocktail->Write("PCM_v2_Cocktail");
  graphDirect1->Write("PCM_v2_Direct_TotalError");
  graphRGamma->Write("comb_RGamma");
  output_File->Close();
  
}