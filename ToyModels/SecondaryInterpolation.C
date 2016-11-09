/****************************************************************************************************************************
******          provided by Gamma Conversion Group, PWG-GA*****
******          Daniel Mühlheim, d.muehlheim@cern.ch
*****************************************************************************************************************************/

#include <Riostream.h>
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
#include "TGaxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/CombinationFunctions.h"
#include "TFitResultPtr.h"

TGraphErrors *GetInterpolSpectrum2D(TGraphErrors *g1, TGraphErrors *g2, Double_t d1, Double_t d2,Double_t dSqrts);
TH1D* ConvertYieldHisto(TH1D* input, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt);
void PlotInterpolationPtBins(TGraphErrors** gPtvSqrts,TGraphErrors** gPtvsEnergies, TF1** fPowerlaw, TGraphErrors* gRpPb,Int_t fColumnPlot, Int_t fRowPlot,TString namePlot);
void PlotAlphavsPt(TGraphErrors* gAlpha, TString method, TString thesisPlotLabel, TString namePlot);

void RemoveZerosAndPrint(TGraphErrors* graph, TString d){
  cout << d.Data() << endl;
  while (graph->GetY()[0] == 0. ) graph->RemovePoint(0);
  graph->Print();
}

void SetStyleTGraphErrorsForGraphs(TGraphErrors* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset = 1, Float_t yTitleOffset = 1, Int_t xNDivisions = 510, Int_t yNDivisions = 510){
    graph->GetXaxis()->SetTitle(XTitle);
    graph->GetYaxis()->SetTitle(YTitle);
    graph->SetTitle("");

    graph->GetXaxis()->SetLabelSize(xLableSize);
    graph->GetXaxis()->SetTitleSize(xTitleSize);
    graph->GetXaxis()->SetTitleOffset(xTitleOffset);
    graph->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->SetLabelFont(42);
    graph->GetXaxis()->SetTitleFont(62);
    graph->GetYaxis()->SetTitleFont(62);


    graph->GetYaxis()->SetDecimals();
    graph->GetYaxis()->SetLabelSize(yLableSize);
    graph->GetYaxis()->SetTitleSize(yTitleSize);
    graph->GetYaxis()->SetTitleOffset(yTitleOffset);
    graph->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);
}

void SetStyleTGraphAsymmErrorsForGraphs(TGraphAsymmErrors* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset = 1, Float_t yTitleOffset = 1, Int_t xNDivisions = 510, Int_t yNDivisions = 510){
    graph->GetXaxis()->SetTitle(XTitle);
    graph->GetYaxis()->SetTitle(YTitle);
    graph->SetTitle("");

    graph->GetXaxis()->SetLabelSize(xLableSize);
    graph->GetXaxis()->SetTitleSize(xTitleSize);
    graph->GetXaxis()->SetTitleOffset(xTitleOffset);
    graph->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->SetLabelFont(42);
    graph->GetXaxis()->SetTitleFont(62);
    graph->GetYaxis()->SetTitleFont(62);


    graph->GetYaxis()->SetDecimals();
    graph->GetYaxis()->SetLabelSize(yLableSize);
    graph->GetYaxis()->SetTitleSize(yTitleSize);
    graph->GetYaxis()->SetTitleOffset(yTitleOffset);
    graph->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);
}

TGraphErrors*      graphAlpha = 0x0;
TGraphErrors**     graphPtvsSqrts = 0x0;
TGraphErrors**     gPtvsEnergiesSystem = 0x0;
TF1**              fPowerlawSystem = 0x0;

void SecondaryInterpolation(TString suffix ="eps"){

  //*************************************************************************************************
  //*************************** general settings
  //*************************************************************************************************

  TH1::AddDirectory(kFALSE);
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  StyleSettingsThesis(suffix);
  SetPlotStyle();

  TString optionEnergy = "8TeV";
//  Float_t pT_low        = 0.3;
//  Float_t pT_high       = 20.;

//  cout << "pT low: " << pT_low << endl;
//  cout << "pT_high: " << pT_high << endl;

  TString dateForOutput = ReturnDateStringForOutput();
  TString outputDir = Form("SecondaryInterpolation/%s/%s/%s",dateForOutput.Data(),optionEnergy.Data(),suffix.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  Color_t colorTrigg      [10]                = {kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2};
  Color_t colorTriggShade [10]                = {kGray+1, kGray, kRed-6, kBlue-6, kGreen-8, kCyan-6, kViolet-8, kMagenta-8,  kRed-8, kBlue-8};
  Marker_t markerTrigg    [10]                = {20, 20, 21, 34, 29, 33, 21, 27, 28, 30 };
  Marker_t markerTriggMC  [10]                = {24, 24, 25, 28, 30, 27, 25, 27, 28, 30 };

  Size_t sizeTrigg        [10]                = {1.5, 1.5, 1.5, 2, 2.2, 2., 1.5, 2., 2.5, 1.5 };

  //*************************************************************************************************
  //*************************** read input
  //*************************************************************************************************

  TFile* inputFile            = new TFile("ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root");
  TFile* inputLambda          = new TFile("ExternalInput/OtherParticles/Lambda-pp7TeV-Preliminary.root");
  TFile* inputAntiLambda      = new TFile("ExternalInput/OtherParticles/AntiLambda-pp7TeV-Preliminary.root");

  TH1D* hChargedKaon2760GeV   = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat2760GeV");
  TH1D* hChargedKaon7TeV      = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat7TeV");

  TH1D* hLambda2760GeV        = (TH1D*)inputFile->Get("histoLambda1115SpecStat2760GeV");
  TH1D* hLambda7TeV           = 0x0;

  TH1D* hInput7TeVLambda      = (TH1D*)inputLambda->Get("fHistPtLambdaStatAndSystExceptNormalization");
  TH1D* hInput7TeVAntiLambda  = (TH1D*)inputAntiLambda->Get("fHistPtAntiLambdaStatAndSystExceptNormalization");
  if(hInput7TeVLambda && hInput7TeVAntiLambda){
    hInput7TeVLambda->Scale(0.5);
    hInput7TeVLambda->Add(hInput7TeVAntiLambda,0.5);
    hLambda7TeV = ConvertYieldHisto(hInput7TeVLambda, kTRUE, kTRUE, kFALSE, kFALSE);
  }else{
    cout << "\n\n\n\n\tWARNING: " << __LINE__ << " 7 TeV Lambda reference not found!!\n\n\n\n" << endl;
    return;
  }


  TGraphErrors* graphChargedKaon2760GeV = new TGraphErrors(hChargedKaon2760GeV);
  RemoveZerosAndPrint(graphChargedKaon2760GeV,"graphChargedKaon2760GeV");

  TGraphErrors* graphChargedKaon7TeV    = new TGraphErrors(hChargedKaon7TeV);
  RemoveZerosAndPrint(graphChargedKaon7TeV,"graphChargedKaon7TeV");

  TGraphErrors* graphLambda2760GeV      = new TGraphErrors(hLambda2760GeV);
  graphLambda2760GeV->RemovePoint(graphLambda2760GeV->GetN()-1);
  RemoveZerosAndPrint(graphLambda2760GeV,"graphLambda2760GeV");

  TGraphErrors* graphLambda7TeV         = new TGraphErrors(hLambda7TeV);
  RemoveZerosAndPrint(graphLambda7TeV,"graphLambda7TeV");

  //*************************************************************************************************
  //*************************** Fits
  //*************************************************************************************************

  Double_t paramFitKaon2760[3] = {0.15, 6., 0.13};
  TF1* fitKaon2760      = FitObject("l","fitKaon2760","K",graphChargedKaon2760GeV,graphChargedKaon2760GeV->GetX()[0],graphChargedKaon2760GeV->GetX()[graphChargedKaon2760GeV->GetN()-1],paramFitKaon2760,"QNRMEX0+");

  Double_t paramFitKaon7[3] = {0.15, 6., 0.13};
  TF1* fitKaon7      = FitObject("l","fitKaon7","K",graphChargedKaon7TeV,graphChargedKaon7TeV->GetX()[0],graphChargedKaon7TeV->GetX()[graphChargedKaon7TeV->GetN()-1],paramFitKaon7,"QNRMEX0+");

  Double_t paramFitLambda2760[3] = {0.15, 6., 0.13};
  TF1* fitLambda2760      = FitObject("l","fitLambda2760","Lambda",graphLambda2760GeV,graphLambda2760GeV->GetX()[0],graphLambda2760GeV->GetX()[graphLambda2760GeV->GetN()-1],paramFitLambda2760,"QNRMEX0+");

  Double_t paramFitLambda7[3] = {0.15, 6., 0.13};
  TF1* fitLambda7      = FitObject("l","fitLambda7","Lambda",graphLambda7TeV,graphLambda7TeV->GetX()[0],graphLambda7TeV->GetX()[graphLambda7TeV->GetN()-1],paramFitLambda7,"QNRMEX0+");

//  Double_t paramTCMKaon2760[5]  = { graphChargedKaon2760GeV->GetY()[0],0.1,
//                                    graphChargedKaon2760GeV->GetY()[20],0.6,3.0};
//  TF1* fitKaon2760  = FitObject("tcm","fitKaon2760","Pi0",graphChargedKaon2760GeV,graphChargedKaon2760GeV->GetX()[0],graphChargedKaon2760GeV->GetX()[graphChargedKaon2760GeV->GetN()-1],paramTCMKaon2760,"QNRMEX0+","", kFALSE);

//  Double_t paramTCMKaon7[5]  = { graphChargedKaon7TeV->GetY()[0],0.1,
//                                 graphChargedKaon7TeV->GetY()[15],0.6,3.0};
//  TF1* fitKaon7  = FitObject("tcm","fitKaon7","Pi0",graphChargedKaon7TeV,graphChargedKaon7TeV->GetX()[0],graphChargedKaon7TeV->GetX()[graphChargedKaon7TeV->GetN()-1],paramTCMKaon7,"QNRMEX0+","", kFALSE);

//  Double_t paramTCMLambda2760[5]  = { graphLambda2760GeV->GetY()[0],0.1,
//                                    graphLambda2760GeV->GetY()[15],0.6,3.0};
//  TF1* fitLambda2760  = FitObject("tcm","fitLambda2760","Pi0",graphLambda2760GeV,graphLambda2760GeV->GetX()[0],graphLambda2760GeV->GetX()[graphLambda2760GeV->GetN()-1],paramTCMLambda2760,"QNRMEX0+","", kFALSE);

//  Double_t paramTCMLambda7[5]  = { graphLambda7TeV->GetY()[0],0.1,
//                                 graphLambda7TeV->GetY()[15],0.6,3.0};
//  TF1* fitLambda7  = FitObject("tcm","fitLambda7","Pi0",graphLambda7TeV,graphLambda7TeV->GetX()[0],graphLambda7TeV->GetX()[graphLambda7TeV->GetN()-1],paramTCMLambda7,"QNRMEX0+","", kFALSE);


  TGraphErrors* graphChargedKaon2760GeVrebin = (TGraphErrors*) graphChargedKaon7TeV->Clone("graphChargedKaon2760GeVrebin");
  for(Int_t i = 0; i<graphChargedKaon2760GeVrebin->GetN(); i++){
    graphChargedKaon2760GeVrebin->SetPoint(i,graphChargedKaon7TeV->GetX()[i],fitKaon2760->Eval(graphChargedKaon7TeV->GetX()[i]));
    graphChargedKaon2760GeVrebin->SetPointError(i,graphChargedKaon7TeV->GetEX()[i],hChargedKaon2760GeV->GetBinError(hChargedKaon2760GeV->FindBin(graphChargedKaon7TeV->GetX()[i])));
  }

  TGraphErrors* graphLambda2760GeVrebin = (TGraphErrors*) graphLambda7TeV->Clone("graphLambda7TeV");
  for(Int_t i = 0; i<graphLambda2760GeVrebin->GetN(); i++){
    graphLambda2760GeVrebin->SetPoint(i,graphLambda7TeV->GetX()[i],fitLambda2760->Eval(graphLambda7TeV->GetX()[i]));
    graphLambda2760GeVrebin->SetPointError(i,graphLambda7TeV->GetEX()[i],hLambda2760GeV->GetBinError(hLambda2760GeV->FindBin(graphLambda7TeV->GetX()[i])));
  }

  //*************************************************************************************************
  //*************************** extrapolate spectra
  //*************************************************************************************************

  TGraphErrors*  graphKaon8TeV = GetInterpolSpectrum2D(
                                   graphChargedKaon2760GeVrebin,
                                   graphChargedKaon7TeV,
                                   2760,
                                   7000,
                                   8000);

  if(graphAlpha && graphPtvsSqrts && gPtvsEnergiesSystem && fPowerlawSystem ){
    PlotInterpolationPtBins(graphPtvsSqrts,gPtvsEnergiesSystem,fPowerlawSystem,graphKaon8TeV,9,7,Form("%s/Kaon_Pt_vs_Sqrts.%s",outputDir.Data(),suffix.Data()));
    PlotAlphavsPt(graphAlpha, "pp", "Lambda", Form("%s/Kaon_Alpha_vs_Pt.%s", outputDir.Data(),suffix.Data()));
  }else{
    cout << "ERROR: NULL pointer - returning..." << endl;
    cout << graphAlpha << endl;
    cout << graphPtvsSqrts << endl;
    cout << gPtvsEnergiesSystem << endl;
    cout << fPowerlawSystem << endl;
    return;
  }

  TGraphErrors*  graphLambda8TeV = GetInterpolSpectrum2D(
                                   graphLambda2760GeVrebin,
                                   graphLambda7TeV,
                                   2760,
                                   7000,
                                   8000);

  if(graphAlpha && graphPtvsSqrts && gPtvsEnergiesSystem && fPowerlawSystem ){
    PlotInterpolationPtBins(graphPtvsSqrts,gPtvsEnergiesSystem,fPowerlawSystem,graphLambda8TeV,7,5,Form("%s/Lambda_Pt_vs_Sqrts.%s",outputDir.Data(),suffix.Data()));
    PlotAlphavsPt(graphAlpha, "pp", "Lambda", Form("%s/Lambda_Alpha_vs_Pt.%s", outputDir.Data(),suffix.Data()));
  }else{
    cout << "ERROR: NULL pointer - returning..." << endl;
    cout << graphAlpha << endl;
    cout << graphPtvsSqrts << endl;
    cout << gPtvsEnergiesSystem << endl;
    cout << fPowerlawSystem << endl;
    return;
  }

  Double_t paramFitKaon8[3] = {0.15, 6., 0.13};
  TF1* fitKaon8      = FitObject("l","fitKaon8","K",graphKaon8TeV,graphKaon8TeV->GetX()[0],graphKaon8TeV->GetX()[graphKaon8TeV->GetN()-1],paramFitKaon8,"QNRMEX0+");

  Double_t paramFitLambda8[3] = {0.15, 6., 0.13};
  TF1* fitLambda8      = FitObject("l","fitLambda8","Lambda",graphLambda8TeV,graphLambda8TeV->GetX()[0],graphLambda8TeV->GetX()[graphLambda8TeV->GetN()-1],paramFitLambda8,"QNRMEX0+");

//  Double_t paramTCMKaon8[5]  = { graphKaon8TeV->GetY()[0],0.1,
//                                 graphKaon8TeV->GetY()[18],0.6,3.0};
//  TF1* fitKaon8  = FitObject("tcm","fitKaon8","Pi0",graphKaon8TeV,graphKaon8TeV->GetX()[0],graphKaon8TeV->GetX()[graphKaon8TeV->GetN()-1],paramTCMKaon8,"QNRMEX0+","", kFALSE);

//  Double_t paramTCMLambda8[5]  = { graphLambda8TeV->GetY()[0],0.1,
//                                   graphLambda8TeV->GetY()[15],0.6,3.0};
//  TF1* fitLambda8  = FitObject("tcm","fitLambda8","Pi0",graphLambda8TeV,graphLambda8TeV->GetX()[0],graphLambda8TeV->GetX()[graphLambda8TeV->GetN()-1],paramTCMLambda8,"QNRMEX0+","", kFALSE);


  //*************************************************************************************************
  //*************************** plotting
  //*************************************************************************************************

  TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
  DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
  canvasDummy2->SetLogy();
  canvasDummy2->SetLogx();
  TH2F * histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.1,30.,1000,1e-10,2);
  SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","", 0.032,0.04, 0.04,0.04, 0.8,1.55);
  histo2DDummy2->DrawCopy();

  DrawGammaSetMarkerTGraphErr(graphChargedKaon2760GeV, 20., 1., colorTrigg[0], colorTrigg[0], 1.4, kTRUE);
  graphChargedKaon2760GeV->Draw("pEsame");
  DrawGammaSetMarkerTGraphErr(graphChargedKaon2760GeVrebin, 20., 1., colorTrigg[1], colorTrigg[1], 1.4, kTRUE);
  graphChargedKaon2760GeVrebin->Draw("pEsame");

  fitKaon2760->SetLineColor(colorTrigg[1]);
  fitKaon2760->Draw("same");

  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/inputKaon2760_withFit.%s",outputDir.Data(),suffix.Data()));

  //**************************************************************************************************
  canvasDummy2->Clear();
  histo2DDummy2->DrawCopy();

  DrawGammaSetMarkerTGraphErr(graphChargedKaon7TeV, 20., 1., colorTrigg[0], colorTrigg[0], 1.4, kTRUE);
  graphChargedKaon7TeV->Draw("pEsame");

  fitKaon7->SetLineColor(colorTrigg[0]);
  fitKaon7->Draw("same");

  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/inputKaon7_withFit.%s",outputDir.Data(),suffix.Data()));

  //**************************************************************************************************
  canvasDummy2->Clear();
  histo2DDummy2->DrawCopy();

  DrawGammaSetMarkerTGraphErr(graphLambda2760GeV, 20., 1., colorTrigg[0], colorTrigg[0], 1.4, kTRUE);
  graphLambda2760GeV->Draw("pEsame");
  DrawGammaSetMarkerTGraphErr(graphLambda2760GeVrebin, 20., 1., colorTrigg[1], colorTrigg[1], 1.4, kTRUE);
  graphLambda2760GeVrebin->Draw("pEsame");

  fitLambda2760->SetLineColor(colorTrigg[1]);
  fitLambda2760->Draw("same");

  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/inputLambda2760_withFit.%s",outputDir.Data(),suffix.Data()));

  //**************************************************************************************************
  canvasDummy2->Clear();
  histo2DDummy2->DrawCopy();

  DrawGammaSetMarkerTGraphErr(graphLambda7TeV, 20., 1., colorTrigg[0], colorTrigg[0], 1.4, kTRUE);
  graphLambda7TeV->Draw("pEsame");

  fitLambda7->SetLineColor(colorTrigg[0]);
  fitLambda7->Draw("same");

  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/inputLambda7_withFit.%s",outputDir.Data(),suffix.Data()));


  //**************************************************************************************************
  //**************************************************************************************************
  canvasDummy2->Clear();
  histo2DDummy2->DrawCopy();

  DrawGammaSetMarkerTGraphErr(graphChargedKaon2760GeVrebin, 20., 1., colorTrigg[1], colorTrigg[1], 1.4, kTRUE);
  graphChargedKaon2760GeVrebin->Draw("pEsame");
  DrawGammaSetMarkerTGraphErr(graphChargedKaon7TeV, 20., 1., colorTrigg[0], colorTrigg[0], 1.4, kTRUE);
  graphChargedKaon7TeV->Draw("pEsame");
  DrawGammaSetMarkerTGraphErr(graphKaon8TeV, 21., 1., colorTrigg[2], colorTrigg[2], 1.4, kTRUE);
  graphKaon8TeV->Draw("pEsame");

  fitKaon2760->SetLineColor(colorTrigg[1]);
  fitKaon2760->Draw("same");
  fitKaon7->SetLineColor(colorTrigg[0]);
  fitKaon7->Draw("same");
  fitKaon8->SetLineColor(colorTrigg[2]);
  fitKaon8->Draw("same");

  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/inputKaon_withFit.%s",outputDir.Data(),suffix.Data()));

  canvasDummy2->Clear();
  histo2DDummy2->DrawCopy();

  DrawGammaSetMarkerTGraphErr(graphLambda2760GeVrebin, 20., 1., colorTrigg[1], colorTrigg[1], 1.4, kTRUE);
  graphLambda2760GeVrebin->Draw("pEsame");
  DrawGammaSetMarkerTGraphErr(graphLambda7TeV, 20., 1., colorTrigg[0], colorTrigg[0], 1.4, kTRUE);
  graphLambda7TeV->Draw("pEsame");
  DrawGammaSetMarkerTGraphErr(graphLambda8TeV, 21., 1., colorTrigg[2], colorTrigg[2], 1.4, kTRUE);
  graphLambda8TeV->Draw("pEsame");

  fitLambda2760->SetLineColor(colorTrigg[1]);
  fitLambda2760->Draw("same");
  fitLambda7->SetLineColor(colorTrigg[0]);
  fitLambda7->Draw("same");
  fitLambda8->SetLineColor(colorTrigg[2]);
  fitLambda8->Draw("same");


  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/inputLambda_withFit.%s",outputDir.Data(),suffix.Data()));



  return;
}


//________________________________________________________________________________________________________________________
TGraphErrors *GetInterpolSpectrum2D(TGraphErrors *g1, TGraphErrors *g2, Double_t d1, Double_t d2,Double_t dSqrts)
{
  if(!g1) return 0x0;
  if(!g2) return 0x0;

  TGraphErrors  *gInterpol     = new TGraphErrors(g1->GetN());
  TGraphErrors  *gAlpha        = new TGraphErrors(g1->GetN());
  TGraphErrors** gPtvsSqrts    = new TGraphErrors*[g1->GetN()];
  TGraphErrors** gPtvsEnergies = new TGraphErrors*[g1->GetN()];
  TF1**          fPowerlawFits = new TF1*[g1->GetN()];

  for(Int_t i = 0; i < g1->GetN(); i++){
    TGraphErrors *grint = new TGraphErrors(1);
    grint->SetPoint(0, dSqrts, 0);

    TGraphErrors *gToFit = new TGraphErrors(2);
    gToFit->SetPoint(0, d1, g1->GetY()[i]);
    gToFit->SetPointError(0, 0, g1->GetEY()[i]);

    gToFit->SetPoint(1, d2, g2->GetY()[i]);
    gToFit->SetPointError(1, 0, g2->GetEY()[i]);

    TF1 *fPowerlaw = new TF1("fPowerlaw","[0]*x^([1])", 0,10000);
    if(i==0){
      fPowerlaw->SetParameters(0, 0.005);
      fPowerlaw->SetParameters(1, 0.13);
    }else{
      fPowerlaw->SetParameters(0, 0.1);
      fPowerlaw->SetParameters(1, 2.0);
    }

    for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlaw,"Q");
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
    Double_t alpha  = fPowerlaw->GetParameter(1);
    cout << g1->GetX()[i] << ": " << fPowerlaw->GetParameter(0) << " - " << fPowerlaw->GetParameter(1) << endl;

    gInterpol->SetPoint(i, g1->GetX()[i],fPowerlaw->Eval(dSqrts));
    gInterpol->SetPointError(i, 0, grint->GetEY()[0]);

    gAlpha->SetPoint(i, g1->GetX()[i],alpha);
    //gAlpha->SetPointError(i, 0,alphaE);

    gPtvsSqrts[i]= new TGraphErrors(1);
    gPtvsSqrts[i]->SetPoint(0,dSqrts,gInterpol->GetY()[i]);
    //gPtvsSqrts[i]->Sort();

    fPowerlawFits[i] = fPowerlaw;
    gPtvsEnergies[i] = gToFit;

    gToFit->SetPoint(2, dSqrts, gInterpol->GetY()[i]);
    gToFit->SetPointError(2, 0, gInterpol->GetEY()[i]);

    delete grint;
  }

  graphAlpha          = gAlpha;
  graphPtvsSqrts      = new TGraphErrors*[g1->GetN()];
  fPowerlawSystem	  = new TF1*[g1->GetN()];
  gPtvsEnergiesSystem = new TGraphErrors*[g1->GetN()];

  for ( Int_t i = 0; i < g1->GetN(); i++ ){
    graphPtvsSqrts[i] = gPtvsSqrts[i];
    fPowerlawSystem[i] = fPowerlawFits[i];
    gPtvsEnergiesSystem[i] = gPtvsEnergies[i];
  }

  return gInterpol;
}


//________________________________________________________________________________________________________________________
TH1D* ConvertYieldHisto(TH1D* input, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt){
    if (!input) {
        cout << "Error: Histogram is NULL" << endl;
        return NULL;
    }

    Int_t nBins                 = input->GetNbinsX();
    Double_t newValue           = 0;
    Double_t newErrorValue      = 0;
    Double_t correctionValue    = 1;

    //correct by 2pi if specified
    if (DivideBy2pi) input->Scale(1/(2*TMath::Pi()));
    if (MultiplyBy2pi) input->Scale(2*TMath::Pi());

    for(Int_t i=0;i<nBins;i++){

        //correct by 1/Pt if specified
        if(DivideByPt)    correctionValue  = 1/(input->GetBinCenter(i+1));
        if(MultiplyByPt)  correctionValue  = input->GetBinCenter(i+1);

        //set the value and error of the bin
        input->SetBinContent(i+1,input->GetBinContent(i+1)*correctionValue);
        input->SetBinError(i+1,input->GetBinError(i+1)*correctionValue);
    }

    return input;
}

//________________________________________________________________________________________________________________________
void PlotInterpolationPtBins(TGraphErrors** gPtvSqrts,TGraphErrors** gPtvsEnergies, TF1** fPowerlaw, TGraphErrors* gRpPb,Int_t fColumnPlot, Int_t fRowPlot,TString namePlot){

    TGaxis::SetMaxDigits(3);
    TString nameCanvas = "";
    TString namePad    = "";

    TCanvas * canvasPtvsSqrts     = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasPtvsSqrts->SetTopMargin(0.00);
    canvasPtvsSqrts->SetBottomMargin(0.00);
    canvasPtvsSqrts->SetRightMargin(0.0);
    canvasPtvsSqrts->SetLeftMargin(0.00);

    TPad * padPtvsSqrts           = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padPtvsSqrts->SetFillColor(0);
    padPtvsSqrts->GetFrame()->SetFillColor(0);
    padPtvsSqrts->SetBorderMode(0);
    padPtvsSqrts->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padPtvsSqrts->SetLeftMargin(0.);
    padPtvsSqrts->SetRightMargin(0.);
    padPtvsSqrts->SetTopMargin(0.);
    padPtvsSqrts->SetBottomMargin(0.);
    padPtvsSqrts->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    Int_t place = 0;
    for(Int_t iPt = 0; iPt <gRpPb->GetN(); iPt++){
      place++;
      padPtvsSqrts->cd(place);
      padPtvsSqrts->cd(place)->SetTopMargin(0.12);
      padPtvsSqrts->cd(place)->SetBottomMargin(0.15);
      padPtvsSqrts->cd(place)->SetRightMargin(0.035);
      //padPtvsSqrts->cd(place)->SetLeftMargin(0.25);

      int remaining           = (place-1)%fColumnPlot;
      if (remaining > 0) padPtvsSqrts->cd(place)->SetLeftMargin(0.15);
      else padPtvsSqrts->cd(place)->SetLeftMargin(0.25);
      DrawGammaSetMarkerTGraphErr(gPtvSqrts[iPt],21,1.5, kRed , kRed);
      DrawGammaSetMarkerTGraphErr(gPtvsEnergies[iPt],20,1.5, kBlack , kBlack);

      gPtvsEnergies[iPt]->GetXaxis()->SetTitle("#sqrt{s}");
      gPtvsEnergies[iPt]->GetXaxis()->SetTitleSize(0.065);
      gPtvsEnergies[iPt]->GetXaxis()->SetLabelSize(0.06);
      gPtvsEnergies[iPt]->GetYaxis()->SetTitle("invariant cross section");
      gPtvsEnergies[iPt]->GetYaxis()->SetTitleSize(0.065);
      gPtvsEnergies[iPt]->GetYaxis()->SetLabelSize(0.06);
      gPtvsEnergies[iPt]->GetXaxis()->SetNdivisions(308,kTRUE);
      gPtvsEnergies[iPt]->GetYaxis()->SetNdivisions(304,kTRUE);
      gPtvsEnergies[iPt]->GetXaxis()->SetLabelOffset(0.015);
      gPtvsEnergies[iPt]->GetYaxis()->SetLabelOffset(0.01);

      gPtvsEnergies[iPt]->Draw("ap");
      gPtvSqrts[iPt]->Draw("same, p");
      fPowerlaw[iPt]->SetLineColor(kBlue);
      fPowerlaw[iPt]->SetLineWidth(2);
      fPowerlaw[iPt]->Draw("same");

      TString Title = Form("#it{p}_{T} = %3.2f GeV/#it{c} ",gRpPb->GetX()[iPt]);

      if(Title.Length() > 0){
        gPtvsEnergies[iPt]->SetTitle("");

        Double_t xMin,yMin;
        yMin = 0.78;
        if ( remaining > 0 ) xMin = 0.20;
        else xMin = 0.30;

        TLatex *alice = new TLatex(xMin,yMin,Form("%s",Title.Data())); // Bo: this was
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(0.062);
        alice->Draw();
      }
    }

    canvasPtvsSqrts->Print(namePlot.Data());
    delete padPtvsSqrts;
    delete canvasPtvsSqrts;
}

//________________________________________________________________________________________________________________________
void PlotAlphavsPt(TGraphErrors* gAlpha, TString method, TString thesisPlotLabel, TString namePlot){

    TCanvas * canvasAlphavsPt     = new TCanvas("AlphavsPt","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasAlphavsPt,  0.1, 0.02, 0.03, 0.1);
    SetStyleTGraphErrorsForGraphs(gAlpha,"pT","#alpha", 0.04,0.04, 0.04,0.04, 1.,1.2, 512, 512);
    DrawGammaSetMarkerTGraphErr(gAlpha,21,1.5, kRed , kRed);
    gAlpha->GetYaxis()->CenterTitle();
    gAlpha->Draw("ap");

    TLatex *labelThesis = new TLatex(0.15,0.90,thesisPlotLabel.Data());
    SetStyleTLatex( labelThesis, 0.04,4);
    labelThesis->Draw();

    canvasAlphavsPt->Print(namePlot.Data());
    delete canvasAlphavsPt;
}
