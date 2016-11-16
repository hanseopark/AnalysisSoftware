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
TF1* DoFitWithTsallis(TGraphErrors* graph, TString name, TString particle, Double_t p0, Double_t p1, Double_t p2);
void PlotInterpolationPtBins(TGraphErrors** gPtvSqrts,TGraphErrors** gPtvsEnergies, TF1** fPowerlaw, TGraphErrors* gRpPb,Int_t fColumnPlot, Int_t fRowPlot,TString namePlot);
void PlotAlphavsPt(TGraphErrors* gAlpha, TString method, TString thesisPlotLabel, TString namePlot);
void PlotWithFit(TCanvas *canvas, TH2F hist, TGraphErrors* graph, TGraphErrors* graphRebin, TF1* fit, TString name, TString outputDir, TString suffix, Color_t* colorTrigg);
void SetStyleTGraphErrorsForGraphs(TGraphErrors* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset = 1, Float_t yTitleOffset = 1, Int_t xNDivisions = 510, Int_t yNDivisions = 510);
void SetStyleTGraphAsymmErrorsForGraphs(TGraphAsymmErrors* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset = 1, Float_t yTitleOffset = 1, Int_t xNDivisions = 510, Int_t yNDivisions = 510);

void RemoveZerosAndPrint(TGraphErrors* graph, TString d){
  cout << d.Data() << endl;
  while (graph->GetY()[0] == 0. ) graph->RemovePoint(0);
  graph->Print();
}

TGraphErrors*      graphAlpha = 0x0;
TGraphErrors**     graphPtvsSqrts = 0x0;
TGraphErrors**     gPtvsEnergiesSystem = 0x0;
TF1**              fPowerlawSystem = 0x0;

void CalculateStatPlusSysErrors(TH1D* histStat, TH1D* histSys);

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
  TString dateForOutput = ReturnDateStringForOutput();
  TString outputDir = Form("SecondaryInterpolation/%s/%s/%s",dateForOutput.Data(),optionEnergy.Data(),suffix.Data());
  gSystem->Exec("mkdir -p "+outputDir);

  Color_t colorTrigg      [10]                = {kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2};

  //*************************************************************************************************
  //*************************** read input
  //*************************************************************************************************

  TFile* inputFile            = new TFile("ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root");
  TFile* inputLambda          = new TFile("ExternalInput/OtherParticles/Lambda-pp7TeV-Preliminary.root");
  TFile* inputAntiLambda      = new TFile("ExternalInput/OtherParticles/AntiLambda-pp7TeV-Preliminary.root");

  TH1D* hChargedKaon2760GeV   = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat2760GeV");
  TH1D* hChargedKaon7TeV      = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat7TeV");
  TH1D* hChargedKaon2760GeVSys= (TH1D*)inputFile->Get("histoChargedKaonSpecPubSyst2760GeV");
  TH1D* hChargedKaon7TeVSys   = (TH1D*)inputFile->Get("histoChargedKaonSpecPubSyst7TeV");
  CalculateStatPlusSysErrors(hChargedKaon2760GeV,hChargedKaon2760GeVSys);
  CalculateStatPlusSysErrors(hChargedKaon7TeV,hChargedKaon7TeVSys);

  TH1D* hLambda2760GeV        = (TH1D*)inputFile->Get("histoLambda1115SpecStat2760GeV");
  TH1D* hLambda2760GeVSys     = (TH1D*)inputFile->Get("histoLambda1115SpecSyst2760GeV");
  CalculateStatPlusSysErrors(hLambda2760GeV,hLambda2760GeVSys);
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

  TH1D* hProton2760GeV   = (TH1D*)inputFile->Get("histoProtonSpecPubStat2760GeV");
  TH1D* hProton7TeV      = (TH1D*)inputFile->Get("histoProtonSpecPubStat7TeV");
  TH1D* hProton2760GeVSys= (TH1D*)inputFile->Get("histoProtonSpecPubSyst2760GeV");
  TH1D* hProton7TeVSys   = (TH1D*)inputFile->Get("histoProtonSpecPubSyst7TeV");
  CalculateStatPlusSysErrors(hProton2760GeV,hProton2760GeVSys);
  CalculateStatPlusSysErrors(hProton7TeV,hProton7TeVSys);

  //*************************************************************************************************

  TGraphErrors* graphChargedKaon2760GeV = new TGraphErrors(hChargedKaon2760GeV);
  RemoveZerosAndPrint(graphChargedKaon2760GeV,"graphChargedKaon2760GeV");
  TGraphErrors* graphChargedKaon7TeV    = new TGraphErrors(hChargedKaon7TeV);
  RemoveZerosAndPrint(graphChargedKaon7TeV,"graphChargedKaon7TeV");

  TGraphErrors* graphLambda2760GeV      = new TGraphErrors(hLambda2760GeV);
  graphLambda2760GeV->RemovePoint(graphLambda2760GeV->GetN()-1);
  RemoveZerosAndPrint(graphLambda2760GeV,"graphLambda2760GeV");
  TGraphErrors* graphLambda7TeV         = new TGraphErrors(hLambda7TeV);
  RemoveZerosAndPrint(graphLambda7TeV,"graphLambda7TeV");

  TGraphErrors* graphProton2760GeV = new TGraphErrors(hProton2760GeV);
  RemoveZerosAndPrint(graphProton2760GeV,"graphProton2760GeV");
  TGraphErrors* graphProton7TeV = new TGraphErrors(hProton7TeV);
  RemoveZerosAndPrint(graphProton7TeV,"graphProton7TeV");

  //*************************************************************************************************
  //*************************** Fits
  //*************************************************************************************************

  TF1* fitKaon2760 = DoFitWithTsallis(graphChargedKaon2760GeV,"fitKaon2760","K", 0.2,7.,0.2);
  TF1* fitKaon7 = DoFitWithTsallis(graphChargedKaon7TeV,"fitKaon7","K", 0.2,7.,0.2);

  TF1* fitLambda2760 = DoFitWithTsallis(graphLambda2760GeV,"fitLambda2760","Lambda", 0.1,7.,0.2);
  TF1* fitLambda7 = DoFitWithTsallis(graphLambda7TeV,"fitLambda7","Lambda", 0.1,7.,0.2);

  TF1* fitProton2760 = DoFitWithTsallis(graphProton2760GeV,"fitProton2760","P", 0.15,7.,0.2);
  TF1* fitProton7 = DoFitWithTsallis(graphProton7TeV,"fitProton7","P", 0.15,7.,0.2);

  //*************************************************************************************************
  //*************************** rebin of 2.76 TeV spectra according to 7 TeV
  //*************************************************************************************************

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

  TGraphErrors* graphProton2760GeVrebin = (TGraphErrors*) graphProton7TeV->Clone("graphProton7TeV");
  for(Int_t i = 0; i<graphProton2760GeVrebin->GetN(); i++){
    graphProton2760GeVrebin->SetPoint(i,graphProton7TeV->GetX()[i],fitProton2760->Eval(graphProton7TeV->GetX()[i]));
    graphProton2760GeVrebin->SetPointError(i,graphProton7TeV->GetEX()[i],hProton2760GeV->GetBinError(hProton2760GeV->FindBin(graphProton7TeV->GetX()[i])));
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

  TGraphErrors*  graphProton8TeV = GetInterpolSpectrum2D(
                                   graphProton2760GeVrebin,
                                   graphProton7TeV,
                                   2760,
                                   7000,
                                   8000);

  if(graphAlpha && graphPtvsSqrts && gPtvsEnergiesSystem && fPowerlawSystem ){
    PlotInterpolationPtBins(graphPtvsSqrts,gPtvsEnergiesSystem,fPowerlawSystem,graphProton8TeV,9,7,Form("%s/Proton_Pt_vs_Sqrts.%s",outputDir.Data(),suffix.Data()));
    PlotAlphavsPt(graphAlpha, "pp", "Proton", Form("%s/Proton_Alpha_vs_Pt.%s", outputDir.Data(),suffix.Data()));
  }else{
    cout << "ERROR: NULL pointer - returning..." << endl;
    cout << graphAlpha << endl;
    cout << graphPtvsSqrts << endl;
    cout << gPtvsEnergiesSystem << endl;
    cout << fPowerlawSystem << endl;
    return;
  }

  //*************************************************************************************************
  //*************************** Fit extrapolations
  //*************************************************************************************************

  TF1* fitKaon8 = DoFitWithTsallis(graphKaon8TeV,"fitKaon8","K", 0.3,6.5,0.2);
  TF1* fitLambda8 = DoFitWithTsallis(graphLambda8TeV,"fitLambda8","Lambda", 0.1,6.5,0.2);
  TF1* fitProton8 = DoFitWithTsallis(graphProton8TeV,"fitProton8","P", 0.15,6.5,0.2);

  //*************************************************************************************************
  //*************************** plotting
  //*************************************************************************************************

  TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
  DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
  canvasDummy2->SetLogy();
  canvasDummy2->SetLogx();
  TH2F * histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.1,30.,1000,1e-10,2);
  SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","", 0.032,0.04, 0.04,0.04, 0.8,1.55);

  PlotWithFit(canvasDummy2, *histo2DDummy2, graphChargedKaon2760GeV, graphChargedKaon2760GeVrebin, fitKaon2760, "input_Kaon2760_withFit", outputDir, suffix, colorTrigg);
  PlotWithFit(canvasDummy2, *histo2DDummy2, graphLambda2760GeV, graphLambda2760GeVrebin, fitLambda2760, "input_Lambda2760_withFit", outputDir, suffix, colorTrigg);
  PlotWithFit(canvasDummy2, *histo2DDummy2, graphProton2760GeV, graphProton2760GeVrebin, fitProton2760, "input_Proton2760_withFit", outputDir, suffix, colorTrigg);

  PlotWithFit(canvasDummy2, *histo2DDummy2, graphChargedKaon7TeV, 0x0, fitKaon7, "input_Kaon7_withFit", outputDir, suffix, colorTrigg);
  PlotWithFit(canvasDummy2, *histo2DDummy2, graphLambda7TeV, 0x0, fitLambda7, "input_Lambda7_withFit", outputDir, suffix, colorTrigg);
  PlotWithFit(canvasDummy2, *histo2DDummy2, graphProton7TeV, 0x0, fitProton7, "input_Proton7_withFit", outputDir, suffix, colorTrigg);

  PlotWithFit(canvasDummy2, *histo2DDummy2, graphKaon8TeV, 0x0, fitKaon8, "input_Kaon8_withFit", outputDir, suffix, colorTrigg);
  PlotWithFit(canvasDummy2, *histo2DDummy2, graphLambda8TeV, 0x0, fitLambda8, "input_Lambda8_withFit", outputDir, suffix, colorTrigg);
  PlotWithFit(canvasDummy2, *histo2DDummy2, graphProton8TeV, 0x0, fitProton8, "input_Proton8_withFit", outputDir, suffix, colorTrigg);

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


  canvasDummy2->Clear();
  histo2DDummy2->DrawCopy();

  DrawGammaSetMarkerTGraphErr(graphProton2760GeVrebin, 20., 1., colorTrigg[1], colorTrigg[1], 1.4, kTRUE);
  graphProton2760GeVrebin->Draw("pEsame");
  DrawGammaSetMarkerTGraphErr(graphProton7TeV, 20., 1., colorTrigg[0], colorTrigg[0], 1.4, kTRUE);
  graphProton7TeV->Draw("pEsame");
  DrawGammaSetMarkerTGraphErr(graphProton8TeV, 21., 1., colorTrigg[2], colorTrigg[2], 1.4, kTRUE);
  graphProton8TeV->Draw("pEsame");

  fitProton2760->SetLineColor(colorTrigg[1]);
  fitProton2760->Draw("same");
  fitProton7->SetLineColor(colorTrigg[0]);
  fitProton7->Draw("same");
  fitProton8->SetLineColor(colorTrigg[2]);
  fitProton8->Draw("same");


  canvasDummy2->Update();
  canvasDummy2->Print(Form("%s/inputProton_withFit.%s",outputDir.Data(),suffix.Data()));



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


//________________________________________________________________________________________________________________________
TF1* DoFitWithTsallis(TGraphErrors* graph, TString name, TString particle, Double_t p0, Double_t p1, Double_t p2){

  cout << "-----------------------------------" << endl;
  cout << "fit: '" << name.Data() << "' for '" << particle.Data() << "'" << endl;

  Double_t paramFit[3] = {p0, p1, p2};
  TF1* fit = FitObject("l",name.Data(),particle.Data(),graph,graph->GetX()[0],graph->GetX()[graph->GetN()-1],paramFit,"QNRMEX0+");

  cout << "chi2/ndf: " << fit->GetChisquare()/fit->GetNDF() << endl;
  cout << endl;

  paramFit[0]=fit->GetParameter(0);
  paramFit[1]=fit->GetParameter(1);
  paramFit[2]=fit->GetParameter(2);
  fit = FitObject("l",name.Data(),particle.Data(),graph,graph->GetX()[0],graph->GetX()[graph->GetN()-1],paramFit,"QNRMEX0+");

  cout << "p0: " << fit->GetParameter(0) << " +- " << fit->GetParError(0) << endl;
  cout << "p1: " << fit->GetParameter(1) << " +- " << fit->GetParError(1) << endl;
  cout << "p2: " << fit->GetParameter(2) << " +- " << fit->GetParError(2) << endl;
  cout << "chi2/ndf: " << fit->GetChisquare()/fit->GetNDF() << endl;
  cout << "-----------------------------------" << endl;
  cout << endl;
  return fit;
}

//________________________________________________________________________________________________________________________
void CalculateStatPlusSysErrors(TH1D* histStat, TH1D* histSys){

  for(Int_t i=1; i<=histStat->GetNbinsX(); i++){
    histStat->SetBinError(i,TMath::Sqrt(TMath::Power(histStat->GetBinError(i),2)+TMath::Power(histSys->GetBinError(i),2)));
  }

  return;
}

//________________________________________________________________________________________________________________________
void PlotWithFit(TCanvas* canvas, TH2F hist, TGraphErrors* graph, TGraphErrors* graphRebin, TF1* fit, TString name, TString outputDir, TString suffix, Color_t *colorTrigg){

  canvas->Clear();
  hist.DrawCopy();

  DrawGammaSetMarkerTGraphErr(graph, 20., 1., colorTrigg[0], colorTrigg[0], 1.4, kTRUE);
  graph->Draw("pEsame");
  if(graphRebin){
    DrawGammaSetMarkerTGraphErr(graphRebin, 20., 1., colorTrigg[1], colorTrigg[1], 1.4, kTRUE);
    graphRebin->Draw("pEsame");
  }

  fit->SetLineColor(colorTrigg[1]);
  fit->Draw("same");

  canvas->Update();
  canvas->Print(Form("%s/%s.%s",outputDir.Data(),name.Data(),suffix.Data()));

  canvas->SetLogy(kFALSE);
  if(graphRebin){
    canvas->Clear();
    hist.GetYaxis()->SetRangeUser(0.5,1.5);
    hist.DrawCopy();

    TGraphErrors* graphRatio = CalculateGraphErrRatioToFit (graph, fit);
    TGraphErrors* graphRebinRatio = CalculateGraphErrRatioToFit (graphRebin, fit);

    DrawGammaSetMarkerTGraphErr(graphRatio, 20., 1., colorTrigg[0], colorTrigg[0], 1.4, kTRUE);
    graphRatio->Draw("pEsame");
    DrawGammaSetMarkerTGraphErr(graphRebinRatio, 20., 1., colorTrigg[1], colorTrigg[1], 1.4, kTRUE);
    graphRebinRatio->Draw("pEsame");

    canvas->Update();
    canvas->Print(Form("%s/%s_ratio.%s",outputDir.Data(),name.Data(),suffix.Data()));
  }else{
    canvas->Clear();
    hist.GetYaxis()->SetRangeUser(0,2.);
    hist.DrawCopy();

    TGraphErrors* graphRatio = CalculateGraphErrRatioToFit (graph, fit);

    DrawGammaSetMarkerTGraphErr(graphRatio, 20., 1., colorTrigg[0], colorTrigg[0], 1.4, kTRUE);
    graphRatio->Draw("pEsame");

    canvas->Update();
    canvas->Print(Form("%s/%s_ratio.%s",outputDir.Data(),name.Data(),suffix.Data()));
  }
  canvas->SetLogy();

  return;
}

//________________________________________________________________________________________________________________________
void SetStyleTGraphErrorsForGraphs(TGraphErrors* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset, Float_t yTitleOffset, Int_t xNDivisions, Int_t yNDivisions){
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
//________________________________________________________________________________________________________________________
void SetStyleTGraphAsymmErrorsForGraphs(TGraphAsymmErrors* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset, Float_t yTitleOffset, Int_t xNDivisions, Int_t yNDivisions){
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
