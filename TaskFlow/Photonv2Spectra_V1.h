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
#include "TPaveText.h"
#include "/home/mike/alicesw/aliphysics/master/src/PWG/FLOW/Base/AliFlowCommonHistResults.h"
#include "/home/mike/alicesw/aliphysics/master/src/PWG/FLOW/Base/AliFlowCommonHist.h"
#include "/home/mike/afterburner_v3/AnalysisSoftware/CommonHeaders/PlottingGammaConversionAdditional.h"
#include "/home/mike/afterburner_v3/AnalysisSoftware/CommonHeaders/PlottingGammaConversionHistos.h"

void SetStyle(Bool_t graypalette=kTRUE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
TH1D *IncluSPv0(TString filename, TString directnameQa, TString directnameQb, TString detectdir,bool reb = kFALSE);
TH1D *Spectra(TString filename, TString directnameQa, TString directnameQb, TString detectdir,bool reb = kFALSE);
void libs();

void libs()
{
  // load necessary libraries
  cout << " --- Loading libraries --- " << endl;
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libProof");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
}



//------------------------------------------------------------------
TH1D *PionSubtr(TH1D *v2Incl, TH1D *v2Pions,TH1D *Purity)
{
  TH1D *rat1 = (TH1D*)v2Incl->Clone();
  rat1->Divide(Purity);

  TH1D *rat2 = (TH1D*)v2Pions->Clone();
  rat2->Divide(Purity);
  rat2->Add(v2Pions,-1);

  TH1D *Final = (TH1D*)rat1->Clone();
  Final->Add(rat2,-1);

  return Final;
}
//-------------------------------------------------------------------
TH1D* RebinInPt(Int_t nPtIntervals, Double_t *ptInterval, Int_t *nMergedBins, TH1D *hist)
{
  // Rebin original histograms
  cout << " > Integrating differential flow " << endl;
  if (!hist){
    cout << endl;
    cout << " WARNING: hist is NULL in RebinInPt() !!!!" << endl;
    cout << endl;
    exit(0);
  }
  
  Double_t binWidthOld = hist->GetXaxis()->GetBinWidth(4);
  Int_t nBinsOld = hist->GetXaxis()->GetNbins();
  for (Int_t b = 1; b <= nBinsOld; b++){
    if (TMath::Abs(hist->GetXaxis()->GetBinWidth(b) - binWidthOld) > 1.e-44){
      cout << endl;
      cout << Form(" WARNING: %s have bins of unequal width !!!!", hist->GetName()) << endl;
      cout << "               Do not trust rebinning for high pt." << endl;
      cout << endl;
    }
  }
  
  if (binWidthOld < 1.e-44){
    cout << endl;
    cout << Form(" WARNING: in %s bin width is 0 !!!!", hist->GetName()) << endl;
    cout << "               Cannot do rebinning in pt." << endl;
    cout << endl;
    exit(0);
  }
  
  // Book rebinned histogram
  Int_t nBinsNew = 0;
  for (Int_t i = 0; i < nPtIntervals; i++){
    if (0 == nMergedBins[i]){
      continue;
    }
    Double_t xMin = TMath::Nint(ptInterval[i] / binWidthOld) * binWidthOld;
    Double_t xMax = TMath::Nint(ptInterval[i + 1] / binWidthOld) * binWidthOld;
    Int_t nBins = TMath::Nint((xMax - xMin) / binWidthOld);
    if (nBins <= 0){
      cout << endl;
      cout << Form(" WARNING: nBins <=0 when rebinning %s !!!!", hist->GetName()) << endl;
      cout << "               Check entries in array ptInterval." << endl;
      cout << endl;
      exit(0);
    }
    if (nBins % nMergedBins[i] == 0){
      nBinsNew += nBins / nMergedBins[i];
    }
    else{
      nBinsNew += (nBins / nMergedBins[i] + 1);
    }
  }
  
  const Int_t nBinsRebinned = nBinsNew;
  Double_t binEdges[nBinsRebinned + 1];
  Int_t counterForRebinnedBins = 0;
  for (Int_t i = 0; i < nPtIntervals; i++){
    if (0 == nMergedBins[i]){
      continue;
    }
    Double_t xMin = TMath::Nint(ptInterval[i] / binWidthOld) * binWidthOld;
    Double_t xMax = TMath::Nint(ptInterval[i + 1] / binWidthOld) * binWidthOld;
    Int_t nBins = TMath::Nint((xMax - xMin) / binWidthOld);
    if (nBins % nMergedBins[i] == 0){
      nBins = nBins / nMergedBins[i];
    }
    else{
      nBins = (nBins / nMergedBins[i] + 1);
    }
    for (Int_t b = 0; b < nBins; b++){
      binEdges[counterForRebinnedBins] = xMin + b * binWidthOld * nMergedBins[i];
      counterForRebinnedBins++;
    }
  }
  
  // Last bin edge:
  binEdges[counterForRebinnedBins] = hist->GetXaxis()->GetXmax();
  TH1D *temp = new TH1D("", "", nBinsRebinned, binEdges); // rebinned histogram
  for (Int_t br = 0; br < nBinsRebinned; br++){ // bins in rebinned histogram
    Double_t value = 0.;
    Double_t error = 0.;
    Double_t dSum1 = 0.; // sum value_i/(error_i)^2
    Double_t dSum2 = 0.; // sum 1/(error_i)^2
    Int_t startingBin = hist->FindBin(binEdges[br]);
    Int_t endingBin = hist->FindBin(binEdges[br + 1]);
    for (Int_t bo = startingBin; bo < endingBin; bo++){ // bins in original histogram
      value = hist->GetBinContent(bo);
      error = hist->GetBinError(bo);
      if (error > 0.){
        dSum1 += value / (error * error);
        dSum2 += 1. / (error * error);
      }
    }
    if (dSum2 > 0.){
      temp->SetBinContent(br + 1, dSum1 / dSum2);
      temp->SetBinError(br + 1, pow(1. / dSum2, 0.5));
    }
  }
  
  return temp;
}
//_____________________________________________________________________________
TGraphErrors* Change2Graph(TH1D * original)
{
    // cast TH1D* to TGraphErrors to enable shift in X
    cout << " > (TGraphErrors*)TH2D called " << endl;
    TGraphErrors *toReturn = new TGraphErrors(original);
    toReturn->GetXaxis()->SetTitle(original->GetXaxis()->GetTitle());
    toReturn->GetYaxis()->SetTitle(original->GetYaxis()->GetTitle());
    return toReturn;
}
//_____________________________________________________________________________
void ShiftGraphInX(TGraphErrors * Graph, double Shift)
{
    // shift graph in X for beatiful plots
    cout << " > Shifting plot in X " << endl;
    TString titleX = Graph->GetXaxis()->GetTitle();
    TString titleY = Graph->GetYaxis()->GetTitle();
    int ns = Graph->GetN();
    Double_t *xs = Graph->GetX();
    Double_t *ys = Graph->GetY();
    for (int i = 0; i != ns; ++i){
      Graph->SetPoint(i, xs[i] + Shift, ys[i]);
    }
    Graph->GetXaxis()->SetTitle(titleX.Data());
    Graph->GetYaxis()->SetTitle(titleY.Data());
}
//_____________________________________________________________________________

TH1D *IncluSPv0(TString filename, TString directnameQa, TString directnameQb, TString detectdir, bool reb){
  
  //Define the binning here
  const Int_t nPtIntervals = 200; // set some input values for RebininPt which will be used to calculate integrated flow
  Double_t ptInterval2[nPtIntervals + 1] = {0.0,0.9,2.7,3.3,4.1,4.6,7.0,8.0,14.0}; // in GeV
  Int_t nMergedBins2[nPtIntervals] = {1, 2, 3, 4, 5, 8, 10, 30};
  
  //------------Inclusive electron flow VZERO SP Qa----------------------------------
  TFile *FileName = new TFile(filename);
  TDirectoryFile *dirSP = (TDirectoryFile*)FileName->Get(detectdir);
  
  TList *listSPQa = (TList*)dirSP->Get(directnameQa);
  AliFlowCommonHistResults *lBB1Qa = (AliFlowCommonHistResults*) listSPQa->FindObject(Form("AliFlowCommonHistResults_SP"));
  TH1D *fFlowPOISPQa, *rebinnedFlowSPQa;
  fFlowPOISPQa = (TH1D*) lBB1Qa->GetHistDiffFlowPtPOI();
  fFlowPOISPQa->GetYaxis()->SetTitle("v_{2}{EP, |#Delta#eta| > 0.9}");
  fFlowPOISPQa->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  
  TList *listSPQb = (TList*)dirSP->Get(directnameQb);
  AliFlowCommonHistResults *lBB1Qb = (AliFlowCommonHistResults*) listSPQb->FindObject(Form("AliFlowCommonHistResults_SP"));
  TH1D *fFlowPOISPQb, *rebinnedFlowSPQb;
  fFlowPOISPQb = (TH1D*) lBB1Qb->GetHistDiffFlowPtPOI();
  fFlowPOISPQb->GetYaxis()->SetTitle("v_{2}{EP, |#Delta#eta| > 0.9}");
  fFlowPOISPQb->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  if(reb){
    rebinnedFlowSPQa = (TH1D*)(RebinInPt(nPtIntervals, ptInterval2, nMergedBins2, fFlowPOISPQa));
    rebinnedFlowSPQb = (TH1D*)(RebinInPt(nPtIntervals, ptInterval2, nMergedBins2, fFlowPOISPQb));
  }


  TH1D *fFlowPOISPQaclon =(TH1D*)fFlowPOISPQa->Clone();
  TH1D *v2IncSPVzero;
  
  if(reb){
    TH1D *rebinnedFlowSPQaclon =(TH1D*)rebinnedFlowSPQa->Clone();
    for(Int_t o = 0; o < rebinnedFlowSPQaclon->GetNbinsX(); o++) {
        (rebinnedFlowSPQaclon->GetBinContent(o+1)*rebinnedFlowSPQb->GetBinContent(o+1) > 0) ? rebinnedFlowSPQaclon->SetBinContent(o+1, TMath::Sqrt(rebinnedFlowSPQaclon->GetBinContent(o+1)*rebinnedFlowSPQb->GetBinContent(o+1))) : rebinnedFlowSPQaclon->SetBinContent(o+1, 0);
        //rebinnedFlowSPQaclon->SetBinError(o+1, TMath::Sqrt(TMath::Power(rebinnedFlowSPQaclon->GetBinError(o+1), 2)+TMath::Power(rebinnedFlowSPQb->GetBinError(o+1), 2)));
        rebinnedFlowSPQa->SetBinError(o+1, rebinnedFlowSPQb->GetBinError(o+1));
    }
    v2IncSPVzero = (TH1D*)rebinnedFlowSPQaclon->Clone();
  }else{
    for(Int_t o = 0; o < fFlowPOISPQaclon->GetNbinsX(); o++){
      (fFlowPOISPQaclon->GetBinContent(o+1)*fFlowPOISPQb->GetBinContent(o+1) > 0) ? fFlowPOISPQaclon->SetBinContent(o+1, TMath::Sqrt(fFlowPOISPQaclon->GetBinContent(o+1)*fFlowPOISPQb->GetBinContent(o+1))) : fFlowPOISPQaclon->SetBinContent(o+1, 0);
      //fFlowPOISPQaclon->SetBinError(o+1, TMath::Sqrt(TMath::Power(fFlowPOISPQaclon->GetBinError(o+1), 2)+TMath::Power(fFlowPOISPQb->GetBinError(o+1), 2)));
      fFlowPOISPQa->SetBinError(o+1, fFlowPOISPQb->GetBinError(o+1));
    }
    v2IncSPVzero = (TH1D*)fFlowPOISPQaclon->Clone();
  }
  
  return v2IncSPVzero;
}



//========================================================================================================

//========================================================================================================
TH1D *Spectra(TString filename, TString directnameQa, TString directnameQb, TString detectdir,bool reb){
  
  //------------Inclusive electron flow VZERO SP Qa----------------------------------
  TFile *FileName = new TFile(filename);
  TDirectoryFile *dirSP = (TDirectoryFile*)FileName->Get(detectdir);
  TList *listSPQa = (TList*)dirSP->Get(directnameQa);
  AliFlowCommonHist *lBB1Qa = (AliFlowCommonHist*) listSPQa->FindObject(Form("AliFlowCommonHist_SP"));

  TH1D *HistPtPOI, *Nev;
  HistPtPOI = (TH1D*) lBB1Qa->GetHistPtPOI();
  HistPtPOI->Sumw2();
  HistPtPOI->GetYaxis()->SetTitle("Counts/Nev");
  HistPtPOI->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");

  Nev = (TH1D*)lBB1Qa->GetHistAngleQ();
  cout << "Nev == " << Nev->GetEntries() << endl;

  if(reb == kTRUE){
    Double_t binsPt[28] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.3, 3.7, 4.1, 4.6, 5.4, 6.2, 7.0, 8.0};
    TH1D* HistPtPOIReb = (TH1D*)HistPtPOI->Rebin(27, "", binsPt);
    HistPtPOIReb->Scale(1./Nev->GetEntries(),"width");
    return HistPtPOIReb;
  }else{
    HistPtPOI->Scale(1./Nev->GetEntries(),"width");
    return HistPtPOI;
  }
}

//==========================
void SetStyle(Bool_t graypalette){
  cout << "Setting style!" << endl;
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}