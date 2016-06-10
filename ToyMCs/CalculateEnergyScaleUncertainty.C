#include <stdlib.h>
#include <iostream>
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
#include "TProfile2D.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
//#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
//#include "../CommonHeaders/FittingGammaConversion.h"
//#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
//#include "../CommonHeaders/ConversionFunctions.h"

#include "CalculateEnergyScaleUncertainty.h"


void CalculateEnergyScaleUncertainty(const Int_t ntotal=1e8){

  gRandom = new TRandom3(time(0));
  
  const Int_t nbins = 500;
  const Double_t histmin = 0;
  const Double_t histmax = 50;
  const Double_t dnmin = 0.5;
  const Double_t dnmax = 50;
  
  // function to describe energy uncertainty
  TF1* funcertainty = new TF1("funcertainty","pol0",dnmin,dnmax);
  funcertainty->SetParameter(0,1.01);
  
  // function for pT sampling (flat for good statistics)
  TF1* fpt = new TF1("fpt","pol0",dnmin,dnmax);
  fpt->SetParameter(0,1.0);
  
  // function for input spectrum shape
  TF1* fdndpt = new TF1("fdndpt","1/x^[0]",dnmin,dnmax);
  fdndpt->SetParameter(0,7);
  
  // unshifted histogram
  TH1D* fHdN = new TH1D("hdN","",nbins,histmin,histmax);

  // shifted histogram
  TH1D* fHdNShift = new TH1D("hdNShift","",nbins,histmin,histmax);
  
  Int_t counter = 0;
  
  for (int i = 0; i<ntotal; i++) { // loop over events
    if(counter % 100000 == 0)
      cout << counter << endl;
    
    // get random pt
    Double_t pt = fpt->GetRandom();
    
    // get weight
    Double_t weight = fdndpt->Eval(pt);
    
    // fill histograms
    // unshifted
    fHdN->Fill(pt,weight);
    
    // get shift
    Double_t ptshift = funcertainty->Eval(pt);
    
    // fill shifted histogram
    fHdNShift->Fill(pt*ptshift,weight);

    counter++;
  }
  
  // draw both histos and ratio
  
  // canvas
  TCanvas* c1 = new TCanvas("c1","",900,900);
  
  // top panel
  TPad* ptop = new TPad("ptop","",0,0.3,1,1);
  ptop->Draw();
  ptop->cd();
  gPad->SetLogy();
  
  // draw before and after
  fHdN->SetMarkerStyle(20);
  fHdN->SetMarkerColor(1);
  fHdN->SetLineColor(1);
  fHdN->SetStats(kFALSE);
  fHdN->SetXTitle("p_{T} (GeV/c)");
  fHdN->SetYTitle("N^{#pi^{0}}");
  fHdN->Draw("p");

  fHdNShift->SetMarkerStyle(21);
  fHdNShift->SetMarkerColor(2);
  fHdNShift->SetLineColor(2);
  fHdNShift->SetStats(kFALSE);
  fHdNShift->Draw("p,same");

  c1->cd();
  // bottom panel
  TPad* pbottom = new TPad("pbottom","",0,0,1,0.3);
  pbottom->Draw();
  pbottom->cd();

  // calculate ratio
  TH1D* fHdNRatio = (TH1D*) fHdNShift->Clone("fHdNRatio");
  fHdNRatio->Divide(fHdN);
  fHdNRatio->SetYTitle("shifted / regular");
  fHdNRatio->SetMarkerColor(1);
  fHdNRatio->SetLineColor(1);
  fHdNRatio->SetAxisRange(0.8,1.2,"Y");
  fHdNRatio->Draw("p");
  
}
