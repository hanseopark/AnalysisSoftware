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
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

#include "CalculateEnergyScaleUncertainty.h"


void CalculateEnergyScaleUncertainty(TString fileName = "data_EMCAL-EMCALResultsFullCorrection_PP_2016_06_08_RetunedThresholds.root",
                                     Int_t ntotal=5e7, // number of "particles"
                                     const Int_t nfuncdndpt = 2, // switch dn/dpt function
                                     const Int_t nfuncdpt = 2,   // switch pt dependence of energy uncertainty
                                     const Float_t mass = 0.1349766, // mass of pi0
                                     Bool_t variablebinning = 1){

  gRandom = new TRandom3(time(0));

  const Int_t nbinsfinal = 33;
  Float_t fBinsPi02760GeVPtMerged[nbinsfinal+1]            = { 0.0, 0.4, 0.6, 0.8, 1.0,
    1.2, 1.4, 1.6, 1.8, 2.0,
    2.2, 2.4, 2.6, 3.0, 3.5,
    4.0, 5.0, 6.0, 7.0, 8.0,
    9.0, 10.0, 11.0, 12.0, 14.0,
    16.0, 18., 22., 26., 30.,
    35.0, 40., 50., 60.};

  
  const Int_t nbins = 500;
  const Double_t histmin = 0;
  const Double_t histmax = 50;
  const Double_t dnmin = 0.5;
  const Double_t dnmax = 50;
  
  // function to describe energy uncertainty
  TF1* funcertainty;
  
  cout << nfuncdpt << endl;

  TString dNName;
  TString uncName;
  
  switch(nfuncdpt){

    case 1:
      funcertainty = new TF1("funcertainty","pol0",dnmin,dnmax);
      funcertainty->SetParameter(0,1.01);
      break;

    case 2:
      funcertainty = new TF1("funcertainty","[0] + [1]/(1+[2]*x)",dnmin,dnmax);
      funcertainty->SetParameters(1.008,0.015,0.3);
      break;
    default:
      cout << "no input function for energy scale uncertainty defined, switching to constant, running only 1 event" << endl;
      funcertainty = new TF1("funcertainty","pol0",dnmin,dnmax);
      funcertainty->SetParameter(0,1);
      ntotal = 1;

  }
  
  // function for input spectrum shape
  TF1* fdndpt;
  cout << nfuncdndpt << endl;
  switch(nfuncdndpt){
    case 1:    // power law
      fdndpt = new TF1("fdndpt","([0]/x)^[1]",dnmin,dnmax);
      fdndpt->SetParameters(1,6);
      dNName = "power law";
      break;

    case 2:    // tsallis
      fdndpt = new TF1("fdndpt",
                       Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * x* pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mass,mass,mass,mass),dnmin,dnmax);
      fdndpt->SetParameters(2e11,7., 0.137);
      dNName = "Tsallis";
      break;

    case 3:    // bylinkin
      fdndpt = new TF1("fdndpt",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mass,mass,mass),dnmin,dnmax);
      fdndpt->SetParameters(450e8,0.3,1,0.3,8.);
      dNName = "Bylinkin";
      break;
      
    default:
      cout << "no input function for spectral shape defined, switching to constant, running only 1 event" << endl;
      fdndpt = new TF1("fdndpt","pol0",dnmin,dnmax);
      fdndpt->SetParameter(0,1);
      ntotal = 1;
  }

  // function for pT sampling (flat for good statistics)
  TF1* fpt = new TF1("fpt","pol0",dnmin,dnmax);
  fpt->SetParameter(0,1.0);
  
  TH1D* fHdN;
  TH1D* fHdNShift;

  if(variablebinning){
    // unshifted histogram
    fHdN = new TH1D("hdN","",nbinsfinal,fBinsPi02760GeVPtMerged);
    // shifted histogram
    fHdNShift = new TH1D("hdNShift","",nbinsfinal,fBinsPi02760GeVPtMerged);
  }
  else{
    // unshifted histogram
    fHdN = new TH1D("hdN","",nbins,histmin,histmax);
    // shifted histogram
    fHdNShift = new TH1D("hdNShift","",nbins,histmin,histmax);
  }
  
  fHdN->Sumw2();
  fHdNShift->Sumw2();
  
  // open file and get spectrum
  TFile* fileInput  = new TFile(fileName.Data());
  TDirectory* directoryEMCALPi0                           = (TDirectory*)fileInput->Get("Pi02.76TeV");
  TH1D* hpi0spectrum = (TH1D*) directoryEMCALPi0->Get("InvCrossSectionPi0")->Clone("hpi0spectrum");
  
  // make dN/dpt!
  for(Int_t ibin = 1;ibin<hpi0spectrum->GetNbinsX();ibin++){
    hpi0spectrum->SetBinContent(ibin,hpi0spectrum->GetBinContent(ibin)*hpi0spectrum->GetBinCenter(ibin));
    hpi0spectrum->SetBinError(ibin,hpi0spectrum->GetBinError(ibin)*hpi0spectrum->GetBinCenter(ibin));
  }
  
  // fit spectrum!
  hpi0spectrum->Fit(fdndpt,"RE","",2.5,16);
  
  // draw it
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  c1->SetLogy();
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.02);
  c1->SetTopMargin(0.02);
  hpi0spectrum->SetAxisRange(0,15.5,"X");
  hpi0spectrum->SetAxisRange(3e3,3e9,"Y");
  hpi0spectrum->SetXTitle("#it{p}_{T} (GeV/c)");
  hpi0spectrum->SetYTitle("#it{E} #frac{d^{3}#it{N}}{d#it{p}^{3}}(pb GeV^{-2}c^{3}");
  hpi0spectrum->SetTitleSize(0.03,"X");
  hpi0spectrum->SetTitleSize(0.03,"Y");
  hpi0spectrum->SetTitleOffset(1.2,"X");
  hpi0spectrum->SetTitleOffset(1.5,"Y");
  hpi0spectrum->SetMarkerStyle(20);
  hpi0spectrum->Draw("p");
  fdndpt->Draw("same");
  
  // add legend
  TLegend* linput = new TLegend(0.6,0.85,0.9,0.97);
  linput->SetHeader("ALICE, #sqrt{s} = 2.76 TeV");
  linput->AddEntry(hpi0spectrum,"pp #rightarrow #pi^{0}+X");
  linput->AddEntry(fdndpt,Form("fit with %s function",dNName.Data()));
  linput->SetBorderSize(0);
  linput->Draw();
  
  Int_t counter = 0;
  
  for (int i = 0; i<ntotal; i++) { // loop over events
    if(counter % 1000000 == 0)
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
  
  // divide by number of iterations
  fHdN->Scale(1./ntotal);
  fHdNShift->Scale(1./ntotal);
  
  // draw both histos and ratio
  
  // canvas
  TCanvas* c2 = new TCanvas("c2","",900,900);
  c2->cd(0);
  
  // top panel
  TPad* ptop = new TPad("ptop","",0,0.3,1,1);
  ptop->SetBottomMargin(0);
  ptop->SetTopMargin(0.02);
  ptop->SetRightMargin(0.02);
  ptop->Draw();
  ptop->cd();
  ptop->SetLogy();
  
  // divide by bin width
  for(Int_t ibin=1;ibin<fHdN->GetNbinsX();ibin++){
    fHdN->SetBinContent(ibin,fHdN->GetBinContent(ibin)/fHdN->GetBinWidth(ibin));
    fHdNShift->SetBinContent(ibin,fHdNShift->GetBinContent(ibin)/fHdNShift->GetBinWidth(ibin));
    fHdN->SetBinError(ibin,fHdN->GetBinError(ibin)/fHdN->GetBinWidth(ibin));
    fHdNShift->SetBinError(ibin,fHdNShift->GetBinError(ibin)/fHdNShift->GetBinWidth(ibin));
  }
  
  // draw before and after
  fHdN->SetMarkerStyle(20);
  fHdN->SetMarkerColor(1);
  fHdN->SetLineColor(1);
  fHdN->SetStats(kFALSE);
  fHdN->SetXTitle("#it{p}_{T} (GeV/c)");
  fHdN->SetYTitle("d#it{N}^{#pi^{0}}/d#it{p}_{T} (a.u.)");
  fHdN->SetAxisRange(histmin,histmax-0.01,"X");
  fHdN->SetAxisRange(2e-1,8e9,"Y");
  fHdN->Draw("p");

  fHdNShift->SetMarkerStyle(21);
  fHdNShift->SetMarkerColor(2);
  fHdNShift->SetLineColor(2);
  fHdNShift->SetStats(kFALSE);
  fHdNShift->Draw("p,same");

  TLegend* legendSpectra = new TLegend(0.4, 0.8, 0.95, 0.95);
  legendSpectra->SetMargin(0.12);
  legendSpectra->SetNColumns(1);
  legendSpectra->AddEntry(fHdN,Form("unshifted spectrum, %s function",dNName.Data()));
  legendSpectra->AddEntry(fHdNShift,Form("spectrum shifted with energy uncertainty %s",uncName.Data()));
  legendSpectra->SetBorderSize(0);
  legendSpectra->Draw();
  
  c2->cd(0);
  // bottom panel
  TPad* pbottom = new TPad("pbottom","",0,0,1,0.3);
  pbottom->SetTopMargin(0.01);
  pbottom->SetRightMargin(0.02);
  pbottom->SetBottomMargin(0.25);
  pbottom->Draw();
  pbottom->cd();

  // calculate ratio
  TH1D* fHdNRatio = (TH1D*) fHdNShift->Clone("fHdNRatio");
  fHdNRatio->SetAxisRange(histmin,histmax-0.01,"X");
  fHdNRatio->Divide(fHdN);
  fHdNRatio->SetYTitle("shifted / regular");
  fHdNRatio->SetXTitle("#it{p}_{T} (GeV/c)");
  fHdNRatio->SetMarkerColor(1);
  fHdNRatio->SetLineColor(1);
  fHdNRatio->SetAxisRange(0.81,1.19,"Y");
  fHdNRatio->SetLabelSize(0.07,"X");
  fHdNRatio->SetTitleSize(0.1,"X");
  fHdNRatio->SetLabelSize(0.07,"Y");
  fHdNRatio->SetTitleSize(0.08,"Y");
  fHdNRatio->SetTitleOffset(0.5,"Y");
  fHdNRatio->Draw("p");
  
  funcertainty->SetLineStyle(2);
  funcertainty->SetLineColor(1);
  funcertainty->Draw("same");

  TLegend* legendRatio = new TLegend(0.14, 0.29, 0.42, 0.48);
  legendRatio->SetMargin(0.12);
  legendRatio->SetNColumns(1);
  legendRatio->SetBorderSize(0);
  legendRatio->AddEntry(fHdNRatio,"ratio of shifted and regular spectrum");
  legendRatio->AddEntry(funcertainty,"energy uncertainty function");
  legendRatio->Draw();

}
