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
                                     Bool_t variablebinning = 1,
                                     const Int_t ipoln = 5){
  
  gRandom = new TRandom3(time(0));
  
  TString outputDir ="./plots/";
  TString suffix = "eps";
  
  gSystem->Exec("mkdir -p "+outputDir);

  const Int_t nbinsfinal = 33;
  Float_t fBinsPi02760GeVPtMerged[nbinsfinal+1]            = { 0.0, 0.4, 0.6, 0.8, 1.0,
    1.2, 1.4, 1.6, 1.8, 2.0,
    2.2, 2.4, 2.6, 3.0, 3.5,
    4.0, 5.0, 6.0, 7.0, 8.0,
    9.0, 10.0, 11.0, 12.0, 14.0,
    16.0, 18., 22., 26., 30.,
    35.0, 40., 50., 60.};
  
  const Int_t ntriggers = 6;
  const char* cTrigName[ntriggers] = {"INT1","INT7","EMC1","EMC7","EG2","EG1"};
  const Float_t minpttrig[ntriggers] = {1.4,1.4,6. ,3. ,6. ,11.};
  const Float_t maxpttrig[ntriggers] = {6. ,6. ,20.,11.,16.,20.};
  
  const Int_t nbins = 500;
  const Double_t histmin = 0;
  const Double_t histmax = 50;
  const Double_t dnmin = 0.5;
  const Double_t dnmax = 50;
  
  TString dNName;
  TString uncName;
  TString ratiofitname;
  
  // switch function for fit to mass ratio
  TF1* fratiofit;
  switch(ipoln){
    case 0: // pol0
      fratiofit = new TF1("fratiofit","pol0",dnmin,dnmax);
      ratiofitname = "pol0";
      break;
      
    case 1: // pol1
      fratiofit = new TF1("fratiofit","pol1",dnmin,dnmax);
      ratiofitname = "pol1";
      break;
      
    case 2: // pol2
      fratiofit = new TF1("fratiofit","pol2",dnmin,dnmax);
      ratiofitname = "pol2";
      break;
      
    case 3: // pol3
      fratiofit = new TF1("fratiofit","pol3",dnmin,dnmax);
      ratiofitname = "pol3";
      break;
      
    case 4: // pol4
      fratiofit = new TF1("fratiofit","pol4",dnmin,dnmax);
      ratiofitname = "pol4";
      break;
      
    case 5: // pol5
      fratiofit = new TF1("fratiofit","pol5",dnmin,dnmax);
      ratiofitname = "pol5";
      break;
      
    case 6: // pol6
      fratiofit = new TF1("fratiofit","pol6",dnmin,dnmax);
      ratiofitname = "pol6";
      break;
      
    case 7: // pol7
      fratiofit = new TF1("fratiofit","pol7",dnmin,dnmax);
      ratiofitname = "pol7";
      break;
      
    case 8: // pol8
      fratiofit = new TF1("fratiofit","pol8",dnmin,dnmax);
      ratiofitname = "pol8";
      break;
      
    case 9: // pol9
      fratiofit = new TF1("fratiofit","pol9",dnmin,dnmax);
      ratiofitname = "pol9";
      break;
      
    default:
      fratiofit = new TF1("fratiofit","pol0",dnmin,dnmax);
      ratiofitname = "pol0";
      ntotal = 1;
      cout << "maximum pol9! setting back to pol0, running only 1 event" << endl;
      
  }
  
  // function to describe energy uncertainty
  TF1* funcertainty;
  
  // switch different cases of energy uncertainty
  switch(nfuncdpt){
      
    case 1: // 1% constant
      funcertainty = new TF1("funcertainty","pol0",dnmin,dnmax);
      funcertainty->SetParameter(0,1.01);
      break;
      
    case 2: // 0.5% constant
      funcertainty = new TF1("funcertainty","pol0",dnmin,dnmax);
      funcertainty->SetParameter(0,1.005);
      break;
      
    case 8: // %larger at low pT, going to constant at high pT
      funcertainty = new TF1("funcertainty","[0] + [1]/(1+[2]*x)",dnmin,dnmax);
      funcertainty->SetParameters(1.008,0.015,0.3);
      break;
      
    case 9:
      funcertainty = (TF1*) fratiofit;
      break;
      
    default:
      cout << "no input function for energy scale uncertainty defined, switching to constant, running only 1 event" << endl;
      funcertainty = new TF1("funcertainty","pol0",dnmin,dnmax);
      funcertainty->SetParameter(0,1);
      ntotal = 1;
      
  }
  
  // switch function for input spectrum shape
  TF1* fdndpt;
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
  
  // histograms for simulated spectra
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
  
  // make dN/dpt! (multiply with pt)
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
  hpi0spectrum->SetYTitle("#it{E} #frac{d^{ 3}#it{N}}{d#it{p}^{ 3}} (pb GeV^{ -2}c^{ 3})");
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
  
  // TGraphErrors for mass ratio data/mc vs pT
  TGraphErrors* gMassRatios = new TGraphErrors(60);
  TGraphErrors* gMassRatio[ntriggers];
  
  for(Int_t i=0;i<ntriggers;i++){
    gMassRatio[i] = new TGraphErrors(30);
  }
  
  // get masses for data and MC
  TH1D* hPi0MassData[ntriggers];
  TH1D* hPi0MassMC[ntriggers];
  TH1D* hPi0MassRatio[ntriggers];
  
  for(Int_t i=0;i<ntriggers;i++){
    hPi0MassData[i] = (TH1D*) directoryEMCALPi0->Get(Form("Pi0_Mass_data_%s",cTrigName[i]))->Clone(Form("hPi0MassData%s",cTrigName[i]));
    hPi0MassMC[i] = (TH1D*) directoryEMCALPi0->Get(Form("Pi0_Mass_MC_%s",cTrigName[i]))->Clone(Form("hPi0MassMC%s",cTrigName[i]));
    hPi0MassRatio[i] = (TH1D*) hPi0MassData[i]->Clone(Form("hPi0MassRatio%s",cTrigName[i]));
    hPi0MassRatio[i]->Divide(hPi0MassMC[i]);
  }
  // fill all used points in a graph
  Int_t ipoint = 0;
  Int_t ipoints[ntriggers] = {0,0,0,0,0,0};
  for(Int_t i=0;i<ntriggers;i++){
    for(Int_t ibin=1;ibin<hPi0MassRatio[i]->GetNbinsX();ibin++){
      if(hPi0MassRatio[i]->GetBinCenter(ibin) > minpttrig[i] && hPi0MassRatio[i]->GetBinCenter(ibin) < maxpttrig[i]){
        
        gMassRatio[i]->SetPoint(ipoints[i],hPi0MassRatio[i]->GetBinCenter(ibin),hPi0MassRatio[i]->GetBinContent(ibin));
        gMassRatio[i]->SetPointError(ipoints[i],0.,hPi0MassRatio[i]->GetBinError(ibin));
        
        gMassRatios->SetPoint(ipoint,hPi0MassRatio[i]->GetBinCenter(ibin),hPi0MassRatio[i]->GetBinContent(ibin));
        gMassRatios->SetPointError(ipoint,0.,hPi0MassRatio[i]->GetBinError(ibin));
        
        ipoints[i]++;
        ipoint++;
      }
    }
  }
  
  gMassRatios->Fit(fratiofit,"RE","",dnmin,16.);
  
  TCanvas* c3 = new TCanvas("c3","",800,800);
  DrawGammaCanvasSettings( c3, 0.12, 0.015, 0.015, 0.07);
  TH2F * histoFrameRatio;
  histoFrameRatio = new TH2F(" histoFrameRatio"," histoFrameRatio",1000,0., 16.,10000,0.9,1.1);
  histoFrameRatio->SetStats(kFALSE);
  SetStyleHistoTH2ForGraphs(histoFrameRatio, "#it{p}_{T} (GeV/#it{c})","#it{M}^{ #pi^{0}}_{data} / #it{M}^{ #pi^{0}}_{MC}",
                            0.032,0.033, 0.032,0.035, 1.,1.5);
  histoFrameRatio->Draw();

  fratiofit->SetLineColor(kRed);
  fratiofit->Draw("same");

  TLegend* lmassrat = new TLegend(0.6,0.7,0.9,0.9);
  for(Int_t i=0;i<ntriggers;i++){
    gMassRatio[i]->SetMarkerStyle(20+i);
    gMassRatio[i]->SetFillColor(0);
    gMassRatio[i]->Draw("p");
    lmassrat->AddEntry(gMassRatio[i],cTrigName[i]);
  }
  lmassrat->AddEntry(ratiofitname,Form("%s fit to all triggers",ratiofitname.Data()));
  lmassrat->SetBorderSize(0);
  lmassrat->Draw();
  
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
  fHdN->SetYTitle("d#it{N}^{ #pi^{0}}/d#it{p}_{T} (a.u.)");
  fHdN->SetTitleOffset(1.1,"Y");
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
  
  c1->SaveAs(Form("%s/pi0_crosssection_fit%s.%s",outputDir.Data(),dNName.Data(),suffix.Data()));
  c2->SaveAs(Form("%s/pi0_massratio_data_mc_fit_%s.%s",outputDir.Data(),ratiofitname.Data(),suffix.Data()));
  c3->SaveAs(Form("%s/pi0_dndpt_shift_dnfit%s_uncfit%s.%s",outputDir.Data(),dNName.Data(),uncName.Data(),suffix.Data()));
  
}
