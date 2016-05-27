// provided to run with output from PWGGA/GammaConv/macros/AddTask_GammaPureMC.C
// takes pT hard binned pythia MC output, calculates weights, adds spectra from different bins and plots particle spectra and acceptance

// run for example as
// root -b -l -q 'AnalysePythiaMB.C+("grid/Legotrain_vAN-20150825-Pythia6/PythiaAnalysisResults","Pythia8","Perugia2011","Pi0")'

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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

#include "AnalysePythiaMB.h"

TString acceptanceOf = "";
TString ptHardRange = "";
TString fDetectionProcess ="";

//**********************************************************************************
//******************* return minimum for 2 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestEntryIn2D(TH2* histo){
  Double_t minimum = 1;
  for (Int_t i = 2; i<histo->GetNbinsX(); i++){
    for (Int_t j = 2; j<histo->GetNbinsY(); j++){
      if (histo->GetBinContent(i,j) < minimum && histo->GetBinContent(i,j) > 0){
        minimum = histo->GetBinContent(i,j);
      }
    }
  }
  return minimum;
}

//**********************************************************************************
//******************* return maximum for 2 D histo  ********************************
//**********************************************************************************
Double_t FindLargestEntryIn2D(TH2* histo){
  Double_t maximum = 1;
  for (Int_t i = 2; i<histo->GetNbinsX(); i++){
    for (Int_t j = 2; j<histo->GetNbinsY(); j++){
      if (histo->GetBinContent(i,j) > maximum && histo->GetBinContent(i,j) > 0){
        maximum = histo->GetBinContent(i,j);
      }
    }
  }
  return maximum;
}

//**********************************************************************************
//******************* return minimum for 1 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestEntryIn1D(TH1D* histo){
  Double_t minimum = 1;
  for (Int_t i = 1; i<histo->GetNbinsX(); i++){
    if (histo->GetBinContent(i) < minimum && histo->GetBinContent(i) > 0){
      minimum = histo->GetBinContent(i);
    }
  }
  return minimum;
}

//**********************************************************************************
//******************* return maximum for 1 D histo  ********************************
//**********************************************************************************
Double_t FindLargestEntryIn1D(TH1D* histo){
  Double_t maximum = 1;
  for (Int_t i = 1; i<histo->GetNbinsX(); i++){
    if (histo->GetBinContent(i) > maximum ){
      maximum = histo->GetBinContent(i);
    }
  }
  return maximum;
}

//**********************************************************************************
//******************* Standardized plotting of 2D plots ****************************
//**********************************************************************************
void PlotStandard2D( TH2* histo2D, TString nameOutput, TString title, TString xTitle, TString yTitle, Bool_t kRangeY, Double_t startY, Double_t endY, Bool_t kRangeX, Double_t startX, Double_t endX, Int_t logX, Int_t logY, Int_t logZ, Float_t* floatLogo, Int_t canvasSizeX = 500, Int_t canvasSizeY = 500, TString generator ="" , TString period ="",Bool_t bPlotPtHard = 1){
  TCanvas * canvasStandard = new TCanvas("canvasStandard","",10,10,canvasSizeX,canvasSizeY);  // gives the page size
  canvasStandard->SetLogx(logX);
  canvasStandard->SetLogy(logY);
  canvasStandard->SetLogz(logZ);
  canvasStandard->SetRightMargin(0.12);
  canvasStandard->SetLeftMargin(0.12);
  canvasStandard->SetBottomMargin(0.1);
  canvasStandard->SetTopMargin(0.04);
  canvasStandard->cd();
  histo2D->SetTitle("");
  DrawAutoGammaHisto2D(   histo2D,
                       title.Data(), xTitle.Data(), yTitle.Data(),"",kRangeY, startY, endY, kRangeX, startX, endX);
  histo2D->GetXaxis()->SetTitleOffset(1.05);
  //    cout << histo2D->GetYaxis()->GetTitleOffset() << endl;
  histo2D->GetYaxis()->SetTitleOffset(1.35);
  if (logX==1){
    //       cout << histo2D->GetXaxis()->GetLabelOffset() << endl;
    histo2D->GetXaxis()->SetLabelOffset(0.);
  }
  Float_t maximumZ = FindLargestEntryIn2D(histo2D)*5;
  Float_t minimumZ = FindSmallestEntryIn2D(histo2D);
  histo2D->SetAxisRange(minimumZ,maximumZ,"Z");
  
  histo2D->Draw("colz");
  DrawLabelsEvents(floatLogo[0],floatLogo[1],floatLogo[2], 0.00, "pp", "Pythia", period);
  TLatex *detprocess = 	new TLatex(floatLogo[0], floatLogo[1] - 3.2*floatLogo[2], fDetectionProcess);
  detprocess->SetNDC();
  detprocess->SetTextColor(1);
  detprocess->SetTextSize(floatLogo[2]);
  detprocess->Draw();
  TLatex *ptHardBin = 	new TLatex(floatLogo[0], floatLogo[1] - 4.25*floatLogo[2], ptHardRange);
  ptHardBin->SetNDC();
  ptHardBin->SetTextColor(1);
  ptHardBin->SetTextSize(floatLogo[2]);
  if(bPlotPtHard){
    ptHardBin->Draw();
  }
  
  canvasStandard->Update();
  canvasStandard->SaveAs(nameOutput.Data());
  delete canvasStandard;
}


//**********************************************************************************
//******************* Main function ************************************************
//**********************************************************************************

void AnalysePythiaMB(TString filemask="",  // should read something like "$PATHTODIRECTORY/PythiaAnalysisResultsBin"
                         TString generator="",
                         TString tune="",
                         TString particle="Pi0",
                         TString energy="2760GeV"){
  
  TH1::AddDirectory(kFALSE);
  
  const Int_t nbinsfinal = 33;
  Double_t fBinsPi02760GeVPtMerged[nbinsfinal+1]            = { 0.0, 0.4, 0.6, 0.8, 1.0,
    1.2, 1.4, 1.6, 1.8, 2.0,
    2.2, 2.4, 2.6, 3.0, 3.5,
    4.0, 5.0, 6.0, 7.0, 8.0,
    9.0, 10.0, 11.0, 12.0, 14.0,
    16.0, 18., 22., 26., 30.,
    35.0, 40., 50., 60.};
  
  //***************************************************************************************************************
  //************************************ Layouting preparations & general setup ***********************************
  //***************************************************************************************************************
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  
  StyleSettingsThesis();
  SetPlotStyle();
  
  const Int_t nreb = 4;
  
  const Int_t MaxNumberOfFiles = 12;
  
  Color_t colorBins[12]		= {kBlack, kRed+2, kBlue+2, kGreen+2, kCyan+2, kViolet, kMagenta+2, kGray+1, kRed-2, kBlue-2, kViolet+2, kCyan-2};
  Color_t colorBinsShade[12]	= {kGray+1, kRed-6, kBlue-6, kGreen-8, kCyan-6, kViolet-8, kMagenta-8, kGray, kRed-8, kBlue-8, kViolet-6, kCyan-8};
  Marker_t markerBins[12]		= {20, 21, 33, 34, 29, 24, 25, 27, 28, 30, 20, 21};
  
  
  Float_t floatLocationRightUp2D[4] = {0.45,0.95,0.035, 0.02};
  Float_t floatLocationLeftDown2D[4] = {0.15,0.25,0.035, 0.02};
  Float_t floatLocationRightDown2D[4] = {0.45,0.25,0.035, 0.02};
  
  TString outputDir ="./plots/";
  TString suffix = "eps";
  
  gSystem->Exec("mkdir -p "+outputDir);
  
  //fDetectionProcess 			= ReturnFullTextReconstructionProcess(mode);
  
  TString particlelegend = "";
  // translate particle string to plotting string
  if(particle.CompareTo("Pi0") == 0){
    particlelegend = "#pi^{0}";
  }
  else if(particle.CompareTo("Eta") == 0){
    particlelegend = "#eta";
  }
  else if(particle.CompareTo("PiPl") == 0){
    particlelegend = "#pi^{+}";
  }
  else if(particle.CompareTo("PiMi") == 0){
    particlelegend = "#pi^{-}";
  }
  else if(particle.CompareTo("EtaPrim") == 0){
    particlelegend = "#eta^{l}";
  }
  else if(particle.CompareTo("Omega") == 0){
    particlelegend = "#omega";
  }
  else{
    cout << particle.Data() << " unknown! Aborting";
    return;
  }
  
  TString ptcut="";
  TString pthard_legend="";
  TString mbias="";
  
  
  TString detector="PCM";
  TString binning="";
  // histograms for MB
  TH2F* fHPt_Y_ParticleMB;
  TH2F* fHPt_Y_ParticleAccMB;
  TH2F* fHPt_Alpha_ParticleAccMB;
  
  TH1D* fHPt_ParticleMB;
  TH1D* fHPt_ParticleMBUB;
  TH1D* fHPt_ParticleAccMB;
  TH1D* fHPtFrac_ParticleMB;

  TH1D* fHNEventsMB;
  TH1D* fHXSectionMB;
  
  TString nameMainDir = "GammaPureMC";
  
   // open MB stuff
   sprintf(name,"%sMB%s.root",filemask.Data(),tune.Data());
   TFile f(name);
   // get tdirectoryfile
   TDirectoryFile* directoryAll = (TDirectoryFile*)f.Get(nameMainDir.Data());
  
  // get top list
    TList* HistosAll =(TList*)directoryAll->Get(nameMainDir.Data())->Clone();
    if(HistosAll == NULL){
        cout<<"ERROR: List not found"<<endl;
        return;
      }
    fHNEventsMB = (TH1D*) HistosAll->FindObject("NEvents")->Clone("fHNEventsMB");
    fHXSectionMB = (TH1D*) HistosAll->FindObject("XSection")->Clone("fHXSectionMB");
  
    // get particle histograms
    fHPt_Y_ParticleMB = (TH2F*) HistosAll->FindObject(Form("Pt_Y_%s",particle.Data()))->Clone(Form("Pt_Y_%s_MB",particle.Data()));
    fHPt_Y_ParticleAccMB = (TH2F*) HistosAll->FindObject(Form("Pt_Y_%sGG%sAcc",particle.Data(),"PCM"))->Clone(Form("Pt_Y_%sGG%sAcc_MB",particle.Data(),detector.Data()));
    fHPt_Alpha_ParticleAccMB = (TH2F*) HistosAll->FindObject(Form("Pt_Alpha_%sGG%sAcc",particle.Data(),detector.Data()))->Clone(Form("Pt_Alpha_%sGG%sAcc_MB",particle.Data(),detector.Data()));
  
    //scale
    fHPt_Y_ParticleMB->Sumw2();
    fHPt_Y_ParticleAccMB->Sumw2();
   fHPt_Alpha_ParticleAccMB->Sumw2();
  
  
    cout << "these are " << fHNEventsMB->GetBinContent(1) << " events with a cross section of " << fHXSectionMB->GetMean() << " mb" << endl;
  
  
    fHPt_Y_ParticleMB->Scale(1./fHNEventsMB->GetBinContent(1));
    fHPt_Y_ParticleAccMB->Scale(1./fHNEventsMB->GetBinContent(1));
    fHPt_Alpha_ParticleAccMB->Scale(1./fHNEventsMB->GetBinContent(1));
  
    // now projections on pT axis
    fHPt_ParticleMBUB = (TH1D*)fHPt_Y_ParticleMB->ProjectionX(Form("hPt_%s_MB",particle.Data()),1,160)->Clone(Form("hPt_%s_MB",particle.Data()));
    fHPt_ParticleMB  = (TH1D*) fHPt_ParticleMBUB->Rebin(nbinsfinal,"hPt_ParticleMB",fBinsPi02760GeVPtMerged);
    fHPt_ParticleAccMB = (TH1D*)fHPt_Y_ParticleAccMB->ProjectionX(Form("hPt_%sAcc_MB",particle.Data()),1,160)->Clone(Form("hPt_%sAcc_MB",particle.Data()));
  
    fHPt_ParticleMB->SetName(Form("hPt%s%s%s",particle.Data(),generator.Data(),tune.Data()));
    fHPt_ParticleMB->SetTitle(Form("pp, %s, %s, %s %s",energy.Data(),particle.Data(),generator.Data(),tune.Data()));
    fHPt_ParticleMB->Sumw2();
    fHPt_ParticleAccMB->Sumw2();
    //fHPt_ParticleMB->Rebin(nreb);
    fHPt_ParticleAccMB->Rebin(nreb);
  
  // divide MB by bin width
  for(int ibin=1;ibin<fHPt_ParticleMB->GetNbinsX();ibin++){
    double dpt = fHPt_ParticleMB->GetBinWidth(ibin);
    if(dpt == 0){
      continue;
    }
    fHPt_ParticleMB->SetBinContent(ibin,fHPt_ParticleMB->GetBinContent(ibin)/dpt);
    fHPt_ParticleMB->SetBinError(ibin,fHPt_ParticleMB->GetBinError(ibin)/dpt);
  }
  
  // divide MB by pT
  for(int ibin=1;ibin<fHPt_ParticleMB->GetNbinsX();ibin++){
    double binpt = fHPt_ParticleMB->GetBinCenter(ibin);
    if(binpt == 0){
      continue;
    }
    fHPt_ParticleMB->SetBinContent(ibin,fHPt_ParticleMB->GetBinContent(ibin)/binpt);
    fHPt_ParticleMB->SetBinError(ibin,fHPt_ParticleMB->GetBinError(ibin)/binpt);
  }

  // scale with 1/2pi y
  fHPt_ParticleMB->Scale(1./(2.*3.14159265*1.6));
  
  fHPt_ParticleMB->Print("all");
    // now make plots like no other
    
    // now make plots like no other
  
    //***************************************************************************************************************
    //*************************************Plotting *****************************************************************
    //***************************************************************************************************************
  
    // input
    TCanvas* canvasInputScaled = new TCanvas("canvasInputScaled","",0,0,1350,1350);  // gives the page size
    DrawGammaCanvasSettings( canvasInputScaled, 0.13, 0.015, 0.015, 0.07);
    canvasInputScaled->SetLogy();
  
    Float_t maximumPi0Scaled = FindLargestEntryIn1D(fHPt_ParticleMB)*10;
    Float_t minimumPi0Scaled = FindSmallestEntryIn1D(fHPt_ParticleMB)*0.1;
  
    TH2F * histo2DInputScaledPi0;
    histo2DInputScaledPi0 = new TH2F("histo2DInputScaledPi0","histo2DInputScaledPi0",1000,0., 50,10000,minimumPi0Scaled,maximumPi0Scaled);
    SetStyleHistoTH2ForGraphs(histo2DInputScaledPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",
                                                           0.028,0.035, 0.028,0.035, 0.8,1.6);
  
    histo2DInputScaledPi0->DrawCopy();
  
    //    TLegend* legendScaled = GetAndSetLegend2(0.25, 0.99-(1.15*2/2*0.040), 0.95, 0.99,25);
    //    legendScaled->SetMargin(0.12);
    //    legendScaled->SetNColumns(2);
  
    DrawGammaSetMarker(fHPt_ParticleMB, markerBins[0], 1., colorBins[0], colorBins[0]);
    fHPt_ParticleMB->DrawCopy("e1,same");
  
    //legendScaled->AddEntry(fHPt_ParticleMB,"MB","p");
  
    //legendScaled->Draw();
    //    labelMCName->Draw();
    TLatex *labelSimulation = new TLatex(0.38,0.99-(1.15*(1/2+2)*0.04),Form("%s, %s, %s",particlelegend.Data(),generator.Data(),tune.Data()));
    SetStyleTLatex( labelSimulation, 32,4);
    labelSimulation->SetTextFont(43);
    labelSimulation->Draw();
  
    canvasInputScaled->Update();
    canvasInputScaled->SaveAs(Form("%s/%s_MC_InputUnscaled_%s%s_%s_%s%s.%s",outputDir.Data(),particle.Data(),generator.Data(),tune.Data(),binning.Data(),mbias.Data(),ptcut.Data(),suffix.Data()));
  
    // write out
    gSystem->Exec("mkdir -p ../ExternalInput/Theory/Pythia/");
    TFile* fout = new TFile(Form("../ExternalInput/Theory/Pythia/hist_%s_%s_%s.root",generator.Data(),tune.Data(),energy.Data()),"RECREATE");
    fHPt_ParticleMB->Write();
    fout->Close();
  
  
}
