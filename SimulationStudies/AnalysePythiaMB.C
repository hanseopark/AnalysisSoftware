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
//******************* Main function ************************************************
//**********************************************************************************

void AnalysePythiaMB(TString fileName="Legotrain_vAN-20150825-Pythia8/PythiaAnalysisResults",  // should read something like "$PATHTODIRECTORY/PythiaAnalysisResultsBin" // use MC_gen output here!
                         TString generator="Pythia",  // e.g. Pythia
                         TString tune="Monash",    // e.g. Monash
                         TString energy="2760GeV"){
  
  TH1::AddDirectory(kFALSE);
  
//   const Int_t nbinsfinal                    = 93;
//   Double_t fBinsFine[nbinsfinal+1]          = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
//                                                 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
//                                                 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
//                                                 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 
//                                                 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8,
//                                                 7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 9.4, 9.8, 10.2, 10.6,
//                                                 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5,
//                                                 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 
//                                                 26.0, 27.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0,
//                                                 44.0, 46.0, 48.0, 50.0};
//   const Int_t nbinsfinal                    = 84;
//   Double_t fBinsFine[nbinsfinal+1]          = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
//                                                 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
//                                                 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
//                                                 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 
//                                                 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8,
//                                                 7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 9.4, 9.8, 10.2, 10.6,
//                                                 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5,
//                                                 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 24.0, 26.0, 28.0, 
//                                                 32.0, 36.0, 40.0, 45.0, 50.0};

  const Int_t nbinsfinal                    = 76;
  Double_t fBinsFine[nbinsfinal+1]          = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
                                                1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                                2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                                                3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 
                                                5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8,
                                                7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 9.4, 9.8, 10.2, 10.6,
                                                11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 
                                                26.0, 28.0, 32.0, 36.0, 40.0, 45.0, 50.0};
                                                
  //***************************************************************************************************************
  //************************************ Layouting preparations & general setup ***********************************
  //***************************************************************************************************************
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  
  StyleSettingsThesis();
  SetPlotStyle();
  
  Color_t colorBins[12]		= {kBlack, kRed+2, kBlue+2, kGreen+2, kCyan+2, kViolet, kMagenta+2, kGray+1, kRed-2, kBlue-2, kViolet+2, kCyan-2};
  Color_t colorBinsShade[12]	= {kGray+1, kRed-6, kBlue-6, kGreen-8, kCyan-6, kViolet-8, kMagenta-8, kGray, kRed-8, kBlue-8, kViolet-6, kCyan-8};
  Marker_t markerBins[12]		= {20, 21, 33, 34, 29, 24, 25, 27, 28, 30, 20, 21};
  
  
  TString outputDir             ="./plots/";
  TString suffix                = "eps";
  
  gSystem->Exec("mkdir -p "+outputDir);
  
  //fDetectionProcess 			= ReturnFullTextReconstructionProcess(mode);
  TString particle[6]           = {"Pi0", "Eta", "EtaPrim", "Omega", "PiPl", "PiMi"};  
  TString particlelegend[6]     = {"#pi^{0}", "#eta", "#eta^{l}", "#omega", "#pi^{+}", "#pi^{-}"}; 
  // translate particle string to plotting string
  
  TString ptcut                 ="";
  TString pthard_legend         ="";
  TString mbias                 ="";
  
  TString binning="";
  // histograms for MB
  TH2F* fHPt_Y_ParticleMB[6];  
  TH1D* fHPt_ParticleMB[6];
  TH1D* fHPt_ParticleMBUB[6];
  TH1D* fHPt_ParticleMB_XSec[6];

  TH1D* fHNEventsMB;
  TH1D* fHXSectionMB;
  
  TString nameMainDir           = "GammaPureMC";
  
   // open MB stuff
   sprintf(name,"%s",fileName.Data());
   TFile f(name);
   // get tdirectoryfile
   TDirectoryFile* directoryAll = (TDirectoryFile*)f.Get(nameMainDir.Data());
  
  // get top list
    TList* HistosAll =(TList*)directoryAll->Get(nameMainDir.Data())->Clone();
    if(HistosAll == NULL){
        cout<<"ERROR: List not found"<<endl;
        return;
      }
    fHNEventsMB         = (TH1D*) HistosAll->FindObject("NEvents")->Clone("fHNEventsMB");
    fHXSectionMB        = (TH1D*) HistosAll->FindObject("XSection")->Clone("fHXSectionMB");
    cout << "these are " << fHNEventsMB->GetBinContent(1) << " events with a cross section of " << fHXSectionMB->GetMean() << " mb" << endl;
    cout << "total inel. cross section #sigma_{inel} = " << fHXSectionMB->GetMean() * 1e9 << " pb" << endl;
    Int_t eventInMio    = (Int_t)(fHNEventsMB->GetBinContent(1)/1e6);
    
    for (Int_t i = 0; i < 6; i++){
        // get particle histograms
        fHPt_Y_ParticleMB[i]    = (TH2F*) HistosAll->FindObject(Form("Pt_Y_%s",particle[i].Data()))->Clone(Form("Pt_Y_%s_MB",particle[i].Data()));
    
        //scale
        fHPt_Y_ParticleMB[i]->Sumw2();
        fHPt_Y_ParticleMB[i]->Scale(1./fHNEventsMB->GetBinContent(1));
    
        // now projections on pT axis
        fHPt_ParticleMBUB[i]    = (TH1D*)fHPt_Y_ParticleMB[i]->ProjectionX(Form("hPt_%s_MB",particle[i].Data()),1,160)->Clone(Form("hPt_%s_MB",particle[i].Data()));
        fHPt_ParticleMB[i]      = (TH1D*)fHPt_ParticleMBUB[i]->Rebin(nbinsfinal,Form("hPt_%s_MB", particle[i].Data()),fBinsFine);
    
        fHPt_ParticleMB[i]->SetName(Form("hPt%s%s%s",particle[i].Data(),generator.Data(),tune.Data()));
        fHPt_ParticleMB[i]->SetTitle(Form("pp, %s, %s, %s %s",energy.Data(),particle[i].Data(),generator.Data(),tune.Data()));
        fHPt_ParticleMB[i]->Sumw2();
    
        // divide MB by bin width
        for(int ibin=1;ibin<fHPt_ParticleMB[i]->GetNbinsX()+1;ibin++){
            double dpt = fHPt_ParticleMB[i]->GetBinWidth(ibin);
            if(dpt == 0){
                continue;
            }
            fHPt_ParticleMB[i]->SetBinContent(ibin,fHPt_ParticleMB[i]->GetBinContent(ibin)/dpt);
            fHPt_ParticleMB[i]->SetBinError(ibin,fHPt_ParticleMB[i]->GetBinError(ibin)/dpt);
        }
        
        // divide MB by pT
        for(int ibin=1;ibin<fHPt_ParticleMB[i]->GetNbinsX()+1;ibin++){
            double binpt = fHPt_ParticleMB[i]->GetBinCenter(ibin);
            if(binpt == 0){
            continue;
            }
            fHPt_ParticleMB[i]->SetBinContent(ibin,fHPt_ParticleMB[i]->GetBinContent(ibin)/binpt);
            fHPt_ParticleMB[i]->SetBinError(ibin,fHPt_ParticleMB[i]->GetBinError(ibin)/binpt);
        }

        // scale with 1/2pi y
        fHPt_ParticleMB[i]->Scale(1./(2.*TMath::Pi()*1.6));
        
        // print total cross section
        
        // make histogram for x section
        fHPt_ParticleMB_XSec[i] = (TH1D*) fHPt_ParticleMB[i]->Clone(Form("hPt_%s_MB_XSec",particle[i].Data()));
        fHPt_ParticleMB_XSec[i]->Scale(fHXSectionMB->GetMean()); // get xsec in mb from pythia
        fHPt_ParticleMB_XSec[i]->Scale(1e9);                     // pb!
    }  
    
//   fHPt_ParticleMB_XSec->Print("all");
  
    //***************************************************************************************************************
    //*************************************Plotting *****************************************************************
    //***************************************************************************************************************
  
    // input
    TCanvas* canvasInputScaled = new TCanvas("canvasInputScaled","",0,0,1350,1350);  // gives the page size
    DrawGammaCanvasSettings( canvasInputScaled, 0.13, 0.015, 0.015, 0.07);
    canvasInputScaled->SetLogy();
  
    Float_t maximumPi0Scaled = FindLargestEntryIn1D(fHPt_ParticleMB[0])*10;
    Float_t minimumPi0Scaled = FindSmallestEntryIn1D(fHPt_ParticleMB[2])*0.1;
  
    TH2F * histo2DInputScaledPi0;
    histo2DInputScaledPi0 = new TH2F("histo2DInputScaledPi0","histo2DInputScaledPi0",1000,0., 50,10000,minimumPi0Scaled,maximumPi0Scaled);
    SetStyleHistoTH2ForGraphs(histo2DInputScaledPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",
                                                           0.028,0.035, 0.028,0.035, 0.8,1.6);
  
    histo2DInputScaledPi0->DrawCopy();
  
    TLegend* legendScaled = GetAndSetLegend2(0.85, 0.95-(1.15*6/1.5*0.040), 0.95, 0.95,35);
    legendScaled->SetMargin(0.12);
    //    legendScaled->SetNColumns(2);
  
    for (Int_t i = 0; i< 6; i++){
        DrawGammaSetMarker(fHPt_ParticleMB[i], markerBins[i], 1., colorBins[i], colorBins[i]);
        DrawGammaSetMarker(fHPt_ParticleMB_XSec[i], markerBins[i], 1., colorBins[i], colorBins[i]);
        fHPt_ParticleMB[i]->DrawCopy("e1,same");
        legendScaled->AddEntry(fHPt_ParticleMB[i],particlelegend[i],"p");
    }
    legendScaled->Draw();
    
    TLatex *labelSimulation = new TLatex(0.38,0.99-(1.15*(1/2+2)*0.04),Form("%s, %s",generator.Data(),tune.Data()));
    SetStyleTLatex( labelSimulation, 32,4);
    labelSimulation->SetTextFont(43);
    labelSimulation->Draw();

    canvasInputScaled->Update();
    canvasInputScaled->SaveAs(Form("%s/%s_MC_InputUnscaled_%s_%dMio.%s",outputDir.Data(),generator.Data(),tune.Data(),eventInMio,suffix.Data()));
  
    // write out
//     gSystem->Exec("mkdir -p ../ExternalInput/Theory/Pythia/");
    TFile* fout = new TFile(Form("%s_%s_%s_%dMio.root",generator.Data(),tune.Data(),energy.Data(),eventInMio),"RECREATE");
    fHNEventsMB->Write();
    fHXSectionMB->Write();
    for (Int_t i = 0; i< 6; i++){
        fHPt_ParticleMB[i]->Write();
        fHPt_ParticleMB_XSec[i]->Write();
    }    
    fout->Close();
  
}
