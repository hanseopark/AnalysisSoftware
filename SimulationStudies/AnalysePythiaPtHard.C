// provided to run with output from PWGGA/GammaConv/macros/AddTask_GammaPureMC.C
// takes pT hard binned pythia MC output, calculates weights, adds spectra from different bins and plots particle spectra and acceptance

// run for example as
// root -b -l -q 'AnalysePythiaPtHard.C+("grid/Legotrain_vAN-20150825-Pythia6/PythiaAnalysisResults",5,"Pythia8","Perugia2011","Pi0",1,1,1,1,0.,4)'

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

#include "AnalysePythiaPtHard.h"

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

void AnalysePythiaPtHard(TString filemask="",  // should read something like "$PATHTODIRECTORY/PythiaAnalysisResultsBin"
                         const Int_t nbins=0,    // run MB if nbins is 0!
                         TString generator="",
                         TString tune="",
                         TString particle="Pi0",
                         Bool_t bscale = 0,
                         Bool_t bexternalinput = 1,
                         Bool_t pthard_mode = 1, // 1 for narrow binning, other for wide binning (2 is wide full)
                         Bool_t bDoMb = 1,
                         Double_t CutPtH = 0., // if 0.: no cut, if 0.5: cut in middle of bin, if any number > 1: cut at number*pthard_max
                         Int_t mode = 4){
  
  TH1::AddDirectory(kFALSE);
  
  //***************************************************************************************************************
  //************************************ Layouting preparations & general setup ***********************************
  //***************************************************************************************************************
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  
  StyleSettingsThesis();
  SetPlotStyle();
  
  const Int_t nreb = 1;
  
  const Int_t MaxNumberOfFiles = 12;
  
  Color_t colorBins[12]		= {kBlack, kRed+2, kBlue+2, kGreen+2, kCyan+2, kViolet, kMagenta+2, kGray+1, kRed-2, kBlue-2, kViolet+2, kCyan-2};
  Color_t colorBinsShade[12]	= {kGray+1, kRed-6, kBlue-6, kGreen-8, kCyan-6, kViolet-8, kMagenta-8, kGray, kRed-8, kBlue-8, kViolet-6, kCyan-8};
  Marker_t markerBins[12]		= {20, 21, 33, 34, 29, 24, 25, 27, 28, 30, 20, 21};
  
  // some numbers
  Float_t minPtHardSmall[MaxNumberOfFiles] = {5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85};
  Float_t maxPtHardSmall[MaxNumberOfFiles] = {7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 1000};
  
  Float_t minPtHardWide[MaxNumberOfFiles] = {5, 11, 21, 36, 57, 70, 128, 136, 145, 157, 170, 185};
  Float_t maxPtHardWide[MaxNumberOfFiles] = {11, 21, 36, 57, 70, 85, 136, 145, 157, 170, 185, 1000};
  
  Float_t minPtHard[MaxNumberOfFiles];
  Float_t maxPtHard[MaxNumberOfFiles];
  if(pthard_mode == 1){
    for(int i=0;i<MaxNumberOfFiles;i++){
      minPtHard[i] = minPtHardSmall[i];
      maxPtHard[i] = maxPtHardSmall[i];
    }
  }
  // for the wide
  else{
    for(int i=0;i<MaxNumberOfFiles;i++){
      minPtHard[i] = minPtHardWide[i];
      maxPtHard[i] = maxPtHardWide[i];
    }
  }
  
  Float_t floatLocationRightUp2D[4] = {0.45,0.95,0.035, 0.02};
  Float_t floatLocationLeftDown2D[4] = {0.15,0.25,0.035, 0.02};
  Float_t floatLocationRightDown2D[4] = {0.45,0.25,0.035, 0.02};
  
  TString outputDir ="./plots/";
  TString suffix = "eps";
  
  gSystem->Exec("mkdir -p "+outputDir);
  
  fDetectionProcess 			= ReturnFullTextReconstructionProcess(mode);
  
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
  
  TString detector="";
  
  if(mode == 0){
    detector="PCM";
  }
  else if(mode == 2){
    detector="PCMEMC";
  }
  else if(mode == 3){
    detector="PCMPHO";
  }
  else if(mode == 4){
    detector="EMC";
  }
  else if(mode == 5){
    detector="PHO";
  }
  else{
    cout << "wrong mode! Aborting" << endl;
    return;
  }
  
  TString binning="";
  if(pthard_mode == 1){
    binning="narrow";
  }
  else if(pthard_mode==2){
    binning="widefull";
  }
  else{
    binning="wide";
  }
  
  TString mbias="";
  if(bDoMb==1){
    mbias="withMB";
  }
  else{
    mbias="noMB";
  }
  
  TString ptcut="";
  TString pthard_legend="";
  if(CutPtH == 0.5){
    ptcut="_CutPtHCenter";
    pthard_legend=Form(", #it{p}^{%s}_{T} <  (#it{p}^{hard}_{T, min} + #it{p}^{hard}_{T, max})/2",particlelegend.Data());
  }
  else if(CutPtH>=1.){
    ptcut=Form("_CutPtH%.2fPtHMax",CutPtH);
    pthard_legend=Form(", #it{p}^{%s}_{T} <  %.2f #it{p}^{hard}_{T, max}",particlelegend.Data(),CutPtH);
  }
    
  // histograms for MB comparison
  TH2F* fHPt_Y_ParticleMB;
  TH2F* fHPt_Y_ParticleAccMB;
  TH2F* fHPt_Alpha_ParticleAccMB;
  
  TH1D* fHPt_ParticleMB;
  TH1D* fHPt_ParticleAccMB;
  TH1D* fHPtFrac_ParticleMB;

  TH1D* fHNEventsMB;
  TH1D* fHXSectionMB;
  
  // histograms for all pT hard bins
  // control histograms
  TH1D* fHNEvents[nbins];
  TH1D* fHXSection[nbins];
  
  // original 2D histograms
  TH2F* fHPt_Y_Particle[nbins]; // particle pt vs y
  TH2F* fHPt_Y_ParticleAcc[nbins]; // accepted particle pt vs y
  TH2F* fHPt_Alpha_ParticleAcc[nbins]; // accepted particle pt vs alpha
  // pT projections
  TH1D* fHPt_Particle[nbins]; // particle pt
  TH1D* fHPt_ParticleAcc[nbins]; // accepted particle pt
  // fraction to total
  TH1D* fHPtFrac_Particle[nbins]; // particle pt
  TH1D* fHPtFrac_ParticleAcc[nbins]; // accepted particle pt
  
  
  // inclusive (all bins)
  TH1D* fHPt_ParticleAll; // sum, particle pt
  TH1D* fHPt_ParticleAccAll; // sum, accepted particle pt

  TH1D* fHPt_ParticleRelStatAll; // sum, particle pt
  TH1D* fHPt_ParticleRelStatAccAll; // sum, accepted particle pt

  TH2F* fHPt_Alpha_ParticleAccAll;
  //TH1D* fHPt_ParticleAcceptanceAll; // acceptance vs pt
  
  // number of trials
  Double_t number_of_trials[MaxNumberOfFiles];
  Double_t number_of_trials_small[MaxNumberOfFiles];
  Double_t number_of_trials_wide[MaxNumberOfFiles];
  
  TString nameMainDir = "GammaPureMC";
  
  if(nbins < 0){
    cout << "does not run with negative number of bins!";
    return;
  }
  
  else if(nbins > 0){
    
    // start by setting the variable to 1
    for(int i = 0; i < MaxNumberOfFiles; i++){
      number_of_trials[i] = 1.;
    }
    
    cout << __LINE__ << endl;
    
    // read trials from external file, for small bins
    if(bexternalinput){
      ifstream in1(Form("./ExternalInput/scale%s",tune.Data()));
      Int_t i1= 0;
      while(!in1.eof() && i1 < MaxNumberOfFiles){
        in1 >> number_of_trials_small[i1];
        i1++;
      }
      for(int i=i1;i<MaxNumberOfFiles;i++){
        number_of_trials_small[i] = 1.;
      }
      in1.close();
      
      // read trials from external file, for wide bins
      ifstream in2(Form("./ExternalInput/scale%sWide",tune.Data()));
      Int_t i2 = 0;
      while(!in2.eof() && i2 < MaxNumberOfFiles){
        in2 >> number_of_trials_wide[i2];
        i2++;
      }
      for(int i=i2;i<MaxNumberOfFiles;i++){
        number_of_trials_wide[i] = 1.;
      }
      in2.close();
      
      // now fill used variable accordingly!
      if(pthard_mode == 1){
        for(int i=0;i<MaxNumberOfFiles;i++){
          number_of_trials[i] = number_of_trials_small[i];
        }
      }
      else{
        for(int i=0;i<MaxNumberOfFiles;i++){
          number_of_trials[i] = number_of_trials_wide[i];
        }
      }
      
    }
    
    // look at number_of_trials
    for(int i=0;i<nbins;i++){
      cout << "number of trials bin " << i << " " << number_of_trials[i] << endl;
    }
    
    
    // open MB stuff
    sprintf(name,"%sMB%s.root",filemask.Data(),tune.Data());
    TFile f(name);
    // get tdirectoryfile
    TDirectoryFile* directoryAll = (TDirectoryFile*)f.Get(nameMainDir.Data());
    
    // get top list
    TList* HistosAll =(TList*)directoryAll->Get(nameMainDir.Data())->Clone();
    if(HistosAll == NULL){
      cout<<"ERROR: List not found"<<endl;
      bDoMb = 0;
    }
    if(bDoMb){
      fHNEventsMB = (TH1D*) HistosAll->FindObject("NEvents")->Clone("fHNEventsMB");
      fHXSectionMB = (TH1D*) HistosAll->FindObject("XSection")->Clone("fHXSectionMB");
      
      // get particle histograms
      fHPt_Y_ParticleMB = (TH2F*) HistosAll->FindObject(Form("Pt_Y_%s",particle.Data()))->Clone(Form("Pt_Y_%s_MB",particle.Data()));
      fHPt_Y_ParticleAccMB = (TH2F*) HistosAll->FindObject(Form("Pt_Y_%sGG%sAcc",particle.Data(),detector.Data()))->Clone(Form("Pt_Y_%sGG%sAcc_MB",particle.Data(),detector.Data()));
      fHPt_Alpha_ParticleAccMB = (TH2F*) HistosAll->FindObject(Form("Pt_Alpha_%sGG%sAcc",particle.Data(),detector.Data()))->Clone(Form("Pt_Alpha_%sGG%sAcc_MB",particle.Data(),detector.Data()));
      
      //scale
      // ATTENTION! ATTENTION! ATTENTION! ATTENTION! ATTENTION!
      // LOOK OUT, THERE IS A FUDGE FACTOR OF 0.75 RIGHT NOW!!!
      // ATTENTION! ATTENTION! ATTENTION! ATTENTION! ATTENTION!
      fHPt_Y_ParticleMB->Sumw2();
      fHPt_Y_ParticleAccMB->Sumw2();
      fHPt_Alpha_ParticleAccMB->Sumw2();
      
      fHPt_Y_ParticleMB->Scale(0.75*fHXSectionMB->GetMean()/fHNEventsMB->GetBinContent(1));
      fHPt_Y_ParticleAccMB->Scale(0.75*fHXSectionMB->GetMean()/fHNEventsMB->GetBinContent(1));
      fHPt_Alpha_ParticleAccMB->Scale(0.75*fHXSectionMB->GetMean()/fHNEventsMB->GetBinContent(1));
      
      // now projections on pT axis
      fHPt_ParticleMB = (TH1D*)fHPt_Y_ParticleMB->ProjectionX(Form("hPt_%s_MB",particle.Data()),1,160)->Clone(Form("hPt_%s_MB",particle.Data()));
      fHPt_ParticleAccMB = (TH1D*)fHPt_Y_ParticleAccMB->ProjectionX(Form("hPt_%sAcc_MB",particle.Data()),1,160)->Clone(Form("hPt_%sAcc_MB",particle.Data()));
      fHPt_ParticleMB->Sumw2();
      fHPt_ParticleAccMB->Sumw2();
      fHPt_ParticleMB->Rebin(nreb);
      fHPt_ParticleAccMB->Rebin(nreb);
    }
    
    //***************************************************************************************************************
    //**************************************** loop through pT hard bins ********************************************
    //***************************************************************************************************************
    
    for(int ibin = 0; ibin < nbins; ibin++){
      sprintf(name,"%sBin%d%s.root",filemask.Data(),ibin+1,tune.Data());
      cout << "opening file " << name << endl;
      TFile f(name);
      
      // get tdirectoryfile
      TDirectoryFile* directoryAll = (TDirectoryFile*)f.Get(nameMainDir.Data());
      
      // get top list
      TList* HistosAll =(TList*)directoryAll->Get(nameMainDir.Data())->Clone();
      if(HistosAll == NULL){
        cout<<"ERROR: List not found"<<endl;
        return;
      }
      
      // now get the necessary histograms
      fHNEvents[ibin] = (TH1D*) HistosAll->FindObject("NEvents")->Clone(Form("fHNEvents%d",ibin+1));
      fHXSection[ibin] = (TH1D*) HistosAll->FindObject("XSection")->Clone(Form("fHXSections%d",ibin+1));
      
      Double_t scalefactor = 1.;
      if(bscale){
        if(bexternalinput){
          Double_t xsection = fHXSection[ibin]->GetMean();
          Double_t nevents = fHNEvents[ibin]->GetBinContent(1);
          
          scalefactor = xsection/number_of_trials[ibin]*1/nevents;
        }
        else{
          // calculate scaling of this bin
          Double_t xsection = fHXSection[ibin]->GetMean();
          Double_t ntrials = fHNEvents[ibin]->GetBinContent(3); // nevents cancels in the end
          
          scalefactor = xsection/ntrials;
        }
        cout << "scaling bin " << ibin+1 << " by " << scalefactor << endl;
      }
      
      // get particle histograms
      fHPt_Y_Particle[ibin] = (TH2F*) HistosAll->FindObject(Form("Pt_Y_%s",particle.Data()))->Clone(Form("Pt_Y_%s_Bin%d",particle.Data(),ibin+1));
      fHPt_Y_ParticleAcc[ibin] = (TH2F*) HistosAll->FindObject(Form("Pt_Y_%sGG%sAcc",particle.Data(),detector.Data()))->Clone(Form("Pt_Y_%sGG%sAcc_Bin%d",particle.Data(),detector.Data(),ibin+1));
      fHPt_Alpha_ParticleAcc[ibin] = (TH2F*) HistosAll->FindObject(Form("Pt_Alpha_%sGG%sAcc",particle.Data(),detector.Data()))->Clone(Form("Pt_Alpha_%sGG%sAcc_Bin%d",particle.Data(),detector.Data(),ibin+1));
      
      fHPt_Y_Particle[ibin]->Sumw2();
      fHPt_Y_ParticleAcc[ibin]->Sumw2();
      fHPt_Alpha_ParticleAcc[ibin]->Sumw2();
      
      // scale them
      fHPt_Y_Particle[ibin]->Scale(scalefactor);
      fHPt_Y_ParticleAcc[ibin]->Scale(scalefactor);
      fHPt_Alpha_ParticleAcc[ibin]->Scale(scalefactor);
      
      
      // now projections on pT axis
      fHPt_Particle[ibin] = (TH1D*)fHPt_Y_Particle[ibin]->ProjectionX(Form("hPt_%s_Bin%d",particle.Data(),ibin+1),1,160)->Clone(Form("hPt_%s_Bin%d",particle.Data(),ibin+1));
      fHPt_ParticleAcc[ibin] = (TH1D*)fHPt_Y_ParticleAcc[ibin]->ProjectionX(Form("hPt_%sAcc_Bin%d",particle.Data(),ibin+1),1,160)->Clone(Form("hPt_%sAcc_Bin%d",particle.Data(),ibin+1));
      fHPt_Particle[ibin]->Sumw2();
      fHPt_ParticleAcc[ibin]->Sumw2();

      fHPt_Particle[ibin]->Rebin(nreb);
      fHPt_ParticleAcc[ibin]->Rebin(nreb);

      // here we can cut away the high pT tails of each spectrum
      if(CutPtH>0){
        Double_t pthard_cut = 0.;
        if(CutPtH == 0.5){
          pthard_cut = (minPtHard[ibin]+maxPtHard[ibin])/2;
        }
        else if(CutPtH>=1.){
          pthard_cut = CutPtH*maxPtHard[ibin];
        }
        else{
          cout << "wrong value of pthard cut. Aborting" << endl;
          return;
        }
        for(int hbin=1;hbin<fHPt_Particle[ibin]->GetNbinsX();hbin++){
          Double_t pt = fHPt_Particle[ibin]->GetBinCenter(hbin);
          if(pt > pthard_cut){
            fHPt_Particle[ibin]->SetBinContent(hbin,0.);
            fHPt_Particle[ibin]->SetBinError(hbin,0.);
          }
        }
        for(int hbin=1;hbin<fHPt_ParticleAcc[ibin]->GetNbinsX();hbin++){
          Double_t pt = fHPt_ParticleAcc[ibin]->GetBinCenter(hbin);
          if(pt > pthard_cut){
            fHPt_ParticleAcc[ibin]->SetBinContent(hbin,0.);
            fHPt_ParticleAcc[ibin]->SetBinError(hbin,0.);
          }
        }
      }
      // and create composed histogram
      if(ibin == 0){
        fHPt_ParticleAll = (TH1D*)fHPt_Particle[ibin]->Clone(Form("hPt_%s",particle.Data()));
        fHPt_ParticleAccAll = (TH1D*)fHPt_ParticleAcc[ibin]->Clone(Form("hPt_%sAcc",particle.Data()));
        fHPt_ParticleAll->Sumw2();
        fHPt_ParticleAccAll->Sumw2();
        fHPt_Alpha_ParticleAccAll = (TH2F*)fHPt_Alpha_ParticleAcc[ibin]->Clone(Form("hPt_Alpha_%s",particle.Data()));
        fHPt_Alpha_ParticleAccAll->Sumw2();
      }
      else{
        fHPt_ParticleAll->Add(fHPt_Particle[ibin]);
        fHPt_ParticleAccAll->Add(fHPt_ParticleAcc[ibin]);
        fHPt_Alpha_ParticleAccAll->Add(fHPt_Alpha_ParticleAcc[ibin]);
      }
      cout << " done" << endl;
    }
    
    // now we can make ratios
    if(bDoMb){
      fHPtFrac_ParticleMB = (TH1D*)fHPt_ParticleMB->Clone(Form("hPtFrac_%MB",particle.Data()));
      fHPtFrac_ParticleMB->Divide(fHPt_ParticleAll);
    }
    for(int ibin=0;ibin<nbins;ibin++){
      fHPtFrac_Particle[ibin] = (TH1D*)fHPt_Particle[ibin]->Clone(Form("hPtFrac_%s_Bin%d",particle.Data(),ibin+1));
      fHPtFrac_ParticleAcc[ibin] = (TH1D*)fHPt_ParticleAcc[ibin]->Clone(Form("hPtFrac_%sAcc_Bin%d",particle.Data(),ibin+1));
      
      fHPtFrac_Particle[ibin]->Divide(fHPt_ParticleAll);
      fHPtFrac_ParticleAcc[ibin]->Divide(fHPt_ParticleAccAll);
    }
    
    // relative statistical errors
    fHPt_ParticleRelStatAll = (TH1D*)fHPt_ParticleAll->Clone(Form("hPt_%sRelStat",particle.Data()));
    fHPt_ParticleRelStatAccAll = (TH1D*)fHPt_ParticleAccAll->Clone(Form("hPt_%sAccRelStat",particle.Data()));
    fHPt_ParticleRelStatAll->Sumw2();
    fHPt_ParticleRelStatAccAll->Sumw2();

    for(int hbin=1;hbin<fHPt_ParticleRelStatAll->GetNbinsX();hbin++){
      if(fHPt_ParticleRelStatAll->GetBinContent(hbin) != 0){
        fHPt_ParticleRelStatAll->SetBinContent(hbin,fHPt_ParticleRelStatAll->GetBinError(hbin)/fHPt_ParticleRelStatAll->GetBinContent(hbin));
      }
      else{
        fHPt_ParticleRelStatAll->SetBinContent(hbin,0.);
      }
      fHPt_ParticleRelStatAll->SetBinError(hbin,0);
    }
    for(int hbin=1;hbin<fHPt_ParticleRelStatAccAll->GetNbinsX();hbin++){
      if(fHPt_ParticleRelStatAccAll->GetBinContent(hbin) != 0){
        fHPt_ParticleRelStatAccAll->SetBinContent(hbin,fHPt_ParticleRelStatAccAll->GetBinError(hbin)/fHPt_ParticleRelStatAccAll->GetBinContent(hbin));
      }
      else{
        fHPt_ParticleRelStatAccAll->SetBinContent(hbin,0.);
      }
      fHPt_ParticleRelStatAccAll->SetBinError(hbin,0);
    }
    
    // now make plots like no other
    
    //***************************************************************************************************************
    //*************************************Plotting *****************************************************************
    //***************************************************************************************************************
    
    // inputs
    TCanvas* canvasInputScaled = new TCanvas("canvasInputScaled","",0,0,1000,1350);  // gives the page size
    DrawGammaCanvasSettings( canvasInputScaled, 0.1, 0.015, 0.015, 0.07);
    canvasInputScaled->SetLogy();
    
    Float_t maximumPi0Scaled = FindLargestEntryIn1D(fHPt_Particle[0])*10;
    Float_t minimumPi0Scaled = FindSmallestEntryIn1D(fHPt_Particle[nbins-1]);
    
    TH2F * histo2DInputScaledPi0;
    histo2DInputScaledPi0 = new TH2F("histo2DInputScaledPi0","histo2DInputScaledPi0",1000,0., 50,10000,minimumPi0Scaled,maximumPi0Scaled);
    if(bscale){
      SetStyleHistoTH2ForGraphs(histo2DInputScaledPi0, "#it{p}_{T} (GeV/#it{c})",Form("N_{%s} reweighted",particlelegend.Data()),
                                0.032,0.04, 0.032,0.04, 0.8,1);
    }
    else{
      SetStyleHistoTH2ForGraphs(histo2DInputScaledPi0, "#it{p}_{T} (GeV/#it{c})",Form("N_{%s}",particlelegend.Data()),
                                0.032,0.04, 0.032,0.04, 0.8,1);
    }
    histo2DInputScaledPi0->DrawCopy();
    
    TLegend* legendScaled = GetAndSetLegend2(0.25, 0.99-(1.15*nbins/2*0.040), 0.95, 0.99,25);
    legendScaled->SetMargin(0.12);
    legendScaled->SetNColumns(2);
    DrawGammaSetMarker(fHPt_ParticleAll, markerBins[0], 1., colorBins[0], colorBins[0]);
    fHPt_ParticleAll->DrawCopy("e1,same");
    legendScaled->AddEntry(fHPt_ParticleAll,"summed","p");
    if(bscale && bDoMb){
      DrawGammaSetMarker(fHPt_ParticleMB, markerBins[5], 1., colorBins[0], colorBins[0]);
      fHPt_ParticleMB->DrawCopy("e1,same");
      legendScaled->AddEntry(fHPt_ParticleMB,"MB","p");
    }
    for (Int_t i = 0; i< nbins; i++){
      DrawGammaSetMarker(fHPt_Particle[i], markerBins[i+1], 1., colorBins[i+1], colorBins[i+1]);
      fHPt_Particle[i]->DrawCopy("e1,same");
      legendScaled->AddEntry(fHPt_Particle[i],Form("%1.0f GeV/#it{c} < #it{p}^{hard}_{T} < %1.0f GeV/#it{c}", minPtHard[i], maxPtHard[i]),"p");
    }
    legendScaled->Draw();
    //    labelMCName->Draw();
    TLatex *labelSimulation = new TLatex(0.38,0.99-(1.15*(nbins/2+2)*0.04),Form("%s, %s, %s%s",particlelegend.Data(),generator.Data(),tune.Data(),pthard_legend.Data()));
    SetStyleTLatex( labelSimulation, 32,4);
    labelSimulation->SetTextFont(43);
    labelSimulation->Draw();
    
    canvasInputScaled->Update();
    if(bscale){
      canvasInputScaled->SaveAs(Form("%s/%s_MC_InputScaled_%s%s_%s_%s%s.%s",outputDir.Data(),particle.Data(),generator.Data(),tune.Data(),binning.Data(),mbias.Data(),ptcut.Data(),suffix.Data()));
    }
    else{
      canvasInputScaled->SaveAs(Form("%s/%s_MC_InputUnscaled_%s%s_%s_%s%s.%s",outputDir.Data(),particle.Data(),generator.Data(),tune.Data(),binning.Data(),mbias.Data(),ptcut.Data(),suffix.Data()));
    }
    
    // draw accepted for scaled case
    if(bscale){
      SetStyleHistoTH2ForGraphs(histo2DInputScaledPi0, "#it{p}_{T} (GeV/#it{c})",Form("N_{%s}^{accepted} reweighted ",particlelegend.Data()),
                                0.032,0.04, 0.032,0.04, 0.8,1);
      histo2DInputScaledPi0->DrawCopy();
      DrawGammaSetMarker(fHPt_ParticleAccAll, markerBins[0], 1., colorBins[0], colorBins[0]);
      fHPt_ParticleAccAll->DrawCopy("e1,same");
      if(bDoMb){
        DrawGammaSetMarker(fHPt_ParticleAccMB, markerBins[5], 1., colorBins[0], colorBins[0]);
        fHPt_ParticleAccMB->DrawCopy("e1,same");
      }
      for (Int_t i = 0; i< nbins; i++){
        DrawGammaSetMarker(fHPt_ParticleAcc[i], markerBins[i+1], 1., colorBins[i+1], colorBins[i+1]);
        fHPt_ParticleAcc[i]->DrawCopy("e1,same");
      }
      legendScaled->Draw();
      legendScaled->Draw();
      
      TLatex *labelAcceptance = new TLatex(0.45,0.99-(1.15*(nbins/2+2)*0.047),Form("%s",fDetectionProcess.Data()));
      SetStyleTLatex( labelAcceptance, 32,4);
      labelAcceptance->SetTextFont(43);
      labelAcceptance->Draw();
      
      labelSimulation->Draw();
      
      canvasInputScaled->Update();
      canvasInputScaled->SaveAs(Form("%s/%s_MC_InputScaledInAcceptance%s_%s%s_%s_%s%s.%s",outputDir.Data(),particle.Data(),detector.Data(),generator.Data(),tune.Data(),binning.Data(),mbias.Data(),ptcut.Data(),suffix.Data()));
      
    // fractions to total for scaled case

      canvasInputScaled->SetLogx();
      canvasInputScaled->SetLogy(0);
      
      TH2F * histo2DFractionScaledPi0;
      if(bDoMb){
        histo2DFractionScaledPi0 = new TH2F("histo2DFractionScaledPi0","histo2DFractionScaledPi0",1000,1e0, 50,10000,0.,2.);
      }
      else{
        histo2DFractionScaledPi0 = new TH2F("histo2DFractionScaledPi0","histo2DFractionScaledPi0",1000,1e0, 50,10000,0.,1.2);
      }
      SetStyleHistoTH2ForGraphs(histo2DFractionScaledPi0, "#it{p}_{T} (GeV/#it{c})",Form("N_{%s}^{bin}/N_{%s}^{all} reweighted",particlelegend.Data(),particlelegend.Data()),
                                0.032,0.04, 0.032,0.04, 0.8,1);
      histo2DFractionScaledPi0->DrawCopy();
      TLegend* legendFracScaled = GetAndSetLegend2(0.25, 0.99-(1.15*nbins/2*0.040), 0.95, 0.99,25);
      legendFracScaled->SetMargin(0.12);
      legendFracScaled->SetNColumns(2);
      if(bDoMb){
        DrawGammaSetMarker(fHPtFrac_ParticleMB, markerBins[0], 1., colorBins[0], colorBins[0]);
        fHPtFrac_ParticleMB->DrawCopy("e1,same");
        legendFracScaled->AddEntry(fHPtFrac_ParticleMB,"MB Comparison","p");
      }
      for (Int_t i = 0; i< nbins; i++){
        DrawGammaSetMarker(fHPtFrac_Particle[i], markerBins[i+1], 1., colorBins[i+1], colorBins[i+1]);
        fHPtFrac_Particle[i]->DrawCopy("e1,same");
        legendFracScaled->AddEntry(fHPtFrac_Particle[i],Form("%1.0f GeV/#it{c} < #it{p}^{hard}_{T} < %1.0f GeV/#it{c}", minPtHard[i], maxPtHard[i]),"p");
      }
      legendFracScaled->Draw();
      labelSimulation->Draw();
      
      canvasInputScaled->Update();
      canvasInputScaled->SaveAs(Form("%s/%s_MC_InputFractionsLogx_%s%s_%s_%s%s.%s",outputDir.Data(),particle.Data(),generator.Data(),tune.Data(),binning.Data(),mbias.Data(),ptcut.Data(),suffix.Data()));
      canvasInputScaled->SetLogx(0);
      
      canvasInputScaled->Update();
      canvasInputScaled->SaveAs(Form("%s/%s_MC_InputFractionsLinx_%s%s_%s_%s%s.%s",outputDir.Data(),particle.Data(),generator.Data(),tune.Data(),binning.Data(),mbias.Data(),ptcut.Data(),suffix.Data()));
      
      
      
      canvasInputScaled->SetLogy();
    
    // relative statistical error for scaled case

      canvasInputScaled->SetLogy(0);
      TH2F * histo2DRelStatPi0;
      histo2DRelStatPi0 = new TH2F("histo2DRelStatPi0","histo2DRelStatPi0",1000,1e0, 50,10000,0.,0.7);
      SetStyleHistoTH2ForGraphs(histo2DRelStatPi0, "#it{p}_{T} (GeV/#it{c})",Form("%s relative statistical uncertainty",particlelegend.Data()),
                                0.032,0.04, 0.032,0.04, 0.8,1);
      histo2DRelStatPi0->DrawCopy();

      DrawGammaSetMarker(fHPt_ParticleRelStatAll, markerBins[0], 1., colorBins[0], colorBins[0]);
      fHPt_ParticleRelStatAll->DrawCopy("e1,same");

      DrawGammaSetMarker(fHPt_ParticleRelStatAccAll, markerBins[5], 1., colorBins[1], colorBins[1]);
      fHPt_ParticleRelStatAccAll->DrawCopy("p,e1,same");
      
      TLegend* legendRelStat = GetAndSetLegend2(0.25, 0.99-(1.15*0.040), 0.95, 0.99,25);
      legendRelStat->SetMargin(0.12);
      legendRelStat->SetNColumns(2);

      legendRelStat->AddEntry(fHPt_ParticleRelStatAll,Form("rel. stat. err., %s, input",particlelegend.Data()));
      legendRelStat->AddEntry(fHPt_ParticleRelStatAccAll,Form("rel. stat. err., %s, accepted",particlelegend.Data()));
      
      legendRelStat->Draw();
      labelAcceptance->Draw();
      labelSimulation->Draw();
      
      canvasInputScaled->Update();
      canvasInputScaled->SaveAs(Form("%s/%s_MC_RelStatErrorssLinx_%s_%s%s_%s_%s%s.%s",outputDir.Data(),particle.Data(),detector.Data(),generator.Data(),tune.Data(),binning.Data(),mbias.Data(),ptcut.Data(),suffix.Data()));

      canvasInputScaled->SetLogy();
    }
    
    // mode 0: PCM
    // mode 2: PCM-EMC
    // mode 3: PCM-PHOS
    // mode 4: EMC
    // mode 5: PHOS

    // plot pt vs alpha for all pT hard bins
    if(bscale){
      Double_t minimumX = 0.1;
//      if(mode == 2){
//        minimumX = 0.3;
//      }
//      else if(mode == 3){
//        minimumX = 0.25;
//      }
//      else if(mode == 4){
//        minimumX = 0.7;
//      }
//      else if(mode == 5){
//        minimumX = 0.5;
//      }

      if(bDoMb){
        PlotStandard2D( fHPt_Alpha_ParticleAccMB ,
                       Form("%s/MC_%sPt_Alpha_MB_%s_%s%s_%s%s.%s",outputDir.Data(),particle.Data(),detector.Data(),generator.Data(),tune.Data(),binning.Data(),ptcut.Data(),suffix.Data()),
                       "", Form("#it{p}_{%s,T} (GeV/#it{c})",particlelegend.Data()), "#alpha",
                       kFALSE, 0., 1., kTRUE, minimumX, 40., 1, 0, 1,
                       floatLocationRightUp2D,500,500,"MC", "",0);
      }
      PlotStandard2D( fHPt_Alpha_ParticleAccAll ,
                     Form("%s/MC_%sPt_Alpha_Sum_%s_%s%s_%s%s.%s",outputDir.Data(),particle.Data(),detector.Data(),generator.Data(),tune.Data(),binning.Data(),ptcut.Data(),suffix.Data()),
                     "", Form("#it{p}_{%s,T} (GeV/#it{c})",particlelegend.Data()), "#alpha",
                     kFALSE, 0., 1., kTRUE, minimumX, 40., 1, 0, 1,
                     floatLocationRightUp2D,500,500,"MC", "",0);
      for(int ibin=0;ibin<nbins;ibin++){
        
        ptHardRange = Form("%1.0f GeV/#it{c} < #it{p}^{hard}_{T} < %1.0f GeV/#it{c}", minPtHard[ibin], maxPtHard[ibin]);
        
        PlotStandard2D( fHPt_Alpha_ParticleAcc[ibin] ,
                       Form("%s/MC_%sPt_Alpha_Bin%d_%s_%s%s_%s%s.%s",outputDir.Data(),particle.Data(),ibin+1,detector.Data(),generator.Data(),tune.Data(),binning.Data(),ptcut.Data(),suffix.Data()),
                       "", Form("#it{p}_{%s,T} (GeV/#it{c})",particlelegend.Data()), "#alpha",
                       kFALSE, 0., 1., kTRUE, minimumX, 40., 1, 0, 1,
                       floatLocationRightUp2D,500,500,"MC", "",1);
      }
    }
  }
  
  else{
    cout << "this will be implemented later" << endl;
    
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
    // get particle histograms
    fHPt_Y_ParticleMB = (TH2F*) HistosAll->FindObject(Form("Pt_Y_%s",particle.Data()))->Clone(Form("Pt_Y_%s_MB",particle.Data()));
    fHPt_Y_ParticleAccMB = (TH2F*) HistosAll->FindObject(Form("Pt_Y_%sGG%sAcc",particle.Data(),detector.Data()))->Clone(Form("Pt_Y_%sGG%sAcc_MB",particle.Data(),detector.Data()));
    fHPt_Alpha_ParticleAccMB = (TH2F*) HistosAll->FindObject(Form("Pt_Alpha_%sGG%sAcc",particle.Data(),detector.Data()))->Clone(Form("Pt_Alpha_%sGG%sAcc_MB",particle.Data(),detector.Data()));
    // now projections on pT axis
    fHPt_ParticleMB = (TH1D*)fHPt_Y_ParticleMB->ProjectionX(Form("hPt_%s_MB",particle.Data()),1,160)->Clone(Form("hPt_%s_MB",particle.Data()));
    fHPt_ParticleAccMB = (TH1D*)fHPt_Y_ParticleAccMB->ProjectionX(Form("hPt_%sAcc_MB",particle.Data()),1,160)->Clone(Form("hPt_%sAcc_MB",particle.Data()));
    
    return;
  }
  
  
}
