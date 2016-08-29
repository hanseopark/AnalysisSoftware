// provided by Gamma Conversion Group, $ALICE_ROOT/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
// ***************************************************************************************************************
// **   Friederike Bock, friederike.bock@cern.ch                                                                **
// **   Lucas Altenkaemper, lucas.altenkaemper@cern.ch                                                          **
// ***************************************************************************************************************

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
#include "TGaxis.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
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
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "PrepareCocktail.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
//#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"

void PrepareCocktail(   TString nameFileCocktail    = "",
                        TString nameFilePi0         = "",
                        TString suffix              = "eps",
                        TString cutSelection        = "",
                        TString option              = "",
                        Double_t rapidity           = 0.80,
                        TString quantity            = "",
                        TString period              = "",
                        Int_t numberOfBins          = 30,
                        Int_t mode                  = 0
                     ) {
    
    gROOT->Reset();
    
    //************************** Set general style settings *********************************************************
    StyleSettingsThesis();
    SetPlotStyle();
    
    //************************** Set output directory ***************************************************************
    TString outputDir                                           = Form("%s/%s/%s/PrepareCocktail",cutSelection.Data(),option.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    //************************** Set global variables ***************************************************************
    //fDate                                                     = ReturnDateString();
    fEnergyFlag                                                 = option;
    fPeriodFlag                                                 = period;
    fSuffix                                                     = suffix;
    fMode                                                       = mode;
    fRapidity                                                   = rapidity;
    cout << "Pictures are saved as " << suffix.Data() << endl;
    
    //***************************** Separate cutstrings *************************************************************
    if(cutSelection.Length() == 0){
        cout<<"ERROR: Cut selection is not set, please do!"<<endl;
        return;
    }
    fCutSelection                                               = cutSelection;
    ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, fMode);
    
    //**************************** Determine Centrality *************************************************************
    centralityString                                            = GetCentralityString(fEventCutSelection);
    if (centralityString.CompareTo("pp")==0){
        fTextCent                                               = "MinBias";
    } else {
        fTextCent                                               = Form("%s central", centralityString.Data());
    }
    if (centralityString.CompareTo("pp")!=0 && !centralityString.Contains("0-100%") ){
        fCollisionSystem                                        = Form("%s %s", centralityString.Data(), fCollisionSystem.Data());
    }
    
    //***************************** Load binning for spectrum *******************************************************
    Initialize(fEnergyFlag, numberOfBins);

    //***************************** Spectra quantity ****************************************************************
    if (quantity.CompareTo("") == 0) {
        cout << "ERROR: Quantity of cocktail spectra not specified, returning!" << endl;
        return;
    }
    TString quantityLatex                                       = "";
    if (quantity.CompareTo("dNdydpT") == 0 || quantity.CompareTo("dNdydpt") == 0) {
        quantityLatex                                           = "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})";
    } else if (quantity.CompareTo("invYield") == 0) {
        quantityLatex                                           = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})";
    }
    
    //***************************** Cocktail file *******************************************************************
    TFile fileCocktail(nameFileCocktail.Data());
    TDirectoryFile* topDirCocktail                              = (TDirectoryFile*)fileCocktail.Get("GammaCocktailMC");
    if (!topDirCocktail) {
        cout << "ERROR: TopDirCocktail not found!" << endl;
        return;
    }
    TList* histoListCocktail                                    = (TList*)topDirCocktail->Get(Form("GammaCocktailMC_%.2f", rapidity));
    cout << "searching for " << Form("GammaCocktailMC_%.2f", rapidity) << endl;
    if (!histoListCocktail) {
        cout << "ERROR: Folder with rapidity " << rapidity << " not contained in file cocktail file!" << endl;
        return;
    }
    delete topDirCocktail;

    //***************************** Pi0 file ************************************************************************
    if (nameFilePi0.CompareTo("") != 0) {
        TFile* filePi0                                          = new TFile(nameFilePi0);
        histoPi0YieldData                                       = (TH1D*)filePi0->Get("CorrectedYieldTrueEff");
    } else {
        cout << "WARNING: No pi0 file specified!" << endl;
        histoPi0YieldData                                       = NULL;
    }
    
    //***************************** Get number of events ************************************************************
    histoNEvents                                                = (TH1F*)histoListCocktail->FindObject("NEvents");
    nEvents                                                     = histoNEvents->GetEntries();
    cout << nEvents << " events" << endl;

    //***************************** Read histograms from cocktail file **********************************************
    histoGammaSumPtY                                            = (TH2F*)histoListCocktail->FindObject("Pt_Y_Gamma");
    histoGammaSumPtPhi                                          = (TH2F*)histoListCocktail->FindObject("Pt_Phi_Gamma");
    if (histoGammaSumPtY) {
        histoGammaSumPtY->SetName("Gamma_Pt_Y_OrBin");
        histoGammaSumPtY->Sumw2();
    } else {
        cout << "ERROR: Gamma histo not found!" << endl;
        return;
    }
    
    histoDecayChannels                                          = new TH1F*[nMotherParticles];
    histoGammaPtY                                               = new TH2F*[nMotherParticles];
    histoGammaPtPhi                                             = new TH2F*[nMotherParticles];
    histoGammaMotherPtY                                         = new TH2F*[nMotherParticles];
    histoGammaMotherPtPhi                                       = new TH2F*[nMotherParticles];
//    histoPtGammaPtMother                                      = new TH2F*[nMotherParticles];
//    histoPhiGammaPhiMother                                    = new TH2F*[nMotherParticles];
//    histoGammaMotherPtDeltaPhi                                = new TH2F*[nMotherParticles];
//    histoGammaMotherPtAlpha                                   = new TH2F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        histoDecayChannels[i]                                   = (TH1F*)histoListCocktail->FindObject(Form("DecayChannels_%s", motherParticles[i].Data()));
        
        histoGammaPtY[i]                                        = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_Gamma_From_%s", motherParticles[i].Data()));
        if (histoGammaPtY[i]->GetEntries()) {
            histoGammaPtY[i]->SetName(Form("Gamma_From_%s_Pt_Y_OrBin", motherParticles[i].Data()));
            histoGammaPtY[i]->Sumw2();
            hasGammasFromSource[i]                              = kTRUE;
        } else
            histoGammaPtY[i]                                    = NULL;
      
        histoGammaPtPhi[i]                                      = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_Gamma_From_%s", motherParticles[i].Data()));
        if (histoGammaPtPhi[i]->GetEntries() && hasGammasFromSource[i]) {
            histoGammaPtPhi[i]->SetName(Form("Gamma_From_%s_Pt_Phi_OrBin", motherParticles[i].Data()));
            histoGammaPtPhi[i]->Sumw2();
        } else
            histoGammaPtPhi[i]                                  = NULL;
        
        histoGammaMotherPtY[i]                                  = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_%s", motherParticles[i].Data()));
        if (histoGammaMotherPtY[i]->GetEntries() && hasGammasFromSource[i]) {
            histoGammaMotherPtY[i]->SetName(Form("%s_Pt_Y_OrBin", motherParticles[i].Data()));
            histoGammaMotherPtY[i]->Sumw2();
        } else
            histoGammaMotherPtY[i]                              = NULL;

        histoGammaMotherPtPhi[i]                                = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_%s", motherParticles[i].Data()));
        if (histoGammaMotherPtPhi[i]->GetEntries() && hasGammasFromSource[i]) {
            histoGammaMotherPtPhi[i]->SetName(Form("%s_Pt_Phi_OrBin", motherParticles[i].Data()));
            histoGammaMotherPtPhi[i]->Sumw2();
        } else
            histoGammaMotherPtPhi[i]                            = NULL;
    
//        histoPtGammaPtMother[i]                               = (TH2F*)histoListCocktail->FindObject(Form("PtGamma_PtMother_%s", motherParticles[i].Data()));
//        if (histoPtGammaPtMother[i]->GetEntries() && hasGammasFromSource[i]) {
//            histoPtGammaPtMother[i]->SetName(Form("%s_PtGamma_PtMother_OrBin", motherParticles[i].Data()));
//            histoPtGammaPtMother[i]->Sumw2();
//        } else
//            histoPtGammaPtMother[i]                           = NULL;
//        
//        histoPhiGammaPhiMother[i]                             = (TH2F*)histoListCocktail->FindObject(Form("PhiGamma_PhiMother_%s", motherParticles[i].Data()));
//        if (histoPhiGammaPhiMother[i]->GetEntries() && hasGammasFromSource[i]) {
//            histoPhiGammaPhiMother[i]->SetName(Form("%s_PhiGamma_PhiMother_OrBin", motherParticles[i].Data()));
//            histoPhiGammaPhiMother[i]->Sumw2();
//        } else
//            histoPhiGammaPhiMother[i]                         = NULL;
        
        // histos affected bei lightweight output mode
//        histoGammaMotherPtDeltaPhi[i]                         = (TH2F*)histoListCocktail->FindObject(Form("Pt_DeltaPhi_%s", motherParticles[i].Data()));
//        if (histoGammaMotherPtDeltaPhi[i] && hasGammasFromSource[i]) {
//            if (histoGammaMotherPtDeltaPhi[i]->GetEntries()) {
//                histoGammaMotherPtDeltaPhi[i]->SetName(Form("%s_Pt_DeltaPhi_OrBin", motherParticles[i].Data()));
//                histoGammaMotherPtDeltaPhi[i]->Sumw2();
//            } else
//                histoGammaMotherPtDeltaPhi[i]                 = NULL;
//        }
//        
//        histoGammaMotherPtAlpha[i]                            = (TH2F*)histoListCocktail->FindObject(Form("Pt_Alpha_%s", motherParticles[i].Data()));
//        if (histoGammaMotherPtAlpha[i] && hasGammasFromSource[i]) {
//            if (histoGammaMotherPtAlpha[i]->GetEntries()) {
//                histoGammaMotherPtAlpha[i]->SetName(Form("%s_Pt_Alpha_OrBin", motherParticles[i].Data()));
//                histoGammaMotherPtAlpha[i]->Sumw2();
//            } else
//                histoGammaMotherPtAlpha[i]                    = NULL;
//        }
    }
    
    //***************************** Read params from cocktail input file (or mT scaling) ****************************  // <<== this is not complete, which parametrization should be taken?
    cocktailInputList                                           = GetCocktailInput(fEnergyFlag, centralityString);
    cocktailInputParametrizations                               = new TF1*[nMotherParticles];
    cocktailInputParametrizationsMtScaled                       = new TF1*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (cocktailInputList->FindObject(Form("%sCombStat_Fit", motherParticlesCocktailInput[i].Data()))) {
            cocktailInputParametrizations[i]                    = (TF1*)cocktailInputList->FindObject(Form("%sCombStat_Fit", motherParticlesCocktailInput[i].Data()));
            cocktailInputParametrizations[i]->SetName(Form("%s_Pt_Param", motherParticles[i].Data()));
            if (i>0) {
                cocktailInputParametrizationsMtScaled[i]        = (TF1*)MtScaledParam(cocktailInputParametrizations[0], i);
                cocktailInputParametrizationsMtScaled[i]->SetName(Form("%s_Pt_Param_mTscaled", motherParticles[i].Data()));
            }
        } else if (cocktailInputList->FindObject(Form("%sStat_Fit", motherParticlesCocktailInput[i].Data()))) {
            cocktailInputParametrizations[i]                    = (TF1*)cocktailInputList->FindObject(Form("%sStat_Fit", motherParticlesCocktailInput[i].Data()));
            cocktailInputParametrizations[i]->SetName(Form("%s_Pt_Param", motherParticles[i].Data()));
            if (i>0) {
                cocktailInputParametrizationsMtScaled[i]        = (TF1*)MtScaledParam(cocktailInputParametrizations[0], i);
                cocktailInputParametrizationsMtScaled[i]->SetName(Form("%s_Pt_Param_mTscaled", motherParticles[i].Data()));
            }
        } else {
            cocktailInputParametrizations[i]                    = NULL;
            cocktailInputParametrizationsMtScaled[i]            = (TF1*)MtScaledParam(cocktailInputParametrizations[0], i);
            cocktailInputParametrizationsMtScaled[i]->SetName(Form("%s_Pt_Param_mTscaled", motherParticles[i].Data()));
        }
    }
    
    //***************************** Get number of spectra ***********************************************************
    Int_t nSpectra                                              = 0;
    for (Int_t i=0; i<nMotherParticles; i++)
        if (hasGammasFromSource[i]) nSpectra++;
    Int_t nRows                                                 = 0;
    if (nSpectra%2 == 0) nRows                                  = nSpectra/2 + 1;
    else nRows                                                  = (nSpectra+1)/2 + 1;

    //***************************** Get decay channels  *************************************************************
    TString tempBinLabel                                        = "";
    Int_t counter                                               = 0;
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoDecayChannels[i]) {
            for (Int_t bin=2; bin<20; bin++) {
                if (histoDecayChannels[i]->GetBinContent(bin)) {
                    tempBinLabel                                = histoDecayChannels[i]->GetXaxis()->GetBinLabel(bin);
                    if (tempBinLabel.Contains("#gamma")) {
                        if (counter==0)
                            decayChannelsLatex[i]               = tempBinLabel.Data();
                        else if (counter==1)
                            decayChannelsLatex[i]               = decayChannelsLatex[i] + " (" + tempBinLabel;
                        else
                            decayChannelsLatex[i]               = decayChannelsLatex[i] + ", " + tempBinLabel;
                        counter++;
                    }
                }
            }
            if (counter>1) decayChannelsLatex[i]                = Form("%s)", decayChannelsLatex[i].Data());
            counter                                             = 0;
        } else {
            decayChannelsLatex[i]                               = "";
        }
    }
    
    //***************************** Project from 2D histograms ******************************************************
    histoGammaSumPtOrBin                                        = (TH1F*)histoGammaSumPtY->ProjectionX("Gamma_Pt_OrBin", histoGammaSumPtY->GetYaxis()->FindBin(-rapidity), histoGammaSumPtY->GetYaxis()->FindBin(rapidity), "e");
    histoGammaSumYOrBin                                         = (TH1F*)histoGammaSumPtY->ProjectionY("Gamma_Y_OrBin", 1, histoGammaSumPtY->GetNbinsX(), "e");
    histoGammaSumPhiOrBin                                       = (TH1F*)histoGammaSumPtPhi->ProjectionY("Gamma_Phi_OrBin", 1, histoGammaSumPtPhi->GetNbinsX(), "e");
    
    histoGammaPtOrBin                                           = new TH1F*[nMotherParticles];
    histoGammaYOrBin                                            = new TH1F*[nMotherParticles];
    histoGammaPhiOrBin                                          = new TH1F*[nMotherParticles];
    histoGammaMotherPtOrBin                                     = new TH1F*[nMotherParticles];
    histoGammaMotherYOrBin                                      = new TH1F*[nMotherParticles];
    histoGammaMotherPhiOrBin                                    = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtY[i]) {
            histoGammaPtOrBin[i]                                = (TH1F*)histoGammaPtY[i]->ProjectionX(Form("Gamma_From_%s_Pt_OrBin", motherParticles[i].Data()), histoGammaPtY[i]->GetYaxis()->FindBin(-rapidity), histoGammaPtY[i]->GetYaxis()->FindBin(rapidity), "e");
            histoGammaYOrBin[i]                                 = (TH1F*)histoGammaPtY[i]->ProjectionY(Form("Gamma_From_%s_Y_OrBin", motherParticles[i].Data()), 1, histoGammaPtY[i]->GetNbinsX(), "e");
        } else {
            histoGammaPtOrBin[i]                                = NULL;
            histoGammaYOrBin[i]                                 = NULL;
        }
        if (histoGammaPtPhi[i])
            histoGammaPhiOrBin[i]                               = (TH1F*)histoGammaPtPhi[i]->ProjectionY(Form("Gamma_From_%s_Phi_OrBin", motherParticles[i].Data()), 1, histoGammaPtPhi[i]->GetNbinsX(), "e");
        else
            histoGammaPhiOrBin[i]                               = NULL;
        
        if (histoGammaMotherPtY[i]) {
            histoGammaMotherPtOrBin[i]                          = (TH1F*)histoGammaMotherPtY[i]->ProjectionX(Form("%s_Pt_OrBin", motherParticles[i].Data()), histoGammaMotherPtY[i]->GetYaxis()->FindBin(-rapidity), histoGammaMotherPtY[i]->GetYaxis()->FindBin(rapidity), "e");
            histoGammaMotherYOrBin[i]                           = (TH1F*)histoGammaMotherPtY[i]->ProjectionY(Form("%s_Y_OrBin", motherParticles[i].Data()), 1, histoGammaMotherPtY[i]->GetNbinsX(), "e");
        } else {
            histoGammaMotherPtOrBin[i]                          = NULL;
            histoGammaMotherYOrBin[i]                           = NULL;
        }
        if (histoGammaMotherPtPhi[i])
            histoGammaMotherPhiOrBin[i]                         = (TH1F*)histoGammaMotherPtPhi[i]->ProjectionY(Form("%s_Phi_OrBin", motherParticles[i].Data()), 1, histoGammaMotherPtPhi[i]->GetNbinsX(), "e");
        else
            histoGammaMotherPhiOrBin[i]                         = NULL;
    }
    
    //***************************** Rebin pt spectra ****************************************************************
    histoGammaSumPt                                             = (TH1F*)histoGammaSumPtOrBin->Clone("Gamma_Pt");
    histoGammaSumPt->Sumw2();
    histoGammaSumPt->GetXaxis()->SetRangeUser(ptMin, ptMax);
    RebinSpectrum(histoGammaSumPt,"");
    
    histoGammaPt                                                = new TH1F*[nMotherParticles];
    histoGammaMotherPt                                          = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i]) {
            histoGammaPt[i]                                     = (TH1F*)histoGammaPtOrBin[i]->Clone(Form("Gamma_From_%s_Pt", motherParticles[i].Data()));
            histoGammaPt[i]->Sumw2();
            histoGammaPt[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            RebinSpectrum(histoGammaPt[i],"");
        } else {
            histoGammaPt[i]                                     = NULL;
        }
        if (histoGammaMotherPtOrBin[i]) {
            histoGammaMotherPt[i]                               = (TH1F*)histoGammaMotherPtOrBin[i]->Clone(Form("%s_Pt", motherParticles[i].Data()));
            histoGammaMotherPt[i]->Sumw2();
            histoGammaMotherPt[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            RebinSpectrum(histoGammaMotherPt[i],"");
        } else {
            histoGammaMotherPt[i]                               = NULL;
        }
    }

    //***************************** Scale spectra with nEvents ******************************************************
    histoGammaSumPt->Scale(1./nEvents);
    histoGammaSumPtOrBin->Scale(1./nEvents);
    histoGammaSumYOrBin->Scale(1./nEvents);
    histoGammaSumPhiOrBin->Scale(1./nEvents);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPt[i])                histoGammaPt[i]->Scale(1./nEvents);
        if (histoGammaPtOrBin[i])           histoGammaPtOrBin[i]->Scale(1./nEvents);
        if (histoGammaYOrBin[i])            histoGammaYOrBin[i]->Scale(1./nEvents);
        if (histoGammaMotherPt[i])          histoGammaMotherPt[i]->Scale(1./nEvents);
        if (histoGammaMotherPtOrBin[i])     histoGammaMotherPtOrBin[i]->Scale(1./nEvents);
        if (histoGammaMotherYOrBin[i])      histoGammaMotherYOrBin[i]->Scale(1./nEvents);
        if (histoGammaMotherPhiOrBin[i])    histoGammaMotherPhiOrBin[i]->Scale(1./nEvents);
    }

    //***************************** Calculate ratio to input parametrizations ***************************************
//    histoGammaMotherPtOrBinRatioToParam                         = new TH1F*[nMotherParticles];
//    for (Int_t i=0; i<nMotherParticles; i++) {
//        if (histoGammaMotherPtOrBin[i] && cocktailInputParametrizations[i]) {
//            histoGammaMotherPtOrBinRatioToParam[i]              = (TH1F*)histoGammaMotherPtOrBin[i]->Clone(Form("%s_RatioToParam", histoGammaMotherPtOrBin[i]->GetName()));
//            histoGammaMotherPtOrBinRatioToParam[i]->Sumw2();
//            histoGammaMotherPtOrBinRatioToParam[i]->Divide(cocktailInputParametrizations[i]);
//        } else {
//            histoGammaMotherPtOrBinRatioToParam[i]              = NULL;
//        }
//    }

    //***************************** Transform yields ****************************************************************
    if (quantity.CompareTo("dNdydpT") == 0 || quantity.CompareTo("dNdydpt") == 0) {
        histoGammaSumPt                                         = ConvertYieldHisto(histoGammaSumPt,            kTRUE, kTRUE, kFALSE, kFALSE);
        for (Int_t i=0; i<nMotherParticles; i++) {
            if (histoGammaPt[i])        histoGammaPt[i]         = ConvertYieldHisto(histoGammaPt[i],            kTRUE, kTRUE, kFALSE, kFALSE);
            if (histoGammaMotherPt[i])  histoGammaMotherPt[i]   = ConvertYieldHisto(histoGammaMotherPt[i],      kTRUE, kTRUE, kFALSE, kFALSE);
        }
    }
    
    TH1D* dummyHist                                             = NULL;
    //***************************** Plot cocktail mothers (pt) ******************************************************
    TCanvas *canvasMothers                                      = new TCanvas("canvasMothers","",1100,1200);
    DrawGammaCanvasSettings(canvasMothers, 0.165, 0.02, 0.02, 0.09);
    canvasMothers->SetLogy();
        
    TLegend* legendMothers                                      = GetAndSetLegend2(0.5, 0.92-(0.045*nRows), 0.9, 0.92, 40, 2);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, histoGammaMotherPtOrBin[0]->GetXaxis()->GetXmin(), histoGammaMotherPtOrBin[0]->GetXaxis()->GetXmax());
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", quantityLatex, 1e-10, 2, 1.0, 1.8);
    dummyHist->Draw();
        
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherPtOrBin[i]) {
            DrawGammaSetMarker(         histoGammaMotherPtOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothers->AddEntry(    histoGammaMotherPtOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "l");
            histoGammaMotherPtOrBin[i]->Draw("csamehist");
        }
    }
    legendMothers->Draw("same");
        
    canvasMothers->SaveAs(Form("%s/CocktailMothers_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothers;
    delete canvasMothers;

    //***************************** Plot cocktail mothers + input param (pt) ****************************************
    TCanvas *canvasMothersParam                                 = new TCanvas("canvasMothersParam","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersParam, 0.165, 0.02, 0.02, 0.09);
    canvasMothersParam->SetLogy();
    
    TLegend* legendMothersParam                                 = GetAndSetLegend2(0.5, 0.92-(0.045*nRows), 0.9, 0.92, 40, 2);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, histoGammaMotherPtOrBin[0]->GetXaxis()->GetXmin(), histoGammaMotherPtOrBin[0]->GetXaxis()->GetXmax());
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", quantityLatex, 1e-10, 2, 1.0, 1.8);
    dummyHist->Draw();
    
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherPtOrBin[i]) {
            if (cocktailInputParametrizations[i]) {
                cocktailInputParametrizations[i]->SetLineColor(cocktailColor[i]);
                cocktailInputParametrizations[i]->SetLineStyle(2);
                cocktailInputParametrizations[i]->Draw("same");
            }
            legendMothersParam->AddEntry(histoGammaMotherPtOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "l");
            histoGammaMotherPtOrBin[i]->Draw("csamehist");
        }
    }
    legendMothersParam->Draw("same");
    
    canvasMothersParam->SaveAs(Form("%s/CocktailMothersInclParam_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothersParam;
    delete canvasMothersParam;
    
    //***************************** Plot ratio cocktail mothers to pi0 (pt) *****************************************
    TCanvas *canvasMothersRatio                                 = new TCanvas("canvasMothersRatio","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersRatio, 0.1, 0.02, 0.02, 0.09);
    canvasMothersRatio->SetLogy();
    
    TLegend* legendMothersRatio                                 = GetAndSetLegend2(0.5, 0.92-(0.045*nRows), 0.9, 0.92, 40, 2);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, histoGammaMotherPtOrBin[0]->GetXaxis()->GetXmin(), histoGammaMotherPtOrBin[0]->GetXaxis()->GetXmax());
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "ratio X / #pi^{0}", 1e-4, 10, 1.0, 1.2);
    dummyHist->Draw();
    
    TH1F* tempRatio                                             = NULL;
    for (Int_t i=1; i<nMotherParticles; i++) {
        if (histoGammaMotherPtOrBin[0] && histoGammaMotherPtOrBin[i]) {
            tempRatio                                           = (TH1F*)histoGammaMotherPtOrBin[i]->Clone("tempRatio");
            tempRatio->Sumw2();
            tempRatio->Divide(histoGammaMotherPtOrBin[0]);
            DrawGammaSetMarker(             tempRatio, cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothersRatio->AddEntry(   tempRatio, Form("%s / #pi^{0}", motherParticlesLatex[i].Data()), "l");
            tempRatio->Draw("csamehist");
        }
    }
    legendMothersRatio->Draw("same");
    
    canvasMothersRatio->SaveAs(Form("%s/CocktailMothersRatioToPi0_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete tempRatio;
    delete legendMothersRatio;
    delete canvasMothersRatio;
    
    //***************************** Plot cocktail gammas (pt) *******************************************************
    TCanvas *CocktailGammas                                     = new TCanvas("CocktailGammas","",1100,1200);
    DrawGammaCanvasSettings(CocktailGammas, 0.165, 0.02, 0.02, 0.09);
    CocktailGammas->SetLogy();
    
    TLegend* legendGammas                                       = GetAndSetLegend2(0.4, 0.92-(0.045*nRows*1.6), 0.75, 0.95, 40);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, histoGammaPtOrBin[0]->GetXaxis()->GetXmin(), histoGammaPtOrBin[0]->GetXaxis()->GetXmax());
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", quantityLatex, 1e-10, 2, 1.0, 1.8);
    dummyHist->Draw();
    
    DrawGammaSetMarker(     histoGammaSumPtOrBin, 20, 1, kBlack,  kBlack);
    legendGammas->AddEntry( histoGammaSumPtOrBin, "all #gamma", "l");
    histoGammaSumPtOrBin->Draw("chistsame");
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i]) {
            DrawGammaSetMarker(         histoGammaPtOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendGammas->AddEntry(     histoGammaPtOrBin[i], Form("%s #rightarrow %s", motherParticlesLatex[i].Data(), decayChannelsLatex[i].Data()), "l");
            histoGammaPtOrBin[i]->Draw("csamehist");
        }
    }
    legendGammas->Draw("same");
    
    CocktailGammas->SaveAs(Form("%s/CocktailGammas_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendGammas;
    delete CocktailGammas;

    //***************************** Plot cocktail gammas to pi0 ratio ***********************************************
    if (histoGammaMotherPtOrBin[0]) {
        TCanvas *canvasGammasRatio                              = new TCanvas("canvasGammasRatio","",1100,1200);
        DrawGammaCanvasSettings(canvasGammasRatio, 0.1, 0.02, 0.02, 0.09);
        canvasGammasRatio->SetLogy();
        
        TLegend* legendGammasRatio                              = GetAndSetLegend2(0.4, 0.92-(0.045*nRows*1.6), 0.75, 0.95, 40);
        dummyHist                                               = new TH1D("dummyHist", "", 1000, histoGammaSumPtOrBin->GetXaxis()->GetXmin(), histoGammaSumPtOrBin->GetXaxis()->GetXmax());
        SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "ratio #gamma / #pi^{0} from different sources", 1e-4, 10, 1.0, 1.2);
        dummyHist->Draw();
        
        TH1F* tempRatioGammas                                   = (TH1F*)histoGammaSumPtOrBin->Clone("tempRatioGammas");
        tempRatioGammas->Sumw2();
        DrawGammaSetMarker(             tempRatioGammas, 20, 1, kBlack,  kBlack);
        legendGammasRatio->AddEntry(    tempRatioGammas, "all #gamma", "l");
        tempRatioGammas->Draw("chistsame");
        for (Int_t i=0; i<nMotherParticles; i++) {
            if (histoGammaPtOrBin[i]) {
                tempRatioGammas                                 = (TH1F*)histoGammaMotherPtOrBin[i]->Clone("tempRatioGammas");
                tempRatioGammas->Sumw2();
                tempRatioGammas->Divide(histoGammaMotherPtOrBin[0]);
                DrawGammaSetMarker(             tempRatioGammas, cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
                legendGammasRatio->AddEntry(    tempRatioGammas, Form("%s #rightarrow %s", motherParticlesLatex[i].Data(), decayChannelsLatex[i].Data()), "l");
                tempRatioGammas->Draw("csamehist");
            }
        }
        legendGammasRatio->Draw("same");
        
        canvasGammasRatio->SaveAs(Form("%s/CocktailGammasRatioToPi0_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
        delete tempRatioGammas;
        delete legendGammasRatio;
        delete canvasGammasRatio;
    }
    
    //***************************** Plot cocktail gammas to all gammas ratio ****************************************
    TCanvas *canvasGammasRatio2                                 = new TCanvas("canvasGammasRatio2","",1100,1200);
    DrawGammaCanvasSettings(canvasGammasRatio2, 0.1, 0.02, 0.02, 0.09);
    canvasGammasRatio2->SetLogy();
        
    TLegend* legendGammasRatio2                                 = GetAndSetLegend2(0.4, 0.92-(0.045*nRows*1.6), 0.75, 0.95, 40);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, histoGammaSumPtOrBin->GetXaxis()->GetXmin(), histoGammaSumPtOrBin->GetXaxis()->GetXmax());
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "ratio #gamma from X / all", 1e-4, 10, 1.0, 1.2);
    dummyHist->Draw();
    
    TH1F* tempRatioGammas2                                      = NULL;
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i]) {
            tempRatioGammas2                                    = (TH1F*)histoGammaPtOrBin[i]->Clone("tempRatioGammas2");
            tempRatioGammas2->Sumw2();
            tempRatioGammas2->Divide(histoGammaSumPtOrBin);
            DrawGammaSetMarker(             tempRatioGammas2, cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendGammasRatio2->AddEntry(   tempRatioGammas2, Form("%s #rightarrow %s", motherParticlesLatex[i].Data(), decayChannelsLatex[i].Data()), "l");
            tempRatioGammas2->Draw("csamehist");
        }
    }
    legendGammasRatio2->Draw("same");
        
    canvasGammasRatio2->SaveAs(Form("%s/CocktailGammasRatioToAll_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete tempRatioGammas2;
    delete legendGammasRatio2;
    delete canvasGammasRatio2;
    
    //***************************** Plot cocktail mothers + mT scaled (if input param available) ********************
    TCanvas *canvasMothersParamMt                               = new TCanvas("canvasMothersParamMt","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersParamMt, 0.165, 0.02, 0.02, 0.09);
    canvasMothersParamMt->SetLogy();
    
    Int_t numberOfSpectra                                       = 0;
    for (Int_t i=1; i<nMotherParticles; i++) {
        if (cocktailInputParametrizations[i]) numberOfSpectra++;
    }
    
    TLegend* legendMothersParamMt                               = GetAndSetLegend2(0.5, 0.92-(0.045*numberOfSpectra*3), 0.9, 0.92, 40);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, histoGammaMotherPtOrBin[0]->GetXaxis()->GetXmin(), histoGammaMotherPtOrBin[0]->GetXaxis()->GetXmax());
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", quantityLatex, 1e-10, 2, 1.0, 1.8);
    dummyHist->Draw();
    
    for (Int_t i=1; i<nMotherParticles; i++) {
        if (histoGammaMotherPtOrBin[i]) {
            if (cocktailInputParametrizations[i]) {
                cocktailInputParametrizations[i]->SetLineColor(cocktailColor[i]);
                cocktailInputParametrizations[i]->SetLineStyle(2);
                cocktailInputParametrizations[i]->Draw("same");
                
                cocktailInputParametrizationsMtScaled[i]->SetLineColor(cocktailColor[i]);
                cocktailInputParametrizationsMtScaled[i]->SetLineStyle(4);
                cocktailInputParametrizationsMtScaled[i]->Draw("same");
                
                legendMothersParamMt->AddEntry(histoGammaMotherPtOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "l");
                legendMothersParamMt->AddEntry(cocktailInputParametrizations[i], Form("%s input param.", motherParticlesLatex[i].Data()), "l");
                legendMothersParamMt->AddEntry(cocktailInputParametrizationsMtScaled[i], Form("%s m_{T} scaled param.", motherParticlesLatex[i].Data()), "l");
                histoGammaMotherPtOrBin[i]->Draw("csamehist");
            }
        }
    }
    legendMothersParamMt->Draw("same");
    
    canvasMothersParamMt->SaveAs(Form("%s/CocktailMothersInclParamMt_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothersParamMt;
    delete canvasMothersParamMt;
    
    //***************************** Plot mT scaling cross check *****************************************************
    TCanvas *canvasMtCrossCheck                                 = new TCanvas("canvasMtCrossCheck","",1100,1200);
    DrawGammaCanvasSettings(canvasMtCrossCheck, 0.165, 0.02, 0.02, 0.09);
    canvasMtCrossCheck->SetLogy();
    
    TLegend* legendMtCrossCheck                                 = GetAndSetLegend2(0.5, 0.92-(0.045*2), 0.9, 0.92, 40);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, histoGammaMotherPtOrBin[3]->GetXaxis()->GetXmin(), histoGammaMotherPtOrBin[3]->GetXaxis()->GetXmax());
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", quantityLatex, 1e-10, 2, 1.0, 1.8);
    dummyHist->Draw();
    
    legendMothersParamMt->AddEntry(histoGammaMotherPtOrBin[3], Form("%s", motherParticlesLatex[3].Data()), "l");
    legendMothersParamMt->AddEntry(cocktailInputParametrizationsMtScaled[3], Form("%s m_{T} scaled param.", motherParticlesLatex[3].Data()), "l");

    cocktailInputParametrizationsMtScaled[3]->Draw("same");
    histoGammaMotherPtOrBin[3]->Draw("csamehist");
    legendMtCrossCheck->Draw("same");
    
    canvasMtCrossCheck->SaveAs(Form("%s/MtScalingOmega_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMtCrossCheck;
    delete canvasMtCrossCheck;
    
    //***************************** Plot pi0 from data vs. cocktail *************************************************
    if (histoPi0YieldData) {
        
        TCanvas *canvasPi0                                          = new TCanvas("canvasPi0","",1100,1200);
        DrawGammaCanvasSettings(canvasPi0, 0.165, 0.02, 0.02, 0.09);
        canvasPi0->SetLogy();

        TLegend* legendPi0                                          = GetAndSetLegend2(0.5, 0.92-(0.045*2), 0.9, 0.92, 40);
        dummyHist                                                   = new TH1D("dummyHist", "", 1000, histoGammaMotherPt[0]->GetXaxis()->GetXmin(), histoGammaMotherPt[0]->GetXaxis()->GetXmax());
        SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", histoPi0YieldData->GetMinimum(0)*0.1, histoPi0YieldData->GetMaximum()*2, 1.0, 1.8);
        dummyHist->Draw();

        DrawGammaSetMarker(histoPi0YieldData, 24, 1, kBlack,  kBlack);

        legendPi0->AddEntry(histoPi0YieldData,      Form("%s data", motherParticlesLatex[0].Data()), "p");
        legendPi0->AddEntry(histoGammaMotherPt[0],  Form("%s cocktail", motherParticlesLatex[0].Data()), "l");

        histoPi0YieldData->Draw("same");
        histoGammaMotherPt[0]->Draw("csamehist");
        legendPi0->Draw("same");

        canvasPi0->SaveAs(Form("%s/Pi0DataCocktail_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
        delete legendPi0;
        delete canvasPi0;
    }
    delete dummyHist;
    
    //***************************** Save histograms *****************************************************************
    SaveHistos();
    
    //***************************** Delete objects ******************************************************************
    DeleteObjects();

}

//************************** Initialize binning *********************************************************************
void Initialize(TString energy, Int_t numberOfBins){
    
    InitializeBinning("Pi0", numberOfBins, energy, "directPhoton", fMode, fEventCutSelection, fClusterCutSelection);
    
    fDeltaPt                                    = new TH1F("deltaPt","",fNBinsPt,fBinsPt);
    for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPt+1;iPt++){
        fDeltaPt->SetBinContent(iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
        fDeltaPt->SetBinError(iPt,0);
    }
    
    ptMin                                       = fBinsPt[0];
    ptMax                                       = fBinsPt[fNBinsPt];
    cout << "Using " << ptMin << " <= p_T <= " << ptMax << " GeV/c" << endl;
}

//************************** Rebin spectrum *************************************************************************
void RebinSpectrum(TH1F *Spectrum, TString NewName){
    if(NewName.CompareTo(""))
        NewName = Spectrum->GetName();
    
    *Spectrum = *((TH1F*)Spectrum->Rebin(fNBinsPt,NewName,fBinsPt));
    Spectrum->Divide(fDeltaPt);
}

//************************** Convert yield histo ********************************************************************
TH1F* ConvertYieldHisto(TH1F* input, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt){
    
    if (!input) {
        cout << "Error: Histogram is NULL" << endl;
        return NULL;
    }
    
    Int_t nBins                             = input->GetNbinsX();
    Double_t newValue                       = 0;
    Double_t newErrorValue                  = 0;
    Double_t correctionValue                = 1;
    
    //correct by 2pi if specified
    if (DivideBy2pi) input->Scale(1/(2*TMath::Pi()));
    if (MultiplyBy2pi) input->Scale(2*TMath::Pi());
    
    for(Int_t i=0;i<nBins;i++){
        
        //correct by 1/Pt if specified
        if(DivideByPt)    correctionValue  = 1/(input->GetBinCenter(i+1));
        if(MultiplyByPt)  correctionValue  = input->GetBinCenter(i+1);
        
        //set the value and error of the bin
        input->SetBinContent(i+1,   input->GetBinContent(i+1)*correctionValue);
        input->SetBinError(i+1,     input->GetBinError(i+1)*correctionValue);
    }
    
    return input;
}

//************************** Routine for saving histograms **********************************************************
void SaveHistos() {
    TString nameOutput                  = Form("%s/%s/GammaCocktail%s_%.2f_%s.root", fCutSelection.Data(), fEnergyFlag.Data(), fPeriodFlag.Data(), fRapidity, fCutSelection.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *outputFile                   = new TFile(nameOutput,"UPDATE");
    
    // write number of events
    histoNEvents->Write("NEvents", TObject::kOverwrite);
    
    // write binning histogram
    fDeltaPt->Write("deltaPt", TObject::kOverwrite);
    
    // write original histograms
    histoGammaSumPtY->Write(histoGammaSumPtY->GetName(), TObject::kOverwrite);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtY[i])           histoGammaPtY[i]->Write(histoGammaPtY[i]->GetName(), TObject::kOverwrite);
        if (histoGammaMotherPtY[i])     histoGammaMotherPtY[i]->Write(histoGammaMotherPtY[i]->GetName(), TObject::kOverwrite);
        if (histoGammaMotherPtPhi[i])   histoGammaMotherPtPhi[i]->Write(histoGammaMotherPtPhi[i]->GetName(), TObject::kOverwrite);
    }
    
    // write projections
    histoGammaSumPtOrBin->Write(histoGammaSumPtOrBin->GetName(), TObject::kOverwrite);
    histoGammaSumYOrBin->Write(histoGammaSumYOrBin->GetName(), TObject::kOverwrite);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i])           histoGammaPtOrBin[i]->Write(histoGammaPtOrBin[i]->GetName(), TObject::kOverwrite);
        if (histoGammaYOrBin[i])            histoGammaYOrBin[i]->Write(histoGammaYOrBin[i]->GetName(), TObject::kOverwrite);
        if (histoGammaMotherPtOrBin[i])     histoGammaMotherPtOrBin[i]->Write(histoGammaMotherPtOrBin[i]->GetName(), TObject::kOverwrite);
        if (histoGammaMotherYOrBin[i])      histoGammaMotherYOrBin[i]->Write(histoGammaMotherYOrBin[i]->GetName(), TObject::kOverwrite);
        if (histoGammaMotherPhiOrBin[i])    histoGammaMotherPhiOrBin[i]->Write(histoGammaMotherPhiOrBin[i]->GetName(), TObject::kOverwrite);
    }

    // write rebinned histograms
    histoGammaSumPt->Write(histoGammaSumPt->GetName(), TObject::kOverwrite);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPt[i])                histoGammaPt[i]->Write(histoGammaPt[i]->GetName(), TObject::kOverwrite);
        if (histoGammaMotherPt[i])          histoGammaMotherPt[i]->Write(histoGammaMotherPt[i]->GetName(), TObject::kOverwrite);
    }
    
    // write input parametrizations and mt scaled ones
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (cocktailInputParametrizations[i])           cocktailInputParametrizations[i]->Write(cocktailInputParametrizations[i]->GetName(), TObject::kOverwrite);
        if (cocktailInputParametrizationsMtScaled[i])   cocktailInputParametrizationsMtScaled[i]->Write(cocktailInputParametrizationsMtScaled[i]->GetName(), TObject::kOverwrite);
    }
}

//************************** Routine to get cocktail input file *****************************************************
TList* GetCocktailInput(TString energy, TString centrality) {
    
    if (energy.Contains("pPb")) {
        // pPb
        cocktailInputFile                       = new TFile("CocktailInput/Parametrization/CocktailInputPPb_Param.root");
    
        return NULL;

    } else if (energy.Contains("PbPb")) {
        //PbPb
        cocktailInputFile                       = new TFile("CocktailInput/Parametrization/CocktailInputPbPb_Param.root");
    
        return NULL;
        
    } else {
        // pp
        cocktailInputFile                       = new TFile("CocktailInput/Parametrization/CocktailInputPP_Param.root");
        
        if (energy.CompareTo("900GeV"))
            energy                              = "0.9TeV";
        if (energy.CompareTo("13TeV"))
            energy                              = "7TeV";
        
        if (cocktailInputFile->GetListOfKeys()->Contains(Form("pp_%s", energy.Data()))) {
            cocktailInputList                   = (TList*)cocktailInputFile->Get(Form("pp_%s", energy.Data()));
        } else {
            cout << "ERROR: Energy/collision system not found in cocktail input file!" << endl;
            cocktailInputList                   = NULL;
        }
    }

    return cocktailInputList;
}

//************************** Routine to calculate mt scaled params **************************************************
TF1* MtScaledParam(TF1* param, Int_t particleNumber) {

    if (!param)
        return NULL;
    
    Int_t collSysNumber;
    if (fEnergyFlag.Contains("PbPb"))
        collSysNumber                           = 2;
    else if (fEnergyFlag.Contains("pPb"))
        collSysNumber                           = 1;
    else
        collSysNumber                           = 0;
    
    Double_t scaleFactor                        = mtScaleFactor[collSysNumber][particleNumber];
    Double_t mass                               = GetMass(motherParticles[particleNumber]);
    Double_t massPi0                            = GetMass("Pi0");
    
    Double_t xMin, xMax;
    param->GetRange(xMin, xMax);
    
    paramScaleBase                              = param;
    
    TF1* scaledParam                            = new TF1("scaledParam",
                                                          [&](double*x, double *p)
                                                          {return p[2] * paramScaleBase->Eval(5.) / paramScaleBase->Eval(TMath::Sqrt(25. + p[0]*p[0] - p[1]*p[1])) * x[0]/TMath::Sqrt(x[0]*x[0] + p[0]*p[0] - p[1]*p[1])*paramScaleBase->Eval(TMath::Sqrt(x[0]*x[0] + p[0]*p[0] - p[1]*p[1]));},
                                                          xMin, xMax, 3);
    scaledParam->SetParameters(mass, massPi0, scaleFactor);

    return scaledParam;
}

//************************** Return partile mass ********************************************************************
Double_t GetMass(TString particleName) {
    
    TDatabasePDG* pdg                       = new TDatabasePDG;
    TParticlePDG* particle                  = NULL;
    Double_t mass                           = 0;
    
    if (particleName.CompareTo("Pi0") == 0) {
        particle                            = pdg->GetParticle(111);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("Eta") == 0) {
        particle                            = pdg->GetParticle(221);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("EtaPrim") == 0) {
        particle                            = pdg->GetParticle(331);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("omega") == 0) {
        particle                            = pdg->GetParticle(223);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("rho0") == 0) {
        particle                            = pdg->GetParticle(113);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("rho+") == 0) {
        particle                            = pdg->GetParticle(213);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("rho-") == 0) {
        particle                            = pdg->GetParticle(-213);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("phi") == 0) {
        particle                            = pdg->GetParticle(333);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("J/psi") == 0) {
        particle                            = pdg->GetParticle(443);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("Delta-") == 0) {
        particle                            = pdg->GetParticle(1114);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("Delta0") == 0) {
        particle                            = pdg->GetParticle(2114);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("Delta+") == 0) {
        particle                            = pdg->GetParticle(2214);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("Delta++") == 0) {
        particle                            = pdg->GetParticle(2224);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("Sigma0") == 0) {
        particle                            = pdg->GetParticle(3212);
        mass                                = particle->Mass();
    } else {
        mass                                = -1;
    }

    delete pdg;
    return mass;
}

//************************** Free pointer ***************************************************************************
void DeleteObjects() {
    
    delete cocktailInputFile;
    delete cocktailInputList;

    delete fDeltaPt;

    delete[] histoDecayChannels;
    delete[] histoGammaPtPhi;
    delete[] histoGammaPtOrBin;
    delete[] histoGammaPt;
    delete[] histoGammaYOrBin;
    delete[] histoGammaPhiOrBin;
    delete[] histoGammaMotherPtY;
    delete[] histoGammaMotherPtOrBin;
    delete[] histoGammaMotherPt;
    delete[] histoGammaMotherYOrBin;
    delete[] histoGammaMotherPtPhi;
    delete[] histoGammaMotherPhiOrBin;
    delete[] cocktailInputParametrizations;
    delete[] cocktailInputParametrizationsMtScaled;
}

