// provided by Gamma Conversion Group, $ALICE_ROOT/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
// ***************************************************************************************************************
// **   Friederike Bock, friederike.bock@cern.ch                                                                **
// **   Lucas Altenkaemper, lucas.altenkamper@cern.ch                                                           **
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
#include "TTree.h"
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
                        Double_t rapidity           = 0.85,
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
    
    Bool_t isPCM                                                = 0;
    Bool_t isCalo                                               = 0;
    if ( mode == 0 || mode == 2 || mode == 3)
        isPCM                                                   = 1;
    if ( mode == 2 || mode == 3 || mode == 4 || mode == 5)
        isCalo                                                  = 1;
    
    //**************************** Determine Centrality *************************************************************
    TString centrality                                          = GetCentralityString(fEventCutSelection);
    TString collisionSystem                                     = ReturnFullCollisionsSystem(option);
    TString cent                                                = "";
    TString textMeasurement                                     = ""; //"#gamma";
    TString detectionProcess                                    = ReturnFullTextReconstructionProcess(mode);
    //TString detectionProcess1                                   = "";
    //TString detectionProcess2                                   = "";
    //if (isPCM && isCalo){
    //    detectionProcess1                                       = ReturnFullTextReconstructionProcess(mode,1);
    //    detectionProcess                                        = detectionProcess1;
    //    detectionProcess2                                       = ReturnFullTextReconstructionProcess(mode,2);
    //}
    detectionProcess                                            = "";
    if(option.Contains("PbPb")){
        cent                                                    = Form("%s %s", centrality.Data(), collisionSystem.Data());
    } else {
        cent                                                    = collisionSystem;
    }
    
    //***************************** Load binning for spectrum *******************************************************
    Initialize(fEnergyFlag, numberOfBins);
    
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
        cout << "ERROR: Folder with rapidity " << rapidity << " not contained in cocktail file!" << endl;
        return;
    }
    TTree* cocktailSettingsTree                                 = (TTree*)histoListCocktail->FindObject("cocktailSettings");
    TList* cocktailSettingsList                                 = NULL;
    if (cocktailSettingsTree) cocktailSettingsList              = (TList*)cocktailSettingsTree->GetUserInfo();
    if (!cocktailSettingsTree || !cocktailSettingsList) {
        cout << "ERROR: Cocktail settings not contained in cocktail file!" << endl;
        return;
    }
    
    //***************************** read cocktail settings **********************************************************
    histMtScalingFactors                                        = (TH1F*)cocktailSettingsList->FindObject("histoMtScaleFactor");
    TString tempBinLabel                                        = "";
    for (Int_t i=1; i<histMtScalingFactors->GetNbinsX()+1; i++) {
        tempBinLabel                                            = (TString)histMtScalingFactors->GetXaxis()->GetBinLabel(i);
        for (Int_t j=0; j<nMotherParticles; j++) {
            if (tempBinLabel.CompareTo(motherParticlesPDG[j]) == 0)
                mtScaleFactor[j]                                = histMtScalingFactors->GetBinContent(i);
        }
    }
    
    TObject* tempObject                                         = NULL;
    TString tempObjectName                                      = "";
    TObjArray* arr                                              = NULL;
    for (Int_t i=0; i<cocktailSettingsList->GetEntries(); i++) {
        tempObject                                              = (TObject*)cocktailSettingsList->At(i);
        tempObjectName                                          = (TString)tempObject->GetName();
        
        if (tempObjectName.BeginsWith("selectMothers")) {
            arr                                                 = tempObjectName.Tokenize("_");
            selectedMothers                                     = (((TObjString*)arr->At(1))->GetString()).Atoi();
            for (Int_t j=0; j<nMotherParticles; j++) {
                if (selectedMothers&motherParticleDec[j]) {
                    TH2F* tempHist                              = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_Gamma_From_%s", motherParticles[j].Data()));
                    if (tempHist) {
                        if (tempHist->GetEntries())
                            hasMother[j]                        = kTRUE;
                    } else {
                        hasMother[j]                            = kFALSE;
                    }
                }
            }
        }
        if (tempObjectName.BeginsWith("nParticles")) {
            arr                                                 = tempObjectName.Tokenize("_");
            nParticles                                          = (((TObjString*)arr->At(1))->GetString()).Atoi();
        }
        if (tempObjectName.BeginsWith("ptMin")) {
            arr                                                 = tempObjectName.Tokenize("_");
            ptGenMin                                            = (((TObjString*)arr->At(1))->GetString()).Atof();
        }
        if (tempObjectName.BeginsWith("ptMax")) {
            arr                                                 = tempObjectName.Tokenize("_");
            ptGenMax                                            = (((TObjString*)arr->At(1))->GetString()).Atof();
        }
    }
    
    cocktailInputParametrizations                               = new TF1*[nMotherParticles];
    cocktailInputParametrizationsMtScaled                       = new TF1*[nMotherParticles];
    TF1* paramTemp                                              = NULL;
    TF1* paramMtTemp                                            = NULL;
    for (Int_t i=0; i<nMotherParticles; i++) {
        paramTemp                                               = (TF1*)cocktailSettingsList->FindObject(Form("%s_pt", motherParticlesPDG[i].Data()));
        paramMtTemp                                             = (TF1*)cocktailSettingsList->FindObject(Form("%s_pt_mtScaled", motherParticlesPDG[i].Data()));
        if (paramTemp)  cocktailInputParametrizations[i]        = new TF1(*paramTemp);
        else            cocktailInputParametrizations[i]        = NULL;
        if (paramMtTemp)
            cocktailInputParametrizationsMtScaled[i]            = new TF1(*paramMtTemp);
        else
            cocktailInputParametrizationsMtScaled[i]            = (TF1*)MtScaledParam(cocktailInputParametrizations[0], i);
    }
    

    for (Int_t i=0; i<nMotherParticles; i++) {
        for (Int_t j=0; j<18; j++) {
            decayChannelsBR[i][j]                               = 0.;
        }
    }
    histoDecayChannelsBR                                        = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (!hasMother[i]) continue;
        histoDecayChannelsBR[i]                                 = (TH1F*)cocktailSettingsList->FindObject(Form("PythiaBR_%s", motherParticles[i].Data()));
        for (Int_t j=0; j<18; j++) {
            if (histoDecayChannelsBR[i]->GetBinContent(j+2))
                decayChannelsBR[i][j]                           = histoDecayChannelsBR[i]->GetBinContent(j+2);
        }
    }

    //***************************** ranges **************************************************************************
    Double_t deltaRap                                           = 2*rapidity;
    Double_t deltaEta                                           = 2*0.9;                    // this must be taken from cut or sth., also if rap = 0.8 and eta = 0.9 -> we are missing photons
    Double_t deltaPtGen                                         = ptGenMax-ptGenMin;
    Double_t deltaPhi                                           = 2*TMath::Pi();
    
    //***************************** Get number of spectra ***********************************************************
    Int_t nSpectra                                              = 0;
    for (Int_t i=0; i<nMotherParticles; i++)
        if (hasMother[i]) nSpectra++;
    Int_t nRows                                                 = 0;
    if (nSpectra%6 == 0) nRows                                  = nSpectra/6 + 1;
    else nRows                                                  = (nSpectra+1)/6 + 1;

    //***************************** Pi0 file ************************************************************************
    TFile* filePi0                                              = new TFile(nameFilePi0);
    if (filePi0 && !filePi0->IsZombie()) {
        histoPi0YieldData                                       = (TH1D*)filePi0->Get("CorrectedYieldTrueEff");
    } else {
        cout << "WARNING: No pi0 file specified!" << endl;
        histoPi0YieldData                                       = NULL;
    }

    //***************************** Get number of events (cocktail) *************************************************
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
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (hasMother[i]) {
            histoDecayChannels[i]                               = (TH1F*)histoListCocktail->FindObject(Form("DecayChannels_%s", motherParticles[i].Data()));
            
            histoGammaPtY[i]                                    = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_Gamma_From_%s", motherParticles[i].Data()));
            histoGammaPtY[i]->SetName(Form("Gamma_From_%s_Pt_Y_OrBin", motherParticles[i].Data()));
            histoGammaPtY[i]->Sumw2();
            
            histoGammaPtPhi[i]                                  = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_Gamma_From_%s", motherParticles[i].Data()));
            histoGammaPtPhi[i]->SetName(Form("Gamma_From_%s_Pt_Phi_OrBin", motherParticles[i].Data()));
            histoGammaPtPhi[i]->Sumw2();
            
            histoGammaMotherPtY[i]                              = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_%s", motherParticles[i].Data()));
            histoGammaMotherPtY[i]->SetName(Form("%s_Pt_Y_OrBin", motherParticles[i].Data()));
            histoGammaMotherPtY[i]->Sumw2();
            
            histoGammaMotherPtPhi[i]                            = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_%s", motherParticles[i].Data()));
            histoGammaMotherPtPhi[i]->SetName(Form("%s_Pt_Phi_OrBin", motherParticles[i].Data()));
            histoGammaMotherPtPhi[i]->Sumw2();
        } else {
            histoDecayChannels[i]                               = NULL;
            histoGammaPtY[i]                                    = NULL;
            histoGammaPtPhi[i]                                  = NULL;
            histoGammaMotherPtY[i]                              = NULL;
            histoGammaMotherPtPhi[i]                            = NULL;
        }
        
    }

    //***************************** Get decay channels  *************************************************************
    for (Int_t i=0; i<nMotherParticles; i++) {
        for (Int_t j=0; j<18; j++) {
            decayChannelsLatex[i][j]                            = "";
        }
    }
    
    tempBinLabel                                                = "";
    Int_t counter                                               = 0;
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoDecayChannels[i]) {
            for (Int_t bin=2; bin<20; bin++) {
                if (histoDecayChannels[i]->GetBinContent(bin)) {
                    tempBinLabel                                = histoDecayChannels[i]->GetXaxis()->GetBinLabel(bin);
                    if (tempBinLabel.Contains("#gamma"))
                        decayChannelsLatex[i][bin-2]            = tempBinLabel;
                }
            }
        }
    }
    
    //***************************** Project from 2D histograms ******************************************************
    histoGammaSumPtOrBin                                        = (TH1F*)histoGammaSumPtY->ProjectionX("Gamma_Pt_OrBin", 1, histoGammaSumPtY->GetNbinsY(), "e");
    SetHistogramTitles(histoGammaSumPtOrBin,"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
    histoGammaSumPtOrBin->Sumw2();
    histoGammaSumPtOrBin->GetXaxis()->SetRangeUser(ptGenMin, ptGenMax);
    histoGammaSumPtOrBin->Scale(deltaRap*deltaPhi/deltaEta);
    histoGammaSumYOrBin                                         = (TH1F*)histoGammaSumPtY->ProjectionY("Gamma_Y_OrBin", 1, histoGammaSumPtY->GetNbinsX(), "e");
    SetHistogramTitles(histoGammaSumYOrBin,"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
    histoGammaSumYOrBin->Scale(deltaPtGen*deltaPhi/deltaEta);
    histoGammaSumPhiOrBin                                       = (TH1F*)histoGammaSumPtPhi->ProjectionY("Gamma_Phi_OrBin", 1, histoGammaSumPtPhi->GetNbinsX(), "e");
    SetHistogramTitles(histoGammaSumPhiOrBin,"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
    histoGammaSumPhiOrBin->Scale(deltaPtGen*deltaRap/deltaEta);
    
    histoGammaPtOrBin                                           = new TH1F*[nMotherParticles];
    histoGammaYOrBin                                            = new TH1F*[nMotherParticles];
    histoGammaPhiOrBin                                          = new TH1F*[nMotherParticles];
    histoGammaMotherPtOrBin                                     = new TH1F*[nMotherParticles];
    histoGammaMotherYOrBin                                      = new TH1F*[nMotherParticles];
    histoGammaMotherPhiOrBin                                    = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtY[i]) {
            histoGammaPtOrBin[i]                                = (TH1F*)histoGammaPtY[i]->ProjectionX(Form("Gamma_From_%s_Pt_OrBin", motherParticles[i].Data()), 1, histoGammaPtY[i]->GetNbinsY(), "e");
            SetHistogramTitles(histoGammaPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaPtOrBin[i]->Sumw2();
            histoGammaPtOrBin[i]->GetXaxis()->SetRangeUser(ptGenMin, ptGenMax);
            histoGammaPtOrBin[i]->Scale(deltaRap*deltaPhi/deltaEta);
            histoGammaYOrBin[i]                                 = (TH1F*)histoGammaPtY[i]->ProjectionY(Form("Gamma_From_%s_Y_OrBin", motherParticles[i].Data()), 1, histoGammaPtY[i]->GetNbinsX(), "e");
            SetHistogramTitles(histoGammaYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaYOrBin[i]->Sumw2();
            histoGammaYOrBin[i]->Scale(deltaPtGen*deltaPhi/deltaEta);
        } else {
            histoGammaPtOrBin[i]                                = NULL;
            histoGammaYOrBin[i]                                 = NULL;
        }
        if (histoGammaPtPhi[i]) {
            histoGammaPhiOrBin[i]                               = (TH1F*)histoGammaPtPhi[i]->ProjectionY(Form("Gamma_From_%s_Phi_OrBin", motherParticles[i].Data()), 1, histoGammaPtPhi[i]->GetNbinsX(), "e");
            SetHistogramTitles(histoGammaPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaPhiOrBin[i]->Sumw2();
            histoGammaPhiOrBin[i]->Scale(deltaPtGen*deltaRap/deltaEta);
        } else
            histoGammaPhiOrBin[i]                               = NULL;
        
        if (histoGammaMotherPtY[i]) {
            histoGammaMotherPtOrBin[i]                          = (TH1F*)histoGammaMotherPtY[i]->ProjectionX(Form("%s_Pt_OrBin", motherParticles[i].Data()), 1, histoGammaMotherPtY[i]->GetNbinsY(), "e");
            SetHistogramTitles(histoGammaMotherPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaMotherPtOrBin[i]->Sumw2();
            histoGammaMotherPtOrBin[i]->GetXaxis()->SetRangeUser(ptGenMin, ptGenMax);
            histoGammaMotherPtOrBin[i]->Scale(deltaRap*deltaPhi/deltaEta);
            histoGammaMotherYOrBin[i]                           = (TH1F*)histoGammaMotherPtY[i]->ProjectionY(Form("%s_Y_OrBin", motherParticles[i].Data()), 1, histoGammaMotherPtY[i]->GetNbinsX(), "e");
            SetHistogramTitles(histoGammaMotherYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaMotherYOrBin[i]->Sumw2();
            histoGammaMotherYOrBin[i]->Scale(deltaPtGen*deltaPhi/deltaEta);
        } else {
            histoGammaMotherPtOrBin[i]                          = NULL;
            histoGammaMotherYOrBin[i]                           = NULL;
        }
        if (histoGammaMotherPtPhi[i]) {
            histoGammaMotherPhiOrBin[i]                         = (TH1F*)histoGammaMotherPtPhi[i]->ProjectionY(Form("%s_Phi_OrBin", motherParticles[i].Data()), 1, histoGammaMotherPtPhi[i]->GetNbinsX(), "e");
            SetHistogramTitles(histoGammaMotherPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaMotherPhiOrBin[i]->Sumw2();
            histoGammaMotherPhiOrBin[i]->Scale(deltaPtGen*deltaRap/deltaEta);
        } else
            histoGammaMotherPhiOrBin[i]                         = NULL;
    }
    
    if (!histoGammaMotherPtOrBin[0]) {
        cout << "ERROR: Didn't get pi0 pt spectrum, returning!" << endl;
        return;
    }

    //***************************** Scale spectra *******************************************************************
    histoGammaSumPtOrBin->Scale(1./nEvents);
    histoGammaSumYOrBin->Scale(1./nEvents);
    histoGammaSumPhiOrBin->Scale(1./nEvents);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i])           histoGammaPtOrBin[i]->Scale(1./nEvents);
        if (histoGammaYOrBin[i])            histoGammaYOrBin[i]->Scale(1./nEvents);
        if (histoGammaPhiOrBin[i])          histoGammaPhiOrBin[i]->Scale(1./nEvents);
        if (histoGammaMotherPtOrBin[i])     histoGammaMotherPtOrBin[i]->Scale(1./nEvents);
        if (histoGammaMotherYOrBin[i])      histoGammaMotherYOrBin[i]->Scale(1./nEvents);
        if (histoGammaMotherPhiOrBin[i])    histoGammaMotherPhiOrBin[i]->Scale(1./nEvents);
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

    //***************************** Transform yields ****************************************************************
    histoGammaSumPt                                             = ConvertYieldHisto(histoGammaSumPt);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPt[i])        histoGammaPt[i]             = ConvertYieldHisto(histoGammaPt[i]);
        if (histoGammaMotherPt[i])  histoGammaMotherPt[i]       = ConvertYieldHisto(histoGammaMotherPt[i]);
    }
    
    TH1D* dummyHist                                             = NULL;
    //***************************** Plot cocktail mothers (pt) ******************************************************
    TCanvas *canvasMothers                                      = new TCanvas("canvasMothers","",1100,1200);
    DrawGammaCanvasSettings(canvasMothers, 0.165, 0.02, 0.02, 0.09);
    canvasMothers->SetLogy();
    canvasMothers->SetLogx();
    
    TLegend* legendMothers                                      = GetAndSetLegend2(0.2, 0.95-(0.04*nRows), 0.95, 0.95, 40, 6);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, 5e-2, ptGenMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", 1e-7, 6e2, 1.0, 1.8);
    dummyHist->Draw();
        
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherPtOrBin[i]) {
            DrawGammaSetMarker(         histoGammaMotherPtOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothers->AddEntry(    histoGammaMotherPtOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "l");
            histoGammaMotherPtOrBin[i]->SetLineWidth(2);
            histoGammaMotherPtOrBin[i]->Draw("csamehist");
        }
    }
    legendMothers->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(0.25, 0.22, 0.032, cent, textMeasurement, detectionProcess, 42, 0.03);
    
    canvasMothers->SaveAs(Form("%s/CocktailMothers_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothers;
    delete canvasMothers;

    //***************************** Plot cocktail mothers + input param (pt) ****************************************
    TCanvas *canvasMothersParam                                 = new TCanvas("canvasMothersParam","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersParam, 0.165, 0.02, 0.02, 0.09);
    canvasMothersParam->SetLogy();
    canvasMothersParam->SetLogx();
    
    TLegend* legendMothersParam                                 = GetAndSetLegend2(0.2, 0.95-(0.04*nRows), 0.95, 0.95, 40, 6);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, 5e-2, ptGenMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", 1e-7, 6e2, 1.0, 1.8);
    dummyHist->Draw();
    
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherPtOrBin[i]) {
            if (cocktailInputParametrizations[i]) {
                cocktailInputParametrizations[i]->SetLineColor(cocktailColor[i]);
                cocktailInputParametrizations[i]->SetLineStyle(2);
                cocktailInputParametrizations[i]->SetLineWidth(2);
                cocktailInputParametrizations[i]->Draw("same");
            }
            legendMothersParam->AddEntry(histoGammaMotherPtOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "l");
            histoGammaMotherPtOrBin[i]->Draw("csamehist");
        }
    }
    legendMothersParam->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(0.25, 0.22, 0.032, cent, textMeasurement, detectionProcess, 42, 0.03);
    
    canvasMothersParam->SaveAs(Form("%s/CocktailMothersInclParam_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothersParam;
    delete canvasMothersParam;
    
    //***************************** Plot ratio cocktail mothers to pi0 (pt) *****************************************
    TCanvas *canvasMothersRatio                                 = new TCanvas("canvasMothersRatio","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersRatio, 0.1, 0.02, 0.02, 0.09);
    canvasMothersRatio->SetLogy();
    canvasMothersRatio->SetLogx();
    
    TLegend* legendMothersRatio                                 = GetAndSetLegend2(0.13, 0.95-(0.04*nRows), 0.95, 0.95, 40, 5);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, 5e-2, ptGenMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "X / #pi^{0}", 1e-3, 10, 1.0, 1.2);
    dummyHist->Draw();

    TH1F* tempRatio                                             = NULL;
    for (Int_t i=1; i<nMotherParticles; i++) {
        if (histoGammaMotherPtOrBin[0] && histoGammaMotherPtOrBin[i]) {
            tempRatio                                           = (TH1F*)histoGammaMotherPtOrBin[i]->Clone("tempRatio");
            tempRatio->Sumw2();
            tempRatio->Divide(histoGammaMotherPtOrBin[0]);
            DrawGammaSetMarker(             tempRatio, cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothersRatio->AddEntry(   tempRatio, Form("%s / #pi^{0}", motherParticlesLatex[i].Data()), "l");
            tempRatio->SetLineWidth(2);
            tempRatio->Draw("csamehist");
        }
    }
    legendMothersRatio->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(0.25, 0.22, 0.032, cent, textMeasurement, detectionProcess, 42, 0.03);

    canvasMothersRatio->SaveAs(Form("%s/CocktailMothersRatioToPi0_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete tempRatio;
    delete legendMothersRatio;
    delete canvasMothersRatio;

    //***************************** Plot cocktail mothers (y) *******************************************************
    TCanvas *canvasMothersY                                     = new TCanvas("canvasMothersY","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersY, 0.165, 0.02, 0.02, 0.09);
    canvasMothersY->SetLogy();
    
    TLegend* legendMothersY                                     = GetAndSetLegend2(0.2, 0.95-(0.04*nRows), 0.95, 0.95, 40, 6);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, -fRapidity, fRapidity);
    SetHistogramm(dummyHist, "y", "#frac{1}{N_{ev}} #frac{d#it{N}}{dy}", 5e-3, 1e2, 1.0, 1.8);
    dummyHist->Draw();
    
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherYOrBin[i]) {
            DrawGammaSetMarker(         histoGammaMotherYOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothersY->AddEntry(   histoGammaMotherYOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "p");
            histoGammaMotherYOrBin[i]->Draw("same");
        }
    }
    legendMothersY->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(0.25, 0.22, 0.032, cent, textMeasurement, detectionProcess, 42, 0.03);

    canvasMothersY->SaveAs(Form("%s/CocktailMothersY_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothersY;
    delete canvasMothersY;

    //***************************** Plot cocktail mothers (phi) *****************************************************
    TCanvas *canvasMothersPhi                                   = new TCanvas("canvasMothersPhi","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersPhi, 0.165, 0.02, 0.02, 0.09);
    canvasMothersPhi->SetLogy();

    TLegend* legendMothersPhi                                   = GetAndSetLegend2(0.2, 0.95-(0.04*nRows), 0.95, 0.95, 40, 6);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, 0, 7.);
    SetHistogramm(dummyHist, "#phi", "#frac{1}{N_{ev}} #frac{d#it{N}}{d#phi}", 5e-3, 1e2, 1.0, 1.8);
    dummyHist->Draw();
    
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherPhiOrBin[i]) {
            DrawGammaSetMarker(         histoGammaMotherPhiOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothersPhi->AddEntry( histoGammaMotherPhiOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "p");
            histoGammaMotherPhiOrBin[i]->Draw("same");
        }
    }
    legendMothersPhi->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(0.25, 0.22, 0.032, cent, textMeasurement, detectionProcess, 42, 0.03);

    canvasMothersPhi->SaveAs(Form("%s/CocktailMothersPhi_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothersPhi;
    delete canvasMothersPhi;
    
    //***************************** Plot cocktail gammas (pt) *******************************************************
    TCanvas *CocktailGammas                                     = new TCanvas("CocktailGammas","",1100,1200);
    DrawGammaCanvasSettings(CocktailGammas, 0.165, 0.02, 0.02, 0.09);
    CocktailGammas->SetLogy();
    CocktailGammas->SetLogx();
    
    TLegend* legendGammas                                       = GetAndSetLegend2(0.2, 0.95-(0.05*nRows), 0.95, 0.95, 40, 6);
    legendGammas->SetHeader("#gamma from");
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, 5e-2, ptGenMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", 1e-7, 8e3, 1.0, 1.8);
    dummyHist->Draw();
    
    DrawGammaSetMarker(     histoGammaSumPtOrBin, 20, 1, kBlack,  kBlack);
    legendGammas->AddEntry( histoGammaSumPtOrBin, "all #gamma", "l");
    histoGammaSumPtOrBin->SetLineWidth(2);
    histoGammaSumPtOrBin->Draw("chistsame");
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i]) {
            DrawGammaSetMarker(         histoGammaPtOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendGammas->AddEntry(     histoGammaPtOrBin[i], motherParticlesLatex[i].Data(), "l");
            histoGammaPtOrBin[i]->SetLineWidth(2);
            histoGammaPtOrBin[i]->Draw("csamehist");
        }
    }
    legendGammas->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(0.25, 0.22, 0.032, cent, textMeasurement, detectionProcess, 42, 0.03);

    CocktailGammas->SaveAs(Form("%s/CocktailGammas_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendGammas;
    delete CocktailGammas;

    //***************************** Plot cocktail gammas to pi0 ratio ***********************************************
    TCanvas *canvasGammasRatio                                  = new TCanvas("canvasGammasRatio","",1100,1200);
    DrawGammaCanvasSettings(canvasGammasRatio, 0.12, 0.02, 0.02, 0.09);
    canvasGammasRatio->SetLogy();
        
    TLegend* legendGammasRatio                                  = GetAndSetLegend2(0.2, 0.95-(0.05*nRows), 0.95, 0.95, 40, 6);
    legendGammasRatio->SetHeader("#gamma from");
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, 5e-2, ptGenMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#gamma_{decay} / #pi^{0}", 1e-6, 5e2, 1.0, 1.3);
    dummyHist->Draw();

    TH1F* tempRatioGammas                                       = (TH1F*)histoGammaSumPtOrBin->Clone("tempRatioGammas");
    tempRatioGammas->Sumw2();
    tempRatioGammas->Divide(histoGammaMotherPtOrBin[0]);
    DrawGammaSetMarker(             tempRatioGammas, 20, 1, kBlack,  kBlack);
    legendGammasRatio->AddEntry(    tempRatioGammas, "all #gamma", "l");
    tempRatioGammas->SetLineWidth(2);
    tempRatioGammas->Draw("chistsame");
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i]) {
            tempRatioGammas                                     = (TH1F*)histoGammaPtOrBin[i]->Clone("tempRatioGammas");
            tempRatioGammas->Sumw2();
            tempRatioGammas->Divide(histoGammaMotherPtOrBin[0]);
            DrawGammaSetMarker(             tempRatioGammas, cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendGammasRatio->AddEntry(    tempRatioGammas, motherParticlesLatex[i].Data(), "l");
            tempRatioGammas->SetLineWidth(2);
            tempRatioGammas->Draw("csamehist");
        }
    }
    legendGammasRatio->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(0.25, 0.22, 0.032, cent, textMeasurement, detectionProcess, 42, 0.03);

    canvasGammasRatio->SaveAs(Form("%s/CocktailGammasRatioToPi0_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete tempRatioGammas;
    delete legendGammasRatio;
    delete canvasGammasRatio;

    //***************************** Plot cocktail gammas to all gammas ratio ****************************************
    TCanvas *canvasGammasRatio2                                 = new TCanvas("canvasGammasRatio2","",1100,1200);
    DrawGammaCanvasSettings(canvasGammasRatio2, 0.12, 0.02, 0.02, 0.09);
    canvasGammasRatio2->SetLogy();
        
    TLegend* legendGammasRatio2                                 = GetAndSetLegend2(0.2, 0.95-(0.05*nRows), 0.95, 0.95, 40, 6);
    legendGammasRatio2->SetHeader("#gamma from");
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, 5e-2, ptGenMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#gamma_{source} / #gamma_{decay}", 1e-6, 5e2, 1.0, 1.3);
    dummyHist->Draw();

    TH1F* tempRatioGammas2                                      = NULL;
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i]) {
            tempRatioGammas2                                    = (TH1F*)histoGammaPtOrBin[i]->Clone("tempRatioGammas2");
            tempRatioGammas2->Sumw2();
            tempRatioGammas2->Divide(histoGammaSumPtOrBin);
            DrawGammaSetMarker(             tempRatioGammas2, cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendGammasRatio2->AddEntry(   tempRatioGammas2, motherParticlesLatex[i].Data(), "l");
            tempRatioGammas2->SetLineWidth(2);
            tempRatioGammas2->Draw("csamehist");
        }
    }
    legendGammasRatio2->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(0.25, 0.22, 0.032, cent, textMeasurement, detectionProcess, 42, 0.03);

    canvasGammasRatio2->SaveAs(Form("%s/CocktailGammasRatioToAll_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete tempRatioGammas2;
    delete legendGammasRatio2;
    delete canvasGammasRatio2;

    //***************************** Plot cocktail gammas (y) ********************************************************
    TCanvas *canvasGammasY                                      = new TCanvas("canvasGammasY","",1100,1200);
    DrawGammaCanvasSettings(canvasGammasY, 0.165, 0.02, 0.02, 0.09);
    canvasGammasY->SetLogy();
    
    TLegend* legendGammasY                                      = GetAndSetLegend2(0.2, 0.95-(0.05*nRows), 0.95, 0.95, 40, 6);
    legendGammasY->SetHeader("#gamma from");
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, -fRapidity, fRapidity);
    SetHistogramm(dummyHist, "y", "#frac{1}{N_{ev}} #frac{d#it{N}}{dy}", 1e-5, 5e3, 1.0, 1.8);
    dummyHist->Draw();
    
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaYOrBin[i]) {
            DrawGammaSetMarker(         histoGammaYOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendGammasY->AddEntry(    histoGammaYOrBin[i], motherParticlesLatex[i].Data(), "p");
            histoGammaYOrBin[i]->Draw("same");
        }
    }
    legendGammasY->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(0.25, 0.22, 0.032, cent, textMeasurement, detectionProcess, 42, 0.03);

    canvasGammasY->SaveAs(Form("%s/CocktailGammasY_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendGammasY;
    delete canvasGammasY;
    
    //***************************** Plot cocktail mothers (phi) *****************************************************
    TCanvas *canvasGammasPhi                                    = new TCanvas("canvasGammasPhi","",1100,1200);
    DrawGammaCanvasSettings(canvasGammasPhi, 0.165, 0.02, 0.02, 0.09);
    canvasGammasPhi->SetLogy();
    
    TLegend* legendGammasPhi                                    = GetAndSetLegend2(0.2, 0.95-(0.05*nRows), 0.95, 0.95, 40, 6);
    legendGammasPhi->SetHeader("#gamma from");
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, 0, 7.);
    SetHistogramm(dummyHist, "#phi", "#frac{1}{N_{ev}} #frac{d#it{N}}{d#phi}", 1e-5, 5e3, 1.0, 1.8);
    dummyHist->Draw();
    
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherPhiOrBin[i]) {
            DrawGammaSetMarker(         histoGammaPhiOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendGammasPhi->AddEntry(  histoGammaPhiOrBin[i], motherParticlesLatex[i].Data(), "p");
            histoGammaPhiOrBin[i]->Draw("same");
        }
    }
    legendGammasPhi->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(0.25, 0.22, 0.032, cent, textMeasurement, detectionProcess, 42, 0.03);

    canvasGammasPhi->SaveAs(Form("%s/CocktailGammasPhi_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendGammasPhi;
    delete canvasGammasPhi;
    
    //***************************** Plot mT scaling cross check *****************************************************
    TCanvas* canvasMtCrossCheck                                     = NULL;
    TLegend* legendMtCrossCheck                                     = NULL;
    TPad* padMtCrossCheck                                           = NULL;
    TPad* padMtCrossCheckRatio                                      = NULL;
    dummyHist                                                       = NULL;
    TH1D* dummyHistRatio                                            = NULL;
    TH1D* tempRatio1                                                = NULL;
    TH1D* tempRatio2                                                = NULL;
    for (Int_t particle=0; particle<nMotherParticles; particle++) {
        if (histoGammaMotherPtOrBin[particle] && cocktailInputParametrizations[particle] && cocktailInputParametrizationsMtScaled[particle]) {
            canvasMtCrossCheck                                  = new TCanvas("canvasMtCrossCheck","",1100,1200);
            padMtCrossCheck                                     = new TPad("padMtCrossCheck", "", 0., 0.25, 1., 1.,-1, -1, -2);
            padMtCrossCheckRatio                                = new TPad("padMtCrossCheckRatio", "", 0., 0., 1., 0.25,-1, -1, -2);
            legendMtCrossCheck                                  = GetAndSetLegend2(0.55, 0.87-(0.045*3), 0.9, 0.87, 40);
            
            DrawGammaCanvasSettings(canvasMtCrossCheck, 0.165, 0.015, 0.025, 0.25);
            DrawGammaPadSettings(padMtCrossCheck,       0.165, 0.015, 0.025, 0.);
            DrawGammaPadSettings(padMtCrossCheckRatio,  0.165, 0.015, 0.0, 0.25);
            
            padMtCrossCheck->Draw();
            padMtCrossCheck->SetLogy();
            padMtCrossCheck->SetLogx();
            
            padMtCrossCheckRatio->Draw();
            padMtCrossCheckRatio->SetLogx();

            // dummy hist
            dummyHist                                           = new TH1D("dummyHist", "", 1000, 5e-2, ptGenMax);
            dummyHistRatio                                      = new TH1D("dummyHist", "", 1000, 5e-2, ptGenMax);
            
            // spectrum + parametrizations
            padMtCrossCheck->cd();
            SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", 5e-8, 1e2, 1.0, 1.8);
            dummyHist->Draw();

            cocktailInputParametrizationsMtScaled[particle]->SetLineColor(cocktailColor[particle]);
            cocktailInputParametrizationsMtScaled[particle]->SetLineStyle(4);
            cocktailInputParametrizationsMtScaled[particle]->SetLineWidth(2);
            
            legendMtCrossCheck->AddEntry(histoGammaMotherPtOrBin[particle], Form("%s", motherParticlesLatex[particle].Data()), "l");
            legendMtCrossCheck->AddEntry(cocktailInputParametrizations[particle], Form("%s param.", motherParticlesLatex[particle].Data()), "l");
            legendMtCrossCheck->AddEntry(cocktailInputParametrizationsMtScaled[particle], Form("%s m_{T} scaled param.", motherParticlesLatex[particle].Data()), "l");
            
            histoGammaMotherPtOrBin[particle]->Draw("csamehist");
            cocktailInputParametrizations[particle]->Draw("same");
            cocktailInputParametrizationsMtScaled[particle]->Draw("same");
            legendMtCrossCheck->Draw("same");
            
            PutProcessLabelAndEnergyOnPlot(0.55, 0.95, 0.045, cent, textMeasurement, detectionProcess, 42, 0.03);
            
            // ratios of parametrizations to spectrum
            padMtCrossCheckRatio->cd();
            SetStyleHistoTH1ForGraphs(dummyHistRatio, "#it{p}_{T} (GeV/#it{c})","#frac{spec}{param}", 0.12, 0.1, 0.12, 0.1, 1.1, 0.6, 510, 505);
            dummyHistRatio->GetXaxis()->SetLabelOffset(-0.025);
            dummyHistRatio->GetYaxis()->SetRangeUser(0,2.3);
            
            tempRatio1                                          = (TH1D*)CalculateRatioToTF1((TH1D*)histoGammaMotherPtOrBin[particle], cocktailInputParametrizations[particle]);
            tempRatio2                                          = (TH1D*)CalculateRatioToTF1((TH1D*)histoGammaMotherPtOrBin[particle], cocktailInputParametrizationsMtScaled[particle]);
            
            tempRatio1->SetLineColor(cocktailColor[particle]);
            tempRatio2->SetLineColor(cocktailColor[particle]);
            tempRatio1->SetLineStyle(2);
            tempRatio2->SetLineStyle(4);
            
            dummyHistRatio->Draw();
            tempRatio1->Draw("csamehist");
            tempRatio2->Draw("csamehist ");
            
            canvasMtCrossCheck->SaveAs(Form("%s/MtScaling%s_%.2f_%s.%s",outputDir.Data(), motherParticles[particle].Data(),fRapidity,cutSelection.Data(),suffix.Data()));
        }
    }

  
  
    //***************************** Plot pi0 from data vs. cocktail *************************************************
    if (histoPi0YieldData) {
        
        TCanvas *canvasPi0                                          = new TCanvas("canvasPi0","",1100,1200);
        DrawGammaCanvasSettings(canvasPi0, 0.165, 0.02, 0.02, 0.09);
        canvasPi0->SetLogy();
        canvasPi0->SetLogx();
        
        TLegend* legendPi0                                          = GetAndSetLegend2(0.65, 0.87-(0.045*2), 0.9, 0.87, 40);
        dummyHist                                                   = new TH1D("dummyHist", "", 1000, 1e-1, ptMax);
        SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", histoPi0YieldData->GetMinimum(0)*0.1, histoPi0YieldData->GetMaximum()*2, 1.0, 1.8);
        dummyHist->Draw();
        
        DrawGammaSetMarker(histoPi0YieldData, 24, 1, kBlack,  kBlack);
        DrawGammaSetMarker(histoGammaMotherPt[0], 20, 1, kBlue,  kBlue);

        legendPi0->AddEntry(histoPi0YieldData,      Form("%s data", motherParticlesLatex[0].Data()), "p");
        legendPi0->AddEntry(histoGammaMotherPt[0],  Form("%s cocktail", motherParticlesLatex[0].Data()), "p");

        histoPi0YieldData->Draw("same");
        histoGammaMotherPt[0]->Draw("same");
        legendPi0->Draw("same");

        PutProcessLabelAndEnergyOnPlot(0.65, 0.92, 0.03, cent, textMeasurement, detectionProcess, 42, 0.03);

        canvasPi0->SaveAs(Form("%s/Pi0DataCocktail_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
        delete legendPi0;
        delete canvasPi0;
    }
    delete dummyHist;
    
    //***************************** Save histograms *****************************************************************
    SaveHistos();
    CreateBRTableLatex();
    
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
        NewName             = Spectrum->GetName();

    Double_t newBinContent  = 0.;
    for (Int_t i=1; i<Spectrum->GetNbinsX()+1; i++) {
        newBinContent       = Spectrum->GetBinContent(i) * Spectrum->GetBinWidth(i);
        Spectrum->SetBinContent(i, newBinContent);
    }
    
    *Spectrum = *((TH1F*)Spectrum->Rebin(fNBinsPt,NewName,fBinsPt));
    Spectrum->Divide(fDeltaPt);

}

//************************** Convert yield histo ********************************************************************
TH1F* ConvertYieldHisto(TH1F* input){
    
    if (!input) {
        cout << "Error: Histogram is NULL" << endl;
        return NULL;
    }
    
    Int_t nBins                     = input->GetNbinsX();
    Double_t newValue               = 0;
    Double_t newErrorValue          = 0;
    Double_t correctionValue        = 1;
    
    // divide py 2*pi
    input->Scale(1/(2*TMath::Pi()));
    
    // divide by pT
    for(Int_t i=0;i<nBins;i++){
        correctionValue             = 1/(input->GetBinCenter(i+1));
        input->SetBinContent(i+1,   input->GetBinContent(i+1)*correctionValue);
        input->SetBinError(i+1,     input->GetBinError(i+1)*correctionValue);
    }
    
    SetHistogramTitles(input,"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
    
    return input;
}

//************************** Routine to calculate mt scaled params **************************************************
TF1* MtScaledParam(TF1* param, Int_t particleNumber) {

    if (!param || particleNumber==0)
        return NULL;
    
    Double_t scaleFactor                        = mtScaleFactor[particleNumber];
    Double_t mass                               = GetMass(motherParticles[particleNumber]);
    Double_t massPi0                            = GetMass("Pi0");
    
    if (!scaleFactor || !mass || !massPi0)
        return NULL;
    
    Double_t xMin, xMax;
    param->GetRange(xMin, xMax);
    
    paramScaleBase                              = param;
    
    TF1* scaledParam                            = new TF1("scaledParam",
                                                          [&](double*x, double *p)
                                                          {return p[2] * paramScaleBase->Eval(5.) / paramScaleBase->Eval(TMath::Sqrt(25. + p[0]*p[0] - p[1]*p[1])) * x[0]/TMath::Sqrt(x[0]*x[0] + p[0]*p[0] - p[1]*p[1])*paramScaleBase->Eval(TMath::Sqrt(x[0]*x[0] + p[0]*p[0] - p[1]*p[1]));},
                                                          xMin, xMax, 3);
    scaledParam->SetParameters(mass, massPi0, scaleFactor);
    scaledParam->SetName(Form("%s_pt_mtScaled", motherParticlesPDG[particleNumber].Data()));

    return scaledParam;
}

//************************** Calculate ratio of TF1 to TH1D *********************************************************
TH1D* CalculateRatioToTF1(TH1D* hist, TF1* func) {
    
    if (!hist) return NULL;
    
    TH1D* resultHist            = (TH1D*)hist->Clone(Form("%s_to_%s", hist->GetName(), func->GetName()));
    
    Double_t tempVal            = 0.;
    Double_t tempErr            = 0.;
    
    for (Int_t i=1; i<hist->GetNbinsX()+1; i++) {
        if (hist->GetBinContent(i)) {
            
            tempVal             = func->Integral(hist->GetXaxis()->GetBinLowEdge(i), hist->GetXaxis()->GetBinUpEdge(i));
            tempVal             = tempVal/(hist->GetBinContent(i)*hist->GetBinWidth(i));
            tempErr             = tempVal*hist->GetBinError(i)/hist->GetBinContent(i);
            
            resultHist->SetBinContent(  i, tempVal);
            resultHist->SetBinError(    i, tempErr);
            
        } else {
            resultHist->SetBinContent(  i, 0);
            resultHist->SetBinError(    i, 0);
        }
    }
    
    return resultHist;
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

//************************** Routine to set histogram titles *********************************************************
void SetHistogramTitles(TH1F* input, TString title, TString xTitle, TString yTitle) {
    
    if (!input) return;
    
    input->SetTitle(title);
    input->GetXaxis()->SetTitle(xTitle);
    input->GetYaxis()->SetTitle(yTitle);
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

//************************** Routine for saving histograms **********************************************************
void SaveHistos() {
    TString nameOutput                  = Form("%s/%s/GammaCocktail%s_%.2f_%s.root", fCutSelection.Data(), fEnergyFlag.Data(), fPeriodFlag.Data(), fRapidity, fCutSelection.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *outputFile                   = new TFile(nameOutput,"UPDATE");
    
    // write number of events
    histoNEvents->Write("NEvents", TObject::kOverwrite);
    
    // write binning histogram
    fDeltaPt->Write("deltaPt", TObject::kOverwrite);
    
    // write projections
    histoGammaSumPtOrBin->Write(histoGammaSumPtOrBin->GetName(), TObject::kOverwrite);
    histoGammaSumYOrBin->Write(histoGammaSumYOrBin->GetName(), TObject::kOverwrite);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i])           histoGammaPtOrBin[i]->Write(histoGammaPtOrBin[i]->GetName(), TObject::kOverwrite);
        if (histoGammaYOrBin[i])            histoGammaYOrBin[i]->Write(histoGammaYOrBin[i]->GetName(), TObject::kOverwrite);
        if (histoGammaPhiOrBin[i])          histoGammaPhiOrBin[i]->Write(histoGammaPhiOrBin[i]->GetName(), TObject::kOverwrite);
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

//************************** Create tex file containing BR table ****************************************************
void CreateBRTableLatex() {

    // tex file
    TString texFileName                  = Form("%s/%s/GammaCocktail%s_%.2f_%s_BranchingRatioTable.tex", fCutSelection.Data(), fEnergyFlag.Data(), fPeriodFlag.Data(), fRapidity, fCutSelection.Data());
    ofstream texFile;
    texFile.open(texFileName);

    // table
    texFile << "{\\def\\arraystretch{1.2}\\tabcolsep=5pt" << endl;
    texFile << "  \\begin{tabular}{c | c | c | c}" << endl;
    texFile << "    particle & mass (MeV) & decay & branching ratio \\\\" << endl;
    texFile << "    \\hline" << endl;
    
    TString tempMother                  = "";
    TString tempChannel                 = "";
    TString tempBR                      = "";
    Int_t counter                       = 0;

    // last mother
    Int_t lastMother                    = 0;
    for (Int_t i=nMotherParticles-1; i>=0; i--) {
        if(hasMother[i]) {
            lastMother                  = i;
            break;
        }
    }
    
    // loop over particles
    for(Int_t particle=0; particle<nMotherParticles; particle++) {
        if(!hasMother[particle]) continue;
        
        tempMother                      = motherParticlesLatex[particle];
        tempMother.ReplaceAll("#", "\\");
        
        // loop over channels
        counter                         = 0;
        for (Int_t channel=0; channel<18; channel++) {
            if (decayChannelsLatex[particle][channel]=="") continue;
            
            tempChannel                 = decayChannelsLatex[particle][channel];
            tempChannel.ReplaceAll("#", "\\");
            
            tempBR                      = Form("%.3e", decayChannelsBR[particle][channel]);
            tempBR.ReplaceAll("e", " \\times 10^{");
            tempBR                      = tempBR + "}";
            
            if (counter==0) {
                if (cocktailInputParametrizations[particle])
                    texFile << "    $" << Form("%s", tempMother.Data()) << "$ & $" << Form("%.2f", GetMass(motherParticles[particle].Data())*1e3) << "$ & ";
                else
                    texFile << "    $" << Form("%s (%.2f)", tempMother.Data(), mtScaleFactor[particle]) << "$ & $" << Form("%.2f", GetMass(motherParticles[particle].Data())*1e3) << "$ & ";
            } else
                texFile << "    &  & ";
            
            texFile << "$" << tempChannel.Data() << "$ & ";
            texFile << "$" << tempBR.Data() << "$ \\\\" << endl;
            
            counter++;
        }
    
        if (counter && particle<lastMother) texFile << "    \\hline" << endl;
    }
    
    texFile << "  \\end{tabular}" << endl;
    texFile << "}";
    texFile.close();
}

